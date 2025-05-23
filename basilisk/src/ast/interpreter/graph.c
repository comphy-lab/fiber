typedef struct _Node Node;

typedef struct {
  Node * n;
  float c;
} Edge;

struct _Node {
  Key * k;
  Edge * parent, * child;
};

KHASH_MAP_INIT_INT(GRAPH, Node *)
  
typedef struct {
  khash_t(GRAPH) * hash;
} Graph;

/**
 * @brief Checks if two node pointers refer to the same node.
 *
 * @param a First node pointer.
 * @param b Second node pointer.
 * @return true if both pointers are identical; false otherwise.
 */
static inline
bool identical_nodes (const Node * a, const Node * b)
{
  return a == b;
}

/**
 * @brief Retrieves a node from the graph by key, creating it if it does not exist.
 *
 * If a node with the specified key is present in the graph, returns it. Otherwise, allocates and initializes a new node with the given key, inserts it into the graph, and returns the new node.
 *
 * @param g Pointer to the graph.
 * @param key Pointer to the key identifying the node.
 * @return Pointer to the retrieved or newly created node.
 */
static
Node * get_node (Graph * g, Key * key)
{
  khiter_t k = kh_get (GRAPH, g->hash, key->j);
  if (k != kh_end (g->hash))
    return kh_value (g->hash, k);
  Node * n = calloc (1, sizeof (Node));
  n->k = key;
  int ret;
  k = kh_put (GRAPH, g->hash, key->j, &ret);
  assert (ret > 0);
  kh_value (g->hash, k) = n;
  return n;
}

/**
 * @brief Counts the number of edges in a null-terminated edge list.
 *
 * @param list Pointer to the first edge in the list.
 * @return int The number of edges in the list, or 0 if the list is NULL.
 */
static
int nedge (const Edge * list)
{
  int n = 0;
  if (list) for (const Edge * e = list; e->n; e++, n++);
  return n;
}

/**
 * @brief Removes the edge in the list that points to the specified node.
 *
 * Searches the edge list for an edge targeting node `b`, removes it by shifting subsequent edges, and returns its cost.
 * Returns 0 if no such edge is found.
 *
 * @param list Null-terminated list of edges.
 * @param b Node to remove from the edge list.
 * @return float Cost of the removed edge, or 0 if not found.
 */
static
float remove_edge (Edge * list, Node * b)
{
  if (list)
    for (Edge * e = list; e->n; e++)
      if (identical_nodes (e->n, b)) {
	float c = e->c;
	for (; e->n; e++)
	  *e = *(e + 1);
	return c;
      }
  return 0.;
}

/**
 * @brief Adds or updates an edge to a node in a null-terminated edge list.
 *
 * If an edge to node @p b exists, increments its cost by @p c and removes the edge if the resulting cost is zero. If no such edge exists, appends a new edge to @p b with cost @p c. The edge list remains null-terminated.
 *
 * @param list Null-terminated list of edges.
 * @param b Target node for the edge.
 * @param c Cost to add to the edge.
 * @return Updated null-terminated edge list.
 */
static
Edge * add_mono_edge (Edge * list, Node * b, float c)
{
  if (list) {
    for (Edge * e = list; e->n; e++)
      if (identical_nodes (e->n, b)) {
	e->c += c;
	if (!e->c)
	  for (; e->n; e++)
	    *e = *(e + 1);
	return list;
      }
  }
  int n = nedge (list);
  list = realloc (list, (n + 2)*sizeof(Edge));
  list[n + 1].n = NULL;
  list[n].n = b;
  list[n].c = c;
  return list;
}

/**
 * @brief Adds a directed edge from node a to node b with the specified cost.
 *
 * Updates both the child edge list of node a and the parent edge list of node b to maintain bidirectional consistency.
 * If an edge already exists, its cost is incremented by c; if the resulting cost is zero, the edge is removed.
 *
 * @param a Source node.
 * @param b Destination node.
 * @param c Edge cost; must satisfy |c| < 100.
 */
static
void add_edge (Node * a, Node * b, float c)
{
  assert (fabs (c) < 100);
  a->child = add_mono_edge (a->child, b, c);
  b->parent = add_mono_edge (b->parent, a, c);
}

/**
 * @brief Retrieves the key at the specified row and column in the system's matrix.
 *
 * Searches the given row of the system for a key with column index `j`. Returns the matching key if found, or NULL if the row is empty or no such key exists.
 *
 * @param s Pointer to the system containing the matrix.
 * @param i Row index.
 * @param j Column index.
 * @return Key* Pointer to the matching key, or NULL if not found.
 */
static Key * matrix_key (const System * s, int i, int j)
{
  Dimension ** r = s->r + i;
  if (!r[0]->c)
    return NULL;
  foreach_key (r[0], c)
    if (c->j == j)
      return c;
  return NULL;
}

/**
 * @brief Returns the smallest column index among keys in the given dimension.
 *
 * Iterates through all keys in the dimension and finds the minimum value of the `j` index.
 *
 * @param d Pointer to the dimension containing a null-terminated array of keys.
 * @return int The smallest column index (`j`) found, or `INT_MAX` if the dimension is empty.
 */
static
int leftmost (const Dimension * d)
{
  int left = INT_MAX;
  for (Key ** c = d->c; c[0]; c++)
    if (c[0]->j < left)
      left = c[0]->j;
  return left;
}

/**
 * @brief Constructs a directed weighted graph from a system of constraints.
 *
 * For each constraint in the system, creates a node for the leftmost key and adds directed edges from other keys in the constraint to this node. Edge weights are determined by the ratio of matrix coefficients.
 *
 * @param s Pointer to the system of constraints to convert.
 * @return Pointer to the constructed graph representing the system.
 */
Graph * system_to_graph (const System * s)
{
  Graph * g = malloc (sizeof (Graph));
  g->hash = kh_init (GRAPH);
  int n = 0;
  foreach_constraint (s, i) {
    if (i->c) {
      int left = leftmost (i);
      Key * kleft = matrix_key (s, n, left);
      float coef = - matrix (s, n, left);
      Node * node = get_node (g, kleft);
      foreach_key (i, c)
	if (c->j != left)
	  add_edge (get_node (g, c), node, matrix (s, n, c->j)/coef);
    }
    n++;
  }
  return g;
}

/**
 * @brief Verifies the consistency of a node's parent and child edge lists.
 *
 * Asserts that for each child of the node, the node appears exactly once in the child's parent edge list, and for each parent, the node appears exactly once in the parent's child edge list.
 */
void check_node (const Node * n)
{
  if (n->child)
    for (const Edge * e = n->child; e->n; e++) {
      assert (e->n->parent);
      int np = 0;
      for (const Edge * f = e->n->parent; f->n; f++)
	if (identical_nodes (f->n, n))
	  np++;
      assert (np == 1);
    }
  if (n->parent)
    for (const Edge * e = n->parent; e->n; e++) {
      assert (e->n->child);
      int np = 0;
      for (const Edge * f = e->n->child; f->n; f++)
	if (identical_nodes (f->n, n))
	  np++;
      assert (np == 1);
    }
}

/**
 * @brief Removes a node from the graph and frees its memory.
 *
 * Deletes the node from the graph's hash map and releases all memory associated with the node, including its child and parent edge lists.
 */
static
void remove_node_internal (Graph * g, Node * node)
{
  khiter_t k = kh_get (GRAPH, g->hash, node->k->j);
  kh_del (GRAPH, g->hash, k);
  free (node->child);
  free (node->parent);
  free (node);  
}

/**
 * @brief Collapses a node into its parent, redirecting all child edges to the parent.
 *
 * Removes the edge from the parent to the node, then for each child of the node,
 * removes the corresponding parent edge and adds a new edge from the parent to the child
 * with an adjusted cost. Checks the consistency of the parent and deletes the node from the graph.
 */
static
void edge_collapse_parent (Graph * g, Node * node, Node * parent)
{
  assert (remove_edge (parent->child, node));
  if (node->child)
    for (Edge * e = node->child; e->n; e++) {
      assert (remove_edge (e->n->parent, node));
      add_edge (parent, e->n, e->c*node->parent->c);
    }
  check_node (parent);
  remove_node_internal (g, node);
}

/**
 * @brief Collapses a node into its child, redirecting edges and updating costs.
 *
 * Removes the edge from the specified child to the node, then redirects all parent edges of the node to the child, adjusting their costs accordingly. For each child of the node (except the specified child), redirects edges from the child to those nodes with updated costs. Ensures the child's edge consistency and removes the original node from the graph.
 */
static
void edge_collapse_child (Graph * g, Node * node, Node * child)
{
  float c = remove_edge (child->parent, node);
  assert (c);
  if (node->parent)
    for (Edge * e = node->parent; e->n; e++) {
      assert (remove_edge (e->n->child, node));
      add_edge (e->n, child, e->c*c);
    }
  if (node->child)
    for (Edge * e = node->child; e->n; e++)
      if (!identical_nodes (e->n, child)) {
	assert (remove_edge (e->n->parent, node));
	add_edge (child, e->n, e->c*c);
      }
  check_node (child);
  remove_node_internal (g, node);
}

/**
 * @brief Removes a node from the graph and deletes all edges connected to it.
 *
 * For the specified node, removes all references from its child and parent nodes, then deletes the node from the graph.
 */
static
void remove_node (Graph * g, Node * node)
{
  if (node->child)
    for (Edge * e = node->child; e->n; e++)
      assert (remove_edge (e->n->parent, node));
  if (node->parent)
    for (Edge * e = node->parent; e->n; e++)
      assert (remove_edge (e->n->child, node));
  remove_node_internal (g, node);
}

/**
 * @brief Simplifies the graph by removing or collapsing nodes based on structural and terminal value conditions.
 *
 * Iterates through all nodes in the graph and performs the following simplifications:
 * - Removes nodes with no children if their key's right terminal value is zero.
 * - Collapses nodes with exactly one parent if their key's right terminal value is zero.
 * - Collapses nodes with exactly one parent if the parent's key's right terminal value is zero.
 *
 * Each simplification updates the graph structure and maintains edge consistency.
 *
 * @param g Pointer to the graph to be simplified.
 * @return int The number of simplifications performed.
 */
int graph_simplify (Graph * g)
{
  int ns = 0;
  for (khiter_t k = kh_begin (g->hash); k != kh_end (g->hash); ++k)
    if (kh_exist (g->hash, k)) {
      Node * node = kh_value (g->hash, k);
      check_node (node);
      if (!atof (ast_right_terminal (node->k->parent)->start)) {
#if 1	
	if (nedge (node->child) == 0) {
	  remove_node (g, node);
	  ns++;
	}
#endif
#if 1
	else if (nedge (node->parent) == 1) {
	  edge_collapse_parent (g, node, node->parent->n);
	  ns++;
	}
#endif
      }
#if 1      
      else if (nedge (node->parent) == 1 &&
	       !atof (ast_right_terminal (node->parent->n->k->parent)->start)) {
	edge_collapse_child (g, node->parent->n, node);
	ns++;	  	
      }
#endif
    }
  return ns;
}

/**
 * @brief Prints a graph node label in DOT format using the key's index and label.
 *
 * Outputs the node's identifier and label to the specified file stream, formatted for DOT graph visualization.
 */
static
void print_node_label (const Key * k, FILE * fp)
{
  fprintf (fp, "  %d [label=\"", k->j);
  print_key_label (k, '\n', fp, LINENO);
#if 0  
  if (should_be_dimensionless (k))
    fputs ("\" style=filled fillcolor=\"aquamarine", fp);
#endif
  fprintf (fp, " %d\"]\n", k->j);
}

/**
 * @brief Outputs the graph in DOT format for visualization.
 *
 * Prints all nodes with at least one parent or child edge, and outputs directed edges with their weights as labels to the specified file stream.
 *
 * @param g Pointer to the graph to be visualized.
 * @param fp File stream to which the DOT representation is written.
 */
void graph_dot (const Graph * g, FILE * fp)
{
  fputs ("digraph mygraph {\n", fp);
  int key;
  Node * node;
  kh_foreach (g->hash, key, node, {
      if (node->child || node->parent) // do not print "orphaned" nodes
	print_node_label (node->k, fp);
      if (node->child)
	for (Edge * e = node->child; e->n; e++)
	  fprintf (fp, "  %d -> %d [label=\"%g\"]\n",
		   key, e->n->k->j, e->c);
#if 0  
      if (node->parent)
	for (Edge * e = node->parent; e->n; e++)
	  fprintf (fp, "  %d -> %d [label=\"p%g\"]\n",
		   e->n->k->index, key, e->c);
#endif
    });
  fputs ("}\n", fp);
}
