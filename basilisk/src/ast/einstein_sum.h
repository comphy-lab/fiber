/**
# Einstein summation notation

This header file contains the functions related to the Einstein summation macro implemented by the $\texttt{Basilisk C}$ preprocessor. This macro allows the user to write tensor algebra in index notation using $\texttt{Basilisk C}$ code. The syntax follows the famous [Einstein summation convention](https://en.wikipedia.org/wiki/Einstein_notation).

## Mathematical definitions

Let us give the basis of the Einstein summation convention through the example of a simple scalar product between two second-order tensors in $3D$ space. In mathematics and physics, it is customary to write second-order tensors as $\textbf{A} = \sum_{i,j = 1}^3 A_{ij} \bm e_i\bm e_j$ where $\bm e_i$ and $\bm e_j$ are basis vectors in $3D$ space and $A_{ij}$ is the component of the tensor $\textbf A$ projected on the vectors $\bm e_i\bm e_j$. In tensor notation the inner product between two arbitrary second-order tensors: $\textbf{A}$ and $\textbf{B}$, can be written as $\textbf{C} =  \textbf{A}\cdot\textbf{B}$. If the basis is orthonormal, the components of the tensor $\textbf C$ are given by
$$ C_{ij} = \sum_{k=1}^3 A_{ik} B_{kj}.$$
Since the summation operator can be quite cumbersome in the Einstein summation convention we rewrite this operation as
$$ C_{ij} = A_{ik} B_{kj} $$
where the summation over $k$ from $1$ to $3$ is implicit.
We took the example of a simple scalar product, nevertheless note that this notation extends to any tensor operations which can be much more complex.   
Therefore, to remove any ambiguity in the Einstein summation notation, one has to follow the following rules:
1. each index can appear at most twice in any term,
2. repeated indices are implicitly summed over,
3. each term must contain identical non-repeated indices.

For example the expression: $M_{ij}v_iv_j$ is valid and indicate a double summation over $i$ and $j$, 
while the rexpression $M_{ij}v_ju_j$ is ambigous since the index $j$ appears $3$ times. 
For a more complete and rigourous definition of the Einstein summation convention we recommand the Wikipedia page: [Einstein summation convention](https://en.wikipedia.org/wiki/Einstein_notation).

## User interface

To write a scalar product in $\texttt{Basilisk C}$ one needs first to define vector and second-order tensor-like structures. The default vector structure is the `coord` struct. The second or higher-order structures must be defined by the user. 

~~~literatec
typedef struct {
  coord x,y,z;
} mytensor;

mytensor A,B,C; 
~~~

Note that the name `mytensor` is arbitrary. However, the name of the structure members must be $x$, $y$ and $z$. Following the Einstein summation convention the scalar product introduced in the previous paragraph can be written in $\texttt{Basilisk C}$ as

~~~literatec
einstein_sum(i,j,k)
  C.i.j = A.i.k*B.k.j;
~~~

With this macro, the preprocessor of $\texttt{Basilisk C}$ interprets all lines of code within the braces as tensor operations. The letters given within the parenthesis indicate the indices on which the Einstein summation takes place. In this case a summation will be applied on the index $k$ and permutations will be performed on the indices $i$ and $j$. 
To verify that the $\texttt{Basilisk C}$ preprocessor gives the desired results, one can precompile the $\texttt{Basilisk C}$ file with the command line

~~~bash
qcc -source my_file.c
~~~

The results of the precompilation stored in `_my_file.c` will be

~~~literatec
{
  C.x.x = A.x.x*B.x.x+ A.x.y*B.y.x+ A.x.z*B.z.x;
  C.x.y = A.x.x*B.x.y+ A.x.y*B.y.y+ A.x.z*B.z.y;
  C.x.z = A.x.x*B.x.z+ A.x.y*B.y.z+ A.x.z*B.z.z;
  C.y.x = A.y.x*B.x.x+ A.y.y*B.y.x+ A.y.z*B.z.x;
  C.y.y = A.y.x*B.x.y+ A.y.y*B.y.y+ A.y.z*B.z.y;
  C.y.z = A.y.x*B.x.z+ A.y.y*B.y.z+ A.y.z*B.z.z;
  C.z.x = A.z.x*B.x.x+ A.z.y*B.y.x+ A.z.z*B.z.x;
  C.z.y = A.z.x*B.x.y+ A.z.y*B.y.y+ A.z.z*B.z.y;
  C.z.z = A.z.x*B.x.z+ A.z.y*B.y.z+ A.z.z*B.z.z;
}
~~~

Note that the preprocessor duplicated the lines for each permutation of $i$ and $j$ and applied a summation over the index $k$. This macro also applies to $\texttt{Basilisk C}$ vector fields (such as the velocity field `u.x[]`) and any other user-defined structures with members named $x$, $y$ and $z$.

For higher-order rank tensors one can note that these expressions become quite cumbersome.
It is for that reason that the `einstein_sum` macro has been created. 
For more examples one can check the test file [einstein_sum.c](/src/test/einstein_sum.c).

## Specific cases 

As long as the user follows the rules mentioned above, the desired result will be obtained. 
However, the `einstein_sum` macro is somewhat more permissive than the original rules of the convention. 
For example the code

~~~literatec
U.i = M.i.j*u.j*v.j;
~~~

which is **not defined correctly** will still be compiled, and in $2D$ will give

~~~literatec
U.x = M.x.x*u.x*v.x + M.x.y*u.y*v.y;
U.y = M.y.x*u.x*v.x + M.y.y*u.y*v.y;
~~~

The use of C functions within the macro is also allowed. 
Let `function(args)` be an arbitrary C function. 
Then the expression, 

~~~literatec
U.i = function (M.i.j*u.j*v.j);
~~~

will be expanded to

~~~literatec
U.x = function (M.x.x*u.x*v.x + M.x.y*u.y*v.y);
U.y = function (M.y.x*u.x*v.x + M.y.y*u.y*v.y);
~~~

The summation is applied at the lowest level of the expression where all the indices that must be summed over are present. 
Here the summation takes place inside the parenthesis since all indices `j` are contained within them. 

A useful example of the use of a function is the computation of the L2 norm of the second-order tensor $\textbf C$

~~~literatec
L2 = sqrt (C.i.j*C.i.j);
~~~

Lastly, if no equality sign is identified, the preprocessor will not perform any summation operations. It will only carry the permutation on the current line of code. For example if one wants to print the content of a tensor one may write 

~~~literatec
einstein_sum (i,j)
  fprintf (stderr, "%g\n", C.i.j);
~~~

which gives in $2D$

~~~literatec
  fprintf (stderr, "%g\n", C.x.x);
  fprintf (stderr, "%g\n", C.x.y);
  fprintf (stderr, "%g\n", C.y.x);
  fprintf (stderr, "%g\n", C.y.y);
~~~

## Conclusion

In summary the `einstein_sum` macro follows the following steps:

1. It Identifies if a given line of code is an assignment expression with a '=' sign or not. 
    1. If it is, it identifies the indices on the right and left-hand side of the expression. 
       It then carries a summation over the indices that appear only on the right-hand side. 
    2. If it is not an equality, it does nothing at this step. 
2. It then copies the lines of codes and carries permutations on the indices present on the line (excluding the one that has been summed over). 
3. Steps 1. and 2. are repeated for all lines within the macro. 
  
<div class="message">
<div id="msg_logo"><img src="/img/warning.png"></div>
Note that these steps do not follow the original Einstein summation convention rules. 
Therefore some expressions may be ambiguous and the user is advised to always check for the pre-compiler results using `qcc -source my_file.c` before using the code.</div>

# Implementation

This structure carries information about the indices present in the expression. */

typedef struct {
  int dimension;
  char * current_id, * forbiden_id;
  int * values_sum, id_N;
} Einstein_sumData;

/**
 * @brief Appends an AST node to an additive expression list with a '+' separator.
 *
 * If the appended item is already an additive expression of the specified type, it is added directly; otherwise, it is wrapped before appending. Returns the new additive expression node.
 *
 * @param list The existing additive expression AST node.
 * @param item_sym The symbol representing the additive expression type (e.g., addition).
 * @param item The AST node to append.
 * @return Ast* The updated additive expression AST node.
 */

Ast * ast_add_list_append (Ast * list, int item_sym, Ast * item)
{
  ast_set_line (item, ast_right_terminal (list), false);
  Ast * parent = list->parent;
  int index = ast_child_index (list);
  Ast * l;
  if (item->sym == item_sym)
    l = ast_new_children (ast_new (parent, list->sym),
                          list,
                          ast_terminal_new_char (item, "+"),
                          item);
  else {
    l =  ast_new_children (ast_new (parent, list->sym),
                           list,
                           ast_terminal_new_char (item, "+"),
                           ast_new (item, item_sym));
    ast_attach (l->child[2], item);
  }
  ast_set_child (parent, index, l);
  return l;
}

/**
 * @brief Collects single-character identifiers from the macro argument list.
 *
 * Traverses the AST node representing macro arguments, appending each single-character identifier to the provided buffer. Emits an error and exits if any identifier is longer than one character.
 *
 * @param n AST node representing a macro argument.
 * @param stack Unused.
 * @param data Pointer to a character buffer where collected identifiers are appended.
 */

static void einstein_sum_id_list (Ast * n, Stack * stack, void * data)
{
  if (n->sym == sym_IDENTIFIER) {
    AstTerminal * t = ast_terminal (n);
    if (strlen (t->start) > 1){
      fprintf (stderr,
	       "%s:%d: error: the args of einstein_sum(...,%s,...) must be of length one\n",
	       t->file, t->line, t->start);
      exit (1);
    }
    char * id_list = data, first_char[2] = {*t->start, '\0'};
    strcat (id_list, first_char);
  }
}

/**
 * @brief Extracts the argument list (summation indices) from an `einstein_sum` macro invocation in the AST.
 *
 * Traverses the AST upwards from the given node to locate the enclosing `einstein_sum` macro statement, then collects its argument identifiers into a newly allocated string.
 *
 * @return Pointer to a heap-allocated string containing the concatenated summation indices. Caller is responsible for freeing the memory.
 */
static char * get_einstein_sum_args (Ast * n, Stack * stack)
{
  // stop when the macro einstein sum is found
  Ast * ein_macro = n;
  char * identifier = "init";
  while (strcmp (identifier, "einstein_sum")) {
    ein_macro = ast_parent (n, sym_macro_statement);
    identifier = ast_terminal(ast_schema (ein_macro, sym_macro_statement,
					  0, sym_MACRO))->start;
  }
  // gather the args in the buffer
  char buffer[1000] = {0}; // fixme: buffer overflows??
  ast_traverse (ast_schema (ein_macro, sym_macro_statement,
			    2, sym_argument_expression_list),
		stack, einstein_sum_id_list, buffer);
  int length = strlen (buffer);
  buffer[length] = '\0';
  char * indices = malloc ((length + 1)*sizeof(char));
  strcpy (indices, buffer);
  return indices;
}

/**
 * @brief Appends member identifier characters to a buffer if they match summation indices.
 *
 * Traverses an AST node representing a member identifier and, if the identifier's first character
 * matches any of the summation indices specified in the macro arguments, appends it to the provided buffer.
 *
 * @param n AST node to examine.
 * @param stack AST traversal stack (unused in this function).
 * @param data Pointer to a character buffer where matching identifier characters are appended.
 */
static void einstein_sum_get_member_id (Ast * n, Stack * stack, void * data)
{
  if (n->sym == sym_member_identifier) {
    AstTerminal * t = ast_terminal (ast_schema (n, 
                                                sym_member_identifier,
                                                0, sym_generic_identifier,
                                                0, sym_IDENTIFIER));
    char * sub_id_list = data;
    char * id_list = get_einstein_sum_args (n, stack);
    char first_char[2] = {*t->start, '\0'};
    if (strchr (id_list, *t->start))
      strcat (sub_id_list, first_char);
  }
}

/**
 * @brief Finds the enclosing expression statement AST node for a given node.
 *
 * Traverses the AST upwards from the given node until it reaches a node representing an expression statement or until the traversal condition fails.
 *
 * @param n The starting AST node.
 * @return Ast* Pointer to the enclosing expression statement node, or the last traversed node if not found.
 */
static Ast * get_expression_statement (Ast * n)
{
  while (n->sym != sym_expression_statement && n->sym == sym_expression)
    n = n->parent;  
  return n;
}

/**
 * @brief Returns a string containing all unique index identifiers found in the given AST subtree.
 *
 * Traverses the AST node `n` and collects single-character member identifiers (e.g., `.i`, `.j`) into a dynamically allocated string. The returned string contains the indices present in the expression, in the order they are encountered.
 *
 * @return Pointer to a newly allocated string of indices. Caller is responsible for freeing the memory.
 */

static char * get_expression_id (Ast * n, Stack * stack)
{
  char buffer[1000] = {0}; // fixme: buffer overflows??
  ast_traverse (n, stack, einstein_sum_get_member_id, buffer);
  int length = strlen (buffer);
  buffer[length + 1] = '\0';
  char * sub_id_list = malloc ((length + 2)*sizeof(char));
  strcpy (sub_id_list, buffer);
  return sub_id_list;
}

/**
 * @brief Retrieves the right-hand side expression of the nearest enclosing assignment.
 *
 * Traverses the AST upwards from the given node to locate the closest assignment expression,
 * then returns its right-hand side child node.
 *
 * @return Ast* The AST node representing the right-hand side of the assignment, or NULL if not found.
 */

static Ast * get_right_hand_side (Ast * n, Stack * stack)
{
  while (!ast_child (n, sym_assignment_operator))
    n = ast_parent (n, sym_assignment_expression);
  return ast_child (n, sym_assignment_expression);
}

/**
 * @brief Appends the suffix "_x" to member identifiers matching summation indices.
 *
 * Traverses the AST node and, for each member identifier whose name matches a single-character summation index from the macro arguments, modifies the identifier by appending "_x". This marks the index for later permutation during Einstein summation expansion.
 */

static void einstein_sum_replace_id (Ast * n, Stack * stack, void * data)
{
  if (n->sym == sym_member_identifier) {
    AstTerminal * t = ast_terminal (ast_schema (n, 
                                                sym_member_identifier,
                                                0, sym_generic_identifier,
                                                0, sym_IDENTIFIER));
    char * id_list = get_einstein_sum_args (n, stack);
    if (t->start[1] == '\0' && strchr (id_list, *t->start))
      // char * id_old
      strcat (t->start, "_x");
  }
}

/**
 * @brief Updates the suffix of member identifiers to reflect the current permutation index.
 *
 * For member identifiers ending with a suffix such as `_x`, `_y`, or `_z`, replaces the suffix with the appropriate character (`x`, `y`, or `z`) based on the current permutation values in the provided Einstein_sumData. This is used to permute tensor component indices during Einstein summation expansion.
 *
 * @param n AST node representing a member identifier.
 * @param stack Unused stack parameter.
 * @param data Pointer to Einstein_sumData containing permutation state.
 */

static void einstein_sum_rotate (Ast * n, Stack * stack, void * data)
{
  if (n->sym == sym_member_identifier) {
    AstTerminal * t = ast_terminal (ast_schema (n, 
                                                sym_member_identifier,
                                                0, sym_generic_identifier,
                                                0, sym_IDENTIFIER));
    int len = strlen (t->start);
    Einstein_sumData * d = data;
    if (len >= 2 && t->start[len - 2] == '_' && strchr ("xyz", t->start[len - 1]))
      for (int i = 0; i < d->id_N; i++)
        if(t->start[0] == d->current_id[i]){
          t->start[len - 1] = 'x' + d->values_sum[i];
        }
  }
}

/**
 * @brief Generates all possible index permutations for tensor components.
 *
 * Fills the provided array with every combination of index values for the specified number of indices and dimension, where 0, 1, and 2 correspond to x, y, and z components.
 *
 * @param LOP Output array to store the permutations. Must have size at least pow(dim, num_of_index) * num_of_index.
 * @param num_of_index Number of indices to permute.
 * @param dim Dimension of the tensor space (e.g., 2 for x/y, 3 for x/y/z).
 */

static void generates_list_of_permutations(int * LOP,int num_of_index,int dim){
    int length = pow(dim,num_of_index);
    int bitval = dim - 1;
    for(int i=0; i<num_of_index;i++)
      for(int j=0;j<length;j++){
        if(!(j % (int) pow(dim,i))){
          if(bitval == (dim - 1)) bitval =0;
          else bitval++;
        }
        LOP[j*num_of_index + i] = bitval;
      }
}

/**
 * @brief Expands summation indices in an expression according to Einstein notation.
 *
 * Identifies indices in the expression that should be summed over, generates all possible permutations for these indices, and duplicates the expression for each permutation, appending the resulting terms as an additive sum. Updates the forbidden index list to prevent repeated summation over the same indices.
 *
 * @param n The AST node representing the expression to process.
 * @param stack The AST traversal stack.
 * @param data Pointer to an Einstein_sumData struct containing summation state.
 */

static void einstein_sum_sum (Ast * n, Stack * stack, void * data)
{
  if (
    (n->sym == sym_multiplicative_expression 
    && n->parent->sym == sym_additive_expression) 
    ||n->sym == sym_additive_expression 
    || n == get_right_hand_side (n, stack)
  ) {
    Einstein_sumData * d = data;

    if (n == get_right_hand_side (n, stack)){
      while (n->sym != sym_additive_expression)
        n = ast_last_child(n);      
    }
    Ast * body = n;

    char * current_id_list = get_expression_id (body, stack);
    // indentify the list of the indicies on the right hand side of the expr 
    Ast * right_hand_side = get_right_hand_side (body, stack);
    char * id_list_after = get_expression_id (right_hand_side, stack);
    // list that contain index to be summed
    char sum_id[80];
    int j = 0;
    for(int i = 0; i < strlen(current_id_list); i++) {
      int k = 0, m = 0;
      // check how many times the index appears in the expr
      for (int l = 0; l < strlen (current_id_list); l++)
        if (current_id_list[l] == current_id_list[i])
          k++;
      // check how many times the index appears in the original expr
      for (int l = 0; l < strlen (id_list_after); l++)
        if (id_list_after[l] == current_id_list[i])
          m++;
      
      if (!strchr (d->forbiden_id, current_id_list[i]) && 
	  !strchr (sum_id, current_id_list[i]) && k == m)
        sum_id[j++] = current_id_list[i];
      
      sum_id[j] = '\0';
    }
    d->id_N = strlen (sum_id);
    if (!j) {
      d->id_N = 0;
      sum_id[0] = '\0';
    }
    d->current_id = malloc (sizeof(char)*(d->id_N+1));
    for(int i=0;i<d->id_N;i++) d->current_id[i] = sum_id[i];
    d->current_id[d->id_N+1] = '\0';


    int number_of_terms = pow(d->dimension,d->id_N);
    int * list_of_permutations = malloc(sizeof(int) * number_of_terms * d->id_N);
    generates_list_of_permutations(list_of_permutations,
                                  d->id_N,
                                  d->dimension);

    // we must permute and add alll multiplicative expr in the additive expr 
    int mult=0;
    if(body->sym==sym_multiplicative_expression){
      body = body->parent;
      mult=1;
    }
    Ast * b = body; 
    while (b->sym==sym_additive_expression 
      &&  ast_last_child(b)->sym==sym_multiplicative_expression){
      Ast * copy = ast_last_child(ast_copy(b)); //first multiplicative 
      for(int j=1;j<number_of_terms;j++){//we start at j=1 
        d->values_sum = &list_of_permutations[j * d->id_N];
        copy = ast_copy (copy);
        stack_push (stack, &copy);
        ast_traverse (copy, stack, einstein_sum_rotate, d);
        ast_pop_scope (stack, copy);  
        body = ast_add_list_append (body, copy->sym, copy);
      }
      b = b->child[0];
      if(mult) break;
    }
    free(list_of_permutations);
    strcat (d->forbiden_id, sum_id);
  }
}

/**
 * @brief Restores original index names in member identifiers after permutation expansion.
 *
 * Reverts temporary suffixes (e.g., `_x`, `_y`, `_z`) on member identifiers back to their original single-character index names within the AST subtree rooted at `n`.
 */
static void einstein_sum_replace_id_back (Ast * n, Stack * stack, void * data)
{
  if (n->sym == sym_member_identifier) {
    AstTerminal * t = ast_terminal (ast_schema (n, 
                                                sym_member_identifier,
                                                0, sym_generic_identifier,
                                                0, sym_IDENTIFIER));
    char * id_list = get_einstein_sum_args (n, stack);
    int len = strlen (t->start);
    if (t->start[len - 2] == '_' &&
	strchr ("xyz", t->start[len - 1]) &&
	strchr (id_list, t->start[len - 3]))
      t->start[1] = '\0', t->start[0] = t->start[len - 1] ;
  }
}

/**
 * @brief Expands an assignment expression within the Einstein summation macro into explicit summations and permutations.
 *
 * For each assignment expression inside the macro block, this function:
 * - Marks tensor indices for permutation by appending suffixes (e.g., `.i` to `.i_x`).
 * - Identifies summation indices (those appearing only on the right-hand side) and expands the expression into a sum over all possible values of these indices.
 * - Removes duplicate indices from the left-hand side.
 * - Generates all permutations of remaining free indices and duplicates the expression accordingly.
 * - Restores original member names after expansion.
 *
 * The expanded expressions are inserted into the AST, replacing the original macro statement.
 */
static void einstein_sum_expression (Ast * n, Stack * stack, void * data)
{ 
  if (n->sym == sym_assignment_expression &&
      get_expression_statement(n->parent)->sym == sym_expression_statement &&
      get_expression_statement(n->parent)->parent->sym == sym_statement &&
      n->parent->sym == sym_expression){
    Einstein_sumData * d = data;
    
    // transform into block item 
    Ast * statement = get_expression_statement (n->parent)->parent;
    Ast * item = ast_block_list_get_item (statement);
    if (!item) item = get_expression_statement (n->parent);
    Ast * item_list = item->parent;

    /**
    This is the first step, it just permute every indices such as .i into 
    .i_x and so on. */
    
    stack_push (stack, &item);
    ast_traverse (item, stack, einstein_sum_replace_id, d);
    ast_pop_scope (stack, item);
    
    /**
    The second step identifies and performs the summation. */
    
    char * id_list_before = get_expression_id (n, stack);
    int id_N = strlen (id_list_before);
    // check if the expression is an equallity  
    if (ast_child (n, sym_assignment_operator)) {
      Ast * right_hand_side = ast_child (n, sym_assignment_expression);
      Ast * left_hand_side = ast_child (n, sym_unary_expression);
      // indices appearing on the left hand side of the expression
      id_list_before = get_expression_id (left_hand_side, stack);
      // save the indices appearing on the left hand side
      d->forbiden_id = malloc (sizeof(char)*80);
      strcpy (d->forbiden_id, id_list_before);
      // traverse over the whole expression and performs the sum
      stack_push (stack, &right_hand_side);
      ast_traverse (right_hand_side, stack, einstein_sum_sum, d);
      ast_pop_scope (stack, right_hand_side);  
      id_N = strlen (id_list_before);
    }

    // Remove any duplicates in "id_list_before"
    for (int i = 0; i < id_N; i++)
      for (int j = i + 1; j < id_N;){
        if (id_list_before[i] == id_list_before[j]) {
          for (int k = j; k < id_N; k++){
            id_list_before[k] = id_list_before[k + 1];
            id_N--;
          }
        }
	else
          j++;
      }
    d->id_N = id_N; 
    d->current_id = malloc (sizeof(char)*(d->id_N +1));
    for(int i=0;i<d->id_N;i++) d->current_id[i] = id_list_before[id_N - i - 1];
    d->current_id[d->id_N+1] = '\0';

    /**
    The final step is to perform permutations on the remaining indices.*/

    int number_of_terms = pow(d->dimension,d->id_N);
    int * list_of_permutations = malloc(sizeof(int) * number_of_terms * d->id_N);
    generates_list_of_permutations(list_of_permutations,
                                  d->id_N,
                                  d->dimension);
    
    Ast * copy = statement;
    for(int j=1;j<number_of_terms;j++){//we start at j=1 
      d->values_sum = &list_of_permutations[j * d->id_N];
      copy = ast_copy (copy);
      stack_push (stack, &copy);
      ast_traverse (copy, stack, einstein_sum_rotate, d);
      ast_pop_scope (stack, copy);  
      item_list = ast_block_list_append (item_list, sym_block_item, copy);
    }
    free(list_of_permutations);
    // replace the i_x indicators by the actual member name x,yz...
    stack_push (stack, &item_list);
    ast_traverse (item_list, stack, einstein_sum_replace_id_back, d);
    ast_pop_scope (stack, item_list);  
  }
}

/**
 * @brief Expands the `einstein_sum` macro in the AST to explicit tensor summation and permutation expressions.
 *
 * Processes the `einstein_sum` macro invocation by extracting summation indices, ensuring the macro body is a compound statement, and traversing its AST to expand index notation into explicit component-wise summations and permutations according to the specified tensor dimension. The macro body is replaced with the expanded AST reflecting the full summation and permutation logic.
 *
 * @param n AST node representing the macro statement.
 * @param stack Stack used for AST traversal and scope management.
 * @param dimension Dimension of the tensor space (e.g., 2 for 2D, 3 for 3D).
 */
void einstein_sum_global (Ast * n, Stack * stack, int dimension)
{
  Ast * item = ast_block_list_get_item (ast_ancestor (n, 2));
  Einstein_sumData data = {dimension, NULL, NULL}; 
  if (!ast_schema (n, sym_macro_statement,
		   2, sym_argument_expression_list)) {
    AstTerminal * t = ast_left_terminal (n);
    fprintf (stderr,
	     "%s:%d: error: missing summation indices in macro einstein_sum(...)\n",
	     t->file, t->line);
    exit (1);
  }
  Ast * body = ast_child (n, sym_statement);
  if (!ast_schema (body, sym_statement,
		   0, sym_compound_statement)) {
    AstTerminal * open = NCB(body, "{"), * close = NCA(body, "}");
    int index = ast_child_index (body);
    body = NN(item, sym_statement,
	      NN(item, sym_compound_statement,
		 open,
		 NN(item, sym_block_item_list,
		    NN(item, sym_block_item,
		       body)),
		 close));
    ast_set_child (n, index, body);
  }
  stack_push (stack, &body);
  ast_traverse (body, stack, einstein_sum_expression, &data);
  ast_pop_scope (stack, body);
  ast_replace_child (item, 0, body);
}
