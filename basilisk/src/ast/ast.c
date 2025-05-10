#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <stdarg.h>

#include "parser.h"
#include "basilisk.h"
#include "symbols.h"

Ast * const ast_placeholder = (Ast *) 128;

/**
 * @brief Attaches multiple child AST nodes to a parent node.
 *
 * Allocates and sets the parent's child array to include all provided child nodes, updating each child's parent pointer. The list of children is specified as variadic arguments and must be terminated with NULL.
 *
 * @param parent The parent AST node to which children will be attached.
 * @return Ast* The parent node with its children attached.
 */
Ast * ast_new_children_internal (Ast * parent, ...)
{
  va_list ap;
  va_start (ap, parent);
  int n = 1;
  Ast * c = va_arg (ap, Ast *);
  while (c) {
    n++;
    c = va_arg (ap, Ast *);
  }
  va_end (ap);

  if (n > 1) {
    parent->child = allocate (ast_get_root (parent)->alloc, n*sizeof (Ast *));
    parent->child[n - 1] = NULL;
    va_start (ap, parent);
    n = 0;
    c = va_arg (ap, Ast *);
    while (c) {
      ast_set_child (parent, n++, c);
      c = va_arg (ap, Ast *);
    }
    va_end (ap);
  }
  
  return parent;
}

/**
 * @brief Recursively clears file and line information from all terminal nodes in an AST subtree.
 *
 * Sets the file pointer to NULL and the line number to 0 for each terminal node under the given AST node.
 */
static void cancel_file_line (Ast * n)
{
  AstTerminal * t = ast_terminal(n);
  if (t)
    t->file = NULL, t->line = 0;
  else
    for (Ast ** c = n->child; *c; c++)
      cancel_file_line (*c);
}

/**
 * @brief Attaches a single AST node as the last descendant child of a parent.
 *
 * Traverses to the deepest child of the given parent and attaches the node `n` as its only child, allocating space for two child pointers. Returns the parent node where the attachment occurred.
 *
 * @param parent The AST node to which the child will be attached.
 * @param n The AST node to attach as a child.
 * @return Ast* The parent node after the child has been attached.
 */
static Ast * ast_attach_single (Ast * parent, Ast * n)
{
  AstRoot * root = ast_get_root (parent);
  while (parent->child)
    parent = parent->child[0];
  parent->child = allocate (root->alloc, 2*sizeof (Ast *));
  parent->child[0] = NULL;
  parent->child[1] = NULL;
  ast_set_child (parent, 0, n);
  return parent;
}

/**
 * @brief Attaches a sequence of AST nodes as nested children to a parent node.
 *
 * Each node in the variadic argument list is attached as the sole child of the previous node, forming a linear chain of parent-child relationships starting from the given parent. The argument list must be terminated with NULL.
 *
 * @param parent The root AST node to which the chain will be attached.
 * @return Ast* The original parent node.
 */
Ast * ast_attach_internal (Ast * parent, ...)
{
  va_list ap;
  va_start (ap, parent);
  Ast * p = parent, * n = va_arg (ap, Ast *);
  while (n) {
    ast_attach_single (p, n);
    p = n;
    n = va_arg (ap, Ast *);
  }
  va_end (ap);
  return parent;
}

/**
 * @brief Recursively frees memory associated with an AST node and its descendants.
 *
 * Releases all memory held by the AST subtree rooted at the given node, including terminal strings and root-level resources such as stacks and allocators.
 */
static void ast_destroy_internal (Ast * n)
{
  if (n->child) {
    for (Ast ** c = n->child; *c; c++)
      if (*c != ast_placeholder)
	ast_destroy_internal (*c);
  }
  else {
    AstTerminal * t = ast_terminal (n);
    free (t->before);
    free (t->start);
    free (t->after);
    t->before = t->start = t->after = NULL;
  }
  AstRoot * r = ast_root (n);
  if (r) {
    free (r->before);
    free (r->after);
    if (r->stack)
      stack_destroy (r->stack);
    if (r->alloc)
      free_allocator (r->alloc);
  }
}

/**
 * @brief Recursively destroys an AST node and its subtree.
 *
 * If the node has a parent, replaces its entry in the parent's child array with a placeholder before freeing all associated memory for the node and its descendants.
 */
void ast_destroy (Ast * n)
{
  if (!n || n == ast_placeholder)
    return;
  if (n->parent && n->parent->child) {
    Ast ** c;
    for (c = n->parent->child; *c && *c != n; c++);
    if (*c == n)
      *c = ast_placeholder;
  }
  ast_destroy_internal (n);
}

/**
 * @brief Returns the line number of a terminal node as a string.
 *
 * Converts the line number of the given AST terminal node to a static string buffer.
 *
 * @param t Pointer to the AST terminal node.
 * @return Pointer to a static string containing the line number.
 *
 * @note The returned string is stored in a static buffer and will be overwritten by subsequent calls.
 */
char * ast_line (AstTerminal * t)
{
  static char s[20];
  snprintf (s, 19, "%d", t->line);
  return s;
}

/**
 * @brief Returns the root node of the AST containing the given node.
 *
 * Traverses parent pointers from the specified AST node up to the root and returns the corresponding AstRoot structure.
 *
 * @param n Pointer to an AST node.
 * @return AstRoot* Pointer to the root of the AST tree containing the node.
 */
AstRoot * ast_get_root (const Ast * n)
{
  const Ast * root = n;
  while (root->parent) {
    assert (root != root->parent);
    root = root->parent;
  }
  return ast_root (root);
}

/**
 * @brief Counts the number of newline characters in a string.
 *
 * @param s Input string to scan.
 * @return int Number of lines, determined by the count of '\n' characters.
 */
static int count_lines (const char * s)
{
  int line = 0;
  while (*s != '\0') {
    if (*s == '\n')
      line++;
    s++;
  }
  return line;
}

/**
 * @brief Sets the specified child of an AST node, updating parent relationships.
 *
 * Replaces the child at the given index in the parent's child array with the provided child node.
 * If the child node already has a parent, it is removed from its previous parent's child array and replaced with a placeholder.
 * Updates the child's parent pointer to the new parent.
 *
 * @param parent The AST node whose child is to be set.
 * @param index The index in the parent's child array to set.
 * @param child The AST node to set as the child.
 */
void ast_set_child (Ast * parent, int index, Ast * child)
{
  if (!child)
    return;
  if (child == ast_placeholder) {
    parent->child[index] = child;
    return;
  }
  if (child->parent) {
    Ast * oldparent = child->parent;
    if (oldparent->child) {
      Ast ** c;
      for (c = oldparent->child; *c && *c != child; c++);
      if (*c == child)
	*c = ast_placeholder;
    }
  }
  parent->child[index] = child;
  child->parent = parent;
}

/**
 * @brief Returns the leftmost terminal node with file and line information in an AST subtree.
 *
 * Recursively searches the subtree rooted at @p n for the first terminal node that contains valid file and line metadata. Returns NULL if no such terminal exists.
 *
 * @param n The root of the AST subtree to search.
 * @return AstTerminal* Pointer to the leftmost terminal with file and line info, or NULL if not found.
 */
static AstTerminal * ast_left_line (Ast * n)
{
  AstTerminal * t = ast_terminal (n);
  if (t)
    return t->file ? t : NULL;
  else
    for (Ast ** c = n->child; *c; c++) {
      t = ast_left_line (*c);
      if (t)
	return t;
    }
  return NULL;
}

/**
 * @brief Replaces a child node of an AST parent with a new child, transferring file and line information.
 *
 * Replaces the child at the specified index in the parent's child array with the given child node. If the old child is not a placeholder, its file and line information, as well as its 'before' string, are transferred to the new child's leftmost terminal if missing. The old child is then destroyed.
 *
 * @param parent The AST parent node whose child will be replaced.
 * @param index The index of the child to replace.
 * @param child The new AST node to insert as the child.
 */
void ast_replace_child (Ast * parent, int index, Ast * child)
{
  Ast * oldchild = parent->child[index];
  ast_set_child (parent, index, child);
  AstTerminal * left = ast_left_terminal (child);
  if (oldchild != ast_placeholder) {
    AstTerminal * oldleft = ast_left_terminal (oldchild);
    if (oldleft && left) {
      if (!left->file) {
	if (left->before)
	  free (left->before), left->before = NULL;
	assert (!left->before);
	left->file = oldleft->file;
	left->line = oldleft->line;
	ast_set_line (child, left, false);
      }
      if (left->before)
	str_append (left->before, oldleft->before);
      else
	left->before = oldleft->before, oldleft->before = NULL;
    }
    ast_destroy (oldchild);
  }
  if (left && !left->file) {
    AstTerminal * line = ast_left_line (child);
    if (line) {
      if (left->before)
	str_append (left->before, line->before);
      else
	left->before = line->before, line->before = NULL;      
      left->file = line->file;
      left->line = line->line;
      ast_set_line (child, left, false);
    }      
  }
}

/**
 * @brief Replaces a terminal node with a specified AST subtree within an AST.
 *
 * Searches the AST subtree rooted at @p n for a terminal node whose start string matches @p terminal. If found, replaces the nearest ancestor node with the same symbol as @p with by @p with, destroys the replaced node, and returns the rightmost terminal of @p with. If not found, recursively searches child nodes. Updates line information for affected nodes.
 *
 * @param n The root of the AST subtree to search.
 * @param terminal The start string of the terminal node to replace.
 * @param with The AST subtree to insert in place of the matched node.
 * @return AstTerminal* The rightmost terminal of the inserted subtree if replacement occurs, or NULL if no match is found.
 */
AstTerminal * ast_replace (Ast * n, const char * terminal, Ast * with)
{
  AstTerminal * t = ast_terminal (n);
  if (t && !strcmp (t->start, terminal)) {
    while (n && n->sym != with->sym)
      n = n->parent;
    if (n) {
      ast_set_child (n->parent, ast_child_index (n), with);
      ast_destroy (n);
      return ast_right_terminal (with);
    }
    return NULL;
  }
  if (n->child)
    for (Ast ** c = n->child; *c; c++) {
      AstTerminal * right = ast_replace (*c, terminal, with);
      if (right) {
	int after = count_lines (right->start);
	right->line += after;
	for (c++; *c; c++)
	  ast_set_line (*c, right, false);
	right->line -= after;
	return right;
      }
    }
  return NULL;
}

typedef struct {
  const char * name;
  int line, macro;
} File;

/**
 * @brief Updates file name and line number in a File struct based on preprocessor directives.
 *
 * Scans the given preprocessor string for `#line` directives and updates the file's line number and name accordingly. Increments the line number for each newline encountered outside of directives.
 *
 * @param preproc Preprocessor string to scan for line and file information.
 * @param file Pointer to the File struct to update.
 */
static void update_file_line (const char * preproc, File * file)
{
  const char * s = preproc;
  while (*s != '\0') {
    if (*s == '#') {
      const char * u = s > preproc ? s - 1 : preproc;
      while (u != preproc && !strchr ("\n\r", *u) && strchr (" \t", *u))
	u--;
      if (u == preproc || strchr ("\n\r", *u)) {
	s++;
	while (*s && strchr (" \t", *s)) s++;
	if (strchr ("0123456789", *s)) {
	  int line = atoi (s) - 1;
	  while (*s && !strchr (" \t", *s)) s++;
	  while (*s && strchr (" \t", *s)) s++;
	  if (*s == '"') {
	    file->line = line;
	    file->name = s + 1;
	  }
	}
      }
    }
    else if (*s == '\n')
      file->line++;
    s++;
  }  
}

typedef struct {
  char * s;
  int len, size;
} String;

/**
 * @brief Appends multiple C strings to a dynamic String structure.
 *
 * Concatenates all provided null-terminated strings to the end of the String's buffer,
 * reallocating memory as needed to accommodate the new content.
 *
 * @param a Pointer to a String structure to which the strings will be appended.
 *          The argument list must be terminated with a NULL pointer.
 */
void string_append (void * a, ...)
{
  va_list ap;
  va_start (ap, a);
  char * s = va_arg(ap, char *);
  while (s) {
    String * b = a;
    b->len += strlen (s);
    if (b->len >= b->size) {
      b->size = b->len + 1024;
      char * n = realloc (b->s, b->size);
      if (!b->s) n[0] = '\0';
      b->s = n;
    }
    strcat (b->s, s);
    s = va_arg(ap, char *);
  }
  va_end (ap);
}

/**
 * @brief Checks if a string contains only whitespace and determines the last newline position if the file line number exceeds the terminal's line.
 *
 * @param s Input string to check.
 * @param file File structure containing the starting line number.
 * @param t AST terminal node providing a reference line number.
 * @return Pointer to the character after the last newline if the string contains only whitespace and the file's line number is greater than the terminal's line; otherwise, NULL.
 */
static const char * only_spaces (const char * s, const File * file, const AstTerminal * t)
{
  const char * spaces = NULL;
  int line = file->line;
  while (*s != '\0') {
    if (!strchr (" \t\n\r", *s))
      return NULL;
    if (strchr ("\n\r", *s)) {
      spaces = s + 1;
      line++;
    }
    s++;
  }
  return line > t->line ? spaces : NULL;
}

/**
 * @brief Recursively prints an AST subtree using a custom output function.
 *
 * Traverses the AST rooted at `n`, formatting and emitting its textual representation via the provided `output` callback. Handles terminals, macros, symbol annotations, file and line directives, and special cases such as macro definitions and array accesses. Updates the `file` structure to track current file and line state during output.
 *
 * @param n The AST node to print.
 * @param sym Controls symbol and annotation output modes.
 * @param real If nonzero, annotates real number nodes.
 * @param file Pointer to a File struct tracking current file and line for line directives.
 * @param output Callback function to receive formatted output fragments.
 * @param data Opaque pointer passed to the output callback.
 */
static void str_print_internal (const Ast * n, int sym, int real, File * file,
				void (* output) (void *, ...), void * data)
{
  if (n == ast_placeholder) {
    output (data, "$", NULL);
    return;
  }

  //  ast_print_file_line (n, stderr);
  AstTerminal * t = ast_terminal (n);
  if (t) {
    if (t->before) {
      const char * spaces = only_spaces (t->before, file, t);
      if (spaces)
	output (data, spaces, NULL);
      else {
	output (data, t->before, NULL);
	update_file_line (t->before, file);
      }
    }
    if (t->file) {
      int len = strlen (t->file);
      if (!file->name || strncmp (t->file, file->name, len) ||
	  (file->name[len] != '"' && file->name[len] != '\0')) {
	char line[20]; snprintf (line, 19, "%d", t->line);
	if (sym == 3)
	  output (data, "\n", t->file, ":", line, " ", NULL);
	else
	  output (data, "\n#line ", line, " \"", t->file, "\"\n", NULL);
	file->line = t->line;
	file->name = t->file;      
      }
      else if (t->line != file->line) {
	if (file->line > t->line || t->line - file->line > 10) {
	  char line[20]; snprintf (line, 19, "%d", t->line);
	  if (sym == 3)
	    output (data, "\n", t->file, ":", line, " ", NULL);
	  else
	    output (data, "\n#line ", line, "\n", NULL);
	  file->line = t->line;
	}
	else
	  while (file->line < t->line) {
	    output (data, "\n", NULL);
	    file->line++;
	  }
      }
    }
    if (sym == 1)
      output (data, "|", symbol_name (n->sym), "|", NULL);
    else if (sym == 2) {
      if (t->value)
	output (data, "^", NULL);
    }
    if (real && (n->sym == sym_DOUBLE || n->sym == sym_FLOAT))
      output (data, "real", NULL);
    else if (n->sym == sym_MACRODEF) {
      if (ast_schema (ast_ancestor (n, 2), sym_declaration_specifiers,
		      1, sym_declaration_specifiers))
	output (data, "", NULL);
      else
	output (data, "void", NULL);
    }
    else if (n->sym == sym_AUTO) {
      if (ast_schema (ast_ancestor (n, 2), sym_declaration_specifiers,
		      1, sym_declaration_specifiers,
		      0, sym_storage_class_specifier,
		      0, sym_MACRODEF) ||
	  ast_schema (ast_ancestor (n, 3), sym_parameter_declaration))
	output (data, "    ", NULL);
      else
	output (data, t->start, NULL);
    }
    else if (n->sym == sym_ELLIPSIS_MACRO)
      output (data, "{}", NULL);
    else if (n->sym == sym_IDENTIFIER) {
      if (ast_find (ast_schema (ast_ancestor (n, 5), sym_function_declaration,
				0, sym_declaration_specifiers), sym_MACRODEF)) {
	char s[20]; snprintf (s, 19, "%d", file->macro++);
	output (data, "macro", s, "_", NULL);
      }
      output (data, t->start, NULL);
    }
    else
      output (data, t->start, NULL);
    file->line += count_lines (t->start);
    if (sym == 1)
      output (data, "/", NULL);
    if (t->after) {
      output (data, t->after, NULL);
      file->line += count_lines (t->after);
    }
  }
  else { // !terminal

    /**
    Do not output macro definitions. */

    if (ast_is_macro_declaration (ast_schema (n, sym_function_definition,
					      0, sym_function_declaration))) {
      AstTerminal * t = ast_left_terminal (n);
      if (t->before) {
	const char * spaces = only_spaces (t->before, file, t);
	if (spaces)
	  output (data, spaces, NULL);
	else {
	  output (data, t->before, NULL);
	  update_file_line (t->before, file);
	}
      }
      return;
    }
    
    /**
    Ignore 'break =' macro parameters. */
    
    if (ast_schema (ast_child (n, sym_parameter_declaration), sym_parameter_declaration,
		    0, sym_BREAK)) {
      Ast * list = ast_child (n, sym_parameter_list);
      if (list)
	str_print_internal (list, sym, real, file, output, data);
      return;
    }

    AstRoot * r = ast_root (n);
    if (r && r->before) {
      output (data, r->before, NULL);
      update_file_line (r->before, file);
    }

    /** 
    Ignore dimensioned constant expressions e.g. `1.23
    [1,-1]`, `(2.*3 + 4.)[0,1]` etc. */
    
    if (n->sym == sym_array_access && n->child[2] &&
	ast_evaluate_constant_expression (n->child[0]) < DBL_MAX)
      str_print_internal (n->child[0], sym, real, file, output, data);

    /**
    Ignore the values of optional parameters. */

    else if (ast_schema (n, sym_parameter_declaration,
			 3, sym_initializer))
      for (int i = 0; i < 2; i++)
	str_print_internal (n->child[i], sym, real, file, output, data);
    
    else
      for (Ast ** c = n->child; *c; c++)
	str_print_internal (*c, sym, real, file, output, data);
    
    if (r && r->after) {
      output (data, r->after, NULL);
      file->line += count_lines (r->after);
    }
  }
}

/**
 * @brief Converts an AST subtree to a string representation.
 *
 * Generates a string representation of the AST subtree rooted at `n`, optionally including symbol and real number annotations. If `s` is provided, appends the result to `s` and returns it; otherwise, returns a newly allocated string containing the result.
 *
 * @param n The root of the AST subtree to convert.
 * @param s Optional string to append the result to. If NULL, a new string is allocated.
 * @param sym If nonzero, includes symbol names in the output.
 * @param real If nonzero, includes real number annotations in the output.
 * @return char* The resulting string, either appended to `s` or newly allocated.
 */
char * ast_str_print (const Ast * n, char * s, int sym, int real)
{
  String a = {0};
  str_print_internal (n, sym, real, &(File){0}, string_append, &a);
  if (s) {
    str_append (s, a.s);
    free (a.s);
    return s;
  }
  return a.s;
}

/**
 * @brief Writes multiple strings to a file stream.
 *
 * Appends each provided null-terminated string to the given file stream in order. The argument list must be terminated with a NULL pointer.
 *
 * @param a Pointer to a FILE stream where the strings will be written.
 */
static void file_append (void * a, ...)
{
  va_list ap;
  va_start (ap, a);
  char * s = va_arg(ap, char *);
  while (s) {
    fputs (s, a);
    s = va_arg(ap, char *);
  }
  va_end (ap);
}

/**
 * @brief Prints the textual representation of an AST subtree to a file.
 *
 * Outputs the AST rooted at @p n to the specified file stream @p fp. If @p sym is nonzero, symbol names are included in the output.
 *
 * @param n Root of the AST subtree to print.
 * @param fp Output file stream.
 * @param sym If nonzero, include symbol names in the output.
 */
void ast_print (const Ast * n, FILE * fp, int sym)
{
  str_print_internal (n, sym, 0, &(File){0}, file_append, fp);
}

/**
 * @brief Appends the textual content of an AST subtree to a string.
 *
 * Traverses the AST rooted at @p n and appends the concatenated terminal and root strings to @p s in order.
 *
 * @param n The root of the AST subtree to convert to text.
 * @param s The destination string to append to.
 * @return The resulting string with the AST's textual content appended.
 */
char * ast_str_append (const Ast * n, char * s)
{
  AstTerminal * t = ast_terminal (n);
  if (t) {
    if (t->before)
      str_append (s, t->before);
    str_append (s, t->start);
    if (t->after)
      str_append (s, t->after);
  }
  else {
    AstRoot * r = ast_root (n);
    if (r && r->before)
      str_append (s, r->before);
    for (Ast ** c = n->child; *c; c++)
      s = ast_str_append (*c, s);
    if (r && r->after)
      str_append (s, r->after);
  }
  return s;
}

/**
 * @brief Prints the file name, line number, and identifier or symbol name of an AST node to a file.
 *
 * Outputs the file and line information from the leftmost terminal of the AST node `n` to the specified file stream `fp`, followed by either the terminal's start string or the node's symbol name.
 */
void ast_print_file_line (Ast * n, FILE * fp)
{
  assert (n);
  AstTerminal * t = ast_left_terminal (n);
  assert (t);
  fprintf (fp, "%s:%d: %s\n", t->file, t->line,
	   ast_terminal(n) ? t->start : symbol_name (n->sym));
}

/**
 * @brief Prints a child AST node in a tree-like structure with indentation and branch characters.
 *
 * Prints the given AST node `n` as a child in a tree diagram to the specified file stream, using the provided indentation and branch formatting. Adjusts indentation and branch symbols based on whether the node is the last child, and delegates further printing to `ast_print_tree`.
 *
 * @param n The AST node to print as a child.
 * @param fp The file stream to print to.
 * @param indent The current indentation string.
 * @param isLast Indicates if this node is the last child at its level.
 * @param compress If true, compresses linear chains in the output.
 * @param maxdepth Maximum depth to print; deeper nodes are omitted.
 */
static void print_child_tree (Ast * n, FILE * fp,
			      const char * indent, bool isLast,
			      bool compress, int maxdepth)
{
  char * ind;
  if (indent) {
    fputs (indent, fp);
    ind = strdup (indent);
  }
  else
    ind = strdup ("");
    
  if (isLast) {
    fputs ("└─", fp);
    str_append (ind, "  ");
  }
  else {
    fputs ("├─", fp);
    str_append (ind, "│ ");
  }
  ast_print_tree (n, fp, ind, compress, maxdepth);
  free (ind);
}

/**
 * @brief Prints the structure of an AST subtree in a tree-like format.
 *
 * Recursively prints the AST rooted at node `n` to the given file stream, using indentation to represent tree hierarchy. Optionally compresses linear chains of nodes and limits the printed depth.
 *
 * @param n Root of the AST subtree to print.
 * @param fp Output file stream.
 * @param indent Indentation string for tree levels.
 * @param compress If true, compresses linear chains of single-child nodes.
 * @param maxdepth Maximum depth to print; printing stops with "..." if exceeded.
 */
void ast_print_tree (Ast * n, FILE * fp, const char * indent,
		     bool compress, int maxdepth)
{
  if (n == ast_placeholder) {
    fputs ("_placeholder_\n", fp);
    return;
  }
  if (!maxdepth) {
    fputs ("...\n", fp);
    return;
  }
  if (compress)
    while (n->child && !n->child[1] && n->child[0] != ast_placeholder)
      n = n->child[0];
  fprintf (fp, "%s", symbol_name (n->sym));
  AstTerminal * t = ast_terminal (n);
  if (t)
    fprintf (fp, " %s %s:%d\n", t->start, t->file, t->line);
  else {
    fputc ('\n', fp);
    for (Ast **c = n->child; *c; c++)
      print_child_tree (*c, fp, indent, *(c + 1) == NULL,
			compress, maxdepth - 1);
  }
}

/**
 * @brief Prints C constructor macros to recreate the AST subtree.
 *
 * Outputs C code to the given file stream that, when compiled and executed, would reconstruct the AST subtree rooted at the specified node using constructor macros. Indentation is applied for readability.
 *
 * @param n Root of the AST subtree to print.
 * @param fp Output file stream.
 * @param indent Optional indentation string for formatting.
 */
void ast_print_constructor (Ast * n, FILE * fp, const char * indent)
{
  if (indent)
    fputs (indent, fp);
  AstTerminal * t = ast_terminal (n);
  if (t) {
    if (t->start[1] == '\0')
      fprintf (fp, "NCA(n, \"%s\")", t->start);
    else
      fprintf (fp, "NA(n, sym_%s, \"%s\")", symbol_name (n->sym), t->start);
  }
  else {
    fprintf (fp, "NN(n, sym_%s,\n", symbol_name (n->sym));
    char * ind = indent ? strdup (indent) : strdup ("");
    str_append (ind, "   ");
    for (Ast **c = n->child; *c; c++) {
      ast_print_constructor (*c, fp, ind);
      if (*(c + 1))
	fputs (",\n", fp);
    }
    fputs (")", fp);
    free (ind);
  }
}

/**
 * @brief Creates a new terminal AST node with the specified symbol and source text.
 *
 * The new terminal node inherits file and line information from the leftmost terminal of the parent node.
 *
 * @param parent The parent AST node to which the terminal will be attached.
 * @param symbol The symbol representing the terminal's token type.
 * @param start The source text associated with the terminal.
 * @return Pointer to the newly created AstTerminal node.
 */
AstTerminal * ast_terminal_new (Ast * parent, int symbol, const char * start)
{
  AstTerminal * t = allocate (ast_get_root (parent)->alloc,
			      sizeof (AstTerminal));
  memset (t, 0, sizeof (AstTerminal));
  ((Ast *)t)->sym = symbol;
  ((Ast *)t)->parent = parent;
  t->start = strdup (start);
  AstTerminal * r = ast_left_terminal (parent);
  t->file = r->file;
  t->line = r->line;
  return t;
}

/**
 * @brief Creates a nested chain of AST nodes from a variadic list of symbols.
 *
 * Constructs a sequence of AST nodes, each as the sole child of the previous node, using symbols provided via a variadic argument list. The chain starts from the given parent node and terminates when a negative symbol is encountered.
 *
 * @param parent The parent AST node to which the chain is attached.
 * @param ap Variadic argument list of integer symbols, ending with a negative value.
 * @return Ast* Pointer to the root of the newly created AST chain, or NULL if the first symbol is negative.
 */
static Ast * vast_new_internal (Ast * parent, va_list ap)
{
  int sym = va_arg (ap, int);
  if (sym < 0)
    return NULL;
  AstRoot * root = ast_get_root (parent);
  assert (root->alloc);
  Ast * n = allocate (root->alloc, sizeof (Ast));
  n->sym = sym;
  n->parent = parent;
  n->child = NULL;
  Ast * m = n;
  sym = va_arg (ap, int);
  while (sym >= 0) {
    Ast * c = allocate (root->alloc, sizeof (Ast));
    c->sym = sym;
    c->parent = m;
    m->child = allocate (root->alloc, 2*sizeof (Ast *));
    m->child[0] = c;
    m->child[1] = NULL;
    m = c;
    sym = va_arg (ap, int);
  }
  return n;
}

/**
 * @brief Creates a nested chain of AST nodes as children of a parent, using a variadic list of symbols.
 *
 * Constructs a sequence of AST nodes where each symbol in the argument list becomes a child of the previous node, starting from the given parent.
 *
 * @return Pointer to the deepest (last) AST node created in the chain.
 */
Ast * ast_new_internal (Ast * parent, ...)
{
  va_list ap;
  va_start (ap, parent);
  Ast * n = vast_new_internal (parent, ap);
  va_end (ap);
  return n;
}

/**
 * @brief Checks if an AST node matches a variadic schema of symbols and child indices.
 *
 * Traverses the AST starting from node `n`, following a sequence of (symbol, child index) pairs provided via a variadic argument list. At each step, verifies that the current node matches the expected symbol and that the specified child exists. Returns the final matched node if the schema is satisfied, or NULL otherwise.
 *
 * @param n The starting AST node.
 * @param ap Variadic argument list specifying alternating symbol and child index pairs, ending with a negative child index.
 * @return Ast* The matched AST node if the schema is satisfied, or NULL if not.
 */
static Ast * vast_schema_internal (const Ast * n, va_list ap)
{
  int sym = va_arg(ap, int);
  if (n == ast_placeholder || n->sym != sym)
    return NULL;
  int c = va_arg(ap, int);
  while (c >= 0) {
    if (!n->child)
      return NULL;
    int i;
    for (i = 0; i <= c && n->child[i]; i++);
    if (i <= c)
      return NULL;
    n = n->child[c];
    int sym = va_arg(ap, int);
    if (n == ast_placeholder || n->sym != sym)
      return NULL;
    c = va_arg(ap, int);
  }
  return (Ast *) n;
}

/**
 * @brief Checks if an AST node matches a specified schema.
 *
 * Traverses the AST node `n` according to a variadic schema of symbol and child index pairs.
 * Returns the node matching the schema, or NULL if no match is found.
 *
 * @param n The root AST node to check.
 * @return Ast* The matched AST node, or NULL if the schema does not match.
 */
Ast * ast_schema_internal (const Ast * n, ...)
{
  if (!n)
    return NULL;
  va_list ap;
  va_start (ap, n);
  n = vast_schema_internal (n, ap);
  va_end (ap);
  return (Ast *) n;
}

/**
 * @brief Recursively searches an AST subtree for a node matching a schema and optional identifier.
 *
 * Traverses the AST rooted at @p n, looking for a node that matches the schema specified by the variadic argument list. If @p identifier is provided, the node's terminal string must also match it. Returns the first matching node found, or NULL if none is found.
 *
 * @param n The root of the AST subtree to search.
 * @param identifier Optional terminal string to match; if NULL, only schema is matched.
 * @param ap Variadic argument list specifying the schema to match.
 * @return Ast* Pointer to the matching AST node, or NULL if not found.
 */
static Ast * vast_find_internal (const Ast * n, const char * identifier, va_list ap)
{
  if (n == ast_placeholder)
    return NULL;
  va_list bp;
  va_copy (bp, ap);
  Ast * found = vast_schema_internal (n, bp);
  va_end (bp);
  if (found && (!identifier ||
		(ast_terminal (found) && !strcmp (ast_terminal (found)->start, identifier))))
    return found;
  if (!ast_terminal(n))
    for (Ast ** c = n->child; *c; c++)
      if ((found = vast_find_internal (*c, identifier, ap)))
	return found;
  return NULL;
}

/**
 * @brief Searches an AST subtree for a node matching a schema and optional identifier.
 *
 * Recursively searches the AST subtree rooted at `n` for a node that matches a variadic schema of symbol/index pairs and, if provided, a terminal node with the specified identifier string.
 *
 * @param n Root of the AST subtree to search.
 * @param identifier Optional identifier string to match in terminal nodes. If NULL, only the schema is matched.
 * @return Ast* Pointer to the first matching AST node, or NULL if no match is found.
 */
Ast * ast_find_internal (const Ast * n, const char * identifier, ...)
{
  if (!n)
    return NULL;
  va_list ap;
  va_start (ap, identifier);
  n = vast_find_internal (n, identifier, ap);
  va_end (ap);
  return (Ast *) n;
}

/**
 * @brief Copies basic AST node fields from a source node to a destination node.
 *
 * Sets the symbol of the destination node to match the source, clears its child pointer, and sets its parent pointer to itself.
 *
 * @param dst Destination AST node to populate.
 * @param src Source AST node to copy from.
 * @return Ast* The destination AST node after copying.
 */
static Ast * copy_ast (Ast * dst, const Ast * src)
{
  dst->sym = src->sym;
  dst->child = NULL;
  dst->parent = dst;
  return dst;
}

/**
 * @brief Creates a deep copy of an AST terminal node, duplicating its strings and metadata.
 *
 * Copies the terminal node's string fields (`before`, `start`, `after`) and file/line information into a newly allocated node associated with the destination AST root.
 *
 * @param src The source AST terminal node to copy.
 * @param dst_root The destination AST root for memory allocation.
 * @param src_root The source AST root (unused in this function, but may be relevant for context).
 * @return Ast* Pointer to the newly copied AST terminal node.
 */
static Ast * terminal_copy (const AstTerminal * src,
			    const AstRoot * dst_root, const AstRoot * src_root)
{
  AstTerminal * dst =
    (AstTerminal *) copy_ast (allocate (dst_root->alloc, sizeof (AstTerminal)),
			      (Ast *)src);
  dst->before = src->before ? strdup (src->before) : NULL;
  dst->start = strdup (src->start);
  dst->after = src->after ? strdup (src->after) : NULL;
  dst->line = src->line;
  dst->file = src->file;
  return (Ast *)dst;
}

/**
 * @brief Creates a deep copy of an AST root node, including its metadata and allocation context.
 *
 * Allocates a new AST root node, copies its basic fields from the source, duplicates the `before` and `after` strings if present, and initializes a new allocator. The copied root's stack and parent pointers are set to NULL.
 *
 * @param src The source AST root node to copy.
 * @return Ast* Pointer to the newly allocated AST root node.
 */
static Ast * root_copy (const AstRoot * src)
{
  Allocator * alloc = new_allocator();
  AstRoot * dst = (AstRoot *) copy_ast (allocate (alloc, sizeof (AstRoot)),
					(Ast *)src);
  dst->alloc = alloc;
  dst->before = src->before ? strdup (src->before) : NULL;
  dst->after = src->after ? strdup (src->after) : NULL;
  dst->stack = NULL;
  ((Ast *)dst)->parent = NULL;
  return (Ast *)dst;
}

/**
 * @brief Creates a shallow copy of a single AST node, duplicating terminal or root data as needed.
 *
 * If the node is a root, a new root is created and pointers to the source and destination roots are updated.
 * If the node is a terminal, a new terminal node is created with duplicated string and file/line information.
 * For non-terminal, non-root nodes, a shallow copy is made and an empty child array is allocated.
 *
 * @param n The AST node to copy.
 * @param dst_root Pointer to receive the destination root if a new root is created.
 * @param src_root Pointer to receive the source root if applicable.
 * @return Ast* The copied AST node.
 */
Ast * ast_copy_single (const Ast * n,
		       AstRoot ** dst_root, AstRoot ** src_root)
{
  Ast * c = NULL;
  AstRoot * r = ast_root (n);
  if (r) {
    c = root_copy (r);
    *src_root = r;
    *dst_root = ast_root (c);
  }
  AstTerminal * t = ast_terminal (n);
  if (t)
    c = terminal_copy (t, *dst_root, *src_root);
  else if (!c)
    c = copy_ast (allocate ((*dst_root)->alloc, sizeof (Ast)), n);
  if (!t) {
    int len = 0;
    for (Ast ** i = n->child; *i; i++, len++);
    c->child = allocate ((*dst_root)->alloc, (len + 1)*sizeof (Ast *));
    c->child[len] = NULL;
  }
  return c;
}

/**
 * @brief Recursively copies an AST subtree, marking if a schema match is found.
 *
 * Copies the AST node `n` and its descendants, using a variadic schema to determine if any node matches the specified pattern. Sets `*found` to true if a matching node is encountered during the copy. Child nodes that are placeholders are preserved as placeholders in the copy.
 *
 * @param n The AST node to copy.
 * @param ap Variadic argument list specifying the schema to match.
 * @param found Pointer to a boolean set to true if a schema match is found.
 * @param dst_root The destination AST root for the copy.
 * @param src_root The source AST root.
 * @return Ast* The root of the copied AST subtree.
 */
static Ast * vast_copy_internal (const Ast * n, va_list ap, bool * found,
				 AstRoot * dst_root, AstRoot * src_root)
{
  Ast * c = ast_copy_single (n, &dst_root, &src_root);
  va_list cp;
  va_copy (cp, ap);
  if (vast_schema_internal (n, cp))
    *found = true;
  va_end (cp);

  if (!ast_terminal (n)) {
    int len = 0;
    for (Ast ** i = n->child, ** j = c->child;
	 *i && !(*found); i++, j++, len++) {
      if (*i == ast_placeholder)
	*j = ast_placeholder;
      else {
	*j = vast_copy_internal (*i, ap, found, dst_root, src_root);
	(*j)->parent = c;
      }
    }
    c->child[len] = NULL;
  }
  return c;
}

/**
 * @brief Creates a deep copy of an AST subtree matching a specified schema.
 *
 * Recursively copies the AST subtree rooted at `n` that matches a variadic schema of symbol and child index pairs. The copied subtree's parent pointer is set to the original node.
 *
 * @param n The root of the AST subtree to copy.
 * @return Ast* Pointer to the root of the copied AST subtree, or NULL if `n` is NULL.
 */
Ast * ast_copy_internal (const Ast * n, ...)
{
  if (!n) return NULL;
  va_list ap;
  va_start (ap, n);
  bool found = false;
  AstRoot * src_root = ast_get_root (n);
  Ast * c = vast_copy_internal (n, ap, &found, src_root, src_root);
  c->parent = (Ast *) n;
  va_end (ap);
  return c;
}

/**
 * @brief Parses a C expression string into an AST node.
 *
 * Wraps the given expression in a dummy function, parses it, extracts the corresponding statement or declaration node, removes file and line information, and returns the resulting AST subtree.
 *
 * @param expr The C expression to parse.
 * @return Ast* The AST node representing the parsed expression, or NULL on failure.
 */
Ast * ast_parse_expression (const char * expr, AstRoot * parent)
{
  char * s = NULL;
  str_append (s, "void main() {", expr, "}");
  Ast * n = (Ast *) ast_parse (s, parent);
  free (s);
  if (n) {
    ast_pop_scope (parent->stack, n);
    Ast * c = ast_find (n, sym_statement);
    if (c)
      c = c->child[0];
    else
      c = ast_find (n, sym_declaration); 
    cancel_file_line (c);
    Ast ** i;
    for (i = c->parent->child; *i && *i != c; i++);
    *i = ast_placeholder;
    ast_destroy (n);
    n = c;
  }
  return n;
}

/**
 * @brief Parses a C external declaration string into an AST node.
 *
 * Parses the given external declaration string, extracts the corresponding external declaration AST node, removes file and line information, detaches it from its temporary parent, and returns it. Returns NULL if parsing fails or the external declaration is not found.
 *
 * @param decl The C external declaration as a string.
 * @param parent The AST root to associate with the parsed node.
 * @return Ast* The extracted external declaration AST node, or NULL on failure.
 */
Ast * ast_parse_external_declaration (const char * decl, AstRoot * parent)
{
  Ast * def = (Ast *) ast_parse (decl, parent);
  Ast * n = ast_find (def, sym_external_declaration);
  if (!n) {
    ast_destroy (def);
    return NULL;
  }
  cancel_file_line (n);
  Ast ** i;
  for (i = n->parent->child; *i && *i != n; i++);
  assert (*i == n);
  *i = ast_placeholder;
  ast_destroy (def);
  return n;
}

/**
 * @brief Parses the contents of a file into an AST root node.
 *
 * Reads the entire file stream, parses its contents as source code, and returns the resulting AST root node. The returned AST root is associated with the given parent, if provided.
 *
 * @param fp Input file stream to parse.
 * @param parent Optional parent AST root for context.
 * @return AstRoot* The parsed AST root node, or NULL on failure.
 */
AstRoot * ast_parse_file (FILE * fp, AstRoot * parent)
{
  char * buffer = NULL;
  size_t len = 0, maxlen = 0;
  int c;
  while ((c = fgetc (fp)) != EOF) {
    if (len >= maxlen) {
      maxlen += 4096;
      buffer = realloc (buffer, maxlen);      
    }
    buffer[len++] = c;
  }
  if (len >= maxlen) {
    maxlen++;
    buffer = realloc (buffer, maxlen);      
  }
  buffer[len++] = '\0';
  AstRoot * root = ast_parse (buffer, parent);
  free (buffer);
  return root;
}

/**
 * @brief Returns the identifier string with known typedef prefixes removed.
 *
 * Removes the prefixes "face ", "vertex ", or "symmetric " from the beginning of the identifier string if present.
 *
 * @param identifier The input identifier string.
 * @return const char* Pointer to the identifier string after any recognized prefix, or the original string if no prefix is found.
 */

static const char * ignore_prefixes (const char * identifier)
{
  if (!strncmp (identifier, "face ", 5))
    return identifier + 5;
  else if (!strncmp (identifier, "vertex ", 7))
    return identifier + 7;
  else if (!strncmp (identifier, "symmetric ", 10))
    return identifier + 10;
  return identifier;
}

/**
 * @brief Searches a stack for an identifier declaration between two AST nodes.
 *
 * Traverses the stack from the node after `start` up to (but not including) `end`, returning the first AST node representing a declaration of the specified identifier. Ignores known prefixes in the identifier string. Handles both standard and lexer-specific identifier representations.
 *
 * @param identifier The identifier string to search for (prefixes ignored).
 * @param start The AST node to start searching after (inclusive if NULL).
 * @param end The AST node at which to stop searching (exclusive; may be NULL).
 * @return Pointer to the AST node declaring the identifier, or NULL if not found.
 */
Ast * ast_identifier_declaration_from_to (Stack * stack, const char * identifier,
					  const Ast * start, const Ast * end)
{
  if (!identifier) return NULL;
  identifier = ignore_prefixes (identifier);
  
  Ast ** d;
  int i = 0;
  if (start) {
    for (; (d = stack_index (stack, i)) && *d != start; i++);
    if (d) {
      assert (*d == start);
      i++;
    }
    else {
      assert (false); // start not found!!
      return NULL;
    }
  }
  for (; (d = stack_index (stack, i)); i++)
    if (end && *d == end)
      break;
    else if (*d && (*d)->sym == sym_IDENTIFIER) {

      /**
      WARNING: this assumes that the "after" string is never modified
      for the declaration identifiers stored in the stack. */

      if (!ast_terminal(*d)->after) {
	if (ast_terminal(*d)->start &&
	    !strcmp (ast_terminal(*d)->start, identifier))
	  return *d;
      }

      /**
      If 'after' is defined, we assume that the function is called
      from the [lexer](tokens.lex) and that 'after' is the last
      character of the declaration identifier (as set by the
      lexer). */
      
      else {
	char * s = ast_terminal(*d)->start, * end = ast_terminal(*d)->after;
	const char * i = identifier;
	for (; *i != '\0' && s <= end && *s == *i; s++, i++);
	if (*i == '\0' && s == end + 1)
	  return *d;
      }
    }
  return NULL;
}

/**
 * @brief Searches a stack for an AST node declaring the specified identifier.
 *
 * Looks up an identifier declaration in the given stack, first attempting an exact match, then retrying after removing known prefixes ("face ", "vertex ", "symmetric ") from the identifier.
 *
 * @param identifier The identifier string to search for, with or without known prefixes.
 * @return Ast* Pointer to the AST node declaring the identifier, or NULL if not found.
 */
Ast * ast_identifier_declaration (Stack * stack, const char * identifier)
{
  Ast ** n = fast_stack_find (stack, identifier);
  if (!n)
    n = fast_stack_find (stack, ignore_prefixes (identifier));
  return n ? *n : NULL;
}

/**
 * @brief Appends multiple strings to a dynamically allocated string, reallocating as needed.
 *
 * Concatenates all provided strings to the end of `src`, reallocating memory to fit the result.
 * The argument list must be terminated with a NULL pointer.
 *
 * @param src The original string to append to, or NULL to start a new string.
 * @return char* Pointer to the newly allocated string containing the concatenation, or NULL if the result is empty.
 */
char * str_append_realloc (char * src, ...)
{
  va_list ap;
  va_start (ap, src);
  char * i = va_arg(ap, char *);
  int len = src ? strlen(src) : 0;
  while (i) {
    len += strlen (i);
    i = va_arg(ap, char *);
  }
  va_end (ap);

  if (len < 1)
    return NULL;
  
  char * dst = realloc (src, len + 1);
  if (!src)
    dst[0] = '\0';
  
  va_start (ap, src);
  i = va_arg(ap, char *);
  while (i) {
    strcat (dst, i);
    i = va_arg(ap, char *);
  }
  va_end (ap);

  return dst;
}

/**
 * @brief Prepends one or more strings to a source string, reallocating memory as needed.
 *
 * Concatenates all provided strings in order before the original `src` string, allocating a new buffer for the result. The original `src` is freed. The argument list must be terminated with a NULL pointer.
 *
 * @param src The original string to prepend to. May be NULL.
 * @return char* Newly allocated string with all prepended content, or NULL if the result is empty.
 */
char * str_prepend_realloc (char * src, ...)
{
  va_list ap;
  va_start (ap, src);
  char * i = va_arg(ap, char *);
  int len = src ? strlen(src) : 0;
  while (i) {
    len += strlen (i);
    i = va_arg(ap, char *);
  }
  va_end (ap);

  if (len < 1)
    return NULL;
  
  char * dst = malloc (len + 1);
  dst[0] = '\0';

  va_start (ap, src);
  i = va_arg(ap, char *);
  while (i) {
    strcat (dst, i);
    i = va_arg(ap, char *);
  }
  va_end (ap);

  if (src)
    strcat (dst, src);

  free (src);
  return dst;
}

/**
 * @brief Prints the file, line, pointer, and identifier string of an AST terminal node.
 *
 * Outputs the source file name, line number, memory address, and the identifier string associated with the given AST node to the specified file stream.
 */
void ast_identifier_print (Ast * identifier, FILE * fp)
{
  AstTerminal * t = ast_terminal (identifier);
  fprintf (fp, "%s:%d: %p ", t->file, t->line, t);
  if (!t->after)
    fprintf (fp, "* %s", t->start);
  else {
    char * s = t->start, * end = t->after;
    for (; s <= end; s++)
      fputc (*s, fp);
  }
  fputc ('\n', fp);
}

/**
 * @brief Prints the contents of an AST node stack to a file.
 *
 * Each stack entry is printed as its identifier (if terminal), as "_placeholder_" for placeholders, or as the symbol name and pointer for other nodes.
 *
 * @param stack Stack of AST node pointers to print.
 * @param fp Output file stream.
 */
void ast_stack_print (Stack * stack, FILE * fp)
{
  Ast ** n;
  for (int i = 0; (n = stack_index (stack, i)); i++)
    if (*n) {
      if (ast_terminal (*n))
	ast_identifier_print (*n, fp);
      else if (*n == ast_placeholder)
	fprintf (fp, "_placeholder_\n");
      else
	fprintf (fp, "%s %p\n", symbol_name ((*n)->sym), *n);
    }
}

/**
 * @brief Sets the symbol and start character of an AST node.
 *
 * Updates the AST node's symbol to correspond to the given character and sets the first character of its terminal's start string to that character.
 *
 * @param n The AST node to update.
 * @param c The character to assign as the node's symbol and terminal start.
 */
void ast_set_char (Ast * n, int c)
{
  n->sym = token_symbol(c), ast_terminal (n)->start[0] = c;
}

/**
 * @brief Recursively removes the content of an AST node and its descendants, merging terminal 'before' strings.
 *
 * For each terminal node in the subtree rooted at `n`, prepends its 'before' string to the provided `before` terminal and frees associated strings.
 */
void ast_remove_internal (Ast * n, AstTerminal * before)
{
  if (n->child) {
    for (Ast ** c = n->child; *c; c++)
      if (*c != ast_placeholder)
	ast_remove_internal (*c, before);
  }
  else {
    AstTerminal * t = ast_terminal (n);
    if (t->before) {
      str_prepend (before->before, t->before);
      free (t->before), t->before = NULL;
    }
    free (t->start), t->start = NULL;
    free (t->after), t->after = NULL;
  }
}

/**
 * @brief Removes an AST node from its parent's child list and merges its 'before' string.
 *
 * If the node has a parent, it is removed from the parent's child array by shifting subsequent children left.
 * Then, recursively removes the node's content, prepending its 'before' string to the provided terminal.
 *
 * @param n The AST node to remove.
 * @param before The terminal to which the removed node's 'before' string will be prepended.
 */
void ast_remove (Ast * n, AstTerminal * before)
{
  if (n->parent) {
    Ast ** c;
    for (c = n->parent->child; *c && *c != n; c++);
    if (*c == n)
      for (; *c; c++)
	*c = *(c + 1);
  }
  ast_remove_internal (n, before);
}

/**
 * @brief Returns the next terminal node in an in-order traversal of the AST.
 *
 * Searches for the next terminal node following the given node `n` within the AST, skipping placeholder nodes. Returns NULL if there is no subsequent terminal.
 *
 * @param n The current AST node.
 * @return AstTerminal* Pointer to the next terminal node, or NULL if none exists.
 */
AstTerminal * ast_next_terminal (const Ast * n)
{
  if (!n || n == ast_placeholder)
    return NULL;
  Ast * parent = n->parent;
  while (parent) {
    int index = ast_child_index (n) + 1;
    if (index <= 0)
      return NULL;
    while ((n = parent->child[index])) {
      if (n != ast_placeholder)
	return ast_left_terminal (n);
      index++;
    }
    n = parent;
    parent = n->parent;
  }
  return NULL;
}

/**
 * @brief Removes the content of an AST node and destroys it.
 *
 * Recursively removes the content of the AST node `n`, merging its `before` string with the next terminal node if available, and then frees all associated memory.
 */
void ast_erase (Ast * n)
{
  AstTerminal * t = ast_next_terminal (n);
  if (t)
    ast_remove_internal (n, t);
  ast_destroy (n);
}

/**
 * @brief Recursively verifies that each child's parent pointer references its actual parent node.
 *
 * Asserts that for every child of the given AST node, the child's parent pointer is set to the current node, and recursively checks all descendants.
 */
static void ast_check_children (Ast * n)
{
  if (n->child)
    for (Ast ** c = n->child; *c; c++) {
      assert ((*c)->parent == n);
      ast_check_children (*c);
    }
}

/**
 * @brief Verifies the integrity of the AST subtree rooted at the given node.
 *
 * Recursively checks that each child's parent pointer is correct and that the node is properly linked in its parent's child array, traversing up to the root. Asserts are triggered if inconsistencies are found.
 */
void ast_check (Ast * n)
{
  ast_check_children (n);
  while (n->parent) {
    assert (n->parent->child);
    Ast ** c;
    for (c = n->parent->child; *c && *c != n; c++);
    assert (*c == n);
    assert (n->parent != n);
    n = n->parent;
  }
}

/**
 * @brief Recursively sets file and line information for an AST node and its descendants.
 *
 * If the node or its terminals lack file information, or if overwrite is true, assigns the file and line from the provided terminal to each terminal node in the subtree.
 *
 * @param n The root AST node whose subtree will be updated.
 * @param l The terminal node providing the file and line information.
 * @param overwrite If true, existing file and line information will be overwritten.
 */
void ast_set_line (Ast * n, AstTerminal * l, bool overwrite)
{
  if (n == ast_placeholder)
    return;
  AstTerminal * t = ast_terminal (n);
  if (t && (overwrite || !t->file)) {
    t->file = l->file;
    t->line = l->line;
  }
  if (n->child)
    for (Ast ** c = n->child; *c; c++)
      ast_set_line (*c, l, overwrite);
}

/**
 * @brief Normalizes whitespace and file/line information in an AST subtree.
 *
 * Recursively replaces the `before` and `after` strings of all terminal nodes in the subtree rooted at `n` with a single space and NULL, respectively, and sets their file and line information to match those of terminal `t`.
 *
 * @param n The root of the AST subtree to flatten.
 * @param t The terminal node whose file and line information will be propagated.
 * @return The root of the flattened AST subtree.
 */
Ast * ast_flatten (Ast * n, AstTerminal * t)
{
  AstTerminal * r = ast_terminal (n);
  if (r) {
    if (r->before) {
      free (r->before);
      r->before = strdup (" ");
    }
    if (r->after) {
      free (r->after);
      r->after = NULL;
    }
    r->file = t->file;
    r->line = t->line;
  }
  if (n->child)
    for (Ast ** c = n->child; *c; c++)
      ast_flatten (*c, t);
  return n;
}

/**
 * @brief Recursively checks whether two AST subtrees are structurally and textually identical.
 *
 * Compares the structure, symbols, and terminal strings of both AST nodes and their descendants.
 *
 * @param a First AST node to compare.
 * @param b Second AST node to compare.
 * @return true if the ASTs are identical in structure and terminal content, false otherwise.
 */
bool ast_are_identical (const Ast * a, const Ast * b)
{
  if (a == NULL && b == NULL)
    return true;
  if (a == NULL || b == NULL)
    return false;
  if (a->sym != b->sym)
    return false;
  AstTerminal * ta = ast_terminal (a), * tb = ast_terminal (b);
  if (ta) {
    if (!tb || strcmp (ta->start, tb->start))
      return false;
  }
  else if (tb)
    return false;
  else {
    Ast ** c, ** d;
    for (c = a->child, d = b->child; *c && *d; c++, d++)
      if (!ast_are_identical (*c, *d))
	return false;
    if (*c || *d)
      return false;
  }
  return true;
}
