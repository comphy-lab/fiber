#include <stdio.h>
#include <stdlib.h>
#include "ast.h"
#include "symbols.h"

/**
 * @brief Performs a grammar check on the given AST node.
 *
 * This function is intended to validate the grammar of the specified AST node.
 * The actual grammar rules are included from "grammar.h". Currently, the function
 * always returns false, indicating that no grammar check passes.
 *
 * @param n Pointer to the AST node to check.
 * @return true if the node passes the grammar check; false otherwise.
 */
static bool check (Ast * n)
{
  int sym_error = sym_YYerror;
  #include "grammar.h"
  return false;
}

/**
 * @brief Checks the grammar correctness of an AST node and its subtree.
 *
 * Validates the grammar of the given AST node. If the node is non-terminal, checks its structure using grammar rules and, if `recursive` is true, recursively checks all child nodes. If `stencils` is true, allows a special recursive macro statement pattern used for stencil construction. For terminal nodes, asserts the presence of source location metadata. Aborts the program and prints the subtree if a grammatical error is detected.
 *
 * @param n The AST node to check.
 * @param recursive Whether to recursively check all child nodes.
 * @param stencils Enables acceptance of the stencil-specific macro statement pattern.
 * @return The original AST node if grammar is correct; otherwise, the program aborts.
 */

Ast * ast_check_grammar (Ast * n, bool recursive, bool stencils)
{
  if (!n)
    return NULL;
  if (n->child) {
    if (!check (n) &&
	!(stencils && n->sym == sym_macro_statement &&
	  ((n->child[0] && n->child[0]->sym == sym_MACRO &&
	    n->child[1] && n->child[1]->sym == token_symbol('(') &&
	    n->child[2] && n->child[2]->sym == token_symbol(')') &&
	    n->child[3] && n->child[3]->sym == sym_statement &&
	    n->child[4] && n->child[4]->sym == sym_macro_statement && !n->child[5]) ||
	   (n->child[0] && n->child[0]->sym == sym_MACRO &&
	    n->child[1] && n->child[1]->sym == token_symbol('(') &&
	    n->child[2] && n->child[2]->sym == sym_argument_expression_list &&
	    n->child[3] && n->child[3]->sym == token_symbol(')') &&
	    n->child[4] && n->child[4]->sym == sym_statement &&
	    n->child[5] && n->child[5]->sym == sym_macro_statement && !n->child[6])))) {
      fprintf (stderr, "\ngrammatical error:\n");
      ast_print_tree (n, stderr, 0, 0, 10);
      abort ();
    }
    if (recursive)
      for (Ast ** c = n->child; *c; c++) {
	assert ((*c)->parent == n);
	ast_check_grammar (*c, true, stencils);
      }
  }
  else {
    AstTerminal * t = ast_terminal (n);
    assert (t->file && t->line);
  }
  return n;
}
