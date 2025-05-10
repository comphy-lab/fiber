#include "ast.h"
#include "symbols.h"

/**
 * @brief Entry point for the AST parser and printer command-line tool.
 *
 * Parses the first command-line argument as code or an expression, constructs its abstract syntax tree (AST), and prints the AST in different formats depending on the presence of additional arguments. Exits with an error code if parsing fails.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return int Returns 0 on success, or 1 if usage is incorrect or parsing fails.
 */
int main (int argc, char * argv[])
{
  if (argc < 2) {
    fprintf (stderr, "usage: %s 'code' [OPTIONS]\n", argv[0]);
    return 1;
  }
  Ast * n = (Ast *) ast_parse (argv[1], NULL);
  if (!n)
    n = ast_parse_expression (argv[1], NULL);
  if (!n) {
    fprintf (stderr, "%s: error: could not parse code\n", argv[0]);
    return 1;
  }
  if (argc > 2) {
    ast_print_constructor (n, stderr, 0);
    fputc ('\n', stderr);
  }
  else
    ast_print_tree (n, stderr, 0, false, -1);
  ast_destroy (n);
}
