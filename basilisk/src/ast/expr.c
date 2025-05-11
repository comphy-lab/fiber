#include "ast.h"
#include "symbols.h"

/**
 * @brief Parses input code from the command line and prints its abstract syntax tree (AST).
 *
 * Expects at least one argument containing code to parse. Attempts to parse the input as a full AST, falling back to expression parsing if necessary. On success, prints the AST in constructor form if additional arguments are provided; otherwise, prints a formatted tree representation. Prints usage or error messages to stderr and returns a nonzero exit code on failure.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return int Exit status code (0 on success, 1 on error).
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
