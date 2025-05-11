/**
# Stencils accesses for automatic boundary conditions

This file defines the `ast_stencil()` function which returns an
[AST](README) containing code obtained by transforming the input point
function or foreach loop into a simplified version containing only
calls to functions recording stencil read/write accesses to different
fields.

This code is typically added by the [Basilisk C
translator](translate.c#stencils) before each foreach loop body to
detect which fields will be modified by the loop and which fields
require updates to boundary conditions (to ensure consistent stencil
accesses). */

#include <stdlib.h>
#include <string.h>
#include "ast.h"
#include "symbols.h"

#if 0 // for debugging
# define CHECK(x) ast_check_grammar(x, true, false)
#else
# define CHECK(x) (x)
#endif

/**
 * @brief Determines if a function AST node represents a stencil function.
 *
 * Checks whether the function's identifier starts with the "_stencil_" prefix.
 *
 * @param n The AST node representing a function.
 * @return true if the function is a stencil function, false otherwise.
 */
bool ast_is_stencil_function (Ast * n)
{
  Ast * identifier = ast_function_identifier (n);
  if (!identifier)
    return false;
  int len = strlen (ast_terminal (identifier)->start) - 9;
  return len > 0 && !strncmp (ast_terminal (identifier)->start, "_stencil_", 9);
}

/**
 * @brief Determines if an AST node represents a non-constant scalar or vertex scalar.
 *
 * Checks whether the given AST node has a type of "scalar" or "vertex scalar" and is not a constant postfix expression.
 *
 * @return true if the node is a non-constant scalar or vertex scalar; false otherwise.
 */
static bool is_scalar (Ast * n, Stack * stack)
{
  const char * typename = ast_typedef_name (ast_expression_type (n, stack, true));  
  return typename &&
    (!strcmp (typename, "scalar") || !strcmp (typename, "vertex scalar")) &&
    !ast_constant_postfix_expression (n, stack);
}

/**
 * @brief Checks if an AST node represents a point function call.
 *
 * Scans the argument list of a function call AST node and returns the argument node if it is an identifier named "point" or a function call to "neighborp". Returns NULL if no such argument is found.
 *
 * @param n The AST node to check.
 * @return Ast* The argument node corresponding to a point function call, or NULL if not found.
 */
static Ast * is_point_function_call (Ast * n)
{  
  Ast * arguments = ast_child (n, sym_argument_expression_list);
  if (!arguments)
    return NULL;
  foreach_item (arguments, 2, argument)
    if (argument != ast_placeholder) {
      Ast * identifier = ast_is_identifier_expression (argument->child[0]);
      if (identifier) {
	if (!strcmp (ast_terminal(identifier)->start, "point"))
	  return argument;
      }
      else if ((identifier = ast_function_call_identifier
		(ast_schema (ast_is_unary_expression (argument->child[0]), sym_unary_expression,
			     0, sym_postfix_expression,
			     0, sym_function_call))) &&
	       !strcmp (ast_terminal (identifier)->start, "neighborp"))
	return argument;
    }
  return NULL;
}

/**
 * @brief Determines if an AST node represents a stencil field access.
 *
 * Returns true if the node is an array access of a scalar field, a function call to a recognized stencil-related function (such as "val", "neighbor", "_assign", etc.), or a point function call. Otherwise, returns false.
 *
 * @param n The AST node to check.
 * @param stack The current symbol stack for type resolution.
 * @return true if the node is a stencil field access; false otherwise.
 */
static bool is_field_access (Ast * n, Stack * stack)
{  
  switch (n->sym) {

  case sym_array_access: {
    if (is_scalar (n->child[0], stack))
      return true;
    break;
  }

  case sym_function_call: {
    Ast * identifier = ast_function_call_identifier (n);
    AstTerminal * t;
    if (identifier && (t = ast_terminal (identifier)) &&
	(!strcmp (t->start, "val") ||
	 !strcmp (t->start, "val_diagonal") ||
	 !strcmp (t->start, "fine") ||
	 !strcmp (t->start, "coarse") ||
	 !strcmp (t->start, "neighbor") ||
	 !strcmp (t->start, "aparent") ||
	 !strcmp (t->start, "child") ||
	 !strcmp (t->start, "_assign") ||
	 !strcmp (t->start, "_overflow") ||
	 !strcmp (t->start, "r_assign")))
      return true;
    if (is_point_function_call (n))
      return true;
    break;
  }

  }

  return false;
}

/**
 * @brief Appends an AST item to a block list, creating a new list if necessary.
 *
 * If the list is empty, returns a new list containing the item. Otherwise, returns a new list with the item appended.
 *
 * @param list The AST block list to append to.
 * @param item The AST node to append.
 * @return Ast* The resulting AST block list with the item appended.
 */
static Ast * block_list_append (Ast * list, Ast * item)
{
  return !list->child ? ast_new_children (list, item) :
    ast_new_children (ast_new (item, list->sym), list, item);
}

/**
 * @brief Retrieves the local variable declaration for an identifier within a given scope.
 *
 * Searches for the declaration of the identifier referenced by the AST node `n` in the specified scope. Returns NULL if the identifier is part of a function call or cannot be resolved, or returns the node itself if the identifier is "point".
 *
 * @param n AST node containing the identifier reference.
 * @param stack Stack used for symbol resolution.
 * @param scope AST node representing the scope to search within.
 * @return Ast* The AST node for the local variable declaration, the original node if the identifier is "point", or NULL if not found.
 */
static
Ast * get_local_variable_reference (Ast * n, Stack * stack, Ast * scope)
{
  Ast * identifier = ast_child (ast_find (n, sym_postfix_expression,
					  0, sym_primary_expression),
				sym_IDENTIFIER);
  if (!identifier || ast_ancestor (identifier, 3)->sym == sym_function_call)
    return NULL;
  if (!strcmp (ast_terminal (identifier)->start, "point"))
    return n;
  return ast_identifier_declaration_from_to (stack, ast_terminal (identifier)->start, NULL, scope);
}

/**
 * @brief Cleans up jump statements and transforms certain postfix operations in the AST.
 *
 * Removes jump statements, converts postfix increment/decrement operations that do not reference local variables into assignment expressions, and erases assignment right-hand sides for non-local variable assignments.
 */
static
void set_conditionals (Ast * n, Stack * stack, void * scope)
{
  switch (n->sym) {

  case sym_jump_statement:
    ast_erase (n);
    break;

  case sym_postfix_expression:
    if (n->child[1] && (n->child[1]->sym == sym_INC_OP ||
			n->child[1]->sym == sym_DEC_OP) &&
	!get_local_variable_reference (n, stack, scope)) {
      Ast * parent = n->parent;
      int index = ast_child_index (n);
      AstTerminal * t = ast_terminal_new (n, sym_ADD_ASSIGN, "+=");
      Ast * assign = NN (n, sym_assignment_expression,
			 NN (n, sym_unary_expression, n->child[0]),
			 NN (n, sym_assignment_operator, t),
			 ast_placeholder);
      ast_erase (n->child[1]);
      ast_set_child (parent, index, assign);
    }
    break;
    
  case sym_assignment_expression:
    if (n->child[1] && n->child[2] != ast_placeholder &&
	!get_local_variable_reference (n, stack, scope))
      ast_erase (n->child[2]);
    break;
    
  }
}

/**
 * @brief Wraps an AST node as a stencil access using the "o_stencil" identifier.
 *
 * Creates a unary expression node representing a stencil access and attaches it to the given AST node, marking it for stencil processing.
 *
 * @param n The AST node to be wrapped as a stencil access.
 * @return Ast* The new AST node representing the stencil access.
 */
static inline Ast * o_stencil (Ast * n)
{
  return ast_attach (ast_new_unary_expression (n),
		     NN (n, sym_postfix_expression,
			 NN (n, sym_primary_expression,
			     NB (n, sym_IDENTIFIER, "o_stencil"))));
}

/**
 * @brief Flattens nested AST nodes with identical symbols by shifting their children up.
 *
 * Recursively replaces a node's single child with its own children if the child has the same symbol,
 * effectively collapsing chains of symbolically identical parent-child nodes into a single level.
 */

void ast_shift_children (Ast * n)
{
  while (n->child && !n->child[1] && n->child[0] != ast_placeholder &&
	 n->child[0]->sym == n->sym) {
    n->child = n->child[0]->child;
    for (Ast ** c = n->child; *c; c++)
      if (*c != ast_placeholder)
	(*c)->parent = n;
  }
}

/**
 * @brief Determines if an AST node represents a field type.
 *
 * Checks whether the given AST node corresponds to a field type such as scalar, vertex scalar, vector, face vector, tensor, or symmetric tensor. Returns the type AST node if it matches a field type, otherwise returns NULL.
 *
 * @param n The AST node to check.
 * @param stack The stack used for type resolution.
 * @return Ast* The type node if it is a recognized field type, or NULL otherwise.
 */
static Ast * is_field (Ast * n, Stack * stack)
{
  Ast * type = ast_expression_type (n, stack, true);
  const char * typename = ast_typedef_name (type);
  return typename && (!strcmp (typename, "scalar") ||
		      !strcmp (typename, "vertex scalar") ||
		      !strcmp (typename, "vector") ||
		      !strcmp (typename, "face vector") ||
		      !strcmp (typename, "tensor") ||
		      !strcmp (typename, "symmetric tensor")) ?
    type : NULL;
}

/**
 * @brief Determines if a function call has field-type arguments or is a "neighborp" call.
 *
 * Checks whether any argument of the given function call is a field, or if the function being called is "neighborp".
 *
 * @param function_call The AST node representing the function call.
 * @param stack The stack used for type/context resolution.
 * @return true if the function call has at least one field argument or is "neighborp"; false otherwise.
 */
static
bool has_field_arguments (Ast * function_call, Stack * stack)
{
  Ast * identifier = ast_function_call_identifier (function_call);
  if (identifier && !strcmp (ast_terminal (identifier)->start, "neighborp"))
    return true;
  Ast * arguments = ast_child (function_call, sym_argument_expression_list);
  foreach_item (arguments, 2, argument)
    if (is_field (argument, stack))
      return true;
  return false;
}

/**
 * @brief Cleans up the AST by removing placeholders, invalid, or unreachable nodes.
 *
 * Traverses the AST node and its children, removing or transforming nodes that are placeholders, incomplete, or no longer valid after previous transformation passes. Handles special cases for statements, expressions, and declarations, ensuring the resulting AST is well-formed and contains only valid constructs. Also removes unreachable statements (such as returns inside foreach loops), empty statements, and nodes that do not reference fields or are otherwise unnecessary for stencil analysis.
 *
 * @param n The AST node to clean up.
 * @param stack The traversal stack used for context.
 * @param scope The current scope node.
 * @param init_declarator Indicates whether to erase incomplete init declarators.
 */

void ast_cleanup (Ast * n, Stack * stack, Ast * scope, bool init_declarator)
{
  if (!n->child) // terminals are clean
    return;
  
  int nc = 0, np = 0;
  for (Ast ** c = n->child; *c; c++)
    if (*c != ast_placeholder)
      nc++;
    else
      np++;
  if (!nc) {
    // all children are destroyed, destroy parent
    if (n->sym != sym_argument_expression_list)
      ast_erase (n);
    else {
      Ast * function_call = ast_parent (n, sym_function_call);
      if (function_call &&
	  !is_point_function_call (function_call) &&
	  is_field_access (function_call, stack))
	ast_set_child (n, 0, NN (n, sym_argument_expression_list_item,
				 o_stencil (function_call)));
    }
    return;
  }
  
  if (!np) { // clean
    ast_shift_children (n);
    
    switch (n->sym) {
    
    /**
    Remove function calls which do not reference fields. */

    case sym_function_call:
      if (n->sym == sym_function_call && !is_point_function_call (n) &&
	  !is_field_access (n, stack) && !has_field_arguments (n, stack))
	ast_erase (n);
      break;

    case sym_jump_statement:
      if (n->child[0]->sym == sym_RETURN) {
	if (ast_is_foreach_statement (scope)) {
	  fprintf (stderr, "%s:%d: error: cannot return from a foreach() loop\n",
		   ast_terminal (n->child[0])->file,
		   ast_terminal (n->child[0])->line);
	  exit (1);
	}
	
	/**
	Remove return from compound statements. */

	else if (scope->sym == sym_compound_statement)
	  ast_erase (n);
	
	/**
	Remove values from return statements in point functions. */

	else if (n->child[2] && scope->sym == sym_function_definition) {
	    ast_erase (n->child[1]);
	    n->child[1] = n->child[2];
	    n->child[2] = NULL;
	  }
      }
      
      /**
      Remove goto statements. */

      else if (n->child[0]->sym == sym_GOTO)
	ast_erase (n);
      
      break;

    /**
    Remove empty expression statements. */

    case sym_expression_statement:
      if (!n->child[1] && ast_ancestor (n, 2)->sym == sym_block_item)
	ast_erase (n);
      break;
      
    /**
    Remove labeled (goto) statements. */

    case sym_labeled_statement:
      if (n->child[0]->sym == sym_generic_identifier)
	ast_erase (n);
      break;
	 
    }
      
    return;
  }

  /**
  Fix placeholders. */
  
  switch (n->sym) {

  case sym_selection_statement:
    if (n->child[0]->sym == sym_IF) {
      if (n->child[4] == ast_placeholder &&
	  (!n->child[5] || n->child[6] == ast_placeholder)) {
	ast_erase (n);
	return;
      }
      else if (n->child[2] == ast_placeholder) {
	Ast * list = ast_new (n, sym_block_item_list);
	int ns = 0;
	for (int i = 4; i <= 6 && n->child[i-1]; i += 2)
	  if (n->child[i] != ast_placeholder) {
	    Ast * item = NN (n, sym_block_item, n->child[i]);
	    list = block_list_append (list, item);
	    ns++;
	  }
	if (!ns) {
	  ast_erase (list);
	  ast_erase (n);
	  return;
	}
	ast_erase (n->child[1]);
	ast_erase (n->child[3]);
	if (n->child[5])
	  ast_erase (n->child[5]);
	n->sym = sym_statement;
	if (ns == 1) {
	  ast_erase (n->child[0]);
	  ast_set_child (n, 0, ast_schema (list, sym_block_item_list,
					   0, sym_block_item,
					   0, sym_statement)->child[0]);
	}
	else {
	  Ast * open = ast_terminal_new_char (n, "{");
	  Ast * close = ast_terminal_new_char (n, "}");
	  ast_right_terminal(close)->line = ast_right_terminal(list)->line;
	  ast_erase (n->child[0]);
	  ast_set_child (n, 0, NN (n, sym_compound_statement,
				   open, list, close));
	}
	n->child[1] = NULL;
	stack_push (stack, &n);
	ast_traverse (n, stack, set_conditionals, n);
	ast_pop_scope (stack, n);
	return;
      }
      else if (n->child[4] != ast_placeholder &&
	       n->child[5] &&
	       n->child[6] == ast_placeholder) {
	// IF '(' expression_error ')' statement ELSE _placeholder_
	ast_erase (n->child[5]);
	n->child[5] = NULL;
	return;
      }
      else if (n->child[4] == ast_placeholder) {
	n->child[4] = NN (n, sym_statement,
			  NN (n, sym_expression_statement,
			      NCB (n->child[5], ";")));
	return;
      }
      else
	// fixme: deal with other cases
	assert (false);
    }
    // fixme: deal with other selection statements (switch etc.)
    assert (false);    
    break;
    
  case sym_iteration_statement:
  case sym_for_declaration_statement:
    if (!ast_child (n, sym_statement)) {
      ast_erase (n);
      return;
    }
    break;

  case sym_forin_declaration_statement:
  case sym_forin_statement:
    ast_erase (n);
    break;
    
  case sym_argument_expression_list: {
    Ast * function_call = ast_parent (n, sym_function_call);
    if (is_point_function_call (function_call)) {
      // fixme: undefined function call argument
    }
    else if (is_field_access (function_call, stack)) {
      if (n->child[0] == ast_placeholder) {
	if (!n->child[1])
	  ast_set_child (n, 0, NN (n, sym_argument_expression_list_item,
				   o_stencil (n)));
	else
	  ast_set_child (n, 0, NN (n, sym_argument_expression_list,
				   NN (n, sym_argument_expression_list_item,
				       o_stencil (n->child[1]))));
      }
      if (n->child[1] && n->child[2] == ast_placeholder)
	ast_set_child (n, 2, NN (n, sym_argument_expression_list_item,
				 o_stencil (n->child[1])));
    }
    else
      ast_erase (n);
    return;
  }

  case sym_jump_statement:
    assert (n->child[0]->sym == sym_RETURN);
    n->child[1] = n->child[2];
    n->child[2] = NULL;
    break;
    
  case sym_expression: {
    Ast * parent = n->parent;
    while (parent->sym == sym_expression)
      parent = parent->parent;
    if (is_field_access (parent, stack)) {
      assert (n->child[1]);
      if (n->child[0] == ast_placeholder)
	ast_set_child (n, 0, NN (n, sym_expression, o_stencil (n->child[1])));
      if (n->child[2] == ast_placeholder)
	ast_set_child (n, 2, o_stencil (n->child[1]));
      return;
    }
    // fall through
  }

  case sym_init_declarator_list: {
    if (n->child[0] == ast_placeholder) {
      if (!n->child[2] || n->child[2] == ast_placeholder) {
	ast_erase (n);
	return;
      }
      n->child[0] = n->child[2];
    }
    ast_erase (n->child[1]);
    n->child[1] = NULL;
    ast_shift_children (n);
    break;
  }
    
  case sym_block_item_list:
    if (n->child[0] == ast_placeholder)
      n->child[0] = n->child[1];
    n->child[1] = NULL;
    ast_shift_children (n);    
    break;

  case sym_array_access:
    if (!is_field_access (n, stack))
      ast_erase (n);
    else {
      assert (n->child[2] == ast_placeholder);
      ast_set_child (n, 2, NN (n, sym_expression, o_stencil (n)));
    }      
    break;

  case sym_function_call:
    if (n->child[0] != ast_placeholder)
      assert (!is_field_access (n, stack));
    ast_erase (n);
    break;

  case sym_declarator:
  case sym_direct_declarator:
  case sym_declaration:
  case sym_labeled_statement:
  case sym_cast_expression:
  case sym_compound_statement:
  case sym_initializer:
  case sym_postfix_initializer:
  case sym_initializer_list:
  case sym_foreach_dimension_statement:
  case sym_postfix_expression:
  case sym_conditional_expression:
  case sym_relational_expression:
  case sym_equality_expression:
  case sym_and_expression:
  case sym_exclusive_or_expression:
  case sym_inclusive_or_expression:
  case sym_logical_and_expression:
  case sym_logical_or_expression:    
  case sym_primary_expression:
  case sym_unary_expression:
  case sym_expression_statement:
  case sym_additive_expression:
  case sym_shift_expression:
  case sym_multiplicative_expression:
    assert (!is_field_access (n, stack));
    ast_erase (n);
    return;

  case sym_foreach_statement:
    ast_erase (n);
    return;
    
  case sym_function_definition:
    // keep incomplete function definitions and foreach statements
    return;
    
  /**
  We need to keep incomplete init declarators for the
  undefined_variables() function below. */
    
  case sym_init_declarator:
    if (init_declarator)
      ast_erase (n);
    break;
    
  case sym_assignment_expression:
    if (n->child[0] == ast_placeholder || init_declarator) {
      ast_erase (n);
      return;
    }

    /**
    Same as above for incomplete assignments. */
    
    break;

  case sym_macro_statement:
    if (n->child[0] == ast_placeholder) {
      Ast * statement = ast_ancestor (n, 2);
      ast_set_child (statement, 0, n->child[1]);
    }
    else
      ast_erase (n);
    break;
    
  default:
#if 1   
    ast_print_tree (n, stderr, 0, 0, -1);
    abort();
#endif
    //    ast_print (n, stderr, 0);
    ;

  }
}

/**
 * @brief Moves a field access into its own expression statement and records access type.
 *
 * Isolates a field access AST node as a separate expression statement, wrapping it as an assignment, read, or overflow access as appropriate. The new statement is inserted before or after the parent statement based on the `after` flag.
 *
 * @param parent The parent AST node containing the field access.
 * @param n The AST node representing the field access.
 * @param after If true, inserts the new statement after the parent item; otherwise, inserts before.
 * @param overflow If true, treats the access as an overflow operation.
 * @return Ast* The newly created statement node, or NULL if the field access is already isolated.
 */

static
Ast * move_field_access (Ast * parent, Ast * n, bool after, bool overflow)
{
  Ast * assignment = ast_parent (n, sym_assignment_expression);
  Ast * op = ast_child (assignment, sym_assignment_operator);
  if (!op && ast_ancestor (assignment, 2)->sym == sym_expression_statement)
    return NULL; // already on its own in an expression statement
  Ast * item = ast_block_list_get_item (parent);
  Ast * postfix = n->parent;
  if (op || overflow) {
    Ast * assign = ast_attach (ast_new_unary_expression (parent),
			       postfix);
    AstTerminal * func = NB(assign, sym_IDENTIFIER,
			    op ? (op->child[0]->sym == token_symbol('=') ?
				  "_assign" : "r_assign") :
			    "_overflow");
    postfix = NN(parent, sym_postfix_expression,
		 NN(parent, sym_function_call,
		    NN(parent, sym_postfix_expression,
		       NN(parent, sym_primary_expression,
			  func)),
		    NCB(assign, "("),
		    NN(parent, sym_argument_expression_list,
		       NN(parent, sym_argument_expression_list_item,
			  assign)),
		    NCA(assign, ")")));
  }
  Ast * assign = ast_attach (ast_new_unary_expression (parent),
			     postfix);
  Ast * statement = NN(parent, sym_statement,
		       NN(parent, sym_expression_statement,
			  NN(parent, sym_expression,
			     assign),
			  NCA(assign, ";")));
  if (after)
    return ast_block_list_append (item->parent, item->sym, statement);
  else
    return ast_block_list_insert_before2 (item, statement);
}

typedef struct {
  Ast * scope;
  bool parallel, overflow, nowarning;
  bool undefined;
  int index;
} Undefined;

/**
 * @brief Moves field accesses into separate expression statements.
 *
 * Recursively traverses the AST and isolates each field access by moving it into its own expression statement, placed before its parent statement or declaration. This ensures that field reads and writes are explicitly separated for further stencil analysis and transformation.
 */

void move_field_accesses (Ast * n, Stack * stack, Undefined * u)
{
  if (n == ast_placeholder)
    return;
  
  Ast * scope = ast_push_declarations (n, stack);

  switch (n->sym) {

  case sym_assignment_expression: {
    /**
    In assignment expressions the right-hand-side must be evaluated
    before the left-hand-side. */
    
    Ast ** c;
    for (c = n->child; *c; c++);
    for (c--; c >= n->child; c--)
      move_field_accesses (*c, stack, u);
    break;
  }

  default:
    if (n->child)
      for (Ast ** c = n->child; *c; c++)
	move_field_accesses (*c, stack, u);
  }

  if (is_field_access (n, stack)) {
    Ast * parent = n->parent;
    while (parent &&
	   parent->sym != sym_statement &&
	   parent->sym != sym_declaration)
      parent = parent->parent;
    assert (parent);
    if (parent->sym == sym_statement) {
      switch (parent->child[0]->sym) {

      case sym_jump_statement:
      case sym_selection_statement:
      case sym_expression_statement: {
	move_field_access (parent, n, false, u->overflow);
	break;
      }
	
      default:  // fixme: deal with other statements
	ast_print_tree (parent->child[0], stderr, 0, 0, 1);
	abort();
      }
    }
    else {
      assert (parent->sym == sym_declaration);
      move_field_access (parent, n, true, u->overflow);
    }
  }
  else
    ast_cleanup (n, stack, u->scope, false);
  
  ast_pop_scope (stack, scope);
}

/**
 * @brief Marks an AST node as undefined within a given scope.
 *
 * Sets the value of the terminal AST node to indicate it is undefined in the specified scope.
 *
 * @param n The AST node to mark as undefined.
 * @param scope The scope in which the node is considered undefined.
 */

static inline void set_undefined (Ast * n, Ast * scope)
{
  ast_terminal (n)->value = scope;
}

/**
 * @brief Checks if an AST node is marked as undefined within a given scope.
 *
 * @param n The AST node to check.
 * @param scope The scope used to identify undefined status.
 * @return true if the node's value matches the scope, indicating it is undefined; false otherwise.
 */
static inline bool is_undefined (Ast * n, Ast * scope)
{
  return ast_terminal (n)->value == scope;
}

/**
 * @brief Retrieves the variable declaration AST node for an identifier within a given scope.
 *
 * Searches for the identifier referenced by the AST node, skipping placeholders, function calls, and the special "point" identifier. Returns the corresponding declaration node if found; otherwise, prints an error and exits.
 *
 * @param n AST node potentially referencing a variable.
 * @param stack Stack used for symbol resolution.
 * @param scope Scope in which to search for the variable declaration.
 * @return Ast* Declaration node for the referenced variable, or NULL if not applicable.
 */
static Ast * get_variable_reference (Ast * n, Stack * stack, Ast * scope)
{
  Ast * identifier = ast_find (n, sym_postfix_expression,
			       0, sym_primary_expression);
  if (identifier->child[0] == ast_placeholder)
    return NULL;
  identifier = ast_find (n, sym_postfix_expression,
			 0, sym_primary_expression,
			 0, sym_IDENTIFIER);
  if (!identifier || ast_ancestor (identifier, 3)->sym == sym_function_call)
    return NULL;
  if (!strcmp (ast_terminal (identifier)->start, "point"))
    return NULL;
  Ast * ref = ast_identifier_declaration_from_to
    (stack, ast_terminal (identifier)->start, NULL, scope);
  if (!ref) {
    fprintf (stderr, "%s:%d: error: undeclared identifier '%s'\n",
	     ast_terminal (identifier)->file,
	     ast_terminal (identifier)->line,
	     ast_terminal (identifier)->start);
    exit (1);
  }
  return ref;
}

/**
 * @brief Marks iterators as undefined if they are modified in unsupported ways.
 *
 * Traverses AST nodes representing postfix increment/decrement or assignment expressions and marks the corresponding iterator variable as undefined within the given scope if it is modified in a way that is not a simple assignment.
 *
 * @param n The AST node to analyze.
 * @param stack The stack used for variable resolution.
 * @param data Pointer to an Undefined structure containing the current scope.
 */
static
void undefined_iterators (Ast * n, Stack * stack, void * data)
{
  Undefined * undef = data;

  switch (n->sym) {

  case sym_postfix_expression:
    if (n->child[1] && (n->child[1]->sym == sym_INC_OP ||
			n->child[1]->sym == sym_DEC_OP)) {
      Ast * ref = get_variable_reference (n, stack, NULL);
      if (!ref)
	return;
      set_undefined (ref, undef->scope);
    }
    break;

  case sym_assignment_expression:
    if (n->child[1]) {
      Ast * ref = get_variable_reference (n, stack, NULL);
      if (!ref)
	return;
      if (n->child[1]->child[0]->sym != token_symbol('='))
	set_undefined (ref, undef->scope);
    }
    return;
    
  }
}

/**
 * @brief Checks if a parameter is marked as undefined for stencil analysis.
 *
 * If the given AST node represents a parameter declaration with type `_stencil_undefined`, returns the corresponding type identifier node; otherwise, returns NULL.
 *
 * @param n AST node to check.
 * @return Ast* The type identifier node if the parameter is undefined, or NULL.
 */
static Ast * is_undefined_parameter (const Ast * n)
{
  if (n->sym == sym_IDENTIFIER) n = ast_parent (n, sym_parameter_declaration);
  Ast * identifier;
  if ((identifier = ast_schema (n, sym_parameter_declaration,
				0, sym_declaration_specifiers,
				0, sym_type_specifier,
				0, sym_types,
				0, sym_TYPEDEF_NAME)) &&
      !strcmp (ast_terminal (identifier)->start, "_stencil_undefined"))
    return identifier;
  return NULL;
}

/**
 * @brief Determines if an AST node represents a local variable declaration within a given scope.
 *
 * Checks if the AST node corresponds to the special "point" identifier or is present in the stack before reaching the specified scope.
 *
 * @param n The AST node to check.
 * @param stack The stack representing the current scope chain.
 * @param scope The scope boundary for the search.
 * @return true if the node is a local declaration within the scope, false otherwise.
 */
static
bool is_local_declaration (Ast * n, Stack * stack, Ast * scope)
{
  if (ast_terminal (n)->start && !strcmp (ast_terminal (n)->start, "point"))
    return true;
  Ast ** d;
  for (int i = 0; (d = stack_index (stack, i)); i++)
    if (*d == n)
      return true;
    else if (*d == scope)
      return false;
  return false;
}

/**
 * @brief Finds the nearest enclosing foreach statement in the stack.
 *
 * Traverses the stack from top to bottom and returns the first AST node that represents a foreach statement.
 *
 * @return Ast* Pointer to the nearest foreach statement AST node.
 * @note This function asserts if no foreach statement is found in the stack.
 */
static
Ast * calling_foreach (Stack * stack)
{
  Ast ** d;
  for (int i = 0; (d = stack_index (stack, i)); i++)
    if (ast_is_foreach_statement (*d))
      return *d;
  assert (false);
  return NULL;
}

/**
 * @brief Reports an error if a non-local variable modified in a foreach loop or point function is not listed as a reduction.
 *
 * Checks whether the given variable is included in the reduction list of the enclosing foreach loop or point function. If not, prints an error message and terminates the program.
 *
 * @param n AST node representing the variable being modified.
 * @param stack Call stack for context.
 * @param scope AST node representing the enclosing foreach loop or function definition.
 */
static
void check_missing_reductions (Ast * n, Stack * stack, Ast * scope)
{
  Ast * list = ast_find (ast_schema (scope, sym_foreach_statement,
				     2, sym_argument_expression_list),
			 sym_reduction_list);
  foreach_item (list, 1, reduction) {
    Ast * identifier = ast_schema (reduction, sym_reduction,
				   4, sym_reduction_array,
				   0, sym_generic_identifier,
				   0, sym_IDENTIFIER);
    if (!strcmp (ast_terminal (identifier)->start,
		 ast_terminal (n)->start))
      return;
  }
  AstTerminal * t = ast_left_terminal (scope);
  if (ast_is_foreach_statement (scope))
    fprintf (stderr,
	     "%s:%d: error: non-local variable '%s' is modified by "
	     "this foreach loop:\n"
	     "%s:%d: error: use a loop-local variable, a reduction operation\n"
	     "%s:%d: error: or a serial loop to get rid of this error\n",
	     t->file, t->line, ast_terminal (n)->start,
	     t->file, t->line, t->file, t->line);
  else if (scope->sym == sym_function_definition) {
    AstTerminal * f = ast_left_terminal (calling_foreach (stack));
    fprintf (stderr,
	     "%s:%d: error: non-local variable '%s' is modified by this "
	     "point function\n"
	     "%s:%d: error: use a local variable or\n"
	     "%s:%d: error: a serial loop (here) to get rid of this error\n",
	     t->file, t->line, ast_terminal (n)->start,
	     t->file, t->line,
	     f->file, f->line);
  }
  else
    assert (false);
  exit (1);
}

/**
 * @brief Determines if the given AST node refers to the POINT_VARIABLES identifier.
 *
 * Traverses the AST from the provided reference node to locate an identifier and checks if it matches "POINT_VARIABLES".
 *
 * @param ref The AST node reference to check.
 * @return true if the node refers to POINT_VARIABLES, false otherwise.
 */
static
bool is_point_variable (const Ast * ref)
{
  Ast * identifier =
    ast_schema (ast_parent (ref, sym_function_definition), sym_function_definition,
		0, sym_function_declaration,
		1, sym_declarator,
		0, sym_direct_declarator,
		0, sym_direct_declarator,
		0, sym_generic_identifier,
		0, sym_IDENTIFIER);
  return identifier && !strcmp (ast_terminal (identifier)->start, "POINT_VARIABLES");
}

/**
 * @brief Determines if an AST node is an argument to a foreach statement.
 *
 * Traverses the AST upwards from the given node to check if it is part of the argument list of a foreach statement.
 *
 * @param n The AST node to check.
 * @return true if the node is an argument to a foreach statement, false otherwise.
 */
bool ast_is_foreach_parameter (Ast * n)
{
  n = ast_parent (n, sym_argument_expression_list);
  while (n && n->sym == sym_argument_expression_list)
    n = n->parent;
  return ast_is_foreach_statement (n);
}

/**
 * @brief Propagates and removes undefined variable references in the AST.
 *
 * Traverses the AST to mark variables as undefined if they are assigned undefined values, are non-local and illegally modified in parallel contexts, or are referenced before declaration. Removes or erases AST nodes corresponding to undefined or illegal variable usages, including increment/decrement operations and assignments. Handles point function calls that may modify variables by address, and cleans up incomplete or invalid loop constructs. Ensures the AST does not contain references to undefined or illegal variables after transformation.
 */
static
void undefined_variables (Ast * n, Stack * stack, void * data)
{
  Undefined * undef = data;
  
  switch (n->sym) {

  case sym_IDENTIFIER: {

    /**
    Variable "point" is always defined. */
    
    if (!ast_terminal (n)->start || !strcmp (ast_terminal (n)->start, "point"))
      break;
    
    Ast * ref = ast_identifier_declaration (stack, ast_terminal (n)->start);
    if (ref && (is_point_variable (ref) || is_undefined_parameter (ref)))
      ref = NULL;
    
    /**
    Reset state when the variable is declared. */
    
    if (ref == n) {
      set_undefined (ref, NULL);
      break;
    }

    /**
    Only consider variable identifiers i.e. not struct members,
    function identifiers or foreach parameters. */
    
    if (n->parent->sym != sym_primary_expression ||
	ast_ancestor (n, 3)->sym == sym_function_call ||
	ast_is_foreach_parameter (n))
      break;

    /**
    If the variable is undeclared or marked as undefined, we remove it. */
    
    if (!ref || is_undefined (ref, undef->scope)) {
      ast_erase (n);
      undef->undefined = true;
      return;
    }

    /**
    Check for global variables. */

    if (!is_local_declaration (ref, stack, undef->scope)) {

      /**
      Check whether this is an (illegal parallel) `for (IDENTIFIER in ...)` construct. */

      if (undef->parallel) {
	Ast * assignment = ast_parent (n, sym_assignment_expression);
	if (ast_is_identifier_expression (assignment)) {
	  while (assignment->parent->sym == sym_expression)
	    assignment = assignment->parent;
	  if (ast_schema (assignment->parent, sym_forin_statement,
			  2, sym_expression))
	    check_missing_reductions (ref, stack, undef->scope);
	}
      }
    }
      
    break;
  }
    
  case sym_init_declarator: {
    if (n->child[1] && n->child[2] == ast_placeholder) {
      Ast * ref = ast_find (n->child[0], sym_direct_declarator,
			    0, sym_generic_identifier,
			    0, sym_IDENTIFIER);
      if (ref)
	set_undefined (ref, undef->scope);
      ast_erase (n->child[1]);
      n->child[1] = NULL;
    }
    break;
  }

  /**
  Incremented non-local variables must be removed. */
    
  case sym_postfix_expression: {
    Ast * ref;
    if (n->child[1] && n->child[0] != ast_placeholder &&
	(n->child[1]->sym == sym_INC_OP ||
	 n->child[1]->sym == sym_DEC_OP) &&
	(ref = get_variable_reference (n, stack, NULL)) &&
	!is_local_declaration (ref, stack, undef->scope)) {
      if (undef->parallel)
	check_missing_reductions (ref, stack, undef->scope);
      set_undefined (ref, undef->scope);
      ast_erase (n);
      undef->undefined = true;
      return;
    }
    break;
  }

  case sym_assignment_expression:
    if (n->child[1] && n->child[0] != ast_placeholder) {
      Ast * ref = get_variable_reference (n, stack, NULL);
      if (!ref)
	return;
      bool local = is_local_declaration (ref, stack, undef->scope);
      if (!local || n->child[2] == ast_placeholder) {
	if (!local && undef->parallel)
	  check_missing_reductions (ref, stack, undef->scope);
	set_undefined (ref, undef->scope);
	ast_erase (n);
	undef->undefined = true;
	return;
      }
      else if (is_undefined (ref, undef->scope) &&
	       n->child[1]->child[0]->sym == token_symbol('='))
	set_undefined (ref, NULL);
    }
    break;

  /**
  Point function calls which take the address of a variable as
  argument: we assume the point function modifies the variable. */

  case sym_function_call:
    if (is_point_function_call (n)) {
      Ast * arguments = ast_child (n, sym_argument_expression_list);
      foreach_item (arguments, 2, argument)
	if (ast_schema (argument, sym_argument_expression_list_item,
			0, sym_assignment_expression,
			0, sym_conditional_expression) &&
	    ast_find (argument, sym_unary_expression,
		      0, sym_unary_operator,
		      0, token_symbol ('&'))) {
	  Ast * ref = get_variable_reference (argument, stack, NULL);
	  if (!ref)
	    break;
#if 0 // we assume non-local variables are not modified	  
	  if (undef->parallel &&
	      !is_local_declaration (ref, stack, undef->scope))
	    check_missing_reductions (ref, stack, undef->scope);
#endif
	  set_undefined (ref, undef->scope);
	  ast_erase (argument);
	  undef->undefined = true;
	}
    }
    break;
    
    
  /**
  For loops where the condition is undefined. We need to set all the
  corresponding iterative variables to undefined. */
    
  case sym_expression:
    if ((n->parent->sym == sym_for_declaration_statement ||
	 n->parent->sym == sym_iteration_statement) &&
	n->parent->child[3] == ast_placeholder)
      ast_traverse (n, stack, undefined_iterators, undef);
    break;

  /**
  Cleanup of incomplete for loops. */

  case sym_for_declaration_statement:
  case sym_iteration_statement: {
    Ast ** c;
    for (c = n->child; *c; c++)
      if (*c == ast_placeholder)
	break;
    if (!*c)
      break;
    Ast * statement = ast_child (n, sym_statement);
    for (c = n->child; *c; c++)
      if (*c != ast_placeholder && *c != statement)
	ast_erase (*c);
    if (n->sym == sym_for_declaration_statement)
      n = n->parent;
    if (statement)
      ast_set_child (n->parent, ast_child_index (n), statement->child[0]);
    break;
  }
    
  }

  ast_cleanup (n, stack, undef->scope, false);
}

/**
 * @brief Finds a stencil function corresponding to a given function definition.
 *
 * Searches for a previously defined stencil function whose name matches the given function definition, using the "_stencil_" prefix. Returns the matching stencil function AST node if found, or NULL otherwise.
 *
 * @param function_definition The AST node representing the original function definition.
 * @return Ast* The AST node of the matching stencil function, or NULL if not found.
 */

static Ast * get_stencil_function (Ast * function_definition)
{
  Ast * identifier = ast_function_identifier (function_definition);
  if (!identifier)
    return NULL;
  char * stencil_name = NULL;
  str_append (stencil_name, "_stencil_", ast_terminal (identifier)->start);
  Ast * stencil;
  if (((stencil = ast_schema (ast_ancestor (function_definition, 3),
			      sym_translation_unit,
			      1, sym_external_declaration,
			      0, sym_function_definition)) ||
       ((stencil = ast_schema (ast_ancestor (function_definition, 4),
			       sym_translation_unit,
			       1, sym_external_declaration,
			       0, sym_external_foreach_dimension)) &&
	(stencil = ast_child (stencil, sym_function_definition)))) &&
      (identifier = ast_function_identifier (stencil)) &&
      !strcmp (ast_terminal (identifier)->start, stencil_name)) {
    free (stencil_name);
    return stencil;
  }
  free (stencil_name);
  return NULL;
}

/**
 * @brief Finds a function declaration by name, considering possible _x, _y, or _z suffix rotations.
 *
 * Searches for a function declaration matching the given name within the specified AST range. If not found and the name ends with an underscore followed by 'x', 'y', or 'z', attempts to find a declaration by rotating the suffix among 'x', 'y', and 'z'.
 *
 * @param name The function name, which may be modified temporarily during the search.
 * @return Ast* The matching function declaration, or NULL if none is found.
 */

static
Ast * identifier_function_declaration (Stack * stack, char * name,
				       Ast * from, Ast * to)
{
  Ast * n = ast_identifier_declaration_from_to (stack, name, from, to);
  int len;
  if (!n && (len = strlen(name)) > 2 &&
      name[len - 2] == '_' && strchr ("xyz", name[len - 1])) {
    char o = name[len - 1], c = 'x';
    while (!n && c <= 'z') {
      name[len - 1] = c++;
      n = ast_identifier_declaration_from_to (stack, name, from, to);
    }
    name[len - 1] = o;
  }
  return n;
}

/**
 * @brief Retrieves the function definition AST node for a given identifier.
 *
 * Searches for the function definition corresponding to the provided identifier, starting from the given declaration and traversing parent nodes as needed. Returns NULL if no matching function definition is found.
 *
 * @param identifier The AST node representing the function identifier.
 * @param declaration The starting declaration node for the search.
 * @return Ast* The AST node for the function definition, or NULL if not found.
 */
Ast * ast_get_function_definition (Stack * stack, Ast * identifier, Ast * declaration)
{
  if (!identifier)
    return NULL;
  Ast * declaration1 = identifier_function_declaration
    (stack, ast_terminal (identifier)->start, declaration, NULL);
  Ast * function_definition = declaration1;
  while (function_definition &&
	 function_definition->sym != sym_declaration &&
	 function_definition->sym != sym_function_definition)
    function_definition = function_definition->parent;
  if (!function_definition)
    return NULL;
  if (function_definition->sym == sym_function_definition)
    return function_definition;
  if (!ast_schema (function_definition, sym_declaration,
		   1, sym_init_declarator_list,
		   0, sym_init_declarator,
		   0, sym_declarator,
		   0, sym_direct_declarator,
		   0, sym_direct_declarator,
		   0, sym_generic_identifier,
		   0, sym_IDENTIFIER))
    return NULL;
  if (declaration1 == declaration)
    return NULL;
  return ast_get_function_definition (stack, identifier, declaration1);
}

/**
 * @brief Appends a function declaration to the appropriate location in the AST.
 *
 * Depending on the parent node's type, inserts the function declaration either directly into the external declarations block or wraps it within a foreach dimension node before appending. Asserts if the parent context is invalid.
 */
static void append_function_declaration (Ast * parent, Ast * declaration)
{
  if (parent->parent->sym == sym_external_declaration)
    ast_block_list_append (ast_ancestor (parent, 2),
			   sym_external_declaration, declaration);
  else if (parent->parent->sym ==
	   sym_external_foreach_dimension) {
    int index = ast_child_index (parent);
    parent->parent->child[index] = ast_placeholder;
    Ast * foreach_dimension = ast_copy (parent->parent);
    parent->parent->child[index] = parent;
    ast_set_child (foreach_dimension, index, declaration);
    ast_block_list_append (ast_ancestor (parent, 3),
			   sym_external_declaration, foreach_dimension);
  }
  else
    assert (false);
}

/**
 * @brief Replaces a function call AST node with a default stencil call.
 *
 * Transforms the given AST node into a call to the "default_stencil" function, constructing its arguments and initializer list from the original node. If no field arguments are present, the node is destroyed.
 */
static void default_stencil (Ast * n, Stack * stack, void * scope)
{
  Ast * initializer = NN (n, sym_postfix_initializer,
			  NCA (n, "{"), ast_placeholder, NCA (n, "}"));
  ast_replace_child (n, 0,
		     NN (n, sym_postfix_expression,
			 NN (n, sym_primary_expression,
			     NB (n, sym_IDENTIFIER, "default_stencil"))));
  Ast * arguments = ast_child (n, sym_argument_expression_list);
  Ast * list = ast_new (n, sym_initializer_list);
  Ast * args = NN (n, sym_argument_expression_list,
		   NN (n, sym_argument_expression_list_item,
		       ast_attach (ast_new_unary_expression (n),
				   NN (n, sym_postfix_expression,
				       NN (n, sym_primary_expression,
					   NB (n, sym_IDENTIFIER, "point"))))));
  ast_set_child (n, 2, args);
  args = ast_list_append (args, sym_argument_expression_list_item, initializer, ",");
  ast_set_child (initializer, 1, list);
  foreach_item (arguments, 2, argument)
    if (argument != ast_placeholder && is_field (argument, stack)) {
      if (!list->child)
	ast_attach (list, NN (list, sym_initializer, argument->child[0]));
      else
	list = ast_list_append (list, sym_initializer, argument->child[0], ",");
    }
  ast_destroy (arguments);
  if (!list->child) { // no field arguments
    list->child = allocate (ast_get_root(n)->alloc, sizeof (Ast *));
    list->child[0] = NULL;
    ast_destroy (n);
  }
}

/**
 * @brief Creates an AST node representing a NULL expression.
 *
 * Constructs a unary expression AST node with the identifier "NULL" as its operand.
 *
 * @return Ast* The newly created NULL expression AST node.
 */
static
Ast * null_expression (Ast * n)
{		     
  return ast_attach (ast_new_unary_expression (n),
		     NN (n, sym_postfix_expression,
			 NN (n, sym_primary_expression,
			     NA (n, sym_IDENTIFIER, "NULL"))));
}

/**
 * @brief Marks a function parameter as undefined by setting its type to `_stencil_undefined`.
 *
 * Replaces the parameter's type with `_stencil_undefined` and updates its declarator accordingly.
 */
static
void set_undefined_parameter (Ast * parameter)
{
  AstTerminal * pointer = NCB (parameter->child[1], "*");
  Ast * undefined = 
    NN (parameter, sym_declaration_specifiers,
	NN (parameter, sym_type_specifier,
	    NN (parameter, sym_types,
		NB (parameter->child[1], sym_TYPEDEF_NAME,
		    "_stencil_undefined"))));
  ast_erase (parameter->child[0]);
  ast_set_child (parameter, 0, undefined);
  ast_replace_child (parameter, 1,
		     NN (parameter, sym_declarator,
			 NN (parameter, sym_pointer,
			     pointer),
			 NN (parameter, sym_direct_declarator,
			     NN (parameter, sym_generic_identifier,
				 ast_find (parameter->child[1],
					   sym_direct_declarator,
					   0, sym_generic_identifier,
					   0, sym_IDENTIFIER)))));
  ast_flatten (parameter, ast_left_terminal (parameter));
}

/**
 * @brief Marks point function parameters as undefined if corresponding call arguments are undefined.
 *
 * For each argument in a function call, sets the corresponding parameter in the point function as undefined if the argument is a placeholder. Exits with an error if there are more arguments than parameters.
 */

static
void set_undefined_parameters (Ast * point_function, const Ast * function_call)
{
  Ast * arguments = ast_schema (function_call, sym_function_call,
				2, sym_argument_expression_list);
  Ast * parameters = ast_find (point_function, sym_direct_declarator,
			       2, sym_parameter_type_list,
			       0, sym_parameter_list);
  Ast * parameter = parameters ? (parameters->child[1] ? parameters->child[2] :
				  parameters->child[0]) : NULL;
  foreach_item (arguments, 2, argument) {
    if (!parameter) {
      fprintf (stderr, "%s:%d: error: too many arguments in function call\n",
	       ast_left_terminal (function_call)->file, ast_left_terminal (function_call)->line);
      exit (1);
    }
    if (argument == ast_placeholder)
      set_undefined_parameter (parameter);
    parameters = parameters && parameters != ast_placeholder &&
      parameters->child[1] ? parameters->child[0] : NULL;
    parameter = parameters ? (parameters->child[1] ? parameters->child[2] :
			      parameters->child[0]) : NULL;
    arguments = _list;
  }
}

/**
 * @brief Replaces point function calls in the AST with corresponding stencil function calls.
 *
 * For each point function call found in the AST node, this function locates or generates the corresponding stencil function (prefixed with "_stencil_"), ensuring that undefined arguments and parameters are handled consistently. If the stencil function does not exist, it is created by transforming a copy of the original function. The function also replaces undefined arguments with NULL expressions and issues warnings or errors for unresolved or mismatched cases.
 */
static void point_function_calls (Ast * n, Stack * stack, void * data)
{
  Undefined * undef = data;
  ast_cleanup (n, stack, undef->scope, false);
  
  if (n->sym != sym_function_call || !is_point_function_call (n))
    return;

  Ast * identifier = ast_function_call_identifier (n);
  if (!identifier) {
    if (!undef->nowarning)
      fprintf (stderr,
	       "%s:%d: warning: stencils: "
	       "cannot analyze point function pointers\n",
	       ast_left_terminal (n)->file, ast_left_terminal (n)->line);
    default_stencil (n, stack, undef->scope);
    return;
  }

  Ast * function_declaration = NULL;
  Ast * function_definition = ast_get_function_definition (stack, identifier, NULL);
  if (function_definition) {
    function_declaration =
      identifier_function_declaration (stack, ast_terminal (identifier)->start,
				       NULL, NULL);
    while (function_declaration &&
	   function_declaration->sym != sym_declaration &&
	   function_declaration->sym != sym_function_definition)
    function_declaration = function_declaration->parent;
    if (function_declaration->sym != sym_declaration)
      function_declaration = NULL;
  }
  else {
    if (!undef->nowarning)
      fprintf (stderr,
	       "%s:%d: warning: stencils: point function '%s' is not defined\n",
	       ast_left_terminal (n)->file, ast_left_terminal (n)->line,
	       ast_terminal (identifier)->start);
    default_stencil (n, stack, undef->scope);
    return;
  }

  assert (!ast_is_stencil_function (function_definition));

  str_prepend (ast_terminal (identifier)->start, "_stencil_");
  
  if (function_definition == undef->scope)
    return; // recursive function call

  /**
  Look for a previously-defined stencil function. */

  Ast * stencil = get_stencil_function (function_definition);
  if (!stencil) {
    
    /**
    We need to create the new stencil function. */  

    Ast * new_stencil = ast_copy (function_definition);
    set_undefined_parameters (new_stencil, n);
    stencil = ast_stencil (new_stencil, undef->parallel, undef->overflow, undef->nowarning);
    if (!stencil) {
      ast_destroy (new_stencil);
      ast_erase (n);
      return;
    }
    else {
      Ast * m = function_definition;
      AstTerminal * t = ast_terminal_new (m, sym_VOID, "void");
      str_append (t->before, " ");
      Ast * specifiers = NN (m, sym_declaration_specifiers,
			     NN (m, sym_storage_class_specifier,
				 ast_terminal_new (m, sym_STATIC, "static")),
			     NN (m, sym_declaration_specifiers,
				 NN (m, sym_type_specifier,
				     NN (m, sym_types, t))));
      ast_replace_child (ast_schema (stencil, sym_function_definition,
				     0, sym_function_declaration),
			 0,
			 specifiers);
      str_prepend (ast_terminal (ast_function_identifier (stencil))->start,
		   "_stencil_");
      append_function_declaration (function_definition, stencil);
    
      /**
      We also create the corresponding declaration if necessary. */

      if (function_declaration) {
	Ast * semicolumn =
	  ast_terminal_new_char ((Ast *) ast_right_terminal (stencil->child[0]),
				 ";");
	Ast * specifiers = ast_copy (ast_schema (stencil, sym_function_definition,
						 0, sym_function_declaration,
						 0, sym_declaration_specifiers));
	Ast * declarator = ast_copy (ast_schema (stencil, sym_function_definition,
						 0, sym_function_declaration,
						 1, sym_declarator));
	Ast * declaration = NN (n, sym_declaration,
				specifiers,
				NN (n, sym_init_declarator_list,
				    NN (n, sym_init_declarator,
					declarator)),
				semicolumn);
	append_function_declaration (function_declaration, declaration);
      }
    }   
  }
    
  /**
  We check that the undefined function call arguments match undefined
  function parameters. */
  
  Ast * arguments = ast_schema (n, sym_function_call,
				2, sym_argument_expression_list);
  Ast * parameters = ast_find (stencil, sym_direct_declarator,
			       2, sym_parameter_type_list,
			       0, sym_parameter_list);
  Ast * parameter = parameters ?
    (parameters->child[1] ? parameters->child[2] :
     parameters->child[0]) : NULL;
  foreach_item (arguments, 2, argument) {
    if (!parameter) {
      fprintf (stderr, "%s:%d: error: too many arguments in function call\n",
	       ast_left_terminal (n)->file, ast_left_terminal (n)->line);
      exit (1);
    }
    if (argument == ast_placeholder) {

      /**
      We check that the undefined argument corresponds with an undefined parameter. */

      if (!is_undefined_parameter (parameter)) {
	fprintf (stderr, "%s:%d: error: stencils: not expecting an undefined "
		 "argument for '%s'\n",
		 ast_left_terminal (n)->file, ast_left_terminal (n)->line,
		 ast_terminal (ast_find (parameter->child[1],
					 sym_IDENTIFIER))->start);
	exit (1);
      }

      /**
      We replace the undefined argument with a NULL value. */

      int index = arguments->child[1] ? 2 : 0;
      ast_replace_child (arguments, index,
			 NN (n, sym_argument_expression_list_item,
			     null_expression (n)));
    }

    /**
    If the parameter is undefined we replace the argument with a
    NULL value. */

    else if (is_undefined_parameter (parameter))
      ast_replace_child (argument, 0, null_expression (argument));

    parameters = parameters && parameters != ast_placeholder &&
      parameters->child[1] ? parameters->child[0] : NULL;
    parameter = parameters ? (parameters->child[1] ? parameters->child[2] :
			      parameters->child[0]) : NULL;
    arguments = _list;
  }
}

/**
 * @brief Finds the nearest ancestor node with a specific symbol within a given scope.
 *
 * Traverses the parent chain of the AST node `n` to locate the closest ancestor whose symbol matches `sym`, stopping if a node with symbol `scope` is encountered.
 *
 * @param n The starting AST node.
 * @param sym The symbol to search for among ancestor nodes.
 * @param scope The symbol indicating the scope boundary for the search.
 * @return Ast* Pointer to the ancestor node with symbol `sym`, or NULL if not found before reaching the scope boundary.
 */

Ast * ast_scope_parent (Ast * n, int sym, int scope)
{
  n = n->parent;
  while (n && n->sym != scope) {
    if (n->sym == sym)
      return n;
    n = n->parent;
  }
  return NULL;
}

/**
 * @brief Marks an AST node as initialized by setting its terminal value to itself.
 *
 * This function is used to indicate that the given AST node has been initialized,
 * typically during AST analysis or transformation passes.
 */
static inline void initialize (Ast * n)
{
#if 0 
  fprintf (stderr, "%s:%d: initialize %s\n",
	   ast_terminal (n)->file,
	   ast_terminal (n)->line,
	   ast_terminal (n)->start);
#endif
  ast_terminal (n)->value = n;
}

/**
 * @brief Checks if an AST node is marked as initialized.
 *
 * @param n The AST node to check.
 * @return true if the node is initialized; false otherwise.
 */
static inline bool is_initialized (Ast * n)
{
  return ast_terminal (n)->value == n;
}
  
/**
 * @brief Marks an AST node as declared within a given scope.
 *
 * Associates the provided scope with the AST node to indicate its declaration context.
 */
static inline void declare (Ast * n, Ast * scope)
{
#if 0 
  fprintf (stderr, "%s:%d: declare %s\n",
	   ast_terminal (n)->file,
	   ast_terminal (n)->line,
	   ast_terminal (n)->start);
#endif
  ast_terminal (n)->value = scope;
}

/**
 * @brief Checks if an AST node is declared within a given scope.
 *
 * @param n The AST node to check.
 * @param scope The scope to compare against.
 * @return true if the node's value matches the scope, false otherwise.
 */
static inline bool is_declared (Ast * n, Ast * scope)
{
  return ast_terminal (n)->value == scope;
}

/**
 * @brief Marks an AST node as used by clearing its value.
 *
 * Sets the value of the terminal AST node to NULL to indicate usage.
 */
static inline void use (Ast * n)
{
  ast_terminal (n)->value = NULL;
}

/**
 * @brief Marks variables in the AST as used, declared, or initialized, and removes unused identifiers.
 *
 * Traverses the AST node to update the usage state of variables based on their context. Identifiers are marked as declared or initialized depending on their position in declarations or assignments. Unused identifiers that are declared but not initialized are removed from the AST.
 */
static
void mark_unused (Ast * n, Stack * stack, void * scope)
{
  switch (n->sym) {

  case sym_IDENTIFIER: {
    if (ast_ancestor (n, 3)->sym == sym_function_call || ast_parent (n, sym_forin_arguments))
      return;
    if (ast_ancestor (n, 2)->sym == sym_direct_declarator) {      
      if (ast_scope_parent (n, sym_struct_declarator, sym_declaration))
	return;
      Ast * parent;
      if (ast_ancestor (n, 4)->sym == sym_forin_declaration_statement ||
	  ast_parent (n, sym_function_declaration) ||
	  ((parent = ast_scope_parent (n, sym_init_declarator, sym_declaration)) && parent->child[1]))
	initialize (n);
      else
	declare (n, scope);
      return;
    }
    if (n->parent->sym != sym_primary_expression)
      return;
    Ast * ref = ast_identifier_declaration_from_to (stack, ast_terminal (n)->start, NULL, scope);
    if (ref) {
      Ast * expr = ast_scope_parent (n, sym_expression, sym_statement);
      if (expr) {
	while (expr->parent->sym == sym_expression)
	  expr = expr->parent;
	if (expr->parent->sym == sym_forin_statement) {
	  initialize (ref);
	  return;
	}
      }
      Ast * assign = ast_parent (n, sym_assignment_expression);
      if (ast_schema (assign, sym_assignment_expression,
		      1, sym_assignment_operator,
		      0, token_symbol ('='))) {
	if (is_declared (ref, scope))
	  initialize (ref);
	return;
      }
      if (is_declared (ref, scope))
	// declared but not initialized
	ast_erase (n);
      else {
#if 0
	fprintf (stderr, "%s:%d: use %s\n",
		 ast_terminal (n)->file,
		 ast_terminal (n)->line,
		 ast_terminal (n)->start);
#endif	
	use (ref);
      }
    }
    break;
  }
    
  }
}

/**
 * @brief Removes references to undefined identifiers from the AST.
 *
 * Recursively traverses the AST node and erases identifier nodes that are not declared or initialized in the current scope, except in certain contexts such as function calls, declarations, or foreach arguments. Cleans up remaining nodes as needed.
 */
static
void remove_undefined (Ast * n, Stack * stack, void * scope)
{
  Ast * ref;
  if (n->sym == sym_IDENTIFIER &&
      ast_ancestor (n, 3)->sym != sym_function_call &&
      ast_ancestor (n, 2)->sym != sym_direct_declarator &&
      n->parent->sym == sym_primary_expression &&
      !ast_parent (n, sym_forin_arguments) &&
      strcmp (ast_terminal (n)->start, "point") &&
      (ref = ast_identifier_declaration (stack, ast_terminal (n)->start)) &&
      (is_declared (ref, scope) || is_initialized (ref))) {
#if 0
    fprintf (stderr, "%s:%d: '%s' undefined %d %d\n",
	     ast_terminal (n)->file,
	     ast_terminal (n)->line,
	     ast_terminal (n)->start,
	     is_declared (ref, scope), is_initialized (ref));
#endif
    ast_erase (n);
  }
  else
    ast_cleanup (n, stack, scope, true);
}

/**
 * @brief Removes unused or undefined identifiers and marks unused parameters as undefined.
 *
 * Traverses the AST to erase identifiers that are declared or initialized but unused within the given scope, and replaces unused function parameters with an undefined type. Also recursively cleans up other nodes as needed.
 */
static
void remove_unused (Ast * n, Stack * stack, void * data)
{
  Undefined * undef = data;
  if (n->sym == sym_IDENTIFIER &&
      (is_declared (n, undef->scope) || is_initialized (n))) {
    Ast * decl = ast_parent (n, sym_function_declaration);
    if (!decl || decl->parent != undef->scope) {
#if 0
      fprintf (stderr, "%s:%d: '%s' unused %d %d\n",
	       ast_terminal (n)->file,
	       ast_terminal (n)->line,
	       ast_terminal (n)->start,
	       is_declared (n, undef->scope), is_initialized (n));
#endif
      ast_erase (n);
      undef->undefined = true;
    }
    else { // (unused) Point function parameters
      Ast * parameter = ast_parent (n, sym_parameter_declaration);
      if (!parameter) { // the function name
	use (n);
	return;
      }
      Ast * name;
      if (!strcmp (ast_terminal (n)->start, "point") &&
	  (name = ast_schema (ast_ancestor (n, 4), sym_parameter_declaration,
			      0, sym_declaration_specifiers,
			      0, sym_type_specifier,
			      0, sym_types,
			      0, sym_TYPEDEF_NAME)) &&
	  !strcmp (ast_terminal (name)->start, "Point")) { // 'Point point' parameter
	use (n);
	return;
      }
#if 0
      fprintf (stderr, "%s:%d: '%s' parameter unused %d %d\n",
	       ast_terminal (n)->file,
	       ast_terminal (n)->line,
	       ast_terminal (n)->start,
	       is_declared (n, undef->scope), is_initialized (n));
#endif
      /**
      We replace the unused parameter with an undefined type. */

      set_undefined_parameter (parameter);
    }
  }
  else
    ast_cleanup (n, stack, undef->scope, true);
}

/**
 * @brief Transforms a foreach loop or point function AST into a stencil AST for boundary condition handling.
 *
 * Recursively analyzes and rewrites the given AST node to isolate stencil field accesses, propagate undefined values, replace point function calls with stencil equivalents, and remove unused or undefined variables. Returns a simplified AST that records stencil read/write accesses, or NULL if no stencil accesses remain.
 *
 * @param n The AST node representing a foreach loop or point function.
 * @param parallel Indicates if the transformation should consider parallel execution.
 * @param overflow Enables overflow handling for stencil accesses.
 * @param nowarning Suppresses warnings during transformation.
 * @return Ast* The transformed stencil AST, or NULL if no stencil accesses are found.
 */

Ast * ast_stencil (Ast * n, bool parallel, bool overflow, bool nowarning)
{
  if (ast_is_foreach_statement (n))
    n->sym = sym_foreach_statement;
  else
    assert (n->sym == sym_function_definition);
  AstRoot * root = ast_get_root (n);
  Stack * stack = root->stack;
  stack_push (stack, &n);
  Undefined u = {n, parallel, overflow, nowarning};
  move_field_accesses (n, stack, &u);
  Ast * m = ast_is_foreach_statement (n) ? ast_child (n, sym_statement) : n;
  do {
    u.undefined = false;
    ast_traverse (m, stack, undefined_variables, &u);
  } while (u.undefined);
  ast_traverse (m, stack, point_function_calls, &u);

  Ast * statement =  (ast_is_foreach_statement (n) ? ast_child (n, sym_statement) :
		      ast_child (n, sym_compound_statement));
  do {
    ast_traverse (n, stack, mark_unused, n);
    ast_traverse (statement, stack, remove_undefined, n);
    u.undefined = false;
    ast_traverse (n, stack, remove_unused, &u);
    ast_traverse (n, stack, undefined_variables, &u);
  } while (u.undefined);
  
  ast_pop_scope (stack, n);
  if (ast_is_foreach_statement (n) && !ast_child (n, sym_statement))
    return NULL;
  if (n->sym == sym_function_definition && !ast_child (n, sym_compound_statement))
    return NULL;
  if (n->sym == sym_foreach_statement)
    n->sym = sym_macro_statement;
  return CHECK (n);
}
