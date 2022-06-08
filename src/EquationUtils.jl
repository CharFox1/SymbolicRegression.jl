module EquationUtilsModule

import ..CoreModule: CONST_TYPE, Node, copy_node, Options

# for symbolic constraints
using PyCall

# Count the operators, constants, variables in an equation
function count_nodes(tree::Node)::Int
    if tree.degree == 0
        return 1
    elseif tree.degree == 1
        return 1 + count_nodes(tree.l)
    else
        return 1 + count_nodes(tree.l) + count_nodes(tree.r)
    end
end

# Count the max depth of a tree
function count_depth(tree::Node)::Int
    if tree.degree == 0
        return 1
    elseif tree.degree == 1
        return 1 + count_depth(tree.l)
    else
        return 1 + max(count_depth(tree.l), count_depth(tree.r))
    end
end

# Count the number of unary operators in the equation
function count_unary_operators(tree::Node)::Int
    if tree.degree == 0
        return 0
    elseif tree.degree == 1
        return 1 + count_unary_operators(tree.l)
    else
        return 0 + count_unary_operators(tree.l) + count_unary_operators(tree.r)
    end
end

# Count the number of binary operators in the equation
function count_binary_operators(tree::Node)::Int
    if tree.degree == 0
        return 0
    elseif tree.degree == 1
        return 0 + count_binary_operators(tree.l)
    else
        return 1 + count_binary_operators(tree.l) + count_binary_operators(tree.r)
    end
end

# Count the number of operators in the equation
function count_operators(tree::Node)::Int
    return count_unary_operators(tree) + count_binary_operators(tree)
end

# Count the number of constants in an equation
function count_constants(tree::Node)::Int
    if tree.degree == 0
        if tree.constant
            return 1
        else
            return 0
        end
    elseif tree.degree == 1
        return 0 + count_constants(tree.l)
    else
        return 0 + count_constants(tree.l) + count_constants(tree.r)
    end
end

"""
Compute the complexity of a tree.

By default, this is the number of nodes in a tree.
However, it could use the custom settings in options.complexity_mapping
if these are defined.
"""
function compute_complexity(tree::Node, options::Options)::Int
    if options.complexity_mapping.use
        return round(Int, _compute_complexity(tree, options))
    else
        return count_nodes(tree)
    end
end

function _compute_complexity(
    tree::Node, options::Options{A,B,dA,dB,C,complexity_type}
)::complexity_type where {A,B,dA,dB,C,complexity_type<:Real}
    if tree.degree == 0
        if tree.constant
            return options.complexity_mapping.constant_complexity
        else
            return options.complexity_mapping.variable_complexity
        end
    elseif tree.degree == 1
        return (
            options.complexity_mapping.unaop_complexities[tree.op] +
            _compute_complexity(tree.l, options)
        )
    else # tree.degree == 2
        return (
            options.complexity_mapping.binop_complexities[tree.op] +
            _compute_complexity(tree.l, options) +
            _compute_complexity(tree.r, options)
        )
    end
end

# Get all the constants from a tree
function get_constants(tree::Node)::AbstractVector{CONST_TYPE}
    if tree.degree == 0
        if tree.constant
            return [tree.val]
        else
            return CONST_TYPE[]
        end
    elseif tree.degree == 1
        return get_constants(tree.l)
    else
        both = [get_constants(tree.l), get_constants(tree.r)]
        return [constant for subtree in both for constant in subtree]
    end
end

# Set all the constants inside a tree
function set_constants(tree::Node, constants::AbstractVector{T}) where {T<:Real}
    if tree.degree == 0
        if tree.constant
            tree.val = convert(CONST_TYPE, constants[1])
        end
    elseif tree.degree == 1
        set_constants(tree.l, constants)
    else
        numberLeft = count_constants(tree.l)
        set_constants(tree.l, constants)
        set_constants(tree.r, constants[(numberLeft + 1):end])
    end
end

## Assign index to nodes of a tree
# This will mirror a Node struct, rather
# than adding a new attribute to Node.
mutable struct NodeIndex
    constant_index::Int  # Index of this constant (if a constant exists here)
    l::NodeIndex
    r::NodeIndex

    NodeIndex() = new()
end

function index_constants(tree::Node)::NodeIndex
    return index_constants(tree, 0)
end

function index_constants(tree::Node, left_index::Int)::NodeIndex
    index_tree = NodeIndex()
    index_constants(tree, index_tree, left_index)
    return index_tree
end

# Count how many constants to the left of this node, and put them in a tree
function index_constants(tree::Node, index_tree::NodeIndex, left_index::Int)
    if tree.degree == 0
        if tree.constant
            index_tree.constant_index = left_index + 1
        end
    elseif tree.degree == 1
        index_tree.constant_index = count_constants(tree.l)
        index_tree.l = NodeIndex()
        index_constants(tree.l, index_tree.l, left_index)
    else
        index_tree.l = NodeIndex()
        index_tree.r = NodeIndex()
        index_constants(tree.l, index_tree.l, left_index)
        index_tree.constant_index = count_constants(tree.l)
        left_index_here = left_index + index_tree.constant_index
        index_constants(tree.r, index_tree.r, left_index_here)
    end
end

# python function to check if expression is monontonically increasing
py"""
import sympy
def is_monotonic_increasing_test(expr, interval, var):

    # constant value never decreases    
    if expr.is_constant():
        return True

    # get critical points as list
    turning_points = list(sympy.solveset(expr.diff(var), var, interval))
    turning_points.sort()
    # failed to find critical points
    # there could be 0 or infinite...
    if (turning_points == []):
        # fall back to simpler increasing function
        return bool(1 if (expr.limit(var, interval.end) - expr.limit(var, interval.start)) >= 0 else 0)
    increasing = 1
    # turn to false if interval from start of main interval to first critical point not increasing
    increasing = min(increasing, (1 if (expr.limit(var, turning_points[0]) - expr.limit(var, interval.start)) >= 0 else 0))
    # check intervals between all critical points
    for i in range(len(turning_points)-1):
        thisPoint = turning_points[i]
        nextPoint = turning_points[i+1]
        increasing = min(increasing, (1 if (expr.limit(var, nextPoint) - expr.limit(var, thisPoint)) >= 0 else 0))
        #increasing = min(increasing, sympy.is_increasing(expr, sympy.Interval(thisPoint, nextPoint, false, false), var))
    # check last interval
    increasing = min(increasing, (1 if (expr.limit(var, interval.end) - expr.limit(var, turning_points[-1])) >= 0 else 0))
    #increasing = min(increasing, sympy.is_increasing(expr, sympy.Interval(turning_points[-1], interval.end, false, false), var))
    return bool(increasing)
"""

function thermoConstraints(expr::PyObject, var::PyObject)

    results = [true, true, true]
    
    # Axiom 1: the expr needs to pass through the origin
    try
        if sympy.limit(expr, var, 0, "+") != 0
            #println("constraint 1")
            results[0] = false
        end
    catch error
        #println(error)
        #println("SymPy cannot evaluate Axiom 1")
        results[0] = false
    end
        # Axiom 2: the expr needs to converge to Henry's Law at zero pressure
    try
        if (sympy.limit(sympy.diff(expr, var), var, 0) == sympy.oo 
            || sympy.limit(sympy.diff(expr, var), var, 0) == -sympy.oo 
            || sympy.limit(sympy.diff(expr, var), var, 0) == 0)
            #println("constraint 2")
            results[1] = false
        end
    catch error
        #println(error)
        #println("SymPy cannot evaluate Axiom 2")
        results[1] = false
    end

    # Axiom 3: the expr must be strictly increasing as pressure increases
    try
        # use custom function because sympy doesn't work as expected
        if not(is_monotonic_increasing_test(expr, sympy.Interval(0,sympy.oo), var))
            #println("constraint 3")
            results[2] = false
        end
    catch error
        #println("SymPy cannot evaluate Axiom 3")
        results[2] = false
    end

    return results
end

end