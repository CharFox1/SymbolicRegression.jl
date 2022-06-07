using FromFile
@from "Core.jl" import CONST_TYPE, Node, copyNode, Options

# for symbolic constraints
using PyCall
#sympy = pyimport("sympy")

# Count the operators, constants, variables in an equation
function countNodes(tree::Node)::Int
    if tree.degree == 0
        return 1
    elseif tree.degree == 1
        return 1 + countNodes(tree.l)
    else
        return 1 + countNodes(tree.l) + countNodes(tree.r)
    end
end

# Count the max depth of a tree
function countDepth(tree::Node)::Int
    if tree.degree == 0
        return 1
    elseif tree.degree == 1
        return 1 + countDepth(tree.l)
    else
        return 1 + max(countDepth(tree.l), countDepth(tree.r))
    end
end

# Count the number of unary operators in the equation
function countUnaryOperators(tree::Node)::Int
    if tree.degree == 0
        return 0
    elseif tree.degree == 1
        return 1 + countUnaryOperators(tree.l)
    else
        return 0 + countUnaryOperators(tree.l) + countUnaryOperators(tree.r)
    end
end

# Count the number of binary operators in the equation
function countBinaryOperators(tree::Node)::Int
    if tree.degree == 0
        return 0
    elseif tree.degree == 1
        return 0 + countBinaryOperators(tree.l)
    else
        return 1 + countBinaryOperators(tree.l) + countBinaryOperators(tree.r)
    end
end

# Count the number of operators in the equation
function countOperators(tree::Node)::Int
    return countUnaryOperators(tree) + countBinaryOperators(tree)
end

# Count the number of constants in an equation
function countConstants(tree::Node)::Int
    if tree.degree == 0
        if tree.constant
            return 1
        else
            return 0
        end
    elseif tree.degree == 1
        return 0 + countConstants(tree.l)
    else
        return 0 + countConstants(tree.l) + countConstants(tree.r)
    end
end

# Get all the constants from a tree
function getConstants(tree::Node)::AbstractVector{CONST_TYPE}
    if tree.degree == 0
        if tree.constant
            return [tree.val]
        else
            return CONST_TYPE[]
        end
    elseif tree.degree == 1
        return getConstants(tree.l)
    else
        both = [getConstants(tree.l), getConstants(tree.r)]
        return [constant for subtree in both for constant in subtree]
    end
end

# Set all the constants inside a tree
function setConstants(tree::Node, constants::AbstractVector{T}) where {T<:Real}
    if tree.degree == 0
        if tree.constant
            tree.val = convert(CONST_TYPE, constants[1])
        end
    elseif tree.degree == 1
        setConstants(tree.l, constants)
    else
        numberLeft = countConstants(tree.l)
        setConstants(tree.l, constants)
        setConstants(tree.r, constants[(numberLeft + 1):end])
    end
end

# custom function to check if expression is monontonically increasing
py"""
import sympy
def is_monotonic_increasing(expr, interval, var):
    # get critical points as list
    turning_points = list(sympy.solveset(expr.diff(var), var, interval))
    turning_points.sort()
    print("turning points =", turning_points)

    # failed to find critical points
    # there could be 0 or infinite...
    if (turning_points == []):
        # fall back to simpler increasing function
        return sympy.is_increasing(expr, sympy.Interval(0,sympy.oo), var)

    increasing = 1

    #print("turning points = ", turning_points)
    # turn to false if interval from start of main interval to first critical point not increasing
    increasing = min(increasing, sympy.is_increasing(expr, sympy.Interval(interval.start, turning_points[0]), var))
    # check intervals between all critical points
    for i in range(len(turning_points)-1):
        thisPoint = turning_points[i]
        nextPoint = turning_points[i+1]
        increasing = min(increasing, sympy.is_increasing(expr, sympy.Interval(thisPoint, nextPoint), var))
    # check last interval
    increasing = min(increasing, sympy.is_increasing(expr, sympy.Interval(turning_points[-1], interval.end), var))
    return bool(increasing)
"""

# no parameters to fill variant
function thermoConstraints(expr::PyObject, var::PyObject)

    # don't do this here :(
    #sympy = pyimport("sympy")
    #var = sympy.symbols("p")
    #expr = sympy.parse_expr(expr)

    #println("expression:", expr)
    results = [true, true, true]

    # Axiom 1: the expr needs to pass through the origin
    try
        if sympy.limit(expr, var, 0, "+") != 0
            #println("constraint 1")
            results[1] = false
        end
    catch error
        #println(error)
        #println("SymPy cannot evaluate Axiom 1")
        results[1] = false
    end
    # Axiom 2: the expr needs to converge to Henry's Law at zero pressure
    try
        if (
            sympy.limit(sympy.diff(expr, var), var, 0) == sympy.oo ||
            sympy.limit(sympy.diff(expr, var), var, 0) == -sympy.oo ||
            sympy.limit(sympy.diff(expr, var), var, 0) == 0
        )
            #println("constraint 2")
            results[2] = false
        end
    catch error
        #println(error)
        #println("SymPy cannot evaluate Axiom 2")
        results[2] = false
    end

    # Axiom 3: the expr must be strictly increasing as pressure increases
    try
        # use custom function because sympy doesn't work as expected
        if !(py"is_monotonic_increasing"(expr, sympy.Interval(0, sympy.oo), var))
            #println("constraint 3")
            results[3] = false
        end
    catch error
        #println("SymPy cannot evaluate Axiom 3")
        #println(error)
        results[3] = false
    end

    return results
end
