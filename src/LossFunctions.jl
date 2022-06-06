using FromFile
using Random: randperm
using LossFunctions
@from "Core.jl" import Options, Dataset, Node
@from "Equation.jl" import stringTree
@from "EquationUtils.jl" import countNodes, thermoConstraints
@from "EvaluateEquation.jl" import evalTreeArray, differentiableEvalTreeArray

# for symbolic constraints
using PyCall
const _sympy = PyNULL() 

# test init function so sympy can load
#function __init__()
    
#end

#using PyCall
#sympy = pyimport("sympy")

function Loss(x::AbstractArray{T}, y::AbstractArray{T}, options::Options{A,B,C})::T where {T<:Real,A,B,C<:SupervisedLoss}
    value(options.loss, y, x, AggMode.Mean())
end
function Loss(x::AbstractArray{T}, y::AbstractArray{T}, options::Options{A,B,C})::T where {T<:Real,A,B,C<:Function}
    sum(options.loss.(x, y))/length(y)
end

function Loss(x::AbstractArray{T}, y::AbstractArray{T}, w::AbstractArray{T}, options::Options{A,B,C})::T where {T<:Real,A,B,C<:SupervisedLoss}
    value(options.loss, y, x, AggMode.WeightedMean(w))
end
function Loss(x::AbstractArray{T}, y::AbstractArray{T}, w::AbstractArray{T}, options::Options{A,B,C})::T where {T<:Real,A,B,C<:Function}
    sum(options.loss.(x, y, w))/sum(w)
end

# Loss function. Only MSE implemented right now. TODO
# Also need to put actual loss function in scoreFuncBatch!
function EvalLoss(tree::Node, dataset::Dataset{T}, options::Options;
                  allow_diff=false)::T where {T<:Real}
    if !allow_diff
        (prediction, completion) = evalTreeArray(tree, dataset.X, options)
    else
        (prediction, completion) = differentiableEvalTreeArray(tree, dataset.X, options)
    end
    if !completion
        return T(1000000000)
    end

    if dataset.weighted
        return Loss(prediction, dataset.y, dataset.weights, options)
    else
        return Loss(prediction, dataset.y, options)
    end
end

# Score an equation
function scoreFunc(dataset::Dataset{T},
                   baseline::T, tree::Node,
                   options::Options; allow_diff=false)::T where {T<:Real}
    mse = EvalLoss(tree, dataset, options; allow_diff=allow_diff)

    # try to pass python libary through options
    #if !haskey(options.pyLibs, "sympy")
    #    copy!(sympy, pyimport("sympy"))
    #    options.pyLibs["sympy"] = sympy
    #end

    #sympyVar = options.pyLibs["sympy"]

    if options.penalties !== nothing
        if ispynull(_sympy)
            println("loading lib")
        end
        # good way to check if lib needs to be loaded from PyCall/serialize.jl
        sympy = ispynull(_sympy) ? copy!(_sympy, pyimport_conda("sympy", "sympy")) : _sympy
        
        if length(dataset.varMap) >= 1 && ("p" in dataset.varMap)
            # turn tree to string for sympy to use
            # get boolean result for each constraint
            str = stringTree(tree, options)
            expr = sympy.parse_expr(str)
            #expr = sympyVar.parse_expr("p*2 + 1")
            var = sympy.symbols("p")
            results = thermoConstraints(expr, var)
            #results = [true, true, true]
            # for each constraint, apply penalty if expression fails
            for i in 1:length(options.penalties)
                if results[i] == false
                    mse *= options.penalties[i]
                end
            end
        else
            println("Need pressure variable \"p\"!")
        end
    end

    # generic symbolic constraint code (not working yet)
    #if options.sym_constraints !== nothing
    #    # check symbolic constraints
    #    for item in options.sym_constraints
    #        # if constraint is not passed, apply penalty
    #        result = iszero(item[2](dataset, baseline, tree, options::Options, allow_diff=allow_diff))
    #        if result
    #            mse = mse * item[1] # scale by penalty 
    #        end
    #    end 
    #end

    return mse / baseline + countNodes(tree)*options.parsimony
end

# Score an equation with a small batch
function scoreFuncBatch(dataset::Dataset{T}, baseline::T,
                        tree::Node, options::Options)::T where {T<:Real}
    batchSize = options.batchSize
    batch_idx = randperm(dataset.n)[1:options.batchSize]
    batch_X = dataset.X[:, batch_idx]
    batch_y = dataset.y[batch_idx]
    (prediction, completion) = evalTreeArray(tree, batch_X, options)
    if !completion
        return T(1000000000)
    end

    if !dataset.weighted
        mse = Loss(prediction, batch_y, options)
    else
        batch_w = dataset.weights[batch_idx]
        mse = Loss(prediction, batch_y, batch_w, options)
    end
    return mse / baseline + countNodes(tree) * options.parsimony
end
