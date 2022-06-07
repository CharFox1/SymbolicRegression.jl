module EvaluateEquationDerivativeModule

using LinearAlgebra
import ..CoreModule: Node, Options
import ..UtilsModule: @return_on_false2, is_bad_array
import ..EquationUtilsModule: count_constants, index_constants, NodeIndex
import ..EvaluateEquationModule: deg0_eval

"""
    eval_diff_tree_array(tree::Node, cX::AbstractMatrix{T}, options::Options, direction::Int)

Compute the forward derivative of an expression, using a similar
structure and optimization to eval_tree_array. `direction` is the index of a particular
variable in the expression. e.g., `direction=1` would indicate derivative with
respect to `x1`.

# Returns

- `(evaluation, derivative, complete)::Tuple{AbstractVector{T}, AbstractVector{T}, Bool}`: the normal evaluation,
    the derivative, and whether the evaluation completed as normal (or encountered a nan or inf).
"""
function eval_diff_tree_array(
    tree::Node, cX::AbstractMatrix{T}, options::Options, direction::Int
)::Tuple{AbstractVector{T},AbstractVector{T},Bool} where {T<:Real}
    # TODO: Implement quick check for whether the variable is actually used
    # in this tree. Otherwise, return zero.
    evaluation, derivative, complete = _eval_diff_tree_array(tree, cX, options, direction)
    @return_on_false2 complete evaluation derivative
    return evaluation, derivative, !(is_bad_array(evaluation) || is_bad_array(derivative))
end

function _eval_diff_tree_array(
    tree::Node, cX::AbstractMatrix{T}, options::Options, direction::Int
)::Tuple{AbstractVector{T},AbstractVector{T},Bool} where {T<:Real}
    if tree.degree == 0
        diff_deg0_eval(tree, cX, options, direction)
    elseif tree.degree == 1
        diff_deg1_eval(tree, cX, Val(tree.op), options, direction)
    else
        diff_deg2_eval(tree, cX, Val(tree.op), options, direction)
    end
end

function diff_deg0_eval(
    tree::Node, cX::AbstractMatrix{T}, options::Options, direction::Int
)::Tuple{AbstractVector{T},AbstractVector{T},Bool} where {T<:Real}
    n = size(cX, 2)
    const_part = deg0_eval(tree, cX, options)[1]
    derivative_part = (tree.feature == direction) ? ones(T, n) : zeros(T, n)
    return (const_part, derivative_part, true)
end

function diff_deg1_eval(
    tree::Node, cX::AbstractMatrix{T}, ::Val{op_idx}, options::Options, direction::Int
)::Tuple{AbstractVector{T},AbstractVector{T},Bool} where {T<:Real,op_idx}
    n = size(cX, 2)
    (cumulator, dcumulator, complete) = eval_diff_tree_array(tree.l, cX, options, direction)
    @return_on_false2 complete cumulator dcumulator

    op = options.unaops[op_idx]
    diff_op = options.diff_unaops[op_idx]

    @inbounds @simd for j in 1:n
        x = op(cumulator[j])::T
        dx = diff_op(cumulator[j])::T * dcumulator[j]

        cumulator[j] = x
        dcumulator[j] = dx
    end
    return (cumulator, dcumulator, true)
end

function diff_deg2_eval(
    tree::Node, cX::AbstractMatrix{T}, ::Val{op_idx}, options::Options, direction::Int
)::Tuple{AbstractVector{T},AbstractVector{T},Bool} where {T<:Real,op_idx}
    n = size(cX, 2)
    (cumulator, dcumulator, complete) = eval_diff_tree_array(tree.l, cX, options, direction)
    @return_on_false2 complete cumulator dcumulator
    (array2, dcumulator2, complete2) = eval_diff_tree_array(tree.r, cX, options, direction)
    @return_on_false2 complete2 array2 dcumulator2

    op = options.binops[op_idx]
    diff_op = options.diff_binops[op_idx]

    @inbounds @simd for j in 1:n
        x = op(cumulator[j], array2[j])

        dx = dot(diff_op(cumulator[j], array2[j]), [dcumulator[j], dcumulator2[j]])

        cumulator[j] = x
        dcumulator[j] = dx
    end
    return (cumulator, dcumulator, true)
end

"""
    eval_grad_tree_array(tree::Node, cX::AbstractMatrix{T}, options::Options; variable::Bool=false)

Compute the forward-mode derivative of an expression, using a similar
structure and optimization to eval_tree_array. `variable` specifies whether
we should take derivatives with respect to features (i.e., cX), or with respect
to every constant in the expression.

# Returns

- `(evaluation, gradient, complete)::Tuple{AbstractVector{T}, AbstractMatrix{T}, Bool}`: the normal evaluation,
    the gradient, and whether the evaluation completed as normal (or encountered a nan or inf).
"""
function eval_grad_tree_array(
    tree::Node, cX::AbstractMatrix{T}, options::Options; variable::Bool=false
)::Tuple{AbstractVector{T},AbstractMatrix{T},Bool} where {T<:Real}
    n = size(cX, 2)
    if variable
        n_gradients = size(cX, 1)
    else
        n_gradients = count_constants(tree)
    end
    index_tree = index_constants(tree, 0)
    return eval_grad_tree_array(
        tree, n, n_gradients, index_tree, cX, options, Val(variable)
    )
end

function eval_grad_tree_array(
    tree::Node,
    n::Int,
    n_gradients::Int,
    index_tree::NodeIndex,
    cX::AbstractMatrix{T},
    options::Options,
    ::Val{variable},
)::Tuple{AbstractVector{T},AbstractMatrix{T},Bool} where {T<:Real,variable}
    evaluation, gradient, complete = _eval_grad_tree_array(
        tree, n, n_gradients, index_tree, cX, options, Val(variable)
    )
    @return_on_false2 complete evaluation gradient
    return evaluation, gradient, !(is_bad_array(evaluation) || is_bad_array(gradient))
end

function _eval_grad_tree_array(
    tree::Node,
    n::Int,
    n_gradients::Int,
    index_tree::NodeIndex,
    cX::AbstractMatrix{T},
    options::Options,
    ::Val{variable},
)::Tuple{AbstractVector{T},AbstractMatrix{T},Bool} where {T<:Real,variable}
    if tree.degree == 0
        grad_deg0_eval(tree, n, n_gradients, index_tree, cX, options, Val(variable))
    elseif tree.degree == 1
        grad_deg1_eval(
            tree, n, n_gradients, index_tree, cX, Val(tree.op), options, Val(variable)
        )
    else
        grad_deg2_eval(
            tree, n, n_gradients, index_tree, cX, Val(tree.op), options, Val(variable)
        )
    end
end

function grad_deg0_eval(
    tree::Node,
    n::Int,
    n_gradients::Int,
    index_tree::NodeIndex,
    cX::AbstractMatrix{T},
    options::Options,
    ::Val{variable},
)::Tuple{AbstractVector{T},AbstractMatrix{T},Bool} where {T<:Real,variable}
    const_part = deg0_eval(tree, cX, options)[1]

    if variable == tree.constant
        return (const_part, zeros(T, n_gradients, n), true)
    end

    index = variable ? tree.feature : index_tree.constant_index
    derivative_part = zeros(T, n_gradients, n)
    derivative_part[index, :] .= T(1)
    return (const_part, derivative_part, true)
end

function grad_deg1_eval(
    tree::Node,
    n::Int,
    n_gradients::Int,
    index_tree::NodeIndex,
    cX::AbstractMatrix{T},
    ::Val{op_idx},
    options::Options,
    ::Val{variable},
)::Tuple{AbstractVector{T},AbstractMatrix{T},Bool} where {T<:Real,op_idx,variable}
    (cumulator, dcumulator, complete) = eval_grad_tree_array(
        tree.l, n, n_gradients, index_tree.l, cX, options, Val(variable)
    )
    @return_on_false2 complete cumulator dcumulator

    op = options.unaops[op_idx]
    diff_op = options.diff_unaops[op_idx]

    @inbounds @simd for j in 1:n
        x = op(cumulator[j])::T
        dx = diff_op(cumulator[j])

        cumulator[j] = x
        for k in 1:n_gradients
            dcumulator[k, j] = dx * dcumulator[k, j]
        end
    end
    return (cumulator, dcumulator, true)
end

function grad_deg2_eval(
    tree::Node,
    n::Int,
    n_gradients::Int,
    index_tree::NodeIndex,
    cX::AbstractMatrix{T},
    ::Val{op_idx},
    options::Options,
    ::Val{variable},
)::Tuple{AbstractVector{T},AbstractMatrix{T},Bool} where {T<:Real,op_idx,variable}
    derivative_part = Array{T,2}(undef, n_gradients, n)
    (cumulator1, dcumulator1, complete) = eval_grad_tree_array(
        tree.l, n, n_gradients, index_tree.l, cX, options, Val(variable)
    )
    @return_on_false2 complete cumulator1 dcumulator1
    (cumulator2, dcumulator2, complete2) = eval_grad_tree_array(
        tree.r, n, n_gradients, index_tree.r, cX, options, Val(variable)
    )
    @return_on_false2 complete2 cumulator1 dcumulator1

    op = options.binops[op_idx]
    diff_op = options.diff_binops[op_idx]

    @inbounds @simd for j in 1:n
        x = op(cumulator1[j], cumulator2[j])
        dx = diff_op(cumulator1[j], cumulator2[j])
        cumulator1[j] = x
        for k in 1:n_gradients
            derivative_part[k, j] = dx[1] * dcumulator1[k, j] + dx[2] * dcumulator2[k, j]
        end
    end

    return (cumulator1, derivative_part, true)
end

end
