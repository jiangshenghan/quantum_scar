module Utilities

fib(n::Int)=n<3?1:fib(n-1)+fib(n-2)

"""
transform from a spin-1/2 state, which stored as int, to binary qubits, which stored as array of Bool
e.g. for L=6, 11 --> 1011 --> [1,1,0,1,0,0] ~ |dduduu>

return ψ_qb (BitArray)
"""
function qubits_from_state(ψ::Int,L::Int)
    ψ_qb=BitArray(map(x->parse(Bool,x),collect(reverse(bin(ψ)))))
    append!(ψ_qb,zeros(Bool,L-length(ψ_qb)))
    return ψ_qb
end

"""
transform from binary qubit, which stored as Bool array or string, to a spin-1/2 chain state, which represented as int
e.g. |dduduu> ~ [1,1,0,1,0,0]/"110100" --> 001011 --> 11

return state as int
"""
function state_from_qubits(ψ_qb::BitArray)
    return parse(Int,reverse(*(bin.(ψ_qb)...)),2)
end
function state_from_qubits(ψ_qb::String)
    return parse(Int,reverse(ψ_qb),2)
end

"""
obtain the translation little group for a particular state ψ
"""
function Rψ_from_qubits(ψ_qb::BitArray)
    φ_qb=copy(ψ_qb)
    L=length(ψ_qb)
    j=0
    for j=1:L
        φ_qb=rol(φ_qb,-1)
        if φ_qb==ψ_qb break end
    end
    return j
end
Rψ_from_state(ψ::Int,L::Int)=Rψ_from_qubits(qubits_from_state(ψ,L))

"""
check if the spin state is allowed in Fib chain Hilbert space
PBC is assumed
"""
function is_fib_state(ψ::Int,L::Int)
    ψ_qb=qubits_from_state(ψ,L)
    #@show Array{Int}(ψ_qb)
    for j=1:L
        if ψ_qb[j]==1 && ψ_qb[mod(j,L)+1]==1
            return false
        end
    end
    return true
end

"""
check if the spin state is allowed in Fib chain Hilbert space for open chain
"""
function is_fib_chain_obc()
    ψ_qb=qubits_from_state(ψ,L)
    for j=1:L-1
        if ψ_qb[j]==1 && ψ_qb[j+1]==1
            return false
        end
    end
    return true
end

"""
check if the real space state is the represent state of a momentum state. Namely, check if the real space state have the smallest integer among all states related by translate this state.
"""
function is_rep_state(ψ::Int,L::Int)
    ψ_qb=qubits_from_state(ψ,L)
    φ_qb=copy(ψ_qb)
    for j=1:L-1
        φ_qb=rol(φ_qb,-1)
        φ=state_from_qubits(φ_qb)
        #@show φ,Array{Int}(φ_qb)
        if φ<ψ return false end
        if φ==ψ return true end
    end
    return true
end

"""
given a real space state, find the corresponding rep state
return the rep state and the step to rep state
|ψ>=T^{r}|rep>

return |rep> (stored as int), r (distance from ref state)
"""
function find_rep_state(ψ::Int,L::Int)
    rep_state=ψ
    φ_qb=qubits_from_state(ψ,L)
    r=0
    for j=1:L-1
        #perform T^{-1}
        φ_qb=rol(φ_qb,1)
        φ=state_from_qubits(φ_qb)
        if φ<rep_state rep_state=φ; r=j end
        if φ==ψ break end
    end
    return rep_state,r
end

end
