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
function is_fib_state_obc(ψ::Int,L::Int)
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
function find_krep_state(ψ::Int,L::Int)
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


"""
Given a real space spin state |s>, obtain <s|nk> for the eigenstates {|1k>,|2k>,...}'s, where |nk>=a_nk|sk>+... is eigenstate of both Hamiltonian and momentum 
<s|nk>=1/sqrt(N_s0)*a_nk*exp(-ik*r_0)
where |s_0> is the rep state and |s>=T^{r_0}|s_0>. N_s0=L^2/Rs and T^Rs|s>=|s>

eigenstates are stored as matrix, where each col rep one state

return {<s|1k>,<s|2k>,...}
"""
function spin_state_nk_eigenstate_overlap(s::Int,eigvecs,nk,L,basis)
    Rs=Rψ_from_state(s,L)
    s0,r0=find_krep_state(s,L)
    s0_ind=searchsorted(basis[mod(nk,L)+1],s0)
    #@show s0,r0,s0_ind
    if length(s0_ind)==0 return zeros(size(eigvecs,2)) end
    s0_ind=s0_ind[1]
    return exp(-im*2pi/L*nk*r0)*sqrt(1/Rs)*eigvecs[s0_ind,:]
end


