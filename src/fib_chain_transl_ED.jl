
import Utilities

"""
obtain basis in momentum space for fib chain, with momentums {k1,...}, where kn=(n-1)*2pi/L

return basis
"""
function fib_chain_transl_basis(L)
    #label momentum state
    #identify rep states with k=0(2pi), which include rep states for other cases
    rep_states=[]
    for s=0:2^L-1
        if Utilities.is_fib_state(s,L)==false continue end
        if Utilities.is_rep_state(s,L) push!(rep_states,s) end
    end

    #basis[n] stores rep states with momentum k=2pi.mks[n]/L, mks given by users (default value 0,...L-1)
    basis=[]
    for n=1:L
        nk_basis=[]
        for s in rep_states
            Rs=Utilities.Rψ_from_state(s,L)
            if mod(n-1,div(L,Rs))==0 
                push!(nk_basis,s)
            end
        end
        push!(basis,nk_basis)
    end

    return basis
end

"""
obtain quantum scar Hamiltonian matrix for transl inv fib chain with a fixed momentum k=2pi*nk/L. The real space Hamiltonian reads
H=\sum P_{i-1}X_i P_{i+1}

return hmat
"""
function fib_chain_qs_transl_hmat(basis,nk,L)
    #construct Hamiltonian in momentum basis
    hdim=length(basis[mod(nk,L)+1])
    hmat=zeros(Complex128,hdim,hdim)
    for s_ind=1:length(basis[mod(nk,L)+1])
        s_qb=Utilities.qubits_from_state(basis[mod(nk,L)+1][s_ind],L)
        Rs=Utilities.Rψ_from_qubits(s_qb)
        for j=1:L
            if s_qb[mod(j-2,L)+1]==1 || s_qb[mod(j,L)+1]==1 continue end
            t_qb=copy(s_qb)
            t_qb[j]=!t_qb[j]
            Rt=Utilities.Rψ_from_qubits(t_qb)
            t,rt=Utilities.find_rep_state(Utilities.state_from_qubits(t_qb),L)
            t_ind=searchsorted(basis[mod(nk,L)+1],t)
            if length(t_ind)==0 continue end #t not consistent with mks[n]
            t_ind=t_ind[1]
            #@show t_ind,s_ind
            hmat[t_ind,s_ind]+=exp(im*2π/L*nk*rt)*sqrt(Rs/Rt)
        end
        hmat=(hmat+hmat')/2
    end

    return hmat
end


"""
obtain XZZ Hamiltonian matrix for transl inv fib chain with a fixed momentum k=2pi*nk/L. The real space Hamiltonian reads
H=\sum P_{i-1}X_i P_{i+1} + g\sum Z_i Z_{i+1}

return hmat
"""
function fib_chain_XZZ_transl_hmat(basis,nk,L;g=0.3)
    #construct Hamiltonian in momentum basis
    hdim=length(basis[mod(nk,L)+1])
    hmat=zeros(Complex128,hdim,hdim)
    for s_ind=1:length(basis[mod(nk,L)+1])
        s_qb=Utilities.qubits_from_state(basis[mod(nk,L)+1][s_ind],L)
        Rs=Utilities.Rψ_from_qubits(s_qb)
        hmat[s_ind,s_ind]+=g*sum(i->(-2*s_qb[i]+1)*(-2*s_qb[mod(i,L)+1]+1),collect(1:L))
        for j=1:L
            if s_qb[mod(j-2,L)+1]==1 || s_qb[mod(j,L)+1]==1 continue end
            t_qb=copy(s_qb)
            t_qb[j]=!t_qb[j]
            Rt=Utilities.Rψ_from_qubits(t_qb)
            t,rt=Utilities.find_rep_state(Utilities.state_from_qubits(t_qb),L)
            t_ind=searchsorted(basis[mod(nk,L)+1],t)
            if length(t_ind)==0 continue end #t not consistent with mks[n]
            t_ind=t_ind[1]
            #@show t_ind,s_ind
            hmat[t_ind,s_ind]+=exp(im*2π/L*nk*rt)*sqrt(Rs/Rt)
        end
        hmat=(hmat+hmat')/2
    end

    return hmat
end


"""
Given a spin state |s>, obtain <s|nk> for the eigenstates {|1k>,|2k>,...}'s, where |nk>=a_nk|sk>+...
<s|nk>=1/sqrt(N_s0)*a_nk*exp(-ik*r_0)
where |s_0> is the rep state and |s>=T^{r_0}|s_0>. N_s0=L^2/Rs and T^Rs|s>=|s>

eigenstates are stored as matrix, where each col rep one state

return {<s|1k>,<s|2k>,...}
"""
function spin_state_overlap(s::Int,eigvecs,nk,L,basis)
    Rs=Utilities.Rψ_from_state(s,L)
    s0,r0=Utilities.find_rep_state(s,L)
    s0_ind=searchsorted(basis[mod(nk,L)+1],s0)
    #@show s0,r0,s0_ind
    if length(s0_ind)==0 return zeros(size(eigvecs,2)) end
    s0_ind=s0_ind[1]
    return exp(-im*2pi/L*nk*r0)*sqrt(1/Rs)*eigvecs[s0_ind,:]
end


