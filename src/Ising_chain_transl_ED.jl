
import Utilities

"""
obtain basis in momentum space for Ising chain, with momentums {k1,...}, where kn=(n-1)*2pi/L
"""
function Ising_chain_transl_basis(L)
    #label momentum state
    #identify rep states with k=0(2pi), which include rep states for other cases
    rep_states=[]
    for s=0:2^L-1
        if Utilities.is_rep_state(s,L) push!(rep_states,s) end
    end

    #basis[n] stores rep states with momentum k=2pi.n/L
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
obtain Hamiltonian matrix for transl inv Ising chain with a fixed momentum k=2pi*nk/L. The real space Hamiltonian reads
H=-\sum_i Z_i\otimes Z_{i+1} - g\sum X_i - h\sum Z_i

return hmat
"""
function Ising_chain_transl_hmat(basis,nk,L;g=-1.05,h=0.5,elemtype=Complex128)
    #construct Hamiltonian in momentum basis
    hdim=length(basis[mod(nk,L)+1])
    hmat=zeros(elemtype,hdim,hdim)
    for s_ind=1:length(basis[mod(nk,L)+1])
        s_qb=Utilities.qubits_from_state(basis[mod(nk,L)+1][s_ind],L)
        Rs=Utilities.Rψ_from_qubits(s_qb)

        # -\sum_j Z_j\otimes Z_{j+1} -h\sum_j Z_j term
        hmat[s_ind,s_ind]+=-sum(j->(1-2*s_qb[j])*(1-2*s_qb[mod(j,L)+1]),1:L)-h*sum(x->(1-2x),s_qb)

        # -g X_j term
        for j=1:L
            t_qb=copy(s_qb)
            t_qb[j]=!t_qb[j]
            Rt=Utilities.Rψ_from_qubits(t_qb)
            t,rt=Utilities.find_krep_state(Utilities.state_from_qubits(t_qb),L)
            t_ind=searchsorted(basis[mod(nk,L)+1],t)
            if length(t_ind)==0 continue end #t not consistent with nk
            t_ind=t_ind[1]
            #@show t_ind,s_ind
            hmat[t_ind,s_ind]+=-g*exp(im*2π/L*nk*rt)*sqrt(Rs/Rt)
        end
        hmat=(hmat+hmat')/2
    end

    return hmat
end
