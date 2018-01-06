
import Utilities

"""
obtain basis for fib chain Hilbert space, no neighbour 11 is allowed

return basis
"""
function fib_chain_basis(L)
    #label states
    ind=0
    basis=Array{Int,1}()
    for s=0:2^L-1
        if Utilities.is_fib_state(s,L)
            ind+=1
            push!(basis,s)
        end
    end
    return basis
end


"""
obtain the fib chain quantum scar Hamitonian matrix 
H=\sum P_{i-1}X_i P_{i+1}
PBC is assumed

return hmat
"""
function fib_chain_qs_hmat(basis,L)
    hdim=length(basis)
    #construct Hamiltonian
    hmat=zeros(hdim,hdim)
    for ind=1:hdim
        ψ_qb=Utilities.qubits_from_state(basis[ind],L)
        #Hj=P_{j-1}X_jP_{j+1}
        for j=1:L
            if ψ_qb[mod(j-2,L)+1]==1 || ψ_qb[mod(j,L)+1]==1 continue end
            φ_qb=copy(ψ_qb)
            φ_qb[j]=!φ_qb[j]
            φ=Utilities.state_from_qubits(φ_qb)
            φ_ind=searchsorted(basis,φ)
            if length(φ_ind)>0 φ_ind=φ_ind[1] end
            hmat[φ_ind,ind]+=1
        end
    end

    return hmat
end

"""
obtain random Hamiltonian for fib chain

return hmat
"""
function fib_chain_rand_hmat(basis,L)
    hdim=length(basis)
    hmat=zeros(hdim,hdim)
    for j=1:hdim,k=j+1:hdim
        hmat[j,k]=2*rand()-1 #matrix element from -1 to 1
    end
    hmat=hmat+hmat'
    foreach(j->hmat[j,j]=rand(),1:hdim)

    return hmat
end

