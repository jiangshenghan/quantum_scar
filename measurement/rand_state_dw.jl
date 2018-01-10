#measure <Z1.Z2> for a random state

include("../src/utilities.jl")
include("../src/fib_chain_ED.jl")

L=20
basis=fib_chain_basis(L)

#=
hmat=fib_chain_qs_hmat(basis,L)
heig=eigfact(hmat)
heigvals=heig[:values]
heigvecs=heig[:vectors]
#store heig into a file
writedlm("/home/jiangsh/quantum_scar/result/qs_heigvals_L=$L.txt",heig[:values])
writedlm("/home/jiangsh/quantum_scar/result/qs_heigvecs_L=$L.txt",heig[:vectors])
# =#

#read heigs from file
heigvals=readdlm("/home/jiangsh/quantum_scar/result/qs_heigvals_L=$L.txt")
heigvecs=readdlm("/home/jiangsh/quantum_scar/result/qs_heigvecs_L=$L.txt")
print("successfully read eigs!\n")

@show length(heigvals)

#measure domain wall <Z1.Z2> / single <Z1> and write to file
zz_nm=zeros(length(basis),length(basis))
z_nm=zeros(length(basis),length(basis))
for n=1:length(basis), m=1:length(basis)
    zz_nm[n,m]=sum(1:length(basis)) do j
        s_qb=Utilities.qubits_from_state(basis[j],L)
        return conj(heigvecs[j,n])*heigvecs[j,m]*(-2*s_qb[1]+1)*(-2*s_qb[2]+1)
    end
    z_nm[n,m]=sum(1:length(basis)) do j
        zj=(-2*Utilities.qubits_from_state(basis[j],L)[1]+1)
        return conj(heigvecs[j,n])*heigvecs[j,m]*zj
    end
    if mod(n,10)==0 && m==1 @show n; end
end
writedlm("/home/jiangsh/quantum_scar/result/qs_z1z2nm_L=$L.txt",zz_nm)
writedlm("/home/jiangsh/quantum_scar/result/qs_z1nm_L=$L.txt",z_nm)
