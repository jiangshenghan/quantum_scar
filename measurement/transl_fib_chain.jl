
include("../src/utilities.jl")
include("../src/fib_chain_transl_ED.jl")

L=26
nk=0
@show L,nk
basis=fib_chain_transl_basis(L)
@show length(basis[mod(nk,L)+1])
#=
hmat=fib_chain_qs_transl_hmat(basis,nk,L)
heig=eigfact(hmat)
heigvals=heig[:values]
heigvecs=heig[:vectors]
@show norm(imag(heigvecs))
#store heig into a file
writecsv("/home/jiangsh/quantum_scar/result/transl_qs_heigvals_L=$L\_nk=$nk.txt",heigvals)
if norm(imag(heigvecs))<1e-10
    writecsv("/home/jiangsh/quantum_scar/result/transl_qs_heigvecs_L=$L\_nk=$nk.txt",real(heigvecs))
else
    writecsv("/home/jiangsh/quantum_scar/result/transl_qs_heigvecs_L=$L\_nk=$nk\_real.txt",real(heigvecs))
    writecsv("/home/jiangsh/quantum_scar/result/transl_qs_heigvecs_L=$L\_nk=$nk\_imag.txt",imag(heigvecs))
end
# =#

#read heigs from file
heigvals=readcsv("/home/jiangsh/quantum_scar/result/transl_qs_heigvals_L=$L\_nk=$nk.txt")
heigvecs=readcsv("/home/jiangsh/quantum_scar/result/transl_qs_heigvecs_L=$L\_nk=$nk.txt")
print("read eigs successfully!\n")

#measure domain wall <Z1.Z2> / single <Z1> and write to file
zz_nm=zeros(length(basis[nk+1]),length(basis[nk+1]))
z_nm=zeros(length(basis[nk+1]),length(basis[nk+1]))
for n=1:length(basis[nk+1]), m=1:length(basis[nk+1])
    for j=1:length(basis[nk+1])
        s_qb=Utilities.qubits_from_state(basis[nk+1][j],L)
        zzs=sum(i->(-2*s_qb[i]+1)*(-2*s_qb[mod(i,L)+1]+1),collect(1:L))/L
        zz_nm[n,m]+=conj(heigvecs[j,n])*heigvecs[j,m]*zzs
        zs=sum(i->(-2*s_qb[i]+1),collect(1:L))/L
        z_nm[n,m]+=conj(heigvecs[j,n])*heigvecs[j,m]*zs
    end
    if mod(n,10)==0 && m==1 @show n; end
end
writecsv("/home/jiangsh/quantum_scar/result/transl_qs_z1z2nm_L=$L\_nk=$nk.txt",zz_nm)
writecsv("/home/jiangsh/quantum_scar/result/transl_qs_z1nm_L=$L\_nk=$nk.txt",z_nm)
