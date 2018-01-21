#measure observables using tebd method. Fib chain Hilbert space are imposed by soft constraint term
# H=\sum_j X_j + V \sum_j 1/4(I_j-Z_j).(I_{j+1}-Z_{j+1}), with infinite V limit

include("/home/jiangsh/JTensor.jl/src/JTensor.jl")
#include("/mnt/c/Users/sheng/OneDrive\ -\ California\ Institute\ of\ Technology/Code/JTensor.jl/src/JTensor.jl")
using JTensor

d=2
L=24
tf=100
δt=0.02
chi=200
@show L,tf,δt,chi

I=[1 0; 0 1]
X=[0 1; 1 0]
Y=[0 -im; im 0]
Z=[1 0; 0 -1]

#construct Hamiltonian
V=200
@show V
H=[]
push!(H,0.25*V*jcontract([I-Z,I-Z],[[-1,-3],[-2,-4]])+jcontract([X,I],[[-1,-3],[-2,-4]])+0.5*jcontract([I,X],[[-1,-3],[-2,-4]]))
for j=2:L-2
    push!(H,0.25*V*jcontract([I-Z,I-Z],[[-1,-3],[-2,-4]])+0.5*jcontract([X,I],[[-1,-3],[-2,-4]])+0.5*jcontract([I,X],[[-1,-3],[-2,-4]]))
end
push!(H,0.25*V*jcontract([I-Z,I-Z],[[-1,-3],[-2,-4]])+0.5*jcontract([X,I],[[-1,-3],[-2,-4]])+jcontract([I,X],[[-1,-3],[-2,-4]]))
Hmat=[]
for j=1:L-1
    push!(Hmat,reshape(H[j],d^2,d^2))
end

#construct time evolution operator
Ue=[]
for k=1:div(L-1,2)
    j=2k
    Heig=eigfact(Hmat[j])
    push!(Ue,Heig[:vectors]*diagm(exp.(-im*δt*Heig[:values]))*Heig[:vectors]')
    Ue[k]=reshape(transpose(Ue[k]),d,d,d,d)
end
Uo=[]
for k=1:div(L,2)
    j=2k-1
    Heig=eigfact(Hmat[j])
    push!(Uo,Heig[:vectors]*diagm(exp.(-im*0.5δt*Heig[:values]))*Heig[:vectors]')
    Uo[k]=reshape(transpose(Uo[k]),d,d,d,d)
end

#initilize state with both site and bond tensor
#|g>=|000...>
A=[]

for j=1:L
    push!(A,zeros(1,1,d))
    A[j][1]=1
end
B=[]
for j=1:L-1 push!(B,ones(1,1)) end

#tebd
t=0
cind=div(L,2)
while t<tf
    #The measurement assumes canonical form for MPS
    #zzt=jcontract([A[1],B[1],A[2],B[2],Z,Z,conj(A[1]),conj(B[1]),conj(A[2]),conj(B[2])],[[1,2,6],[2,3],[3,4,7],[4,5],[6,8],[7,9],[1,10,8],[10,11],[11,12,9],[12,5]])
    zzt=jcontract([B[cind-1],A[cind],B[cind],A[cind+1],B[cind+1],Z,Z,conj(B[cind-1]),conj(A[cind]),conj(B[cind]),conj(A[cind+1]),conj(B[cind+1])],[[1,2],[2,3,7],[3,4],[4,5,8],[5,6],[7,9],[8,10],[1,11],[11,12,9],[12,13],[13,14,10],[14,6]])

    #measure energy of tebd, which serves as a criteria for accuracy
    if abs(t-round(t))<1e-7
        wf_norm=jcontract([A[1],conj(A[1]),B[1],conj(B[1])],[[1,3,2],[1,4,2],[3,-1],[4,-2]])
        for j=2:L-1 
            wf_norm=jcontract([wf_norm,A[j],conj(A[j]),B[j],conj(B[j])],[[1,2],[1,4,3],[2,5,3],[4,-1],[5,-2]])
        end
        wf_norm=jcontract([wf_norm,A[L],conj(A[L])],[[1,2],[1,4,3],[2,4,3]])

        bond_energy=zeros(Complex128,L-1)
        bond_energy[1]=jcontract([A[1],B[1],A[2],B[2],H[1],conj(A[1]),conj(B[1]),conj(A[2]),conj(B[2])],[[1,2,5],[2,3],[3,4,6],[4,9],[5,6,7,8],[1,10,7],[10,11],[11,12,8],[12,9]])
        for j=2:L-2
            bond_energy[j]=jcontract([B[j-1],A[j],B[j],A[j+1],B[j+1],H[j],conj(B[j-1]),conj(A[j]),conj(B[j]),conj(A[j+1]),conj(B[j+1])],[[1,2],[2,3,7],[3,4],[4,5,8],[5,6],[7,8,9,10],[1,11],[11,12,9],[12,13],[13,14,10],[14,6]])
        end
        bond_energy[end]=jcontract([B[L-2],A[L-1],B[L-1],A[L],H[L-1],conj(B[L-2]),conj(A[L-1]),conj(B[L-1]),conj(A[L])],[[1,2],[2,3,6],[3,4],[4,5,7],[6,7,8,9],[1,10],[10,11,8],[11,12],[12,5,9]])

        energy=sum(bond_energy)
        @show t,wf_norm,energy
    end

    @show t,zzt
    A,B=tebd_even_odd_one_step(A,B,Ue,Uo,eps=1e-6,chi=chi)
    t=t+δt
end
