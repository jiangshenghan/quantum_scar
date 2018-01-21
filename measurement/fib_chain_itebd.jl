#measure observables using itebd method. Fib chain Hilbert space are imposed by soft constraint term
# H=\sum_j X_j + V \sum_j 1/4(I_j-Z_j).(I_{j+1}-Z_{j+1}), with infinite V limit

include("/home/jiangsh/JTensor.jl/src/JTensor.jl")
#include("/mnt/c/Users/sheng/OneDrive\ -\ California\ Institute\ of\ Technology/Code/JTensor.jl/src/JTensor.jl")
using JTensor

d=2
tf=100
δt=0.02
chi=200
@show tf,δt,chi

I=[1 0; 0 1]
X=[0 1; 1 0]
Y=[0 -im; im 0]
Z=[1 0; 0 -1]

#construct Hamiltonian
V=400
@show V
H=0.25*V*jcontract([I-Z,I-Z],[[-1,-3],[-2,-4]])+0.5*jcontract([X,I],[[-1,-3],[-2,-4]])+0.5*jcontract([I,X],[[-1,-3],[-2,-4]])
Hmat=reshape(H,d^2,d^2)

#construct time evolution operator
Heig=eigfact(Hmat)
Ue=Heig[:vectors]*diagm(exp.(-im*δt*Heig[:values]))*Heig[:vectors]'
Ue=reshape(transpose(Ue),d,d,d,d)
Uo=Heig[:vectors]*diagm(exp.(-im*0.5δt*Heig[:values]))*Heig[:vectors]'
Uo=reshape(transpose(Uo),d,d,d,d)

#initilize state with both site and bond tensor
#|g>=|000...>
A=[]
for j=1:2 push!(A,zeros(1,1,d)); A[j][1]=1 end
B=[]
for j=1:2 push!(B,ones(1,1)) end

#itebd
t=0
while t<tf
    #The measurement assumes canonical form for MPS
    zzt_odd=jcontract([B[2],A[1],B[1],A[2],B[2],Z,Z,conj(B[2]),conj(A[1]),conj(B[1]),conj(A[2]),conj(B[2])],[[1,2],[2,3,7],[3,4],[4,5,8],[5,6],[7,9],[8,10],[1,11],[11,12,9],[12,13],[13,14,10],[14,6]])
    zzt_even=jcontract([B[1],A[2],B[2],A[1],B[1],Z,Z,conj(B[1]),conj(A[2]),conj(B[2]),conj(A[1]),conj(B[1])],[[1,2],[2,3,7],[3,4],[4,5,8],[5,6],[7,9],[8,10],[1,11],[11,12,9],[12,13],[13,14,10],[14,6]])
    zzt=(zzt_odd+zzt_even)/2

    @show t,zzt_odd,zzt_even
    @show t,zzt
    A,B=itebd_even_odd_one_step(A,B,Ue,Uo;eps=1e-8,chi=chi)
    t=t+δt
end

