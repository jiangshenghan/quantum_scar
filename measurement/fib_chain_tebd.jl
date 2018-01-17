#measure observables using tebd method. Fib chain Hilbert space are imposed by soft constraint term
# H=\sum_j X_j + V \sum_j 1/4(I_j-Z_j).(I_{j+1}-Z_{j+1}), with infinite V limit

include("/home/jiangsh/JTensor.jl/src/JTensor.jl")
using JTensor

d=2
L=20
tf=100
δt=0.05
@show L,tf,δt

I=[1 0; 0 1]
X=[0 1; 1 0]
Y=[0 -im; im 0]
Z=[1 0; 0 -1]

#construct Hamiltonian
V=100
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
while t<tf
    #The measurement assumes canonical form for MPS
    #zzt=jcontract([A[1],B[1],A[2],B[2],Z,Z,conj(A[1]),conj(B[1]),conj(A[2]),conj(B[2])],[[1,2,6],[2,3],[3,4,7],[4,5],[6,8],[7,9],[1,10,8],[10,11],[11,12,9],[12,5]])
    cind=div(L,2)
    zzt=jcontract([B[cind-1],A[cind],B[cind],A[cind+1],B[cind+1],Z,Z,conj(B[cind-1]),conj(A[cind]),conj(B[cind]),conj(A[cind+1]),conj(B[cind+1])],[[1,2],[2,3,7],[3,4],[4,5,8],[5,6],[7,9],[8,10],[1,11],[11,12,9],[12,13],[13,14,10],[14,6]])

    @show t,zzt
    A,B=tebd_even_odd_one_step(A,B,Ue,Uo;eps=1e-6)
    t=t+δt
end
