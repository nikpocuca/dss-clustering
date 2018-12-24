
#Pkg.add("Distributions")
#Pkg.add("JLD")

function matnormpdf(X,M,Sigmastar,Psistar,detSig,detPsi,n,p)
    logdens=-n*p*log(2*pi)/2+(p/2)*(detSig)+(n/2)*(detPsi)-(1/2)*tr(Sigmastar*(X-M)*Psistar*(X-M)')

    return logdens
end

function zupdate_init(X,pig,M,LoadA,LoadB,detA,detB,n,p,N,G)
    #print(LoadA)
    zmat=zeros(N,G)
    zmat2=zeros(N,G)
    diff=zeros(N)
    #print(N)
    for g=1:G
        for i=1:N
            zmat[i,g]=log(pig[g])+matnormpdf(X[:,:,i],M[:,:,g],LoadA[:,:,g],LoadB[:,:,g],detA[g],detB[g],n,p)
            if isnan(zmat[i,g])
                return [1]
            end
        end
    end
    #

    for i=1:N
        if maximum(zmat[i,:])< -700
            difer=-700-maximum(zmat[i,:])
            zmat2[i,:]=zmat[i,:].+difer
            diff[i]=difer
        end
    end


    zmat2=exp.(zmat2)
    return [zmat2,diff,zmat]
end

function zup_init(zdens,N,G)
    z=zeros(N,G)

    for i in 1:N
            z[i,:]=zdens[i,:]/sum(zdens[i,:])
    end
    return z
end


function pigupdate_init(z,G,N)
    pig=zeros(G)
    for g=1:G
        pig[g]=sum(z[:,g])/N
    end
    return pig
end

function Mupdate_init(X,z,G,n,p,N)
    M=zeros(n,p,G)
    for g=1:G
        denom=sum(z[:,g])
        numer=zeros(n,p)
        for i=1:N
            numer=numer+z[i,g]*X[:,:,i]
        end
        M[:,:,g]=numer/denom
    end
    return M
end

function WAG_update_init(A,q,Sigmainv,G,Tol_Sing)
    prob=0
    iden=Matrix(I,q,q)
    WAG=zeros(q,q,G)
    for g=1:G
        WAGtemp=iden+A[:,:,g]'*Sigmainv[:,:,g]*A[:,:,g]
        rcond=1/cond(WAGtemp,1)
        if rcond<Tol_Sing
            return [1]
        end
        WAG[:,:,g]=inv(WAGtemp)
    end
    return [0,WAG]
end

function WBG_update_init(B,r,Psiinv,G,Tol_Sing)
    iden=Matrix(I,r,r)
    WBG=zeros(r,r,G)
    for g=1:G
        WBGtemp=iden+B[:,:,g]'*Psiinv[:,:,g]*B[:,:,g]
        rcond=1/cond(WBGtemp,1)
        if rcond<Tol_Sing
            return [1]
        end
        WBG[:,:,g]=inv(WBGtemp)
    end
    return [0,WBG]
end

function Estep2_init(M,Sigmainv,A,X,LoadB,G,p,N,q)
    WAG=WAG_update_init(A,q,Sigmainv,G,1e-9)
    if WAG[1]==1
        return [1]
    end
    WAG=WAG[2]
    aB=zeros(q,p,N,G)
    bB=zeros(q,q,N,G)

    for g=1:G
        WAGtemp=WAG[:,:,g]
        Atemp=A[:,:,g]
        Sigmatemp=Sigmainv[:,:,g]
        Mtemp=M[:,:,g]
        for i=1:N
            aB[:,:,i,g]=WAGtemp*Atemp'*Sigmatemp*(X[:,:,i]-Mtemp)
            bB[:,:,i,g]=p*WAGtemp+aB[:,:,i,g]*LoadB[:,:,g]*aB[:,:,i,g]'
        end
    end
    return [0,aB,bB]
end

function Estep3_init(M,Psiinv,B,X,LoadA,G,n,N,r)
    WBG=WBG_update_init(B,r,Psiinv,G,1e-9)
    if WBG[1]==1
        return [1]
    end
    WBG=WBG[2]
    aA=zeros(n,r,N,G)
    bA=zeros(r,r,N,G)

    for g=1:G
        WBGtemp=WBG[:,:,g]
        Btemp=B[:,:,g]
        Psitemp=Psiinv[:,:,g]
        Mtemp=M[:,:,g]
        for i=1:N
            aA[:,:,i,g]=(X[:,:,i]-Mtemp)*Psitemp*Btemp*WBGtemp
            bA[:,:,i,g]=n*WBGtemp+aA[:,:,i,g]'*LoadA[:,:,g]*aA[:,:,i,g]
        end
    end
    return [0,aA,bA]
end

function A_update_init(z,X,M,LoadB,aB,bB,Tol_Sing,N,G,n,q)
    Anew=zeros(n,q,G)
    for g=1:G
        Term1=zeros(n,q)
        Term2=zeros(q,q)
        Mtemp=M[:,:,g]
        aBtemp=aB[:,:,:,g]
        bBtemp=bB[:,:,:,g]
        LoadBtemp=LoadB[:,:,g]
        for i=1:N
            Term1=Term1+z[i,g]*(X[:,:,i]-Mtemp)*LoadBtemp*aBtemp[:,:,i]'
            Term2=Term2+z[i,g]*bBtemp[:,:,i]
        end
        rcond=cond(Term2,1)
        if rcond<Tol_Sing
            return [1]
        end
        Term2inv=inv(Term2)
        Anew[:,:,g]=Term1*Term2inv
    end
    return [0,Anew]
end


function Sigmaup_init(z,X,M,A,LoadB,aB,bB,n,p,N,G,Tol_Sing)
    Sigmat=zeros(n,n,G)
    Sigmainv=zeros(n,n,G)
    Sigmadet=zeros(G)
    for g=1:G
        denom=sum(z[:,g])*p
        Term1=zeros(n,n)
        Term2=zeros(n,n)
        for i=1:N
            Term1=Term1+z[i,g]*(X[:,:,i]-M[:,:,g])*LoadB[:,:,g]*(X[:,:,i]-M[:,:,g])'
            Term2=Term2+z[i,g]*A[:,:,g]*aB[:,:,i,g]*LoadB[:,:,g]*(X[:,:,i]-M[:,:,g])'
        end
        Sigmatemp=(Term1-Term2)/denom
        Sigma=Diagonal(Sigmatemp)
        if minimum(diag(Sigma))<Tol_Sing
            #print(z)
            #print(diag(Sigma))
            return[1]
        end
        Sigmat[:,:,g]=Sigma
        Sigmainv[:,:,g]=inv(Sigma)
        Sigmadet[g]=logdet(Sigmainv[:,:,g])
    end
    return [0,Sigmat,Sigmainv,Sigmadet]
end

function B_update_init(z,X,M,LoadA,aB,bB,Tol_Sing,N,G,p,r)
    Bnew=zeros(p,r,G)
    for g=1:G
        Term1=zeros(p,r)
        Term2=zeros(r,r)
        Mtemp=M[:,:,g]
        aBtemp=aB[:,:,:,g]
        bBtemp=bB[:,:,:,g]
        LoadAtemp=LoadA[:,:,g]
        for i=1:N
            Term1=Term1+z[i,g]*(X[:,:,i]-Mtemp)'*LoadAtemp*aBtemp[:,:,i]
            Term2=Term2+z[i,g]*bBtemp[:,:,i]
        end
        rcond=cond(Term2,1)
        if rcond<Tol_Sing
            return [1]
        end
        Term2=inv(Term2)
        Bnew[:,:,g]=Term1*Term2
    end
    return [0,Bnew]
end


function Psiup_init(z,X,M,B,LoadA,aB,bB,n,p,N,G,Tol_Sing)
    Psiinv=zeros(p,p,G)
    Psit=zeros(p,p,G)
    Psidet=zeros(G)
    for g=1:G
        denom=sum(z[:,g])*n
        Term1=zeros(p,p)
        Term2=zeros(p,p)
        for i=1:N
            Term1=Term1+z[i,g]*(X[:,:,i]-M[:,:,g])'*LoadA[:,:,g]*(X[:,:,i]-M[:,:,g])
            Term2=Term2+z[i,g]*B[:,:,g]*aB[:,:,i,g]'*LoadA[:,:,g]*(X[:,:,i]-M[:,:,g])
        end
        Psitemp=Term1-Term2
        Psi=Diagonal(Psitemp)/denom
        #print(Psi)

        if minimum(diag(Psi))<Tol_Sing
            return[1]
        end
        Psit[:,:,g]=Psi
        Psiinv[:,:,g]=inv(Psi)
        Psidet[g]=logdet(Psiinv[:,:,g])
    end
    return [0,Psit,Psiinv,Psidet]
end

function AECM1_init(z,X,G,n,p,N)
    pig=pigupdate_init(z,G,N)
    M=Mupdate_init(X,z,G,n,p,N)
    return [pig,M]
end

function AECM2_init(z,X,M,A,Sigmainv,LoadB,N,n,p,q,G,Tol_Sing)
    Estep=Estep2_init(M,Sigmainv,A,X,LoadB,G,p,N,q)

    if Estep[1]==1
        #print("Probelm with sedond E Step")
        return [1]
    end

    aB=Estep[2]
    bB=Estep[3]
    Anew=A_update_init(z,X,M,LoadB,aB,bB,Tol_Sing,N,G,n,q)
    if Anew[1]==1
        print("Problem with A update")
        return [1]
    end
    Anew=Anew[2]
    Sigmanew=Sigmaup_init(z,X,M,A,LoadB,aB,bB,n,p,N,G,Tol_Sing)
    if Sigmanew[1]==1
        print("Problem with Sigma Update")
        return [1]
    end
    LoadA=zeros(n,n,G)
    detA=zeros(G)
    for g=1:G
        LoadAtemp=Sigmanew[2][:,:,g]+Anew[:,:,g]*Anew[:,:,g]'
        rcond=cond(LoadAtemp)
        if rcond<Tol_Sing
            print("Singular LoadA")
            return [1]
        end
        LoadA[:,:,g]=inv(LoadAtemp)
        detA[g]=logdet(LoadA[:,:,g])
    end
    return [0,Anew,Sigmanew[3],detA,LoadA]
end

function AECM3_init(z,X,M,B,Psiinv,LoadA,N,n,p,r,G,Tol_Sing)
    Estep=Estep3_init(M,Psiinv,B,X,LoadA,G,n,N,r)

    if Estep[1]==1
        #print("Probelm with third E Step")
        return [1]
    end

    aB=Estep[2]
    bB=Estep[3]
    Bnew=B_update_init(z,X,M,LoadA,aB,bB,Tol_Sing,N,G,p,r)
    if Bnew[1]==1
        print("Problem with B update")
        return [1]
    end
    Bnew=Bnew[2]
    Psinew=Psiup_init(z,X,M,B,LoadA,aB,bB,n,p,N,G,Tol_Sing)
    if Psinew[1]==1
        print("Problem with Psi Update")
        return [1]
    end
    LoadB=zeros(p,p,G)
    detB=zeros(G)
    for g=1:G
        LoadBtemp=Psinew[2][:,:,g]+Bnew[:,:,g]*Bnew[:,:,g]'
        rcond=cond(LoadBtemp)
        if rcond<Tol_Sing
            print("Singular LoadB")
            return [1]
        end
        LoadB[:,:,g]=inv(LoadBtemp)
        detB[g]=logdet(LoadB[:,:,g])
    end
    return [0,Bnew,Psinew[3],detB,LoadB]
end

function likcalc(zdens,N,class,z)
    likeclass=0
    for i=1:N
        like=0
        if class[i]==0
            like=log(sum(zdens[1][i,:]))-zdens[2][i]
        else
            like=sum(z[i,:].*zdens[3][i,:])
        end
        likeclass=likeclass+like
    end
    return likeclass
end


# X is an array of matrices n x n x N 
# q is the column factors 
# r is the row factors 
# G number of groups 
# sets random seed. 
# max iterations is the EM steps for iterations 
# class is 0 if unknown, and group number if class is known , equal to the number of observations 
# BIC is in the thirteenth position of the list 

function EM_Main_init(X,q,r,G,seedno,maxiter,class)
    #print(class)
    Random.seed!(seedno)
    #d=Uniform(-1,1)
    #dA=Normal(0,1)
    class=convert(Array{Int,1},class)
    N=size(X)[3]
    n=size(X)[1]
    p=size(X)[2]
    zinit=zeros(N,G)
    for i=1:N
        temp=rand(1:100,G)
        zinit[i,:]=temp/sum(temp)
    end
    for i=1:N
        if class[i]>0
            temp=zeros(G)
            temp[class[i]]=1
            zinit[i,:]=temp
        end
    end
    #print(zinit)
    #zinit=mymap(zinit)
    piginit=pigupdate_init(zinit,G,N)
    #print(piginit)
    Minit=zeros(n,p,G)
    Sigmastarinit=zeros(n,n,G)
    Sigmainit=zeros(n,n,G)
    Psistarinit=zeros(p,p,G)
    Psiinit=zeros(p,p,G)
    Loada=zeros(n,q,G)
    Loadb=zeros(p,r,G)
    detSig=zeros(G)
    detPsi=zeros(G)
    for g=1:G
        for i=1:N
            Minit[:,:,g]=Minit[:,:,g]+zinit[i,g]*X[:,:,i]/sum(zinit[:,g])
        end
        Sigmatemp=zeros(n,n)
        Psitemp=zeros(p,p)
        for i=1:N
            Sigmatemp=Sigmatemp+(zinit[i,g]*(X[:,:,i]-Minit[:,:,g])*(X[:,:,i]-Minit[:,:,g])'/(sum(zinit[:,g])*p))
            Psitemp=Psitemp+(zinit[i,g]*(X[:,:,i]-Minit[:,:,g])'*(X[:,:,i]-Minit[:,:,g])/(n*sum(zinit[:,g])))
        end
        #print(1/cond(Sigmatemp))
        Sigmainit[:,:,g]=inv(Diagonal(Sigmatemp))
        Psiinit[:,:,g]=inv(Diagonal(Psitemp))
        Loada[:,:,g]=reshape(rand(-1:0.1:1,n*q),n,q)
        Loadb[:,:,g]=reshape(rand(-1:0.1:1,p*r),p,r)
        Sigmastarinit[:,:,g]=inv(inv(Sigmainit[:,:,g])+Loada[:,:,g]*Loada[:,:,g]')
        #print(1/cond(Sigmastarinit[:,:,g]))
        Psistarinit[:,:,g]=inv(inv(Psiinit[:,:,g])+Loadb[:,:,g]*Loadb[:,:,g]')
        detSig[g]=logdet(Sigmastarinit[:,:,g])
        detPsi[g]=logdet(Psistarinit[:,:,g])
    end
    pig=piginit
    M=Minit
    #print(Minit)
    Sigmastar=Sigmastarinit
    Sigma=Sigmainit
    Psistar=Psistarinit
    Psi=Psiinit
    #print(detSig)


    #zdens=zupdate_init(X,pig,M,Sigmastar,Psistar,detSig,detPsi,n,p,N,G)
    #z=zup_init(zdens[1],N,G)
    #print(z)
    #print("This is Thursday")
    #print(z)
    z=zinit
    #print(z)
    conv=0
    tol=0
    lik=zeros(maxiter)
    iter=1
    for i=1:maxiter
        ECM1=AECM1_init(z,X,G,n,p,N)
        pig=ECM1[1]
        #print(pig[1])
        if pig[1]==NaN
            return [1]
        end
        #print(pig)
        M=ECM1[2]

        ECM2=AECM2_init(z,X,M,Loada,Sigma,Psistar,N,n,p,q,G,1e-32)
        if ECM2[1]==1
            return [1]
        end

        Loada=ECM2[2]
        Sigma=ECM2[3]
        detSig=ECM2[4]
        Sigmastar=ECM2[5]
        #print(detSig)

        ECM3=AECM3_init(z,X,M,Loadb,Psi,Sigmastar,N,n,p,r,G,1e-32)

        if ECM3[1]==1
            return [1]
        end

        Loadb=ECM3[2]
        Psi=ECM3[3]
        detPsi=ECM3[4]
        Psistar=ECM3[5]

        zdens=zupdate_init(X,pig,M,Sigmastar,Psistar,detSig,detPsi,n,p,N,G)
        if zdens[1]==1
            return [1]

        end
        z=zup_init(zdens[1],N,G)
        for i=1:N
            if class[i]>0
                temp=zeros(G)
                temp[class[i]]=1
                z[i,:]=temp
            end
        end
        classpred=mymap2(z)
        #print(z)
        #print(z)
        lik[i]=likcalc(zdens,N,class,z)
        if iter==5
            tol=sign(lik[i])*lik[i]/10^4
        end


        if iter>1
            if (lik[iter]-lik[iter-1])<0
                print("Decreasing Likelihood Q=$q R=$r G=$G")
                print(lik[iter]-lik[iter-1])
                return [1,pig,z,iter,lik[1:maxiter],M,Sigmastar,Psistar]
            end
        end
        if (iter>5)
            if lik[iter-1]-lik[iter-2]==0
                BIC=BICcalc(lik[1:iter],n,p,q,r,G,N)
                return [0,pig,M,Sigma,Psi,Loada,Loadb,Sigmastar,Psistar,z,classpred,lik[1:iter],BIC]
            end

            ak=(lik[iter]-lik[iter-1])/(lik[iter-1]-lik[iter-2])
            linf=lik[iter-1]+(lik[iter]-lik[iter-1])/(1-ak)

            if (abs(linf-lik[iter-1]))<tol
                print("converged")
                print([q,r])
                BIC=BICcalc(lik[1:iter],n,p,q,r,G,N)
                return [0,pig,M,Sigma,Psi,Loada,Loadb,Sigmastar,Psistar,z,classpred,lik[1:iter],BIC,class,[q,r]]
            end
        end

        iter=iter+1
        #print(iter)
    end
    #print(z)
    print("Maximum Iterations Reached")
    return [1,pig,z,iter,lik[1:maxiter],M,Sigmastar,Psistar,Sigma,Psi,detSig,detPsi,Loada,Loadb,-Inf]
end


function Simulx(M,Sigma,Psi,n,p,G,NG)
    x=matnormal(NG[1],n,p,M[:,:,1],Sigma[:,:,1],Psi[:,:,1])
    for g=2:G
        x=cat(3,x,matnormal(NG[g],n,p,M[:,:,g],Sigma[:,:,g],Psi[:,:,g]))
    end
    return x
end

function BICcalc(ll,n,p,q,r,G,N)
    free=(G-1)+G*(n*p+n*q+n-0.5*q*(q-1)+p*r+p-0.5*r*(r-1))
    BIC=2*maximum(ll)-free*log(N)
    return BIC
end

function mymap(z,N)
    class=zeros(N)
    for i=1:N
        class[i]=indmax(z[i,:])
    end
    return class
end
function mymap2(z)
    N=size(z)[1]
    class=zeros(N)
    for i=1:N
        class[i]=argmax(z[i,:])
    end
    return class
end
