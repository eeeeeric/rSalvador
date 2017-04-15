# July 31, 2016: adding Fisher information to rsalvador 1.6

#---- June 18, 2016: computing expected Fisher info for m of the MK distribution ---

fisher.MK=function(m,w,n) {

p=prob.MK(m=m,w,n=n)

r=1/w

B=( 2:n -1)/(2:n +r)

B=Reduce(prod, B, init=1,accumulate=TRUE)/(1+r)

eta=1/(1:n +r)

eta=Reduce(sum,eta[-1],init=eta[1],accumulate=TRUE)

h10=append(r*B,-1,after=0)

q.1=seq.convolute(p,h10)

U=0 

for (i in 1:(n+1) ) {

U=U+(q.1[i])^2/p[i]

}

return( U )

}



# --------- June 19, 2016: fisher info with plating -------------

fisher.LD.plating=function(m,e,n){

k=e/(1-e)

eta.seq=etaSeq(e,n)

p0=prob.LD.plating(m,e,n)

p1=seq.convolute(eta.seq,p0)*k

U=0

for (i in 1:(n+1) ) {

 U=U+(p1[i])^2/p0[i] }

return(U)

}



### April 12, 2017: adding sample size calculation to rsalvador 1.7

### MK model ###

samp.size.MK=function(m=1.2,w=0.9,psi=0.25,trunc=3000) {

fisher=fisher.MK(m,w,trunc)

n=(1.96/m/psi)^2 / fisher

return( ceiling(n) )

}

### LD model with partial plating ###

samp.size.LD.plating=function(m=10.2,e=0.2,psi=0.25,trunc=3000) {

fisher=fisher.LD.plating(m,e,trunc)

n=(1.96/m/psi)^2 / fisher

return( ceiling(n) )

}


## P0 method with partial plating, April 12, 2017


p0.LD.plating=function (data, e) {
    p0 = min(data)
    if (p0 == 0) {
        return( p0.plating(data=data, e=e) )
    }
    else {
        message("The P0 method is not applicable to the given data.")
    }
}
