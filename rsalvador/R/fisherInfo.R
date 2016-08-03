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

