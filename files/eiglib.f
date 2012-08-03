CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   package EIGLIB   -- from numerical recipes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine eig(A,n,np,e,vec,work)
C
C  takes in symmetric n x n matrix (phys dimension np x np)
C  computes the eigenvalues in e  and vectors in vec 
C   work is a dummy array 
C
      implicit none
      
      integer n,np
      real a(np,np)		! array to be diagonalized 
      real e(np)		! eigenvalues out 
      real vec(np,np)		! eigenvectors out
      real work(np)		!
      integer i,j
      
      do i =1,n
      	do j = 1,n
      	     vec(i,j)=a(i,j)
      	enddo
      enddo
      call tred2(vec,n,np,e,work)
      call tqli(e,work,n,np,vec)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE tred2(a,n,np,d,e)
      INTEGER n,np
      REAL a(np,np),d(np),e(np)
      INTEGER i,j,k,l
      REAL f,g,h,hh,scale
      do 18 i=n,2,-1
        l=i-1
        h=0.
        scale=0.
        if(l.gt.1)then
          do 11 k=1,l
            scale=scale+abs(a(i,k))
11        continue
          if(scale.eq.0.)then
            e(i)=a(i,l)
          else
            do 12 k=1,l
              a(i,k)=a(i,k)/scale
              h=h+a(i,k)**2
12          continue
            f=a(i,l)
            g=-sign(sqrt(h),f)
            e(i)=scale*g
            h=h-f*g
            a(i,l)=f-g
            f=0.
            do 15 j=1,l
C     Omit following line if finding only eigenvalues
              a(j,i)=a(i,j)/h
              g=0.
              do 13 k=1,j
                g=g+a(j,k)*a(i,k)
13            continue
              do 14 k=j+1,l
                g=g+a(k,j)*a(i,k)
14            continue
              e(j)=g/h
              f=f+e(j)*a(i,j)
15          continue
            hh=f/(h+h)
            do 17 j=1,l
              f=a(i,j)
              g=e(j)-hh*f
              e(j)=g
              do 16 k=1,j
                a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
16            continue
17          continue
          endif
        else
          e(i)=a(i,l)
        endif
        d(i)=h
18    continue
C     Omit following line if finding only eigenvalues.
      d(1)=0.
      e(1)=0.
      do 24 i=1,n
C     Delete lines from here ...
        l=i-1
        if(d(i).ne.0.)then
          do 22 j=1,l
            g=0.
            do 19 k=1,l
              g=g+a(i,k)*a(k,j)
19          continue
            do 21 k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
21          continue
22        continue
        endif
C     ... to here when finding only eigenvalues.
        d(i)=a(i,i)
C     Also delete lines from here ...
        a(i,i)=1.
        do 23 j=1,l
          a(i,j)=0.
          a(j,i)=0.
23      continue
C     ... to here when finding only eigenvalues.
24    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #<6V*1(.31..


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE tqli(d,e,n,np,z)
      INTEGER n,np
      REAL d(np),e(np),z(np,np)
CU    USES pythag
      INTEGER i,iter,k,l,m
      REAL b,c,dd,f,g,p,r,s,pythag
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.
      do 15 l=1,n
        iter=0
1       do 12 m=l,n-1
          dd=abs(d(m))+abs(d(m+1))
          if (abs(e(m))+dd.eq.dd) goto 2
12      continue
        m=n
2       if(m.ne.l)then
          if(iter.eq.30)pause 'too many iterations in tqli'
          iter=iter+1
          g=(d(l+1)-d(l))/(2.*e(l))
          r=pythag(g,1.)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.
          c=1.
          p=0.
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythag(f,g)
            e(i+1)=r
            if(r.eq.0.)then
              d(i+1)=d(i+1)-p
              e(m)=0.
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
C     Omit lines from here ...
            do 13 k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
13          continue
C     ... to here when finding only eigenvalues.
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.
          goto 1
        endif
15    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #<6V*1(.31..
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION pythag(a,b)
      REAL a,b,pythag
      REAL absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
          pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #<6V*1(.31..
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE eigsrt(d,v,n,np)
      INTEGER n,np
      REAL d(np),v(np,np)
      INTEGER i,j,k
      REAL p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #<6V*1(.31..
