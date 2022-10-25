import numpy as np
A=list()
B=list()
C=list()
D=list()

countA=0
countB=0
countC=0
countD=0
#Heron's Formula to find Areas
def heron(a,b,c):
    s=(a+b+c)/2
    area=np.sqrt(s*(s-a)*(s-b)*(s-c))
    return area

def distance(x1,y1,z1,x2,y2,z2):
    a=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2
    d=np.sqrt(a)
    return d

def areaQuad(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4):
    a=distance(x1,y1,z1,x2,y2,z2)
    b=distance(x2,y2,z2,x3,y3,z3)
    c=distance(x3,y3,z3,x1,y1,z1)
    Atriangle1=heron(a,b,c)
    p=distance(x1,y1,z1,x3,y3,z3)
    q=distance(x3,y3,z3,x4,y4,z4)
    r=distance(x4,y4,z4,x1,y1,z1)
    Atriangle2=heron(p,q,r)
    Area=Atriangle1+Atriangle2
    return Area
#Entering Coordinates of the Quadrilateral
print('For 2D diagrams in XY plane use Z=0')
while countA<3:
    Acoord=input('Enter X,Y,Z coordinates of point A:')
    A.append(float(Acoord))
    countA=countA+1
print('A:',A)
while countB<3:
    Bcoord=input('Enter X,Y,Z coordinates of point B:')
    B.append(float(Bcoord))
    countB=countB+1
print('B:',B)
while countC<3:
    Ccoord=input('Enter X,Y,Z coordinates of point C:')
    C.append(float(Ccoord))
    countC=countC+1
print('C:',C)
while countD<3:
    Dcoord=input('Enter X,Y,Z coordinates of point D:')
    D.append(float(Dcoord))
    countD=countD+1
print('D:',D)
print('\nQuadrilateral ABCD with Coordinates')
print('A:',A)
print('B:',B)
print('C:',C)
print('D:',D)

#Checking if the entered coordinates are coplanar (AB X AD) . AC = 0
CheckPlanei=B[1]*D[2]-B[1]*A[2]-A[1]*D[2]-B[2]*D[1]+B[2]*A[1]+A[2]*D[1]
CheckPlanej=B[2]*D[0]-B[2]*A[0]-A[2]*D[0]-B[0]*D[2]+B[0]*A[2]+A[0]*D[2]
CheckPlanek=B[0]*D[1]-B[0]*A[1]-A[0]*D[1]-B[1]*D[0]+B[1]*A[0]+A[1]*D[0]
AC=[C[0]-A[0],C[1]-A[1],C[2]-A[2]]
CheckPlane=CheckPlanei*AC[0]+CheckPlanej*AC[1]+CheckPlanek*AC[2]

if CheckPlane==0:
    print('Points are coplanar')
else:
    print('\nPlease Enter Valid Coordinates ( Enter Coplanar coordinates)')
    quit()

#Area of Triangle ABC=0.5*(AB x AC) cross product method
AreaABCi=B[1]*C[2]-B[1]*A[2]-A[1]*C[2]-B[2]*C[1]+B[2]*A[1]+A[2]*C[1]
AreaABCj=B[2]*C[0]-B[2]*A[0]-A[2]*C[0]-B[0]*C[2]+B[0]*A[2]+A[0]*C[2]
AreaABCk=B[0]*C[1]-B[0]*A[1]-A[0]*C[1]-B[1]*C[0]+B[1]*A[0]+A[1]*C[0]
AreaABC=(np.sqrt(AreaABCi**2+AreaABCj**2+AreaABCk**2))*0.5

#Area of Triangle ACD=0.5*(AC x AD) cross product method
AreaACDi=C[1]*D[2]-C[1]*A[2]-A[1]*D[2]-C[2]*D[1]+C[2]*A[1]+A[2]*D[1]
AreaACDj=C[2]*D[0]-C[2]*A[0]-A[2]*D[0]-C[0]*D[2]+C[0]*A[2]+A[0]*D[2]
AreaACDk=C[0]*D[1]-C[0]*A[1]-A[0]*D[1]-C[1]*D[0]+C[1]*A[0]+A[1]*D[0]
AreaACD=(np.sqrt(AreaACDi**2+AreaACDj**2+AreaACDk**2))*0.5
AreaABCD=AreaABC+AreaACD
#print('Area ABCD (by cross product):',AreaABCD)
heronarea1=areaQuad(A[0],A[1],A[2],B[0],B[1],B[2],C[0],C[1],C[2],D[0],D[1],D[2])
#print("Area ABCD (by Heron's Formula ):",format(heronarea1,".3f"))

#mutually perpenducular unit vectors e1,e2,m
magAB=np.sqrt(((B[0]-A[0])**2)+((B[1]-A[1])**2)+((B[2]-A[2])**2))
e1=[(B[0]-A[0])/magAB,(B[1]-A[1])/magAB,(B[2]-A[2])/magAB]
print('\ne1:',e1)
magAC =np.sqrt(((C[0]-A[0])**2)+((C[1]-A[1])**2)+((C[2]-A[2])**2))
b1=[(C[0]-A[0])/magAC,(C[1]-A[1])/magAC,(C[2]-A[2])/magAC]
print('b1:',b1)
# M=e1 x b1
M=[e1[1]*b1[2]-b1[1]*e1[2],b1[0]*e1[2]-e1[0]*b1[2],e1[0]*b1[1]-e1[1]*b1[0]]
MagM=np.sqrt(M[0]**2+M[1]**2+M[2]**2)
m=[M[0]/MagM,M[1]/MagM,M[2]/MagM]
print('m:',m)
# E2=m x e1
E2=[m[1]*e1[2]-e1[1]*m[2],e1[0]*m[2]-m[0]*e1[2],m[0]*e1[1]-e1[0]*m[1]]
MagE2=np.sqrt(E2[0]**2+E2[1]**2+E2[2]**2)
e2=[E2[0]/MagE2,E2[1]/MagE2,E2[2]/MagE2]
print('e2:',e2)

#center of gravity or centroid
cg=[(A[0]+B[0]+C[0]+D[0])/4,(A[1]+B[1]+C[1]+D[1])/4,(A[2]+B[2]+C[2]+D[2])/4]
print('\ncentroid cg:',cg)

Acg=[cg[0]-A[0],cg[1]-A[1],cg[2]-A[2]]
Bcg=[cg[0]-B[0],cg[1]-B[1],cg[2]-B[2]]
Ccg=[cg[0]-C[0],cg[1]-C[1],cg[2]-C[2]]
Dcg=[cg[0]-D[0],cg[1]-D[1],cg[2]-D[2]]

xx1=(Acg[0]*e1[0])+(Acg[1]*e1[1])+(Acg[2]*e1[2]) #Acg . e1
yy1=(Acg[0]*e2[0])+(Acg[1]*e2[1])+(Acg[2]*e2[2]) #Acg . e2
xx2=(Bcg[0]*e1[0])+(Bcg[1]*e1[1])+(Bcg[2]*e1[2]) #Bcg . e1
yy2=(Bcg[0]*e2[0])+(Bcg[1]*e2[1])+(Bcg[2]*e2[2]) #Bcg . e2
xx3=(Ccg[0]*e1[0])+(Ccg[1]*e1[1])+(Ccg[2]*e1[2]) #Ccg . e1
yy3=(Ccg[0]*e2[0])+(Ccg[1]*e2[1])+(Ccg[2]*e2[2]) #Ccg . e2
xx4=(Dcg[0]*e1[0])+(Dcg[1]*e1[1])+(Dcg[2]*e1[2]) #Dcg . e1
yy4=(Dcg[0]*e2[0])+(Dcg[1]*e2[1])+(Dcg[2]*e2[2]) #Dcg . e2

P=[xx1,yy1]
Q=[xx2,yy2]
R=[xx3,yy3]
S=[xx4,yy4]
print('\nArea ABCD (by cross product):',AreaABCD)
print("Area ABCD (by Heron's Formula ):",format(heronarea1,".3f"))
# Area of PQRS
#AreaPQRS=0.5*(P[0]*Q[1]+Q[0]*R[1]+R[0]*S[1]+S[0]*P[1]-P[1]*Q[0]-Q[1]*R[0]-R[1]*S[0]-S[1]*P[0])
#AreaPQRS=0.5*((Q[0]-P[0])*(Q[1]+P[1])+(R[0]-Q[0])*(R[1]+Q[1])+(R[0]-S[0])*(R[1]+S[1])+(S[0]-P[0])*(S[1]+P[1]))
Areapqrs=0.5*((xx2-xx1)*(yy2+yy1)+(xx3-xx2)*(yy3+yy2)+(xx4-xx3)*(yy4+yy3)+(xx1-xx4)*(yy1+yy4))
AreaPQRS=abs(Areapqrs)
#AreaPQRS=0.5*((xx1*yy2+xx2*yy3+xx3*yy4+xx4*yy1)-(xx2*yy1+xx3*yy2+xx4*yy3+xx1*yy4))
print('\nArea PQRS :',round(AreaPQRS,4))
heronarea2=areaQuad(P[0],P[1],0,Q[0],Q[1],0,R[0],R[1],0,S[0],S[1],0)
print("Area PQRS (by Heron's Formula ):",format(heronarea2,".3f"))

 # some inputs try (−4,−2,0),(−3,−5,0),(3,−2,0),(2,3,0). giving very small error in order of  10^-15 so rounded off to 4 digits
if (abs(AreaABCD-AreaPQRS)<0.001):
    print('\nconfirmed that Areas are same')
#else:
    #if heronarea1==heronarea2:
        #print('\nconfirmed that Areas are same (areas by herons formula)')
else :
     print("\nAreas don't match")
