##CRANK-NICOLSON
def CkNicolson(func1,func2,funcSol1,funcSol2,a,b,n,y1,y2):
  h = abs(b-a)/n
  X, Y1a, Y2a = BkEuler(func1,func2,funcSol1,funcSol2,a,b,n,y1,y2)
  X, Y1b, Y2b = FrEuler(func1,func2,funcSol1,funcSol2,a,b,n,y1,y2)
  Y = np.zeros(n+1)
  prom1 = (Y1a+Y1b)/2
  prom2 = (Y2a+Y2b)/2
  return X, prom1, prom2
CkNicolson(func1,func2,funcSol1,funcSol2,0,1,2,4/3,2/3)

# 100 linearly spaced numbers
x = np.linspace(0,1,1000)
real1 = np.zeros(len(x))
real2 = np.zeros(len(x))
for i in range(len(x)):
  real1[i] = funcSol1(x[i])
  real2[i] = funcSol2(x[i])

x1,y11,y21=CkNicolson(func1,func2,funcSol1,funcSol2,0,1,2,(4/3),(2/3))
x2,y12,y22=CkNicolson(func1,func2,funcSol1,funcSol2,0,1,4,(4/3),(2/3))
x3,y13,y23=CkNicolson(func1,func2,funcSol1,funcSol2,0,1,8,(4/3),(2/3))
x4,y14,y24=CkNicolson(func1,func2,funcSol1,funcSol2,0,1,16,(4/3),(2/3))
x5,y15,y25=CkNicolson(func1,func2,funcSol1,funcSol2,0,1,32,(4/3),(2/3))
x6,y16,y26=CkNicolson(func1,func2,funcSol1,funcSol2,0,1,64,(4/3),(2/3))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.plot(x,real2, 'b', label='Exacta')
plt.plot(x1,y11, 'g', label='h=1/2')
plt.plot(x2,y12, 'r', label='h=1/4')
plt.plot(x3,y13, 'c', label='h=1/8')
plt.plot(x4,y14, 'm', label='h=1/16')
plt.plot(x5,y15, 'y', label='h=1/32')
plt.plot(x6,y16, 'k', label='h=1/64')
plt.ylim(top = 1, bottom = -2)
plt.legend(loc='upper rigth')
# show the plot
