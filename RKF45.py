##Runge Kutta

def RKF(func1, func2, a, b, n, y1, y2):
  h = abs(b-a)/n
  
  x = [a]
  re1 = [y1]
  re2 = [y2]
  K = lambda a1, a2, c: (func1(x[i-1] + h*a1, re1[i-1] + h*a2, re2[i-1] + h*c),
                             func2(x[i-1] + h*a1,re1[i-1] + h*a2, re2[i-1] + h*c))
  
  auxReal1 = [y1]
  auxReal2 = [y2]
  
  s2 = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]
  s1 = [[0, 0, 0, 0, 0, 0],
      [1/4, 0, 0, 0, 0, 0],
      [3/32, 9/32, 0, 0, 0, 0],
      [1932/2197, -7200/2197, 7296/2197, 0, 0, 0],
      [439/216, -8, 3680/513, -845/4104, 0, 0],
      [-8/27, 2, -3544/2565, 1859/4104, -11/40, 0]]
  s3 = [0, 1/4, 3/8, 12/13, 1, 1/2]
  i = 1
  
  for i in range(1, n+1):
    k1Y1, k1Y2 = K(0, 0, 0)
    k2Y1, k2Y2 = K(s1[1][0]*k1Y1,
                    s1[1][0]*k1Y2, s3[1])
    k3Y1, k3Y2 = K(s1[2][0]*k1Y1 + s1[2][1]*k2Y1,
                    s1[2][0]*k1Y2 + s1[2][1]*k2Y2, s3[2])
    k4Y1, k4Y2 = K(s1[3][0]*k1Y1 + s1[3][1]*k2Y1 + s1[3][2]*k3Y1,
                    s1[3][0]*k1Y2 + s1[3][1]*k2Y2 + s1[3][2]*k3Y2, s3[3])
    k5Y1, k5Y2 = K(s1[4][0]*k1Y1 + s1[4][1]*k2Y1 + s1[4][2]*k3Y1 + s1[4][3]*k4Y1,
                    s1[4][0]*k1Y2 + s1[4][1]*k2Y2 + s1[4][2]*k3Y2 + s1[4][3]*k4Y2, s3[4])
    k6Y1, k6Y2 = K(s1[5][0]*k1Y1 + s1[5][1]*k2Y1 + s1[5][2]*k3Y1 + s1[5][3]*k4Y1 + s1[5][4]*k5Y1,
                    s1[5][0]*k1Y2 + s1[5][1]*k2Y2 + s1[5][2]*k3Y2 + s1[5][3]*k4Y2 + s1[5][4]*k5Y2, s3[5])

    re1 = np.append(re1, re1[i-1] + h*(s2[0]*k1Y1 + s2[1]*k2Y1 + s2[2]*k3Y1 + s2[3]*k4Y1 + s2[4]*k5Y1 + s2[5]*k6Y1))
    re2 = np.append(re2, re2[i-1] + h*(s2[0]*k1Y2 + s2[1]*k2Y2 + s2[2]*k3Y2 + s2[3]*k4Y2 + s2[4]*k5Y2 + s2[5]*k6Y2))
    

    x = np.append(x, x[i-1] + h)
      

    auxReal1 = np.append(auxReal1, re1[i])
    auxReal2 = np.append(auxReal2, re2[i])

    

  return x, auxReal1, auxReal2
    
RKF(func1, func2, 0,1,2,4/3, 2/3)
