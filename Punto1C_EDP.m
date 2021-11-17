x = 0:0.001:2*pi;
z = (exp(3*1i*x)-exp(2*1i*x))./(5/12-4/3*exp(1i*x)+23/12*exp(2*1i*x));
% (720 * (exp(4*1i*x)-exp(3*1i*x)))./(-19+106*exp(1i*x)-264*exp(2*1i*x) +646*exp(3*1i*x) +251*exp(4*1i*x));
i = imag(z);
r = real(z);
plot(r, i, '.')
fill(r, i, 'magenta')