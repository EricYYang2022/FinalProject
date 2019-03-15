syms r;
syms t;

s = 0.7211; a = 0.0827;
flux = 1/pi*((-2*r^2*s*sin(t)*sin(a)+r*cos(a))^3)/(4*r^4*s^2+r);

fuckonce = int(flux, r, 0, 2);
fucktwice = int(fuckonce, t, 0, 2*pi);
disp(double(fucktwice));
fuckflux = matlabFunction(flux);
f = integral2(fuckflux, 0, 2, 0, 2*pi);
disp(f);