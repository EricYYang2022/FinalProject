clear; close all; clc;

dimens = 20;
a = linspace(0, pi/2, dimens);
s = linspace(0.1, 6, dimens);

syms r;
syms t;

%Set up cross products (n) of s_u cross s_v for the rocket head based on
%values of s. The cross products are predetermined through calculations.
cross_s1 = cell(dimens, 1);
cross_s2 = [-2*cos(t), -2*sin(t), 0];
cross_s3 = [0, 0 , r];
for ii = 1:dimens
    cross_s1{ii} = ...
        [-2*r^2*s(ii)*cos(t), ...
         -2*r^2*s(ii)*sin(t), ...
         r];
end

%Calculating cos(phi) using definition of dot product. Since the rocket is
%symmetrical about its center line, we just need an angle from this center
%line to determine our angle of direction, as all directions with this same
%angle alpha should generate the same angle phi.
unit_v = cell(dimens, 1);
for jj = 1:dimens
    unit_v{jj} = [0, sin(a(jj)), cos(a(jj))];
end

%Now we loop through for all values of s and a.
cosphi = cell(1,3);
M = zeros(dimens, dimens);
N = zeros(dimens, dimens);
for kk = 1:dimens
    %angle first
    u = unit_v{kk};
        cosphi{1,2} = dot(u, cross_s2)/norm(cross_s2, 2);
        cosphi{1,3} = dot(u, cross_s3)/norm(cross_s3, 2);
       
        %sometimes this stuff generates 0's so we can't use matlabfunction
        integraled2 = 1/pi * (cosphi{1, 2})^3 * norm(cross_s2, 2);
        int_once2 = int(integraled2, r, 0, 10);
        int_twice2 = int(int_once2, t, 0, 2*pi);
        
        integraled3 = 1/pi * (cosphi{1, 3})^3 * norm(cross_s3, 2);
        int_once3 = int(integraled3, r, 0, 2);
        int_twice3 = int(int_once3, t, 0, 2*pi);
        
        sumflux = int_twice2 + int_twice3;
    for ll = 1:dimens

        n1 = cross_s1{ll};
        cosphi{1,1} = dot(u, n1)/norm(n1, 2); %will be scaled to be F
        normm = norm(cross_s1{ll}, 2); %||ru x rv||
        cosphii = matlabFunction(1./pi.*(cosphi{1,1}).^3.*normm);
        int_twice = integral2(cosphii, 0, 2, 0, 2*pi);
%         integraled = 1/pi * (cosphi{1,1})^3;
%         int_once = int(integraled, r, 0, 5);
%         int_twice = int(int_once, t, 0, 2*pi);

%         M(ll, kk) = double(sumflux + int_twice);
        M(ll, kk) = sumflux + int_twice;
        volumee = 10*pi*2^2 + ...
        integral2(@(rr,tt)(-s(ll).*(rr.^2) + 4.*s(ll)), 0, 2, 0, 2*pi)
        N(ll, kk) = M(ll, kk)/volumee;
%         disp(M(ll, kk));
%         disp(N(ll, kk));
            
    end
end
hold on; 
surf(s, a*180/pi, M); 
view([60, 60, 60])
title('Friction Coefficient')
xlabel('s')
ylabel(['\alpha (' 176 ')'])

figure;
surf(s, a*180/pi, N);
view([60, 60, 60])
title('Friction Coefficient per mass')
xlabel('s')
ylabel(['\alpha (' 176 ')'])