N = 5;
d   = 0.5;
deg = 0;

ph0 = -pi/2:0.01:pi/2;           % sudut azimuth (degree)
deg = deg *(pi/180);       % perbedaan phase tiap arus (rad)
a   = -2*pi*d*cos(deg);
ps0 = 2*pi*d*cos(ph0) ;  % phase dalam domain psi
AF   = sin(N.*ps0/2) ./ (N.*sin(ps0/2)); % array factor UE, ESLA

g   = 1;                    % uniform 
f   = abs(AF.*g);
f   = f / max(f);

figure(1);
ph1 = ph0*180/pi;
plot(ph1,f);
xlabel('(\theta)');
ylabel('Normalisasi Array Factor');
title('Array Respons BTS (M = 5)')

figure(2)
polar(ph0, f)
title('Pola Radiasi BTS (M = 5)')

