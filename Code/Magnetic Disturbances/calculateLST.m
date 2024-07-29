function LST = calculateLST(date, lon)

GMT = mod(date2mjd2000(date)-0.5,1)*2*pi;

LST = GMT + lon;

%{
jd = date2jd(date);

date0 = [date(1) date(2) date(3) 0 0 0];

jd0 = date2jd(date0);

H = (jd - jd0) * 24;

Dtt = jd - 2451545;

Dut = jd0 - 2451545;

T = Dtt/36525;
    
GMST = mod(6.697375 + 0.065709824279*Dut + 1.0027379*H + 0.0000258*T^2, 24)*360/24;

LST = GMST + lon;
%}
end