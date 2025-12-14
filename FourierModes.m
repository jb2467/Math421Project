function [Field] = Fourier_Modes(sx,sy,sz,lx,ly,lz,m,N_modes,init,Amp,isVector)

% Create spatial grid
x = linspace(0, lx, sx);
y = linspace(0, ly, sy);
z = linspace(0, lz, sz);
[X, Y, Z] = ndgrid(x, y, z);

% Grid spacing
delta = lx / (sx - 1);

% Wave numbers
alpha = 1.453;
L_t = 0.1 * lx;
p = 2;
k_e = alpha * 9 * pi / (55 * L_t);
k_min = k_e / p;
k_max = 2 * pi / (2 * delta);
k_values = linspace(k_min, k_max, N_modes);

% Random phases for all modes
psi = 2*pi*rand(1, N_modes);

if isVector
    Field = zeros(sx, sy, sz, m, 3);  % 3 for Ux, Uy, Uz
else
    Field = zeros(sx, sy, sz, m);
end

for n = 1:N_modes
    k = k_values(n);

    theta = 2*pi*rand();
    phi = acos(2*rand()-1);
    dir = [sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)];

    phase = k*(dir(1)*X + dir(2)*Y + dir(3)*Z);

    for t_i = 1:m
        if isVector
            Field(:,:,:,t_i,1) = Field(:,:,:,t_i,1) + Amp * dir(1) * cos(phase + psi(n) - k*t_i);
            Field(:,:,:,t_i,2) = Field(:,:,:,t_i,2) + Amp * dir(2) * cos(phase + psi(n) - k*t_i);
            Field(:,:,:,t_i,3) = Field(:,:,:,t_i,3) + Amp * dir(3) * cos(phase + psi(n) - k*t_i);
        else
            Field(:,:,:,t_i) = Field(:,:,:,t_i) + Amp * cos(phase + psi(n) - k*t_i);
        end
    end
end

if isVector
    for d = 1:3
        Field(:,:,:,:,d) = Field(:,:,:,:,d) + init(d);
    end
else
    Field = Field + init;
end

end





sx = 100; sy = 1; sz = 1;
lx = 10; ly = 1; lz = 1;
m = 2000;
N_modes = 1200;

T0 = 300;
AmpT = 1;
P0 = 101325;
AmpP = 500;
Rho0 = 1.225;
AmpRho = 0.05;

T = Fourier_Modes(sx, sy, sz, lx, ly, lz, m, N_modes, T0, AmpT, false);
P = Fourier_Modes(sx, sy, sz, lx, ly, lz, m, N_modes, P0, AmpP, false);
Rho = Fourier_Modes(sx, sy, sz, lx, ly, lz, m, N_modes, Rho0, AmpRho, false);

BaselineVel = [0, 0, 0];
AmpVel = 1.0;             % velocity fluctuation amplitude
U = Fourier_Modes(sx, sy, sz, lx, ly, lz, m, N_modes, BaselineVel, AmpVel, true);

save("/MATLAB Drive/MATH 421/Data/T.mat", "T");
save("/MATLAB Drive/MATH 421/Data/P.mat", "P");
save("/MATLAB Drive/MATH 421/Data/Rho.mat", "Rho");
save("/MATLAB Drive/MATH 421/Data/U.mat", "U");
