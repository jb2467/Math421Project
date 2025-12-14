function f = computeBGKf(TemperatureDim, Ux, Rho, tau, v_grid, dt, dx)
m_air = 4.81e-26;
kB = 1.38e-23;

Nx = size(TemperatureDim,1);
Nt = size(TemperatureDim,4);
Nv = numel(v_grid);

dv = v_grid(2) - v_grid(1);

f = zeros(Nx, Nv, Nt);

for t = 1:Nt

    T   = squeeze(TemperatureDim(:,:,:,t));
    u   = squeeze(Ux(:,:,:,t));
    rho = squeeze(Rho(:,:,:,t));

    dTdx = zeros(Nx,1);
    dudx = zeros(Nx,1);
    drhodx = zeros(Nx,1);

    dTdx(2:Nx-1)   = (T(3:Nx)   - T(1:Nx-2))   / (2*dx);
    dudx(2:Nx-1)   = (u(3:Nx)   - u(1:Nx-2))   / (2*dx);
    drhodx(2:Nx-1) = (rho(3:Nx) - rho(1:Nx-2)) / (2*dx);

    dTdx(1) = (T(2)-T(1))/dx;       dTdx(Nx) = (T(Nx)-T(Nx-1))/dx;
    dudx(1) = (u(2)-u(1))/dx;       dudx(Nx) = (u(Nx)-u(Nx-1))/dx;
    drhodx(1) = (rho(2)-rho(1))/dx; drhodx(Nx) = (rho(Nx)-rho(Nx-1))/dx;

    if t < Nt
        dTdt   = (squeeze(TemperatureDim(:,:,:,t+1)) - T)   / dt;
        dudt   = (squeeze(Ux(:,:,:,t+1))            - u)   / dt;
        drhodt = (squeeze(Rho(:,:,:,t+1))           - rho) / dt;
    else
        dTdt   = (T   - squeeze(TemperatureDim(:,:,:,t-1))) / dt;
        dudt   = (u   - squeeze(Ux(:,:,:,t-1)))             / dt;
        drhodt = (rho - squeeze(Rho(:,:,:,t-1)))            / dt;
    end

    for ix = 1:Nx

        rho_loc = rho(ix);
        u_loc   = u(ix);
        T_loc   = T(ix);

        c = v_grid - u_loc;

        f0 = rho_loc ./ sqrt(2*pi*kB*T_loc/m_air) .* ...
            exp(-m_air*c.^2 / (2*kB*T_loc));

        df0dt = f0 .* ( ...
            drhodt(ix)/rho_loc + ...
            (m_air*c.^2/(2*kB*T_loc^2) - 1/(2*T_loc)) * dTdt(ix) + ...
            (m_air*c/(kB*T_loc)) * dudt(ix) );

        df0dx = f0 .* ( ...
            drhodx(ix)/rho_loc + ...
            (m_air*c.^2/(2*kB*T_loc^2) - 1/(2*T_loc)) * dTdx(ix) + ...
            (m_air*c/(kB*T_loc)) * dudx(ix) );

        f1 = -tau * (df0dt + v_grid .* df0dx);

        mass_f1 = sum(f1) * dv;
        f1 = f1 - mass_f1 * f0 / rho_loc;

        f(ix,:,t) = f0 + f1;
    end
end

end
