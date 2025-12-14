function T = solveHeat1D(T0, rho_mass, q, k, cv, dx, dt)

[Nx, Nt] = size(q);
x = (0:Nx-1)' * dx;
t = (0:Nt-1) * dt;

dqdx = zeros(Nx, Nt);
for n = 1:Nt
    dqdx(2:Nx-1,n) = (q(3:Nx,n) - q(1:Nx-2,n)) / (2*dx);
    dqdx(1,n)      = (q(2,n) - q(1,n)) / dx;
    dqdx(Nx,n)     = (q(Nx,n) - q(Nx-1,n)) / dx;
end

rho_itp  = griddedInterpolant({x,t}, rho_mass, 'linear', 'nearest');
dqdx_itp = griddedInterpolant({x,t}, dqdx,     'linear', 'nearest');

m = 0;
pde = @(xpos,tval,u,dudx) pdefun(xpos,tval,u,dudx,...
    rho_itp,dqdx_itp,k,cv);
ic  = @(xpos) interp1(x, T0, xpos, 'linear', T0(1));
bc  = @(xl,ul,xr,ur,t) deal(0,1,0,1);

sol = pdepe(m, pde, ic, bc, x, t);

T = squeeze(sol)';

end


function [c,f,s] = pdefun(xpos,tval,u,dudx,rho_itp,dqdx_itp,k,cv)
rho = rho_itp(xpos, tval);
c = 1;
f = (k / (rho * cv + eps)) * dudx;
s = dqdx_itp({xpos, tval}) / (rho * cv + eps);
end
