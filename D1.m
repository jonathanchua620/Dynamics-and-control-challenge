open('MeasuredSignals-1.mat');

t = MeasuredSignals.t;

x=MeasuredSignals.x;

varphi=MeasuredSignals.varphi;

figure
plot(t, x);

figure
plot(t, varphi);



% Question 1a: Position and velocity vectors

r3 = [x + L23*cos(phi-phi_h3), ...
      y + L23*sin(phi-phi_h3)];

rdot3 = [xdot - L23*sin(phi-phi_h3)*phi_dot, ...
         ydot + L23*cos(phi-phi_h3)*phi_dot];

r4 = [x + L24*cos(phi), ...
      y + L24*sin(phi)];

rdot4 = [xdot - L24*sin(phi)*phi_dot, ...
         ydot + L24*cos(phi)*phi_dot];


% Question 1b: Kinetic energy

J3 = Jh2 - mh*L23^2;

T = 0.5*mh*(rdot3(1)^2 + rdot3(2)^2) + 0.5*J3*phi_dot^2 + 0.5*mb*ydot^2 + 0.5*mc*(xdot^2 + ydot^2);


% Question 1c: Potential energy

V = mb*g*y + mc*g*y + mh*g*(y + L23*sin(phi-phi_h3)) + 0.5*ky*(y-yref)^2 + 0.5*kPD*(x-xPDref)^2 + 0.5*k_phi*(phi-phi_ref)^2;


% Question 1d: Generalized non-potential forces

f_np = [-(dPD+dx)*xdot;
         M_phi-d_phi*phi_dot];

% Question 1e:

