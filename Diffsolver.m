function xdot = Diffsolver(t,x)
global R L v p J Tl k frame f
tee = t;

%Phase Voltage: No initial phase for Vas. A-B-C sequence for Voltage
vas = sqrt(2/3)*v*cos(2*pi*f*tee);
vbs = sqrt(2/3)*v*cos(2*pi*f*tee-2*pi/3);
vcs = sqrt(2/3)*v*cos(2*pi*f*tee+2*pi/3);
if frame == 1
    rot = exp(-1i*2*pi*f*tee);
    a = exp(1i*(2*pi/3));
    b = exp(1i*(-2*pi/3));
    vqs = real((2/3)*rot*(vas + a*vbs + b*vcs));
    vds = -imag((2/3)*rot*(vas + a*vbs + b*vcs));
    V = [vds;vqs;0;0];
    W = zeros(4,4); W(1,2) = -2*pi*f; W(2,1) = 2*pi*f; W(3,4) = -(2*pi*f-x(5,1));  W(4,3) = 2*pi*f -x(5,1);
end

if frame == 2
    a = exp(1i*(2*pi/3));
    b = exp(1i*(-2*pi/3));
    vqs = real((2/3)*(vas + a*vbs + b*vcs));
    vds = -imag((2/3)*(vas + a*vbs + b*vcs));
    V = [vds;vqs;0;0];
    W = zeros(4,4);  W(3,4) = x(5,1);  W(4,3) = -x(5,1);
end

if frame == 3
    rot = exp(-1i*x(6,1));
    a = exp(1i*(2*pi/3));
    b = exp(1i*(-2*pi/3));
    vqs = real((2/3)*rot*(vas + a*vbs + b*vcs));
    vds = -imag((2/3)*rot*(vas + a*vbs + b*vcs));
    V = [vds;vqs;0;0];
    W = zeros(4,4); W(1,2) = -x(5,1); W(2,1) = x(5,1); 
end

xdot = zeros(6,1);
xdot(1:4,1) = V-R/L*x(1:4,1)-W*x(1:4,1);
Te = 0.75*p*k*(x(2,1)*x(3,1)-x(4,1)*x(1,1));
xdot(5,1) = (p/(2*J))*(Te-Tl);
xdot(6,1) = x(5,1);