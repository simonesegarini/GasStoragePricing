function Vs = simulationV(params, dt, nSim)

% Sistemare

% Assign parameters.
alpha = params(1); b = params(2);
a = exp(-b*dt);

% Preallocate for memory saving.
Vs = zeros(nSim, 1);

% Define the PDF of W.
fW = @(w) log(a^(-alpha)) * (exp(w * log(a^(-alpha))) - 1) ...
    \ (a^(-alpha) - 1 - log(a^(-alpha)));

% Discretization of the w values for fW.
L = 100;
ws = linspace(0, 1, L+1);
fW_evaluated = fW(ws);

% Computation of the parts needed for the acceptance/rejection method.
qls = 0.5 * ( fW_evaluated(1:end-1) + fW_evaluated(2:end) ) / L;
GL = sum(qls);

pls = qls/GL;
pls_cum = cumsum(pls);

for it = 1:length(Vs)

    acceptance = false;
    % check better this part
    while ~acceptance
        
        s = find(pls_cum > rand(1), 1, 'first');

        coeff_ang = (fW_evaluated(s+1) - fW_evaluated(s)) / (ws(s+1) - ws(s));

        y = ws(s) - fW_evaluated(s) / coeff_ang + ...
            sqrt(fW_evaluated(s)^2 / coeff_ang^2 + 2 / coeff_ang * qls(s) * rand(1));

        u = rand(1);

        fw_y = fW(y);
        gl_y = (coeff_ang * (y-ws(s))) + fW_evaluated(s);

        acceptance = (u <= (fw_y/gl_y));
    end

    Vs(it) = a^(-y);
end