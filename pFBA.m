function [ fluxes, exitflag, opt_obj ] = pFBA( in_model, minimizeObj )
exitflag = 0;
fluxes = [];
opt_obj = NaN;
f = zeros(size(in_model.S,2),1);
obj_rxn_i = find(in_model.c ~= 0);
f(obj_rxn_i)=-1;
beq = zeros(size(in_model.S,1),1);

options = [];optimoptions('linprog','Display','none');


model.obj = f;
model.A = sparse(in_model.S); % A must be sparse
n = size(model.A, 2);
model.vtype = repmat('C', n, 1);
model.sense = repmat('=',size(in_model.S,1),1);
model.rhs = full(beq(:)); % rhs must be dense
model.lb = in_model.lb;
model.ub = in_model.ub;


% Extract relevant Gurobi parameters from (subset of) options
params = struct();
params.OutputFlag = 0;
%params.MIPGap = options.RelativeGapTolerance;
%params.MIPGapAbs = options.AbsoluteGapTolerance;
params.FeasibilityTol = 1E-9;
params.OptimalityTol = 1E-9;
params.NumericFocus=3;
%params.IntFeasTol = options.IntFeasTol;
params.Quad = 1;
% Solve model with Gurobi
result = gurobi(model, params);

if ~strcmp(result.status, 'OPTIMAL')
    exitflag = -1;
    return;
end

fluxes = result.x;
val = result.objval;
opt_obj = -val;
minFluxes = fluxes;
maxFluxes = fluxes;

model.lb = in_model.lb;
model.lb(obj_rxn_i) = opt_obj;
if isempty(minimizeObj)
    model.obj = ones(size(in_model.S,2),1);
else
    model.obj = minimizeObj;
end
result = gurobi(model, params);
if ~strcmp(result.status, 'OPTIMAL')
    exitflag = -1;
    return;
end
fluxes = result.x;



