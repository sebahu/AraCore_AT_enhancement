function [ minFluxes, maxFluxes, exitflag, opt_obj ] = slowFVA( in_model, objectivePart )
exitflag = 0;
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
model.lb(obj_rxn_i) = objectivePart*opt_obj;

for rxn_i = 1:size(in_model.S,2)
    f = zeros(size(in_model.S,2),1);
    if model.lb(rxn_i) < minFluxes(rxn_i)
        f(rxn_i)=1;
        model.obj = f;
        result = gurobi(model, params);

        if ~strcmp(result.status, 'OPTIMAL')
            exitflag = -1;
            return;
        end
        fluxes = result.x;
        minFluxes = min(minFluxes, fluxes);
    end
    if model.ub(rxn_i) > maxFluxes(rxn_i)
        f(rxn_i)=-1;
        model.obj = f;
        result = gurobi(model, params);

        if ~strcmp(result.status, 'OPTIMAL')
            exitflag = -1;
            return;
        end
        fluxes = result.x;
        maxFluxes = max(maxFluxes, fluxes);
    end    
end



