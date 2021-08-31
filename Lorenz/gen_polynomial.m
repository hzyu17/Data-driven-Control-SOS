function [C_a_Psi, C_cx_Psi, C_ab_Psi, C_bc_Psi, Psi] = gen_poly(x , bx, deg_a, deg_c, order, opt)
% Generate coefficients of polynomials in the basis Psi.

% Method: generate a common polynomial basis M(x), namely xmon in the code,
% for all the polynomials. For each polynomial, its coefficient in basis Psi is 
% calculated firstly in xmon, then we use the transition from xmon to Psi
% to express the polynomial in Psi.

% -----------------------------------------------
% HISTORY: 
% 12.2: start writing C_Psi_xmon an C_ax_Psi
% -----------------------------------------------

Psi = monomials(x, 0 : order);
nPsi = length(Psi);
nx = length(x);

% ----------------------------------------------------------
% Reorganize Psi in M(x): Find matrix C_Psi_xmon, s.t.:
%      Psi = C_Psi_xmon' * xmon
% ----------------------------------------------------------

% 1. Create a M_Psi(x) and find C_Psi_Psimon, s.t.:
%     Psi = C_Psi_Psimon' * M_Psi(x)
nPsimon = size(Psi.degmat, 1);
tmp = repmat(x', nPsimon, 1).^ Psi.degmat(:, end-nx+1:end);
Psimon = ones(nPsimon, 1);
for i1 = 1 : size(tmp, 2)
    Psimon = Psimon .* tmp(:, i1); % M_Psi(x)
end
C_Psi_Psimon = Psi.coefficient; % C_Psi_Psimon

%% 2. Reorganize Psi in M(x): Find C_Psi_xmon s.t.:
%        Psi = C_Psi_xmon' * xmon

% ---- a. Generate common polynomial basis, M(x) 
deg_xmon = 0:Psi.maxdeg;
xmon = monomials(x, deg_xmon); % M(x)
nxmon = length(xmon);

% ---- b. Reorganize
[~, ivarPsi] = ismember(xmon.varname, Psimon.varname);

c_mat_Psi_xmon = zeros(nxmon, nPsimon);
[ixmon, iPsimon] = ismember(xmon.degmat, Psimon.degmat(:, ivarPsi), 'row');
c_mat_Psi_xmon(find(ixmon==1),:) = Psimon.coefficient(nonzeros(iPsimon), :);
% Psimon = c_mat_Psi_xmon' * xmon
C_Psi_xmon = c_mat_Psi_xmon * C_Psi_Psimon;

%%
% --------------------------------------------------------------
% Reorganize a(x) in Psi by firstly reorganize a(x) in M(x)
% --------------------------------------------------------------
% ---- a. amon and C_a_amon
amon = monomials(x, 0:deg_a);
namon = length(amon);
C_a_amon = mpvar('a', namon, 1);
ax = C_a_amon' * amon;

% ---- b. Create matrix c_mat_ax_xmon, s.t.: 
%            amon = c_mat_ax_xmon' * xmon, and:
%     C_a_xmon = c_mat_ax_xmon * C_a_amon, and:
%                  ax = C_a_xmon' * xmon
c_mat_ax_xmon = zeros(nxmon, namon);
[ixmon_a, iamon] = ismember(xmon.degmat, amon.degmat, 'row');
c_mat_ax_xmon(find(ixmon_a==1),: ) = amon.coefficient(nonzeros(iamon), :);
C_a_xmon = c_mat_ax_xmon * C_a_amon;

% ---- c. Create C_a_Psi, s.t: 
%     ax = C_a_Psi' * Psi.

%     ax = C_a_Psi' * Psi = (C_a_Psi' *  C_Psi_xmon') * xmon 
%          = C_a_xmon' * xmon
C_a_Psi = s2p(C_Psi_xmon \ p2s(C_a_xmon));
C_a_Psi.coefficient(find(abs(C_a_Psi.coefficient)<=opt.round)) = 0;

%% 
% --------------------------
% Reorganize c(x) in Psi
% --------------------------

% 1. Create c_monomials and coefficients
cmon = monomials(x, 0:deg_c);
ncmon = length(cmon);
C_cx_cmon = mpvar('c', ncmon, 1);
cx = C_cx_cmon' * cmon;

% 2. Reorganize cx in xmon: find C_cx_xmon s.t.:
%       cx = C_cx_xmon' * xmon
c_mat_cx_xmon = zeros(nxmon, ncmon);
[ixmon_c, icmon] = ismember(xmon.degmat, cmon.degmat, 'row');
c_mat_cx_xmon(find(ixmon_c==1), :) = cmon.coefficient(nonzeros(icmon), :);
C_cx_xmon = c_mat_cx_xmon * C_cx_cmon; % cx = C_cx_xmon' * xmon

% 3. Reorganize cx in Psi
%   cx = C_cx_xmon' * xmon
%        = C_cx_Psi' * Psi
%        = C_cx_Psi' * C_Psi_xmon' * xmon
C_cx_Psi = s2p(C_Psi_xmon \ p2s(C_cx_xmon)); % cx = C_cx_Psi' * Psi
C_cx_Psi.coefficient(find(abs(C_cx_Psi.coefficient)<=opt.round)) = 0;

%% 
% ---------------------------
% Reorganize ab(x) in Psi
% ---------------------------
abx = bx * ax;
% 1. Create M_ab (abmon)
nabmon = size(abx.degmat, 1); % number of monomials in M_ab
[ivarxmon, ivarabx] = ismember(xmon.varname, abx.varname);

temp = repmat(x', nabmon, 1) .^ abx.degmat(:, ivarabx);
abmon = ones(nabmon, 1);
for i1 = 1:size(temp, 2)
    abmon = abmon .* temp(:, i1);
end

% 2. Find the coefficient variables for abx in abmon
coeff_vars = [];
for i1 = 1 : ivarabx(1)-1
    coeff_vars = [coeff_vars; s2p(str2sym(abx.varname{i1}))];
end
for i1 = ivarabx(end)+1 : length(abx.varname)
    coeff_vars = [coeff_vars; s2p(str2sym(abx.varname{i1}))];
end
degmat_var = [abx.degmat(:, 1:ivarabx(1)-1), abx.degmat(:, ivarabx(end)+1 : end)];

% C_abx_abmon = [ 1       [ 0 1 0       [ a2
%                              1   .*    2 0 0   *    a1
%                              1]        0 0 3]        a3]
%                         = diag([a2, a1, a3])
C_abx_abmon = abx.coefficient .* degmat_var * coeff_vars;
% abx = C_abx_abmon' * abmon

% 3. Reorganize abx in xmon
% reorganize abmon in xmon
[idegxmon, idegabx] = ismember(xmon.degmat(:, find(ivarxmon==1)), abmon.degmat, 'row');
c_mat_ab_xmon = zeros(nxmon, nabmon);
c_mat_ab_xmon(find(idegxmon==1), :) = abmon.coefficient(nonzeros(idegabx), :);
C_ab_xmon = c_mat_ab_xmon * C_abx_abmon; % abx = C_abx_xmon' * xmon

% 4. From C_abx_xmon to C_abx_Psi 
C_ab_Psi = s2p(C_Psi_xmon \ p2s(C_ab_xmon));
C_ab_Psi.coefficient(find(abs(C_ab_Psi.coefficient)<=opt.round)) = 0;
% abx = C_abx_Psi' * Psi;

%%
% ----------------------------------------------------------
% Reorganize bc(x) in Psi: similar to the case with a(x)
% ----------------------------------------------------------
bcx = bx * cx;
% 1. Create M_ab (abmon)
nbcmon = size(bcx.degmat, 1); % number of monomials in M_ab
[ivarxmon, ivarbcx] = ismember(xmon.varname, bcx.varname);

temp = repmat(x', nbcmon, 1) .^ bcx.degmat(:, ivarbcx);
bcmon = ones(nbcmon, 1);
for i1 = 1:size(temp, 2)
    bcmon = bcmon .* temp(:, i1);
end

% 2. Find the coefficient variables for abx in abmon
coeff_vars = [];
for i1 = 1 : ivarbcx(1)-1
    coeff_vars = [coeff_vars; s2p(str2sym(bcx.varname{i1}))];
end
for i1 = ivarbcx(end)+1 : length(bcx.varname)
    coeff_vars = [coeff_vars; s2p(str2sym(bcx.varname{i1}))];
end
degmat_var = [bcx.degmat(:, 1:ivarbcx(1)-1), bcx.degmat(:, ivarbcx(end)+1 : end)];

C_bc_bcmon = bcx.coefficient .* degmat_var * coeff_vars;
% bcx - C_bc_bcmon' * bcmon

% 3. Reorganize abx in xmon
% reorganize abmon in xmon
[idegxmon, idegbcx] = ismember(xmon.degmat(:, find(ivarxmon==1)), bcmon.degmat, 'row');
c_mat_bc_xmon = zeros(nxmon, nbcmon);
c_mat_bc_xmon(find(idegxmon==1), :) = bcmon.coefficient(nonzeros(idegbcx), :);
C_bc_xmon = c_mat_bc_xmon * C_bc_bcmon; % abx = C_abx_xmon' * xmon

% 4. From C_abx_xmon to C_abx_Psi 
C_bc_Psi = s2p(C_Psi_xmon \ p2s(C_bc_xmon));
C_bc_Psi.coefficient(find(abs(C_bc_Psi.coefficient)<=opt.round)) = 0;
% bcx - C_bc_Psi' * Psi

%% 
% ----------------------
% Reorganize x in Psi
% ----------------------
cmat_x_xmon = zeros(nxmon, nx);
[~, ivarinx] = ismember(xmon.varname, x.varname);
[ixmon_x, ix] = ismember(xmon.degmat, x.degmat(:, ivarinx), 'row');
cmat_x_xmon(find(ixmon_x==1), :) = x.coefficient(nonzeros(ix), :);
C_x_Psi = C_Psi_xmon\cmat_x_xmon;
C_x_Psi(find(C_x_Psi<=opt.round)) = 0;

end
