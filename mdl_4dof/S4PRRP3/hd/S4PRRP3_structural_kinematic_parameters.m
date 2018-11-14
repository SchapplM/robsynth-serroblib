% Return Structural Kinematic Parameters of the Robot 
% S4PRRP3
%
% Output:
% v_mdh [4x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [4x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [4x1]
%   Aktives Gelenk = 1, Passiv = 0

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh] = S4PRRP3_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2; 3;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [1; 0; 0; 1;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 1; 1; 1;];
mu_mdh = t1;
