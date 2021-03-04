% Return Structural Kinematic Parameters of the Robot 
% S2PP1
%
% Output:
% v_mdh [2x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [2x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [2x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x1) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2021-03-03 18:41
% Revision: 33b345ae0dd6ec4aa15499ab3d43edbbded0bea5 (2021-02-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = S2PP1_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [1; 1;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 1;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 3;
% Anzahl der Kinematikparameter
NKP = 1;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 2;
% Namen der Kinematikparameter
pkin_names = {'a2'};
