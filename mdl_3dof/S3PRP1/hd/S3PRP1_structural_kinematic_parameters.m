% Return Structural Kinematic Parameters of the Robot 
% S3PRP1
%
% Output:
% v_mdh [3x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [3x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [3x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x3) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = S3PRP1_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [1; 0; 1;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 1; 1;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 4;
% Anzahl der Kinematikparameter
NKP = 3;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 3;
% Namen der Kinematikparameter
pkin_names = {'a2', 'a3', 'd2'};
