% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S1P1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% Tc_mdh [4x4x(1+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   2:  mdh base (link 0) -> mdh frame (2-1), link (2-1)
%   ...
%   1+1:  mdh base (link 0) -> mdh frame (1)
% T_c_stack [(1+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S1P1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:45
% EndTime: 2021-01-14 12:21:46
% DurationCPUTime: 0.10s
% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, qJ(1) + 0;];
Tc_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,1+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,1+1]); end % symbolisch
for i = 1:1+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
