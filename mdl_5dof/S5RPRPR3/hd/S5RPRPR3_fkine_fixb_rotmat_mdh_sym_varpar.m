% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)
% T_c_stack [(5+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPRPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:29
% EndTime: 2022-01-23 09:20:30
% DurationCPUTime: 0.11s
% Computational Cost: add. (126->30), mult. (52->30), div. (0->0), fcn. (85->10), ass. (0->27)
t15 = sin(pkin(9));
t14 = qJ(1) + pkin(8);
t11 = qJ(3) + t14;
t6 = sin(t11);
t33 = t6 * t15;
t7 = cos(t11);
t32 = t7 * t15;
t16 = cos(pkin(9));
t17 = sin(qJ(5));
t31 = t16 * t17;
t19 = cos(qJ(5));
t30 = t16 * t19;
t29 = pkin(5) + 0;
t18 = sin(qJ(1));
t28 = t18 * pkin(1) + 0;
t20 = cos(qJ(1));
t27 = t20 * pkin(1) + 0;
t9 = sin(t14);
t26 = pkin(2) * t9 + t28;
t10 = cos(t14);
t25 = pkin(2) * t10 + t27;
t24 = qJ(2) + t29;
t8 = pkin(6) + t24;
t23 = pkin(4) * t16 + pkin(7) * t15;
t22 = t7 * pkin(3) + t6 * qJ(4) + t25;
t21 = t6 * pkin(3) - t7 * qJ(4) + t26;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t20, -t18, 0, 0; t18, t20, 0, 0; 0, 0, 1, t29; t10, -t9, 0, t27; t9, t10, 0, t28; 0, 0, 1, t24; t7, -t6, 0, t25; t6, t7, 0, t26; 0, 0, 1, t8; t7 * t16, -t32, t6, t22; t6 * t16, -t33, -t7, t21; t15, t16, 0, t8; t6 * t17 + t30 * t7, t6 * t19 - t31 * t7, t32, t23 * t7 + t22; -t7 * t17 + t30 * t6, -t7 * t19 - t31 * t6, t33, t23 * t6 + t21; t15 * t19, -t15 * t17, -t16, t15 * pkin(4) - t16 * pkin(7) + t8;];
Tc_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
