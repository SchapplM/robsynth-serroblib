% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR2
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RPRPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:42
% EndTime: 2022-01-23 09:18:43
% DurationCPUTime: 0.09s
% Computational Cost: add. (112->30), mult. (33->22), div. (0->0), fcn. (61->10), ass. (0->23)
t26 = pkin(5) + 0;
t15 = qJ(1) + pkin(8);
t19 = sin(qJ(1));
t25 = t19 * pkin(1) + 0;
t20 = cos(qJ(1));
t24 = t20 * pkin(1) + 0;
t8 = sin(t15);
t23 = pkin(2) * t8 + t25;
t10 = cos(t15);
t22 = pkin(2) * t10 + t24;
t21 = qJ(2) + t26;
t6 = pkin(6) + t21;
t18 = -pkin(7) - qJ(4);
t17 = cos(pkin(9));
t16 = sin(pkin(9));
t14 = pkin(9) + qJ(5);
t11 = qJ(3) + t15;
t9 = cos(t14);
t7 = sin(t14);
t5 = pkin(4) * t17 + pkin(3);
t4 = cos(t11);
t3 = sin(t11);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t20, -t19, 0, 0; t19, t20, 0, 0; 0, 0, 1, t26; t10, -t8, 0, t24; t8, t10, 0, t25; 0, 0, 1, t21; t4, -t3, 0, t22; t3, t4, 0, t23; 0, 0, 1, t6; t4 * t17, -t4 * t16, t3, pkin(3) * t4 + qJ(4) * t3 + t22; t3 * t17, -t3 * t16, -t4, pkin(3) * t3 - qJ(4) * t4 + t23; t16, t17, 0, t6; t4 * t9, -t4 * t7, t3, -t18 * t3 + t4 * t5 + t22; t3 * t9, -t3 * t7, -t4, t18 * t4 + t3 * t5 + t23; t7, t9, 0, pkin(4) * t16 + t6;];
Tc_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
