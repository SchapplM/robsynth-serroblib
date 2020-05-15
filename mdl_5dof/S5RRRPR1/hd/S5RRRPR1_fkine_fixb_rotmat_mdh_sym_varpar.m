% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:41
% EndTime: 2019-12-05 18:37:41
% DurationCPUTime: 0.09s
% Computational Cost: add. (112->42), mult. (50->38), div. (0->0), fcn. (86->10), ass. (0->23)
t24 = -pkin(7) - pkin(6);
t22 = cos(qJ(2));
t10 = t22 * pkin(2) + pkin(1);
t19 = qJ(2) + qJ(3);
t18 = pkin(5) + 0;
t13 = cos(t19);
t2 = pkin(3) * t13 + t10;
t17 = -qJ(4) + t24;
t20 = sin(qJ(2));
t26 = t20 * pkin(2) + t18;
t11 = pkin(9) + t19;
t12 = sin(t19);
t25 = pkin(3) * t12 + t26;
t23 = cos(qJ(1));
t21 = sin(qJ(1));
t14 = -pkin(8) + t17;
t9 = qJ(5) + t11;
t6 = cos(t11);
t5 = sin(t11);
t4 = cos(t9);
t3 = sin(t9);
t1 = pkin(4) * t6 + t2;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t21, 0, 0; t21, t23, 0, 0; 0, 0, 1, t18; 0, 0, 0, 1; t23 * t22, -t23 * t20, t21, t23 * pkin(1) + t21 * pkin(6) + 0; t21 * t22, -t21 * t20, -t23, t21 * pkin(1) - t23 * pkin(6) + 0; t20, t22, 0, t18; 0, 0, 0, 1; t23 * t13, -t23 * t12, t21, t23 * t10 - t21 * t24 + 0; t21 * t13, -t21 * t12, -t23, t21 * t10 + t23 * t24 + 0; t12, t13, 0, t26; 0, 0, 0, 1; t23 * t6, -t23 * t5, t21, -t21 * t17 + t23 * t2 + 0; t21 * t6, -t21 * t5, -t23, t23 * t17 + t21 * t2 + 0; t5, t6, 0, t25; 0, 0, 0, 1; t23 * t4, -t23 * t3, t21, t23 * t1 - t21 * t14 + 0; t21 * t4, -t21 * t3, -t23, t21 * t1 + t23 * t14 + 0; t3, t4, 0, pkin(4) * t5 + t25; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
