% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR12
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
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRPR12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:19
% EndTime: 2019-12-31 18:29:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (121->48), mult. (92->54), div. (0->0), fcn. (138->10), ass. (0->32)
t14 = pkin(8) + qJ(3);
t11 = cos(t14);
t22 = sin(qJ(1));
t35 = t22 * t11;
t16 = sin(pkin(9));
t34 = t22 * t16;
t18 = cos(pkin(9));
t33 = t22 * t18;
t23 = cos(qJ(1));
t32 = t23 * t11;
t31 = t23 * t16;
t30 = t23 * t18;
t15 = pkin(5) + 0;
t19 = cos(pkin(8));
t6 = t19 * pkin(2) + pkin(1);
t29 = t23 * t6 + 0;
t21 = -pkin(6) - qJ(2);
t28 = t23 * t21 + t22 * t6 + 0;
t17 = sin(pkin(8));
t27 = t17 * pkin(2) + t15;
t20 = -pkin(7) - qJ(4);
t5 = t18 * pkin(4) + pkin(3);
t9 = sin(t14);
t26 = t11 * t5 - t20 * t9;
t25 = pkin(3) * t11 + qJ(4) * t9;
t24 = -t22 * t21 + t29;
t13 = pkin(9) + qJ(5);
t10 = cos(t13);
t8 = sin(t13);
t4 = t23 * t9;
t3 = t22 * t9;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t23, -t22, 0, 0; t22, t23, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t23 * t19, -t23 * t17, t22, t23 * pkin(1) + t22 * qJ(2) + 0; t22 * t19, -t22 * t17, -t23, t22 * pkin(1) - t23 * qJ(2) + 0; t17, t19, 0, t15; 0, 0, 0, 1; t32, -t4, t22, t24; t35, -t3, -t23, t28; t9, t11, 0, t27; 0, 0, 0, 1; t11 * t30 + t34, -t11 * t31 + t33, t4, t25 * t23 + t24; t11 * t33 - t31, -t11 * t34 - t30, t3, t25 * t22 + t28; t9 * t18, -t9 * t16, -t11, t9 * pkin(3) - t11 * qJ(4) + t27; 0, 0, 0, 1; t10 * t32 + t22 * t8, t22 * t10 - t8 * t32, t4, t26 * t23 + (pkin(4) * t16 - t21) * t22 + t29; t10 * t35 - t23 * t8, -t23 * t10 - t8 * t35, t3, -pkin(4) * t31 + t26 * t22 + t28; t9 * t10, -t9 * t8, -t11, t11 * t20 + t9 * t5 + t27; 0, 0, 0, 1;];
T_ges = t1;
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
