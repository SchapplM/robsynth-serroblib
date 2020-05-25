% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:05
% EndTime: 2020-01-03 11:49:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (114->56), mult. (117->54), div. (0->0), fcn. (168->8), ass. (0->37)
t19 = sin(qJ(3));
t29 = -t19 * pkin(3) - qJ(2);
t23 = -pkin(7) - pkin(6);
t21 = cos(qJ(3));
t9 = t21 * pkin(3) + pkin(2);
t16 = qJ(3) + qJ(4);
t10 = sin(t16);
t17 = sin(pkin(8));
t39 = t17 * t10;
t18 = cos(pkin(8));
t20 = sin(qJ(1));
t38 = t20 * t18;
t37 = t20 * t19;
t36 = t20 * t21;
t22 = cos(qJ(1));
t35 = t22 * t17;
t34 = t22 * t18;
t33 = t22 * t19;
t32 = t22 * t21;
t31 = -pkin(4) * t10 + t29;
t15 = pkin(5) + 0;
t30 = t20 * pkin(1) + 0;
t28 = -t20 * qJ(2) + 0;
t27 = pkin(2) * t18 + pkin(6) * t17;
t14 = -qJ(5) + t23;
t11 = cos(t16);
t5 = pkin(4) * t11 + t9;
t26 = -t14 * t17 + t18 * t5;
t25 = -t17 * t23 + t18 * t9;
t24 = -t22 * qJ(2) + t30;
t8 = t20 * t17;
t7 = t17 * t11;
t4 = -t20 * t10 - t11 * t34;
t3 = t10 * t34 - t20 * t11;
t2 = -t22 * t10 + t11 * t38;
t1 = -t10 * t38 - t22 * t11;
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t15; t20, t22, 0, 0; -t22, t20, 0, 0; 0, 0, 0, 1; t17, t18, 0, t15; t38, -t8, -t22, t24; -t34, t35, -t20, -t22 * pkin(1) + t28; 0, 0, 0, 1; t17 * t21, -t17 * t19, -t18, t17 * pkin(2) - t18 * pkin(6) + t15; t18 * t36 - t33, -t18 * t37 - t32, t8, t27 * t20 + t24; -t18 * t32 - t37, t18 * t33 - t36, -t35, (-pkin(1) - t27) * t22 + t28; 0, 0, 0, 1; t7, -t39, -t18, t17 * t9 + t18 * t23 + t15; t2, t1, t8, t25 * t20 + t29 * t22 + t30; t4, t3, -t35, 0 + t29 * t20 + (-pkin(1) - t25) * t22; 0, 0, 0, 1; t7, -t39, -t18, t18 * t14 + t17 * t5 + t15; t2, t1, t8, t26 * t20 + t31 * t22 + t30; t4, t3, -t35, 0 + t31 * t20 + (-pkin(1) - t26) * t22; 0, 0, 0, 1;];
T_ges = t6;
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
