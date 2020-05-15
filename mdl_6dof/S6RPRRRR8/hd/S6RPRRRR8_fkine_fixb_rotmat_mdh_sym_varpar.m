% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRRRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:37:17
% EndTime: 2018-11-23 16:37:18
% DurationCPUTime: 0.14s
% Computational Cost: add. (141->56), mult. (110->56), div. (0->0), fcn. (160->10), ass. (0->40)
t17 = sin(qJ(1));
t14 = qJ(3) + qJ(4);
t5 = sin(t14);
t44 = t17 * t5;
t13 = qJ(5) + qJ(6);
t6 = cos(t13);
t43 = t17 * t6;
t7 = cos(t14);
t42 = t17 * t7;
t20 = cos(qJ(1));
t41 = t20 * t5;
t40 = t20 * t6;
t15 = sin(qJ(5));
t39 = t17 * t15;
t16 = sin(qJ(3));
t38 = t17 * t16;
t18 = cos(qJ(5));
t37 = t17 * t18;
t36 = t20 * t15;
t35 = t20 * t18;
t12 = pkin(6) + 0;
t34 = t17 * pkin(1) + 0;
t33 = pkin(2) + t12;
t22 = -pkin(8) - pkin(7);
t32 = pkin(5) * t15 - t22;
t31 = t20 * pkin(1) + t17 * qJ(2) + 0;
t30 = -pkin(3) * t16 - qJ(2);
t19 = cos(qJ(3));
t29 = t19 * pkin(3) + t33;
t28 = pkin(3) * t38 + t31;
t27 = pkin(4) * t5 - pkin(9) * t7;
t21 = -pkin(10) - pkin(9);
t3 = t18 * pkin(5) + pkin(4);
t26 = t21 * t7 + t3 * t5;
t25 = -t17 * t22 + t34;
t24 = -t20 * qJ(2) + t34;
t23 = -t20 * t22 + t28;
t4 = sin(t13);
t1 = t20 * t7;
t2 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t20, -t17, 0, 0; t17, t20, 0, 0; 0, 0, 1, t12; 0, 0, 0, 1; 0, -t20, t17, t31; 0, -t17, -t20, t24; 1, 0, 0, t12; 0, 0, 0, 1; t38, t17 * t19, t20, t20 * pkin(7) + t31; -t20 * t16, -t20 * t19, t17, t17 * pkin(7) + t24; t19, -t16, 0, t33; 0, 0, 0, 1; t44, t42, t20, t23; -t41, -t1, t17, t30 * t20 + t25; t7, -t5, 0, t29; 0, 0, 0, 1; t5 * t37 + t36, -t5 * t39 + t35, -t42, t27 * t17 + t23; -t5 * t35 + t39, t5 * t36 + t37, t1 (-t27 + t30) * t20 + t25; t7 * t18, -t7 * t15, t5, t7 * pkin(4) + t5 * pkin(9) + t29; 0, 0, 0, 1; t20 * t4 + t5 * t43, -t4 * t44 + t40, -t42, t26 * t17 + t32 * t20 + t28; t17 * t4 - t5 * t40, t4 * t41 + t43, t1, t32 * t17 + (-t26 + t30) * t20 + t34; t7 * t6, -t7 * t4, t5, -t5 * t21 + t7 * t3 + t29; 0, 0, 0, 1;];
T_ges = t2;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
