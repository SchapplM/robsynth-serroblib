% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 15:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPRRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:51:37
% EndTime: 2018-11-23 15:51:37
% DurationCPUTime: 0.14s
% Computational Cost: add. (141->56), mult. (110->56), div. (0->0), fcn. (160->10), ass. (0->40)
t19 = sin(qJ(1));
t12 = pkin(10) + qJ(4);
t5 = cos(t12);
t44 = t19 * t5;
t14 = qJ(5) + qJ(6);
t6 = sin(t14);
t43 = t19 * t6;
t7 = cos(t14);
t42 = t19 * t7;
t21 = cos(qJ(1));
t41 = t21 * t6;
t40 = t21 * t7;
t15 = sin(pkin(10));
t39 = t19 * t15;
t18 = sin(qJ(5));
t38 = t19 * t18;
t20 = cos(qJ(5));
t37 = t19 * t20;
t36 = t21 * t18;
t35 = t21 * t20;
t13 = pkin(6) + 0;
t34 = t19 * pkin(1) + 0;
t33 = pkin(2) + t13;
t17 = -pkin(7) - qJ(3);
t32 = pkin(5) * t18 - t17;
t31 = t21 * pkin(1) + t19 * qJ(2) + 0;
t30 = -pkin(3) * t15 - qJ(2);
t16 = cos(pkin(10));
t29 = t16 * pkin(3) + t33;
t28 = pkin(3) * t39 + t31;
t4 = sin(t12);
t27 = pkin(4) * t4 - pkin(8) * t5;
t22 = -pkin(9) - pkin(8);
t3 = t20 * pkin(5) + pkin(4);
t26 = t22 * t5 + t3 * t4;
t25 = -t19 * t17 + t34;
t24 = -t21 * qJ(2) + t34;
t23 = -t21 * t17 + t28;
t1 = t21 * t5;
t2 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t19, 0, 0; t19, t21, 0, 0; 0, 0, 1, t13; 0, 0, 0, 1; 0, -t21, t19, t31; 0, -t19, -t21, t24; 1, 0, 0, t13; 0, 0, 0, 1; t39, t19 * t16, t21, t21 * qJ(3) + t31; -t21 * t15, -t21 * t16, t19, t19 * qJ(3) + t24; t16, -t15, 0, t33; 0, 0, 0, 1; t19 * t4, t44, t21, t23; -t21 * t4, -t1, t19, t30 * t21 + t25; t5, -t4, 0, t29; 0, 0, 0, 1; t4 * t37 + t36, -t4 * t38 + t35, -t44, t27 * t19 + t23; -t4 * t35 + t38, t4 * t36 + t37, t1 (-t27 + t30) * t21 + t25; t5 * t20, -t5 * t18, t4, t5 * pkin(4) + t4 * pkin(8) + t29; 0, 0, 0, 1; t4 * t42 + t41, -t4 * t43 + t40, -t44, t26 * t19 + t32 * t21 + t28; -t4 * t40 + t43, t4 * t41 + t42, t1, t32 * t19 + (-t26 + t30) * t21 + t34; t5 * t7, -t5 * t6, t4, -t4 * t22 + t5 * t3 + t29; 0, 0, 0, 1;];
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
