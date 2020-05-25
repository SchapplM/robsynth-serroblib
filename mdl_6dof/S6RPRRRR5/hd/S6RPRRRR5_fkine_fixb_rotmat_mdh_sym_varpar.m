% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:35
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRRRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:34:55
% EndTime: 2018-11-23 16:34:55
% DurationCPUTime: 0.14s
% Computational Cost: add. (200->59), mult. (112->64), div. (0->0), fcn. (166->12), ass. (0->40)
t24 = cos(pkin(11));
t10 = t24 * pkin(2) + pkin(1);
t22 = qJ(5) + qJ(6);
t15 = sin(t22);
t27 = sin(qJ(1));
t45 = t27 * t15;
t16 = cos(t22);
t44 = t27 * t16;
t26 = sin(qJ(5));
t43 = t27 * t26;
t28 = cos(qJ(5));
t42 = t27 * t28;
t29 = cos(qJ(1));
t41 = t29 * t15;
t40 = t29 * t16;
t39 = t29 * t26;
t38 = t29 * t28;
t25 = -pkin(7) - qJ(2);
t21 = pkin(6) + 0;
t20 = pkin(11) + qJ(3);
t13 = cos(t20);
t3 = pkin(3) * t13 + t10;
t37 = t29 * t3 + 0;
t19 = -pkin(8) + t25;
t36 = t29 * t19 + t27 * t3 + 0;
t23 = sin(pkin(11));
t35 = t23 * pkin(2) + t21;
t12 = sin(t20);
t34 = pkin(3) * t12 + t35;
t14 = qJ(4) + t20;
t8 = sin(t14);
t9 = cos(t14);
t33 = pkin(4) * t9 + pkin(9) * t8;
t11 = t28 * pkin(5) + pkin(4);
t30 = -pkin(10) - pkin(9);
t32 = t11 * t9 - t30 * t8;
t31 = -t27 * t19 + t37;
t5 = t29 * t8;
t4 = t27 * t8;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t27, 0, 0; t27, t29, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t29 * t24, -t29 * t23, t27, t29 * pkin(1) + t27 * qJ(2) + 0; t27 * t24, -t27 * t23, -t29, t27 * pkin(1) - t29 * qJ(2) + 0; t23, t24, 0, t21; 0, 0, 0, 1; t29 * t13, -t29 * t12, t27, t29 * t10 - t27 * t25 + 0; t27 * t13, -t27 * t12, -t29, t27 * t10 + t29 * t25 + 0; t12, t13, 0, t35; 0, 0, 0, 1; t29 * t9, -t5, t27, t31; t27 * t9, -t4, -t29, t36; t8, t9, 0, t34; 0, 0, 0, 1; t9 * t38 + t43, -t9 * t39 + t42, t5, t29 * t33 + t31; t9 * t42 - t39, -t9 * t43 - t38, t4, t27 * t33 + t36; t8 * t28, -t8 * t26, -t9, t8 * pkin(4) - t9 * pkin(9) + t34; 0, 0, 0, 1; t9 * t40 + t45, -t9 * t41 + t44, t5, t32 * t29 + (pkin(5) * t26 - t19) * t27 + t37; t9 * t44 - t41, -t9 * t45 - t40, t4, -pkin(5) * t39 + t27 * t32 + t36; t8 * t16, -t8 * t15, -t9, t8 * t11 + t9 * t30 + t34; 0, 0, 0, 1;];
T_ges = t1;
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
