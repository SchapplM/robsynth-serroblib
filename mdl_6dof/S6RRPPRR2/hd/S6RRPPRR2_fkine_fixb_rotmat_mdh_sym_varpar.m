% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2018-11-23 16:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPPRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:49:03
% EndTime: 2018-11-23 16:49:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (199->66), mult. (137->74), div. (0->0), fcn. (196->12), ass. (0->40)
t23 = sin(pkin(11));
t44 = t23 * pkin(4);
t24 = cos(pkin(11));
t9 = t24 * pkin(4) + pkin(3);
t21 = qJ(2) + pkin(10);
t15 = cos(t21);
t28 = sin(qJ(1));
t43 = t28 * t15;
t42 = t28 * t23;
t41 = t28 * t24;
t30 = cos(qJ(1));
t40 = t30 * t15;
t39 = t30 * t23;
t38 = t30 * t24;
t26 = -pkin(8) - qJ(4);
t22 = pkin(6) + 0;
t29 = cos(qJ(2));
t11 = t29 * pkin(2) + pkin(1);
t37 = t30 * t11 + 0;
t20 = pkin(11) + qJ(5);
t27 = sin(qJ(2));
t36 = t27 * pkin(2) + t22;
t25 = -qJ(3) - pkin(7);
t35 = t28 * t11 + t30 * t25 + 0;
t14 = cos(t20);
t1 = pkin(5) * t14 + t9;
t13 = sin(t21);
t19 = -pkin(9) + t26;
t34 = t1 * t15 - t13 * t19;
t33 = -t13 * t26 + t15 * t9;
t32 = pkin(3) * t15 + qJ(4) * t13;
t31 = -t28 * t25 + t37;
t16 = qJ(6) + t20;
t12 = sin(t20);
t8 = cos(t16);
t7 = sin(t16);
t6 = t30 * t13;
t5 = t28 * t13;
t2 = pkin(5) * t12 + t44;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t28, 0, 0; t28, t30, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t30 * t29, -t30 * t27, t28, t30 * pkin(1) + t28 * pkin(7) + 0; t28 * t29, -t28 * t27, -t30, t28 * pkin(1) - t30 * pkin(7) + 0; t27, t29, 0, t22; 0, 0, 0, 1; t40, -t6, t28, t31; t43, -t5, -t30, t35; t13, t15, 0, t36; 0, 0, 0, 1; t15 * t38 + t42, -t15 * t39 + t41, t6, t32 * t30 + t31; t15 * t41 - t39, -t15 * t42 - t38, t5, t32 * t28 + t35; t13 * t24, -t13 * t23, -t15, t13 * pkin(3) - t15 * qJ(4) + t36; 0, 0, 0, 1; t28 * t12 + t14 * t40, -t12 * t40 + t28 * t14, t6, t33 * t30 + (-t25 + t44) * t28 + t37; -t30 * t12 + t14 * t43, -t12 * t43 - t30 * t14, t5, -pkin(4) * t39 + t33 * t28 + t35; t13 * t14, -t13 * t12, -t15, t13 * t9 + t15 * t26 + t36; 0, 0, 0, 1; t28 * t7 + t8 * t40, t28 * t8 - t7 * t40, t6, t34 * t30 + (t2 - t25) * t28 + t37; -t30 * t7 + t8 * t43, -t30 * t8 - t7 * t43, t5, -t30 * t2 + t34 * t28 + t35; t13 * t8, -t13 * t7, -t15, t13 * t1 + t15 * t19 + t36; 0, 0, 0, 1;];
T_ges = t3;
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
