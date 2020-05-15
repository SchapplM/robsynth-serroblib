% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 17:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:00:56
% EndTime: 2018-11-23 17:00:56
% DurationCPUTime: 0.15s
% Computational Cost: add. (199->66), mult. (137->74), div. (0->0), fcn. (196->12), ass. (0->40)
t25 = sin(qJ(4));
t44 = t25 * pkin(4);
t28 = cos(qJ(4));
t10 = t28 * pkin(4) + pkin(3);
t21 = qJ(2) + pkin(10);
t15 = cos(t21);
t27 = sin(qJ(1));
t43 = t27 * t15;
t42 = t27 * t25;
t41 = t27 * t28;
t30 = cos(qJ(1));
t40 = t30 * t15;
t39 = t30 * t25;
t38 = t30 * t28;
t23 = -qJ(5) - pkin(8);
t22 = pkin(6) + 0;
t20 = qJ(4) + pkin(11);
t29 = cos(qJ(2));
t11 = t29 * pkin(2) + pkin(1);
t37 = t30 * t11 + 0;
t24 = -qJ(3) - pkin(7);
t36 = t27 * t11 + t30 * t24 + 0;
t26 = sin(qJ(2));
t35 = t26 * pkin(2) + t22;
t13 = sin(t21);
t34 = pkin(3) * t15 + pkin(8) * t13;
t14 = cos(t20);
t1 = pkin(5) * t14 + t10;
t19 = -pkin(9) + t23;
t33 = t1 * t15 - t13 * t19;
t32 = t10 * t15 - t13 * t23;
t31 = -t27 * t24 + t37;
t16 = qJ(6) + t20;
t12 = sin(t20);
t8 = cos(t16);
t7 = sin(t16);
t6 = t30 * t13;
t5 = t27 * t13;
t2 = pkin(5) * t12 + t44;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t27, 0, 0; t27, t30, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t30 * t29, -t30 * t26, t27, t30 * pkin(1) + t27 * pkin(7) + 0; t27 * t29, -t27 * t26, -t30, t27 * pkin(1) - t30 * pkin(7) + 0; t26, t29, 0, t22; 0, 0, 0, 1; t40, -t6, t27, t31; t43, -t5, -t30, t36; t13, t15, 0, t35; 0, 0, 0, 1; t15 * t38 + t42, -t15 * t39 + t41, t6, t34 * t30 + t31; t15 * t41 - t39, -t15 * t42 - t38, t5, t34 * t27 + t36; t13 * t28, -t13 * t25, -t15, t13 * pkin(3) - t15 * pkin(8) + t35; 0, 0, 0, 1; t27 * t12 + t14 * t40, -t12 * t40 + t27 * t14, t6, t32 * t30 + (-t24 + t44) * t27 + t37; -t30 * t12 + t14 * t43, -t12 * t43 - t30 * t14, t5, -pkin(4) * t39 + t32 * t27 + t36; t13 * t14, -t13 * t12, -t15, t13 * t10 + t15 * t23 + t35; 0, 0, 0, 1; t27 * t7 + t8 * t40, t27 * t8 - t7 * t40, t6, t33 * t30 + (t2 - t24) * t27 + t37; -t30 * t7 + t8 * t43, -t30 * t8 - t7 * t43, t5, -t30 * t2 + t33 * t27 + t36; t13 * t8, -t13 * t7, -t15, t13 * t1 + t15 * t19 + t35; 0, 0, 0, 1;];
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
