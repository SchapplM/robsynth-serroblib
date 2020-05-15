% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2018-11-23 17:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:05:04
% EndTime: 2018-11-23 17:05:04
% DurationCPUTime: 0.17s
% Computational Cost: add. (204->64), mult. (230->71), div. (0->0), fcn. (314->10), ass. (0->39)
t24 = pkin(10) + qJ(4);
t19 = sin(t24);
t20 = cos(t24);
t34 = cos(qJ(1));
t31 = sin(qJ(1));
t33 = cos(qJ(2));
t48 = t31 * t33;
t3 = t19 * t48 + t34 * t20;
t4 = -t34 * t19 + t20 * t48;
t53 = t4 * pkin(4) + t3 * qJ(5);
t47 = t34 * t33;
t5 = t19 * t47 - t31 * t20;
t6 = t31 * t19 + t20 * t47;
t52 = t6 * pkin(4) + t5 * qJ(5);
t30 = sin(qJ(2));
t50 = t30 * t19;
t12 = t30 * t20;
t26 = sin(pkin(10));
t49 = t31 * t26;
t17 = t31 * t30;
t18 = t34 * t30;
t25 = pkin(6) + 0;
t45 = t31 * pkin(1) + 0;
t28 = -pkin(8) - qJ(3);
t44 = t30 * (-pkin(9) - t28);
t43 = t34 * pkin(1) + t31 * pkin(7) + 0;
t27 = cos(pkin(10));
t15 = t27 * pkin(3) + pkin(2);
t42 = t30 * t15 + t33 * t28 + t25;
t41 = pkin(2) * t33 + qJ(3) * t30;
t40 = -t34 * pkin(7) + t45;
t39 = pkin(3) * t49 + t15 * t47 + t43;
t38 = pkin(4) * t12 + qJ(5) * t50 + t42;
t37 = t15 * t48 + (-pkin(3) * t26 - pkin(7)) * t34 + t45;
t36 = -t28 * t18 + t39;
t35 = -t28 * t17 + t37;
t32 = cos(qJ(6));
t29 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t31, 0, 0; t31, t34, 0, 0; 0, 0, 1, t25; 0, 0, 0, 1; t47, -t18, t31, t43; t48, -t17, -t34, t40; t30, t33, 0, t25; 0, 0, 0, 1; t27 * t47 + t49, -t26 * t47 + t31 * t27, t18, t41 * t34 + t43; -t34 * t26 + t27 * t48, -t26 * t48 - t34 * t27, t17, t41 * t31 + t40; t30 * t27, -t30 * t26, -t33, t30 * pkin(2) - t33 * qJ(3) + t25; 0, 0, 0, 1; t6, -t5, t18, t36; t4, -t3, t17, t35; t12, -t50, -t33, t42; 0, 0, 0, 1; t6, t18, t5, t36 + t52; t4, t17, t3, t35 + t53; t12, -t33, t50, t38; 0, 0, 0, 1; t5 * t29 + t6 * t32, -t6 * t29 + t5 * t32, -t18, t6 * pkin(5) + t34 * t44 + t39 + t52; t3 * t29 + t4 * t32, -t4 * t29 + t3 * t32, -t17, t4 * pkin(5) + t31 * t44 + t37 + t53; (t19 * t29 + t20 * t32) * t30 (t19 * t32 - t20 * t29) * t30, t33, pkin(5) * t12 + t33 * pkin(9) + t38; 0, 0, 0, 1;];
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
