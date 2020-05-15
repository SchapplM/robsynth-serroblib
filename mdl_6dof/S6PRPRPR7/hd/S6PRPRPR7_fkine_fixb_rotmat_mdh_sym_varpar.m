% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2018-11-23 14:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:58:51
% EndTime: 2018-11-23 14:58:51
% DurationCPUTime: 0.20s
% Computational Cost: add. (488->74), mult. (533->69), div. (0->0), fcn. (592->14), ass. (0->51)
t69 = pkin(5) + pkin(8);
t33 = sin(pkin(10));
t35 = cos(pkin(10));
t42 = cos(qJ(2));
t63 = pkin(6) + qJ(2);
t54 = sin(t63) / 0.2e1;
t64 = pkin(6) - qJ(2);
t56 = sin(t64);
t47 = t54 - t56 / 0.2e1;
t13 = t33 * t42 + t35 * t47;
t68 = t13 * pkin(8);
t15 = -t33 * t47 + t35 * t42;
t67 = t15 * pkin(8);
t55 = cos(t64) / 0.2e1;
t57 = cos(t63);
t23 = t55 - t57 / 0.2e1;
t66 = t23 * pkin(8);
t34 = sin(pkin(6));
t26 = t33 * t34;
t65 = t35 * t34;
t62 = pkin(7) * t65;
t61 = t33 * pkin(1) + 0;
t60 = qJ(1) + 0;
t59 = t35 * pkin(1) + pkin(7) * t26 + 0;
t36 = cos(pkin(6));
t58 = t36 * pkin(7) + t60;
t53 = t55 + t57 / 0.2e1;
t39 = sin(qJ(2));
t12 = t33 * t39 - t35 * t53;
t52 = t13 * pkin(2) + t12 * qJ(3) + t61;
t14 = t33 * t53 + t35 * t39;
t51 = t15 * pkin(2) + t14 * qJ(3) + t59;
t22 = t54 + t56 / 0.2e1;
t50 = t23 * pkin(2) - t22 * qJ(3) + t58;
t49 = pkin(3) * t26 + t51;
t48 = t36 * pkin(3) + t50;
t46 = (-pkin(3) - pkin(7)) * t65 + t52;
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t3 = -t14 * t41 + t38 * t26;
t4 = t14 * t38 + t41 * t26;
t45 = t4 * pkin(4) + t3 * qJ(5) + t49;
t16 = t22 * t41 + t36 * t38;
t17 = -t22 * t38 + t36 * t41;
t44 = t17 * pkin(4) + t16 * qJ(5) + t48;
t5 = t12 * t41 + t38 * t65;
t6 = -t12 * t38 + t41 * t65;
t43 = -t6 * pkin(4) - t5 * qJ(5) + t46;
t40 = cos(qJ(6));
t37 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t33, 0, 0; t33, t35, 0, 0; 0, 0, 1, t60; 0, 0, 0, 1; t15, -t14, t26, t59; t13, -t12, -t65, t61 - t62; t23, t22, t36, t58; 0, 0, 0, 1; t26, -t15, t14, t51; -t65, -t13, t12, t52 - t62; t36, -t23, -t22, t50; 0, 0, 0, 1; t4, -t3, t15, t49 + t67; -t6, t5, t13, t46 + t68; t17, -t16, t23, t48 + t66; 0, 0, 0, 1; t15, -t4, t3, t45 + t67; t13, t6, -t5, t43 + t68; t23, -t17, t16, t44 + t66; 0, 0, 0, 1; t15 * t40 + t3 * t37, -t15 * t37 + t3 * t40, t4, t4 * pkin(9) + t69 * t15 + t45; t13 * t40 - t5 * t37, -t13 * t37 - t5 * t40, -t6, -t6 * pkin(9) + t69 * t13 + t43; t16 * t37 + t23 * t40, t16 * t40 - t23 * t37, t17, t17 * pkin(9) + t69 * t23 + t44; 0, 0, 0, 1;];
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
