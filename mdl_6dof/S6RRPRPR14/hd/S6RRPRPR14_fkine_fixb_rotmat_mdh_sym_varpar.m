% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2018-11-23 17:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRPR14_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:09:47
% EndTime: 2018-11-23 17:09:47
% DurationCPUTime: 0.17s
% Computational Cost: add. (488->74), mult. (533->69), div. (0->0), fcn. (592->14), ass. (0->51)
t69 = pkin(5) + pkin(9);
t38 = sin(qJ(1));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t63 = pkin(6) + qJ(2);
t54 = sin(t63) / 0.2e1;
t64 = pkin(6) - qJ(2);
t56 = sin(t64);
t47 = t54 - t56 / 0.2e1;
t15 = t38 * t41 + t42 * t47;
t68 = t15 * pkin(9);
t17 = -t38 * t47 + t42 * t41;
t67 = t17 * pkin(9);
t55 = cos(t64) / 0.2e1;
t57 = cos(t63);
t23 = t55 - t57 / 0.2e1;
t66 = t23 * pkin(9);
t33 = sin(pkin(6));
t28 = t38 * t33;
t65 = t42 * t33;
t62 = pkin(7) + 0;
t61 = pkin(8) * t65;
t60 = t38 * pkin(1) + 0;
t34 = cos(pkin(6));
t59 = t34 * pkin(8) + t62;
t58 = t42 * pkin(1) + pkin(8) * t28 + 0;
t53 = t55 + t57 / 0.2e1;
t37 = sin(qJ(2));
t14 = t38 * t37 - t42 * t53;
t52 = t15 * pkin(2) + t14 * qJ(3) + t60;
t22 = t54 + t56 / 0.2e1;
t51 = t23 * pkin(2) - t22 * qJ(3) + t59;
t16 = t42 * t37 + t38 * t53;
t50 = t17 * pkin(2) + t16 * qJ(3) + t58;
t49 = t34 * pkin(3) + t51;
t48 = pkin(3) * t28 + t50;
t46 = (-pkin(3) - pkin(8)) * t65 + t52;
t36 = sin(qJ(4));
t40 = cos(qJ(4));
t12 = t22 * t40 + t34 * t36;
t13 = -t22 * t36 + t34 * t40;
t45 = t13 * pkin(4) + t12 * qJ(5) + t49;
t3 = -t16 * t40 + t36 * t28;
t4 = t16 * t36 + t40 * t28;
t44 = t4 * pkin(4) + t3 * qJ(5) + t48;
t5 = t14 * t40 + t36 * t65;
t6 = -t14 * t36 + t40 * t65;
t43 = -t6 * pkin(4) - t5 * qJ(5) + t46;
t39 = cos(qJ(6));
t35 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t42, -t38, 0, 0; t38, t42, 0, 0; 0, 0, 1, t62; 0, 0, 0, 1; t17, -t16, t28, t58; t15, -t14, -t65, t60 - t61; t23, t22, t34, t59; 0, 0, 0, 1; t28, -t17, t16, t50; -t65, -t15, t14, t52 - t61; t34, -t23, -t22, t51; 0, 0, 0, 1; t4, -t3, t17, t48 + t67; -t6, t5, t15, t46 + t68; t13, -t12, t23, t49 + t66; 0, 0, 0, 1; t17, -t4, t3, t44 + t67; t15, t6, -t5, t43 + t68; t23, -t13, t12, t45 + t66; 0, 0, 0, 1; t17 * t39 + t3 * t35, -t17 * t35 + t3 * t39, t4, t4 * pkin(10) + t69 * t17 + t44; t15 * t39 - t5 * t35, -t15 * t35 - t5 * t39, -t6, -t6 * pkin(10) + t69 * t15 + t43; t12 * t35 + t23 * t39, t12 * t39 - t23 * t35, t13, t13 * pkin(10) + t69 * t23 + t45; 0, 0, 0, 1;];
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
