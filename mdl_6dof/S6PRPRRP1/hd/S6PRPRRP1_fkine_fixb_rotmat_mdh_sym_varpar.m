% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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

function [T_c_mdh, Tc_stack] = S6PRPRRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:59:29
% EndTime: 2018-11-23 14:59:29
% DurationCPUTime: 0.23s
% Computational Cost: add. (705->81), mult. (534->91), div. (0->0), fcn. (579->20), ass. (0->67)
t50 = pkin(6) - qJ(2);
t38 = cos(t50) / 0.2e1;
t49 = pkin(6) + qJ(2);
t46 = cos(t49);
t84 = t38 - t46 / 0.2e1;
t37 = sin(t49) / 0.2e1;
t44 = sin(t50);
t26 = t37 - t44 / 0.2e1;
t51 = sin(pkin(10));
t52 = sin(pkin(6));
t35 = t51 * t52;
t53 = cos(pkin(10));
t81 = t53 * t52;
t48 = qJ(2) + pkin(11);
t80 = qJ(1) + 0;
t57 = sin(qJ(5));
t79 = pkin(5) * t57 + pkin(8);
t56 = pkin(7) + qJ(3);
t20 = pkin(2) * t26 - t52 * t56;
t62 = cos(qJ(2));
t41 = pkin(2) * t62 + pkin(1);
t78 = t20 * t53 + t41 * t51 + 0;
t77 = pkin(6) - t48;
t76 = pkin(6) + t48;
t68 = sin(t76) / 0.2e1;
t72 = sin(t77);
t23 = t68 - t72 / 0.2e1;
t45 = cos(t48);
t14 = t23 * t53 + t45 * t51;
t75 = pkin(3) * t14 + t78;
t74 = -t51 * t20 + t41 * t53 + 0;
t73 = cos(t76);
t54 = cos(pkin(6));
t71 = pkin(2) * t84 + t54 * t56 + t80;
t16 = -t23 * t51 + t45 * t53;
t70 = pkin(3) * t16 + t74;
t69 = cos(t77) / 0.2e1;
t25 = t69 - t73 / 0.2e1;
t67 = pkin(3) * t25 + t71;
t42 = sin(t48);
t63 = t73 / 0.2e1 + t69;
t13 = t42 * t51 - t53 * t63;
t66 = pkin(8) * t13 + t75;
t15 = t42 * t53 + t51 * t63;
t65 = pkin(8) * t15 + t70;
t24 = t72 / 0.2e1 + t68;
t64 = -pkin(8) * t24 + t67;
t61 = cos(qJ(4));
t60 = cos(qJ(5));
t59 = sin(qJ(2));
t58 = sin(qJ(4));
t55 = -qJ(6) - pkin(9);
t40 = pkin(5) * t60 + pkin(4);
t27 = t38 + t46 / 0.2e1;
t18 = t25 * t61 + t54 * t58;
t17 = t25 * t58 - t54 * t61;
t10 = t16 * t61 + t35 * t58;
t9 = t16 * t58 - t35 * t61;
t8 = t14 * t61 - t58 * t81;
t7 = t14 * t58 + t61 * t81;
t6 = t18 * t60 - t24 * t57;
t5 = -t18 * t57 - t24 * t60;
t4 = t10 * t60 + t15 * t57;
t3 = -t10 * t57 + t15 * t60;
t2 = t13 * t57 + t60 * t8;
t1 = t13 * t60 - t57 * t8;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t53, -t51, 0, 0; t51, t53, 0, 0; 0, 0, 1, t80; 0, 0, 0, 1; -t26 * t51 + t53 * t62, -t27 * t51 - t53 * t59, t35, pkin(1) * t53 + pkin(7) * t35 + 0; t26 * t53 + t51 * t62, t27 * t53 - t51 * t59, -t81, pkin(1) * t51 - pkin(7) * t81 + 0; t84, t37 + t44 / 0.2e1, t54, pkin(7) * t54 + t80; 0, 0, 0, 1; t16, -t15, t35, t74; t14, -t13, -t81, t78; t25, t24, t54, t71; 0, 0, 0, 1; t10, -t9, t15, t65; t8, -t7, t13, t66; t18, -t17, -t24, t64; 0, 0, 0, 1; t4, t3, t9, pkin(4) * t10 + pkin(9) * t9 + t65; t2, t1, t7, pkin(4) * t8 + pkin(9) * t7 + t66; t6, t5, t17, pkin(4) * t18 + pkin(9) * t17 + t64; 0, 0, 0, 1; t4, t3, t9, t10 * t40 + t15 * t79 - t9 * t55 + t70; t2, t1, t7, t13 * t79 + t8 * t40 - t7 * t55 + t75; t6, t5, t17, -t17 * t55 + t18 * t40 - t24 * t79 + t67; 0, 0, 0, 1;];
T_ges = t11;
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
