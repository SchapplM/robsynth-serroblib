% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR14_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t51 = sin(pkin(5));
t50 = sin(pkin(6));
t65 = cos(pkin(5));
t61 = t65 * t50;
t52 = cos(pkin(11));
t53 = cos(pkin(6));
t77 = t52 * t53;
t93 = t51 * t77 + t61;
t49 = sin(pkin(11));
t62 = pkin(1) * t65;
t66 = qJ(2) * t51;
t34 = t49 * t62 + t52 * t66;
t20 = t93 * pkin(8) + t34;
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t42 = t52 * t62;
t81 = t49 * t51;
t23 = t65 * pkin(2) + t42 + (-pkin(8) * t53 - qJ(2)) * t81;
t27 = (-pkin(8) * t49 * t50 - pkin(2) * t52 - pkin(1)) * t51;
t60 = t23 * t53 + t27 * t50;
t11 = -t56 * t20 + t59 * t60;
t22 = t56 * t61 + (t49 * t59 + t56 * t77) * t51;
t78 = t51 * t52;
t31 = t50 * t78 - t53 * t65;
t55 = sin(qJ(4));
t58 = cos(qJ(4));
t16 = t22 * t55 + t31 * t58;
t92 = -0.2e1 * t16;
t21 = t56 * t81 - t93 * t59;
t91 = 0.2e1 * t21;
t90 = -0.2e1 * t22;
t89 = -0.2e1 * t55;
t88 = 0.2e1 * t58;
t45 = t51 ^ 2;
t87 = pkin(1) * t45;
t57 = cos(qJ(5));
t86 = pkin(4) * t57;
t54 = sin(qJ(5));
t85 = pkin(9) * t54;
t12 = t59 * t20 + t56 * t60;
t10 = -t31 * pkin(9) + t12;
t15 = -t50 * t23 + t53 * t27;
t8 = t21 * pkin(3) - t22 * pkin(9) + t15;
t5 = -t55 * t10 + t58 * t8;
t3 = -t21 * pkin(4) - t5;
t84 = t3 * t54;
t83 = t3 * t57;
t17 = t22 * t58 - t31 * t55;
t14 = t17 * t57 + t21 * t54;
t82 = t14 * t54;
t80 = t50 * t56;
t79 = t50 * t59;
t76 = t54 * t16;
t75 = t54 * t55;
t74 = t54 * t57;
t73 = t54 * t58;
t72 = t55 * t16;
t71 = t55 * t21;
t70 = t57 * t16;
t69 = t57 * t55;
t68 = t57 * t58;
t67 = t58 * t21;
t64 = t55 * t88;
t6 = t58 * t10 + t55 * t8;
t9 = t31 * pkin(3) - t11;
t48 = t57 ^ 2;
t47 = t55 ^ 2;
t46 = t54 ^ 2;
t38 = -t58 * pkin(4) - t55 * pkin(10) - pkin(3);
t36 = t55 * t53 + t58 * t80;
t35 = -t58 * t53 + t55 * t80;
t33 = -t49 * t66 + t42;
t29 = pkin(9) * t68 + t54 * t38;
t28 = -pkin(9) * t73 + t57 * t38;
t25 = t57 * t36 - t54 * t79;
t24 = -t54 * t36 - t57 * t79;
t13 = t17 * t54 - t21 * t57;
t7 = t16 * pkin(4) - t17 * pkin(10) + t9;
t4 = t21 * pkin(10) + t6;
t2 = t57 * t4 + t54 * t7;
t1 = -t54 * t4 + t57 * t7;
t18 = [1, 0, 0, 0.2e1 * t33 * t65 + 0.2e1 * t52 * t87, -0.2e1 * t34 * t65 - 0.2e1 * t49 * t87, 0.2e1 * (-t33 * t49 + t34 * t52) * t51, t45 * pkin(1) ^ 2 + t33 ^ 2 + t34 ^ 2, t22 ^ 2, t21 * t90, t31 * t90, t31 * t91, t31 ^ 2, -0.2e1 * t11 * t31 + 0.2e1 * t15 * t21, 0.2e1 * t12 * t31 + 0.2e1 * t15 * t22, t17 ^ 2, t17 * t92, t17 * t91, t21 * t92, t21 ^ 2, 0.2e1 * t9 * t16 + 0.2e1 * t5 * t21, 0.2e1 * t9 * t17 - 0.2e1 * t6 * t21, t14 ^ 2, -0.2e1 * t14 * t13, 0.2e1 * t14 * t16, t13 * t92, t16 ^ 2, 0.2e1 * t1 * t16 + 0.2e1 * t3 * t13, 0.2e1 * t3 * t14 - 0.2e1 * t2 * t16; 0, 0, 0, -t78, t81, 0, -t51 * pkin(1), 0, 0, 0, 0, 0, t53 * t21 - t31 * t79, t53 * t22 + t31 * t80, 0, 0, 0, 0, 0, -t16 * t79 - t35 * t21, -t17 * t79 - t36 * t21, 0, 0, 0, 0, 0, t35 * t13 + t24 * t16, t35 * t14 - t25 * t16; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, -t31, t11, -t12, t17 * t55, t17 * t58 - t72, t71, t67, 0, -pkin(3) * t16 - pkin(9) * t71 - t9 * t58, -pkin(3) * t17 - pkin(9) * t67 + t9 * t55, t14 * t69, (-t13 * t57 - t82) * t55, -t14 * t58 + t16 * t69, t13 * t58 - t54 * t72, -t16 * t58, -t1 * t58 + t28 * t16 + (pkin(9) * t13 + t84) * t55, -t29 * t16 + t2 * t58 + (pkin(9) * t14 + t83) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t80, 0, 0, 0, 0, 0, t58 * t79, -t55 * t79, 0, 0, 0, 0, 0, -t24 * t58 + t35 * t75, t25 * t58 + t35 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t47, t64, 0, 0, 0, pkin(3) * t88, pkin(3) * t89, t48 * t47, -0.2e1 * t47 * t74, t68 * t89, t54 * t64, t58 ^ 2, -0.2e1 * t28 * t58 + 0.2e1 * t47 * t85, 0.2e1 * t47 * pkin(9) * t57 + 0.2e1 * t29 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t21, t5, -t6, t82, -t54 * t13 + t14 * t57, t76, t70, 0, -pkin(4) * t13 - pkin(10) * t76 - t83, -pkin(4) * t14 - pkin(10) * t70 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, 0, 0, 0, 0, -t35 * t57, t35 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t58, 0, -t55 * pkin(9), -t58 * pkin(9), t54 * t69, (-t46 + t48) * t55, -t73, -t68, 0, -pkin(9) * t69 + (-pkin(4) * t55 + pkin(10) * t58) * t54, pkin(10) * t68 + (t85 - t86) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t46, 0.2e1 * t74, 0, 0, 0, 0.2e1 * t86, -0.2e1 * pkin(4) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t75, -t58, t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t57, 0, -t54 * pkin(10), -t57 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t18;
