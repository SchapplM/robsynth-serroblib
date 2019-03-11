% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t98 = -2 * pkin(3);
t51 = cos(pkin(7));
t54 = sin(qJ(5));
t58 = cos(qJ(5));
t49 = sin(pkin(7));
t59 = cos(qJ(3));
t87 = t49 * t59;
t27 = t54 * t51 + t58 * t87;
t97 = -0.2e1 * t27;
t96 = 0.2e1 * t49;
t57 = cos(qJ(6));
t95 = 0.2e1 * t57;
t94 = 2 * qJ(4);
t61 = -pkin(3) - pkin(10);
t55 = sin(qJ(3));
t93 = pkin(2) * t55;
t92 = pkin(2) * t59;
t41 = t49 * t55;
t36 = pkin(9) * t41;
t66 = -pkin(3) - t92;
t13 = pkin(4) * t41 + t36 + (-pkin(10) + t66) * t51;
t65 = -qJ(4) * t55 - pkin(2);
t20 = (t61 * t59 + t65) * t49;
t7 = t58 * t13 - t54 * t20;
t5 = -pkin(5) * t41 - t7;
t53 = sin(qJ(6));
t91 = t5 * t53;
t90 = t5 * t57;
t28 = t58 * t51 - t54 * t87;
t18 = t57 * t28 + t53 * t41;
t89 = t18 * t53;
t44 = t49 ^ 2;
t88 = t44 * t59;
t50 = sin(pkin(6));
t56 = sin(qJ(2));
t86 = t50 * t56;
t60 = cos(qJ(2));
t85 = t50 * t60;
t84 = t51 * t59;
t83 = t51 * t60;
t82 = t53 * t27;
t81 = t53 * t54;
t80 = t53 * t57;
t79 = t53 * t58;
t78 = t54 * t61;
t77 = t57 * t27;
t76 = t57 * t54;
t42 = t57 * t58;
t75 = t57 * t61;
t74 = t58 * t18;
t73 = t58 * t27;
t72 = t58 * t54;
t71 = t58 * t61;
t30 = pkin(9) * t87 + t51 * t93;
t46 = t54 ^ 2;
t48 = t58 ^ 2;
t70 = -t46 - t48;
t69 = 0.2e1 * t41;
t68 = -0.2e1 * t72;
t67 = t54 * t41;
t35 = t58 * t41;
t43 = t51 * qJ(4);
t23 = -t43 - t30;
t19 = pkin(4) * t87 - t23;
t64 = -pkin(5) * t58 - pkin(11) * t54;
t8 = t54 * t13 + t58 * t20;
t52 = cos(pkin(6));
t14 = -t50 * t59 * t83 - t52 * t87 + t55 * t86;
t26 = -t49 * t85 + t52 * t51;
t63 = t14 * t51 + t26 * t87;
t15 = t52 * t41 - (-t55 * t83 - t56 * t59) * t50;
t62 = -t15 * t51 + t26 * t41;
t47 = t57 ^ 2;
t45 = t53 ^ 2;
t40 = t44 * t55 ^ 2;
t32 = t54 * pkin(5) - t58 * pkin(11) + qJ(4);
t29 = pkin(2) * t84 - t36;
t25 = t66 * t51 + t36;
t24 = (-pkin(3) * t59 + t65) * t49;
t22 = t53 * t32 + t54 * t75;
t21 = t57 * t32 - t53 * t78;
t17 = t53 * t28 - t57 * t41;
t11 = t14 * t54 + t26 * t58;
t10 = -t14 * t58 + t26 * t54;
t9 = t27 * pkin(5) - t28 * pkin(11) + t19;
t6 = pkin(11) * t41 + t8;
t4 = t11 * t57 + t15 * t53;
t3 = -t11 * t53 + t15 * t57;
t2 = t53 * t9 + t57 * t6;
t1 = -t53 * t6 + t57 * t9;
t12 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t15 ^ 2 + t26 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t85, -t86, 0, 0, 0, 0, 0, -t63, t62 (t14 * t55 + t15 * t59) * t49, t63, -t62, t14 * t25 - t15 * t23 + t26 * t24, 0, 0, 0, 0, 0, -t10 * t41 + t15 * t27, -t11 * t41 + t15 * t28, 0, 0, 0, 0, 0, t10 * t17 + t3 * t27, t10 * t18 - t4 * t27; 0, 1, 0, 0, t40, 0.2e1 * t55 * t88, t51 * t69, t84 * t96, t51 ^ 2, 0.2e1 * pkin(2) * t88 + 0.2e1 * t29 * t51, -0.2e1 * t30 * t51 - 0.2e1 * t44 * t93 (-t23 * t59 + t25 * t55) * t96, 0.2e1 * t24 * t87 + 0.2e1 * t25 * t51, -0.2e1 * t23 * t51 - 0.2e1 * t24 * t41, t23 ^ 2 + t24 ^ 2 + t25 ^ 2, t28 ^ 2, t28 * t97, t28 * t69, t41 * t97, t40, 0.2e1 * t19 * t27 + 0.2e1 * t7 * t41, 0.2e1 * t19 * t28 - 0.2e1 * t8 * t41, t18 ^ 2, -0.2e1 * t18 * t17, 0.2e1 * t18 * t27, t17 * t97, t27 ^ 2, 0.2e1 * t1 * t27 + 0.2e1 * t5 * t17, 0.2e1 * t5 * t18 - 0.2e1 * t2 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, t14, t15, -t14 * pkin(3) + t15 * qJ(4), 0, 0, 0, 0, 0, t15 * t54, t15 * t58, 0, 0, 0, 0, 0, t10 * t79 + t3 * t54, t10 * t42 - t4 * t54; 0, 0, 0, 0, 0, 0, t41, t87, t51, t29, -t30 (-pkin(3) * t55 + qJ(4) * t59) * t49, t36 + (t98 - t92) * t51, 0.2e1 * t43 + t30, -t25 * pkin(3) - t23 * qJ(4), t28 * t58, -t28 * t54 - t73, t35, -t67, 0, qJ(4) * t27 + t19 * t54 + t61 * t35, qJ(4) * t28 + t19 * t58 - t61 * t67, t57 * t74 (-t17 * t57 - t89) * t58, t18 * t54 + t57 * t73, -t17 * t54 - t53 * t73, t27 * t54, t1 * t54 + t21 * t27 + (-t17 * t61 + t91) * t58, -t2 * t54 - t22 * t27 + (-t18 * t61 + t90) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t98, t94, pkin(3) ^ 2 + qJ(4) ^ 2, t48, t68, 0, 0, 0, t54 * t94, t58 * t94, t47 * t48, -0.2e1 * t48 * t80, t72 * t95, t53 * t68, t46, -0.2e1 * t48 * t61 * t53 + 0.2e1 * t21 * t54, -0.2e1 * t22 * t54 - 0.2e1 * t48 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t51, 0, t25, 0, 0, 0, 0, 0, t35, -t67, 0, 0, 0, 0, 0, -t58 * t17 - t27 * t81, -t27 * t76 - t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t53, t70 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t10 * t57, t10 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, t41, t7, -t8, t89, -t53 * t17 + t18 * t57, t82, t77, 0, -pkin(5) * t17 - pkin(11) * t82 - t90, -pkin(5) * t18 - pkin(11) * t77 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t54, 0, t71, -t78, t53 * t42 (-t45 + t47) * t58, t81, t76, 0, t64 * t53 + t57 * t71, -t53 * t71 + t64 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t54, 0, 0, 0, 0, 0, t42, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t45, 0.2e1 * t80, 0, 0, 0, pkin(5) * t95, -0.2e1 * pkin(5) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, t27, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t79, t54, t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t57, 0, -t53 * pkin(11), -t57 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t12;
