% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t51 = cos(pkin(6));
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t50 = sin(pkin(6));
t58 = cos(qJ(2));
t83 = t50 * t58;
t23 = -t51 * t54 - t57 * t83;
t55 = sin(qJ(2));
t36 = t50 * t55;
t53 = sin(qJ(6));
t56 = cos(qJ(6));
t13 = t23 * t56 + t53 * t36;
t98 = -0.2e1 * t13;
t21 = -t51 * t57 + t54 * t83;
t97 = 0.2e1 * t21;
t96 = 0.2e1 * t50;
t95 = -0.2e1 * t54;
t94 = 0.2e1 * t57;
t59 = (-pkin(2) - pkin(3));
t93 = 2 * t59;
t92 = pkin(1) * t55;
t91 = pkin(5) * t56;
t19 = -t50 * pkin(1) - pkin(2) * t83 - qJ(3) * t36;
t15 = pkin(3) * t83 - t19;
t8 = (pkin(4) * t55 + pkin(9) * t58) * t50 + t15;
t45 = -pkin(9) + t59;
t70 = qJ(4) * t50;
t82 = t51 * t58;
t73 = -pkin(1) * t82 + pkin(8) * t36;
t63 = -t55 * t70 + t73;
t9 = t45 * t51 + t63;
t5 = -t54 * t9 + t57 * t8;
t3 = -pkin(5) * t36 - t5;
t90 = t3 * t53;
t89 = t3 * t56;
t88 = t51 * pkin(2);
t87 = t13 * t53;
t86 = t13 * t57;
t44 = t50 ^ 2;
t85 = t44 * t58;
t84 = t45 * t53;
t81 = t53 * t21;
t80 = t53 * t56;
t79 = t53 * t57;
t78 = t54 * t21;
t77 = t56 * t21;
t76 = t56 * t54;
t39 = t56 * t57;
t12 = t23 * t53 - t56 * t36;
t75 = t57 * t12;
t74 = t57 * t45;
t52 = qJ(3) + pkin(4);
t25 = pkin(8) * t83 + t51 * t92;
t47 = t54 ^ 2;
t49 = t57 ^ 2;
t72 = -t47 - t49;
t71 = qJ(3) * t58;
t69 = 0.2e1 * t36;
t68 = t54 * t94;
t67 = t54 * t36;
t66 = t57 * t36;
t65 = t53 * t78;
t64 = t21 * t76;
t6 = t54 * t8 + t57 * t9;
t62 = -t58 * t70 + t25;
t41 = t51 * qJ(3);
t14 = -t41 - t62;
t10 = t51 * pkin(4) - t14;
t61 = qJ(3) ^ 2;
t60 = 0.2e1 * qJ(3);
t48 = t56 ^ 2;
t46 = t53 ^ 2;
t40 = 0.2e1 * t41;
t38 = t53 * t54;
t35 = t44 * t55 ^ 2;
t26 = t57 * pkin(5) + t54 * pkin(10) + t52;
t20 = t73 - t88;
t18 = t41 + t25;
t17 = t53 * t26 + t56 * t74;
t16 = t56 * t26 - t53 * t74;
t11 = t59 * t51 + t63;
t7 = -t21 * pkin(5) - t23 * pkin(10) + t10;
t4 = pkin(10) * t36 + t6;
t2 = t56 * t4 + t53 * t7;
t1 = -t53 * t4 + t56 * t7;
t22 = [1, 0, 0, t35, 0.2e1 * t55 * t85, t51 * t69, t82 * t96, t51 ^ 2, 0.2e1 * pkin(1) * t85 - 0.2e1 * t51 * t73, -0.2e1 * t25 * t51 - 0.2e1 * t44 * t92, -0.2e1 * t19 * t83 - 0.2e1 * t20 * t51 (t18 * t58 + t20 * t55) * t96, 0.2e1 * t18 * t51 - 0.2e1 * t19 * t36, t18 ^ 2 + t19 ^ 2 + t20 ^ 2, -0.2e1 * t14 * t51 + 0.2e1 * t15 * t36, 0.2e1 * t11 * t51 - 0.2e1 * t15 * t83 (-t11 * t55 + t14 * t58) * t96, t11 ^ 2 + t14 ^ 2 + t15 ^ 2, t23 ^ 2, t23 * t97, t23 * t69, t21 * t69, t35, -0.2e1 * t10 * t21 + 0.2e1 * t5 * t36, 0.2e1 * t10 * t23 - 0.2e1 * t6 * t36, t13 ^ 2, t12 * t98, t21 * t98, t12 * t97, t21 ^ 2, -0.2e1 * t1 * t21 + 0.2e1 * t3 * t12, 0.2e1 * t3 * t13 + 0.2e1 * t2 * t21; 0, 0, 0, 0, 0, t36, t83, t51, -t73, -t25, -t73 + 0.2e1 * t88 (-pkin(2) * t55 + t71) * t50, t40 + t25, -t20 * pkin(2) + t18 * qJ(3), t40 + t62, t51 * t93 + t63 (-t55 * t59 - t71) * t50, -t14 * qJ(3) + t11 * t59, -t23 * t54, -t23 * t57 - t78, -t67, -t66, 0, t10 * t57 - t52 * t21 - t45 * t67, -t10 * t54 + t52 * t23 - t45 * t66, -t13 * t76 (t12 * t56 + t87) * t54, t64 + t86, -t65 - t75, -t21 * t57, t1 * t57 - t16 * t21 + (t12 * t45 - t90) * t54, t17 * t21 - t2 * t57 + (t13 * t45 - t89) * t54; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, t60, pkin(2) ^ 2 + t61, t60, t93, 0 (t59 ^ 2) + t61, t47, t68, 0, 0, 0, t52 * t94, t52 * t95, t48 * t47, -0.2e1 * t47 * t80, t39 * t95, t53 * t68, t49, 0.2e1 * t16 * t57 - 0.2e1 * t47 * t84, -0.2e1 * t47 * t45 * t56 - 0.2e1 * t17 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t36, 0, t20, 0, t51, -t36, t11, 0, 0, 0, 0, 0, -t67, -t66, 0, 0, 0, 0, 0, t54 * t12 + t21 * t79, t54 * t13 + t21 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 1, 0, t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t53, t72 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t83, 0, t15, 0, 0, 0, 0, 0, t66, -t67, 0, 0, 0, 0, 0, t65 - t75, t64 - t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t21, t36, t5, -t6, t87, -t53 * t12 + t13 * t56, -t81, -t77, 0, -pkin(5) * t12 + pkin(10) * t81 - t89, -pkin(5) * t13 + pkin(10) * t77 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t57, 0, -t54 * t45, -t74, -t53 * t76 (t46 - t48) * t54, t79, t39, 0, -t45 * t76 + (pkin(5) * t54 - pkin(10) * t57) * t53, -pkin(10) * t39 + (t84 + t91) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t57, 0, 0, 0, 0, 0, -t76, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t54, 0, 0, 0, 0, 0, t39, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t46, 0.2e1 * t80, 0, 0, 0, 0.2e1 * t91, -0.2e1 * pkin(5) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, -t21, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t38, t57, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t56, 0, -t53 * pkin(10), -t56 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t22;
