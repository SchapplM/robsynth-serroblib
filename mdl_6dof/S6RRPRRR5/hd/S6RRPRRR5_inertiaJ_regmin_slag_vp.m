% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t73 = sin(pkin(12));
t74 = sin(pkin(6));
t75 = cos(pkin(12));
t80 = sin(qJ(2));
t84 = cos(qJ(2));
t47 = (t73 * t84 + t75 * t80) * t74;
t76 = cos(pkin(6));
t79 = sin(qJ(4));
t83 = cos(qJ(4));
t36 = t47 * t79 - t76 * t83;
t132 = -0.2e1 * t36;
t131 = 0.2e1 * t36;
t82 = cos(qJ(5));
t67 = -t82 * pkin(5) - pkin(4);
t130 = 0.2e1 * t67;
t129 = -0.2e1 * t83;
t128 = 0.2e1 * t83;
t127 = pkin(10) + pkin(11);
t126 = pkin(1) * t80;
t125 = pkin(4) * t82;
t124 = t36 * pkin(5);
t77 = sin(qJ(6));
t123 = t77 * pkin(5);
t81 = cos(qJ(6));
t122 = t81 * pkin(5);
t37 = t47 * t83 + t76 * t79;
t108 = t74 * t84;
t109 = t74 * t80;
t46 = -t75 * t108 + t73 * t109;
t78 = sin(qJ(5));
t25 = t37 * t78 - t46 * t82;
t62 = t76 * t84 * pkin(1);
t92 = pkin(8) + qJ(3);
t38 = t76 * pkin(2) - t92 * t109 + t62;
t89 = t76 * t126;
t43 = t92 * t108 + t89;
t24 = t73 * t38 + t75 * t43;
t21 = t76 * pkin(9) + t24;
t57 = (-pkin(2) * t84 - pkin(1)) * t74;
t27 = t46 * pkin(3) - t47 * pkin(9) + t57;
t14 = t83 * t21 + t79 * t27;
t10 = t46 * pkin(10) + t14;
t23 = t75 * t38 - t73 * t43;
t20 = -t76 * pkin(3) - t23;
t12 = t36 * pkin(4) - t37 * pkin(10) + t20;
t7 = t82 * t10 + t78 * t12;
t5 = -t25 * pkin(11) + t7;
t121 = t81 * t5;
t120 = t83 * pkin(5);
t13 = -t79 * t21 + t83 * t27;
t9 = -t46 * pkin(4) - t13;
t119 = t9 * t78;
t118 = t9 * t82;
t26 = t37 * t82 + t46 * t78;
t16 = -t77 * t25 + t81 * t26;
t117 = t16 * t83;
t116 = t26 * t78;
t115 = t26 * t83;
t114 = t36 * t83;
t106 = t78 * t79;
t98 = t82 * t79;
t49 = -t77 * t106 + t81 * t98;
t113 = t49 * t36;
t56 = t77 * t82 + t81 * t78;
t112 = t56 * t83;
t64 = t73 * pkin(2) + pkin(9);
t111 = t64 * t78;
t68 = t74 ^ 2;
t110 = t68 * t84;
t107 = t78 * t36;
t105 = t78 * t82;
t104 = t78 * t83;
t103 = t79 * t36;
t102 = t79 * t46;
t101 = t79 * t64;
t65 = -t75 * pkin(2) - pkin(3);
t54 = -t83 * pkin(4) - t79 * pkin(10) + t65;
t93 = t83 * t64;
t86 = t82 * t93;
t31 = t86 + (-pkin(11) * t79 + t54) * t78;
t100 = t81 * t31;
t99 = t82 * t36;
t97 = t82 * t83;
t15 = t81 * t25 + t77 * t26;
t96 = t83 * t15;
t95 = t83 * t25;
t55 = t77 * t78 - t81 * t82;
t94 = t83 * t55;
t91 = 0.2e1 * t74 * t76;
t90 = t79 * t128;
t88 = t78 * t103;
t87 = t36 * t98;
t6 = -t78 * t10 + t82 * t12;
t4 = -t26 * pkin(11) + t124 + t6;
t1 = t81 * t4 - t77 * t5;
t50 = t82 * t54;
t30 = -pkin(11) * t98 + t50 + (-pkin(5) - t111) * t83;
t17 = t81 * t30 - t77 * t31;
t72 = t83 ^ 2;
t71 = t82 ^ 2;
t70 = t79 ^ 2;
t69 = t78 ^ 2;
t60 = t127 * t82;
t59 = t127 * t78;
t53 = pkin(8) * t108 + t89;
t52 = -pkin(8) * t109 + t62;
t51 = (pkin(5) * t78 + t64) * t79;
t48 = t56 * t79;
t44 = t83 * t46;
t42 = -t77 * t59 + t81 * t60;
t41 = -t81 * t59 - t77 * t60;
t35 = t36 ^ 2;
t34 = t78 * t54 + t86;
t33 = -t78 * t93 + t50;
t28 = t48 * t36;
t18 = t77 * t30 + t100;
t8 = t25 * pkin(5) + t9;
t2 = t77 * t4 + t121;
t3 = [1, 0, 0, t68 * t80 ^ 2, 0.2e1 * t80 * t110, t80 * t91, t84 * t91, t76 ^ 2, 0.2e1 * pkin(1) * t110 + 0.2e1 * t52 * t76, -0.2e1 * t68 * t126 - 0.2e1 * t53 * t76, -0.2e1 * t23 * t47 - 0.2e1 * t24 * t46, t23 ^ 2 + t24 ^ 2 + t57 ^ 2, t37 ^ 2, t37 * t132, 0.2e1 * t37 * t46, t46 * t132, t46 ^ 2, 0.2e1 * t13 * t46 + 0.2e1 * t20 * t36, -0.2e1 * t14 * t46 + 0.2e1 * t20 * t37, t26 ^ 2, -0.2e1 * t26 * t25, t26 * t131, t25 * t132, t35, 0.2e1 * t9 * t25 + 0.2e1 * t6 * t36, 0.2e1 * t9 * t26 - 0.2e1 * t7 * t36, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t131, t15 * t132, t35, 0.2e1 * t1 * t36 + 0.2e1 * t8 * t15, 0.2e1 * t8 * t16 - 0.2e1 * t2 * t36; 0, 0, 0, 0, 0, t109, t108, t76, t52, -t53 (-t46 * t73 - t47 * t75) * pkin(2) (t23 * t75 + t24 * t73) * pkin(2), t37 * t79, t37 * t83 - t103, t102, t44, 0, -t46 * t101 - t20 * t83 + t65 * t36, t20 * t79 + t65 * t37 - t46 * t93, t26 * t98 (-t25 * t82 - t116) * t79, t87 - t115, -t88 + t95, -t114, t33 * t36 - t6 * t83 + (t25 * t64 + t119) * t79, -t34 * t36 + t7 * t83 + (t26 * t64 + t118) * t79, t16 * t49, -t49 * t15 - t16 * t48, t113 - t117, -t28 + t96, -t114, -t1 * t83 + t51 * t15 + t17 * t36 + t8 * t48, t51 * t16 - t18 * t36 + t2 * t83 + t8 * t49; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t73 ^ 2 + t75 ^ 2) * pkin(2) ^ 2, t70, t90, 0, 0, 0, t65 * t129, 0.2e1 * t65 * t79, t71 * t70, -0.2e1 * t70 * t105, -0.2e1 * t79 * t97, t78 * t90, t72, 0.2e1 * t111 * t70 - 0.2e1 * t33 * t83, 0.2e1 * t70 * t64 * t82 + 0.2e1 * t34 * t83, t49 ^ 2, -0.2e1 * t49 * t48, t49 * t129, t48 * t128, t72, -0.2e1 * t17 * t83 + 0.2e1 * t51 * t48, 0.2e1 * t18 * t83 + 0.2e1 * t51 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, t44, -t102, 0, 0, 0, 0, 0, -t88 - t95, -t87 - t115, 0, 0, 0, 0, 0, -t28 - t96, -t113 - t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, t46, t13, -t14, t116, -t78 * t25 + t26 * t82, t107, t99, 0, -pkin(4) * t25 - pkin(10) * t107 - t118, -pkin(4) * t26 - pkin(10) * t99 + t119, t16 * t56, -t56 * t15 - t16 * t55, t56 * t36, -t55 * t36, 0, t67 * t15 + t41 * t36 + t8 * t55, t67 * t16 - t42 * t36 + t8 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t83, 0, -t101, -t93, t78 * t98 (-t69 + t71) * t79, -t104, -t97, 0, -t64 * t98 + (-pkin(4) * t79 + pkin(10) * t83) * t78, pkin(10) * t97 + (t111 - t125) * t79, t49 * t56, -t56 * t48 - t49 * t55, -t112, t94, 0, -t41 * t83 + t67 * t48 + t51 * t55, t42 * t83 + t67 * t49 + t51 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t79, 0, 0, 0, 0, 0, t97, -t104, 0, 0, 0, 0, 0, -t94, -t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t69, 0.2e1 * t105, 0, 0, 0, 0.2e1 * t125, -0.2e1 * pkin(4) * t78, t56 ^ 2, -0.2e1 * t56 * t55, 0, 0, 0, t55 * t130, t56 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t36, t6, -t7, 0, 0, t16, -t15, t36, t36 * t122 + t1, -t121 + (-t4 - t124) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -t106, -t83, t33, -t34, 0, 0, t49, -t48, -t83, -t81 * t120 + t17, -t100 + (-t30 + t120) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t98, 0, 0, 0, 0, 0, -t48, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t82, 0, -t78 * pkin(10), -t82 * pkin(10), 0, 0, t56, -t55, 0, t41, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t122, -0.2e1 * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t36, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t48, -t83, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t55, 0, t41, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t122, -t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
