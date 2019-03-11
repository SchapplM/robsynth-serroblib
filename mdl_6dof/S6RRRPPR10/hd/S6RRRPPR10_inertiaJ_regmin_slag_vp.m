% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t78 = sin(pkin(11));
t80 = cos(pkin(11));
t83 = sin(qJ(6));
t86 = cos(qJ(6));
t55 = -t83 * t78 + t86 * t80;
t79 = sin(pkin(6));
t85 = sin(qJ(2));
t112 = t79 * t85;
t81 = cos(pkin(6));
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t49 = t87 * t112 + t81 * t84;
t127 = -0.2e1 * t49;
t67 = t78 * pkin(5) + qJ(4);
t126 = 0.2e1 * t67;
t125 = -0.2e1 * t84;
t124 = 0.2e1 * t84;
t123 = 0.2e1 * t87;
t122 = 0.2e1 * qJ(4);
t121 = pkin(1) * t85;
t88 = cos(qJ(2));
t120 = pkin(1) * t88;
t106 = pkin(3) + qJ(5);
t119 = -pkin(10) - t106;
t118 = t49 * t80;
t41 = t49 * t84;
t54 = t86 * t78 + t83 * t80;
t117 = t54 * t49;
t116 = t54 * t84;
t73 = t79 ^ 2;
t115 = t73 * t88;
t114 = t78 * t49;
t113 = t78 * t87;
t111 = t79 * t88;
t110 = t80 * t87;
t109 = t81 * t85;
t102 = pkin(8) * t111;
t38 = t102 + (pkin(9) + t121) * t81;
t39 = (-pkin(2) * t88 - pkin(9) * t85 - pkin(1)) * t79;
t22 = -t84 * t38 + t87 * t39;
t66 = pkin(3) * t111;
t21 = -t22 + t66;
t12 = t49 * pkin(4) + qJ(5) * t111 + t21;
t48 = t84 * t112 - t81 * t87;
t65 = pkin(8) * t112;
t37 = t65 + (-pkin(2) - t120) * t81;
t91 = -t49 * qJ(4) + t37;
t14 = t106 * t48 + t91;
t6 = t78 * t12 + t80 * t14;
t23 = t87 * t38 + t84 * t39;
t98 = -t84 * qJ(4) - pkin(2);
t52 = -t106 * t87 + t98;
t69 = t84 * pkin(9);
t61 = t84 * pkin(4) + t69;
t28 = t80 * t52 + t78 * t61;
t70 = t87 * pkin(9);
t62 = t87 * pkin(4) + t70;
t64 = t78 ^ 2 + t80 ^ 2;
t76 = t84 ^ 2;
t105 = t87 ^ 2 + t76;
t104 = qJ(4) * t87;
t103 = 0.2e1 * t111;
t101 = t84 * t111;
t100 = t87 * t111;
t99 = qJ(4) * t111;
t5 = t80 * t12 - t78 * t14;
t97 = pkin(9) * t101;
t96 = pkin(9) * t100;
t95 = t5 * t80 + t6 * t78;
t94 = -pkin(3) * t84 + t104;
t20 = t99 - t23;
t93 = -t20 * t87 + t21 * t84;
t57 = t80 * t61;
t27 = -t78 * t52 + t57;
t15 = t27 * t80 + t28 * t78;
t92 = -t106 * t84 + t104;
t16 = -t48 * pkin(4) - t20;
t89 = qJ(4) ^ 2;
t60 = -t87 * pkin(3) + t98;
t59 = t119 * t80;
t58 = t119 * t78;
t53 = t64 * t106;
t51 = pkin(1) * t109 + t102;
t50 = t81 * t120 - t65;
t47 = pkin(5) * t110 + t62;
t46 = t49 ^ 2;
t45 = t55 * t84;
t43 = t54 * t87;
t42 = t55 * t87;
t33 = t80 * t111 - t48 * t78;
t32 = t78 * t111 + t48 * t80;
t31 = t55 * t49;
t30 = t86 * t58 + t83 * t59;
t29 = -t83 * t58 + t86 * t59;
t25 = -pkin(10) * t110 + t28;
t24 = t84 * pkin(5) + t57 + (pkin(10) * t87 - t52) * t78;
t19 = t48 * pkin(3) + t91;
t18 = t83 * t32 - t86 * t33;
t17 = -t86 * t32 - t83 * t33;
t9 = t83 * t24 + t86 * t25;
t8 = t86 * t24 - t83 * t25;
t7 = -t32 * pkin(5) + t16;
t4 = t32 * pkin(10) + t6;
t3 = t49 * pkin(5) + t33 * pkin(10) + t5;
t2 = t83 * t3 + t86 * t4;
t1 = t86 * t3 - t83 * t4;
t10 = [1, 0, 0, t73 * t85 ^ 2, 0.2e1 * t85 * t115, 0.2e1 * t79 * t109, t81 * t103, t81 ^ 2, 0.2e1 * pkin(1) * t115 + 0.2e1 * t50 * t81, -0.2e1 * t73 * t121 - 0.2e1 * t51 * t81, t46, t48 * t127, t111 * t127, t48 * t103, t73 * t88 ^ 2, -0.2e1 * t22 * t111 + 0.2e1 * t37 * t48, 0.2e1 * t23 * t111 + 0.2e1 * t37 * t49, 0.2e1 * t20 * t48 + 0.2e1 * t21 * t49, -0.2e1 * t21 * t111 - 0.2e1 * t19 * t48, 0.2e1 * t20 * t111 - 0.2e1 * t19 * t49, t19 ^ 2 + t20 ^ 2 + t21 ^ 2, -0.2e1 * t16 * t32 + 0.2e1 * t5 * t49, -0.2e1 * t16 * t33 - 0.2e1 * t6 * t49, 0.2e1 * t6 * t32 + 0.2e1 * t5 * t33, t16 ^ 2 + t5 ^ 2 + t6 ^ 2, t18 ^ 2, -0.2e1 * t18 * t17, 0.2e1 * t18 * t49, t17 * t127, t46, 0.2e1 * t1 * t49 + 0.2e1 * t7 * t17, 0.2e1 * t7 * t18 - 0.2e1 * t2 * t49; 0, 0, 0, 0, 0, t112, t111, t81, t50, -t51, t41, -t84 * t48 + t49 * t87, -t101, -t100, 0, -pkin(2) * t48 - t37 * t87 + t97, -pkin(2) * t49 + t37 * t84 + t96 (-t48 * t87 + t41) * pkin(9) + t93, t19 * t87 - t60 * t48 - t97, -t19 * t84 - t60 * t49 - t96, t93 * pkin(9) + t19 * t60, t110 * t16 + t27 * t49 - t62 * t32 + t5 * t84, -t113 * t16 - t28 * t49 - t62 * t33 - t6 * t84, t27 * t33 + t28 * t32 + (t5 * t78 - t6 * t80) * t87, t16 * t62 + t5 * t27 + t6 * t28, -t18 * t43, t43 * t17 - t18 * t42, t18 * t84 - t43 * t49, -t17 * t84 - t42 * t49, t41, t1 * t84 + t47 * t17 + t7 * t42 + t8 * t49, t47 * t18 - t2 * t84 - t7 * t43 - t9 * t49; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t76, t84 * t123, 0, 0, 0, pkin(2) * t123, pkin(2) * t125, 0.2e1 * t105 * pkin(9), t60 * t123, t60 * t125, t105 * pkin(9) ^ 2 + t60 ^ 2, 0.2e1 * t110 * t62 + 0.2e1 * t27 * t84, -0.2e1 * t113 * t62 - 0.2e1 * t28 * t84 (t27 * t78 - t28 * t80) * t123, t27 ^ 2 + t28 ^ 2 + t62 ^ 2, t43 ^ 2, 0.2e1 * t43 * t42, -t43 * t124, -t42 * t124, t76, 0.2e1 * t47 * t42 + 0.2e1 * t8 * t84, -0.2e1 * t47 * t43 - 0.2e1 * t9 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t48, -t111, t22, -t23, -t49 * pkin(3) - qJ(4) * t48, -t22 + 0.2e1 * t66, -0.2e1 * t99 + t23, -t21 * pkin(3) - t20 * qJ(4), -qJ(4) * t32 - t106 * t118 + t16 * t78, -qJ(4) * t33 + t106 * t114 + t16 * t80 (-t106 * t33 - t5) * t80 + (-t106 * t32 - t6) * t78, t16 * qJ(4) - t106 * t95, t18 * t55, -t55 * t17 - t18 * t54, t31, -t117, 0, t67 * t17 + t29 * t49 + t7 * t54, t67 * t18 - t30 * t49 + t7 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t87, 0, -t69, -t70, t94, t69, t70, t94 * pkin(9), t62 * t78 + t80 * t92, t62 * t80 - t78 * t92, -t15, t62 * qJ(4) - t106 * t15, -t43 * t55, -t55 * t42 + t43 * t54, t45, -t116, 0, t29 * t84 + t67 * t42 + t47 * t54, -t30 * t84 - t67 * t43 + t47 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t122, pkin(3) ^ 2 + t89, t78 * t122, t80 * t122, 0.2e1 * t53, t106 ^ 2 * t64 + t89, t55 ^ 2, -0.2e1 * t55 * t54, 0, 0, 0, t54 * t126, t55 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t111, 0, t21, t118, -t114, t78 * t32 + t80 * t33, t95, 0, 0, 0, 0, 0, t31, -t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, t69, t80 * t84, -t78 * t84, 0, t15, 0, 0, 0, 0, 0, t45, -t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, -t64, -t53, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, t16, 0, 0, 0, 0, 0, t17, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, -t113, 0, t62, 0, 0, 0, 0, 0, t42, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t80, 0, qJ(4), 0, 0, 0, 0, 0, t54, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, t49, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t42, t84, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t54, 0, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
