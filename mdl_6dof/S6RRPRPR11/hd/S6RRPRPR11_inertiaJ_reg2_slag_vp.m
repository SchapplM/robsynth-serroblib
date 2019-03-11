% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t79 = sin(qJ(2));
t71 = t79 ^ 2;
t81 = cos(qJ(2));
t73 = t81 ^ 2;
t131 = t71 + t73;
t108 = cos(qJ(6));
t75 = sin(pkin(10));
t76 = cos(pkin(10));
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t44 = t75 * t80 + t76 * t78;
t46 = -t75 * t78 + t76 * t80;
t77 = sin(qJ(6));
t86 = t108 * t44 + t77 * t46;
t130 = t86 ^ 2;
t87 = t108 * t46 - t77 * t44;
t129 = t87 ^ 2;
t128 = (t44 * t75 + t46 * t76) * pkin(4);
t102 = t78 * t81;
t99 = t80 * t81;
t33 = t75 * t102 - t76 * t99;
t34 = t44 * t81;
t15 = -t108 * t33 - t77 * t34;
t127 = t15 * t86;
t17 = -t108 * t34 + t77 * t33;
t126 = t17 * t87;
t125 = t86 * t79;
t124 = t87 * t79;
t118 = t44 ^ 2;
t43 = t46 ^ 2;
t123 = t43 + t118;
t122 = t129 + t130;
t109 = t79 * pkin(4);
t66 = t79 * pkin(7);
t53 = t79 * pkin(3) + t66;
t47 = t80 * t53;
t82 = -pkin(2) - pkin(8);
t94 = -t79 * qJ(3) - pkin(1);
t42 = t82 * t81 + t94;
t93 = qJ(5) * t81 - t42;
t18 = t93 * t78 + t109 + t47;
t104 = t78 * t53;
t20 = -t93 * t80 + t104;
t8 = t76 * t18 - t75 * t20;
t4 = t79 * pkin(5) + t34 * pkin(9) + t8;
t9 = t75 * t18 + t76 * t20;
t5 = t33 * pkin(9) + t9;
t1 = t108 * t4 - t77 * t5;
t2 = t108 * t5 + t77 * t4;
t121 = t1 * t87 + t2 * t86;
t50 = (-qJ(5) + t82) * t78;
t63 = t80 * t82;
t51 = -t80 * qJ(5) + t63;
t29 = -t75 * t50 + t76 * t51;
t13 = -t46 * pkin(9) + t29;
t30 = t76 * t50 + t75 * t51;
t14 = -t44 * pkin(9) + t30;
t6 = t108 * t13 - t77 * t14;
t7 = t108 * t14 + t77 * t13;
t120 = t6 * t87 + t7 * t86;
t111 = t75 * pkin(4);
t110 = t76 * pkin(4);
t60 = pkin(5) + t110;
t36 = t108 * t60 - t77 * t111;
t37 = t108 * t111 + t77 * t60;
t119 = t36 * t87 + t37 * t86;
t61 = t78 * pkin(4) + qJ(3);
t32 = t44 * pkin(5) + t61;
t117 = 0.2e1 * t32;
t116 = 0.2e1 * t61;
t115 = -0.2e1 * t79;
t114 = 0.2e1 * t79;
t113 = 0.2e1 * t81;
t112 = 0.2e1 * qJ(3);
t107 = t44 * t33;
t106 = t44 * t79;
t105 = t46 * t34;
t103 = t78 * t79;
t101 = t79 * t81;
t100 = t80 * t78;
t98 = t131 * pkin(7) ^ 2;
t68 = t81 * pkin(7);
t54 = t81 * pkin(3) + t68;
t70 = t78 ^ 2;
t72 = t80 ^ 2;
t56 = t70 + t72;
t97 = t81 * qJ(3);
t96 = -0.2e1 * t101;
t95 = t78 * t99;
t40 = pkin(4) * t99 + t54;
t92 = t9 * t44 + t8 * t46;
t91 = -t79 * pkin(2) + t97;
t27 = -t78 * t42 + t47;
t28 = t80 * t42 + t104;
t10 = t27 * t80 + t28 * t78;
t90 = t29 * t46 + t30 * t44;
t88 = t79 * t82 + t97;
t83 = qJ(3) ^ 2;
t62 = t80 * t79;
t58 = 0.2e1 * t101;
t52 = -t81 * pkin(2) + t94;
t49 = 0.2e1 * t131 * pkin(7);
t48 = t56 * t82;
t39 = t46 * t79;
t21 = -t33 * pkin(5) + t40;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t71, t58, 0, t73, 0, 0, pkin(1) * t113, pkin(1) * t115, t49, pkin(1) ^ 2 + t98, 0, 0, 0, t71, t58, t73, t49, t52 * t113, t52 * t115, t52 ^ 2 + t98, t70 * t73, 0.2e1 * t73 * t100, t78 * t96, t72 * t73, t80 * t96, t71, 0.2e1 * t27 * t79 + 0.2e1 * t54 * t99, -0.2e1 * t54 * t102 - 0.2e1 * t28 * t79 (t27 * t78 - t28 * t80) * t113, t27 ^ 2 + t28 ^ 2 + t54 ^ 2, t34 ^ 2, -0.2e1 * t34 * t33, -t34 * t114, t33 ^ 2, t33 * t114, t71, -0.2e1 * t40 * t33 + 0.2e1 * t8 * t79, -0.2e1 * t40 * t34 - 0.2e1 * t9 * t79, 0.2e1 * t9 * t33 + 0.2e1 * t8 * t34, t40 ^ 2 + t8 ^ 2 + t9 ^ 2, t17 ^ 2, -0.2e1 * t17 * t15, t17 * t114, t15 ^ 2, t15 * t115, t71, 0.2e1 * t1 * t79 + 0.2e1 * t21 * t15, 0.2e1 * t21 * t17 - 0.2e1 * t2 * t79, -0.2e1 * t1 * t17 - 0.2e1 * t2 * t15, t1 ^ 2 + t2 ^ 2 + t21 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, t81, 0, -t66, -t68, 0, 0, 0, -t79, -t81, 0, 0, 0, t91, t66, t68, t91 * pkin(7), -t95 (t70 - t72) * t81, t62, t95, -t103, 0, t54 * t78 + t88 * t80, t54 * t80 - t88 * t78, -t10, t54 * qJ(3) + t10 * t82, -t105, t46 * t33 + t34 * t44, t39, -t107, -t106, 0, t29 * t79 - t61 * t33 + t40 * t44, -t30 * t79 - t61 * t34 + t40 * t46, t29 * t34 + t30 * t33 - t92, t8 * t29 + t9 * t30 + t40 * t61, t126, -t15 * t87 - t17 * t86, t124, t127, -t125, 0, t32 * t15 + t21 * t86 + t6 * t79, t32 * t17 + t21 * t87 - t7 * t79, -t7 * t15 - t6 * t17 - t121, t1 * t6 + t2 * t7 + t21 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(2), t112, pkin(2) ^ 2 + t83, t72, -0.2e1 * t100, 0, t70, 0, 0, t78 * t112, t80 * t112, -0.2e1 * t48, t56 * t82 ^ 2 + t83, t43, -0.2e1 * t46 * t44, 0, t118, 0, 0, t44 * t116, t46 * t116, -0.2e1 * t90, t29 ^ 2 + t30 ^ 2 + t61 ^ 2, t129, -0.2e1 * t87 * t86, 0, t130, 0, 0, t86 * t117, t87 * t117, -0.2e1 * t120, t32 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, t66, 0, 0, 0, 0, 0, 0, t62, -t103, 0, t10, 0, 0, 0, 0, 0, 0, t39, -t106, t105 + t107, t92, 0, 0, 0, 0, 0, 0, t124, -t125, -t126 - t127, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t56, t48, 0, 0, 0, 0, 0, 0, 0, 0, -t123, t90, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, 0, -t99, t79, t27, -t28, 0, 0, 0, 0, -t34, 0, t33, t79, t76 * t109 + t8, -t75 * t109 - t9 (t33 * t75 + t34 * t76) * pkin(4) (t75 * t9 + t76 * t8) * pkin(4), 0, 0, t17, 0, -t15, t79, t36 * t79 + t1, -t37 * t79 - t2, -t37 * t15 - t36 * t17, t1 * t36 + t2 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, -t78, 0, t63, -t78 * t82, 0, 0, 0, 0, t46, 0, -t44, 0, t29, -t30, -t128 (t29 * t76 + t30 * t75) * pkin(4), 0, 0, t87, 0, -t86, 0, t6, -t7, -t119, t6 * t36 + t7 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t78, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t44, 0, t128, 0, 0, 0, 0, 0, 0, t87, -t86, 0, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t110, -0.2e1 * t111, 0 (t75 ^ 2 + t76 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t37, 0, t36 ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t34, 0, t40, 0, 0, 0, 0, 0, 0, t15, t17, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t46, 0, t61, 0, 0, 0, 0, 0, 0, t86, t87, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t15, t79, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, -t86, 0, t6, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t86, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, -t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t3;
