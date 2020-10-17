% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:27:25
% EndTime: 2019-05-07 10:27:29
% DurationCPUTime: 1.15s
% Computational Cost: add. (1507->138), mult. (2998->231), div. (0->0), fcn. (3689->10), ass. (0->111)
t92 = sin(pkin(11));
t93 = cos(pkin(11));
t106 = t92 ^ 2 + t93 ^ 2;
t107 = t106 * qJ(4);
t95 = sin(qJ(5));
t99 = cos(qJ(5));
t67 = t95 * t92 - t99 * t93;
t124 = t67 * pkin(5);
t121 = t93 * pkin(4);
t114 = cos(qJ(3));
t103 = t114 * pkin(2);
t85 = -t103 - pkin(3);
t73 = t85 - t121;
t57 = t73 + t124;
t131 = 0.2e1 * t57;
t83 = -pkin(3) - t121;
t58 = t83 + t124;
t130 = 0.2e1 * t58;
t115 = cos(qJ(2));
t96 = sin(qJ(3));
t97 = sin(qJ(2));
t70 = -t114 * t115 + t96 * t97;
t129 = -0.2e1 * t70;
t128 = 0.2e1 * t70;
t127 = 0.2e1 * t73;
t126 = 0.2e1 * t83;
t86 = -t115 * pkin(2) - pkin(1);
t125 = 0.2e1 * t86;
t68 = t99 * t92 + t95 * t93;
t123 = t68 * pkin(10);
t122 = t70 * pkin(5);
t94 = sin(qJ(6));
t120 = t94 * pkin(5);
t119 = t96 * pkin(2);
t98 = cos(qJ(6));
t118 = t98 * pkin(5);
t71 = t114 * t97 + t96 * t115;
t38 = t68 * t71;
t111 = t93 * t71;
t45 = t70 * pkin(3) - t71 * qJ(4) + t86;
t78 = (-pkin(8) - pkin(7)) * t97;
t104 = t115 * pkin(7);
t79 = t115 * pkin(8) + t104;
t56 = t114 * t79 + t96 * t78;
t25 = t93 * t45 - t92 * t56;
t17 = t70 * pkin(4) - pkin(9) * t111 + t25;
t112 = t92 * t71;
t26 = t92 * t45 + t93 * t56;
t20 = -pkin(9) * t112 + t26;
t8 = t95 * t17 + t99 * t20;
t6 = -t38 * pkin(10) + t8;
t117 = t98 * t6;
t116 = pkin(3) - t85;
t55 = -t114 * t78 + t96 * t79;
t113 = t55 * t93;
t110 = t57 + t58;
t109 = t73 + t83;
t82 = qJ(4) + t119;
t108 = t106 * t82;
t105 = 0.2e1 * t115;
t39 = t67 * t71;
t7 = t99 * t17 - t95 * t20;
t4 = t39 * pkin(10) + t122 + t7;
t1 = t98 * t4 - t94 * t6;
t62 = (-pkin(9) - t82) * t92;
t89 = t93 * pkin(9);
t63 = t93 * t82 + t89;
t42 = t99 * t62 - t95 * t63;
t74 = (-pkin(9) - qJ(4)) * t92;
t75 = t93 * qJ(4) + t89;
t51 = t99 * t74 - t95 * t75;
t34 = pkin(4) * t112 + t55;
t102 = -pkin(3) * t71 - qJ(4) * t70;
t9 = -t25 * t92 + t26 * t93;
t43 = t95 * t62 + t99 * t63;
t101 = -t70 * t82 + t71 * t85;
t52 = t95 * t74 + t99 * t75;
t69 = t70 ^ 2;
t65 = t68 ^ 2;
t64 = t67 * pkin(10);
t54 = t68 * t70;
t53 = t67 * t70;
t50 = -0.2e1 * t68 * t67;
t49 = t55 * t92;
t47 = -t94 * t67 + t98 * t68;
t46 = t98 * t67 + t94 * t68;
t44 = t47 ^ 2;
t36 = t52 - t64;
t35 = t51 - t123;
t33 = t47 * t70;
t32 = t46 * t70;
t31 = t43 - t64;
t30 = t42 - t123;
t29 = t39 * t68;
t28 = t34 * t68;
t27 = t34 * t67;
t24 = -0.2e1 * t47 * t46;
t23 = t38 * pkin(5) + t34;
t22 = -t94 * t38 - t98 * t39;
t21 = t98 * t38 - t94 * t39;
t19 = t94 * t35 + t98 * t36;
t18 = t98 * t35 - t94 * t36;
t15 = -t68 * t38 + t39 * t67;
t14 = t94 * t30 + t98 * t31;
t13 = t98 * t30 - t94 * t31;
t12 = t23 * t47;
t11 = t23 * t46;
t10 = t22 * t47;
t5 = -t47 * t21 - t22 * t46;
t2 = t94 * t4 + t117;
t3 = [1, 0, 0, t97 ^ 2, t97 * t105, 0, 0, 0, pkin(1) * t105, -0.2e1 * pkin(1) * t97, t71 ^ 2, t71 * t129, 0, 0, 0, t70 * t125, t71 * t125, 0.2e1 * t55 * t112 + 0.2e1 * t25 * t70, 0.2e1 * t55 * t111 - 0.2e1 * t26 * t70, 0.2e1 * (-t25 * t93 - t26 * t92) * t71, t25 ^ 2 + t26 ^ 2 + t55 ^ 2, t39 ^ 2, 0.2e1 * t39 * t38, -t39 * t128, t38 * t129, t69, 0.2e1 * t34 * t38 + 0.2e1 * t7 * t70, -0.2e1 * t34 * t39 - 0.2e1 * t8 * t70, t22 ^ 2, -0.2e1 * t22 * t21, t22 * t128, t21 * t129, t69, 0.2e1 * t1 * t70 + 0.2e1 * t23 * t21, -0.2e1 * t2 * t70 + 0.2e1 * t23 * t22; 0, 0, 0, 0, 0, t97, t115, 0, -t97 * pkin(7), -t104, 0, 0, t71, -t70, 0, -t55, -t56, t101 * t92 - t113, t101 * t93 + t49, t9, t55 * t85 + t9 * t82, -t29, t15, t54, -t53, 0, t73 * t38 + t42 * t70 + t27, -t73 * t39 - t43 * t70 + t28, t10, t5, t33, -t32, 0, t13 * t70 + t57 * t21 + t11, -t14 * t70 + t57 * t22 + t12; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t103, -0.2e1 * t119, -0.2e1 * t85 * t93, 0.2e1 * t85 * t92, 0.2e1 * t108, t106 * t82 ^ 2 + t85 ^ 2, t65, t50, 0, 0, 0, t67 * t127, t68 * t127, t44, t24, 0, 0, 0, t46 * t131, t47 * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t70, 0, -t55, -t56, t102 * t92 - t113, t102 * t93 + t49, t9, -t55 * pkin(3) + t9 * qJ(4), -t29, t15, t54, -t53, 0, t83 * t38 + t51 * t70 + t27, -t83 * t39 - t52 * t70 + t28, t10, t5, t33, -t32, 0, t18 * t70 + t58 * t21 + t11, -t19 * t70 + t58 * t22 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t103, -t119, t116 * t93, -t116 * t92, t107 + t108, -t85 * pkin(3) + t82 * t107, t65, t50, 0, 0, 0, t109 * t67, t109 * t68, t44, t24, 0, 0, 0, t110 * t46, t110 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t93, -0.2e1 * pkin(3) * t92, 0.2e1 * t107, t106 * qJ(4) ^ 2 + pkin(3) ^ 2, t65, t50, 0, 0, 0, t67 * t126, t68 * t126, t44, t24, 0, 0, 0, t46 * t130, t47 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t111, 0, t55, 0, 0, 0, 0, 0, t38, -t39, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t92, 0, t85, 0, 0, 0, 0, 0, t67, t68, 0, 0, 0, 0, 0, t46, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t92, 0, -pkin(3), 0, 0, 0, 0, 0, t67, t68, 0, 0, 0, 0, 0, t46, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t38, t70, t7, -t8, 0, 0, t22, -t21, t70, t70 * t118 + t1, -t117 + (-t4 - t122) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t67, 0, t42, -t43, 0, 0, t47, -t46, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t67, 0, t51, -t52, 0, 0, t47, -t46, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t118, -0.2e1 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, t70, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t46, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t46, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t118, -t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
