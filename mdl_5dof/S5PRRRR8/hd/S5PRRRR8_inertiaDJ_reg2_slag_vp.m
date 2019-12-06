% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:41
% EndTime: 2019-12-05 17:16:49
% DurationCPUTime: 1.63s
% Computational Cost: add. (1968->187), mult. (5092->356), div. (0->0), fcn. (4913->10), ass. (0->118)
t68 = sin(qJ(5));
t64 = t68 ^ 2;
t72 = cos(qJ(5));
t65 = t72 ^ 2;
t122 = t64 - t65;
t141 = t122 * qJD(5);
t140 = qJD(3) + qJD(4);
t70 = sin(qJ(3));
t139 = t70 ^ 2;
t138 = -pkin(8) - pkin(7);
t137 = cos(qJ(4));
t69 = sin(qJ(4));
t126 = t69 * t70;
t73 = cos(qJ(3));
t52 = t138 * t73;
t36 = t138 * t126 - t137 * t52;
t98 = qJD(3) * t138;
t48 = t70 * t98;
t99 = t137 * t73;
t93 = qJD(3) * t99;
t20 = t36 * qJD(4) - t138 * t93 + t69 * t48;
t100 = t137 * t70;
t35 = -t138 * t100 - t69 * t52;
t136 = t20 * t35;
t66 = sin(pkin(5));
t71 = sin(qJ(2));
t128 = t66 * t71;
t67 = cos(pkin(5));
t40 = -t70 * t128 + t67 * t73;
t41 = t73 * t128 + t67 * t70;
t26 = t137 * t41 + t69 * t40;
t74 = cos(qJ(2));
t119 = qJD(2) * t74;
t105 = t66 * t119;
t34 = t40 * qJD(3) + t73 * t105;
t77 = t41 * qJD(3) + t70 * t105;
t11 = t26 * qJD(4) + t137 * t77 + t69 * t34;
t25 = -t137 * t40 + t69 * t41;
t135 = t25 * t11;
t134 = t25 * t69;
t97 = qJD(4) * t137;
t32 = t140 * t126 - t73 * t97 - t93;
t133 = t32 * t72;
t132 = t35 * t69;
t125 = t69 * t73;
t46 = t100 + t125;
t131 = t46 * t32;
t130 = t46 * t68;
t129 = t46 * t72;
t127 = t66 * t74;
t62 = qJD(5) * t72;
t124 = t20 * t68 + t35 * t62;
t121 = pkin(3) * qJD(4);
t108 = t69 * t121;
t60 = -t137 * pkin(3) - pkin(4);
t123 = t68 * t108 + t60 * t62;
t120 = qJD(2) * t71;
t117 = qJD(5) * t68;
t116 = t70 * qJD(3);
t115 = t73 * qJD(3);
t33 = t140 * t46;
t45 = -t99 + t126;
t114 = 0.2e1 * t45 * t33;
t113 = -0.2e1 * pkin(2) * qJD(3);
t112 = t68 * t133;
t111 = pkin(3) * t116;
t110 = pkin(4) * t117;
t109 = pkin(4) * t62;
t107 = t74 * t116;
t106 = t66 * t120;
t104 = t68 * t62;
t103 = t70 * t115;
t61 = -pkin(3) * t73 - pkin(2);
t102 = t68 * t137;
t101 = t72 * t137;
t43 = t46 ^ 2;
t96 = t43 * t104;
t95 = t66 ^ 2 * t71 * t119;
t94 = pkin(3) * t97;
t92 = t11 * t35 + t25 * t20;
t81 = -pkin(4) * t45 + pkin(9) * t46 - t61;
t78 = t72 * t81;
t15 = -t68 * t36 - t78;
t16 = t72 * t36 - t68 * t81;
t91 = t15 * t72 + t16 * t68;
t86 = t68 * t127 - t72 * t26;
t87 = t72 * t127 + t68 * t26;
t90 = -t68 * t86 - t72 * t87;
t59 = pkin(3) * t69 + pkin(9);
t89 = t45 * t59 - t46 * t60;
t88 = -t72 * t108 + t60 * t117;
t85 = -t32 * t68 + t46 * t62;
t84 = t46 * t117 + t133;
t83 = t45 * t117 - t33 * t72;
t82 = (t64 + t65) * t137;
t80 = (-t137 * t45 + t46 * t69) * qJD(4);
t79 = pkin(4) * t33 + pkin(9) * t32 + t111;
t19 = t35 * qJD(4) - t98 * t125 - t137 * t48;
t3 = qJD(5) * t78 + t36 * t117 + t72 * t19 - t68 * t79;
t4 = -qJD(5) * t16 + t68 * t19 + t72 * t79;
t1 = -t91 * qJD(5) - t3 * t72 - t4 * t68;
t10 = t25 * qJD(4) - t137 * t34 + t69 * t77;
t5 = t87 * qJD(5) + t72 * t10 - t68 * t106;
t6 = t86 * qJD(5) + t68 * t10 + t72 * t106;
t2 = -t90 * qJD(5) - t5 * t72 - t6 * t68;
t76 = pkin(3) * t80 - t32 * t60 - t33 * t59;
t75 = t139 * t105 - t40 * t115 + t34 * t73;
t55 = -0.2e1 * t104;
t54 = 0.2e1 * t104;
t44 = -0.2e1 * t141;
t39 = t82 * t121;
t30 = t35 * t117;
t24 = t33 * t68 + t45 * t62;
t14 = t46 * t141 + t112;
t12 = -0.4e1 * t46 * t104 + t122 * t32;
t8 = -t11 * t72 + t25 * t117;
t7 = t11 * t68 + t25 * t62;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t41 * t34 - 0.2e1 * t40 * t77 - 0.2e1 * t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t26 * t10 + 0.2e1 * t135 - 0.2e1 * t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t5 * t86 - 0.2e1 * t6 * t87 + 0.2e1 * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105, 0, 0, 0, 0, 0, 0, 0, 0, (-t73 * t120 - t107) * t66, (-t74 * t115 + t70 * t120) * t66, t75, -pkin(2) * t106 + t75 * pkin(7), 0, 0, 0, 0, 0, 0, (t45 * t120 - t33 * t74) * t66, (t46 * t120 + t32 * t74) * t66, t10 * t45 + t11 * t46 - t25 * t32 - t26 * t33, -t10 * t36 - t26 * t19 + (-pkin(3) * t107 + t61 * t120) * t66 + t92, 0, 0, 0, 0, 0, 0, t11 * t130 + t25 * t85 - t33 * t87 + t45 * t6, t11 * t129 - t25 * t84 + t33 * t86 + t45 * t5, t90 * t32 + (t5 * t68 - t6 * t72 + (-t68 * t87 + t72 * t86) * qJD(5)) * t46, t15 * t6 - t16 * t5 + t3 * t86 - t4 * t87 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t103, 0.2e1 * (t73 ^ 2 - t139) * qJD(3), 0, -0.2e1 * t103, 0, 0, t70 * t113, t73 * t113, 0, 0, -0.2e1 * t131, 0.2e1 * t32 * t45 - 0.2e1 * t33 * t46, 0, t114, 0, 0, 0.2e1 * t45 * t111 + 0.2e1 * t33 * t61, 0.2e1 * t46 * t111 - 0.2e1 * t32 * t61, 0.2e1 * t19 * t45 + 0.2e1 * t20 * t46 - 0.2e1 * t32 * t35 - 0.2e1 * t33 * t36, 0.2e1 * t61 * t111 - 0.2e1 * t19 * t36 + 0.2e1 * t136, -0.2e1 * t65 * t131 - 0.2e1 * t96, 0.4e1 * t46 * t112 + 0.2e1 * t43 * t141, 0.2e1 * t33 * t129 - 0.2e1 * t45 * t84, -0.2e1 * t64 * t131 + 0.2e1 * t96, -0.2e1 * t33 * t130 - 0.2e1 * t45 * t85, t114, 0.2e1 * t20 * t130 + 0.2e1 * t15 * t33 + 0.2e1 * t35 * t85 + 0.2e1 * t4 * t45, 0.2e1 * t20 * t129 - 0.2e1 * t16 * t33 + 0.2e1 * t3 * t45 - 0.2e1 * t35 * t84, 0.2e1 * t91 * t32 + 0.2e1 * (t3 * t68 - t4 * t72 + (t15 * t68 - t16 * t72) * qJD(5)) * t46, 0.2e1 * t15 * t4 - 0.2e1 * t16 * t3 + 0.2e1 * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t34, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, (-t137 * t11 - t10 * t69 + (t137 * t26 + t134) * qJD(4)) * pkin(3), 0, 0, 0, 0, 0, 0, t8, t7, t2, t11 * t60 + (-t101 * t86 + t102 * t87 + t134) * t121 + t2 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, 0, -t116, 0, -pkin(7) * t115, pkin(7) * t116, 0, 0, 0, 0, -t32, 0, -t33, 0, -t20, t19, (t137 * t32 - t33 * t69 + t80) * pkin(3), (-t137 * t20 - t19 * t69 + (t137 * t36 + t132) * qJD(4)) * pkin(3), -t14, t12, t24, t14, -t83, 0, t30 + (-qJD(5) * t89 - t20) * t72 + t76 * t68, t89 * t117 + t72 * t76 + t124, t1, t20 * t60 + (t101 * t16 - t102 * t15 + t132) * t121 + t1 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t108, -0.2e1 * t94, 0, 0, t54, t44, 0, t55, 0, 0, 0.2e1 * t88, 0.2e1 * t123, 0.2e1 * t39, 0.2e1 * (t59 * t82 + t60 * t69) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, t2, -pkin(4) * t11 + pkin(9) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, -t33, 0, -t20, t19, 0, 0, -t14, t12, t24, t14, -t83, 0, t30 + (pkin(4) * t32 - pkin(9) * t33) * t68 + (-t20 + (-pkin(4) * t46 - pkin(9) * t45) * qJD(5)) * t72, pkin(4) * t84 + pkin(9) * t83 + t124, t1, -pkin(4) * t20 + pkin(9) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t94, 0, 0, t54, t44, 0, t55, 0, 0, t88 - t110, -t109 + t123, t39, (-pkin(4) * t69 + pkin(9) * t82) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t44, 0, t55, 0, 0, -0.2e1 * t110, -0.2e1 * t109, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, 0, -t85, t33, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t117, 0, -t59 * t62 - t68 * t94, t59 * t117 - t72 * t94, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t117, 0, -pkin(9) * t62, pkin(9) * t117, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
