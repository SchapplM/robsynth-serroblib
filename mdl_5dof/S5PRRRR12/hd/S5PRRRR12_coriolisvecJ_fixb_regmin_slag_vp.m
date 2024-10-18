% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:20
% EndTime: 2024-09-28 18:08:21
% DurationCPUTime: 0.73s
% Computational Cost: add. (1413->147), mult. (3605->236), div. (0->0), fcn. (2816->12), ass. (0->114)
t57 = qJD(2) + qJD(3);
t54 = qJD(4) + t57;
t62 = cos(pkin(6));
t127 = t62 * t54;
t49 = qJD(5) + t127;
t107 = qJD(5) - t49;
t61 = sin(pkin(5));
t66 = sin(qJ(3));
t67 = sin(qJ(2));
t70 = cos(qJ(3));
t71 = cos(qJ(2));
t41 = (-t66 * t67 + t70 * t71) * t61;
t139 = t107 * t54;
t63 = cos(pkin(5));
t113 = qJD(1) * t63;
t114 = qJD(1) * t61;
t102 = t67 * t114;
t99 = t71 * t114;
t48 = qJD(2) * pkin(2) + t99;
t37 = t102 * t70 + t48 * t66;
t65 = sin(qJ(4));
t124 = t65 * t37;
t90 = t66 * t102;
t36 = t48 * t70 - t90;
t33 = pkin(3) * t57 + t36;
t69 = cos(qJ(4));
t14 = t33 * t69 - t124;
t13 = pkin(4) * t54 + t14;
t60 = sin(pkin(6));
t11 = t113 * t62 - t13 * t60;
t133 = t54 * t60;
t112 = qJD(3) * t48;
t116 = t57 * t90;
t89 = qJD(2) * t99;
t20 = (t89 + t112) * t70 - t116;
t111 = qJD(4) * t65;
t121 = t67 * t70;
t85 = t66 * t71 + t121;
t81 = t85 * qJD(2);
t73 = (-qJD(3) * t121 - t81) * t114;
t21 = -t112 * t66 + t73;
t94 = -t111 * t37 + t65 * t21;
t5 = (qJD(4) * t33 + t20) * t69 + t94;
t82 = t113 * t60 + t13 * t62;
t138 = -t107 * t82 - t11 * t133 - t5;
t42 = t85 * t61;
t24 = t41 * t69 - t42 * t65;
t137 = -(t24 * t62 + t60 * t63) * t49 + (-t24 * t60 + t62 * t63) * t133;
t119 = t69 * t37;
t15 = -t33 * t65 - t119;
t95 = -t65 * t20 + t21 * t69;
t6 = qJD(4) * t15 + t95;
t64 = sin(qJ(5));
t136 = t6 * t64;
t25 = t41 * t65 + t42 * t69;
t135 = t25 * t49;
t56 = t60 ^ 2;
t134 = t54 * t56;
t132 = t56 * t64;
t68 = cos(qJ(5));
t131 = t56 * t68;
t130 = t60 * t64;
t129 = t60 * t68;
t128 = t61 * qJD(2) ^ 2;
t126 = t62 * t64;
t125 = t62 * t68;
t123 = t65 * t66;
t122 = t66 * t69;
t12 = pkin(10) * t133 - t15;
t120 = t68 * t12;
t110 = qJD(4) * t69;
t38 = t85 * t114;
t39 = qJD(1) * t41;
t52 = pkin(2) * t70 + pkin(3);
t118 = -t69 * t38 - t65 * t39 + t52 * t111 - (-t66 * t110 + (-t65 * t70 - t122) * qJD(3)) * pkin(2);
t117 = -t65 * t38 + t69 * t39 - t52 * t110 - (-t66 * t111 + (t69 * t70 - t123) * qJD(3)) * pkin(2);
t115 = t64 ^ 2 - t68 ^ 2;
t109 = qJD(5) * t64;
t108 = qJD(5) * t68;
t101 = t60 * t109;
t4 = t6 * t125;
t106 = (-t64 * t5 + t4 + (-t64 * t82 - t120) * qJD(5)) * t62 + t6 * t131 + t11 * t101;
t104 = t54 * t132;
t103 = t54 * t131;
t100 = t60 * t108;
t98 = -pkin(3) * t54 - t33;
t97 = -(t6 * t126 + t68 * t5 + (-t64 * t12 + t68 * t82) * qJD(5)) * t62 + t11 * t100;
t96 = t15 * t54 - t6;
t93 = t49 + t127;
t92 = 0.2e1 * qJD(5) * t134;
t16 = -t36 * t65 - t119;
t88 = pkin(3) * t111 + t16;
t17 = t36 * t69 - t124;
t87 = -pkin(3) * t110 + t17;
t83 = (-pkin(2) * t57 - t48) * qJD(3);
t43 = -pkin(2) * t123 + t52 * t69 + pkin(4);
t80 = -t109 * t43 - t118 * t68;
t79 = -t108 * t43 + t118 * t64;
t51 = pkin(3) * t69 + pkin(4);
t76 = -t109 * t51 - t68 * t88;
t75 = -t108 * t51 + t64 * t88;
t55 = t60 * pkin(10);
t53 = t54 ^ 2;
t50 = pkin(3) * t65 + t55;
t44 = t64 * t68 * t92;
t40 = pkin(2) * t122 + t52 * t65 + t55;
t35 = t115 * t92;
t31 = t93 * t100;
t30 = t93 * t101;
t27 = (-qJD(3) * t85 - t81) * t61;
t26 = t57 * t41;
t8 = -qJD(4) * t25 - t65 * t26 + t69 * t27;
t7 = qJD(4) * t24 + t69 * t26 + t65 * t27;
t1 = [0, 0, -t67 * t128, -t71 * t128, 0, t27 * t57, -t26 * t57, 0, t8 * t54, -t7 * t54, 0, 0, 0, 0, 0, (t125 * t8 - t64 * t7) * t49 + t8 * t103 + (-t68 * t135 + t137 * t64) * qJD(5), -(t126 * t8 + t68 * t7) * t49 - t8 * t104 + (t64 * t135 + t137 * t68) * qJD(5); 0, 0, 0, 0, 0, t38 * t57 + t66 * t83 + t73, t39 * t57 + (t83 - t89) * t70 + t116, 0, -t118 * t54 + t6, t117 * t54 - t5, t44, -t35, t31, -t30, 0, (-t108 * t40 + t117 * t64 + t62 * t80) * t49 + t80 * t134 + t106, (t109 * t40 + t117 * t68 + t62 * t79) * t49 + (t54 * t79 - t136) * t56 + t97; 0, 0, 0, 0, 0, t37 * t57 + t21, t36 * t57 - t20, 0, -t16 * t54 + (t65 * t98 - t119) * qJD(4) + t95, t17 * t54 + (qJD(4) * t98 - t20) * t69 - t94, t44, -t35, t31, -t30, 0, (-t108 * t50 + t62 * t76 + t64 * t87) * t49 + t76 * t134 + t106, (t109 * t50 + t62 * t75 + t68 * t87) * t49 + (t54 * t75 - t136) * t56 + t97; 0, 0, 0, 0, 0, 0, 0, 0, -t96, t14 * t54 - t5, t44, -t35, t31, -t30, 0, -(t15 * t125 - t64 * t14) * t49 - t15 * t103 + ((-pkin(4) * t126 - pkin(10) * t129) * t49 - pkin(4) * t104) * qJD(5) + t106, (t15 * t126 + t68 * t14) * t49 + t96 * t132 + (-(pkin(4) * t125 - pkin(10) * t130) * t49 - pkin(4) * t103) * qJD(5) + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t53 * t131, t115 * t56 * t53, t129 * t139, -t130 * t139, 0, -t107 * t120 + t138 * t64 + t4, (t107 * t12 - t62 * t6) * t64 + t138 * t68;];
tauc_reg = t1;
