% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:41
% EndTime: 2019-12-31 20:53:43
% DurationCPUTime: 0.71s
% Computational Cost: add. (1008->196), mult. (1698->228), div. (0->0), fcn. (767->4), ass. (0->128)
t76 = qJD(1) + qJD(2);
t83 = cos(qJ(3));
t134 = t76 * t83;
t125 = pkin(1) * qJD(1);
t82 = sin(qJ(2));
t110 = t82 * t125;
t49 = t76 * pkin(7) + t110;
t40 = t83 * t49;
t27 = pkin(4) * t134 + t40;
t114 = qJD(3) * qJ(4);
t99 = -qJD(5) - t114;
t17 = -t99 + t27;
t81 = sin(qJ(3));
t121 = qJD(3) * t81;
t124 = pkin(1) * qJD(2);
t105 = qJD(1) * t124;
t84 = cos(qJ(2));
t100 = t84 * t105;
t55 = t83 * t100;
t77 = qJD(3) * qJD(4);
t129 = t55 + t77;
t15 = t49 * t121 - t129;
t120 = qJD(3) * t83;
t54 = t81 * t100;
t20 = t49 * t120 + t54;
t144 = -t15 * t83 + t20 * t81;
t80 = -pkin(3) - qJ(5);
t143 = qJD(3) * t80;
t39 = t81 * t49;
t142 = -qJD(4) - t39;
t85 = qJD(3) ^ 2;
t141 = pkin(7) * t85;
t140 = t76 * pkin(2);
t139 = t84 * pkin(1);
t122 = t81 * qJ(4);
t38 = t80 * t83 - pkin(2) - t122;
t137 = t38 * t76;
t96 = -t83 * pkin(3) - t122;
t51 = -pkin(2) + t96;
t136 = t51 * t76;
t135 = t76 * t81;
t133 = t85 * t81;
t71 = t85 * t83;
t109 = t84 * t125;
t50 = -t109 - t140;
t64 = t82 * t105;
t132 = t50 * t120 + t81 * t64;
t101 = t76 * t110;
t119 = qJD(3) * t84;
t106 = t81 * t119;
t131 = t83 * t101 + t106 * t125;
t130 = t55 + 0.2e1 * t77;
t108 = t76 * t121;
t128 = pkin(3) * t108 + t64;
t78 = t81 ^ 2;
t79 = t83 ^ 2;
t127 = t78 - t79;
t126 = t78 + t79;
t123 = qJ(4) * t83;
t32 = -t40 - t114;
t118 = t32 * qJD(3);
t117 = t81 * qJD(4);
t116 = -qJD(1) - t76;
t26 = -pkin(4) * t135 - t39;
t115 = qJD(4) - t26;
t75 = t76 ^ 2;
t113 = t81 * t75 * t83;
t107 = t76 * t120;
t112 = pkin(4) * t107 + t20;
t111 = t84 * t124;
t70 = t82 * t124;
t92 = qJ(5) * t81 - t123;
t1 = (t92 * qJD(3) - t83 * qJD(5) - t117) * t76 + t128;
t68 = pkin(3) * t121;
t18 = qJ(5) * t121 + t99 * t83 - t117 + t68;
t13 = t70 + t18;
t104 = -t13 * t76 - t1;
t103 = -t18 * t76 - t1;
t19 = -t109 + t136;
t102 = -t19 - t136;
t98 = (-pkin(4) * t76 - t49) * t81;
t88 = -t83 * t114 - t117;
t33 = t68 + t88;
t97 = t33 * t76 + t141;
t28 = t33 + t70;
t66 = t82 * pkin(1) + pkin(7);
t95 = t28 * t76 + t66 * t85;
t29 = -qJD(3) * pkin(3) - t142;
t94 = t29 * t83 + t32 * t81;
t93 = t29 * t81 - t32 * t83;
t91 = t81 * t118 + t29 * t120 + t144;
t90 = t126 * t109;
t14 = t115 + t143;
t6 = qJD(3) * t98 + t129;
t7 = -qJD(3) * qJD(5) + t112;
t89 = t14 * t120 - t17 * t121 + t6 * t83 + t7 * t81;
t41 = t51 - t139;
t87 = qJD(3) * (-t41 * t76 + t111 - t19);
t86 = t94 * qJD(3) + t144;
t74 = t83 * pkin(4);
t73 = t81 * pkin(4);
t69 = pkin(4) * t120;
t67 = -pkin(2) - t139;
t63 = pkin(3) * t135;
t61 = qJD(4) * t134;
t58 = t83 * pkin(7) + t74;
t57 = t81 * pkin(7) + t73;
t52 = -t78 * t75 - t85;
t48 = 0.2e1 * t81 * t107;
t46 = t81 * t101;
t45 = pkin(7) * t120 + t69;
t44 = (-pkin(4) - pkin(7)) * t121;
t43 = t83 * t66 + t74;
t42 = t81 * t66 + t73;
t36 = t50 * t121;
t34 = -t76 * t123 + t63;
t31 = t38 - t139;
t30 = -0.2e1 * t127 * t76 * qJD(3);
t23 = t92 * t76 + t63;
t22 = t81 * t111 + t66 * t120 + t69;
t21 = t83 * t111 + (-pkin(4) - t66) * t121;
t12 = t19 * t135;
t11 = t88 * t76 + t128;
t10 = -t109 + t137;
t8 = t11 * t83;
t5 = t10 * t121;
t4 = t10 * t134;
t2 = [0, 0, 0, 0, -t76 * t70 - t64, t116 * t111, t48, t30, t71, -t133, 0, t67 * t108 - t66 * t71 + t36 + (t116 * t83 * t82 - t106) * t124, t67 * t107 + t66 * t133 + (-t83 * t119 + t82 * t135) * t124 + t132, t126 * t76 * t111 + t91, t81 * t87 + t95 * t83 + t8, (-t11 - t95) * t81 + t83 * t87, t11 * t41 + t93 * t111 + t19 * t28 + t86 * t66, (t21 * t83 + t22 * t81 + (t42 * t83 - t43 * t81) * qJD(3)) * t76 + t89, t104 * t81 + (t21 + (-t31 * t76 - t10) * t83) * qJD(3), t5 + t104 * t83 + (t31 * t135 - t22) * qJD(3), t1 * t31 + t10 * t13 + t14 * t22 + t17 * t21 + t7 * t42 + t6 * t43; 0, 0, 0, 0, -t64 + t101, (-qJD(2) + t76) * t109, t48, t30, t71, -t133, 0, -pkin(2) * t108 + t36 + (-t64 - t141) * t83 + t131, pkin(7) * t133 - t46 + (t109 - t140) * t120 + t132, -t76 * t90 + t91, t102 * t121 + t97 * t83 - t131 + t8, t46 + (-t11 - t97) * t81 + (t102 - t109) * t120, t11 * t51 + t19 * t33 + (-t19 * t82 - t93 * t84) * t125 + t86 * pkin(7), (t44 * t83 + t45 * t81 + (t57 * t83 - t58 * t81) * qJD(3) - t90) * t76 + t89, t46 + t103 * t81 + (t44 + (-t10 - t109 - t137) * t83) * qJD(3), t5 + t103 * t83 + (t38 * t135 - t45) * qJD(3) + t131, t1 * t38 + t10 * t18 + t14 * t45 + t17 * t44 + t7 * t57 + t6 * t58 + (-t10 * t82 + (-t14 * t81 - t17 * t83) * t84) * t125; 0, 0, 0, 0, 0, 0, -t113, t127 * t75, 0, 0, 0, -t50 * t135 - t54, -t50 * t134 - t55, t61 + (t96 * qJD(3) - t94) * t76, -t34 * t134 + t12 + t54, (t19 * t83 + t34 * t81) * t76 + t130, -t20 * pkin(3) - t15 * qJ(4) + t142 * t32 - t19 * t34 - t29 * t40, t61 + (-t14 - t26 + t143) * t134, t23 * t135 + t4 + (-t26 + t98) * qJD(3) + t130, (-t10 * t81 + t23 * t83) * t76 + (0.2e1 * qJD(5) + t27) * qJD(3) - t112, t6 * qJ(4) - t10 * t23 + t7 * t80 + t115 * t17 + (-qJD(5) - t27) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t52, t12 + t118 + t20, 0, t52, -t113, t10 * t135 + (-qJD(5) - t17) * qJD(3) + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t79 * t75 - t85, t4 + (t14 + t98) * qJD(3) + t129;];
tauc_reg = t2;
