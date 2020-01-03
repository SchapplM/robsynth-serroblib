% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:25
% EndTime: 2019-12-31 21:49:28
% DurationCPUTime: 0.61s
% Computational Cost: add. (1356->156), mult. (2398->187), div. (0->0), fcn. (1241->6), ass. (0->111)
t63 = sin(qJ(2));
t104 = pkin(1) * qJD(2);
t89 = qJD(1) * t104;
t105 = pkin(1) * qJD(1);
t95 = t63 * t105;
t136 = qJD(3) * t95 + t63 * t89;
t58 = qJD(1) + qJD(2);
t66 = cos(qJ(2));
t92 = t66 * t105;
t42 = t58 * pkin(2) + t92;
t62 = sin(qJ(3));
t65 = cos(qJ(3));
t29 = t62 * t42 + t65 * t95;
t56 = qJD(3) + t58;
t24 = t56 * pkin(8) + t29;
t61 = sin(qJ(4));
t119 = t61 * t24;
t108 = t136 * t62;
t83 = t66 * t89;
t12 = (qJD(3) * t42 + t83) * t65 - t108;
t64 = cos(qJ(4));
t9 = t64 * t12;
t2 = t9 + (qJD(5) - t119) * qJD(4);
t100 = qJD(4) * t64;
t8 = t61 * t12;
t4 = t24 * t100 + t8;
t135 = t2 * t64 + t4 * t61;
t117 = t63 * t65;
t76 = t62 * t66 + t117;
t36 = t76 * t105;
t103 = qJD(3) * t62;
t93 = pkin(2) * t103;
t134 = t56 * (-t36 + t93);
t49 = t62 * t95;
t37 = t65 * t92 - t49;
t102 = qJD(3) * t65;
t94 = pkin(2) * t102;
t133 = -t37 + t94;
t118 = t62 * t63;
t53 = t66 * pkin(1) + pkin(2);
t132 = pkin(1) * t118 - t65 * t53;
t67 = qJD(4) ^ 2;
t131 = pkin(8) * t67;
t130 = t56 * pkin(3);
t129 = t65 * pkin(2);
t17 = t53 * t102 + (-t63 * t103 + (t65 * t66 - t118) * qJD(2)) * pkin(1);
t128 = t17 * t56;
t18 = t53 * t103 + (t76 * qJD(2) + t63 * t102) * pkin(1);
t127 = t18 * t56;
t28 = t65 * t42 - t49;
t126 = t28 * t56;
t125 = t29 * t56;
t35 = pkin(1) * t117 + t62 * t53 + pkin(8);
t124 = t35 * t67;
t43 = -t64 * pkin(4) - t61 * qJ(5) - pkin(3);
t123 = t43 * t56;
t51 = t62 * pkin(2) + pkin(8);
t122 = t51 * t67;
t121 = t56 * t61;
t120 = t56 * t64;
t116 = t64 * t24;
t114 = t67 * t61;
t13 = t42 * t103 + t136 * t65 + t62 * t83;
t23 = -t28 - t130;
t113 = t23 * t100 + t13 * t61;
t101 = qJD(4) * t61;
t112 = t28 * t101 + t29 * t120;
t111 = t37 * t101 + t36 * t120;
t97 = qJD(4) * qJ(5);
t98 = t61 * qJD(5);
t38 = pkin(4) * t101 - t64 * t97 - t98;
t110 = t29 - t38;
t30 = t38 + t93;
t109 = t30 - t36;
t59 = t61 ^ 2;
t60 = t64 ^ 2;
t107 = t59 - t60;
t106 = t59 + t60;
t16 = t97 + t116;
t99 = t16 * qJD(4);
t55 = t56 ^ 2;
t96 = t61 * t55 * t64;
t91 = t56 * t101;
t79 = pkin(4) * t61 - qJ(5) * t64;
t5 = (t79 * qJD(4) - t98) * t56 + t13;
t90 = -t5 - t131;
t88 = -t5 - t122;
t26 = t43 + t132;
t87 = t26 * t56 - t17;
t86 = -qJD(4) * pkin(4) + qJD(5);
t81 = (-qJD(2) + t58) * t105;
t80 = (-qJD(1) - t58) * t104;
t15 = t86 + t119;
t78 = t15 * t61 + t16 * t64;
t77 = t124 + t127;
t75 = qJD(4) * ((-pkin(3) + t132) * t56 - t17);
t7 = t18 + t38;
t74 = -t56 * t7 - t124 - t5;
t40 = t43 - t129;
t73 = t40 * t56 - t94;
t72 = (-pkin(3) - t129) * t56 - t94;
t71 = t15 * t100 - t61 * t99 + t135;
t68 = (t15 * t64 - t16 * t61) * qJD(4) + t135;
t57 = t67 * t64;
t41 = 0.2e1 * t64 * t91;
t33 = t79 * t56;
t31 = -0.2e1 * t107 * t56 * qJD(4);
t20 = t23 * t101;
t11 = -t28 + t123;
t6 = t11 * t101;
t1 = [0, 0, 0, 0, t63 * t80, t66 * t80, 0, -t13 - t127, -t12 - t128, t41, t31, t57, -t114, 0, t20 + t61 * t75 + (-t13 - t77) * t64, t61 * t77 + t64 * t75 + t113, t87 * t101 + t74 * t64 + t6, t106 * t128 + t71, t74 * t61 + (-t11 - t87) * t100, t11 * t7 + t17 * t78 + t5 * t26 + t35 * t68; 0, 0, 0, 0, t63 * t81, t66 * t81, 0, -t13 - t134, t37 * t56 + (-t83 + (-pkin(2) * t56 - t42) * qJD(3)) * t65 + t108, t41, t31, t57, -t114, 0, t20 + t72 * t101 + (-t56 * t93 - t122 - t13) * t64 + t111, (t122 + t134) * t61 + (t37 + t72) * t100 + t113, t6 + t73 * t101 + (-t30 * t56 + t88) * t64 + t111, t133 * t56 * t106 + t71, (-t109 * t56 + t88) * t61 + (-t11 - t37 - t73) * t100, t109 * t11 + t133 * t78 + t5 * t40 + t68 * t51; 0, 0, 0, 0, 0, 0, 0, -t13 + t125, -t12 + t126, t41, t31, t57, -t114, 0, -pkin(3) * t91 + t20 + (-t13 - t131) * t64 + t112, (-t125 + t131) * t61 + (t28 - t130) * t100 + t113, t43 * t91 + t6 + (-t38 * t56 + t90) * t64 + t112, -t106 * t126 + t71, (-t11 - t28 - t123) * t100 + (t110 * t56 + t90) * t61, t68 * pkin(8) - t110 * t11 - t78 * t28 + t5 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, t107 * t55, 0, 0, 0, -t23 * t121 - t8, -t23 * t120 - t9, -t8 + (-t11 * t61 + t33 * t64) * t56, ((t16 - t97) * t61 + (-t15 + t86) * t64) * t56, 0.2e1 * qJD(4) * qJD(5) + t9 + (t11 * t64 + t33 * t61) * t56, -t15 * t116 - t4 * pkin(4) + t2 * qJ(5) - t11 * t33 + (qJD(5) + t119) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, 0, -t59 * t55 - t67, t11 * t121 + t4 - t99;];
tauc_reg = t1;
