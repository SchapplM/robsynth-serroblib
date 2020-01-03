% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:16
% EndTime: 2019-12-31 17:14:18
% DurationCPUTime: 0.61s
% Computational Cost: add. (736->161), mult. (1107->187), div. (0->0), fcn. (552->8), ass. (0->109)
t58 = qJ(1) + qJ(2);
t52 = cos(t58);
t128 = g(2) * t52;
t51 = sin(t58);
t45 = g(1) * t51;
t134 = t128 - t45;
t59 = sin(qJ(3));
t62 = cos(qJ(3));
t79 = t62 * pkin(3) + t59 * qJ(4);
t28 = -pkin(2) - t79;
t132 = g(1) * t52 + g(2) * t51;
t55 = qJD(1) + qJD(2);
t109 = pkin(1) * qJD(1);
t60 = sin(qJ(2));
t90 = t60 * t109;
t26 = pkin(6) * t55 + t90;
t117 = t59 * t26;
t63 = cos(qJ(2));
t104 = qJD(2) * t63;
t54 = qJDD(1) + qJDD(2);
t96 = qJDD(1) * t60;
t16 = t54 * pkin(6) + (qJD(1) * t104 + t96) * pkin(1);
t12 = t62 * t16;
t94 = qJDD(3) * qJ(4);
t4 = t94 + t12 + (qJD(4) - t117) * qJD(3);
t11 = t59 * t16;
t100 = qJDD(3) * pkin(3);
t101 = qJD(3) * t62;
t130 = t26 * t101 - t100;
t5 = qJDD(4) + t11 + t130;
t133 = t4 * t62 + t5 * t59;
t56 = t59 ^ 2;
t57 = t62 ^ 2;
t110 = t56 + t57;
t131 = t110 * t55;
t65 = qJD(3) ^ 2;
t129 = pkin(6) * t65;
t127 = t54 * pkin(2);
t126 = t55 * pkin(2);
t124 = t63 * pkin(1);
t123 = t28 * t54;
t122 = t28 * t55;
t47 = pkin(1) * t60 + pkin(6);
t121 = t47 * t65;
t120 = t51 * t59;
t119 = t52 * t59;
t118 = t55 * t59;
t116 = t62 * t26;
t115 = t62 * t54;
t114 = g(1) * t120 - g(2) * t119;
t112 = -qJD(2) * t90 + qJDD(1) * t124;
t111 = t56 - t57;
t107 = pkin(6) * qJDD(3);
t106 = qJD(1) * t63;
t105 = qJD(2) * t60;
t103 = qJD(3) * t55;
t102 = qJD(3) * t59;
t97 = qJD(3) * qJ(4);
t18 = t97 + t116;
t99 = t18 * qJD(3);
t98 = t59 * qJD(4);
t95 = qJDD(3) * t47;
t53 = t55 ^ 2;
t93 = t59 * t53 * t62;
t37 = t62 * t45;
t84 = t55 * t90;
t89 = pkin(1) * t106;
t92 = t102 * t89 + t62 * t84 + t37;
t91 = pkin(1) * t104;
t88 = t55 * t105;
t15 = -t112 - t127;
t87 = -t15 - t128;
t86 = t110 * t54;
t27 = -t89 - t126;
t85 = t101 * t27 + t15 * t59 - t114;
t83 = -qJD(3) * pkin(3) + qJD(4);
t82 = -t127 + t129;
t61 = sin(qJ(1));
t64 = cos(qJ(1));
t80 = g(1) * t61 - g(2) * t64;
t78 = pkin(3) * t59 - qJ(4) * t62;
t14 = t83 + t117;
t77 = t14 * t59 + t18 * t62;
t76 = t112 - t134;
t75 = g(1) * t119 + g(2) * t120 - g(3) * t62 - t11;
t1 = (qJD(3) * t78 - t98) * t55 + t123 - t112;
t74 = -t1 - t123 - t129;
t24 = t28 - t124;
t73 = t24 * t55 - t91;
t72 = -qJDD(4) + t75;
t19 = pkin(3) * t102 - t62 * t97 - t98;
t71 = t101 * t14 - t59 * t99 - t132 + t133;
t10 = pkin(1) * t105 + t19;
t70 = -t10 * t55 - t24 * t54 - t1 - t121;
t48 = -pkin(2) - t124;
t69 = pkin(1) * t88 + t48 * t54 + t121;
t68 = -t95 + (t48 * t55 - t91) * qJD(3);
t67 = (t14 * t62 - t18 * t59) * qJD(3) + t133;
t66 = -t132 * pkin(6) + t134 * t28;
t39 = t59 * t54;
t30 = qJDD(3) * t62 - t59 * t65;
t29 = qJDD(3) * t59 + t62 * t65;
t22 = t27 * t102;
t20 = t78 * t55;
t17 = 0.2e1 * t101 * t118 + t54 * t56;
t8 = -t89 + t122;
t7 = -0.2e1 * t103 * t111 + 0.2e1 * t115 * t59;
t6 = t8 * t102;
t2 = [qJDD(1), t80, g(1) * t64 + g(2) * t61, t54, (t54 * t63 - t88) * pkin(1) + t76, ((-qJDD(1) - t54) * t60 + (-qJD(1) - t55) * t104) * pkin(1) + t132, t17, t7, t29, t30, 0, t22 + t37 + t68 * t59 + (-t69 + t87) * t62, t59 * t69 + t62 * t68 + t85, t37 + t6 + (qJD(3) * t73 - t95) * t59 + (t70 - t128) * t62, t131 * t91 + t47 * t86 + t71, (t95 + (-t73 - t8) * qJD(3)) * t62 + t70 * t59 + t114, t1 * t24 + t8 * t10 + (t104 * t77 + t80) * pkin(1) + t67 * t47 + t66; 0, 0, 0, t54, t76 + t84, (-t96 + (-qJD(2) + t55) * t106) * pkin(1) + t132, t17, t7, t29, t30, 0, t22 + (-pkin(2) * t103 - t107) * t59 + (-t82 + t87) * t62 + t92, (-t107 + (t89 - t126) * qJD(3)) * t62 + (t82 - t84) * t59 + t85, t6 + (t103 * t28 - t107) * t59 + (-t19 * t55 - t128 + t74) * t62 + t92, pkin(6) * t86 - t131 * t89 + t71, (t107 + (-t8 - t89 - t122) * qJD(3)) * t62 + ((-t19 + t90) * t55 + t74) * t59 + t114, t1 * t28 + t8 * t19 + (-t60 * t8 - t63 * t77) * t109 + t67 * pkin(6) + t66; 0, 0, 0, 0, 0, 0, -t93, t111 * t53, t39, t115, qJDD(3), -t118 * t27 + t75, g(3) * t59 - t12 + (-t27 * t55 + t132) * t62, 0.2e1 * t100 + (t20 * t62 - t59 * t8) * t55 + t72, -t78 * t54 + ((t18 - t97) * t59 + (-t14 + t83) * t62) * t55, 0.2e1 * t94 + 0.2e1 * qJD(3) * qJD(4) + t12 + (t20 * t55 - g(3)) * t59 + (t55 * t8 - t132) * t62, t4 * qJ(4) - t5 * pkin(3) - t8 * t20 - t14 * t116 - g(3) * t79 + (qJD(4) + t117) * t18 + t132 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t93, t39, -t53 * t56 - t65, t118 * t8 + t130 - t72 - t99;];
tau_reg = t2;
