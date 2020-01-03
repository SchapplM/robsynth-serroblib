% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR10
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:57
% EndTime: 2019-12-31 17:12:00
% DurationCPUTime: 1.07s
% Computational Cost: add. (653->204), mult. (1460->271), div. (0->0), fcn. (864->6), ass. (0->127)
t149 = 2 * qJD(2);
t62 = cos(qJ(2));
t106 = t62 * qJDD(1);
t104 = qJD(1) * qJD(2);
t59 = sin(qJ(2));
t96 = t59 * t104;
t148 = -t106 + t96;
t123 = t59 * qJ(3);
t87 = pkin(2) * t62 + t123;
t81 = pkin(1) + t87;
t20 = t81 * qJD(1);
t147 = t81 * qJDD(1);
t143 = pkin(3) + pkin(5);
t146 = t143 * t59;
t121 = qJD(1) * t59;
t46 = qJD(4) + t121;
t113 = qJD(4) * t46;
t107 = t59 * qJDD(1);
t95 = t62 * t104;
t76 = t95 + t107;
t26 = qJDD(4) + t76;
t61 = cos(qJ(4));
t18 = t61 * t26;
t58 = sin(qJ(4));
t145 = -t58 * t113 + t18;
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t89 = g(1) * t60 - g(2) * t63;
t90 = g(1) * t63 + g(2) * t60;
t108 = qJD(2) * qJ(3);
t120 = qJD(1) * t62;
t50 = pkin(5) * t120;
t37 = -t50 - t108;
t51 = pkin(3) * t120;
t19 = -t37 + t51;
t64 = -pkin(2) - pkin(6);
t144 = t64 * t26 + (t19 - t50 - t51) * t46;
t99 = t58 * t120;
t6 = -qJD(4) * t99 + t58 * qJDD(2) + (qJD(2) * qJD(4) - t148) * t61;
t138 = g(3) * t59;
t55 = g(3) * t62;
t117 = qJD(2) * t58;
t27 = t61 * t120 + t117;
t5 = -qJD(4) * t27 + t61 * qJDD(2) + t148 * t58;
t137 = t5 * t61;
t136 = t27 * t46;
t115 = qJD(2) * t61;
t29 = -t99 + t115;
t135 = t29 * t46;
t134 = t46 * t58;
t133 = t58 * t26;
t132 = t58 * t60;
t131 = t58 * t62;
t130 = t58 * t63;
t129 = t59 * t61;
t128 = t60 * t61;
t127 = t61 * t62;
t126 = t61 * t63;
t56 = t59 ^ 2;
t57 = t62 ^ 2;
t125 = t56 - t57;
t124 = qJ(3) * t62;
t122 = pkin(5) * qJDD(2);
t119 = qJD(2) * t27;
t118 = qJD(2) * t29;
t116 = qJD(2) * t59;
t114 = qJD(2) * t62;
t112 = qJD(4) * t62;
t111 = qJDD(2) * pkin(2);
t110 = t59 * qJD(3);
t49 = pkin(5) * t121;
t109 = pkin(3) * t121 + qJD(3) + t49;
t47 = pkin(5) * t107;
t105 = qJDD(3) + t47;
t103 = qJDD(2) * qJ(3);
t66 = qJD(1) ^ 2;
t102 = t59 * t66 * t62;
t39 = t143 * t62;
t25 = t64 * t62 - pkin(1) - t123;
t45 = pkin(2) * t96;
t85 = pkin(6) * t59 - t124;
t71 = t85 * qJD(2) - t110;
t1 = t71 * qJD(1) + t25 * qJDD(1) + t45;
t44 = pkin(5) * t95;
t97 = t44 + t105;
t8 = t76 * pkin(3) + t64 * qJDD(2) + t97;
t98 = -t58 * t1 + t61 * t8;
t14 = t64 * qJD(2) + t109;
t94 = qJD(4) * t14 + t1;
t11 = t25 * qJD(1);
t93 = qJD(4) * t11 - t8;
t92 = -qJD(2) * pkin(2) + qJD(3);
t48 = pkin(5) * t106;
t91 = -t48 - t103;
t86 = pkin(2) * t59 - t124;
t3 = t11 * t61 + t14 * t58;
t35 = t49 + t92;
t83 = t35 * t62 + t37 * t59;
t79 = -t61 * t113 - t133;
t78 = -0.2e1 * pkin(1) * t104 - t122;
t77 = -t62 * t108 - t110;
t75 = pkin(1) * t66 + t90;
t65 = qJD(2) ^ 2;
t74 = pkin(5) * t65 - t89;
t73 = t20 * t149 + t122;
t72 = 0.2e1 * qJDD(1) * pkin(1) - t74;
t70 = -t20 * t121 - t90 * t59 + t105 + t55;
t52 = pkin(2) * t116;
t17 = t52 + t77;
t7 = t77 * qJD(1) - t147 + t45;
t69 = qJD(1) * t17 - t147 + t7 + t74;
t53 = pkin(2) * t121;
t9 = pkin(3) * t106 + (-qJD(1) * t146 + qJD(3)) * qJD(2) - t91;
t68 = -t138 + t9 - t90 * t62 + (t85 * qJD(1) - qJD(4) * t64 + t53) * t46;
t12 = (-qJD(3) + t49) * qJD(2) + t91;
t15 = t97 - t111;
t67 = t83 * qJD(2) - t12 * t62 + t15 * t59 - t90;
t34 = qJD(2) * t39;
t32 = qJD(2) * t146;
t30 = -qJ(3) * t120 + t53;
t24 = -t59 * t132 + t126;
t23 = t59 * t128 + t130;
t22 = t59 * t130 + t128;
t21 = t59 * t126 - t132;
t10 = t52 + t71;
t2 = -t11 * t58 + t14 * t61;
t4 = [qJDD(1), t89, t90, qJDD(1) * t56 + 0.2e1 * t59 * t95, -0.2e1 * t125 * t104 + 0.2e1 * t59 * t106, qJDD(2) * t59 + t62 * t65, qJDD(2) * t62 - t59 * t65, 0, t78 * t59 + t72 * t62, -t72 * t59 + t78 * t62, (t56 + t57) * qJDD(1) * pkin(5) + t67, t73 * t59 + t69 * t62, -t69 * t59 + t73 * t62, t67 * pkin(5) - t20 * t17 + (-t7 + t89) * t81, -t5 * t131 + (-t61 * t112 + t58 * t116) * t29, (-t27 * t58 + t29 * t61) * t116 + (-t137 + t58 * t6 + (t27 * t61 + t29 * t58) * qJD(4)) * t62, (t46 * t117 + t5) * t59 + (t79 + t118) * t62, (t46 * t115 - t6) * t59 + (-t119 - t145) * t62, t46 * t114 + t26 * t59, (-t58 * t10 + t61 * t34) * t46 + (t146 * t61 - t25 * t58) * t26 + t98 * t59 - t32 * t27 + t39 * t6 + t9 * t127 - g(1) * t24 - g(2) * t22 + (-t129 * t19 + t2 * t62) * qJD(2) + ((-t146 * t58 - t25 * t61) * t46 - t3 * t59 - t19 * t131) * qJD(4), -t3 * t114 + g(1) * t23 - g(2) * t21 - t32 * t29 + t39 * t5 + (-(qJD(4) * t146 + t10) * t46 - t25 * t26 - t94 * t59 - t19 * t112) * t61 + (-(-qJD(4) * t25 + t34) * t46 - t146 * t26 - t9 * t62 + (qJD(2) * t19 + t93) * t59) * t58; 0, 0, 0, -t102, t125 * t66, t107, t106, qJDD(2), t75 * t59 - t47 - t55, t75 * t62 + t138 - t48, -t86 * qJDD(1) + ((-t37 - t108) * t59 + (-t35 + t92) * t62) * qJD(1), -t30 * t120 - 0.2e1 * t111 + t70, 0.2e1 * t103 + qJD(3) * t149 + t48 + (qJD(1) * t30 - g(3)) * t59 + (-qJD(1) * t20 - t90) * t62, -t83 * qJD(1) * pkin(5) - t15 * pkin(2) - g(3) * t87 - t12 * qJ(3) - t37 * qJD(3) + t20 * t30 + t90 * t86, -t134 * t29 + t137, (-t6 - t135) * t61 + (-t5 + t136) * t58, (-t59 * t134 - t29 * t62) * qJD(1) + t145, (-t46 * t129 + t27 * t62) * qJD(1) + t79, -t46 * t120, qJ(3) * t6 + t109 * t27 - t2 * t120 + t144 * t61 + t68 * t58, qJ(3) * t5 + t109 * t29 + t3 * t120 - t144 * t58 + t68 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, qJDD(2) + t102, -t56 * t66 - t65, qJD(2) * t37 - t111 + t44 + t70, 0, 0, 0, 0, 0, -t134 * t46 - t119 + t18, -t46 ^ 2 * t61 - t118 - t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t27, -t27 ^ 2 + t29 ^ 2, t5 + t136, t135 - t6, t26, -g(1) * t21 - g(2) * t23 + g(3) * t127 - t19 * t29 + t98 + (-qJD(4) + t46) * t3, g(1) * t22 - g(2) * t24 + t19 * t27 + t2 * t46 - t94 * t61 + (t93 - t55) * t58;];
tau_reg = t4;
