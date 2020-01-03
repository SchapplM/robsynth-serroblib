% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:14
% EndTime: 2019-12-31 17:41:17
% DurationCPUTime: 0.86s
% Computational Cost: add. (700->205), mult. (1288->232), div. (0->0), fcn. (666->4), ass. (0->123)
t80 = cos(qJ(3));
t133 = qJ(4) * t80;
t146 = pkin(3) + pkin(4);
t79 = sin(qJ(3));
t149 = t146 * t79;
t93 = t133 - t149;
t66 = t79 * qJ(4);
t109 = pkin(2) + t66;
t154 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t117 = qJD(2) * qJD(3);
t110 = t80 * t117;
t60 = t79 * qJDD(2);
t153 = t110 + t60;
t130 = qJD(2) * t79;
t56 = pkin(6) * t130;
t22 = t80 * qJD(1) - t56;
t152 = qJD(4) - t22;
t112 = t146 * qJD(3);
t121 = qJ(5) * qJD(2);
t12 = t79 * t121 + t22;
t122 = qJD(4) - t12;
t7 = -t112 + t122;
t108 = t146 * qJDD(3);
t23 = t80 * qJD(2) * pkin(6) + t79 * qJD(1);
t14 = -t80 * t121 + t23;
t75 = qJD(3) * qJ(4);
t10 = t14 + t75;
t72 = pkin(7) + qJ(2);
t58 = sin(t72);
t59 = cos(t72);
t138 = g(1) * t59 + g(2) * t58;
t151 = 0.2e1 * t154;
t20 = t75 + t23;
t9 = qJD(5) + (t146 * t80 + t109) * qJD(2);
t132 = qJD(5) + t9;
t150 = t132 * t79;
t17 = -qJD(3) * pkin(3) + t152;
t142 = t59 * t79;
t143 = t58 * t79;
t148 = -g(1) * t142 - g(2) * t143 + g(3) * t80;
t131 = pkin(6) * qJDD(3);
t70 = t80 * pkin(3);
t94 = -t109 - t70;
t18 = t94 * qJD(2);
t137 = t70 + t66;
t24 = -pkin(2) - t137;
t147 = (qJD(2) * t24 + t18) * qJD(3) - t131;
t46 = g(1) * t58;
t145 = g(2) * t59;
t69 = t80 * pkin(4);
t144 = t10 * t80;
t141 = t59 * t80;
t83 = qJD(2) ^ 2;
t140 = t79 * t83;
t139 = pkin(6) - qJ(5);
t76 = t79 ^ 2;
t77 = t80 ^ 2;
t136 = t76 - t77;
t135 = t76 + t77;
t134 = pkin(6) * qJD(3);
t29 = t139 * t80;
t129 = qJD(3) * t29;
t128 = qJD(3) * t79;
t78 = qJDD(2) * pkin(2);
t127 = qJDD(3) * pkin(3);
t126 = t10 * qJD(3);
t125 = t79 * qJD(4);
t124 = t79 * qJD(5);
t123 = t80 * qJD(5);
t120 = qJ(5) * qJD(3);
t62 = t80 * qJDD(2);
t119 = qJ(5) * qJDD(2);
t118 = qJD(1) * qJD(3);
t116 = t80 * t140;
t115 = pkin(6) * t62 + t79 * qJDD(1) + t80 * t118;
t114 = t69 + t137;
t113 = t46 - t145;
t111 = t79 * t117;
t19 = pkin(2) + t114;
t107 = qJD(2) * t19 + t9;
t106 = pkin(3) * t141 + t58 * pkin(6) + t109 * t59;
t105 = t153 * pkin(6) - t80 * qJDD(1) + t79 * t118;
t102 = 0.2e1 * t110;
t101 = pkin(3) * t62 + t153 * qJ(4) + qJD(2) * t125 + t78;
t100 = t115 + t154;
t82 = qJD(3) ^ 2;
t99 = pkin(6) * t82 + t145;
t97 = pkin(3) * t79 - t133;
t96 = -qJDD(4) - t105;
t95 = g(3) * t79 - t115;
t92 = t105 + t148;
t91 = -0.2e1 * pkin(2) * t117 - t131;
t90 = pkin(4) * t62 + qJDD(5) + t101;
t89 = -t99 + 0.2e1 * t78;
t6 = -t96 - t127;
t3 = -t146 * t111 + t90;
t8 = t93 * qJD(3) + t125;
t88 = qJD(2) * t8 + qJDD(2) * t19 - t145 + t3;
t87 = t23 * qJD(3) - t92;
t16 = t97 * qJD(3) - t125;
t4 = pkin(3) * t111 - t101;
t86 = -qJD(2) * t16 - qJDD(2) * t24 - t4 - t99;
t85 = -t79 * t119 + qJDD(4) + t92;
t5 = -pkin(6) * t111 + t100;
t84 = t5 * t80 + t6 * t79 + (t17 * t80 - t20 * t79) * qJD(3);
t43 = t59 * pkin(6);
t39 = qJ(5) * t111;
t38 = -t76 * t83 - t82;
t36 = t80 * t46;
t35 = g(1) * t143;
t32 = t59 * t133;
t30 = t58 * t133;
t28 = t139 * t79;
t27 = -qJDD(3) - t116;
t26 = qJDD(3) * t80 - t82 * t79;
t25 = qJDD(3) * t79 + t82 * t80;
t21 = t97 * qJD(2);
t15 = -t124 + t129;
t13 = -t139 * t128 - t123;
t11 = t93 * qJD(2);
t2 = -t80 * t119 + t39 + (-pkin(6) * t128 - t123) * qJD(2) + t100;
t1 = -t153 * qJ(5) - qJD(2) * t124 - t108 - t96;
t31 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t26, 0, t25, t5 * t79 - t6 * t80 - g(3) + (t17 * t79 + t20 * t80) * qJD(3), t26, t25, 0, -t1 * t80 + t2 * t79 - g(3) + (t7 * t79 + t144) * qJD(3); 0, qJDD(2), t113, t138, t76 * qJDD(2) + t79 * t102, -0.2e1 * t136 * t117 + 0.2e1 * t79 * t62, t25, t26, 0, t91 * t79 + t89 * t80 + t36, -t89 * t79 + t91 * t80 - t35, t147 * t79 + t86 * t80 + t36, t135 * qJDD(2) * pkin(6) - t138 + t84, -t147 * t80 + t86 * t79 + t35, t84 * pkin(6) - g(1) * t43 - g(2) * t106 + t18 * t16 + t4 * t24 - t94 * t46, -t28 * qJDD(3) + t36 + (-t107 * t79 - t15) * qJD(3) + t88 * t80, t29 * qJDD(3) + t35 + (t107 * t80 + t13) * qJD(3) + t88 * t79, (-qJD(3) * t7 - qJDD(2) * t29 - t2 + (-qJD(3) * t28 - t13) * qJD(2)) * t80 + (t126 - qJDD(2) * t28 - t1 + (-t15 + t129) * qJD(2)) * t79 + t138, t2 * t29 + t10 * t13 + t1 * t28 + t7 * t15 + t3 * t19 + t9 * t8 - g(1) * (-t59 * qJ(5) + t43) - g(2) * (pkin(4) * t141 + t106) + (-g(1) * (t94 - t69) + g(2) * qJ(5)) * t58; 0, 0, 0, 0, -t116, t136 * t83, t60, t62, qJDD(3), pkin(2) * t140 + t87, (t22 + t56) * qJD(3) + (pkin(2) * t83 + t138) * t80 + t95, 0.2e1 * t127 - qJDD(4) + (-t18 * t79 + t21 * t80) * qJD(2) + t87, -t97 * qJDD(2), -t22 * qJD(3) - t138 * t80 + (t18 * t80 + (t21 - t134) * t79) * qJD(2) - t95 + t151, t5 * qJ(4) - t6 * pkin(3) - t18 * t21 - t17 * t23 - g(1) * (-pkin(3) * t142 + t32) - g(2) * (-pkin(3) * t143 + t30) - g(3) * t137 + t152 * t20, t14 * qJD(3) + 0.2e1 * t108 + ((-t11 + t120) * t80 + t150) * qJD(2) - t85, -t12 * qJD(3) + t39 + (-g(3) + (-t11 - t134) * qJD(2)) * t79 + (-t132 * qJD(2) - t119 - t138) * t80 + t115 + t151, -t93 * qJDD(2), -g(1) * t32 - g(2) * t30 - g(3) * t114 + t2 * qJ(4) - t1 * t146 + t122 * t10 - t9 * t11 + t138 * t149 - t7 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t60, t38, -t20 * qJD(3) + t18 * t130 + t148 + t6, t27, t38, -t60, -t126 - t108 + (-t80 * t120 - t150) * qJD(2) + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 - 0.2e1 * t111, t60 + t102, -t135 * t83, (t144 + (t7 - t112) * t79) * qJD(2) + t90 + t113;];
tau_reg = t31;
