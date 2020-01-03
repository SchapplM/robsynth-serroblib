% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR9
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:50
% EndTime: 2019-12-31 18:02:53
% DurationCPUTime: 0.97s
% Computational Cost: add. (963->230), mult. (1796->310), div. (0->0), fcn. (1173->8), ass. (0->129)
t115 = qJD(1) * qJD(4);
t73 = sin(qJ(4));
t108 = t73 * t115;
t75 = cos(qJ(4));
t59 = t75 * qJDD(1);
t157 = -t108 + t59;
t121 = qJ(2) * qJD(1);
t76 = -pkin(1) - pkin(2);
t42 = t76 * qJD(1) + qJD(2);
t69 = sin(pkin(8));
t70 = cos(pkin(8));
t24 = t70 * t121 + t69 * t42;
t21 = -qJD(1) * pkin(6) + t24;
t15 = qJD(3) * t73 + t21 * t75;
t132 = qJD(4) * t15;
t117 = qJ(2) * qJDD(1);
t41 = t76 * qJDD(1) + qJDD(2);
t141 = t70 * t117 + t69 * t41;
t116 = qJD(1) * qJD(2);
t47 = t70 * t116;
t19 = t47 + t141;
t17 = -qJDD(1) * pkin(6) + t19;
t2 = -qJDD(4) * pkin(4) - t75 * qJDD(3) + t73 * t17 + t132;
t43 = qJD(1) * t75 + qJD(5);
t150 = sin(qJ(1));
t151 = cos(qJ(1));
t27 = -t150 * t69 - t151 * t70;
t155 = g(1) * t27;
t28 = -t150 * t70 + t151 * t69;
t94 = g(2) * t28 + t155;
t96 = -pkin(4) * t73 + pkin(7) * t75;
t156 = (pkin(7) * qJD(5) + t96 * qJD(1)) * t43 + t73 * t94 - g(3) * t75 + t2;
t154 = g(1) * t28;
t153 = g(2) * t27;
t72 = sin(qJ(5));
t125 = qJD(5) * t72;
t111 = t73 * t125;
t114 = qJD(4) * qJD(5);
t74 = cos(qJ(5));
t107 = t75 * t115;
t118 = t73 * qJDD(1);
t88 = t107 + t118;
t7 = qJD(1) * t111 + t72 * qJDD(4) + (t114 - t88) * t74;
t152 = t7 * t72;
t149 = t28 * t75;
t127 = qJD(4) * t74;
t135 = qJD(1) * t73;
t29 = t72 * t135 + t127;
t148 = t29 * t43;
t123 = t72 * qJD(4);
t30 = t74 * t135 - t123;
t147 = t30 * t43;
t146 = t43 * t74;
t145 = t70 * t73;
t26 = -qJDD(5) - t157;
t144 = t72 * t26;
t143 = t72 * t75;
t142 = t74 * t75;
t36 = t70 * qJ(2) + t69 * t76;
t140 = t151 * pkin(1) + t150 * qJ(2);
t139 = g(1) * t150 - g(2) * t151;
t66 = t73 ^ 2;
t138 = -t75 ^ 2 + t66;
t77 = qJD(4) ^ 2;
t78 = qJD(1) ^ 2;
t137 = t77 + t78;
t136 = pkin(1) * qJDD(1);
t134 = qJD(2) * t70;
t14 = qJD(3) * t75 - t21 * t73;
t133 = qJD(4) * t14;
t131 = qJD(4) * t29;
t130 = qJD(4) * t30;
t32 = -pkin(6) + t36;
t129 = qJD(4) * t32;
t128 = qJD(4) * t73;
t126 = qJD(4) * t75;
t124 = qJD(5) * t74;
t122 = qJDD(3) + g(3);
t120 = qJDD(4) * t73;
t119 = qJDD(4) * t75;
t113 = t43 * t123;
t112 = t43 * t127;
t110 = 0.2e1 * t116;
t45 = t69 * t116;
t12 = qJD(4) * pkin(7) + t15;
t106 = t32 * t43 + t12;
t105 = -t69 * t117 + t41 * t70;
t23 = -t69 * t121 + t70 * t42;
t35 = -t69 * qJ(2) + t70 * t76;
t20 = qJD(1) * pkin(3) - t23;
t97 = pkin(4) * t75 + pkin(7) * t73;
t13 = t97 * qJD(1) + t20;
t103 = -qJDD(4) * pkin(7) - qJD(5) * t13 - t73 * qJDD(3) - t75 * t17 - t133;
t102 = 0.2e1 * t107;
t101 = qJDD(2) - t136;
t31 = pkin(3) - t35;
t18 = t105 - t45;
t100 = t7 + t112;
t8 = t74 * qJDD(4) + (t75 * t123 + t73 * t124) * qJD(1) + (t118 - t114) * t72;
t99 = -t8 + t113;
t98 = -t150 * pkin(1) + t151 * qJ(2);
t95 = t153 - t154;
t93 = t23 * t69 - t24 * t70;
t92 = qJD(4) * t21 - t122;
t91 = -qJD(5) * t12 + t153;
t16 = qJDD(1) * pkin(3) - t18;
t90 = t43 * t124 - t144;
t89 = t43 * t125 + t74 * t26;
t87 = g(1) * t151 + g(2) * t150;
t86 = -t90 - t131;
t85 = t89 - t130;
t84 = -g(2) * t149 - g(3) * t73 + t103;
t11 = -qJD(4) * pkin(4) - t14;
t83 = pkin(7) * t26 + (t11 + t14) * t43;
t82 = qJD(1) * t20 - qJD(3) * qJD(4) - t17 - t94;
t81 = -qJDD(4) * t32 + (-qJD(1) * t31 - t134 - t20) * qJD(4);
t80 = -qJD(4) * t11 - t43 * t134 + t26 * t32 + t103;
t79 = qJDD(1) * t31 - t32 * t77 + t16 + t45 + t95;
t39 = -t73 * t77 + t119;
t38 = -t75 * t77 - t120;
t25 = t69 * qJD(2) + t96 * qJD(4);
t22 = t31 + t97;
t10 = -t27 * t142 + t28 * t72;
t9 = t27 * t143 + t28 * t74;
t6 = t157 * pkin(4) + t88 * pkin(7) + t16;
t5 = t74 * t6;
t4 = t12 * t74 + t13 * t72;
t3 = -t12 * t72 + t13 * t74;
t1 = [qJDD(1), t139, t87, -qJDD(2) + 0.2e1 * t136 + t139, t110 - t87 + 0.2e1 * t117, -t101 * pkin(1) - g(1) * t98 - g(2) * t140 + (t110 + t117) * qJ(2), -qJDD(1) * t35 - t105 + 0.2e1 * t45 + t95, qJDD(1) * t36 + t141 + 0.2e1 * t47 + t94, t19 * t36 + t18 * t35 - g(1) * (-t150 * pkin(2) + t98) - g(2) * (t151 * pkin(2) + t140) - t93 * qJD(2), qJDD(1) * t66 + t73 * t102, -0.2e1 * t138 * t115 + 0.2e1 * t73 * t59, t38, -t39, 0, t81 * t73 + t79 * t75, -t79 * t73 + t81 * t75, -t7 * t74 * t73 + (t126 * t74 - t111) * t30, (-t29 * t74 - t30 * t72) * t126 + (t152 - t74 * t8 + (t29 * t72 - t30 * t74) * qJD(5)) * t73, (t7 - t112) * t75 + (t89 + t130) * t73, (t8 + t113) * t75 + (t90 - t131) * t73, -t128 * t43 - t26 * t75, -g(2) * t10 + (-t29 * t129 + t5) * t75 + (-qJD(4) * t3 - t29 * t134 - t32 * t8) * t73 + (-g(1) * t149 - t22 * t26 + t25 * t43 + (-t106 * t75 - t11 * t73) * qJD(5)) * t74 + ((-qJD(5) * t22 + t128 * t32) * t43 - t2 * t73 - t155 + t80 * t75) * t72, -(t124 * t22 + t72 * t25) * t43 + t22 * t144 - t74 * t155 - g(2) * t9 + (-t30 * t129 + (qJD(5) * t106 + t154 - t6) * t72 + t80 * t74) * t75 + (-t30 * t134 + t11 * t125 - t2 * t74 + t32 * t7 + (t32 * t146 + t4) * qJD(4)) * t73; 0, 0, 0, -qJDD(1), -t78, -qJ(2) * t78 + t101 - t139, -qJDD(1) * t70 - t69 * t78, qJDD(1) * t69 - t70 * t78, t93 * qJD(1) + t18 * t70 + t19 * t69 - t139, 0, 0, 0, 0, 0, (0.2e1 * t108 - t59) * t70 + (-t137 * t75 - t120) * t69, (t102 + t118) * t70 + (t137 * t73 - t119) * t69, 0, 0, 0, 0, 0, t89 * t70 + (t73 * t99 + t75 * t86) * t69 + (-(-t70 * t143 + t69 * t74) * t43 + t29 * t145) * qJD(1), t90 * t70 + (t100 * t73 + t75 * t85) * t69 + ((t70 * t142 + t69 * t72) * t43 + t30 * t145) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, t122, 0, 0, 0, 0, 0, t39, t38, 0, 0, 0, 0, 0, t73 * t86 - t75 * t99, -t100 * t75 + t73 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 * t78 * t75, t138 * t78, -t118, -t59, qJDD(4), t82 * t73 - t92 * t75 + t132, t92 * t73 + t75 * t82 + t133, -t30 * t146 + t152, (t7 + t148) * t74 + (t8 + t147) * t72, (t43 * t142 - t30 * t73) * qJD(1) + t90, (-t43 * t143 + t29 * t73) * qJD(1) - t89, t43 * t135, pkin(4) * t8 + t3 * t135 + t15 * t29 - t156 * t74 + t83 * t72, -pkin(4) * t7 - t4 * t135 + t15 * t30 + t156 * t72 + t83 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t29, -t29 ^ 2 + t30 ^ 2, t7 - t148, t8 - t147, -t26, -g(1) * t9 + t11 * t30 + t4 * t43 + t72 * t84 + t74 * t91 + t5, g(1) * t10 - t11 * t29 + t3 * t43 + (-t6 - t91) * t72 + t84 * t74;];
tau_reg = t1;
