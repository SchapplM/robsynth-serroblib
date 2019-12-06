% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:44
% EndTime: 2019-12-05 15:54:49
% DurationCPUTime: 1.42s
% Computational Cost: add. (1247->218), mult. (2856->302), div. (0->0), fcn. (2313->14), ass. (0->128)
t100 = cos(qJ(4));
t92 = sin(pkin(9));
t94 = cos(pkin(9));
t97 = sin(qJ(4));
t66 = t100 * t92 + t94 * t97;
t172 = t66 * qJD(2);
t158 = t100 * t94;
t141 = qJD(2) * t158;
t163 = t92 * t97;
t142 = qJD(2) * t163;
t55 = -t141 + t142;
t96 = sin(qJ(5));
t99 = cos(qJ(5));
t20 = t172 * t96 + t99 * t55;
t91 = qJD(4) + qJD(5);
t165 = t20 * t91;
t152 = qJD(5) * t99;
t153 = qJD(5) * t96;
t144 = qJDD(2) * t100;
t147 = t94 * qJDD(2);
t143 = qJD(4) * t141 + t92 * t144 + t97 * t147;
t27 = -qJD(4) * t142 + t143;
t148 = t92 * qJDD(2);
t128 = -t94 * t144 + t97 * t148;
t59 = t66 * qJD(4);
t28 = qJD(2) * t59 + t128;
t4 = -t55 * t152 - t153 * t172 + t99 * t27 - t96 * t28;
t184 = t4 + t165;
t124 = -t172 * t99 + t55 * t96;
t107 = t124 * qJD(5) - t96 * t27 - t99 * t28;
t166 = t124 * t91;
t183 = t107 - t166;
t101 = cos(qJ(2));
t149 = t101 * qJD(1);
t132 = qJD(3) - t149;
t182 = t124 * t20;
t168 = g(3) * t101;
t93 = sin(pkin(8));
t95 = cos(pkin(8));
t131 = g(1) * t95 + g(2) * t93;
t98 = sin(qJ(2));
t177 = t131 * t98;
t110 = t177 - t168;
t150 = qJDD(1) - g(3);
t181 = t150 * t101 + t177;
t180 = t124 ^ 2 - t20 ^ 2;
t155 = qJD(1) * t98;
t73 = qJD(2) * qJ(3) + t155;
t135 = pkin(6) * qJD(2) + t73;
t48 = t135 * t92;
t49 = t135 * t94;
t121 = -t100 * t49 + t48 * t97;
t13 = -pkin(7) * t55 - t121;
t156 = t101 * t95;
t157 = t101 * t93;
t169 = g(3) * t98;
t50 = qJDD(2) * qJ(3) + t98 * qJDD(1) + (qJD(3) + t149) * qJD(2);
t134 = pkin(6) * qJDD(2) + t50;
t33 = t134 * t92;
t34 = t134 * t94;
t138 = -t100 * t33 - t97 * t34;
t2 = qJDD(4) * pkin(4) - pkin(7) * t27 + t121 * qJD(4) + t138;
t82 = -pkin(3) * t94 - pkin(2);
t60 = t82 * qJD(2) + t132;
t31 = t55 * pkin(4) + t60;
t90 = pkin(9) + qJ(4);
t86 = qJ(5) + t90;
t80 = sin(t86);
t81 = cos(t86);
t179 = t31 * t20 + t81 * t169 + t13 * t153 + (-t13 * t91 - t2) * t96 - g(1) * (-t81 * t156 - t80 * t93) - g(2) * (-t81 * t157 + t80 * t95);
t167 = t13 * t99;
t176 = -t100 * t48 - t49 * t97;
t12 = -pkin(7) * t172 + t176;
t9 = qJD(4) * pkin(4) + t12;
t127 = -t9 * t96 - t167;
t122 = t100 * t34 - t97 * t33;
t3 = -t28 * pkin(7) + t176 * qJD(4) + t122;
t178 = t31 * t124 - g(1) * (-t80 * t156 + t81 * t93) - g(2) * (-t80 * t157 - t81 * t95) + t127 * qJD(5) + t99 * t2 - t96 * t3 + t80 * t169;
t162 = pkin(6) + qJ(3);
t69 = t162 * t92;
t70 = t162 * t94;
t160 = t100 * t70 - t97 * t69;
t175 = -t160 * qJD(4) - t132 * t66;
t174 = t131 * t101;
t120 = t158 - t163;
t58 = t120 * qJD(4);
t145 = qJD(1) * qJD(2);
t151 = qJDD(2) * pkin(2);
t146 = t98 * t145 + qJDD(3);
t119 = -t101 * qJDD(1) + t146;
t54 = t119 - t151;
t171 = (t131 + t145) * t98 + t151 - t54 - t168;
t113 = t120 * t101;
t63 = t100 * t69;
t170 = qJD(1) * t113 + (qJD(3) * t92 + qJD(4) * t70) * t97 - qJD(3) * t158 + qJD(4) * t63;
t159 = t92 ^ 2 + t94 ^ 2;
t154 = qJD(2) * t98;
t140 = t159 * t50;
t137 = -t70 * t97 - t63;
t136 = t101 * t159;
t133 = t159 * qJDD(2);
t15 = -pkin(7) * t66 + t137;
t130 = t59 * pkin(7) - qJD(5) * t15 + t170;
t16 = pkin(7) * t120 + t160;
t129 = pkin(7) * t58 + qJD(5) * t16 - t175;
t51 = t66 * t98;
t52 = t120 * t98;
t126 = -t51 * t99 - t52 * t96;
t125 = -t51 * t96 + t52 * t99;
t29 = -t120 * t99 + t66 * t96;
t30 = t120 * t96 + t66 * t99;
t123 = pkin(4) * t59 - t155;
t102 = qJD(2) ^ 2;
t117 = qJDD(2) * t101 - t102 * t98;
t108 = t132 * t159;
t39 = t82 * qJDD(2) + t119;
t105 = t140 - t169 - t174;
t87 = qJDD(4) + qJDD(5);
t85 = cos(t90);
t84 = sin(t90);
t71 = -qJD(2) * pkin(2) + t132;
t43 = -pkin(4) * t120 + t82;
t18 = -t101 * t172 - t98 * t58;
t17 = qJD(2) * t113 - qJD(4) * t51;
t14 = pkin(4) * t28 + t39;
t7 = t30 * qJD(5) + t96 * t58 + t99 * t59;
t6 = -t29 * qJD(5) + t99 * t58 - t96 * t59;
t1 = [t150, 0, t117, -qJDD(2) * t98 - t101 * t102, t117 * t94, -t117 * t92, t102 * t136 + t98 * t133, -t101 * t54 - g(3) + t98 * t140 + (t73 * t136 + t71 * t98) * qJD(2), 0, 0, 0, 0, 0, qJD(4) * t18 - qJDD(4) * t51 - t101 * t28 + t55 * t154, -qJD(4) * t17 - qJDD(4) * t52 - t101 * t27 + t154 * t172, 0, 0, 0, 0, 0, (-qJD(5) * t125 - t96 * t17 + t99 * t18) * t91 + t126 * t87 + t20 * t154 + t101 * t107, -(qJD(5) * t126 + t99 * t17 + t96 * t18) * t91 - t125 * t87 - t124 * t154 - t101 * t4; 0, qJDD(2), t181, -t150 * t98 + t174, t171 * t94, -t171 * t92, qJ(3) * t133 + t108 * qJD(2) + t105, -t71 * t155 + (-t54 + t110) * pkin(2) + t105 * qJ(3) + t108 * t73, t172 * t58 + t27 * t66, t120 * t27 - t172 * t59 - t28 * t66 - t55 * t58, qJD(4) * t58 + qJDD(4) * t66, -qJD(4) * t59 + qJDD(4) * t120, 0, t175 * qJD(4) + t137 * qJDD(4) + t110 * t85 - t120 * t39 - t55 * t155 + t82 * t28 + t60 * t59, t170 * qJD(4) - t160 * qJDD(4) - t110 * t84 - t155 * t172 + t82 * t27 + t39 * t66 + t60 * t58, -t124 * t6 + t30 * t4, t107 * t30 + t124 * t7 - t20 * t6 - t29 * t4, t30 * t87 + t6 * t91, -t29 * t87 - t7 * t91, 0, (t15 * t99 - t16 * t96) * t87 - t43 * t107 + t14 * t29 + t31 * t7 + (-t129 * t99 + t130 * t96) * t91 + t123 * t20 + t110 * t81, -(t15 * t96 + t16 * t99) * t87 + t43 * t4 + t14 * t30 + t31 * t6 + (t129 * t96 + t130 * t99) * t91 - t123 * t124 - t110 * t80; 0, 0, 0, 0, -t147, t148, -t159 * t102, -t159 * t73 * qJD(2) + t146 - t151 - t181, 0, 0, 0, 0, 0, 0.2e1 * qJD(4) * t172 + t128, (-t55 - t142) * qJD(4) + t143, 0, 0, 0, 0, 0, -t107 - t166, t4 - t165; 0, 0, 0, 0, 0, 0, 0, 0, t172 * t55, t172 ^ 2 - t55 ^ 2, (t55 - t142) * qJD(4) + t143, -t128, qJDD(4), -t60 * t172 - g(1) * (-t84 * t156 + t85 * t93) - g(2) * (-t84 * t157 - t85 * t95) + t84 * t169 + t138, t60 * t55 - g(1) * (-t85 * t156 - t84 * t93) - g(2) * (-t85 * t157 + t84 * t95) + t85 * t169 - t122, -t182, t180, t184, t183, t87, -(-t12 * t96 - t167) * t91 + (-t91 * t153 - t172 * t20 + t87 * t99) * pkin(4) + t178, (-qJD(5) * t9 + t12 * t91 - t3) * t99 + (t124 * t172 - t91 * t152 - t87 * t96) * pkin(4) + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t180, t184, t183, t87, -t127 * t91 + t178, (-t3 + (-qJD(5) + t91) * t9) * t99 + t179;];
tau_reg = t1;
