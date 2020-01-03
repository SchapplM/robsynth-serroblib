% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:23
% EndTime: 2019-12-31 18:57:27
% DurationCPUTime: 1.44s
% Computational Cost: add. (1430->266), mult. (2776->361), div. (0->0), fcn. (1631->6), ass. (0->144)
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t182 = -g(1) * t78 + g(2) * t81;
t77 = sin(qJ(3));
t79 = cos(qJ(4));
t158 = t77 * t79;
t80 = cos(qJ(3));
t53 = pkin(3) * t77 - pkin(7) * t80 + qJ(2);
t76 = sin(qJ(4));
t82 = -pkin(1) - pkin(6);
t148 = t82 * t158 + t76 * t53;
t141 = qJD(1) * t77;
t62 = qJD(4) + t141;
t30 = t53 * qJD(1);
t60 = t82 * qJD(1) + qJD(2);
t52 = t77 * t60;
t37 = qJD(3) * pkin(7) + t52;
t11 = t79 * t30 - t37 * t76;
t137 = qJD(3) * t76;
t140 = qJD(1) * t80;
t48 = t140 * t79 + t137;
t7 = -qJ(5) * t48 + t11;
t6 = pkin(4) * t62 + t7;
t12 = t30 * t76 + t37 * t79;
t131 = t79 * qJD(3);
t46 = t140 * t76 - t131;
t8 = -qJ(5) * t46 + t12;
t104 = t6 * t79 + t76 * t8;
t183 = -qJD(1) * t104 + t182;
t115 = pkin(4) * t76 - t82;
t153 = t79 * t81;
t161 = t76 * t78;
t31 = -t161 * t77 + t153;
t156 = t78 * t79;
t157 = t77 * t81;
t33 = t157 * t76 + t156;
t181 = -g(1) * t31 - g(2) * t33;
t172 = g(3) * t77;
t136 = qJD(3) * t77;
t110 = -qJDD(3) * pkin(3) + t60 * t136;
t58 = t82 * qJDD(1) + qJDD(2);
t165 = t58 * t80;
t22 = t110 - t165;
t180 = qJD(4) * pkin(7) * t62 - t182 * t80 - t172 + t22;
t120 = t76 * t141;
t127 = t80 * qJDD(1);
t14 = -qJD(3) * t120 + qJD(4) * t48 - t79 * qJDD(3) + t127 * t76;
t119 = t77 * t131;
t133 = qJD(4) * t76;
t92 = t133 * t80 + t119;
t13 = qJD(1) * t92 - qJD(4) * t131 - t76 * qJDD(3) - t79 * t127;
t179 = t48 ^ 2;
t178 = -t7 + t6;
t174 = g(1) * t81;
t171 = g(3) * t80;
t151 = qJ(5) + pkin(7);
t111 = qJD(4) * t151;
t160 = t76 * t80;
t107 = pkin(3) * t80 + pkin(7) * t77;
t51 = t107 * qJD(1);
t36 = t79 * t51;
t170 = -t76 * qJD(5) - t111 * t79 + t60 * t160 - t36 - (pkin(4) * t80 + qJ(5) * t158) * qJD(1);
t169 = t13 * t76;
t168 = t46 * t62;
t167 = t48 * t62;
t166 = t48 * t79;
t164 = t60 * t80;
t163 = t62 * t76;
t126 = qJD(1) * qJD(3);
t113 = t80 * t126;
t128 = t77 * qJDD(1);
t44 = qJDD(4) + t113 + t128;
t162 = t76 * t44;
t159 = t76 * t81;
t155 = t79 * t44;
t154 = t79 * t80;
t152 = t80 * t13;
t130 = t79 * qJD(5);
t149 = t60 * t154 + t76 * t51;
t150 = -qJ(5) * t120 - t111 * t76 + t130 - t149;
t147 = t81 * pkin(1) + t78 * qJ(2);
t74 = t80 ^ 2;
t146 = t77 ^ 2 - t74;
t83 = qJD(3) ^ 2;
t84 = qJD(1) ^ 2;
t145 = -t83 - t84;
t144 = qJ(2) * t84;
t143 = qJ(5) * t80;
t142 = pkin(1) * qJDD(1);
t139 = qJD(3) * t46;
t138 = qJD(3) * t48;
t135 = qJD(3) * t80;
t134 = qJD(3) * t82;
t132 = qJD(4) * t79;
t129 = qJDD(3) * t77;
t125 = qJDD(1) * qJ(2);
t45 = qJD(3) * t107 + qJD(2);
t18 = qJD(1) * t45 + qJDD(1) * t53;
t23 = qJDD(3) * pkin(7) + t135 * t60 + t58 * t77;
t124 = -t30 * t132 - t76 * t18 - t79 * t23;
t118 = t80 * t134;
t123 = t79 * t118 + t53 * t132 + t76 * t45;
t121 = t79 * t143;
t117 = 0.2e1 * qJD(1) * qJD(2);
t116 = -t76 * t82 + pkin(4);
t112 = t62 * t82 + t37;
t109 = -qJD(4) * t30 - t23;
t108 = qJD(4) * t77 + qJD(1);
t106 = g(2) * t78 + t174;
t103 = t6 * t76 - t79 * t8;
t102 = qJDD(2) + t182;
t64 = pkin(4) * t79 + pkin(3);
t100 = -t151 * t80 + t64 * t77;
t38 = -qJD(3) * pkin(3) - t164;
t99 = -t182 - t58;
t97 = t62 * t132 + t162;
t96 = -t62 * t133 + t155;
t95 = -t133 * t37 - t124;
t93 = 0.2e1 * qJ(2) * t126 + qJDD(3) * t82;
t91 = pkin(4) * t14 + qJDD(5) + t110;
t90 = t99 + t144;
t89 = -pkin(7) * t44 + t62 * t38;
t87 = -t106 + t117 + 0.2e1 * t125;
t16 = t79 * t18;
t1 = pkin(4) * t44 + qJ(5) * t13 - qJD(4) * t12 - qJD(5) * t48 - t76 * t23 + t16;
t2 = -qJ(5) * t14 - qJD(5) * t46 + t95;
t86 = -qJD(4) * t104 - t1 * t76 + t2 * t79;
t85 = -t82 * t83 + t87;
t70 = t81 * qJ(2);
t67 = qJDD(3) * t80;
t56 = t151 * t79;
t55 = t151 * t76;
t43 = t46 ^ 2;
t42 = t79 * t53;
t34 = t153 * t77 - t161;
t32 = t156 * t77 + t159;
t29 = t79 * t45;
t21 = -t143 * t76 + t148;
t20 = pkin(4) * t46 + qJD(5) + t38;
t17 = t116 * t77 - t121 + t42;
t5 = t91 - t165;
t4 = -qJD(4) * t121 + (-qJD(5) * t80 + (qJ(5) * qJD(3) - qJD(4) * t82) * t77) * t76 + t123;
t3 = qJ(5) * t119 + t29 - t148 * qJD(4) + (qJ(5) * t133 + qJD(3) * t116 - t130) * t80;
t9 = [qJDD(1), -t182, t106, t102 - 0.2e1 * t142, t87, -(qJDD(2) - t142) * pkin(1) - g(1) * (-t78 * pkin(1) + t70) - g(2) * t147 + (t117 + t125) * qJ(2), qJDD(1) * t74 - 0.2e1 * t113 * t77, 0.2e1 * t126 * t146 - 0.2e1 * t127 * t77, -t77 * t83 + t67, -t80 * t83 - t129, 0, t77 * t85 + t80 * t93, -t77 * t93 + t80 * t85, -t152 * t79 - t48 * t92, (t46 * t79 + t48 * t76) * t136 + (t169 - t14 * t79 + (t46 * t76 - t166) * qJD(4)) * t80, (-t62 * t131 - t13) * t77 + (t96 + t138) * t80, (t62 * t137 - t14) * t77 + (-t97 - t139) * t80, t62 * t135 + t44 * t77, -g(1) * t34 - g(2) * t32 + t29 * t62 + t42 * t44 + (-t112 * t132 + t134 * t46 + t16) * t77 + (qJD(3) * t11 + t132 * t38 - t14 * t82) * t80 + ((-qJD(4) * t53 - t118) * t62 + t22 * t80 + (-qJD(3) * t38 - t44 * t82 + t109) * t77) * t76, -t123 * t62 - t148 * t44 + g(1) * t33 - g(2) * t31 + (t112 * t133 + (-t38 * t79 + t48 * t82) * qJD(3) + t124) * t77 + (-qJD(3) * t12 + t13 * t82 - t133 * t38 + t22 * t79) * t80, t13 * t17 - t14 * t21 - t3 * t48 - t4 * t46 + t104 * t136 + (qJD(4) * t103 - t1 * t79 - t2 * t76 + t106) * t80, t2 * t21 + t8 * t4 + t1 * t17 + t6 * t3 - g(1) * (t157 * t64 + t70) - g(2) * (pkin(4) * t159 + t81 * pkin(6) + t147) - t20 * t115 * t136 + (pkin(4) * t132 * t20 + t115 * t5 + t151 * t174) * t80 + (g(1) * t115 - g(2) * t100) * t78; 0, 0, 0, qJDD(1), -t84, t102 - t142 - t144, 0, 0, 0, 0, 0, t145 * t77 + t67, t145 * t80 - t129, 0, 0, 0, 0, 0, -t80 * t14 + (t139 - t162) * t77 + (-t108 * t79 - t135 * t76) * t62, t152 + (t138 - t155) * t77 + (t108 * t76 - t131 * t80) * t62, (t108 * t48 - t135 * t46 - t14 * t77) * t79 + (t108 * t46 - t13 * t77 + t135 * t48) * t76, (-qJD(3) * t103 - t5) * t80 + (qJD(3) * t20 + t86) * t77 + t183; 0, 0, 0, 0, 0, 0, t80 * t84 * t77, -t146 * t84, t127, -t128, qJDD(3), -t80 * t90 + t172, t77 * t90 + t171, t62 * t166 - t169, (-t13 - t168) * t79 + (-t14 - t167) * t76, (t62 * t158 - t48 * t80) * qJD(1) + t97, (-t77 * t163 + t46 * t80) * qJD(1) + t96, -t62 * t140, -t11 * t140 - t46 * t52 - pkin(3) * t14 - t36 * t62 + (t62 * t164 + t89) * t76 - t180 * t79, pkin(3) * t13 + t12 * t140 + t149 * t62 + t180 * t76 - t48 * t52 + t89 * t79, -t13 * t55 - t14 * t56 - t150 * t46 - t170 * t48 + t183 * t77 - t171 + t86, t2 * t56 - t1 * t55 - t5 * t64 + g(3) * t100 + t150 * t8 + t170 * t6 + (pkin(4) * t163 - t52) * t20 + t182 * (t151 * t77 + t64 * t80); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t46, -t43 + t179, -t13 + t168, -t14 + t167, t44, -t37 * t132 + t12 * t62 - t38 * t48 + t16 + (t109 + t171) * t76 + t181, g(1) * t32 - g(2) * t34 + g(3) * t154 + t11 * t62 + t38 * t46 - t95, pkin(4) * t13 - t178 * t46, t178 * t8 + (g(3) * t160 - t20 * t48 + t1 + t181) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 - t179, t46 * t8 + t48 * t6 + t80 * t99 - t172 + t91;];
tau_reg = t9;
