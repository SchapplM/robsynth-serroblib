% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:45:01
% EndTime: 2019-03-09 01:45:04
% DurationCPUTime: 1.76s
% Computational Cost: add. (2180->269), mult. (4094->346), div. (0->0), fcn. (2765->14), ass. (0->153)
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t187 = cos(pkin(10));
t159 = t187 * t118;
t186 = sin(pkin(10));
t67 = -t186 * t115 + t159;
t111 = sin(pkin(9));
t81 = pkin(1) * t111 + qJ(3);
t188 = qJDD(1) * t81;
t157 = t186 * t118;
t127 = t115 * t187 + t157;
t57 = t127 * qJD(1);
t209 = qJD(6) + t57;
t114 = sin(qJ(6));
t117 = cos(qJ(6));
t175 = t117 * qJD(4);
t59 = t67 * qJD(1);
t44 = t114 * t59 - t175;
t212 = t209 * t44;
t105 = qJ(1) + pkin(9);
t96 = sin(t105);
t98 = cos(t105);
t165 = -g(1) * t96 + g(2) * t98;
t152 = t117 * t209;
t149 = qJDD(1) * t186;
t150 = qJDD(1) * t187;
t36 = -t59 * qJD(4) - t115 * t150 - t118 * t149;
t34 = -qJDD(6) + t36;
t194 = t114 * t34;
t211 = t152 * t209 - t194;
t70 = qJD(1) * t81;
t130 = -t70 * qJD(1) + t165;
t112 = cos(pkin(9));
t89 = -pkin(1) * t112 - pkin(2);
t77 = -pkin(7) + t89;
t189 = qJ(5) - t77;
t62 = t189 * t115;
t208 = qJDD(1) * t89;
t64 = qJD(1) * t77 + qJD(3);
t151 = qJ(5) * qJD(1) - t64;
t207 = -qJD(2) * t118 + t115 * t151;
t61 = t127 * qJD(4);
t138 = t209 * t61 + t34 * t67;
t178 = qJD(6) * t114;
t166 = t67 * t178;
t206 = -t117 * t138 - t166 * t209;
t205 = 0.2e1 * qJD(4) * t70 + qJDD(4) * t77;
t134 = -qJ(5) * qJDD(1) - qJD(1) * qJD(5);
t63 = qJDD(1) * t77 + qJDD(3);
t53 = t118 * t63;
t14 = qJDD(4) * pkin(4) + qJD(4) * t207 - t115 * qJDD(2) + t134 * t118 + t53;
t19 = (-qJD(4) * t151 + qJDD(2)) * t118 + (-qJD(2) * qJD(4) + t134 + t63) * t115;
t3 = t14 * t187 - t186 * t19;
t1 = -qJDD(4) * pkin(5) - t3;
t163 = t186 * t207;
t180 = qJD(1) * t118;
t55 = t118 * t64;
t41 = -qJ(5) * t180 - t115 * qJD(2) + t55;
t40 = qJD(4) * pkin(4) + t41;
t17 = t187 * t40 + t163;
t12 = -qJD(4) * pkin(5) - t17;
t102 = t115 * pkin(4);
t56 = qJD(1) * t102 + qJD(5) + t70;
t24 = pkin(5) * t57 - pkin(8) * t59 + t56;
t4 = t186 * t14 + t187 * t19;
t160 = qJDD(4) * pkin(8) + qJD(6) * t24 + t4;
t128 = qJD(4) * t62 - qJD(5) * t118;
t174 = t118 * qJD(4);
t47 = -qJD(5) * t115 - t174 * t189;
t23 = t128 * t186 + t187 * t47;
t144 = t102 + t81;
t28 = pkin(5) * t127 - pkin(8) * t67 + t144;
t31 = -t157 * t189 - t187 * t62;
t204 = t1 * t67 - t12 * t61 - (qJD(6) * t28 + t23) * t209 - t160 * t127 + t31 * t34;
t87 = pkin(4) * t186 + pkin(8);
t104 = qJ(4) + pkin(10);
t95 = sin(t104);
t97 = cos(t104);
t203 = -t165 * t97 + (pkin(4) * t180 + pkin(5) * t59 + pkin(8) * t57 + qJD(6) * t87) * t209 - g(3) * t95 + t1;
t37 = qJD(1) * t61 + t115 * t149 - t118 * t150;
t15 = qJD(6) * t175 + t114 * qJDD(4) - t117 * t37 - t178 * t59;
t46 = qJD(4) * t114 + t117 * t59;
t58 = t67 * qJD(4);
t201 = t127 * t15 + t46 * t58;
t199 = t12 * t67;
t198 = t28 * t34;
t197 = t44 * t59;
t196 = t46 * t59;
t38 = t187 * t207;
t18 = t186 * t40 - t38;
t195 = t114 * t15;
t193 = t114 * t96;
t192 = t114 * t98;
t29 = t117 * t34;
t191 = t117 * t96;
t190 = t117 * t98;
t185 = qJD(1) * t56;
t184 = pkin(4) * t174 + qJD(3);
t110 = qJDD(2) - g(3);
t109 = t118 ^ 2;
t182 = t115 ^ 2 - t109;
t120 = qJD(4) ^ 2;
t121 = qJD(1) ^ 2;
t181 = -t120 - t121;
t177 = qJD(6) * t117;
t173 = qJD(1) * qJD(4);
t107 = qJD(3) * qJD(1);
t172 = qJDD(4) * t115;
t171 = qJDD(4) * t118;
t170 = t115 * qJDD(1);
t169 = t118 * qJDD(1);
t119 = cos(qJ(1));
t168 = t119 * pkin(1) + t98 * pkin(2) + t96 * qJ(3);
t65 = t107 + t188;
t116 = sin(qJ(1));
t164 = -pkin(1) * t116 + t98 * qJ(3);
t162 = t118 * t173;
t13 = qJD(4) * pkin(8) + t18;
t133 = qJDD(5) + t65 + (t162 + t170) * pkin(4);
t8 = -pkin(5) * t36 + pkin(8) * t37 + t133;
t161 = qJD(6) * t13 - t8;
t153 = -t117 * qJDD(4) - t114 * t37;
t148 = -qJD(6) * t127 - qJD(1);
t146 = -g(1) * t98 - g(2) * t96;
t141 = g(1) * t116 - g(2) * t119;
t140 = t15 * t67 - t46 * t61;
t16 = qJD(6) * t46 + t153;
t139 = t127 * t16 + t44 * t58;
t136 = -t29 + (-t114 * t57 - t178) * t209;
t135 = g(3) * t97 - t160;
t132 = qJDD(3) + t208;
t129 = t146 + t188;
t21 = t187 * t41 + t163;
t126 = t87 * t34 + (t12 + t21) * t209;
t125 = -t177 * t209 * t67 + t114 * t138;
t124 = t127 * t4 - t17 * t61 + t18 * t58 + t3 * t67 + t165;
t123 = -t120 * t77 + t107 + t129 + t65;
t113 = -qJ(5) - pkin(7);
t88 = -pkin(4) * t187 - pkin(5);
t69 = -t115 * t120 + t171;
t68 = -t118 * t120 - t172;
t51 = t95 * t190 - t193;
t50 = t95 * t192 + t191;
t49 = t95 * t191 + t192;
t48 = -t95 * t193 + t190;
t30 = t159 * t189 - t186 * t62;
t25 = pkin(5) * t58 + pkin(8) * t61 + t184;
t22 = -t128 * t187 + t186 * t47;
t20 = t186 * t41 - t38;
t7 = t117 * t8;
t6 = t114 * t24 + t117 * t13;
t5 = -t114 * t13 + t117 * t24;
t2 = [qJDD(1), t141, g(1) * t119 + g(2) * t116 (t141 + (t111 ^ 2 + t112 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(3) + t165 + 0.2e1 * t208, 0.2e1 * t107 + t129 + t188, t65 * t81 + t70 * qJD(3) + t132 * t89 - g(1) * (-pkin(2) * t96 + t164) - g(2) * t168, qJDD(1) * t109 - 0.2e1 * t115 * t162, -0.2e1 * t115 * t169 + 0.2e1 * t173 * t182, t69, t68, 0, t123 * t115 + t118 * t205, -t115 * t205 + t123 * t118, t22 * t59 - t23 * t57 - t30 * t37 + t31 * t36 - t124, t4 * t31 + t18 * t23 - t3 * t30 - t17 * t22 + t133 * t144 + t56 * t184 - g(1) * (t98 * t102 + (-pkin(2) + t113) * t96 + t164) - g(2) * (t96 * t102 - t113 * t98 + t168) t117 * t140 - t166 * t46 (t114 * t46 + t117 * t44) * t61 + (-t195 - t117 * t16 + (t114 * t44 - t117 * t46) * qJD(6)) * t67, t201 + t206, t125 - t139, -t127 * t34 + t209 * t58, -g(1) * t51 - g(2) * t49 + t30 * t16 + t22 * t44 + t5 * t58 + t7 * t127 + (t25 * t209 - t198 + (-t127 * t13 - t209 * t31 + t199) * qJD(6)) * t117 + t204 * t114, g(1) * t50 - g(2) * t48 + t30 * t15 + t22 * t46 - t6 * t58 + (-(-qJD(6) * t31 + t25) * t209 + t198 + t161 * t127 - qJD(6) * t199) * t114 + t204 * t117; 0, 0, 0, t110, 0, 0, t110, 0, 0, 0, 0, 0, t68, -t69, -t127 * t37 + t36 * t67 + t57 * t61 + t58 * t59, -t127 * t3 - t17 * t58 - t18 * t61 + t4 * t67 - g(3), 0, 0, 0, 0, 0, t125 + t139, t201 - t206; 0, 0, 0, 0, qJDD(1), -t121, t132 + t130, 0, 0, 0, 0, 0, t115 * t181 + t171, t118 * t181 - t172, t127 * t36 + t37 * t67 - t57 * t58 + t59 * t61, t124 - t185, 0, 0, 0, 0, 0, t127 * t194 - t16 * t67 + t44 * t61 + (-t114 * t58 + t117 * t148) * t209, t127 * t29 + (-t114 * t148 - t117 * t58) * t209 - t140; 0, 0, 0, 0, 0, 0, 0, t118 * t121 * t115, -t182 * t121, t169, -t170, qJDD(4), -t110 * t115 + t118 * t130 + t53, t55 * qJD(4) + (-qJD(4) * t64 - t110) * t118 + (-t130 - t63) * t115 (t18 - t20) * t59 - (t17 - t21) * t57 + (t186 * t36 + t187 * t37) * pkin(4), t17 * t20 - t18 * t21 + (t187 * t3 + t186 * t4 + g(3) * t115 + (t165 - t185) * t118) * pkin(4), t152 * t46 + t195 (t15 - t212) * t117 + (-t209 * t46 - t16) * t114, -t196 + t211, t136 + t197, -t209 * t59, t126 * t114 - t117 * t203 + t88 * t16 - t20 * t44 - t5 * t59, t114 * t203 + t126 * t117 + t88 * t15 - t20 * t46 + t6 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 ^ 2 - t59 ^ 2, t17 * t59 + t18 * t57 + t133 + t146, 0, 0, 0, 0, 0, t136 - t197, -t196 - t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t44, -t44 ^ 2 + t46 ^ 2, t15 + t212, -t153 + (-qJD(6) + t209) * t46, -t34, -g(1) * t48 - g(2) * t50 + t114 * t135 - t12 * t46 - t13 * t177 + t209 * t6 + t7, g(1) * t49 - g(2) * t51 + t114 * t161 + t117 * t135 + t12 * t44 + t209 * t5;];
tau_reg  = t2;
