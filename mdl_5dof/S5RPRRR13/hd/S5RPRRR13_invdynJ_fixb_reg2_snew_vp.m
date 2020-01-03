% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR13_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:28
% EndTime: 2019-12-31 19:15:35
% DurationCPUTime: 2.10s
% Computational Cost: add. (8045->278), mult. (15790->375), div. (0->0), fcn. (10182->8), ass. (0->192)
t153 = sin(qJ(5));
t159 = cos(qJ(3));
t186 = qJD(1) * qJD(3);
t143 = t159 * t186;
t155 = sin(qJ(3));
t144 = t155 * qJDD(1);
t133 = -t144 - t143;
t128 = qJDD(4) - t133;
t125 = qJDD(5) + t128;
t154 = sin(qJ(4));
t158 = cos(qJ(4));
t190 = qJD(1) * t159;
t129 = -t158 * qJD(3) + t154 * t190;
t131 = qJD(3) * t154 + t158 * t190;
t157 = cos(qJ(5));
t106 = t157 * t129 + t131 * t153;
t108 = -t129 * t153 + t131 * t157;
t76 = t108 * t106;
t217 = -t76 + t125;
t222 = t153 * t217;
t113 = t131 * t129;
t215 = -t113 + t128;
t221 = t154 * t215;
t220 = t157 * t217;
t219 = t158 * t215;
t146 = t159 * qJDD(1);
t181 = t155 * t186;
t134 = t146 - t181;
t178 = -t158 * qJDD(3) + t134 * t154;
t101 = -qJD(4) * t131 - t178;
t170 = -qJDD(3) * t154 - t134 * t158;
t102 = -qJD(4) * t129 - t170;
t63 = -qJD(5) * t106 + t101 * t153 + t102 * t157;
t142 = qJD(1) * t155 + qJD(4);
t140 = qJD(5) + t142;
t93 = t140 * t106;
t218 = t63 - t93;
t162 = qJD(1) ^ 2;
t213 = pkin(6) + pkin(1);
t216 = t213 * t162;
t121 = t142 * t129;
t85 = t102 + t121;
t179 = -t157 * t101 + t102 * t153;
t43 = (qJD(5) - t140) * t108 + t179;
t81 = (qJD(4) - t142) * t131 + t178;
t104 = t106 ^ 2;
t105 = t108 ^ 2;
t126 = t129 ^ 2;
t127 = t131 ^ 2;
t139 = t140 ^ 2;
t141 = t142 ^ 2;
t185 = qJD(2) * qJD(1);
t148 = 0.2e1 * t185;
t150 = qJDD(1) * qJ(2);
t156 = sin(qJ(1));
t160 = cos(qJ(1));
t175 = g(1) * t160 + g(2) * t156;
t167 = -t150 + t175;
t165 = t148 - t167;
t172 = -t134 + t181;
t173 = -t133 + t143;
t80 = t173 * pkin(3) + t172 * pkin(7) + t165 - t216;
t180 = g(1) * t156 - t160 * g(2);
t174 = qJDD(2) - t180;
t201 = qJ(2) * t162;
t163 = t174 - t201;
t120 = -t213 * qJDD(1) + t163;
t110 = t159 * g(3) - t155 * t120;
t161 = qJD(3) ^ 2;
t177 = pkin(3) * t155 - pkin(7) * t159;
t166 = t162 * t177;
t90 = -t161 * pkin(3) + qJDD(3) * pkin(7) - t155 * t166 - t110;
t58 = t154 * t90 - t158 * t80;
t33 = t215 * pkin(4) - t85 * pkin(8) - t58;
t176 = pkin(4) * t142 - pkin(8) * t131;
t59 = t154 * t80 + t158 * t90;
t34 = -t126 * pkin(4) + t101 * pkin(8) - t142 * t176 + t59;
t17 = t153 * t34 - t157 * t33;
t18 = t153 * t33 + t157 * t34;
t7 = t153 * t18 - t157 * t17;
t214 = pkin(4) * t7;
t46 = t63 + t93;
t22 = -t153 * t43 - t157 * t46;
t212 = pkin(4) * t22;
t211 = t154 * t7;
t210 = t158 * t7;
t109 = g(3) * t155 + t120 * t159;
t89 = qJDD(3) * pkin(3) + pkin(7) * t161 - t159 * t166 + t109;
t51 = pkin(4) * t101 + pkin(8) * t126 - t131 * t176 + t89;
t209 = t153 * t51;
t68 = t76 + t125;
t208 = t153 * t68;
t207 = t154 * t89;
t97 = t113 + t128;
t206 = t154 * t97;
t205 = t157 * t51;
t204 = t157 * t68;
t203 = t158 * t89;
t202 = t158 * t97;
t200 = qJDD(1) * pkin(1);
t199 = t140 * t153;
t198 = t140 * t157;
t197 = t142 * t154;
t196 = t142 * t158;
t151 = t155 ^ 2;
t195 = t151 * t162;
t152 = t159 ^ 2;
t194 = t152 * t162;
t182 = t155 * t162 * t159;
t193 = t155 * (qJDD(3) + t182);
t192 = t159 * (qJDD(3) - t182);
t191 = t151 + t152;
t188 = qJD(4) + t142;
t184 = t155 * t76;
t183 = t155 * t113;
t8 = t153 * t17 + t157 * t18;
t31 = t154 * t58 + t158 * t59;
t171 = t154 * t59 - t158 * t58;
t79 = t159 * t109 - t155 * t110;
t72 = -t139 - t104;
t36 = t153 * t72 + t220;
t169 = pkin(4) * t36 - t17;
t168 = qJ(2) + t177;
t88 = -t105 - t139;
t49 = t157 * t88 - t208;
t164 = pkin(4) * t49 - t18;
t136 = t191 * qJDD(1);
t135 = t146 - 0.2e1 * t181;
t132 = t144 + 0.2e1 * t143;
t122 = -t163 + t200;
t119 = -t127 + t141;
t118 = t126 - t141;
t117 = t167 - 0.2e1 * t185 + t216;
t115 = -t193 + t159 * (-t161 - t194);
t114 = t155 * (-t161 - t195) + t192;
t112 = t127 - t126;
t111 = -t127 - t141;
t103 = -t141 - t126;
t95 = t126 + t127;
t92 = -t105 + t139;
t91 = t104 - t139;
t86 = t188 * t129 + t170;
t84 = t102 - t121;
t82 = -t188 * t131 - t178;
t75 = t105 - t104;
t74 = -t111 * t154 - t202;
t73 = t111 * t158 - t206;
t71 = t103 * t158 - t221;
t70 = t103 * t154 + t219;
t66 = (-t106 * t157 + t108 * t153) * t140;
t65 = (-t106 * t153 - t108 * t157) * t140;
t64 = -t104 - t105;
t62 = -qJD(5) * t108 - t179;
t61 = t154 * t85 - t158 * t81;
t56 = t157 * t91 - t208;
t55 = -t153 * t92 + t220;
t54 = t153 * t91 + t204;
t53 = t157 * t92 + t222;
t52 = t155 * t74 + t159 * t86;
t50 = -t153 * t88 - t204;
t48 = t155 * t71 + t159 * t82;
t42 = (qJD(5) + t140) * t108 + t179;
t41 = -t108 * t199 + t157 * t63;
t40 = t108 * t198 + t153 * t63;
t39 = t106 * t198 - t153 * t62;
t38 = t106 * t199 + t157 * t62;
t37 = t157 * t72 - t222;
t35 = t155 * t61 + t159 * t95;
t29 = -pkin(8) * t49 - t205;
t28 = -t154 * t49 + t158 * t50;
t27 = t154 * t50 + t158 * t49;
t26 = t155 * t31 + t159 * t89;
t25 = -pkin(8) * t36 - t209;
t24 = t153 * t46 - t157 * t43;
t23 = -t153 * t218 - t157 * t42;
t21 = -t153 * t42 + t157 * t218;
t20 = -t154 * t36 + t158 * t37;
t19 = t154 * t37 + t158 * t36;
t15 = -pkin(4) * t218 + pkin(8) * t50 - t209;
t14 = t155 * t28 - t159 * t218;
t13 = -pkin(4) * t42 + pkin(8) * t37 + t205;
t12 = t155 * t20 - t159 * t42;
t11 = -t154 * t22 + t158 * t24;
t10 = t154 * t24 + t158 * t22;
t9 = t11 * t155 - t159 * t64;
t6 = pkin(4) * t51 + pkin(8) * t8;
t5 = -pkin(8) * t22 - t7;
t4 = -pkin(4) * t64 + pkin(8) * t24 + t8;
t3 = t158 * t8 - t211;
t2 = t154 * t8 + t210;
t1 = t155 * t3 + t159 * t51;
t16 = [0, 0, 0, 0, 0, qJDD(1), t180, t175, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t174 - 0.2e1 * t200, t148 + 0.2e1 * t150 - t175, pkin(1) * t122 + qJ(2) * (-pkin(1) * t162 + t165), -t172 * t159, -t132 * t159 - t135 * t155, t192 - t155 * (t161 - t194), t173 * t155, t159 * (-t161 + t195) - t193, 0, qJ(2) * t132 - t213 * t114 - t117 * t155, qJ(2) * t135 - t213 * t115 - t117 * t159, t213 * t136 - t191 * t201 - t79, -qJ(2) * t117 - t213 * t79, t159 * (t102 * t158 - t131 * t197) + t183, t159 * (-t154 * t84 + t158 * t82) + t155 * t112, t159 * (-t119 * t154 + t219) + t155 * t85, t159 * (-t101 * t154 + t129 * t196) - t183, t159 * (t118 * t158 - t206) - t155 * t81, t155 * t128 + t159 * (-t129 * t158 + t131 * t154) * t142, t159 * (-pkin(7) * t70 - t207) - t155 * (-pkin(3) * t70 + t58) + qJ(2) * t70 - t213 * t48, t159 * (-pkin(7) * t73 - t203) - t155 * (-pkin(3) * t73 + t59) + qJ(2) * t73 - t213 * t52, -t159 * t171 + t168 * (-t154 * t81 - t158 * t85) - t213 * t35, t168 * t171 - t213 * t26, t159 * (-t154 * t40 + t158 * t41) + t184, t159 * (-t154 * t21 + t158 * t23) + t155 * t75, t159 * (-t154 * t53 + t158 * t55) + t155 * t46, t159 * (-t154 * t38 + t158 * t39) - t184, t159 * (-t154 * t54 + t158 * t56) - t155 * t43, t159 * (-t154 * t65 + t158 * t66) + t155 * t125, t159 * (-pkin(7) * t19 - t13 * t154 + t158 * t25) - t155 * (-pkin(3) * t19 - t169) + qJ(2) * t19 - t213 * t12, t159 * (-pkin(7) * t27 - t15 * t154 + t158 * t29) - t155 * (-pkin(3) * t27 - t164) + qJ(2) * t27 - t213 * t14, t159 * (-pkin(7) * t10 - t154 * t4 + t158 * t5) - t155 * (-pkin(3) * t10 - t212) + qJ(2) * t10 - t213 * t9, t159 * (-pkin(7) * t2 - pkin(8) * t210 - t154 * t6) - t155 * (-pkin(3) * t2 - t214) + qJ(2) * t2 - t213 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t162, -t122, 0, 0, 0, 0, 0, 0, t114, t115, -t136, t79, 0, 0, 0, 0, 0, 0, t48, t52, t35, t26, 0, 0, 0, 0, 0, 0, t12, t14, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, (-t151 + t152) * t162, t146, -t182, -t144, qJDD(3), t109, t110, 0, 0, t102 * t154 + t131 * t196, t154 * t82 + t158 * t84, t119 * t158 + t221, t101 * t158 + t129 * t197, t118 * t154 + t202, (-t129 * t154 - t131 * t158) * t142, pkin(3) * t82 + pkin(7) * t71 + t203, pkin(3) * t86 + pkin(7) * t74 - t207, pkin(3) * t95 + pkin(7) * t61 + t31, pkin(3) * t89 + pkin(7) * t31, t154 * t41 + t158 * t40, t154 * t23 + t158 * t21, t154 * t55 + t158 * t53, t154 * t39 + t158 * t38, t154 * t56 + t158 * t54, t154 * t66 + t158 * t65, -pkin(3) * t42 + pkin(7) * t20 + t13 * t158 + t154 * t25, -pkin(3) * t218 + pkin(7) * t28 + t15 * t158 + t154 * t29, -pkin(3) * t64 + pkin(7) * t11 + t154 * t5 + t158 * t4, pkin(3) * t51 + pkin(7) * t3 - pkin(8) * t211 + t158 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t112, t85, -t113, -t81, t128, -t58, -t59, 0, 0, t76, t75, t46, -t76, -t43, t125, t169, t164, t212, t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t75, t46, -t76, -t43, t125, -t17, -t18, 0, 0;];
tauJ_reg = t16;
