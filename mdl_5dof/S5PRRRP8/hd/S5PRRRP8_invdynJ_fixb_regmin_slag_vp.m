% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:38
% EndTime: 2019-12-05 17:00:47
% DurationCPUTime: 2.79s
% Computational Cost: add. (2281->354), mult. (5328->486), div. (0->0), fcn. (4119->10), ass. (0->175)
t122 = cos(qJ(3));
t189 = qJD(2) * t122;
t245 = qJD(4) - t189;
t119 = sin(qJ(3));
t190 = qJD(2) * t119;
t244 = qJD(4) * t190 - qJDD(3);
t118 = sin(qJ(4));
t121 = cos(qJ(4));
t178 = t119 * qJDD(2);
t36 = t118 * ((qJD(4) + t189) * qJD(3) + t178) + t244 * t121;
t151 = pkin(3) * t122 + pkin(8) * t119 + pkin(2);
t185 = qJD(4) * t118;
t116 = sin(pkin(5));
t193 = qJD(1) * t116;
t120 = sin(qJ(2));
t123 = cos(qJ(2));
t197 = t122 * t123;
t240 = t118 * t197 - t120 * t121;
t156 = pkin(3) * t119 - pkin(8) * t122;
t89 = t156 * qJD(3);
t243 = -t121 * t89 - t151 * t185 - t240 * t193;
t184 = qJD(4) * t121;
t61 = (t118 * t120 + t121 * t197) * t116;
t242 = -qJD(1) * t61 + t118 * t89 - t151 * t184;
t117 = cos(pkin(5));
t201 = t117 * t122;
t90 = qJD(2) * pkin(7) + t120 * t193;
t241 = qJD(1) * t201 - t119 * t90;
t115 = sin(pkin(9));
t210 = cos(pkin(9));
t163 = t116 * t210;
t203 = t116 * t122;
t162 = t210 * t120;
t205 = t115 * t123;
t71 = t117 * t162 + t205;
t161 = t210 * t123;
t206 = t115 * t120;
t73 = -t117 * t206 + t161;
t204 = t116 * t120;
t75 = t119 * t204 - t201;
t239 = g(3) * t75 - g(2) * (-t71 * t119 - t122 * t163) - g(1) * (t115 * t203 - t119 * t73);
t50 = -qJD(3) * pkin(3) - t241;
t182 = t121 * qJD(3);
t84 = t118 * t190 - t182;
t188 = qJD(3) * t118;
t86 = t121 * t190 + t188;
t16 = pkin(4) * t84 - qJ(5) * t86 + t50;
t111 = t122 * qJDD(2);
t180 = qJD(2) * qJD(3);
t82 = t119 * t180 + qJDD(4) - t111;
t233 = pkin(8) * t82;
t238 = -t16 * t245 + t233;
t124 = qJD(3) ^ 2;
t181 = qJD(1) * qJD(2);
t169 = t120 * t181;
t202 = t116 * t123;
t154 = -qJDD(1) * t202 + t116 * t169;
t70 = -t117 * t161 + t206;
t229 = g(2) * t70;
t72 = t117 * t205 + t162;
t231 = g(1) * t72;
t157 = t229 + t231;
t237 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t124 + (-g(3) * t123 + t169) * t116 - t154 + t157;
t235 = t86 ^ 2;
t234 = pkin(4) * t82;
t192 = qJD(1) * t119;
t102 = t117 * t192;
t58 = t122 * t90 + t102;
t51 = qJD(3) * pkin(8) + t58;
t176 = t123 * t193;
t59 = -t151 * qJD(2) - t176;
t15 = t118 * t59 + t121 * t51;
t7 = qJ(5) * t245 + t15;
t227 = t245 * t7;
t226 = t86 * t84;
t183 = qJD(4) * t122;
t187 = qJD(3) * t119;
t225 = qJ(5) * t187 - qJD(5) * t122 + (-t118 * t183 - t119 * t182) * pkin(7) + t242;
t224 = -pkin(4) * t187 + (-t118 * t187 + t121 * t183) * pkin(7) + t243;
t88 = t156 * qJD(2);
t223 = t118 * t88 + t121 * t241;
t152 = pkin(4) * t118 - qJ(5) * t121;
t222 = -qJD(5) * t118 + t245 * t152 - t58;
t220 = pkin(8) * qJD(4);
t219 = qJ(5) * t82;
t218 = qJD(2) * pkin(2);
t217 = t245 * t15;
t216 = t245 * t84;
t215 = t245 * t86;
t168 = t122 * t180;
t35 = -qJD(4) * t182 + (-t168 - t178) * t121 + t244 * t118;
t214 = t118 * t35;
t213 = t121 * t86;
t212 = t121 * t151;
t198 = t121 * t122;
t211 = pkin(7) * t198 - t118 * t151;
t209 = qJD(3) * t84;
t208 = qJD(3) * t86;
t200 = t118 * t122;
t14 = -t118 * t51 + t121 * t59;
t196 = qJD(5) - t14;
t195 = qJDD(1) - g(3);
t113 = t119 ^ 2;
t194 = -t122 ^ 2 + t113;
t191 = qJD(2) * t116;
t186 = qJD(3) * t122;
t179 = qJDD(1) * t117;
t174 = t120 * t191;
t173 = t123 * t191;
t172 = t245 * t188;
t171 = t245 * t182;
t170 = t245 * t185;
t167 = t123 * t180;
t165 = t119 * t179;
t63 = qJDD(2) * pkin(7) + (qJDD(1) * t120 + t123 * t181) * t116;
t10 = qJDD(3) * pkin(8) + qJD(3) * t241 + t122 * t63 + t165;
t28 = qJD(2) * t89 - t151 * qJDD(2) + t154;
t164 = t118 * t10 - t121 * t28 + t51 * t184 + t59 * t185;
t159 = t84 * t176;
t158 = t86 * t176;
t6 = -pkin(4) * t245 + t196;
t155 = -t118 * t7 + t121 * t6;
t153 = pkin(4) * t121 + qJ(5) * t118;
t125 = qJD(2) ^ 2;
t150 = qJDD(2) * t123 - t120 * t125;
t149 = pkin(3) + t153;
t148 = pkin(7) + t152;
t76 = t117 * t119 + t120 * t203;
t45 = t118 * t76 + t121 * t202;
t46 = -t118 * t202 + t121 * t76;
t145 = t121 * t10 + t118 * t28 + t59 * t184 - t51 * t185;
t142 = -t118 * t82 - t184 * t245;
t141 = t121 * t82 - t170;
t43 = -t75 * qJD(3) + t122 * t173;
t4 = t46 * qJD(4) + t118 * t43 - t121 * t174;
t44 = t76 * qJD(3) + t119 * t173;
t140 = -t245 * t4 + t75 * t36 + t44 * t84 - t45 * t82;
t29 = -t71 * t121 - t70 * t200;
t31 = -t73 * t121 - t72 * t200;
t60 = t240 * t116;
t139 = g(1) * t31 + g(2) * t29 + g(3) * t60;
t30 = t118 * t71 - t70 * t198;
t32 = t118 * t73 - t72 * t198;
t138 = -g(1) * t32 - g(2) * t30 - g(3) * t61;
t40 = -t119 * t163 + t122 * t71;
t42 = t115 * t116 * t119 + t122 * t73;
t136 = g(1) * t42 + g(2) * t40 + g(3) * t76;
t135 = qJD(3) * t102 + t119 * t63 - t122 * t179 + t90 * t186;
t134 = t245 * t50 - t233;
t11 = -qJDD(3) * pkin(3) + t135;
t5 = -t45 * qJD(4) + t118 * t174 + t121 * t43;
t133 = -t245 * t5 - t35 * t75 + t44 * t86 - t46 * t82;
t18 = t118 * t40 - t70 * t121;
t20 = t118 * t42 - t72 * t121;
t132 = g(1) * t20 + g(2) * t18 + g(3) * t45 - t164;
t131 = t220 * t245 - t239;
t3 = pkin(4) * t36 + qJ(5) * t35 - qJD(5) * t86 + t11;
t129 = -t131 - t3;
t91 = -t176 - t218;
t128 = -pkin(7) * qJDD(3) + (t176 + t91 - t218) * qJD(3);
t19 = t118 * t70 + t121 * t40;
t21 = t118 * t72 + t121 * t42;
t127 = -g(1) * t21 - g(2) * t19 - g(3) * t46 + t145;
t126 = t16 * t86 + qJDD(5) - t132;
t68 = t148 * t119;
t53 = t212 + (pkin(7) * t118 + pkin(4)) * t122;
t52 = -qJ(5) * t122 + t211;
t38 = pkin(4) * t86 + qJ(5) * t84;
t27 = (t153 * qJD(4) - qJD(5) * t121) * t119 + t148 * t186;
t26 = -pkin(4) * t190 + t118 * t241 - t121 * t88;
t25 = qJ(5) * t190 + t223;
t13 = -t35 + t216;
t2 = qJDD(5) + t164 - t234;
t1 = qJD(5) * t245 + t145 + t219;
t8 = [t195, 0, t150 * t116, (-qJDD(2) * t120 - t123 * t125) * t116, 0, 0, 0, 0, 0, -qJD(3) * t44 - qJDD(3) * t75 + (-t119 * t167 + t150 * t122) * t116, -qJD(3) * t43 - qJDD(3) * t76 + (-t150 * t119 - t122 * t167) * t116, 0, 0, 0, 0, 0, t140, t133, t140, -t35 * t45 - t36 * t46 + t4 * t86 - t5 * t84, -t133, t1 * t46 + t16 * t44 + t2 * t45 + t3 * t75 + t4 * t6 + t5 * t7 - g(3); 0, qJDD(2), t195 * t202 + t157, g(1) * t73 + g(2) * t71 - t195 * t204, qJDD(2) * t113 + 0.2e1 * t119 * t168, 0.2e1 * t119 * t111 - 0.2e1 * t194 * t180, qJDD(3) * t119 + t122 * t124, qJDD(3) * t122 - t119 * t124, 0, t128 * t119 + t122 * t237, -t119 * t237 + t128 * t122, t122 * t86 * t182 + (-t121 * t35 - t86 * t185) * t119, (-t118 * t86 - t121 * t84) * t186 + (t214 - t121 * t36 + (t118 * t84 - t213) * qJD(4)) * t119, (t35 + t171) * t122 + (t141 + t208) * t119, (t36 - t172) * t122 + (t142 - t209) * t119, -t122 * t82 + t187 * t245, -t82 * t212 - t243 * t245 + (t50 * t188 + (t142 + t209) * pkin(7) + t164) * t122 + (-t159 + t50 * t184 + t14 * qJD(3) + t11 * t118 + (t36 + t172) * pkin(7)) * t119 + t138, -t211 * t82 - t242 * t245 + (t50 * t182 + (t170 + t208) * pkin(7) + t145) * t122 + (-t158 - t50 * t185 - t15 * qJD(3) + t11 * t121 + (-t35 + t171) * pkin(7)) * t119 + t139, t27 * t84 + t36 * t68 - t53 * t82 + (t16 * t188 + t2) * t122 - t224 * t245 + (-qJD(3) * t6 + t118 * t3 + t16 * t184 - t159) * t119 + t138, -t35 * t53 - t36 * t52 + t224 * t86 - t225 * t84 + t155 * t186 + (-g(3) * t202 - t1 * t118 + t121 * t2 + (-t118 * t6 - t121 * t7) * qJD(4) + t157) * t119, -t27 * t86 + t35 * t68 + t52 * t82 + (-t16 * t182 - t1) * t122 + t225 * t245 + (qJD(3) * t7 - t121 * t3 + t16 * t185 + t158) * t119 - t139, t1 * t52 + t3 * t68 + t16 * t27 + t2 * t53 - g(1) * (pkin(4) * t32 + pkin(7) * t73 + qJ(5) * t31) - g(2) * (pkin(4) * t30 + pkin(7) * t71 + qJ(5) * t29) - g(3) * (pkin(4) * t61 + qJ(5) * t60) + t151 * t231 + t151 * t229 + t225 * t7 + t224 * t6 + (-g(3) * pkin(7) * t120 + (-g(3) * t151 - t16 * t192) * t123) * t116; 0, 0, 0, 0, -t119 * t125 * t122, t194 * t125, t178, t111, qJDD(3), qJD(3) * t58 - t91 * t190 - t135 + t239, -t165 + (-qJD(2) * t91 - t63) * t122 + t136, t213 * t245 - t214, (-t35 - t216) * t121 + (-t215 - t36) * t118, (-t119 * t86 - t198 * t245) * qJD(2) - t142, (t119 * t84 + t200 * t245) * qJD(2) + t141, -t245 * t190, -t14 * t190 - pkin(3) * t36 - t58 * t84 + (t241 * t245 + t134) * t118 + (-t11 - (t88 + t220) * t245 + t239) * t121, pkin(3) * t35 + t223 * t245 + t15 * t190 - t58 * t86 + t134 * t121 + (t11 + t131) * t118, -t118 * t238 + t129 * t121 - t149 * t36 + t6 * t190 + t222 * t84 + t245 * t26, t25 * t84 - t26 * t86 + (t1 + t245 * t6 + (qJD(4) * t86 - t36) * pkin(8)) * t121 + (t2 - t227 + (qJD(4) * t84 - t35) * pkin(8)) * t118 - t136, t129 * t118 + t121 * t238 - t149 * t35 - t7 * t190 - t222 * t86 - t245 * t25, -t7 * t25 - t6 * t26 + t222 * t16 + (qJD(4) * t155 + t1 * t121 + t2 * t118 - t136) * pkin(8) + (-t3 + t239) * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, -t84 ^ 2 + t235, t13, -t36 + t215, t82, -t50 * t86 + t132 + t217, t14 * t245 + t50 * t84 - t127, -t38 * t84 - t126 + t217 + 0.2e1 * t234, pkin(4) * t35 - qJ(5) * t36 + (-t15 + t7) * t86 + (t6 - t196) * t84, 0.2e1 * t219 - t16 * t84 + t38 * t86 - (-0.2e1 * qJD(5) + t14) * t245 + t127, t1 * qJ(5) - t2 * pkin(4) - t16 * t38 - t6 * t15 - g(1) * (-pkin(4) * t20 + qJ(5) * t21) - g(2) * (-pkin(4) * t18 + qJ(5) * t19) - g(3) * (-pkin(4) * t45 + qJ(5) * t46) + t196 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82 + t226, t13, -t245 ^ 2 - t235, t126 - t227 - t234;];
tau_reg = t8;
