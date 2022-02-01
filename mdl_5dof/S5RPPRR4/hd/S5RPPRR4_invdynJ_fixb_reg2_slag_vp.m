% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:59
% EndTime: 2022-01-23 09:17:06
% DurationCPUTime: 3.82s
% Computational Cost: add. (5009->361), mult. (12538->496), div. (0->0), fcn. (9478->14), ass. (0->206)
t175 = cos(pkin(8));
t226 = t175 * qJD(1);
t286 = qJD(4) - t226;
t176 = sin(qJ(5));
t179 = cos(qJ(5));
t228 = qJD(5) * t179;
t229 = qJD(5) * t176;
t173 = sin(pkin(8));
t177 = sin(qJ(4));
t174 = cos(pkin(9));
t180 = cos(qJ(4));
t242 = t180 * t174;
t218 = t173 * t242;
t172 = sin(pkin(9));
t223 = qJDD(1) * t172;
t124 = t180 * t172 + t177 * t174;
t279 = t124 * qJD(4);
t56 = -qJDD(1) * t218 + t173 * (qJD(1) * t279 + t177 * t223);
t234 = qJD(1) * t173;
t216 = t172 * t234;
t203 = t177 * t216;
t230 = qJD(4) * t180;
t215 = t174 * t230;
t57 = -qJD(4) * t203 + (qJD(1) * t215 + qJDD(1) * t124) * t173;
t192 = qJD(1) * t124;
t94 = t173 * t192;
t97 = qJD(1) * t218 - t203;
t16 = t176 * t57 + t179 * t56 + t94 * t228 + t97 * t229;
t141 = -qJD(5) - t286;
t43 = t176 * t97 + t179 * t94;
t257 = t43 * t141;
t285 = -t16 - t257;
t198 = -t176 * t94 + t179 * t97;
t267 = t198 ^ 2;
t268 = t43 ^ 2;
t284 = t267 - t268;
t266 = t43 * t198;
t123 = -t177 * t172 + t242;
t254 = t286 * t123;
t253 = t175 * t192 - t279;
t17 = qJD(5) * t198 - t176 * t56 + t179 * t57;
t258 = t198 * t141;
t283 = -t17 - t258;
t129 = pkin(2) * t175 + t173 * qJ(3) + pkin(1);
t114 = -qJD(1) * t129 + qJD(2);
t103 = t174 * t114;
t190 = -t174 * t173 * pkin(6) + (-qJ(2) * t172 - pkin(3)) * t175;
t58 = qJD(1) * t190 + t103;
t217 = qJ(2) * t226;
t72 = t172 * t114 + t174 * t217;
t63 = -pkin(6) * t216 + t72;
t30 = -t177 * t63 + t180 * t58;
t22 = -t97 * pkin(7) + t30;
t20 = pkin(4) * t286 + t22;
t31 = t177 * t58 + t180 * t63;
t23 = -t94 * pkin(7) + t31;
t224 = qJD(1) * qJD(2);
t213 = t175 * t224;
t232 = qJD(3) * t173;
t83 = -qJD(1) * t232 - qJDD(1) * t129 + qJDD(2);
t76 = t174 * t83;
t39 = qJDD(1) * t190 - t172 * t213 + t76;
t221 = t173 * qJDD(1);
t212 = t172 * t221;
t220 = t175 * qJDD(1);
t210 = t174 * t220;
t55 = qJ(2) * t210 + t172 * t83 + t174 * t213;
t41 = -pkin(6) * t212 + t55;
t13 = -qJD(4) * t31 - t177 * t41 + t180 * t39;
t144 = -qJDD(4) + t220;
t6 = -t144 * pkin(4) + t56 * pkin(7) + t13;
t231 = qJD(4) * t177;
t12 = t177 * t39 + t180 * t41 + t58 * t230 - t63 * t231;
t7 = -t57 * pkin(7) + t12;
t1 = (qJD(5) * t20 + t7) * t179 + t176 * t6 - t23 * t229;
t171 = pkin(9) + qJ(4);
t159 = qJ(5) + t171;
t153 = cos(t159);
t270 = g(3) * t173;
t142 = qJ(2) * t234 + qJD(3);
t115 = pkin(3) * t216 + t142;
t61 = t94 * pkin(4) + t115;
t152 = sin(t159);
t181 = cos(qJ(1));
t241 = t181 * t152;
t178 = sin(qJ(1));
t245 = t178 * t153;
t90 = -t175 * t245 + t241;
t240 = t181 * t153;
t246 = t178 * t152;
t92 = t175 * t240 + t246;
t282 = g(1) * t92 - g(2) * t90 + t153 * t270 + t61 * t43 - t1;
t165 = g(2) * t181;
t272 = g(1) * t178;
t280 = -t165 + t272;
t225 = qJ(2) * qJDD(1);
t157 = cos(t171);
t238 = t181 * t157;
t156 = sin(t171);
t244 = t178 * t156;
t104 = t175 * t244 + t238;
t239 = t181 * t156;
t243 = t178 * t157;
t106 = -t175 * t239 + t243;
t278 = -g(1) * t106 + g(2) * t104 + t156 * t270;
t259 = t179 * t23;
t9 = t176 * t20 + t259;
t2 = -qJD(5) * t9 - t176 * t7 + t179 * t6;
t89 = t175 * t246 + t240;
t91 = -t175 * t241 + t245;
t277 = -g(1) * t91 + g(2) * t89 + t152 * t270 - t61 * t198 + t2;
t275 = t97 ^ 2;
t274 = t57 * pkin(4);
t269 = t172 * pkin(3);
t265 = t97 * t94;
t264 = qJ(3) + pkin(6);
t263 = t174 * pkin(3) + pkin(2);
t69 = t179 * t123 - t176 * t124;
t262 = qJD(5) * t69 + t253 * t176 + t254 * t179;
t70 = t176 * t123 + t179 * t124;
t261 = -qJD(5) * t70 - t254 * t176 + t253 * t179;
t122 = t174 * t129;
t66 = -t122 + t190;
t250 = t172 * t173;
t87 = t174 * t175 * qJ(2) - t172 * t129;
t73 = -pkin(6) * t250 + t87;
t35 = t177 * t66 + t180 * t73;
t260 = t176 * t23;
t256 = t94 * t286;
t255 = t97 * t286;
t252 = qJDD(1) * pkin(1);
t251 = (pkin(4) * t157 + t263) * t175;
t249 = t172 * t175;
t248 = t173 * t181;
t182 = qJD(1) ^ 2;
t247 = t175 * t182;
t125 = pkin(3) * t250 + t173 * qJ(2);
t160 = t178 * qJ(2);
t237 = t181 * pkin(1) + t160;
t167 = t172 ^ 2;
t169 = t174 ^ 2;
t236 = -t167 - t169;
t168 = t173 ^ 2;
t170 = t175 ^ 2;
t235 = t168 + t170;
t233 = qJD(2) * t175;
t227 = t173 * qJD(2);
t222 = t168 * qJDD(1);
t120 = qJ(2) * t221 + t173 * t224 + qJDD(3);
t34 = -t177 * t73 + t180 * t66;
t206 = t235 * t182;
t205 = 0.2e1 * t173 * t220;
t154 = qJDD(2) - t252;
t204 = 0.2e1 * t235;
t84 = pkin(3) * t212 + t120;
t147 = t173 * t272;
t201 = -g(2) * t248 + t147;
t200 = g(1) * t181 + g(2) * t178;
t112 = t123 * t173;
t28 = -t175 * pkin(4) - t112 * pkin(7) + t34;
t111 = t124 * t173;
t29 = -t111 * pkin(7) + t35;
t14 = -t176 * t29 + t179 * t28;
t15 = t176 * t28 + t179 * t29;
t60 = -t176 * t111 + t179 * t112;
t197 = t224 + t225;
t196 = t154 - t252 + t165;
t116 = -t172 * t233 - t174 * t232;
t54 = -t197 * t249 + t76;
t86 = -qJ(2) * t249 - t122;
t195 = -t116 * qJD(1) - t86 * qJDD(1) - t54;
t117 = -t172 * t232 + t174 * t233;
t194 = t117 * qJD(1) + t87 * qJDD(1) + t55;
t193 = t204 * t224;
t24 = t177 * t116 + t180 * t117 + t66 * t230 - t73 * t231;
t25 = -t35 * qJD(4) + t180 * t116 - t177 * t117;
t163 = g(3) * t175;
t188 = -t173 * t200 + t163 + t84;
t186 = t120 * t173 + t168 * t197 - t200;
t166 = -pkin(7) - t264;
t162 = t181 * qJ(2);
t155 = t170 * qJDD(1);
t148 = qJ(2) + t269;
t138 = -qJDD(5) + t144;
t128 = pkin(4) * t156 + t269;
t113 = t264 * t173 + t263 * t175 + pkin(1);
t107 = t175 * t238 + t244;
t105 = -t175 * t243 + t239;
t101 = t173 * t215 - t231 * t250;
t100 = t173 * t279;
t88 = t94 ^ 2;
t74 = t101 * pkin(4) + t227;
t71 = -t172 * t217 + t103;
t68 = t111 * pkin(4) + t125;
t59 = t179 * t111 + t176 * t112;
t36 = t84 + t274;
t27 = qJD(5) * t60 - t176 * t100 + t179 * t101;
t26 = t179 * t100 + t176 * t101 + t111 * t228 + t112 * t229;
t19 = t100 * pkin(7) + t25;
t18 = -t101 * pkin(7) + t24;
t11 = t179 * t22 - t260;
t10 = -t176 * t22 - t259;
t8 = t179 * t20 - t260;
t4 = -qJD(5) * t15 - t176 * t18 + t179 * t19;
t3 = qJD(5) * t14 + t176 * t19 + t179 * t18;
t5 = [0, 0, 0, 0, 0, qJDD(1), t280, t200, 0, 0, t222, t205, 0, t155, 0, 0, (-t196 + t272) * t175, t173 * t196 - t147, t204 * t225 + t193 - t200, -t154 * pkin(1) - g(1) * (-t178 * pkin(1) + t162) - g(2) * t237 + (t235 * t225 + t193) * qJ(2), t169 * t222, -0.2e1 * t174 * t172 * t222, -0.2e1 * t173 * t210, t167 * t222, t172 * t205, t155, (t174 * t280 + t195) * t175 + t186 * t172, (-t172 * t280 + t194) * t175 + t186 * t174, t147 + (-t172 * t194 + t174 * t195 - t165) * t173, t55 * t87 + t72 * t117 + t54 * t86 + t71 * t116 - g(1) * (-t129 * t178 + t162) - g(2) * (t129 * t181 + t160) + (t120 * qJ(2) + t142 * qJD(2)) * t173, -t97 * t100 - t56 * t112, t100 * t94 - t97 * t101 + t56 * t111 - t112 * t57, -t100 * t286 - t112 * t144 + t56 * t175, t94 * t101 + t57 * t111, -t101 * t286 + t111 * t144 + t57 * t175, t144 * t175, -g(1) * t105 - g(2) * t107 + t115 * t101 + t84 * t111 + t125 * t57 - t13 * t175 - t34 * t144 + t227 * t94 + t25 * t286, -g(1) * t104 - g(2) * t106 - t115 * t100 + t84 * t112 + t12 * t175 - t125 * t56 + t35 * t144 + t227 * t97 - t24 * t286, t30 * t100 - t31 * t101 - t12 * t111 - t13 * t112 - t24 * t94 - t25 * t97 + t34 * t56 - t35 * t57 + t201, t12 * t35 + t31 * t24 + t13 * t34 + t30 * t25 + t84 * t125 + t115 * t227 - g(1) * (-t113 * t178 + t148 * t181) - g(2) * (t113 * t181 + t148 * t178), -t16 * t60 - t198 * t26, t16 * t59 - t60 * t17 - t198 * t27 + t26 * t43, -t60 * t138 + t26 * t141 + t16 * t175, t17 * t59 + t43 * t27, t59 * t138 + t27 * t141 + t17 * t175, t138 * t175, -g(1) * t90 - g(2) * t92 - t14 * t138 - t4 * t141 + t68 * t17 - t2 * t175 + t61 * t27 + t36 * t59 + t74 * t43, -g(1) * t89 - g(2) * t91 + t1 * t175 + t15 * t138 + t3 * t141 - t68 * t16 + t198 * t74 - t61 * t26 + t36 * t60, -t1 * t59 + t14 * t16 - t15 * t17 - t198 * t4 - t2 * t60 + t8 * t26 - t9 * t27 - t3 * t43 + t201, t1 * t15 + t9 * t3 + t2 * t14 + t8 * t4 + t36 * t68 + t61 * t74 - g(1) * (t181 * t128 + t162) - g(2) * (-t166 * t248 + t181 * t251 + t237) + (-g(1) * (t166 * t173 - pkin(1) - t251) - g(2) * t128) * t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t221, -t206, -qJ(2) * t206 + t154 - t280, 0, 0, 0, 0, 0, 0, -t172 * t206 - t210, t172 * t220 - t174 * t206, t236 * t221, t55 * t172 + t54 * t174 + (-t142 * t173 + (t172 * t71 - t174 * t72) * t175) * qJD(1) - t280, 0, 0, 0, 0, 0, 0, -t123 * t144 - t234 * t94 + t253 * t286, t124 * t144 - t234 * t97 - t254 * t286, t123 * t56 - t124 * t57 - t253 * t97 - t254 * t94, -t115 * t234 + t12 * t124 + t13 * t123 + t253 * t30 + t254 * t31 - t280, 0, 0, 0, 0, 0, 0, -t69 * t138 - t141 * t261 - t234 * t43, t70 * t138 + t141 * t262 - t198 * t234, t69 * t16 - t70 * t17 - t198 * t261 - t262 * t43, t1 * t70 + t2 * t69 - t234 * t61 + t261 * t8 + t262 * t9 - t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t174 * t247 + t223) * t173, (qJDD(1) * t174 + t172 * t247) * t173, t236 * t182 * t168, t163 + ((t172 * t72 + t174 * t71) * qJD(1) - t200) * t173 + t120, 0, 0, 0, 0, 0, 0, t255 + t57, -t56 - t256, -t88 - t275, t30 * t97 + t31 * t94 + t188, 0, 0, 0, 0, 0, 0, t17 - t258, -t16 + t257, -t267 - t268, t198 * t8 + t43 * t9 + t188 + t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, -t88 + t275, -t56 + t256, -t265, -t57 + t255, -t144, -t115 * t97 + t286 * t31 + t13 + t278, g(1) * t107 - g(2) * t105 + t115 * t94 + t157 * t270 + t286 * t30 - t12, 0, 0, t266, t284, t285, -t266, t283, -t138, t10 * t141 + (-t138 * t179 + t141 * t229 - t43 * t97) * pkin(4) + t277, -t11 * t141 + (t138 * t176 + t141 * t228 - t198 * t97) * pkin(4) + t282, t10 * t198 + t11 * t43 + t9 * t198 - t8 * t43 + (t16 * t179 - t17 * t176 + (t176 * t198 - t179 * t43) * qJD(5)) * pkin(4), -t8 * t10 - t9 * t11 + (t1 * t176 + t2 * t179 - t61 * t97 + (-t8 * t176 + t9 * t179) * qJD(5) + t278) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266, t284, t285, -t266, t283, -t138, -t9 * t141 + t277, -t8 * t141 + t282, 0, 0;];
tau_reg = t5;
