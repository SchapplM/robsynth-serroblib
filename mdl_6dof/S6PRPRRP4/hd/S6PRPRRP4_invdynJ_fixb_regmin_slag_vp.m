% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRP4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:22
% EndTime: 2019-03-08 20:12:31
% DurationCPUTime: 4.06s
% Computational Cost: add. (5111->443), mult. (11891->583), div. (0->0), fcn. (9816->14), ass. (0->222)
t164 = sin(pkin(11));
t170 = sin(qJ(4));
t167 = cos(pkin(11));
t292 = cos(qJ(4));
t239 = t292 * t167;
t198 = -t170 * t164 + t239;
t166 = sin(pkin(6));
t173 = cos(qJ(2));
t258 = t166 * t173;
t184 = t198 * t258;
t100 = qJD(1) * t184;
t288 = pkin(8) + qJ(3);
t139 = t288 * t164;
t140 = t288 * t167;
t199 = -t292 * t139 - t170 * t140;
t53 = qJD(3) * t198 + qJD(4) * t199;
t314 = t100 - t53;
t171 = sin(qJ(2));
t251 = qJD(1) * t166;
t238 = t171 * t251;
t136 = qJD(2) * qJ(3) + t238;
t274 = cos(pkin(6));
t230 = qJD(1) * t274;
t149 = t167 * t230;
t282 = pkin(8) * qJD(2);
t92 = t149 + (-t136 - t282) * t164;
t107 = t167 * t136 + t164 * t230;
t93 = t167 * t282 + t107;
t39 = t170 * t92 + t292 * t93;
t313 = t39 * qJD(4);
t312 = t198 * qJD(2);
t130 = t198 * qJD(4);
t134 = t292 * t164 + t170 * t167;
t187 = t134 * qJDD(2);
t175 = qJD(2) * t130 + t187;
t311 = -qJD(4) * qJD(5) - t175;
t305 = t134 * qJD(2);
t310 = qJD(4) * t305;
t119 = qJD(5) - t312;
t169 = sin(qJ(5));
t172 = cos(qJ(5));
t247 = qJD(5) * t172;
t248 = qJD(5) * t169;
t245 = qJDD(1) * t166;
t234 = t173 * t245;
t250 = qJD(2) * t171;
t235 = qJD(1) * t250;
t304 = t166 * t235 + qJDD(3);
t203 = -t234 + t304;
t272 = qJDD(2) * pkin(2);
t110 = t203 - t272;
t242 = t167 * qJDD(2);
t131 = t134 * qJD(4);
t243 = t164 * qJDD(2);
t218 = -qJDD(2) * t239 + t170 * t243;
t81 = qJD(2) * t131 + t218;
t25 = -pkin(3) * t242 + t81 * pkin(4) - pkin(9) * t175 + t110;
t35 = qJD(4) * pkin(9) + t39;
t154 = t167 * pkin(3) + pkin(2);
t236 = t173 * t251;
t219 = qJD(3) - t236;
t118 = -t154 * qJD(2) + t219;
t47 = -pkin(4) * t312 - pkin(9) * t305 + t118;
t244 = qJDD(2) * qJ(3);
t105 = t171 * t245 + t244 + (qJD(3) + t236) * qJD(2);
t229 = qJDD(1) * t274;
t147 = t167 * t229;
t59 = t147 + (-pkin(8) * qJDD(2) - t105) * t164;
t73 = t167 * t105 + t164 * t229;
t60 = pkin(8) * t242 + t73;
t208 = t170 * t59 + t292 * t60;
t307 = -t170 * t93 + t292 * t92;
t8 = qJDD(4) * pkin(9) + qJD(4) * t307 + t208;
t233 = t169 * t8 - t172 * t25 + t35 * t247 + t47 * t248;
t75 = qJDD(5) + t81;
t298 = pkin(5) * t75;
t2 = qJDD(6) + t233 - t298;
t16 = t169 * t47 + t172 * t35;
t11 = qJ(6) * t119 + t16;
t279 = t11 * t119;
t309 = -t2 + t279;
t185 = t134 * t258;
t98 = -t170 * t139 + t292 * t140;
t284 = -qJD(1) * t185 + qJD(3) * t134 + qJD(4) * t98;
t77 = pkin(4) * t131 - pkin(9) * t130;
t78 = -pkin(4) * t198 - pkin(9) * t134 - t154;
t308 = -t78 * t247 + t98 * t248 + t314 * t172 + (t238 - t77) * t169;
t163 = pkin(11) + qJ(4);
t156 = sin(t163);
t165 = sin(pkin(10));
t273 = cos(pkin(10));
t223 = t274 * t273;
t123 = t165 * t171 - t173 * t223;
t232 = t165 * t274;
t125 = t273 * t171 + t173 * t232;
t226 = g(1) * t125 + g(2) * t123;
t190 = -g(3) * t258 + t226;
t306 = t190 * t156;
t174 = qJD(2) ^ 2;
t191 = (qJDD(2) * t173 - t171 * t174) * t166;
t145 = t169 * t258;
t259 = t166 * t171;
t121 = -t164 * t259 + t274 * t167;
t122 = t274 * t164 + t167 * t259;
t64 = t170 * t121 + t292 * t122;
t50 = t172 * t64 - t145;
t124 = t165 * t173 + t171 * t223;
t126 = -t171 * t232 + t273 * t173;
t157 = cos(t163);
t231 = t166 * t273;
t260 = t165 * t166;
t193 = -g(3) * (-t156 * t259 + t274 * t157) - g(2) * (-t124 * t156 - t157 * t231) - g(1) * (-t126 * t156 + t157 * t260);
t101 = -t172 * qJD(4) + t169 * t305;
t249 = qJD(4) * t169;
t103 = t172 * t305 + t249;
t34 = -qJD(4) * pkin(4) - t307;
t19 = t101 * pkin(5) - t103 * qJ(6) + t34;
t297 = pkin(9) * t75;
t303 = t119 * t19 - t297;
t210 = -t110 + t226;
t302 = t166 * (-g(3) * t173 + t235) + t210 + t272;
t216 = (-t136 * t164 + t149) * t164 - t107 * t167;
t301 = t173 * t216 - (-qJD(2) * pkin(2) + t219) * t171;
t300 = t103 ^ 2;
t299 = t119 ^ 2;
t294 = qJ(6) * t131 - qJD(6) * t198 - t308;
t283 = t169 * t78 + t172 * t98;
t61 = t100 * t169 - t172 * t238;
t293 = -t131 * pkin(5) + t283 * qJD(5) + t169 * t53 - t172 * t77 - t61;
t220 = pkin(5) * t169 - qJ(6) * t172;
t221 = pkin(5) * t172 + qJ(6) * t169;
t287 = t220 * t130 + (t221 * qJD(5) - qJD(6) * t172) * t134 + t284;
t224 = -t172 * qJDD(4) + t247 * t305;
t33 = t169 * t187 + (qJD(5) + t312) * t249 + t224;
t286 = -t101 * t247 - t169 * t33;
t76 = pkin(4) * t305 - pkin(9) * t312;
t285 = t169 * t76 + t172 * t307;
t281 = pkin(9) * qJD(5);
t280 = qJ(6) * t75;
t278 = t119 * t16;
t66 = t169 * t75;
t67 = t172 * t75;
t32 = -t169 * qJDD(4) + t172 * t311 + t305 * t248;
t276 = t32 * t169;
t275 = -t169 * qJD(6) + t119 * t220 - t39;
t271 = t101 * t312;
t270 = t103 * t101;
t228 = t103 * t119;
t269 = t103 * t305;
t268 = t119 * t169;
t267 = t305 * t101;
t266 = t130 * t169;
t265 = t130 * t172;
t264 = t134 * t172;
t262 = t157 * t169;
t261 = t157 * t172;
t256 = t172 * t173;
t255 = t173 * t174;
t15 = -t169 * t35 + t172 * t47;
t254 = qJD(6) - t15;
t253 = qJDD(1) - g(3);
t252 = t164 ^ 2 + t167 ^ 2;
t240 = t166 * t256;
t237 = t166 * t250;
t225 = g(1) * t126 + g(2) * t124;
t206 = t169 * t25 + t172 * t8 + t47 * t247 - t35 * t248;
t1 = qJD(6) * t119 + t206 + t280;
t10 = -pkin(5) * t119 + t254;
t222 = t119 * t10 + t1;
t217 = t10 * t172 - t11 * t169;
t214 = -t119 * t248 + t268 * t312 + t67;
t213 = t66 + (-t172 * t312 + t247) * t119;
t212 = pkin(4) + t221;
t211 = pkin(4) * t157 + pkin(9) * t156 + t154;
t49 = t169 * t64 + t240;
t205 = t170 * t60 - t292 * t59 + t313;
t201 = t119 * t34 - t297;
t200 = t292 * t121 - t170 * t122;
t197 = -t134 * t248 + t265;
t36 = qJD(2) * t184 + qJD(4) * t200;
t14 = t50 * qJD(5) + t169 * t36 - t172 * t237;
t37 = qJD(2) * t185 + qJD(4) * t64;
t196 = t37 * t101 - t119 * t14 - t200 * t33 - t49 * t75;
t108 = t157 * t145 - t172 * t259;
t55 = -t123 * t262 - t124 * t172;
t57 = -t125 * t262 - t126 * t172;
t195 = g(1) * t57 + g(2) * t55 + g(3) * t108;
t109 = (t157 * t256 + t169 * t171) * t166;
t56 = -t123 * t261 + t124 * t169;
t58 = -t125 * t261 + t126 * t169;
t194 = -g(1) * t58 - g(2) * t56 - g(3) * t109;
t113 = t274 * t156 + t157 * t259;
t83 = t124 * t157 - t156 * t231;
t85 = t126 * t157 + t156 * t260;
t192 = -g(1) * t85 - g(2) * t83 - g(3) * t113;
t9 = -qJDD(4) * pkin(4) + t205;
t183 = t190 + t234;
t72 = -t105 * t164 + t147;
t182 = -t72 * t164 + t73 * t167 - t225;
t13 = -qJD(5) * t49 + t169 * t237 + t172 * t36;
t181 = t103 * t37 - t119 * t13 + t200 * t32 - t50 * t75;
t43 = -t123 * t172 + t169 * t83;
t45 = -t125 * t172 + t169 * t85;
t86 = t113 * t169 + t240;
t180 = g(1) * t45 + g(2) * t43 + g(3) * t86 - t233;
t179 = -t119 * t281 + t193;
t3 = t33 * pkin(5) + t32 * qJ(6) - t103 * qJD(6) + t9;
t178 = t179 - t3;
t44 = t123 * t169 + t172 * t83;
t46 = t125 * t169 + t172 * t85;
t87 = t113 * t172 - t145;
t177 = -g(1) * t46 - g(2) * t44 - g(3) * t87 + t206;
t176 = t103 * t19 + qJDD(6) - t180;
t95 = -t154 * qJDD(2) + t203;
t48 = pkin(5) * t103 + qJ(6) * t101;
t40 = t220 * t134 - t199;
t28 = pkin(5) * t198 + t169 * t98 - t172 * t78;
t27 = -qJ(6) * t198 + t283;
t20 = t101 * t119 - t32;
t18 = -pkin(5) * t305 + t169 * t307 - t172 * t76;
t17 = qJ(6) * t305 + t285;
t4 = [t253, 0, t191 (-qJDD(2) * t171 - t255) * t166, t167 * t191, -t164 * t191, t252 * t166 * t255 + (-t121 * t164 + t122 * t167) * qJDD(2), t72 * t121 + t73 * t122 - g(3) + (-t301 * qJD(2) - t110 * t173) * t166, 0, 0, 0, 0, 0, -t37 * qJD(4) + t200 * qJDD(4) + (-t173 * t81 - t250 * t312) * t166, -t36 * qJD(4) - t64 * qJDD(4) + (-t173 * t187 + (-t130 * t173 + t171 * t305) * qJD(2)) * t166, 0, 0, 0, 0, 0, t196, t181, t196, -t101 * t13 + t103 * t14 - t32 * t49 - t33 * t50, -t181, t1 * t50 + t10 * t14 + t11 * t13 + t19 * t37 + t2 * t49 - t200 * t3 - g(3); 0, qJDD(2), t183, -t253 * t259 + t225, t302 * t167, -t302 * t164, -g(3) * t259 + t182 + (t219 * qJD(2) + t244) * t252, -t216 * qJD(3) + t210 * pkin(2) + t182 * qJ(3) + (-g(3) * (pkin(2) * t173 + qJ(3) * t171) + t301 * qJD(1)) * t166, t130 * t305 + t134 * t175, t130 * t312 - t131 * t305 - t134 * t81 + t175 * t198, qJD(4) * t130 + qJDD(4) * t134, -qJD(4) * t131 + qJDD(4) * t198, 0, -t284 * qJD(4) + qJDD(4) * t199 + t118 * t131 - t154 * t81 + t190 * t157 - t198 * t95 + t238 * t312, qJD(4) * t314 - t98 * qJDD(4) + t118 * t130 + t95 * t134 - t154 * t175 - t305 * t238 - t306, t103 * t197 - t32 * t264 (-t101 * t172 - t103 * t169) * t130 + (t276 - t172 * t33 + (t101 * t169 - t103 * t172) * qJD(5)) * t134, t103 * t131 + t119 * t197 + t198 * t32 + t75 * t264, -t134 * t66 - t101 * t131 + t33 * t198 + (-t134 * t247 - t266) * t119, t119 * t131 - t198 * t75, t233 * t198 + t15 * t131 - t199 * t33 + t61 * t119 + t284 * t101 + ((-qJD(5) * t98 + t77) * t119 + t78 * t75 + t34 * qJD(5) * t134) * t172 + ((-qJD(5) * t78 - t53) * t119 - t98 * t75 + t9 * t134 + t34 * t130) * t169 + t194, -t283 * t75 + t206 * t198 - t16 * t131 + t199 * t32 + t34 * t265 + (t9 * t172 - t34 * t248) * t134 + t308 * t119 + t284 * t103 + t195, t19 * t266 - t10 * t131 + t198 * t2 - t28 * t75 + t33 * t40 + (t169 * t3 + t19 * t247) * t134 - t293 * t119 + t287 * t101 + t194, -t27 * t33 - t28 * t32 + t217 * t130 + t293 * t103 - t294 * t101 + t306 + (-t1 * t169 + t2 * t172 + (-t10 * t169 - t11 * t172) * qJD(5)) * t134, -t19 * t265 - t1 * t198 + t11 * t131 + t27 * t75 + t32 * t40 + (-t172 * t3 + t19 * t248) * t134 + t294 * t119 - t287 * t103 - t195, t1 * t27 + t3 * t40 + t2 * t28 - g(1) * (t58 * pkin(5) + t57 * qJ(6) + t126 * t288) - g(2) * (t56 * pkin(5) + t55 * qJ(6) + t124 * t288) + t287 * t19 + t294 * t11 + t293 * t10 + t226 * t211 + (-t109 * pkin(5) - t108 * qJ(6) - (t171 * t288 + t173 * t211) * t166) * g(3); 0, 0, 0, 0, -t242, t243, -t252 * t174, qJD(2) * t216 - t183 - t272 + t304, 0, 0, 0, 0, 0, t218 + 0.2e1 * t310, 0.2e1 * qJD(4) * t312 + t187, 0, 0, 0, 0, 0, t214 - t267, -t299 * t172 - t269 - t66, -t119 * t268 - t267 + t67 (t32 + t271) * t172 + t169 * t228 + t286, t213 + t269, t222 * t169 + t309 * t172 - t19 * t305 - t190; 0, 0, 0, 0, 0, 0, 0, 0, -t305 * t312, t305 ^ 2 - t312 ^ 2, t187, -t218, qJDD(4), -t118 * t305 + t193 - t205 + t313, -t118 * t312 - t192 - t208, t172 * t228 - t276 (-t32 + t271) * t172 - t103 * t268 + t286, t213 - t269, t214 + t267, -t119 * t305, -pkin(4) * t33 - t39 * t101 - t15 * t305 + (t119 * t307 + t201) * t169 + (-t9 + (-t76 - t281) * t119 + t193) * t172, pkin(4) * t32 + t285 * t119 + t16 * t305 - t39 * t103 + t201 * t172 + (-t179 + t9) * t169, t10 * t305 + t275 * t101 + t119 * t18 + t303 * t169 + t178 * t172 - t212 * t33, t101 * t17 - t103 * t18 + ((qJD(5) * t103 - t33) * pkin(9) + t222) * t172 + ((qJD(5) * t101 - t32) * pkin(9) - t309) * t169 + t192, -t275 * t103 - t11 * t305 - t119 * t17 + t178 * t169 - t303 * t172 - t212 * t32, -t10 * t18 - t11 * t17 + t275 * t19 + (qJD(5) * t217 + t1 * t172 + t2 * t169 + t192) * pkin(9) + (-t3 + t193) * t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270, -t101 ^ 2 + t300, t20, t169 * t311 - t224 + t228, t75, -t103 * t34 + t180 + t278, t101 * t34 + t119 * t15 - t177, -t101 * t48 - t176 + t278 + 0.2e1 * t298, pkin(5) * t32 - t33 * qJ(6) + (t11 - t16) * t103 + (t10 - t254) * t101, 0.2e1 * t280 - t101 * t19 + t103 * t48 + (0.2e1 * qJD(6) - t15) * t119 + t177, t1 * qJ(6) - t2 * pkin(5) - t19 * t48 - t10 * t16 - g(1) * (-pkin(5) * t45 + qJ(6) * t46) - g(2) * (-pkin(5) * t43 + qJ(6) * t44) - g(3) * (-pkin(5) * t86 + qJ(6) * t87) + t254 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(5) - t218 + t270 - t310, t20, -t299 - t300, t176 - t279 - t298;];
tau_reg  = t4;
