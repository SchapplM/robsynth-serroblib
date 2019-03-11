% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:55:56
% EndTime: 2019-03-09 04:56:05
% DurationCPUTime: 4.01s
% Computational Cost: add. (3903->498), mult. (7243->587), div. (0->0), fcn. (4315->6), ass. (0->237)
t149 = sin(qJ(4));
t330 = t149 * qJ(5) + pkin(3);
t306 = pkin(4) + qJ(6);
t153 = cos(qJ(3));
t242 = qJD(1) * qJD(3);
t219 = t153 * t242;
t150 = sin(qJ(3));
t239 = t150 * qJDD(1);
t83 = qJDD(4) + t219 + t239;
t229 = t306 * t83;
t245 = qJD(4) * t153;
t223 = t149 * t245;
t152 = cos(qJ(4));
t243 = t152 * qJD(3);
t329 = t150 * t243 + t223;
t142 = t153 * pkin(8);
t323 = t150 * pkin(3) + qJ(2) - t142;
t67 = t323 * qJD(1);
t155 = -pkin(1) - pkin(7);
t114 = t155 * qJD(1) + qJD(2);
t93 = t150 * t114;
t74 = qJD(3) * pkin(8) + t93;
t27 = t149 * t74 - t152 * t67;
t251 = qJD(3) * t149;
t252 = qJD(1) * t153;
t87 = t152 * t252 + t251;
t191 = t87 * pkin(5) + t27;
t264 = qJD(5) + t191;
t254 = qJD(1) * t150;
t119 = qJD(4) + t254;
t205 = qJD(4) * t150 + qJD(1);
t249 = qJD(3) * t153;
t227 = t149 * t249;
t250 = qJD(3) * t150;
t228 = t149 * t250;
t238 = t153 * qJDD(1);
t280 = qJD(4) * t87;
t31 = -qJD(1) * t228 - t152 * qJDD(3) + t149 * t238 + t280;
t328 = (t152 * t205 + t227) * t119 + t153 * t31;
t105 = t119 * qJD(5);
t73 = t83 * qJ(5);
t327 = -t105 - t73;
t30 = qJD(1) * t329 - qJD(4) * t243 - t149 * qJDD(3) - t152 * t238;
t293 = t153 * t30;
t326 = -t87 * t250 - t293;
t116 = t119 ^ 2;
t82 = t87 ^ 2;
t325 = -t82 - t116;
t171 = t149 * t205 - t153 * t243;
t324 = t171 * t119;
t247 = qJD(4) * t149;
t253 = qJD(1) * t152;
t284 = qJ(5) * t150;
t322 = pkin(4) * t247 - t149 * qJD(5) - t253 * t284 - t93;
t321 = -pkin(5) * t31 + qJDD(6);
t85 = t149 * t252 - t243;
t319 = t85 ^ 2;
t318 = 0.2e1 * t73;
t317 = pkin(5) + pkin(8);
t315 = pkin(5) * t85;
t314 = pkin(8) * t83;
t313 = t30 * pkin(5);
t312 = t83 * pkin(4);
t311 = pkin(4) * t152;
t151 = sin(qJ(1));
t310 = g(1) * t151;
t154 = cos(qJ(1));
t145 = g(2) * t154;
t309 = g(3) * t150;
t308 = t153 * pkin(5);
t307 = t87 * t85;
t220 = t306 * t150;
t283 = qJ(6) * t149;
t305 = qJD(1) * t149 * t220 - t152 * qJD(6) + (-qJ(5) * t152 + t283) * qJD(4) + t322;
t28 = t149 * t67 + t152 * t74;
t246 = qJD(4) * t152;
t304 = pkin(4) * t149 * t254 - qJ(5) * t246 + t322;
t268 = t153 * t114;
t200 = pkin(3) * t153 + pkin(8) * t150;
t90 = t200 * qJD(1);
t303 = t149 * t90 + t152 * t268;
t278 = t149 * t150;
t302 = -t317 * t247 - (pkin(5) * t278 + qJ(5) * t153) * qJD(1) - t303;
t104 = t317 * t152;
t275 = t150 * t152;
t180 = -pkin(5) * t275 - t306 * t153;
t79 = t149 * t268;
t212 = -t152 * t90 + t79;
t301 = -qJD(1) * t180 + qJD(4) * t104 - t212;
t300 = pkin(8) * qJD(4);
t299 = qJ(5) * t31;
t10 = -t306 * t119 + t264;
t298 = t10 * t119;
t21 = -qJ(5) * t119 - t28;
t297 = t119 * t21;
t296 = t119 * t28;
t295 = t149 * t83;
t294 = t152 * t83;
t291 = t30 * t149;
t290 = t85 * qJ(5);
t289 = t85 * t119;
t288 = t87 * t119;
t273 = t150 * t155;
t286 = t149 * t323 + t152 * t273;
t285 = pkin(1) * qJDD(1);
t282 = qJD(3) * t85;
t281 = qJD(3) * t87;
t277 = t149 * t153;
t276 = t150 * t151;
t274 = t150 * t154;
t272 = t151 * t149;
t271 = t151 * t152;
t270 = t151 * t153;
t269 = t152 * t153;
t267 = t153 * t154;
t266 = t154 * t152;
t156 = qJD(3) ^ 2;
t265 = t155 * t156;
t263 = qJD(5) + t27;
t17 = t28 - t315;
t262 = -qJD(6) - t17;
t235 = g(2) * t267;
t261 = g(3) * t278 + t149 * t235;
t260 = g(3) * t275 + t152 * t235;
t259 = pkin(4) * t277 - qJ(5) * t269;
t258 = g(1) * t267 + g(2) * t270;
t257 = t154 * pkin(1) + t151 * qJ(2);
t147 = t153 ^ 2;
t256 = t150 ^ 2 - t147;
t157 = qJD(1) ^ 2;
t255 = -t156 - t157;
t248 = qJD(3) * t155;
t244 = qJD(4) * t155;
t241 = qJDD(1) * qJ(2);
t240 = qJDD(3) * t150;
t224 = t153 * t248;
t84 = qJD(3) * t200 + qJD(2);
t237 = t149 * t84 + t152 * t224 + t246 * t323;
t236 = g(1) * t270;
t130 = g(2) * t274;
t234 = 0.2e1 * qJD(1) * qJD(2);
t75 = -qJD(3) * pkin(3) - t268;
t177 = -t87 * qJ(5) + t75;
t13 = t306 * t85 + t177;
t233 = t13 * t247;
t232 = t13 * t246;
t231 = t151 * t269;
t222 = t149 * t244;
t221 = t150 * t246;
t216 = -t145 + t310;
t215 = g(3) * t153 - t130;
t68 = t150 * t272 - t266;
t69 = t149 * t154 + t150 * t271;
t214 = -t68 * pkin(4) + qJ(5) * t69;
t70 = t149 * t274 + t271;
t71 = t150 * t266 - t272;
t213 = t70 * pkin(4) - qJ(5) * t71;
t37 = qJD(1) * t84 + qJDD(1) * t323;
t107 = t155 * qJDD(1) + qJDD(2);
t43 = qJDD(3) * pkin(8) + t107 * t150 + t114 * t249;
t211 = -t149 * t37 - t152 * t43 - t67 * t246 + t74 * t247;
t210 = -t149 * t43 + t152 * t37 - t74 * t246 - t67 * t247;
t111 = t149 * t273;
t209 = t152 * t323 - t111;
t207 = pkin(4) * t231 + pkin(8) * t276 + t330 * t270;
t206 = qJDD(2) - t285;
t204 = t329 * qJ(5) + t150 * t248 + t245 * t311;
t42 = -qJDD(3) * pkin(3) - t153 * t107 + t114 * t250;
t167 = t30 * qJ(5) - t87 * qJD(5) + t42;
t2 = t85 * qJD(6) + t306 * t31 + t167;
t203 = -t2 - t236;
t202 = g(1) * t70 + g(2) * t68;
t201 = -g(1) * t71 - g(2) * t69;
t46 = -t284 - t286;
t199 = g(1) * t154 + g(2) * t151;
t198 = -qJDD(5) + t210;
t197 = qJ(2) * t157 + t310;
t12 = qJD(6) - t21 - t315;
t196 = t10 * t152 - t12 * t149;
t195 = t10 * t149 + t12 * t152;
t20 = -pkin(4) * t119 + t263;
t194 = t149 * t21 + t152 * t20;
t193 = t149 * t20 - t152 * t21;
t192 = -t149 * t224 + t152 * t84 - t155 * t221 - t247 * t323;
t4 = t211 + t327;
t190 = t330 + t311;
t188 = t234 + 0.2e1 * t241;
t187 = -t107 + t197;
t185 = t119 * t246 + t295;
t184 = -t119 * t247 + t294;
t183 = -g(1) * t276 - t215;
t182 = 0.2e1 * qJ(2) * t242 + qJDD(3) * t155;
t181 = t119 * t300 + t236;
t78 = -t306 * t152 - t330;
t179 = pkin(3) * t276 + t69 * pkin(4) + t154 * pkin(7) + t68 * qJ(5) + t257;
t6 = pkin(4) * t31 + t167;
t176 = t181 + t6;
t175 = t119 * t75 - t314;
t24 = pkin(4) * t85 + t177;
t174 = -t119 * t24 + t314;
t141 = t154 * qJ(2);
t173 = pkin(3) * t274 + t71 * pkin(4) - pkin(8) * t267 + t70 * qJ(5) + t141;
t172 = g(1) * t68 - g(2) * t70 + g(3) * t277 + t210;
t170 = -qJDD(5) + t172;
t169 = t188 - t199;
t168 = t83 - t307;
t166 = g(1) * t69 - g(2) * t71 + g(3) * t269 + t211;
t5 = -t198 - t312;
t165 = qJD(4) * t194 + t5 * t149 - t4 * t152;
t164 = t24 * t87 - t170;
t163 = t13 * t87 - t170 - t313;
t15 = t289 - t30;
t162 = -t13 * t85 - t166 + t321;
t161 = t288 - t31;
t160 = t85 * t250 - t83 * t278 - t328;
t159 = t83 * t275 - t324 + t326;
t158 = t171 * t85 - t275 * t31 - t30 * t278 + (t221 + t227 + t253) * t87;
t138 = qJDD(3) * t153;
t103 = t317 * t149;
t48 = -t153 * t155 + t259;
t47 = -t150 * pkin(4) - t209;
t45 = (-t155 + t283) * t153 + t259;
t40 = pkin(4) * t87 + t290;
t36 = -pkin(5) * t277 - t46;
t35 = -pkin(4) * t252 + t212;
t34 = -qJ(5) * t252 - t303;
t26 = t111 + (-t323 + t308) * t152 - t220;
t22 = t306 * t87 + t290;
t19 = -pkin(4) * t228 - qJD(5) * t269 + t204;
t14 = -pkin(4) * t249 - t192;
t11 = -qJ(5) * t249 + (-qJD(5) + t222) * t150 - t237;
t9 = (qJ(6) * qJD(4) - qJD(5)) * t269 + (-qJD(3) * t220 + qJD(6) * t153) * t149 + t204;
t8 = (-pkin(5) * t246 + qJ(5) * qJD(3)) * t153 + (qJD(5) + (pkin(5) * qJD(3) - t244) * t149) * t150 + t237;
t7 = -pkin(5) * t223 + qJD(3) * t180 - t150 * qJD(6) - t192;
t3 = -t4 + t321;
t1 = -t119 * qJD(6) - t198 - t229 - t313;
t16 = [qJDD(1), t216, t199, qJDD(2) - t216 - 0.2e1 * t285, t169, -t206 * pkin(1) - g(1) * (-t151 * pkin(1) + t141) - g(2) * t257 + (t234 + t241) * qJ(2), qJDD(1) * t147 - 0.2e1 * t150 * t219, -0.2e1 * t150 * t238 + 0.2e1 * t242 * t256, -t150 * t156 + t138, -t153 * t156 - t240, 0, t182 * t153 + (t169 - t265) * t150, -t182 * t150 + (t188 - t265) * t153 - t258, t326 * t152 - t87 * t223 (t149 * t87 + t152 * t85) * t250 + (t291 - t152 * t31 + (t149 * t85 - t152 * t87) * qJD(4)) * t153 (-t119 * t243 - t30) * t150 + (t184 + t281) * t153 (t119 * t251 - t31) * t150 + (-t185 - t282) * t153, t119 * t249 + t150 * t83, t192 * t119 + t209 * t83 + ((-t149 * t75 + t155 * t85) * qJD(3) + t210) * t150 + (-qJD(3) * t27 + t42 * t149 - t155 * t31 + t246 * t75) * t153 + t201, -t237 * t119 - t286 * t83 + (t119 * t222 + (-t152 * t75 + t155 * t87) * qJD(3) + t211) * t150 + (-qJD(3) * t28 + t42 * t152 + t155 * t30 - t247 * t75) * t153 + t202, t11 * t85 + t14 * t87 - t47 * t30 + t46 * t31 - t194 * t250 + (-qJD(4) * t193 + t149 * t4 + t152 * t5) * t153 + t258, t14 * t119 - t19 * t85 - t48 * t31 + t47 * t83 + (t24 * t251 + t5) * t150 + (qJD(3) * t20 - t6 * t149 - t24 * t246) * t153 - t201, -t11 * t119 - t19 * t87 + t48 * t30 - t46 * t83 + (t24 * t243 - t4) * t150 + (-qJD(3) * t21 - t6 * t152 + t24 * t247) * t153 - t202, t6 * t48 + t24 * t19 + t4 * t46 + t21 * t11 + t5 * t47 + t20 * t14 - g(1) * (t155 * t151 + t173) - g(2) * (-pkin(8) * t270 + t179) -t26 * t30 - t36 * t31 + t7 * t87 - t8 * t85 - t196 * t250 + (-qJD(4) * t195 + t1 * t152 - t149 * t3) * t153 + t258, t8 * t119 + t45 * t30 + t36 * t83 - t9 * t87 + (t13 * t243 + t3) * t150 + (qJD(3) * t12 - t2 * t152 + t233) * t153 - t202, -t7 * t119 - t26 * t83 + t45 * t31 + t9 * t85 + (-t13 * t251 - t1) * t150 + (-qJD(3) * t10 + t2 * t149 + t232) * t153 + t201, t2 * t45 + t13 * t9 + t1 * t26 + t10 * t7 + t3 * t36 + t12 * t8 - g(1) * (-pkin(5) * t267 + t71 * qJ(6) + t173) - g(2) * (t69 * qJ(6) + t179) + (g(2) * t317 * t153 - g(1) * t155) * t151; 0, 0, 0, qJDD(1), -t157, t145 - t197 + t206, 0, 0, 0, 0, 0, t150 * t255 + t138, t153 * t255 - t240, 0, 0, 0, 0, 0, t160, t293 + (t281 - t294) * t150 + t324, t158 (-t282 + t295) * t150 + t328, t159, t194 * qJD(1) + (qJD(3) * t193 - t6) * t153 + (qJD(3) * t24 + t165) * t150 - t216, t158, t159, t160, t196 * qJD(1) + (qJD(3) * t195 - t2) * t153 + (qJD(3) * t13 + qJD(4) * t196 + t1 * t149 + t3 * t152) * t150 - t216; 0, 0, 0, 0, 0, 0, t153 * t157 * t150, -t256 * t157, t238, -t239, qJDD(3), t309 + (-t187 + t145) * t153, t150 * t187 + t215, t152 * t288 - t291 (-t30 - t289) * t152 + (-t31 - t288) * t149 (t119 * t275 - t153 * t87) * qJD(1) + t185 (-t119 * t278 + t153 * t85) * qJD(1) + t184, -t119 * t252, t27 * t252 - t85 * t93 - pkin(3) * t31 + t79 * t119 + (-t236 - t42 + (-t90 - t300) * t119) * t152 + t175 * t149 + t260, pkin(3) * t30 - t87 * t93 + t303 * t119 + t28 * t252 + t175 * t152 + (t181 + t42) * t149 - t261, -t34 * t85 - t35 * t87 + (-t4 + t119 * t20 + (-t31 + t280) * pkin(8)) * t152 + (t5 + t297 + (qJD(4) * t85 - t30) * pkin(8)) * t149 + t183, -t35 * t119 + t149 * t174 + t152 * t176 + t190 * t31 - t20 * t252 - t304 * t85 - t260, t34 * t119 - t149 * t176 + t152 * t174 - t190 * t30 + t21 * t252 - t304 * t87 + t261, -t21 * t34 - t20 * t35 - g(1) * t207 - g(3) * t142 + t304 * t24 + (t165 + t130) * pkin(8) + (-t6 + t309 + t235) * t190, -t103 * t30 - t104 * t31 + t301 * t87 - t302 * t85 + (t3 + t298) * t152 + (-t12 * t119 + t1) * t149 + t183, -t232 + t104 * t83 + t78 * t30 - t305 * t87 + t203 * t149 + t302 * t119 + (-t12 * t153 - t13 * t275) * qJD(1) + t261, t233 - t103 * t83 + t78 * t31 + t305 * t85 + t203 * t152 - t301 * t119 + (t10 * t153 + t13 * t278) * qJD(1) + t260, t2 * t78 + t1 * t103 + t3 * t104 - g(1) * (qJ(6) * t231 + t207) - g(3) * (t142 + t308) + t305 * t13 + t302 * t12 + t301 * t10 + (-pkin(5) * t310 - g(3) * t78) * t150 - (-t317 * t150 + t78 * t153) * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t307, t82 - t319, t15, t161, t83, -t75 * t87 + t172 + t296, -t119 * t27 + t75 * t85 + t166, pkin(4) * t30 - t299 + (-t21 - t28) * t87 + (t20 - t263) * t85, t40 * t85 + t164 - t296 - 0.2e1 * t312, t119 * t263 - t24 * t85 + t40 * t87 + t105 - t166 + t318, -t5 * pkin(4) - g(1) * t214 - g(2) * t213 + g(3) * t259 - t4 * qJ(5) - t20 * t28 - t21 * t263 - t24 * t40, -t299 + t306 * t30 + (t12 + t262) * t87 + (t10 - t264) * t85, t119 * t191 + t22 * t87 + 0.2e1 * t105 + t162 + t318, -t22 * t85 + (0.2e1 * qJD(6) + t17) * t119 + 0.2e1 * t229 - t163, -t1 * t306 + t3 * qJ(5) - t13 * t22 - g(1) * (-qJ(6) * t68 + t214) - g(2) * (qJ(6) * t70 + t213) - g(3) * (-qJ(6) * t277 - t259) + t264 * t12 + t262 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t168, t325, t164 + t297 - t312, t15, t325, -t168, -t229 + (-qJD(6) - t12) * t119 + t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, t83 + t307, -t116 - t319, t162 + t298 - t327;];
tau_reg  = t16;
