% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR15
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR15_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:14
% EndTime: 2019-12-31 20:43:23
% DurationCPUTime: 4.77s
% Computational Cost: add. (5529->525), mult. (11476->660), div. (0->0), fcn. (7036->10), ass. (0->266)
t332 = pkin(3) + pkin(6);
t184 = sin(qJ(5));
t185 = sin(qJ(4));
t188 = cos(qJ(5));
t189 = cos(qJ(4));
t116 = t184 * t189 + t185 * t188;
t186 = sin(qJ(2));
t212 = t116 * t186;
t338 = qJD(4) + qJD(5);
t60 = t338 * t116;
t317 = qJD(1) * t212 + t60;
t276 = qJD(2) * t185;
t190 = cos(qJ(2));
t279 = qJD(1) * t190;
t113 = t189 * t279 + t276;
t250 = t185 * t279;
t274 = qJD(2) * t189;
t115 = -t250 + t274;
t225 = t113 * t184 - t115 * t188;
t54 = t113 * t188 + t115 * t184;
t322 = t54 * t225;
t298 = t185 * t186;
t220 = pkin(4) * t190 - pkin(8) * t298;
t271 = qJD(4) * t185;
t333 = pkin(2) + pkin(7);
t320 = pkin(8) + t333;
t162 = pkin(6) * t279;
t123 = pkin(3) * t279 + t162;
t280 = qJD(1) * t186;
t165 = pkin(2) * t280;
t305 = qJ(3) * t190;
t226 = pkin(7) * t186 - t305;
t85 = qJD(1) * t226 + t165;
t47 = t123 * t189 - t185 * t85;
t348 = -qJD(1) * t220 + t320 * t271 - t47;
t127 = t320 * t189;
t254 = t189 * t280;
t48 = t123 * t185 + t189 * t85;
t347 = pkin(8) * t254 + qJD(4) * t127 + t48;
t346 = t225 ^ 2 - t54 ^ 2;
t150 = qJD(4) + t280;
t138 = qJD(5) + t150;
t269 = qJD(4) * t190;
t253 = t185 * t269;
t208 = t186 * t274 + t253;
t263 = t190 * qJDD(1);
t270 = qJD(4) * t189;
t257 = qJD(2) * t270 + qJDD(2) * t185 + t189 * t263;
t199 = qJD(1) * t208 - t257;
t267 = qJD(5) * t188;
t268 = qJD(5) * t184;
t265 = qJD(1) * qJD(2);
t248 = t186 * t265;
t218 = -t189 * qJDD(2) + (-t248 + t263) * t185;
t49 = qJD(4) * t113 + t218;
t14 = t113 * t267 + t115 * t268 - t184 * t199 + t188 * t49;
t345 = t138 * t54 - t14;
t170 = t186 * qJ(3);
t245 = -pkin(1) - t170;
t210 = -t190 * t333 + t245;
t73 = t210 * qJD(1);
t161 = pkin(6) * t280;
t339 = qJD(3) + t161;
t266 = pkin(3) * t280 + t339;
t79 = -qJD(2) * t333 + t266;
t36 = -t185 * t73 + t189 * t79;
t29 = -pkin(8) * t115 + t36;
t25 = pkin(4) * t150 + t29;
t37 = t185 * t79 + t189 * t73;
t30 = -pkin(8) * t113 + t37;
t247 = t190 * t265;
t264 = t186 * qJDD(1);
t211 = t247 + t264;
t112 = qJDD(4) + t211;
t149 = pkin(2) * t248;
t272 = qJD(3) * t186;
t203 = qJD(2) * t226 - t272;
t35 = qJD(1) * t203 + qJDD(1) * t210 + t149;
t148 = pkin(6) * t247;
t158 = pkin(6) * t264;
t246 = qJDD(3) + t148 + t158;
t51 = pkin(3) * t211 - qJDD(2) * t333 + t246;
t9 = -qJD(4) * t37 - t185 * t35 + t189 * t51;
t6 = pkin(4) * t112 + pkin(8) * t49 + t9;
t262 = t185 * t51 + t189 * t35 + t270 * t79;
t8 = -t271 * t73 + t262;
t7 = pkin(8) * t199 + t8;
t1 = (qJD(5) * t25 + t7) * t188 + t184 * t6 - t30 * t268;
t183 = qJ(4) + qJ(5);
t168 = sin(t183);
t177 = g(3) * t190;
t180 = qJD(2) * qJ(3);
t95 = t180 + t123;
t62 = pkin(4) * t113 + t95;
t169 = cos(t183);
t187 = sin(qJ(1));
t191 = cos(qJ(1));
t294 = t186 * t191;
t81 = t168 * t294 + t169 * t187;
t296 = t186 * t187;
t83 = -t168 * t296 + t169 * t191;
t344 = g(1) * t81 - g(2) * t83 - t168 * t177 + t54 * t62 - t1;
t309 = t188 * t30;
t11 = t184 * t25 + t309;
t2 = -qJD(5) * t11 - t184 * t7 + t188 * t6;
t80 = -t168 * t187 + t169 * t294;
t82 = t168 * t191 + t169 * t296;
t343 = -g(1) * t80 - g(2) * t82 + t169 * t177 + t225 * t62 + t2;
t202 = qJD(5) * t225 + t184 * t49 + t188 * t199;
t342 = -t138 * t225 + t202;
t341 = t150 * t37 + t9;
t275 = qJD(2) * t186;
t122 = t332 * t275;
t159 = pkin(6) * t263;
t178 = qJDD(2) * qJ(3);
t179 = qJD(2) * qJD(3);
t256 = t159 + t178 + t179;
t234 = pkin(3) * t263 + t256;
t52 = -qJD(1) * t122 + t234;
t340 = qJD(4) * t150 * t333 + t52;
t290 = t188 * t189;
t300 = t184 * t185;
t222 = -t290 + t300;
t86 = t222 * t190;
t227 = g(1) * t191 + g(2) * t187;
t289 = t189 * t190;
t288 = t189 * t191;
t97 = -t185 * t187 + t186 * t288;
t292 = t187 * t189;
t99 = t185 * t191 + t186 * t292;
t337 = -g(1) * t97 - g(2) * t99 + g(3) * t289;
t336 = t14 * t222 + t225 * t317;
t316 = -t184 * t271 - t185 * t268 + t188 * t254 - t280 * t300 + t290 * t338;
t335 = -t116 * t202 + t316 * t54;
t311 = t184 * t30;
t10 = t188 * t25 - t311;
t255 = -g(1) * t294 - g(2) * t296 + t177;
t334 = t1 * t116 - t10 * t317 + t11 * t316 - t2 * t222 + t255;
t329 = pkin(4) * t185;
t328 = pkin(7) * t190;
t327 = g(1) * t187;
t324 = g(2) * t191;
t323 = g(3) * t186;
t174 = t190 * pkin(2);
t126 = t320 * t185;
t63 = t126 * t184 - t127 * t188;
t319 = qJD(5) * t63 + t184 * t348 - t188 * t347;
t64 = -t126 * t188 - t127 * t184;
t318 = -qJD(5) * t64 + t184 * t347 + t188 * t348;
t313 = t150 * t36;
t310 = t185 * t36;
t308 = t189 * t49;
t157 = pkin(4) * t189 + pkin(3);
t307 = pkin(4) * t270 + t157 * t280 + t339;
t284 = t174 + t170;
t231 = t284 + t328;
t106 = -pkin(1) - t231;
t134 = t332 * t186;
t118 = t185 * t134;
t59 = t106 * t189 + t118;
t306 = pkin(6) * qJDD(2);
t304 = qJDD(2) * pkin(2);
t303 = t112 * t185;
t302 = t113 * t189;
t301 = t115 * t113;
t299 = t185 * t113;
t297 = t185 * t190;
t295 = t186 * t189;
t195 = qJD(1) ^ 2;
t293 = t186 * t195;
t291 = t187 * t190;
t91 = t189 * t112;
t287 = t190 * t191;
t192 = -pkin(8) - pkin(7);
t286 = t190 * t192;
t285 = t333 * t112;
t135 = t332 * t190;
t283 = pkin(1) * t191 + pkin(6) * t187;
t181 = t186 ^ 2;
t182 = t190 ^ 2;
t282 = t181 - t182;
t281 = t181 + t182;
t278 = qJD(2) * t113;
t277 = qJD(2) * t115;
t273 = qJD(2) * t190;
t261 = pkin(4) * t297;
t259 = t113 * t295;
t258 = t185 * t294;
t252 = t189 * t269;
t244 = pkin(8) * t190 - t106;
t241 = qJD(1) * t59 + t37;
t124 = t332 * t273;
t164 = pkin(2) * t275;
t70 = t164 + t203;
t240 = t124 * t189 - t185 * t70;
t239 = -qJD(2) * pkin(2) + qJD(3);
t238 = -qJD(1) * t135 - t95;
t237 = t150 + t280;
t236 = pkin(2) * t287 + qJ(3) * t294 + t283;
t235 = -t158 - t255;
t105 = qJDD(5) + t112;
t233 = -t222 * t105 - t138 * t317;
t232 = t186 * t247;
t230 = t257 * t185;
t229 = t281 * qJDD(1) * pkin(6);
t194 = qJD(2) ^ 2;
t228 = pkin(6) * t194 + t324;
t119 = t189 * t134;
t39 = pkin(4) * t186 + t185 * t244 + t119;
t44 = -pkin(8) * t289 + t59;
t20 = -t184 * t44 + t188 * t39;
t21 = t184 * t39 + t188 * t44;
t224 = -t115 * t189 + t299;
t125 = t161 + t239;
t133 = -t162 - t180;
t223 = t125 * t190 + t133 * t186;
t221 = t150 * t185;
t219 = t245 - t174;
t96 = t219 * qJD(1);
t217 = t280 * t96 + qJDD(3) - t235;
t216 = -0.2e1 * pkin(1) * t265 - t306;
t215 = -t150 * t270 - t303;
t214 = -t105 * t116 - t138 * t316;
t213 = -qJ(3) * t273 - t272;
t22 = -t106 * t271 + t124 * t185 + t134 * t270 + t189 * t70;
t209 = 0.2e1 * qJDD(1) * pkin(1) - t228;
t206 = t189 * t248 - t257;
t128 = -pkin(1) - t284;
t205 = t306 + (-qJD(1) * t128 - t96) * qJD(2);
t204 = -t190 * t227 - t323;
t50 = qJD(1) * t213 + qJDD(1) * t219 + t149;
t90 = t164 + t213;
t201 = qJD(1) * t90 + qJDD(1) * t128 + t228 + t50;
t65 = -pkin(4) * t253 + (-pkin(6) - t157) * t275;
t76 = pkin(6) * t248 - t256;
t84 = t246 - t304;
t200 = qJD(2) * t223 + t84 * t186 - t76 * t190;
t175 = t191 * pkin(6);
t155 = qJ(3) + t329;
t153 = g(1) * t291;
t147 = qJ(3) * t287;
t145 = qJ(3) * t291;
t143 = t190 * t293;
t131 = t282 * t195;
t130 = qJDD(2) * t190 - t186 * t194;
t129 = qJDD(2) * t186 + t190 * t194;
t120 = -qJ(3) * t279 + t165;
t108 = qJDD(1) * t182 - 0.2e1 * t232;
t107 = qJDD(1) * t181 + 0.2e1 * t232;
t100 = -t185 * t296 + t288;
t98 = t258 + t292;
t94 = pkin(4) * t289 + t135;
t87 = t116 * t190;
t69 = 0.2e1 * t186 * t263 - 0.2e1 * t265 * t282;
t58 = -t106 * t185 + t119;
t32 = t190 * t60 - t222 * t275;
t31 = qJD(2) * t212 + t338 * t86;
t24 = pkin(4) * t257 + qJD(1) * t65 + t234;
t23 = -qJD(4) * t59 + t240;
t19 = pkin(8) * t208 + t22;
t18 = t220 * qJD(2) + (t189 * t244 - t118) * qJD(4) + t240;
t13 = t188 * t29 - t311;
t12 = -t184 * t29 - t309;
t4 = -qJD(5) * t21 + t18 * t188 - t184 * t19;
t3 = qJD(5) * t20 + t18 * t184 + t188 * t19;
t5 = [0, 0, 0, 0, 0, qJDD(1), -t324 + t327, t227, 0, 0, t107, t69, t129, t108, t130, 0, t186 * t216 + t190 * t209 + t153, t216 * t190 + (-t209 - t327) * t186, -t227 + 0.2e1 * t229, -g(1) * (-pkin(1) * t187 + t175) - g(2) * t283 + (pkin(6) ^ 2 * t281 + pkin(1) ^ 2) * qJDD(1), 0, -t129, -t130, t107, t69, t108, t229 + t200 - t227, t186 * t205 + t190 * t201 - t153, t205 * t190 + (-t201 + t327) * t186, pkin(6) * t200 - g(1) * t175 - g(2) * t236 + t50 * t128 - t219 * t327 + t96 * t90, t49 * t297 + (t185 * t275 - t252) * t115, -t224 * t275 + (-t206 * t185 + t308 + (t302 + (t115 - t250) * t185) * qJD(4)) * t190, (t150 * t276 - t49) * t186 + (t215 + t277) * t190, -t113 * t208 - t199 * t289, (t237 * t274 - t257) * t186 + (t237 * t271 - t278 - t91) * t190, t112 * t186 + t150 * t273, t23 * t150 + t58 * t112 - t122 * t113 + t135 * t257 - g(1) * t100 - g(2) * t98 + (t238 * t274 + t9) * t186 + (t36 * qJD(2) + t52 * t189 + t238 * t271) * t190, g(1) * t99 - g(2) * t97 - t112 * t59 - t115 * t122 - t135 * t49 - t150 * t22 + (t276 * t95 - t8) * t186 + (-qJD(2) * t37 - t52 * t185 - t270 * t95) * t190, -t22 * t113 - t59 * t257 - t23 * t115 + t58 * t49 + t153 + (t189 * t241 - t310) * t275 + (-t324 + t9 * t185 - t8 * t189 + (t185 * t241 + t189 * t36) * qJD(4)) * t190, t8 * t59 + t37 * t22 + t9 * t58 + t36 * t23 + t52 * t135 - t95 * t122 - g(1) * (pkin(3) * t191 + t175) - g(2) * (pkin(7) * t287 + t236) + (-g(1) * (t219 - t328) - g(2) * pkin(3)) * t187, t14 * t87 - t225 * t31, -t14 * t86 - t202 * t87 - t225 * t32 - t31 * t54, -t105 * t87 + t138 * t31 - t14 * t186 - t225 * t273, t202 * t86 - t32 * t54, t105 * t86 + t138 * t32 + t186 * t202 - t273 * t54, t105 * t186 + t138 * t273, -g(1) * t83 - g(2) * t81 + t10 * t273 + t105 * t20 + t138 * t4 + t186 * t2 - t202 * t94 - t24 * t86 - t32 * t62 + t54 * t65, g(1) * t82 - g(2) * t80 - t1 * t186 - t105 * t21 - t11 * t273 - t138 * t3 - t14 * t94 - t225 * t65 - t24 * t87 + t31 * t62, -g(2) * t287 + t1 * t86 - t10 * t31 + t11 * t32 + t14 * t20 + t2 * t87 + t202 * t21 + t225 * t4 - t3 * t54 + t153, t1 * t21 + t11 * t3 + t2 * t20 + t10 * t4 + t24 * t94 + t62 * t65 - g(1) * (t157 * t191 + t175) - g(2) * (pkin(4) * t258 - t191 * t286 + t236) + (-g(1) * (-pkin(4) * t298 + t219 + t286) - g(2) * t157) * t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t131, t264, t143, t263, qJDD(2), pkin(1) * t293 + t235, t323 - t159 + (pkin(1) * t195 + t227) * t190, 0, 0, qJDD(2), -t264, -t263, -t143, t131, t143, (-pkin(2) * t186 + t305) * qJDD(1) + ((-t133 - t180) * t186 + (-t125 + t239) * t190) * qJD(1), -t120 * t279 + t217 - 0.2e1 * t304, t159 + 0.2e1 * t178 + 0.2e1 * t179 + (qJD(1) * t120 - g(3)) * t186 + (qJD(1) * t96 - t227) * t190, -t76 * qJ(3) - t133 * qJD(3) - t84 * pkin(2) - t96 * t120 - g(1) * (-pkin(2) * t294 + t147) - g(2) * (-pkin(2) * t296 + t145) - g(3) * t284 - t223 * qJD(1) * pkin(6), -t115 * t221 - t308, -t257 * t189 + t49 * t185 + t224 * qJD(4) + (t185 * t252 + (t299 + (-t115 + t274) * t189) * t186) * qJD(1), -t150 * t271 + t91 + (-t115 * t190 - t150 * t298) * qJD(1), t230 + t113 * t270 + (-t185 * t208 + t259) * qJD(1), (t113 * t190 - t150 * t295) * qJD(1) + t215, -t150 * t279, qJ(3) * t257 - t47 * t150 - t36 * t279 + t266 * t113 + (t95 * qJD(4) - t285 + (t95 - t180) * t280) * t189 + (-t323 + (-qJ(3) * qJD(1) * qJD(4) - t227) * t190 + t340) * t185, t37 * t279 - qJ(3) * t49 + t150 * t48 + t266 * t115 + (-t150 * t95 + t285) * t185 + (t204 + t340) * t189, t48 * t113 + t47 * t115 + (-t37 * t280 - t333 * t49 - t9 + (t113 * t333 - t37) * qJD(4)) * t189 + (-t333 * t206 - t8 + t36 * t280 + (t36 - (t115 + t250) * t333) * qJD(4)) * t185 - t255, -g(1) * t147 - g(2) * t145 - g(3) * t231 + t52 * qJ(3) + t266 * t95 - t36 * t47 - t37 * t48 + (-t8 * t185 - t9 * t189 - (t189 * t37 - t310) * qJD(4) + t227 * t186) * t333, t336, t116 * t14 - t202 * t222 + t225 * t316 + t317 * t54, t225 * t279 + t233, t335, t279 * t54 + t214, -t138 * t279, -t10 * t279 + t105 * t63 + t116 * t24 + t138 * t318 - t155 * t202 + t168 * t204 + t307 * t54 + t316 * t62, -t105 * t64 + t11 * t279 - t138 * t319 - t14 * t155 + t169 * t204 - t222 * t24 - t225 * t307 - t317 * t62, t14 * t63 + t202 * t64 + t225 * t318 - t319 * t54 - t334, t1 * t64 + t2 * t63 + t24 * t155 - g(1) * (t191 * t261 + t147) - g(2) * (t187 * t261 + t145) - g(3) * (t284 - t286) + t307 * t62 + (-g(3) * t329 + t227 * (pkin(2) - t192)) * t186 + t319 * t11 + t318 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, qJDD(2) + t143, -t181 * t195 - t194, qJD(2) * t133 + t148 + t217 - t304, 0, 0, 0, 0, 0, 0, -t150 * t221 - t278 + t91, -t150 ^ 2 * t189 - t277 - t303, -t230 + t308 + (t115 * t185 - t302) * qJD(4) + (-t259 + (t253 + (t115 + t274) * t186) * t185) * qJD(1), -qJD(2) * t95 + t341 * t189 + (t8 - t313) * t185 + t255, 0, 0, 0, 0, 0, 0, -qJD(2) * t54 + t233, qJD(2) * t225 + t214, -t335 - t336, -qJD(2) * t62 + t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301, -t113 ^ 2 + t115 ^ 2, -t218 + (-qJD(4) + t150) * t113, -t301, t115 * t150 + t199, t112, -t115 * t95 + t337 + t341, g(1) * t98 - g(2) * t100 + t113 * t95 + t313 + (qJD(4) * t73 - t177) * t185 - t262, 0, 0, -t322, t346, t345, t322, t342, t105, -t12 * t138 + (t105 * t188 - t115 * t54 - t138 * t268) * pkin(4) + t343, t13 * t138 + (-t105 * t184 + t115 * t225 - t138 * t267) * pkin(4) + t344, -t10 * t54 - t11 * t225 - t12 * t225 + t13 * t54 + (t14 * t188 + t202 * t184 + (-t184 * t225 - t188 * t54) * qJD(5)) * pkin(4), -t10 * t12 - t11 * t13 + (t1 * t184 + t2 * t188 - t62 * t115 + (-t10 * t184 + t11 * t188) * qJD(5) + t337) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, t346, t345, t322, t342, t105, t11 * t138 + t343, t10 * t138 + t344, 0, 0;];
tau_reg = t5;
