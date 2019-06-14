% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPPRP5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPPRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:27:55
% EndTime: 2019-05-06 09:28:18
% DurationCPUTime: 7.90s
% Computational Cost: add. (19455->428), mult. (42983->520), div. (0->0), fcn. (26915->8), ass. (0->267)
t251 = sin(qJ(2));
t247 = sin(pkin(9));
t248 = cos(pkin(9));
t253 = cos(qJ(5));
t254 = cos(qJ(2));
t301 = qJD(1) * qJD(2);
t235 = t254 * t301;
t237 = t251 * qJDD(1);
t214 = t237 + t235;
t207 = qJDD(5) + t214;
t305 = qJD(1) * t254;
t204 = qJD(2) * t247 + t248 * t305;
t206 = qJD(2) * t248 - t247 * t305;
t250 = sin(qJ(5));
t169 = t204 * t253 + t206 * t250;
t171 = -t204 * t250 + t206 * t253;
t322 = t171 * t169;
t127 = -t322 - t207;
t329 = t127 * t250;
t168 = t171 ^ 2;
t306 = qJD(1) * t251;
t232 = qJD(5) + t306;
t352 = t232 ^ 2;
t366 = -t168 - t352;
t71 = t253 * t366 + t329;
t328 = t127 * t253;
t73 = -t250 * t366 + t328;
t46 = t247 * t73 + t248 * t71;
t416 = t251 * t46;
t346 = -pkin(2) - qJ(4);
t415 = t346 * t46;
t414 = -pkin(3) * t46 - pkin(4) * t71;
t313 = t251 * qJ(3);
t266 = t254 * t346 - pkin(1) - t313;
t413 = t266 * (t247 * t71 - t248 * t73);
t353 = t169 ^ 2;
t150 = t353 - t352;
t84 = -t150 * t250 + t328;
t88 = -t150 * t253 - t329;
t157 = t232 * t171;
t238 = t254 * qJDD(1);
t293 = t251 * t301;
t215 = t238 - t293;
t182 = -qJDD(2) * t247 - t215 * t248;
t183 = qJDD(2) * t248 - t215 * t247;
t287 = -t182 * t253 + t183 * t250;
t269 = qJD(5) * t171 + t287;
t93 = -t157 + t269;
t412 = -t251 * t93 + t254 * (t247 * t88 + t248 * t84);
t410 = pkin(8) * t71;
t409 = pkin(8) * t73;
t135 = t168 - t353;
t275 = t182 * t250 + t183 * t253;
t116 = -qJD(5) * t169 + t275;
t323 = t169 * t232;
t368 = t116 - t323;
t92 = t157 + t269;
t58 = -t250 * t92 + t253 * t368;
t339 = t250 * t368;
t60 = t253 * t92 + t339;
t408 = t251 * t135 + t254 * (t247 * t60 - t248 * t58);
t405 = t247 * t84 - t248 * t88;
t403 = t247 * t58 + t248 * t60;
t362 = -t322 + t207;
t327 = t362 * t250;
t361 = -t352 - t353;
t370 = t253 * t361 - t327;
t117 = t253 * t362;
t371 = t250 * t361 + t117;
t381 = t247 * t370 + t248 * t371;
t400 = t251 * t381;
t399 = t346 * t381;
t398 = pkin(3) * t381 + pkin(4) * t371;
t152 = -t168 + t352;
t382 = t152 * t253 + t327;
t383 = -t152 * t250 + t117;
t397 = -t247 * t382 + t248 * t383;
t367 = t116 + t323;
t396 = t251 * t367 + t254 * (-t247 * t383 - t248 * t382);
t395 = 2 * qJD(3);
t102 = -t353 - t168;
t394 = pkin(3) * t102;
t393 = pkin(4) * t102;
t391 = pkin(8) * t370;
t390 = pkin(8) * t371;
t389 = qJ(3) * t102;
t388 = qJ(6) * t368;
t387 = t102 * t254;
t380 = t266 * (-t247 * t371 + t248 * t370);
t379 = -2 * qJD(4);
t244 = t251 ^ 2;
t257 = qJD(1) ^ 2;
t239 = t244 * t257;
t256 = qJD(2) ^ 2;
t227 = -t239 - t256;
t311 = t251 * t257;
t295 = t254 * t311;
t221 = -qJDD(2) + t295;
t310 = t254 * t221;
t378 = pkin(7) * (-t227 * t251 + t310);
t174 = t206 * t204;
t365 = -t174 + t214;
t377 = t247 * t365;
t376 = t248 * t365;
t252 = sin(qJ(1));
t255 = cos(qJ(1));
t284 = g(1) * t255 + g(2) * t252;
t334 = qJDD(1) * pkin(7);
t196 = -pkin(1) * t257 - t284 + t334;
t347 = t254 * pkin(2);
t279 = -t313 - t347;
t211 = t279 * qJD(1);
t272 = (qJD(1) * t211 + t196) * t254;
t372 = qJD(2) * t395 + t272;
t278 = t214 + t235;
t369 = t278 * qJ(3);
t188 = t204 * t306;
t161 = t183 + t188;
t222 = pkin(3) * t306 - qJD(2) * qJ(4);
t231 = pkin(2) * t293;
t292 = qJD(3) * t306;
t234 = -0.2e1 * t292;
t245 = t254 ^ 2;
t290 = g(1) * t252 - g(2) * t255;
t271 = -qJDD(1) * pkin(1) - t290;
t105 = -t222 * t306 + t231 + t234 + (-pkin(3) * t245 - pkin(7)) * t257 + t346 * t215 - t369 + t271;
t320 = t196 * t251;
t267 = -qJDD(2) * pkin(2) - qJ(3) * t256 + t211 * t306 + qJDD(3) + t320;
t128 = pkin(3) * t214 - qJDD(2) * qJ(4) + (-pkin(3) * t301 - qJ(4) * t311 + g(3)) * t254 + t267;
t281 = -t105 * t247 + t128 * t248 + t206 * t379;
t65 = t105 * t248 + t128 * t247 + t204 * t379;
t314 = t245 * t257;
t364 = t310 - (-t256 + t314) * t251;
t360 = pkin(5) * t269 - t388;
t184 = pkin(4) * t306 - pkin(8) * t206;
t202 = t204 ^ 2;
t359 = -pkin(4) * t182 - pkin(8) * t202 + t184 * t206;
t216 = t238 - 0.2e1 * t293;
t228 = t256 + t314;
t220 = qJDD(2) + t295;
t319 = t220 * t251;
t358 = pkin(7) * (t228 * t254 + t319) - pkin(1) * t216;
t268 = (-t169 * t250 - t171 * t253) * t232;
t316 = t232 * t250;
t149 = t171 * t316;
t315 = t232 * t253;
t298 = t169 * t315;
t282 = t149 - t298;
t357 = -t247 * t268 + t248 * t282;
t270 = t250 * t269 + t298;
t283 = t169 * t316 - t253 * t269;
t356 = -t247 * t283 + t248 * t270;
t355 = t251 * t207 + t254 * (-t247 * t282 - t248 * t268);
t297 = t251 * t322;
t354 = -t297 + t254 * (-t247 * t270 - t248 * t283);
t203 = t206 ^ 2;
t351 = 2 * qJD(6);
t349 = pkin(5) * t253;
t348 = g(3) * t254;
t241 = t251 * g(3);
t134 = pkin(5) * t169 - qJ(6) * t171;
t259 = pkin(4) * t365 - pkin(8) * t161 + t281;
t55 = -pkin(4) * t202 + pkin(8) * t182 - t184 * t306 + t65;
t30 = t250 * t259 + t253 * t55;
t280 = qJ(6) * t207 - t134 * t169 + t232 * t351 + t30;
t23 = -pkin(5) * t352 + t280;
t29 = t250 * t55 - t253 * t259;
t25 = -pkin(5) * t207 - qJ(6) * t352 + t134 * t171 + qJDD(6) + t29;
t345 = -pkin(5) * t25 + qJ(6) * t23;
t344 = -pkin(5) * t367 - qJ(6) * t93;
t14 = t250 * t30 - t253 * t29;
t343 = t14 * t247;
t342 = t14 * t248;
t299 = qJDD(2) * qJ(3);
t307 = -pkin(2) * t256 - t241;
t265 = pkin(3) * t215 - qJ(4) * t314 + qJDD(4) + t299 + t307;
t122 = t272 + (t395 + t222) * qJD(2) + t265;
t70 = t122 + t359;
t341 = t250 * t70;
t337 = t253 * t70;
t335 = qJ(6) * t253;
t333 = t367 * t250;
t294 = t206 * t306;
t159 = t182 + t294;
t113 = t159 * t247 - t161 * t248;
t332 = t113 * t251;
t331 = t122 * t247;
t330 = t122 * t248;
t163 = t174 + t214;
t325 = t163 * t247;
t324 = t163 * t248;
t312 = t251 * t216;
t218 = t239 + t314;
t308 = pkin(1) * t218 + (t244 + t245) * t334;
t302 = qJD(5) + t232;
t296 = t251 * t174;
t291 = -qJ(6) * t250 - pkin(4);
t15 = t250 * t29 + t253 * t30;
t179 = t320 + t348;
t180 = t196 * t254 - t241;
t288 = t179 * t251 + t180 * t254;
t79 = t116 * t250 + t171 * t315;
t80 = t116 * t253 - t149;
t285 = t254 * (-t247 * t80 - t248 * t79) + t297;
t38 = t247 * t65 + t248 * t281;
t277 = -t247 * t281 + t248 * t65;
t274 = (-t239 + t256) * t254 + t319;
t195 = pkin(7) * t257 - t271;
t264 = -pkin(5) * t366 - qJ(6) * t127 + t23;
t147 = t267 + t348;
t263 = t307 + t372;
t262 = pkin(5) * t362 + qJ(6) * t361 - t25;
t261 = pkin(2) * t215 + t195 - t231;
t144 = t263 + t299;
t260 = t261 + 0.2e1 * t292;
t258 = -qJD(2) * t222 + t171 * t351 - t265 - t359 - t360 - t372;
t219 = t239 - t314;
t213 = t237 + 0.2e1 * t235;
t198 = t251 * t214;
t187 = -t203 - t239;
t186 = -t203 + t239;
t185 = t202 - t239;
t178 = t235 * t251 + t198;
t177 = (t215 - t293) * t254;
t173 = t213 * t254 + t312;
t172 = -t239 - t202;
t160 = t183 - t188;
t158 = -t182 + t294;
t154 = -t202 - t203;
t140 = t187 * t248 - t325;
t129 = t172 * t247 + t376;
t98 = -t169 * t302 + t275;
t95 = (-qJD(5) + t232) * t171 - t287;
t94 = t171 * t302 + t287;
t89 = t253 * t367;
t63 = t253 * t95 + t333;
t61 = -t253 * t93 + t333;
t59 = t250 * t95 - t89;
t57 = -t250 * t93 - t89;
t53 = -t247 * t79 + t248 * t80;
t45 = t337 - t410;
t44 = t341 - t390;
t37 = -pkin(4) * t98 + t341 + t409;
t36 = (pkin(5) * t232 - (2 * qJD(6))) * t171 + t70 + t360;
t35 = -pkin(4) * t92 - t337 + t391;
t32 = t247 * t63 + t248 * t59;
t31 = t247 * t61 + t248 * t57;
t27 = t258 + (-t94 - t157) * pkin(5);
t26 = -pkin(5) * t157 + t258 + t388;
t21 = -qJ(6) * t102 + t25;
t20 = (-t102 - t352) * pkin(5) + t280;
t19 = -t250 * t27 - t335 * t94 - t390;
t18 = -pkin(5) * t339 + t253 * t26 + t410;
t17 = t253 * t27 + t291 * t94 + t391;
t16 = -t409 + t250 * t26 + (pkin(4) + t349) * t368;
t13 = -pkin(4) * t70 + pkin(8) * t15;
t12 = -pkin(8) * t59 - t14;
t11 = t23 * t253 + t25 * t250;
t10 = t23 * t250 - t25 * t253;
t9 = pkin(8) * t63 + t15 - t393;
t8 = -pkin(8) * t57 - t20 * t250 + t21 * t253;
t7 = pkin(8) * t61 + t20 * t253 + t21 * t250 - t393;
t5 = t15 * t247 + t342;
t4 = -pkin(8) * t10 + (pkin(5) * t250 - t335) * t36;
t3 = pkin(8) * t11 + (t291 - t349) * t36;
t1 = t10 * t248 + t11 * t247;
t2 = [0, 0, 0, 0, 0, qJDD(1), t290, t284, 0, 0, t178, t173, t274, t177, -t364, 0, t254 * t195 - t358, -pkin(1) * t213 - t251 * t195 + t378, t288 + t308, pkin(1) * t195 + pkin(7) * t288, 0, -t274, t364, t178, t173, t177, t251 * (qJ(3) * t218 + t267) + (pkin(2) * t218 + t144 + t241) * t254 + t308, t254 * (-pkin(2) * t216 + t234 - t261) + (-t254 * t278 - t312) * qJ(3) + t358, t251 * t260 - t378 + (pkin(1) + t347) * t213 + (t213 + t278) * t313, pkin(7) * (t144 * t254 + t147 * t251) + (pkin(1) - t279) * (t260 + t369), t296 + t254 * (-t183 * t247 - t248 * t294), t251 * (t203 - t202) + t254 * (t158 * t247 - t160 * t248), t251 * t161 + t254 * (-t186 * t248 - t377), -t296 + t254 * (-t182 * t248 - t188 * t247), t251 * t159 + t254 * (-t185 * t247 - t324), t198 + (t204 * t247 + t206 * t248) * t251 * t305, t251 * (pkin(3) * t129 + t281) + t254 * (pkin(3) * t158 + t330) + pkin(7) * (t129 * t251 + t158 * t254) + t266 * (t172 * t248 - t377), t251 * (pkin(3) * t140 - t65) + t254 * (pkin(3) * t160 - t331) + pkin(7) * (t140 * t251 + t160 * t254) + t266 * (-t187 * t247 - t324), pkin(3) * t332 + t254 * (pkin(3) * t154 - t277) + pkin(7) * (t154 * t254 + t332) + t266 * (t159 * t248 + t161 * t247), t266 * t277 + (pkin(3) + pkin(7)) * (t122 * t254 + t251 * t38), t285, t408, t396, t354, t412, t355, t251 * (-t29 + t398) + t254 * (pkin(3) * t92 - t247 * t44 - t248 * t35) + pkin(7) * (t254 * t92 + t400) + t380, t251 * (-t30 - t414) + t254 * (pkin(3) * t98 - t247 * t45 - t248 * t37) + pkin(7) * (t254 * t98 + t416) - t413, t251 * (pkin(3) * t32 + pkin(4) * t59) + t254 * (-t12 * t247 - t248 * t9 + t394) + pkin(7) * (t251 * t32 + t387) + t266 * (-t247 * t59 + t248 * t63), t251 * (pkin(3) * t5 + pkin(4) * t14) + t254 * (pkin(3) * t70 + pkin(8) * t343 - t13 * t248) + pkin(7) * (t251 * t5 + t254 * t70) + t266 * (t15 * t248 - t343), t285, t396, -t408, t355, -t412, t354, t251 * (t262 + t398) + t254 * (pkin(3) * t94 - t17 * t248 - t19 * t247) + pkin(7) * (t254 * t94 + t400) + t380, t251 * (pkin(3) * t31 + pkin(4) * t57 + t344) + t254 * (-t247 * t8 - t248 * t7 + t394) + pkin(7) * (t251 * t31 + t387) + t266 * (-t247 * t57 + t248 * t61), t251 * (t264 + t414) + t254 * (-pkin(3) * t368 - t16 * t248 - t18 * t247) + pkin(7) * (-t254 * t368 - t416) + t413, t251 * (pkin(3) * t1 + pkin(4) * t10 + t345) + t254 * (pkin(3) * t36 - t247 * t4 - t248 * t3) + pkin(7) * (t1 * t251 + t254 * t36) + t266 * (-t10 * t247 + t11 * t248); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t295, t219, t237, t295, t238, qJDD(2), -t179, -t180, 0, 0, qJDD(2), -t237, -t238, -t295, t219, t295, (-pkin(2) * t251 + qJ(3) * t254) * qJDD(1), -pkin(2) * t220 + qJ(3) * t228 + t147, -pkin(2) * t227 + (qJDD(2) - t221) * qJ(3) + t263, -pkin(2) * t147 + qJ(3) * t144, t183 * t248 - t247 * t294, -t158 * t248 - t160 * t247, -t186 * t247 + t376, -t182 * t247 + t188 * t248, t185 * t248 - t325, (-t204 * t248 + t206 * t247) * t306, qJ(3) * t158 + t129 * t346 + t331, qJ(3) * t160 + t140 * t346 + t330, qJ(3) * t154 + t113 * t346 - t38, qJ(3) * t122 + t346 * t38, t53, -t403, t397, t356, t405, t357, qJ(3) * t92 - t247 * t35 + t248 * t44 + t399, qJ(3) * t98 - t247 * t37 + t248 * t45 + t415, t12 * t248 - t247 * t9 + t32 * t346 + t389, -pkin(8) * t342 + qJ(3) * t70 - t13 * t247 + t346 * t5, t53, t397, t403, t357, -t405, t356, qJ(3) * t94 - t17 * t247 + t19 * t248 + t399, -t247 * t7 + t248 * t8 + t31 * t346 + t389, -qJ(3) * t368 - t16 * t247 + t18 * t248 - t415, qJ(3) * t36 + t1 * t346 - t247 * t3 + t248 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, t220, t227, t147, 0, 0, 0, 0, 0, 0, t129, t140, t113, t38, 0, 0, 0, 0, 0, 0, t381, t46, t32, t5, 0, 0, 0, 0, 0, 0, t381, t31, -t46, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t160, t154, t122, 0, 0, 0, 0, 0, 0, t92, t98, t102, t70, 0, 0, 0, 0, 0, 0, t94, t102, -t368, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, t135, t367, -t322, -t93, t207, -t29, -t30, 0, 0, t322, t367, -t135, t207, t93, -t322, t262, t344, t264, t345; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t362, t367, t366, t25;];
tauJ_reg  = t2;
