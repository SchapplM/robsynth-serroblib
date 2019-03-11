% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:38
% EndTime: 2019-03-09 07:14:56
% DurationCPUTime: 8.16s
% Computational Cost: add. (11040->530), mult. (26357->691), div. (0->0), fcn. (21458->18), ass. (0->275)
t258 = cos(pkin(11));
t267 = cos(qJ(3));
t347 = t267 * t258;
t235 = qJD(1) * t347;
t257 = sin(pkin(11));
t262 = sin(qJ(3));
t357 = t257 * t262;
t323 = qJD(1) * t357;
t196 = t235 - t323;
t399 = qJD(4) + qJD(5);
t426 = t196 - t399;
t209 = t257 * t267 + t258 * t262;
t197 = t209 * qJD(1);
t261 = sin(qJ(4));
t266 = cos(qJ(4));
t332 = t266 * qJD(3);
t163 = t197 * t261 - t332;
t165 = qJD(3) * t261 + t197 * t266;
t260 = sin(qJ(5));
t265 = cos(qJ(5));
t106 = t265 * t163 + t165 * t260;
t259 = sin(qJ(6));
t264 = cos(qJ(6));
t295 = t163 * t260 - t265 * t165;
t296 = t106 * t259 + t264 * t295;
t51 = t264 * t106 - t259 * t295;
t418 = t296 * t51;
t212 = t260 * t261 - t265 * t266;
t344 = t426 * t212;
t353 = t260 * t266;
t213 = t261 * t265 + t353;
t343 = t426 * t213;
t414 = t296 ^ 2 - t51 ^ 2;
t187 = qJD(4) - t196;
t183 = qJD(5) + t187;
t177 = qJD(6) + t183;
t329 = t258 * qJDD(1);
t330 = t257 * qJDD(1);
t326 = qJD(3) * t235 + t262 * t329 + t267 * t330;
t148 = -qJD(3) * t323 + t326;
t338 = qJD(4) * t261;
t98 = qJD(4) * t332 + t261 * qJDD(3) + t266 * t148 - t197 * t338;
t99 = qJD(4) * t165 - t266 * qJDD(3) + t261 * t148;
t274 = qJD(5) * t295 - t260 * t98 - t265 * t99;
t333 = qJD(6) * t264;
t334 = qJD(6) * t259;
t335 = qJD(5) * t265;
t336 = qJD(5) * t260;
t37 = -t163 * t335 - t165 * t336 - t260 * t99 + t265 * t98;
t8 = -t106 * t333 + t259 * t274 + t264 * t37 + t295 * t334;
t412 = t177 * t51 + t8;
t256 = qJ(4) + qJ(5);
t252 = qJ(6) + t256;
t240 = sin(t252);
t241 = cos(t252);
t268 = cos(qJ(1));
t255 = pkin(11) + qJ(3);
t246 = cos(t255);
t263 = sin(qJ(1));
t363 = t246 * t263;
t167 = t240 * t268 - t241 * t363;
t362 = t246 * t268;
t169 = t240 * t263 + t241 * t362;
t237 = -pkin(2) * t258 - pkin(1);
t222 = qJD(1) * t237 + qJD(2);
t117 = -t196 * pkin(3) - t197 * pkin(8) + t222;
t392 = pkin(7) + qJ(2);
t225 = t392 * t257;
t210 = qJD(1) * t225;
t226 = t392 * t258;
t211 = qJD(1) * t226;
t152 = -t262 * t210 + t267 * t211;
t144 = qJD(3) * pkin(8) + t152;
t80 = t266 * t117 - t144 * t261;
t63 = -pkin(9) * t165 + t80;
t44 = pkin(4) * t187 + t63;
t81 = t117 * t261 + t144 * t266;
t64 = -pkin(9) * t163 + t81;
t60 = t265 * t64;
t27 = t260 * t44 + t60;
t420 = pkin(10) * t106;
t21 = t27 - t420;
t19 = t21 * t334;
t245 = sin(t255);
t394 = g(3) * t245;
t402 = -t210 * t267 - t262 * t211;
t143 = -qJD(3) * pkin(3) - t402;
t100 = pkin(4) * t163 + t143;
t56 = pkin(5) * t106 + t100;
t411 = g(1) * t169 - g(2) * t167 + t241 * t394 + t51 * t56 + t19;
t166 = t240 * t363 + t241 * t268;
t168 = -t240 * t362 + t241 * t263;
t199 = t209 * qJD(3);
t298 = t262 * t330 - t267 * t329;
t149 = qJD(1) * t199 + t298;
t142 = qJDD(4) + t149;
t139 = qJDD(5) + t142;
t221 = qJDD(1) * t237 + qJDD(2);
t93 = t149 * pkin(3) - t148 * pkin(8) + t221;
t85 = t266 * t93;
t331 = qJD(1) * qJD(2);
t396 = qJDD(1) * t392 + t331;
t181 = t396 * t257;
t182 = t396 * t258;
t294 = -t262 * t181 + t267 * t182;
t90 = qJDD(3) * pkin(8) + qJD(3) * t402 + t294;
t12 = pkin(4) * t142 - pkin(9) * t98 - qJD(4) * t81 - t261 * t90 + t85;
t337 = qJD(4) * t266;
t286 = t117 * t337 - t144 * t338 + t261 * t93 + t266 * t90;
t17 = -pkin(9) * t99 + t286;
t319 = t265 * t12 - t260 * t17;
t276 = -qJD(5) * t27 + t319;
t2 = pkin(5) * t139 - pkin(10) * t37 + t276;
t313 = -t260 * t12 - t265 * t17 - t44 * t335 + t64 * t336;
t3 = pkin(10) * t274 - t313;
t325 = t264 * t2 - t259 * t3;
t425 = -g(1) * t168 + g(2) * t166 + t240 * t394 + t56 * t296 + t325;
t275 = qJD(6) * t296 - t259 * t37 + t264 * t274;
t406 = -t177 * t296 + t275;
t395 = pkin(8) + pkin(9);
t324 = qJD(4) * t395;
t145 = pkin(3) * t197 - pkin(8) * t196;
t346 = t261 * t145 + t266 * t402;
t372 = t196 * t261;
t424 = -pkin(9) * t372 + t261 * t324 + t346;
t131 = t266 * t145;
t423 = pkin(4) * t197 - t261 * t402 + t131 + (-pkin(9) * t196 + t324) * t266;
t302 = g(1) * t268 + g(2) * t263;
t279 = -g(3) * t246 + t245 * t302;
t339 = qJD(3) * t267;
t340 = qJD(3) * t262;
t292 = -t181 * t267 - t262 * t182 + t210 * t340 - t211 * t339;
t91 = -qJDD(3) * pkin(3) - t292;
t422 = -qJD(4) * pkin(8) * t187 + t279 - t91;
t421 = t338 - t372;
t419 = pkin(10) * t295;
t154 = -t212 * t259 + t213 * t264;
t389 = qJD(6) * t154 + t344 * t259 - t343 * t264;
t416 = t295 * t106;
t321 = t209 * t337;
t208 = -t347 + t357;
t198 = t208 * qJD(3);
t371 = t198 * t261;
t415 = t321 - t371;
t413 = -t106 ^ 2 + t295 ^ 2;
t410 = t106 * t183 + t37;
t58 = t260 * t64;
t26 = t265 * t44 - t58;
t20 = t26 + t419;
t18 = pkin(5) * t183 + t20;
t385 = t264 * t21;
t7 = t259 * t18 + t385;
t409 = -qJD(6) * t7 + t425;
t251 = cos(t256);
t359 = t251 * t263;
t250 = sin(t256);
t360 = t250 * t268;
t172 = -t246 * t359 + t360;
t358 = t251 * t268;
t361 = t250 * t263;
t174 = t246 * t358 + t361;
t408 = g(1) * t174 - g(2) * t172 + t100 * t106 + t251 * t394 + t313;
t171 = t246 * t361 + t358;
t173 = -t246 * t360 + t359;
t407 = -g(1) * t173 + g(2) * t171 + t100 * t295 + t250 * t394 + t276;
t405 = -t183 * t295 + t274;
t404 = t423 * t265;
t133 = t213 * t209;
t305 = pkin(4) * t421 - t152;
t161 = t225 * t267 + t262 * t226;
t227 = t395 * t261;
t228 = t395 * t266;
t342 = -t260 * t227 + t265 * t228;
t401 = t227 * t335 + t228 * t336 + t423 * t260 + t265 * t424;
t400 = qJ(2) * qJDD(1);
t128 = qJDD(6) + t139;
t153 = t264 * t212 + t213 * t259;
t390 = -qJD(6) * t153 + t259 * t343 + t264 * t344;
t398 = -t154 * t128 - t177 * t390;
t397 = -t213 * t139 - t183 * t344;
t391 = t265 * t63 - t58;
t147 = pkin(3) * t208 - pkin(8) * t209 + t237;
t136 = t266 * t147;
t162 = -t225 * t262 + t226 * t267;
t368 = t209 * t266;
t69 = pkin(4) * t208 - pkin(9) * t368 - t162 * t261 + t136;
t155 = t266 * t162;
t345 = t261 * t147 + t155;
t369 = t209 * t261;
t83 = -pkin(9) * t369 + t345;
t387 = t260 * t69 + t265 * t83;
t386 = t197 * t51;
t384 = t296 * t197;
t383 = t98 * t261;
t382 = -pkin(5) * t343 + t305;
t381 = qJDD(1) * pkin(1);
t380 = t106 * t197;
t379 = t295 * t197;
t378 = t128 * t260;
t376 = t163 * t187;
t375 = t163 * t197;
t374 = t165 * t187;
t373 = t165 * t197;
t370 = t198 * t266;
t242 = pkin(4) * t265 + pkin(5);
t364 = t242 * t128;
t354 = t260 * t264;
t352 = t261 * t142;
t351 = t261 * t263;
t350 = t261 * t268;
t349 = t263 * t266;
t127 = t266 * t142;
t348 = t266 * t268;
t341 = t257 ^ 2 + t258 ^ 2;
t243 = -pkin(4) * t266 - pkin(3);
t322 = t209 * t338;
t320 = qJD(6) * t18 + t3;
t118 = -t208 * qJD(2) - qJD(3) * t161;
t146 = pkin(3) * t199 + pkin(8) * t198;
t132 = t266 * t146;
t32 = pkin(9) * t370 + pkin(4) * t199 - t261 * t118 + t132 + (-t155 + (pkin(9) * t209 - t147) * t261) * qJD(4);
t284 = t266 * t118 + t261 * t146 + t147 * t337 - t162 * t338;
t40 = -pkin(9) * t415 + t284;
t317 = -t260 * t40 + t265 * t32;
t316 = -t260 * t63 - t60;
t315 = -t260 * t83 + t265 * t69;
t312 = t341 * qJD(1) ^ 2;
t311 = -qJD(4) * t117 - t90;
t309 = -t265 * t227 - t228 * t260;
t308 = t187 * t266;
t119 = qJD(2) * t209 - t225 * t340 + t226 * t339;
t307 = 0.2e1 * t341;
t306 = -t153 * t128 - t177 * t389;
t125 = -pkin(10) * t213 + t309;
t304 = -pkin(10) * t343 - qJD(6) * t125 + t401;
t126 = -pkin(10) * t212 + t342;
t303 = pkin(5) * t197 + pkin(10) * t344 + t342 * qJD(5) + qJD(6) * t126 - t260 * t424 + t404;
t301 = g(1) * t263 - g(2) * t268;
t300 = -t144 * t337 + t85;
t299 = -t212 * t139 + t183 * t343;
t120 = pkin(4) * t369 + t161;
t134 = t212 * t209;
t87 = t264 * t133 - t134 * t259;
t88 = -t133 * t259 - t134 * t264;
t291 = -t187 * t421 + t127;
t92 = pkin(4) * t415 + t119;
t290 = t260 * t32 + t265 * t40 + t69 * t335 - t336 * t83;
t288 = -t322 - t370;
t285 = -pkin(8) * t142 + t143 * t187;
t282 = -t301 - t381;
t244 = qJDD(2) - t381;
t281 = -t244 - t282;
t41 = pkin(4) * t99 + t91;
t273 = t307 * t331 - t302;
t191 = t246 * t348 + t351;
t190 = -t246 * t350 + t349;
t189 = -t246 * t349 + t350;
t188 = t246 * t351 + t348;
t180 = pkin(5) * t212 + t243;
t89 = pkin(4) * t165 - pkin(5) * t295;
t86 = pkin(5) * t133 + t120;
t47 = -t198 * t353 - t260 * t322 - t336 * t369 + (t368 * t399 - t371) * t265;
t46 = -t133 * t399 + t212 * t198;
t33 = pkin(5) * t47 + t92;
t29 = -pkin(10) * t133 + t387;
t28 = pkin(5) * t208 + pkin(10) * t134 + t315;
t23 = t391 + t419;
t22 = t316 + t420;
t15 = qJD(6) * t88 + t259 * t46 + t264 * t47;
t14 = -qJD(6) * t87 - t259 * t47 + t264 * t46;
t13 = -pkin(5) * t274 + t41;
t6 = t264 * t18 - t21 * t259;
t5 = -pkin(10) * t47 + t290;
t4 = pkin(5) * t199 - pkin(10) * t46 - qJD(5) * t387 + t317;
t1 = [qJDD(1), t301, t302, t281 * t258, -t281 * t257, t307 * t400 + t273 (-t244 + t301) * pkin(1) + (t341 * t400 + t273) * qJ(2), t148 * t209 - t197 * t198, -t148 * t208 - t149 * t209 - t196 * t198 - t197 * t199, -qJD(3) * t198 + qJDD(3) * t209, -qJD(3) * t199 - qJDD(3) * t208, 0, -qJD(3) * t119 - qJDD(3) * t161 + t149 * t237 + t199 * t222 + t208 * t221 + t246 * t301, -qJD(3) * t118 - qJDD(3) * t162 + t148 * t237 - t198 * t222 + t209 * t221 - t245 * t301, t165 * t288 + t368 * t98 -(-t163 * t266 - t165 * t261) * t198 + (-t383 - t266 * t99 + (t163 * t261 - t165 * t266) * qJD(4)) * t209, t127 * t209 + t165 * t199 + t187 * t288 + t98 * t208, -t163 * t199 - t187 * t415 - t99 * t208 - t209 * t352, t142 * t208 + t187 * t199 (-t162 * t337 + t132) * t187 + t136 * t142 + t300 * t208 + t80 * t199 + t119 * t163 + t161 * t99 + t143 * t321 - g(1) * t189 - g(2) * t191 + ((-qJD(4) * t147 - t118) * t187 - t162 * t142 + t311 * t208 + t91 * t209 - t143 * t198) * t261, -g(1) * t188 - g(2) * t190 + t119 * t165 - t142 * t345 + t143 * t288 + t161 * t98 - t187 * t284 - t81 * t199 - t208 * t286 + t368 * t91, -t134 * t37 - t295 * t46, -t106 * t46 - t133 * t37 - t134 * t274 + t295 * t47, -t134 * t139 + t183 * t46 - t199 * t295 + t208 * t37, -t106 * t199 - t133 * t139 - t183 * t47 + t208 * t274, t139 * t208 + t183 * t199, t317 * t183 + t315 * t139 + t319 * t208 + t26 * t199 + t92 * t106 - t120 * t274 + t41 * t133 + t100 * t47 - g(1) * t172 - g(2) * t174 + (-t183 * t387 - t208 * t27) * qJD(5), -g(1) * t171 - g(2) * t173 + t100 * t46 + t120 * t37 - t41 * t134 - t139 * t387 - t183 * t290 - t27 * t199 + t208 * t313 - t295 * t92, -t14 * t296 + t8 * t88, -t14 * t51 + t15 * t296 + t275 * t88 - t8 * t87, t128 * t88 + t14 * t177 - t199 * t296 + t208 * t8, -t128 * t87 - t15 * t177 - t199 * t51 + t208 * t275, t128 * t208 + t177 * t199 (-t259 * t5 + t264 * t4) * t177 + (-t259 * t29 + t264 * t28) * t128 + t325 * t208 + t6 * t199 + t33 * t51 - t86 * t275 + t13 * t87 + t56 * t15 - g(1) * t167 - g(2) * t169 + ((-t259 * t28 - t264 * t29) * t177 - t7 * t208) * qJD(6), -g(1) * t166 - g(2) * t168 + t13 * t88 + t56 * t14 + t19 * t208 - t7 * t199 - t33 * t296 + t86 * t8 + (-(-qJD(6) * t29 + t4) * t177 - t28 * t128 - t2 * t208) * t259 + (-(qJD(6) * t28 + t5) * t177 - t29 * t128 - t320 * t208) * t264; 0, 0, 0, -t329, t330, -t312, -qJ(2) * t312 + qJDD(2) + t282, 0, 0, 0, 0, 0, 0.2e1 * t197 * qJD(3) + t298 (t196 - t323) * qJD(3) + t326, 0, 0, 0, 0, 0, t291 - t375, -t187 ^ 2 * t266 - t352 - t373, 0, 0, 0, 0, 0, t299 - t380, t379 + t397, 0, 0, 0, 0, 0, t306 - t386, t384 + t398; 0, 0, 0, 0, 0, 0, 0, -t197 * t196, -t196 ^ 2 + t197 ^ 2 (-t196 - t323) * qJD(3) + t326, -t298, qJDD(3), qJD(3) * t152 - t197 * t222 + t279 + t292, -t196 * t222 + t302 * t246 - t294 + t394, t165 * t308 + t383 (t98 - t376) * t266 + (-t99 - t374) * t261, t187 * t308 + t352 - t373, t291 + t375, -t187 * t197, -pkin(3) * t99 - t131 * t187 - t152 * t163 - t80 * t197 + (t187 * t402 + t285) * t261 + t422 * t266, -pkin(3) * t98 - t152 * t165 + t346 * t187 + t81 * t197 - t261 * t422 + t285 * t266, t37 * t213 - t295 * t344, -t106 * t344 - t212 * t37 + t213 * t274 - t295 * t343, t379 - t397, t299 + t380, -t183 * t197, t309 * t139 - t243 * t274 + t41 * t212 - t26 * t197 + (-t228 * t335 + (qJD(5) * t227 + t424) * t260 - t404) * t183 + t305 * t106 - t343 * t100 + t279 * t251, t344 * t100 - t342 * t139 + t183 * t401 + t27 * t197 + t41 * t213 + t243 * t37 - t250 * t279 - t295 * t305, t8 * t154 - t296 * t390, -t153 * t8 + t154 * t275 + t296 * t389 - t390 * t51, t384 - t398, t306 + t386, -t177 * t197 (t125 * t264 - t126 * t259) * t128 - t180 * t275 + t13 * t153 - t6 * t197 + t389 * t56 + t382 * t51 + (t259 * t304 - t264 * t303) * t177 + t279 * t241 -(t125 * t259 + t126 * t264) * t128 + t180 * t8 + t13 * t154 + t7 * t197 + t390 * t56 - t382 * t296 + (t259 * t303 + t264 * t304) * t177 - t279 * t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165 * t163, -t163 ^ 2 + t165 ^ 2, t98 + t376, t374 - t99, t142, -g(1) * t190 + g(2) * t188 - t143 * t165 + t187 * t81 + (t311 + t394) * t261 + t300, g(1) * t191 - g(2) * t189 + t143 * t163 + t187 * t80 + t266 * t394 - t286, -t416, t413, t410, t405, t139, -t316 * t183 + (-t106 * t165 + t139 * t265 - t183 * t336) * pkin(4) + t407, t391 * t183 + (-t139 * t260 + t165 * t295 - t183 * t335) * pkin(4) + t408, -t418, t414, t412, t406, t128, t264 * t364 - (t22 * t264 - t23 * t259) * t177 - t89 * t51 + (-t259 * t378 + (-t259 * t265 - t354) * t177 * qJD(5)) * pkin(4) + ((-pkin(4) * t354 - t242 * t259) * t177 - t7) * qJD(6) + t425, t89 * t296 + (-t364 - t2 + (t22 - (-qJD(5) - qJD(6)) * t260 * pkin(4)) * t177) * t259 + (-pkin(4) * t378 + (-pkin(4) * t335 - qJD(6) * t242 + t23) * t177 - t320) * t264 + t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t416, t413, t410, t405, t139, t183 * t27 + t407, t183 * t26 + t408, -t418, t414, t412, t406, t128 -(-t20 * t259 - t385) * t177 + (t128 * t264 - t177 * t334 + t295 * t51) * pkin(5) + t409 (-t177 * t21 - t2) * t259 + (t177 * t20 - t320) * t264 + (-t128 * t259 - t177 * t333 - t295 * t296) * pkin(5) + t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t418, t414, t412, t406, t128, t177 * t7 + t409, t177 * t6 - t259 * t2 - t264 * t320 + t411;];
tau_reg  = t1;
