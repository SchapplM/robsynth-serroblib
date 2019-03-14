% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:23
% EndTime: 2019-03-09 18:23:35
% DurationCPUTime: 7.39s
% Computational Cost: add. (7387->525), mult. (16029->663), div. (0->0), fcn. (11939->14), ass. (0->291)
t250 = sin(qJ(3));
t251 = sin(qJ(2));
t255 = cos(qJ(2));
t392 = cos(qJ(3));
t187 = t250 * t255 + t251 * t392;
t171 = t187 * qJD(1);
t248 = sin(qJ(6));
t249 = sin(qJ(5));
t253 = cos(qJ(6));
t254 = cos(qJ(5));
t184 = t248 * t254 + t249 * t253;
t398 = qJD(5) + qJD(6);
t350 = (t171 + t398) * t184;
t410 = qJD(5) + t171;
t156 = qJD(6) + t410;
t339 = qJD(6) * t248;
t342 = qJD(5) * t249;
t354 = t253 * t254;
t399 = -t248 * t249 + t354;
t411 = t399 * t171 - t248 * t342 - t249 * t339 + t354 * t398;
t242 = qJD(2) + qJD(3);
t356 = t250 * t251;
t298 = t242 * t356;
t327 = t392 * t255;
t308 = qJD(1) * t327;
t317 = qJDD(1) * t392;
t333 = t255 * qJDD(1);
t309 = t242 * t308 + t250 * t333 + t251 * t317;
t98 = qJD(1) * t298 - t309;
t91 = -qJDD(5) + t98;
t90 = -qJDD(6) + t91;
t416 = -t156 * t411 + t184 * t90;
t345 = qJD(1) * t251;
t326 = t250 * t345;
t169 = -t308 + t326;
t138 = -t254 * t169 + t242 * t249;
t415 = t138 * t410;
t185 = -t327 + t356;
t116 = t399 * t185;
t414 = -t156 * t350 - t399 * t90;
t140 = t169 * t249 + t242 * t254;
t293 = t138 * t248 - t253 * t140;
t75 = t253 * t138 + t140 * t248;
t413 = t293 * t75;
t394 = pkin(3) + pkin(9);
t340 = qJD(5) * t394;
t119 = pkin(3) * t171 + qJ(4) * t169;
t164 = t171 * pkin(9);
t92 = t119 + t164;
t412 = t92 + t340;
t409 = t293 ^ 2 - t75 ^ 2;
t338 = qJD(6) * t253;
t241 = qJDD(2) + qJDD(3);
t341 = qJD(5) * t254;
t131 = t242 * t187;
t334 = t251 * qJDD(1);
t297 = t250 * t334 - t255 * t317;
t99 = qJD(1) * t131 + t297;
t58 = t169 * t341 + t254 * t241 - t242 * t342 + t249 * t99;
t316 = t249 * t241 - t254 * t99;
t59 = qJD(5) * t140 + t316;
t16 = -t138 * t338 - t140 * t339 - t248 * t59 + t253 * t58;
t408 = t156 * t75 + t16;
t246 = qJ(5) + qJ(6);
t236 = sin(t246);
t238 = cos(t246);
t252 = sin(qJ(1));
t363 = t238 * t252;
t247 = qJ(2) + qJ(3);
t237 = sin(t247);
t256 = cos(qJ(1));
t364 = t237 * t256;
t147 = t236 * t364 + t363;
t362 = t238 * t256;
t365 = t237 * t252;
t149 = -t236 * t365 + t362;
t239 = cos(t247);
t225 = g(3) * t239;
t393 = pkin(8) + pkin(7);
t198 = t393 * t255;
t190 = qJD(1) * t198;
t173 = t250 * t190;
t197 = t393 * t251;
t188 = qJD(1) * t197;
t382 = qJD(2) * pkin(2);
t177 = -t188 + t382;
t121 = -t392 * t177 + t173;
t337 = qJD(4) + t121;
t389 = t171 * pkin(4);
t352 = t389 + t337;
t70 = -t242 * t394 + t352;
t387 = t255 * pkin(2);
t228 = pkin(1) + t387;
t196 = t228 * qJD(1);
t275 = -t171 * qJ(4) - t196;
t73 = t169 * t394 + t275;
t41 = t249 * t70 + t254 * t73;
t30 = -pkin(10) * t138 + t41;
t27 = t30 * t339;
t235 = t242 * qJ(4);
t176 = t392 * t190;
t122 = t250 * t177 + t176;
t390 = t169 * pkin(4);
t97 = t122 - t390;
t82 = t235 + t97;
t61 = pkin(5) * t138 + t82;
t407 = g(1) * t147 - g(2) * t149 - t236 * t225 + t61 * t75 + t27;
t146 = -t236 * t252 + t237 * t362;
t148 = t236 * t256 + t237 * t363;
t335 = qJD(1) * qJD(2);
t324 = t251 * t335;
t166 = pkin(2) * t324 - t228 * qJDD(1);
t265 = t98 * qJ(4) - t171 * qJD(4) + t166;
t19 = t394 * t99 + t265;
t323 = t255 * t335;
t136 = qJDD(2) * pkin(2) + t393 * (-t323 - t334);
t137 = t393 * (-t324 + t333);
t325 = qJD(3) * t392;
t344 = qJD(3) * t250;
t311 = -t392 * t136 + t250 * t137 + t177 * t344 + t190 * t325;
t291 = qJDD(4) + t311;
t26 = -t98 * pkin(4) - t241 * t394 + t291;
t319 = -t249 * t19 + t254 * t26;
t268 = -qJD(5) * t41 + t319;
t3 = -t91 * pkin(5) - t58 * pkin(10) + t268;
t332 = -t254 * t19 - t249 * t26 - t70 * t341;
t286 = -t73 * t342 - t332;
t4 = -pkin(10) * t59 + t286;
t329 = -t248 * t4 + t253 * t3;
t40 = -t249 * t73 + t254 * t70;
t29 = -pkin(10) * t140 + t40;
t20 = pkin(5) * t410 + t29;
t380 = t253 * t30;
t9 = t20 * t248 + t380;
t406 = -g(1) * t146 - g(2) * t148 - qJD(6) * t9 + t238 * t225 + t61 * t293 + t329;
t267 = qJD(6) * t293 - t248 * t58 - t253 * t59;
t405 = -t156 * t293 + t267;
t310 = -t250 * t136 - t392 * t137 - t177 * t325 + t190 * t344;
t229 = t241 * qJ(4);
t401 = -t242 * qJD(4) - t229;
t49 = t310 + t401;
t28 = -pkin(4) * t99 - t49;
t404 = t28 * t249 + t82 * t341;
t125 = -t250 * t188 + t176;
t103 = t125 - t390;
t331 = pkin(2) * t344;
t403 = (-t103 + t331) * t254;
t117 = t184 * t185;
t402 = -t169 * pkin(5) - pkin(10) * t342;
t126 = -t188 * t392 - t173;
t349 = -pkin(2) * t325 - qJD(4) + t126;
t347 = t239 * pkin(3) + t237 * qJ(4);
t330 = -pkin(5) * t254 - pkin(4);
t400 = pkin(5) * t341 - t330 * t171;
t144 = -t250 * t197 + t198 * t392;
t303 = g(1) * t252 - g(2) * t256;
t328 = qJD(2) * t393;
t189 = t251 * t328;
t191 = t255 * t328;
t79 = t392 * t189 + t250 * t191 + t197 * t325 + t198 * t344;
t397 = t144 * t241 + t237 * t303 - t79 * t242;
t143 = t197 * t392 + t250 * t198;
t80 = t144 * qJD(3) - t250 * t189 + t191 * t392;
t396 = t143 * t241 - t239 * t303 + t80 * t242;
t395 = t171 ^ 2;
t391 = pkin(10) * t171;
t224 = g(3) * t237;
t388 = t241 * pkin(3);
t227 = -pkin(2) * t392 - pkin(3);
t219 = -pkin(9) + t227;
t386 = -pkin(10) + t219;
t385 = -pkin(10) - t394;
t112 = pkin(2) * t345 + t119;
t81 = t112 + t164;
t384 = t249 * t103 + t254 * t81;
t383 = t249 * t97 + t254 * t92;
t84 = t254 * t91;
t379 = t394 * t91;
t378 = t58 * t254;
t376 = t400 - t349;
t375 = t337 + t400;
t290 = -t187 * qJ(4) - t228;
t101 = t185 * t394 + t290;
t113 = t187 * pkin(4) + t143;
t109 = t249 * t113;
t374 = t254 * t101 + t109;
t373 = t122 * t242;
t372 = t131 * t249;
t371 = t156 * t169;
t370 = t410 * t169;
t369 = t171 * t169;
t368 = t171 * t254;
t367 = t185 * t249;
t366 = t185 * t254;
t361 = t239 * t252;
t360 = t239 * t256;
t358 = t249 * t252;
t357 = t249 * t256;
t355 = t252 * t254;
t353 = t254 * t256;
t348 = t389 - t349;
t244 = t251 ^ 2;
t346 = -t255 ^ 2 + t244;
t343 = qJD(5) * t219;
t233 = t251 * t382;
t179 = t386 * t254;
t193 = t385 * t254;
t221 = pkin(2) * t250 + qJ(4);
t322 = -pkin(10) * t185 - t101;
t321 = qJD(6) * t20 + t4;
t320 = -t41 * t169 + t28 * t254;
t315 = t410 * t82;
t314 = t410 ^ 2;
t313 = t410 * t140;
t307 = -t125 + t331;
t305 = -pkin(2) * t251 - pkin(3) * t237;
t304 = g(1) * t256 + g(2) * t252;
t178 = t386 * t249;
t302 = qJD(6) * t178 + t219 * t342 + (-t81 - t391) * t249 + t402 - t403;
t151 = pkin(10) * t368;
t301 = -t179 * t398 - t249 * t331 + t151 + t384;
t192 = t385 * t249;
t87 = t254 * t97;
t300 = qJD(6) * t192 + t402 + t87 + (-t391 - t412) * t249;
t299 = -t193 * t398 + t151 + t383;
t294 = t40 * t169 + t82 * t368 + t404;
t289 = t228 + t347;
t288 = -t249 * t314 - t84;
t287 = -0.2e1 * pkin(1) * t335 - pkin(7) * qJDD(2);
t285 = t185 * t341 + t372;
t284 = -t131 * t254 + t185 * t342;
t130 = -qJD(2) * t327 - t255 * t325 + t298;
t283 = t130 * qJ(4) - t187 * qJD(4) + t233;
t45 = t131 * t394 + t283;
t57 = -t130 * pkin(4) + t80;
t282 = -t101 * t342 + t113 * t341 + t249 * t57 + t254 * t45;
t278 = -g(1) * t364 - g(2) * t365 + t225 + t311;
t277 = g(1) * t360 + g(2) * t361 + t224 + t310;
t276 = t249 * t91 - t254 * t314;
t273 = -t239 * t304 - t224;
t56 = -t131 * pkin(4) - t79;
t259 = qJD(2) ^ 2;
t272 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t259 + t303;
t260 = qJD(1) ^ 2;
t271 = pkin(1) * t260 - pkin(7) * qJDD(1) + t304;
t270 = -t196 * t169 + t277;
t269 = t196 * t171 - t278;
t106 = t169 * pkin(3) + t275;
t266 = t106 * t171 + qJDD(4) + t278;
t264 = -t106 * t169 - t277 - t401;
t64 = t309 + (t169 - t326) * t242;
t15 = pkin(5) * t59 + t28;
t8 = t20 * t253 - t248 * t30;
t262 = t15 * t184 + t8 * t169 + t236 * t273 + t411 * t61;
t261 = t15 * t399 - t9 * t169 + t238 * t273 - t350 * t61;
t240 = t249 * pkin(5);
t222 = qJ(4) + t240;
t200 = qJ(4) * t360;
t199 = qJ(4) * t361;
t195 = t221 + t240;
t162 = -t237 * t358 + t353;
t161 = t237 * t355 + t357;
t160 = t237 * t357 + t355;
t159 = t237 * t353 - t358;
t120 = pkin(3) * t185 + t290;
t115 = -t235 - t122;
t114 = -t185 * pkin(4) + t144;
t111 = -pkin(3) * t242 + t337;
t110 = t254 * t113;
t102 = -t169 ^ 2 + t395;
t83 = t330 * t185 + t144;
t60 = pkin(3) * t131 + t283;
t52 = t254 * t57;
t50 = t291 - t388;
t48 = pkin(10) * t366 + t374;
t43 = t187 * pkin(5) + t322 * t249 + t110;
t38 = t99 * pkin(3) + t265;
t37 = pkin(5) * t284 + t56;
t34 = t117 * t398 - t131 * t399;
t33 = t116 * t398 + t184 * t131;
t32 = -t138 * t169 + t276;
t31 = t140 * t169 + t288;
t21 = -t249 * t313 + t378;
t14 = -t75 * t169 + t416;
t13 = -t169 * t293 + t414;
t10 = (-t59 - t313) * t254 + (-t58 + t415) * t249;
t7 = -pkin(10) * t284 + t282;
t6 = -t130 * pkin(5) + t52 + (-pkin(10) * t131 - t45) * t249 + (t322 * t254 - t109) * qJD(5);
t5 = t16 * t399 + t293 * t350;
t1 = -t16 * t184 + t267 * t399 + t293 * t411 + t350 * t75;
t2 = [qJDD(1), t303, t304, qJDD(1) * t244 + 0.2e1 * t251 * t323, 0.2e1 * t251 * t333 - 0.2e1 * t346 * t335, qJDD(2) * t251 + t255 * t259, qJDD(2) * t255 - t251 * t259, 0, t251 * t287 + t255 * t272, -t251 * t272 + t255 * t287, -t130 * t171 - t187 * t98, t130 * t169 - t131 * t171 + t185 * t98 - t187 * t99, -t130 * t242 + t187 * t241, -t131 * t242 - t185 * t241, 0, -t196 * t131 + t166 * t185 + t169 * t233 - t228 * t99 - t396, t196 * t130 + t166 * t187 + t171 * t233 + t228 * t98 - t397, -t111 * t130 + t115 * t131 - t143 * t98 - t144 * t99 + t169 * t79 + t171 * t80 + t185 * t49 + t187 * t50 - t304, -t106 * t131 - t120 * t99 - t60 * t169 - t38 * t185 + t396, t106 * t130 + t120 * t98 - t60 * t171 - t38 * t187 + t397, t106 * t60 + t111 * t80 + t115 * t79 + t38 * t120 + t50 * t143 - t49 * t144 + (-g(1) * t393 - g(2) * t289) * t256 + (g(1) * t289 - g(2) * t393) * t252, t140 * t285 + t58 * t367 (-t138 * t249 + t140 * t254) * t131 + (-t249 * t59 + t378 + (-t138 * t254 - t140 * t249) * qJD(5)) * t185, -t140 * t130 + t58 * t187 + t285 * t410 - t91 * t367, t138 * t130 - t59 * t187 - t284 * t410 - t91 * t366, -t130 * t410 - t187 * t91 (-t249 * t45 + t52) * t410 - (-t101 * t249 + t110) * t91 + t319 * t187 - t40 * t130 + t56 * t138 + t114 * t59 - g(1) * t162 - g(2) * t160 + (-t131 * t82 - t185 * t28) * t254 + (-t41 * t187 + t82 * t367 - t374 * t410) * qJD(5), g(1) * t161 - g(2) * t159 + t114 * t58 + t41 * t130 + t56 * t140 + t185 * t404 - t286 * t187 - t282 * t410 + t82 * t372 + t374 * t91, t117 * t16 - t293 * t33, t116 * t16 + t117 * t267 + t293 * t34 - t33 * t75, -t117 * t90 + t130 * t293 + t156 * t33 + t16 * t187, -t116 * t90 + t130 * t75 - t156 * t34 + t187 * t267, -t130 * t156 - t187 * t90 (-t248 * t7 + t253 * t6) * t156 - (-t248 * t48 + t253 * t43) * t90 + t329 * t187 - t8 * t130 + t37 * t75 - t83 * t267 - t15 * t116 + t61 * t34 - g(1) * t149 - g(2) * t147 + ((-t248 * t43 - t253 * t48) * t156 - t9 * t187) * qJD(6), g(1) * t148 - g(2) * t146 + t15 * t117 + t9 * t130 + t83 * t16 + t27 * t187 + t61 * t33 - t37 * t293 + (-(-qJD(6) * t48 + t6) * t156 + t43 * t90 - t3 * t187) * t248 + (-(qJD(6) * t43 + t7) * t156 + t48 * t90 - t321 * t187) * t253; 0, 0, 0, -t251 * t260 * t255, t346 * t260, t334, t333, qJDD(2), -g(3) * t255 + t251 * t271, g(3) * t251 + t255 * t271, t369, t102, t64, -t297, t241, t125 * t242 + (-t169 * t345 + t241 * t392 - t242 * t344) * pkin(2) + t269, t126 * t242 + (-t171 * t345 - t241 * t250 - t242 * t325) * pkin(2) + t270, -t221 * t99 - t227 * t98 + (-t115 + t307) * t171 + (t111 + t349) * t169, t112 * t169 + t307 * t242 + (-pkin(3) + t227) * t241 + t266, t112 * t171 + t221 * t241 - t349 * t242 + t264, -t49 * t221 + t50 * t227 - t106 * t112 - g(1) * (t305 * t256 + t200) - g(2) * (t305 * t252 + t199) - g(3) * (t347 + t387) + t349 * t115 + t307 * t111, t21, t10, t31, t32, t370, -t219 * t84 + t221 * t59 + t403 * t410 + t348 * t138 + ((t81 - t343) * t410 + t273) * t249 + t294, t221 * t58 + t384 * t410 + t348 * t140 + (-t343 * t410 + t273) * t254 + (t219 * t91 - t331 * t410 - t315) * t249 + t320, t5, t1, t13, t14, t371 -(-t178 * t248 + t179 * t253) * t90 - t195 * t267 + t376 * t75 + (t248 * t301 - t253 * t302) * t156 + t262 (t178 * t253 + t179 * t248) * t90 + t195 * t16 - t376 * t293 + (t248 * t302 + t253 * t301) * t156 + t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t369, t102, t64, -t297, t241, t269 + t373, -t121 * t242 + t270, pkin(3) * t98 - qJ(4) * t99 + (-t115 - t122) * t171 + (t111 - t337) * t169, t119 * t169 + t266 - t373 - 0.2e1 * t388, t119 * t171 + t337 * t242 + t229 + t264, -t49 * qJ(4) - t50 * pkin(3) - t106 * t119 - t111 * t122 - g(1) * (-pkin(3) * t364 + t200) - g(2) * (-pkin(3) * t365 + t199) - g(3) * t347 - t337 * t115, t21, t10, t31, t32, t370, t254 * t379 + qJ(4) * t59 - t87 * t410 + t352 * t138 + (t410 * t412 + t273) * t249 + t294, qJ(4) * t58 + t383 * t410 + t352 * t140 + (-t315 - t379) * t249 + (t340 * t410 + t273) * t254 + t320, t5, t1, t13, t14, t371 -(-t192 * t248 + t193 * t253) * t90 - t222 * t267 + t375 * t75 + (t248 * t299 - t253 * t300) * t156 + t262 (t192 * t253 + t193 * t248) * t90 + t222 * t16 - t375 * t293 + (t248 * t300 + t253 * t299) * t156 + t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t241 - t369, -t242 ^ 2 - t395, t115 * t242 + t266 - t388, 0, 0, 0, 0, 0, -t242 * t138 + t288, -t242 * t140 + t276, 0, 0, 0, 0, 0, -t242 * t75 + t414, t242 * t293 + t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140 * t138, -t138 ^ 2 + t140 ^ 2, t58 + t415, -t316 + (-qJD(5) + t410) * t140, -t91, -g(1) * t159 - g(2) * t161 - t82 * t140 + t254 * t225 + t41 * t410 + t268, g(1) * t160 - g(2) * t162 + t82 * t138 + t40 * t410 + (qJD(5) * t73 - t225) * t249 + t332, -t413, t409, t408, t405, -t90 -(-t248 * t29 - t380) * t156 + (-t140 * t75 - t156 * t339 - t253 * t90) * pkin(5) + t406 (-t156 * t30 - t3) * t248 + (t156 * t29 - t321) * t253 + (t140 * t293 - t156 * t338 + t248 * t90) * pkin(5) + t407; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t413, t409, t408, t405, -t90, t9 * t156 + t406, t8 * t156 - t248 * t3 - t253 * t321 + t407;];
tau_reg  = t2;