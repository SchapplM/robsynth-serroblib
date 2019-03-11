% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPPR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:09:34
% EndTime: 2019-03-09 16:09:50
% DurationCPUTime: 6.46s
% Computational Cost: add. (6792->633), mult. (16705->794), div. (0->0), fcn. (13003->10), ass. (0->296)
t242 = sin(qJ(3));
t449 = -t242 * qJ(4) - pkin(2);
t239 = sin(pkin(6));
t246 = cos(qJ(2));
t361 = qJDD(1) * t246;
t213 = t239 * t361;
t243 = sin(qJ(2));
t363 = qJD(1) * qJD(2);
t338 = t243 * t363;
t448 = t239 * t338 - t213;
t150 = qJDD(3) + t448;
t426 = pkin(3) + pkin(4);
t358 = t426 * t150;
t362 = qJDD(1) * t243;
t337 = t239 * t362;
t339 = t246 * t363;
t447 = t239 * t339 + t337;
t373 = qJD(1) * t239;
t345 = t246 * t373;
t319 = t242 * t345;
t370 = qJD(3) * t242;
t446 = t319 - t370;
t191 = -qJD(3) + t345;
t396 = cos(pkin(6));
t332 = t396 * qJD(1);
t317 = pkin(1) * t332;
t155 = pkin(8) * t345 + t243 * t317;
t296 = t332 + qJD(2);
t108 = pkin(9) * t296 + t155;
t294 = -pkin(2) * t246 - pkin(9) * t243 - pkin(1);
t114 = t294 * t373;
t245 = cos(qJ(3));
t397 = t242 * t108 - t245 * t114;
t445 = qJD(4) + t397;
t423 = cos(qJ(1));
t311 = t396 * t423;
t422 = sin(qJ(1));
t169 = t422 * t243 - t246 * t311;
t241 = sin(qJ(6));
t244 = cos(qJ(6));
t170 = t243 * t311 + t422 * t246;
t348 = t239 * t423;
t97 = t170 * t242 + t245 * t348;
t444 = t169 * t244 + t241 * t97;
t443 = -t169 * t241 + t244 * t97;
t280 = qJD(3) * t296;
t323 = t396 * qJDD(1);
t292 = t323 + qJDD(2);
t371 = qJD(2) * t246;
t342 = t242 * t371;
t368 = qJD(3) * t245;
t61 = t239 * (qJD(1) * (t243 * t368 + t342) + t242 * t362) + t242 * t280 - t245 * t292;
t235 = t239 ^ 2;
t442 = 0.2e1 * t235;
t418 = pkin(9) - qJ(5);
t346 = t243 * t373;
t132 = t242 * t346 - t245 * t296;
t441 = -t61 * qJ(5) - t132 * qJD(5);
t152 = -pkin(8) * t346 + t246 * t317;
t308 = pkin(2) * t243 - pkin(9) * t246;
t153 = t308 * t373;
t378 = t245 * t152 + t242 * t153;
t65 = qJ(4) * t346 + t378;
t440 = -qJD(5) * t245 - t418 * t370 - t65;
t394 = t132 * t191;
t386 = t239 * t243;
t354 = t242 * t386;
t318 = qJD(3) * t354;
t60 = qJD(1) * t318 - t242 * t292 + (-t280 - t447) * t245;
t439 = -t60 + t394;
t134 = t242 * t296 + t245 * t346;
t392 = t134 * t191;
t438 = t61 - t392;
t437 = -t134 * qJ(5) + t445;
t131 = t134 ^ 2;
t436 = -t191 ^ 2 - t131;
t349 = pkin(1) * t396;
t385 = t239 * t246;
t271 = pkin(8) * t385 + t243 * t349;
t157 = t271 * qJD(2);
t168 = t242 * t396 + t245 * t386;
t93 = qJD(3) * t168 + t239 * t342;
t343 = t239 * t371;
t94 = -t318 + (qJD(3) * t396 + t343) * t245;
t32 = t93 * pkin(3) - t94 * qJ(4) - t168 * qJD(4) + t157;
t167 = -t245 * t396 + t354;
t329 = -t167 * pkin(3) + t168 * qJ(4);
t136 = t242 * t152;
t326 = -t153 * t245 + t136;
t435 = -qJD(5) * t242 + t418 * t368 - t326;
t434 = (qJDD(2) + 0.2e1 * t323) * t239;
t433 = t242 * qJD(4) + t155 + (-t245 * t345 + t368) * qJ(4);
t432 = t61 * pkin(3) + t60 * qJ(4) - t134 * qJD(4);
t419 = pkin(9) * t150;
t107 = -t296 * pkin(2) - t152;
t44 = t132 * pkin(3) - t134 * qJ(4) + t107;
t431 = t191 * t44 + t419;
t81 = t132 * t244 + t191 * t241;
t23 = qJD(6) * t81 + t244 * t150 + t241 * t61;
t123 = qJD(6) + t134;
t56 = -qJDD(6) + t60;
t403 = t244 * t56;
t274 = t123 ^ 2 * t241 + t403;
t310 = t396 * t422;
t172 = -t243 * t310 + t423 * t246;
t347 = t239 * t422;
t102 = t172 * t245 + t242 * t347;
t98 = t170 * t245 - t242 * t348;
t282 = g(1) * t102 + g(2) * t98 + g(3) * t168;
t425 = -pkin(4) - pkin(10);
t364 = pkin(3) - t425;
t142 = t150 * qJ(4);
t177 = qJD(4) * t191;
t293 = qJD(2) * t317;
t314 = pkin(1) * t323;
t270 = t448 * pkin(8) - t243 * t314 - t246 * t293;
t74 = pkin(9) * t292 - t270;
t286 = t308 * qJD(2);
t77 = (qJD(1) * t286 + qJDD(1) * t294) * t239;
t331 = t108 * t370 - t114 * t368 - t242 * t77 - t245 * t74;
t12 = t142 - t177 - t331;
t8 = -t12 + t441;
t6 = pkin(5) * t150 - t8;
t76 = t134 * pkin(3) + t132 * qJ(4);
t430 = t123 * (-pkin(5) * t132 - qJD(6) * t364 + t425 * t134 - t76) + t282 - t6;
t427 = 0.2e1 * t142;
t248 = qJD(1) ^ 2;
t424 = pkin(4) * t61;
t421 = pkin(3) * t150;
t420 = pkin(4) * t167;
t240 = qJ(4) + pkin(5);
t417 = pkin(9) * qJD(3);
t416 = qJ(4) * t61;
t415 = qJ(5) * t60;
t79 = t132 * t241 - t244 * t191;
t414 = t123 * t79;
t413 = t123 * t81;
t412 = t132 * t79;
t411 = t132 * t81;
t179 = t191 * qJ(4);
t59 = t245 * t108 + t242 * t114;
t37 = qJ(5) * t132 + t59;
t33 = t179 - t37;
t410 = t191 * t33;
t51 = -t179 + t59;
t409 = t191 * t51;
t408 = t191 * t59;
t407 = t191 * t79;
t406 = t191 * t81;
t365 = qJD(6) * t244;
t366 = qJD(6) * t241;
t22 = -t132 * t366 - t241 * t150 + t191 * t365 + t244 * t61;
t405 = t22 * t241;
t404 = t241 * t56;
t402 = t446 * t426 + t433;
t384 = t242 * t246;
t401 = -(pkin(5) * t243 + qJ(5) * t384) * t373 + t440;
t400 = qJ(5) * t319 - t440;
t399 = -t446 * pkin(3) - t433;
t382 = t245 * t246;
t356 = qJ(5) * t382;
t398 = -(-t426 * t243 - t356) * t373 + t435;
t395 = qJ(5) * t167;
t393 = t134 * t132;
t389 = t169 * t245;
t171 = t423 * t243 + t246 * t310;
t388 = t171 * t245;
t387 = t235 * t248;
t383 = t244 * t246;
t148 = pkin(9) * t396 + t271;
t376 = pkin(2) * t385 + pkin(9) * t386;
t149 = -pkin(1) * t239 - t376;
t379 = t245 * t148 + t242 * t149;
t377 = -pkin(8) * t386 + t246 * t349;
t237 = t243 ^ 2;
t374 = -t246 ^ 2 + t237;
t372 = qJD(2) * t243;
t369 = qJD(3) * t244;
t367 = qJD(4) * t246;
t218 = pkin(3) * t385;
t185 = -t245 * pkin(3) + t449;
t357 = qJ(4) * t385;
t355 = t246 * t387;
t352 = -pkin(3) * t389 + t449 * t169;
t351 = -pkin(3) * t388 + t449 * t171;
t147 = -t396 * pkin(2) - t377;
t344 = t239 * t372;
t341 = t239 * t367;
t340 = pkin(1) * t442;
t335 = -t97 * pkin(3) + qJ(4) * t98;
t101 = t172 * t242 - t245 * t347;
t334 = -t101 * pkin(3) + qJ(4) * t102;
t330 = t108 * t368 + t114 * t370 + t242 * t74 - t245 * t77;
t327 = -t242 * t148 + t149 * t245;
t325 = t123 * t244;
t324 = t364 * t243;
t173 = t245 * pkin(4) - t185;
t321 = t447 * pkin(8) + t243 * t293 - t246 * t314;
t320 = t245 * t218 + t242 * t357 + t376;
t313 = g(1) * t97 - g(2) * t101;
t194 = t418 * t242;
t312 = qJD(6) * t194 - t433 + t191 * (pkin(5) * t245 - t242 * t364);
t309 = t239 * t248 * t396;
t307 = -g(1) * t98 + g(2) * t102;
t306 = g(1) * t171 + g(2) * t169;
t305 = g(1) * t169 - g(2) * t171;
t304 = g(1) * t172 + g(2) * t170;
t196 = t292 * pkin(2);
t75 = -t196 + t321;
t14 = t75 + t432;
t266 = qJDD(5) - t14;
t4 = -pkin(5) * t60 + t425 * t61 + t266;
t300 = qJDD(4) + t330;
t262 = -qJD(5) * t134 + t300 + t415;
t5 = -t150 * t364 + t262;
t303 = t241 * t4 + t244 * t5;
t122 = pkin(5) * t242 + pkin(10) * t245 + t173;
t302 = -qJD(6) * t122 + (-t324 - t356) * t373 - t435;
t64 = t218 - t327;
t276 = (-t241 * t243 + t242 * t383) * t239;
t113 = qJD(1) * t276;
t301 = t242 * t369 - t113;
t289 = qJD(5) - t44;
t20 = pkin(5) * t134 + t425 * t132 + t289;
t24 = t191 * t364 + t437;
t9 = t20 * t244 - t24 * t241;
t10 = t20 * t241 + t24 * t244;
t62 = t147 - t329;
t27 = pkin(5) * t168 + t425 * t167 - t62;
t40 = pkin(4) * t385 - qJ(5) * t168 + t64;
t35 = pkin(10) * t385 + t40;
t299 = -t241 * t35 + t244 * t27;
t298 = t241 * t27 + t244 * t35;
t295 = 0.2e1 * t332 + qJD(2);
t291 = g(3) * t386 + t304;
t63 = -t357 + t379;
t154 = t239 * t286;
t156 = t377 * qJD(2);
t290 = -t148 * t368 - t149 * t370 + t154 * t245 - t242 * t156;
t288 = -t167 * t241 + t239 * t383;
t96 = t167 * t244 + t241 * t385;
t287 = t241 * t384 + t243 * t244;
t285 = t191 * t245;
t284 = t423 * pkin(1) + t172 * pkin(2) + t102 * pkin(3) + pkin(8) * t347 + qJ(4) * t101;
t281 = -t148 * t370 + t149 * t368 + t242 * t154 + t245 * t156;
t279 = t194 * t56 + t304;
t275 = -t123 * t325 + t404;
t273 = qJ(4) * t344 + t281;
t272 = g(3) * t385 - t306;
t269 = -t422 * pkin(1) - t170 * pkin(2) - pkin(3) * t98 + pkin(8) * t348 - qJ(4) * t97;
t11 = t266 - t424;
t268 = -t11 + t272;
t29 = -pkin(5) * t191 - t33;
t267 = -t364 * t56 + (-t29 + t37) * t123;
t265 = -t107 * t191 - t419;
t263 = t150 - t393;
t2 = -qJD(6) * t10 - t241 * t5 + t244 * t4;
t261 = g(1) * t101 + g(2) * t97 + g(3) * t167 - t330;
t260 = -t282 - t331;
t259 = -qJ(5) * t94 - qJD(5) * t168 - t290;
t258 = t191 * t417 - t272;
t257 = qJDD(4) - t261;
t256 = -t272 - t321;
t255 = -t14 + t258;
t254 = t60 + t394;
t253 = t191 * t397 - t260;
t252 = qJ(5) * t93 + qJD(5) * t167 + t273;
t250 = t134 * t44 + t257;
t31 = -pkin(4) * t132 + t289;
t249 = t415 + (-qJD(5) - t31) * t134 + t257;
t195 = t418 * t245;
t130 = t132 ^ 2;
t112 = t287 * t373;
t69 = t101 * t244 - t171 * t241;
t68 = -t101 * t241 - t171 * t244;
t67 = -pkin(3) * t346 + t326;
t50 = -pkin(4) * t134 - t76;
t49 = pkin(3) * t191 + t445;
t47 = -t63 - t395;
t45 = -t62 - t420;
t42 = qJD(6) * t288 - t241 * t344 + t244 * t93;
t41 = qJD(6) * t96 + t241 * t93 + t244 * t344;
t39 = -t240 * t385 + t379 + t395;
t30 = -pkin(3) * t344 - t290;
t28 = t426 * t191 + t437;
t26 = t273 - t341;
t21 = -t93 * pkin(4) - t32;
t19 = -t252 + t341;
t18 = -t426 * t344 + t259;
t17 = (pkin(5) * t372 - t367) * t239 + t252;
t16 = -qJD(2) * t239 * t324 + t259;
t15 = t300 - t421;
t13 = t94 * pkin(5) + t425 * t93 - t32;
t7 = -t358 + t262;
t1 = t9 * qJD(6) + t303;
t3 = [qJDD(1), g(1) * t422 - g(2) * t423, g(1) * t423 + g(2) * t422 (qJDD(1) * t237 + 0.2e1 * t246 * t338) * t235 (t243 * t361 - t363 * t374) * t442, t434 * t243 + t295 * t343, t434 * t246 - t295 * t344, t292 * t396, -t157 * t296 + t377 * t292 - t321 * t396 + g(1) * t170 - g(2) * t172 + (-t338 + t361) * t340, -t156 * t296 - t271 * t292 + t270 * t396 + (-t339 - t362) * t340 - t305, t134 * t94 - t168 * t60, -t132 * t94 - t134 * t93 + t167 * t60 - t168 * t61, t150 * t168 - t191 * t94 + (t134 * t372 + t246 * t60) * t239, -t150 * t167 + t191 * t93 + (-t132 * t372 + t246 * t61) * t239 (-t150 * t246 - t191 * t372) * t239, -t290 * t191 + t327 * t150 + t157 * t132 + t147 * t61 + t75 * t167 + t107 * t93 + (t246 * t330 - t372 * t397) * t239 - t307, t281 * t191 - t379 * t150 + t157 * t134 - t147 * t60 + t75 * t168 + t107 * t94 + (-t246 * t331 - t372 * t59) * t239 - t313, t132 * t32 + t14 * t167 - t150 * t64 + t191 * t30 + t44 * t93 + t61 * t62 + (t15 * t246 - t372 * t49) * t239 - t307, -t12 * t167 - t132 * t26 + t134 * t30 + t15 * t168 + t49 * t94 - t51 * t93 - t60 * t64 - t61 * t63 + t305, -t134 * t32 - t14 * t168 + t150 * t63 - t191 * t26 - t44 * t94 + t60 * t62 + (-t12 * t246 + t372 * t51) * t239 + t313, t12 * t63 + t51 * t26 + t14 * t62 + t44 * t32 + t15 * t64 + t49 * t30 - g(1) * (-t169 * pkin(9) + t269) - g(2) * (pkin(9) * t171 + t284) t11 * t168 + t134 * t21 - t150 * t47 + t19 * t191 + t31 * t94 - t45 * t60 + (t246 * t8 - t33 * t372) * t239 + t313, t11 * t167 + t132 * t21 + t150 * t40 - t18 * t191 + t31 * t93 + t45 * t61 + (-t246 * t7 + t28 * t372) * t239 + t307, -t132 * t19 - t134 * t18 - t167 * t8 - t168 * t7 - t28 * t94 - t33 * t93 + t40 * t60 - t47 * t61 - t305, t7 * t40 + t28 * t18 + t8 * t47 + t33 * t19 + t11 * t45 + t31 * t21 - g(1) * (-pkin(4) * t98 - t418 * t169 + t269) - g(2) * (pkin(4) * t102 + t418 * t171 + t284) t22 * t96 + t42 * t81, t22 * t288 - t23 * t96 - t41 * t81 - t42 * t79, t123 * t42 + t168 * t22 - t56 * t96 + t81 * t94, -t123 * t41 - t168 * t23 - t288 * t56 - t79 * t94, t123 * t94 - t168 * t56 (-qJD(6) * t298 + t13 * t244 - t16 * t241) * t123 - t299 * t56 + t2 * t168 + t9 * t94 + t17 * t79 + t39 * t23 - t6 * t288 + t29 * t41 + g(1) * t443 - g(2) * t69 -(qJD(6) * t299 + t13 * t241 + t16 * t244) * t123 + t298 * t56 - t1 * t168 - t10 * t94 + t17 * t81 + t39 * t22 + t6 * t96 + t29 * t42 - g(1) * t444 - g(2) * t68; 0, 0, 0, -t243 * t355, t374 * t387, -t246 * t309 + t337, t243 * t309 + t213, t292, pkin(1) * t243 * t387 + t155 * t296 + t256, pkin(1) * t355 + t152 * t296 + t270 + t291, -t134 * t285 - t242 * t60, -t438 * t242 + t439 * t245, -t191 * t368 + t242 * t150 + (-t134 * t243 + t191 * t382) * t373, t191 * t370 + t245 * t150 + (t132 * t243 - t191 * t384) * t373, t191 * t346, t397 * t346 - pkin(2) * t61 - t155 * t132 - t136 * t191 + t265 * t242 + (-t75 + (t153 + t417) * t191 - t272) * t245, pkin(2) * t60 - t378 * t191 + t59 * t346 - t155 * t134 + t265 * t245 + (-t258 + t75) * t242, t399 * t132 + t185 * t61 - t191 * t67 - t242 * t431 + t255 * t245 + t49 * t346, t132 * t65 - t134 * t67 + (t12 - t191 * t49 + (qJD(3) * t134 - t61) * pkin(9)) * t245 + (t15 + t409 + (qJD(3) * t132 - t60) * pkin(9)) * t242 - t291, -t399 * t134 + t185 * t60 + t191 * t65 + t255 * t242 + t245 * t431 - t51 * t346, t14 * t185 - t51 * t65 - t49 * t67 - g(1) * t351 - g(2) * t352 - g(3) * t320 + t399 * t44 + (t12 * t245 + t15 * t242 + (-t242 * t51 + t245 * t49) * qJD(3) - t304) * pkin(9), t31 * t368 + t150 * t195 - t173 * t60 + t400 * t191 + t402 * t134 + (t243 * t33 - t31 * t382) * t373 - t268 * t242, t31 * t370 + t150 * t194 + t173 * t61 - t398 * t191 + t402 * t132 + (-t243 * t28 - t31 * t384) * t373 + t268 * t245, t194 * t60 + t195 * t61 - t398 * t134 - t400 * t132 + (t191 * t28 + t8) * t245 + (-t7 + t410) * t242 + t291, t7 * t194 - t8 * t195 + t11 * t173 - g(1) * (-pkin(4) * t388 + t418 * t172 + t351) - g(2) * (-pkin(4) * t389 + t418 * t170 + t352) - g(3) * ((pkin(4) * t382 - qJ(5) * t243) * t239 + t320) + t400 * t33 + t402 * t31 + t398 * t28, -t22 * t244 * t245 + (t245 * t366 + t301) * t81, t112 * t81 + t113 * t79 + (-t241 * t81 - t244 * t79) * t370 + (t405 + t23 * t244 + (-t241 * t79 + t244 * t81) * qJD(6)) * t245, t22 * t242 + t301 * t123 + (t123 * t366 + t403 - t406) * t245, -t23 * t242 + (-t241 * t370 + t112) * t123 + (t123 * t365 - t404 + t407) * t245, -t123 * t285 - t242 * t56, -t122 * t403 - t29 * t112 + t195 * t23 + t401 * t79 + t279 * t241 + (t29 * t241 * qJD(3) + t244 * t306 + t2) * t242 - g(3) * t276 + (t241 * t302 - t244 * t312) * t123 + (-t191 * t9 - t6 * t241 - t29 * t365) * t245, t122 * t404 - t29 * t113 + t195 * t22 + t401 * t81 + t279 * t244 + (-t241 * t306 + t29 * t369 - t1) * t242 + g(3) * t287 * t239 + (t241 * t312 + t244 * t302) * t123 + (t10 * t191 - t6 * t244 + t29 * t366) * t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t393, -t130 + t131, -t254, -t392 - t61, t150, -t107 * t134 + t261 - t408, t107 * t132 + t253, -t132 * t76 - t250 - t408 + 0.2e1 * t421, pkin(3) * t60 - t416 + (t51 - t59) * t134 + (t49 - t445) * t132, -t132 * t44 + t134 * t76 - 0.2e1 * t177 - t253 + t427, -t15 * pkin(3) - g(1) * t334 - g(2) * t335 - g(3) * t329 + t12 * qJ(4) - t44 * t76 + t445 * t51 - t49 * t59, t132 * t31 - t134 * t50 - t191 * t437 - t177 + t260 + t427 - t441, -t132 * t50 + t191 * t37 + t249 - 0.2e1 * t358, t416 - t426 * t60 + (t33 + t37) * t134 + (-t28 + t437) * t132, -t7 * t426 - t8 * qJ(4) - t28 * t37 - t31 * t50 - g(1) * (-pkin(4) * t101 + t334) - g(2) * (-pkin(4) * t97 + t335) - g(3) * (t329 - t420) - t437 * t33, -t325 * t81 - t405 (-t22 + t414) * t244 + (t23 + t413) * t241, t275 + t411, t274 - t412, t123 * t132, t9 * t132 + t240 * t23 + t267 * t241 - t244 * t430 + t437 * t79, -t10 * t132 + t240 * t22 + t241 * t430 + t267 * t244 + t437 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t263, -t254, t436, t250 + t409 - t421, t436, t263, t254, t249 - t358 - t410, 0, 0, 0, 0, 0, t275 + t407, t274 + t406; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t439, t438, -t131 - t130, t132 * t33 + t134 * t28 + qJDD(5) + t196 + t256 - t424 - t432, 0, 0, 0, 0, 0, -t274 - t412, t275 - t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t79, -t79 ^ 2 + t81 ^ 2, t22 + t414, -t23 + t413, -t56, -g(1) * t68 + g(2) * t444 - g(3) * t288 + t10 * t123 - t29 * t81 + t2, t29 * t79 + g(1) * t69 + g(2) * t443 + g(3) * t96 - t303 + (t123 - qJD(6)) * t9;];
tau_reg  = t3;
