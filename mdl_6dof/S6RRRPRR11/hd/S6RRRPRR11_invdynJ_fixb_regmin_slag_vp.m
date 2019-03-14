% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:34:19
% EndTime: 2019-03-09 19:34:46
% DurationCPUTime: 11.50s
% Computational Cost: add. (10964->659), mult. (26494->883), div. (0->0), fcn. (21495->12), ass. (0->299)
t284 = cos(qJ(3));
t423 = cos(pkin(6));
t362 = t423 * qJD(1);
t322 = t362 + qJD(2);
t280 = sin(qJ(3));
t281 = sin(qJ(2));
t277 = sin(pkin(6));
t399 = qJD(1) * t277;
t376 = t281 * t399;
t352 = t280 * t376;
t177 = -t284 * t322 + t352;
t179 = t280 * t322 + t284 * t376;
t279 = sin(qJ(5));
t283 = cos(qJ(5));
t117 = t177 * t279 + t179 * t283;
t285 = cos(qJ(2));
t387 = qJDD(1) * t285;
t257 = t277 * t387;
t388 = qJD(1) * qJD(2);
t370 = t281 * t388;
t348 = t277 * t370;
t192 = qJDD(3) - t257 + t348;
t188 = -qJDD(5) + t192;
t375 = t285 * t399;
t236 = -qJD(3) + t375;
t227 = qJD(5) + t236;
t278 = sin(qJ(6));
t282 = cos(qJ(6));
t306 = qJD(3) * t322;
t355 = t423 * qJDD(1);
t317 = t355 + qJDD(2);
t385 = t281 * qJDD(1);
t368 = t277 * t385;
t369 = t285 * t388;
t473 = t277 * t369 + t368;
t100 = qJD(3) * t352 - t280 * t317 + (-t306 - t473) * t284;
t396 = qJD(2) * t285;
t372 = t280 * t396;
t394 = qJD(3) * t284;
t101 = t277 * (qJD(1) * (t281 * t394 + t372) + t280 * t385) + t280 * t306 - t284 * t317;
t392 = qJD(5) * t283;
t393 = qJD(5) * t279;
t33 = -t283 * t100 + t279 * t101 + t177 * t392 - t179 * t393;
t390 = qJD(6) * t282;
t391 = qJD(6) * t278;
t12 = -t117 * t391 - t278 * t188 + t227 * t390 + t282 * t33;
t81 = t117 * t282 + t227 * t278;
t13 = qJD(6) * t81 + t282 * t188 + t278 * t33;
t459 = -t283 * t177 + t179 * t279;
t470 = qJD(6) + t459;
t431 = t470 * t81;
t79 = t117 * t278 - t282 * t227;
t432 = t470 * t79;
t494 = (t13 + t431) * t278 - (t12 - t432) * t282;
t413 = t277 * t281;
t206 = t280 * t413 - t284 * t423;
t207 = t280 * t423 + t284 * t413;
t326 = t206 * t283 - t207 * t279;
t442 = cos(qJ(1));
t343 = t423 * t442;
t441 = sin(qJ(1));
t210 = t281 * t343 + t285 * t441;
t378 = t277 * t442;
t143 = t210 * t280 + t284 * t378;
t144 = t210 * t284 - t280 * t378;
t463 = t143 * t283 - t144 * t279;
t342 = t423 * t441;
t212 = -t281 * t342 + t285 * t442;
t377 = t277 * t441;
t147 = t212 * t280 - t284 * t377;
t148 = t212 * t284 + t280 * t377;
t88 = t147 * t283 - t148 * t279;
t309 = g(1) * t88 + g(2) * t463 + g(3) * t326;
t349 = pkin(1) * t362;
t318 = qJD(2) * t349;
t346 = pkin(1) * t355;
t380 = pkin(8) * t257 + t281 * t346 + t285 * t318;
t297 = -pkin(8) * t348 + t380;
t120 = pkin(9) * t317 + t297;
t340 = pkin(2) * t281 - pkin(9) * t285;
t311 = t340 * qJD(2);
t319 = -pkin(2) * t285 - pkin(9) * t281 - pkin(1);
t123 = (qJD(1) * t311 + qJDD(1) * t319) * t277;
t200 = pkin(8) * t375 + t281 * t349;
t154 = pkin(9) * t322 + t200;
t191 = t319 * t277;
t166 = qJD(1) * t191;
t395 = qJD(3) * t280;
t354 = t280 * t120 - t284 * t123 + t154 * t394 + t166 * t395;
t320 = qJDD(4) + t354;
t446 = pkin(3) + pkin(4);
t16 = pkin(10) * t100 - t192 * t446 + t320;
t187 = t192 * qJ(4);
t223 = t236 * qJD(4);
t308 = t284 * t120 + t280 * t123 - t154 * t395 + t166 * t394;
t26 = t187 - t223 + t308;
t19 = pkin(10) * t101 + t26;
t98 = -t280 * t154 + t284 * t166;
t410 = qJD(4) - t98;
t474 = -pkin(10) * t179 + t410;
t54 = t236 * t446 + t474;
t225 = t236 * qJ(4);
t99 = t284 * t154 + t280 * t166;
t68 = pkin(10) * t177 + t99;
t60 = -t225 + t68;
t365 = -t283 * t16 + t279 * t19 + t60 * t392 + t54 * t393;
t4 = pkin(5) * t188 + t365;
t493 = t309 + t4;
t490 = t284 * t375 - t394;
t351 = t280 * t375;
t481 = t351 - t395;
t197 = -pkin(8) * t376 + t285 * t349;
t153 = -t322 * pkin(2) - t197;
t72 = t177 * pkin(3) - t179 * qJ(4) + t153;
t59 = -pkin(4) * t177 - t72;
t489 = -t117 * t59 - t309 - t365;
t209 = t281 * t441 - t285 * t343;
t85 = t143 * t279 + t144 * t283;
t488 = t209 * t282 + t278 * t85;
t487 = -t209 * t278 + t282 * t85;
t430 = t12 * t278;
t482 = t470 * t282;
t486 = t482 * t81 + t430;
t364 = -t100 * t279 - t283 * t101;
t34 = qJD(5) * t117 + t364;
t30 = qJDD(6) + t34;
t427 = t278 * t30;
t485 = -t117 * t81 + t470 * t482 + t427;
t198 = t340 * t399;
t405 = t284 * t197 + t280 * t198;
t107 = qJ(4) * t376 + t405;
t445 = pkin(9) - pkin(10);
t484 = pkin(10) * t351 + t445 * t395 + t107;
t239 = t445 * t284;
t181 = t280 * t197;
t359 = -t198 * t284 + t181;
t411 = t284 * t285;
t483 = qJD(3) * t239 - (-pkin(10) * t411 - t281 * t446) * t399 - t359;
t426 = t282 * t30;
t480 = t278 * t470 ^ 2 - t117 * t79 - t426;
t479 = pkin(5) * t117;
t24 = t279 * t54 + t283 * t60;
t22 = pkin(11) * t227 + t24;
t25 = pkin(5) * t459 - pkin(11) * t117 + t59;
t10 = t22 * t282 + t25 * t278;
t478 = t10 * t117;
t477 = t470 * t117;
t334 = t22 * t278 - t25 * t282;
t476 = t117 * t334;
t475 = t117 * t459;
t216 = t279 * t280 + t283 * t284;
t137 = (qJD(3) - qJD(5)) * t216;
t412 = t277 * t285;
t167 = t216 * t412;
t156 = qJD(1) * t167;
t409 = t137 - t156;
t408 = -t490 * t279 + t280 * t392 + t481 * t283 - t284 * t393;
t472 = t490 * qJ(4) - t280 * qJD(4) - t200;
t211 = t281 * t442 + t285 * t342;
t471 = -g(1) * t211 - g(2) * t209 + g(3) * t412;
t469 = t117 ^ 2 - t459 ^ 2;
t465 = t227 * t459 + t33;
t132 = t206 * t279 + t207 * t283;
t313 = t279 * t16 + t283 * t19 + t54 * t392 - t393 * t60;
t89 = t147 * t279 + t148 * t283;
t464 = g(1) * t89 + g(2) * t85 + g(3) * t132 + t459 * t59 - t313;
t274 = t277 ^ 2;
t462 = 0.2e1 * t274;
t238 = t445 * t280;
t159 = t238 * t279 + t239 * t283;
t461 = -qJD(5) * t159 + t484 * t279 + t483 * t283;
t325 = t238 * t283 - t239 * t279;
t460 = -qJD(5) * t325 - t483 * t279 + t484 * t283;
t406 = t481 * t446 - t472;
t141 = qJD(3) * t207 + t277 * t372;
t373 = t277 * t396;
t142 = -qJD(3) * t206 + t284 * t373;
t379 = pkin(1) * t423;
t298 = pkin(8) * t412 + t281 * t379;
t202 = t298 * qJD(2);
t57 = t141 * pkin(3) - t142 * qJ(4) - t207 * qJD(4) + t202;
t458 = -t206 * pkin(3) + t207 * qJ(4);
t401 = t283 * qJ(4) - t279 * t446;
t457 = (qJDD(2) + 0.2e1 * t355) * t277;
t455 = t284 * pkin(3) + t280 * qJ(4) + pkin(2);
t439 = pkin(9) * t192;
t454 = t236 * t72 + t439;
t453 = t470 * (pkin(11) * t470 + t479) + t493;
t219 = -pkin(11) + t401;
t122 = t179 * pkin(3) + t177 * qJ(4);
t77 = -pkin(4) * t179 - t122;
t452 = t470 * (-pkin(11) * t459 + qJD(6) * t219 - t479 + t77) - t493;
t190 = pkin(9) * t423 + t298;
t360 = -t280 * t190 + t191 * t284;
t104 = pkin(3) * t412 - t360;
t69 = pkin(4) * t412 - pkin(10) * t207 + t104;
t407 = t284 * t190 + t280 * t191;
t103 = -qJ(4) * t412 + t407;
t74 = pkin(10) * t206 + t103;
t329 = t279 * t69 + t283 * t74;
t199 = t277 * t311;
t403 = -pkin(8) * t413 + t285 * t379;
t201 = t403 * qJD(2);
t315 = -t190 * t394 - t191 * t395 + t199 * t284 - t280 * t201;
t397 = qJD(2) * t281;
t374 = t277 * t397;
t45 = -pkin(10) * t142 - t374 * t446 - t315;
t307 = -t190 * t395 + t191 * t394 + t280 * t199 + t284 * t201;
t50 = qJ(4) * t374 - qJD(4) * t412 + t307;
t46 = pkin(10) * t141 + t50;
t450 = -qJD(5) * t329 - t279 * t46 + t283 * t45;
t3 = -pkin(11) * t188 + t313;
t353 = pkin(8) * t473 + t281 * t318 - t285 * t346;
t121 = -t317 * pkin(2) + t353;
t27 = t101 * pkin(3) + t100 * qJ(4) - t179 * qJD(4) + t121;
t20 = -pkin(4) * t101 - t27;
t6 = pkin(5) * t34 - pkin(11) * t33 + t20;
t1 = -t334 * qJD(6) + t278 * t6 + t282 * t3;
t447 = t179 ^ 2;
t287 = qJD(1) ^ 2;
t440 = pkin(3) * t192;
t434 = -pkin(5) * t376 - t461;
t433 = pkin(9) * qJD(3);
t78 = -t225 + t99;
t429 = t236 * t78;
t428 = t236 * t99;
t327 = -qJ(4) * t279 - t283 * t446;
t425 = -qJD(5) * t327 + t279 * t68 - t283 * t474;
t424 = t401 * qJD(5) + t279 * t474 + t283 * t68;
t422 = t177 * t236;
t421 = t179 * t177;
t420 = t179 * t236;
t324 = t279 * t284 - t280 * t283;
t416 = t324 * t278;
t415 = t324 * t282;
t414 = t274 * t287;
t404 = -t481 * pkin(3) + t472;
t275 = t281 ^ 2;
t400 = -t285 ^ 2 + t275;
t383 = t285 * t414;
t189 = -t423 * pkin(2) - t403;
t371 = pkin(1) * t462;
t358 = t283 * t227;
t213 = t284 * pkin(4) + t455;
t126 = pkin(5) * t216 + pkin(11) * t324 + t213;
t345 = -pkin(11) * t376 - qJD(6) * t126 + t460;
t344 = -pkin(5) * t408 + pkin(11) * t409 + qJD(6) * t159 - t406;
t341 = t277 * t287 * t423;
t339 = g(1) * t143 - g(2) * t147;
t338 = g(1) * t144 - g(2) * t148;
t337 = -g(1) * t209 + g(2) * t211;
t336 = -g(1) * t212 - g(2) * t210;
t39 = pkin(11) * t412 + t329;
t102 = t189 - t458;
t73 = -pkin(4) * t206 - t102;
t40 = -pkin(5) * t326 - pkin(11) * t132 + t73;
t333 = t278 * t40 + t282 * t39;
t332 = -t278 * t39 + t282 * t40;
t23 = -t279 * t60 + t283 * t54;
t330 = -t279 * t74 + t283 * t69;
t321 = 0.2e1 * t362 + qJD(2);
t314 = -t132 * t278 + t282 * t412;
t106 = t132 * t282 + t278 * t412;
t312 = t279 * t45 + t283 * t46 + t69 * t392 - t393 * t74;
t127 = t156 * t278 + t282 * t376;
t305 = t137 * t278 - t324 * t390 - t127;
t128 = t156 * t282 - t278 * t376;
t304 = t137 * t282 + t324 * t391 - t128;
t299 = -g(3) * t413 + t336;
t21 = -pkin(5) * t227 - t23;
t296 = -pkin(11) * t30 + (t21 + t23) * t470;
t295 = -t153 * t236 - t439;
t2 = -qJD(6) * t10 - t278 * t3 + t282 * t6;
t293 = -t219 * t30 + (-t21 + t425) * t470;
t292 = g(1) * t147 + g(2) * t143 + g(3) * t206 - t354;
t291 = -t236 * t433 + t471;
t290 = -t27 - t291;
t289 = t179 * t72 + qJDD(4) - t292;
t47 = -t141 * pkin(4) - t57;
t288 = g(1) * t148 + g(2) * t144 + g(3) * t207 - t236 * t98 - t308;
t218 = pkin(5) - t327;
t125 = t216 * t211;
t124 = t216 * t209;
t110 = -pkin(3) * t376 + t359;
t76 = pkin(3) * t236 + t410;
t63 = -t211 * t278 + t282 * t89;
t62 = -t211 * t282 - t278 * t89;
t61 = -t100 - t422;
t56 = -pkin(3) * t374 - t315;
t53 = qJD(5) * t326 + t141 * t279 + t142 * t283;
t52 = qJD(5) * t132 - t141 * t283 + t142 * t279;
t38 = -pkin(5) * t412 - t330;
t32 = qJD(6) * t314 - t278 * t374 + t282 * t53;
t31 = qJD(6) * t106 + t278 * t53 + t282 * t374;
t28 = t320 - t440;
t11 = t52 * pkin(5) - t53 * pkin(11) + t47;
t8 = pkin(5) * t374 - t450;
t7 = -pkin(11) * t374 + t312;
t5 = [qJDD(1), g(1) * t441 - g(2) * t442, g(1) * t442 + g(2) * t441 (qJDD(1) * t275 + 0.2e1 * t281 * t369) * t274 (t285 * t385 - t388 * t400) * t462, t281 * t457 + t321 * t373, t285 * t457 - t321 * t374, t317 * t423, -t202 * t322 + t403 * t317 - t353 * t423 + g(1) * t210 - g(2) * t212 + (-t370 + t387) * t371, -t201 * t322 - t298 * t317 - t297 * t423 + (-t369 - t385) * t371 + t337, -t100 * t207 + t142 * t179, t100 * t206 - t101 * t207 - t141 * t179 - t142 * t177, -t142 * t236 + t192 * t207 + (t100 * t285 + t179 * t397) * t277, t141 * t236 - t192 * t206 + (t101 * t285 - t177 * t397) * t277 (-t192 * t285 - t236 * t397) * t277, -t315 * t236 + t360 * t192 + t202 * t177 + t189 * t101 + t121 * t206 + t153 * t141 + (t285 * t354 + t397 * t98) * t277 + t338, t307 * t236 - t407 * t192 + t202 * t179 - t189 * t100 + t121 * t207 + t153 * t142 + (t285 * t308 - t397 * t99) * t277 - t339, t101 * t102 - t104 * t192 + t141 * t72 + t177 * t57 + t206 * t27 + t236 * t56 + (t28 * t285 - t397 * t76) * t277 + t338, -t100 * t104 - t101 * t103 - t141 * t78 + t142 * t76 - t177 * t50 + t179 * t56 - t206 * t26 + t207 * t28 - t337, t100 * t102 + t103 * t192 - t142 * t72 - t179 * t57 - t207 * t27 - t236 * t50 + (-t26 * t285 + t397 * t78) * t277 + t339, t26 * t103 + t78 * t50 + t27 * t102 + t72 * t57 + t28 * t104 + t76 * t56 - g(1) * (-pkin(1) * t441 - t210 * pkin(2) - pkin(3) * t144 + pkin(8) * t378 - t209 * pkin(9) - qJ(4) * t143) - g(2) * (pkin(1) * t442 + t212 * pkin(2) + t148 * pkin(3) + pkin(8) * t377 + t211 * pkin(9) + t147 * qJ(4)) t117 * t53 + t132 * t33, -t117 * t52 - t132 * t34 + t326 * t33 - t459 * t53, -t132 * t188 + t227 * t53 + (-t117 * t397 + t285 * t33) * t277, -t326 * t188 - t227 * t52 + (-t285 * t34 + t397 * t459) * t277 (-t188 * t285 - t227 * t397) * t277, t450 * t227 - t330 * t188 + t47 * t459 + t73 * t34 - t20 * t326 + t59 * t52 + g(1) * t85 - g(2) * t89 + (-t23 * t397 - t285 * t365) * t277, -t312 * t227 + t329 * t188 + t47 * t117 + t73 * t33 + t20 * t132 + t59 * t53 + g(1) * t463 - g(2) * t88 + (t24 * t397 - t285 * t313) * t277, t106 * t12 + t32 * t81, -t106 * t13 + t12 * t314 - t31 * t81 - t32 * t79, t106 * t30 - t12 * t326 + t32 * t470 + t52 * t81, t13 * t326 + t30 * t314 - t31 * t470 - t52 * t79, -t30 * t326 + t470 * t52 (-qJD(6) * t333 + t11 * t282 - t278 * t7) * t470 + t332 * t30 - t2 * t326 - t334 * t52 + t8 * t79 + t38 * t13 - t4 * t314 + t21 * t31 + g(1) * t487 - g(2) * t63 -(qJD(6) * t332 + t11 * t278 + t282 * t7) * t470 - t333 * t30 + t1 * t326 - t10 * t52 + t8 * t81 + t38 * t12 + t4 * t106 + t21 * t32 - g(1) * t488 - g(2) * t62; 0, 0, 0, -t281 * t383, t400 * t414, -t285 * t341 + t368, t281 * t341 + t257, t317, pkin(1) * t281 * t414 + t200 * t322 - t353 - t471, pkin(1) * t383 + t197 * t322 + (pkin(8) * t388 + g(3)) * t413 - t336 - t380, -t100 * t280 - t284 * t420 (-t100 + t422) * t284 + (-t101 + t420) * t280, -t236 * t394 + t192 * t280 + (-t179 * t281 + t236 * t411) * t399, t236 * t395 + t192 * t284 + (-t236 * t280 * t285 + t177 * t281) * t399, t236 * t376, -t98 * t376 - pkin(2) * t101 - t200 * t177 - t181 * t236 + t295 * t280 + (-t121 + (t198 + t433) * t236 - t471) * t284, pkin(2) * t100 - t405 * t236 + t99 * t376 - t200 * t179 + t295 * t284 + (t121 + t291) * t280, -t101 * t455 - t110 * t236 + t404 * t177 - t280 * t454 + t290 * t284 + t76 * t376, t107 * t177 - t110 * t179 + (t26 - t236 * t76 + (qJD(3) * t179 - t101) * pkin(9)) * t284 + (t28 + t429 + (qJD(3) * t177 - t100) * pkin(9)) * t280 + t299, -t100 * t455 + t107 * t236 - t404 * t179 + t290 * t280 + t284 * t454 - t78 * t376, -t78 * t107 - t76 * t110 + t404 * t72 + (t26 * t284 + t28 * t280 + (-t280 * t78 + t284 * t76) * qJD(3) + t299) * pkin(9) + (-t27 - t471) * t455, t117 * t409 - t324 * t33, -t117 * t408 - t216 * t33 + t324 * t34 - t409 * t459, t117 * t376 + t188 * t324 + t227 * t409, t188 * t216 - t227 * t408 - t376 * t459, t227 * t376, g(1) * t125 + g(2) * t124 - g(3) * t167 - t188 * t325 + t20 * t216 + t213 * t34 + t227 * t461 + t23 * t376 + t406 * t459 + t408 * t59, t406 * t117 + t159 * t188 + t213 * t33 + t460 * t227 - t24 * t376 + t409 * t59 + (-t20 + t471) * t324, -t12 * t415 + t304 * t81, t127 * t81 + t128 * t79 + (-t278 * t81 - t282 * t79) * t137 - (-t430 - t13 * t282 + (t278 * t79 - t282 * t81) * qJD(6)) * t324, t12 * t216 - t30 * t415 + t304 * t470 + t408 * t81, -t13 * t216 + t30 * t416 - t305 * t470 - t408 * t79, t216 * t30 + t408 * t470 (t126 * t282 - t159 * t278) * t30 + t2 * t216 - t325 * t13 - t4 * t416 - g(1) * (-t125 * t282 - t212 * t278) - g(2) * (-t124 * t282 - t210 * t278) - g(3) * (t167 * t282 - t278 * t413) - t408 * t334 + t434 * t79 + (t278 * t345 - t282 * t344) * t470 + t305 * t21 -(t126 * t278 + t159 * t282) * t30 - t1 * t216 - t325 * t12 - t4 * t415 - g(1) * (t125 * t278 - t212 * t282) - g(2) * (t124 * t278 - t210 * t282) - g(3) * (-t167 * t278 - t282 * t413) + t434 * t81 + (t278 * t344 + t282 * t345) * t470 - t408 * t10 + t304 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t421, -t177 ^ 2 + t447, t61, -t101 - t420, t192, -t153 * t179 + t292 - t428, t153 * t177 + t288, -t122 * t177 - t289 - t428 + 0.2e1 * t440, pkin(3) * t100 - qJ(4) * t101 + (t78 - t99) * t179 + (t76 - t410) * t177, t122 * t179 - t177 * t72 + 0.2e1 * t187 - 0.2e1 * t223 - t288, t26 * qJ(4) - t28 * pkin(3) - t72 * t122 - t76 * t99 - g(1) * (-pkin(3) * t147 + qJ(4) * t148) - g(2) * (-pkin(3) * t143 + qJ(4) * t144) - g(3) * t458 + t410 * t78, -t475, -t469, -t465, -t117 * t227 + t34, t188, -t188 * t327 - t227 * t424 - t459 * t77 - t489, -t77 * t117 + t188 * t401 + t227 * t425 - t464, -t486, t494, -t485, t480, t477, t218 * t13 + t293 * t278 - t282 * t452 + t424 * t79 - t476, t218 * t12 + t278 * t452 + t293 * t282 + t424 * t81 - t478; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192 + t421, t61, -t236 ^ 2 - t447, t289 + t429 - t440, 0, 0, 0, 0, 0, -t227 ^ 2 * t279 - t179 * t459 - t188 * t283, -t117 * t179 + t188 * t279 - t227 * t358, 0, 0, 0, 0, 0, -t283 * t13 + (-t282 * t179 - t278 * t358) * t470 + (t227 * t79 - t390 * t470 - t427) * t279, -t283 * t12 + (t278 * t179 - t282 * t358) * t470 + (t227 * t81 + t391 * t470 - t426) * t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t475, t469, t465, -t364 + (-qJD(5) + t227) * t117, -t188, t227 * t24 + t489, t227 * t23 + t464, t486, -t494, t485, -t480, -t477, -pkin(5) * t13 - t24 * t79 + t296 * t278 - t282 * t453 + t476, -pkin(5) * t12 - t24 * t81 + t278 * t453 + t296 * t282 + t478; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t79, -t79 ^ 2 + t81 ^ 2, t12 + t432, -t13 + t431, t30, -g(1) * t62 + g(2) * t488 - g(3) * t314 + t10 * t470 - t21 * t81 + t2, g(1) * t63 + g(2) * t487 + g(3) * t106 + t21 * t79 - t334 * t470 - t1;];
tau_reg  = t5;