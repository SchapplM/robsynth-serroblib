% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:35:48
% EndTime: 2019-03-10 03:36:18
% DurationCPUTime: 13.48s
% Computational Cost: add. (33427->657), mult. (82801->856), div. (0->0), fcn. (63003->10), ass. (0->326)
t336 = cos(qJ(2));
t497 = cos(qJ(3));
t410 = t497 * t336;
t385 = qJD(1) * t410;
t333 = sin(qJ(2));
t495 = sin(qJ(3));
t408 = t495 * t333;
t271 = -qJD(1) * t408 + t385;
t287 = t333 * t497 + t336 * t495;
t272 = qJD(1) * t287;
t332 = sin(qJ(4));
t496 = cos(qJ(4));
t222 = -t496 * t271 + t272 * t332;
t334 = cos(qJ(6));
t335 = cos(qJ(5));
t422 = qJD(5) + qJD(6);
t425 = qJD(6) * t334;
t427 = qJD(5) * t335;
t330 = sin(qJ(6));
t331 = sin(qJ(5));
t445 = t330 * t331;
t507 = t334 * t335 - t445;
t545 = -t507 * t222 - t334 * t427 - t335 * t425 + t422 * t445;
t444 = t330 * t335;
t286 = t331 * t334 + t444;
t242 = t422 * t286;
t544 = t286 * t222 + t242;
t526 = t222 * t335;
t543 = t427 + t526;
t327 = qJD(2) + qJD(3);
t325 = qJD(4) + t327;
t368 = t332 * t271 + t272 * t496;
t201 = -t335 * t325 + t331 * t368;
t499 = -pkin(8) - pkin(7);
t305 = t499 * t336;
t293 = qJD(1) * t305;
t273 = t495 * t293;
t303 = t499 * t333;
t291 = qJD(1) * t303;
t484 = qJD(2) * pkin(2);
t279 = t291 + t484;
t232 = t497 * t279 + t273;
t266 = t272 * pkin(9);
t199 = -t266 + t232;
t187 = t327 * pkin(3) + t199;
t277 = t497 * t293;
t233 = t279 * t495 - t277;
t492 = t271 * pkin(9);
t200 = t233 + t492;
t195 = t496 * t200;
t128 = t332 * t187 + t195;
t124 = t325 * pkin(10) + t128;
t321 = -t336 * pkin(2) - pkin(1);
t301 = qJD(1) * t321;
t243 = -pkin(3) * t271 + t301;
t146 = pkin(4) * t222 - pkin(10) * t368 + t243;
t80 = t124 * t335 + t146 * t331;
t65 = -pkin(11) * t201 + t80;
t480 = t330 * t65;
t522 = qJD(5) + t222;
t203 = t325 * t331 + t335 * t368;
t79 = -t124 * t331 + t335 * t146;
t64 = -pkin(11) * t203 + t79;
t55 = pkin(5) * t522 + t64;
t23 = t334 * t55 - t480;
t478 = t334 * t65;
t24 = t330 * t55 + t478;
t412 = qJD(2) * t499;
t386 = qJD(1) * t412;
t280 = t333 * t386;
t281 = t336 * t386;
t402 = t495 * qJD(3);
t404 = qJD(3) * t497;
t175 = t279 * t404 + t497 * t280 + t495 * t281 + t293 * t402;
t346 = t327 * t287;
t343 = t346 * qJD(1);
t132 = -pkin(9) * t343 + t175;
t379 = -t280 * t495 + t497 * t281;
t176 = -qJD(3) * t233 + t379;
t430 = qJD(1) * t333;
t229 = (qJD(2) * t495 + t402) * t430 - t327 * t385;
t133 = t229 * pkin(9) + t176;
t403 = qJD(4) * t496;
t429 = qJD(4) * t332;
t349 = -t132 * t496 - t332 * t133 - t187 * t403 + t200 * t429;
t428 = qJD(5) * t331;
t342 = t496 * t346;
t125 = qJD(1) * t342 + qJD(4) * t368 - t332 * t229;
t423 = qJD(1) * qJD(2);
t401 = t333 * t423;
t313 = pkin(2) * t401;
t206 = pkin(3) * t343 + t313;
t340 = -t229 * t496 + t271 * t403 - t272 * t429 - t332 * t343;
t63 = t125 * pkin(4) - pkin(10) * t340 + t206;
t13 = -t124 * t428 + t146 * t427 + t331 * t63 - t335 * t349;
t413 = t325 * t428 + t331 * t340 + t368 * t427;
t10 = -pkin(11) * t413 + t13;
t426 = qJD(6) * t330;
t101 = -t325 * t427 - t335 * t340 + t368 * t428;
t14 = -qJD(5) * t80 + t331 * t349 + t335 * t63;
t9 = pkin(5) * t125 + pkin(11) * t101 + t14;
t3 = (qJD(6) * t55 + t10) * t334 + t330 * t9 - t65 * t426;
t4 = -qJD(6) * t24 - t10 * t330 + t334 * t9;
t359 = t23 * t545 - t544 * t24 - t4 * t286 + t3 * t507;
t194 = t332 * t200;
t127 = t187 * t496 - t194;
t123 = -t325 * pkin(4) - t127;
t102 = t201 * pkin(5) + t123;
t393 = t332 * t132 - t496 * t133;
t52 = qJD(4) * t128 + t393;
t35 = pkin(5) * t413 + t52;
t363 = -t102 * t545 + t24 * t368 + t35 * t286;
t138 = t334 * t201 + t203 * t330;
t370 = t201 * t330 - t334 * t203;
t353 = qJD(6) * t370 + t101 * t330 - t334 * t413;
t44 = t334 * t101 + t201 * t425 + t203 * t426 + t330 * t413;
t7 = t138 * t545 + t286 * t353 + t544 * t370 - t44 * t507;
t537 = t522 * t331;
t95 = t331 * t413;
t20 = -t101 * t335 - t201 * t543 - t203 * t537 - t95;
t213 = qJD(6) + t522;
t36 = t286 * t125 - t213 * t545 + t368 * t370;
t16 = -t44 * t286 + t370 * t545;
t122 = t335 * t125;
t53 = t201 * t368 - t522 * t537 + t122;
t527 = t222 * t331;
t420 = pkin(11) * t527;
t388 = pkin(3) * t403;
t135 = t199 * t496 - t194;
t172 = pkin(4) * t368 + pkin(10) * t222;
t494 = pkin(3) * t272;
t154 = t172 + t494;
t82 = -t135 * t331 + t335 * t154;
t542 = -t331 * t388 - t82;
t418 = t497 * pkin(2);
t319 = t418 + pkin(3);
t409 = t495 * t332;
t227 = t319 * t403 + (-qJD(4) * t409 + (t496 * t497 - t409) * qJD(3)) * pkin(2);
t240 = t497 * t291 + t273;
t205 = -t266 + t240;
t239 = -t291 * t495 + t277;
t362 = t239 - t492;
t148 = t205 * t496 + t332 * t362;
t322 = pkin(2) * t430;
t149 = t154 + t322;
t84 = -t148 * t331 + t335 * t149;
t541 = -t227 * t331 - t84;
t83 = t335 * t135 + t331 * t154;
t540 = -t335 * t388 + t83;
t85 = t335 * t148 + t331 * t149;
t539 = -t227 * t335 + t85;
t508 = (t428 + t527) * pkin(5);
t383 = pkin(5) * t368 + pkin(11) * t526;
t364 = t102 * t544 - t23 * t368 - t35 * t507;
t17 = t138 * t544 + t353 * t507;
t98 = t101 * t331;
t56 = t203 * t543 - t98;
t37 = t125 * t507 + t138 * t368 - t213 * t544;
t120 = t331 * t125;
t509 = -t427 * t522 - t120;
t54 = -t203 * t368 + t522 * t526 - t509;
t387 = t495 * t496;
t268 = pkin(2) * t387 + t332 * t319;
t264 = pkin(10) + t268;
t490 = -pkin(11) - t264;
t399 = qJD(5) * t490;
t533 = t331 * t399 - t420 - t539;
t532 = -t335 * t399 + t383 - t541;
t493 = pkin(3) * t332;
t317 = pkin(10) + t493;
t489 = -pkin(11) - t317;
t398 = qJD(5) * t489;
t531 = t331 * t398 - t420 - t540;
t530 = t335 * t398 - t383 + t542;
t498 = -pkin(11) - pkin(10);
t411 = qJD(5) * t498;
t91 = t335 * t127 + t331 * t172;
t529 = t331 * t411 - t420 - t91;
t90 = -t127 * t331 + t335 * t172;
t528 = t335 * t411 - t383 - t90;
t467 = t138 * t370;
t525 = t222 * t368;
t448 = t243 * t368;
t523 = -t52 - t448;
t106 = -t222 ^ 2 + t368 ^ 2;
t521 = -t138 ^ 2 + t370 ^ 2;
t520 = t102 * t370 + t4;
t519 = t138 * t213 - t44;
t518 = t102 * t138 - t3;
t517 = -t213 * t370 + t353;
t348 = t243 * t222 + t349;
t103 = t222 * t325 + t340;
t515 = 0.2e1 * t301;
t514 = -0.2e1 * t423;
t511 = -t522 * t80 - t14;
t462 = t213 * t368;
t461 = t522 * t368;
t134 = t199 * t332 + t195;
t381 = pkin(3) * t429 - t134;
t247 = t497 * t303 + t305 * t495;
t216 = -t287 * pkin(9) + t247;
t249 = t495 * t303 - t497 * t305;
t361 = t408 - t410;
t217 = -pkin(9) * t361 + t249;
t169 = t332 * t216 + t217 * t496;
t162 = t335 * t169;
t354 = t496 * t361;
t237 = t287 * t332 + t354;
t357 = t332 * t361;
t238 = t287 * t496 - t357;
t255 = pkin(3) * t361 + t321;
t167 = t237 * pkin(4) - t238 * pkin(10) + t255;
t93 = t331 * t167 + t162;
t434 = -t227 + t148;
t433 = -t205 * t332 + t496 * t362 + t319 * t429 + (qJD(4) * t387 + (t332 * t497 + t387) * qJD(3)) * pkin(2);
t506 = -t331 * t79 + t335 * t80;
t369 = t123 * t428 - t52 * t335 - t79 * t368;
t382 = t123 * t427 + t52 * t331 + t368 * t80;
t104 = t325 * t368 - t125;
t491 = t335 * pkin(5);
t244 = t490 * t331;
t326 = t335 * pkin(11);
t245 = t264 * t335 + t326;
t186 = t244 * t330 + t245 * t334;
t488 = qJD(6) * t186 + t330 * t533 + t532 * t334;
t185 = t244 * t334 - t245 * t330;
t487 = -qJD(6) * t185 + t532 * t330 - t334 * t533;
t12 = t13 * t335;
t168 = -t496 * t216 + t217 * t332;
t483 = t168 * t52;
t482 = t522 * t79;
t282 = t489 * t331;
t283 = t317 * t335 + t326;
t230 = t282 * t334 - t283 * t330;
t473 = qJD(6) * t230 + t330 * t530 + t334 * t531;
t231 = t282 * t330 + t283 * t334;
t472 = -qJD(6) * t231 - t330 * t531 + t334 * t530;
t302 = t498 * t331;
t304 = pkin(10) * t335 + t326;
t246 = t302 * t334 - t304 * t330;
t471 = qJD(6) * t246 + t330 * t528 + t334 * t529;
t248 = t302 * t330 + t304 * t334;
t470 = -qJD(6) * t248 - t330 * t529 + t334 * t528;
t468 = t123 * t222;
t107 = t125 * t237;
t344 = t327 * t361;
t150 = qJD(4) * t354 + t287 * t429 + t332 * t346 + t344 * t496;
t466 = t150 * t331;
t465 = t150 * t335;
t464 = t201 * t331;
t463 = t203 * t201;
t450 = t238 * t331;
t449 = t238 * t335;
t447 = t272 * t271;
t446 = t301 * t272;
t338 = qJD(1) ^ 2;
t442 = t336 * t338;
t337 = qJD(2) ^ 2;
t441 = t337 * t333;
t440 = t337 * t336;
t439 = t508 + t381;
t438 = t508 + t433;
t431 = t333 ^ 2 - t336 ^ 2;
t421 = -t526 * t79 - t527 * t80 + t12;
t417 = t496 * pkin(3);
t416 = t495 * pkin(2);
t324 = t333 * t484;
t414 = t333 * t442;
t407 = t238 * t428;
t292 = t333 * t412;
t294 = t336 * t412;
t183 = t497 * t292 + t495 * t294 + t303 * t404 + t305 * t402;
t161 = -pkin(9) * t346 + t183;
t378 = -t292 * t495 + t497 * t294;
t339 = pkin(9) * t344 - t303 * t402 + t305 * t404 + t378;
t68 = t161 * t496 + t216 * t403 - t217 * t429 + t332 * t339;
t151 = -qJD(4) * t357 + t287 * t403 - t332 * t344 + t342;
t226 = pkin(3) * t346 + t324;
t77 = t151 * pkin(4) + t150 * pkin(10) + t226;
t400 = -t331 * t68 + t335 * t77;
t394 = pkin(1) * t514;
t92 = t335 * t167 - t169 * t331;
t318 = -t417 - pkin(4);
t384 = t336 * t401;
t380 = -t128 + t508;
t100 = t413 * t335;
t267 = -pkin(2) * t409 + t319 * t496;
t74 = pkin(5) * t237 - pkin(11) * t449 + t92;
t81 = -pkin(11) * t450 + t93;
t42 = -t330 * t81 + t334 * t74;
t43 = t330 * t74 + t334 * t81;
t377 = t331 * t80 + t335 * t79;
t373 = -t125 * t264 + t468;
t372 = -t125 * t317 + t468;
t371 = -t127 * t222 + t128 * t368;
t263 = -pkin(4) - t267;
t367 = t238 * t427 - t466;
t366 = -t407 - t465;
t18 = t167 * t427 - t169 * t428 + t331 * t77 + t335 * t68;
t365 = -t301 * t271 - t175;
t356 = t361 * qJD(2);
t352 = -qJD(5) * t377 - t14 * t331 + t12;
t69 = qJD(4) * t169 + t332 * t161 - t496 * t339;
t320 = -pkin(4) - t491;
t298 = t318 - t491;
t253 = t263 - t491;
t252 = t322 + t494;
t204 = -t271 ^ 2 + t272 ^ 2;
t190 = t272 * t327 - t343;
t189 = -t271 * t327 - t229;
t184 = -qJD(3) * t249 + t378;
t174 = t507 * t238;
t173 = t286 * t238;
t118 = pkin(5) * t450 + t168;
t57 = t464 * t522 - t100;
t48 = -t150 * t444 - t330 * t407 - t426 * t450 + (t422 * t449 - t466) * t334;
t47 = t150 * t507 + t238 * t242;
t46 = pkin(5) * t367 + t69;
t26 = t334 * t64 - t480;
t25 = -t330 * t64 - t478;
t19 = -qJD(5) * t93 + t400;
t15 = -pkin(11) * t367 + t18;
t11 = pkin(11) * t465 + pkin(5) * t151 + (-t162 + (pkin(11) * t238 - t167) * t331) * qJD(5) + t400;
t6 = -qJD(6) * t43 + t11 * t334 - t15 * t330;
t5 = qJD(6) * t42 + t11 * t330 + t15 * t334;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t384, t431 * t514, t440, -0.2e1 * t384, -t441, 0, -pkin(7) * t440 + t333 * t394, pkin(7) * t441 + t336 * t394, 0, 0, -t229 * t287 - t272 * t344, t229 * t361 - t287 * t343 + t327 * (-t271 * t361 - t272 * t287) -t344 * t327, -t271 * t346 + t343 * t361, -t346 * t327, 0, -t271 * t324 + t322 * t356 + (t287 * t515 + t184) * t327, -t183 * t327 - t321 * t229 + t272 * t324 + t287 * t313 - t301 * t344, t183 * t271 - t249 * t343 - t175 * t361 - t233 * t346 - t184 * t272 + t247 * t229 - t176 * t287 + t232 * (qJD(3) * t361 + t356) t175 * t249 + t176 * t247 + t183 * t233 + t184 * t232 + t324 * t515, -t368 * t150 + t238 * t340, -t238 * t125 + t150 * t222 - t151 * t368 - t237 * t340, -t150 * t325, t151 * t222 + t107, -t151 * t325, 0, t125 * t255 + t151 * t243 + t206 * t237 + t222 * t226 - t325 * t69, -t243 * t150 + t206 * t238 + t226 * t368 + t255 * t340 - t68 * t325, -t169 * t125 + t127 * t150 - t128 * t151 + t168 * t340 - t68 * t222 + t237 * t349 + t52 * t238 + t368 * t69, -t127 * t69 + t128 * t68 - t169 * t349 + t206 * t255 + t226 * t243 + t483, -t101 * t449 + t203 * t366 (t201 * t335 + t203 * t331) * t150 + (-t100 + t98 + (-t203 * t335 + t464) * qJD(5)) * t238, -t101 * t237 + t122 * t238 + t151 * t203 + t366 * t522, t201 * t367 + t238 * t95, -t120 * t238 - t201 * t151 - t237 * t413 - t367 * t522, t151 * t522 + t107, t123 * t367 + t92 * t125 + t14 * t237 + t79 * t151 + t168 * t413 + t19 * t522 + t69 * t201 + t450 * t52, -t101 * t168 + t123 * t366 - t125 * t93 - t13 * t237 - t151 * t80 - t18 * t522 + t203 * t69 + t449 * t52, -t18 * t201 - t93 * t413 - t19 * t203 + t92 * t101 + t377 * t150 + (-qJD(5) * t506 - t13 * t331 - t14 * t335) * t238, t123 * t69 + t13 * t93 + t14 * t92 + t18 * t80 + t19 * t79 + t483, -t174 * t44 + t370 * t47, t138 * t47 + t173 * t44 + t174 * t353 + t370 * t48, t125 * t174 - t151 * t370 - t213 * t47 - t237 * t44, t138 * t48 - t173 * t353, -t125 * t173 - t138 * t151 - t213 * t48 + t237 * t353, t151 * t213 + t107, t102 * t48 - t118 * t353 + t125 * t42 + t138 * t46 + t151 * t23 + t173 * t35 + t213 * t6 + t237 * t4, -t102 * t47 - t118 * t44 - t125 * t43 - t151 * t24 + t174 * t35 - t213 * t5 - t237 * t3 - t370 * t46, -t138 * t5 - t173 * t3 - t174 * t4 + t23 * t47 - t24 * t48 + t353 * t43 + t370 * t6 + t42 * t44, t102 * t46 + t118 * t35 + t23 * t6 + t24 * t5 + t3 * t43 + t4 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t414, t431 * t338, 0, t414, 0, 0, t338 * pkin(1) * t333, pkin(1) * t442, 0, 0, -t447, t204, t189, t447, t190, 0, t271 * t322 - t239 * t327 - t446 + (-t327 * t416 - t233) * qJD(3) + t379, t240 * t327 + (-t272 * t430 - t327 * t404) * pkin(2) + t365, t229 * t418 - t343 * t416 + (pkin(2) * t402 + t233 + t239) * t272 + (pkin(2) * t404 + t232 - t240) * t271, -t232 * t239 - t233 * t240 + (-t301 * t430 + t497 * t176 + t495 * t175 + (-t232 * t495 + t233 * t497) * qJD(3)) * pkin(2), t525, t106, t103, -t525, t104, 0, -t252 * t222 - t325 * t433 + t523, -t252 * t368 + t325 * t434 + t348, -t268 * t125 + t222 * t434 - t267 * t340 + t368 * t433 + t371, -t127 * t433 - t128 * t434 - t243 * t252 - t267 * t52 - t268 * t349, t56, t20, t54, t57, t53, -t461, t263 * t413 + t373 * t331 + t433 * t201 + (-t264 * t427 + t541) * t522 + t369, -t101 * t263 + t373 * t335 + t433 * t203 + (t264 * t428 + t539) * t522 + t382, t85 * t201 + t84 * t203 + (-t227 * t201 - t264 * t413 + (t203 * t264 - t79) * qJD(5)) * t335 + (-t264 * t101 + t227 * t203 - t14 + (t201 * t264 - t80) * qJD(5)) * t331 + t421, t123 * t433 + t227 * t506 + t263 * t52 + t264 * t352 - t79 * t84 - t80 * t85, t16, t7, t36, t17, t37, -t462, t125 * t185 + t138 * t438 - t213 * t488 - t253 * t353 + t364, -t125 * t186 + t213 * t487 - t253 * t44 - t370 * t438 + t363, t138 * t487 + t185 * t44 + t186 * t353 - t370 * t488 + t359, t102 * t438 + t185 * t4 + t186 * t3 - t23 * t488 - t24 * t487 + t253 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t447, t204, t189, t447, t190, 0, t233 * t327 + t176 - t446, t232 * t327 + t365, 0, 0, t525, t106, t103, -t525, t104, 0, -t222 * t494 + t134 * t325 - t448 + (-t195 + (-pkin(3) * t325 - t187) * t332) * qJD(4) - t393, t135 * t325 + (-t272 * t368 - t325 * t403) * pkin(3) + t348, -t125 * t493 - t340 * t417 + t371 + t381 * t368 + (t135 - t388) * t222, t127 * t134 - t128 * t135 + (-t496 * t52 - t243 * t272 - t332 * t349 + (-t127 * t332 + t128 * t496) * qJD(4)) * pkin(3), t56, t20, t54, t57, t53, -t461, t318 * t413 + t372 * t331 + t381 * t201 + (-t317 * t427 + t542) * t522 + t369, -t318 * t101 + t372 * t335 + t381 * t203 + (t317 * t428 + t540) * t522 + t382, t83 * t201 + t82 * t203 + (-t201 * t388 - t317 * t413 + (t203 * t317 - t79) * qJD(5)) * t335 + (t203 * t388 - t317 * t101 - t14 + (t201 * t317 - t80) * qJD(5)) * t331 + t421, -t123 * t134 + t52 * t318 - t79 * t82 - t80 * t83 + (t123 * t332 + t496 * t506) * qJD(4) * pkin(3) + t352 * t317, t16, t7, t36, t17, t37, -t462, t125 * t230 + t138 * t439 + t213 * t472 - t298 * t353 + t364, -t125 * t231 - t213 * t473 - t298 * t44 - t370 * t439 + t363, -t138 * t473 + t230 * t44 + t231 * t353 + t370 * t472 + t359, t102 * t439 + t23 * t472 + t230 * t4 + t231 * t3 + t24 * t473 + t298 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t525, t106, t103, -t525, t104, 0, t128 * t325 + t523, t127 * t325 + t348, 0, 0, t56, t20, t54, t201 * t537 - t100, t53, -t461, -pkin(4) * t413 + pkin(10) * t509 + t123 * t527 - t128 * t201 - t522 * t90 + t369, t123 * t526 + pkin(4) * t101 - t128 * t203 + t522 * t91 + (t428 * t522 - t122) * pkin(10) + t382, t91 * t201 + t90 * t203 + t12 + (-t482 + (qJD(5) * t203 - t413) * pkin(10)) * t335 + ((qJD(5) * t201 - t101) * pkin(10) + t511) * t331, -pkin(4) * t52 + pkin(10) * t352 - t123 * t128 - t79 * t90 - t80 * t91, t16, t7, t36, t17, t37, -t462, t125 * t246 + t138 * t380 + t213 * t470 - t320 * t353 + t364, -t125 * t248 - t213 * t471 - t320 * t44 - t370 * t380 + t363, -t138 * t471 + t246 * t44 + t248 * t353 + t370 * t470 + t359, t102 * t380 + t23 * t470 + t24 * t471 + t246 * t4 + t248 * t3 + t320 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t463, -t201 ^ 2 + t203 ^ 2, t201 * t522 - t101, -t463, t203 * t522 - t413, t125, -t123 * t203 - t511, t123 * t201 - t13 + t482, 0, 0, -t467, t521, t519, t467, t517, t125, -t213 * t25 + (t125 * t334 - t138 * t203 - t213 * t426) * pkin(5) + t520, t213 * t26 + (-t125 * t330 + t203 * t370 - t213 * t425) * pkin(5) + t518, -t370 * t24 + t138 * t26 - t138 * t23 - t370 * t25 + (t330 * t353 + t334 * t44 + (-t138 * t334 - t330 * t370) * qJD(6)) * pkin(5), -t23 * t25 - t24 * t26 + (-t102 * t203 + t3 * t330 + t334 * t4 + (-t23 * t330 + t24 * t334) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t467, t521, t519, t467, t517, t125, t213 * t24 + t520, t213 * t23 + t518, 0, 0;];
tauc_reg  = t1;
