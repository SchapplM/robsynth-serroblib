% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRRRP11
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 07:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRRRP11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 07:08:44
% EndTime: 2019-05-08 07:10:03
% DurationCPUTime: 30.97s
% Computational Cost: add. (299292->798), mult. (739511->1214), div. (0->0), fcn. (628066->14), ass. (0->531)
t376 = cos(pkin(6));
t524 = qJD(1) * t376;
t368 = qJD(2) + t524;
t516 = qJD(2) + t368;
t379 = sin(qJ(3));
t380 = sin(qJ(2));
t383 = cos(qJ(3));
t374 = sin(pkin(6));
t525 = qJD(1) * t374;
t375 = cos(pkin(7));
t384 = cos(qJ(2));
t542 = t375 * t384;
t373 = sin(pkin(7));
t548 = t373 * t379;
t343 = t368 * t548 + (t379 * t542 + t380 * t383) * t525;
t522 = qJD(1) * t384;
t499 = t374 * t522;
t353 = t368 * t375 - t373 * t499 + qJD(3);
t378 = sin(qJ(4));
t382 = cos(qJ(4));
t330 = t343 * t382 + t353 * t378;
t377 = sin(qJ(5));
t381 = cos(qJ(5));
t489 = t375 * t499;
t523 = qJD(1) * t380;
t500 = t374 * t523;
t546 = t373 * t383;
t341 = -t368 * t546 + t379 * t500 - t383 * t489;
t420 = qJD(4) + t341;
t308 = t330 * t377 - t381 * t420;
t310 = t330 * t381 + t377 * t420;
t270 = t310 * t308;
t511 = qJDD(1) * t384;
t513 = qJD(1) * qJD(2);
t433 = -t380 * t513 + t511;
t414 = t433 * t375;
t494 = t384 * t513;
t512 = qJDD(1) * t380;
t434 = t494 + t512;
t488 = qJDD(1) * t376 + qJDD(2);
t442 = t373 * t488;
t391 = (t379 * t414 + t383 * t434) * t374 + t379 * t442;
t320 = -t341 * qJD(3) + t391;
t415 = t433 * t374;
t401 = t373 * t415 - t375 * t488 - qJDD(3);
t492 = t378 * t320 + t382 * t401;
t276 = -t330 * qJD(4) - t492;
t275 = qJDD(5) - t276;
t613 = -t270 + t275;
t623 = pkin(5) * t613;
t328 = t343 * t378 - t353 * t382;
t448 = -t320 * t382 + t378 * t401;
t277 = -qJD(4) * t328 - t448;
t416 = t434 * t374;
t402 = t379 * t416 + (-t374 * t414 - t442) * t383;
t397 = t343 * qJD(3) + t402;
t395 = qJDD(4) + t397;
t227 = -t308 * qJD(5) + t381 * t277 + t377 * t395;
t325 = qJD(5) + t328;
t284 = t325 * t308;
t207 = t284 + t227;
t622 = qJ(6) * t207;
t299 = t330 * t328;
t611 = -t299 + t395;
t621 = t378 * t611;
t620 = t382 * t611;
t619 = t516 * pkin(2) * t380;
t561 = t613 * t377;
t560 = t613 * t381;
t367 = t368 ^ 2;
t603 = sin(qJ(1));
t604 = cos(qJ(1));
t430 = g(1) * t603 - g(2) * t604;
t410 = qJDD(1) * pkin(1) + t430;
t403 = t376 * t410;
t431 = g(1) * t604 + g(2) * t603;
t606 = qJD(1) ^ 2;
t412 = pkin(1) * t606 + t431;
t446 = t488 * pkin(2);
t589 = pkin(10) * t373;
t390 = t367 * t589 + t380 * t412 + t384 * t403 + t446;
t547 = t373 * t380;
t601 = pkin(2) * t384;
t445 = -pkin(10) * t547 - t601;
t427 = t445 * t374;
t588 = pkin(10) * t375;
t486 = (-qJD(2) + t368) * t588;
t501 = pkin(9) + t588;
t590 = pkin(9) * t376;
t509 = t384 * t590;
t578 = t384 * g(3);
t387 = (-t578 - t501 * t512 + ((-t380 * t427 + t509) * qJD(1) + t384 * t486) * qJD(1)) * t374 + t390;
t436 = -pkin(1) + t445;
t491 = t384 * t516;
t579 = t376 * g(3);
t371 = t380 ^ 2;
t372 = t384 ^ 2;
t610 = -t371 - t372;
t389 = -t579 + (t436 * qJDD(1) + (-pkin(9) * t525 + (t375 * t525 * t610 - t373 * t491) * pkin(10) + t619) * qJD(1) - t430) * t374;
t618 = t373 * t389 + t375 * t387;
t306 = t310 ^ 2;
t324 = t325 ^ 2;
t262 = -t306 - t324;
t305 = t308 ^ 2;
t322 = pkin(3) * t341 - pkin(11) * t343;
t545 = t374 * t380;
t369 = g(3) * t545;
t296 = -t384 * t412 + t380 * t403 - t369 + pkin(10) * t442 - pkin(2) * t367 + (t501 * t511 + ((t380 * t590 + t384 * t427) * qJD(1) + t380 * t486) * qJD(1)) * t374;
t535 = t383 * t296;
t607 = t353 ^ 2;
t211 = -t607 * pkin(3) - pkin(11) * t401 - t341 * t322 + t379 * t618 + t535;
t333 = t353 * t341;
t292 = -t333 + t320;
t550 = t353 * t343;
t385 = -t373 * t387 + t375 * t389 - t292 * pkin(11) + (t397 + t550) * pkin(3);
t159 = t211 * t382 + t378 * t385;
t297 = pkin(4) * t328 - pkin(12) * t330;
t413 = t420 ^ 2;
t128 = -pkin(4) * t413 + pkin(12) * t395 - t328 * t297 + t159;
t228 = t379 * t296 - t383 * t618;
t210 = pkin(3) * t401 - pkin(11) * t607 + t343 * t322 + t228;
t316 = t420 * t328;
t257 = t277 - t316;
t411 = t420 * t330;
t150 = -t257 * pkin(12) + (-t276 + t411) * pkin(4) + t210;
t102 = t128 * t381 + t150 * t377;
t493 = t277 * t377 - t381 * t395;
t226 = -qJD(5) * t310 - t493;
t280 = pkin(5) * t325 - qJ(6) * t310;
t435 = qJ(6) * t226 - 0.2e1 * qJD(6) * t308 - t280 * t325 + t102;
t617 = -t435 + (t262 + t305) * pkin(5);
t551 = t343 * t341;
t396 = -t401 - t551;
t615 = t379 * t396;
t614 = t383 * t396;
t612 = -t284 + t227;
t517 = t374 * t606;
t399 = pkin(9) * t517 + t410;
t398 = t376 * t399;
t580 = t374 * pkin(9);
t400 = qJDD(1) * t580 - t412;
t334 = t380 * t400 + (t374 * g(3) - t398) * t384;
t335 = t380 * t398 + t384 * t400 - t369;
t609 = t380 * t334 + t384 * t335;
t204 = t310 * (qJD(5) - t325) + t493;
t370 = t374 ^ 2;
t608 = (t524 - t516) * t370;
t326 = t328 ^ 2;
t327 = t330 ^ 2;
t338 = t341 ^ 2;
t339 = t343 ^ 2;
t519 = qJD(6) * t310;
t302 = -0.2e1 * t519;
t101 = t128 * t377 - t150 * t381;
t421 = -t101 - t622 + t623;
t80 = t302 + t421;
t605 = pkin(5) * t80;
t156 = -t204 * t381 + t207 * t377;
t248 = -t305 - t306;
t123 = t156 * t378 - t248 * t382;
t600 = pkin(3) * t123;
t260 = -t324 - t305;
t182 = t260 * t381 - t561;
t203 = (qJD(5) + t325) * t310 + t493;
t135 = t182 * t378 - t203 * t382;
t599 = pkin(3) * t135;
t220 = t270 + t275;
t562 = t220 * t381;
t188 = -t262 * t377 - t562;
t140 = t188 * t378 - t382 * t612;
t598 = pkin(3) * t140;
t597 = pkin(3) * t379;
t596 = pkin(3) * t383;
t181 = t260 * t377 + t560;
t595 = pkin(4) * t181;
t563 = t220 * t377;
t187 = t262 * t381 - t563;
t594 = pkin(4) * t187;
t593 = pkin(4) * t378;
t592 = pkin(4) * t382;
t591 = pkin(5) * t207;
t586 = pkin(11) * t123;
t585 = pkin(11) * t135;
t584 = pkin(11) * t140;
t154 = -t204 * t377 - t207 * t381;
t583 = pkin(12) * t154;
t582 = pkin(12) * t181;
t581 = pkin(12) * t187;
t124 = t156 * t382 + t248 * t378;
t470 = t124 * t379 - t154 * t383;
t74 = t123 * t375 + t373 * t470;
t75 = -t123 * t373 + t375 * t470;
t99 = t124 * t383 + t154 * t379;
t577 = pkin(1) * (-t374 * t74 + (t380 * t99 + t384 * t75) * t376) + (-t380 * t75 + t384 * t99) * t580;
t136 = t182 * t382 + t203 * t378;
t113 = t136 * t383 + t181 * t379;
t468 = t136 * t379 - t181 * t383;
t85 = t135 * t375 + t373 * t468;
t86 = -t135 * t373 + t375 * t468;
t576 = pkin(1) * (-t374 * t85 + (t113 * t380 + t384 * t86) * t376) + (t113 * t384 - t380 * t86) * t580;
t141 = t188 * t382 + t378 * t612;
t117 = t141 * t383 + t187 * t379;
t467 = t141 * t379 - t187 * t383;
t88 = t140 * t375 + t373 * t467;
t89 = -t140 * t373 + t375 * t467;
t575 = pkin(1) * (-t374 * t88 + (t117 * t380 + t384 * t89) * t376) + (t117 * t384 - t380 * t89) * t580;
t574 = pkin(2) * t75 + t589 * t99;
t573 = t377 * t80;
t158 = t211 * t378 - t382 * t385;
t127 = -pkin(4) * t395 - pkin(12) * t413 + t297 * t330 + t158;
t106 = -pkin(5) * t226 - qJ(6) * t305 + t280 * t310 + qJDD(6) + t127;
t81 = -pkin(5) * t305 + t435;
t47 = t381 * t81 - t573;
t35 = -t106 * t382 + t378 * t47;
t36 = t106 * t378 + t382 * t47;
t570 = t381 * t80;
t46 = t377 * t81 + t570;
t476 = t36 * t379 - t383 * t46;
t11 = -t35 * t373 + t375 * t476;
t572 = t380 * t11;
t62 = t101 * t377 + t102 * t381;
t51 = t127 * t378 + t382 * t62;
t61 = -t101 * t381 + t102 * t377;
t473 = t379 * t51 - t383 * t61;
t50 = -t127 * t382 + t378 * t62;
t16 = -t373 * t50 + t375 * t473;
t571 = t380 * t16;
t569 = pkin(2) * t86 + t113 * t589;
t568 = pkin(2) * t89 + t117 * t589;
t567 = t127 * t377;
t566 = t127 * t381;
t565 = t210 * t378;
t564 = t210 * t382;
t507 = pkin(2) * t511;
t261 = t375 * t579 + (-t375 * (-t410 - t507) + ((-t445 * t547 - t375 * (t588 * t610 - pkin(9))) * t525 - t375 * t619) * qJD(1)) * t374 + (t390 + (-pkin(9) * t512 - t578 + (0.2e1 * pkin(10) * t368 * t542 + t509 * qJD(1)) * qJD(1)) * t374) * t373;
t559 = t261 * t379;
t558 = t261 * t383;
t266 = t299 + t395;
t557 = t266 * t378;
t556 = t266 * t382;
t311 = t401 - t551;
t555 = t311 * t379;
t554 = t311 * t383;
t553 = t325 * t377;
t552 = t325 * t381;
t549 = t353 * t383;
t544 = t374 * t384;
t543 = t375 * t379;
t254 = -t330 * t341 + t492;
t258 = t277 + t316;
t190 = -t254 * t378 - t258 * t382;
t192 = -t254 * t382 + t258 * t378;
t278 = t326 + t327;
t461 = t192 * t379 + t278 * t383;
t130 = -t190 * t373 + t375 * t461;
t541 = t380 * t130;
t286 = -t413 - t326;
t230 = t286 * t378 + t620;
t231 = t286 * t382 - t621;
t419 = 0.2e1 * qJD(4) + t341;
t255 = -t330 * t419 - t492;
t458 = t231 * t379 + t255 * t383;
t161 = -t230 * t373 + t375 * t458;
t540 = t380 * t161;
t295 = -t327 - t413;
t233 = t295 * t382 - t557;
t234 = -t295 * t378 - t556;
t259 = t328 * t419 + t448;
t457 = t234 * t379 + t259 * t383;
t165 = -t233 * t373 + t375 * t457;
t539 = t380 * t165;
t307 = -t338 - t339;
t293 = t333 + t320;
t393 = (-qJD(3) + t353) * t343 - t402;
t453 = -t293 * t383 + t379 * t393;
t218 = -t307 * t373 + t375 * t453;
t538 = t380 * t218;
t518 = t370 * t606;
t366 = t384 * t380 * t518;
t355 = t366 + t488;
t536 = t380 * t355;
t356 = -t366 + t488;
t533 = t384 * t356;
t532 = -pkin(3) * t154 + pkin(11) * t124;
t531 = -pkin(3) * t181 + pkin(11) * t136;
t530 = -pkin(3) * t187 + pkin(11) * t141;
t529 = -pkin(4) * t248 + pkin(12) * t156;
t528 = pkin(4) * t203 - pkin(12) * t182;
t527 = pkin(4) * t612 - pkin(12) * t188;
t515 = qJD(3) + t353;
t506 = t378 * t270;
t505 = t382 * t270;
t504 = t379 * t299;
t503 = t383 * t299;
t502 = -pkin(2) * t74 + t588 * t99;
t498 = t371 * t518;
t497 = t372 * t518;
t496 = -pkin(2) * t85 + t113 * t588;
t495 = -pkin(2) * t88 + t117 * t588;
t105 = t158 * t378 + t159 * t382;
t487 = t374 * t494;
t485 = -pkin(4) * t127 + pkin(12) * t62;
t484 = -pkin(11) * t379 - t596;
t483 = -t528 - t566;
t482 = -t527 + t567;
t319 = -t607 - t338;
t273 = t319 * t383 - t615;
t481 = pkin(10) * t273 + t558;
t321 = -t339 - t607;
t279 = -t321 * t379 + t554;
t480 = pkin(10) * t279 - t559;
t118 = -pkin(4) * t154 + t591;
t68 = -qJ(6) * t204 + (-t248 - t305) * pkin(5) + t435;
t303 = 0.2e1 * t519;
t70 = t303 - t421 + t622;
t33 = -t377 * t68 + t381 * t70 - t583;
t23 = -t118 * t378 + t33 * t382 - t586;
t429 = t377 * t70 + t381 * t68 + t529;
t25 = -t429 - t600;
t479 = t23 * t379 + t25 * t383;
t409 = t421 + t623;
t67 = t303 - t409 - t595;
t90 = -pkin(5) * t203 + qJ(6) * t260 - t106;
t71 = -qJ(6) * t560 - t377 * t90 - t582;
t31 = -t378 * t67 + t382 * t71 - t585;
t417 = -qJ(6) * t561 + t381 * t90 - t528;
t54 = -t417 - t599;
t478 = t31 * t379 + t383 * t54;
t69 = -t594 - t617;
t103 = -qJ(6) * t262 + t106;
t166 = -pkin(5) * t612 - qJ(6) * t220;
t72 = t103 * t381 - t166 * t377 - t581;
t34 = -t378 * t69 + t382 * t72 - t584;
t426 = t103 * t377 + t166 * t381 - t527;
t55 = -t426 - t598;
t477 = t34 * t379 + t383 * t55;
t53 = -t61 - t583;
t37 = t154 * t593 + t382 * t53 - t586;
t444 = t62 + t529;
t40 = -t444 - t600;
t475 = t37 * t379 + t383 * t40;
t109 = t567 - t582;
t82 = t101 - t595;
t49 = t109 * t382 - t378 * t82 - t585;
t76 = -t483 - t599;
t474 = t379 * t49 + t383 * t76;
t110 = t566 - t581;
t83 = t102 - t594;
t52 = t110 * t382 - t378 * t83 - t584;
t77 = -t482 - t598;
t472 = t379 * t52 + t383 * t77;
t471 = t105 * t379 - t210 * t383;
t155 = -t203 * t381 - t377 * t612;
t268 = t306 - t305;
t132 = t155 * t382 + t268 * t378;
t153 = -t203 * t377 + t381 * t612;
t469 = t132 * t379 - t153 * t383;
t282 = -t306 + t324;
t195 = -t282 * t377 + t560;
t146 = t195 * t382 + t207 * t378;
t193 = t282 * t381 + t561;
t466 = t146 * t379 - t193 * t383;
t281 = t305 - t324;
t196 = t281 * t381 - t563;
t147 = t196 * t382 - t204 * t378;
t194 = t281 * t377 + t562;
t465 = t147 * t379 - t194 * t383;
t104 = -t158 * t382 + t159 * t378;
t198 = -t226 * t377 + t308 * t552;
t174 = t198 * t382 - t506;
t197 = -t226 * t381 - t308 * t553;
t464 = t174 * t379 + t197 * t383;
t200 = t227 * t381 - t310 * t553;
t175 = t200 * t382 + t506;
t199 = t227 * t377 + t310 * t552;
t463 = t175 * t379 - t199 * t383;
t191 = t255 * t382 - t257 * t378;
t298 = t327 - t326;
t462 = t191 * t379 - t298 * t383;
t238 = (-t308 * t381 + t310 * t377) * t325;
t215 = t238 * t382 + t275 * t378;
t237 = (-t308 * t377 - t310 * t381) * t325;
t460 = t215 * t379 - t237 * t383;
t349 = t374 * t399 + t579;
t425 = t368 * t373 + t489;
t229 = t535 + (t375 * (-g(3) * t544 + t384 * t398 + t446) + t373 * (-t374 * t507 - t349) + (t375 * (t368 * t425 - t375 * t487) + t373 * (-t373 * t487 - t425 * t499)) * pkin(10) + (t375 * t412 + ((-pkin(10) * t373 ^ 2 - t375 * t501) * qJDD(1) + (t373 * t516 + t489) * qJD(1) * pkin(2)) * t374) * t380) * t379;
t459 = t228 * t383 - t229 * t379;
t171 = t228 * t379 + t229 * t383;
t315 = -t327 + t413;
t243 = -t315 * t378 + t620;
t456 = t243 * t379 - t258 * t383;
t314 = t326 - t413;
t244 = t314 * t382 - t557;
t455 = t244 * t379 + t254 * t383;
t288 = t343 * t515 + t402;
t454 = -t288 * t379 + t292 * t383;
t452 = t321 * t383 + t555;
t331 = t338 - t607;
t451 = t331 * t379 - t554;
t332 = -t339 + t607;
t450 = t332 * t383 + t615;
t449 = t319 * t379 + t614;
t447 = -t341 * t379 - t343 * t383;
t441 = (-t373 * t74 - t375 * t75) * pkin(10);
t440 = (-t373 * t85 - t375 * t86) * pkin(10);
t439 = (-t373 * t88 - t375 * t89) * pkin(10);
t405 = t382 * t316;
t251 = -t276 * t378 + t405;
t438 = t251 * t379 + t503;
t406 = t378 * t411;
t253 = t277 * t382 - t406;
t437 = t253 * t379 - t503;
t18 = t36 * t383 + t379 * t46;
t58 = -pkin(5) * t106 + qJ(6) * t81;
t19 = -pkin(12) * t46 - qJ(6) * t570 - t377 * t58;
t27 = -pkin(4) * t46 - t605;
t5 = -pkin(11) * t35 + t19 * t382 - t27 * t378;
t408 = -pkin(4) * t106 + pkin(12) * t47 - qJ(6) * t573 + t381 * t58;
t9 = -pkin(3) * t35 - t408;
t432 = pkin(10) * t18 + t379 * t5 + t383 * t9;
t17 = -pkin(11) * t50 + (-pkin(12) * t382 + t593) * t61;
t22 = -pkin(3) * t50 - t485;
t24 = t379 * t61 + t383 * t51;
t428 = pkin(10) * t24 + t17 * t379 + t22 * t383;
t120 = -pkin(3) * t230 + t158;
t169 = -pkin(11) * t230 + t565;
t183 = t231 * t383 - t255 * t379;
t424 = pkin(10) * t183 + t120 * t383 + t169 * t379;
t121 = -pkin(3) * t233 + t159;
t170 = -pkin(11) * t233 + t564;
t184 = t234 * t383 - t259 * t379;
t423 = pkin(10) * t184 + t121 * t383 + t170 * t379;
t249 = t293 * t379 + t383 * t393;
t422 = pkin(10) * t249 + t171;
t176 = t192 * t383 - t278 * t379;
t94 = -pkin(11) * t190 - t104;
t418 = pkin(10) * t176 - t190 * t596 + t379 * t94;
t407 = t378 * t316;
t404 = t382 * t411;
t394 = t383 * t397;
t272 = -t405 + t406;
t392 = t379 * t272 - t383 * t395;
t360 = t368 * t499;
t359 = t368 * t500;
t358 = (t371 - t372) * t518;
t357 = -t367 - t497;
t354 = -t498 - t367;
t348 = -t359 + t415;
t347 = t359 + t415;
t346 = -t360 + t416;
t323 = t339 - t338;
t291 = -t341 * t515 + t391;
t287 = t353 * t373 * t447 - t375 * t401;
t271 = -t407 - t404;
t264 = t320 * t548 + (t341 * t375 + t353 * t546) * t343;
t263 = -t373 * t394 + (-t343 * t375 + t353 * t548) * t341;
t252 = t277 * t378 + t404;
t250 = t276 * t382 + t407;
t246 = t373 * t451 + t375 * t393;
t245 = t293 * t375 + t373 * t450;
t242 = t314 * t378 + t556;
t241 = t315 * t382 + t621;
t240 = -t291 * t373 + t375 * t452;
t239 = t291 * t375 + t373 * t452;
t236 = -t288 * t373 + t375 * t449;
t235 = t288 * t375 + t373 * t449;
t232 = t323 * t375 + t373 * t454;
t216 = t375 * t271 + t373 * t392;
t214 = t238 * t378 - t275 * t382;
t189 = t255 * t378 + t257 * t382;
t178 = t252 * t375 + t373 * t437;
t177 = t250 * t375 + t373 * t438;
t173 = t200 * t378 - t505;
t172 = t198 * t378 + t505;
t168 = t242 * t375 + t373 * t455;
t167 = t241 * t375 + t373 * t456;
t164 = t233 * t375 + t373 * t457;
t163 = t261 * t373 - t375 * t459;
t162 = -t261 * t375 - t373 * t459;
t160 = t230 * t375 + t373 * t458;
t145 = t196 * t378 + t204 * t382;
t144 = t195 * t378 - t207 * t382;
t143 = pkin(2) * t240 - t229 * t375 + t373 * t480;
t142 = pkin(3) * t259 + pkin(11) * t234 + t565;
t138 = pkin(2) * t236 - t228 * t375 + t373 * t481;
t137 = pkin(3) * t255 + pkin(11) * t231 - t564;
t133 = t189 * t375 + t373 * t462;
t131 = t155 * t378 - t268 * t382;
t129 = t190 * t375 + t373 * t461;
t125 = t214 * t375 + t373 * t460;
t119 = pkin(2) * t218 + t373 * t422;
t116 = pkin(2) * t163 + t171 * t589;
t108 = t173 * t375 + t373 * t463;
t107 = t172 * t375 + t373 * t464;
t96 = -pkin(3) * t210 + pkin(11) * t105;
t95 = t105 * t383 + t210 * t379;
t93 = t145 * t375 + t373 * t465;
t92 = t144 * t375 + t373 * t466;
t91 = pkin(3) * t278 + pkin(11) * t192 + t105;
t79 = t376 * t125 + (t380 * (t215 * t383 + t237 * t379) + t384 * (-t214 * t373 + t375 * t460)) * t374;
t78 = t131 * t375 + t373 * t469;
t65 = -t104 * t373 + t375 * t471;
t64 = t104 * t375 + t373 * t471;
t60 = pkin(2) * t165 + t142 * t375 + t373 * t423;
t59 = pkin(2) * t161 + t137 * t375 + t373 * t424;
t57 = t376 * t108 + (t380 * (t175 * t383 + t199 * t379) + t384 * (-t173 * t373 + t375 * t463)) * t374;
t56 = t376 * t107 + (t380 * (t174 * t383 - t197 * t379) + t384 * (-t172 * t373 + t375 * t464)) * t374;
t45 = t110 * t378 + t382 * t83 + t530;
t44 = t109 * t378 + t382 * t82 + t531;
t43 = t376 * t93 + (t380 * (t147 * t383 + t194 * t379) + t384 * (-t145 * t373 + t375 * t465)) * t374;
t42 = t376 * t92 + (t380 * (t146 * t383 + t193 * t379) + t384 * (-t144 * t373 + t375 * t466)) * t374;
t38 = pkin(2) * t130 + t373 * t418 + t375 * t91;
t32 = -t154 * t592 + t378 * t53 + t532;
t30 = t376 * t78 + (t380 * (t132 * t383 + t153 * t379) + t384 * (-t131 * t373 + t375 * t469)) * t374;
t29 = t378 * t72 + t382 * t69 + t530;
t26 = t378 * t71 + t382 * t67 + t531;
t21 = t118 * t382 + t33 * t378 + t532;
t20 = pkin(2) * t65 + t375 * t96 + (pkin(10) * t95 + t104 * t484) * t373;
t15 = t373 * t473 + t375 * t50;
t14 = pkin(11) * t51 + (-pkin(12) * t378 - pkin(3) - t592) * t61;
t13 = t373 * t472 + t375 * t45 + t568;
t12 = t373 * t474 + t375 * t44 + t569;
t10 = t35 * t375 + t373 * t476;
t8 = t29 * t375 + t373 * t477 + t568;
t7 = t26 * t375 + t373 * t478 + t569;
t6 = t32 * t375 + t373 * t475 + t574;
t4 = -pkin(3) * t46 + pkin(11) * t36 + t19 * t378 + t27 * t382;
t3 = t21 * t375 + t373 * t479 + t574;
t2 = pkin(2) * t16 + t14 * t375 + t373 * t428;
t1 = pkin(2) * t11 + t373 * t432 + t375 * t4;
t28 = [0, 0, 0, 0, 0, qJDD(1), t430, t431, 0, 0, (t370 * t512 - t522 * t608) * t380, t376 * t358 + (t380 * t348 + (t360 + t416) * t384) * t374, t376 * t346 + (t536 + t384 * (t367 - t498)) * t374, (t370 * t511 + t523 * t608) * t384, t376 * t347 + (t380 * (-t367 + t497) + t533) * t374, t376 * t488, (-t334 + pkin(1) * (t355 * t384 + t357 * t380)) * t376 + (t384 * t349 + pkin(1) * t348 + pkin(9) * (t357 * t384 - t536)) * t374, -t349 * t545 - t376 * t335 + pkin(1) * ((t354 * t384 - t356 * t380) * t376 + (-qJD(1) * t491 - t512) * t370) + (-t354 * t380 - t533) * t580, pkin(1) * ((-t346 * t384 + t347 * t380) * t376 - t610 * t370 * t517) + (t346 * t380 + t347 * t384) * t580 + t609 * t374, pkin(1) * (t374 * t349 + (-t334 * t384 + t335 * t380) * t376) + t609 * t580, t376 * t264 + (t380 * (t320 * t383 - t379 * t550) + t384 * (t320 * t543 + (-t341 * t373 + t375 * t549) * t343)) * t374, t376 * t232 + (t380 * (-t288 * t383 - t292 * t379) + t384 * (-t323 * t373 + t375 * t454)) * t374, t376 * t245 + (t380 * (-t332 * t379 + t614) + t384 * (-t293 * t373 + t375 * t450)) * t374, (t341 * t549 + t379 * t397) * t545 + (-t375 * t394 + (t343 * t373 + t353 * t543) * t341) * t544 + t376 * t263, t376 * t246 + (t380 * (t331 * t383 + t555) + t384 * (-t373 * t393 + t375 * t451)) * t374, t401 * t373 * t544 + t376 * t287 + (t380 * (-t341 * t383 + t343 * t379) + t447 * t542) * t374 * t353, (t138 + pkin(1) * (t236 * t384 + t273 * t380)) * t376 + (t380 * (-t559 + (-t235 * t373 - t236 * t375) * pkin(10)) + t384 * (-pkin(2) * t235 + t228 * t373 + t375 * t481) - pkin(1) * t235 + pkin(9) * (-t236 * t380 + t273 * t384)) * t374, (t143 + pkin(1) * (t240 * t384 + t279 * t380)) * t376 + (t380 * (-t558 + (-t239 * t373 - t240 * t375) * pkin(10)) + t384 * (-pkin(2) * t239 + t229 * t373 + t375 * t480) - pkin(1) * t239 + pkin(9) * (-t240 * t380 + t279 * t384)) * t374, (t119 + pkin(1) * (t218 * t384 + t249 * t380)) * t376 + (t380 * t459 + pkin(9) * (t249 * t384 - t538) + t436 * (t307 * t375 + t373 * t453) + (-pkin(10) * t538 + t384 * t422) * t375) * t374, (t116 + pkin(1) * (t163 * t384 + t171 * t380)) * t376 + (pkin(9) * (-t163 * t380 + t171 * t384) + (-pkin(1) - t601) * t162 + (t380 * (-t162 * t373 - t163 * t375) + t171 * t542) * pkin(10)) * t374, t376 * t178 + (t380 * (t253 * t383 + t504) + t384 * (-t252 * t373 + t375 * t437)) * t374, t376 * t133 + (t380 * (t191 * t383 + t298 * t379) + t384 * (-t189 * t373 + t375 * t462)) * t374, t376 * t167 + (t380 * (t243 * t383 + t258 * t379) + t384 * (-t241 * t373 + t375 * t456)) * t374, t376 * t177 + (t380 * (t251 * t383 - t504) + t384 * (-t250 * t373 + t375 * t438)) * t374, t376 * t168 + (t380 * (t244 * t383 - t254 * t379) + t384 * (-t242 * t373 + t375 * t455)) * t374, (t383 * t272 + t379 * t395) * t545 + (-t271 * t373 + t375 * t392) * t544 + t376 * t216, (t59 + pkin(1) * (t161 * t384 + t183 * t380)) * t376 + (t380 * (-t120 * t379 - t160 * t589 + t169 * t383) + t384 * (-pkin(2) * t160 - t137 * t373) - pkin(1) * t160 + pkin(9) * (t183 * t384 - t540) + (-pkin(10) * t540 + t384 * t424) * t375) * t374, (t60 + pkin(1) * (t165 * t384 + t184 * t380)) * t376 + (t380 * (-t121 * t379 - t164 * t589 + t170 * t383) + t384 * (-pkin(2) * t164 - t142 * t373) - pkin(1) * t164 + pkin(9) * (t184 * t384 - t539) + (-pkin(10) * t539 + t384 * t423) * t375) * t374, (t38 + pkin(1) * (t130 * t384 + t176 * t380)) * t376 + (t380 * (-t129 * t589 + t190 * t597 + t383 * t94) + t384 * (-pkin(2) * t129 - t373 * t91) - pkin(1) * t129 + pkin(9) * (t176 * t384 - t541) + (-pkin(10) * t541 + t384 * t418) * t375) * t374, (t20 + pkin(1) * (t380 * t95 + t384 * t65)) * t376 + (t384 * (-pkin(2) * t64 - t373 * t96) - pkin(1) * t64 + pkin(9) * (-t380 * t65 + t384 * t95) + (t380 * (-t373 * t64 - t375 * t65) + t95 * t542) * pkin(10) + (t380 * (-pkin(11) * t383 + t597) + t484 * t542) * t104) * t374, t57, t30, t42, t56, t43, t79, t376 * t12 + (t380 * (-t379 * t76 + t383 * t49 + t440) + t384 * (-t373 * t44 + t375 * t474 + t496)) * t374 + t576, t376 * t13 + (t380 * (-t379 * t77 + t383 * t52 + t439) + t384 * (-t373 * t45 + t375 * t472 + t495)) * t374 + t575, t376 * t6 + (t380 * (t37 * t383 - t379 * t40 + t441) + t384 * (-t32 * t373 + t375 * t475 + t502)) * t374 + t577, (t2 + pkin(1) * (t16 * t384 + t24 * t380)) * t376 + (t380 * (-t15 * t589 + t17 * t383 - t22 * t379) + t384 * (-pkin(2) * t15 - t14 * t373) - pkin(1) * t15 + pkin(9) * (t24 * t384 - t571) + (-pkin(10) * t571 + t384 * t428) * t375) * t374, t57, t30, t42, t56, t43, t79, t376 * t7 + (t380 * (t31 * t383 - t379 * t54 + t440) + t384 * (-t26 * t373 + t375 * t478 + t496)) * t374 + t576, t376 * t8 + (t380 * (t34 * t383 - t379 * t55 + t439) + t384 * (-t29 * t373 + t375 * t477 + t495)) * t374 + t575, t376 * t3 + (t380 * (t23 * t383 - t25 * t379 + t441) + t384 * (-t21 * t373 + t375 * t479 + t502)) * t374 + t577, (t1 + pkin(1) * (t11 * t384 + t18 * t380)) * t376 + (t380 * (-t10 * t589 - t379 * t9 + t383 * t5) + t384 * (-pkin(2) * t10 - t373 * t4) - pkin(1) * t10 + pkin(9) * (t18 * t384 - t572) + (-pkin(10) * t572 + t384 * t432) * t375) * t374; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t366, t358, t346, t366, t347, t488, -t334, -t335, 0, 0, t264, t232, t245, t263, t246, t287, t138, t143, t119, t116, t178, t133, t167, t177, t168, t216, t59, t60, t38, t20, t108, t78, t92, t107, t93, t125, t12, t13, t6, t2, t108, t78, t92, t107, t93, t125, t7, t8, t3, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t551, t323, t293, -t551, t393, -t401, -t228, -t229, 0, 0, t252, t189, t241, t250, t242, t271, t137, t142, t91, t96, t173, t131, t144, t172, t145, t214, t44, t45, t32, t14, t173, t131, t144, t172, t145, t214, t26, t29, t21, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, t298, t258, -t299, -t254, t395, -t158, -t159, 0, 0, t199, t153, t193, -t197, t194, t237, t483, t482, t444, t485, t199, t153, t193, -t197, t194, t237, t417, t426, t429, t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270, t268, t207, -t270, -t204, t275, -t101, -t102, 0, 0, t270, t268, t207, -t270, -t204, t275, t302 + t409, t617, -t591, t605; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, t612, t248, t106;];
tauJ_reg  = t28;
