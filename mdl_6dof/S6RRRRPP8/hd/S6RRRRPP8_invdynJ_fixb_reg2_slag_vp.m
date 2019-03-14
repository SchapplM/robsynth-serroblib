% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPP8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:38:37
% EndTime: 2019-03-09 21:39:10
% DurationCPUTime: 16.76s
% Computational Cost: add. (17789->931), mult. (43884->1128), div. (0->0), fcn. (34688->10), ass. (0->396)
t344 = sin(qJ(3));
t556 = cos(pkin(6));
t467 = t556 * qJD(1);
t417 = t467 + qJD(2);
t342 = sin(pkin(6));
t345 = sin(qJ(2));
t522 = qJD(1) * t345;
t480 = t342 * t522;
t584 = cos(qJ(3));
t226 = t344 * t480 - t584 * t417;
t217 = qJD(4) + t226;
t398 = t344 * t417;
t459 = t556 * qJDD(1);
t412 = t459 + qJDD(2);
t476 = qJD(3) * t584;
t347 = cos(qJ(2));
t520 = qJD(2) * t347;
t478 = t344 * t520;
t512 = qJDD(1) * t345;
t119 = qJD(3) * t398 + t342 * (qJD(1) * (t345 * t476 + t478) + t344 * t512) - t584 * t412;
t113 = qJDD(4) + t119;
t483 = t347 * t584;
t446 = qJD(2) * t483;
t485 = t345 * t584;
t620 = t226 * qJD(3) - t344 * t412 - (qJD(1) * t446 + qJDD(1) * t485) * t342;
t623 = t344 * t620;
t523 = qJD(1) * t342;
t324 = t347 * t523;
t449 = t344 * t324;
t519 = qJD(3) * t344;
t603 = t449 - t519;
t448 = pkin(1) * t467;
t261 = -pkin(8) * t480 + t347 * t448;
t437 = pkin(2) * t345 - pkin(9) * t347;
t262 = t437 * t523;
t156 = t584 * t261 + t344 * t262;
t133 = pkin(10) * t480 + t156;
t493 = pkin(1) * t556;
t332 = t345 * t493;
t411 = pkin(3) * t344 - pkin(10) * t584;
t538 = t342 * t347;
t165 = (t332 + (pkin(8) + t411) * t538) * qJD(1);
t301 = -pkin(3) * t584 - t344 * pkin(10) - pkin(2);
t343 = sin(qJ(4));
t346 = cos(qJ(4));
t475 = t584 * qJD(4);
t457 = pkin(9) * t475;
t518 = qJD(4) * t343;
t622 = -t343 * t133 + t301 * t518 + (t165 + t457) * t346;
t621 = t113 * qJ(5) + t217 * qJD(5);
t539 = t342 * t345;
t283 = -pkin(8) * t539 + t347 * t493;
t265 = qJD(2) * t283;
t585 = cos(qJ(1));
t440 = t556 * t585;
t583 = sin(qJ(1));
t279 = t345 * t440 + t347 * t583;
t489 = t342 * t585;
t193 = t279 * t584 - t344 * t489;
t278 = t345 * t583 - t347 * t440;
t134 = t193 * t343 - t278 * t346;
t135 = t193 * t346 + t278 * t343;
t488 = t342 * t584;
t447 = qJD(1) * t488;
t228 = t345 * t447 + t398;
t425 = t324 - qJD(3);
t168 = t228 * t343 + t346 * t425;
t461 = t217 * t343;
t533 = t346 * t113;
t619 = -t168 * t228 + t217 * t461 - t533;
t170 = t346 * t228 - t343 * t425;
t419 = t347 * t447;
t212 = t343 * t419 - t346 * t480;
t230 = (t343 * t345 + t346 * t483) * t342;
t213 = qJD(1) * t230;
t484 = t346 * t584;
t486 = t343 * t584;
t513 = qJD(1) * qJD(2);
t472 = t345 * t513;
t443 = t342 * t472;
t511 = qJDD(1) * t347;
t323 = t342 * t511;
t509 = qJDD(3) - t323;
t392 = t443 + t509;
t352 = -t343 * t620 - t346 * t392;
t68 = qJD(4) * t170 + t352;
t567 = t68 * t346;
t458 = qJD(4) * t425;
t67 = t228 * t518 - t343 * t392 + (t620 + t458) * t346;
t568 = t67 * t343;
t618 = t344 * (qJD(4) * (t168 * t343 - t170 * t346) - t567 + t568) - (t168 * t484 + t170 * t486) * qJD(3) + t213 * t168 + t170 * t212;
t410 = t343 * t476 - t212;
t517 = qJD(4) * t346;
t536 = t343 * t113;
t356 = -t410 * t217 + t344 * (t168 * t425 - t217 * t517 - t536) + t584 * t68;
t409 = t346 * t476 - t213;
t617 = -t409 * t217 + t344 * (t170 * t425 + t217 * t518 - t533) - t67 * t584;
t550 = t170 * t217;
t552 = t168 * t217;
t613 = t67 + t552;
t12 = t343 * (t68 + t550) + t346 * t613;
t339 = t342 ^ 2;
t616 = 0.2e1 * t339;
t615 = t68 * qJ(6) + t168 * qJD(6);
t291 = t411 * qJD(3);
t614 = t346 * t291 - t622;
t201 = -pkin(2) * t417 - t261;
t100 = t226 * pkin(3) - t228 * pkin(10) + t201;
t526 = pkin(8) * t538 + t332;
t252 = t556 * pkin(9) + t526;
t202 = qJD(2) * pkin(9) + qJD(1) * t252;
t415 = -pkin(2) * t347 - pkin(9) * t345 - pkin(1);
t214 = t415 * t523;
t118 = t202 * t584 + t344 * t214;
t103 = -pkin(10) * t425 + t118;
t51 = t346 * t100 - t343 * t103;
t531 = qJD(5) - t51;
t611 = 0.2e1 * t621;
t610 = qJD(5) * t343 + t118;
t506 = pkin(9) * t519;
t609 = t156 + t506;
t276 = t344 * t539 - t556 * t584;
t189 = -qJD(3) * t276 + t342 * t446;
t608 = -qJD(4) * t538 + t189;
t215 = t217 ^ 2;
t587 = t170 ^ 2;
t607 = -t215 - t587;
t337 = t584 * pkin(4);
t606 = -t584 * pkin(5) - t337;
t605 = (qJDD(2) + 0.2e1 * t459) * t342;
t604 = t419 - t476;
t588 = t168 ^ 2;
t602 = t587 - t588;
t266 = t526 * qJD(2);
t110 = t113 * pkin(4);
t601 = t110 - qJDD(5);
t86 = t346 * t133 + t343 * t165;
t600 = -qJ(5) * t603 - t584 * qJD(5) - t86;
t117 = -t344 * t202 + t584 * t214;
t102 = pkin(3) * t425 - t117;
t383 = -t170 * qJ(5) + t102;
t50 = t168 * pkin(4) + t383;
t580 = pkin(10) * t113;
t599 = t217 * t50 - t580;
t586 = pkin(4) + pkin(5);
t31 = -t168 * t586 + qJD(6) - t383;
t439 = t556 * t583;
t281 = -t345 * t439 + t347 * t585;
t487 = t342 * t583;
t197 = t281 * t584 + t344 * t487;
t280 = t345 * t585 + t347 * t439;
t138 = t197 * t343 - t280 * t346;
t277 = t342 * t485 + t344 * t556;
t190 = t277 * t343 + t346 * t538;
t413 = qJD(2) * t448;
t442 = pkin(1) * t459;
t496 = pkin(8) * t323 + t345 * t442 + t347 * t413;
t172 = -pkin(8) * t443 + t496;
t146 = pkin(9) * t412 + t172;
t406 = t437 * qJD(2);
t153 = (qJD(1) * t406 + qJDD(1) * t415) * t342;
t399 = -t584 * t146 - t344 * t153 + t202 * t519 - t214 * t476;
t39 = pkin(10) * t392 - t399;
t470 = t342 * t512;
t471 = t347 * t513;
t454 = t345 * t413 - t347 * t442 + (t342 * t471 + t470) * pkin(8);
t147 = -pkin(2) * t412 + t454;
t45 = t119 * pkin(3) + pkin(10) * t620 + t147;
t468 = t100 * t518 + t103 * t517 + t343 * t39 - t346 * t45;
t372 = g(1) * t138 + g(2) * t134 + g(3) * t190 - t468;
t364 = t372 + t601;
t572 = qJ(6) * t67;
t597 = (qJD(6) + t31) * t170 + t364 - t572;
t554 = qJ(5) * t346;
t593 = t343 * t586 - t554;
t537 = t343 * qJ(5);
t592 = -t586 * t346 - t537;
t349 = qJD(1) ^ 2;
t582 = pkin(5) * t113;
t581 = pkin(5) * t190;
t579 = pkin(10) * t170;
t196 = t281 * t344 - t487 * t584;
t578 = pkin(10) * t196;
t463 = -t279 * t344 - t585 * t488;
t577 = t463 * pkin(10);
t576 = pkin(10) - qJ(6);
t445 = qJ(6) * t476;
t494 = -pkin(9) * t343 - pkin(4);
t514 = t346 * qJD(6);
t575 = qJ(6) * t213 + t586 * t449 + (-t445 - t291) * t346 + (qJ(6) * t518 - t514 + (-pkin(5) + t494) * qJD(3)) * t344 + t622;
t529 = t343 * t291 + t301 * t517;
t534 = t344 * t346;
t574 = -qJ(6) * t212 + (-pkin(9) * qJD(3) + qJ(6) * qJD(4)) * t534 + (qJD(6) * t344 + t445 - t457) * t343 + t529 + t600;
t573 = qJ(5) * t68;
t208 = t217 * qJ(5);
t52 = t100 * t343 + t103 * t346;
t33 = qJ(6) * t168 + t52;
t29 = t208 + t33;
t571 = t217 * t29;
t42 = t208 + t52;
t570 = t217 * t42;
t569 = t217 * t52;
t155 = -t344 * t261 + t584 * t262;
t401 = pkin(3) * t480 + t155;
t379 = qJ(5) * t213 + t401;
t481 = t584 * qJ(5);
t507 = t584 * pkin(9);
t386 = -t346 * t481 + t507;
t515 = qJD(5) * t346;
t566 = t586 * t212 + (qJD(4) * t592 + t515) * t344 + (t343 * t606 - t386) * qJD(3) - t379;
t148 = (-t343 * t475 - t346 * t519) * pkin(9) + t529;
t565 = t148 + t600;
t564 = pkin(4) * t449 + t494 * t519 - t614;
t427 = pkin(4) * t346 + t537;
t563 = -pkin(4) * t212 + (qJD(4) * t427 - t515) * t344 + (pkin(4) * t486 + t386) * qJD(3) + t379;
t562 = t148 - t86;
t561 = t343 * t506 + t614;
t560 = -t217 * t593 + t610;
t553 = qJ(6) * t226;
t150 = pkin(3) * t228 + pkin(10) * t226;
t76 = t346 * t117 + t343 * t150;
t56 = t228 * qJ(5) + t76;
t559 = t343 * t553 - t518 * t576 - t514 - t56;
t426 = pkin(4) * t343 - t554;
t558 = t217 * t426 - t610;
t106 = t343 * t117;
t309 = t576 * t346;
t557 = qJD(4) * t309 - qJD(6) * t343 - t106 - (-t150 + t553) * t346 + t586 * t228;
t555 = qJ(5) * t168;
t551 = t170 * t168;
t549 = t463 * t346;
t548 = t196 * t346;
t547 = t217 * t228;
t546 = t228 * t226;
t545 = t276 * t346;
t543 = t278 * t344;
t541 = t280 * t344;
t540 = t339 * t349;
t535 = t343 * t344;
t32 = qJ(6) * t170 + t51;
t532 = qJD(5) - t32;
t251 = -pkin(2) * t556 - t283;
t268 = t276 * pkin(3);
t124 = -t277 * pkin(10) + t251 + t268;
t527 = pkin(2) * t538 + pkin(9) * t539;
t253 = -pkin(1) * t342 - t527;
t152 = t584 * t252 + t344 * t253;
t126 = -pkin(10) * t538 + t152;
t73 = t343 * t124 + t346 * t126;
t151 = -t344 * t252 + t584 * t253;
t240 = pkin(9) * t484 + t343 * t301;
t525 = t585 * pkin(1) + pkin(8) * t487;
t340 = t345 ^ 2;
t341 = t347 ^ 2;
t524 = t340 - t341;
t521 = qJD(2) * t345;
t503 = t347 * t540;
t502 = t344 * t538;
t59 = t276 * qJ(5) + t73;
t183 = t463 * pkin(3);
t501 = pkin(4) * t549 + t463 * t537 + t183;
t185 = t196 * pkin(3);
t500 = -pkin(4) * t548 - t196 * t537 - t185;
t125 = pkin(3) * t538 - t151;
t499 = -pkin(4) * t545 - t276 * t537 - t268;
t269 = t278 * pkin(2);
t491 = t278 * t584;
t498 = -pkin(3) * t491 - pkin(10) * t543 - t269;
t271 = t280 * pkin(2);
t490 = t280 * t584;
t497 = -pkin(3) * t490 - pkin(10) * t541 - t271;
t495 = -g(1) * t541 - g(2) * t543 + g(3) * t502;
t492 = t119 * t584;
t479 = t342 * t521;
t473 = pkin(1) * t616;
t466 = -t134 * pkin(4) + qJ(5) * t135;
t139 = t197 * t346 + t280 * t343;
t465 = -t138 * pkin(4) + qJ(5) * t139;
t191 = t277 * t346 - t343 * t538;
t464 = -t190 * pkin(4) + qJ(5) * t191;
t75 = t150 * t346 - t106;
t72 = t124 * t346 - t343 * t126;
t333 = pkin(9) * t486;
t239 = t346 * t301 - t333;
t460 = t217 * t346;
t456 = t345 * t503;
t455 = t344 * t146 - t584 * t153 + t202 * t476 + t214 * t519;
t263 = t342 * t406;
t88 = -t252 * t476 - t253 * t519 + t584 * t263 - t344 * t265;
t452 = t342 * t483;
t453 = pkin(3) * t452 + pkin(10) * t502 + t527;
t444 = t345 * t471;
t441 = -pkin(1) * t583 + pkin(8) * t489;
t438 = t342 * t349 * t556;
t436 = -g(1) * t134 + g(2) * t138;
t435 = g(1) * t135 - g(2) * t139;
t434 = g(1) * t463 + g(2) * t196;
t433 = -g(1) * t278 + g(2) * t280;
t432 = g(1) * t281 + g(2) * t279;
t430 = t425 * t584;
t429 = (qJD(4) * t168 - t67) * pkin(10);
t428 = t425 * t523;
t97 = t277 * t517 + t343 * t608 - t346 * t479;
t25 = t168 * t97 + t190 * t68;
t40 = -t392 * pkin(3) + t455;
t422 = t170 * t228 - t536;
t188 = qJD(3) * t277 + t342 * t478;
t61 = t113 * t276 + t188 * t217;
t416 = 0.2e1 * t467 + qJD(2);
t414 = t281 * pkin(2) + pkin(9) * t280 + t525;
t87 = -t252 * t519 + t253 * t476 + t344 * t263 + t584 * t265;
t83 = pkin(10) * t479 + t87;
t93 = t188 * pkin(3) - t189 * pkin(10) + t266;
t20 = -t124 * t518 - t126 * t517 - t343 * t83 + t346 * t93;
t205 = -t481 + t240;
t6 = t468 - t601;
t408 = qJD(3) * t430;
t8 = t100 * t517 - t103 * t518 + t343 * t45 + t346 * t39;
t403 = t197 * pkin(3) + t414;
t402 = t425 * t344;
t19 = t124 * t517 - t126 * t518 + t343 * t93 + t346 * t83;
t400 = t102 * t217 - t580;
t74 = -t464 + t125;
t397 = g(1) * t585 + g(2) * t583;
t159 = -t278 * t486 - t279 * t346;
t161 = -t280 * t486 - t281 * t346;
t229 = t343 * t452 - t346 * t539;
t396 = -g(1) * t161 - g(2) * t159 - g(3) * t229;
t160 = -t278 * t484 + t279 * t343;
t162 = -t280 * t484 + t281 * t343;
t395 = -g(1) * t162 - g(2) * t160 - g(3) * t230;
t394 = g(1) * t196 - g(2) * t463 + g(3) * t276;
t393 = g(1) * t197 + g(2) * t193 + g(3) * t277;
t10 = t68 * pkin(4) + t67 * qJ(5) - t170 * qJD(5) + t40;
t3 = -pkin(5) * t68 + qJDD(6) - t10;
t391 = t3 + t394;
t390 = t230 * pkin(4) + qJ(5) * t229 + t453;
t28 = t168 * t461 - t567;
t389 = -t279 * pkin(2) - t278 * pkin(9) + t441;
t388 = -pkin(10) * t567 - t393;
t4 = t8 + t621;
t385 = pkin(3) * t479 + t88;
t17 = t188 * qJ(5) + t276 * qJD(5) + t19;
t382 = -pkin(3) * t193 + t389;
t381 = t160 * pkin(4) + pkin(9) * t279 + qJ(5) * t159 + t498;
t380 = t162 * pkin(4) + pkin(9) * t281 + qJ(5) * t161 + t497;
t378 = t67 - t552;
t98 = -t277 * t518 + t343 * t479 + t346 * t608;
t7 = t168 * t98 + t170 * t97 - t190 * t67 + t191 * t68;
t375 = t139 * pkin(4) + qJ(5) * t138 + t403;
t374 = t113 * t190 + t168 * t188 + t217 * t97 + t276 * t68;
t14 = t113 * t191 + t170 * t188 + t217 * t98 - t276 * t67;
t373 = -pkin(10) * qJD(4) * t217 + t394;
t370 = t344 * t455 - t399 * t584 - t432;
t368 = -t10 + t373;
t367 = -t373 + t40;
t69 = -t113 * t584 - t217 * t402;
t366 = -pkin(4) * t135 - qJ(5) * t134 + t382;
t363 = g(1) * t139 + g(2) * t135 + g(3) * t191 - t8;
t362 = qJ(5) * t98 + qJD(5) * t191 + t385;
t360 = t170 * t50 - t364;
t22 = t68 * t535 + (t344 * t517 + t410) * t168;
t359 = t217 * t51 + t363;
t355 = t551 - t113;
t351 = -t228 * t517 + t343 * t458 - t352;
t350 = -t351 - t550;
t308 = t576 * t343;
t296 = -pkin(3) - t427;
t287 = pkin(3) - t592;
t264 = t526 * qJD(1);
t254 = (pkin(9) + t426) * t344;
t206 = -t239 + t337;
t203 = (-pkin(9) - t593) * t344;
t179 = qJ(6) * t535 + t205;
t174 = t333 + (-qJ(6) * t344 - t301) * t346 - t606;
t92 = pkin(4) * t170 + t555;
t70 = -t170 * t586 - t555;
t60 = -pkin(4) * t276 - t72;
t58 = -pkin(4) * t228 - t75;
t54 = -t74 - t581;
t46 = qJ(6) * t190 + t59;
t41 = -pkin(4) * t217 + t531;
t38 = -qJ(6) * t191 - t276 * t586 - t72;
t30 = t217 * t460 - t422;
t27 = t170 * t460 - t568;
t26 = -t217 * t586 + t532;
t24 = t170 * t98 - t191 * t67;
t23 = pkin(4) * t97 - t362;
t21 = -t67 * t534 + (-t344 * t518 + t409) * t170;
t18 = -pkin(4) * t188 - t20;
t16 = -t586 * t97 + t362;
t13 = qJ(6) * t97 + qJD(6) * t190 + t17;
t11 = -qJ(6) * t98 - qJD(6) * t191 - t188 * t586 - t20;
t2 = t4 + t615;
t1 = -qJD(6) * t170 + t572 - t582 + t6;
t5 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t583 - g(2) * t585, t397, 0, 0 (qJDD(1) * t340 + 0.2e1 * t444) * t339 (t345 * t511 - t513 * t524) * t616, t342 * t416 * t520 + t345 * t605 (qJDD(1) * t341 - 0.2e1 * t444) * t339, t347 * t605 - t416 * t479, t412 * t556, -t266 * t417 + t283 * t412 - t454 * t556 + g(1) * t279 - g(2) * t281 + (-t472 + t511) * t473, -t265 * t417 - t526 * t412 - t172 * t556 + (-t471 - t512) * t473 + t433 ((-t261 * qJD(2) + qJDD(1) * t526 + t172) * t347 + (-qJD(2) * t264 - qJDD(1) * t283 + t454) * t345 - t397) * t342, t339 * qJDD(1) * pkin(1) ^ 2 - g(1) * t441 - g(2) * t525 + t172 * t526 - t261 * t266 + t264 * t265 - t283 * t454, t228 * t189 - t277 * t620, -t277 * t119 - t228 * t188 - t189 * t226 + t276 * t620, -t189 * t425 + t277 * t509 + (t228 * t521 + t277 * t472 + t347 * t620) * t342, t119 * t276 + t188 * t226, t188 * t425 - t276 * t509 + (t119 * t347 + (-qJD(1) * t276 - t226) * t521) * t342 (-t509 * t347 + (-t324 - t425) * t521) * t342, -t88 * t425 + t151 * t509 + t266 * t226 + t251 * t119 + t147 * t276 + t201 * t188 + g(1) * t193 - g(2) * t197 + (t455 * t347 + (qJD(1) * t151 + t117) * t521) * t342, -t118 * t479 + t147 * t277 - t152 * t392 + t201 * t189 + t266 * t228 - t251 * t620 - t399 * t538 + t425 * t87 + t434, -t117 * t189 - t118 * t188 - t152 * t119 + t151 * t620 - t87 * t226 - t88 * t228 + t276 * t399 + t277 * t455 - t433, -g(1) * t389 - g(2) * t414 + t117 * t88 + t118 * t87 + t147 * t251 - t151 * t455 - t152 * t399 + t201 * t266, t24, -t7, t14, t25, -t374, t61, t102 * t97 + t113 * t72 + t125 * t68 - t168 * t385 + t188 * t51 + t190 * t40 + t20 * t217 - t276 * t468 + t435, t102 * t98 - t113 * t73 - t125 * t67 - t170 * t385 - t188 * t52 - t19 * t217 + t191 * t40 - t276 * t8 + t436, -t168 * t19 - t170 * t20 - t190 * t8 + t191 * t468 - t51 * t98 - t52 * t97 + t67 * t72 - t68 * t73 - t434, t8 * t73 + t52 * t19 - t468 * t72 + t51 * t20 + t40 * t125 - t102 * t385 - g(1) * (t382 + t577) - g(2) * (t403 + t578) t24, t14, t7, t61, t374, t25, t10 * t190 - t113 * t60 + t168 * t23 - t18 * t217 - t188 * t41 - t276 * t6 + t50 * t97 + t68 * t74 + t435, -t168 * t17 + t170 * t18 - t190 * t4 + t191 * t6 + t41 * t98 - t42 * t97 - t59 * t68 - t60 * t67 - t434, -t10 * t191 + t113 * t59 + t17 * t217 - t170 * t23 + t188 * t42 + t276 * t4 - t50 * t98 + t67 * t74 - t436, t4 * t59 + t42 * t17 + t10 * t74 + t50 * t23 + t6 * t60 + t41 * t18 - g(1) * (t366 + t577) - g(2) * (t375 + t578) t24, t7, -t14, t25, -t374, t61, -t1 * t276 - t11 * t217 - t113 * t38 - t16 * t168 - t188 * t26 - t190 * t3 - t31 * t97 - t54 * t68 + t435, t113 * t46 + t13 * t217 + t16 * t170 + t188 * t29 + t191 * t3 + t2 * t276 + t31 * t98 - t54 * t67 - t436, -t1 * t191 - t11 * t170 + t13 * t168 + t190 * t2 - t26 * t98 + t29 * t97 + t38 * t67 + t46 * t68 + t434, t2 * t46 + t29 * t13 + t1 * t38 + t26 * t11 + t3 * t54 + t31 * t16 - g(1) * (-pkin(5) * t135 + t463 * t576 + t366) - g(2) * (pkin(5) * t139 + t196 * t576 + t375); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t456, t524 * t540, -t347 * t438 + t470, t456, t345 * t438 + t323, t412, pkin(1) * t345 * t540 + g(1) * t280 + g(2) * t278 - g(3) * t538 + t264 * t417 - t454, pkin(1) * t503 + t261 * t417 + (pkin(8) * t513 + g(3)) * t539 + t432 - t496, 0, 0, -t228 * t604 - t623, -t344 * t119 + t226 * t604 + t228 * t603 - t584 * t620, -t408 + t344 * t509 + (t347 * t430 + (qJD(2) * t344 - t228) * t345) * t523, -t226 * t402 - t492, t584 * t392 + t226 * t480 + (qJD(3) * t425 - t347 * t428) * t344, t345 * t428, pkin(9) * t408 - pkin(2) * t119 - t147 * t584 + t155 * t425 - t264 * t226 + g(1) * t490 + g(2) * t491 + (-g(3) * t483 - t117 * t522) * t342 + (-pkin(9) * t392 - t201 * t425) * t344, pkin(2) * t620 + t118 * t480 + t147 * t344 - t201 * t604 - t264 * t228 - t392 * t507 - t425 * t609 + t495, -g(3) * t539 + t155 * t228 + t370 + t609 * t226 + t603 * t118 + t604 * t117 + (t228 * t476 - t492 - t623) * pkin(9), -t147 * pkin(2) - t118 * t156 - t117 * t155 - t201 * t264 + g(1) * t271 + g(2) * t269 - g(3) * t527 + ((-t117 * t584 - t118 * t344) * qJD(3) + t370) * pkin(9), t21, t618, -t617, t22, t356, t69, t468 * t584 - t102 * t212 + t239 * t113 + t401 * t168 + t561 * t217 + (t102 * t486 + t168 * t507) * qJD(3) + (pkin(9) * t68 + t102 * t517 + t343 * t40 - t425 * t51) * t344 + t395, t8 * t584 - t102 * t213 - t240 * t113 + t401 * t170 - t562 * t217 + (t102 * t484 + t170 * t507) * qJD(3) + (-pkin(9) * t67 - t102 * t518 + t346 * t40 + t425 * t52) * t344 - t396, t52 * t212 + t51 * t213 + t239 * t67 - t240 * t68 - t561 * t170 - t562 * t168 + (-t484 * t51 - t486 * t52) * qJD(3) + (-t343 * t8 + t346 * t468 + (t343 * t51 - t346 * t52) * qJD(4)) * t344 - t495, t8 * t240 - t468 * t239 + t102 * t401 - g(1) * t497 - g(2) * t498 - g(3) * t453 + t562 * t52 + t561 * t51 + (t102 * t476 + t40 * t344 - t432) * pkin(9), t21, -t617, -t618, t69, -t356, t22, t6 * t584 - t206 * t113 + t254 * t68 + t410 * t50 - t564 * t217 + t563 * t168 + (t10 * t343 + t41 * t425 + t50 * t517) * t344 + t395, -t205 * t68 - t206 * t67 + t42 * t212 - t41 * t213 + t564 * t170 - t565 * t168 + (t41 * t484 - t42 * t486) * qJD(3) + (-t343 * t4 + t346 * t6 + (-t343 * t41 - t346 * t42) * qJD(4)) * t344 - t495, -t4 * t584 + t205 * t113 + t254 * t67 - t409 * t50 + t565 * t217 - t563 * t170 + (-t10 * t346 - t42 * t425 + t50 * t518) * t344 + t396, -g(1) * t380 - g(2) * t381 - g(3) * t390 + t10 * t254 + t4 * t205 + t6 * t206 + t41 * t564 + t42 * t565 + t50 * t563, t21, -t618, t617, t22, t356, t69, t1 * t584 - t174 * t113 - t203 * t68 - t410 * t31 - t575 * t217 - t566 * t168 + (t26 * t425 - t3 * t343 - t31 * t517) * t344 + t395, -t2 * t584 + t179 * t113 - t203 * t67 + t409 * t31 + t574 * t217 + t566 * t170 + (-t29 * t425 + t3 * t346 - t31 * t518) * t344 + t396, t174 * t67 + t179 * t68 - t29 * t212 + t26 * t213 - t575 * t170 + t574 * t168 + (-t26 * t484 + t29 * t486) * qJD(3) + (-t1 * t346 + t2 * t343 + (t26 * t343 + t29 * t346) * qJD(4)) * t344 + t495, t2 * t179 + t1 * t174 + t3 * t203 - g(1) * (pkin(5) * t162 + qJ(6) * t541 + t380) - g(2) * (pkin(5) * t160 + qJ(6) * t543 + t381) - g(3) * (pkin(5) * t230 - qJ(6) * t502 + t390) + t566 * t31 + t574 * t29 + t575 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t546, -t226 ^ 2 + t228 ^ 2, -t226 * t425 - t620, -t546, -t228 * t425 - t119, t392, -t118 * t425 - t201 * t228 + t394 - t455, -t117 * t425 + t201 * t226 + t393 + t399, 0, 0, t27, -t12, t30, t28, -t619, -t547, -pkin(3) * t68 - t118 * t168 - t217 * t75 - t228 * t51 + t343 * t400 - t346 * t367, pkin(3) * t67 - t118 * t170 + t217 * t76 + t228 * t52 + t343 * t367 + t346 * t400, t168 * t76 + t170 * t75 + (-t226 * t51 + t8 + (-t51 + t579) * qJD(4)) * t346 + (t429 + t468 - t569) * t343 + t388, -t40 * pkin(3) + g(1) * t185 - g(2) * t183 + g(3) * t268 - t102 * t118 - t51 * t75 - t52 * t76 + (t468 * t343 + t8 * t346 + (-t343 * t52 - t346 * t51) * qJD(4) - t393) * pkin(10), t27, t30, t12, -t547, t619, t28, t558 * t168 + t217 * t58 + t228 * t41 + t296 * t68 + t343 * t599 + t368 * t346, t168 * t56 - t170 * t58 + (t226 * t41 + t4 + (t41 + t579) * qJD(4)) * t346 + (t429 + t6 - t570) * t343 + t388, -t558 * t170 - t217 * t56 - t228 * t42 + t296 * t67 + t368 * t343 - t346 * t599, t10 * t296 - t42 * t56 - t41 * t58 - g(1) * t500 - g(2) * t501 - g(3) * t499 + t558 * t50 + (t6 * t343 + t4 * t346 + (-t343 * t42 + t346 * t41) * qJD(4) - t393) * pkin(10), t27, t12, -t215 * t346 + t422, t28, -t619, -t547, -t113 * t308 - t168 * t560 - t217 * t557 + t228 * t26 - t287 * t68 - t31 * t461 + t346 * t391, t113 * t309 + t170 * t560 + t217 * t559 - t228 * t29 - t287 * t67 + t31 * t460 + t343 * t391, t308 * t67 + t309 * t68 - t557 * t170 + t559 * t168 + (-t217 * t26 - t2) * t346 + (-t1 + t571) * t343 + t393, t2 * t309 + t1 * t308 + t3 * t287 - g(1) * (-pkin(5) * t548 + t197 * t576 + t500) - g(2) * (pkin(5) * t549 + t193 * t576 + t501) - g(3) * (-pkin(5) * t545 + t277 * t576 + t499) + t560 * t31 + t559 * t29 + t557 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t551, t602, -t378, -t551, -t350, t113, -t102 * t170 + t372 + t569, t102 * t168 + t359, 0, 0, t551, -t378, -t602, t113, t350, -t551, -t168 * t92 + t110 - t360 + t569, pkin(4) * t67 - t573 + (t42 - t52) * t170 + (t41 - t531) * t168, -t168 * t50 + t170 * t92 - t359 + t611, -t6 * pkin(4) - g(1) * t465 - g(2) * t466 - g(3) * t464 + t4 * qJ(5) - t41 * t52 + t42 * t531 - t50 * t92, t551, -t602, t378, -t551, -t350, t113, t168 * t70 + t217 * t33 + (pkin(5) + t586) * t113 + t597, t168 * t31 - t170 * t70 - t217 * t32 - t363 + t611 + t615, t573 - t586 * t67 + (-t29 + t33) * t170 + (-t26 + t532) * t168, t2 * qJ(5) - t1 * t586 - t26 * t33 - t31 * t70 - g(1) * (-pkin(5) * t138 + t465) - g(2) * (-pkin(5) * t134 + t466) - g(3) * (t464 - t581) + t532 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t355, -t378, t607, t360 - t570, 0, 0, 0, 0, 0, 0, t355, t607, t378, -t571 - t582 - t597; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t351 - t550, -t613, -t587 - t588, -t168 * t29 + t170 * t26 + t391;];
tau_reg  = t5;