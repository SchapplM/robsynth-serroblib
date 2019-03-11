% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:23:21
% EndTime: 2019-03-09 22:24:02
% DurationCPUTime: 25.76s
% Computational Cost: add. (34453->888), mult. (78494->1138), div. (0->0), fcn. (57990->18), ass. (0->393)
t407 = sin(qJ(2));
t410 = cos(qJ(2));
t456 = pkin(2) * t407 - pkin(8) * t410;
t329 = t456 * qJD(1);
t409 = cos(qJ(3));
t406 = sin(qJ(3));
t503 = qJD(1) * t407;
t475 = t406 * t503;
t258 = pkin(7) * t475 + t409 * t329;
t518 = t409 * t410;
t445 = pkin(3) * t407 - pkin(9) * t518;
t412 = -pkin(9) - pkin(8);
t482 = qJD(3) * t412;
t627 = -qJD(1) * t445 + t409 * t482 - t258;
t305 = t406 * t329;
t522 = t407 * t409;
t524 = t406 * t410;
t626 = t305 + (-pkin(7) * t522 - pkin(9) * t524) * qJD(1) - t406 * t482;
t499 = qJD(2) * t409;
t321 = -t475 + t499;
t481 = t409 * t503;
t501 = qJD(2) * t406;
t322 = t481 + t501;
t405 = sin(qJ(4));
t575 = cos(qJ(4));
t234 = -t575 * t321 + t322 * t405;
t402 = sin(pkin(11));
t403 = cos(pkin(11));
t441 = t405 * t321 + t322 * t575;
t161 = t234 * t403 + t402 * t441;
t404 = sin(qJ(6));
t502 = qJD(1) * t410;
t464 = qJD(3) + t502;
t446 = t464 * qJD(2);
t485 = t407 * qJDD(1);
t430 = t446 + t485;
t495 = qJD(3) * t407;
t469 = qJD(1) * t495;
t351 = t406 * t469;
t486 = t406 * qJDD(2);
t463 = -t351 + t486;
t422 = t430 * t409 + t463;
t362 = t406 * t485;
t488 = qJD(1) * qJD(2);
t470 = t410 * t488;
t601 = qJD(2) * qJD(3) + t470;
t462 = t406 * t601 + t409 * t469 + t362;
t437 = t409 * qJDD(2) - t462;
t472 = t575 * qJD(4);
t493 = qJD(4) * t405;
t117 = -t321 * t472 + t322 * t493 - t405 * t437 - t575 * t422;
t118 = qJD(4) * t441 + t405 * t422 - t575 * t437;
t467 = t117 * t402 - t403 * t118;
t574 = cos(qJ(6));
t474 = qJD(6) * t574;
t492 = qJD(6) * t404;
t595 = -t234 * t402 + t403 * t441;
t66 = -t117 * t403 - t118 * t402;
t20 = t161 * t474 - t404 * t467 + t492 * t595 - t574 * t66;
t369 = -qJD(3) + t502;
t355 = -qJD(4) + t369;
t344 = -qJD(6) + t355;
t85 = t161 * t574 + t404 * t595;
t551 = t344 * t85;
t625 = -t20 - t551;
t526 = t405 * t406;
t440 = t409 * t575 - t526;
t581 = qJD(3) + qJD(4);
t584 = t575 * qJD(3) + t472;
t509 = -t409 * t584 + t440 * t502 + t526 * t581;
t324 = t405 * t409 + t406 * t575;
t248 = t581 * t324;
t508 = -t324 * t502 + t248;
t558 = t85 ^ 2;
t615 = -t404 * t161 + t574 * t595;
t559 = t615 ^ 2;
t620 = -t558 + t559;
t557 = t85 * t615;
t348 = t412 * t406;
t349 = t412 * t409;
t512 = t348 * t472 + t349 * t493 + t405 * t627 - t626 * t575;
t252 = t405 * t348 - t575 * t349;
t511 = -qJD(4) * t252 + t626 * t405 + t575 * t627;
t457 = pkin(2) * t410 + pkin(8) * t407;
t340 = -pkin(1) - t457;
t312 = t340 * qJD(1);
t384 = pkin(7) * t502;
t347 = qJD(2) * pkin(8) + t384;
t243 = t409 * t312 - t347 * t406;
t202 = -pkin(9) * t322 + t243;
t189 = -pkin(3) * t369 + t202;
t244 = t312 * t406 + t347 * t409;
t203 = pkin(9) * t321 + t244;
t195 = t405 * t203;
t126 = t575 * t189 - t195;
t591 = qJ(5) * t441;
t97 = t126 - t591;
t90 = -pkin(4) * t355 + t97;
t197 = t575 * t203;
t127 = t405 * t189 + t197;
t602 = t234 * qJ(5);
t98 = t127 - t602;
t91 = t402 * t98;
t45 = t403 * t90 - t91;
t592 = pkin(10) * t595;
t34 = -pkin(5) * t355 + t45 - t592;
t550 = t403 * t98;
t46 = t402 * t90 + t550;
t608 = pkin(10) * t161;
t38 = t46 - t608;
t387 = t410 * qJDD(1);
t585 = -t407 * t488 + t387;
t318 = qJDD(3) - t585;
t311 = qJDD(4) + t318;
t332 = t456 * qJD(2);
t249 = qJD(1) * t332 + qJDD(1) * t340;
t294 = pkin(7) * t585 + qJDD(2) * pkin(8);
t494 = qJD(3) * t409;
t496 = qJD(3) * t406;
t148 = t406 * t249 + t409 * t294 + t312 * t494 - t347 * t496;
t112 = pkin(9) * t437 + t148;
t239 = t409 * t249;
t427 = t485 + t601;
t466 = -qJD(3) * t312 - t294;
t99 = t318 * pkin(3) + t351 * pkin(9) + t239 + (-pkin(9) * qJDD(2) + t466) * t406 + (-pkin(9) * t427 - qJD(3) * t347) * t409;
t37 = -qJD(4) * t127 - t405 * t112 + t575 * t99;
t26 = t311 * pkin(4) + t117 * qJ(5) - qJD(5) * t441 + t37;
t465 = -t575 * t112 - t189 * t472 + t203 * t493 - t405 * t99;
t30 = -qJ(5) * t118 - qJD(5) * t234 - t465;
t8 = t403 * t26 - t30 * t402;
t6 = pkin(5) * t311 - pkin(10) * t66 + t8;
t9 = t402 * t26 + t403 * t30;
t7 = pkin(10) * t467 + t9;
t1 = t34 * t474 - t38 * t492 + t404 * t6 + t574 * t7;
t346 = -qJD(2) * pkin(2) + pkin(7) * t503;
t260 = -pkin(3) * t321 + t346;
t182 = pkin(4) * t234 + qJD(5) + t260;
t103 = pkin(5) * t161 + t182;
t401 = qJ(3) + qJ(4);
t388 = pkin(11) + t401;
t378 = qJ(6) + t388;
t363 = sin(t378);
t364 = cos(t378);
t411 = cos(qJ(1));
t408 = sin(qJ(1));
t519 = t408 * t410;
t254 = t363 * t411 - t364 * t519;
t516 = t410 * t411;
t256 = t363 * t408 + t364 * t516;
t560 = g(3) * t407;
t617 = g(1) * t256 - g(2) * t254 + t103 * t85 + t364 * t560 - t1;
t21 = qJD(6) * t615 + t404 * t66 - t467 * t574;
t549 = t615 * t344;
t613 = -t21 - t549;
t624 = -pkin(4) * t503 + qJ(5) * t509 - t324 * qJD(5) + t511;
t623 = -qJ(5) * t508 + qJD(5) * t440 + t512;
t382 = pkin(7) * t485;
t539 = qJDD(2) * pkin(2);
t295 = pkin(7) * t470 + t382 - t539;
t397 = g(3) * t410;
t455 = g(1) * t411 + g(2) * t408;
t431 = -t455 * t407 + t397;
t426 = -t295 - t431;
t622 = qJD(3) * pkin(8) * t369 + t426 + t539;
t536 = t161 ^ 2;
t538 = t595 ^ 2;
t621 = -t536 + t538;
t13 = t404 * t34 + t38 * t574;
t2 = -qJD(6) * t13 - t404 * t7 + t574 * t6;
t253 = t363 * t519 + t364 * t411;
t255 = -t363 * t516 + t364 * t408;
t610 = -g(1) * t255 + g(2) * t253 - t103 * t615 + t363 * t560 + t2;
t535 = t161 * t355;
t619 = t66 - t535;
t603 = t161 * t595;
t514 = -t402 * t509 + t403 * t508;
t513 = -t402 * t508 - t403 * t509;
t537 = t595 * t355;
t618 = t467 - t537;
t373 = sin(t388);
t374 = cos(t388);
t263 = t373 * t411 - t374 * t519;
t265 = t373 * t408 + t374 * t516;
t616 = g(1) * t265 - g(2) * t263 + t161 * t182 + t374 * t560 - t9;
t614 = pkin(5) * t595;
t554 = -t402 * t623 + t403 * t624;
t553 = t402 * t624 + t403 * t623;
t135 = -t202 * t405 - t197;
t104 = t135 + t602;
t136 = t575 * t202 - t195;
t105 = t136 - t591;
t527 = t403 * t405;
t552 = pkin(3) * qJD(4);
t542 = -t403 * t104 + t105 * t402 + (-t402 * t575 - t527) * t552;
t528 = t402 * t405;
t541 = -t402 * t104 - t403 * t105 + (t403 * t575 - t528) * t552;
t262 = t373 * t519 + t374 * t411;
t264 = -t373 * t516 + t374 * t408;
t609 = -g(1) * t264 + g(2) * t262 - t595 * t182 + t373 * t560 + t8;
t607 = -pkin(5) * t503 - pkin(10) * t513 + t554;
t606 = -pkin(10) * t514 + t553;
t605 = t542 - t608;
t604 = t541 + t592;
t534 = t234 * t441;
t573 = pkin(3) * t406;
t458 = pkin(3) * t496 - t502 * t573 - t384;
t498 = qJD(2) * t410;
t480 = t406 * t498;
t600 = t407 * t494 + t480;
t599 = -t234 ^ 2 + t441 ^ 2;
t598 = -t234 * t355 - t117;
t389 = sin(t401);
t390 = cos(t401);
t278 = t389 * t411 - t390 * t519;
t280 = t389 * t408 + t390 * t516;
t596 = g(1) * t280 - g(2) * t278 + t234 * t260 + t390 * t560 + t465;
t593 = pkin(4) * t441;
t586 = t243 * t369 + t148;
t510 = pkin(4) * t508 + t458;
t320 = t409 * t340;
t570 = pkin(7) * t406;
t241 = -pkin(9) * t522 + t320 + (-pkin(3) - t570) * t410;
t371 = pkin(7) * t518;
t268 = t406 * t340 + t371;
t525 = t406 * t407;
t250 = -pkin(9) * t525 + t268;
t176 = t405 * t241 + t575 * t250;
t300 = t406 * t519 + t409 * t411;
t302 = -t406 * t516 + t408 * t409;
t582 = -g(1) * t302 + g(2) * t300;
t277 = t389 * t519 + t390 * t411;
t279 = -t389 * t516 + t390 * t408;
t580 = -g(1) * t279 + g(2) * t277 + t389 * t560;
t579 = -t260 * t441 + t37 + t580;
t578 = -t355 * t441 - t118;
t577 = -0.2e1 * pkin(1);
t576 = t467 * pkin(5);
t572 = pkin(4) * t389;
t571 = pkin(4) * t402;
t569 = pkin(8) * t318;
t566 = g(1) * t408;
t561 = g(2) * t411;
t392 = t407 * pkin(7);
t251 = t575 * t348 + t349 * t405;
t213 = -qJ(5) * t324 + t251;
t214 = qJ(5) * t440 + t252;
t143 = t403 * t213 - t214 * t402;
t232 = t324 * t403 + t402 * t440;
t114 = -pkin(10) * t232 + t143;
t144 = t402 * t213 + t403 * t214;
t231 = t324 * t402 - t403 * t440;
t115 = -pkin(10) * t231 + t144;
t60 = t114 * t574 - t404 * t115;
t556 = qJD(6) * t60 + t404 * t607 + t574 * t606;
t61 = t404 * t114 + t115 * t574;
t555 = -qJD(6) * t61 - t404 * t606 + t574 * t607;
t394 = t409 * pkin(3);
t381 = t394 + pkin(2);
t461 = t575 * t498;
t180 = t248 * t407 + t405 * t480 - t409 * t461;
t288 = t440 * t407;
t500 = qJD(2) * t407;
t507 = t409 * t332 + t500 * t570;
t172 = t445 * qJD(2) + (-t371 + (pkin(9) * t407 - t340) * t406) * qJD(3) + t507;
t199 = t406 * t332 + t340 * t494 + (-t407 * t499 - t410 * t496) * pkin(7);
t179 = -pkin(9) * t600 + t199;
t77 = -qJD(4) * t176 + t575 * t172 - t405 * t179;
t49 = pkin(4) * t500 + t180 * qJ(5) - t288 * qJD(5) + t77;
t478 = t406 * t495;
t181 = t406 * t461 - t405 * t478 - t493 * t525 + (t405 * t498 + t407 * t584) * t409;
t287 = t324 * t407;
t76 = t405 * t172 + t575 * t179 + t241 * t472 - t250 * t493;
t53 = -qJ(5) * t181 - qJD(5) * t287 + t76;
t23 = t402 * t49 + t403 * t53;
t52 = t403 * t97 - t91;
t166 = -t404 * t231 + t232 * t574;
t548 = -qJD(6) * t166 - t404 * t513 - t514 * t574;
t547 = t231 * t474 + t232 * t492 + t404 * t514 - t513 * t574;
t380 = pkin(3) * t575 + pkin(4);
t290 = -pkin(3) * t528 + t403 * t380;
t283 = pkin(5) + t290;
t292 = pkin(3) * t527 + t380 * t402;
t209 = t283 * t574 - t404 * t292;
t546 = qJD(6) * t209 + t404 * t605 + t574 * t604;
t210 = t404 * t283 + t292 * t574;
t545 = -qJD(6) * t210 - t404 * t604 + t574 * t605;
t375 = pkin(4) * t403 + pkin(5);
t291 = t375 * t574 - t404 * t571;
t51 = -t402 * t97 - t550;
t39 = t51 + t608;
t40 = t52 - t592;
t544 = t291 * qJD(6) - t404 * t39 - t40 * t574;
t293 = t404 * t375 + t571 * t574;
t543 = -t293 * qJD(6) - t39 * t574 + t404 * t40;
t540 = pkin(7) * qJDD(1);
t532 = t244 * t369;
t531 = t321 * t369;
t530 = t322 * t321;
t529 = t322 * t369;
t523 = t406 * t411;
t521 = t407 * t411;
t520 = t407 * t412;
t517 = t410 * t369;
t515 = pkin(5) * t514 + t510;
t175 = t575 * t241 - t250 * t405;
t141 = -pkin(4) * t410 - qJ(5) * t288 + t175;
t150 = -qJ(5) * t287 + t176;
t79 = t402 * t141 + t403 * t150;
t377 = pkin(4) * t390;
t337 = t377 + t394;
t333 = pkin(3) * t525 + t392;
t506 = t411 * pkin(1) + t408 * pkin(7);
t399 = t407 ^ 2;
t400 = t410 ^ 2;
t505 = t399 - t400;
t504 = t399 + t400;
t497 = qJD(3) * t244;
t491 = t321 * qJD(3);
t490 = t322 * qJD(3);
t489 = t346 * qJD(3);
t398 = -qJ(5) + t412;
t414 = qJD(1) ^ 2;
t483 = t407 * t414 * t410;
t361 = pkin(5) * t374;
t282 = t361 + t337;
t261 = pkin(3) * t600 + pkin(7) * t498;
t479 = t369 * t496;
t476 = t355 * t503;
t22 = -t402 * t53 + t403 * t49;
t78 = t403 * t141 - t150 * t402;
t273 = -pkin(4) * t440 - t381;
t460 = t407 * t470;
t372 = t407 * t566;
t459 = -g(2) * t521 + t372;
t245 = pkin(4) * t287 + t333;
t192 = pkin(3) * t322 + t593;
t315 = -pkin(5) * t373 - t572;
t454 = pkin(7) * t321 - t346 * t406;
t275 = pkin(2) + t282;
t391 = -pkin(10) + t398;
t452 = t275 * t410 - t391 * t407;
t328 = pkin(2) + t337;
t450 = t328 * t410 - t398 * t407;
t448 = t381 * t410 - t520;
t151 = pkin(4) * t181 + t261;
t444 = -t410 * pkin(7) * t322 + t244 * t407;
t205 = -t287 * t402 + t288 * t403;
t67 = -pkin(5) * t410 - pkin(10) * t205 + t78;
t204 = t403 * t287 + t288 * t402;
t68 = -pkin(10) * t204 + t79;
t31 = -t404 * t68 + t574 * t67;
t32 = t404 * t67 + t574 * t68;
t442 = -pkin(7) * qJDD(2) + t488 * t577;
t138 = -t404 * t204 + t205 * t574;
t439 = t406 * t318 - t369 * t494;
t438 = pkin(1) * t414 + t455;
t413 = qJD(2) ^ 2;
t436 = pkin(7) * t413 + qJDD(1) * t577 + t561;
t429 = -t410 * t455 - t560;
t424 = t427 * t409;
t187 = -pkin(3) * t437 + t295;
t80 = t118 * pkin(4) + qJDD(5) + t187;
t417 = t80 + t431;
t395 = t411 * pkin(7);
t336 = t572 + t573;
t316 = t361 + t377;
t303 = t406 * t408 + t409 * t516;
t301 = -t408 * t518 + t523;
t296 = qJDD(6) + t311;
t281 = -t315 + t573;
t267 = -pkin(7) * t524 + t320;
t259 = -pkin(7) * t481 + t305;
t242 = -t311 * t410 - t355 * t500;
t200 = -qJD(3) * t268 + t507;
t188 = pkin(5) * t231 + t273;
t162 = t231 * t574 + t232 * t404;
t158 = pkin(5) * t204 + t245;
t149 = -t406 * t294 + t239 - t497;
t137 = t204 * t574 + t205 * t404;
t119 = t593 + t614;
t113 = t192 + t614;
t110 = -t180 * t403 - t181 * t402;
t109 = -t180 * t402 + t403 * t181;
t73 = pkin(5) * t109 + t151;
t42 = qJD(6) * t138 + t109 * t574 + t404 * t110;
t41 = t404 * t109 - t110 * t574 + t204 * t474 + t205 * t492;
t33 = t80 - t576;
t17 = -pkin(10) * t109 + t23;
t16 = pkin(5) * t500 - pkin(10) * t110 + t22;
t12 = t34 * t574 - t404 * t38;
t4 = -qJD(6) * t32 + t16 * t574 - t404 * t17;
t3 = qJD(6) * t31 + t404 * t16 + t17 * t574;
t5 = [0, 0, 0, 0, 0, qJDD(1), -t561 + t566, t455, 0, 0, qJDD(1) * t399 + 0.2e1 * t460, 0.2e1 * t387 * t407 - 0.2e1 * t488 * t505, qJDD(2) * t407 + t410 * t413, qJDD(1) * t400 - 0.2e1 * t460, qJDD(2) * t410 - t407 * t413, 0, t442 * t407 + (-t436 + t566) * t410, t407 * t436 + t410 * t442 - t372, 0.2e1 * t504 * t540 - t455, -g(1) * (-pkin(1) * t408 + t395) - g(2) * t506 + (pkin(7) ^ 2 * t504 + pkin(1) ^ 2) * qJDD(1), -t322 * t478 + (t322 * t498 + (t424 + t463) * t407) * t409 (t321 * t409 - t322 * t406) * t498 + ((-t463 - t491) * t406 + (-t406 * t446 - t362 + t437 - t490) * t409) * t407, -t463 * t410 + (qJD(2) * t322 + t479) * t407 + ((t318 - t387) * t407 + (-t369 - t464) * t498) * t409, -t321 * t600 - t437 * t525 (t369 * t501 - t437) * t410 + (qJD(2) * t321 - t439) * t407, -t318 * t410 - t369 * t500, -g(1) * t301 - g(2) * t303 - t200 * t369 + t267 * t318 + (-qJD(2) * t454 - t149) * t410 + (-pkin(7) * t437 + t243 * qJD(2) + t295 * t406 + t409 * t489) * t407, -g(1) * t300 - g(2) * t302 + t148 * t410 + t199 * t369 - t268 * t318 + (-t406 * t489 + t295 * t409 + (t409 * t485 + t463) * pkin(7)) * t407 + ((t346 * t410 + t392 * t464) * t409 - t444) * qJD(2), t199 * t321 + t268 * t437 - t200 * t322 - t267 * t463 + t372 + (-t244 * t524 + (-t243 * t410 - t267 * t464) * t409) * qJD(2) + (-t561 + (qJD(3) * t243 - t148) * t406 + (-qJDD(1) * t267 - t149 - t497) * t409) * t407, t148 * t268 + t244 * t199 + t149 * t267 + t243 * t200 - g(1) * t395 - g(2) * (t411 * t457 + t506) - t340 * t566 + (t295 * t407 + t346 * t498) * pkin(7), -t117 * t288 - t180 * t441, t117 * t287 - t118 * t288 + t180 * t234 - t181 * t441, t117 * t410 + t180 * t355 + t288 * t311 + t441 * t500, t118 * t287 + t181 * t234, t118 * t410 + t181 * t355 - t234 * t500 - t287 * t311, t242, -g(1) * t278 - g(2) * t280 + t118 * t333 + t126 * t500 + t175 * t311 + t181 * t260 + t187 * t287 + t234 * t261 - t355 * t77 - t37 * t410, -g(1) * t277 - g(2) * t279 - t117 * t333 - t127 * t500 - t176 * t311 - t180 * t260 + t187 * t288 + t261 * t441 + t355 * t76 - t410 * t465, t117 * t175 - t118 * t176 + t126 * t180 - t127 * t181 - t234 * t76 + t287 * t465 - t288 * t37 - t441 * t77 + t459, -t465 * t176 + t127 * t76 + t37 * t175 + t126 * t77 + t187 * t333 + t260 * t261 - g(1) * (pkin(3) * t523 + t395) - g(2) * (t381 * t516 - t411 * t520 + t506) + (-g(1) * (-pkin(1) - t448) - g(2) * t573) * t408, t110 * t595 + t205 * t66, -t109 * t595 - t110 * t161 - t204 * t66 + t205 * t467, -t110 * t355 + t205 * t311 - t410 * t66 + t500 * t595, t109 * t161 - t204 * t467, t109 * t355 - t161 * t500 - t204 * t311 - t410 * t467, t242, -g(1) * t263 - g(2) * t265 + t109 * t182 + t151 * t161 + t204 * t80 - t22 * t355 - t245 * t467 + t311 * t78 - t410 * t8 + t45 * t500, -g(1) * t262 - g(2) * t264 + t110 * t182 + t151 * t595 + t205 * t80 + t23 * t355 + t245 * t66 - t311 * t79 + t410 * t9 - t46 * t500, -t109 * t46 - t110 * t45 - t161 * t23 - t204 * t9 - t205 * t8 - t22 * t595 + t467 * t79 - t66 * t78 + t459, t9 * t79 + t46 * t23 + t8 * t78 + t45 * t22 + t80 * t245 + t182 * t151 - g(1) * (t336 * t411 + t395) - g(2) * (t328 * t516 - t398 * t521 + t506) + (-g(1) * (-pkin(1) - t450) - g(2) * t336) * t408, -t138 * t20 - t41 * t615, t137 * t20 - t138 * t21 + t41 * t85 - t42 * t615, t138 * t296 + t20 * t410 + t344 * t41 + t500 * t615, t137 * t21 + t42 * t85, -t137 * t296 + t21 * t410 + t344 * t42 - t500 * t85, -t296 * t410 - t344 * t500, -g(1) * t254 - g(2) * t256 + t103 * t42 + t12 * t500 + t137 * t33 + t158 * t21 - t2 * t410 + t296 * t31 - t344 * t4 + t73 * t85, -g(1) * t253 - g(2) * t255 + t1 * t410 - t103 * t41 - t13 * t500 + t138 * t33 - t158 * t20 - t296 * t32 + t3 * t344 + t615 * t73, -t1 * t137 + t12 * t41 - t13 * t42 - t138 * t2 + t20 * t31 - t21 * t32 - t3 * t85 - t4 * t615 + t459, t1 * t32 + t13 * t3 + t2 * t31 + t12 * t4 + t33 * t158 + t103 * t73 - g(1) * (t281 * t411 + t395) - g(2) * (t275 * t516 - t391 * t521 + t506) + (-g(1) * (-pkin(1) - t452) - g(2) * t281) * t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t483, t505 * t414, t485, t483, t387, qJDD(2), t407 * t438 - t382 - t397, t560 + (t438 - t540) * t410, 0, 0, t463 * t406 + (t406 * t430 - t529) * t409 (-t462 + t529) * t406 + (-t351 + t424 + 0.2e1 * t486 - t531) * t409 (-t322 * t407 + t409 * t517) * qJD(1) + t439, t406 * t531 + t409 * t437, t479 + t409 * t318 + (-t321 * t407 - t406 * t517) * qJD(1), t369 * t503, -pkin(2) * t462 + t258 * t369 + (t489 - t569) * t406 + (-t243 * t407 + t410 * t454) * qJD(1) + t622 * t409, pkin(2) * t351 - t259 * t369 + t444 * qJD(1) + (-pkin(2) * t427 - t346 * t369 - t569) * t409 - t622 * t406, t258 * t322 - t259 * t321 + (-t149 + t532 + (t463 - t491) * pkin(8)) * t406 + ((t406 * t427 + t437 + t490) * pkin(8) + t586) * t409 + t429, -t346 * t384 - t243 * t258 - t244 * t259 + t426 * pkin(2) + (t148 * t409 - t149 * t406 + (-t243 * t409 - t244 * t406) * qJD(3) + t429) * pkin(8), -t117 * t324 - t441 * t509, -t117 * t440 - t118 * t324 + t234 * t509 - t441 * t508, t311 * t324 + t355 * t509 - t441 * t503, -t118 * t440 + t234 * t508, t234 * t503 + t311 * t440 + t355 * t508, t476, -t118 * t381 - t126 * t503 - t187 * t440 + t234 * t458 + t251 * t311 + t260 * t508 - t355 * t511 - t390 * t431, t117 * t381 + t127 * t503 + t187 * t324 - t252 * t311 - t260 * t509 + t355 * t512 + t389 * t431 + t441 * t458, t117 * t251 - t118 * t252 + t126 * t509 - t127 * t508 - t234 * t512 - t324 * t37 - t440 * t465 - t441 * t511 + t429, -g(3) * t448 + t126 * t511 + t127 * t512 - t187 * t381 + t37 * t251 - t465 * t252 + t260 * t458 + t455 * (t381 * t407 + t410 * t412) t232 * t66 + t513 * t595, -t161 * t513 - t231 * t66 + t232 * t467 - t514 * t595, t232 * t311 - t355 * t513 - t503 * t595, t161 * t514 - t231 * t467, t161 * t503 - t231 * t311 + t355 * t514, t476, t143 * t311 + t161 * t510 + t182 * t514 + t231 * t80 - t273 * t467 - t355 * t554 - t374 * t431 - t45 * t503, -t144 * t311 + t182 * t513 + t232 * t80 + t273 * t66 + t355 * t553 + t373 * t431 + t46 * t503 + t510 * t595, -t143 * t66 + t144 * t467 - t161 * t553 - t231 * t9 - t232 * t8 - t45 * t513 - t46 * t514 - t554 * t595 + t429, -g(3) * t450 + t8 * t143 + t9 * t144 + t182 * t510 + t80 * t273 + t45 * t554 + t46 * t553 + t455 * (t328 * t407 + t398 * t410) -t166 * t20 - t547 * t615, t162 * t20 - t166 * t21 + t547 * t85 + t548 * t615, t166 * t296 + t344 * t547 - t503 * t615, t162 * t21 - t548 * t85, -t162 * t296 - t344 * t548 + t503 * t85, t344 * t503, -t103 * t548 - t12 * t503 + t162 * t33 + t188 * t21 + t296 * t60 - t344 * t555 - t364 * t431 + t515 * t85, -t103 * t547 + t13 * t503 + t166 * t33 - t188 * t20 - t296 * t61 + t344 * t556 + t363 * t431 + t515 * t615, -t1 * t162 + t12 * t547 + t13 * t548 - t166 * t2 + t20 * t60 - t21 * t61 - t555 * t615 - t556 * t85 + t429, -g(3) * t452 + t1 * t61 + t103 * t515 + t12 * t555 + t13 * t556 + t33 * t188 + t2 * t60 + t455 * (t275 * t407 + t391 * t410); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t530, -t321 ^ 2 + t322 ^ 2, t422 + t531, t530, t437 - t529, t318, -t347 * t494 - t532 - t322 * t346 + t239 + (t466 + t560) * t406 + t582, g(1) * t303 - g(2) * t301 + g(3) * t522 - t321 * t346 - t586, 0, 0, t534, t599, t598, -t534, t578, t311, t135 * t355 + (-t234 * t322 + t311 * t575 + t355 * t493) * pkin(3) + t579, -t136 * t355 + (-t311 * t405 - t322 * t441 + t355 * t472) * pkin(3) + t596, -t126 * t234 + t127 * t441 + t135 * t441 + t136 * t234 + (t575 * t117 - t118 * t405 + (-t234 * t575 + t405 * t441) * qJD(4)) * pkin(3), -t126 * t135 - t127 * t136 + (-t465 * t405 + t37 * t575 - t260 * t322 + g(3) * t525 + (-t126 * t405 + t127 * t575) * qJD(4) + t582) * pkin(3), t603, t621, t619, -t603, t618, t311, -t161 * t192 + t290 * t311 - t355 * t542 + t609, -t192 * t595 - t292 * t311 + t355 * t541 + t616, -t290 * t66 + t292 * t467 + (t46 - t542) * t595 + (-t541 - t45) * t161, t9 * t292 + t8 * t290 - t182 * t192 - g(1) * (-t336 * t516 + t337 * t408) - g(2) * (-t336 * t519 - t337 * t411) + t336 * t560 + t541 * t46 + t542 * t45, t557, t620, t625, -t557, t613, t296, -t113 * t85 + t209 * t296 - t344 * t545 + t610, -t113 * t615 - t210 * t296 + t344 * t546 + t617, t20 * t209 - t21 * t210 + (t13 - t545) * t615 + (-t12 - t546) * t85, t1 * t210 + t2 * t209 - t103 * t113 - g(1) * (-t281 * t516 + t282 * t408) - g(2) * (-t281 * t519 - t282 * t411) + t281 * t560 + t546 * t13 + t545 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t534, t599, t598, -t534, t578, t311, -t127 * t355 + t579, -t126 * t355 + t596, 0, 0, t603, t621, t619, -t603, t618, t311, t355 * t51 + (-t161 * t441 + t311 * t403) * pkin(4) + t609, -t355 * t52 + (-t311 * t402 - t441 * t595) * pkin(4) + t616 (t402 * t467 - t403 * t66) * pkin(4) + (t46 + t51) * t595 + (t52 - t45) * t161, -t45 * t51 - t46 * t52 + (-t182 * t441 + t9 * t402 + t8 * t403 + t580) * pkin(4), t557, t620, t625, -t557, t613, t296, -t119 * t85 + t291 * t296 - t344 * t543 + t610, -t119 * t615 - t293 * t296 + t344 * t544 + t617, t20 * t291 - t21 * t293 + (t13 - t543) * t615 + (-t12 - t544) * t85, t1 * t293 + t2 * t291 - t103 * t119 - g(1) * (t315 * t516 + t316 * t408) - g(2) * (t315 * t519 - t316 * t411) - t315 * t560 + t544 * t13 + t543 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t467 - t537, t66 + t535, -t536 - t538, t161 * t46 + t45 * t595 + t417, 0, 0, 0, 0, 0, 0, t21 - t549, -t20 + t551, -t558 - t559, t12 * t615 + t13 * t85 + t417 - t576; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t557, t620, t625, -t557, t613, t296, -t13 * t344 + t610, -t12 * t344 + t617, 0, 0;];
tau_reg  = t5;
