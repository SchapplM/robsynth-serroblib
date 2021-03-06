% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRR13
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR13_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR13_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_invdynJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:05:20
% EndTime: 2019-03-09 20:06:13
% DurationCPUTime: 25.62s
% Computational Cost: add. (32330->878), mult. (92599->1240), div. (0->0), fcn. (79451->18), ass. (0->386)
t390 = sin(qJ(6));
t394 = cos(qJ(6));
t396 = cos(qJ(3));
t389 = cos(pkin(7));
t397 = cos(qJ(2));
t387 = sin(pkin(6));
t523 = qJD(1) * t387;
t500 = t397 * t523;
t472 = t389 * t500;
t344 = t396 * t472;
t559 = cos(pkin(6));
t479 = t559 * qJDD(1);
t368 = t479 + qJDD(2);
t392 = sin(qJ(3));
t386 = sin(pkin(7));
t484 = t559 * qJD(1);
t443 = t484 + qJD(2);
t429 = t386 * t443;
t414 = qJD(3) * t429;
t512 = qJDD(1) * t397;
t491 = t387 * t512;
t470 = t389 * t491;
t393 = sin(qJ(2));
t501 = t393 * t523;
t475 = t392 * t501;
t545 = t386 * t392;
t513 = qJDD(1) * t393;
t492 = t387 * t513;
t514 = qJD(1) * qJD(2);
t493 = t397 * t514;
t588 = -t387 * t493 - t492;
t155 = qJD(3) * t344 + t368 * t545 + t392 * t470 + (-qJD(2) * t389 - qJD(3)) * t475 + (t414 - t588) * t396;
t494 = t393 * t514;
t469 = t387 * t494;
t257 = t368 * t389 + qJDD(3) + (t469 - t491) * t386;
t385 = sin(pkin(13));
t388 = cos(pkin(13));
t116 = t385 * t155 - t257 * t388;
t117 = t155 * t388 + t257 * t385;
t391 = sin(qJ(5));
t395 = cos(qJ(5));
t483 = t395 * t116 + t391 * t117;
t536 = t393 * t396;
t537 = t392 * t397;
t435 = t389 * t537 + t536;
t418 = t435 * t387;
t253 = qJD(1) * t418 + t392 * t429;
t304 = t386 * t500 - t389 * t443 - qJD(3);
t193 = t253 * t388 - t304 * t385;
t243 = t385 * t253;
t482 = t304 * t388 + t243;
t597 = t395 * t193 - t391 * t482;
t45 = qJD(5) * t597 + t483;
t43 = qJDD(6) + t45;
t184 = t395 * t482;
t122 = t193 * t391 + t184;
t120 = qJD(6) + t122;
t608 = t120 ^ 2;
t610 = -t390 * t608 + t394 * t43;
t503 = pkin(1) * t559;
t374 = t397 * t503;
t365 = qJD(1) * t374;
t572 = pkin(10) * t389;
t465 = t387 * (-pkin(9) - t572);
t440 = t393 * t465;
t282 = qJD(1) * t440 + t365;
t373 = t393 * t503;
t405 = t397 * t465 - t373;
t283 = t405 * qJD(1);
t573 = pkin(10) * t386;
t439 = pkin(2) * t393 - t397 * t573;
t317 = t439 * t523;
t520 = qJD(3) * t396;
t495 = t389 * t520;
t540 = t389 * t392;
t609 = pkin(2) * t495 - t396 * t282 - t283 * t540 - t317 * t545;
t251 = -t396 * t429 - t344 + t475;
t248 = qJD(5) + t251;
t607 = t122 * t248;
t204 = -t283 * t386 + t389 * t317;
t436 = t389 * t536 + t537;
t315 = t436 * t387;
t302 = qJD(1) * t315;
t535 = t396 * t397;
t538 = t392 * t393;
t434 = -t389 * t538 + t535;
t316 = t434 * t387;
t303 = qJD(1) * t316;
t606 = -pkin(3) * t302 + qJ(4) * t303 - t204 + (-qJD(4) * t392 + (pkin(3) * t392 - qJ(4) * t396) * qJD(3)) * t386;
t476 = t386 * t501;
t521 = qJD(3) * t392;
t498 = t386 * t521;
t605 = pkin(10) * t498 + qJ(4) * t476 - qJD(4) * t389 - t609;
t235 = t303 * t385 - t388 * t476;
t497 = t386 * t520;
t604 = -t385 * t497 + t235;
t96 = -t394 * t248 + t390 * t597;
t603 = t597 * t96;
t98 = t248 * t390 + t394 * t597;
t602 = t597 * t98;
t562 = t385 * t605 + t388 * t606;
t561 = t385 * t606 - t388 * t605;
t398 = cos(qJ(1));
t487 = t398 * t559;
t576 = sin(qJ(1));
t330 = t393 * t576 - t397 * t487;
t331 = t393 * t487 + t397 * t576;
t541 = t387 * t398;
t509 = t386 * t541;
t218 = t330 * t540 - t331 * t396 + t392 * t509;
t287 = -t330 * t386 + t389 * t541;
t382 = pkin(13) + qJ(5);
t379 = sin(t382);
t380 = cos(t382);
t159 = t218 * t380 + t287 * t379;
t539 = t389 * t396;
t215 = t330 * t539 + t331 * t392 + t396 * t509;
t601 = t159 * t390 + t215 * t394;
t600 = t159 * t394 - t215 * t390;
t381 = t387 ^ 2;
t599 = 0.2e1 * t381;
t236 = t303 * t388 + t385 * t476;
t571 = pkin(11) * t388;
t595 = pkin(4) * t302 - pkin(11) * t236 - (pkin(4) * t392 - t396 * t571) * t386 * qJD(3) - t562;
t594 = pkin(11) * t604 + t561;
t496 = t389 * t521;
t598 = pkin(2) * t496 + pkin(10) * t497 - t392 * t282 + t283 * t539;
t596 = pkin(1) * t599;
t325 = t385 * t545 - t388 * t389;
t327 = t385 * t389 + t388 * t545;
t246 = t395 * t325 + t327 * t391;
t342 = t385 * t391 - t395 * t388;
t532 = -qJD(5) * t246 + t235 * t391 - t236 * t395 - t342 * t497;
t247 = -t325 * t391 + t327 * t395;
t343 = t385 * t395 + t388 * t391;
t531 = qJD(5) * t247 - t395 * t235 - t236 * t391 + t343 * t497;
t589 = t302 - t498;
t528 = -(-pkin(3) * t501 - t317 * t396) * t386 + t598;
t165 = t342 * t251;
t328 = t342 * qJD(5);
t527 = -t328 - t165;
t526 = t248 * t343;
t400 = qJD(2) * t436 + qJD(3) * t435;
t544 = t386 * t396;
t156 = -t368 * t544 + t387 * (qJD(1) * t400 + t392 * t513) + t392 * t414 - t396 * t470;
t587 = t218 * t379 - t287 * t380;
t525 = pkin(2) * t540 + pkin(10) * t544;
t320 = qJ(4) * t389 + t525;
t321 = (-pkin(3) * t396 - qJ(4) * t392 - pkin(2)) * t386;
t222 = -t320 * t385 + t388 * t321;
t182 = -pkin(4) * t544 - pkin(11) * t327 + t222;
t223 = t388 * t320 + t385 * t321;
t196 = -pkin(11) * t325 + t223;
t530 = t391 * t182 + t395 * t196;
t570 = pkin(11) + qJ(4);
t355 = t570 * t385;
t356 = t570 * t388;
t444 = -t355 * t395 - t356 * t391;
t542 = t387 * t397;
t412 = pkin(9) * t542 + t373;
t244 = t412 * qJD(1) + (t429 + t472) * pkin(10);
t407 = pkin(2) * t559 + t440;
t249 = qJD(2) * pkin(2) + qJD(1) * t407 + t365;
t438 = pkin(2) * t397 + t393 * t573;
t432 = -pkin(1) - t438;
t309 = t432 * t387;
t296 = qJD(1) * t309;
t138 = -t392 * t244 + t396 * (t249 * t389 + t296 * t386);
t170 = pkin(3) * t253 + qJ(4) * t251;
t94 = -t138 * t385 + t388 * t170;
t75 = pkin(4) * t253 + t251 * t571 + t94;
t555 = t251 * t385;
t95 = t388 * t138 + t385 * t170;
t84 = pkin(11) * t555 + t95;
t586 = qJD(4) * t342 - qJD(5) * t444 + t391 * t75 + t395 * t84;
t292 = -t355 * t391 + t356 * t395;
t585 = -qJD(4) * t343 - qJD(5) * t292 + t391 * t84 - t395 * t75;
t486 = t559 * t386;
t275 = t392 * t486 + t418;
t326 = t386 * t542 - t389 * t559;
t210 = t275 * t385 + t326 * t388;
t211 = t275 * t388 - t326 * t385;
t144 = -t210 * t391 + t211 * t395;
t581 = t389 * t535 - t538;
t274 = -t581 * t387 - t396 * t486;
t267 = t274 * t394;
t584 = -t144 * t390 + t267;
t529 = -pkin(4) * t604 + t528;
t518 = qJD(5) * t395;
t519 = qJD(5) * t391;
t583 = -t182 * t518 + t196 * t519 + t595 * t391 - t594 * t395;
t582 = -t494 + t512;
t580 = -pkin(3) * t257 + qJDD(4);
t471 = qJD(2) * t503;
t441 = qJD(1) * t471;
t467 = pkin(1) * t479;
t504 = pkin(9) * t491 + t393 * t467 + t397 * t441;
t411 = -pkin(9) * t469 + t504;
t180 = (t582 * t389 * t387 + t368 * t386) * pkin(10) + t411;
t413 = -t393 * t441 + t397 * t467;
t431 = -t493 - t513;
t415 = t431 * pkin(9);
t185 = t368 * pkin(2) + (t431 * t572 + t415) * t387 + t413;
t422 = t439 * qJD(2);
t221 = (qJD(1) * t422 + qJDD(1) * t432) * t387;
t62 = t396 * t180 + t185 * t540 + t221 * t545 - t244 * t521 + t249 * t495 + t296 * t497;
t52 = qJ(4) * t257 - qJD(4) * t304 + t62;
t130 = -t185 * t386 + t389 * t221;
t55 = pkin(3) * t156 - qJ(4) * t155 - qJD(4) * t253 + t130;
t23 = -t385 * t52 + t388 * t55;
t15 = pkin(4) * t156 - pkin(11) * t117 + t23;
t24 = t385 * t55 + t388 * t52;
t17 = -pkin(11) * t116 + t24;
t183 = -t249 * t386 + t389 * t296;
t113 = pkin(3) * t251 - qJ(4) * t253 + t183;
t139 = t396 * t244 + t249 * t540 + t296 * t545;
t119 = -qJ(4) * t304 + t139;
t69 = t388 * t113 - t119 * t385;
t49 = pkin(4) * t251 - pkin(11) * t193 + t69;
t70 = t385 * t113 + t388 * t119;
t57 = -pkin(11) * t482 + t70;
t22 = t391 * t49 + t395 * t57;
t6 = -qJD(5) * t22 + t395 * t15 - t17 * t391;
t153 = qJDD(5) + t156;
t44 = -qJD(5) * t184 - t391 * t116 + t395 * t117 - t193 * t519;
t26 = qJD(6) * t98 - t394 * t153 + t390 * t44;
t4 = -pkin(5) * t153 - t6;
t464 = t559 * t576;
t332 = -t393 * t464 + t397 * t398;
t410 = t398 * t393 + t397 * t464;
t502 = t387 * t576;
t477 = t386 * t502;
t220 = t332 * t396 + (-t389 * t410 + t477) * t392;
t289 = t386 * t410 + t389 * t502;
t160 = -t220 * t379 + t289 * t380;
t427 = g(1) * t160 + g(2) * t587 + g(3) * (-t275 * t379 - t326 * t380);
t579 = t120 * (pkin(5) * t597 + t120 * pkin(12)) + t4 + t427;
t463 = qJD(3) * t486;
t203 = t396 * t463 + (t434 * qJD(2) + t581 * qJD(3)) * t387;
t522 = qJD(2) * t387;
t499 = t393 * t522;
t474 = t386 * t499;
t178 = t203 * t388 + t385 * t474;
t202 = t387 * t400 + t392 * t463;
t269 = (t389 * t542 + t486) * pkin(10) + t412;
t279 = t374 + t407;
t366 = t397 * t471;
t285 = qJD(2) * t440 + t366;
t286 = t405 * qJD(2);
t318 = t387 * t422;
t408 = -t269 * t521 + t279 * t495 + t396 * t285 + t286 * t540 + t309 * t497 + t318 * t545;
t86 = qJ(4) * t474 - qJD(4) * t326 + t408;
t205 = -t286 * t386 + t389 * t318;
t93 = pkin(3) * t202 - qJ(4) * t203 - qJD(4) * t275 + t205;
t46 = -t385 * t86 + t388 * t93;
t34 = pkin(4) * t202 - pkin(11) * t178 + t46;
t177 = t203 * t385 - t388 * t474;
t47 = t385 * t93 + t388 * t86;
t37 = -pkin(11) * t177 + t47;
t198 = -t279 * t386 + t389 * t309;
t134 = pkin(3) * t274 - qJ(4) * t275 + t198;
t506 = t396 * t269 + t279 * t540 + t309 * t545;
t141 = -qJ(4) * t326 + t506;
t82 = t388 * t134 - t141 * t385;
t61 = pkin(4) * t274 - pkin(11) * t211 + t82;
t83 = t385 * t134 + t388 * t141;
t66 = -pkin(11) * t210 + t83;
t449 = t391 * t61 + t395 * t66;
t577 = -qJD(5) * t449 + t34 * t395 - t37 * t391;
t437 = t392 * t180 - t185 * t539 - t221 * t544 + t244 * t520 + t249 * t496 + t296 * t498;
t56 = t437 + t580;
t36 = pkin(4) * t116 + t56;
t10 = pkin(5) * t45 - pkin(12) * t44 + t36;
t5 = t391 * t15 + t395 * t17 + t49 * t518 - t519 * t57;
t3 = pkin(12) * t153 + t5;
t20 = pkin(12) * t248 + t22;
t118 = t304 * pkin(3) + qJD(4) - t138;
t89 = pkin(4) * t482 + t118;
t38 = t122 * pkin(5) - pkin(12) * t597 + t89;
t454 = t20 * t390 - t38 * t394;
t1 = -t454 * qJD(6) + t390 * t10 + t394 * t3;
t399 = qJD(1) ^ 2;
t569 = t589 * pkin(5) + t530 * qJD(5) + t594 * t391 + t595 * t395;
t566 = t120 * t96;
t565 = t120 * t98;
t516 = qJD(6) * t394;
t517 = qJD(6) * t390;
t25 = t390 * t153 + t248 * t516 + t394 * t44 - t517 * t597;
t564 = t25 * t390;
t563 = t390 * t43;
t560 = pkin(5) * t253 - t585;
t558 = qJ(4) * t156;
t557 = t139 * t304;
t554 = t274 * t390;
t552 = t343 * t390;
t551 = t343 * t394;
t550 = t379 * t386;
t549 = t380 * t386;
t548 = t380 * t390;
t547 = t380 * t394;
t546 = t381 * t399;
t543 = t387 * t393;
t213 = t247 * t390 + t394 * t544;
t534 = -qJD(6) * t213 - t589 * t390 + t532 * t394;
t508 = t390 * t544;
t533 = -qJD(6) * t508 + t247 * t516 + t532 * t390 + t589 * t394;
t383 = t393 ^ 2;
t524 = -t397 ^ 2 + t383;
t515 = -qJD(4) + t118;
t511 = t397 * t546;
t510 = t386 * t543;
t378 = -pkin(4) * t388 - pkin(3);
t480 = t120 * t394;
t101 = -pkin(4) * t555 + t139;
t478 = -t269 * t520 - t279 * t496 - t392 * t285 - t309 * t498;
t369 = pkin(10) * t545;
t323 = t369 + (-pkin(2) * t396 - pkin(3)) * t389;
t250 = pkin(4) * t325 + t323;
t140 = pkin(5) * t246 - pkin(12) * t247 + t250;
t466 = t589 * pkin(12) - qJD(6) * t140 + t583;
t462 = t387 * t399 * t559;
t406 = t410 * t396;
t219 = t332 * t392 + t389 * t406 - t396 * t477;
t461 = g(1) * t215 - g(2) * t219;
t109 = -pkin(12) * t544 + t530;
t459 = -t531 * pkin(5) + t532 * pkin(12) + qJD(6) * t109 - t529;
t262 = pkin(5) * t342 - pkin(12) * t343 + t378;
t458 = pkin(12) * t253 - qJD(6) * t262 + t586;
t457 = -t526 * pkin(5) + t527 * pkin(12) + qJD(6) * t292 + t101;
t12 = t20 * t394 + t38 * t390;
t28 = pkin(12) * t274 + t449;
t143 = t395 * t210 + t211 * t391;
t260 = t392 * t269;
t445 = t279 * t389 + t309 * t386;
t142 = pkin(3) * t326 - t396 * t445 + t260;
t99 = pkin(4) * t210 + t142;
t50 = pkin(5) * t143 - pkin(12) * t144 + t99;
t453 = t28 * t394 + t390 * t50;
t452 = -t28 * t390 + t394 * t50;
t21 = -t391 * t57 + t395 * t49;
t450 = -t391 * t66 + t395 * t61;
t103 = t144 * t394 + t554;
t447 = t182 * t395 - t196 * t391;
t442 = 0.2e1 * t484 + qJD(2);
t433 = t391 * t34 + t395 * t37 + t61 * t518 - t519 * t66;
t428 = t387 * (t479 + t368);
t426 = g(1) * t219 + g(2) * t215 + g(3) * t274;
t425 = -g(1) * t220 + g(2) * t218 - g(3) * t275;
t239 = -t330 * t392 + t331 * t539;
t241 = t332 * t539 - t392 * t410;
t424 = g(1) * t241 + g(2) * t239 + g(3) * t315;
t240 = -t330 * t396 - t331 * t540;
t242 = -t332 * t540 - t406;
t423 = g(1) * t242 + g(2) * t240 + g(3) * t316;
t128 = t165 * t390 - t394 * t253;
t421 = -t328 * t390 + t343 * t516 - t128;
t129 = t165 * t394 + t253 * t390;
t420 = -t328 * t394 - t343 * t517 - t129;
t417 = t426 - t56;
t416 = g(1) * t332 + g(2) * t331 + g(3) * t543;
t19 = -pkin(5) * t248 - t21;
t409 = -pkin(12) * t43 + (t19 + t21) * t120;
t2 = -qJD(6) * t12 + t394 * t10 - t390 * t3;
t403 = t118 * t520 - t416;
t402 = t426 - t437;
t401 = t443 * t412;
t92 = -t286 * t539 + (-pkin(3) * t499 - t318 * t396) * t386 - t478;
t76 = pkin(4) * t177 + t92;
t238 = t316 * t380 + t379 * t510;
t214 = t247 * t394 - t508;
t200 = t275 * t380 - t326 * t379;
t187 = t242 * t380 + t332 * t550;
t186 = t240 * t380 + t331 * t550;
t161 = t220 * t380 + t289 * t379;
t108 = pkin(5) * t544 - t447;
t107 = t161 * t394 + t219 * t390;
t106 = -t161 * t390 + t219 * t394;
t78 = qJD(5) * t144 + t395 * t177 + t178 * t391;
t77 = -qJD(5) * t143 - t177 * t391 + t178 * t395;
t40 = qJD(6) * t103 - t202 * t394 + t390 * t77;
t39 = t584 * qJD(6) + t202 * t390 + t394 * t77;
t27 = -pkin(5) * t274 - t450;
t18 = pkin(5) * t78 - pkin(12) * t77 + t76;
t8 = -pkin(5) * t202 - t577;
t7 = pkin(12) * t202 + t433;
t9 = [qJDD(1), g(1) * t576 - g(2) * t398, g(1) * t398 + g(2) * t576 (qJDD(1) * t383 + 0.2e1 * t393 * t493) * t381 (t393 * t512 - t514 * t524) * t599, t397 * t442 * t522 + t393 * t428, t397 * t428 - t442 * t499, t368 * t559, -qJD(2) * t401 + (-pkin(9) * t543 + t374) * t368 + (t387 * t415 + t413) * t559 + g(1) * t331 - g(2) * t332 + t582 * t596 -(-pkin(9) * t499 + t366) * t443 - t412 * t368 - t411 * t559 - g(1) * t330 + g(2) * t410 + t431 * t596, t155 * t275 + t203 * t253, -t155 * t274 - t156 * t275 - t202 * t253 - t203 * t251, -t155 * t326 - t203 * t304 + t253 * t474 + t257 * t275, t156 * t326 + t202 * t304 - t251 * t474 - t257 * t274, -t257 * t326 - t304 * t474, -t478 * t304 - t260 * t257 + t437 * t326 + t138 * t474 + t205 * t251 + t198 * t156 + t130 * t274 + t183 * t202 - g(1) * t218 - g(2) * t220 + (-(t286 * t389 + t318 * t386) * t304 + t445 * t257) * t396, t130 * t275 - t139 * t474 + t198 * t155 + t183 * t203 + t205 * t253 - t257 * t506 + t304 * t408 + t62 * t326 - t461, t46 * t251 + t82 * t156 + t23 * t274 + t69 * t202 + t92 * t482 + t142 * t116 + t56 * t210 + t118 * t177 - g(1) * (t218 * t388 + t287 * t385) - g(2) * (t220 * t388 + t289 * t385) -t47 * t251 - t83 * t156 - t24 * t274 - t70 * t202 + t92 * t193 + t142 * t117 + t56 * t211 + t118 * t178 - g(1) * (-t218 * t385 + t287 * t388) - g(2) * (-t220 * t385 + t289 * t388) -t83 * t116 - t82 * t117 - t70 * t177 - t69 * t178 - t46 * t193 - t24 * t210 - t23 * t211 - t47 * t482 + t461, t24 * t83 + t70 * t47 + t23 * t82 + t69 * t46 + t56 * t142 + t118 * t92 - g(1) * (-pkin(1) * t576 - t331 * pkin(2) + t218 * pkin(3) + pkin(9) * t541 + pkin(10) * t287 - qJ(4) * t215) - g(2) * (t398 * pkin(1) + t332 * pkin(2) + t220 * pkin(3) + pkin(9) * t502 + pkin(10) * t289 + t219 * qJ(4)) t144 * t44 + t597 * t77, -t122 * t77 - t143 * t44 - t144 * t45 - t597 * t78, t144 * t153 + t202 * t597 + t248 * t77 + t274 * t44, -t122 * t202 - t143 * t153 - t248 * t78 - t274 * t45, t153 * t274 + t202 * t248, -g(1) * t159 - g(2) * t161 + t76 * t122 + t36 * t143 + t450 * t153 + t21 * t202 + t248 * t577 + t6 * t274 + t99 * t45 + t89 * t78, g(1) * t587 - g(2) * t160 + t36 * t144 - t449 * t153 - t22 * t202 - t433 * t248 - t5 * t274 + t99 * t44 + t597 * t76 + t89 * t77, t103 * t25 + t39 * t98, -t103 * t26 + t25 * t584 - t39 * t96 - t40 * t98, t103 * t43 + t120 * t39 + t143 * t25 + t78 * t98, -t120 * t40 - t143 * t26 + t43 * t584 - t78 * t96, t120 * t78 + t143 * t43 (-qJD(6) * t453 + t18 * t394 - t390 * t7) * t120 + t452 * t43 + t2 * t143 - t454 * t78 + t8 * t96 + t27 * t26 - t4 * t584 + t19 * t40 - g(1) * t600 - g(2) * t107 -(qJD(6) * t452 + t18 * t390 + t394 * t7) * t120 - t453 * t43 - t1 * t143 - t12 * t78 + t8 * t98 + t27 * t25 + t4 * t103 + t19 * t39 + g(1) * t601 - g(2) * t106; 0, 0, 0, -t393 * t511, t524 * t546, -t397 * t462 + t492, t393 * t462 + t491, t368, pkin(1) * t393 * t546 + t588 * pkin(9) + g(1) * t410 + g(2) * t330 - g(3) * t542 + qJD(1) * t401 + t413, pkin(1) * t511 + (-pkin(9) * t501 + t365) * t484 + t365 * qJD(2) + t416 - t504, t155 * t545 + (-t303 + t497) * t253, t251 * t303 + t253 * t302 + (t155 * t396 - t156 * t392 + (-t251 * t396 - t253 * t392) * qJD(3)) * t386, t155 * t389 + t303 * t304 + (-t253 * t501 + t257 * t392 - t304 * t520) * t386, -t156 * t389 - t302 * t304 + (t251 * t501 + t257 * t396 + t304 * t521) * t386, t257 * t389 + t304 * t476 (pkin(2) * t539 - t369) * t257 - t437 * t389 - t204 * t251 - t183 * t302 + t598 * t304 + (-t138 * t501 + t183 * t521 - pkin(2) * t156 + (t304 * t317 - t130) * t396) * t386 - t423, -t525 * t257 - t62 * t389 - t204 * t253 - t183 * t303 + t609 * t304 + (t139 * t501 - pkin(2) * t155 + t130 * t392 + (-pkin(10) * t304 * t392 + t183 * t396) * qJD(3)) * t386 + t424, t323 * t116 - t118 * t235 + t222 * t156 - t69 * t302 + t56 * t325 + t562 * t251 + t528 * t243 + (t304 * t528 - t423) * t388 + (-t23 * t396 + t385 * t403 + t521 * t69) * t386, t323 * t117 - t118 * t236 - t223 * t156 + t70 * t302 + t56 * t327 - t561 * t251 + t528 * t193 + t423 * t385 + (t24 * t396 + t388 * t403 - t521 * t70) * t386, -t223 * t116 - t24 * t325 - t222 * t117 - t23 * t327 + t70 * t235 + t69 * t236 - t562 * t193 + (-t385 * t70 - t388 * t69) * t497 - t424 - t561 * t482, t24 * t223 + t23 * t222 + t56 * t323 - g(1) * (-pkin(2) * t410 + t242 * pkin(3) + t241 * qJ(4) + t332 * t573) - g(2) * (-pkin(2) * t330 + pkin(3) * t240 + qJ(4) * t239 + t331 * t573) - g(3) * (pkin(3) * t316 + qJ(4) * t315 + t387 * t438) + t561 * t70 + t562 * t69 + t528 * t118, t247 * t44 + t532 * t597, -t122 * t532 - t246 * t44 - t247 * t45 - t531 * t597, t153 * t247 + t248 * t532 - t44 * t544 - t589 * t597, t122 * t589 - t153 * t246 - t248 * t531 + t45 * t544, -t153 * t544 - t248 * t589, t447 * t153 - t6 * t544 + t250 * t45 + t36 * t246 - g(1) * t187 - g(2) * t186 - g(3) * t238 + t531 * t89 + ((-qJD(5) * t196 - t595) * t395 + (-qJD(5) * t182 - t594) * t391) * t248 - t589 * t21 + t529 * t122, -t530 * t153 + t5 * t544 + t250 * t44 + t36 * t247 - g(1) * (-t242 * t379 + t332 * t549) - g(2) * (-t240 * t379 + t331 * t549) - g(3) * (-t316 * t379 + t380 * t510) + t532 * t89 + t583 * t248 + t589 * t22 + t529 * t597, t214 * t25 + t534 * t98, -t213 * t25 - t214 * t26 - t533 * t98 - t534 * t96, t120 * t534 + t214 * t43 + t246 * t25 + t531 * t98, -t120 * t533 - t213 * t43 - t246 * t26 - t531 * t96, t120 * t531 + t246 * t43 (-t109 * t390 + t140 * t394) * t43 + t2 * t246 + t108 * t26 + t4 * t213 - g(1) * (t187 * t394 + t241 * t390) - g(2) * (t186 * t394 + t239 * t390) - g(3) * (t238 * t394 + t315 * t390) + t569 * t96 + t533 * t19 + (t390 * t466 - t394 * t459) * t120 - t531 * t454 -(t109 * t394 + t140 * t390) * t43 - t1 * t246 + t108 * t25 + t4 * t214 - g(1) * (-t187 * t390 + t241 * t394) - g(2) * (-t186 * t390 + t239 * t394) - g(3) * (-t238 * t390 + t315 * t394) + t569 * t98 + t534 * t19 + (t390 * t459 + t394 * t466) * t120 - t531 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253 * t251, -t251 ^ 2 + t253 ^ 2, -t251 * t304 + t155, -t253 * t304 - t156, t257, -t183 * t253 + t402 - t557, -t138 * t304 + t183 * t251 - t425 - t62, -t385 * t558 - pkin(3) * t116 - t139 * t243 - t69 * t253 + (t385 * t515 - t94) * t251 + (t417 - t557) * t388, -t388 * t558 - pkin(3) * t117 - t139 * t193 + t253 * t70 + (t388 * t515 + t95) * t251 - t417 * t385, t94 * t193 + t95 * t243 + (-qJ(4) * t116 - qJD(4) * t482 - t69 * t251 + t95 * t304 + t24) * t388 + (qJ(4) * t117 + qJD(4) * t193 - t251 * t70 - t23) * t385 + t425, -t118 * t139 - t69 * t94 - t70 * t95 + (-t385 * t69 + t388 * t70) * qJD(4) + t417 * pkin(3) + (-t23 * t385 + t24 * t388 + t425) * qJ(4), t343 * t44 + t527 * t597, -t122 * t527 - t342 * t44 - t343 * t45 - t526 * t597, t153 * t343 + t248 * t527 - t253 * t597, t122 * t253 - t153 * t342 - t248 * t526, -t248 * t253, -t101 * t122 + t153 * t444 - t21 * t253 + t585 * t248 + t36 * t342 + t378 * t45 + t426 * t380 + t526 * t89, -t101 * t597 - t292 * t153 + t22 * t253 + t586 * t248 + t36 * t343 + t378 * t44 - t426 * t379 + t527 * t89, t25 * t551 + t420 * t98, t128 * t98 + t129 * t96 - (-t390 * t98 - t394 * t96) * t328 + (-t564 - t26 * t394 + (t390 * t96 - t394 * t98) * qJD(6)) * t343, t120 * t420 + t25 * t342 + t43 * t551 + t526 * t98, -t120 * t421 - t26 * t342 - t43 * t552 - t526 * t96, t120 * t526 + t342 * t43 (t262 * t394 - t292 * t390) * t43 + t2 * t342 - t444 * t26 + t4 * t552 - g(1) * (-t219 * t547 + t220 * t390) - g(2) * (-t215 * t547 - t218 * t390) - g(3) * (-t274 * t547 + t275 * t390) + t560 * t96 + (t390 * t458 - t394 * t457) * t120 - t526 * t454 + t421 * t19 -(t262 * t390 + t292 * t394) * t43 - t1 * t342 - t444 * t25 + t4 * t551 - g(1) * (t219 * t548 + t220 * t394) - g(2) * (t215 * t548 - t218 * t394) - g(3) * (t274 * t548 + t275 * t394) + t560 * t98 + (t390 * t457 + t394 * t458) * t120 - t526 * t12 + t420 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193 * t251 + t116, -t251 * t482 + t117, -t193 ^ 2 - t482 ^ 2, t193 * t69 + t482 * t70 - t402 + t580, 0, 0, 0, 0, 0, t248 * t597 + t45, t44 - t607, 0, 0, 0, 0, 0, -t603 + t610, -t394 * t608 - t563 - t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t597 * t122, -t122 ^ 2 + t597 ^ 2, t44 + t607, -t483 + (-qJD(5) + t248) * t597, t153, t22 * t248 - t597 * t89 - t427 + t6, g(1) * t161 - g(2) * t159 + g(3) * t200 + t122 * t89 + t21 * t248 - t5, t480 * t98 + t564 (t25 - t566) * t394 + (-t26 - t565) * t390, t120 * t480 + t563 - t602, t603 + t610, -t120 * t597, -pkin(5) * t26 - t22 * t96 + t409 * t390 - t394 * t579 + t454 * t597, -pkin(5) * t25 + t12 * t597 - t22 * t98 + t390 * t579 + t409 * t394; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98 * t96, -t96 ^ 2 + t98 ^ 2, t25 + t566, -t26 + t565, t43, t12 * t120 - t19 * t98 - g(1) * t106 - g(2) * t601 - g(3) * (-t200 * t390 + t267) + t2, -t454 * t120 + t19 * t96 + g(1) * t107 - g(2) * t600 - g(3) * (-t200 * t394 - t554) - t1;];
tau_reg  = t9;
