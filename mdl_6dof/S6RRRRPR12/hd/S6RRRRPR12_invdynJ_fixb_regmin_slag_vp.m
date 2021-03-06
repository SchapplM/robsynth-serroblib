% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x33]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPR12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_invdynJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:44:51
% EndTime: 2019-03-09 23:45:38
% DurationCPUTime: 22.62s
% Computational Cost: add. (33905->828), mult. (96892->1181), div. (0->0), fcn. (82285->18), ass. (0->389)
t394 = sin(pkin(6));
t396 = cos(pkin(7));
t404 = cos(qJ(3));
t405 = cos(qJ(2));
t552 = t404 * t405;
t400 = sin(qJ(3));
t401 = sin(qJ(2));
t555 = t400 * t401;
t446 = -t396 * t555 + t552;
t314 = t446 * t394;
t301 = qJD(1) * t314;
t393 = sin(pkin(7));
t535 = qJD(3) * t404;
t509 = t393 * t535;
t633 = t301 - t509;
t579 = cos(pkin(6));
t515 = pkin(1) * t579;
t379 = t405 * t515;
t370 = qJD(1) * t379;
t593 = pkin(10) * t396;
t472 = t394 * (-pkin(9) - t593);
t451 = t401 * t472;
t279 = qJD(1) * t451 + t370;
t378 = t401 * t515;
t416 = t405 * t472 - t378;
t280 = t416 * qJD(1);
t594 = pkin(10) * t393;
t450 = pkin(2) * t401 - t405 * t594;
t538 = qJD(1) * t394;
t318 = t450 * t538;
t557 = t396 * t400;
t562 = t393 * t400;
t374 = pkin(10) * t562;
t556 = t396 * t404;
t606 = pkin(2) * t556 - t374;
t632 = t606 * qJD(3) - t404 * t279 - t280 * t557 - t318 * t562;
t206 = -t280 * t393 + t396 * t318;
t553 = t401 * t404;
t554 = t400 * t405;
t448 = t396 * t553 + t554;
t313 = t448 * t394;
t300 = qJD(1) * t313;
t631 = pkin(3) * t300 - pkin(11) * t301 + t206 - (pkin(3) * t400 - pkin(11) * t404) * t393 * qJD(3);
t513 = t401 * t538;
t482 = t393 * t513;
t630 = pkin(11) * t482 - t632;
t399 = sin(qJ(4));
t403 = cos(qJ(4));
t334 = t396 * t399 + t403 * t562;
t545 = -qJD(4) * t334 + t399 * t633 - t403 * t482;
t333 = -t403 * t396 + t399 * t562;
t544 = -qJD(4) * t333 - t399 * t482 - t403 * t633;
t536 = qJD(3) * t400;
t510 = t393 * t536;
t613 = t300 - t510;
t629 = t631 * t403;
t561 = t393 * t404;
t540 = pkin(2) * t557 + pkin(10) * t561;
t323 = pkin(11) * t396 + t540;
t324 = (-pkin(3) * t404 - pkin(11) * t400 - pkin(2)) * t393;
t533 = qJD(4) * t403;
t534 = qJD(4) * t399;
t628 = t323 * t534 - t324 * t533 + t399 * t631 + t630 * t403;
t512 = t405 * t538;
t479 = t396 * t512;
t350 = t404 * t479;
t496 = t579 * qJD(1);
t455 = t496 + qJD(2);
t438 = t393 * t455;
t481 = t400 * t513;
t250 = -t404 * t438 - t350 + t481;
t248 = qJD(4) + t250;
t406 = cos(qJ(1));
t499 = t406 * t579;
t597 = sin(qJ(1));
t335 = t401 * t597 - t405 * t499;
t336 = t401 * t499 + t405 * t597;
t558 = t394 * t406;
t522 = t393 * t558;
t216 = t335 * t557 - t336 * t404 + t400 * t522;
t283 = -t335 * t393 + t396 * t558;
t389 = qJ(4) + pkin(13);
t386 = sin(t389);
t387 = cos(t389);
t167 = t216 * t387 + t283 * t386;
t213 = t335 * t556 + t336 * t400 + t404 * t522;
t398 = sin(qJ(6));
t402 = cos(qJ(6));
t627 = t167 * t398 + t213 * t402;
t626 = t167 * t402 - t213 * t398;
t543 = t403 * t323 + t399 * t324;
t625 = -t613 * pkin(4) - t544 * qJ(5) - qJD(4) * t543 - qJD(5) * t334 + t399 * t630 - t629;
t624 = qJ(5) * t545 - qJD(5) * t333 - t628;
t397 = -qJ(5) - pkin(11);
t500 = qJD(4) * t397;
t559 = t394 * t405;
t422 = pkin(9) * t559 + t378;
t242 = t422 * qJD(1) + (t438 + t479) * pkin(10);
t417 = pkin(2) * t579 + t451;
t249 = qJD(2) * pkin(2) + qJD(1) * t417 + t370;
t445 = -pkin(2) * t405 - t401 * t594 - pkin(1);
t307 = t445 * t394;
t296 = qJD(1) * t307;
t146 = -t400 * t242 + (t249 * t396 + t296 * t393) * t404;
t447 = t396 * t554 + t553;
t428 = t447 * t394;
t252 = qJD(1) * t428 + t400 * t438;
t178 = pkin(3) * t252 + pkin(11) * t250;
t549 = t403 * t146 + t399 * t178;
t573 = t250 * t399;
t623 = -qJ(5) * t573 + qJD(5) * t403 + t399 * t500 - t549;
t177 = t403 * t178;
t622 = -pkin(4) * t252 - t177 + (-qJ(5) * t250 + t500) * t403 + (-qJD(5) + t146) * t399;
t571 = t283 * t399;
t621 = t216 * t403 + t571;
t620 = pkin(4) * t545;
t303 = t393 * t512 - t396 * t455 - qJD(3);
t192 = t252 * t399 + t403 * t303;
t194 = t252 * t403 - t303 * t399;
t392 = sin(pkin(13));
t395 = cos(pkin(13));
t458 = -t192 * t392 + t395 * t194;
t101 = -t402 * t248 + t398 * t458;
t490 = -t395 * t192 - t194 * t392;
t603 = qJD(6) - t490;
t619 = t101 * t603;
t615 = -t392 * t545 - t544 * t395;
t346 = t392 * t399 - t395 * t403;
t614 = t248 * t346;
t546 = t544 * t392 - t395 * t545;
t508 = t396 * t536;
t607 = pkin(2) * t508 + pkin(10) * t509 - t400 * t279 + t280 * t556;
t347 = t392 * t403 + t395 * t399;
t541 = t248 * t347;
t528 = qJDD(1) * t401;
t503 = t394 * t528;
t529 = qJD(1) * qJD(2);
t504 = t405 * t529;
t612 = -t394 * t504 - t503;
t487 = t603 * t402;
t486 = t579 * qJDD(1);
t373 = t486 + qJDD(2);
t424 = qJD(3) * t438;
t527 = qJDD(1) * t405;
t502 = t394 * t527;
t477 = t396 * t502;
t163 = qJD(3) * t350 + t373 * t562 + t400 * t477 + (-qJD(2) * t396 - qJD(3)) * t481 + (t424 - t612) * t404;
t505 = t401 * t529;
t476 = t394 * t505;
t256 = t373 * t396 + qJDD(3) + (t476 - t502) * t393;
t95 = t403 * t163 - t252 * t534 + t399 * t256 - t303 * t533;
t96 = qJD(4) * t194 + t163 * t399 - t403 * t256;
t53 = -t392 * t95 - t395 * t96;
t52 = qJDD(6) - t53;
t611 = -t398 * t52 - t603 * t487;
t409 = qJD(2) * t448 + qJD(3) * t447;
t164 = -t373 * t561 + t394 * (qJD(1) * t409 + t400 * t528) + t400 * t424 - t404 * t477;
t610 = t216 * t399 - t283 * t403;
t587 = -t624 * t392 + t395 * t625;
t586 = t392 * t625 + t624 * t395;
t581 = t392 * t623 - t395 * t622;
t580 = t392 * t622 + t395 * t623;
t278 = t379 + t417;
t200 = -t278 * t393 + t396 * t307;
t498 = t579 * t393;
t604 = t396 * t552 - t555;
t273 = -t394 * t604 - t404 * t498;
t274 = t400 * t498 + t428;
t142 = pkin(3) * t273 - pkin(11) * t274 + t200;
t330 = t393 * t559 - t396 * t579;
t267 = (t396 * t559 + t498) * pkin(10) + t422;
t519 = t404 * t267 + t278 * t557 + t307 * t562;
t150 = -pkin(11) * t330 + t519;
t548 = t399 * t142 + t403 * t150;
t210 = t274 * t399 + t330 * t403;
t211 = t274 * t403 - t330 * t399;
t152 = -t210 * t392 + t211 * t395;
t265 = t273 * t402;
t608 = -t152 * t398 + t265;
t569 = t318 * t404;
t542 = -(-pkin(3) * t513 - t569) * t393 + t607;
t147 = t404 * t242 + t249 * t557 + t296 * t562;
t605 = -t147 + (t534 + t573) * pkin(4);
t602 = pkin(4) * t96 + qJDD(5);
t471 = t579 * t597;
t337 = -t406 * t401 - t405 * t471;
t338 = -t401 * t471 + t405 * t406;
t514 = t394 * t597;
t483 = t393 * t514;
t218 = t338 * t404 + (t337 * t396 + t483) * t400;
t285 = -t337 * t393 + t396 * t514;
t172 = -t218 * t399 + t285 * t403;
t601 = -g(1) * t172 - g(2) * t610 + g(3) * t210;
t160 = qJDD(4) + t164;
t478 = qJD(2) * t515;
t453 = qJD(1) * t478;
t474 = pkin(1) * t486;
t517 = pkin(9) * t502 + t401 * t474 + t405 * t453;
t420 = -pkin(9) * t476 + t517;
t184 = (t373 * t393 + (-t505 + t527) * t396 * t394) * pkin(10) + t420;
t423 = -t401 * t453 + t405 * t474;
t440 = -t504 - t528;
t425 = t440 * pkin(9);
t188 = t373 * pkin(2) + (t440 * t593 + t425) * t394 + t423;
t432 = t450 * qJD(2);
t221 = (qJD(1) * t432 + qJDD(1) * t445) * t394;
t507 = t396 * t535;
t72 = t404 * t184 + t188 * t557 + t221 * t562 - t242 * t536 + t249 * t507 + t296 * t509;
t61 = pkin(11) * t256 + t72;
t137 = -t188 * t393 + t396 * t221;
t68 = pkin(3) * t164 - pkin(11) * t163 + t137;
t185 = -t249 * t393 + t396 * t296;
t115 = pkin(3) * t250 - pkin(11) * t252 + t185;
t120 = -pkin(11) * t303 + t147;
t78 = t115 * t399 + t120 * t403;
t20 = -qJD(4) * t78 - t399 * t61 + t403 * t68;
t14 = pkin(4) * t160 - qJ(5) * t95 - qJD(5) * t194 + t20;
t443 = -t115 * t533 + t120 * t534 - t399 * t68 - t403 * t61;
t18 = -qJ(5) * t96 - qJD(5) * t192 - t443;
t5 = t14 * t395 - t18 * t392;
t3 = -pkin(5) * t160 - t5;
t383 = pkin(4) * t392 + pkin(12);
t600 = t603 * (pkin(4) * t194 + pkin(5) * t458 - pkin(12) * t490 + qJD(6) * t383) + g(1) * (-t218 * t386 + t285 * t387) + g(2) * (t216 * t386 - t283 * t387) + g(3) * (-t274 * t386 - t330 * t387) + t3;
t449 = t400 * t184 - t188 * t556 - t221 * t561 + t242 * t535 + t249 * t508 + t296 * t510;
t595 = pkin(3) * t256;
t62 = t449 - t595;
t38 = t62 + t602;
t54 = -t392 * t96 + t395 * t95;
t12 = -pkin(5) * t53 - pkin(12) * t54 + t38;
t6 = t392 * t14 + t395 * t18;
t4 = pkin(12) * t160 + t6;
t77 = t403 * t115 - t120 * t399;
t63 = -qJ(5) * t194 + t77;
t56 = pkin(4) * t248 + t63;
t64 = -qJ(5) * t192 + t78;
t58 = t395 * t64;
t29 = t392 * t56 + t58;
t27 = pkin(12) * t248 + t29;
t119 = pkin(3) * t303 - t146;
t97 = pkin(4) * t192 + qJD(5) + t119;
t47 = -pkin(5) * t490 - pkin(12) * t458 + t97;
t461 = t27 * t398 - t402 * t47;
t1 = -t461 * qJD(6) + t398 * t12 + t402 * t4;
t560 = t394 * t401;
t426 = g(1) * t338 + g(2) * t336 + g(3) * t560;
t407 = qJD(1) ^ 2;
t388 = t394 ^ 2;
t596 = pkin(1) * t388;
t470 = qJD(3) * t498;
t204 = t404 * t470 + (t446 * qJD(2) + qJD(3) * t604) * t394;
t537 = qJD(2) * t394;
t511 = t401 * t537;
t480 = t393 * t511;
t136 = -qJD(4) * t210 + t204 * t403 + t399 * t480;
t203 = t394 * t409 + t400 * t470;
t282 = t416 * qJD(2);
t319 = t394 * t432;
t207 = -t282 * t393 + t396 * t319;
t107 = pkin(3) * t203 - pkin(11) * t204 + t207;
t371 = t405 * t478;
t281 = qJD(2) * t451 + t371;
t419 = -t267 * t536 + t278 * t507 + t404 * t281 + t282 * t557 + t307 * t509 + t319 * t562;
t99 = pkin(11) * t480 + t419;
t415 = -qJD(4) * t548 + t403 * t107 - t399 * t99;
t22 = pkin(4) * t203 - qJ(5) * t136 - qJD(5) * t211 + t415;
t135 = qJD(4) * t211 + t204 * t399 - t403 * t480;
t441 = t399 * t107 + t142 * t533 - t150 * t534 + t403 * t99;
t25 = -qJ(5) * t135 - qJD(5) * t210 + t441;
t10 = t392 * t22 + t395 * t25;
t492 = t403 * t142 - t150 * t399;
t71 = pkin(4) * t273 - qJ(5) * t211 + t492;
t76 = -qJ(5) * t210 + t548;
t37 = t392 * t71 + t395 * t76;
t588 = pkin(5) * t613 - t587;
t531 = qJD(6) * t402;
t532 = qJD(6) * t398;
t32 = t398 * t160 + t248 * t531 + t402 * t54 - t458 * t532;
t585 = t32 * t398;
t584 = t392 * t64;
t582 = pkin(5) * t252 + t581;
t578 = t101 * t458;
t103 = t248 * t398 + t402 * t458;
t577 = t103 * t458;
t575 = t192 * t248;
t574 = t194 * t248;
t572 = t273 * t398;
t570 = t285 * t399;
t568 = t347 * t398;
t567 = t347 * t402;
t566 = t386 * t393;
t565 = t387 * t398;
t564 = t387 * t402;
t563 = t388 * t407;
t245 = -t333 * t392 + t334 * t395;
t219 = t245 * t398 + t402 * t561;
t551 = -qJD(6) * t219 - t398 * t613 - t402 * t615;
t521 = t398 * t561;
t550 = -qJD(6) * t521 + t245 * t531 - t398 * t615 + t402 * t613;
t489 = -t323 * t399 + t403 * t324;
t191 = -pkin(4) * t561 - qJ(5) * t334 + t489;
t198 = -qJ(5) * t333 + t543;
t128 = t392 * t191 + t395 * t198;
t390 = t401 ^ 2;
t539 = -t405 ^ 2 + t390;
t526 = t405 * t596;
t523 = t401 * t563;
t385 = pkin(4) * t403 + pkin(3);
t506 = t397 * t399;
t497 = -t402 * t160 + t398 * t54;
t494 = -t402 * t252 + t398 * t614;
t493 = t252 * t398 + t402 * t614;
t488 = t403 * t248;
t484 = -t267 * t535 - t278 * t508 - t400 * t281 - t307 * t510;
t244 = t395 * t333 + t334 * t392;
t322 = t374 + (-pkin(2) * t404 - pkin(3)) * t396;
t421 = pkin(4) * t333 + t322;
t148 = pkin(5) * t244 - pkin(12) * t245 + t421;
t473 = pkin(12) * t613 - qJD(6) * t148 - t586;
t469 = t394 * t407 * t579;
t217 = -t337 * t556 + t338 * t400 - t404 * t483;
t468 = -g(1) * t213 + g(2) * t217;
t117 = -pkin(12) * t561 + t128;
t466 = -pkin(5) * t546 - pkin(12) * t615 + qJD(6) * t117 - t542 + t620;
t262 = pkin(5) * t346 - pkin(12) * t347 - t385;
t465 = pkin(12) * t252 - qJD(6) * t262 - t580;
t362 = t397 * t403;
t294 = -t395 * t362 + t392 * t506;
t464 = -pkin(5) * t541 - pkin(12) * t614 + qJD(6) * t294 - t605;
t9 = t22 * t395 - t25 * t392;
t16 = t27 * t402 + t398 * t47;
t35 = pkin(12) * t273 + t37;
t151 = t395 * t210 + t211 * t392;
t258 = t400 * t267;
t456 = t278 * t396 + t307 * t393;
t149 = pkin(3) * t330 - t404 * t456 + t258;
t412 = pkin(4) * t210 + t149;
t57 = pkin(5) * t151 - pkin(12) * t152 + t412;
t460 = t35 * t402 + t398 * t57;
t459 = -t35 * t398 + t402 * t57;
t28 = t395 * t56 - t584;
t36 = -t392 * t76 + t395 * t71;
t110 = t152 * t402 + t572;
t127 = t191 * t395 - t198 * t392;
t454 = 0.2e1 * t496 + qJD(2);
t452 = t402 * t52 + (t398 * t490 - t532) * t603;
t442 = -pkin(11) * t160 + t119 * t248;
t437 = t394 * (t486 + t373);
t436 = g(1) * t217 + g(2) * t213 + g(3) * t273;
t435 = -g(1) * t218 + g(2) * t216 - g(3) * t274;
t238 = -t335 * t400 + t336 * t556;
t240 = t337 * t400 + t338 * t556;
t434 = g(1) * t240 + g(2) * t238 + g(3) * t313;
t239 = -t335 * t404 - t336 * t557;
t241 = t337 * t404 - t338 * t557;
t433 = g(1) * t241 + g(2) * t239 + g(3) * t314;
t431 = t347 * t531 - t494;
t430 = -t347 * t532 - t493;
t26 = -pkin(5) * t248 - t28;
t31 = t395 * t63 - t584;
t418 = -t383 * t52 + (t26 + t31) * t603;
t2 = -qJD(6) * t16 + t402 * t12 - t398 * t4;
t414 = pkin(11) * qJD(4) * t248 - t436 + t62;
t411 = t436 - t449;
t410 = t455 * t422;
t100 = -t282 * t556 + (-pkin(3) * t511 - t319 * t404) * t393 - t484;
t408 = pkin(4) * t135 + t100;
t384 = -pkin(4) * t395 - pkin(5);
t293 = -t362 * t392 - t395 * t506;
t233 = t314 * t387 + t560 * t566;
t220 = t245 * t402 - t521;
t202 = t274 * t387 - t330 * t386;
t190 = t241 * t387 + t338 * t566;
t189 = t239 * t387 + t336 * t566;
t173 = t218 * t403 + t570;
t169 = t218 * t387 + t285 * t386;
t116 = pkin(5) * t561 - t127;
t112 = t169 * t402 + t217 * t398;
t111 = -t169 * t398 + t217 * t402;
t89 = -t135 * t392 + t136 * t395;
t88 = t395 * t135 + t136 * t392;
t49 = qJD(6) * t110 - t203 * t402 + t398 * t89;
t48 = qJD(6) * t608 + t203 * t398 + t402 * t89;
t34 = -pkin(5) * t273 - t36;
t33 = qJD(6) * t103 + t497;
t30 = t392 * t63 + t58;
t23 = pkin(5) * t88 - pkin(12) * t89 + t408;
t8 = pkin(12) * t203 + t10;
t7 = -pkin(5) * t203 - t9;
t11 = [qJDD(1), g(1) * t597 - g(2) * t406, g(1) * t406 + g(2) * t597 (qJDD(1) * t390 + 0.2e1 * t401 * t504) * t388, 0.2e1 * (t401 * t527 - t529 * t539) * t388, t405 * t454 * t537 + t401 * t437, t405 * t437 - t454 * t511, t373 * t579, 0.2e1 * qJDD(1) * t526 - 0.2e1 * t505 * t596 - qJD(2) * t410 + (-pkin(9) * t560 + t379) * t373 + (t394 * t425 + t423) * t579 + g(1) * t336 - g(2) * t338 -(-pkin(9) * t511 + t371) * t455 - t422 * t373 - t420 * t579 - g(1) * t335 - g(2) * t337 + 0.2e1 * t440 * t596, t163 * t274 + t204 * t252, -t163 * t273 - t164 * t274 - t203 * t252 - t204 * t250, -t163 * t330 - t204 * t303 + t252 * t480 + t256 * t274, t164 * t330 + t203 * t303 - t250 * t480 - t256 * t273, -t256 * t330 - t303 * t480, -t484 * t303 - t258 * t256 + t449 * t330 + t146 * t480 + t207 * t250 + t200 * t164 + t137 * t273 + t185 * t203 - g(1) * t216 - g(2) * t218 + (-(t282 * t396 + t319 * t393) * t303 + t456 * t256) * t404, t137 * t274 - t147 * t480 + t200 * t163 + t185 * t204 + t207 * t252 - t256 * t519 + t303 * t419 + t72 * t330 + t468, t136 * t194 + t211 * t95, -t135 * t194 - t136 * t192 - t210 * t95 - t211 * t96, t136 * t248 + t160 * t211 + t194 * t203 + t273 * t95, -t135 * t248 - t160 * t210 - t192 * t203 - t273 * t96, t160 * t273 + t203 * t248, -g(1) * t621 - g(2) * t173 + t100 * t192 + t119 * t135 + t149 * t96 + t492 * t160 + t20 * t273 + t77 * t203 + t62 * t210 + t415 * t248, g(1) * t610 - g(2) * t172 + t100 * t194 + t119 * t136 + t149 * t95 - t548 * t160 - t78 * t203 + t62 * t211 - t441 * t248 + t443 * t273, t10 * t490 - t151 * t6 - t152 * t5 - t28 * t89 - t29 * t88 - t36 * t54 + t37 * t53 - t458 * t9 - t468, t6 * t37 + t29 * t10 + t5 * t36 + t28 * t9 + t38 * t412 + t97 * t408 - g(1) * (-pkin(1) * t597 - t336 * pkin(2) + pkin(4) * t571 + pkin(9) * t558 + pkin(10) * t283 + t213 * t397 + t216 * t385) - g(2) * (t406 * pkin(1) + t338 * pkin(2) + pkin(4) * t570 + pkin(9) * t514 + pkin(10) * t285 - t217 * t397 + t218 * t385) t103 * t48 + t110 * t32, -t101 * t48 - t103 * t49 - t110 * t33 + t32 * t608, t103 * t88 + t110 * t52 + t151 * t32 + t48 * t603, -t101 * t88 - t151 * t33 - t49 * t603 + t52 * t608, t151 * t52 + t603 * t88 (-qJD(6) * t460 + t23 * t402 - t398 * t8) * t603 + t459 * t52 + t2 * t151 - t461 * t88 + t7 * t101 + t34 * t33 - t3 * t608 + t26 * t49 - g(1) * t626 - g(2) * t112 -(qJD(6) * t459 + t23 * t398 + t402 * t8) * t603 - t460 * t52 - t1 * t151 - t16 * t88 + t7 * t103 + t34 * t32 + t3 * t110 + t26 * t48 + g(1) * t627 - g(2) * t111; 0, 0, 0, -t405 * t523, t539 * t563, -t405 * t469 + t503, t401 * t469 + t502, t373, pkin(1) * t523 + pkin(9) * t612 - g(1) * t337 + g(2) * t335 - g(3) * t559 + qJD(1) * t410 + t423, t407 * t526 + (-pkin(9) * t513 + t370) * t496 + t370 * qJD(2) + t426 - t517, t163 * t562 - t252 * t633, t250 * t301 + t252 * t300 + (t163 * t404 - t164 * t400 + (-t250 * t404 - t252 * t400) * qJD(3)) * t393, t163 * t396 + t301 * t303 + (-t252 * t513 + t256 * t400 - t303 * t535) * t393, -t164 * t396 - t300 * t303 + (t250 * t513 + t256 * t404 + t303 * t536) * t393, t256 * t396 + t303 * t482, t606 * t256 - t449 * t396 - t206 * t250 - t185 * t300 + t607 * t303 + (-t146 * t513 + t185 * t536 - pkin(2) * t164 + (t303 * t318 - t137) * t404) * t393 - t433, -t540 * t256 - t72 * t396 - t206 * t252 - t185 * t301 + t632 * t303 + (-pkin(2) * t163 + t137 * t400 + t147 * t513 + t185 * t535) * t393 + t434, t194 * t544 + t334 * t95, -t192 * t544 + t194 * t545 - t333 * t95 - t334 * t96, t160 * t334 - t194 * t613 + t248 * t544 - t561 * t95, -t160 * t333 + t192 * t613 + t248 * t545 + t561 * t96, -t160 * t561 - t248 * t613, t489 * t160 + t322 * t96 + t62 * t333 - t77 * t300 - t433 * t403 + (-t20 * t404 - t399 * t426 + t536 * t77) * t393 + (-t323 * t533 + (-qJD(4) * t324 + t630) * t399 - t629) * t248 + t542 * t192 - t545 * t119, -t543 * t160 + t322 * t95 + t62 * t334 + t78 * t300 + t433 * t399 + (-t403 * t426 - t404 * t443 - t536 * t78) * t393 + t628 * t248 + t542 * t194 + t544 * t119, -t127 * t54 + t128 * t53 - t244 * t6 - t245 * t5 + t28 * t615 - t546 * t29 - t587 * t458 + t586 * t490 - t434, t6 * t128 + t5 * t127 + t38 * t421 - g(1) * (pkin(2) * t337 - t240 * t397 + t241 * t385) - g(2) * (-pkin(2) * t335 - t238 * t397 + t239 * t385) - g(3) * (pkin(2) * t559 - t313 * t397 + t314 * t385) + (t607 - t620) * t97 + t586 * t29 + t587 * t28 + ((pkin(3) * qJD(1) * t560 + t569) * t97 - t426 * (pkin(4) * t399 + pkin(10))) * t393, t103 * t551 + t220 * t32, -t101 * t551 - t103 * t550 - t219 * t32 - t220 * t33, t103 * t546 + t220 * t52 + t244 * t32 + t551 * t603, -t101 * t546 - t219 * t52 - t244 * t33 - t550 * t603, t244 * t52 + t546 * t603 (-t117 * t398 + t148 * t402) * t52 + t2 * t244 + t116 * t33 + t3 * t219 - g(1) * (t190 * t402 + t240 * t398) - g(2) * (t189 * t402 + t238 * t398) - g(3) * (t233 * t402 + t313 * t398) + t550 * t26 - t546 * t461 + (t398 * t473 - t402 * t466) * t603 + t588 * t101 -(t117 * t402 + t148 * t398) * t52 - t1 * t244 + t116 * t32 + t3 * t220 - g(1) * (-t190 * t398 + t240 * t402) - g(2) * (-t189 * t398 + t238 * t402) - g(3) * (-t233 * t398 + t313 * t402) + t551 * t26 - t546 * t16 + (t398 * t466 + t402 * t473) * t603 + t588 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252 * t250, -t250 ^ 2 + t252 ^ 2, -t250 * t303 + t163, -t252 * t303 - t164, t256, -t147 * t303 - t185 * t252 + t411, -t146 * t303 + t185 * t250 - t435 - t72, t194 * t488 + t399 * t95 (t95 - t575) * t403 + (-t96 - t574) * t399, t160 * t399 - t194 * t252 + t248 * t488, -t248 ^ 2 * t399 + t160 * t403 + t192 * t252, -t248 * t252, -pkin(3) * t96 - t147 * t192 - t177 * t248 - t77 * t252 + (t146 * t248 + t442) * t399 - t414 * t403, -pkin(3) * t95 - t147 * t194 + t248 * t549 + t78 * t252 + t399 * t414 + t403 * t442, t28 * t614 - t541 * t29 + t293 * t54 + t294 * t53 - t346 * t6 - t347 * t5 + t581 * t458 + t580 * t490 + t435, t6 * t294 - t5 * t293 - t38 * t385 - g(1) * (-t217 * t385 - t218 * t397) - g(2) * (-t213 * t385 + t216 * t397) - g(3) * (-t273 * t385 - t274 * t397) + t605 * t97 + t580 * t29 - t581 * t28, t103 * t430 + t32 * t567, t494 * t103 + t493 * t101 + (-t585 - t33 * t402 + (t101 * t398 - t103 * t402) * qJD(6)) * t347, t103 * t541 + t32 * t346 + t430 * t603 + t52 * t567, -t101 * t541 - t33 * t346 - t431 * t603 - t52 * t568, t346 * t52 + t541 * t603 (t262 * t402 - t294 * t398) * t52 + t2 * t346 + t293 * t33 + t3 * t568 - g(1) * (-t217 * t564 + t218 * t398) - g(2) * (-t213 * t564 - t216 * t398) - g(3) * (-t273 * t564 + t274 * t398) - t541 * t461 + (t398 * t465 - t402 * t464) * t603 + t582 * t101 + t431 * t26 -(t262 * t398 + t294 * t402) * t52 - t1 * t346 + t293 * t32 + t3 * t567 - g(1) * (t217 * t565 + t218 * t402) - g(2) * (t213 * t565 - t216 * t402) - g(3) * (t273 * t565 + t274 * t402) - t541 * t16 + (t398 * t464 + t402 * t465) * t603 + t582 * t103 + t430 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194 * t192, -t192 ^ 2 + t194 ^ 2, t95 + t575, t574 - t96, t160, -t119 * t194 + t78 * t248 + t20 + t601, g(1) * t173 - g(2) * t621 + g(3) * t211 + t119 * t192 + t77 * t248 + t443 (t392 * t53 - t395 * t54) * pkin(4) + (-t31 + t28) * t490 + (t29 - t30) * t458, t28 * t30 - t29 * t31 + (-t97 * t194 + t6 * t392 + t5 * t395 + t601) * pkin(4), t103 * t487 + t585 (t32 - t619) * t402 + (-t103 * t603 - t33) * t398, -t577 - t611, t452 + t578, -t603 * t458, -t30 * t101 + t384 * t33 + t418 * t398 - t402 * t600 + t458 * t461, -t30 * t103 + t16 * t458 + t384 * t32 + t398 * t600 + t418 * t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t458 ^ 2 - t490 ^ 2, t28 * t458 - t29 * t490 - t411 - t595 + t602, 0, 0, 0, 0, 0, t452 - t578, -t577 + t611; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103 * t101, -t101 ^ 2 + t103 ^ 2, t32 + t619, -t497 + (-qJD(6) + t603) * t103, t52, t16 * t603 - t26 * t103 - g(1) * t111 - g(2) * t627 - g(3) * (-t202 * t398 + t265) + t2, -t461 * t603 + t26 * t101 + g(1) * t112 - g(2) * t626 - g(3) * (-t202 * t402 - t572) - t1;];
tau_reg  = t11;
