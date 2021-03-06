% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRRP12
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
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRP12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:20:44
% EndTime: 2019-03-10 03:21:50
% DurationCPUTime: 30.74s
% Computational Cost: add. (37263->967), mult. (106118->1314), div. (0->0), fcn. (89759->14), ass. (0->424)
t385 = cos(qJ(2));
t612 = cos(pkin(6));
t557 = pkin(1) * t612;
t371 = t385 * t557;
t363 = qJD(1) * t371;
t382 = sin(qJ(2));
t378 = sin(pkin(6));
t611 = cos(pkin(7));
t470 = t378 * (-pkin(10) * t611 - pkin(9));
t448 = t382 * t470;
t274 = qJD(1) * t448 + t363;
t370 = t382 * t557;
t413 = t385 * t470 - t370;
t275 = t413 * qJD(1);
t610 = sin(pkin(7));
t535 = t385 * t610;
t445 = pkin(2) * t382 - pkin(10) * t535;
t578 = qJD(1) * t378;
t314 = t445 * t578;
t641 = cos(qJ(3));
t513 = t611 * t641;
t381 = sin(qJ(3));
t539 = t381 * t610;
t457 = pkin(2) * t513 - pkin(10) * t539;
t540 = t381 * t611;
t675 = t457 * qJD(3) - t641 * t274 - t275 * t540 - t314 * t539;
t433 = -t382 * t540 + t385 * t641;
t310 = t433 * t378;
t296 = qJD(1) * t310;
t512 = t610 * t641;
t472 = qJD(3) * t512;
t444 = t472 - t296;
t195 = -t275 * t610 + t611 * t314;
t431 = t381 * t385 + t382 * t513;
t309 = t431 * t378;
t295 = qJD(1) * t309;
t674 = t295 * pkin(3) - t296 * pkin(11) + t195 - (pkin(3) * t539 - pkin(11) * t512) * qJD(3);
t542 = t378 * t610;
t511 = qJD(1) * t542;
t479 = t382 * t511;
t673 = pkin(11) * t479 - t675;
t421 = t431 * qJD(2);
t536 = t385 * t611;
t432 = t381 * t536 + t382 * t641;
t395 = qJD(3) * t432 + t421;
t565 = qJDD(1) * t382;
t546 = t381 * t565;
t392 = (qJD(1) * t395 + t546) * t378;
t533 = t612 * qJD(1);
t487 = t533 + qJD(2);
t438 = t610 * t487;
t428 = t381 * t438;
t523 = t612 * qJDD(1);
t477 = t523 + qJDD(2);
t437 = t477 * t610;
t481 = t385 * t513;
t462 = t378 * t481;
t583 = -qJDD(1) * t462 - t641 * t437;
t411 = qJD(3) * t428 + t583;
t149 = t392 + t411;
t672 = -qJDD(4) - t149;
t380 = sin(qJ(4));
t384 = cos(qJ(4));
t325 = t380 * t539 - t384 * t611;
t656 = -qJD(4) * t325 - t380 * t479 + t384 * t444;
t537 = t384 * t610;
t326 = t380 * t611 + t381 * t537;
t586 = qJD(4) * t326 + t380 * t444 + t384 * t479;
t581 = pkin(2) * t540 + pkin(10) * t512;
t671 = t581 * qJD(3) - t381 * t274 + t275 * t513 + t314 * t512;
t508 = qJD(3) * t539;
t468 = t508 - t295;
t655 = pkin(3) * t479 + t671;
t318 = pkin(11) * t611 + t581;
t555 = t610 * pkin(2);
t319 = -pkin(3) * t512 - pkin(11) * t539 - t555;
t572 = qJD(4) * t384;
t574 = qJD(4) * t380;
t670 = t318 * t574 - t319 * t572 + t380 * t674 + t673 * t384;
t554 = t382 * t578;
t521 = t381 * t554;
t582 = -qJD(1) * t462 - t641 * t438;
t242 = t521 + t582;
t442 = qJD(4) + t242;
t669 = -pkin(12) * t468 + t670;
t668 = t586 * pkin(4) - pkin(12) * t656 + t655;
t386 = cos(qJ(1));
t534 = t386 * t612;
t640 = sin(qJ(1));
t327 = t382 * t640 - t385 * t534;
t328 = t382 * t534 + t385 * t640;
t598 = t378 * t386;
t214 = -t327 * t540 + t328 * t641 - t539 * t598;
t543 = t378 * t611;
t278 = t327 * t610 - t386 * t543;
t156 = t214 * t384 + t278 * t380;
t213 = t327 * t513 + t328 * t381 + t512 * t598;
t379 = sin(qJ(5));
t383 = cos(qJ(5));
t99 = t156 * t379 - t213 * t383;
t100 = t156 * t383 + t213 * t379;
t599 = t378 * t385;
t439 = pkin(9) * t599 + t370;
t515 = t378 * t536;
t238 = t439 * qJD(1) + (qJD(1) * t515 + t438) * pkin(10);
t417 = pkin(2) * t612 + t448;
t241 = qJD(2) * pkin(2) + qJD(1) * t417 + t363;
t538 = t382 * t610;
t446 = pkin(2) * t385 + pkin(10) * t538;
t443 = -pkin(1) - t446;
t303 = t443 * t378;
t288 = qJD(1) * t303;
t124 = t641 * t238 + t241 * t540 + t288 * t539;
t664 = -t124 + t442 * (pkin(4) * t380 - pkin(12) * t384);
t281 = t379 * t326 + t383 * t512;
t589 = qJD(5) * t281 - t468 * t379 - t383 * t656;
t482 = t379 * t512;
t570 = qJD(5) * t383;
t588 = -qJD(5) * t482 + t326 * t570 + t379 * t656 - t468 * t383;
t667 = -t318 * t572 - t319 * t574 + t673 * t380 - t384 * t674;
t665 = pkin(11) * t574;
t618 = -pkin(4) * t468 - t667;
t423 = t432 * t378;
t244 = qJD(1) * t423 + t428;
t567 = t385 * t511 - qJD(3);
t419 = -t487 * t611 + t567;
t181 = t380 * t244 + t384 * t419;
t516 = t378 * t535;
t416 = -qJDD(1) * t516 + t477 * t611 + qJDD(3);
t517 = t378 * t538;
t478 = qJD(2) * t517;
t398 = qJD(1) * t478 + t416;
t525 = qJDD(1) * t611;
t502 = t385 * t525;
t532 = t611 * qJD(2);
t548 = t378 * t565;
t566 = qJD(1) * qJD(2);
t549 = t385 * t566;
t659 = t378 * t549 + t548;
t148 = -t532 * t521 + t659 * t641 + (t378 * t502 + t437) * t381 - t242 * qJD(3);
t594 = t384 * t148;
t394 = t380 * t398 + t594;
t391 = -qJD(4) * t181 + t394;
t661 = qJD(5) * t442 + t391;
t317 = -pkin(3) * t611 - t457;
t209 = t325 * pkin(4) - t326 * pkin(12) + t317;
t585 = t384 * t318 + t380 * t319;
t212 = -pkin(12) * t512 + t585;
t571 = qJD(5) * t379;
t660 = -t209 * t570 + t212 * t571 - t668 * t379 + t383 * t669;
t606 = t242 * t380;
t658 = t606 + t574;
t180 = qJD(5) + t181;
t157 = -t214 * t380 + t278 * t384;
t375 = t378 ^ 2;
t639 = pkin(1) * t375;
t657 = 0.2e1 * t639;
t587 = t379 * t209 + t383 * t212;
t630 = t586 * pkin(5) - qJD(5) * t587 + t669 * t379 + t668 * t383;
t629 = qJ(6) * t586 + qJD(6) * t325 - t660;
t282 = t383 * t326 - t482;
t628 = pkin(5) * t588 + qJ(6) * t589 - t282 * qJD(6) + t618;
t405 = -qJDD(4) - t411;
t654 = t392 - t405;
t514 = t612 * t640;
t434 = t386 * t382 + t385 * t514;
t653 = t434 * t611 - t640 * t542;
t652 = (qJDD(2) + 0.2e1 * t523) * t378;
t485 = pkin(4) * t384 + pkin(12) * t380 + pkin(3);
t123 = -t381 * t238 + t241 * t513 + t288 * t512;
t165 = pkin(3) * t244 + pkin(11) * t242;
t592 = t384 * t123 + t380 * t165;
t74 = pkin(12) * t244 + t592;
t651 = -t379 * t664 + t383 * t74 + t485 * t570;
t183 = t384 * t244 - t380 * t419;
t527 = t380 * t148 - t384 * t398;
t81 = t183 * qJD(4) + t527;
t80 = qJDD(5) + t81;
t649 = qJD(3) * t244 + t378 * (qJD(1) * t421 + t546) + qJDD(4) + t583;
t329 = -t382 * t514 + t385 * t386;
t217 = t329 * t381 + t641 * t653;
t498 = t612 * t610;
t447 = t641 * t498;
t596 = t381 * t382;
t268 = t378 * t596 - t447 - t462;
t449 = g(1) * t217 + g(2) * t213 + g(3) * t268;
t218 = t329 * t641 - t381 * t653;
t399 = -t434 * t610 - t543 * t640;
t159 = t218 * t380 + t384 * t399;
t471 = t381 * t498;
t269 = t471 + t423;
t499 = t612 * t611;
t420 = t516 - t499;
t207 = t269 * t380 + t384 * t420;
t452 = g(1) * t159 - g(2) * t157 + g(3) * t207;
t135 = t183 * t379 - t383 * t442;
t137 = t383 * t183 + t379 * t442;
t106 = -pkin(11) * t419 + t124;
t176 = -t241 * t610 + t611 * t288;
t96 = t242 * pkin(3) - t244 * pkin(11) + t176;
t58 = -t380 * t106 + t384 * t96;
t51 = -pkin(4) * t442 - t58;
t26 = t135 * pkin(5) - t137 * qJ(6) + t51;
t642 = pkin(12) * t80;
t648 = t180 * t26 - t642;
t193 = qJD(3) * t471 + t378 * t395;
t273 = t371 + t417;
t191 = -t273 * t610 + t611 * t303;
t631 = t269 * pkin(11);
t119 = t268 * pkin(3) + t191 - t631;
t260 = (t515 + t498) * pkin(10) + t439;
t560 = t641 * t260 + t273 * t540 + t303 * t539;
t126 = -pkin(11) * t420 + t560;
t520 = qJD(2) * t557;
t364 = t385 * t520;
t276 = qJD(2) * t448 + t364;
t277 = t413 * qJD(2);
t435 = qJD(2) * t445;
t315 = t378 * t435;
t473 = qJD(3) * t513;
t576 = qJD(3) * t381;
t429 = -t260 * t576 + t273 * t473 + t641 * t276 + t277 * t540 + t303 * t472 + t315 * t539;
t84 = pkin(11) * t478 + t429;
t194 = qJD(3) * t447 + ((t481 - t596) * qJD(3) + t433 * qJD(2)) * t378;
t196 = -t277 * t610 + t611 * t315;
t91 = t193 * pkin(3) - t194 * pkin(11) + t196;
t463 = t119 * t572 - t126 * t574 + t380 * t91 + t384 * t84;
t21 = pkin(12) * t193 + t463;
t208 = t269 * t384 - t380 * t420;
t111 = qJD(4) * t208 + t194 * t380 - t384 * t478;
t112 = -qJD(4) * t207 + t194 * t384 + t380 * t478;
t531 = t611 * qJD(3);
t510 = t381 * t531;
t550 = t641 * qJD(3);
t401 = -t260 * t550 - t273 * t510 - t381 * t276 + t277 * t513 - t303 * t508 + t315 * t512;
t85 = -pkin(3) * t478 - t401;
t36 = t111 * pkin(4) - t112 * pkin(12) + t85;
t591 = t380 * t119 + t384 * t126;
t67 = pkin(12) * t268 + t591;
t409 = -t381 * t260 + t273 * t513 + t303 * t512;
t125 = pkin(3) * t420 - t409;
t78 = t207 * pkin(4) - t208 * pkin(12) + t125;
t492 = t379 * t78 + t383 * t67;
t647 = -qJD(5) * t492 - t21 * t379 + t36 * t383;
t522 = pkin(9) * t554;
t483 = qJD(1) * t520;
t518 = pkin(1) * t523;
t564 = qJDD(1) * t385;
t547 = t378 * t564;
t558 = pkin(9) * t547 + t382 * t518 + t385 * t483;
t436 = -qJD(2) * t522 + t558;
t503 = qJD(1) * t532;
t171 = (t437 + (-t382 * t503 + t502) * t378) * pkin(10) + t436;
t440 = -t382 * t483 + t385 * t518;
t461 = -t549 - t565;
t441 = t461 * pkin(9);
t177 = t477 * pkin(2) + ((-t382 * t525 - t385 * t503) * pkin(10) + t441) * t378 + t440;
t220 = (qJD(1) * t435 + qJDD(1) * t443) * t378;
t469 = t381 * t171 - t177 * t513 - t220 * t512 + t238 * t550 + t241 * t510 + t288 * t508;
t46 = -pkin(3) * t398 + t469;
t646 = t46 - t449;
t645 = t137 ^ 2;
t644 = t180 ^ 2;
t387 = qJD(1) ^ 2;
t643 = pkin(5) * t80;
t638 = pkin(11) * t384;
t59 = t384 * t106 + t380 * t96;
t625 = pkin(12) * qJD(5);
t624 = qJ(6) * t80;
t52 = pkin(12) * t442 + t59;
t105 = pkin(3) * t419 - t123;
t64 = t181 * pkin(4) - t183 * pkin(12) + t105;
t24 = t379 * t64 + t383 * t52;
t20 = qJ(6) * t180 + t24;
t623 = t180 * t20;
t622 = t180 * t24;
t38 = t183 * t571 + t379 * t672 - t383 * t661;
t621 = t379 * t38;
t620 = t379 * t80;
t619 = t383 * t80;
t569 = qJD(5) * t384;
t573 = qJD(4) * t383;
t617 = -qJD(6) * t384 + (-t379 * t569 - t380 * t573) * pkin(11) - t651 + t658 * qJ(6);
t616 = -t485 * t571 + (-t74 - t665) * t379 + (pkin(11) * t569 - t664) * t383 - t658 * pkin(5);
t597 = t379 * t384;
t161 = -t242 * t597 - t383 * t244;
t595 = t383 * t384;
t162 = -t242 * t595 + t244 * t379;
t495 = pkin(5) * t379 - qJ(6) * t383;
t475 = pkin(11) + t495;
t496 = pkin(5) * t383 + qJ(6) * t379;
t530 = -t380 * t123 + t165 * t384;
t73 = -pkin(4) * t244 - t530;
t615 = -pkin(5) * t161 + qJ(6) * t162 + (qJD(5) * t496 - qJD(6) * t383) * t380 + t475 * t572 - t73;
t614 = -qJD(6) * t379 + t180 * t495 - t59;
t108 = pkin(4) * t183 + pkin(12) * t181;
t613 = t379 * t108 + t383 * t58;
t609 = t135 * t180;
t608 = t137 * t135;
t607 = t137 * t180;
t605 = t242 * t384;
t602 = t485 * t383;
t601 = t375 * t387;
t600 = t378 * t382;
t23 = -t379 * t52 + t383 * t64;
t593 = qJD(6) - t23;
t580 = pkin(11) * t595 - t379 * t485;
t376 = t382 ^ 2;
t579 = -t385 ^ 2 + t376;
t577 = qJD(2) * t378;
t575 = qJD(4) * t379;
t568 = t105 * qJD(4);
t562 = t382 * t639;
t561 = t385 * t601;
t556 = pkin(10) * t610;
t553 = t382 * t577;
t552 = t180 * t571;
t551 = qJD(4) + t582;
t430 = -t641 * t171 - t177 * t540 - t220 * t539 + t238 * t576 - t241 * t473 - t288 * t472;
t45 = pkin(11) * t398 - t430;
t113 = -t177 * t610 + t611 * t220;
t53 = t149 * pkin(3) - t148 * pkin(11) + t113;
t464 = t106 * t574 - t380 * t53 - t384 * t45 - t96 * t572;
t11 = pkin(12) * t654 - t464;
t17 = t81 * pkin(4) - pkin(12) * t391 + t46;
t544 = t379 * t11 - t383 * t17 + t52 * t570 + t64 * t571;
t541 = t380 * t610;
t14 = -t106 * t572 - t380 * t45 + t384 * t53 - t96 * t574;
t529 = t119 * t384 - t380 * t126;
t526 = -t380 * t318 + t319 * t384;
t524 = t180 * t383;
t509 = t378 * t387 * t612;
t160 = t218 * t384 - t380 * t399;
t103 = t160 * t379 - t217 * t383;
t506 = -g(1) * t99 + g(2) * t103;
t104 = t160 * t383 + t217 * t379;
t505 = g(1) * t100 - g(2) * t104;
t504 = g(1) * t157 + g(2) * t159;
t211 = pkin(4) * t512 - t526;
t501 = t379 * t572 - t161;
t500 = t383 * t572 - t162;
t153 = t208 * t379 - t268 * t383;
t154 = t208 * t383 + t268 * t379;
t497 = -pkin(5) * t153 + qJ(6) * t154;
t19 = -pkin(5) * t180 + t593;
t494 = t19 * t383 - t20 * t379;
t491 = -t379 * t67 + t383 * t78;
t490 = t209 * t383 - t212 * t379;
t486 = 0.2e1 * t533 + qJD(2);
t484 = -t119 * t574 - t126 * t572 - t380 * t84 + t384 * t91;
t476 = pkin(4) + t496;
t66 = -pkin(4) * t268 - t529;
t467 = -t180 * t570 - t620;
t3 = t383 * t11 + t379 * t17 - t52 * t571 + t64 * t570;
t466 = t383 * t21 + t379 * t36 + t78 * t570 - t571 * t67;
t465 = t180 * t51 - t642;
t458 = t385 * (t531 + qJD(2));
t127 = -t213 * t597 - t214 * t383;
t129 = -t217 * t597 - t218 * t383;
t178 = -t268 * t597 - t269 * t383;
t456 = g(1) * t129 + g(2) * t127 + g(3) * t178;
t128 = -t213 * t595 + t214 * t379;
t130 = -t217 * t595 + t218 * t379;
t179 = -t268 * t595 + t269 * t379;
t455 = -g(1) * t130 - g(2) * t128 - g(3) * t179;
t235 = -t327 * t641 - t328 * t540;
t185 = t235 * t384 + t328 * t541;
t234 = -t327 * t381 + t328 * t513;
t131 = t185 * t379 - t234 * t383;
t237 = -t329 * t540 - t434 * t641;
t187 = t237 * t384 + t329 * t541;
t236 = t329 * t513 - t381 * t434;
t133 = t187 * t379 - t236 * t383;
t246 = t310 * t384 + t380 * t517;
t188 = t246 * t379 - t309 * t383;
t454 = -g(1) * t133 - g(2) * t131 - g(3) * t188;
t132 = t185 * t383 + t234 * t379;
t134 = t187 * t383 + t236 * t379;
t189 = t246 * t383 + t309 * t379;
t453 = -g(1) * t134 - g(2) * t132 - g(3) * t189;
t451 = -g(1) * t160 - g(2) * t156 - g(3) * t208;
t184 = t235 * t380 - t328 * t537;
t186 = t237 * t380 - t329 * t537;
t245 = t310 * t380 - t384 * t517;
t450 = -g(1) * t186 - g(2) * t184 - g(3) * t245;
t22 = -pkin(4) * t193 - t484;
t427 = qJD(4) * t442;
t425 = g(1) * t103 + g(2) * t99 + g(3) * t153 - t544;
t424 = qJD(2) * t513 + t550;
t422 = t180 * t625 - t452;
t12 = -pkin(4) * t654 - t14;
t39 = t183 * t570 + t379 * t661 + t383 * t672;
t5 = t39 * pkin(5) + t38 * qJ(6) - t137 * qJD(6) + t12;
t418 = -t422 - t5;
t415 = t382 * t424;
t412 = qJD(4) * t419;
t410 = -g(1) * t104 - g(2) * t100 - g(3) * t154 + t3;
t407 = t137 * t26 + qJDD(6) - t425;
t406 = t419 * t610;
t404 = t487 * t439;
t403 = qJD(3) * t406;
t397 = t242 * t442 + t427;
t396 = t398 * t610;
t320 = t475 * t380;
t291 = t602 + (pkin(11) * t379 + pkin(5)) * t384;
t290 = -qJ(6) * t384 + t580;
t139 = pkin(5) * t281 - qJ(6) * t282 + t211;
t110 = -pkin(5) * t325 - t490;
t109 = qJ(6) * t325 + t587;
t77 = pkin(5) * t137 + qJ(6) * t135;
t61 = -qJD(5) * t153 + t112 * t383 + t193 * t379;
t60 = qJD(5) * t154 + t112 * t379 - t193 * t383;
t37 = -t497 + t66;
t30 = -pkin(5) * t183 - t108 * t383 + t379 * t58;
t29 = qJ(6) * t183 + t613;
t28 = -pkin(5) * t207 - t491;
t27 = qJ(6) * t207 + t492;
t25 = -t38 + t609;
t8 = pkin(5) * t60 - qJ(6) * t61 - qJD(6) * t154 + t22;
t7 = -pkin(5) * t111 - t647;
t6 = qJ(6) * t111 + qJD(6) * t207 + t466;
t2 = qJDD(6) + t544 - t643;
t1 = qJD(6) * t180 + t3 + t624;
t4 = [qJDD(1), g(1) * t640 - g(2) * t386, g(1) * t386 + g(2) * t640 (qJDD(1) * t376 + 0.2e1 * t382 * t549) * t375, 0.2e1 * (t382 * t564 - t566 * t579) * t375, t385 * t486 * t577 + t382 * t652, t385 * t652 - t486 * t553, t477 * t612, t564 * t657 - 0.2e1 * t562 * t566 - qJD(2) * t404 + (-pkin(9) * t600 + t371) * t477 + (t378 * t441 + t440) * t612 + g(1) * t328 - g(2) * t329 -(-pkin(9) * t553 + t364) * t487 - t439 * t477 - t436 * t612 - g(1) * t327 + g(2) * t434 + t461 * t657, t148 * t269 + t194 * t244, -t148 * t268 - t149 * t269 - t193 * t244 - t194 * t242, -t148 * t420 - t194 * t419 + t244 * t478 + t269 * t398, t149 * t420 + t193 * t419 - t242 * t478 - t268 * t398, -t398 * t420 - t406 * t553, g(1) * t214 - g(2) * t218 + t113 * t268 + t123 * t478 + t191 * t149 + t176 * t193 + t196 * t242 + t398 * t409 - t401 * t419 + t420 * t469, -g(1) * t213 + g(2) * t217 + t113 * t269 - t124 * t478 + t191 * t148 + t176 * t194 + t196 * t244 - t398 * t560 + t419 * t429 - t420 * t430, t183 * t112 + t208 * t391, -t183 * t111 - t112 * t181 - t207 * t391 - t208 * t81, t112 * t551 - t208 * t405 + (-t244 * t574 + t380 * t416 - t384 * t412 + t594) * t268 + t183 * t193 + (t208 * t546 + (t208 * t381 * t458 + (qJD(2) * t268 * t541 + t112 * t381 + t208 * t424) * t382) * qJD(1)) * t378, -t111 * t551 + t207 * t405 - t81 * t268 - t181 * t193 + (-t207 * t546 + (-t207 * t415 + (-t111 * t382 - t207 * t458) * t381) * qJD(1)) * t378, -t405 * t268 + t551 * t193 + (t268 * t546 + (t268 * t415 + (t382 * t193 + t268 * t458) * t381) * qJD(1)) * t378, g(1) * t156 - g(2) * t160 + t105 * t111 + t125 * t81 + t14 * t268 + t85 * t181 + t58 * t193 + t46 * t207 + t442 * t484 - t529 * t672, t105 * t112 + t125 * t391 + t85 * t183 - t59 * t193 + t46 * t208 + t268 * t464 - t442 * t463 + t591 * t672 + t504, t137 * t61 - t154 * t38, -t135 * t61 - t137 * t60 + t153 * t38 - t154 * t39, t111 * t137 + t154 * t80 + t180 * t61 - t207 * t38, -t111 * t135 - t153 * t80 - t180 * t60 - t207 * t39, t111 * t180 + t207 * t80, t23 * t111 + t12 * t153 + t22 * t135 + t180 * t647 - t207 * t544 + t66 * t39 + t491 * t80 + t51 * t60 + t505, -t24 * t111 + t12 * t154 + t22 * t137 - t180 * t466 - t3 * t207 - t66 * t38 - t492 * t80 + t51 * t61 + t506, -t111 * t19 + t135 * t8 + t153 * t5 - t180 * t7 - t2 * t207 + t26 * t60 - t28 * t80 + t37 * t39 + t505, -t1 * t153 - t135 * t6 + t137 * t7 + t154 * t2 + t19 * t61 - t20 * t60 - t27 * t39 - t28 * t38 - t504, t1 * t207 + t111 * t20 - t137 * t8 - t154 * t5 + t180 * t6 - t26 * t61 + t27 * t80 + t37 * t38 - t506, t1 * t27 + t20 * t6 + t5 * t37 + t26 * t8 + t2 * t28 + t19 * t7 - g(1) * (-pkin(1) * t640 - t328 * pkin(2) - pkin(3) * t214 - pkin(4) * t156 - pkin(5) * t100 + pkin(9) * t598 - pkin(11) * t213 + t157 * pkin(12) - qJ(6) * t99) - g(2) * (pkin(9) * t378 * t640 + t386 * pkin(1) + t329 * pkin(2) + t218 * pkin(3) + t160 * pkin(4) + t104 * pkin(5) + t217 * pkin(11) + t159 * pkin(12) + t103 * qJ(6)) + (g(1) * t278 + g(2) * t399) * pkin(10); 0, 0, 0, -t382 * t561, t579 * t601, -t385 * t509 + t548, t382 * t509 + t547, t477, -pkin(9) * t659 + g(1) * t434 + g(2) * t327 - g(3) * t599 + qJD(1) * t404 + t387 * t562 + t440, pkin(1) * t561 + (t363 - t522) * t533 + g(1) * t329 + g(2) * t328 + g(3) * t600 + t363 * qJD(2) - t558, t148 * t539 + t244 * t444, t148 * t512 - t149 * t539 + t296 * t242 + t244 * t295 + (-t242 * t512 - t244 * t539) * qJD(3), t148 * t611 - t244 * t479 + t296 * t419 + t381 * t396 - t403 * t641, -t149 * t611 + t242 * t479 - t295 * t419 + t381 * t403 + t396 * t641, t416 * t611 - (qJD(1) * t499 - t567) * t479, -g(1) * t237 - g(2) * t235 - g(3) * t310 - t113 * t512 - t123 * t479 - t149 * t555 + t176 * t468 - t195 * t242 + t398 * t457 + t419 * t671 - t469 * t611, g(1) * t236 + g(2) * t234 + g(3) * t309 + t113 * t539 + t124 * t479 - t148 * t555 + t176 * t444 - t195 * t244 - t398 * t581 + t419 * t675 + t430 * t611, t183 * t656 + t391 * t326, -t181 * t656 - t183 * t586 - t325 * t391 - t326 * t81, t183 * t468 - t326 * t672 - t391 * t512 + t442 * t656, t325 * t405 + t81 * t512 - t468 * t181 + (-t325 * t546 + (-t325 * t415 + (-t325 * t458 - t382 * t586) * t381) * qJD(1)) * t378 - t586 * t551, t442 * t468 + t512 * t672, -g(1) * t187 - g(2) * t185 - g(3) * t246 + t586 * t105 - t14 * t512 + t655 * t181 + t317 * t81 + t46 * t325 + t442 * t667 + t468 * t58 - t526 * t672, t656 * t105 + t655 * t183 + t317 * t391 + t46 * t326 + t442 * t670 - t464 * t512 - t468 * t59 + t585 * t672 - t450, -t137 * t589 - t282 * t38, t135 * t589 - t137 * t588 + t281 * t38 - t282 * t39, t137 * t586 - t180 * t589 + t282 * t80 - t325 * t38, -t135 * t586 - t180 * t588 - t281 * t80 - t325 * t39, t180 * t586 + t325 * t80, t490 * t80 - t544 * t325 + t211 * t39 + t12 * t281 + t588 * t51 + t586 * t23 + ((-qJD(5) * t212 + t668) * t383 + (-qJD(5) * t209 + t669) * t379) * t180 + t618 * t135 + t453, t12 * t282 + t618 * t137 + t180 * t660 - t211 * t38 - t586 * t24 - t3 * t325 - t589 * t51 - t587 * t80 - t454, -t110 * t80 + t135 * t628 + t139 * t39 + t180 * t630 - t19 * t586 - t2 * t325 + t26 * t588 + t281 * t5 + t453, -t1 * t281 - t109 * t39 - t110 * t38 - t135 * t629 - t137 * t630 - t19 * t589 + t2 * t282 - t20 * t588 + t450, t1 * t325 + t109 * t80 - t137 * t628 + t139 * t38 + t180 * t629 + t20 * t586 + t26 * t589 - t282 * t5 + t454, t1 * t109 + t5 * t139 + t2 * t110 - g(1) * (-pkin(2) * t434 + t237 * pkin(3) + t187 * pkin(4) + t134 * pkin(5) + t236 * pkin(11) + t186 * pkin(12) + t133 * qJ(6) + t329 * t556) - g(2) * (-t327 * pkin(2) + t235 * pkin(3) + t185 * pkin(4) + t132 * pkin(5) + t234 * pkin(11) + t184 * pkin(12) + t131 * qJ(6) + t328 * t556) - g(3) * (t310 * pkin(3) + t246 * pkin(4) + t189 * pkin(5) + t309 * pkin(11) + t245 * pkin(12) + t188 * qJ(6) + t378 * t446) + t628 * t26 + t629 * t20 - t630 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244 * t242, -t242 ^ 2 + t244 ^ 2, -t242 * t419 + t148, -t244 * t419 - t149, t398, -t124 * t419 - t176 * t244 + t449 - t469, g(1) * t218 + g(2) * t214 + g(3) * t269 - t123 * t419 + t176 * t242 + t430 (-qJD(4) * t244 + t398) * t380 ^ 2 + ((t148 - t412) * t380 + t442 * t183) * t384, -t380 * t81 + t384 * t391 - t658 * t183 + (-t572 - t605) * t181, -t183 * t244 + t380 * t649 + t384 * t397, t181 * t244 - t380 * t397 + t384 * t649, -t442 * t244, -pkin(3) * t81 + t105 * t606 - t124 * t181 - t58 * t244 - t427 * t638 - t442 * t530 + (pkin(11) * t672 + t568) * t380 - t646 * t384, -pkin(3) * t391 + t105 * t605 - t124 * t183 + t59 * t244 + t384 * t568 + t672 * t638 + (t592 + t665) * t442 + t646 * t380, -t38 * t380 * t383 + (-t380 * t571 + t500) * t137, t135 * t162 + t137 * t161 + (-t135 * t383 - t137 * t379) * t572 + (t621 - t383 * t39 + (t135 * t379 - t137 * t383) * qJD(5)) * t380, t38 * t384 + t500 * t180 + (t137 * t442 - t552 + t619) * t380, t384 * t39 - t501 * t180 + (-t135 * t442 + t467) * t380, t180 * t380 * t442 - t384 * t80, -t80 * t602 - t73 * t135 - t51 * t161 + (t664 * t383 + (qJD(5) * t485 + t74) * t379) * t180 + (t51 * t575 + t544 + (qJD(4) * t135 + t467) * pkin(11)) * t384 + (t51 * t570 + t12 * t379 + t442 * t23 + (t180 * t575 + t39) * pkin(11)) * t380 + t455, -t580 * t80 - t73 * t137 - t51 * t162 + t651 * t180 + (t51 * t573 + t3 + (qJD(4) * t137 + t552) * pkin(11)) * t384 + (-t51 * t571 + t12 * t383 - t442 * t24 + (t180 * t573 - t38) * pkin(11)) * t380 + t456, t2 * t384 - t291 * t80 + t320 * t39 + t501 * t26 - t616 * t180 + t615 * t135 + (-t19 * t442 + t26 * t570 + t379 * t5) * t380 + t455, t161 * t20 - t162 * t19 - t290 * t39 - t291 * t38 + t616 * t137 - t617 * t135 + t494 * t572 + (-t1 * t379 + t2 * t383 + (-t19 * t379 - t20 * t383) * qJD(5) + t449) * t380, -t1 * t384 + t290 * t80 + t320 * t38 - t500 * t26 + t617 * t180 - t615 * t137 + (t20 * t442 + t26 * t571 - t383 * t5) * t380 - t456, t1 * t290 + t5 * t320 + t2 * t291 - g(1) * (pkin(5) * t130 + pkin(11) * t218 + qJ(6) * t129) - g(2) * (pkin(5) * t128 + pkin(11) * t214 + qJ(6) * t127) - g(3) * (pkin(5) * t179 + qJ(6) * t178 + t631) + t615 * t26 + t617 * t20 + t616 * t19 + t449 * t485; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183 * t181, -t181 ^ 2 + t183 ^ 2, t242 * t181 + t394, t183 * t242 - t527, -t672, -t105 * t183 + t442 * t59 + t14 + t452, t105 * t181 + t442 * t58 - t451 + t464, t137 * t524 - t621 (-t38 - t609) * t383 + (-t39 - t607) * t379, -t137 * t183 + t180 * t524 + t620, t135 * t183 - t379 * t644 + t619, -t180 * t183, -pkin(4) * t39 - t59 * t135 - t23 * t183 + (t58 * t180 + t465) * t379 + (-t12 + (-t108 - t625) * t180 + t452) * t383, pkin(4) * t38 + t613 * t180 + t24 * t183 - t59 * t137 + t465 * t383 + (t12 + t422) * t379, t614 * t135 + t180 * t30 + t183 * t19 + t379 * t648 + t418 * t383 - t39 * t476, t135 * t29 - t137 * t30 + (t1 + t180 * t19 + (qJD(5) * t137 - t39) * pkin(12)) * t383 + (t2 - t623 + (qJD(5) * t135 - t38) * pkin(12)) * t379 + t451, -t614 * t137 - t180 * t29 - t183 * t20 + t418 * t379 - t38 * t476 - t383 * t648, -t19 * t30 - t20 * t29 + t614 * t26 + (qJD(5) * t494 + t1 * t383 + t2 * t379 + t451) * pkin(12) + (-t5 + t452) * t476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t608, -t135 ^ 2 + t645, t25, -t39 + t607, t80, -t137 * t51 + t425 + t622, t135 * t51 + t180 * t23 - t410, -t135 * t77 - t407 + t622 + 0.2e1 * t643, pkin(5) * t38 - qJ(6) * t39 + (t20 - t24) * t137 + (t19 - t593) * t135, 0.2e1 * t624 - t135 * t26 + t137 * t77 + (0.2e1 * qJD(6) - t23) * t180 + t410, t1 * qJ(6) - t2 * pkin(5) - t26 * t77 - t19 * t24 - g(1) * (-pkin(5) * t103 + qJ(6) * t104) - g(2) * (-pkin(5) * t99 + qJ(6) * t100) - g(3) * t497 + t593 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t608 - t80, t25, -t644 - t645, t407 - t623 - t643;];
tau_reg  = t4;
