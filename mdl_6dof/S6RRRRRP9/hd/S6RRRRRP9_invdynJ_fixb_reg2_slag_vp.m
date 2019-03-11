% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRP9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:11:48
% EndTime: 2019-03-10 02:12:36
% DurationCPUTime: 27.34s
% Computational Cost: add. (34162->996), mult. (82671->1271), div. (0->0), fcn. (66457->14), ass. (0->422)
t415 = sin(qJ(2));
t610 = cos(pkin(6));
t531 = pkin(1) * t610;
t393 = t415 * t531;
t414 = sin(qJ(3));
t417 = cos(qJ(3));
t481 = pkin(3) * t414 - pkin(10) * t417;
t411 = sin(pkin(6));
t418 = cos(qJ(2));
t581 = t411 * t418;
t686 = -(t393 + (pkin(8) + t481) * t581) * qJD(1) + t481 * qJD(3);
t413 = sin(qJ(4));
t416 = cos(qJ(4));
t557 = qJD(1) * t411;
t524 = t418 * t557;
t494 = t417 * t524;
t525 = t415 * t557;
t286 = t413 * t494 - t416 * t525;
t553 = qJD(3) * t417;
t684 = -t413 * t553 + t286;
t551 = qJD(4) * t416;
t683 = -t414 * t551 + t684;
t640 = cos(qJ(5));
t513 = t640 * qJD(5);
t685 = t416 * (t640 * qJD(4) + t513);
t505 = t610 * qJD(1);
t493 = pkin(1) * t505;
t312 = -pkin(8) * t525 + t418 * t493;
t483 = pkin(2) * t415 - pkin(9) * t418;
t313 = t483 * t557;
t231 = t417 * t312 + t414 * t313;
t206 = pkin(10) * t525 + t231;
t638 = pkin(3) * t417;
t482 = pkin(10) * t414 + t638;
t359 = -pkin(2) - t482;
t552 = qJD(4) * t413;
t554 = qJD(3) * t414;
t566 = -t416 * t206 + t359 * t551 + (-t416 * t554 - t417 * t552) * pkin(9) + t686 * t413;
t572 = t417 * t418;
t576 = t413 * t415;
t287 = (t416 * t572 + t576) * t557;
t475 = t416 * t553 - t287;
t682 = -t414 * t552 + t475;
t540 = pkin(9) * t554;
t681 = t686 * t416 + (t206 + t540) * t413;
t560 = pkin(8) * t581 + t393;
t309 = t610 * pkin(9) + t560;
t274 = qJD(2) * pkin(9) + qJD(1) * t309;
t468 = -pkin(2) * t418 - pkin(9) * t415 - pkin(1);
t288 = t468 * t557;
t187 = -t414 * t274 + t288 * t417;
t470 = t505 + qJD(2);
t496 = t414 * t525;
t296 = -t417 * t470 + t496;
t298 = t414 * t470 + t417 * t525;
t219 = pkin(3) * t298 + pkin(10) * t296;
t126 = -t187 * t413 + t416 * t219;
t419 = -pkin(11) - pkin(10);
t530 = qJD(4) * t419;
t680 = -pkin(4) * t298 - t126 + (-pkin(11) * t296 + t530) * t416;
t127 = t416 * t187 + t413 * t219;
t595 = t296 * t413;
t679 = pkin(11) * t595 - t413 * t530 + t127;
t574 = t416 * t417;
t395 = pkin(9) * t574;
t495 = t414 * t524;
t678 = -pkin(4) * t495 + pkin(11) * t287 + (pkin(4) * t414 - pkin(11) * t574) * qJD(3) + (-t395 + (pkin(11) * t414 - t359) * t413) * qJD(4) + t681;
t677 = -t683 * pkin(11) - t566;
t366 = -qJD(3) + t524;
t237 = -t413 * t298 - t366 * t416;
t238 = t298 * t416 - t366 * t413;
t412 = sin(qJ(5));
t144 = -t640 * t237 + t238 * t412;
t141 = t144 ^ 2;
t459 = t412 * t237 + t238 * t640;
t642 = t459 ^ 2;
t676 = -t141 + t642;
t550 = qJD(5) * t412;
t675 = pkin(4) * t550;
t674 = pkin(4) * t640;
t582 = t411 * t415;
t325 = t414 * t610 + t417 * t582;
t461 = -t325 * t416 + t413 * t581;
t673 = pkin(11) * t461;
t545 = qJDD(1) * t418;
t383 = t411 * t545;
t547 = qJD(1) * qJD(2);
t511 = t415 * t547;
t490 = t411 * t511;
t311 = qJDD(3) - t383 + t490;
t448 = qJD(3) * t470;
t510 = t418 * t547;
t546 = qJDD(1) * t415;
t454 = t510 + t546;
t426 = t411 * t454 + t448;
t501 = t610 * qJDD(1);
t464 = t501 + qJDD(2);
t437 = qJD(3) * t496 - t414 * t464;
t655 = t426 * t417 - t437;
t122 = t298 * t552 - t413 * t311 + t366 * t551 - t416 * t655;
t536 = t298 * t551 - t366 * t552 + t413 * t655;
t463 = t311 * t416 - t536;
t54 = t640 * t122 - t237 * t513 + t238 * t550 - t412 * t463;
t289 = qJD(4) + t296;
t282 = qJD(5) + t289;
t605 = t144 * t282;
t672 = -t54 + t605;
t671 = t144 * qJ(6);
t670 = t144 * t459;
t344 = t412 * t416 + t413 * t640;
t648 = qJD(4) + qJD(5);
t256 = t648 * t344;
t492 = t640 * t553;
t570 = t256 * t414 + t287 * t640 - t684 * t412 - t416 * t492;
t577 = t413 * t414;
t538 = t412 * t577;
t569 = -qJD(5) * t538 - t286 * t640 + t682 * t412 + t413 * t492 + t414 * t685;
t669 = t344 * t296 + t256;
t526 = t640 * t416;
t580 = t412 * t413;
t343 = -t526 + t580;
t567 = t343 * t296 + t580 * t648 - t685;
t55 = qJD(5) * t459 - t412 * t122 - t640 * t463;
t602 = t459 * t282;
t668 = -t55 + t602;
t273 = -pkin(2) * t470 - t312;
t161 = t296 * pkin(3) - t298 * pkin(10) + t273;
t188 = t274 * t417 + t288 * t414;
t166 = -pkin(10) * t366 + t188;
t105 = t416 * t161 - t166 * t413;
t466 = qJD(2) * t493;
t489 = pkin(1) * t501;
t533 = pkin(8) * t383 + t415 * t489 + t418 * t466;
t239 = -pkin(8) * t490 + t533;
t215 = pkin(9) * t464 + t239;
t460 = t483 * qJD(2);
t226 = (qJD(1) * t460 + qJDD(1) * t468) * t411;
t101 = t417 * t215 + t414 * t226 - t274 * t554 + t288 * t553;
t88 = pkin(10) * t311 + t101;
t555 = qJD(2) * t418;
t521 = t414 * t555;
t189 = t411 * (qJD(1) * (t415 * t553 + t521) + t414 * t546) + t414 * t448 - t417 * t464;
t509 = t411 * t546;
t498 = t415 * t466 - t418 * t489 + (t411 * t510 + t509) * pkin(8);
t216 = -pkin(2) * t464 + t498;
t96 = t189 * pkin(3) - pkin(10) * t655 + t216;
t456 = -t161 * t551 + t166 * t552 - t413 * t96 - t416 * t88;
t667 = -t105 * t289 - t456;
t651 = t495 - t554;
t335 = -pkin(8) * t582 + t418 * t531;
t316 = qJD(2) * t335;
t641 = cos(qJ(1));
t486 = t610 * t641;
t639 = sin(qJ(1));
t327 = t415 * t486 + t418 * t639;
t529 = t411 * t641;
t262 = t327 * t414 + t417 * t529;
t485 = t610 * t639;
t329 = -t415 * t485 + t418 * t641;
t528 = t411 * t639;
t266 = t329 * t414 - t417 * t528;
t324 = t414 * t582 - t417 * t610;
t450 = g(1) * t266 + g(2) * t262 + g(3) * t324;
t165 = pkin(3) * t366 - t187;
t128 = -pkin(4) * t237 + t165;
t267 = t329 * t417 + t414 * t528;
t328 = t415 * t641 + t418 * t485;
t410 = qJ(4) + qJ(5);
t402 = sin(t410);
t403 = cos(t410);
t200 = t267 * t403 + t328 * t402;
t182 = qJDD(4) + t189;
t106 = t161 * t413 + t166 * t416;
t30 = -qJD(4) * t106 - t413 * t88 + t416 * t96;
t20 = pkin(4) * t182 + pkin(11) * t122 + t30;
t25 = pkin(11) * t463 - t456;
t85 = -pkin(11) * t238 + t105;
t77 = pkin(4) * t289 + t85;
t86 = pkin(11) * t237 + t106;
t3 = t412 * t20 + t640 * t25 + t77 * t513 - t86 * t550;
t263 = t327 * t417 - t414 * t529;
t326 = t415 * t639 - t418 * t486;
t663 = t263 * t403 + t326 * t402;
t432 = g(1) * t200 + g(2) * t663 - g(3) * (-t325 * t403 + t402 * t581) - t3;
t665 = t128 * t144 + t432;
t664 = t263 * t402 - t326 * t403;
t662 = t263 * t413 - t326 * t416;
t661 = t263 * t416 + t326 * t413;
t407 = t411 ^ 2;
t660 = 0.2e1 * t407;
t625 = g(3) * t411;
t342 = t416 * t359;
t575 = t414 * t416;
t633 = pkin(9) * t413;
t251 = -pkin(11) * t575 + t342 + (-pkin(4) - t633) * t417;
t306 = t413 * t359 + t395;
t270 = -pkin(11) * t577 + t306;
t619 = t251 * t513 - t270 * t550 + t678 * t412 - t677 * t640;
t168 = t412 * t251 + t640 * t270;
t618 = -qJD(5) * t168 + t677 * t412 + t678 * t640;
t308 = -pkin(2) * t610 - t335;
t202 = t324 * pkin(3) - t325 * pkin(10) + t308;
t561 = pkin(2) * t581 + pkin(9) * t582;
t310 = -pkin(1) * t411 - t561;
t221 = t417 * t309 + t414 * t310;
t204 = -pkin(10) * t581 + t221;
t125 = t413 * t202 + t416 * t204;
t462 = t325 * t413 + t416 * t581;
t107 = -pkin(11) * t462 + t125;
t578 = t413 * t204;
t124 = t416 * t202 - t578;
t99 = pkin(4) * t324 + t124 + t673;
t62 = t640 * t107 + t412 * t99;
t367 = t419 * t413;
t368 = t419 * t416;
t612 = t367 * t513 + t368 * t550 + t680 * t412 - t679 * t640;
t278 = t412 * t367 - t640 * t368;
t611 = -qJD(5) * t278 + t679 * t412 + t680 * t640;
t659 = -t106 * t289 - t30;
t658 = qJ(6) * t459;
t487 = -t188 + (t552 + t595) * pkin(4);
t178 = qJDD(5) + t182;
t500 = pkin(4) * t513;
t636 = pkin(4) * t412;
t656 = -t178 * t636 - t282 * t500;
t230 = -t414 * t312 + t313 * t417;
t205 = -pkin(3) * t525 - t230;
t401 = pkin(9) * t553;
t563 = -t683 * pkin(4) - t205 + t401;
t207 = -t267 * t413 + t328 * t416;
t654 = -g(1) * t207 + g(2) * t662 + g(3) * t462;
t653 = (qJDD(2) + 0.2e1 * t501) * t411;
t652 = t494 - t553;
t649 = -t311 * t417 - t366 * t554;
t317 = t560 * qJD(2);
t174 = t178 * pkin(5);
t616 = t54 * qJ(6);
t647 = -t459 * qJD(6) + t174 + t616;
t199 = -t267 * t402 + t328 * t403;
t646 = -g(3) * (-t325 * t402 - t403 * t581) + g(2) * t664 - g(1) * t199;
t84 = t640 * t86;
t43 = t412 * t77 + t84;
t507 = t640 * t20 - t412 * t25;
t4 = -qJD(5) * t43 + t507;
t423 = t4 + t646;
t645 = -t128 * t459 + t423;
t556 = qJD(2) * t415;
t520 = t416 * t556;
t522 = t411 * t555;
t260 = -qJD(3) * t324 + t417 * t522;
t597 = t260 * t413;
t446 = t411 * t520 - t597;
t427 = qJD(4) * t461 + t446;
t506 = -pkin(5) * t144 - qJD(6);
t80 = t128 - t506;
t643 = -t80 * t459 + t646;
t420 = qJD(1) ^ 2;
t637 = pkin(4) * t238;
t635 = pkin(4) * t413;
t634 = pkin(4) * t416;
t624 = t311 * pkin(3);
t404 = t414 * pkin(9);
t354 = pkin(5) * t402 + t635;
t623 = pkin(9) + t354;
t82 = t412 * t86;
t42 = t640 * t77 - t82;
t31 = t42 - t658;
t28 = pkin(5) * t282 + t31;
t622 = -t31 + t28;
t319 = t414 * t526 - t538;
t621 = -t651 * pkin(5) + t570 * qJ(6) - t319 * qJD(6) + t618;
t318 = t344 * t414;
t620 = -t569 * qJ(6) - qJD(6) * t318 + t619;
t48 = t640 * t85 - t82;
t617 = qJ(6) * t55;
t615 = -t669 * qJ(6) - qJD(6) * t343 + t612;
t614 = -pkin(5) * t298 + t567 * qJ(6) - t344 * qJD(6) + t611;
t613 = -t144 * t500 - t55 * t636;
t607 = t122 * t413;
t601 = t182 * t416;
t600 = t237 * t289;
t599 = t238 * t237;
t598 = t238 * t289;
t596 = t282 * t298;
t594 = t298 * t296;
t593 = t298 * t366;
t588 = t326 * t414;
t586 = t328 * t414;
t457 = t366 * t414;
t585 = t402 * t417;
t584 = t403 * t417;
t583 = t407 * t420;
t579 = t413 * t182;
t573 = t417 * t189;
t571 = t569 * pkin(5) + t563;
t565 = -qJD(4) * t306 + t681;
t564 = t669 * pkin(5) + t487;
t352 = pkin(4) * t577 + t404;
t559 = t641 * pkin(1) + pkin(8) * t528;
t408 = t415 ^ 2;
t409 = t418 ^ 2;
t558 = t408 - t409;
t549 = qJD(6) * t144;
t543 = g(3) * t582;
t542 = g(3) * t581;
t539 = t418 * t583;
t314 = t411 * t460;
t136 = -t309 * t554 + t310 * t553 + t414 * t314 + t417 * t316;
t523 = t411 * t556;
t132 = pkin(10) * t523 + t136;
t259 = qJD(3) * t325 + t411 * t521;
t153 = t259 * pkin(3) - t260 * pkin(10) + t317;
t537 = t416 * t132 + t413 * t153 + t202 * t551;
t535 = t450 * t402;
t534 = t329 * pkin(2) + t559;
t400 = pkin(3) + t634;
t532 = pkin(9) + t635;
t515 = g(3) * t561;
t512 = pkin(1) * t660;
t355 = pkin(5) * t403 + t634;
t47 = -t412 * t85 - t84;
t61 = -t107 * t412 + t640 * t99;
t167 = t640 * t251 - t270 * t412;
t220 = -t414 * t309 + t417 * t310;
t277 = t640 * t367 + t368 * t412;
t503 = t289 * t413;
t502 = t289 * t416;
t499 = t415 * t539;
t102 = -t414 * t215 + t417 * t226 - t274 * t553 - t288 * t554;
t491 = t415 * t510;
t488 = -pkin(1) * t639 + pkin(8) * t529;
t484 = t411 * t420 * t610;
t480 = -g(1) * t664 - g(2) * t199;
t479 = g(1) * t663 - g(2) * t200;
t478 = -g(1) * t262 + g(2) * t266;
t477 = -g(1) * t326 + g(2) * t328;
t476 = g(1) * t329 + g(2) * t327;
t203 = pkin(3) * t581 - t220;
t474 = -t105 * t416 - t106 * t413;
t348 = pkin(3) + t355;
t406 = -qJ(6) + t419;
t473 = t348 * t417 - t406 * t414;
t472 = t400 * t417 - t414 * t419;
t469 = 0.2e1 * t505 + qJD(2);
t467 = pkin(9) * t328 + t534;
t465 = -t327 * pkin(2) + t488;
t137 = -t309 * t553 - t310 * t554 + t417 * t314 - t414 * t316;
t160 = -qJD(4) * t462 + t260 * t416 + t413 * t523;
t67 = -qJD(4) * t125 - t132 * t413 + t416 * t153;
t41 = pkin(4) * t259 - pkin(11) * t160 + t67;
t51 = t446 * pkin(11) + (-t578 + t673) * qJD(4) + t537;
t10 = -t107 * t550 + t412 * t41 + t640 * t51 + t99 * t513;
t89 = -t102 - t624;
t455 = -pkin(10) * t182 + t165 * t289;
t453 = g(1) * t641 + g(2) * t639;
t452 = -g(1) * (t328 * t585 + t329 * t403) - g(2) * (t326 * t585 + t327 * t403) - (-t402 * t572 + t403 * t415) * t625;
t451 = -g(1) * (-t328 * t584 + t329 * t402) - g(2) * (-t326 * t584 + t327 * t402) - (t402 * t415 + t403 * t572) * t625;
t449 = -g(1) * t267 - g(2) * t263 - g(3) * t325;
t447 = t463 * t416;
t443 = -t326 * pkin(9) + t465;
t442 = t450 - t89;
t441 = -g(1) * t328 - g(2) * t326 + t542;
t440 = -t476 - t543;
t438 = t450 * t403;
t436 = t640 * t462;
t435 = t441 * t414;
t433 = t101 * t417 - t102 * t414 - t476;
t431 = pkin(10) * qJD(4) * t289 - t442;
t11 = -qJD(5) * t62 + t640 * t41 - t412 * t51;
t171 = -t412 * t462 - t461 * t640;
t149 = pkin(4) * t462 + t203;
t429 = t432 + t617;
t63 = -pkin(4) * t463 + t89;
t21 = t55 * pkin(5) + qJDD(6) + t63;
t90 = -(-t325 * t551 - t597) * pkin(4) + (-pkin(3) * t556 - (t418 * t552 + t520) * pkin(4)) * t411 - t137;
t399 = pkin(5) + t674;
t322 = t328 * pkin(2);
t320 = t326 * pkin(2);
t315 = t560 * qJD(1);
t307 = pkin(5) * t343 - t400;
t305 = -t417 * t633 + t342;
t253 = pkin(5) * t318 + t352;
t236 = -qJ(6) * t343 + t278;
t235 = -qJ(6) * t344 + t277;
t208 = t267 * t416 + t328 * t413;
t170 = -t412 * t461 + t436;
t140 = -qJ(6) * t318 + t168;
t138 = -pkin(5) * t417 - qJ(6) * t319 + t167;
t133 = -pkin(3) * t523 - t137;
t121 = pkin(5) * t459 + t637;
t116 = -t178 * t417 - t282 * t457;
t115 = t178 * t324 + t259 * t282;
t103 = t170 * pkin(5) + t149;
t74 = qJD(5) * t171 + t412 * t160 - t427 * t640;
t73 = qJD(5) * t436 - t160 * t640 - t412 * t427 - t461 * t550;
t66 = -t204 * t552 + t537;
t60 = t144 * t298 - t178 * t343 - t282 * t669;
t59 = t178 * t344 - t282 * t567 - t298 * t459;
t49 = -qJ(6) * t170 + t62;
t45 = t74 * pkin(5) + t90;
t44 = pkin(5) * t324 - qJ(6) * t171 + t61;
t34 = t48 - t658;
t33 = t47 + t671;
t32 = t43 - t671;
t27 = t144 * t669 + t343 * t55;
t26 = -t344 * t54 - t459 * t567;
t23 = t144 * t569 + t318 * t55;
t22 = -t319 * t54 - t459 * t570;
t17 = t144 * t74 + t170 * t55;
t16 = -t171 * t54 - t459 * t73;
t15 = t144 * t457 - t178 * t318 - t282 * t569 + t417 * t55;
t14 = t178 * t319 - t282 * t570 + t417 * t54 - t457 * t459;
t13 = -t144 * t259 - t170 * t178 - t282 * t74 - t324 * t55;
t12 = t171 * t178 + t259 * t459 - t282 * t73 - t324 * t54;
t9 = t144 * t567 + t343 * t54 - t344 * t55 - t459 * t669;
t8 = t144 * t570 + t318 * t54 - t319 * t55 - t459 * t569;
t7 = -qJ(6) * t74 - qJD(6) * t170 + t10;
t6 = t144 * t73 + t170 * t54 - t171 * t55 - t459 * t74;
t5 = t259 * pkin(5) + t73 * qJ(6) - t171 * qJD(6) + t11;
t2 = t3 - t549 - t617;
t1 = t4 + t647;
t18 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t639 - g(2) * t641, t453, 0, 0 (qJDD(1) * t408 + 0.2e1 * t491) * t407 (t415 * t545 - t547 * t558) * t660, t415 * t653 + t469 * t522 (qJDD(1) * t409 - 0.2e1 * t491) * t407, t418 * t653 - t469 * t523, t464 * t610, -t317 * t470 + t335 * t464 - t498 * t610 + g(1) * t327 - g(2) * t329 + (-t511 + t545) * t512, -t239 * t610 - t316 * t470 - t454 * t512 - t464 * t560 + t477 ((-t312 * qJD(2) + qJDD(1) * t560 + t239) * t418 + (-qJD(2) * t315 - qJDD(1) * t335 + t498) * t415 - t453) * t411, t407 * qJDD(1) * pkin(1) ^ 2 - g(1) * t488 - g(2) * t559 + t239 * t560 - t312 * t317 + t315 * t316 - t335 * t498, t298 * t260 + t325 * t655, -t325 * t189 - t298 * t259 - t260 * t296 - t324 * t655, -t260 * t366 + t325 * t311 + (t298 * t556 - t418 * t655) * t411, t189 * t324 + t259 * t296, t259 * t366 - t311 * t324 + (t189 * t418 - t296 * t556) * t411 (-t311 * t418 - t366 * t556) * t411, g(1) * t263 - g(2) * t267 - t137 * t366 + t189 * t308 + t216 * t324 + t220 * t311 + t259 * t273 + t296 * t317 + (-t102 * t418 + t187 * t556) * t411, t101 * t581 + t136 * t366 - t188 * t523 + t216 * t325 - t221 * t311 + t273 * t260 + t317 * t298 + t308 * t655 + t478, -t101 * t324 - t102 * t325 - t136 * t296 - t137 * t298 - t187 * t260 - t188 * t259 - t221 * t189 - t220 * t655 - t477, -g(1) * t443 - g(2) * t467 + t101 * t221 + t102 * t220 + t188 * t136 + t187 * t137 + t216 * t308 + t273 * t317, t122 * t461 + t160 * t238, t122 * t462 + t160 * t237 + t238 * t427 - t461 * t463, -t122 * t324 + t160 * t289 - t182 * t461 + t238 * t259, t237 * t427 - t462 * t463, -t182 * t462 + t237 * t259 + t289 * t427 + t324 * t463, t182 * t324 + t259 * t289, g(1) * t661 - g(2) * t208 + t105 * t259 + t124 * t182 - t133 * t237 - t165 * t427 - t203 * t463 + t67 * t289 + t30 * t324 + t89 * t462, -g(1) * t662 - g(2) * t207 - t106 * t259 - t203 * t122 - t125 * t182 + t133 * t238 + t165 * t160 - t66 * t289 + t456 * t324 - t89 * t461, -t105 * t160 + t106 * t427 + t124 * t122 + t125 * t463 + t66 * t237 - t67 * t238 + t30 * t461 + t456 * t462 - t478, -t456 * t125 + t106 * t66 + t30 * t124 + t105 * t67 + t89 * t203 + t165 * t133 - g(1) * (-pkin(3) * t263 - pkin(10) * t262 + t443) - g(2) * (pkin(3) * t267 + pkin(10) * t266 + t467) t16, t6, t12, t17, t13, t115, t11 * t282 + t128 * t74 + t144 * t90 + t149 * t55 + t170 * t63 + t178 * t61 + t259 * t42 + t324 * t4 + t479, -t10 * t282 - t128 * t73 - t149 * t54 + t171 * t63 - t178 * t62 - t259 * t43 - t3 * t324 + t459 * t90 + t480, -t10 * t144 - t11 * t459 - t170 * t3 - t171 * t4 + t42 * t73 - t43 * t74 + t54 * t61 - t55 * t62 - t478, t3 * t62 + t43 * t10 + t4 * t61 + t42 * t11 + t63 * t149 + t128 * t90 - g(1) * (t262 * t419 - t263 * t400 - t326 * t532 + t465) - g(2) * (-t266 * t419 + t267 * t400 + t328 * t532 + t534) t16, t6, t12, t17, t13, t115, t1 * t324 + t103 * t55 + t144 * t45 + t170 * t21 + t178 * t44 + t259 * t28 + t282 * t5 + t74 * t80 + t479, -t103 * t54 + t171 * t21 - t178 * t49 - t2 * t324 - t259 * t32 - t282 * t7 + t45 * t459 - t73 * t80 + t480, -t1 * t171 - t144 * t7 - t170 * t2 + t28 * t73 - t32 * t74 + t44 * t54 - t459 * t5 - t49 * t55 - t478, t2 * t49 + t32 * t7 + t1 * t44 + t28 * t5 + t21 * t103 + t80 * t45 - g(1) * (t262 * t406 - t263 * t348 - t326 * t623 + t465) - g(2) * (-t266 * t406 + t267 * t348 + t328 * t623 + t534); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t499, t558 * t583, -t418 * t484 + t509, t499, t415 * t484 + t383, t464, pkin(1) * t415 * t583 + t315 * t470 - t441 - t498, pkin(1) * t539 + t312 * t470 + (pkin(8) * t547 + g(3)) * t582 + t476 - t533, 0, 0, -t437 * t414 + (t414 * t426 - t593) * t417, -t414 * t189 + t296 * t652 + t298 * t651 + t417 * t655, -t366 * t553 + t311 * t414 + (-t298 * t415 + t366 * t572) * t557, -t296 * t457 - t573 (t296 * t415 - t418 * t457) * t557 - t649, t366 * t525, -t187 * t525 - pkin(2) * t189 + t230 * t366 - t296 * t315 + (-pkin(9) * t311 - t273 * t366) * t414 + (qJD(3) * pkin(9) * t366 - t216 - t441) * t417, -pkin(2) * t655 - g(1) * t586 - g(2) * t588 + t188 * t525 - t231 * t366 - t315 * t298 + (t216 + t542) * t414 - t652 * t273 + t649 * pkin(9), -pkin(9) * t573 + t655 * t404 + t433 - t543 + (t230 + t401) * t298 + (t231 + t540) * t296 + t651 * t188 + t652 * t187, -t216 * pkin(2) - t188 * t231 - t187 * t230 - t273 * t315 + g(1) * t322 + g(2) * t320 - t515 + ((-t187 * t417 - t188 * t414) * qJD(3) + t433) * pkin(9), -t122 * t575 + t682 * t238, -t287 * t237 + t238 * t286 + (t237 * t416 - t238 * t413) * t553 + (t447 + t607 + (-t237 * t413 - t238 * t416) * qJD(4)) * t414, t122 * t417 + t475 * t289 + (-t238 * t366 - t289 * t552 + t601) * t414, t683 * t237 - t463 * t577, -t463 * t417 + t684 * t289 + (-t237 * t366 - t289 * t551 - t579) * t414, -t182 * t417 - t289 * t457, -t165 * t286 + t305 * t182 + t205 * t237 + t565 * t289 + t440 * t413 + (-t30 + (-pkin(9) * t237 + t165 * t413) * qJD(3) - t441 * t416) * t417 + (-pkin(9) * t463 - t105 * t366 + t165 * t551 + t89 * t413) * t414, -t165 * t287 - t306 * t182 - t205 * t238 - t566 * t289 + t440 * t416 + (-t456 + (pkin(9) * t238 + t165 * t416) * qJD(3) + t441 * t413) * t417 + (-pkin(9) * t122 + t106 * t366 - t165 * t552 + t89 * t416) * t414, t306 * t463 + t305 * t122 + t106 * t286 + t105 * t287 - t565 * t238 + t566 * t237 + t474 * t553 + (t456 * t413 - t30 * t416 + (t105 * t413 - t106 * t416) * qJD(4) - t441) * t414, -t456 * t306 + t30 * t305 - t165 * t205 - g(1) * (-pkin(10) * t586 - t328 * t638 - t322) - g(2) * (-pkin(10) * t588 - t326 * t638 - t320) - g(3) * (t482 * t581 + t561) + t566 * t106 + t565 * t105 + (t165 * t553 + t414 * t89 - t476) * pkin(9), t22, t8, t14, t23, t15, t116, t128 * t569 + t144 * t563 + t167 * t178 + t282 * t618 + t318 * t63 + t352 * t55 - t4 * t417 - t42 * t457 + t451, -t128 * t570 - t168 * t178 - t282 * t619 + t3 * t417 + t319 * t63 - t352 * t54 + t43 * t457 + t459 * t563 + t452, -t144 * t619 + t167 * t54 - t168 * t55 - t3 * t318 - t319 * t4 + t42 * t570 - t43 * t569 - t459 * t618 - t435, t3 * t168 + t4 * t167 + t63 * t352 - g(1) * (-t328 * t472 + t329 * t532 - t322) - g(2) * (-t326 * t472 + t327 * t532 - t320) - t515 + t619 * t43 + t618 * t42 - (pkin(4) * t576 + t418 * t472) * t625 + t563 * t128, t22, t8, t14, t23, t15, t116, -t1 * t417 + t138 * t178 + t144 * t571 + t21 * t318 + t253 * t55 - t28 * t457 + t282 * t621 + t569 * t80 + t451, -t140 * t178 + t2 * t417 + t21 * t319 - t253 * t54 - t282 * t620 + t32 * t457 + t459 * t571 - t570 * t80 + t452, -t1 * t319 + t138 * t54 - t140 * t55 - t144 * t620 - t2 * t318 + t28 * t570 - t32 * t569 - t459 * t621 - t435, t2 * t140 + t1 * t138 + t21 * t253 - g(1) * (-t328 * t473 + t329 * t623 - t322) - g(2) * (-t326 * t473 + t327 * t623 - t320) - t515 + t571 * t80 - (t354 * t415 + t418 * t473) * t625 + t620 * t32 + t621 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t594, -t296 ^ 2 + t298 ^ 2, -t296 * t366 + t655, -t594, -t189 - t593, t311, -t188 * t366 - t273 * t298 + t102 + t450, -t187 * t366 + t273 * t296 - t101 - t449, 0, 0, t238 * t502 - t607 (-t122 + t600) * t416 + (t463 - t598) * t413, -t238 * t298 + t289 * t502 + t579, -t237 * t503 + t447, -t237 * t298 - t289 * t503 + t601, -t289 * t298, -pkin(3) * t536 - t126 * t289 - t105 * t298 + t188 * t237 + t455 * t413 + (-t431 + t624) * t416, pkin(3) * t122 + t106 * t298 + t127 * t289 - t188 * t238 + t413 * t431 + t416 * t455, t126 * t238 - t127 * t237 + ((qJD(4) * t238 + t463) * pkin(10) + t667) * t416 + ((-qJD(4) * t237 - t122) * pkin(10) + t659) * t413 + t449, -t105 * t126 - t106 * t127 - t165 * t188 + t442 * pkin(3) + (qJD(4) * t474 - t30 * t413 - t416 * t456 + t449) * pkin(10), t26, t9, t59, t27, t60, -t596, t128 * t669 + t144 * t487 + t178 * t277 + t282 * t611 - t298 * t42 + t343 * t63 - t400 * t55 + t438, -t128 * t567 - t178 * t278 - t282 * t612 + t298 * t43 + t344 * t63 + t400 * t54 + t459 * t487 - t535, -t144 * t612 + t277 * t54 - t278 * t55 - t3 * t343 - t344 * t4 + t42 * t567 - t43 * t669 - t459 * t611 + t449, t3 * t278 + t4 * t277 - t63 * t400 - g(1) * (-t266 * t400 - t267 * t419) - g(2) * (-t262 * t400 - t263 * t419) - g(3) * (-t324 * t400 - t325 * t419) + t612 * t43 + t611 * t42 + t487 * t128, t26, t9, t59, t27, t60, -t596, t144 * t564 + t178 * t235 + t21 * t343 - t28 * t298 + t282 * t614 + t307 * t55 + t669 * t80 + t438, -t178 * t236 + t21 * t344 - t282 * t615 + t298 * t32 - t307 * t54 + t459 * t564 - t567 * t80 - t535, -t1 * t344 - t144 * t615 - t2 * t343 + t235 * t54 - t236 * t55 + t28 * t567 - t32 * t669 - t459 * t614 + t449, t2 * t236 + t1 * t235 + t21 * t307 - g(1) * (-t266 * t348 - t267 * t406) - g(2) * (-t262 * t348 - t263 * t406) - g(3) * (-t324 * t348 - t325 * t406) + t564 * t80 + t615 * t32 + t614 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t599, -t237 ^ 2 + t238 ^ 2, -t122 - t600, t599, t463 + t598, t182, -t165 * t238 + t654 - t659, g(1) * t208 + g(2) * t661 - g(3) * t461 - t165 * t237 - t667, 0, 0, t670, t676, t672, -t670, t668, t178, -t47 * t282 + (-t144 * t238 + t178 * t640 - t282 * t550) * pkin(4) + t645, t282 * t48 - t459 * t637 + t656 + t665, t54 * t674 + t613 + (t43 + t47 + t675) * t459 + (t48 - t42) * t144, -t42 * t47 - t43 * t48 + (t3 * t412 + t4 * t640 - t128 * t238 + (-t42 * t412 + t43 * t640) * qJD(5) + t654) * pkin(4), t670, t676, t672, -t670, t668, t178, -t121 * t144 + t399 * t178 - t33 * t282 + (-t84 + (-pkin(4) * t282 - t77) * t412) * qJD(5) + t507 + t643 + t647, -t121 * t459 + t144 * t80 + t282 * t34 + t429 + t549 + t656, t399 * t54 + t613 + (t32 + t33 + t675) * t459 + (t34 - t28) * t144, t1 * t399 - t32 * t34 - t28 * t33 - t80 * t121 - g(1) * (-t267 * t354 + t328 * t355) - g(2) * (-t263 * t354 + t326 * t355) - g(3) * (-t325 * t354 - t355 * t581) + (t2 * t412 + (-t28 * t412 + t32 * t640) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t670, t676, t672, -t670, t668, t178, t43 * t282 + t645, t282 * t42 + t665, 0, 0, t670, t676, t672, -t670, t668, t178, t616 + t32 * t282 + 0.2e1 * t174 + (t506 - t80) * t459 + t423, -pkin(5) * t642 + t282 * t31 + (qJD(6) + t80) * t144 + t429, pkin(5) * t54 - t144 * t622, t622 * t32 + (t1 + t643) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 + t602, -t54 - t605, -t141 - t642, t32 * t144 + t28 * t459 + t21 - t450;];
tau_reg  = t18;
