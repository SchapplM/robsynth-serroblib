% Calculate inertial parameters regressor of coriolis joint torque vector for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% tauc_reg [7x(7*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S7RRRRRRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 07:42:54
% EndTime: 2019-03-10 07:45:34
% DurationCPUTime: 70.63s
% Computational Cost: add. (37035->1210), mult. (90317->1724), div. (0->0), fcn. (72186->12), ass. (0->495)
t299 = sin(qJ(4));
t301 = sin(qJ(2));
t304 = cos(qJ(3));
t303 = cos(qJ(4));
t305 = cos(qJ(2));
t579 = t303 * t305;
t235 = (t299 * t301 + t304 * t579) * qJD(1);
t300 = sin(qJ(3));
t554 = qJD(4) * t299;
t699 = -t300 * t554 + t235;
t563 = qJD(1) * t305;
t494 = t304 * t563;
t557 = qJD(3) * t304;
t698 = t557 + t494;
t490 = t303 * t557;
t697 = t490 + t699;
t696 = t698 * pkin(2);
t518 = -pkin(2) * t303 - pkin(3);
t430 = t300 * t518;
t559 = qJD(3) * t300;
t481 = t303 * t559;
t552 = qJD(4) * t304;
t695 = pkin(3) * t559 + (t299 * t552 + t481) * pkin(2) - t430 * t563;
t694 = pkin(3) * t697 + t696;
t298 = sin(qJ(5));
t653 = cos(qJ(5));
t474 = t653 * qJD(4);
t477 = qJD(5) * t653;
t480 = qJD(3) * t653;
t496 = t300 * t563;
t574 = (-t303 * t480 - t477) * t304 + (t299 * t474 + (qJD(5) * t303 + qJD(3)) * t298) * t300 - t235 * t653 + t298 * t496;
t565 = qJD(1) * t301;
t234 = t299 * t494 - t303 * t565;
t553 = qJD(4) * t303;
t487 = t300 * t553;
t360 = t299 * t557 + t234 + t487;
t256 = t518 * t304;
t672 = pkin(3) * t303 + pkin(2);
t258 = t672 * t300;
t662 = -t298 * t256 + t653 * t258;
t578 = qJD(5) * t662 + t694 * t298 + t653 * t695;
t449 = qJD(3) + t563;
t415 = t303 * t449;
t291 = t300 * qJD(2);
t564 = qJD(1) * t304;
t392 = -t301 * t564 + t291;
t670 = t299 * t392;
t209 = t415 - t670;
t400 = t303 * t392;
t211 = -t299 * t449 - t400;
t297 = sin(qJ(6));
t652 = cos(qJ(6));
t440 = t652 * t653;
t133 = -t209 * t440 + t297 * t211;
t411 = qJD(5) * t440;
t384 = t411 - t133;
t548 = qJD(6) * t297;
t693 = t298 * t548 - t384;
t503 = t653 * t300;
t248 = t298 * t304 + t303 * t503;
t475 = qJD(6) * t652;
t432 = t299 * t475;
t418 = t300 * t432;
t615 = t418 - t574 * t652 + (-qJD(6) * t248 + t360) * t297;
t497 = t300 * t565;
t561 = qJD(2) * t304;
t253 = t497 + t561;
t510 = t303 * t653;
t178 = -t253 * t510 + t298 * t392;
t476 = qJD(5) * t652;
t433 = t298 * t476;
t473 = t653 * qJD(6);
t592 = t297 * t299;
t575 = t178 * t652 - t253 * t592 - (-qJD(4) * t440 + t475) * t303 - (t433 + (t473 - qJD(4)) * t297) * t299;
t446 = t305 * t503;
t580 = t303 * t304;
t573 = qJD(1) * t446 + t503 * qJD(3) + t248 * qJD(5) + (qJD(3) * t580 + t699) * t298;
t587 = t298 * t303;
t177 = -t253 * t587 - t392 * t653;
t438 = t299 * t477;
t350 = t298 * t553 + t177 + t438;
t196 = t256 * t653 + t298 * t258;
t452 = pkin(2) * t496;
t427 = t299 * t452;
t585 = t299 * t300;
t523 = t297 * t585;
t456 = pkin(2) * t523;
t651 = pkin(2) * t304;
t617 = -t297 * t427 + t196 * t548 - qJD(3) * t456 - (-t297 * t553 - t432) * t651 - t578 * t652;
t577 = qJD(5) * t196 + t298 * t695 - t694 * t653;
t549 = qJD(5) * t299;
t660 = -t298 * t549 + t303 * t474;
t692 = t178 + t660;
t482 = t291 * t305;
t372 = t301 * t557 + t482;
t540 = qJD(2) * qJD(3);
t470 = t300 * t540;
t207 = (qJD(1) * t372 - t470) * pkin(2);
t250 = t392 * pkin(2);
t691 = t449 * t250 + t207;
t541 = qJD(1) * qJD(2);
t472 = t305 * t541;
t286 = t304 * t540;
t558 = qJD(3) * t301;
t471 = qJD(1) * t558;
t568 = t300 * t471 + t286;
t386 = t304 * t472 - t568;
t364 = t299 * t386;
t469 = t301 * t541;
t276 = t303 * t469;
t408 = qJD(4) * t449;
t570 = t299 * t408 + t276;
t340 = t364 - t570;
t394 = qJD(4) * t392;
t376 = t303 * t394;
t690 = -0.2e1 * t376 + t340;
t689 = -pkin(4) * t573 + t617;
t688 = -t615 * pkin(4) + t577;
t543 = qJD(4) - t253;
t165 = t211 * t653 + t298 * t543;
t332 = qJD(5) + t209;
t190 = t652 * t332;
t114 = t165 * t297 - t190;
t542 = qJD(7) - t114;
t116 = t165 * t652 + t297 * t332;
t237 = t653 * t543;
t163 = t211 * t298 - t237;
t153 = qJD(6) + t163;
t296 = sin(qJ(7));
t302 = cos(qJ(7));
t68 = t116 * t296 + t153 * t302;
t687 = t542 * t68;
t70 = t116 * t302 - t296 * t153;
t686 = t542 * t70;
t478 = qJD(4) * t652;
t434 = t303 * t478;
t479 = qJD(3) * t652;
t618 = -t427 * t652 + t196 * t475 - (-t304 * t434 + (t300 * t479 + t304 * t548) * t299) * pkin(2) + t578 * t297;
t183 = pkin(3) * t392 + t250 * t303;
t337 = t672 * t253;
t123 = t183 * t298 - t653 * t337;
t685 = -pkin(3) * t660 - pkin(4) * t575 - t123;
t124 = t183 * t653 + t298 * t337;
t102 = t652 * t124 + t250 * t592;
t484 = t299 * t548;
t684 = -t102 + (t299 * t411 + (t434 - t484) * t298) * pkin(3) + t350 * pkin(4);
t251 = t253 * pkin(2);
t514 = t299 * t653;
t600 = t209 * t298;
t147 = pkin(3) * t600 + t251 * t514;
t581 = t303 * t251;
t118 = t652 * t147 - t297 * t581;
t435 = t297 * t473;
t550 = qJD(5) * t298;
t683 = -t118 + (-t435 - t433) * pkin(3) + (-t600 - t550) * pkin(4);
t504 = t653 * t209;
t589 = t298 * t299;
t146 = pkin(3) * t504 - t251 * t589;
t682 = -pkin(3) * t477 + t693 * pkin(4) - t146;
t681 = t209 * t543;
t680 = t211 * t543;
t678 = -t400 - t211;
t492 = t300 * t558;
t560 = qJD(2) * t305;
t373 = -t304 * t560 + t492;
t293 = t301 ^ 2;
t413 = t305 * t449;
t677 = t293 * qJD(1) - t413;
t566 = -t305 ^ 2 + t293;
t656 = t566 * qJD(1);
t547 = qJD(7) * t296;
t342 = pkin(3) * t211 - t250;
t166 = t653 * t342;
t180 = pkin(3) * t543 - t581;
t112 = t298 * t180 + t166;
t61 = pkin(4) * t116 + t112;
t113 = t180 * t653 - t298 * t342;
t380 = -t113 * t652 + t251 * t592;
t62 = -t153 * pkin(4) - t380;
t208 = (qJD(1) * t373 + t286) * pkin(2);
t571 = -t299 * t394 + t303 * t408;
t387 = t299 * t469 - t571;
t667 = t303 * t386;
t324 = -t387 - t667;
t322 = t298 * t324;
t569 = t300 * t472 + t304 * t471;
t404 = t470 - t569;
t582 = t303 * t207;
t119 = pkin(3) * t404 + t251 * t554 - t582;
t507 = t653 * t119;
t551 = qJD(5) * t112;
t37 = pkin(3) * t322 + t208 * t298 + t507 - t551;
t488 = t251 * t553;
t586 = t299 * t207;
t13 = -t113 * t548 - t251 * t432 + t652 * t37 + (-t488 - t586) * t297;
t318 = -t653 * t404 - t322;
t67 = qJD(5) * t165 + t318;
t7 = -pkin(4) * t67 + t13;
t317 = -t376 + t340;
t66 = -qJD(5) * t237 + t211 * t550 - t298 * t404 + t653 * t324;
t27 = -qJD(6) * t190 + t165 * t548 - t297 * t317 + t652 * t66;
t310 = pkin(3) * t324 + t208;
t466 = t298 * t119 - t653 * t310;
t38 = qJD(5) * t113 + t466;
t8 = -t27 * pkin(4) + t38;
t2 = -t296 * t7 + t61 * t547 + (-qJD(7) * t62 - t8) * t302;
t22 = -t296 * t61 + t302 * t62;
t676 = t542 * t22 + t2;
t425 = t296 * t62 + t302 * t61;
t1 = -qJD(7) * t425 - t296 * t8 + t302 * t7;
t675 = -t542 * t425 - t1;
t379 = t300 * t404;
t674 = t304 * ((-qJD(3) + t563) * t561 + qJD(3) * (qJD(4) - t497) - t568) + t379;
t673 = 0.2e1 * t299;
t583 = t301 * t304;
t244 = t299 * t583 + t579;
t671 = t244 * t317;
t669 = t300 * t449;
t668 = t303 * t317;
t385 = t251 * t449 - t208;
t666 = t385 * t300;
t352 = t386 * t583;
t665 = t392 * t253;
t314 = t653 * t317;
t245 = -t305 * t299 + t301 * t580;
t215 = -pkin(2) * t583 - pkin(3) * t245;
t239 = t301 * t430;
t663 = t653 * t215 - t298 * t239;
t395 = qJD(3) * t392;
t661 = -t386 * t300 + t304 * t395;
t612 = qJD(7) * t68;
t5 = t302 * t27 + t296 * t67 + t612;
t409 = qJD(3) * t449;
t659 = qJD(5) * t332;
t513 = t299 * t652;
t14 = qJD(6) * t380 - t207 * t513 - t251 * t434 - t297 * t37;
t658 = -t112 * t116 - t153 * t380 + t14;
t6 = t302 * (qJD(7) * t116 + t67) - t153 * t547 - t27 * t296;
t28 = qJD(6) * t116 - t297 * t66 - t652 * t317;
t176 = t652 * t196 - t592 * t651;
t243 = t300 * t587 - t304 * t653;
t134 = pkin(4) * t243 + t176;
t201 = -t248 * t652 - t523;
t135 = pkin(4) * t201 - t662;
t75 = -t134 * t296 - t135 * t302;
t654 = -qJD(7) * t75 + t296 * t688 + t302 * t689;
t649 = pkin(4) * t114;
t648 = t296 * t5;
t647 = t296 * t6;
t646 = t302 * t5;
t645 = t302 * t6;
t644 = t70 * t68;
t76 = t134 * t302 - t135 * t296;
t643 = -qJD(7) * t76 + t296 * t689 - t302 * t688;
t159 = t298 * t215 + t239 * t653;
t443 = t300 * t513;
t129 = -pkin(2) * t301 * t443 - t297 * t159;
t417 = t300 * t434;
t439 = t304 * t479;
t331 = -t301 * t490 + (t301 * t554 - t303 * t560) * t300;
t155 = pkin(2) * t331 - pkin(3) * t372;
t448 = -qJD(4) + t561;
t174 = t448 * t579 + (-t481 + (qJD(2) - t552) * t299) * t301;
t338 = pkin(2) * t373 - pkin(3) * t174;
t59 = qJD(5) * t663 + t653 * t155 + t298 * t338;
t32 = -t159 * t475 - t297 * t59 + (-t301 * t417 + (-t301 * t439 + (t301 * t548 - t560 * t652) * t300) * t299) * pkin(2);
t80 = -t297 * t113 - t251 * t513;
t642 = t14 * t129 + t80 * t32;
t247 = -t297 * t303 + t299 * t440;
t214 = -pkin(3) * t514 - t247 * pkin(4);
t238 = (pkin(3) * t652 + pkin(4)) * t589;
t422 = t214 * t296 - t238 * t302;
t641 = -qJD(7) * t422 + t296 * t684 + t302 * t685;
t156 = -t214 * t302 - t238 * t296;
t640 = -qJD(7) * t156 + t296 * t685 - t302 * t684;
t639 = t243 * t547 + (qJD(7) * t201 - t573) * t302 - t615 * t296;
t200 = -t247 * t302 + t296 * t589;
t638 = -qJD(7) * t200 + t296 * t575 + t302 * t350;
t588 = t298 * t302;
t520 = t299 * t588;
t637 = -qJD(7) * t520 - t247 * t547 - t296 * t350 + t302 * t575;
t636 = t14 * t296;
t635 = t14 * t302;
t445 = t301 * t503;
t202 = t245 * t298 + t445;
t634 = t202 * t67;
t633 = t243 * t67;
t632 = t27 * t297;
t631 = t28 * t296;
t630 = t28 * t302;
t629 = t297 * t67;
t628 = t297 * t80;
t627 = t302 * t70;
t626 = t38 * t297;
t625 = t38 * t298;
t624 = t66 * t298;
t101 = -t297 * t124 + t250 * t513;
t623 = t80 * t101;
t149 = t201 * t296 - t243 * t302;
t621 = qJD(7) * t149 - t296 * t573 + t302 * t615;
t424 = pkin(3) * t440;
t255 = pkin(4) * t653 + t424;
t538 = t652 * pkin(4);
t257 = (t538 + pkin(3)) * t298;
t193 = -t255 * t296 - t257 * t302;
t620 = qJD(7) * t193 + t296 * t682 + t302 * t683;
t195 = t255 * t302 - t257 * t296;
t619 = -qJD(7) * t195 - t296 * t683 + t302 * t682;
t616 = -t234 * t652 + t248 * t475 - t297 * t574 - t299 * t439 + t300 * t484 - t417;
t512 = t302 * t653;
t614 = -qJD(7) * t512 + (t600 + (qJD(7) * t652 + qJD(5)) * t298) * t296 + t693 * t302;
t511 = t302 * t652;
t246 = t296 * t653 + t298 * t511;
t613 = t246 * qJD(7) - t296 * t693 + t332 * t588;
t610 = t112 * t123;
t609 = t112 * t146;
t608 = t114 * t297;
t607 = t116 * t114;
t606 = t116 * t153;
t458 = t153 * t297;
t605 = t163 * t298;
t604 = t165 * t163;
t603 = t165 * t298;
t602 = t207 * t300;
t601 = t207 * t304;
t599 = t209 * t303;
t598 = t211 * t209;
t597 = t211 * t299;
t596 = t250 * t304;
t292 = t299 ^ 2;
t595 = t251 * t292;
t594 = t251 * t299;
t593 = t297 * t298;
t584 = t300 * t301;
t412 = qJD(6) * t440;
t576 = t253 * t513 - t303 * t548 + (t412 - t478) * t299 + t692 * t297;
t572 = (-pkin(2) * t559 - t452) * t595;
t294 = t303 ^ 2;
t567 = t292 + t294;
t562 = qJD(2) * t301;
t556 = qJD(3) * t305;
t555 = qJD(4) * t209;
t546 = qJD(7) * t302;
t544 = t305 * qJD(4);
t539 = pkin(2) * t584;
t536 = t28 * t593;
t535 = t67 * t589;
t534 = t13 * t652;
t533 = t27 * t653;
t532 = t28 * t653;
t26 = t28 * t652;
t531 = t37 * t653;
t530 = t38 * t652;
t529 = t66 * t653;
t528 = t67 * t653;
t527 = t653 * t38;
t526 = t652 * t67;
t525 = t652 * t380;
t524 = t207 * t584;
t522 = t298 * t584;
t307 = qJD(1) ^ 2;
t521 = t301 * t307 * t305;
t519 = t208 * t300 * pkin(2) + t250 * t696;
t517 = t296 * t652;
t516 = t297 * t653;
t515 = t298 * t652;
t509 = t653 * t112;
t508 = t653 * t116;
t506 = t653 * t163;
t505 = t653 * t165;
t502 = t652 * t114;
t501 = t652 * t116;
t500 = t652 * t153;
t499 = t652 * t163;
t498 = t652 * t303;
t117 = -t297 * t147 - t251 * t498;
t468 = pkin(3) * t550 * t628 - t80 * t117;
t463 = t567 * t207;
t462 = t567 * t251;
t460 = t298 * t332;
t459 = t299 * t543;
t457 = t542 * t296;
t454 = -t253 + t561;
t453 = t542 * t538;
t447 = t14 * t516;
t444 = t297 * t508;
t442 = t301 * t487;
t175 = -t297 * t196 - t513 * t651;
t441 = t14 * t175 - t618 * t80;
t436 = t297 * t477;
t431 = t305 * t469;
t60 = qJD(5) * t159 + t298 * t155 - t653 * t338;
t426 = t112 * t60 - t38 * t663;
t130 = t652 * t159 - t301 * t456;
t103 = -pkin(4) * t202 + t130;
t203 = t245 * t653 - t522;
t152 = t203 * t652 + t244 * t297;
t91 = pkin(4) * t152 - t663;
t49 = t103 * t302 - t296 * t91;
t48 = -t103 * t296 - t302 * t91;
t110 = t152 * t296 + t202 * t302;
t423 = -t597 + t599;
t421 = t250 * t300 + t251 * t304;
t420 = t114 * t440;
t419 = t298 * t432;
t410 = qJD(1) * t449;
t83 = t302 * t165 - t296 * t499;
t407 = t296 * t475 - t83;
t84 = -t296 * t165 - t302 * t499;
t406 = -t302 * t475 + t84;
t132 = -t211 * t652 - t297 * t504;
t405 = t436 - t132;
t403 = -t525 - t628;
t402 = t542 * t546 - t631;
t401 = t542 * t547 + t630;
t399 = t305 * t392;
t398 = -pkin(4) * t28 - t542 * t80;
t397 = -t203 * t297 + t244 * t652;
t396 = -t113 * t298 + t509;
t393 = t449 * t583;
t391 = t301 * t409;
t390 = t112 * t577 - t38 * t662;
t383 = t415 + qJD(5);
t382 = -t153 * t475 - t629;
t381 = -t153 * t548 + t526;
t378 = t404 * t304;
t377 = qJD(2) * t399;
t374 = qJD(1) * t399;
t369 = t300 * t391;
t368 = t477 + t504;
t367 = t475 + t499;
t365 = (t251 * t300 - t596) * t560;
t363 = t387 * t304;
t361 = t301 * t379;
t356 = -t22 * t517 + t425 * t511;
t354 = qJD(3) * t543;
t353 = qJD(4) * t543;
t351 = -t251 * t669 + t601;
t348 = t383 * t553;
t347 = -t376 - t570;
t346 = t300 * t305 * t543;
t345 = t543 * t554;
t344 = t300 * t353;
t343 = t380 * t516 - t440 * t80;
t341 = t298 * t475 + t405;
t339 = t303 * t344;
t336 = -t1 * t296 - t2 * t302 + (-t22 * t302 - t296 * t425) * qJD(7);
t335 = (-t400 + t211) * qJD(4) - t570;
t333 = t305 * t410 + t409;
t330 = -t299 * t372 - t442;
t329 = -t527 + t298 * t37 + (t112 * t298 + t113 * t653) * qJD(5);
t328 = (-qJD(3) * t250 - t207) * t304 + t666;
t327 = t299 * t332;
t323 = t332 * t653;
t321 = t300 * t324;
t320 = qJD(5) * t323;
t319 = t253 * t372 - t361;
t316 = -t253 * t543 + t353;
t315 = t300 * t317;
t313 = (-t387 - t555) * t299 + (qJD(4) * t678 - t570) * t303;
t31 = -t159 * t548 + t652 * t59 + (t297 * t330 - t301 * t418) * pkin(2);
t312 = (-qJD(3) * qJD(4) + (t301 * t669 - t544) * qJD(1) - t568) * t300 + (0.2e1 * qJD(2) * t669 - t569) * t304;
t311 = -qJD(2) * t346 - t354 * t583 - t361;
t309 = -t298 * t317 - t320;
t308 = t299 * t314 + t323 * t553 - t589 * t659;
t306 = qJD(2) ^ 2;
t242 = t297 * t514 + t498;
t241 = t296 * t515 - t512;
t221 = t250 * pkin(2) * t492;
t199 = -t297 * t248 + t443;
t198 = -t247 * t296 - t520;
t173 = t553 * t583 - t303 * t562 + (-t544 - t373) * t299;
t150 = t201 * t302 + t243 * t296;
t111 = t152 * t302 - t202 * t296;
t97 = -qJD(5) * t445 + t174 * t653 + (-qJD(5) * t245 - t372) * t298;
t96 = qJD(2) * t446 - qJD(5) * t522 + t174 * t298 + t245 * t477 + t480 * t583;
t72 = -pkin(4) * t499 + t113;
t71 = -t165 * pkin(4) - t112 * t652;
t52 = t112 * t628;
t51 = t296 * t649 + t302 * t80;
t50 = -t296 * t80 + t302 * t649;
t47 = qJD(6) * t397 + t173 * t297 + t652 * t97;
t46 = qJD(6) * t152 - t173 * t652 + t97 * t297;
t30 = -t296 * t72 + t302 * t71;
t29 = -t296 * t71 - t302 * t72;
t24 = t28 * t242;
t23 = t28 * t199;
t20 = t28 * t397;
t19 = t47 * pkin(4) + t60;
t18 = -t96 * pkin(4) + t31;
t16 = -qJD(7) * t110 - t296 * t96 + t302 * t47;
t15 = t296 * t47 - t202 * t547 + (qJD(7) * t152 + t96) * t302;
t4 = -qJD(7) * t49 - t18 * t296 - t19 * t302;
t3 = qJD(7) * t48 + t18 * t302 - t19 * t296;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t431, -0.2e1 * t566 * t541, t306 * t305, -0.2e1 * t431, -t306 * t301, 0, 0, 0, 0, 0, -t304 * t377 + t395 * t584 + t352, t253 * t373 + t300 * t377 + (t378 + t661) * t301, -t369 - t568 * t305 + (t291 * t301 + (t556 - 0.2e1 * t656) * t304) * qJD(2), t319, -t304 * t391 - t569 * t305 + (t253 * t301 + t300 * t656) * qJD(2) (-qJD(3) - 0.2e1 * t563) * t562, -t250 * t562 + t208 * t305 + (t561 * t677 + t369) * pkin(2), -t251 * t562 + t207 * t305 + (qJD(3) * t393 - t291 * t677) * pkin(2), t365 + (qJD(3) * t421 - t208 * t304 + t602) * t301 + (t373 * t392 + t319 + t352) * pkin(2), t221 + (t365 + (t602 + (qJD(3) * t251 - t208) * t304) * t301) * pkin(2), t211 * t174 - t245 * t324, -t211 * t173 - t174 * t209 + t244 * t324 - t245 * t317, t174 * t543 - t211 * t372 + t245 * t404 + t301 * t321, t209 * t173 + t671, -t173 * t543 + t209 * t372 - t244 * t404 + t301 * t315, t311, -t299 * t524 + t250 * t173 + t208 * t244 + t330 * t251 + (t301 * t339 - t347 * t583 + (-t352 - t311) * t299 + t373 * t209) * pkin(2), -t303 * t524 + t250 * t174 + t208 * t245 + t331 * t251 + (-t345 * t584 - t301 * t363 + (t372 * t543 + t404 * t584 - t352) * t303 + t373 * t211) * pkin(2) -(-t244 * t303 + t245 * t299) * t207 + (t303 * t173 - t299 * t174 + (-t244 * t299 - t245 * t303) * qJD(4)) * t251 + (t313 * t584 + t372 * t423) * pkin(2), t221 + ((t300 * t462 - t596) * t560 + (t300 * t463 + (qJD(3) * t462 - t208) * t304) * t301) * pkin(2), t165 * t97 - t203 * t66, -t163 * t97 - t165 * t96 + t202 * t66 - t203 * t67, t165 * t173 + t203 * t317 - t66 * t244 + t332 * t97, t163 * t96 + t634, -t163 * t173 - t202 * t317 - t67 * t244 - t332 * t96, t173 * t332 + t671, -t60 * t383 - t663 * t570 - t38 * t244 - t112 * t173 + (-t163 * t539 - t251 * t202 - t392 * t663) * t553 + (-t207 * t202 - t251 * t96 + t60 * t392 + t663 * t386 + (-t163 * t372 - t584 * t67) * pkin(2)) * t299, -t59 * t383 + t159 * t570 - t37 * t244 - t113 * t173 + (t159 * t392 - t165 * t539 - t251 * t203) * t553 + (-t207 * t203 - t251 * t97 + t59 * t392 - t159 * t386 + (-t165 * t372 + t584 * t66) * pkin(2)) * t299, t112 * t97 - t113 * t96 - t159 * t67 - t163 * t59 + t165 * t60 - t202 * t37 + t203 * t38 + t66 * t663, t113 * t59 + t159 * t37 + (t292 * t524 + (t292 * t372 + t442 * t673) * t251) * pkin(2) + t426, t116 * t47 - t152 * t27, -t114 * t47 - t116 * t46 - t152 * t28 - t27 * t397, t116 * t96 + t152 * t67 + t153 * t47 - t202 * t27, t114 * t46 - t20, -t114 * t96 - t153 * t46 - t202 * t28 + t397 * t67, t153 * t96 + t634, t112 * t46 + t114 * t60 + t129 * t67 + t14 * t202 + t153 * t32 - t28 * t663 - t38 * t397 + t80 * t96, t112 * t47 + t116 * t60 - t13 * t202 - t130 * t67 + t152 * t38 - t153 * t31 + t27 * t663 + t380 * t96, -t114 * t31 - t116 * t32 + t129 * t27 + t13 * t397 - t130 * t28 - t14 * t152 + t380 * t46 - t47 * t80, t13 * t130 - t31 * t380 + t426 + t642, -t111 * t5 + t16 * t70, t110 * t5 - t111 * t6 - t15 * t70 - t16 * t68, -t111 * t28 + t16 * t542 - t397 * t5 - t46 * t70, t110 * t6 + t15 * t68, t110 * t28 - t15 * t542 - t397 * t6 + t46 * t68, -t46 * t542 - t20, t110 * t14 + t129 * t6 + t15 * t80 + t2 * t397 - t28 * t48 + t32 * t68 + t4 * t542 + t425 * t46, -t1 * t397 + t111 * t14 - t129 * t5 + t16 * t80 + t22 * t46 + t28 * t49 - t3 * t542 + t32 * t70, -t1 * t110 - t111 * t2 - t15 * t22 + t16 * t425 - t3 * t68 - t4 * t70 + t48 * t5 - t49 * t6, t1 * t49 + t2 * t48 + t22 * t3 - t4 * t425 + t642; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t521, t566 * t307, 0, t521, 0, 0, 0, 0, 0, 0, t304 * t374 + t661 (t253 * qJD(3) - t454 * t563 + t568) * t304 + (-t374 + (-t291 - t392) * qJD(3) + t569) * t300, -t304 * qJD(3) ^ 2 + (-0.2e1 * t556 + t656) * t564, -t253 * t669 - t378, t300 * t333 + t454 * t565, t301 * t410, t250 * t565 + (-t300 * t469 + t304 * t333) * pkin(2), t251 * t565 + (-t304 * t469 + (-qJD(1) * t413 - t409) * t300) * pkin(2), t691 * t304 - t666 + ((t291 * qJD(3) + (-t393 + t482) * qJD(1) + t569) * t304 + (t568 + t449 * (-t253 - t561)) * t300) * pkin(2), pkin(2) * t351 + t519, -t211 * t697 + t303 * t321, t235 * t209 + t211 * t234 + (t597 + t599) * t557 + ((t387 - t555) * t299 + (0.2e1 * t364 + t335) * t303) * t300, t211 * t669 - t235 * t543 + t299 * t344 - t303 * t674 - t363, -t209 * t360 - t299 * t315, -t209 * t669 + t234 * t543 + t299 * t674 + t347 * t304 + t339, qJD(1) * t346 + t300 * t354 - t378, -t250 * t234 - t421 * t553 + t328 * t299 + (-t303 * t448 * t552 + t300 * (-t291 * t553 - t570) + t312 * t299 + t698 * t209) * pkin(2), -t250 * t235 + t421 * t554 + t328 * t303 + (t300 * t387 + t303 * t312 + (t211 * t449 - t345) * t304) * pkin(2) (-t234 * t303 + t235 * t299) * t251 + (t304 * t313 - t423 * t669) * pkin(2) (t292 * t601 + t294 * t351) * pkin(2) + t519 + t572, t165 * t574 + t248 * t66, -t163 * t574 + t165 * t573 - t243 * t66 + t248 * t67, -t165 * t360 - t248 * t317 + t332 * t574 + t585 * t66, -t163 * t573 - t633, t163 * t360 + t243 * t317 + t332 * t573 + t585 * t67, -t327 * t557 - t332 * t234 + (-t690 * t299 - t348) * t300, -t662 * t570 + t112 * t234 + (t112 * t300 - t163 * t651 + t251 * t243 - t392 * t662) * t553 + (t207 * t243 + t662 * t386 + t38 * t300 + t112 * t557 + t573 * t251 + (t163 * t669 - t304 * t67) * pkin(2) + t577 * t392) * t299 - t577 * t383, t196 * t570 + t113 * t234 + (t113 * t300 - t165 * t651 + t196 * t392 + t251 * t248) * t553 + (t207 * t248 - t196 * t386 + t37 * t300 + t113 * t557 - t574 * t251 + (t165 * t669 + t304 * t66) * pkin(2) + t578 * t392) * t299 - t578 * t383, t112 * t574 + t113 * t573 - t163 * t578 + t165 * t577 - t196 * t67 + t243 * t37 - t248 * t38 + t66 * t662, t196 * t37 + t578 * t113 + (t207 * t292 + t488 * t673) * t651 + t390 + t572, -t116 * t615 - t201 * t27, t114 * t615 + t116 * t616 + t199 * t27 - t201 * t28, -t116 * t573 - t153 * t615 + t201 * t67 + t243 * t27, -t114 * t616 + t23, t114 * t573 + t153 * t616 - t199 * t67 + t243 * t28, -t153 * t573 - t633, -t112 * t616 + t114 * t577 - t14 * t243 - t153 * t618 + t175 * t67 + t199 * t38 - t28 * t662 - t573 * t80, -t112 * t615 + t116 * t577 + t13 * t243 + t153 * t617 - t176 * t67 + t201 * t38 + t27 * t662 - t380 * t573, t114 * t617 + t116 * t618 - t13 * t199 - t14 * t201 + t175 * t27 - t176 * t28 - t380 * t616 + t615 * t80, t13 * t176 + t380 * t617 + t390 + t441, -t150 * t5 - t621 * t70, t149 * t5 - t150 * t6 + t621 * t68 - t639 * t70, -t150 * t28 + t199 * t5 - t542 * t621 + t616 * t70, t149 * t6 + t639 * t68, t149 * t28 + t199 * t6 - t542 * t639 - t616 * t68, t542 * t616 + t23, t14 * t149 + t175 * t6 - t199 * t2 - t28 * t75 - t425 * t616 + t542 * t643 - t618 * t68 + t639 * t80, t1 * t199 + t14 * t150 - t175 * t5 - t22 * t616 + t28 * t76 + t542 * t654 - t618 * t70 - t621 * t80, -t1 * t149 - t150 * t2 - t22 * t639 - t425 * t621 + t5 * t75 - t6 * t76 - t643 * t70 + t654 * t68, t1 * t76 + t2 * t75 - t22 * t654 - t425 * t643 + t441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t665, -t253 ^ 2 + t392 ^ 2, t253 * t449 + t386, t665, -t392 * t449 + t404, -t469, -t385, t691, 0, 0, t299 * t324 - t303 * t680 (t571 - t667 + t681) * t303 + (-t211 * t253 - t276 + t335 + t364) * t299, -t211 * t392 - t299 * t404 - t303 * t316, -t209 * t459 + t668, t209 * t392 + t299 * t316 - t303 * t404, -t543 * t392, t208 * t303 + (-t209 - t670) * t251, -t208 * t299 + t251 * t678, t250 * t423 + t463 (-0.1e1 + t567) * t251 * t250, -t692 * t165 + t66 * t514, t178 * t163 + t165 * t177 + (t506 + t603) * t553 + (t528 - t624 + (t505 - t605) * qJD(5)) * t299, -t165 * t459 - t178 * t332 - t66 * t303 - t308, -t163 * t350 - t535, t177 * t332 + t348 * t298 - t67 * t303 + (t543 * t163 + t690 * t298 + t320) * t299, -t327 * t543 + t668, -t38 * t303 + t123 * t383 + (t298 * t207 + t251 * t477) * t292 + (-t250 * t163 + t251 * t177 - t123 * t392 - t112 * t253 + (0.2e1 * t298 * t581 + t112) * qJD(4)) * t299 + t308 * pkin(3), -t37 * t303 + t124 * t383 + (t207 * t653 - t251 * t550) * t292 + (-t250 * t165 + t251 * t178 - t124 * t392 - t113 * t253 + (0.2e1 * t251 * t510 + t113) * qJD(4)) * t299 + (-t659 * t514 + (-t299 * t317 - t332 * t553) * t298) * pkin(3), -t112 * t178 + t113 * t177 - t123 * t165 + t124 * t163 + ((-t505 - t605) * pkin(3) - t396) * t553 + ((t529 - t298 * t67 + (-t506 + t603) * qJD(5)) * pkin(3) + t329) * t299, t250 * t595 - t610 - t113 * t124 + (t299 * t329 - t396 * t553) * pkin(3), -t116 * t575 + t247 * t27, t114 * t575 + t116 * t576 - t242 * t27 + t247 * t28, -t116 * t350 - t153 * t575 - t247 * t67 + t27 * t589, -t114 * t576 - t24, t114 * t350 + t153 * t576 + t242 * t67 + t28 * t589, -t153 * t350 - t535, -t14 * t589 - t101 * t153 - t123 * t114 - t38 * t242 - t576 * t112 - t350 * t80 + ((-t114 * t653 - t153 * t593) * t553 + (-t153 * t436 - t532 + (qJD(5) * t114 + t382) * t298) * t299) * pkin(3), t13 * t589 + t102 * t153 - t123 * t116 - t38 * t247 - t575 * t112 - t350 * t380 + ((-t298 * t500 - t508) * t553 + (-t153 * t411 + t533 + (qJD(5) * t116 - t381) * t298) * t299) * pkin(3), t101 * t116 + t102 * t114 + t13 * t242 + t14 * t247 - t576 * t380 + t575 * t80 + ((-t420 + t444) * t549 + ((t116 * t297 - t502) * t553 + (-t26 - t632 + (t501 + t608) * qJD(6)) * t299) * t298) * pkin(3), -t623 + t380 * t102 - t610 + ((t298 * t403 - t509) * t553 + (-t527 + (-t380 * t440 - t516 * t80) * qJD(5) + (t534 + t551 - t14 * t297 + (t297 * t380 - t652 * t80) * qJD(6)) * t298) * t299) * pkin(3), -t200 * t5 - t637 * t70, t198 * t5 - t200 * t6 + t637 * t68 + t638 * t70, -t200 * t28 - t242 * t5 - t542 * t637 + t576 * t70, t198 * t6 - t638 * t68, t198 * t28 - t242 * t6 + t542 * t638 - t576 * t68, t542 * t576 - t24, -t101 * t68 + t14 * t198 - t156 * t28 + t2 * t242 - t638 * t80 - t576 * t425 - t641 * t542 + (-t68 * t419 + (-t68 * t438 + (-t299 * t6 - t553 * t68) * t298) * t297) * pkin(3), -t1 * t242 - t101 * t70 + t14 * t200 - t422 * t28 - t637 * t80 - t576 * t22 + t640 * t542 + (-t70 * t419 + (-t70 * t438 + (t299 * t5 - t553 * t70) * t298) * t297) * pkin(3), -t1 * t198 + t156 * t5 - t2 * t200 + t22 * t638 + t422 * t6 - t425 * t637 + t640 * t68 + t641 * t70, -t1 * t422 - t623 + t2 * t156 - t640 * t22 + t641 * t425 + (-t80 * t419 + (-t80 * t438 + (-t14 * t299 - t553 * t80) * t298) * t297) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t598, -t209 ^ 2 + t211 ^ 2, -t324 + t681, -t598, -t317 + t680, t404, -t250 * t211 + t253 * t581 + t586, t250 * t209 - t253 * t594 + t582, 0, 0, t165 * t368 - t624, -t529 - t368 * t163 + (-t165 * t332 - t67) * t298, -t165 * t211 + t209 * t323 - t309, t163 * t460 - t528, t314 + t163 * t211 + (-t209 * t332 - t659) * t298, -t332 * t211, t207 * t514 - t146 * t332 + t112 * t211 + t309 * pkin(3) + ((t474 + t163) * t303 - t332 * t589) * t251, -t298 * t586 + t147 * t332 + t113 * t211 + (t298 * t659 - t314) * pkin(3) + ((-qJD(4) * t298 + t165) * t303 - t368 * t299) * t251, t531 + t146 * t165 + t147 * t163 + t368 * t112 + (t165 * t477 - t528) * pkin(3) + (t38 - t332 * t113 + (qJD(5) * t163 - t66) * pkin(3)) * t298, -t303 * t251 ^ 2 * t299 + t609 - t113 * t147 + (qJD(5) * t396 + t531 + t625) * pkin(3), -t116 * t693 - t27 * t515, t133 * t114 + t116 * t132 + (-t420 - t444) * qJD(5) + (-t26 + t632 + (-t501 + t608) * qJD(6)) * t298, t533 + t384 * t153 + (t116 * t332 + t381) * t298, t114 * t341 + t536, t532 - t405 * t153 + (-t114 * t332 + t382) * t298, t153 * t460 - t528, -t14 * t653 + t146 * t114 - t117 * t153 + t405 * t112 + (t114 * t477 - t153 * t412 - t516 * t67) * pkin(3) + (t112 * t475 + t626 + t332 * t80 + (qJD(5) * t458 + t28) * pkin(3)) * t298, t13 * t653 + t146 * t116 + t118 * t153 + t384 * t112 + (t116 * t477 + t153 * t435 - t440 * t67) * pkin(3) + (-t112 * t548 + t530 + t332 * t380 + (t153 * t476 - t27) * pkin(3)) * t298, t118 * t114 + t117 * t116 - t380 * t132 + t80 * t133 + t343 * qJD(5) + (-qJD(6) * t403 - t13 * t297 - t14 * t652) * t298 + (t116 * t412 + t114 * t433 - t28 * t440 + (t114 * t473 - t116 * t550 - t533) * t297) * pkin(3), t609 + t380 * t118 + (t13 * t440 - t447 + t625 + t343 * qJD(6) + (t380 * t515 + t509) * qJD(5)) * pkin(3) + t468, -t246 * t5 - t614 * t70, t241 * t5 - t246 * t6 - t613 * t70 + t614 * t68, -t246 * t28 - t341 * t70 + t5 * t593 - t542 * t614, t241 * t6 + t613 * t68, t241 * t28 + t341 * t68 - t542 * t613 + t593 * t6, -t341 * t542 + t536, -t117 * t68 - t425 * t132 + t14 * t241 - t193 * t28 + t613 * t80 + t619 * t542 + (-t424 * t68 + t425 * t515) * qJD(6) + (t425 * t477 - t2 * t298 + (t550 * t68 - t6 * t653) * pkin(3)) * t297, -t117 * t70 - t22 * t132 + t14 * t246 + t195 * t28 - t614 * t80 - t620 * t542 + (t22 * t515 - t424 * t70) * qJD(6) + (t22 * t477 + t1 * t298 + (t5 * t653 + t550 * t70) * pkin(3)) * t297, -t1 * t241 + t193 * t5 - t195 * t6 - t2 * t246 - t22 * t613 - t425 * t614 - t619 * t70 - t620 * t68, t1 * t195 + t2 * t193 + t620 * t22 - t619 * t425 + (-t412 * t80 - t447) * pkin(3) + t468; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t604, -t163 ^ 2 + t165 ^ 2, t163 * t332 - t66, -t604, t165 * t209 - t318, t317, t113 * t209 + t165 * t594 - t466, qJD(5) * t166 - t112 * t332 - t163 * t594 + t180 * t550 - t298 * t310 - t507, 0, 0, t116 * t367 - t632, -t27 * t652 - t367 * t114 + (-t28 - t606) * t297, -t116 * t165 + t153 * t367 + t629, t114 * t458 - t26, -t153 ^ 2 * t297 + t114 * t165 + t526, -t153 * t165, -t113 * t114 - t80 * t165 - t530, -t113 * t116 - t380 * t165 + t626 + (-t500 + t367) * t112, -t112 * t502 - t297 * t658 - t367 * t80 + t534, -t52 + (-t525 - t113) * t112, -t297 * t646 + (-t297 * t547 - t406) * t70, t84 * t68 + t70 * t83 + (-t511 * t68 - t517 * t70) * qJD(6) + (t648 - t645 + (t296 * t68 - t627) * qJD(7)) * t297, -t5 * t652 - t406 * t542 + (-t153 * t70 - t401) * t297, t297 * t647 + (t297 * t546 + t407) * t68, -t6 * t652 - t407 * t542 + (t153 * t68 - t402) * t297, -t458 * t542 - t26, t2 * t652 - t29 * t542 - t80 * t83 + (-t302 * t453 + t517 * t80) * qJD(6) + (pkin(4) * t401 - t112 * t68 + t153 * t425 + t546 * t80 + t636) * t297, -t1 * t652 + t30 * t542 - t80 * t84 + (t296 * t453 + t511 * t80) * qJD(6) + (pkin(4) * t402 - t112 * t70 + t153 * t22 - t547 * t80 + t635) * t297, -t425 * t84 + t22 * t83 + t29 * t70 + t30 * t68 + ((t511 * t70 + t517 * t68) * pkin(4) + t356) * qJD(6) + ((t647 - t646 + (-t296 * t70 + t302 * t68) * qJD(7)) * pkin(4) + t336) * t297, t425 * t29 - t22 * t30 - t52 + (qJD(6) * t356 + t297 * t336) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t607, -t114 ^ 2 + t116 ^ 2, t114 * t153 - t27, -t607, -t28 + t606, t67, t658, t112 * t114 + t153 * t80 - t13, 0, 0, -t542 * t627 + t648 (t5 + t687) * t302 + (t6 + t686) * t296, -t302 * t542 ^ 2 + t116 * t70 + t631, -t457 * t68 + t645, -t116 * t68 + t457 * t542 + t630, t542 * t116, -t116 * t425 + t635 - t68 * t380 - (-pkin(4) * t546 + t50) * t542 + t398 * t296, -t116 * t22 - t636 - t70 * t380 - (pkin(4) * t547 - t51) * t542 + t398 * t302, t50 * t70 + t51 * t68 + ((-qJD(7) * t70 + t6) * pkin(4) + t675) * t302 + ((t5 - t612) * pkin(4) + t676) * t296, t425 * t50 - t22 * t51 - t80 * t380 + (-t1 * t302 + t2 * t296 + (t22 * t296 - t302 * t425) * qJD(7)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t644, -t68 ^ 2 + t70 ^ 2, -t5 + t687, -t644, -t6 + t686, -t28, -t70 * t80 + t676, t68 * t80 + t675, 0, 0;];
tauc_reg  = t9;
