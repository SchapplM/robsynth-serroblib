% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:03:27
% EndTime: 2019-03-10 01:03:41
% DurationCPUTime: 9.34s
% Computational Cost: add. (14101->545), mult. (35328->686), div. (0->0), fcn. (26301->8), ass. (0->247)
t573 = cos(qJ(5));
t662 = qJD(5) * t573;
t571 = sin(qJ(3));
t572 = sin(qJ(2));
t670 = qJD(1) * t572;
t643 = t571 * t670;
t574 = cos(qJ(3));
t575 = cos(qJ(2));
t669 = qJD(1) * t575;
t644 = t574 * t669;
t526 = -t643 + t644;
t527 = -t571 * t669 - t574 * t670;
t570 = sin(qJ(4));
t713 = cos(qJ(4));
t721 = t713 * t526 + t570 * t527;
t735 = t721 * t573;
t743 = t662 - t735;
t569 = sin(qJ(5));
t714 = pkin(7) + pkin(8);
t549 = t714 * t575;
t543 = qJD(1) * t549;
t532 = t574 * t543;
t548 = t714 * t572;
t541 = qJD(1) * t548;
t629 = t541 * t571 - t532;
t709 = pkin(9) * t526;
t486 = t629 - t709;
t523 = t527 * pkin(9);
t528 = t571 * t543;
t672 = -t574 * t541 - t528;
t487 = t523 + t672;
t451 = t570 * t486 + t487 * t713;
t561 = pkin(2) * t574 + pkin(3);
t642 = qJD(4) * t713;
t665 = qJD(4) * t570;
t688 = t570 * t571;
t503 = t561 * t642 + (-t571 * t665 + (t574 * t713 - t688) * qJD(3)) * pkin(2);
t734 = t503 - t451;
t742 = t569 * t734;
t605 = -t570 * t526 + t527 * t713;
t663 = qJD(5) * t569;
t566 = qJD(2) + qJD(3);
t687 = t571 * t575;
t540 = t572 * t574 + t687;
t716 = qJD(1) * t540;
t582 = t566 * t716;
t581 = t570 * t582;
t656 = qJD(1) * qJD(2);
t641 = t575 * t656;
t505 = qJD(3) * t644 - t566 * t643 + t574 * t641;
t647 = t713 * t505;
t579 = t647 - t581;
t578 = qJD(4) * t721 + t579;
t639 = qJD(4) + t566;
t733 = -qJD(5) * t639 - t578;
t419 = t573 * t733 - t605 * t663;
t417 = t419 * t569;
t483 = t569 * t639 - t573 * t605;
t420 = t483 * qJD(5) + t578 * t569;
t633 = t570 * t505 + t713 * t582;
t433 = -t605 * qJD(4) + t633;
t430 = t573 * t433;
t481 = -t569 * t605 - t573 * t639;
t428 = t569 * t433;
t732 = -qJD(5) + t721;
t679 = -t662 * t732 + t428;
t740 = t732 * t569;
t741 = -t721 ^ 2 * MDP(19) + (-t721 * t566 + t579) * MDP(20) - t633 * MDP(21) + (MDP(18) * t721 + MDP(19) * t605 - MDP(21) * t566 - MDP(29) * t732) * t605 + (t743 * t483 - t417) * MDP(25) + (t483 * t605 + t732 * t735 + t679) * MDP(27) + (-t419 * t573 - t569 * t420 - t743 * t481 + t483 * t740) * MDP(26) + (-t481 * t605 - t732 * t740 + t430) * MDP(28);
t534 = qJD(2) * pkin(2) - t541;
t632 = t574 * t534 - t528;
t479 = t523 + t632;
t471 = pkin(3) * t566 + t479;
t612 = -t534 * t571 - t532;
t480 = -t612 + t709;
t477 = t713 * t480;
t436 = t570 * t471 + t477;
t432 = pkin(10) * t639 + t436;
t562 = -pkin(2) * t575 - pkin(1);
t547 = t562 * qJD(1);
t512 = -pkin(3) * t526 + t547;
t448 = -pkin(4) * t721 + pkin(10) * t605 + t512;
t407 = -t432 * t569 + t448 * t573;
t657 = qJD(6) - t407;
t395 = pkin(5) * t732 + t657;
t563 = pkin(2) * t670;
t556 = qJD(2) * t563;
t587 = t526 * t642 + t527 * t665 + t647;
t666 = qJD(3) * t574;
t667 = qJD(3) * t571;
t668 = qJD(2) * t572;
t390 = t433 * pkin(4) - t587 * pkin(10) + t556 + (pkin(10) * t570 + pkin(3)) * qJD(1) * (qJD(2) * t687 + t572 * t666 + t574 * t668 + t575 * t667);
t651 = qJD(2) * t714;
t623 = qJD(1) * t651;
t536 = t575 * t623;
t630 = -t571 * t536 - t543 * t667;
t535 = t572 * t623;
t724 = t574 * (qJD(3) * t534 - t535);
t442 = -pkin(9) * t582 + t630 + t724;
t631 = t571 * t535 - t574 * t536;
t591 = qJD(3) * t612 + t631;
t443 = -pkin(9) * t505 + t591;
t588 = -t442 * t713 - t570 * t443 - t471 * t642 + t480 * t665;
t597 = t569 * t390 - t432 * t663 + t448 * t662 - t573 * t588;
t707 = qJ(6) * t433;
t371 = -qJD(6) * t732 + t597 + t707;
t624 = -t573 * t390 + t432 * t662 + t448 * t663 - t569 * t588;
t711 = pkin(5) * t433;
t373 = t624 - t711;
t723 = t371 * t573 + t373 * t569;
t739 = t743 * t395 + t723;
t710 = pkin(5) * t605;
t492 = t605 * qJ(6);
t736 = t721 * t569;
t383 = t570 * t442 - t713 * t443 + t471 * t665 + t480 * t642;
t599 = t512 * t605 - t383;
t376 = pkin(5) * t420 + qJ(6) * t419 - qJD(6) * t483 + t383;
t476 = t570 * t480;
t435 = t471 * t713 - t476;
t431 = -pkin(4) * t639 - t435;
t410 = t481 * pkin(5) - t483 * qJ(6) + t431;
t610 = -t376 * t573 - t395 * t605 + t410 * t663;
t408 = t432 * t573 + t448 * t569;
t396 = -qJ(6) * t732 + t408;
t706 = t376 * t569;
t731 = t396 * t605 - t706;
t608 = -t383 * t573 + t407 * t605 + t431 * t663;
t620 = t383 * t569 - t408 * t605 + t431 * t662;
t584 = -t512 * t721 + t588;
t466 = -pkin(4) * t605 - pkin(10) * t721;
t727 = -0.2e1 * t656;
t726 = MDP(4) * t572;
t725 = MDP(5) * (t572 ^ 2 - t575 ^ 2);
t444 = t570 * t479 + t477;
t622 = pkin(3) * t665 - t444;
t539 = t571 * t572 - t574 * t575;
t508 = -t570 * t539 + t540 * t713;
t515 = pkin(3) * t539 + t562;
t604 = -t539 * t713 - t570 * t540;
t462 = -pkin(4) * t604 - pkin(10) * t508 + t515;
t493 = -pkin(9) * t540 - t548 * t574 - t549 * t571;
t611 = t548 * t571 - t549 * t574;
t494 = -pkin(9) * t539 - t611;
t464 = t570 * t493 + t494 * t713;
t675 = t569 * t462 + t573 * t464;
t645 = t713 * t571;
t673 = t486 * t713 - t570 * t487 + t561 * t665 + (t571 * t642 + (t570 * t574 + t645) * qJD(3)) * pkin(2);
t720 = pkin(2) * t688 - t713 * t561;
t719 = -t581 + t587;
t522 = pkin(5) * t663 - qJ(6) * t662 - t569 * qJD(6);
t718 = pkin(5) * t736 - qJ(6) * t735 - t522;
t717 = t713 * t493 - t570 * t494;
t715 = t483 ^ 2;
t712 = pkin(3) * t527;
t708 = pkin(3) * qJD(4);
t705 = t408 * t732;
t704 = t410 * t721;
t703 = t420 * t573;
t702 = t431 * t721;
t521 = pkin(2) * t645 + t570 * t561 + pkin(10);
t701 = t433 * t521;
t559 = pkin(3) * t570 + pkin(10);
t700 = t433 * t559;
t699 = t481 * t569;
t698 = t483 * t481;
t697 = t483 * t573;
t692 = t503 * t573;
t691 = t508 * t573;
t690 = t547 * t527;
t576 = qJD(2) ^ 2;
t686 = t572 * t576;
t685 = t575 * t576;
t577 = qJD(1) ^ 2;
t684 = t575 * t577;
t618 = pkin(5) * t569 - qJ(6) * t573;
t682 = t618 * t721 + t436 - t522;
t681 = t718 - t622;
t680 = t718 - t673;
t678 = t573 * t435 + t569 * t466;
t445 = t479 * t713 - t476;
t456 = t466 - t712;
t677 = t573 * t445 + t569 * t456;
t452 = t456 + t563;
t676 = t573 * t451 + t569 * t452;
t664 = qJD(5) * t410;
t655 = t713 * pkin(3);
t565 = pkin(2) * t668;
t650 = t569 * t713;
t649 = t573 * t713;
t640 = -pkin(2) * t566 - t534;
t595 = t540 * qJD(3);
t510 = t540 * qJD(2) + t595;
t502 = pkin(3) * t510 + t565;
t638 = pkin(1) * t727;
t401 = -t492 + t676;
t635 = -t401 + t692;
t634 = t452 * t573 - t710 + t742;
t625 = pkin(3) * t642;
t619 = t573 * pkin(5) + t569 * qJ(6);
t616 = t395 * t573 - t396 * t569;
t615 = -t701 - t702;
t614 = -t700 - t702;
t613 = -t435 * t569 + t466 * t573;
t545 = -pkin(4) - t619;
t609 = t410 * t735 + t731;
t607 = -t547 * t526 - t630;
t509 = t566 * t539;
t453 = qJD(4) * t604 - t509 * t713 - t570 * t510;
t603 = t453 * t569 + t508 * t662;
t602 = -t453 * t573 + t508 * t663;
t601 = -t521 * t663 + t692;
t600 = t410 * t483 + t624;
t542 = t572 * t651;
t544 = t575 * t651;
t598 = -t574 * t542 - t571 * t544 - t548 * t666 - t549 * t667;
t459 = -pkin(9) * t510 + t598;
t590 = qJD(3) * t611 + t542 * t571 - t574 * t544;
t460 = pkin(9) * t509 + t590;
t393 = t717 * qJD(4) + t713 * t459 + t570 * t460;
t454 = qJD(4) * t508 - t570 * t509 + t510 * t713;
t404 = pkin(4) * t454 - pkin(10) * t453 + t502;
t596 = t573 * t393 + t569 * t404 + t462 * t662 - t464 * t663;
t594 = t679 * pkin(10);
t592 = -t559 * t663 + t573 * t625;
t589 = (-t663 + t736) * t396 + t739;
t586 = qJD(5) * t616 + t723;
t585 = -t417 - t703 + (t697 + t699) * qJD(5);
t394 = qJD(4) * t464 + t570 * t459 - t460 * t713;
t580 = t527 * t526 * MDP(11) + (-t526 * t566 + t505) * MDP(13) + (-t527 * t566 - t582) * MDP(14) + (-t526 ^ 2 + t527 ^ 2) * MDP(12) + t741;
t560 = -t655 - pkin(4);
t537 = -t655 + t545;
t520 = -pkin(4) + t720;
t513 = t563 - t712;
t511 = t545 + t720;
t488 = pkin(3) * t582 + t556;
t449 = pkin(5) * t483 + qJ(6) * t481;
t423 = t508 * t618 - t717;
t412 = pkin(5) * t604 - t462 * t573 + t464 * t569;
t411 = -qJ(6) * t604 + t675;
t406 = -t613 + t710;
t405 = t678 - t492;
t400 = t445 * t569 - t456 * t573 + t710;
t399 = -t492 + t677;
t398 = -t481 * t732 - t419;
t378 = t618 * t453 + (qJD(5) * t619 - qJD(6) * t573) * t508 + t394;
t375 = -pkin(5) * t454 + t675 * qJD(5) + t393 * t569 - t404 * t573;
t374 = qJ(6) * t454 - qJD(6) * t604 + t596;
t1 = [(-t509 * MDP(13) - t510 * MDP(14) + MDP(16) * t590 - MDP(17) * t598) * t566 + (t562 * t505 - t547 * t509 + (-t527 + t716) * t565) * MDP(17) + (t505 * t540 + t509 * t527) * MDP(11) + (-t505 * t539 - t509 * t526 + t527 * t510 - t540 * t582) * MDP(12) + (t512 * t453 + t488 * t508 - t502 * t605 + t515 * t578) * MDP(24) + (-t453 * t605 + t719 * t508) * MDP(18) + (t515 * t433 + t512 * t454 - t488 * t604 - t502 * t721) * MDP(23) + (-t508 * t433 + t453 * t721 + t454 * t605 + t604 * t719) * MDP(19) + (-t526 * t565 + t547 * t510 + (t562 * t595 + (t572 * pkin(2) * t539 + t540 * t562) * qJD(2)) * qJD(1)) * MDP(16) + (t373 * t604 + t375 * t732 + t378 * t481 - t395 * t454 + t410 * t603 - t412 * t433 + t420 * t423 + t508 * t706) * MDP(32) + (-t371 * t604 - t374 * t732 - t376 * t691 - t378 * t483 + t396 * t454 + t410 * t602 + t411 * t433 + t419 * t423) * MDP(34) + (t420 * t604 - t428 * t508 - t454 * t481 + t603 * t732) * MDP(28) + (t419 * t604 + t430 * t508 + t454 * t483 + t602 * t732) * MDP(27) + (-t433 * t604 - t454 * t732) * MDP(29) + (t624 * t604 + t407 * t454 + t394 * t481 - t717 * t420 + (-(-qJD(5) * t464 + t404) * t732 + t462 * t433 + t431 * qJD(5) * t508) * t573 + (-(-qJD(5) * t462 - t393) * t732 - t464 * t433 + t383 * t508 + t431 * t453) * t569) * MDP(30) + (t383 * t691 + t394 * t483 - t408 * t454 + t419 * t717 - t431 * t602 - t433 * t675 + t596 * t732 + t597 * t604) * MDP(31) + (t453 * MDP(20) - t454 * MDP(21) - t394 * MDP(23) - t393 * MDP(24)) * t639 + ((-t481 * t573 - t483 * t569) * t453 + (t417 - t703 + (-t697 + t699) * qJD(5)) * t508) * MDP(26) + (-t419 * t691 - t483 * t602) * MDP(25) + (pkin(7) * t686 + t575 * t638) * MDP(10) + (-pkin(7) * t685 + t572 * t638) * MDP(9) - MDP(7) * t686 + 0.2e1 * t641 * t726 + t725 * t727 + MDP(6) * t685 + (t371 * t411 + t373 * t412 + t374 * t396 + t375 * t395 + t376 * t423 + t378 * t410) * MDP(35) + (-t374 * t481 + t375 * t483 - t411 * t420 - t412 * t419 + t616 * t453 + (-t371 * t569 + t373 * t573 + (-t395 * t569 - t396 * t573) * qJD(5)) * t508) * MDP(33); t580 + (t420 * t511 + (-t701 - t704) * t569 - t680 * t481 - (-t521 * t662 - t634) * t732 + t610) * MDP(32) + (t419 * t511 + (-t664 + t701) * t573 + t680 * t483 - (-t401 + t601) * t732 + t609) * MDP(34) + (t376 * t511 + t395 * t634 + t396 * t635 - t410 * t680 + t521 * t586) * MDP(35) + (-t520 * t419 + t615 * t573 + t673 * t483 - (-t601 + t676) * t732 + t620) * MDP(31) + (t520 * t420 + t615 * t569 + t673 * t481 - ((-qJD(5) * t521 - t452) * t573 - t742) * t732 + t608) * MDP(30) + (-t481 * t635 + t483 * t634 + t521 * t585 + t589) * MDP(33) + (t527 * t563 + t672 * t566 + (qJD(3) * t640 + t535) * t574 + t607) * MDP(17) + (t526 * t563 + t690 - t629 * t566 + (t571 * t640 - t532) * qJD(3) + t631) * MDP(16) + t577 * t725 - t684 * t726 + (t513 * t721 + t599) * MDP(23) + (t513 * t605 + t584) * MDP(24) + (-t673 * MDP(23) - MDP(24) * t734) * t639 + (MDP(9) * t572 * t577 + MDP(10) * t684) * pkin(1); t580 + (t376 * t537 - t395 * t400 - t396 * t399 - t681 * t410 + (t395 * t650 + t396 * t649) * t708 + t586 * t559) * MDP(35) + (-t560 * t419 + t614 * t573 + t622 * t483 - (-t592 + t677) * t732 + t620) * MDP(31) + (t399 * t481 - t400 * t483 + (-t481 * t649 + t483 * t650) * t708 + t585 * t559 + t589) * MDP(33) + (t560 * t420 + t614 * t569 + t622 * t481 - ((-qJD(5) * t559 - t456) * t573 + (-t625 + t445) * t569) * t732 + t608) * MDP(30) + (t444 * t639 + (-t527 * t721 - t639 * t665) * pkin(3) + t599) * MDP(23) + (t537 * t419 + (-t664 + t700) * t573 + t681 * t483 - (-t399 + t592) * t732 + t609) * MDP(34) + (t537 * t420 + (-t700 - t704) * t569 - t681 * t481 - (-t559 * t662 - t569 * t625 + t400) * t732 + t610) * MDP(32) + (-t566 * t612 + t591 + t690) * MDP(16) + (t566 * t632 + t607 - t724) * MDP(17) + (t445 * t639 + (-t527 * t605 - t639 * t642) * pkin(3) + t584) * MDP(24); (t436 * t639 + t599) * MDP(23) + (t435 * t639 + t584) * MDP(24) + (-pkin(4) * t420 - t431 * t736 - t436 * t481 + t613 * t732 - t594 + t608) * MDP(30) + (pkin(4) * t419 - t678 * t732 - t436 * t483 - t431 * t735 + (-t663 * t732 - t430) * pkin(10) + t620) * MDP(31) + (-t406 * t732 - t410 * t736 + t420 * t545 - t481 * t682 - t594 + t610) * MDP(32) + (pkin(10) * t585 + t396 * t740 + t405 * t481 - t406 * t483 + t739) * MDP(33) + (t419 * t545 - (-pkin(10) * t663 - t405) * t732 + t682 * t483 + (pkin(10) * t433 + t410 * t732) * t573 + t731) * MDP(34) + (pkin(10) * t586 + t376 * t545 - t395 * t406 - t396 * t405 - t410 * t682) * MDP(35) + t741; MDP(25) * t698 + (-t481 ^ 2 + t715) * MDP(26) + t398 * MDP(27) + (-t483 * t732 + t569 * t733 + t605 * t662) * MDP(28) + t433 * MDP(29) + (-t431 * t483 - t624 - t705) * MDP(30) + (-t407 * t732 + t431 * t481 - t597) * MDP(31) + (-t449 * t481 - t600 - t705 + 0.2e1 * t711) * MDP(32) + (pkin(5) * t419 - qJ(6) * t420 + (t396 - t408) * t483 + (t395 - t657) * t481) * MDP(33) + (0.2e1 * t707 - t410 * t481 + t449 * t483 - (0.2e1 * qJD(6) - t407) * t732 + t597) * MDP(34) + (-pkin(5) * t373 + qJ(6) * t371 - t395 * t408 + t396 * t657 - t410 * t449) * MDP(35); (-t433 + t698) * MDP(32) + t398 * MDP(33) + (-t732 ^ 2 - t715) * MDP(34) + (t396 * t732 + t600 - t711) * MDP(35);];
tauc  = t1;
