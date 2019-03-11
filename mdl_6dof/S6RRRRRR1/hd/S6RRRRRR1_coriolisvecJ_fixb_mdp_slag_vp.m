% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRRR1
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
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:30:05
% EndTime: 2019-03-10 03:30:21
% DurationCPUTime: 9.73s
% Computational Cost: add. (14759->487), mult. (38247->635), div. (0->0), fcn. (30100->10), ass. (0->241)
t614 = cos(qJ(6));
t693 = qJD(6) * t614;
t610 = sin(qJ(5));
t615 = cos(qJ(5));
t617 = cos(qJ(3));
t618 = cos(qJ(2));
t700 = qJD(1) * t618;
t686 = t617 * t700;
t612 = sin(qJ(3));
t613 = sin(qJ(2));
t701 = qJD(1) * t613;
t687 = t612 * t701;
t566 = -t686 + t687;
t568 = -t612 * t700 - t617 * t701;
t611 = sin(qJ(4));
t616 = cos(qJ(4));
t640 = t566 * t611 + t616 * t568;
t641 = -t566 * t616 + t611 * t568;
t496 = t610 * t640 + t615 * t641;
t786 = t496 * t614;
t789 = t693 - t786;
t606 = qJD(2) + qJD(3);
t690 = qJD(1) * qJD(2);
t685 = t618 * t690;
t544 = qJD(3) * t686 - t606 * t687 + t617 * t685;
t580 = t612 * t618 + t613 * t617;
t550 = t606 * t580;
t545 = t550 * qJD(1);
t697 = qJD(4) * t616;
t698 = qJD(4) * t611;
t472 = t616 * t544 - t611 * t545 - t566 * t697 + t568 * t698;
t605 = qJD(4) + t606;
t670 = t544 * t611 + t616 * t545;
t766 = -qJD(4) * t640 + t670;
t609 = sin(qJ(6));
t694 = qJD(6) * t609;
t695 = qJD(5) * t615;
t696 = qJD(5) * t610;
t430 = t615 * t472 - t610 * t766 + t640 * t696 + t641 * t695;
t596 = qJD(5) + t605;
t712 = t614 * t430 + t596 * t693;
t765 = t610 * t641 - t615 * t640;
t418 = -t694 * t765 + t712;
t416 = t418 * t609;
t417 = t418 * t614;
t483 = t596 * t609 + t614 * t765;
t736 = t430 * t609;
t419 = t483 * qJD(6) + t736;
t431 = qJD(5) * t765 + t472 * t610 + t615 * t766;
t429 = t614 * t431;
t727 = t765 * t609;
t481 = -t614 * t596 + t727;
t691 = -qJD(6) + t496;
t427 = t609 * t431;
t713 = -t691 * t693 + t427;
t785 = t691 * t609;
t774 = (t483 * t789 + t416) * MDP(32) + (-t483 * t765 + t691 * t786 + t713) * MDP(34) + (t481 * t765 - t691 * t785 + t429) * MDP(35) + (-t609 * t419 - t481 * t789 + t483 * t785 + t417) * MDP(33) - t431 * MDP(28) - t496 ^ 2 * MDP(26) + (-t496 * t596 + t430) * MDP(27) + (-MDP(25) * t496 + MDP(26) * t765 + MDP(28) * t596 + MDP(36) * t691) * t765;
t788 = t641 * MDP(18) * t640 + (-t605 * t641 + t472) * MDP(20) + (t640 ^ 2 - t641 ^ 2) * MDP(19) + (-t605 * t640 - t766) * MDP(21) + t774;
t742 = pkin(7) + pkin(8);
t688 = qJD(2) * t742;
t655 = qJD(1) * t688;
t577 = t618 * t655;
t589 = t742 * t618;
t583 = qJD(1) * t589;
t699 = qJD(3) * t612;
t667 = -t612 * t577 - t583 * t699;
t588 = t742 * t613;
t581 = qJD(1) * t588;
t737 = qJD(2) * pkin(2);
t575 = -t581 + t737;
t576 = t613 * t655;
t754 = t617 * (qJD(3) * t575 - t576);
t474 = -pkin(9) * t545 + t667 + t754;
t562 = t568 * pkin(9);
t569 = t612 * t583;
t669 = t617 * t575 - t569;
t515 = t562 + t669;
t506 = pkin(3) * t606 + t515;
t502 = t616 * t506;
t530 = pkin(10) * t640;
t573 = t617 * t583;
t639 = -t575 * t612 - t573;
t668 = t612 * t576 - t617 * t577;
t626 = qJD(3) * t639 + t668;
t475 = -pkin(9) * t544 + t626;
t739 = pkin(9) * t566;
t516 = -t639 - t739;
t675 = t611 * t475 - t516 * t698;
t413 = t616 * t474 - t670 * pkin(10) + (t502 + t530) * qJD(4) + t675;
t508 = t611 * t516;
t673 = t502 - t508;
t458 = t673 + t530;
t455 = pkin(4) * t605 + t458;
t510 = t616 * t516;
t646 = -t506 * t611 - t510;
t676 = -t611 * t474 + t616 * t475;
t627 = qJD(4) * t646 + t676;
t414 = -pkin(10) * t472 + t627;
t738 = pkin(10) * t641;
t459 = -t646 + t738;
t678 = t414 * t610 - t459 * t696;
t400 = t615 * (qJD(5) * t455 + t413) + t678;
t602 = -pkin(2) * t618 - pkin(1);
t587 = t602 * qJD(1);
t551 = pkin(3) * t566 + t587;
t500 = -pkin(4) * t641 + t551;
t772 = t496 * t500;
t623 = -t400 - t772;
t732 = t459 * t610;
t434 = t455 * t615 - t732;
t432 = -pkin(5) * t596 - t434;
t787 = t432 * t496;
t731 = t459 * t615;
t435 = t455 * t610 + t731;
t679 = t413 * t610 - t615 * t414;
t401 = qJD(5) * t435 + t679;
t753 = t765 * t500;
t622 = -t401 - t753;
t451 = pkin(5) * t765 - pkin(11) * t496;
t666 = t581 * t612 - t573;
t518 = t666 + t739;
t705 = -t617 * t581 - t569;
t519 = t562 + t705;
t601 = pkin(2) * t617 + pkin(3);
t721 = t611 * t612;
t781 = -t601 * t697 - (-t612 * t698 + (t616 * t617 - t721) * qJD(3)) * pkin(2) + t611 * t518 + t616 * t519;
t719 = t612 * t616;
t780 = -t601 * t698 + (-t612 * t697 + (-t611 * t617 - t719) * qJD(3)) * pkin(2) - t616 * t518 + t519 * t611;
t726 = t551 * t640;
t779 = t627 + t726;
t778 = -t551 * t641 - t675;
t433 = pkin(11) * t596 + t435;
t446 = -pkin(5) * t496 - pkin(11) * t765 + t500;
t411 = t433 * t614 + t446 * t609;
t777 = t401 * t609 + t411 * t765 + t432 * t693;
t648 = t433 * t609 - t446 * t614;
t776 = -t401 * t614 + t432 * t694 + t648 * t765;
t740 = pkin(4) * t640;
t769 = -t738 - t780;
t768 = t530 + t781;
t767 = (-qJD(4) * t506 - t474) * t616 + t778;
t760 = -0.2e1 * t690;
t758 = MDP(4) * t613;
t757 = MDP(5) * (t613 ^ 2 - t618 ^ 2);
t752 = qJD(1) * t580;
t582 = t613 * t688;
t584 = t618 * t688;
t723 = t588 * t617;
t634 = -qJD(3) * t723 - t617 * t582 - t612 * t584 - t589 * t699;
t487 = -pkin(9) * t550 + t634;
t579 = t612 * t613 - t617 * t618;
t549 = t606 * t579;
t638 = t588 * t612 - t589 * t617;
t624 = qJD(3) * t638 + t582 * t612 - t617 * t584;
t488 = pkin(9) * t549 + t624;
t528 = -pkin(9) * t580 - t589 * t612 - t723;
t529 = -pkin(9) * t579 - t638;
t645 = -t528 * t611 - t529 * t616;
t746 = qJD(4) * t645 - t487 * t611 + t616 * t488;
t548 = -t579 * t611 + t580 * t616;
t478 = qJD(4) * t548 - t549 * t611 + t616 * t550;
t633 = t616 * t487 + t611 * t488 + t528 * t697 - t529 * t698;
t422 = -pkin(10) * t478 + t633;
t547 = t616 * t579 + t580 * t611;
t477 = -qJD(4) * t547 - t549 * t616 - t550 * t611;
t423 = -pkin(10) * t477 + t746;
t470 = -pkin(10) * t548 + t528 * t616 - t529 * t611;
t471 = -pkin(10) * t547 - t645;
t647 = t470 * t615 - t471 * t610;
t402 = qJD(5) * t647 + t422 * t615 + t423 * t610;
t498 = t615 * t547 + t548 * t610;
t442 = -qJD(5) * t498 + t477 * t615 - t478 * t610;
t445 = t470 * t610 + t471 * t615;
t499 = -t547 * t610 + t548 * t615;
t554 = pkin(3) * t579 + t602;
t511 = pkin(4) * t547 + t554;
t450 = pkin(5) * t498 - pkin(11) * t499 + t511;
t745 = t401 * t499 - t445 * t431 + t432 * t442 + t691 * (qJD(6) * t450 + t402) - t498 * (qJD(6) * t446 + t400);
t741 = pkin(3) * t568;
t734 = t432 * t499;
t733 = t450 * t431;
t724 = t587 * t568;
t722 = t610 * t611;
t720 = t611 * t615;
t619 = qJD(2) ^ 2;
t718 = t613 * t619;
t717 = t618 * t619;
t620 = qJD(1) ^ 2;
t716 = t618 * t620;
t672 = -t515 * t611 - t510;
t460 = t672 - t738;
t708 = t616 * t515 - t508;
t461 = t530 + t708;
t600 = pkin(3) * t616 + pkin(4);
t711 = t460 * t610 + t461 * t615 - t600 * t695 - (-t611 * t696 + (t615 * t616 - t722) * qJD(4)) * pkin(3);
t560 = -pkin(2) * t721 + t601 * t616 + pkin(4);
t563 = pkin(2) * t719 + t601 * t611;
t643 = t560 * t615 - t563 * t610;
t710 = -qJD(5) * t643 + t610 * t769 + t615 * t768;
t642 = t560 * t610 + t563 * t615;
t709 = qJD(5) * t642 - t610 * t768 + t615 * t769;
t706 = t460 * t615 - t461 * t610 + t600 * t696 + (t611 * t695 + (t610 * t616 + t720) * qJD(4)) * pkin(3);
t604 = t613 * t737;
t603 = pkin(2) * t701;
t684 = -pkin(2) * t606 - t575;
t683 = -pkin(3) * t605 - t506;
t520 = t545 * pkin(3) + qJD(2) * t603;
t538 = pkin(3) * t550 + t604;
t682 = -pkin(4) * t596 - t455;
t681 = pkin(1) * t760;
t503 = -t741 - t740;
t448 = t451 + t503;
t524 = pkin(11) + t642;
t658 = qJD(6) * t524 + t448 + t603;
t561 = pkin(3) * t720 + t600 * t610 + pkin(11);
t657 = qJD(6) * t561 + t448;
t598 = pkin(4) * t610 + pkin(11);
t656 = qJD(6) * t598 + t451 - t740;
t436 = t458 * t610 + t731;
t654 = pkin(4) * t696 - t436;
t437 = t458 * t615 - t732;
t653 = -pkin(4) * t695 + t437;
t651 = -t431 * t598 - t787;
t650 = -t431 * t524 - t787;
t649 = -t431 * t561 - t787;
t464 = pkin(4) * t478 + t538;
t636 = t587 * t566 - t667;
t635 = t442 * t614 - t499 * t694;
t452 = pkin(4) * t766 + t520;
t621 = -t568 * t566 * MDP(11) + t544 * MDP(13) + (-t566 ^ 2 + t568 ^ 2) * MDP(12) + (t566 * MDP(13) + (-t568 - t752) * MDP(14)) * t606 + t788;
t599 = -pkin(4) * t615 - pkin(5);
t559 = pkin(3) * t722 - t600 * t615 - pkin(5);
t552 = t603 - t741;
t523 = -pkin(5) - t643;
t501 = t503 + t603;
t443 = qJD(5) * t499 + t477 * t610 + t615 * t478;
t409 = pkin(5) * t443 - pkin(11) * t442 + t464;
t405 = t431 * pkin(5) - t430 * pkin(11) + t452;
t404 = t614 * t405;
t403 = qJD(5) * t445 + t422 * t610 - t423 * t615;
t1 = [(t417 * t499 + t483 * t635) * MDP(32) + ((-t481 * t614 - t483 * t609) * t442 + (-t416 - t419 * t614 + (t481 * t609 - t483 * t614) * qJD(6)) * t499) * MDP(33) - MDP(7) * t718 + (pkin(7) * t718 + t618 * t681) * MDP(10) + (-pkin(7) * t717 + t613 * t681) * MDP(9) + (t431 * t511 + t443 * t500 + t452 * t498 - t464 * t496) * MDP(30) + (-t430 * t498 - t431 * t499 + t442 * t496 - t443 * t765) * MDP(26) + (-t472 * t547 + t477 * t641 + t478 * t640 - t548 * t766) * MDP(19) + (t430 * t511 + t442 * t500 + t452 * t499 + t464 * t765) * MDP(31) + (t430 * t499 + t442 * t765) * MDP(25) + (t477 * MDP(20) - t478 * MDP(21) + MDP(23) * t746 - t633 * MDP(24)) * t605 + (t472 * t548 - t477 * t640) * MDP(18) + (t554 * t472 + t551 * t477 + t520 * t548 - t538 * t640) * MDP(24) + (-t544 * t579 - t545 * t580 + t549 * t566 + t550 * t568) * MDP(12) + (t544 * t580 + t549 * t568) * MDP(11) + (t602 * t544 - t587 * t549 + (-t568 + t752) * t604) * MDP(17) + (-t549 * MDP(13) - t550 * MDP(14) + MDP(16) * t624 - MDP(17) * t634) * t606 + 0.2e1 * t685 * t758 + (t551 * t478 + t520 * t547 - t538 * t641 + t554 * t766) * MDP(23) + t757 * t760 + MDP(6) * t717 + (-t499 * t427 - t419 * t498 - t443 * t481 - (-t442 * t609 - t499 * t693) * t691) * MDP(35) + (t418 * t498 + t429 * t499 + t443 * t483 - t635 * t691) * MDP(34) + (t431 * t498 - t443 * t691) * MDP(36) + (t403 * t481 + t404 * t498 - t648 * t443 - t647 * t419 + (-t409 * t691 + t733 + (-t433 * t498 + t445 * t691 + t734) * qJD(6)) * t614 + t745 * t609) * MDP(37) + (t403 * t483 - t411 * t443 - t647 * t418 + ((-qJD(6) * t445 + t409) * t691 - t733 - (-qJD(6) * t433 + t405) * t498 - qJD(6) * t734) * t609 + t745 * t614) * MDP(38) + (MDP(27) * t442 - MDP(28) * t443 - MDP(30) * t403 - MDP(31) * t402) * t596 + (t602 * t545 + t587 * t550 + (qJD(1) * t579 + t566) * t604) * MDP(16); (-t501 * t765 + t596 * t710 + t623) * MDP(31) + (t496 * t501 - t596 * t709 + t622) * MDP(30) + t620 * t757 + (t552 * t641 + t605 * t780 + t779) * MDP(23) + (t552 * t640 + t605 * t781 + t767) * MDP(24) + t621 + (t523 * t419 + t650 * t609 + t709 * t481 - (t609 * t710 - t614 * t658) * t691 + t776) * MDP(37) + (t568 * t603 + t705 * t606 + (qJD(3) * t684 + t576) * t617 + t636) * MDP(17) + (-t566 * t603 + t724 - t666 * t606 + (t612 * t684 - t573) * qJD(3) + t668) * MDP(16) + (t523 * t418 + t650 * t614 + t709 * t483 - (t609 * t658 + t614 * t710) * t691 + t777) * MDP(38) - t716 * t758 + (MDP(9) * t613 * t620 + MDP(10) * t716) * pkin(1); (-t503 * t765 + t596 * t711 + t623) * MDP(31) + (t496 * t503 - t596 * t706 + t622) * MDP(30) + (t606 * t669 + t636 - t754) * MDP(17) + (-t606 * t639 + t626 + t724) * MDP(16) + t621 + (t559 * t419 + t649 * t609 + t706 * t481 - (t609 * t711 - t614 * t657) * t691 + t776) * MDP(37) + (t559 * t418 + t649 * t614 + t706 * t483 - (t609 * t657 + t614 * t711) * t691 + t777) * MDP(38) + (-t640 * t741 + t708 * t605 + (qJD(4) * t683 - t474) * t616 + t778) * MDP(24) + (-t641 * t741 + t726 - t672 * t605 + (t611 * t683 - t510) * qJD(4) + t676) * MDP(23); (-t605 * t646 + t779) * MDP(23) + (t605 * t673 + t767) * MDP(24) + (-t496 * t740 + t436 * t596 - t753 + (t610 * t682 - t731) * qJD(5) - t679) * MDP(30) + (t765 * t740 + t437 * t596 - t772 + (qJD(5) * t682 - t413) * t615 - t678) * MDP(31) + (t599 * t419 + t651 * t609 + t654 * t481 - (t609 * t653 - t614 * t656) * t691 + t776) * MDP(37) + (t599 * t418 + t651 * t614 + t654 * t483 - (t609 * t656 + t614 * t653) * t691 + t777) * MDP(38) + t788; (t435 * t596 + t622) * MDP(30) + (t434 * t596 + t623) * MDP(31) + (-pkin(5) * t419 + (-t434 * t609 + t451 * t614) * t691 - t435 * t481 - t609 * t787 - t713 * pkin(11) + t776) * MDP(37) + (-pkin(5) * t418 - (t434 * t614 + t451 * t609) * t691 - t435 * t483 - t432 * t786 + (-t691 * t694 - t429) * pkin(11) + t777) * MDP(38) + t774; t483 * t481 * MDP(32) + (-t481 ^ 2 + t483 ^ 2) * MDP(33) + (-t481 * t691 + t712) * MDP(34) + (-t483 * t691 - t736) * MDP(35) + t431 * MDP(36) + (-t400 * t609 - t411 * t691 - t432 * t483 + t404) * MDP(37) + (-t400 * t614 - t405 * t609 + t432 * t481 + t648 * t691) * MDP(38) + (-MDP(34) * t727 - MDP(35) * t483 - MDP(37) * t411 + MDP(38) * t648) * qJD(6);];
tauc  = t1;
