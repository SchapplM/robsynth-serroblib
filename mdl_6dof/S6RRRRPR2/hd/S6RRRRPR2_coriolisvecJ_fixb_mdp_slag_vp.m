% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:59:30
% EndTime: 2019-03-09 21:59:48
% DurationCPUTime: 10.39s
% Computational Cost: add. (13003->523), mult. (32955->692), div. (0->0), fcn. (25066->10), ass. (0->248)
t647 = sin(qJ(3));
t648 = sin(qJ(2));
t725 = qJD(1) * t648;
t708 = t647 * t725;
t651 = cos(qJ(3));
t652 = cos(qJ(2));
t724 = qJD(1) * t652;
t709 = t651 * t724;
t590 = -t708 + t709;
t591 = -t647 * t724 - t651 * t725;
t646 = sin(qJ(4));
t650 = cos(qJ(4));
t561 = t650 * t590 + t591 * t646;
t554 = qJD(6) - t561;
t644 = cos(pkin(11));
t649 = cos(qJ(6));
t740 = t649 * t644;
t643 = sin(pkin(11));
t645 = sin(qJ(6));
t743 = t643 * t645;
t604 = -t740 + t743;
t795 = t554 * t604;
t605 = t643 * t649 + t644 * t645;
t794 = t554 * t605;
t790 = t561 * t643;
t549 = pkin(5) * t790;
t712 = pkin(10) * t790;
t585 = t591 * pkin(9);
t768 = pkin(7) + pkin(8);
t618 = t768 * t652;
t610 = qJD(1) * t618;
t592 = t647 * t610;
t617 = t768 * t648;
t608 = qJD(1) * t617;
t763 = qJD(2) * pkin(2);
t600 = -t608 + t763;
t697 = t651 * t600 - t592;
t538 = t585 + t697;
t640 = qJD(2) + qJD(3);
t531 = pkin(3) * t640 + t538;
t596 = t651 * t610;
t680 = -t600 * t647 - t596;
t764 = pkin(9) * t590;
t539 = -t680 + t764;
t537 = t650 * t539;
t483 = t646 * t531 + t537;
t713 = qJD(3) + qJD(4);
t638 = qJD(2) + t713;
t478 = qJ(5) * t638 + t483;
t635 = -pkin(2) * t652 - pkin(1);
t616 = t635 * qJD(1);
t574 = -pkin(3) * t590 + t616;
t681 = t590 * t646 - t650 * t591;
t499 = -pkin(4) * t561 - qJ(5) * t681 + t574;
t449 = -t478 * t643 + t644 * t499;
t789 = t561 * t644;
t793 = t449 * t789;
t607 = t647 * t652 + t648 * t651;
t775 = qJD(1) * t607;
t656 = t640 * t775;
t710 = qJD(2) * t768;
t689 = qJD(1) * t710;
t602 = t652 * t689;
t723 = qJD(3) * t647;
t695 = -t647 * t602 - t610 * t723;
t601 = t648 * t689;
t777 = t651 * (qJD(3) * t600 - t601);
t487 = -pkin(9) * t656 + t695 + t777;
t714 = qJD(1) * qJD(2);
t707 = t652 * t714;
t566 = qJD(3) * t709 - t640 * t708 + t651 * t707;
t696 = t647 * t601 - t651 * t602;
t658 = qJD(3) * t680 + t696;
t488 = -pkin(9) * t566 + t658;
t721 = qJD(4) * t646;
t699 = t646 * t488 - t539 * t721;
t662 = -(qJD(4) * t531 + t487) * t650 - t699;
t747 = t574 * t561;
t792 = t662 - t747;
t688 = pkin(5) * t681 - pkin(10) * t789;
t543 = t638 * t643 + t644 * t681;
t429 = -pkin(5) * t561 - pkin(10) * t543 + t449;
t450 = t644 * t478 + t643 * t499;
t541 = -t644 * t638 + t643 * t681;
t435 = -pkin(10) * t541 + t450;
t409 = t429 * t645 + t435 * t649;
t720 = qJD(4) * t650;
t426 = t646 * t487 - t650 * t488 + t531 * t721 + t539 * t720;
t479 = t650 * t566 + t590 * t720 + t591 * t721 - t646 * t656;
t759 = t479 * t643;
t416 = pkin(5) * t759 + t426;
t535 = t646 * t539;
t482 = t531 * t650 - t535;
t477 = -pkin(4) * t638 + qJD(5) - t482;
t462 = pkin(5) * t541 + t477;
t665 = t409 * t681 + t416 * t605 - t795 * t462;
t408 = t429 * t649 - t435 * t645;
t666 = -t408 * t681 + t416 * t604 + t794 * t462;
t678 = -qJD(6) * t543 - t759;
t716 = qJD(6) * t649;
t732 = t479 * t740 - t541 * t716;
t421 = t678 * t645 + t732;
t785 = t541 * t645 - t543 * t649;
t422 = -qJD(6) * t785 + t479 * t605;
t480 = qJD(4) * t681 + t566 * t646 + t650 * t656;
t786 = t649 * t541;
t491 = t543 * t645 + t786;
t791 = -t561 ^ 2 * MDP(19) + (-t561 * t638 + t479) * MDP(20) + (-MDP(18) * t561 + MDP(19) * t681 + MDP(21) * t638) * t681 + t421 * t605 * MDP(29) + (-t421 * t604 - t605 * t422) * MDP(30) + (t795 * MDP(30) + t681 * MDP(32)) * t491 + (MDP(31) * t605 - MDP(32) * t604 - MDP(21)) * t480 + (t795 * MDP(29) + t794 * MDP(30) + t681 * MDP(31)) * t785 + (-t795 * MDP(31) - t794 * MDP(32) - MDP(33) * t681) * t554;
t522 = pkin(4) * t681 - qJ(5) * t561;
t782 = -0.2e1 * t714;
t779 = MDP(4) * t648;
t778 = MDP(5) * (t648 ^ 2 - t652 ^ 2);
t424 = qJD(5) * t638 - t662;
t636 = pkin(2) * t725;
t547 = pkin(3) * t656 + qJD(2) * t636;
t430 = t480 * pkin(4) - t479 * qJ(5) - qJD(5) * t681 + t547;
t406 = t644 * t424 + t643 * t430;
t404 = t406 * t644;
t405 = -t424 * t643 + t644 * t430;
t685 = -t405 * t643 + t404;
t495 = t538 * t650 - t535;
t767 = pkin(3) * t591;
t507 = t522 - t767;
t452 = -t495 * t643 + t644 * t507;
t623 = pkin(3) * t720 + qJD(5);
t703 = t623 * t643 + t452;
t453 = t644 * t495 + t643 * t507;
t702 = t623 * t644 - t453;
t694 = t608 * t647 - t596;
t545 = t694 - t764;
t729 = -t651 * t608 - t592;
t546 = t585 + t729;
t501 = t545 * t646 + t546 * t650;
t502 = t507 + t636;
t454 = -t501 * t643 + t644 * t502;
t628 = t646 * t647 * pkin(2);
t634 = pkin(2) * t651 + pkin(3);
t722 = qJD(3) * t651;
t770 = t650 * pkin(2) * t722 - t628 * t713 + t634 * t720;
t564 = qJD(5) + t770;
t701 = t564 * t643 + t454;
t455 = t644 * t501 + t643 * t502;
t700 = t564 * t644 - t455;
t457 = t644 * t482 + t643 * t522;
t776 = -qJD(5) * t644 + t457;
t493 = t538 * t646 + t537;
t686 = pkin(3) * t721 - t493;
t742 = t647 * t650;
t731 = -t545 * t650 + t546 * t646 - t634 * t721 - (t647 * t720 + (t646 * t651 + t742) * qJD(3)) * pkin(2);
t551 = -pkin(9) * t607 - t617 * t651 - t618 * t647;
t606 = t647 * t648 - t651 * t652;
t679 = t617 * t647 - t618 * t651;
t552 = -pkin(9) * t606 - t679;
t518 = -t551 * t650 + t646 * t552;
t704 = t426 * t643 + t450 * t681;
t669 = -t574 * t681 - t426;
t684 = -t426 * t644 - t449 * t681;
t766 = pkin(3) * t650;
t765 = pkin(5) * t644;
t639 = t644 * pkin(10);
t760 = t477 * t561;
t758 = t479 * t644;
t568 = t650 * t606 + t607 * t646;
t570 = t640 * t606;
t667 = t607 * qJD(3);
t571 = t607 * qJD(2) + t667;
t503 = -qJD(4) * t568 - t570 * t650 - t571 * t646;
t757 = t503 * t643;
t569 = -t606 * t646 + t607 * t650;
t749 = t569 * t643;
t748 = t569 * t644;
t746 = t616 * t591;
t653 = qJD(2) ^ 2;
t741 = t648 * t653;
t739 = t652 * t653;
t654 = qJD(1) ^ 2;
t738 = t652 * t654;
t504 = qJD(4) * t569 - t570 * t646 + t650 * t571;
t637 = t648 * t763;
t563 = pkin(3) * t571 + t637;
t440 = pkin(4) * t504 - qJ(5) * t503 - qJD(5) * t569 + t563;
t609 = t648 * t710;
t611 = t652 * t710;
t668 = -t651 * t609 - t647 * t611 - t617 * t722 - t618 * t723;
t510 = -pkin(9) * t571 + t668;
t657 = qJD(3) * t679 + t609 * t647 - t651 * t611;
t511 = pkin(9) * t570 + t657;
t443 = -qJD(4) * t518 + t510 * t650 + t511 * t646;
t411 = t643 * t440 + t644 * t443;
t578 = pkin(3) * t606 + t635;
t517 = pkin(4) * t568 - qJ(5) * t569 + t578;
t519 = t551 * t646 + t552 * t650;
t461 = t643 * t517 + t644 * t519;
t730 = -t549 - t731;
t718 = qJD(6) * t435;
t717 = qJD(6) * t569;
t631 = -pkin(4) - t765;
t706 = -pkin(2) * t640 - t600;
t705 = pkin(1) * t782;
t410 = t644 * t440 - t443 * t643;
t456 = -t482 * t643 + t644 * t522;
t460 = t644 * t517 - t519 * t643;
t401 = -pkin(10) * t759 + t406;
t691 = -qJD(6) * t429 - t401;
t584 = -t634 * t650 - pkin(4) + t628;
t687 = -t549 + t686;
t447 = pkin(5) * t568 - pkin(10) * t748 + t460;
t451 = -pkin(10) * t749 + t461;
t683 = t447 * t649 - t451 * t645;
t682 = t447 * t645 + t451 * t649;
t677 = t450 * t790 + t685 + t793;
t676 = -t616 * t590 - t695;
t630 = pkin(3) * t646 + qJ(5);
t599 = t630 * t644 + t639;
t675 = qJD(6) * t599 + t688 + t703;
t583 = pkin(2) * t742 + t634 * t646 + qJ(5);
t573 = t583 * t644 + t639;
t674 = qJD(6) * t573 + t688 + t701;
t598 = (-pkin(10) - t630) * t643;
t673 = -qJD(6) * t598 - t702 - t712;
t572 = (-pkin(10) - t583) * t643;
t672 = -qJD(6) * t572 - t700 - t712;
t615 = qJ(5) * t644 + t639;
t671 = qJD(5) * t643 + qJD(6) * t615 + t456 + t688;
t614 = (-pkin(10) - qJ(5)) * t643;
t670 = -qJD(6) * t614 - t712 + t776;
t664 = t426 * t569 + t477 * t503 + t479 * t518;
t663 = -pkin(4) * t479 - qJ(5) * t480 - (-qJD(5) + t477) * t561;
t661 = t479 * t584 - t480 * t583 + t561 * t564 - t760;
t633 = -pkin(4) - t766;
t660 = t479 * t633 - t480 * t630 + t561 * t623 - t760;
t444 = qJD(4) * t519 + t510 * t646 - t511 * t650;
t655 = t591 * t590 * MDP(11) + (-t590 * t640 + t566) * MDP(13) + (-t591 * t640 - t656) * MDP(14) + (-t590 ^ 2 + t591 ^ 2) * MDP(12) + t791;
t613 = t631 - t766;
t576 = t584 - t765;
t575 = t636 - t767;
t524 = t604 * t569;
t523 = t605 * t569;
t476 = pkin(5) * t749 + t518;
t465 = t483 + t549;
t432 = t503 * t605 + t716 * t748 - t717 * t743;
t431 = -t503 * t604 - t605 * t717;
t417 = pkin(5) * t757 + t444;
t407 = -pkin(10) * t757 + t411;
t402 = pkin(5) * t504 - t503 * t639 + t410;
t400 = pkin(5) * t480 - pkin(10) * t758 + t405;
t399 = t649 * t400;
t1 = [(MDP(20) * t503 - MDP(21) * t504 - MDP(23) * t444 - MDP(24) * t443) * t638 + (t421 * t568 + t431 * t554 - t480 * t524 - t504 * t785) * MDP(31) + (-t421 * t523 + t422 * t524 - t431 * t491 + t432 * t785) * MDP(30) + (-t421 * t524 - t431 * t785) * MDP(29) + (-(t402 * t645 + t407 * t649) * t554 - t682 * t480 - (t400 * t645 + t401 * t649) * t568 - t409 * t504 - t417 * t785 + t476 * t421 - t416 * t524 + t462 * t431 + (-t408 * t568 - t554 * t683) * qJD(6)) * MDP(35) + (-t590 * t637 + t616 * t571 + (t635 * t667 + (t648 * pkin(2) * t606 + t607 * t635) * qJD(2)) * qJD(1)) * MDP(16) + (t479 * t578 + t503 * t574 + t547 * t569 + t563 * t681) * MDP(24) + (t479 * t569 + t503 * t681) * MDP(18) + (-t410 * t543 - t411 * t541 + (-t405 * t569 - t449 * t503 - t460 * t479) * t644 + (-t406 * t569 - t450 * t503 - t461 * t479) * t643) * MDP(27) + 0.2e1 * t707 * t779 + (-t570 * MDP(13) - t571 * MDP(14) + MDP(16) * t657 - MDP(17) * t668) * t640 + (t635 * t566 - t616 * t570 + (-t591 + t775) * t637) * MDP(17) + (t566 * t607 + t570 * t591) * MDP(11) + (-t566 * t606 - t570 * t590 + t591 * t571 - t607 * t656) * MDP(12) + (-t479 * t568 - t480 * t569 + t503 * t561 - t504 * t681) * MDP(19) + (t480 * t578 + t504 * t574 + t547 * t568 - t561 * t563) * MDP(23) + (-t406 * t568 + t411 * t561 + t444 * t543 - t450 * t504 - t461 * t480 + t644 * t664) * MDP(26) + (t405 * t568 - t410 * t561 + t444 * t541 + t449 * t504 + t460 * t480 + t643 * t664) * MDP(25) + (-t422 * t568 - t432 * t554 - t480 * t523 - t491 * t504) * MDP(32) + (t480 * t568 + t504 * t554) * MDP(33) + (t405 * t460 + t406 * t461 + t410 * t449 + t411 * t450 + t426 * t518 + t444 * t477) * MDP(28) + MDP(6) * t739 + t778 * t782 + ((t402 * t649 - t407 * t645) * t554 + t683 * t480 + (-t401 * t645 + t399) * t568 + t408 * t504 + t417 * t491 + t476 * t422 + t416 * t523 + t462 * t432 + (-t409 * t568 - t554 * t682) * qJD(6)) * MDP(34) + (-pkin(7) * t739 + t648 * t705) * MDP(9) - MDP(7) * t741 + (pkin(7) * t741 + t652 * t705) * MDP(10); t655 + (t591 * t636 + t729 * t640 + (qJD(3) * t706 + t601) * t651 + t676) * MDP(17) + (t590 * t636 + t746 - t694 * t640 + (t647 * t706 - t596) * qJD(3) + t696) * MDP(16) + (t454 * t561 - t541 * t731 + t643 * t661 + t684) * MDP(25) + (-t455 * t561 - t543 * t731 + t644 * t661 + t704) * MDP(26) + (t426 * t584 - t449 * t701 + t450 * t700 - t477 * t731 + t583 * t685) * MDP(28) + (t561 * t575 + t638 * t731 + t669) * MDP(23) + (-t541 * t700 + t543 * t701 + t677) * MDP(27) - t738 * t779 + (-t575 * t681 + (t501 - t770) * t638 + t792) * MDP(24) + t654 * t778 + ((t572 * t649 - t573 * t645) * t480 + t576 * t422 + (t645 * t672 - t649 * t674) * t554 + t730 * t491 + t666) * MDP(34) + (-(t572 * t645 + t573 * t649) * t480 + t576 * t421 + (t645 * t674 + t649 * t672) * t554 - t730 * t785 + t665) * MDP(35) + (MDP(9) * t648 * t654 + MDP(10) * t738) * pkin(1); t655 + (t452 * t561 + t541 * t686 + t643 * t660 + t684) * MDP(25) + (-t453 * t561 + t543 * t686 + t644 * t660 + t704) * MDP(26) + (t426 * t633 - t449 * t703 + t450 * t702 + t477 * t686 + t630 * t685) * MDP(28) + (t681 * t767 + t495 * t638 - t747 + (-t487 + (-pkin(3) * t638 - t531) * qJD(4)) * t650 - t699) * MDP(24) + (t493 * t638 + (-t561 * t591 - t638 * t721) * pkin(3) + t669) * MDP(23) + (-t541 * t702 + t543 * t703 + t677) * MDP(27) + (t640 * t697 + t676 - t777) * MDP(17) + (-t640 * t680 + t658 + t746) * MDP(16) + ((t598 * t649 - t599 * t645) * t480 + t613 * t422 + (t645 * t673 - t649 * t675) * t554 + t687 * t491 + t666) * MDP(34) + (-(t598 * t645 + t599 * t649) * t480 + t613 * t421 + (t645 * t675 + t649 * t673) * t554 - t687 * t785 + t665) * MDP(35); (t483 * t638 + t669) * MDP(23) + (t482 * t638 + t792) * MDP(24) + (t456 * t561 - t483 * t541 + t643 * t663 + t684) * MDP(25) + (-t457 * t561 - t483 * t543 + t644 * t663 + t704) * MDP(26) + (t793 + t456 * t543 + t404 + t776 * t541 + (qJD(5) * t543 + t450 * t561 - t405) * t643) * MDP(27) + (-pkin(4) * t426 - t449 * t456 - t450 * t457 - t477 * t483 + (-t449 * t643 + t450 * t644) * qJD(5) + t685 * qJ(5)) * MDP(28) + ((t614 * t649 - t615 * t645) * t480 + t631 * t422 - t465 * t491 + (t645 * t670 - t649 * t671) * t554 + t666) * MDP(34) + (-(t614 * t645 + t615 * t649) * t480 + t631 * t421 + t465 * t785 + (t645 * t671 + t649 * t670) * t554 + t665) * MDP(35) + t791; (-t543 * t561 + t759) * MDP(25) + (t541 * t561 + t758) * MDP(26) + (-t541 ^ 2 - t543 ^ 2) * MDP(27) + (t449 * t543 + t450 * t541 + t426) * MDP(28) + (-t785 * t554 + t422) * MDP(34) + (-t554 * t786 + (-t543 * t554 + t678) * t645 + t732) * MDP(35); -t491 ^ 2 * MDP(30) + (t491 * t554 + t732) * MDP(31) + t480 * MDP(33) + (t409 * t554 + t399) * MDP(34) + (t408 * t554 + t462 * t491) * MDP(35) - (MDP(29) * t491 - MDP(30) * t785 + MDP(32) * t554 - MDP(34) * t462) * t785 + (MDP(32) * t678 - MDP(34) * t718 + MDP(35) * t691) * t649 + (t678 * MDP(31) + (qJD(6) * t541 - t758) * MDP(32) + t691 * MDP(34) + (-t400 + t718) * MDP(35)) * t645;];
tauc  = t1;
