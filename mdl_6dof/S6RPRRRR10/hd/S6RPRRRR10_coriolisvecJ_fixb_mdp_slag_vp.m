% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:34:32
% EndTime: 2019-03-09 07:34:55
% DurationCPUTime: 13.94s
% Computational Cost: add. (20018->595), mult. (65244->831), div. (0->0), fcn. (57405->14), ass. (0->246)
t664 = sin(qJ(3));
t657 = sin(pkin(13));
t659 = sin(pkin(6));
t762 = qJD(1) * t659;
t745 = t657 * t762;
t660 = cos(pkin(13));
t793 = cos(pkin(7));
t795 = cos(qJ(3));
t713 = t793 * t795;
t691 = t660 * t713;
t658 = sin(pkin(7));
t794 = cos(pkin(6));
t740 = t794 * t658;
t693 = t795 * t740;
t833 = t659 * t691 + t693;
t820 = qJD(1) * t833 - t664 * t745;
t590 = qJD(4) - t820;
t588 = qJD(5) + t590;
t741 = t664 * t793;
t678 = t659 * (t657 * t741 - t660 * t795);
t742 = qJD(3) * t795;
t835 = qJD(1) * t678 + t658 * t742;
t665 = cos(qJ(6));
t755 = qJD(6) * t665;
t604 = t659 * (t657 * t795 + t660 * t741) + t664 * t740;
t597 = qJD(1) * t604;
t708 = t794 * t793;
t744 = t660 * t762;
t619 = -qJD(1) * t708 + t658 * t744 - qJD(3);
t663 = sin(qJ(4));
t667 = cos(qJ(4));
t556 = -t663 * t597 - t619 * t667;
t557 = t597 * t667 - t619 * t663;
t662 = sin(qJ(5));
t666 = cos(qJ(5));
t507 = -t666 * t556 + t557 * t662;
t824 = t507 * t665;
t834 = t755 + t824;
t832 = t657 * t713 + t660 * t664;
t549 = pkin(3) * t597 - pkin(10) * t820;
t548 = t667 * t549;
t796 = pkin(10) + pkin(11);
t747 = qJD(4) * t796;
t739 = qJD(1) * t794;
t720 = pkin(1) * t739;
t622 = qJ(2) * t744 + t657 * t720;
t777 = t659 * t660;
t677 = (t777 * t793 + t740) * pkin(9);
t580 = qJD(1) * t677 + t622;
t644 = t660 * t720;
t780 = t657 * t659;
t670 = t794 * pkin(2) + (-pkin(9) * t793 - qJ(2)) * t780;
t591 = qJD(1) * t670 + t644;
t781 = t657 * t658;
t618 = (-pkin(2) * t660 - pkin(9) * t781 - pkin(1)) * t659;
t611 = qJD(1) * t618 + qJD(2);
t746 = t658 * t795;
t807 = -t664 * t580 + t591 * t713 + t611 * t746;
t830 = pkin(4) * t597 - t663 * t807 + t548 + (-pkin(11) * t820 + t747) * t667;
t769 = t663 * t549 + t667 * t807;
t783 = t820 * t663;
t829 = -pkin(11) * t783 + t663 * t747 + t769;
t778 = t658 * t664;
t625 = t663 * t793 + t667 * t778;
t723 = t658 * t745;
t828 = qJD(4) * t625 + t835 * t663 + t667 * t723;
t624 = -t663 * t778 + t667 * t793;
t827 = -qJD(4) * t624 + t663 * t723 - t835 * t667;
t697 = t556 * t662 + t666 * t557;
t577 = t820 * qJD(3);
t759 = qJD(4) * t667;
t760 = qJD(4) * t663;
t530 = t667 * t577 - t597 * t760 - t619 * t759;
t531 = qJD(4) * t557 + t577 * t663;
t757 = qJD(5) * t666;
t758 = qJD(5) * t662;
t465 = t666 * t530 - t662 * t531 + t556 * t757 - t557 * t758;
t596 = t604 * qJD(3);
t578 = qJD(1) * t596;
t661 = sin(qJ(6));
t749 = t665 * t465 + t661 * t578 + t588 * t755;
t756 = qJD(6) * t661;
t447 = -t697 * t756 + t749;
t445 = t447 * t661;
t490 = t588 * t661 + t665 * t697;
t735 = t465 * t661 - t665 * t578;
t448 = qJD(6) * t490 + t735;
t466 = qJD(5) * t697 + t530 * t662 + t666 * t531;
t462 = t665 * t466;
t787 = t697 * t661;
t488 = -t665 * t588 + t787;
t753 = -qJD(6) - t507;
t460 = t661 * t466;
t773 = -t753 * t755 + t460;
t823 = t753 * t661;
t826 = t578 * MDP(26) - t466 * MDP(25) - t507 ^ 2 * MDP(23) + (t507 * t588 + t465) * MDP(24) + (MDP(22) * t507 + MDP(23) * t697 + MDP(25) * t588 + MDP(33) * t753) * t697 + (t834 * t490 + t445) * MDP(29) + (-t490 * t697 - t753 * t824 + t773) * MDP(31) + (t488 * t697 - t753 * t823 + t462) * MDP(32) + (t447 * t665 - t661 * t448 - t834 * t488 + t490 * t823) * MDP(30);
t550 = -t591 * t658 + t793 * t611;
t500 = -pkin(3) * t820 - pkin(10) * t597 + t550;
t524 = t795 * t580 + t591 * t741 + t611 * t778;
t502 = -pkin(10) * t619 + t524;
t476 = t667 * t500 - t502 * t663;
t463 = -pkin(11) * t557 + t476;
t458 = pkin(4) * t590 + t463;
t477 = t500 * t663 + t502 * t667;
t464 = pkin(11) * t556 + t477;
t791 = t464 * t662;
t431 = t458 * t666 - t791;
t429 = -pkin(5) * t588 - t431;
t825 = t429 * t507;
t630 = t662 * t663 - t666 * t667;
t767 = t588 * t630;
t631 = t662 * t667 + t663 * t666;
t766 = t588 * t631;
t761 = qJD(2) * t659;
t822 = t832 * t761;
t821 = MDP(4) * t657 + MDP(5) * t660;
t480 = pkin(5) * t697 + pkin(12) * t507;
t501 = t619 * pkin(3) - t807;
t485 = -t556 * pkin(4) + t501;
t671 = qJD(2) * t678;
t804 = qJD(3) * t807;
t496 = -qJD(1) * t671 + t804;
t721 = t761 * t781;
t539 = pkin(3) * t578 - pkin(10) * t577 + qJD(1) * t721;
t734 = -t496 * t663 + t667 * t539;
t673 = -qJD(4) * t477 + t734;
t437 = pkin(4) * t578 - pkin(11) * t530 + t673;
t687 = t667 * t496 + t500 * t759 - t502 * t760 + t663 * t539;
t441 = -pkin(11) * t531 + t687;
t724 = -t662 * t437 - t666 * t441 - t458 * t757 + t464 * t758;
t819 = t485 * t507 + t724;
t813 = (t657 ^ 2 + t660 ^ 2) * MDP(6) * t659 ^ 2;
t623 = t658 * t777 - t708;
t563 = t604 * t667 - t623 * t663;
t779 = t657 * t664;
t603 = t659 * t779 - t833;
t748 = pkin(1) * t794;
t648 = t660 * t748;
t605 = t648 + t670;
t558 = -t605 * t658 + t793 * t618;
t516 = pkin(3) * t603 - pkin(10) * t604 + t558;
t765 = qJ(2) * t777 + t657 * t748;
t600 = t677 + t765;
t522 = t795 * t600 - t623 * pkin(10) + (t605 * t793 + t618 * t658) * t664;
t732 = t667 * t516 - t522 * t663;
t469 = pkin(4) * t603 - pkin(11) * t563 + t732;
t562 = t604 * t663 + t623 * t667;
t770 = t663 * t516 + t667 * t522;
t475 = -pkin(11) * t562 + t770;
t812 = t662 * t469 + t666 * t475;
t716 = -t524 + (t760 - t783) * pkin(4);
t695 = t624 * t666 - t625 * t662;
t811 = -qJD(5) * t695 + t662 * t828 + t666 * t827;
t582 = t624 * t662 + t625 * t666;
t810 = qJD(5) * t582 - t662 * t827 + t666 * t828;
t641 = t796 * t663;
t642 = t796 * t667;
t694 = -t641 * t666 - t642 * t662;
t809 = -qJD(5) * t694 + t662 * t830 + t829 * t666;
t615 = -t641 * t662 + t642 * t666;
t808 = -qJD(5) * t615 + t829 * t662 - t666 * t830;
t743 = qJD(3) * t778;
t806 = t762 * t832 - t743;
t805 = -t664 * t600 + t605 * t713 + t618 * t746;
t790 = t464 * t666;
t432 = t458 * t662 + t790;
t736 = -t666 * t437 + t441 * t662;
t674 = -qJD(5) * t432 - t736;
t419 = -pkin(5) * t578 - t674;
t430 = pkin(12) * t588 + t432;
t455 = t507 * pkin(5) - pkin(12) * t697 + t485;
t705 = t430 * t661 - t455 * t665;
t802 = -t419 * t665 + t429 * t756 + t697 * t705;
t417 = t419 * t661;
t423 = t430 * t665 + t455 * t661;
t801 = t423 * t697 + t429 * t755 + t417;
t799 = -t485 * t697 + t674;
t786 = t556 * t590;
t785 = t557 * t590;
t782 = t631 * t665;
t768 = pkin(5) * t597 - t808;
t752 = qJD(1) * qJD(2);
t653 = -pkin(4) * t667 - pkin(3);
t418 = pkin(12) * t578 - t724;
t712 = qJD(3) * t741;
t497 = qJD(1) * t822 + t580 * t742 + t591 * t712 + t611 * t743;
t481 = pkin(4) * t531 + t497;
t425 = pkin(5) * t466 - pkin(12) * t465 + t481;
t738 = -t418 * t661 + t665 * t425;
t511 = qJD(3) * t805 - t671;
t595 = (t693 + (t691 - t779) * t659) * qJD(3);
t544 = pkin(3) * t596 - pkin(10) * t595 + t721;
t733 = -t511 * t663 + t667 * t544;
t731 = -t665 * t597 + t661 * t767;
t730 = t597 * t661 + t665 * t767;
t728 = t590 * t667;
t651 = pkin(4) * t662 + pkin(12);
t725 = pkin(4) * t557 + qJD(6) * t651 + t480;
t433 = t463 * t662 + t790;
t715 = pkin(4) * t758 - t433;
t434 = t463 * t666 - t791;
t714 = -pkin(4) * t757 + t434;
t601 = pkin(5) * t630 - pkin(12) * t631 + t653;
t710 = pkin(12) * t597 - qJD(6) * t601 + t809;
t709 = -pkin(5) * t766 - pkin(12) * t767 + qJD(6) * t615 - t716;
t512 = t600 * t742 + t605 * t712 + t618 * t743 + t822;
t707 = t418 * t665 + t425 * t661;
t706 = -t466 * t651 + t825;
t439 = pkin(12) * t603 + t812;
t521 = t623 * pkin(3) - t805;
t491 = t562 * pkin(4) + t521;
t525 = t666 * t562 + t563 * t662;
t526 = -t562 * t662 + t563 * t666;
t456 = t525 * pkin(5) - t526 * pkin(12) + t491;
t704 = t439 * t665 + t456 * t661;
t703 = -t439 * t661 + t456 * t665;
t538 = -qJD(4) * t562 + t595 * t667;
t443 = pkin(4) * t596 - pkin(11) * t538 - qJD(4) * t770 + t733;
t537 = qJD(4) * t563 + t595 * t663;
t686 = t667 * t511 + t516 * t759 - t522 * t760 + t663 * t544;
t450 = -pkin(11) * t537 + t686;
t702 = t443 * t666 - t450 * t662;
t701 = t469 * t666 - t475 * t662;
t494 = t526 * t665 + t603 * t661;
t493 = t526 * t661 - t603 * t665;
t696 = (-qJ(2) * t745 + t644) * t657 - t622 * t660;
t688 = -pkin(10) * t578 + t501 * t590;
t685 = t662 * t443 + t666 * t450 + t469 * t757 - t475 * t758;
t684 = -t661 * t582 - t665 * t746;
t683 = -t665 * t582 + t661 * t746;
t682 = t631 * t755 - t731;
t681 = -t631 * t756 - t730;
t484 = pkin(4) * t537 + t512;
t652 = -pkin(4) * t666 - pkin(5);
t551 = t578 * t603;
t473 = qJD(5) * t526 + t666 * t537 + t538 * t662;
t472 = -qJD(5) * t525 - t537 * t662 + t538 * t666;
t454 = qJD(6) * t494 + t472 * t661 - t596 * t665;
t453 = -qJD(6) * t493 + t472 * t665 + t596 * t661;
t438 = -pkin(5) * t603 - t701;
t426 = pkin(5) * t473 - pkin(12) * t472 + t484;
t421 = -pkin(5) * t596 + qJD(5) * t812 - t702;
t420 = pkin(12) * t596 + t685;
t416 = -qJD(6) * t423 + t738;
t415 = -qJD(6) * t705 + t707;
t1 = [(t733 * t590 + t732 * t578 + t734 * t603 + t476 * t596 - t512 * t556 + t521 * t531 + t497 * t562 + t501 * t537 + (-t477 * t603 - t590 * t770) * qJD(4)) * MDP(20) + (((t660 * t765 + (qJ(2) * t780 - t648) * t657) * qJD(1) - t696) * MDP(7) - 0.2e1 * t821 * t739) * t761 + (-t477 * t596 + t497 * t563 + t501 * t538 + t512 * t557 + t521 * t530 - t578 * t770 - t590 * t686 - t603 * t687) * MDP(21) + (-t465 * t525 - t466 * t526 - t472 * t507 - t473 * t697) * MDP(23) + (t465 * t526 + t472 * t697) * MDP(22) + (t465 * t603 + t472 * t588 + t526 * t578 + t596 * t697) * MDP(24) + 0.2e1 * t752 * t813 + (-t530 * t562 - t531 * t563 - t537 * t557 + t538 * t556) * MDP(16) + (t530 * t563 + t538 * t557) * MDP(15) + (-t577 * t623 - t595 * t619) * MDP(10) + (t578 * t623 + t596 * t619) * MDP(11) + (-t577 * t603 - t578 * t604 + t595 * t820 - t596 * t597) * MDP(9) + (t497 * t623 + t512 * t619 + t550 * t596 + t558 * t578 + (qJD(1) * t603 - t820) * t721) * MDP(13) + (t577 * t604 + t595 * t597) * MDP(8) + (-t466 * t603 - t473 * t588 - t507 * t596 - t525 * t578) * MDP(25) + (-t531 * t603 - t537 * t590 + t556 * t596 - t562 * t578) * MDP(18) + (t530 * t603 + t538 * t590 + t557 * t596 + t563 * t578) * MDP(17) + (t588 * t596 + t551) * MDP(26) + (t590 * t596 + t551) * MDP(19) + (t496 * t623 + t511 * t619 + t550 * t595 + t558 * t577 + 0.2e1 * t597 * t721) * MDP(14) + ((qJD(6) * t703 + t420 * t665 + t426 * t661) * t753 - t704 * t466 - t415 * t525 - t423 * t473 + t421 * t490 + t438 * t447 + t419 * t494 + t429 * t453) * MDP(35) + (t466 * t525 - t473 * t753) * MDP(33) + (-t448 * t525 + t454 * t753 - t466 * t493 - t473 * t488) * MDP(32) + (t447 * t525 - t453 * t753 + t466 * t494 + t473 * t490) * MDP(31) + (-(-qJD(6) * t704 - t420 * t661 + t426 * t665) * t753 + t703 * t466 + t416 * t525 - t705 * t473 + t421 * t488 + t438 * t448 + t419 * t493 + t429 * t454) * MDP(34) + (t702 * t588 + t701 * t578 - t736 * t603 + t431 * t596 + t484 * t507 + t491 * t466 + t481 * t525 + t485 * t473 + (-t432 * t603 - t588 * t812) * qJD(5)) * MDP(27) + (-t432 * t596 + t491 * t465 + t485 * t472 + t481 * t526 + t484 * t697 - t578 * t812 - t588 * t685 + t603 * t724) * MDP(28) + (t447 * t494 + t453 * t490) * MDP(29) + (-t447 * t493 - t448 * t494 - t453 * t488 - t454 * t490) * MDP(30); t696 * MDP(7) * t762 + (t578 * t793 - t619 * t806 + t723 * t820) * MDP(13) + (t793 * t577 - t597 * t723 + t835 * t619) * MDP(14) + (-t531 * t746 + t806 * t556 + t624 * t578 - t590 * t828) * MDP(20) + (-t530 * t746 - t806 * t557 - t625 * t578 + t590 * t827) * MDP(21) + (-t466 * t746 - t507 * t806 + t578 * t695 - t588 * t810) * MDP(27) + (-t465 * t746 - t582 * t578 + t588 * t811 - t697 * t806) * MDP(28) + (-t695 * t448 + t684 * t466 - (qJD(6) * t683 + t661 * t811 - t665 * t806) * t753 + t810 * t488) * MDP(34) + (-t695 * t447 + t683 * t466 - (-qJD(6) * t684 + t661 * t806 + t665 * t811) * t753 + t810 * t490) * MDP(35) + (t659 * t794 * t821 - t813) * qJD(1) ^ 2; (-pkin(3) * t530 + t497 * t663 - t524 * t557 + (pkin(10) * t760 + t769) * t590 + t688 * t667) * MDP(21) + (t530 * t663 + t557 * t728) * MDP(15) + (t465 * t631 - t697 * t767) * MDP(22) + (-t465 * t630 - t466 * t631 + t507 * t767 - t697 * t766) * MDP(23) + (t578 * t631 - t588 * t767) * MDP(24) + (t653 * t465 + t481 * t631 - t767 * t485 - t615 * t578 + t588 * t809 + t716 * t697) * MDP(28) + (t578 * t663 + t590 * t728) * MDP(17) + (-t578 * t630 - t588 * t766) * MDP(25) + (-t590 ^ 2 * t663 + t578 * t667) * MDP(18) + (-t524 * t619 - t497) * MDP(13) - MDP(11) * t578 + (-MDP(11) * t619 - MDP(13) * t550 - MDP(17) * t557 - MDP(18) * t556 - MDP(19) * t590 - MDP(20) * t476 + MDP(21) * t477 - MDP(24) * t697 + MDP(25) * t507 - MDP(26) * t588 - MDP(27) * t431 + MDP(28) * t432 - MDP(8) * t820 + MDP(9) * t597) * t597 + (t619 * t820 + t577) * MDP(10) + (-t550 * t820 - t619 * t807 + t678 * t752 - t804) * MDP(14) - t820 ^ 2 * MDP(9) + (t466 * t630 - t753 * t766) * MDP(33) + (t447 * t630 + t462 * t631 + t490 * t766 - t681 * t753) * MDP(31) + (-t448 * t630 - t460 * t631 - t488 * t766 + t682 * t753) * MDP(32) + (-pkin(3) * t531 - t497 * t667 + t524 * t556 + (-pkin(10) * t759 - t548) * t590 + (t590 * t807 + t688) * t663) * MDP(20) + (t653 * t466 + t481 * t630 + t766 * t485 + t716 * t507 + t578 * t694 + t588 * t808) * MDP(27) + (-(t601 * t661 + t615 * t665) * t466 - t415 * t630 - t694 * t447 + t419 * t782 - (t661 * t709 + t665 * t710) * t753 + t768 * t490 - t766 * t423 + t681 * t429) * MDP(35) + ((t601 * t665 - t615 * t661) * t466 + t416 * t630 - t694 * t448 + t631 * t417 - (t661 * t710 - t665 * t709) * t753 + t768 * t488 - t766 * t705 + t682 * t429) * MDP(34) + (t447 * t782 + t490 * t681) * MDP(29) + ((t530 + t786) * t667 + (-t531 - t785) * t663) * MDP(16) + (t731 * t490 + t730 * t488 + (-t445 - t448 * t665 + (t488 * t661 - t490 * t665) * qJD(6)) * t631) * MDP(30); (t434 * t588 + (-t557 * t697 - t578 * t662 - t588 * t757) * pkin(4) + t819) * MDP(28) + (t652 * t448 + t706 * t661 + t715 * t488 - (t661 * t714 - t665 * t725) * t753 + t802) * MDP(34) + (t652 * t447 + t706 * t665 + t715 * t490 - (t661 * t725 + t665 * t714) * t753 + t801) * MDP(35) + (t433 * t588 + (-t507 * t557 + t578 * t666 - t588 * t758) * pkin(4) + t799) * MDP(27) + (t477 * t590 - t501 * t557 + t673) * MDP(20) + (t476 * t590 - t501 * t556 - t687) * MDP(21) + t578 * MDP(19) + (-t556 ^ 2 + t557 ^ 2) * MDP(16) - t557 * t556 * MDP(15) + (-t531 + t785) * MDP(18) + (t530 - t786) * MDP(17) + t826; (t432 * t588 + t799) * MDP(27) + (t431 * t588 + t819) * MDP(28) + (-pkin(5) * t448 + (-t431 * t661 + t480 * t665) * t753 - t432 * t488 + t661 * t825 - t773 * pkin(12) + t802) * MDP(34) + (-pkin(5) * t447 - (t431 * t665 + t480 * t661) * t753 - t432 * t490 + t429 * t824 + (-t753 * t756 - t462) * pkin(12) + t801) * MDP(35) + t826; t490 * t488 * MDP(29) + (-t488 ^ 2 + t490 ^ 2) * MDP(30) + (-t488 * t753 + t749) * MDP(31) + (-t490 * t753 - t735) * MDP(32) + t466 * MDP(33) + (-t423 * t753 - t429 * t490 + t738) * MDP(34) + (t429 * t488 + t705 * t753 - t707) * MDP(35) + (-MDP(31) * t787 - MDP(32) * t490 - MDP(34) * t423 + MDP(35) * t705) * qJD(6);];
tauc  = t1;
