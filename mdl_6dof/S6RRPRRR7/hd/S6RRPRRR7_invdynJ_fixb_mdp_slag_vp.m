% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:00:05
% EndTime: 2019-03-09 14:00:20
% DurationCPUTime: 11.39s
% Computational Cost: add. (6734->629), mult. (14060->798), div. (0->0), fcn. (10197->12), ass. (0->283)
t654 = sin(qJ(4));
t659 = cos(qJ(4));
t660 = cos(qJ(2));
t762 = qJD(1) * t660;
t655 = sin(qJ(2));
t763 = qJD(1) * t655;
t572 = -t654 * t762 + t659 * t763;
t643 = qJD(2) - qJD(4);
t653 = sin(qJ(5));
t658 = cos(qJ(5));
t536 = t572 * t658 - t643 * t653;
t748 = qJD(1) * qJD(2);
t733 = t655 * t748;
t746 = qJDD(1) * t660;
t840 = t733 - t746;
t732 = t660 * t748;
t747 = qJDD(1) * t655;
t841 = t732 + t747;
t677 = t654 * t840 + t659 * t841;
t577 = t654 * t655 + t659 * t660;
t681 = t577 * qJD(4);
t486 = -qJD(1) * t681 + t677;
t642 = qJDD(2) - qJDD(4);
t753 = qJD(5) * t658;
t754 = qJD(5) * t653;
t461 = t658 * t486 - t572 * t754 - t653 * t642 - t643 * t753;
t726 = t486 * t653 + t658 * t642;
t462 = qJD(5) * t536 + t726;
t534 = t572 * t653 + t658 * t643;
t652 = sin(qJ(6));
t657 = cos(qJ(6));
t751 = qJD(6) * t657;
t740 = t657 * t461 - t652 * t462 - t534 * t751;
t752 = qJD(6) * t652;
t434 = -t536 * t752 + t740;
t698 = t534 * t652 - t657 * t536;
t727 = t461 * t652 + t657 * t462;
t435 = -qJD(6) * t698 + t727;
t794 = t536 * t652;
t470 = t657 * t534 + t794;
t756 = qJD(4) * t659;
t757 = qJD(4) * t654;
t759 = qJD(2) * t660;
t855 = t654 * t759 + t655 * t756 - t660 * t757;
t487 = qJD(1) * t855 + qJDD(1) * t577 - t659 * t733;
t823 = t577 * qJD(1);
t839 = qJD(5) + t823;
t559 = qJD(6) + t839;
t788 = t652 * t653;
t576 = -t657 * t658 + t788;
t787 = t652 * t658;
t579 = t653 * t657 + t787;
t745 = qJD(5) + qJD(6);
t774 = (t745 + t823) * t579;
t825 = -qJD(6) * t658 - t753;
t675 = t825 * t657;
t776 = t576 * t823 + t745 * t788 + t675;
t484 = qJDD(5) + t487;
t786 = t653 * t484;
t797 = t484 * t658;
t483 = qJDD(6) + t484;
t798 = t483 * t579;
t799 = t483 * t576;
t800 = t461 * t653;
t829 = t839 * t536;
t830 = t839 ^ 2;
t856 = t534 * t839;
t857 = -(t653 * (t462 + t829) + (-t461 + t856) * t658) * MDP(23) + (t434 * t579 + t698 * t776) * MDP(29) + (-t434 * t576 - t435 * t579 + t470 * t776 + t698 * t774) * MDP(30) - (t559 * t776 - t572 * t698 - t798) * MDP(31) - (-t470 * t572 + t559 * t774 + t799) * MDP(32) + (t658 * t829 + t800) * MDP(22) + (-t536 * t572 + t658 * t830 + t786) * MDP(24) - (-t534 * t572 + t653 * t830 - t797) * MDP(25) - (t572 * t643 + t487) * MDP(18) - (-t572 ^ 2 + t823 ^ 2) * MDP(16) - t642 * MDP(19) + (MDP(15) * t823 - MDP(26) * t839 - MDP(33) * t559) * t572;
t630 = pkin(7) * t763;
t583 = pkin(8) * t763 - t630;
t854 = qJD(3) - t583;
t574 = -qJD(1) * pkin(1) - pkin(2) * t762 - qJ(3) * t763;
t544 = pkin(3) * t762 - t574;
t481 = pkin(4) * t823 - pkin(9) * t572 + t544;
t662 = -pkin(2) - pkin(3);
t738 = t662 * qJD(2);
t551 = t738 + t854;
t631 = pkin(7) * t762;
t585 = -pkin(8) * t762 + t631;
t646 = qJD(2) * qJ(3);
t573 = t585 + t646;
t505 = t654 * t551 + t659 * t573;
t492 = -pkin(9) * t643 + t505;
t454 = t481 * t653 + t492 * t658;
t445 = -pkin(10) * t534 + t454;
t453 = t658 * t481 - t492 * t653;
t444 = -pkin(10) * t536 + t453;
t439 = pkin(5) * t839 + t444;
t804 = t439 * t657;
t431 = -t445 * t652 + t804;
t613 = pkin(7) * t732;
t625 = pkin(7) * t747;
t731 = qJDD(3) + t613 + t625;
t516 = -pkin(8) * t841 + t662 * qJDD(2) + t731;
t626 = pkin(7) * t746;
t644 = qJDD(2) * qJ(3);
t645 = qJD(2) * qJD(3);
t545 = -pkin(7) * t733 + t626 + t644 + t645;
t517 = pkin(8) * t840 + t545;
t691 = -t659 * t516 + t654 * t517 + t551 * t757 + t573 * t756;
t449 = pkin(4) * t642 + t691;
t436 = pkin(5) * t462 + t449;
t504 = t551 * t659 - t654 * t573;
t491 = pkin(4) * t643 - t504;
t466 = pkin(5) * t534 + t491;
t649 = qJ(5) + qJ(6);
t636 = cos(t649);
t656 = sin(qJ(1));
t781 = t656 * t660;
t784 = t655 * t659;
t561 = t654 * t781 - t656 * t784;
t661 = cos(qJ(1));
t783 = t655 * t661;
t785 = t654 * t660;
t563 = -t659 * t783 + t661 * t785;
t678 = g(1) * t563 + g(2) * t561 + g(3) * t577;
t853 = t431 * t572 - t436 * t576 - t774 * t466 - t678 * t636;
t803 = t445 * t657;
t432 = t439 * t652 + t803;
t635 = sin(t649);
t852 = -t432 * t572 - t436 * t579 + t776 * t466 + t678 * t635;
t845 = t470 * t559;
t844 = t559 * t698;
t523 = t583 * t659 + t585 * t654;
t725 = -t654 * qJ(3) + t659 * t662;
t554 = qJD(3) * t659 + qJD(4) * t725;
t771 = t554 - t523;
t443 = t445 * t752;
t562 = t577 * t656;
t519 = -t562 * t636 - t635 * t661;
t564 = t577 * t661;
t521 = t564 * t636 - t656 * t635;
t578 = t784 - t785;
t809 = g(3) * t578;
t838 = g(1) * t521 - g(2) * t519 + t466 * t470 + t636 * t809 + t443;
t520 = -t564 * t635 - t636 * t656;
t697 = t562 * t635 - t636 * t661;
t634 = t655 * qJD(3);
t650 = qJDD(1) * pkin(1);
t706 = pkin(2) * t746 + qJ(3) * t841 + qJD(1) * t634 + t650;
t716 = t655 * t738;
t489 = pkin(3) * t746 + qJD(1) * t716 + t706;
t440 = pkin(4) * t487 - pkin(9) * t486 + t489;
t438 = t658 * t440;
t684 = t654 * t516 + t659 * t517 + t551 * t756 - t573 * t757;
t448 = -pkin(9) * t642 + t684;
t428 = pkin(5) * t484 - pkin(10) * t461 - qJD(5) * t454 - t653 * t448 + t438;
t683 = -t653 * t440 - t658 * t448 - t481 * t753 + t492 * t754;
t429 = -pkin(10) * t462 - t683;
t728 = t657 * t428 - t652 * t429;
t837 = -g(1) * t520 + g(2) * t697 + t466 * t698 + t635 * t809 + t728;
t832 = t483 * MDP(33) + (-t470 ^ 2 + t698 ^ 2) * MDP(30) - t470 * MDP(29) * t698;
t766 = t660 * pkin(2) + t655 * qJ(3);
t590 = -pkin(1) - t766;
t831 = t491 * t839;
t508 = t579 * t578;
t767 = t659 * qJ(3) + t654 * t662;
t770 = qJD(4) * t767 + t585 * t659 + t654 * t854;
t827 = -t652 * t754 - t653 * t752;
t792 = t823 * t653;
t826 = (t754 + t792) * pkin(5);
t810 = g(2) * t661;
t813 = g(1) * t656;
t824 = -t810 + t813;
t811 = g(2) * t656;
t812 = g(1) * t661;
t713 = t811 + t812;
t807 = pkin(7) * qJDD(2);
t818 = (qJD(1) * t590 + t574) * qJD(2) - t807;
t817 = pkin(7) - pkin(8);
t816 = pkin(9) + pkin(10);
t815 = pkin(5) * t658;
t814 = g(1) * t562;
t582 = -pkin(9) + t767;
t808 = pkin(10) - t582;
t806 = qJDD(2) * pkin(2);
t802 = t453 * t572;
t801 = t454 * t572;
t528 = qJD(2) * t577 - t681;
t796 = t528 * t653;
t795 = t528 * t658;
t793 = t823 * t643;
t790 = t578 * t653;
t789 = t578 * t658;
t664 = qJD(1) ^ 2;
t782 = t655 * t664;
t593 = t817 * t655;
t595 = t817 * t660;
t538 = t593 * t654 + t595 * t659;
t531 = t658 * t538;
t514 = pkin(4) * t572 + pkin(9) * t823;
t779 = t658 * t504 + t653 * t514;
t621 = qJ(3) * t762;
t560 = t662 * t763 + t621;
t485 = -t514 + t560;
t778 = t653 * t485 + t658 * t523;
t575 = t660 * pkin(3) - t590;
t502 = pkin(4) * t577 - pkin(9) * t578 + t575;
t773 = t653 * t502 + t531;
t772 = -t826 + t770;
t768 = qJ(3) * t759 + t634;
t647 = t655 ^ 2;
t648 = t660 ^ 2;
t765 = t647 - t648;
t761 = qJD(2) * t655;
t760 = qJD(2) * t659;
t758 = qJD(4) * t559;
t755 = qJD(5) * t839;
t743 = pkin(10) * t792;
t741 = t660 * t782;
t739 = qJD(5) * t816;
t730 = qJD(5) * t808;
t729 = -qJD(2) * pkin(2) + qJD(3);
t719 = -qJD(5) * t481 - t448;
t718 = qJD(6) * t439 + t429;
t717 = t643 ^ 2;
t581 = pkin(4) - t725;
t715 = -t505 + t826;
t663 = qJD(2) ^ 2;
t714 = pkin(7) * t663 + t810;
t712 = g(2) * t562 + t809;
t480 = t658 * t485;
t547 = t808 * t658;
t693 = -pkin(10) * t658 * t823 - pkin(5) * t572;
t711 = -qJD(6) * t547 + t653 * t771 - t658 * t730 + t480 + t693;
t507 = t658 * t514;
t594 = t816 * t658;
t710 = qJD(6) * t594 - t504 * t653 + t658 * t739 + t507 - t693;
t546 = t808 * t653;
t709 = -qJD(6) * t546 - t554 * t658 - t653 * t730 - t743 + t778;
t592 = t816 * t653;
t708 = qJD(6) * t592 + t653 * t739 + t743 + t779;
t707 = pkin(2) * t655 - qJ(3) * t660;
t589 = t630 + t729;
t591 = t631 + t646;
t696 = t589 * t660 - t591 * t655;
t695 = t593 * t659 - t595 * t654;
t694 = g(1) * t783 - g(3) * t660 + t655 * t811 - t625;
t690 = -0.2e1 * pkin(1) * t748 - t807;
t689 = t578 * t753 + t796;
t688 = -t578 * t754 + t795;
t687 = -pkin(9) * t484 + t831;
t685 = -qJDD(3) + t694;
t542 = t716 + t768;
t527 = -t655 * t760 + t855;
t465 = pkin(4) * t527 - pkin(9) * t528 + t542;
t584 = t817 * t761;
t586 = qJD(2) * t595;
t475 = qJD(4) * t695 - t584 * t659 + t586 * t654;
t682 = t653 * t465 + t658 * t475 + t502 * t753 - t538 * t754;
t679 = -t582 * t484 - t831;
t676 = -t714 + 0.2e1 * t650;
t673 = -t449 + t678;
t513 = pkin(2) * t733 - t706;
t565 = pkin(2) * t761 - t768;
t671 = -qJD(1) * t565 - qJDD(1) * t590 - t513 - t714;
t670 = pkin(9) * t755 - t673;
t669 = t582 * t755 + t673;
t476 = qJD(4) * t538 - t584 * t654 - t586 * t659;
t667 = t544 * t572 - t678 + t691;
t666 = -g(1) * t564 - t544 * t823 + t684 - t712;
t558 = t731 - t806;
t665 = qJD(2) * t696 + t545 * t660 + t558 * t655 - t713;
t624 = -pkin(4) - t815;
t615 = g(1) * t781;
t580 = pkin(2) * t763 - t621;
t571 = t653 * t763 + t658 * t760;
t568 = -t653 * t760 + t658 * t763;
t567 = t581 + t815;
t533 = t564 * t658 - t653 * t656;
t532 = -t564 * t653 - t656 * t658;
t509 = t576 * t578;
t503 = pkin(5) * t790 - t695;
t494 = t658 * t502;
t464 = t658 * t465;
t460 = -pkin(10) * t790 + t773;
t455 = pkin(5) * t577 - pkin(10) * t789 - t538 * t653 + t494;
t451 = pkin(5) * t689 + t476;
t442 = t528 * t787 + (t745 * t789 + t796) * t657 + t827 * t578;
t441 = -t508 * t745 - t576 * t528;
t433 = -pkin(10) * t689 + t682;
t430 = -pkin(10) * t795 + pkin(5) * t527 - t475 * t653 + t464 + (-t531 + (pkin(10) * t578 - t502) * t653) * qJD(5);
t1 = [(t461 * t789 + t536 * t688) * MDP(22) + (qJDD(1) * t647 + 0.2e1 * t655 * t732) * MDP(4) + 0.2e1 * (t655 * t746 - t748 * t765) * MDP(5) + ((t430 * t657 - t433 * t652) * t559 + (t455 * t657 - t460 * t652) * t483 + t728 * t577 + t431 * t527 + t451 * t470 + t503 * t435 + t436 * t508 + t466 * t442 - g(1) * t519 - g(2) * t521 + ((-t455 * t652 - t460 * t657) * t559 - t432 * t577) * qJD(6)) * MDP(34) + t713 * MDP(3) + (t655 * t690 + t660 * t676 + t615) * MDP(9) + ((t647 + t648) * qJDD(1) * pkin(7) + t665) * MDP(12) + (t690 * t660 + (-t676 - t813) * t655) * MDP(10) + (-g(1) * t561 + g(2) * t563 + t475 * t643 + t486 * t575 + t489 * t578 + t528 * t544 + t538 * t642 + t542 * t572) * MDP(21) + (t527 * t643 + t577 * t642) * MDP(18) + (-t528 * t643 - t578 * t642) * MDP(17) + (-t435 * t577 - t442 * t559 - t470 * t527 - t483 * t508) * MDP(32) + (t483 * t577 + t527 * t559) * MDP(33) + (t486 * t578 + t528 * t572) * MDP(15) + ((-t534 * t658 - t536 * t653) * t528 + (-t800 - t462 * t658 + (t534 * t653 - t536 * t658) * qJD(5)) * t578) * MDP(23) + (t461 * t577 + t484 * t789 + t527 * t536 + t688 * t839) * MDP(24) + (-t462 * t577 - t527 * t534 - t578 * t786 - t689 * t839) * MDP(25) + (t484 * t577 + t527 * t839) * MDP(26) + (-t682 * t839 - t773 * t484 + t683 * t577 - t454 * t527 + t476 * t536 - t695 * t461 + t449 * t789 - g(1) * (t562 * t653 - t658 * t661) - g(2) * t532 + t688 * t491) * MDP(28) + (-g(2) * t533 + t438 * t577 + t453 * t527 - t695 * t462 + t464 * t839 + t476 * t534 + t494 * t484 + (t814 + (t491 * t578 - t492 * t577 - t538 * t839) * qJD(5)) * t658 + ((-qJD(5) * t502 - t475) * t839 - t538 * t484 + t719 * t577 + t449 * t578 + t491 * t528 + t812) * t653) * MDP(27) + qJDD(1) * MDP(1) + (t443 * t577 - t432 * t527 - t451 * t698 + t503 * t434 - t436 * t509 + t466 * t441 - g(1) * t697 - g(2) * t520 + (-(-qJD(6) * t460 + t430) * t559 - t455 * t483 - t428 * t577) * t652 + (-(qJD(6) * t455 + t433) * t559 - t460 * t483 - t718 * t577) * t657) * MDP(35) + (t434 * t577 + t441 * t559 - t483 * t509 - t527 * t698) * MDP(31) + (-t434 * t508 + t435 * t509 - t441 * t470 + t442 * t698) * MDP(30) + (-t434 * t509 - t441 * t698) * MDP(29) + (-t486 * t577 - t487 * t578 - t527 * t572 - t528 * t823) * MDP(16) + (-g(2) * t564 + t476 * t643 + t487 * t575 + t489 * t577 + t527 * t544 + t542 * t823 - t642 * t695 + t814) * MDP(20) + (t665 * pkin(7) + t574 * t565 + (t513 - t824) * t590) * MDP(14) + t824 * MDP(2) + (qJDD(2) * t655 + t660 * t663) * MDP(6) + (qJDD(2) * t660 - t655 * t663) * MDP(7) + (t655 * t818 + t671 * t660 + t615) * MDP(11) + (-t818 * t660 + (t671 + t813) * t655) * MDP(13); (pkin(1) * t782 + t694) * MDP(9) + (-t560 * t572 + t642 * t767 + t643 * t771 + t666) * MDP(21) + (-t707 * qJDD(1) + ((t591 - t646) * t655 + (-t589 + t729) * t660) * qJD(1)) * MDP(12) + (t626 + 0.2e1 * t644 + 0.2e1 * t645 + (qJD(1) * t580 - g(3)) * t655 + (qJD(1) * t574 - t713) * t660) * MDP(13) + (g(3) * t655 - t626 + (pkin(1) * t664 + t713) * t660) * MDP(10) - t857 + (-(t546 * t652 - t547 * t657) * t483 + t567 * t434 + (t652 * t711 + t657 * t709) * t559 - t772 * t698 + t852) * MDP(35) + ((t546 * t657 + t547 * t652) * t483 + t567 * t435 + (t652 * t709 - t657 * t711) * t559 + t772 * t470 + t853) * MDP(34) - MDP(4) * t741 + (0.2e1 * t806 + (-t574 * t655 + t580 * t660) * qJD(1) + t685) * MDP(11) + (t581 * t461 + t778 * t839 - t801 + t770 * t536 + (-t554 * t839 + t679) * t658 + t669 * t653) * MDP(28) + (t802 + t581 * t462 - t480 * t839 + t770 * t534 + (-t771 * t839 + t679) * t653 - t669 * t658) * MDP(27) + t765 * MDP(5) * t664 + qJDD(2) * MDP(8) + (-t560 * t823 - t642 * t725 + t643 * t770 + t667) * MDP(20) + (qJD(4) * t823 - t677 + t793) * MDP(17) + (-t696 * qJD(1) * pkin(7) - t558 * pkin(2) - g(3) * t766 + t545 * qJ(3) + t591 * qJD(3) - t574 * t580 + t707 * t713) * MDP(14) + MDP(7) * t746 + MDP(6) * t747; (-qJDD(2) - t741) * MDP(11) + MDP(12) * t747 + (-t647 * t664 - t663) * MDP(13) + (-qJD(2) * t591 + t574 * t763 + t613 - t685 - t806) * MDP(14) + (-t642 * t659 - t654 * t717 - t763 * t823) * MDP(20) + (-t572 * t763 + t642 * t654 - t659 * t717) * MDP(21) + (-t462 * t659 + (-t653 * t756 - t568) * t839 + (-t534 * t643 - t753 * t839 - t786) * t654) * MDP(27) + (-t461 * t659 + (-t658 * t756 + t571) * t839 + (-t536 * t643 + t754 * t839 - t797) * t654) * MDP(28) + (-(t568 * t657 - t571 * t652) * t559 + (-t579 * t758 - t435) * t659 + ((t675 - t827) * t559 - t798 - t643 * t470) * t654) * MDP(34) + ((t568 * t652 + t571 * t657) * t559 + (t576 * t758 - t434) * t659 + (-(t652 * t825 - t653 * t751 - t657 * t754) * t559 + t799 + t643 * t698) * t654) * MDP(35); ((-t592 * t657 - t594 * t652) * t483 + t624 * t435 + (t652 * t708 - t657 * t710) * t559 + t715 * t470 - t853) * MDP(34) + (-(-t592 * t652 + t594 * t657) * t483 + t624 * t434 + (t652 * t710 + t657 * t708) * t559 - t715 * t698 - t852) * MDP(35) + (t486 - t793) * MDP(17) + (-t505 * t643 - t667) * MDP(20) + (-t504 * t643 - t666) * MDP(21) + (-pkin(4) * t462 - t802 - t505 * t534 - t507 * t839 + (t504 * t839 + t687) * t653 - t670 * t658) * MDP(27) + (-pkin(4) * t461 - t505 * t536 + t653 * t670 + t658 * t687 + t779 * t839 + t801) * MDP(28) + t857; t536 * t534 * MDP(22) + (-t534 ^ 2 + t536 ^ 2) * MDP(23) + (t461 + t856) * MDP(24) + (-t726 + (-qJD(5) + t839) * t536) * MDP(25) + t484 * MDP(26) + (-g(1) * t532 + t454 * t839 - t491 * t536 + t438 + (-qJD(5) * t492 - t810) * t658 + (t712 + t719) * t653) * MDP(27) + (t453 * t839 + t491 * t534 + g(1) * t533 - g(2) * (-t562 * t658 - t653 * t661) + g(3) * t789 + t683) * MDP(28) + (t434 + t845) * MDP(31) + (-t435 - t844) * MDP(32) + (-(-t444 * t652 - t803) * t559 - t432 * qJD(6) + (-t470 * t536 + t483 * t657 - t559 * t752) * pkin(5) + t837) * MDP(34) + ((-t445 * t559 - t428) * t652 + (t444 * t559 - t718) * t657 + (-t483 * t652 + t536 * t698 - t559 * t751) * pkin(5) + t838) * MDP(35) + t832; (t740 + t845) * MDP(31) + (-t727 - t844) * MDP(32) + (t432 * t559 + t837) * MDP(34) + (-t652 * t428 - t657 * t429 + t431 * t559 + t838) * MDP(35) + (-MDP(31) * t794 + MDP(32) * t698 - MDP(34) * t432 - MDP(35) * t804) * qJD(6) + t832;];
tau  = t1;
