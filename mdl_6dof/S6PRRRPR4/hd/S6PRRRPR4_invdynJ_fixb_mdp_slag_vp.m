% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:20:46
% EndTime: 2019-03-08 23:21:03
% DurationCPUTime: 12.16s
% Computational Cost: add. (6413->590), mult. (14839->814), div. (0->0), fcn. (11645->16), ass. (0->264)
t665 = cos(qJ(6));
t662 = sin(qJ(4));
t666 = cos(qJ(4));
t746 = t666 * qJD(3);
t663 = sin(qJ(3));
t758 = qJD(2) * t663;
t613 = t662 * t758 - t746;
t755 = qJD(3) * t662;
t615 = t666 * t758 + t755;
t656 = sin(pkin(12));
t658 = cos(pkin(12));
t698 = -t613 * t658 - t615 * t656;
t782 = t665 * t698;
t541 = t613 * t656 - t615 * t658;
t661 = sin(qJ(6));
t795 = t541 * t661;
t479 = t782 + t795;
t667 = cos(qJ(3));
t757 = qJD(2) * t667;
t638 = -qJD(4) + t757;
t631 = -qJD(6) + t638;
t797 = t479 * t631;
t750 = qJD(4) * t663;
t825 = -qJD(2) * t750 + qJDD(3);
t744 = qJD(2) * qJD(3);
t724 = t667 * t744;
t742 = qJDD(2) * t663;
t826 = t724 + t742;
t533 = qJD(4) * t746 + t662 * t825 + t666 * t826;
t650 = t667 * qJDD(2);
t818 = -t663 * t744 + t650;
t605 = qJDD(4) - t818;
t657 = sin(pkin(6));
t664 = sin(qJ(2));
t668 = cos(qJ(2));
t745 = qJD(1) * qJD(2);
t576 = qJDD(2) * pkin(8) + (qJDD(1) * t664 + t668 * t745) * t657;
t659 = cos(pkin(6));
t743 = qJDD(1) * t659;
t722 = t663 * t743;
t762 = qJD(1) * t657;
t618 = qJD(2) * pkin(8) + t664 * t762;
t786 = t659 * t667;
t820 = qJD(1) * t786 - t663 * t618;
t496 = qJDD(3) * pkin(9) + qJD(3) * t820 + t576 * t667 + t722;
t761 = qJD(1) * t663;
t636 = t659 * t761;
t570 = t667 * t618 + t636;
t562 = qJD(3) * pkin(9) + t570;
t621 = -pkin(3) * t667 - pkin(9) * t663 - pkin(2);
t760 = qJD(1) * t668;
t735 = t657 * t760;
t571 = qJD(2) * t621 - t735;
t510 = t562 * t666 + t571 * t662;
t709 = pkin(3) * t663 - pkin(9) * t667;
t617 = t709 * qJD(3);
t726 = t664 * t745;
t788 = t657 * t668;
t704 = -qJDD(1) * t788 + t657 * t726;
t519 = qJD(2) * t617 + qJDD(2) * t621 + t704;
t518 = t666 * t519;
t676 = -qJD(4) * t510 - t662 * t496 + t518;
t433 = pkin(4) * t605 - qJ(5) * t533 - qJD(5) * t615 + t676;
t534 = t662 * (qJD(3) * (qJD(4) + t757) + t742) - t825 * t666;
t749 = qJD(4) * t666;
t739 = t666 * t496 + t662 * t519 + t571 * t749;
t751 = qJD(4) * t662;
t685 = t562 * t751 - t739;
t435 = -qJ(5) * t534 - qJD(5) * t613 - t685;
t426 = t656 * t433 + t658 * t435;
t466 = -t533 * t656 - t534 * t658;
t424 = pkin(10) * t466 + t426;
t561 = -qJD(3) * pkin(3) - t820;
t527 = pkin(4) * t613 + qJD(5) + t561;
t472 = -pkin(5) * t698 + t527;
t801 = sin(pkin(11));
t714 = t801 * t668;
t802 = cos(pkin(11));
t717 = t802 * t664;
t589 = t659 * t717 + t714;
t719 = t657 * t802;
t549 = t589 * t667 - t663 * t719;
t715 = t801 * t664;
t716 = t802 * t668;
t591 = -t659 * t715 + t716;
t718 = t657 * t801;
t551 = t591 * t667 + t663 * t718;
t588 = -t659 * t716 + t715;
t590 = t659 * t714 + t717;
t789 = t657 * t664;
t596 = t659 * t663 + t667 * t789;
t651 = qJ(4) + pkin(12) + qJ(6);
t642 = sin(t651);
t643 = cos(t651);
t425 = t658 * t433 - t435 * t656;
t467 = t533 * t658 - t534 * t656;
t423 = pkin(5) * t605 - pkin(10) * t467 + t425;
t508 = -t562 * t662 + t666 * t571;
t488 = -qJ(5) * t615 + t508;
t470 = -pkin(4) * t638 + t488;
t489 = -qJ(5) * t613 + t510;
t787 = t658 * t489;
t445 = t656 * t470 + t787;
t823 = pkin(10) * t698;
t440 = t445 + t823;
t747 = qJD(6) * t661;
t712 = t661 * t423 - t440 * t747;
t841 = -t472 * t479 - g(1) * (-t551 * t643 - t590 * t642) - g(2) * (-t549 * t643 - t588 * t642) - g(3) * (-t596 * t643 + t642 * t788) - t665 * t424 - t712;
t601 = qJDD(6) + t605;
t824 = -t665 * t541 + t661 * t698;
t840 = t601 * MDP(25) + (-t479 ^ 2 + t824 ^ 2) * MDP(22) - t479 * MDP(21) * t824;
t754 = qJD(3) * t663;
t780 = t667 * t668;
t809 = pkin(8) * t662;
t839 = (-t662 * t780 + t664 * t666) * t762 - t666 * t617 - t754 * t809;
t838 = -(t662 * t664 + t666 * t780) * t762 + t662 * t617 + t621 * t749;
t798 = t824 * t631;
t781 = t666 * t667;
t640 = pkin(8) * t781;
t810 = pkin(4) * t663;
t696 = -qJ(5) * t781 + t810;
t748 = qJD(5) * t666;
t836 = -t663 * t748 + t696 * qJD(3) + (-t640 + (qJ(5) * t663 - t621) * t662) * qJD(4) - t839;
t784 = t663 * t666;
t835 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t784 + (-qJD(5) * t663 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t667) * t662 + t838;
t616 = t709 * qJD(2);
t599 = t666 * t616;
t660 = -qJ(5) - pkin(9);
t720 = qJD(4) * t660;
t834 = -qJD(2) * t696 + t666 * t720 - t599 + (-qJD(5) + t820) * t662;
t730 = t662 * t757;
t770 = t662 * t616 + t666 * t820;
t833 = -qJ(5) * t730 - t662 * t720 - t748 + t770;
t708 = g(1) * t590 + g(2) * t588;
t832 = -g(3) * t788 + t708;
t713 = t665 * t423 - t661 * t424;
t831 = -t472 * t824 - g(1) * (-t551 * t642 + t590 * t643) - g(2) * (-t549 * t642 + t588 * t643) - g(3) * (-t596 * t642 - t643 * t788) + t713;
t830 = pkin(10) * t541;
t607 = t656 * t666 + t658 * t662;
t686 = t607 * t667;
t828 = qJD(2) * t686 - t607 * qJD(4);
t697 = t656 * t662 - t658 * t666;
t827 = t638 * t697;
t707 = g(1) * t591 + g(2) * t589;
t680 = -g(3) * t789 - t707;
t776 = -t656 * t835 + t658 * t836;
t775 = t656 * t836 + t658 * t835;
t772 = t656 * t833 + t658 * t834;
t771 = t656 * t834 - t658 * t833;
t554 = -t596 * t662 - t666 * t788;
t819 = -g(1) * (-t551 * t662 + t590 * t666) - g(2) * (-t549 * t662 + t588 * t666) - g(3) * t554;
t817 = -t570 + (-t730 + t751) * pkin(4);
t816 = pkin(4) * t534 + qJDD(5);
t669 = qJD(3) ^ 2;
t815 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t669 + t657 * (-g(3) * t668 + t726) - t704 + t708;
t813 = qJD(4) * (pkin(8) * t638 + t562) + t832;
t711 = -t665 * t466 + t467 * t661;
t430 = qJD(6) * t824 + t711;
t811 = pkin(4) * t656;
t803 = qJD(2) * pkin(2);
t799 = qJDD(3) * pkin(3);
t796 = t533 * t662;
t794 = t613 * t638;
t793 = t615 * t638;
t792 = t615 * t666;
t791 = t642 * t667;
t790 = t643 * t667;
t483 = t656 * t489;
t785 = t662 * t663;
t444 = t658 * t470 - t483;
t437 = -pkin(5) * t638 + t444 + t830;
t783 = t665 * t437;
t779 = qJDD(1) - g(3);
t727 = t667 * t746;
t753 = qJD(3) * t667;
t728 = t662 * t753;
t524 = t607 * t750 + t656 * t728 - t658 * t727;
t778 = -pkin(5) * t754 - pkin(10) * t524 - t776;
t523 = -qJD(3) * t686 + t697 * t750;
t777 = pkin(10) * t523 + t775;
t699 = -t607 * t661 - t665 * t697;
t774 = qJD(6) * t699 + t661 * t828 + t665 * t827;
t540 = t607 * t665 - t661 * t697;
t773 = qJD(6) * t540 + t661 * t827 - t665 * t828;
t451 = t658 * t488 - t483;
t609 = t666 * t621;
t545 = -qJ(5) * t784 + t609 + (-pkin(4) - t809) * t667;
t766 = t662 * t621 + t640;
t556 = -qJ(5) * t785 + t766;
t491 = t656 * t545 + t658 * t556;
t769 = -pkin(5) * t828 + t817;
t622 = t660 * t662;
t623 = t660 * t666;
t559 = t656 * t622 - t658 * t623;
t765 = pkin(4) * t785 + t663 * pkin(8);
t654 = t663 ^ 2;
t764 = -t667 ^ 2 + t654;
t759 = qJD(2) * t657;
t756 = qJD(3) * t613;
t752 = qJD(4) * t638;
t740 = qJD(6) * t782 + t661 * t466 + t665 * t467;
t737 = pkin(4) * t728 + pkin(8) * t753 + t749 * t810;
t646 = pkin(4) * t666 + pkin(3);
t733 = t663 * t760;
t732 = t664 * t759;
t731 = t668 * t759;
t729 = t638 * t746;
t450 = -t488 * t656 - t787;
t490 = t658 * t545 - t556 * t656;
t558 = t658 * t622 + t623 * t656;
t529 = -pkin(10) * t697 + t559;
t706 = pkin(5) * t758 + pkin(10) * t827 + qJD(6) * t529 - t772;
t528 = -pkin(10) * t607 + t558;
t705 = pkin(10) * t828 + qJD(6) * t528 + t771;
t428 = t661 * t437 + t665 * t440;
t583 = t697 * t663;
t459 = -pkin(5) * t667 + pkin(10) * t583 + t490;
t582 = t607 * t663;
t460 = -pkin(10) * t582 + t491;
t703 = t459 * t661 + t460 * t665;
t690 = -t596 * t666 + t662 * t788;
t493 = t554 * t658 + t656 * t690;
t494 = t554 * t656 - t658 * t690;
t702 = t493 * t665 - t494 * t661;
t701 = t493 * t661 + t494 * t665;
t700 = -t665 * t582 + t583 * t661;
t522 = -t582 * t661 - t583 * t665;
t644 = pkin(4) * t658 + pkin(5);
t695 = t644 * t661 + t665 * t811;
t694 = t644 * t665 - t661 * t811;
t595 = t663 * t789 - t786;
t688 = t605 * t662 - t638 * t749;
t687 = t605 * t666 + t638 * t751;
t429 = t541 * t747 + t740;
t548 = t589 * t663 + t667 * t719;
t550 = t591 * t663 - t667 * t718;
t684 = g(1) * t550 + g(2) * t548 + g(3) * t595;
t683 = g(1) * t551 + g(2) * t549 + g(3) * t596;
t682 = qJD(3) * t636 + t663 * t576 + t618 * t753 - t667 * t743;
t679 = -pkin(9) * t605 - t561 * t638;
t497 = t682 - t799;
t675 = pkin(9) * t752 - t497 + t684;
t619 = -t735 - t803;
t674 = -pkin(8) * qJDD(3) + (t619 + t735 - t803) * qJD(3);
t458 = t497 + t816;
t672 = -t682 + t684;
t670 = qJD(2) ^ 2;
t577 = pkin(5) * t697 - t646;
t553 = qJD(3) * t596 + t663 * t731;
t552 = -qJD(3) * t595 + t667 * t731;
t546 = pkin(5) * t582 + t765;
t514 = pkin(4) * t615 - pkin(5) * t541;
t492 = -pkin(5) * t523 + t737;
t482 = qJD(4) * t554 + t552 * t666 + t662 * t732;
t481 = qJD(4) * t690 - t552 * t662 + t666 * t732;
t455 = qJD(6) * t522 - t665 * t523 - t524 * t661;
t454 = qJD(6) * t700 + t523 * t661 - t524 * t665;
t449 = t481 * t656 + t482 * t658;
t447 = t481 * t658 - t482 * t656;
t443 = -pkin(5) * t466 + t458;
t442 = t451 + t830;
t441 = t450 - t823;
t427 = -t440 * t661 + t783;
t1 = [t779 * MDP(1) + (-qJD(3) * t553 - qJDD(3) * t595) * MDP(10) + (-qJD(3) * t552 - qJDD(3) * t596) * MDP(11) + (-t481 * t638 + t534 * t595 + t553 * t613 + t554 * t605) * MDP(17) + (t482 * t638 + t533 * t595 + t553 * t615 + t605 * t690) * MDP(18) + (t447 * t541 + t449 * t698 + t466 * t494 - t467 * t493) * MDP(19) + (t425 * t493 + t426 * t494 + t444 * t447 + t445 * t449 + t458 * t595 + t527 * t553 - g(3)) * MDP(20) + (-(-qJD(6) * t701 + t447 * t665 - t449 * t661) * t631 + t702 * t601 - t553 * t479 + t595 * t430) * MDP(26) + ((qJD(6) * t702 + t447 * t661 + t449 * t665) * t631 - t701 * t601 + t553 * t824 + t595 * t429) * MDP(27) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t667 + MDP(11) * t663 - MDP(3)) * t670) * t664 + (t818 * MDP(10) - MDP(11) * t826 + qJDD(2) * MDP(3) - t670 * MDP(4)) * t668) * t657; (t674 * t663 + t667 * t815) * MDP(10) + (-t663 * t815 + t674 * t667) * MDP(11) + (-t703 * t601 + ((qJD(6) * t437 + t424) * t665 + t712) * t667 - t428 * t754 + t492 * t824 + t546 * t429 + t443 * t522 + t472 * t454 - g(1) * (t590 * t791 + t591 * t643) - g(2) * (t588 * t791 + t589 * t643) + ((qJD(6) * t459 + t777) * t665 + (-qJD(6) * t460 - t778) * t661) * t631 + (-t824 * t733 - g(3) * (-t642 * t780 + t643 * t664)) * t657) * MDP(27) + (-t429 * t667 - t454 * t631 + t522 * t601 + t754 * t824) * MDP(23) + (t429 * t522 + t454 * t824) * MDP(21) + (t426 * t491 + t425 * t490 + t458 * t765 + t775 * t445 + t776 * t444 + t680 * (pkin(4) * t662 + pkin(8)) + (-t761 * t788 + t737) * t527 + t832 * (t646 * t667 - t660 * t663 + pkin(2))) * MDP(20) + ((-t613 * t666 - t615 * t662) * t753 + (-t796 - t534 * t666 + (t613 * t662 - t792) * qJD(4)) * t663) * MDP(13) + (qJDD(2) * t654 + 0.2e1 * t663 * t724) * MDP(5) + (t533 * t784 + (-t662 * t750 + t727) * t615) * MDP(12) + (t779 * t788 + t708) * MDP(3) + (-t766 * t605 + t838 * t638 + t680 * t666 + ((pkin(8) * t615 + t561 * t666) * qJD(3) - t813 * t662 + t739) * t667 + (-t615 * t735 - t561 * t751 - t510 * qJD(3) + t497 * t666 + (t533 - t729) * pkin(8)) * t663) * MDP(18) + (t609 * t605 + t839 * t638 + (t621 * t752 + t680) * t662 + (pkin(8) * t756 - t518 + (-pkin(8) * t605 + qJD(3) * t561 + qJD(4) * t571 + t496) * t662 + t813 * t666) * t667 + (pkin(8) * t534 + qJD(3) * t508 + t497 * t662 + t561 * t749 - t613 * t735) * t663) * MDP(17) + 0.2e1 * (t650 * t663 - t744 * t764) * MDP(6) + ((t459 * t665 - t460 * t661) * t601 - t713 * t667 + t427 * t754 - t492 * t479 + t546 * t430 - t443 * t700 + t472 * t455 - g(1) * (-t590 * t790 + t591 * t642) - g(2) * (-t588 * t790 + t589 * t642) + (t661 * t777 + t665 * t778) * t631 + (t428 * t667 + t631 * t703) * qJD(6) + (t479 * t733 - g(3) * (t642 * t664 + t643 * t780)) * t657) * MDP(26) + (t430 * t667 + t455 * t631 + t479 * t754 + t601 * t700) * MDP(24) + (t429 * t700 - t430 * t522 + t454 * t479 - t455 * t824) * MDP(22) + (t425 * t583 - t426 * t582 + t444 * t524 + t445 * t523 + t466 * t491 - t467 * t490 + t541 * t776 + t663 * t832 + t698 * t775) * MDP(19) + (-t779 * t789 + t707) * MDP(4) + (-t601 * t667 - t631 * t754) * MDP(25) + ((t638 * t755 + t534) * t667 + (-t688 - t756) * t663) * MDP(15) + (-t605 * t667 - t638 * t754) * MDP(16) + qJDD(2) * MDP(2) + (qJDD(3) * t663 + t667 * t669) * MDP(7) + (qJDD(3) * t667 - t663 * t669) * MDP(8) + ((-t533 - t729) * t667 + (qJD(3) * t615 + t687) * t663) * MDP(14); MDP(7) * t742 + MDP(8) * t650 + qJDD(3) * MDP(9) + (qJD(3) * t570 + t672) * MDP(10) + (-t722 + (-qJD(2) * t619 - t576) * t667 + t683) * MDP(11) + (-t638 * t792 + t796) * MDP(12) + ((t533 + t794) * t666 + (-t534 + t793) * t662) * MDP(13) + ((-t615 * t663 + t638 * t781) * qJD(2) + t688) * MDP(14) + ((-t638 * t662 * t667 + t613 * t663) * qJD(2) + t687) * MDP(15) + (-pkin(3) * t534 - t570 * t613 + t599 * t638 + (-t638 * t820 + t679) * t662 + t675 * t666) * MDP(17) + (-pkin(3) * t533 - t570 * t615 - t638 * t770 - t662 * t675 + t666 * t679) * MDP(18) + (-t425 * t607 - t426 * t697 - t444 * t827 + t445 * t828 + t466 * t559 - t467 * t558 + t772 * t541 + t771 * t698 - t683) * MDP(19) + (t426 * t559 + t425 * t558 - t458 * t646 - g(1) * (-t550 * t646 - t551 * t660) - g(2) * (-t548 * t646 - t549 * t660) - g(3) * (-t595 * t646 - t596 * t660) + t817 * t527 + t771 * t445 + t772 * t444) * MDP(20) + (t429 * t540 + t774 * t824) * MDP(21) + (t429 * t699 - t430 * t540 + t479 * t774 - t773 * t824) * MDP(22) + (t540 * t601 - t631 * t774) * MDP(23) + (t601 * t699 + t631 * t773) * MDP(24) + ((t528 * t665 - t529 * t661) * t601 + t577 * t430 - t443 * t699 + (t661 * t705 + t665 * t706) * t631 - t769 * t479 + t773 * t472 + t684 * t643) * MDP(26) + (-(t528 * t661 + t529 * t665) * t601 + t577 * t429 + t443 * t540 + (-t661 * t706 + t665 * t705) * t631 + t769 * t824 + t774 * t472 - t684 * t642) * MDP(27) + (-MDP(10) * t619 + t638 * MDP(16) - t508 * MDP(17) + t510 * MDP(18) - MDP(23) * t824 - MDP(24) * t479 + t631 * MDP(25) - t427 * MDP(26) + t428 * MDP(27)) * t758 + (-MDP(5) * t663 * t667 + MDP(6) * t764) * t670; t615 * t613 * MDP(12) + (-t613 ^ 2 + t615 ^ 2) * MDP(13) + (t533 - t794) * MDP(14) + (-t534 - t793) * MDP(15) + t605 * MDP(16) + (-t510 * t638 - t561 * t615 + t676 + t819) * MDP(17) + (-t508 * t638 + t561 * t613 - g(1) * (-t551 * t666 - t590 * t662) - g(2) * (-t549 * t666 - t588 * t662) - g(3) * t690 + t685) * MDP(18) + ((t466 * t656 - t467 * t658) * pkin(4) + (t444 - t451) * t698 + (-t445 - t450) * t541) * MDP(19) + (-t444 * t450 - t445 * t451 + (t425 * t658 + t426 * t656 - t527 * t615 + t819) * pkin(4)) * MDP(20) + (t429 + t797) * MDP(23) + (-t430 - t798) * MDP(24) + (t694 * t601 + (t441 * t665 - t442 * t661) * t631 + t514 * t479 + (t631 * t695 - t428) * qJD(6) + t831) * MDP(26) + (-t695 * t601 - (t441 * t661 + t442 * t665) * t631 - t514 * t824 + (t631 * t694 - t783) * qJD(6) + t841) * MDP(27) + t840; (-t541 ^ 2 - t698 ^ 2) * MDP(19) + (-t444 * t541 - t445 * t698 - t672 - t799 + t816) * MDP(20) + (t430 - t798) * MDP(26) + (t429 - t797) * MDP(27); (t740 + t797) * MDP(23) + (-t711 - t798) * MDP(24) + (-t428 * t631 + t831) * MDP(26) + (-t427 * t631 + t841) * MDP(27) + (MDP(23) * t795 - MDP(24) * t824 - MDP(26) * t428 - MDP(27) * t783) * qJD(6) + t840;];
tau  = t1;
