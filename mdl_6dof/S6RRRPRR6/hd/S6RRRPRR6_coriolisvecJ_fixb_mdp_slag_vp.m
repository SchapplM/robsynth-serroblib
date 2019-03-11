% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:29:52
% EndTime: 2019-03-09 18:30:15
% DurationCPUTime: 14.36s
% Computational Cost: add. (11640->576), mult. (29102->786), div. (0->0), fcn. (22187->10), ass. (0->259)
t707 = sin(qJ(2));
t711 = cos(qJ(2));
t671 = -pkin(2) * t711 - pkin(8) * t707 - pkin(1);
t655 = t671 * qJD(1);
t784 = qJD(1) * t711;
t695 = pkin(7) * t784;
t678 = qJD(2) * pkin(8) + t695;
t706 = sin(qJ(3));
t710 = cos(qJ(3));
t618 = t710 * t655 - t678 * t706;
t783 = qJD(2) * t706;
t785 = qJD(1) * t707;
t667 = t710 * t785 + t783;
t582 = -qJ(4) * t667 + t618;
t689 = -qJD(3) + t784;
t571 = -pkin(3) * t689 + t582;
t820 = t706 * t655;
t619 = t678 * t710 + t820;
t765 = t706 * t785;
t772 = t710 * qJD(2);
t665 = t765 - t772;
t583 = -qJ(4) * t665 + t619;
t702 = sin(pkin(11));
t576 = t702 * t583;
t703 = cos(pkin(11));
t510 = t703 * t571 - t576;
t610 = t665 * t702 - t667 * t703;
t855 = pkin(9) * t610;
t490 = -pkin(4) * t689 + t510 + t855;
t821 = t703 * t583;
t511 = t702 * t571 + t821;
t728 = -t665 * t703 - t667 * t702;
t844 = pkin(9) * t728;
t492 = t511 + t844;
t705 = sin(qJ(5));
t709 = cos(qJ(5));
t451 = t490 * t705 + t492 * t709;
t551 = t610 * t705 + t709 * t728;
t865 = pkin(10) * t551;
t447 = t451 + t865;
t704 = sin(qJ(6));
t774 = qJD(6) * t704;
t443 = t447 * t774;
t708 = cos(qJ(6));
t847 = -t709 * t610 + t705 * t728;
t828 = t847 * t704;
t487 = t551 * t708 - t828;
t677 = -qJD(2) * pkin(2) + pkin(7) * t785;
t623 = pkin(3) * t665 + qJD(4) + t677;
t560 = -pkin(4) * t728 + t623;
t499 = -pkin(5) * t551 + t560;
t874 = -t499 * t487 + t443;
t483 = t551 * t704 + t708 * t847;
t758 = MDP(31) * t785;
t873 = qJD(2) * t758 + (t483 ^ 2 - t487 ^ 2) * MDP(28) - t487 * MDP(27) * t483;
t771 = qJD(1) * qJD(2);
t759 = t711 * t771;
t779 = qJD(3) * t707;
t762 = t706 * t779;
t770 = qJD(2) * qJD(3);
t625 = -qJD(1) * t762 + (t759 + t770) * t710;
t778 = qJD(3) * t710;
t761 = t707 * t778;
t781 = qJD(2) * t711;
t764 = t706 * t781;
t849 = t761 + t764;
t626 = qJD(1) * t849 + t706 * t770;
t568 = -t625 * t702 - t626 * t703;
t569 = t625 * t703 - t626 * t702;
t775 = qJD(5) * t709;
t776 = qJD(5) * t705;
t477 = t705 * t568 + t709 * t569 + t610 * t776 + t728 * t775;
t738 = pkin(2) * t707 - pkin(8) * t711;
t669 = t738 * qJD(2);
t656 = qJD(1) * t669;
t760 = t707 * t771;
t741 = pkin(7) * t760;
t796 = -t710 * t656 - t706 * t741;
t718 = -qJD(3) * t619 - t796;
t509 = pkin(3) * t760 - qJ(4) * t625 - qJD(4) * t667 + t718;
t780 = qJD(3) * t706;
t725 = t655 * t778 + t706 * t656 - t678 * t780;
t714 = -t710 * t741 + t725;
t517 = -qJ(4) * t626 - qJD(4) * t665 + t714;
t466 = t703 * t509 - t517 * t702;
t457 = pkin(4) * t760 - pkin(9) * t569 + t466;
t467 = t702 * t509 + t703 * t517;
t459 = pkin(9) * t568 + t467;
t751 = t709 * t457 - t705 * t459;
t716 = -qJD(5) * t451 + t751;
t434 = pkin(5) * t760 - pkin(10) * t477 + t716;
t478 = qJD(5) * t847 - t709 * t568 + t569 * t705;
t740 = -t705 * t457 - t709 * t459 - t490 * t775 + t492 * t776;
t435 = -pkin(10) * t478 - t740;
t872 = -t704 * t434 - t708 * t435 + t874;
t773 = qJD(6) * t708;
t768 = t708 * t477 - t704 * t478 + t551 * t773;
t440 = -t774 * t847 + t768;
t749 = t477 * t704 + t708 * t478;
t715 = -qJD(6) * t483 - t749;
t756 = MDP(24) * t785;
t680 = -qJD(5) + t689;
t829 = t551 * t680;
t830 = t847 * t680;
t672 = -qJD(6) + t680;
t869 = t487 * t672;
t870 = t483 * t672;
t871 = qJD(2) * t756 + (-t478 - t830) * MDP(23) + (-t551 ^ 2 + t847 ^ 2) * MDP(21) - t551 * MDP(20) * t847 + (t477 + t829) * MDP(22) + (t440 + t869) * MDP(29) + (t715 - t870) * MDP(30) + t873;
t813 = t710 * t711;
t724 = pkin(3) * t707 - qJ(4) * t813;
t833 = -qJ(4) - pkin(8);
t754 = qJD(3) * t833;
t668 = t738 * qJD(1);
t793 = pkin(7) * t765 + t710 * t668;
t868 = -qJD(1) * t724 - qJD(4) * t706 + t710 * t754 - t793;
t651 = t706 * t668;
t777 = qJD(4) * t710;
t817 = t707 * t710;
t818 = t706 * t711;
t867 = t651 + (-pkin(7) * t817 - qJ(4) * t818) * qJD(1) - t706 * t754 - t777;
t752 = t708 * t434 - t704 * t435;
t848 = -t499 * t483 + t752;
t659 = t702 * t710 + t703 * t706;
t722 = t659 * t711;
t857 = qJD(1) * t722 - t659 * qJD(3);
t727 = t702 * t706 - t703 * t710;
t864 = t689 * t727;
t799 = t702 * t867 + t703 * t868;
t798 = t702 * t868 - t703 * t867;
t863 = -t560 * t551 + t740;
t860 = pkin(10) * t847;
t858 = pkin(4) * t785 + pkin(9) * t864 - t799;
t851 = pkin(9) * t857 + t798;
t856 = -t560 * t847 + t716;
t729 = -t659 * t705 - t709 * t727;
t802 = qJD(5) * t729 + t705 * t857 + t709 * t864;
t608 = t659 * t709 - t705 * t727;
t801 = qJD(5) * t608 + t705 * t864 - t709 * t857;
t850 = -t695 + (-t706 * t784 + t780) * pkin(3);
t845 = -0.2e1 * t771;
t843 = MDP(4) * t707;
t700 = t707 ^ 2;
t842 = MDP(5) * (-t711 ^ 2 + t700);
t839 = t858 * t709;
t661 = t710 * t671;
t834 = pkin(7) * t706;
t615 = -qJ(4) * t817 + t661 + (-pkin(3) - t834) * t711;
t691 = pkin(7) * t813;
t791 = t706 * t671 + t691;
t819 = t706 * t707;
t620 = -qJ(4) * t819 + t791;
t554 = t703 * t615 - t620 * t702;
t640 = t727 * t707;
t527 = -pkin(4) * t711 + pkin(9) * t640 + t554;
t555 = t702 * t615 + t703 * t620;
t639 = t659 * t707;
t531 = -pkin(9) * t639 + t555;
t803 = t705 * t527 + t709 * t531;
t674 = t833 * t706;
t675 = t833 * t710;
t621 = t703 * t674 + t675 * t702;
t596 = -pkin(9) * t659 + t621;
t622 = t702 * t674 - t703 * t675;
t597 = -pkin(9) * t727 + t622;
t800 = t705 * t596 + t709 * t597;
t797 = -pkin(4) * t857 + t850;
t520 = -t582 * t702 - t821;
t500 = t520 - t844;
t521 = t703 * t582 - t576;
t501 = t521 + t855;
t692 = pkin(3) * t703 + pkin(4);
t835 = pkin(3) * t702;
t739 = t709 * t692 - t705 * t835;
t838 = t739 * qJD(5) - t705 * t500 - t709 * t501;
t645 = t692 * t705 + t709 * t835;
t837 = t645 * qJD(5) + t709 * t500 - t501 * t705;
t836 = t596 * t775 - t597 * t776 - t705 * t858 + t709 * t851;
t827 = t625 * t706;
t826 = t665 * t689;
t825 = t667 * t689;
t824 = t677 * t706;
t823 = t677 * t710;
t822 = t689 * t710;
t712 = qJD(2) ^ 2;
t816 = t707 * t712;
t450 = t709 * t490 - t492 * t705;
t446 = t450 - t860;
t442 = -pkin(5) * t680 + t446;
t815 = t708 * t442;
t814 = t708 * t447;
t812 = t711 * t712;
t713 = qJD(1) ^ 2;
t811 = t711 * t713;
t810 = t837 - t865;
t809 = t838 + t860;
t545 = t608 * t704 - t708 * t729;
t808 = -qJD(6) * t545 - t704 * t801 + t708 * t802;
t546 = t608 * t708 + t704 * t729;
t807 = qJD(6) * t546 + t704 * t802 + t708 * t801;
t804 = pkin(5) * t801 + t797;
t782 = qJD(2) * t707;
t794 = t710 * t669 + t782 * t834;
t542 = -t707 * t777 + t724 * qJD(2) + (-t691 + (qJ(4) * t707 - t671) * t706) * qJD(3) + t794;
t795 = t706 * t669 + t671 * t778;
t553 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t817 + (-qJD(4) * t707 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t711) * t706 + t795;
t481 = t702 * t542 + t703 * t553;
t789 = pkin(3) * t819 + t707 * pkin(7);
t767 = pkin(3) * t849 + pkin(7) * t781;
t766 = -pkin(3) * t710 - pkin(2);
t763 = t711 * t772;
t755 = MDP(15) * t782;
t609 = pkin(3) * t626 + pkin(7) * t759;
t753 = pkin(1) * t845;
t480 = t703 * t542 - t553 * t702;
t591 = t659 * t779 + t702 * t764 - t703 * t763;
t470 = pkin(4) * t782 + pkin(9) * t591 + t480;
t590 = -qJD(2) * t722 + t727 * t779;
t476 = pkin(9) * t590 + t481;
t750 = t709 * t470 - t476 * t705;
t747 = t709 * t527 - t531 * t705;
t745 = t709 * t596 - t597 * t705;
t744 = t665 + t772;
t743 = -t667 + t783;
t742 = qJD(6) * t442 + t435;
t616 = pkin(4) * t639 + t789;
t573 = pkin(3) * t667 - pkin(4) * t610;
t502 = -pkin(10) * t608 + t745;
t737 = -pkin(10) * t801 + qJD(6) * t502 + t836;
t503 = pkin(10) * t729 + t800;
t736 = pkin(5) * t785 + pkin(10) * t802 + t800 * qJD(5) + qJD(6) * t503 + t705 * t851 + t839;
t437 = t704 * t442 + t814;
t589 = -t639 * t705 - t640 * t709;
t462 = -pkin(5) * t711 - pkin(10) * t589 + t747;
t732 = -t709 * t639 + t640 * t705;
t465 = pkin(10) * t732 + t803;
t735 = t462 * t704 + t465 * t708;
t523 = t589 * t704 - t708 * t732;
t524 = t589 * t708 + t704 * t732;
t641 = pkin(5) + t739;
t731 = t641 * t708 - t645 * t704;
t730 = t641 * t704 + t645 * t708;
t726 = qJD(1) * t700 - t689 * t711;
t631 = pkin(4) * t727 + t766;
t556 = -pkin(4) * t590 + t767;
t528 = -pkin(4) * t568 + t609;
t721 = t705 * t470 + t709 * t476 + t527 * t775 - t531 * t776;
t570 = -pkin(5) * t729 + t631;
t538 = -pkin(5) * t732 + t616;
t504 = pkin(5) * t847 + t573;
t498 = qJD(5) * t589 - t709 * t590 - t591 * t705;
t497 = qJD(5) * t732 + t590 * t705 - t591 * t709;
t471 = pkin(5) * t498 + t556;
t454 = pkin(5) * t478 + t528;
t445 = qJD(6) * t524 + t497 * t704 + t708 * t498;
t444 = -qJD(6) * t523 + t497 * t708 - t498 * t704;
t439 = -pkin(10) * t498 + t721;
t438 = pkin(5) * t782 - pkin(10) * t497 - qJD(5) * t803 + t750;
t436 = -t447 * t704 + t815;
t1 = [t842 * t845 + (t477 * t732 - t478 * t589 + t497 * t551 - t498 * t847) * MDP(21) + (-t750 * t680 - t751 * t711 - t556 * t551 + t616 * t478 - t528 * t732 + t560 * t498 + (t451 * t711 + t680 * t803) * qJD(5)) * MDP(25) + (t445 * t672 - t711 * t715) * MDP(30) + (t616 * t477 + t560 * t497 + t528 * t589 + t556 * t847 + t721 * t680 - t740 * t711) * MDP(26) + (t477 * t589 + t497 * t847) * MDP(20) + 0.2e1 * t759 * t843 + (t625 * t817 + (-t762 + t763) * t667) * MDP(11) - MDP(7) * t816 + (pkin(7) * t816 + t711 * t753) * MDP(10) + (-pkin(7) * t812 + t707 * t753) * MDP(9) + (t466 * t554 + t467 * t555 + t510 * t480 + t511 * t481 + t609 * t789 + t623 * t767) * MDP(19) + (-t689 - t784) * t755 + (t466 * t640 - t467 * t639 + t480 * t610 + t481 * t728 + t510 * t591 + t511 * t590 - t554 * t569 + t555 * t568) * MDP(18) + MDP(6) * t812 + (-(t438 * t708 - t439 * t704) * t672 - t752 * t711 - t471 * t487 - t538 * t715 + t454 * t523 + t499 * t445 + (t437 * t711 + t672 * t735) * qJD(6)) * MDP(32) + (t538 * t440 - t443 * t711 + t499 * t444 + t454 * t524 + t471 * t483 + ((-qJD(6) * t465 + t438) * t672 + t434 * t711) * t704 + ((qJD(6) * t462 + t439) * t672 + t742 * t711) * t708) * MDP(33) + (t440 * t524 + t444 * t483) * MDP(27) + ((-qJD(1) * t803 - t451) * MDP(26) + ((t462 * t708 - t465 * t704) * qJD(1) + t436) * MDP(32) + (qJD(1) * t589 + t847) * MDP(22) + (qJD(1) * t732 + t551) * MDP(23) + (qJD(1) * t524 + t483) * MDP(29) + (-qJD(1) * t523 + t487) * MDP(30) + (qJD(1) * t747 + t450) * MDP(25) + (-qJD(1) * t735 - t437) * MDP(33) + (-t680 - t784) * MDP(24) + (-t672 - t784) * MDP(31)) * t782 + (-t440 * t523 + t444 * t487 - t445 * t483 + t524 * t715) * MDP(28) + ((-pkin(7) * t711 * t780 + t795) * t689 + t725 * t711 + (pkin(7) * t625 - t677 * t780) * t707 + ((pkin(7) * t667 + t823) * t711 + (-pkin(7) * t822 - qJD(1) * t791 - t619) * t707) * qJD(2)) * MDP(17) + (-(-t671 * t780 + t794) * t689 + (t677 * t778 + pkin(7) * t626 + (qJD(1) * t661 + t618) * qJD(2)) * t707 + ((pkin(7) * t665 + t824) * qJD(2) + (t820 + (pkin(7) * t689 + t678) * t710) * qJD(3) + t796) * t711) * MDP(16) + (t478 * t711 + t498 * t680) * MDP(23) + (-t477 * t711 - t497 * t680) * MDP(22) + (-t440 * t711 - t444 * t672) * MDP(29) + (t689 * t761 + t626 * t711 + (-t665 * t707 - t706 * t726) * qJD(2)) * MDP(14) + (t689 * t762 - t625 * t711 + (t667 * t707 + t710 * t726) * qJD(2)) * MDP(13) + ((-t665 * t710 - t667 * t706) * t781 + (-t827 - t626 * t710 + (t665 * t706 - t667 * t710) * qJD(3)) * t707) * MDP(12); (MDP(9) * t707 * t713 + MDP(10) * t811) * pkin(1) + t713 * t842 + t689 * MDP(15) * t785 + (t631 * t478 - t528 * t729 + (t597 * t775 + (qJD(5) * t596 + t851) * t705 + t839) * t680 + t801 * t560 - t797 * t551 + (qJD(2) * t745 - t450) * t785) * MDP(25) + (t477 * t729 - t478 * t608 + t551 * t802 - t801 * t847) * MDP(21) + (t801 * t680 + (qJD(2) * t729 - t551) * t785) * MDP(23) + (t631 * t477 + t528 * t608 + t836 * t680 + t802 * t560 + t797 * t847 + (-qJD(2) * t800 + t451) * t785) * MDP(26) + (t477 * t608 + t802 * t847) * MDP(20) + (-t802 * t680 + (qJD(2) * t608 - t847) * t785) * MDP(22) - t811 * t843 + (t689 * t780 + (-t689 * t818 + t707 * t744) * qJD(1)) * MDP(14) + (-t689 * t778 + (t689 * t813 + t707 * t743) * qJD(1)) * MDP(13) + t680 * t756 + t672 * t758 + (t466 * t621 + t467 * t622 + t799 * t510 + t798 * t511 + t609 * t766 + t623 * t850) * MDP(19) + (-t466 * t659 - t467 * t727 - t510 * t864 + t857 * t511 + t568 * t622 - t569 * t621 + t799 * t610 + t798 * t728) * MDP(18) + (-t570 * t715 + t454 * t545 + (t704 * t737 + t708 * t736) * t672 + t807 * t499 - t804 * t487 + ((t502 * t708 - t503 * t704) * qJD(2) - t436) * t785) * MDP(32) + (t807 * t672 + (-qJD(2) * t545 - t487) * t785) * MDP(30) + (t440 * t546 + t483 * t808) * MDP(27) + (-t808 * t672 + (qJD(2) * t546 - t483) * t785) * MDP(29) + (t570 * t440 + t454 * t546 + (-t704 * t736 + t708 * t737) * t672 + t808 * t499 + t804 * t483 + (-(t502 * t704 + t503 * t708) * qJD(2) + t437) * t785) * MDP(33) + (-t440 * t545 - t483 * t807 + t487 * t808 + t546 * t715) * MDP(28) + (-pkin(2) * t625 - t651 * t689 + (-pkin(8) * t689 * t706 + t823) * qJD(3) + (-t677 * t813 + (-pkin(8) * t772 + t619) * t707 + (t689 * t817 + t711 * t743) * pkin(7)) * qJD(1)) * MDP(17) + (-pkin(2) * t626 + t793 * t689 + (pkin(8) * t822 + t824) * qJD(3) + ((-pkin(8) * t783 - t618) * t707 + (-pkin(7) * t744 - t824) * t711) * qJD(1)) * MDP(16) + ((t625 + t826) * t710 + (-t626 + t825) * t706) * MDP(12) + (-t667 * t822 + t827) * MDP(11); (t551 * t573 + t680 * t837 + t739 * t760 + t856) * MDP(25) + (-t573 * t847 - t645 * t760 + t680 * t838 + t863) * MDP(26) + (-t626 - t825) * MDP(14) + (t731 * t760 + t504 * t487 + (t704 * t809 + t708 * t810) * t672 + (t672 * t730 - t437) * qJD(6) + t848) * MDP(32) + qJD(1) * t755 + (-t730 * t760 - t504 * t483 + (-t704 * t810 + t708 * t809) * t672 + (t672 * t731 - t815) * qJD(6) + t872) * MDP(33) + (-t665 ^ 2 + t667 ^ 2) * MDP(12) + (-t619 * t689 - t667 * t677 + t718) * MDP(16) + (-t618 * t689 + t665 * t677 - t714) * MDP(17) + ((t568 * t702 - t569 * t703) * pkin(3) + (t510 - t521) * t728 + (-t511 - t520) * t610) * MDP(18) + t667 * t665 * MDP(11) + (-t510 * t520 - t511 * t521 + (t466 * t703 + t467 * t702 - t623 * t667) * pkin(3)) * MDP(19) + (t625 - t826) * MDP(13) + t871; (-t610 ^ 2 - t728 ^ 2) * MDP(18) + (-t510 * t610 - t511 * t728 + t609) * MDP(19) + (t478 - t830) * MDP(25) + (t477 - t829) * MDP(26) + (-t715 - t870) * MDP(32) + (t440 - t869) * MDP(33); (-t451 * t680 + t856) * MDP(25) + (-t450 * t680 + t863) * MDP(26) + ((-t446 * t704 - t814) * t672 - t437 * qJD(6) + (t487 * t847 + t672 * t774 + t708 * t760) * pkin(5) + t848) * MDP(32) + ((t447 * t672 - t434) * t704 + (-t446 * t672 - t742) * t708 + (-t483 * t847 + t672 * t773 - t704 * t760) * pkin(5) + t874) * MDP(33) + t871; (t768 + t869) * MDP(29) + (-t749 - t870) * MDP(30) + (-t437 * t672 + t848) * MDP(32) + (-t436 * t672 + t872) * MDP(33) + (-MDP(29) * t828 - t483 * MDP(30) - MDP(32) * t437 - MDP(33) * t815) * qJD(6) + t873;];
tauc  = t1;
