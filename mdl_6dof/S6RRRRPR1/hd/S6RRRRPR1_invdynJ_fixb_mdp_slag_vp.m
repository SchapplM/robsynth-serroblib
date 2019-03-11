% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:55:16
% EndTime: 2019-03-09 21:55:31
% DurationCPUTime: 10.95s
% Computational Cost: add. (14954->557), mult. (36703->717), div. (0->0), fcn. (28426->18), ass. (0->271)
t714 = cos(qJ(6));
t805 = qJD(6) * t714;
t712 = sin(qJ(2));
t857 = cos(qJ(3));
t792 = qJD(1) * t857;
t711 = sin(qJ(3));
t716 = cos(qJ(2));
t828 = t711 * t716;
t627 = -qJD(1) * t828 - t712 * t792;
t710 = sin(qJ(4));
t715 = cos(qJ(4));
t810 = qJD(1) * t712;
t875 = -t711 * t810 + t716 * t792;
t584 = t627 * t710 + t715 * t875;
t707 = sin(pkin(11));
t708 = cos(pkin(11));
t756 = t715 * t627 - t710 * t875;
t864 = t584 * t708 + t707 * t756;
t882 = -t864 * t714 + t805;
t703 = qJD(2) + qJD(3);
t881 = t703 * t875;
t787 = qJDD(1) * t857;
t801 = qJDD(1) * t716;
t561 = t711 * t801 + t712 * t787 + t881;
t640 = t712 * t857 + t828;
t595 = t703 * t640;
t802 = qJDD(1) * t712;
t762 = t711 * t802 - t716 * t787;
t562 = qJD(1) * t595 + t762;
t807 = qJD(4) * t715;
t808 = qJD(4) * t710;
t512 = t715 * t561 - t710 * t562 + t627 * t808 + t807 * t875;
t734 = qJD(4) * t756 - t561 * t710 - t715 * t562;
t484 = t512 * t708 + t707 * t734;
t702 = qJDD(2) + qJDD(3);
t695 = qJDD(4) + t702;
t696 = qJD(4) + t703;
t709 = sin(qJ(6));
t797 = t714 * t484 + t709 * t695 + t696 * t805;
t806 = qJD(6) * t709;
t863 = t584 * t707 - t708 * t756;
t466 = -t806 * t863 + t797;
t464 = t466 * t714;
t525 = t696 * t709 + t714 * t863;
t666 = t714 * t695;
t467 = qJD(6) * t525 + t484 * t709 - t666;
t523 = -t714 * t696 + t709 * t863;
t880 = -t709 * t467 - t523 * t882 + t464;
t463 = t466 * t709;
t483 = -t512 * t707 + t708 * t734;
t480 = qJDD(6) - t483;
t477 = t709 * t480;
t532 = qJD(6) - t864;
t840 = t525 * t863;
t879 = t695 * MDP(22) + t584 * MDP(18) * t756 + (-t584 ^ 2 + t756 ^ 2) * MDP(19) + (-t584 * t696 + t512) * MDP(20) + (-t696 * t756 + t734) * MDP(21) + (t525 * t882 + t463) * MDP(27) + (t477 - t840) * MDP(29) + (MDP(29) * t882 - t863 * MDP(31)) * t532;
t577 = t756 * qJ(5);
t620 = t627 * pkin(9);
t858 = pkin(7) + pkin(8);
t657 = t858 * t716;
t646 = qJD(1) * t657;
t628 = t711 * t646;
t656 = t858 * t712;
t644 = qJD(1) * t656;
t847 = qJD(2) * pkin(2);
t634 = -t644 + t847;
t779 = t857 * t634 - t628;
t559 = t620 + t779;
t548 = pkin(3) * t703 + t559;
t632 = t857 * t646;
t748 = -t711 * t634 - t632;
t850 = t875 * pkin(9);
t560 = -t748 + t850;
t550 = t710 * t560;
t785 = t715 * t548 - t550;
t506 = t785 + t577;
t500 = pkin(4) * t696 + t506;
t552 = t715 * t560;
t758 = -t548 * t710 - t552;
t846 = qJ(5) * t584;
t507 = -t758 + t846;
t842 = t507 * t707;
t473 = t500 * t708 - t842;
t471 = -pkin(5) * t696 - t473;
t845 = t471 * t864;
t841 = t523 * t863;
t778 = t644 * t711 - t632;
t564 = t778 - t850;
t816 = -t857 * t644 - t628;
t565 = t620 + t816;
t691 = pkin(2) * t857 + pkin(3);
t830 = t710 * t711;
t877 = -t691 * t807 - (-t711 * t808 + (t715 * t857 - t830) * qJD(3)) * pkin(2) + t710 * t564 + t715 * t565;
t829 = t711 * t715;
t876 = -t691 * t808 + (-t711 * t807 + (-t710 * t857 - t829) * qJD(3)) * pkin(2) - t715 * t564 + t565 * t710;
t706 = qJ(2) + qJ(3);
t700 = qJ(4) + t706;
t688 = sin(t700);
t689 = cos(t700);
t713 = sin(qJ(1));
t717 = cos(qJ(1));
t764 = g(1) * t717 + g(2) * t713;
t874 = -g(3) * t689 + t688 * t764;
t502 = t708 * t507;
t474 = t707 * t500 + t502;
t472 = pkin(10) * t696 + t474;
t701 = t716 * pkin(2);
t849 = pkin(1) + t701;
t655 = t849 * qJD(1);
t598 = -pkin(3) * t875 - t655;
t546 = -pkin(4) * t584 + qJD(5) + t598;
t489 = -pkin(5) * t864 - pkin(10) * t863 + t546;
t459 = -t472 * t709 + t489 * t714;
t685 = pkin(11) + t700;
t672 = sin(t685);
t826 = t714 * t717;
t827 = t713 * t714;
t873 = -t459 * t863 + t471 * t806 + (g(1) * t826 + g(2) * t827) * t672;
t460 = t472 * t714 + t489 * t709;
t803 = qJD(1) * qJD(2);
t789 = t716 * t803;
t596 = qJDD(2) * pkin(2) + t858 * (-t789 - t802);
t790 = t712 * t803;
t597 = t858 * (-t790 + t801);
t732 = qJD(3) * t748 + t857 * t596 - t711 * t597;
t504 = t702 * pkin(3) - t561 * pkin(9) + t732;
t791 = qJD(3) * t857;
t809 = qJD(3) * t711;
t730 = t711 * t596 + t597 * t857 + t634 * t791 - t646 * t809;
t510 = -t562 * pkin(9) + t730;
t736 = qJD(4) * t758 + t715 * t504 - t710 * t510;
t456 = pkin(4) * t695 - qJ(5) * t512 + qJD(5) * t756 + t736;
t860 = t715 * (qJD(4) * t548 + t510) + t710 * t504 - t560 * t808;
t458 = qJ(5) * t734 + qJD(5) * t584 + t860;
t446 = t456 * t708 - t458 * t707;
t444 = -pkin(5) * t695 - t446;
t673 = cos(t685);
t867 = g(3) * t673 + t444;
t872 = t460 * t863 + t471 * t805 + t709 * t867;
t870 = pkin(5) * t863 - pkin(10) * t864;
t869 = pkin(4) * t756;
t866 = -t577 - t877;
t865 = t846 + t876;
t726 = g(3) * t688 - t598 * t584 + t689 * t764 - t860;
t722 = t598 * t756 + t736 + t874;
t681 = pkin(4) * t707 + pkin(10);
t861 = t532 * (qJD(6) * t681 - t869 + t870);
t777 = -t857 * t656 - t657 * t711;
t578 = -pkin(9) * t640 + t777;
t747 = -t711 * t712 + t716 * t857;
t815 = -t711 * t656 + t857 * t657;
t579 = pkin(9) * t747 + t815;
t819 = t710 * t578 + t715 * t579;
t593 = t640 * t715 + t710 * t747;
t594 = t703 * t747;
t520 = qJD(4) * t593 + t594 * t710 + t715 * t595;
t592 = t640 * t710 - t715 * t747;
t794 = qJD(2) * t858;
t645 = t712 * t794;
t647 = t716 * t794;
t742 = -t857 * t645 - t711 * t647 - t656 * t791 - t657 * t809;
t529 = -pkin(9) * t595 + t742;
t731 = -qJD(3) * t815 + t711 * t645 - t857 * t647;
t530 = -t594 * pkin(9) + t731;
t741 = t715 * t529 + t710 * t530 + t578 * t807 - t579 * t808;
t468 = -qJ(5) * t520 - qJD(5) * t592 + t741;
t519 = -qJD(4) * t592 + t594 * t715 - t595 * t710;
t735 = -qJD(4) * t819 - t529 * t710 + t715 * t530;
t724 = -qJ(5) * t519 - qJD(5) * t593 + t735;
t450 = t708 * t468 + t707 * t724;
t517 = -qJ(5) * t592 + t819;
t781 = t715 * t578 - t579 * t710;
t749 = -qJ(5) * t593 + t781;
t488 = t708 * t517 + t707 * t749;
t493 = t519 * t708 - t520 * t707;
t541 = t708 * t592 + t593 * t707;
t542 = -t592 * t707 + t593 * t708;
t603 = -pkin(3) * t747 - t849;
t745 = pkin(4) * t592 + t603;
t495 = pkin(5) * t541 - pkin(10) * t542 + t745;
t447 = t707 * t456 + t708 * t458;
t445 = pkin(10) * t695 + t447;
t774 = qJD(6) * t489 + t445;
t859 = t444 * t542 + t471 * t493 - t488 * t480 - t532 * (qJD(6) * t495 + t450) - t541 * t774;
t856 = pkin(3) * t627;
t853 = g(3) * t672;
t848 = pkin(3) * qJD(4);
t844 = t471 * t542;
t843 = t495 * t480;
t839 = t525 * t709;
t698 = cos(t706);
t836 = t698 * t713;
t835 = t698 * t717;
t834 = t707 * t710;
t833 = t708 * t710;
t832 = t709 * t713;
t831 = t709 * t717;
t478 = t714 * t480;
t823 = t707 * t866 - t708 * t865;
t822 = t707 * t865 + t708 * t866;
t821 = t715 * t559 - t550;
t767 = -pkin(2) * t830 + t715 * t691;
t617 = pkin(4) + t767;
t622 = pkin(2) * t829 + t691 * t710;
t573 = t707 * t617 + t708 * t622;
t511 = t577 + t821;
t784 = -t559 * t710 - t552;
t751 = t784 - t846;
t818 = -t511 * t707 + t708 * t751 + (t707 * t715 + t833) * t848;
t817 = -t708 * t511 - t707 * t751 + (t708 * t715 - t834) * t848;
t814 = pkin(3) * t698 + pkin(4) * t689;
t690 = pkin(3) * t715 + pkin(4);
t619 = pkin(3) * t833 + t707 * t690;
t704 = t712 ^ 2;
t813 = -t716 ^ 2 + t704;
t694 = t712 * t847;
t796 = t701 + t814;
t586 = pkin(3) * t595 + t694;
t776 = t532 * t709;
t621 = pkin(2) * t790 - qJDD(1) * t849;
t545 = pkin(3) * t562 + t621;
t729 = -pkin(4) * t734 + qJDD(5) + t545;
t451 = -pkin(5) * t483 - pkin(10) * t484 + t729;
t772 = qJD(6) * t472 - t451;
t766 = -t856 - t869;
t491 = t766 + t870;
t568 = pkin(10) + t573;
t693 = pkin(2) * t810;
t771 = qJD(6) * t568 + t491 + t693;
t613 = pkin(10) + t619;
t770 = qJD(6) * t613 + t491;
t697 = sin(t706);
t765 = -pkin(3) * t697 - pkin(4) * t688;
t763 = g(1) * t713 - g(2) * t717;
t761 = -t480 * t568 - t845;
t760 = -t480 * t613 - t845;
t759 = t473 * t864 + t474 * t863;
t572 = t617 * t708 - t622 * t707;
t755 = pkin(4) * t520 + t586;
t754 = t478 + (t709 * t864 - t806) * t532;
t618 = -pkin(3) * t834 + t690 * t708;
t753 = t764 * t672;
t752 = -0.2e1 * pkin(1) * t803 - pkin(7) * qJDD(2);
t746 = t493 * t714 - t542 * t806;
t476 = t506 * t708 - t842;
t739 = t476 * t532 - t480 * t681 - t845;
t718 = qJD(2) ^ 2;
t738 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t718 + t763;
t719 = qJD(1) ^ 2;
t737 = pkin(1) * t719 - pkin(7) * qJDD(1) + t764;
t733 = -t714 * t867 + t873;
t727 = -t709 * t753 + t872;
t723 = g(1) * t835 + g(2) * t836 + g(3) * t697 + t655 * t875 - t730;
t721 = -g(3) * t698 - t655 * t627 + t697 * t764 + t732;
t720 = t627 * t875 * MDP(11) + (-t525 * t776 + t880) * MDP(28) + (-t532 * t776 + t478 + t841) * MDP(30) + (t561 - t881) * MDP(13) + (-t762 + (-qJD(1) * t640 - t627) * t703) * MDP(14) + (t627 ^ 2 - t875 ^ 2) * MDP(12) + t702 * MDP(15) + t879;
t699 = -qJ(5) - pkin(9) - t858;
t682 = -pkin(4) * t708 - pkin(5);
t612 = -pkin(5) - t618;
t608 = pkin(1) + t796;
t607 = t673 * t826 + t832;
t606 = -t673 * t831 + t827;
t605 = -t673 * t827 + t831;
t604 = t673 * t832 + t826;
t599 = t693 - t856;
t567 = -pkin(5) - t572;
t492 = t519 * t707 + t708 * t520;
t487 = t517 * t707 - t708 * t749;
t475 = t506 * t707 + t502;
t461 = pkin(5) * t492 - pkin(10) * t493 + t755;
t449 = t468 * t707 - t708 * t724;
t448 = t714 * t451;
t1 = [t763 * MDP(2) + t764 * MDP(3) + (t712 * t752 + t716 * t738) * MDP(9) + (-t712 * t738 + t716 * t752) * MDP(10) + (t598 * t520 + t545 * t592 - t584 * t586 - t603 * t734 + t689 * t763 + t695 * t781 + t696 * t735) * MDP(23) + (-t512 * t592 + t519 * t584 + t520 * t756 + t593 * t734) * MDP(19) + (-t446 * t542 - t447 * t541 + t449 * t863 + t450 * t864 - t473 * t493 - t474 * t492 + t483 * t488 + t484 * t487 - t764) * MDP(25) + (-t561 * t849 - t655 * t594 + t621 * t640 - t627 * t694 - t697 * t763 - t702 * t815 - t703 * t742) * MDP(17) + (t480 * t541 + t492 * t532) * MDP(31) + (qJDD(2) * t716 - t712 * t718) * MDP(7) + (qJDD(2) * t712 + t716 * t718) * MDP(6) + (t512 * t593 - t519 * t756) * MDP(18) + (t603 * t512 + t598 * t519 + t545 * t593 - t586 * t756 - t688 * t763 - t695 * t819 - t696 * t741) * MDP(24) + (-t595 * t703 + t702 * t747) * MDP(14) + (t561 * t640 - t594 * t627) * MDP(11) + (t594 * t703 + t640 * t702) * MDP(13) + (-t520 * t696 - t592 * t695) * MDP(21) + (t519 * t696 + t593 * t695) * MDP(20) + (-g(1) * t604 - g(2) * t606 + t449 * t525 - t460 * t492 + t487 * t466 + (-(-qJD(6) * t488 + t461) * t532 - t843 + t772 * t541 - qJD(6) * t844) * t709 + t859 * t714) * MDP(33) + (-g(1) * t605 - g(2) * t607 + t448 * t541 + t449 * t523 + t459 * t492 + t487 * t467 + (t461 * t532 + t843 + (-t472 * t541 - t488 * t532 + t844) * qJD(6)) * t714 + t859 * t709) * MDP(32) + (g(1) * t836 - g(2) * t835 - t562 * t849 - t655 * t595 - t621 * t747 - t694 * t875 + t702 * t777 + t703 * t731) * MDP(16) + (t561 * t747 - t562 * t640 + t594 * t875 + t595 * t627) * MDP(12) + (qJDD(1) * t704 + 0.2e1 * t712 * t789) * MDP(4) + 0.2e1 * (t712 * t801 - t803 * t813) * MDP(5) + (t466 * t541 + t478 * t542 + t492 * t525 + t532 * t746) * MDP(29) + (-t542 * t477 - t467 * t541 - t492 * t523 + (-t493 * t709 - t542 * t805) * t532) * MDP(30) + (t464 * t542 + t525 * t746) * MDP(27) + ((-t523 * t714 - t839) * t493 + (-t463 - t467 * t714 + (t523 * t709 - t525 * t714) * qJD(6)) * t542) * MDP(28) + qJDD(1) * MDP(1) + (t447 * t488 + t474 * t450 - t446 * t487 - t473 * t449 + t729 * t745 + t546 * t755 - g(1) * (-t608 * t713 - t699 * t717) - g(2) * (t608 * t717 - t699 * t713)) * MDP(26); (g(3) * t712 + t716 * t737) * MDP(10) + t720 + (t483 * t573 - t484 * t572 + t822 * t864 + t823 * t863 + t759) * MDP(25) + MDP(6) * t802 + (t567 * t466 + t761 * t714 + t823 * t525 + (t709 * t771 - t714 * t822) * t532 + t727) * MDP(33) + (-g(3) * t716 + t712 * t737) * MDP(9) + (t447 * t573 + t446 * t572 - t546 * (t693 + t766) - g(3) * t796 - t764 * (-pkin(2) * t712 + t765) + t822 * t474 - t823 * t473) * MDP(26) + (t816 * t703 + (t627 * t810 - t711 * t702 - t703 * t791) * pkin(2) + t723) * MDP(17) + (-t778 * t703 + (t702 * t857 - t703 * t809 + t810 * t875) * pkin(2) + t721) * MDP(16) + (t567 * t467 + t761 * t709 + t823 * t523 + (-t709 * t822 - t714 * t771) * t532 + t733) * MDP(32) + MDP(7) * t801 + (t599 * t584 + t767 * t695 + t696 * t876 + t722) * MDP(23) + (t599 * t756 - t622 * t695 + t696 * t877 + t726) * MDP(24) + qJDD(2) * MDP(8) + (-t712 * t716 * MDP(4) + MDP(5) * t813) * t719; t720 + (t612 * t467 + t760 * t709 + t818 * t523 + (-t709 * t817 - t714 * t770) * t532 + t733) * MDP(32) + (-g(3) * t814 + t446 * t618 + t447 * t619 - t473 * t818 + t474 * t817 - t546 * t766 - t764 * t765) * MDP(26) + (t821 * t696 + (-t627 * t756 - t695 * t710 - t696 * t807) * pkin(3) + t726) * MDP(24) + (-t784 * t696 + (-t584 * t627 + t695 * t715 - t696 * t808) * pkin(3) + t722) * MDP(23) + (t612 * t466 + t760 * t714 + t818 * t525 + (t709 * t770 - t714 * t817) * t532 + t727) * MDP(33) + (t703 * t779 + t723) * MDP(17) + (-t703 * t748 + t721) * MDP(16) + (t483 * t619 - t484 * t618 + t817 * t864 + t818 * t863 + t759) * MDP(25); (-t696 * t758 + t722) * MDP(23) + (t696 * t785 + t726) * MDP(24) + ((t483 * t707 - t484 * t708) * pkin(4) + (t473 - t476) * t864 + (t474 - t475) * t863) * MDP(25) + (t473 * t475 - t474 * t476 + (t446 * t708 + t447 * t707 + t546 * t756 + t874) * pkin(4)) * MDP(26) + (-t532 * t839 + t880) * MDP(28) + (t754 + t841) * MDP(30) + (t682 * t467 - t475 * t523 + t739 * t709 + (-t867 - t861) * t714 + t873) * MDP(32) + (t682 * t466 - t475 * t525 + t739 * t714 + (-t753 + t861) * t709 + t872) * MDP(33) + t879; (-t863 ^ 2 - t864 ^ 2) * MDP(25) + (t473 * t863 - t474 * t864 + t729 - t763) * MDP(26) + (t754 - t841) * MDP(32) + (-t532 ^ 2 * t714 - t477 - t840) * MDP(33); t525 * t523 * MDP(27) + (-t523 ^ 2 + t525 ^ 2) * MDP(28) + (t523 * t532 + t797) * MDP(29) + (t525 * t532 + t666) * MDP(30) + t480 * MDP(31) + (-g(1) * t606 + g(2) * t604 + t460 * t532 - t471 * t525 + t448) * MDP(32) + (g(1) * t607 - g(2) * t605 + t459 * t532 + t471 * t523) * MDP(33) + ((-t445 + t853) * MDP(33) + (-MDP(30) * t863 - MDP(32) * t472 - MDP(33) * t489) * qJD(6)) * t714 + (-qJD(6) * t863 * MDP(29) + (-qJD(6) * t696 - t484) * MDP(30) + (-t774 + t853) * MDP(32) + t772 * MDP(33)) * t709;];
tau  = t1;
