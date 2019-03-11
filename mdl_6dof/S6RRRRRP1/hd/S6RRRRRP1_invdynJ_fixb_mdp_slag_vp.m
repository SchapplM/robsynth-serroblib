% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:58:34
% EndTime: 2019-03-10 00:58:49
% DurationCPUTime: 11.19s
% Computational Cost: add. (12628->577), mult. (29355->714), div. (0->0), fcn. (22028->14), ass. (0->273)
t694 = cos(qJ(5));
t793 = qJD(5) * t694;
t691 = sin(qJ(3));
t692 = sin(qJ(2));
t800 = qJD(1) * t692;
t773 = t691 * t800;
t696 = cos(qJ(3));
t697 = cos(qJ(2));
t799 = qJD(1) * t697;
t775 = t696 * t799;
t604 = -t773 + t775;
t605 = -t691 * t799 - t696 * t800;
t690 = sin(qJ(4));
t695 = cos(qJ(4));
t868 = t695 * t604 + t605 * t690;
t879 = t868 * t694;
t885 = t793 - t879;
t689 = sin(qJ(5));
t683 = qJD(2) + qJD(3);
t786 = qJDD(1) * t697;
t788 = qJD(1) * qJD(2);
t771 = t697 * t788;
t787 = qJDD(1) * t692;
t875 = t771 + t787;
t551 = qJD(3) * t775 - t683 * t773 + t691 * t786 + t696 * t875;
t617 = t691 * t697 + t692 * t696;
t581 = t683 * t617;
t740 = t691 * t787 - t696 * t786;
t552 = qJD(1) * t581 + t740;
t735 = t695 * t551 - t690 * t552;
t708 = qJD(4) * t868 + t735;
t732 = t604 * t690 - t695 * t605;
t767 = qJD(4) + t683;
t739 = t694 * t767;
t682 = qJDD(2) + qJDD(3);
t750 = qJDD(4) + t682;
t794 = qJD(5) * t689;
t482 = -qJD(5) * t739 - t689 * t750 - t694 * t708 + t732 * t794;
t480 = t482 * t689;
t555 = t689 * t767 + t694 * t732;
t704 = t689 * t708 - t694 * t750;
t483 = qJD(5) * t555 + t704;
t765 = t690 * t551 + t695 * t552;
t496 = t732 * qJD(4) + t765;
t495 = qJDD(5) + t496;
t491 = t689 * t495;
t492 = t694 * t495;
t553 = t689 * t732 - t739;
t789 = qJD(5) - t868;
t882 = t789 * t689;
t884 = t750 * MDP(22) - t765 * MDP(21) - t868 ^ 2 * MDP(19) + (-t868 * t683 + t735) * MDP(20) + (-MDP(18) * t868 + MDP(19) * t732 + MDP(21) * t683 - MDP(29) * t789) * t732 + (t555 * t885 - t480) * MDP(25) + (-t555 * t732 + t789 * t885 + t491) * MDP(27) + (-t482 * t694 - t689 * t483 - t553 * t885 - t555 * t882) * MDP(26) + (t553 * t732 - t789 * t882 + t492) * MDP(28);
t687 = qJ(2) + qJ(3);
t680 = qJ(4) + t687;
t666 = sin(t680);
t693 = sin(qJ(1));
t698 = cos(qJ(1));
t746 = g(1) * t698 + g(2) * t693;
t883 = t746 * t666;
t853 = pkin(7) + pkin(8);
t640 = t853 * t697;
t627 = qJD(1) * t640;
t610 = t696 * t627;
t639 = t853 * t692;
t625 = qJD(1) * t639;
t839 = qJD(2) * pkin(2);
t612 = -t625 + t839;
t731 = -t612 * t691 - t610;
t847 = pkin(9) * t604;
t549 = -t731 + t847;
t545 = t690 * t549;
t600 = t605 * pkin(9);
t606 = t691 * t627;
t763 = t696 * t612 - t606;
t548 = t600 + t763;
t507 = t548 * t695 - t545;
t795 = qJD(4) * t695;
t878 = -pkin(3) * t795 + t507;
t762 = t625 * t691 - t610;
t558 = t762 - t847;
t805 = -t696 * t625 - t606;
t559 = t600 + t805;
t671 = pkin(2) * t696 + pkin(3);
t796 = qJD(4) * t690;
t825 = t690 * t691;
t869 = t558 * t690 + t559 * t695 - t671 * t795 - (-t691 * t796 + (t695 * t696 - t825) * qJD(3)) * pkin(2);
t540 = pkin(3) * t683 + t548;
t502 = t695 * t540 - t545;
t499 = -pkin(4) * t767 - t502;
t837 = t499 * t868;
t678 = cos(t687);
t667 = cos(t680);
t848 = pkin(5) * t694;
t669 = pkin(4) + t848;
t688 = -qJ(6) - pkin(10);
t760 = -t666 * t688 + t667 * t669;
t743 = pkin(3) * t678 + t760;
t824 = t691 * t695;
t876 = (t691 * t795 + (t690 * t696 + t824) * qJD(3)) * pkin(2) + t558 * t695;
t532 = pkin(4) * t732 - pkin(10) * t868;
t681 = t697 * pkin(2);
t840 = pkin(1) + t681;
t638 = t840 * qJD(1);
t584 = -pkin(3) * t604 - t638;
t658 = g(3) * t666;
t778 = t667 * t746 + t658;
t582 = qJDD(2) * pkin(2) - t853 * t875;
t772 = t692 * t788;
t583 = t853 * (-t772 + t786);
t713 = qJD(3) * t731 + t696 * t582 - t691 * t583;
t488 = pkin(3) * t682 - pkin(9) * t551 + t713;
t798 = qJD(3) * t691;
t856 = (qJD(3) * t612 + t583) * t696 + t691 * t582 - t627 * t798;
t494 = -pkin(9) * t552 + t856;
t857 = -t695 * (qJD(4) * t540 + t494) - t690 * t488 + t549 * t796;
t706 = -t584 * t868 + t778 + t857;
t679 = t694 * qJ(6);
t744 = pkin(5) * t732 - t868 * t679;
t745 = g(1) * t693 - g(2) * t698;
t870 = t745 * t666;
t810 = -t559 * t690 + t671 * t796 + t876;
t804 = -t691 * t639 + t696 * t640;
t676 = t694 * qJD(6);
t834 = t868 * t689;
t867 = qJ(6) * t834 + t676;
t677 = sin(t687);
t730 = t666 * t669 + t667 * t688;
t866 = pkin(3) * t677 + t730;
t852 = pkin(3) * t605;
t520 = t532 - t852;
t865 = t689 * t520 + t694 * t878;
t673 = pkin(2) * t800;
t515 = t520 + t673;
t864 = t689 * t515 + t694 * t869;
t822 = t694 * t698;
t827 = t689 * t693;
t591 = t667 * t827 + t822;
t823 = t693 * t694;
t826 = t689 * t698;
t593 = -t667 * t826 + t823;
t863 = -g(1) * t593 + g(2) * t591;
t842 = g(3) * t667;
t862 = t842 - t883;
t546 = t695 * t549;
t503 = t690 * t540 + t546;
t500 = pkin(10) * t767 + t503;
t512 = -pkin(4) * t868 - pkin(10) * t732 + t584;
t476 = -t500 * t689 + t694 * t512;
t861 = -t476 * t732 + t499 * t794 + t694 * t883;
t748 = -t695 * t488 + t690 * t494 + t540 * t796 + t549 * t795;
t451 = -pkin(4) * t750 + t748;
t477 = t500 * t694 + t512 * t689;
t498 = t499 * t793;
t841 = g(3) * t689;
t860 = t451 * t689 + t477 * t732 + t667 * t841 + t498;
t714 = -t584 * t732 - t748 - t862;
t854 = t555 ^ 2;
t850 = pkin(3) * t695;
t849 = pkin(5) * t689;
t463 = -qJ(6) * t555 + t476;
t460 = pkin(5) * t789 + t463;
t838 = t460 * t694;
t616 = t691 * t692 - t696 * t697;
t577 = -t616 * t690 + t617 * t695;
t831 = t577 * t689;
t830 = t577 * t694;
t761 = -t696 * t639 - t640 * t691;
t562 = -pkin(9) * t617 + t761;
t563 = -pkin(9) * t616 + t804;
t529 = t562 * t690 + t563 * t695;
t526 = t694 * t529;
t803 = pkin(2) * t824 + t690 * t671;
t597 = pkin(10) + t803;
t821 = -qJ(6) - t597;
t668 = pkin(3) * t690 + pkin(10);
t820 = -qJ(6) - t668;
t819 = -t463 + t460;
t817 = t694 * t502 + t689 * t532;
t576 = t695 * t616 + t617 * t690;
t588 = pkin(3) * t616 - t840;
t527 = pkin(4) * t576 - pkin(10) * t577 + t588;
t814 = t689 * t527 + t526;
t758 = qJD(5) * t821;
t813 = t689 * t758 - t864 + t867;
t511 = t694 * t515;
t812 = t694 * t758 - t511 - t744 + (-qJD(6) + t869) * t689;
t757 = qJD(5) * t820;
t809 = t689 * t757 - t865 + t867;
t519 = t694 * t520;
t808 = t694 * t757 - t519 - t744 + (-qJD(6) + t878) * t689;
t768 = qJD(5) * t688;
t807 = t676 - t817 + (qJ(6) * t868 + t768) * t689;
t531 = t694 * t532;
t806 = t694 * t768 - t531 + (-qJD(6) + t502) * t689 - t744;
t802 = (t794 - t834) * pkin(5);
t685 = t692 ^ 2;
t801 = -t697 ^ 2 + t685;
t797 = qJD(3) * t696;
t675 = t692 * t839;
t784 = qJD(5) * pkin(10) * t789;
t776 = qJD(2) * t853;
t626 = t692 * t776;
t628 = t697 * t776;
t720 = -t696 * t626 - t691 * t628 - t639 * t797 - t640 * t798;
t524 = -pkin(9) * t581 + t720;
t580 = t683 * t616;
t712 = -qJD(3) * t804 + t626 * t691 - t696 * t628;
t525 = pkin(9) * t580 + t712;
t734 = t562 * t695 - t563 * t690;
t467 = qJD(4) * t734 + t524 * t695 + t525 * t690;
t516 = -qJD(4) * t576 - t580 * t695 - t581 * t690;
t517 = qJD(4) * t577 - t580 * t690 + t695 * t581;
t571 = pkin(3) * t581 + t675;
t474 = pkin(4) * t517 - pkin(10) * t516 + t571;
t781 = t694 * t467 + t689 * t474 + t527 * t793;
t774 = t577 * t793;
t770 = pkin(9) + t853 + t849;
t769 = -t451 - t842;
t759 = -pkin(2) * t825 + t671 * t695;
t450 = pkin(10) * t750 - t857;
t752 = -qJD(5) * t512 - t450;
t601 = pkin(2) * t772 - qJDD(1) * t840;
t537 = t552 * pkin(3) + t601;
t456 = t496 * pkin(4) - pkin(10) * t708 + t537;
t719 = t694 * t450 + t689 * t456 - t500 * t794 + t512 * t793;
t443 = -qJ(6) * t483 - qJD(6) * t553 + t719;
t749 = t443 * t694 - t778;
t596 = -pkin(4) - t759;
t506 = t548 * t690 + t546;
t747 = pkin(3) * t796 - t506;
t453 = t694 * t456;
t742 = -t500 * t793 + t453;
t741 = -pkin(10) * t495 - t837;
t464 = -qJ(6) * t553 + t477;
t738 = -t464 * t689 - t838;
t737 = -t495 * t597 - t837;
t736 = -t495 * t668 - t837;
t729 = -qJ(6) * t516 - qJD(6) * t577;
t727 = -0.2e1 * pkin(1) * t788 - pkin(7) * qJDD(2);
t726 = -t840 - t743;
t724 = t516 * t689 + t774;
t723 = t516 * t694 - t577 * t794;
t699 = qJD(2) ^ 2;
t716 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t699 + t745;
t700 = qJD(1) ^ 2;
t715 = pkin(1) * t700 - pkin(7) * qJDD(1) + t746;
t711 = t769 * t694 + t861;
t468 = qJD(4) * t529 + t524 * t690 - t525 * t695;
t707 = -t689 * t883 + t860;
t447 = t483 * pkin(5) + qJDD(6) + t451;
t705 = g(3) * t677 + t638 * t604 + t678 * t746 - t856;
t441 = pkin(5) * t495 + qJ(6) * t482 - qJD(5) * t477 - qJD(6) * t555 - t450 * t689 + t453;
t703 = qJD(5) * t738 - t441 * t689 + t460 * t879 + t464 * t834 + t749;
t702 = -g(3) * t678 - t638 * t605 + t677 * t746 + t713;
t701 = t605 * t604 * MDP(11) + (-t604 * t683 + t551) * MDP(13) + (-t740 + (-qJD(1) * t617 - t605) * t683) * MDP(14) + (-t604 ^ 2 + t605 ^ 2) * MDP(12) + t682 * MDP(15) + t884;
t670 = -pkin(4) - t850;
t637 = pkin(10) * t694 + t679;
t636 = t688 * t689;
t614 = t668 * t694 + t679;
t613 = t820 * t689;
t594 = t667 * t822 + t827;
t592 = -t667 * t823 + t826;
t587 = t673 - t852;
t586 = t597 * t694 + t679;
t585 = t821 * t689;
t550 = t553 ^ 2;
t523 = t694 * t527;
t486 = t553 * pkin(5) + qJD(6) + t499;
t478 = -qJ(6) * t831 + t814;
t473 = t694 * t474;
t471 = pkin(5) * t576 - t529 * t689 - t577 * t679 + t523;
t445 = -qJ(6) * t774 + (-qJD(5) * t529 + t729) * t689 + t781;
t444 = pkin(5) * t517 - t467 * t689 + t473 + t729 * t694 + (-t526 + (qJ(6) * t577 - t527) * t689) * qJD(5);
t1 = [(-t468 * t767 + t588 * t496 + t584 * t517 + t537 * t576 - t571 * t868 + t667 * t745 + t734 * t750) * MDP(23) + (-t577 * t496 + t516 * t868 - t517 * t732 - t576 * t708) * MDP(19) + (-t551 * t616 - t552 * t617 - t580 * t604 + t581 * t605) * MDP(12) + (t551 * t617 + t580 * t605) * MDP(11) + (-t580 * t683 + t617 * t682) * MDP(13) + ((-t529 * t793 + t473) * t789 + t523 * t495 + t742 * t576 + t476 * t517 + t468 * t553 - t734 * t483 + t577 * t498 - g(1) * t592 - g(2) * t594 + ((-qJD(5) * t527 - t467) * t789 - t529 * t495 + t752 * t576 + t451 * t577 + t499 * t516) * t689) * MDP(30) + (-(-t529 * t794 + t781) * t789 - t814 * t495 - t719 * t576 - t477 * t517 + t468 * t555 + t734 * t482 + t451 * t830 - g(1) * t591 - g(2) * t593 + t723 * t499) * MDP(31) + (t495 * t576 + t517 * t789) * MDP(29) + (-t482 * t576 + t492 * t577 + t517 * t555 + t723 * t789) * MDP(27) + (-t483 * t576 - t491 * t577 - t517 * t553 - t724 * t789) * MDP(28) + (qJDD(2) * t692 + t697 * t699) * MDP(6) + (qJDD(2) * t697 - t692 * t699) * MDP(7) + t745 * MDP(2) + t746 * MDP(3) + (t516 * t732 + t577 * t708) * MDP(18) + (t692 * t727 + t697 * t716) * MDP(9) + (-t692 * t716 + t697 * t727) * MDP(10) + (-t552 * t840 - t638 * t581 + t601 * t616 - t604 * t675 + t678 * t745 + t682 * t761 + t683 * t712) * MDP(16) + (-t551 * t840 + t638 * t580 + t601 * t617 - t605 * t675 - t677 * t745 - t682 * t804 - t683 * t720) * MDP(17) + (-t581 * t683 - t616 * t682) * MDP(14) + qJDD(1) * MDP(1) + (-t467 * t767 + t584 * t516 - t529 * t750 + t537 * t577 + t571 * t732 + t588 * t708 - t870) * MDP(24) + (-t444 * t555 - t445 * t553 + t471 * t482 - t478 * t483 + t870 + t738 * t516 + (-t441 * t694 - t443 * t689 + (t460 * t689 - t464 * t694) * qJD(5)) * t577) * MDP(32) + (-t517 * t767 - t576 * t750) * MDP(21) + (t516 * t767 + t577 * t750) * MDP(20) + (qJDD(1) * t685 + 0.2e1 * t692 * t771) * MDP(4) + 0.2e1 * (t692 * t786 - t788 * t801) * MDP(5) + (-t482 * t830 + t555 * t723) * MDP(25) + (t443 * t478 + t464 * t445 + t441 * t471 + t460 * t444 + t447 * (pkin(5) * t831 - t734) + t486 * (pkin(5) * t724 + t468) + (-g(1) * t770 + g(2) * t726) * t698 + (-g(1) * t726 - g(2) * t770) * t693) * MDP(33) + ((-t553 * t694 - t555 * t689) * t516 + (t480 - t483 * t694 + (t553 * t689 - t555 * t694) * qJD(5)) * t577) * MDP(26); (t482 * t585 - t483 * t586 - t553 * t813 - t555 * t812 + t703) * MDP(32) + (t805 * t683 + (t605 * t800 - t682 * t691 - t683 * t797) * pkin(2) + t705) * MDP(17) + (-t762 * t683 + (t604 * t800 + t682 * t696 - t683 * t798) * pkin(2) + t702) * MDP(16) + t701 + (-g(3) * t697 + t692 * t715) * MDP(9) + (g(3) * t692 + t697 * t715) * MDP(10) + (t443 * t586 + t441 * t585 + t447 * (t596 - t848) - g(3) * (t681 + t743) + t813 * t464 + t812 * t460 + t746 * (pkin(2) * t692 + t866) + ((qJD(4) * t671 - t559) * t690 + t802 + t876) * t486) * MDP(33) + (t596 * t483 + t737 * t689 + t810 * t553 + (-t597 * t793 + t689 * t869 - t511) * t789 + t711) * MDP(30) + (t587 * t868 + t750 * t759 - t767 * t810 + t714) * MDP(23) + (-t587 * t732 - t750 * t803 + t767 * t869 + t706) * MDP(24) + MDP(7) * t786 + MDP(6) * t787 + qJDD(2) * MDP(8) + (-t596 * t482 + t737 * t694 + t810 * t555 + (t597 * t794 + t864) * t789 + t707) * MDP(31) + (-MDP(4) * t692 * t697 + MDP(5) * t801) * t700; (t443 * t614 + t441 * t613 + t447 * (-t669 - t850) - g(3) * t743 + (-t546 + (pkin(3) * qJD(4) - t548) * t690 + t802) * t486 + t809 * t464 + t808 * t460 + t746 * t866) * MDP(33) + t701 + (t482 * t613 - t483 * t614 - t553 * t809 - t555 * t808 + t703) * MDP(32) + (t670 * t483 + t736 * t689 + t747 * t553 + (-t668 * t793 + t689 * t878 - t519) * t789 + t711) * MDP(30) + (t507 * t767 + (t605 * t732 - t690 * t750 - t767 * t795) * pkin(3) + t706) * MDP(24) + (-t670 * t482 + t736 * t694 + t747 * t555 + (t668 * t794 + t865) * t789 + t707) * MDP(31) + (-t683 * t731 + t702) * MDP(16) + (t683 * t763 + t705) * MDP(17) + (t506 * t767 + (-t605 * t868 + t695 * t750 - t767 * t796) * pkin(3) + t714) * MDP(23); (t503 * t767 + t714) * MDP(23) + (t502 * t767 + t706) * MDP(24) + (-pkin(4) * t483 - t503 * t553 - t531 * t789 + (t502 * t789 + t741) * t689 + (t769 - t784) * t694 + t861) * MDP(30) + (pkin(4) * t482 + t817 * t789 - t503 * t555 + t741 * t694 + (-t883 + t784) * t689 + t860) * MDP(31) + (t482 * t636 - t483 * t637 - t806 * t555 - t807 * t553 - t789 * t838 + (-t464 * t789 - t441) * t689 + t749) * MDP(32) + (t443 * t637 + t441 * t636 - t447 * t669 - g(3) * t760 + (t789 * t849 - t503) * t486 + t807 * t464 + t806 * t460 + t746 * t730) * MDP(33) + t884; t555 * t553 * MDP(25) + (-t550 + t854) * MDP(26) + (t553 * t789 - t482) * MDP(27) + (-t704 + (-qJD(5) + t789) * t555) * MDP(28) + t495 * MDP(29) + (t477 * t789 - t499 * t555 + (t752 + t658) * t689 + t742 + t863) * MDP(30) + (g(1) * t594 - g(2) * t592 + t476 * t789 + t499 * t553 + t658 * t694 - t719) * MDP(31) + (pkin(5) * t482 - t553 * t819) * MDP(32) + (t819 * t464 + (-t486 * t555 + t666 * t841 + t441 + t863) * pkin(5)) * MDP(33); (-t550 - t854) * MDP(32) + (t460 * t555 + t464 * t553 + t447 + t862) * MDP(33);];
tau  = t1;
