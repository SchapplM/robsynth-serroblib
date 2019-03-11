% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR13_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR13_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR13_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR13_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:25:46
% EndTime: 2019-03-09 04:26:02
% DurationCPUTime: 13.33s
% Computational Cost: add. (11225->692), mult. (35746->921), div. (0->0), fcn. (31387->14), ass. (0->293)
t671 = sin(qJ(3));
t826 = cos(pkin(7));
t832 = cos(qJ(3));
t743 = t826 * t832;
t667 = sin(pkin(6));
t668 = cos(pkin(12));
t813 = t667 * t668;
t713 = t743 * t813;
t700 = qJD(1) * t713;
t666 = sin(pkin(7));
t827 = cos(pkin(6));
t765 = t827 * t666;
t726 = t832 * t765;
t665 = sin(pkin(12));
t799 = qJD(1) * t667;
t778 = t665 * t799;
t572 = -qJD(1) * t726 + t671 * t778 - t700;
t736 = t827 * t826;
t777 = t668 * t799;
t611 = -qJD(1) * t736 + t666 * t777 - qJD(3);
t670 = sin(qJ(5));
t674 = cos(qJ(5));
t526 = -t674 * t572 - t611 * t670;
t769 = t671 * t826;
t744 = t668 * t769;
t745 = t671 * t765;
t586 = t667 * (t665 * t832 + t744) + t745;
t841 = qJD(1) * t586;
t852 = qJD(5) + t841;
t856 = t526 * t852;
t672 = sin(qJ(1));
t675 = cos(qJ(1));
t767 = t675 * t827;
t617 = t672 * t665 - t668 * t767;
t618 = t665 * t767 + t672 * t668;
t779 = t667 * t832;
t753 = t666 * t779;
t532 = t617 * t743 + t618 * t671 + t675 * t753;
t770 = t667 * t826;
t590 = t617 * t666 - t675 * t770;
t509 = t532 * t670 + t590 * t674;
t811 = t667 * t675;
t814 = t666 * t671;
t533 = -t617 * t769 + t618 * t832 - t811 * t814;
t669 = sin(qJ(6));
t673 = cos(qJ(6));
t855 = t509 * t669 - t533 * t673;
t854 = t509 * t673 + t533 * t669;
t691 = t667 * (t665 * t743 + t668 * t671);
t600 = qJD(1) * t691;
t798 = qJD(3) * t671;
t776 = t666 * t798;
t853 = t600 - t776;
t616 = t666 * t813 - t736;
t781 = pkin(1) * t827;
t621 = qJ(2) * t813 + t665 * t781;
t690 = (t668 * t770 + t765) * pkin(9);
t582 = t690 + t621;
t658 = t668 * t781;
t815 = t665 * t667;
t687 = t827 * pkin(2) + (-pkin(9) * t826 - qJ(2)) * t815;
t587 = t658 + t687;
t604 = (-pkin(9) * t665 * t666 - pkin(2) * t668 - pkin(1)) * t667;
t725 = qJD(3) * t743;
t746 = t665 * t770;
t775 = qJD(3) * t832;
t748 = t666 * t775;
t752 = t668 * t779;
t679 = -t671 * (qJD(2) * t746 + qJD(3) * t582) + qJD(2) * t752 + t587 * t725 + t604 * t748;
t486 = t616 * qJD(4) - t679;
t785 = t671 * t815;
t849 = (-t726 + t785) * qJD(3);
t848 = -t713 - t726;
t525 = qJD(6) + t526;
t847 = t532 * t674 - t590 * t670;
t528 = t572 * t670 - t611 * t674;
t819 = t528 * t669;
t500 = -t673 * t852 + t819;
t846 = t500 * t852;
t663 = t667 ^ 2;
t845 = t663 * (t665 ^ 2 + t668 ^ 2);
t757 = t674 * t852;
t780 = t666 * t832;
t686 = t671 * t582 - t587 * t743 - t604 * t780;
t833 = pkin(3) + pkin(10);
t469 = t586 * pkin(4) + t616 * t833 + t686;
t585 = t785 + t848;
t529 = -t587 * t666 + t826 * t604;
t824 = qJ(4) * t586;
t718 = t529 - t824;
t478 = t585 * t833 + t718;
t844 = t670 * t469 + t674 * t478;
t696 = -t670 * t826 - t674 * t780;
t751 = t666 * t778;
t843 = -qJD(5) * t696 + t670 * t853 + t674 * t751;
t623 = -t670 * t780 + t674 * t826;
t842 = qJD(5) * t623 - t670 * t751 + t674 * t853;
t601 = (-t671 * t746 + t752) * qJD(1);
t724 = t748 - t601;
t762 = qJD(1) * t827;
t749 = pkin(1) * t762;
t614 = qJ(2) * t777 + t665 * t749;
t557 = qJD(1) * t690 + t614;
t650 = t668 * t749;
t567 = qJD(1) * t687 + t650;
t596 = qJD(1) * t604 + qJD(2);
t495 = t671 * t557 - t567 * t743 - t596 * t780;
t791 = -qJD(4) - t495;
t548 = t567 * t769;
t496 = t832 * t557 + t596 * t814 + t548;
t481 = -pkin(4) * t572 + t496;
t606 = t611 * qJ(4);
t466 = t481 - t606;
t787 = qJDD(1) * t667;
t773 = t665 * t787;
t522 = qJD(1) * t849 - qJD(3) * t700 - qJDD(1) * t745 - t744 * t787 - t832 * t773;
t520 = -qJDD(5) + t522;
t840 = -t466 * t852 - t520 * t833;
t790 = pkin(4) * t841 - t791;
t463 = t611 * t833 + t790;
t524 = -t567 * t666 + t826 * t596;
t719 = -qJ(4) * t841 + t524;
t465 = t572 * t833 + t719;
t440 = t463 * t670 + t465 * t674;
t786 = qJDD(1) * t668;
t772 = t667 * t786;
t610 = -qJDD(1) * t736 + t666 * t772 - qJDD(3);
t747 = qJDD(1) * t781;
t788 = qJD(1) * qJD(2);
t774 = t667 * t788;
t598 = qJ(2) * t772 + t665 * t747 + t668 * t774;
t542 = qJDD(1) * t690 + t598;
t648 = t668 * t747;
t543 = qJDD(1) * t687 - t665 * t774 + t648;
t592 = qJDD(1) * t604 + qJDD(2);
t722 = qJD(3) * t548 + t671 * t542 - t543 * t743 + t557 * t775 - t592 * t780 + t596 * t776;
t702 = qJDD(4) + t722;
t443 = -pkin(4) * t522 + t610 * t833 + t702;
t575 = t586 * qJD(3);
t704 = qJDD(1) * t848 + t671 * t773;
t523 = qJD(1) * t575 + t704;
t512 = -t543 * t666 + t826 * t592;
t693 = qJ(4) * t522 - qJD(4) * t841 + t512;
t451 = t523 * t833 + t693;
t760 = -t674 * t443 + t451 * t670;
t835 = t440 * qJD(5) + t760;
t430 = pkin(5) * t520 + t835;
t768 = t672 * t827;
t619 = -t665 * t768 + t675 * t668;
t701 = t675 * t665 + t668 * t768;
t689 = t701 * t826;
t536 = t619 * t671 - t672 * t753 + t689 * t832;
t591 = t666 * t701 + t672 * t770;
t510 = t536 * t674 - t591 * t670;
t728 = t585 * t674 + t616 * t670;
t707 = g(1) * t510 + g(2) * t847 + g(3) * t728;
t839 = t525 * (pkin(5) * t528 + pkin(11) * t525) + t430 + t707;
t758 = t674 * t523 + t610 * t670;
t475 = qJD(5) * t528 - t758;
t472 = qJDD(6) + t475;
t642 = pkin(5) * t670 - pkin(11) * t674 + qJ(4);
t742 = pkin(5) * t674 + pkin(11) * t670;
t838 = t525 * ((-pkin(4) - t742) * t841 - qJD(5) * t742 + t791) - t642 * t472;
t834 = t841 ^ 2;
t831 = pkin(1) * t663;
t830 = pkin(3) * t585;
t829 = pkin(3) * t610;
t828 = g(1) * t672;
t825 = qJ(4) * t572;
t603 = qJ(4) * t610;
t796 = qJD(5) * t674;
t797 = qJD(5) * t670;
t474 = t670 * t523 + t572 * t796 - t674 * t610 + t611 * t797;
t793 = qJD(6) * t673;
t784 = t673 * t474 - t669 * t520 + t793 * t852;
t794 = qJD(6) * t669;
t448 = -t528 * t794 + t784;
t823 = t448 * t669;
t822 = t496 * t611;
t821 = t500 * t525;
t502 = t528 * t673 + t669 * t852;
t820 = t502 * t525;
t818 = t572 * t841;
t817 = t841 * t670;
t812 = t667 * t672;
t810 = t669 * t472;
t809 = t669 * t833;
t808 = t673 * t472;
t807 = t673 * t833;
t806 = t674 * t448;
t499 = t833 * t841 + t825;
t802 = t670 * t481 + t674 * t499;
t801 = t675 * pkin(1) + qJ(2) * t812;
t795 = qJD(5) * t833;
t792 = t610 * MDP(12);
t565 = t832 * t582;
t764 = t826 * t587;
t782 = t604 * t814 + t671 * t764 + t565;
t771 = -t672 * pkin(1) + qJ(2) * t811;
t677 = qJD(1) ^ 2;
t766 = t677 * t827;
t709 = t670 * t443 + t674 * t451 + t463 * t796 - t465 * t797;
t429 = -pkin(11) * t520 + t709;
t605 = qJD(4) * t611;
t723 = -t832 * t542 - t543 * t769 + t557 * t798 - t567 * t725 - t592 * t814 - t596 * t748;
t454 = t603 + t605 + t723;
t444 = -pkin(4) * t523 - t454;
t434 = pkin(5) * t475 - pkin(11) * t474 + t444;
t761 = -t669 * t429 + t673 * t434;
t759 = t474 * t669 + t673 * t520;
t756 = t852 * t528;
t755 = t852 * t502;
t754 = t525 * t673;
t493 = t616 * qJ(4) - t782;
t635 = t666 * qJD(2) * t815;
t741 = -g(1) * t532 + g(2) * t536;
t740 = g(1) * t675 + g(2) * t672;
t739 = qJD(2) * t762;
t738 = g(2) * t811 - g(3) * t827;
t514 = -t572 * t669 + t673 * t817;
t737 = -t673 * t797 - t514;
t735 = t673 * t429 + t669 * t434;
t438 = pkin(11) * t852 + t440;
t455 = pkin(5) * t526 - pkin(11) * t528 + t466;
t432 = t438 * t673 + t455 * t669;
t734 = t438 * t669 - t455 * t673;
t447 = pkin(11) * t586 + t844;
t483 = -pkin(4) * t585 - t493;
t531 = t585 * t670 - t616 * t674;
t458 = -pkin(5) * t728 - pkin(11) * t531 + t483;
t733 = t447 * t673 + t458 * t669;
t732 = -t447 * t669 + t458 * t673;
t439 = t463 * t674 - t465 * t670;
t731 = t469 * t674 - t478 * t670;
t490 = qJD(2) * t691 + (t565 + (t604 * t666 + t764) * t671) * qJD(3);
t574 = -qJD(3) * t713 + t849;
t473 = -t574 * pkin(4) + t490;
t714 = qJ(4) * t574 - qJD(4) * t586 + t635;
t482 = t575 * t833 + t714;
t729 = t473 * t674 - t482 * t670;
t506 = t531 * t673 + t586 * t669;
t505 = t531 * t669 - t586 * t673;
t727 = (-qJ(2) * t778 + t650) * t665 - t614 * t668;
t721 = -t623 * t669 + t673 * t814;
t720 = t623 * t673 + t669 * t814;
t717 = -t670 * t852 ^ 2 - t674 * t520;
t716 = -t525 * t793 - t810;
t715 = -t525 * t794 + t808;
t708 = t469 * t796 + t670 * t473 - t478 * t797 + t674 * t482;
t706 = g(1) * t536 + g(2) * t532 + g(3) * t585;
t537 = t619 * t832 + (t666 * t812 - t689) * t671;
t705 = g(1) * t537 + g(2) * t533 + g(3) * t586;
t703 = g(1) * t533 - g(2) * t537 + t490 * t611;
t699 = t670 * t520 - t757 * t852;
t697 = t444 - t705;
t437 = -pkin(5) * t852 - t439;
t692 = -pkin(11) * t472 + (t437 + t439) * t525;
t688 = -qJD(6) * t525 * t833 + t705;
t683 = (-pkin(11) * t572 - qJD(6) * t642 + t802) * t525 + t706;
t682 = t706 - t722;
t681 = -t705 - t723;
t680 = -t572 * t611 - t522;
t487 = pkin(3) * t572 + t719;
t678 = t487 * t841 + qJDD(4) - t682;
t464 = -t575 * pkin(4) - t486;
t651 = -pkin(1) * t787 + qJDD(2);
t620 = -qJ(2) * t815 + t658;
t597 = t648 + (-qJ(2) * qJDD(1) - t788) * t815;
t521 = pkin(3) * t841 + t825;
t513 = t673 * t572 + t669 * t817;
t511 = t536 * t670 + t591 * t674;
t504 = qJD(5) * t531 - t575 * t674;
t503 = qJD(5) * t728 + t575 * t670;
t497 = pkin(3) * t575 + t714;
t494 = t616 * pkin(3) + t686;
t492 = t718 + t830;
t489 = t606 - t496;
t488 = pkin(3) * t611 - t791;
t485 = t511 * t673 + t537 * t669;
t484 = -t511 * t669 + t537 * t673;
t461 = -qJD(6) * t505 + t503 * t673 - t574 * t669;
t460 = qJD(6) * t506 + t503 * t669 + t574 * t673;
t457 = pkin(3) * t523 + t693;
t456 = t702 + t829;
t452 = pkin(5) * t572 - t481 * t674 + t499 * t670;
t449 = qJD(6) * t502 + t759;
t446 = -pkin(5) * t586 - t731;
t445 = t504 * pkin(5) - t503 * pkin(11) + t464;
t436 = pkin(5) * t574 + qJD(5) * t844 - t729;
t435 = -pkin(11) * t574 + t708;
t428 = -t432 * qJD(6) + t761;
t427 = -t734 * qJD(6) + t735;
t1 = [(t788 * t845 + (-t597 * t665 + t598 * t668 + (-t620 * t665 + t621 * t668) * qJDD(1) - t740) * t667) * MDP(6) + ((-qJD(6) * t733 - t435 * t669 + t445 * t673) * t525 + t732 * t472 - t428 * t728 - t734 * t504 + t436 * t500 + t446 * t449 + t430 * t505 + t437 * t460 + g(1) * t854 - g(2) * t485) * MDP(31) + (-(qJD(6) * t732 + t435 * t673 + t445 * t669) * t525 - t733 * t472 + t427 * t728 - t432 * t504 + t436 * t502 + t446 * t448 + t430 * t506 + t437 * t461 - g(1) * t855 - g(2) * t484) * MDP(32) + (g(1) * t590 - g(2) * t591 + t454 * t585 + t456 * t586 + t486 * t572 - t488 * t574 + t489 * t575 + t490 * t841 + t493 * t523 - t494 * t522) * MDP(15) + (t512 * t586 - t529 * t522 - t524 * t574 + t610 * t782 + t611 * t679 - t616 * t723 + t635 * t841 + t741) * MDP(14) + (t454 * t616 - t457 * t586 + t486 * t611 + t487 * t574 + t492 * t522 + t493 * t610 - t497 * t841 - t741) * MDP(17) + (t522 * t585 - t523 * t586 + t572 * t574 - t575 * t841) * MDP(9) + (-t522 * t586 - t574 * t841) * MDP(8) + (t729 * t852 - t731 * t520 - t760 * t586 - t439 * t574 + t464 * t526 + t483 * t475 - t444 * t728 + t466 * t504 + g(1) * t509 - g(2) * t511 + (-t440 * t586 - t844 * t852) * qJD(5)) * MDP(24) + (g(1) * t847 - g(2) * t510 + t440 * t574 + t444 * t531 + t464 * t528 + t466 * t503 + t483 * t474 + t844 * t520 - t709 * t586 - t708 * t852) * MDP(25) + (-t475 * t586 - t504 * t852 - t520 * t728 + t526 * t574) * MDP(22) + (t474 * t586 + t503 * t852 - t520 * t531 - t528 * t574) * MDP(21) + (-t520 * t586 - t574 * t852) * MDP(23) + (t474 * t728 - t475 * t531 - t503 * t526 - t504 * t528) * MDP(20) + (-t472 * t728 + t504 * t525) * MDP(30) + (t449 * t728 - t460 * t525 - t472 * t505 - t500 * t504) * MDP(29) + (-t448 * t728 + t461 * t525 + t472 * t506 + t502 * t504) * MDP(28) + (t457 * t492 + t487 * t497 + t454 * t493 + t489 * t486 + t456 * t494 + t488 * t490 - g(1) * (-t618 * pkin(2) - pkin(3) * t533 - pkin(9) * t590 - qJ(4) * t532 + t771) - g(2) * (t619 * pkin(2) + t537 * pkin(3) + pkin(9) * t591 + t536 * qJ(4) + t801)) * MDP(18) + (t522 * t616 + t574 * t611 - t586 * t610) * MDP(10) + (-t448 * t505 - t449 * t506 - t460 * t502 - t461 * t500) * MDP(27) + (t448 * t506 + t461 * t502) * MDP(26) + (t474 * t531 + t503 * t528) * MDP(19) + (t598 * t621 + t597 * t620 - g(1) * t771 - g(2) * t801 + (-t651 * pkin(1) - qJD(2) * t727) * t667) * MDP(7) + (t523 * t616 + t575 * t611 + t585 * t610) * MDP(11) + (-t598 * t827 - g(1) * t617 + g(2) * t701 + (t651 * t665 - t668 * t739) * t667 + (-t621 * t827 - t665 * t831) * qJDD(1)) * MDP(5) + (t597 * t827 + g(1) * t618 - g(2) * t619 + (-t651 * t668 - t665 * t739) * t667 + (t620 * t827 + t668 * t831) * qJDD(1)) * MDP(4) + (-t456 * t616 - t457 * t585 - t487 * t575 - t492 * t523 - t494 * t610 - t497 * t572 - t703) * MDP(16) + qJDD(1) * MDP(1) + (t512 * t585 + t529 * t523 + t524 * t575 + t572 * t635 + t610 * t686 + t616 * t722 + t703) * MDP(13) + (-g(2) * t675 + t828) * MDP(2) + t740 * MDP(3) + t616 * t792; -t677 * MDP(6) * t845 + (qJDD(2) + t738) * MDP(7) + (t601 * t572 - t600 * t841 + (t832 * t522 - t523 * t671 + (-t572 * t832 + t671 * t841) * qJD(3)) * t666) * MDP(15) + (-g(1) * t812 + t457 * t826 - t488 * t600 + t489 * t601 + (-t487 * t778 - t832 * t456 - t454 * t671 + (t488 * t671 - t489 * t832) * qJD(3)) * t666 + t738) * MDP(18) + (t475 * t814 - t520 * t696 + t724 * t526 - t842 * t852) * MDP(24) + (t474 * t814 + t623 * t520 + t724 * t528 + t843 * t852) * MDP(25) + (-t696 * t449 + t721 * t472 + (-t720 * qJD(6) + t669 * t843 + t673 * t724) * t525 + t842 * t500) * MDP(31) + (-t696 * t448 - t720 * t472 + (-t721 * qJD(6) - t669 * t724 + t673 * t843) * t525 + t842 * t502) * MDP(32) + (MDP(14) - MDP(17)) * (t666 * (t610 * t671 + t611 * t775 - t778 * t841) - t522 * t826 - t601 * t611) + (-MDP(13) + MDP(16)) * (t666 * (t572 * t778 + t610 * t832 - t611 * t798) - t523 * t826 + t600 * t611) + ((t665 * t766 - t786) * MDP(4) + (qJDD(1) * t665 + t668 * t766) * MDP(5) + (-pkin(1) * qJDD(1) + qJD(1) * t727 - t828) * MDP(7)) * t667; t852 * t572 * MDP(23) + (t474 * t674 - t670 * t756) * MDP(19) + (-t449 * t670 + (t669 * t797 + t513) * t525 + (t716 - t846) * t674) * MDP(29) + (t472 * t670 + t525 * t757) * MDP(30) + (-t704 - (qJD(3) + t611) * t841) * MDP(11) + (-t437 * t514 - t452 * t502 + t838 * t669 + t683 * t673 + (t472 * t807 - t427 + (-t437 * t673 - t502 * t833) * qJD(5) + t688 * t669) * t670 + (-t437 * t794 + t430 * t673 - t432 * t841 + t833 * t448 + (t525 * t807 - t432) * qJD(5)) * t674) * MDP(32) + (-t437 * t513 - t452 * t500 - t838 * t673 + t683 * t669 + (t472 * t809 + t428 + (-t437 * t669 - t500 * t833) * qJD(5) - t688 * t673) * t670 + (t437 * t793 + t430 * t669 - t734 * t841 + t833 * t449 + (t525 * t809 - t734) * qJD(5)) * t674) * MDP(31) + (-t487 * t572 + t521 * t841 + t611 * t791 - 0.2e1 * t603 - t605 + t681) * MDP(17) + (-t524 * t841 + t682 - t822) * MDP(13) + (pkin(3) * t522 - qJ(4) * t523 + (-t489 - t496) * t841 + (t488 + t791) * t572) * MDP(15) + (qJ(4) * t475 + t439 * t572 + t790 * t526 + (-t481 * t852 - t840) * t674 + ((t499 + t795) * t852 + t697) * t670) * MDP(24) + (qJ(4) * t474 + t802 * t852 - t440 * t572 + t790 * t528 + t840 * t670 + (t795 * t852 + t697) * t674) * MDP(25) + ((-t475 - t756) * t674 + (-t474 + t856) * t670) * MDP(20) + (t673 * t806 + (-t674 * t794 + t737) * t502) * MDP(26) + MDP(8) * t818 + (-t572 ^ 2 + t834) * MDP(9) + (t500 * t514 + t502 * t513 + (t500 * t673 + t502 * t669) * t797 + (-t823 - t449 * t673 + (t500 * t669 - t502 * t673) * qJD(6)) * t674) * MDP(27) + (t495 * t611 + t524 * t572 - t681) * MDP(14) + t680 * MDP(10) + (t521 * t572 + t678 + t822 + 0.2e1 * t829) * MDP(16) + (-t454 * qJ(4) - t456 * pkin(3) - t487 * t521 - t488 * t496 - g(1) * (-pkin(3) * t536 + qJ(4) * t537) - g(2) * (-pkin(3) * t532 + qJ(4) * t533) - g(3) * (t824 - t830) + t791 * t489) * MDP(18) + (-t526 * t572 + t699) * MDP(22) + (t528 * t572 + t717) * MDP(21) - t792 + (t448 * t670 + t737 * t525 + (t755 + t715) * t674) * MDP(28); t680 * MDP(15) + (-t610 - t818) * MDP(16) + (-t611 ^ 2 - t834) * MDP(17) + (-t489 * t611 + t678 + t829) * MDP(18) + (t526 * t611 + t717) * MDP(24) + (t528 * t611 + t699) * MDP(25) + (-t674 * t449 + (t673 * t611 - t669 * t757) * t525 + (t716 + t846) * t670) * MDP(31) + (-t806 + (-t669 * t611 - t673 * t757) * t525 + (t755 - t715) * t670) * MDP(32); -t526 ^ 2 * MDP(20) + (t474 + t856) * MDP(21) + t758 * MDP(22) - t520 * MDP(23) + (t440 * t852 - t707 - t835) * MDP(24) + (g(1) * t511 + g(2) * t509 + g(3) * t531 + t439 * t852 + t466 * t526 - t709) * MDP(25) + (t502 * t754 + t823) * MDP(26) + ((t448 - t821) * t673 + (-t449 - t820) * t669) * MDP(27) + (t525 * t754 + t810) * MDP(28) + (-t525 ^ 2 * t669 + t808) * MDP(29) + (-pkin(5) * t449 - t440 * t500 + t692 * t669 - t673 * t839) * MDP(31) + (-pkin(5) * t448 - t440 * t502 + t669 * t839 + t692 * t673) * MDP(32) + (t526 * MDP(19) + (-qJD(5) + t852) * MDP(22) - t466 * MDP(24) - t502 * MDP(28) + t500 * MDP(29) - t525 * MDP(30) + t734 * MDP(31) + t432 * MDP(32) + t528 * MDP(20)) * t528; t502 * t500 * MDP(26) + (-t500 ^ 2 + t502 ^ 2) * MDP(27) + (t784 + t821) * MDP(28) + (-t759 + t820) * MDP(29) + t472 * MDP(30) + (-g(1) * t484 + g(2) * t855 + g(3) * t505 + t432 * t525 - t437 * t502 + t761) * MDP(31) + (g(1) * t485 + g(2) * t854 + g(3) * t506 + t437 * t500 - t734 * t525 - t735) * MDP(32) + (-MDP(28) * t819 - MDP(29) * t502 - MDP(31) * t432 + MDP(32) * t734) * qJD(6);];
tau  = t1;
