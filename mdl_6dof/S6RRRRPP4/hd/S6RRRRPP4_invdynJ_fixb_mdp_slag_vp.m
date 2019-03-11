% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRRPP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:02:59
% EndTime: 2019-03-09 21:03:17
% DurationCPUTime: 13.94s
% Computational Cost: add. (14610->718), mult. (32604->890), div. (0->0), fcn. (23013->14), ass. (0->311)
t775 = sin(qJ(4));
t776 = sin(qJ(3));
t779 = cos(qJ(4));
t780 = cos(qJ(3));
t696 = t775 * t776 - t779 * t780;
t928 = qJD(3) + qJD(4);
t626 = t928 * t696;
t781 = cos(qJ(2));
t802 = t696 * t781;
t641 = qJD(1) * t802;
t884 = -t626 + t641;
t697 = t775 * t780 + t776 * t779;
t627 = t928 * t697;
t870 = qJD(1) * t781;
t640 = t697 * t870;
t883 = t627 - t640;
t777 = sin(qJ(2));
t871 = qJD(1) * t777;
t850 = t776 * t871;
t859 = t780 * qJD(2);
t693 = t850 - t859;
t867 = qJD(2) * t776;
t695 = t780 * t871 + t867;
t615 = t779 * t693 + t695 * t775;
t773 = sin(pkin(10));
t774 = cos(pkin(10));
t813 = t693 * t775 - t779 * t695;
t564 = t774 * t615 - t773 * t813;
t935 = -t615 * t773 - t774 * t813;
t951 = pkin(5) * t935 + qJ(6) * t564;
t719 = -qJD(2) * pkin(2) + pkin(7) * t871;
t633 = pkin(3) * t693 + t719;
t579 = pkin(4) * t615 + qJD(5) + t633;
t519 = pkin(5) * t564 - qJ(6) * t935 + t579;
t772 = qJ(3) + qJ(4);
t759 = pkin(10) + t772;
t744 = cos(t759);
t743 = sin(t759);
t782 = cos(qJ(1));
t894 = t782 * t743;
t778 = sin(qJ(1));
t897 = t778 * t781;
t636 = t744 * t897 - t894;
t895 = t781 * t782;
t638 = t743 * t778 + t744 * t895;
t857 = qJD(1) * qJD(2);
t841 = t781 * t857;
t856 = qJDD(1) * t777;
t863 = qJD(3) * t777;
t938 = -qJD(1) * t863 + qJDD(2);
t602 = qJD(3) * t859 + (t841 + t856) * t780 + t938 * t776;
t603 = t776 * (qJD(2) * (qJD(3) + t870) + t856) - t938 * t780;
t860 = qJD(4) * t779;
t861 = qJD(4) * t775;
t545 = t779 * t602 - t775 * t603 - t693 * t860 - t695 * t861;
t758 = t781 * qJDD(1);
t930 = -t777 * t857 + t758;
t690 = qJDD(3) - t930;
t682 = qJDD(4) + t690;
t712 = -pkin(2) * t781 - pkin(8) * t777 - pkin(1);
t683 = t712 * qJD(1);
t753 = pkin(7) * t870;
t720 = qJD(2) * pkin(8) + t753;
t623 = t780 * t683 - t720 * t776;
t590 = -pkin(9) * t695 + t623;
t738 = -qJD(3) + t870;
t583 = -pkin(3) * t738 + t590;
t624 = t683 * t776 + t720 * t780;
t591 = -pkin(9) * t693 + t624;
t589 = t779 * t591;
t550 = t583 * t775 + t589;
t822 = pkin(2) * t777 - pkin(8) * t781;
t705 = t822 * qJD(2);
t628 = qJD(1) * t705 + qJDD(1) * t712;
t620 = t780 * t628;
t664 = pkin(7) * t930 + qJDD(2) * pkin(8);
t535 = pkin(3) * t690 - pkin(9) * t602 - qJD(3) * t624 - t664 * t776 + t620;
t862 = qJD(3) * t780;
t864 = qJD(3) * t776;
t801 = t776 * t628 + t780 * t664 + t683 * t862 - t720 * t864;
t542 = -pkin(9) * t603 + t801;
t836 = t779 * t535 - t775 * t542;
t789 = -qJD(4) * t550 + t836;
t491 = pkin(4) * t682 - qJ(5) * t545 + qJD(5) * t813 + t789;
t788 = qJD(4) * t813 - t602 * t775 - t779 * t603;
t827 = -t775 * t535 - t779 * t542 - t583 * t860 + t591 * t861;
t494 = qJ(5) * t788 - qJD(5) * t615 - t827;
t483 = t773 * t491 + t774 * t494;
t852 = t682 * qJ(6) + t483;
t906 = t744 * t777;
t950 = -g(1) * t638 - g(2) * t636 - g(3) * t906 - t519 * t564 + t852;
t587 = t775 * t591;
t549 = t779 * t583 - t587;
t941 = qJ(5) * t813;
t533 = t549 + t941;
t728 = -qJD(4) + t738;
t528 = -pkin(4) * t728 + t533;
t942 = qJ(5) * t615;
t534 = t550 - t942;
t529 = t774 * t534;
t500 = t773 * t528 + t529;
t498 = -qJ(6) * t728 + t500;
t948 = t498 * t935;
t947 = t500 * t935;
t751 = pkin(7) * t856;
t665 = -qJDD(2) * pkin(2) + pkin(7) * t841 + t751;
t768 = g(3) * t781;
t820 = g(1) * t782 + g(2) * t778;
t932 = t777 * t820;
t794 = -t768 + t932;
t946 = qJD(3) * pkin(8) * t738 - t665 + t794;
t922 = pkin(3) * t776;
t823 = pkin(3) * t864 - t870 * t922 - t753;
t945 = pkin(4) * t883 + t823;
t482 = t774 * t491 - t773 * t494;
t481 = -t682 * pkin(5) + qJDD(6) - t482;
t635 = t743 * t897 + t744 * t782;
t637 = -t778 * t744 + t781 * t894;
t914 = g(3) * t777;
t927 = -g(1) * t637 - g(2) * t635 + t519 * t935 - t743 * t914 + t481;
t944 = pkin(4) * t813;
t940 = t564 * t728;
t889 = t884 * t773 + t883 * t774;
t888 = -t883 * t773 + t884 * t774;
t865 = qJD(2) * t781;
t848 = t776 * t865;
t939 = t777 * t862 + t848;
t761 = sin(t772);
t762 = cos(t772);
t646 = t761 * t782 - t762 * t897;
t648 = t761 * t778 + t762 * t895;
t937 = g(1) * t648 - g(2) * t646 + t615 * t633 + t762 * t914 + t827;
t904 = t762 * t782;
t645 = t761 * t897 + t904;
t853 = t761 * t895;
t905 = t762 * t778;
t647 = -t853 + t905;
t936 = -g(1) * t647 + g(2) * t645 + t633 * t813 + t761 * t914 + t789;
t934 = t682 * MDP(22) + (-t615 ^ 2 + t813 ^ 2) * MDP(19) + (-t615 * t728 + t545) * MDP(20) + (t728 * t813 + t788) * MDP(21) - t615 * MDP(18) * t813;
t933 = t935 ^ 2;
t654 = t697 * t777;
t896 = t780 * t781;
t810 = pkin(3) * t777 - pkin(9) * t896;
t702 = t822 * qJD(1);
t875 = pkin(7) * t850 + t780 * t702;
t601 = qJD(1) * t810 + t875;
t596 = t779 * t601;
t675 = t776 * t702;
t899 = t777 * t780;
t901 = t776 * t781;
t621 = t675 + (-pkin(7) * t899 - pkin(9) * t901) * qJD(1);
t551 = pkin(4) * t871 + qJ(5) * t641 - t621 * t775 + t596;
t923 = pkin(8) + pkin(9);
t851 = qJD(3) * t923;
t703 = t776 * t851;
t704 = t780 * t851;
t721 = t923 * t776;
t722 = t923 * t780;
t799 = -t779 * t703 - t775 * t704 - t721 * t860 - t722 * t861;
t552 = -qJ(5) * t627 - qJD(5) * t696 + t799;
t886 = t775 * t601 + t779 * t621;
t554 = -qJ(5) * t640 + t886;
t678 = t779 * t704;
t876 = -t775 * t721 + t779 * t722;
t786 = qJ(5) * t626 - qJD(4) * t876 - qJD(5) * t697 + t703 * t775 - t678;
t892 = (t551 - t786) * t774 + (t552 - t554) * t773;
t692 = t780 * t712;
t919 = pkin(7) * t776;
t622 = -pkin(9) * t899 + t692 + (-pkin(3) - t919) * t781;
t740 = pkin(7) * t896;
t874 = t776 * t712 + t740;
t902 = t776 * t777;
t629 = -pkin(9) * t902 + t874;
t885 = t775 * t622 + t779 * t629;
t887 = t779 * t590 - t587;
t538 = t887 + t941;
t834 = -t590 * t775 - t589;
t806 = t834 + t942;
t855 = pkin(3) * t773 * t775;
t880 = -qJD(4) * t855 - t773 * t806 + (pkin(3) * t860 - t538) * t774;
t929 = -t776 * t863 + t781 * t859;
t926 = t781 * t820 + t914;
t924 = -0.2e1 * pkin(1);
t921 = pkin(4) * t761;
t920 = pkin(5) * t743;
t918 = g(1) * t778;
t915 = g(2) * t782;
t765 = t780 * pkin(3);
t913 = pkin(2) + t765;
t505 = t533 * t773 + t529;
t912 = t505 * t935;
t911 = t534 * t773;
t910 = t602 * t776;
t909 = t693 * t738;
t908 = t695 * t738;
t907 = t695 * t780;
t903 = t774 * t775;
t900 = t777 * t778;
t898 = t777 * t782;
t708 = pkin(4) * t762 + t765;
t701 = pkin(2) + t708;
t681 = t781 * t701;
t707 = t921 + t922;
t688 = t782 * t707;
t577 = -qJD(2) * t802 - t654 * t928;
t655 = t696 * t777;
t866 = qJD(2) * t777;
t877 = t780 * t705 + t866 * t919;
t571 = t810 * qJD(2) + (-t740 + (pkin(9) * t777 - t712) * t776) * qJD(3) + t877;
t878 = t776 * t705 + t712 * t862;
t576 = -t939 * pkin(9) + (-t777 * t859 - t781 * t864) * pkin(7) + t878;
t835 = t779 * t571 - t576 * t775;
t502 = pkin(4) * t866 - qJ(5) * t577 - qJD(4) * t885 + qJD(5) * t655 + t835;
t578 = -t861 * t902 + (t899 * t928 + t848) * t779 + t929 * t775;
t800 = t775 * t571 + t779 * t576 + t622 * t860 - t629 * t861;
t507 = -qJ(5) * t578 - qJD(5) * t654 + t800;
t488 = t773 * t502 + t774 * t507;
t893 = pkin(5) * t871 + t892;
t517 = t773 * t551 + t774 * t554;
t512 = qJ(6) * t871 + t517;
t515 = t774 * t552 + t773 * t786;
t891 = t515 - t512;
t613 = -t696 * t773 + t697 * t774;
t890 = pkin(5) * t889 - qJ(6) * t888 - qJD(6) * t613 + t945;
t832 = t779 * t622 - t629 * t775;
t556 = -pkin(4) * t781 + qJ(5) * t655 + t832;
t561 = -qJ(5) * t654 + t885;
t525 = t773 * t556 + t774 * t561;
t882 = qJD(6) + t880;
t881 = -t538 * t773 + t774 * t806 + (t773 * t779 + t903) * qJD(4) * pkin(3);
t879 = -t781 * t688 + t778 * t708;
t749 = pkin(3) * t779 + pkin(4);
t659 = pkin(3) * t903 + t773 * t749;
t706 = pkin(3) * t902 + t777 * pkin(7);
t770 = t777 ^ 2;
t873 = -t781 ^ 2 + t770;
t869 = qJD(2) * t693;
t868 = qJD(2) * t695;
t506 = t533 * t774 - t911;
t858 = qJD(6) - t506;
t634 = pkin(3) * t939 + pkin(7) * t865;
t849 = t738 * t859;
t846 = t738 * t864;
t845 = t738 * t862;
t839 = t881 * t935;
t838 = -t635 * pkin(5) + qJ(6) * t636;
t837 = -t637 * pkin(5) + qJ(6) * t638;
t510 = t545 * t773 - t774 * t788;
t831 = -t707 * t897 - t708 * t782;
t830 = -t779 * t721 - t722 * t775;
t829 = -qJD(3) * t683 - t664;
t826 = pkin(4) * t696 - t913;
t741 = g(1) * t900;
t825 = -g(2) * t898 + t741;
t824 = pkin(4) * t654 + t706;
t821 = pkin(3) * t695 - t944;
t819 = t720 * t862 - t620;
t818 = pkin(5) * t744 + qJ(6) * t743;
t817 = -pkin(8) * t690 + qJD(3) * t719;
t487 = t502 * t774 - t507 * t773;
t499 = t528 * t774 - t911;
t511 = t545 * t774 + t773 * t788;
t524 = t556 * t774 - t561 * t773;
t811 = pkin(4) * t578 + t634;
t658 = t749 * t774 - t855;
t807 = -pkin(7) * qJDD(2) + t857 * t924;
t805 = -qJ(5) * t697 + t830;
t804 = t690 * t776 - t845;
t803 = t690 * t780 + t846;
t784 = qJD(1) ^ 2;
t798 = pkin(1) * t784 + t820;
t582 = pkin(3) * t603 + t665;
t783 = qJD(2) ^ 2;
t797 = pkin(7) * t783 + qJDD(1) * t924 + t915;
t769 = -qJ(5) - t923;
t796 = t782 * pkin(1) + t701 * t895 - t769 * t898 + (pkin(7) + t707) * t778;
t795 = t688 + t769 * t900 + t782 * pkin(7) + (-pkin(1) - t681) * t778;
t526 = -pkin(4) * t788 + qJDD(5) + t582;
t597 = -qJ(5) * t696 + t876;
t559 = t597 * t773 - t774 * t805;
t560 = t774 * t597 + t773 * t805;
t787 = -t560 * t510 + t511 * t559 - t515 * t564 - t926;
t484 = pkin(5) * t510 - qJ(6) * t511 - qJD(6) * t935 + t526;
t745 = -pkin(4) * t774 - pkin(5);
t742 = pkin(4) * t773 + qJ(6);
t724 = pkin(4) * t905;
t711 = t728 * qJD(6);
t709 = qJ(6) * t906;
t673 = t776 * t778 + t780 * t895;
t672 = -t776 * t895 + t778 * t780;
t671 = t776 * t782 - t778 * t896;
t670 = t776 * t897 + t780 * t782;
t650 = -pkin(5) - t658;
t649 = qJ(6) + t659;
t612 = t774 * t696 + t697 * t773;
t593 = -t654 * t773 - t655 * t774;
t592 = t774 * t654 - t655 * t773;
t558 = pkin(5) * t612 - qJ(6) * t613 + t826;
t543 = pkin(5) * t592 - qJ(6) * t593 + t824;
t540 = t577 * t774 - t578 * t773;
t539 = t577 * t773 + t774 * t578;
t523 = -t944 + t951;
t522 = pkin(5) * t781 - t524;
t521 = -qJ(6) * t781 + t525;
t520 = t821 + t951;
t497 = pkin(5) * t728 + qJD(6) - t499;
t495 = pkin(5) * t539 - qJ(6) * t540 - qJD(6) * t593 + t811;
t486 = -pkin(5) * t866 - t487;
t485 = qJ(6) * t866 - qJD(6) * t781 + t488;
t480 = -t711 + t852;
t1 = [(-g(1) * t645 - g(2) * t647 + t706 * t545 - t550 * t866 + t633 * t577 - t582 * t655 - t634 * t813 - t682 * t885 + t728 * t800 - t781 * t827) * MDP(24) + (-t545 * t781 - t577 * t728 - t655 * t682 - t813 * t866) * MDP(20) + (-t545 * t655 - t577 * t813) * MDP(18) + (t602 * t899 + t695 * t929) * MDP(11) + (qJDD(2) * t777 + t781 * t783) * MDP(6) + (qJDD(2) * t781 - t777 * t783) * MDP(7) + (t480 * t521 + t498 * t485 + t484 * t543 + t519 * t495 + t481 * t522 + t497 * t486 - g(1) * (-pkin(5) * t636 - qJ(6) * t635 + t795) - g(2) * (pkin(5) * t638 + qJ(6) * t637 + t796)) * MDP(30) + ((t738 * t867 + t603) * t781 + (-t804 - t869) * t777) * MDP(14) + 0.2e1 * (t758 * t777 - t857 * t873) * MDP(5) + (-(-t712 * t864 + t877) * t738 + t692 * t690 - g(1) * t671 - g(2) * t673 + ((t845 + t869) * pkin(7) + (-pkin(7) * t690 + qJD(2) * t719 - t829) * t776 + t819) * t781 + (pkin(7) * t603 + qJD(2) * t623 + t665 * t776 + t719 * t862) * t777) * MDP(16) + (t878 * t738 - t874 * t690 - g(1) * t670 - g(2) * t672 + (t719 * t859 + (-t846 + t868) * pkin(7) + t801) * t781 + (-t719 * t864 - t624 * qJD(2) + t665 * t780 + (t602 - t849) * pkin(7)) * t777) * MDP(17) + (-t835 * t728 + t832 * t682 - t836 * t781 + t549 * t866 + t634 * t615 - t706 * t788 + t582 * t654 + t633 * t578 - g(1) * t646 - g(2) * t648 + (t550 * t781 + t728 * t885) * qJD(4)) * MDP(23) + (g(1) * t635 - g(2) * t637 - t480 * t781 - t484 * t593 - t485 * t728 - t495 * t935 + t498 * t866 - t511 * t543 - t519 * t540 + t521 * t682) * MDP(29) + (-t482 * t593 - t483 * t592 - t487 * t935 - t488 * t564 - t499 * t540 - t500 * t539 - t510 * t525 - t511 * t524 + t825) * MDP(25) + (-t480 * t592 + t481 * t593 - t485 * t564 + t486 * t935 + t497 * t540 - t498 * t539 - t510 * t521 + t511 * t522 + t825) * MDP(28) + (qJDD(1) * t770 + 0.2e1 * t777 * t841) * MDP(4) + (t777 * t797 + t781 * t807 - t741) * MDP(10) + ((-t693 * t780 - t695 * t776) * t865 + (-t910 - t603 * t780 + (t693 * t776 - t907) * qJD(3)) * t777) * MDP(12) + (-t690 * t781 - t738 * t866) * MDP(15) + (-t682 * t781 - t728 * t866) * MDP(22) + (g(1) * t636 - g(2) * t638 + t481 * t781 + t484 * t592 + t486 * t728 + t495 * t564 - t497 * t866 + t510 * t543 + t519 * t539 - t522 * t682) * MDP(27) + ((-t602 - t849) * t781 + (t803 + t868) * t777) * MDP(13) + (-t545 * t654 - t577 * t615 + t578 * t813 - t655 * t788) * MDP(19) + (t578 * t728 - t615 * t866 - t654 * t682 - t781 * t788) * MDP(21) + (t807 * t777 + (-t797 + t918) * t781) * MDP(9) + (-t915 + t918) * MDP(2) + t820 * MDP(3) + (-g(1) * t795 - g(2) * t796 + t482 * t524 + t483 * t525 + t499 * t487 + t500 * t488 + t526 * t824 + t579 * t811) * MDP(26) + qJDD(1) * MDP(1); (t483 * t560 - t482 * t559 + t526 * t826 - g(3) * (-t769 * t777 + t681) + t945 * t579 + (t515 - t517) * t500 - t892 * t499 + t820 * (t701 * t777 + t769 * t781)) * MDP(26) + (-pkin(2) * t603 + t875 * t738 + t817 * t776 + (-t623 * t777 + (-pkin(7) * t693 - t719 * t776) * t781) * qJD(1) + t946 * t780) * MDP(16) + (-pkin(2) * t602 - t675 * t738 + t817 * t780 + (-t719 * t896 + t624 * t777 + (-t695 * t781 + t738 * t899) * pkin(7)) * qJD(1) - t946 * t776) * MDP(17) + (t545 * t697 - t813 * t884) * MDP(18) + (t738 * MDP(15) + MDP(20) * t813 + t615 * MDP(21) + t728 * MDP(22) - t549 * MDP(23) + t550 * MDP(24) + t497 * MDP(27) - t498 * MDP(29)) * t871 + (-g(3) * t681 + t480 * t560 + t481 * t559 + t484 * t558 + t890 * t519 + t891 * t498 + t893 * t497 + (-g(3) * t818 + t769 * t820) * t781 + (g(3) * t769 + t820 * (t701 + t818)) * t777) * MDP(30) + (t777 * t798 - t751 - t768) * MDP(9) + ((-t695 * t777 + t738 * t896) * qJD(1) + t804) * MDP(13) + (-t876 * t682 - t913 * t545 + t582 * t697 + (t799 - t886) * t728 + t884 * t633 - t823 * t813 - t794 * t761) * MDP(24) + (-t484 * t613 - t511 * t558 - t519 * t888 + t560 * t682 - t728 * t891 + t743 * t794 - t890 * t935) * MDP(29) + (-t482 * t613 - t483 * t612 - t499 * t888 - t500 * t889 + t517 * t564 + t892 * t935 + t787) * MDP(25) + (-t480 * t612 + t481 * t613 + t497 * t888 - t498 * t889 + t512 * t564 + t893 * t935 + t787) * MDP(28) + (t830 * t682 + t913 * t788 + t582 * t696 + (t722 * t860 + t596 + t678 + (-qJD(4) * t721 - t621 - t703) * t775) * t728 + t883 * t633 + t823 * t615 + t794 * t762) * MDP(23) + (-t738 * t907 + t910) * MDP(11) + (t682 * t697 - t728 * t884) * MDP(20) + (-t545 * t696 - t615 * t884 + t697 * t788 + t813 * t883) * MDP(19) + (-MDP(4) * t777 * t781 + MDP(5) * t873) * t784 + qJDD(2) * MDP(8) + (t484 * t612 + t510 * t558 + t519 * t889 - t559 * t682 + t564 * t890 + t728 * t893 + t744 * t794) * MDP(27) + ((t602 + t909) * t780 + (-t603 + t908) * t776) * MDP(12) + MDP(7) * t758 + MDP(6) * t856 + (t914 + (-pkin(7) * qJDD(1) + t798) * t781) * MDP(10) + (-t682 * t696 + t728 * t883) * MDP(21) + ((t693 * t777 - t738 * t901) * qJD(1) + t803) * MDP(14); t934 + (t834 * t728 + (-t615 * t695 + t682 * t779 + t728 * t861) * pkin(3) + t936) * MDP(23) + (-t887 * t728 + (-t682 * t775 + t695 * t813 + t728 * t860) * pkin(3) + t937) * MDP(24) + (-t603 - t908) * MDP(14) + (t520 * t935 + t649 * t682 - t728 * t882 - t711 + t950) * MDP(29) + (-t510 * t649 + t511 * t650 + t839 + t948) * MDP(28) + (-t510 * t659 - t511 * t658 + t839 + t947) * MDP(25) + (-t650 * t682 + t728 * t881 - t927) * MDP(27) + (t480 * t649 + t481 * t650 - t519 * t520 - g(1) * (t837 + t879) - g(2) * (t831 + t838) - g(3) * (t709 + (-t707 - t920) * t777) + t882 * t498 + t881 * t497) * MDP(30) + (t602 - t909) * MDP(13) + (-g(1) * t879 - g(2) * t831 + t482 * t658 + t483 * t659 - t499 * t881 + t500 * t880 - t579 * t821 + t707 * t914) * MDP(26) + (-g(1) * t672 + g(2) * t670 - t624 * t738 - t695 * t719 + (t829 + t914) * t776 - t819) * MDP(16) + t695 * t693 * MDP(11) + (-t693 ^ 2 + t695 ^ 2) * MDP(12) + t690 * MDP(15) + (g(1) * t673 - g(2) * t671 + g(3) * t899 - t623 * t738 + t693 * t719 - t801) * MDP(17) + ((t497 - t882) * MDP(28) + (-t499 - t880) * MDP(25) - t520 * MDP(27)) * t564; (-t550 * t728 + t936) * MDP(23) + (-t549 * t728 + t937) * MDP(24) + (t947 - t912 + (-t510 * t773 - t511 * t774) * pkin(4) + (-t499 + t506) * t564) * MDP(25) + (-g(1) * t724 + t499 * t505 - t500 * t506 + (g(2) * t904 + t482 * t774 + t483 * t773 + t579 * t813 + t761 * t926) * pkin(4)) * MDP(26) + (-t505 * t728 - t523 * t564 - t682 * t745 - t927) * MDP(27) + (t948 - t510 * t742 + t511 * t745 - t912 + (t497 - t858) * t564) * MDP(28) + (t506 * t728 + t523 * t935 + t682 * t742 - 0.2e1 * t711 + t950) * MDP(29) + (t480 * t742 + t481 * t745 - t519 * t523 - t497 * t505 - g(1) * (-pkin(4) * t853 + t724 + t837) - g(2) * (-pkin(4) * t645 + t838) - g(3) * (t709 + (-t920 - t921) * t777) + t858 * t498) * MDP(30) + t934; (t499 * t935 + t500 * t564 + t526 + t768) * MDP(26) + (-t728 * t935 + t510) * MDP(27) + (-t511 - t940) * MDP(29) + (-t497 * t935 + t498 * t564 + t484 + t768) * MDP(30) + (-MDP(26) - MDP(30)) * t932 + (MDP(25) + MDP(28)) * (-t564 ^ 2 - t933); (t564 * t935 - t682) * MDP(27) + (t511 - t940) * MDP(28) + (-t728 ^ 2 - t933) * MDP(29) + (t498 * t728 + t927) * MDP(30);];
tau  = t1;
