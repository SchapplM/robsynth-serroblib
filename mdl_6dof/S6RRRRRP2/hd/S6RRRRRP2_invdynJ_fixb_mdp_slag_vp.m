% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP2
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
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:03:36
% EndTime: 2019-03-10 01:03:55
% DurationCPUTime: 14.36s
% Computational Cost: add. (16874->671), mult. (39188->818), div. (0->0), fcn. (29442->14), ass. (0->303)
t717 = cos(qJ(5));
t832 = qJD(5) * t717;
t718 = cos(qJ(3));
t714 = sin(qJ(3));
t715 = sin(qJ(2));
t839 = qJD(1) * t715;
t809 = t714 * t839;
t719 = cos(qJ(2));
t838 = qJD(1) * t719;
t624 = t718 * t838 - t809;
t625 = -t714 * t838 - t718 * t839;
t713 = sin(qJ(4));
t897 = cos(qJ(4));
t903 = t897 * t624 + t713 * t625;
t918 = t903 * t717;
t926 = t832 - t918;
t712 = sin(qJ(5));
t764 = -t713 * t624 + t625 * t897;
t706 = qJDD(2) + qJDD(3);
t792 = qJDD(4) + t706;
t833 = qJD(5) * t712;
t635 = t714 * t715 - t718 * t719;
t754 = t635 * qJDD(1);
t707 = qJD(2) + qJD(3);
t823 = qJDD(1) * t719;
t826 = qJD(1) * qJD(2);
t806 = t719 * t826;
t824 = qJDD(1) * t715;
t759 = -t806 - t824;
t825 = qJD(1) * qJD(3);
t913 = -t719 * t825 + t759;
t575 = -t707 * t809 + t714 * t823 - t718 * t913;
t636 = t714 * t719 + t715 * t718;
t602 = t707 * t636;
t733 = qJD(1) * t602;
t914 = t897 * t575 - t713 * t733;
t724 = -t713 * t754 + t914;
t723 = qJD(4) * t903 + t724;
t802 = qJD(4) + t707;
t915 = qJD(5) * t802 + t723;
t507 = -t712 * t792 - t717 * t915 - t764 * t833;
t503 = t507 * t712;
t508 = t712 * t915 - t717 * t792 - t764 * t832;
t725 = -t754 - t733;
t799 = t713 * t575 - t897 * t725;
t524 = -t764 * qJD(4) + t799;
t523 = qJDD(5) + t524;
t519 = t712 * t523;
t520 = t717 * t523;
t576 = -t712 * t764 - t717 * t802;
t578 = t712 * t802 - t717 * t764;
t912 = -qJD(5) + t903;
t924 = t912 * t712;
t925 = t792 * MDP(22) - t903 ^ 2 * MDP(19) + (-t903 * t707 + t724) * MDP(20) - t799 * MDP(21) + (MDP(18) * t903 + MDP(19) * t764 - MDP(21) * t707 - MDP(29) * t912) * t764 + (t578 * t926 - t503) * MDP(25) + (t578 * t764 - t912 * t926 + t519) * MDP(27) + (-t507 * t717 - t712 * t508 - t576 * t926 + t578 * t924) * MDP(26) + (-t576 * t764 - t912 * t924 + t520) * MDP(28);
t705 = t719 * pkin(2);
t885 = pkin(1) + t705;
t609 = pkin(3) * t635 - t885;
t662 = qJD(1) * t885;
t923 = qJDD(1) * t885;
t620 = t625 * pkin(9);
t898 = pkin(7) + pkin(8);
t664 = t898 * t719;
t647 = qJD(1) * t664;
t626 = t714 * t647;
t663 = t898 * t715;
t645 = qJD(1) * t663;
t883 = qJD(2) * pkin(2);
t632 = -t645 + t883;
t798 = t718 * t632 - t626;
t573 = t620 + t798;
t564 = pkin(3) * t707 + t573;
t630 = t718 * t647;
t772 = -t632 * t714 - t630;
t892 = pkin(9) * t624;
t574 = -t772 + t892;
t570 = t897 * t574;
t532 = t713 * t564 + t570;
t529 = pkin(10) * t802 + t532;
t606 = -pkin(3) * t624 - t662;
t540 = -pkin(4) * t903 + pkin(10) * t764 + t606;
t498 = -t529 * t712 + t540 * t717;
t827 = qJD(6) - t498;
t487 = pkin(5) * t912 + t827;
t711 = qJ(2) + qJ(3);
t704 = qJ(4) + t711;
t693 = sin(t704);
t683 = g(3) * t693;
t694 = cos(t704);
t720 = cos(qJ(1));
t862 = t694 * t720;
t716 = sin(qJ(1));
t863 = t694 * t716;
t816 = g(1) * t862 + g(2) * t863 + t683;
t603 = qJDD(2) * pkin(2) + t759 * t898;
t807 = t715 * t826;
t758 = -t807 + t823;
t605 = t898 * t758;
t745 = qJD(3) * t772 + t718 * t603 - t714 * t605;
t516 = pkin(3) * t706 - pkin(9) * t575 + t745;
t837 = qJD(3) * t714;
t622 = t647 * t837;
t793 = -qJD(3) * t632 - t605;
t522 = -t622 + (pkin(9) * t913 + t603) * t714 + ((-t715 * t825 + t758) * pkin(9) - t793) * t718;
t808 = qJD(4) * t897;
t835 = qJD(4) * t713;
t741 = t713 * t516 + t522 * t897 + t564 * t808 - t574 * t835;
t473 = pkin(10) * t792 + t741;
t684 = pkin(2) * t807;
t881 = qJDD(1) * pkin(1);
t479 = -pkin(2) * t823 - pkin(3) * t725 + t524 * pkin(4) - pkin(10) * t723 + t684 - t881;
t756 = t717 * t473 + t712 * t479 - t529 * t833 + t540 * t832;
t882 = qJ(6) * t523;
t462 = -qJD(6) * t912 + t756 + t882;
t789 = -t712 * t473 + t717 * t479 - t529 * t832 - t540 * t833;
t894 = pkin(5) * t523;
t464 = qJDD(6) - t789 - t894;
t905 = t462 * t717 + t464 * t712;
t922 = t487 * t926 - t816 + t905;
t893 = pkin(5) * t764;
t569 = t713 * t574;
t531 = t564 * t897 - t569;
t528 = -pkin(4) * t802 - t531;
t501 = t576 * pkin(5) - t578 * qJ(6) + t528;
t878 = t501 * t903;
t874 = t528 * t903;
t585 = t764 * qJ(6);
t702 = sin(t711);
t784 = g(1) * t720 + g(2) * t716;
t919 = t784 * t702;
t797 = t645 * t714 - t630;
t581 = t797 - t892;
t844 = -t718 * t645 - t626;
t582 = t620 + t844;
t543 = t713 * t581 + t582 * t897;
t697 = pkin(2) * t718 + pkin(3);
t858 = t713 * t714;
t595 = t697 * t808 + (-t714 * t835 + (t718 * t897 - t858) * qJD(3)) * pkin(2);
t904 = -t595 + t543;
t917 = t904 * t712;
t864 = t693 * t720;
t865 = t693 * t716;
t916 = g(1) * t864 + g(2) * t865;
t780 = t717 * pkin(5) + t712 * qJ(6);
t843 = t916 * t717;
t911 = -t487 * t764 + t501 * t833 + t843;
t499 = t529 * t717 + t540 * t712;
t488 = -qJ(6) * t912 + t499;
t887 = g(3) * t712;
t817 = -t694 * t887 + t712 * t916;
t790 = -t897 * t516 + t713 * t522 + t564 * t835 + t574 * t808;
t474 = -pkin(4) * t792 + t790;
t465 = t508 * pkin(5) + t507 * qJ(6) - t578 * qJD(6) + t474;
t880 = t465 * t712;
t910 = t488 * t764 + t817 - t880;
t909 = t498 * t764 + t528 * t833 + t843;
t753 = t474 * t712 - t499 * t764 + t528 * t832 - t817;
t731 = -t606 * t903 - t741 + t816;
t888 = g(3) * t694;
t746 = t606 * t764 - t790 - t888 + t916;
t558 = -pkin(4) * t764 - pkin(10) * t903;
t536 = t713 * t573 + t570;
t787 = pkin(3) * t835 - t536;
t600 = -t713 * t635 + t636 * t897;
t763 = -t635 * t897 - t713 * t636;
t554 = -pkin(4) * t763 - pkin(10) * t600 + t609;
t796 = -t718 * t663 - t664 * t714;
t586 = -pkin(9) * t636 + t796;
t842 = -t714 * t663 + t718 * t664;
t587 = -pkin(9) * t635 + t842;
t556 = t713 * t586 + t587 * t897;
t849 = t712 * t554 + t717 * t556;
t810 = t897 * t714;
t846 = t581 * t897 - t713 * t582 + t697 * t835 + (t714 * t808 + (t713 * t718 + t810) * qJD(3)) * pkin(2);
t902 = t694 * pkin(4) + t693 * pkin(10);
t901 = t625 * t835 + t624 * t808 + t713 * (-t714 * t824 + t718 * t823) + t914;
t619 = pkin(5) * t833 - qJ(6) * t832 - t712 * qJD(6);
t869 = t903 * t712;
t900 = pkin(5) * t869 - qJ(6) * t918 - t619;
t899 = t578 ^ 2;
t896 = pkin(3) * t625;
t891 = pkin(10) * t523;
t884 = pkin(10) * qJD(5);
t879 = t499 * t912;
t877 = t508 * t717;
t841 = pkin(2) * t810 + t713 * t697;
t618 = pkin(10) + t841;
t876 = t523 * t618;
t695 = pkin(3) * t713 + pkin(10);
t875 = t523 * t695;
t873 = t576 * t712;
t872 = t578 * t576;
t871 = t578 * t717;
t867 = t595 * t717;
t866 = t600 * t717;
t860 = t712 * t716;
t857 = t716 * t717;
t856 = t717 * t720;
t855 = t720 * t712;
t854 = t900 - t787;
t852 = t717 * t531 + t712 * t558;
t537 = t573 * t897 - t569;
t548 = t558 - t896;
t851 = t717 * t537 + t712 * t548;
t699 = pkin(2) * t839;
t544 = t548 + t699;
t850 = t717 * t543 + t712 * t544;
t847 = -t900 + t846;
t779 = pkin(5) * t712 - qJ(6) * t717;
t845 = -t779 * t903 - t532 + t619;
t709 = t715 ^ 2;
t840 = -t719 ^ 2 + t709;
t836 = qJD(3) * t718;
t834 = qJD(5) * t501;
t822 = t897 * pkin(3);
t701 = t715 * t883;
t814 = qJD(2) * t898;
t813 = t712 * t897;
t812 = t717 * t897;
t594 = pkin(3) * t602 + t701;
t804 = -t465 - t888;
t803 = -t474 - t888;
t492 = -t585 + t850;
t801 = -t492 + t867;
t800 = t544 * t717 - t893 - t917;
t791 = pkin(3) * t808;
t788 = t694 * t780 + t902;
t612 = t694 * t860 + t856;
t614 = t694 * t855 - t857;
t786 = -g(1) * t612 + g(2) * t614;
t613 = t694 * t857 - t855;
t615 = t694 * t856 + t860;
t785 = g(1) * t613 - g(2) * t615;
t783 = g(1) * t716 - g(2) * t720;
t781 = -pkin(2) * t858 + t697 * t897;
t703 = cos(t711);
t688 = pkin(3) * t703;
t777 = t688 + t788;
t775 = t487 * t717 - t488 * t712;
t774 = -t874 - t876;
t773 = -t874 - t875;
t771 = pkin(4) + t780;
t770 = t688 + t885 + t902;
t767 = -0.2e1 * pkin(1) * t826 - pkin(7) * qJDD(2);
t766 = t586 * t897 - t713 * t587;
t601 = t707 * t635;
t545 = qJD(4) * t763 - t601 * t897 - t713 * t602;
t762 = t545 * t712 + t600 * t832;
t761 = -t545 * t717 + t600 * t833;
t760 = -t618 * t833 + t867;
t646 = t715 * t814;
t648 = t719 * t814;
t757 = -t718 * t646 - t714 * t648 - t663 * t836 - t664 * t837;
t551 = -pkin(9) * t602 + t757;
t744 = -qJD(3) * t842 + t646 * t714 - t718 * t648;
t552 = pkin(9) * t601 + t744;
t485 = qJD(4) * t766 + t551 * t897 + t713 * t552;
t546 = qJD(4) * t600 - t713 * t601 + t602 * t897;
t495 = pkin(4) * t546 - pkin(10) * t545 + t594;
t755 = t717 * t485 + t712 * t495 + t554 * t832 - t556 * t833;
t751 = t501 * t918 + t910;
t750 = -t695 * t833 + t717 * t791;
t721 = qJD(2) ^ 2;
t749 = -pkin(7) * t721 + t783 + 0.2e1 * t881;
t722 = qJD(1) ^ 2;
t748 = pkin(1) * t722 - pkin(7) * qJDD(1) + t784;
t747 = g(1) * t614 + g(2) * t612 + t693 * t887 + t789;
t743 = t804 * t717 + t911;
t742 = t803 * t717 + t909;
t740 = qJD(5) * t775 + t905;
t739 = -t503 - t877 + (t871 + t873) * qJD(5);
t738 = t501 * t578 + qJDD(6) - t747;
t737 = (-t833 + t869) * t488 + t922;
t736 = -g(1) * t615 - g(2) * t613 - t683 * t717 + t756;
t486 = qJD(4) * t556 + t713 * t551 - t552 * t897;
t735 = g(3) * t702 - t714 * t603 + t662 * t624 + t703 * t784 + t793 * t718 + t622;
t730 = -g(3) * t703 - t662 * t625 + t745 + t919;
t729 = t784 * t771 * t693;
t727 = t625 * t624 * MDP(11) + (-t624 * t707 + t575) * MDP(13) + (-t625 * t707 + t725) * MDP(14) + (-t624 ^ 2 + t625 ^ 2) * MDP(12) + t706 * MDP(15) + t925;
t652 = pkin(10) * t863;
t654 = pkin(10) * t862;
t726 = -g(1) * t654 - g(2) * t652 + t729;
t708 = -pkin(9) - t898;
t696 = -t822 - pkin(4);
t649 = -pkin(2) * t715 - pkin(3) * t702;
t633 = -t822 - t771;
t621 = t684 - t923;
t617 = -pkin(4) - t781;
t607 = t699 - t896;
t604 = -t771 - t781;
t561 = pkin(3) * t733 + qJDD(1) * t609 + t684;
t541 = pkin(5) * t578 + qJ(6) * t576;
t514 = t600 * t779 - t766;
t506 = pkin(5) * t763 - t554 * t717 + t556 * t712;
t505 = -qJ(6) * t763 + t849;
t497 = t531 * t712 - t558 * t717 + t893;
t496 = t852 - t585;
t491 = t537 * t712 - t548 * t717 + t893;
t490 = -t585 + t851;
t482 = -t576 * t912 - t507;
t469 = t779 * t545 + (qJD(5) * t780 - qJD(6) * t717) * t600 + t486;
t467 = -pkin(5) * t546 + qJD(5) * t849 + t485 * t712 - t495 * t717;
t466 = qJ(6) * t546 - qJD(6) * t763 + t755;
t1 = [(-g(1) * t865 + g(2) * t864 - t485 * t802 + t606 * t545 - t556 * t792 + t561 * t600 - t594 * t764 + t609 * t723) * MDP(24) + (-t545 * t764 + t600 * t901) * MDP(18) + (-t486 * t802 + t609 * t524 + t606 * t546 - t561 * t763 - t594 * t903 + t694 * t783 + t766 * t792) * MDP(23) + (-t600 * t524 + t545 * t903 + t546 * t764 + t763 * t901) * MDP(19) + (-t546 * t802 + t763 * t792) * MDP(21) + (-t601 * t707 + t636 * t706) * MDP(13) + (t575 * t636 + t601 * t625) * MDP(11) + (-t575 * t635 - t601 * t624 + t625 * t602 + t636 * t725) * MDP(12) + (t507 * t763 + t520 * t600 + t546 * t578 + t761 * t912) * MDP(27) + (t464 * t763 + t467 * t912 + t469 * t576 - t487 * t546 + t501 * t762 - t506 * t523 + t508 * t514 + t600 * t880 + t785) * MDP(32) + (-t523 * t763 - t546 * t912) * MDP(29) + (-t462 * t763 - t465 * t866 - t466 * t912 - t469 * t578 + t488 * t546 + t501 * t761 + t505 * t523 + t507 * t514 - t786) * MDP(34) + (t508 * t763 - t519 * t600 - t546 * t576 + t762 * t912) * MDP(28) + (t474 * t866 + t486 * t578 - t499 * t546 + t507 * t766 - t523 * t849 - t528 * t761 + t755 * t912 + t756 * t763 + t786) * MDP(31) + (-t789 * t763 + t498 * t546 + t486 * t576 - t766 * t508 + (-(-qJD(5) * t556 + t495) * t912 + t554 * t523 + t528 * qJD(5) * t600) * t717 + (-(-qJD(5) * t554 - t485) * t912 - t556 * t523 + t474 * t600 + t528 * t545) * t712 + t785) * MDP(30) + (-t575 * t885 + t662 * t601 + t621 * t636 - t625 * t701 - t702 * t783 - t706 * t842 - t707 * t757) * MDP(17) + (-t624 * t701 + t703 * t783 + t706 * t796 + t707 * t744 + (t621 - t923) * t635 - 0.2e1 * t662 * t602) * MDP(16) + t783 * MDP(2) + (-t466 * t576 + t467 * t578 - t505 * t508 - t506 * t507 + t783 * t693 + t775 * t545 + (-t462 * t712 + t464 * t717 + (-t487 * t712 - t488 * t717) * qJD(5)) * t600) * MDP(33) + t784 * MDP(3) + (t715 * t767 + t719 * t749) * MDP(9) + (-t715 * t749 + t719 * t767) * MDP(10) + qJDD(1) * MDP(1) + 0.2e1 * (t715 * t823 - t826 * t840) * MDP(5) + (t545 * t802 + t600 * t792) * MDP(20) + (qJDD(1) * t709 + 0.2e1 * t715 * t806) * MDP(4) + (qJDD(2) * t715 + t719 * t721) * MDP(6) + (qJDD(2) * t719 - t715 * t721) * MDP(7) + (-t602 * t707 - t635 * t706) * MDP(14) + ((-t576 * t717 - t578 * t712) * t545 + (t503 - t877 + (-t871 + t873) * qJD(5)) * t600) * MDP(26) + (-t507 * t866 - t578 * t761) * MDP(25) + (t462 * t505 + t488 * t466 + t465 * t514 + t501 * t469 + t464 * t506 + t487 * t467 - g(1) * (-pkin(5) * t613 - qJ(6) * t612) - g(2) * (pkin(5) * t615 + qJ(6) * t614) + (g(1) * t708 - g(2) * t770) * t720 + (g(1) * t770 + g(2) * t708) * t716) * MDP(35); (t617 * t508 + t774 * t712 + t846 * t576 - ((-qJD(5) * t618 - t544) * t717 + t917) * t912 + t742) * MDP(30) + t727 + MDP(6) * t824 + qJDD(2) * MDP(8) + (t508 * t604 + (-t876 - t878) * t712 + t847 * t576 - (-t618 * t832 - t800) * t912 + t743) * MDP(32) + (-t576 * t801 + t578 * t800 + t618 * t739 + t737) * MDP(33) + (t844 * t707 + (t625 * t839 - t706 * t714 - t707 * t836) * pkin(2) + t735) * MDP(17) + (-t797 * t707 + (t624 * t839 + t706 * t718 - t707 * t837) * pkin(2) + t730) * MDP(16) + (-t617 * t507 + t774 * t717 + t846 * t578 - (-t760 + t850) * t912 + t753) * MDP(31) + (-g(3) * t719 + t715 * t748) * MDP(9) + (g(3) * t715 + t719 * t748) * MDP(10) + (t465 * t604 - g(1) * (t649 * t720 + t654) - g(2) * (t649 * t716 + t652) - g(3) * (t705 + t777) + t847 * t501 + t801 * t488 + t800 * t487 + t729 + t740 * t618) * MDP(35) + (t607 * t764 - t792 * t841 + t731) * MDP(24) + (t607 * t903 + t781 * t792 + t746) * MDP(23) + (t507 * t604 + (-t834 + t876) * t717 - t847 * t578 - (-t492 + t760) * t912 + t751) * MDP(34) + MDP(7) * t823 + (-MDP(23) * t846 + MDP(24) * t904) * t802 + (-t715 * t719 * MDP(4) + MDP(5) * t840) * t722; (t696 * t508 + t773 * t712 + t787 * t576 - ((-qJD(5) * t695 - t548) * t717 + (-t791 + t537) * t712) * t912 + t742) * MDP(30) + (-t696 * t507 + t773 * t717 + t787 * t578 - (-t750 + t851) * t912 + t753) * MDP(31) + (-t707 * t772 + t730) * MDP(16) + (t707 * t798 + t735) * MDP(17) + (t633 * t507 + (-t834 + t875) * t717 + t854 * t578 - (-t490 + t750) * t912 + t751) * MDP(34) + t727 + (t465 * t633 - t488 * t490 - t487 * t491 - g(3) * t777 - t854 * t501 + (t919 + (t487 * t813 + t488 * t812) * qJD(4)) * pkin(3) + t740 * t695 + t726) * MDP(35) + (t633 * t508 + (-t875 - t878) * t712 - t854 * t576 - (-t695 * t832 - t712 * t791 + t491) * t912 + t743) * MDP(32) + (t536 * t802 + (-t625 * t903 + t792 * t897 - t802 * t835) * pkin(3) + t746) * MDP(23) + (t537 * t802 + (-t625 * t764 - t713 * t792 - t802 * t808) * pkin(3) + t731) * MDP(24) + (t490 * t576 - t491 * t578 + (-t576 * t812 + t578 * t813) * qJD(4) * pkin(3) + t739 * t695 + t737) * MDP(33); (t531 * t802 + t731) * MDP(24) + (t532 * t802 + t746) * MDP(23) + (-pkin(4) * t508 - t532 * t576 + (-t531 * t912 - t874 - t891) * t712 + (-(-t558 - t884) * t912 + t803) * t717 + t909) * MDP(30) + (pkin(4) * t507 - t852 * t912 - t532 * t578 - t528 * t918 + (-t833 * t912 - t520) * pkin(10) + t753) * MDP(31) + (-t497 * t912 - t508 * t771 + (-t891 - t878) * t712 + t845 * t576 + (t884 * t912 + t804) * t717 + t911) * MDP(32) + (-t507 * t771 - (-pkin(10) * t833 - t496) * t912 - t845 * t578 + (t501 * t912 + t891) * t717 + t910) * MDP(34) + (pkin(10) * t740 - g(3) * t788 - t465 * t771 - t487 * t497 - t488 * t496 + t501 * t845 + t726) * MDP(35) + (pkin(10) * t739 + t488 * t924 + t496 * t576 - t497 * t578 + t922) * MDP(33) + t925; MDP(25) * t872 + (-t576 ^ 2 + t899) * MDP(26) + t482 * MDP(27) + (-t578 * t912 - t508) * MDP(28) + t523 * MDP(29) + (-t528 * t578 + t747 - t879) * MDP(30) + (-t498 * t912 + t528 * t576 - t736) * MDP(31) + (-t541 * t576 - t738 - t879 + 0.2e1 * t894) * MDP(32) + (pkin(5) * t507 - qJ(6) * t508 + (t488 - t499) * t578 + (t487 - t827) * t576) * MDP(33) + (0.2e1 * t882 - t501 * t576 + t541 * t578 - (0.2e1 * qJD(6) - t498) * t912 + t736) * MDP(34) + (t462 * qJ(6) - t464 * pkin(5) - t501 * t541 - t487 * t499 - g(1) * (-pkin(5) * t614 + qJ(6) * t615) - g(2) * (-pkin(5) * t612 + qJ(6) * t613) + t779 * t683 + t827 * t488) * MDP(35); (t872 - t523) * MDP(32) + t482 * MDP(33) + (-t912 ^ 2 - t899) * MDP(34) + (t488 * t912 + t738 - t894) * MDP(35);];
tau  = t1;
