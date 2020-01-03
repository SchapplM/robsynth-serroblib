% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRRPR2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:49
% EndTime: 2020-01-03 12:07:58
% DurationCPUTime: 10.02s
% Computational Cost: add. (59094->416), mult. (74873->563), div. (0->0), fcn. (44118->10), ass. (0->285)
t831 = qJD(1) + qJD(2);
t824 = qJD(3) + t831;
t822 = t824 ^ 2;
t830 = qJDD(1) + qJDD(2);
t823 = qJDD(3) + t830;
t835 = sin(pkin(9));
t836 = cos(pkin(9));
t782 = t822 * t836 + t823 * t835;
t785 = t822 * t835 - t823 * t836;
t838 = sin(qJ(3));
t842 = cos(qJ(3));
t717 = t782 * t842 - t785 * t838;
t834 = g(1) - qJDD(4);
t758 = qJ(4) * t782 - t834 * t836;
t911 = qJ(4) * t785 - t834 * t835;
t661 = pkin(7) * t717 + t758 * t842 - t838 * t911;
t721 = t782 * t838 + t785 * t842;
t839 = sin(qJ(2));
t843 = cos(qJ(2));
t671 = t717 * t843 - t721 * t839;
t928 = pkin(7) * t721 + t758 * t838 + t842 * t911;
t601 = pkin(6) * t671 + t661 * t843 - t839 * t928;
t840 = sin(qJ(1));
t844 = cos(qJ(1));
t675 = t717 * t839 + t721 * t843;
t932 = t671 * t840 + t675 * t844;
t946 = pkin(6) * t675 + t661 * t839 + t843 * t928;
t958 = pkin(5) * t932 + t840 * t601 + t844 * t946;
t628 = t671 * t844 - t675 * t840;
t957 = -pkin(5) * t628 - t844 * t601 + t840 * t946;
t816 = g(2) * t844 + g(3) * t840;
t804 = qJDD(1) * pkin(1) - t816;
t815 = g(2) * t840 - g(3) * t844;
t847 = qJD(1) ^ 2;
t805 = -pkin(1) * t847 - t815;
t743 = -t804 * t843 + t805 * t839;
t732 = pkin(2) * t830 - t743;
t744 = t804 * t839 + t805 * t843;
t829 = t831 ^ 2;
t733 = -pkin(2) * t829 + t744;
t685 = t732 * t838 + t733 * t842;
t682 = -pkin(3) * t822 + t685;
t684 = -t732 * t842 + t733 * t838;
t852 = pkin(3) * t823 - t684;
t634 = t682 * t835 - t836 * t852;
t635 = t836 * t682 + t835 * t852;
t872 = t634 * t835 + t635 * t836;
t590 = t634 * t836 - t635 * t835;
t884 = t842 * t590;
t566 = -t838 * t872 + t884;
t888 = t838 * t590;
t918 = t842 * t872 + t888;
t550 = t566 * t843 - t839 * t918;
t943 = t566 * t839 + t843 * t918;
t956 = t550 * t840 + t844 * t943;
t955 = -t550 * t844 + t840 * t943;
t788 = t822 * t842 + t823 * t838;
t791 = t822 * t838 - t823 * t842;
t725 = t788 * t843 - t791 * t839;
t762 = pkin(7) * t788 - g(1) * t842;
t912 = pkin(7) * t791 - g(1) * t838;
t669 = pkin(6) * t725 + t762 * t843 - t839 * t912;
t729 = t788 * t839 + t791 * t843;
t909 = t725 * t840 + t729 * t844;
t929 = pkin(6) * t729 + t762 * t839 + t843 * t912;
t948 = pkin(5) * t909 + t840 * t669 + t844 * t929;
t680 = t725 * t844 - t729 * t840;
t947 = -pkin(5) * t680 - t844 * t669 + t840 * t929;
t871 = t684 * t838 + t685 * t842;
t640 = t684 * t842 - t685 * t838;
t883 = t843 * t640;
t594 = -t839 * t871 + t883;
t887 = t839 * t640;
t919 = t843 * t871 + t887;
t942 = t594 * t840 + t844 * t919;
t941 = -t594 * t844 + t840 * t919;
t797 = t829 * t843 + t830 * t839;
t767 = pkin(6) * t797 - g(1) * t843;
t800 = t829 * t839 - t830 * t843;
t855 = t797 * t840 + t800 * t844;
t913 = pkin(6) * t800 - g(1) * t839;
t931 = pkin(5) * t855 + t840 * t767 + t844 * t913;
t735 = t797 * t844 - t800 * t840;
t930 = -pkin(5) * t735 - t844 * t767 + t840 * t913;
t870 = t743 * t839 + t744 * t843;
t692 = t743 * t843 - t744 * t839;
t882 = t844 * t692;
t921 = t840 * t870 - t882;
t886 = t840 * t692;
t920 = t844 * t870 + t886;
t630 = -t822 * pkin(4) + t823 * pkin(8) + t635;
t837 = sin(qJ(5));
t841 = cos(qJ(5));
t622 = t630 * t837 + t834 * t841;
t623 = t630 * t841 - t834 * t837;
t585 = t622 * t837 + t623 * t841;
t832 = t837 ^ 2;
t894 = t832 * t822;
t629 = -pkin(4) * t823 - pkin(8) * t822 + t634;
t624 = t837 * t629;
t810 = t841 * t822 * t837;
t801 = qJDD(5) + t810;
t891 = t837 * t801;
t802 = qJDD(5) - t810;
t890 = t837 * t802;
t889 = t837 * t823;
t625 = t841 * t629;
t885 = t841 * t802;
t817 = t841 * t823;
t587 = pkin(3) * t590;
t881 = -pkin(2) * t566 - t587;
t880 = -pkin(4) * t629 + pkin(8) * t585;
t833 = t841 ^ 2;
t879 = t832 + t833;
t878 = qJD(5) * t824;
t846 = qJD(5) ^ 2;
t807 = -t846 - t894;
t752 = -t807 * t837 - t885;
t814 = t841 * t878;
t777 = 0.2e1 * t814 + t889;
t877 = -pkin(4) * t777 + pkin(8) * t752 + t624;
t818 = t833 * t822;
t809 = -t818 - t846;
t750 = t809 * t841 - t891;
t874 = t837 * t878;
t780 = t817 - 0.2e1 * t874;
t876 = pkin(4) * t780 + pkin(8) * t750 - t625;
t573 = t585 * t835 - t629 * t836;
t875 = pkin(3) * t573 + t880;
t811 = -qJDD(1) * t840 - t844 * t847;
t873 = pkin(5) * t811 + g(1) * t844;
t867 = -t815 * t840 - t816 * t844;
t866 = t835 * t810;
t865 = t836 * t810;
t574 = t585 * t836 + t629 * t835;
t556 = t573 * t842 + t574 * t838;
t864 = pkin(2) * t556 + t875;
t786 = t879 * t823;
t792 = t818 + t894;
t863 = pkin(4) * t792 + pkin(8) * t786 + t585;
t699 = t752 * t835 - t777 * t836;
t862 = pkin(3) * t699 + t877;
t698 = t750 * t835 + t780 * t836;
t861 = pkin(3) * t698 + t876;
t860 = -pkin(3) * t785 - t634;
t859 = -pkin(2) * t791 - t684;
t723 = t786 * t835 + t792 * t836;
t858 = pkin(3) * t723 + t863;
t700 = t750 * t836 - t780 * t835;
t650 = t698 * t842 + t700 * t838;
t857 = pkin(2) * t650 + t861;
t701 = t752 * t836 + t777 * t835;
t651 = t699 * t842 + t701 * t838;
t856 = pkin(2) * t651 + t862;
t584 = t622 * t841 - t623 * t837;
t854 = t815 * t844 - t816 * t840;
t853 = -pkin(2) * t721 + t860;
t724 = t786 * t836 - t792 * t835;
t677 = t723 * t842 + t724 * t838;
t851 = pkin(2) * t677 + t858;
t850 = -pkin(2) * t788 - t685;
t849 = -pkin(3) * t782 - t635;
t848 = -pkin(2) * t717 + t849;
t845 = pkin(1) * g(1);
t812 = qJDD(1) * t844 - t840 * t847;
t808 = t818 - t846;
t806 = t846 - t894;
t795 = t841 * t801;
t794 = pkin(5) * t812 + g(1) * t840;
t793 = -t818 + t894;
t779 = t817 - t874;
t778 = t814 + t889;
t771 = t879 * t878;
t754 = qJDD(5) * t835 + t771 * t836;
t753 = -qJDD(5) * t836 + t771 * t835;
t751 = -t806 * t837 + t795;
t749 = t808 * t841 - t890;
t748 = t807 * t841 - t890;
t747 = t806 * t841 + t891;
t746 = t809 * t837 + t795;
t745 = t808 * t837 + t885;
t742 = t778 * t841 - t832 * t878;
t741 = -t779 * t837 - t833 * t878;
t737 = (t778 + t814) * t837;
t736 = (t779 - t874) * t841;
t715 = -t777 * t837 + t780 * t841;
t714 = t777 * t841 + t780 * t837;
t711 = t751 * t836 + t835 * t889;
t710 = t749 * t836 + t817 * t835;
t709 = t751 * t835 - t836 * t889;
t708 = t749 * t835 - t817 * t836;
t707 = -pkin(1) * t800 - t743;
t706 = -pkin(1) * t797 - t744;
t705 = t742 * t836 - t866;
t704 = t741 * t836 + t866;
t703 = t742 * t835 + t865;
t702 = t741 * t835 - t865;
t695 = -t753 * t838 + t754 * t842;
t694 = t753 * t842 + t754 * t838;
t691 = t715 * t836 + t793 * t835;
t688 = t715 * t835 - t793 * t836;
t687 = pkin(1) * t692;
t686 = pkin(6) * t870 + t845;
t678 = -t723 * t838 + t724 * t842;
t665 = -t709 * t838 + t711 * t842;
t664 = -t708 * t838 + t710 * t842;
t663 = t709 * t842 + t711 * t838;
t662 = t708 * t842 + t710 * t838;
t657 = -t703 * t838 + t705 * t842;
t656 = -t702 * t838 + t704 * t842;
t655 = t703 * t842 + t705 * t838;
t654 = t702 * t842 + t704 * t838;
t653 = -t699 * t838 + t701 * t842;
t652 = -t698 * t838 + t700 * t842;
t647 = -t694 * t839 + t695 * t843;
t646 = t694 * t843 + t695 * t839;
t645 = -pkin(1) * t729 + t859;
t644 = -pkin(1) * t725 + t850;
t643 = -t688 * t838 + t691 * t842;
t642 = t688 * t842 + t691 * t838;
t637 = pkin(2) * t640;
t636 = pkin(2) * g(1) + pkin(7) * t871;
t632 = -t677 * t839 + t678 * t843;
t631 = t677 * t843 + t678 * t839;
t619 = -t663 * t839 + t665 * t843;
t618 = -t662 * t839 + t664 * t843;
t617 = t663 * t843 + t665 * t839;
t616 = t662 * t843 + t664 * t839;
t615 = -t655 * t839 + t657 * t843;
t614 = -t654 * t839 + t656 * t843;
t613 = t655 * t843 + t657 * t839;
t612 = t654 * t843 + t656 * t839;
t611 = -pkin(8) * t748 + t625;
t610 = -pkin(8) * t746 + t624;
t609 = -pkin(4) * t748 + t623;
t608 = -pkin(4) * t746 + t622;
t607 = -t651 * t839 + t653 * t843;
t606 = -t650 * t839 + t652 * t843;
t605 = t651 * t843 + t653 * t839;
t604 = t650 * t843 + t652 * t839;
t603 = -t642 * t839 + t643 * t843;
t602 = t642 * t843 + t643 * t839;
t597 = -pkin(1) * t675 + t853;
t596 = -pkin(1) * t671 + t848;
t586 = pkin(3) * t834 + qJ(4) * t872;
t581 = -qJ(4) * t723 + t584 * t836;
t580 = qJ(4) * t724 + t584 * t835;
t579 = -pkin(1) * t594 - t637;
t578 = -qJ(4) * t699 - t609 * t835 + t611 * t836;
t577 = -qJ(4) * t698 - t608 * t835 + t610 * t836;
t576 = -pkin(3) * t748 + qJ(4) * t701 + t609 * t836 + t611 * t835;
t575 = -pkin(3) * t746 + qJ(4) * t700 + t608 * t836 + t610 * t835;
t571 = pkin(1) * t605 + t856;
t570 = pkin(1) * t604 + t857;
t569 = pkin(6) * t594 + pkin(7) * t883 - t636 * t839;
t568 = pkin(6) * t919 + pkin(7) * t887 + t636 * t843 + t845;
t562 = pkin(1) * t631 + t851;
t561 = -pkin(7) * t677 - t580 * t838 + t581 * t842;
t560 = pkin(7) * t678 + t580 * t842 + t581 * t838;
t559 = -pkin(7) * t651 - t576 * t838 + t578 * t842;
t558 = -pkin(7) * t650 - t575 * t838 + t577 * t842;
t557 = -t573 * t838 + t574 * t842;
t554 = -pkin(2) * t748 + pkin(7) * t653 + t576 * t842 + t578 * t838;
t553 = -pkin(2) * t746 + pkin(7) * t652 + t575 * t842 + t577 * t838;
t552 = -qJ(4) * t573 - (pkin(4) * t835 - pkin(8) * t836) * t584;
t547 = pkin(7) * t566 + qJ(4) * t884 - t586 * t838;
t546 = pkin(2) * t834 + pkin(7) * t918 + qJ(4) * t888 + t586 * t842;
t545 = qJ(4) * t574 - (-pkin(4) * t836 - pkin(8) * t835 - pkin(3)) * t584;
t544 = -pkin(6) * t631 - t560 * t839 + t561 * t843;
t543 = pkin(6) * t632 + t560 * t843 + t561 * t839;
t542 = -t556 * t839 + t557 * t843;
t541 = t556 * t843 + t557 * t839;
t540 = -pkin(1) * t550 + t881;
t539 = -pkin(6) * t605 - t554 * t839 + t559 * t843;
t538 = -pkin(6) * t604 - t553 * t839 + t558 * t843;
t537 = -pkin(1) * t748 + pkin(6) * t607 + t554 * t843 + t559 * t839;
t536 = -pkin(1) * t746 + pkin(6) * t606 + t553 * t843 + t558 * t839;
t535 = pkin(6) * t550 - t546 * t839 + t547 * t843;
t534 = pkin(1) * t834 + pkin(6) * t943 + t546 * t843 + t547 * t839;
t533 = -pkin(7) * t556 - t545 * t838 + t552 * t842;
t532 = pkin(1) * t541 + t864;
t531 = pkin(2) * t584 + pkin(7) * t557 + t545 * t842 + t552 * t838;
t530 = -pkin(6) * t541 - t531 * t839 + t533 * t843;
t529 = pkin(1) * t584 + pkin(6) * t542 + t531 * t843 + t533 * t839;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t816, t815, 0, 0, 0, 0, 0, 0, 0, t830, t707, t706, 0, -t687, 0, 0, 0, 0, 0, t823, t645, t644, 0, t579, 0, 0, 0, 0, 0, t823, t597, t596, 0, t540, t737, t714, t747, t736, t745, 0, t570, t571, t562, t532; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t811, 0, t812, 0, t873, -t794, -t854, -pkin(5) * t854, 0, 0, t735, 0, -t855, 0, t930, t931, t920, pkin(5) * t920 + pkin(6) * t886 + t844 * t686, 0, 0, t680, 0, -t909, 0, t947, t948, t942, pkin(5) * t942 + t844 * t568 + t840 * t569, 0, 0, t628, 0, -t932, 0, t957, t958, t956, pkin(5) * t956 + t844 * t534 + t840 * t535, t613 * t844 + t615 * t840, t602 * t844 + t603 * t840, t617 * t844 + t619 * t840, t612 * t844 + t614 * t840, t616 * t844 + t618 * t840, t646 * t844 + t647 * t840, t840 * t538 + t844 * t536 - pkin(5) * (t604 * t840 - t606 * t844), t840 * t539 + t844 * t537 - pkin(5) * (t605 * t840 - t607 * t844), t840 * t544 + t844 * t543 - pkin(5) * (t631 * t840 - t632 * t844), t840 * t530 + t844 * t529 - pkin(5) * (t541 * t840 - t542 * t844); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t812, 0, -t811, 0, t794, t873, t867, pkin(5) * t867, 0, 0, t855, 0, t735, 0, -t931, t930, t921, pkin(5) * t921 - pkin(6) * t882 + t840 * t686, 0, 0, t909, 0, t680, 0, -t948, t947, t941, pkin(5) * t941 + t840 * t568 - t844 * t569, 0, 0, t932, 0, t628, 0, -t958, t957, t955, pkin(5) * t955 + t840 * t534 - t844 * t535, t613 * t840 - t615 * t844, t602 * t840 - t603 * t844, t617 * t840 - t619 * t844, t612 * t840 - t614 * t844, t616 * t840 - t618 * t844, t646 * t840 - t647 * t844, -t844 * t538 + t840 * t536 + pkin(5) * (t604 * t844 + t606 * t840), -t844 * t539 + t840 * t537 + pkin(5) * (t605 * t844 + t607 * t840), -t844 * t544 + t840 * t543 + pkin(5) * (t631 * t844 + t632 * t840), -t844 * t530 + t840 * t529 + pkin(5) * (t541 * t844 + t542 * t840); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t847, 0, 0, -g(1), t816, 0, 0, 0, -t800, 0, -t797, 0, t913, t767, t692, pkin(6) * t692, 0, 0, -t729, 0, -t725, 0, t929, t669, t594, t569, 0, 0, -t675, 0, -t671, 0, t946, t601, t550, t535, t615, t603, t619, t614, t618, t647, t538, t539, t544, t530; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t847, 0, qJDD(1), 0, g(1), 0, -t815, 0, 0, 0, t797, 0, -t800, 0, -t767, t913, t870, t686, 0, 0, t725, 0, -t729, 0, -t669, t929, t919, t568, 0, 0, t671, 0, -t675, 0, -t601, t946, t943, t534, t613, t602, t617, t612, t616, t646, t536, t537, t543, t529; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t816, t815, 0, 0, 0, 0, 0, 0, 0, t830, t707, t706, 0, -t687, 0, 0, 0, 0, 0, t823, t645, t644, 0, t579, 0, 0, 0, 0, 0, t823, t597, t596, 0, t540, t737, t714, t747, t736, t745, 0, t570, t571, t562, t532; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t830, 0, -t829, 0, 0, -g(1), t743, 0, 0, 0, -t791, 0, -t788, 0, t912, t762, t640, pkin(7) * t640, 0, 0, -t721, 0, -t717, 0, t928, t661, t566, t547, t657, t643, t665, t656, t664, t695, t558, t559, t561, t533; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t829, 0, t830, 0, g(1), 0, t744, 0, 0, 0, t788, 0, -t791, 0, -t762, t912, t871, t636, 0, 0, t717, 0, -t721, 0, -t661, t928, t918, t546, t655, t642, t663, t654, t662, t694, t553, t554, t560, t531; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t830, -t743, -t744, 0, 0, 0, 0, 0, 0, 0, t823, t859, t850, 0, -t637, 0, 0, 0, 0, 0, t823, t853, t848, 0, t881, t737, t714, t747, t736, t745, 0, t857, t856, t851, t864; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t823, 0, -t822, 0, 0, -g(1), t684, 0, 0, 0, -t785, 0, -t782, 0, t911, t758, t590, qJ(4) * t590, t705, t691, t711, t704, t710, t754, t577, t578, t581, t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t822, 0, t823, 0, g(1), 0, t685, 0, 0, 0, t782, 0, -t785, 0, -t758, t911, t872, t586, t703, t688, t709, t702, t708, t753, t575, t576, t580, t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t823, -t684, -t685, 0, 0, 0, 0, 0, 0, 0, t823, t860, t849, 0, -t587, t737, t714, t747, t736, t745, 0, t861, t862, t858, t875; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t823, 0, -t822, 0, 0, -t834, t634, 0, t742, t715, t751, t741, t749, t771, t610, t611, t584, pkin(8) * t584; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t822, 0, t823, 0, t834, 0, t635, 0, t810, -t793, -t889, -t810, -t817, -qJDD(5), t608, t609, 0, pkin(4) * t584; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t823, -t634, -t635, 0, 0, t737, t714, t747, t736, t745, 0, t876, t877, t863, t880; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t778, t780, t801, -t814, t808, t814, 0, t629, t622, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t874, t777, t806, t779, t802, -t874, -t629, 0, t623, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t810, t793, t889, t810, t817, qJDD(5), -t622, -t623, 0, 0;];
m_new_reg = t1;