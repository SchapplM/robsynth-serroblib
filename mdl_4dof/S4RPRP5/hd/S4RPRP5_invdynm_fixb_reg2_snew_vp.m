% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RPRP5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:09
% EndTime: 2019-12-31 16:45:14
% DurationCPUTime: 5.66s
% Computational Cost: add. (10034->372), mult. (25098->442), div. (0->0), fcn. (16550->6), ass. (0->259)
t778 = sin(pkin(6));
t779 = cos(pkin(6));
t781 = sin(qJ(3));
t783 = cos(qJ(3));
t800 = t778 * t783 + t779 * t781;
t746 = t800 * qJD(1);
t741 = t746 ^ 2;
t785 = qJD(3) ^ 2;
t698 = t785 + t741;
t744 = (t778 * t781 - t779 * t783) * qJD(1);
t849 = t746 * t744;
t864 = qJDD(3) + t849;
t841 = t781 * t864;
t642 = t783 * t698 + t841;
t832 = t783 * t864;
t663 = t781 * t698 - t832;
t593 = t779 * t642 - t778 * t663;
t919 = pkin(1) * t593;
t918 = qJ(2) * t593;
t617 = t778 * t642 + t779 * t663;
t917 = qJ(2) * t617;
t782 = sin(qJ(1));
t916 = t782 * t617;
t784 = cos(qJ(1));
t915 = t784 * t617;
t858 = t744 ^ 2;
t725 = t858 - t785;
t652 = t781 * t725 + t832;
t658 = t783 * t725 - t841;
t613 = t778 * t652 - t779 * t658;
t765 = t778 * qJDD(1);
t767 = t779 * qJDD(1);
t742 = t781 * t765 - t783 * t767;
t914 = t782 * t613 - t784 * t742;
t913 = t784 * t613 + t782 * t742;
t912 = pkin(2) * t642;
t911 = pkin(5) * t642;
t910 = pkin(5) * t663;
t726 = t741 - t785;
t865 = qJDD(3) - t849;
t840 = t781 * t865;
t887 = -t783 * t726 + t840;
t689 = t783 * t865;
t888 = t781 * t726 + t689;
t896 = -t778 * t887 + t779 * t888;
t909 = t782 * t896;
t908 = t784 * t896;
t606 = t779 * t652 + t778 * t658;
t743 = t800 * qJDD(1);
t875 = -t783 * t742 + t781 * t743;
t876 = -t781 * t742 - t783 * t743;
t882 = t778 * t875 + t779 * t876;
t907 = pkin(1) * t882;
t866 = -t858 - t785;
t871 = t783 * t866 - t840;
t877 = t781 * t866 + t689;
t884 = t778 * t871 + t779 * t877;
t906 = pkin(1) * t884;
t905 = qJ(2) * t882;
t904 = qJ(2) * t884;
t826 = t746 * qJD(3);
t705 = t742 + 0.2e1 * t826;
t885 = -t778 * t877 + t779 * t871;
t903 = -pkin(1) * t705 + qJ(2) * t885;
t674 = -t858 - t741;
t883 = -t778 * t876 + t779 * t875;
t902 = -pkin(1) * t674 + qJ(2) * t883;
t901 = pkin(4) * (-t784 * t705 + t782 * t885);
t900 = pkin(4) * (-t784 * t674 + t782 * t883);
t899 = pkin(4) * (t782 * t705 + t784 * t885);
t898 = pkin(4) * (t782 * t674 + t784 * t883);
t897 = t778 * t888 + t779 * t887;
t856 = pkin(2) * t876;
t894 = pkin(2) * t877;
t893 = pkin(5) * t871;
t892 = pkin(5) * t876;
t891 = pkin(5) * t877;
t886 = -pkin(2) * t674 + pkin(5) * t875;
t868 = t741 - t858;
t881 = t782 * t868;
t880 = t784 * t868;
t735 = qJD(3) * t744;
t708 = t743 - t735;
t869 = -t735 + t708;
t879 = t869 * qJ(4);
t758 = t784 * g(1) + t782 * g(2);
t786 = qJD(1) ^ 2;
t812 = -t786 * pkin(1) + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t758;
t878 = pkin(5) * qJDD(1) + t812;
t795 = (-t744 * t781 - t746 * t783) * qJD(3);
t723 = t781 * t826;
t818 = t783 * t735;
t803 = t723 - t818;
t861 = -t778 * t795 + t779 * t803;
t874 = t782 * qJDD(3) + t784 * t861;
t819 = t784 * t849;
t706 = t742 + t826;
t798 = t781 * t706 + t818;
t804 = -t783 * t706 + t781 * t735;
t860 = -t778 * t804 + t779 * t798;
t873 = t782 * t860 + t819;
t872 = -t784 * qJDD(3) + t782 * t861;
t820 = t782 * t849;
t870 = t784 * t860 - t820;
t773 = t778 ^ 2;
t774 = t779 ^ 2;
t867 = t773 + t774;
t863 = t786 * t867;
t862 = t778 * t803 + t779 * t795;
t859 = t778 * t798 + t779 * t804;
t844 = t779 * t786;
t852 = t779 * g(3);
t679 = -t852 + (pkin(2) * t844 - t878) * t778;
t768 = t774 * t786;
t853 = t778 * g(3);
t682 = -pkin(2) * t768 + t878 * t779 - t853;
t629 = -t783 * t679 + t781 * t682;
t630 = t781 * t679 + t783 * t682;
t586 = -t783 * t629 + t781 * t630;
t857 = pkin(2) * t586;
t855 = pkin(3) * t783;
t854 = t706 * pkin(3);
t850 = qJDD(1) * pkin(1);
t848 = t773 * t786;
t847 = t778 * t586;
t846 = t778 * t779;
t845 = t779 * t586;
t843 = t781 * t869;
t757 = t782 * g(1) - t784 * g(2);
t805 = -qJDD(2) + t757;
t699 = (pkin(2) * t779 + pkin(1)) * qJDD(1) + (t867 * pkin(5) + qJ(2)) * t786 + t805;
t842 = t781 * t699;
t839 = t781 * t705;
t737 = t786 * qJ(2) + t805 + t850;
t835 = t782 * t737;
t833 = t783 * t699;
t831 = t783 * t705;
t830 = t784 * t737;
t824 = qJD(4) * qJD(3);
t770 = 0.2e1 * t824;
t695 = t744 * pkin(3) - t746 * qJ(4);
t802 = -t785 * pkin(3) + qJDD(3) * qJ(4) - t744 * t695 + t630;
t601 = t770 + t802;
t603 = -qJDD(3) * pkin(3) - t785 * qJ(4) + t746 * t695 + qJDD(4) + t629;
t828 = -pkin(3) * t603 + qJ(4) * t601;
t827 = -pkin(3) * t743 - qJ(4) * t742;
t823 = t782 * qJDD(1);
t822 = t784 * qJDD(1);
t817 = t778 * t767;
t816 = -qJ(4) * t781 - pkin(2);
t815 = t737 + t850;
t587 = t781 * t629 + t783 * t630;
t715 = t812 * t778 + t852;
t716 = t812 * t779 - t853;
t651 = t778 * t715 + t779 * t716;
t814 = -t782 * t757 - t784 * t758;
t575 = t781 * t601 - t783 * t603;
t811 = pkin(2) * t575 + t828;
t810 = -t630 - t912;
t809 = t827 + t856;
t756 = -t782 * t786 + t822;
t808 = -pkin(4) * t756 - t782 * g(3);
t668 = t781 * t708 + t783 * t826;
t669 = t783 * t708 - t723;
t625 = -t778 * t668 + t779 * t669;
t807 = t782 * t625 - t819;
t806 = t784 * t625 + t820;
t650 = t779 * t715 - t778 * t716;
t801 = t784 * t757 - t782 * t758;
t755 = t784 * t786 + t823;
t799 = -t629 + t894;
t749 = t779 * t863;
t797 = -t782 * t749 + t779 * t822;
t796 = t784 * t749 + t779 * t823;
t794 = pkin(3) * t698 + qJ(4) * t864 + t802;
t793 = t794 + t912;
t792 = pkin(3) * t865 + qJ(4) * t866 - t603;
t791 = t792 + t894;
t790 = -pkin(3) * t826 + 0.2e1 * qJD(4) * t746 + t699;
t789 = t790 + t879;
t766 = t774 * qJDD(1);
t764 = t773 * qJDD(1);
t760 = t778 * t844;
t759 = 0.2e1 * t817;
t754 = -t768 + t848;
t753 = t768 + t848;
t752 = t766 - t764;
t751 = t766 + t764;
t748 = t778 * t863;
t738 = -pkin(4) * t755 + t784 * g(3);
t721 = t756 * t846;
t720 = t755 * t846;
t718 = t784 * t748 + t778 * t823;
t717 = t782 * t748 - t778 * t822;
t707 = t743 - 0.2e1 * t735;
t681 = -qJ(2) * t749 + t815 * t779;
t680 = qJ(2) * t748 - t815 * t778;
t672 = t735 + t708;
t648 = -t781 * t707 - t831;
t646 = t783 * t707 - t839;
t638 = pkin(1) * t737 + qJ(2) * t651;
t633 = t831 + t843;
t632 = -t783 * t869 + t839;
t631 = pkin(1) * t753 + qJ(2) * t751 + t651;
t627 = -t833 + t911;
t626 = -t842 - t891;
t622 = t779 * t668 + t778 * t669;
t605 = -pkin(2) * t707 - t842 + t910;
t604 = -pkin(2) * t705 + t833 + t893;
t600 = -t778 * t646 + t779 * t648;
t598 = t779 * t646 + t778 * t648;
t591 = t789 - t854;
t590 = -t809 - t907;
t589 = -t778 * t632 + t779 * t633;
t588 = t779 * t632 + t778 * t633;
t585 = (-t705 - t706) * pkin(3) + t789;
t584 = -qJ(4) * t674 + t603;
t583 = -pkin(3) * t674 + t601;
t582 = t790 - t854 + 0.2e1 * t879;
t581 = pkin(2) * t699 + pkin(5) * t587;
t580 = -t856 - t907;
t579 = -t586 - t892;
t578 = -t810 + t919;
t577 = -qJ(4) * t831 - t781 * t585 - t891;
t576 = t783 * t601 + t781 * t603;
t574 = t587 + t886;
t573 = -t799 - t906;
t572 = t783 * t585 + t816 * t705 + t893;
t571 = -pkin(3) * t843 + t783 * t582 - t911;
t570 = -t778 * t605 + t779 * t627 + t918;
t569 = t779 * t587 - t847;
t568 = t778 * t587 + t845;
t567 = -t791 - t906;
t566 = -t910 + t781 * t582 + (pkin(2) + t855) * t869;
t565 = -pkin(1) * t707 + t779 * t605 + t778 * t627 + t917;
t564 = -t778 * t604 + t779 * t626 - t904;
t563 = -t781 * t583 + t783 * t584 - t892;
t562 = -t793 - 0.2e1 * t824 - t919;
t561 = t779 * t604 + t778 * t626 + t903;
t560 = t783 * t583 + t781 * t584 + t886;
t559 = -pkin(1) * t568 - t857;
t558 = -t778 * t575 + t779 * t576;
t557 = t779 * t575 + t778 * t576;
t556 = -pkin(5) * t575 + (-pkin(3) * t781 + qJ(4) * t783) * t591;
t555 = -t778 * t574 + t779 * t579 - t905;
t554 = -t778 * t572 + t779 * t577 - t904;
t553 = t779 * t572 + t778 * t577 + t903;
t552 = pkin(5) * t576 + (-t816 + t855) * t591;
t551 = t779 * t574 + t778 * t579 + t902;
t550 = -pkin(5) * t845 - qJ(2) * t568 - t778 * t581;
t549 = pkin(1) * t699 - pkin(5) * t847 + qJ(2) * t569 + t779 * t581;
t548 = -t778 * t566 + t779 * t571 - t918;
t547 = pkin(1) * t869 + t779 * t566 + t778 * t571 - t917;
t546 = -t778 * t560 + t779 * t563 - t905;
t545 = t779 * t560 + t778 * t563 + t902;
t544 = -pkin(1) * t557 - t811;
t543 = -qJ(2) * t557 - t778 * t552 + t779 * t556;
t542 = pkin(1) * t591 + qJ(2) * t558 + t779 * t552 + t778 * t556;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t756, 0, -t755, 0, t808, -t738, -t801, -pkin(4) * t801, t721, t784 * t752 + t782 * t754, t718, -t721, t796, 0, -pkin(4) * t797 - t782 * t715 - t778 * t830, -pkin(4) * t717 - t782 * t716 - t779 * t830, t784 * t650 - pkin(4) * (t782 * t751 + t784 * t753), -pkin(4) * (t782 * t651 + t830) - (t782 * pkin(1) - t784 * qJ(2)) * t650, t806, t784 * t600 + t881, t782 * t743 + t908, t870, -t913, t874, t784 * t564 - t782 * t573 - t901, t784 * t570 - t782 * t578 - pkin(4) * (-t784 * t707 + t916), t784 * t555 - t782 * t580 - t900, t784 * t550 - t782 * t559 - pkin(4) * (t782 * t569 + t784 * t699), t806, t782 * t672 + t908, t784 * t589 - t881, t874, t913, t870, t784 * t554 - t782 * t567 - t901, t784 * t546 - t782 * t590 - t900, t784 * t548 - t782 * t562 - pkin(4) * (t784 * t869 - t916), t784 * t543 - t782 * t544 - pkin(4) * (t782 * t558 + t784 * t591); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t755, 0, t756, 0, t738, t808, t814, pkin(4) * t814, t720, t782 * t752 - t784 * t754, t717, -t720, -t797, 0, -pkin(4) * t796 + t784 * t715 - t778 * t835, pkin(4) * t718 + t784 * t716 - t779 * t835, t782 * t650 + pkin(4) * (t784 * t751 - t782 * t753), pkin(4) * (t784 * t651 - t835) - (-t784 * pkin(1) - t782 * qJ(2)) * t650, t807, t782 * t600 - t880, -t784 * t743 + t909, t873, -t914, t872, t782 * t564 + t784 * t573 + t899, t782 * t570 + t784 * t578 + pkin(4) * (t782 * t707 + t915), t782 * t555 + t784 * t580 + t898, t782 * t550 + t784 * t559 + pkin(4) * (t784 * t569 - t782 * t699), t807, -t784 * t672 + t909, t782 * t589 + t880, t872, t914, t873, t782 * t554 + t784 * t567 + t899, t782 * t546 + t784 * t590 + t898, t782 * t548 + t784 * t562 + pkin(4) * (-t782 * t869 - t915), t782 * t543 + t784 * t544 + pkin(4) * (t784 * t558 - t782 * t591); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t757, t758, 0, 0, t764, t759, 0, t766, 0, 0, t681, t680, t631, t638, t622, t598, t897, t859, t606, t862, t561, t565, t551, t549, t622, t897, t588, t862, -t606, t859, t553, t545, t547, t542; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t786, 0, 0, -g(3), -t757, 0, t817, t752, t748, -t817, t749, 0, -t778 * t737, -t779 * t737, t650, qJ(2) * t650, t625, t600, t896, t860, -t613, t861, t564, t570, t555, t550, t625, t896, t589, t861, t613, t860, t554, t546, t548, t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t786, 0, qJDD(1), 0, g(3), 0, -t758, 0, t760, -t754, -t765, -t760, -t767, 0, t715, t716, 0, pkin(1) * t650, -t849, -t868, -t743, t849, t742, -qJDD(3), t573, t578, t580, t559, -t849, -t672, t868, -qJDD(3), -t742, t849, t567, t590, t562, t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t757, t758, 0, 0, t764, t759, 0, t766, 0, 0, t681, t680, t631, t638, t622, t598, t897, t859, t606, t862, t561, t565, t551, t549, t622, t897, t588, t862, -t606, t859, t553, t545, t547, t542; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t765, t767, t760, 0, t768, 0, 0, -t737, t715, 0, t669, t648, t888, t798, t658, t803, t626, t627, t579, -pkin(5) * t586, t669, t888, t633, t803, -t658, t798, t577, t563, t571, t556; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t765, -t848, t767, -t760, 0, t737, 0, t716, 0, t668, t646, t887, t804, t652, t795, t604, t605, t574, t581, t668, t887, t632, t795, -t652, t804, t572, t560, t566, t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t760, t754, t765, t760, t767, 0, -t715, -t716, 0, 0, t849, t868, t743, -t849, -t742, qJDD(3), t799, t810, t856, t857, t849, t672, -t868, qJDD(3), t742, -t849, t791, t809, t770 + t793, t811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t708, -t705, t865, t735, t725, -t735, 0, -t699, t629, 0, t708, t865, t705, -t735, -t725, t735, -qJ(4) * t705, t584, t582, qJ(4) * t591; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t826, t707, -t726, -t706, t864, -t826, t699, 0, t630, 0, t826, -t726, -t869, -t826, -t864, -t706, t585, t583, pkin(3) * t869, pkin(3) * t591; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t849, t868, t743, -t849, -t742, qJDD(3), -t629, -t630, 0, 0, t849, t672, -t868, qJDD(3), t742, -t849, t792, t827, t770 + t794, t828; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t708, t865, t705, -t735, -t725, t735, 0, t603, t591, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t849, t672, -t868, qJDD(3), t742, -t849, -t603, 0, t601, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t826, t726, t869, t826, t864, t706, -t591, -t601, 0, 0;];
m_new_reg = t1;
