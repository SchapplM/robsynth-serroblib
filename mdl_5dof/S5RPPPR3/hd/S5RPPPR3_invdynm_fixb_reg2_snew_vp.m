% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPPPR3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPPPR3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:10
% EndTime: 2019-12-31 17:44:18
% DurationCPUTime: 7.89s
% Computational Cost: add. (19525->498), mult. (42889->612), div. (0->0), fcn. (26061->8), ass. (0->323)
t808 = sin(pkin(8));
t803 = t808 ^ 2;
t810 = cos(pkin(8));
t804 = t810 ^ 2;
t817 = qJD(1) ^ 2;
t897 = (t803 + t804) * t817;
t757 = t808 * t897;
t809 = sin(pkin(7));
t811 = cos(pkin(7));
t858 = t811 * qJDD(1);
t723 = -t809 * t757 + t808 * t858;
t859 = t809 * qJDD(1);
t726 = t811 * t757 + t808 * t859;
t813 = sin(qJ(1));
t815 = cos(qJ(1));
t906 = t723 * t815 - t813 * t726;
t920 = pkin(5) * t906;
t812 = sin(qJ(5));
t814 = cos(qJ(5));
t827 = t808 * t812 + t810 * t814;
t748 = t827 * qJD(1);
t863 = qJD(1) * t810;
t864 = qJD(1) * t808;
t750 = -t812 * t863 + t814 * t864;
t877 = t750 * t748;
t915 = qJDD(5) - t877;
t919 = t812 * t915;
t918 = t814 * t915;
t776 = t815 * g(1) + t813 * g(2);
t764 = -t817 * pkin(1) - t776;
t775 = t813 * g(1) - t815 * g(2);
t825 = qJDD(1) * pkin(1) + t775;
t708 = t809 * t764 - t811 * t825;
t709 = t811 * t764 + t809 * t825;
t846 = t809 * t708 + t811 * t709;
t647 = t811 * t708 - t809 * t709;
t879 = t647 * t813;
t917 = t815 * t846 + t879;
t878 = t647 * t815;
t916 = -t813 * t846 + t878;
t914 = pkin(1) * t723;
t665 = t723 * t813 + t815 * t726;
t913 = pkin(5) * t665;
t912 = qJ(2) * t723;
t911 = qJ(2) * t726;
t767 = t811 * t817 + t859;
t806 = g(3) - qJDD(2);
t735 = qJ(2) * t767 - t811 * t806;
t768 = -t809 * t817 + t858;
t829 = -qJ(2) * t768 - t809 * t806;
t900 = t815 * t767 + t813 * t768;
t910 = pkin(5) * t900 + t815 * t735 - t813 * t829;
t718 = -t813 * t767 + t815 * t768;
t909 = -pkin(5) * t718 + t813 * t735 + t815 * t829;
t796 = t803 * qJDD(1);
t798 = t804 * qJDD(1);
t766 = t798 - t796;
t801 = t804 * t817;
t876 = t803 * t817;
t771 = -t801 + t876;
t713 = t809 * t766 - t811 * t771;
t716 = t811 * t766 + t809 * t771;
t908 = t713 * t815 + t716 * t813;
t907 = t713 * t813 - t716 * t815;
t758 = t810 * t897;
t848 = t810 * t858;
t725 = -t809 * t758 + t848;
t799 = t810 * qJDD(1);
t728 = t811 * t758 + t799 * t809;
t664 = t815 * t725 - t728 * t813;
t905 = t813 * t725 + t728 * t815;
t694 = -t817 * pkin(2) + qJDD(1) * qJ(3) + t709;
t833 = -pkin(3) * t810 - qJ(4) * t808;
t761 = t833 * qJD(1);
t826 = t694 + ((2 * qJD(3)) + t761) * qJD(1);
t789 = t810 * t806;
t857 = qJDD(4) + t789;
t874 = t810 * t817;
t625 = (-pkin(4) * t874 - pkin(6) * qJDD(1) + t826) * t808 + t857;
t860 = qJD(1) * qJD(3);
t850 = t810 * t860;
t786 = 0.2e1 * t850;
t868 = t810 * t694 - t808 * t806;
t851 = -t761 * t863 - t868;
t643 = t786 - t851;
t631 = -pkin(4) * t801 - pkin(6) * t799 + t643;
t564 = -t814 * t625 + t812 * t631;
t565 = t812 * t625 + t814 * t631;
t535 = -t814 * t564 + t565 * t812;
t676 = t789 + (t694 + 0.2e1 * t860) * t808;
t677 = t786 + t868;
t614 = t808 * t676 + t810 * t677;
t797 = t808 * qJDD(1);
t899 = -pkin(2) * t797 + qJ(3) * t757;
t898 = t827 * qJDD(1);
t536 = t812 * t564 + t814 * t565;
t891 = pkin(3) + pkin(4);
t896 = qJ(4) * t536 - t891 * t535;
t744 = t748 ^ 2;
t816 = qJD(5) ^ 2;
t698 = -t816 - t744;
t635 = t812 * t698 + t918;
t636 = t814 * t698 - t919;
t895 = qJ(4) * t636 - t891 * t635 + t564;
t747 = t797 * t814 - t812 * t799;
t638 = -t814 * t747 - t812 * t898;
t640 = t812 * t747 - t814 * t898;
t894 = qJ(4) * t640 - t891 * t638;
t745 = t750 ^ 2;
t738 = -t745 - t816;
t700 = qJDD(5) + t877;
t872 = t812 * t700;
t657 = t814 * t738 - t872;
t869 = t814 * t700;
t660 = -t812 * t738 - t869;
t893 = qJ(4) * t660 - t891 * t657 + t565;
t765 = t798 + t796;
t770 = t801 + t876;
t712 = t809 * t765 + t811 * t770;
t715 = t811 * t765 - t809 * t770;
t890 = pkin(5) * (t815 * t712 + t813 * t715);
t889 = pkin(5) * t664;
t887 = pkin(6) * t535;
t886 = pkin(6) * t536;
t885 = qJ(2) * t712;
t884 = qJ(2) * t725;
t693 = -qJDD(1) * pkin(2) - t817 * qJ(3) + qJDD(3) + t708;
t687 = t808 * t693;
t875 = t808 * t810;
t688 = t810 * t693;
t790 = qJ(4) * t797;
t793 = pkin(3) * t799;
t821 = t693 - t793;
t849 = qJD(4) * t864;
t669 = -t790 + t821 - 0.2e1 * t849;
t637 = -pkin(4) * t799 + pkin(6) * t897 + t669;
t873 = t812 * t637;
t632 = t814 * t637;
t867 = -pkin(2) * t693 + qJ(3) * t614;
t866 = pkin(2) * t770 + qJ(3) * t765;
t865 = pkin(2) * t799 - qJ(3) * t758;
t862 = t748 * qJD(5);
t861 = t750 * qJD(5);
t856 = pkin(3) * t797;
t855 = t809 * t877;
t854 = t811 * t877;
t781 = t808 * t874;
t853 = t687 + t899;
t852 = -t688 + t865;
t780 = t808 * t799;
t847 = -pkin(6) * t657 - t632;
t845 = -t775 * t813 - t815 * t776;
t684 = -t744 - t745;
t823 = -pkin(6) * t640 - t536;
t523 = t891 * t684 + t823;
t824 = -pkin(6) * t638 - t535;
t530 = qJ(4) * t684 + t824;
t583 = t808 * t638 + t810 * t640;
t844 = pkin(2) * t684 + qJ(3) * t583 + t810 * t523 + t808 * t530;
t517 = t535 * t808 + t536 * t810;
t527 = -t891 * t637 - t886;
t532 = -qJ(4) * t637 - t887;
t843 = -pkin(2) * t637 + qJ(3) * t517 + t810 * t527 + t808 * t532;
t702 = t898 + 0.2e1 * t861;
t831 = -pkin(6) * t636 - t632;
t544 = t891 * t702 + t831;
t832 = -pkin(6) * t635 - t873;
t558 = qJ(4) * t702 + t832;
t579 = t808 * t635 + t810 * t636;
t842 = pkin(2) * t702 + qJ(3) * t579 + t810 * t544 + t808 * t558;
t704 = -0.2e1 * t862 + t747;
t830 = -pkin(6) * t660 + t873;
t550 = t891 * t704 + t830;
t567 = qJ(4) * t704 + t847;
t604 = t808 * t657 + t810 * t660;
t841 = pkin(2) * t704 + qJ(3) * t604 + t810 * t550 + t808 * t567;
t642 = t808 * t826 + t857;
t629 = qJ(4) * t770 + t642;
t630 = pkin(3) * t770 + t643;
t840 = t808 * t629 + t810 * t630 + t866;
t784 = 0.2e1 * t849;
t670 = t784 + 0.2e1 * t790 - t821;
t839 = pkin(3) * t780 + t808 * t670 - t899;
t838 = t866 + t614;
t671 = -t693 + t784 + t790 + 0.2e1 * t793;
t837 = qJ(4) * t780 + t810 * t671 + t865;
t773 = t815 * qJDD(1) - t813 * t817;
t836 = -pkin(5) * t773 - g(3) * t813;
t834 = -pkin(3) * t642 + qJ(4) * t643;
t613 = t810 * t676 - t808 * t677;
t730 = t767 * t875;
t731 = -t781 * t809 + t808 * t848;
t678 = t815 * t730 + t813 * t731;
t681 = t813 * t730 - t815 * t731;
t828 = t775 * t815 - t776 * t813;
t587 = t808 * t642 + t810 * t643;
t820 = qJ(3) * t587 + (-pkin(2) + t833) * t669;
t791 = qJ(4) * t799;
t779 = -0.2e1 * t780;
t778 = 0.2e1 * t780;
t772 = t813 * qJDD(1) + t815 * t817;
t759 = -t791 + t856;
t743 = -pkin(5) * t772 + t815 * g(3);
t737 = -t745 + t816;
t736 = t744 - t816;
t721 = pkin(1) * t725;
t720 = qJ(2) * t728;
t711 = pkin(1) * t712;
t710 = qJ(2) * t715;
t706 = t745 - t744;
t705 = -t862 + t747;
t703 = -t898 - t861;
t691 = (-t748 * t814 + t750 * t812) * qJD(5);
t690 = (-t748 * t812 - t750 * t814) * qJD(5);
t686 = -pkin(1) * t767 - t709;
t685 = pkin(1) * t768 - t708;
t675 = t814 * t705 - t812 * t861;
t674 = t812 * t705 + t814 * t861;
t673 = -t812 * t703 + t814 * t862;
t672 = -t814 * t703 - t812 * t862;
t661 = pkin(5) * t905;
t659 = -t812 * t737 + t918;
t658 = t814 * t736 - t872;
t656 = t814 * t737 + t919;
t655 = t812 * t736 + t869;
t651 = pkin(5) * (-t813 * t712 + t815 * t715);
t650 = -pkin(3) * t796 + t810 * t670;
t649 = qJ(4) * t798 - t808 * t671;
t644 = pkin(1) * t647;
t641 = -t814 * t702 - t812 * t704;
t639 = -t812 * t702 + t814 * t704;
t633 = pkin(1) * t806 + qJ(2) * t846;
t628 = -0.2e1 * t850 + (-pkin(3) * t803 + qJ(4) * t875) * t817 + t851;
t624 = qJ(4) * t801 + (-pkin(3) * t874 + t826) * t808 + t857;
t622 = t808 * t690 + t810 * t691;
t621 = -t810 * t690 + t808 * t691;
t618 = t721 + t852;
t617 = t853 - t914;
t616 = -t809 * qJDD(5) + t811 * t622;
t615 = t811 * qJDD(5) + t809 * t622;
t611 = t839 + t914;
t609 = t721 + t837;
t608 = t808 * t674 + t810 * t675;
t607 = -t808 * t672 + t810 * t673;
t606 = -t810 * t674 + t808 * t675;
t605 = t810 * t672 + t808 * t673;
t603 = t808 * t656 + t810 * t659;
t602 = t808 * t655 + t810 * t658;
t601 = -t810 * t657 + t808 * t660;
t600 = -t810 * t656 + t808 * t659;
t599 = -t810 * t655 + t808 * t658;
t597 = -t809 * t677 + t811 * t688 + t912;
t596 = -t809 * t676 + t811 * t687 - t884;
t595 = t811 * t677 + t809 * t688 + t911;
t594 = t811 * t676 + t809 * t687 - t720;
t593 = t811 * t603 - t809 * t747;
t592 = t811 * t602 + t809 * t898;
t591 = t809 * t603 + t811 * t747;
t590 = t809 * t602 - t811 * t898;
t589 = t811 * t613 - t885;
t588 = t809 * t613 + t710;
t586 = -t810 * t642 + t808 * t643;
t584 = t808 * t639 + t810 * t641;
t582 = -t810 * t639 + t808 * t641;
t581 = -t810 * t638 + t808 * t640;
t578 = -t810 * t635 + t808 * t636;
t576 = t608 * t811 - t855;
t575 = t607 * t811 + t855;
t574 = t608 * t809 + t854;
t573 = t607 * t809 - t854;
t572 = t614 * t811 + t693 * t809;
t571 = t614 * t809 - t693 * t811;
t570 = t604 * t811 - t704 * t809;
t569 = t604 * t809 + t704 * t811;
t568 = t711 + t838;
t566 = t629 * t810 - t630 * t808;
t560 = -t628 * t809 + t650 * t811 - t912;
t559 = t628 * t811 + t650 * t809 - t911;
t557 = -t624 * t809 + t649 * t811 - t884;
t556 = t624 * t811 + t649 * t809 - t720;
t554 = t584 * t811 - t706 * t809;
t553 = t584 * t809 + t706 * t811;
t552 = t579 * t811 - t702 * t809;
t551 = t579 * t809 + t702 * t811;
t549 = t583 * t811 - t684 * t809;
t548 = t583 * t809 + t684 * t811;
t546 = t587 * t811 + t669 * t809;
t545 = t587 * t809 - t669 * t811;
t543 = t566 * t811 - t759 * t809 - t885;
t542 = t566 * t809 + t759 * t811 + t710;
t540 = t711 + t840;
t539 = -qJ(3) * t586 + (pkin(3) * t808 - qJ(4) * t810) * t669;
t538 = -pkin(2) * t586 - t834;
t537 = pkin(1) * t571 + t867;
t534 = -pkin(2) * t581 - t894;
t533 = -qJ(2) * t571 - (pkin(2) * t809 - qJ(3) * t811) * t613;
t529 = -pkin(2) * t601 - t893;
t525 = -qJ(3) * t601 - t550 * t808 + t567 * t810;
t524 = qJ(2) * t572 - (-pkin(2) * t811 - qJ(3) * t809 - pkin(1)) * t613;
t521 = pkin(1) * t545 + t820;
t520 = -pkin(2) * t578 - t895;
t519 = -qJ(3) * t578 - t544 * t808 + t558 * t810;
t518 = pkin(1) * t569 + t841;
t516 = -t535 * t810 + t536 * t808;
t514 = -qJ(2) * t545 - t538 * t809 + t539 * t811;
t513 = t517 * t811 + t637 * t809;
t512 = t517 * t809 - t637 * t811;
t511 = pkin(1) * t551 + t842;
t510 = -pkin(1) * t586 + qJ(2) * t546 + t538 * t811 + t539 * t809;
t509 = -qJ(2) * t569 + t525 * t811 - t529 * t809;
t508 = -qJ(3) * t581 - t523 * t808 + t530 * t810;
t507 = -pkin(1) * t601 + qJ(2) * t570 + t525 * t809 + t529 * t811;
t506 = -qJ(2) * t551 + t519 * t811 - t520 * t809;
t505 = pkin(1) * t548 + t844;
t504 = -pkin(1) * t578 + qJ(2) * t552 + t519 * t809 + t520 * t811;
t503 = -qJ(3) * t516 - t527 * t808 + t532 * t810;
t502 = -qJ(2) * t548 + t508 * t811 - t534 * t809;
t501 = -pkin(2) * t516 - t896;
t500 = -pkin(1) * t581 + qJ(2) * t549 + t508 * t809 + t534 * t811;
t499 = pkin(1) * t512 + t843;
t498 = -qJ(2) * t512 - t501 * t809 + t503 * t811;
t497 = -pkin(1) * t516 + qJ(2) * t513 + t501 * t811 + t503 * t809;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t773, 0, -t772, 0, t836, -t743, -t828, -pkin(5) * t828, 0, 0, t718, 0, -t900, 0, t909, t910, t916, pkin(5) * t916 + qJ(2) * t878 - t813 * t633, -t681, -t907, t665, t681, t905, 0, -t594 * t813 + t596 * t815 - t889, -t595 * t813 + t597 * t815 + t920, -t588 * t813 + t589 * t815 - t890, t815 * t533 - t813 * t524 - pkin(5) * (t571 * t815 + t572 * t813), -t681, t665, t907, 0, -t905, t681, -t556 * t813 + t557 * t815 - t889, -t542 * t813 + t543 * t815 - t890, -t813 * t559 + t815 * t560 - t920, t815 * t514 - t813 * t510 - pkin(5) * (t545 * t815 + t546 * t813), -t574 * t813 + t576 * t815, -t553 * t813 + t554 * t815, -t591 * t813 + t593 * t815, -t573 * t813 + t575 * t815, -t590 * t813 + t592 * t815, -t615 * t813 + t616 * t815, t815 * t506 - t813 * t504 - pkin(5) * (t551 * t815 + t552 * t813), t815 * t509 - t813 * t507 - pkin(5) * (t569 * t815 + t570 * t813), t815 * t502 - t813 * t500 - pkin(5) * (t548 * t815 + t549 * t813), t815 * t498 - t813 * t497 - pkin(5) * (t512 * t815 + t513 * t813); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t772, 0, t773, 0, t743, t836, t845, pkin(5) * t845, 0, 0, t900, 0, t718, 0, -t910, t909, t917, pkin(5) * t917 + qJ(2) * t879 + t815 * t633, t678, t908, -t906, -t678, -t664, 0, t594 * t815 + t596 * t813 - t661, t595 * t815 + t597 * t813 + t913, t588 * t815 + t589 * t813 + t651, t813 * t533 + t815 * t524 + pkin(5) * (-t571 * t813 + t572 * t815), t678, -t906, -t908, 0, t664, -t678, t556 * t815 + t557 * t813 - t661, t542 * t815 + t543 * t813 + t651, t815 * t559 + t813 * t560 - t913, t813 * t514 + t815 * t510 + pkin(5) * (-t545 * t813 + t546 * t815), t574 * t815 + t576 * t813, t553 * t815 + t554 * t813, t591 * t815 + t593 * t813, t573 * t815 + t575 * t813, t590 * t815 + t592 * t813, t615 * t815 + t616 * t813, t813 * t506 + t815 * t504 + pkin(5) * (-t551 * t813 + t552 * t815), t813 * t509 + t815 * t507 + pkin(5) * (-t569 * t813 + t570 * t815), t813 * t502 + t815 * t500 + pkin(5) * (-t548 * t813 + t549 * t815), t813 * t498 + t815 * t497 + pkin(5) * (-t512 * t813 + t513 * t815); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t775, t776, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t685, t686, 0, -t644, t796, t778, 0, t798, 0, 0, t618, t617, t568, t537, t796, 0, t779, 0, 0, t798, t609, t540, t611, t521, t606, t582, t600, t605, t599, t621, t511, t518, t505, t499; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t817, 0, 0, -g(3), -t775, 0, 0, 0, t768, 0, -t767, 0, t829, t735, t647, qJ(2) * t647, t731, t716, t726, -t731, t728, 0, t596, t597, t589, t533, t731, t726, -t716, 0, -t728, -t731, t557, t543, t560, t514, t576, t554, t593, t575, t592, t616, t506, t509, t502, t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t817, 0, qJDD(1), 0, g(3), 0, -t776, 0, 0, 0, t767, 0, t768, 0, -t735, t829, t846, t633, t730, t713, -t723, -t730, -t725, 0, t594, t595, t588, t524, t730, -t723, -t713, 0, t725, -t730, t556, t542, t559, t510, t574, t553, t591, t573, t590, t615, t504, t507, t500, t497; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t775, t776, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t685, t686, 0, -t644, t796, t778, 0, t798, 0, 0, t618, t617, t568, t537, t796, 0, t779, 0, 0, t798, t609, t540, t611, t521, t606, t582, t600, t605, t599, t621, t511, t518, t505, t499; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t817, 0, 0, -t806, t708, 0, t780, t766, t757, -t780, t758, 0, t687, t688, t613, qJ(3) * t613, t780, t757, -t766, 0, -t758, -t780, t649, t566, t650, t539, t608, t584, t603, t607, t602, t622, t519, t525, t508, t503; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t817, 0, qJDD(1), 0, t806, 0, t709, 0, t781, -t771, -t797, -t781, -t799, 0, t676, t677, 0, pkin(2) * t613, t781, -t797, t771, 0, t799, -t781, t624, t759, t628, t538, t877, t706, t747, -t877, -t898, qJDD(5), t520, t529, t534, t501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t708, -t709, 0, 0, t796, t778, 0, t798, 0, 0, t852, t853, t838, t867, t796, 0, t779, 0, 0, t798, t837, t840, t839, t820, t606, t582, t600, t605, t599, t621, t842, t841, t844, t843; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t797, t799, t781, 0, t801, 0, 0, t693, t676, 0, t797, t781, -t799, 0, -t801, 0, t791, t629, t670, -qJ(4) * t669, t675, t641, t659, t673, t658, t691, t558, t567, t530, t532; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t797, -t876, t799, -t781, 0, -t693, 0, t677, 0, 0, -t876, -t797, 0, t781, t799, t671, t630, t856, -pkin(3) * t669, -t674, -t639, -t656, t672, -t655, -t690, t544, t550, t523, t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t781, t771, t797, t781, t799, 0, -t676, -t677, 0, 0, -t781, t797, -t771, 0, -t799, t781, -t624, -t759, -t628, t834, -t877, -t706, -t747, t877, t898, -qJDD(5), t895, t893, t894, t896; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t797, t781, -t799, 0, -t801, 0, 0, t642, -t669, 0, t675, t641, t659, t673, t658, t691, t832, t847, t824, -t887; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t781, t797, -t771, 0, -t799, t781, -t642, 0, t643, 0, -t877, -t706, -t747, t877, t898, -qJDD(5), -pkin(4) * t635 + t564, -pkin(4) * t657 + t565, -pkin(4) * t638, -pkin(4) * t535; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t876, t797, 0, -t781, -t799, t669, -t643, 0, 0, t674, t639, t656, -t672, t655, t690, -pkin(4) * t702 - t831, -pkin(4) * t704 - t830, -pkin(4) * t684 - t823, pkin(4) * t637 + t886; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t705, -t702, t915, t862, t736, -t862, 0, -t637, t564, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t861, t704, t737, t703, t700, -t861, t637, 0, t565, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t877, t706, t747, -t877, -t898, qJDD(5), -t564, -t565, 0, 0;];
m_new_reg = t1;
