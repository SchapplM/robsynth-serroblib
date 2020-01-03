% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPP4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:59
% EndTime: 2019-12-31 20:56:17
% DurationCPUTime: 18.82s
% Computational Cost: add. (50006->563), mult. (113797->786), div. (0->0), fcn. (80422->8), ass. (0->416)
t773 = sin(qJ(3));
t776 = cos(qJ(3));
t777 = cos(qJ(2));
t774 = sin(qJ(2));
t819 = qJD(1) * t774;
t728 = -t776 * t777 * qJD(1) + t773 * t819;
t730 = (t773 * t777 + t774 * t776) * qJD(1);
t771 = sin(pkin(8));
t772 = cos(pkin(8));
t695 = -t771 * t728 + t772 * t730;
t691 = t695 ^ 2;
t768 = qJD(2) + qJD(3);
t862 = t768 ^ 2;
t619 = t862 + t691;
t693 = t772 * t728 + t771 * t730;
t631 = t695 * t693;
t767 = qJDD(2) + qJDD(3);
t872 = t631 + t767;
t847 = t771 * t872;
t555 = t772 * t619 + t847;
t842 = t772 * t872;
t588 = t771 * t619 - t842;
t498 = t776 * t555 - t773 * t588;
t527 = t773 * t555 + t776 * t588;
t474 = t774 * t498 + t777 * t527;
t775 = sin(qJ(1));
t778 = cos(qJ(1));
t814 = qJD(1) * qJD(2);
t804 = t777 * t814;
t813 = t774 * qJDD(1);
t739 = t804 + t813;
t762 = t777 * qJDD(1);
t805 = t774 * t814;
t740 = t762 - t805;
t801 = t773 * t739 - t776 * t740;
t661 = -t730 * qJD(3) - t801;
t783 = t728 * qJD(3) - t776 * t739 - t773 * t740;
t611 = t771 * t661 - t772 * t783;
t858 = t768 * t693;
t874 = t611 - t858;
t446 = t775 * t474 - t778 * t874;
t967 = pkin(5) * t446;
t448 = t778 * t474 + t775 * t874;
t966 = pkin(5) * t448;
t450 = t777 * t498 - t774 * t527;
t965 = pkin(6) * t450;
t964 = -pkin(1) * t450 - pkin(2) * t498 - pkin(3) * t555;
t963 = pkin(1) * t874 - pkin(6) * t474;
t802 = -t772 * t661 - t771 * t783;
t857 = t768 * t695;
t786 = t802 + t857;
t514 = -t771 * t786 + t772 * t874;
t849 = t771 * t874;
t516 = t772 * t786 + t849;
t464 = -t776 * t514 + t773 * t516;
t468 = t773 * t514 + t776 * t516;
t425 = t774 * t464 - t777 * t468;
t863 = t693 ^ 2;
t626 = t691 - t863;
t962 = t775 * t425 - t778 * t626;
t961 = t778 * t425 + t775 * t626;
t673 = t863 - t862;
t592 = t771 * t673 + t842;
t596 = t772 * t673 - t847;
t531 = t776 * t592 + t773 * t596;
t536 = t773 * t592 - t776 * t596;
t479 = t774 * t531 + t777 * t536;
t573 = t802 - t857;
t960 = t775 * t479 - t778 * t573;
t959 = t778 * t479 + t775 * t573;
t958 = pkin(7) * t498;
t956 = pkin(2) * t874 - pkin(7) * t527;
t955 = t777 * t464 + t774 * t468;
t954 = t777 * t531 - t774 * t536;
t605 = -t863 - t691;
t875 = -t611 - t858;
t886 = -t772 * t573 - t771 * t875;
t887 = -t771 * t573 + t772 * t875;
t898 = t773 * t886 + t776 * t887;
t899 = -t773 * t887 + t776 * t886;
t922 = -t774 * t898 + t777 * t899;
t936 = t775 * t605 + t778 * t922;
t953 = pkin(5) * t936;
t873 = -t631 + t767;
t846 = t771 * t873;
t868 = -t862 - t863;
t884 = t772 * t868 - t846;
t617 = t772 * t873;
t885 = t771 * t868 + t617;
t900 = t773 * t884 + t776 * t885;
t901 = -t773 * t885 + t776 * t884;
t920 = -t774 * t900 + t777 * t901;
t937 = t775 * t786 + t778 * t920;
t952 = pkin(5) * t937;
t938 = -t778 * t605 + t775 * t922;
t951 = pkin(5) * t938;
t939 = t775 * t920 - t778 * t786;
t950 = pkin(5) * t939;
t674 = -t691 + t862;
t907 = t772 * t674 + t846;
t908 = -t771 * t674 + t617;
t924 = -t773 * t907 + t776 * t908;
t925 = t773 * t908 + t776 * t907;
t935 = -t774 * t925 + t777 * t924;
t949 = t775 * t935 + t778 * t875;
t948 = -t775 * t875 + t778 * t935;
t921 = t774 * t901 + t777 * t900;
t946 = pkin(6) * t921;
t923 = t774 * t899 + t777 * t898;
t945 = pkin(6) * t923;
t944 = qJ(4) * t555;
t943 = qJ(4) * t588;
t942 = -pkin(1) * t921 - pkin(2) * t900 - pkin(3) * t885;
t403 = -pkin(1) * t923 - pkin(2) * t898 - pkin(3) * t887;
t941 = -pkin(1) * t786 + pkin(6) * t920;
t940 = -pkin(1) * t605 + pkin(6) * t922;
t934 = t774 * t924 + t777 * t925;
t933 = pkin(7) * t898;
t932 = pkin(7) * t900;
t927 = -pkin(2) * t786 + pkin(7) * t901;
t926 = -pkin(2) * t605 + pkin(7) * t899;
t915 = qJ(4) * t884;
t914 = qJ(4) * t885;
t913 = qJ(4) * t887;
t758 = t778 * t767;
t784 = (-t693 * t771 - t695 * t772) * t768;
t856 = t768 * t771;
t668 = t695 * t856;
t855 = t768 * t772;
t810 = t693 * t855;
t792 = t668 - t810;
t866 = -t773 * t784 + t776 * t792;
t867 = t773 * t792 + t776 * t784;
t881 = -t774 * t867 + t777 * t866;
t906 = t775 * t881 - t758;
t807 = t778 * t631;
t787 = t771 * t802 + t810;
t793 = t693 * t856 - t772 * t802;
t864 = t773 * t787 + t776 * t793;
t865 = -t773 * t793 + t776 * t787;
t882 = -t774 * t864 + t777 * t865;
t905 = t775 * t882 + t807;
t809 = t775 * t631;
t904 = t778 * t882 - t809;
t831 = t775 * t767;
t903 = t778 * t881 + t831;
t902 = -pkin(3) * t605 + qJ(4) * t886;
t897 = 2 * qJD(5);
t701 = t730 * t728;
t871 = -t701 + t767;
t894 = t773 * t871;
t891 = t776 * t871;
t888 = t874 * qJ(5);
t722 = t768 * t728;
t643 = -t722 + t783;
t883 = t774 * t865 + t777 * t864;
t880 = t774 * t866 + t777 * t867;
t870 = -t722 - t783;
t625 = t693 * pkin(4) - t695 * qJ(5);
t780 = qJD(1) ^ 2;
t833 = t774 * t780;
t750 = t778 * g(1) + t775 * g(2);
t732 = -t780 * pkin(1) + qJDD(1) * pkin(6) - t750;
t836 = t774 * t732;
t653 = qJDD(2) * pkin(2) - t739 * pkin(7) - t836 + (pkin(2) * t833 + pkin(7) * t814 - g(3)) * t777;
t716 = -t774 * g(3) + t777 * t732;
t770 = t777 ^ 2;
t764 = t770 * t780;
t790 = qJD(2) * pkin(2) - pkin(7) * t819;
t657 = -pkin(2) * t764 + t740 * pkin(7) - qJD(2) * t790 + t716;
t607 = -t776 * t653 + t773 * t657;
t543 = t871 * pkin(3) + t643 * qJ(4) - t607;
t608 = t773 * t653 + t776 * t657;
t726 = t728 ^ 2;
t791 = t768 * pkin(3) - t730 * qJ(4);
t548 = -t726 * pkin(3) + t661 * qJ(4) - t768 * t791 + t608;
t823 = t771 * t543 + t772 * t548;
t869 = t767 * qJ(5) - t693 * t625 + t768 * t897 + t823;
t639 = (qJD(3) - t768) * t730 + t801;
t727 = t730 ^ 2;
t861 = pkin(4) * t772;
t860 = t802 * pkin(4);
t859 = qJ(5) * t772;
t854 = t768 * t773;
t853 = t768 * t776;
t769 = t774 ^ 2;
t852 = t769 * t780;
t749 = t775 * g(1) - t778 * g(2);
t789 = qJDD(1) * pkin(1) + t749;
t664 = t740 * pkin(2) - t790 * t819 + (pkin(7) * t770 + pkin(6)) * t780 + t789;
t563 = t661 * pkin(3) + t726 * qJ(4) - t730 * t791 - qJDD(4) + t664;
t851 = t771 * t563;
t844 = t772 * t563;
t817 = qJD(4) * t695;
t688 = 0.2e1 * t817;
t822 = -t772 * t543 + t771 * t548;
t488 = t688 + t822;
t818 = qJD(4) * t693;
t686 = -0.2e1 * t818;
t489 = t686 + t823;
t436 = -t772 * t488 + t771 * t489;
t841 = t773 * t436;
t840 = t773 * t664;
t698 = t701 + t767;
t839 = t773 * t698;
t539 = -t776 * t607 + t773 * t608;
t838 = t774 * t539;
t731 = t780 * pkin(6) + t789;
t837 = t774 * t731;
t757 = t777 * t833;
t747 = qJDD(2) + t757;
t835 = t774 * t747;
t748 = qJDD(2) - t757;
t834 = t774 * t748;
t830 = t776 * t436;
t829 = t776 * t664;
t828 = t776 * t698;
t827 = t777 * t539;
t826 = t777 * t731;
t825 = t777 * t748;
t821 = -t605 - t862;
t820 = t769 + t770;
t812 = t775 * qJDD(1);
t811 = t778 * qJDD(1);
t808 = t775 * t701;
t806 = t778 * t701;
t803 = -qJ(5) * t771 - pkin(3);
t437 = t771 * t488 + t772 * t489;
t540 = t773 * t607 + t776 * t608;
t715 = t777 * g(3) + t836;
t656 = t774 * t715 + t777 * t716;
t707 = -t775 * t749 - t778 * t750;
t800 = t775 * t757;
t799 = t778 * t757;
t797 = t695 * t625 + qJDD(5) + t822;
t744 = -t775 * t780 + t811;
t796 = -pkin(5) * t744 - t775 * g(3);
t568 = t771 * t611 + t695 * t855;
t569 = t772 * t611 - t668;
t508 = t776 * t568 + t773 * t569;
t511 = -t773 * t568 + t776 * t569;
t461 = -t774 * t508 + t777 * t511;
t795 = t775 * t461 - t807;
t794 = t778 * t461 + t809;
t655 = t777 * t715 - t774 * t716;
t706 = t778 * t749 - t775 * t750;
t788 = t686 + t869;
t785 = -t767 * pkin(4) + t797;
t782 = -pkin(4) * t857 + t695 * t897 + t563;
t781 = t782 + t888;
t779 = qJD(2) ^ 2;
t755 = -t764 - t779;
t754 = t764 - t779;
t753 = -t779 - t852;
t752 = t779 - t852;
t746 = t764 - t852;
t745 = t764 + t852;
t743 = t778 * t780 + t812;
t742 = t820 * qJDD(1);
t741 = t762 - 0.2e1 * t805;
t738 = 0.2e1 * t804 + t813;
t736 = t777 * t747;
t735 = t820 * t814;
t725 = -pkin(5) * t743 + t778 * g(3);
t720 = -t727 + t862;
t719 = t726 - t862;
t718 = t777 * t739 - t769 * t814;
t717 = -t774 * t740 - t770 * t814;
t714 = -t727 - t862;
t713 = -t774 * t753 - t825;
t712 = -t774 * t752 + t736;
t711 = t777 * t755 - t835;
t710 = t777 * t754 - t834;
t709 = t777 * t753 - t834;
t708 = t774 * t755 + t736;
t704 = t778 * t742 - t775 * t745;
t703 = t775 * t742 + t778 * t745;
t702 = -t774 * t738 + t777 * t741;
t700 = -t727 + t726;
t696 = -t862 - t726;
t680 = t778 * t713 + t775 * t738;
t679 = t778 * t711 - t775 * t741;
t678 = t775 * t713 - t778 * t738;
t677 = t775 * t711 + t778 * t741;
t672 = (-t728 * t776 + t730 * t773) * t768;
t671 = (-t728 * t773 - t730 * t776) * t768;
t670 = -pkin(6) * t709 - t826;
t669 = -pkin(6) * t708 - t837;
t663 = -t726 - t727;
t660 = -pkin(1) * t709 + t716;
t659 = -pkin(1) * t708 + t715;
t649 = t776 * t719 - t839;
t648 = -t773 * t720 + t891;
t647 = t773 * t719 + t828;
t646 = t776 * t720 + t894;
t645 = -t773 * t714 - t828;
t644 = t776 * t714 - t839;
t638 = (qJD(3) + t768) * t730 + t801;
t637 = -t730 * t854 - t776 * t783;
t636 = t730 * t853 - t773 * t783;
t635 = -t773 * t661 + t728 * t853;
t634 = t776 * t661 + t728 * t854;
t633 = t778 * t656 - t775 * t731;
t632 = t775 * t656 + t778 * t731;
t629 = t776 * t696 - t894;
t628 = t773 * t696 + t891;
t612 = -t774 * t671 + t777 * t672;
t603 = -pkin(7) * t644 - t829;
t600 = -pkin(7) * t628 - t840;
t599 = -t774 * t647 + t777 * t649;
t598 = -t774 * t646 + t777 * t648;
t585 = -t774 * t644 + t777 * t645;
t584 = t777 * t644 + t774 * t645;
t583 = -t639 * t776 - t773 * t643;
t582 = -t776 * t638 - t773 * t870;
t581 = -t639 * t773 + t776 * t643;
t580 = -t773 * t638 + t776 * t870;
t562 = -t774 * t636 + t777 * t637;
t561 = -t774 * t634 + t777 * t635;
t560 = -t774 * t628 + t777 * t629;
t559 = t777 * t628 + t774 * t629;
t554 = -pkin(2) * t870 + pkin(7) * t645 - t840;
t549 = -pkin(2) * t638 + pkin(7) * t629 + t829;
t545 = t778 * t585 + t775 * t870;
t544 = t775 * t585 - t778 * t870;
t538 = t778 * t560 + t775 * t638;
t537 = t775 * t560 - t778 * t638;
t524 = pkin(2) * t664 + pkin(7) * t540;
t523 = -t774 * t581 + t777 * t583;
t522 = -t774 * t580 + t777 * t582;
t521 = t777 * t581 + t774 * t583;
t520 = -t844 + t944;
t505 = -pkin(1) * t584 - pkin(2) * t644 + t608;
t504 = -t851 - t914;
t503 = t778 * t523 + t775 * t663;
t502 = t775 * t523 - t778 * t663;
t497 = -pkin(1) * t559 - pkin(2) * t628 + t607;
t496 = -pkin(7) * t581 - t539;
t493 = -pkin(2) * t663 + pkin(7) * t583 + t540;
t492 = -pkin(1) * t521 - pkin(2) * t581;
t491 = -pkin(3) * t874 - t851 + t943;
t490 = -pkin(6) * t584 - t774 * t554 + t777 * t603;
t487 = -pkin(3) * t786 + t844 + t915;
t486 = t777 * t540 - t838;
t485 = t774 * t540 + t827;
t483 = t781 - t860;
t482 = -pkin(6) * t559 - t774 * t549 + t777 * t600;
t481 = t778 * t486 - t775 * t664;
t480 = t775 * t486 + t778 * t664;
t463 = qJ(5) * t862 - t785 - 0.2e1 * t817;
t462 = -pkin(4) * t862 + t788;
t458 = t777 * t508 + t774 * t511;
t455 = (-t786 - t802) * pkin(4) + t781;
t454 = t782 - t860 + 0.2e1 * t888;
t445 = -pkin(1) * t485 - pkin(2) * t539;
t444 = t821 * qJ(5) + t688 + t785;
t443 = t821 * pkin(4) + t788;
t438 = -t771 * t455 - t786 * t859 - t914;
t435 = -pkin(4) * t849 + t772 * t454 - t944;
t434 = -pkin(6) * t485 - pkin(7) * t827 - t774 * t524;
t433 = -t773 * t491 + t776 * t520 + t958;
t432 = t772 * t455 + t786 * t803 + t915;
t431 = -pkin(6) * t521 - t774 * t493 + t777 * t496;
t430 = -t943 + t771 * t454 + (pkin(3) + t861) * t874;
t429 = pkin(3) * t563 + qJ(4) * t437;
t428 = -t773 * t487 + t776 * t504 - t932;
t427 = t776 * t491 + t773 * t520 - t956;
t420 = t772 * t462 - t771 * t463;
t419 = t771 * t462 + t772 * t463;
t418 = -t436 - t913;
t417 = t776 * t487 + t773 * t504 + t927;
t412 = t437 + t902;
t411 = t489 - t964;
t410 = t488 + t942;
t409 = t688 + (-t868 - t862) * qJ(5) + (-t873 - t767) * pkin(4) + t797 + t942;
t408 = -t771 * t443 + t772 * t444 - t913;
t407 = t772 * t443 + t771 * t444 + t902;
t406 = t776 * t437 - t841;
t405 = t773 * t437 + t830;
t404 = -qJ(5) * t872 + 0.2e1 * t818 + (-t619 + t862) * pkin(4) - t869 + t964;
t402 = -qJ(4) * t419 + (-pkin(4) * t771 + t859) * t483;
t401 = -t773 * t432 + t776 * t438 - t932;
t400 = -pkin(4) * t875 + qJ(5) * t573 + t403;
t399 = -t773 * t430 + t776 * t435 - t958;
t398 = t776 * t432 + t773 * t438 + t927;
t397 = -t773 * t419 + t776 * t420;
t396 = t776 * t419 + t773 * t420;
t395 = t776 * t430 + t773 * t435 + t956;
t394 = qJ(4) * t420 + (-t803 + t861) * t483;
t393 = -t774 * t427 + t777 * t433 + t965;
t392 = -t774 * t417 + t777 * t428 - t946;
t391 = -t773 * t412 + t776 * t418 - t933;
t390 = t776 * t412 + t773 * t418 + t926;
t389 = -t774 * t405 + t777 * t406;
t388 = t777 * t405 + t774 * t406;
t387 = -pkin(7) * t405 - qJ(4) * t830 - t773 * t429;
t386 = t778 * t389 - t775 * t563;
t385 = t775 * t389 + t778 * t563;
t384 = pkin(2) * t563 + pkin(7) * t406 - qJ(4) * t841 + t776 * t429;
t383 = -t773 * t407 + t776 * t408 - t933;
t382 = t776 * t407 + t773 * t408 + t926;
t381 = -t774 * t396 + t777 * t397;
t380 = t777 * t396 + t774 * t397;
t379 = -t774 * t398 + t777 * t401 - t946;
t378 = t778 * t381 - t775 * t483;
t377 = t775 * t381 + t778 * t483;
t376 = -t774 * t395 + t777 * t399 - t965;
t375 = -pkin(1) * t388 - pkin(2) * t405 - pkin(3) * t436;
t374 = -pkin(7) * t396 - t773 * t394 + t776 * t402;
t373 = -t774 * t390 + t777 * t391 - t945;
t372 = pkin(2) * t483 + pkin(7) * t397 + t776 * t394 + t773 * t402;
t371 = -pkin(1) * t380 - pkin(2) * t396 - pkin(3) * t419 - pkin(4) * t463 - qJ(5) * t462;
t370 = -t774 * t382 + t777 * t383 - t945;
t369 = -pkin(6) * t388 - t774 * t384 + t777 * t387;
t368 = -pkin(6) * t380 - t774 * t372 + t777 * t374;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t743, -t744, 0, t707, 0, 0, 0, 0, 0, 0, t679, t680, t704, t633, 0, 0, 0, 0, 0, 0, t538, t545, t503, t481, 0, 0, 0, 0, 0, 0, t937, t448, t936, t386, 0, 0, 0, 0, 0, 0, t937, t936, -t448, t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t744, -t743, 0, t706, 0, 0, 0, 0, 0, 0, t677, t678, t703, t632, 0, 0, 0, 0, 0, 0, t537, t544, t502, t480, 0, 0, 0, 0, 0, 0, t939, t446, t938, t385, 0, 0, 0, 0, 0, 0, t939, t938, -t446, t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t708, t709, 0, -t655, 0, 0, 0, 0, 0, 0, t559, t584, t521, t485, 0, 0, 0, 0, 0, 0, t921, -t450, t923, t388, 0, 0, 0, 0, 0, 0, t921, t923, t450, t380; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t744, 0, -t743, 0, t796, -t725, -t706, -pkin(5) * t706, t778 * t718 - t800, t778 * t702 - t775 * t746, t778 * t712 + t774 * t812, t778 * t717 + t800, t778 * t710 + t775 * t762, t775 * qJDD(2) + t778 * t735, -pkin(5) * t677 - t775 * t659 + t778 * t669, -pkin(5) * t678 - t775 * t660 + t778 * t670, -pkin(5) * t703 + t778 * t655, -pkin(5) * t632 - (pkin(1) * t775 - pkin(6) * t778) * t655, t778 * t562 + t808, t778 * t522 - t775 * t700, t778 * t598 - t775 * t643, t778 * t561 - t808, t778 * t599 - t775 * t639, t778 * t612 + t831, -pkin(5) * t537 + t778 * t482 - t775 * t497, -pkin(5) * t544 + t778 * t490 - t775 * t505, -pkin(5) * t502 + t778 * t431 - t775 * t492, -pkin(5) * t480 + t778 * t434 - t775 * t445, t794, t961, t948, t904, -t959, t903, t778 * t392 - t775 * t410 - t950, t778 * t393 - t775 * t411 - t967, t778 * t373 - t775 * t403 - t951, -pkin(5) * t385 + t778 * t369 - t775 * t375, t794, t948, -t961, t903, t959, t904, t778 * t379 - t775 * t409 - t950, t778 * t370 - t775 * t400 - t951, t778 * t376 - t775 * t404 + t967, -pkin(5) * t377 + t778 * t368 - t775 * t371; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t743, 0, t744, 0, t725, t796, t707, pkin(5) * t707, t775 * t718 + t799, t775 * t702 + t778 * t746, t775 * t712 - t774 * t811, t775 * t717 - t799, t775 * t710 - t777 * t811, -t778 * qJDD(2) + t775 * t735, pkin(5) * t679 + t778 * t659 + t775 * t669, pkin(5) * t680 + t778 * t660 + t775 * t670, pkin(5) * t704 + t775 * t655, pkin(5) * t633 - (-pkin(1) * t778 - pkin(6) * t775) * t655, t775 * t562 - t806, t775 * t522 + t778 * t700, t775 * t598 + t778 * t643, t775 * t561 + t806, t775 * t599 + t778 * t639, t775 * t612 - t758, pkin(5) * t538 + t775 * t482 + t778 * t497, pkin(5) * t545 + t775 * t490 + t778 * t505, pkin(5) * t503 + t775 * t431 + t778 * t492, pkin(5) * t481 + t775 * t434 + t778 * t445, t795, t962, t949, t905, -t960, t906, t775 * t392 + t778 * t410 + t952, t775 * t393 + t778 * t411 + t966, t775 * t373 + t778 * t403 + t953, pkin(5) * t386 + t775 * t369 + t778 * t375, t795, t949, -t962, t906, t960, t905, t775 * t379 + t778 * t409 + t952, t775 * t370 + t778 * t400 + t953, t775 * t376 + t778 * t404 - t966, pkin(5) * t378 + t775 * t368 + t778 * t371; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t749, t750, 0, 0, (t739 + t804) * t774, t777 * t738 + t774 * t741, t777 * t752 + t835, (t740 - t805) * t777, t774 * t754 + t825, 0, pkin(1) * t741 + pkin(6) * t711 + t826, -pkin(1) * t738 + pkin(6) * t713 - t837, pkin(1) * t745 + pkin(6) * t742 + t656, pkin(1) * t731 + pkin(6) * t656, t777 * t636 + t774 * t637, t777 * t580 + t774 * t582, t777 * t646 + t774 * t648, t777 * t634 + t774 * t635, t777 * t647 + t774 * t649, t777 * t671 + t774 * t672, -pkin(1) * t638 + pkin(6) * t560 + t777 * t549 + t774 * t600, -pkin(1) * t870 + pkin(6) * t585 + t777 * t554 + t774 * t603, -pkin(1) * t663 + pkin(6) * t523 + t777 * t493 + t774 * t496, pkin(1) * t664 + pkin(6) * t486 - pkin(7) * t838 + t777 * t524, t458, -t955, t934, t883, t954, t880, t777 * t417 + t774 * t428 + t941, t777 * t427 + t774 * t433 - t963, t777 * t390 + t774 * t391 + t940, pkin(1) * t563 + pkin(6) * t389 + t777 * t384 + t774 * t387, t458, t934, t955, t880, -t954, t883, t777 * t398 + t774 * t401 + t941, t777 * t382 + t774 * t383 + t940, t777 * t395 + t774 * t399 + t963, pkin(1) * t483 + pkin(6) * t381 + t777 * t372 + t774 * t374;];
tauB_reg = t1;
