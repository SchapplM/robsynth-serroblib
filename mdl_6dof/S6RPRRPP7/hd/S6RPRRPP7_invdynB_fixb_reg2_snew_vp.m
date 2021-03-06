% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RPRRPP7_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP7_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:48:41
% EndTime: 2019-05-05 21:49:05
% DurationCPUTime: 16.28s
% Computational Cost: add. (19613->531), mult. (38073->653), div. (0->0), fcn. (23347->6), ass. (0->379)
t701 = sin(qJ(4));
t704 = cos(qJ(4));
t705 = cos(qJ(3));
t763 = qJD(1) * t705;
t663 = -t704 * qJD(3) + t701 * t763;
t702 = sin(qJ(3));
t758 = qJD(1) * qJD(3);
t747 = t702 * t758;
t754 = qJDD(1) * t705;
t669 = -t747 + t754;
t729 = -t701 * qJDD(3) - t704 * t669;
t720 = qJD(4) * t663 + t729;
t690 = qJD(1) * t702 + qJD(4);
t779 = t663 * t690;
t561 = t779 - t720;
t665 = qJD(3) * t701 + t704 * t763;
t743 = -t704 * qJDD(3) + t701 * t669;
t721 = (qJD(4) - t690) * t665 + t743;
t488 = -t561 * t701 + t704 * t721;
t660 = t665 ^ 2;
t819 = t663 ^ 2;
t583 = -t819 - t660;
t445 = t488 * t702 + t583 * t705;
t482 = t561 * t704 + t701 * t721;
t703 = sin(qJ(1));
t706 = cos(qJ(1));
t417 = t445 * t703 + t482 * t706;
t942 = pkin(6) * t417;
t938 = t445 * t706 - t482 * t703;
t941 = pkin(6) * t938;
t815 = pkin(7) + pkin(1);
t940 = t445 * t815;
t448 = t488 * t705 - t583 * t702;
t939 = t448 * t815;
t937 = qJ(2) * t482 - t940;
t936 = pkin(2) * t482 - t939;
t759 = qJD(4) + t690;
t555 = t665 * t759 + t743;
t565 = t779 + t720;
t802 = t565 * t701;
t490 = -t555 * t704 + t802;
t611 = -t819 + t660;
t784 = t611 * t705;
t454 = t490 * t702 - t784;
t800 = t565 * t704;
t484 = -t555 * t701 - t800;
t935 = t454 * t703 + t484 * t706;
t934 = t454 * t706 - t484 * t703;
t818 = t690 ^ 2;
t827 = t819 - t818;
t746 = t705 * t758;
t756 = qJDD(1) * t702;
t668 = -t746 - t756;
t661 = qJDD(4) - t668;
t780 = t663 * t665;
t715 = t661 + t780;
t834 = t715 * t701;
t846 = t704 * t827 - t834;
t874 = t702 * t846;
t476 = t721 * t705 + t874;
t833 = t715 * t704;
t850 = t701 * t827 + t833;
t867 = t706 * t850;
t933 = t476 * t703 + t867;
t872 = t703 * t850;
t932 = t476 * t706 - t872;
t931 = pkin(2) * t445 + pkin(3) * t583 + pkin(8) * t488 - qJ(2) * t448;
t632 = t660 - t818;
t716 = t661 - t780;
t831 = t716 * t704;
t853 = t632 * t701 + t831;
t473 = -t561 * t705 + t702 * t853;
t832 = t716 * t701;
t520 = t632 * t704 - t832;
t894 = t473 * t703 - t520 * t706;
t930 = t473 * t706 + t520 * t703;
t829 = t660 + t818;
t864 = -t829 * t701 + t833;
t901 = t702 * t864;
t895 = -t565 * t705 + t901;
t865 = t829 * t704 + t834;
t900 = t703 * t865;
t911 = -t706 * t895 + t900;
t928 = pkin(6) * t911;
t898 = t706 * t865;
t912 = t703 * t895 + t898;
t927 = pkin(6) * t912;
t926 = pkin(8) * t482;
t869 = t705 * t846;
t920 = t721 * t702 - t869;
t919 = pkin(3) * t482 + qJ(5) * t721;
t785 = t611 * t702;
t917 = t490 * t705 + t785;
t596 = qJD(4) * t665 + t743;
t639 = t690 * t665;
t554 = t596 + t639;
t828 = -t818 - t819;
t845 = t704 * t828 - t832;
t875 = t702 * t845;
t861 = -t554 * t705 + t875;
t849 = t701 * t828 + t831;
t873 = t703 * t849;
t890 = -t706 * t861 + t873;
t916 = pkin(6) * t890;
t868 = t706 * t849;
t893 = t703 * t861 + t868;
t915 = pkin(6) * t893;
t858 = t561 * t702 + t705 * t853;
t899 = t705 * t864;
t892 = t565 * t702 + t899;
t907 = pkin(2) * t865;
t910 = -t815 * t892 + t907;
t903 = qJ(2) * t865;
t909 = -t815 * t895 + t903;
t905 = pkin(8) * t864;
t908 = pkin(2) * t895 - pkin(3) * t565 - qJ(2) * t892 + t905;
t906 = pkin(3) * t865;
t904 = pkin(8) * t865;
t902 = t565 * qJ(5);
t762 = qJD(5) * t690;
t679 = -0.2e1 * t762;
t897 = -qJ(5) * t715 + t679 - t906;
t870 = t705 * t845;
t860 = t554 * t702 + t870;
t886 = pkin(2) * t849;
t889 = -t815 * t860 + t886;
t881 = qJ(2) * t849;
t888 = -t815 * t861 + t881;
t883 = pkin(8) * t845;
t887 = pkin(2) * t861 - pkin(3) * t554 - qJ(2) * t860 + t883;
t884 = pkin(3) * t849;
t882 = pkin(8) * t849;
t880 = qJ(5) * t583;
t863 = -qJ(5) * t828 - t884;
t771 = t690 * t704;
t736 = t665 * t771 - t701 * t720;
t604 = t705 * t780;
t772 = t690 * t701;
t625 = t665 * t772;
t765 = -t704 * t720 - t625;
t855 = -t765 * t702 + t604;
t852 = t703 * t736 + t706 * t855;
t851 = -t703 * t855 + t706 * t736;
t856 = qJ(6) * t561;
t602 = pkin(4) * t663 - qJ(5) * t665;
t757 = qJD(2) * qJD(1);
t696 = 0.2e1 * t757;
t681 = t706 * g(1) + t703 * g(2);
t698 = qJDD(1) * qJ(2);
t728 = t681 - t698;
t725 = t696 - t728;
t731 = -t669 + t747;
t732 = -t668 + t746;
t708 = qJD(1) ^ 2;
t835 = t708 * t815;
t540 = pkin(3) * t732 + pkin(8) * t731 + t725 - t835;
t680 = t703 * g(1) - t706 * g(2);
t733 = qJDD(2) - t680;
t722 = -t708 * qJ(2) + t733;
t637 = -qJDD(1) * t815 + t722;
t606 = -t705 * g(3) + t702 * t637;
t707 = qJD(3) ^ 2;
t814 = pkin(3) * t702;
t739 = -pkin(8) * t705 + t814;
t726 = t708 * t739;
t568 = -t707 * pkin(3) + qJDD(3) * pkin(8) - t702 * t726 + t606;
t478 = -t704 * t540 + t701 * t568;
t723 = -t661 * pkin(4) - qJ(5) * t818 + qJDD(5) + t478;
t411 = -t661 * pkin(5) - t665 * (-pkin(5) * t663 + (2 * qJD(6)) - t602) + t723 - t856;
t751 = t663 * t771;
t734 = t625 - t751;
t781 = t661 * t705;
t822 = -t702 * t734 + t781;
t724 = (-t663 * t701 - t665 * t704) * t690;
t837 = t706 * t724;
t848 = -t703 * t822 + t837;
t735 = -t704 * t596 + t663 * t772;
t727 = t596 * t701 + t751;
t823 = -t702 * t727 - t604;
t847 = -t703 * t823 + t706 * t735;
t840 = t703 * t724;
t844 = t706 * t822 + t840;
t843 = t703 * t735 + t706 * t823;
t603 = t702 * t780;
t824 = t705 * t727 - t603;
t766 = t705 * t765 + t603;
t782 = t661 * t702;
t821 = t705 * t734 + t782;
t623 = -pkin(5) * t690 - qJ(6) * t665;
t817 = 0.2e1 * qJD(5);
t820 = -t596 * pkin(5) - t819 * qJ(6) + qJDD(6) + (t817 + t623) * t665;
t816 = pkin(4) + pkin(5);
t813 = t596 * pkin(4);
t812 = qJDD(1) * pkin(1);
t810 = t554 * t701;
t808 = t554 * t704;
t605 = t702 * g(3) + t705 * t637;
t567 = qJDD(3) * pkin(3) + t707 * pkin(8) - t705 * t726 + t605;
t796 = t567 * t701;
t795 = t567 * t704;
t699 = t702 ^ 2;
t700 = t705 ^ 2;
t764 = t699 + t700;
t671 = t764 * qJDD(1);
t778 = t671 * t703;
t777 = t671 * t706;
t750 = t702 * t705 * t708;
t676 = qJDD(3) + t750;
t776 = t676 * t702;
t775 = t676 * t705;
t677 = qJDD(3) - t750;
t774 = t677 * t702;
t773 = t677 * t705;
t770 = t699 * t708;
t769 = t700 * t708;
t627 = t728 - 0.2e1 * t757 + t835;
t768 = t702 * t627;
t767 = t705 * t627;
t479 = t701 * t540 + t704 * t568;
t761 = qJD(6) * t663;
t755 = qJDD(1) * t703;
t753 = qJDD(1) * t706;
t752 = t665 * t817;
t749 = pkin(3) * t705 + pkin(2);
t748 = (-t554 - t596) * pkin(4);
t422 = t478 * t701 + t704 * t479;
t642 = -t708 * pkin(1) + t725;
t643 = -t722 + t812;
t578 = t706 * t642 - t643 * t703;
t615 = -t680 * t703 - t706 * t681;
t741 = t703 * t750;
t740 = t706 * t750;
t672 = -t703 * t708 + t753;
t738 = pkin(6) * t672 + g(3) * t703;
t673 = t706 * t708 + t755;
t737 = -pkin(6) * t673 + g(3) * t706;
t730 = pkin(4) * t818 - t661 * qJ(5) + t663 * t602 - t479;
t421 = -t478 * t704 + t479 * t701;
t538 = t705 * t605 + t702 * t606;
t539 = -t605 * t702 + t606 * t705;
t575 = t642 * t703 + t643 * t706;
t614 = t680 * t706 - t681 * t703;
t678 = 0.2e1 * t762;
t442 = t678 - t730;
t719 = pkin(5) * t819 - t690 * t623 + t730;
t718 = t665 * t602 + t723;
t714 = -t596 * qJ(6) + t719;
t713 = -pkin(4) * t639 + t567;
t712 = t713 - t902;
t711 = t712 + t752;
t710 = t713 - t813 - 0.2e1 * t902;
t709 = t712 + t820;
t687 = -t707 - t769;
t686 = t707 - t769;
t685 = -t707 - t770;
t684 = -t707 + t770;
t675 = (-t699 + t700) * t708;
t674 = t764 * t708;
t670 = -0.2e1 * t747 + t754;
t667 = 0.2e1 * t746 + t756;
t662 = t764 * t758;
t649 = -0.2e1 * t761;
t629 = -t669 * t702 - t700 * t758;
t628 = -t668 * t705 - t699 * t758;
t622 = -t687 * t702 - t775;
t621 = t685 * t705 - t774;
t620 = t687 * t705 - t776;
t619 = -t686 * t705 - t774;
t618 = t685 * t702 + t773;
t617 = -t684 * t702 - t775;
t610 = -t674 * t706 - t778;
t609 = -t674 * t703 + t777;
t601 = t667 * t702 - t670 * t705;
t580 = t620 * t703 + t670 * t706;
t579 = t618 * t703 + t667 * t706;
t577 = -t620 * t706 + t670 * t703;
t576 = -t618 * t706 + t667 * t703;
t574 = (-t663 * t704 + t665 * t701) * t690;
t563 = t663 * t759 + t729;
t553 = t596 - t639;
t532 = -t574 * t702 + t781;
t511 = -pkin(2) * t674 - t539;
t510 = pkin(2) * t620 - qJ(2) * t622 - t606;
t509 = pkin(2) * t618 - qJ(2) * t621 + t605;
t502 = pkin(2) * t667 - t621 * t815 - t767;
t501 = pkin(2) * t670 - t622 * t815 + t768;
t500 = t538 * t703 - t627 * t706;
t499 = -t538 * t706 - t627 * t703;
t492 = -qJ(5) * t554 + qJ(6) * t716;
t486 = t808 - t802;
t480 = t810 + t800;
t474 = -t553 * t705 - t874;
t470 = -t795 + t904;
t469 = pkin(2) * t538 - qJ(2) * t539;
t468 = -t796 - t882;
t466 = -t563 * t702 - t899;
t463 = t563 * t705 - t901;
t460 = t555 * t702 + t870;
t457 = -t555 * t705 + t875;
t453 = -t486 * t702 - t784;
t452 = -pkin(2) * t627 - t539 * t815;
t451 = -qJ(6) * t715 - t565 * t816;
t444 = t711 - t813;
t441 = t479 + t906;
t440 = t478 - t884;
t438 = t463 * t703 - t898;
t435 = -t463 * t706 - t900;
t433 = t718 - t880;
t431 = t457 * t703 + t868;
t428 = -t457 * t706 + t873;
t426 = pkin(4) * t561 + t919;
t425 = -pkin(4) * t583 + t442;
t424 = t748 + t711;
t423 = t710 + t752;
t420 = t709 - t813;
t413 = -t561 * t816 - t919;
t412 = t678 - t714 + 0.2e1 * t761;
t410 = t422 * t705 - t567 * t702;
t409 = t422 * t702 + t567 * t705;
t408 = -pkin(4) * t716 + t718 + t863;
t407 = -pkin(4) * t829 + t730 + t897;
t406 = qJ(6) * t829 + t710 + t820;
t405 = -t421 + t926;
t404 = -qJ(5) * t808 - t424 * t701 - t882;
t403 = pkin(4) * t802 + t423 * t704 - t904;
t402 = t442 * t704 + t701 * t718;
t401 = t442 * t701 - t704 * t718;
t400 = -t411 + t856 + t880;
t399 = -pkin(5) * t554 - qJ(6) * t828 + t709 + t748;
t398 = t649 + t679 + t816 * t583 + (-t721 - t596) * qJ(6) + t719;
t397 = pkin(2) * t463 + pkin(3) * t563 - qJ(2) * t466 - t796 - t905;
t396 = -t816 * t829 + t649 + t714 + t897;
t395 = pkin(2) * t457 - pkin(3) * t555 - qJ(2) * t460 + t795 + t883;
t394 = -t716 * t816 + t411 + t863;
t393 = t409 * t703 + t421 * t706;
t392 = -t409 * t706 + t421 * t703;
t391 = qJ(5) * t420 - qJ(6) * t411;
t390 = t402 * t705 - t444 * t702;
t389 = t402 * t702 + t444 * t705;
t388 = -t425 * t701 + t433 * t704 + t926;
t387 = t406 * t704 - t451 * t701 - t904;
t386 = t411 * t701 + t412 * t704;
t385 = -t411 * t704 + t412 * t701;
t384 = -t399 * t701 + t492 * t704 - t882;
t383 = -t705 * t441 - t466 * t815 - t702 * t470 - t907;
t382 = -qJ(5) * t810 + t424 * t704 + t887;
t381 = -t705 * t440 - t460 * t815 - t702 * t468 + t886;
t380 = -pkin(4) * t800 + t423 * t701 + t908;
t379 = t422 - t931;
t378 = -pkin(8) * t401 + (-pkin(4) * t701 + qJ(5) * t704) * t444;
t377 = -pkin(3) * t401 + pkin(4) * t718 - qJ(5) * t442;
t376 = -qJ(6) * t412 + t420 * t816;
t375 = t386 * t705 - t420 * t702;
t374 = t386 * t702 + t420 * t705;
t373 = pkin(2) * t409 + pkin(3) * t567 + pkin(8) * t422 - qJ(2) * t410;
t372 = t406 * t701 + t451 * t704 + t908;
t371 = -t702 * t405 - t482 * t749 + t939;
t370 = t399 * t704 + t492 * t701 + t887;
t369 = -t398 * t701 + t400 * t704 - t926;
t368 = t425 * t704 + t433 * t701 - t931;
t367 = t389 * t703 + t401 * t706;
t366 = -t389 * t706 + t401 * t703;
t365 = -t702 * t404 - t705 * t408 + t889;
t364 = -t702 * t403 - t705 * t407 + t910;
t363 = -t702 * t388 - t705 * t426 - t936;
t362 = t398 * t704 + t400 * t701 + t931;
t361 = -t815 * t410 + (pkin(8) * t702 + t749) * t421;
t360 = -t702 * t387 - t705 * t396 + t910;
t359 = t374 * t703 + t385 * t706;
t358 = -t374 * t706 + t385 * t703;
t357 = -t702 * t384 - t705 * t394 + t889;
t356 = -pkin(3) * t385 - qJ(5) * t412 + t411 * t816;
t355 = -t702 * t369 - t705 * t413 + t936;
t354 = -pkin(8) * t385 - t376 * t701 + t391 * t704;
t353 = pkin(2) * t389 + pkin(8) * t402 - qJ(2) * t390 + (pkin(4) * t704 + qJ(5) * t701 + pkin(3)) * t444;
t352 = pkin(2) * t401 - t705 * t377 - t702 * t378 - t390 * t815;
t351 = pkin(2) * t374 + pkin(3) * t420 + pkin(8) * t386 - qJ(2) * t375 + t376 * t704 + t391 * t701;
t350 = pkin(2) * t385 - t702 * t354 - t705 * t356 - t375 * t815;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t673, -t672, 0, t615, 0, 0, 0, 0, 0, 0, 0, t673, t672, t578, 0, 0, 0, 0, 0, 0, t579, t580, t610, t500, 0, 0, 0, 0, 0, 0, t431, t438, -t417, t393, 0, 0, 0, 0, 0, 0, t893, -t417, t912, t367, 0, 0, 0, 0, 0, 0, t893, t912, t417, t359; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t672, -t673, 0, t614, 0, 0, 0, 0, 0, 0, 0, -t672, t673, t575, 0, 0, 0, 0, 0, 0, t576, t577, t609, t499, 0, 0, 0, 0, 0, 0, t428, t435, t938, t392, 0, 0, 0, 0, 0, 0, t890, t938, t911, t366, 0, 0, 0, 0, 0, 0, t890, t911, -t938, t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t621, t622, 0, t539, 0, 0, 0, 0, 0, 0, t460, t466, -t448, t410, 0, 0, 0, 0, 0, 0, t860, -t448, t892, t390, 0, 0, 0, 0, 0, 0, t860, t892, t448, t375; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t672, 0, -t673, 0, -t738, -t737, -t614, -pkin(6) * t614, 0, -t672, t673, 0, 0, 0, -t575, t738, t737, -pkin(6) * t575 + (-pkin(1) * t703 + qJ(2) * t706) * g(3), -t629 * t703 + t740, -t601 * t703 + t675 * t706, -t619 * t703 + t705 * t753, -t628 * t703 - t740, -t617 * t703 - t702 * t753, qJDD(3) * t706 - t662 * t703, -pkin(6) * t576 - t502 * t703 + t509 * t706, -pkin(6) * t577 - t501 * t703 + t510 * t706, -pkin(2) * t777 - pkin(6) * t609 - t511 * t703, -pkin(6) * t499 - t452 * t703 + t469 * t706, t851, t935, t894, t847, t933, t848, -pkin(6) * t428 - t381 * t703 + t395 * t706, -pkin(6) * t435 - t383 * t703 + t397 * t706, -t371 * t703 + t379 * t706 - t941, -pkin(6) * t392 - t361 * t703 + t373 * t706, t851, t894, -t935, t848, -t933, t847, -t365 * t703 + t382 * t706 - t916, -t363 * t703 + t368 * t706 - t941, -t364 * t703 + t380 * t706 - t928, -pkin(6) * t366 - t352 * t703 + t353 * t706, t851, -t453 * t703 + t480 * t706, -t894, t847, -t474 * t703 + t867, -t532 * t703 + t837, -t357 * t703 + t370 * t706 - t916, -t360 * t703 + t372 * t706 - t928, -t355 * t703 + t362 * t706 + t941, -pkin(6) * t358 - t350 * t703 + t351 * t706; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t673, 0, t672, 0, t737, -t738, t615, pkin(6) * t615, 0, -t673, -t672, 0, 0, 0, t578, -t737, t738, pkin(6) * t578 + (pkin(1) * t706 + qJ(2) * t703) * g(3), t629 * t706 + t741, t601 * t706 + t675 * t703, t619 * t706 + t703 * t754, t628 * t706 - t741, t617 * t706 - t702 * t755, qJDD(3) * t703 + t662 * t706, pkin(6) * t579 + t502 * t706 + t509 * t703, pkin(6) * t580 + t501 * t706 + t510 * t703, -pkin(2) * t778 + pkin(6) * t610 + t511 * t706, pkin(6) * t500 + t452 * t706 + t469 * t703, t852, -t934, -t930, t843, -t932, t844, pkin(6) * t431 + t381 * t706 + t395 * t703, pkin(6) * t438 + t383 * t706 + t397 * t703, t371 * t706 + t379 * t703 - t942, pkin(6) * t393 + t361 * t706 + t373 * t703, t852, -t930, t934, t844, t932, t843, t365 * t706 + t382 * t703 + t915, t363 * t706 + t368 * t703 - t942, t364 * t706 + t380 * t703 + t927, pkin(6) * t367 + t352 * t706 + t353 * t703, t852, t453 * t706 + t480 * t703, t930, t843, t474 * t706 + t872, t532 * t706 + t840, t357 * t706 + t370 * t703 + t915, t360 * t706 + t372 * t703 + t927, t355 * t706 + t362 * t703 + t942, pkin(6) * t359 + t350 * t706 + t351 * t703; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t680, t681, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t733 - 0.2e1 * t812, -t681 + t696 + 0.2e1 * t698, pkin(1) * t643 + qJ(2) * t642, -t731 * t705, -t667 * t705 - t670 * t702, -t686 * t702 + t773, t732 * t702, t684 * t705 - t776, 0, qJ(2) * t667 - t618 * t815 - t768, qJ(2) * t670 - t620 * t815 - t767, -qJ(2) * t674 + t671 * t815 - t538, -qJ(2) * t627 - t538 * t815, t766, t917, t858, t824, -t920, t821, -t702 * t440 - t457 * t815 + t705 * t468 + t881, -t702 * t441 - t463 * t815 + t705 * t470 - t903, t705 * t405 - (qJ(2) + t814) * t482 + t940, -t815 * t409 + (qJ(2) + t739) * t421, t766, t858, -t917, t821, t920, t824, t705 * t404 - t702 * t408 + t888, t705 * t388 - t702 * t426 - t937, t705 * t403 - t702 * t407 + t909, qJ(2) * t401 - t702 * t377 + t705 * t378 - t389 * t815, t766, t486 * t705 - t785, -t858, t824, -t553 * t702 + t869, t574 * t705 + t782, t705 * t384 - t702 * t394 + t888, t705 * t387 - t702 * t396 + t909, t705 * t369 - t702 * t413 + t937, qJ(2) * t385 + t705 * t354 - t702 * t356 - t374 * t815;];
tauB_reg  = t1;
