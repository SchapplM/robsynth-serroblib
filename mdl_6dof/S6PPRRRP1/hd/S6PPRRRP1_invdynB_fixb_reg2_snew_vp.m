% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 20:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6PPRRRP1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynB_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:29:12
% EndTime: 2019-05-04 20:29:43
% DurationCPUTime: 28.09s
% Computational Cost: add. (180178->831), mult. (332704->1337), div. (0->0), fcn. (267143->14), ass. (0->633)
t764 = sin(pkin(11));
t768 = cos(pkin(11));
t740 = g(1) * t768 + g(2) * t764;
t763 = sin(pkin(12));
t767 = cos(pkin(12));
t739 = g(1) * t764 - t768 * g(2);
t760 = g(3) - qJDD(1);
t766 = sin(pkin(6));
t770 = cos(pkin(6));
t815 = t739 * t770 - t760 * t766;
t662 = -t767 * t740 + t763 * t815;
t775 = sin(qJ(3));
t778 = cos(qJ(3));
t661 = -t763 * t740 - t767 * t815;
t704 = t739 * t766 + t760 * t770 - qJDD(2);
t765 = sin(pkin(7));
t769 = cos(pkin(7));
t963 = t661 * t769 + t704 * t765;
t581 = t778 * t662 - t775 * t963;
t779 = qJD(3) ^ 2;
t575 = -t779 * pkin(3) + qJDD(3) * pkin(9) + t581;
t623 = -t661 * t765 + t769 * t704;
t777 = cos(qJ(4));
t620 = t777 * t623;
t774 = sin(qJ(4));
t955 = pkin(4) * t777;
t874 = -pkin(10) * t774 - t955;
t729 = t874 * qJD(3);
t959 = qJD(4) ^ 2;
t520 = (qJD(3) * t729 + t575) * t774 - qJDD(4) * pkin(4) - t959 * pkin(10) + t620;
t894 = qJDD(3) * t775;
t813 = t778 * t779 + t894;
t713 = t813 * t765;
t715 = t813 * t769;
t892 = qJDD(3) * t778;
t735 = -t775 * t779 + t892;
t818 = t715 * t767 + t735 * t763;
t645 = -t766 * t713 + t770 * t818;
t677 = t715 * t763 - t735 * t767;
t596 = t645 * t764 + t677 * t768;
t966 = t645 * t768 - t677 * t764;
t896 = qJD(3) * qJD(4);
t752 = t774 * t896;
t893 = qJDD(3) * t777;
t732 = -t752 + t893;
t723 = -qJDD(5) + t732;
t773 = sin(qJ(5));
t776 = cos(qJ(5));
t900 = qJD(3) * t774;
t726 = -t776 * qJD(4) + t773 * t900;
t728 = qJD(4) * t773 + t776 * t900;
t927 = t726 * t728;
t879 = t723 + t927;
t962 = t879 * pkin(5);
t886 = t777 * t896;
t895 = qJDD(3) * t774;
t731 = t886 + t895;
t878 = -t776 * qJDD(4) + t773 * t731;
t674 = -qJD(5) * t728 - t878;
t899 = qJD(3) * t777;
t749 = -qJD(5) + t899;
t703 = -pkin(5) * t749 - qJ(6) * t728;
t721 = t726 ^ 2;
t492 = -t674 * pkin(5) - t721 * qJ(6) + t728 * t703 + qJDD(6) + t520;
t932 = t879 * t773;
t931 = t879 * t776;
t711 = t726 * t749;
t782 = qJD(5) * t726 - qJDD(4) * t773 - t731 * t776;
t961 = t711 - t782;
t650 = (qJD(5) + t749) * t728 + t878;
t960 = t770 * t713 + t766 * t818;
t722 = t728 ^ 2;
t747 = t749 ^ 2;
t679 = -t747 - t721;
t610 = t679 * t773 - t931;
t958 = pkin(4) * t610;
t685 = -t722 - t747;
t664 = t723 - t927;
t934 = t664 * t773;
t617 = t685 * t776 + t934;
t957 = pkin(4) * t617;
t956 = pkin(4) * t774;
t954 = pkin(8) * t765;
t953 = pkin(8) * t769;
t654 = t711 + t782;
t585 = -t650 * t776 - t654 * t773;
t663 = -t721 - t722;
t552 = t585 * t774 - t663 * t777;
t952 = pkin(9) * t552;
t611 = t679 * t776 + t932;
t649 = (qJD(5) - t749) * t728 + t878;
t562 = t611 * t774 - t649 * t777;
t951 = pkin(9) * t562;
t933 = t664 * t776;
t618 = -t685 * t773 + t933;
t565 = t618 * t774 - t777 * t961;
t950 = pkin(9) * t565;
t583 = -t650 * t773 + t654 * t776;
t949 = pkin(10) * t583;
t948 = pkin(10) * t610;
t947 = pkin(10) * t617;
t553 = t585 * t777 + t663 * t774;
t839 = t553 * t775 - t583 * t778;
t461 = t769 * t552 + t765 * t839;
t462 = -t765 * t552 + t769 * t839;
t497 = t553 * t778 + t583 * t775;
t851 = t462 * t767 + t497 * t763;
t370 = -t766 * t461 + t770 * t851;
t412 = -t462 * t763 + t497 * t767;
t335 = t370 * t768 + t412 * t764;
t946 = qJ(1) * t335;
t563 = t611 * t777 + t649 * t774;
t837 = t563 * t775 - t610 * t778;
t475 = t769 * t562 + t765 * t837;
t476 = -t765 * t562 + t769 * t837;
t515 = t563 * t778 + t610 * t775;
t848 = t476 * t767 + t515 * t763;
t384 = -t766 * t475 + t770 * t848;
t430 = -t476 * t763 + t515 * t767;
t345 = t384 * t768 + t430 * t764;
t945 = qJ(1) * t345;
t566 = t618 * t777 + t774 * t961;
t836 = t566 * t775 - t617 * t778;
t478 = t769 * t565 + t765 * t836;
t479 = -t765 * t565 + t769 * t836;
t524 = t566 * t778 + t617 * t775;
t847 = t479 * t767 + t524 * t763;
t390 = -t766 * t478 + t770 * t847;
t435 = -t479 * t763 + t524 * t767;
t349 = t390 * t768 + t435 * t764;
t944 = qJ(1) * t349;
t943 = qJ(2) * t766;
t942 = qJ(2) * t770;
t536 = t777 * t575 - t774 * t623;
t521 = -pkin(4) * t959 + qJDD(4) * pkin(10) + t729 * t899 + t536;
t580 = t775 * t662 + t778 * t963;
t574 = -qJDD(3) * pkin(3) - t779 * pkin(9) + t580;
t869 = -t732 + t752;
t870 = t731 + t886;
t551 = pkin(4) * t869 - pkin(10) * t870 + t574;
t905 = -t773 * t521 + t776 * t551;
t889 = -qJ(6) * t782 - t905;
t806 = qJ(6) * t711 - t889;
t898 = qJD(6) * t728;
t451 = t806 - 0.2e1 * t898 - t962;
t941 = t451 * t773;
t940 = t451 * t776;
t939 = t520 * t773;
t938 = t520 * t776;
t937 = t574 * t774;
t936 = t574 * t777;
t930 = t704 * t764;
t928 = t704 * t768;
t748 = t774 * t779 * t777;
t741 = qJDD(4) + t748;
t926 = t741 * t774;
t742 = qJDD(4) - t748;
t925 = t742 * t774;
t924 = t742 * t777;
t923 = t749 * t773;
t922 = t749 * t776;
t758 = t774 ^ 2;
t921 = t758 * t779;
t920 = t764 * t760;
t919 = t766 * t704;
t917 = t766 * t767;
t916 = t768 * t760;
t915 = t770 * t704;
t913 = t775 * t623;
t912 = t778 * t623;
t911 = pkin(1) * t370 + t412 * t943;
t910 = pkin(1) * t384 + t430 * t943;
t909 = pkin(1) * t390 + t435 * t943;
t908 = pkin(2) * t462 + t497 * t954;
t907 = pkin(2) * t476 + t515 * t954;
t906 = pkin(2) * t479 + t524 * t954;
t473 = t776 * t521 + t773 * t551;
t904 = -pkin(3) * t583 + pkin(9) * t553;
t903 = -pkin(3) * t610 + pkin(9) * t563;
t902 = -pkin(3) * t617 + pkin(9) * t566;
t759 = t777 ^ 2;
t901 = t758 + t759;
t891 = t774 * t927;
t890 = t777 * t927;
t887 = qJDD(3) * t763 * t765;
t369 = t770 * t461 + t766 * t851;
t885 = -pkin(1) * t369 + t412 * t942;
t383 = t770 * t475 + t766 * t848;
t884 = -pkin(1) * t383 + t430 * t942;
t389 = t770 * t478 + t766 * t847;
t883 = -pkin(1) * t389 + t435 * t942;
t882 = -pkin(2) * t461 + t497 * t953;
t881 = -pkin(2) * t475 + t515 * t953;
t880 = -pkin(2) * t478 + t524 * t953;
t535 = t575 * t774 + t620;
t470 = t535 * t774 + t777 * t536;
t690 = -t739 * t764 - t768 * t740;
t876 = t775 * t748;
t875 = t778 * t748;
t469 = t535 * t777 - t536 * t774;
t734 = t901 * qJDD(3);
t756 = t759 * t779;
t737 = t756 + t921;
t686 = t734 * t778 - t737 * t775;
t873 = pkin(8) * t686 + t469 * t775;
t872 = -pkin(8) * t735 - t913;
t871 = -pkin(8) * t813 + t912;
t805 = t674 * qJ(6) - 0.2e1 * qJD(6) * t726 + t749 * t703 + t473;
t434 = -qJ(6) * t650 + (-t663 - t721) * pkin(5) + t805;
t718 = 0.2e1 * t898;
t439 = t718 + (-t654 - t711) * qJ(6) + t962 + t889;
t366 = -t434 * t773 + t439 * t776 - t949;
t545 = -pkin(4) * t583 - pkin(5) * t654;
t340 = t366 * t774 + t545 * t777 + t904;
t347 = t366 * t777 - t545 * t774 - t952;
t801 = -pkin(3) * t552 + pkin(4) * t663 - pkin(10) * t585;
t355 = -t434 * t776 - t439 * t773 + t801;
t860 = t347 * t775 + t355 * t778;
t293 = -t765 * t340 + t769 * t860 + t882;
t811 = (-t461 * t765 - t462 * t769) * pkin(8);
t300 = t778 * t347 - t775 * t355 + t811;
t868 = t293 * t767 + t300 * t763;
t407 = t473 * t773 + t776 * t905;
t397 = -t407 - t949;
t364 = t397 * t774 - t583 * t955 + t904;
t371 = t397 * t777 + t583 * t956 - t952;
t408 = t473 * t776 - t773 * t905;
t378 = -t408 + t801;
t856 = t371 * t775 + t378 * t778;
t297 = -t765 * t364 + t769 * t856 + t882;
t310 = t778 * t371 - t775 * t378 + t811;
t867 = t297 * t767 + t310 * t763;
t427 = t718 - t806 - t958 + 0.2e1 * t962;
t471 = -pkin(5) * t649 + qJ(6) * t679 - t492;
t446 = qJ(6) * t931 - t471 * t773 - t948;
t360 = t427 * t777 + t446 * t774 + t903;
t365 = -t427 * t774 + t446 * t777 - t951;
t800 = -pkin(3) * t562 + pkin(4) * t649 - pkin(10) * t611;
t413 = -qJ(6) * t932 - t471 * t776 + t800;
t858 = t365 * t775 + t413 * t778;
t299 = -t765 * t360 + t769 * t858 + t881;
t810 = (-t475 * t765 - t476 * t769) * pkin(8);
t320 = t778 * t365 - t775 * t413 + t810;
t866 = t299 * t767 + t320 * t763;
t431 = -t957 + (-t685 - t721) * pkin(5) + t805;
t488 = -qJ(6) * t685 + t492;
t598 = -pkin(5) * t961 + qJ(6) * t664;
t449 = t488 * t776 - t598 * t773 - t947;
t361 = t431 * t777 + t449 * t774 + t902;
t367 = -t431 * t774 + t449 * t777 - t950;
t799 = -pkin(3) * t565 + pkin(4) * t961 - pkin(10) * t618;
t415 = -t488 * t773 - t598 * t776 + t799;
t857 = t367 * t775 + t415 * t778;
t305 = -t765 * t361 + t769 * t857 + t880;
t809 = (-t478 * t765 - t479 * t769) * pkin(8);
t322 = t778 * t367 - t775 * t415 + t809;
t865 = t305 * t767 + t322 * t763;
t454 = -pkin(5) * t721 + t805;
t386 = t454 * t776 - t941;
t362 = t386 * t774 - t492 * t777;
t363 = t386 * t777 + t492 * t774;
t385 = t454 * t773 + t940;
t859 = t363 * t775 - t385 * t778;
t307 = -t765 * t362 + t769 * t859;
t327 = t363 * t778 + t385 * t775;
t864 = t307 * t767 + t327 * t763;
t458 = -t905 - t958;
t493 = t939 - t948;
t381 = t458 * t777 + t493 * t774 + t903;
t393 = -t458 * t774 + t493 * t777 - t951;
t450 = t800 + t938;
t853 = t393 * t775 + t450 * t778;
t312 = -t765 * t381 + t769 * t853 + t881;
t331 = t778 * t393 - t775 * t450 + t810;
t863 = t312 * t767 + t331 * t763;
t459 = t473 - t957;
t494 = t938 - t947;
t387 = t459 * t777 + t494 * t774 + t902;
t396 = -t459 * t774 + t494 * t777 - t950;
t453 = t799 - t939;
t852 = t396 * t775 + t453 * t778;
t315 = -t765 * t387 + t769 * t852 + t880;
t332 = t778 * t396 - t775 * t453 + t809;
t862 = t315 * t767 + t332 * t763;
t379 = t408 * t774 - t520 * t777;
t380 = t408 * t777 + t520 * t774;
t855 = t380 * t775 - t407 * t778;
t319 = -t765 * t379 + t769 * t855;
t339 = t380 * t778 + t407 * t775;
t861 = t319 * t767 + t339 * t763;
t849 = t470 * t775 - t574 * t778;
t392 = t469 * t765 + t769 * t849;
t444 = t470 * t778 + t574 * t775;
t854 = t392 * t767 + t444 * t763;
t584 = -t649 * t776 - t773 * t961;
t691 = -t722 + t721;
t558 = t584 * t774 + t691 * t777;
t559 = t584 * t777 - t691 * t774;
t582 = t649 * t773 - t776 * t961;
t838 = t559 * t775 + t582 * t778;
t464 = -t765 * t558 + t769 * t838;
t498 = t559 * t778 - t582 * t775;
t850 = t464 * t767 + t498 * t763;
t709 = -t722 + t747;
t635 = -t709 * t773 - t931;
t567 = t635 * t774 + t654 * t777;
t569 = t635 * t777 - t654 * t774;
t633 = -t709 * t776 + t932;
t835 = t569 * t775 + t633 * t778;
t486 = -t765 * t567 + t769 * t835;
t533 = t569 * t778 - t633 * t775;
t846 = t486 * t767 + t533 * t763;
t708 = t721 - t747;
t636 = t708 * t776 + t934;
t568 = t636 * t774 + t650 * t777;
t570 = t636 * t777 - t650 * t774;
t634 = -t708 * t773 + t933;
t834 = t570 * t775 + t634 * t778;
t487 = -t765 * t568 + t769 * t834;
t534 = t570 * t778 - t634 * t775;
t845 = t487 * t767 + t534 * t763;
t506 = t580 * t778 - t581 * t775;
t490 = -t506 * t769 + t765 * t623;
t507 = t580 * t775 + t581 * t778;
t844 = t490 * t767 + t507 * t763;
t640 = -t674 * t773 - t726 * t922;
t599 = t640 * t774 + t890;
t601 = t640 * t777 - t891;
t639 = -t674 * t776 + t726 * t923;
t833 = t601 * t775 + t639 * t778;
t501 = -t765 * t599 + t769 * t833;
t554 = t601 * t778 - t639 * t775;
t843 = t501 * t767 + t554 * t763;
t642 = t728 * t923 - t776 * t782;
t600 = t642 * t774 - t890;
t602 = t642 * t777 + t891;
t641 = t728 * t922 + t773 * t782;
t832 = t602 * t775 + t641 * t778;
t502 = -t765 * t600 + t769 * t832;
t555 = t602 * t778 - t641 * t775;
t842 = t502 * t767 + t555 * t763;
t504 = t507 * t769;
t841 = t504 * t767 + t506 * t763;
t658 = (t726 * t776 - t728 * t773) * t749;
t637 = t658 * t774 + t723 * t777;
t638 = t658 * t777 - t723 * t774;
t657 = (-t726 * t773 - t728 * t776) * t749;
t824 = t638 * t775 + t657 * t778;
t538 = -t765 * t637 + t769 * t824;
t578 = t638 * t778 - t657 * t775;
t840 = t538 * t767 + t578 * t763;
t730 = 0.2e1 * t886 + t895;
t733 = -0.2e1 * t752 + t893;
t682 = t730 * t777 + t733 * t774;
t683 = -t730 * t774 + t733 * t777;
t738 = t756 - t921;
t822 = t683 * t775 + t738 * t778;
t608 = -t765 * t682 + t769 * t822;
t655 = t683 * t778 - t738 * t775;
t831 = t608 * t767 + t655 * t763;
t725 = t777 * t741;
t746 = -t756 - t959;
t696 = t746 * t774 + t725;
t700 = t746 * t777 - t926;
t820 = t700 * t775 + t733 * t778;
t614 = -t765 * t696 + t769 * t820;
t659 = t700 * t778 - t733 * t775;
t830 = t614 * t767 + t659 * t763;
t744 = -t921 - t959;
t698 = t744 * t777 - t925;
t702 = -t744 * t774 - t924;
t819 = t702 * t775 - t730 * t778;
t615 = -t765 * t698 + t769 * t819;
t660 = t702 * t778 + t730 * t775;
t829 = t615 * t767 + t660 * t763;
t745 = t756 - t959;
t695 = t745 * t774 + t924;
t699 = t745 * t777 - t925;
t808 = t699 * t775 - t777 * t892;
t626 = -t765 * t695 + t769 * t808;
t670 = t699 * t778 + t775 * t893;
t828 = t626 * t767 + t670 * t763;
t743 = -t921 + t959;
t697 = t743 * t777 + t926;
t701 = -t743 * t774 + t725;
t807 = t701 * t775 - t774 * t892;
t627 = -t765 * t697 + t769 * t807;
t671 = t701 * t778 + t774 * t894;
t827 = t627 * t767 + t671 * t763;
t693 = t869 * t777;
t705 = -t732 * t774 - t759 * t896;
t798 = t705 * t775 - t875;
t630 = t765 * t693 + t769 * t798;
t672 = t705 * t778 + t876;
t826 = t630 * t767 + t672 * t763;
t694 = t870 * t774;
t706 = t731 * t777 - t758 * t896;
t797 = t706 * t775 + t875;
t631 = -t765 * t694 + t769 * t797;
t673 = t706 * t778 - t876;
t825 = t631 * t767 + t673 * t763;
t605 = t661 * t767 - t662 * t763;
t606 = t661 * t763 + t662 * t767;
t816 = t734 * t775 + t737 * t778;
t681 = t816 * t769;
t823 = t681 * t767 + t686 * t763;
t724 = t901 * t896;
t812 = -qJDD(4) * t778 + t724 * t775;
t688 = t812 * t769;
t707 = qJDD(4) * t775 + t724 * t778;
t821 = t688 * t767 + t707 * t763;
t716 = t735 * t769;
t817 = t716 * t767 - t763 * t813;
t689 = t739 * t768 - t740 * t764;
t804 = (-t369 * t766 - t370 * t770) * qJ(2);
t803 = (-t383 * t766 - t384 * t770) * qJ(2);
t802 = (-t389 * t766 - t390 * t770) * qJ(2);
t406 = -pkin(5) * t492 + qJ(6) * t454;
t328 = -pkin(10) * t385 - qJ(6) * t940 - t406 * t773;
t351 = -pkin(4) * t385 - pkin(5) * t451;
t294 = -pkin(9) * t362 + t328 * t777 - t351 * t774;
t303 = -pkin(3) * t362 + pkin(4) * t492 - pkin(10) * t386 + qJ(6) * t941 - t406 * t776;
t796 = pkin(8) * t327 + t294 * t775 + t303 * t778;
t323 = -pkin(9) * t379 + (-pkin(10) * t777 + t956) * t407;
t333 = -pkin(3) * t379 + pkin(4) * t520 - pkin(10) * t408;
t795 = pkin(8) * t339 + t323 * t775 + t333 * t778;
t516 = -pkin(3) * t696 + t535;
t556 = -pkin(9) * t696 + t937;
t794 = pkin(8) * t659 + t516 * t778 + t556 * t775;
t517 = -pkin(3) * t698 + t536;
t557 = -pkin(9) * t698 + t936;
t793 = pkin(8) * t660 + t517 * t778 + t557 * t775;
t290 = -pkin(3) * t385 + pkin(9) * t363 + t328 * t774 + t351 * t777;
t306 = t769 * t362 + t765 * t859;
t263 = -pkin(2) * t306 - t765 * t290 + t769 * t796;
t265 = t778 * t294 - t775 * t303 + (-t306 * t765 - t307 * t769) * pkin(8);
t291 = -t307 * t763 + t327 * t767;
t792 = qJ(2) * t291 + t263 * t767 + t265 * t763;
t308 = pkin(9) * t380 + (-pkin(3) + t874) * t407;
t318 = t769 * t379 + t765 * t855;
t270 = -pkin(2) * t318 - t765 * t308 + t769 * t795;
t278 = t778 * t323 - t775 * t333 + (-t318 * t765 - t319 * t769) * pkin(8);
t295 = -t319 * t763 + t339 * t767;
t791 = qJ(2) * t295 + t270 * t767 + t278 * t763;
t391 = -t469 * t769 + t765 * t849;
t445 = -pkin(3) * t574 + pkin(9) * t470;
t780 = pkin(8) * t444 - (-pkin(3) * t778 - pkin(9) * t775) * t469;
t317 = -pkin(2) * t391 - t765 * t445 + t769 * t780;
t321 = -(pkin(3) * t775 - pkin(9) * t778) * t469 + (-t391 * t765 - t392 * t769) * pkin(8);
t352 = -t392 * t763 + t444 * t767;
t790 = qJ(2) * t352 + t317 * t767 + t321 * t763;
t467 = pkin(3) * t737 + pkin(9) * t734 + t470;
t680 = t816 * t765;
t395 = -pkin(2) * t680 - t765 * t467 + t769 * t873;
t455 = t778 * t469 + (-t680 * t765 - t681 * t769) * pkin(8);
t632 = -t681 * t763 + t686 * t767;
t789 = qJ(2) * t632 + t395 * t767 + t455 * t763;
t549 = pkin(3) * t733 + pkin(9) * t700 - t936;
t612 = t769 * t696 + t765 * t820;
t424 = -pkin(2) * t612 - t765 * t549 + t769 * t794;
t442 = -t775 * t516 + t778 * t556 + (-t612 * t765 - t614 * t769) * pkin(8);
t572 = -t614 * t763 + t659 * t767;
t788 = qJ(2) * t572 + t424 * t767 + t442 * t763;
t548 = -pkin(3) * t730 + pkin(9) * t702 + t937;
t613 = t769 * t698 + t765 * t819;
t425 = -pkin(2) * t613 - t765 * t548 + t769 * t793;
t443 = -t775 * t517 + t778 * t557 + (-t613 * t765 - t615 * t769) * pkin(8);
t573 = -t615 * t763 + t660 * t767;
t787 = qJ(2) * t573 + t425 * t767 + t443 * t763;
t489 = -t506 * t765 - t769 * t623;
t426 = (-t489 * t765 - t490 * t769) * pkin(8);
t436 = -t490 * t763 + t507 * t767;
t437 = -pkin(2) * t489 + t507 * t953;
t786 = qJ(2) * t436 + t426 * t763 + t437 * t767;
t513 = pkin(2) * t713 + t765 * t581 + t769 * t872;
t587 = -t912 + (t713 * t765 + t715 * t769) * pkin(8);
t785 = qJ(2) * t677 + t513 * t767 + t587 * t763;
t714 = t735 * t765;
t514 = -pkin(2) * t714 + t765 * t580 + t769 * t871;
t586 = -t913 + (-t714 * t765 - t716 * t769) * pkin(8);
t678 = -t716 * t763 - t767 * t813;
t784 = qJ(2) * t678 + t514 * t767 + t586 * t763;
t712 = (-t765 * t767 * t770 - t766 * t769) * qJDD(3);
t687 = t812 * t765;
t648 = -t688 * t763 + t707 * t767;
t647 = -t766 * t714 + t770 * t817;
t644 = t770 * t714 + t766 * t817;
t629 = t769 * t694 + t765 * t797;
t628 = -t769 * t693 + t765 * t798;
t625 = t769 * t697 + t765 * t807;
t624 = t769 * t695 + t765 * t808;
t622 = t661 * t766 + t767 * t915;
t621 = t662 * t766 - t763 * t915;
t607 = t769 * t682 + t765 * t822;
t604 = t606 * t770;
t603 = -t766 * t687 + t770 * t821;
t597 = -t647 * t764 + t678 * t768;
t595 = t647 * t768 + t678 * t764;
t593 = -t766 * t680 + t770 * t823;
t592 = t770 * t680 + t766 * t823;
t591 = -t631 * t763 + t673 * t767;
t590 = -t630 * t763 + t672 * t767;
t589 = -t627 * t763 + t671 * t767;
t588 = -t626 * t763 + t670 * t767;
t577 = -t605 * t770 + t919;
t576 = -t605 * t766 - t915;
t561 = -t608 * t763 + t655 * t767;
t544 = -t593 * t764 + t632 * t768;
t543 = t593 * t768 + t632 * t764;
t542 = -t766 * t629 + t770 * t825;
t541 = -t766 * t628 + t770 * t826;
t540 = -t766 * t625 + t770 * t827;
t539 = -t766 * t624 + t770 * t828;
t537 = t769 * t637 + t765 * t824;
t531 = -t766 * t613 + t770 * t829;
t530 = -t766 * t612 + t770 * t830;
t529 = t770 * t613 + t766 * t829;
t528 = t770 * t612 + t766 * t830;
t527 = -pkin(1) * t576 + t606 * t942;
t526 = -t577 * t764 + t606 * t768;
t525 = t577 * t768 + t606 * t764;
t512 = pkin(2) * t716 - t769 * t580 + t765 * t871;
t511 = -pkin(2) * t715 - t769 * t581 + t765 * t872;
t508 = -t766 * t607 + t770 * t831;
t505 = (-t576 * t766 - t577 * t770) * qJ(2);
t503 = t507 * t765;
t500 = t769 * t600 + t765 * t832;
t499 = t769 * t599 + t765 * t833;
t491 = -t538 * t763 + t578 * t767;
t485 = t769 * t568 + t765 * t834;
t484 = t769 * t567 + t765 * t835;
t483 = -t531 * t764 + t573 * t768;
t482 = -t530 * t764 + t572 * t768;
t481 = t531 * t768 + t573 * t764;
t480 = t530 * t768 + t572 * t764;
t466 = -t502 * t763 + t555 * t767;
t465 = -t501 * t763 + t554 * t767;
t463 = t769 * t558 + t765 * t838;
t457 = -t763 * t513 + t767 * t587 + (t645 * t770 + t766 * t960) * qJ(2);
t456 = -t763 * t514 + t767 * t586 + (-t644 * t766 - t647 * t770) * qJ(2);
t452 = -t504 * t763 + t506 * t767;
t448 = -t766 * t537 + t770 * t840;
t447 = t770 * t537 + t766 * t840;
t441 = -t487 * t763 + t534 * t767;
t440 = -t486 * t763 + t533 * t767;
t438 = pkin(2) * t490 + t507 * t954;
t423 = pkin(2) * t615 + t769 * t548 + t765 * t793;
t422 = pkin(2) * t614 + t769 * t549 + t765 * t794;
t421 = pkin(1) * t960 - t766 * t511 + t770 * t785;
t420 = -pkin(1) * t644 - t766 * t512 + t770 * t784;
t419 = -t766 * t500 + t770 * t842;
t418 = -t766 * t499 + t770 * t843;
t417 = t770 * t500 + t766 * t842;
t416 = t770 * t499 + t766 * t843;
t414 = -t464 * t763 + t498 * t767;
t411 = -t766 * t503 + t770 * t841;
t405 = -t448 * t764 + t491 * t768;
t404 = t448 * t768 + t491 * t764;
t403 = -t766 * t489 + t770 * t844;
t402 = t770 * t489 + t766 * t844;
t401 = -t766 * t485 + t770 * t845;
t400 = -t766 * t484 + t770 * t846;
t399 = t770 * t485 + t766 * t845;
t398 = t770 * t484 + t766 * t846;
t394 = pkin(2) * t681 + t769 * t467 + t765 * t873;
t377 = -t419 * t764 + t466 * t768;
t376 = -t418 * t764 + t465 * t768;
t375 = t419 * t768 + t466 * t764;
t374 = t418 * t768 + t465 * t764;
t373 = -t766 * t463 + t770 * t850;
t372 = t770 * t463 + t766 * t850;
t359 = -t401 * t764 + t441 * t768;
t358 = -t400 * t764 + t440 * t768;
t357 = t401 * t768 + t441 * t764;
t356 = t400 * t768 + t440 * t764;
t354 = -t403 * t764 + t436 * t768;
t353 = t403 * t768 + t436 * t764;
t350 = -t390 * t764 + t435 * t768;
t348 = qJ(1) * t350;
t346 = -t384 * t764 + t430 * t768;
t344 = qJ(1) * t346;
t343 = -t763 * t395 + t767 * t455 + (-t592 * t766 - t593 * t770) * qJ(2);
t342 = -t763 * t425 + t767 * t443 + (-t529 * t766 - t531 * t770) * qJ(2);
t341 = -t763 * t424 + t767 * t442 + (-t528 * t766 - t530 * t770) * qJ(2);
t338 = -t373 * t764 + t414 * t768;
t337 = t373 * t768 + t414 * t764;
t336 = -t370 * t764 + t412 * t768;
t334 = qJ(1) * t336;
t330 = -pkin(1) * t529 - t766 * t423 + t770 * t787;
t329 = -pkin(1) * t528 - t766 * t422 + t770 * t788;
t326 = -t766 * t391 + t770 * t854;
t325 = t770 * t391 + t766 * t854;
t324 = -pkin(1) * t592 - t766 * t394 + t770 * t789;
t316 = pkin(2) * t392 + t769 * t445 + t765 * t780;
t314 = t769 * t387 + t765 * t852 + t906;
t313 = t767 * t426 - t763 * t437 + (-t402 * t766 - t403 * t770) * qJ(2);
t311 = t769 * t381 + t765 * t853 + t907;
t309 = -pkin(1) * t402 - t766 * t438 + t770 * t786;
t304 = t769 * t361 + t765 * t857 + t906;
t302 = -t326 * t764 + t352 * t768;
t301 = t326 * t768 + t352 * t764;
t298 = t769 * t360 + t765 * t858 + t907;
t296 = t769 * t364 + t765 * t856 + t908;
t292 = t769 * t340 + t765 * t860 + t908;
t289 = -t763 * t315 + t767 * t332 + t802;
t288 = -t763 * t312 + t767 * t331 + t803;
t287 = -t766 * t318 + t770 * t861;
t286 = t770 * t318 + t766 * t861;
t285 = -t763 * t305 + t767 * t322 + t802;
t284 = -t763 * t299 + t767 * t320 + t803;
t283 = -t766 * t306 + t770 * t864;
t282 = t770 * t306 + t766 * t864;
t281 = -t766 * t314 + t770 * t862 + t883;
t280 = -t766 * t311 + t770 * t863 + t884;
t279 = -t763 * t297 + t767 * t310 + t804;
t277 = -t763 * t317 + t767 * t321 + (-t325 * t766 - t326 * t770) * qJ(2);
t276 = -t766 * t304 + t770 * t865 + t883;
t275 = -t763 * t293 + t767 * t300 + t804;
t274 = -t287 * t764 + t295 * t768;
t273 = t287 * t768 + t295 * t764;
t272 = -pkin(1) * t325 - t766 * t316 + t770 * t790;
t271 = -t766 * t298 + t770 * t866 + t884;
t269 = pkin(2) * t319 + t769 * t308 + t765 * t795;
t268 = -t766 * t296 + t770 * t867 + t885;
t267 = -t283 * t764 + t291 * t768;
t266 = t283 * t768 + t291 * t764;
t264 = -t766 * t292 + t770 * t868 + t885;
t262 = pkin(2) * t307 + t769 * t290 + t765 * t796;
t261 = -t763 * t270 + t767 * t278 + (-t286 * t766 - t287 * t770) * qJ(2);
t260 = -pkin(1) * t286 - t766 * t269 + t770 * t791;
t259 = -t763 * t263 + t767 * t265 + (-t282 * t766 - t283 * t770) * qJ(2);
t258 = -pkin(1) * t282 - t766 * t262 + t770 * t792;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t690, 0, 0, 0, 0, 0, 0, 0, 0, 0, t526, 0, 0, 0, 0, 0, 0, t597, t596, 0, t354, 0, 0, 0, 0, 0, 0, t482, t483, t544, t302, 0, 0, 0, 0, 0, 0, t346, t350, t336, t274, 0, 0, 0, 0, 0, 0, t346, t350, t336, t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t689, 0, 0, 0, 0, 0, 0, 0, 0, 0, t525, 0, 0, 0, 0, 0, 0, t595, -t966, 0, t353, 0, 0, 0, 0, 0, 0, t480, t481, t543, t301, 0, 0, 0, 0, 0, 0, t345, t349, t335, t273, 0, 0, 0, 0, 0, 0, t345, t349, t335, t266; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t760, 0, 0, 0, 0, 0, 0, 0, 0, 0, t576, 0, 0, 0, 0, 0, 0, t644, -t960, 0, t402, 0, 0, 0, 0, 0, 0, t528, t529, t592, t325, 0, 0, 0, 0, 0, 0, t383, t389, t369, t286, 0, 0, 0, 0, 0, 0, t383, t389, t369, t282; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t920, -t916, -t689, -qJ(1) * t689, 0, 0, 0, 0, 0, 0, -t622 * t764 - t763 * t928, -t621 * t764 - t767 * t928, -t604 * t764 + t605 * t768, -qJ(1) * t525 + t505 * t768 - t527 * t764, 0, 0, -t596, 0, t597, -t712 * t764 + t768 * t887, -qJ(1) * t595 - t420 * t764 + t456 * t768, qJ(1) * t966 - t421 * t764 + t457 * t768, -t411 * t764 + t452 * t768, -qJ(1) * t353 - t309 * t764 + t313 * t768, -t542 * t764 + t591 * t768, -t508 * t764 + t561 * t768, -t540 * t764 + t589 * t768, -t541 * t764 + t590 * t768, -t539 * t764 + t588 * t768, -t603 * t764 + t648 * t768, -qJ(1) * t480 - t329 * t764 + t341 * t768, -qJ(1) * t481 - t330 * t764 + t342 * t768, -qJ(1) * t543 - t324 * t764 + t343 * t768, -qJ(1) * t301 - t272 * t764 + t277 * t768, t377, t338, t358, t376, t359, t405, -t280 * t764 + t288 * t768 - t945, -t281 * t764 + t289 * t768 - t944, -t268 * t764 + t279 * t768 - t946, -qJ(1) * t273 - t260 * t764 + t261 * t768, t377, t338, t358, t376, t359, t405, -t271 * t764 + t284 * t768 - t945, -t276 * t764 + t285 * t768 - t944, -t264 * t764 + t275 * t768 - t946, -qJ(1) * t266 - t258 * t764 + t259 * t768; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t916, -t920, t690, qJ(1) * t690, 0, 0, 0, 0, 0, 0, t622 * t768 - t763 * t930, t621 * t768 - t767 * t930, t604 * t768 + t605 * t764, qJ(1) * t526 + t505 * t764 + t527 * t768, 0, 0, t966, 0, t595, t712 * t768 + t764 * t887, qJ(1) * t597 + t420 * t768 + t456 * t764, qJ(1) * t596 + t421 * t768 + t457 * t764, t411 * t768 + t452 * t764, qJ(1) * t354 + t309 * t768 + t313 * t764, t542 * t768 + t591 * t764, t508 * t768 + t561 * t764, t540 * t768 + t589 * t764, t541 * t768 + t590 * t764, t539 * t768 + t588 * t764, t603 * t768 + t648 * t764, qJ(1) * t482 + t329 * t768 + t341 * t764, qJ(1) * t483 + t330 * t768 + t342 * t764, qJ(1) * t544 + t324 * t768 + t343 * t764, qJ(1) * t302 + t272 * t768 + t277 * t764, t375, t337, t356, t374, t357, t404, t280 * t768 + t288 * t764 + t344, t281 * t768 + t289 * t764 + t348, t268 * t768 + t279 * t764 + t334, qJ(1) * t274 + t260 * t768 + t261 * t764, t375, t337, t356, t374, t357, t404, t271 * t768 + t284 * t764 + t344, t276 * t768 + t285 * t764 + t348, t264 * t768 + t275 * t764 + t334, qJ(1) * t267 + t258 * t768 + t259 * t764; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t739, t740, 0, 0, 0, 0, 0, 0, 0, 0, -t661 * t770 + t704 * t917, -t662 * t770 - t763 * t919, t606 * t766, pkin(1) * t577 + t606 * t943, 0, 0, t960, 0, t644, (-t765 * t917 + t769 * t770) * qJDD(3), pkin(1) * t647 + t770 * t512 + t766 * t784, -pkin(1) * t645 + t770 * t511 + t766 * t785, t770 * t503 + t766 * t841, pkin(1) * t403 + t770 * t438 + t766 * t786, t770 * t629 + t766 * t825, t770 * t607 + t766 * t831, t770 * t625 + t766 * t827, t770 * t628 + t766 * t826, t770 * t624 + t766 * t828, t770 * t687 + t766 * t821, pkin(1) * t530 + t770 * t422 + t766 * t788, pkin(1) * t531 + t770 * t423 + t766 * t787, pkin(1) * t593 + t770 * t394 + t766 * t789, pkin(1) * t326 + t770 * t316 + t766 * t790, t417, t372, t398, t416, t399, t447, t770 * t311 + t766 * t863 + t910, t770 * t314 + t766 * t862 + t909, t770 * t296 + t766 * t867 + t911, pkin(1) * t287 + t770 * t269 + t766 * t791, t417, t372, t398, t416, t399, t447, t770 * t298 + t766 * t866 + t910, t770 * t304 + t766 * t865 + t909, t770 * t292 + t766 * t868 + t911, pkin(1) * t283 + t770 * t262 + t766 * t792;];
tauB_reg  = t1;
