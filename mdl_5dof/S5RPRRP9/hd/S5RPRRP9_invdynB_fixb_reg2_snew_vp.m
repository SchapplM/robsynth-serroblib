% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRRP9_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:22
% EndTime: 2019-12-31 18:49:39
% DurationCPUTime: 16.99s
% Computational Cost: add. (41672->525), mult. (102606->738), div. (0->0), fcn. (76742->8), ass. (0->385)
t713 = cos(pkin(8));
t718 = cos(qJ(3));
t712 = sin(pkin(8));
t715 = sin(qJ(3));
t784 = t712 * t715;
t680 = (t713 * t718 - t784) * qJD(1);
t732 = t712 * t718 + t713 * t715;
t681 = t732 * qJD(1);
t714 = sin(qJ(4));
t717 = cos(qJ(4));
t640 = t714 * t680 + t717 * t681;
t637 = t640 ^ 2;
t711 = qJD(3) + qJD(4);
t798 = t711 ^ 2;
t583 = t798 + t637;
t638 = -t717 * t680 + t714 * t681;
t593 = t640 * t638;
t708 = qJDD(3) + qJDD(4);
t811 = t593 + t708;
t778 = t714 * t811;
t513 = t717 * t583 + t778;
t767 = t717 * t811;
t544 = t714 * t583 - t767;
t460 = t718 * t513 - t715 * t544;
t484 = t715 * t513 + t718 * t544;
t438 = t712 * t460 + t713 * t484;
t716 = sin(qJ(1));
t719 = cos(qJ(1));
t755 = qJDD(1) * t713;
t678 = qJDD(1) * t784 - t718 * t755;
t758 = t681 * qJD(3);
t649 = -t678 - t758;
t759 = t680 * qJD(3);
t806 = t732 * qJDD(1);
t651 = t806 + t759;
t726 = t638 * qJD(4) - t714 * t649 - t717 * t651;
t790 = t711 * t638;
t819 = -t726 - t790;
t407 = t716 * t438 - t719 * t819;
t901 = pkin(5) * t407;
t409 = t719 * t438 + t716 * t819;
t900 = pkin(5) * t409;
t413 = t713 * t460 - t712 * t484;
t899 = qJ(2) * t413;
t898 = -pkin(1) * t413 - pkin(2) * t460 - pkin(3) * t513;
t897 = pkin(1) * t819 - qJ(2) * t438;
t743 = -t717 * t649 + t714 * t651;
t530 = (qJD(4) + t711) * t640 + t743;
t473 = -t714 * t530 + t717 * t819;
t781 = t714 * t819;
t475 = t717 * t530 + t781;
t423 = -t718 * t473 + t715 * t475;
t427 = t715 * t473 + t718 * t475;
t386 = t712 * t423 - t713 * t427;
t799 = t638 ^ 2;
t590 = -t799 + t637;
t896 = t716 * t386 - t719 * t590;
t895 = t719 * t386 + t716 * t590;
t623 = t799 - t798;
t549 = t714 * t623 + t767;
t553 = t717 * t623 - t778;
t488 = t718 * t549 + t715 * t553;
t493 = t715 * t549 - t718 * t553;
t443 = t712 * t488 + t713 * t493;
t531 = (qJD(4) - t711) * t640 + t743;
t894 = t716 * t443 - t719 * t531;
t893 = t719 * t443 + t716 * t531;
t565 = -t799 - t637;
t808 = -t790 + t726;
t836 = -t717 * t531 - t714 * t808;
t837 = -t714 * t531 + t717 * t808;
t853 = t715 * t836 + t718 * t837;
t854 = -t715 * t837 + t718 * t836;
t865 = -t712 * t853 + t713 * t854;
t880 = t716 * t565 + t719 * t865;
t892 = pkin(5) * t880;
t882 = -t719 * t565 + t716 * t865;
t891 = pkin(5) * t882;
t890 = pkin(6) * t460;
t888 = pkin(2) * t819 - pkin(6) * t484;
t887 = t713 * t423 + t712 * t427;
t886 = t713 * t488 - t712 * t493;
t866 = t712 * t854 + t713 * t853;
t885 = qJ(2) * t866;
t366 = -pkin(1) * t866 - pkin(2) * t853 - pkin(3) * t837;
t884 = -pkin(1) * t565 + qJ(2) * t865;
t624 = -t637 + t798;
t812 = -t593 + t708;
t777 = t714 * t812;
t838 = t717 * t624 + t777;
t581 = t717 * t812;
t839 = -t714 * t624 + t581;
t851 = -t715 * t838 + t718 * t839;
t852 = t715 * t839 + t718 * t838;
t868 = -t712 * t852 + t713 * t851;
t883 = t716 * t868 + t719 * t808;
t881 = -t716 * t808 + t719 * t868;
t878 = pkin(6) * t853;
t877 = pkin(7) * t513;
t876 = pkin(7) * t544;
t807 = -t798 - t799;
t817 = t717 * t807 - t777;
t818 = t714 * t807 + t581;
t830 = t715 * t817 + t718 * t818;
t831 = -t715 * t818 + t718 * t817;
t855 = -t712 * t830 + t713 * t831;
t875 = qJ(2) * t855;
t856 = t712 * t831 + t713 * t830;
t874 = qJ(2) * t856;
t873 = t716 * t855;
t872 = t719 * t855;
t871 = -pkin(1) * t856 - pkin(2) * t830 - pkin(3) * t818;
t869 = -pkin(2) * t565 + pkin(6) * t854;
t867 = t712 * t851 + t713 * t852;
t863 = pkin(6) * t830;
t862 = pkin(6) * t831;
t861 = pkin(7) * t837;
t857 = -pkin(3) * t565 + pkin(7) * t836;
t846 = pkin(7) * t817;
t845 = pkin(7) * t818;
t840 = t819 * qJ(5);
t693 = t719 * g(1) + t716 * g(2);
t721 = qJD(1) ^ 2;
t820 = -t721 * pkin(1) + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t693;
t728 = (-t638 * t714 - t640 * t717) * t711;
t788 = t711 * t714;
t620 = t640 * t788;
t787 = t711 * t717;
t749 = t638 * t787;
t736 = t620 - t749;
t802 = -t715 * t728 + t718 * t736;
t803 = t715 * t736 + t718 * t728;
t814 = -t712 * t803 + t713 * t802;
t835 = -t719 * t708 + t716 * t814;
t748 = t719 * t593;
t568 = -t640 * qJD(4) - t743;
t731 = -t714 * t568 + t749;
t737 = t717 * t568 + t638 * t788;
t800 = t715 * t731 + t718 * t737;
t801 = -t715 * t737 + t718 * t731;
t815 = -t712 * t800 + t713 * t801;
t834 = t716 * t815 + t748;
t751 = t716 * t593;
t833 = t719 * t815 - t751;
t832 = t716 * t708 + t719 * t814;
t829 = 2 * qJD(5);
t654 = t680 * t681;
t805 = qJDD(3) + t654;
t826 = t715 * t805;
t823 = t718 * t805;
t816 = t712 * t801 + t713 * t800;
t813 = t712 * t802 + t713 * t803;
t709 = t712 ^ 2;
t710 = t713 ^ 2;
t810 = t709 + t710;
t794 = t713 * g(3);
t797 = pkin(2) * t713;
t622 = -t794 + (-pkin(6) * qJDD(1) + t721 * t797 - t820) * t712;
t656 = -t712 * g(3) + t820 * t713;
t704 = t710 * t721;
t631 = -pkin(2) * t704 + pkin(6) * t755 + t656;
t572 = -t718 * t622 + t715 * t631;
t510 = (-t651 + t759) * pkin(7) + t805 * pkin(3) - t572;
t573 = t715 * t622 + t718 * t631;
t676 = t680 ^ 2;
t734 = qJD(3) * pkin(3) - t681 * pkin(7);
t519 = -t676 * pkin(3) + t649 * pkin(7) - qJD(3) * t734 + t573;
t459 = t714 * t510 + t717 * t519;
t589 = t638 * pkin(4) - t640 * qJ(5);
t733 = t708 * qJ(5) - t638 * t589 + t711 * t829 + t459;
t804 = t721 * t810;
t677 = t681 ^ 2;
t796 = pkin(4) * t717;
t795 = t568 * pkin(4);
t793 = qJ(5) * t717;
t792 = qJDD(1) * pkin(1);
t791 = t709 * t721;
t789 = t711 * t640;
t502 = -t718 * t572 + t715 * t573;
t786 = t712 * t502;
t785 = t712 * t713;
t783 = t713 * t502;
t692 = t716 * g(1) - t719 * g(2);
t735 = -qJDD(2) + t692;
t644 = (pkin(1) + t797) * qJDD(1) + (t810 * pkin(6) + qJ(2)) * t721 + t735;
t559 = t649 * pkin(3) + t676 * pkin(7) - t681 * t734 + t644;
t779 = t714 * t559;
t458 = -t717 * t510 + t714 * t519;
t411 = -t717 * t458 + t714 * t459;
t775 = t715 * t411;
t774 = t715 * t644;
t646 = qJDD(3) - t654;
t773 = t715 * t646;
t674 = t721 * qJ(2) + t735 + t792;
t771 = t716 * t674;
t768 = t717 * t559;
t766 = t718 * t411;
t765 = t718 * t644;
t764 = t718 * t646;
t762 = t719 * t674;
t761 = -t565 - t798;
t754 = t716 * qJDD(1);
t753 = t719 * qJDD(1);
t750 = t716 * t654;
t747 = t719 * t654;
t746 = -qJ(5) * t714 - pkin(3);
t745 = t674 + t792;
t412 = t714 * t458 + t717 * t459;
t503 = t715 * t572 + t718 * t573;
t655 = t820 * t712 + t794;
t602 = t712 * t655 + t713 * t656;
t664 = -t716 * t692 - t719 * t693;
t741 = t640 * t589 + qJDD(5) + t458;
t691 = -t716 * t721 + t753;
t740 = -pkin(5) * t691 - t716 * g(3);
t524 = t640 * t787 - t714 * t726;
t525 = -t717 * t726 - t620;
t466 = t718 * t524 + t715 * t525;
t469 = -t715 * t524 + t718 * t525;
t422 = -t712 * t466 + t713 * t469;
t739 = t716 * t422 - t748;
t738 = t719 * t422 + t751;
t601 = t713 * t655 - t712 * t656;
t663 = t719 * t692 - t716 * t693;
t690 = t719 * t721 + t754;
t684 = t713 * t804;
t660 = -t716 * t684 + t713 * t753;
t730 = t719 * t684 + t713 * t754;
t729 = -t708 * pkin(4) + t741;
t725 = -pkin(4) * t789 + t640 * t829 + t559;
t724 = t725 + t840;
t720 = qJD(3) ^ 2;
t703 = t710 * qJDD(1);
t702 = t709 * qJDD(1);
t689 = t704 - t791;
t688 = t704 + t791;
t687 = t703 - t702;
t686 = t703 + t702;
t683 = t712 * t804;
t675 = -pkin(5) * t690 + t719 * g(3);
t669 = -t677 - t720;
t668 = -t677 + t720;
t667 = t676 - t720;
t666 = t691 * t785;
t665 = t690 * t785;
t661 = t719 * t683 + t712 * t754;
t659 = t716 * t683 - t712 * t753;
t658 = t719 * t686 - t716 * t688;
t657 = t716 * t686 + t719 * t688;
t653 = -t677 + t676;
t650 = t806 + 0.2e1 * t759;
t648 = t678 + 0.2e1 * t758;
t643 = -t720 - t676;
t633 = (t680 * t718 + t681 * t715) * qJD(3);
t632 = (t680 * t715 - t681 * t718) * qJD(3);
t614 = -t676 - t677;
t612 = t718 * t651 - t715 * t758;
t611 = t715 * t651 + t718 * t758;
t610 = -t715 * t649 - t718 * t759;
t609 = t718 * t649 - t715 * t759;
t608 = -t715 * t669 - t764;
t607 = -t715 * t668 + t823;
t606 = t718 * t667 - t773;
t605 = t718 * t669 - t773;
t604 = t718 * t668 + t826;
t603 = t715 * t667 + t764;
t599 = -t718 * t648 - t715 * t650;
t598 = -t678 * t718 + t715 * t806;
t597 = -t715 * t648 + t718 * t650;
t596 = -t678 * t715 - t718 * t806;
t595 = t718 * t643 - t826;
t594 = t715 * t643 + t823;
t580 = -t712 * t632 + t713 * t633;
t575 = t719 * t602 - t771;
t574 = t716 * t602 + t762;
t570 = -pkin(6) * t605 - t765;
t567 = -pkin(6) * t594 - t774;
t561 = -t712 * t611 + t713 * t612;
t560 = -t712 * t609 + t713 * t610;
t558 = -t712 * t605 + t713 * t608;
t557 = -t712 * t604 + t713 * t607;
t556 = -t712 * t603 + t713 * t606;
t555 = t713 * t605 + t712 * t608;
t546 = -pkin(2) * t650 + pkin(6) * t608 - t774;
t541 = -pkin(2) * t648 + pkin(6) * t595 + t765;
t540 = -t712 * t597 + t713 * t599;
t539 = -t712 * t596 + t713 * t598;
t538 = t713 * t596 + t712 * t598;
t529 = -t568 + t789;
t528 = -t712 * t594 + t713 * t595;
t527 = t713 * t594 + t712 * t595;
t518 = t719 * t558 + t716 * t650;
t517 = t716 * t558 - t719 * t650;
t501 = t719 * t528 + t716 * t648;
t500 = t716 * t528 - t719 * t648;
t499 = t719 * t539 + t716 * t614;
t498 = t716 * t539 - t719 * t614;
t497 = pkin(2) * t644 + pkin(6) * t503;
t496 = -pkin(1) * t538 - pkin(2) * t596;
t495 = -t768 + t877;
t494 = -pkin(6) * t596 - t502;
t481 = -pkin(1) * t555 - pkin(2) * t605 + t573;
t480 = -t779 - t845;
t479 = -pkin(2) * t614 + pkin(6) * t598 + t503;
t470 = -pkin(1) * t527 - pkin(2) * t594 + t572;
t456 = -qJ(2) * t555 - t712 * t546 + t713 * t570;
t453 = t713 * t503 - t786;
t452 = t712 * t503 + t783;
t451 = -pkin(3) * t819 - t779 + t876;
t450 = -qJ(2) * t527 - t712 * t541 + t713 * t567;
t449 = t719 * t453 - t716 * t644;
t448 = t716 * t453 + t719 * t644;
t447 = -pkin(3) * t530 + t768 + t846;
t446 = t724 + t795;
t445 = qJ(5) * t798 - t729;
t444 = -pkin(4) * t798 + t733;
t435 = -pkin(1) * t452 - pkin(2) * t502;
t434 = t761 * qJ(5) + t729;
t433 = t761 * pkin(4) + t733;
t432 = (-t529 + t568) * pkin(4) + t724;
t431 = t725 + t795 + 0.2e1 * t840;
t419 = t713 * t466 + t712 * t469;
t410 = t716 * t529 + t872;
t408 = -t719 * t529 + t873;
t406 = -qJ(2) * t538 - t712 * t479 + t713 * t494;
t405 = -pkin(6) * t783 - qJ(2) * t452 - t712 * t497;
t404 = pkin(3) * t559 + pkin(7) * t412;
t403 = t716 * t530 + t872;
t401 = -t719 * t530 + t873;
t399 = -t714 * t432 - t529 * t793 - t845;
t398 = -t715 * t451 + t718 * t495 + t890;
t397 = -pkin(4) * t781 + t717 * t431 - t877;
t396 = t717 * t444 - t714 * t445;
t395 = t714 * t444 + t717 * t445;
t394 = -t411 - t861;
t393 = t717 * t432 + t746 * t529 + t846;
t392 = -t715 * t447 + t718 * t480 - t863;
t391 = t718 * t451 + t715 * t495 - t888;
t390 = -t876 + t714 * t431 + (pkin(3) + t796) * t819;
t389 = t412 + t857;
t388 = -pkin(2) * t530 + t718 * t447 + t715 * t480 + t862;
t381 = t459 - t898;
t376 = t718 * t412 - t775;
t375 = t715 * t412 + t766;
t374 = -t714 * t433 + t717 * t434 - t861;
t373 = t458 + t871;
t372 = (-t807 - t798) * qJ(5) + (-t812 - t708) * pkin(4) + t741 + t871;
t371 = t717 * t433 + t714 * t434 + t857;
t370 = -qJ(5) * t811 + (-t583 + t798) * pkin(4) - t733 + t898;
t369 = -pkin(7) * t395 + (-pkin(4) * t714 + t793) * t446;
t368 = -t715 * t395 + t718 * t396;
t367 = t718 * t395 + t715 * t396;
t365 = -t715 * t393 + t718 * t399 - t863;
t364 = pkin(7) * t396 + (-t746 + t796) * t446;
t363 = -t715 * t390 + t718 * t397 - t890;
t362 = -pkin(2) * t529 + t718 * t393 + t715 * t399 + t862;
t361 = -pkin(4) * t808 + qJ(5) * t531 + t366;
t360 = t718 * t390 + t715 * t397 + t888;
t359 = -t712 * t391 + t713 * t398 + t899;
t358 = -t715 * t389 + t718 * t394 - t878;
t357 = -t712 * t375 + t713 * t376;
t356 = t713 * t375 + t712 * t376;
t355 = t718 * t389 + t715 * t394 + t869;
t354 = -pkin(6) * t375 - pkin(7) * t766 - t715 * t404;
t353 = t719 * t357 - t716 * t559;
t352 = t716 * t357 + t719 * t559;
t351 = -t712 * t388 + t713 * t392 - t874;
t350 = pkin(2) * t559 + pkin(6) * t376 - pkin(7) * t775 + t718 * t404;
t349 = -t715 * t371 + t718 * t374 - t878;
t348 = -t712 * t367 + t713 * t368;
t347 = t713 * t367 + t712 * t368;
t346 = t718 * t371 + t715 * t374 + t869;
t345 = t719 * t348 - t716 * t446;
t344 = t716 * t348 + t719 * t446;
t343 = -pkin(1) * t356 - pkin(2) * t375 - pkin(3) * t411;
t342 = -t712 * t362 + t713 * t365 - t874;
t341 = -t712 * t360 + t713 * t363 - t899;
t340 = -pkin(6) * t367 - t715 * t364 + t718 * t369;
t339 = pkin(2) * t446 + pkin(6) * t368 + t718 * t364 + t715 * t369;
t338 = -t712 * t355 + t713 * t358 - t885;
t337 = -pkin(1) * t347 - pkin(2) * t367 - pkin(3) * t395 - pkin(4) * t445 - qJ(5) * t444;
t336 = -qJ(2) * t356 - t712 * t350 + t713 * t354;
t335 = -t712 * t346 + t713 * t349 - t885;
t334 = -qJ(2) * t347 - t712 * t339 + t713 * t340;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t690, -t691, 0, t664, 0, 0, 0, 0, 0, 0, -t730, t661, t658, t575, 0, 0, 0, 0, 0, 0, t501, t518, t499, t449, 0, 0, 0, 0, 0, 0, t403, t409, t880, t353, 0, 0, 0, 0, 0, 0, t410, t880, -t409, t345; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t691, -t690, 0, t663, 0, 0, 0, 0, 0, 0, t660, t659, t657, t574, 0, 0, 0, 0, 0, 0, t500, t517, t498, t448, 0, 0, 0, 0, 0, 0, t401, t407, t882, t352, 0, 0, 0, 0, 0, 0, t408, t882, -t407, t344; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t601, 0, 0, 0, 0, 0, 0, t527, t555, t538, t452, 0, 0, 0, 0, 0, 0, t856, -t413, t866, t356, 0, 0, 0, 0, 0, 0, t856, t866, t413, t347; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t691, 0, -t690, 0, t740, -t675, -t663, -pkin(5) * t663, t666, t719 * t687 - t716 * t689, t661, -t666, t730, 0, -pkin(5) * t660 - t716 * t655 - t712 * t762, -pkin(5) * t659 - t716 * t656 - t713 * t762, -pkin(5) * t657 + t719 * t601, -pkin(5) * t574 - (pkin(1) * t716 - qJ(2) * t719) * t601, t719 * t561 - t750, t719 * t540 - t716 * t653, t719 * t557 + t716 * t806, t719 * t560 + t750, t719 * t556 - t716 * t678, t716 * qJDD(3) + t719 * t580, -pkin(5) * t500 + t719 * t450 - t716 * t470, -pkin(5) * t517 + t719 * t456 - t716 * t481, -pkin(5) * t498 + t719 * t406 - t716 * t496, -pkin(5) * t448 + t719 * t405 - t716 * t435, t738, t895, t881, t833, -t893, t832, -pkin(5) * t401 + t719 * t351 - t716 * t373, t719 * t359 - t716 * t381 - t901, t719 * t338 - t716 * t366 - t891, -pkin(5) * t352 + t719 * t336 - t716 * t343, t738, t881, -t895, t832, t893, t833, -pkin(5) * t408 + t719 * t342 - t716 * t372, t719 * t335 - t716 * t361 - t891, t719 * t341 - t716 * t370 + t901, -pkin(5) * t344 + t719 * t334 - t716 * t337; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t690, 0, t691, 0, t675, t740, t664, pkin(5) * t664, t665, t716 * t687 + t719 * t689, t659, -t665, -t660, 0, -pkin(5) * t730 + t719 * t655 - t712 * t771, pkin(5) * t661 + t719 * t656 - t713 * t771, pkin(5) * t658 + t716 * t601, pkin(5) * t575 - (-pkin(1) * t719 - qJ(2) * t716) * t601, t716 * t561 + t747, t716 * t540 + t719 * t653, t716 * t557 - t719 * t806, t716 * t560 - t747, t716 * t556 + t719 * t678, -t719 * qJDD(3) + t716 * t580, pkin(5) * t501 + t716 * t450 + t719 * t470, pkin(5) * t518 + t716 * t456 + t719 * t481, pkin(5) * t499 + t716 * t406 + t719 * t496, pkin(5) * t449 + t716 * t405 + t719 * t435, t739, t896, t883, t834, -t894, t835, pkin(5) * t403 + t716 * t351 + t719 * t373, t716 * t359 + t719 * t381 + t900, t716 * t338 + t719 * t366 + t892, pkin(5) * t353 + t716 * t336 + t719 * t343, t739, t883, -t896, t835, t894, t834, pkin(5) * t410 + t716 * t342 + t719 * t372, t716 * t335 + t719 * t361 + t892, t716 * t341 + t719 * t370 - t900, pkin(5) * t345 + t716 * t334 + t719 * t337; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t692, t693, 0, 0, t702, 0.2e1 * t712 * t755, 0, t703, 0, 0, -qJ(2) * t684 + t713 * t745, qJ(2) * t683 - t712 * t745, pkin(1) * t688 + qJ(2) * t686 + t602, pkin(1) * t674 + qJ(2) * t602, t713 * t611 + t712 * t612, t713 * t597 + t712 * t599, t713 * t604 + t712 * t607, t713 * t609 + t712 * t610, t713 * t603 + t712 * t606, t713 * t632 + t712 * t633, -pkin(1) * t648 + qJ(2) * t528 + t713 * t541 + t712 * t567, -pkin(1) * t650 + qJ(2) * t558 + t713 * t546 + t712 * t570, -pkin(1) * t614 + qJ(2) * t539 + t713 * t479 + t712 * t494, pkin(1) * t644 - pkin(6) * t786 + qJ(2) * t453 + t713 * t497, t419, -t887, t867, t816, t886, t813, -pkin(1) * t530 + t713 * t388 + t712 * t392 + t875, t713 * t391 + t712 * t398 - t897, t713 * t355 + t712 * t358 + t884, pkin(1) * t559 + qJ(2) * t357 + t713 * t350 + t712 * t354, t419, t867, t887, t813, -t886, t816, -pkin(1) * t529 + t713 * t362 + t712 * t365 + t875, t713 * t346 + t712 * t349 + t884, t713 * t360 + t712 * t363 + t897, pkin(1) * t446 + qJ(2) * t348 + t713 * t339 + t712 * t340;];
tauB_reg = t1;