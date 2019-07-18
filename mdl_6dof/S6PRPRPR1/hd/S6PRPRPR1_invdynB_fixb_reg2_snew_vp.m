% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 22:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6PRPRPR1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynB_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:11:07
% EndTime: 2019-05-04 22:11:37
% DurationCPUTime: 29.85s
% Computational Cost: add. (178415->908), mult. (363372->1481), div. (0->0), fcn. (269099->14), ass. (0->645)
t831 = cos(qJ(2));
t818 = sin(pkin(10));
t822 = cos(pkin(10));
t786 = g(1) * t822 + g(2) * t818;
t828 = sin(qJ(2));
t785 = g(1) * t818 - t822 * g(2);
t819 = sin(pkin(6));
t823 = cos(pkin(6));
t920 = g(3) - qJDD(1);
t853 = t823 * t785 - t819 * t920;
t692 = -t831 * t786 + t828 * t853;
t833 = qJD(2) ^ 2;
t686 = -t833 * pkin(2) + t692;
t817 = sin(pkin(11));
t821 = cos(pkin(11));
t691 = -t828 * t786 - t831 * t853;
t835 = qJDD(2) * pkin(2) - t691;
t622 = t686 * t817 - t821 * t835;
t623 = t821 * t686 + t817 * t835;
t898 = t622 * t817 + t821 * t623;
t527 = t622 * t821 - t623 * t817;
t940 = t527 * t828;
t964 = t831 * t898 + t940;
t939 = t527 * t831;
t448 = -t828 * t898 + t939;
t816 = sin(pkin(12));
t820 = cos(pkin(12));
t830 = cos(qJ(4));
t917 = qJD(2) * t830;
t827 = sin(qJ(4));
t918 = qJD(2) * t827;
t766 = t816 * t918 - t820 * t917;
t768 = t816 * t917 + t820 * t918;
t718 = t768 * t766;
t949 = qJDD(4) - t718;
t963 = t816 * t949;
t962 = t820 * t949;
t826 = sin(qJ(6));
t829 = cos(qJ(6));
t734 = -t829 * qJD(4) + t768 * t826;
t736 = qJD(4) * t826 + t768 * t829;
t666 = t736 * t734;
t911 = qJD(2) * qJD(4);
t900 = t830 * t911;
t909 = qJDD(2) * t827;
t775 = t900 + t909;
t807 = t830 * qJDD(2);
t902 = t827 * t911;
t776 = t807 - t902;
t896 = t775 * t816 - t820 * t776;
t890 = qJDD(6) + t896;
t950 = -t666 + t890;
t961 = t826 * t950;
t960 = t829 * t950;
t915 = qJD(4) * t768;
t667 = t896 + t915;
t752 = t785 * t819 + t823 * t920;
t747 = -qJDD(3) + t752;
t778 = qJDD(2) * t817 + t821 * t833;
t690 = qJ(3) * t778 - t747 * t821;
t779 = qJDD(2) * t821 - t817 * t833;
t889 = -qJ(3) * t779 - t747 * t817;
t959 = t690 * t828 + t889 * t831;
t958 = t690 * t831 - t889 * t828;
t724 = t778 * t828 - t779 * t831;
t855 = t778 * t831 + t779 * t828;
t953 = t855 * t823;
t641 = t724 * t822 + t818 * t953;
t639 = t724 * t818 - t822 * t953;
t957 = t818 * t920;
t956 = t822 * t920;
t703 = t855 * t819;
t720 = t775 * t820 + t776 * t816;
t647 = -t734 * qJD(6) + t826 * qJDD(4) + t829 * t720;
t759 = qJD(6) + t766;
t685 = t759 * t734;
t601 = -t685 + t647;
t897 = -t829 * qJDD(4) + t720 * t826;
t598 = (qJD(6) - t759) * t736 + t897;
t732 = t734 ^ 2;
t733 = t736 ^ 2;
t758 = t759 ^ 2;
t760 = t766 ^ 2;
t761 = t768 ^ 2;
t948 = 2 * qJD(5);
t947 = pkin(2) * t527;
t946 = pkin(5) * t816;
t628 = t691 * t828 + t692 * t831;
t945 = pkin(7) * t628;
t852 = -pkin(3) * t833 + qJDD(2) * pkin(8) + t623;
t583 = -t827 * t747 + t830 * t852;
t787 = qJD(4) * pkin(4) - qJ(5) * t918;
t813 = t830 ^ 2;
t810 = t813 * t833;
t554 = -pkin(4) * t810 + qJ(5) * t776 - qJD(4) * t787 + t583;
t837 = t827 * t852;
t924 = t827 * t833;
t834 = -t837 - t775 * qJ(5) + qJDD(4) * pkin(4) + (pkin(4) * t924 + qJ(5) * t911 - t747) * t830;
t899 = t554 * t816 - t820 * t834;
t464 = t768 * t948 + t899;
t465 = -0.2e1 * qJD(5) * t766 + t820 * t554 + t816 * t834;
t401 = -t464 * t820 + t465 * t816;
t944 = t401 * t827;
t943 = t401 * t830;
t701 = pkin(5) * t766 - pkin(9) * t768;
t832 = qJD(4) ^ 2;
t445 = -qJDD(4) * pkin(5) - t832 * pkin(9) + (t948 + t701) * t768 + t899;
t942 = t445 * t826;
t941 = t445 * t829;
t616 = -qJDD(2) * pkin(3) - t833 * pkin(8) + t622;
t566 = -t776 * pkin(4) - qJ(5) * t810 + t787 * t918 + qJDD(5) + t616;
t938 = t566 * t816;
t937 = t566 * t820;
t936 = t616 * t827;
t935 = t616 * t830;
t630 = t666 + t890;
t934 = t630 * t826;
t933 = t630 * t829;
t709 = qJDD(4) + t718;
t932 = t709 * t816;
t931 = t709 * t820;
t930 = t759 * t826;
t929 = t759 * t829;
t798 = t830 * t924;
t788 = qJDD(4) + t798;
t928 = t788 * t827;
t789 = qJDD(4) - t798;
t927 = t789 * t827;
t926 = t789 * t830;
t812 = t827 ^ 2;
t925 = t812 * t833;
t923 = t828 * t752;
t921 = t831 * t752;
t446 = -pkin(5) * t832 + qJDD(4) * pkin(9) - t701 * t766 + t465;
t916 = qJD(4) * t766;
t895 = -t720 + t916;
t503 = pkin(5) * t667 + t895 * pkin(9) + t566;
t409 = t829 * t446 + t826 * t503;
t919 = t812 + t813;
t914 = qJD(4) * t816;
t913 = qJD(4) * t820;
t910 = qJDD(2) * t819;
t908 = qJDD(4) * t821;
t907 = t816 * t666;
t906 = t820 * t666;
t905 = t817 * t718;
t904 = t821 * t718;
t903 = -pkin(5) * t820 - pkin(4);
t901 = t822 * t910;
t408 = t446 * t826 - t829 * t503;
t402 = t464 * t816 + t820 * t465;
t582 = t830 * t747 + t837;
t492 = t582 * t827 + t830 * t583;
t729 = -t785 * t818 - t822 * t786;
t894 = t817 * t798;
t893 = t821 * t798;
t781 = qJDD(2) * t831 - t828 * t833;
t892 = -pkin(7) * t781 - t923;
t854 = qJDD(2) * t828 + t831 * t833;
t891 = -pkin(7) * t854 + t921;
t359 = t408 * t826 + t409 * t829;
t333 = t359 * t816 - t445 * t820;
t334 = t359 * t820 + t445 * t816;
t306 = -t333 * t827 + t334 * t830;
t358 = -t408 * t829 + t409 * t826;
t288 = t306 * t817 - t358 * t821;
t289 = t306 * t821 + t358 * t817;
t888 = t288 * t831 + t289 * t828;
t350 = t402 * t830 - t944;
t339 = t350 * t817 - t566 * t821;
t340 = t350 * t821 + t566 * t817;
t887 = t339 * t831 + t340 * t828;
t602 = -t685 - t647;
t510 = -t598 * t829 - t602 * t826;
t638 = t732 + t733;
t473 = t510 * t816 + t638 * t820;
t474 = t510 * t820 - t638 * t816;
t411 = -t473 * t827 + t474 * t830;
t508 = -t598 * t826 + t602 * t829;
t378 = t411 * t817 - t508 * t821;
t379 = t411 * t821 + t508 * t817;
t886 = t378 * t831 + t379 * t828;
t600 = (-qJD(6) - t759) * t736 - t897;
t511 = t600 * t829 - t601 * t826;
t664 = -t733 + t732;
t477 = t511 * t816 + t664 * t820;
t478 = t511 * t820 - t664 * t816;
t417 = -t477 * t827 + t478 * t830;
t509 = -t600 * t826 - t601 * t829;
t385 = t417 * t817 + t509 * t821;
t386 = t417 * t821 - t509 * t817;
t885 = t385 * t831 + t386 * t828;
t645 = -t758 - t732;
t553 = t645 * t829 - t961;
t481 = t553 * t816 + t600 * t820;
t482 = t553 * t820 - t600 * t816;
t422 = -t481 * t827 + t482 * t830;
t552 = t645 * t826 + t960;
t397 = t422 * t817 - t552 * t821;
t398 = t422 * t821 + t552 * t817;
t884 = t397 * t831 + t398 * t828;
t657 = -t733 - t758;
t558 = -t657 * t826 - t933;
t487 = t558 * t816 - t601 * t820;
t488 = t558 * t820 + t601 * t816;
t424 = -t487 * t827 + t488 * t830;
t557 = t657 * t829 - t934;
t399 = t424 * t817 - t557 * t821;
t400 = t424 * t821 + t557 * t817;
t883 = t399 * t831 + t400 * t828;
t678 = -t733 + t758;
t564 = -t678 * t826 + t960;
t496 = t564 * t816 + t602 * t820;
t498 = t564 * t820 - t602 * t816;
t433 = -t496 * t827 + t498 * t830;
t562 = -t678 * t829 - t961;
t404 = t433 * t817 + t562 * t821;
t406 = t433 * t821 - t562 * t817;
t882 = t404 * t831 + t406 * t828;
t677 = t732 - t758;
t565 = t677 * t829 - t934;
t497 = t565 * t816 + t598 * t820;
t499 = t565 * t820 - t598 * t816;
t434 = -t497 * t827 + t499 * t830;
t563 = -t677 * t826 - t933;
t405 = t434 * t817 + t563 * t821;
t407 = t434 * t821 - t563 * t817;
t881 = t405 * t831 + t407 * t828;
t646 = -qJD(6) * t736 - t897;
t589 = -t646 * t826 + t734 * t929;
t537 = t589 * t816 + t906;
t539 = t589 * t820 - t907;
t461 = -t537 * t827 + t539 * t830;
t588 = -t646 * t829 - t734 * t930;
t435 = t461 * t817 + t588 * t821;
t437 = t461 * t821 - t588 * t817;
t880 = t435 * t831 + t437 * t828;
t591 = t647 * t829 - t736 * t930;
t538 = t591 * t816 - t906;
t540 = t591 * t820 + t907;
t462 = -t538 * t827 + t540 * t830;
t590 = -t647 * t826 - t736 * t929;
t436 = t462 * t817 + t590 * t821;
t438 = t462 * t821 - t590 * t817;
t879 = t436 * t831 + t438 * t828;
t625 = (-t734 * t829 + t736 * t826) * t759;
t567 = t625 * t816 - t820 * t890;
t568 = t625 * t820 + t816 * t890;
t484 = -t567 * t827 + t568 * t830;
t624 = (t734 * t826 + t736 * t829) * t759;
t453 = t484 * t817 + t624 * t821;
t454 = t484 * t821 - t624 * t817;
t878 = t453 * t831 + t454 * t828;
t455 = t492 * t817 - t616 * t821;
t456 = t492 * t821 + t616 * t817;
t877 = t455 * t831 + t456 * t828;
t669 = -t896 + t915;
t671 = -t720 - t916;
t605 = t669 * t816 + t671 * t820;
t607 = t669 * t820 - t671 * t816;
t517 = -t605 * t827 + t607 * t830;
t665 = -t760 - t761;
t479 = t517 * t817 - t665 * t821;
t480 = t517 * t821 + t665 * t817;
t876 = t479 * t831 + t480 * t828;
t604 = -t667 * t816 - t820 * t895;
t606 = -t667 * t820 + t816 * t895;
t516 = -t604 * t827 + t606 * t830;
t716 = -t761 + t760;
t493 = t516 * t817 + t716 * t821;
t494 = t516 * t821 - t716 * t817;
t875 = t493 * t831 + t494 * t828;
t702 = -t832 - t760;
t636 = t702 * t816 + t962;
t637 = t702 * t820 - t963;
t544 = -t636 * t827 + t637 * t830;
t512 = t544 * t817 - t667 * t821;
t513 = t544 * t821 + t667 * t817;
t874 = t512 * t831 + t513 * t828;
t753 = t760 - t832;
t651 = t753 * t816 + t931;
t654 = t753 * t820 - t932;
t572 = -t651 * t827 + t654 * t830;
t529 = t572 * t817 - t669 * t821;
t532 = t572 * t821 + t669 * t817;
t871 = t529 * t831 + t532 * t828;
t754 = -t761 + t832;
t652 = t754 * t820 + t963;
t655 = -t754 * t816 + t962;
t573 = -t652 * t827 + t655 * t830;
t530 = t573 * t817 + t671 * t821;
t533 = t573 * t821 - t671 * t817;
t870 = t530 * t831 + t533 * t828;
t755 = -t761 - t832;
t653 = t755 * t820 - t932;
t656 = -t755 * t816 - t931;
t574 = -t653 * t827 + t656 * t830;
t531 = t574 * t817 + t821 * t895;
t534 = t574 * t821 - t817 * t895;
t869 = t531 * t831 + t534 * t828;
t659 = t766 * t914 - t820 * t896;
t660 = t766 * t913 + t816 * t896;
t586 = -t659 * t827 + t660 * t830;
t548 = t586 * t817 + t904;
t550 = t586 * t821 - t905;
t868 = t548 * t831 + t550 * t828;
t661 = t720 * t816 + t768 * t913;
t662 = t720 * t820 - t768 * t914;
t587 = -t661 * t827 + t662 * t830;
t549 = t587 * t817 - t904;
t551 = t587 * t821 + t905;
t867 = t549 * t831 + t551 * t828;
t491 = t582 * t830 - t583 * t827;
t683 = (-t766 * t816 - t768 * t820) * qJD(4);
t684 = (-t766 * t820 + t768 * t816) * qJD(4);
t621 = -t683 * t827 + t684 * t830;
t613 = t621 * t817 - t908;
t802 = t817 * qJDD(4);
t614 = t621 * t821 + t802;
t866 = t613 * t831 + t614 * t828;
t774 = 0.2e1 * t900 + t909;
t777 = t807 - 0.2e1 * t902;
t722 = -t774 * t827 + t777 * t830;
t784 = t810 - t925;
t672 = t722 * t817 + t784 * t821;
t673 = t722 * t821 - t784 * t817;
t865 = t672 * t831 + t673 * t828;
t797 = -t810 - t832;
t744 = t797 * t830 - t928;
t679 = t744 * t817 + t777 * t821;
t681 = t744 * t821 - t777 * t817;
t864 = t679 * t831 + t681 * t828;
t795 = -t832 - t925;
t746 = -t795 * t827 - t926;
t680 = t746 * t817 - t774 * t821;
t682 = t746 * t821 + t774 * t817;
t863 = t680 * t831 + t682 * t828;
t627 = t691 * t831 - t692 * t828;
t796 = t810 - t832;
t743 = t796 * t830 - t927;
t693 = t743 * t817 - t807 * t821;
t695 = t743 * t821 + t807 * t817;
t862 = t693 * t831 + t695 * t828;
t773 = t830 * t788;
t794 = t832 - t925;
t745 = -t794 * t827 + t773;
t694 = t745 * t817 - t821 * t909;
t696 = t745 * t821 + t817 * t909;
t861 = t694 * t831 + t696 * t828;
t750 = -t776 * t827 - t813 * t911;
t697 = t750 * t817 - t893;
t699 = t750 * t821 + t894;
t860 = t697 * t831 + t699 * t828;
t751 = t775 * t830 - t812 * t911;
t698 = t751 * t817 + t893;
t700 = t751 * t821 - t894;
t859 = t698 * t831 + t700 * t828;
t780 = t919 * qJDD(2);
t783 = t810 + t925;
t726 = t780 * t817 + t783 * t821;
t727 = t780 * t821 - t783 * t817;
t858 = t726 * t831 + t727 * t828;
t772 = t919 * t911;
t748 = t772 * t817 - t908;
t749 = t772 * t821 + t802;
t857 = t748 * t831 + t749 * t828;
t764 = t854 * t823;
t856 = t764 * t822 + t781 * t818;
t714 = t764 * t818 - t781 * t822;
t728 = t785 * t822 - t786 * t818;
t283 = qJ(5) * t334 + (-pkin(9) * t816 + t903) * t358;
t296 = -qJ(5) * t333 + (-pkin(9) * t820 + t946) * t358;
t305 = t333 * t830 + t334 * t827;
t265 = -pkin(8) * t305 - t283 * t827 + t296 * t830;
t276 = -pkin(3) * t305 - pkin(4) * t333 + pkin(5) * t445 - pkin(9) * t359;
t252 = -pkin(2) * t305 + qJ(3) * t289 + t265 * t817 + t276 * t821;
t255 = -qJ(3) * t288 + t265 * t821 - t276 * t817;
t268 = -t288 * t828 + t289 * t831;
t851 = pkin(7) * t268 + t252 * t831 + t255 * t828;
t341 = -pkin(9) * t508 - t358;
t316 = qJ(5) * t474 + t341 * t816 + t508 * t903;
t320 = -qJ(5) * t473 + t341 * t820 + t508 * t946;
t410 = t473 * t830 + t474 * t827;
t285 = -pkin(8) * t410 - t316 * t827 + t320 * t830;
t312 = -pkin(3) * t410 - pkin(4) * t473 - pkin(5) * t638 - pkin(9) * t510 - t359;
t269 = -pkin(2) * t410 + qJ(3) * t379 + t285 * t817 + t312 * t821;
t271 = -qJ(3) * t378 + t285 * t821 - t312 * t817;
t332 = -t378 * t828 + t379 * t831;
t850 = pkin(7) * t332 + t269 * t831 + t271 * t828;
t349 = t402 * t827 + t943;
t382 = -pkin(4) * t566 + qJ(5) * t402;
t310 = -pkin(8) * t349 - qJ(5) * t943 - t382 * t827;
t319 = -pkin(3) * t349 - pkin(4) * t401;
t272 = -pkin(2) * t349 + qJ(3) * t340 + t310 * t817 + t319 * t821;
t279 = -qJ(3) * t339 + t310 * t821 - t319 * t817;
t309 = -t339 * t828 + t340 * t831;
t849 = pkin(7) * t309 + t272 * t831 + t279 * t828;
t387 = -pkin(5) * t552 + t408;
t418 = -pkin(9) * t552 + t942;
t328 = -pkin(4) * t552 + qJ(5) * t482 + t387 * t820 + t418 * t816;
t335 = -qJ(5) * t481 - t387 * t816 + t418 * t820;
t421 = t481 * t830 + t482 * t827;
t300 = -pkin(8) * t421 - t328 * t827 + t335 * t830;
t342 = -pkin(3) * t421 - pkin(4) * t481 - pkin(5) * t600 - pkin(9) * t553 + t941;
t275 = -pkin(2) * t421 + qJ(3) * t398 + t300 * t817 + t342 * t821;
t280 = -qJ(3) * t397 + t300 * t821 - t342 * t817;
t346 = -t397 * t828 + t398 * t831;
t848 = pkin(7) * t346 + t275 * t831 + t280 * t828;
t389 = -pkin(5) * t557 + t409;
t420 = -pkin(9) * t557 + t941;
t329 = -pkin(4) * t557 + qJ(5) * t488 + t389 * t820 + t420 * t816;
t337 = -qJ(5) * t487 - t389 * t816 + t420 * t820;
t423 = t487 * t830 + t488 * t827;
t301 = -pkin(8) * t423 - t329 * t827 + t337 * t830;
t343 = -pkin(3) * t423 - pkin(4) * t487 + pkin(5) * t601 - pkin(9) * t558 - t942;
t278 = -pkin(2) * t423 + qJ(3) * t400 + t301 * t817 + t343 * t821;
t281 = -qJ(3) * t399 + t301 * t821 - t343 * t817;
t348 = -t399 * t828 + t400 * t831;
t847 = pkin(7) * t348 + t278 * t831 + t281 * t828;
t380 = -pkin(4) * t665 + qJ(5) * t607 + t402;
t388 = -qJ(5) * t605 - t401;
t515 = t605 * t830 + t607 * t827;
t325 = -pkin(8) * t515 - t380 * t827 + t388 * t830;
t468 = -pkin(3) * t515 - pkin(4) * t605;
t311 = -pkin(2) * t515 + qJ(3) * t480 + t325 * t817 + t468 * t821;
t313 = -qJ(3) * t479 + t325 * t821 - t468 * t817;
t419 = -t479 * t828 + t480 * t831;
t846 = pkin(7) * t419 + t311 * t831 + t313 * t828;
t489 = -pkin(4) * t667 + qJ(5) * t637 - t937;
t519 = -qJ(5) * t636 + t938;
t543 = t636 * t830 + t637 * t827;
t403 = -pkin(8) * t543 - t489 * t827 + t519 * t830;
t415 = -pkin(3) * t543 - pkin(4) * t636 + t464;
t336 = -pkin(2) * t543 + qJ(3) * t513 + t403 * t817 + t415 * t821;
t344 = -qJ(3) * t512 + t403 * t821 - t415 * t817;
t441 = -t512 * t828 + t513 * t831;
t845 = pkin(7) * t441 + t336 * t831 + t344 * t828;
t495 = pkin(4) * t895 + qJ(5) * t656 + t938;
t524 = -qJ(5) * t653 + t937;
t571 = t653 * t830 + t656 * t827;
t414 = -pkin(8) * t571 - t495 * t827 + t524 * t830;
t425 = -pkin(3) * t571 - pkin(4) * t653 + t465;
t345 = -pkin(2) * t571 + qJ(3) * t534 + t414 * t817 + t425 * t821;
t353 = -qJ(3) * t531 + t414 * t821 - t425 * t817;
t452 = -t531 * t828 + t534 * t831;
t844 = pkin(7) * t452 + t345 * t831 + t353 * t828;
t361 = qJ(3) * t456 - (-pkin(3) * t821 - pkin(8) * t817 - pkin(2)) * t491;
t377 = -qJ(3) * t455 - (pkin(3) * t817 - pkin(8) * t821) * t491;
t395 = -t455 * t828 + t456 * t831;
t843 = pkin(7) * t395 + t361 * t831 + t377 * t828;
t740 = t797 * t827 + t773;
t555 = -pkin(3) * t740 + t582;
t580 = -pkin(8) * t740 + t936;
t457 = -pkin(2) * t740 + qJ(3) * t681 + t555 * t821 + t580 * t817;
t466 = -qJ(3) * t679 - t555 * t817 + t580 * t821;
t618 = -t679 * t828 + t681 * t831;
t842 = pkin(7) * t618 + t457 * t831 + t466 * t828;
t742 = t795 * t830 - t927;
t556 = -pkin(3) * t742 + t583;
t581 = -pkin(8) * t742 + t935;
t458 = -pkin(2) * t742 + qJ(3) * t682 + t556 * t821 + t581 * t817;
t467 = -qJ(3) * t680 - t556 * t817 + t581 * t821;
t619 = -t680 * t828 + t682 * t831;
t841 = pkin(7) * t619 + t458 * t831 + t467 * t828;
t475 = qJ(3) * t727 + t491 * t817;
t476 = -qJ(3) * t726 + t491 * t821;
t650 = -t726 * t828 + t727 * t831;
t840 = pkin(7) * t650 + t475 * t831 + t476 * t828;
t839 = -pkin(7) * t855 - t958;
t838 = pkin(7) * t724 + t959;
t518 = pkin(2) * t747 + qJ(3) * t898;
t836 = pkin(7) * t964 + qJ(3) * t940 + t518 * t831;
t803 = t823 * qJDD(2);
t790 = t818 * t910;
t765 = t781 * t823;
t763 = t781 * t819;
t762 = t854 * t819;
t741 = t794 * t830 + t928;
t739 = t796 * t827 + t926;
t738 = (t775 + t900) * t827;
t737 = (t776 - t902) * t830;
t721 = t774 * t830 + t777 * t827;
t715 = -t765 * t818 - t822 * t854;
t713 = t765 * t822 - t818 * t854;
t707 = t724 * t823;
t704 = t724 * t819;
t663 = -t748 * t828 + t749 * t831;
t658 = t857 * t823;
t649 = -t921 + (t762 * t819 + t764 * t823) * pkin(7);
t648 = -t923 + (-t763 * t819 - t765 * t823) * pkin(7);
t644 = t858 * t823;
t643 = t858 * t819;
t642 = t707 * t818 - t822 * t855;
t640 = -t707 * t822 - t818 * t855;
t635 = -t698 * t828 + t700 * t831;
t634 = -t697 * t828 + t699 * t831;
t633 = -t694 * t828 + t696 * t831;
t632 = -t693 * t828 + t695 * t831;
t626 = t628 * t823;
t620 = t683 * t830 + t684 * t827;
t612 = -pkin(1) * t763 + t691 * t819 + t823 * t891;
t611 = pkin(1) * t762 + t692 * t819 + t823 * t892;
t610 = -t672 * t828 + t673 * t831;
t609 = -pkin(2) * t778 - t623;
t608 = pkin(2) * t779 - t622;
t597 = -t738 * t819 + t823 * t859;
t596 = -t737 * t819 + t823 * t860;
t595 = -t741 * t819 + t823 * t861;
t594 = -t739 * t819 + t823 * t862;
t593 = -t627 * t823 + t752 * t819;
t592 = -t627 * t819 - t752 * t823;
t585 = t661 * t830 + t662 * t827;
t584 = t659 * t830 + t660 * t827;
t579 = -t742 * t819 + t823 * t863;
t578 = -t740 * t819 + t823 * t864;
t577 = t742 * t823 + t819 * t863;
t576 = t740 * t823 + t819 * t864;
t570 = t652 * t830 + t655 * t827;
t569 = t651 * t830 + t654 * t827;
t561 = -t644 * t818 + t650 * t822;
t560 = t644 * t822 + t650 * t818;
t559 = -t721 * t819 + t823 * t865;
t542 = (t703 * t819 + t823 * t953) * pkin(7) + t958;
t541 = (t704 * t819 + t707 * t823) * pkin(7) + t959;
t536 = pkin(2) * t679 + pkin(3) * t777 + pkin(8) * t744 - t935;
t535 = pkin(2) * t680 - pkin(3) * t774 + pkin(8) * t746 + t936;
t523 = -pkin(1) * t592 + t823 * t945;
t522 = -t593 * t818 + t628 * t822;
t521 = t593 * t822 + t628 * t818;
t520 = -t613 * t828 + t614 * t831;
t514 = t604 * t830 + t606 * t827;
t507 = -t579 * t818 + t619 * t822;
t506 = -t578 * t818 + t618 * t822;
t505 = t579 * t822 + t619 * t818;
t504 = t578 * t822 + t618 * t818;
t500 = (-t592 * t819 - t593 * t823) * pkin(7);
t486 = pkin(1) * t703 - t609 * t819 + t823 * t838;
t485 = pkin(1) * t704 - t608 * t819 + t823 * t839;
t483 = t567 * t830 + t568 * t827;
t472 = -t620 * t819 + t823 * t866;
t471 = pkin(2) * t726 + pkin(3) * t783 + pkin(8) * t780 + t492;
t470 = -t549 * t828 + t551 * t831;
t469 = -t548 * t828 + t550 * t831;
t460 = t538 * t830 + t540 * t827;
t459 = t537 * t830 + t539 * t827;
t451 = -t530 * t828 + t533 * t831;
t450 = -t529 * t828 + t532 * t831;
t447 = t964 * t823;
t443 = -t448 * t823 + t747 * t819;
t442 = -t448 * t819 - t747 * t823;
t440 = -t585 * t819 + t823 * t867;
t439 = -t584 * t819 + t823 * t868;
t432 = t497 * t830 + t499 * t827;
t431 = t496 * t830 + t498 * t827;
t430 = -t571 * t819 + t823 * t869;
t429 = -t570 * t819 + t823 * t870;
t428 = -t569 * t819 + t823 * t871;
t427 = t571 * t823 + t819 * t869;
t426 = -t493 * t828 + t494 * t831;
t416 = t477 * t830 + t478 * t827;
t413 = -t543 * t819 + t823 * t874;
t412 = t543 * t823 + t819 * t874;
t396 = pkin(2) * t455 - pkin(3) * t616 + pkin(8) * t492;
t394 = -t453 * t828 + t454 * t831;
t393 = -t475 * t828 + t476 * t831 + (-t643 * t819 - t644 * t823) * pkin(7);
t392 = -t514 * t819 + t823 * t875;
t391 = -t515 * t819 + t823 * t876;
t390 = t515 * t823 + t819 * t876;
t384 = -t443 * t818 + t822 * t964;
t383 = t443 * t822 + t818 * t964;
t381 = pkin(2) * t531 + pkin(3) * t895 + pkin(8) * t574 + t495 * t830 + t524 * t827;
t376 = -t430 * t818 + t452 * t822;
t375 = t430 * t822 + t452 * t818;
t374 = pkin(2) * t512 - pkin(3) * t667 + pkin(8) * t544 + t489 * t830 + t519 * t827;
t373 = -t436 * t828 + t438 * t831;
t372 = -t435 * t828 + t437 * t831;
t371 = -t458 * t828 + t467 * t831 + (-t577 * t819 - t579 * t823) * pkin(7);
t370 = -t457 * t828 + t466 * t831 + (-t576 * t819 - t578 * t823) * pkin(7);
t369 = t491 * t819 + t823 * t877;
t368 = -t491 * t823 + t819 * t877;
t367 = -t483 * t819 + t823 * t878;
t366 = -pkin(1) * t643 - t471 * t819 + t823 * t840;
t365 = -t413 * t818 + t441 * t822;
t364 = t413 * t822 + t441 * t818;
t363 = -pkin(1) * t577 - t535 * t819 + t823 * t841;
t362 = -pkin(1) * t576 - t536 * t819 + t823 * t842;
t360 = qJ(3) * t939 - t518 * t828 + (-t442 * t819 - t443 * t823) * pkin(7);
t357 = -t405 * t828 + t407 * t831;
t356 = -t404 * t828 + t406 * t831;
t355 = -t460 * t819 + t823 * t879;
t354 = -t459 * t819 + t823 * t880;
t352 = -t391 * t818 + t419 * t822;
t351 = t391 * t822 + t419 * t818;
t347 = -pkin(1) * t442 + t819 * t947 + t823 * t836;
t338 = -t385 * t828 + t386 * t831;
t331 = -t369 * t818 + t395 * t822;
t330 = t369 * t822 + t395 * t818;
t327 = -t432 * t819 + t823 * t881;
t326 = -t431 * t819 + t823 * t882;
t324 = -t423 * t819 + t823 * t883;
t323 = t423 * t823 + t819 * t883;
t322 = -t421 * t819 + t823 * t884;
t321 = t421 * t823 + t819 * t884;
t318 = -t819 * t416 + t823 * t885;
t317 = pkin(2) * t479 - pkin(3) * t665 + pkin(8) * t517 + t380 * t830 + t388 * t827;
t315 = -t819 * t410 + t823 * t886;
t314 = t823 * t410 + t819 * t886;
t308 = -t324 * t818 + t348 * t822;
t307 = t324 * t822 + t348 * t818;
t304 = -t322 * t818 + t346 * t822;
t303 = t322 * t822 + t346 * t818;
t302 = -t345 * t828 + t353 * t831 + (-t427 * t819 - t430 * t823) * pkin(7);
t299 = -t315 * t818 + t332 * t822;
t298 = t315 * t822 + t332 * t818;
t297 = -t361 * t828 + t377 * t831 + (-t368 * t819 - t369 * t823) * pkin(7);
t295 = -t336 * t828 + t344 * t831 + (-t412 * t819 - t413 * t823) * pkin(7);
t294 = -t819 * t349 + t823 * t887;
t293 = t823 * t349 + t819 * t887;
t292 = -pkin(1) * t427 - t381 * t819 + t823 * t844;
t291 = -pkin(1) * t368 - t396 * t819 + t823 * t843;
t290 = pkin(2) * t339 - pkin(3) * t566 + pkin(8) * t350 - qJ(5) * t944 + t382 * t830;
t287 = pkin(2) * t399 - pkin(3) * t557 + pkin(8) * t424 + t329 * t830 + t337 * t827;
t286 = pkin(2) * t397 - pkin(3) * t552 + pkin(8) * t422 + t328 * t830 + t335 * t827;
t284 = -pkin(1) * t412 - t374 * t819 + t823 * t845;
t282 = pkin(2) * t378 - pkin(3) * t508 + pkin(8) * t411 + t316 * t830 + t320 * t827;
t277 = -t311 * t828 + t313 * t831 + (-t390 * t819 - t391 * t823) * pkin(7);
t274 = -t294 * t818 + t309 * t822;
t273 = t294 * t822 + t309 * t818;
t270 = -pkin(1) * t390 - t317 * t819 + t823 * t846;
t267 = -t819 * t305 + t823 * t888;
t266 = t823 * t305 + t819 * t888;
t264 = -t278 * t828 + t281 * t831 + (-t323 * t819 - t324 * t823) * pkin(7);
t263 = -t275 * t828 + t280 * t831 + (-t321 * t819 - t322 * t823) * pkin(7);
t262 = pkin(2) * t288 - pkin(3) * t358 + pkin(8) * t306 + t283 * t830 + t296 * t827;
t261 = -pkin(1) * t323 - t287 * t819 + t823 * t847;
t260 = -pkin(1) * t321 - t286 * t819 + t823 * t848;
t259 = -t269 * t828 + t271 * t831 + (-t314 * t819 - t315 * t823) * pkin(7);
t258 = -t267 * t818 + t268 * t822;
t257 = t267 * t822 + t268 * t818;
t256 = -t272 * t828 + t279 * t831 + (-t293 * t819 - t294 * t823) * pkin(7);
t254 = -pkin(1) * t293 - t290 * t819 + t823 * t849;
t253 = -pkin(1) * t314 - t282 * t819 + t823 * t850;
t251 = -t252 * t828 + t255 * t831 + (-t266 * t819 - t267 * t823) * pkin(7);
t250 = -pkin(1) * t266 - t262 * t819 + t823 * t851;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t729, 0, 0, 0, 0, 0, 0, t715, t714, 0, t522, 0, 0, 0, 0, 0, 0, t642, t641, 0, t384, 0, 0, 0, 0, 0, 0, t506, t507, t561, t331, 0, 0, 0, 0, 0, 0, t365, t376, t352, t274, 0, 0, 0, 0, 0, 0, t304, t308, t299, t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t728, 0, 0, 0, 0, 0, 0, t713, -t856, 0, t521, 0, 0, 0, 0, 0, 0, t640, t639, 0, t383, 0, 0, 0, 0, 0, 0, t504, t505, t560, t330, 0, 0, 0, 0, 0, 0, t364, t375, t351, t273, 0, 0, 0, 0, 0, 0, t303, t307, t298, t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t920, 0, 0, 0, 0, 0, 0, t763, -t762, 0, t592, 0, 0, 0, 0, 0, 0, -t704, -t703, 0, t442, 0, 0, 0, 0, 0, 0, t576, t577, t643, t368, 0, 0, 0, 0, 0, 0, t412, t427, t390, t293, 0, 0, 0, 0, 0, 0, t321, t323, t314, t266; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t957, -t956, -t728, -qJ(1) * t728, 0, 0, -t714, 0, t715, t790, -qJ(1) * t713 - t612 * t818 + t648 * t822, qJ(1) * t856 - t611 * t818 + t649 * t822, -t626 * t818 + t627 * t822, -qJ(1) * t521 + t500 * t822 - t523 * t818, 0, 0, -t641, 0, t642, t790, -qJ(1) * t640 - t485 * t818 + t541 * t822, -qJ(1) * t639 - t486 * t818 + t542 * t822, -t447 * t818 + t448 * t822, -qJ(1) * t383 - t347 * t818 + t360 * t822, -t597 * t818 + t635 * t822, -t559 * t818 + t610 * t822, -t595 * t818 + t633 * t822, -t596 * t818 + t634 * t822, -t594 * t818 + t632 * t822, -t658 * t818 + t663 * t822, -qJ(1) * t504 - t362 * t818 + t370 * t822, -qJ(1) * t505 - t363 * t818 + t371 * t822, -qJ(1) * t560 - t366 * t818 + t393 * t822, -qJ(1) * t330 - t291 * t818 + t297 * t822, -t440 * t818 + t470 * t822, -t392 * t818 + t426 * t822, -t429 * t818 + t451 * t822, -t439 * t818 + t469 * t822, -t428 * t818 + t450 * t822, -t472 * t818 + t520 * t822, -qJ(1) * t364 - t284 * t818 + t295 * t822, -qJ(1) * t375 - t292 * t818 + t302 * t822, -qJ(1) * t351 - t270 * t818 + t277 * t822, -qJ(1) * t273 - t254 * t818 + t256 * t822, -t355 * t818 + t373 * t822, -t318 * t818 + t338 * t822, -t326 * t818 + t356 * t822, -t354 * t818 + t372 * t822, -t327 * t818 + t357 * t822, -t367 * t818 + t394 * t822, -qJ(1) * t303 - t260 * t818 + t263 * t822, -qJ(1) * t307 - t261 * t818 + t264 * t822, -qJ(1) * t298 - t253 * t818 + t259 * t822, -qJ(1) * t257 - t250 * t818 + t251 * t822; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t956, -t957, t729, qJ(1) * t729, 0, 0, t856, 0, t713, -t901, qJ(1) * t715 + t612 * t822 + t648 * t818, qJ(1) * t714 + t611 * t822 + t649 * t818, t626 * t822 + t627 * t818, qJ(1) * t522 + t500 * t818 + t523 * t822, 0, 0, -t639, 0, t640, -t901, qJ(1) * t642 + t485 * t822 + t541 * t818, qJ(1) * t641 + t486 * t822 + t542 * t818, t447 * t822 + t448 * t818, qJ(1) * t384 + t347 * t822 + t360 * t818, t597 * t822 + t635 * t818, t559 * t822 + t610 * t818, t595 * t822 + t633 * t818, t596 * t822 + t634 * t818, t594 * t822 + t632 * t818, t658 * t822 + t663 * t818, qJ(1) * t506 + t362 * t822 + t370 * t818, qJ(1) * t507 + t363 * t822 + t371 * t818, qJ(1) * t561 + t366 * t822 + t393 * t818, qJ(1) * t331 + t291 * t822 + t297 * t818, t440 * t822 + t470 * t818, t392 * t822 + t426 * t818, t429 * t822 + t451 * t818, t439 * t822 + t469 * t818, t428 * t822 + t450 * t818, t472 * t822 + t520 * t818, qJ(1) * t365 + t284 * t822 + t295 * t818, qJ(1) * t376 + t292 * t822 + t302 * t818, qJ(1) * t352 + t270 * t822 + t277 * t818, qJ(1) * t274 + t254 * t822 + t256 * t818, t355 * t822 + t373 * t818, t318 * t822 + t338 * t818, t326 * t822 + t356 * t818, t354 * t822 + t372 * t818, t327 * t822 + t357 * t818, t367 * t822 + t394 * t818, qJ(1) * t304 + t260 * t822 + t263 * t818, qJ(1) * t308 + t261 * t822 + t264 * t818, qJ(1) * t299 + t253 * t822 + t259 * t818, qJ(1) * t258 + t250 * t822 + t251 * t818; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t785, t786, 0, 0, 0, 0, t762, 0, t763, t803, pkin(1) * t765 - t691 * t823 + t819 * t891, -pkin(1) * t764 - t692 * t823 + t819 * t892, t628 * t819, pkin(1) * t593 + t819 * t945, 0, 0, t703, 0, -t704, t803, -pkin(1) * t707 + t608 * t823 + t819 * t839, -pkin(1) * t953 + t609 * t823 + t819 * t838, t964 * t819, pkin(1) * t443 + t819 * t836 - t823 * t947, t738 * t823 + t819 * t859, t721 * t823 + t819 * t865, t741 * t823 + t819 * t861, t737 * t823 + t819 * t860, t739 * t823 + t819 * t862, t857 * t819, pkin(1) * t578 + t536 * t823 + t819 * t842, pkin(1) * t579 + t535 * t823 + t819 * t841, pkin(1) * t644 + t471 * t823 + t819 * t840, pkin(1) * t369 + t396 * t823 + t819 * t843, t585 * t823 + t819 * t867, t514 * t823 + t819 * t875, t570 * t823 + t819 * t870, t584 * t823 + t819 * t868, t569 * t823 + t819 * t871, t620 * t823 + t819 * t866, pkin(1) * t413 + t374 * t823 + t819 * t845, pkin(1) * t430 + t381 * t823 + t819 * t844, pkin(1) * t391 + t317 * t823 + t819 * t846, pkin(1) * t294 + t290 * t823 + t819 * t849, t460 * t823 + t819 * t879, t823 * t416 + t819 * t885, t431 * t823 + t819 * t882, t459 * t823 + t819 * t880, t432 * t823 + t819 * t881, t483 * t823 + t819 * t878, pkin(1) * t322 + t286 * t823 + t819 * t848, pkin(1) * t324 + t287 * t823 + t819 * t847, pkin(1) * t315 + t282 * t823 + t819 * t850, pkin(1) * t267 + t262 * t823 + t819 * t851;];
tauB_reg  = t1;