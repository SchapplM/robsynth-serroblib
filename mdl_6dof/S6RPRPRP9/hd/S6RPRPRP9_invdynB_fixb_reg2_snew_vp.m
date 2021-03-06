% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RPRPRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RPRPRP9_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP9_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:07:16
% EndTime: 2019-05-05 18:07:43
% DurationCPUTime: 23.79s
% Computational Cost: add. (50166->614), mult. (105978->835), div. (0->0), fcn. (69954->8), ass. (0->435)
t782 = sin(pkin(9));
t783 = cos(pkin(9));
t788 = cos(qJ(3));
t843 = qJD(1) * t788;
t747 = qJD(3) * t783 - t782 * t843;
t748 = qJD(3) * t782 + t783 * t843;
t784 = sin(qJ(5));
t787 = cos(qJ(5));
t704 = t747 * t784 + t748 * t787;
t701 = t704 ^ 2;
t785 = sin(qJ(3));
t841 = t785 * qJD(1);
t772 = qJD(5) + t841;
t888 = t772 ^ 2;
t642 = t888 + t701;
t702 = -t787 * t747 + t748 * t784;
t648 = t704 * t702;
t839 = qJD(1) * qJD(3);
t821 = t788 * t839;
t837 = qJDD(1) * t785;
t753 = t821 + t837;
t749 = qJDD(5) + t753;
t900 = t648 + t749;
t870 = t900 * t784;
t561 = t642 * t787 + t870;
t869 = t900 * t787;
t570 = t642 * t784 - t869;
t511 = t561 * t782 + t570 * t783;
t822 = t785 * t839;
t835 = qJDD(1) * t788;
t754 = -t822 + t835;
t727 = t783 * qJDD(3) - t754 * t782;
t728 = qJDD(3) * t782 + t754 * t783;
t794 = qJD(5) * t702 - t727 * t784 - t728 * t787;
t862 = t702 * t772;
t909 = -t794 - t862;
t487 = t511 * t785 - t788 * t909;
t504 = t561 * t783 - t570 * t782;
t786 = sin(qJ(1));
t789 = cos(qJ(1));
t443 = t487 * t789 + t504 * t786;
t1005 = pkin(6) * t443;
t450 = t487 * t786 - t504 * t789;
t1004 = pkin(6) * t450;
t818 = -t787 * t727 + t784 * t728;
t592 = (qJD(5) + t772) * t704 + t818;
t535 = -t592 * t784 + t787 * t909;
t874 = t909 * t784;
t539 = -t592 * t787 - t874;
t475 = -t535 * t782 + t539 * t783;
t889 = t702 ^ 2;
t645 = t701 - t889;
t463 = t475 * t785 - t645 * t788;
t471 = t535 * t783 + t539 * t782;
t1003 = t463 * t786 + t471 * t789;
t1002 = t463 * t789 - t471 * t786;
t886 = pkin(7) + pkin(1);
t1001 = qJ(2) * t504 + t487 * t886;
t489 = t511 * t788 + t785 * t909;
t1000 = pkin(2) * t504 + t489 * t886;
t678 = t889 - t888;
t584 = t678 * t784 + t869;
t588 = t678 * t787 - t870;
t531 = t584 * t782 - t588 * t783;
t593 = (qJD(5) - t772) * t704 + t818;
t493 = t531 * t785 - t593 * t788;
t526 = t584 * t783 + t588 * t782;
t999 = t493 * t786 - t526 * t789;
t998 = t493 * t789 + t526 * t786;
t997 = -pkin(2) * t487 + pkin(3) * t909 + qJ(2) * t489 - qJ(4) * t511;
t897 = -t862 + t794;
t928 = -t593 * t787 - t897 * t784;
t929 = -t593 * t784 + t787 * t897;
t944 = t782 * t928 + t783 * t929;
t612 = -t889 - t701;
t945 = -t782 * t929 + t783 * t928;
t962 = -t612 * t788 + t785 * t945;
t979 = t786 * t944 - t789 * t962;
t995 = pkin(6) * t979;
t981 = t786 * t962 + t789 * t944;
t994 = pkin(6) * t981;
t992 = qJ(4) * t504;
t985 = -pkin(3) * t504 - pkin(4) * t561;
t983 = t475 * t788 + t645 * t785;
t982 = t531 * t788 + t593 * t785;
t680 = -t701 + t888;
t901 = -t648 + t749;
t868 = t901 * t784;
t930 = t787 * t680 + t868;
t629 = t787 * t901;
t931 = -t680 * t784 + t629;
t943 = t782 * t931 + t783 * t930;
t942 = -t782 * t930 + t783 * t931;
t963 = -t785 * t942 - t788 * t897;
t980 = -t786 * t963 + t789 * t943;
t978 = t786 * t943 + t789 * t963;
t960 = t612 * t785 + t788 * t945;
t977 = pkin(2) * t944 - t886 * t960;
t976 = qJ(2) * t944 - t886 * t962;
t975 = pkin(2) * t962 - pkin(3) * t612 - qJ(2) * t960 + qJ(4) * t945;
t972 = pkin(8) * t561;
t971 = pkin(8) * t570;
t969 = qJ(4) * t944;
t447 = -pkin(3) * t944 - pkin(4) * t929;
t961 = -t785 * t897 + t788 * t942;
t896 = -t888 - t889;
t905 = t787 * t896 - t868;
t908 = t784 * t896 + t629;
t922 = t782 * t905 + t783 * t908;
t959 = pkin(2) * t922;
t957 = pkin(8) * t929;
t956 = qJ(2) * t922;
t955 = qJ(4) * t922;
t923 = -t782 * t908 + t783 * t905;
t954 = qJ(4) * t923;
t951 = t785 * t923;
t950 = t786 * t922;
t949 = t788 * t923;
t948 = t789 * t922;
t947 = -pkin(3) * t922 - pkin(4) * t908;
t946 = -pkin(4) * t612 + pkin(8) * t928;
t938 = pkin(8) * t905;
t937 = pkin(8) * t908;
t936 = qJ(6) * t909;
t796 = (-t702 * t784 - t704 * t787) * t772;
t851 = t772 * t784;
t671 = t704 * t851;
t850 = t772 * t787;
t831 = t702 * t850;
t809 = t671 - t831;
t893 = t782 * t809 + t783 * t796;
t892 = -t782 * t796 + t783 * t809;
t906 = t749 * t788 - t785 * t892;
t927 = -t786 * t906 + t789 * t893;
t627 = -qJD(5) * t704 - t818;
t800 = -t627 * t784 + t831;
t810 = t787 * t627 + t702 * t851;
t890 = t782 * t800 + t783 * t810;
t832 = t788 * t648;
t891 = -t782 * t810 + t783 * t800;
t907 = -t785 * t891 - t832;
t926 = -t786 * t907 + t789 * t890;
t925 = t786 * t893 + t789 * t906;
t924 = t786 * t890 + t789 * t907;
t921 = 2 * qJD(6);
t712 = t747 * t748;
t899 = t712 + t753;
t919 = t782 * t899;
t918 = t783 * t899;
t833 = t785 * t648;
t904 = t788 * t891 - t833;
t903 = t749 * t785 + t788 * t892;
t791 = qJD(1) ^ 2;
t902 = t791 * t886;
t838 = qJD(2) * qJD(1);
t777 = 0.2e1 * t838;
t765 = t789 * g(1) + t786 * g(2);
t779 = qJDD(1) * qJ(2);
t803 = t765 - t779;
t798 = t777 - t803;
t805 = -t754 + t822;
t806 = t753 + t821;
t659 = pkin(3) * t806 + qJ(4) * t805 + t798 - t902;
t764 = t786 * g(1) - t789 * g(2);
t808 = qJDD(2) - t764;
t795 = -t791 * qJ(2) + t808;
t729 = -qJDD(1) * t886 + t795;
t709 = -t788 * g(3) + t785 * t729;
t790 = qJD(3) ^ 2;
t885 = pkin(3) * t785;
t807 = -qJ(4) * t788 + t885;
t799 = t791 * t807;
t663 = -t790 * pkin(3) + qJDD(3) * qJ(4) - t785 * t799 + t709;
t887 = 2 * qJD(4);
t572 = -t783 * t659 + t782 * t663 + t748 * t887;
t826 = t747 * t841;
t693 = -t728 + t826;
t545 = pkin(4) * t899 + pkin(8) * t693 - t572;
t573 = t782 * t659 + t783 * t663 + t747 * t887;
t745 = t747 ^ 2;
t802 = pkin(4) * t841 - pkin(8) * t748;
t554 = -t745 * pkin(4) + t727 * pkin(8) - t802 * t841 + t573;
t484 = t784 * t545 + t787 * t554;
t644 = pkin(5) * t702 - qJ(6) * t704;
t804 = t749 * qJ(6) - t702 * t644 + t772 * t921 + t484;
t580 = t704 * t850 - t784 * t794;
t581 = -t787 * t794 - t671;
t520 = t580 * t783 + t581 * t782;
t523 = -t580 * t782 + t581 * t783;
t801 = -t523 * t785 + t832;
t895 = t786 * t520 + t789 * t801;
t894 = t789 * t520 - t786 * t801;
t746 = t748 ^ 2;
t884 = pkin(5) * t787;
t883 = t627 * pkin(5);
t882 = qJ(6) * t787;
t881 = qJDD(1) * pkin(1);
t483 = -t787 * t545 + t784 * t554;
t437 = -t483 * t787 + t484 * t784;
t880 = t437 * t782;
t879 = t437 * t783;
t708 = t785 * g(3) + t788 * t729;
t662 = qJDD(3) * pkin(3) + t790 * qJ(4) - t788 * t799 - qJDD(4) + t708;
t603 = t727 * pkin(4) + t745 * pkin(8) - t748 * t802 + t662;
t872 = t603 * t784;
t871 = t603 * t787;
t867 = t662 * t782;
t866 = t662 * t783;
t695 = -t712 + t753;
t864 = t695 * t782;
t863 = t695 * t783;
t859 = t753 * t788;
t780 = t785 ^ 2;
t781 = t788 ^ 2;
t844 = t780 + t781;
t756 = t844 * qJDD(1);
t858 = t756 * t786;
t857 = t756 * t789;
t828 = t785 * t788 * t791;
t761 = qJDD(3) + t828;
t856 = t761 * t785;
t855 = t761 * t788;
t762 = qJDD(3) - t828;
t854 = t762 * t785;
t853 = t762 * t788;
t852 = t772 * t704;
t849 = t780 * t791;
t848 = t781 * t791;
t724 = t803 - 0.2e1 * t838 + t902;
t847 = t785 * t724;
t846 = t788 * t724;
t845 = -t612 - t888;
t836 = qJDD(1) * t786;
t834 = qJDD(1) * t789;
t830 = t785 * t712;
t829 = t788 * t712;
t827 = pkin(3) * t788 + pkin(2);
t825 = t748 * t841;
t824 = t782 * t841;
t823 = t783 * t841;
t820 = -qJ(6) * t784 - pkin(4);
t438 = t483 * t784 + t787 * t484;
t514 = t572 * t782 + t783 * t573;
t733 = -t791 * pkin(1) + t798;
t734 = -t795 + t881;
t675 = t789 * t733 - t734 * t786;
t714 = -t764 * t786 - t789 * t765;
t817 = t786 * t828;
t816 = t789 * t828;
t814 = t704 * t644 + qJDD(6) + t483;
t757 = -t786 * t791 + t834;
t813 = pkin(6) * t757 + g(3) * t786;
t758 = t789 * t791 + t836;
t812 = -pkin(6) * t758 + g(3) * t789;
t811 = t788 * t523 + t833;
t513 = -t572 * t783 + t573 * t782;
t657 = t788 * t708 + t785 * t709;
t658 = -t708 * t785 + t709 * t788;
t672 = t733 * t786 + t734 * t789;
t713 = t764 * t789 - t765 * t786;
t797 = -t749 * pkin(5) + t814;
t691 = t727 + t825;
t793 = -pkin(5) * t852 + t704 * t921 + t603;
t792 = t793 + t936;
t770 = -t790 - t848;
t769 = t790 - t848;
t768 = -t790 - t849;
t767 = -t790 + t849;
t760 = (-t780 + t781) * t791;
t759 = t844 * t791;
t755 = -0.2e1 * t822 + t835;
t752 = 0.2e1 * t821 + t837;
t750 = t844 * t839;
t732 = -t746 - t849;
t731 = -t746 + t849;
t730 = t745 - t849;
t726 = -t754 * t785 - t781 * t839;
t725 = -t780 * t839 + t859;
t721 = -t770 * t785 - t855;
t720 = t768 * t788 - t854;
t719 = t770 * t788 - t856;
t718 = -t769 * t788 - t854;
t717 = t768 * t785 + t853;
t716 = -t767 * t785 - t855;
t711 = -t759 * t789 - t858;
t710 = -t759 * t786 + t857;
t707 = -t746 + t745;
t706 = t752 * t785 - t755 * t788;
t705 = -t849 - t745;
t692 = t728 + t826;
t690 = t727 - t825;
t684 = t745 + t746;
t682 = (t747 * t783 + t748 * t782) * t841;
t681 = (t747 * t782 - t748 * t783) * t841;
t677 = t719 * t786 + t755 * t789;
t676 = t717 * t786 + t752 * t789;
t674 = -t719 * t789 + t755 * t786;
t673 = -t717 * t789 + t752 * t786;
t669 = t728 * t783 - t748 * t824;
t668 = t728 * t782 + t748 * t823;
t667 = -t727 * t782 - t747 * t823;
t666 = t727 * t783 - t747 * t824;
t660 = -t682 * t785 + t859;
t654 = t730 * t783 - t864;
t653 = -t732 * t782 - t863;
t652 = -t731 * t782 + t918;
t651 = t730 * t782 + t863;
t650 = t732 * t783 - t864;
t649 = t731 * t783 + t919;
t643 = -pkin(2) * t759 - t658;
t640 = t705 * t783 - t919;
t639 = t705 * t782 + t918;
t638 = pkin(2) * t719 - qJ(2) * t721 - t709;
t637 = pkin(2) * t717 - qJ(2) * t720 + t708;
t636 = -t669 * t785 - t829;
t635 = -t667 * t785 + t829;
t626 = t691 * t783 - t693 * t782;
t625 = t690 * t783 - t692 * t782;
t624 = t691 * t782 + t693 * t783;
t623 = t690 * t782 + t692 * t783;
t617 = pkin(2) * t752 - t720 * t886 - t846;
t616 = pkin(2) * t755 - t721 * t886 + t847;
t615 = t657 * t786 - t724 * t789;
t614 = -t657 * t789 - t724 * t786;
t609 = t653 * t788 + t692 * t785;
t608 = -t654 * t785 + t691 * t788;
t607 = t653 * t785 - t692 * t788;
t606 = -t652 * t785 - t693 * t788;
t605 = -t625 * t785 - t707 * t788;
t604 = -qJ(4) * t650 - t866;
t602 = t640 * t788 - t690 * t785;
t601 = t640 * t785 + t690 * t788;
t594 = -t627 + t852;
t590 = pkin(2) * t657 - qJ(2) * t658;
t575 = t626 * t788 - t684 * t785;
t574 = t626 * t785 + t684 * t788;
t567 = -qJ(4) * t639 - t867;
t565 = -pkin(2) * t724 - t658 * t886;
t556 = t607 * t786 + t650 * t789;
t555 = -t607 * t789 + t650 * t786;
t551 = -pkin(3) * t650 + t573;
t548 = t601 * t786 + t639 * t789;
t547 = -t601 * t789 + t639 * t786;
t546 = -pkin(3) * t639 + t572;
t542 = t574 * t786 + t624 * t789;
t541 = -t574 * t789 + t624 * t786;
t532 = -t871 + t972;
t508 = -t872 - t937;
t503 = t514 * t788 - t662 * t785;
t502 = t514 * t785 + t662 * t788;
t497 = -qJ(4) * t624 - t513;
t496 = pkin(2) * t607 - pkin(3) * t692 - qJ(2) * t609 + qJ(4) * t653 - t867;
t495 = -pkin(4) * t909 - t872 + t971;
t490 = t594 * t785 + t949;
t488 = -t594 * t788 + t951;
t486 = t792 + t883;
t485 = -pkin(4) * t592 + t871 + t938;
t481 = pkin(2) * t601 + pkin(3) * t690 - qJ(2) * t602 + qJ(4) * t640 + t866;
t480 = t592 * t785 + t949;
t478 = -t592 * t788 + t951;
t468 = (-t594 + t627) * pkin(5) + t792;
t467 = t793 + t883 + 0.2e1 * t936;
t466 = qJ(6) * t888 - t797;
t465 = -pkin(5) * t888 + t804;
t462 = t502 * t786 + t513 * t789;
t461 = -t502 * t789 + t513 * t786;
t460 = pkin(2) * t650 - t788 * t551 - t785 * t604 - t609 * t886;
t455 = qJ(6) * t845 + t797;
t454 = pkin(5) * t845 + t804;
t453 = pkin(2) * t574 + pkin(3) * t684 - qJ(2) * t575 + qJ(4) * t626 + t514;
t452 = pkin(2) * t639 - t788 * t546 - t785 * t567 - t602 * t886;
t451 = t488 * t786 + t948;
t449 = -t488 * t789 + t950;
t446 = t478 * t786 + t948;
t444 = -t478 * t789 + t950;
t442 = -t468 * t784 - t594 * t882 - t937;
t441 = -t785 * t497 - t575 * t886 + t624 * t827;
t440 = -pkin(5) * t874 + t467 * t787 - t972;
t439 = t484 - t985;
t436 = t483 + t947;
t435 = t787 * t468 + t594 * t820 + t938;
t434 = -t495 * t782 + t532 * t783 + t992;
t433 = -t971 + t784 * t467 + (pkin(4) + t884) * t909;
t432 = pkin(4) * t603 + pkin(8) * t438;
t431 = pkin(2) * t502 + pkin(3) * t662 - qJ(2) * t503 + qJ(4) * t514;
t430 = -pkin(5) * t897 + qJ(6) * t593 + t447;
t429 = -t485 * t782 + t508 * t783 - t955;
t428 = (-t896 - t888) * qJ(6) + (-t901 - t749) * pkin(5) + t814 + t947;
t427 = -t437 - t957;
t426 = t465 * t787 - t466 * t784;
t425 = t465 * t784 + t466 * t787;
t420 = t438 + t946;
t419 = -qJ(6) * t900 + (-t642 + t888) * pkin(5) - t804 + t985;
t418 = -t886 * t503 + (qJ(4) * t785 + t827) * t513;
t417 = -t454 * t784 + t455 * t787 - t957;
t416 = t454 * t787 + t455 * t784 + t946;
t415 = t438 * t783 - t880;
t414 = t438 * t782 + t879;
t413 = t415 * t788 - t603 * t785;
t412 = t415 * t785 + t603 * t788;
t411 = t495 * t783 + t532 * t782 - t997;
t410 = -pkin(8) * t425 + (-pkin(5) * t784 + t882) * t486;
t409 = -t435 * t782 + t442 * t783 - t955;
t408 = pkin(2) * t478 - pkin(3) * t592 - qJ(2) * t480 + t485 * t783 + t508 * t782 + t954;
t407 = -t433 * t782 + t440 * t783 - t992;
t406 = -t425 * t782 + t426 * t783;
t405 = t425 * t783 + t426 * t782;
t404 = pkin(8) * t426 + (-t820 + t884) * t486;
t403 = -pkin(3) * t414 - pkin(4) * t437;
t402 = t406 * t788 - t486 * t785;
t401 = t406 * t785 + t486 * t788;
t400 = -t420 * t782 + t427 * t783 - t969;
t399 = -t785 * t434 - t788 * t439 - t1000;
t398 = pkin(2) * t488 - pkin(3) * t594 - qJ(2) * t490 + t435 * t783 + t442 * t782 + t954;
t397 = t433 * t783 + t440 * t782 + t997;
t396 = -t785 * t429 - t788 * t436 - t480 * t886 + t959;
t395 = -pkin(8) * t879 - qJ(4) * t414 - t432 * t782;
t394 = t412 * t786 + t414 * t789;
t393 = -t412 * t789 + t414 * t786;
t392 = -t416 * t782 + t417 * t783 - t969;
t391 = t420 * t783 + t427 * t782 + t975;
t390 = -t785 * t409 - t788 * t428 - t490 * t886 + t959;
t389 = -pkin(3) * t405 - pkin(4) * t425 - pkin(5) * t466 - qJ(6) * t465;
t388 = -t785 * t407 - t788 * t419 + t1000;
t387 = t401 * t786 + t405 * t789;
t386 = -t401 * t789 + t405 * t786;
t385 = t416 * t783 + t417 * t782 + t975;
t384 = -t785 * t400 - t788 * t447 + t977;
t383 = -qJ(4) * t405 - t404 * t782 + t410 * t783;
t382 = -t785 * t392 - t788 * t430 + t977;
t381 = pkin(2) * t412 + pkin(3) * t603 - pkin(8) * t880 - qJ(2) * t413 + qJ(4) * t415 + t432 * t783;
t380 = pkin(2) * t414 - t785 * t395 - t788 * t403 - t413 * t886;
t379 = pkin(2) * t401 + pkin(3) * t486 - qJ(2) * t402 + qJ(4) * t406 + t404 * t783 + t410 * t782;
t378 = pkin(2) * t405 - t785 * t383 - t788 * t389 - t402 * t886;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t758, -t757, 0, t714, 0, 0, 0, 0, 0, 0, 0, t758, t757, t675, 0, 0, 0, 0, 0, 0, t676, t677, t711, t615, 0, 0, 0, 0, 0, 0, t548, t556, t542, t462, 0, 0, 0, 0, 0, 0, t446, t450, t981, t394, 0, 0, 0, 0, 0, 0, t451, t981, -t450, t387; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t757, -t758, 0, t713, 0, 0, 0, 0, 0, 0, 0, -t757, t758, t672, 0, 0, 0, 0, 0, 0, t673, t674, t710, t614, 0, 0, 0, 0, 0, 0, t547, t555, t541, t461, 0, 0, 0, 0, 0, 0, t444, -t443, t979, t393, 0, 0, 0, 0, 0, 0, t449, t979, t443, t386; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t720, t721, 0, t658, 0, 0, 0, 0, 0, 0, t602, t609, t575, t503, 0, 0, 0, 0, 0, 0, t480, t489, t960, t413, 0, 0, 0, 0, 0, 0, t490, t960, -t489, t402; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t757, 0, -t758, 0, -t813, -t812, -t713, -pkin(6) * t713, 0, -t757, t758, 0, 0, 0, -t672, t813, t812, -pkin(6) * t672 + (-pkin(1) * t786 + qJ(2) * t789) * g(3), -t726 * t786 + t816, -t706 * t786 + t760 * t789, -t718 * t786 + t788 * t834, -t725 * t786 - t816, -t716 * t786 - t785 * t834, qJDD(3) * t789 - t750 * t786, -pkin(6) * t673 - t617 * t786 + t637 * t789, -pkin(6) * t674 - t616 * t786 + t638 * t789, -pkin(2) * t857 - pkin(6) * t710 - t643 * t786, -pkin(6) * t614 - t565 * t786 + t590 * t789, -t636 * t786 + t668 * t789, -t605 * t786 + t623 * t789, -t606 * t786 + t649 * t789, -t635 * t786 + t666 * t789, -t608 * t786 + t651 * t789, -t660 * t786 + t681 * t789, -pkin(6) * t547 - t452 * t786 + t481 * t789, -pkin(6) * t555 - t460 * t786 + t496 * t789, -pkin(6) * t541 - t441 * t786 + t453 * t789, -pkin(6) * t461 - t418 * t786 + t431 * t789, t894, t1003, t980, t926, -t999, t927, -pkin(6) * t444 - t396 * t786 + t408 * t789, -t399 * t786 + t411 * t789 + t1005, -t384 * t786 + t391 * t789 - t995, -pkin(6) * t393 - t380 * t786 + t381 * t789, t894, t980, -t1003, t927, t999, t926, -pkin(6) * t449 - t390 * t786 + t398 * t789, -t382 * t786 + t385 * t789 - t995, -t388 * t786 + t397 * t789 - t1005, -pkin(6) * t386 - t378 * t786 + t379 * t789; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t758, 0, t757, 0, t812, -t813, t714, pkin(6) * t714, 0, -t758, -t757, 0, 0, 0, t675, -t812, t813, pkin(6) * t675 + (pkin(1) * t789 + qJ(2) * t786) * g(3), t726 * t789 + t817, t706 * t789 + t760 * t786, t718 * t789 + t786 * t835, t725 * t789 - t817, t716 * t789 - t785 * t836, qJDD(3) * t786 + t750 * t789, pkin(6) * t676 + t617 * t789 + t637 * t786, pkin(6) * t677 + t616 * t789 + t638 * t786, -pkin(2) * t858 + pkin(6) * t711 + t643 * t789, pkin(6) * t615 + t565 * t789 + t590 * t786, t636 * t789 + t668 * t786, t605 * t789 + t623 * t786, t606 * t789 + t649 * t786, t635 * t789 + t666 * t786, t608 * t789 + t651 * t786, t660 * t789 + t681 * t786, pkin(6) * t548 + t452 * t789 + t481 * t786, pkin(6) * t556 + t460 * t789 + t496 * t786, pkin(6) * t542 + t441 * t789 + t453 * t786, pkin(6) * t462 + t418 * t789 + t431 * t786, t895, -t1002, t978, t924, t998, t925, pkin(6) * t446 + t396 * t789 + t408 * t786, t399 * t789 + t411 * t786 + t1004, t384 * t789 + t391 * t786 + t994, pkin(6) * t394 + t380 * t789 + t381 * t786, t895, t978, t1002, t925, -t998, t924, pkin(6) * t451 + t390 * t789 + t398 * t786, t382 * t789 + t385 * t786 + t994, t388 * t789 + t397 * t786 - t1004, pkin(6) * t387 + t378 * t789 + t379 * t786; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t764, t765, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t808 - 0.2e1 * t881, -t765 + t777 + 0.2e1 * t779, pkin(1) * t734 + qJ(2) * t733, -t805 * t788, -t752 * t788 - t755 * t785, -t769 * t785 + t853, t806 * t785, t767 * t788 - t856, 0, qJ(2) * t752 - t717 * t886 - t847, qJ(2) * t755 - t719 * t886 - t846, -qJ(2) * t759 + t756 * t886 - t657, -qJ(2) * t724 - t657 * t886, t669 * t788 - t830, t625 * t788 - t707 * t785, t652 * t788 - t693 * t785, t667 * t788 + t830, t654 * t788 + t691 * t785, t682 * t788 + t753 * t785, qJ(2) * t639 - t785 * t546 + t788 * t567 - t601 * t886, qJ(2) * t650 - t785 * t551 + t788 * t604 - t607 * t886, t788 * t497 + (qJ(2) + t885) * t624 - t886 * t574, -t886 * t502 + (qJ(2) + t807) * t513, t811, t983, t961, t904, -t982, t903, t788 * t429 - t785 * t436 - t478 * t886 + t956, t788 * t434 - t785 * t439 - t1001, t788 * t400 - t785 * t447 + t976, qJ(2) * t414 + t788 * t395 - t785 * t403 - t412 * t886, t811, t961, -t983, t903, t982, t904, t788 * t409 - t785 * t428 - t488 * t886 + t956, t788 * t392 - t785 * t430 + t976, t788 * t407 - t785 * t419 + t1001, qJ(2) * t405 + t788 * t383 - t785 * t389 - t401 * t886;];
tauB_reg  = t1;
