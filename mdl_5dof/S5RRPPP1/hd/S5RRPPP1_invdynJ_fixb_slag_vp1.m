% Calculate vector of inverse dynamics joint torques for
% S5RRPPP1
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:21
% EndTime: 2019-12-31 19:24:29
% DurationCPUTime: 53.76s
% Computational Cost: add. (17486->903), mult. (51357->1110), div. (0->0), fcn. (54081->8), ass. (0->392)
t877 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t876 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t875 = Icges(4,4) + Icges(5,6) - Icges(6,6);
t874 = Icges(5,4) - Icges(4,5) - Icges(6,5);
t873 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t872 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t554 = sin(pkin(5));
t555 = sin(qJ(2));
t557 = cos(qJ(2));
t778 = cos(pkin(8));
t779 = cos(pkin(5));
t633 = t779 * t778;
t777 = sin(pkin(8));
t448 = t555 * t633 + t557 * t777;
t558 = cos(qJ(1));
t566 = t448 * t558;
t556 = sin(qJ(1));
t664 = t556 * t778;
t352 = -t554 * t664 + t566;
t661 = t558 * t778;
t663 = t556 * t777;
t567 = t554 * t663 + t557 * t661;
t632 = t779 * t777;
t595 = t555 * t632;
t353 = -t558 * t595 + t567;
t754 = t555 * t558;
t694 = t554 * t754;
t451 = t556 * t779 + t694;
t350 = t448 * t556 + t554 * t661;
t449 = t557 * t778 - t595;
t660 = t558 * t777;
t641 = t554 * t660;
t351 = t449 * t556 - t641;
t755 = t555 * t556;
t450 = t554 * t755 - t558 * t779;
t826 = t872 * t350 - t875 * t351 + t873 * t450;
t828 = -t875 * t350 + t877 * t351 - t874 * t450;
t854 = t873 * t350 - t874 * t351 + t876 * t450;
t899 = -t352 * t826 - t353 * t828 - t451 * t854;
t858 = t350 * t826 + t351 * t828 + t450 * t854;
t827 = -t875 * t352 + t877 * t353 - t874 * t451;
t825 = t872 * t352 - t875 * t353 + t873 * t451;
t823 = t873 * t352 - t874 * t353 + t876 * t451;
t594 = t557 * t633;
t446 = t555 * t777 - t594;
t447 = t555 * t778 + t557 * t632;
t756 = t554 * t557;
t853 = -t873 * t446 + t874 * t447 + t876 * t756;
t814 = t872 * t446 - t875 * t447 - t873 * t756;
t812 = -t875 * t446 + t877 * t447 + t874 * t756;
t547 = Icges(3,4) * t557;
t612 = -Icges(3,2) * t555 + t547;
t412 = Icges(3,6) * t556 + t558 * t612;
t775 = Icges(3,4) * t555;
t508 = Icges(3,1) * t557 - t775;
t414 = Icges(3,5) * t556 + t508 * t558;
t753 = t556 * t557;
t376 = t414 * t753;
t504 = Icges(3,5) * t557 - Icges(3,6) * t555;
t410 = Icges(3,3) * t556 + t504 * t558;
t655 = t410 * t558 - t376;
t190 = -t412 * t755 - t655;
t857 = t825 * t350 + t827 * t351 + t823 * t450;
t807 = t190 + t857;
t752 = t557 * t558;
t727 = t556 * t410 + t414 * t752;
t192 = -t412 * t754 + t727;
t856 = t825 * t352 + t827 * t353 + t823 * t451;
t805 = t192 + t856;
t868 = rSges(6,1) + pkin(4);
t757 = t554 * t555;
t833 = t557 * pkin(2) + qJ(3) * t757;
t587 = -pkin(1) - t833;
t865 = t556 * t587;
t898 = t899 * t558;
t897 = t858 * t558;
t505 = Icges(3,2) * t557 + t775;
t507 = Icges(3,1) * t555 + t547;
t603 = t505 * t555 - t507 * t557;
t503 = Icges(3,5) * t555 + Icges(3,6) * t557;
t758 = t503 * t558;
t846 = t814 * t350 + t812 * t351 - t853 * t450 - t556 * t603 - t758;
t759 = t503 * t556;
t845 = t814 * t352 + t812 * t353 - t853 * t451 - t558 * t603 + t759;
t713 = qJD(1) * t350;
t832 = t446 * qJD(2);
t225 = t558 * t832 + t713;
t406 = t447 * t558;
t583 = qJD(1) * t595;
t639 = qJD(1) * t664;
t226 = -qJD(1) * t641 + qJD(2) * t406 - t556 * t583 + t557 * t639;
t706 = qJD(2) * t558;
t678 = t557 * t706;
t491 = t554 * t678;
t652 = -qJD(1) * t450 + t491;
t883 = -t225 * t872 + t226 * t875 + t652 * t873;
t227 = qJD(1) * t566 - t554 * t639 - t556 * t832;
t592 = qJD(2) * t632;
t657 = qJD(2) * t778;
t228 = qJD(1) * t567 - t558 * t583 - t592 * t753 - t657 * t755;
t696 = t554 * t753;
t834 = qJD(1) * t451 + qJD(2) * t696;
t882 = -t227 * t872 + t228 * t875 - t834 * t873;
t881 = t225 * t875 - t226 * t877 - t652 * t874;
t880 = t227 * t875 - t228 * t877 + t834 * t874;
t879 = -t225 * t873 + t226 * t874 + t652 * t876;
t893 = t227 * t873 - t228 * t874 + t834 * t876;
t417 = t448 * qJD(2);
t418 = -t555 * t592 + t557 * t657;
t708 = qJD(2) * t555;
t679 = t554 * t708;
t892 = t417 * t873 - t418 * t874 + t679 * t876;
t891 = t417 * t872 - t418 * t875 + t679 * t873;
t890 = -t417 * t875 + t418 * t877 - t679 * t874;
t772 = Icges(3,6) * t558;
t411 = Icges(3,4) * t753 - Icges(3,2) * t755 - t772;
t771 = Icges(3,3) * t558;
t409 = Icges(3,5) * t753 - Icges(3,6) * t755 - t771;
t537 = Icges(3,4) * t755;
t774 = Icges(3,5) * t558;
t413 = Icges(3,1) * t753 - t537 - t774;
t728 = -t556 * t409 - t413 * t752;
t191 = -t411 * t754 - t728;
t889 = -t191 * t558 + t805 * t556 + t898;
t761 = t411 * t555;
t607 = -t413 * t557 + t761;
t189 = -t409 * t558 - t556 * t607;
t888 = -t189 * t558 + t556 * t807 - t897;
t866 = rSges(6,3) + qJ(5);
t886 = -t350 * rSges(6,2) - t868 * t450;
t658 = t779 * qJ(3);
t530 = t558 * t658;
t379 = t556 * t833 - t530;
t432 = qJD(3) * t451;
t552 = t558 * pkin(7);
t519 = pkin(1) * t556 - t552;
t482 = qJD(1) * t519;
t709 = qJD(1) * t558;
t543 = pkin(7) * t709;
t704 = qJD(3) * t555;
t527 = t554 * t704;
t656 = qJD(3) * t779;
t643 = qJ(3) * t491 + qJD(1) * t530 + t558 * t527 + t556 * t656;
t675 = t555 * t706;
t885 = -pkin(2) * t675 + qJD(1) * t379 - t432 + t482 + t543 + t643;
t540 = pkin(2) * t752;
t380 = qJ(3) * t694 + t556 * t658 + t540;
t468 = pkin(2) * t555 - qJ(3) * t756;
t707 = qJD(2) * t556;
t520 = t558 * pkin(1) + t556 * pkin(7);
t859 = qJD(1) * t520;
t884 = -qJD(1) * t380 - qJD(3) * t450 + t468 * t707 - t859;
t871 = t845 * qJD(1);
t870 = t846 * qJD(1);
t700 = qJD(1) * qJD(2);
t478 = qJDD(2) * t556 + t558 * t700;
t213 = t353 * pkin(3) + qJ(4) * t352;
t710 = qJD(1) * t556;
t581 = -t557 * t710 - t675;
t681 = t555 * t710;
t209 = -qJ(3) * t554 * t681 + pkin(2) * t581 + t643;
t720 = qJD(1) * (-pkin(1) * t710 + t543) + qJDD(1) * t520;
t591 = qJD(1) * t209 + qJD(3) * t834 + qJDD(1) * t380 + qJDD(3) * t450 + t720;
t332 = qJD(4) * t352;
t89 = -t226 * pkin(3) - qJ(4) * t225 + t332;
t569 = qJD(1) * t89 + qJD(4) * t227 + qJDD(1) * t213 + qJDD(4) * t350 + t591;
t702 = qJD(5) * t447;
t676 = qJD(3) * t756;
t398 = qJD(2) * t833 - t676;
t420 = qJD(4) * t446;
t743 = -pkin(3) * t418 - qJ(4) * t417 - t398 - t420;
t647 = -rSges(6,2) * t417 - t418 * t866 - t679 * t868 - t702 + t743;
t599 = qJD(2) * t647;
t322 = pkin(3) * t447 + qJ(4) * t446;
t733 = -t322 - t468;
t737 = rSges(6,2) * t446 + t447 * t866 - t756 * t868;
t644 = t733 - t737;
t746 = t352 * rSges(6,2) + t353 * t866 + t451 * t868;
t330 = qJD(5) * t353;
t751 = -t225 * rSges(6,2) - t226 * t866 + t652 * t868 + t330;
t3 = qJD(1) * t751 + qJD(5) * t228 + qJDD(1) * t746 + qJDD(5) * t351 + t478 * t644 + t556 * t599 + t569;
t869 = -g(2) + t3;
t333 = t350 * qJ(4);
t212 = -t351 * pkin(3) - t333;
t864 = -qJD(1) * t212 - t332 + t885 + t89;
t703 = qJD(4) * t350;
t863 = -qJD(1) * t213 + t322 * t707 - t703 + t884;
t732 = -t449 * pkin(3) - qJ(4) * t448 - t833;
t862 = rSges(5,1) * t757 - t449 * rSges(5,2) + t448 * rSges(5,3) - t732;
t645 = -t448 * rSges(6,2) - t866 * t449 - t757 * t868 + t732;
t403 = t555 * t663 - t556 * t594;
t404 = t447 * t556;
t248 = rSges(5,1) * t696 + t404 * rSges(5,2) - t403 * rSges(5,3);
t279 = -t404 * pkin(3) - qJ(4) * t403;
t698 = pkin(2) * t755;
t434 = qJ(3) * t696 - t698;
t738 = -t279 - t434;
t861 = -t738 + t248;
t741 = -t403 * rSges(6,2) - t404 * t866 + t696 * t868;
t860 = t738 - t741;
t747 = -t351 * t866 + t886;
t852 = qJD(2) * t888 + t870;
t851 = qJD(2) * t889 + t871;
t585 = qJD(2) * t505;
t316 = qJD(1) * t412 - t556 * t585;
t586 = qJD(2) * t507;
t318 = qJD(1) * t414 - t556 * t586;
t850 = t607 * qJD(2) - t316 * t557 - t318 * t555 + (t557 * t893 - t708 * t854) * t554 + t880 * t447 + t882 * t446 - t828 * t418 - t826 * t417;
t315 = -t558 * t585 + (-t556 * t612 + t772) * qJD(1);
t317 = -t558 * t586 + (-t508 * t556 + t774) * qJD(1);
t760 = t412 * t555;
t606 = -t414 * t557 + t760;
t849 = -qJD(2) * t606 + t315 * t557 + t317 * t555 + (-t557 * t879 + t708 * t823) * t554 + t881 * t447 + t883 * t446 + t827 * t418 + t825 * t417;
t466 = t612 * qJD(2);
t467 = t508 * qJD(2);
t604 = t505 * t557 + t507 * t555;
t562 = qJD(1) * t503 - qJD(2) * t604 - t466 * t555 + t467 * t557;
t800 = t603 * qJD(1) + t504 * qJD(2);
t848 = -t225 * t814 - t226 * t812 + t352 * t891 + t353 * t890 + t451 * t892 + t556 * t800 + t562 * t558 - t652 * t853;
t847 = t227 * t814 + t228 * t812 + t350 * t891 + t351 * t890 + t450 * t892 + t562 * t556 - t558 * t800 - t834 * t853;
t806 = t191 - t899;
t281 = t411 * t557 + t413 * t555;
t804 = t446 * t826 + t447 * t828 - t756 * t854 + t281;
t282 = t412 * t557 + t414 * t555;
t803 = t446 * t825 + t447 * t827 - t756 * t823 + t282;
t329 = qJD(5) * t351;
t844 = t746 * qJD(1) + t329 - t863;
t715 = t530 + t552;
t843 = t715 + t865;
t716 = qJD(2) * t698 + t558 * t656;
t842 = t587 * t709 + t716 + ((-t658 - pkin(7)) * qJD(1) + (-qJ(3) * qJD(2) * t557 - t704) * t554) * t556;
t840 = -(t350 * t556 + t352 * t558) * qJD(2) + t417;
t839 = t446 * t706 - t225 + t713;
t838 = -qJD(1) * t352 + t446 * t707 + t227;
t251 = -t404 * rSges(4,1) + t403 * rSges(4,2) + rSges(4,3) * t696;
t835 = t434 + t251;
t831 = -t227 * rSges(6,2) - t834 * t868 - t329;
t180 = t353 * rSges(4,1) - t352 * rSges(4,2) + t451 * rSges(4,3);
t829 = -qJD(1) * t180 + t884;
t822 = -t403 * t875 + t404 * t877 + t696 * t874;
t405 = t555 * t660 - t558 * t594;
t695 = t554 * t752;
t821 = t405 * t875 - t406 * t877 - t695 * t874;
t820 = t403 * t872 - t404 * t875 - t696 * t873;
t819 = -t405 * t872 + t406 * t875 + t695 * t873;
t818 = -t403 * t873 + t404 * t874 + t696 * t876;
t817 = -t405 * t873 + t406 * t874 + t695 * t876;
t816 = t448 * t873 - t449 * t874 + t757 * t876;
t813 = t448 * t872 - t449 * t875 + t757 * t873;
t811 = -t448 * t875 + t449 * t877 - t757 * t874;
t178 = t351 * rSges(4,1) - t350 * rSges(4,2) + t450 * rSges(4,3);
t184 = -t450 * rSges(5,1) + t351 * rSges(5,2) - t350 * rSges(5,3);
t186 = t451 * rSges(5,1) - t353 * rSges(5,2) + t352 * rSges(5,3);
t809 = -qJD(1) * t186 + t863;
t808 = t189 + t858;
t584 = qJD(2) * t503;
t802 = -t558 * t584 + (-t504 * t556 + t606 + t771) * qJD(1);
t712 = qJD(1) * t410;
t801 = qJD(1) * t607 - t556 * t584 + t712;
t563 = -qJD(2) * t282 - t315 * t555 + t317 * t557 + t712;
t564 = qJD(1) * t409 - qJD(2) * t281 - t316 * t555 + t318 * t557;
t799 = (-t826 * t227 - t828 * t228 + t882 * t350 + t880 * t351 - t450 * t893 + t558 * t801 - t854 * t834) * t558 + (t563 * t556 + (-t564 - t802) * t558 + t823 * t834 + t879 * t450 + t881 * t351 + t883 * t350 + t827 * t228 + t825 * t227) * t556;
t798 = (t826 * t225 + t828 * t226 + t882 * t352 + t880 * t353 - t451 * t893 - t564 * t558 - t854 * t652) * t558 + (t556 * t802 + (t563 - t801) * t558 + t823 * t652 + t879 * t451 + t881 * t353 + t883 * t352 - t827 * t226 - t825 * t225) * t556;
t796 = -m(5) - m(6);
t795 = t478 / 0.2e1;
t479 = -qJDD(2) * t558 + t556 * t700;
t794 = t479 / 0.2e1;
t789 = rSges(3,1) * t557;
t306 = -rSges(5,1) * t756 - rSges(5,2) * t447 + rSges(5,3) * t446;
t726 = -t379 - t519;
t690 = t212 + t726;
t651 = t184 + t690;
t687 = -t306 + t733;
t731 = t332 + t432;
t44 = qJD(1) * t651 + t687 * t706 + t731;
t788 = t44 * t306;
t548 = t556 * rSges(3,3);
t307 = rSges(4,1) * t447 - rSges(4,2) * t446 - rSges(4,3) * t756;
t693 = -t178 + t726;
t736 = -t307 - t468;
t79 = qJD(1) * t693 + t706 * t736 + t432;
t781 = t79 * t307;
t511 = rSges(3,1) * t555 + rSges(3,2) * t557;
t680 = t511 * t706;
t714 = rSges(3,2) * t755 + t558 * rSges(3,3);
t415 = rSges(3,1) * t753 - t714;
t721 = -t415 - t519;
t277 = qJD(1) * t721 - t680;
t765 = t277 * t556;
t764 = t277 * t558;
t416 = rSges(3,1) * t752 - rSges(3,2) * t754 + t548;
t278 = qJD(1) * t416 - t511 * t707 + t859;
t459 = t511 * t558;
t763 = t278 * t459;
t750 = t228 * t866 - t831;
t749 = -t180 - t380;
t210 = qJ(3) * t834 + qJD(1) * t540 + t556 * t527 - t716;
t745 = -t210 - t859;
t744 = -t213 - t380;
t742 = t227 * qJ(4) + t703;
t740 = -t405 * rSges(6,2) - t406 * t866 + t695 * t868;
t739 = -rSges(4,1) * t418 + rSges(4,2) * t417 - rSges(4,3) * t679 - t398;
t735 = -t449 * rSges(4,1) + t448 * rSges(4,2) - rSges(4,3) * t757 - t833;
t445 = t468 * t710;
t734 = t322 * t710 + t445;
t729 = t556 * t379 + t558 * t380;
t493 = t558 * t676;
t725 = -qJD(4) * t405 + t493;
t435 = -pkin(2) * t754 + qJ(3) * t695;
t724 = qJD(1) * t435 + t556 * t676;
t723 = -Icges(3,2) * t753 + t413 - t537;
t722 = -t505 * t558 + t414;
t719 = -t505 + t508;
t718 = t507 + t612;
t717 = rSges(3,2) * t681 + rSges(3,3) * t709;
t711 = qJD(1) * t504;
t699 = -pkin(3) - t866;
t90 = t228 * pkin(3) + t742;
t697 = -t90 + t745;
t692 = -t186 + t744;
t691 = t558 * t209 + t556 * t210 + t379 * t709;
t689 = -rSges(5,1) * t679 + rSges(5,2) * t418 - rSges(5,3) * t417 + t743;
t109 = -t226 * rSges(4,1) + t225 * rSges(4,2) + rSges(4,3) * t652;
t686 = qJD(3) * t652 + qJDD(3) * t451 + t479 * t468;
t112 = rSges(5,1) * t652 + t226 * rSges(5,2) - t225 * rSges(5,3);
t252 = -t406 * rSges(4,1) + t405 * rSges(4,2) + rSges(4,3) * t695;
t685 = t434 * t707 + t435 * t706 + t527;
t250 = rSges(5,1) * t695 + t406 * rSges(5,2) - t405 * rSges(5,3);
t674 = -pkin(1) - t789;
t671 = -t707 / 0.2e1;
t670 = t707 / 0.2e1;
t669 = -t706 / 0.2e1;
t668 = t706 / 0.2e1;
t617 = t690 + t747;
t41 = qJD(1) * t617 + t644 * t706 + t330 + t731;
t666 = t41 * t737;
t80 = -t307 * t707 - t829;
t659 = t80 * t736;
t280 = -t406 * pkin(3) - qJ(4) * t405;
t654 = -t409 + t760;
t653 = qJD(2) * t739;
t650 = t744 - t746;
t648 = -t212 * t556 + t558 * t213 + t729;
t646 = qJD(1) * t280 - qJD(4) * t403 + t724;
t45 = -t306 * t707 - t809;
t642 = t45 * t687;
t636 = qJD(2) * t689;
t514 = rSges(2,1) * t558 - rSges(2,2) * t556;
t512 = rSges(2,1) * t556 + rSges(2,2) * t558;
t513 = -rSges(3,2) * t555 + t789;
t614 = t380 + t520;
t609 = -t278 * t556 - t764;
t319 = rSges(3,1) * t581 - rSges(3,2) * t678 + t717;
t458 = t511 * t556;
t320 = -qJD(2) * t458 + (t513 * t558 + t548) * qJD(1);
t608 = t319 * t558 + t320 * t556;
t605 = t415 * t556 + t416 * t558;
t42 = -t707 * t737 + t844;
t602 = t42 * t644;
t601 = t212 + t715;
t598 = -t212 * t709 + t556 * t90 + t558 * t89 + t691;
t597 = t379 * t707 + t380 * t706 - t676;
t590 = qJD(4) * t448 + t279 * t707 + t280 * t706 + t685;
t589 = -qJD(4) * t225 + qJDD(4) * t352 + t479 * t322 + t686;
t110 = t228 * rSges(4,1) - t227 * rSges(4,2) + rSges(4,3) * t834;
t114 = rSges(5,1) * t834 - t228 * rSges(5,2) + t227 * rSges(5,3);
t578 = t280 + t435;
t577 = qJD(2) * t527 - qJDD(3) * t756 + t209 * t706 + t210 * t707 + t478 * t379;
t576 = t411 * t558 - t412 * t556;
t575 = t213 + t614;
t574 = -t212 * t707 + t213 * t706 + t420 + t597;
t573 = g(1) * t865;
t568 = (-t555 * t718 + t557 * t719) * qJD(1);
t565 = qJD(4) * t417 + qJDD(4) * t446 - t212 * t478 + t89 * t706 + t90 * t707 + t577;
t560 = t576 * t557 + (-t556 * t722 + t558 * t723) * t555;
t470 = t513 * qJD(2);
t310 = (t450 * t556 + t451 * t558) * qJD(2);
t275 = t605 * qJD(2);
t120 = qJD(1) * t319 + qJDD(1) * t416 - t470 * t707 - t478 * t511 + t720;
t119 = -t470 * t706 + t479 * t511 + t721 * qJDD(1) + (-t320 - t859) * qJD(1);
t68 = (t178 * t556 + t180 * t558) * qJD(2) + t597;
t43 = (-t184 * t556 + t186 * t558) * qJD(2) + t574;
t34 = t702 + (-t556 * t747 + t558 * t746) * qJD(2) + t574;
t33 = qJD(1) * t109 + qJDD(1) * t180 + t478 * t736 + t556 * t653 + t591;
t32 = t307 * t479 + t558 * t653 + t693 * qJDD(1) + (-t110 + t745) * qJD(1) + t686;
t31 = t178 * t478 + t749 * t479 + (t109 * t558 + t110 * t556) * qJD(2) + t577;
t6 = qJD(1) * t112 + qJDD(1) * t186 + t478 * t687 + t556 * t636 + t569;
t5 = t306 * t479 + t558 * t636 + t651 * qJDD(1) + (-t114 + t697) * qJD(1) + t589;
t4 = -t184 * t478 + (t112 * t558 + t114 * t556) * qJD(2) + t692 * t479 + t565;
t2 = -qJD(5) * t226 + qJDD(5) * t353 + t737 * t479 + t558 * t599 + t617 * qJDD(1) + (t697 - t750) * qJD(1) + t589;
t1 = (t556 * t750 + t558 * t751) * qJD(2) - t747 * t478 + t650 * t479 + t565 + qJD(5) * t418 + qJDD(5) * t447;
t7 = [-m(2) * (-g(1) * t512 + g(2) * t514) + (-qJD(2) * t603 + t466 * t557 + t467 * t555 + (-t557 * t892 - t853 * t708) * t554 + t890 * t447 + t891 * t446 + t812 * t418 + t814 * t417) * qJD(1) + (-(t556 * t666 + t602 * t558) * qJD(2) + (t699 * t351 - t333 + t843 + t886) * t2 - g(1) * (t601 + t747) - t573 + t869 * (t575 + t746) + (-t747 * qJD(1) + t587 * t710 - t330 + t751 + t864) * t42 + (t699 * t228 - t742 + t831 + t842 + t844) * t41) * m(6) + (-g(1) * (t601 + t184) - t573 - (t556 * t788 + t642 * t558) * qJD(2) + (t6 - g(2)) * (t575 + t186) + (t184 + t212 + t843) * t5 + (t112 + (-t184 + t865) * qJD(1) + t864) * t45 + (-t114 - t90 - t809 + t842) * t44) * m(5) + (-(t556 * t781 + t659 * t558) * qJD(2) - g(1) * (-t178 + t715) - t573 + (t109 + (t178 + t865) * qJD(1) + t885) * t80 + (-t178 + t843) * t32 + (-t110 - t829 + t842) * t79 + (t33 - g(2)) * (t614 + t180)) * m(4) + (-(-qJD(1) * t415 - t277 - t482 - t680) * t278 + t278 * (t543 + t717) + (t511 * t765 - t763) * qJD(2) + ((-pkin(1) - t513) * t764 + (t277 * (-rSges(3,3) - pkin(7)) + t278 * t674) * t556) * qJD(1) + (-g(2) + t120) * (t416 + t520) + (-g(1) + t119) * (t674 * t556 + t552 + t714)) * m(3) + (((t190 - t376 + (t410 + t761) * t558 + t728) * t558 + (t727 + t856) * t556 + t898) * qJD(2) + t871) * t668 + (t604 + m(2) * (t512 ^ 2 + t514 ^ 2) + Icges(2,3) + t853 * t756 + t812 * t447 + t814 * t446) * qJDD(1) + (t845 + t803) * t795 + (t846 + t804) * t794 + (t848 + t849) * t670 + (t847 - t850 + t851) * t669 + (((t558 * t654 + t192 - t727) * t558 + (t556 * t654 + t655 + t806 - t857 + t899) * t556 + t897) * qJD(2) + t852 - t870) * t671; (t41 * t734 + t1 * t648 + t34 * t598 + (t2 * t644 + t41 * t647 + t1 * t746 + t34 * t751 + (-t34 * t747 + t602) * qJD(1)) * t558 + (t3 * t644 + t42 * t647 - t1 * t747 + t34 * t750 + (t34 * t650 + t666) * qJD(1)) * t556 - g(1) * (t578 + t740) + g(2) * t860 + g(3) * t645 - t41 * (-qJD(5) * t406 + t725) - t42 * (-qJD(5) * t404 + t646) - t34 * (qJD(5) * t449 + t590) - (t41 * t860 + t42 * t740) * qJD(1) - ((t34 * t740 + t41 * t645) * t558 + (t34 * t741 + t42 * t645) * t556) * qJD(2)) * m(6) + (-g(1) * (t578 + t250) - g(2) * t861 - g(3) * t862 - t44 * t725 - t45 * t646 - t43 * t590 - (t45 * t250 - t44 * t861) * qJD(1) - ((t43 * t250 - t44 * t862) * t558 + (t43 * t248 - t45 * t862) * t556) * qJD(2) + t44 * t734 + t4 * t648 + t43 * t598 + (t5 * t687 + t44 * t689 + t4 * t186 + t43 * t112 + (-t184 * t43 + t642) * qJD(1)) * t558 + (t6 * t687 + t45 * t689 - t4 * t184 + t43 * t114 + (t43 * t692 + t788) * qJD(1)) * t556) * m(5) + (t79 * t445 + t31 * t729 + t68 * t691 + (t32 * t736 + t79 * t739 + t31 * t180 + t68 * t109 + (t68 * t178 + t659) * qJD(1)) * t558 + (t33 * t736 + t80 * t739 + t31 * t178 + t68 * t110 + (t68 * t749 + t781) * qJD(1)) * t556 - g(1) * (t435 + t252) - g(2) * t835 + g(3) * t735 - t79 * (-qJD(1) * t835 + t493) - t80 * (qJD(1) * t252 + t724) - t68 * t685 - ((t68 * t252 + t735 * t79) * t558 + (t68 * t251 + t735 * t80) * t556) * qJD(2)) * m(4) + (g(1) * t459 + g(2) * t458 - g(3) * t513 + (qJD(2) * t608 + t415 * t478 - t416 * t479) * t605 + t275 * ((t415 * t558 - t416 * t556) * qJD(1) + t608) + t609 * t470 + (-t119 * t558 - t120 * t556 + (-t278 * t558 + t765) * qJD(1)) * t511 - (t277 * t458 - t763) * qJD(1) - (t275 * (-t458 * t556 - t459 * t558) + t609 * t513) * qJD(2)) * m(3) + t889 * t795 + t888 * t794 + (qJD(1) * t848 + t798 * qJD(2) + qJDD(1) * t845 + t805 * t478 + t806 * t479) * t556 / 0.2e1 - (qJD(1) * t847 + t799 * qJD(2) + qJDD(1) * t846 + t807 * t478 + t808 * t479) * t558 / 0.2e1 - ((t576 * t555 + (t446 * t820 + t447 * t822 - t448 * t826 - t449 * t828 - t723 * t557) * t558 + (t446 * t819 + t447 * t821 + t448 * t825 + t449 * t827 + t722 * t557) * t556 + ((-t555 * t854 + t557 * t818) * t558 + (t555 * t823 - t557 * t817) * t556) * t554) * qJD(2) + (t555 * t719 + t557 * t718 + (-t555 * t853 - t557 * t816) * t554 + t812 * t449 + t814 * t448 + t811 * t447 + t813 * t446) * qJD(1)) * qJD(1) / 0.2e1 + (t850 * t558 + t849 * t556 + (t556 * t804 + t558 * t803) * qJD(1)) * qJD(1) / 0.2e1 + (t556 * t803 - t558 * t804) * qJDD(1) / 0.2e1 + t852 * t710 / 0.2e1 + t851 * t709 / 0.2e1 + ((-t707 * t758 + t711) * t556 + t568 * t558 + ((t352 * t819 + t353 * t821 - t405 * t825 - t406 * t827 + t451 * t817 + t695 * t823) * t556 + (t352 * t820 + t353 * t822 + t405 * t826 + t406 * t828 - t451 * t818 + t556 * t759 - t695 * t854 + t560) * t558) * qJD(2) + (t352 * t813 + t353 * t811 - t405 * t814 - t406 * t812 + t451 * t816 - t695 * t853) * qJD(1)) * t671 + ((t556 * t806 + t558 * t805) * qJD(1) + t798) * t670 + ((t556 * t808 + t558 * t807) * qJD(1) + t799) * t669 + ((-t706 * t759 - t711) * t558 + t568 * t556 + ((t350 * t820 + t351 * t822 + t403 * t826 + t404 * t828 - t450 * t818 - t696 * t854) * t558 + (t350 * t819 + t351 * t821 - t403 * t825 - t404 * t827 + t450 * t817 + t558 * t758 + t696 * t823 + t560) * t556) * qJD(2) + (t350 * t813 + t351 * t811 - t403 * t814 - t404 * t812 + t450 * t816 - t696 * t853) * qJD(1)) * t668; (-m(4) + t796) * (g(1) * t451 + g(2) * t450 - g(3) * t756) + m(4) * (t32 * t451 + t33 * t450 + (-t31 * t557 + t68 * t708) * t554) - m(4) * t310 * t68 + (t2 * t451 + t3 * t450 + (-t1 * t557 + t34 * t708) * t554 - t310 * t34) * m(6) + (t450 * t6 + t451 * t5 + (-t4 * t557 + t43 * t708) * t554 - t310 * t43) * m(5); t796 * (g(1) * t352 + g(2) * t350 + g(3) * t446) + (t1 * t446 + t2 * t352 + t3 * t350 + t34 * t840 + t41 * t839 + t42 * t838) * m(6) + (t350 * t6 + t352 * t5 + t4 * t446 + t43 * t840 + t44 * t839 + t45 * t838) * m(5); (-t226 * t41 + t228 * t42 + (-t42 * qJD(1) - g(1) + t2) * t353 + (qJD(1) * t41 + t869) * t351 + (t418 - (t351 * t556 + t353 * t558) * qJD(2)) * t34 + (t1 - (-t41 * t558 - t42 * t556) * qJD(2) - g(3)) * t447) * m(6);];
tau = t7;
