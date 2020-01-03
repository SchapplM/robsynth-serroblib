% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:32:09
% EndTime: 2019-12-31 22:32:49
% DurationCPUTime: 19.30s
% Computational Cost: add. (36873->839), mult. (88298->1188), div. (0->0), fcn. (96858->10), ass. (0->457)
t530 = sin(qJ(5));
t532 = cos(qJ(5));
t529 = sin(pkin(5));
t531 = sin(qJ(2));
t702 = t529 * t531;
t735 = cos(pkin(5));
t789 = sin(qJ(3));
t791 = cos(qJ(3));
t479 = t702 * t791 + t735 * t789;
t563 = t702 * t789 - t735 * t791;
t788 = sin(qJ(4));
t790 = cos(qJ(4));
t549 = t479 * t790 - t563 * t788;
t533 = cos(qJ(2));
t701 = t529 * t533;
t308 = -t530 * t549 - t532 * t701;
t297 = Ifges(6,4) * t308;
t309 = -t530 * t701 + t532 * t549;
t159 = -Ifges(6,2) * t309 + t297;
t777 = Ifges(6,4) * t309;
t160 = Ifges(6,1) * t308 - t777;
t501 = mrSges(6,1) * t530 + mrSges(6,2) * t532;
t776 = Ifges(6,4) * t530;
t504 = Ifges(6,2) * t532 + t776;
t607 = Ifges(6,1) * t532 - t776;
t792 = t532 / 0.4e1;
t796 = t530 / 0.4e1;
t523 = Ifges(6,4) * t532;
t861 = -Ifges(6,2) * t530 + t523;
t375 = t479 * t788 + t790 * t563;
t107 = Ifges(6,1) * t309 + Ifges(6,5) * t375 + t297;
t693 = t532 * t107;
t106 = Ifges(6,2) * t308 + Ifges(6,6) * t375 + t777;
t700 = t530 * t106;
t874 = t693 / 0.4e1 - t700 / 0.4e1;
t506 = Ifges(6,1) * t530 + t523;
t897 = t506 / 0.4e1;
t900 = t309 / 0.4e1;
t901 = t308 / 0.4e1;
t654 = pkin(1) * t735;
t484 = pkin(7) * t701 + t531 * t654;
t458 = pkin(8) * t735 + t484;
t459 = (-pkin(2) * t533 - pkin(8) * t531 - pkin(1)) * t529;
t332 = -t458 * t789 + t791 * t459;
t274 = -t479 * pkin(9) + t332;
t238 = -pkin(3) * t701 + t274;
t333 = t458 * t791 + t459 * t789;
t275 = -pkin(9) * t563 + t333;
t646 = t788 * t275;
t118 = t238 * t790 - t646;
t103 = pkin(4) * t701 - t118;
t902 = t103 / 0.2e1;
t922 = t607 * t900 + t861 * t901 + t159 * t792 + t160 * t796 - t309 * t504 / 0.4e1 + t308 * t897 + t501 * t902 + t874;
t105 = Ifges(6,5) * t309 + Ifges(6,6) * t308 + t375 * Ifges(6,3);
t778 = Ifges(5,4) * t549;
t205 = -Ifges(5,2) * t375 - Ifges(5,6) * t701 + t778;
t362 = Ifges(5,4) * t375;
t206 = Ifges(5,1) * t549 - Ifges(5,5) * t701 - t362;
t492 = t788 * t789 - t790 * t791;
t493 = -t788 * t791 - t789 * t790;
t766 = Ifges(6,3) * t493;
t522 = Ifges(6,5) * t532;
t770 = Ifges(6,6) * t530;
t862 = t522 - t770;
t318 = -t492 * t862 - t766;
t320 = -Ifges(6,6) * t493 - t492 * t861;
t322 = -Ifges(6,5) * t493 - t492 * t607;
t482 = -pkin(7) * t702 + t533 * t654;
t457 = -pkin(2) * t735 - t482;
t389 = pkin(3) * t563 + t457;
t398 = t501 * t492;
t741 = t530 * mrSges(6,3);
t406 = mrSges(6,2) * t493 + t492 * t741;
t740 = t532 * mrSges(6,3);
t408 = -mrSges(6,1) * t493 + t492 * t740;
t413 = -mrSges(5,1) * t493 - mrSges(5,2) * t492;
t415 = -Ifges(5,5) * t492 + Ifges(5,6) * t493;
t487 = Ifges(5,4) * t492;
t416 = Ifges(5,2) * t493 - t487;
t744 = t493 * Ifges(5,4);
t418 = -Ifges(5,1) * t492 + t744;
t648 = t790 * t275;
t119 = t238 * t788 + t648;
t104 = -pkin(10) * t701 + t119;
t142 = t375 * pkin(4) - pkin(10) * t549 + t389;
t49 = -t104 * t530 + t142 * t532;
t50 = t104 * t532 + t142 * t530;
t898 = t493 / 0.4e1;
t899 = -t492 / 0.4e1;
t921 = -t492 * t874 - t415 * t701 / 0.4e1 + t205 * t898 - t493 * t105 / 0.4e1 + t206 * t899 + t389 * t413 / 0.2e1 + t549 * t418 / 0.4e1 + t322 * t900 + t320 * t901 - t103 * t398 / 0.2e1 + t50 * t406 / 0.2e1 + t49 * t408 / 0.2e1 + (-t416 / 0.4e1 + t318 / 0.4e1) * t375;
t502 = Ifges(6,5) * t530 + Ifges(6,6) * t532;
t687 = t532 * t506;
t694 = t530 * t504;
t582 = -t687 / 0.2e1 + t694 / 0.2e1;
t793 = t532 / 0.2e1;
t797 = t530 / 0.2e1;
t804 = -t493 / 0.2e1;
t577 = t320 * t793 + t322 * t797 + t582 * t492 + t502 * t804 + t415;
t500 = -mrSges(6,1) * t532 + mrSges(6,2) * t530;
t677 = t789 * pkin(8);
t578 = -pkin(9) * t789 - t677;
t526 = t791 * pkin(8);
t683 = t791 * pkin(9) + t526;
t851 = t788 * t578 + t683 * t790;
t887 = t851 * t500;
t893 = t851 * mrSges(5,1);
t443 = -t790 * t578 + t683 * t788;
t916 = t443 * mrSges(5,2);
t920 = t577 + t887 - t893 + t916;
t919 = -t887 / 0.2e1 + t893 / 0.2e1 - t916 / 0.2e1;
t864 = t549 * t502;
t881 = Ifges(6,6) * t549 - t375 * t861;
t910 = t532 * t881;
t882 = Ifges(6,5) * t549 - t375 * t607;
t911 = t530 * t882;
t918 = t864 / 0.2e1 + t910 / 0.2e1 + t911 / 0.2e1;
t868 = Ifges(5,6) * t549;
t894 = Ifges(5,5) * t375;
t917 = -t864 / 0.4e1 + t868 / 0.2e1 - (t694 / 0.4e1 - t687 / 0.4e1) * t375 + t894 / 0.2e1 - t910 / 0.4e1 - t911 / 0.4e1;
t915 = t119 * t443;
t914 = t443 * t530;
t913 = t443 * t532;
t912 = t443 * t788;
t712 = t443 * t851;
t909 = -Ifges(5,2) * t549 - t362;
t686 = -t894 - t868;
t521 = -pkin(3) * t791 - pkin(2);
t402 = t492 * pkin(4) + t493 * pkin(10) + t521;
t258 = t402 * t532 - t530 * t851;
t259 = t402 * t530 + t532 * t851;
t319 = t492 * Ifges(6,3) - t493 * t862;
t399 = t501 * t493;
t417 = -t492 * Ifges(5,2) - t744;
t419 = -t493 * Ifges(5,1) - t487;
t747 = t492 * Ifges(6,5);
t323 = -t493 * t607 + t747;
t689 = t532 * t323;
t746 = t492 * Ifges(6,6);
t321 = -t493 * t861 + t746;
t696 = t530 * t321;
t795 = -t532 / 0.2e1;
t905 = (-t689 / 0.2e1 + t696 / 0.2e1 + t318 / 0.2e1 - t419 / 0.2e1 - t416 / 0.2e1) * t492 - t851 * t399 + (t322 * t795 + t320 * t797 - t319 / 0.2e1 - t418 / 0.2e1 + t417 / 0.2e1) * t493 + t258 * t408 + t259 * t406 - t443 * t398 + t521 * t413;
t811 = -t443 / 0.2e1;
t884 = mrSges(5,1) * t549 - mrSges(5,2) * t375;
t904 = (-t419 / 0.4e1 - t689 / 0.4e1 + t696 / 0.4e1 + mrSges(5,3) * t811) * t375 + t521 * t884 / 0.2e1;
t824 = t309 / 0.2e1;
t825 = t308 / 0.2e1;
t891 = t501 * t375;
t903 = -t103 * t891 + t389 * t884 + t882 * t824 + t881 * t825;
t895 = pkin(4) * t851;
t892 = t103 * t851;
t679 = t790 * pkin(3);
t520 = -t679 - pkin(4);
t890 = t520 * t851;
t665 = t770 / 0.2e1;
t590 = -t522 / 0.2e1 + t665;
t889 = t590 * t375;
t750 = t479 * mrSges(4,3);
t441 = -mrSges(4,1) * t701 - t750;
t885 = t750 + t441;
t222 = pkin(4) * t549 + pkin(10) * t375;
t883 = -Ifges(5,1) * t375 - t778;
t866 = Ifges(6,3) * t549;
t880 = -t375 * t862 + t866;
t877 = -t206 / 0.2e1;
t876 = t549 / 0.2e1;
t641 = -t701 / 0.2e1;
t872 = mrSges(6,1) * t549;
t871 = mrSges(6,2) * t549;
t865 = t482 * mrSges(3,2);
t157 = -mrSges(6,1) * t308 + mrSges(6,2) * t309;
t756 = t549 * mrSges(5,3);
t325 = -mrSges(5,1) * t701 - t756;
t863 = -t325 + t157;
t524 = Ifges(4,4) * t791;
t610 = -Ifges(4,2) * t789 + t524;
t527 = t530 ^ 2;
t528 = t532 ^ 2;
t860 = t527 + t528;
t452 = t492 * t701;
t394 = t452 * t530 + t532 * t702;
t451 = t493 * t701;
t272 = mrSges(6,2) * t451 + mrSges(6,3) * t394;
t395 = -t452 * t532 + t530 * t702;
t273 = -mrSges(6,1) * t451 - mrSges(6,3) * t395;
t798 = -t530 / 0.2e1;
t859 = t272 * t793 + t273 * t798;
t191 = -mrSges(6,2) * t375 + mrSges(6,3) * t308;
t192 = mrSges(6,1) * t375 - t309 * mrSges(6,3);
t584 = t191 * t798 + t192 * t795;
t856 = -t791 * mrSges(4,1) + t789 * mrSges(4,2);
t855 = -Ifges(4,5) * t563 - Ifges(4,6) * t479 + t686;
t420 = -t493 * pkin(4) + t492 * pkin(10);
t678 = t789 * pkin(3);
t405 = t678 + t420;
t266 = t405 * t532 + t914;
t267 = t405 * t530 - t913;
t594 = -t266 * t530 + t267 * t532;
t126 = t274 * t790 - t646;
t787 = pkin(3) * t479;
t181 = t222 + t787;
t58 = -t126 * t530 + t181 * t532;
t59 = t126 * t532 + t181 * t530;
t603 = -t58 * t530 + t59 * t532;
t662 = -t747 / 0.2e1;
t854 = t492 * t665 + t532 * t662 - t766 / 0.2e1;
t853 = t417 / 0.4e1 - t319 / 0.4e1;
t852 = mrSges(6,3) * (t308 * t798 + t309 * t793) - t584;
t213 = t375 * t741 - t871;
t215 = t375 * t740 + t872;
t284 = t420 * t532 + t914;
t285 = t420 * t530 - t913;
t706 = t493 * t532;
t409 = mrSges(6,1) * t492 + mrSges(6,3) * t706;
t745 = t493 * mrSges(5,3);
t659 = t745 / 0.2e1;
t550 = t119 * t659 + t921;
t64 = -t118 * t530 + t222 * t532;
t65 = t118 * t532 + t222 * t530;
t763 = t119 * mrSges(5,3);
t664 = -t763 / 0.2e1;
t751 = t851 * mrSges(5,3);
t794 = -t532 / 0.4e1;
t707 = t493 * t530;
t675 = mrSges(6,3) * t707;
t407 = -mrSges(6,2) * t492 + t675;
t814 = t407 / 0.2e1;
t826 = t285 / 0.2e1;
t829 = t258 / 0.2e1;
t835 = t64 / 0.2e1;
t843 = m(6) / 0.2e1;
t849 = t65 * t814 + t409 * t835 - t119 * t399 / 0.2e1 + t284 * t192 / 0.2e1 + t191 * t826 + t259 * t213 / 0.2e1 + t215 * t829 + t550 + (t258 * t64 + t259 * t65 + t284 * t49 + t285 * t50 + t892 + t915) * t843 + (-t883 / 0.4e1 + t882 * t794 + t881 * t796 + t664) * t493 + (-t909 / 0.4e1 + t880 / 0.4e1) * t492 + (-t751 / 0.2e1 - t853) * t549 + (-t325 / 0.2e1 + t157 / 0.2e1) * t851 + t904;
t483 = (pkin(2) * t531 - pkin(8) * t533) * t529;
t387 = -t482 * t789 + t791 * t483;
t651 = t533 * t791;
t294 = (pkin(3) * t531 - pkin(9) * t651) * t529 + t387;
t388 = t791 * t482 + t789 * t483;
t650 = t533 * t789;
t614 = t529 * t650;
t331 = -pkin(9) * t614 + t388;
t154 = t294 * t790 - t331 * t788;
t144 = -pkin(4) * t702 - t154;
t155 = t788 * t294 + t790 * t331;
t775 = Ifges(5,5) * t452;
t848 = t395 * t897 + t394 * t504 / 0.4e1 + t144 * t500 / 0.2e1 - t155 * mrSges(5,2) / 0.2e1 + t154 * mrSges(5,1) / 0.2e1 - t775 / 0.2e1;
t844 = -m(6) / 0.2e1;
t842 = pkin(4) / 0.2e1;
t841 = m(5) * pkin(3);
t840 = m(6) * pkin(3);
t839 = -mrSges(6,1) / 0.2e1;
t838 = -mrSges(6,2) / 0.2e1;
t837 = mrSges(6,2) / 0.2e1;
t836 = Ifges(6,3) / 0.2e1;
t145 = pkin(10) * t702 + t155;
t446 = pkin(3) * t614 + t484;
t231 = -t451 * pkin(4) + t452 * pkin(10) + t446;
t71 = -t145 * t530 + t231 * t532;
t834 = t71 / 0.2e1;
t72 = t145 * t532 + t231 * t530;
t833 = -t72 / 0.2e1;
t831 = -t205 / 0.2e1;
t830 = -t258 / 0.2e1;
t828 = -t259 / 0.2e1;
t827 = -t266 / 0.2e1;
t757 = t375 * mrSges(5,3);
t324 = mrSges(5,2) * t701 - t757;
t823 = -t324 / 0.2e1;
t820 = -t375 / 0.2e1;
t819 = t375 / 0.2e1;
t817 = t394 / 0.2e1;
t816 = t395 / 0.2e1;
t815 = -t407 / 0.2e1;
t813 = -t409 / 0.2e1;
t810 = t443 / 0.2e1;
t809 = t451 / 0.2e1;
t808 = -t451 / 0.2e1;
t807 = -t452 / 0.2e1;
t806 = t479 / 0.2e1;
t805 = t492 / 0.4e1;
t803 = t502 / 0.4e1;
t802 = t862 / 0.4e1;
t676 = t788 * pkin(3);
t519 = t676 + pkin(10);
t801 = -t519 / 0.2e1;
t800 = t520 / 0.2e1;
t156 = mrSges(6,1) * t309 + mrSges(6,2) * t308;
t786 = pkin(4) * t156;
t785 = pkin(4) * t398;
t784 = pkin(4) * t501;
t783 = pkin(10) * t530;
t782 = pkin(10) * t532;
t780 = Ifges(3,4) * t531;
t779 = Ifges(4,4) * t479;
t774 = Ifges(6,5) * t395;
t772 = Ifges(5,6) * t451;
t771 = Ifges(6,6) * t394;
t767 = Ifges(6,3) * t451;
t765 = t118 * mrSges(5,2);
t764 = t119 * mrSges(5,1);
t125 = t274 * t788 + t648;
t762 = t125 * mrSges(5,1);
t761 = t126 * mrSges(5,2);
t758 = t259 * mrSges(6,3);
t748 = t492 * mrSges(5,3);
t695 = t530 * t375;
t214 = mrSges(6,3) * t695 - t871;
t688 = t532 * t375;
t216 = mrSges(6,3) * t688 + t872;
t219 = mrSges(5,1) * t375 + mrSges(5,2) * t549;
t344 = -Ifges(4,2) * t563 - Ifges(4,6) * t701 + t779;
t475 = Ifges(4,4) * t563;
t345 = Ifges(4,1) * t479 - Ifges(4,5) * t701 - t475;
t558 = t563 * mrSges(4,3);
t440 = mrSges(4,2) * t701 - t558;
t552 = t479 * mrSges(4,1) - mrSges(4,2) * t563;
t553 = -Ifges(4,1) * t563 - t779;
t623 = -Ifges(4,2) * t479 - t475;
t628 = -t688 / 0.2e1;
t635 = t695 / 0.2e1;
t5 = t903 + t863 * t125 + t855 * t641 + m(6) * (t103 * t125 + t49 * t58 + t50 * t59) + m(5) * (-t118 * t125 + t119 * t126 + t389 * t787) + t457 * t552 + t375 * t877 + (t118 * t375 - t119 * t549) * mrSges(5,3) - t479 * t344 / 0.2e1 + t332 * t440 + t126 * t324 + t49 * t216 + t50 * t214 + t59 * t191 + t58 * t192 + t880 * t819 + (t883 + t105) * t876 + t553 * t806 + t219 * t787 + t107 * t628 + t549 * t831 + (-t623 / 0.2e1 - t345 / 0.2e1 + t332 * mrSges(4,3)) * t563 + t909 * t820 + t106 * t635 - t885 * t333;
t743 = t5 * qJD(1);
t170 = -t767 + t771 + t774;
t171 = Ifges(6,4) * t395 + Ifges(6,2) * t394 - Ifges(6,6) * t451;
t172 = Ifges(6,1) * t395 + Ifges(6,4) * t394 - Ifges(6,5) * t451;
t225 = -mrSges(6,1) * t394 + mrSges(6,2) * t395;
t287 = -Ifges(5,4) * t452 + Ifges(5,2) * t451 + Ifges(5,6) * t702;
t288 = -Ifges(5,1) * t452 + Ifges(5,4) * t451 + Ifges(5,5) * t702;
t312 = -mrSges(5,1) * t451 - mrSges(5,2) * t452;
t403 = -mrSges(5,2) * t702 + mrSges(5,3) * t451;
t404 = mrSges(5,1) * t702 + mrSges(5,3) * t452;
t425 = (Ifges(4,6) * t531 + t533 * t610) * t529;
t670 = Ifges(4,4) * t789;
t576 = Ifges(4,1) * t791 - t670;
t426 = (Ifges(4,5) * t531 + t533 * t576) * t529;
t574 = mrSges(4,1) * t789 + mrSges(4,2) * t791;
t456 = t574 * t701;
t471 = (-mrSges(4,2) * t531 - mrSges(4,3) * t650) * t529;
t472 = (mrSges(4,1) * t531 - mrSges(4,3) * t651) * t529;
t575 = Ifges(4,5) * t791 - Ifges(4,6) * t789;
t570 = t575 * t533;
t591 = -t614 / 0.2e1;
t640 = t701 / 0.2e1;
t592 = t791 * t640;
t612 = Ifges(3,5) * t701 - Ifges(3,6) * t702;
t642 = t702 / 0.2e1;
t643 = -t702 / 0.2e1;
t6 = t345 * t592 + t344 * t591 + t288 * t876 + t333 * t471 + t332 * t472 + t457 * t456 + t446 * t219 + t388 * t440 + t387 * t441 + t119 * t403 + t118 * t404 + t389 * t312 + t155 * t324 + t154 * t325 + t287 * t820 + t172 * t824 + t171 * t825 + t107 * t816 + t106 * t817 + t170 * t819 + t50 * t272 + t49 * t273 + t103 * t225 + t72 * t191 + t71 * t192 + t144 * t157 + ((Ifges(3,2) * t533 + t780) * t643 + (-pkin(1) * (mrSges(3,1) * t531 + mrSges(3,2) * t533) + t531 * (Ifges(3,1) * t533 - t780) / 0.2e1 - t533 * (Ifges(4,3) * t531 + t570) / 0.2e1) * t529) * t529 + t426 * t806 + t206 * t807 + t105 * t808 + t205 * t809 + t484 * (mrSges(4,1) * t563 + t479 * mrSges(4,2)) - t563 * t425 / 0.2e1 + (0.2e1 * Ifges(3,4) * t701 + (Ifges(3,1) - Ifges(3,2)) * t702) * t640 + (Ifges(4,5) * t479 + Ifges(5,5) * t549 - Ifges(4,6) * t563 - Ifges(5,6) * t375 + (-Ifges(4,3) - Ifges(5,3)) * t701) * t642 + (Ifges(5,3) * t702 + t772 - t775) * t641 + (Ifges(3,6) * t643 + Ifges(3,5) * t640 - t484 * mrSges(3,1) - t865 + t612 / 0.2e1) * t735 + m(4) * (t332 * t387 + t333 * t388 + t457 * t484) + m(5) * (t118 * t154 + t119 * t155 + t389 * t446) + m(6) * (t103 * t144 + t49 * t71 + t50 * t72);
t737 = t6 * qJD(1);
t7 = m(6) * (t49 * t64 + t50 * t65) + t65 * t191 + t50 * t213 + t64 * t192 + t49 * t215 + t118 * t324 + t686 * t641 + (t105 / 0.2e1 + t883 / 0.2e1 + t831 - t763) * t549 + (t880 / 0.2e1 + t877 - t909 / 0.2e1 - t693 / 0.2e1 + t700 / 0.2e1 + t118 * mrSges(5,3)) * t375 + (m(6) * t103 + t863) * t119 + t903;
t736 = t7 * qJD(1);
t732 = t119 * t500;
t158 = Ifges(6,5) * t308 - Ifges(6,6) * t309;
t12 = t158 * t819 + t49 * t191 - t50 * t192 + t103 * t156 + (t160 / 0.2e1 - t106 / 0.2e1 - t50 * mrSges(6,3)) * t309 + (t107 / 0.2e1 + t159 / 0.2e1 - t49 * mrSges(6,3)) * t308;
t731 = t12 * qJD(1);
t730 = t125 * t500;
t713 = t443 * t125;
t705 = t519 * t530;
t704 = t519 * t532;
t703 = t520 * t398;
t681 = t791 / 0.2e1;
t680 = t789 / 0.2e1;
t672 = t783 / 0.2e1;
t671 = -t782 / 0.2e1;
t661 = t747 / 0.2e1;
t660 = -t746 / 0.2e1;
t658 = -t741 / 0.2e1;
t657 = t741 / 0.2e1;
t656 = -t740 / 0.2e1;
t655 = t740 / 0.2e1;
t653 = t530 * t790;
t652 = t532 * t790;
t645 = -t705 / 0.2e1;
t644 = t704 / 0.2e1;
t401 = t506 * t493;
t625 = -t321 / 0.4e1 + t401 / 0.4e1;
t400 = t504 * t493;
t624 = t323 / 0.4e1 + t400 / 0.4e1;
t621 = mrSges(5,3) * t679;
t620 = mrSges(5,3) * t676;
t617 = t679 / 0.2e1;
t615 = t676 / 0.2e1;
t613 = -t653 / 0.2e1;
t611 = mrSges(6,3) * (t528 / 0.2e1 + t527 / 0.2e1);
t602 = -t64 * t530 + t65 * t532;
t601 = -t71 * t530 + t72 * t532;
t599 = pkin(3) * t613;
t598 = t532 * t617;
t414 = mrSges(5,1) * t492 - mrSges(5,2) * t493;
t505 = Ifges(4,2) * t791 + t670;
t507 = Ifges(4,1) * t789 + t524;
t534 = -m(5) * (-t915 + t713 + (t389 * t789 + t479 * t521) * pkin(3) + (-t118 + t126) * t851) / 0.2e1 - t904 + t853 * t549 - t457 * t574 / 0.2e1 - t479 * t576 / 0.4e1 + t529 * t570 / 0.4e1 - t789 * t553 / 0.4e1 + t789 * t344 / 0.4e1 - t414 * t787 / 0.2e1 + pkin(2) * t552 / 0.2e1 - t219 * t678 / 0.2e1 + t751 * t876 - t851 * t157 / 0.2e1 + t851 * t325 / 0.2e1 + t479 * t505 / 0.4e1 + (t126 / 0.2e1 - t118 / 0.2e1) * t748 + t216 * t830 + t192 * t827 + t214 * t828 + t59 * t815 - t267 * t191 / 0.2e1 + (t507 / 0.4e1 + t610 / 0.4e1) * t563 - t881 * t707 / 0.4e1 + t882 * t706 / 0.4e1 - t891 * t811 - t443 * t823 + t880 * t899 + t883 * t898 + t58 * t813 - (t623 + t345) * t791 / 0.4e1 + (t399 / 0.2e1 + t659) * t125 + (t558 / 0.2e1 + t440 / 0.2e1) * t677 + t909 * t805 + t885 * t526 / 0.2e1 + (t258 * t58 + t259 * t59 + t266 * t49 + t267 * t50 + t713 + t892) * t844;
t555 = t72 * t655 + t71 * t658 + t172 * t796 + t171 * t792 - t451 * t803 + t772 / 0.2e1 + Ifges(5,3) * t642 + t848;
t537 = t144 * t520 * t843 + t555 + t225 * t800 + t387 * mrSges(4,1) / 0.2e1 - t388 * mrSges(4,2) / 0.2e1 + Ifges(4,5) * t592 + Ifges(4,6) * t591 + (t154 * t790 + t155 * t788) * t841 / 0.2e1 + Ifges(4,3) * t642 + t403 * t615 + t404 * t617 + (t601 * t843 + t859) * t519;
t1 = t493 * t664 + t534 + t537 - t921;
t17 = t267 * t407 + t266 * t409 + t576 * t680 - pkin(2) * t574 - t789 * t505 / 0.2e1 + m(6) * (t258 * t266 + t259 * t267 + t712) + (t610 + t507) * t681 + t905 + (m(5) * t521 + t414) * t678;
t597 = -t1 * qJD(1) + t17 * qJD(2);
t21 = m(6) * (t258 * t284 + t259 * t285 + t712) + t285 * t407 + t284 * t409 + t905;
t554 = (-pkin(4) * t144 + pkin(10) * t601) * t844 + t225 * t842;
t3 = (-t171 / 0.4e1 + mrSges(6,3) * t833 - pkin(10) * t272 / 0.2e1) * t532 + (-t172 / 0.4e1 + mrSges(6,3) * t834 + pkin(10) * t273 / 0.2e1) * t530 + (t823 - t891 / 0.2e1) * t443 - (-t502 / 0.4e1 + Ifges(5,6) / 0.2e1) * t451 + Ifges(5,3) * t643 + t554 - t848 + t849;
t596 = t3 * qJD(1) + t21 * qJD(2);
t397 = t500 * t493;
t539 = (mrSges(6,3) * t830 + t624) * t308 + (-t758 / 0.2e1 + t625) * t309 + (t375 * t803 + t160 * t794 + t106 * t792 + (t49 * t798 + t50 * t793) * mrSges(6,3) + (t159 + t107) * t796) * t493 + t397 * t902 + t191 * t829 + t192 * t828 + t156 * t810 + t49 * t814 + t158 * t805 + t50 * t813;
t551 = t774 / 0.2e1 + t771 / 0.2e1 - t767 / 0.2e1 + mrSges(6,1) * t834 + mrSges(6,2) * t833;
t10 = t539 - t551;
t27 = -t443 * t397 + t259 * t409 + ((-t323 / 0.2e1 - t400 / 0.2e1 + t662) * t530 + (t401 / 0.2e1 - t321 / 0.2e1 - t758 + t660) * t532) * t493 + (-t407 + t675) * t258;
t595 = t10 * qJD(1) - t27 * qJD(2);
t593 = -t284 * t530 + t285 * t532;
t589 = mrSges(6,1) * t827 + t267 * t837;
t588 = mrSges(6,2) * t826 + t284 * t839;
t587 = pkin(10) * t815 + t625;
t586 = pkin(10) * t813 + t624;
t585 = t492 * t802 + t501 * t810;
t581 = t407 * t801 + t625;
t580 = t409 * t801 + t624;
t573 = t860 * t790;
t538 = (t890 + t593 * t519 + (-t258 * t653 + t259 * t652 + t912) * pkin(3)) * t843 - t703 / 0.2e1 + t284 * t658 + t285 * t655 + t408 * t645 + t406 * t644 - t399 * t615 + t409 * t599 + t407 * t598 - t919;
t541 = (pkin(10) * t594 - t895) * t844 - t785 / 0.2e1 + t266 * t657 + t267 * t656 + t408 * t672 + t406 * t671 + t919;
t22 = t538 + t541;
t546 = (-mrSges(5,1) + t500) * t676 + (mrSges(6,3) * t860 - mrSges(5,2)) * t679;
t223 = (t519 * t573 + t520 * t788) * t840 + t546;
t535 = (t520 * t119 + t602 * t519 + (t103 * t788 - t49 * t653 + t50 * t652) * pkin(3)) * t843 - t765 / 0.2e1 - t764 / 0.2e1 + t732 / 0.2e1 - t891 * t800 + t215 * t645 + t213 * t644 + t64 * t658 + t65 * t655 + t157 * t615 + t192 * t599 + t191 * t598 + (t324 + t757) * t617 - (t325 + t756) * t676 / 0.2e1 - t917;
t540 = (-pkin(4) * t125 + pkin(10) * t603) * t844 - t891 * t842 + t762 / 0.2e1 - t730 / 0.2e1 + t761 / 0.2e1 + t216 * t672 + t214 * t671 + t58 * t657 + t59 * t656 + t917;
t8 = t540 + t535;
t572 = t8 * qJD(1) + t22 * qJD(2) + t223 * qJD(3);
t557 = t375 * t802 + (t657 + t658) * t50 + t922;
t542 = t156 * t800 - t519 * t852 + t557;
t569 = -t866 / 0.2e1 + t58 * t839 + t59 * t837;
t14 = t542 + t569 - t889;
t480 = t520 * t501;
t548 = t607 * t798 + t795 * t861 + t582;
t245 = -t480 + t548;
t547 = t504 * t792 + t607 * t794 + (t506 + t861) * t796;
t543 = t519 * t611 + t547;
t565 = t397 * t800 + t585;
t26 = (t661 + t580) * t532 + (t660 + t581) * t530 + (t836 + t543) * t493 + t565 + t589;
t571 = t14 * qJD(1) + t26 * qJD(2) - t245 * qJD(3);
t568 = t866 / 0.2e1 + mrSges(6,1) * t835 + t65 * t838;
t567 = -pkin(4) * t397 / 0.2e1 + t585;
t15 = t786 / 0.2e1 + (-t862 / 0.4e1 + t590) * t375 + t852 * pkin(10) + t568 - t922;
t545 = t480 / 0.2e1 - t784 / 0.2e1 - t548;
t556 = (mrSges(6,1) * t613 + t652 * t838) * pkin(3);
t179 = t556 - t545;
t268 = t548 + t784;
t544 = pkin(10) * t611 + t547;
t29 = (t661 + t586) * t532 + (t660 + t587) * t530 + (t836 + t544) * t493 + t567 + t588;
t561 = t15 * qJD(1) - t29 * qJD(2) + t179 * qJD(3) + t268 * qJD(4);
t180 = t556 + t545;
t28 = t493 * t544 + t530 * t587 + t532 * t586 + t567 - t588 + t854;
t25 = t493 * t543 + t530 * t581 + t532 * t580 + t565 - t589 + t854;
t18 = t538 - t541 + t577;
t16 = t557 - t786 / 0.2e1 + t889 + t568 + (t308 * t657 + t309 * t656 + t584) * pkin(10);
t13 = Ifges(6,5) * t628 + Ifges(6,6) * t635 + t542 - t569;
t11 = t539 + t551;
t9 = -t540 + t535;
t4 = pkin(10) * t859 + t324 * t811 - t810 * t891 - t554 + t555 + t849;
t2 = t537 - t534 + t550;
t19 = [qJD(2) * t6 + qJD(3) * t5 + qJD(4) * t7 + qJD(5) * t12, t737 + ((-m(4) * pkin(2) - mrSges(3,1) + t856) * t484 + t505 * t591 + t507 * t592 + m(6) * (t258 * t71 + t259 * t72) - t155 * t748 - t172 * t706 / 0.2e1 + t171 * t707 / 0.2e1 - t472 * t677 + m(5) * (t155 * t851 + t446 * t521) + t851 * t403 + (m(4) * pkin(8) + mrSges(4,3)) * (-t789 * t387 + t791 * t388) + t521 * t312 - pkin(2) * t456 + t446 * t414 + t72 * t407 + t71 * t409 - t144 * t399 + t323 * t816 + t321 * t817 + t259 * t272 + t258 * t273 + (Ifges(4,5) * t789 - Ifges(5,5) * t493 + Ifges(4,6) * t791) * t642 + t288 * t804 + t419 * t807 + t319 * t808 + t417 * t809 + t154 * t745 + (-t287 / 0.2e1 + t170 / 0.2e1 - Ifges(5,6) * t642) * t492 + t426 * t680 + t425 * t681 + t471 * t526 + (-m(5) * t154 + m(6) * t144 + t225 - t404) * t443 + t612 - t865) * qJD(2) + t2 * qJD(3) + t4 * qJD(4) + t11 * qJD(5), t743 + t2 * qJD(2) + (-t520 * t891 + t730 - t333 * mrSges(4,1) - t332 * mrSges(4,2) - t762 - t761 + m(6) * (t125 * t520 + t519 * t603) - t549 * t620 + t375 * t621 + t506 * t628 + t504 * t635 + t214 * t704 - t216 * t705 + (-t125 * t790 + t126 * t788) * t841 + t603 * mrSges(6,3) + t855 + t918) * qJD(3) + t9 * qJD(4) + t13 * qJD(5), t736 + t4 * qJD(2) + t9 * qJD(3) + (m(6) * (-pkin(4) * t119 + pkin(10) * t602) + t213 * t782 - t215 * t783 + t732 + pkin(4) * t891 - t764 - t765 + t582 * t375 + t602 * mrSges(6,3) + t686 + t918) * qJD(4) + t16 * qJD(5), t731 + t11 * qJD(2) + t13 * qJD(3) + t16 * qJD(4) + (-mrSges(6,1) * t50 - mrSges(6,2) * t49 + t158) * qJD(5); -qJD(3) * t1 + qJD(4) * t3 + qJD(5) * t10 - t737, qJD(3) * t17 + qJD(4) * t21 - qJD(5) * t27, (m(6) * (t519 * t594 + t890) + t575 - t703 + t492 * t621 + t493 * t620 + t406 * t704 - t408 * t705 + (-t790 * t851 - t912) * t841 + t856 * pkin(8) + t594 * mrSges(6,3) + t920) * qJD(3) + t18 * qJD(4) + t25 * qJD(5) + t597, t18 * qJD(3) + (m(6) * (pkin(10) * t593 - t895) + t406 * t782 - t408 * t783 + t785 + t593 * mrSges(6,3) + t920) * qJD(4) + t28 * qJD(5) + t596, t25 * qJD(3) + t28 * qJD(4) + (-mrSges(6,1) * t259 - mrSges(6,2) * t258 + t493 * t502) * qJD(5) + t595; qJD(2) * t1 + qJD(4) * t8 + qJD(5) * t14 - t743, qJD(4) * t22 + qJD(5) * t26 - t597, qJD(4) * t223 - qJD(5) * t245, ((-pkin(4) * t788 + pkin(10) * t573) * t840 + t546) * qJD(4) + t180 * qJD(5) + t572, t180 * qJD(4) + (t500 * t519 + t862) * qJD(5) + t571; -qJD(2) * t3 - qJD(3) * t8 - qJD(5) * t15 - t736, -qJD(3) * t22 + qJD(5) * t29 - t596, -qJD(5) * t179 - t572, -t268 * qJD(5), (pkin(10) * t500 + t862) * qJD(5) - t561; -qJD(2) * t10 - qJD(3) * t14 + qJD(4) * t15 - t731, -qJD(3) * t26 - qJD(4) * t29 - t595, qJD(4) * t179 - t571, t561, 0;];
Cq = t19;
