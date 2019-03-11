% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:21:59
% EndTime: 2019-03-09 13:22:46
% DurationCPUTime: 27.71s
% Computational Cost: add. (69210->838), mult. (134396->1123), div. (0->0), fcn. (157455->10), ass. (0->443)
t540 = sin(qJ(4));
t705 = sin(pkin(11));
t624 = t705 * pkin(2);
t529 = t624 + pkin(8);
t761 = pkin(9) + t529;
t499 = t761 * t540;
t543 = cos(qJ(4));
t500 = t761 * t543;
t539 = sin(qJ(5));
t542 = cos(qJ(5));
t424 = -t499 * t539 + t542 * t500;
t861 = -t539 * t540 + t542 * t543;
t365 = pkin(10) * t861 + t424;
t538 = sin(qJ(6));
t541 = cos(qJ(6));
t590 = t539 * t543 + t542 * t540;
t864 = -t542 * t499 - t539 * t500;
t888 = -t590 * pkin(10) + t864;
t925 = -t365 * t538 + t541 * t888;
t830 = -t925 / 0.2e1;
t886 = qJD(4) + qJD(5);
t448 = t538 * t861 + t541 * t590;
t611 = -t538 * t590 + t541 * t861;
t649 = Ifges(7,5) * t611 - Ifges(7,6) * t448;
t931 = t925 * mrSges(7,2);
t233 = t365 * t541 + t538 * t888;
t935 = t233 * mrSges(7,1);
t939 = -t935 / 0.2e1 - t931 / 0.2e1;
t944 = 0.2e1 * t939 + t649;
t946 = t944 * qJD(6);
t752 = Ifges(7,4) * t448;
t332 = Ifges(7,2) * t611 + t752;
t440 = Ifges(7,4) * t611;
t335 = Ifges(7,1) * t448 + t440;
t450 = mrSges(6,1) * t590 + mrSges(6,2) * t861;
t504 = Ifges(6,4) * t861;
t452 = -Ifges(6,2) * t590 + t504;
t706 = cos(pkin(11));
t625 = t706 * pkin(2);
t530 = -t625 - pkin(3);
t515 = -t543 * pkin(4) + t530;
t455 = Ifges(6,1) * t590 + t504;
t796 = t455 / 0.2e1;
t800 = t448 / 0.2e1;
t875 = t611 / 0.2e1;
t878 = -Ifges(7,2) * t448 + t440;
t896 = -t448 / 0.2e1;
t904 = Ifges(7,1) * t611 - t752;
t468 = -pkin(5) * t861 + t515;
t869 = t611 * mrSges(7,2);
t893 = t448 * mrSges(7,1);
t902 = t893 + t869;
t922 = t468 * t902;
t945 = t922 + t515 * t450 + (t796 + t452 / 0.2e1) * t861 + t332 * t896 + t904 * t800 + (t878 + t335) * t875;
t329 = -mrSges(7,1) * t611 + mrSges(7,2) * t448;
t943 = m(7) * t468 + t329;
t778 = sin(qJ(2));
t779 = cos(qJ(2));
t506 = t705 * t778 - t706 * t779;
t508 = -t705 * t779 - t706 * t778;
t395 = t861 * t508;
t397 = t590 * t508;
t612 = t395 * t538 + t541 * t397;
t290 = t395 * t541 - t397 * t538;
t753 = Ifges(7,4) * t290;
t156 = Ifges(7,2) * t612 + t506 * Ifges(7,6) - t753;
t285 = Ifges(7,4) * t612;
t158 = -Ifges(7,1) * t290 + Ifges(7,5) * t506 + t285;
t640 = t778 * pkin(7);
t516 = -qJ(3) * t778 - t640;
t642 = t779 * pkin(7);
t518 = qJ(3) * t779 + t642;
t461 = -t706 * t516 + t518 * t705;
t678 = t508 * t540;
t375 = -pkin(4) * t678 + t461;
t297 = -t397 * pkin(5) + t375;
t246 = -mrSges(7,2) * t506 + mrSges(7,3) * t612;
t828 = t246 / 0.2e1;
t879 = Ifges(7,2) * t290 + t285;
t880 = -mrSges(7,1) * t290 + mrSges(7,2) * t612;
t903 = Ifges(7,1) * t612 + t753;
t942 = t925 * t828 + t880 * t468 / 0.2e1 + (t158 / 0.4e1 + t879 / 0.4e1) * t611 + t902 * t297 / 0.2e1 + (t903 / 0.4e1 - t156 / 0.4e1) * t448;
t605 = Ifges(6,5) * t861 - Ifges(6,6) * t590 + t649;
t920 = -t424 * mrSges(6,1) - t864 * mrSges(6,2) + t605;
t926 = -t931 - t935;
t941 = t920 + t926;
t829 = -t233 / 0.2e1;
t932 = (t830 + t925 / 0.2e1) * mrSges(7,2) + (t829 + t233 / 0.2e1) * mrSges(7,1);
t940 = qJD(6) * t932;
t626 = -t893 / 0.2e1;
t804 = -t611 / 0.2e1;
t572 = t626 + t893 / 0.2e1 + (t804 + t875) * mrSges(7,2);
t938 = qJD(2) * t932 + qJD(3) * t572;
t934 = -t233 * t541 + t538 * t925;
t391 = Ifges(6,4) * t397;
t255 = -Ifges(6,1) * t395 + Ifges(6,5) * t506 + t391;
t301 = -mrSges(6,1) * t395 + mrSges(6,2) * t397;
t304 = Ifges(6,2) * t395 + t391;
t789 = t506 / 0.4e1;
t248 = mrSges(7,1) * t506 + t290 * mrSges(7,3);
t827 = t248 / 0.2e1;
t901 = t878 / 0.4e1 + t335 / 0.4e1;
t918 = t332 / 0.4e1 - t904 / 0.4e1;
t933 = -t233 * t827 + (t304 / 0.4e1 + t255 / 0.4e1) * t861 + t375 * t450 / 0.2e1 + t515 * t301 / 0.2e1 + t605 * t789 + t901 * t612 + t918 * t290 + t942;
t836 = m(7) * pkin(5);
t915 = -t836 / 0.2e1;
t787 = t861 / 0.2e1;
t534 = -pkin(2) * t779 - pkin(1);
t433 = t506 * pkin(3) + t508 * pkin(8) + t534;
t863 = t705 * t516 + t706 * t518;
t325 = t543 * t433 - t540 * t863;
t677 = t508 * t543;
t286 = pkin(9) * t677 + t325;
t249 = t506 * pkin(4) + t286;
t326 = t540 * t433 + t543 * t863;
t287 = pkin(9) * t678 + t326;
t268 = t539 * t287;
t141 = t542 * t249 - t268;
t163 = t542 * t286 - t268;
t921 = t141 - t163;
t396 = t590 * t506;
t398 = t861 * t506;
t291 = t396 * t541 + t398 * t538;
t294 = t396 * t538 - t398 * t541;
t730 = t294 * mrSges(7,2);
t731 = t291 * mrSges(7,1);
t655 = t731 / 0.2e1 - t730 / 0.2e1;
t726 = t398 * mrSges(6,2);
t728 = t396 * mrSges(6,1);
t607 = t655 + t726 / 0.2e1 + t728 / 0.2e1;
t919 = (t291 * t541 + t294 * t538) * t915 - t607;
t806 = t424 / 0.2e1;
t788 = -t508 / 0.2e1;
t916 = -t612 / 0.2e1;
t641 = t778 * pkin(2);
t914 = m(4) * t641;
t911 = t297 * t880;
t270 = t542 * t287;
t142 = t249 * t539 + t270;
t162 = -t286 * t539 - t270;
t906 = t142 + t162;
t821 = t294 / 0.2e1;
t823 = t291 / 0.2e1;
t587 = Ifges(7,5) * t821 + Ifges(7,6) * t823;
t434 = -t508 * pkin(3) + t506 * pkin(8) + t641;
t336 = t543 * t434 + t461 * t540;
t679 = t506 * t543;
t250 = -t508 * pkin(4) + pkin(9) * t679 + t336;
t337 = t540 * t434 - t461 * t543;
t680 = t506 * t540;
t288 = pkin(9) * t680 + t337;
t150 = t542 * t250 - t288 * t539;
t100 = -pkin(5) * t508 + pkin(10) * t398 + t150;
t151 = t539 * t250 + t542 * t288;
t104 = pkin(10) * t396 + t151;
t75 = t100 * t541 - t104 * t538;
t76 = t100 * t538 + t104 * t541;
t862 = t76 * mrSges(7,2) / 0.2e1 - t75 * mrSges(7,1) / 0.2e1;
t600 = Ifges(7,3) * t788 + t587 - t862;
t870 = Ifges(7,5) * t612;
t895 = Ifges(7,6) * t290;
t654 = t870 + t895;
t617 = t895 / 0.2e1 + t870 / 0.2e1;
t900 = -t508 * mrSges(4,1) - t506 * mrSges(4,2);
t727 = t397 * mrSges(6,3);
t344 = -mrSges(6,2) * t506 + t727;
t729 = t395 * mrSges(6,3);
t346 = mrSges(6,1) * t506 + t729;
t786 = t590 / 0.2e1;
t899 = (t290 * t800 + t611 * t916) * mrSges(7,3) - t346 * t786 + t344 * t787 - t246 * t804 - t248 * t800;
t897 = t290 / 0.2e1;
t814 = -t346 / 0.2e1;
t711 = t590 * mrSges(6,3);
t891 = t506 * (t540 ^ 2 + t543 ^ 2);
t773 = pkin(10) * t397;
t103 = t142 + t773;
t703 = t103 * t538;
t394 = t395 * pkin(10);
t102 = t141 + t394;
t99 = pkin(5) * t506 + t102;
t71 = t541 * t99 - t703;
t889 = -t71 * mrSges(7,3) + t158 / 0.2e1;
t755 = Ifges(6,4) * t395;
t253 = Ifges(6,2) * t397 + t506 * Ifges(6,6) - t755;
t305 = Ifges(6,1) * t397 + t755;
t601 = (t253 / 0.4e1 - t305 / 0.4e1) * t590;
t754 = Ifges(6,4) * t590;
t453 = Ifges(6,2) * t861 + t754;
t454 = Ifges(6,1) * t861 - t754;
t614 = t453 / 0.4e1 - t454 / 0.4e1;
t887 = t614 * t395 - t601;
t885 = t915 * t934 - t939;
t111 = t162 - t773;
t112 = t394 + t163;
t86 = t111 * t538 + t112 * t541;
t762 = t86 * mrSges(7,2);
t85 = t111 * t541 - t112 * t538;
t763 = t85 * mrSges(7,1);
t759 = t763 / 0.2e1 - t762 / 0.2e1;
t884 = (t538 * t86 + t541 * t85) * t915 - t759;
t881 = t151 * mrSges(6,2) / 0.2e1 - t150 * mrSges(6,1) / 0.2e1;
t535 = Ifges(5,5) * t543;
t748 = Ifges(5,6) * t540;
t868 = Ifges(4,4) - t535 / 0.2e1 + t748 / 0.2e1;
t810 = -t398 / 0.2e1;
t811 = t396 / 0.2e1;
t866 = Ifges(6,5) * t810 + Ifges(6,6) * t811;
t536 = Ifges(5,4) * t543;
t521 = Ifges(5,1) * t540 + t536;
t606 = Ifges(6,5) * t397 + Ifges(6,6) * t395 + t654;
t598 = Ifges(5,2) * t540 - t536;
t702 = t103 * t541;
t72 = t538 * t99 + t702;
t834 = -t71 / 0.2e1;
t564 = t233 * t897 + t611 * t834 + t612 * t830 + t72 * t896;
t615 = t452 / 0.4e1 + t455 / 0.4e1;
t809 = -t864 / 0.2e1;
t81 = -t102 * t538 - t702;
t815 = t344 / 0.2e1;
t82 = t102 * t541 - t703;
t841 = m(7) / 0.2e1;
t855 = ((t72 + t81) * t925 + (-t71 + t82) * t233) * t841 + (mrSges(6,3) * t809 + t615) * t397 + (t81 * t896 + t82 * t875 + t564) * mrSges(7,3) + t864 * t815 + t933;
t854 = qJD(6) * t572;
t575 = -t869 + 0.2e1 * t626;
t853 = t575 * qJD(6);
t768 = t72 * mrSges(7,1);
t770 = t71 * mrSges(7,2);
t852 = -t768 / 0.2e1 - t770 / 0.2e1 + t617;
t168 = -mrSges(7,1) * t612 - mrSges(7,2) * t290;
t813 = -t395 / 0.2e1;
t851 = (t297 * t590 - t395 * t468) * t841 + t168 * t786 + t329 * t813;
t565 = t587 + t866;
t850 = -t565 + t862 + t881;
t849 = Ifges(6,3) * t788 + t600 + t866 - t881;
t845 = m(5) / 0.2e1;
t844 = -m(6) / 0.2e1;
t843 = m(6) / 0.2e1;
t842 = -m(7) / 0.2e1;
t840 = pkin(4) / 0.2e1;
t838 = m(4) * pkin(2);
t837 = m(6) * pkin(4);
t833 = -t72 / 0.2e1;
t822 = t612 / 0.2e1;
t820 = -t290 / 0.2e1;
t817 = t335 / 0.2e1;
t812 = t395 / 0.2e1;
t808 = t864 / 0.2e1;
t807 = -t424 / 0.2e1;
t799 = t453 / 0.2e1;
t533 = pkin(4) * t542 + pkin(5);
t670 = t538 * t539;
t489 = -pkin(4) * t670 + t533 * t541;
t794 = -t489 / 0.2e1;
t666 = t539 * t541;
t490 = pkin(4) * t666 + t533 * t538;
t793 = -t490 / 0.2e1;
t497 = (-t538 * t542 - t666) * pkin(4);
t792 = -t497 / 0.2e1;
t498 = (t541 * t542 - t670) * pkin(4);
t791 = t498 / 0.2e1;
t790 = t506 / 0.2e1;
t785 = -t529 / 0.2e1;
t784 = -t538 / 0.2e1;
t783 = t540 / 0.2e1;
t782 = -t541 / 0.2e1;
t781 = -t543 / 0.2e1;
t780 = t543 / 0.2e1;
t777 = (-t448 * t541 + t538 * t611) * t836;
t776 = pkin(4) * t540;
t775 = pkin(5) * t590;
t774 = pkin(5) * t538;
t771 = t395 * pkin(5);
t765 = t81 * mrSges(7,1);
t764 = t82 * mrSges(7,2);
t760 = t765 / 0.2e1 - t764 / 0.2e1;
t758 = mrSges(4,3) * t508;
t756 = Ifges(5,4) * t540;
t751 = Ifges(5,5) * t506;
t749 = Ifges(5,6) * t506;
t746 = pkin(4) * qJD(4);
t745 = pkin(5) * qJD(5);
t744 = t141 * mrSges(6,2);
t743 = t142 * mrSges(6,1);
t740 = t162 * mrSges(6,1);
t739 = t163 * mrSges(6,2);
t720 = t611 * mrSges(7,3);
t717 = t448 * mrSges(7,3);
t716 = t497 * mrSges(7,1);
t715 = t498 * mrSges(7,2);
t155 = Ifges(7,4) * t294 + Ifges(7,2) * t291 - Ifges(7,6) * t508;
t157 = Ifges(7,1) * t294 + Ifges(7,4) * t291 - Ifges(7,5) * t508;
t167 = t730 - t731;
t245 = mrSges(7,2) * t508 + mrSges(7,3) * t291;
t247 = -mrSges(7,1) * t508 - mrSges(7,3) * t294;
t252 = -Ifges(6,4) * t398 + Ifges(6,2) * t396 - Ifges(6,6) * t508;
t254 = -Ifges(6,1) * t398 + Ifges(6,4) * t396 - Ifges(6,5) * t508;
t374 = -pkin(4) * t680 + t863;
t296 = -t396 * pkin(5) + t374;
t302 = -t726 - t728;
t303 = -mrSges(6,1) * t397 - mrSges(6,2) * t395;
t343 = mrSges(6,2) * t508 + mrSges(6,3) * t396;
t345 = -mrSges(6,1) * t508 + mrSges(6,3) * t398;
t357 = -Ifges(5,6) * t508 + t506 * t598;
t358 = t508 * t598 + t749;
t522 = Ifges(5,1) * t543 - t756;
t359 = -Ifges(5,5) * t508 - t506 * t522;
t360 = -t508 * t522 + t751;
t709 = t543 * mrSges(5,2);
t710 = t540 * mrSges(5,1);
t517 = t709 + t710;
t417 = t517 * t506;
t418 = t517 * t508;
t441 = t508 * mrSges(5,2) + mrSges(5,3) * t680;
t442 = -t506 * mrSges(5,2) + mrSges(5,3) * t678;
t443 = -t508 * mrSges(5,1) + mrSges(5,3) * t679;
t444 = t506 * mrSges(5,1) + mrSges(5,3) * t677;
t5 = (Ifges(3,1) - Ifges(3,2)) * t779 * t778 + t254 * t813 + t157 * t820 + t255 * t810 + t253 * t811 + m(6) * (t141 * t150 + t142 * t151 + t374 * t375) + m(7) * (t296 * t297 + t71 * t75 + t72 * t76) + m(5) * (t325 * t336 + t326 * t337 + t461 * t863) - pkin(1) * (mrSges(3,1) * t778 + mrSges(3,2) * t779) + (-t778 ^ 2 + t779 ^ 2) * Ifges(3,4) - t461 * t417 + t326 * t441 + t337 * t442 + t325 * t443 + t336 * t444 + t397 * t252 / 0.2e1 + t374 * t303 + t375 * t302 + t151 * t344 + t141 * t345 + t150 * t346 + t142 * t343 + t296 * t168 + t297 * t167 + t71 * t247 + t75 * t248 + t72 * t245 + t76 * t246 + (mrSges(4,1) * t641 + t358 * t783 + t360 * t781 + t868 * t506 + t565) * t506 - t863 * t418 + t158 * t821 + t155 * t822 + t156 * t823 + (t900 + t914) * t534 + (t359 * t781 + t357 * t783 + Ifges(6,5) * t812 - Ifges(6,6) * t397 / 0.2e1 + Ifges(7,5) * t897 + Ifges(7,6) * t916 - mrSges(4,2) * t641 - t868 * t508 + (-Ifges(5,3) - Ifges(4,2) + Ifges(4,1) - Ifges(6,3) - Ifges(7,3)) * t506) * t508;
t714 = t5 * qJD(1);
t713 = t506 * mrSges(4,3);
t712 = t861 * mrSges(6,3);
t342 = -pkin(4) * t677 - t771;
t599 = -mrSges(5,1) * t543 + t540 * mrSges(5,2);
t416 = t508 * t599;
t519 = Ifges(5,2) * t543 + t756;
t419 = t508 * t519;
t420 = t508 * t521;
t608 = t156 / 0.2e1 + t72 * mrSges(7,3);
t550 = t911 + t375 * t301 + t608 * t290 + (t304 / 0.2e1 + t255 / 0.2e1) * t397 - t141 * t727 + t879 * t822 + t903 * t820 + t606 * t790 + t889 * t612;
t629 = t751 / 0.2e1;
t6 = t550 + t461 * t416 + t325 * t442 - t326 * t444 + t163 * t344 + t162 * t346 + t342 * t168 + t85 * t248 + t86 * t246 + ((t419 / 0.2e1 + t360 / 0.2e1 + t629 - t325 * mrSges(5,3)) * t540 + (-t420 / 0.2e1 + t358 / 0.2e1 + t749 / 0.2e1 + t326 * mrSges(5,3) + (-m(6) * t375 - t303) * pkin(4)) * t543) * t508 + m(6) * (t141 * t162 + t142 * t163) + (t142 * mrSges(6,3) + t253 / 0.2e1 - t305 / 0.2e1) * t395 + m(7) * (t297 * t342 + t71 * t85 + t72 * t86);
t708 = t6 * qJD(1);
t7 = t550 + (-t346 + t729) * t142 - t168 * t771 + t253 * t812 + t305 * t813 + t141 * t344 + t81 * t248 + t82 * t246 + m(7) * (-t297 * t771 + t71 * t81 + t72 * t82);
t707 = t7 * qJD(1);
t451 = -mrSges(6,1) * t861 + mrSges(6,2) * t590;
t546 = (-t530 * t508 - t529 * t891) * t845 + (t396 * t864 - t398 * t424 - t508 * t515) * t843 + (t233 * t294 + t291 * t925 - t468 * t508) * t841 - t416 / 0.2e1 + (-t506 * t705 + t508 * t706) * t838 / 0.2e1 - t291 * t717 / 0.2e1 + t720 * t821 - t711 * t811 + t712 * t810 + (t329 + t451) * t788 - mrSges(5,3) * t891 / 0.2e1;
t548 = (t336 * t543 + t540 * t337) * t845 + (t150 * t861 + t151 * t590) * t843 + (t448 * t76 + t611 * t75) * t841 + t247 * t875 + t245 * t800 + t345 * t787 + t343 * t786 + t441 * t783 + t443 * t780 + t914 / 0.2e1;
t21 = -t546 + t548 + t900;
t704 = qJD(1) * t21;
t12 = t911 + t654 * t790 + t71 * t246 - t72 * t248 - (t903 / 0.2e1 - t608) * t290 + (t879 / 0.2e1 + t889) * t612;
t701 = t12 * qJD(1);
t658 = t543 * t442;
t665 = t540 * t444;
t686 = t461 * t508;
t23 = t294 * t246 + t291 * t248 - t398 * t344 + t396 * t346 + (-t658 + t665 + t713) * t506 + (-t168 - t303 + t418 + t758) * t508 + m(7) * (t291 * t71 + t294 * t72 - t297 * t508) + m(6) * (t141 * t396 - t142 * t398 - t375 * t508) + m(5) * (-t686 + (t325 * t540 - t326 * t543) * t506) + m(4) * (-t506 * t863 - t686);
t700 = t23 * qJD(1);
t561 = (-t290 * t896 + t612 * t804) * mrSges(7,3) + t246 * t875 + t248 * t896;
t25 = t561 - t655;
t698 = t25 * qJD(1);
t691 = t395 * t590;
t690 = t397 * t861;
t685 = t489 * t612;
t684 = t489 * t611;
t683 = t490 * t290;
t682 = t490 * t448;
t674 = t529 * t540;
t673 = t529 * t543;
t672 = t538 * t290;
t671 = t538 * t448;
t669 = t539 * t343;
t668 = t539 * t395;
t664 = t541 * t247;
t663 = t541 * t612;
t662 = t541 * t611;
t661 = t542 * t345;
t660 = t542 * t397;
t223 = -t490 * mrSges(7,1) - t489 * mrSges(7,2);
t646 = qJD(6) * t223;
t644 = t837 / 0.2e1;
t639 = t777 / 0.2e1;
t638 = t774 / 0.2e1;
t637 = pkin(5) * t782;
t632 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t630 = mrSges(5,3) * t529 / 0.2e1;
t622 = t680 / 0.2e1;
t621 = -t679 / 0.2e1;
t619 = -t141 / 0.2e1 + t163 / 0.2e1;
t618 = t142 / 0.2e1 + t162 / 0.2e1;
t610 = mrSges(7,3) * t638;
t609 = mrSges(7,3) * t637;
t603 = -t774 / 0.2e1 + t793;
t602 = t637 + t794;
t597 = -t448 * t71 + t611 * t72;
t38 = t922 + (t904 / 0.2e1 - t332 / 0.2e1) * t448 + (t817 + t878 / 0.2e1) * t611;
t549 = -(mrSges(7,3) * t829 - t918) * t290 + (mrSges(7,3) * t830 + t901) * t612 + t248 * t829 + t649 * t789 + t942;
t9 = t549 - t600;
t596 = t9 * qJD(1) + t38 * qJD(2);
t552 = (-t590 * t921 + t906 * t861) * t844 + (t448 * t86 + t611 * t85 + t597) * t842 + t665 / 0.2e1 - t658 / 0.2e1;
t553 = (t710 / 0.2e1 + t709 / 0.2e1) * t506 + (t291 * t489 + t294 * t490) * t841 + (t396 * t542 - t398 * t539) * t644 + t607;
t15 = (-t691 / 0.2e1 + t690 / 0.2e1) * mrSges(6,3) + t552 + t553 - t899;
t595 = t15 * qJD(1);
t574 = -(-t691 + t690) * mrSges(6,3) / 0.2e1 + t899;
t557 = (t448 * t82 + t611 * t81 + t597) * t841 + t574;
t18 = t557 + t919;
t594 = t18 * qJD(1);
t591 = -t448 * t489 + t490 * t611;
t585 = m(7) * (t538 * t76 + t541 * t75);
t584 = t683 / 0.2e1 - t685 / 0.2e1;
t583 = t444 * t785 + t419 / 0.4e1 + t360 / 0.4e1;
t580 = -t902 - t450;
t474 = t775 + t776;
t544 = (t395 * t806 + t397 * t809 - t590 * t618 + t619 * t861) * mrSges(6,3) + (t375 * t776 + t906 * t864) * t843 + t615 * t397 + t461 * t517 / 0.2e1 + t530 * t416 / 0.2e1 + t474 * t168 / 0.2e1 + t344 * t808 + t342 * t329 / 0.2e1 + (t297 * t474 + t342 * t468 + (t72 + t85) * t925 + (-t71 + t86) * t233) * t841 + t535 * t789 + (t85 * t896 + t86 * t875 + t564) * mrSges(7,3) + t887 + (-t843 * t921 + t814) * t424 + t933;
t554 = (t489 * t75 + t490 * t76) * t842 - t336 * mrSges(5,1) / 0.2e1 + t337 * mrSges(5,2) / 0.2e1 + t247 * t794 + t245 * t793 - (t150 * t542 + t151 * t539) * t837 / 0.2e1;
t555 = (t521 / 0.4e1 - t598 / 0.4e1 + t540 * t630) * t540 + (-t522 / 0.4e1 + t519 / 0.4e1 + t543 * t630 + (t515 * t844 - t451 / 0.2e1) * pkin(4)) * t543;
t573 = -t358 / 0.4e1 + t420 / 0.4e1 + t303 * t840 + t442 * t785;
t2 = t554 + t544 + (t629 + t583) * t543 + (Ifges(5,3) / 0.2e1 + t555 + t632) * t508 + (-0.3e1 / 0.4e1 * t749 + t573) * t540 + (-t661 / 0.2e1 - t669 / 0.2e1) * pkin(4) + t850;
t29 = -(t799 - t454 / 0.2e1) * t590 + m(6) * t515 * t776 + (-t598 / 0.2e1 + t521 / 0.2e1) * t543 + t530 * t517 + (pkin(4) * t451 - t519 / 0.2e1 + t522 / 0.2e1) * t540 + t943 * t474 + t945;
t577 = t2 * qJD(1) + t29 * qJD(2);
t3 = t632 * t508 + t346 * t807 + (-t664 / 0.2e1 + t245 * t784 - t585 / 0.2e1 + t851) * pkin(5) - t601 + (mrSges(6,3) * t806 + t614) * t395 + t850 + t855;
t31 = t454 * t786 - t590 * t799 + t775 * t943 + t945;
t576 = t3 * qJD(1) + t31 * qJD(2);
t551 = (t542 * t815 + t539 * t814 + (t668 / 0.2e1 - t660 / 0.2e1) * mrSges(6,3)) * pkin(4) + (t489 * t81 + t490 * t82 + t497 * t71 + t498 * t72) * t841 + t497 * t827 + t246 * t791 + t760;
t11 = t619 * mrSges(6,2) - t618 * mrSges(6,1) + ((-t672 / 0.2e1 + t663 / 0.2e1) * pkin(5) + t584) * mrSges(7,3) + t551 + t884;
t566 = t716 + (-mrSges(6,1) * t539 - mrSges(6,2) * t542) * pkin(4) - t715;
t215 = m(7) * (t489 * t497 + t490 * t498) + t566;
t560 = ((t490 + t497) * t925 + (-t489 + t498) * t233) * t841 + t939;
t563 = t611 * t791 + t448 * t792 - t682 / 0.2e1 - t684 / 0.2e1;
t37 = (t808 + t809) * mrSges(6,2) + (t806 + t807) * mrSges(6,1) + ((t671 / 0.2e1 + t662 / 0.2e1) * pkin(5) + t563) * mrSges(7,3) + t560 + t885;
t567 = m(7) * (t448 * t498 + t497 * t611 + t591);
t94 = t639 - t567 / 0.2e1;
t570 = t11 * qJD(1) + t37 * qJD(2) - t94 * qJD(3) + t215 * qJD(4);
t559 = (-t290 * t793 + t612 * t794) * mrSges(7,3) + t489 * t828 + t248 * t793 + t617;
t14 = (t834 + t86 / 0.2e1) * mrSges(7,2) + (t833 - t85 / 0.2e1) * mrSges(7,1) + t559 - t617;
t569 = t14 * qJD(1) + t223 * qJD(4) + t938;
t196 = (t791 + t602) * mrSges(7,2) + (t792 + t603) * mrSges(7,1);
t558 = (t541 * t828 + t248 * t784 + (-t290 * t784 + t612 * t782) * mrSges(7,3)) * pkin(5) + t617;
t20 = (t834 + t82 / 0.2e1) * mrSges(7,2) + (t833 - t81 / 0.2e1) * mrSges(7,1) + t558 - t617;
t514 = (mrSges(7,1) * t538 + mrSges(7,2) * t541) * pkin(5);
t562 = -qJD(1) * t20 - qJD(4) * t196 + qJD(5) * t514 - t938;
t509 = t514 * qJD(6);
t197 = -t715 / 0.2e1 + t716 / 0.2e1 + t602 * mrSges(7,2) + t603 * mrSges(7,1);
t89 = t567 / 0.2e1 + t639 + t580;
t36 = t563 * mrSges(7,3) - t448 * t610 + t611 * t609 + t560 - t885 + t920;
t24 = t561 + t655;
t22 = t546 + t548;
t19 = t558 + t760 + t852;
t17 = t557 - t919;
t16 = -t552 + t553 + t574;
t13 = t559 + t759 + t852;
t10 = t551 + t584 * mrSges(7,3) + t740 / 0.2e1 - t739 / 0.2e1 - t743 / 0.2e1 - t744 / 0.2e1 + t612 * t609 + t290 * t610 + t606 - t884;
t8 = t549 + t600;
t4 = (t814 + t729 / 0.2e1) * t424 + t245 * t638 + (t664 + t585) * pkin(5) / 0.2e1 + t851 * pkin(5) + t849 + t855 + t887;
t1 = -t554 + t555 * t508 + t544 + Ifges(5,5) * t621 + Ifges(5,6) * t622 + Ifges(5,3) * t788 + t583 * t543 + (-t749 / 0.4e1 + t573) * t540 + (t661 + t669) * t840 + t849;
t26 = [qJD(2) * t5 + qJD(3) * t23 + qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t12, t22 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + t8 * qJD(6) + t714 + (t521 * t621 + t519 * t622 + m(7) * (t233 * t76 + t296 * t468 + t75 * t925) + t925 * t247 + t294 * t817 + t157 * t800 - t398 * t796 + t396 * t799 + t357 * t780 + t359 * t783 + t254 * t786 + t252 * t787 + t624 * t758 + t76 * t720 + (Ifges(5,5) * t540 + Ifges(6,5) * t590 + Ifges(7,5) * t448 + Ifges(5,6) * t543 + Ifges(6,6) * t861 + Ifges(7,6) * t611) * t788 + (m(5) * t529 + mrSges(5,3)) * (-t336 * t540 + t337 * t543) - Ifges(3,6) * t778 + Ifges(3,5) * t779 - t75 * t717 - t443 * t674 - t150 * t711 + t151 * t712 + t625 * t713 + t441 * t673 + mrSges(3,2) * t640 - mrSges(3,1) * t642 + t515 * t302 - t530 * t417 - Ifges(4,5) * t506 + Ifges(4,6) * t508 + t468 * t167 + t374 * t451 + t424 * t343 + t296 * t329 + t233 * t245 + t155 * t875 + (m(5) * t530 - t706 * t838 - mrSges(4,1) + t599) * t863 + (-t705 * t838 + mrSges(4,2)) * t461 + m(6) * (t150 * t864 + t151 * t424 + t374 * t515) + t864 * t345 + t332 * t823) * qJD(2), t700 + t22 * qJD(2) + 0.2e1 * ((t291 * t611 + t294 * t448) * t841 + (t396 * t861 - t398 * t590) * t843) * qJD(3) + t16 * qJD(4) + t17 * qJD(5) + t24 * qJD(6), t708 + t1 * qJD(2) + t16 * qJD(3) + (Ifges(5,5) * t678 + Ifges(5,6) * t677 + m(7) * (t489 * t85 + t490 * t86) - t762 + t763 + t740 - t739 - t325 * mrSges(5,2) - t326 * mrSges(5,1) + t606 + (t683 - t685) * mrSges(7,3)) * qJD(4) + t10 * qJD(5) + t13 * qJD(6) + (m(6) * (t162 * t542 + t163 * t539) + (-t660 + t668) * mrSges(6,3)) * t746, t707 + t4 * qJD(2) + t17 * qJD(3) + t10 * qJD(4) + (t606 - t743 - t744 - t764 + t765) * qJD(5) + t19 * qJD(6) + (m(7) * (t538 * t82 + t541 * t81) + (-t663 + t672) * mrSges(7,3)) * t745, t701 + t8 * qJD(2) + t24 * qJD(3) + t13 * qJD(4) + t19 * qJD(5) + (t654 - t768 - t770) * qJD(6); -qJD(3) * t21 + qJD(4) * t2 + qJD(5) * t3 + qJD(6) * t9 - t714, qJD(4) * t29 + qJD(5) * t31 + qJD(6) * t38, -t704 (m(7) * (-t233 * t489 + t490 * t925) + t535 - t748 - mrSges(5,1) * t673 + mrSges(5,2) * t674 + (-t682 - t684) * mrSges(7,3) + t941) * qJD(4) + t36 * qJD(5) + t946 + (m(6) * (-t424 * t542 + t539 * t864) + (-t539 * t590 - t542 * t861) * mrSges(6,3)) * t746 + t577, t36 * qJD(4) + t941 * qJD(5) + t946 + (m(7) * t934 + (-t662 - t671) * mrSges(7,3)) * t745 + t576 (t649 + t926) * qJD(6) + t596 + t886 * t944; qJD(2) * t21 - qJD(4) * t15 + qJD(5) * t18 + qJD(6) * t25 - t700, t704, 0 (-t517 + t580) * qJD(4) + t89 * qJD(5) + 0.2e1 * (t591 * t841 + (t539 * t861 - t542 * t590) * t644) * qJD(4) - t595 + t853, t89 * qJD(4) + (t580 + t777) * qJD(5) + t594 + t853, -qJD(6) * t902 + t575 * t886 + t698; -qJD(2) * t2 + qJD(3) * t15 + qJD(5) * t11 + qJD(6) * t14 - t708, qJD(5) * t37 - t577 + t940, -qJD(5) * t94 + t595 + t854, qJD(5) * t215 + t646 ((t497 * t541 + t498 * t538) * t836 + t566) * qJD(5) + t197 * qJD(6) + t570, t197 * qJD(5) + t569 + t646; -qJD(2) * t3 - qJD(3) * t18 - qJD(4) * t11 + qJD(6) * t20 - t707, -qJD(4) * t37 - t576 + t940, qJD(4) * t94 - t594 + t854, qJD(6) * t196 - t570, -t509, -t509 - t562; -qJD(2) * t9 - qJD(3) * t25 - qJD(4) * t14 - qJD(5) * t20 - t701, -t886 * t932 - t596, -t572 * t886 - t698, -qJD(5) * t196 - t569, t562, 0;];
Cq  = t26;
