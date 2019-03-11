% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:23
% EndTime: 2019-03-09 07:08:54
% DurationCPUTime: 19.62s
% Computational Cost: add. (61879->745), mult. (120451->974), div. (0->0), fcn. (146710->10), ass. (0->421)
t476 = cos(qJ(5));
t738 = t476 * pkin(5);
t465 = -pkin(4) - t738;
t472 = sin(pkin(11));
t473 = cos(pkin(11));
t745 = sin(qJ(3));
t748 = cos(qJ(3));
t441 = -t472 * t745 + t473 * t748;
t442 = t472 * t748 + t473 * t745;
t475 = sin(qJ(4));
t747 = cos(qJ(4));
t394 = -t747 * t441 + t442 * t475;
t474 = sin(qJ(5));
t744 = sin(qJ(6));
t746 = cos(qJ(6));
t823 = t744 * t474 - t746 * t476;
t875 = t823 * t394;
t703 = t875 * mrSges(7,2);
t513 = t474 * t746 + t476 * t744;
t872 = t513 * t394;
t705 = t872 * mrSges(7,1);
t934 = t703 - t705;
t942 = t465 * t934;
t399 = mrSges(7,1) * t823 + mrSges(7,2) * t513;
t722 = Ifges(6,4) * t474;
t456 = Ifges(6,1) * t476 - t722;
t467 = Ifges(6,4) * t476;
t455 = Ifges(6,1) * t474 + t467;
t625 = t476 * t455;
t453 = Ifges(6,2) * t476 + t722;
t635 = t474 * t453;
t520 = t635 / 0.2e1 - t625 / 0.2e1;
t739 = pkin(5) * t474;
t750 = t476 / 0.2e1;
t829 = -Ifges(6,2) * t474 + t467;
t910 = t474 / 0.2e1;
t941 = t399 * t739 + t456 * t910 + t750 * t829 - t520;
t578 = -pkin(2) * t473 - pkin(1);
t417 = -pkin(3) * t441 + t578;
t519 = t475 * t441 + t442 * t747;
t238 = pkin(4) * t394 - pkin(9) * t519 + t417;
t728 = pkin(7) + qJ(2);
t550 = t728 * t473;
t551 = t728 * t472;
t409 = -t550 * t745 - t551 * t748;
t345 = -t442 * pkin(8) + t409;
t410 = t550 * t748 - t551 * t745;
t346 = t441 * pkin(8) + t410;
t832 = t345 * t475 + t747 * t346;
t123 = t476 * t238 - t474 * t832;
t124 = t474 * t238 + t476 * t832;
t266 = t823 * t519;
t269 = t513 * t519;
t147 = mrSges(7,1) * t269 - mrSges(7,2) * t266;
t236 = t747 * t345 - t475 * t346;
t657 = t519 * t474;
t180 = pkin(5) * t657 - t236;
t687 = t476 * mrSges(6,2);
t451 = t474 * mrSges(6,1) + t687;
t287 = t451 * t519;
t96 = -pkin(10) * t657 + t124;
t596 = t744 * t96;
t656 = t519 * t476;
t95 = -pkin(10) * t656 + t123;
t87 = t394 * pkin(5) + t95;
t54 = t746 * t87 - t596;
t597 = t746 * t96;
t55 = t744 * t87 + t597;
t784 = -t269 / 0.2e1;
t786 = -t266 / 0.2e1;
t874 = t451 * t394;
t893 = -mrSges(7,2) * t519 + mrSges(7,3) * t872;
t873 = t476 * t394;
t894 = mrSges(6,1) * t519 + mrSges(6,3) * t873;
t876 = t394 * t474;
t895 = -mrSges(6,2) * t519 + mrSges(6,3) * t876;
t897 = -pkin(5) * t876 + t832;
t918 = Ifges(7,4) * t875 + Ifges(7,2) * t872 + Ifges(7,6) * t519;
t919 = Ifges(7,1) * t875 + Ifges(7,4) * t872 + Ifges(7,5) * t519;
t920 = mrSges(7,1) * t519 - mrSges(7,3) * t875;
t940 = t123 * t894 + t124 * t895 + t897 * t147 + t180 * t934 + t236 * t874 + t832 * t287 + t54 * t920 + t55 * t893 + t918 * t784 + t919 * t786;
t821 = -mrSges(6,1) * t476 + t474 * mrSges(6,2);
t868 = t832 * t821;
t877 = t832 * mrSges(5,1);
t720 = Ifges(7,4) * t513;
t402 = -Ifges(7,2) * t823 + t720;
t901 = t872 * t402;
t907 = t236 * mrSges(5,2);
t925 = t897 * t399;
t434 = Ifges(7,4) * t823;
t404 = Ifges(7,1) * t513 - t434;
t926 = t875 * t404;
t936 = t823 * t918;
t937 = t513 * t919;
t939 = -t868 / 0.2e1 - (-t625 / 0.4e1 + t635 / 0.4e1) * t394 + t877 / 0.2e1 - t901 / 0.4e1 + t907 / 0.2e1 - t925 / 0.2e1 - t926 / 0.4e1 + t936 / 0.4e1 - t937 / 0.4e1;
t938 = t868 - t877 + t901 / 0.2e1 - t907 + t925 + t926 / 0.2e1 - t936 / 0.2e1 + t937 / 0.2e1;
t470 = t474 ^ 2;
t471 = t476 ^ 2;
t899 = (-t470 / 0.2e1 - t471 / 0.2e1) * mrSges(6,3);
t523 = t519 * t899;
t611 = t747 * pkin(3);
t464 = -t611 - pkin(4);
t449 = t464 - t738;
t793 = -pkin(10) - pkin(9);
t459 = t793 * t474;
t460 = t793 * t476;
t531 = t746 * t459 + t460 * t744;
t742 = pkin(3) * t475;
t463 = pkin(9) + t742;
t727 = pkin(10) + t463;
t439 = t727 * t474;
t440 = t727 * t476;
t532 = -t746 * t439 - t440 * t744;
t690 = t823 * mrSges(7,3);
t584 = -t690 / 0.2e1;
t585 = t690 / 0.2e1;
t642 = t464 * t451;
t740 = pkin(4) * t451;
t801 = m(7) / 0.2e1;
t689 = t513 * mrSges(7,3);
t582 = -t689 / 0.2e1;
t388 = t439 * t744 - t440 * t746;
t411 = -t459 * t744 + t460 * t746;
t867 = -t411 - t388;
t902 = t867 * t582;
t935 = (t449 + t465) * t739 * t801 - t740 / 0.2e1 + t642 / 0.2e1 - t902 + (t585 + t584) * (t532 + t531) + t941;
t613 = pkin(5) * t746;
t544 = t613 / 0.2e1;
t612 = pkin(5) * t744;
t794 = Ifges(6,3) / 0.2e1;
t932 = -t920 * t544 - t893 * t612 / 0.2e1 - t519 * t794;
t929 = t180 * t897;
t928 = t449 * t897;
t927 = t465 * t897;
t741 = pkin(4) * t519;
t301 = pkin(9) * t394 + t741;
t905 = t474 * t236;
t140 = t476 * t301 - t905;
t904 = t476 * t236;
t141 = t474 * t301 + t904;
t797 = m(7) * pkin(5);
t615 = t797 / 0.2e1;
t909 = pkin(10) * t876;
t105 = t141 + t909;
t896 = pkin(5) * t519 + pkin(10) * t873;
t93 = t140 + t896;
t71 = -t105 * t744 + t746 * t93;
t72 = t105 * t746 + t744 * t93;
t795 = -mrSges(6,2) / 0.2e1;
t796 = mrSges(6,1) / 0.2e1;
t783 = t875 / 0.2e1;
t785 = t872 / 0.2e1;
t522 = Ifges(7,5) * t783 + Ifges(7,6) * t785;
t841 = Ifges(7,3) * t519;
t911 = -mrSges(7,2) / 0.2e1;
t912 = mrSges(7,1) / 0.2e1;
t809 = t841 / 0.2e1 + t71 * t912 + t72 * t911 + t522;
t466 = Ifges(6,5) * t476;
t713 = Ifges(6,6) * t474;
t815 = (-t466 / 0.2e1 + t713 / 0.2e1) * t394;
t924 = -t140 * t796 - t141 * t795 - (t71 * t746 + t72 * t744) * t615 - t815 - t809 + t932;
t743 = pkin(3) * t442;
t243 = t301 + t743;
t127 = t476 * t243 - t905;
t128 = t474 * t243 + t904;
t560 = t876 / 0.2e1;
t88 = t127 + t896;
t97 = t128 + t909;
t58 = -t744 * t97 + t746 * t88;
t59 = t744 * t88 + t746 * t97;
t712 = Ifges(7,6) * t872;
t717 = Ifges(7,5) * t875;
t807 = -t841 / 0.2e1 - t712 / 0.2e1 - t717 / 0.2e1 + t59 * mrSges(7,2) / 0.2e1 - t58 * mrSges(7,1) / 0.2e1;
t923 = t807 - t127 * t796 - t128 * t795 - (t58 * t746 + t59 * t744) * t615 + Ifges(6,5) * t873 / 0.2e1 - Ifges(6,6) * t560 + t932;
t922 = t411 * mrSges(7,1) - t531 * mrSges(7,2);
t921 = t388 * mrSges(7,1) - t532 * mrSges(7,2);
t913 = t442 ^ 2;
t68 = t746 * t95 - t596;
t908 = t54 - t68;
t906 = t236 * t475;
t671 = t236 * t832;
t618 = -Ifges(7,5) * t823 - Ifges(7,6) * t513;
t81 = t618 + t921;
t903 = t81 * qJD(6);
t94 = t618 + t922;
t900 = t94 * qJD(6);
t753 = -t474 / 0.2e1;
t773 = -t519 / 0.2e1;
t776 = t394 / 0.2e1;
t830 = t466 - t713;
t855 = t519 / 0.2e1;
t857 = Ifges(6,6) * t519 - t394 * t829;
t858 = Ifges(6,5) * t519 - t394 * t456;
t888 = t858 * t750 + t857 * t753 + Ifges(7,5) * t786 + Ifges(7,6) * t784 + t830 * t855 + 0.2e1 * Ifges(5,4) * t773 + (Ifges(5,2) + Ifges(6,3) + Ifges(7,3)) * t776;
t887 = -pkin(4) / 0.2e1;
t778 = -t388 / 0.2e1;
t886 = -t394 / 0.2e1;
t763 = -t411 / 0.2e1;
t885 = -t823 / 0.2e1;
t884 = -t513 / 0.2e1;
t883 = t513 / 0.2e1;
t882 = m(7) * t739;
t881 = pkin(4) * t832;
t67 = -t744 * t95 - t597;
t851 = t55 + t67;
t828 = t470 + t471;
t880 = mrSges(6,3) * t828;
t871 = t821 * t519;
t401 = -Ifges(7,2) * t513 - t434;
t403 = -Ifges(7,1) * t823 - t720;
t538 = t402 * t884 + t403 * t883 + (t401 + t404) * t885;
t398 = mrSges(7,1) * t513 - mrSges(7,2) * t823;
t641 = t465 * t398;
t645 = t449 * t398;
t865 = t538 + t641 / 0.2e1 + t645 / 0.2e1;
t862 = t398 * qJD(6);
t553 = t402 / 0.4e1 - t403 / 0.4e1;
t554 = t401 / 0.4e1 + t404 / 0.4e1;
t861 = -t236 * t451 / 0.2e1 + t553 * t266 - t554 * t269;
t144 = -mrSges(7,1) * t266 - mrSges(7,2) * t269;
t704 = t269 * mrSges(7,3);
t192 = -mrSges(7,2) * t394 - t704;
t761 = t449 / 0.2e1;
t854 = t532 / 0.2e1;
t859 = t144 * t761 + t192 * t854;
t846 = mrSges(5,3) * t519;
t833 = t147 + t287;
t766 = -t531 / 0.2e1;
t825 = t266 * t763 - t269 * t766;
t782 = -t532 / 0.2e1;
t824 = t266 * t778 - t269 * t782;
t822 = -t123 * t474 + t124 * t476;
t552 = -t453 / 0.4e1 + t456 / 0.4e1;
t718 = Ifges(6,5) * t394;
t209 = t456 * t519 + t718;
t631 = t476 * t209;
t714 = Ifges(6,6) * t394;
t206 = t519 * t829 + t714;
t640 = t474 * t206;
t820 = -Ifges(5,4) * t394 + t417 * mrSges(5,2) + Ifges(5,1) * t855 - t640 / 0.2e1 + t631 / 0.2e1;
t297 = t394 * mrSges(6,1) - mrSges(6,3) * t656;
t751 = -t476 / 0.2e1;
t819 = t297 * t751 + t523;
t195 = mrSges(7,1) * t394 + t266 * mrSges(7,3);
t817 = t192 * t885 + t195 * t884;
t816 = t871 / 0.2e1 + t399 * t855;
t536 = t747 * t744;
t537 = t747 * t746;
t418 = (-t474 * t537 - t476 * t536) * pkin(3);
t419 = (-t474 * t536 + t476 * t537) * pkin(3);
t619 = t418 * t912 + t419 * t911;
t623 = t705 / 0.2e1 - t703 / 0.2e1;
t583 = t689 / 0.2e1;
t539 = t266 * t583 - t269 * t585 + t817;
t22 = t539 - t623;
t814 = t22 * qJD(1);
t756 = t465 / 0.2e1;
t561 = t144 * t756;
t765 = t531 / 0.2e1;
t787 = t195 / 0.2e1;
t812 = t192 * t765 + t411 * t787 - t871 * t887 + t561;
t757 = t464 / 0.2e1;
t811 = t388 * t787 - t757 * t871 + t859;
t721 = Ifges(7,4) * t266;
t112 = -Ifges(7,2) * t269 + t394 * Ifges(7,6) - t721;
t149 = -Ifges(7,1) * t269 + t721;
t775 = t394 / 0.4e1;
t806 = t180 * t398 / 0.2e1 + t618 * t775 - (t112 / 0.4e1 - t149 / 0.4e1) * t513;
t804 = -m(6) / 0.2e1;
t803 = m(6) / 0.2e1;
t802 = -m(7) / 0.2e1;
t800 = -pkin(9) / 0.2e1;
t799 = m(5) * pkin(3);
t798 = m(6) * pkin(3);
t792 = t112 / 0.2e1;
t262 = Ifges(7,4) * t269;
t115 = -Ifges(7,1) * t266 + Ifges(7,5) * t394 - t262;
t790 = t115 / 0.2e1;
t788 = t192 / 0.2e1;
t779 = t388 / 0.2e1;
t771 = t399 / 0.2e1;
t764 = t411 / 0.2e1;
t758 = -t463 / 0.2e1;
t752 = t474 / 0.4e1;
t749 = t476 / 0.4e1;
t737 = t54 * mrSges(7,2);
t736 = t55 * mrSges(7,1);
t733 = t67 * mrSges(7,1);
t732 = t68 * mrSges(7,2);
t726 = m(7) * qJD(2);
t719 = Ifges(5,5) * t394;
t715 = Ifges(5,6) * t519;
t294 = -mrSges(6,2) * t394 - mrSges(6,3) * t657;
t386 = t519 * mrSges(5,1);
t429 = t442 * mrSges(4,1);
t535 = t794 + Ifges(5,2) / 0.2e1 + Ifges(7,3) / 0.2e1;
t3 = (t578 * mrSges(4,2) + Ifges(4,4) * t441 + (-Ifges(4,2) + Ifges(4,1)) * t442) * t441 + t578 * t429 + (m(5) * t743 + t386) * t417 + m(6) * (t123 * t127 + t124 * t128 - t671) + (t712 + t717) * t776 + t872 * t792 + m(7) * (t54 * t58 + t55 * t59 + t929) - t913 * Ifges(4,4) + (mrSges(5,1) * t743 + t535 * t519 - t776 * t830 - t820) * t394 + t128 * t294 + t127 * t297 + t875 * t790 + (mrSges(5,2) * t743 + Ifges(5,1) * t886 + t888) * t519 + t59 * t192 + t58 * t195 + t940;
t702 = t3 * qJD(1);
t697 = t394 * mrSges(5,2);
t8 = (t815 + t522 - t820) * t394 + t115 * t783 + t112 * t785 + m(6) * (t123 * t140 + t124 * t141 - t671) + t141 * t294 + t140 * t297 + m(7) * (t54 * t71 + t55 * t72 + t929) + (t417 * mrSges(5,1) + (-Ifges(5,1) / 0.2e1 + t535) * t394 + t888) * t519 + t72 * t192 + t71 * t195 + t940;
t686 = t8 * qJD(1);
t288 = t519 * t453;
t289 = t455 * t519;
t148 = Ifges(7,2) * t266 - t262;
t622 = -Ifges(7,5) * t269 + Ifges(7,6) * t266;
t492 = t180 * t144 - (t148 / 0.2e1 + t790) * t269 + (t55 * mrSges(7,3) + t792 - t149 / 0.2e1) * t266 + t54 * t704 + t622 * t776;
t9 = t67 * t195 + m(7) * (t54 * t67 + t55 * t68) + t68 * t192 + t123 * t294 - t124 * t297 + t236 * t871 + ((-t718 / 0.2e1 + t123 * mrSges(6,3) + t288 / 0.2e1 - t209 / 0.2e1) * t474 + (-t714 / 0.2e1 - t124 * mrSges(6,3) - t206 / 0.2e1 - t289 / 0.2e1 + (m(7) * t180 + t147) * pkin(5)) * t476) * t519 + t492;
t685 = t9 * qJD(1);
t616 = t799 / 0.2e1;
t487 = (t127 * t476 + t474 * t128) * t803 + (t513 * t59 - t58 * t823) * t801 + t920 * t885 + t893 * t883 + t895 * t910 + t894 * t750 + t442 * t616;
t548 = t828 * t394;
t489 = (-t463 * t548 + t464 * t519) * t803 + (-t388 * t875 + t449 * t519 + t532 * t872) * t801 + (-t394 * t475 - t519 * t747) * t616;
t20 = t441 * mrSges(4,2) - t394 * t899 + t872 * t583 + t585 * t875 + t386 + t429 + t487 - t489 - t697 - t816;
t684 = qJD(1) * t20;
t681 = t127 * t474;
t680 = t128 * t476;
t14 = t54 * t192 - t55 * t195 + t492;
t679 = t14 * qJD(1);
t678 = t140 * t474;
t677 = t141 * t476;
t493 = (t744 * t875 + t746 * t872) * t615 + mrSges(6,1) * t560 - t873 * t795 + t623;
t629 = t476 * t294;
t558 = t629 / 0.2e1;
t637 = t474 * t297;
t495 = (-t513 * t908 - t851 * t823) * t801 - t637 / 0.2e1 + t558;
t16 = t266 * t582 - t269 * t584 + t493 - t495 - t817;
t676 = t16 * qJD(1);
t669 = t236 * t519;
t18 = t875 * t192 + t872 * t195 + (t441 ^ 2 + t913) * mrSges(4,3) - (-mrSges(5,3) * t394 + t629 - t637) * t394 + (t833 + t846) * t519 + m(7) * (t180 * t519 + t54 * t872 + t55 * t875) + m(6) * (-t394 * t822 - t669) + m(5) * (-t394 * t832 - t669) + m(4) * (-t409 * t442 + t410 * t441) + (m(3) * qJ(2) + mrSges(3,3)) * (t472 ^ 2 + t473 ^ 2);
t674 = t18 * qJD(1);
t400 = Ifges(7,5) * t513 - Ifges(7,6) * t823;
t661 = t519 * t400;
t452 = Ifges(6,5) * t474 + Ifges(6,6) * t476;
t659 = t519 * t452;
t644 = t463 * t474;
t643 = t463 * t476;
t639 = t474 * t858;
t638 = t474 * t894;
t632 = t476 * t857;
t630 = t476 * t895;
t614 = pkin(5) * t656;
t610 = mrSges(6,3) * t681;
t609 = mrSges(6,3) * t680;
t304 = -t418 * t823 + t419 * t513;
t608 = t304 * t801;
t605 = pkin(5) * t147 / 0.2e1;
t604 = t294 * t800;
t599 = t67 / 0.2e1 + t55 / 0.2e1;
t598 = t68 / 0.2e1 - t54 / 0.2e1;
t579 = qJD(4) * t608;
t576 = t744 * t195;
t563 = -t644 / 0.2e1;
t555 = t209 / 0.4e1 - t288 / 0.4e1;
t546 = mrSges(7,3) * t613;
t545 = mrSges(7,3) * t612;
t542 = -t611 / 0.2e1;
t534 = t269 * t546;
t533 = t266 * t545;
t526 = t474 * t542;
t525 = t680 - t681;
t524 = t677 - t678;
t521 = t582 * t872 + t584 * t875 + t880 * t886 + t816;
t488 = (t140 * t476 + t474 * t141) * t804 + (t513 * t72 - t71 * t823) * t802 + mrSges(5,1) * t773 + t823 * t920 / 0.2e1 + t893 * t884 + t895 * t753 + t894 * t751;
t490 = -t386 / 0.2e1 + (-pkin(9) * t548 - t741) * t803 + (-t411 * t875 + t465 * t519 + t531 * t872) * t801 + t521;
t24 = 0.2e1 * mrSges(5,2) * t776 + t488 + t490;
t516 = -t24 * qJD(1) + qJD(3) * t608;
t514 = t828 * t747;
t511 = t865 + t902;
t477 = (t452 + t400) * t519 / 0.4e1 - t939 + (t180 * t742 - t388 * t72 + t418 * t54 + t419 * t55 + t532 * t71 + t928) * t801 + t833 * t742 / 0.2e1 + t297 * t526 + (t677 / 0.2e1 - t678 / 0.2e1) * mrSges(6,3) + t419 * t788 + Ifges(5,6) * t773 + t418 * t787 - t874 * t757 + Ifges(5,5) * t886 + t72 * t584 + t857 * t749 + t858 * t752 + t71 * t582 + t920 * t854 + t934 * t761 + (t464 * t832 + t524 * t463 + (t747 * t822 - t906) * pkin(3)) * t803 + t558 * t611 + t893 * t778 + t895 * t643 / 0.2e1 + t894 * t563;
t478 = -t661 / 0.4e1 - t659 / 0.4e1 - t639 / 0.4e1 - t632 / 0.4e1 + (pkin(9) * t525 - t881) * t804 + t939 + t610 / 0.2e1 + (-t411 * t59 + t531 * t58 + t927) * t802 + t715 / 0.2e1 + t719 / 0.2e1 + t630 * t800 + t874 * t887 + t59 * t585 + t58 * t583 + t920 * t766 - t942 / 0.2e1 - t609 / 0.2e1 + t893 * t764 + pkin(9) * t638 / 0.2e1;
t2 = t478 + t477;
t486 = -t418 * t689 - t419 * t690 + (-mrSges(5,1) + t399 + t821) * t742 + (-mrSges(5,2) + t880) * t611;
t80 = m(7) * (-t388 * t419 + t418 * t532 + t449 * t742) + (t463 * t514 + t464 * t475) * t798 + t486;
t510 = -t2 * qJD(1) - t80 * qJD(3) - t304 * t726 / 0.2e1;
t508 = t538 + t941;
t44 = t449 * t882 + t508 + t642 + t645;
t479 = -t640 / 0.4e1 + t631 / 0.4e1 + t474 * t605 + t830 * t775 + t614 * t771 - t288 * t749 - t289 * t752 + t68 * t584 + t54 * t585 - (t148 + t115) * t823 / 0.4e1 + t552 * t656 + t851 * t582 - (t829 + t455) * t657 / 0.4e1 + t806 + t861;
t171 = t180 * t739;
t503 = t388 * t908 + t851 * t532 + t171;
t5 = (t449 * t614 + t503) * t801 + t479 + t294 * t563 + t824 * mrSges(7,3) + t819 * t463 + t811 + t923;
t509 = -t5 * qJD(1) - t44 * qJD(3);
t494 = -(t148 / 0.4e1 + t115 / 0.4e1) * t823 + t806;
t481 = -(mrSges(7,3) * t782 + t554) * t269 + (mrSges(7,3) * t778 + t553) * t266 + t195 * t779 + t494 + t859;
t11 = t481 + t807;
t73 = t538 + t645;
t507 = -t11 * qJD(1) - t73 * qJD(3);
t505 = -t513 * t599 - t598 * t823;
t502 = t411 * t908 + t851 * t531 + t171;
t491 = (t418 * t746 + t419 * t744) * t615 + mrSges(6,1) * t526 + t542 * t687 + t619;
t34 = t583 * t867 + t491 - t865 - t935;
t45 = t465 * t882 + t508 + t641 - t740;
t6 = (t465 * t614 + t502) * t801 + t479 + t474 * t604 + t825 * mrSges(7,3) + t819 * pkin(9) + t812 + t924;
t501 = -t6 * qJD(1) + t34 * qJD(3) - t45 * qJD(4);
t482 = (mrSges(7,3) * t763 + t553) * t266 - (mrSges(7,3) * t766 + t554) * t269 + t531 * t788 + t195 * t764 + t561 + t494;
t13 = t482 - t809;
t498 = -(t779 + t764) * t689 + t511;
t42 = t498 - t619;
t75 = t538 + t641;
t500 = -t13 * qJD(1) - t42 * qJD(3) - t75 * qJD(4);
t499 = -t513 * t545 + t546 * t823 + t618 + t830;
t497 = (-t455 / 0.4e1 - t829 / 0.4e1) * t519 - t206 / 0.4e1 - t289 / 0.4e1 - t714 / 0.4e1 + t605;
t104 = (t766 + t765) * mrSges(7,2) + (t764 + t763) * mrSges(7,1);
t19 = -t598 * mrSges(7,2) + t599 * mrSges(7,1) + (-t746 * t192 / 0.2e1 + t576 / 0.2e1 + (t744 * t786 + t746 * t784) * mrSges(7,3)) * pkin(5);
t448 = (mrSges(7,1) * t744 + mrSges(7,2) * t746) * pkin(5);
t86 = (t782 + t854) * mrSges(7,2) + (t779 + t778) * mrSges(7,1);
t496 = t19 * qJD(1) - t86 * qJD(3) - t104 * qJD(4) + t448 * qJD(5);
t485 = t466 * t775 + t494 + t861;
t443 = t448 * qJD(6);
t43 = t498 + t619;
t35 = t491 + t511 + t935;
t25 = -t488 + t697 / 0.2e1 + mrSges(5,2) * t886 + t490;
t23 = t539 + t623;
t21 = t487 + t489 + t521;
t17 = t493 + t495 + t539;
t15 = -t736 / 0.2e1 + t192 * t544 + t533 / 0.2e1 - pkin(5) * t576 / 0.2e1 + t534 / 0.2e1 - t737 / 0.2e1 + t733 / 0.2e1 - t732 / 0.2e1 + t622;
t12 = t482 + t809;
t10 = t481 - t807;
t7 = (t505 + t825) * mrSges(7,3) + (t604 + t497) * t474 + (t297 * t800 + ((m(7) * t756 + t771) * pkin(5) + t552) * t519 + t555) * t476 + t502 * t801 + t485 + pkin(9) * t523 + t812 - t924;
t4 = (t505 + t824) * mrSges(7,3) + (t294 * t758 + t497) * t474 + (t297 * t758 + ((m(7) * t761 + t771) * pkin(5) + t552) * t519 + t555) * t476 + t463 * t523 + t485 + t503 * t801 + t811 - t923;
t1 = -t478 + t477;
t26 = [qJD(2) * t18 + qJD(3) * t3 + qJD(4) * t8 + qJD(5) * t9 + qJD(6) * t14, t674 + (t513 * t875 - t823 * t872) * t726 + t21 * qJD(3) + t25 * qJD(4) + t17 * qJD(5) + t23 * qJD(6), t21 * qJD(2) + t1 * qJD(4) + t4 * qJD(5) + t10 * qJD(6) + t702 + (t661 / 0.2e1 + t659 / 0.2e1 + t639 / 0.2e1 + t632 / 0.2e1 + t938 - t610 - t742 * t846 + m(7) * (-t388 * t59 + t532 * t58 + t928) - t715 - t719 + (m(6) * t525 + t630 - t638) * t463 + (m(6) * t832 - t874) * t464 - t58 * t689 - t59 * t690 + t532 * t920 + t449 * t934 - Ifges(4,6) * t442 + t609 + Ifges(4,5) * t441 - t409 * mrSges(4,2) - t410 * mrSges(4,1) + (-t747 * t832 + t906) * t799 - t388 * t893 - (-mrSges(5,3) * t611 - t520) * t394) * qJD(3), t25 * qJD(2) + t1 * qJD(3) + t7 * qJD(5) + t12 * qJD(6) + t686 + (-t71 * t689 - t72 * t690 + t942 - t411 * t893 + t531 * t920 + pkin(4) * t874 + (pkin(9) * t895 + t141 * mrSges(6,3) + t857 / 0.2e1) * t476 + (-pkin(9) * t894 - t140 * mrSges(6,3) + t858 / 0.2e1) * t474 + 0.2e1 * (pkin(9) * t524 - t881) * t803 + 0.2e1 * (-t411 * t72 + t531 * t71 + t927) * t801 + (t452 / 0.2e1 + t400 / 0.2e1 - Ifges(5,6)) * t519 + (-Ifges(5,5) + t520) * t394 + t938) * qJD(4), t685 + t17 * qJD(2) + t4 * qJD(3) + t7 * qJD(4) + (-Ifges(6,6) * t656 - Ifges(6,5) * t657 + t733 + (t67 * t746 + t68 * t744) * t797 + t533 + t534 - t732 - t123 * mrSges(6,2) - t124 * mrSges(6,1) + t622) * qJD(5) + t15 * qJD(6), t679 + t23 * qJD(2) + t10 * qJD(3) + t12 * qJD(4) + t15 * qJD(5) + (t622 - t736 - t737) * qJD(6); qJD(3) * t20 - qJD(4) * t24 - qJD(5) * t16 + qJD(6) * t22 - t674, 0, t579 + t684, t516, -t676 - t862 + (-t398 - t451 + (-t513 * t613 - t612 * t823) * m(7)) * qJD(5), -qJD(5) * t398 + t814 - t862; -qJD(2) * t20 + qJD(4) * t2 + qJD(5) * t5 + qJD(6) * t11 - t702, t579 - t684, qJD(4) * t80 + qJD(5) * t44 + qJD(6) * t73 (m(7) * (-t411 * t419 + t418 * t531 + t465 * t742) + (-pkin(4) * t475 + pkin(9) * t514) * t798 + t486) * qJD(4) + t35 * qJD(5) + t43 * qJD(6) - t510, t35 * qJD(4) + ((t388 * t746 + t532 * t744) * t797 + mrSges(6,2) * t644 - mrSges(6,1) * t643 + t499 + t921) * qJD(5) + t903 - t509, t43 * qJD(4) + t81 * qJD(5) - t507 + t903; qJD(2) * t24 - qJD(3) * t2 + qJD(5) * t6 + qJD(6) * t13 - t686, -t516, -t34 * qJD(5) + t42 * qJD(6) + t510, qJD(5) * t45 + qJD(6) * t75 ((t411 * t746 + t531 * t744) * t797 + t499 + t821 * pkin(9) + t922) * qJD(5) + t900 - t501, t94 * qJD(5) - t500 + t900; qJD(2) * t16 - qJD(3) * t5 - qJD(4) * t6 - qJD(6) * t19 - t685, t676, t34 * qJD(4) + t86 * qJD(6) + t509, t104 * qJD(6) + t501, -t443, -t443 - t496; -qJD(2) * t22 - qJD(3) * t11 - qJD(4) * t13 + qJD(5) * t19 - t679, -t814, -t42 * qJD(4) - t86 * qJD(5) + t507, -t104 * qJD(5) + t500, t496, 0;];
Cq  = t26;
