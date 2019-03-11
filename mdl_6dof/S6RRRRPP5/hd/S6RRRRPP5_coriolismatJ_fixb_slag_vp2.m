% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRRPP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:21
% EndTime: 2019-03-09 21:06:01
% DurationCPUTime: 21.52s
% Computational Cost: add. (23220->883), mult. (49652->1110), div. (0->0), fcn. (48403->6), ass. (0->438)
t963 = -Ifges(7,4) - Ifges(6,5);
t951 = Ifges(6,4) + Ifges(5,5);
t967 = Ifges(7,2) + Ifges(6,3);
t950 = Ifges(5,6) - Ifges(6,6);
t837 = cos(qJ(3));
t686 = t837 * pkin(8);
t516 = pkin(9) * t837 + t686;
t835 = sin(qJ(3));
t683 = t835 * pkin(8);
t601 = -pkin(9) * t835 - t683;
t834 = sin(qJ(4));
t836 = cos(qJ(4));
t899 = t836 * t516 + t834 * t601;
t917 = t899 * mrSges(6,1);
t918 = t899 * mrSges(5,1);
t380 = -t834 * t516 + t836 * t601;
t919 = t380 * mrSges(6,3);
t920 = t380 * mrSges(5,2);
t500 = t834 * t835 - t836 * t837;
t911 = t500 * qJ(6) + t899;
t949 = t911 * mrSges(7,1);
t501 = -t834 * t837 - t835 * t836;
t269 = t501 * qJ(6) - t380;
t962 = t269 * mrSges(7,2);
t966 = t919 - t917 - t918 - t920 - t949 - t962;
t965 = t917 / 0.2e1 + t918 / 0.2e1 - t919 / 0.2e1 + t920 / 0.2e1 + t949 / 0.2e1 + t962 / 0.2e1;
t964 = -Ifges(5,1) - Ifges(6,1) - Ifges(7,1);
t550 = sin(qJ(2));
t460 = t501 * t550;
t850 = -t460 / 0.2e1;
t862 = t911 / 0.2e1;
t867 = pkin(4) + pkin(5);
t708 = t501 * t867;
t748 = qJ(5) * t500;
t262 = -t748 + t708;
t684 = t835 * pkin(3);
t576 = -t684 + t262;
t571 = m(7) * t576;
t961 = t963 * t501;
t960 = t963 * t500;
t459 = t500 * t550;
t440 = Ifges(7,4) * t459;
t286 = Ifges(7,1) * t460 - t440;
t438 = Ifges(6,5) * t459;
t287 = Ifges(6,1) * t460 - t438;
t822 = Ifges(5,4) * t459;
t288 = Ifges(5,1) * t460 + t822;
t551 = cos(qJ(2));
t238 = -t551 * Ifges(6,6) - Ifges(6,3) * t460 - t438;
t753 = t551 * Ifges(7,6);
t240 = -Ifges(7,2) * t460 - t440 + t753;
t903 = t240 + t238;
t959 = t288 + t287 + t286 + t903;
t821 = Ifges(5,4) * t501;
t914 = t500 * t967 + t961;
t958 = t500 * t964 + t821 + t914 + t961;
t685 = t836 * pkin(3);
t534 = -t685 - pkin(4);
t527 = -pkin(5) + t534;
t682 = t834 * pkin(3);
t530 = t682 + qJ(5);
t955 = -t530 * t269 + t527 * t911;
t756 = t551 * mrSges(7,2);
t773 = t460 * mrSges(7,3);
t384 = -t756 - t773;
t777 = t459 * mrSges(7,3);
t388 = mrSges(7,1) * t551 + t777;
t954 = -t269 * t384 / 0.2e1 + t388 * t862;
t547 = t550 * pkin(7);
t657 = t550 * t835;
t511 = pkin(3) * t657 + t547;
t749 = qJ(5) * t459;
t623 = -t511 - t749;
t129 = t460 * t867 + t623;
t831 = pkin(8) * t550;
t513 = -pkin(2) * t551 - pkin(1) - t831;
t499 = t837 * t513;
t658 = t550 * t837;
t639 = pkin(9) * t658;
t313 = -t639 + t499 + (-pkin(7) * t835 - pkin(3)) * t551;
t656 = t551 * t837;
t430 = pkin(7) * t656 + t513 * t835;
t367 = -pkin(9) * t657 + t430;
t320 = t836 * t367;
t141 = t313 * t834 + t320;
t203 = -pkin(4) * t460 - t623;
t535 = -pkin(3) * t837 - pkin(2);
t602 = t501 * qJ(5) + t535;
t219 = -t500 * t867 - t602;
t242 = Ifges(5,2) * t460 - t551 * Ifges(5,6) - t822;
t274 = -mrSges(6,1) * t459 - mrSges(6,3) * t460;
t275 = t459 * mrSges(7,1) + t460 * mrSges(7,2);
t276 = -mrSges(5,1) * t459 + mrSges(5,2) * t460;
t299 = t500 * pkin(4) + t602;
t335 = -mrSges(6,1) * t501 + mrSges(6,3) * t500;
t336 = t501 * mrSges(7,1) - t500 * mrSges(7,2);
t337 = -mrSges(5,1) * t501 - mrSges(5,2) * t500;
t340 = -Ifges(7,5) * t500 - Ifges(7,6) * t501;
t348 = -Ifges(5,2) * t500 - t821;
t632 = -t500 * t951 + t950 * t501;
t761 = t501 * mrSges(5,3);
t921 = mrSges(6,2) + mrSges(5,3);
t493 = Ifges(5,4) * t500;
t888 = t501 * t964 - t493 - t960;
t942 = Ifges(5,2) * t501 - t493 + t888;
t943 = -t501 * t967 + t960;
t442 = Ifges(5,4) * t460;
t285 = Ifges(5,2) * t459 + t442;
t754 = t551 * Ifges(7,5);
t819 = Ifges(7,4) * t460;
t244 = -Ifges(7,1) * t459 + t754 - t819;
t816 = Ifges(6,5) * t460;
t246 = -Ifges(6,1) * t459 - t551 * Ifges(6,4) - t816;
t248 = -Ifges(5,1) * t459 - t551 * Ifges(5,5) + t442;
t890 = t248 + t246 + t244;
t944 = t285 + t890;
t283 = -Ifges(6,3) * t459 + t816;
t284 = -Ifges(7,2) * t459 + t819;
t945 = t284 + t283;
t953 = t921 * (t899 * t459 / 0.2e1 + t380 * t850) + t141 * t761 / 0.2e1 + t535 * t276 / 0.2e1 + t511 * t337 / 0.2e1 + t203 * t335 / 0.2e1 + t129 * t336 / 0.2e1 + t299 * t274 / 0.2e1 + t219 * t275 / 0.2e1 + (t340 / 0.4e1 - t632 / 0.4e1) * t551 + (t945 / 0.4e1 - t944 / 0.4e1) * t500 + (-t943 / 0.4e1 + t942 / 0.4e1) * t460 + (t242 / 0.4e1 - t959 / 0.4e1) * t501 + (t348 / 0.4e1 - t958 / 0.4e1) * t459;
t824 = mrSges(7,3) * t501;
t947 = t911 * t824;
t338 = mrSges(6,1) * t500 + mrSges(6,3) * t501;
t941 = -m(6) * t299 - t338;
t825 = mrSges(6,3) + mrSges(7,2);
t868 = m(6) + m(7);
t521 = qJ(5) * t868 + t825;
t940 = qJD(4) * t521;
t939 = t521 * qJD(5);
t764 = t500 * mrSges(7,3);
t475 = t764 / 0.2e1;
t767 = t500 * mrSges(6,2);
t477 = -t767 / 0.2e1;
t870 = -mrSges(7,3) / 0.2e1;
t871 = mrSges(6,2) / 0.2e1;
t680 = t871 + t870;
t875 = m(7) / 0.2e1;
t935 = m(7) * t862 - t500 * t680 + t911 * t875 + t475 + t477;
t882 = 0.2e1 * m(7);
t931 = mrSges(6,1) / 0.2e1;
t930 = mrSges(7,1) / 0.2e1;
t922 = m(5) * t535;
t739 = t380 * t500;
t736 = t899 * t501;
t461 = t501 * t551;
t698 = t551 * t500;
t916 = t963 * t698 + (Ifges(6,6) - Ifges(7,6)) * t550 - t967 * t461;
t339 = -mrSges(7,1) * t500 - mrSges(7,2) * t501;
t910 = m(7) * t219 + t339;
t757 = t550 * mrSges(7,1);
t769 = t698 * mrSges(7,3);
t391 = -t757 + t769;
t909 = t391 / 0.2e1;
t634 = t682 / 0.2e1;
t906 = mrSges(6,2) - mrSges(7,3);
t755 = t551 * mrSges(6,3);
t776 = t460 * mrSges(6,2);
t383 = -t755 + t776;
t774 = t460 * mrSges(5,3);
t385 = mrSges(5,2) * t551 + t774;
t902 = t385 + t383;
t382 = mrSges(6,2) * t461 + mrSges(6,3) * t550;
t386 = mrSges(7,2) * t550 - mrSges(7,3) * t461;
t901 = t386 + t382;
t778 = t459 * mrSges(5,3);
t389 = -mrSges(5,1) * t551 + t778;
t779 = t459 * mrSges(6,2);
t390 = mrSges(6,1) * t551 - t779;
t900 = -t389 + t390;
t842 = t530 / 0.2e1;
t896 = mrSges(6,2) * t842 + t530 * t870;
t759 = t534 * mrSges(6,2);
t760 = t527 * mrSges(7,3);
t895 = t760 / 0.2e1 - t759 / 0.2e1;
t666 = -t773 / 0.2e1;
t894 = t666 + t776 / 0.2e1;
t592 = -t837 * mrSges(4,1) + t835 * mrSges(4,2);
t633 = t459 * t950 + t460 * t951;
t893 = t736 + t739;
t892 = t760 - t759;
t891 = m(6) / 0.4e1 + m(7) / 0.4e1;
t889 = t964 * t698 + (-Ifges(7,5) + t951) * t550 + (Ifges(5,4) + t963) * t461;
t887 = t530 * t906;
t433 = t460 * qJ(6);
t702 = t551 * qJ(5);
t573 = t141 - 0.2e1 * t702;
t876 = -m(7) / 0.2e1;
t879 = -m(6) / 0.2e1;
t886 = -t551 * t825 - t573 * t879 - (-t433 + t573) * t876;
t885 = -t836 * mrSges(5,2) + (-mrSges(5,1) - mrSges(6,1) - mrSges(7,1)) * t834;
t832 = pkin(2) * t550;
t517 = -pkin(8) * t551 + t832;
t431 = pkin(7) * t657 + t837 * t517;
t316 = t550 * pkin(3) - pkin(9) * t656 + t431;
t432 = -pkin(7) * t658 + t835 * t517;
t655 = t551 * t835;
t372 = -pkin(9) * t655 + t432;
t152 = t834 * t316 + t836 * t372;
t135 = t550 * qJ(5) + t152;
t619 = -t316 * t836 + t834 * t372;
t136 = -t550 * pkin(4) + t619;
t77 = qJ(6) * t698 - t550 * t867 + t619;
t89 = -qJ(6) * t461 + t135;
t884 = t152 * mrSges(5,2) / 0.2e1 + t619 * mrSges(5,1) / 0.2e1 + t136 * t931 - t135 * mrSges(6,3) / 0.2e1 - t89 * mrSges(7,2) / 0.2e1 + t77 * t930;
t880 = 0.2e1 * qJD(4);
t878 = m(6) / 0.2e1;
t873 = m(5) * pkin(3);
t872 = -mrSges(6,2) / 0.2e1;
t94 = t141 - t433;
t869 = t94 / 0.2e1;
t866 = pkin(4) * mrSges(6,2);
t864 = t242 / 0.2e1;
t277 = -mrSges(6,1) * t460 + mrSges(6,3) * t459;
t861 = t277 / 0.2e1;
t278 = mrSges(7,1) * t460 - mrSges(7,2) * t459;
t860 = t278 / 0.2e1;
t859 = t338 / 0.2e1;
t858 = t339 / 0.2e1;
t855 = -t389 / 0.2e1;
t854 = t390 / 0.2e1;
t758 = t550 * mrSges(6,1);
t771 = t698 * mrSges(6,2);
t393 = -t758 - t771;
t853 = t393 / 0.2e1;
t852 = -t459 / 0.2e1;
t851 = t460 / 0.2e1;
t849 = t461 / 0.2e1;
t848 = -t461 / 0.2e1;
t847 = -t698 / 0.2e1;
t846 = -t500 / 0.2e1;
t845 = t500 / 0.2e1;
t844 = -t501 / 0.2e1;
t843 = t501 / 0.2e1;
t841 = -t550 / 0.2e1;
t840 = t550 / 0.2e1;
t839 = -t551 / 0.2e1;
t838 = t551 / 0.2e1;
t833 = m(6) * t899;
t830 = t459 * pkin(5);
t548 = t551 * pkin(7);
t651 = t834 * t367;
t140 = t836 * t313 - t651;
t435 = t459 * qJ(6);
t93 = -t435 + t140;
t827 = t93 * mrSges(7,2);
t826 = t94 * mrSges(7,1);
t823 = Ifges(3,4) * t550;
t820 = Ifges(6,4) * t698;
t817 = Ifges(5,5) * t698;
t814 = Ifges(7,5) * t698;
t813 = Ifges(6,2) * t550;
t812 = Ifges(5,6) * t461;
t811 = Ifges(6,6) * t461;
t810 = Ifges(7,6) * t461;
t809 = Ifges(4,3) * t550;
t808 = Ifges(5,3) * t550;
t807 = Ifges(7,3) * t550;
t429 = -pkin(7) * t655 + t499;
t366 = -t639 + t429;
t163 = t366 * t834 + t320;
t106 = t163 - t433;
t806 = t106 * mrSges(7,1);
t164 = t836 * t366 - t651;
t107 = -t435 + t164;
t805 = t107 * mrSges(7,2);
t802 = t140 * mrSges(5,2);
t801 = t140 * mrSges(6,3);
t800 = t141 * mrSges(5,1);
t799 = t141 * mrSges(6,1);
t796 = t163 * mrSges(5,1);
t795 = t163 * mrSges(6,1);
t794 = t164 * mrSges(5,2);
t793 = t164 * mrSges(6,3);
t772 = t461 * mrSges(7,1);
t770 = t698 * mrSges(7,2);
t124 = t141 - t702;
t127 = t551 * pkin(4) - t140;
t512 = pkin(3) * t655 + t548;
t622 = -qJ(5) * t698 - t512;
t130 = t461 * t867 + t622;
t204 = -pkin(4) * t461 - t622;
t243 = -Ifges(5,4) * t698 + Ifges(5,2) * t461 + t550 * Ifges(5,6);
t279 = -mrSges(5,1) * t460 - mrSges(5,2) * t459;
t280 = -mrSges(6,1) * t461 + mrSges(6,3) * t698;
t281 = -t770 + t772;
t282 = -mrSges(5,1) * t461 - mrSges(5,2) * t698;
t387 = -mrSges(5,2) * t550 + mrSges(5,3) * t461;
t392 = mrSges(5,1) * t550 + mrSges(5,3) * t698;
t679 = Ifges(4,4) * t837;
t597 = -Ifges(4,2) * t835 + t679;
t581 = t597 * t550;
t452 = -Ifges(4,6) * t551 + t581;
t453 = Ifges(4,6) * t550 + t551 * t597;
t678 = Ifges(4,4) * t835;
t599 = Ifges(4,1) * t837 - t678;
t582 = t599 * t550;
t454 = -Ifges(4,5) * t551 + t582;
t455 = Ifges(4,5) * t550 + t551 * t599;
t593 = mrSges(4,1) * t835 + mrSges(4,2) * t837;
t484 = t593 * t551;
t637 = mrSges(4,3) * t657;
t507 = t551 * mrSges(4,2) - t637;
t508 = -t550 * mrSges(4,2) - mrSges(4,3) * t655;
t638 = mrSges(4,3) * t658;
t509 = -t551 * mrSges(4,1) - t638;
t510 = t550 * mrSges(4,1) - mrSges(4,3) * t656;
t578 = t593 * t547;
t595 = Ifges(4,5) * t837 - Ifges(4,6) * t835;
t583 = t551 * t595;
t625 = -t655 / 0.2e1;
t626 = t656 / 0.2e1;
t627 = -t657 / 0.2e1;
t628 = t658 / 0.2e1;
t703 = t550 * t551;
t76 = pkin(5) * t551 + t127 + t435;
t81 = -t433 + t124;
t5 = t889 * t852 + t890 * t847 + t903 * t848 + t916 * t850 + m(5) * (-t140 * t619 + t141 * t152 + t511 * t512) - t619 * t389 + t551 * t578 + m(4) * (pkin(7) ^ 2 * t703 + t429 * t431 + t430 * t432) + t243 * t851 + t242 * t849 + (t550 * t595 - t823 + t950 * t460 - t951 * t459 + (Ifges(3,1) - Ifges(6,2) - Ifges(4,3) - Ifges(5,3)) * t551) * t840 + t484 * t547 + (-Ifges(7,5) * t459 - Ifges(7,6) * t460 + t823 + (Ifges(3,2) + Ifges(7,3)) * t551) * t841 - pkin(1) * (mrSges(3,1) * t550 + mrSges(3,2) * t551) + t429 * t510 + t511 * t282 + t512 * t279 + t432 * t507 + t430 * t508 + t431 * t509 + t152 * t385 + t81 * t386 + t141 * t387 + t77 * t388 + t136 * t390 + t76 * t391 + t140 * t392 + t127 * t393 + t124 * t382 + t135 * t383 + t89 * t384 + m(6) * (t124 * t135 + t127 * t136 + t203 * t204) + m(7) * (t129 * t130 + t76 * t77 + t81 * t89) + t130 * t278 + t203 * t280 + t129 * t281 + t204 * t277 + (0.2e1 * Ifges(3,4) * t551 - t807 - t810 - t814 + (Ifges(3,1) - Ifges(3,2)) * t550) * t838 + (-t811 + t813 - t820 + t808 + t812 - t817 + t583 + t809) * t839 + t452 * t625 + t454 * t626 + t453 * t627 + t455 * t628;
t768 = t5 * qJD(1);
t765 = t500 * mrSges(5,3);
t762 = t501 * mrSges(6,2);
t752 = t867 * mrSges(7,3);
t434 = qJ(5) * t460;
t273 = -t459 * pkin(4) - t434;
t640 = pkin(3) * t658;
t222 = t640 + t273;
t146 = -t222 + t830;
t567 = t124 * t779 + t129 * t275 - t140 * t774 + t203 * t274 + t511 * t276 + t633 * t839 - t76 * t773;
t594 = Ifges(4,5) * t835 + Ifges(4,6) * t837;
t596 = Ifges(4,2) * t837 + t678;
t598 = Ifges(4,1) * t835 + t679;
t616 = Ifges(7,5) * t460 - Ifges(7,6) * t459;
t629 = -t658 / 0.2e1;
t688 = t835 / 0.2e1;
t691 = -t837 / 0.2e1;
t6 = t567 + t459 * t864 - t81 * t777 + t127 * t776 + t594 * t703 / 0.2e1 + t616 * t838 + t141 * t778 + t106 * t388 + t107 * t384 + t146 * t278 + m(7) * (t106 * t76 + t107 * t81 + t129 * t146) + t454 * t627 + t452 * t629 + t945 * t850 + (m(5) * t511 + t279) * t640 + (-t509 - t638) * t430 + (t507 + t637) * t429 + (m(6) * t203 + t277) * t222 + (m(5) * t141 + m(6) * t124 + t902) * t164 + (-m(5) * t140 + m(6) * t127 + t900) * t163 + (-pkin(7) * t592 + t596 * t688 + t598 * t691) * t550 ^ 2 + t944 * t851 + t959 * t852;
t751 = t6 * qJD(1);
t174 = -t273 + t830;
t7 = t567 + t94 * t388 + t93 * t384 + t174 * t278 + t273 * t277 + t900 * t141 + t902 * t140 + m(6) * (t124 * t140 + t127 * t141 + t203 * t273) + m(7) * (t129 * t174 + t76 * t94 + t81 * t93) - (-t127 * mrSges(6,2) - t754 / 0.2e1 - t285 / 0.2e1 - t248 / 0.2e1 - t246 / 0.2e1 - t244 / 0.2e1 + t284 / 0.2e1 + t283 / 0.2e1) * t460 + (-t81 * mrSges(7,3) + t141 * mrSges(5,3) - t753 / 0.2e1 - t288 / 0.2e1 - t287 / 0.2e1 - t286 / 0.2e1 - t240 / 0.2e1 - t238 / 0.2e1 + t864) * t459;
t750 = t7 * qJD(1);
t22 = (-t383 - t384) * t551 + (t277 - t278) * t459 + m(6) * (-t124 * t551 + t203 * t459) + m(7) * (-t129 * t459 - t551 * t81);
t747 = qJD(1) * t22;
t31 = m(7) * (t459 * t76 - t460 * t81) + t459 * t388 - t460 * t384;
t746 = qJD(1) * t31;
t729 = t459 * t527;
t722 = t460 * t530;
t697 = t867 * t459;
t692 = t825 * t685;
t290 = mrSges(6,2) * t739;
t671 = t824 / 0.2e1;
t670 = -t777 / 0.2e1;
t668 = -t776 / 0.2e1;
t667 = t773 / 0.2e1;
t665 = -t764 / 0.2e1;
t664 = -t762 / 0.2e1;
t659 = -t752 / 0.2e1;
t653 = t836 * t530;
t643 = mrSges(5,3) * t682;
t641 = t500 * t685;
t636 = t685 / 0.2e1;
t635 = -t682 / 0.2e1;
t621 = t685 * t774;
t620 = t459 * t643;
t617 = mrSges(5,1) * t500 - mrSges(5,2) * t501;
t334 = -t501 * pkin(4) + t748;
t559 = mrSges(5,3) * t736 + t219 * t336 + t299 * t335 + t535 * t337 + t348 * t843 + t844 * t958 + t943 * t845 + t942 * t846 + t290 - t947;
t574 = t835 * t596;
t575 = t837 * t598;
t584 = t684 + t334;
t604 = pkin(3) * t617;
t10 = -t559 - t684 * t922 - t219 * t571 - t575 / 0.2e1 + t597 * t691 + t574 / 0.2e1 + pkin(2) * t593 - t576 * t339 - t947 + (-t599 / 0.2e1 - t604) * t835 + (-t739 + t893) * mrSges(5,3) + (-t736 + t893) * mrSges(6,2) + t941 * t584;
t609 = -t163 * t380 + t164 * t899;
t553 = t902 * t380 / 0.2e1 - (t835 ^ 2 + t837 ^ 2) * mrSges(4,3) * t831 / 0.2e1 + t954 + m(5) * (-t899 * t140 + t380 * t141 + (t511 * t835 + t535 * t658) * pkin(3) + t609) / 0.2e1 + t899 * t854 + t899 * t855 + (t380 * t124 + t127 * t899 + t203 * t584 + t299 * t222 + t609) * t878 - (t581 + t452) * t835 / 0.4e1 + (t582 + t454) * t837 / 0.4e1 - t583 / 0.4e1 + t592 * t832 / 0.2e1 + t140 * t765 / 0.2e1 - t509 * t686 / 0.2e1 + t279 * t684 / 0.2e1 - t507 * t683 / 0.2e1 + t584 * t861 + t146 * t858 + t222 * t859 + t576 * t860 + t106 * t671 + t107 * t475 + t578 / 0.2e1 + t921 * (t163 * t844 + t164 * t846) + t598 * t627 + t604 * t628 + t596 * t629 + (t129 * t576 + t219 * t146 + (t107 + t76) * t911 + (t106 - t81) * t269) * t875;
t569 = t810 / 0.2e1 - t811 / 0.2e1 + t812 / 0.2e1 + t814 / 0.2e1 - t817 / 0.2e1 - t820 / 0.2e1 + t807 / 0.2e1 + t808 / 0.2e1 + t813 / 0.2e1 - t884;
t557 = t387 * t634 + t392 * t636 + t809 / 0.2e1 + t527 * t909 + t534 * t853 + t431 * mrSges(4,1) / 0.2e1 - t432 * mrSges(4,2) / 0.2e1 + t569 + (t152 * t834 - t619 * t836) * t873 / 0.2e1 + (t527 * t77 + t530 * t89) * t875 + (t135 * t530 + t136 * t534) * t878 + Ifges(4,5) * t626 + Ifges(4,6) * t625 + t901 * t842;
t2 = -t553 + t557 + t127 * t767 / 0.2e1 + t777 * t862 + t81 * t671 + t124 * t664 + t76 * t665 + t269 * t667 - t953;
t615 = -t2 * qJD(1) - t10 * qJD(2);
t11 = t559 + (-mrSges(5,3) * t899 + mrSges(7,3) * t911) * t501 - t290 - t941 * t334 + t910 * t262;
t560 = t911 * t670 - t81 * t824 / 0.2e1 + t124 * t762 / 0.2e1 + t76 * t475 + t127 * t477 + t269 * t666 + t953;
t554 = t560 + (t203 * t334 + t273 * t299 - (-t124 + t141) * t380 + (t127 + t140) * t899) * t878 + (t140 * t872 + t93 * mrSges(7,3) / 0.2e1) * t500 + t273 * t859 + t174 * t858 + t334 * t861 + t262 * t860 + (mrSges(7,3) * t869 + (t872 - mrSges(5,3) / 0.2e1) * t141) * t501 + (t129 * t262 + t174 * t219 + (-t81 + t94) * t269 + (t93 + t76) * t911) * t875 - (-t385 / 0.2e1 - t383 / 0.2e1) * t380 + (t855 + t854) * t899 + t954;
t566 = (-pkin(4) * t136 + qJ(5) * t135) * t879 + (qJ(5) * t89 - t77 * t867) * t876 + pkin(4) * t853 + t867 * t909;
t4 = t554 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(7,3) / 0.2e1) * t550 + (-t386 / 0.2e1 - t382 / 0.2e1) * qJ(5) - (-Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t461 - (-Ifges(5,5) / 0.2e1 + Ifges(7,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * t698 + t566 + t884;
t614 = t4 * qJD(1) + t11 * qJD(2);
t562 = (t861 - t278 / 0.2e1) * t501 + (t859 - t339 / 0.2e1) * t459 + (t203 * t501 + t299 * t459 - t551 * t899) * t878 + (-t129 * t501 - t219 * t459 - t551 * t911) * t875;
t606 = t136 * t879 + t77 * t876;
t15 = (t930 + t931) * t550 + t562 + t606 + t906 * t698;
t40 = (-t910 - t941) * t501;
t613 = qJD(1) * t15 + qJD(2) * t40;
t563 = (t459 * t843 - t460 * t845) * mrSges(7,3) + (t269 * t459 - t460 * t911 + t500 * t81 + t501 * t76) * t875 + t384 * t845 + t388 * t843;
t580 = t130 * t875 + t772 / 0.2e1 - t770 / 0.2e1;
t19 = t563 - t580;
t51 = -m(7) * (t269 * t501 + t500 * t911) + (-t500 ^ 2 - t501 ^ 2) * mrSges(7,3);
t612 = qJD(1) * t19 - qJD(2) * t51;
t37 = (t729 / 0.4e1 - t722 / 0.4e1 - t146 / 0.4e1) * t882 - t275;
t591 = (t500 * t530 + t501 * t527) * t875;
t50 = t591 - t571 / 0.2e1 - t336;
t611 = qJD(1) * t37 + qJD(2) * t50;
t43 = (-t434 / 0.4e1 - t697 / 0.4e1 - t174 / 0.4e1) * t882 - t275;
t58 = (t748 / 0.4e1 - t708 / 0.4e1 - t262 / 0.4e1) * t882 - t336;
t610 = qJD(1) * t43 + qJD(2) * t58;
t234 = m(7) * t459;
t301 = m(7) * t501;
t607 = qJD(1) * t234 + qJD(2) * t301;
t605 = t106 * t876 + t163 * t879;
t603 = qJD(4) * (qJ(5) * t906 + Ifges(7,6));
t600 = qJD(4) * (Ifges(7,5) - t752 + t866);
t590 = m(7) * t869 + t141 * t878 + t894;
t100 = -t692 + (-m(6) * (t534 * t834 + t653) - m(7) * (t527 * t834 + t653) - t885) * pkin(3);
t556 = (-(t682 - t530) * t380 + (t534 + t685) * t899) * t878 + ((t269 * t834 + t836 * t911) * pkin(3) + t955) * t875 + t641 * t872 + t636 * t764 + t635 * t762 + t634 * t824 + t896 * t501 + t895 * t500 - t965;
t558 = -(t876 * t911 + t665) * t867 + (-t879 * t899 + t477) * pkin(4) + (-t269 * t876 + t380 * t879 + t664 + t671) * qJ(5) + t965;
t13 = t556 + t558;
t555 = (t530 * t140 + t534 * t141 + (t124 * t836 + t127 * t834) * pkin(3)) * t878 + (t527 * t94 + t530 * t93 + (t76 * t834 + t81 * t836) * pkin(3)) * t875 - t802 / 0.2e1 + t801 / 0.2e1 - t800 / 0.2e1 - t799 / 0.2e1 + t827 / 0.2e1 - t826 / 0.2e1 + t389 * t635 - t621 / 0.2e1 + t620 / 0.2e1 + (t388 + t390) * t634 - t895 * t460 + t896 * t459 + (t384 + t902) * t636;
t565 = (-pkin(4) * t163 + qJ(5) * t164) * t878 + (qJ(5) * t107 - t106 * t867) * t875 - t806 / 0.2e1 + t805 / 0.2e1;
t9 = -t555 + pkin(4) * t668 + t749 * t871 + qJ(5) * t670 - t794 / 0.2e1 + t793 / 0.2e1 - t795 / 0.2e1 - t796 / 0.2e1 - t460 * t659 + t565;
t589 = -t9 * qJD(1) + t13 * qJD(2) - t100 * qJD(3);
t570 = (-qJ(5) - t530) * t551 + t141;
t561 = t570 * t878 + (-t433 + t570) * t875 - t756 - t755 + t894;
t24 = t667 + t668 + t561 + t605;
t398 = 0.4e1 * t530 * t891 + t825;
t588 = -qJD(1) * t24 - qJD(3) * t398;
t579 = -t340 + t632;
t25 = -t460 * t680 + t590 - t886;
t572 = qJD(1) * t25 - qJD(3) * t521 - t940;
t542 = qJ(5) * t685;
t315 = t825 + t868 * t634 + t891 * (0.2e1 * t682 + 0.4e1 * qJ(5));
t49 = t591 + t571 / 0.2e1;
t45 = t833 + t935;
t42 = (-t697 - t434 + t174) * t875;
t36 = (-t722 + t729 + t146) * t875;
t35 = t833 / 0.2e1 + t899 * t878 + t935;
t26 = t590 + t894 + t886;
t23 = t561 - t605 + t894;
t18 = t563 + t580;
t14 = t769 / 0.2e1 - t757 / 0.2e1 - t771 / 0.2e1 - t758 / 0.2e1 + t680 * t698 + t562 - t606;
t12 = t556 - t558 + t579;
t8 = t555 + (qJ(5) * t680 + Ifges(7,6)) * t459 - (t866 / 0.2e1 + Ifges(7,5) + t659) * t460 + (-mrSges(6,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t163 + (-mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t164 + t565 + t633;
t3 = t554 - t566 + t569 + t901 * qJ(5) / 0.2e1;
t1 = t553 + t557 + t560;
t16 = [qJD(2) * t5 + qJD(3) * t6 + qJD(4) * t7 + qJD(5) * t22 + qJD(6) * t31, t1 * qJD(3) + t3 * qJD(4) + t14 * qJD(5) + t18 * qJD(6) + t768 + ((-m(4) * pkin(2) - mrSges(3,1) + t592) * t548 - t941 * t204 + t914 * t848 + t916 * t845 + (m(5) * t152 + m(6) * t135 + t382 + t387) * t899 - t619 * t761 + (-m(5) * t619 - m(6) * t136 + t392 - t393) * t380 + t911 * t386 + m(7) * (t130 * t219 + t269 * t77 + t89 * t911) + (m(4) * pkin(8) + mrSges(4,3)) * (-t835 * t431 + t837 * t432) + t837 * t453 / 0.2e1 - t152 * t765 - t135 * t767 - t136 * t762 - t510 * t683 + t348 * t849 + (-t500 * t950 - t501 * t951 + t594) * t840 + t243 * t846 + t575 * t838 + t574 * t839 + (-Ifges(7,5) * t501 + Ifges(7,6) * t500) * t841 + t77 * t824 + mrSges(3,2) * t547 + t89 * t764 + t508 * t686 + t455 * t688 + Ifges(3,5) * t551 - Ifges(3,6) * t550 + t535 * t282 - pkin(2) * t484 + t269 * t391 + t130 * t339 + t299 * t280 + t219 * t281 + t888 * t847 + t889 * t844 + (t617 + t922) * t512) * qJD(2), t751 + t1 * qJD(2) + (-t616 + m(7) * (t106 * t527 + t107 * t530) + m(6) * (t163 * t534 + t164 * t530) + t620 - t621 + t633 - t429 * mrSges(4,2) - t430 * mrSges(4,1) - t794 + t793 - t795 - t796 - t806 + t805 + (-t163 * t836 + t164 * t834) * t873 - Ifges(4,5) * t657 - Ifges(4,6) * t658 + t459 * t887 - t892 * t460) * qJD(3) + t8 * qJD(4) + t23 * qJD(5) + t36 * qJD(6), t750 + t3 * qJD(2) + t8 * qJD(3) + (t633 - t799 - t800 + t801 - t802 - t826 + t827) * qJD(4) + t26 * qJD(5) + t42 * qJD(6) + ((-pkin(4) * t141 + qJ(5) * t140) * t878 + (qJ(5) * t93 - t867 * t94) * t875) * t880 - t460 * t600 + t459 * t603, qJD(2) * t14 + qJD(3) * t23 + qJD(4) * t26 + t747, qJD(2) * t18 + qJD(3) * t36 + qJD(4) * t42 + t746; -qJD(3) * t2 + qJD(4) * t4 + qJD(5) * t15 + qJD(6) * t19 - t768, -qJD(3) * t10 + qJD(4) * t11 + qJD(5) * t40 - qJD(6) * t51 (m(7) * t955 + m(6) * (t380 * t530 + t534 * t899) + t579 + (t380 * t834 - t836 * t899) * t873 + mrSges(5,3) * t641 + t595 + t892 * t500 + t592 * pkin(8) + (t643 + t887) * t501 + t966) * qJD(3) + t12 * qJD(4) + t35 * qJD(5) + t49 * qJD(6) + t615, t12 * qJD(3) + (t632 + t966) * qJD(4) + t45 * qJD(5) + ((-pkin(4) * t899 + qJ(5) * t380) * t878 + (-qJ(5) * t269 - t867 * t911) * t875) * t880 + t501 * t603 + t500 * t600 + t614, qJD(3) * t35 + qJD(4) * t45 + t613, qJD(3) * t49 + t612; qJD(2) * t2 - qJD(4) * t9 + qJD(5) * t24 + qJD(6) * t37 - t751, qJD(4) * t13 + qJD(6) * t50 - t615, -qJD(4) * t100 + qJD(5) * t398, t315 * qJD(5) + ((-pkin(4) * t682 + t542) * t878 + (-t682 * t867 + t542) * t875) * t880 + t589 + (pkin(3) * t885 + t692) * qJD(4), qJD(4) * t315 - t588, t611; -qJD(2) * t4 + qJD(3) * t9 - qJD(5) * t25 + qJD(6) * t43 - t750, -qJD(3) * t13 + qJD(6) * t58 - t614, -t589 + t939, t939, -t572, t610; -qJD(2) * t15 - qJD(3) * t24 + qJD(4) * t25 + qJD(6) * t234 - t747, qJD(6) * t301 - t613, t588 - t940, t572, 0, t607; -qJD(2) * t19 - qJD(3) * t37 - qJD(4) * t43 - qJD(5) * t234 - t746, -qJD(3) * t50 - qJD(4) * t58 - qJD(5) * t301 - t612, -t611, -t610, -t607, 0;];
Cq  = t16;
