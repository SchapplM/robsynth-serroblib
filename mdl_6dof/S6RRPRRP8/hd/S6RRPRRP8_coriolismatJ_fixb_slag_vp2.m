% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:22:05
% EndTime: 2019-03-09 12:22:42
% DurationCPUTime: 19.11s
% Computational Cost: add. (36221->827), mult. (77340->1095), div. (0->0), fcn. (84383->8), ass. (0->404)
t518 = sin(pkin(10));
t519 = cos(pkin(10));
t742 = sin(qJ(4));
t744 = cos(qJ(4));
t481 = -t742 * t518 + t744 * t519;
t504 = -pkin(3) * t519 - pkin(2);
t458 = -pkin(4) * t481 + t504;
t520 = sin(qJ(5));
t566 = -t518 * t744 - t519 * t742;
t743 = cos(qJ(5));
t430 = t743 * t481 + t520 * t566;
t557 = t743 * t566;
t840 = t520 * t481 - t557;
t658 = qJ(6) * t840;
t582 = -pkin(5) * t430 - t658;
t232 = t458 + t582;
t860 = t840 * mrSges(7,1);
t881 = -t430 * mrSges(7,3) + t860;
t894 = t232 * t881;
t421 = Ifges(7,5) * t840;
t272 = -Ifges(7,3) * t430 + t421;
t885 = -Ifges(6,4) + Ifges(7,5);
t886 = Ifges(6,1) + Ifges(7,1);
t893 = t430 * t886 + t885 * t840 + t272;
t672 = Ifges(7,5) * t430;
t882 = Ifges(7,3) * t840 + t672;
t424 = Ifges(6,4) * t430;
t883 = -Ifges(6,2) * t840 + t424;
t892 = -t882 / 0.4e1 + t883 / 0.4e1;
t521 = sin(qJ(2));
t514 = t521 * pkin(7);
t651 = t518 * t521;
t489 = pkin(3) * t651 + t514;
t561 = t521 * t566;
t434 = -pkin(4) * t561 + t489;
t465 = t481 * t521;
t374 = t465 * t520 - t521 * t557;
t538 = t465 * t743 + t520 * t561;
t584 = t374 * pkin(5) - qJ(6) * t538;
t163 = t434 + t584;
t845 = t374 * mrSges(6,2);
t856 = t538 * mrSges(6,1) - t845;
t821 = t538 * mrSges(7,1);
t857 = t374 * mrSges(7,3) + t821;
t864 = t430 * mrSges(6,2);
t880 = t840 * mrSges(6,1) + t864;
t891 = t232 * t857 / 0.2e1 + t458 * t856 / 0.2e1 + t163 * t881 / 0.2e1 + t434 * t880 / 0.2e1;
t890 = 0.2e1 * mrSges(6,1);
t889 = 0.2e1 * mrSges(7,3);
t863 = t430 * mrSges(7,2);
t397 = t863 / 0.2e1;
t888 = 0.2e1 * t397;
t799 = 2 * m(7);
t887 = 2 * t799;
t630 = -t863 / 0.2e1;
t829 = Ifges(7,4) + Ifges(6,5);
t363 = Ifges(7,5) * t538;
t522 = cos(qJ(2));
t196 = -t522 * Ifges(7,6) + Ifges(7,3) * t374 + t363;
t229 = -Ifges(7,1) * t374 + t363;
t723 = Ifges(6,4) * t538;
t230 = -Ifges(6,1) * t374 - t723;
t884 = t196 + t229 + t230;
t269 = -mrSges(7,1) * t430 - mrSges(7,3) * t840;
t841 = m(7) * t232 + t269;
t654 = t430 * qJ(6);
t737 = pkin(5) * t840;
t583 = t654 - t737;
t715 = Ifges(7,5) * t374;
t858 = Ifges(7,3) * t538 - t715;
t366 = Ifges(6,4) * t374;
t859 = -Ifges(6,2) * t538 - t366;
t879 = t858 / 0.4e1 - t859 / 0.4e1;
t866 = Ifges(7,6) * t840;
t867 = Ifges(6,6) * t840;
t869 = Ifges(6,5) * t430;
t870 = Ifges(7,4) * t430;
t811 = t869 - t867 + t870 + t866;
t876 = t866 / 0.2e1 - t867 / 0.2e1 + t869 / 0.2e1 + t870 / 0.2e1;
t875 = t430 / 0.2e1;
t874 = -t840 / 0.2e1;
t872 = t840 / 0.2e1;
t767 = t840 / 0.4e1;
t844 = t374 * mrSges(7,2);
t631 = -t844 / 0.2e1;
t676 = t840 * mrSges(7,2);
t659 = qJ(6) * t374;
t738 = pkin(5) * t538;
t220 = t738 + t659;
t824 = Ifges(7,6) * t538;
t825 = Ifges(6,6) * t538;
t846 = Ifges(6,5) * t374;
t847 = Ifges(7,4) * t374;
t810 = -t847 - t846 + t824 - t825;
t854 = -t824 / 0.2e1 + t825 / 0.2e1 + t846 / 0.2e1 + t847 / 0.2e1;
t853 = -mrSges(6,1) / 0.2e1;
t852 = mrSges(6,2) / 0.2e1;
t851 = -mrSges(7,3) / 0.2e1;
t779 = -t374 / 0.2e1;
t850 = t374 / 0.2e1;
t728 = pkin(8) + qJ(3);
t601 = t742 * t728;
t602 = t744 * t728;
t381 = (-pkin(9) * t742 - t601) * t519 + (-pkin(9) * t744 - t602) * t518;
t442 = t481 * t728;
t542 = t481 * pkin(9) + t442;
t806 = t381 * t520 + t542 * t743;
t787 = -t806 / 0.2e1;
t830 = mrSges(6,3) + mrSges(7,2);
t464 = t566 * t522;
t466 = t481 * t522;
t375 = -t464 * t743 + t466 * t520;
t378 = t520 * t464 + t466 * t743;
t842 = t885 * t375 + t378 * t886 + t829 * t521;
t276 = Ifges(7,1) * t840 - t672;
t278 = Ifges(6,1) * t840 + t424;
t818 = t278 + t276;
t492 = t521 * pkin(2) - qJ(3) * t522;
t459 = pkin(7) * t651 + t519 * t492;
t649 = t519 * t521;
t460 = -pkin(7) * t649 + t518 * t492;
t839 = -t459 * t518 + t460 * t519;
t511 = m(7) * qJ(6) + mrSges(7,3);
t837 = qJD(5) * t511;
t836 = t511 * qJD(6);
t749 = t521 / 0.2e1;
t746 = t522 / 0.2e1;
t834 = -t538 / 0.2e1;
t781 = t538 / 0.2e1;
t780 = t538 / 0.4e1;
t739 = pkin(4) * t520;
t634 = t739 / 0.2e1;
t831 = mrSges(7,1) + mrSges(6,1);
t828 = -Ifges(6,6) + Ifges(7,6);
t562 = mrSges(5,3) * t566;
t827 = Ifges(5,4) * t566;
t820 = t538 * mrSges(7,2);
t200 = Ifges(7,1) * t538 - t522 * Ifges(7,4) + t715;
t202 = Ifges(6,1) * t538 - t522 * Ifges(6,5) - t366;
t819 = t202 + t200;
t512 = t522 * mrSges(7,3);
t314 = -t512 - t844;
t690 = t374 * mrSges(6,3);
t315 = mrSges(6,2) * t522 - t690;
t817 = t314 + t315;
t477 = Ifges(5,4) * t481;
t438 = -Ifges(5,1) * t566 + t477;
t816 = Ifges(5,2) * t566 + t438 + t477;
t559 = Ifges(5,4) * t561;
t537 = Ifges(5,1) * t465 - Ifges(5,5) * t522 + t559;
t815 = -Ifges(5,2) * t465 + t537 + t559;
t809 = -t314 / 0.2e1 - t315 / 0.2e1;
t808 = Ifges(5,5) * t481 + Ifges(5,6) * t566 + t811;
t807 = Ifges(5,5) * t561 - Ifges(5,6) * t465 + t810;
t234 = -t743 * t381 + t520 * t542;
t644 = t522 * qJ(6);
t491 = -pkin(2) * t522 - t521 * qJ(3) - pkin(1);
t480 = t519 * t491;
t432 = -pkin(8) * t649 + t480 + (-pkin(7) * t518 - pkin(3)) * t522;
t648 = t519 * t522;
t450 = pkin(7) * t648 + t518 * t491;
t439 = -pkin(8) * t651 + t450;
t282 = t432 * t742 + t439 * t744;
t242 = pkin(9) * t561 + t282;
t628 = t743 * t242;
t281 = t744 * t432 - t439 * t742;
t241 = -t465 * pkin(9) + t281;
t218 = -t522 * pkin(4) + t241;
t647 = t520 * t218;
t98 = t628 + t647;
t85 = t98 - t644;
t804 = m(7) * t85 + t817;
t681 = t538 * mrSges(6,3);
t317 = -mrSges(6,1) * t522 - t681;
t318 = mrSges(7,1) * t522 + t820;
t646 = t520 * t242;
t97 = t218 * t743 - t646;
t86 = t522 * pkin(5) - t97;
t803 = m(7) * t86 - t317 + t318;
t433 = t521 * pkin(3) - pkin(8) * t648 + t459;
t650 = t518 * t522;
t440 = -pkin(8) * t650 + t460;
t283 = t744 * t433 - t440 * t742;
t237 = t521 * pkin(4) - t466 * pkin(9) + t283;
t284 = t742 * t433 + t744 * t440;
t243 = pkin(9) * t464 + t284;
t105 = t237 * t743 - t520 * t243;
t106 = t520 * t237 + t743 * t243;
t94 = qJ(6) * t521 + t106;
t95 = -t521 * pkin(5) - t105;
t802 = t106 * t852 + t105 * t853 + t95 * mrSges(7,1) / 0.2e1 + t94 * t851;
t801 = m(6) * t98 + t804;
t800 = m(6) * t97 - t803;
t798 = m(4) / 0.2e1;
t797 = m(5) / 0.2e1;
t796 = -m(6) / 0.2e1;
t795 = m(6) / 0.2e1;
t794 = -m(7) / 0.2e1;
t793 = m(7) / 0.2e1;
t792 = m(6) * pkin(4);
t791 = m(7) * pkin(4);
t790 = t98 / 0.2e1;
t789 = t220 / 0.2e1;
t788 = t806 / 0.2e1;
t786 = -t583 / 0.2e1;
t662 = t521 * mrSges(7,1);
t679 = t378 * mrSges(7,2);
t320 = -t662 + t679;
t783 = t320 / 0.2e1;
t777 = -t375 / 0.2e1;
t776 = t375 / 0.2e1;
t770 = t378 / 0.2e1;
t765 = -t430 / 0.2e1;
t665 = t465 * mrSges(5,3);
t446 = -mrSges(5,1) * t522 - t665;
t759 = t446 / 0.2e1;
t726 = Ifges(4,4) * t518;
t600 = Ifges(4,1) * t519 - t726;
t758 = Ifges(4,5) * t749 + t600 * t746;
t757 = t464 / 0.2e1;
t756 = t465 / 0.2e1;
t755 = t466 / 0.2e1;
t753 = t481 / 0.2e1;
t503 = qJ(6) + t739;
t752 = -t503 / 0.2e1;
t751 = -t518 / 0.2e1;
t750 = t519 / 0.2e1;
t748 = -t522 / 0.2e1;
t741 = m(7) * t806;
t740 = pkin(4) * t465;
t515 = t522 * pkin(7);
t734 = t97 * mrSges(6,2);
t733 = t97 * mrSges(6,3);
t732 = t97 * mrSges(7,3);
t731 = t98 * mrSges(6,1);
t730 = t98 * mrSges(7,1);
t729 = t98 * mrSges(6,3);
t727 = Ifges(3,4) * t521;
t725 = Ifges(4,4) * t519;
t724 = Ifges(5,4) * t465;
t721 = Ifges(7,4) * t378;
t719 = Ifges(5,5) * t466;
t717 = Ifges(6,5) * t378;
t714 = Ifges(7,2) * t521;
t713 = Ifges(5,6) * t464;
t711 = Ifges(6,6) * t375;
t708 = Ifges(7,6) * t375;
t706 = Ifges(5,3) * t521;
t705 = Ifges(6,3) * t521;
t110 = t241 * t520 + t628;
t702 = t110 * mrSges(6,1);
t701 = t110 * mrSges(7,1);
t111 = t241 * t743 - t646;
t700 = t111 * mrSges(6,2);
t699 = t111 * mrSges(7,3);
t698 = t806 * mrSges(6,1);
t697 = t806 * mrSges(7,1);
t696 = t234 * mrSges(6,2);
t695 = t234 * mrSges(7,3);
t688 = t375 * mrSges(6,1);
t687 = t375 * mrSges(7,1);
t680 = t378 * mrSges(6,2);
t678 = t378 * mrSges(7,3);
t667 = t840 * Ifges(6,4);
t666 = t464 * mrSges(5,1);
t664 = t466 * mrSges(5,2);
t490 = pkin(3) * t650 + t515;
t435 = -pkin(4) * t464 + t490;
t164 = pkin(5) * t375 - qJ(6) * t378 + t435;
t197 = Ifges(7,5) * t378 + t521 * Ifges(7,6) + Ifges(7,3) * t375;
t198 = -Ifges(6,2) * t374 - t522 * Ifges(6,6) + t723;
t199 = Ifges(6,4) * t378 - Ifges(6,2) * t375 + t521 * Ifges(6,6);
t223 = mrSges(7,1) * t374 - mrSges(7,3) * t538;
t224 = mrSges(6,1) * t374 + mrSges(6,2) * t538;
t225 = -t678 + t687;
t226 = t680 + t688;
t313 = -mrSges(7,2) * t375 + mrSges(7,3) * t521;
t316 = -mrSges(6,2) * t521 - mrSges(6,3) * t375;
t319 = mrSges(6,1) * t521 - mrSges(6,3) * t378;
t370 = Ifges(5,4) * t466 + Ifges(5,2) * t464 + t521 * Ifges(5,6);
t371 = t466 * Ifges(5,1) + t464 * Ifges(5,4) + t521 * Ifges(5,5);
t382 = t664 - t666;
t558 = mrSges(5,3) * t561;
t444 = t522 * mrSges(5,2) + t558;
t445 = -mrSges(5,2) * t521 + mrSges(5,3) * t464;
t447 = mrSges(5,1) * t521 - mrSges(5,3) * t466;
t449 = -pkin(7) * t650 + t480;
t595 = -Ifges(4,2) * t518 + t725;
t462 = Ifges(4,6) * t521 + t522 * t595;
t475 = (mrSges(4,1) * t518 + mrSges(4,2) * t519) * t522;
t484 = t522 * mrSges(4,2) - mrSges(4,3) * t651;
t485 = -t521 * mrSges(4,2) - mrSges(4,3) * t650;
t486 = -mrSges(4,1) * t522 - mrSges(4,3) * t649;
t487 = t521 * mrSges(4,1) - mrSges(4,3) * t648;
t536 = Ifges(5,2) * t561 - Ifges(5,6) * t522 + t724;
t554 = t561 / 0.2e1;
t590 = Ifges(4,5) * t519 - Ifges(4,6) * t518;
t619 = -t650 / 0.2e1;
t5 = t842 * t781 + (0.2e1 * Ifges(3,4) * t522 + (Ifges(3,1) - Ifges(3,2)) * t521) * t746 + (Ifges(4,3) * t521 + t522 * t590 + t705 + t706 + t708 - t711 + t713 + t714 + t717 + t719 + t721) * t748 + (-Ifges(4,5) * t522 + t521 * t600) * t648 / 0.2e1 - t521 * (Ifges(3,2) * t522 + t727) / 0.2e1 + t819 * t770 + t197 * t850 - pkin(1) * (t521 * mrSges(3,1) + mrSges(3,2) * t522) + t460 * t484 + t450 * t485 + t459 * t486 + t449 * t487 + t489 * t382 + t283 * t446 + t281 * t447 + t284 * t444 + t282 * t445 + t434 * t226 + t435 * t224 + t106 * t315 + t98 * t316 + t105 * t317 + t95 * t318 + t97 * t319 + t86 * t320 + t85 * t313 + t94 * t314 + t164 * t223 + t163 * t225 + m(4) * (pkin(7) ^ 2 * t521 * t522 + t449 * t459 + t450 * t460) + m(5) * (t281 * t283 + t282 * t284 + t489 * t490) + m(6) * (t105 * t97 + t106 * t98 + t434 * t435) + m(7) * (t163 * t164 + t85 * t94 + t86 * t95) + t370 * t554 + t196 * t776 + t198 * t777 + t199 * t779 + t537 * t755 + t371 * t756 + t536 * t757 + t649 * t758 + 0.2e1 * t475 * t514 - t462 * t651 / 0.2e1 + t490 * (-mrSges(5,1) * t561 + t465 * mrSges(5,2)) + (Ifges(5,5) * t465 + Ifges(5,6) * t561 + t521 * t590 - t727 + t829 * t538 + t828 * t374 + (Ifges(3,1) - Ifges(7,2) - Ifges(4,3) - Ifges(5,3) - Ifges(6,3)) * t522) * t749 + (-Ifges(4,6) * t522 + t521 * t595) * t619;
t663 = t5 * qJD(1);
t545 = Ifges(5,1) * t561 - t724;
t567 = t220 + t740;
t613 = t465 * mrSges(5,1) + mrSges(5,2) * t561;
t6 = t86 * t844 - t224 * t740 - t489 * t613 + t858 * t779 + t198 * t781 + t536 * t756 + t538 * t729 - t374 * t733 - t567 * t223 + t85 * t820 - t465 * t545 / 0.2e1 - t815 * t561 / 0.2e1 + (-m(6) * t740 - t856) * t434 + (t446 + t665) * t282 + (t558 - t444) * t281 + (-m(7) * t567 - t857) * t163 - t801 * t111 + t800 * t110 + t807 * t746 + t884 * t834 + (t859 + t819) * t850;
t661 = t6 * qJD(1);
t7 = t434 * t856 + t220 * t223 + (-t85 * mrSges(7,2) - t729 + t230 / 0.2e1 + t229 / 0.2e1 + t196 / 0.2e1 - t198 / 0.2e1) * t538 + (-t86 * mrSges(7,2) + t733 - t202 / 0.2e1 - t200 / 0.2e1 + t858 / 0.2e1 - t859 / 0.2e1) * t374 + (m(7) * t220 + t857) * t163 + t803 * t98 + t804 * t97 + t810 * t748;
t660 = t7 * qJD(1);
t16 = t444 * t561 - t465 * t446 + m(5) * (-t281 * t465 + t282 * t561) - t484 * t651 - t486 * t649 + m(4) * (-t449 * t519 - t450 * t518) * t521 - t801 * t374 - t800 * t538;
t657 = qJD(1) * t16;
t32 = m(7) * (-t163 * t538 - t522 * t85) - t538 * t223 - t522 * t314;
t656 = qJD(1) * t32;
t638 = mrSges(6,3) * t739;
t637 = t792 / 0.2e1;
t636 = t743 * pkin(4);
t633 = -Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1;
t632 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t618 = t430 * t746;
t606 = t234 * t538 - t374 * t806;
t604 = mrSges(6,3) * t636;
t603 = t636 / 0.2e1;
t563 = t566 * pkin(4);
t248 = -t563 - t583;
t270 = -mrSges(6,1) * t430 + mrSges(6,2) * t840;
t274 = Ifges(6,2) * t430 + t667;
t437 = Ifges(5,2) * t481 - t827;
t547 = Ifges(5,1) * t481 + t827;
t552 = -mrSges(5,1) * t566 + t481 * mrSges(5,2);
t10 = t882 * t765 + t894 + t274 * t874 - t270 * t563 + t504 * t552 + t816 * t753 + (-t547 / 0.2e1 + t437 / 0.2e1) * t566 + (-m(6) * t563 + t880) * t458 + t841 * t248 + (t883 + t818) * t875 + t893 * t872;
t441 = -t518 * t602 - t519 * t601;
t549 = mrSges(5,3) * t554;
t575 = t234 * t110 + t111 * t806;
t523 = t318 * t787 + (-t819 / 0.4e1 + t879) * t430 - t884 * t840 / 0.4e1 + (t248 * t163 + t232 * t567 + t806 * t86 + t575) * t794 + (-t806 * t97 + (-t434 * t566 + t458 * t465) * pkin(4) + t575) * t796 - t891 - t270 * t740 / 0.2e1 - t504 * t613 / 0.2e1 + t830 * (t110 * t874 + t111 * t765 + t234 * t850 + t781 * t806) + t733 * t875 + t729 * t872 - t815 * t481 / 0.4e1 - t816 * t561 / 0.4e1 - (t794 * t85 + t796 * t98 + t809) * t234 + t465 * t437 / 0.4e1 + t808 * t522 / 0.4e1 - t248 * t223 / 0.2e1 - t893 * t538 / 0.4e1 + (t818 / 0.4e1 + t892) * t374 + t274 * t780 + t317 * t788 + t198 * t767 - t567 * t269 / 0.2e1 + t85 * t676 / 0.2e1 - t566 * t536 / 0.4e1 + t566 * t545 / 0.4e1 + t224 * t563 / 0.2e1 - t489 * t552 / 0.2e1 + (t665 / 0.2e1 + t759) * t442 + (-t444 / 0.2e1 + t549) * t441 + t86 * t630 - t465 * t547 / 0.4e1;
t505 = -t636 - pkin(5);
t560 = t708 / 0.2e1 - t711 / 0.2e1 + t717 / 0.2e1 + t721 / 0.2e1 + t705 / 0.2e1 + t714 / 0.2e1 - t802;
t527 = (t105 * t743 + t106 * t520) * t637 + t560 + t706 / 0.2e1 + t503 * t313 / 0.2e1 + t505 * t783 + t719 / 0.2e1 + t713 / 0.2e1 + t283 * mrSges(5,1) / 0.2e1 - t284 * mrSges(5,2) / 0.2e1 + (t503 * t94 + t505 * t95) * t793 + t316 * t634 + t319 * t603;
t2 = t527 + t523;
t581 = -t2 * qJD(1) + t10 * qJD(2);
t275 = Ifges(7,1) * t430 + t421;
t277 = Ifges(6,1) * t430 - t667;
t13 = t894 + t458 * t880 + (t275 / 0.2e1 + t272 / 0.2e1 + t277 / 0.2e1 - t274 / 0.2e1) * t840 - (-t276 / 0.2e1 + t882 / 0.2e1 - t278 / 0.2e1 - t883 / 0.2e1) * t430 - t841 * t583;
t525 = t809 * t234 + (-t317 / 0.2e1 + t318 / 0.2e1) * t806 + (-t198 / 0.4e1 + t230 / 0.4e1 + t229 / 0.4e1 + t196 / 0.4e1) * t840 - (-t202 / 0.4e1 - t200 / 0.4e1 + t879) * t430 + (t806 * t834 + t234 * t779 + (t790 - t85 / 0.2e1) * t840 - (-t86 / 0.2e1 - t97 / 0.2e1) * t430) * mrSges(7,2) + (t277 / 0.4e1 + t275 / 0.4e1 + t272 / 0.4e1 - t274 / 0.4e1 + mrSges(6,3) * t787) * t538 + (-t278 / 0.4e1 - t276 / 0.4e1 - t234 * mrSges(6,3) / 0.2e1 - t892) * t374 + (-t163 * t583 + t220 * t232 + (-t85 + t98) * t234 + (t86 + t97) * t806) * t793 + t269 * t789 + t223 * t786 - t811 * t522 / 0.4e1 + t891;
t541 = (-pkin(5) * t95 + qJ(6) * t94) * t794 + pkin(5) * t783 - qJ(6) * t313 / 0.2e1;
t3 = t633 * t378 + t525 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t521 + t632 * t375 + t541 + t802;
t580 = t3 * qJD(1) + t13 * qJD(2);
t524 = (-t449 * t518 + t450 * t519) * t798 + (t281 * t566 + t481 * t282 - t441 * t465 + t442 * t561) * t797 + (t430 * t98 - t840 * t97 + t606) * t795 + (t430 * t85 + t840 * t86 + t606) * t793 + t317 * t874 + t318 * t872 + t444 * t753 + t486 * t751 + t484 * t750 + t566 * t759 - t562 * t756 + t481 * t549 + t817 * t875 + t830 * (-t374 * t875 + t538 * t872);
t529 = -m(5) * t490 / 0.2e1 + t435 * t796 + t164 * t794 - t688 / 0.2e1 - t687 / 0.2e1 - t680 / 0.2e1 + t678 / 0.2e1 + t666 / 0.2e1 - t664 / 0.2e1 - m(4) * t515 / 0.2e1 + mrSges(4,1) * t619 - mrSges(4,2) * t648 / 0.2e1;
t12 = t529 + t524;
t22 = (t481 ^ 2 + t566 ^ 2) * mrSges(5,3) + m(5) * (t441 * t566 + t442 * t481) + (m(4) * qJ(3) + mrSges(4,3)) * (t518 ^ 2 + t519 ^ 2) + (m(7) + m(6)) * (t234 * t840 + t430 * t806) + t830 * (t430 ^ 2 + t840 ^ 2);
t579 = -qJD(1) * t12 - qJD(2) * t22;
t534 = (-t163 * t840 - t232 * t538 - t522 * t806) * t793 + t269 * t834 + t223 * t874;
t572 = t95 * t794 + t662 / 0.2e1;
t24 = (-t618 - t378 / 0.2e1) * mrSges(7,2) + t534 + t572;
t61 = t841 * t840;
t578 = qJD(1) * t24 - qJD(2) * t61;
t532 = (-t374 * t503 + t505 * t538) * t793 + (-t374 * t520 - t538 * t743) * t637;
t540 = t465 * t637 + t567 * t793;
t29 = t532 - t540 - t857 - t856 - t613;
t531 = (t430 * t503 + t505 * t840) * t793 + (t430 * t520 - t743 * t840) * t637;
t539 = t248 * t793 + t563 * t796;
t37 = t531 - t539 - t552 - t881 - t880;
t577 = qJD(1) * t29 + qJD(2) * t37;
t30 = t779 * t889 + t834 * t890 + (-t738 / 0.4e1 - t659 / 0.4e1 - t220 / 0.4e1) * t799 - t821 + t845;
t38 = t875 * t889 + t874 * t890 + (-t737 / 0.4e1 + t654 / 0.4e1 + t583 / 0.4e1) * t799 - t860 - t864;
t576 = qJD(1) * t30 + qJD(2) * t38;
t190 = t780 * t887;
t257 = t767 * t887;
t573 = qJD(1) * t190 + qJD(2) * t257;
t571 = t110 * t794 + t844 / 0.2e1;
t568 = m(7) * (-pkin(5) * t806 - qJ(6) * t234);
t530 = ((-t503 + t739) * t234 + (t505 + t636) * t806) * t793 + t876;
t535 = t505 * t875 + t840 * t752 + (t520 * t872 + t743 * t875) * pkin(4);
t14 = t633 * t430 + t632 * t840 - t568 / 0.2e1 + (pkin(5) * t875 + t658 / 0.2e1 + t535) * mrSges(7,2) + t530 + t831 * (t787 + t788);
t548 = -(-mrSges(7,3) + mrSges(6,2)) * t636 - t739 * t831;
t358 = -(t503 * t743 + t505 * t520) * t791 - t548;
t526 = (t503 * t97 + t505 * t98 + (t520 * t86 + t743 * t85) * pkin(4)) * t793 - t734 / 0.2e1 + t732 / 0.2e1 - t731 / 0.2e1 - t730 / 0.2e1 + t820 * t752 + t505 * t631 + t318 * t634 - (t317 + t681) * t739 / 0.2e1 + (t690 + t817) * t603 - t854;
t528 = (-pkin(5) * t110 + qJ(6) * t111) * t794 + t702 / 0.2e1 + t701 / 0.2e1 + t700 / 0.2e1 - t699 / 0.2e1 + pkin(5) * t631 + qJ(6) * t820 / 0.2e1 + t854;
t8 = t526 + t528;
t565 = t8 * qJD(1) + t14 * qJD(2) - t358 * qJD(4);
t533 = -t512 + ((-qJ(6) - t503) * t522 + t98) * t793 + t631;
t34 = t533 + t571;
t488 = m(7) * t503 + mrSges(7,3);
t52 = (t875 + t765) * mrSges(7,2);
t564 = -qJD(1) * t34 + qJD(2) * t52 - qJD(4) * t488;
t265 = t397 + t630;
t36 = -t512 + (t628 / 0.4e1 - t644 / 0.2e1 + t647 / 0.4e1 - t98 / 0.4e1) * t799;
t550 = qJD(1) * t36 + qJD(2) * t265 + qJD(4) * t511 + t837;
t474 = mrSges(7,3) + (qJ(6) + 0.2e1 * t634) * m(7);
t258 = (t874 + t872) * m(7);
t191 = (t834 + t781) * m(7);
t115 = t888 + t741;
t53 = t741 / 0.2e1 + t888 + m(7) * t788;
t51 = t531 + t539;
t48 = t532 + t540;
t39 = m(7) * t786 + t583 * t793;
t35 = (t98 - 0.2e1 * t644) * t793 + m(7) * t790 + t314;
t33 = t533 - t571;
t31 = m(7) * t789 - t220 * t793;
t23 = -mrSges(7,2) * t618 + t679 / 0.2e1 + t534 - t572;
t15 = t530 + (t852 + t851) * t234 + t568 / 0.2e1 + t696 / 0.2e1 - t695 / 0.2e1 - t697 / 0.2e1 - t698 / 0.2e1 + (t853 - mrSges(7,1) / 0.2e1) * t806 + pkin(5) * t630 + (t535 - t658 / 0.2e1) * mrSges(7,2) + t876;
t11 = -t529 + t524;
t9 = t526 - t528;
t4 = t525 - t541 + t560;
t1 = t527 - t523;
t17 = [qJD(2) * t5 + qJD(3) * t16 - qJD(4) * t6 + qJD(5) * t7 + qJD(6) * t32, t663 + (t839 * mrSges(4,3) + t818 * t770 + (t316 + t313) * t806 + t94 * t863 + t842 * t872 + t199 * t875 + (-t105 * t840 + t106 * t430) * mrSges(6,3) - Ifges(3,6) * t521 + t504 * t382 - pkin(2) * t475 + t458 * t226 + t441 * t447 + t442 * t445 + t435 * t270 + t95 * t676 + t164 * t269 + t232 * t225 + (-t319 + t320) * t234 + t272 * t776 + t274 * t777 + t197 * t765 + (t485 * t519 - t487 * t518) * qJ(3) + t462 * t750 + t370 * t753 + t438 * t755 + t437 * t757 + t518 * t758 + mrSges(3,2) * t514 + t490 * (-t481 * mrSges(5,1) - mrSges(5,2) * t566) - t566 * t371 / 0.2e1 + t283 * t562 + t284 * t481 * mrSges(5,3)) * qJD(2) + t11 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + t23 * qJD(6) + (Ifges(3,5) + (Ifges(4,1) * t518 + t725) * t750 + (Ifges(4,2) * t519 + t726) * t751 + (-mrSges(4,1) * t519 + mrSges(4,2) * t518 - mrSges(3,1)) * pkin(7)) * qJD(2) * t522 + 0.2e1 * ((-pkin(2) * t515 + t839 * qJ(3)) * t798 + (t283 * t441 + t284 * t442 + t490 * t504) * t797 + (-t105 * t234 + t106 * t806 + t435 * t458) * t795 + (t164 * t232 + t234 * t95 + t806 * t94) * t793) * qJD(2) + (Ifges(4,5) * t518 - Ifges(5,5) * t566 + Ifges(4,6) * t519 + Ifges(5,6) * t481 - t430 * t828 + t829 * t840) * qJD(2) * t749, qJD(2) * t11 + qJD(4) * t48 + qJD(5) * t31 + qJD(6) * t191 + t657, -t661 + t1 * qJD(2) + t48 * qJD(3) + (t374 * t604 - t538 * t638 - t505 * t844 + m(7) * (t110 * t505 + t111 * t503) - t503 * t820 + t699 - t700 - t702 - t701 + (-t110 * t743 + t111 * t520) * t792 - t281 * mrSges(5,2) - t282 * mrSges(5,1) + t807) * qJD(4) + t9 * qJD(5) + t33 * qJD(6), t660 + t4 * qJD(2) + t31 * qJD(3) + t9 * qJD(4) + (m(7) * (-pkin(5) * t98 + qJ(6) * t97) + t732 - t734 - t731 - t730 + t584 * mrSges(7,2) + t810) * qJD(5) + t35 * qJD(6), qJD(2) * t23 + qJD(3) * t191 + qJD(4) * t33 + qJD(5) * t35 + t656; qJD(3) * t12 - qJD(4) * t2 + qJD(5) * t3 + qJD(6) * t24 - t663, qJD(3) * t22 + qJD(4) * t10 + qJD(5) * t13 - qJD(6) * t61, qJD(4) * t51 + qJD(5) * t39 + qJD(6) * t258 - t579, t51 * qJD(3) + (-t697 - t695 + m(7) * (-t234 * t503 + t505 * t806) + t696 - t698 + (-t234 * t520 - t743 * t806) * t792 - t441 * mrSges(5,2) - t442 * mrSges(5,1) - t840 * t638 - t430 * t604 + t505 * t863 - t503 * t676 + t808) * qJD(4) + t15 * qJD(5) + t53 * qJD(6) + t581, t39 * qJD(3) + t15 * qJD(4) + (t582 * mrSges(7,2) + (-m(7) * pkin(5) - t831) * t806 + (mrSges(6,2) - t511) * t234 + t811) * qJD(5) + t115 * qJD(6) + t580, qJD(3) * t258 + qJD(4) * t53 + qJD(5) * t115 + t578; -qJD(2) * t12 - qJD(4) * t29 - qJD(5) * t30 - qJD(6) * t190 - t657, -qJD(4) * t37 - qJD(5) * t38 - qJD(6) * t257 + t579, 0, -t577, -t576, -t573; qJD(2) * t2 + qJD(3) * t29 + qJD(5) * t8 + qJD(6) * t34 + t661, qJD(3) * t37 + qJD(5) * t14 - qJD(6) * t52 - t581, t577, -qJD(5) * t358 + qJD(6) * t488 ((-pkin(5) * t520 + qJ(6) * t743) * t791 + t548) * qJD(5) + t474 * qJD(6) + t565, qJD(5) * t474 - t564; -qJD(2) * t3 + qJD(3) * t30 - qJD(4) * t8 + qJD(6) * t36 - t660, qJD(3) * t38 - qJD(4) * t14 + qJD(6) * t265 - t580, t576, -t565 + t836, t836, t550; -qJD(2) * t24 + qJD(3) * t190 - qJD(4) * t34 - qJD(5) * t36 - t656, qJD(3) * t257 + qJD(4) * t52 - qJD(5) * t265 - t578, t573, t564 - t837, -t550, 0;];
Cq  = t17;
