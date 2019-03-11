% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRP9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP9_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP9_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:58
% EndTime: 2019-03-09 17:22:32
% DurationCPUTime: 18.96s
% Computational Cost: add. (19725->897), mult. (39981->1126), div. (0->0), fcn. (37276->6), ass. (0->457)
t920 = Ifges(5,4) + Ifges(4,5);
t898 = -Ifges(7,4) - Ifges(6,5);
t789 = sin(qJ(3));
t505 = t789 * qJ(4);
t520 = cos(qJ(3));
t569 = t520 * pkin(3) + t505;
t465 = -pkin(2) - t569;
t432 = t520 * pkin(4) - t465;
t518 = sin(qJ(5));
t790 = cos(qJ(5));
t663 = t790 * t520;
t449 = t518 * t789 + t663;
t618 = t790 * t789;
t450 = t518 * t520 - t618;
t591 = pkin(5) * t449 + qJ(6) * t450;
t200 = t591 + t432;
t252 = mrSges(7,1) * t449 + mrSges(7,3) * t450;
t892 = m(7) * t200 + t252;
t253 = mrSges(6,1) * t449 - mrSges(6,2) * t450;
t905 = m(6) * t432 + t253;
t682 = Ifges(4,4) * t789;
t582 = Ifges(4,1) * t520 - t682;
t511 = Ifges(5,5) * t789;
t886 = Ifges(5,1) * t520 + t511;
t924 = t886 + t582;
t507 = t520 * qJ(4);
t519 = sin(qJ(2));
t491 = t519 * t507;
t693 = t789 * pkin(3);
t564 = -pkin(4) * t789 - t693;
t555 = -pkin(7) + t564;
t298 = t519 * t555 + t491;
t709 = t519 * t520;
t392 = t518 * t709 - t519 * t618;
t673 = t519 * t789;
t391 = -t518 * t673 - t519 * t663;
t717 = t391 * qJ(6);
t592 = t392 * pkin(5) + t717;
t115 = t298 + t592;
t521 = cos(qJ(2));
t599 = mrSges(7,1) * t450 - mrSges(7,3) * t449;
t600 = mrSges(7,1) * t391 - mrSges(7,3) * t392;
t601 = mrSges(6,1) * t450 + mrSges(6,2) * t449;
t602 = mrSges(6,1) * t391 + mrSges(6,2) * t392;
t435 = Ifges(7,6) * t450;
t436 = Ifges(6,6) * t450;
t439 = Ifges(6,5) * t449;
t440 = Ifges(7,4) * t449;
t875 = -t439 + t436 - t440 - t435;
t923 = t875 * t521 / 0.4e1 - t599 * t115 / 0.2e1 - t600 * t200 / 0.2e1 - t601 * t298 / 0.2e1 - t602 * t432 / 0.2e1;
t824 = -t391 / 0.2e1;
t477 = (pkin(8) - pkin(9)) * t520;
t692 = t789 * pkin(8);
t605 = -t789 * pkin(9) + t692;
t866 = t518 * t477 - t605 * t790;
t922 = t866 / 0.2e1;
t304 = t790 * t477 + t518 * t605;
t891 = t790 * t304 + t518 * t866;
t901 = mrSges(6,1) + mrSges(7,1);
t900 = mrSges(7,2) + mrSges(6,3);
t899 = -mrSges(7,3) + mrSges(6,2);
t897 = -Ifges(6,6) + Ifges(7,6);
t394 = t521 * t449;
t740 = t394 * mrSges(7,2);
t792 = t521 / 0.2e1;
t358 = Ifges(7,6) * t391;
t359 = Ifges(6,6) * t391;
t362 = Ifges(6,5) * t392;
t363 = Ifges(7,4) * t392;
t874 = -t363 - t358 - t362 + t359;
t917 = t874 * t792;
t706 = t520 * t521;
t393 = t518 * t706 - t521 * t618;
t916 = t898 * t519 + (Ifges(6,1) + Ifges(7,1)) * t394 + (-Ifges(6,4) + Ifges(7,5)) * t393;
t437 = Ifges(7,5) * t449;
t259 = -t450 * Ifges(7,1) + t437;
t441 = Ifges(6,4) * t449;
t261 = -t450 * Ifges(6,1) - t441;
t915 = t261 + t259;
t307 = -mrSges(7,2) * t393 - mrSges(7,3) * t519;
t312 = mrSges(6,2) * t519 - mrSges(6,3) * t393;
t914 = t307 + t312;
t913 = t519 * t920 + t924 * t521;
t512 = Ifges(4,4) * t520;
t472 = Ifges(4,1) * t789 + t512;
t769 = Ifges(5,5) * t520;
t912 = Ifges(5,1) * t789 + t472 - t769;
t578 = Ifges(4,5) * t789 + Ifges(4,6) * t520;
t911 = -Ifges(5,4) * t789 - t578;
t780 = pkin(8) * t521;
t784 = pkin(2) * t519;
t478 = -t780 + t784;
t356 = -pkin(7) * t709 + t789 * t478;
t315 = t519 * qJ(4) + t356;
t619 = -pkin(7) * t789 - pkin(3);
t707 = t520 * t478;
t316 = t519 * t619 - t707;
t909 = t315 * t520 + t789 * t316;
t669 = t790 * qJ(4);
t522 = -pkin(3) - pkin(4);
t710 = t518 * t522;
t464 = t669 + t710;
t453 = -qJ(6) + t464;
t732 = t464 * mrSges(6,3);
t908 = -t453 * mrSges(7,2) - t732;
t463 = -t518 * qJ(4) + t522 * t790;
t462 = pkin(5) - t463;
t733 = t463 * mrSges(6,3);
t907 = t462 * mrSges(7,2) - t733;
t906 = t304 * t824 + t392 * t922;
t544 = t518 * t901 + t899 * t790;
t904 = -t200 * t599 - t432 * t601;
t782 = pkin(8) * t519;
t466 = -pkin(2) * t521 - pkin(1) - t782;
t658 = t789 * t466;
t783 = pkin(7) * t520;
t305 = t658 + (-qJ(4) + t783) * t521;
t627 = pkin(9) * t673;
t237 = t627 + t305;
t668 = t790 * t237;
t670 = t521 * t789;
t699 = pkin(7) * t670 - t520 * t466;
t282 = pkin(9) * t709 - t699;
t517 = t521 * pkin(3);
t225 = pkin(4) * t521 - t282 + t517;
t715 = t518 * t225;
t103 = t668 + t715;
t705 = t521 * qJ(6);
t87 = t103 + t705;
t902 = m(7) * t87;
t561 = m(7) * (-pkin(5) * t518 + qJ(6) * t790);
t896 = -t740 / 0.2e1;
t743 = t391 * mrSges(6,3);
t741 = t392 * mrSges(6,3);
t102 = t225 * t790 - t518 * t237;
t89 = t790 * t103;
t895 = t518 * t102 - t89;
t697 = t790 / 0.2e1;
t617 = t392 * t697;
t360 = Ifges(7,5) * t392;
t196 = -t391 * Ifges(7,1) + t521 * Ifges(7,4) + t360;
t364 = Ifges(6,4) * t392;
t198 = -t391 * Ifges(6,1) + t521 * Ifges(6,5) - t364;
t894 = t198 + t196;
t205 = mrSges(7,1) * t392 + mrSges(7,3) * t391;
t206 = mrSges(6,1) * t392 - mrSges(6,2) * t391;
t893 = t205 + t206;
t306 = t517 + t699;
t890 = t306 - t699;
t309 = -mrSges(6,2) * t521 - t741;
t509 = t521 * mrSges(7,3);
t742 = t392 * mrSges(7,2);
t637 = -t509 + t742;
t889 = t309 - t637;
t310 = t521 * mrSges(6,1) + t743;
t744 = t391 * mrSges(7,2);
t311 = -t521 * mrSges(7,1) - t744;
t888 = t310 - t311;
t361 = Ifges(7,5) * t391;
t598 = Ifges(7,1) * t392 + t361;
t365 = Ifges(6,4) * t391;
t212 = -Ifges(6,1) * t392 + t365;
t557 = t582 * t519;
t887 = t519 * t886 - t521 * t920 + t557;
t702 = t453 * t790 + t462 * t518;
t438 = Ifges(7,5) * t450;
t597 = Ifges(7,1) * t449 + t438;
t442 = Ifges(6,4) * t450;
t260 = -Ifges(6,1) * t449 + t442;
t468 = -t507 + t693;
t469 = -Ifges(5,3) * t520 + t511;
t885 = -Ifges(4,2) * t789 + t512;
t679 = Ifges(4,6) * t789;
t579 = Ifges(4,5) * t520 - t679;
t678 = Ifges(5,6) * t789;
t580 = Ifges(5,4) * t520 + t678;
t884 = t579 + t580;
t495 = Ifges(5,5) * t709;
t371 = -Ifges(5,6) * t521 + Ifges(5,3) * t673 + t495;
t883 = -Ifges(5,1) * t673 + t371 + t495;
t882 = -(t520 ^ 2 + t789 ^ 2) * t782 / 0.2e1;
t684 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t686 = -Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1;
t880 = t684 * t393 + t686 * t394;
t470 = Ifges(4,2) * t520 + t682;
t695 = t789 / 0.2e1;
t696 = -t789 / 0.2e1;
t879 = t469 * t695 + t470 * t696;
t467 = -t520 * mrSges(5,1) - t789 * mrSges(5,3);
t570 = -t520 * mrSges(4,1) + t789 * mrSges(4,2);
t254 = -Ifges(7,3) * t450 - t437;
t209 = -Ifges(7,3) * t391 - t360;
t595 = -Ifges(6,2) * t450 + t441;
t596 = -Ifges(6,2) * t391 + t364;
t831 = -t310 / 0.2e1;
t640 = t831 + t311 / 0.2e1;
t641 = t309 / 0.2e1 - t637 / 0.2e1;
t854 = -m(7) / 0.2e1;
t873 = -m(6) / 0.2e1 + t854;
t872 = t885 + t912;
t871 = -t254 + t915;
t870 = -t209 + t894;
t257 = -t449 * Ifges(6,2) - t442;
t869 = -t260 + t597 + t257;
t194 = -t392 * Ifges(6,2) + t521 * Ifges(6,6) - t365;
t868 = -t212 + t598 + t194;
t504 = m(7) * qJ(6) + mrSges(7,3);
t510 = t521 * mrSges(5,3);
t624 = mrSges(5,2) * t673;
t461 = -t510 - t624;
t867 = -m(5) * t305 - t461;
t864 = t440 / 0.2e1 + t439 / 0.2e1 - t436 / 0.2e1 + t435 / 0.2e1;
t863 = t363 / 0.2e1 + t362 / 0.2e1 - t359 / 0.2e1 + t358 / 0.2e1;
t88 = -t521 * pkin(5) - t102;
t862 = m(7) * t88 - t888;
t861 = t889 + t741;
t759 = Ifges(6,3) * t519;
t762 = Ifges(7,6) * t393;
t763 = Ifges(6,6) * t393;
t764 = Ifges(7,2) * t519;
t768 = Ifges(6,5) * t394;
t770 = Ifges(7,4) * t394;
t860 = -t759 / 0.2e1 + t762 / 0.2e1 - t763 / 0.2e1 - t764 / 0.2e1 + t768 / 0.2e1 + t770 / 0.2e1;
t859 = 0.2e1 * m(7);
t857 = m(5) / 0.2e1;
t855 = m(6) / 0.2e1;
t853 = m(7) / 0.2e1;
t852 = pkin(5) / 0.2e1;
t851 = m(5) * pkin(8);
t850 = mrSges(6,1) / 0.2e1;
t849 = -mrSges(6,2) / 0.2e1;
t848 = -mrSges(7,2) / 0.2e1;
t847 = mrSges(7,3) / 0.2e1;
t226 = (-pkin(9) * t521 - t478) * t520 + (-pkin(4) + t619) * t519;
t240 = pkin(9) * t670 + t315;
t105 = t226 * t790 - t518 * t240;
t99 = t519 * pkin(5) - t105;
t846 = -t99 / 0.2e1;
t845 = qJ(6) / 0.2e1;
t844 = -t102 / 0.2e1;
t843 = -t103 / 0.2e1;
t842 = t103 / 0.2e1;
t839 = t205 / 0.2e1;
t838 = t252 / 0.2e1;
t836 = -t304 / 0.2e1;
t835 = -t866 / 0.2e1;
t834 = t304 / 0.2e1;
t729 = t519 * mrSges(7,1);
t314 = t729 + t740;
t829 = t314 / 0.2e1;
t822 = t391 / 0.2e1;
t821 = -t392 / 0.2e1;
t820 = t392 / 0.2e1;
t818 = -t393 / 0.2e1;
t817 = t393 / 0.2e1;
t816 = t394 / 0.2e1;
t810 = -t449 / 0.2e1;
t809 = t449 / 0.2e1;
t807 = -t450 / 0.2e1;
t805 = t450 / 0.2e1;
t804 = t453 / 0.2e1;
t803 = -t462 / 0.2e1;
t802 = t463 / 0.2e1;
t801 = -t464 / 0.2e1;
t800 = t464 / 0.2e1;
t799 = t518 / 0.2e1;
t798 = -t519 / 0.2e1;
t797 = t519 / 0.2e1;
t796 = -t520 / 0.2e1;
t795 = t520 / 0.2e1;
t793 = -t521 / 0.2e1;
t354 = pkin(7) * t706 + t658;
t542 = t627 + t354;
t123 = t282 * t518 - t542 * t790;
t788 = m(7) * t123;
t786 = m(7) * t304;
t785 = m(7) * t518;
t781 = pkin(8) * t520;
t779 = t519 * pkin(7);
t777 = t518 * t88 + t790 * t87;
t771 = Ifges(3,4) * t519;
t767 = Ifges(5,2) * t519;
t761 = Ifges(4,3) * t519;
t756 = t102 * mrSges(6,2);
t755 = t102 * mrSges(7,3);
t754 = t103 * mrSges(6,1);
t753 = t103 * mrSges(7,1);
t752 = t123 * mrSges(6,1);
t751 = t123 * mrSges(7,1);
t124 = t282 * t790 + t518 * t542;
t750 = t124 * mrSges(6,2);
t749 = t124 * mrSges(7,3);
t748 = t304 * mrSges(6,1);
t747 = t304 * mrSges(7,1);
t746 = t866 * mrSges(6,2);
t745 = t866 * mrSges(7,3);
t739 = t449 * mrSges(7,2);
t738 = t449 * mrSges(6,3);
t737 = t450 * mrSges(7,2);
t736 = t450 * mrSges(6,3);
t106 = t518 * t226 + t790 * t240;
t492 = qJ(4) * t706;
t299 = t521 * t555 + t492;
t116 = t393 * pkin(5) - t394 * qJ(6) + t299;
t192 = t521 * Ifges(7,6) + t392 * Ifges(7,3) - t361;
t193 = Ifges(7,5) * t394 - t519 * Ifges(7,6) + Ifges(7,3) * t393;
t195 = Ifges(6,4) * t394 - Ifges(6,2) * t393 - t519 * Ifges(6,6);
t207 = mrSges(7,1) * t393 - mrSges(7,3) * t394;
t208 = mrSges(6,1) * t393 + mrSges(6,2) * t394;
t313 = -mrSges(6,1) * t519 - mrSges(6,3) * t394;
t355 = pkin(7) * t673 + t707;
t577 = Ifges(5,3) * t789 + t769;
t372 = Ifges(5,6) * t519 + t521 * t577;
t373 = -Ifges(4,6) * t521 + t519 * t885;
t374 = Ifges(4,6) * t519 + t521 * t885;
t620 = t693 + pkin(7);
t379 = t519 * t620 - t491;
t380 = t521 * t620 - t492;
t571 = mrSges(5,1) * t789 - t520 * mrSges(5,3);
t429 = t571 * t519;
t430 = t571 * t521;
t572 = mrSges(4,1) * t789 + t520 * mrSges(4,2);
t431 = t572 * t521;
t623 = mrSges(4,3) * t673;
t454 = t521 * mrSges(4,2) - t623;
t455 = -t519 * mrSges(4,2) - mrSges(4,3) * t670;
t688 = mrSges(4,3) * t709;
t456 = -mrSges(4,1) * t521 - t688;
t689 = mrSges(5,2) * t709;
t457 = mrSges(5,1) * t521 + t689;
t458 = mrSges(4,1) * t519 - mrSges(4,3) * t706;
t493 = mrSges(5,2) * t706;
t730 = t519 * mrSges(5,1);
t459 = t493 - t730;
t460 = -mrSges(5,2) * t670 + t519 * mrSges(5,3);
t546 = t572 * t779;
t558 = t521 * t579;
t559 = t521 * t580;
t611 = t670 / 0.2e1;
t612 = -t670 / 0.2e1;
t616 = -t673 / 0.2e1;
t647 = t706 / 0.2e1;
t649 = t709 / 0.2e1;
t708 = t519 * t521;
t94 = -qJ(6) * t519 + t106;
t5 = t916 * t824 + t913 * t649 + (0.2e1 * Ifges(3,4) * t521 - t759 + t762 - t763 - t764 + t768 + t770 + (Ifges(3,1) - Ifges(3,2)) * t519) * t792 + (t559 + t767 + t558 + t761) * t793 + m(5) * (t305 * t315 + t306 * t316 + t379 * t380) + m(6) * (t102 * t105 + t103 * t106 + t298 * t299) + m(7) * (t115 * t116 + t87 * t94 + t88 * t99) + t371 * t611 + t373 * t612 + t374 * t616 + t193 * t820 + t195 * t821 + t192 * t817 + t194 * t818 + t372 * t673 / 0.2e1 - pkin(1) * (mrSges(3,1) * t519 + mrSges(3,2) * t521) + t316 * t457 + t306 * t459 + t305 * t460 + t315 * t461 + t356 * t454 + t354 * t455 + t355 * t456 + t380 * t429 + t379 * t430 + t102 * t313 + t88 * t314 + t87 * t307 + t106 * t309 + t105 * t310 + t99 * t311 + t103 * t312 + t298 * t208 + t299 * t206 + t116 * t205 + t115 * t207 + t521 * t546 + (-t771 + t884 * t519 + (Ifges(3,1) - Ifges(5,2) - Ifges(4,3)) * t521) * t797 + t887 * t647 + t894 * t816 + t431 * t779 + (t771 + t897 * t392 + t898 * t391 + (Ifges(3,2) + Ifges(7,2) + Ifges(6,3)) * t521) * t798 - t94 * t637 + m(4) * (pkin(7) ^ 2 * t708 + t354 * t356 - t355 * t699) - t699 * t458;
t731 = t5 * qJD(1);
t726 = t520 * mrSges(5,2);
t202 = -t391 * pkin(5) + t392 * qJ(6);
t333 = (t520 * t522 - t505) * t519;
t125 = -t202 + t333;
t428 = t569 * t519;
t494 = Ifges(5,6) * t709;
t560 = t519 * t467;
t650 = -t709 / 0.2e1;
t6 = (m(5) * t306 - t456 + t457 - t688) * t354 + (m(6) * t333 + t602) * t298 + (-m(6) * t102 + t862) * t123 + (m(7) * t125 + t600) * t115 - t917 + t373 * t650 + ((-t469 / 0.2e1 + t470 / 0.2e1) * t789 - pkin(7) * t570 - t472 * t795) * t519 ^ 2 + t596 * t821 + t192 * t822 - t103 * t743 - t87 * t744 - t102 * t741 + t578 * t708 / 0.2e1 - t305 * t689 + (m(5) * t428 - t560) * t379 + t428 * t429 + t333 * t206 + t125 * t205 + t868 * t824 + t870 * t820 + t883 * t649 + t887 * t616 + t88 * t742 + (-Ifges(5,4) * t673 + t494) * t793 + (m(6) * t103 + t889 + t902) * t124 - t306 * t624 - (t454 + t623 - t867) * t699;
t725 = t6 * qJD(1);
t7 = t202 * t205 - t298 * t602 + (t209 / 0.2e1 - t198 / 0.2e1 + t596 / 0.2e1 - t88 * mrSges(7,2) - t196 / 0.2e1) * t392 + (t87 * mrSges(7,2) - t192 / 0.2e1 - t212 / 0.2e1 + t194 / 0.2e1 + t598 / 0.2e1) * t391 + t917 + (m(7) * t202 - t600) * t115 + (t862 + t743) * t103 + (t861 + t902) * t102;
t724 = t7 * qJD(1);
t666 = t790 * t309;
t667 = t790 * t637;
t711 = t518 * t521;
t18 = t888 * t711 + (-m(5) * t379 + m(6) * t298 + m(7) * t115 - t429 + t893) * t709 + (m(6) * t895 - m(7) * t777 - t666 + t667 + t867) * t521;
t723 = qJD(1) * t18;
t32 = m(7) * (t115 * t391 + t521 * t87) - t521 * t637 + t391 * t205;
t722 = qJD(1) * t32;
t698 = pkin(5) * t848;
t691 = qJD(4) * t785;
t690 = qJD(6) * t785;
t685 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t675 = -t737 / 0.2e1;
t671 = t521 * t790;
t662 = t789 * t305;
t660 = t789 * t354;
t651 = t711 / 0.2e1;
t639 = -t463 / 0.2e1 + t803;
t638 = t801 + t804;
t628 = -t463 * t518 + t464 * t790;
t625 = -qJD(1) * t521 + qJD(3);
t613 = t671 / 0.2e1;
t249 = -t450 * pkin(5) + t449 * qJ(6);
t445 = t507 + t564;
t201 = -t249 + t445;
t255 = t449 * Ifges(7,3) - t438;
t586 = t123 * t866 + t124 * t304;
t523 = (-t373 / 0.4e1 + t883 / 0.4e1) * t789 + (t192 / 0.4e1 - t868 / 0.4e1) * t450 + (-t596 / 0.4e1 + t870 / 0.4e1) * t449 + t900 * (t123 * t807 + t124 * t810 + t906) + (-t595 / 0.4e1 + t871 / 0.4e1) * t392 + (t255 / 0.4e1 - t869 / 0.4e1) * t391 - t923 + (t298 * t445 + t333 * t432 + t586) * t855 + (t115 * t201 + t125 * t200 + t586) * t853 + t546 / 0.2e1 + t470 * t650 - (-t102 * t855 + t853 * t88 + t640) * t304 + (t886 + 0.2e1 * t469) * t709 / 0.4e1 + (-t472 + t577) * t673 / 0.4e1 - t558 / 0.4e1 + t201 * t839 + t736 * t843 + t738 * t844 + t125 * t838 + t570 * t784 / 0.2e1 + t88 * t739 / 0.2e1 - (t461 + t454) * t692 / 0.2e1 + t379 * t571 / 0.2e1 + t428 * t467 / 0.2e1 + t468 * t429 / 0.2e1 + t445 * t206 / 0.2e1 + t333 * t253 / 0.2e1 - t465 * t560 / 0.2e1 + (t457 / 0.2e1 - t456 / 0.2e1) * t781 - t872 * t673 / 0.4e1 + (t660 / 0.2e1 - t662 / 0.2e1 + t882) * mrSges(5,2) + t882 * mrSges(4,3) + (t557 + t887) * t520 / 0.4e1 + (t468 * t379 + t465 * t428 + (t520 * t890 + t660 - t662) * pkin(8)) * t857 + t890 * t726 / 0.2e1 + t87 * t675 + (t103 * t855 + t853 * t87 + t641) * t866 - t559 / 0.4e1;
t541 = mrSges(7,1) * t846 + t105 * t850 + t106 * t849 + t847 * t94;
t527 = (-pkin(3) * t316 + qJ(4) * t315) * t857 + (t105 * t463 + t106 * t464) * t855 + (t453 * t94 + t462 * t99) * t853 - pkin(3) * t459 / 0.2e1 + qJ(4) * t460 / 0.2e1 + t315 * mrSges(5,3) / 0.2e1 - t316 * mrSges(5,1) / 0.2e1 + t355 * mrSges(4,1) / 0.2e1 - t356 * mrSges(4,2) / 0.2e1 + t307 * t804 + t462 * t829 + t313 * t802 + t312 * t800 - t541;
t1 = -t523 + t767 / 0.2e1 + t761 / 0.2e1 + Ifges(4,6) * t612 + Ifges(5,6) * t611 + t527 + t920 * t647 - t860;
t10 = t595 * t810 + t577 * t796 + t255 * t805 - pkin(2) * t572 + (m(5) * t468 + t571) * t465 + t468 * t467 + t924 * t695 + t871 * t809 + t869 * t807 + t872 * t795 + t879 - t904 + t905 * t445 + t892 * t201;
t589 = -t1 * qJD(1) + t10 * qJD(2);
t550 = t891 * t521;
t525 = (-t789 * t379 + (-t465 * t519 - t780) * t520) * t857 + (t298 * t789 + t432 * t709 - t550) * t855 + (t115 * t789 + t200 * t709 - t550) * t853 + t429 * t696 + t467 * t650 + t893 * t695 + (t252 + t253) * t649 + t900 * (t449 * t613 + t450 * t651);
t529 = t316 * t857 + (t105 * t790 + t518 * t106) * t855 + (t518 * t94 - t790 * t99) * t853 - t730 / 0.2e1 + t313 * t697 - t790 * t314 / 0.2e1 + t914 * t799;
t11 = -t493 + t525 - t529;
t50 = (-m(5) * t465 - t467 + t892 + t905) * t789;
t588 = qJD(1) * t11 + qJD(2) * t50;
t537 = (t115 * t450 + t200 * t391 + t304 * t521) * t853 + t252 * t822 + t205 * t805;
t583 = m(7) * t846 - t729 / 0.2e1;
t26 = t537 + t583 + 0.2e1 * t896;
t61 = t892 * t450;
t587 = -qJD(1) * t26 - qJD(2) * t61;
t568 = m(7) * (-pkin(5) * t123 + qJ(6) * t124);
t567 = m(7) * (pkin(5) * t304 + qJ(6) * t866);
t563 = -t463 * t899 - t464 * t901;
t556 = (t777 + t895) * t853;
t532 = ((-t453 + t464) * t866 + (t462 + t463) * t304) * t853 + t864;
t16 = t684 * t450 + t686 * t449 - t567 / 0.2e1 + ((t845 + t638) * t450 + (t852 + t639) * t449) * mrSges(7,2) + t532 + t899 * (t922 + t835) + t901 * (t834 + t836);
t52 = m(7) * (t453 * t463 + t462 * t464) - t563;
t533 = (t102 * t453 + t103 * t462 + t463 * t87 + t464 * t88) * t853 + t863;
t9 = t640 * t464 + t641 * t463 - t568 / 0.2e1 + (t733 / 0.2e1 + (t803 + t852) * mrSges(7,2) + t686) * t392 + (t732 / 0.2e1 + (t804 + t845) * mrSges(7,2) + t684) * t391 + t533 + t901 * (t123 / 0.2e1 + t842) + t899 * (t124 / 0.2e1 + t102 / 0.2e1);
t554 = t9 * qJD(1) + t16 * qJD(2) + t52 * qJD(3);
t14 = (t257 / 0.2e1 - t260 / 0.2e1 - t255 / 0.2e1 + t597 / 0.2e1) * t450 + (t595 / 0.2e1 - t261 / 0.2e1 + t254 / 0.2e1 - t259 / 0.2e1) * t449 + t892 * t249 + t904;
t524 = -t641 * t866 + t640 * t304 + (t596 / 0.4e1 + t209 / 0.4e1 - t198 / 0.4e1 - t196 / 0.4e1) * t449 + (-t212 / 0.4e1 + t598 / 0.4e1 + t194 / 0.4e1 - t192 / 0.4e1) * t450 + ((t843 + t87 / 0.2e1) * t450 + (t844 - t88 / 0.2e1) * t449 - t906) * mrSges(7,2) + (t595 / 0.4e1 - t261 / 0.4e1 - t259 / 0.4e1 + t254 / 0.4e1 + mrSges(6,3) * t835) * t392 + (t257 / 0.4e1 - t255 / 0.4e1 - t260 / 0.4e1 + t597 / 0.4e1 + mrSges(6,3) * t834) * t391 + (t115 * t249 + t200 * t202 + (t103 - t87) * t866 + (t88 + t102) * t304) * t853 + t202 * t838 + t249 * t839 + t923;
t531 = (-pkin(5) * t99 + qJ(6) * t94) * t854 + pkin(5) * t829 - qJ(6) * t307 / 0.2e1 - t541;
t4 = t519 * t685 + t524 + t531 + t880;
t553 = t4 * qJD(1) + t14 * qJD(2);
t534 = mrSges(6,2) * t613 - mrSges(7,3) * t671 / 0.2e1 + t901 * t651;
t543 = -t667 / 0.2e1 + t666 / 0.2e1;
t530 = t311 * t799 + t518 * t831 - t534 + t543;
t526 = -t510 + (t658 + (-0.2e1 * qJ(4) + t783) * t521) * t857 + (-t464 * t671 + t89 + (t463 * t521 - t102) * t518) * t855 + (-t521 * t702 + t777) * t853 + t530;
t535 = t900 * (t391 * t799 + t617);
t528 = -m(5) * t354 / 0.2e1 + t535 + t873 * (-t123 * t790 + t518 * t124);
t15 = t526 + t528;
t536 = t873 * t891;
t540 = (t855 + t853) * t891;
t31 = t536 + t540;
t70 = m(5) * qJ(4) + m(6) * t628 + m(7) * t702 + mrSges(5,3) + t544;
t552 = qJD(1) * t15 + qJD(2) * t31 + qJD(3) * t70;
t547 = t521 * t561;
t20 = t556 + t547 / 0.2e1 + t530 + t535;
t538 = (-t628 + t702) * t853;
t57 = -t561 / 0.2e1 + t538 + t544;
t551 = t20 * qJD(1) + t57 * qJD(3);
t390 = -m(7) * t453 + mrSges(7,3);
t539 = -t509 + ((-qJ(6) + t453) * t521 - t103) * t853;
t44 = -t788 / 0.2e1 + t539;
t549 = qJD(1) * t44 + qJD(3) * t390;
t233 = mrSges(7,3) + (t464 / 0.4e1 - t669 / 0.4e1 - t710 / 0.4e1 + t845) * t859;
t34 = t509 + (t668 / 0.4e1 + t705 / 0.2e1 + t715 / 0.4e1 - t103 / 0.4e1) * t859;
t548 = qJD(1) * t34 - qJD(3) * t233 + qJD(5) * t504;
t234 = -mrSges(7,3) + (-0.2e1 * qJ(6) + t464) * t853 + m(7) * t800;
t126 = -t739 + t786;
t104 = t739 - t786 / 0.2e1 + m(7) * t836;
t59 = t561 / 0.2e1 + t538;
t39 = t742 + t788 / 0.2e1 + t539;
t33 = (t103 + 0.2e1 * t705) * t853 + m(7) * t842 - t637;
t28 = (mrSges(5,2) + t851) * t520 - t536 + t540 + t900 * (-t790 * t449 - t518 * t450);
t25 = t896 + t740 / 0.2e1 + t537 - t583;
t19 = -t547 / 0.2e1 + ((mrSges(7,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t391 + t640) * t518 + t556 + t534 + t543 + t900 * t617;
t17 = qJ(6) * t675 + t449 * t698 - t746 / 0.2e1 + t745 / 0.2e1 + t747 / 0.2e1 + t748 / 0.2e1 + t532 + (t849 + t847) * t866 + (mrSges(7,1) / 0.2e1 + t850) * t304 + (t449 * t639 + t450 * t638) * mrSges(7,2) + t567 / 0.2e1 + t864;
t13 = -t624 + t526 - t528;
t12 = t525 + t529;
t8 = t863 + t861 * t802 - t755 / 0.2e1 + t756 / 0.2e1 + t753 / 0.2e1 + t754 / 0.2e1 + (t311 + t743) * t800 + t533 - t752 / 0.2e1 - t750 / 0.2e1 - t751 / 0.2e1 + (t391 * t804 + t392 * t803) * mrSges(7,2) + t717 * t848 + t392 * t698 + t310 * t801 + t568 / 0.2e1 + t749 / 0.2e1;
t3 = t524 - t531 + t860;
t2 = (t678 / 0.2e1 - t679 / 0.2e1 + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t520) * t521 + t523 + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + t685) * t519 + t527 + t880;
t21 = [qJD(2) * t5 + qJD(3) * t6 + qJD(4) * t18 + qJD(5) * t7 + qJD(6) * t32, t2 * qJD(3) + t12 * qJD(4) + t3 * qJD(5) + t25 * qJD(6) + t731 + ((t459 - t458) * t692 + t916 * t807 + t912 * t647 + t913 * t695 + (m(6) * t106 + m(7) * t94 + t914) * t304 + t915 * t816 + m(5) * (pkin(8) * t909 + t465 * t380) + t909 * mrSges(5,2) + (-Ifges(5,6) * t520 - t911) * t797 + t905 * t299 + (t449 * t897 + t450 * t898) * t798 + (m(4) * pkin(8) + mrSges(4,3)) * (-t789 * t355 + t356 * t520) + t892 * t116 + (t460 + t455) * t781 - t99 * t737 - t106 * t738 - t94 * t739 + t193 * t809 + t195 * t810 + t255 * t817 + t257 * t818 + t372 * t796 - Ifges(3,6) * t519 + t465 * t430 + t380 * t467 - pkin(2) * t431 + t432 * t208 + ((-m(4) * pkin(2) - mrSges(3,1) + t570) * pkin(7) + Ifges(3,5) + t879) * t521 + t200 * t207 + t105 * t736 + mrSges(3,2) * t779 + t374 * t795 + (-m(6) * t105 + m(7) * t99 - t313 + t314) * t866) * qJD(2), t2 * qJD(2) + t13 * qJD(4) + t8 * qJD(5) + t39 * qJD(6) + t725 + (-t354 * mrSges(4,1) - t354 * mrSges(5,1) + t699 * mrSges(4,2) - t699 * mrSges(5,3) + t494 - t749 + t750 + t751 + t752 + 0.2e1 * (-t123 * t463 + t124 * t464) * t855 + 0.2e1 * (t123 * t462 + t124 * t453) * t853 + 0.2e1 * (-pkin(3) * t354 - qJ(4) * t699) * t857 + (t898 + t907) * t392 + (-t897 + t908) * t391 + (mrSges(5,2) * t468 + t911) * t519) * qJD(3), qJD(2) * t12 + qJD(3) * t13 + qJD(5) * t19 + t723, t724 + t3 * qJD(2) + t8 * qJD(3) + t19 * qJD(4) + (m(7) * (-pkin(5) * t103 + qJ(6) * t102) + t755 - t754 - t753 - t756 + t592 * mrSges(7,2) + t874) * qJD(5) + t33 * qJD(6), qJD(2) * t25 + qJD(3) * t39 + qJD(5) * t33 + t722; -qJD(3) * t1 + qJD(4) * t11 + qJD(5) * t4 + qJD(6) * t26 - t731, qJD(3) * t10 + qJD(4) * t50 + qJD(5) * t14 + qJD(6) * t61, t28 * qJD(4) + t17 * qJD(5) + t104 * qJD(6) + t589 + (-mrSges(5,2) * t505 + m(6) * (t463 * t304 + t464 * t866) + m(7) * (-t304 * t462 + t453 * t866) - pkin(3) * t726 + t746 - t745 - t747 - t748 - t569 * t851 + t875 + t884 + (t570 + t467) * pkin(8) + t908 * t450 + t907 * t449) * qJD(3), qJD(3) * t28 + t588, t17 * qJD(3) + (t591 * mrSges(7,2) + (-m(7) * pkin(5) - t901) * t304 + (mrSges(6,2) - t504) * t866 + t875) * qJD(5) + t126 * qJD(6) + t553, qJD(3) * t104 + qJD(5) * t126 - t587; qJD(2) * t1 + qJD(4) * t15 + qJD(5) * t9 + qJD(6) * t44 - t725, qJD(4) * t31 + qJD(5) * t16 - t589, qJD(4) * t70 + qJD(5) * t52 + qJD(6) * t390, qJD(5) * t59 + t552, t59 * qJD(4) + (m(7) * (-pkin(5) * t464 + qJ(6) * t463) + t563) * qJD(5) + t234 * qJD(6) + t554, qJD(5) * t234 + t549; -qJD(2) * t11 - qJD(3) * t15 + qJD(5) * t20 + t521 * t690 - t723, -qJD(3) * t31 - t588, qJD(5) * t57 - t552 - t690, 0 (-t544 + t561) * qJD(5) + t690 + t551 (qJD(5) - t625) * t785; -qJD(2) * t4 - qJD(3) * t9 - qJD(4) * t20 + qJD(6) * t34 - t724, -qJD(3) * t16 - t553, -qJD(4) * t57 - qJD(6) * t233 - t554, -t551, t504 * qJD(6), t548; -qJD(2) * t26 - qJD(3) * t44 - qJD(5) * t34 - t521 * t691 - t722, t587, qJD(5) * t233 - t549 + t691, t625 * t785, -t548, 0;];
Cq  = t21;
