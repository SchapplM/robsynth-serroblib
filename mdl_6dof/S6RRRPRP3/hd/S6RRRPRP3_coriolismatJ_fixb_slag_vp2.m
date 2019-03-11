% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:39:14
% EndTime: 2019-03-09 16:39:41
% DurationCPUTime: 14.52s
% Computational Cost: add. (34945->739), mult. (70467->949), div. (0->0), fcn. (78688->8), ass. (0->413)
t730 = sin(qJ(3));
t731 = sin(qJ(2));
t733 = cos(qJ(3));
t734 = cos(qJ(2));
t472 = t730 * t731 - t733 * t734;
t500 = sin(pkin(10));
t501 = cos(pkin(10));
t729 = sin(qJ(5));
t732 = cos(qJ(5));
t532 = t500 * t732 + t501 * t729;
t353 = t532 * t472;
t468 = t729 * t500 - t732 * t501;
t356 = t468 * t472;
t763 = -mrSges(7,3) / 0.2e1;
t764 = mrSges(6,2) / 0.2e1;
t615 = t764 + t763;
t765 = mrSges(7,1) / 0.2e1;
t616 = t765 + mrSges(6,1) / 0.2e1;
t699 = t501 * mrSges(5,2);
t701 = t500 * mrSges(5,1);
t767 = m(7) / 0.2e1;
t769 = m(6) / 0.2e1;
t771 = m(5) / 0.2e1;
t628 = t731 * pkin(7);
t479 = -pkin(8) * t731 - t628;
t631 = t734 * pkin(7);
t480 = pkin(8) * t734 + t631;
t791 = t730 * t479 + t733 * t480;
t658 = t472 * t500;
t814 = -pkin(4) * t658 + t791;
t793 = -t353 * pkin(5) - t356 * qJ(6) + t814;
t847 = -t616 * t353 + t615 * t356 + (-t699 / 0.2e1 - t701 / 0.2e1) * t472 + t767 * t793 + t769 * t814 + t771 * t791;
t473 = -t730 * t734 - t731 * t733;
t525 = t468 * t473;
t492 = -pkin(2) * t734 - pkin(1);
t381 = t472 * pkin(3) + t473 * qJ(4) + t492;
t243 = t501 * t381 - t500 * t791;
t244 = t500 * t381 + t501 * t791;
t539 = m(5) * (-t500 * t243 + t244 * t501);
t796 = t532 * t473;
t797 = t796 * mrSges(6,3);
t265 = -mrSges(6,2) * t472 + t797;
t703 = t472 * mrSges(7,3);
t270 = mrSges(7,2) * t796 + t703;
t583 = -t270 / 0.2e1 - t265 / 0.2e1;
t269 = -mrSges(7,1) * t472 + mrSges(7,2) * t525;
t798 = mrSges(6,3) * t525;
t268 = mrSges(6,1) * t472 - t798;
t801 = -t268 / 0.2e1;
t584 = t269 / 0.2e1 + t801;
t656 = t473 * t500;
t388 = -mrSges(5,2) * t472 + mrSges(5,3) * t656;
t648 = t501 * t388;
t587 = -t648 / 0.2e1;
t655 = t473 * t501;
t390 = t472 * mrSges(5,1) + mrSges(5,3) * t655;
t650 = t500 * t390;
t589 = t650 / 0.2e1;
t614 = mrSges(7,2) / 0.2e1 + mrSges(6,3) / 0.2e1;
t846 = -(t525 * t614 + t584) * t532 + (t614 * t796 - t583) * t468 - t539 / 0.2e1 + t589 + t587 + t847;
t741 = t532 / 0.2e1;
t742 = -t468 / 0.2e1;
t795 = t265 + t270;
t841 = mrSges(6,3) + mrSges(7,2);
t845 = t539 / 0.2e1 + t532 * t801 + t269 * t741 - t650 / 0.2e1 + t648 / 0.2e1 + t795 * t742 + t841 * (t525 * t741 + t796 * t742) + t847;
t150 = Ifges(6,4) * t356 + Ifges(6,2) * t353 - t473 * Ifges(6,6);
t152 = Ifges(7,1) * t356 - Ifges(7,4) * t473 - Ifges(7,5) * t353;
t154 = Ifges(6,1) * t356 + Ifges(6,4) * t353 - Ifges(6,5) * t473;
t715 = Ifges(5,2) * t500;
t720 = Ifges(5,4) * t501;
t283 = -Ifges(5,6) * t473 + (t715 - t720) * t472;
t721 = Ifges(5,4) * t500;
t284 = -Ifges(5,5) * t473 + (-Ifges(5,1) * t501 + t721) * t472;
t456 = Ifges(7,5) * t532;
t398 = Ifges(7,3) * t468 + t456;
t718 = Ifges(6,4) * t532;
t400 = -Ifges(6,2) * t468 + t718;
t716 = Ifges(7,5) * t468;
t402 = Ifges(7,1) * t532 + t716;
t459 = Ifges(6,4) * t468;
t404 = Ifges(6,1) * t532 - t459;
t657 = t472 * t501;
t735 = t500 / 0.2e1;
t738 = -t473 / 0.2e1;
t753 = t356 / 0.2e1;
t755 = -t353 / 0.2e1;
t756 = t353 / 0.2e1;
t760 = Ifges(7,5) * t753 + Ifges(7,6) * t738 + Ifges(7,3) * t755;
t838 = -Ifges(6,6) + Ifges(7,6);
t840 = Ifges(7,4) + Ifges(6,5);
t520 = t150 * t742 + t284 * t735 + t398 * t755 + t400 * t756 + t468 * t760 + t501 * t283 / 0.2e1 + (Ifges(5,2) * t501 + t721) * t658 / 0.2e1 - (Ifges(5,1) * t500 + t720) * t657 / 0.2e1 + Ifges(4,6) * t473 - Ifges(4,5) * t472 + (t402 + t404) * t753 + (t152 + t154) * t741 + (Ifges(5,5) * t500 + Ifges(5,6) * t501 + t468 * t838 + t532 * t840) * t738;
t395 = mrSges(7,1) * t468 - mrSges(7,3) * t532;
t816 = t793 * t395;
t477 = -mrSges(5,1) * t501 + t500 * mrSges(5,2);
t817 = t791 * t477;
t824 = t791 * mrSges(4,1);
t422 = -t733 * t479 + t480 * t730;
t825 = t422 * mrSges(4,2);
t396 = mrSges(6,1) * t468 + mrSges(6,2) * t532;
t834 = t814 * t396;
t844 = t520 + t816 + t817 + t825 - t824 + t834;
t843 = t816 / 0.2e1 + t825 / 0.2e1 + t834 / 0.2e1;
t768 = -m(7) / 0.2e1;
t702 = t473 * mrSges(7,1);
t711 = t356 * mrSges(7,2);
t267 = t702 + t711;
t842 = t267 / 0.2e1;
t839 = Ifges(7,2) + Ifges(6,3);
t317 = -pkin(4) * t656 + t422;
t837 = t317 * t814;
t630 = t733 * pkin(2);
t491 = -t630 - pkin(3);
t723 = t501 * pkin(4);
t476 = t491 - t723;
t836 = t476 * t814;
t487 = -pkin(3) - t723;
t835 = t487 * t814;
t690 = qJ(6) * t796;
t726 = pkin(5) * t525;
t833 = t690 - t726;
t831 = t817 / 0.2e1 - t824 / 0.2e1;
t567 = t733 * t729;
t568 = t733 * t732;
t431 = (t500 * t568 + t501 * t567) * pkin(2);
t432 = (-t500 * t567 + t501 * t568) * pkin(2);
t830 = t841 * (t431 * t532 - t432 * t468);
t558 = -pkin(5) * t796 - qJ(6) * t525;
t118 = t317 + t558;
t190 = -t353 * mrSges(7,1) - t356 * mrSges(7,3);
t191 = -t353 * mrSges(6,1) + t356 * mrSges(6,2);
t192 = -mrSges(7,1) * t796 - mrSges(7,3) * t525;
t193 = -mrSges(6,1) * t796 + mrSges(6,2) * t525;
t264 = mrSges(6,2) * t473 + mrSges(6,3) * t353;
t266 = -mrSges(6,1) * t473 - mrSges(6,3) * t356;
t271 = mrSges(7,2) * t353 - mrSges(7,3) * t473;
t559 = -t699 - t701;
t378 = t559 * t472;
t379 = t559 * t473;
t387 = mrSges(5,2) * t473 + mrSges(5,3) * t658;
t389 = -t473 * mrSges(5,1) + mrSges(5,3) * t657;
t717 = Ifges(7,5) * t796;
t153 = Ifges(7,1) * t525 + Ifges(7,4) * t472 - t717;
t348 = Ifges(6,4) * t796;
t155 = Ifges(6,1) * t525 + Ifges(6,5) * t472 + t348;
t585 = t155 / 0.2e1 + t153 / 0.2e1;
t345 = Ifges(7,5) * t525;
t149 = Ifges(7,6) * t472 - Ifges(7,3) * t796 + t345;
t719 = Ifges(6,4) * t525;
t151 = Ifges(6,2) * t796 + Ifges(6,6) * t472 + t719;
t586 = t149 / 0.2e1 - t151 / 0.2e1;
t659 = t472 * qJ(6);
t139 = t472 * pkin(4) + pkin(9) * t655 + t243;
t597 = t729 * t139;
t176 = pkin(9) * t656 + t244;
t603 = t732 * t176;
t85 = t603 + t597;
t70 = t85 + t659;
t84 = t139 * t732 - t176 * t729;
t71 = -t472 * pkin(5) - t84;
t829 = t814 * t193 + t793 * t192 + t791 * t379 + t118 * t190 + t243 * t389 + t244 * t387 + t317 * t191 + t422 * t378 + t492 * (-mrSges(4,1) * t473 - mrSges(4,2) * t472) + t70 * t271 + t71 * t267 + t84 * t266 + t85 * t264 - (t760 - t150 / 0.2e1) * t796 - t586 * t353 + (t154 / 0.2e1 + t152 / 0.2e1) * t525 + t585 * t356;
t800 = -t525 / 0.2e1;
t754 = t796 / 0.2e1;
t826 = pkin(3) * t791;
t823 = t118 * t793;
t556 = t468 * pkin(5) - qJ(6) * t532;
t373 = t487 + t556;
t342 = -t630 + t373;
t822 = t342 * t793;
t821 = t373 * t793;
t664 = t422 * t791;
t820 = t491 * t791;
t819 = t500 * t422;
t818 = t501 * t422;
t799 = mrSges(7,1) + mrSges(6,1);
t815 = t799 * t800;
t499 = t501 ^ 2;
t792 = t500 ^ 2 + t499;
t612 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t613 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t813 = -t612 * t353 - t613 * t356;
t645 = t711 / 0.2e1 + t702 / 0.2e1;
t724 = t473 * pkin(5);
t407 = -t473 * pkin(3) + t472 * qJ(4);
t629 = t731 * pkin(2);
t384 = t629 + t407;
t250 = t501 * t384 + t819;
t566 = -t473 * pkin(4) + pkin(9) * t657;
t140 = t250 + t566;
t251 = t500 * t384 - t818;
t636 = pkin(9) * t658;
t185 = t636 + t251;
t88 = t140 * t732 - t185 * t729;
t73 = -t88 + t724;
t809 = t73 * t768 - t645;
t627 = t730 * pkin(2);
t802 = m(6) + m(7);
t808 = (m(5) + t802) * t627 / 0.2e1;
t807 = -qJ(6) * t271 / 0.2e1 + pkin(5) * t842;
t806 = t792 * mrSges(5,3);
t805 = (t71 + t84) * t767 + t584;
t794 = -t268 + t269;
t790 = -t468 * t840 + t532 * t838;
t789 = t525 * t838 + t796 * t840;
t681 = t251 * t501;
t682 = t250 * t500;
t788 = t681 - t682;
t784 = -t532 * t192 / 0.2e1 + t395 * t800;
t495 = m(7) * qJ(6) + mrSges(7,3);
t562 = t627 + qJ(4);
t783 = t792 * t562;
t393 = mrSges(7,1) * t532 + t468 * mrSges(7,3);
t394 = mrSges(6,1) * t532 - t468 * mrSges(6,2);
t689 = qJ(6) * t468;
t725 = pkin(5) * t532;
t557 = -t689 - t725;
t774 = 0.2e1 * m(7);
t112 = (-t725 / 0.4e1 - t689 / 0.4e1 + t557 / 0.4e1) * t774 - t393 - t394;
t40 = (-t726 / 0.4e1 + t690 / 0.4e1 + t833 / 0.4e1) * t774 + 0.2e1 * t815 + 0.2e1 * (mrSges(7,3) - mrSges(6,2)) * t754;
t637 = qJD(2) + qJD(3);
t782 = qJD(1) * t40 + t112 * t637;
t172 = t525 * t774 / 0.2e1;
t375 = m(7) * t532;
t781 = -qJD(1) * t172 - t375 * t637;
t780 = -m(7) * pkin(5) - t799;
t779 = -mrSges(6,2) + t495;
t778 = (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t473 + t807 - t813;
t777 = -t733 * mrSges(4,2) + (-mrSges(4,1) + t395 + t396 + t477) * t730;
t776 = Ifges(6,6) * t756 + Ifges(7,6) * t755 + t738 * t839 + t753 * t840 - t807;
t496 = t501 * pkin(9);
t478 = t501 * qJ(4) + t496;
t580 = (-pkin(9) - qJ(4)) * t500;
t419 = t478 * t729 - t580 * t732;
t420 = t478 * t732 + t580 * t729;
t649 = t501 * t387;
t651 = t500 * t389;
t653 = t487 * t191;
t670 = t373 * t190;
t89 = t729 * t140 + t732 * t185;
t693 = t89 * t468;
t694 = t88 * t532;
t695 = t73 * t532;
t460 = t473 * qJ(6);
t72 = -t460 + t89;
t696 = t72 * t468;
t728 = pkin(3) * t378;
t775 = (qJ(4) * t788 - t826) * t771 + (-t419 * t88 + t420 * t89 + t835) * t769 + (t419 * t73 + t420 * t72 + t821) * t767 - t728 / 0.2e1 + t670 / 0.2e1 + t653 / 0.2e1 + (t649 / 0.2e1 - t651 / 0.2e1) * qJ(4) + (t681 / 0.2e1 - t682 / 0.2e1) * mrSges(5,3) + (-t694 / 0.2e1 - t693 / 0.2e1) * mrSges(6,3) + (t695 / 0.2e1 - t696 / 0.2e1) * mrSges(7,2) + t843;
t770 = -m(6) / 0.2e1;
t766 = -mrSges(6,1) / 0.2e1;
t258 = t501 * t407 + t819;
t164 = t258 + t566;
t259 = t500 * t407 - t818;
t198 = t636 + t259;
t92 = t164 * t732 - t198 * t729;
t81 = -t92 + t724;
t762 = t81 / 0.2e1;
t759 = -t833 / 0.2e1;
t758 = t342 / 0.2e1;
t544 = t501 * t562;
t449 = t544 + t496;
t537 = t500 * (-pkin(9) - t562);
t367 = t449 * t729 - t537 * t732;
t751 = -t367 / 0.2e1;
t368 = t449 * t732 + t537 * t729;
t750 = -t368 / 0.2e1;
t749 = t373 / 0.2e1;
t748 = -t557 / 0.2e1;
t747 = -t419 / 0.2e1;
t746 = -t420 / 0.2e1;
t745 = -t431 / 0.2e1;
t744 = t431 / 0.2e1;
t743 = -t432 / 0.2e1;
t737 = t476 / 0.2e1;
t736 = t487 / 0.2e1;
t722 = -t70 + t85;
t698 = t501 * Ifges(5,5);
t700 = t500 * Ifges(5,6);
t511 = (t499 * Ifges(5,1) / 0.2e1 - Ifges(5,3) - Ifges(4,2) + Ifges(4,1) + (-t720 + t715 / 0.2e1) * t500 - t839) * t473 + (Ifges(4,4) - t698 + t700) * t472 + t813;
t513 = -t501 * t284 / 0.2e1 + t283 * t735 + (t698 / 0.2e1 - t700 / 0.2e1 - Ifges(4,4)) * t473 + t613 * t525 + t612 * t796;
t3 = (mrSges(4,1) * t629 + t511) * t472 - pkin(1) * (mrSges(3,1) * t731 + mrSges(3,2) * t734) + (-mrSges(4,2) * t629 + t513) * t473 + m(4) * t492 * t629 + m(5) * (t243 * t250 + t244 * t251 + t664) + m(6) * (t84 * t88 + t85 * t89 + t837) + m(7) * (t70 * t72 + t71 * t73 + t823) + t251 * t388 + t250 * t390 + t89 * t265 + t88 * t268 + t73 * t269 + t72 * t270 + (Ifges(3,1) - Ifges(3,2)) * t734 * t731 + (-t731 ^ 2 + t734 ^ 2) * Ifges(3,4) + t829;
t714 = t3 * qJD(1);
t712 = t796 * mrSges(7,3);
t710 = t525 * mrSges(6,1);
t709 = t525 * mrSges(7,1);
t704 = t468 * mrSges(7,2);
t93 = t729 * t164 + t732 * t198;
t80 = -t460 + t93;
t6 = t511 * t472 + t513 * t473 + m(5) * (t243 * t258 + t244 * t259 + t664) + m(6) * (t84 * t92 + t85 * t93 + t837) + m(7) * (t70 * t80 + t71 * t81 + t823) + t259 * t388 + t258 * t390 + t93 * t265 + t92 * t268 + t81 * t269 + t80 * t270 + t829;
t697 = t6 * qJD(1);
t188 = t709 - t712;
t189 = t796 * mrSges(6,2) + t710;
t194 = Ifges(7,3) * t525 + t717;
t195 = -Ifges(6,2) * t525 + t348;
t196 = Ifges(7,1) * t796 + t345;
t197 = Ifges(6,1) * t796 - t719;
t9 = t317 * t189 - t833 * t192 + (t196 / 0.2e1 + t197 / 0.2e1 - t70 * mrSges(7,2) + t586) * t525 - (t194 / 0.2e1 - t195 / 0.2e1 - t71 * mrSges(7,2) - t585) * t796 + (-m(7) * t833 + t188) * t118 + (m(7) * t71 + t794 - t798) * t85 + (m(7) * t70 + t795 - t797) * t84 + t789 * t472 / 0.2e1;
t692 = t9 * qJD(1);
t16 = t795 * t796 + t794 * t525 + m(7) * (t525 * t71 + t70 * t796) + m(6) * (-t525 * t84 + t796 * t85) + (m(5) * (t243 * t501 + t244 * t500) + t500 * t388 + t501 * t390) * t473;
t688 = qJD(1) * t16;
t33 = m(7) * (-t118 * t525 + t472 * t70) + t472 * t270 - t525 * t192;
t687 = qJD(1) * t33;
t684 = t118 * t532;
t680 = t258 * t500;
t679 = t259 * t501;
t676 = t342 * t190;
t674 = t367 * t266;
t673 = t367 * t267;
t672 = t368 * t264;
t671 = t368 * t271;
t669 = t419 * t266;
t668 = t419 * t267;
t667 = t420 * t264;
t666 = t420 * t271;
t654 = t476 * t191;
t652 = t491 * t378;
t647 = t112 * qJD(4);
t231 = m(7) * t748 + t557 * t767;
t646 = t231 * qJD(5);
t642 = -t112 * qJD(5) - t375 * qJD(6);
t641 = t375 * qJD(4);
t635 = t81 * t532 * mrSges(7,2);
t634 = t92 * t532 * mrSges(6,3);
t633 = t93 * t468 * mrSges(6,3);
t626 = mrSges(5,3) * t680;
t625 = mrSges(5,3) * t679;
t618 = m(7) * t744;
t445 = -t704 / 0.2e1;
t609 = t704 / 0.2e1;
t604 = t472 * t445 + t784;
t578 = t367 * t431 + t368 * t432;
t575 = t419 * t431 + t420 * t432;
t570 = 0.2e1 * t445;
t560 = m(7) * t762 + t645;
t555 = -t468 * t70 + t532 * t71;
t554 = -t468 * t85 - t532 * t84;
t535 = t387 * t544;
t545 = t500 * t562;
t536 = t389 * t545;
t502 = t80 * t609 - t843 - t831 - t676 / 0.2e1 - t673 / 0.2e1 + t674 / 0.2e1 - t672 / 0.2e1 - t671 / 0.2e1 - t654 / 0.2e1 - t652 / 0.2e1 + t633 / 0.2e1 + t634 / 0.2e1 + (t317 * t627 - t367 * t92 + t368 * t93 - t431 * t84 + t432 * t85 + t836) * t770 + t626 / 0.2e1 + t269 * t745 + t268 * t744 - t625 / 0.2e1 + t536 / 0.2e1 + t589 * t630 + t587 * t630 + t795 * t743 - t535 / 0.2e1 - m(5) * (t422 * t627 + t820 + (t244 * t630 + t259 * t562) * t501 + (-t243 * t630 - t258 * t562) * t500) / 0.2e1 + (t118 * t627 + t367 * t81 + t368 * t80 + t431 * t71 + t432 * t70 + t822) * t768 - (t379 + t193 + t192) * t627 / 0.2e1 - t635 / 0.2e1;
t2 = t502 - t669 / 0.2e1 + t666 / 0.2e1 + t667 / 0.2e1 + t668 / 0.2e1 + t775 + t831;
t533 = t792 * t733;
t31 = (mrSges(5,3) * t533 + t777) * pkin(2) + m(6) * (t476 * t627 + t578) + m(7) * (t342 * t627 + t578) + m(5) * (t491 * t627 + t630 * t783) + t830;
t553 = -t2 * qJD(1) + t31 * qJD(2);
t217 = t342 * t557;
t397 = Ifges(7,3) * t532 - t716;
t399 = -Ifges(6,2) * t532 - t459;
t401 = -Ifges(7,1) * t468 + t456;
t403 = -Ifges(6,1) * t468 - t718;
t516 = -t557 * t395 + (t397 / 0.2e1 - t404 / 0.2e1 - t402 / 0.2e1 - t399 / 0.2e1) * t468 - (-t398 / 0.2e1 - t403 / 0.2e1 - t401 / 0.2e1 + t400 / 0.2e1) * t532;
t20 = -m(7) * t217 + t342 * t393 + t476 * t394 + t516;
t100 = t557 * t118;
t510 = (-t195 / 0.4e1 - t155 / 0.4e1 - t153 / 0.4e1 + t194 / 0.4e1) * t468 - (-t197 / 0.4e1 - t196 / 0.4e1 - t149 / 0.4e1 + t151 / 0.4e1) * t532 + t118 * t393 / 0.2e1 + t395 * t759 + t317 * t394 / 0.2e1 + t192 * t748 + t790 * t472 / 0.4e1;
t523 = (-t71 / 0.2e1 - t84 / 0.2e1) * t468 - (-t85 / 0.2e1 + t70 / 0.2e1) * t532;
t540 = t398 / 0.4e1 - t400 / 0.4e1 + t401 / 0.4e1 + t403 / 0.4e1;
t541 = t397 / 0.4e1 - t399 / 0.4e1 - t402 / 0.4e1 - t404 / 0.4e1;
t504 = (-t342 * t833 + t367 * t722 - t100) * t767 + t510 + (t525 * t750 - t751 * t796 + t523) * mrSges(7,2) + t583 * t367 + t189 * t737 + t188 * t758 + (mrSges(6,3) * t750 + t540) * t525 - (mrSges(6,3) * t751 + t541) * t796 + t805 * t368;
t518 = (-pkin(5) * t73 + qJ(6) * t72) * t768 + t72 * t763 + t73 * t765 + t88 * t766 + t89 * t764;
t4 = t504 + t518 + t778;
t552 = t4 * qJD(1) + t20 * qJD(2);
t550 = t367 * t525 + t368 * t796;
t515 = (t550 + t554) * t770 + (t550 + t555) * t768;
t10 = t515 + t846;
t282 = t368 * t468;
t527 = t806 + t841 * (t468 ^ 2 + t532 ^ 2);
t47 = t527 + t802 * (t367 * t532 - t282) + m(5) * t783;
t551 = qJD(1) * t10 - qJD(2) * t47;
t549 = t419 * t525 + t420 * t796;
t301 = t532 * t395;
t138 = -t342 * t375 - t301;
t521 = (-t342 * t525 + t368 * t472 - t684) * t767 + t604;
t26 = t521 + t809;
t547 = -qJD(1) * t26 - qJD(2) * t138;
t39 = t703 + (t603 / 0.4e1 + t597 / 0.4e1 + t659 / 0.2e1 - t85 / 0.4e1) * t774;
t546 = qJD(1) * t39 + qJD(5) * t495;
t543 = m(7) * (-pkin(5) * t431 + qJ(6) * t432);
t538 = t301 + (t342 + t373) * t375 / 0.2e1;
t245 = t373 * t557;
t506 = (t749 + t758) * t393 + (t736 + t737) * t394 + (-t217 - t245) * t767 + t516;
t14 = t615 * t432 + t616 * t431 - t543 / 0.2e1 + t506;
t21 = -m(7) * t245 + t373 * t393 + t487 * t394 + t516;
t503 = (-t373 * t833 + t419 * t722 - t100) * t767 + t510 + (t525 * t746 - t747 * t796 + t523) * mrSges(7,2) + t583 * t419 + (mrSges(6,3) * t746 + t540) * t525 - (mrSges(6,3) * t747 + t541) * t796 + t189 * t736 + t188 * t749 + t805 * t420;
t517 = (-pkin(5) * t81 + qJ(6) * t80) * t768 + t80 * t763 + mrSges(7,1) * t762 + t92 * t766 + t93 * t764;
t8 = t503 + t517 + t778;
t531 = t8 * qJD(1) + t14 * qJD(2) + t21 * qJD(3);
t514 = (t549 + t554) * t770 + (t549 + t555) * t768;
t12 = t514 + t846;
t338 = t420 * t468;
t509 = t527 + t792 * t771 * (t627 + 0.2e1 * qJ(4)) + (t769 + t767) * (-t282 - t338 - (-t367 - t419) * t532);
t36 = -t509 + t808;
t48 = m(5) * t792 * qJ(4) + t527 + t802 * (t419 * t532 - t338);
t530 = qJD(1) * t12 + qJD(2) * t36 - qJD(3) * t48;
t142 = -t373 * t375 - t301;
t526 = m(7) * (-t373 * t525 + t420 * t472 - t684);
t27 = t472 * t609 - t526 / 0.2e1 + t560 - t784;
t94 = t618 + t538;
t528 = qJD(1) * t27 + qJD(2) * t94 - qJD(3) * t142;
t522 = mrSges(7,2) * t556 + t790;
t257 = m(7) * t420 + t570;
t216 = t231 * qJD(4);
t211 = m(7) * t368 + t570;
t173 = m(7) * t800 + t525 * t767;
t95 = t618 - t538;
t41 = t710 / 0.2e1 + m(7) * t759 - t712 / 0.2e1 + t709 / 0.2e1 + t833 * t767 + mrSges(7,3) * t754 + t815;
t38 = 0.2e1 * t70 * t767 + t270;
t37 = t509 + t808;
t28 = t526 / 0.2e1 + t560 + t604;
t25 = t521 - t809;
t15 = t543 / 0.2e1 + t432 * mrSges(7,3) / 0.2e1 + mrSges(6,2) * t743 + t506 + t799 * t745;
t13 = -t514 + t845;
t11 = -t515 + t845;
t7 = t503 - t517 + t776;
t5 = t504 - t518 + t776;
t1 = -t502 + (t477 / 0.2e1 - mrSges(4,1) / 0.2e1) * t791 + (t271 / 0.2e1 + t264 / 0.2e1) * t420 + (-t266 / 0.2e1 + t842) * t419 + t520 + t775;
t17 = [qJD(2) * t3 + qJD(3) * t6 + qJD(4) * t16 + qJD(5) * t9 + qJD(6) * t33, t714 + (t844 + t676 + (t472 * t630 + t473 * t627) * mrSges(4,3) + t673 - t674 + t672 + t671 + t654 + t652 + m(4) * (-t422 * t730 - t733 * t791) * pkin(2) + t788 * mrSges(5,3) + m(6) * (-t367 * t88 + t368 * t89 + t836) + Ifges(3,5) * t734 - Ifges(3,6) * t731 - mrSges(3,1) * t631 - t536 + mrSges(3,2) * t628 + t535 + m(5) * (-t250 * t545 + t251 * t544 + t820) + m(7) * (t367 * t73 + t368 * t72 + t822) + (-t693 - t694) * mrSges(6,3) + (t695 - t696) * mrSges(7,2)) * qJD(2) + t1 * qJD(3) + t11 * qJD(4) + t5 * qJD(5) + t25 * qJD(6), t697 + t1 * qJD(2) + (t670 - t669 + t668 + t667 + t666 + t653 - t633 - t634 - t626 - t80 * t704 + t625 - t728 + (t649 - t651) * qJ(4) + t635 + t844) * qJD(3) + t13 * qJD(4) + t7 * qJD(5) + t28 * qJD(6) + 0.2e1 * ((-t826 + (t679 - t680) * qJ(4)) * t771 + (-t419 * t92 + t420 * t93 + t835) * t769 + (t419 * t81 + t420 * t80 + t821) * t767) * qJD(3), qJD(2) * t11 + qJD(3) * t13 + qJD(5) * t41 + qJD(6) * t173 + t688, t692 + t5 * qJD(2) + t7 * qJD(3) + t41 * qJD(4) + (t558 * mrSges(7,2) + t779 * t84 + t780 * t85 + t789) * qJD(5) + t38 * qJD(6), qJD(2) * t25 + qJD(3) * t28 + qJD(4) * t173 + qJD(5) * t38 + t687; -qJD(3) * t2 - qJD(4) * t10 + qJD(5) * t4 + qJD(6) * t26 - t714, qJD(3) * t31 + qJD(4) * t47 + qJD(5) * t20 + qJD(6) * t138, t37 * qJD(4) + t15 * qJD(5) + t95 * qJD(6) + t553 + (m(6) * (t487 * t627 + t575) + m(7) * (t373 * t627 + t575) + (m(5) * (-pkin(3) * t730 + qJ(4) * t533) + t777) * pkin(2) + t630 * t806 + t830) * qJD(3), qJD(3) * t37 - t551 + t646, t15 * qJD(3) + t216 + (-t367 * t779 + t368 * t780 + t522) * qJD(5) + t211 * qJD(6) + t552, qJD(3) * t95 + qJD(5) * t211 - t547; qJD(2) * t2 - qJD(4) * t12 + qJD(5) * t8 - qJD(6) * t27 - t697, -qJD(4) * t36 + qJD(5) * t14 - qJD(6) * t94 - t553, qJD(4) * t48 + qJD(5) * t21 + qJD(6) * t142, -t530 + t646, t216 + (-t419 * t779 + t420 * t780 + t522) * qJD(5) + t257 * qJD(6) + t531, qJD(5) * t257 - t528; qJD(2) * t10 + qJD(3) * t12 - qJD(5) * t40 - qJD(6) * t172 - t688, qJD(3) * t36 + t551 + t642, t530 + t642, 0, -t782, t781; -qJD(2) * t4 - qJD(3) * t8 + qJD(4) * t40 + qJD(6) * t39 - t692, -qJD(3) * t14 - t552 + t647, -t531 + t647, t782, t495 * qJD(6), t546; -qJD(2) * t26 + qJD(3) * t27 + qJD(4) * t172 - qJD(5) * t39 - t687, qJD(3) * t94 + t547 + t641, t528 + t641, -t781, -t546, 0;];
Cq  = t17;
