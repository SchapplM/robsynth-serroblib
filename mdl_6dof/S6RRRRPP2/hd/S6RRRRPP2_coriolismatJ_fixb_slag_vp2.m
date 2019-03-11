% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRRPP2
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
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRRPP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:49:02
% EndTime: 2019-03-09 20:49:32
% DurationCPUTime: 16.59s
% Computational Cost: add. (21638->831), mult. (43392->1056), div. (0->0), fcn. (43813->6), ass. (0->425)
t501 = sin(qJ(4));
t499 = t501 ^ 2;
t502 = cos(qJ(4));
t500 = t502 ^ 2;
t822 = t499 + t500;
t740 = sin(qJ(3));
t741 = sin(qJ(2));
t742 = cos(qJ(3));
t743 = cos(qJ(2));
t429 = t740 * t741 - t742 * t743;
t488 = t502 * mrSges(7,2);
t709 = t501 * mrSges(7,1);
t766 = m(7) / 0.2e1;
t665 = t429 * t501;
t650 = t502 * qJ(5);
t634 = t743 * pkin(7);
t465 = pkin(8) * t743 + t634;
t631 = t741 * pkin(7);
t541 = -pkin(8) * t741 - t631;
t783 = t742 * t465 + t740 * t541;
t798 = -pkin(4) * t665 + t429 * t650 + t783;
t820 = pkin(5) * t665 - t798;
t841 = (-t488 / 0.2e1 + t709 / 0.2e1) * t429 + t820 * t766;
t430 = -t740 * t743 - t741 * t742;
t663 = t430 * t501;
t388 = mrSges(7,3) * t663;
t293 = t429 * mrSges(7,2) - t388;
t662 = t430 * t502;
t298 = -mrSges(7,1) * t429 + mrSges(7,3) * t662;
t745 = t502 / 0.2e1;
t749 = t501 / 0.2e1;
t780 = t293 * t745 + t298 * t749;
t840 = t780 + t841;
t839 = -t780 + t841;
t720 = Ifges(6,5) * t502;
t446 = Ifges(6,3) * t501 + t720;
t181 = -Ifges(6,6) * t430 - t429 * t446;
t721 = Ifges(7,4) * t502;
t449 = Ifges(7,2) * t501 + t721;
t183 = Ifges(7,6) * t430 - t429 * t449;
t495 = Ifges(5,4) * t502;
t562 = Ifges(5,2) * t501 - t495;
t185 = -Ifges(5,6) * t430 + t429 * t562;
t493 = Ifges(7,4) * t501;
t454 = Ifges(7,1) * t502 + t493;
t187 = Ifges(7,5) * t430 - t429 * t454;
t491 = Ifges(6,5) * t501;
t456 = Ifges(6,1) * t502 + t491;
t189 = -Ifges(6,4) * t430 - t429 * t456;
t722 = Ifges(5,4) * t501;
t458 = Ifges(5,1) * t502 - t722;
t191 = -Ifges(5,5) * t430 - t429 * t458;
t451 = Ifges(5,2) * t502 + t722;
t664 = t429 * t502;
t604 = -t664 / 0.2e1;
t605 = t665 / 0.2e1;
t606 = -t665 / 0.2e1;
t747 = -t502 / 0.2e1;
t754 = -t430 / 0.2e1;
t453 = Ifges(7,1) * t501 - t721;
t455 = Ifges(6,1) * t501 - t720;
t457 = Ifges(5,1) * t501 + t495;
t796 = t453 + t455 + t457;
t445 = -Ifges(6,3) * t502 + t491;
t448 = -Ifges(7,2) * t502 + t493;
t801 = t445 + t448;
t811 = Ifges(6,4) + Ifges(5,5);
t523 = t185 * t745 + t451 * t605 + t430 * (Ifges(7,5) * t501 - Ifges(7,6) * t502) / 0.2e1 + Ifges(4,6) * t430 - Ifges(4,5) * t429 + ((Ifges(5,6) - Ifges(6,6)) * t502 + t811 * t501) * t754 + (t181 + t183) * t747 + t801 * t606 + (t187 + t189 + t191) * t749 + t796 * t604;
t438 = -t502 * mrSges(5,1) + t501 * mrSges(5,2);
t803 = t783 * t438;
t809 = t783 * mrSges(4,1);
t343 = t465 * t740 - t742 * t541;
t810 = t343 * mrSges(4,2);
t436 = -t502 * mrSges(6,1) - t501 * mrSges(6,3);
t823 = t798 * t436;
t437 = mrSges(7,1) * t502 + mrSges(7,2) * t501;
t833 = t820 * t437;
t838 = t523 + t803 + t810 - t809 + t823 + t833;
t837 = -t803 / 0.2e1 + t809 / 0.2e1 - t810 / 0.2e1 - t823 / 0.2e1 - t833 / 0.2e1;
t761 = pkin(4) + pkin(5);
t647 = t501 * t761;
t775 = -t650 + t647;
t105 = t430 * t775 - t343;
t836 = t105 * t820;
t633 = t742 * pkin(2);
t480 = -t633 - pkin(3);
t485 = t501 * qJ(5);
t785 = t502 * pkin(4) + t485;
t404 = t480 - t785;
t497 = t502 * pkin(5);
t368 = t497 - t404;
t835 = t368 * t820;
t635 = pkin(3) + t785;
t405 = t497 + t635;
t834 = t405 * t820;
t831 = t822 * t742;
t481 = -pkin(2) * t743 - pkin(1);
t732 = t430 * pkin(9);
t284 = t429 * pkin(3) + t481 + t732;
t125 = t284 * t501 + t502 * t783;
t439 = pkin(4) * t501 - t650;
t155 = -t430 * t439 + t343;
t441 = mrSges(6,1) * t501 - mrSges(6,3) * t502;
t272 = t441 * t429;
t442 = t488 - t709;
t273 = t442 * t429;
t726 = mrSges(5,2) * t502;
t443 = mrSges(5,1) * t501 + t726;
t274 = t443 * t429;
t275 = t441 * t430;
t276 = t442 * t430;
t277 = t443 * t430;
t291 = -mrSges(7,2) * t430 - mrSges(7,3) * t665;
t292 = mrSges(5,2) * t430 + mrSges(5,3) * t665;
t710 = t430 * mrSges(7,1);
t295 = mrSges(7,3) * t664 + t710;
t296 = -mrSges(5,1) * t430 + mrSges(5,3) * t664;
t711 = t430 * mrSges(6,1);
t297 = -mrSges(6,2) * t664 + t711;
t302 = mrSges(6,2) * t665 - mrSges(6,3) * t430;
t646 = -t502 * t284 + t501 * t783;
t543 = qJ(6) * t662 + t646;
t57 = -t429 * t761 + t543;
t484 = t501 * qJ(6);
t385 = t430 * t484;
t406 = t429 * qJ(5);
t91 = t125 + t406;
t68 = -t385 + t91;
t92 = -pkin(4) * t429 + t646;
t828 = -t820 * t276 - t798 * t275 - t783 * t277 - t105 * t273 - t646 * t296 + t125 * t292 - t155 * t272 - t343 * t274 + t68 * t291 + t57 * t295 + t92 * t297 + t91 * t302 + t481 * (-mrSges(4,1) * t430 - mrSges(4,2) * t429);
t789 = mrSges(6,2) + mrSges(5,3);
t827 = t789 / 0.2e1;
t731 = mrSges(6,2) - mrSges(7,3);
t826 = t155 * t798;
t825 = t404 * t798;
t824 = t635 * t798;
t294 = -mrSges(5,2) * t429 + mrSges(5,3) * t663;
t626 = mrSges(6,2) * t663;
t301 = t429 * mrSges(6,3) + t626;
t787 = t301 + t294;
t620 = t761 * t430;
t806 = t501 * t343;
t821 = t620 - t806;
t816 = t275 / 0.2e1;
t815 = -t276 / 0.2e1;
t814 = -t437 / 0.2e1;
t813 = pkin(3) * t783;
t812 = mrSges(5,1) + mrSges(6,1);
t729 = Ifges(5,6) + Ifges(7,6);
t728 = t57 - t543;
t702 = -t646 + t92;
t808 = t343 * t740;
t679 = t343 * t783;
t807 = t480 * t783;
t805 = t502 * t343;
t804 = t502 * t702;
t529 = t789 * (t499 / 0.2e1 + t500 / 0.2e1);
t802 = t436 + t438;
t800 = t449 + t446;
t799 = t831 * pkin(2);
t797 = Ifges(6,2) + Ifges(5,3) + Ifges(7,3);
t794 = -t562 + t796;
t793 = t458 + t456 + t454 + t801;
t792 = 0.2e1 * t766;
t768 = m(6) / 0.2e1;
t791 = 0.2e1 * t768;
t790 = m(6) + m(7);
t579 = t633 / 0.2e1;
t730 = mrSges(6,3) + mrSges(7,2);
t270 = t430 * t437;
t788 = t292 + t302;
t645 = t404 - t635;
t591 = t436 / 0.2e1 + t814;
t784 = (-mrSges(7,1) / 0.2e1 - mrSges(6,1) / 0.2e1 + t591 * t502) * t430 + t731 * t664;
t781 = t761 * t295 / 0.2e1 + pkin(4) * t297 / 0.2e1;
t779 = t788 * pkin(9);
t630 = t740 * pkin(2);
t479 = t630 + pkin(9);
t778 = t788 * t479 * t502;
t771 = 0.2e1 * m(7);
t212 = (-t650 / 0.4e1 + t647 / 0.4e1 + t775 / 0.4e1) * t771 - t442;
t386 = t430 * t485;
t172 = t502 * t620 + t386;
t599 = t662 / 0.4e1;
t48 = -t270 + (-t761 * t599 - t386 / 0.4e1 - t172 / 0.4e1) * t771;
t641 = qJD(2) + qJD(3);
t777 = qJD(1) * t48 + t212 * t641;
t188 = -t429 * Ifges(7,5) - t430 * t454;
t190 = t429 * Ifges(6,4) - t430 * t456;
t192 = t429 * Ifges(5,5) - t430 * t458;
t572 = Ifges(7,5) / 0.2e1 - Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t776 = t572 * t429 - t188 / 0.2e1 - t190 / 0.2e1 - t192 / 0.2e1;
t299 = mrSges(5,1) * t429 + mrSges(5,3) * t662;
t300 = -mrSges(6,1) * t429 - mrSges(6,2) * t662;
t750 = -t501 / 0.2e1;
t774 = t299 * t747 + t300 * t745 + t787 * t750;
t603 = t664 / 0.2e1;
t763 = mrSges(6,2) / 0.2e1;
t773 = mrSges(6,2) * t604 + mrSges(7,3) * t603 + t710 / 0.2e1 + t711 / 0.2e1 + (t591 * t430 + (t763 - mrSges(7,3) / 0.2e1) * t429) * t502;
t268 = -pkin(4) * t662 - t386;
t772 = t439 * t816 + t775 * t815 - t343 * t443 / 0.2e1 - t268 * t436 / 0.2e1 + t172 * t814 - t155 * t441 / 0.2e1 - t105 * t442 / 0.2e1;
t770 = m(5) / 0.2e1;
t769 = -m(6) / 0.2e1;
t767 = -m(7) / 0.2e1;
t765 = m(6) * pkin(2);
t764 = m(7) * pkin(2);
t733 = t429 * pkin(9);
t306 = -t430 * pkin(3) + t733;
t700 = qJ(6) * t429;
t62 = (-t306 + t700) * t502 + t821;
t762 = -t62 / 0.2e1;
t151 = t306 * t502 + t806;
t734 = pkin(4) * t430;
t107 = -t151 + t734;
t759 = -t107 / 0.2e1;
t758 = -t125 / 0.2e1;
t757 = t368 / 0.2e1;
t756 = t404 / 0.2e1;
t435 = pkin(9) * t501 - t484;
t753 = t435 / 0.2e1;
t751 = t480 / 0.2e1;
t739 = m(6) * t479;
t738 = m(7) * t501;
t737 = pkin(3) * t274;
t736 = pkin(3) * t443;
t86 = t125 - t385;
t727 = -t68 + t86;
t632 = t741 * pkin(2);
t286 = t632 + t306;
t139 = t286 * t502 + t806;
t140 = t501 * t286 - t805;
t621 = Ifges(7,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t571 = -Ifges(6,6) / 0.2e1 + t621;
t510 = -Ifges(4,4) * t430 + (-t191 / 0.2e1 - t189 / 0.2e1 - t187 / 0.2e1 - t572 * t430) * t502 + (t185 / 0.2e1 - t183 / 0.2e1 - t181 / 0.2e1 - t571 * t430) * t501 + (Ifges(4,1) - Ifges(4,2) - t797) * t429;
t390 = Ifges(6,5) * t662;
t182 = Ifges(6,6) * t429 - Ifges(6,3) * t663 - t390;
t391 = Ifges(7,4) * t662;
t184 = -Ifges(7,2) * t663 - t429 * Ifges(7,6) - t391;
t186 = t429 * Ifges(5,6) + t430 * t562;
t567 = -t182 / 0.2e1 - t184 / 0.2e1 + t186 / 0.2e1;
t515 = Ifges(4,4) * t429 + t776 * t502 + (t429 * t571 + t567) * t501;
t59 = (-t286 + t700) * t502 + t821;
t384 = t429 * t484;
t407 = t430 * qJ(5);
t98 = -t407 + t140;
t69 = -t384 + t98;
t99 = -t139 + t734;
t4 = (mrSges(4,1) * t632 + t515) * t429 - pkin(1) * (mrSges(3,1) * t741 + mrSges(3,2) * t743) + (-mrSges(4,2) * t632 + t510) * t430 + t59 * t298 + t139 * t299 + t99 * t300 + t98 * t301 + t69 * t293 + t140 * t294 + m(5) * (t125 * t140 - t139 * t646 + t679) + m(6) * (t91 * t98 + t92 * t99 + t826) + m(7) * (t57 * t59 + t68 * t69 + t836) + m(4) * t481 * t632 + (-Ifges(3,2) + Ifges(3,1)) * t743 * t741 + (-t741 ^ 2 + t743 ^ 2) * Ifges(3,4) + t828;
t712 = t4 * qJD(1);
t708 = t501 * mrSges(7,3);
t707 = t501 * t99;
t489 = t502 * mrSges(6,2);
t706 = t502 * mrSges(7,3);
t705 = t502 * t98;
t152 = t501 * t306 - t805;
t103 = -t407 + t152;
t72 = -t384 + t103;
t8 = t515 * t429 + t510 * t430 + t62 * t298 + t151 * t299 + t107 * t300 + t103 * t301 + t72 * t293 + t152 * t294 + m(5) * (t125 * t152 - t151 * t646 + t679) + m(6) * (t103 * t91 + t107 * t92 + t826) + m(7) * (t57 * t62 + t68 * t72 + t836) + t828;
t704 = t8 * qJD(1);
t269 = t430 * t436;
t271 = t438 * t430;
t278 = t430 * t445;
t279 = t430 * t448;
t280 = t430 * t451;
t281 = Ifges(7,1) * t663 - t391;
t282 = Ifges(6,1) * t663 - t390;
t283 = t457 * t430;
t389 = Ifges(6,6) * t662;
t9 = -t429 * t389 / 0.2e1 + t343 * t271 + t86 * t298 - t543 * t293 - t268 * t275 - t172 * t276 + t155 * t269 + t105 * t270 + (-t299 + t300) * t125 - t787 * t646 + m(6) * (t125 * t92 + t155 * t268 - t646 * t91) + m(7) * (t105 * t172 - t543 * t68 + t57 * t86) + ((-t281 / 0.2e1 - t282 / 0.2e1 - t283 / 0.2e1 + t125 * mrSges(5,3) + t91 * mrSges(6,2) - t68 * mrSges(7,3) + t621 * t429 + t567) * t502 + (-t278 / 0.2e1 - t279 / 0.2e1 + t280 / 0.2e1 + t646 * mrSges(5,3) + t92 * mrSges(6,2) - t57 * mrSges(7,3) - t776) * t501) * t430;
t703 = t9 * qJD(1);
t701 = t125 - t91;
t19 = (-t275 + t276) * t662 + (t293 + t301) * t429 + m(6) * (t155 * t662 + t429 * t91) + m(7) * (-t105 * t662 + t429 * t68);
t699 = qJD(1) * t19;
t32 = (-t501 * t293 + m(7) * (-t501 * t68 + t57 * t502) + t502 * t298) * t430;
t698 = qJD(1) * t32;
t695 = t105 * t775;
t692 = t107 * t501;
t691 = t125 * t501;
t690 = t139 * t501;
t689 = t140 * t502;
t688 = t152 * t502;
t686 = t155 * t439;
t684 = t155 * t501;
t676 = t368 * t273;
t675 = t368 * t442;
t674 = t404 * t272;
t673 = t405 * t273;
t672 = t405 * t442;
t408 = t479 * t501 - t484;
t671 = t408 * t295;
t670 = t408 * t501;
t409 = (-qJ(6) + t479) * t502;
t669 = t409 * t291;
t668 = t409 * t502;
t666 = t775 * t437;
t661 = t635 * t272;
t660 = t435 * t295;
t658 = t439 * t436;
t440 = (pkin(9) - qJ(6)) * t502;
t657 = t440 * t291;
t655 = t480 * t274;
t654 = t480 * t443;
t653 = t501 * t296;
t652 = t501 * t297;
t644 = t799 * pkin(9);
t643 = -t212 * qJD(4) + qJD(5) * t738;
t642 = t822 * mrSges(7,3);
t450 = Ifges(6,4) * t502 + Ifges(6,6) * t501;
t640 = m(7) * t662;
t638 = t69 * t706;
t637 = mrSges(6,2) * t705;
t629 = qJD(6) * t738;
t628 = mrSges(5,3) * t690;
t627 = mrSges(5,3) * t689;
t619 = t479 * t653;
t618 = t479 * t652;
t617 = mrSges(6,2) * t750;
t616 = mrSges(5,3) * t750;
t615 = -t708 / 0.2e1;
t614 = t708 / 0.2e1;
t613 = t489 / 0.2e1;
t612 = -t706 / 0.2e1;
t611 = t501 * t742;
t610 = t502 * t742;
t593 = t300 / 0.2e1 - t299 / 0.2e1;
t592 = -t301 / 0.2e1 - t294 / 0.2e1;
t589 = m(6) * t645;
t588 = m(7) * (t368 + t405);
t587 = m(6) * t439 + t441;
t586 = (-t408 - t435) * t501;
t585 = (-t409 - t440) * t502;
t580 = -t633 / 0.2e1;
t578 = -t630 / 0.2e1;
t492 = Ifges(5,5) * t502;
t576 = mrSges(7,3) * t485 + t450 + t492;
t575 = qJ(5) * t610;
t570 = t645 * t768;
t287 = t368 * t775;
t312 = t405 * t775;
t569 = (-t287 - t312) * t766;
t568 = (t816 + t815) * t501;
t561 = -t501 * t57 - t502 * t68;
t560 = t705 + t707;
t559 = t501 * t580;
t557 = t502 * t579;
t518 = -t666 + (t453 / 0.2e1 - t449 / 0.2e1 + t455 / 0.2e1 - t446 / 0.2e1 + t457 / 0.2e1 - t562 / 0.2e1) * t502 + (t454 / 0.2e1 + t448 / 0.2e1 + t456 / 0.2e1 + t445 / 0.2e1 + t458 / 0.2e1 - t451 / 0.2e1) * t501;
t514 = t518 + t658;
t26 = -m(7) * t287 + t404 * t587 + t514 + t654 + t675;
t444 = Ifges(7,5) * t502 + Ifges(7,6) * t501;
t447 = -Ifges(5,6) * t501 + t492;
t509 = -t501 * t186 / 0.4e1 - t429 * t444 / 0.4e1 - t772 + t125 * t616 + t91 * t617 + t68 * t614 + t86 * t615 + t451 * t599 + (t450 + t447) * t429 / 0.4e1 - (t279 + t278) * t502 / 0.4e1 + t702 * t613 + t728 * t612 + t691 * t827 + (t280 + t192 + t190 + t188) * t502 / 0.4e1 + (t283 + t282 + t281 + t184 + t182) * t501 / 0.4e1 - t793 * t662 / 0.4e1 + (-t800 / 0.4e1 + t794 / 0.4e1) * t663;
t506 = t509 + t271 * t751 + t269 * t756 - t408 * t293 / 0.2e1 + t409 * t298 / 0.2e1 + t270 * t757 + (t268 * t404 + t686) * t768 + (t172 * t368 + t408 * t727 + t409 * t728 - t695) * t766;
t511 = (-pkin(4) * t99 + qJ(5) * t98) * t769 + (qJ(5) * t69 - t59 * t761) * t767 - t139 * mrSges(5,1) / 0.2e1 + t140 * mrSges(5,2) / 0.2e1 + t59 * mrSges(7,1) / 0.2e1 - t69 * mrSges(7,2) / 0.2e1 - t98 * mrSges(6,3) / 0.2e1 + t99 * mrSges(6,1) / 0.2e1;
t535 = t501 * t701 + t804;
t538 = (-t670 / 0.2e1 - t668 / 0.2e1) * mrSges(7,3);
t3 = t506 + t511 + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + t538 + t529 * t479) * t430 + (-t501 * t571 - t502 * t572) * t429 + (-t302 / 0.2e1 - t291 / 0.2e1) * qJ(5) + (t501 * t592 + t502 * t593) * t479 + t535 * t739 / 0.2e1 + t781;
t556 = t3 * qJD(1) + t26 * qJD(2);
t513 = -mrSges(4,2) * t633 + (-mrSges(4,1) - t437 + t802) * t630 + (mrSges(5,3) + t731) * t799;
t534 = t831 * t479;
t36 = -(t404 * t740 + t534) * t765 - (-t368 * t740 + t408 * t611 + t409 * t610) * t764 - m(5) * (t480 * t740 + t534) * pkin(2) - t513;
t552 = -t151 * t501 + t688;
t554 = t103 * t502 + t692;
t504 = (t807 + t552 * t479 + (t125 * t610 + t611 * t646 + t808) * pkin(2)) * t770 + mrSges(5,3) * t688 / 0.2e1 - t676 / 0.2e1 - t674 / 0.2e1 + t671 / 0.2e1 + t669 / 0.2e1 - t655 / 0.2e1 + t778 / 0.2e1 + t299 * t559 + t618 / 0.2e1 - t837 + (-t277 - t275) * t630 / 0.2e1 + t62 * t615 + t151 * t616 - t276 * t578 + t72 * t612 + t103 * t613 + (t825 + t554 * t479 + (t155 * t740 + t610 * t91 + t611 * t92) * pkin(2)) * t768 - t619 / 0.2e1 + t692 * t763 + (t300 + t298) * t501 * t579 + (t835 + t408 * t62 + t409 * t72 + (-t105 * t740 + t57 * t611 + t610 * t68) * pkin(2)) * t766 + (t293 + t787) * t557;
t553 = t689 - t690;
t505 = t673 / 0.2e1 - t660 / 0.2e1 - t661 / 0.2e1 - t657 / 0.2e1 + t628 / 0.2e1 - t627 / 0.2e1 + t638 / 0.2e1 - t637 / 0.2e1 + t99 * t617 + t747 * t779 + t59 * t614 + (pkin(9) * t560 - t824) * t769 + (t435 * t59 + t440 * t69 + t834) * t767 - t737 / 0.2e1 - pkin(9) * t652 / 0.2e1 + pkin(9) * t653 / 0.2e1 - m(5) * (pkin(9) * t553 - t813) / 0.2e1 + t837;
t7 = t504 + t505;
t555 = t7 * qJD(1) - t36 * qJD(2);
t551 = -t268 * t635 + t686;
t428 = t501 * t437;
t121 = -t368 * t738 - t428 + (m(6) * t404 + t436) * t501;
t101 = t501 * t105;
t516 = t568 + (-t684 + (t404 * t430 + t429 * t479) * t502) * t768 + (-t368 * t662 + t409 * t429 + t101) * t766;
t547 = t59 * t767 + t769 * t99;
t13 = t516 + t547 + t784;
t550 = qJD(1) * t13 - qJD(2) * t121;
t179 = m(7) * (-t668 - t670) + t642;
t526 = m(7) * ((t408 * t502 - t409 * t501) * t430 + t561);
t21 = -t526 / 0.2e1 + t840;
t549 = -qJD(1) * t21 + qJD(2) * t179;
t542 = t125 + 0.2e1 * t406;
t520 = t730 * t429 + t542 * t768 + (-t385 + t542) * t766;
t544 = m(6) * t758 + t767 * t86;
t29 = t520 + t544;
t468 = qJ(5) * t790 + t730;
t548 = qJD(1) * t29 + qJD(4) * t468;
t546 = -mrSges(6,2) * pkin(4) + mrSges(7,3) * t761 - Ifges(7,5);
t545 = m(6) * t759 + m(7) * t762;
t507 = (-pkin(4) * t611 + t575) * t765 / 0.2e1 + (-t611 * t761 + t575) * t764 / 0.2e1 + t580 * t726 + t730 * t557 + (mrSges(7,1) + t812) * t559;
t10 = t439 * t570 - t507 + t451 * t750 + t654 / 0.2e1 + t672 / 0.2e1 + t675 / 0.2e1 - t736 / 0.2e1 - t666 + t658 + t569 + t645 * t441 / 0.2e1 + t800 * t747 + t794 * t745 + t793 * t749;
t30 = -m(7) * t312 - t587 * t635 + t514 + t672 - t736;
t527 = Ifges(7,5) * t603 + Ifges(6,6) * t606 - t781 + (t291 + t302) * qJ(5) / 0.2e1 + t729 * t605 + t811 * t604 + t797 * t754;
t508 = t151 * mrSges(5,1) / 0.2e1 - t152 * mrSges(5,2) / 0.2e1 + t103 * mrSges(6,3) / 0.2e1 + mrSges(6,1) * t759 + t72 * mrSges(7,2) / 0.2e1 + mrSges(7,1) * t762 + t527 + (qJ(5) * t72 - t62 * t761) * t766 + (-pkin(4) * t107 + qJ(5) * t103) * t768;
t512 = (t172 * t405 + t435 * t727 + t440 * t728 - t695) * t767 + pkin(3) * t271 / 0.2e1 - t405 * t270 / 0.2e1 + t635 * t269 / 0.2e1 + t293 * t753 - t440 * t298 / 0.2e1;
t5 = (-t283 / 0.4e1 - t282 / 0.4e1 - t281 / 0.4e1 + t186 / 0.4e1 - t184 / 0.4e1 - t182 / 0.4e1 + (t86 / 0.2e1 - t68 / 0.2e1) * mrSges(7,3) + (t758 + t91 / 0.2e1) * mrSges(6,2) + (t701 * t769 - t592) * pkin(9)) * t501 + (-t192 / 0.4e1 - t190 / 0.4e1 - t188 / 0.4e1 - t280 / 0.4e1 + t279 / 0.4e1 + t278 / 0.4e1 + (-t543 / 0.2e1 + t57 / 0.2e1) * mrSges(7,3) + (t646 / 0.2e1 - t92 / 0.2e1) * mrSges(6,2) + (t702 * t769 - t593) * pkin(9)) * t502 + t508 + t512 + ((t458 / 0.4e1 + t456 / 0.4e1 + t454 / 0.4e1 - t451 / 0.4e1 + t448 / 0.4e1 + t445 / 0.4e1 + t440 * mrSges(7,3) / 0.2e1) * t502 + (t562 / 0.4e1 + t449 / 0.4e1 + t446 / 0.4e1 - t457 / 0.4e1 - t455 / 0.4e1 - t453 / 0.4e1 + mrSges(7,3) * t753) * t501 - pkin(9) * t529) * t430 + t551 * t769 + (t444 / 0.4e1 - t447 / 0.4e1 - t450 / 0.4e1) * t429 + t772;
t539 = -t5 * qJD(1) + t10 * qJD(2) + t30 * qJD(3);
t517 = t568 + (-t684 + (-t430 * t635 + t733) * t502) * t768 + (-t405 * t662 + t429 * t440 + t101) * t766;
t15 = t517 + t545 + t784;
t153 = -t405 * t738 - t428 + (-m(6) * t635 + t436) * t501;
t533 = t790 * t579;
t49 = -t428 + (t436 + t533 + t589 / 0.2e1 - t588 / 0.2e1) * t501;
t537 = qJD(1) * t15 - qJD(2) * t49 - qJD(3) * t153;
t525 = m(7) * ((t435 * t502 - t440 * t501) * t430 + t561);
t23 = -t525 / 0.2e1 + t840;
t234 = m(7) * (-t435 * t501 - t440 * t502) + t642;
t95 = (t578 - t585 / 0.2e1 - t586 / 0.2e1) * m(7) - t642;
t536 = -qJD(1) * t23 - qJD(2) * t95 + qJD(3) * t234;
t522 = ((-mrSges(6,2) * qJ(5) - t729) * t501 + t546 * t502) * qJD(4);
t521 = qJD(4) * (-m(6) * t785 + t802);
t317 = m(7) * t440 + t489 + (m(6) * pkin(9) - mrSges(7,3)) * t502;
t304 = (-qJD(1) * t662 + t501 * t641) * m(7);
t243 = m(7) * t409 + t489 + (-mrSges(7,3) + t739) * t502;
t202 = t212 * qJD(6);
t96 = (t585 + t586) * t766 + m(7) * t578 + t642;
t70 = (-t662 * t761 + t172 - t386) * t766;
t50 = t428 + (-t436 + t533) * t501 + (-t589 + t588) * t749;
t28 = -t388 + t520 - t544 + t626;
t24 = t525 / 0.2e1 + t839;
t22 = t526 / 0.2e1 + t839;
t14 = t517 - t545 + t773;
t12 = t516 - t547 + t773;
t11 = (t436 + t570) * t439 + t507 + t569 + (t751 - pkin(3) / 0.2e1) * t443 + (t405 / 0.2e1 + t757) * t442 + (-t635 / 0.2e1 + t756) * t441 + t518;
t6 = t551 * t768 + t508 + t509 - t512 + (t435 * t615 + t440 * t612) * t430 + (t535 * t768 + t774) * pkin(9) + t822 * t732 * t827;
t2 = ((-t91 * t501 + t691 + t804) * t768 + t529 * t430 + t774) * t479 + t506 - t511 + t527 + t430 * t538;
t1 = t504 - t505 + t523;
t16 = [qJD(2) * t4 + qJD(3) * t8 + qJD(4) * t9 + qJD(5) * t19 + qJD(6) * t32, t712 + (t838 + m(4) * (-t742 * t783 - t808) * pkin(2) + m(5) * (t479 * t553 + t807) - t59 * t708 - t676 - t674 + t671 + t669 - t655 - Ifges(3,6) * t741 + t778 - t628 + t627 - t638 + t637 + t618 + (t429 * t633 + t430 * t630) * mrSges(4,3) + m(6) * (t479 * t560 + t825) - mrSges(3,1) * t634 - t619 + mrSges(3,2) * t631 + Ifges(3,5) * t743 + m(7) * (t408 * t59 + t409 * t69 + t835) + mrSges(6,2) * t707) * qJD(2) + t1 * qJD(3) + t2 * qJD(4) + t12 * qJD(5) + t22 * qJD(6), t1 * qJD(2) + t6 * qJD(4) + t14 * qJD(5) + t24 * qJD(6) + t704 + (t657 + t660 + t661 - t673 + t737 + (t103 * mrSges(6,2) + t152 * mrSges(5,3) - t72 * mrSges(7,3) + t779) * t502 + (t107 * mrSges(6,2) - t151 * mrSges(5,3) - t62 * mrSges(7,3) + (-t296 + t297) * pkin(9)) * t501 + (pkin(9) * t554 - t824) * t791 + 0.2e1 * (pkin(9) * t552 - t813) * t770 + (t435 * t62 + t440 * t72 + t834) * t792 + t838) * qJD(3), t2 * qJD(2) + t6 * qJD(3) + t28 * qJD(5) + t70 * qJD(6) + t703 + (-t86 * mrSges(7,1) - t543 * mrSges(7,2) - t389 + (-qJ(5) * t543 - t761 * t86) * t792 + ((qJ(5) * t731 + t729) * t502 + (t546 + t811) * t501) * t430 + (-pkin(4) * t791 - t812) * t125 + (-qJ(5) * t791 + mrSges(5,2) - mrSges(6,3)) * t646) * qJD(4), qJD(2) * t12 + qJD(3) * t14 + qJD(4) * t28 + t699, qJD(2) * t22 + qJD(3) * t24 + qJD(4) * t70 + t698; qJD(3) * t7 + qJD(4) * t3 + qJD(5) * t13 - qJD(6) * t21 - t712, -qJD(3) * t36 + qJD(4) * t26 - qJD(5) * t121 + qJD(6) * t179 (m(6) * (-t630 * t635 + t644) + m(5) * (-pkin(3) * t630 + t644) + (-t405 * t740 + t435 * t611 + t440 * t610) * t764 + t513) * qJD(3) + t11 * qJD(4) + t50 * qJD(5) + t96 * qJD(6) + t555, t11 * qJD(3) + (m(7) * (-qJ(5) * t408 - t409 * t761) - t409 * mrSges(7,1) - t408 * mrSges(7,2) + t576) * qJD(4) + t243 * qJD(5) + t479 * t521 + t522 + t556, qJD(3) * t50 + qJD(4) * t243 + t550, qJD(3) * t96 + t549; -qJD(2) * t7 - qJD(4) * t5 + qJD(5) * t15 - qJD(6) * t23 - t704, qJD(4) * t10 - qJD(5) * t49 - qJD(6) * t95 - t555, qJD(4) * t30 - qJD(5) * t153 + qJD(6) * t234 (m(7) * (-qJ(5) * t435 - t440 * t761) - t440 * mrSges(7,1) - t435 * mrSges(7,2) + t576) * qJD(4) + t317 * qJD(5) + pkin(9) * t521 + t522 + t539, qJD(4) * t317 + t537, t536; -qJD(2) * t3 + qJD(3) * t5 + qJD(5) * t29 + qJD(6) * t48 - t703, -qJD(3) * t10 + t202 - t556, t202 - t539, t468 * qJD(5), t548, t777; -qJD(2) * t13 - qJD(3) * t15 - qJD(4) * t29 + qJD(6) * t640 - t699, qJD(3) * t49 - t550 - t629, -t537 - t629, -t548, 0, -t304; qJD(2) * t21 + qJD(3) * t23 - qJD(4) * t48 - qJD(5) * t640 - t698, qJD(3) * t95 - t549 + t643, -t536 + t643, -t777, t304, 0;];
Cq  = t16;
