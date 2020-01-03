% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:16
% EndTime: 2019-12-31 22:08:52
% DurationCPUTime: 16.70s
% Computational Cost: add. (19078->926), mult. (47629->1281), div. (0->0), fcn. (48239->8), ass. (0->427)
t550 = sin(qJ(4));
t724 = -qJ(5) - pkin(9);
t511 = t724 * t550;
t761 = -t511 / 0.2e1;
t553 = cos(qJ(4));
t715 = Ifges(6,4) * t553;
t521 = t550 * Ifges(6,1) + t715;
t717 = Ifges(5,4) * t553;
t522 = t550 * Ifges(5,1) + t717;
t823 = t521 / 0.4e1 + t522 / 0.4e1;
t839 = mrSges(6,3) * t761 + t823;
t555 = cos(qJ(2));
t679 = cos(pkin(5));
t626 = pkin(1) * t679;
t549 = sin(pkin(5));
t552 = sin(qJ(2));
t670 = t549 * t552;
t468 = -pkin(7) * t670 + t555 * t626;
t469 = (pkin(2) * t552 - pkin(8) * t555) * t549;
t551 = sin(qJ(3));
t554 = cos(qJ(3));
t287 = -t551 * t468 + t469 * t554;
t240 = -pkin(3) * t670 - t287;
t662 = t554 * t555;
t405 = (-t550 * t662 + t552 * t553) * t549;
t142 = -pkin(4) * t405 + t240;
t406 = (t550 * t552 + t553 * t662) * t549;
t696 = t406 * mrSges(6,2);
t697 = t405 * mrSges(6,1);
t230 = t696 - t697;
t231 = -mrSges(5,1) * t405 + mrSges(5,2) * t406;
t288 = t554 * t468 + t551 * t469;
t664 = t551 * t555;
t633 = t549 * t664;
t308 = -mrSges(6,2) * t633 + t405 * mrSges(6,3);
t310 = mrSges(6,1) * t633 - t406 * mrSges(6,3);
t514 = t724 * t553;
t729 = pkin(4) * t553;
t544 = -pkin(3) - t729;
t591 = mrSges(5,1) * t553 - mrSges(5,2) * t550;
t716 = Ifges(6,4) * t550;
t518 = t553 * Ifges(6,2) + t716;
t718 = Ifges(5,4) * t550;
t519 = t553 * Ifges(5,2) + t718;
t607 = t518 / 0.4e1 + t519 / 0.4e1;
t241 = pkin(9) * t670 + t288;
t731 = pkin(3) * t551;
t526 = -pkin(9) * t554 + t731;
t600 = t552 * t626;
t669 = t549 * t555;
t301 = t600 + (pkin(7) + t526) * t669;
t99 = -t550 * t241 + t553 * t301;
t64 = pkin(4) * t633 - t406 * qJ(5) + t99;
t100 = t553 * t241 + t550 * t301;
t677 = t100 * t553;
t680 = t99 * t550;
t72 = qJ(5) * t405 + t100;
t746 = -t544 / 0.2e1;
t759 = t514 / 0.2e1;
t512 = -mrSges(6,1) * t553 + mrSges(6,2) * t550;
t760 = -t512 / 0.2e1;
t809 = -m(6) / 0.2e1;
t838 = -t823 * t406 - m(5) * (-pkin(3) * t240 + (t677 - t680) * pkin(9)) / 0.2e1 + (t142 * t544 + t511 * t64 - t514 * t72) * t809 + pkin(3) * t231 / 0.2e1 + t142 * t760 + t240 * t591 / 0.2e1 - t287 * mrSges(4,1) / 0.2e1 + t288 * mrSges(4,2) / 0.2e1 + t310 * t761 + t308 * t759 + t230 * t746 - t607 * t405;
t516 = Ifges(6,5) * t550 + Ifges(6,6) * t553;
t517 = Ifges(5,5) * t550 + Ifges(5,6) * t553;
t836 = -t517 / 0.2e1 - t516 / 0.2e1;
t590 = mrSges(5,1) * t550 + mrSges(5,2) * t553;
t474 = t590 * t551;
t835 = pkin(8) * t474;
t635 = m(6) * pkin(4) + mrSges(6,1);
t466 = t679 * t551 + t554 * t670;
t337 = -t466 * t550 - t553 * t669;
t465 = t551 * t670 - t679 * t554;
t338 = t466 * t553 - t550 * t669;
t704 = t338 * Ifges(6,4);
t112 = t337 * Ifges(6,2) + t465 * Ifges(6,6) + t704;
t705 = t338 * Ifges(5,4);
t113 = t337 * Ifges(5,2) + t465 * Ifges(5,6) + t705;
t833 = t113 + t112;
t334 = Ifges(6,4) * t337;
t114 = t338 * Ifges(6,1) + t465 * Ifges(6,5) + t334;
t335 = Ifges(5,4) * t337;
t115 = t338 * Ifges(5,1) + t465 * Ifges(5,5) + t335;
t832 = t115 + t114;
t166 = -Ifges(6,2) * t338 + t334;
t167 = -Ifges(5,2) * t338 + t335;
t831 = t167 + t166;
t168 = Ifges(6,1) * t337 - t704;
t169 = Ifges(5,1) * t337 - t705;
t830 = t169 + t168;
t188 = Ifges(6,4) * t406 + Ifges(6,2) * t405 + Ifges(6,6) * t633;
t189 = Ifges(5,4) * t406 + Ifges(5,2) * t405 + Ifges(5,6) * t633;
t829 = t189 + t188;
t190 = Ifges(6,1) * t406 + Ifges(6,4) * t405 + Ifges(6,5) * t633;
t191 = Ifges(5,1) * t406 + Ifges(5,4) * t405 + Ifges(5,5) * t633;
t828 = t191 + t190;
t585 = -Ifges(6,2) * t550 + t715;
t574 = t585 * t551;
t426 = -t554 * Ifges(6,6) + t574;
t586 = -Ifges(5,2) * t550 + t717;
t575 = t586 * t551;
t428 = -t554 * Ifges(5,6) + t575;
t827 = t426 + t428;
t587 = Ifges(6,1) * t553 - t716;
t576 = t587 * t551;
t430 = -t554 * Ifges(6,5) + t576;
t588 = Ifges(5,1) * t553 - t718;
t577 = t588 * t551;
t432 = -t554 * Ifges(5,5) + t577;
t826 = t430 + t432;
t825 = t518 + t519;
t583 = Ifges(6,5) * t553 - Ifges(6,6) * t550;
t584 = Ifges(5,5) * t553 - Ifges(5,6) * t550;
t824 = t584 + t583;
t821 = -t517 / 0.4e1 - t516 / 0.4e1;
t819 = t428 / 0.2e1 + t426 / 0.2e1;
t818 = t550 ^ 2 + t553 ^ 2;
t605 = t522 / 0.2e1 + t521 / 0.2e1;
t612 = t428 / 0.4e1 + t426 / 0.4e1;
t730 = pkin(4) * t550;
t627 = pkin(8) + t730;
t504 = t627 * t551;
t733 = m(6) * t504;
t685 = t553 * mrSges(6,2);
t687 = t550 * mrSges(6,1);
t589 = t685 + t687;
t473 = t589 * t551;
t771 = t473 / 0.2e1;
t817 = (t771 + t733 / 0.2e1) * pkin(4) - t612;
t470 = pkin(7) * t669 + t600;
t816 = -t470 * mrSges(3,1) - t468 * mrSges(3,2);
t692 = t465 * mrSges(6,2);
t708 = t337 * mrSges(6,3);
t220 = -t692 + t708;
t693 = t465 * mrSges(6,1);
t706 = t338 * mrSges(6,3);
t222 = t693 - t706;
t727 = pkin(9) * t551;
t510 = -pkin(3) * t554 - pkin(2) - t727;
t493 = t553 * t510;
t665 = t551 * t553;
t593 = -qJ(5) * t665 + t493;
t324 = (-pkin(8) * t550 - pkin(4)) * t554 + t593;
t663 = t553 * t554;
t408 = pkin(8) * t663 + t510 * t550;
t667 = t550 * t551;
t342 = -qJ(5) * t667 + t408;
t415 = -t679 * pkin(2) - t468;
t198 = t465 * pkin(3) - t466 * pkin(9) + t415;
t416 = pkin(8) * t679 + t470;
t417 = (-pkin(2) * t555 - pkin(8) * t552 - pkin(1)) * t549;
t229 = t554 * t416 + t551 * t417;
t203 = -pkin(9) * t669 + t229;
t66 = t198 * t550 + t203 * t553;
t46 = qJ(5) * t337 + t66;
t65 = t553 * t198 - t203 * t550;
t45 = -qJ(5) * t338 + t65;
t36 = pkin(4) * t465 + t45;
t698 = t36 * t553;
t740 = -t553 / 0.2e1;
t744 = -t550 / 0.2e1;
t684 = t554 * mrSges(6,1);
t500 = -mrSges(6,3) * t665 - t684;
t768 = -t500 / 0.2e1;
t652 = mrSges(6,3) * t667;
t683 = t554 * mrSges(6,2);
t496 = -t652 + t683;
t770 = t496 / 0.2e1;
t808 = m(6) / 0.2e1;
t815 = (-t324 * t338 + t337 * t342 + (-t46 * t550 - t698) * t551) * t808 + t337 * t770 + t338 * t768 + (t220 * t744 + t222 * t740) * t551;
t638 = Ifges(5,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t641 = Ifges(5,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t814 = t550 * t638 - t553 * t641;
t810 = m(5) / 0.2e1;
t806 = -mrSges(5,1) / 0.2e1;
t805 = -mrSges(5,2) / 0.2e1;
t804 = -mrSges(6,2) / 0.2e1;
t803 = -mrSges(5,3) / 0.2e1;
t802 = mrSges(5,3) / 0.2e1;
t801 = -mrSges(6,3) / 0.2e1;
t673 = t465 * t553;
t228 = -t551 * t416 + t417 * t554;
t297 = pkin(3) * t466 + pkin(9) * t465;
t89 = -t228 * t550 + t553 * t297;
t55 = pkin(4) * t466 + qJ(5) * t673 + t89;
t800 = t55 / 0.2e1;
t799 = -t64 / 0.2e1;
t798 = pkin(3) * mrSges(5,1);
t797 = pkin(3) * mrSges(5,2);
t796 = pkin(4) * t64;
t202 = pkin(3) * t669 - t228;
t107 = -pkin(4) * t337 + t202;
t795 = -t107 / 0.2e1;
t794 = t107 / 0.2e1;
t161 = mrSges(5,1) * t338 + mrSges(5,2) * t337;
t793 = t161 / 0.2e1;
t792 = -t202 / 0.2e1;
t791 = t202 / 0.2e1;
t790 = t220 / 0.2e1;
t789 = t222 / 0.2e1;
t788 = t229 / 0.2e1;
t282 = mrSges(6,1) * t466 + mrSges(6,3) * t673;
t787 = t282 / 0.2e1;
t786 = t324 / 0.2e1;
t411 = pkin(8) * t667 + t553 * t526;
t325 = pkin(4) * t551 - qJ(5) * t663 + t411;
t785 = t325 / 0.2e1;
t784 = t337 / 0.2e1;
t782 = t338 / 0.2e1;
t719 = Ifges(4,4) * t551;
t780 = (Ifges(4,5) * t552 + (Ifges(4,1) * t554 - t719) * t555) * t549 / 0.2e1;
t778 = t406 / 0.2e1;
t666 = t550 * t554;
t656 = pkin(8) * t666;
t407 = t493 - t656;
t777 = t407 / 0.2e1;
t776 = t411 / 0.2e1;
t651 = mrSges(5,3) * t667;
t497 = mrSges(5,2) * t554 - t651;
t769 = t497 / 0.2e1;
t767 = t500 / 0.2e1;
t501 = -mrSges(5,1) * t554 - mrSges(5,3) * t665;
t766 = t501 / 0.2e1;
t502 = mrSges(6,1) * t551 - mrSges(6,3) * t663;
t765 = t502 / 0.2e1;
t764 = -t504 / 0.2e1;
t763 = t504 / 0.2e1;
t505 = t627 * t554;
t762 = t505 / 0.2e1;
t745 = t544 / 0.2e1;
t743 = -t550 / 0.4e1;
t742 = t550 / 0.2e1;
t741 = t550 / 0.4e1;
t738 = t553 / 0.2e1;
t737 = t553 / 0.4e1;
t735 = t554 / 0.2e1;
t734 = m(6) * t107;
t732 = m(6) * t544;
t728 = pkin(8) * t551;
t726 = pkin(9) * t553;
t725 = -Ifges(5,6) - Ifges(6,6);
t723 = -t36 + t45;
t722 = Ifges(4,1) * t466;
t721 = Ifges(3,4) * t552;
t720 = Ifges(4,4) * t466;
t548 = Ifges(4,4) * t554;
t547 = Ifges(4,5) * t554;
t714 = Ifges(5,5) * t406;
t713 = Ifges(6,5) * t406;
t712 = Ifges(5,6) * t405;
t711 = Ifges(6,6) * t405;
t710 = Ifges(4,3) * t552;
t709 = t324 * t46;
t707 = t338 * mrSges(5,3);
t341 = t593 - t656;
t703 = t341 * t46;
t702 = t342 * mrSges(6,3);
t701 = t342 * t36;
t700 = t342 * t45;
t699 = t36 * t550;
t695 = t408 * mrSges(5,3);
t694 = t46 * t553;
t691 = t465 * mrSges(4,3);
t690 = t466 * mrSges(4,3);
t110 = Ifges(6,5) * t338 + Ifges(6,6) * t337 + t465 * Ifges(6,3);
t111 = Ifges(5,5) * t338 + Ifges(5,6) * t337 + t465 * Ifges(5,3);
t162 = -mrSges(6,1) * t337 + mrSges(6,2) * t338;
t163 = -mrSges(5,1) * t337 + mrSges(5,2) * t338;
t186 = Ifges(6,3) * t633 + t711 + t713;
t187 = Ifges(5,3) * t633 + t712 + t714;
t221 = -mrSges(5,2) * t465 + t337 * mrSges(5,3);
t223 = mrSges(5,1) * t465 - t707;
t243 = -Ifges(4,2) * t465 - Ifges(4,6) * t669 + t720;
t452 = Ifges(4,4) * t465;
t244 = -Ifges(4,5) * t669 - t452 + t722;
t309 = -mrSges(5,2) * t633 + t405 * mrSges(5,3);
t311 = mrSges(5,1) * t633 - t406 * mrSges(5,3);
t339 = (Ifges(4,6) * t552 + (-Ifges(4,2) * t551 + t548) * t555) * t549;
t368 = mrSges(4,2) * t669 - t691;
t369 = -mrSges(4,1) * t669 - t690;
t515 = mrSges(4,1) * t551 + mrSges(4,2) * t554;
t409 = t515 * t669;
t442 = (-mrSges(4,2) * t552 - mrSges(4,3) * t664) * t549;
t443 = (mrSges(4,1) * t552 - mrSges(4,3) * t662) * t549;
t536 = Ifges(3,5) * t669;
t622 = t669 / 0.2e1;
t596 = t554 * t622;
t598 = t551 * t622;
t623 = -t669 / 0.2e1;
t599 = t551 * t623;
t624 = t670 / 0.2e1;
t5 = (0.2e1 * Ifges(3,4) * t669 + (Ifges(3,1) - Ifges(3,2)) * t670) * t622 + (t187 + t186) * t465 / 0.2e1 + t828 * t782 + t829 * t784 + t832 * t778 + t833 * t405 / 0.2e1 + (-Ifges(3,6) * t670 + t536 / 0.2e1 + Ifges(3,5) * t622 + t816) * t679 + (t111 + t110) * t598 + (-(Ifges(3,2) * t555 + t721) * t670 / 0.2e1 + (-pkin(1) * (mrSges(3,1) * t552 + mrSges(3,2) * t555) - t555 * (t710 + t555 * (-Ifges(4,6) * t551 + t547)) / 0.2e1 + t552 * (Ifges(3,1) * t555 - t721) / 0.2e1) * t549) * t549 + t466 * t780 + (Ifges(4,5) * t466 - Ifges(4,6) * t465 - Ifges(4,3) * t669) * t624 + t470 * (mrSges(4,1) * t465 + mrSges(4,2) * t466) - t465 * t339 / 0.2e1 + t229 * t442 + t228 * t443 + t415 * t409 + t288 * t368 + t287 * t369 + t46 * t308 + t66 * t309 + t36 * t310 + t65 * t311 + t240 * t163 + t107 * t230 + t202 * t231 + t72 * t220 + t100 * t221 + t64 * t222 + t99 * t223 + m(4) * (t228 * t287 + t229 * t288 + t415 * t470) + m(5) * (t100 * t66 + t202 * t240 + t65 * t99) + m(6) * (t107 * t142 + t36 * t64 + t46 * t72) + t142 * t162 + t243 * t599 + t244 * t596;
t689 = t5 * qJD(1);
t688 = t514 * mrSges(6,3);
t686 = t550 * mrSges(6,3);
t674 = t465 * t550;
t143 = -pkin(4) * t674 + t229;
t179 = Ifges(6,3) * t466 - t465 * t583;
t180 = Ifges(5,3) * t466 - t465 * t584;
t250 = t589 * t465;
t251 = t590 * t465;
t280 = -mrSges(6,2) * t466 + mrSges(6,3) * t674;
t281 = -mrSges(5,2) * t466 + mrSges(5,3) * t674;
t283 = mrSges(5,1) * t466 + mrSges(5,3) * t673;
t294 = mrSges(4,1) * t466 - mrSges(4,2) * t465;
t295 = -Ifges(4,2) * t466 - t452;
t296 = -Ifges(4,1) * t465 - t720;
t183 = Ifges(6,5) * t466 - t465 * t587;
t184 = Ifges(5,5) * t466 - t465 * t588;
t613 = t183 / 0.2e1 + t184 / 0.2e1;
t181 = Ifges(6,6) * t466 - t465 * t585;
t182 = Ifges(5,6) * t466 - t465 * t586;
t614 = t181 / 0.2e1 + t182 / 0.2e1;
t615 = -t114 / 0.2e1 - t115 / 0.2e1;
t616 = t112 / 0.2e1 + t113 / 0.2e1;
t660 = -Ifges(4,5) * t465 - Ifges(4,6) * t466;
t90 = t553 * t228 + t550 * t297;
t71 = qJ(5) * t674 + t90;
t6 = t660 * t623 + t415 * t294 + t228 * t368 + t46 * t280 + t66 * t281 + t36 * t282 + t65 * t283 - t202 * t251 - t107 * t250 + t71 * t220 + t90 * t221 + t55 * t222 + t89 * t223 + t143 * t162 + t613 * t338 + t614 * t337 + (-t369 + t163) * t229 + m(5) * (t202 * t229 + t65 * t89 + t66 * t90) + m(6) * (t107 * t143 + t36 * t55 + t46 * t71) + (-t229 * mrSges(4,3) - t243 / 0.2e1 + t296 / 0.2e1 + t110 / 0.2e1 + t111 / 0.2e1) * t466 + (t228 * mrSges(4,3) + t179 / 0.2e1 + t180 / 0.2e1 - t295 / 0.2e1 - t244 / 0.2e1 + t615 * t553 + t616 * t550) * t465;
t682 = t6 * qJD(1);
t329 = t337 * mrSges(6,2);
t160 = mrSges(6,1) * t338 + t329;
t164 = Ifges(6,5) * t337 - Ifges(6,6) * t338;
t165 = Ifges(5,5) * t337 - Ifges(5,6) * t338;
t634 = m(6) * t723;
t9 = t107 * t160 + t202 * t161 + t45 * t220 + t65 * t221 - t66 * t223 + (t164 / 0.2e1 + t165 / 0.2e1) * t465 + (t634 - t222) * t46 + (t169 / 0.2e1 - t46 * mrSges(6,3) + t168 / 0.2e1 - t66 * mrSges(5,3) + (t162 + t734) * pkin(4) - t616) * t338 + (t167 / 0.2e1 - t36 * mrSges(6,3) + t166 / 0.2e1 - t65 * mrSges(5,3) - t615) * t337;
t681 = t9 * qJD(1);
t23 = m(6) * (t337 * t46 - t338 * t36) - t338 * t222 + t337 * t220;
t678 = qJD(1) * t23;
t676 = t107 * t550;
t675 = t338 * t544;
t672 = t511 * t496;
t671 = t511 * t551;
t668 = t550 * t337;
t661 = -t324 + t341;
t659 = qJD(3) * t550;
t658 = qJD(3) * t553;
t657 = pkin(4) * t665;
t655 = pkin(9) * t802;
t650 = t143 * t808;
t649 = m(6) * t762;
t648 = pkin(4) * t512 / 0.2e1;
t647 = t730 / 0.2e1;
t646 = pkin(9) * t744;
t645 = -t726 / 0.2e1;
t644 = -Ifges(5,1) / 0.4e1 - Ifges(6,1) / 0.4e1;
t643 = Ifges(5,4) / 0.4e1 + Ifges(6,4) / 0.4e1;
t642 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t640 = 0.3e1 / 0.4e1 * Ifges(5,5) + 0.3e1 / 0.4e1 * Ifges(6,5);
t639 = Ifges(5,2) / 0.4e1 + Ifges(6,2) / 0.4e1;
t637 = 0.3e1 / 0.4e1 * Ifges(5,6) + 0.3e1 / 0.4e1 * Ifges(6,6);
t636 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t632 = mrSges(6,3) * t738;
t631 = t324 * t801;
t630 = -t706 / 0.2e1;
t628 = t686 / 0.2e1;
t625 = -t676 / 0.2e1;
t617 = -t555 * t547 / 0.4e1;
t427 = Ifges(6,6) * t551 + t554 * t585;
t429 = Ifges(5,6) * t551 + t554 * t586;
t611 = t429 / 0.2e1 + t427 / 0.2e1;
t610 = -t432 / 0.4e1 - t430 / 0.4e1;
t431 = Ifges(6,5) * t551 + t554 * t587;
t433 = Ifges(5,5) * t551 + t554 * t588;
t609 = -t433 / 0.2e1 - t431 / 0.2e1;
t608 = t518 / 0.2e1 + t519 / 0.2e1;
t604 = m(6) * t661;
t602 = -m(5) * pkin(3) - mrSges(4,1) - t591;
t601 = mrSges(6,3) * pkin(4) - Ifges(5,5) - Ifges(6,5);
t412 = -pkin(8) * t665 + t550 * t526;
t541 = mrSges(6,1) * t665;
t594 = -mrSges(6,2) * t667 + t541;
t424 = -Ifges(6,3) * t554 + t551 * t583;
t425 = -Ifges(5,3) * t554 + t551 * t584;
t520 = Ifges(4,2) * t554 + t719;
t592 = t424 / 0.2e1 + t425 / 0.2e1 - t520 / 0.2e1;
t350 = -qJ(5) * t666 + t412;
t475 = t589 * t554;
t476 = t590 * t554;
t498 = -mrSges(6,2) * t551 - mrSges(6,3) * t666;
t499 = -mrSges(5,2) * t551 - mrSges(5,3) * t666;
t503 = mrSges(5,1) * t551 - mrSges(5,3) * t663;
t523 = t551 * Ifges(4,1) + t548;
t12 = pkin(2) * t515 - t324 * t502 - t325 * t500 - t342 * t498 - t350 * t496 - t407 * t503 - t408 * t499 - t411 * t501 - t412 * t497 - t504 * t475 - t505 * t473 - t476 * t728 - t554 * t835 - m(5) * (t407 * t411 + t408 * t412) - m(6) * (t324 * t325 + t342 * t350 + t504 * t505) + (t719 / 0.2e1 + t609 * t553 + t611 * t550 - t592) * t551 + (-t523 / 0.2e1 - t548 / 0.2e1 + (-t432 / 0.2e1 - t430 / 0.2e1 + t641 * t554) * t553 + (-t638 * t554 + t819) * t550 + (-m(5) * pkin(8) ^ 2 + Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1 + t636) * t551) * t554;
t556 = (-t520 / 0.4e1 + t425 / 0.4e1 + t424 / 0.4e1) * t466 + (t107 * t505 + t143 * t504 + t324 * t55 + t325 * t36 + t342 * t71 + t350 * t46) * t808 + (t407 * t89 + t408 * t90 + t411 * t65 + t412 * t66) * t810 + t282 * t786 + t474 * t788 + t350 * t790 + t476 * t791 + t475 * t794 + t223 * t776 + t283 * t777 + t222 * t785 + t415 * t515 / 0.2e1 + t65 * t503 / 0.2e1 + t46 * t498 / 0.2e1 + t66 * t499 / 0.2e1 + t412 * t221 / 0.2e1 + t408 * t281 / 0.2e1 + t342 * t280 / 0.2e1 - pkin(2) * t294 / 0.2e1 + t162 * t762 - t250 * t763 + t36 * t765 + t89 * t766 + t55 * t767 + t90 * t769 + t71 * t770 + t143 * t771 + (-t523 / 0.4e1 - t548 / 0.4e1 + t610 * t553 + t612 * t550) * t465 + (t431 / 0.4e1 + t433 / 0.4e1) * t338 + (t427 / 0.4e1 + t429 / 0.4e1) * t337;
t562 = t824 * t465 / 0.4e1 + t833 * t743 + t832 * t737;
t559 = t244 / 0.4e1 - t180 / 0.4e1 - t179 / 0.4e1 + t722 / 0.4e1 + t295 / 0.4e1 + (-t369 / 0.2e1 + t163 / 0.2e1 + m(5) * t791 - t690 / 0.2e1) * pkin(8) + t562;
t561 = (-t182 / 0.4e1 - t181 / 0.4e1) * t550 + (t184 / 0.4e1 + t183 / 0.4e1) * t553 + (Ifges(5,3) / 0.4e1 + Ifges(6,3) / 0.4e1 + Ifges(4,2) / 0.4e1) * t465 + (m(5) * t788 - t368 / 0.2e1 - t251 / 0.2e1 - t691 / 0.2e1) * pkin(8) + t110 / 0.4e1 + t111 / 0.4e1 - t243 / 0.4e1 + t296 / 0.4e1 - t720 / 0.4e1;
t2 = t556 + ((0.3e1 / 0.4e1 * Ifges(4,6) + t821) * t669 + t561) * t551 + (Ifges(4,5) * t623 + t559) * t554 + (-pkin(9) * t309 / 0.2e1 + t100 * t803 + t72 * t801 - t189 / 0.4e1 - t188 / 0.4e1) * t553 + (t99 * t802 + t64 * mrSges(6,3) / 0.2e1 + pkin(9) * t311 / 0.2e1 - t191 / 0.4e1 - t190 / 0.4e1) * t550 + (-t710 / 0.2e1 + t617) * t549 + t838;
t582 = t2 * qJD(1) - t12 * qJD(2);
t15 = -t324 * t652 - t341 * t496 + t408 * t501 - t473 * t657 - t504 * t594 + (-t497 - t651) * t407 + (t500 - t604) * t342 + (-t729 * t733 + t826 * t742 + t836 * t554 + (t695 + t702) * t553 + t827 * t738 + (-pkin(8) * t591 + t825 * t744 + (t521 + t522) * t738) * t551) * t551;
t557 = pkin(8) * t793 + t591 * t791 + mrSges(6,2) * t625 + (t734 / 0.2e1 + t162 / 0.2e1) * t729 + (t699 / 0.2e1 - t694 / 0.2e1) * mrSges(6,3) + (t65 * t742 + t66 * t740) * mrSges(5,3) - t833 * t553 / 0.4e1 + t830 * t737 + t821 * t465 - t823 * t338 - t607 * t337 + (t831 + t832) * t743;
t558 = (-t165 / 0.4e1 - t164 / 0.4e1) * t554 + (t407 * t803 - t610 + t631) * t337 + (-t695 / 0.2e1 - t702 / 0.2e1 + t817) * t338 + t541 * t794 + t341 * t790 - t342 * t222 / 0.2e1 + t221 * t777 - t408 * t223 / 0.2e1 + t45 * t770 + t46 * t768 + t160 * t763 + t65 * t769 - t66 * t501 / 0.2e1;
t566 = -pkin(4) * t310 / 0.2e1 + t100 * mrSges(5,2) / 0.2e1 + mrSges(6,1) * t799 + t72 * mrSges(6,2) / 0.2e1 + t99 * t806;
t4 = t558 - t641 * t406 - t638 * t405 + (-t796 / 0.2e1 - t709 / 0.2e1 + t703 / 0.2e1 - t701 / 0.2e1 + t700 / 0.2e1) * m(6) + (-t636 * t669 + t557) * t551 + t566;
t581 = t4 * qJD(1) - t15 * qJD(2);
t570 = t142 * t809 + t697 / 0.2e1 - t696 / 0.2e1;
t18 = t570 + t815;
t68 = (t550 * t496 - m(6) * (-t324 * t553 - t342 * t550) + t553 * t500) * t551;
t580 = qJD(1) * t18 - qJD(2) * t68;
t579 = t639 + t644;
t578 = t685 / 0.2e1 + t687 / 0.2e1;
t563 = t636 * t551 + (m(6) * t785 + t765) * pkin(4) + mrSges(6,1) * t785 + t350 * t804 + mrSges(5,1) * t776 + t412 * t805;
t10 = t541 * t746 - t672 / 0.2e1 - (t767 - t604 / 0.2e1) * t514 + (pkin(9) * t769 + mrSges(6,1) * t764 - t637 * t554 + (mrSges(6,2) * t745 + pkin(8) * t806 - t797 / 0.2e1 + (t655 - t579) * t550 + t839) * t551 - t817) * t550 + (pkin(9) * t766 + mrSges(6,2) * t764 + t640 * t554 + (t786 - t341 / 0.2e1) * mrSges(6,3) + (pkin(8) * t805 + t798 / 0.2e1 - t688 / 0.2e1 + (t655 + t579) * t553 + (Ifges(5,4) + Ifges(6,4)) * t550 + (t760 - t732 / 0.2e1) * pkin(4) + t607) * t551 + t610) * t553 + t563;
t29 = -t512 * t730 - t544 * t589 + (-t553 * t642 - t605 + t797) * t553 + (-pkin(4) * t732 + t798 + t642 * t550 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t553 + t608) * t550;
t565 = t636 * t466 + mrSges(6,1) * t800 + t71 * t804 + t89 * mrSges(5,1) / 0.2e1 + t90 * t805;
t569 = pkin(3) * t793 + t160 * t746 + t220 * t761;
t8 = (-t115 / 0.4e1 - t114 / 0.4e1 + mrSges(6,2) * t795 + mrSges(5,2) * t792 - t167 / 0.4e1 - t166 / 0.4e1) * t553 + (mrSges(6,1) * t795 + mrSges(5,1) * t792 - t169 / 0.4e1 - t168 / 0.4e1 + t113 / 0.4e1 + t112 / 0.4e1) * t550 + (t550 * t643 + t553 * t644 + t607) * t338 + (t550 * t639 - t553 * t643 - t823) * t337 - (-t634 / 0.2e1 + t789) * t514 + (t162 * t744 + t338 * t760 + t787 + (t625 - t675 / 0.2e1 + t800) * m(6)) * pkin(4) + (t223 * t738 + t221 * t742 + (-t668 / 0.2e1 + t338 * t738) * mrSges(5,3)) * pkin(9) + (t511 * t784 - t514 * t782 + (-t45 / 0.2e1 + t36 / 0.2e1) * t553) * mrSges(6,3) + (t550 * t637 - t553 * t640) * t465 + t565 + t569;
t573 = -t8 * qJD(1) - t10 * qJD(2) - t29 * qJD(3);
t568 = m(6) * (-t337 * t514 - t338 * t511 + t694 - t699);
t19 = (-t692 / 0.2e1 - t708 / 0.2e1 - t220 / 0.2e1) * t553 + (-t693 / 0.2e1 + t630 + t789) * t550 + t650 - t568 / 0.2e1;
t268 = m(6) * (-t511 * t550 - t514 * t553) + t818 * mrSges(6,3);
t567 = m(6) * ((t342 - t671) * t553 + (t514 * t551 - t324) * t550);
t34 = (t683 / 0.2e1 - t496 / 0.2e1) * t553 + (t684 / 0.2e1 + t767) * t550 + t649 - t567 / 0.2e1;
t572 = -qJD(1) * t19 - qJD(2) * t34 + qJD(3) * t268;
t393 = -m(6) * t657 - t594;
t477 = t550 * t635 + t685;
t81 = -t338 * t635 - t329;
t571 = qJD(1) * t81 + qJD(2) * t393 - qJD(3) * t477;
t35 = t496 * t738 + t567 / 0.2e1 + t500 * t744 + t649 + t578 * t554;
t20 = t337 * t632 + t338 * t628 + t568 / 0.2e1 + t222 * t744 + t220 * t738 + t650 - t578 * t465;
t17 = -t570 + t815;
t11 = t501 * t645 + t497 * t646 + t473 * t647 + t553 * t631 + t341 * t632 + t672 / 0.2e1 + t835 / 0.2e1 - t591 * t731 / 0.2e1 + (-t661 * t514 + (t504 * t550 + t544 * t665) * pkin(4)) * t808 + t589 * t763 + t594 * t745 + t500 * t759 + t648 * t665 + t628 * t671 + t563 - t605 * t667 + t818 * t727 * t803 - (-t688 + t825) * t665 / 0.2e1 + (t575 + t574 + t827) * t743 + (t577 + t576 + t826) * t737 + (-t814 - t824 / 0.4e1) * t554;
t7 = t221 * t646 + t162 * t647 - t514 * t630 + t45 * t632 + t565 - t569 + (m(6) * t800 + t787) * pkin(4) + t698 * t801 + (-t723 * t514 + (t675 + t676) * pkin(4)) * t808 + t590 * t791 + t589 * t794 + t814 * t465 + t222 * t759 + t562 + t655 * t668 + (t586 + t585) * t337 / 0.4e1 + (t588 + t587) * t338 / 0.4e1 + t830 * t741 + t831 * t737 + (t223 + t707) * t645 + (t648 - t607) * t338 + t839 * t337;
t3 = t558 + t557 * t551 + t712 / 0.2e1 + t711 / 0.2e1 + t714 / 0.2e1 + t713 / 0.2e1 - t566 + (t700 - t701 + t703 - t709 + t796) * t808 + (Ifges(5,3) + Ifges(6,3)) * t598;
t1 = t556 + t309 * t726 / 0.2e1 + t677 * t802 + t72 * t632 + Ifges(4,3) * t624 + t680 * t803 + t686 * t799 + t311 * t646 + t559 * t554 + t549 * t617 + Ifges(4,5) * t596 + Ifges(4,6) * t599 + t828 * t741 + t829 * t737 + (t561 + (Ifges(4,6) + t517 + t516) * t669 / 0.4e1) * t551 - t838;
t13 = [qJD(2) * t5 + qJD(3) * t6 + qJD(4) * t9 + qJD(5) * t23, t1 * qJD(3) + t3 * qJD(4) + t17 * qJD(5) + t689 + (t816 + t536 + m(4) * (-pkin(2) * t470 + (-t287 * t551 + t288 * t554) * pkin(8)) + 0.2e1 * (t100 * t408 + t240 * t728 + t407 * t99) * t810 + 0.2e1 * (t142 * t504 + t324 * t64 + t342 * t72) * t808 + t504 * t230 + t72 * t496 + t100 * t497 + t64 * t500 + t99 * t501 + t240 * t474 + t142 * t473 + t407 * t311 + t408 * t309 - pkin(2) * t409 + t342 * t308 + t324 * t310 + t826 * t778 + t819 * t405 + (pkin(8) * t442 + t288 * mrSges(4,3) - t470 * mrSges(4,1) + t339 / 0.2e1 - t186 / 0.2e1 - t187 / 0.2e1) * t554 + (-t287 * mrSges(4,3) + t470 * mrSges(4,2) + t780 + (t190 / 0.2e1 + t191 / 0.2e1) * t553 + (-t188 / 0.2e1 - t189 / 0.2e1) * t550 + (-t443 + t231) * pkin(8)) * t551 + ((Ifges(4,5) * t551 / 0.2e1 + Ifges(4,6) * t735 - Ifges(3,6)) * t552 + (t523 * t735 + t551 * t592) * t555) * t549) * qJD(2), t682 + t1 * qJD(2) + t7 * qJD(4) + t20 * qJD(5) + (t90 * mrSges(5,3) + t71 * mrSges(6,3) - t605 * t465 + (m(5) * t90 + t281) * pkin(9) + t614) * t658 + (-t89 * mrSges(5,3) - t55 * mrSges(6,3) + t608 * t465 + (-m(5) * t89 - t283) * pkin(9) + t613) * t659 + (m(6) * (t143 * t544 + t511 * t55 - t514 * t71) - t544 * t250 + t511 * t282 + t143 * t512 - t514 * t280 + pkin(3) * t251 - t228 * mrSges(4,2) + t660 + t602 * t229 - t836 * t466) * qJD(3), t681 + t3 * qJD(2) + t7 * qJD(3) + (-mrSges(5,1) * t66 - mrSges(6,1) * t46 - mrSges(5,2) * t65 - mrSges(6,2) * t45 + (-m(6) * t46 - t708) * pkin(4) + t165 + t164) * qJD(4), qJD(2) * t17 + qJD(3) * t20 + t678; qJD(3) * t2 + qJD(4) * t4 + qJD(5) * t18 - t689, -qJD(3) * t12 - qJD(4) * t15 - qJD(5) * t68, t11 * qJD(4) + t35 * qJD(5) + (t412 * mrSges(5,3) + t350 * mrSges(6,3) + t605 * t554 + (m(5) * t412 + t499) * pkin(9) + t611) * t658 + (-t411 * mrSges(5,3) - t325 * mrSges(6,3) - t608 * t554 + (-m(5) * t411 - t503) * pkin(9) - t609) * t659 + t582 + (t547 + m(6) * (t325 * t511 - t350 * t514 + t505 * t544) + t544 * t475 + t511 * t502 + t505 * t512 - t514 * t498 - pkin(3) * t476 + (-Ifges(4,6) - t836) * t551 + (t551 * mrSges(4,2) + t554 * t602) * pkin(8)) * qJD(3), t11 * qJD(3) + t581 + (-mrSges(5,1) * t408 - mrSges(5,2) * t407 - mrSges(6,2) * t341 + (t550 * t601 + t553 * t725) * t551 - t635 * t342) * qJD(4), qJD(3) * t35 + t580; -qJD(2) * t2 - qJD(4) * t8 - qJD(5) * t19 - t682, -qJD(4) * t10 - qJD(5) * t34 - t582, -qJD(4) * t29 + qJD(5) * t268, t573 + (-mrSges(6,2) * t511 + (mrSges(5,2) * pkin(9) + t725) * t550 + (-mrSges(5,1) * pkin(9) - t601) * t553 + t635 * t514) * qJD(4), t572; -qJD(2) * t4 + qJD(3) * t8 + qJD(5) * t81 - t681, qJD(3) * t10 + qJD(5) * t393 - t581, -qJD(5) * t477 - t573, 0, t571; -qJD(2) * t18 + qJD(3) * t19 - qJD(4) * t81 - t678, qJD(3) * t34 - qJD(4) * t393 - t580, qJD(4) * t477 - t572, -t571, 0;];
Cq = t13;
