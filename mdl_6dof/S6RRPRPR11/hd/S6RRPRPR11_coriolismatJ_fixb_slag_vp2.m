% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPR11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:49
% EndTime: 2019-03-09 11:12:15
% DurationCPUTime: 15.33s
% Computational Cost: add. (29170->731), mult. (55494->979), div. (0->0), fcn. (58834->8), ass. (0->377)
t489 = sin(pkin(10));
t491 = sin(qJ(4));
t494 = cos(qJ(4));
t622 = cos(pkin(10));
t557 = t622 * t494;
t438 = t489 * t491 - t557;
t490 = sin(qJ(6));
t493 = cos(qJ(6));
t520 = t489 * t494 + t491 * t622;
t338 = t438 * t493 + t490 * t520;
t550 = t438 * t490 - t493 * t520;
t599 = Ifges(7,5) * t550 + Ifges(7,6) * t338;
t496 = -pkin(2) - pkin(8);
t605 = t491 * t496;
t447 = -t491 * qJ(5) + t605;
t448 = (-qJ(5) + t496) * t494;
t360 = t622 * t447 + t489 * t448;
t273 = -pkin(9) * t520 + t360;
t757 = -t489 * t447 + t622 * t448;
t777 = pkin(9) * t438 + t757;
t156 = t273 * t493 + t490 * t777;
t789 = -t273 * t490 + t493 * t777;
t819 = -t156 * mrSges(7,1) - t789 * mrSges(7,2);
t24 = t599 + t819;
t821 = t24 * qJD(6);
t495 = cos(qJ(2));
t492 = sin(qJ(2));
t621 = qJ(3) * t492;
t434 = t495 * t496 - pkin(1) - t621;
t740 = pkin(3) + pkin(7);
t459 = t740 * t492;
t349 = -t434 * t491 + t494 * t459;
t606 = t491 * t495;
t303 = qJ(5) * t606 + t349;
t285 = pkin(4) * t492 + t303;
t350 = t434 * t494 + t459 * t491;
t601 = t494 * t495;
t304 = -qJ(5) * t601 + t350;
t286 = t489 * t304;
t165 = t622 * t285 - t286;
t406 = t520 * t495;
t682 = pkin(9) * t406;
t128 = pkin(5) * t492 + t165 + t682;
t559 = t622 * t304;
t166 = t489 * t285 + t559;
t405 = t438 * t495;
t683 = pkin(9) * t405;
t136 = t166 + t683;
t63 = t128 * t493 - t136 * t490;
t177 = -t303 * t489 - t559;
t140 = t177 - t683;
t178 = t622 * t303 - t286;
t141 = t178 + t682;
t73 = t140 * t490 + t141 * t493;
t820 = t73 - t63;
t818 = t789 / 0.2e1;
t485 = t491 * pkin(4);
t479 = qJ(3) + t485;
t388 = pkin(5) * t520 + t479;
t817 = m(7) * t388;
t64 = t128 * t490 + t136 * t493;
t72 = t140 * t493 - t141 * t490;
t814 = t64 + t72;
t552 = t493 * t405 + t406 * t490;
t281 = t405 * t490 - t406 * t493;
t649 = t281 * Ifges(7,4);
t151 = Ifges(7,2) * t552 + t492 * Ifges(7,6) + t649;
t261 = Ifges(7,4) * t552;
t153 = t281 * Ifges(7,1) + t492 * Ifges(7,5) + t261;
t641 = Ifges(7,4) * t338;
t198 = Ifges(7,2) * t550 - t641;
t768 = mrSges(7,3) * t552;
t236 = -mrSges(7,2) * t492 + t768;
t460 = t740 * t495;
t422 = pkin(4) * t601 + t460;
t323 = -pkin(5) * t405 + t422;
t658 = mrSges(7,3) * t281;
t238 = mrSges(7,1) * t492 - t658;
t731 = -t238 / 0.2e1;
t735 = -t156 / 0.2e1;
t776 = -Ifges(7,2) * t281 + t261;
t782 = t281 * mrSges(7,1);
t796 = t552 * mrSges(7,2) + t782;
t781 = t338 * mrSges(7,1);
t797 = t550 * mrSges(7,2) - t781;
t798 = Ifges(7,1) * t552 - t649;
t799 = Ifges(7,1) * t550 + t641;
t815 = t156 * t731 + (-t798 / 0.4e1 + t151 / 0.4e1) * t338 + (t776 / 0.4e1 + t153 / 0.4e1) * t550 + t236 * t818 + t796 * t388 / 0.2e1 + t797 * t323 / 0.2e1 + (t799 / 0.4e1 + mrSges(7,3) * t735 - t198 / 0.4e1) * t281;
t801 = t550 * mrSges(7,1);
t582 = t801 / 0.2e1;
t812 = t323 * t796;
t811 = t388 * t797;
t810 = t165 - t178;
t553 = t492 * pkin(2) - qJ(3) * t495;
t437 = pkin(8) * t492 + t553;
t442 = t494 * t460;
t290 = pkin(4) * t495 + t442 + (-qJ(5) * t492 - t437) * t491;
t362 = t494 * t437 + t491 * t460;
t604 = t492 * t494;
t307 = qJ(5) * t604 + t362;
t175 = t622 * t290 - t307 * t489;
t176 = t489 * t290 + t622 * t307;
t607 = t491 * t492;
t404 = -t489 * t607 + t492 * t557;
t407 = t520 * t492;
t276 = t404 * t493 - t407 * t490;
t235 = -mrSges(7,2) * t495 + mrSges(7,3) * t276;
t551 = t404 * t490 + t493 * t407;
t237 = mrSges(7,1) * t495 - mrSges(7,3) * t551;
t359 = -t437 * t491 + t442;
t368 = -mrSges(6,2) * t495 + mrSges(6,3) * t404;
t370 = mrSges(6,1) * t495 - mrSges(6,3) * t407;
t573 = t622 * pkin(4);
t478 = t573 + pkin(5);
t685 = pkin(4) * t489;
t409 = t478 * t493 - t490 * t685;
t410 = t478 * t490 + t493 * t685;
t705 = t407 / 0.2e1;
t711 = t404 / 0.2e1;
t527 = Ifges(6,5) * t705 + Ifges(6,6) * t711;
t727 = t551 / 0.2e1;
t730 = t276 / 0.2e1;
t525 = Ifges(7,5) * t727 + Ifges(7,6) * t730;
t129 = pkin(5) * t495 - pkin(9) * t407 + t175;
t137 = pkin(9) * t404 + t176;
t67 = t129 * t493 - t137 * t490;
t68 = t129 * t490 + t137 * t493;
t804 = t495 / 0.2e1;
t546 = Ifges(7,3) * t804 - t68 * mrSges(7,2) / 0.2e1 + t67 * mrSges(7,1) / 0.2e1 + t525;
t567 = t604 / 0.2e1;
t568 = t607 / 0.2e1;
t741 = m(6) * pkin(4);
t592 = t741 / 0.2e1;
t703 = t410 / 0.2e1;
t704 = t409 / 0.2e1;
t742 = m(7) / 0.2e1;
t769 = Ifges(5,3) + Ifges(6,3);
t809 = (t409 * t67 + t410 * t68) * t742 + t175 * mrSges(6,1) / 0.2e1 - t176 * mrSges(6,2) / 0.2e1 + t359 * mrSges(5,1) / 0.2e1 - t362 * mrSges(5,2) / 0.2e1 + t237 * t704 + t235 * t703 + (t175 * t622 + t176 * t489) * t592 + Ifges(5,5) * t568 + Ifges(5,6) * t567 + t368 * t685 / 0.2e1 + t370 * t573 / 0.2e1 + t527 + t769 * t804 + t546;
t780 = t338 * mrSges(7,2);
t808 = t801 + t780;
t807 = -t798 / 0.2e1 + t151 / 0.2e1;
t803 = -t550 / 0.2e1;
t802 = t550 / 0.2e1;
t581 = -t781 / 0.2e1;
t583 = t782 / 0.2e1;
t686 = m(6) * t479;
t785 = t338 / 0.2e1;
t800 = t236 * t785;
t758 = t177 + t166;
t764 = Ifges(7,5) * t552;
t784 = Ifges(7,6) * t281;
t600 = t764 - t784;
t330 = Ifges(7,4) * t550;
t201 = -Ifges(7,1) * t338 + t330;
t775 = Ifges(7,2) * t338 + t330;
t795 = t775 / 0.4e1 + t201 / 0.4e1;
t794 = t338 ^ 2 + t550 ^ 2;
t526 = t764 / 0.2e1 - t784 / 0.2e1;
t792 = t338 * t67 + t550 * t68;
t788 = 0.2e1 * mrSges(7,2);
t786 = -t338 / 0.2e1;
t778 = t802 * t658;
t774 = t338 * t410 - t409 * t550;
t635 = t406 * mrSges(6,3);
t371 = mrSges(6,1) * t492 + t635;
t771 = t371 / 0.2e1;
t728 = -t552 / 0.2e1;
t729 = t552 / 0.2e1;
t665 = Ifges(5,6) * t494;
t669 = Ifges(5,5) * t491;
t761 = Ifges(3,4) + Ifges(4,6) - t669 / 0.2e1 - t665 / 0.2e1;
t756 = Ifges(6,5) * t405 + Ifges(6,6) * t406 + t600;
t392 = Ifges(6,4) * t405;
t268 = -t406 * Ifges(6,1) + t492 * Ifges(6,5) + t392;
t755 = Ifges(6,2) * t406 + t268 + t392;
t706 = t406 / 0.2e1;
t710 = -t405 / 0.2e1;
t754 = t360 * t706 + t710 * t757;
t753 = -Ifges(6,5) * t520 + Ifges(6,6) * t438 + t599;
t752 = t175 * t438 - t176 * t520;
t751 = (t491 ^ 2 + t494 ^ 2) * t495;
t650 = t281 * mrSges(7,2);
t656 = t552 * mrSges(7,1);
t160 = t650 - t656;
t670 = Ifges(6,4) * t406;
t266 = t405 * Ifges(6,2) + t492 * Ifges(6,6) - t670;
t351 = -t438 * mrSges(6,1) - mrSges(6,2) * t520;
t427 = Ifges(6,4) * t520;
t353 = Ifges(6,2) * t438 - t427;
t632 = t438 * Ifges(6,4);
t354 = -Ifges(6,2) * t520 - t632;
t355 = -Ifges(6,1) * t520 + t632;
t356 = -t438 * Ifges(6,1) - t427;
t590 = pkin(4) * t606;
t365 = -pkin(5) * t406 - t590;
t637 = t405 * mrSges(6,3);
t369 = -mrSges(6,2) * t492 + t637;
t684 = pkin(4) * t494;
t393 = -pkin(5) * t438 + t684;
t627 = t494 * mrSges(5,1);
t630 = t491 * mrSges(5,2);
t452 = t627 - t630;
t543 = Ifges(6,1) * t405 + t670;
t556 = -t406 * mrSges(6,1) + t405 * mrSges(6,2);
t691 = t492 / 0.4e1;
t695 = t460 / 0.2e1;
t712 = t393 / 0.2e1;
t714 = t365 / 0.2e1;
t748 = (-t543 / 0.4e1 + t266 / 0.4e1) * t438 + (t354 / 0.4e1 - t355 / 0.4e1) * t406 + (t356 + t353) * t405 / 0.4e1 + t479 * t556 / 0.2e1 - t360 * t771 + t795 * t552 - t808 * t714 + t160 * t712 + t452 * t695 + (t156 * t820 + t323 * t393 + t365 * t388 + t814 * t789) * t742 + (t63 * t803 + t789 * t728 + t73 * t802 + t785 * t814) * mrSges(7,3) + t422 * t351 / 0.2e1 + t757 * t369 / 0.2e1 - t755 * t520 / 0.4e1 + t753 * t691 + t815;
t745 = -m(6) / 0.2e1;
t744 = m(6) / 0.2e1;
t743 = -m(7) / 0.2e1;
t739 = -qJ(3) / 0.2e1;
t736 = -t789 / 0.2e1;
t734 = -t165 / 0.2e1;
t726 = t281 / 0.2e1;
t709 = t405 / 0.2e1;
t707 = -t406 / 0.2e1;
t702 = -t438 / 0.2e1;
t701 = t438 / 0.2e1;
t700 = -t520 / 0.2e1;
t672 = Ifges(5,4) * t491;
t541 = t494 * Ifges(5,2) + t672;
t698 = t541 / 0.4e1;
t671 = Ifges(5,4) * t494;
t544 = Ifges(5,1) * t491 + t671;
t697 = t544 / 0.4e1;
t458 = Ifges(5,1) * t494 - t672;
t696 = -t458 / 0.2e1;
t694 = -t491 / 0.2e1;
t693 = t491 / 0.2e1;
t690 = -t494 / 0.2e1;
t689 = t494 / 0.2e1;
t687 = t496 / 0.2e1;
t680 = t63 * mrSges(7,2);
t679 = t64 * mrSges(7,1);
t676 = t72 * mrSges(7,1);
t675 = t73 * mrSges(7,2);
t668 = Ifges(5,5) * t494;
t666 = Ifges(5,6) * t491;
t657 = t276 * mrSges(7,1);
t652 = t551 * mrSges(7,2);
t150 = Ifges(7,4) * t551 + Ifges(7,2) * t276 + Ifges(7,6) * t495;
t152 = Ifges(7,1) * t551 + Ifges(7,4) * t276 + Ifges(7,5) * t495;
t159 = t652 - t657;
t265 = Ifges(6,4) * t407 + Ifges(6,2) * t404 + Ifges(6,6) * t495;
t267 = Ifges(6,1) * t407 + Ifges(6,4) * t404 + Ifges(6,5) * t495;
t634 = t407 * mrSges(6,2);
t639 = t404 * mrSges(6,1);
t291 = t634 - t639;
t636 = t406 * mrSges(6,2);
t638 = t405 * mrSges(6,1);
t292 = -t636 - t638;
t421 = (-t684 - t740) * t492;
t322 = -pkin(5) * t404 + t421;
t399 = Ifges(5,6) * t495 + t492 * t541;
t401 = Ifges(5,5) * t495 + t492 * t544;
t428 = t452 * t492;
t625 = t495 * mrSges(5,1);
t443 = -mrSges(5,3) * t607 + t625;
t589 = mrSges(5,3) * t606;
t444 = mrSges(5,1) * t492 + t589;
t624 = t495 * mrSges(5,2);
t445 = mrSges(5,3) * t604 - t624;
t588 = mrSges(5,3) * t601;
t446 = -mrSges(5,2) * t492 - t588;
t536 = -pkin(2) * t495 - t621;
t451 = -pkin(1) + t536;
t453 = t495 * mrSges(4,2) - t492 * mrSges(4,3);
t561 = m(4) * t451 + t453;
t628 = t492 * Ifges(5,6);
t400 = -t495 * t541 + t628;
t603 = t494 * t400;
t629 = t492 * Ifges(5,5);
t402 = -t495 * t544 + t629;
t609 = t491 * t402;
t3 = t322 * t160 + t323 * t159 + t64 * t235 + t68 * t236 + t63 * t237 + t67 * t238 + t153 * t727 + t150 * t729 + t151 * t730 + t152 * t726 + t268 * t705 + t267 * t707 + t265 * t709 + t266 * t711 + m(7) * (t322 * t323 + t63 * t67 + t64 * t68) + m(6) * (t165 * t175 + t166 * t176 + t421 * t422) + m(5) * (t349 * t359 + t350 * t362 - t459 * t460) + (-pkin(1) * mrSges(3,2) - t451 * mrSges(4,3) - t459 * t452 + t401 * t694 + t399 * t690 + Ifges(7,5) * t726 + Ifges(7,6) * t729 + Ifges(6,5) * t707 + Ifges(6,6) * t709 + (Ifges(3,1) - Ifges(3,2) + Ifges(7,3) - Ifges(4,3) + Ifges(4,2) + t769) * t492 + t761 * t495) * t495 + t166 * t368 + t176 * t369 + t165 * t370 + t175 * t371 + t421 * t292 + t422 * t291 + t349 * t443 + t359 * t444 + t350 * t445 + t362 * t446 - t460 * t428 + t561 * t553 + (-pkin(1) * mrSges(3,1) - t451 * mrSges(4,2) + t609 / 0.2e1 + t603 / 0.2e1 - t761 * t492 + t525 + t527) * t492;
t648 = t3 * qJD(1);
t456 = -Ifges(5,2) * t491 + t671;
t626 = t494 * mrSges(5,2);
t545 = -t491 * mrSges(5,1) - t626;
t522 = t545 * t495;
t4 = (t402 * t689 + t400 * t694 + (t456 * t690 + t458 * t694) * t495 + (m(6) * t422 + t292) * t485) * t495 - t812 - t422 * t556 - m(6) * (t165 * t177 + t166 * t178) - t73 * t236 - t72 * t238 + (-t588 - t446) * t349 - t460 * t522 + (-t589 + t444) * t350 + t543 * t706 + t266 * t707 - m(7) * (t323 * t365 + t63 * t72 + t64 * t73) + t63 * t768 + t165 * t637 + t64 * t658 - t166 * t635 - t365 * t160 - t178 * t369 - t177 * t371 + t755 * t710 - ((t666 - t668) * t495 + t756) * t492 / 0.2e1 + (t776 + t153) * t728 + t807 * t281;
t640 = t4 * qJD(1);
t633 = t438 * mrSges(6,3);
t631 = t520 * mrSges(6,3);
t7 = -t64 * t238 + t63 * t236 + t492 * t600 / 0.2e1 + t812 + (-t64 * mrSges(7,3) - t807) * t281 + (-t63 * mrSges(7,3) + t153 / 0.2e1 + t776 / 0.2e1) * t552;
t623 = t7 * qJD(1);
t523 = t238 * t803 + t800;
t524 = mrSges(7,1) * t727 + mrSges(7,2) * t730;
t17 = (t281 * t803 + t552 * t786) * mrSges(7,3) + t523 + t524;
t620 = t17 * qJD(1);
t602 = t494 * t446;
t608 = t491 * t444;
t21 = -t276 * t236 + t551 * t238 - t404 * t369 + t407 * t371 + m(7) * (-t276 * t64 + t551 * t63) + m(6) * (t165 * t407 - t166 * t404) + (-t602 + m(5) * (t349 * t491 - t350 * t494) + t608 - t561) * t492;
t617 = t21 * qJD(1);
t616 = t404 * t520;
t615 = t405 * t520;
t614 = t406 * t438;
t613 = t407 * t438;
t612 = t410 * t281;
t611 = t422 * t494;
t127 = t410 * mrSges(7,1) + mrSges(7,2) * t409;
t594 = t127 * qJD(6);
t593 = mrSges(6,3) * t685;
t591 = t745 + t743;
t587 = t684 / 0.2e1;
t585 = -t768 / 0.2e1;
t579 = -t631 / 0.2e1;
t578 = mrSges(5,3) * t687;
t566 = -t496 * t444 / 0.2e1;
t565 = t446 * t687;
t562 = t438 ^ 2 + t520 ^ 2;
t352 = mrSges(6,1) * t520 - t438 * mrSges(6,2);
t549 = -t590 / 0.2e1;
t548 = mrSges(6,3) * t573;
t539 = -t665 - t669;
t22 = t811 - (t799 / 0.2e1 - t198 / 0.2e1) * t338 + (t201 / 0.2e1 + t775 / 0.2e1) * t550;
t504 = (mrSges(7,3) * t736 + t795) * t552 + t599 * t691 + t815;
t6 = t504 - t546;
t535 = t6 * qJD(1) + t22 * qJD(2);
t502 = (-t276 * t802 + t551 * t785) * mrSges(7,3) + (t613 / 0.2e1 + t616 / 0.2e1) * mrSges(6,3) + m(5) * t695 + (-t360 * t404 + t407 * t757 + t422) * t744 + (-t156 * t276 + t551 * t789 + t323) * t742 - t656 / 0.2e1 + t650 / 0.2e1 - t638 / 0.2e1 - t636 / 0.2e1;
t528 = t359 * t494 + t362 * t491;
t506 = -m(5) * t528 / 0.2e1 - t752 * t745 - t792 * t743 + t237 * t785 + t235 * t802 + t370 * t701 + t368 * t700;
t14 = (-t443 / 0.2e1 + t625 / 0.2e1) * t494 + (-t445 / 0.2e1 - t624 / 0.2e1) * t491 + t502 + t506;
t512 = t808 - t352 + t545;
t85 = mrSges(4,3) + (m(5) + m(4)) * qJ(3) + t686 + t817 - t512;
t534 = qJD(1) * t14 + qJD(2) * t85;
t503 = (-t281 * t785 + t552 * t802) * mrSges(7,3) + (t614 / 0.2e1 - t615 / 0.2e1) * mrSges(6,3) + (t165 * t438 - t166 * t520 + t360 * t405 + t406 * t757) * t744 + (t156 * t552 - t281 * t789 + t338 * t63 + t550 * t64) * t742 + t238 * t785 + t236 * t802 + t371 * t701 + t369 * t700;
t507 = t421 * t744 + t322 * t742 - t657 / 0.2e1 + t652 / 0.2e1 - t639 / 0.2e1 + t634 / 0.2e1;
t16 = t503 - t507;
t35 = t794 * mrSges(7,3) + t562 * mrSges(6,3) + m(7) * (t156 * t550 + t338 * t789) + m(6) * (-t360 * t520 + t438 * t757);
t533 = -qJD(1) * t16 - qJD(2) * t35;
t509 = (-t281 * t409 + t410 * t552) * t742 + (t405 * t489 + t406 * t622) * t592;
t518 = m(6) * t549 + m(7) * t714;
t46 = -t509 + t518 + t796 + t556;
t508 = (t338 * t409 + t410 * t550) * t742 + (t438 * t622 - t489 * t520) * t592;
t521 = m(6) * t587 + m(7) * t712;
t49 = t797 + t351 - t508 + t521;
t532 = qJD(1) * t46 + qJD(2) * t49;
t60 = (t281 * t338 - t550 * t552) * t742 + (-t614 + t615) * t744;
t511 = t562 * t745 - t742 * t794;
t77 = t511 + t591;
t531 = qJD(1) * t60 + qJD(2) * t77;
t79 = t729 * t788 + 0.2e1 * t583;
t90 = t788 * t802 + 0.2e1 * t581;
t530 = qJD(1) * t79 + qJD(2) * t90;
t26 = t552 * t236 - t281 * t238 + t405 * t369 + t406 * t371 + m(7) * (-t281 * t63 + t552 * t64) + m(6) * (t165 * t406 + t166 * t405);
t529 = qJD(1) * t26 + qJD(3) * t60;
t12 = t198 * t785 + t799 * t786 + t811 + qJ(3) * t452 + t479 * t351 + (t541 / 0.2e1 + t696) * t491 - (t353 / 0.2e1 + t356 / 0.2e1) * t520 + (t354 / 0.2e1 - t355 / 0.2e1) * t438 + (-t456 / 0.2e1 - t544 / 0.2e1 + pkin(4) * t352) * t494 + t684 * t686 + (t775 + t201) * t802 + (-t808 + t817) * t393;
t514 = -t360 * t810 + t758 * t757;
t2 = -t609 / 0.4e1 - t603 / 0.4e1 + t494 * t565 + t491 * t566 + qJ(3) * t522 / 0.2e1 + t352 * t549 - t631 * t734 + t539 * t691 + t606 * t697 + t456 * t606 / 0.2e1 + ((-t479 * t606 + t611) * pkin(4) + t514) * t744 + t178 * t579 + t292 * t587 + (t698 + t696) * t601 + t758 * t633 / 0.2e1 + t578 * t751 + t754 * mrSges(6,3) + t748 - t809;
t517 = t2 * qJD(1) + t12 * qJD(2);
t500 = (-t758 * t438 - t520 * t810) * t745 + (-t814 * t338 - t550 * t820) * t743 + t550 * t731 + t800 + t369 * t701 + t520 * t771 + t608 / 0.2e1 - t602 / 0.2e1 + t768 * t786 - t778 + t633 * t710 + t406 * t579 - mrSges(5,3) * t751 / 0.2e1;
t505 = (-t276 * t410 + t409 * t551) * t742 + mrSges(6,2) * t711 + mrSges(6,1) * t705 + (-t404 * t489 + t407 * t622) * t592 + mrSges(5,1) * t568 + mrSges(5,2) * t567 + t524;
t8 = t500 + t505;
t516 = t8 * qJD(1);
t515 = -t409 * t236 / 0.2e1 + t238 * t703 - t526;
t10 = (t552 * t704 + t612 / 0.2e1) * mrSges(7,3) + (-t73 / 0.2e1 + t63 / 0.2e1) * mrSges(7,2) + (t72 / 0.2e1 + t64 / 0.2e1) * mrSges(7,1) + t515 + t526;
t23 = (t736 + t818) * mrSges(7,2) + (t735 + t156 / 0.2e1) * mrSges(7,1);
t44 = t582 - t801 / 0.2e1 + (t785 + t786) * mrSges(7,2);
t513 = t10 * qJD(1) - t23 * qJD(2) - t44 * qJD(3) + t127 * qJD(4);
t91 = t581 + t781 / 0.2e1;
t83 = t508 + t521;
t80 = t583 - t782 / 0.2e1;
t76 = t511 - t591;
t69 = t509 + t518;
t59 = t60 * qJD(5);
t45 = t780 + 0.2e1 * t582;
t18 = -t338 * t585 - t523 + t524 + t778;
t15 = t503 + t507;
t13 = (m(4) * pkin(7) + t627 / 0.2e1 - t630 / 0.2e1 + mrSges(4,1)) * t495 + t502 + t445 * t693 + t443 * t689 - t506;
t11 = -t679 / 0.2e1 - t680 / 0.2e1 + t409 * t585 - mrSges(7,3) * t612 / 0.2e1 + t676 / 0.2e1 - t675 / 0.2e1 - t515 + t526;
t9 = -t500 + t505;
t5 = t504 + t546;
t1 = (-(t178 / 0.2e1 + t734) * t520 + (t166 / 0.2e1 + t177 / 0.2e1) * t438 + t754) * mrSges(6,3) + ((mrSges(5,2) * t739 - t458 / 0.4e1 + t698 + (t578 - Ifges(5,1) / 0.4e1) * t494) * t494 + (mrSges(5,1) * t739 + t671 / 0.2e1 + t697 + t456 / 0.4e1 + (t578 - Ifges(5,2) / 0.4e1) * t491 + (-t352 / 0.2e1 - t686 / 0.2e1) * pkin(4)) * t491) * t495 + (pkin(4) * t611 + t514) * t744 + (-t402 / 0.4e1 - t629 / 0.4e1 + t566) * t491 + (-t628 / 0.4e1 - t400 / 0.4e1 + pkin(4) * t292 / 0.2e1 + t565) * t494 + t748 + t809;
t19 = [qJD(2) * t3 + qJD(3) * t21 - qJD(4) * t4 + qJD(5) * t26 + qJD(6) * t7, t13 * qJD(3) + t1 * qJD(4) + t15 * qJD(5) + t5 * qJD(6) + t648 + (t156 * t235 + t152 * t786 + m(5) * (-qJ(3) * t459 + t496 * t528) + t792 * mrSges(7,3) + (-t399 / 0.2e1 - t459 * mrSges(5,1) - t362 * mrSges(5,3) + t496 * t445) * t491 + t198 * t730 + t201 * t727 + t265 * t700 + t267 * t702 + t356 * t705 + t354 * t711 + (-qJ(3) * mrSges(4,1) + t456 * t689 + t458 * t693 + Ifges(4,5) - Ifges(3,6)) * t492 + (m(4) * t536 - t495 * mrSges(3,1) + t492 * mrSges(3,2) + t453) * pkin(7) + (t401 / 0.2e1 - t459 * mrSges(5,2) - t359 * mrSges(5,3) + t496 * t443) * t494 + (-Ifges(4,4) + Ifges(3,5) + Ifges(7,5) * t786 + Ifges(7,6) * t802 + Ifges(6,5) * t702 + Ifges(6,6) * t700 + t668 / 0.2e1 - t666 / 0.2e1 - pkin(2) * mrSges(4,1)) * t495 + t150 * t802 + t789 * t237 + 0.2e1 * (t156 * t68 + t322 * t388 + t67 * t789) * t742 - t322 * t808 + t360 * t368 + t388 * t159 + t421 * t352 - qJ(3) * t428 + t479 * t291 + 0.2e1 * (t175 * t757 + t176 * t360 + t421 * t479) * t744 + t757 * t370 + t752 * mrSges(6,3)) * qJD(2), t617 + t13 * qJD(2) + 0.2e1 * ((t276 * t550 - t338 * t551) * t742 + (-t613 - t616) * t744) * qJD(3) + t9 * qJD(4) + t59 + t18 * qJD(6), -t640 + t1 * qJD(2) + t9 * qJD(3) + (-t405 * t548 + t406 * t593 + t676 - t675 - t409 * t768 - t410 * t658 + m(7) * (t409 * t72 + t410 * t73) + t177 * mrSges(6,1) - t178 * mrSges(6,2) + (t177 * t622 + t178 * t489) * t741 - t349 * mrSges(5,2) - t350 * mrSges(5,1) - Ifges(5,5) * t601 + Ifges(5,6) * t606 + t756) * qJD(4) + t69 * qJD(5) + t11 * qJD(6), qJD(2) * t15 + qJD(4) * t69 + qJD(6) * t80 + t529, t623 + t5 * qJD(2) + t18 * qJD(3) + t11 * qJD(4) + t80 * qJD(5) + (t600 - t679 - t680) * qJD(6); qJD(3) * t14 + qJD(4) * t2 + qJD(5) * t16 + qJD(6) * t6 - t648, qJD(3) * t85 + qJD(4) * t12 + qJD(5) * t35 + qJD(6) * t22, qJD(5) * t76 + t534 (m(7) * (-t156 * t409 + t410 * t789) + (-t360 * t622 + t489 * t757) * t741 - t757 * mrSges(6,2) - t360 * mrSges(6,1) + t438 * t593 + t520 * t548 - t496 * t626 - mrSges(5,1) * t605 + t539 + t774 * mrSges(7,3) + t753 + t819) * qJD(4) + t83 * qJD(5) + t821 + t517, qJD(3) * t76 + qJD(4) * t83 + qJD(6) * t91 - t533, t24 * qJD(4) + t91 * qJD(5) + t535 + t821; -qJD(2) * t14 - qJD(4) * t8 - qJD(6) * t17 + t59 - t617, qJD(5) * t77 - t534, 0 ((-t438 * t489 - t520 * t622) * t741 - m(7) * t774 + t512) * qJD(4) + t45 * qJD(6) - t516, t531, t45 * qJD(4) + qJD(6) * t808 - t620; -qJD(2) * t2 + qJD(3) * t8 - qJD(5) * t46 - qJD(6) * t10 + t640, -qJD(5) * t49 + qJD(6) * t23 - t517, qJD(6) * t44 + t516, -t594, -t532, -t513 - t594; -qJD(2) * t16 + qJD(4) * t46 + qJD(6) * t79 - t529, -qJD(3) * t77 + qJD(4) * t49 + qJD(6) * t90 + t533, -t531, t532, 0, t530; -qJD(2) * t6 + qJD(3) * t17 + qJD(4) * t10 - qJD(5) * t79 - t623, -qJD(4) * t23 - qJD(5) * t90 - t535, -t44 * qJD(4) + t620, t513, -t530, 0;];
Cq  = t19;
