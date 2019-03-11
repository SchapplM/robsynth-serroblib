% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:32:50
% EndTime: 2019-03-09 18:34:30
% DurationCPUTime: 55.89s
% Computational Cost: add. (39237->1092), mult. (94662->1463), div. (0->0), fcn. (77887->18), ass. (0->498)
t442 = sin(pkin(12));
t444 = cos(pkin(12));
t449 = sin(qJ(3));
t454 = cos(qJ(3));
t388 = t442 * t454 + t444 * t449;
t455 = cos(qJ(2));
t443 = sin(pkin(6));
t559 = qJD(1) * t443;
t527 = t455 * t559;
t300 = t388 * t527;
t375 = t388 * qJD(3);
t765 = -t300 + t375;
t404 = qJD(3) - t527;
t398 = qJD(5) + t404;
t641 = t398 / 0.2e1;
t445 = cos(pkin(6));
t558 = qJD(1) * t445;
t423 = qJD(2) + t558;
t450 = sin(qJ(2));
t528 = t450 * t559;
t325 = t423 * t454 - t449 * t528;
t326 = t423 * t449 + t454 * t528;
t241 = t325 * t442 + t326 * t444;
t448 = sin(qJ(5));
t453 = cos(qJ(5));
t515 = t444 * t325 - t326 * t442;
t740 = t241 * t453 + t448 * t515;
t655 = t740 / 0.2e1;
t764 = Ifges(6,4) * t655 + Ifges(6,6) * t641;
t741 = -t241 * t448 + t453 * t515;
t171 = qJD(6) - t741;
t660 = -t171 / 0.2e1;
t447 = sin(qJ(6));
t452 = cos(qJ(6));
t154 = t398 * t447 + t452 * t740;
t668 = -t154 / 0.2e1;
t153 = t398 * t452 - t447 * t740;
t670 = -t153 / 0.2e1;
t763 = Ifges(7,5) * t668 + Ifges(7,6) * t670 + Ifges(7,3) * t660;
t757 = mrSges(6,2) - mrSges(7,3);
t657 = t741 / 0.2e1;
t762 = Ifges(6,2) * t657;
t502 = -mrSges(7,1) * t452 + mrSges(7,2) * t447;
t756 = mrSges(6,1) - t502;
t418 = pkin(8) * t527;
t544 = pkin(1) * t558;
t361 = t450 * t544 + t418;
t508 = t449 * t527;
t554 = qJD(3) * t449;
t715 = -t361 + (-t508 + t554) * pkin(3);
t761 = Ifges(4,3) + Ifges(5,3);
t695 = m(7) * pkin(5);
t760 = -t695 - t756;
t694 = m(7) * pkin(11);
t512 = -t694 + t757;
t358 = -pkin(8) * t528 + t455 * t544;
t476 = (pkin(2) * t450 - pkin(9) * t455) * t443;
t359 = qJD(1) * t476;
t254 = -t449 * t358 + t454 * t359;
t568 = t454 * t455;
t217 = (pkin(3) * t450 - qJ(4) * t568) * t559 + t254;
t255 = t454 * t358 + t449 * t359;
t235 = -qJ(4) * t508 + t255;
t155 = t444 * t217 - t235 * t442;
t387 = -t442 * t449 + t444 * t454;
t301 = t387 * t527;
t134 = pkin(4) * t528 - pkin(10) * t301 + t155;
t156 = t442 * t217 + t444 * t235;
t138 = -pkin(10) * t300 + t156;
t446 = -qJ(4) - pkin(9);
t516 = qJD(3) * t446;
t372 = qJD(4) * t454 + t449 * t516;
t373 = -qJD(4) * t449 + t454 * t516;
t257 = t444 * t372 + t442 * t373;
t221 = -pkin(10) * t375 + t257;
t256 = -t372 * t442 + t444 * t373;
t376 = t387 * qJD(3);
t477 = -pkin(10) * t376 + t256;
t405 = t446 * t449;
t406 = t446 * t454;
t306 = t444 * t405 + t406 * t442;
t258 = -pkin(10) * t388 + t306;
t307 = t442 * t405 - t444 * t406;
t259 = pkin(10) * t387 + t307;
t485 = t453 * t258 - t259 * t448;
t727 = qJD(5) * t485 + (-t138 + t221) * t453 + (-t134 + t477) * t448;
t717 = pkin(4) * t765 + t715;
t170 = Ifges(6,4) * t741;
t590 = t398 * Ifges(6,5);
t109 = Ifges(6,1) * t740 + t170 + t590;
t308 = -t423 * pkin(2) - t358;
t237 = -t325 * pkin(3) + qJD(4) + t308;
t178 = -pkin(4) * t515 + t237;
t309 = pkin(9) * t423 + t361;
t348 = (-pkin(2) * t455 - pkin(9) * t450 - pkin(1)) * t443;
t317 = qJD(1) * t348;
t231 = -t309 * t449 + t454 * t317;
t195 = -qJ(4) * t326 + t231;
t183 = pkin(3) * t404 + t195;
t232 = t309 * t454 + t317 * t449;
t196 = qJ(4) * t325 + t232;
t191 = t442 * t196;
t128 = t444 * t183 - t191;
t749 = pkin(10) * t241;
t101 = pkin(4) * t404 + t128 - t749;
t572 = t444 * t196;
t129 = t442 * t183 + t572;
t732 = pkin(10) * t515;
t110 = t129 + t732;
t56 = t101 * t453 - t110 * t448;
t623 = t56 * mrSges(6,3);
t706 = t178 * mrSges(6,2) - t623;
t759 = Ifges(6,1) * t655 + Ifges(6,4) * t657 + Ifges(6,5) * t641 + t109 / 0.2e1 + t706;
t57 = t101 * t448 + t110 * t453;
t55 = pkin(11) * t398 + t57;
t87 = -pkin(5) * t741 - pkin(11) * t740 + t178;
t25 = -t447 * t55 + t452 * t87;
t26 = t447 * t87 + t452 * t55;
t758 = t763 + t764;
t753 = -t178 * mrSges(6,1) - t25 * mrSges(7,1) + t26 * mrSges(7,2) + t57 * mrSges(6,3) + t758 + t762;
t755 = -pkin(11) * t528 + t727;
t483 = t453 * t387 - t388 * t448;
t203 = qJD(5) * t483 - t375 * t448 + t376 * t453;
t284 = t387 * t448 + t388 * t453;
t204 = qJD(5) * t284 + t453 * t375 + t376 * t448;
t219 = t453 * t300 + t301 * t448;
t220 = -t300 * t448 + t301 * t453;
t754 = t717 + (-t203 + t220) * pkin(11) + (t204 - t219) * pkin(5);
t555 = qJD(2) * t455;
t365 = (qJD(1) * t555 + qJDD(1) * t450) * t443;
t547 = qJDD(1) * t445;
t422 = qJDD(2) + t547;
t233 = qJD(3) * t325 + t365 * t454 + t422 * t449;
t234 = -qJD(3) * t326 - t365 * t449 + t422 * t454;
t164 = t233 * t444 + t234 * t442;
t556 = qJD(2) * t443;
t526 = t450 * t556;
t548 = qJDD(1) * t443;
t364 = -qJD(1) * t526 + t455 * t548;
t352 = qJDD(3) - t364;
t635 = pkin(1) * t445;
t543 = qJD(2) * t635;
t509 = qJD(1) * t543;
t540 = pkin(1) * t547;
t264 = pkin(8) * t364 + t450 * t540 + t455 * t509;
t246 = pkin(9) * t422 + t264;
t251 = -pkin(1) * t548 - pkin(2) * t364 - pkin(9) * t365;
t553 = qJD(3) * t454;
t135 = t454 * t246 + t449 * t251 - t309 * t554 + t317 * t553;
t106 = qJ(4) * t234 + qJD(4) * t325 + t135;
t136 = -qJD(3) * t232 - t246 * t449 + t454 * t251;
t97 = pkin(3) * t352 - qJ(4) * t233 - qJD(4) * t326 + t136;
t52 = -t106 * t442 + t444 * t97;
t37 = pkin(4) * t352 - pkin(10) * t164 + t52;
t163 = -t233 * t442 + t234 * t444;
t53 = t444 * t106 + t442 * t97;
t39 = pkin(10) * t163 + t53;
t551 = qJD(5) * t453;
t552 = qJD(5) * t448;
t10 = t101 * t551 - t110 * t552 + t448 * t37 + t453 * t39;
t576 = t443 * t450;
t424 = pkin(8) * t576;
t265 = -qJD(2) * t418 - qJDD(1) * t424 - t450 * t509 + t455 * t540;
t247 = -t422 * pkin(2) - t265;
t176 = -t234 * pkin(3) + qJDD(4) + t247;
t113 = -t163 * pkin(4) + t176;
t335 = qJDD(5) + t352;
t81 = qJD(5) * t741 + t163 * t448 + t164 * t453;
t49 = qJD(6) * t153 + t335 * t447 + t452 * t81;
t50 = -qJD(6) * t154 + t335 * t452 - t447 * t81;
t82 = -qJD(5) * t740 + t163 * t453 - t164 * t448;
t80 = qJDD(6) - t82;
t14 = Ifges(7,5) * t49 + Ifges(7,6) * t50 + Ifges(7,3) * t80;
t644 = t335 / 0.2e1;
t681 = t82 / 0.2e1;
t682 = t81 / 0.2e1;
t683 = t80 / 0.2e1;
t687 = t50 / 0.2e1;
t688 = t49 / 0.2e1;
t22 = -t82 * pkin(5) - t81 * pkin(11) + t113;
t7 = pkin(11) * t335 + t10;
t2 = qJD(6) * t25 + t22 * t447 + t452 * t7;
t3 = -qJD(6) * t26 + t22 * t452 - t447 * t7;
t705 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t752 = t705 + mrSges(6,1) * t113 - mrSges(6,3) * t10 + Ifges(7,5) * t688 + Ifges(7,6) * t687 + Ifges(7,3) * t683 + t14 / 0.2e1 + (-t644 - t335 / 0.2e1) * Ifges(6,6) + (-t681 - t82 / 0.2e1) * Ifges(6,2) + (-t682 - t81 / 0.2e1) * Ifges(6,4);
t659 = t171 / 0.2e1;
t667 = t154 / 0.2e1;
t669 = t153 / 0.2e1;
t750 = Ifges(7,5) * t667 + Ifges(7,6) * t669 + Ifges(7,3) * t659 - t753 - t762 - t764;
t117 = pkin(5) * t740 - pkin(11) * t741;
t11 = -qJD(5) * t57 + t37 * t453 - t39 * t448;
t748 = t11 * mrSges(6,1) - t10 * mrSges(6,2);
t746 = t264 * mrSges(3,2);
t745 = -m(4) * pkin(9) + m(5) * t446 - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t159 = mrSges(6,1) * t398 - mrSges(6,3) * t740;
t98 = -mrSges(7,1) * t153 + mrSges(7,2) * t154;
t726 = t98 - t159;
t744 = Ifges(4,5) * t233 + Ifges(5,5) * t164 + Ifges(4,6) * t234 + Ifges(5,6) * t163 + t352 * t761;
t642 = -t398 / 0.2e1;
t656 = -t740 / 0.2e1;
t658 = -t741 / 0.2e1;
t743 = -Ifges(6,4) * t656 - Ifges(6,2) * t658 - Ifges(6,6) * t642 + t753 + t763;
t441 = qJ(3) + pkin(12);
t438 = qJ(5) + t441;
t431 = sin(t438);
t432 = cos(t438);
t629 = pkin(3) * t454;
t434 = pkin(2) + t629;
t436 = sin(t441);
t437 = cos(t441);
t504 = -t454 * mrSges(4,1) + t449 * mrSges(4,2);
t742 = -m(4) * pkin(2) - m(5) * t434 - t437 * mrSges(5,1) + t436 * mrSges(5,2) + t431 * t512 + t432 * t760 + t504;
t501 = mrSges(7,1) * t447 + mrSges(7,2) * t452;
t54 = -pkin(5) * t398 - t56;
t474 = t54 * t501;
t152 = Ifges(7,4) * t153;
t72 = t154 * Ifges(7,1) + t171 * Ifges(7,5) + t152;
t584 = t452 * t72;
t637 = t447 / 0.2e1;
t608 = Ifges(7,4) * t154;
t71 = Ifges(7,2) * t153 + Ifges(7,6) * t171 + t608;
t739 = -t584 / 0.2e1 + t71 * t637 - t474 - t706;
t738 = -t136 * mrSges(4,1) - t52 * mrSges(5,1) + t135 * mrSges(4,2) + t53 * mrSges(5,2);
t736 = t113 * mrSges(6,2) - mrSges(6,3) * t11 + 0.2e1 * Ifges(6,1) * t682 + 0.2e1 * Ifges(6,4) * t681 + 0.2e1 * Ifges(6,5) * t644;
t733 = m(6) + m(7);
t666 = t163 / 0.2e1;
t665 = t164 / 0.2e1;
t652 = t233 / 0.2e1;
t651 = t234 / 0.2e1;
t643 = t352 / 0.2e1;
t336 = -pkin(4) * t387 - t434;
t184 = -pkin(5) * t483 - pkin(11) * t284 + t336;
t188 = t258 * t448 + t259 * t453;
t124 = t184 * t447 + t188 * t452;
t731 = -qJD(6) * t124 - t447 * t755 + t452 * t754;
t123 = t184 * t452 - t188 * t447;
t730 = qJD(6) * t123 + t447 * t754 + t452 * t755;
t696 = m(5) * pkin(3);
t728 = -t696 - mrSges(4,1);
t573 = t443 * t455;
t385 = pkin(8) * t573 + t450 * t635;
t347 = pkin(9) * t445 + t385;
t249 = -t449 * t347 + t454 * t348;
t574 = t443 * t454;
t378 = t445 * t449 + t450 * t574;
t207 = -pkin(3) * t573 - t378 * qJ(4) + t249;
t250 = t454 * t347 + t449 * t348;
t377 = t445 * t454 - t449 * t576;
t216 = qJ(4) * t377 + t250;
t143 = t444 * t207 - t442 * t216;
t261 = t377 * t442 + t378 * t444;
t122 = -pkin(4) * t573 - t261 * pkin(10) + t143;
t144 = t442 * t207 + t444 * t216;
t260 = t377 * t444 - t378 * t442;
t126 = pkin(10) * t260 + t144;
t723 = t448 * t122 + t453 * t126;
t722 = -t155 + t256;
t721 = -t156 + t257;
t158 = -mrSges(6,2) * t398 + mrSges(6,3) * t741;
t104 = -mrSges(7,2) * t171 + mrSges(7,3) * t153;
t105 = mrSges(7,1) * t171 - mrSges(7,3) * t154;
t488 = t104 * t452 - t105 * t447;
t720 = -t158 - t488;
t201 = -t220 * t447 + t452 * t528;
t549 = qJD(6) * t452;
t473 = t447 * t203 + t284 * t549;
t719 = t201 + t473;
t202 = t220 * t452 + t447 * t528;
t550 = qJD(6) * t447;
t472 = -t452 * t203 + t284 * t550;
t718 = t202 + t472;
t511 = mrSges(3,3) * t528;
t716 = -mrSges(3,1) * t423 - mrSges(4,1) * t325 + mrSges(4,2) * t326 + t511;
t634 = pkin(1) * t455;
t384 = t445 * t634 - t424;
t631 = pkin(3) * t444;
t433 = pkin(4) + t631;
t632 = pkin(3) * t442;
t371 = t448 * t433 + t453 * t632;
t461 = mrSges(3,2) + t745;
t714 = -t461 + t501;
t319 = -t431 * t576 + t432 * t445;
t320 = t431 * t445 + t432 * t576;
t712 = -t756 * t319 + t320 * t757;
t636 = cos(qJ(1));
t529 = t636 * t455;
t451 = sin(qJ(1));
t570 = t450 * t451;
t382 = -t445 * t570 + t529;
t575 = t443 * t451;
t281 = t382 * t431 - t432 * t575;
t282 = t382 * t432 + t431 * t575;
t711 = t756 * t281 + t282 * t757;
t530 = t636 * t450;
t569 = t451 * t455;
t380 = t445 * t530 + t569;
t531 = t443 * t636;
t277 = -t380 * t431 - t432 * t531;
t278 = t380 * t432 - t431 * t531;
t710 = -t756 * t277 + t278 * t757;
t708 = t135 * t454 - t136 * t449;
t23 = mrSges(7,1) * t80 - mrSges(7,3) * t49;
t24 = -mrSges(7,2) * t80 + mrSges(7,3) * t50;
t707 = -t447 * t23 + t452 * t24;
t525 = t443 * t555;
t288 = -qJD(3) * t378 - t449 * t525;
t289 = qJD(3) * t377 + t454 * t525;
t213 = t288 * t442 + t289 * t444;
t360 = qJD(2) * t476;
t362 = t384 * qJD(2);
t186 = -qJD(3) * t250 + t454 * t360 - t362 * t449;
t141 = pkin(3) * t526 - qJ(4) * t289 - qJD(4) * t378 + t186;
t185 = -t347 * t554 + t348 * t553 + t449 * t360 + t454 * t362;
t149 = qJ(4) * t288 + qJD(4) * t377 + t185;
t91 = t444 * t141 - t149 * t442;
t69 = pkin(4) * t526 - pkin(10) * t213 + t91;
t212 = t288 * t444 - t289 * t442;
t92 = t442 * t141 + t444 * t149;
t75 = pkin(10) * t212 + t92;
t20 = -qJD(5) * t723 - t448 * t75 + t453 * t69;
t703 = mrSges(3,1) - t742;
t674 = -t109 / 0.2e1;
t701 = Ifges(6,1) * t656 + Ifges(6,4) * t658 + Ifges(6,5) * t642 + t674;
t613 = Ifges(3,4) * t450;
t498 = Ifges(3,2) * t455 + t613;
t605 = Ifges(3,6) * t423;
t700 = t129 * mrSges(5,2) + t498 * t559 / 0.2e1 + t605 / 0.2e1 - t128 * mrSges(5,1) - t56 * mrSges(6,1);
t697 = t443 ^ 2;
t15 = t49 * Ifges(7,4) + t50 * Ifges(7,2) + t80 * Ifges(7,6);
t692 = t15 / 0.2e1;
t16 = t49 * Ifges(7,1) + t50 * Ifges(7,4) + t80 * Ifges(7,5);
t691 = t16 / 0.2e1;
t680 = Ifges(5,4) * t665 + Ifges(5,2) * t666 + Ifges(5,6) * t643;
t679 = Ifges(5,1) * t665 + Ifges(5,4) * t666 + Ifges(5,5) * t643;
t678 = pkin(1) * mrSges(3,1);
t677 = pkin(1) * mrSges(3,2);
t672 = Ifges(4,4) * t652 + Ifges(4,2) * t651 + Ifges(4,6) * t643;
t671 = Ifges(4,1) * t652 + Ifges(4,4) * t651 + Ifges(4,5) * t643;
t166 = t241 * Ifges(5,4) + Ifges(5,2) * t515 + t404 * Ifges(5,6);
t664 = -t166 / 0.2e1;
t663 = t166 / 0.2e1;
t167 = t241 * Ifges(5,1) + Ifges(5,4) * t515 + t404 * Ifges(5,5);
t662 = -t167 / 0.2e1;
t661 = t167 / 0.2e1;
t592 = t326 * Ifges(4,4);
t223 = t325 * Ifges(4,2) + t404 * Ifges(4,6) + t592;
t654 = t223 / 0.2e1;
t318 = Ifges(4,4) * t325;
t224 = t326 * Ifges(4,1) + t404 * Ifges(4,5) + t318;
t653 = t224 / 0.2e1;
t650 = -t515 / 0.2e1;
t649 = t515 / 0.2e1;
t648 = -t241 / 0.2e1;
t647 = t241 / 0.2e1;
t645 = t326 / 0.2e1;
t640 = -t404 / 0.2e1;
t639 = t404 / 0.2e1;
t638 = t445 / 0.2e1;
t633 = pkin(3) * t326;
t630 = pkin(3) * t449;
t626 = t377 * pkin(3);
t617 = mrSges(5,3) * t128;
t616 = mrSges(5,3) * t129;
t615 = mrSges(7,3) * t447;
t614 = mrSges(7,3) * t452;
t612 = Ifges(3,4) * t455;
t611 = Ifges(4,4) * t449;
t610 = Ifges(4,4) * t454;
t607 = Ifges(7,4) * t447;
t606 = Ifges(7,4) * t452;
t599 = t741 * Ifges(6,6);
t598 = t740 * Ifges(6,5);
t597 = t231 * mrSges(4,3);
t596 = t232 * mrSges(4,3);
t595 = t515 * Ifges(5,6);
t594 = t241 * Ifges(5,5);
t593 = t325 * Ifges(4,6);
t591 = t326 * Ifges(4,5);
t588 = t398 * Ifges(6,3);
t587 = t423 * Ifges(3,5);
t578 = t284 * t447;
t577 = t284 * t452;
t132 = t444 * t195 - t191;
t396 = pkin(4) * t436 + t630;
t397 = pkin(4) * t437 + t629;
t563 = -t382 * t396 + t397 * t575;
t561 = -t396 * t576 + t445 * t397;
t363 = pkin(8) * t525 + t450 * t543;
t560 = t636 * pkin(1) + pkin(8) * t575;
t545 = Ifges(6,5) * t81 + Ifges(6,6) * t82 + Ifges(6,3) * t335;
t539 = t447 * t573;
t538 = t452 * t573;
t535 = t584 / 0.2e1;
t532 = Ifges(3,5) * t365 + Ifges(3,6) * t364 + Ifges(3,3) * t422;
t31 = -t82 * mrSges(6,1) + t81 * mrSges(6,2);
t519 = -t550 / 0.2e1;
t518 = -t451 * pkin(1) + pkin(8) * t531;
t517 = -t281 * pkin(5) + pkin(11) * t282;
t103 = -t163 * mrSges(5,1) + t164 * mrSges(5,2);
t131 = -t195 * t442 - t572;
t514 = -t380 * t437 + t436 * t531;
t411 = t449 * t531;
t513 = -t380 * t454 + t411;
t510 = mrSges(3,3) * t527;
t245 = -pkin(3) * t288 + t363;
t200 = pkin(4) * t241 + t633;
t505 = mrSges(4,1) * t377 - mrSges(4,2) * t378;
t500 = Ifges(4,1) * t454 - t611;
t499 = Ifges(7,1) * t452 - t607;
t497 = -Ifges(4,2) * t449 + t610;
t496 = -Ifges(7,2) * t447 + t606;
t495 = Ifges(4,5) * t454 - Ifges(4,6) * t449;
t494 = Ifges(7,5) * t452 - Ifges(7,6) * t447;
t493 = t25 * t452 + t26 * t447;
t492 = -t25 * t447 + t26 * t452;
t61 = -pkin(11) * t573 + t723;
t190 = t260 * t448 + t261 * t453;
t346 = t424 + (-pkin(2) - t634) * t445;
t272 = t346 - t626;
t199 = -t260 * pkin(4) + t272;
t484 = t453 * t260 - t261 * t448;
t99 = -pkin(5) * t484 - t190 * pkin(11) + t199;
t36 = t447 * t99 + t452 * t61;
t35 = -t447 * t61 + t452 * t99;
t381 = t445 * t569 + t530;
t395 = pkin(2) + t397;
t440 = -pkin(10) + t446;
t489 = -t381 * t440 + t382 * t395 + t396 * t575 + t560;
t62 = t453 * t122 - t448 * t126;
t88 = t134 * t453 - t138 * t448;
t479 = -t380 * t396 - t397 * t531;
t370 = t433 * t453 - t448 * t632;
t478 = t131 - t732;
t179 = -t447 * t190 - t538;
t475 = -t452 * t190 + t539;
t291 = -t382 * t449 + t451 * t574;
t470 = t153 * t496;
t469 = t154 * t499;
t468 = t171 * t494;
t19 = t122 * t551 - t126 * t552 + t448 * t69 + t453 * t75;
t169 = -pkin(4) * t212 + t245;
t463 = t380 * t436 + t437 * t531;
t462 = t380 * t449 + t454 * t531;
t379 = -t445 * t529 + t570;
t460 = -g(1) * t381 - g(2) * t379 + g(3) * t573;
t459 = -qJD(6) * t493 - t3 * t447;
t458 = t2 * t452 + t459;
t8 = -pkin(5) * t335 - t11;
t456 = t16 * t637 + t2 * t614 + t452 * t692 + t8 * t502 + (Ifges(7,1) * t447 + t606) * t688 + (Ifges(7,2) * t452 + t607) * t687 + t545 + t71 * t519 + (Ifges(7,5) * t447 + Ifges(7,6) * t452) * t683 + (t474 + t535) * qJD(6) + (t470 + t469 + t468) * qJD(6) / 0.2e1 + t748;
t416 = Ifges(3,4) * t527;
t383 = (-mrSges(3,1) * t455 + mrSges(3,2) * t450) * t443;
t356 = -pkin(5) - t370;
t355 = -t423 * mrSges(3,2) + t510;
t314 = t319 * pkin(5);
t304 = Ifges(3,1) * t528 + t416 + t587;
t292 = t382 * t454 + t449 * t575;
t286 = t382 * t437 + t436 * t575;
t285 = -t382 * t436 + t437 * t575;
t276 = mrSges(4,1) * t404 - mrSges(4,3) * t326;
t275 = -mrSges(4,2) * t404 + mrSges(4,3) * t325;
t273 = t277 * pkin(5);
t226 = t282 * t452 + t381 * t447;
t225 = -t282 * t447 + t381 * t452;
t222 = t404 * Ifges(4,3) + t591 + t593;
t215 = mrSges(5,1) * t404 - mrSges(5,3) * t241;
t214 = -mrSges(5,2) * t404 + mrSges(5,3) * t515;
t198 = -mrSges(4,2) * t352 + mrSges(4,3) * t234;
t197 = mrSges(4,1) * t352 - mrSges(4,3) * t233;
t177 = -mrSges(5,1) * t515 + mrSges(5,2) * t241;
t168 = -mrSges(4,1) * t234 + mrSges(4,2) * t233;
t165 = t404 * Ifges(5,3) + t594 + t595;
t146 = mrSges(5,1) * t352 - mrSges(5,3) * t164;
t145 = -mrSges(5,2) * t352 + mrSges(5,3) * t163;
t119 = qJD(5) * t188 + t221 * t448 - t453 * t477;
t116 = t132 - t749;
t115 = -mrSges(6,1) * t741 + mrSges(6,2) * t740;
t112 = qJD(5) * t190 - t453 * t212 + t213 * t448;
t111 = qJD(5) * t484 + t212 * t448 + t213 * t453;
t107 = t588 + t598 + t599;
t90 = t117 + t200;
t85 = -pkin(5) * t528 - t88;
t84 = qJD(6) * t475 - t447 * t111 + t452 * t526;
t83 = qJD(6) * t179 + t452 * t111 + t447 * t526;
t65 = -mrSges(6,2) * t335 + mrSges(6,3) * t82;
t64 = mrSges(6,1) * t335 - mrSges(6,3) * t81;
t60 = pkin(5) * t573 - t62;
t59 = t453 * t116 + t448 * t478;
t42 = pkin(5) * t112 - pkin(11) * t111 + t169;
t33 = t117 * t447 + t452 * t56;
t32 = t117 * t452 - t447 * t56;
t28 = t447 * t90 + t452 * t59;
t27 = -t447 * t59 + t452 * t90;
t21 = -mrSges(7,1) * t50 + mrSges(7,2) * t49;
t18 = -pkin(5) * t526 - t20;
t17 = pkin(11) * t526 + t19;
t5 = -qJD(6) * t36 - t17 * t447 + t42 * t452;
t4 = qJD(6) * t35 + t17 * t452 + t42 * t447;
t1 = [(Ifges(4,1) * t378 + Ifges(4,4) * t377) * t652 + (Ifges(3,5) * t638 - t384 * mrSges(3,3) + (t450 * Ifges(3,1) + t612 - t677) * t443) * t365 + (Ifges(3,6) * t638 + t385 * mrSges(3,3) + (t498 + t678) * t443) * t364 + (mrSges(3,3) * t264 - Ifges(4,5) * t652 - Ifges(5,5) * t665 - Ifges(6,5) * t682 - Ifges(4,6) * t651 - Ifges(5,6) * t666 - Ifges(6,6) * t681 - Ifges(6,3) * t644 - t643 * t761 + t738 - t748) * t573 + (Ifges(7,5) * t83 + Ifges(7,6) * t84) * t659 - t752 * t484 - t247 * t505 + (Ifges(4,4) * t378 + Ifges(4,2) * t377) * t651 + (-pkin(1) * t383 * t443 + Ifges(2,3)) * qJDD(1) - t445 * t746 + (-m(4) * (-pkin(2) * t380 + t518) - t513 * mrSges(4,1) - t462 * mrSges(4,2) - m(5) * (pkin(3) * t411 - t380 * t434 + t518) - t514 * mrSges(5,1) - t463 * mrSges(5,2) - m(3) * t518 + t380 * mrSges(3,1) - mrSges(3,3) * t531 + t451 * mrSges(2,1) + t636 * mrSges(2,2) + t512 * t277 - t760 * t278 + t714 * t379 + t733 * (-t379 * t440 + t380 * t395 - t396 * t531 - t518)) * g(1) + t759 * t111 + m(6) * (t10 * t723 + t11 * t62 + t113 * t199 + t169 * t178 + t19 * t57 + t20 * t56) + t723 * t65 + (t179 * t2 - t25 * t83 + t26 * t84 + t3 * t475) * mrSges(7,3) + (-Ifges(7,5) * t475 + Ifges(7,6) * t179) * t683 + (-Ifges(7,4) * t475 + Ifges(7,2) * t179) * t687 + (-Ifges(7,1) * t475 + Ifges(7,4) * t179) * t688 + t8 * (-mrSges(7,1) * t179 - mrSges(7,2) * t475) - t475 * t691 + (Ifges(4,5) * t289 + Ifges(5,5) * t213 + Ifges(4,6) * t288 + Ifges(5,6) * t212) * t639 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t697 + t264 * t385 + t265 * t384 - t358 * t363 + t361 * t362) + (-t128 * t213 + t129 * t212 + t260 * t53 - t261 * t52) * mrSges(5,3) + (Ifges(5,1) * t261 + Ifges(5,4) * t260) * t665 + t716 * t363 + t750 * t112 + t736 * t190 + (Ifges(4,5) * t378 + Ifges(5,5) * t261 + Ifges(4,6) * t377 + Ifges(5,6) * t260) * t643 + (-t636 * mrSges(2,1) - m(5) * (t382 * t434 + t560) - t286 * mrSges(5,1) - t285 * mrSges(5,2) - m(3) * t560 - t382 * mrSges(3,1) - m(4) * (pkin(2) * t382 + t560) - t292 * mrSges(4,1) - t291 * mrSges(4,2) - m(6) * t489 - t282 * mrSges(6,1) - m(7) * (pkin(5) * t282 + t489) - t226 * mrSges(7,1) - t225 * mrSges(7,2) + (mrSges(2,2) + (-m(5) * t630 - mrSges(3,3)) * t443) * t451 + t512 * t281 + t461 * t381) * g(2) + (Ifges(3,3) * t638 - t385 * mrSges(3,2) + t384 * mrSges(3,1) + (Ifges(3,5) * t450 + Ifges(3,6) * t455) * t443) * t422 + (Ifges(7,4) * t83 + Ifges(7,2) * t84) * t669 + t325 * (Ifges(4,4) * t289 + Ifges(4,2) * t288) / 0.2e1 + t4 * t104 + t5 * t105 + t18 * t98 + (Ifges(7,1) * t83 + Ifges(7,4) * t84) * t667 + ((-t358 * mrSges(3,3) + t587 / 0.2e1 + t304 / 0.2e1 + (-t677 + t612 / 0.2e1) * t559) * t455 + (t588 / 0.2e1 + (-t678 - t613 / 0.2e1 + (-Ifges(3,2) / 0.2e1 + Ifges(3,1) / 0.2e1) * t455) * t559 + t594 / 0.2e1 + t595 / 0.2e1 - t605 / 0.2e1 + t591 / 0.2e1 + t593 / 0.2e1 - t700 + t598 / 0.2e1 + t599 / 0.2e1 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t404 + t222 / 0.2e1 - t232 * mrSges(4,2) + t231 * mrSges(4,1) - t57 * mrSges(6,2) - t361 * mrSges(3,3) + t165 / 0.2e1 + t107 / 0.2e1) * t450) * t556 + m(7) * (t18 * t54 + t2 * t36 + t25 * t5 + t26 * t4 + t3 * t35 + t60 * t8) + m(5) * (t128 * t91 + t129 * t92 + t143 * t52 + t144 * t53 + t176 * t272 + t237 * t245) + m(4) * (t135 * t250 + t136 * t249 + t185 * t232 + t186 * t231 + t247 * t346 + t308 * t363) + t83 * t72 / 0.2e1 + t54 * (-mrSges(7,1) * t84 + mrSges(7,2) * t83) + t84 * t71 / 0.2e1 - (t545 + t744) * t573 / 0.2e1 + (t135 * t377 - t136 * t378 - t231 * t289 + t232 * t288) * mrSges(4,3) + t265 * (mrSges(3,1) * t445 - mrSges(3,3) * t576) + (Ifges(5,4) * t261 + Ifges(5,2) * t260) * t666 + t62 * t64 + t60 * t21 + t36 * t24 + t35 * t23 + t179 * t692 + t261 * t679 + t260 * t680 + t378 * t671 + t377 * t672 + t289 * t653 + t288 * t654 + t213 * t661 + t212 * t663 + t362 * t355 + t346 * t168 + t144 * t145 + t143 * t146 + t19 * t158 + t20 * t159 + t169 * t115 + t532 * t638 + (Ifges(4,1) * t289 + Ifges(4,4) * t288) * t645 + (Ifges(5,1) * t213 + Ifges(5,4) * t212) * t647 + (Ifges(5,4) * t213 + Ifges(5,2) * t212) * t649 + t199 * t31 + t92 * t214 + t91 * t215 + t237 * (-mrSges(5,1) * t212 + mrSges(5,2) * t213) + t245 * t177 + t249 * t197 + t250 * t198 + t176 * (-mrSges(5,1) * t260 + mrSges(5,2) * t261) + t272 * t103 + t185 * t275 + t186 * t276 + t308 * (-mrSges(4,1) * t288 + mrSges(4,2) * t289); (-pkin(2) * t247 - t231 * t254 - t232 * t255) * m(4) - t752 * t483 + t247 * t504 + (Ifges(5,1) * t647 + Ifges(5,4) * t649 + Ifges(5,5) * t639 - t617 + t661) * t376 - t15 * t578 / 0.2e1 + (t535 + t759) * t203 - ((-Ifges(3,2) * t528 + t454 * t224 + t304 + t416) * t455 + (t222 + t165 + t107) * t450 + t404 * (Ifges(4,3) * t450 + t455 * t495) + t326 * (Ifges(4,5) * t450 + t455 * t500) + t325 * (Ifges(4,6) * t450 + t455 * t497) + t423 * (Ifges(3,5) * t455 - Ifges(3,6) * t450)) * t559 / 0.2e1 - (t21 - t64) * t485 + (t123 * t3 + t124 * t2 - t485 * t8 + (-t85 + t119) * t54 + t730 * t26 + t731 * t25) * m(7) + (t10 * t188 + t11 * t485 + t113 * t336 + t727 * t57 + (-t88 - t119) * t56 + t717 * t178) * m(6) + (-t178 * t220 + t528 * t57) * mrSges(6,2) + t404 * t308 * (mrSges(4,1) * t449 + mrSges(4,2) * t454) + (m(4) * ((-t231 * t454 - t232 * t449) * qJD(3) + t708) - t276 * t553 - t275 * t554 + t454 * t198 - t449 * t197) * pkin(9) + t708 * mrSges(4,3) + t715 * t177 + (-m(4) * t308 + t511 - t716) * t361 + t717 * t115 + (t325 * t497 + t326 * t500 + t404 * t495) * qJD(3) / 0.2e1 + t750 * t204 + (-Ifges(7,5) * t472 - Ifges(7,6) * t473) * t659 + (Ifges(7,5) * t202 + Ifges(7,6) * t201) * t660 + (t494 * t683 + t496 * t687 + t499 * t688 + t501 * t8 + t519 * t72 + t736) * t284 + (Ifges(5,5) * t301 - Ifges(5,6) * t300) * t640 + (Ifges(5,1) * t301 - Ifges(5,4) * t300) * t648 + (Ifges(5,4) * t301 - Ifges(5,2) * t300) * t650 + (t128 * t301 + t129 * t300 + t387 * t53 - t388 * t52) * mrSges(5,3) + (-t733 * (-t381 * t395 - t382 * t440) - t714 * t382 + t703 * t381) * g(1) + (-t733 * (-t379 * t395 - t380 * t440) - t714 * t380 + t703 * t379) * g(2) + ((-t301 + t376) * mrSges(5,2) + t765 * mrSges(5,1)) * t237 + t730 * t104 + t731 * t105 + (-t450 * (Ifges(3,1) * t455 - t613) / 0.2e1 + pkin(1) * (mrSges(3,1) * t450 + mrSges(3,2) * t455)) * qJD(1) ^ 2 * t697 - t746 + t532 + (-t231 * (mrSges(4,1) * t450 - mrSges(4,3) * t568) - t232 * (-mrSges(4,3) * t449 * t455 - mrSges(4,2) * t450)) * t559 - (Ifges(5,4) * t647 + Ifges(5,2) * t649 + Ifges(5,6) * t639 + t616 + t663) * t375 - t85 * t98 + (Ifges(5,5) * t648 + Ifges(6,5) * t656 + Ifges(5,6) * t650 + Ifges(6,6) * t658 + Ifges(5,3) * t640 + Ifges(6,3) * t642 + t700) * t528 + (t510 - t355) * t358 + (-Ifges(7,4) * t472 - Ifges(7,2) * t473) * t669 + (Ifges(7,4) * t202 + Ifges(7,2) * t201) * t670 + (Ifges(7,1) * t202 + Ifges(7,4) * t201) * t668 + (-Ifges(7,1) * t472 - Ifges(7,4) * t473) * t667 + (t701 + t623) * t220 + (-t733 * t395 * t573 + t383 + (t742 * t455 + (t440 * t733 - t501 + t745) * t450) * t443) * g(3) + t743 * t219 + t726 * t119 + t727 * t158 - t719 * t71 / 0.2e1 + (mrSges(7,1) * t719 - mrSges(7,2) * t718) * t54 + (-t2 * t578 + t25 * t718 - t26 * t719 - t3 * t577) * mrSges(7,3) + t721 * t214 + t722 * t215 + (t128 * t722 + t129 * t721 - t176 * t434 + t237 * t715 + t306 * t52 + t307 * t53) * m(5) + (-t596 - t223 / 0.2e1) * t554 + (-t597 + t653) * t553 + t577 * t691 + t388 * t679 + t387 * t680 + t449 * t671 + t454 * t672 + (Ifges(4,1) * t449 + t610) * t652 + t508 * t654 + t301 * t662 - t300 * t664 + (Ifges(5,1) * t388 + Ifges(5,4) * t387) * t665 + (Ifges(5,4) * t388 + Ifges(5,2) * t387) * t666 - t434 * t103 + t176 * (-mrSges(5,1) * t387 + mrSges(5,2) * t388) + t336 * t31 + t123 * t23 + t124 * t24 - t88 * t159 - pkin(2) * t168 + (Ifges(4,2) * t454 + t611) * t651 + t188 * t65 - t202 * t72 / 0.2e1 + t265 * mrSges(3,1) + (Ifges(4,5) * t449 + Ifges(5,5) * t388 + Ifges(4,6) * t454 + Ifges(5,6) * t387) * t643 - t255 * t275 - t254 * t276 + t306 * t146 + t307 * t145; -t738 - t3 * t615 - t177 * t633 - m(5) * (t128 * t131 + t129 * t132 + t237 * t633) + t744 - (mrSges(5,1) * t237 + Ifges(5,4) * t648 + Ifges(5,2) * t650 + Ifges(5,6) * t640 - t616 + t664) * t241 + (m(7) * t458 - t104 * t550 - t105 * t549 + t707) * (pkin(11) + t371) + (-t505 - (-t436 * t576 + t437 * t445) * mrSges(5,1) - (-t436 * t445 - t437 * t576) * mrSges(5,2) - m(5) * t626 - m(7) * (pkin(11) * t320 + t314 + t561) - m(6) * t561 + t712) * g(3) - (-Ifges(4,2) * t326 + t224 + t318) * t325 / 0.2e1 + t743 * t740 + (t25 * t614 + t26 * t615 + t494 * t660 + t496 * t670 + t499 * t668 + t701 + t739) * t741 + (-t25 * t549 - t26 * t550) * mrSges(7,3) - t28 * t104 - t27 * t105 + (-t25 * t27 - t26 * t28 + t356 * t8) * m(7) + (-m(6) * t56 + m(7) * t54 + t726) * (t371 * qJD(5) - t116 * t448 + t453 * t478) + (t10 * t371 - t178 * t200 - t57 * t59) * m(6) - t326 * (Ifges(4,1) * t325 - t592) / 0.2e1 + (-t285 * mrSges(5,1) + t286 * mrSges(5,2) - m(7) * (t517 + t563) - m(6) * t563 + mrSges(4,2) * t292 + t728 * t291 + t711) * g(1) + (-t513 * mrSges(4,2) + t463 * mrSges(5,1) - t514 * mrSges(5,2) - m(7) * (t278 * pkin(11) + t273 + t479) - m(6) * t479 - t728 * t462 + t710) * g(2) + (-mrSges(5,2) * t237 + Ifges(5,1) * t648 + Ifges(5,4) * t650 + Ifges(5,5) * t640 + t617 + t662) * t515 + (t442 * t53 + t444 * t52) * t696 + t371 * t65 + t356 * t21 + t326 * t596 + t325 * t597 + t456 - t59 * t158 + t146 * t631 + t145 * t632 + (Ifges(4,5) * t325 - Ifges(4,6) * t326) * t640 + t223 * t645 - t200 * t115 - t132 * t214 - t131 * t215 + ((m(6) * t57 + m(7) * t492 - t720) * qJD(5) + t11 * m(6) + t64) * t370 - t231 * t275 + t232 * t276 - t308 * (mrSges(4,1) * t326 + mrSges(4,2) * t325); -t515 * t214 + t241 * t215 + t452 * t23 + t447 * t24 - t726 * t740 + t488 * qJD(6) + t720 * t741 + t103 + t31 + (t171 * t492 + t2 * t447 + t3 * t452 - t740 * t54 + t460) * m(7) + (t56 * t740 - t57 * t741 + t113 + t460) * m(6) + (t128 * t241 - t129 * t515 + t176 + t460) * m(5); (-m(7) * t273 + t710) * g(2) + (-m(7) * t517 + t711) * g(1) + (-t470 / 0.2e1 - t469 / 0.2e1 - t468 / 0.2e1 - t590 / 0.2e1 + t674 - t170 / 0.2e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t740 + t493 * mrSges(7,3) + t739) * t741 + t459 * mrSges(7,3) + ((-t104 * t447 - t105 * t452) * qJD(6) + (-g(2) * t278 - g(3) * t320) * m(7) + t707) * pkin(11) + t458 * t694 - t8 * t695 - t726 * t57 + (t753 + t758) * t740 - t33 * t104 - t32 * t105 + (-m(7) * t314 + t712) * g(3) - m(7) * (t25 * t32 + t26 * t33 + t54 * t57) - pkin(5) * t21 + t456 - t56 * t158; -t54 * (mrSges(7,1) * t154 + mrSges(7,2) * t153) + (Ifges(7,1) * t153 - t608) * t668 + t71 * t667 + (Ifges(7,5) * t153 - Ifges(7,6) * t154) * t660 - t25 * t104 + t26 * t105 - g(1) * (mrSges(7,1) * t225 - mrSges(7,2) * t226) - g(2) * ((-t278 * t447 + t379 * t452) * mrSges(7,1) + (-t278 * t452 - t379 * t447) * mrSges(7,2)) - g(3) * ((-t320 * t447 - t538) * mrSges(7,1) + (-t320 * t452 + t539) * mrSges(7,2)) + (t153 * t25 + t154 * t26) * mrSges(7,3) + t14 + (-Ifges(7,2) * t154 + t152 + t72) * t670 + t705;];
tau  = t1;
