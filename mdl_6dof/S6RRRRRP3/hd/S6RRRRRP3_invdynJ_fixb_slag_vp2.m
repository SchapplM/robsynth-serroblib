% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:52
% EndTime: 2019-03-10 01:06:47
% DurationCPUTime: 30.58s
% Computational Cost: add. (22474->885), mult. (48276->1122), div. (0->0), fcn. (35325->14), ass. (0->413)
t731 = Ifges(6,4) + Ifges(7,4);
t732 = Ifges(6,1) + Ifges(7,1);
t730 = Ifges(6,5) + Ifges(7,5);
t729 = Ifges(6,2) + Ifges(7,2);
t728 = Ifges(6,6) + Ifges(7,6);
t450 = qJ(4) + qJ(5);
t443 = cos(t450);
t734 = -mrSges(7,1) - mrSges(6,1);
t771 = t443 * t734;
t454 = sin(qJ(3));
t459 = cos(qJ(3));
t460 = cos(qJ(2));
t548 = qJD(1) * t460;
t455 = sin(qJ(2));
t549 = qJD(1) * t455;
t354 = -t454 * t549 + t459 * t548;
t379 = t454 * t460 + t455 * t459;
t355 = t379 * qJD(1);
t281 = pkin(3) * t355 - pkin(9) * t354;
t531 = pkin(2) * t549;
t258 = t281 + t531;
t463 = -pkin(8) - pkin(7);
t409 = t463 * t460;
t386 = qJD(1) * t409;
t357 = t454 * t386;
t407 = t463 * t455;
t385 = qJD(1) * t407;
t290 = t385 * t459 + t357;
t453 = sin(qJ(4));
t458 = cos(qJ(4));
t192 = t453 * t258 + t458 * t290;
t545 = qJD(3) * t459;
t529 = pkin(2) * t545;
t770 = t458 * t529 - t192;
t191 = t458 * t258 - t290 * t453;
t769 = -t453 * t529 - t191;
t727 = Ifges(7,3) + Ifges(6,3);
t441 = sin(t450);
t733 = -mrSges(7,2) - mrSges(6,2);
t758 = t441 * t733;
t588 = t354 * t458;
t499 = t355 * pkin(4) - pkin(10) * t588;
t633 = pkin(2) * t454;
t431 = pkin(9) + t633;
t617 = -pkin(10) - t431;
t510 = qJD(4) * t617;
t768 = t458 * t510 - t499 + t769;
t364 = qJD(2) * pkin(2) + t385;
t285 = t364 * t459 + t357;
t194 = t458 * t281 - t285 * t453;
t462 = -pkin(10) - pkin(9);
t518 = qJD(4) * t462;
t767 = t458 * t518 - t194 - t499;
t589 = t354 * t453;
t537 = pkin(10) * t589;
t766 = -t453 * t510 - t537 - t770;
t195 = t453 * t281 + t458 * t285;
t765 = -t453 * t518 + t195 - t537;
t449 = qJD(2) + qJD(3);
t306 = -t355 * t453 + t449 * t458;
t307 = t355 * t458 + t449 * t453;
t452 = sin(qJ(5));
t457 = cos(qJ(5));
t505 = t457 * t306 - t307 * t452;
t764 = t731 * t505;
t378 = t452 * t458 + t453 * t457;
t251 = t378 * t354;
t689 = qJD(4) + qJD(5);
t293 = t689 * t378;
t763 = t251 - t293;
t483 = t452 * t453 - t457 * t458;
t252 = t483 * t354;
t292 = t689 * t483;
t748 = t252 - t292;
t216 = t306 * t452 + t307 * t457;
t762 = t731 * t216;
t759 = -mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t540 = qJD(1) * qJD(2);
t392 = qJDD(1) * t460 - t455 * t540;
t393 = qJDD(1) * t455 + t460 * t540;
t242 = -qJD(3) * t355 + t392 * t459 - t454 * t393;
t238 = qJDD(4) - t242;
t233 = qJDD(5) + t238;
t377 = t454 * t455 - t459 * t460;
t474 = t377 * qJD(3);
t241 = -qJD(1) * t474 + t392 * t454 + t393 * t459;
t447 = qJDD(2) + qJDD(3);
t174 = qJD(4) * t306 + t241 * t458 + t447 * t453;
t175 = -qJD(4) * t307 - t241 * t453 + t447 * t458;
t65 = qJD(5) * t505 + t174 * t457 + t175 * t452;
t66 = -qJD(5) * t216 - t174 * t452 + t175 * t457;
t761 = -t729 * t66 / 0.2e1 - t731 * t65 / 0.2e1 - t728 * t233 / 0.2e1;
t349 = qJD(4) - t354;
t334 = qJD(5) + t349;
t710 = t334 * t728 + t505 * t729 + t762;
t709 = t216 * t732 + t730 * t334 + t764;
t562 = t459 * t386;
t289 = t385 * t454 - t562;
t546 = qJD(3) * t454;
t760 = pkin(2) * t546 - t289;
t757 = t710 / 0.2e1;
t625 = pkin(5) * t216;
t755 = t233 * t730 + t65 * t732 + t66 * t731;
t446 = t460 * pkin(2);
t435 = t446 + pkin(1);
t405 = t435 * qJD(1);
t253 = -pkin(3) * t354 - pkin(9) * t355 - t405;
t286 = t454 * t364 - t562;
t265 = pkin(9) * t449 + t286;
t170 = t458 * t253 - t265 * t453;
t754 = t170 * mrSges(5,1);
t171 = t253 * t453 + t265 * t458;
t753 = t171 * mrSges(5,2);
t367 = t617 * t453;
t445 = t458 * pkin(10);
t368 = t431 * t458 + t445;
t278 = t452 * t367 + t457 * t368;
t716 = -qJD(5) * t278 + t452 * t766 + t457 * t768;
t541 = qJD(5) * t457;
t542 = qJD(5) * t452;
t715 = t367 * t541 - t368 * t542 + t452 * t768 - t457 * t766;
t406 = t462 * t453;
t408 = pkin(9) * t458 + t445;
t310 = t452 * t406 + t457 * t408;
t714 = -qJD(5) * t310 + t452 * t765 + t457 * t767;
t713 = t406 * t541 - t408 * t542 + t452 * t767 - t457 * t765;
t752 = qJ(6) * t216;
t598 = t355 * mrSges(4,3);
t751 = -mrSges(4,1) * t449 - mrSges(5,1) * t306 + mrSges(5,2) * t307 + t598;
t750 = qJ(6) * t763 - qJD(6) * t483;
t544 = qJD(4) * t453;
t438 = pkin(4) * t544;
t747 = -pkin(5) * t763 + t438;
t330 = pkin(4) * t589;
t746 = -t330 + t760;
t745 = -t355 * pkin(5) - qJ(6) * t748 - qJD(6) * t378;
t294 = -qJD(2) * t377 - t474;
t543 = qJD(4) * t458;
t517 = t379 * t543;
t477 = t294 * t453 + t517;
t595 = t453 * mrSges(5,2);
t616 = mrSges(5,1) * t458;
t744 = t595 - t616;
t593 = qJDD(1) * pkin(1);
t344 = -pkin(2) * t392 - t593;
t135 = -pkin(3) * t242 - pkin(9) * t241 + t344;
t371 = t393 * pkin(7);
t298 = qJDD(2) * pkin(2) - pkin(8) * t393 - t371;
t370 = t392 * pkin(7);
t305 = pkin(8) * t392 + t370;
t152 = t454 * t298 + t459 * t305 + t364 * t545 + t386 * t546;
t143 = pkin(9) * t447 + t152;
t35 = t453 * t135 + t458 * t143 + t253 * t543 - t265 * t544;
t36 = -qJD(4) * t171 + t458 * t135 - t143 * t453;
t743 = t35 * t458 - t36 * t453;
t133 = -pkin(10) * t307 + t170;
t115 = pkin(4) * t349 + t133;
t134 = pkin(10) * t306 + t171;
t128 = t457 * t134;
t57 = t115 * t452 + t128;
t712 = qJ(6) * t505;
t34 = t57 + t712;
t742 = mrSges(6,3) * t57 + mrSges(7,3) * t34;
t725 = t349 * Ifges(5,3);
t726 = t306 * Ifges(5,6);
t739 = t307 * Ifges(5,5) + t216 * t730 + t334 * t727 + t505 * t728 + t725 + t726;
t451 = qJ(2) + qJ(3);
t442 = sin(t451);
t738 = (t616 - t771) * t442;
t444 = cos(t451);
t737 = -t444 * mrSges(4,1) + (mrSges(4,2) + t759) * t442;
t700 = t444 * pkin(3) + t442 * pkin(9);
t626 = pkin(4) * t458;
t433 = pkin(3) + t626;
t702 = t444 * t433 - t442 * t462;
t397 = pkin(5) * t443 + t626;
t384 = pkin(3) + t397;
t448 = -qJ(6) + t462;
t705 = t444 * t384 - t442 * t448;
t736 = -m(5) * t700 - m(6) * t702 - m(7) * t705;
t675 = m(6) * pkin(4);
t668 = t174 / 0.2e1;
t667 = t175 / 0.2e1;
t659 = t238 / 0.2e1;
t735 = t392 / 0.2e1;
t638 = t460 / 0.2e1;
t665 = -t505 / 0.2e1;
t724 = t449 * Ifges(4,5);
t723 = t449 * Ifges(4,6);
t722 = t460 * Ifges(3,2);
t721 = t750 + t715;
t720 = t745 + t716;
t719 = t750 + t713;
t718 = t745 + t714;
t717 = t675 + mrSges(5,1);
t264 = -pkin(3) * t449 - t285;
t494 = mrSges(5,1) * t453 + mrSges(5,2) * t458;
t711 = t264 * t494;
t271 = t483 * t379;
t226 = t330 + t286;
t708 = -t226 + t747;
t707 = t746 + t747;
t284 = pkin(3) * t377 - pkin(9) * t379 - t435;
t311 = t407 * t454 - t409 * t459;
t203 = t458 * t284 - t311 * t453;
t584 = t379 * t458;
t159 = pkin(4) * t377 - pkin(10) * t584 + t203;
t300 = t458 * t311;
t204 = t453 * t284 + t300;
t585 = t379 * t453;
t181 = -pkin(10) * t585 + t204;
t89 = t452 * t159 + t457 * t181;
t704 = t459 * t407 + t409 * t454;
t703 = t438 + t746;
t699 = t438 - t226;
t697 = t233 * t727 + t65 * t730 + t66 * t728;
t622 = pkin(7) * t460;
t623 = pkin(7) * t455;
t696 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t549) * t622 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t548) * t623;
t456 = sin(qJ(1));
t461 = cos(qJ(1));
t574 = t444 * t461;
t327 = -t441 * t574 + t443 * t456;
t328 = t441 * t456 + t443 * t574;
t695 = t327 * t734 - t328 * t733;
t575 = t444 * t456;
t325 = t441 * t575 + t443 * t461;
t326 = t441 * t461 - t443 * t575;
t694 = -t325 * t734 + t326 * t733;
t693 = t370 * t460 + t371 * t455;
t692 = mrSges(6,1) * t441 - t443 * t733;
t691 = g(1) * t461 + g(2) * t456;
t126 = t452 * t134;
t56 = t457 * t115 - t126;
t33 = t56 - t752;
t31 = pkin(5) * t334 + t33;
t690 = -mrSges(6,3) * t56 - mrSges(7,3) * t31;
t688 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t687 = m(7) + m(5) + m(6) + m(4);
t484 = -t433 * t442 - t444 * t462;
t486 = -t384 * t442 - t444 * t448;
t632 = pkin(2) * t455;
t686 = -m(7) * (t486 - t632) - m(6) * (t484 - t632) - m(5) * (-pkin(3) * t442 - t632) + t738;
t685 = -m(6) * t484 - m(7) * t486 + t738;
t532 = t442 * t595;
t578 = t442 * t461;
t684 = -t461 * t532 + t574 * t759 + t578 * t758;
t579 = t442 * t456;
t683 = -t456 * t532 + t575 * t759 + t579 * t758;
t682 = t737 + (t744 - t758 + t771) * t444;
t628 = pkin(4) * t453;
t396 = pkin(5) * t441 + t628;
t681 = -m(6) * t628 - m(7) * t396;
t680 = t36 * mrSges(5,1) - t35 * mrSges(5,2);
t404 = -mrSges(3,1) * t460 + mrSges(3,2) * t455;
t679 = m(3) * pkin(1) + mrSges(2,1) - t404 - t737;
t109 = mrSges(5,1) * t238 - mrSges(5,3) * t174;
t110 = -mrSges(5,2) * t238 + mrSges(5,3) * t175;
t613 = mrSges(5,3) * t306;
t244 = -mrSges(5,2) * t349 + t613;
t612 = mrSges(5,3) * t307;
t245 = mrSges(5,1) * t349 - t612;
t678 = m(5) * ((-t170 * t458 - t171 * t453) * qJD(4) + t743) - t245 * t543 - t244 * t544 + t458 * t110 - t453 * t109;
t26 = pkin(4) * t238 - pkin(10) * t174 + t36;
t29 = pkin(10) * t175 + t35;
t8 = -qJD(5) * t57 + t457 * t26 - t29 * t452;
t2 = pkin(5) * t233 - qJ(6) * t65 - qJD(6) * t216 + t8;
t7 = t115 * t541 - t134 * t542 + t452 * t26 + t457 * t29;
t4 = qJ(6) * t66 + qJD(6) * t505 + t7;
t677 = t8 * mrSges(6,1) + t2 * mrSges(7,1) - t7 * mrSges(6,2) - t4 * mrSges(7,2);
t674 = m(7) * pkin(5);
t673 = t65 / 0.2e1;
t672 = t66 / 0.2e1;
t671 = Ifges(5,1) * t668 + Ifges(5,4) * t667 + Ifges(5,5) * t659;
t664 = t505 / 0.2e1;
t662 = -t216 / 0.2e1;
t661 = t216 / 0.2e1;
t660 = t233 / 0.2e1;
t651 = -t306 / 0.2e1;
t650 = -t307 / 0.2e1;
t649 = t307 / 0.2e1;
t648 = -t334 / 0.2e1;
t647 = t334 / 0.2e1;
t646 = -t349 / 0.2e1;
t644 = t354 / 0.2e1;
t642 = t355 / 0.2e1;
t631 = pkin(2) * t459;
t630 = pkin(4) * t307;
t627 = pkin(4) * t457;
t619 = g(3) * t442;
t611 = mrSges(6,3) * t505;
t610 = mrSges(6,3) * t216;
t609 = mrSges(7,3) * t505;
t608 = mrSges(7,3) * t216;
t607 = Ifges(3,4) * t455;
t606 = Ifges(3,4) * t460;
t605 = Ifges(5,4) * t307;
t604 = Ifges(5,4) * t453;
t603 = Ifges(5,4) * t458;
t600 = t285 * mrSges(4,3);
t597 = t355 * Ifges(4,4);
t594 = qJ(6) * t378;
t591 = t294 * t458;
t197 = Ifges(5,2) * t306 + Ifges(5,6) * t349 + t605;
t572 = t453 * t197;
t571 = t453 * t456;
t570 = t453 * t461;
t567 = t456 * t396;
t566 = t456 * t458;
t304 = Ifges(5,4) * t306;
t198 = Ifges(5,1) * t307 + Ifges(5,5) * t349 + t304;
t564 = t458 * t198;
t563 = t458 * t461;
t69 = t457 * t133 - t126;
t547 = qJD(2) * t455;
t530 = pkin(2) * t547;
t528 = pkin(4) * t542;
t527 = pkin(4) * t541;
t520 = Ifges(5,5) * t174 + Ifges(5,6) * t175 + Ifges(5,3) * t238;
t519 = qJD(2) * t463;
t514 = t564 / 0.2e1;
t21 = -t66 * mrSges(7,1) + t65 * mrSges(7,2);
t511 = -t544 / 0.2e1;
t509 = t540 / 0.2e1;
t68 = -t133 * t452 - t128;
t88 = t457 * t159 - t181 * t452;
t295 = t449 * t379;
t201 = pkin(3) * t295 - pkin(9) * t294 + t530;
t390 = t455 * t519;
t391 = t460 * t519;
t219 = qJD(3) * t704 + t390 * t459 + t391 * t454;
t506 = t458 * t201 - t219 * t453;
t277 = t457 * t367 - t368 * t452;
t308 = t457 * t406 - t408 * t452;
t254 = pkin(4) * t585 - t704;
t497 = mrSges(3,1) * t455 + mrSges(3,2) * t460;
t495 = mrSges(4,1) * t442 + mrSges(4,2) * t444;
t492 = Ifges(5,1) * t458 - t604;
t491 = t607 + t722;
t490 = -Ifges(5,2) * t453 + t603;
t489 = Ifges(3,5) * t460 - Ifges(3,6) * t455;
t488 = Ifges(5,5) * t458 - Ifges(5,6) * t453;
t329 = pkin(5) * t483 - t433;
t153 = t298 * t459 - t454 * t305 - t364 * t546 + t386 * t545;
t478 = pkin(1) * t497;
t342 = -t444 * t570 + t566;
t340 = t444 * t571 + t563;
t476 = t379 * t544 - t591;
t475 = t455 * (Ifges(3,1) * t460 - t607);
t46 = -pkin(10) * t591 + pkin(4) * t295 + (-t300 + (pkin(10) * t379 - t284) * t453) * qJD(4) + t506;
t80 = t453 * t201 + t458 * t219 + t284 * t543 - t311 * t544;
t59 = -pkin(10) * t477 + t80;
t11 = t159 * t541 - t181 * t542 + t452 * t46 + t457 * t59;
t144 = -pkin(3) * t447 - t153;
t202 = -pkin(4) * t306 + t264;
t82 = -pkin(4) * t175 + t144;
t12 = -qJD(5) * t89 - t452 * t59 + t457 * t46;
t220 = qJD(3) * t311 + t390 * t454 - t459 * t391;
t468 = t677 + t697;
t141 = pkin(4) * t477 + t220;
t117 = -pkin(5) * t505 + qJD(6) + t202;
t259 = t354 * Ifges(4,2) + t597 + t723;
t346 = Ifges(4,4) * t354;
t260 = t355 * Ifges(4,1) + t346 + t724;
t27 = -pkin(5) * t66 + qJDD(6) + t82;
t72 = Ifges(5,4) * t174 + Ifges(5,2) * t175 + Ifges(5,6) * t238;
t464 = -t31 * (mrSges(7,1) * t355 + mrSges(7,3) * t252) - t56 * (mrSges(6,1) * t355 + mrSges(6,3) * t252) - t34 * (-mrSges(7,2) * t355 - mrSges(7,3) * t251) - t57 * (-mrSges(6,2) * t355 - mrSges(6,3) * t251) + t405 * (mrSges(4,1) * t355 + mrSges(4,2) * t354) + (-mrSges(6,1) * t763 + mrSges(6,2) * t748) * t202 + (-mrSges(7,1) * t763 + mrSges(7,2) * t748) * t117 + (t514 + t711) * qJD(4) + (t755 / 0.2e1 + t82 * mrSges(6,2) + t27 * mrSges(7,2) - mrSges(6,3) * t8 - mrSges(7,3) * t2 + t730 * t660 + t731 * t672 + t732 * t673) * t378 - t354 * t711 + t355 * t753 + t286 * t598 + t354 * t600 + ((-t544 + t589) * t171 + (-t543 + t588) * t170 + t743) * mrSges(5,3) + t144 * t744 + (t306 * t490 + t307 * t492 + t349 * t488) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t355 + t260 + t346 + t564) * t354 / 0.2e1 - (Ifges(4,1) * t354 - t597 + t739) * t355 / 0.2e1 - t152 * mrSges(4,2) + t458 * t72 / 0.2e1 - t449 * (Ifges(4,5) * t354 - Ifges(4,6) * t355) / 0.2e1 + Ifges(4,3) * t447 + (mrSges(6,1) * t82 + mrSges(7,1) * t27 - mrSges(6,3) * t7 - mrSges(7,3) * t4 - t728 * t660 - t729 * t672 - t731 * t673 + t761) * t483 + t259 * t642 + t572 * t644 + (Ifges(5,3) * t355 + t354 * t488) * t646 - t355 * t754 + Ifges(4,5) * t241 + Ifges(4,6) * t242 + (-t251 * t728 - t252 * t730 + t355 * t727) * t648 - (t647 * t728 + t661 * t731 + t664 * t729 + t742) * t293 + (-t251 * t729 - t252 * t731 + t355 * t728) * t665 - (t647 * t730 + t661 * t732 + t731 * t664 + t690) * t292 + (-t251 * t731 - t252 * t732 + t355 * t730) * t662 + t197 * t511 + t709 * (-t292 / 0.2e1 + t252 / 0.2e1) + t710 * (-t293 / 0.2e1 + t251 / 0.2e1) + (Ifges(5,5) * t355 + t354 * t492) * t650 + (Ifges(5,6) * t355 + t354 * t490) * t651 + (Ifges(5,5) * t453 + Ifges(5,6) * t458) * t659 + (Ifges(5,2) * t458 + t604) * t667 + (Ifges(5,1) * t453 + t603) * t668 + t453 * t671 + t153 * mrSges(4,1);
t437 = Ifges(3,4) * t548;
t434 = -pkin(3) - t631;
t432 = pkin(5) + t627;
t421 = pkin(9) * t574;
t420 = pkin(9) * t575;
t401 = -t433 - t631;
t366 = t483 * qJ(6);
t353 = Ifges(3,1) * t549 + Ifges(3,5) * qJD(2) + t437;
t352 = Ifges(3,6) * qJD(2) + qJD(1) * t491;
t343 = t444 * t563 + t571;
t341 = -t444 * t566 + t570;
t315 = -mrSges(4,2) * t449 + mrSges(4,3) * t354;
t314 = t329 - t631;
t280 = -mrSges(4,1) * t354 + mrSges(4,2) * t355;
t270 = t378 * t379;
t263 = -t366 + t310;
t262 = t308 - t594;
t240 = -t366 + t278;
t239 = t277 - t594;
t225 = -mrSges(4,2) * t447 + mrSges(4,3) * t242;
t224 = mrSges(4,1) * t447 - mrSges(4,3) * t241;
t188 = pkin(5) * t270 + t254;
t187 = mrSges(6,1) * t334 - t610;
t186 = mrSges(7,1) * t334 - t608;
t185 = -mrSges(6,2) * t334 + t611;
t184 = -mrSges(7,2) * t334 + t609;
t169 = t630 + t625;
t125 = -mrSges(6,1) * t505 + mrSges(6,2) * t216;
t124 = -mrSges(7,1) * t505 + mrSges(7,2) * t216;
t102 = t271 * t689 - t378 * t294;
t101 = -t293 * t379 - t294 * t483;
t90 = -mrSges(5,1) * t175 + mrSges(5,2) * t174;
t81 = -qJD(4) * t204 + t506;
t74 = -qJ(6) * t270 + t89;
t67 = pkin(5) * t377 + qJ(6) * t271 + t88;
t54 = -pkin(5) * t102 + t141;
t42 = -mrSges(6,2) * t233 + mrSges(6,3) * t66;
t41 = -mrSges(7,2) * t233 + mrSges(7,3) * t66;
t40 = mrSges(6,1) * t233 - mrSges(6,3) * t65;
t39 = mrSges(7,1) * t233 - mrSges(7,3) * t65;
t38 = t69 - t752;
t37 = t68 - t712;
t22 = -mrSges(6,1) * t66 + mrSges(6,2) * t65;
t10 = qJ(6) * t102 - qJD(6) * t270 + t11;
t9 = pkin(5) * t295 - qJ(6) * t101 + qJD(6) * t271 + t12;
t1 = [(t344 * mrSges(4,2) - t153 * mrSges(4,3) + Ifges(4,1) * t241 + Ifges(4,4) * t242 + Ifges(4,5) * t447 + t144 * t494 + t198 * t511 + t488 * t659 + t490 * t667 + t492 * t668) * t379 + t306 * (-Ifges(5,4) * t476 - Ifges(5,2) * t477) / 0.2e1 + (t392 * t622 + t393 * t623 + t693) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t693) + (t353 * t638 + t489 * qJD(2) / 0.2e1 - t696) * qJD(2) + (t520 + t697) * t377 / 0.2e1 + (t460 * t606 + t475) * t509 + t27 * (mrSges(7,1) * t270 - mrSges(7,2) * t271) + (-t101 * t56 + t102 * t57 - t270 * t7 + t271 * t8) * mrSges(6,3) + (-t101 * t31 + t102 * t34 + t2 * t271 - t270 * t4) * mrSges(7,3) + t82 * (mrSges(6,1) * t270 - mrSges(6,2) * t271) - (-m(4) * t153 + m(5) * t144 - t224 + t90) * t704 + m(5) * (t170 * t81 + t171 * t80 + t203 * t36 + t204 * t35) + (t170 * t476 - t171 * t477 - t35 * t585 - t36 * t584) * mrSges(5,3) + (Ifges(3,4) * t393 + Ifges(3,2) * t392) * t638 + (t739 / 0.2e1 - t286 * mrSges(4,3) + t754 - t753 + t725 / 0.2e1 + t726 / 0.2e1 - t723 / 0.2e1 - t405 * mrSges(4,1) - Ifges(4,4) * t642 - Ifges(4,2) * t644 + Ifges(5,5) * t649 + t56 * mrSges(6,1) + t31 * mrSges(7,1) - t57 * mrSges(6,2) - t34 * mrSges(7,2) - t259 / 0.2e1 + t728 * t664 + t730 * t661 + t727 * t647) * t295 - t352 * t547 / 0.2e1 + t102 * t757 + (-m(4) * t285 + m(5) * t264 + t751) * t220 + t349 * (-Ifges(5,5) * t476 - Ifges(5,6) * t477) / 0.2e1 + t10 * t184 + t393 * t606 / 0.2e1 - t72 * t585 / 0.2e1 - t197 * t517 / 0.2e1 + (t514 - t600 - t572 / 0.2e1 + t724 / 0.2e1 - t405 * mrSges(4,2) + Ifges(4,1) * t642 + Ifges(4,4) * t644 + t260 / 0.2e1) * t294 + m(4) * (t152 * t311 + t219 * t286 - t344 * t435 - t405 * t530) + (t101 * t731 + t102 * t729) * t664 + (-t270 * t729 - t271 * t731) * t672 + (-t270 * t731 - t271 * t732) * t673 + (t101 * t732 + t102 * t731) * t661 + (-t341 * mrSges(5,1) - t340 * mrSges(5,2) + t734 * t326 + t733 * t325 + (t463 * t687 + t681 + t688) * t461 + (-m(7) * (-t435 - t705) - m(5) * (-t435 - t700) - m(6) * (-t435 - t702) + m(4) * t435 + t679) * t456) * g(1) + t491 * t735 - t404 * t593 + (-mrSges(3,1) * t623 - mrSges(3,2) * t622 + 0.2e1 * Ifges(3,6) * t638) * qJDD(2) + (Ifges(3,1) * t393 + Ifges(3,4) * t735 + Ifges(3,5) * qJDD(2) - t509 * t722) * t455 + t280 * t530 + t270 * t761 + t709 * t101 / 0.2e1 - t478 * t540 + (-Ifges(5,1) * t476 - Ifges(5,4) * t477) * t649 + Ifges(2,3) * qJDD(1) + t141 * t125 + t54 * t124 + t117 * (-mrSges(7,1) * t102 + mrSges(7,2) * t101) + t88 * t40 + t89 * t42 + t74 * t41 + t67 * t39 + (-t571 * t675 - m(7) * t567 - t343 * mrSges(5,1) - t342 * mrSges(5,2) - t687 * (t461 * t435 - t456 * t463) + t734 * t328 + t733 * t327 + t688 * t456 + (-t679 + t736) * t461) * g(2) + t264 * (mrSges(5,1) * t477 - mrSges(5,2) * t476) - t435 * (-mrSges(4,1) * t242 + mrSges(4,2) * t241) - pkin(1) * (-mrSges(3,1) * t392 + mrSges(3,2) * t393) + t11 * t185 + t9 * t186 + t12 * t187 + t188 * t21 + t202 * (-mrSges(6,1) * t102 + mrSges(6,2) * t101) + t203 * t109 + t204 * t110 + m(7) * (t10 * t34 + t117 * t54 + t188 * t27 + t2 * t67 + t31 * t9 + t4 * t74) + m(6) * (t11 * t57 + t12 * t56 + t141 * t202 + t254 * t82 + t7 * t89 + t8 * t88) + (-t270 * t728 - t271 * t730) * t660 + (t101 * t730 + t102 * t728) * t647 + (t344 * mrSges(4,1) - t152 * mrSges(4,3) - Ifges(4,4) * t241 + Ifges(5,5) * t668 - Ifges(4,2) * t242 - Ifges(4,6) * t447 + Ifges(5,6) * t667 + Ifges(5,3) * t659 + t660 * t727 + t672 * t728 + t673 * t730 + t677 + t680) * t377 + t80 * t244 + t81 * t245 - t755 * t271 / 0.2e1 + t254 * t22 + t311 * t225 + t219 * t315 + t584 * t671; t720 * t186 + (t117 * t707 + t2 * t239 + t240 * t4 + t27 * t314 + t31 * t720 + t34 * t721) * m(7) + t721 * t184 + t715 * t185 + (t202 * t703 + t277 * t8 + t278 * t7 + t401 * t82 + t56 * t716 + t57 * t715) * m(6) + t716 * t187 - (-Ifges(3,2) * t549 + t353 + t437) * t548 / 0.2e1 + (m(4) * t632 + t495 + t497) * t691 + (t529 - t290) * t315 + t703 * t125 + (-m(5) * t420 + t456 * t686 + t683) * g(2) + (-m(5) * t421 + t461 * t686 + t684) * g(1) + t464 + ((t152 * t454 + t153 * t459 + (-t285 * t454 + t286 * t459) * qJD(3)) * pkin(2) + t285 * t289 - t286 * t290 + t405 * t531) * m(4) + (t696 + (t478 - t475 / 0.2e1) * qJD(1)) * qJD(1) + t352 * t549 / 0.2e1 + (-m(6) * (t446 + t702) - m(7) * (t446 + t705) - m(5) * (t446 + t700) - m(4) * t446 + t404 + t682) * g(3) - t280 * t531 + (-t170 * t191 - t171 * t192 - t264 * t289 + t144 * t434 + (t264 * t454 + (-t170 * t453 + t171 * t458) * t459) * qJD(3) * pkin(2)) * m(5) + t678 * t431 + t707 * t124 - t489 * t540 / 0.2e1 + Ifges(3,3) * qJDD(2) + t434 * t90 + t401 * t22 + Ifges(3,6) * t392 + Ifges(3,5) * t393 - t370 * mrSges(3,2) - t371 * mrSges(3,1) + t239 * t39 + t240 * t41 + t277 * t40 + t278 * t42 + t769 * t245 + t770 * t244 + t751 * t760 + t314 * t21 + t224 * t631 + t225 * t633; t718 * t186 + (t117 * t708 + t2 * t262 + t263 * t4 + t27 * t329 + t31 * t718 + t34 * t719) * m(7) + t719 * t184 + t713 * t185 + (t202 * t699 + t308 * t8 + t310 * t7 - t433 * t82 + t56 * t714 + t57 * t713) * m(6) + t714 * t187 + t691 * t495 + (-pkin(3) * t144 - t170 * t194 - t171 * t195 - t264 * t286) * m(5) - t751 * t286 + t699 * t125 + (-m(5) * (-pkin(3) * t579 + t420) + t685 * t456 + t683) * g(2) + (-m(5) * (-pkin(3) * t578 + t421) + t685 * t461 + t684) * g(1) + t464 + t678 * pkin(9) + t708 * t124 - pkin(3) * t90 + (t682 + t736) * g(3) - t433 * t22 - t195 * t244 - t194 * t245 + t262 * t39 + t263 * t41 + t308 * t40 + t310 * t42 - t285 * t315 + t329 * t21; (-m(7) * (-t397 * t461 - t444 * t567) - mrSges(5,2) * t341 + t717 * t340 + t694) * g(2) + (-m(7) * (-t396 * t574 + t397 * t456) + mrSges(5,2) * t343 - t717 * t342 + t695) * g(1) + (mrSges(7,1) * t441 + t494 - t681 + t692) * t619 + (t42 + t41) * pkin(4) * t452 + t680 + (-Ifges(5,2) * t307 + t198 + t304) * t651 + (t613 - t244) * t170 + (-t528 - t68) * t187 + (-t38 + t527) * t184 - m(6) * (t202 * t630 + t56 * t68 + t57 * t69) + (-t528 - t37) * t186 + (t612 + t245) * t171 + (-t202 * mrSges(6,2) - t117 * mrSges(7,2) + t648 * t730 + t662 * t732 + t665 * t731 - t690) * t505 + (t527 - t69) * t185 + t709 * t665 - t125 * t630 + t468 + t520 + t432 * t39 + (Ifges(5,5) * t306 - Ifges(5,6) * t307) * t646 - t169 * t124 + (-t202 * mrSges(6,1) - t117 * mrSges(7,1) - t728 * t648 - t731 * t662 - t729 * t665 + t742 + t757) * t216 - t264 * (mrSges(5,1) * t307 + mrSges(5,2) * t306) + t40 * t627 + t197 * t649 + (Ifges(5,1) * t306 - t605) * t650 + (t452 * t7 + t457 * t8 + (-t452 * t56 + t457 * t57) * qJD(5)) * t675 + (-t117 * t169 - t31 * t37 - t34 * t38 + t2 * t432 + (t4 * t452 + (-t31 * t452 + t34 * t457) * qJD(5)) * pkin(4)) * m(7); t31 * t609 - t33 * t184 - t124 * t625 + t468 + pkin(5) * t39 - t202 * (mrSges(6,1) * t216 + mrSges(6,2) * t505) + t2 * t674 + (t505 * t732 - t762) * t662 + t710 * t661 + (-t216 * t728 + t505 * t730) * t648 + (-(-mrSges(7,1) - t674) * t441 + t692) * t619 + (t610 + t187) * t57 + (t611 - t185) * t56 + (-m(7) * t625 - mrSges(7,1) * t216 - mrSges(7,2) * t505) * t117 + (t608 - m(7) * (-t31 + t33) + t186) * t34 + (t325 * t674 + t694) * g(2) + (-t327 * t674 + t695) * g(1) + (-t216 * t729 + t709 + t764) * t665; -t505 * t184 + t216 * t186 + (g(3) * t444 + t31 * t216 - t34 * t505 - t442 * t691 + t27) * m(7) + t21;];
tau  = t1;
