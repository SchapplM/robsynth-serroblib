% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:12:24
% EndTime: 2019-03-09 23:14:27
% DurationCPUTime: 77.54s
% Computational Cost: add. (39162->1255), mult. (93291->1703), div. (0->0), fcn. (75889->18), ass. (0->508)
t477 = sin(qJ(2));
t610 = cos(pkin(6));
t562 = pkin(1) * t610;
t453 = t477 * t562;
t476 = sin(qJ(3));
t480 = cos(qJ(3));
t528 = pkin(3) * t476 - pkin(10) * t480;
t471 = sin(pkin(6));
t481 = cos(qJ(2));
t597 = t471 * t481;
t791 = -(t453 + (pkin(8) + t528) * t597) * qJD(1) + t528 * qJD(3);
t543 = t610 * qJD(1);
t534 = pkin(1) * t543;
t581 = qJD(1) * t471;
t559 = t477 * t581;
t379 = -pkin(8) * t559 + t481 * t534;
t504 = (pkin(2) * t477 - pkin(9) * t481) * t471;
t380 = qJD(1) * t504;
t288 = t480 * t379 + t476 * t380;
t266 = pkin(10) * t559 + t288;
t475 = sin(qJ(4));
t479 = cos(qJ(4));
t579 = qJD(3) * t476;
t570 = pkin(9) * t579;
t790 = t791 * t479 + (t266 + t570) * t475;
t424 = -pkin(3) * t480 - pkin(10) * t476 - pkin(2);
t575 = qJD(4) * t479;
t789 = -t479 * t266 + t424 * t575 + t475 * t791;
t583 = pkin(8) * t597 + t453;
t382 = t583 * qJD(1);
t445 = t543 + qJD(2);
t335 = t445 * pkin(9) + t382;
t344 = (-pkin(2) * t481 - pkin(9) * t477 - pkin(1)) * t581;
t250 = -t476 * t335 + t344 * t480;
t351 = t445 * t480 - t476 * t559;
t535 = t480 * t559;
t352 = t445 * t476 + t535;
t279 = pkin(3) * t352 - pkin(10) * t351;
t177 = -t250 * t475 + t479 * t279;
t473 = -qJ(5) - pkin(10);
t544 = qJD(4) * t473;
t605 = t351 * t479;
t788 = -pkin(4) * t352 + qJ(5) * t605 - qJD(5) * t475 + t479 * t544 - t177;
t178 = t479 * t250 + t475 * t279;
t574 = qJD(5) * t479;
t606 = t351 * t475;
t787 = -qJ(5) * t606 - t475 * t544 + t178 - t574;
t590 = t480 * t481;
t342 = (t475 * t477 + t479 * t590) * t581;
t591 = t479 * t480;
t455 = pkin(9) * t591;
t558 = t481 * t581;
t536 = t476 * t558;
t786 = -pkin(4) * t536 + t342 * qJ(5) - t476 * t574 + (pkin(4) * t476 - qJ(5) * t591) * qJD(3) + (-t455 + (qJ(5) * t476 - t424) * t475) * qJD(4) + t790;
t341 = (-t475 * t590 + t477 * t479) * t581;
t593 = t476 * t479;
t785 = qJ(5) * t341 - (-pkin(9) * qJD(3) - qJ(5) * qJD(4)) * t593 - (-qJD(5) * t476 + (-pkin(9) * qJD(4) - qJ(5) * qJD(3)) * t480) * t475 - t789;
t470 = sin(pkin(12));
t472 = cos(pkin(12));
t410 = t470 * t479 + t472 * t475;
t253 = t410 * t351;
t395 = t410 * qJD(4);
t784 = t253 - t395;
t510 = t470 * t475 - t472 * t479;
t254 = t510 * t351;
t396 = t510 * qJD(4);
t783 = t254 - t396;
t257 = t341 * t472 - t342 * t470;
t576 = qJD(4) * t476;
t578 = qJD(3) * t480;
t290 = -t410 * t578 + t510 * t576;
t720 = t257 - t290;
t258 = t341 * t470 + t342 * t472;
t291 = -t395 * t476 - t510 * t578;
t719 = t258 - t291;
t716 = t536 - t579;
t734 = t787 * t470 + t472 * t788;
t733 = t470 * t788 - t787 * t472;
t732 = t470 * t785 + t472 * t786;
t731 = t470 * t786 - t472 * t785;
t718 = t475 * t578 + t476 * t575 + t341;
t623 = Ifges(4,4) * t352;
t428 = qJD(3) - t558;
t738 = t428 * Ifges(4,6);
t241 = t351 * Ifges(4,2) + t623 + t738;
t346 = qJD(4) - t351;
t337 = qJD(6) + t346;
t652 = t337 / 0.2e1;
t298 = -t352 * t475 + t428 * t479;
t299 = t352 * t479 + t428 * t475;
t197 = t298 * t470 + t299 * t472;
t664 = t197 / 0.2e1;
t542 = t472 * t298 - t299 * t470;
t666 = t542 / 0.2e1;
t474 = sin(qJ(6));
t478 = cos(qJ(6));
t125 = t197 * t478 + t474 * t542;
t672 = t125 / 0.2e1;
t762 = -t197 * t474 + t478 * t542;
t674 = t762 / 0.2e1;
t334 = -t445 * pkin(2) - t379;
t217 = -t351 * pkin(3) - t352 * pkin(10) + t334;
t251 = t480 * t335 + t476 * t344;
t220 = pkin(10) * t428 + t251;
t146 = t217 * t475 + t220 * t479;
t113 = qJ(5) * t298 + t146;
t108 = t470 * t113;
t145 = t479 * t217 - t220 * t475;
t112 = -qJ(5) * t299 + t145;
t99 = pkin(4) * t346 + t112;
t59 = t472 * t99 - t108;
t596 = t472 * t113;
t60 = t470 * t99 + t596;
t702 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t745 = Ifges(5,3) + Ifges(6,3);
t710 = t299 * Ifges(5,5) + t197 * Ifges(6,5) + t125 * Ifges(7,5) + t298 * Ifges(5,6) + Ifges(6,6) * t542 + t762 * Ifges(7,6) + t337 * Ifges(7,3) + t346 * t745;
t782 = -t241 / 0.2e1 + Ifges(6,5) * t664 + Ifges(7,5) * t672 + Ifges(6,6) * t666 + Ifges(7,6) * t674 + Ifges(7,3) * t652 + t710 / 0.2e1 - t702;
t526 = t480 * mrSges(4,1) - t476 * mrSges(4,2);
t707 = m(7) * (-pkin(11) + t473) - mrSges(7,3) + m(6) * t473 - mrSges(6,3) - m(5) * pkin(10) - mrSges(5,3);
t469 = qJ(4) + pkin(12);
t464 = cos(t469);
t630 = pkin(4) * t479;
t422 = pkin(5) * t464 + t630;
t416 = pkin(3) + t422;
t465 = qJ(6) + t469;
t458 = sin(t465);
t459 = cos(t465);
t461 = pkin(3) + t630;
t463 = sin(t469);
t524 = -mrSges(5,1) * t479 + mrSges(5,2) * t475;
t708 = m(5) * pkin(3) + m(6) * t461 + m(7) * t416 + mrSges(6,1) * t464 + mrSges(7,1) * t459 - mrSges(6,2) * t463 - mrSges(7,2) * t458 - t524;
t781 = t476 * t707 - t480 * t708 - t526;
t780 = -m(5) - m(4);
t779 = -pkin(5) * t716 + pkin(11) * t719 + t732;
t778 = -pkin(11) * t720 + t731;
t777 = -pkin(5) * t352 - pkin(11) * t783 + t734;
t776 = pkin(11) * t784 + t733;
t287 = -t476 * t379 + t380 * t480;
t265 = -pkin(3) * t559 - t287;
t462 = pkin(9) * t578;
t772 = t462 - t265;
t541 = t610 * qJDD(1);
t573 = qJD(1) * qJD(2);
t766 = -pkin(8) * t471 * t573 + pkin(1) * t541;
t572 = qJDD(1) * t471;
t767 = pkin(8) * t572 + qJD(2) * t534;
t300 = t477 * t766 + t481 * t767;
t444 = t541 + qJDD(2);
t274 = pkin(9) * t444 + t300;
t385 = (-qJDD(1) * t481 + t477 * t573) * t471;
t386 = (qJDD(1) * t477 + t481 * t573) * t471;
t568 = pkin(1) * t572;
t282 = pkin(2) * t385 - pkin(9) * t386 - t568;
t137 = t480 * t274 + t476 * t282 - t335 * t579 + t344 * t578;
t374 = qJDD(3) + t385;
t119 = pkin(10) * t374 + t137;
t255 = qJD(3) * t351 + t386 * t480 + t444 * t476;
t256 = -qJD(3) * t535 - t476 * t386 + t444 * t480 - t445 * t579;
t301 = -t477 * t767 + t481 * t766;
t275 = -t444 * pkin(2) - t301;
t135 = -t256 * pkin(3) - t255 * pkin(10) + t275;
t577 = qJD(4) * t475;
t48 = t479 * t119 + t475 * t135 + t217 * t575 - t220 * t577;
t49 = -qJD(4) * t146 - t119 * t475 + t479 * t135;
t166 = qJD(4) * t298 + t255 * t479 + t374 * t475;
t249 = qJDD(4) - t256;
t30 = pkin(4) * t249 - qJ(5) * t166 - qJD(5) * t299 + t49;
t167 = -qJD(4) * t299 - t255 * t475 + t374 * t479;
t39 = qJ(5) * t167 + qJD(5) * t298 + t48;
t12 = t470 * t30 + t472 * t39;
t747 = t12 * mrSges(6,2);
t11 = t472 * t30 - t39 * t470;
t748 = t11 * mrSges(6,1);
t771 = t49 * mrSges(5,1) - t48 * mrSges(5,2) - t747 + t748;
t769 = pkin(11) * t197;
t46 = pkin(5) * t346 + t59 - t769;
t749 = pkin(11) * t542;
t50 = t60 + t749;
t15 = t46 * t478 - t474 * t50;
t16 = t46 * t474 + t478 * t50;
t614 = t251 * mrSges(4,3);
t770 = t614 - t334 * mrSges(4,1) - t145 * mrSges(5,1) - t15 * mrSges(7,1) + t146 * mrSges(5,2) + t16 * mrSges(7,2) + t738 / 0.2e1;
t768 = mrSges(4,2) + t707;
t722 = pkin(4) * t718 + t772;
t715 = -t251 + (t577 - t606) * pkin(4);
t709 = -m(7) - m(6) + t780;
t764 = -pkin(2) * t709 + mrSges(3,1) - t781;
t523 = t475 * mrSges(5,1) + t479 * mrSges(5,2);
t763 = -t463 * mrSges(6,1) - t458 * mrSges(7,1) - t464 * mrSges(6,2) - t459 * mrSges(7,2) - t523;
t631 = pkin(4) * t475;
t421 = pkin(5) * t463 + t631;
t700 = mrSges(4,3) - mrSges(3,2) + m(7) * (pkin(9) + t421) + m(6) * (pkin(9) + t631);
t642 = cos(qJ(1));
t530 = t610 * t642;
t641 = sin(qJ(1));
t400 = t477 * t530 + t481 * t641;
t561 = t471 * t642;
t321 = t400 * t480 - t476 * t561;
t399 = t477 * t641 - t481 * t530;
t761 = t321 * t458 - t399 * t459;
t760 = -t321 * t459 - t399 * t458;
t759 = t321 * t463 - t399 * t464;
t758 = t321 * t464 + t399 * t463;
t757 = t321 * t475 - t399 * t479;
t756 = t321 * t479 + t399 * t475;
t93 = t166 * t472 + t167 * t470;
t6 = pkin(5) * t249 - pkin(11) * t93 + t11;
t92 = -t166 * t470 + t167 * t472;
t7 = pkin(11) * t92 + t12;
t2 = qJD(6) * t15 + t474 * t6 + t478 * t7;
t3 = -qJD(6) * t16 - t474 * t7 + t478 * t6;
t755 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t34 = qJD(6) * t762 + t474 * t92 + t478 * t93;
t691 = t34 / 0.2e1;
t35 = -qJD(6) * t125 - t474 * t93 + t478 * t92;
t690 = t35 / 0.2e1;
t681 = t92 / 0.2e1;
t680 = t93 / 0.2e1;
t670 = t166 / 0.2e1;
t669 = t167 / 0.2e1;
t238 = qJDD(6) + t249;
t663 = t238 / 0.2e1;
t661 = t249 / 0.2e1;
t660 = t255 / 0.2e1;
t659 = t256 / 0.2e1;
t645 = t374 / 0.2e1;
t412 = t479 * t424;
t310 = -qJ(5) * t593 + t412 + (-pkin(9) * t475 - pkin(4)) * t480;
t362 = t475 * t424 + t455;
t594 = t475 * t476;
t328 = -qJ(5) * t594 + t362;
t221 = t472 * t310 - t328 * t470;
t376 = t510 * t476;
t188 = -pkin(5) * t480 + pkin(11) * t376 + t221;
t222 = t470 * t310 + t472 * t328;
t375 = t410 * t476;
t193 = -pkin(11) * t375 + t222;
t114 = t188 * t478 - t193 * t474;
t744 = qJD(6) * t114 + t474 * t779 + t478 * t778;
t115 = t188 * t474 + t193 * t478;
t743 = -qJD(6) * t115 - t474 * t778 + t478 * t779;
t742 = Ifges(5,5) * t166 + Ifges(6,5) * t93 + Ifges(5,6) * t167 + Ifges(6,6) * t92 + t249 * t745;
t429 = t473 * t475;
t430 = t473 * t479;
t332 = t472 * t429 + t430 * t470;
t294 = -pkin(11) * t410 + t332;
t333 = t470 * t429 - t472 * t430;
t295 = -pkin(11) * t510 + t333;
t191 = t294 * t478 - t295 * t474;
t741 = qJD(6) * t191 + t474 * t777 + t478 * t776;
t192 = t294 * t474 + t295 * t478;
t740 = -qJD(6) * t192 - t474 * t776 + t478 * t777;
t739 = t428 * Ifges(4,5);
t693 = m(6) * pkin(4);
t737 = -mrSges(5,1) - t693;
t632 = pkin(4) * t472;
t460 = pkin(5) + t632;
t633 = pkin(4) * t470;
t388 = t460 * t474 + t478 * t633;
t64 = -t112 * t470 - t596;
t51 = t64 - t749;
t65 = t472 * t112 - t108;
t52 = t65 - t769;
t736 = -t388 * qJD(6) + t474 * t52 - t478 * t51;
t387 = t460 * t478 - t474 * t633;
t735 = t387 * qJD(6) - t474 * t51 - t478 * t52;
t728 = -pkin(5) * t784 + t715;
t283 = -t375 * t478 + t376 * t474;
t155 = qJD(6) * t283 + t290 * t474 + t291 * t478;
t171 = t257 * t474 + t258 * t478;
t727 = t155 - t171;
t284 = -t375 * t474 - t376 * t478;
t156 = -qJD(6) * t284 + t290 * t478 - t291 * t474;
t170 = t257 * t478 - t258 * t474;
t726 = t156 - t170;
t725 = pkin(5) * t720 + t722;
t724 = -qJD(4) * t362 + t790;
t723 = (-t479 * t579 - t480 * t577) * pkin(9) + t789;
t296 = Ifges(5,4) * t298;
t182 = t299 * Ifges(5,1) + t346 * Ifges(5,5) + t296;
t345 = Ifges(4,4) * t351;
t242 = t352 * Ifges(4,1) + t345 + t739;
t721 = t479 * t182 + t242;
t598 = t471 * t477;
t406 = -pkin(8) * t598 + t481 * t562;
t369 = -pkin(2) * t610 - t406;
t397 = t476 * t598 - t480 * t610;
t398 = t476 * t610 + t480 * t598;
t262 = t397 * pkin(3) - t398 * pkin(10) + t369;
t370 = pkin(9) * t610 + t583;
t584 = pkin(2) * t597 + pkin(9) * t598;
t636 = pkin(1) * t471;
t371 = -t584 - t636;
t281 = t480 * t370 + t476 * t371;
t264 = -pkin(10) * t597 + t281;
t174 = t475 * t262 + t479 * t264;
t717 = t475 * t576 - t479 * t578 + t342;
t138 = -t476 * t274 + t282 * t480 - t335 * t578 - t344 * t579;
t714 = t137 * t480 - t138 * t476;
t713 = -t475 * t49 + t479 * t48;
t8 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t238;
t712 = t8 + t742;
t625 = Ifges(3,4) * t477;
t696 = t471 ^ 2;
t711 = (t477 * (Ifges(3,1) * t481 - t625) / 0.2e1 - pkin(1) * (mrSges(3,1) * t477 + mrSges(3,2) * t481)) * t696;
t706 = mrSges(4,1) + t708;
t540 = mrSges(3,3) * t559;
t703 = -m(4) * t334 + mrSges(3,1) * t445 + mrSges(4,1) * t351 - mrSges(4,2) * t352 - t540;
t699 = pkin(9) * t780 - t700 + t763;
t651 = -t346 / 0.2e1;
t653 = -t337 / 0.2e1;
t656 = -t299 / 0.2e1;
t658 = -t298 / 0.2e1;
t665 = -t197 / 0.2e1;
t667 = -t542 / 0.2e1;
t673 = -t125 / 0.2e1;
t675 = -t762 / 0.2e1;
t697 = Ifges(5,5) * t656 + Ifges(6,5) * t665 + Ifges(7,5) * t673 + Ifges(5,6) * t658 + Ifges(6,6) * t667 + Ifges(7,6) * t675 + Ifges(7,3) * t653 + t651 * t745 + t702;
t694 = Ifges(7,4) * t691 + Ifges(7,2) * t690 + Ifges(7,6) * t663;
t692 = Ifges(7,1) * t691 + Ifges(7,4) * t690 + Ifges(7,5) * t663;
t689 = Ifges(6,4) * t680 + Ifges(6,2) * t681 + Ifges(6,6) * t661;
t688 = Ifges(6,1) * t680 + Ifges(6,4) * t681 + Ifges(6,5) * t661;
t618 = Ifges(7,4) * t125;
t67 = Ifges(7,2) * t762 + Ifges(7,6) * t337 + t618;
t687 = -t67 / 0.2e1;
t686 = t67 / 0.2e1;
t121 = Ifges(7,4) * t762;
t68 = Ifges(7,1) * t125 + Ifges(7,5) * t337 + t121;
t685 = -t68 / 0.2e1;
t684 = t68 / 0.2e1;
t80 = Ifges(5,4) * t166 + Ifges(5,2) * t167 + Ifges(5,6) * t249;
t683 = t80 / 0.2e1;
t682 = Ifges(5,1) * t670 + Ifges(5,4) * t669 + Ifges(5,5) * t661;
t105 = Ifges(6,4) * t197 + Ifges(6,2) * t542 + Ifges(6,6) * t346;
t679 = -t105 / 0.2e1;
t678 = t105 / 0.2e1;
t106 = Ifges(6,1) * t197 + Ifges(6,4) * t542 + Ifges(6,5) * t346;
t677 = -t106 / 0.2e1;
t676 = t106 / 0.2e1;
t671 = Ifges(4,1) * t660 + Ifges(4,4) * t659 + Ifges(4,5) * t645;
t668 = t182 / 0.2e1;
t657 = t298 / 0.2e1;
t655 = t299 / 0.2e1;
t650 = t346 / 0.2e1;
t648 = t351 / 0.2e1;
t646 = t352 / 0.2e1;
t640 = mrSges(6,3) * t59;
t639 = mrSges(6,3) * t60;
t638 = mrSges(7,3) * t15;
t637 = mrSges(7,3) * t16;
t635 = pkin(4) * t299;
t318 = -t398 * t475 - t479 * t597;
t634 = pkin(4) * t318;
t629 = pkin(9) * t480;
t466 = t476 * pkin(9);
t580 = qJD(2) * t481;
t556 = t471 * t580;
t317 = -qJD(3) * t397 + t480 * t556;
t557 = qJD(2) * t598;
t216 = qJD(4) * t318 + t317 * t479 + t475 * t557;
t316 = qJD(3) * t398 + t476 * t556;
t503 = -t398 * t479 + t475 * t597;
t381 = qJD(2) * t504;
t383 = t406 * qJD(2);
t189 = -t370 * t579 + t371 * t578 + t476 * t381 + t480 * t383;
t184 = pkin(10) * t557 + t189;
t384 = t583 * qJD(2);
t207 = t316 * pkin(3) - t317 * pkin(10) + t384;
t86 = -qJD(4) * t174 - t184 * t475 + t479 * t207;
t57 = pkin(4) * t316 - qJ(5) * t216 + qJD(5) * t503 + t86;
t215 = qJD(4) * t503 - t317 * t475 + t479 * t557;
t85 = t479 * t184 + t475 * t207 + t262 * t575 - t264 * t577;
t61 = qJ(5) * t215 + qJD(5) * t318 + t85;
t21 = t470 * t57 + t472 * t61;
t624 = Ifges(3,4) * t481;
t622 = Ifges(4,4) * t476;
t621 = Ifges(4,4) * t480;
t620 = Ifges(5,4) * t475;
t619 = Ifges(5,4) * t479;
t617 = t145 * mrSges(5,3);
t616 = t146 * mrSges(5,3);
t615 = t250 * mrSges(4,3);
t613 = t299 * Ifges(5,4);
t120 = -pkin(3) * t374 - t138;
t609 = t120 * t476;
t181 = t298 * Ifges(5,2) + t346 * Ifges(5,6) + t613;
t595 = t475 * t181;
t173 = t479 * t262 - t264 * t475;
t134 = pkin(4) * t397 + qJ(5) * t503 + t173;
t147 = qJ(5) * t318 + t174;
t76 = t470 * t134 + t472 * t147;
t589 = -mrSges(7,1) * t761 + mrSges(7,2) * t760;
t529 = t610 * t641;
t402 = -t477 * t529 + t481 * t642;
t560 = t471 * t641;
t325 = t402 * t480 + t476 * t560;
t401 = t477 * t642 + t481 * t529;
t243 = -t325 * t458 + t401 * t459;
t244 = t325 * t459 + t401 * t458;
t588 = t243 * mrSges(7,1) - t244 * mrSges(7,2);
t587 = (-t398 * t458 - t459 * t597) * mrSges(7,1) + (-t398 * t459 + t458 * t597) * mrSges(7,2);
t419 = pkin(4) * t594 + t466;
t582 = t642 * pkin(1) + pkin(8) * t560;
t566 = Ifges(4,5) * t255 + Ifges(4,6) * t256 + Ifges(4,3) * t374;
t565 = Ifges(3,5) * t386 - Ifges(3,6) * t385 + Ifges(3,3) * t444;
t564 = t402 * pkin(2) + t582;
t553 = t598 / 0.2e1;
t552 = -t595 / 0.2e1;
t47 = -t92 * mrSges(6,1) + t93 * mrSges(6,2);
t13 = -t35 * mrSges(7,1) + t34 * mrSges(7,2);
t545 = -t576 / 0.2e1;
t20 = -t470 * t61 + t472 * t57;
t75 = t472 * t134 - t147 * t470;
t280 = -t476 * t370 + t371 * t480;
t539 = mrSges(3,3) * t558;
t531 = -pkin(1) * t641 + pkin(8) * t561;
t263 = pkin(3) * t597 - t280;
t527 = mrSges(4,1) * t397 + mrSges(4,2) * t398;
t525 = mrSges(5,1) * t318 + mrSges(5,2) * t503;
t520 = Ifges(4,1) * t480 - t622;
t519 = Ifges(5,1) * t479 - t620;
t518 = Ifges(5,1) * t475 + t619;
t517 = -Ifges(4,2) * t476 + t621;
t516 = -Ifges(5,2) * t475 + t619;
t515 = Ifges(5,2) * t479 + t620;
t514 = Ifges(4,5) * t480 - Ifges(4,6) * t476;
t513 = Ifges(5,5) * t479 - Ifges(5,6) * t475;
t512 = Ifges(5,5) * t475 + Ifges(5,6) * t479;
t225 = t318 * t470 - t472 * t503;
t62 = pkin(5) * t397 - pkin(11) * t225 + t75;
t224 = t318 * t472 + t470 * t503;
t69 = pkin(11) * t224 + t76;
t22 = -t474 * t69 + t478 * t62;
t23 = t474 * t62 + t478 * t69;
t152 = t224 * t478 - t225 * t474;
t153 = t224 * t474 + t225 * t478;
t267 = -t325 * t475 + t401 * t479;
t308 = -t410 * t474 - t478 * t510;
t309 = t410 * t478 - t474 * t510;
t508 = pkin(9) * t401 + t564;
t506 = -t400 * pkin(2) + t531;
t190 = -t370 * t578 - t371 * t579 + t381 * t480 - t476 * t383;
t505 = t755 + t8;
t219 = -pkin(3) * t428 - t250;
t502 = t219 * t523;
t500 = (t481 * Ifges(3,2) + t625) * t471;
t320 = t400 * t476 + t480 * t561;
t324 = t402 * t476 - t480 * t560;
t496 = -g(1) * t324 - g(2) * t320 - g(3) * t397;
t198 = t263 - t634;
t492 = t445 * t471 * (Ifges(3,5) * t481 - Ifges(3,6) * t477);
t488 = -t399 * pkin(9) + t506;
t179 = -pkin(4) * t298 + qJD(5) + t219;
t185 = -pkin(3) * t557 - t190;
t82 = -pkin(4) * t167 + qJDD(5) + t120;
t126 = -pkin(4) * t215 + t185;
t439 = Ifges(3,4) * t558;
t403 = (-mrSges(3,1) * t481 + mrSges(3,2) * t477) * t471;
t378 = -t445 * mrSges(3,2) + t539;
t363 = pkin(5) * t510 - t461;
t361 = -t475 * t629 + t412;
t330 = Ifges(3,1) * t559 + t445 * Ifges(3,5) + t439;
t329 = Ifges(3,6) * t445 + qJD(1) * t500;
t311 = pkin(5) * t375 + t419;
t307 = mrSges(4,1) * t428 - mrSges(4,3) * t352;
t306 = -mrSges(4,2) * t428 + mrSges(4,3) * t351;
t268 = t325 * t479 + t401 * t475;
t260 = t325 * t464 + t401 * t463;
t259 = -t325 * t463 + t401 * t464;
t240 = t352 * Ifges(4,5) + t351 * Ifges(4,6) + t428 * Ifges(4,3);
t227 = mrSges(5,1) * t346 - mrSges(5,3) * t299;
t226 = -mrSges(5,2) * t346 + mrSges(5,3) * t298;
t213 = -qJD(6) * t309 - t395 * t478 + t396 * t474;
t212 = qJD(6) * t308 - t395 * t474 - t396 * t478;
t206 = -mrSges(4,2) * t374 + mrSges(4,3) * t256;
t205 = mrSges(4,1) * t374 - mrSges(4,3) * t255;
t202 = -mrSges(5,1) * t298 + mrSges(5,2) * t299;
t176 = mrSges(6,1) * t346 - mrSges(6,3) * t197;
t175 = -mrSges(6,2) * t346 + mrSges(6,3) * t542;
t172 = -mrSges(4,1) * t256 + mrSges(4,2) * t255;
t169 = -t253 * t474 - t254 * t478;
t168 = -t253 * t478 + t254 * t474;
t162 = pkin(5) * t197 + t635;
t157 = t255 * Ifges(4,4) + t256 * Ifges(4,2) + t374 * Ifges(4,6);
t143 = -pkin(5) * t224 + t198;
t142 = t215 * t470 + t216 * t472;
t140 = t215 * t472 - t216 * t470;
t127 = -mrSges(6,1) * t542 + mrSges(6,2) * t197;
t117 = -mrSges(5,2) * t249 + mrSges(5,3) * t167;
t116 = mrSges(5,1) * t249 - mrSges(5,3) * t166;
t107 = -pkin(5) * t542 + t179;
t101 = mrSges(7,1) * t337 - mrSges(7,3) * t125;
t100 = -mrSges(7,2) * t337 + mrSges(7,3) * t762;
t94 = -mrSges(5,1) * t167 + mrSges(5,2) * t166;
t78 = mrSges(6,1) * t249 - mrSges(6,3) * t93;
t77 = -mrSges(6,2) * t249 + mrSges(6,3) * t92;
t73 = -pkin(5) * t140 + t126;
t71 = -mrSges(7,1) * t762 + mrSges(7,2) * t125;
t54 = -qJD(6) * t153 + t140 * t478 - t142 * t474;
t53 = qJD(6) * t152 + t140 * t474 + t142 * t478;
t45 = -pkin(5) * t92 + t82;
t25 = -mrSges(7,2) * t238 + mrSges(7,3) * t35;
t24 = mrSges(7,1) * t238 - mrSges(7,3) * t34;
t17 = pkin(11) * t140 + t21;
t14 = pkin(5) * t316 - pkin(11) * t142 + t20;
t5 = -qJD(6) * t23 + t14 * t478 - t17 * t474;
t4 = qJD(6) * t22 + t14 * t474 + t17 * t478;
t1 = [(-Ifges(4,4) * t646 + Ifges(5,5) * t655 - Ifges(4,2) * t648 + Ifges(5,6) * t657 + t745 * t650 - t770 + t782) * t316 + (Ifges(5,6) * t669 + Ifges(5,5) * t670 + Ifges(7,6) * t690 + Ifges(7,5) * t691 + Ifges(6,5) * t680 + Ifges(6,6) * t681 + t745 * t661 - t137 * mrSges(4,3) + t712 / 0.2e1 - t157 / 0.2e1 + Ifges(7,3) * t663 - Ifges(4,2) * t659 - Ifges(4,4) * t660 - Ifges(4,6) * t645 + t755 + t771) * t397 + (Ifges(7,5) * t153 + Ifges(7,6) * t152) * t663 + (Ifges(7,5) * t53 + Ifges(7,6) * t54) * t652 + (Ifges(4,1) * t398 - Ifges(4,5) * t597) * t660 + (Ifges(4,1) * t317 + Ifges(4,5) * t557) * t646 + m(4) * (t137 * t281 + t138 * t280 + t189 * t251 + t190 * t250 + t275 * t369) + (t300 * t597 - t301 * t598 - t379 * t556 - t382 * t557 - t385 * t583 - t386 * t406) * mrSges(3,3) + (Ifges(3,4) * t386 - Ifges(3,2) * t385 + Ifges(3,6) * t444) * t597 / 0.2e1 + (t301 * t610 - t385 * t636 + t406 * t444) * mrSges(3,1) + (Ifges(3,1) * t386 - Ifges(3,4) * t385 + Ifges(3,5) * t444) * t553 + (t137 * t597 - t251 * t557 + t317 * t334) * mrSges(4,2) + (-m(7) * (t325 * t416 + t564) - t244 * mrSges(7,1) - t243 * mrSges(7,2) - m(3) * t582 - t402 * mrSges(3,1) - mrSges(3,3) * t560 - m(6) * (t325 * t461 + t564) - t260 * mrSges(6,1) - t259 * mrSges(6,2) - m(5) * (pkin(3) * t325 + t508) - t268 * mrSges(5,1) - t267 * mrSges(5,2) - m(4) * t508 - t325 * mrSges(4,1) - mrSges(2,1) * t642 + mrSges(2,2) * t641 - t700 * t401 + t768 * t324) * g(2) + (-m(7) * (-t321 * t416 + t506) - t760 * mrSges(7,1) - t761 * mrSges(7,2) - m(5) * (-pkin(3) * t321 + t488) + t756 * mrSges(5,1) - t757 * mrSges(5,2) - m(6) * (-t321 * t461 + t506) + t758 * mrSges(6,1) - t759 * mrSges(6,2) - m(3) * t531 + t400 * mrSges(3,1) - mrSges(3,3) * t561 + mrSges(2,1) * t641 + mrSges(2,2) * t642 - m(4) * t488 + t321 * mrSges(4,1) + t700 * t399 - t768 * t320) * g(1) + (-t11 * t225 + t12 * t224 + t140 * t60 - t142 * t59) * mrSges(6,3) + (Ifges(5,1) * t216 + Ifges(5,4) * t215) * t655 + (Ifges(7,1) * t153 + Ifges(7,4) * t152) * t691 + (Ifges(7,1) * t53 + Ifges(7,4) * t54) * t672 + (Ifges(6,4) * t225 + Ifges(6,2) * t224) * t681 + (Ifges(6,4) * t142 + Ifges(6,2) * t140) * t666 + (Ifges(6,1) * t225 + Ifges(6,4) * t224) * t680 + (Ifges(6,1) * t142 + Ifges(6,4) * t140) * t664 + (-t15 * t53 + t152 * t2 - t153 * t3 + t16 * t54) * mrSges(7,3) + (t471 * t330 + t696 * qJD(1) * (-Ifges(3,2) * t477 + t624)) * t580 / 0.2e1 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t696 + t300 * t583 + t301 * t406 + t382 * t383) + (-t300 * t610 - t386 * t636 - t444 * t583) * mrSges(3,2) - t329 * t557 / 0.2e1 + (Ifges(5,4) * t216 + Ifges(5,2) * t215) * t657 + t224 * t689 + t153 * t692 + t152 * t694 + t140 * t678 + t318 * t683 + t53 * t684 + t54 * t686 + t225 * t688 + (Ifges(5,5) * t216 + Ifges(6,5) * t142 + Ifges(5,6) * t215 + Ifges(6,6) * t140) * t650 - t566 * t597 / 0.2e1 + t138 * (-mrSges(4,1) * t597 - t398 * mrSges(4,3)) + t4 * t100 + t5 * t101 + (t240 * t553 + t492 / 0.2e1) * qJD(2) + t250 * (mrSges(4,1) * t557 - mrSges(4,3) * t317) + t386 * (Ifges(3,5) * t610 + (t477 * Ifges(3,1) + t624) * t471) / 0.2e1 + t610 * t565 / 0.2e1 - t385 * (Ifges(3,6) * t610 + t500) / 0.2e1 + t444 * (Ifges(3,3) * t610 + (Ifges(3,5) * t477 + Ifges(3,6) * t481) * t471) / 0.2e1 + (Ifges(7,4) * t153 + Ifges(7,2) * t152) * t690 + (Ifges(7,4) * t53 + Ifges(7,2) * t54) * t674 + (Ifges(4,4) * t398 - Ifges(4,6) * t597) * t659 + (Ifges(4,4) * t317 + Ifges(4,6) * t557) * t648 + t428 * (Ifges(4,5) * t317 + Ifges(4,3) * t557) / 0.2e1 + (Ifges(4,5) * t398 - Ifges(4,3) * t597) * t645 - t403 * t568 + t711 * t573 + (-m(3) * t379 - t703) * t384 + t107 * (-mrSges(7,1) * t54 + mrSges(7,2) * t53) + t76 * t77 + t75 * t78 + t73 * t71 + Ifges(2,3) * qJDD(1) + t23 * t25 + t22 * t24 + t383 * t378 + t369 * t172 + t216 * t668 + t398 * t671 + t142 * t676 - t120 * t525 + t275 * t527 + t126 * t127 + m(7) * (t107 * t73 + t143 * t45 + t15 * t5 + t16 * t4 + t2 * t23 + t22 * t3) + m(6) * (t11 * t75 + t12 * t76 + t126 * t179 + t198 * t82 + t20 * t59 + t21 * t60) + m(5) * (t120 * t263 + t145 * t86 + t146 * t85 + t173 * t49 + t174 * t48 + t185 * t219) + t143 * t13 + t45 * (-mrSges(7,1) * t152 + mrSges(7,2) * t153) + t173 * t116 + t174 * t117 + t21 * t175 + t20 * t176 + t179 * (-mrSges(6,1) * t140 + mrSges(6,2) * t142) + t198 * t47 + t185 * t202 + t215 * t181 / 0.2e1 + t219 * (-mrSges(5,1) * t215 + mrSges(5,2) * t216) + t82 * (-mrSges(6,1) * t224 + mrSges(6,2) * t225) + t85 * t226 + t86 * t227 + t263 * t94 + t280 * t205 + t281 * t206 + (-Ifges(5,5) * t503 + Ifges(6,5) * t225 + Ifges(5,6) * t318 + Ifges(6,6) * t224) * t661 + (-t145 * t216 + t146 * t215 + t318 * t48 + t49 * t503) * mrSges(5,3) + (-Ifges(5,1) * t503 + Ifges(5,4) * t318) * t670 + (-Ifges(5,4) * t503 + Ifges(5,2) * t318) * t669 - t503 * t682 + t189 * t306 + t190 * t307 + t317 * t242 / 0.2e1; (Ifges(6,3) * t650 - t614 + t782) * t579 + (t403 + t709 * t584 + (t781 * t481 + (-m(6) * t631 - m(7) * t421 - mrSges(4,3) + t763) * t477) * t471) * g(3) + t772 * t202 + (t329 * t553 - t492 / 0.2e1 - t711 * qJD(1)) * qJD(1) + (-t251 * (-mrSges(4,3) * t476 * t481 - mrSges(4,2) * t477) - t250 * (mrSges(4,1) * t477 - mrSges(4,3) * t590)) * t581 + (t479 * t545 - t341 / 0.2e1) * t181 + (Ifges(6,5) * t291 + Ifges(6,6) * t290 - t512 * t576 + (Ifges(5,3) * t476 + t480 * t513) * qJD(3)) * t650 + (Ifges(5,1) * t342 + Ifges(5,4) * t341) * t656 + (-Ifges(6,1) * t376 - Ifges(6,4) * t375 - Ifges(6,5) * t480) * t680 + (-Ifges(6,4) * t376 - Ifges(6,2) * t375 - Ifges(6,6) * t480) * t681 + t82 * (mrSges(6,1) * t375 - mrSges(6,2) * t376) + (Ifges(5,5) * t342 + Ifges(6,5) * t258 + Ifges(5,6) * t341 + Ifges(6,6) * t257) * t651 + (-t205 + t94) * t466 + (Ifges(7,1) * t155 + Ifges(7,4) * t156) * t672 + (Ifges(7,1) * t171 + Ifges(7,4) * t170) * t673 + (-t615 + t552) * t578 + (-t378 + t539) * t379 + (t351 * t517 + t352 * t520 + t428 * t514) * qJD(3) / 0.2e1 + (-t570 - t288) * t306 - t375 * t689 + (Ifges(7,4) * t284 + Ifges(7,2) * t283 - Ifges(7,6) * t480) * t690 + (Ifges(7,1) * t284 + Ifges(7,4) * t283 - Ifges(7,5) * t480) * t691 + t284 * t692 + t283 * t694 + t290 * t678 + t257 * t679 + t593 * t682 + t155 * t684 + t171 * t685 + t156 * t686 + t170 * t687 - t376 * t688 + t49 * (-mrSges(5,1) * t480 - mrSges(5,3) * t593) + (-Ifges(6,5) * t376 - Ifges(6,6) * t375 + t476 * t513 - t480 * t745) * t661 + t731 * t175 + (t11 * t221 + t12 * t222 + t179 * t722 + t419 * t82 + t59 * t732 + t60 * t731) * m(6) - t480 * t748 - ((-Ifges(3,2) * t559 + t480 * t242 + t476 * t710 + t330 + t439) * t481 + t352 * (Ifges(4,5) * t477 + t481 * t520) + t351 * (Ifges(4,6) * t477 + t481 * t517) + t477 * t240 + t428 * (Ifges(4,3) * t477 + t481 * t514)) * t581 / 0.2e1 - t80 * t594 / 0.2e1 + t48 * (mrSges(5,2) * t480 - mrSges(5,3) * t594) + (-mrSges(5,1) * t716 + mrSges(5,3) * t717) * t145 + (mrSges(5,1) * t718 - mrSges(5,2) * t717) * t219 + (mrSges(5,2) * t716 - mrSges(5,3) * t718) * t146 + t565 + (Ifges(6,4) * t291 + Ifges(6,2) * t290) * t666 + (Ifges(6,4) * t258 + Ifges(6,2) * t257) * t667 + t723 * t226 + t724 * t227 + (t361 * t49 + t362 * t48 + (t219 * t578 + t609) * pkin(9) - t219 * t265 + t723 * t146 + t724 * t145) * m(5) + t725 * t71 + (mrSges(7,2) * t716 + mrSges(7,3) * t726) * t16 + (-mrSges(7,1) * t726 + mrSges(7,2) * t727) * t107 + (-mrSges(7,1) * t716 - mrSges(7,3) * t727) * t15 + (-t250 * t287 - t251 * t288 - pkin(2) * t275 + ((-t250 * t480 - t251 * t476) * qJD(3) + t714) * pkin(9)) * m(4) + t714 * mrSges(4,3) + t743 * t101 + t744 * t100 + (t107 * t725 + t114 * t3 + t115 * t2 + t15 * t743 + t16 * t744 + t311 * t45) * m(7) + (Ifges(6,1) * t291 + Ifges(6,4) * t290) * t664 + (Ifges(6,1) * t258 + Ifges(6,4) * t257) * t665 + (t241 / 0.2e1 + t697) * t536 + t114 * t24 + t115 * t25 + t732 * t176 + t428 * t334 * (mrSges(4,1) * t476 + mrSges(4,2) * t480) + (Ifges(5,4) * t342 + Ifges(5,2) * t341) * t658 - t712 * t480 / 0.2e1 + (mrSges(6,1) * t720 - mrSges(6,2) * t719) * t179 + (t11 * t376 - t12 * t375 + t59 * t719 - t60 * t720) * mrSges(6,3) + t721 * t578 / 0.2e1 + t722 * t127 + (t475 * t545 - t342 / 0.2e1) * t182 + (Ifges(7,4) * t155 + Ifges(7,2) * t156) * t674 + (Ifges(7,4) * t171 + Ifges(7,2) * t170) * t675 + (Ifges(7,5) * t155 + Ifges(7,6) * t156) * t652 + (Ifges(7,5) * t171 + Ifges(7,6) * t170) * t653 + (t540 + t703) * t382 + t480 * t157 / 0.2e1 + t2 * (mrSges(7,2) * t480 + mrSges(7,3) * t283) + t3 * (-mrSges(7,1) * t480 - mrSges(7,3) * t284) + t419 * t47 + t361 * t116 + t362 * t117 + (-Ifges(5,6) * t480 + t476 * t516) * t669 + (-Ifges(5,5) * t480 + t476 * t519) * t670 + t476 * t671 + t291 * t676 + t258 * t677 + (Ifges(7,5) * t284 + Ifges(7,6) * t283 - Ifges(7,3) * t480) * t663 + (Ifges(4,2) * t480 + t622) * t659 + (Ifges(4,1) * t476 + t621) * t660 + (-t518 * t576 + (Ifges(5,5) * t476 + t480 * t519) * qJD(3)) * t655 + (-t515 * t576 + (Ifges(5,6) * t476 + t480 * t516) * qJD(3)) * t657 + (Ifges(4,5) * t476 + Ifges(4,6) * t480) * t645 + t206 * t629 + t523 * t609 - t275 * t526 + (-t462 - t287) * t307 - pkin(2) * t172 + (t399 * t764 + t400 * t699) * g(2) + (t401 * t764 + t402 * t699) * g(1) + t480 * t747 + t221 * t78 + t222 * t77 + t45 * (-mrSges(7,1) * t283 + mrSges(7,2) * t284) - t300 * mrSges(3,2) + t301 * mrSges(3,1) + t311 * t13; (t697 + t770) * t352 + (-t202 + t307) * t251 + (Ifges(7,4) * t169 + Ifges(7,2) * t168) * t675 + (Ifges(7,1) * t169 + Ifges(7,4) * t168) * t673 + (t15 * t169 - t16 * t168 + t2 * t308 - t3 * t309) * mrSges(7,3) + (t552 + t502) * qJD(4) - (Ifges(6,4) * t664 + Ifges(6,2) * t666 + Ifges(6,6) * t650 + t639 + t678) * t395 + (-pkin(3) * t120 - t145 * t177 - t146 * t178 - t219 * t251) * m(5) + (-t617 + t668) * t575 + (t320 * t706 + t321 * t768) * g(2) + (t324 * t706 + t325 * t768) * g(1) + (-mrSges(6,1) * t784 + mrSges(6,2) * t783) * t179 + (t298 * t516 + t299 * t519 + t346 * t513) * qJD(4) / 0.2e1 + (Ifges(7,4) * t309 + Ifges(7,2) * t308) * t690 + (Ifges(7,1) * t309 + Ifges(7,4) * t308) * t691 + t309 * t692 + t308 * t694 - t254 * t677 - t253 * t679 + t475 * t682 + t479 * t683 + t169 * t685 + t168 * t687 + t410 * t688 + (-t502 - t739 / 0.2e1 - t334 * mrSges(4,2) + t516 * t658 + t513 * t651 + t519 * t656 + t615) * t351 + (t397 * t708 + t398 * t707 + t527) * g(3) + t566 + (-Ifges(6,4) * t254 - Ifges(6,2) * t253) * t667 + (-Ifges(6,1) * t254 - Ifges(6,4) * t253) * t665 + (-Ifges(6,5) * t254 - Ifges(6,6) * t253) * t651 + t728 * t71 + (m(5) * ((-t145 * t479 - t146 * t475) * qJD(4) + t713) - t475 * t116 + t479 * t117 - t227 * t575 - t226 * t577) * pkin(10) + (t145 * t605 + t146 * t606 + t713) * mrSges(5,3) + t715 * t127 + (Ifges(7,4) * t672 + Ifges(7,2) * t674 + Ifges(7,6) * t652 + t637 + t686) * t213 + t740 * t101 + t741 * t100 + (t107 * t728 + t15 * t740 + t16 * t741 + t191 * t3 + t192 * t2 + t363 * t45) * m(7) + t733 * t175 + t734 * t176 + (t11 * t332 + t12 * t333 + t179 * t715 - t461 * t82 + t59 * t734 + t60 * t733) * m(6) + (Ifges(7,1) * t672 + Ifges(7,4) * t674 + Ifges(7,5) * t652 - t638 + t684) * t212 - pkin(3) * t94 - (Ifges(4,1) * t351 - t623 + t710) * t352 / 0.2e1 - t577 * t616 - (-Ifges(4,2) * t352 + t345 + t721) * t351 / 0.2e1 - t461 * t47 + t363 * t13 + t515 * t669 + t518 * t670 + (Ifges(7,5) * t309 + Ifges(7,6) * t308) * t663 + t241 * t646 + t595 * t648 + t120 * t524 + ((-t169 + t212) * mrSges(7,2) + (t168 - t213) * mrSges(7,1)) * t107 - t137 * mrSges(4,2) + t138 * mrSges(4,1) + t191 * t24 + t192 * t25 - (Ifges(6,1) * t664 + Ifges(6,4) * t666 + Ifges(6,5) * t650 - t640 + t676) * t396 - t178 * t226 - t177 * t227 + (Ifges(6,1) * t410 - Ifges(6,4) * t510) * t680 + (Ifges(6,4) * t410 - Ifges(6,2) * t510) * t681 + (-t11 * t410 - t12 * t510 + t253 * t60 - t254 * t59) * mrSges(6,3) + (Ifges(6,5) * t410 - Ifges(6,6) * t510 + t512) * t661 + t82 * (mrSges(6,1) * t510 + mrSges(6,2) * t410) - t510 * t689 - t250 * t306 + t45 * (-mrSges(7,1) * t308 + mrSges(7,2) * t309) + (Ifges(7,5) * t169 + Ifges(7,6) * t168) * t653 + t332 * t78 + t333 * t77; t771 + (-mrSges(7,2) * t107 + Ifges(7,1) * t673 + Ifges(7,4) * t675 + Ifges(7,5) * t653 + t638 + t685) * t762 + (-mrSges(6,2) * t179 + Ifges(6,1) * t665 + Ifges(6,4) * t667 + Ifges(6,5) * t651 + t640 + t677) * t542 + (t11 * t472 + t12 * t470) * t693 - t127 * t635 - m(6) * (t179 * t635 + t59 * t64 + t60 * t65) + t742 + (-Ifges(5,2) * t299 + t182 + t296) * t658 + t735 * t100 + t736 * t101 + (-t107 * t162 + t15 * t736 + t16 * t735 + t2 * t388 + t3 * t387) * m(7) + (-m(7) * (-t325 * t421 + t401 * t422) - t588 - t259 * mrSges(6,1) + t260 * mrSges(6,2) + mrSges(5,2) * t268 + t737 * t267) * g(1) + t505 - (mrSges(7,1) * t107 + Ifges(7,4) * t673 + Ifges(7,2) * t675 + Ifges(7,6) * t653 - t637 + t687) * t125 + t387 * t24 + t388 * t25 + (Ifges(5,5) * t298 - Ifges(5,6) * t299) * t651 + t181 * t655 + (Ifges(5,1) * t298 - t613) * t656 + t78 * t632 + t77 * t633 + t299 * t616 + t298 * t617 + (-m(7) * (-t321 * t421 + t399 * t422) - t589 + t756 * mrSges(5,2) + t759 * mrSges(6,1) + t758 * mrSges(6,2) - t737 * t757) * g(2) - t162 * t71 - t65 * t175 - t64 * t176 + (-(-t398 * t463 - t464 * t597) * mrSges(6,1) - (-t398 * t464 + t463 * t597) * mrSges(6,2) - m(6) * t634 - m(7) * (-t398 * t421 - t422 * t597) - t587 - t525) * g(3) - t145 * t226 + t146 * t227 - (mrSges(6,1) * t179 + Ifges(6,4) * t665 + Ifges(6,2) * t667 + Ifges(6,6) * t651 - t639 + t679) * t197 - t219 * (mrSges(5,1) * t299 + mrSges(5,2) * t298); -t762 * t100 + t125 * t101 - t542 * t175 + t197 * t176 + t13 + t47 + (t125 * t15 - t16 * t762 + t45 + t496) * m(7) + (t197 * t59 - t542 * t60 + t496 + t82) * m(6); -t107 * (mrSges(7,1) * t125 + mrSges(7,2) * t762) + (Ifges(7,1) * t762 - t618) * t673 + t67 * t672 + (Ifges(7,5) * t762 - Ifges(7,6) * t125) * t653 - t15 * t100 + t16 * t101 - g(1) * t588 - g(2) * t589 - g(3) * t587 + (t125 * t16 + t15 * t762) * mrSges(7,3) + t505 + (-Ifges(7,2) * t125 + t121 + t68) * t675;];
tau  = t1;
