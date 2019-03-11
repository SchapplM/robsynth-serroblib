% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP9
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:23
% EndTime: 2019-03-10 02:02:27
% DurationCPUTime: 74.13s
% Computational Cost: add. (29970->1124), mult. (71371->1503), div. (0->0), fcn. (57340->14), ass. (0->461)
t443 = sin(qJ(2));
t581 = cos(pkin(6));
t528 = pkin(1) * t581;
t425 = t443 * t528;
t442 = sin(qJ(3));
t446 = cos(qJ(3));
t492 = pkin(3) * t442 - pkin(10) * t446;
t439 = sin(pkin(6));
t447 = cos(qJ(2));
t568 = t439 * t447;
t776 = -(t425 + (pkin(8) + t492) * t568) * qJD(1) + t492 * qJD(3);
t441 = sin(qJ(4));
t445 = cos(qJ(4));
t550 = qJD(1) * t439;
t561 = t446 * t447;
t319 = (-t441 * t561 + t443 * t445) * t550;
t544 = qJD(4) * t445;
t547 = qJD(3) * t446;
t704 = t441 * t547 + t442 * t544 + t319;
t507 = t581 * qJD(1);
t498 = pkin(1) * t507;
t524 = t443 * t550;
t349 = -pkin(8) * t524 + t447 * t498;
t469 = (pkin(2) * t443 - pkin(9) * t447) * t439;
t350 = qJD(1) * t469;
t262 = t446 * t349 + t442 * t350;
t238 = pkin(10) * t524 + t262;
t396 = -pkin(3) * t446 - pkin(10) * t442 - pkin(2);
t546 = qJD(4) * t441;
t548 = qJD(3) * t442;
t707 = -t445 * t238 + t396 * t544 + (-t445 * t548 - t446 * t546) * pkin(9) + t776 * t441;
t536 = pkin(9) * t548;
t775 = t776 * t445 + (t238 + t536) * t441;
t733 = Ifges(6,4) + Ifges(7,4);
t552 = pkin(8) * t568 + t425;
t352 = t552 * qJD(1);
t417 = t507 + qJD(2);
t309 = t417 * pkin(9) + t352;
t321 = (-pkin(2) * t447 - pkin(9) * t443 - pkin(1)) * t550;
t210 = -t442 * t309 + t321 * t446;
t328 = t417 * t446 - t442 * t524;
t499 = t446 * t524;
t329 = t417 * t442 + t499;
t251 = pkin(3) * t329 - pkin(10) * t328;
t141 = -t210 * t441 + t445 * t251;
t448 = -pkin(11) - pkin(10);
t527 = qJD(4) * t448;
t576 = t328 * t445;
t774 = -pkin(4) * t329 + pkin(11) * t576 + t445 * t527 - t141;
t142 = t445 * t210 + t441 * t251;
t577 = t328 * t441;
t773 = -pkin(11) * t577 - t441 * t527 + t142;
t565 = t441 * t443;
t320 = (t445 * t561 + t565) * t550;
t562 = t445 * t446;
t427 = pkin(9) * t562;
t523 = t447 * t550;
t500 = t442 * t523;
t772 = -pkin(4) * t500 + pkin(11) * t320 + (pkin(4) * t442 - pkin(11) * t562) * qJD(3) + (-t427 + (pkin(11) * t442 - t396) * t441) * qJD(4) + t775;
t771 = pkin(11) * t704 - t707;
t598 = Ifges(4,4) * t329;
t400 = qJD(3) - t523;
t718 = t400 * Ifges(4,6);
t203 = t328 * Ifges(4,2) + t598 + t718;
t323 = qJD(4) - t328;
t315 = qJD(5) + t323;
t636 = t315 / 0.2e1;
t268 = -t329 * t441 + t400 * t445;
t269 = t329 * t445 + t400 * t441;
t440 = sin(qJ(5));
t444 = cos(qJ(5));
t164 = t268 * t440 + t269 * t444;
t658 = t164 / 0.2e1;
t506 = t444 * t268 - t269 * t440;
t661 = t506 / 0.2e1;
t729 = Ifges(7,3) + Ifges(6,3);
t730 = Ifges(6,6) + Ifges(7,6);
t732 = Ifges(6,5) + Ifges(7,5);
t690 = t269 * Ifges(5,5) + t268 * Ifges(5,6) + t323 * Ifges(5,3) + t164 * t732 + t315 * t729 + t506 * t730;
t770 = t636 * t729 + t658 * t732 + t661 * t730 - t203 / 0.2e1 + t690 / 0.2e1;
t438 = qJ(4) + qJ(5);
t434 = cos(t438);
t610 = pkin(4) * t445;
t394 = pkin(5) * t434 + t610;
t386 = pkin(3) + t394;
t431 = pkin(3) + t610;
t488 = -mrSges(5,1) * t445 + mrSges(5,2) * t441;
t769 = -m(5) * pkin(3) - m(6) * t431 - m(7) * t386 + t488;
t688 = m(6) * t448 + m(7) * (-qJ(6) + t448) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(10) - mrSges(5,3);
t768 = mrSges(6,1) + mrSges(7,1);
t767 = mrSges(6,2) + mrSges(7,2);
t734 = Ifges(6,1) + Ifges(7,1);
t731 = Ifges(6,2) + Ifges(7,2);
t382 = t440 * t445 + t441 * t444;
t693 = qJD(4) + qJD(5);
t289 = t693 * t382;
t474 = t440 * t441 - t444 * t445;
t195 = -t289 * t442 - t474 * t547;
t226 = t319 * t440 + t320 * t444;
t710 = t195 - t226;
t356 = t474 * t442;
t196 = t356 * t693 - t382 * t547;
t225 = t319 * t444 - t320 * t440;
t709 = t196 - t225;
t223 = t382 * t328;
t766 = t223 - t289;
t224 = t474 * t328;
t288 = t693 * t474;
t765 = t224 - t288;
t261 = -t442 * t349 + t350 * t446;
t237 = -pkin(3) * t524 - t261;
t432 = pkin(9) * t547;
t764 = t432 - t237;
t701 = t500 - t548;
t182 = -pkin(3) * t400 - t210;
t143 = -pkin(4) * t268 + t182;
t308 = -t417 * pkin(2) - t349;
t178 = -t328 * pkin(3) - t329 * pkin(10) + t308;
t211 = t446 * t309 + t442 * t321;
t183 = pkin(10) * t400 + t211;
t115 = t445 * t178 - t183 * t441;
t91 = -pkin(11) * t269 + t115;
t78 = pkin(4) * t323 + t91;
t116 = t178 * t441 + t183 * t445;
t92 = pkin(11) * t268 + t116;
t88 = t440 * t92;
t35 = t444 * t78 - t88;
t753 = qJ(6) * t164;
t28 = t35 - t753;
t25 = pkin(5) * t315 + t28;
t87 = -pkin(5) * t506 + qJD(6) + t143;
t682 = t87 * mrSges(7,2) - mrSges(6,3) * t35 - mrSges(7,3) * t25;
t763 = t143 * mrSges(6,2) + t682;
t90 = t444 * t92;
t36 = t440 * t78 + t90;
t712 = qJ(6) * t506;
t29 = t36 + t712;
t683 = -t87 * mrSges(7,1) + mrSges(6,3) * t36 + mrSges(7,3) * t29;
t748 = t733 * t164;
t721 = t315 * t730 + t506 * t731 + t748;
t758 = t721 / 0.2e1;
t762 = -t143 * mrSges(6,1) + t683 + t758;
t587 = t211 * mrSges(4,3);
t761 = t587 - t308 * mrSges(4,1) - t115 * mrSges(5,1) - t35 * mrSges(6,1) - t25 * mrSges(7,1) + t116 * mrSges(5,2) + t36 * mrSges(6,2) + t29 * mrSges(7,2) + t718 / 0.2e1;
t752 = t733 * t506;
t720 = t164 * t734 + t315 * t732 + t752;
t759 = t720 / 0.2e1;
t541 = qJD(1) * qJD(2);
t358 = (qJDD(1) * t443 + t447 * t541) * t439;
t505 = t581 * qJDD(1);
t416 = t505 + qJDD(2);
t222 = -qJD(3) * t499 - t442 * t358 + t416 * t446 - t417 * t548;
t209 = qJDD(4) - t222;
t198 = qJDD(5) + t209;
t221 = qJD(3) * t328 + t358 * t446 + t416 * t442;
t357 = (-qJDD(1) * t447 + t443 * t541) * t439;
t346 = qJDD(3) + t357;
t132 = qJD(4) * t268 + t221 * t445 + t346 * t441;
t133 = -qJD(4) * t269 - t221 * t441 + t346 * t445;
t55 = qJD(5) * t506 + t132 * t444 + t133 * t440;
t56 = -qJD(5) * t164 - t132 * t440 + t133 * t444;
t726 = t198 * t732 + t55 * t734 + t56 * t733;
t757 = t726 / 0.2e1;
t727 = t198 * t730 + t55 * t733 + t56 * t731;
t756 = -t727 / 0.2e1;
t380 = t445 * t396;
t564 = t442 * t445;
t284 = -pkin(11) * t564 + t380 + (-pkin(9) * t441 - pkin(4)) * t446;
t339 = t441 * t396 + t427;
t566 = t441 * t442;
t304 = -pkin(11) * t566 + t339;
t185 = t440 * t284 + t444 * t304;
t723 = -qJD(5) * t185 + t440 * t771 + t444 * t772;
t542 = qJD(5) * t444;
t543 = qJD(5) * t440;
t722 = t284 * t542 - t304 * t543 + t440 * t772 - t444 * t771;
t754 = mrSges(4,2) + t688;
t402 = t448 * t441;
t403 = t448 * t445;
t312 = t440 * t402 - t444 * t403;
t714 = -qJD(5) * t312 + t440 * t773 + t444 * t774;
t713 = t402 * t542 + t403 * t543 + t440 * t774 - t444 * t773;
t702 = pkin(4) * t704 + t764;
t540 = qJDD(1) * t439;
t751 = pkin(8) * t540 + qJD(2) * t498;
t433 = sin(t438);
t612 = pkin(4) * t441;
t393 = pkin(5) * t433 + t612;
t487 = t441 * mrSges(5,1) + t445 * mrSges(5,2);
t750 = -m(7) * t393 - mrSges(4,3) - t487;
t749 = -pkin(8) * t439 * t541 + pkin(1) * t505;
t700 = -t211 + (t546 - t577) * pkin(4);
t622 = cos(qJ(1));
t494 = t581 * t622;
t621 = sin(qJ(1));
t366 = t443 * t494 + t447 * t621;
t526 = t439 * t622;
t297 = t366 * t446 - t442 * t526;
t365 = t443 * t621 - t447 * t494;
t747 = t297 * t433 - t365 * t434;
t228 = -t297 * t434 - t365 * t433;
t746 = t297 * t441 - t365 * t445;
t745 = t297 * t445 + t365 * t441;
t490 = t446 * mrSges(4,1) - t442 * mrSges(4,2);
t743 = t442 * t688 + t446 * t769 - t490;
t271 = -t443 * t751 + t447 * t749;
t247 = -t416 * pkin(2) - t271;
t109 = -t222 * pkin(3) - t221 * pkin(10) + t247;
t270 = t443 * t749 + t447 * t751;
t246 = pkin(9) * t416 + t270;
t534 = pkin(1) * t540;
t254 = pkin(2) * t357 - pkin(9) * t358 - t534;
t111 = t446 * t246 + t442 * t254 - t309 * t548 + t321 * t547;
t96 = pkin(10) * t346 + t111;
t26 = t441 * t109 + t178 * t544 - t183 * t546 + t445 * t96;
t27 = -qJD(4) * t116 + t445 * t109 - t441 * t96;
t742 = t27 * mrSges(5,1) - t26 * mrSges(5,2);
t741 = t636 * t730 + t658 * t733 + t661 * t731;
t740 = t636 * t732 + t658 * t734 + t661 * t733;
t21 = pkin(4) * t209 - pkin(11) * t132 + t27;
t24 = pkin(11) * t133 + t26;
t6 = -qJD(5) * t36 + t444 * t21 - t24 * t440;
t2 = pkin(5) * t198 - qJ(6) * t55 - qJD(6) * t164 + t6;
t5 = t440 * t21 + t444 * t24 + t78 * t542 - t543 * t92;
t3 = qJ(6) * t56 + qJD(6) * t506 + t5;
t739 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t666 = t132 / 0.2e1;
t665 = t133 / 0.2e1;
t651 = t209 / 0.2e1;
t650 = t221 / 0.2e1;
t649 = t222 / 0.2e1;
t629 = t346 / 0.2e1;
t662 = -t506 / 0.2e1;
t728 = t198 * t729 + t55 * t732 + t56 * t730;
t725 = -pkin(5) * t701 - qJ(6) * t710 + qJD(6) * t356 + t723;
t355 = t382 * t442;
t724 = qJ(6) * t709 - qJD(6) * t355 + t722;
t719 = t400 * Ifges(4,5);
t675 = m(6) * pkin(4);
t717 = -mrSges(5,1) - t675;
t716 = -pkin(5) * t329 - qJ(6) * t765 - qJD(6) * t382 + t714;
t715 = qJ(6) * t766 - qJD(6) * t474 + t713;
t569 = t439 * t443;
t375 = -pkin(8) * t569 + t447 * t528;
t341 = -pkin(2) * t581 - t375;
t363 = t442 * t569 - t446 * t581;
t364 = t442 * t581 + t446 * t569;
t234 = t363 * pkin(3) - t364 * pkin(10) + t341;
t342 = pkin(9) * t581 + t552;
t553 = pkin(2) * t568 + pkin(9) * t569;
t616 = pkin(1) * t439;
t343 = -t553 - t616;
t253 = t446 * t342 + t442 * t343;
t236 = -pkin(10) * t568 + t253;
t135 = t445 * t234 - t236 * t441;
t468 = -t364 * t445 + t441 * t568;
t108 = pkin(4) * t363 + pkin(11) * t468 + t135;
t136 = t441 * t234 + t445 * t236;
t294 = -t364 * t441 - t445 * t568;
t117 = pkin(11) * t294 + t136;
t60 = t440 * t108 + t444 * t117;
t711 = -pkin(5) * t709 + t702;
t264 = Ifges(5,4) * t268;
t146 = t269 * Ifges(5,1) + t323 * Ifges(5,5) + t264;
t322 = Ifges(4,4) * t328;
t204 = t329 * Ifges(4,1) + t322 + t719;
t708 = t445 * t146 + t204;
t706 = -qJD(4) * t339 + t775;
t705 = -pkin(5) * t766 + t700;
t545 = qJD(4) * t442;
t703 = t441 * t545 - t445 * t547 + t320;
t282 = -t364 * t433 - t434 * t568;
t699 = -t767 * (-t364 * t434 + t433 * t568) - t768 * t282;
t493 = t581 * t621;
t368 = -t443 * t493 + t447 * t622;
t525 = t439 * t621;
t301 = t368 * t446 + t442 * t525;
t367 = t443 * t622 + t447 * t493;
t231 = -t301 * t433 + t367 * t434;
t232 = t301 * t434 + t367 * t433;
t698 = -t231 * t768 + t767 * t232;
t697 = -t767 * t228 + t747 * t768;
t112 = -t442 * t246 + t254 * t446 - t309 * t547 - t321 * t548;
t696 = t111 * t446 - t112 * t442;
t695 = t26 * t445 - t27 * t441;
t694 = -m(6) - m(7) - m(4);
t600 = Ifges(3,4) * t443;
t677 = t439 ^ 2;
t692 = (t443 * (Ifges(3,1) * t447 - t600) / 0.2e1 - pkin(1) * (mrSges(3,1) * t443 + mrSges(3,2) * t447)) * t677;
t61 = Ifges(5,5) * t132 + Ifges(5,6) * t133 + Ifges(5,3) * t209;
t691 = t61 + t728;
t689 = -t767 * t433 + t434 * t768 - t769;
t686 = mrSges(4,1) + t689;
t504 = mrSges(3,3) * t524;
t684 = -m(4) * t308 + mrSges(3,1) * t417 + mrSges(4,1) * t328 - mrSges(4,2) * t329 - t504;
t681 = mrSges(3,1) - t743;
t680 = m(6) * (pkin(9) + t612) + m(7) * (pkin(9) + t393) - mrSges(3,2) + mrSges(4,3);
t679 = -m(5) * pkin(9) - m(6) * t612 + mrSges(3,2) + t750;
t674 = m(7) * pkin(5);
t673 = t55 / 0.2e1;
t672 = t56 / 0.2e1;
t62 = t132 * Ifges(5,4) + t133 * Ifges(5,2) + t209 * Ifges(5,6);
t671 = t62 / 0.2e1;
t670 = Ifges(5,1) * t666 + Ifges(5,4) * t665 + Ifges(5,5) * t651;
t667 = Ifges(4,1) * t650 + Ifges(4,4) * t649 + Ifges(4,5) * t629;
t664 = t146 / 0.2e1;
t659 = -t164 / 0.2e1;
t653 = t198 / 0.2e1;
t644 = -t268 / 0.2e1;
t643 = t268 / 0.2e1;
t642 = -t269 / 0.2e1;
t641 = t269 / 0.2e1;
t637 = -t315 / 0.2e1;
t635 = -t323 / 0.2e1;
t634 = t323 / 0.2e1;
t632 = t328 / 0.2e1;
t630 = t329 / 0.2e1;
t615 = pkin(4) * t269;
t614 = pkin(4) * t294;
t611 = pkin(4) * t444;
t608 = pkin(9) * t446;
t435 = t442 * pkin(9);
t45 = t444 * t91 - t88;
t604 = mrSges(6,3) * t506;
t603 = mrSges(6,3) * t164;
t602 = mrSges(7,3) * t506;
t601 = mrSges(7,3) * t164;
t599 = Ifges(3,4) * t447;
t597 = Ifges(4,4) * t442;
t596 = Ifges(4,4) * t446;
t595 = Ifges(5,4) * t269;
t594 = Ifges(5,4) * t441;
t593 = Ifges(5,4) * t445;
t590 = t115 * mrSges(5,3);
t589 = t116 * mrSges(5,3);
t588 = t210 * mrSges(4,3);
t97 = -pkin(3) * t346 - t112;
t582 = t442 * t97;
t571 = t433 * t446;
t570 = t434 * t446;
t145 = t268 * Ifges(5,2) + t323 * Ifges(5,6) + t595;
t567 = t441 * t145;
t391 = pkin(4) * t566 + t435;
t551 = t622 * pkin(1) + pkin(8) * t525;
t549 = qJD(2) * t447;
t538 = pkin(4) * t543;
t537 = pkin(4) * t542;
t532 = Ifges(4,5) * t221 + Ifges(4,6) * t222 + Ifges(4,3) * t346;
t531 = Ifges(3,5) * t358 - Ifges(3,6) * t357 + Ifges(3,3) * t416;
t530 = t368 * pkin(2) + t551;
t522 = qJD(2) * t569;
t521 = t439 * t549;
t518 = t569 / 0.2e1;
t517 = -t567 / 0.2e1;
t19 = -t56 * mrSges(7,1) + t55 * mrSges(7,2);
t510 = -t545 / 0.2e1;
t44 = -t440 * t91 - t90;
t59 = t444 * t108 - t117 * t440;
t184 = t444 * t284 - t304 * t440;
t252 = -t442 * t342 + t343 * t446;
t311 = t444 * t402 + t403 * t440;
t503 = mrSges(3,3) * t523;
t495 = -pkin(1) * t621 + pkin(8) * t526;
t235 = pkin(3) * t568 - t252;
t491 = mrSges(4,1) * t363 + mrSges(4,2) * t364;
t489 = mrSges(5,1) * t294 + mrSges(5,2) * t468;
t486 = Ifges(4,1) * t446 - t597;
t485 = Ifges(5,1) * t445 - t594;
t484 = Ifges(5,1) * t441 + t593;
t483 = -Ifges(4,2) * t442 + t596;
t482 = -Ifges(5,2) * t441 + t593;
t481 = Ifges(5,2) * t445 + t594;
t480 = Ifges(4,5) * t446 - Ifges(4,6) * t442;
t479 = Ifges(5,5) * t445 - Ifges(5,6) * t441;
t478 = Ifges(5,5) * t441 + Ifges(5,6) * t445;
t189 = t294 * t444 + t440 * t468;
t190 = t294 * t440 - t444 * t468;
t239 = -t301 * t441 + t367 * t445;
t472 = t367 * pkin(9) + t530;
t470 = -t366 * pkin(2) + t495;
t351 = qJD(2) * t469;
t353 = t375 * qJD(2);
t155 = -t342 * t547 - t343 * t548 + t351 * t446 - t442 * t353;
t467 = t182 * t487;
t465 = (t447 * Ifges(3,2) + t600) * t439;
t293 = -qJD(3) * t363 + t446 * t521;
t177 = qJD(4) * t294 + t293 * t445 + t441 * t522;
t292 = qJD(3) * t364 + t442 * t521;
t154 = -t342 * t548 + t343 * t547 + t442 * t351 + t446 * t353;
t150 = pkin(10) * t522 + t154;
t354 = t552 * qJD(2);
t172 = t292 * pkin(3) - t293 * pkin(10) + t354;
t68 = -qJD(4) * t136 - t150 * t441 + t445 * t172;
t34 = pkin(4) * t292 - pkin(11) * t177 + t68;
t176 = qJD(4) * t468 - t293 * t441 + t445 * t522;
t67 = t445 * t150 + t441 * t172 + t234 * t544 - t236 * t546;
t48 = pkin(11) * t176 + t67;
t9 = t108 * t542 - t117 * t543 + t440 * t34 + t444 * t48;
t165 = t235 - t614;
t296 = t366 * t442 + t446 * t526;
t458 = t417 * t439 * (Ifges(3,5) * t447 - Ifges(3,6) * t443);
t454 = -t365 * pkin(9) + t470;
t151 = -pkin(3) * t522 - t155;
t66 = -pkin(4) * t133 + t97;
t10 = -qJD(5) * t60 + t444 * t34 - t440 * t48;
t450 = t728 + t739;
t98 = -pkin(4) * t176 + t151;
t430 = pkin(5) + t611;
t411 = Ifges(3,4) * t523;
t369 = (-mrSges(3,1) * t447 + mrSges(3,2) * t443) * t439;
t361 = t367 * pkin(2);
t359 = t365 * pkin(2);
t348 = -mrSges(3,2) * t417 + t503;
t340 = pkin(5) * t474 - t431;
t338 = -t441 * t608 + t380;
t306 = Ifges(3,1) * t524 + t417 * Ifges(3,5) + t411;
t305 = t417 * Ifges(3,6) + qJD(1) * t465;
t300 = t368 * t442 - t446 * t525;
t286 = pkin(5) * t355 + t391;
t275 = mrSges(4,1) * t400 - mrSges(4,3) * t329;
t274 = -mrSges(4,2) * t400 + mrSges(4,3) * t328;
t267 = -qJ(6) * t474 + t312;
t266 = -qJ(6) * t382 + t311;
t240 = t301 * t445 + t367 * t441;
t202 = t329 * Ifges(4,5) + t328 * Ifges(4,6) + t400 * Ifges(4,3);
t188 = mrSges(5,1) * t323 - mrSges(5,3) * t269;
t187 = -mrSges(5,2) * t323 + mrSges(5,3) * t268;
t171 = -mrSges(4,2) * t346 + mrSges(4,3) * t222;
t170 = mrSges(4,1) * t346 - mrSges(4,3) * t221;
t169 = -mrSges(5,1) * t268 + mrSges(5,2) * t269;
t157 = -qJ(6) * t355 + t185;
t156 = -pkin(5) * t446 + qJ(6) * t356 + t184;
t140 = mrSges(6,1) * t315 - t603;
t139 = mrSges(7,1) * t315 - t601;
t138 = -mrSges(6,2) * t315 + t604;
t137 = -mrSges(7,2) * t315 + t602;
t134 = -mrSges(4,1) * t222 + mrSges(4,2) * t221;
t128 = pkin(5) * t164 + t615;
t123 = t221 * Ifges(4,4) + t222 * Ifges(4,2) + t346 * Ifges(4,6);
t113 = -pkin(5) * t189 + t165;
t100 = -mrSges(6,1) * t506 + mrSges(6,2) * t164;
t99 = -mrSges(7,1) * t506 + mrSges(7,2) * t164;
t94 = -mrSges(5,2) * t209 + mrSges(5,3) * t133;
t93 = mrSges(5,1) * t209 - mrSges(5,3) * t132;
t75 = -qJD(5) * t190 + t176 * t444 - t177 * t440;
t74 = qJD(5) * t189 + t176 * t440 + t177 * t444;
t73 = -mrSges(5,1) * t133 + mrSges(5,2) * t132;
t46 = qJ(6) * t189 + t60;
t42 = -pkin(5) * t75 + t98;
t41 = pkin(5) * t363 - qJ(6) * t190 + t59;
t40 = -mrSges(6,2) * t198 + mrSges(6,3) * t56;
t39 = -mrSges(7,2) * t198 + mrSges(7,3) * t56;
t38 = mrSges(6,1) * t198 - mrSges(6,3) * t55;
t37 = mrSges(7,1) * t198 - mrSges(7,3) * t55;
t31 = t45 - t753;
t30 = t44 - t712;
t22 = -pkin(5) * t56 + qJDD(6) + t66;
t20 = -mrSges(6,1) * t56 + mrSges(6,2) * t55;
t8 = qJ(6) * t75 + qJD(6) * t189 + t9;
t7 = pkin(5) * t292 - qJ(6) * t74 - qJD(6) * t190 + t10;
t1 = [(Ifges(4,1) * t364 - Ifges(4,5) * t568) * t650 + (Ifges(4,1) * t293 + Ifges(4,5) * t522) * t630 + (Ifges(5,5) * t177 + Ifges(5,6) * t176) * t634 + (-t111 * mrSges(4,3) - Ifges(4,2) * t649 - Ifges(4,4) * t650 - Ifges(4,6) * t629 + Ifges(5,6) * t665 + Ifges(5,5) * t666 + Ifges(5,3) * t651 - t123 / 0.2e1 + t730 * t672 + t732 * t673 + t729 * t653 + t691 / 0.2e1 + t739 + t742) * t363 + (Ifges(5,1) * t177 + Ifges(5,4) * t176) * t641 + m(4) * (t111 * t253 + t112 * t252 + t154 * t211 + t155 * t210 + t247 * t341) + (t439 * t306 + t677 * qJD(1) * (-Ifges(3,2) * t443 + t599)) * t549 / 0.2e1 + (-t270 * t581 - t358 * t616 - t416 * t552) * mrSges(3,2) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t677 + t270 * t552 + t271 * t375 + t352 * t353) + (t111 * t568 - t211 * t522 + t293 * t308) * mrSges(4,2) + (t202 * t518 + t458 / 0.2e1) * qJD(2) + (Ifges(4,4) * t364 - Ifges(4,6) * t568) * t649 + (Ifges(4,4) * t293 + Ifges(4,6) * t522) * t632 - t357 * (Ifges(3,6) * t581 + t465) / 0.2e1 + t416 * (Ifges(3,3) * t581 + (Ifges(3,5) * t443 + Ifges(3,6) * t447) * t439) / 0.2e1 + t581 * t531 / 0.2e1 - t532 * t568 / 0.2e1 + t112 * (-mrSges(4,1) * t568 - mrSges(4,3) * t364) + t294 * t671 - t305 * t522 / 0.2e1 + t210 * (mrSges(4,1) * t522 - mrSges(4,3) * t293) + (-Ifges(5,4) * t468 + Ifges(5,2) * t294) * t665 + (t741 + t762) * t75 + (t759 + t740 + t763) * t74 + m(7) * (t113 * t22 + t2 * t41 + t25 * t7 + t29 * t8 + t3 * t46 + t42 * t87) + m(6) * (t10 * t35 + t143 * t98 + t165 * t66 + t36 * t9 + t5 * t60 + t59 * t6) + m(5) * (t115 * t68 + t116 * t67 + t135 * t27 + t136 * t26 + t151 * t182 + t235 * t97) + (Ifges(5,4) * t177 + Ifges(5,2) * t176) * t643 + (-m(6) * (-t297 * t431 + t470) - m(3) * t495 + t366 * mrSges(3,1) - mrSges(3,3) * t526 + mrSges(2,1) * t621 + mrSges(2,2) * t622 - m(5) * (-pkin(3) * t297 + t454) + t745 * mrSges(5,1) - t746 * mrSges(5,2) - m(7) * (-t297 * t386 + t470) - m(4) * t454 + t297 * mrSges(4,1) + t680 * t365 - t768 * t228 - t767 * t747 - t754 * t296) * g(1) + (-m(6) * (t301 * t431 + t530) - m(4) * t472 - t301 * mrSges(4,1) - mrSges(2,1) * t622 + mrSges(2,2) * t621 - m(5) * (pkin(3) * t301 + t472) - t240 * mrSges(5,1) - t239 * mrSges(5,2) - m(7) * (t301 * t386 + t530) - m(3) * t551 - t368 * mrSges(3,1) - mrSges(3,3) * t525 - t680 * t367 - t768 * t232 - t767 * t231 + t754 * t300) * g(2) + t247 * t491 + t177 * t664 + t364 * t667 + t358 * (Ifges(3,5) * t581 + (t443 * Ifges(3,1) + t599) * t439) / 0.2e1 - t369 * t534 - t97 * t489 + t400 * (Ifges(4,5) * t293 + Ifges(4,3) * t522) / 0.2e1 + (Ifges(4,5) * t364 - Ifges(4,3) * t568) * t629 + (-Ifges(5,5) * t468 + Ifges(5,6) * t294) * t651 + (-Ifges(5,1) * t468 + Ifges(5,4) * t294) * t666 + (-t115 * t177 + t116 * t176 + t26 * t294 + t27 * t468) * mrSges(5,3) - t468 * t670 + t60 * t40 + t59 * t38 + t46 * t39 + t41 * t37 + Ifges(2,3) * qJDD(1) + t353 * t348 + t341 * t134 + t293 * t204 / 0.2e1 + t154 * t274 + t155 * t275 + t252 * t170 + t253 * t171 + t235 * t73 + (-Ifges(4,4) * t630 + Ifges(5,5) * t641 - Ifges(4,2) * t632 + Ifges(5,6) * t643 + Ifges(5,3) * t634 - t761 + t770) * t292 + (t270 * t568 - t271 * t569 - t349 * t521 - t352 * t522 - t357 * t552 - t358 * t375) * mrSges(3,3) + (t271 * t581 - t357 * t616 + t375 * t416) * mrSges(3,1) + (Ifges(3,4) * t358 - Ifges(3,2) * t357 + Ifges(3,6) * t416) * t568 / 0.2e1 + (Ifges(3,1) * t358 - Ifges(3,4) * t357 + Ifges(3,5) * t416) * t518 + t42 * t99 + t98 * t100 + t113 * t19 + t135 * t93 + t136 * t94 + (-m(3) * t349 - t684) * t354 + t8 * t137 + t9 * t138 + t7 * t139 + t10 * t140 + (t5 * mrSges(6,3) + t3 * mrSges(7,3) + t731 * t672 + t733 * t673 + t730 * t653 + t727 / 0.2e1 - t66 * mrSges(6,1) - t22 * mrSges(7,1)) * t189 + (t66 * mrSges(6,2) + t22 * mrSges(7,2) - t6 * mrSges(6,3) - t2 * mrSges(7,3) + t732 * t653 + t733 * t672 + t734 * t673 + t757) * t190 + t165 * t20 + t151 * t169 + t176 * t145 / 0.2e1 + t182 * (-mrSges(5,1) * t176 + mrSges(5,2) * t177) + t692 * t541 + t67 * t187 + t68 * t188; t764 * t169 + (t328 * t483 + t329 * t486 + t400 * t480) * qJD(3) / 0.2e1 + (t445 * t510 - t319 / 0.2e1) * t145 + (-t432 - t261) * t275 + (t503 - t348) * t349 + t400 * t308 * (mrSges(4,1) * t442 + mrSges(4,2) * t446) + t203 * t500 / 0.2e1 + t27 * (-mrSges(5,1) * t446 - mrSges(5,3) * t564) + (-t170 + t73) * t435 - ((-Ifges(3,2) * t524 + t446 * t204 + t442 * t690 + t306 + t411) * t447 + t400 * (Ifges(4,3) * t443 + t447 * t480) + t329 * (Ifges(4,5) * t443 + t447 * t486) + t328 * (Ifges(4,6) * t443 + t447 * t483) + t443 * t202) * t550 / 0.2e1 + (-t536 - t262) * t274 + t564 * t670 + (-t210 * (mrSges(4,1) * t443 - mrSges(4,3) * t561) - t211 * (-mrSges(4,3) * t442 * t447 - mrSges(4,2) * t443)) * t550 + (t441 * t510 - t320 / 0.2e1) * t146 + (t369 + (-m(5) + t694) * t553 + (-t768 * (t433 * t443 + t434 * t561) - t767 * (-t433 * t561 + t434 * t443) - t565 * t675 + t750 * t443 + t743 * t447) * t439) * g(3) + (m(5) * t361 - t768 * (-t367 * t570 + t368 * t433) - t767 * (t367 * t571 + t368 * t434) + t694 * (t368 * pkin(9) - t361) + t679 * t368 + t681 * t367) * g(1) + (m(5) * t359 - t768 * (-t365 * t570 + t366 * t433) - t767 * (t365 * t571 + t366 * t434) + t694 * (t366 * pkin(9) - t359) + t679 * t366 + t681 * t365) * g(2) + (-Ifges(5,6) * t446 + t442 * t482) * t665 + (-Ifges(5,5) * t446 + t442 * t485) * t666 + t442 * t667 + (Ifges(4,2) * t446 + t597) * t649 + (Ifges(4,1) * t442 + t596) * t650 + (-Ifges(5,3) * t446 + t442 * t479) * t651 + (-t484 * t545 + (Ifges(5,5) * t442 + t446 * t485) * qJD(3)) * t641 + (Ifges(5,1) * t320 + Ifges(5,4) * t319 + Ifges(5,5) * t500) * t642 + (-t481 * t545 + (Ifges(5,6) * t442 + t446 * t482) * qJD(3)) * t643 + (Ifges(5,4) * t320 + Ifges(5,2) * t319 + Ifges(5,6) * t500) * t644 + (Ifges(4,5) * t442 + Ifges(4,6) * t446) * t629 + (-t478 * t545 + (Ifges(5,3) * t442 + t446 * t479) * qJD(3)) * t634 + (Ifges(5,5) * t320 + Ifges(5,6) * t319 + Ifges(5,3) * t500) * t635 + t531 - t247 * t490 + (t305 * t518 - t458 / 0.2e1 - t692 * qJD(1)) * qJD(1) + t446 * t123 / 0.2e1 + t391 * t20 + t338 * t93 + t339 * t94 + (t195 * t733 + t196 * t731) * t661 + (-t355 * t731 - t356 * t733 - t446 * t730) * t672 + (t225 * t731 + t226 * t733 + t500 * t730) * t662 + (t195 * t734 + t196 * t733) * t658 + (-t355 * t733 - t356 * t734 - t446 * t732) * t673 + (t225 * t733 + t226 * t734 + t500 * t732) * t659 + (t195 * t732 + t196 * t730) * t636 + (t225 * t730 + t226 * t732 + t500 * t729) * t637 + (-t355 * t730 - t356 * t732 - t446 * t729) * t653 - t720 * t226 / 0.2e1 - t721 * t225 / 0.2e1 + t722 * t138 + t723 * t140 + (t143 * t702 + t184 * t6 + t185 * t5 + t35 * t723 + t36 * t722 + t391 * t66) * m(6) + t724 * t137 + t725 * t139 + (t156 * t2 + t157 * t3 + t22 * t286 + t25 * t725 + t29 * t724 + t711 * t87) * m(7) - t726 * t356 / 0.2e1 + t286 * t19 + t271 * mrSges(3,1) - t270 * mrSges(3,2) + (-t588 + t517) * t547 + (-mrSges(6,1) * t701 - mrSges(6,3) * t710) * t35 + (-mrSges(7,1) * t701 - mrSges(7,3) * t710) * t25 + t711 * t99 + t708 * t547 / 0.2e1 + (mrSges(6,2) * t701 + mrSges(6,3) * t709) * t36 + (mrSges(7,2) * t701 + mrSges(7,3) * t709) * t29 + (-mrSges(7,1) * t709 + mrSges(7,2) * t710) * t87 + (-mrSges(6,1) * t709 + mrSges(6,2) * t710) * t143 + (-t587 + t770) * t548 + t5 * (mrSges(6,2) * t446 - mrSges(6,3) * t355) + t3 * (mrSges(7,2) * t446 - mrSges(7,3) * t355) + t6 * (-mrSges(6,1) * t446 + mrSges(6,3) * t356) + t2 * (-mrSges(7,1) * t446 + mrSges(7,3) * t356) + t66 * (mrSges(6,1) * t355 - mrSges(6,2) * t356) + t22 * (mrSges(7,1) * t355 - mrSges(7,2) * t356) - pkin(2) * t134 + (t504 + t684) * t352 + t355 * t756 + t196 * t758 + t195 * t759 + t156 * t37 + t157 * t39 - t62 * t566 / 0.2e1 + t26 * (mrSges(5,2) * t446 - mrSges(5,3) * t566) + t184 * t38 + t185 * t40 - t691 * t446 / 0.2e1 + (-pkin(2) * t247 + ((-t210 * t446 - t211 * t442) * qJD(3) + t696) * pkin(9) - t210 * t261 - t211 * t262) * m(4) + t696 * mrSges(4,3) + t702 * t100 + (-mrSges(5,1) * t701 + mrSges(5,3) * t703) * t115 + (mrSges(5,1) * t704 - mrSges(5,2) * t703) * t182 + (mrSges(5,2) * t701 - mrSges(5,3) * t704) * t116 + t706 * t188 + t707 * t187 + (t26 * t339 + t27 * t338 + (t182 * t547 + t582) * pkin(9) - t182 * t237 + t707 * t116 + t706 * t115) * m(5) + t487 * t582 + t171 * t608; (-mrSges(6,1) * t766 + mrSges(6,2) * t765) * t143 + (Ifges(5,5) * t642 + Ifges(5,6) * t644 + Ifges(5,3) * t635 + t729 * t637 + t732 * t659 + t730 * t662 + t761) * t329 + (-pkin(3) * t97 - t115 * t141 - t116 * t142 - t182 * t211) * m(5) + (t268 * t482 + t269 * t485 + t323 * t479) * qJD(4) / 0.2e1 + (-t590 + t664) * t544 - (t682 + t740) * t288 - (t683 + t741) * t289 + (t275 - t169) * t211 - pkin(3) * t73 - t546 * t589 + t441 * t670 + t445 * t671 + t532 + (t382 * t733 - t474 * t731) * t672 + (t382 * t734 - t474 * t733) * t673 + (t382 * t732 - t474 * t730) * t653 + (t223 * t36 - t224 * t35 - t382 * t6 - t474 * t5) * mrSges(6,3) + (-t2 * t382 + t223 * t29 - t224 * t25 - t3 * t474) * mrSges(7,3) + t66 * (mrSges(6,1) * t474 + mrSges(6,2) * t382) + t22 * (mrSges(7,1) * t474 + mrSges(7,2) * t382) + (-t288 / 0.2e1 + t224 / 0.2e1) * t720 + (-t289 / 0.2e1 + t223 / 0.2e1) * t721 + t481 * t665 + t484 * t666 + t478 * t651 + t203 * t630 + t567 * t632 + t97 * t488 + (t296 * t686 + t297 * t754) * g(2) + (t300 * t686 + t301 * t754) * g(1) + (t517 + t467) * qJD(4) - t87 * (mrSges(7,1) * t223 - mrSges(7,2) * t224) - t431 * t20 + t340 * t19 + (-t223 * t731 - t224 * t733) * t662 + t311 * t38 + t312 * t40 + (-t223 * t733 - t224 * t734) * t659 + (-t223 * t730 - t224 * t732) * t637 - t210 * t274 + (t485 * t642 + t482 * t644 + t479 * t635 - t467 - t719 / 0.2e1 - t308 * mrSges(4,2) + t588) * t328 + t713 * t138 + t266 * t37 + t267 * t39 + (t143 * t700 + t311 * t6 + t312 * t5 + t35 * t714 + t36 * t713 - t431 * t66) * m(6) + t714 * t140 + t715 * t137 + t716 * t139 + (t2 * t266 + t22 * t340 + t25 * t716 + t267 * t3 + t29 * t715 + t705 * t87) * m(7) - (-Ifges(4,2) * t329 + t322 + t708) * t328 / 0.2e1 - t111 * mrSges(4,2) + t112 * mrSges(4,1) + t474 * t756 + t382 * t757 + (t363 * t689 + t364 * t688 + t491) * g(3) - (Ifges(4,1) * t328 - t598 + t690) * t329 / 0.2e1 - t142 * t187 - t141 * t188 + (m(5) * ((-t115 * t445 - t116 * t441) * qJD(4) + t695) + t445 * t94 - t441 * t93 - t188 * t544 - t187 * t546) * pkin(10) + (t115 * t576 + t116 * t577 + t695) * mrSges(5,3) + t700 * t100 + t705 * t99; t742 + (-Ifges(5,2) * t269 + t146 + t264) * t644 + (-t538 - t30) * t139 + (t537 - t31) * t137 + (t537 - t45) * t138 + (-t538 - t44) * t140 + t61 + (t40 + t39) * pkin(4) * t440 + (t440 * t5 + t444 * t6 + (-t35 * t440 + t36 * t444) * qJD(5)) * t675 + (-t730 * t637 - t659 * t733 - t731 * t662 + t762) * t164 + (t637 * t732 + t659 * t734 + t662 * t733 - t763) * t506 - t100 * t615 - m(6) * (t143 * t615 + t35 * t44 + t36 * t45) + t720 * t662 + t145 * t641 + (Ifges(5,1) * t268 - t595) * t642 + (Ifges(5,5) * t268 - Ifges(5,6) * t269) * t635 + t450 + (t745 * mrSges(5,2) - m(7) * (-t297 * t393 + t365 * t394) - t717 * t746 + t697) * g(2) + t430 * t37 + (-m(7) * (-t301 * t393 + t367 * t394) + mrSges(5,2) * t240 + t717 * t239 + t698) * g(1) - t182 * (mrSges(5,1) * t269 + mrSges(5,2) * t268) - t128 * t99 - t115 * t187 + t116 * t188 + (-m(7) * (-t364 * t393 - t394 * t568) - m(6) * t614 - t489 + t699) * g(3) + (-t128 * t87 - t25 * t30 - t29 * t31 + t2 * t430 + (t3 * t440 + (-t25 * t440 + t29 * t444) * qJD(5)) * pkin(4)) * m(7) + t269 * t589 + t268 * t590 + t38 * t611; t2 * t674 + t450 - t28 * t137 - t87 * (mrSges(7,1) * t164 + mrSges(7,2) * t506) - t143 * (mrSges(6,1) * t164 + mrSges(6,2) * t506) + t25 * t602 + (t506 * t734 - t748) * t659 + t721 * t658 + (-t164 * t730 + t506 * t732) * t637 + (t140 + t603) * t36 + (-t138 + t604) * t35 + (-m(7) * (-t25 + t28) + t139 + t601) * t29 + (-t282 * t674 + t699) * g(3) + (t674 * t747 + t697) * g(2) + (-t231 * t674 + t698) * g(1) + (-t164 * t731 + t720 + t752) * t662 + (t37 + (-m(7) * t87 - t99) * t164) * pkin(5); -t506 * t137 + t164 * t139 + (-g(1) * t300 - g(2) * t296 - g(3) * t363 + t164 * t25 - t29 * t506 + t22) * m(7) + t19;];
tau  = t1;
