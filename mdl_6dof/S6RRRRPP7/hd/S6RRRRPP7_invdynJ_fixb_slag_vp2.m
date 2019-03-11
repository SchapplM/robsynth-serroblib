% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:17:10
% EndTime: 2019-03-09 21:18:54
% DurationCPUTime: 67.60s
% Computational Cost: add. (23971->1038), mult. (57950->1373), div. (0->0), fcn. (46194->14), ass. (0->444)
t431 = sin(qJ(2));
t578 = cos(pkin(6));
t520 = pkin(1) * t578;
t412 = t431 * t520;
t430 = sin(qJ(3));
t433 = cos(qJ(3));
t481 = pkin(3) * t430 - pkin(10) * t433;
t427 = sin(pkin(6));
t434 = cos(qJ(2));
t557 = t427 * t434;
t780 = -(t412 + (pkin(8) + t481) * t557) * qJD(1) + t481 * qJD(3);
t496 = t578 * qJD(1);
t404 = t496 + qJD(2);
t539 = qJD(1) * t427;
t516 = t431 * t539;
t309 = t404 * t433 - t430 * t516;
t303 = qJD(4) - t309;
t618 = t303 / 0.2e1;
t310 = t404 * t430 + t433 * t516;
t515 = t434 * t539;
t387 = qJD(3) - t515;
t429 = sin(qJ(4));
t432 = cos(qJ(4));
t246 = -t310 * t429 + t387 * t432;
t247 = t310 * t432 + t387 * t429;
t426 = sin(pkin(11));
t577 = cos(pkin(11));
t456 = t426 * t246 + t247 * t577;
t633 = t456 / 0.2e1;
t146 = -t577 * t246 + t247 * t426;
t636 = t146 / 0.2e1;
t637 = -t146 / 0.2e1;
t712 = Ifges(7,4) + Ifges(6,5);
t714 = Ifges(6,1) + Ifges(7,1);
t713 = Ifges(6,4) - Ifges(7,5);
t708 = -t146 * t713 + t303 * t712 + t456 * t714;
t747 = t708 / 0.2e1;
t489 = pkin(1) * t496;
t343 = -pkin(8) * t516 + t434 * t489;
t289 = -t404 * pkin(2) - t343;
t166 = -t309 * pkin(3) - t310 * pkin(10) + t289;
t541 = pkin(8) * t557 + t412;
t346 = t541 * qJD(1);
t290 = t404 * pkin(9) + t346;
t301 = (-pkin(2) * t434 - pkin(9) * t431 - pkin(1)) * t539;
t192 = t433 * t290 + t430 * t301;
t171 = pkin(10) * t387 + t192;
t109 = t166 * t429 + t171 * t432;
t83 = qJ(5) * t246 + t109;
t580 = t426 * t83;
t108 = t432 * t166 - t171 * t429;
t82 = -qJ(5) * t247 + t108;
t71 = pkin(4) * t303 + t82;
t30 = t577 * t71 - t580;
t27 = -t303 * pkin(5) + qJD(6) - t30;
t191 = -t430 * t290 + t301 * t433;
t170 = -pkin(3) * t387 - t191;
t134 = -pkin(4) * t246 + qJD(5) + t170;
t53 = pkin(5) * t146 - qJ(6) * t456 + t134;
t754 = mrSges(7,2) * t27 - mrSges(6,3) * t30 - mrSges(7,3) * t53;
t779 = Ifges(6,4) * t637 + Ifges(7,5) * t636 + t712 * t618 + t714 * t633 + t747 + t754;
t461 = (pkin(2) * t431 - pkin(9) * t434) * t427;
t344 = qJD(1) * t461;
t237 = t433 * t343 + t430 * t344;
t213 = pkin(10) * t516 + t237;
t537 = qJD(3) * t430;
t528 = pkin(9) * t537;
t778 = t780 * t432 + (t213 + t528) * t429;
t383 = -pkin(3) * t433 - pkin(10) * t430 - pkin(2);
t533 = qJD(4) * t432;
t777 = -t432 * t213 + t383 * t533 + t429 * t780;
t705 = mrSges(6,2) * t134;
t776 = t705 + t779;
t79 = t577 * t83;
t31 = t426 * t71 + t79;
t28 = qJ(6) * t303 + t31;
t619 = -t303 / 0.2e1;
t634 = -t456 / 0.2e1;
t775 = Ifges(6,4) * t633 + Ifges(7,5) * t634 + Ifges(6,6) * t618 + Ifges(7,6) * t619 + (Ifges(6,2) + Ifges(7,3)) * t637 - t53 * mrSges(7,1) + mrSges(7,2) * t28 + mrSges(6,3) * t31;
t425 = qJ(4) + pkin(11);
t421 = sin(t425);
t422 = cos(t425);
t670 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t750 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t774 = t421 * t750 + t422 * t670;
t710 = Ifges(6,3) + Ifges(7,2);
t773 = Ifges(5,3) + t710;
t226 = pkin(3) * t310 - pkin(10) * t309;
t132 = -t191 * t429 + t432 * t226;
t428 = -qJ(5) - pkin(10);
t499 = qJD(4) * t428;
t571 = t309 * t432;
t771 = -pkin(4) * t310 + qJ(5) * t571 - qJD(5) * t429 + t432 * t499 - t132;
t133 = t432 * t191 + t429 * t226;
t532 = qJD(5) * t432;
t572 = t309 * t429;
t770 = -qJ(5) * t572 - t429 * t499 + t133 - t532;
t548 = t433 * t434;
t553 = t429 * t431;
t300 = (t432 * t548 + t553) * t539;
t549 = t432 * t433;
t414 = pkin(9) * t549;
t490 = t430 * t515;
t769 = -pkin(4) * t490 + t300 * qJ(5) - t430 * t532 + (pkin(4) * t430 - qJ(5) * t549) * qJD(3) + (-t414 + (qJ(5) * t430 - t383) * t429) * qJD(4) + t778;
t299 = (-t429 * t548 + t431 * t432) * t539;
t552 = t430 * t432;
t768 = qJ(5) * t299 - (-pkin(9) * qJD(3) - qJ(5) * qJD(4)) * t552 - (-qJD(5) * t430 + (-pkin(9) * qJD(4) - qJ(5) * qJD(3)) * t433) * t429 - t777;
t536 = qJD(3) * t433;
t688 = t429 * t536 + t430 * t533 + t299;
t767 = Ifges(6,2) * t637 - Ifges(7,3) * t636 + t633 * t713 + t775;
t702 = t387 * Ifges(4,6);
t765 = -t289 * mrSges(4,1) - t108 * mrSges(5,1) + t109 * mrSges(5,2) + t702 / 0.2e1;
t763 = -Ifges(6,4) * t636 - Ifges(7,5) * t637 - t712 * t619 - t634 * t714 + t754;
t531 = qJD(1) * qJD(2);
t350 = (qJDD(1) * t431 + t434 * t531) * t427;
t495 = t578 * qJDD(1);
t403 = t495 + qJDD(2);
t198 = -qJD(3) * t310 - t350 * t430 + t403 * t433;
t190 = qJDD(4) - t198;
t630 = t190 / 0.2e1;
t197 = qJD(3) * t309 + t350 * t433 + t403 * t430;
t349 = (-qJDD(1) * t434 + t431 * t531) * t427;
t338 = qJDD(3) + t349;
t123 = qJD(4) * t246 + t197 * t432 + t338 * t429;
t124 = -qJD(4) * t247 - t197 * t429 + t338 * t432;
t66 = t123 * t577 + t426 * t124;
t647 = t66 / 0.2e1;
t762 = t712 * t630 + t714 * t647;
t65 = t123 * t426 - t124 * t577;
t649 = -t65 / 0.2e1;
t719 = -m(7) - m(6);
t745 = pkin(4) * t719;
t711 = Ifges(6,6) - Ifges(7,6);
t761 = -mrSges(5,1) + t745;
t199 = -t299 * t577 + t300 * t426;
t497 = t577 * t432;
t498 = t577 * t429;
t534 = qJD(4) * t430;
t684 = t429 * t534 - t432 * t536;
t240 = t426 * t684 - t497 * t534 - t498 * t536;
t760 = -t199 - t240;
t200 = t426 * t299 + t300 * t577;
t372 = t426 * t432 + t498;
t559 = t426 * t429;
t455 = t497 - t559;
t241 = -t372 * t534 + t455 * t536;
t759 = t200 - t241;
t236 = -t430 * t343 + t344 * t433;
t212 = -pkin(3) * t516 - t236;
t420 = pkin(9) * t536;
t758 = -t212 + t420;
t685 = t490 - t537;
t607 = cos(qJ(1));
t483 = t578 * t607;
t606 = sin(qJ(1));
t361 = t431 * t483 + t434 * t606;
t519 = t427 * t607;
t275 = t361 * t433 - t430 * t519;
t360 = t431 * t606 - t434 * t483;
t757 = -t275 * t429 + t360 * t432;
t756 = t705 + t763;
t592 = Ifges(4,4) * t310;
t184 = t309 * Ifges(4,2) + t592 + t702;
t584 = t192 * mrSges(4,3);
t664 = t30 * mrSges(6,1) - t27 * mrSges(7,1) - t31 * mrSges(6,2) + t28 * mrSges(7,3);
t675 = t247 * Ifges(5,5) + t246 * Ifges(5,6) - t146 * t711 + t303 * t773 + t456 * t712;
t751 = t664 + Ifges(6,6) * t637 + Ifges(7,6) * t636 + t633 * t712 - t184 / 0.2e1 + t675 / 0.2e1 - t584;
t748 = -t708 / 0.2e1;
t746 = t649 * t713 + t762;
t744 = -m(4) + t719;
t527 = m(5) * pkin(10) + mrSges(5,3);
t715 = mrSges(6,3) + mrSges(7,2);
t677 = -t527 - t715;
t743 = mrSges(4,2) + t677;
t701 = t426 * t768 + t577 * t769;
t700 = t426 * t769 - t577 * t768;
t697 = t426 * t770 + t577 * t771;
t694 = t426 * t771 - t577 * t770;
t535 = qJD(4) * t429;
t691 = -t192 + (t535 - t572) * pkin(4);
t686 = pkin(4) * t688 + t758;
t477 = -mrSges(5,1) * t432 + mrSges(5,2) * t429;
t453 = m(5) * pkin(3) - t477;
t479 = t433 * mrSges(4,1) - t430 * mrSges(4,2);
t742 = -t430 * t527 - t433 * t453 - t479;
t530 = qJDD(1) * t427;
t741 = pkin(8) * t530 + qJD(2) * t489;
t740 = -pkin(8) * t427 * t531 + pkin(1) * t495;
t482 = t578 * t606;
t363 = -t431 * t482 + t434 * t607;
t518 = t427 * t606;
t279 = t363 * t433 + t430 * t518;
t362 = t431 * t607 + t434 * t482;
t214 = -t279 * t429 + t362 * t432;
t658 = -Ifges(6,2) * t636 + Ifges(7,3) * t637 - t634 * t713 + t775;
t704 = t134 * mrSges(6,1);
t734 = t658 - t704;
t733 = t704 - t767;
t201 = t275 * t421 - t360 * t422;
t731 = t275 * t422 + t360 * t421;
t568 = t360 * t429;
t730 = t275 * t432 + t568;
t476 = t429 * mrSges(5,1) + t432 * mrSges(5,2);
t716 = mrSges(4,3) - mrSges(3,2);
t723 = -m(5) * pkin(9) - t421 * t670 + t422 * t750 + t429 * t745 - t476 - t716;
t419 = pkin(4) * t432 + pkin(3);
t722 = mrSges(3,1) - t742 + (t428 * t719 + t715) * t430 + (-t419 * t719 + t774) * t433;
t248 = t431 * t740 + t434 * t741;
t221 = pkin(9) * t403 + t248;
t526 = pkin(1) * t530;
t233 = pkin(2) * t349 - pkin(9) * t350 - t526;
t101 = t433 * t221 + t430 * t233 - t290 * t537 + t301 * t536;
t87 = pkin(10) * t338 + t101;
t249 = -t431 * t741 + t434 * t740;
t222 = -t403 * pkin(2) - t249;
t98 = -t198 * pkin(3) - t197 * pkin(10) + t222;
t24 = -qJD(4) * t109 - t429 * t87 + t432 * t98;
t11 = pkin(4) * t190 - qJ(5) * t123 - qJD(5) * t247 + t24;
t23 = t166 * t533 - t171 * t535 + t429 * t98 + t432 * t87;
t13 = qJ(5) * t124 + qJD(5) * t246 + t23;
t4 = t426 * t11 + t577 * t13;
t1 = qJ(6) * t190 + qJD(6) * t303 + t4;
t3 = t11 * t577 - t426 * t13;
t2 = -t190 * pkin(5) + qJDD(6) - t3;
t721 = t24 * mrSges(5,1) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t23 * mrSges(5,2) - t4 * mrSges(6,2) + t1 * mrSges(7,3);
t102 = -t430 * t221 + t233 * t433 - t290 * t536 - t301 * t537;
t88 = -pkin(3) * t338 - t102;
t48 = -pkin(4) * t124 + qJDD(5) + t88;
t6 = pkin(5) * t65 - qJ(6) * t66 - qJD(6) * t456 + t48;
t648 = t65 / 0.2e1;
t720 = mrSges(6,2) * t48 + mrSges(7,2) * t2 - mrSges(6,3) * t3 - mrSges(7,3) * t6 + Ifges(6,4) * t649 + Ifges(7,5) * t648 + t746 + t762;
t640 = t123 / 0.2e1;
t639 = t124 / 0.2e1;
t628 = t197 / 0.2e1;
t627 = t198 / 0.2e1;
t613 = t338 / 0.2e1;
t707 = -qJ(6) * t685 - qJD(6) * t433 + t700;
t706 = pkin(5) * t685 - t701;
t703 = t387 * Ifges(4,5);
t340 = t455 * t430;
t699 = pkin(5) * t760 + qJ(6) * t759 - qJD(6) * t340 + t686;
t195 = t372 * t309;
t196 = t455 * t309;
t356 = t372 * qJD(4);
t357 = t455 * qJD(4);
t698 = -qJD(6) * t372 + t691 + (t196 - t357) * qJ(6) + (-t195 + t356) * pkin(5);
t696 = t310 * pkin(5) - t697;
t695 = -qJ(6) * t310 + t694;
t34 = t577 * t82 - t580;
t693 = qJD(6) - t34;
t130 = mrSges(6,1) * t303 - mrSges(6,3) * t456;
t131 = -mrSges(7,1) * t303 + mrSges(7,2) * t456;
t692 = t130 - t131;
t558 = t427 * t431;
t367 = -pkin(8) * t558 + t434 * t520;
t332 = -pkin(2) * t578 - t367;
t358 = t430 * t558 - t433 * t578;
t359 = t430 * t578 + t433 * t558;
t208 = t358 * pkin(3) - t359 * pkin(10) + t332;
t333 = pkin(9) * t578 + t541;
t542 = pkin(2) * t557 + pkin(9) * t558;
t601 = pkin(1) * t427;
t334 = -t542 - t601;
t228 = t433 * t333 + t430 * t334;
t210 = -pkin(10) * t557 + t228;
t127 = t429 * t208 + t432 * t210;
t690 = (-t432 * t537 - t433 * t535) * pkin(9) + t777;
t324 = t429 * t383 + t414;
t689 = -qJD(4) * t324 + t778;
t687 = t300 + t684;
t244 = Ifges(5,4) * t246;
t137 = t247 * Ifges(5,1) + t303 * Ifges(5,5) + t244;
t302 = Ifges(4,4) * t309;
t185 = t310 * Ifges(4,1) + t302 + t703;
t683 = t432 * t137 + t185;
t682 = t101 * t433 - t102 * t430;
t681 = t23 * t432 - t24 * t429;
t594 = Ifges(3,4) * t431;
t655 = t427 ^ 2;
t680 = (pkin(1) * (mrSges(3,1) * t431 + mrSges(3,2) * t434) - t431 * (Ifges(3,1) * t434 - t594) / 0.2e1) * t655;
t678 = Ifges(5,5) * t123 + Ifges(5,6) * t124 + t190 * t773 - t65 * t711 + t66 * t712;
t673 = t453 + t774;
t669 = mrSges(4,1) + t673;
t494 = mrSges(3,3) * t516;
t667 = -m(4) * t289 + mrSges(3,1) * t404 + mrSges(4,1) * t309 - mrSges(4,2) * t310 - t494;
t657 = mrSges(6,1) * t48 + mrSges(7,1) * t6 - mrSges(7,2) * t1 - mrSges(6,3) * t4 + 0.2e1 * Ifges(7,3) * t648 - t66 * Ifges(6,4) / 0.2e1 - t190 * Ifges(6,6) / 0.2e1 + Ifges(7,6) * t630 + (-t713 + Ifges(7,5)) * t647 + (-t649 + t648) * Ifges(6,2);
t622 = -t247 / 0.2e1;
t624 = -t246 / 0.2e1;
t656 = Ifges(5,5) * t622 + Ifges(5,6) * t624 + Ifges(6,6) * t636 + Ifges(7,6) * t637 + t634 * t712 - t664;
t46 = t123 * Ifges(5,4) + t124 * Ifges(5,2) + t190 * Ifges(5,6);
t651 = t46 / 0.2e1;
t650 = Ifges(5,1) * t640 + Ifges(5,4) * t639 + Ifges(5,5) * t630;
t641 = Ifges(4,1) * t628 + Ifges(4,4) * t627 + Ifges(4,5) * t613;
t638 = t137 / 0.2e1;
t623 = t246 / 0.2e1;
t621 = t247 / 0.2e1;
t616 = t309 / 0.2e1;
t614 = t310 / 0.2e1;
t600 = pkin(4) * t247;
t272 = -t359 * t429 - t432 * t557;
t449 = pkin(4) * t272;
t599 = pkin(4) * t426;
t597 = pkin(9) * t433;
t423 = t430 * pkin(9);
t538 = qJD(2) * t434;
t513 = t427 * t538;
t271 = -qJD(3) * t358 + t433 * t513;
t514 = qJD(2) * t558;
t165 = qJD(4) * t272 + t271 * t432 + t429 * t514;
t270 = qJD(3) * t359 + t430 * t513;
t460 = -t359 * t432 + t429 * t557;
t345 = qJD(2) * t461;
t347 = t367 * qJD(2);
t143 = -t333 * t537 + t334 * t536 + t430 * t345 + t433 * t347;
t139 = pkin(10) * t514 + t143;
t348 = t541 * qJD(2);
t157 = t270 * pkin(3) - t271 * pkin(10) + t348;
t52 = -qJD(4) * t127 - t139 * t429 + t432 * t157;
t26 = pkin(4) * t270 - qJ(5) * t165 + qJD(5) * t460 + t52;
t164 = qJD(4) * t460 - t271 * t429 + t432 * t514;
t51 = t432 * t139 + t429 * t157 + t208 * t533 - t210 * t535;
t32 = qJ(5) * t164 + qJD(5) * t272 + t51;
t9 = t426 * t26 + t577 * t32;
t593 = Ifges(3,4) * t434;
t591 = Ifges(4,4) * t430;
t590 = Ifges(4,4) * t433;
t589 = Ifges(5,4) * t429;
t588 = Ifges(5,4) * t432;
t587 = t108 * mrSges(5,3);
t586 = t109 * mrSges(5,3);
t585 = t191 * mrSges(4,3);
t581 = t247 * Ifges(5,4);
t579 = t430 * t88;
t110 = qJ(5) * t272 + t127;
t126 = t432 * t208 - t210 * t429;
t97 = pkin(4) * t358 + qJ(5) * t460 + t126;
t40 = t577 * t110 + t426 * t97;
t565 = t362 * t429;
t136 = t246 * Ifges(5,2) + t303 * Ifges(5,6) + t581;
t555 = t429 * t136;
t554 = t429 * t430;
t551 = t430 * t434;
t374 = t432 * t383;
t261 = -qJ(5) * t552 + t374 + (-pkin(9) * t429 - pkin(4)) * t433;
t282 = -qJ(5) * t554 + t324;
t173 = t426 * t261 + t577 * t282;
t377 = pkin(4) * t554 + t423;
t540 = t607 * pkin(1) + pkin(8) * t518;
t524 = t427 * t551;
t523 = t427 * t548;
t522 = Ifges(4,5) * t197 + Ifges(4,6) * t198 + Ifges(4,3) * t338;
t521 = Ifges(3,5) * t350 - Ifges(3,6) * t349 + Ifges(3,3) * t403;
t517 = t577 * pkin(4);
t508 = t558 / 0.2e1;
t507 = -t555 / 0.2e1;
t22 = t65 * mrSges(6,1) + t66 * mrSges(6,2);
t21 = t65 * mrSges(7,1) - t66 * mrSges(7,3);
t502 = -t534 / 0.2e1;
t44 = -t190 * mrSges(7,1) + t66 * mrSges(7,2);
t227 = -t430 * t333 + t334 * t433;
t493 = mrSges(3,3) * t515;
t484 = -pkin(1) * t606 + pkin(8) * t519;
t209 = pkin(3) * t557 - t227;
t480 = mrSges(4,1) * t358 + mrSges(4,2) * t359;
t478 = t272 * mrSges(5,1) + mrSges(5,2) * t460;
t474 = Ifges(4,1) * t433 - t591;
t473 = Ifges(5,1) * t432 - t589;
t472 = Ifges(5,1) * t429 + t588;
t471 = -Ifges(4,2) * t430 + t590;
t470 = -Ifges(5,2) * t429 + t588;
t469 = Ifges(5,2) * t432 + t589;
t468 = Ifges(4,5) * t433 - Ifges(4,6) * t430;
t467 = Ifges(5,5) * t432 - Ifges(5,6) * t429;
t466 = Ifges(5,5) * t429 + Ifges(5,6) * t432;
t463 = t363 * pkin(2) + pkin(9) * t362 + t540;
t144 = -t333 * t536 - t334 * t537 + t345 * t433 - t430 * t347;
t8 = t26 * t577 - t426 * t32;
t459 = t170 * t476;
t457 = (t434 * Ifges(3,2) + t594) * t427;
t39 = -t426 * t110 + t577 * t97;
t172 = t261 * t577 - t426 * t282;
t274 = t361 * t430 + t433 * t519;
t278 = t363 * t430 - t433 * t518;
t451 = -g(1) * t278 - g(2) * t274 - g(3) * t358;
t149 = t209 - t449;
t448 = t404 * t427 * (Ifges(3,5) * t434 - Ifges(3,6) * t431);
t443 = -t361 * pkin(2) - t360 * pkin(9) + t484;
t140 = -pkin(3) * t514 - t144;
t89 = -pkin(4) * t164 + t140;
t418 = -t517 - pkin(5);
t415 = qJ(6) + t599;
t398 = Ifges(3,4) * t515;
t389 = t428 * t432;
t364 = (-mrSges(3,1) * t434 + mrSges(3,2) * t431) * t427;
t353 = t362 * pkin(2);
t351 = t360 * pkin(2);
t342 = -t404 * mrSges(3,2) + t493;
t339 = t372 * t430;
t323 = -t429 * t597 + t374;
t288 = -t389 * t577 + t428 * t559;
t287 = -t389 * t426 - t428 * t498;
t285 = Ifges(3,1) * t516 + t404 * Ifges(3,5) + t398;
t284 = Ifges(3,6) * t404 + qJD(1) * t457;
t258 = t359 * t421 + t422 * t557;
t253 = mrSges(4,1) * t387 - mrSges(4,3) * t310;
t252 = -mrSges(4,2) * t387 + mrSges(4,3) * t309;
t250 = -pkin(5) * t455 - qJ(6) * t372 - t419;
t215 = t279 * t432 + t565;
t211 = pkin(5) * t339 - qJ(6) * t340 + t377;
t206 = t279 * t422 + t362 * t421;
t205 = t279 * t421 - t362 * t422;
t183 = t310 * Ifges(4,5) + t309 * Ifges(4,6) + t387 * Ifges(4,3);
t177 = mrSges(5,1) * t303 - mrSges(5,3) * t247;
t176 = -mrSges(5,2) * t303 + mrSges(5,3) * t246;
t175 = t426 * t272 - t460 * t577;
t174 = -t272 * t577 - t426 * t460;
t168 = t433 * pkin(5) - t172;
t167 = -qJ(6) * t433 + t173;
t156 = -mrSges(4,2) * t338 + mrSges(4,3) * t198;
t155 = mrSges(4,1) * t338 - mrSges(4,3) * t197;
t153 = -mrSges(5,1) * t246 + mrSges(5,2) * t247;
t129 = -mrSges(6,2) * t303 - mrSges(6,3) * t146;
t128 = -mrSges(7,2) * t146 + mrSges(7,3) * t303;
t125 = -mrSges(4,1) * t198 + mrSges(4,2) * t197;
t114 = t197 * Ifges(4,4) + t198 * Ifges(4,2) + t338 * Ifges(4,6);
t106 = t426 * t164 + t165 * t577;
t104 = -t164 * t577 + t165 * t426;
t91 = mrSges(6,1) * t146 + mrSges(6,2) * t456;
t90 = mrSges(7,1) * t146 - mrSges(7,3) * t456;
t85 = -mrSges(5,2) * t190 + mrSges(5,3) * t124;
t84 = mrSges(5,1) * t190 - mrSges(5,3) * t123;
t70 = pkin(5) * t456 + qJ(6) * t146 + t600;
t68 = pkin(5) * t174 - qJ(6) * t175 + t149;
t67 = -mrSges(5,1) * t124 + mrSges(5,2) * t123;
t43 = mrSges(6,1) * t190 - mrSges(6,3) * t66;
t42 = -mrSges(6,2) * t190 - mrSges(6,3) * t65;
t41 = -mrSges(7,2) * t65 + mrSges(7,3) * t190;
t36 = -t358 * pkin(5) - t39;
t35 = qJ(6) * t358 + t40;
t33 = t426 * t82 + t79;
t20 = pkin(5) * t104 - qJ(6) * t106 - qJD(6) * t175 + t89;
t7 = -t270 * pkin(5) - t8;
t5 = qJ(6) * t270 + qJD(6) * t358 + t9;
t10 = [(Ifges(5,4) * t165 + Ifges(5,2) * t164) * t623 + (Ifges(5,5) * t165 + Ifges(5,6) * t164 - t104 * t711) * t618 + (-m(3) * t540 - t363 * mrSges(3,1) - mrSges(3,3) * t518 - mrSges(2,1) * t607 + mrSges(2,2) * t606 - m(5) * (pkin(3) * t279 + t463) - t215 * mrSges(5,1) - t214 * mrSges(5,2) - m(4) * t463 - t279 * mrSges(4,1) + t719 * (pkin(4) * t565 - t278 * t428 + t279 * t419 + t463) - t716 * t362 - t670 * t206 - t750 * t205 + t743 * t278) * g(2) + (-m(5) * (-pkin(3) * t275 + t443) + t730 * mrSges(5,1) + t757 * mrSges(5,2) + mrSges(2,1) * t606 + mrSges(2,2) * t607 - m(4) * t443 + t275 * mrSges(4,1) - m(3) * t484 + t361 * mrSges(3,1) - mrSges(3,3) * t519 + t719 * (-pkin(4) * t568 + t274 * t428 - t275 * t419 + t443) + t716 * t360 + t670 * t731 + t750 * t201 - t743 * t274) * g(1) + (-Ifges(4,6) * t613 + Ifges(7,6) * t648 + Ifges(6,6) * t649 + Ifges(5,6) * t639 + Ifges(5,5) * t640 - Ifges(4,2) * t627 + t721 - t101 * mrSges(4,3) - Ifges(4,4) * t628 - t114 / 0.2e1 + t678 / 0.2e1 + t773 * t630 + t712 * t647) * t358 + (-Ifges(4,4) * t614 + Ifges(5,5) * t621 - Ifges(4,2) * t616 + Ifges(5,6) * t623 + t618 * t773 + t751 - t765) * t270 + t387 * (Ifges(4,5) * t271 + Ifges(4,3) * t514) / 0.2e1 + (Ifges(4,5) * t359 - Ifges(4,3) * t557) * t613 + (Ifges(5,1) * t165 + Ifges(5,4) * t164) * t621 + (Ifges(4,4) * t271 + Ifges(4,6) * t514) * t616 + (Ifges(4,4) * t359 - Ifges(4,6) * t557) * t627 + (-Ifges(5,4) * t460 + Ifges(5,2) * t272) * t639 + (-Ifges(5,1) * t460 + Ifges(5,4) * t272) * t640 + (-Ifges(5,5) * t460 + Ifges(5,6) * t272 - t174 * t711) * t630 + (-t108 * t165 + t109 * t164 + t23 * t272 + t24 * t460) * mrSges(5,3) - t460 * t650 + (Ifges(4,1) * t271 + Ifges(4,5) * t514) * t614 + (Ifges(4,1) * t359 - Ifges(4,5) * t557) * t628 + (t655 * qJD(1) * (-Ifges(3,2) * t431 + t593) + t427 * t285) * t538 / 0.2e1 + (-t248 * t578 - t350 * t601 - t403 * t541) * mrSges(3,2) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t655 + t248 * t541 + t249 * t367 + t346 * t347) + m(7) * (t1 * t35 + t2 * t36 + t20 * t53 + t27 * t7 + t28 * t5 + t6 * t68) + m(6) * (t134 * t89 + t149 * t48 + t3 * t39 + t30 * t8 + t31 * t9 + t4 * t40) + m(5) * (t108 * t52 + t109 * t51 + t126 * t24 + t127 * t23 + t140 * t170 + t209 * t88) + t350 * (Ifges(3,5) * t578 + (t431 * Ifges(3,1) + t593) * t427) / 0.2e1 + t657 * t174 - t284 * t514 / 0.2e1 + t191 * (mrSges(4,1) * t514 - mrSges(4,3) * t271) + t733 * t104 + (t448 / 0.2e1 + t183 * t508) * qJD(2) + (t248 * t557 - t249 * t558 - t343 * t513 - t346 * t514 - t349 * t541 - t350 * t367) * mrSges(3,3) + (Ifges(3,1) * t350 - Ifges(3,4) * t349 + Ifges(3,5) * t403) * t508 + (t249 * t578 - t349 * t601 + t367 * t403) * mrSges(3,1) + (Ifges(3,4) * t350 - Ifges(3,2) * t349 + Ifges(3,6) * t403) * t557 / 0.2e1 - t364 * t526 + t720 * t175 - t522 * t557 / 0.2e1 + t102 * (-mrSges(4,1) * t557 - t359 * mrSges(4,3)) + m(4) * (t101 * t228 + t102 * t227 + t143 * t192 + t144 * t191 + t222 * t332) + t776 * t106 + Ifges(2,3) * qJDD(1) + t347 * t342 + t332 * t125 + t271 * t185 / 0.2e1 + t143 * t252 + t144 * t253 + t227 * t155 + t228 * t156 + t209 * t67 + t51 * t176 + t52 * t177 + t170 * (-mrSges(5,1) * t164 + mrSges(5,2) * t165) + t164 * t136 / 0.2e1 + t149 * t22 + t140 * t153 + t8 * t130 + t7 * t131 + t126 * t84 + t127 * t85 + t5 * t128 + t9 * t129 + t20 * t90 + t89 * t91 + t68 * t21 + t40 * t42 + t39 * t43 + t36 * t44 + t35 * t41 + t272 * t651 + t359 * t641 - t88 * t478 + t222 * t480 + (t101 * t557 - t192 * t514 + t271 * t289) * mrSges(4,2) + t578 * t521 / 0.2e1 - t349 * (Ifges(3,6) * t578 + t457) / 0.2e1 + t403 * (Ifges(3,3) * t578 + (Ifges(3,5) * t431 + Ifges(3,6) * t434) * t427) / 0.2e1 + (-m(3) * t343 - t667) * t348 + t165 * t638 - t680 * t531; (t507 - t585) * t536 + (-t528 - t237) * t252 + (t364 + (-m(5) - m(4)) * t542 - t715 * t524 + t719 * (t419 * t523 - t428 * t524 + t542) - t750 * (t421 * t523 - t422 * t558) + (t553 * t745 + t742 * t434 + (-t476 - mrSges(4,3)) * t431 - t670 * (t421 * t431 + t422 * t548)) * t427) * g(3) + (Ifges(5,5) * t300 + Ifges(5,6) * t299 - t199 * t711 + t490 * t773) * t619 + (-t339 * t711 + t340 * t712 + t430 * t467 - t433 * t773) * t630 + (-t448 / 0.2e1 + t284 * t508 + t680 * qJD(1)) * qJD(1) + (-t1 * t433 - t340 * t6) * mrSges(7,3) + (-t420 - t236) * t253 + (Ifges(6,4) * t340 - Ifges(6,6) * t433) * t649 - (t387 * (Ifges(4,3) * t431 + t434 * t468) + t310 * (Ifges(4,5) * t431 + t434 * t474) + t309 * (Ifges(4,6) * t431 + t434 * t471) + t431 * t183 + (-Ifges(3,2) * t516 + t433 * t185 + t430 * t675 + t285 + t398) * t434) * t539 / 0.2e1 + (-t466 * t534 + (Ifges(5,3) * t430 + t433 * t467) * qJD(3) + t711 * t240) * t618 + (-t155 + t67) * t423 + (t309 * t471 + t310 * t474 + t387 * t468) * qJD(3) / 0.2e1 + (Ifges(5,1) * t300 + Ifges(5,4) * t299) * t622 + (t748 - t763) * t200 + (t340 * t48 + t4 * t433) * mrSges(6,2) + t340 * t746 + t476 * t579 + (Ifges(7,5) * t340 - Ifges(7,6) * t433) * t648 + (-t342 + t493) * t343 + t657 * t339 + t658 * t199 + (-t191 * (mrSges(4,1) * t431 - mrSges(4,3) * t548) - t192 * (-mrSges(4,2) * t431 - mrSges(4,3) * t551)) * t539 + (t184 / 0.2e1 + t656) * t490 + t156 * t597 + (Ifges(4,5) * t430 + Ifges(4,6) * t433) * t613 + t521 + (-t299 / 0.2e1 + t432 * t502) * t136 + (-t300 / 0.2e1 + t429 * t502) * t137 + t324 * t85 + t323 * t84 - t248 * mrSges(3,2) + t249 * mrSges(3,1) + t211 * t21 + t172 * t43 + t173 * t42 + t167 * t41 + t168 * t44 - pkin(2) * t125 + (t340 * t714 - t433 * t712) * t647 + t706 * t131 + t707 * t128 + (t1 * t167 + t168 * t2 + t211 * t6 + t27 * t706 + t28 * t707 + t53 * t699) * m(7) + t552 * t650 + (-Ifges(5,6) * t433 + t430 * t470) * t639 + (-Ifges(5,5) * t433 + t430 * t473) * t640 + t430 * t641 + t699 * t90 + t700 * t129 + t701 * t130 + (t134 * t686 + t172 * t3 + t173 * t4 + t30 * t701 + t31 * t700 + t377 * t48) * m(6) + (-t472 * t534 + (Ifges(5,5) * t430 + t433 * t473) * qJD(3)) * t621 + (-t469 * t534 + (Ifges(5,6) * t430 + t433 * t470) * qJD(3)) * t623 + (Ifges(4,2) * t433 + t591) * t627 + (Ifges(4,1) * t430 + t590) * t628 + t387 * t289 * (mrSges(4,1) * t430 + mrSges(4,2) * t433) - t222 * t479 + t767 * t240 + t377 * t22 + (t618 * t710 + t751) * t537 + (t494 + t667) * t346 + t779 * t241 + (Ifges(5,4) * t300 + Ifges(5,2) * t299) * t624 + t24 * (-mrSges(5,1) * t433 - mrSges(5,3) * t552) - t46 * t554 / 0.2e1 + t23 * (mrSges(5,2) * t433 - mrSges(5,3) * t554) + t758 * t153 + (mrSges(6,1) * t760 - mrSges(6,2) * t759) * t134 + (m(5) * t353 + t744 * (pkin(9) * t363 - t353) + t723 * t363 + t722 * t362) * g(1) + (m(5) * t351 + t744 * (pkin(9) * t361 - t351) + t723 * t361 + t722 * t360) * g(2) - t678 * t433 / 0.2e1 + t2 * (mrSges(7,1) * t433 + mrSges(7,2) * t340) + t3 * (-mrSges(6,1) * t433 - mrSges(6,3) * t340) + t433 * t114 / 0.2e1 + (-pkin(2) * t222 + ((-t191 * t433 - t192 * t430) * qJD(3) + t682) * pkin(9) - t191 * t236 - t192 * t237) * m(4) + t682 * mrSges(4,3) + t683 * t536 / 0.2e1 + t686 * t91 + (-mrSges(5,1) * t685 + mrSges(5,3) * t687) * t108 + (mrSges(5,1) * t688 - mrSges(5,2) * t687) * t170 + (mrSges(5,2) * t685 - mrSges(5,3) * t688) * t109 + t689 * t177 + t690 * t176 + (-t170 * t212 + t23 * t324 + t24 * t323 + (t170 * t536 + t579) * pkin(9) + t690 * t109 + t689 * t108) * m(5); (t42 + t41) * t288 + (t44 - t43) * t287 + (-t195 * t711 + t309 * t467 + t310 * t773) * t619 + (t507 + t459) * qJD(4) + (t719 * (-t278 * t419 - t279 * t428) + t743 * t279 + t669 * t278) * g(1) + (t719 * (-t274 * t419 - t275 * t428) + t743 * t275 + t669 * t274) * g(2) - (-t630 * t711 + t657) * t455 + t734 * t195 + (t246 * t470 + t247 * t473 + t303 * t467) * qJD(4) / 0.2e1 + (t253 - t153) * t192 + (-t587 + t638) * t533 + (-t618 * t711 + t733) * t356 + t184 * t614 + t555 * t616 - t535 * t586 + t522 + (-pkin(3) * t88 - t108 * t132 - t109 * t133 - t170 * t192) * m(5) + t720 * t372 + t776 * t357 + t250 * t21 - t191 * t252 - t133 * t176 - t132 * t177 + t102 * mrSges(4,1) - t101 * mrSges(4,2) - pkin(3) * t67 + (t480 + t719 * (-t358 * t419 - t359 * t428) + t677 * t359 + t673 * t358) * g(3) + t429 * t650 + t432 * t651 + t469 * t639 + t472 * t640 + t697 * t130 + (t134 * t691 - t287 * t3 + t288 * t4 + t30 * t697 + t31 * t694 - t419 * t48) * m(6) + t698 * t90 + (t1 * t288 + t2 * t287 + t250 * t6 + t27 * t696 + t28 * t695 + t53 * t698) * m(7) + (t585 - t459 - t289 * mrSges(4,2) + t473 * t622 + t470 * t624 - t703 / 0.2e1) * t309 + t694 * t129 + t695 * t128 + t696 * t131 + t691 * t91 + t88 * t477 + (t584 + t656 + t765) * t310 - t419 * t22 + (t748 - t756) * t196 + t466 * t630 - (Ifges(4,1) * t309 - t592 + t675) * t310 / 0.2e1 + (-t177 * t533 - t176 * t535 + m(5) * ((-t108 * t432 - t109 * t429) * qJD(4) + t681) - t429 * t84 + t432 * t85) * pkin(10) + (t108 * t571 + t109 * t572 + t681) * mrSges(5,3) - (-Ifges(4,2) * t310 + t302 + t683) * t309 / 0.2e1; (-t478 + t719 * t449 - t750 * (t359 * t422 - t421 * t557) + t670 * t258) * g(3) + (mrSges(5,2) * t215 + t670 * t205 - t206 * t750 + t214 * t761) * g(1) + (mrSges(5,2) * t730 + t670 * t201 - t731 * t750 + t757 * t761) * g(2) + t721 + (-t619 * t711 + t734) * t456 + ((t3 * t577 + t4 * t426) * pkin(4) - t134 * t600 + t30 * t33 - t31 * t34) * m(6) + t247 * t586 + t246 * t587 + t42 * t599 + t678 + (-Ifges(5,2) * t247 + t137 + t244) * t624 + t43 * t517 - t170 * (mrSges(5,1) * t247 + mrSges(5,2) * t246) - t108 * t176 + t109 * t177 - t34 * t129 - t70 * t90 + t692 * t33 + (t1 * t415 + t2 * t418 - t27 * t33 + t28 * t693 - t53 * t70) * m(7) + t693 * t128 - t91 * t600 + (Ifges(5,5) * t246 - Ifges(5,6) * t247) * t619 + t136 * t621 + (Ifges(5,1) * t246 - t581) * t622 + t415 * t41 + t418 * t44 + (t747 + t756) * t146; -(-t128 - t129) * t146 + t692 * t456 + t21 + t22 + (t146 * t28 - t27 * t456 + t451 + t6) * m(7) + (t146 * t31 + t30 * t456 + t451 + t48) * m(6); -t303 * t128 + t456 * t90 + (-g(1) * t205 - g(2) * t201 - g(3) * t258 - t28 * t303 + t456 * t53 + t2) * m(7) + t44;];
tau  = t10;
