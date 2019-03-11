% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR8
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
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:44:09
% EndTime: 2019-03-09 18:45:58
% DurationCPUTime: 67.28s
% Computational Cost: add. (36593->1123), mult. (87590->1513), div. (0->0), fcn. (71721->18), ass. (0->483)
t439 = cos(qJ(2));
t434 = sin(qJ(2));
t427 = sin(pkin(6));
t537 = qJD(1) * t427;
t506 = t434 * t537;
t429 = cos(pkin(6));
t536 = qJD(1) * t429;
t524 = pkin(1) * t536;
t349 = -pkin(8) * t506 + t439 * t524;
t459 = (pkin(2) * t434 - pkin(9) * t439) * t427;
t350 = qJD(1) * t459;
t433 = sin(qJ(3));
t438 = cos(qJ(3));
t261 = -t349 * t433 + t438 * t350;
t430 = -qJ(4) - pkin(9);
t495 = qJD(3) * t430;
t544 = t438 * t439;
t753 = -(pkin(3) * t434 - qJ(4) * t544) * t537 - t261 - qJD(4) * t433 + t438 * t495;
t262 = t438 * t349 + t433 * t350;
t505 = t439 * t537;
t485 = t433 * t505;
t752 = -qJ(4) * t485 - qJD(4) * t438 - t433 * t495 + t262;
t426 = sin(pkin(12));
t428 = cos(pkin(12));
t375 = t426 * t438 + t428 * t433;
t301 = t375 * t505;
t361 = t375 * qJD(3);
t751 = -t301 + t361;
t388 = qJD(3) - t505;
t610 = -t388 / 0.2e1;
t406 = qJD(2) + t536;
t323 = t406 * t438 - t433 * t506;
t324 = t406 * t433 + t438 * t506;
t463 = t323 * t426 + t428 * t324;
t619 = -t463 / 0.2e1;
t690 = t428 * t323 - t426 * t324;
t622 = -t690 / 0.2e1;
t750 = Ifges(5,4) * t619 + Ifges(5,2) * t622 + Ifges(5,6) * t610;
t437 = cos(qJ(5));
t417 = pkin(5) * t437 + pkin(4);
t425 = qJ(5) + qJ(6);
t421 = sin(t425);
t422 = cos(t425);
t432 = sin(qJ(5));
t477 = -mrSges(6,1) * t437 + mrSges(6,2) * t432;
t680 = m(6) * pkin(4) + m(7) * t417 + mrSges(7,1) * t422 - mrSges(7,2) * t421 + mrSges(5,1) - t477;
t672 = mrSges(5,2) + m(7) * (-pkin(11) - pkin(10)) - mrSges(7,3) - m(6) * pkin(10);
t695 = t426 * t753 - t752 * t428;
t352 = pkin(8) * t505 + t434 * t524;
t533 = qJD(3) * t433;
t686 = -t352 + (-t485 + t533) * pkin(3);
t448 = -mrSges(6,3) + t672;
t749 = pkin(10) * t506 - t695;
t374 = t426 * t433 - t428 * t438;
t302 = t374 * t505;
t362 = t374 * qJD(3);
t748 = t686 + (-t302 + t362) * pkin(10) + t751 * pkin(4);
t309 = pkin(9) * t406 + t352;
t338 = (-pkin(2) * t439 - pkin(9) * t434 - pkin(1)) * t427;
t314 = qJD(1) * t338;
t231 = -t309 * t433 + t438 * t314;
t194 = -qJ(4) * t324 + t231;
t181 = pkin(3) * t388 + t194;
t232 = t309 * t438 + t314 * t433;
t195 = qJ(4) * t323 + t232;
t550 = t428 * t195;
t110 = t426 * t181 + t550;
t436 = cos(qJ(6));
t241 = qJD(5) - t690;
t207 = t388 * t432 + t437 * t463;
t105 = pkin(10) * t388 + t110;
t308 = -pkin(2) * t406 - t349;
t242 = -pkin(3) * t323 + qJD(4) + t308;
t121 = -pkin(4) * t690 - pkin(10) * t463 + t242;
t55 = -t105 * t432 + t437 * t121;
t48 = -pkin(11) * t207 + t55;
t45 = pkin(5) * t241 + t48;
t431 = sin(qJ(6));
t206 = t388 * t437 - t432 * t463;
t56 = t105 * t437 + t121 * t432;
t49 = pkin(11) * t206 + t56;
t567 = t431 * t49;
t16 = t436 * t45 - t567;
t566 = t436 * t49;
t17 = t431 * t45 + t566;
t747 = mrSges(5,1) * t242 + mrSges(6,1) * t55 + mrSges(7,1) * t16 - mrSges(6,2) * t56 - mrSges(7,2) * t17 - mrSges(5,3) * t110 + t750;
t266 = t302 * t432 + t437 * t506;
t530 = qJD(5) * t437;
t502 = t375 * t530;
t456 = -t362 * t432 + t502;
t694 = t266 + t456;
t609 = t388 / 0.2e1;
t618 = t463 / 0.2e1;
t621 = t690 / 0.2e1;
t624 = t241 / 0.2e1;
t238 = qJD(6) + t241;
t626 = t238 / 0.2e1;
t633 = t207 / 0.2e1;
t635 = t206 / 0.2e1;
t135 = t206 * t431 + t207 * t436;
t646 = t135 / 0.2e1;
t493 = t436 * t206 - t207 * t431;
t648 = t493 / 0.2e1;
t746 = -Ifges(5,4) * t618 + Ifges(6,5) * t633 + Ifges(7,5) * t646 - Ifges(5,2) * t621 - Ifges(5,6) * t609 + Ifges(6,6) * t635 + Ifges(7,6) * t648 + Ifges(6,3) * t624 + Ifges(7,3) * t626 + t747;
t424 = qJ(3) + pkin(12);
t419 = sin(t424);
t420 = cos(t424);
t479 = -mrSges(4,1) * t438 + mrSges(4,2) * t433;
t745 = -m(4) * pkin(2) + t448 * t419 - t420 * t680 + t479;
t744 = Ifges(4,3) + Ifges(5,3);
t743 = t432 * t749 + t437 * t748;
t418 = pkin(3) * t438 + pkin(2);
t282 = pkin(4) * t374 - pkin(10) * t375 - t418;
t392 = t430 * t433;
t393 = t430 * t438;
t307 = t392 * t426 - t393 * t428;
t531 = qJD(5) * t432;
t703 = t282 * t530 - t307 * t531 + t432 * t748 - t437 * t749;
t529 = qJD(1) * qJD(2);
t356 = (qJDD(1) * t434 + t439 * t529) * t427;
t527 = qJDD(1) * t429;
t405 = qJDD(2) + t527;
t234 = -qJD(3) * t324 - t356 * t433 + t405 * t438;
t528 = qJDD(1) * t427;
t718 = pkin(8) * t528 + qJD(2) * t524;
t719 = -pkin(8) * t427 * t529 + pkin(1) * t527;
t278 = -t434 * t718 + t439 * t719;
t252 = -pkin(2) * t405 - t278;
t177 = -pkin(3) * t234 + qJDD(4) + t252;
t277 = t434 * t719 + t439 * t718;
t251 = pkin(9) * t405 + t277;
t355 = (-qJDD(1) * t439 + t434 * t529) * t427;
t256 = -pkin(1) * t528 + pkin(2) * t355 - pkin(9) * t356;
t119 = -t232 * qJD(3) - t251 * t433 + t438 * t256;
t233 = qJD(3) * t323 + t356 * t438 + t405 * t433;
t342 = qJDD(3) + t355;
t82 = pkin(3) * t342 - qJ(4) * t233 - qJD(4) * t324 + t119;
t532 = qJD(3) * t438;
t118 = t438 * t251 + t433 * t256 - t309 * t533 + t314 * t532;
t86 = qJ(4) * t234 + qJD(4) * t323 + t118;
t43 = -t426 * t86 + t428 * t82;
t613 = t342 / 0.2e1;
t167 = t233 * t428 + t234 * t426;
t640 = t167 / 0.2e1;
t166 = -t426 * t233 + t234 * t428;
t641 = t166 / 0.2e1;
t740 = t177 * mrSges(5,2) - mrSges(5,3) * t43 + 0.2e1 * Ifges(5,1) * t640 + 0.2e1 * Ifges(5,4) * t641 + 0.2e1 * Ifges(5,5) * t613;
t165 = qJDD(5) - t166;
t93 = qJD(5) * t206 + t167 * t437 + t342 * t432;
t94 = -qJD(5) * t207 - t167 * t432 + t342 * t437;
t37 = Ifges(6,5) * t93 + Ifges(6,6) * t94 + Ifges(6,3) * t165;
t44 = t426 * t82 + t428 * t86;
t642 = t165 / 0.2e1;
t159 = qJDD(6) + t165;
t643 = t159 / 0.2e1;
t653 = t94 / 0.2e1;
t654 = t93 / 0.2e1;
t34 = -qJD(6) * t135 - t431 * t93 + t436 * t94;
t663 = t34 / 0.2e1;
t33 = qJD(6) * t493 + t431 * t94 + t436 * t93;
t664 = t33 / 0.2e1;
t42 = pkin(10) * t342 + t44;
t59 = -pkin(4) * t166 - pkin(10) * t167 + t177;
t11 = -t105 * t531 + t121 * t530 + t437 * t42 + t432 * t59;
t12 = -qJD(5) * t56 - t42 * t432 + t437 * t59;
t676 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t6 = pkin(5) * t165 - pkin(11) * t93 + t12;
t7 = pkin(11) * t94 + t11;
t2 = qJD(6) * t16 + t431 * t6 + t436 * t7;
t3 = -qJD(6) * t17 - t431 * t7 + t436 * t6;
t677 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t8 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t159;
t739 = mrSges(5,1) * t177 + Ifges(6,5) * t654 + Ifges(7,5) * t664 + Ifges(6,6) * t653 + Ifges(7,6) * t663 + Ifges(6,3) * t642 + Ifges(7,3) * t643 + t37 / 0.2e1 + t8 / 0.2e1 - mrSges(5,3) * t44 + (-t342 / 0.2e1 - t613) * Ifges(5,6) + (-t166 / 0.2e1 - t641) * Ifges(5,2) + (-t167 / 0.2e1 - t640) * Ifges(5,4) + t677 + t676;
t267 = -t302 * t437 + t432 * t506;
t292 = t437 * t307;
t559 = t362 * t437;
t737 = pkin(11) * t267 + pkin(11) * t559 + (-t292 + (pkin(11) * t375 - t282) * t432) * qJD(5) + t743 + t751 * pkin(5);
t736 = pkin(11) * t694 - t703;
t602 = pkin(3) * t426;
t415 = pkin(10) + t602;
t594 = pkin(11) + t415;
t494 = qJD(5) * t594;
t561 = t690 * t432;
t190 = t426 * t195;
t114 = t194 * t428 - t190;
t603 = pkin(3) * t324;
t149 = pkin(4) * t463 - pkin(10) * t690 + t603;
t68 = t437 * t114 + t432 * t149;
t735 = -pkin(11) * t561 + t432 * t494 + t68;
t560 = t690 * t437;
t67 = -t114 * t432 + t437 * t149;
t734 = -pkin(5) * t463 + pkin(11) * t560 - t437 * t494 - t67;
t697 = t752 * t426 + t428 * t753;
t731 = t531 - t561;
t109 = t181 * t428 - t190;
t728 = t242 * mrSges(5,2) - t109 * mrSges(5,3);
t171 = Ifges(5,1) * t463 + Ifges(5,4) * t690 + t388 * Ifges(5,5);
t637 = t171 / 0.2e1;
t727 = Ifges(5,1) * t618 + Ifges(5,4) * t621 + Ifges(5,5) * t609 + t637;
t666 = m(7) * pkin(5);
t585 = Ifges(6,4) * t207;
t102 = t206 * Ifges(6,2) + t241 * Ifges(6,6) + t585;
t724 = -t102 / 0.2e1;
t712 = -m(6) - m(7);
t683 = m(5) - t712;
t722 = pkin(3) * t683;
t721 = t277 * mrSges(3,2);
t720 = t207 * Ifges(6,5) + t135 * Ifges(7,5) + t206 * Ifges(6,6) + Ifges(7,6) * t493 + t241 * Ifges(6,3) + t238 * Ifges(7,3);
t696 = pkin(4) * t506 - t697;
t717 = Ifges(4,5) * t233 + Ifges(5,5) * t167 + Ifges(4,6) * t234 + Ifges(5,6) * t166 + t342 * t744;
t607 = cos(qJ(1));
t507 = t607 * t439;
t435 = sin(qJ(1));
t547 = t434 * t435;
t368 = -t429 * t547 + t507;
t552 = t427 * t438;
t298 = -t368 * t433 + t435 * t552;
t554 = t427 * t434;
t363 = t429 * t438 - t433 * t554;
t625 = -t241 / 0.2e1;
t627 = -t238 / 0.2e1;
t634 = -t207 / 0.2e1;
t636 = -t206 / 0.2e1;
t647 = -t135 / 0.2e1;
t649 = -t493 / 0.2e1;
t716 = Ifges(6,5) * t634 + Ifges(7,5) * t647 + Ifges(6,6) * t636 + Ifges(7,6) * t649 + Ifges(6,3) * t625 + Ifges(7,3) * t627 - t747 - t750;
t475 = t421 * mrSges(7,1) + t422 * mrSges(7,2);
t565 = t437 * mrSges(6,2);
t476 = t432 * mrSges(6,1) + t565;
t526 = t432 * t666;
t709 = m(4) * pkin(9) + mrSges(4,3) + mrSges(5,3);
t715 = -t475 - t526 - t476 - t709;
t714 = -t119 * mrSges(4,1) - t43 * mrSges(5,1) + t118 * mrSges(4,2) + t44 * mrSges(5,2);
t629 = t233 / 0.2e1;
t628 = t234 / 0.2e1;
t204 = t437 * t282 - t307 * t432;
t555 = t375 * t437;
t179 = pkin(5) * t374 - pkin(11) * t555 + t204;
t205 = t432 * t282 + t292;
t556 = t375 * t432;
t188 = -pkin(11) * t556 + t205;
t99 = t179 * t436 - t188 * t431;
t711 = qJD(6) * t99 + t431 * t737 - t436 * t736;
t100 = t179 * t431 + t188 * t436;
t710 = -qJD(6) * t100 + t431 * t736 + t436 * t737;
t372 = t594 * t432;
t373 = t594 * t437;
t280 = -t372 * t436 - t373 * t431;
t708 = qJD(6) * t280 + t431 * t734 - t436 * t735;
t281 = -t372 * t431 + t373 * t436;
t707 = -qJD(6) * t281 + t431 * t735 + t436 * t734;
t113 = t194 * t426 + t550;
t706 = pkin(5) * t731 - t113;
t705 = -t666 - mrSges(6,1);
t704 = -qJD(5) * t205 + t743;
t462 = t431 * t432 - t436 * t437;
t272 = t462 * t375;
t551 = t427 * t439;
t371 = t429 * t434 * pkin(1) + pkin(8) * t551;
t337 = pkin(9) * t429 + t371;
t254 = -t337 * t433 + t438 * t338;
t364 = t429 * t433 + t434 * t552;
t201 = -pkin(3) * t551 - qJ(4) * t364 + t254;
t255 = t438 * t337 + t433 * t338;
t215 = qJ(4) * t363 + t255;
t131 = t426 * t201 + t428 * t215;
t125 = -pkin(10) * t551 + t131;
t269 = -t428 * t363 + t364 * t426;
t270 = t363 * t426 + t364 * t428;
t407 = pkin(8) * t554;
t604 = pkin(1) * t439;
t336 = t407 + (-pkin(2) - t604) * t429;
t279 = -pkin(3) * t363 + t336;
t168 = pkin(4) * t269 - pkin(10) * t270 + t279;
t76 = t437 * t125 + t432 * t168;
t378 = t431 * t437 + t432 * t436;
t682 = qJD(5) + qJD(6);
t294 = t682 * t378;
t150 = -t294 * t375 + t362 * t462;
t185 = t266 * t431 + t267 * t436;
t702 = t150 - t185;
t151 = t272 * t682 + t378 * t362;
t184 = t266 * t436 - t267 * t431;
t701 = t151 - t184;
t700 = pkin(5) * t694 + t696;
t139 = -mrSges(6,1) * t206 + mrSges(6,2) * t207;
t214 = mrSges(5,1) * t388 - mrSges(5,3) * t463;
t699 = t214 - t139;
t490 = mrSges(3,3) * t506;
t698 = -mrSges(3,1) * t406 - mrSges(4,1) * t323 + mrSges(4,2) * t324 + t490;
t455 = t375 * t531 + t559;
t693 = t267 + t455;
t174 = t462 * t690;
t293 = t682 * t462;
t692 = -t293 + t174;
t173 = t378 * t690;
t691 = -t294 + t173;
t370 = t429 * t604 - t407;
t508 = t607 * t434;
t546 = t435 * t439;
t366 = t429 * t508 + t546;
t509 = t427 * t607;
t446 = t366 * t433 + t438 * t509;
t685 = t11 * t437 - t12 * t432;
t104 = -pkin(4) * t388 - t109;
t81 = -pkin(5) * t206 + t104;
t679 = -t81 * mrSges(7,1) + mrSges(7,3) * t17;
t678 = t81 * mrSges(7,2) - mrSges(7,3) * t16;
t675 = mrSges(3,2) + t715;
t674 = mrSges(3,1) - t745;
t673 = t432 * t724 + t728;
t669 = t427 ^ 2;
t668 = Ifges(7,4) * t664 + Ifges(7,2) * t663 + Ifges(7,6) * t643;
t667 = m(5) * pkin(3);
t665 = Ifges(7,1) * t664 + Ifges(7,4) * t663 + Ifges(7,5) * t643;
t38 = t93 * Ifges(6,4) + t94 * Ifges(6,2) + t165 * Ifges(6,6);
t662 = t38 / 0.2e1;
t661 = Ifges(6,1) * t654 + Ifges(6,4) * t653 + Ifges(6,5) * t642;
t582 = Ifges(7,4) * t135;
t61 = Ifges(7,2) * t493 + Ifges(7,6) * t238 + t582;
t660 = -t61 / 0.2e1;
t659 = t61 / 0.2e1;
t127 = Ifges(7,4) * t493;
t62 = Ifges(7,1) * t135 + Ifges(7,5) * t238 + t127;
t658 = -t62 / 0.2e1;
t657 = t62 / 0.2e1;
t652 = pkin(1) * mrSges(3,1);
t651 = pkin(1) * mrSges(3,2);
t645 = Ifges(4,4) * t629 + Ifges(4,2) * t628 + Ifges(4,6) * t613;
t644 = Ifges(4,1) * t629 + Ifges(4,4) * t628 + Ifges(4,5) * t613;
t571 = t324 * Ifges(4,4);
t225 = t323 * Ifges(4,2) + t388 * Ifges(4,6) + t571;
t631 = t225 / 0.2e1;
t315 = Ifges(4,4) * t323;
t226 = t324 * Ifges(4,1) + t388 * Ifges(4,5) + t315;
t630 = t226 / 0.2e1;
t614 = t324 / 0.2e1;
t608 = t429 / 0.2e1;
t601 = pkin(3) * t428;
t600 = pkin(5) * t207;
t593 = mrSges(4,3) * t323;
t591 = mrSges(6,3) * t206;
t590 = mrSges(6,3) * t207;
t589 = Ifges(3,4) * t434;
t588 = Ifges(3,4) * t439;
t587 = Ifges(4,4) * t433;
t586 = Ifges(4,4) * t438;
t584 = Ifges(6,4) * t432;
t583 = Ifges(6,4) * t437;
t575 = t232 * mrSges(4,3);
t574 = t690 * Ifges(5,6);
t573 = t463 * Ifges(5,5);
t572 = t323 * Ifges(4,6);
t570 = t324 * Ifges(4,5);
t569 = t406 * Ifges(3,5);
t568 = t406 * Ifges(3,6);
t564 = t118 * t438;
t563 = t119 * t433;
t553 = t427 * t435;
t203 = Ifges(6,4) * t206;
t103 = t207 * Ifges(6,1) + t241 * Ifges(6,5) + t203;
t545 = t437 * t103;
t351 = qJD(2) * t459;
t353 = t370 * qJD(2);
t183 = -qJD(3) * t255 + t438 * t351 - t353 * t433;
t534 = qJD(2) * t427;
t503 = t439 * t534;
t296 = qJD(3) * t363 + t438 * t503;
t504 = t434 * t534;
t126 = pkin(3) * t504 - qJ(4) * t296 - qJD(4) * t364 + t183;
t182 = -t337 * t533 + t338 * t532 + t433 * t351 + t438 * t353;
t295 = -qJD(3) * t364 - t433 * t503;
t138 = qJ(4) * t295 + qJD(4) * t363 + t182;
t73 = t426 * t126 + t428 * t138;
t286 = t366 * t420 - t419 * t509;
t365 = -t429 * t507 + t547;
t543 = (-t286 * t421 + t365 * t422) * mrSges(7,1) + (-t286 * t422 - t365 * t421) * mrSges(7,2);
t290 = t368 * t420 + t419 * t553;
t367 = t429 * t546 + t508;
t221 = -t290 * t421 + t367 * t422;
t222 = t290 * t422 + t367 * t421;
t542 = t221 * mrSges(7,1) - t222 * mrSges(7,2);
t335 = t419 * t429 + t420 * t554;
t541 = (-t335 * t421 - t422 * t551) * mrSges(7,1) + (-t335 * t422 + t421 * t551) * mrSges(7,2);
t538 = t607 * pkin(1) + pkin(8) * t553;
t517 = t432 * t551;
t515 = t433 * t553;
t513 = t437 * t551;
t510 = Ifges(3,5) * t356 - Ifges(3,6) * t355 + Ifges(3,3) * t405;
t500 = t545 / 0.2e1;
t497 = -t531 / 0.2e1;
t496 = -t435 * pkin(1) + pkin(8) * t509;
t85 = -t166 * mrSges(5,1) + t167 * mrSges(5,2);
t75 = -t125 * t432 + t437 * t168;
t72 = t126 * t428 - t426 * t138;
t130 = t201 * t428 - t426 * t215;
t285 = -t366 * t419 - t420 * t509;
t396 = t433 * t509;
t491 = -t366 * t438 + t396;
t306 = -t428 * t392 - t393 * t426;
t489 = mrSges(3,3) * t505;
t124 = pkin(4) * t551 - t130;
t481 = -mrSges(3,2) + t709;
t480 = mrSges(4,1) * t363 - mrSges(4,2) * t364;
t474 = Ifges(4,1) * t438 - t587;
t473 = Ifges(6,1) * t437 - t584;
t472 = Ifges(3,2) * t439 + t589;
t471 = -Ifges(4,2) * t433 + t586;
t470 = -Ifges(6,2) * t432 + t583;
t469 = Ifges(4,5) * t438 - Ifges(4,6) * t433;
t468 = Ifges(6,5) * t437 - Ifges(6,6) * t432;
t458 = -t270 * t437 + t517;
t53 = pkin(5) * t269 + pkin(11) * t458 + t75;
t239 = -t270 * t432 - t513;
t65 = pkin(11) * t239 + t76;
t27 = -t431 * t65 + t436 * t53;
t28 = t431 * t53 + t436 * t65;
t466 = pkin(3) * t515 - t367 * t430 + t368 * t418 + t538;
t155 = -mrSges(6,2) * t241 + t591;
t156 = mrSges(6,1) * t241 - t590;
t465 = t155 * t437 - t156 * t432;
t175 = t239 * t436 + t431 * t458;
t176 = t239 * t431 - t436 * t458;
t235 = -t290 * t432 + t367 * t437;
t41 = -pkin(4) * t342 - t43;
t460 = t677 + t8;
t453 = t104 * t476;
t211 = -t428 * t295 + t296 * t426;
t212 = t295 * t426 + t296 * t428;
t354 = t371 * qJD(2);
t250 = -pkin(3) * t295 + t354;
t108 = pkin(4) * t211 - pkin(10) * t212 + t250;
t70 = pkin(10) * t504 + t73;
t22 = t432 * t108 - t125 * t531 + t168 * t530 + t437 * t70;
t444 = -g(1) * t367 - g(2) * t365 + g(3) * t551;
t69 = -pkin(4) * t504 - t72;
t23 = -qJD(5) * t76 + t437 * t108 - t432 * t70;
t416 = -pkin(4) - t601;
t401 = Ifges(3,4) * t505;
t384 = -t417 - t601;
t369 = (-mrSges(3,1) * t439 + mrSges(3,2) * t434) * t427;
t348 = -mrSges(3,2) * t406 + t489;
t304 = Ifges(3,1) * t506 + t401 + t569;
t303 = t472 * t537 + t568;
t299 = t368 * t438 + t515;
t289 = t368 * t419 - t420 * t553;
t284 = mrSges(4,1) * t388 - mrSges(4,3) * t324;
t283 = -mrSges(4,2) * t388 + t593;
t271 = t378 * t375;
t260 = pkin(5) * t556 + t306;
t236 = t290 * t437 + t367 * t432;
t224 = t388 * Ifges(4,3) + t570 + t572;
t213 = -mrSges(5,2) * t388 + mrSges(5,3) * t690;
t199 = -mrSges(4,2) * t342 + mrSges(4,3) * t234;
t198 = mrSges(4,1) * t342 - mrSges(4,3) * t233;
t178 = -mrSges(5,1) * t690 + mrSges(5,2) * t463;
t172 = -mrSges(4,1) * t234 + mrSges(4,2) * t233;
t169 = t388 * Ifges(5,3) + t573 + t574;
t144 = qJD(5) * t458 - t212 * t432 + t437 * t504;
t143 = qJD(5) * t239 + t212 * t437 + t432 * t504;
t137 = mrSges(5,1) * t342 - mrSges(5,3) * t167;
t136 = -mrSges(5,2) * t342 + mrSges(5,3) * t166;
t98 = mrSges(7,1) * t238 - mrSges(7,3) * t135;
t97 = -mrSges(7,2) * t238 + mrSges(7,3) * t493;
t95 = -pkin(5) * t239 + t124;
t74 = -mrSges(7,1) * t493 + mrSges(7,2) * t135;
t64 = -mrSges(6,2) * t165 + mrSges(6,3) * t94;
t63 = mrSges(6,1) * t165 - mrSges(6,3) * t93;
t52 = -qJD(6) * t176 - t143 * t431 + t144 * t436;
t51 = qJD(6) * t175 + t143 * t436 + t144 * t431;
t47 = -mrSges(6,1) * t94 + mrSges(6,2) * t93;
t46 = -pkin(5) * t144 + t69;
t26 = -mrSges(7,2) * t159 + mrSges(7,3) * t34;
t25 = mrSges(7,1) * t159 - mrSges(7,3) * t33;
t24 = -pkin(5) * t94 + t41;
t19 = t436 * t48 - t567;
t18 = -t431 * t48 - t566;
t15 = pkin(11) * t144 + t22;
t14 = pkin(5) * t211 - pkin(11) * t143 + t23;
t13 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t5 = -qJD(6) * t28 + t14 * t436 - t15 * t431;
t4 = qJD(6) * t27 + t14 * t431 + t15 * t436;
t1 = [(t11 * t239 + t12 * t458 - t143 * t55 + t144 * t56) * mrSges(6,3) + (-t16 * t51 + t17 * t52 + t175 * t2 - t176 * t3) * mrSges(7,3) + (-Ifges(6,4) * t458 + Ifges(6,2) * t239) * t653 + (-Ifges(6,5) * t458 + Ifges(6,6) * t239) * t642 + (Ifges(6,5) * t143 + Ifges(6,6) * t144) * t624 + (Ifges(4,1) * t364 + Ifges(4,4) * t363) * t629 + (Ifges(7,5) * t176 + Ifges(7,6) * t175) * t643 + (Ifges(7,5) * t51 + Ifges(7,6) * t52) * t626 + (-Ifges(6,1) * t458 + Ifges(6,4) * t239) * t654 + (Ifges(6,4) * t143 + Ifges(6,2) * t144) * t635 + (Ifges(7,4) * t176 + Ifges(7,2) * t175) * t663 + (Ifges(7,4) * t51 + Ifges(7,2) * t52) * t648 + (t720 / 0.2e1 + t746) * t211 + (Ifges(4,5) * t364 + Ifges(4,6) * t363) * t613 - t429 * t721 + (Ifges(6,1) * t143 + Ifges(6,4) * t144) * t633 + (Ifges(4,4) * t364 + Ifges(4,2) * t363) * t628 + t41 * (-mrSges(6,1) * t239 - mrSges(6,2) * t458) - t458 * t661 + t278 * (mrSges(3,1) * t429 - mrSges(3,3) * t554) + t51 * t657 + t52 * t659 + t239 * t662 + (t727 + t728) * t212 + (-t717 / 0.2e1 + mrSges(3,3) * t277 - Ifges(4,5) * t629 - Ifges(5,5) * t640 - Ifges(4,6) * t628 - Ifges(5,6) * t641 - t744 * t613 + t714) * t551 + (t118 * t363 - t119 * t364 - t231 * t296 + t232 * t295) * mrSges(4,3) + t739 * t269 + t28 * t26 - (t371 * mrSges(3,3) + Ifges(3,6) * t608 + (t472 + t652) * t427) * t355 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t669 + t277 * t371 + t278 * t370 - t349 * t354 + t352 * t353) + (Ifges(7,1) * t176 + Ifges(7,4) * t175) * t664 + (Ifges(7,1) * t51 + Ifges(7,4) * t52) * t646 + t740 * t270 + m(7) * (t16 * t5 + t17 * t4 + t2 * t28 + t24 * t95 + t27 * t3 + t46 * t81) + m(6) * (t104 * t69 + t11 * t76 + t12 * t75 + t124 * t41 + t22 * t56 + t23 * t55) + m(5) * (t109 * t72 + t110 * t73 + t130 * t43 + t131 * t44 + t177 * t279 + t242 * t250) + m(4) * (t118 * t255 + t119 * t254 + t182 * t232 + t183 * t231 + t252 * t336 + t308 * t354) + t176 * t665 + t175 * t668 + t27 * t25 - t252 * t480 + (-t607 * mrSges(2,1) - m(5) * t466 - t290 * mrSges(5,1) - m(6) * (pkin(4) * t290 + t466) - t236 * mrSges(6,1) - t235 * mrSges(6,2) - m(4) * (pkin(2) * t368 + t538) - t299 * mrSges(4,1) - t298 * mrSges(4,2) - m(3) * t538 - t368 * mrSges(3,1) - m(7) * (t290 * t417 + t466) - t222 * mrSges(7,1) - t221 * mrSges(7,2) + (-mrSges(3,3) * t427 + mrSges(2,2)) * t435 + (-t481 - t526) * t367 + t448 * t289) * g(2) + (t370 * mrSges(3,1) - t371 * mrSges(3,2) + Ifges(3,3) * t608 + (Ifges(3,5) * t434 + Ifges(3,6) * t439) * t427) * t405 + (-t370 * mrSges(3,3) + Ifges(3,5) * t608 + (t434 * Ifges(3,1) + t588 - t651) * t427) * t356 + ((-t349 * mrSges(3,3) + t304 / 0.2e1 + t569 / 0.2e1 + (-t651 + t588 / 0.2e1) * t537) * t439 + (t573 / 0.2e1 + t574 / 0.2e1 - t110 * mrSges(5,2) + t109 * mrSges(5,1) + t572 / 0.2e1 + t570 / 0.2e1 - t232 * mrSges(4,2) + t231 * mrSges(4,1) - t303 / 0.2e1 + t169 / 0.2e1 + t224 / 0.2e1 - t568 / 0.2e1 - t352 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t388 + (-t652 - t589 / 0.2e1 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t439) * t537) * t434) * t534 + t323 * (Ifges(4,4) * t296 + Ifges(4,2) * t295) / 0.2e1 + (-m(3) * t496 + t366 * mrSges(3,1) - mrSges(3,3) * t509 + t435 * mrSges(2,1) + t607 * mrSges(2,2) - m(4) * (-pkin(2) * t366 + t496) - t491 * mrSges(4,1) - t446 * mrSges(4,2) + t680 * t286 + (-t432 * t705 + t475 + t481 + t565) * t365 + t448 * t285 + t683 * (-pkin(3) * t396 - t365 * t430 + t366 * t418 - t496)) * g(1) + t353 * t348 + t336 * t172 + t698 * t354 + t308 * (-mrSges(4,1) * t295 + mrSges(4,2) * t296) + t182 * t283 + t183 * t284 + t279 * t85 + t255 * t199 + t250 * t178 + t254 * t198 + t73 * t213 + t72 * t214 + t24 * (-mrSges(7,1) * t175 + mrSges(7,2) * t176) + t364 * t644 + t363 * t645 + t296 * t630 + t295 * t631 + (Ifges(4,1) * t296 + Ifges(4,4) * t295) * t614 + t510 * t608 + (Ifges(4,5) * t296 + Ifges(4,6) * t295) * t609 + (-pkin(1) * t369 * t427 + Ifges(2,3)) * qJDD(1) + t46 * t74 + t75 * t63 + t76 * t64 + t81 * (-mrSges(7,1) * t52 + mrSges(7,2) * t51) + t95 * t13 + t4 * t97 + t5 * t98 + t124 * t47 + t131 * t136 + t130 * t137 + t69 * t139 + t143 * t103 / 0.2e1 + t144 * t102 / 0.2e1 + t104 * (-mrSges(6,1) * t144 + mrSges(6,2) * t143) + t22 * t155 + t23 * t156; (Ifges(6,1) * t267 + Ifges(6,4) * t266) * t634 + (-Ifges(6,1) * t455 - Ifges(6,4) * t456) * t633 + (-Ifges(7,5) * t272 - Ifges(7,6) * t271) * t643 + t24 * (mrSges(7,1) * t271 - mrSges(7,2) * t272) + (t110 * t506 + t242 * t302) * mrSges(5,2) + (-Ifges(5,1) * t302 + Ifges(5,5) * t506) * t619 + (-Ifges(5,5) * t302 + Ifges(5,3) * t506) * t610 + (-Ifges(5,4) * t302 + Ifges(5,6) * t506) * t622 - t109 * (mrSges(5,1) * t506 + mrSges(5,3) * t302) + (-Ifges(7,1) * t272 - Ifges(7,4) * t271) * t664 + (-Ifges(7,4) * t272 - Ifges(7,2) * t271) * t663 + (t323 * t471 + t324 * t474 + t388 * t469) * qJD(3) / 0.2e1 + (-mrSges(4,3) * t532 - (mrSges(4,1) * t434 - mrSges(4,3) * t544) * t537) * t231 + (-t434 * (Ifges(3,1) * t439 - t589) / 0.2e1 + pkin(1) * (mrSges(3,1) * t434 + mrSges(3,2) * t439)) * qJD(1) ^ 2 * t669 + t720 * (-t301 / 0.2e1 + t361 / 0.2e1) + (-t284 * t532 - t283 * t533 + t438 * t199 + m(4) * (t564 - t563 + (-t231 * t438 - t232 * t433) * qJD(3)) - t433 * t198) * pkin(9) + (-Ifges(6,5) * t455 - Ifges(6,6) * t456) * t624 + (Ifges(6,5) * t267 + Ifges(6,6) * t266) * t625 + (t369 - t683 * t418 * t551 + (t745 * t439 + (t430 * t683 + t715) * t434) * t427) * g(3) + t716 * t301 + (t489 - t348) * t349 + t746 * t361 + (Ifges(7,1) * t150 + Ifges(7,4) * t151) * t646 + (Ifges(7,1) * t185 + Ifges(7,4) * t184) * t647 + (-Ifges(6,4) * t455 - Ifges(6,2) * t456) * t635 + (Ifges(6,4) * t267 + Ifges(6,2) * t266) * t636 + (Ifges(7,4) * t150 + Ifges(7,2) * t151) * t648 + (Ifges(7,4) * t185 + Ifges(7,2) * t184) * t649 + (t47 - t137) * t306 + (-pkin(2) * t252 - t231 * t261 - t232 * t262) * m(4) - ((t224 + t169) * t434 + (-Ifges(3,2) * t506 + t438 * t226 + t304 + t401) * t439 + t388 * (Ifges(4,3) * t434 + t439 * t469) + t324 * (Ifges(4,5) * t434 + t439 * t474) + t323 * (Ifges(4,6) * t434 + t439 * t471) + t406 * (Ifges(3,5) * t439 - Ifges(3,6) * t434)) * t537 / 0.2e1 + t150 * t657 + t185 * t658 + t151 * t659 + t184 * t660 + t555 * t661 + (Ifges(7,5) * t150 + Ifges(7,6) * t151) * t626 + (Ifges(7,5) * t185 + Ifges(7,6) * t184) * t627 + t510 - t721 - t38 * t556 / 0.2e1 + t303 * t506 / 0.2e1 + t739 * t374 - mrSges(4,3) * t563 + (t103 * t497 + t41 * t476 + t468 * t642 + t470 * t653 + t473 * t654 + t740) * t375 - (t500 + t673 + t727) * t362 + (-t225 / 0.2e1 - t575) * t533 - t272 * t665 - t271 * t668 - t232 * (-mrSges(4,3) * t433 * t439 - mrSges(4,2) * t434) * t537 + t252 * t479 - t418 * t85 + t710 * t98 + t711 * t97 + (t100 * t2 + t16 * t710 + t17 * t711 + t24 * t260 + t3 * t99 + t700 * t81) * m(7) + t703 * t155 + t704 * t156 + (t104 * t696 + t11 * t205 + t12 * t204 + t306 * t41 + t55 * t704 + t56 * t703) * m(6) + (-m(4) * t308 + t490 - t698) * t352 + t700 * t74 + (-mrSges(7,1) * t701 + mrSges(7,2) * t702) * t81 + (-t16 * t702 + t17 * t701 - t2 * t271 + t272 * t3) * mrSges(7,3) + (mrSges(6,1) * t694 - mrSges(6,2) * t693) * t104 + (-t11 * t556 - t12 * t555 + t55 * t693 - t56 * t694) * mrSges(6,3) + t695 * t213 + t696 * t139 + t697 * t214 + (t109 * t697 + t110 * t695 - t177 * t418 + t242 * t686 - t306 * t43 + t307 * t44) * m(5) + t307 * t136 + t686 * t178 + (-t683 * (-t367 * t418 - t368 * t430) + t675 * t368 + t674 * t367) * g(1) + (-t683 * (-t365 * t418 - t366 * t430) + t675 * t366 + t674 * t365) * g(2) - t262 * t283 - t261 * t284 + t278 * mrSges(3,1) - t267 * t103 / 0.2e1 + t260 * t13 + t204 * t63 + t205 * t64 + t433 * t644 + t438 * t645 + (Ifges(4,2) * t438 + t587) * t628 + (Ifges(4,1) * t433 + t586) * t629 + t532 * t630 + t485 * t631 + (Ifges(4,5) * t433 + Ifges(4,6) * t438) * t613 + mrSges(4,3) * t564 + (t266 + t502) * t724 + t302 * t637 + t99 * t25 + t100 * t26 + t388 * t308 * (mrSges(4,1) * t433 + mrSges(4,2) * t438) - pkin(2) * t172; (-Ifges(7,5) * t174 - Ifges(7,6) * t173) * t627 + (-Ifges(7,4) * t174 - Ifges(7,2) * t173) * t649 + (-Ifges(7,1) * t174 - Ifges(7,4) * t173) * t647 - t81 * (mrSges(7,1) * t173 - mrSges(7,2) * t174) + (t206 * t470 + t207 * t473 + t241 * t468) * qJD(5) / 0.2e1 - (-Ifges(4,2) * t324 + t226 + t315) * t323 / 0.2e1 + (t500 + t453) * qJD(5) + (-t104 * t113 + t41 * t416 - t55 * t67 - t56 * t68) * m(6) + (-t283 + t593) * t231 + t720 * t619 + (mrSges(4,2) * t299 + t680 * t289 + t672 * t290 + (-mrSges(4,1) - t722) * t298) * g(1) + t717 + t716 * t463 - m(5) * (t110 * t114 + t242 * t603) + (Ifges(5,1) * t619 + Ifges(5,4) * t622 + Ifges(5,5) * t610 + t468 * t625 + t470 * t636 + t473 * t634 - t453 - t673) * t690 + (Ifges(7,1) * t378 - Ifges(7,4) * t462) * t664 + t24 * (mrSges(7,1) * t462 + mrSges(7,2) * t378) + (Ifges(7,5) * t378 - Ifges(7,6) * t462) * t643 + (-t16 * t174 + t17 * t173 - t2 * t462 - t3 * t378) * mrSges(7,3) + (Ifges(7,4) * t378 - Ifges(7,2) * t462) * t663 - t462 * t668 + (Ifges(6,2) * t437 + t584) * t653 + (Ifges(6,1) * t432 + t583) * t654 - t174 * t658 - t173 * t660 + t432 * t661 + t437 * t662 + t378 * t665 + (-g(1) * t290 - g(2) * t286 - g(3) * t335 - t731 * t56 + (-t530 + t560) * t55 + t685) * mrSges(6,3) - t714 - t178 * t603 - t324 * (Ifges(4,1) * t323 - t571) / 0.2e1 + (t426 * t44 + t428 * t43) * t667 + t41 * t477 + (-mrSges(4,2) * t491 - t680 * t285 + t672 * t286 + (-pkin(3) * t712 + mrSges(4,1) + t667) * t446) * g(2) + t416 * t47 + t384 * t13 + t706 * t74 + t707 * t98 + t708 * t97 + (t16 * t707 + t17 * t708 + t2 * t281 + t24 * t384 + t280 * t3 + t706 * t81) * m(7) + (m(5) * t109 + t699) * t113 - t308 * (mrSges(4,1) * t324 + mrSges(4,2) * t323) + (m(6) * ((-t432 * t56 - t437 * t55) * qJD(5) + t685) - t155 * t531 - t156 * t530 + t437 * t64 - t432 * t63) * t415 - (Ifges(7,4) * t646 + Ifges(7,2) * t648 + Ifges(7,6) * t626 + t659 + t679) * t294 + t280 * t25 + t281 * t26 + t232 * t284 - (Ifges(7,1) * t646 + Ifges(7,4) * t648 + Ifges(7,5) * t626 + t657 + t678) * t293 - t114 * t213 + (Ifges(6,5) * t432 + Ifges(6,6) * t437) * t642 + t225 * t614 + t137 * t601 + t136 * t602 + (Ifges(4,5) * t323 - Ifges(4,6) * t324) * t610 + t324 * t575 + (-t480 + t672 * t335 - t680 * (-t419 * t554 + t420 * t429) - t363 * t722) * g(3) + t102 * t497 + (t545 + t171) * t622 - t68 * t155 - t67 * t156; -t462 * t25 + t378 * t26 + t432 * t64 + t437 * t63 + t691 * t98 + t692 * t97 + t465 * qJD(5) + (-t213 - t465) * t690 + (-t74 + t699) * t463 + t85 + (t16 * t691 + t17 * t692 + t2 * t378 - t3 * t462 - t463 * t81 + t444) * m(7) + (-t104 * t463 + t11 * t432 + t12 * t437 + t444 + t241 * (-t432 * t55 + t437 * t56)) * m(6) + (t109 * t463 - t110 * t690 + t177 + t444) * m(5); (mrSges(6,2) * t236 + t235 * t705 - t542) * g(1) + (-t543 - (-t286 * t437 - t365 * t432) * mrSges(6,2) + t705 * (-t286 * t432 + t365 * t437)) * g(2) - (Ifges(7,4) * t647 + Ifges(7,2) * t649 + Ifges(7,6) * t627 + t660 - t679) * t135 + (Ifges(7,1) * t647 + Ifges(7,4) * t649 + Ifges(7,5) * t627 + t658 - t678) * t493 + (-t541 - (-t335 * t437 + t517) * mrSges(6,2) + t705 * (-t335 * t432 - t513)) * g(3) + t460 + t37 - t74 * t600 - m(7) * (t16 * t18 + t17 * t19 + t600 * t81) + (t2 * t431 + t3 * t436 + (-t16 * t431 + t17 * t436) * qJD(6)) * t666 + t676 - t104 * (mrSges(6,1) * t207 + mrSges(6,2) * t206) + (Ifges(6,1) * t206 - t585) * t634 + t102 * t633 + (Ifges(6,5) * t206 - Ifges(6,6) * t207) * t625 + (t591 - t155) * t55 - t19 * t97 - t18 * t98 + (t590 + t156) * t56 + (-Ifges(6,2) * t207 + t103 + t203) * t636 + ((-t431 * t98 + t436 * t97) * qJD(6) + t25 * t436 + t26 * t431) * pkin(5); -t81 * (mrSges(7,1) * t135 + mrSges(7,2) * t493) + (Ifges(7,1) * t493 - t582) * t647 + t61 * t646 + (Ifges(7,5) * t493 - Ifges(7,6) * t135) * t627 - t16 * t97 + t17 * t98 - g(1) * t542 - g(2) * t543 - g(3) * t541 + (t135 * t17 + t16 * t493) * mrSges(7,3) + t460 + (-Ifges(7,2) * t135 + t127 + t62) * t649;];
tau  = t1;
