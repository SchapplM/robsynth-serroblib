% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:48:10
% EndTime: 2019-03-09 05:49:12
% DurationCPUTime: 38.20s
% Computational Cost: add. (24503->955), mult. (77974->1234), div. (0->0), fcn. (67310->14), ass. (0->454)
t329 = sin(qJ(4));
t474 = qJD(4) * t329;
t326 = sin(pkin(7));
t330 = sin(qJ(3));
t327 = sin(pkin(6));
t495 = cos(pkin(12));
t496 = cos(pkin(7));
t415 = t496 * t495;
t539 = cos(qJ(3));
t379 = t539 * t415;
t372 = t327 * t379;
t497 = cos(pkin(6));
t421 = t497 * t539;
t386 = qJD(1) * t421;
t325 = sin(pkin(12));
t478 = qJD(1) * t327;
t455 = t325 * t478;
t236 = qJD(1) * t372 + t326 * t386 - t330 * t455;
t489 = t236 * t329;
t609 = t474 - t489;
t437 = t330 * t496;
t359 = t327 * (-t325 * t437 + t495 * t539);
t265 = qJD(1) * t359;
t450 = qJD(3) * t539;
t423 = t326 * t450;
t669 = t423 - t265;
t229 = qJD(4) - t236;
t670 = Ifges(5,4) + Ifges(6,6);
t416 = t497 * t495;
t538 = sin(qJ(1));
t540 = cos(qJ(1));
t371 = t538 * t325 - t416 * t540;
t458 = t327 * t540;
t668 = t326 * t458 + t371 * t496;
t320 = pkin(1) * t416;
t313 = qJD(1) * t320;
t456 = t497 * pkin(2);
t486 = t325 * t327;
t360 = t456 + (-pkin(9) * t496 - qJ(2)) * t486;
t230 = qJD(1) * t360 + t313;
t218 = t230 * t437;
t440 = t325 * t497;
t319 = pkin(1) * t440;
t438 = t327 * t495;
t419 = qJD(1) * t438;
t277 = qJ(2) * t419 + qJD(1) * t319;
t435 = t497 * t326;
t357 = (t327 * t415 + t435) * pkin(9);
t224 = qJD(1) * t357 + t277;
t537 = pkin(9) * t325;
t270 = (-pkin(2) * t495 - t326 * t537 - pkin(1)) * t327;
t257 = qJD(1) * t270 + qJD(2);
t484 = t326 * t330;
t134 = t539 * t224 + t257 * t484 + t218;
t667 = pkin(4) * t609 - qJD(5) * t329 - t134;
t283 = t440 * t540 + t495 * t538;
t207 = -t283 * t539 + t330 * t668;
t332 = cos(qJ(4));
t439 = t327 * t496;
t346 = t371 * t326 - t540 * t439;
t159 = t207 * t332 - t329 * t346;
t158 = t207 * t329 + t332 * t346;
t367 = t325 * t539 + t330 * t415;
t249 = t327 * t367 + t330 * t435;
t239 = t249 * qJD(1);
t282 = -t326 * t438 + t496 * t497;
t350 = qJDD(1) * t282 + qJDD(3);
t485 = t325 * t330;
t366 = t379 - t485;
t469 = qJD(1) * qJD(3);
t336 = (qJDD(1) * t367 + t366 * t469) * t327;
t432 = qJDD(1) * t497;
t354 = qJD(3) * t386 + t330 * t432;
t334 = t354 * t326 + t336;
t351 = qJD(1) * t282 + qJD(3);
t658 = qJD(4) * t351 + t334;
t104 = t239 * t474 - t329 * t350 - t332 * t658;
t565 = -t104 / 0.2e1;
t473 = qJD(4) * t332;
t105 = t239 * t473 + t329 * t658 - t332 * t350;
t563 = -t105 / 0.2e1;
t176 = t327 * (qJDD(1) * t366 - t367 * t469) - (t330 * t469 * t497 - qJDD(1) * t421) * t326;
t174 = qJDD(4) - t176;
t556 = t174 / 0.2e1;
t644 = Ifges(6,4) - Ifges(5,5);
t643 = Ifges(5,6) - Ifges(6,5);
t642 = Ifges(5,3) + Ifges(6,1);
t666 = t667 + t229 * (pkin(11) * t329 - qJ(5) * t332);
t420 = t496 * t539;
t389 = t230 * t420;
t459 = t326 * t539;
t133 = -t330 * t224 + t257 * t459 + t389;
t129 = t329 * t133;
t177 = pkin(3) * t239 - pkin(10) * t236;
t568 = pkin(5) + pkin(10);
t309 = t568 * t332;
t569 = pkin(4) + pkin(11);
t665 = qJD(4) * t309 - t129 - (pkin(5) * t236 - t177) * t332 + t569 * t239;
t290 = t329 * t496 + t332 * t484;
t425 = t326 * t455;
t664 = qJD(4) * t290 + t669 * t329 + t332 * t425;
t358 = t327 * (t325 * t420 + t330 * t495);
t264 = qJD(1) * t358;
t477 = qJD(3) * t330;
t453 = t326 * t477;
t663 = t453 - t264;
t488 = t236 * t332;
t610 = t473 - t488;
t188 = t239 * t329 - t332 * t351;
t328 = sin(qJ(6));
t331 = cos(qJ(6));
t136 = t188 * t331 - t229 * t328;
t46 = qJD(6) * t136 + t105 * t328 + t174 * t331;
t97 = qJDD(6) - t104;
t28 = mrSges(7,1) * t97 - mrSges(7,3) * t46;
t137 = t188 * t328 + t229 * t331;
t47 = -qJD(6) * t137 + t105 * t331 - t174 * t328;
t29 = -mrSges(7,2) * t97 + mrSges(7,3) * t47;
t189 = t332 * t239 + t329 * t351;
t187 = qJD(6) + t189;
t90 = -mrSges(7,2) * t187 + mrSges(7,3) * t136;
t91 = mrSges(7,1) * t187 - mrSges(7,3) * t137;
t391 = -t328 * t91 + t331 * t90;
t662 = t391 * qJD(6) + t331 * t28 + t328 * t29;
t564 = t104 / 0.2e1;
t562 = t105 / 0.2e1;
t661 = -t174 / 0.2e1;
t660 = -mrSges(6,1) - mrSges(5,3);
t186 = Ifges(5,4) * t188;
t634 = t189 * Ifges(5,1) + t229 * Ifges(5,5) + t137 * Ifges(7,5) + t136 * Ifges(7,6) + t187 * Ifges(7,3) - t186;
t184 = -t326 * t230 + t257 * t496;
t364 = -t236 * pkin(3) - t239 * pkin(10) + t184;
t292 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t327;
t261 = qJDD(1) * t319 + t495 * t292;
t213 = qJDD(1) * t357 + t261;
t260 = qJDD(1) * t320 - t325 * t292;
t214 = (-t439 * t537 + t456) * qJDD(1) + t260;
t253 = qJDD(1) * t270 + qJDD(2);
t69 = qJD(3) * t389 + t539 * t213 + t214 * t437 - t224 * t477 + t253 * t484 + t257 * t423;
t659 = pkin(10) * t350 + qJD(4) * t364 + t69;
t289 = t329 * t484 - t332 * t496;
t617 = qJD(4) * t289 + t329 * t425 - t332 * t669;
t363 = t325 * t540 + t416 * t538;
t340 = t363 * t326 + t538 * t439;
t598 = m(7) * pkin(11) + mrSges(5,1) - mrSges(6,2);
t657 = pkin(4) * (m(7) + m(6)) + t598;
t113 = -pkin(3) * t351 - t133;
t68 = t188 * pkin(4) - t189 * qJ(5) + t113;
t656 = t113 * mrSges(5,2) - t68 * mrSges(6,3);
t114 = pkin(10) * t351 + t134;
t55 = t114 * t329 - t332 * t364;
t385 = pkin(5) * t189 + t55;
t654 = qJD(5) + t385;
t570 = t97 / 0.2e1;
t573 = t47 / 0.2e1;
t574 = t46 / 0.2e1;
t70 = -qJD(3) * t218 - t330 * t213 + t214 * t420 - t224 * t450 + t253 * t459 - t257 * t453;
t66 = -pkin(3) * t350 - t70;
t18 = t105 * pkin(4) + t104 * qJ(5) - t189 * qJD(5) + t66;
t12 = t105 * pkin(11) + t18;
t34 = -t229 * t569 + t654;
t51 = t188 * pkin(11) + t68;
t15 = -t328 * t51 + t331 * t34;
t243 = t496 * t253;
t84 = -t176 * pkin(3) + t243 + (-pkin(10) * t354 - t214) * t326 - pkin(10) * t336;
t14 = -t114 * t473 - t329 * t659 + t332 * t84;
t381 = qJDD(5) - t14;
t5 = -pkin(5) * t104 - t174 * t569 + t381;
t1 = qJD(6) * t15 + t12 * t331 + t328 * t5;
t16 = t328 * t34 + t331 * t51;
t2 = -qJD(6) * t16 - t12 * t328 + t331 * t5;
t597 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t9 = Ifges(7,5) * t46 + Ifges(7,6) * t47 + Ifges(7,3) * t97;
t653 = 0.2e1 * Ifges(5,1) * t565 + Ifges(6,4) * t661 + Ifges(7,5) * t574 + Ifges(7,6) * t573 + Ifges(7,3) * t570 + t597 + (-t564 + t565) * Ifges(6,2) + t9 / 0.2e1 + t670 * t563 + (-t644 + Ifges(5,5)) * t556;
t204 = t283 * t330 + t539 * t668;
t652 = Ifges(5,6) * t661 + 0.2e1 * Ifges(6,3) * t562 + t670 * t564 + (t562 - t563) * Ifges(5,2) + (-t643 + Ifges(6,5)) * t556;
t525 = mrSges(6,1) * t189;
t139 = mrSges(6,2) * t229 + t525;
t523 = mrSges(5,3) * t189;
t141 = mrSges(5,1) * t229 - t523;
t52 = -pkin(4) * t229 + qJD(5) + t55;
t651 = -m(6) * t52 - t139 + t141;
t650 = t249 / 0.2e1;
t648 = m(7) * t568;
t646 = t69 * mrSges(4,2);
t645 = mrSges(5,2) - mrSges(6,3);
t73 = mrSges(6,1) * t105 - mrSges(6,3) * t174;
t76 = -mrSges(5,2) * t174 - mrSges(5,3) * t105;
t641 = -t73 + t76;
t74 = -t104 * mrSges(6,1) + t174 * mrSges(6,2);
t75 = mrSges(5,1) * t174 + mrSges(5,3) * t104;
t640 = t74 - t75;
t228 = Ifges(4,4) * t236;
t639 = Ifges(4,2) * t236;
t637 = t351 * Ifges(4,5);
t636 = t351 * Ifges(4,6);
t635 = -t188 * t643 - t189 * t644 + t229 * t642;
t526 = mrSges(6,1) * t188;
t138 = -mrSges(6,3) * t229 + t526;
t85 = -mrSges(7,1) * t136 + mrSges(7,2) * t137;
t633 = t138 - t85;
t493 = qJ(5) * t329;
t444 = -pkin(3) - t493;
t291 = -t332 * t569 + t444;
t308 = t568 * t329;
t251 = -t291 * t328 + t308 * t331;
t632 = qJD(6) * t251 + t328 * t665 + t331 * t666;
t252 = t291 * t331 + t308 * t328;
t631 = -qJD(6) * t252 - t328 * t666 + t331 * t665;
t630 = -qJ(5) * t610 + t667;
t87 = t332 * t133 + t329 * t177;
t77 = -qJ(5) * t239 - t87;
t629 = pkin(5) * t489 - t568 * t474 + t77;
t524 = mrSges(5,3) * t188;
t140 = -mrSges(5,2) * t229 - t524;
t626 = -t140 + t138;
t481 = t329 * t331;
t165 = t236 * t481 - t239 * t328;
t471 = qJD(6) * t332;
t452 = t328 * t471;
t625 = -t331 * t474 + t165 - t452;
t483 = t328 * t329;
t166 = t236 * t483 + t239 * t331;
t624 = -t328 * t474 + t331 * t471 + t166;
t375 = -t328 * t289 + t331 * t459;
t623 = qJD(6) * t375 - t328 * t663 + t331 * t664;
t258 = t331 * t289 + t328 * t459;
t622 = qJD(6) * t258 + t328 * t664 + t331 * t663;
t503 = t239 * mrSges(4,3);
t621 = mrSges(4,1) * t351 - mrSges(5,1) * t188 - mrSges(5,2) * t189 - t503;
t491 = t204 * t332;
t620 = -pkin(4) * t491 - t204 * t493;
t284 = -t440 * t538 + t495 * t540;
t457 = t327 * t538;
t615 = -t326 * t457 + t363 * t496;
t208 = t284 * t330 + t539 * t615;
t490 = t208 * t332;
t619 = -pkin(4) * t490 - t208 * t493;
t286 = qJ(2) * t438 + t319;
t242 = t357 + t286;
t250 = t320 + t360;
t142 = -t330 * t242 + t250 * t420 + t270 * t459;
t390 = t326 * t421;
t248 = t327 * t485 - t372 - t390;
t487 = t248 * t332;
t618 = -pkin(4) * t487 - t248 * t493;
t614 = t113 * (mrSges(5,1) * t329 + mrSges(5,2) * t332) + t68 * (-mrSges(6,2) * t329 - mrSges(6,3) * t332);
t613 = -t329 * t643 - t332 * t644;
t611 = t104 * t644 - t105 * t643 + t174 * t642;
t13 = -t114 * t474 + t329 * t84 + t332 * t659;
t607 = t13 * t332 - t14 * t329;
t606 = Ifges(4,4) * t249 + Ifges(4,6) * t282;
t605 = Ifges(4,1) * t249 + Ifges(4,5) * t282;
t7 = -qJ(5) * t174 - qJD(5) * t229 - t13;
t8 = -pkin(4) * t174 + t381;
t604 = t329 * t8 - t332 * t7;
t600 = -g(1) * t457 + g(2) * t458 - g(3) * t497;
t56 = t332 * t114 + t329 * t364;
t53 = -t229 * qJ(5) - t56;
t534 = t188 * pkin(5);
t39 = -t53 - t534;
t409 = t331 * mrSges(7,1) - t328 * mrSges(7,2);
t135 = Ifges(7,4) * t136;
t64 = Ifges(7,1) * t137 + Ifges(7,5) * t187 + t135;
t501 = t328 * t64;
t599 = t39 * t409 - t501 / 0.2e1;
t596 = -m(5) * t113 + t621;
t595 = t15 * mrSges(7,1) - t16 * mrSges(7,2);
t594 = -mrSges(7,3) - t598;
t593 = -t113 * mrSges(5,1) + t68 * mrSges(6,2);
t408 = mrSges(7,1) * t328 + mrSges(7,2) * t331;
t377 = m(7) * qJ(5) + t408;
t592 = -m(6) * qJ(5) - t377 + t645;
t591 = t184 * mrSges(4,2) + t637 / 0.2e1;
t185 = Ifges(6,6) * t188;
t100 = t229 * Ifges(6,4) - t189 * Ifges(6,2) + t185;
t589 = -t100 / 0.2e1 + t595;
t209 = t284 * t539 - t330 * t615;
t160 = t209 * t329 - t332 * t340;
t202 = t249 * t329 - t282 * t332;
t588 = g(1) * t160 - g(2) * t158 + g(3) * t202;
t587 = -g(1) * t209 + g(2) * t207 - g(3) * t249;
t406 = mrSges(6,2) * t332 - mrSges(6,3) * t329;
t411 = -mrSges(5,1) * t332 + mrSges(5,2) * t329;
t586 = mrSges(7,1) * t483 + mrSges(7,2) * t481 - t406 - t411;
t585 = -t409 - t648;
t583 = mrSges(4,1) + t586;
t582 = mrSges(4,2) + t585;
t581 = t14 * mrSges(5,1) - t13 * mrSges(5,2) + t8 * mrSges(6,2) - t7 * mrSges(6,3);
t555 = -t187 / 0.2e1;
t559 = -t137 / 0.2e1;
t561 = -t136 / 0.2e1;
t580 = Ifges(7,5) * t559 + Ifges(7,6) * t561 + Ifges(7,3) * t555 - t595;
t579 = t53 * mrSges(6,3) + t55 * mrSges(5,1) + t56 * mrSges(5,2) - t184 * mrSges(4,1) + t636 / 0.2e1 - t52 * mrSges(6,2);
t10 = Ifges(7,4) * t46 + Ifges(7,2) * t47 + Ifges(7,6) * t97;
t577 = -t10 / 0.2e1;
t11 = Ifges(7,1) * t46 + Ifges(7,4) * t47 + Ifges(7,5) * t97;
t576 = t11 / 0.2e1;
t516 = Ifges(7,4) * t137;
t63 = Ifges(7,2) * t136 + Ifges(7,6) * t187 + t516;
t572 = -t63 / 0.2e1;
t571 = t63 / 0.2e1;
t506 = t189 * Ifges(5,4);
t101 = -t188 * Ifges(5,2) + t229 * Ifges(5,6) + t506;
t566 = -t101 / 0.2e1;
t560 = t136 / 0.2e1;
t558 = t137 / 0.2e1;
t554 = t187 / 0.2e1;
t553 = -t188 / 0.2e1;
t552 = t188 / 0.2e1;
t551 = -t189 / 0.2e1;
t550 = t189 / 0.2e1;
t548 = -t229 / 0.2e1;
t547 = t229 / 0.2e1;
t546 = -t236 / 0.2e1;
t544 = -t239 / 0.2e1;
t543 = t239 / 0.2e1;
t533 = t2 * t331;
t532 = t202 * pkin(11);
t531 = t204 * pkin(10);
t530 = t208 * pkin(10);
t522 = mrSges(7,3) * t328;
t520 = Ifges(4,4) * t239;
t518 = Ifges(5,4) * t329;
t517 = Ifges(5,4) * t332;
t515 = Ifges(7,4) * t328;
t514 = Ifges(7,4) * t331;
t511 = Ifges(6,6) * t329;
t510 = Ifges(6,6) * t332;
t507 = t15 * t328;
t505 = t189 * Ifges(6,6);
t504 = t236 * mrSges(4,3);
t494 = qJ(5) * t188;
t492 = t189 * t331;
t482 = t328 * t332;
t480 = t331 * t332;
t190 = -t250 * t326 + t496 * t270;
t247 = t248 * pkin(3);
t441 = pkin(10) * t249 - t247;
t125 = t190 - t441;
t227 = t539 * t242;
t434 = t496 * t250;
t143 = t270 * t484 + t330 * t434 + t227;
t132 = pkin(10) * t282 + t143;
t72 = t329 * t125 + t332 * t132;
t479 = t540 * pkin(1) + qJ(2) * t457;
t472 = qJD(6) * t331;
t468 = qJDD(1) * t327;
t465 = pkin(10) * t474;
t464 = pkin(10) * t473;
t461 = Ifges(4,5) * t334 + Ifges(4,6) * t176 + Ifges(4,3) * t350;
t454 = qJD(2) * t486;
t449 = t325 * t468;
t445 = -t472 / 0.2e1;
t198 = t204 * pkin(3);
t443 = -pkin(10) * t207 - t198;
t200 = t208 * pkin(3);
t442 = pkin(10) * t209 - t200;
t197 = t202 * pkin(4);
t203 = t249 * t332 + t282 * t329;
t433 = t203 * qJ(5) - t197;
t71 = t125 * t332 - t329 * t132;
t86 = t177 * t332 - t129;
t424 = t326 * t454;
t422 = -pkin(1) * t538 + qJ(2) * t458;
t58 = -qJ(5) * t248 - t72;
t418 = qJDD(1) * t438;
t414 = mrSges(4,1) * t248 + mrSges(4,2) * t249;
t413 = mrSges(4,1) * t282 - mrSges(4,3) * t249;
t412 = mrSges(5,1) * t202 + mrSges(5,2) * t203;
t150 = t202 * t331 - t248 * t328;
t151 = t202 * t328 + t248 * t331;
t410 = mrSges(7,1) * t150 - mrSges(7,2) * t151;
t407 = -t202 * mrSges(6,2) - t203 * mrSges(6,3);
t405 = Ifges(5,1) * t332 - t518;
t404 = Ifges(7,1) * t331 - t515;
t403 = Ifges(7,1) * t328 + t514;
t402 = -Ifges(5,2) * t329 + t517;
t400 = -Ifges(7,2) * t328 + t514;
t399 = Ifges(7,2) * t331 + t515;
t397 = Ifges(7,5) * t331 - Ifges(7,6) * t328;
t396 = Ifges(7,5) * t328 + Ifges(7,6) * t331;
t395 = -Ifges(6,2) * t332 + t511;
t394 = Ifges(6,3) * t329 - t510;
t392 = -t16 * t331 + t507;
t41 = pkin(5) * t203 - t248 * t569 - t71;
t131 = -t282 * pkin(3) - t142;
t81 = t131 - t433;
t57 = t81 + t532;
t19 = -t328 * t57 + t331 * t41;
t20 = t328 * t41 + t331 * t57;
t117 = qJD(2) * t359 + qJD(3) * t142;
t237 = (t327 * t366 + t390) * qJD(3);
t238 = t249 * qJD(3);
t163 = pkin(3) * t238 - pkin(10) * t237 + t424;
t33 = -t329 * t117 - t125 * t474 - t132 * t473 + t163 * t332;
t382 = -(-qJ(2) * t455 + t313) * t325 + t277 * t495;
t32 = t332 * t117 + t125 * t473 - t132 * t474 + t329 * t163;
t378 = -mrSges(3,1) * t418 + mrSges(3,2) * t449;
t376 = mrSges(3,1) * t497 - mrSges(3,3) * t486;
t370 = -mrSges(3,2) * t497 + mrSges(3,3) * t438;
t356 = -qJD(6) * t392 + t1 * t328 + t533;
t27 = -qJ(5) * t238 - qJD(5) * t248 - t32;
t344 = -t283 * pkin(2) - pkin(9) * t346 + t422;
t343 = t284 * pkin(2) + pkin(9) * t340 + t479;
t342 = t207 * pkin(3) + t344;
t341 = t209 * pkin(3) + t343;
t339 = qJD(1) * t237;
t338 = t159 * pkin(4) + t158 * qJ(5) + t342;
t161 = t209 * t332 + t329 * t340;
t337 = t161 * pkin(4) + t160 * qJ(5) + t341;
t118 = qJD(2) * t358 + (t227 + (t270 * t326 + t434) * t330) * qJD(3);
t149 = -t249 * t474 + (qJD(4) * t282 + t237) * t332;
t335 = -t149 * qJ(5) - t203 * qJD(5) + t118;
t314 = -pkin(1) * t468 + qJDD(2);
t301 = -pkin(4) * t332 + t444;
t288 = t370 * qJD(1);
t287 = t376 * qJD(1);
t285 = -qJ(2) * t486 + t320;
t191 = -mrSges(4,2) * t351 + t504;
t175 = -mrSges(4,1) * t236 + mrSges(4,2) * t239;
t164 = -t326 * t214 + t243;
t148 = qJD(4) * t203 + t237 * t329;
t147 = Ifges(4,1) * t239 + t228 + t637;
t146 = t520 + t636 + t639;
t145 = -mrSges(4,2) * t350 + t176 * mrSges(4,3);
t144 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t339 + qJDD(1) * t413;
t123 = -mrSges(6,2) * t188 - mrSges(6,3) * t189;
t121 = pkin(4) * t189 + t494;
t108 = t160 * t328 + t208 * t331;
t107 = t160 * t331 - t208 * t328;
t106 = -t176 * mrSges(4,1) + mrSges(4,2) * t334;
t98 = t229 * Ifges(6,5) + t188 * Ifges(6,3) - t505;
t89 = t189 * t569 + t494;
t83 = qJD(6) * t150 + t148 * t328 + t238 * t331;
t82 = -qJD(6) * t151 + t148 * t331 - t238 * t328;
t79 = -pkin(4) * t239 - t86;
t59 = -pkin(4) * t248 - t71;
t50 = mrSges(5,1) * t105 - mrSges(5,2) * t104;
t49 = -mrSges(6,2) * t105 + mrSges(6,3) * t104;
t48 = -pkin(5) * t202 - t58;
t45 = t56 - t534;
t40 = t148 * pkin(4) + t335;
t31 = t148 * t569 + t335;
t30 = -pkin(4) * t238 - t33;
t24 = t328 * t45 + t331 * t89;
t23 = -t328 * t89 + t331 * t45;
t22 = -pkin(5) * t148 - t27;
t21 = pkin(5) * t149 - t238 * t569 - t33;
t17 = -mrSges(7,1) * t47 + mrSges(7,2) * t46;
t6 = -pkin(5) * t105 - t7;
t4 = -qJD(6) * t20 + t21 * t331 - t31 * t328;
t3 = qJD(6) * t19 + t21 * t328 + t31 * t331;
t25 = [(-m(5) * (t342 - t531) - m(6) * (t338 - t531) - m(7) * t338 + mrSges(2,1) * t538 + mrSges(2,2) * t540 - m(3) * t422 + t283 * mrSges(3,1) - mrSges(3,2) * t371 - mrSges(3,3) * t458 - m(4) * t344 - t207 * mrSges(4,1) + mrSges(4,3) * t346 + (-t408 + t645) * t158 - (t582 + t660) * t204 + t594 * t159) * g(1) + (t7 * mrSges(6,1) - t13 * mrSges(5,3) - Ifges(5,4) * t565 + Ifges(6,6) * t564 + t652) * t202 + t334 * t605 / 0.2e1 + t176 * t606 / 0.2e1 + (t611 / 0.2e1 - t69 * mrSges(4,3) - Ifges(4,2) * t176 - t606 * qJDD(1) / 0.2e1 + Ifges(6,5) * t562 + Ifges(5,6) * t563 + Ifges(6,4) * t564 + Ifges(5,5) * t565 + t642 * t556 + (-qJDD(3) / 0.2e1 - t350 / 0.2e1) * Ifges(4,6) + (-t339 / 0.2e1 - t334 / 0.2e1) * Ifges(4,4) + t581) * t248 + (mrSges(6,1) * t8 - mrSges(5,3) * t14 + Ifges(5,4) * t563 - Ifges(6,6) * t562 + t653) * t203 + (t635 / 0.2e1 - t134 * mrSges(4,3) - t639 / 0.2e1 - t146 / 0.2e1 - Ifges(4,4) * t543 + Ifges(5,5) * t550 + Ifges(6,4) * t551 + Ifges(6,5) * t552 + Ifges(5,6) * t553 + t642 * t547 - t579) * t238 - t287 * t454 + m(4) * (t117 * t134 + t142 * t70 + t143 * t69 + t164 * t190 + t184 * t424) + m(7) * (t1 * t20 + t15 * t4 + t16 * t3 + t19 * t2 + t22 * t39 + t48 * t6) + m(6) * (t18 * t81 + t27 * t53 + t30 * t52 + t40 * t68 + t58 * t7 + t59 * t8) + (-m(4) * t133 - t596) * t118 + qJD(2) * t288 * t438 + (Ifges(3,5) * t449 + Ifges(3,6) * t418 + Ifges(3,3) * t432) * t497 + (-t133 * mrSges(4,3) + t228 / 0.2e1 + t147 / 0.2e1 + Ifges(4,1) * t543 + t591) * t237 + (Ifges(7,4) * t83 + Ifges(7,2) * t82) * t560 + (Ifges(7,4) * t151 + Ifges(7,2) * t150) * t573 + t18 * t407 + t282 * t461 / 0.2e1 + m(5) * (t13 * t72 + t131 * t66 + t14 * t71 + t32 * t56 - t33 * t55) - t6 * t410 + t66 * t412 + t70 * t413 + t164 * t414 + t350 * (Ifges(4,5) * t249 + Ifges(4,3) * t282) / 0.2e1 + m(3) * (t260 * t285 + t261 * t286) + ((Ifges(3,1) * t325 + Ifges(3,4) * t495) * t449 + (Ifges(3,4) * t325 + Ifges(3,2) * t495) * t418 + (Ifges(3,5) * t325 + Ifges(3,6) * t495) * t432 + t314 * (-mrSges(3,1) * t495 + mrSges(3,2) * t325) + m(3) * (-pkin(1) * t314 + qJD(2) * t382) - pkin(1) * t378) * t327 + (Ifges(7,5) * t83 + Ifges(7,6) * t82) * t554 + (Ifges(7,5) * t151 + Ifges(7,6) * t150) * t570 + t175 * t424 + (Ifges(7,1) * t83 + Ifges(7,4) * t82) * t558 + (Ifges(7,1) * t151 + Ifges(7,4) * t150) * t574 + (Ifges(4,1) * t339 + Ifges(4,4) * t176 + Ifges(4,5) * qJDD(3)) * t650 + t260 * t376 + (t1 * t150 - t15 * t83 - t151 * t2 + t16 * t82) * mrSges(7,3) + t261 * t370 + t190 * t106 + t117 * t191 + t150 * t10 / 0.2e1 + t32 * t140 + t33 * t141 + t142 * t144 + t143 * t145 + t27 * t138 + t30 * t139 + t131 * t50 + t40 * t123 + t3 * t90 + t4 * t91 + t39 * (-mrSges(7,1) * t82 + mrSges(7,2) * t83) + t83 * t64 / 0.2e1 + t22 * t85 + t81 * t49 + t58 * t73 + t59 * t74 + t71 * t75 + t72 * t76 + t48 * t17 + t19 * t28 + t20 * t29 + (t634 / 0.2e1 + t52 * mrSges(6,1) + t55 * mrSges(5,3) + Ifges(5,1) * t550 + Ifges(5,4) * t553 + Ifges(7,5) * t558 - Ifges(6,2) * t551 - Ifges(6,6) * t552 + Ifges(7,6) * t560 + Ifges(7,3) * t554 - t547 * t644 + t589 + t656) * t149 + (-t56 * mrSges(5,3) + t53 * mrSges(6,1) + t98 / 0.2e1 - Ifges(5,4) * t550 + Ifges(6,6) * t551 + Ifges(6,3) * t552 - Ifges(5,2) * t553 + t566 - t643 * t547 - t593) * t148 - t282 * t646 + (t285 * t376 + t286 * t370 + t605 * t650 + Ifges(2,3)) * qJDD(1) + (-m(5) * (t341 + t530) - m(6) * (t337 + t530) - m(7) * t337 - t108 * mrSges(7,1) - t107 * mrSges(7,2) - m(3) * t479 - t284 * mrSges(3,1) + mrSges(3,2) * t363 - mrSges(3,3) * t457 - mrSges(2,1) * t540 + mrSges(2,2) * t538 - m(4) * t343 - t209 * mrSges(4,1) - mrSges(4,3) * t340 + t645 * t160 + (mrSges(4,2) - t648 + t660) * t208 + t594 * t161) * g(2) + t82 * t571 + t151 * t576; (t133 * t264 - t134 * t265 - t184 * t425 + t164 * t496 + (t539 * t70 + t330 * t69 + (-t133 * t330 + t134 * t539) * qJD(3)) * t326 + t600) * m(4) + (-t382 * t478 + t314 + t600) * m(3) - t175 * t425 + t496 * t106 + t623 * t91 + (-t50 - t49 + t144) * t459 - t375 * t29 - t288 * t419 + (-t264 * t68 + t8 * t289 - t7 * t290 + (-t18 * t539 + t477 * t68) * t326 + t617 * t53 + t600) * m(6) + (-t113 * t264 + t13 * t290 - t14 * t289 + (t113 * t477 - t539 * t66) * t326 - t617 * t56 + t600) * m(5) + t622 * t90 + (-t1 * t375 + t15 * t623 + t16 * t622 + t2 * t258 + t290 * t6 - t39 * t617 + t600) * m(7) + t378 + t669 * t191 + t258 * t28 + t287 * t455 + t640 * t289 + (t17 + t641) * t290 + t145 * t484 + (m(5) * t55 - t651) * t664 + t617 * (-t85 + t626) + t663 * (t123 - t621); (-m(7) * (-pkin(11) * t491 - t198 + t620) + mrSges(7,3) * t491 - m(5) * t443 - m(6) * (t443 + t620) - t582 * t207 + t583 * t204) * g(2) + (-pkin(3) * t66 + t55 * t86 - t56 * t87) * m(5) + (t18 * t301 - t52 * t79 - t53 * t77 + t630 * t68) * m(6) + (((-t329 * t56 + t332 * t55) * qJD(4) + t607) * m(5) + t640 * t329 + t641 * t332 + ((t329 * t53 + t332 * t52) * qJD(4) + t604) * m(6)) * pkin(10) + (-t396 * t570 - t399 * t573 - t403 * t574 + t6 * t409 + t64 * t445 - t652) * t332 + (t52 * t610 + t53 * t609 + t587 + t604) * mrSges(6,1) + (t55 * t610 - t56 * t609 + t587 + t607) * mrSges(5,3) + (Ifges(4,1) * t544 + t394 * t553 + t395 * t550 + t402 * t552 + t405 * t551 + t548 * t613 - t591 - t614) * t236 + t653 * t329 + (Ifges(7,5) * t166 + Ifges(7,6) * t165) * t555 + (t464 - t79) * t139 - t510 * t564 + t517 * t565 + (t503 + t596) * t134 + ((-t402 / 0.2e1 + t394 / 0.2e1) * t188 + (t405 / 0.2e1 - t395 / 0.2e1) * t189 + (Ifges(7,3) * t332 + t329 * t396) * t554 + (Ifges(7,5) * t332 + t329 * t403) * t558 + (Ifges(7,6) * t332 + t329 * t399) * t560 + t614 + t613 * t547) * qJD(4) + (mrSges(7,1) * t625 - mrSges(7,2) * t624) * t39 + (-t1 * t480 + t15 * t624 - t16 * t625 + t2 * t482) * mrSges(7,3) + (Ifges(7,4) * t166 + Ifges(7,2) * t165) * t561 + (-t520 + t635) * t544 + t589 * t473 + (-t465 - t87) * t140 + (t100 / 0.2e1 + t580) * t488 + (t228 + t147) * t546 + (t101 / 0.2e1 - t98 / 0.2e1) * t489 + t629 * t85 + t630 * t123 + t631 * t91 + t632 * t90 + (t1 * t252 + t631 * t15 + t632 * t16 + t2 * t251 + t309 * t6 + t629 * t39) * m(7) + (Ifges(7,1) * t166 + Ifges(7,4) * t165) * t559 - t511 * t562 + t518 * t563 + t301 * t49 + t18 * t406 + (t331 * t63 + t501 + t98) * t474 / 0.2e1 + t66 * t411 + (-m(6) * (t441 + t618) - m(7) * (-pkin(11) * t487 - t247 + t618) + mrSges(7,3) * t487 + t414 - m(5) * t441 + t585 * t249 + t586 * t248) * g(3) + (-m(7) * (-pkin(11) * t490 - t200 + t619) + mrSges(7,3) * t490 - m(5) * t442 - m(6) * (t442 + t619) + t582 * t209 + t583 * t208) * g(1) + (-t397 * t554 - t400 * t560 - t404 * t558) * t471 + (-t464 - t86) * t141 - t11 * t482 / 0.2e1 - t646 + t251 * t28 + t252 * t29 + (-t191 + t504) * t133 + (t465 - t77) * t138 - t166 * t64 / 0.2e1 + t70 * mrSges(4,1) - pkin(3) * t50 + t461 + t309 * t17 + (Ifges(6,4) * t550 + Ifges(5,5) * t551 + Ifges(6,5) * t553 - Ifges(4,2) * t546 + Ifges(5,6) * t552 + t548 * t642 + t579) * t239 + t634 * (t473 / 0.2e1 - t488 / 0.2e1) + t146 * t543 + t474 * t566 + t452 * t571 + t165 * t572 + t480 * t577; (-m(7) * t356 - t662) * t569 + (qJD(6) * t507 - t533 + (-t472 - t492) * t16 + t588) * mrSges(7,3) + t611 + (t523 + t651) * t56 + (-t158 * t657 - t159 * t592) * g(2) + (-pkin(4) * t8 - qJ(5) * t7 - qJD(5) * t53 - t121 * t68) * m(6) + (-m(6) * t53 + t524 - t626) * t55 + t63 * t445 + t599 * qJD(6) + (-t506 + t98) * t551 + (t505 + t101) * t550 + (-m(6) * t433 + t407 + t412 - m(7) * (-t197 - t532) - t377 * t203) * g(3) + t385 * t85 + (t185 + t100) * t553 - t633 * qJD(5) + (-Ifges(5,2) * t189 - t186 + t634) * t552 + (t6 * qJ(5) - t15 * t23 - t16 * t24 + t654 * t39) * m(7) - t53 * t525 + t6 * t408 - (t136 * t399 + t137 * t403 + t187 * t396) * qJD(6) / 0.2e1 - t1 * t522 + t581 + (-t73 + t17) * qJ(5) - t121 * t123 - t24 * t90 - t23 * t91 - pkin(4) * t74 + (-Ifges(5,1) * t551 + Ifges(6,2) * t550 + t548 * t644 - t580 + t656) * t188 + (t160 * t657 + t161 * t592) * g(1) + (Ifges(6,3) * t553 + t15 * t522 + t396 * t555 + t399 * t561 + t403 * t559 - t548 * t643 + t593 + t599) * t189 + t52 * t526 + t397 * t570 + t492 * t572 + t400 * t573 + t404 * t574 + t331 * t576 + t328 * t577; t633 * t229 + (t123 + t391) * t189 + t74 + (-t189 * t392 - t229 * t39 + t356 - t588) * m(7) + (t189 * t68 + t229 * t53 - t588 + t8) * m(6) + t662; -t39 * (mrSges(7,1) * t137 + mrSges(7,2) * t136) + (Ifges(7,1) * t136 - t516) * t559 + t63 * t558 + (Ifges(7,5) * t136 - Ifges(7,6) * t137) * t555 - t15 * t90 + t16 * t91 - g(1) * (mrSges(7,1) * t107 - mrSges(7,2) * t108) - g(2) * ((-t158 * t331 - t204 * t328) * mrSges(7,1) + (t158 * t328 - t204 * t331) * mrSges(7,2)) - g(3) * t410 + (t136 * t15 + t137 * t16) * mrSges(7,3) + t9 + (-Ifges(7,2) * t137 + t135 + t64) * t561 + t597;];
tau  = t25;
