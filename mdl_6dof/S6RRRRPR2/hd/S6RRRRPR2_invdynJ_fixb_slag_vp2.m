% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:25
% EndTime: 2019-03-09 21:58:08
% DurationCPUTime: 25.46s
% Computational Cost: add. (28308->857), mult. (65068->1113), div. (0->0), fcn. (49333->18), ass. (0->417)
t410 = sin(pkin(11));
t411 = cos(pkin(11));
t413 = sin(qJ(6));
t418 = cos(qJ(6));
t453 = t410 * t413 - t411 * t418;
t311 = t453 * qJD(6);
t415 = sin(qJ(3));
t416 = sin(qJ(2));
t419 = cos(qJ(3));
t420 = cos(qJ(2));
t452 = t415 * t416 - t419 * t420;
t443 = t452 * qJD(1);
t576 = cos(qJ(4));
t301 = t576 * t443;
t339 = t415 * t420 + t416 * t419;
t316 = t339 * qJD(1);
t414 = sin(qJ(4));
t259 = t316 * t414 + t301;
t662 = t453 * t259;
t674 = t311 + t662;
t335 = t410 * t418 + t411 * t413;
t312 = t335 * qJD(6);
t663 = t335 * t259;
t673 = t312 + t663;
t672 = -mrSges(6,3) - mrSges(7,3);
t409 = qJ(2) + qJ(3);
t403 = qJ(4) + t409;
t386 = sin(t403);
t406 = pkin(11) + qJ(6);
t398 = cos(t406);
t561 = mrSges(7,1) * t398;
t562 = mrSges(6,1) * t411;
t671 = (t561 + t562) * t386;
t397 = sin(t406);
t558 = mrSges(7,2) * t397;
t559 = mrSges(6,2) * t410;
t670 = (-t558 - t559) * t386;
t665 = t259 * t410;
t249 = pkin(5) * t665;
t502 = pkin(10) * t665;
t404 = t420 * pkin(2);
t390 = t404 + pkin(1);
t299 = pkin(3) * t452 - t390;
t563 = t411 * pkin(5);
t384 = pkin(4) + t563;
t387 = cos(t403);
t412 = -pkin(10) - qJ(5);
t454 = -t384 * t386 - t387 * t412;
t669 = -m(7) * t454 + t671;
t434 = t414 * t443;
t428 = t316 * t576 - t434;
t664 = t259 * t411;
t471 = pkin(5) * t428 + pkin(10) * t664;
t422 = -pkin(8) - pkin(7);
t369 = t422 * t420;
t347 = qJD(1) * t369;
t317 = t415 * t347;
t368 = t422 * t416;
t346 = qJD(1) * t368;
t328 = qJD(2) * pkin(2) + t346;
t270 = t419 * t328 + t317;
t308 = t316 * pkin(9);
t232 = t270 - t308;
t407 = qJD(2) + qJD(3);
t223 = pkin(3) * t407 + t232;
t320 = t419 * t347;
t271 = t415 * t328 - t320;
t436 = pkin(9) * t443;
t233 = t271 - t436;
t225 = t414 * t233;
t153 = t223 * t576 - t225;
t503 = -qJD(3) - qJD(4);
t399 = qJD(2) - t503;
t148 = -t399 * pkin(4) + qJD(5) - t153;
t477 = t411 * t399 - t410 * t428;
t108 = -pkin(5) * t477 + t148;
t505 = qJD(1) * qJD(2);
t350 = qJDD(1) * t420 - t416 * t505;
t351 = qJDD(1) * t416 + t420 * t505;
t442 = t452 * qJD(3);
t240 = -qJD(1) * t442 + t350 * t415 + t351 * t419;
t444 = t339 * qJD(3);
t241 = -qJD(1) * t444 + t350 * t419 - t351 * t415;
t140 = t576 * t240 + (-qJD(4) * t316 + t241) * t414 - qJD(4) * t301;
t483 = qJD(4) * t576;
t141 = -qJD(4) * t434 + t414 * t240 - t576 * t241 + t316 * t483;
t554 = Ifges(5,4) * t428;
t644 = Ifges(5,6) * t399;
t195 = -Ifges(5,2) * t259 + t554 + t644;
t405 = qJDD(2) + qJDD(3);
t396 = qJDD(4) + t405;
t114 = -t140 * t410 + t396 * t411;
t333 = t351 * pkin(7);
t282 = qJDD(2) * pkin(2) - pkin(8) * t351 - t333;
t332 = t350 * pkin(7);
t285 = pkin(8) * t350 + t332;
t179 = -qJD(3) * t271 + t419 * t282 - t285 * t415;
t118 = pkin(3) * t405 - pkin(9) * t240 + t179;
t507 = qJD(3) * t419;
t508 = qJD(3) * t415;
t178 = t415 * t282 + t419 * t285 + t328 * t507 + t347 * t508;
t129 = pkin(9) * t241 + t178;
t506 = qJD(4) * t414;
t51 = t414 * t118 + t576 * t129 + t223 * t483 - t233 * t506;
t43 = qJ(5) * t396 + qJD(5) * t399 + t51;
t542 = qJDD(1) * pkin(1);
t306 = -pkin(2) * t350 - t542;
t212 = -pkin(3) * t241 + t306;
t55 = pkin(4) * t141 - qJ(5) * t140 - qJD(5) * t428 + t212;
t15 = t410 * t55 + t411 * t43;
t11 = pkin(10) * t114 + t15;
t237 = t399 * t410 + t411 * t428;
t226 = t576 * t233;
t154 = t414 * t223 + t226;
t149 = t399 * qJ(5) + t154;
t288 = t299 * qJD(1);
t169 = t259 * pkin(4) - qJ(5) * t428 + t288;
t88 = -t149 * t410 + t411 * t169;
t61 = pkin(5) * t259 - pkin(10) * t237 + t88;
t89 = t411 * t149 + t410 * t169;
t70 = pkin(10) * t477 + t89;
t20 = -t413 * t70 + t418 * t61;
t115 = t140 * t411 + t396 * t410;
t14 = -t410 * t43 + t411 * t55;
t6 = pkin(5) * t141 - pkin(10) * t115 + t14;
t2 = qJD(6) * t20 + t11 * t418 + t413 * t6;
t21 = t413 * t61 + t418 * t70;
t52 = t118 * t576 - t414 * t129 - t223 * t506 - t233 * t483;
t44 = -t396 * pkin(4) + qJDD(5) - t52;
t22 = -t114 * pkin(5) + t44;
t3 = -qJD(6) * t21 - t11 * t413 + t418 * t6;
t39 = t115 * Ifges(6,4) + t114 * Ifges(6,2) + t141 * Ifges(6,6);
t457 = -t14 * t410 + t15 * t411;
t464 = t559 - t562;
t552 = Ifges(6,4) * t411;
t553 = Ifges(6,4) * t410;
t578 = -t399 / 0.2e1;
t581 = t428 / 0.2e1;
t584 = t259 / 0.2e1;
t585 = -t259 / 0.2e1;
t257 = qJD(6) + t259;
t587 = t257 / 0.2e1;
t588 = -t257 / 0.2e1;
t589 = -t237 / 0.2e1;
t590 = -t477 / 0.2e1;
t163 = t237 * t418 + t413 * t477;
t592 = t163 / 0.2e1;
t593 = -t163 / 0.2e1;
t656 = -t237 * t413 + t418 * t477;
t594 = t656 / 0.2e1;
t595 = -t656 / 0.2e1;
t596 = t141 / 0.2e1;
t139 = qJDD(6) + t141;
t597 = t139 / 0.2e1;
t598 = t115 / 0.2e1;
t599 = t114 / 0.2e1;
t156 = Ifges(7,4) * t656;
t83 = Ifges(7,1) * t163 + Ifges(7,5) * t257 + t156;
t600 = t83 / 0.2e1;
t551 = Ifges(7,4) * t163;
t82 = Ifges(7,2) * t656 + Ifges(7,6) * t257 + t551;
t602 = t82 / 0.2e1;
t49 = -qJD(6) * t163 + t114 * t418 - t115 * t413;
t604 = t49 / 0.2e1;
t48 = qJD(6) * t656 + t114 * t413 + t115 * t418;
t605 = t48 / 0.2e1;
t606 = Ifges(6,1) * t598 + Ifges(6,4) * t599 + Ifges(6,5) * t596;
t607 = Ifges(7,1) * t605 + Ifges(7,4) * t604 + Ifges(7,5) * t597;
t608 = Ifges(7,4) * t605 + Ifges(7,2) * t604 + Ifges(7,6) * t597;
t610 = t288 * mrSges(5,1) + t88 * mrSges(6,1) + t20 * mrSges(7,1) - t89 * mrSges(6,2) - t21 * mrSges(7,2) - t154 * mrSges(5,3);
t643 = Ifges(6,6) * t477;
t645 = Ifges(6,5) * t237;
t641 = Ifges(7,5) * t163 + Ifges(7,6) * t656 + Ifges(6,3) * t259 + Ifges(7,3) * t257 + t643 + t645;
t650 = -t428 / 0.2e1;
t668 = -t662 * t83 / 0.2e1 - t663 * t82 / 0.2e1 + (-Ifges(7,5) * t311 - Ifges(7,6) * t312) * t587 + (-Ifges(7,1) * t311 - Ifges(7,4) * t312) * t592 + (-Ifges(7,4) * t311 - Ifges(7,2) * t312) * t594 + (Ifges(6,5) * t589 + Ifges(7,5) * t593 - Ifges(5,2) * t584 - Ifges(5,6) * t578 + Ifges(6,6) * t590 + Ifges(7,6) * t595 + Ifges(6,3) * t585 + Ifges(7,3) * t588 - t610) * t428 + (Ifges(7,5) * t662 + Ifges(7,6) * t663) * t588 + (Ifges(7,4) * t662 + Ifges(7,2) * t663) * t595 + (Ifges(7,1) * t662 + Ifges(7,4) * t663) * t593 + (-t2 * t453 + t20 * t674 - t673 * t21 - t3 * t335) * mrSges(7,3) + (t673 * mrSges(7,1) - mrSges(7,2) * t674) * t108 + (Ifges(7,5) * t335 - Ifges(7,6) * t453) * t597 + (Ifges(7,4) * t335 - Ifges(7,2) * t453) * t604 + (Ifges(7,1) * t335 - Ifges(7,4) * t453) * t605 + t22 * (mrSges(7,1) * t453 + mrSges(7,2) * t335) - t453 * t608 + t195 * t581 + (-t554 + t641) * t650 + (Ifges(6,5) * t410 + Ifges(6,6) * t411) * t596 + (Ifges(6,1) * t410 + t552) * t598 + (Ifges(6,2) * t411 + t553) * t599 - t311 * t600 - t312 * t602 + t410 * t606 + t335 * t607 + t411 * t39 / 0.2e1 + Ifges(5,3) * t396 + (-t664 * t88 - t665 * t89 + t457) * mrSges(6,3) + Ifges(5,5) * t140 - Ifges(5,6) * t141 + t44 * t464 - t51 * mrSges(5,2) + t52 * mrSges(5,1);
t628 = t387 * t384 - t386 * t412;
t667 = m(7) * t628;
t661 = -t387 * mrSges(5,1) + (mrSges(5,2) - mrSges(7,3)) * t386;
t400 = sin(t409);
t401 = cos(t409);
t465 = mrSges(5,1) * t386 + mrSges(5,2) * t387;
t660 = mrSges(4,1) * t400 + mrSges(4,2) * t401 + t465;
t421 = cos(qJ(1));
t525 = t387 * t421;
t659 = t670 * t421 + t525 * t672;
t417 = sin(qJ(1));
t526 = t387 * t417;
t658 = t670 * t417 + t526 * t672;
t622 = g(1) * t421 + g(2) * t417;
t204 = pkin(4) * t428 + qJ(5) * t259;
t463 = mrSges(6,1) * t410 + mrSges(6,2) * t411;
t655 = t288 * mrSges(5,2) - t410 * (Ifges(6,4) * t237 + Ifges(6,2) * t477 + Ifges(6,6) * t259) / 0.2e1 - t153 * mrSges(5,3) + t148 * t463;
t458 = Ifges(6,5) * t411 - Ifges(6,6) * t410;
t460 = -Ifges(6,2) * t410 + t552;
t462 = Ifges(6,1) * t411 - t553;
t654 = -Ifges(5,1) * t650 - Ifges(5,5) * t578 - t458 * t585 - t460 * t590 - t462 * t589;
t651 = t350 / 0.2e1;
t577 = t420 / 0.2e1;
t649 = m(5) * t153;
t646 = Ifges(5,5) * t399;
t642 = t420 * Ifges(3,2);
t570 = pkin(3) * t414;
t383 = qJ(5) + t570;
t324 = (-pkin(10) - t383) * t410;
t402 = t411 * pkin(10);
t528 = t383 * t411;
t325 = t402 + t528;
t263 = t324 * t418 - t325 * t413;
t474 = pkin(3) * t483;
t373 = t474 + qJD(5);
t164 = t232 * t576 - t225;
t572 = pkin(3) * t316;
t180 = t204 + t572;
t92 = -t164 * t410 + t411 * t180;
t68 = t471 + t92;
t93 = t411 * t164 + t410 * t180;
t84 = t502 + t93;
t640 = qJD(6) * t263 - t373 * t453 - t413 * t68 - t418 * t84;
t264 = t324 * t413 + t325 * t418;
t639 = -qJD(6) * t264 - t335 * t373 + t413 * t84 - t418 * t68;
t358 = t412 * t410;
t543 = qJ(5) * t411;
t359 = t402 + t543;
t283 = t358 * t418 - t359 * t413;
t97 = -t153 * t410 + t411 * t204;
t71 = t471 + t97;
t98 = t411 * t153 + t410 * t204;
t87 = t98 + t502;
t638 = -qJD(5) * t453 + qJD(6) * t283 - t413 * t71 - t418 * t87;
t284 = t358 * t413 + t359 * t418;
t637 = -qJD(5) * t335 - qJD(6) * t284 + t413 * t87 - t418 * t71;
t573 = pkin(2) * t419;
t389 = pkin(3) + t573;
t487 = t576 * t415;
t310 = pkin(2) * t487 + t414 * t389;
t304 = qJ(5) + t310;
t286 = (-pkin(10) - t304) * t410;
t531 = t304 * t411;
t287 = t402 + t531;
t214 = t286 * t413 + t287 * t418;
t575 = pkin(2) * t415;
t380 = t414 * t575;
t495 = pkin(2) * t507;
t266 = t380 * t503 + t389 * t483 + t576 * t495;
t265 = qJD(5) + t266;
t277 = -t415 * t346 + t320;
t242 = t277 + t436;
t278 = t419 * t346 + t317;
t243 = -t308 + t278;
t171 = t414 * t242 + t243 * t576;
t511 = qJD(1) * t416;
t393 = pkin(2) * t511;
t172 = t180 + t393;
t94 = -t171 * t410 + t411 * t172;
t69 = t471 + t94;
t95 = t411 * t171 + t410 * t172;
t85 = t502 + t95;
t636 = -qJD(6) * t214 - t265 * t335 + t413 * t85 - t418 * t69;
t213 = t286 * t418 - t287 * t413;
t635 = qJD(6) * t213 - t265 * t453 - t413 * t69 - t418 * t85;
t631 = mrSges(5,1) * t399 + mrSges(6,1) * t477 - mrSges(6,2) * t237 - mrSges(5,3) * t428;
t289 = t419 * t368 + t369 * t415;
t254 = -pkin(9) * t339 + t289;
t290 = t415 * t368 - t419 * t369;
t255 = -pkin(9) * t452 + t290;
t630 = t576 * t254 - t414 * t255;
t629 = t266 - t171;
t479 = t401 * mrSges(4,1) - mrSges(4,2) * t400;
t569 = pkin(4) * t386;
t571 = pkin(3) * t400;
t626 = -m(7) * (t454 - t571) - m(6) * (-t569 - t571) + t671;
t510 = qJD(1) * t420;
t566 = pkin(7) * t420;
t567 = pkin(7) * t416;
t625 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t511) * t566 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t510) * t567;
t624 = t332 * t420 + t333 * t416;
t182 = -mrSges(6,2) * t259 + mrSges(6,3) * t477;
t183 = mrSges(6,1) * t259 - mrSges(6,3) * t237;
t623 = t182 * t411 - t183 * t410;
t376 = t386 * mrSges(6,3);
t620 = -t376 + t661 + (t558 - t561 + t464) * t387;
t62 = -t114 * mrSges(6,1) + t115 * mrSges(6,2);
t619 = m(6) * t44 + t62;
t618 = -t479 + t620;
t617 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t616 = m(6) * t148 - t631;
t408 = -pkin(9) + t422;
t615 = -m(3) * pkin(7) + m(4) * t422 + m(6) * t408 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - t463;
t90 = -mrSges(7,1) * t656 + mrSges(7,2) * t163;
t614 = m(7) * t108 + t616 + t90;
t366 = -t420 * mrSges(3,1) + t416 * mrSges(3,2);
t613 = -m(4) * t390 - (m(6) * pkin(4) - t464) * t387 - mrSges(2,1) - m(3) * pkin(1) + t366 - t479 + t661;
t612 = (t421 * t669 + t659) * g(1) + (t417 * t669 + t658) * g(2);
t519 = t411 * (t237 * Ifges(6,1) + Ifges(6,4) * t477 + t259 * Ifges(6,5));
t611 = t519 / 0.2e1 + t655;
t579 = t316 / 0.2e1;
t574 = pkin(2) * t416;
t385 = pkin(3) * t401;
t568 = pkin(5) * t410;
t279 = -qJD(2) * t452 - t442;
t280 = -qJD(2) * t339 - t444;
t448 = -t414 * t339 - t452 * t576;
t173 = qJD(4) * t448 + t279 * t576 + t414 * t280;
t276 = t339 * t576 - t414 * t452;
t174 = qJD(4) * t276 + t414 * t279 - t280 * t576;
t509 = qJD(2) * t416;
t394 = pkin(2) * t509;
t262 = -pkin(3) * t280 + t394;
t77 = pkin(4) * t174 - qJ(5) * t173 - qJD(5) * t276 + t262;
t488 = qJD(2) * t422;
t348 = t416 * t488;
t349 = t420 * t488;
t217 = t419 * t348 + t415 * t349 + t368 * t507 + t369 * t508;
t186 = pkin(9) * t280 + t217;
t218 = -qJD(3) * t290 - t348 * t415 + t419 * t349;
t187 = -pkin(9) * t279 + t218;
t78 = qJD(4) * t630 + t576 * t186 + t414 * t187;
t26 = t410 * t77 + t411 * t78;
t557 = Ifges(3,4) * t416;
t556 = Ifges(3,4) * t420;
t555 = Ifges(4,4) * t316;
t546 = t316 * mrSges(4,3);
t64 = mrSges(6,1) * t141 - mrSges(6,3) * t115;
t544 = t410 * t64;
t541 = t173 * t410;
t540 = t173 * t411;
t533 = t276 * t410;
t532 = t276 * t411;
t374 = t386 * qJ(5);
t524 = t397 * t417;
t523 = t397 * t421;
t522 = t398 * t417;
t521 = t398 * t421;
t193 = -pkin(4) * t448 - qJ(5) * t276 + t299;
t198 = t414 * t254 + t255 * t576;
t103 = t410 * t193 + t411 * t198;
t513 = t387 * pkin(4) + t374;
t512 = t385 + t404;
t499 = Ifges(7,5) * t48 + Ifges(7,6) * t49 + Ifges(7,3) * t139;
t498 = t576 * pkin(3);
t490 = t385 + t513;
t12 = -t49 * mrSges(7,1) + t48 * mrSges(7,2);
t25 = -t410 * t78 + t411 * t77;
t480 = t505 / 0.2e1;
t102 = t411 * t193 - t198 * t410;
t161 = t232 * t414 + t226;
t170 = -t576 * t242 + t243 * t414;
t388 = -t498 - pkin(4);
t356 = qJ(5) * t526;
t473 = -t417 * t569 + t356;
t357 = qJ(5) * t525;
t472 = -t421 * t569 + t357;
t469 = t385 + t628;
t309 = t389 * t576 - t380;
t468 = mrSges(3,1) * t416 + mrSges(3,2) * t420;
t461 = t557 + t642;
t459 = Ifges(3,5) * t420 - Ifges(3,6) * t416;
t456 = -t410 * t88 + t411 * t89;
t86 = -pkin(5) * t448 - pkin(10) * t532 + t102;
t91 = -pkin(10) * t533 + t103;
t33 = -t413 * t91 + t418 * t86;
t34 = t413 * t86 + t418 * t91;
t305 = -pkin(4) - t309;
t449 = pkin(1) * t468;
t447 = t416 * (Ifges(3,1) * t420 - t557);
t435 = mrSges(4,3) * t443;
t79 = qJD(4) * t198 + t414 * t186 - t576 * t187;
t253 = Ifges(5,4) * t259;
t196 = Ifges(5,1) * t428 - t253 + t646;
t251 = -Ifges(4,2) * t443 + Ifges(4,6) * t407 + t555;
t307 = Ifges(4,4) * t443;
t252 = t316 * Ifges(4,1) + t407 * Ifges(4,5) - t307;
t367 = t390 * qJD(1);
t423 = (-Ifges(5,4) * t584 + t654 + t655) * t259 + t251 * t579 + (-Ifges(4,2) * t316 + t252 - t307) * t443 / 0.2e1 + t367 * (t316 * mrSges(4,1) - mrSges(4,2) * t443) - t316 * (-Ifges(4,1) * t443 - t555) / 0.2e1 - t407 * (-Ifges(4,5) * t443 - Ifges(4,6) * t316) / 0.2e1 + (t519 + t196) * t584 + t668 + t271 * t546 - t270 * t435 + Ifges(4,3) * t405 + Ifges(4,5) * t240 + Ifges(4,6) * t241 - t178 * mrSges(4,2) + t179 * mrSges(4,1);
t392 = Ifges(3,4) * t510;
t355 = t388 - t563;
t352 = -t571 - t574;
t345 = pkin(1) + t512;
t330 = t421 * t352;
t329 = t417 * t352;
t323 = t421 * t345;
t315 = Ifges(3,1) * t511 + Ifges(3,5) * qJD(2) + t392;
t314 = Ifges(3,6) * qJD(2) + qJD(1) * t461;
t298 = t305 - t563;
t297 = t387 * t521 + t524;
t296 = -t387 * t523 + t522;
t295 = -t387 * t522 + t523;
t294 = t387 * t524 + t521;
t293 = mrSges(4,1) * t407 - t546;
t292 = -t407 * mrSges(4,2) - t435;
t291 = t393 + t572;
t269 = mrSges(4,1) * t443 + t316 * mrSges(4,2);
t244 = -mrSges(5,2) * t399 - mrSges(5,3) * t259;
t222 = -mrSges(4,2) * t405 + mrSges(4,3) * t241;
t221 = mrSges(4,1) * t405 - mrSges(4,3) * t240;
t208 = t453 * t276;
t207 = t335 * t276;
t206 = mrSges(5,1) * t259 + mrSges(5,2) * t428;
t147 = pkin(5) * t533 - t630;
t133 = t170 - t249;
t126 = t161 - t249;
t125 = mrSges(7,1) * t257 - mrSges(7,3) * t163;
t124 = -mrSges(7,2) * t257 + mrSges(7,3) * t656;
t121 = -mrSges(5,2) * t396 - mrSges(5,3) * t141;
t120 = mrSges(5,1) * t396 - mrSges(5,3) * t140;
t116 = t154 - t249;
t67 = -t173 * t335 + t276 * t311;
t66 = -t173 * t453 - t276 * t312;
t63 = -mrSges(6,2) * t141 + mrSges(6,3) * t114;
t56 = pkin(5) * t541 + t79;
t24 = -mrSges(7,2) * t139 + mrSges(7,3) * t49;
t23 = mrSges(7,1) * t139 - mrSges(7,3) * t48;
t17 = -pkin(10) * t541 + t26;
t16 = pkin(5) * t174 - pkin(10) * t540 + t25;
t5 = -qJD(6) * t34 + t16 * t418 - t17 * t413;
t4 = qJD(6) * t33 + t16 * t413 + t17 * t418;
t1 = [t22 * (mrSges(7,1) * t207 - mrSges(7,2) * t208) + (-Ifges(7,1) * t208 - Ifges(7,4) * t207) * t605 + (-Ifges(7,5) * t208 - Ifges(7,6) * t207) * t597 + (-Ifges(7,4) * t208 - Ifges(7,2) * t207) * t604 + (-t2 * t207 - t20 * t66 + t208 * t3 + t21 * t67) * mrSges(7,3) + m(6) * (t102 * t14 + t103 * t15 + t25 * t88 + t26 * t89) + m(5) * (t154 * t78 + t198 * t51 + t212 * t299 + t262 * t288) - (-m(5) * t52 - t120 + t619) * t630 + m(4) * (t178 * t290 + t179 * t289 + t217 * t271 + t218 * t270 - t306 * t390 - t367 * t394) + (Ifges(7,1) * t66 + Ifges(7,4) * t67) * t592 + t351 * t556 / 0.2e1 + t461 * t651 - (Ifges(6,5) * t115 + Ifges(6,6) * t114 + Ifges(6,3) * t141 + t499) * t448 / 0.2e1 - (t212 * mrSges(5,1) + t14 * mrSges(6,1) - t15 * mrSges(6,2) - t51 * mrSges(5,3) - Ifges(5,4) * t140 + Ifges(6,5) * t598 + Ifges(7,5) * t605 + Ifges(5,2) * t141 - Ifges(5,6) * t396 + Ifges(6,6) * t599 + Ifges(7,6) * t604 + Ifges(6,3) * t596 + Ifges(7,3) * t597 + t617) * t448 - (-mrSges(4,1) * t306 + mrSges(4,3) * t178 + Ifges(4,4) * t240 + Ifges(4,2) * t241 + Ifges(4,6) * t405) * t452 + (-t295 * mrSges(7,1) - t294 * mrSges(7,2) + (-m(7) * (-t408 + t568) + m(5) * t408 + t615) * t421 + (-m(7) * (-t345 - t628) - m(6) * (-t345 - t374) + t376 + m(5) * t345 - t613) * t417) * g(1) + (Ifges(4,1) * t279 + Ifges(4,4) * t280) * t579 + (Ifges(5,1) * t581 + t458 * t584 + Ifges(5,4) * t585 + t646 / 0.2e1 + t196 / 0.2e1 + t477 * t460 / 0.2e1 + t237 * t462 / 0.2e1 + t611) * t173 + (t641 / 0.2e1 - Ifges(5,4) * t581 + Ifges(6,3) * t584 - Ifges(5,2) * t585 + Ifges(7,3) * t587 + Ifges(7,5) * t592 + Ifges(7,6) * t594 - t644 / 0.2e1 - t195 / 0.2e1 + t643 / 0.2e1 + t645 / 0.2e1 + t610) * t174 + (Ifges(7,4) * t66 + Ifges(7,2) * t67) * t594 + (t212 * mrSges(5,2) - t52 * mrSges(5,3) + Ifges(5,1) * t140 - Ifges(5,4) * t141 + Ifges(5,5) * t396 + t44 * t463 + t458 * t596 + t460 * t599 + t462 * t598) * t276 + (Ifges(7,5) * t66 + Ifges(7,6) * t67) * t587 + t66 * t600 + t67 * t602 + t532 * t606 - t208 * t607 - t207 * t608 + t269 * t394 + (mrSges(4,2) * t306 - mrSges(4,3) * t179 + Ifges(4,1) * t240 + Ifges(4,4) * t241 + Ifges(4,5) * t405) * t339 + (t350 * t566 + t351 * t567 + t624) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t624) + (t315 * t577 + t459 * qJD(2) / 0.2e1 - t625) * qJD(2) + (t420 * t556 + t447) * t480 + (-t14 * t532 - t15 * t533 - t540 * t88 - t541 * t89) * mrSges(6,3) + m(7) * (t108 * t56 + t147 * t22 + t2 * t34 + t20 * t5 + t21 * t4 + t3 * t33) - t366 * t542 + Ifges(2,3) * qJDD(1) + t407 * (Ifges(4,5) * t279 + Ifges(4,6) * t280) / 0.2e1 - t39 * t533 / 0.2e1 - t390 * (-mrSges(4,1) * t241 + mrSges(4,2) * t240) - t367 * (-mrSges(4,1) * t280 + mrSges(4,2) * t279) - pkin(1) * (-mrSges(3,1) * t350 + mrSges(3,2) * t351) - t314 * t509 / 0.2e1 - t449 * t505 + t299 * (mrSges(5,1) * t141 + mrSges(5,2) * t140) + t290 * t222 + t217 * t292 + t218 * t293 + t289 * t221 + t279 * t252 / 0.2e1 + t280 * t251 / 0.2e1 + t262 * t206 + t78 * t244 + t198 * t121 + t25 * t183 + t26 * t182 + t147 * t12 + t4 * t124 + t5 * t125 + t108 * (-mrSges(7,1) * t67 + mrSges(7,2) * t66) + t102 * t64 + t103 * t63 + t56 * t90 + (-t270 * t279 + t271 * t280) * mrSges(4,3) + (-m(6) * t323 - t297 * mrSges(7,1) - t296 * mrSges(7,2) + (-m(7) - m(5)) * (-t408 * t417 + t323) + (-m(7) * t568 + t615) * t417 + (-t667 - (m(6) * qJ(5) + mrSges(6,3)) * t386 + t613) * t421) * g(2) + (t616 - t649) * t79 + (-mrSges(3,1) * t567 - mrSges(3,2) * t566 + 0.2e1 * Ifges(3,6) * t577) * qJDD(2) + (Ifges(3,1) * t351 + Ifges(3,4) * t651 + Ifges(3,5) * qJDD(2) - t480 * t642) * t416 - (Ifges(4,4) * t279 + Ifges(4,2) * t280) * t443 / 0.2e1 + t33 * t23 + t34 * t24 + (Ifges(3,4) * t351 + Ifges(3,2) * t350) * t577; -(-Ifges(3,2) * t511 + t315 + t392) * t510 / 0.2e1 + ((t178 * t415 + t179 * t419 + (-t270 * t415 + t271 * t419) * qJD(3)) * pkin(2) - t270 * t277 - t271 * t278 + t367 * t393) * m(4) + t635 * t124 + (-g(1) * t330 - g(2) * t329 - t108 * t133 + t2 * t214 + t20 * t636 + t21 * t635 + t213 * t3 + t22 * t298) * m(7) + t636 * t125 + t622 * (m(4) * t574 - m(5) * t352 + t468 + t660) + (-m(4) * t404 + t366 - m(5) * t512 - m(6) * (t404 + t490) - m(7) * (t404 + t469) + t618) * g(3) + t221 * t573 + t222 * t575 + (t625 + (t449 - t447 / 0.2e1) * qJD(1)) * qJD(1) + (t495 - t278) * t292 + (t265 * t456 + t304 * t457 + t305 * t44 - t148 * t170 - t88 * t94 - t89 * t95 - g(1) * (t330 + t472) - g(2) * (t329 + t473)) * m(6) + t623 * t265 + t612 + (-pkin(2) * t508 - t277) * t293 + t63 * t531 - t304 * t544 + Ifges(3,6) * t350 + Ifges(3,5) * t351 - t332 * mrSges(3,2) - t333 * mrSges(3,1) + t314 * t511 / 0.2e1 - t459 * t505 / 0.2e1 + t309 * t120 + t310 * t121 + t305 * t62 + t298 * t12 - t269 * t393 - t291 * t206 + Ifges(3,3) * qJDD(2) + t213 * t23 + t214 * t24 - t94 * t183 - t95 * t182 - t133 * t90 + (t153 * t170 + t154 * t629 - t288 * t291 + t309 * t52 + t310 * t51) * m(5) + t629 * t244 + t631 * t170 + (t614 - t649) * (t389 * t506 + (t415 * t483 + (t414 * t419 + t487) * qJD(3)) * pkin(2)) + t423; t639 * t125 + (-t108 * t126 + t2 * t264 + t20 * t639 + t21 * t640 + t22 * t355 + t263 * t3) * m(7) + t640 * t124 + (t421 * t626 + t659) * g(1) + (m(5) * t571 + t660) * t622 + (t417 * t626 + t658) * g(2) + (-m(5) * t385 - m(6) * t490 - m(7) * t469 + t618) * g(3) - t206 * t572 + t121 * t570 + t614 * pkin(3) * t506 + t623 * t373 + t63 * t528 + t120 * t498 - t383 * t544 + t388 * t62 + t355 * t12 - t270 * t292 + t271 * t293 + t263 * t23 + t264 * t24 - t92 * t183 - t93 * t182 - t126 * t90 + (-g(1) * t357 - g(2) * t356 - t148 * t161 + t373 * t456 + t383 * t457 + t388 * t44 - t88 * t92 - t89 * t93) * m(6) + t631 * t161 + (t153 * t161 - t154 * t164 - t288 * t572 + (t576 * t52 + t414 * t51 + (-t153 * t414 + t154 * t576) * qJD(4)) * pkin(3)) * m(5) + (-t164 + t474) * t244 + t423; t637 * t125 + (-t108 * t116 + t2 * t284 + t20 * t637 + t21 * t638 - t22 * t384 + t283 * t3) * m(7) + t638 * t124 + (t611 + t654) * t259 + (-pkin(4) * t44 - g(1) * t472 - g(2) * t473 + qJ(5) * t457 + qJD(5) * t456 - t148 * t154 - t88 * t97 - t89 * t98) * m(6) + t623 * qJD(5) + t612 + t668 + (-t253 + t196) * t584 + t622 * t465 + t63 * t543 - qJ(5) * t544 - t384 * t12 + t283 * t23 + t284 * t24 - t153 * t244 - t97 * t183 - t98 * t182 - t116 * t90 - pkin(4) * t62 + (-m(6) * t513 + t620 - t667) * g(3) + t631 * t154; -t656 * t124 + t163 * t125 - t477 * t182 + t237 * t183 - m(6) * (-t237 * t88 + t477 * t89) + t12 + t619 + (g(3) * t387 - t386 * t622) * (m(6) + m(7)) + (t163 * t20 - t21 * t656 + t22) * m(7); -t108 * (mrSges(7,1) * t163 + mrSges(7,2) * t656) + (Ifges(7,1) * t656 - t551) * t593 + t82 * t592 + (Ifges(7,5) * t656 - Ifges(7,6) * t163) * t588 - t20 * t124 + t21 * t125 - g(1) * (mrSges(7,1) * t296 - mrSges(7,2) * t297) - g(2) * (-mrSges(7,1) * t294 + mrSges(7,2) * t295) - g(3) * (-mrSges(7,1) * t397 - mrSges(7,2) * t398) * t386 + (t163 * t21 + t20 * t656) * mrSges(7,3) + t499 + (-Ifges(7,2) * t163 + t156 + t83) * t595 + t617;];
tau  = t1;
