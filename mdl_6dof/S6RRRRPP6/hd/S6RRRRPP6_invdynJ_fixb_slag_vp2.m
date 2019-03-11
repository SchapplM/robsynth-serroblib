% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:11:16
% EndTime: 2019-03-09 21:12:26
% DurationCPUTime: 44.36s
% Computational Cost: add. (11510->842), mult. (24995->1011), div. (0->0), fcn. (16955->10), ass. (0->411)
t641 = -mrSges(7,2) - mrSges(6,3);
t564 = -mrSges(5,2) - t641;
t642 = mrSges(5,1) - mrSges(6,2);
t654 = Ifges(5,4) + Ifges(6,6);
t320 = cos(qJ(2));
t443 = qJD(1) * t320;
t288 = qJD(3) - t443;
t274 = -qJD(4) - t288;
t518 = t274 / 0.2e1;
t563 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t653 = t563 * t518;
t317 = sin(qJ(2));
t322 = -pkin(9) - pkin(8);
t388 = m(7) * (pkin(5) - t322) + mrSges(7,1);
t640 = -mrSges(5,3) - mrSges(6,1);
t566 = t640 * t317;
t652 = -t388 * t317 + t566;
t444 = qJD(1) * t317;
t302 = pkin(7) * t444;
t261 = -qJD(2) * pkin(2) + t302;
t316 = sin(qJ(3));
t413 = t316 * t444;
t319 = cos(qJ(3));
t441 = qJD(2) * t319;
t357 = t413 - t441;
t177 = pkin(3) * t357 + t261;
t387 = pkin(2) * t320 + pkin(8) * t317;
t253 = -pkin(1) - t387;
t232 = t253 * qJD(1);
t303 = pkin(7) * t443;
t262 = qJD(2) * pkin(8) + t303;
t163 = t316 * t232 + t319 * t262;
t125 = -pkin(9) * t357 + t163;
t516 = cos(qJ(4));
t120 = t516 * t125;
t315 = sin(qJ(4));
t162 = t232 * t319 - t316 * t262;
t417 = t319 * t444;
t243 = qJD(2) * t316 + t417;
t124 = -pkin(9) * t243 + t162;
t336 = pkin(3) * t288 + t124;
t335 = t315 * t336;
t62 = t120 + t335;
t50 = t274 * qJ(5) - t62;
t215 = t516 * t357;
t156 = t243 * t315 + t215;
t612 = t156 * pkin(5) - qJD(6);
t30 = -t50 - t612;
t331 = t243 * t516 - t315 * t357;
t329 = -qJ(5) * t331 + t177;
t501 = pkin(4) + qJ(6);
t46 = t156 * t501 + t329;
t68 = t156 * pkin(4) + t329;
t488 = t331 * Ifges(5,4);
t86 = -Ifges(5,2) * t156 - Ifges(5,6) * t274 + t488;
t651 = -mrSges(5,1) * t177 - mrSges(6,1) * t50 + mrSges(7,1) * t30 + mrSges(6,2) * t68 + mrSges(5,3) * t62 - mrSges(7,3) * t46 + t86 / 0.2e1;
t535 = -t156 / 0.2e1;
t650 = -t243 / 0.2e1;
t649 = -t288 / 0.2e1;
t531 = t331 / 0.2e1;
t603 = t357 / 0.2e1;
t596 = Ifges(5,5) + Ifges(7,5);
t648 = Ifges(6,4) - t596;
t597 = Ifges(7,4) + Ifges(6,5);
t647 = -Ifges(5,6) + t597;
t410 = t516 * qJD(3);
t391 = t319 * t410;
t409 = t516 * qJD(4);
t468 = t315 * t316;
t565 = qJD(3) + qJD(4);
t167 = -t319 * t409 + t468 * t565 - t391;
t418 = t316 * t443;
t419 = t516 * t319;
t194 = -t315 * t418 + t419 * t443;
t646 = t167 + t194;
t519 = -t274 / 0.2e1;
t532 = -t331 / 0.2e1;
t534 = t156 / 0.2e1;
t149 = Ifges(7,6) * t331;
t487 = t331 * Ifges(6,6);
t595 = Ifges(7,2) + Ifges(6,3);
t587 = t156 * t595 - t274 * t597 + t149 - t487;
t598 = -Ifges(5,4) + Ifges(7,6);
t645 = -Ifges(5,2) * t535 + Ifges(6,6) * t532 + t647 * t519 + t598 * t531 + t595 * t534 + t587 / 0.2e1 - t651;
t601 = pkin(5) * t331;
t114 = t516 * t336;
t469 = t315 * t125;
t61 = -t114 + t469;
t363 = t61 + t601;
t613 = qJD(5) + t363;
t29 = t274 * t501 + t613;
t49 = pkin(4) * t274 + qJD(5) + t61;
t151 = Ifges(5,4) * t156;
t489 = Ifges(7,6) * t156;
t599 = Ifges(5,1) + Ifges(7,3);
t586 = -t274 * t596 + t331 * t599 - t151 + t489;
t150 = Ifges(6,6) * t156;
t85 = -Ifges(6,4) * t274 - Ifges(6,2) * t331 + t150;
t644 = -mrSges(6,1) * t49 - mrSges(7,1) * t29 - mrSges(5,2) * t177 + mrSges(7,2) * t46 - mrSges(5,3) * t61 + mrSges(6,3) * t68 + t85 / 0.2e1 - t586 / 0.2e1;
t643 = Ifges(3,2) / 0.2e1;
t435 = qJD(1) * qJD(2);
t354 = qJDD(1) * t317 + t320 * t435;
t332 = t316 * qJDD(2) + t319 * t354;
t343 = t357 * qJD(3);
t324 = -t343 + t332;
t333 = t319 * qJDD(2) - t316 * t354;
t325 = qJD(3) * t243 - t333;
t436 = qJD(4) * t315;
t58 = qJD(4) * t215 + t243 * t436 + t315 * t325 - t324 * t516;
t546 = -t58 / 0.2e1;
t545 = t58 / 0.2e1;
t59 = qJD(4) * t331 + t315 * t324 + t325 * t516;
t543 = t59 / 0.2e1;
t250 = qJDD(1) * t320 - t317 * t435;
t240 = qJDD(3) - t250;
t231 = qJDD(4) + t240;
t524 = t231 / 0.2e1;
t594 = Ifges(6,6) - Ifges(7,6);
t66 = t124 * t516 - t469;
t639 = pkin(3) * t409 - t66;
t467 = t315 * t319;
t245 = t316 * t516 + t467;
t168 = t565 * t245;
t193 = t245 * t443;
t638 = t168 - t193;
t439 = qJD(3) * t316;
t573 = -t303 + (-t418 + t439) * pkin(3);
t637 = Ifges(5,4) * t534 - Ifges(6,2) * t531 - t518 * t648 + t599 * t532 - t594 * t535 + t644;
t636 = -Ifges(5,2) * t534 + Ifges(6,6) * t531 + t647 * t518 + t598 * t532 + t595 * t535 + t651;
t499 = mrSges(6,1) * t331;
t129 = -mrSges(6,2) * t274 + t499;
t495 = mrSges(5,3) * t331;
t131 = -mrSges(5,1) * t274 - t495;
t577 = t129 - t131;
t635 = m(6) * t49 + t577;
t381 = mrSges(3,1) * t317 + mrSges(3,2) * t320;
t494 = Ifges(3,4) * t317;
t634 = t317 * (Ifges(3,1) * t320 - t494) / 0.2e1 - pkin(1) * t381;
t633 = mrSges(7,3) + t642;
t313 = qJ(3) + qJ(4);
t306 = sin(t313);
t307 = cos(t313);
t380 = -mrSges(4,1) * t319 + mrSges(4,2) * t316;
t630 = m(4) * pkin(2) + t564 * t306 + t642 * t307 - t380;
t474 = qJDD(1) * pkin(1);
t166 = -t250 * pkin(2) - pkin(8) * t354 - t474;
t238 = t250 * pkin(7);
t208 = qJDD(2) * pkin(8) + t238;
t437 = qJD(3) * t319;
t74 = t316 * t166 + t319 * t208 + t232 * t437 - t262 * t439;
t396 = t319 * t166 - t316 * t208;
t75 = -qJD(3) * t163 + t396;
t629 = -t75 * mrSges(4,1) + t74 * mrSges(4,2);
t34 = t240 * pkin(3) - pkin(9) * t324 - t232 * t439 - t262 * t437 + t396;
t47 = -pkin(9) * t325 + t74;
t10 = -qJD(4) * t335 - t125 * t409 - t315 * t47 + t34 * t516;
t348 = qJDD(5) - t10;
t1 = -t58 * pkin(5) + t274 * qJD(6) - t231 * t501 + t348;
t9 = qJD(4) * t114 - t125 * t436 + t315 * t34 + t516 * t47;
t5 = -qJ(5) * t231 + qJD(5) * t274 - t9;
t3 = -pkin(5) * t59 + qJDD(6) - t5;
t7 = -t231 * pkin(4) + t348;
t628 = -t10 * mrSges(5,1) + t9 * mrSges(5,2) - t7 * mrSges(6,2) - t3 * mrSges(7,2) + t5 * mrSges(6,3) + t1 * mrSges(7,3);
t627 = -Ifges(5,4) * t535 + Ifges(6,2) * t532 + t519 * t648 - t599 * t531 + t594 * t534 + t644;
t373 = t320 * Ifges(3,2) + t494;
t625 = t29 * mrSges(7,3) + t50 * mrSges(6,3) + t61 * mrSges(5,1) + t62 * mrSges(5,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t373 / 0.2e1 + Ifges(4,5) * t650 + Ifges(4,6) * t603 + Ifges(4,3) * t649 + t647 * t535 + t653 + t648 * t531 - t30 * mrSges(7,2) - t49 * mrSges(6,2);
t239 = t354 * pkin(7);
t209 = -qJDD(2) * pkin(2) + t239;
t115 = pkin(3) * t325 + t209;
t11 = t59 * pkin(4) + t58 * qJ(5) - qJD(5) * t331 + t115;
t4 = t59 * qJ(6) + t156 * qJD(6) + t11;
t544 = -t59 / 0.2e1;
t607 = -t231 / 0.2e1;
t624 = mrSges(5,1) * t115 + mrSges(6,1) * t5 - mrSges(7,1) * t3 - mrSges(6,2) * t11 - mrSges(5,3) * t9 + mrSges(7,3) * t4 + Ifges(5,6) * t607 + 0.2e1 * t595 * t543 + t598 * t546 + (t647 + t597) * t524 + (t543 - t544) * Ifges(5,2) + (t594 + t654) * t545;
t623 = mrSges(6,1) * t7 + mrSges(7,1) * t1 + mrSges(5,2) * t115 - mrSges(7,2) * t4 - mrSges(5,3) * t10 - mrSges(6,3) * t11 + Ifges(6,4) * t607 + 0.2e1 * t599 * t546 + t654 * t544 + (-t594 + t598) * t543 + (-t648 + t596) * t524 + (-t545 + t546) * Ifges(6,2);
t420 = qJD(3) * t322;
t248 = t316 * t420;
t263 = t322 * t316;
t264 = t322 * t319;
t112 = -t516 * t248 - t263 * t409 - t264 * t436 - t420 * t467;
t386 = pkin(2) * t317 - pkin(8) * t320;
t247 = t386 * qJD(1);
t173 = pkin(7) * t413 + t319 * t247;
t456 = t319 * t320;
t365 = pkin(3) * t317 - pkin(9) * t456;
t141 = qJD(1) * t365 + t173;
t224 = t316 * t247;
t462 = t317 * t319;
t464 = t316 * t320;
t159 = t224 + (-pkin(7) * t462 - pkin(9) * t464) * qJD(1);
t91 = t315 * t141 + t516 * t159;
t78 = -qJ(5) * t444 - t91;
t617 = t112 - t78;
t583 = -qJD(5) - t639;
t65 = t124 * t315 + t120;
t616 = pkin(3) * t436 - t65;
t478 = qJ(5) * t156;
t172 = t315 * t263 - t264 * t516;
t113 = qJD(4) * t172 + t315 * t248 - t322 * t391;
t383 = -t141 * t516 + t315 * t159;
t615 = t113 - t383;
t614 = qJ(5) * t646 - qJD(5) * t245 + t573;
t440 = qJD(2) * t320;
t416 = t316 * t440;
t350 = t317 * t437 + t416;
t318 = sin(qJ(1));
t321 = cos(qJ(1));
t455 = t320 * t321;
t221 = -t316 * t455 + t318 * t319;
t609 = -m(3) - m(4);
t608 = -m(7) - m(6);
t523 = t240 / 0.2e1;
t605 = t324 / 0.2e1;
t604 = -t325 / 0.2e1;
t602 = pkin(4) * t331;
t502 = qJD(3) / 0.2e1;
t600 = -mrSges(3,3) + mrSges(2,2);
t40 = mrSges(6,1) * t59 - mrSges(6,3) * t231;
t41 = -t59 * mrSges(7,1) + t231 * mrSges(7,2);
t591 = -t40 + t41;
t590 = -pkin(5) * t638 - t617;
t411 = t317 * t501;
t589 = -pkin(5) * t646 + qJD(1) * t411 + t615;
t244 = -t419 + t468;
t588 = qJD(6) * t244 + t501 * t638 + t614;
t585 = pkin(4) * t638 + t614;
t584 = t616 + t612;
t582 = t601 - t583;
t580 = t331 * t501;
t579 = -t62 + t612;
t500 = mrSges(6,1) * t156;
t127 = mrSges(6,3) * t274 + t500;
t498 = mrSges(7,1) * t156;
t128 = -mrSges(7,2) * t274 - t498;
t578 = -t127 + t128;
t496 = mrSges(5,3) * t156;
t130 = mrSges(5,2) * t274 - t496;
t576 = -t130 + t127;
t575 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t357 + t243 * mrSges(4,2) + mrSges(3,3) * t444;
t310 = t320 * pkin(4);
t477 = qJ(5) * t320;
t574 = t306 * t477 + t307 * t310;
t292 = pkin(7) * t456;
t180 = t316 * t253 + t292;
t572 = t243 * mrSges(4,1) - t357 * mrSges(4,2);
t492 = Ifges(4,4) * t243;
t571 = -Ifges(4,1) * t357 - t492;
t389 = -m(7) * t501 - mrSges(7,3);
t470 = t307 * t317;
t471 = t306 * t317;
t570 = (-mrSges(6,2) - t389) * t471 + t641 * t470;
t236 = Ifges(4,4) * t357;
t140 = t243 * Ifges(4,1) + t288 * Ifges(4,5) - t236;
t301 = Ifges(3,4) * t443;
t569 = Ifges(3,1) * t444 + Ifges(3,5) * qJD(2) + t319 * t140 + t301;
t568 = t238 * t320 + t239 * t317;
t567 = -t316 * t75 + t319 * t74;
t382 = t320 * mrSges(3,1) - t317 * mrSges(3,2);
t560 = t317 * mrSges(4,3) + mrSges(2,1) + t382;
t559 = t231 * t563 + t58 * t648 + t59 * t647;
t454 = t321 * t306;
t203 = -t318 * t307 + t320 * t454;
t204 = t306 * t318 + t307 * t455;
t558 = t203 * t633 - t204 * t564;
t458 = t318 * t320;
t201 = t306 * t458 + t307 * t321;
t202 = t307 * t458 - t454;
t557 = t201 * t633 - t202 * t564;
t556 = m(7) * qJ(6) + mrSges(7,3);
t554 = t556 + t642;
t549 = m(5) * pkin(3);
t542 = Ifges(4,1) * t332 / 0.2e1 + Ifges(4,4) * t333 / 0.2e1 + Ifges(4,5) * t523 + t571 * t502;
t522 = t243 / 0.2e1;
t509 = pkin(3) * t243;
t508 = pkin(3) * t315;
t505 = g(3) * t317;
t308 = t317 * pkin(7);
t503 = -qJD(1) / 0.2e1;
t497 = mrSges(7,1) * t331;
t493 = Ifges(3,4) * t320;
t491 = Ifges(4,4) * t316;
t490 = Ifges(4,4) * t319;
t486 = t163 * mrSges(4,3);
t476 = qJ(6) * t201;
t475 = qJ(6) * t203;
t466 = t316 * t317;
t465 = t316 * t318;
t463 = t316 * t321;
t460 = t317 * t322;
t298 = pkin(3) * t319 + pkin(2);
t267 = t320 * t298;
t242 = t319 * t253;
t161 = -pkin(9) * t462 + t242 + (-pkin(7) * t316 - pkin(3)) * t320;
t170 = -pkin(9) * t466 + t180;
t103 = t315 * t161 + t516 * t170;
t249 = t386 * qJD(2);
t442 = qJD(2) * t317;
t427 = pkin(7) * t442;
t447 = t319 * t249 + t316 * t427;
t289 = pkin(3) * t466;
t251 = t308 + t289;
t445 = t321 * pkin(1) + t318 * pkin(7);
t438 = qJD(3) * t317;
t429 = t516 * pkin(3);
t305 = pkin(7) * t440;
t426 = m(4) * pkin(8) + mrSges(4,3);
t425 = t315 * t466;
t423 = t321 * t460;
t422 = Ifges(4,5) * t324 - Ifges(4,6) * t325 + Ifges(4,3) * t240;
t311 = t321 * pkin(7);
t421 = pkin(3) * t463 + t318 * t460 + t311;
t178 = pkin(3) * t350 + t305;
t415 = t316 * t438;
t139 = -Ifges(4,2) * t357 + Ifges(4,6) * t288 + t492;
t412 = -t316 * t139 / 0.2e1;
t42 = -t58 * mrSges(6,1) + t231 * mrSges(6,2);
t39 = -t58 * mrSges(7,1) - t231 * mrSges(7,3);
t400 = -Ifges(4,2) * t243 - t236;
t399 = -t201 * pkin(4) + qJ(5) * t202;
t398 = -t203 * pkin(4) + qJ(5) * t204;
t397 = -qJ(5) * t306 - t298;
t171 = -t516 * t263 - t264 * t315;
t395 = t267 - t460;
t393 = pkin(3) * t465 + t298 * t455 + t445;
t297 = -t429 - pkin(4);
t392 = t516 * t440;
t265 = qJ(5) * t470;
t390 = -pkin(4) * t471 + t265;
t94 = -t103 + t477;
t206 = t317 * t419 - t425;
t384 = -qJ(5) * t206 + t251;
t102 = t161 * t516 - t315 * t170;
t379 = mrSges(4,1) * t316 + mrSges(4,2) * t319;
t377 = -mrSges(5,1) * t306 - mrSges(5,2) * t307;
t375 = Ifges(4,1) * t319 - t491;
t374 = Ifges(4,1) * t316 + t490;
t372 = -Ifges(4,2) * t316 + t490;
t371 = Ifges(4,2) * t319 + t491;
t370 = Ifges(3,5) * t320 - Ifges(3,6) * t317;
t369 = Ifges(4,5) * t319 - t316 * Ifges(4,6);
t368 = Ifges(4,5) * t316 + t319 * Ifges(4,6);
t367 = t509 + t478;
t366 = t221 * pkin(3);
t364 = -qJ(5) * t245 - t298;
t95 = t310 - t102;
t219 = t316 * t458 + t319 * t321;
t359 = t261 * t379;
t101 = t365 * qJD(2) + (-t292 + (pkin(9) * t317 - t253) * t316) * qJD(3) + t447;
t122 = t316 * t249 + t253 * t437 + (-t317 * t441 - t320 * t439) * pkin(7);
t106 = -pkin(9) * t350 + t122;
t26 = t315 * t101 + t516 * t106 + t161 * t409 - t170 * t436;
t355 = -t101 * t516 + t315 * t106 + t161 * t436 + t170 * t409;
t352 = t219 * pkin(3);
t351 = t319 * t440 - t415;
t349 = t204 * pkin(4) + t203 * qJ(5) + t393;
t345 = t357 * mrSges(4,3);
t344 = -g(1) * t203 - g(2) * t201 - g(3) * t471;
t107 = t168 * t317 + t315 * t416 - t319 * t392;
t342 = qJ(5) * t107 - qJD(5) * t206 + t178;
t340 = t366 + t398;
t339 = Ifges(4,5) * t317 + t320 * t375;
t338 = Ifges(4,6) * t317 + t320 * t372;
t337 = Ifges(4,3) * t317 + t320 * t369;
t334 = -t352 + t399;
t24 = -qJ(5) * t442 + qJD(5) * t320 - t26;
t327 = t559 - t628;
t294 = qJ(5) + t508;
t286 = -qJ(6) + t297;
t259 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t443;
t233 = t379 * t317;
t222 = t319 * t455 + t465;
t220 = -t318 * t456 + t463;
t205 = t245 * t317;
t179 = -pkin(7) * t464 + t242;
t176 = mrSges(4,1) * t288 - mrSges(4,3) * t243;
t175 = -t288 * mrSges(4,2) - t345;
t174 = -pkin(7) * t417 + t224;
t152 = pkin(4) * t244 + t364;
t134 = -t244 * pkin(5) + t172;
t133 = pkin(5) * t245 + t171;
t126 = mrSges(7,3) * t274 + t497;
t123 = -qJD(3) * t180 + t447;
t121 = t244 * t501 + t364;
t118 = pkin(4) * t205 + t384;
t117 = -t240 * mrSges(4,2) - mrSges(4,3) * t325;
t116 = t240 * mrSges(4,1) - mrSges(4,3) * t324;
t108 = t316 * t392 - t315 * t415 - qJD(4) * t425 + (t315 * t440 + (t410 + t409) * t317) * t319;
t100 = -mrSges(6,2) * t156 - mrSges(6,3) * t331;
t99 = mrSges(5,1) * t156 + mrSges(5,2) * t331;
t98 = t478 + t602;
t97 = -mrSges(7,2) * t331 + mrSges(7,3) * t156;
t96 = t205 * t501 + t384;
t92 = -t333 * mrSges(4,1) + t332 * mrSges(4,2) + qJD(3) * t572;
t79 = -pkin(4) * t444 + t383;
t77 = t367 + t602;
t72 = Ifges(4,4) * t332 + Ifges(4,2) * t333 + Ifges(4,6) * t240 + qJD(3) * t400;
t69 = -pkin(5) * t205 - t94;
t67 = t206 * pkin(5) + t320 * qJ(6) + t95;
t64 = t478 + t580;
t48 = t367 + t580;
t38 = -mrSges(5,2) * t231 - mrSges(5,3) * t59;
t37 = mrSges(5,1) * t231 + mrSges(5,3) * t58;
t28 = pkin(4) * t108 + t342;
t25 = -pkin(4) * t442 + t355;
t23 = mrSges(7,2) * t58 + mrSges(7,3) * t59;
t22 = mrSges(5,1) * t59 - mrSges(5,2) * t58;
t21 = -mrSges(6,2) * t59 + mrSges(6,3) * t58;
t14 = qJD(6) * t205 + t108 * t501 + t342;
t13 = -pkin(5) * t108 - t24;
t12 = -t107 * pkin(5) - qJD(2) * t411 + t320 * qJD(6) + t355;
t2 = [(t494 + t373) * t250 / 0.2e1 + t568 * mrSges(3,3) + ((pkin(7) * mrSges(3,3) + t643) * t250 - t597 * t543 - t596 * t546 - Ifges(4,6) * t604 - Ifges(4,5) * t605 + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * qJDD(2) + (-Ifges(3,2) * t317 + t493) * t435 / 0.2e1 - Ifges(4,3) * t523 - Ifges(5,6) * t544 - Ifges(6,4) * t545 - t563 * t524 - t559 / 0.2e1 + t628 - t422 / 0.2e1 + Ifges(3,4) * t354 / 0.2e1 + t629) * t320 + t645 * t108 + t634 * t435 + pkin(1) * t250 * mrSges(3,1) + t382 * t474 + (-m(5) * t421 - t220 * mrSges(4,1) - t219 * mrSges(4,2) + t608 * (-t202 * pkin(4) - qJ(5) * t201 + t421) + t600 * t321 + t609 * t311 + t554 * t202 + t564 * t201 + (m(3) * pkin(1) - m(4) * t253 + (m(7) * pkin(5) + mrSges(7,1)) * t317 + (-m(5) + t608) * (-pkin(1) - t267) + t560 - t566) * t318) * g(1) + (-mrSges(3,1) * t308 + Ifges(3,5) * t317) * qJDD(2) - (t139 * t319 + t140 * t316) * t438 / 0.2e1 + t317 * t372 * t604 + t317 * t375 * t605 + t317 * t369 * t523 + t412 * t440 - t72 * t466 / 0.2e1 + (-t162 * t351 - t163 * t350 - t462 * t75 - t466 * t74) * mrSges(4,3) + m(5) * (t10 * t102 + t103 * t9 + t115 * t251 + t177 * t178 + t26 * t62 + t355 * t61) - t355 * t131 + t623 * t206 + t624 * t205 + (t162 * mrSges(4,1) - t163 * mrSges(4,2) + Ifges(6,4) * t532 + Ifges(5,6) * t535 + t563 * t519 + t596 * t531 + t597 * t534 - t625) * t442 + t627 * t107 + qJD(2) ^ 2 * t370 / 0.2e1 + m(4) * (t122 * t163 + t123 * t162 + t179 * t75 + t180 * t74 + (t209 * t317 + t261 * t440) * pkin(7)) + (-m(7) * t349 - m(5) * (t393 - t423) - m(6) * (t349 - t423) - t222 * mrSges(4,1) - t221 * mrSges(4,2) + t609 * t445 + t600 * t318 - t554 * t204 - t564 * t203 + (-m(4) * t387 - t560 + t652) * t321) * g(2) - t357 * (qJD(2) * t338 - t371 * t438) / 0.2e1 + t288 * (qJD(2) * t337 - t368 * t438) / 0.2e1 + t261 * (mrSges(4,1) * t350 + mrSges(4,2) * t351) + t92 * t308 - t259 * t427 + m(7) * (t1 * t67 + t12 * t29 + t13 * t30 + t14 * t46 + t3 * t69 + t4 * t96) + m(6) * (t11 * t118 + t24 * t50 + t25 * t49 + t28 * t68 + t5 * t94 + t7 * t95) + Ifges(2,3) * qJDD(1) + t251 * t22 + t209 * t233 + t178 * t99 + t179 * t116 + t180 * t117 + t122 * t175 + t123 * t176 + t25 * t129 + t26 * t130 + t12 * t126 + t24 * t127 + t13 * t128 + t118 * t21 + t102 * t37 + t103 * t38 + t28 * t100 + t94 * t40 + t95 * t42 + t96 * t23 + t14 * t97 + t67 * t39 + t69 * t41 + (qJD(2) * t339 - t374 * t438) * t522 + t462 * t542 + (t308 * mrSges(3,3) + t493 / 0.2e1 + t317 * Ifges(3,1) - pkin(1) * mrSges(3,2)) * t354 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t568) + t569 * t440 / 0.2e1 + t575 * t305; (-t10 * t171 - t115 * t298 + t172 * t9 + (-t112 - t91) * t62 + t615 * t61 + t573 * t177) * m(5) + (-t587 / 0.2e1 + t636) * t193 + t637 * t194 - (t301 + t569) * t443 / 0.2e1 + (t38 - t40) * t172 + t371 * t604 + t374 * t605 + (t359 + t412) * qJD(3) + t645 * t168 + (t337 * t503 + t369 * t502) * t288 + (t338 * t603 - t162 * (mrSges(4,1) * t317 - mrSges(4,3) * t456) - t163 * (-mrSges(4,2) * t317 - mrSges(4,3) * t464) - t634 * qJD(1)) * qJD(1) + t383 * t131 + (-pkin(2) * t209 - t162 * t173 - t163 * t174 - t261 * t303) * m(4) + (t339 * t503 + t375 * t502) * t243 - t372 * t343 / 0.2e1 + t259 * t302 - t370 * t435 / 0.2e1 + t140 * t437 / 0.2e1 + t209 * t380 + (g(1) * t321 + g(2) * t318) * (t381 + (-t388 - t426 + (m(5) + m(6)) * t322 + t640) * t320 + (m(5) * t298 - m(6) * (-pkin(4) * t307 + t397) - m(7) * t397 - t307 * t389 + t630) * t317) + t623 * t245 + t624 * t244 + t627 * t167 + (t42 - t37) * t171 - t359 * t443 + (-t317 * t426 - m(5) * t395 - t382 - m(6) * (t395 + t574) - m(7) * (t267 + t574) + (-t307 * t556 - t630) * t320 + t652) * g(3) + (t11 * t152 + t171 * t7 - t172 * t5 + t585 * t68 + t617 * t50 + (t113 - t79) * t49) * m(6) - t439 * t486 + (Ifges(6,4) * t531 + Ifges(5,6) * t534 + t443 * t643 + t596 * t532 + t597 * t535 + t625 + t653) * t444 + Ifges(3,5) * t354 + t139 * t418 / 0.2e1 + t319 * t72 / 0.2e1 - t298 * t22 + Ifges(3,6) * t250 + Ifges(3,3) * qJDD(2) - t238 * mrSges(3,2) - t239 * mrSges(3,1) - t174 * t175 - t173 * t176 + t152 * t21 - t79 * t129 - t91 * t130 + t133 * t39 + t134 * t41 - t78 * t127 + t121 * t23 - pkin(2) * t92 + t368 * t523 + t316 * t542 + (-t176 * t437 - t175 * t439 + m(4) * ((-t162 * t319 - t163 * t316) * qJD(3) + t567) - t316 * t116 + t319 * t117) * pkin(8) + (-t162 * t437 + t567) * mrSges(4,3) + t573 * t99 - t575 * t303 + t576 * t112 + t577 * t113 + t585 * t100 + t588 * t97 + t589 * t126 + t590 * t128 + (t1 * t133 + t121 * t4 + t134 * t3 + t29 * t589 + t30 * t590 + t46 * t588) * m(7); -t629 + t635 * t616 + t636 * t331 - t637 * t156 + t587 * t532 + t639 * t130 + (-t345 - t175) * t162 + t422 + t37 * t429 - (-t316 * t549 + t377) * t505 - t99 * t509 - m(5) * (t177 * t509 + t61 * t65 + t62 * t66) + t38 * t508 + t243 * t486 + (t400 + t140) * t603 + t297 * t42 + t286 * t39 + t163 * t176 - t77 * t100 - t48 * t97 + t139 * t522 + (t516 * t10 + t315 * t9 + (t315 * t61 + t516 * t62) * qJD(4)) * t549 + (-Ifges(4,5) * t357 - Ifges(4,6) * t243) * t649 + t571 * t650 + (-m(7) * (t334 - t476) + m(5) * t352 - m(6) * t334 + mrSges(4,1) * t219 - mrSges(4,2) * t220 + t557) * g(2) + (-m(5) * t366 - m(6) * t340 - m(7) * (t340 - t475) - mrSges(4,1) * t221 + mrSges(4,2) * t222 + t558) * g(1) + (-m(6) * (-t289 + t390) - m(7) * (t265 - t289) + t233 + t570) * g(3) - t261 * t572 + t582 * t128 + t583 * t127 + t584 * t126 + (t1 * t286 + t29 * t584 + t294 * t3 + t30 * t582 - t46 * t48) * m(7) + (-t294 * t5 + t297 * t7 + t50 * t583 - t68 * t77) * m(6) + t327 + t591 * t294; (t495 - t635) * t62 + (-g(3) * t265 + qJ(5) * t3 - t1 * t501 + t579 * t29 + t30 * t613 - t46 * t64) * m(7) - t501 * t39 - t177 * (mrSges(5,1) * t331 - mrSges(5,2) * t156) - t68 * (-mrSges(6,2) * t331 + mrSges(6,3) * t156) - t46 * (mrSges(7,2) * t156 + mrSges(7,3) * t331) - t377 * t505 + t363 * t128 + (-pkin(4) * t7 - g(3) * t390 - qJ(5) * t5 - qJD(5) * t50 - t68 * t98) * m(6) + (t156 * t648 + t331 * t647) * t518 + (Ifges(6,2) * t156 + t487 + t86) * t531 + t49 * t500 + t30 * t497 + t29 * t498 - t50 * t499 - t98 * t100 - t64 * t97 + (-m(7) * (t399 - t476) - m(6) * t399 + t557) * g(2) + (-m(6) * t398 - m(7) * (t398 - t475) + t558) * g(1) + t570 * g(3) + (-m(6) * t50 + t496 - t576) * t61 + t578 * qJD(5) + t579 * t126 + t327 + (-Ifges(5,2) * t331 - t151 + t586) * t534 + t591 * qJ(5) + (t331 * t595 + t150 - t489 + t85) * t535 + (-t156 * t599 + t149 - t488 + t587) * t532 - pkin(4) * t42; t578 * t274 + (t100 + t97) * t331 + t39 + t42 + (t274 * t30 + t331 * t46 + t1 + t344) * m(7) + (-t274 * t50 + t331 * t68 + t344 + t7) * m(6); -t274 * t126 - t156 * t97 + (-g(1) * t204 - g(2) * t202 - g(3) * t470 - t156 * t46 - t29 * t274 + t3) * m(7) + t41;];
tau  = t2;
