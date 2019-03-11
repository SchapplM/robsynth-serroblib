% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:32
% EndTime: 2019-03-09 13:51:18
% DurationCPUTime: 28.16s
% Computational Cost: add. (15710->816), mult. (33194->1049), div. (0->0), fcn. (22953->12), ass. (0->393)
t599 = Ifges(4,4) + Ifges(3,5);
t598 = Ifges(4,6) - Ifges(3,6);
t325 = sin(qJ(4));
t330 = cos(qJ(4));
t331 = cos(qJ(2));
t444 = qJD(1) * t331;
t326 = sin(qJ(2));
t445 = qJD(1) * t326;
t209 = -t325 * t444 + t330 * t445;
t531 = t209 / 0.2e1;
t323 = sin(qJ(6));
t611 = t323 / 0.2e1;
t328 = cos(qJ(6));
t376 = Ifges(7,5) * t328 - Ifges(7,6) * t323;
t610 = t376 / 0.2e1;
t501 = Ifges(7,4) * t328;
t379 = -Ifges(7,2) * t323 + t501;
t609 = t379 / 0.2e1;
t502 = Ifges(7,4) * t323;
t382 = Ifges(7,1) * t328 - t502;
t608 = t382 / 0.2e1;
t607 = mrSges(6,2) - mrSges(7,3);
t317 = -qJD(2) + qJD(4);
t307 = qJD(5) + t317;
t207 = -t325 * t445 - t330 * t444;
t324 = sin(qJ(5));
t329 = cos(qJ(5));
t365 = t207 * t324 + t329 * t209;
t103 = t307 * t328 - t323 * t365;
t104 = t307 * t323 + t328 * t365;
t471 = mrSges(6,1) * t307 + mrSges(7,1) * t103 - mrSges(7,2) * t104 - mrSges(6,3) * t365;
t302 = pkin(7) * t445;
t242 = pkin(8) * t445 - t302;
t334 = -pkin(2) - pkin(3);
t414 = t334 * qJD(2);
t178 = qJD(3) + t414 - t242;
t303 = pkin(7) * t444;
t243 = -pkin(8) * t444 + t303;
t320 = qJD(2) * qJ(3);
t210 = t243 + t320;
t113 = t178 * t325 + t210 * t330;
t524 = pkin(9) * t207;
t96 = t113 + t524;
t482 = t324 * t96;
t112 = t330 * t178 - t210 * t325;
t523 = pkin(9) * t209;
t95 = t112 - t523;
t91 = pkin(4) * t317 + t95;
t60 = t329 * t91 - t482;
t58 = -pkin(5) * t307 - t60;
t580 = -m(6) * t60 + m(7) * t58 - t471;
t125 = -t329 * t207 + t324 * t209;
t119 = qJD(6) + t125;
t85 = pkin(5) * t365 + pkin(10) * t125;
t327 = sin(qJ(1));
t454 = t327 * t331;
t280 = qJ(3) * t454;
t294 = pkin(4) * t330 + pkin(3);
t511 = -pkin(2) - t294;
t411 = t326 * t511;
t606 = t327 * t411 + t280;
t605 = -m(5) - m(4);
t436 = qJD(1) * qJD(2);
t246 = -t331 * qJDD(1) + t326 * t436;
t247 = qJDD(1) * t326 + t331 * t436;
t453 = t330 * t331;
t458 = t325 * t326;
t361 = t453 + t458;
t349 = t361 * qJD(4);
t106 = -qJD(1) * t349 + t246 * t325 + t247 * t330;
t316 = -qJDD(2) + qJDD(4);
t225 = t247 * pkin(7);
t403 = qJDD(3) + t225;
t133 = -pkin(8) * t247 + qJDD(2) * t334 + t403;
t224 = t246 * pkin(7);
t177 = qJDD(2) * qJ(3) + qJD(2) * qJD(3) - t224;
t134 = pkin(8) * t246 + t177;
t65 = -qJD(4) * t113 + t330 * t133 - t134 * t325;
t40 = pkin(4) * t316 - pkin(9) * t106 + t65;
t456 = t326 * t330;
t362 = t325 * t331 - t456;
t350 = t362 * qJD(4);
t107 = qJD(1) * t350 + t246 * t330 - t247 * t325;
t440 = qJD(4) * t330;
t441 = qJD(4) * t325;
t64 = t325 * t133 + t330 * t134 + t178 * t440 - t210 * t441;
t42 = pkin(9) * t107 + t64;
t477 = t329 * t96;
t61 = t324 * t91 + t477;
t15 = -qJD(5) * t61 - t324 * t42 + t329 * t40;
t306 = qJDD(5) + t316;
t12 = -pkin(5) * t306 - t15;
t539 = m(7) * t12;
t554 = qJD(6) * t103;
t56 = -qJD(5) * t125 + t106 * t329 + t107 * t324;
t35 = t306 * t323 + t328 * t56 + t554;
t604 = Ifges(7,5) * t35;
t553 = qJD(6) * t104;
t36 = t306 * t328 - t323 * t56 - t553;
t603 = Ifges(7,6) * t36;
t57 = -qJD(5) * t365 - t106 * t324 + t107 * t329;
t55 = qJDD(6) - t57;
t602 = Ifges(7,3) * t55;
t59 = pkin(10) * t307 + t61;
t211 = -qJD(1) * pkin(1) - pkin(2) * t444 - qJ(3) * t445;
t175 = pkin(3) * t444 - t211;
t128 = -pkin(4) * t207 + t175;
t66 = pkin(5) * t125 - pkin(10) * t365 + t128;
t23 = -t323 * t59 + t328 * t66;
t601 = t23 * mrSges(7,1);
t24 = t323 * t66 + t328 * t59;
t600 = t24 * mrSges(7,2);
t16 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t44 = mrSges(6,1) * t306 - mrSges(6,3) * t56;
t597 = t16 - t44;
t596 = t539 + t16;
t298 = Ifges(3,4) * t444;
t499 = Ifges(4,5) * t331;
t383 = t326 * Ifges(4,1) - t499;
t595 = Ifges(3,1) * t445 + qJD(1) * t383 + qJD(2) * t599 + t298;
t231 = t324 * t330 + t325 * t329;
t551 = qJD(4) + qJD(5);
t594 = (qJD(2) - t551) * t231;
t426 = mrSges(4,2) * t445;
t593 = -mrSges(3,3) * t445 - t426 + (mrSges(3,1) + mrSges(4,1)) * qJD(2);
t18 = mrSges(7,1) * t55 - mrSges(7,3) * t35;
t19 = -mrSges(7,2) * t55 + mrSges(7,3) * t36;
t375 = -t323 * t18 + t328 * t19;
t438 = qJD(6) * t328;
t439 = qJD(6) * t323;
t76 = -mrSges(7,2) * t119 + mrSges(7,3) * t103;
t77 = mrSges(7,1) * t119 - mrSges(7,3) * t104;
t592 = -t77 * t438 - t76 * t439 + t375;
t591 = t598 * t326 + t331 * t599;
t388 = t331 * mrSges(4,1) + t326 * mrSges(4,3);
t390 = mrSges(3,1) * t331 - mrSges(3,2) * t326;
t590 = t388 + t390;
t589 = -t224 * t331 + t225 * t326;
t184 = -qJDD(2) * pkin(2) + t403;
t588 = t177 * t331 + t184 * t326;
t262 = -mrSges(7,1) * t328 + mrSges(7,2) * t323;
t492 = t365 * Ifges(6,1);
t116 = Ifges(6,4) * t125;
t384 = mrSges(7,1) * t323 + mrSges(7,2) * t328;
t355 = t58 * t384;
t496 = t104 * Ifges(7,4);
t51 = t103 * Ifges(7,2) + t119 * Ifges(7,6) + t496;
t418 = t51 * t611;
t101 = Ifges(7,4) * t103;
t52 = Ifges(7,1) * t104 + Ifges(7,5) * t119 + t101;
t480 = t328 * t52;
t487 = t307 * Ifges(6,5);
t513 = t60 * mrSges(6,3);
t534 = -t119 / 0.2e1;
t537 = -t104 / 0.2e1;
t538 = -t103 / 0.2e1;
t81 = -t116 + t487 + t492;
t561 = -t487 / 0.2e1 - t128 * mrSges(6,2) + t376 * t534 + t379 * t538 + t382 * t537 - t355 - t480 / 0.2e1 + t418 + t513 - t81 / 0.2e1 + t116 / 0.2e1;
t587 = -t492 / 0.2e1 + t561;
t586 = -mrSges(2,1) - t590;
t109 = -mrSges(6,2) * t307 - mrSges(6,3) * t125;
t359 = -t323 * t77 + t328 * t76 + t109;
t373 = -t23 * t323 + t24 * t328;
t585 = m(6) * t61 + m(7) * t373 + t359;
t486 = t307 * Ifges(6,6);
t494 = t119 * Ifges(7,3);
t495 = t104 * Ifges(7,5);
t497 = t103 * Ifges(7,6);
t50 = t494 + t495 + t497;
t503 = Ifges(6,4) * t365;
t512 = t61 * mrSges(6,3);
t493 = t125 * Ifges(6,2);
t80 = t486 - t493 + t503;
t584 = -t128 * mrSges(6,1) - t601 + t600 + t512 + t80 / 0.2e1 - t50 / 0.2e1 + t503 / 0.2e1 + t486 / 0.2e1;
t583 = m(7) * pkin(10) - t607;
t14 = qJD(5) * t60 + t324 * t40 + t329 * t42;
t11 = pkin(10) * t306 + t14;
t309 = t326 * qJD(3);
t321 = qJDD(1) * pkin(1);
t130 = t246 * pkin(2) - t247 * qJ(3) - qJD(1) * t309 - t321;
t105 = -pkin(3) * t246 - t130;
t73 = -pkin(4) * t107 + t105;
t17 = -pkin(5) * t57 - pkin(10) * t56 + t73;
t2 = qJD(6) * t23 + t11 * t328 + t17 * t323;
t3 = -qJD(6) * t24 - t11 * t323 + t17 * t328;
t581 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t579 = m(5) * pkin(8) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t197 = Ifges(5,4) * t209;
t114 = t207 * Ifges(5,2) + t317 * Ifges(5,6) + t197;
t196 = Ifges(5,4) * t207;
t115 = t209 * Ifges(5,1) + t317 * Ifges(5,5) + t196;
t507 = mrSges(7,3) * t328;
t351 = t15 * mrSges(6,1) - t14 * mrSges(6,2) + Ifges(6,5) * t56 + Ifges(6,6) * t57 + Ifges(6,3) * t306 + t2 * t507;
t578 = t65 * mrSges(5,1) - t64 * mrSges(5,2) + Ifges(5,5) * t106 + Ifges(5,6) * t107 + Ifges(5,3) * t316 + t351 + (t112 * t207 + t113 * t209) * mrSges(5,3) - t317 * (Ifges(5,5) * t207 - Ifges(5,6) * t209) / 0.2e1 - t175 * (mrSges(5,1) * t209 + mrSges(5,2) * t207) - (-Ifges(5,2) * t209 + t115 + t196) * t207 / 0.2e1 + (-Ifges(5,1) * t207 + t114 + t197) * t531;
t8 = Ifges(7,4) * t35 + Ifges(7,2) * t36 + Ifges(7,6) * t55;
t9 = Ifges(7,1) * t35 + Ifges(7,4) * t36 + Ifges(7,5) * t55;
t577 = t12 * t262 + t553 * t608 + t554 * t609 + t55 * (Ifges(7,5) * t323 + Ifges(7,6) * t328) / 0.2e1 + t36 * (Ifges(7,2) * t328 + t502) / 0.2e1 + t35 * (Ifges(7,1) * t323 + t501) / 0.2e1 + t328 * t8 / 0.2e1 + t9 * t611 + (t119 * t610 + t355) * qJD(6);
t572 = -t493 / 0.2e1;
t143 = t330 * t242 + t325 * t243;
t249 = -qJ(3) * t325 + t330 * t334;
t182 = t330 * qJD(3) + qJD(4) * t249;
t564 = t182 - t143;
t142 = -t242 * t325 + t330 * t243;
t250 = t330 * qJ(3) + t325 * t334;
t183 = -t325 * qJD(3) - qJD(4) * t250;
t563 = t183 - t142;
t241 = -pkin(4) + t249;
t145 = t324 * t241 + t329 * t250;
t397 = mrSges(5,1) * t361 - mrSges(5,2) * t362;
t540 = pkin(7) - pkin(8);
t272 = t540 * t326;
t273 = t540 * t331;
t156 = t325 * t272 + t330 * t273;
t332 = cos(qJ(1));
t314 = t332 * pkin(7);
t333 = -pkin(9) - pkin(8);
t527 = pkin(4) * t325;
t560 = (-pkin(1) + t511 * t331 + (-qJ(3) - t527) * t326) * t327 + t332 * t333 + t314;
t451 = qJ(4) + qJ(5);
t412 = sin(t451);
t413 = cos(t451);
t200 = t326 * t413 - t331 * t412;
t558 = t23 * t438 + t24 * t439;
t555 = g(1) * t332 + g(2) * t327;
t552 = qJD(3) + t302;
t171 = t200 * t327;
t199 = t326 * t412 + t331 * t413;
t172 = t199 * t327;
t550 = t171 * mrSges(6,1) - t172 * t607;
t391 = t332 * t412;
t392 = t332 * t413;
t173 = -t326 * t391 - t331 * t392;
t174 = -t326 * t392 + t331 * t391;
t549 = -t174 * mrSges(6,1) + t173 * t607;
t548 = -t199 * mrSges(6,1) - t200 * t607;
t405 = -t439 / 0.2e1;
t416 = t480 / 0.2e1;
t546 = qJD(6) * t416 + t51 * t405 + t577;
t536 = t104 / 0.2e1;
t529 = t326 / 0.2e1;
t528 = pkin(4) * t209;
t526 = pkin(7) * t326;
t525 = pkin(7) * t331;
t522 = pkin(10) * t173;
t521 = pkin(10) * t200;
t518 = g(3) * t200;
t517 = t2 * t328;
t516 = t3 * t323;
t313 = t331 * pkin(2);
t508 = mrSges(7,3) * t323;
t505 = Ifges(3,4) * t326;
t504 = Ifges(3,4) * t331;
t500 = Ifges(4,5) * t326;
t489 = t23 * t328;
t148 = -qJD(2) * t362 + t350;
t149 = qJD(2) * t361 - t349;
t364 = t324 * t362 - t329 * t361;
t71 = qJD(5) * t364 + t148 * t324 + t149 * t329;
t484 = t323 * t71;
t476 = t331 * mrSges(4,3);
t137 = -t324 * t361 - t329 * t362;
t464 = t137 * t323;
t463 = t137 * t328;
t310 = t326 * qJ(3);
t455 = t326 * t332;
t452 = t331 * t332;
t442 = qJD(2) * t331;
t449 = qJ(3) * t442 + t309;
t447 = t313 + t310;
t446 = t332 * pkin(1) + t327 * pkin(7);
t443 = qJD(2) * t326;
t432 = t602 + t603 + t604;
t431 = g(3) * t453;
t285 = pkin(4) * t458;
t428 = m(6) + m(7) - t605;
t427 = t334 * t326;
t425 = mrSges(4,2) * t444;
t420 = t325 * t452;
t419 = t330 * t455;
t415 = t331 * pkin(3) + t447;
t404 = -t438 / 0.2e1;
t402 = -pkin(1) - t310;
t268 = t454 * t527;
t401 = pkin(10) * t172 - t268;
t271 = pkin(4) * t420;
t400 = -t271 - t522;
t399 = -t436 / 0.2e1;
t155 = t330 * t272 - t273 * t325;
t219 = pkin(1) + t415;
t194 = -qJDD(2) * mrSges(4,1) + t247 * mrSges(4,2);
t396 = t331 * t294 + t285 + t447;
t395 = pkin(2) * t452 + qJ(3) * t455 + t446;
t389 = mrSges(3,1) * t326 + mrSges(3,2) * t331;
t190 = t362 * t327;
t191 = t361 * t327;
t387 = t190 * mrSges(5,1) + t191 * mrSges(5,2);
t192 = -t419 + t420;
t193 = t361 * t332;
t386 = t192 * mrSges(5,1) + t193 * mrSges(5,2);
t381 = t331 * Ifges(3,2) + t505;
t374 = t24 * t323 + t489;
t150 = pkin(4) * t361 + t219;
t75 = -pkin(5) * t364 - pkin(10) * t137 + t150;
t117 = pkin(9) * t362 + t155;
t118 = -pkin(9) * t361 + t156;
t83 = t117 * t324 + t118 * t329;
t37 = -t323 * t83 + t328 * t75;
t38 = t323 * t75 + t328 * t83;
t369 = t329 * t117 - t118 * t324;
t166 = -mrSges(5,2) * t317 + mrSges(5,3) * t207;
t167 = mrSges(5,1) * t317 - mrSges(5,3) * t209;
t368 = t166 * t330 - t167 * t325;
t367 = -t172 * t328 - t323 * t332;
t366 = t172 * t323 - t328 * t332;
t144 = t241 * t329 - t250 * t324;
t251 = -qJD(2) * pkin(2) + t552;
t259 = t303 + t320;
t363 = t251 * t331 - t259 * t326;
t228 = t324 * t325 - t329 * t330;
t289 = qJ(3) * t444;
t189 = qJD(1) * t427 + t289;
t360 = t142 + t524;
t358 = pkin(1) * t389;
t357 = t137 * t438 + t484;
t356 = t137 * t439 - t328 * t71;
t354 = t211 * (t326 * mrSges(4,1) - t476);
t353 = t326 * (Ifges(3,1) * t331 - t505);
t352 = t331 * (Ifges(4,3) * t326 + t499);
t168 = t326 * t414 + t449;
t244 = t540 * t443;
t245 = qJD(2) * t273;
t93 = -t330 * t244 + t325 * t245 + t272 * t440 - t273 * t441;
t348 = m(7) * pkin(5) - t262;
t347 = t332 * t285 + t294 * t452 + t327 * t333 + t395;
t346 = t476 + (-m(4) * pkin(2) - mrSges(4,1)) * t326;
t345 = t171 * t262 - t550;
t344 = -t174 * t262 - t549;
t343 = -t199 * t262 - t548;
t132 = t189 - t528;
t99 = -pkin(4) * t148 + t168;
t342 = -qJD(6) * t374 - t516;
t94 = -qJD(4) * t156 + t244 * t325 + t330 * t245;
t340 = t342 + t517;
t339 = (-t323 * t76 - t328 * t77) * qJD(6) + t375;
t338 = -pkin(9) * t149 + t94;
t336 = -t495 / 0.2e1 - t494 / 0.2e1 - t497 / 0.2e1 + t584;
t297 = Ifges(4,5) * t445;
t281 = qJ(3) * t452;
t270 = pkin(4) * t419;
t267 = t327 * pkin(4) * t456;
t261 = qJD(2) * mrSges(4,3) + t425;
t260 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t444;
t256 = -pkin(1) - t447;
t236 = pkin(2) * t445 - t289;
t235 = t388 * qJD(1);
t208 = t228 * qJD(2);
t203 = Ifges(3,6) * qJD(2) + qJD(1) * t381;
t202 = Ifges(4,6) * qJD(2) - Ifges(4,3) * t444 + t297;
t198 = pkin(2) * t443 - t449;
t195 = -mrSges(4,2) * t246 + qJDD(2) * mrSges(4,3);
t188 = t199 * pkin(5);
t165 = t174 * pkin(5);
t164 = t171 * pkin(5);
t152 = -t208 * t328 + t323 * t445;
t151 = t208 * t323 + t328 * t445;
t146 = t551 * t228;
t140 = pkin(5) - t144;
t139 = -t173 * t328 - t323 * t327;
t138 = t173 * t323 - t327 * t328;
t131 = -mrSges(5,1) * t207 + mrSges(5,2) * t209;
t108 = t143 + t523;
t98 = -mrSges(5,2) * t316 + mrSges(5,3) * t107;
t97 = mrSges(5,1) * t316 - mrSges(5,3) * t106;
t84 = mrSges(6,1) * t125 + mrSges(6,2) * t365;
t79 = pkin(9) * t148 + t93;
t74 = t528 + t85;
t72 = qJD(5) * t137 - t329 * t148 + t149 * t324;
t70 = t329 * t108 + t324 * t360;
t67 = t132 - t85;
t63 = t329 * t95 - t482;
t62 = t324 * t95 + t477;
t45 = -mrSges(6,2) * t306 + mrSges(6,3) * t57;
t32 = t323 * t85 + t328 * t60;
t31 = -t323 * t60 + t328 * t85;
t28 = t323 * t67 + t328 * t70;
t27 = -t323 * t70 + t328 * t67;
t26 = t323 * t74 + t328 * t63;
t25 = -t323 * t63 + t328 * t74;
t22 = pkin(5) * t72 - pkin(10) * t71 + t99;
t20 = qJD(5) * t369 + t324 * t338 + t329 * t79;
t5 = -qJD(6) * t38 - t20 * t323 + t22 * t328;
t4 = qJD(6) * t37 + t20 * t328 + t22 * t323;
t1 = [t71 * t416 + t352 * t399 + (-m(5) * (pkin(3) * t452 + t395) - t193 * mrSges(5,1) + t192 * mrSges(5,2) - m(3) * t446 - m(7) * (-pkin(5) * t173 + t347) - t139 * mrSges(7,1) - t138 * mrSges(7,2) - m(6) * t347 + t173 * mrSges(6,1) - m(4) * t395 - t583 * t174 + t586 * t332 + t579 * t327) * g(2) + t331 * (Ifges(3,4) * t247 + Ifges(3,6) * qJDD(2)) / 0.2e1 + t580 * (qJD(5) * t83 + t324 * t79 - t329 * t338) + t365 * (Ifges(6,1) * t71 - Ifges(6,4) * t72) / 0.2e1 + (t15 * m(6) - t539 - t597) * t369 + (t191 * mrSges(5,1) - t190 * mrSges(5,2) - t367 * mrSges(7,1) - t366 * mrSges(7,2) - (-pkin(5) * t172 + t560) * m(7) - m(6) * t560 + t172 * mrSges(6,1) - t583 * t171 + (-m(3) + t605) * t314 + (-m(4) * (t402 - t313) + m(3) * pkin(1) - m(5) * (t331 * t334 + t402) - t586) * t327 + t579 * t332) * g(1) + (-t603 / 0.2e1 - t604 / 0.2e1 - t602 / 0.2e1 + Ifges(6,6) * t306 - t73 * mrSges(6,1) + Ifges(6,4) * t56 + Ifges(6,2) * t57 - t432 / 0.2e1 + t14 * mrSges(6,3) + t581) * t364 - t72 * t600 + (t599 * t326 - t598 * t331) * qJDD(2) / 0.2e1 + ((Ifges(3,1) + Ifges(4,1)) * t247 + t599 * qJDD(2)) * t529 + ((-t261 - t260) * t443 + m(4) * (qJD(2) * t363 + t588) + m(3) * t589 - t593 * t442) * pkin(7) - t125 * (Ifges(6,4) * t71 - Ifges(6,2) * t72) / 0.2e1 + t595 * t442 / 0.2e1 - t331 * (Ifges(4,5) * t247 + Ifges(4,6) * qJDD(2)) / 0.2e1 + (t353 + t331 * (-Ifges(3,2) * t326 + t504) + t326 * (Ifges(4,1) * t331 + t500)) * t436 / 0.2e1 + (t326 * Ifges(3,1) + t383 + t504) * t247 / 0.2e1 - t256 * mrSges(4,3) * t247 - t358 * t436 + (t500 / 0.2e1 - t381 / 0.2e1 - mrSges(3,3) * t525 + t256 * mrSges(4,1) - pkin(1) * mrSges(3,1) + (Ifges(4,5) - Ifges(3,4)) * t529 + (-Ifges(4,3) - Ifges(3,2) / 0.2e1) * t331) * t246 + (t251 * t442 - t259 * t443 + t588) * mrSges(4,2) + (t247 * t526 + t589) * mrSges(3,3) - t72 * t512 - t71 * t513 - t130 * t388 + (m(3) * pkin(1) ^ 2 + Ifges(2,3)) * qJDD(1) + (-t203 / 0.2e1 + t202 / 0.2e1) * t443 + t103 * (-Ifges(7,4) * t356 - Ifges(7,2) * t357 + Ifges(7,6) * t72) / 0.2e1 + (-t2 * t464 + t23 * t356 - t24 * t357 - t3 * t463) * mrSges(7,3) + t119 * (-Ifges(7,5) * t356 - Ifges(7,6) * t357 + Ifges(7,3) * t72) / 0.2e1 + m(4) * (t130 * t256 + t198 * t211) + m(5) * (t105 * t219 + t112 * t94 + t113 * t93 + t155 * t65 + t156 * t64 + t168 * t175) + (-t112 * t149 + t113 * t148) * mrSges(5,3) - t51 * t484 / 0.2e1 + (-pkin(1) * t247 - qJDD(2) * t525) * mrSges(3,2) + t105 * t397 + (t354 + t591 * qJD(2) / 0.2e1) * qJD(2) + t58 * (mrSges(7,1) * t357 - mrSges(7,2) * t356) + (-qJDD(2) * mrSges(3,1) + t194) * t526 - t8 * t464 / 0.2e1 + t9 * t463 / 0.2e1 + t72 * t601 + (-Ifges(7,1) * t356 - Ifges(7,4) * t357 + Ifges(7,5) * t72) * t536 + t195 * t525 + (Ifges(5,1) * t149 + Ifges(5,4) * t148) * t531 + t317 * (Ifges(5,5) * t149 + Ifges(5,6) * t148) / 0.2e1 + t307 * (Ifges(6,5) * t71 - Ifges(6,6) * t72) / 0.2e1 - t198 * t235 + t219 * (-mrSges(5,1) * t107 + mrSges(5,2) * t106) + t207 * (Ifges(5,4) * t149 + Ifges(5,2) * t148) / 0.2e1 + t175 * (-mrSges(5,1) * t148 + mrSges(5,2) * t149) + t93 * t166 + t94 * t167 + t168 * t131 + t150 * (-mrSges(6,1) * t57 + mrSges(6,2) * t56) + t155 * t97 + t156 * t98 + t148 * t114 / 0.2e1 + t149 * t115 / 0.2e1 + m(7) * (t2 * t38 + t23 * t5 + t24 * t4 + t3 * t37) + m(6) * (t128 * t99 + t14 * t83 + t150 * t73 + t20 * t61) + t128 * (mrSges(6,1) * t72 + mrSges(6,2) * t71) + t20 * t109 + t99 * t84 - t72 * t80 / 0.2e1 + t71 * t81 / 0.2e1 + t83 * t45 + t72 * t50 / 0.2e1 + t4 * t76 + t5 * t77 + t37 * t18 + t38 * t19 + (t73 * mrSges(6,2) - t15 * mrSges(6,3) + Ifges(6,1) * t56 + Ifges(6,4) * t57 + Ifges(6,5) * t306 + t12 * t384 + t35 * t608 + t36 * t609 + t51 * t404 + t52 * t405 + t55 * t610) * t137 + (-mrSges(5,3) * t64 - Ifges(5,4) * t106 - Ifges(5,2) * t107 - Ifges(5,6) * t316) * t361 + (mrSges(5,3) * t65 - Ifges(5,1) * t106 - Ifges(5,4) * t107 - Ifges(5,5) * t316) * t362 + t390 * t321; qJD(6) * t418 + t52 * t404 - t578 - (Ifges(4,1) * t444 + t202 + t297) * t445 / 0.2e1 + t598 * t246 + t599 * t247 + t580 * (qJD(5) * t145 + (-t183 + t360) * t329 + (-t108 + t182) * t324) - (-Ifges(3,2) * t445 + t298 + t595) * t444 / 0.2e1 + (t23 * t507 + t24 * t508 + t587) * t125 - (Ifges(7,5) * t537 + Ifges(7,6) * t538 + Ifges(7,3) * t534 + t572 + t584) * t365 + t585 * (qJD(5) * t144 + t182 * t329 + t183 * t324) + t203 * t445 / 0.2e1 + (-m(7) * (t396 - t521) - t348 * t199 - m(4) * t447 - m(6) * t396 - m(5) * t415 - t397 + t548 - t590) * g(3) + t591 * t399 + (m(7) * t340 + t592) * (-pkin(10) + t145) + t593 * t303 + (-m(5) * (t327 * t427 + t280) - t387 - m(7) * (-t401 + t606) + t348 * t171 - m(4) * t280 - t346 * t327 - m(6) * (t268 + t606) + t550) * g(2) + t555 * t389 + t552 * t261 + (-pkin(2) * t184 + qJ(3) * t177 + qJD(3) * t259 - t211 * t236) * m(4) + (-pkin(7) * t363 * m(4) - t354 + (t352 / 0.2e1 - t353 / 0.2e1 + t358) * qJD(1)) * qJD(1) + t563 * t167 + (t112 * t563 + t113 * t564 - t175 * t189 + t249 * t65 + t250 * t64) * m(5) + t564 * t166 + (Ifges(4,2) + Ifges(3,3)) * qJDD(2) - t251 * t425 + t558 * mrSges(7,3) + (-m(5) * (t332 * t427 + t281) - t386 - m(7) * (t455 * t511 + t281 - t400) - t348 * t174 - m(6) * (t332 * t411 + t271 + t281) - m(4) * t281 - t346 * t332 + t549) * g(1) - t577 + (t12 * t140 - t23 * t27 - t24 * t28) * m(7) + (-t128 * t132 + t14 * t145 + t144 * t15 - t61 * t70) * m(6) + t3 * t508 + t249 * t97 + t250 * t98 + t236 * t235 + t224 * mrSges(3,2) - t225 * mrSges(3,1) - pkin(2) * t194 + qJ(3) * t195 - t184 * mrSges(4,1) - t189 * t131 + t177 * mrSges(4,3) + t140 * t16 + t144 * t44 + t145 * t45 - t132 * t84 - t70 * t109 - t28 * t76 - t27 * t77 + t259 * t426 + t260 * t302; (t45 + t339) * t231 + (-t261 - t368) * qJD(2) + t368 * qJD(4) + t194 + t597 * t228 - t359 * t146 + t330 * t97 + t325 * t98 + t428 * t331 * g(3) + t208 * t109 - t151 * t77 - t152 * t76 + ((-t131 - t235 - t84) * qJD(1) - t555 * t428) * t326 + t471 * t594 + (t12 * t228 - t146 * t373 - t151 * t23 - t152 * t24 + t231 * t340 - t58 * t594) * m(7) + (-t128 * t445 + t14 * t231 - t15 * t228 + (-t146 + t208) * t61 + t594 * t60) * m(6) + (-t175 * t445 + t325 * t64 + t330 * t65 + t317 * (-t112 * t325 + t113 * t330)) * m(5) + (-qJD(2) * t259 + t211 * t445 + t184) * m(4); t578 + t596 * (-pkin(4) * t329 - pkin(5)) + t546 + (-m(7) * (-t165 + t270 + t400) - m(6) * (-t271 + t270) + t344 + t386) * g(1) + (-m(7) * (t164 + t267 + t401) - m(6) * (-t268 + t267) + t345 + t387) * g(2) + (m(7) * t431 - t209 * t84 + t324 * t45 + t329 * t44 + (-t128 * t209 + t14 * t324 + t15 * t329 + t431) * m(6) + (t324 * t580 + t585 * t329) * qJD(5)) * pkin(4) + (t572 + t336) * t365 - t587 * t125 + (m(6) * t285 - m(7) * (-t188 - t285 + t521) + t343 + t397) * g(3) + (-t119 * t489 + (-t119 * t24 - t3) * t323) * mrSges(7,3) + t471 * t62 - m(6) * (-t60 * t62 + t61 * t63) + (m(7) * (-t516 + t517 - t558) + t592) * (pkin(4) * t324 + pkin(10)) - t112 * t166 + t113 * t167 - t63 * t109 - t26 * t76 - t25 * t77 - m(7) * (t23 * t25 + t24 * t26 + t58 * t62); -((-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t365 + t374 * mrSges(7,3) + t561) * t125 + t336 * t365 + t546 + ((t165 + t522) * g(1) - t23 * t31 - t24 * t32 - t58 * t61 + t188 * g(3) - t164 * g(2) + (-g(2) * t172 + t340 - t518) * pkin(10)) * m(7) + t342 * mrSges(7,3) + t471 * t61 + t339 * pkin(10) + t344 * g(1) + t343 * g(3) + t345 * g(2) + t351 - t596 * pkin(5) - t60 * t109 - t32 * t76 - t31 * t77; -t58 * (mrSges(7,1) * t104 + mrSges(7,2) * t103) + (Ifges(7,1) * t103 - t496) * t537 + t51 * t536 + (Ifges(7,5) * t103 - Ifges(7,6) * t104) * t534 - t23 * t76 + t24 * t77 - g(1) * (mrSges(7,1) * t138 - mrSges(7,2) * t139) - g(2) * (-mrSges(7,1) * t366 + mrSges(7,2) * t367) + t384 * t518 + (t103 * t23 + t104 * t24) * mrSges(7,3) + t432 + (-Ifges(7,2) * t104 + t101 + t52) * t538 - t581;];
tau  = t1;
