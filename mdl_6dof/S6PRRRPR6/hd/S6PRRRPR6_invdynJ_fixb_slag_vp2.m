% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:23
% EndTime: 2019-03-08 23:32:03
% DurationCPUTime: 27.85s
% Computational Cost: add. (8017->821), mult. (17978->1089), div. (0->0), fcn. (13321->12), ass. (0->392)
t589 = Ifges(5,1) + Ifges(6,1);
t588 = Ifges(6,4) + Ifges(5,5);
t586 = Ifges(5,6) - Ifges(6,6);
t291 = sin(qJ(3));
t469 = Ifges(4,4) * t291;
t295 = cos(qJ(3));
t479 = t295 / 0.2e1;
t591 = Ifges(4,2) * t479 + t469 / 0.2e1;
t292 = sin(qJ(2));
t287 = sin(pkin(6));
t431 = qJD(1) * t287;
t397 = t292 * t431;
t241 = qJD(2) * pkin(8) + t397;
t288 = cos(pkin(6));
t430 = qJD(1) * t288;
t395 = t291 * t430;
t171 = t241 * t295 + t395;
t154 = qJD(3) * pkin(9) + t171;
t478 = pkin(3) * t295;
t249 = -pkin(9) * t291 - pkin(2) - t478;
t296 = cos(qJ(2));
t396 = t296 * t431;
t173 = qJD(2) * t249 - t396;
t290 = sin(qJ(4));
t294 = cos(qJ(4));
t69 = -t290 * t154 + t294 * t173;
t548 = qJD(5) - t69;
t425 = qJD(2) * t295;
t275 = qJD(4) - t425;
t267 = qJD(6) - t275;
t484 = t267 / 0.2e1;
t416 = t294 * qJD(3);
t427 = qJD(2) * t291;
t230 = t290 * t427 - t416;
t231 = qJD(3) * t290 + t294 * t427;
t289 = sin(qJ(6));
t293 = cos(qJ(6));
t332 = t230 * t289 + t231 * t293;
t492 = t332 / 0.2e1;
t125 = t230 * t293 - t231 * t289;
t494 = t125 / 0.2e1;
t594 = -Ifges(7,5) * t492 - Ifges(7,6) * t494 - Ifges(7,3) * t484;
t415 = qJD(2) * qJD(3);
t240 = qJDD(2) * t291 + t295 * t415;
t422 = qJD(4) * t230;
t117 = qJDD(3) * t290 + t240 * t294 - t422;
t498 = t117 / 0.2e1;
t118 = qJD(4) * t231 - t294 * qJDD(3) + t240 * t290;
t496 = t118 / 0.2e1;
t239 = t295 * qJDD(2) - t291 * t415;
t228 = qJDD(4) - t239;
t490 = t228 / 0.2e1;
t488 = t230 / 0.2e1;
t593 = -t231 / 0.2e1;
t592 = -t275 / 0.2e1;
t590 = qJD(3) / 0.2e1;
t587 = Ifges(6,2) + Ifges(5,3);
t363 = t295 * t396;
t164 = t290 * t363 - t294 * t397;
t357 = pkin(3) * t291 - pkin(9) * t295;
t238 = t357 * qJD(3);
t437 = t294 * t295;
t277 = pkin(8) * t437;
t421 = qJD(4) * t290;
t356 = qJD(4) * t277 - t238 * t294 + t249 * t421;
t420 = qJD(4) * t291;
t391 = t290 * t420;
t398 = -pkin(8) * t290 - pkin(4);
t410 = pkin(10) * t437;
t585 = -t164 + pkin(10) * t391 + (-t410 + (-pkin(5) + t398) * t291) * qJD(3) + t356;
t436 = t295 * t296;
t176 = (t290 * t292 + t294 * t436) * t287;
t165 = qJD(1) * t176;
t424 = qJD(3) * t291;
t280 = qJ(5) * t424;
t419 = qJD(4) * t294;
t433 = t290 * t238 + t249 * t419;
t438 = t291 * t294;
t584 = t165 - t280 - (-pkin(8) * qJD(3) + pkin(10) * qJD(4)) * t438 - (-qJD(5) + (-pkin(8) * qJD(4) + pkin(10) * qJD(3)) * t290) * t295 - t433;
t392 = t290 * t425;
t499 = pkin(9) - pkin(10);
t170 = -t291 * t241 + t295 * t430;
t235 = t357 * qJD(2);
t103 = t294 * t170 + t290 * t235;
t86 = qJ(5) * t427 + t103;
t583 = pkin(10) * t392 + t499 * t421 + t86;
t102 = -t290 * t170 + t235 * t294;
t258 = t499 * t294;
t500 = pkin(4) + pkin(5);
t582 = qJD(4) * t258 - (-t291 * t500 - t410) * qJD(2) + t102;
t581 = -pkin(10) * t231 + t548;
t532 = -t290 * t586 + t294 * t588;
t464 = Ifges(6,5) * t290;
t467 = Ifges(5,4) * t290;
t530 = t294 * t589 + t464 - t467;
t450 = qJ(5) * t290;
t336 = pkin(4) * t294 + t450;
t119 = Ifges(7,4) * t125;
t580 = Ifges(7,2) * t332 - t119;
t34 = -t275 * t500 + t581;
t266 = t275 * qJ(5);
t70 = t294 * t154 + t290 * t173;
t48 = pkin(10) * t230 + t70;
t38 = t266 + t48;
t15 = -t289 * t38 + t293 * t34;
t16 = t289 * t34 + t293 * t38;
t221 = qJDD(6) - t228;
t28 = qJD(6) * t125 + t117 * t293 + t118 * t289;
t29 = -qJD(6) * t332 - t117 * t289 + t118 * t293;
t409 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t221;
t465 = Ifges(7,4) * t332;
t485 = -t267 / 0.2e1;
t493 = -t332 / 0.2e1;
t358 = qJD(3) * pkin(3) + t170;
t315 = qJ(5) * t231 + t358;
t50 = -t230 * t500 + t315;
t413 = qJDD(1) * t288;
t428 = qJD(2) * t287;
t388 = qJD(1) * t428;
t260 = t296 * t388;
t414 = qJDD(1) * t287;
t201 = t292 * t414 + t260;
t576 = qJDD(2) * pkin(8) + qJD(3) * t430 + t201;
t65 = -t241 * t424 + t291 * t413 + t295 * t576;
t59 = qJDD(3) * pkin(9) + t65;
t259 = t292 * t388;
t200 = t296 * t414 - t259;
t177 = -qJDD(2) * pkin(2) - t200;
t89 = -pkin(3) * t239 - pkin(9) * t240 + t177;
t18 = -t154 * t419 - t173 * t421 - t290 * t59 + t294 * t89;
t319 = qJDD(5) - t18;
t5 = -pkin(10) * t117 - t228 * t500 + t319;
t17 = -t154 * t421 + t173 * t419 + t290 * t89 + t294 * t59;
t12 = t228 * qJ(5) + t275 * qJD(5) + t17;
t6 = pkin(10) * t118 + t12;
t1 = qJD(6) * t15 + t289 * t5 + t293 * t6;
t2 = -qJD(6) * t16 - t289 * t6 + t293 * t5;
t519 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t579 = (t125 * t15 + t16 * t332) * mrSges(7,3) + (Ifges(7,5) * t125 - Ifges(7,6) * t332) * t485 + (Ifges(7,1) * t125 - t465) * t493 - t50 * (mrSges(7,1) * t332 + mrSges(7,2) * t125) + t409 + t519;
t497 = -t118 / 0.2e1;
t578 = t588 * t490 + (-Ifges(5,4) + Ifges(6,5)) * t496 + t589 * t498;
t227 = Ifges(5,4) * t230;
t458 = t230 * Ifges(6,5);
t577 = t231 * t589 + t275 * t588 - t227 + t458;
t575 = -qJD(5) * t290 - t395;
t508 = m(7) * pkin(5);
t512 = -mrSges(7,1) * t293 + mrSges(7,2) * t289 - mrSges(5,1) - mrSges(6,1) - t508;
t570 = mrSges(7,1) * t15 - mrSges(7,2) * t16 + Ifges(4,6) * t590 + qJD(2) * t591 + t488 * t586 + t587 * t592 + t588 * t593 - t594;
t330 = t289 * t290 + t293 * t294;
t331 = t289 * t294 - t290 * t293;
t569 = mrSges(7,1) * t330 - mrSges(7,2) * t331;
t507 = t28 / 0.2e1;
t506 = t29 / 0.2e1;
t501 = -m(7) - m(6);
t491 = t221 / 0.2e1;
t439 = t290 * t295;
t276 = pkin(8) * t439;
t285 = t295 * pkin(4);
t116 = pkin(5) * t295 + t276 + t285 + (-pkin(10) * t291 - t249) * t294;
t180 = t290 * t249 + t277;
t160 = -qJ(5) * t295 + t180;
t440 = t290 * t291;
t128 = pkin(10) * t440 + t160;
t40 = t116 * t289 + t128 * t293;
t568 = -qJD(6) * t40 + t289 * t584 + t293 * t585;
t39 = t116 * t293 - t128 * t289;
t567 = qJD(6) * t39 + t289 * t585 - t293 * t584;
t566 = mrSges(3,2) - mrSges(4,3);
t257 = t499 * t290;
t158 = t257 * t293 - t258 * t289;
t563 = qJD(6) * t158 + t289 * t582 - t293 * t583;
t159 = t257 * t289 + t258 * t293;
t562 = -qJD(6) * t159 + t289 * t583 + t293 * t582;
t80 = mrSges(5,1) * t228 - mrSges(5,3) * t117;
t81 = -t228 * mrSges(6,1) + t117 * mrSges(6,2);
t561 = t81 - t80;
t82 = -mrSges(5,2) * t228 - mrSges(5,3) * t118;
t83 = -mrSges(6,2) * t118 + mrSges(6,3) * t228;
t560 = t83 + t82;
t555 = -m(7) * t499 + mrSges(4,2) + mrSges(7,3);
t136 = mrSges(6,1) * t230 - mrSges(6,3) * t231;
t49 = -mrSges(7,1) * t125 + mrSges(7,2) * t332;
t453 = t136 - t49;
t244 = qJ(5) * t293 - t289 * t500;
t554 = -qJD(6) * t244 - t289 * t581 - t293 * t48;
t42 = mrSges(5,1) * t118 + mrSges(5,2) * t117;
t553 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t240 + t42;
t243 = -qJ(5) * t289 - t293 * t500;
t552 = qJD(6) * t243 - t289 * t48 + t293 * t581;
t355 = mrSges(4,1) * t295 - mrSges(4,2) * t291;
t550 = t355 + mrSges(3,1);
t449 = qJ(5) * t294;
t326 = -t290 * t500 + t449;
t549 = -(qJD(2) * t326 - t241) * t295 + qJD(4) * t326 - t575;
t196 = t331 * t291;
t335 = pkin(4) * t290 - t449;
t547 = -(qJD(2) * t335 + t241) * t295 + qJD(4) * t335 + t575;
t452 = cos(pkin(11));
t368 = t452 * t292;
t286 = sin(pkin(11));
t445 = t286 * t296;
t208 = t288 * t368 + t445;
t369 = t287 * t452;
t138 = -t208 * t291 - t295 * t369;
t545 = t336 * t138;
t367 = t452 * t296;
t446 = t286 * t292;
t210 = -t288 * t446 + t367;
t443 = t287 * t295;
t140 = -t210 * t291 + t286 * t443;
t544 = t336 * t140;
t523 = qJD(4) - qJD(6);
t133 = t523 * t330;
t317 = t330 * t295;
t182 = qJD(2) * t317;
t543 = t133 - t182;
t134 = t523 * t331;
t318 = t331 * t295;
t181 = qJD(2) * t318;
t542 = -t134 + t181;
t404 = mrSges(4,3) * t427;
t541 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t230 + mrSges(5,2) * t231 + t404;
t470 = mrSges(5,3) * t231;
t167 = mrSges(5,1) * t275 - t470;
t168 = -mrSges(6,1) * t275 + mrSges(6,2) * t231;
t540 = t167 - t168;
t471 = mrSges(5,3) * t230;
t166 = -mrSges(5,2) * t275 - t471;
t169 = -mrSges(6,2) * t230 + mrSges(6,3) * t275;
t434 = t169 + t166;
t444 = t287 * t292;
t212 = -t288 * t295 + t291 * t444;
t539 = t336 * t212;
t226 = Ifges(6,5) * t231;
t105 = t275 * Ifges(6,6) + t230 * Ifges(6,3) + t226;
t281 = Ifges(4,4) * t425;
t538 = Ifges(4,1) * t427 + Ifges(4,5) * qJD(3) + t290 * t105 + t281;
t537 = t291 * t587 + t295 * t532;
t536 = t291 * t588 + t295 * t530;
t535 = -pkin(4) * t501 - t512;
t351 = mrSges(6,1) * t290 - mrSges(6,3) * t294;
t353 = mrSges(5,1) * t290 + mrSges(5,2) * t294;
t73 = pkin(4) * t230 - t315;
t534 = -t73 * t351 + t353 * t358;
t533 = t290 * t588 + t294 * t586;
t463 = Ifges(6,5) * t294;
t466 = Ifges(5,4) * t294;
t531 = t290 * t589 - t463 + t466;
t527 = t117 * t588 - t118 * t586 + t228 * t587;
t423 = qJD(3) * t295;
t66 = -t241 * t423 - t291 * t576 + t295 * t413;
t526 = -t291 * t66 + t295 * t65;
t525 = t17 * t294 - t18 * t290;
t13 = -pkin(4) * t228 + t319;
t524 = t12 * t294 + t13 * t290;
t522 = -t117 * Ifges(6,5) / 0.2e1 - t228 * Ifges(6,6) / 0.2e1 + Ifges(5,4) * t498 + Ifges(5,6) * t490 + (Ifges(6,3) + Ifges(5,2)) * t497;
t521 = -t294 * t500 - t450;
t352 = mrSges(6,1) * t294 + mrSges(6,3) * t290;
t354 = mrSges(5,1) * t294 - mrSges(5,2) * t290;
t520 = -t294 * t508 - mrSges(4,1) - t352 - t354 - t569;
t518 = m(5) * t358 - t541;
t516 = m(7) * pkin(10) - mrSges(6,2) - mrSges(5,3) + mrSges(7,3);
t514 = -t289 * mrSges(7,1) - t293 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t139 = t208 * t295 - t291 * t369;
t141 = t286 * t287 * t291 + t210 * t295;
t213 = t288 * t291 + t292 * t443;
t513 = -g(1) * t141 - g(2) * t139 - g(3) * t213;
t511 = -t18 * mrSges(5,1) + t13 * mrSges(6,1) + t17 * mrSges(5,2) - t12 * mrSges(6,3);
t298 = qJD(2) ^ 2;
t510 = Ifges(7,4) * t507 + Ifges(7,2) * t506 + Ifges(7,6) * t491;
t509 = Ifges(7,1) * t507 + Ifges(7,4) * t506 + Ifges(7,5) * t491;
t36 = Ifges(7,2) * t125 + Ifges(7,6) * t267 + t465;
t505 = -t36 / 0.2e1;
t504 = t36 / 0.2e1;
t37 = Ifges(7,1) * t332 + Ifges(7,5) * t267 + t119;
t503 = -t37 / 0.2e1;
t502 = t37 / 0.2e1;
t495 = -t125 / 0.2e1;
t489 = -t230 / 0.2e1;
t486 = t231 / 0.2e1;
t482 = t275 / 0.2e1;
t473 = -qJD(2) / 0.2e1;
t472 = qJD(4) / 0.2e1;
t468 = Ifges(4,4) * t295;
t457 = t231 * Ifges(5,4);
t320 = qJDD(3) * pkin(3) + t66;
t456 = t291 * t320;
t451 = qJ(5) * t230;
t207 = -t288 * t367 + t446;
t448 = t207 * t291;
t209 = t288 * t445 + t368;
t447 = t209 * t291;
t442 = t287 * t296;
t432 = pkin(2) * t442 + pkin(8) * t444;
t426 = qJD(2) * t292;
t417 = qJD(5) * t294;
t408 = pkin(8) * t424;
t406 = pkin(9) * t421;
t405 = pkin(9) * t419;
t403 = mrSges(4,3) * t425;
t402 = t291 * t442;
t401 = t290 * t442;
t394 = t287 * t426;
t393 = t296 * t428;
t108 = -t230 * Ifges(5,2) + t275 * Ifges(5,6) + t457;
t390 = -t290 * t108 / 0.2e1;
t382 = -t425 / 0.2e1;
t379 = t423 / 0.2e1;
t376 = -t420 / 0.2e1;
t375 = t419 / 0.2e1;
t374 = -t207 * pkin(2) + pkin(8) * t208;
t373 = -t209 * pkin(2) + pkin(8) * t210;
t129 = t138 * pkin(3);
t372 = pkin(9) * t139 + t129;
t130 = t140 * pkin(3);
t371 = pkin(9) * t141 + t130;
t204 = t212 * pkin(3);
t370 = pkin(9) * t213 - t204;
t366 = t415 / 0.2e1;
t179 = t249 * t294 - t276;
t365 = t287 * pkin(3) * t436 + pkin(9) * t402 + t432;
t345 = -Ifges(5,2) * t290 + t466;
t344 = Ifges(5,2) * t294 + t467;
t341 = Ifges(4,5) * t295 - Ifges(4,6) * t291;
t338 = Ifges(6,3) * t290 + t463;
t337 = -Ifges(6,3) * t294 + t464;
t93 = -mrSges(7,2) * t267 + mrSges(7,3) * t125;
t94 = mrSges(7,1) * t267 - mrSges(7,3) * t332;
t333 = -t289 * t94 + t293 * t93;
t146 = t213 * t290 + t294 * t442;
t147 = t213 * t294 - t401;
t52 = t146 * t293 - t147 * t289;
t53 = t146 * t289 + t147 * t293;
t329 = -pkin(9) * t448 - t207 * t478 + t374;
t328 = -pkin(9) * t447 - t209 * t478 + t373;
t327 = pkin(8) + t335;
t242 = -qJD(2) * pkin(2) - t396;
t323 = t242 * (mrSges(4,1) * t291 + mrSges(4,2) * t295);
t322 = t291 * (Ifges(4,1) * t295 - t469);
t76 = t139 * t290 - t207 * t294;
t78 = t141 * t290 - t209 * t294;
t321 = -g(1) * t78 - g(2) * t76 - g(3) * t146;
t316 = -pkin(8) + t326;
t313 = t295 * t416 - t391;
t312 = t290 * t423 + t291 * t419;
t304 = Ifges(5,6) * t291 + t295 * t345;
t303 = Ifges(6,6) * t291 + t295 * t338;
t301 = qJ(5) * t501 + t514;
t91 = (-t291 * t416 - t295 * t421) * pkin(8) + t433;
t300 = qJ(5) * t117 + qJD(5) * t231 + t320;
t253 = -qJD(3) * mrSges(4,2) + t403;
t245 = -pkin(3) - t336;
t234 = t355 * qJD(2);
t225 = pkin(3) - t521;
t202 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t239;
t197 = t330 * t291;
t193 = t327 * t291;
t175 = -t294 * t444 + t295 * t401;
t161 = -t179 + t285;
t157 = t316 * t291;
t148 = -mrSges(4,1) * t239 + mrSges(4,2) * t240;
t145 = -qJD(3) * t212 + t295 * t393;
t144 = qJD(3) * t213 + t291 * t393;
t135 = pkin(4) * t231 + t451;
t101 = -t209 * t437 + t210 * t290;
t100 = -t209 * t439 - t210 * t294;
t99 = -t207 * t437 + t208 * t290;
t98 = -t207 * t439 - t208 * t294;
t95 = -t231 * t500 - t451;
t92 = t290 * t408 - t356;
t90 = (qJD(4) * t336 - t417) * t291 + t327 * t423;
t88 = -pkin(4) * t427 - t102;
t84 = t398 * t424 + t356;
t79 = t141 * t294 + t209 * t290;
t77 = t139 * t294 + t207 * t290;
t67 = -qJD(5) * t295 + t280 + t91;
t62 = -qJD(3) * t318 + t133 * t291;
t61 = qJD(3) * t317 + t196 * t523;
t57 = (qJD(4) * t521 + t417) * t291 + t316 * t423;
t54 = t266 + t70;
t51 = -pkin(4) * t275 + t548;
t45 = -qJD(4) * t146 + t145 * t294 + t290 * t394;
t44 = -qJD(4) * t401 + t145 * t290 + t213 * t419 - t294 * t394;
t41 = mrSges(6,1) * t118 - mrSges(6,3) * t117;
t23 = -mrSges(7,2) * t221 + mrSges(7,3) * t29;
t22 = mrSges(7,1) * t221 - mrSges(7,3) * t28;
t19 = pkin(4) * t118 - t300;
t14 = -t118 * t500 + t300;
t11 = qJD(6) * t52 + t289 * t44 + t293 * t45;
t10 = -qJD(6) * t53 - t289 * t45 + t293 * t44;
t9 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t3 = [m(2) * qJDD(1) + t10 * t94 + t11 * t93 + t145 * t253 + t213 * t202 + t52 * t22 + t53 * t23 + t434 * t45 - t540 * t44 + t560 * t147 + t561 * t146 + (t41 - t9 + t553) * t212 + (t453 + t541) * t144 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t298 - t148) * t296 + (-mrSges(3,1) * t298 - mrSges(3,2) * qJDD(2) - qJD(2) * t234) * t292) * t287 + (-m(2) - m(3) - m(4) - m(5) + t501) * g(3) + m(4) * (-t144 * t170 + t145 * t171 - t212 * t66 + t213 * t65 + (-t177 * t296 + t242 * t426) * t287) + m(3) * (qJDD(1) * t288 ^ 2 + (t200 * t296 + t201 * t292) * t287) + m(6) * (t12 * t147 + t13 * t146 + t144 * t73 + t19 * t212 + t44 * t51 + t45 * t54) + m(5) * (-t144 * t358 - t146 * t18 + t147 * t17 - t212 * t320 - t44 * t69 + t45 * t70) + m(7) * (t1 * t53 + t10 * t15 + t11 * t16 - t14 * t212 - t144 * t50 + t2 * t52); (Ifges(7,5) * t61 + Ifges(7,6) * t62) * t484 + (t303 * t488 + t304 * t489 + t341 * t590 + t537 * t482 + t536 * t486 + t323) * qJD(3) + (t17 * t180 + t179 * t18 + (t91 - t165) * t70 + (t92 + t164) * t69) * m(5) + (-(t242 * t292 + (-t170 * t291 + t171 * t295) * t296) * t431 - pkin(2) * t177) * m(4) + (t13 * t438 - t312 * t54 + t313 * t51) * mrSges(6,2) + (-t170 * t423 + t526) * mrSges(4,3) + (-t18 * t438 - t312 * t70 - t313 * t69) * mrSges(5,3) + (t260 - t201) * mrSges(3,2) + t438 * t578 + ((-t358 * t423 - t456) * m(5) + t541 * t423 + t553 * t291 + ((-t170 * t295 - t171 * t291) * qJD(3) + t526) * m(4) + t202 * t295) * pkin(8) - t358 * (mrSges(5,1) * t312 + mrSges(5,2) * t313) + (t291 * Ifges(4,1) + t468 / 0.2e1 + Ifges(4,4) * t479) * t240 - t177 * t355 + (Ifges(7,4) * t61 + Ifges(7,2) * t62) * t494 + t73 * (mrSges(6,1) * t312 - mrSges(6,3) * t313) + (mrSges(5,1) * t69 - mrSges(6,1) * t51 - mrSges(5,2) * t70 - t171 * mrSges(4,3) + mrSges(6,3) * t54 - t570 + t594) * t424 + (t12 * t160 + t13 * t161 + t19 * t193 + t73 * t90 + (t67 - t165) * t54 + (t84 - t164) * t51) * m(6) + (t259 + t200) * mrSges(3,1) + t294 * t108 * t376 + (-Ifges(6,6) * t496 - Ifges(5,6) * t497 + Ifges(7,3) * t491 + (-Ifges(4,2) * t291 + t468) * t366 + Ifges(7,6) * t506 + Ifges(7,5) * t507 - t588 * t498 - t587 * t490 + t511 + t519) * t295 + (-t363 - t408) * t253 + (-t1 * t196 - t15 * t61 + t16 * t62 - t197 * t2) * mrSges(7,3) + (Ifges(7,5) * t197 - Ifges(7,6) * t196) * t491 + (Ifges(7,4) * t197 - Ifges(7,2) * t196) * t506 + t14 * (mrSges(7,1) * t196 + mrSges(7,2) * t197) + (Ifges(7,1) * t197 - Ifges(7,4) * t196) * t507 + (-m(6) * t73 + m(7) * t50 - t453 + t518) * t291 * t396 + t90 * t136 + t50 * (-mrSges(7,1) * t62 + mrSges(7,2) * t61) + t57 * t49 + t39 * t22 + t40 * t23 + Ifges(3,3) * qJDD(2) + t61 * t502 + t62 * t504 + t92 * t167 + t84 * t168 + t67 * t169 + t179 * t80 + t180 * t82 + t193 * t41 + t234 * t397 + t577 * (t290 * t376 + t294 * t379) + (-t12 * mrSges(6,2) - t17 * mrSges(5,3) - t522) * t440 - t527 * t295 / 0.2e1 + 0.2e1 * t591 * t239 + t322 * t366 + t197 * t509 - t196 * t510 + (t105 * t375 + t19 * t351 + t338 * t496 + t345 * t497 + t532 * t490 + t530 * t498) * t291 + (-t337 * t488 - t344 * t489 - t482 * t533 - t486 * t531) * t420 + t538 * t379 - t434 * t165 + t540 * t164 + (Ifges(7,1) * t61 + Ifges(7,4) * t62) * t492 + t91 * t166 - pkin(2) * t148 + t157 * t9 + t160 * t83 + t161 * t81 + t567 * t93 + t568 * t94 + (t1 * t40 + t14 * t157 + t15 * t568 + t16 * t567 + t2 * t39 + t50 * t57) * m(7) + (Ifges(4,5) * t291 + 0.2e1 * Ifges(4,6) * t479) * qJDD(3) + (-m(4) * t373 - m(5) * t328 + t501 * (t101 * pkin(4) + qJ(5) * t100 + t328) + t566 * t210 + t550 * t209 + t512 * t101 + t514 * t100 - t516 * t447) * g(1) + (-m(4) * t374 - m(5) * t329 + t501 * (t99 * pkin(4) + qJ(5) * t98 + t329) + t512 * t99 + t514 * t98 + t566 * t208 + t550 * t207 - t516 * t448) * g(2) + (-m(4) * t432 - m(5) * t365 + t501 * (t176 * pkin(4) + qJ(5) * t175 + t365) + (t292 * t566 - t296 * t550) * t287 + t512 * t176 + t514 * t175 + t516 * t402) * g(3) + t390 * t423 - t353 * t456 + t409 * t479; t14 * t569 + (t472 * t532 + t473 * t537) * t275 + t290 * t578 + (t472 * t530 + t473 * t536) * t231 + (Ifges(7,1) * t133 - Ifges(7,4) * t134) * t492 + (Ifges(7,4) * t133 - Ifges(7,2) * t134) * t494 + (Ifges(7,5) * t133 - Ifges(7,6) * t134) * t484 + (Ifges(7,5) * t182 - Ifges(7,6) * t181) * t485 + (Ifges(7,4) * t182 - Ifges(7,2) * t181) * t495 + (Ifges(7,1) * t182 - Ifges(7,4) * t181) * t493 + t320 * t354 + (pkin(3) * t320 - t102 * t69 - t103 * t70) * m(5) - t330 * t510 + (-Ifges(7,5) * t331 - Ifges(7,6) * t330) * t491 + (-Ifges(7,4) * t331 - Ifges(7,2) * t330) * t506 + (-Ifges(7,1) * t331 - Ifges(7,4) * t330) * t507 + (-t1 * t330 - t15 * t543 + t16 * t542 + t2 * t331) * mrSges(7,3) - t331 * t509 - t19 * t352 + (t403 - t253) * t170 + t337 * t496 + t344 * t497 + (-t345 / 0.2e1 + t338 / 0.2e1) * t422 + (-t102 - t405) * t167 + (t19 * t245 - t51 * t88 - t54 * t86 + t547 * t73) * m(6) + (t560 * t294 + t561 * t290 + ((-t290 * t70 - t294 * t69) * qJD(4) + t525) * m(5) + ((-t290 * t54 + t294 * t51) * qJD(4) + t524) * m(6)) * pkin(9) + (t281 + t538) * t382 + (t405 - t88) * t168 + (-m(5) * t370 - m(6) * (t370 - t539) - m(7) * (-t204 - t539) + t555 * t213 - t520 * t212) * g(3) + (-t323 - t69 * (mrSges(5,1) * t291 - mrSges(5,3) * t437) - t51 * (-mrSges(6,1) * t291 + mrSges(6,2) * t437) - t70 * (-mrSges(5,2) * t291 - mrSges(5,3) * t439) - t54 * (-mrSges(6,2) * t439 + mrSges(6,3) * t291) + (t304 / 0.2e1 - t303 / 0.2e1) * t230) * qJD(2) + (-t86 - t406) * t169 + (-t406 - t103) * t166 - t65 * mrSges(4,2) + t66 * mrSges(4,1) - pkin(3) * t42 + Ifges(4,3) * qJDD(3) - t298 * t322 / 0.2e1 + t133 * t502 + t182 * t503 - t134 * t504 - t181 * t505 + (-Ifges(7,5) * t493 - Ifges(4,2) * t382 - Ifges(7,6) * t495 - Ifges(7,3) * t485 + t570) * t427 + t108 * t392 / 0.2e1 + (t404 + t518) * t171 + t225 * t9 + t577 * (t294 * t382 + t375) + t522 * t294 + Ifges(4,6) * t239 + (t419 * t51 - t421 * t54 + t513 + t524) * mrSges(6,2) + Ifges(4,5) * t240 + (-t419 * t69 - t421 * t70 + t513 + t525) * mrSges(5,3) + t245 * t41 + t531 * t498 + t533 * t490 + (t390 - t534) * qJD(4) + t534 * t425 + (-mrSges(7,1) * t542 + mrSges(7,2) * t543) * t50 + t547 * t136 + t549 * t49 + (-m(5) * t371 - m(6) * (t371 + t544) - m(7) * (t130 + t544) + t555 * t141 + t520 * t140) * g(1) + (-m(6) * (t372 + t545) - m(5) * t372 - m(7) * (t129 + t545) + t555 * t139 + t520 * t138) * g(2) - t341 * t415 / 0.2e1 + t158 * t22 + t159 * t23 + t105 * t421 / 0.2e1 + t562 * t94 + t563 * t93 + (t1 * t159 + t14 * t225 + t15 * t562 + t158 * t2 + t16 * t563 + t50 * t549) * m(7); (-t230 * t589 + t105 + t226 - t457) * t593 - t579 - t125 * t503 + t358 * (mrSges(5,1) * t231 - mrSges(5,2) * t230) + t527 + (-t230 * t588 - t231 * t586) * t592 + t580 * t495 + (t230 * t51 + t231 * t54) * mrSges(6,2) + (t540 + t470) * t70 + t332 * t505 - t135 * t136 - t511 - t95 * t49 + qJ(5) * t83 - pkin(4) * t81 + (-t434 - t471) * t69 + qJD(5) * t169 + (-Ifges(5,2) * t231 - t227 + t577) * t488 - t73 * (mrSges(6,1) * t231 + mrSges(6,3) * t230) + t243 * t22 + t244 * t23 + (t301 * t79 + t535 * t78) * g(1) + (t146 * t535 + t147 * t301) * g(3) + (t301 * t77 + t535 * t76) * g(2) + (-pkin(4) * t13 + qJ(5) * t12 - t135 * t73 - t51 * t70 + t54 * t548) * m(6) + t552 * t93 + t554 * t94 + (t1 * t244 + t15 * t554 + t16 * t552 + t2 * t243 - t50 * t95) * m(7) + t108 * t486 + (Ifges(6,3) * t231 - t458) * t489; t293 * t22 + t289 * t23 + t453 * t231 + t333 * qJD(6) + (-t169 - t333) * t275 + t81 + (t1 * t289 + t2 * t293 - t231 * t50 + t321 + t267 * (-t15 * t289 + t16 * t293)) * m(7) + (t231 * t73 - t275 * t54 + t13 + t321) * m(6); t36 * t492 - t15 * t93 + t16 * t94 - g(1) * ((-t289 * t79 + t293 * t78) * mrSges(7,1) + (-t289 * t78 - t293 * t79) * mrSges(7,2)) - g(2) * ((-t289 * t77 + t293 * t76) * mrSges(7,1) + (-t289 * t76 - t293 * t77) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t52 - mrSges(7,2) * t53) + (t37 - t580) * t495 + t579;];
tau  = t3;
