% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:35:08
% EndTime: 2019-03-09 16:35:48
% DurationCPUTime: 25.52s
% Computational Cost: add. (17159->768), mult. (39867->975), div. (0->0), fcn. (29418->14), ass. (0->365)
t617 = Ifges(6,1) + Ifges(7,1);
t588 = Ifges(7,4) + Ifges(6,5);
t622 = -mrSges(6,3) - mrSges(7,2);
t625 = -Ifges(6,4) + Ifges(7,5);
t587 = Ifges(7,2) + Ifges(6,3);
t616 = Ifges(6,6) - Ifges(7,6);
t351 = sin(qJ(5));
t444 = qJD(5) * t351;
t352 = sin(qJ(3));
t353 = sin(qJ(2));
t356 = cos(qJ(3));
t357 = cos(qJ(2));
t283 = -t352 * t353 + t356 * t357;
t265 = t283 * qJD(1);
t284 = t352 * t357 + t353 * t356;
t266 = t284 * qJD(1);
t349 = sin(pkin(10));
t350 = cos(pkin(10));
t408 = t350 * t265 - t266 * t349;
t484 = t408 * t351;
t624 = t444 - t484;
t355 = cos(qJ(5));
t443 = qJD(5) * t355;
t483 = t408 * t355;
t605 = -t483 + t443;
t441 = qJD(1) * qJD(2);
t293 = qJDD(1) * t357 - t353 * t441;
t294 = qJDD(1) * t353 + t357 * t441;
t371 = t283 * qJD(3);
t201 = qJD(1) * t371 + t293 * t352 + t294 * t356;
t372 = t284 * qJD(3);
t202 = -qJD(1) * t372 + t293 * t356 - t294 * t352;
t136 = t201 * t350 + t202 * t349;
t346 = qJDD(2) + qJDD(3);
t381 = t265 * t349 + t350 * t266;
t439 = qJD(2) + qJD(3);
t197 = t351 * t381 - t355 * t439;
t445 = qJD(5) * t197;
t90 = t355 * t136 + t351 * t346 - t445;
t544 = t90 / 0.2e1;
t198 = t351 * t439 + t355 * t381;
t91 = qJD(5) * t198 + t351 * t136 - t355 * t346;
t542 = t91 / 0.2e1;
t135 = -t201 * t349 + t202 * t350;
t130 = qJDD(5) - t135;
t541 = t130 / 0.2e1;
t209 = qJD(5) - t408;
t623 = t209 / 0.2e1;
t591 = mrSges(6,1) + mrSges(7,1);
t590 = mrSges(6,2) - mrSges(7,3);
t193 = Ifges(6,4) * t197;
t492 = t197 * Ifges(7,5);
t621 = t198 * t617 + t588 * t209 - t193 + t492;
t507 = mrSges(7,2) * t197;
t147 = mrSges(7,3) * t209 - t507;
t505 = mrSges(6,3) * t197;
t148 = -mrSges(6,2) * t209 - t505;
t572 = t147 + t148;
t504 = mrSges(6,3) * t198;
t149 = mrSges(6,1) * t209 - t504;
t506 = mrSges(7,2) * t198;
t150 = -mrSges(7,1) * t209 + t506;
t571 = t149 - t150;
t359 = -pkin(8) - pkin(7);
t315 = t359 * t357;
t290 = qJD(1) * t315;
t267 = t352 * t290;
t314 = t359 * t353;
t289 = qJD(1) * t314;
t275 = qJD(2) * pkin(2) + t289;
t222 = t356 * t275 + t267;
t258 = t266 * qJ(4);
t190 = t222 - t258;
t179 = pkin(3) * t439 + t190;
t270 = t356 * t290;
t223 = t275 * t352 - t270;
t487 = qJ(4) * t265;
t191 = t223 + t487;
t182 = t349 * t191;
t113 = t350 * t179 - t182;
t345 = t357 * pkin(2);
t334 = t345 + pkin(1);
t313 = t334 * qJD(1);
t235 = -pkin(3) * t265 + qJD(4) - t313;
t620 = t235 * mrSges(5,2) - t113 * mrSges(5,3);
t464 = t350 * t191;
t114 = t349 * t179 + t464;
t110 = pkin(9) * t439 + t114;
t122 = -pkin(4) * t408 - pkin(9) * t381 + t235;
t45 = -t110 * t351 + t122 * t355;
t35 = -pkin(5) * t209 + qJD(6) - t45;
t46 = t110 * t355 + t122 * t351;
t36 = qJ(6) * t209 + t46;
t619 = t235 * mrSges(5,1) + t45 * mrSges(6,1) - t35 * mrSges(7,1) - t46 * mrSges(6,2) - t114 * mrSges(5,3) + t36 * mrSges(7,3);
t598 = m(6) + m(7);
t618 = m(5) + t598;
t615 = t130 * t588 + t617 * t90 + t625 * t91;
t614 = -t197 * t616 + t198 * t588 + t587 * t209;
t229 = -t289 * t352 + t270;
t199 = t229 - t487;
t230 = t356 * t289 + t267;
t200 = -t258 + t230;
t463 = t350 * t352;
t495 = pkin(2) * qJD(3);
t574 = -t350 * t199 + t200 * t349 - (t349 * t356 + t463) * t495;
t613 = pkin(5) * t624 - qJ(6) * t605 - qJD(6) * t351;
t348 = qJ(2) + qJ(3);
t342 = pkin(10) + t348;
t326 = sin(t342);
t612 = t622 * t326;
t389 = pkin(5) * t355 + qJ(6) * t351;
t379 = -pkin(4) - t389;
t399 = -t355 * mrSges(7,1) - t351 * mrSges(7,3);
t488 = t355 * mrSges(6,1);
t563 = (-m(7) * t379 - t399 + t488) * t326;
t611 = -t351 * t616 + t355 * t588;
t497 = Ifges(7,5) * t351;
t499 = Ifges(6,4) * t351;
t610 = t355 * t617 + t497 - t499;
t327 = cos(t342);
t343 = sin(t348);
t344 = cos(t348);
t609 = mrSges(4,1) * t343 + mrSges(5,1) * t326 + mrSges(4,2) * t344 + mrSges(5,2) * t327;
t608 = -t344 * mrSges(4,1) - t327 * mrSges(5,1) + mrSges(4,2) * t343 + t326 * mrSges(5,2);
t280 = t294 * pkin(7);
t234 = qJDD(2) * pkin(2) - pkin(8) * t294 - t280;
t279 = t293 * pkin(7);
t237 = pkin(8) * t293 + t279;
t144 = -qJD(3) * t223 + t356 * t234 - t237 * t352;
t79 = pkin(3) * t346 - qJ(4) * t201 - qJD(4) * t266 + t144;
t446 = qJD(3) * t356;
t447 = qJD(3) * t352;
t143 = t352 * t234 + t356 * t237 + t275 * t446 + t290 * t447;
t89 = qJ(4) * t202 + qJD(4) * t265 + t143;
t29 = t349 * t79 + t350 * t89;
t26 = pkin(9) * t346 + t29;
t486 = qJDD(1) * pkin(1);
t259 = -pkin(2) * t293 - t486;
t167 = -pkin(3) * t202 + qJDD(4) + t259;
t34 = -pkin(4) * t135 - pkin(9) * t136 + t167;
t6 = -t110 * t444 + t122 * t443 + t355 * t26 + t351 * t34;
t7 = -qJD(5) * t46 - t26 * t351 + t34 * t355;
t604 = -t351 * t7 + t355 * t6;
t2 = qJ(6) * t130 + qJD(6) * t209 + t6;
t4 = -pkin(5) * t130 + qJDD(6) - t7;
t603 = t2 * t355 + t351 * t4;
t354 = sin(qJ(1));
t358 = cos(qJ(1));
t602 = g(1) * t358 + g(2) * t354;
t109 = -pkin(4) * t439 - t113;
t398 = t351 * mrSges(7,1) - t355 * mrSges(7,3);
t400 = mrSges(6,1) * t351 + mrSges(6,2) * t355;
t48 = t197 * pkin(5) - t198 * qJ(6) + t109;
t601 = t109 * t400 + t398 * t48;
t600 = Ifges(7,5) * t544 + Ifges(7,6) * t541 - t90 * Ifges(6,4) / 0.2e1 - t130 * Ifges(6,6) / 0.2e1 + (Ifges(7,3) + Ifges(6,2)) * t542;
t39 = -mrSges(6,2) * t130 - mrSges(6,3) * t91;
t40 = -mrSges(7,2) * t91 + mrSges(7,3) * t130;
t584 = t39 + t40;
t37 = mrSges(6,1) * t130 - mrSges(6,3) * t90;
t38 = -t130 * mrSges(7,1) + t90 * mrSges(7,2);
t585 = -t37 + t38;
t599 = m(6) * ((-t351 * t46 - t355 * t45) * qJD(5) + t604) + m(7) * ((t35 * t355 - t351 * t36) * qJD(5) + t603) - t572 * t444 - t571 * t443 + t355 * t584 + t351 * t585;
t597 = t293 / 0.2e1;
t529 = t346 / 0.2e1;
t527 = t357 / 0.2e1;
t533 = -t381 / 0.2e1;
t596 = t381 / 0.2e1;
t595 = -t408 / 0.2e1;
t594 = t408 / 0.2e1;
t593 = -t439 / 0.2e1;
t592 = t439 / 0.2e1;
t583 = Ifges(5,4) * t381;
t582 = Ifges(5,4) * t408;
t581 = Ifges(4,5) * t284;
t580 = Ifges(3,2) * t357;
t579 = Ifges(4,6) * t283;
t119 = t190 * t349 + t464;
t577 = -t119 + t613;
t576 = t613 - t574;
t134 = t199 * t349 + t200 * t350;
t465 = t349 * t352;
t255 = (t350 * t356 - t465) * t495;
t573 = -t134 + t255;
t225 = -t350 * t283 + t284 * t349;
t226 = t283 * t349 + t284 * t350;
t243 = -pkin(3) * t283 - t334;
t153 = pkin(4) * t225 - pkin(9) * t226 + t243;
t238 = t356 * t314 + t315 * t352;
t212 = -qJ(4) * t284 + t238;
t239 = t352 * t314 - t356 * t315;
t213 = qJ(4) * t283 + t239;
t155 = t212 * t349 + t213 * t350;
t570 = t351 * t153 + t355 * t155;
t569 = mrSges(5,1) * t439 - mrSges(6,1) * t197 - mrSges(6,2) * t198 - mrSges(5,3) * t381;
t474 = t326 * t351;
t434 = mrSges(6,2) * t474;
t471 = t327 * t354;
t568 = -t354 * t434 + t471 * t622;
t469 = t327 * t358;
t567 = -t358 * t434 + t469 * t622;
t565 = -t327 * pkin(4) - t326 * pkin(9);
t520 = pkin(4) * t326;
t523 = pkin(3) * t343;
t564 = m(7) * t523 - m(6) * (-t520 - t523) + t563;
t561 = t130 * t587 + t588 * t90 - t616 * t91;
t449 = qJD(1) * t357;
t450 = qJD(1) * t353;
t517 = pkin(7) * t357;
t518 = pkin(7) * t353;
t560 = (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t449) * t518 + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t450) * t517;
t559 = t279 * t357 + t280 * t353;
t557 = 0.2e1 * t529;
t553 = -m(3) * pkin(7) + m(4) * t359 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t428 = qJD(2) * t359;
t291 = t353 * t428;
t292 = t357 * t428;
t175 = t356 * t291 + t352 * t292 + t314 * t446 + t315 * t447;
t232 = -qJD(2) * t284 - t372;
t131 = qJ(4) * t232 + qJD(4) * t283 + t175;
t176 = -qJD(3) * t239 - t291 * t352 + t356 * t292;
t231 = qJD(2) * t283 + t371;
t132 = -qJ(4) * t231 - qJD(4) * t284 + t176;
t53 = t131 * t350 + t132 * t349;
t168 = t231 * t349 - t350 * t232;
t169 = t231 * t350 + t232 * t349;
t448 = qJD(2) * t353;
t339 = pkin(2) * t448;
t218 = -pkin(3) * t232 + t339;
t77 = pkin(4) * t168 - pkin(9) * t169 + t218;
t13 = -qJD(5) * t570 - t351 * t53 + t355 * t77;
t552 = m(7) * pkin(5) + t591;
t311 = -mrSges(3,1) * t357 + mrSges(3,2) * t353;
t551 = -m(3) * pkin(1) - m(4) * t334 - mrSges(2,1) + t311 + t608;
t470 = t327 * t355;
t472 = t327 * t351;
t550 = -t470 * t591 + t472 * t590 + t608 + t612;
t549 = m(7) * qJ(6) - t590;
t304 = pkin(9) * t471;
t307 = pkin(9) * t469;
t548 = -g(1) * t307 - g(2) * t304;
t547 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t543 = -t91 / 0.2e1;
t539 = -t197 / 0.2e1;
t538 = t197 / 0.2e1;
t537 = -t198 / 0.2e1;
t536 = t198 / 0.2e1;
t535 = -t209 / 0.2e1;
t530 = t266 / 0.2e1;
t528 = t351 / 0.2e1;
t526 = pkin(2) * t353;
t525 = pkin(2) * t356;
t524 = pkin(3) * t266;
t330 = pkin(3) * t344;
t522 = pkin(3) * t349;
t521 = pkin(3) * t350;
t519 = pkin(5) * t381;
t514 = g(3) * t326;
t503 = Ifges(3,4) * t353;
t502 = Ifges(3,4) * t357;
t501 = Ifges(4,4) * t266;
t500 = Ifges(6,4) * t198;
t498 = Ifges(6,4) * t355;
t496 = Ifges(7,5) * t355;
t491 = t265 * mrSges(4,3);
t490 = t266 * mrSges(4,3);
t485 = t169 * t355;
t478 = t255 * t351;
t477 = t255 * t355;
t473 = t326 * t358;
t101 = -Ifges(6,2) * t197 + Ifges(6,6) * t209 + t500;
t462 = t351 * t101;
t461 = t351 * t354;
t460 = t354 * t355;
t459 = t355 * t358;
t458 = t358 * t351;
t120 = t190 * t350 - t182;
t145 = pkin(4) * t381 - pkin(9) * t408 + t524;
t55 = t355 * t120 + t351 * t145;
t337 = pkin(2) * t450;
t137 = t145 + t337;
t57 = t355 * t134 + t351 * t137;
t295 = -t523 - t526;
t455 = t354 * t295 + t304;
t454 = t358 * t295 + t307;
t333 = pkin(3) + t525;
t257 = pkin(2) * t463 + t349 * t333;
t451 = t330 + t345;
t192 = Ifges(7,5) * t198;
t98 = Ifges(7,6) * t209 + Ifges(7,3) * t197 + t192;
t431 = t98 * t528;
t429 = t330 - t565;
t427 = t226 * t443;
t414 = -t444 / 0.2e1;
t413 = t443 / 0.2e1;
t28 = -t349 * t89 + t350 * t79;
t412 = t441 / 0.2e1;
t411 = -t135 * mrSges(5,1) + t136 * mrSges(5,2);
t52 = t131 * t349 - t350 * t132;
t154 = -t350 * t212 + t213 * t349;
t288 = pkin(1) + t451;
t347 = -qJ(4) + t359;
t407 = t358 * t288 - t347 * t354;
t256 = -pkin(2) * t465 + t333 * t350;
t404 = mrSges(3,1) * t353 + mrSges(3,2) * t357;
t395 = t503 + t580;
t394 = -Ifges(6,2) * t351 + t498;
t392 = Ifges(3,5) * t357 - Ifges(3,6) * t353;
t390 = Ifges(7,3) * t351 + t496;
t388 = pkin(5) * t351 - qJ(6) * t355;
t386 = t35 * t351 + t355 * t36;
t385 = -t351 * t45 + t355 * t46;
t383 = pkin(5) * t470 + qJ(6) * t472 + t429;
t54 = -t120 * t351 + t145 * t355;
t56 = -t134 * t351 + t137 * t355;
t66 = t153 * t355 - t155 * t351;
t25 = -pkin(4) * t346 - t28;
t376 = pkin(1) * t404;
t375 = t169 * t351 + t427;
t374 = t226 * t444 - t485;
t373 = t353 * (Ifges(3,1) * t357 - t503);
t12 = t153 * t443 - t155 * t444 + t351 * t77 + t355 * t53;
t157 = Ifges(5,2) * t408 + t439 * Ifges(5,6) + t583;
t158 = Ifges(5,1) * t381 + t439 * Ifges(5,5) + t582;
t210 = Ifges(4,2) * t265 + t439 * Ifges(4,6) + t501;
t260 = Ifges(4,4) * t265;
t211 = Ifges(4,1) * t266 + t439 * Ifges(4,5) + t260;
t9 = pkin(5) * t91 - qJ(6) * t90 - qJD(6) * t198 + t25;
t360 = (-t496 + t498) * t544 - t266 * (Ifges(4,1) * t265 - t501) / 0.2e1 + (Ifges(4,3) + Ifges(5,3)) * t346 + (-t583 + t614) * t533 + t223 * t490 + t222 * t491 + (-Ifges(5,2) * t595 - Ifges(5,6) * t593 + Ifges(6,6) * t538 + Ifges(7,6) * t539 + t535 * t587 + t537 * t588 - t619) * t381 + (t35 * t605 - t36 * t624 + t603) * mrSges(7,2) + (-t45 * t605 - t46 * t624 + t604) * mrSges(6,3) + t621 * (t413 - t483 / 0.2e1) - t98 * t484 / 0.2e1 + (-t394 / 0.2e1 + t390 / 0.2e1) * t445 + (t198 * t610 + t209 * t611) * qJD(5) / 0.2e1 - t25 * t488 + t9 * t399 + (Ifges(4,5) * t265 - Ifges(4,6) * t266) * t593 + t462 * t594 + t157 * t596 + t313 * (mrSges(4,1) * t266 + mrSges(4,2) * t265) - (-Ifges(4,2) * t266 + t211 + t260) * t265 / 0.2e1 + (Ifges(5,1) * t533 + Ifges(5,5) * t593 + t390 * t539 + t394 * t538 + t535 * t611 + t537 * t610 - t601 - t620) * t408 + t497 * t542 + t499 * t543 + t101 * t414 + Ifges(4,6) * t202 + Ifges(4,5) * t201 + t144 * mrSges(4,1) - t143 * mrSges(4,2) + Ifges(5,6) * t135 + Ifges(5,5) * t136 + t28 * mrSges(5,1) - t29 * mrSges(5,2) + (t582 + t158) * t595 + t615 * t528 + (Ifges(6,2) * t543 - Ifges(7,3) * t542 + t541 * t616 - t600) * t355 + (t25 * mrSges(6,2) + t541 * t588 + t544 * t617) * t351 + t210 * t530 + (t431 + t601) * qJD(5);
t336 = Ifges(3,4) * t449;
t329 = -pkin(4) - t521;
t274 = t379 - t521;
t264 = Ifges(3,1) * t450 + Ifges(3,5) * qJD(2) + t336;
t263 = Ifges(3,6) * qJD(2) + qJD(1) * t395;
t251 = -pkin(4) - t256;
t250 = t327 * t459 + t461;
t249 = t327 * t458 - t460;
t248 = t327 * t460 - t458;
t247 = t327 * t461 + t459;
t242 = mrSges(4,1) * t439 - t490;
t241 = -mrSges(4,2) * t439 + t491;
t240 = t337 + t524;
t236 = -t256 + t379;
t220 = -mrSges(4,1) * t265 + mrSges(4,2) * t266;
t207 = t381 * qJ(6);
t203 = -mrSges(5,2) * t439 + mrSges(5,3) * t408;
t181 = -mrSges(4,2) * t346 + mrSges(4,3) * t202;
t180 = mrSges(4,1) * t346 - mrSges(4,3) * t201;
t166 = -mrSges(5,1) * t408 + mrSges(5,2) * t381;
t139 = mrSges(7,1) * t197 - mrSges(7,3) * t198;
t138 = pkin(5) * t198 + qJ(6) * t197;
t118 = mrSges(5,1) * t346 - mrSges(5,3) * t136;
t117 = -mrSges(5,2) * t346 + mrSges(5,3) * t135;
t92 = t226 * t388 + t154;
t51 = -pkin(5) * t225 - t66;
t50 = qJ(6) * t225 + t570;
t44 = -t56 - t519;
t43 = t207 + t57;
t42 = -t54 - t519;
t41 = t207 + t55;
t31 = mrSges(6,1) * t91 + mrSges(6,2) * t90;
t30 = mrSges(7,1) * t91 - mrSges(7,3) * t90;
t14 = t388 * t169 + (qJD(5) * t389 - qJD(6) * t355) * t226 + t52;
t11 = -pkin(5) * t168 - t13;
t10 = qJ(6) * t168 + qJD(6) * t225 + t12;
t1 = [(t579 / 0.2e1 + t581 / 0.2e1) * t346 + (t579 + t581) * t529 + (t621 * t414 + t167 * mrSges(5,2) - t28 * mrSges(5,3) + Ifges(5,1) * t136 + Ifges(5,4) * t135 + t557 * Ifges(5,5) + t25 * t400 + t390 * t542 + t394 * t543 + t9 * t398 + t98 * t413 + t541 * t611 + t544 * t610 + (-t2 * mrSges(7,2) - t6 * mrSges(6,3) + t600) * t351 + (mrSges(7,2) * t4 - mrSges(6,3) * t7 + t615 / 0.2e1) * t355) * t226 + (t357 * t502 + t373) * t412 + (t143 * t283 - t144 * t284 - t222 * t231 + t223 * t232) * mrSges(4,3) + (t293 * t517 + t294 * t518 + t559) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t559) + (t264 * t527 + t392 * qJD(2) / 0.2e1 - t560) * qJD(2) + t621 * t485 / 0.2e1 + (-Ifges(7,5) * t374 + Ifges(7,3) * t375) * t538 + (-Ifges(6,4) * t374 - Ifges(6,2) * t375) * t539 + (mrSges(4,1) * t334 + Ifges(4,4) * t284 + Ifges(4,2) * t283) * t202 + m(4) * (t143 * t239 + t144 * t238 + t175 * t223 + t176 * t222 - t259 * t334 - t313 * t339) - t101 * t427 / 0.2e1 + m(5) * (t114 * t53 + t155 * t29 + t167 * t243 + t218 * t235) + (t552 * t248 + t549 * t247 + (m(5) * t288 - t598 * (-t288 + t565) - t551 - t612) * t354 + (t347 * t618 + t553) * t358) * g(1) + (t374 * t45 - t375 * t46) * mrSges(6,3) + m(7) * (t10 * t36 + t11 * t35 + t14 * t48 + t2 * t50 + t4 * t51 + t9 * t92) + (-t374 * t588 - t375 * t616) * t623 - t311 * t486 + (Ifges(3,4) * t294 + Ifges(3,2) * t293) * t527 + t109 * (mrSges(6,1) * t375 - mrSges(6,2) * t374) + t48 * (mrSges(7,1) * t375 + mrSges(7,3) * t374) + (Ifges(4,5) * t231 + Ifges(4,6) * t232) * t592 + t395 * t597 + (-m(5) * t28 + m(6) * t25 - t118 + t31) * t154 - t376 * t441 + m(6) * (t12 * t46 + t13 * t45 + t570 * t6 + t66 * t7) + t570 * t39 + (Ifges(5,5) * t592 + Ifges(5,4) * t594 + Ifges(5,1) * t596 - t462 / 0.2e1 + t158 / 0.2e1 + t431 + t620) * t169 - t263 * t448 / 0.2e1 + t294 * t502 / 0.2e1 + (-t35 * t374 - t36 * t375) * mrSges(7,2) + (-m(5) * t113 + m(6) * t109 - t569) * t52 + (-mrSges(4,2) * t334 + Ifges(4,1) * t284 + Ifges(4,4) * t283) * t201 + (t561 / 0.2e1 + t167 * mrSges(5,1) - t29 * mrSges(5,3) - Ifges(5,4) * t136 - Ifges(5,2) * t135 - Ifges(5,6) * t557 + Ifges(6,6) * t543 + Ifges(7,6) * t542 + t541 * t587 + t544 * t588 + t547) * t225 + t53 * t203 + t155 * t117 + t10 * t147 + t12 * t148 + t13 * t149 + t11 * t150 + t14 * t139 + Ifges(2,3) * qJDD(1) + t92 * t30 + t66 * t37 + t50 * t40 + t51 * t38 + t243 * t411 + (-mrSges(3,1) * t518 - mrSges(3,2) * t517 + 0.2e1 * Ifges(3,6) * t527) * qJDD(2) + (Ifges(3,1) * t294 + Ifges(3,4) * t597 + Ifges(3,5) * qJDD(2) - t412 * t580) * t353 + t218 * t166 + t231 * t211 / 0.2e1 + t232 * t210 / 0.2e1 + (-t617 * t374 + t375 * t625) * t536 + t238 * t180 + t239 * t181 + t175 * t241 + t220 * t339 + t176 * t242 + (-m(5) * t407 + t622 * t473 - t598 * (pkin(4) * t469 + pkin(9) * t473 + t407) - t552 * t250 - t549 * t249 + t551 * t358 + t553 * t354) * g(2) + t265 * (Ifges(4,4) * t231 + Ifges(4,2) * t232) / 0.2e1 + (t614 / 0.2e1 - Ifges(5,6) * t592 - Ifges(5,2) * t594 - Ifges(5,4) * t596 + t587 * t623 + t588 * t536 - t157 / 0.2e1 + Ifges(7,6) * t538 + Ifges(6,6) * t539 + t619) * t168 + t259 * (-mrSges(4,1) * t283 + mrSges(4,2) * t284) + (Ifges(4,1) * t231 + Ifges(4,4) * t232) * t530 - pkin(1) * (-mrSges(3,1) * t293 + mrSges(3,2) * t294) - t313 * (-mrSges(4,1) * t232 + mrSges(4,2) * t231); t360 - t220 * t337 + (t477 - t43) * t147 + t574 * t569 + (t478 - t44) * t150 + t573 * t203 + (t113 * t574 + t114 * t573 - t235 * t240 + t256 * t28 + t257 * t29) * m(5) + (-g(1) * (-pkin(4) * t473 + t454) - t45 * t56 - t46 * t57 + t25 * t251 + t255 * t385 - g(2) * (-t354 * t520 + t455) - t574 * t109) * m(6) + (m(4) * (t143 * t352 + t144 * t356 + (-t222 * t352 + t223 * t356) * qJD(3)) - t242 * t447 + t241 * t446 + t352 * t181) * pkin(2) + (-g(1) * t454 - g(2) * t455 + t236 * t9 + t255 * t386 - t35 * t44 - t36 * t43 + t48 * t576) * m(7) + t576 * t139 + (-m(7) * (t345 + t383) - m(6) * (t345 + t429) - m(5) * t451 - m(4) * t345 + t311 + t550) * g(3) + t602 * (m(4) * t526 - m(5) * t295 + t404 + t609) + (-t478 - t56) * t149 - (-Ifges(3,2) * t450 + t264 + t336) * t449 / 0.2e1 + (t560 + (-t373 / 0.2e1 + t376) * qJD(1)) * qJD(1) - t392 * t441 / 0.2e1 - m(4) * (t222 * t229 + t223 * t230 - t313 * t337) + t599 * (pkin(9) + t257) + t263 * t450 / 0.2e1 + (t563 * t358 + t567) * g(1) + (t563 * t354 + t568) * g(2) + (t477 - t57) * t148 + Ifges(3,3) * qJDD(2) + t236 * t30 - t240 * t166 - t230 * t241 - t229 * t242 + t251 * t31 + t256 * t118 + t257 * t117 - t279 * mrSges(3,2) - t280 * mrSges(3,1) + Ifges(3,6) * t293 + Ifges(3,5) * t294 + t180 * t525; t360 + ((t28 * t350 + t29 * t349) * pkin(3) + t113 * t119 - t114 * t120 - t235 * t524) * m(5) + (-t109 * t119 + t25 * t329 - t45 * t54 - t46 * t55 + t548) * m(6) - t166 * t524 + (t274 * t9 - t35 * t42 - t36 * t41 + t48 * t577 + t548) * m(7) + t577 * t139 + (-m(5) * t330 - m(6) * t429 - m(7) * t383 + t550) * g(3) + (t358 * t564 + t567) * g(1) + (t354 * t564 + t568) * g(2) + t569 * t119 - t120 * t203 - t41 * t147 - t55 * t148 - t54 * t149 - t42 * t150 - t222 * t241 + t223 * t242 + t274 * t30 + t118 * t521 + t117 * t522 + t329 * t31 + (m(5) * t523 + t609) * t602 + t599 * (pkin(9) + t522); -t408 * t203 + (-t139 + t569) * t381 + (t209 * t572 - t585) * t355 + (-t209 * t571 + t584) * t351 + t411 + (-g(1) * t354 + g(2) * t358) * t618 + (t2 * t351 + t209 * t386 - t355 * t4 - t381 * t48) * m(7) + (-t109 * t381 + t209 * t385 + t351 * t6 + t355 * t7) * m(6) + (t113 * t381 - t114 * t408 + t167) * m(5); (t398 + t400) * t514 + t547 + (-t197 * t588 - t198 * t616) * t535 + (t247 * t591 + t248 * t590) * g(2) + (-m(7) * t36 - t505 - t572) * t45 + t36 * t506 + t35 * t507 + (t249 * t591 + t250 * t590) * g(1) + t561 + (-m(7) * t35 + t504 + t571) * t46 + (-t197 * t617 + t192 - t500 + t98) * t537 - t48 * (mrSges(7,1) * t198 + mrSges(7,3) * t197) - t109 * (mrSges(6,1) * t198 - mrSges(6,2) * t197) + (-Ifges(6,2) * t198 - t193 + t621) * t538 + qJD(6) * t147 - t138 * t139 - pkin(5) * t38 + qJ(6) * t40 + t101 * t536 + (Ifges(7,3) * t198 - t492) * t539 + (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t36 + t388 * t514 - t138 * t48 - g(2) * (-pkin(5) * t247 + qJ(6) * t248) - g(1) * (-pkin(5) * t249 + qJ(6) * t250)) * m(7); t198 * t139 - t209 * t147 + (-g(1) * t249 - g(2) * t247 - g(3) * t474 + t48 * t198 - t36 * t209 + t4) * m(7) + t38;];
tau  = t1;
