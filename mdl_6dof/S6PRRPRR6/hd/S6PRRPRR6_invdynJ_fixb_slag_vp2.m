% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:55
% EndTime: 2019-03-08 22:23:54
% DurationCPUTime: 34.07s
% Computational Cost: add. (16927->899), mult. (42576->1272), div. (0->0), fcn. (36416->18), ass. (0->417)
t341 = sin(pkin(13));
t344 = cos(pkin(13));
t349 = sin(qJ(5));
t353 = cos(qJ(5));
t303 = t341 * t353 + t344 * t349;
t342 = sin(pkin(7));
t354 = cos(qJ(3));
t484 = t342 * t354;
t371 = t303 * t484;
t233 = qJD(2) * t371;
t296 = t303 * qJD(5);
t610 = t233 - t296;
t302 = t341 * t349 - t353 * t344;
t370 = t302 * t484;
t234 = qJD(2) * t370;
t295 = t302 * qJD(5);
t609 = -t234 + t295;
t350 = sin(qJ(3));
t401 = pkin(3) * t350 - qJ(4) * t354;
t469 = qJD(4) * t350;
t241 = (qJD(3) * t401 - t469) * t342;
t345 = cos(pkin(7));
t470 = qJD(3) * t354;
t446 = t345 * t470;
t471 = qJD(3) * t350;
t448 = t342 * t471;
t279 = pkin(2) * t446 - pkin(9) * t448;
t252 = qJD(4) * t345 + t279;
t161 = t344 * t241 - t252 * t341;
t482 = t344 * t354;
t378 = (pkin(4) * t350 - pkin(10) * t482) * t342;
t121 = qJD(3) * t378 + t161;
t481 = t345 * t350;
t298 = pkin(2) * t481 + pkin(9) * t484;
t265 = qJ(4) * t345 + t298;
t402 = -pkin(3) * t354 - qJ(4) * t350;
t266 = (-pkin(2) + t402) * t342;
t184 = -t265 * t341 + t344 * t266;
t485 = t342 * t350;
t292 = t341 * t345 + t344 * t485;
t129 = -pkin(4) * t484 - pkin(10) * t292 + t184;
t162 = t341 * t241 + t344 * t252;
t447 = t342 * t470;
t424 = t341 * t447;
t142 = -pkin(10) * t424 + t162;
t185 = t344 * t265 + t341 * t266;
t290 = -t341 * t485 + t344 * t345;
t143 = pkin(10) * t290 + t185;
t343 = sin(pkin(6));
t521 = cos(qJ(2));
t454 = t521 * t354;
t351 = sin(qJ(2));
t479 = t350 * t351;
t377 = -t345 * t479 + t454;
t262 = t377 * t343;
t246 = qJD(1) * t262;
t483 = t343 * t351;
t452 = qJD(1) * t483;
t427 = t342 * t452;
t199 = -t246 * t341 + t344 * t427;
t200 = t246 * t344 + t341 * t427;
t467 = qJD(5) * t353;
t468 = qJD(5) * t349;
t589 = t129 * t467 - t143 * t468 + (t142 - t200) * t353 + (t121 - t199) * t349;
t472 = qJD(2) * t342;
t299 = pkin(9) * t472 + t452;
t456 = t343 * t521;
t423 = qJD(1) * t456;
t312 = qJD(2) * pkin(2) + t423;
t346 = cos(pkin(6));
t473 = qJD(1) * t346;
t453 = t342 * t473;
t186 = -t350 * t299 + t354 * (t312 * t345 + t453);
t278 = t401 * t472;
t130 = -t186 * t341 + t344 * t278;
t107 = qJD(2) * t378 + t130;
t131 = t344 * t186 + t341 * t278;
t450 = t354 * t472;
t426 = t341 * t450;
t115 = -pkin(10) * t426 + t131;
t514 = pkin(10) + qJ(4);
t315 = t514 * t341;
t316 = t514 * t344;
t390 = -t353 * t315 - t316 * t349;
t585 = -qJD(4) * t302 + qJD(5) * t390 - t349 * t107 - t353 * t115;
t348 = sin(qJ(6));
t352 = cos(qJ(6));
t407 = -mrSges(7,1) * t352 + mrSges(7,2) * t348;
t608 = m(7) * pkin(5) + mrSges(6,1) - t407;
t569 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t607 = -pkin(11) * t448 - t589;
t451 = t350 * t472;
t606 = -pkin(11) * t451 + t585;
t287 = t312 * t481;
t187 = t354 * t299 + t350 * t453 + t287;
t156 = pkin(4) * t426 + t187;
t605 = -pkin(5) * t610 + t609 * pkin(11) - t156;
t391 = t353 * t290 - t292 * t349;
t137 = -qJD(3) * t370 + qJD(5) * t391;
t207 = t290 * t349 + t292 * t353;
t138 = qJD(3) * t371 + qJD(5) * t207;
t280 = t298 * qJD(3);
t227 = pkin(4) * t424 + t280;
t455 = t521 * t350;
t478 = t351 * t354;
t376 = t345 * t478 + t455;
t261 = t376 * t343;
t245 = qJD(1) * t261;
t604 = pkin(5) * t138 - pkin(11) * t137 + t227 - t245;
t331 = qJD(2) * t345 + qJD(3);
t257 = t331 * t344 - t341 * t451;
t258 = t331 * t341 + t344 * t451;
t433 = t353 * t257 - t258 * t349;
t174 = Ifges(6,4) * t433;
t318 = qJD(5) - t450;
t392 = t257 * t349 + t353 * t258;
t92 = Ifges(6,1) * t392 + t318 * Ifges(6,5) + t174;
t541 = t92 / 0.2e1;
t232 = -t315 * t349 + t316 * t353;
t584 = -qJD(4) * t303 - qJD(5) * t232 - t107 * t353 + t115 * t349;
t322 = qJD(2) * t423;
t286 = qJDD(1) * t483 + t322;
t462 = qJDD(2) * t342;
t603 = pkin(9) * t462 + qJD(3) * t453 + t286;
t163 = qJ(4) * t331 + t187;
t326 = t345 * t473;
t205 = t326 + (qJD(2) * t402 - t312) * t342;
t100 = -t163 * t341 + t344 * t205;
t101 = t344 * t163 + t341 * t205;
t160 = -pkin(3) * t331 + qJD(4) - t186;
t248 = -t312 * t342 + t326;
t508 = Ifges(5,4) * t344;
t509 = Ifges(5,4) * t341;
t602 = t160 * (mrSges(5,1) * t341 + mrSges(5,2) * t344) * t484 + (t101 * (-mrSges(5,3) * t341 * t354 - mrSges(5,2) * t350) + t100 * (mrSges(5,1) * t350 - mrSges(5,3) * t482) + t248 * (mrSges(4,1) * t350 + mrSges(4,2) * t354)) * t342 + (t331 * (Ifges(4,5) * t354 - Ifges(4,6) * t350) + t257 * (Ifges(5,6) * t350 + (-Ifges(5,2) * t341 + t508) * t354) + t258 * (Ifges(5,5) * t350 + (Ifges(5,1) * t344 - t509) * t354)) * t342 / 0.2e1;
t140 = t318 * t352 - t348 * t392;
t141 = t318 * t348 + t352 * t392;
t175 = qJD(6) - t433;
t49 = t141 * Ifges(7,5) + t140 * Ifges(7,6) + t175 * Ifges(7,3);
t507 = Ifges(6,4) * t392;
t91 = Ifges(6,2) * t433 + t318 * Ifges(6,6) + t507;
t601 = -t49 / 0.2e1 + t91 / 0.2e1;
t283 = qJD(2) * t448 - t354 * t462;
t525 = t283 / 0.2e1;
t464 = qJD(2) * qJD(3);
t284 = (qJDD(2) * t350 + t354 * t464) * t342;
t330 = qJDD(2) * t345 + qJDD(3);
t220 = t284 * t344 + t330 * t341;
t527 = t220 / 0.2e1;
t600 = Ifges(5,1) * t527 + Ifges(5,5) * t525;
t219 = -t284 * t341 + t330 * t344;
t528 = t219 / 0.2e1;
t599 = Ifges(5,2) * t528 + Ifges(5,6) * t525;
t271 = qJDD(5) + t283;
t526 = t271 / 0.2e1;
t83 = -qJD(5) * t392 + t219 * t353 - t220 * t349;
t543 = t83 / 0.2e1;
t82 = qJD(5) * t433 + t219 * t349 + t220 * t353;
t544 = t82 / 0.2e1;
t552 = Ifges(6,1) * t544 + Ifges(6,4) * t543 + Ifges(6,5) * t526;
t45 = qJD(6) * t140 + t271 * t348 + t352 * t82;
t551 = t45 / 0.2e1;
t46 = -qJD(6) * t141 + t271 * t352 - t348 * t82;
t550 = t46 / 0.2e1;
t80 = qJDD(6) - t83;
t545 = t80 / 0.2e1;
t598 = -m(6) - m(7);
t582 = t349 * t129 + t353 * t143;
t64 = -pkin(11) * t484 + t582;
t332 = pkin(9) * t485;
t268 = t332 + (-pkin(2) * t354 - pkin(3)) * t345;
t209 = -pkin(4) * t290 + t268;
t89 = -pkin(5) * t391 - pkin(11) * t207 + t209;
t30 = -t348 * t64 + t352 * t89;
t597 = qJD(6) * t30 + t604 * t348 - t607 * t352;
t31 = t348 * t89 + t352 * t64;
t596 = -qJD(6) * t31 + t607 * t348 + t604 * t352;
t449 = qJD(2) * t483;
t422 = qJD(1) * t449;
t285 = qJDD(1) * t456 - t422;
t264 = qJDD(2) * pkin(2) + t285;
t463 = qJDD(1) * t346;
t441 = t342 * t463;
t95 = t264 * t481 - t299 * t471 + t312 * t446 + t350 * t441 + t354 * t603;
t595 = t95 * mrSges(4,2);
t96 = t354 * (t264 * t345 + t441) - qJD(3) * t287 - t299 * t470 - t603 * t350;
t594 = t96 * mrSges(4,1);
t16 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t65 = mrSges(6,1) * t271 - mrSges(6,3) * t82;
t593 = t16 - t65;
t337 = pkin(4) * t344 + pkin(3);
t212 = pkin(5) * t302 - pkin(11) * t303 - t337;
t133 = t212 * t352 - t232 * t348;
t592 = qJD(6) * t133 + t605 * t348 + t606 * t352;
t134 = t212 * t348 + t232 * t352;
t591 = -qJD(6) * t134 - t606 * t348 + t605 * t352;
t117 = -pkin(4) * t257 + t160;
t590 = mrSges(6,2) * t117;
t135 = -t219 * mrSges(5,1) + t220 * mrSges(5,2);
t41 = -t83 * mrSges(6,1) + t82 * mrSges(6,2);
t588 = t135 + t41;
t512 = mrSges(6,3) * t392;
t145 = mrSges(6,1) * t318 - t512;
t70 = -mrSges(7,1) * t140 + mrSges(7,2) * t141;
t587 = t145 - t70;
t586 = t258 * Ifges(5,5) + Ifges(6,5) * t392 + t257 * Ifges(5,6) + Ifges(6,6) * t433 - Ifges(5,3) * t450 + t318 * Ifges(6,3);
t583 = pkin(5) * t451 - t584;
t432 = mrSges(4,3) * t451;
t581 = -mrSges(4,1) * t331 - mrSges(5,1) * t257 + mrSges(5,2) * t258 + t432;
t580 = -t199 + t161;
t579 = -t200 + t162;
t203 = t234 * t348 + t352 * t451;
t465 = qJD(6) * t352;
t386 = -t295 * t348 + t303 * t465;
t578 = t203 + t386;
t204 = -t234 * t352 + t348 * t451;
t466 = qJD(6) * t348;
t385 = t295 * t352 + t303 * t466;
t577 = t204 + t385;
t575 = t345 * t454 - t479;
t221 = -t343 * t575 - t346 * t484;
t215 = t221 * t352;
t375 = t345 * t455 + t478;
t222 = t343 * t375 + t346 * t485;
t291 = -t342 * t456 + t346 * t345;
t164 = -t222 * t341 + t291 * t344;
t165 = t222 * t344 + t291 * t341;
t94 = t164 * t349 + t165 * t353;
t61 = -t348 * t94 + t215;
t576 = -t246 + t279;
t325 = t345 * t463;
t116 = pkin(3) * t283 - qJ(4) * t284 + t325 + (-qJD(2) * t469 - t264) * t342;
t75 = qJ(4) * t330 + qJD(4) * t331 + t95;
t47 = t344 * t116 - t341 * t75;
t48 = t341 * t116 + t344 * t75;
t574 = -t341 * t47 + t344 * t48;
t74 = -pkin(4) * t450 - pkin(10) * t258 + t100;
t81 = pkin(10) * t257 + t101;
t39 = t349 * t74 + t353 * t81;
t29 = pkin(11) * t318 + t39;
t54 = -pkin(5) * t433 - pkin(11) * t392 + t117;
t14 = -t29 * t348 + t352 * t54;
t84 = -pkin(3) * t330 + qJDD(4) - t96;
t60 = -pkin(4) * t219 + t84;
t17 = -pkin(5) * t83 - pkin(11) * t82 + t60;
t37 = pkin(4) * t283 - pkin(10) * t220 + t47;
t42 = pkin(10) * t219 + t48;
t5 = t349 * t37 + t353 * t42 + t74 * t467 - t468 * t81;
t3 = pkin(11) * t271 + t5;
t1 = qJD(6) * t14 + t17 * t348 + t3 * t352;
t15 = t29 * t352 + t348 * t54;
t2 = -qJD(6) * t15 + t17 * t352 - t3 * t348;
t573 = t1 * t352 - t2 * t348;
t406 = mrSges(7,1) * t348 + mrSges(7,2) * t352;
t572 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) - t406;
t511 = Ifges(4,4) * t350;
t556 = t342 ^ 2;
t570 = (t350 * (Ifges(4,1) * t354 - t511) / 0.2e1 - t354 * (Ifges(5,3) * t350 + (Ifges(5,5) * t344 - Ifges(5,6) * t341) * t354) / 0.2e1) * t556;
t6 = -qJD(5) * t39 - t349 * t42 + t353 * t37;
t27 = -qJD(5) * t582 + t121 * t353 - t142 * t349;
t340 = pkin(13) + qJ(5);
t338 = sin(t340);
t339 = cos(t340);
t409 = -mrSges(5,1) * t344 + mrSges(5,2) * t341;
t568 = m(5) * pkin(3) - t569 * t338 + t608 * t339 + mrSges(4,1) - t409;
t38 = -t349 * t81 + t353 * t74;
t28 = -pkin(5) * t318 - t38;
t566 = -m(7) * t28 + t587;
t565 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t564 = -m(5) * t160 - t581;
t563 = -t117 * mrSges(6,1) - t14 * mrSges(7,1) + t15 * mrSges(7,2);
t98 = -mrSges(6,1) * t433 + mrSges(6,2) * t392;
t561 = m(4) * t186 - m(6) * t117 + t564 - t98;
t560 = -mrSges(6,3) + t572;
t524 = -t318 / 0.2e1;
t532 = -t433 / 0.2e1;
t534 = -t175 / 0.2e1;
t536 = -t141 / 0.2e1;
t538 = -t140 / 0.2e1;
t559 = Ifges(7,5) * t536 - Ifges(6,2) * t532 - Ifges(6,6) * t524 + Ifges(7,6) * t538 + Ifges(7,3) * t534 + t563;
t9 = Ifges(7,5) * t45 + Ifges(7,6) * t46 + Ifges(7,3) * t80;
t558 = mrSges(6,1) * t60 + Ifges(7,5) * t551 + Ifges(7,6) * t550 + Ifges(7,3) * t545 + t9 / 0.2e1 + t565 + (-t526 - t271 / 0.2e1) * Ifges(6,6) + (-t543 - t83 / 0.2e1) * Ifges(6,2) + (-t544 - t82 / 0.2e1) * Ifges(6,4);
t523 = t318 / 0.2e1;
t529 = t392 / 0.2e1;
t531 = t433 / 0.2e1;
t533 = t175 / 0.2e1;
t535 = t141 / 0.2e1;
t537 = t140 / 0.2e1;
t557 = -Ifges(6,4) * t529 + Ifges(7,5) * t535 - Ifges(6,2) * t531 - Ifges(6,6) * t523 + Ifges(7,6) * t537 + Ifges(7,3) * t533 - t563 - t601;
t355 = qJD(2) ^ 2;
t554 = Ifges(7,1) * t551 + Ifges(7,4) * t550 + Ifges(7,5) * t545;
t506 = Ifges(7,4) * t141;
t50 = t140 * Ifges(7,2) + t175 * Ifges(7,6) + t506;
t547 = t50 / 0.2e1;
t136 = Ifges(7,4) * t140;
t51 = t141 * Ifges(7,1) + t175 * Ifges(7,5) + t136;
t546 = -t51 / 0.2e1;
t540 = Ifges(5,4) * t527 + t599;
t539 = Ifges(5,4) * t528 + t600;
t530 = -t392 / 0.2e1;
t522 = t352 / 0.2e1;
t520 = pkin(2) * t342;
t513 = mrSges(6,3) * t433;
t510 = Ifges(4,4) * t354;
t505 = Ifges(7,4) * t348;
t504 = Ifges(7,4) * t352;
t500 = cos(pkin(12));
t499 = sin(pkin(12));
t498 = t433 * t348;
t497 = t433 * t352;
t496 = t221 * t348;
t410 = t499 * t521;
t435 = t500 * t351;
t293 = t346 * t435 + t410;
t494 = t293 * t342;
t411 = t500 * t521;
t434 = t499 * t351;
t294 = -t346 * t434 + t411;
t493 = t294 * t342;
t491 = t303 * t348;
t490 = t303 * t352;
t489 = t338 * t342;
t488 = t339 * t342;
t487 = t341 * t342;
t486 = t342 * t344;
t480 = t345 * t354;
t459 = t342 * t483;
t474 = pkin(2) * t456 + pkin(9) * t459;
t461 = Ifges(6,5) * t82 + Ifges(6,6) * t83 + Ifges(6,3) * t271;
t458 = t51 * t522;
t457 = Ifges(4,5) * t284 - Ifges(4,6) * t283 + Ifges(4,3) * t330;
t444 = t485 / 0.2e1;
t440 = -t472 / 0.2e1;
t438 = -t466 / 0.2e1;
t437 = t343 * t500;
t436 = t343 * t499;
t431 = mrSges(4,3) * t450;
t430 = t341 * t459;
t425 = t342 * t449;
t419 = t354 * t440;
t363 = -t346 * t411 + t434;
t281 = t363 * pkin(2);
t416 = pkin(9) * t494 - t281;
t364 = t346 * t410 + t435;
t282 = t364 * pkin(2);
t415 = pkin(9) * t493 - t282;
t414 = t342 * t437;
t413 = t342 * t436;
t412 = (pkin(4) * t341 + pkin(9)) * t342;
t405 = Ifges(7,1) * t352 - t505;
t404 = -Ifges(7,2) * t348 + t504;
t403 = Ifges(7,5) * t352 - Ifges(7,6) * t348;
t85 = -mrSges(7,2) * t175 + mrSges(7,3) * t140;
t86 = mrSges(7,1) * t175 - mrSges(7,3) * t141;
t399 = -t348 * t86 + t352 * t85;
t62 = t352 * t94 + t496;
t67 = t129 * t353 - t143 * t349;
t394 = t353 * t164 - t165 * t349;
t172 = -t207 * t348 - t352 * t484;
t389 = -t207 * t352 + t348 * t484;
t388 = t28 * t406;
t384 = t342 * (-mrSges(4,1) * t354 + mrSges(4,2) * t350);
t383 = (Ifges(4,2) * t354 + t511) * t342;
t360 = t363 * t354;
t166 = t293 * t350 + t345 * t360 + t354 * t414;
t362 = t364 * t354;
t168 = t294 * t350 + t345 * t362 - t354 * t413;
t379 = -g(1) * t168 - g(2) * t166 - g(3) * t221;
t361 = t363 * t350;
t323 = Ifges(4,4) * t450;
t297 = pkin(2) * t480 - t332;
t277 = qJD(2) * t384;
t276 = -mrSges(4,2) * t331 + t431;
t236 = Ifges(4,1) * t451 + t331 * Ifges(4,5) + t323;
t235 = t331 * Ifges(4,6) + qJD(2) * t383;
t229 = mrSges(4,1) * t330 - mrSges(4,3) * t284;
t228 = -mrSges(4,2) * t330 - mrSges(4,3) * t283;
t224 = t342 * t364 + t345 * t436;
t223 = t342 * t363 - t345 * t437;
t218 = -mrSges(5,1) * t450 - mrSges(5,3) * t258;
t217 = mrSges(5,2) * t450 + mrSges(5,3) * t257;
t216 = -t264 * t342 + t325;
t208 = mrSges(4,1) * t283 + mrSges(4,2) * t284;
t198 = -t294 * t481 - t362;
t197 = t294 * t480 - t350 * t364;
t196 = -t293 * t481 - t360;
t195 = t293 * t480 - t361;
t171 = mrSges(5,1) * t283 - mrSges(5,3) * t220;
t170 = -mrSges(5,2) * t283 + mrSges(5,3) * t219;
t169 = t294 * t354 + (-t345 * t364 + t413) * t350;
t167 = t293 * t354 - t345 * t361 - t350 * t414;
t159 = t258 * Ifges(5,1) + t257 * Ifges(5,4) - Ifges(5,5) * t450;
t158 = t258 * Ifges(5,4) + t257 * Ifges(5,2) - Ifges(5,6) * t450;
t153 = t346 * t447 + (t377 * qJD(2) + qJD(3) * t575) * t343;
t152 = t346 * t448 + (qJD(2) * t376 + qJD(3) * t375) * t343;
t150 = t222 * t339 + t291 * t338;
t144 = -mrSges(6,2) * t318 + t513;
t123 = t153 * t344 + t341 * t425;
t122 = -t153 * t341 + t344 * t425;
t108 = -t353 * t199 + t200 * t349;
t105 = t169 * t339 + t224 * t338;
t103 = t167 * t339 + t223 * t338;
t99 = pkin(5) * t392 - pkin(11) * t433;
t72 = qJD(6) * t389 - t137 * t348 + t352 * t448;
t71 = qJD(6) * t172 + t137 * t352 + t348 * t448;
t66 = -mrSges(6,2) * t271 + mrSges(6,3) * t83;
t63 = pkin(5) * t484 - t67;
t33 = qJD(5) * t394 + t122 * t349 + t123 * t353;
t25 = -pkin(5) * t448 - t27;
t21 = -mrSges(7,2) * t80 + mrSges(7,3) * t46;
t20 = mrSges(7,1) * t80 - mrSges(7,3) * t45;
t19 = t348 * t99 + t352 * t38;
t18 = -t348 * t38 + t352 * t99;
t13 = -qJD(6) * t62 + t152 * t352 - t33 * t348;
t12 = qJD(6) * t61 + t152 * t348 + t33 * t352;
t10 = t45 * Ifges(7,4) + t46 * Ifges(7,2) + t80 * Ifges(7,6);
t4 = -pkin(5) * t271 - t6;
t7 = [m(3) * (t346 ^ 2 * qJDD(1) + (t285 * t521 + t286 * t351) * t343) + (-m(6) * t38 - t566) * (qJD(5) * t94 - t353 * t122 + t123 * t349) + m(4) * (t153 * t187 + t216 * t291 + t222 * t95 + t248 * t425) + m(5) * (t100 * t122 + t101 * t123 + t164 * t47 + t165 * t48) + m(6) * (t33 * t39 + t5 * t94) + m(7) * (t1 * t62 + t12 * t15 + t13 * t14 + t2 * t61) + t277 * t425 + t153 * t276 - t561 * t152 + t291 * t208 + t222 * t228 + t123 * t217 + t122 * t218 + t165 * t170 + t164 * t171 + t33 * t144 + t94 * t66 + t12 * t85 + t13 * t86 + t61 * t20 + t62 * t21 + (-qJDD(2) * t483 - t355 * t456) * mrSges(3,2) + (qJDD(2) * t456 - t355 * t483) * mrSges(3,1) + (-m(4) * t96 + m(5) * t84 + m(6) * t60 - t229 + t588) * t221 - (-m(6) * t6 + m(7) * t4 + t593) * t394 + (-m(2) - m(3) - m(4) - m(5) + t598) * g(3) + m(2) * qJDD(1); (-m(4) * t416 - t196 * mrSges(4,1) - mrSges(4,3) * t494 - m(5) * (pkin(3) * t196 + t416) - (t196 * t344 + t293 * t487) * mrSges(5,1) - (-t196 * t341 + t293 * t486) * mrSges(5,2) + mrSges(3,1) * t363 + t293 * mrSges(3,2) + t598 * (t195 * t514 + t196 * t337 + t293 * t412 - t281) + t569 * (t196 * t338 - t293 * t488) - t608 * (t196 * t339 + t293 * t489) + t560 * t195) * g(2) + (-m(4) * t415 - t198 * mrSges(4,1) - mrSges(4,3) * t493 - m(5) * (pkin(3) * t198 + t415) - (t198 * t344 + t294 * t487) * mrSges(5,1) - (-t198 * t341 + t294 * t486) * mrSges(5,2) + mrSges(3,1) * t364 + t294 * mrSges(3,2) + t598 * (t197 * t514 + t198 * t337 + t294 * t412 - t282) + t569 * (t198 * t338 - t294 * t488) - t608 * (t198 * t339 + t294 * t489) + t560 * t197) * g(1) + (-m(5) * (pkin(3) * t262 + t474) - (t262 * t344 + t430) * mrSges(5,1) - (-t262 * t341 + t344 * t459) * mrSges(5,2) - m(4) * t474 - t262 * mrSges(4,1) - mrSges(4,3) * t459 - (mrSges(3,1) * t521 - mrSges(3,2) * t351) * t343 + t598 * (pkin(4) * t430 + t261 * t514 + t262 * t337 + t474) + t569 * (t262 * t338 - t339 * t459) - t608 * (t262 * t339 + t338 * t459) + t560 * t261) * g(3) - t208 * t520 + t576 * t276 + (-t186 * t280 + t187 * t576 - t216 * t520 - t248 * t427 + t297 * t96 + t298 * t95) * m(4) + t579 * t217 + (t100 * t580 + t101 * t579 + t160 * t280 + t184 * t47 + t185 * t48 + t268 * t84) * m(5) + t580 * t218 + (Ifges(6,5) * t137 + Ifges(6,3) * t448) * t523 + (Ifges(6,5) * t207 - Ifges(6,3) * t484) * t526 - t158 * t424 / 0.2e1 + t207 * t552 + t292 * t539 + t290 * t540 + t137 * t541 + t72 * t547 + (Ifges(6,1) * t207 - Ifges(6,5) * t484) * t544 + (Ifges(6,1) * t137 + Ifges(6,5) * t448) * t529 - (Ifges(5,5) * t220 + Ifges(5,6) * t219 + Ifges(5,3) * t283 + t461) * t484 / 0.2e1 + t6 * (-mrSges(6,1) * t484 - mrSges(6,3) * t207) + t47 * (-mrSges(5,1) * t484 - mrSges(5,3) * t292) + t48 * (mrSges(5,2) * t484 + mrSges(5,3) * t290) + (Ifges(4,4) * t284 - Ifges(4,2) * t283 + Ifges(4,6) * t330) * t484 / 0.2e1 + t216 * t384 + (Ifges(7,1) * t71 + Ifges(7,4) * t72) * t535 + t345 * t594 + (t117 * t137 + t207 * t60 - t39 * t448 + t484 * t5) * mrSges(6,2) + t284 * (Ifges(4,5) * t345 + (Ifges(4,1) * t350 + t510) * t342) / 0.2e1 - t277 * t427 + (Ifges(4,1) * t284 - Ifges(4,4) * t283 + Ifges(4,5) * t330) * t444 + (Ifges(7,5) * t71 + Ifges(7,6) * t72) * t533 + t330 * (Ifges(4,3) * t345 + (Ifges(4,5) * t350 + Ifges(4,6) * t354) * t342) / 0.2e1 + (-mrSges(6,3) * t39 + t557) * t138 - t283 * (Ifges(4,6) * t345 + t383) / 0.2e1 + (t422 + t285) * mrSges(3,1) + t38 * (mrSges(6,1) * t448 - mrSges(6,3) * t137) + (-t186 * t447 - t187 * t448 + t484 * t95 - t485 * t96) * mrSges(4,3) + (t322 - t286) * mrSges(3,2) + (Ifges(5,5) * t292 + Ifges(5,6) * t290 - Ifges(5,3) * t484) * t525 + (Ifges(5,1) * t292 + Ifges(5,4) * t290 - Ifges(5,5) * t484) * t527 + (Ifges(5,4) * t292 + Ifges(5,2) * t290 - Ifges(5,6) * t484) * t528 + (Ifges(6,4) * t207 - Ifges(6,6) * t484) * t543 + (Ifges(6,4) * t137 + Ifges(6,6) * t448) * t531 + (Ifges(7,4) * t71 + Ifges(7,2) * t72) * t537 + t345 * t457 / 0.2e1 + ((t159 * t344 + t236) * t342 + t556 * qJD(2) * (-Ifges(4,2) * t350 + t510)) * t470 / 0.2e1 + t582 * t66 + (t117 * t227 + t209 * t60 + t5 * t582 + t6 * t67 + t589 * t39 + (t108 + t27) * t38) * m(6) - t389 * t554 + (t1 * t172 - t14 * t71 + t15 * t72 + t2 * t389) * mrSges(7,3) + (-Ifges(7,1) * t389 + Ifges(7,4) * t172) * t551 + (-Ifges(7,5) * t389 + Ifges(7,6) * t172) * t545 + (-Ifges(7,4) * t389 + Ifges(7,2) * t172) * t550 + t4 * (-mrSges(7,1) * t172 - mrSges(7,2) * t389) - (-mrSges(6,3) * t5 + t558) * t391 - t235 * t448 / 0.2e1 + t570 * t464 + t561 * t245 + t268 * t135 + t297 * t229 + t298 * t228 + t84 * (-mrSges(5,1) * t290 + mrSges(5,2) * t292) + t227 * t98 + t209 * t41 + t184 * t171 + t185 * t170 + Ifges(3,3) * qJDD(2) + t172 * t10 / 0.2e1 + t27 * t145 + t28 * (-mrSges(7,1) * t72 + mrSges(7,2) * t71) + t25 * t70 + t71 * t51 / 0.2e1 + t67 * t65 + t63 * t16 + t30 * t20 + t31 * t21 + t581 * t280 + t587 * t108 + t589 * t144 - t345 * t595 + (t586 * t444 + t602) * qJD(3) + t596 * t86 + t597 * t85 + (t1 * t31 + t2 * t30 + t4 * t63 + (-t108 + t25) * t28 + t597 * t15 + t596 * t14) * m(7); (Ifges(7,1) * t204 + Ifges(7,4) * t203) * t536 + (-Ifges(7,1) * t385 - Ifges(7,4) * t386) * t535 + t457 + (-Ifges(6,1) * t234 + Ifges(6,5) * t451) * t530 + (-Ifges(6,5) * t234 + Ifges(6,3) * t451) * t524 + (t117 * t234 + t39 * t451) * mrSges(6,2) + (-pkin(3) * t84 + qJ(4) * t574 - t100 * t130 - t101 * t131) * m(5) + t574 * mrSges(5,3) - t578 * t50 / 0.2e1 + (mrSges(7,1) * t578 - mrSges(7,2) * t577) * t28 + (-t1 * t491 + t14 * t577 - t15 * t578 - t2 * t490) * mrSges(7,3) - t10 * t491 / 0.2e1 + t490 * t554 + t204 * t546 + (-Ifges(4,2) * t451 + t236 + t323) * t419 + (-Ifges(7,4) * t385 - Ifges(7,2) * t386) * t537 + (Ifges(7,4) * t204 + Ifges(7,2) * t203) * t538 + t84 * t409 - t38 * mrSges(6,1) * t451 + t234 * t541 + t158 * t426 / 0.2e1 + t557 * t296 + t558 * t302 + (t431 - t276) * t186 + t509 * t528 + (-g(1) * t169 - g(2) * t167 - g(3) * t222 - t302 * t5 - t303 * t6 + t609 * t38 + t39 * t610) * mrSges(6,3) + (mrSges(6,2) * t60 + t4 * t406 + t403 * t545 + t404 * t550 + t405 * t551 + t438 * t51 + 0.2e1 * t552) * t303 + (t432 + t564) * t187 + t586 * t350 * t440 + (-t117 * t156 + t232 * t5 - t337 * t60 + t38 * t584 + t39 * t585 + t390 * t6) * m(6) + (t1 * t134 + t133 * t2 + t14 * t591 + t15 * t592 + t28 * t583 - t390 * t4) * m(7) - t593 * t390 + (t598 * (-t221 * t337 + t222 * t514) + t572 * t222 + t568 * t221) * g(3) + (t598 * (-t168 * t337 + t169 * t514) + t572 * t169 + t568 * t168) * g(1) + (t598 * (-t166 * t337 + t167 * t514) + t572 * t167 + t568 * t166) * g(2) + (-Ifges(6,4) * t234 + Ifges(6,6) * t451) * t532 - t570 * t355 + t508 * t527 + t594 - t595 + (-Ifges(7,5) * t385 - Ifges(7,6) * t386) * t533 + (Ifges(7,5) * t204 + Ifges(7,6) * t203) * t534 - t337 * t41 + t232 * t66 - t131 * t217 - t130 * t218 - t156 * t98 + t133 * t20 + t134 * t21 - pkin(3) * t135 + (qJ(4) * t170 + t159 * t419 + t540 + (m(5) * t101 + t217) * qJD(4) + t599) * t344 + t583 * t70 + t584 * t145 + t585 * t144 + (-qJ(4) * t171 + t539 + (-m(5) * t100 - t218) * qJD(4) + t600) * t341 - (Ifges(6,1) * t529 + Ifges(6,4) * t531 + Ifges(6,5) * t523 + t458 + t541 + t590) * t295 + t591 * t86 + t592 * t85 + (-Ifges(6,4) * t530 + t559 + t601) * t233 + (t235 * t444 - t602) * qJD(2); t352 * t20 + t348 * t21 - t257 * t217 + t258 * t218 + t587 * t392 + t399 * qJD(6) + (-t144 - t399) * t433 + (t1 * t348 - t392 * t28 + t2 * t352 + t379 + t175 * (-t14 * t348 + t15 * t352)) * m(7) + (t38 * t392 - t39 * t433 + t379 + t60) * m(6) + (t100 * t258 - t101 * t257 + t379 + t84) * m(5) + t588; (t569 * t103 - t608 * (-t167 * t338 + t223 * t339)) * g(2) + (t569 * t105 - t608 * (-t169 * t338 + t224 * t339)) * g(1) + (t569 * t150 - t608 * (-t222 * t338 + t291 * t339)) * g(3) + (t388 + t458) * qJD(6) + (t174 + t92) * t532 + ((-t466 + t498) * t15 + (-t465 + t497) * t14 + t573) * mrSges(7,3) + (-t86 * t465 - t85 * t466 + m(7) * ((-t14 * t352 - t15 * t348) * qJD(6) + t573) + t352 * t21 - t348 * t20) * pkin(11) + (Ifges(7,2) * t352 + t505) * t550 + (Ifges(7,1) * t348 + t504) * t551 + t348 * t554 + (Ifges(7,5) * t348 + Ifges(7,6) * t352) * t545 + t497 * t546 + t498 * t547 + (t140 * t404 + t141 * t405 + t175 * t403) * qJD(6) / 0.2e1 + (-pkin(5) * t4 - t14 * t18 - t15 * t19) * m(7) + t50 * t438 + t4 * t407 + (-t507 + t49) * t530 + t559 * t392 + (t512 + t566) * t39 + t10 * t522 + t91 * t529 + t461 + (t513 - t144) * t38 - t19 * t85 - t18 * t86 - pkin(5) * t16 - t5 * mrSges(6,2) + t6 * mrSges(6,1) + (Ifges(6,1) * t530 + Ifges(6,5) * t524 + t403 * t534 + t404 * t538 + t405 * t536 - t388 - t590) * t433; -t28 * (mrSges(7,1) * t141 + mrSges(7,2) * t140) + (Ifges(7,1) * t140 - t506) * t536 + t50 * t535 + (Ifges(7,5) * t140 - Ifges(7,6) * t141) * t534 - t14 * t85 + t15 * t86 - g(1) * ((-t105 * t348 + t168 * t352) * mrSges(7,1) + (-t105 * t352 - t168 * t348) * mrSges(7,2)) - g(2) * ((-t103 * t348 + t166 * t352) * mrSges(7,1) + (-t103 * t352 - t166 * t348) * mrSges(7,2)) - g(3) * ((-t150 * t348 + t215) * mrSges(7,1) + (-t150 * t352 - t496) * mrSges(7,2)) + (t14 * t140 + t141 * t15) * mrSges(7,3) + t9 + (-Ifges(7,2) * t141 + t136 + t51) * t538 + t565;];
tau  = t7;
