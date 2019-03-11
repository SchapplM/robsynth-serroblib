% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR3
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:20
% EndTime: 2019-03-08 22:03:20
% DurationCPUTime: 35.03s
% Computational Cost: add. (15939->911), mult. (44071->1291), div. (0->0), fcn. (37771->16), ass. (0->435)
t325 = cos(pkin(7));
t333 = cos(qJ(3));
t442 = qJD(3) * t333;
t416 = t325 * t442;
t306 = pkin(2) * t416;
t322 = sin(pkin(7));
t329 = sin(qJ(3));
t506 = -pkin(9) - qJ(4);
t410 = t506 * t329;
t206 = t306 + (qJD(3) * t410 + qJD(4) * t333) * t322;
t463 = t325 * t329;
t314 = pkin(2) * t463;
t411 = t506 * t322;
t441 = qJD(4) * t322;
t207 = -t329 * t441 + (t333 * t411 - t314) * qJD(3);
t321 = sin(pkin(13));
t488 = cos(pkin(13));
t127 = t206 * t488 + t321 * t207;
t330 = sin(qJ(2));
t455 = t330 * t333;
t334 = cos(qJ(2));
t456 = t329 * t334;
t362 = -t325 * t455 - t456;
t323 = sin(pkin(6));
t447 = qJD(1) * t323;
t238 = t362 * t447;
t453 = t333 * t334;
t457 = t329 * t330;
t360 = -t325 * t457 + t453;
t239 = t360 * t447;
t174 = t321 * t238 + t239 * t488;
t580 = -t174 + t127;
t352 = -t329 * t321 + t333 * t488;
t340 = t322 * t352;
t249 = qJD(2) * t340;
t566 = t249 - qJD(5);
t462 = t325 * t333;
t315 = pkin(2) * t462;
t227 = pkin(3) * t325 + t322 * t410 + t315;
t469 = t322 * t333;
t279 = pkin(9) * t469 + t314;
t248 = qJ(4) * t469 + t279;
t169 = t321 * t227 + t488 * t248;
t152 = pkin(10) * t325 + t169;
t402 = t329 * t488;
t283 = -t333 * t321 - t402;
t258 = t283 * t322;
t251 = qJD(3) * t258;
t252 = qJD(3) * t340;
t443 = qJD(3) * t329;
t418 = t322 * t443;
t396 = pkin(3) * t418;
t176 = -pkin(4) * t251 - pkin(10) * t252 + t396;
t515 = pkin(3) * t333;
t320 = pkin(2) + t515;
t288 = t320 * t322;
t184 = -pkin(4) * t340 + pkin(10) * t258 - t288;
t328 = sin(qJ(5));
t332 = cos(qJ(5));
t421 = t330 * t447;
t393 = t322 * t421;
t439 = qJD(5) * t332;
t440 = qJD(5) * t328;
t585 = -t152 * t440 + t184 * t439 + t580 * t332 + (t176 - t393) * t328;
t517 = pkin(3) * t321;
t318 = pkin(10) + t517;
t444 = qJD(2) * t333;
t419 = t322 * t444;
t445 = qJD(2) * t322;
t250 = -t321 * t419 - t402 * t445;
t420 = t329 * t445;
t397 = pkin(3) * t420;
t175 = -pkin(4) * t250 - pkin(10) * t249 + t397;
t293 = qJD(2) * pkin(2) + t334 * t447;
t326 = cos(pkin(6));
t446 = qJD(1) * t326;
t281 = pkin(9) * t445 + t421;
t475 = t281 * t333;
t180 = t293 * t463 + t475 + (qJ(4) * t444 + t329 * t446) * t322;
t160 = t321 * t180;
t422 = t322 * t446;
t297 = t333 * t422;
t448 = t293 * t462 + t297;
t179 = (-qJ(4) * t445 - t281) * t329 + t448;
t91 = t179 * t488 - t160;
t53 = t328 * t175 + t332 * t91;
t608 = -t318 * t440 - t53;
t607 = (t207 - t238) * t488 + (-t206 + t239) * t321;
t606 = m(6) + m(7);
t605 = Ifges(4,3) + Ifges(5,3);
t604 = pkin(11) * t251 - t585;
t222 = -t258 * t328 - t332 * t325;
t166 = -qJD(5) * t222 + t252 * t332;
t223 = -t258 * t332 + t325 * t328;
t167 = qJD(5) * t223 + t252 * t328;
t603 = pkin(5) * t167 - pkin(11) * t166 - t607;
t400 = t488 * t180;
t90 = t179 * t321 + t400;
t602 = -t90 - t566 * (pkin(5) * t328 - pkin(11) * t332);
t601 = -pkin(11) * t250 - t608;
t312 = qJD(2) * t325 + qJD(3);
t204 = t250 * t328 + t312 * t332;
t369 = t250 * t332 - t312 * t328;
t503 = mrSges(5,3) * t250;
t577 = mrSges(5,1) * t312 + mrSges(6,1) * t204 + mrSges(6,2) * t369 + t503;
t363 = t325 * t453 - t457;
t215 = t323 * t363 + t326 * t469;
t600 = t312 / 0.2e1;
t327 = sin(qJ(6));
t331 = cos(qJ(6));
t384 = mrSges(7,1) * t327 + mrSges(7,2) * t331;
t599 = t384 - mrSges(5,2);
t460 = t327 * t332;
t185 = -t249 * t460 - t250 * t331;
t597 = t327 * t439 + t185;
t243 = -t352 * t445 + qJD(5);
t141 = t243 * t331 + t327 * t369;
t142 = t243 * t327 - t331 * t369;
t203 = qJD(6) - t204;
t54 = t142 * Ifges(7,5) + t141 * Ifges(7,6) + t203 * Ifges(7,3);
t497 = Ifges(6,4) * t369;
t586 = t243 * Ifges(6,6);
t587 = t204 * Ifges(6,2);
t97 = -t497 + t586 + t587;
t567 = -t97 / 0.2e1 + t54 / 0.2e1;
t385 = -mrSges(7,1) * t331 + mrSges(7,2) * t327;
t351 = m(7) * pkin(5) - t385;
t551 = -mrSges(6,1) - t351;
t431 = m(7) * pkin(11) + mrSges(7,3);
t562 = mrSges(6,2) - t431;
t596 = mrSges(4,1) * t463 + mrSges(4,2) * t462 + mrSges(3,2) + (-m(4) * pkin(9) + t328 * t551 - t332 * t562 - mrSges(4,3) - mrSges(5,3)) * t322;
t435 = qJD(2) * qJD(3);
t267 = (qJDD(2) * t333 - t329 * t435) * t322;
t268 = (qJDD(2) * t329 + t333 * t435) * t322;
t198 = t267 * t488 - t321 * t268;
t197 = qJDD(5) - t198;
t530 = t197 / 0.2e1;
t199 = t321 * t267 + t268 * t488;
t311 = qJDD(2) * t325 + qJDD(3);
t104 = qJD(5) * t369 - t199 * t328 + t311 * t332;
t535 = t104 / 0.2e1;
t103 = qJD(5) * t204 + t199 * t332 + t311 * t328;
t536 = t103 / 0.2e1;
t546 = Ifges(6,1) * t536 + Ifges(6,4) * t535 + Ifges(6,5) * t530;
t595 = -pkin(4) * t606 + t328 * t562 + t332 * t551 - mrSges(5,1);
t45 = qJD(6) * t141 + t103 * t331 + t197 * t327;
t545 = t45 / 0.2e1;
t46 = -qJD(6) * t142 - t103 * t327 + t197 * t331;
t544 = t46 / 0.2e1;
t102 = qJDD(6) - t104;
t537 = t102 / 0.2e1;
t520 = t311 / 0.2e1;
t582 = t332 * t152 + t328 * t184;
t71 = -pkin(11) * t340 + t582;
t168 = t227 * t488 - t321 * t248;
t151 = -t325 * pkin(4) - t168;
t86 = t222 * pkin(5) - t223 * pkin(11) + t151;
t32 = t327 * t86 + t331 * t71;
t593 = -qJD(6) * t32 + t327 * t604 + t331 * t603;
t31 = -t327 * t71 + t331 * t86;
t592 = qJD(6) * t31 + t327 * t603 - t331 * t604;
t17 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t65 = mrSges(6,1) * t197 - mrSges(6,3) * t103;
t505 = t17 - t65;
t591 = Ifges(4,5) * t268;
t590 = Ifges(5,5) * t199;
t589 = Ifges(4,6) * t267;
t588 = Ifges(5,6) * t198;
t501 = mrSges(6,3) * t369;
t150 = mrSges(6,1) * t243 + t501;
t73 = -mrSges(7,1) * t141 + mrSges(7,2) * t142;
t489 = t150 - t73;
t423 = t488 * pkin(3);
t319 = -t423 - pkin(4);
t280 = -t332 * pkin(5) - t328 * pkin(11) + t319;
t218 = t280 * t331 - t318 * t460;
t584 = qJD(6) * t218 + t327 * t602 - t331 * t601;
t454 = t331 * t332;
t219 = t280 * t327 + t318 * t454;
t583 = -qJD(6) * t219 + t327 * t601 + t331 * t602;
t436 = qJD(6) * t331;
t579 = t328 * t436 + t597;
t186 = t249 * t454 - t250 * t327;
t437 = qJD(6) * t328;
t578 = t327 * t437 - t331 * t439 + t186;
t576 = -t279 * qJD(3) - t238;
t575 = -pkin(9) * t418 - t239 + t306;
t361 = t325 * t456 + t455;
t472 = t322 * t329;
t216 = t323 * t361 + t326 * t472;
t133 = t321 * t215 + t216 * t488;
t466 = t323 * t334;
t273 = -t322 * t466 + t325 * t326;
t254 = t273 * t332;
t574 = -t133 * t328 + t254;
t573 = t396 - t393;
t572 = t311 * t605 + t588 + t589 + t590 + t591;
t483 = t249 * t328;
t571 = t440 - t483;
t20 = mrSges(7,1) * t102 - mrSges(7,3) * t45;
t21 = -mrSges(7,2) * t102 + mrSges(7,3) * t46;
t570 = -t327 * t20 + t331 * t21;
t308 = t325 * t446;
t211 = qJD(4) + t308 + (-pkin(3) * t444 - t293) * t322;
t123 = -pkin(4) * t249 + pkin(10) * t250 + t211;
t409 = qJD(2) * t447;
t304 = t334 * t409;
t434 = qJDD(1) * t323;
t270 = t330 * t434 + t304;
t240 = pkin(9) * qJDD(2) * t322 + t270;
t357 = t293 * t325 + t422;
t408 = qJD(2) * t441;
t391 = t330 * t409;
t269 = t334 * t434 - t391;
t253 = qJDD(2) * pkin(2) + t269;
t433 = qJDD(1) * t326;
t407 = t322 * t433;
t451 = t253 * t462 + t333 * t407;
t61 = -t281 * t442 + pkin(3) * t311 - qJ(4) * t268 + (-qJD(3) * t357 - t240 - t408) * t329 + t451;
t92 = qJD(3) * t297 + t333 * t240 + t253 * t463 - t281 * t443 + t293 * t416 + t329 * t407;
t69 = qJ(4) * t267 + t333 * t408 + t92;
t28 = t321 * t61 + t488 * t69;
t24 = pkin(10) * t311 + t28;
t212 = -t253 * t322 + t325 * t433;
t183 = -pkin(3) * t267 + qJDD(4) + t212;
t72 = -pkin(4) * t198 - pkin(10) * t199 + t183;
t148 = pkin(3) * t312 + t179;
t83 = t321 * t148 + t400;
t77 = pkin(10) * t312 + t83;
t7 = t123 * t439 + t332 * t24 + t328 * t72 - t440 * t77;
t42 = t123 * t328 + t332 * t77;
t8 = -qJD(5) * t42 - t24 * t328 + t332 * t72;
t569 = -t328 * t8 + t332 * t7;
t27 = -t321 * t69 + t488 * t61;
t23 = -t311 * pkin(4) - t27;
t14 = -t104 * pkin(5) - t103 * pkin(11) + t23;
t36 = pkin(11) * t243 + t42;
t82 = t148 * t488 - t160;
t76 = -t312 * pkin(4) - t82;
t47 = -t204 * pkin(5) + pkin(11) * t369 + t76;
t15 = -t327 * t36 + t331 * t47;
t5 = pkin(11) * t197 + t7;
t1 = qJD(6) * t15 + t14 * t327 + t331 * t5;
t16 = t327 * t47 + t331 * t36;
t2 = -qJD(6) * t16 + t14 * t331 - t327 * t5;
t568 = t1 * t331 - t2 * t327;
t241 = -t293 * t322 + t308;
t565 = ((Ifges(4,5) * t333 - Ifges(4,6) * t329) * t600 + t241 * (mrSges(4,1) * t329 + mrSges(4,2) * t333)) * t322;
t564 = 0.2e1 * t520;
t386 = t332 * mrSges(6,1) - t328 * mrSges(6,2);
t563 = -t328 * t431 - t332 * t351 - mrSges(5,1) - t386;
t260 = t283 * t325;
t137 = t323 * (-t260 * t334 + t330 * t352) - t258 * t326;
t190 = t329 * t357 + t475;
t34 = -qJD(5) * t582 - t127 * t328 + t176 * t332;
t561 = pkin(10) * t606 + t599;
t560 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t559 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t558 = -mrSges(4,1) * t333 + mrSges(4,2) * t329;
t555 = mrSges(6,3) + t599;
t554 = t76 * mrSges(6,1) + t15 * mrSges(7,1) - t16 * mrSges(7,2);
t529 = -t203 / 0.2e1;
t532 = -t142 / 0.2e1;
t534 = -t141 / 0.2e1;
t553 = Ifges(7,5) * t532 + Ifges(7,6) * t534 + Ifges(7,3) * t529;
t552 = -m(4) * pkin(2) - mrSges(3,1) + t558;
t93 = -qJD(3) * t190 - t240 * t329 + t451;
t550 = t93 * mrSges(4,1) + t27 * mrSges(5,1) - t92 * mrSges(4,2) - t28 * mrSges(5,2);
t549 = t322 ^ 2;
t335 = qJD(2) ^ 2;
t12 = t45 * Ifges(7,4) + t46 * Ifges(7,2) + t102 * Ifges(7,6);
t548 = t12 / 0.2e1;
t547 = Ifges(7,1) * t545 + Ifges(7,4) * t544 + Ifges(7,5) * t537;
t494 = Ifges(7,4) * t142;
t55 = t141 * Ifges(7,2) + t203 * Ifges(7,6) + t494;
t542 = -t55 / 0.2e1;
t541 = t55 / 0.2e1;
t134 = Ifges(7,4) * t141;
t56 = t142 * Ifges(7,1) + t203 * Ifges(7,5) + t134;
t540 = -t56 / 0.2e1;
t539 = t56 / 0.2e1;
t533 = t141 / 0.2e1;
t531 = t142 / 0.2e1;
t528 = t203 / 0.2e1;
t527 = -t204 / 0.2e1;
t526 = t369 / 0.2e1;
t525 = -t369 / 0.2e1;
t524 = -t243 / 0.2e1;
t522 = -t250 / 0.2e1;
t518 = pkin(2) * t322;
t514 = g(3) * t323;
t6 = -pkin(5) * t197 - t8;
t511 = t328 * t6;
t504 = mrSges(5,3) * t249;
t502 = mrSges(6,3) * t204;
t500 = Ifges(4,4) * t329;
t499 = Ifges(4,4) * t333;
t498 = Ifges(5,4) * t250;
t496 = Ifges(6,4) * t328;
t495 = Ifges(6,4) * t332;
t493 = Ifges(7,4) * t327;
t492 = Ifges(7,4) * t331;
t487 = sin(pkin(12));
t485 = t204 * t327;
t484 = t204 * t331;
t482 = t249 * t332;
t480 = t273 * t328;
t398 = t487 * t334;
t324 = cos(pkin(12));
t465 = t324 * t330;
t275 = t326 * t465 + t398;
t478 = t275 * t329;
t399 = t487 * t330;
t464 = t324 * t334;
t277 = -t326 * t399 + t464;
t476 = t277 * t329;
t471 = t322 * t330;
t468 = t323 * t324;
t467 = t323 * t330;
t461 = t327 * t328;
t459 = t328 * t331;
t272 = pkin(3) * t463 + t411;
t274 = t326 * t464 - t399;
t450 = -t275 * t272 + t274 * t320;
t276 = -t326 * t398 - t465;
t449 = -t277 * t272 + t276 * t320;
t438 = qJD(6) * t327;
t11 = Ifges(7,5) * t45 + Ifges(7,6) * t46 + Ifges(7,3) * t102;
t432 = pkin(3) * t462;
t430 = t322 * t468;
t429 = t322 * t467;
t426 = Ifges(6,5) * t103 + Ifges(6,6) * t104 + Ifges(6,3) * t197;
t417 = t322 * t442;
t414 = t318 * t439;
t412 = t472 / 0.2e1;
t406 = t442 / 0.2e1;
t405 = t439 / 0.2e1;
t404 = -t437 / 0.2e1;
t403 = t323 * t487;
t121 = -t198 * mrSges(5,1) + t199 * mrSges(5,2);
t395 = mrSges(4,3) * t420;
t394 = mrSges(4,3) * t419;
t392 = qJD(2) * t429;
t389 = t322 * t403;
t387 = -t272 * t467 + t320 * t466;
t383 = Ifges(6,1) * t332 - t496;
t382 = Ifges(7,1) * t331 - t493;
t381 = Ifges(7,1) * t327 + t492;
t380 = -Ifges(6,2) * t328 + t495;
t379 = -Ifges(7,2) * t327 + t492;
t378 = Ifges(7,2) * t331 + t493;
t377 = Ifges(6,5) * t332 - Ifges(6,6) * t328;
t376 = Ifges(7,5) * t331 - Ifges(7,6) * t327;
t375 = Ifges(7,5) * t327 + Ifges(7,6) * t331;
t41 = t123 * t332 - t328 * t77;
t52 = t175 * t332 - t328 * t91;
t112 = t133 * t332 + t480;
t132 = -t215 * t488 + t216 * t321;
t59 = t112 * t331 + t132 * t327;
t58 = -t112 * t327 + t132 * t331;
t87 = -t152 * t328 + t184 * t332;
t178 = t223 * t331 - t327 * t340;
t177 = -t223 * t327 - t331 * t340;
t367 = -pkin(3) * t476 + t276 * t432 + t389 * t515;
t364 = -t274 * t325 + t430;
t35 = -pkin(5) * t243 - t41;
t359 = t35 * t384;
t356 = t215 * pkin(3);
t354 = t558 * t322;
t353 = (Ifges(4,2) * t333 + t500) * t322;
t220 = -t274 * t322 - t325 * t468;
t221 = -t276 * t322 + t325 * t403;
t349 = -g(1) * t221 - g(2) * t220 - g(3) * t273;
t345 = t329 * t549 * (Ifges(4,1) * t333 - t500);
t342 = t276 * t325 + t389;
t113 = -t258 * t468 + t260 * t274 - t275 * t352;
t339 = t274 * t432 + (-t333 * t430 - t478) * pkin(3);
t116 = t258 * t403 + t276 * t260 - t277 * t352;
t337 = (-t15 * t331 - t16 * t327) * qJD(6) + t568;
t305 = Ifges(4,4) * t419;
t278 = -pkin(9) * t472 + t315;
t264 = qJD(2) * t354;
t263 = -mrSges(4,2) * t312 + t394;
t262 = mrSges(4,1) * t312 - t395;
t259 = t352 * t325;
t242 = Ifges(5,4) * t249;
t229 = Ifges(4,1) * t420 + t312 * Ifges(4,5) + t305;
t228 = t312 * Ifges(4,6) + qJD(2) * t353;
t225 = mrSges(4,1) * t311 - mrSges(4,3) * t268;
t224 = -mrSges(4,2) * t311 + mrSges(4,3) * t267;
t213 = -mrSges(5,2) * t312 + t504;
t201 = Ifges(6,4) * t204;
t200 = -mrSges(4,1) * t267 + mrSges(4,2) * t268;
t193 = (t260 * t330 + t334 * t352) * t323;
t192 = -t259 * t467 + t283 * t466;
t189 = -t281 * t329 + t448;
t187 = -mrSges(5,1) * t249 - mrSges(5,2) * t250;
t182 = mrSges(5,1) * t311 - mrSges(5,3) * t199;
t181 = -mrSges(5,2) * t311 + mrSges(5,3) * t198;
t172 = -t250 * Ifges(5,1) + t312 * Ifges(5,5) + t242;
t171 = t249 * Ifges(5,2) + t312 * Ifges(5,6) - t498;
t165 = t326 * t417 + (qJD(2) * t360 + qJD(3) * t363) * t323;
t164 = -t326 * t418 + (qJD(2) * t362 - qJD(3) * t361) * t323;
t149 = -mrSges(6,2) * t243 + t502;
t145 = -t259 * t277 + t276 * t283;
t143 = -t259 * t275 + t274 * t283;
t136 = t340 * t326 + (t259 * t334 + t283 * t330) * t323;
t129 = t174 * t328 - t332 * t393;
t128 = -pkin(5) * t369 - pkin(11) * t204;
t120 = t137 * t332 + t480;
t117 = t276 * t259 + t277 * t283 + t340 * t403;
t114 = t259 * t274 + t275 * t283 - t340 * t468;
t98 = -Ifges(6,1) * t369 + t243 * Ifges(6,5) + t201;
t96 = -Ifges(6,5) * t369 + t204 * Ifges(6,6) + t243 * Ifges(6,3);
t95 = mrSges(7,1) * t203 - mrSges(7,3) * t142;
t94 = -mrSges(7,2) * t203 + mrSges(7,3) * t141;
t85 = t321 * t164 + t165 * t488;
t84 = -t164 * t488 + t165 * t321;
t81 = -t116 * t332 + t221 * t328;
t79 = -t113 * t332 + t220 * t328;
t70 = pkin(5) * t340 - t87;
t68 = -qJD(6) * t178 - t166 * t327 - t251 * t331;
t67 = qJD(6) * t177 + t166 * t331 - t251 * t327;
t66 = -mrSges(6,2) * t197 + mrSges(6,3) * t104;
t50 = -mrSges(6,1) * t104 + mrSges(6,2) * t103;
t48 = pkin(5) * t250 - t52;
t40 = qJD(5) * t112 + t328 * t85 - t332 * t392;
t39 = qJD(5) * t574 + t328 * t392 + t332 * t85;
t37 = t103 * Ifges(6,4) + t104 * Ifges(6,2) + t197 * Ifges(6,6);
t30 = pkin(5) * t251 - t34;
t26 = t128 * t327 + t331 * t41;
t25 = t128 * t331 - t327 * t41;
t10 = -qJD(6) * t59 - t327 * t39 + t331 * t84;
t9 = qJD(6) * t58 + t327 * t84 + t331 * t39;
t3 = [-t505 * t574 + m(7) * (t1 * t59 + t10 * t15 + t16 * t9 + t2 * t58 + t35 * t40 - t574 * t6) + m(6) * (t112 * t7 + t132 * t23 + t39 * t42 - t40 * t41 + t574 * t8 + t76 * t84) - t577 * t84 + (-m(2) - m(3) - m(4) - m(5) - t606) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t335) * t334 + (-mrSges(3,1) * t335 - mrSges(3,2) * qJDD(2) + (t187 + t264) * t445) * t330) * t323 - t489 * t40 + (-t182 + t50) * t132 + t133 * t181 + t85 * t213 + (t200 + t121) * t273 + t39 * t149 + t112 * t66 + t9 * t94 + t10 * t95 + t58 * t20 + t59 * t21 + m(3) * (qJDD(1) * t326 ^ 2 + (t269 * t334 + t270 * t330) * t323) + m(4) * (t164 * t189 + t165 * t190 + t212 * t273 + t215 * t93 + t216 * t92 + t241 * t392) + m(5) * (-t132 * t27 + t133 * t28 + t183 * t273 + t211 * t392 - t82 * t84 + t83 * t85) + t216 * t224 + t215 * t225 + t164 * t262 + t165 * t263 + m(2) * qJDD(1); (m(5) * t82 - m(6) * t76 + t577) * t607 + (Ifges(7,4) * t67 + Ifges(7,2) * t68) * t533 + (Ifges(7,4) * t178 + Ifges(7,2) * t177) * t544 + ((Ifges(4,5) * t329 + Ifges(4,6) * t333) * t520 + t268 * (Ifges(4,1) * t329 + t499) / 0.2e1 + t229 * t406) * t322 + (Ifges(7,1) * t67 + Ifges(7,4) * t68) * t531 + (Ifges(7,1) * t178 + Ifges(7,4) * t177) * t545 + (t330 * t514 - t270 + t304) * mrSges(3,2) - (Ifges(6,3) * t530 + Ifges(6,6) * t535 + Ifges(6,5) * t536 + t426 / 0.2e1 + t183 * mrSges(5,1) - Ifges(5,4) * t199 - Ifges(5,2) * t198 - t28 * mrSges(5,3) - t564 * Ifges(5,6) + t560) * t340 + (-m(5) * t450 + t552 * t274 - t606 * (-pkin(10) * t143 + t450) + t555 * t143 + t595 * (t260 * t275 + t274 * t352) + t596 * t275) * g(2) + (-m(5) * t449 + t552 * t276 - t606 * (-pkin(10) * t145 + t449) + t555 * t145 + t595 * (t260 * t277 + t276 * t352) + t596 * t277) * g(1) + (-m(5) * t387 - t193 * mrSges(5,1) - mrSges(5,3) * t429 - t606 * (t193 * pkin(4) - pkin(10) * t192 + t387) + t551 * (t193 * t332 + t328 * t429) + t555 * t192 + t562 * (t193 * t328 - t332 * t429)) * g(3) + (t151 * t23 + t7 * t582 + t8 * t87 + t585 * t42 + (t129 + t34) * t41) * m(6) + t345 * t435 / 0.2e1 + t67 * t539 + t68 * t541 + t178 * t547 - t228 * t418 / 0.2e1 + (-t189 * t417 - t190 * t418 + t469 * t92 - t471 * t514 - t472 * t93) * mrSges(4,3) + (-Ifges(6,6) * t530 - Ifges(6,2) * t535 - Ifges(6,4) * t536 + Ifges(7,3) * t537 + Ifges(7,6) * t544 + Ifges(7,5) * t545 + t11 / 0.2e1 - t37 / 0.2e1 + t23 * mrSges(6,1) - t7 * mrSges(6,3) + t559) * t222 + (t572 / 0.2e1 + t589 / 0.2e1 + t591 / 0.2e1 + t588 / 0.2e1 + t590 / 0.2e1 + t605 * t520 + t550) * t325 + t212 * t354 + t169 * t181 + t168 * t182 + t6 * (-mrSges(7,1) * t177 + mrSges(7,2) * t178) + (t168 * t27 + t169 * t28 - t183 * t288 + t211 * t573 + t580 * t83) * m(5) + t177 * t548 - (t360 * mrSges(4,1) + t362 * mrSges(4,2)) * t514 + t267 * t353 / 0.2e1 + t166 * t98 / 0.2e1 + t34 * t150 + t151 * t50 + (t1 * t177 - t15 * t67 + t16 * t68 - t178 * t2) * mrSges(7,3) + (Ifges(5,5) * t252 + Ifges(5,6) * t251) * t600 + Ifges(3,3) * qJDD(2) + t87 * t65 + t70 * t17 + t30 * t73 + t35 * (-mrSges(7,1) * t68 + mrSges(7,2) * t67) + (t251 * t83 - t252 * t82) * mrSges(5,3) + t204 * (Ifges(6,4) * t166 - Ifges(6,6) * t251) / 0.2e1 + t243 * (Ifges(6,5) * t166 - Ifges(6,3) * t251) / 0.2e1 + (Ifges(6,1) * t166 - Ifges(6,5) * t251) * t525 + (Ifges(5,1) * t252 + Ifges(5,4) * t251) * t522 + (t166 * t76 + t251 * t42) * mrSges(6,2) + t41 * (-mrSges(6,1) * t251 - mrSges(6,3) * t166) + t211 * (-mrSges(5,1) * t251 + mrSges(5,2) * t252) + t249 * (Ifges(5,4) * t252 + Ifges(5,2) * t251) / 0.2e1 + t582 * t66 + t31 * t20 + t32 * t21 + (mrSges(6,2) * t23 - mrSges(6,3) * t8 + 0.2e1 * t546) * t223 - t264 * t393 + (Ifges(4,4) * t268 + Ifges(4,2) * t267 + Ifges(4,6) * t311) * t469 / 0.2e1 - t200 * t518 - (t183 * mrSges(5,2) - t27 * mrSges(5,3) + Ifges(5,1) * t199 + Ifges(5,4) * t198 + Ifges(5,5) * t564) * t258 + t565 * qJD(3) - t251 * t96 / 0.2e1 + t251 * t171 / 0.2e1 + t252 * t172 / 0.2e1 + t573 * t187 + t575 * t263 + t576 * t262 + (-(pkin(2) * t334 + pkin(9) * t471) * t514 - t241 * t393 - t212 * t518 + t278 * t93 + t279 * t92 + t575 * t190 + t576 * t189) * m(4) + t580 * t213 + t278 * t225 + t279 * t224 - t288 * t121 + (Ifges(4,1) * t268 + Ifges(4,4) * t267 + Ifges(4,5) * t311) * t412 + t489 * t129 + t585 * t149 + (Ifges(7,5) * t67 + Ifges(7,6) * t68) * t528 + (Ifges(7,5) * t178 + Ifges(7,6) * t177) * t537 + (-Ifges(6,4) * t525 + Ifges(7,3) * t528 + Ifges(7,5) * t531 + Ifges(7,6) * t533 - t587 / 0.2e1 - mrSges(6,3) * t42 - t586 / 0.2e1 + t554 + t567) * t167 + (-t334 * t514 + t269 + t391) * mrSges(3,1) + t549 * qJD(2) * (-Ifges(4,2) * t329 + t499) * t406 + t592 * t94 + t593 * t95 + (t1 * t32 + t2 * t31 + t6 * t70 + (-t129 + t30) * t35 + t592 * t16 + t593 * t15) * m(7); -t566 * t76 * (mrSges(6,1) * t328 + mrSges(6,2) * t332) + (t23 * t319 - t41 * t52 - t42 * t53 - t76 * t90) * m(6) + (t1 * t219 + t583 * t15 + t584 * t16 + t2 * t218 - t35 * t48) * m(7) + (t505 * t328 + t332 * t66 + ((-t328 * t42 - t332 * t41) * qJD(5) + t569) * m(6) + (t35 * t439 + t511) * m(7)) * t318 + t608 * t149 + (g(1) * t116 + g(2) * t113 - g(3) * t137 - t571 * t42 + (-t439 + t482) * t41 + t569) * mrSges(6,3) + (-m(5) * t356 - mrSges(4,1) * t215 + mrSges(4,2) * t216 - t606 * (t136 * pkin(4) + t356) - t561 * t137 + t563 * t136) * g(3) + (-(-t333 * t364 - t478) * mrSges(4,1) - (-t275 * t333 + t329 * t364) * mrSges(4,2) - m(5) * t339 - t606 * (t114 * pkin(4) + t339) + t561 * t113 + t563 * t114) * g(2) + (-(t333 * t342 - t476) * mrSges(4,1) - (-t277 * t333 - t329 * t342) * mrSges(4,2) - m(5) * t367 - t606 * (t117 * pkin(4) + t367) + t561 * t116 + t563 * t117) * g(1) + t82 * t504 + t384 * t511 + (t395 + t262) * t190 + (-t482 / 0.2e1 + t405) * t98 - t12 * t461 / 0.2e1 + t1 * (mrSges(7,2) * t332 - mrSges(7,3) * t461) + (-t414 - t52) * t150 + (Ifges(7,4) * t186 + Ifges(7,2) * t185) * t534 + (t404 * t55 + t405 * t56) * t331 + t2 * (-mrSges(7,1) * t332 - mrSges(7,3) * t459) + t181 * t517 + t171 * t522 + (-Ifges(6,3) * t250 + t249 * t377) * t524 + (-Ifges(6,5) * t250 + t249 * t383) * t526 + (-Ifges(6,6) * t250 + t249 * t380) * t527 + (-t375 * t437 + (Ifges(7,3) * t328 + t332 * t376) * qJD(5)) * t528 + (Ifges(6,5) * t328 + Ifges(6,6) * t332) * t530 + (-t381 * t437 + (Ifges(7,5) * t328 + t332 * t382) * qJD(5)) * t531 - t42 * mrSges(6,2) * t250 + t41 * mrSges(6,1) * t250 + (Ifges(7,1) * t186 + Ifges(7,4) * t185) * t532 + (Ifges(7,5) * t186 + Ifges(7,6) * t185) * t529 + (-t378 * t437 + (Ifges(7,6) * t328 + t332 * t379) * qJD(5)) * t533 + (Ifges(6,2) * t332 + t496) * t535 + (Ifges(6,1) * t328 + t495) * t536 + (-Ifges(7,3) * t332 + t328 * t376) * t537 + t186 * t540 + (-Ifges(7,6) * t332 + t328 * t379) * t544 + (-Ifges(7,5) * t332 + t328 * t382) * t545 + t328 * t546 + t459 * t547 - (-Ifges(4,2) * t420 + t229 + t305) * t419 / 0.2e1 + (Ifges(5,1) * t249 + t498 + t96) * t250 / 0.2e1 - (Ifges(5,2) * t250 + t172 + t242) * t249 / 0.2e1 + t219 * t21 + (-t48 + t414) * t73 + t182 * t423 + (t394 - t263) * t189 + t218 * t20 - t91 * t213 - t83 * t503 - t335 * t345 / 0.2e1 + t572 - t187 * t397 + (t204 * t380 + t243 * t377 - t369 * t383) * qJD(5) / 0.2e1 + t597 * t542 + (t553 - t567) * t483 + (t228 * t412 - t565) * qJD(2) - t211 * (-mrSges(5,1) * t250 + mrSges(5,2) * t249) + t567 * t440 + t550 + t327 * t56 * t404 + t577 * t90 + (mrSges(7,1) * t571 + mrSges(7,3) * t578) * t15 + (t579 * mrSges(7,1) - t578 * mrSges(7,2)) * t35 + (-mrSges(7,2) * t571 - mrSges(7,3) * t579) * t16 - t312 * (Ifges(5,5) * t249 + Ifges(5,6) * t250) / 0.2e1 + t319 * t50 + t583 * t95 + t584 * t94 - t23 * t386 + ((t27 * t488 + t28 * t321) * pkin(3) - t211 * t397 + t82 * t90 - t83 * t91) * m(5) + t332 * t37 / 0.2e1 - t332 * t11 / 0.2e1; -t185 * t95 - t186 * t94 - t249 * t213 - t577 * t250 + (-t249 * t149 + (-t327 * t95 + t331 * t94 + t149) * qJD(5) - t505) * t332 + (t66 + (-t327 * t94 - t331 * t95) * qJD(6) + t566 * t489 + t570) * t328 + t121 + (t349 + (-t6 + (-t15 * t327 + t16 * t331) * qJD(5)) * t332 + (qJD(5) * t35 + t337) * t328 - t15 * t185 - t16 * t186 - t35 * t483) * m(7) + (t250 * t76 + t328 * t7 + t332 * t8 + t349 - t566 * (-t328 * t41 + t332 * t42)) * m(6) + (-t249 * t83 - t250 * t82 + t183 + t349) * m(5); t426 + (t502 - t149) * t41 + t97 * t525 + (-pkin(5) * t6 - t15 * t25 - t16 * t26) * m(7) + t375 * t537 + t436 * t539 + t484 * t540 + t485 * t541 + t438 * t542 + t378 * t544 + t381 * t545 + t327 * t547 + t331 * t548 + (-t76 * mrSges(6,2) + Ifges(6,1) * t526 + Ifges(6,5) * t524 + t376 * t529 + t379 * t534 + t382 * t532 - t359) * t204 + (t141 * t379 + t142 * t382 + t203 * t376) * qJD(6) / 0.2e1 + qJD(6) * t359 + t560 + (t98 + t201) * t527 - t26 * t94 - t25 * t95 + (t562 * t79 + t551 * (t113 * t328 + t220 * t332)) * g(2) + (t562 * t81 + t551 * (t116 * t328 + t221 * t332)) * g(1) - (-Ifges(6,2) * t527 - Ifges(6,6) * t524 + t553 - t554) * t369 - pkin(5) * t17 + (t497 + t54) * t526 + ((-t438 + t485) * t16 + (-t436 + t484) * t15 + t568) * mrSges(7,3) + (m(7) * t337 - t436 * t95 - t438 * t94 + t570) * pkin(11) + (-m(7) * t35 + t489 - t501) * t42 + (t562 * t120 + t551 * (-t137 * t328 + t254)) * g(3) + t6 * t385; -t35 * (mrSges(7,1) * t142 + mrSges(7,2) * t141) + (Ifges(7,1) * t141 - t494) * t532 + t55 * t531 + (Ifges(7,5) * t141 - Ifges(7,6) * t142) * t529 - t15 * t94 + t16 * t95 - g(1) * ((-t117 * t331 - t327 * t81) * mrSges(7,1) + (t117 * t327 - t331 * t81) * mrSges(7,2)) - g(2) * ((-t114 * t331 - t327 * t79) * mrSges(7,1) + (t114 * t327 - t331 * t79) * mrSges(7,2)) - g(3) * ((-t120 * t327 - t136 * t331) * mrSges(7,1) + (-t120 * t331 + t136 * t327) * mrSges(7,2)) + (t141 * t15 + t142 * t16) * mrSges(7,3) + t11 + (-Ifges(7,2) * t142 + t134 + t56) * t534 + t559;];
tau  = t3;
