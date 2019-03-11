% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:39
% EndTime: 2019-03-09 10:50:30
% DurationCPUTime: 34.85s
% Computational Cost: add. (12586->856), mult. (27790->1094), div. (0->0), fcn. (20052->12), ass. (0->389)
t592 = mrSges(5,1) + mrSges(6,1);
t564 = mrSges(6,3) - mrSges(5,2);
t321 = sin(pkin(10));
t322 = cos(pkin(10));
t325 = sin(qJ(4));
t329 = cos(qJ(4));
t264 = t321 * t329 + t322 * t325;
t330 = cos(qJ(2));
t352 = t264 * t330;
t214 = qJD(1) * t352;
t241 = t264 * qJD(4);
t541 = t214 - t241;
t439 = t321 * t325;
t263 = -t329 * t322 + t439;
t351 = t263 * t330;
t215 = qJD(1) * t351;
t240 = t263 * qJD(4);
t540 = -t215 + t240;
t591 = mrSges(7,3) - mrSges(4,3);
t462 = pkin(8) + qJ(3);
t276 = t462 * t321;
t277 = t462 * t322;
t419 = qJD(4) * t329;
t423 = qJD(3) * t322;
t424 = qJD(3) * t321;
t145 = -t276 * t419 + t329 * t423 + (-qJD(4) * t277 - t424) * t325;
t326 = sin(qJ(2));
t370 = pkin(2) * t326 - qJ(3) * t330;
t266 = t370 * qJD(1);
t428 = qJD(1) * t326;
t405 = t321 * t428;
t201 = pkin(7) * t405 + t322 * t266;
t435 = t322 * t330;
t362 = pkin(3) * t326 - pkin(8) * t435;
t164 = qJD(1) * t362 + t201;
t242 = t321 * t266;
t436 = t322 * t326;
t437 = t321 * t330;
t355 = -pkin(7) * t436 - pkin(8) * t437;
t187 = qJD(1) * t355 + t242;
t101 = t325 * t164 + t329 * t187;
t91 = qJ(5) * t428 + t101;
t590 = t145 - t91;
t562 = -Ifges(5,4) + Ifges(6,5);
t589 = t562 + Ifges(6,5);
t407 = t322 * t428;
t258 = qJD(2) * t321 + t407;
t356 = -t322 * qJD(2) + t405;
t340 = t329 * t258 - t325 * t356;
t487 = t340 / 0.2e1;
t561 = Ifges(6,4) + Ifges(5,5);
t588 = t561 * t487;
t396 = -qJ(3) * t326 - pkin(1);
t274 = -pkin(2) * t330 + t396;
t246 = t274 * qJD(1);
t427 = qJD(1) * t330;
t309 = pkin(7) * t427;
t279 = qJD(2) * qJ(3) + t309;
t190 = t322 * t246 - t279 * t321;
t131 = -pkin(3) * t427 - pkin(8) * t258 + t190;
t191 = t321 * t246 + t322 * t279;
t140 = -pkin(8) * t356 + t191;
t61 = t329 * t131 - t325 * t140;
t587 = qJD(5) - t61;
t100 = t164 * t329 - t325 * t187;
t199 = -t325 * t276 + t329 * t277;
t146 = t264 * qJD(3) + qJD(4) * t199;
t586 = -t146 - t100;
t301 = qJD(4) - t427;
t475 = -t301 / 0.2e1;
t293 = qJD(6) - t301;
t477 = -t293 / 0.2e1;
t182 = t325 * t258 + t329 * t356;
t324 = sin(qJ(6));
t328 = cos(qJ(6));
t574 = t182 * t324 + t328 * t340;
t494 = -t574 / 0.2e1;
t105 = t182 * t328 - t324 * t340;
t496 = -t105 / 0.2e1;
t559 = -Ifges(5,3) - Ifges(6,2);
t585 = Ifges(7,5) * t494 + Ifges(7,6) * t496 + Ifges(7,3) * t477 + t475 * t559;
t417 = qJD(1) * qJD(2);
t270 = qJDD(1) * t326 + t330 * t417;
t216 = qJDD(2) * t322 - t270 * t321;
t217 = qJDD(2) * t321 + t270 * t322;
t89 = -qJD(4) * t182 + t325 * t216 + t329 * t217;
t501 = t89 / 0.2e1;
t90 = qJD(4) * t340 - t329 * t216 + t325 * t217;
t499 = t90 / 0.2e1;
t489 = t182 / 0.2e1;
t313 = t330 * qJDD(1);
t403 = t326 * t417;
t269 = -t313 + t403;
t262 = qJDD(4) + t269;
t480 = t262 / 0.2e1;
t400 = -t427 / 0.2e1;
t584 = -mrSges(5,3) - mrSges(6,2);
t563 = Ifges(5,1) + Ifges(6,1);
t560 = -Ifges(5,6) + Ifges(6,6);
t583 = -pkin(9) * t541 + t590;
t332 = -pkin(4) - pkin(5);
t411 = t332 * t326;
t582 = pkin(9) * t540 - qJD(1) * t411 - t586;
t581 = -pkin(9) * t340 + t587;
t580 = -m(7) * (-pkin(9) + t462) - m(4) * qJ(3) + t591;
t319 = pkin(10) + qJ(4);
t311 = sin(t319);
t312 = cos(t319);
t364 = t311 * t324 + t312 * t328;
t365 = t311 * t328 - t312 * t324;
t579 = t364 * mrSges(7,1) + t365 * mrSges(7,2) + t564 * t311 + t312 * t592;
t331 = cos(qJ(1));
t327 = sin(qJ(1));
t432 = t327 * t330;
t225 = t311 * t432 + t312 * t331;
t430 = t331 * t311;
t226 = t312 * t432 - t430;
t366 = t225 * t324 + t226 * t328;
t533 = t225 * t328 - t226 * t324;
t578 = mrSges(7,1) * t533 - t366 * mrSges(7,2);
t180 = Ifges(5,4) * t182;
t451 = Ifges(6,5) * t182;
t555 = t301 * t561 + t340 * t563 - t180 + t451;
t577 = -t555 / 0.2e1;
t408 = t321 * t427;
t251 = pkin(3) * t408 + t309;
t576 = qJ(5) * t540 - qJD(5) * t264 - t251;
t575 = g(1) * t331 + g(2) * t327;
t308 = pkin(7) * t428;
t386 = qJD(2) * pkin(2) - qJD(3) - t308;
t200 = pkin(3) * t356 - t386;
t336 = qJ(5) * t340 - t200;
t71 = t182 * pkin(4) - t336;
t179 = Ifges(6,5) * t340;
t94 = t301 * Ifges(6,6) + t182 * Ifges(6,3) + t179;
t453 = Ifges(5,4) * t340;
t97 = -t182 * Ifges(5,2) + t301 * Ifges(5,6) + t453;
t573 = mrSges(6,1) * t71 + t94 / 0.2e1 - t97 / 0.2e1;
t46 = t301 * t332 + t581;
t288 = t301 * qJ(5);
t62 = t325 * t131 + t329 * t140;
t52 = pkin(9) * t182 + t62;
t50 = t288 + t52;
t13 = -t324 * t50 + t328 * t46;
t14 = t324 * t46 + t328 * t50;
t457 = Ifges(3,4) * t326;
t554 = t330 * Ifges(3,2);
t374 = t457 + t554;
t57 = -pkin(4) * t301 + t587;
t58 = t288 + t62;
t572 = t61 * mrSges(5,1) - t57 * mrSges(6,1) - t13 * mrSges(7,1) - t62 * mrSges(5,2) + t14 * mrSges(7,2) + t58 * mrSges(6,3) + Ifges(4,5) * t258 / 0.2e1 - Ifges(4,6) * t356 / 0.2e1 + Ifges(4,3) * t400 + t588 + t560 * t489 - Ifges(3,6) * qJD(2) / 0.2e1 - qJD(1) * t374 / 0.2e1 + t585;
t421 = qJD(4) * t325;
t422 = qJD(3) * t326;
t446 = qJDD(1) * pkin(1);
t172 = pkin(2) * t269 - qJ(3) * t270 - qJD(1) * t422 - t446;
t306 = pkin(7) * t313;
t224 = qJDD(2) * qJ(3) + t306 + (qJD(3) - t308) * qJD(2);
t120 = t322 * t172 - t224 * t321;
t69 = pkin(3) * t269 - pkin(8) * t217 + t120;
t121 = t321 * t172 + t322 * t224;
t77 = pkin(8) * t216 + t121;
t16 = -t131 * t421 - t140 * t419 - t325 * t77 + t329 * t69;
t353 = qJDD(5) - t16;
t12 = -pkin(4) * t262 + t353;
t257 = t270 * pkin(7);
t235 = -qJDD(2) * pkin(2) + qJDD(3) + t257;
t161 = -t216 * pkin(3) + t235;
t25 = t90 * pkin(4) - t89 * qJ(5) - qJD(5) * t340 + t161;
t500 = -t90 / 0.2e1;
t571 = mrSges(5,2) * t161 + mrSges(6,2) * t12 - mrSges(5,3) * t16 - mrSges(6,3) * t25 + Ifges(5,4) * t500 + 0.2e1 * t480 * t561 + t499 * t589 + 0.2e1 * t501 * t563;
t23 = qJD(6) * t105 + t324 * t90 + t328 * t89;
t509 = t23 / 0.2e1;
t24 = -qJD(6) * t574 - t324 * t89 + t328 * t90;
t508 = t24 / 0.2e1;
t570 = -m(4) - m(3);
t569 = -m(7) - m(6);
t485 = t216 / 0.2e1;
t484 = t217 / 0.2e1;
t250 = qJDD(6) - t262;
t481 = t250 / 0.2e1;
t568 = -t269 / 0.2e1;
t478 = t269 / 0.2e1;
t566 = mrSges(6,3) * t71;
t565 = qJD(2) / 0.2e1;
t198 = t329 * t276 + t277 * t325;
t155 = -pkin(9) * t264 + t198;
t156 = pkin(9) * t263 + t199;
t74 = t155 * t324 + t156 * t328;
t558 = -qJD(6) * t74 - t324 * t583 + t328 * t582;
t73 = t155 * t328 - t156 * t324;
t557 = qJD(6) * t73 + t324 * t582 + t328 * t583;
t383 = mrSges(3,1) * t330 - mrSges(3,2) * t326;
t552 = -mrSges(2,1) - t383;
t551 = -t332 * t541 - t576;
t272 = t328 * qJ(5) + t324 * t332;
t550 = -qJD(6) * t272 - t324 * t581 - t328 * t52;
t271 = -t324 * qJ(5) + t328 * t332;
t549 = qJD(6) * t271 - t324 * t52 + t328 * t581;
t381 = -mrSges(4,1) * t322 + mrSges(4,2) * t321;
t350 = m(4) * pkin(2) - t381;
t545 = t330 * t350;
t544 = -pkin(4) * t541 + t576;
t459 = mrSges(5,3) * t182;
t147 = -mrSges(5,2) * t301 - t459;
t461 = mrSges(6,2) * t182;
t150 = mrSges(6,3) * t301 - t461;
t543 = t147 + t150;
t458 = mrSges(5,3) * t340;
t148 = mrSges(5,1) * t301 - t458;
t460 = mrSges(6,2) * t340;
t149 = -mrSges(6,1) * t301 + t460;
t542 = -t149 + t148;
t539 = qJD(2) * mrSges(3,1) - mrSges(4,1) * t356 - t258 * mrSges(4,2) - mrSges(3,3) * t428;
t316 = t330 * pkin(4);
t447 = qJ(5) * t330;
t538 = t311 * t447 + t312 * t316;
t537 = t326 * t580 - t545;
t454 = Ifges(4,4) * t322;
t373 = -Ifges(4,2) * t321 + t454;
t455 = Ifges(4,4) * t321;
t375 = Ifges(4,1) * t322 - t455;
t536 = t258 * (Ifges(4,5) * t326 + t330 * t375) - (Ifges(4,6) * t326 + t330 * t373) * t356;
t535 = -t262 * t559 + t560 * t90 + t561 * t89;
t256 = -pkin(7) * t403 + t306;
t534 = t256 * t330 + t257 * t326;
t532 = -t120 * t321 + t121 * t322;
t531 = t584 * t326;
t529 = -m(5) + t569;
t380 = mrSges(4,1) * t321 + mrSges(4,2) * t322;
t528 = -t380 - mrSges(3,3) + mrSges(2,2);
t526 = t386 * t330 * t380 - t191 * (-mrSges(4,2) * t326 - mrSges(4,3) * t437) - t190 * (mrSges(4,1) * t326 - mrSges(4,3) * t435);
t55 = t182 * t332 + t336;
t525 = -mrSges(7,1) * t55 + mrSges(7,3) * t14;
t524 = mrSges(7,2) * t55 - mrSges(7,3) * t13;
t307 = Ifges(3,4) * t427;
t523 = t322 * (Ifges(4,1) * t258 - Ifges(4,4) * t356 - Ifges(4,5) * t427) + Ifges(3,1) * t428 + Ifges(3,5) * qJD(2) + t307;
t522 = m(7) * pkin(5) + t592;
t8 = -pkin(9) * t89 + t262 * t332 + t353;
t15 = t131 * t419 - t140 * t421 + t325 * t69 + t329 * t77;
t11 = t262 * qJ(5) + t301 * qJD(5) + t15;
t9 = pkin(9) * t90 + t11;
t1 = qJD(6) * t13 + t324 * t8 + t328 * t9;
t2 = -qJD(6) * t14 - t324 * t9 + t328 * t8;
t521 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t518 = -t16 * mrSges(5,1) + t12 * mrSges(6,1) + t15 * mrSges(5,2) - t11 * mrSges(6,3);
t474 = t301 / 0.2e1;
t490 = -t182 / 0.2e1;
t516 = Ifges(5,4) * t490 + Ifges(6,5) * t489 + t474 * t561 + t487 * t563 - t566;
t515 = -Ifges(5,2) * t490 + Ifges(6,3) * t489 + t474 * t560 + t487 * t562 + t573;
t513 = mrSges(5,1) * t161 + mrSges(6,1) * t25 + 0.2e1 * Ifges(6,3) * t499 - t89 * Ifges(5,4) / 0.2e1 - t262 * Ifges(5,6) / 0.2e1 + t589 * t501 + (t560 + Ifges(6,6)) * t480 + (-t500 + t499) * Ifges(5,2);
t511 = Ifges(7,4) * t509 + Ifges(7,2) * t508 + Ifges(7,6) * t481;
t510 = Ifges(7,1) * t509 + Ifges(7,4) * t508 + Ifges(7,5) * t481;
t452 = Ifges(7,4) * t574;
t44 = Ifges(7,2) * t105 + Ifges(7,6) * t293 + t452;
t505 = -t44 / 0.2e1;
t504 = t44 / 0.2e1;
t102 = Ifges(7,4) * t105;
t45 = Ifges(7,1) * t574 + Ifges(7,5) * t293 + t102;
t503 = -t45 / 0.2e1;
t502 = t45 / 0.2e1;
t495 = t105 / 0.2e1;
t493 = t574 / 0.2e1;
t492 = Ifges(4,1) * t484 + Ifges(4,4) * t485 + Ifges(4,5) * t478;
t488 = -t340 / 0.2e1;
t476 = t293 / 0.2e1;
t469 = pkin(3) * t321;
t468 = pkin(5) * t312;
t465 = g(3) * t326;
t314 = t326 * pkin(7);
t456 = Ifges(3,4) * t330;
t448 = qJ(5) * t182;
t442 = t235 * t326;
t438 = t321 * t326;
t434 = t462 * t326;
t433 = t326 * t331;
t303 = pkin(3) * t322 + pkin(2);
t284 = t330 * t303;
t431 = t330 * t331;
t255 = t322 * t274;
t189 = -pkin(8) * t436 + t255 + (-pkin(7) * t321 - pkin(3)) * t330;
t211 = pkin(7) * t435 + t321 * t274;
t197 = -pkin(8) * t438 + t211;
t118 = t325 * t189 + t329 * t197;
t236 = qJD(2) * t370 - t422;
t426 = qJD(2) * t326;
t413 = pkin(7) * t426;
t195 = t322 * t236 + t321 * t413;
t425 = qJD(2) * t330;
t310 = pkin(7) * t425;
t406 = t321 * t425;
t252 = pkin(3) * t406 + t310;
t267 = pkin(3) * t438 + t314;
t429 = t331 * pkin(1) + t327 * pkin(7);
t420 = qJD(4) * t326;
t416 = Ifges(7,5) * t23 + Ifges(7,6) * t24 + Ifges(7,3) * t250;
t410 = t462 * t433;
t317 = t331 * pkin(7);
t409 = -t327 * t434 + t331 * t469 + t317;
t39 = t90 * mrSges(5,1) + t89 * mrSges(5,2);
t38 = t90 * mrSges(6,1) - t89 * mrSges(6,3);
t66 = -t262 * mrSges(6,1) + t89 * mrSges(6,2);
t394 = -t417 / 0.2e1;
t393 = t417 / 0.2e1;
t144 = -t216 * mrSges(4,1) + t217 * mrSges(4,2);
t117 = t189 * t329 - t325 * t197;
t389 = t284 + t434;
t388 = t303 * t431 + t327 * t469 + t429;
t7 = -t24 * mrSges(7,1) + t23 * mrSges(7,2);
t111 = t118 - t447;
t234 = t263 * t326;
t384 = -qJ(5) * t234 - t267;
t112 = -t117 + t316;
t382 = mrSges(3,1) * t326 + mrSges(3,2) * t330;
t227 = -t327 * t312 + t330 * t430;
t228 = t311 * t327 + t312 * t431;
t142 = t227 * t328 - t228 * t324;
t143 = t227 * t324 + t228 * t328;
t377 = mrSges(7,1) * t142 - mrSges(7,2) * t143;
t376 = (mrSges(7,1) * t365 - mrSges(7,2) * t364) * t326;
t372 = Ifges(3,5) * t330 - Ifges(3,6) * t326;
t371 = Ifges(4,5) * t322 - Ifges(4,6) * t321;
t70 = pkin(5) * t330 + pkin(9) * t234 + t112;
t233 = t264 * t326;
t72 = pkin(9) * t233 + t111;
t32 = -t324 * t72 + t328 * t70;
t33 = t324 * t70 + t328 * t72;
t78 = -mrSges(7,2) * t293 + mrSges(7,3) * t105;
t79 = mrSges(7,1) * t293 - mrSges(7,3) * t574;
t368 = -t324 * t79 + t328 * t78;
t151 = t233 * t328 + t234 * t324;
t152 = t233 * t324 - t234 * t328;
t185 = t263 * t328 - t264 * t324;
t186 = t263 * t324 + t264 * t328;
t361 = qJ(5) * t264 + t303;
t153 = qJD(2) * t362 + t195;
t222 = t321 * t236;
t165 = qJD(2) * t355 + t222;
t48 = t153 * t329 - t325 * t165 - t189 * t421 - t197 * t419;
t360 = -pkin(4) * t312 - qJ(5) * t311 - t303;
t358 = pkin(1) * t382;
t357 = t326 * (Ifges(3,1) * t330 - t457);
t47 = t325 * t153 + t329 * t165 + t189 * t419 - t197 * t421;
t349 = t416 + t521;
t157 = -qJD(2) * t351 - t264 * t420;
t348 = qJ(5) * t157 - qJD(5) * t234 - t252;
t344 = t228 * pkin(4) + qJ(5) * t227 + t388;
t343 = -g(1) * t227 - g(2) * t225 - t311 * t465;
t40 = qJ(5) * t426 - qJD(5) * t330 + t47;
t341 = t330 * (Ifges(4,3) * t326 + t330 * t371);
t282 = t326 * t312 * qJ(5);
t280 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t427;
t213 = -mrSges(4,1) * t427 - mrSges(4,3) * t258;
t212 = mrSges(4,2) * t427 - mrSges(4,3) * t356;
t210 = -pkin(7) * t437 + t255;
t202 = -pkin(7) * t407 + t242;
t196 = -t322 * t413 + t222;
t178 = pkin(4) * t263 - t361;
t174 = Ifges(4,4) * t258 - Ifges(4,2) * t356 - Ifges(4,6) * t427;
t169 = mrSges(4,1) * t269 - mrSges(4,3) * t217;
t168 = -mrSges(4,2) * t269 + mrSges(4,3) * t216;
t158 = qJD(2) * t352 + t419 * t436 - t420 * t439;
t139 = t214 * t324 - t215 * t328;
t138 = t214 * t328 + t215 * t324;
t137 = t263 * t332 + t361;
t132 = pkin(4) * t233 - t384;
t125 = t217 * Ifges(4,4) + t216 * Ifges(4,2) + t269 * Ifges(4,6);
t116 = mrSges(5,1) * t182 + mrSges(5,2) * t340;
t115 = mrSges(6,1) * t182 - mrSges(6,3) * t340;
t114 = pkin(4) * t340 + t448;
t113 = t233 * t332 + t384;
t92 = -pkin(4) * t428 - t100;
t68 = -mrSges(6,2) * t90 + mrSges(6,3) * t262;
t67 = -mrSges(5,2) * t262 - mrSges(5,3) * t90;
t65 = mrSges(5,1) * t262 - mrSges(5,3) * t89;
t63 = t332 * t340 - t448;
t56 = pkin(4) * t158 - t348;
t54 = -qJD(6) * t152 - t157 * t324 + t158 * t328;
t53 = qJD(6) * t151 + t157 * t328 + t158 * t324;
t49 = -mrSges(7,1) * t105 + mrSges(7,2) * t574;
t42 = t158 * t332 + t348;
t41 = -pkin(4) * t426 - t48;
t31 = pkin(9) * t158 + t40;
t30 = -pkin(9) * t157 + qJD(2) * t411 - t48;
t20 = -mrSges(7,2) * t250 + mrSges(7,3) * t24;
t19 = mrSges(7,1) * t250 - mrSges(7,3) * t23;
t10 = pkin(5) * t90 + t25;
t4 = -qJD(6) * t33 + t30 * t328 - t31 * t324;
t3 = qJD(6) * t32 + t30 * t324 + t31 * t328;
t5 = [(t555 / 0.2e1 + mrSges(5,2) * t200 + mrSges(6,2) * t57 - mrSges(5,3) * t61 + t516) * t157 + (t456 * t393 - t561 * t501 + t559 * t480 - Ifges(4,3) * t478 + Ifges(7,3) * t481 - Ifges(6,6) * t499 - Ifges(5,6) * t500 - Ifges(4,5) * t484 - Ifges(4,6) * t485 + Ifges(7,6) * t508 + Ifges(7,5) * t509 + t121 * mrSges(4,2) - t120 * mrSges(4,1) + pkin(7) * (-qJDD(2) * mrSges(3,2) - mrSges(3,3) * t269) + t518 + t521) * t330 + (Ifges(3,4) * t270 - Ifges(3,2) * t269 + t416) * t330 / 0.2e1 + (Ifges(7,4) * t53 + Ifges(7,2) * t54) * t495 + (Ifges(7,4) * t152 + Ifges(7,2) * t151) * t508 + (t1 * t151 - t13 * t53 + t14 * t54 - t152 * t2) * mrSges(7,3) + (-m(5) * t409 + t366 * mrSges(7,1) + t533 * mrSges(7,2) + t569 * (-t226 * pkin(4) - qJ(5) * t225 + t409) + t570 * t317 + t522 * t226 + t564 * t225 + t528 * t331 + (m(3) * pkin(1) - m(4) * t396 + t545 + t529 * (-pkin(1) - t284) + (-m(7) * pkin(9) - t591) * t326 - t531 - t552) * t327) * g(1) - t571 * t234 + (-Ifges(7,5) * t493 + Ifges(5,6) * t490 + Ifges(6,6) * t489 - Ifges(7,6) * t495 - Ifges(7,3) * t476 - t559 * t474 + t572 + t588) * t426 + (-t120 * t436 - t121 * t438) * mrSges(4,3) + (-qJDD(2) * mrSges(3,1) + t144) * t314 + m(4) * (t120 * t210 + t121 * t211 + t190 * t195 + t191 * t196 + (-t386 * t425 + t442) * pkin(7)) + t383 * t446 + t536 * t565 + t374 * t568 + t380 * t442 + t53 * t502 + t54 * t504 - t125 * t438 / 0.2e1 + qJDD(2) * Ifges(3,6) * t330 + m(7) * (t1 * t33 - t10 * t113 + t13 * t4 + t14 * t3 + t2 * t32 + t42 * t55) + m(6) * (t11 * t111 + t112 * t12 + t132 * t25 + t40 * t58 + t41 * t57 + t56 * t71) + m(5) * (t117 * t16 + t118 * t15 + t161 * t267 + t200 * t252 + t47 * t62 + t48 * t61) + (Ifges(7,5) * t53 + Ifges(7,6) * t54) * t476 + (Ifges(7,5) * t152 + Ifges(7,6) * t151) * t481 + t436 * t492 - t174 * t406 / 0.2e1 + (Ifges(3,1) * t270 + Ifges(3,4) * t568 + Ifges(3,5) * qJDD(2) + t371 * t478 + t373 * t485 + t375 * t484 - t393 * t554) * t326 - t358 * t417 + t210 * t169 + t211 * t168 + t196 * t212 + t195 * t213 + (t372 * t565 - t526) * qJD(2) + (Ifges(7,1) * t53 + Ifges(7,4) * t54) * t493 + (Ifges(7,1) * t152 + Ifges(7,4) * t151) * t509 + t357 * t393 + t341 * t394 + t270 * t456 / 0.2e1 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t534) - (Ifges(4,5) * t217 + Ifges(4,6) * t216 + Ifges(4,3) * t269 + t535) * t330 / 0.2e1 - t539 * t310 + (t270 * t314 + t534) * mrSges(3,3) + t523 * t425 / 0.2e1 - t10 * (-mrSges(7,1) * t151 + mrSges(7,2) * t152) + t47 * t147 + t48 * t148 + t41 * t149 + t40 * t150 + Ifges(2,3) * qJDD(1) + t132 * t38 + t56 * t115 + t117 * t65 + t118 * t67 + t111 * t68 + t112 * t66 + t113 * t7 + (mrSges(5,1) * t200 - mrSges(6,2) * t58 - mrSges(5,3) * t62 + t515) * t158 + t3 * t78 + t4 * t79 + (-mrSges(6,2) * t11 - mrSges(5,3) * t15 + t513) * t233 + t55 * (-mrSges(7,1) * t54 + mrSges(7,2) * t53) + t42 * t49 + t33 * t20 + t32 * t19 + t152 * t510 + t151 * t511 + (-m(5) * (t388 + t410) - m(6) * (t344 + t410) - m(7) * t344 - t143 * mrSges(7,1) - t142 * mrSges(7,2) + t584 * t433 + t570 * t429 - t522 * t228 - t564 * t227 + t528 * t327 + (t537 + t552) * t331) * g(2) + t252 * t116 + t267 * t39 - pkin(1) * (mrSges(3,1) * t269 + mrSges(3,2) * t270) - t280 * t413; (t66 - t65) * t198 + (-t11 * t263 - t540 * t57 + t541 * t58) * mrSges(6,2) + (-t536 / 0.2e1 + t526 + (-t357 / 0.2e1 + t358 + t341 / 0.2e1) * qJD(1)) * qJD(1) + (Ifges(7,5) * t139 + Ifges(7,6) * t138) * t477 + (-Ifges(5,2) * t489 + Ifges(6,3) * t490 + t560 * t475 + t562 * t488 - t573) * t214 + t571 * t264 + (-t424 - t201) * t213 + (t15 * t199 - t16 * t198 - t161 * t303 - t200 * t251 + (-t101 + t145) * t62 + t586 * t61) * m(5) + (t168 * t322 - t169 * t321) * qJ(3) + (-m(5) * t389 - m(6) * (t389 + t538) - m(7) * (t284 + t538) - t383 + (-m(7) * t468 - t579) * t330 + t531 + t537) * g(3) + (t11 * t199 + t12 * t198 + t178 * t25 + t544 * t71 + t590 * t58 + (t146 - t92) * t57) * m(6) + (t423 - t202) * t212 + (-pkin(2) * t235 + (-t190 * t321 + t191 * t322) * qJD(3) + t532 * qJ(3) - t190 * t201 - t191 * t202 + t386 * t309) * m(4) + t551 * t49 + (Ifges(4,5) * t321 + t322 * Ifges(4,6)) * t478 + (Ifges(7,5) * t186 + Ifges(7,6) * t185) * t481 + (-t15 * t263 + t540 * t61 + t541 * t62) * mrSges(5,3) + t139 * t503 + t138 * t505 + t280 * t308 + (t67 + t68) * t199 + (Ifges(7,1) * t139 + Ifges(7,4) * t138) * t494 + (Ifges(4,1) * t321 + t454) * t484 + (Ifges(4,2) * t322 + t455) * t485 + t321 * t492 + (t1 * t185 + t13 * t139 - t138 * t14 - t186 * t2) * mrSges(7,3) + t235 * t381 + (-mrSges(5,1) * t541 - mrSges(5,2) * t540) * t200 - t542 * t146 + t543 * t145 + t544 * t115 + t557 * t78 + t558 * t79 + (t1 * t74 - t10 * t137 + t13 * t558 + t14 * t557 + t2 * t73 + t55 * t551) * m(7) + t372 * t394 + (-Ifges(3,2) * t400 + Ifges(5,6) * t489 + Ifges(6,6) * t490 + t488 * t561 - t572 - t585) * t428 + t539 * t309 + t532 * mrSges(4,3) - t10 * (-mrSges(7,1) * t185 + mrSges(7,2) * t186) + t178 * t38 + (Ifges(7,1) * t493 + Ifges(7,4) * t495 + Ifges(7,5) * t476 + t502 + t524) * (qJD(6) * t185 - t240 * t328 + t241 * t324) + (Ifges(7,4) * t493 + Ifges(7,2) * t495 + Ifges(7,6) * t476 + t504 + t525) * (-qJD(6) * t186 + t240 * t324 + t241 * t328) + (t307 + t523) * t400 - pkin(2) * t144 - t101 * t147 - t100 * t148 - t92 * t149 - t91 * t150 - t55 * (-mrSges(7,1) * t138 + mrSges(7,2) * t139) + t137 * t7 + t515 * t241 - t516 * t240 + t513 * t263 + t73 * t19 + t74 * t20 + Ifges(3,3) * qJDD(2) + t240 * t577 + (Ifges(7,4) * t186 + Ifges(7,2) * t185) * t508 + (Ifges(7,1) * t186 + Ifges(7,4) * t185) * t509 + t186 * t510 + t185 * t511 + t174 * t408 / 0.2e1 + t575 * (t382 + ((-m(5) - m(6)) * t462 + t580 + t584) * t330 + (t350 + m(5) * t303 - m(6) * t360 - m(7) * (t360 - t468) + t579) * t326) - (Ifges(5,4) * t489 + Ifges(6,5) * t490 + t475 * t561 + t488 * t563 + t566 + t577) * t215 - t251 * t116 - t256 * mrSges(3,2) - t257 * mrSges(3,1) - Ifges(3,6) * t269 + Ifges(3,5) * t270 - t303 * t39 + t322 * t125 / 0.2e1 + (Ifges(7,4) * t139 + Ifges(7,2) * t138) * t496; t38 - t7 + t543 * t182 + t542 * t340 + t144 + t105 * t78 - t574 * t79 + t258 * t213 + t39 + t356 * t212 + (t330 * g(3) - t326 * t575) * (m(4) - t529) + (t105 * t14 - t13 * t574 + t10) * m(7) + (t182 * t58 - t340 * t57 + t25) * m(6) + (t182 * t62 + t340 * t61 + t161) * m(5) + (t190 * t258 + t191 * t356 + t235) * m(4); (t569 * (-t225 * pkin(4) + qJ(5) * t226) - t564 * t226 + t522 * t225 + t578) * g(2) + (Ifges(6,3) * t340 - t451) * t490 - t200 * (mrSges(5,1) * t340 - mrSges(5,2) * t182) - t71 * (mrSges(6,1) * t340 + mrSges(6,3) * t182) + t535 - (-mrSges(5,1) * t311 - mrSges(5,2) * t312) * t465 + (t376 - m(6) * t282 - (t312 * mrSges(6,3) + (-m(6) * pkin(4) - mrSges(6,1)) * t311) * t326) * g(3) + t549 * t78 + (-g(3) * (t311 * t411 + t282) - t55 * t63 + t1 * t272 + t2 * t271 + t549 * t14 + t550 * t13) * m(7) + t550 * t79 + t58 * t460 + t57 * t461 + t97 * t487 + (-pkin(4) * t12 + qJ(5) * t11 + qJD(5) * t58 - t114 * t71) * m(6) + (t377 + t569 * (-t227 * pkin(4) + qJ(5) * t228) - t564 * t228 + t522 * t227) * g(1) - t518 - t349 + (-m(6) * t57 + t458 + t542) * t62 + (-m(6) * t58 - t459 - t543) * t61 + (-t182 * t561 + t340 * t560) * t475 + (-t182 * t563 + t179 - t453 + t94) * t488 + (-Ifges(5,2) * t340 - t180 + t555) * t489 + (Ifges(7,4) * t494 + Ifges(7,2) * t496 + Ifges(7,6) * t477 + t505 - t525) * t574 - (Ifges(7,1) * t494 + Ifges(7,4) * t496 + Ifges(7,5) * t477 + t503 - t524) * t105 + qJD(5) * t150 - t114 * t115 - t63 * t49 - pkin(4) * t66 + qJ(5) * t68 + t271 * t19 + t272 * t20; t328 * t19 + t324 * t20 + (t115 - t49) * t340 + t368 * qJD(6) + (-t150 - t368) * t301 + t66 + (t1 * t324 - t340 * t55 + t2 * t328 + t343 + t293 * (-t13 * t324 + t14 * t328)) * m(7) + (-t301 * t58 + t340 * t71 + t12 + t343) * m(6); -t55 * (mrSges(7,1) * t574 + mrSges(7,2) * t105) + (Ifges(7,1) * t105 - t452) * t494 + t44 * t493 + (Ifges(7,5) * t105 - Ifges(7,6) * t574) * t477 - t13 * t78 + t14 * t79 - g(1) * t377 - g(2) * t578 - g(3) * t376 + (t105 * t13 + t14 * t574) * mrSges(7,3) + t349 + (-Ifges(7,2) * t574 + t102 + t45) * t496;];
tau  = t5;
