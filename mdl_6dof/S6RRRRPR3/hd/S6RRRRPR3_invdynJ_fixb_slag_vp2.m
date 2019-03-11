% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:01:57
% EndTime: 2019-03-09 22:02:37
% DurationCPUTime: 21.88s
% Computational Cost: add. (18019->784), mult. (41136->969), div. (0->0), fcn. (30185->14), ass. (0->390)
t589 = mrSges(5,1) - mrSges(6,2);
t324 = sin(qJ(6));
t329 = cos(qJ(6));
t581 = mrSges(7,1) * t324 + mrSges(7,2) * t329;
t326 = sin(qJ(3));
t327 = sin(qJ(2));
t330 = cos(qJ(2));
t506 = cos(qJ(3));
t367 = t326 * t327 - t330 * t506;
t240 = t367 * qJD(1);
t259 = t326 * t330 + t327 * t506;
t241 = t259 * qJD(1);
t325 = sin(qJ(4));
t505 = cos(qJ(4));
t366 = -t325 * t240 + t241 * t505;
t189 = t505 * t240 + t241 * t325;
t459 = qJ(5) * t189;
t136 = pkin(4) * t366 + t459;
t424 = qJD(3) + qJD(4);
t315 = qJD(2) + t424;
t510 = -t315 / 0.2e1;
t515 = t366 / 0.2e1;
t516 = -t366 / 0.2e1;
t187 = qJD(6) + t366;
t520 = -t187 / 0.2e1;
t167 = t189 * t324 + t315 * t329;
t522 = -t167 / 0.2e1;
t166 = t189 * t329 - t315 * t324;
t523 = -t166 / 0.2e1;
t319 = t330 * pkin(2);
t308 = t319 + pkin(1);
t285 = t308 * qJD(1);
t212 = pkin(3) * t240 - t285;
t346 = -qJ(5) * t366 + t212;
t107 = pkin(4) * t189 + t346;
t526 = pkin(4) + pkin(10);
t332 = -pkin(8) - pkin(7);
t287 = t332 * t330;
t265 = qJD(1) * t287;
t242 = t326 * t265;
t286 = t332 * t327;
t264 = qJD(1) * t286;
t249 = qJD(2) * pkin(2) + t264;
t198 = t506 * t249 + t242;
t234 = t241 * pkin(9);
t164 = t198 - t234;
t321 = qJD(2) + qJD(3);
t352 = pkin(3) * t321 + t164;
t146 = t505 * t352;
t245 = t506 * t265;
t199 = t326 * t249 - t245;
t488 = t240 * pkin(9);
t165 = t199 - t488;
t444 = t325 * t165;
t103 = -t146 + t444;
t490 = t366 * pkin(5);
t369 = t103 + t490;
t578 = qJD(5) + t369;
t52 = -t315 * t526 + t578;
t81 = t189 * t526 + t346;
t27 = -t324 * t81 + t329 * t52;
t28 = t324 * t52 + t329 * t81;
t533 = t27 * mrSges(7,1) + t212 * mrSges(5,2) - t28 * mrSges(7,2) - t107 * mrSges(6,3);
t568 = Ifges(6,4) - Ifges(5,5);
t588 = -Ifges(5,1) * t516 - Ifges(7,5) * t522 + Ifges(6,2) * t515 - Ifges(7,6) * t523 - Ifges(7,3) * t520 + t510 * t568 + t533;
t374 = Ifges(7,5) * t324 + Ifges(7,6) * t329;
t468 = Ifges(7,4) * t324;
t376 = Ifges(7,2) * t329 + t468;
t467 = Ifges(7,4) * t329;
t378 = Ifges(7,1) * t324 + t467;
t461 = t329 * mrSges(7,3);
t475 = mrSges(7,3) * t324;
t508 = -t324 / 0.2e1;
t518 = t189 / 0.2e1;
t519 = -t189 / 0.2e1;
t154 = t505 * t165;
t347 = t325 * t352;
t104 = t154 + t347;
t537 = t212 * mrSges(5,1) - t107 * mrSges(6,2) - t104 * mrSges(5,3);
t567 = Ifges(6,5) - Ifges(5,6);
t380 = mrSges(7,1) * t329 - mrSges(7,2) * t324;
t586 = t189 * pkin(5);
t95 = -t315 * qJ(5) - t104;
t55 = -t95 - t586;
t469 = Ifges(7,4) * t167;
t75 = Ifges(7,2) * t166 + Ifges(7,6) * t187 + t469;
t579 = t55 * t380 - t329 * t75 / 0.2e1;
t158 = Ifges(7,4) * t166;
t76 = t167 * Ifges(7,1) + t187 * Ifges(7,5) + t158;
t587 = -Ifges(5,2) * t518 + Ifges(6,3) * t519 + t27 * t475 - t28 * t461 + t374 * t520 + t376 * t523 + t378 * t522 + t76 * t508 + t510 * t567 - t537 + t579;
t183 = Ifges(5,4) * t189;
t564 = Ifges(7,3) * t187;
t565 = Ifges(7,6) * t166;
t561 = Ifges(5,1) * t366 + Ifges(5,5) * t315 + Ifges(7,5) * t167 - t183 + t564 + t565;
t207 = -t264 * t326 + t245;
t171 = t207 + t488;
t208 = t506 * t264 + t242;
t172 = -t234 + t208;
t110 = t325 * t171 + t172 * t505;
t501 = pkin(2) * t326;
t299 = t325 * t501;
t420 = t506 * pkin(2);
t307 = t420 + pkin(3);
t404 = qJD(3) * t506;
t392 = pkin(2) * t404;
t403 = qJD(4) * t505;
t194 = -t299 * t424 + t307 * t403 + t505 * t392;
t193 = -qJD(5) - t194;
t585 = t193 + t110;
t106 = t164 * t505 - t444;
t391 = pkin(3) * t403;
t289 = t391 + qJD(5);
t584 = -t289 + t106;
t323 = qJ(2) + qJ(3);
t316 = sin(t323);
t317 = cos(t323);
t318 = qJ(4) + t323;
t304 = sin(t318);
t305 = cos(t318);
t381 = mrSges(5,1) * t304 + mrSges(5,2) * t305;
t583 = mrSges(4,1) * t316 + mrSges(4,2) * t317 + t381;
t426 = qJD(1) * qJD(2);
t268 = qJDD(1) * t330 - t327 * t426;
t582 = t581 * t305;
t328 = sin(qJ(1));
t331 = cos(qJ(1));
t547 = g(1) * t331 + g(2) * t328;
t580 = -t589 * t305 + (mrSges(5,2) - mrSges(6,3)) * t304;
t182 = Ifges(6,6) * t189;
t128 = Ifges(6,4) * t315 - Ifges(6,2) * t366 + t182;
t94 = -pkin(4) * t315 + qJD(5) + t103;
t577 = -mrSges(6,1) * t94 - mrSges(5,3) * t103 + t128 / 0.2e1;
t458 = qJDD(1) * pkin(1);
t232 = -t268 * pkin(2) - t458;
t269 = qJDD(1) * t327 + t330 * t426;
t353 = t259 * qJD(3);
t343 = -qJD(1) * t353 + t268 * t506 - t326 * t269;
t141 = -pkin(3) * t343 + t232;
t354 = t367 * qJD(3);
t170 = -qJD(1) * t354 + t326 * t268 + t269 * t506;
t429 = qJD(4) * t325;
t83 = -t505 * t170 + t240 * t403 + t241 * t429 - t325 * t343;
t84 = qJD(4) * t366 + t325 * t170 - t343 * t505;
t20 = t84 * pkin(4) + t83 * qJ(5) - qJD(5) * t366 + t141;
t12 = t84 * pkin(10) + t20;
t320 = qJDD(2) + qJDD(3);
t314 = qJDD(4) + t320;
t254 = t269 * pkin(7);
t210 = qJDD(2) * pkin(2) - pkin(8) * t269 - t254;
t253 = t268 * pkin(7);
t211 = pkin(8) * t268 + t253;
t116 = -qJD(3) * t199 + t506 * t210 - t326 * t211;
t58 = t320 * pkin(3) - t170 * pkin(9) + t116;
t430 = qJD(3) * t326;
t115 = t326 * t210 + t506 * t211 + t249 * t404 + t265 * t430;
t70 = pkin(9) * t343 + t115;
t18 = -qJD(4) * t347 - t165 * t403 - t325 * t70 + t505 * t58;
t358 = qJDD(5) - t18;
t6 = -t83 * pkin(5) - t314 * t526 + t358;
t1 = qJD(6) * t27 + t12 * t329 + t324 * t6;
t44 = qJD(6) * t166 + t314 * t329 + t324 * t84;
t45 = -qJD(6) * t167 - t314 * t324 + t329 * t84;
t82 = qJDD(6) - t83;
t10 = t44 * Ifges(7,4) + t45 * Ifges(7,2) + t82 * Ifges(7,6);
t466 = Ifges(6,6) * t366;
t127 = t315 * Ifges(6,5) + t189 * Ifges(6,3) - t466;
t470 = Ifges(5,4) * t366;
t129 = -t189 * Ifges(5,2) + t315 * Ifges(5,6) + t470;
t17 = qJD(4) * t146 - t165 * t429 + t325 * t58 + t505 * t70;
t13 = -qJ(5) * t314 - qJD(5) * t315 - t17;
t15 = -t314 * pkin(4) + t358;
t2 = -qJD(6) * t28 - t12 * t324 + t329 * t6;
t428 = qJD(6) * t324;
t400 = -t428 / 0.2e1;
t427 = qJD(6) * t329;
t527 = t82 / 0.2e1;
t528 = t45 / 0.2e1;
t529 = t44 / 0.2e1;
t530 = Ifges(7,1) * t529 + Ifges(7,4) * t528 + Ifges(7,5) * t527;
t7 = -pkin(5) * t84 - t13;
t576 = -t2 * t461 + t7 * t581 - t1 * t475 + t18 * mrSges(5,1) - t13 * mrSges(6,3) + t15 * mrSges(6,2) - t17 * mrSges(5,2) + (Ifges(7,5) * t329 - Ifges(7,6) * t324) * t527 + (-Ifges(7,2) * t324 + t467) * t528 + (Ifges(7,1) * t329 - t468) * t529 + t329 * t530 + t76 * t400 + t10 * t508 + t567 * t84 + t568 * t83 + (Ifges(6,1) + Ifges(5,3)) * t314 + t579 * qJD(6) + (t27 * t428 - t28 * t427) * mrSges(7,3) - (t166 * t376 + t167 * t378 + t187 * t374) * qJD(6) / 0.2e1 + (t466 + t129) * t515 + (-t470 + t127) * t516;
t573 = -t170 / 0.2e1;
t572 = t268 / 0.2e1;
t511 = t314 / 0.2e1;
t507 = t330 / 0.2e1;
t571 = t343 / 0.2e1;
t570 = m(5) * t103;
t19 = -mrSges(7,1) * t45 + mrSges(7,2) * t44;
t64 = mrSges(6,1) * t84 - mrSges(6,3) * t314;
t566 = t19 - t64;
t563 = t330 * Ifges(3,2);
t560 = t490 - t585;
t559 = t490 - t584;
t108 = -mrSges(7,1) * t166 + mrSges(7,2) * t167;
t484 = mrSges(6,1) * t189;
t173 = -mrSges(6,3) * t315 + t484;
t558 = t108 - t173;
t476 = mrSges(5,3) * t189;
t175 = -mrSges(5,2) * t315 - t476;
t557 = t173 - t175;
t483 = mrSges(6,1) * t366;
t556 = -mrSges(5,3) * t366 + t315 * t589 - t483;
t214 = t326 * t286 - t506 * t287;
t448 = t305 * t331;
t450 = t304 * t331;
t555 = pkin(4) * t448 + qJ(5) * t450;
t398 = t317 * mrSges(4,1) - mrSges(4,2) * t316;
t432 = qJD(1) * t330;
t433 = qJD(1) * t327;
t494 = pkin(7) * t330;
t495 = pkin(7) * t327;
t552 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t433) * t494 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t432) * t495;
t274 = qJ(5) * t448;
t436 = mrSges(6,2) * t450 + mrSges(6,3) * t448;
t551 = -m(7) * t274 - t436;
t449 = t305 * t328;
t272 = qJ(5) * t449;
t451 = t304 * t328;
t437 = mrSges(6,2) * t451 + mrSges(6,3) * t449;
t550 = -m(7) * t272 - t437;
t549 = t253 * t330 + t254 * t327;
t24 = mrSges(7,1) * t82 - mrSges(7,3) * t44;
t25 = -mrSges(7,2) * t82 + mrSges(7,3) * t45;
t548 = t329 * t24 + t324 * t25;
t546 = 0.2e1 * t511;
t544 = -t305 * mrSges(7,3) - t304 * t581 + t580;
t491 = g(3) * t305;
t543 = -t304 * t547 + t491;
t542 = -t398 + t544;
t541 = -m(6) * t95 - t557;
t540 = m(6) * t94 - t556;
t539 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t284 = -mrSges(3,1) * t330 + mrSges(3,2) * t327;
t538 = -m(3) * pkin(1) - m(4) * t308 - mrSges(2,1) + t284 - t398 + t580;
t118 = -mrSges(7,2) * t187 + mrSges(7,3) * t166;
t119 = mrSges(7,1) * t187 - mrSges(7,3) * t167;
t372 = t27 * t324 - t28 * t329;
t345 = -qJD(6) * t372 + t1 * t324 + t2 * t329;
t536 = m(7) * t345 + t118 * t427 - t119 * t428 + t548;
t322 = -pkin(9) + t332;
t535 = -m(3) * pkin(7) + m(4) * t332 - m(7) * (pkin(5) - t322) - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t534 = m(7) * (t27 * t329 + t28 * t324) + t329 * t119 + t324 * t118 + t540;
t359 = Ifges(4,4) * t367;
t471 = Ifges(4,4) * t259;
t512 = t241 / 0.2e1;
t532 = -t285 * (mrSges(4,1) * t259 - mrSges(4,2) * t367) + t321 * (-Ifges(4,5) * t367 - Ifges(4,6) * t259) / 0.2e1 - t240 * (-Ifges(4,2) * t259 - t359) / 0.2e1 + (-Ifges(4,1) * t367 - t471) * t512;
t521 = t167 / 0.2e1;
t509 = t315 / 0.2e1;
t500 = pkin(2) * t327;
t499 = pkin(3) * t241;
t498 = pkin(3) * t316;
t303 = pkin(3) * t317;
t497 = pkin(3) * t325;
t295 = t305 * pkin(4);
t487 = t95 * mrSges(6,1);
t479 = mrSges(4,3) * t240;
t478 = mrSges(4,3) * t241;
t474 = Ifges(3,4) * t327;
t473 = Ifges(3,4) * t330;
t472 = Ifges(4,4) * t241;
t337 = -qJD(2) * t259 - t353;
t338 = -qJD(2) * t367 - t354;
t357 = t325 * t367;
t113 = -qJD(4) * t357 + t259 * t403 + t325 * t338 - t337 * t505;
t457 = t113 * t324;
t456 = t113 * t329;
t351 = t505 * t367;
t205 = t259 * t325 + t351;
t455 = t205 * t324;
t454 = t205 * t329;
t290 = t304 * qJ(5);
t446 = t324 * t328;
t445 = t324 * t331;
t443 = t328 * t329;
t441 = t329 * t331;
t440 = t582 * t328;
t439 = t582 * t331;
t435 = t295 + t290;
t434 = t303 + t319;
t431 = qJD(2) * t327;
t423 = Ifges(7,5) * t44 + Ifges(7,6) * t45 + Ifges(7,3) * t82;
t419 = t505 * pkin(3);
t312 = pkin(2) * t431;
t311 = pkin(2) * t433;
t411 = t303 + t435;
t410 = qJD(2) * t332;
t409 = t505 * t326;
t65 = -t83 * mrSges(6,1) + t314 * mrSges(6,2);
t399 = t426 / 0.2e1;
t263 = pkin(1) + t434;
t395 = -t263 - t290;
t105 = t164 * t325 + t154;
t109 = -t505 * t171 + t172 * t325;
t213 = t506 * t286 + t287 * t326;
t184 = -pkin(9) * t259 + t213;
t185 = -pkin(9) * t367 + t214;
t131 = -t505 * t184 + t185 * t325;
t266 = t327 * t410;
t267 = t330 * t410;
t394 = -t326 * t266 + t506 * t267;
t248 = t331 * t263;
t393 = -t322 * t328 + t248;
t390 = t319 + t411;
t389 = -m(7) * t526 - mrSges(7,3);
t306 = -t419 - pkin(4);
t388 = -pkin(4) * t451 + t272;
t387 = -pkin(4) * t450 + t274;
t386 = -pkin(4) * t304 - t498;
t235 = t307 * t505 - t299;
t384 = mrSges(3,1) * t327 + mrSges(3,2) * t330;
t377 = t474 + t563;
t375 = Ifges(3,5) * t330 - Ifges(3,6) * t327;
t206 = t259 * t505 - t357;
t219 = pkin(3) * t367 - t308;
t344 = -t206 * qJ(5) + t219;
t93 = t205 * t526 + t344;
t96 = pkin(5) * t206 + t131;
t39 = t324 * t96 + t329 * t93;
t38 = -t324 * t93 + t329 * t96;
t371 = t118 * t329 - t119 * t324;
t230 = -pkin(4) - t235;
t370 = t389 * t304;
t368 = pkin(1) * t384;
t132 = t325 * t184 + t185 * t505;
t365 = t205 * t427 + t457;
t364 = t205 * t428 - t456;
t363 = t327 * (Ifges(3,1) * t330 - t474);
t147 = t506 * t266 + t326 * t267 + t286 * t404 + t287 * t430;
t123 = pkin(9) * t337 + t147;
t336 = -pkin(9) * t338 - t286 * t430 + t287 * t404 + t394;
t36 = -t505 * t123 - t184 * t403 + t185 * t429 - t325 * t336;
t236 = pkin(2) * t409 + t325 * t307;
t117 = t499 + t136;
t111 = t117 + t311;
t350 = -m(7) * t498 + t370;
t349 = t328 * t370 + t440;
t348 = t331 * t370 + t439;
t37 = qJD(4) * t132 + t325 * t123 - t505 * t336;
t192 = -pkin(3) * t337 + t312;
t112 = qJD(4) * t351 + t259 * t429 - t325 * t337 - t338 * t505;
t35 = t113 * pkin(4) + t112 * qJ(5) - t206 * qJD(5) + t192;
t180 = -Ifges(4,2) * t240 + Ifges(4,6) * t321 + t472;
t233 = Ifges(4,4) * t240;
t181 = Ifges(4,1) * t241 + Ifges(4,5) * t321 - t233;
t334 = t285 * (mrSges(4,1) * t241 - mrSges(4,2) * t240) - t241 * (-Ifges(4,1) * t240 - t472) / 0.2e1 - t321 * (-Ifges(4,5) * t240 - Ifges(4,6) * t241) / 0.2e1 + t199 * t478 - t198 * t479 + Ifges(4,6) * t343 + t180 * t512 + Ifges(4,3) * t320 - t115 * mrSges(4,2) + t116 * mrSges(4,1) + Ifges(4,5) * t170 + t561 * t518 + (-Ifges(4,2) * t241 + t181 - t233) * t240 / 0.2e1 + (-Ifges(5,4) * t518 + Ifges(6,6) * t519 - t577 + t588) * t189 + (-t487 + t587) * t366 + t576;
t310 = Ifges(3,4) * t432;
t302 = qJ(5) + t497;
t294 = t305 * pkin(10);
t270 = -t498 - t500;
t251 = t331 * t270;
t250 = t328 * t270;
t239 = Ifges(3,1) * t433 + Ifges(3,5) * qJD(2) + t310;
t238 = Ifges(3,6) * qJD(2) + qJD(1) * t377;
t229 = qJ(5) + t236;
t227 = -t304 * t446 + t441;
t226 = t304 * t443 + t445;
t225 = t304 * t445 + t443;
t224 = t304 * t441 - t446;
t217 = mrSges(4,1) * t321 - t478;
t216 = -mrSges(4,2) * t321 - t479;
t215 = t311 + t499;
t197 = mrSges(4,1) * t240 + mrSges(4,2) * t241;
t186 = t366 * pkin(10);
t152 = -t320 * mrSges(4,2) + mrSges(4,3) * t343;
t151 = mrSges(4,1) * t320 - mrSges(4,3) * t170;
t148 = -qJD(3) * t214 + t394;
t140 = -mrSges(6,2) * t189 - mrSges(6,3) * t366;
t139 = mrSges(5,1) * t189 + mrSges(5,2) * t366;
t124 = t205 * pkin(4) + t344;
t102 = t366 * t526 + t459;
t97 = -t205 * pkin(5) + t132;
t88 = t117 + t186;
t87 = t111 + t186;
t85 = t109 - t586;
t71 = t105 - t586;
t63 = -mrSges(5,2) * t314 - mrSges(5,3) * t84;
t62 = mrSges(5,1) * t314 + mrSges(5,3) * t83;
t61 = t104 - t586;
t34 = t102 * t329 + t324 * t61;
t33 = -t102 * t324 + t329 * t61;
t32 = t324 * t85 + t329 * t87;
t31 = -t324 * t87 + t329 * t85;
t30 = t324 * t71 + t329 * t88;
t29 = -t324 * t88 + t329 * t71;
t23 = -t112 * pkin(5) + t37;
t22 = -pkin(5) * t113 - t36;
t21 = t113 * pkin(10) + t35;
t4 = -qJD(6) * t39 - t21 * t324 + t23 * t329;
t3 = qJD(6) * t38 + t21 * t329 + t23 * t324;
t5 = [t187 * (Ifges(7,5) * t365 - Ifges(7,6) * t364) / 0.2e1 + t166 * (Ifges(7,4) * t365 - Ifges(7,2) * t364) / 0.2e1 + (-t198 * t338 + t199 * t337) * mrSges(4,3) + (t15 * mrSges(6,1) + t141 * mrSges(5,2) - t18 * mrSges(5,3) - t20 * mrSges(6,3) + Ifges(5,5) * t511 + Ifges(7,5) * t529 + Ifges(7,6) * t528 + Ifges(7,3) * t527 + (-Ifges(6,6) - Ifges(5,4) / 0.2e1) * t84 - Ifges(6,2) * t83 - t546 * Ifges(6,4) + t539) * t206 + (-Ifges(5,4) * t84 + Ifges(5,5) * t314 + t423) * t206 / 0.2e1 + (-Ifges(5,4) * t519 - Ifges(7,5) * t521 - Ifges(5,1) * t515 + Ifges(6,2) * t516 + Ifges(6,6) * t518 - t565 / 0.2e1 - t564 / 0.2e1 + t568 * t509 - t533 + t577 - t561 / 0.2e1) * t112 + (t330 * t473 + t363) * t399 + (t239 * t507 + t375 * qJD(2) / 0.2e1 + t532 - t552) * qJD(2) + (-m(7) * (pkin(10) * t448 + t248 + t555) - t225 * mrSges(7,1) - t224 * mrSges(7,2) - mrSges(7,3) * t448 - m(6) * (t393 + t555) - m(5) * t393 + t538 * t331 + t535 * t328) * g(2) + (t268 * t494 + t269 * t495 + t549) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t549) + (-m(5) * t104 - t541) * t36 + (-t227 * mrSges(7,1) + t226 * mrSges(7,2) + ((m(5) + m(6)) * t322 + t535) * t331 + (-m(6) * (t395 - t295) - m(7) * t395 - t305 * t389 + m(5) * t263 - t538) * t328) * g(1) + t532 * qJD(3) + (Ifges(7,1) * t365 - Ifges(7,4) * t364) * t521 + (t1 * t454 - t2 * t455 - t27 * t365 - t28 * t364) * mrSges(7,3) + m(4) * (t115 * t214 + t116 * t213 + t147 * t199 + t148 * t198 - t232 * t308 - t285 * t312) + t337 * t180 / 0.2e1 + t338 * t181 / 0.2e1 + (t141 * mrSges(5,1) + t13 * mrSges(6,1) - t20 * mrSges(6,2) - t17 * mrSges(5,3) + t374 * t527 + t376 * t528 + t378 * t529 - t7 * t380 + t75 * t400 + (Ifges(6,3) + Ifges(5,2)) * t84 + (Ifges(6,6) + Ifges(5,4)) * t83 + t567 * t546) * t205 - t308 * (-mrSges(4,1) * t343 + t170 * mrSges(4,2)) + m(6) * (t107 * t35 + t124 * t20) + m(5) * (t141 * t219 + t192 * t212) + (qJD(6) * t76 + t10) * t454 / 0.2e1 + (Ifges(3,4) * t269 + Ifges(3,2) * t268) * t507 - t284 * t458 + t75 * t456 / 0.2e1 + t76 * t457 / 0.2e1 + m(7) * (t1 * t39 + t2 * t38 + t22 * t55 + t27 * t4 + t28 * t3 + t7 * t97) - t238 * t431 / 0.2e1 - t368 * t426 + (m(5) * t17 - m(6) * t13 + t63 - t64) * t132 + t471 * t571 + t377 * t572 + t359 * t573 + (-m(5) * t18 + m(6) * t15 - t62 + t65) * t131 + Ifges(2,3) * qJDD(1) + t455 * t530 + t269 * t473 / 0.2e1 - pkin(1) * (-mrSges(3,1) * t268 + mrSges(3,2) * t269) + t219 * (mrSges(5,1) * t84 - mrSges(5,2) * t83) + t213 * t151 + t214 * t152 + t147 * t216 + t148 * t217 + t192 * t139 - t83 * Ifges(5,1) * t206 + t38 * t24 + t39 * t25 + t197 * t312 + t55 * (mrSges(7,1) * t364 + mrSges(7,2) * t365) + (-Ifges(5,2) * t519 - Ifges(5,4) * t515 + Ifges(6,6) * t516 + Ifges(6,3) * t518 + t487 + t127 / 0.2e1 - t129 / 0.2e1 + t567 * t509 + t537) * t113 + t97 * t19 + (t540 + t570) * t37 + (t232 * mrSges(4,2) - t116 * mrSges(4,3) + Ifges(4,1) * t170 + Ifges(4,4) * t571 + Ifges(4,5) * t320) * t259 + (-mrSges(3,1) * t495 - mrSges(3,2) * t494 + 0.2e1 * Ifges(3,6) * t507) * qJDD(2) + (Ifges(3,1) * t269 + Ifges(3,4) * t572 + Ifges(3,5) * qJDD(2) - t399 * t563) * t327 + (t232 * mrSges(4,1) - t115 * mrSges(4,3) + Ifges(4,4) * t573 - Ifges(4,2) * t343 - Ifges(4,6) * t320) * t367 + t22 * t108 + t3 * t118 + t4 * t119 + t124 * (-mrSges(6,2) * t84 + mrSges(6,3) * t83) + t35 * t140; (-m(4) * t319 - m(5) * t434 - m(7) * (t294 + t390) - m(6) * t390 + t284 + t542) * g(3) + t536 * (-pkin(10) + t230) + (-t103 * t109 + t17 * t236 + t18 * t235 - t212 * t215 + (t194 - t110) * t104) * m(5) + (-t436 - t348) * g(1) + (t552 + (t368 - t363 / 0.2e1) * qJD(1)) * qJD(1) + (-t13 * t229 + t15 * t230 - t107 * t111 - t109 * t94 - g(1) * (t251 + t387) - g(2) * (t250 + t388) + t585 * t95) * m(6) + ((t506 * t116 + t115 * t326 + (-t198 * t326 + t199 * t506) * qJD(3)) * pkin(2) - t198 * t207 - t199 * t208 + t285 * t311) * m(4) + t547 * (m(4) * t500 - m(5) * t270 + t384 + t583) + (-pkin(2) * t430 - t207) * t217 + (t392 - t208) * t216 + (-t437 - t349) * g(2) - (-Ifges(3,2) * t433 + t239 + t310) * t432 / 0.2e1 + t334 + t238 * t433 / 0.2e1 - t375 * t426 / 0.2e1 - t197 * t311 + t151 * t420 + Ifges(3,3) * qJDD(2) + Ifges(3,6) * t268 + Ifges(3,5) * t269 - t253 * mrSges(3,2) - t254 * mrSges(3,1) + t230 * t65 + t235 * t62 + t236 * t63 - t215 * t139 + t194 * t175 + t193 * t173 + t556 * t109 + t557 * t110 + (t229 * t7 - t27 * t31 - t28 * t32 - g(1) * (t251 + t274) - g(2) * (t250 + t272) + t560 * t55) * m(7) + t560 * t108 + t566 * t229 + t152 * t501 + (t534 + t570) * (t307 * t429 + (t326 * t403 + (t325 * t506 + t409) * qJD(3)) * pkin(2)) - t32 * t118 - t31 * t119 - t111 * t140; t534 * pkin(3) * t429 + (-t331 * t350 - t439 + t551) * g(1) + (-t328 * t350 - t440 + t550) * g(2) + (-m(5) * t303 - m(6) * t411 - m(7) * (t294 + t411) + t542) * g(3) + t536 * (-pkin(10) + t306) + (m(5) * t498 + t583) * t547 + (-t105 * t94 - t107 * t117 - g(1) * (t331 * t386 + t274) - g(2) * (t328 * t386 + t272) - t13 * t302 + t15 * t306 + t584 * t95) * m(6) - t139 * t499 + t334 + t175 * t391 + (-t103 * t105 - t104 * t106 - t212 * t499 + (t505 * t18 + t17 * t325 + (t103 * t325 + t104 * t505) * qJD(4)) * pkin(3)) * m(5) + t62 * t419 + t306 * t65 - t289 * t173 - t198 * t216 + t199 * t217 + t556 * t105 + t557 * t106 + t559 * t108 + (-t27 * t29 - t28 * t30 + t302 * t7 + t559 * t55) * m(7) + t566 * t302 + t63 * t497 - t30 * t118 - t29 * t119 - t117 * t140; (-t349 + t550) * g(2) + (-t348 + t551) * g(1) + t547 * t381 + (-m(7) * (t294 + t435) - m(6) * t435 + t544) * g(3) + (t476 + t541) * t103 - t536 * t526 + t588 * t189 + t587 * t366 + (t7 * qJ(5) - t27 * t33 - t28 * t34 + t578 * t55) * m(7) + t369 * t108 + t576 + (-pkin(4) * t15 - g(1) * t387 - g(2) * t388 - qJ(5) * t13 - qJD(5) * t95 - t104 * t94 - t107 * t136) * m(6) - t95 * t483 + (t182 + t128) * t519 + t94 * t484 + t556 * t104 + t558 * qJD(5) + (-t183 + t561) * t518 - pkin(4) * t65 + t566 * qJ(5) - t34 * t118 - t33 * t119 - t136 * t140; -t558 * t315 + t371 * qJD(6) + (t140 + t371) * t366 + t65 + (-t315 * t55 - t366 * t372 + t345 + t543) * m(7) + (t107 * t366 + t315 * t95 + t15 + t543) * m(6) + t548; -t55 * (mrSges(7,1) * t167 + mrSges(7,2) * t166) + (Ifges(7,1) * t166 - t469) * t522 + t75 * t521 + (Ifges(7,5) * t166 - Ifges(7,6) * t167) * t520 - t27 * t118 + t28 * t119 - g(1) * (mrSges(7,1) * t224 - mrSges(7,2) * t225) - g(2) * (mrSges(7,1) * t226 + mrSges(7,2) * t227) + t380 * t491 + (t166 * t27 + t167 * t28) * mrSges(7,3) + t423 + (-Ifges(7,2) * t167 + t158 + t76) * t523 + t539;];
tau  = t5;
