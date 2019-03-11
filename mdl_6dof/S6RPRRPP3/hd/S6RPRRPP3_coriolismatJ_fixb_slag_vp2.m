% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:12
% EndTime: 2019-03-09 04:35:25
% DurationCPUTime: 7.61s
% Computational Cost: add. (6605->537), mult. (13733->688), div. (0->0), fcn. (10785->6), ass. (0->269)
t296 = sin(qJ(4));
t298 = cos(qJ(4));
t457 = Ifges(7,6) * t298;
t458 = Ifges(6,6) * t298;
t533 = t457 - t458 + (Ifges(7,2) + Ifges(6,3)) * t296;
t529 = Ifges(7,4) + Ifges(6,5);
t528 = Ifges(5,5) + Ifges(7,5);
t468 = mrSges(6,1) + mrSges(7,1);
t524 = t296 ^ 2 + t298 ^ 2;
t537 = t524 * (mrSges(5,3) + t468);
t459 = Ifges(6,6) * t296;
t356 = t298 * Ifges(6,3) + t459;
t287 = Ifges(7,6) * t296;
t359 = t298 * Ifges(7,2) - t287;
t519 = t359 + t356;
t288 = Ifges(5,4) * t298;
t238 = Ifges(5,1) * t296 + t288;
t355 = -Ifges(7,3) * t296 + t457;
t518 = -t238 + t355;
t299 = cos(qJ(3));
t430 = t298 * t299;
t270 = mrSges(7,1) * t430;
t297 = sin(qJ(3));
t447 = t297 * mrSges(7,3);
t223 = t270 - t447;
t271 = mrSges(6,1) * t430;
t449 = t297 * mrSges(6,2);
t226 = t271 + t449;
t295 = pkin(4) + qJ(6);
t434 = t296 * t299;
t224 = mrSges(6,1) * t434 - t297 * mrSges(6,3);
t448 = t297 * mrSges(7,2);
t225 = -mrSges(7,1) * t434 + t448;
t387 = -t224 / 0.2e1 + t225 / 0.2e1;
t471 = t297 * pkin(3);
t472 = pkin(8) * t299;
t248 = t471 - t472;
t220 = t296 * t248;
t276 = sin(pkin(9)) * pkin(1) + pkin(7);
t433 = t297 * t298;
t107 = -t276 * t433 + t220;
t453 = t107 * mrSges(5,2);
t239 = t297 * t276;
t431 = t298 * t248;
t106 = t239 * t296 + t431;
t454 = t106 * mrSges(5,1);
t438 = t276 * t296;
t93 = -t431 + (-pkin(4) - t438) * t297;
t469 = t93 * mrSges(6,2);
t382 = qJ(6) + t438;
t375 = -pkin(4) - t382;
t64 = (pkin(5) * t299 - t248) * t298 + t375 * t297;
t496 = -t64 / 0.2e1;
t499 = mrSges(7,2) / 0.2e1;
t503 = m(7) / 0.2e1;
t504 = m(6) / 0.2e1;
t92 = t297 * (t276 * t298 - qJ(5)) - t220;
t71 = -pkin(5) * t434 - t92;
t535 = (-pkin(4) * t93 - qJ(5) * t92) * t504 + (qJ(5) * t71 - t295 * t64) * t503 - pkin(4) * t226 / 0.2e1 + t454 / 0.2e1 - t453 / 0.2e1 - t295 * t223 / 0.2e1 + mrSges(7,3) * t496 + t71 * t499 - t92 * mrSges(6,3) / 0.2e1 + t469 / 0.2e1 + t387 * qJ(5);
t534 = t533 * t297;
t467 = -mrSges(7,2) - mrSges(6,3);
t532 = (mrSges(5,2) + t467) * qJD(4);
t531 = -t297 / 0.2e1;
t484 = t297 / 0.2e1;
t530 = -t299 / 0.2e1;
t478 = t299 / 0.2e1;
t429 = t299 * qJ(5);
t405 = -cos(pkin(9)) * pkin(1) - pkin(2);
t470 = t297 * pkin(8);
t210 = -pkin(3) * t299 + t405 - t470;
t175 = t296 * t210;
t240 = t299 * t276;
t227 = t298 * t240;
t90 = t227 + t175;
t78 = -t90 + t429;
t466 = t78 + t90;
t526 = Ifges(7,3) * t298 + t287;
t525 = -Ifges(5,2) * t296 + t288;
t461 = Ifges(5,4) * t296;
t364 = Ifges(5,1) * t298 - t461;
t328 = t364 * t297;
t523 = t297 * t526 - t299 * t528 + t328;
t522 = -t299 * t529 + t534;
t435 = t296 * t297;
t413 = mrSges(5,3) * t435;
t444 = t299 * mrSges(5,2);
t335 = -t413 + t444;
t410 = mrSges(5,3) * t433;
t337 = -t299 * mrSges(5,1) - t410;
t521 = t296 * t335 + t298 * t337;
t412 = mrSges(6,1) * t433;
t336 = -t299 * mrSges(6,2) + t412;
t285 = t299 * mrSges(6,3);
t415 = mrSges(6,1) * t435;
t373 = t285 + t415;
t520 = t296 * t373 + t298 * t336;
t517 = -t524 * t470 / 0.2e1;
t516 = -mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t515 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t418 = Ifges(6,4) - t528;
t514 = -Ifges(5,6) + t529;
t513 = t525 - t518;
t512 = t526 - t519;
t366 = -t298 * mrSges(7,2) + t296 * mrSges(7,3);
t206 = t299 * t366;
t367 = -t296 * mrSges(6,2) - t298 * mrSges(6,3);
t208 = t299 * t367;
t368 = t296 * mrSges(5,1) + t298 * mrSges(5,2);
t209 = t299 * t368;
t511 = -t209 / 0.2e1 - t208 / 0.2e1 - t206 / 0.2e1;
t436 = t296 * qJ(5);
t323 = t295 * t298 + t436;
t497 = -Ifges(5,6) / 0.2e1;
t378 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1 + t497;
t510 = t92 * mrSges(6,1) - t71 * mrSges(7,1) - t107 * mrSges(5,3) + Ifges(5,6) * t531 + t297 * t378 + t533 * t478 + t484 * t529 + t525 * t530;
t508 = t297 ^ 2;
t507 = 2 * qJD(3);
t506 = 0.2e1 * qJD(4);
t505 = m(5) / 0.2e1;
t500 = -mrSges(6,2) / 0.2e1;
t498 = -Ifges(6,4) / 0.2e1;
t495 = m(6) + m(7);
t494 = pkin(5) + pkin(8);
t487 = t296 / 0.2e1;
t483 = -t298 / 0.2e1;
t481 = t298 / 0.2e1;
t353 = pkin(4) * t298 + t436;
t228 = -pkin(3) - t353;
t477 = m(6) * t228;
t290 = t296 * pkin(4);
t441 = qJ(5) * t298;
t231 = t290 - t441;
t476 = m(6) * t231;
t204 = -pkin(3) - t323;
t475 = m(7) * t204;
t346 = qJ(6) * t296 - t441;
t219 = t290 + t346;
t474 = m(7) * t219;
t473 = m(7) * t295;
t291 = t299 * pkin(4);
t176 = t298 * t210;
t89 = t240 * t296 - t176;
t79 = t291 + t89;
t465 = t79 - t89;
t464 = m(7) * qJD(5);
t463 = m(7) * qJD(6);
t452 = t296 * mrSges(7,1);
t451 = t297 * mrSges(5,1);
t450 = t297 * mrSges(5,2);
t365 = t296 * mrSges(7,2) + t298 * mrSges(7,3);
t329 = t365 * t297;
t315 = t329 / 0.2e1;
t229 = t298 * mrSges(6,2) - t296 * mrSges(6,3);
t330 = t297 * t229;
t316 = -t330 / 0.2e1;
t369 = t298 * mrSges(5,1) - t296 * mrSges(5,2);
t331 = t369 * t297;
t286 = t299 * mrSges(7,2);
t414 = mrSges(7,1) * t435;
t371 = -t286 - t414;
t284 = t299 * mrSges(7,3);
t411 = mrSges(7,1) * t433;
t372 = t284 + t411;
t416 = pkin(5) * t433;
t376 = -t176 + t416;
t394 = -t433 / 0.2e1;
t400 = t435 / 0.2e1;
t420 = t508 / 0.2e1;
t340 = -t291 - t376;
t48 = t299 * t382 - t340;
t7 = t331 * t478 + t371 * t400 + t521 * t484 + (t315 + t316) * t299 + (t372 + m(7) * (t299 * t375 - t376 + t48)) * t394 + (m(6) * ((t465 - t291) * t298 + (-t429 + t466) * t296) + t520) * t531 + t420 * t537;
t443 = t7 * qJD(1);
t442 = -t369 - mrSges(4,1);
t428 = pkin(4) * t435 + t239;
t108 = -qJ(5) * t433 + t428;
t205 = t366 * t297;
t207 = t367 * t297;
t425 = t285 + t286;
t341 = -pkin(5) * t435 + t90;
t61 = t341 - t429;
t86 = t297 * t346 + t428;
t13 = t425 * t299 + ((-t205 - t207) * t298 + t468 * t434) * t297 + m(6) * (-t108 * t433 + t299 * t78) + m(7) * (-t299 * t61 - t433 * t86);
t440 = t13 * qJD(1);
t439 = t228 * t297;
t28 = t205 * t435 + m(7) * (t299 * t48 + t435 * t86) + t299 * t372;
t437 = t28 * qJD(1);
t432 = t297 * t299;
t269 = t298 * t429;
t427 = -pkin(4) * t434 + t269;
t426 = t524 * t472;
t424 = qJD(3) * t299;
t423 = qJD(4) * t297;
t422 = qJD(4) * t298;
t419 = m(6) / 0.4e1 + m(7) / 0.4e1;
t408 = t500 + mrSges(7,3) / 0.2e1;
t109 = t240 - t427;
t402 = -t435 / 0.2e1;
t396 = t365 * t484;
t395 = t239 / 0.2e1;
t388 = -t207 / 0.2e1 - t205 / 0.2e1;
t386 = -t365 + t475;
t377 = Ifges(7,5) / 0.2e1 + t498 + Ifges(5,5) / 0.2e1;
t374 = t427 * t504;
t370 = (-t229 / 0.2e1 + t365 / 0.2e1) * t298;
t237 = t298 * Ifges(5,2) + t461;
t361 = -Ifges(6,2) * t298 + t459;
t360 = Ifges(6,2) * t296 + t458;
t221 = -mrSges(5,3) * t434 - t450;
t222 = -mrSges(5,3) * t430 + t451;
t303 = t93 * mrSges(6,1) + t64 * mrSges(7,1) - t106 * mrSges(5,3) + Ifges(6,4) * t531 + t377 * t297 + t361 * t530 + t528 * t484 + (t526 + t364) * t478;
t87 = qJ(6) * t434 + t109;
t3 = t108 * t208 + t109 * t207 + t87 * t205 + t86 * t206 + t90 * t221 - t89 * t222 + t48 * t223 + t78 * t224 + t61 * t225 + t79 * t226 + t64 * t284 + t92 * t285 - t71 * t286 + m(6) * (t108 * t109 + t78 * t92 + t79 * t93) + m(7) * (t48 * t64 + t61 * t71 + t86 * t87) + m(5) * (-t89 * t106 + t90 * t107) + (t405 * mrSges(4,1) - Ifges(4,4) * t297 + t276 * t209 + t296 * t510 + t303 * t298) * t297 + (Ifges(4,4) * t299 + t405 * mrSges(4,2) + t453 - t469 - t454 + t418 * t430 - t514 * t434 + (m(5) * t276 ^ 2 + Ifges(4,1) - Ifges(4,2) + (t276 * mrSges(5,2) + (Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(7,3) / 0.2e1) * t298) * t298 + (t276 * mrSges(5,1) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t296 + (-Ifges(6,6) + Ifges(7,6) - Ifges(5,4)) * t298) * t296 - t515) * t297) * t299;
t332 = m(5) * (-t106 * t296 + t107 * t298);
t349 = t296 * t93 - t298 * t92;
t350 = -t296 * t48 - t298 * t61;
t6 = ((-t285 / 0.2e1 - t286 / 0.2e1 + t444 / 0.2e1) * t298 + (t284 / 0.2e1 + (mrSges(5,1) / 0.2e1 + t500) * t299) * t296 + (t296 * t79 - t298 * t78 - t109) * t504 + (-t350 - t87) * t503 + (t296 * t89 + t298 * t90) * t505 - m(5) * t240 / 0.2e1 + t511) * t299 + ((t221 / 0.2e1 + t450 / 0.2e1 + t387) * t298 + (-t222 / 0.2e1 + t226 / 0.2e1 + t223 / 0.2e1 + t451 / 0.2e1) * t296 + (t108 + t349) * t504 + (t296 * t64 + t298 * t71 + t86) * t503 + t332 / 0.2e1 + m(5) * t395 - t388) * t297;
t352 = t3 * qJD(1) + t6 * qJD(2);
t306 = -Ifges(5,6) * t299 + t297 * t525;
t327 = t361 * t297;
t311 = -Ifges(6,4) * t299 + t327;
t312 = t297 * t323;
t313 = -t89 - t416;
t324 = t353 * t297;
t5 = -t86 * t329 - t313 * t371 - t341 * t372 - m(7) * (pkin(5) * t350 + t323 * t86) * t297 - t207 * t324 - t205 * t312 + t61 * t411 + t48 * t414 + t79 * t415 + t306 * t433 / 0.2e1 - t78 * t412 + t311 * t402 + t298 * t360 * t420 - t508 * t276 * t369 + t523 * t400 + t522 * t394 + (-m(6) * t324 + t330) * t108 + (-m(6) * t79 - m(7) * t48 - t336 + t337 + t410) * t90 + (-m(6) * t78 + m(7) * t61 + t335 - t373 + t413) * t89 + (t296 * t418 + t298 * t514) * t432 / 0.2e1 - (t518 * t298 + (t237 + t519) * t296) * t508 / 0.2e1;
t351 = -t5 * qJD(1) - t7 * qJD(2);
t32 = 0.4e1 * (m(5) / 0.4e1 + t419) * (-0.1e1 + t524) * t432;
t348 = t6 * qJD(1) + t32 * qJD(2);
t38 = (t229 + t386 + t477) * t296;
t247 = t494 * t298;
t302 = t388 * t296 + (-t296 * t108 + (-t439 - t472) * t298) * t504 + (-t204 * t433 - t247 * t299 - t296 * t86) * t503;
t339 = -m(6) * t93 / 0.2e1 + m(7) * t496;
t9 = -t270 - t271 + (t370 + t408) * t297 + t302 + t339;
t347 = -qJD(1) * t9 + qJD(3) * t38;
t246 = t494 * t296;
t304 = (t204 * t435 + t246 * t299 - t298 * t86) * t503 + t205 * t483;
t338 = -m(7) * t71 / 0.2e1 - t448 / 0.2e1;
t24 = (t299 * mrSges(7,1) - t396) * t296 + t304 + t338;
t83 = t386 * t298;
t345 = -qJD(1) * t24 + qJD(3) * t83;
t344 = t246 * t296 + t247 * t298;
t283 = -0.2e1 * t429;
t23 = 0.2e1 * (t227 / 0.4e1 + t175 / 0.4e1 - t90 / 0.4e1) * m(6) + 0.2e1 * t419 * t283 - t425;
t250 = qJ(5) * t495 - t467;
t343 = qJD(1) * t23 + qJD(4) * t250;
t249 = mrSges(7,3) + t473;
t30 = (-t291 + (-qJ(6) - t295) * t299) * t503 - t284;
t342 = qJD(1) * t30 + qJD(4) * t249;
t333 = m(6) * t353;
t322 = -qJ(5) * t468 + t514;
t321 = pkin(4) * mrSges(6,1) + t295 * mrSges(7,1) + t418;
t300 = -(t327 + t311) * t298 / 0.4e1 - pkin(3) * t331 / 0.2e1 - t246 * t371 / 0.2e1 + t247 * t372 / 0.2e1 + t229 * t324 / 0.2e1 + t237 * t394 + t368 * t395 - t365 * t312 / 0.2e1 + (t341 / 0.2e1 - t61 / 0.2e1) * t452 + t517 * mrSges(5,3) + (t328 + t523) * t298 / 0.4e1 - (t296 * t514 - t298 * t418) * t299 / 0.4e1 + (-t360 / 0.2e1 + t518 / 0.4e1 - t513 / 0.4e1) * t435 + t204 * t315 + t228 * t316 + (t522 + t534) * t296 / 0.4e1 + (-t519 / 0.4e1 + t512 / 0.4e1) * t433 - t296 * t306 / 0.4e1 + (t231 * t108 + t228 * t324) * t504 + (t520 / 0.2e1 - t521 / 0.2e1 + (t296 * t466 + t298 * t465) * t504) * pkin(8) + ((t313 + t48) * t481 + t246 * t402 + t247 * t394) * mrSges(7,1) + (t466 * t487 + t517 + (-t89 / 0.2e1 + t79 / 0.2e1) * t298) * mrSges(6,1) + t86 * t366 / 0.2e1 + t108 * t367 / 0.2e1 + (t219 * t86 + (t48 - t89) * t247 + (-t61 + t90) * t246 + (-pkin(5) * t344 + t204 * t323) * t297) * t503 + t219 * t205 / 0.2e1 + t231 * t207 / 0.2e1;
t1 = t515 * t484 - t300 + (t497 + t529 / 0.2e1) * t434 + (t498 + t528 / 0.2e1) * t430 + t535;
t314 = (t474 + t476) * t478;
t18 = t374 + (-t295 * t434 + t269) * t503 + t314;
t8 = -pkin(3) * t368 + t231 * t229 - t219 * t365 + (t367 + t476) * t228 + (t366 + t474) * t204 - (t361 + t237) * t296 / 0.2e1 + (-t360 + t533) * t483 + t513 * t481 + (t364 + t512) * t487;
t317 = t1 * qJD(1) + t18 * qJD(2) - t8 * qJD(3);
t241 = (qJD(1) * t299 - qJD(4)) * m(7);
t212 = t495 * t433;
t211 = t495 * t434;
t174 = -m(7) * t246 - t452;
t96 = m(7) * t247 + (m(6) * pkin(8) + t468) * t298;
t29 = -t372 + ((-t295 - t382) * t299 + t340 + t313) * t503;
t26 = -t296 * t396 + t304 - t338;
t20 = -t468 * t435 - t425 + (0.2e1 * t90 + t283) * t504 + (0.2e1 * t341 + t283) * t503;
t19 = t374 + t269 * t503 + ((-mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1 + t499) * t298 + (-mrSges(5,1) / 0.2e1 - t473 / 0.2e1 - t408) * t296) * t299 - t314 + t511;
t11 = -t447 / 0.2e1 + t449 / 0.2e1 + t297 * t370 + t302 - t339;
t4 = qJD(3) * t6 - qJD(4) * t7;
t2 = (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1) * t297 + (t296 * t378 + t298 * t377) * t299 + t300 + t535;
t10 = [qJD(3) * t3 - qJD(4) * t5 + qJD(5) * t13 + qJD(6) * t28, t4, t2 * qJD(4) + t11 * qJD(5) + t26 * qJD(6) + (t109 * t477 / 0.2e1 + (t204 * t87 + t246 * t64 + t247 * t71) * t503) * t507 + (Ifges(4,5) + (-t355 / 0.2e1 + t360 / 0.2e1 + t238 / 0.2e1) * t298 + (-t356 / 0.2e1 - t359 / 0.2e1 - t237 / 0.2e1) * t296 + (-m(5) * pkin(3) + t442) * t276) * t424 + t352 + (mrSges(4,2) * t239 - Ifges(4,6) * t297 - pkin(3) * t209 + t109 * t229 + t204 * t206 + t228 * t208 + t246 * t223 + t247 * t225 - t87 * t365 - t510 * t298 + t303 * t296 + ((t221 - t224) * t298 + (-t222 + t226) * t296 + t332 + m(6) * t349) * pkin(8)) * qJD(3), t2 * qJD(3) + t20 * qJD(5) + t29 * qJD(6) + (t322 * t298 + t321 * t296 + (m(7) * (t295 * t296 - t441) + t366) * pkin(5)) * t423 + t351 + (t516 * qJD(4) + (-pkin(4) * t504 - t295 * t503) * t506) * t90 + ((-t504 - t503) * t506 * qJ(5) + t532) * t89, qJD(3) * t11 + qJD(4) * t20 + t440, qJD(3) * t26 + qJD(4) * t29 + t437; t4, t32 * qJD(3), t19 * qJD(4) + t211 * qJD(5) + (t229 - t365 + t442) * qJD(3) * t297 + ((t426 + t439) * t504 + t475 * t484 + (t426 - t471) * t505) * t507 + (t298 * t463 + (m(7) * t344 - mrSges(4,2) + t537) * qJD(3)) * t299 + t348, -t443 + t19 * qJD(3) + t212 * qJD(5) + (t516 * t422 + (t532 - t463) * t296 + (-t333 / 0.2e1 - t323 * t503) * t506) * t297, qJD(3) * t211 + qJD(4) * t212 (-t296 * t423 + t298 * t424) * m(7); -qJD(4) * t1 + qJD(5) * t9 + qJD(6) * t24 - t352, -qJD(4) * t18 - t348, qJD(4) * t8 - qJD(5) * t38 - qJD(6) * t83, t96 * qJD(5) + t174 * qJD(6) - t321 * t422 - t317 + (m(7) * (-qJ(5) * t246 - t247 * t295) - t246 * mrSges(7,2) - t247 * mrSges(7,3) + t322 * t296 + (t229 - t333 - t369) * pkin(8)) * qJD(4), qJD(4) * t96 - t347, qJD(4) * t174 - t345; qJD(3) * t1 + qJD(5) * t23 + qJD(6) * t30 - t351, qJD(3) * t18 + t443, t317, qJD(5) * t250 + qJD(6) * t249, t343, t342; -t9 * qJD(3) - t23 * qJD(4) + t299 * t463 - t440, 0, t347, -t343 - t463, 0, t241; -t24 * qJD(3) - t30 * qJD(4) - t299 * t464 - t437, 0, t345, -t342 + t464, -t241, 0;];
Cq  = t10;
