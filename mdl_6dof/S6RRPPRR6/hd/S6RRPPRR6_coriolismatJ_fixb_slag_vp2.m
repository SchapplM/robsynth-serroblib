% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:13:11
% EndTime: 2019-03-09 09:13:24
% DurationCPUTime: 8.00s
% Computational Cost: add. (21001->554), mult. (39017->731), div. (0->0), fcn. (44567->8), ass. (0->298)
t470 = sin(pkin(10));
t471 = cos(pkin(10));
t510 = sin(qJ(5));
t511 = cos(qJ(5));
t302 = t470 * t510 - t471 * t511;
t346 = sin(qJ(6));
t344 = t346 ^ 2;
t348 = cos(qJ(6));
t345 = t348 ^ 2;
t433 = t344 + t345;
t560 = t433 * t302;
t349 = cos(qJ(2));
t347 = sin(qJ(2));
t406 = t347 * t470;
t304 = -t349 * t471 - t406;
t407 = t347 * t471;
t305 = -t349 * t470 + t407;
t248 = t510 * t304 + t511 * t305;
t485 = t346 * mrSges(7,3);
t431 = t248 * t485;
t549 = -t511 * t304 + t510 * t305;
t158 = -mrSges(7,2) * t549 - t431;
t440 = t348 * t158;
t479 = t348 * mrSges(7,3);
t161 = mrSges(7,1) * t549 - t248 * t479;
t445 = t346 * t161;
t377 = t445 / 0.2e1 - t440 / 0.2e1;
t480 = t348 * mrSges(7,2);
t487 = t346 * mrSges(7,1);
t316 = t480 + t487;
t562 = t316 * t549;
t528 = -t562 / 0.2e1;
t579 = t528 + t377;
t554 = -t349 * pkin(2) - t347 * qJ(3);
t313 = -pkin(1) + t554;
t288 = t349 * pkin(3) - t313;
t261 = -pkin(4) * t304 + t288;
t512 = t348 / 0.2e1;
t515 = -t346 / 0.2e1;
t338 = Ifges(7,6) * t346;
t500 = Ifges(7,5) * t348;
t553 = -t338 + t500;
t570 = -t248 / 0.2e1;
t557 = Ifges(6,4) * t570;
t566 = t549 / 0.2e1;
t569 = t248 / 0.2e1;
t501 = Ifges(7,4) * t348;
t321 = -Ifges(7,2) * t346 + t501;
t91 = Ifges(7,6) * t248 - t321 * t549;
t502 = Ifges(7,4) * t346;
t323 = Ifges(7,1) * t348 - t502;
t94 = Ifges(7,5) * t248 - t323 * t549;
t578 = -t94 * t512 - t91 * t515 - t261 * mrSges(6,1) - t553 * t569 - t557 - (Ifges(6,2) + Ifges(7,3)) * t566;
t577 = -mrSges(7,3) / 0.2e1;
t576 = pkin(5) * t248;
t575 = mrSges(6,2) * t566;
t574 = Ifges(7,3) * t569;
t506 = pkin(7) - qJ(4);
t396 = t506 * t471;
t382 = t349 * t396;
t395 = t506 * t470;
t266 = t347 * t395 + t382;
t508 = t304 * pkin(8);
t215 = t266 + t508;
t307 = t349 * t395;
t264 = t347 * t396 - t307;
t507 = t305 * pkin(8);
t366 = t264 - t507;
t131 = t215 * t510 - t366 * t511;
t573 = t131 * t562;
t404 = t433 * t549;
t481 = t348 * mrSges(7,1);
t486 = t346 * mrSges(7,2);
t394 = t481 - t486;
t572 = -t394 / 0.2e1 - mrSges(6,1) / 0.2e1;
t350 = -pkin(2) - pkin(3);
t311 = -qJ(3) * t470 + t471 * t350;
t309 = -pkin(4) + t311;
t312 = qJ(3) * t471 + t350 * t470;
t258 = t309 * t511 - t312 * t510;
t256 = pkin(5) - t258;
t522 = -t256 / 0.2e1;
t568 = t316 / 0.2e1;
t565 = mrSges(7,1) * t248;
t564 = mrSges(7,2) * t248;
t563 = mrSges(6,1) + t394;
t426 = -t500 / 0.2e1;
t380 = t426 + t338 / 0.2e1;
t561 = t380 * t549;
t403 = t433 * t258;
t168 = pkin(9) * t549 + t576;
t559 = t572 * t248 + t404 * t577 + t575;
t556 = Ifges(3,4) - Ifges(4,5);
t303 = -t470 * t511 - t471 * t510;
t74 = m(7) * (-0.1e1 + t433) * t303 * t302;
t555 = t74 * qJD(3);
t424 = t485 / 0.2e1;
t111 = pkin(5) * t549 - pkin(9) * t248 + t261;
t133 = t215 * t511 + t366 * t510;
t62 = t111 * t346 + t133 * t348;
t484 = t346 * t62;
t552 = t62 * t424 + t484 * t577;
t437 = t348 * t549;
t162 = -mrSges(7,3) * t437 - t565;
t444 = t346 * t549;
t159 = -mrSges(7,3) * t444 + t564;
t439 = t348 * t159;
t551 = t439 / 0.2e1 + t162 * t515;
t491 = t549 * Ifges(7,5);
t95 = t248 * t323 + t491;
t477 = t348 * t95;
t490 = t549 * Ifges(7,6);
t92 = t248 * t321 + t490;
t483 = t346 * t92;
t548 = t261 * mrSges(6,2) + Ifges(6,1) * t569 - Ifges(6,4) * t549 - t483 / 0.2e1 + t477 / 0.2e1;
t320 = Ifges(7,2) * t348 + t502;
t322 = Ifges(7,1) * t346 + t501;
t546 = -t346 * (t322 / 0.4e1 + t321 / 0.4e1) + t348 * (t323 / 0.4e1 - t320 / 0.4e1);
t411 = t321 / 0.2e1 + t322 / 0.2e1;
t412 = t320 / 0.2e1 - t323 / 0.2e1;
t545 = -t412 * t346 + t411 * t348;
t155 = t248 * t320;
t156 = t248 * t322;
t544 = -t348 * t155 / 0.4e1 - t346 * t156 / 0.4e1 + t131 * t568 + t477 / 0.4e1 - t483 / 0.4e1;
t543 = 0.2e1 * m(7);
t542 = -m(5) / 0.2e1;
t541 = m(5) / 0.2e1;
t540 = -m(6) / 0.2e1;
t539 = m(6) / 0.2e1;
t538 = -m(7) / 0.2e1;
t537 = m(7) / 0.2e1;
t536 = -pkin(5) / 0.2e1;
t336 = t349 * qJ(3);
t301 = t347 * t350 + t336;
t262 = -pkin(4) * t305 + t301;
t112 = -t168 + t262;
t265 = -t506 * t407 + t307;
t214 = t265 + t507;
t263 = t406 * t506 + t382;
t362 = t263 + t508;
t132 = t214 * t511 + t362 * t510;
t63 = t112 * t348 - t132 * t346;
t535 = -t63 / 0.2e1;
t64 = t112 * t346 + t132 * t348;
t534 = t64 / 0.2e1;
t70 = t131 * t346 + t168 * t348;
t533 = -t70 / 0.2e1;
t71 = -t131 * t348 + t168 * t346;
t532 = t71 / 0.2e1;
t531 = -t91 / 0.4e1;
t530 = -t94 / 0.4e1;
t153 = t316 * t248;
t529 = -t153 / 0.2e1;
t527 = -t549 / 0.2e1;
t259 = t309 * t510 + t312 * t511;
t257 = -pkin(9) + t259;
t521 = -t257 / 0.2e1;
t520 = t257 / 0.2e1;
t519 = -t258 / 0.2e1;
t518 = t258 / 0.2e1;
t517 = t302 / 0.2e1;
t516 = t303 / 0.2e1;
t514 = t346 / 0.2e1;
t509 = pkin(5) * t316;
t499 = Ifges(6,6) * t248;
t130 = t214 * t510 - t362 * t511;
t165 = mrSges(6,1) * t549 + mrSges(6,2) * t248;
t260 = -mrSges(5,1) * t304 + mrSges(5,2) * t305;
t315 = -t349 * mrSges(4,1) - t347 * mrSges(4,3);
t408 = m(4) * t313 + t315;
t429 = Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t467 = t130 * t131;
t61 = t111 * t348 - t133 * t346;
t1 = t130 * t153 + t573 + t64 * t158 + t62 * t159 + t63 * t161 + t61 * t162 + t262 * t165 + t301 * t260 + t408 * (pkin(2) * t347 - t336) + m(7) * (t61 * t63 + t62 * t64 + t467) + m(6) * (t132 * t133 + t261 * t262 + t467) + m(5) * (t263 * t264 + t265 * t266 + t288 * t301) + (-mrSges(3,2) * pkin(1) - mrSges(4,3) * t313 + t556 * t349) * t349 + (-t288 * mrSges(5,1) + Ifges(5,4) * t305 + (-t263 + t266) * mrSges(5,3)) * t305 + (-mrSges(5,2) * t288 - Ifges(5,4) * t304 + (-Ifges(5,1) + Ifges(5,2)) * t305 + (t264 + t265) * mrSges(5,3)) * t304 + (-mrSges(3,1) * pkin(1) + mrSges(4,1) * t313 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t349 - t556 * t347) * t347 + (t566 * t553 + t548 + (t131 - t132) * mrSges(6,3)) * t549 + (-t429 * t549 + Ifges(6,1) * t566 + Ifges(6,4) * t569 + (t130 + t133) * mrSges(6,3) + t578) * t248;
t498 = t1 * qJD(1);
t497 = t130 * mrSges(6,1);
t496 = t131 * mrSges(6,2);
t495 = t132 * mrSges(6,2);
t494 = t133 * mrSges(6,1);
t492 = t549 * mrSges(6,3);
t488 = t248 * mrSges(6,3);
t482 = t346 * t94;
t478 = t348 * t91;
t151 = t248 * t394;
t425 = t490 / 0.2e1;
t5 = -t131 * t151 + t62 * t161 + ((t95 / 0.2e1 - t155 / 0.2e1 + t491 / 0.2e1) * t346 + (t156 / 0.2e1 + t92 / 0.2e1 + t425 + t62 * mrSges(7,3)) * t348) * t248 + (-t158 - t431) * t61;
t476 = t5 * qJD(1);
t475 = t63 * t346;
t474 = t64 * t348;
t473 = t70 * t346;
t472 = t71 * t348;
t392 = t346 * t61 - t348 * t62;
t464 = t131 * t248;
t14 = (t153 + t488) * t248 + (t304 ^ 2 + t305 ^ 2) * mrSges(5,3) - (t440 - t445 - t492) * t549 + m(7) * (t392 * t549 + t464) + m(6) * (-t133 * t549 + t464) + m(5) * (-t264 * t305 + t266 * t304);
t469 = qJD(1) * t14;
t438 = t348 * t161;
t447 = t346 * t158;
t24 = (m(7) * (t348 * t61 + t484) + t447 + t438 + t165 + m(6) * t261 + m(5) * t288 + t260 - t408) * t347;
t468 = qJD(1) * t24;
t466 = t130 * t302;
t465 = t130 * t394;
t463 = t131 * t303;
t461 = t133 * t394;
t352 = (t304 * t312 - t305 * t311) * t541 + (-t248 * t258 - t259 * t549) * t539 + (t248 * t256 - t257 * t404) * t537 - t559;
t356 = t301 * t541 + t262 * t539 + (t346 * t64 + t348 * t63) * t537 + mrSges(6,1) * t570 + t575 + t159 * t514 + t162 * t512;
t16 = t305 * mrSges(5,1) + t304 * mrSges(5,2) + t352 - t356;
t460 = t16 * qJD(1);
t157 = t485 * t549 - t564;
t160 = t479 * t549 + t565;
t359 = (t346 * t71 + t348 * t70) * t537 + mrSges(6,2) * t527 + mrSges(6,1) * t569 + t157 * t514 + t160 * t512;
t360 = (-pkin(9) * t404 - t576) * t537 + t559;
t19 = -t359 + t360;
t459 = t19 * qJD(1);
t458 = t248 * t302;
t456 = t549 * t553;
t318 = Ifges(7,5) * t346 + Ifges(7,6) * t348;
t455 = t248 * t318;
t454 = t256 * t316;
t379 = -t480 / 0.2e1 - t487 / 0.2e1;
t372 = t379 * t549;
t26 = -t372 + t377;
t453 = t26 * qJD(1);
t452 = t302 * t151;
t204 = t302 * t259;
t451 = t302 * t316;
t446 = t346 * t160;
t443 = t346 * t320;
t441 = t348 * t157;
t436 = t348 * t322;
t405 = t433 * t303;
t355 = (t304 * t470 - t305 * t471) * t541 + (t303 * t549 + t458) * t539 + (t405 * t549 + t458) * t537;
t399 = m(7) * t433 * t347;
t43 = (t540 + t542) * t347 - t399 / 0.2e1 + t355;
t435 = t43 * qJD(1);
t432 = mrSges(7,3) * t472;
t430 = pkin(5) * t528;
t423 = -t256 * t303 - t257 * t560;
t420 = t303 * t514;
t415 = t569 + t570;
t414 = t566 + t527;
t398 = mrSges(7,3) * (-t345 / 0.2e1 - t344 / 0.2e1);
t391 = -t474 + t475;
t390 = t472 - t473;
t363 = -t494 / 0.2e1 - t461 / 0.2e1 - t562 * t522 + t259 * t529;
t385 = t131 * t259 + t133 * t256;
t3 = t385 * t538 + t430 + t415 * Ifges(6,6) + t414 * Ifges(6,5) + (t131 / 0.2e1 - t132 / 0.2e1) * mrSges(6,2) + (m(7) * t536 + t572) * t130 + (t531 + t91 / 0.4e1 + pkin(9) * t159 / 0.2e1 + t158 * t519 + t157 * t521 + (t532 + t534) * mrSges(7,3) + (pkin(9) * t64 / 0.4e1 - t257 * t71 / 0.4e1 - t258 * t62 / 0.4e1) * t543) * t348 + (t530 + t94 / 0.4e1 - pkin(9) * t162 / 0.2e1 + t161 * t518 + t160 * t520 + (t533 + t535) * mrSges(7,3) + (-pkin(9) * t63 / 0.4e1 + t257 * t70 / 0.4e1 + t258 * t61 / 0.4e1) * t543) * t346 + t363;
t370 = t258 * mrSges(6,2) - mrSges(7,3) * t403 + t259 * t563;
t31 = m(7) * (t256 * t259 + t257 * t403) + t370;
t389 = -t3 * qJD(1) + t31 * qJD(2);
t369 = t302 * mrSges(6,2) - mrSges(7,3) * t560 + t303 * t563;
t40 = m(4) * qJ(3) + t470 * mrSges(5,1) + t471 * mrSges(5,2) + mrSges(4,3) + m(7) * t423 + m(6) * (t258 * t303 - t204) + m(5) * (-t311 * t470 + t312 * t471) - t369;
t353 = (t263 * t471 + t265 * t470) * t542 + (-t132 * t303 + t466) * t540 + (t303 * t391 + t466) * t538;
t354 = (-t264 * t470 + t266 * t471) * t541 + (-t133 * t302 - t463) * t539 + (t302 * t392 - t463) * t537;
t9 = (mrSges(6,3) * t415 + t529 + t551) * t303 + (-mrSges(6,3) * t414 + t579) * t302 + t353 + t354;
t388 = -t9 * qJD(1) - t40 * qJD(2);
t13 = (t529 + (-t131 - t390) * t537 - t441 / 0.2e1 + t446 / 0.2e1) * t303 + ((t133 + t392) * t537 + t579) * t302;
t4 = t133 * t153 - t573 + m(7) * (t131 * t133 + t61 * t70 + t62 * t71) + t71 * t158 + t62 * t157 + t70 * t161 + t61 * t160 + (t557 - t578) * t248 + (t561 + (-Ifges(6,1) / 0.2e1 + t429) * t248 - t548) * t549;
t387 = t4 * qJD(1) + t13 * qJD(3);
t113 = -t454 + t545;
t378 = -t447 / 0.2e1 - t438 / 0.2e1;
t351 = t378 * t257 + (t257 * t398 - t546) * t248 - t456 / 0.4e1 - t151 * t522 - t544 + t552;
t367 = mrSges(7,1) * t535 + mrSges(7,2) * t534 + t574;
t7 = t351 + t367 + t561;
t386 = t7 * qJD(1) + t113 * qJD(2);
t175 = (t568 - t379) * t302;
t373 = (-t486 / 0.2e1 + t481 / 0.2e1) * t347;
t21 = -t452 / 0.2e1 + t373 + (t248 * t398 + t378) * t303;
t384 = -t21 * qJD(1) - t175 * qJD(2);
t381 = mrSges(7,1) * t533 + mrSges(7,2) * t532;
t376 = t443 / 0.2e1 - t436 / 0.2e1;
t375 = -Ifges(6,5) + t376;
t374 = t379 * t302;
t357 = (-t258 * t405 + t204 + t423) * t537;
t368 = m(7) * (pkin(5) * t303 - pkin(9) * t560);
t25 = -t368 / 0.2e1 + t357 - t369;
t371 = t13 * qJD(1) + t25 * qJD(2) + t555;
t358 = t378 * pkin(9) + t151 * t536 + t544 + t552;
t361 = pkin(9) * t398 + t546;
t11 = (t553 / 0.4e1 - t380) * t549 + (-Ifges(7,3) / 0.2e1 + t361) * t248 + t358 + t381;
t173 = t509 - t545;
t174 = (-t316 / 0.2e1 - t379) * t302;
t47 = (t522 + t536) * t316 + (mrSges(7,2) * t519 + t411) * t348 + (mrSges(7,1) * t519 - t412) * t346;
t364 = t11 * qJD(1) - t47 * qJD(2) - t174 * qJD(3) - t173 * qJD(5);
t177 = -t451 / 0.2e1 - t374;
t176 = t451 / 0.2e1 - t374;
t48 = t454 / 0.2e1 - t348 * t321 / 0.2e1 + t323 * t515 + t509 / 0.2e1 + t379 * t258 + t376;
t42 = t399 / 0.2e1 + t355 + (m(6) + m(5)) * t347 / 0.2e1;
t28 = t368 / 0.2e1 + t357;
t27 = -t372 - t377;
t22 = t452 / 0.2e1 + t158 * t420 + t373 + (mrSges(7,3) * t248 * t433 + t438) * t516;
t20 = t359 + t360;
t15 = t352 + t356;
t12 = t13 * qJD(5);
t10 = t358 + t361 * t248 + t456 / 0.4e1 + t574 + t346 * t425 + t549 * t426 - t381;
t8 = t162 * t420 + t562 * t517 + (m(4) * pkin(7) + mrSges(4,2)) * t349 + (t492 / 0.2e1 + t377) * t302 - t353 + t354 + (t471 * t304 + t305 * t470) * mrSges(5,3) + (-t439 / 0.2e1 + t529 - t488 / 0.2e1) * t303 + (-t248 * t516 + t517 * t549) * mrSges(6,3);
t6 = t351 - Ifges(7,6) * t444 / 0.2e1 + Ifges(7,5) * t437 / 0.2e1 - t367;
t2 = -t495 / 0.2e1 + (-t391 * t537 + t551) * pkin(9) + t430 + t499 / 0.2e1 - t496 / 0.2e1 - t497 / 0.2e1 + t440 * t518 + t445 * t519 + t441 * t520 + t446 * t521 - t432 / 0.2e1 - t482 / 0.4e1 - t478 / 0.4e1 + (-t475 / 0.2e1 + t474 / 0.2e1) * mrSges(7,3) + 0.2e1 * (-t443 / 0.4e1 + t436 / 0.4e1) * t549 + (-pkin(5) * t130 + t257 * t390 - t258 * t392 + t385) * t537 - t363 + t346 * t530 + t348 * t531 - t455 / 0.4e1 - t465 / 0.2e1 + t70 * t424 - (t318 / 0.4e1 - Ifges(6,6) / 0.2e1) * t248 + 0.2e1 * t566 * Ifges(6,5);
t17 = [qJD(2) * t1 + qJD(3) * t24 + qJD(4) * t14 + qJD(5) * t4 - qJD(6) * t5, t8 * qJD(3) + t15 * qJD(4) + t2 * qJD(5) + t6 * qJD(6) + t498 + (-t263 * mrSges(5,1) + t265 * mrSges(5,2) + Ifges(5,5) * t304 - Ifges(5,6) * t305 + t256 * t562 + t465 + t495 + t497 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t349 + (t91 / 0.2e1 - t64 * mrSges(7,3) + t257 * t159) * t348 + (-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t347 + (t94 / 0.2e1 + t63 * mrSges(7,3) - t257 * t162) * t346 - (-t318 / 0.2e1 + Ifges(6,6) - t259 * mrSges(6,3)) * t248 + 0.2e1 * (-t130 * t258 + t132 * t259) * t539 + 0.2e1 * (t263 * t311 + t265 * t312) * t541 + 0.2e1 * (t130 * t256 - t257 * t391) * t537 + (-t258 * mrSges(6,3) + t375) * t549 + (m(4) * t554 - t349 * mrSges(3,1) + t347 * mrSges(3,2) + t315) * pkin(7) + (t304 * t311 + t305 * t312) * mrSges(5,3)) * qJD(2), qJD(2) * t8 + qJD(4) * t42 + qJD(6) * t22 + t12 + t468, qJD(2) * t15 + qJD(3) * t42 + qJD(5) * t20 + qJD(6) * t27 + t469, t2 * qJD(2) + t20 * qJD(4) + t10 * qJD(6) + t387 + (t455 / 0.2e1 + t482 / 0.2e1 + t478 / 0.2e1 - t461 + t432 - mrSges(7,3) * t473 + t496 - t494 - t499 + t375 * t549 + (-m(7) * t133 + t562) * pkin(5) + (m(7) * t390 + t441 - t446) * pkin(9)) * qJD(5), -t476 + t6 * qJD(2) + t22 * qJD(3) + t27 * qJD(4) + t10 * qJD(5) + (-mrSges(7,1) * t62 - mrSges(7,2) * t61 - t455) * qJD(6); qJD(3) * t9 + qJD(4) * t16 - qJD(5) * t3 + qJD(6) * t7 - t498, qJD(3) * t40 + qJD(5) * t31 + qJD(6) * t113, t28 * qJD(5) + t177 * qJD(6) - t388 + t555, t460, t28 * qJD(3) + (m(7) * (-pkin(5) * t259 + pkin(9) * t403) - t370) * qJD(5) + t48 * qJD(6) + t389, t177 * qJD(3) + t48 * qJD(5) + (-t257 * t394 - t553) * qJD(6) + t386; -qJD(2) * t9 + qJD(4) * t43 - qJD(6) * t21 + t12 - t468, qJD(5) * t25 - qJD(6) * t175 + t388, t74 * qJD(5), t435 (t368 + t369) * qJD(5) + t176 * qJD(6) + t371, qJD(6) * t303 * t394 + t176 * qJD(5) + t384; -qJD(2) * t16 - qJD(3) * t43 - qJD(5) * t19 - qJD(6) * t26 - t469, -t460, -t435, 0, -t459, -qJD(6) * t316 - t453; qJD(2) * t3 + qJD(4) * t19 + qJD(6) * t11 - t387, -qJD(3) * t25 - qJD(6) * t47 - t389, -qJD(6) * t174 - t371, t459, -t173 * qJD(6) (-pkin(9) * t394 + t553) * qJD(6) + t364; -qJD(2) * t7 + qJD(3) * t21 + qJD(4) * t26 - qJD(5) * t11 + t476, qJD(3) * t175 + qJD(5) * t47 - t386, qJD(5) * t174 - t384, t453, -t364, 0;];
Cq  = t17;
