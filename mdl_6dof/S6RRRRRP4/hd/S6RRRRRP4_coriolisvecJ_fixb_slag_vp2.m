% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:12:07
% EndTime: 2019-03-10 01:12:37
% DurationCPUTime: 15.53s
% Computational Cost: add. (18468->696), mult. (44218->908), div. (0->0), fcn. (31825->8), ass. (0->329)
t524 = Ifges(6,1) + Ifges(7,1);
t522 = Ifges(7,4) + Ifges(6,5);
t353 = cos(qJ(2));
t484 = -pkin(8) - pkin(7);
t330 = t484 * t353;
t318 = qJD(1) * t330;
t348 = sin(qJ(3));
t302 = t348 * t318;
t349 = sin(qJ(2));
t328 = t484 * t349;
t317 = qJD(1) * t328;
t307 = qJD(2) * pkin(2) + t317;
t352 = cos(qJ(3));
t263 = t352 * t307 + t302;
t345 = qJD(2) + qJD(3);
t234 = -t345 * pkin(3) - t263;
t314 = t348 * t353 + t352 * t349;
t301 = t314 * qJD(1);
t347 = sin(qJ(4));
t351 = cos(qJ(4));
t364 = t301 * t347 - t345 * t351;
t189 = pkin(4) * t364 + t234;
t402 = qJD(1) * t353;
t403 = qJD(1) * t349;
t300 = -t348 * t403 + t352 * t402;
t296 = qJD(4) - t300;
t288 = qJD(5) + t296;
t340 = -pkin(2) * t353 - pkin(1);
t326 = qJD(1) * t340;
t225 = -t300 * pkin(3) - t301 * pkin(9) + t326;
t407 = t352 * t318;
t264 = t348 * t307 - t407;
t235 = pkin(9) * t345 + t264;
t146 = t351 * t225 - t235 * t347;
t277 = t301 * t351 + t345 * t347;
t123 = -pkin(10) * t277 + t146;
t109 = pkin(4) * t296 + t123;
t350 = cos(qJ(5));
t147 = t347 * t225 + t351 * t235;
t124 = -pkin(10) * t364 + t147;
t346 = sin(qJ(5));
t425 = t124 * t346;
t36 = t109 * t350 - t425;
t508 = qJD(6) - t36;
t34 = -pkin(5) * t288 + t508;
t200 = t346 * t277 + t350 * t364;
t357 = t350 * t277 - t346 * t364;
t73 = t200 * pkin(5) - qJ(6) * t357 + t189;
t566 = t189 * mrSges(6,2) + t34 * mrSges(7,2) - t36 * mrSges(6,3) - t73 * mrSges(7,3);
t435 = Ifges(6,4) * t357;
t102 = -Ifges(6,2) * t200 + Ifges(6,6) * t288 + t435;
t410 = t350 * t124;
t37 = t109 * t346 + t410;
t35 = qJ(6) * t288 + t37;
t197 = Ifges(7,5) * t357;
t99 = Ifges(7,6) * t288 + Ifges(7,3) * t200 + t197;
t565 = t189 * mrSges(6,1) + t73 * mrSges(7,1) - t35 * mrSges(7,2) - t37 * mrSges(6,3) - t102 / 0.2e1 + t99 / 0.2e1;
t523 = -Ifges(6,4) + Ifges(7,5);
t521 = Ifges(7,2) + Ifges(6,3);
t520 = -Ifges(6,6) + Ifges(7,6);
t258 = pkin(3) * t301 - pkin(9) * t300;
t171 = t351 * t258 - t263 * t347;
t418 = t300 * t351;
t378 = t301 * pkin(4) - pkin(10) * t418;
t483 = -pkin(10) - pkin(9);
t388 = qJD(4) * t483;
t564 = t351 * t388 - t171 - t378;
t172 = t347 * t258 + t351 * t263;
t419 = t300 * t347;
t396 = pkin(10) * t419;
t563 = -t347 * t388 + t172 - t396;
t463 = t288 / 0.2e1;
t472 = t357 / 0.2e1;
t475 = t200 / 0.2e1;
t476 = -t200 / 0.2e1;
t562 = -Ifges(6,2) * t476 + Ifges(7,3) * t475 + t463 * t520 + t472 * t523 + t565;
t464 = -t288 / 0.2e1;
t473 = -t357 / 0.2e1;
t561 = -Ifges(6,4) * t475 - Ifges(7,5) * t476 - t522 * t464 - t524 * t473 + t566;
t560 = -Ifges(6,2) * t475 + Ifges(7,3) * t476 + t464 * t520 + t473 * t523 - t565;
t198 = Ifges(6,4) * t200;
t434 = Ifges(7,5) * t200;
t507 = t522 * t288 + t357 * t524 - t198 + t434;
t287 = pkin(4) * t419;
t400 = qJD(4) * t347;
t342 = pkin(4) * t400;
t559 = t342 - t287;
t517 = t345 * Ifges(4,5);
t558 = -t517 / 0.2e1 - t326 * mrSges(4,2);
t389 = qJD(2) * t484;
t379 = qJD(1) * t389;
t308 = t349 * t379;
t361 = t353 * t379;
t178 = qJD(3) * t264 + t308 * t348 - t352 * t361;
t312 = t348 * t349 - t352 * t353;
t270 = t345 * t312;
t252 = t270 * qJD(1);
t186 = -qJD(4) * t277 + t252 * t347;
t105 = -pkin(4) * t186 + t178;
t359 = t364 * qJD(4);
t185 = -t252 * t351 - t359;
t67 = -qJD(5) * t200 + t350 * t185 + t346 * t186;
t68 = qJD(5) * t357 + t346 * t185 - t350 * t186;
t14 = pkin(5) * t68 - qJ(6) * t67 - qJD(6) * t357 + t105;
t271 = t345 * t314;
t253 = t271 * qJD(1);
t401 = qJD(2) * t349;
t393 = pkin(2) * t401;
t158 = pkin(3) * t253 + pkin(9) * t252 + qJD(1) * t393;
t177 = qJD(3) * t263 + t352 * t308 + t348 * t361;
t43 = -qJD(4) * t147 + t351 * t158 - t177 * t347;
t27 = pkin(4) * t253 - pkin(10) * t185 + t43;
t399 = qJD(4) * t351;
t42 = t347 * t158 + t351 * t177 + t225 * t399 - t235 * t400;
t29 = pkin(10) * t186 + t42;
t397 = qJD(5) * t350;
t398 = qJD(5) * t346;
t7 = t109 * t397 - t124 * t398 + t346 * t27 + t350 * t29;
t2 = qJ(6) * t253 + qJD(6) * t288 + t7;
t469 = t253 / 0.2e1;
t488 = t68 / 0.2e1;
t489 = -t68 / 0.2e1;
t490 = t67 / 0.2e1;
t553 = mrSges(6,1) * t105 + mrSges(7,1) * t14 - mrSges(7,2) * t2 - Ifges(6,4) * t67 / 0.2e1 - Ifges(6,6) * t253 / 0.2e1 + 0.2e1 * Ifges(7,3) * t488 - mrSges(6,3) * t7 + (-t489 + t488) * Ifges(6,2) + (t523 + Ifges(7,5)) * t490 + (t520 + Ifges(7,6)) * t469;
t519 = Ifges(4,6) * t345;
t552 = -t326 * mrSges(4,1) - t146 * mrSges(5,1) - t36 * mrSges(6,1) + t34 * mrSges(7,1) + t147 * mrSges(5,2) + t37 * mrSges(6,2) + t519 / 0.2e1 - t35 * mrSges(7,3);
t551 = Ifges(6,4) * t476 + Ifges(7,5) * t475 + t463 * t522 + t472 * t524 + t566;
t550 = t507 / 0.2e1;
t548 = t522 * t253 + t523 * t68 + t524 * t67;
t327 = t483 * t347;
t344 = t351 * pkin(10);
t329 = pkin(9) * t351 + t344;
t280 = t327 * t346 + t329 * t350;
t512 = -qJD(5) * t280 + t346 * t563 + t350 * t564;
t362 = t350 * t327 - t329 * t346;
t511 = qJD(5) * t362 + t346 * t564 - t350 * t563;
t412 = t346 * t351;
t313 = t347 * t350 + t412;
t223 = t313 * t300;
t413 = t346 * t347;
t311 = -t350 * t351 + t413;
t224 = t311 * t300;
t497 = qJD(4) + qJD(5);
t268 = t497 * t311;
t269 = t497 * t313;
t545 = -qJD(6) * t313 + t559 + (-t224 + t268) * qJ(6) + (-t223 + t269) * pkin(5);
t265 = t317 * t348 - t407;
t433 = pkin(2) * qJD(3);
t544 = t348 * t433 - t265;
t276 = Ifges(5,4) * t364;
t175 = Ifges(5,1) * t277 + Ifges(5,5) * t296 - t276;
t409 = t351 * t175;
t438 = Ifges(5,4) * t277;
t174 = -t364 * Ifges(5,2) + Ifges(5,6) * t296 + t438;
t411 = t347 * t174;
t543 = t409 / 0.2e1 - t411 / 0.2e1;
t421 = t270 * t347;
t542 = t314 * t399 - t421;
t541 = -t43 * t347 + t351 * t42;
t516 = t364 * Ifges(5,6);
t518 = Ifges(5,3) * t296;
t539 = Ifges(5,5) * t277 + t200 * t520 + t288 * t521 + t357 * t522 - t516 + t518;
t117 = pkin(5) * t357 + qJ(6) * t200;
t528 = t364 / 0.2e1;
t49 = t123 * t350 - t425;
t515 = pkin(4) * t397 + qJD(6) - t49;
t293 = t301 * qJ(6);
t514 = -t293 + t511;
t450 = pkin(5) * t301;
t513 = t450 - t512;
t510 = -t264 + t545;
t376 = mrSges(5,1) * t347 + mrSges(5,2) * t351;
t509 = t234 * t376;
t506 = t544 + t545;
t262 = t312 * pkin(3) - t314 * pkin(9) + t340;
t281 = t328 * t348 - t330 * t352;
t190 = t351 * t262 - t281 * t347;
t415 = t314 * t351;
t136 = pkin(4) * t312 - pkin(10) * t415 + t190;
t274 = t351 * t281;
t191 = t347 * t262 + t274;
t416 = t314 * t347;
t156 = -pkin(10) * t416 + t191;
t505 = t346 * t136 + t350 * t156;
t504 = Ifges(5,5) * t185 + Ifges(5,6) * t186;
t503 = t544 + t559;
t502 = t352 * t328 + t330 * t348;
t501 = t521 * t253 + t520 * t68 + t522 * t67;
t500 = -t146 * t347 + t147 * t351;
t108 = -mrSges(5,1) * t186 + mrSges(5,2) * t185;
t499 = m(5) * t178 + t108;
t498 = qJD(4) * (-t146 * t351 - t147 * t347);
t8 = -qJD(5) * t37 + t27 * t350 - t29 * t346;
t184 = pkin(3) * t271 + pkin(9) * t270 + t393;
t321 = t349 * t389;
t322 = t353 * t389;
t205 = qJD(3) * t502 + t321 * t352 + t322 * t348;
t381 = t351 * t184 - t205 * t347;
t420 = t270 * t351;
t33 = pkin(10) * t420 + pkin(4) * t271 + (-t274 + (pkin(10) * t314 - t262) * t347) * qJD(4) + t381;
t69 = t347 * t184 + t351 * t205 + t262 * t399 - t281 * t400;
t40 = -pkin(10) * t542 + t69;
t12 = -qJD(5) * t505 + t33 * t350 - t346 * t40;
t496 = t314 * t497;
t493 = m(4) / 0.2e1;
t482 = pkin(1) * mrSges(3,1);
t481 = pkin(1) * mrSges(3,2);
t478 = t185 / 0.2e1;
t477 = t186 / 0.2e1;
t466 = -t277 / 0.2e1;
t465 = t277 / 0.2e1;
t462 = -t296 / 0.2e1;
t460 = t300 / 0.2e1;
t458 = t301 / 0.2e1;
t455 = -t347 / 0.2e1;
t454 = t351 / 0.2e1;
t453 = m(4) * t326;
t451 = pkin(2) * t352;
t449 = t42 * mrSges(5,2);
t448 = t43 * mrSges(5,1);
t336 = pkin(2) * t348 + pkin(9);
t446 = -pkin(10) - t336;
t44 = -mrSges(7,2) * t68 + mrSges(7,3) * t253;
t47 = -mrSges(6,2) * t253 - mrSges(6,3) * t68;
t445 = t44 + t47;
t45 = mrSges(6,1) * t253 - mrSges(6,3) * t67;
t46 = -t253 * mrSges(7,1) + t67 * mrSges(7,2);
t444 = t46 - t45;
t443 = mrSges(4,3) * t300;
t442 = mrSges(6,3) * t200;
t441 = mrSges(6,3) * t357;
t440 = Ifges(3,4) * t349;
t439 = Ifges(4,4) * t301;
t437 = Ifges(5,4) * t347;
t436 = Ifges(5,4) * t351;
t432 = t301 * mrSges(4,3);
t429 = Ifges(3,5) * qJD(2);
t428 = Ifges(3,6) * qJD(2);
t427 = qJD(2) * mrSges(3,1);
t426 = qJD(2) * mrSges(3,2);
t422 = t178 * t502;
t230 = pkin(2) * t403 + t258;
t266 = t317 * t352 + t302;
t167 = t351 * t230 - t266 * t347;
t126 = t167 + t378;
t168 = t347 * t230 + t351 * t266;
t135 = -t396 + t168;
t59 = t346 * t126 + t350 * t135;
t159 = -mrSges(7,2) * t200 + mrSges(7,3) * t288;
t160 = -mrSges(6,2) * t288 - t442;
t406 = t159 + t160;
t161 = mrSges(6,1) * t288 - t441;
t162 = -mrSges(7,1) * t288 + mrSges(7,2) * t357;
t405 = t161 - t162;
t404 = -mrSges(4,1) * t345 + mrSges(5,1) * t364 + t277 * mrSges(5,2) + t432;
t392 = t352 * t433;
t390 = Ifges(5,3) * t253 + t504;
t339 = -pkin(4) * t351 - pkin(3);
t386 = t429 / 0.2e1;
t385 = -t428 / 0.2e1;
t382 = qJD(4) * t446;
t226 = pkin(4) * t416 - t502;
t377 = mrSges(5,1) * t351 - mrSges(5,2) * t347;
t375 = Ifges(5,1) * t351 - t437;
t374 = Ifges(5,1) * t347 + t436;
t373 = -Ifges(5,2) * t347 + t436;
t372 = Ifges(5,2) * t351 + t437;
t371 = Ifges(5,5) * t351 - Ifges(5,6) * t347;
t370 = Ifges(5,5) * t347 + Ifges(5,6) * t351;
t58 = t126 * t350 - t135 * t346;
t129 = mrSges(5,1) * t253 - mrSges(5,3) * t185;
t130 = -mrSges(5,2) * t253 + mrSges(5,3) * t186;
t367 = -t129 * t347 + t130 * t351;
t78 = t136 * t350 - t156 * t346;
t309 = t446 * t347;
t310 = t336 * t351 + t344;
t363 = t350 * t309 - t310 * t346;
t255 = t309 * t346 + t310 * t350;
t11 = t136 * t397 - t156 * t398 + t346 * t33 + t350 * t40;
t259 = pkin(5) * t311 - qJ(6) * t313 + t339;
t358 = -t347 * t392 + t351 * t382;
t206 = qJD(3) * t281 + t321 * t348 - t352 * t322;
t4 = -pkin(5) * t253 - t8;
t356 = t8 * mrSges(6,1) - t4 * mrSges(7,1) - t7 * mrSges(6,2) + t2 * mrSges(7,3) + t501;
t128 = pkin(4) * t542 + t206;
t355 = m(5) * (t498 + t541);
t231 = Ifges(4,2) * t300 + t439 + t519;
t294 = Ifges(4,4) * t300;
t232 = t301 * Ifges(4,1) + t294 + t517;
t82 = t185 * Ifges(5,4) + t186 * Ifges(5,2) + t253 * Ifges(5,6);
t83 = Ifges(5,1) * t185 + Ifges(5,4) * t186 + Ifges(5,5) * t253;
t354 = t560 * t223 + t561 * t224 + t562 * t269 + (t548 / 0.2e1 + t105 * mrSges(6,2) + t4 * mrSges(7,2) - t8 * mrSges(6,3) - t14 * mrSges(7,3) + Ifges(6,4) * t489 + Ifges(7,5) * t488 + t522 * t469 + t524 * t490) * t313 + (t277 * t375 + t296 * t371) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t301 + t232 + t294 + t409) * t300 / 0.2e1 - t551 * t268 + (Ifges(5,5) * t466 + Ifges(5,6) * t528 + Ifges(6,6) * t475 + Ifges(7,6) * t476 + Ifges(5,3) * t462 + t521 * t464 + t522 * t473 + t552) * t301 + (t371 * t462 + t373 * t528 + t375 * t466 - t509 + t558) * t300 - (Ifges(4,1) * t300 - t439 + t539) * t301 / 0.2e1 - t373 * t359 / 0.2e1 + (t146 * t418 + t147 * t419 + t541) * mrSges(5,3) + (t509 + t543) * qJD(4) + t82 * t454 + t507 * (-t268 / 0.2e1 + t224 / 0.2e1) + t553 * t311 + t263 * t443 + (-t377 - mrSges(4,1)) * t178 + t231 * t458 + t411 * t460 - t177 * mrSges(4,2) - Ifges(4,5) * t252 - Ifges(4,6) * t253 + t370 * t469 + t372 * t477 + t374 * t478 + t347 * t83 / 0.2e1;
t341 = Ifges(3,4) * t402;
t337 = -pkin(4) * t350 - pkin(5);
t333 = pkin(4) * t346 + qJ(6);
t325 = mrSges(3,3) * t402 - t426;
t324 = -mrSges(3,3) * t403 + t427;
t323 = t339 - t451;
t299 = Ifges(3,1) * t403 + t341 + t429;
t298 = t428 + (t353 * Ifges(3,2) + t440) * qJD(1);
t283 = -mrSges(4,2) * t345 + t443;
t282 = t347 * t382 + t351 * t392;
t257 = -mrSges(4,1) * t300 + mrSges(4,2) * t301;
t241 = t311 * t314;
t240 = t313 * t314;
t237 = t259 - t451;
t218 = mrSges(5,1) * t296 - mrSges(5,3) * t277;
t217 = -t296 * mrSges(5,2) - mrSges(5,3) * t364;
t210 = t287 + t264;
t141 = qJD(5) * t255 + t282 * t346 - t350 * t358;
t140 = qJD(5) * t363 + t350 * t282 + t346 * t358;
t120 = pkin(5) * t240 + qJ(6) * t241 + t226;
t119 = mrSges(6,1) * t200 + mrSges(6,2) * t357;
t118 = mrSges(7,1) * t200 - mrSges(7,3) * t357;
t97 = pkin(4) * t277 + t117;
t96 = -t270 * t412 + (t415 * t497 - t421) * t350 - t413 * t496;
t95 = t311 * t270 - t313 * t496;
t72 = -pkin(5) * t312 - t78;
t71 = qJ(6) * t312 + t505;
t70 = -qJD(4) * t191 + t381;
t51 = -t58 - t450;
t50 = t293 + t59;
t48 = t123 * t346 + t410;
t25 = mrSges(6,1) * t68 + mrSges(6,2) * t67;
t24 = mrSges(7,1) * t68 - mrSges(7,3) * t67;
t19 = pkin(5) * t96 - qJ(6) * t95 + qJD(6) * t241 + t128;
t10 = -pkin(5) * t271 - t12;
t9 = qJ(6) * t271 + qJD(6) * t312 + t11;
t1 = [(-pkin(7) * t325 - t298 / 0.2e1 + t385 + (-0.2e1 * t482 - 0.3e1 / 0.2e1 * t440 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t353) * qJD(1) + (t257 + qJD(1) * (mrSges(4,1) * t312 + mrSges(4,2) * t314) + 0.2e1 * t453) * pkin(2)) * t401 + t562 * t96 + (t390 + t501 + t504) * t312 / 0.2e1 + t404 * t206 + m(6) * (t105 * t226 + t11 * t37 + t12 * t36 + t128 * t189 + t505 * t7 + t78 * t8) + t505 * t47 - t312 * t449 + (t550 + t551) * t95 + (t539 / 0.2e1 + Ifges(7,6) * t475 + t518 / 0.2e1 - t516 / 0.2e1 + Ifges(6,6) * t476 - t264 * mrSges(4,3) - Ifges(4,4) * t458 - Ifges(4,2) * t460 + Ifges(5,5) * t465 - t231 / 0.2e1 + t522 * t472 + t521 * t463 - t552) * t271 - (mrSges(4,2) * t340 - mrSges(4,3) * t502 - Ifges(4,4) * t312) * t252 - t502 * t108 - (-t263 * mrSges(4,3) + Ifges(4,1) * t458 + Ifges(4,4) * t460 + t232 / 0.2e1 + t543 - t558) * t270 + ((-qJD(4) * t500 - t347 * t42 - t351 * t43) * t314 + t146 * t420 + t147 * t421) * mrSges(5,3) - t364 * (-Ifges(5,4) * t420 + Ifges(5,2) * t421) / 0.2e1 + t234 * (-mrSges(5,1) * t421 - mrSges(5,2) * t420) + m(7) * (t10 * t34 + t120 * t14 + t19 * t73 + t2 * t71 + t35 * t9 + t4 * t72) + t296 * (-Ifges(5,5) * t420 + Ifges(5,6) * t421) / 0.2e1 + (-Ifges(5,1) * t420 + Ifges(5,4) * t421) * t465 + m(5) * (t146 * t70 + t147 * t69 + t190 * t43 + t191 * t42 + t206 * t234 - t422) + m(4) * (t177 * t281 + t205 * t264 - t206 * t263 - t422) - t548 * t241 / 0.2e1 + t19 * t118 + t120 * t24 + ((Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t312 - t281 * mrSges(4,3) - Ifges(4,4) * t314 + t340 * mrSges(4,1)) * t253 + (t299 / 0.2e1 - pkin(7) * t324 + t386 + (-0.2e1 * t481 + 0.3e1 / 0.2e1 * Ifges(3,4) * t353) * qJD(1)) * t353 * qJD(2) + (t83 * t454 - Ifges(4,1) * t252 + t82 * t455 + t375 * t478 + t373 * t477 + (mrSges(4,3) + t376) * t178 + (t175 * t455 - t351 * t174 / 0.2e1 + t234 * t377 + t374 * t466 + t370 * t462 + t372 * t528) * qJD(4)) * t314 + t553 * t240 + t78 * t45 + t71 * t44 + t72 * t46 + t312 * t448 + t128 * t119 + (-t241 * t522 + t312 * t521 + t314 * t371) * t469 + t9 * t159 + t11 * t160 + t12 * t161 + (-t241 * t524 + t312 * t522) * t490 + t10 * t162 + t190 * t129 + t191 * t130 + (-t105 * t241 - t312 * t7) * mrSges(6,2) + (t14 * t241 + t2 * t312) * mrSges(7,3) + t69 * t217 + t70 * t218 - t177 * t312 * mrSges(4,3) + t226 * t25 + (-Ifges(7,5) * t241 + Ifges(7,6) * t312) * t488 + (-Ifges(6,4) * t241 + Ifges(6,6) * t312) * t489 + t8 * (mrSges(6,1) * t312 + mrSges(6,3) * t241) + t4 * (-mrSges(7,1) * t312 - mrSges(7,2) * t241) + t205 * t283; ((-t299 / 0.2e1 - t341 / 0.2e1 + t386 + qJD(1) * t481 + (t324 - t427) * pkin(7)) * t353 + (t385 + t298 / 0.2e1 + (t482 + t440 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t353) * qJD(1) + (t325 + t426) * pkin(7) + (-t257 - t453) * pkin(2)) * t349) * qJD(1) + t506 * t118 + t499 * (-pkin(3) - t451) + 0.2e1 * ((t177 * t348 - t178 * t352) * t493 + (m(5) * (t234 * t348 + t352 * t500) / 0.2e1 + (-t263 * t348 + t264 * t352) * t493) * qJD(3)) * pkin(2) + t503 * t119 - t404 * t265 - t405 * t141 + t406 * t140 + t354 + (t14 * t237 + t2 * t255 - t363 * t4 + t506 * t73 + (-t50 + t140) * t35 + (-t51 + t141) * t34) * m(7) + ((t252 * t352 - t253 * t348) * mrSges(4,3) + (t404 * t348 + (t217 * t351 - t218 * t347 + t283) * t352) * qJD(3)) * pkin(2) + (t105 * t323 + t363 * t8 + t255 * t7 + (-t59 + t140) * t37 + (-t58 - t141) * t36 + t503 * t189) * m(6) - t444 * t363 + t445 * t255 - m(5) * (t146 * t167 + t147 * t168 + t234 * t265) - m(4) * (-t263 * t265 + t264 * t266) + ((-t217 * t347 - t218 * t351) * qJD(4) + t367 + t355) * t336 + mrSges(5,3) * t498 + t264 * t432 - t50 * t159 - t59 * t160 - t58 * t161 - t51 * t162 - t168 * t217 - t167 * t218 + t237 * t24 - t266 * t283 + t323 * t25; t354 + ((-mrSges(5,3) * t146 - pkin(9) * t218) * t351 + (-mrSges(5,3) * t147 + pkin(4) * t119 - pkin(9) * t217) * t347) * qJD(4) + t510 * t118 + (-t404 + t432) * t264 + t513 * t162 + t512 * t161 + t367 * pkin(9) + t511 * t160 - t444 * t362 + t445 * t280 + t514 * t159 - m(5) * (t146 * t171 + t147 * t172 + t234 * t264) - t210 * t119 - t172 * t217 - t171 * t218 + t259 * t24 - t263 * t283 + pkin(9) * t355 + t339 * t25 - t499 * pkin(3) + (t14 * t259 + t2 * t280 + t34 * t513 + t35 * t514 - t362 * t4 + t510 * t73) * m(7) + (t105 * t339 + t362 * t8 + t280 * t7 + t511 * t37 + t512 * t36 + (-t210 + t342) * t189) * m(6); t390 + (-t146 * t364 + t147 * t277) * mrSges(5,3) - t234 * (t277 * mrSges(5,1) - mrSges(5,2) * t364) + t405 * t48 + t356 + t515 * t159 + (-t277 * t119 + t346 * t47 + t350 * t45 + (t160 * t350 - t346 * t405) * qJD(5) + m(7) * t34 * t398 + (0.2e1 * t189 * t466 + t346 * t7 + t350 * t8 - t36 * t398 + t37 * t397) * m(6)) * pkin(4) + (t2 * t333 + t337 * t4 - t34 * t48 + t35 * t515 - t73 * t97) * m(7) - m(6) * (-t36 * t48 + t37 * t49) + t448 - t449 - t97 * t118 + (-Ifges(5,5) * t364 - Ifges(5,6) * t277) * t462 + t560 * t357 - t49 * t160 + t174 * t465 + (-Ifges(5,1) * t364 - t438) * t466 - t146 * t217 + t147 * t218 + t333 * t44 + t337 * t46 + (-Ifges(5,2) * t277 + t175 - t276) * t528 + (t550 + t561) * t200; t356 + (t405 + t441) * t37 + (-t406 - t442) * t36 - t117 * t118 + qJ(6) * t44 - pkin(5) * t46 + (t200 * t34 + t35 * t357) * mrSges(7,2) + qJD(6) * t159 - t73 * (mrSges(7,1) * t357 + mrSges(7,3) * t200) + (Ifges(7,3) * t357 - t434) * t476 + t102 * t472 - t189 * (mrSges(6,1) * t357 - mrSges(6,2) * t200) + (-t200 * t522 + t357 * t520) * t464 + (-pkin(5) * t4 + qJ(6) * t2 - t117 * t73 - t34 * t37 + t35 * t508) * m(7) + (-Ifges(6,2) * t357 - t198 + t507) * t475 + (-t200 * t524 + t197 - t435 + t99) * t473; t357 * t118 - t288 * t159 + 0.2e1 * (t4 / 0.2e1 + t73 * t472 + t35 * t464) * m(7) + t46;];
tauc  = t1(:);
