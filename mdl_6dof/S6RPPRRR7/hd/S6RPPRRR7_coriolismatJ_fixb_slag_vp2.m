% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:12
% EndTime: 2019-03-09 02:33:26
% DurationCPUTime: 8.41s
% Computational Cost: add. (21200->442), mult. (37738->586), div. (0->0), fcn. (44475->8), ass. (0->271)
t275 = sin(qJ(6));
t273 = sin(pkin(10));
t274 = cos(pkin(10));
t450 = sin(qJ(4));
t452 = cos(qJ(4));
t244 = -t273 * t452 - t274 * t450;
t276 = sin(qJ(5));
t259 = t450 * t273;
t353 = t452 * t274;
t329 = t353 - t259;
t451 = cos(qJ(5));
t479 = t244 * t276 + t451 * t329;
t386 = t479 * t275;
t488 = -t451 * t244 + t276 * t329;
t147 = -mrSges(7,2) * t488 - mrSges(7,3) * t386;
t277 = cos(qJ(6));
t372 = t277 * t147;
t385 = t479 * t277;
t150 = mrSges(7,1) * t488 - mrSges(7,3) * t385;
t377 = t275 * t150;
t309 = t377 / 0.2e1 - t372 / 0.2e1;
t258 = t273 * pkin(3) + qJ(2);
t230 = -pkin(4) * t244 + t258;
t446 = pkin(5) * t488;
t111 = -pkin(9) * t479 + t230 + t446;
t444 = pkin(1) + qJ(3);
t365 = pkin(7) + t444;
t317 = t365 * t450;
t318 = t365 * t452;
t227 = -t273 * t318 - t274 * t317;
t294 = t244 * pkin(8) + t227;
t226 = t273 * t317 - t274 * t318;
t480 = -pkin(8) * t329 + t226;
t499 = t276 * t480 + t294 * t451;
t47 = t111 * t277 - t275 * t499;
t48 = t111 * t275 + t277 * t499;
t323 = t275 * t47 - t277 * t48;
t472 = m(7) / 0.2e1;
t552 = -(t499 + t323) * t472 - t309;
t414 = t277 * mrSges(7,1);
t420 = t275 * mrSges(7,2);
t328 = t414 - t420;
t523 = t499 * t328;
t525 = t499 * mrSges(6,1);
t109 = -t276 * t294 + t451 * t480;
t539 = t109 * mrSges(6,2);
t268 = Ifges(7,4) * t277;
t489 = -Ifges(7,2) * t275 + t268;
t529 = Ifges(7,6) * t479 - t488 * t489;
t546 = t277 * t529;
t441 = Ifges(7,4) * t275;
t255 = Ifges(7,1) * t277 - t441;
t530 = Ifges(7,5) * t479 - t255 * t488;
t547 = t275 * t530;
t550 = -t523 / 0.2e1 - t525 / 0.2e1 - t539 / 0.2e1 + t546 / 0.4e1 + t547 / 0.4e1;
t250 = Ifges(7,5) * t275 + Ifges(7,6) * t277;
t503 = t479 * t250;
t510 = Ifges(6,6) * t479;
t549 = t503 / 0.2e1 - t510 - t523 - t525 - t539 + t546 / 0.2e1 + t547 / 0.2e1;
t515 = -t479 / 0.2e1;
t548 = 0.2e1 * t515;
t413 = t277 * mrSges(7,2);
t421 = t275 * mrSges(7,1);
t249 = t413 + t421;
t535 = t249 * t488;
t544 = -m(7) * t499 + t535;
t505 = t249 * t479;
t542 = t109 * t535 + t499 * t505;
t453 = t277 / 0.2e1;
t456 = -t275 / 0.2e1;
t460 = t488 / 0.2e1;
t267 = Ifges(7,5) * t277;
t436 = Ifges(7,6) * t275;
t490 = t267 - t436;
t514 = t479 / 0.2e1;
t541 = t453 * t530 + t456 * t529 + t230 * mrSges(6,1) + t490 * t514 + Ifges(6,4) * t548 + (Ifges(6,2) + Ifges(7,3)) * t460;
t540 = -t488 / 0.2e1;
t449 = m(6) * t230;
t428 = mrSges(6,3) * t488;
t538 = t109 * t275;
t537 = t109 * t276;
t536 = t109 * t277;
t401 = t109 * t499;
t492 = (t267 / 0.2e1 - t436 / 0.2e1) * t488;
t312 = t421 / 0.2e1 + t413 / 0.2e1;
t306 = t312 * t479;
t521 = -t505 / 0.2e1 - t306;
t533 = t521 * qJD(6);
t447 = pkin(5) * t479;
t155 = pkin(9) * t488 + t447;
t254 = Ifges(7,1) * t275 + t268;
t369 = t277 * t254;
t252 = Ifges(7,2) * t277 + t441;
t375 = t275 * t252;
t500 = t375 / 0.4e1 - t369 / 0.4e1;
t528 = t500 * t488 + t503 / 0.4e1 - t510 / 0.2e1;
t527 = Ifges(7,6) / 0.2e1;
t464 = t505 / 0.2e1;
t524 = t472 * t499;
t507 = t488 * mrSges(6,1);
t508 = t479 * mrSges(6,2);
t341 = t507 + t508;
t271 = t275 ^ 2;
t272 = t277 ^ 2;
t366 = t271 + t272;
t504 = t366 * t479;
t518 = pkin(9) * t504 - t446;
t422 = t272 * mrSges(7,3);
t423 = t271 * mrSges(7,3);
t487 = (-t422 / 0.2e1 - t423 / 0.2e1) * t479;
t502 = t488 * t328;
t517 = t487 + t502 / 0.2e1 + t507 / 0.2e1 + t508 / 0.2e1;
t516 = -t249 / 0.2e1;
t513 = mrSges(7,1) * t479;
t512 = mrSges(7,2) * t479;
t430 = t479 * mrSges(6,3);
t392 = t479 * t488;
t501 = -t255 / 0.4e1 + t252 / 0.4e1;
t491 = t369 / 0.2e1 - t375 / 0.2e1;
t330 = mrSges(7,3) * (t271 / 0.2e1 + t272 / 0.2e1);
t367 = t273 ^ 2 + t274 ^ 2;
t371 = t277 * t150;
t380 = t275 * t147;
t485 = -t371 / 0.2e1 - t380 / 0.2e1;
t419 = t275 * mrSges(7,3);
t146 = t419 * t488 - t512;
t373 = t277 * t146;
t427 = t488 * mrSges(7,3);
t149 = t277 * t427 + t513;
t378 = t275 * t149;
t484 = t373 / 0.2e1 - t378 / 0.2e1;
t376 = t275 * t488;
t145 = mrSges(7,3) * t376 - t512;
t374 = t277 * t145;
t370 = t277 * t488;
t148 = mrSges(7,3) * t370 + t513;
t379 = t275 * t148;
t483 = t374 / 0.2e1 - t379 / 0.2e1;
t314 = t329 * pkin(4);
t112 = t155 + t314;
t54 = t112 * t275 + t536;
t408 = t54 * t277;
t53 = t112 * t277 - t538;
t409 = t53 * t275;
t322 = t408 - t409;
t97 = Ifges(7,5) * t488 + t255 * t479;
t410 = t277 * t97;
t94 = Ifges(7,6) * t488 + t479 * t489;
t417 = t275 * t94;
t482 = Ifges(6,1) * t514 - Ifges(6,4) * t488 - t417 / 0.2e1 + t410 / 0.2e1;
t481 = -t488 * t490 / 0.4e1 - t109 * t516;
t471 = -pkin(5) / 0.2e1;
t478 = -pkin(5) * t524 - t535 * t471 + (t322 * t472 + t484) * pkin(9) + (t408 / 0.2e1 - t409 / 0.2e1) * mrSges(7,3) + t550;
t477 = -m(4) / 0.2e1;
t476 = -m(5) / 0.2e1;
t475 = -m(6) / 0.2e1;
t473 = -m(7) / 0.2e1;
t470 = m(6) * pkin(4);
t469 = m(7) * pkin(4);
t468 = -mrSges(7,1) / 0.2e1;
t467 = -mrSges(7,2) / 0.2e1;
t6 = t505 * t460 + (t428 / 0.2e1 + t535 / 0.2e1 + t552) * t479 + (-t430 / 0.2e1 + t484 + (-t109 + t322) * t472) * t488;
t61 = t155 * t275 + t536;
t406 = t61 * t277;
t60 = t155 * t277 - t538;
t407 = t60 * t275;
t321 = t406 - t407;
t465 = -t535 / 0.2e1;
t9 = ((-t109 + t321) * t472 + t464 + t483) * t488 - (t465 - t552) * t479;
t466 = t6 * qJD(4) + t9 * qJD(5);
t448 = pkin(4) * t276;
t262 = pkin(9) + t448;
t457 = t262 / 0.2e1;
t455 = t275 / 0.2e1;
t454 = -t277 / 0.2e1;
t445 = pkin(5) * t249;
t440 = Ifges(6,5) * t488;
t431 = t479 * mrSges(6,1);
t424 = t258 * mrSges(5,1);
t418 = t275 * t48;
t139 = t328 * t479;
t143 = t252 * t479;
t144 = t254 * t479;
t7 = -t109 * t139 + t47 * t147 - t48 * t150 + (t323 * mrSges(7,3) - t144 * t453 + t250 * t540 + t94 * t454 + (-t143 + t97) * t456) * t479;
t405 = t7 * qJD(1);
t295 = t244 ^ 2 + t329 ^ 2;
t400 = t109 * t479;
t14 = (t505 + t430) * t479 + t295 * mrSges(5,3) - (t372 - t377 - t428) * t488 + m(7) * (t323 * t488 - t400) + m(6) * (-t488 * t499 - t400) + m(5) * (-t226 * t329 + t227 * t244) + (m(4) * t444 + mrSges(4,3)) * t367;
t404 = qJD(1) * t14;
t298 = t244 * mrSges(5,1) - mrSges(5,2) * t329 - t341;
t26 = t371 + t380 + t274 * mrSges(4,2) + t273 * mrSges(4,1) + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(7) * (t277 * t47 + t418) + t449 + m(5) * t258 - t298;
t403 = qJD(1) * t26;
t218 = t488 * mrSges(6,2);
t236 = t244 * mrSges(5,2);
t363 = t451 * pkin(4);
t263 = -t363 - pkin(5);
t337 = t366 * t262;
t352 = t328 * t514;
t285 = (t263 * t479 - t337 * t488) * t472 - t352 + (-t276 * t488 - t451 * t479) * t470 / 0.2e1 - t366 * t427 / 0.2e1;
t288 = (t275 * t54 + t277 * t53) * t472 + t146 * t455 + t149 * t453 - t314 * t475;
t17 = -mrSges(5,1) * t329 + t218 - t236 + t285 - t288 - t431;
t397 = t17 * qJD(1);
t296 = -t330 * t479 + t485;
t292 = t139 * t515 + t296 * t488;
t311 = -t420 / 0.2e1 + t414 / 0.2e1;
t19 = t292 - t311;
t396 = t19 * qJD(1);
t339 = t366 * t488;
t290 = -t488 * t330 + t218 / 0.2e1 + (-pkin(9) * t339 - t447) * t472 - t352;
t293 = (t275 * t61 + t277 * t60) * t473 + mrSges(6,2) * t460 + t145 * t456 + t148 * t454;
t21 = mrSges(6,1) * t548 + t290 + t293;
t395 = t21 * qJD(1);
t393 = t479 ^ 2;
t383 = t263 * t139;
t382 = t263 * t249;
t307 = t312 * t488;
t27 = t307 + t309;
t381 = t27 * qJD(1);
t283 = (-t339 * t488 - t393) * t472 + (-t488 ^ 2 - t393) * m(6) / 0.2e1 + t295 * t476 + t367 * t477;
t302 = t366 * t473 + t475 + t476 + t477;
t34 = t283 + t302;
t368 = t34 * qJD(1);
t364 = mrSges(7,3) * t406;
t362 = Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t358 = -t419 / 0.2e1;
t355 = t275 * t451;
t354 = t277 * t451;
t333 = t255 * t455 + t453 * t489 + t491;
t332 = t451 * t467;
t331 = -t355 / 0.2e1;
t1 = t48 * t146 + t54 * t147 + t47 * t149 + t53 * t150 - t230 * t218 + t258 * t236 - (t460 * t490 + t482) * t488 + pkin(4) * t353 * t449 + m(7) * (t47 * t53 + t48 * t54 - t401) + (-t424 + (-t341 - t449) * pkin(4) + (0.2e1 * t353 - t259) * Ifges(5,4)) * t259 + (Ifges(5,4) * t244 + (Ifges(5,1) - Ifges(5,2)) * t329) * t244 + (-t353 * Ifges(5,4) + pkin(4) * t341 + t424) * t353 + (Ifges(6,1) * t540 + t362 * t488 + t541) * t479 + t542;
t325 = t1 * qJD(1) + t6 * qJD(2);
t4 = m(7) * (t47 * t60 + t48 * t61 - t401) + t60 * t150 + t47 * t148 + t61 * t147 + t48 * t145 + t541 * t479 + (-t230 * mrSges(6,2) - t492 + (-Ifges(6,1) / 0.2e1 + t362) * t479 - t482) * t488 + t542;
t324 = t4 * qJD(1) + t9 * qJD(2);
t36 = m(7) * (t488 * t504 - t392);
t320 = t6 * qJD(1) + t36 * qJD(2);
t41 = m(7) * (0.1e1 - t366) * t392;
t319 = t9 * qJD(1) - t41 * qJD(2);
t308 = t366 * t451;
t289 = (-mrSges(6,1) - t328) * t448 + (mrSges(7,3) * t366 - mrSges(6,2)) * t363;
t126 = (t262 * t308 + t263 * t276) * t469 + t289;
t282 = (t263 * t488 + t262 * t504 + (-t276 * t479 + t308 * t488) * pkin(4)) * t472 - t517;
t287 = t473 * t518 + t517;
t25 = t282 + t287;
t281 = Ifges(6,5) * t540 + t448 * t464 + t60 * t358 + t364 / 0.2e1 + (t465 + t524) * t263 + (t321 * t472 + t483) * t262 + ((t354 * t48 - t355 * t47 - t537) * t472 + t150 * t331 + t147 * t354 / 0.2e1) * pkin(4) + t528 + t550;
t3 = -t281 - t440 / 0.2e1 + t478 + t528;
t304 = -t3 * qJD(1) + t25 * qJD(2) + t126 * qJD(4);
t301 = Ifges(7,3) * t514 + t53 * mrSges(7,1) / 0.2e1 + t54 * t467;
t10 = -t383 / 0.2e1 + (Ifges(7,5) * t540 - t97 / 0.4e1 + t143 / 0.4e1 + t150 * t457) * t277 + (t488 * t527 + t144 / 0.4e1 + t94 / 0.4e1 + t147 * t457) * t275 + (t501 * t277 + (t254 / 0.4e1 + t489 / 0.4e1) * t275 + t262 * t330) * t479 + t301 + t481;
t160 = t333 + t382;
t65 = t464 - t306;
t303 = -t10 * qJD(1) - t65 * qJD(2) + t160 * qJD(4);
t300 = Ifges(7,3) * t515 + t60 * t468 + t61 * mrSges(7,2) / 0.2e1;
t105 = (-t263 / 0.2e1 + pkin(5) / 0.2e1) * t249 + (pkin(4) * t332 - t254 / 0.2e1 - t489 / 0.2e1) * t277 + (t363 * t468 - t255 / 0.2e1 + t252 / 0.2e1) * t275;
t297 = t48 * t358 - t275 * t144 / 0.4e1 - t277 * t143 / 0.4e1 + mrSges(7,3) * t418 / 0.2e1 - t417 / 0.4e1 + t410 / 0.4e1 - t481 - t501 * t385 - (t489 + t254) * t386 / 0.4e1;
t286 = pkin(9) * t296 + t139 * t471 + t297;
t13 = t286 + t300 + t492;
t161 = t333 - t445;
t64 = (t516 + t312) * t479;
t299 = t13 * qJD(1) + t64 * qJD(2) - t105 * qJD(4) + t161 * qJD(5);
t106 = t382 / 0.2e1 - t445 / 0.2e1 + (mrSges(7,1) * t331 + t277 * t332) * pkin(4) + t333;
t33 = t283 - t302;
t28 = t307 - t309;
t24 = t282 - t287;
t23 = t285 + t288;
t22 = mrSges(6,1) * t514 - t431 / 0.2e1 + t290 - t293;
t20 = t292 + t311;
t12 = t286 - Ifges(7,5) * t370 / 0.2e1 + t376 * t527 - t300;
t11 = -t492 + t383 / 0.2e1 + t297 + t301 + (t485 + t487) * t262;
t2 = t281 - (Ifges(6,5) / 0.2e1 - t500) * t488 + (t250 / 0.4e1 - Ifges(6,6) / 0.2e1) * t479 + t478;
t5 = [qJD(2) * t26 + qJD(3) * t14 + qJD(4) * t1 + qJD(5) * t4 + qJD(6) * t7, qJD(3) * t33 + qJD(6) * t20 + t403 + t466, qJD(2) * t33 + qJD(4) * t23 + qJD(5) * t22 + qJD(6) * t28 + t404, t23 * qJD(3) + (-t430 * t448 + t363 * t428 + (-t451 * t499 + t537) * t470 - t226 * mrSges(5,2) - t227 * mrSges(5,1) - t440 + Ifges(5,5) * t244 - Ifges(5,6) * t329 - t544 * t263 + (m(7) * t322 + t373 - t378) * t262 - t491 * t488 + t322 * mrSges(7,3) + t549) * qJD(4) + t2 * qJD(5) + t11 * qJD(6) + t325, t22 * qJD(3) + t2 * qJD(4) + t12 * qJD(6) + t324 + (-mrSges(7,3) * t407 + t364 + (-Ifges(6,5) - t491) * t488 + t544 * pkin(5) + (m(7) * t321 + t374 - t379) * pkin(9) + t549) * qJD(5), t405 + t20 * qJD(2) + t28 * qJD(3) + t11 * qJD(4) + t12 * qJD(5) + (-mrSges(7,1) * t48 - mrSges(7,2) * t47 - t503) * qJD(6); qJD(3) * t34 + qJD(6) * t19 - t403 + t466, qJD(4) * t36 - qJD(5) * t41, t368 (t298 - t502 + (m(7) * t263 - t451 * t470) * t488 + (m(7) * t337 + t276 * t470 + t422 + t423) * t479) * qJD(4) + t24 * qJD(5) + t533 + t320, t24 * qJD(4) + (m(7) * t518 + mrSges(7,3) * t504 - t341 - t502) * qJD(5) + t533 + t319, -qJD(6) * t502 + t396 + (qJD(4) + qJD(5)) * t521; -qJD(2) * t34 - qJD(4) * t17 - qJD(5) * t21 - qJD(6) * t27 - t404, -t368, 0, -t397, -t395, -qJD(6) * t249 - t381; qJD(3) * t17 - qJD(5) * t3 - qJD(6) * t10 - t325, qJD(5) * t25 - qJD(6) * t65 - t320, t397, qJD(5) * t126 + qJD(6) * t160 ((-pkin(5) * t276 + pkin(9) * t308) * t469 + t289) * qJD(5) + t106 * qJD(6) + t304, t106 * qJD(5) + (-t262 * t328 + t490) * qJD(6) + t303; qJD(3) * t21 + qJD(4) * t3 + qJD(6) * t13 - t324, -qJD(4) * t25 + qJD(6) * t64 - t319, t395, -qJD(6) * t105 - t304, t161 * qJD(6) (-pkin(9) * t328 + t490) * qJD(6) + t299; -qJD(2) * t19 + qJD(3) * t27 + qJD(4) * t10 - qJD(5) * t13 - t405, qJD(4) * t65 - qJD(5) * t64 - t396, t381, qJD(5) * t105 - t303, -t299, 0;];
Cq  = t5;
