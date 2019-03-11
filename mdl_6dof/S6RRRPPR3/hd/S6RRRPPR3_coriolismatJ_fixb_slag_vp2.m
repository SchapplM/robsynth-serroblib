% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:28
% EndTime: 2019-03-09 15:28:40
% DurationCPUTime: 5.95s
% Computational Cost: add. (11975->467), mult. (22616->592), div. (0->0), fcn. (23084->6), ass. (0->264)
t284 = sin(qJ(6));
t458 = -t284 / 0.2e1;
t287 = cos(qJ(6));
t457 = t287 / 0.2e1;
t285 = sin(qJ(3));
t286 = sin(qJ(2));
t390 = t285 * t286;
t455 = cos(qJ(3));
t456 = cos(qJ(2));
t249 = -t455 * t456 + t390;
t361 = t455 * t286;
t250 = t285 * t456 + t361;
t272 = -t456 * pkin(2) - pkin(1);
t346 = -t249 * pkin(3) + t250 * qJ(4) - t272;
t88 = -pkin(4) * t249 + t346;
t499 = m(6) * t88 + mrSges(6,1) * t250 + mrSges(6,2) * t249;
t446 = Ifges(7,2) * t284;
t449 = Ifges(7,4) * t287;
t333 = -t446 + t449;
t450 = Ifges(7,4) * t284;
t334 = Ifges(7,1) * t287 - t450;
t257 = -Ifges(7,2) * t287 - t450;
t258 = -Ifges(7,1) * t284 - t449;
t502 = t257 * t458 + t258 * t457;
t514 = t284 / 0.2e1;
t340 = t333 * t457 + t334 * t514 - t502;
t375 = t456 * pkin(7);
t261 = t456 * pkin(8) + t375;
t465 = -pkin(8) - pkin(7);
t184 = t261 * t285 - t465 * t361;
t125 = -t250 * qJ(5) + t184;
t277 = t287 * mrSges(7,1);
t424 = t284 * mrSges(7,2);
t336 = t277 - t424;
t507 = t125 * t336;
t491 = t455 * t261 + t465 * t390;
t509 = t491 * mrSges(5,1);
t510 = t491 * mrSges(4,1);
t511 = t184 * mrSges(5,3);
t512 = t184 * mrSges(4,2);
t513 = t125 * mrSges(6,1);
t501 = t249 * qJ(5) + t491;
t525 = t501 * mrSges(6,2);
t529 = t512 - t507 - t509 - t510 - t511 - t513 + t525;
t528 = t507 / 0.2e1 + t509 / 0.2e1 + t510 / 0.2e1 + t511 / 0.2e1 - t512 / 0.2e1 + t513 / 0.2e1 - t525 / 0.2e1;
t281 = t284 ^ 2;
t282 = t287 ^ 2;
t380 = t281 + t282;
t526 = t380 * mrSges(7,3);
t411 = t125 * t501;
t524 = t284 * t501;
t523 = t287 * t501;
t426 = t250 * mrSges(6,3);
t522 = t250 * mrSges(5,2) / 0.2e1 - t426 / 0.2e1;
t500 = -m(5) * t346 + mrSges(5,1) * t249 - mrSges(5,3) * t250;
t444 = Ifges(7,6) * t284;
t447 = Ifges(7,5) * t287;
t332 = t444 - t447;
t469 = Ifges(7,1) / 0.2e1;
t498 = Ifges(4,4) + Ifges(6,4) - Ifges(5,5);
t521 = -mrSges(6,3) * t501 + (t282 * t469 - Ifges(5,1) + Ifges(5,3) - Ifges(4,1) + Ifges(6,1) + Ifges(4,2) - Ifges(6,2) - Ifges(7,3) + (-t449 + t446 / 0.2e1) * t284) * t249 + (-t332 - t498) * t250;
t374 = t455 * pkin(2);
t271 = -t374 - pkin(3);
t266 = -pkin(4) + t271;
t452 = t285 * pkin(2);
t269 = qJ(4) + t452;
t520 = -t269 * t125 + t266 * t501;
t288 = -pkin(3) - pkin(4);
t519 = -qJ(4) * t125 + t288 * t501;
t423 = t284 * mrSges(7,3);
t153 = mrSges(7,2) * t249 - t250 * t423;
t388 = t287 * t153;
t144 = t388 / 0.2e1;
t428 = t249 * mrSges(5,2);
t225 = -t428 / 0.2e1;
t234 = t249 * mrSges(6,3);
t420 = t287 * mrSges(7,2);
t425 = t284 * mrSges(7,1);
t321 = t420 / 0.2e1 + t425 / 0.2e1;
t451 = mrSges(7,3) * t250;
t155 = -mrSges(7,1) * t249 - t287 * t451;
t394 = t284 * t155;
t472 = -mrSges(5,2) / 0.2e1;
t475 = m(7) / 0.4e1;
t478 = m(6) / 0.4e1;
t479 = m(6) / 0.2e1;
t518 = t501 * t479 - t394 / 0.2e1 + t144 + t225 + t234 + (t472 + t321) * t249 + 0.2e1 * (t475 + t478) * t501;
t358 = t394 / 0.2e1;
t464 = t501 / 0.4e1;
t517 = t321 * t249 + t358 - t388 / 0.2e1 + 0.2e1 * (t464 - t501 / 0.4e1) * m(6);
t283 = qJ(4) + pkin(5);
t506 = t125 * t283;
t265 = pkin(5) + t269;
t505 = t265 * t125;
t503 = t269 * t184;
t326 = -mrSges(6,1) - mrSges(5,3) - t336;
t497 = -pkin(3) * t491 - qJ(4) * t184;
t335 = t420 + t425;
t148 = t335 * t249;
t149 = t335 * t250;
t322 = -t444 / 0.2e1 + t447 / 0.2e1;
t466 = -pkin(4) - pkin(9);
t58 = pkin(5) * t250 + t466 * t249 + t346;
t45 = -t125 * t284 + t287 * t58;
t46 = t125 * t287 + t284 * t58;
t86 = -Ifges(7,6) * t249 + t333 * t250;
t87 = -Ifges(7,5) * t249 + t334 * t250;
t496 = -t125 * t148 + (t87 * t457 + t86 * t458 + (-t322 + t498) * t249) * t249 - t346 * (mrSges(5,1) * t250 + mrSges(5,3) * t249) + t46 * t153 + t45 * t155 + t272 * (mrSges(4,1) * t250 - mrSges(4,2) * t249) + t88 * (-t249 * mrSges(6,1) + t250 * mrSges(6,2)) + (t149 + t426) * t501;
t353 = -t282 / 0.2e1 - t281 / 0.2e1;
t337 = t353 * mrSges(7,3);
t492 = t148 + t234;
t230 = t249 * qJ(4);
t160 = t250 * pkin(3) + t230;
t115 = -pkin(4) * t250 - t160;
t476 = m(7) / 0.2e1;
t61 = -pkin(5) * t249 + t466 * t250 - t160;
t51 = t287 * t61 - t524;
t52 = t284 * t61 + t523;
t490 = (t284 * t52 + t287 * t51) * t476 + t115 * t479;
t453 = pkin(2) * t286;
t59 = t61 - t453;
t47 = t287 * t59 - t524;
t48 = t284 * t59 + t523;
t90 = t115 - t453;
t489 = (t284 * t48 + t287 * t47) * t476 + t90 * t479;
t463 = -t249 / 0.2e1;
t488 = t463 * t526;
t487 = t478 + m(5) / 0.4e1;
t403 = t249 * t284;
t152 = -mrSges(7,2) * t250 - mrSges(7,3) * t403;
t389 = t287 * t152;
t402 = t249 * t287;
t154 = t250 * mrSges(7,1) - mrSges(7,3) * t402;
t395 = t284 * t154;
t318 = t395 / 0.2e1 - t389 / 0.2e1;
t484 = 0.2e1 * m(7);
t483 = 0.2e1 * t250;
t482 = m(5) / 0.2e1;
t480 = -m(6) / 0.2e1;
t477 = -m(7) / 0.2e1;
t474 = m(7) * pkin(2);
t473 = mrSges(7,1) / 0.2e1;
t471 = -mrSges(7,2) / 0.2e1;
t470 = mrSges(7,2) / 0.2e1;
t468 = -Ifges(7,3) / 0.2e1;
t467 = m(5) + m(6);
t264 = -pkin(9) + t266;
t462 = t264 / 0.2e1;
t461 = -t265 / 0.2e1;
t280 = -pkin(3) + t466;
t460 = t280 / 0.2e1;
t459 = -t283 / 0.2e1;
t454 = m(5) * t491;
t448 = Ifges(7,5) * t250;
t445 = Ifges(7,6) * t250;
t1 = (mrSges(4,2) * t453 + t521) * t250 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t286 + (m(4) * t272 + mrSges(4,1) * t249) * pkin(2)) * t286 + t48 * t152 + t47 * t154 + m(7) * (t45 * t47 + t46 * t48 - t411) + ((Ifges(3,1) - Ifges(3,2)) * t286 - pkin(1) * mrSges(3,2) + Ifges(3,4) * t456) * t456 + t496 + t500 * (t160 + t453) + t499 * t90;
t443 = t1 * qJD(1);
t2 = t521 * t250 + m(7) * (t45 * t51 + t46 * t52 - t411) + t52 * t152 + t51 * t154 + t500 * t160 + t499 * t115 + t496;
t430 = t2 * qJD(1);
t422 = t284 * t47;
t421 = t284 * t51;
t419 = t287 * t48;
t418 = t287 * t52;
t417 = t288 * mrSges(6,3);
t5 = t45 * t152 - t46 * t154 + ((mrSges(7,1) * t501 - mrSges(7,3) * t46 - Ifges(7,4) * t402 - t445) * t287 + (-t501 * mrSges(7,2) + t45 * mrSges(7,3) - t448 + (t450 + (-Ifges(7,1) + Ifges(7,2)) * t287) * t249) * t284) * t249;
t416 = t5 * qJD(1);
t331 = t284 * t45 - t287 * t46;
t409 = t501 * t249;
t10 = t492 * t249 + (-t389 + t395 + t426) * t250 + m(7) * (t331 * t250 + t409) + m(6) * (-t125 * t250 + t409);
t415 = qJD(1) * t10;
t387 = t287 * t154;
t396 = t284 * t152;
t18 = (t396 + t387 + m(7) * (t284 * t46 + t287 * t45) + t499 - t500) * t250;
t414 = qJD(1) * t18;
t319 = t153 * t514 + t155 * t457;
t294 = (-mrSges(6,1) - t336 / 0.2e1) * t249 + (mrSges(6,2) + t337) * t250 + t319;
t348 = t380 * t250;
t296 = (t249 * t269 - t250 * t266) * t480 + (t249 * t265 - t264 * t348) * t477;
t12 = t294 + t296 + t489;
t413 = t12 * qJD(1);
t408 = t501 * t335;
t299 = (-t250 * t288 + t230) * t480 + (t249 * t283 - t280 * t348) * t477;
t16 = t294 + t299 + t490;
t407 = t16 * qJD(1);
t311 = (-t424 / 0.2e1 + t277 / 0.2e1) * t250;
t320 = t396 / 0.2e1 + t387 / 0.2e1;
t21 = -t249 * t337 + t311 + t320;
t406 = t21 * qJD(1);
t310 = t321 * t250;
t23 = t310 + t318;
t405 = t23 * qJD(1);
t401 = t265 * t149;
t400 = t265 * t335;
t399 = t269 * t250;
t398 = t283 * t149;
t397 = t283 * t335;
t392 = t284 * t264;
t391 = t284 * t280;
t385 = t287 * t264;
t384 = t287 * t280;
t53 = (-t353 * m(7) + t479) * t483;
t383 = t53 * qJD(1);
t378 = mrSges(7,3) * t419;
t377 = mrSges(7,3) * t418;
t373 = t266 * t234;
t372 = mrSges(5,2) * t399;
t371 = mrSges(6,3) * t399;
t365 = -t423 / 0.2e1;
t364 = t423 / 0.2e1;
t363 = t455 * t501;
t362 = t455 * t269;
t360 = t408 / 0.2e1;
t356 = -t392 / 0.2e1;
t354 = -t384 / 0.2e1;
t352 = m(7) * t380;
t347 = t380 * t285;
t345 = t249 * t374;
t330 = t419 - t422;
t329 = t418 - t421;
t292 = (-mrSges(4,1) - mrSges(5,1) + mrSges(6,2) - t526) * t452 + (-mrSges(4,2) - t326) * t374;
t35 = -(t264 * t347 + t455 * t265) * t474 - t292 + (-m(6) * (t266 * t285 + t362) - m(5) * (t271 * t285 + t362)) * pkin(2);
t289 = t492 * t374 / 0.2e1 + (t184 * t482 - t318 + t522) * t452 + t345 * t472 + (-t503 + (t271 + t374) * t491) * t482 + (-t505 + t329 * t264 + (-t331 * t285 + t363) * pkin(2)) * t476 - t377 / 0.2e1 + t401 / 0.2e1 + t371 / 0.2e1 + t373 / 0.2e1 - t372 / 0.2e1 + t264 * t144 + t155 * t356 + t51 * t364 + t271 * t225 + ((t125 * t285 + t363) * pkin(2) + t520) * t479 - t528;
t290 = t417 * t463 + t378 / 0.2e1 - t398 / 0.2e1 + (t330 * t280 - t506) * t477 - m(5) * t497 / 0.2e1 + t522 * qJ(4) + t153 * t354 + t280 * t358 + t47 * t365 + pkin(3) * t225 + t519 * t480 + t528;
t4 = t290 + t289;
t328 = t4 * qJD(1) - t35 * qJD(2);
t304 = (-Ifges(7,1) / 0.4e1 + Ifges(7,2) / 0.2e1) * t284 + t258 / 0.4e1 - 0.3e1 / 0.2e1 * t449;
t314 = (t469 - Ifges(7,2) / 0.4e1) * t287 + t257 / 0.4e1;
t324 = t47 * t473 + t48 * t471;
t6 = t360 + (t154 * t462 + t448) * t287 + (t468 - t264 * t337 + (mrSges(7,1) * t461 + t314) * t287) * t249 + (-t445 + t152 * t462 + (t265 * t470 + t304) * t249) * t284 + t324;
t62 = t340 - t400;
t327 = -t6 * qJD(1) + t62 * qJD(2);
t116 = m(7) * t265 + 0.4e1 * t487 * t269 - t326;
t15 = (t464 + t422 / 0.4e1 - t419 / 0.4e1) * t484 + t517;
t325 = -qJD(1) * t15 - qJD(2) * t116;
t323 = t52 * t471 + t51 * t473;
t317 = t249 * (-Ifges(7,5) * t284 - Ifges(7,6) * t287);
t309 = t336 * t249 / 0.2e1;
t306 = t321 * t452;
t55 = -(t461 + t459) * t335 - t306 - t340;
t65 = t340 - t397;
t8 = t360 + (t154 * t460 + t448) * t287 + (t468 - t280 * t337 + (mrSges(7,1) * t459 + t314) * t287) * t249 + (-t445 + t152 * t460 + (t283 * t470 + t304) * t249) * t284 + t323;
t308 = t8 * qJD(1) + t55 * qJD(2) - t65 * qJD(3);
t187 = m(7) * t283 + t467 * qJ(4) - t326;
t20 = (t464 + t421 / 0.4e1 - t418 / 0.4e1) * t484 + t517;
t63 = (-pkin(5) + (-0.1e1 / 0.2e1 - t353) * t452) * m(7) + (-m(7) - t467) * qJ(4) + t326;
t307 = qJD(1) * t20 - qJD(2) * t63 + qJD(3) * t187;
t305 = t87 * t458 - t317 / 0.2e1 - t287 * t86 / 0.2e1 + (-Ifges(4,5) - Ifges(5,4)) * t249 + (Ifges(5,6) - Ifges(4,6) + t502) * t250;
t302 = -Ifges(6,5) * t250 - Ifges(6,6) * t249 + t305;
t297 = -t336 * t463 + t319 + t380 * t451 / 0.2e1;
t293 = Ifges(7,3) * t463 + t284 * (t333 * t249 + t445) / 0.4e1 - t287 * (t334 * t249 + t448) / 0.4e1 - t408 / 0.2e1 + (t364 + t365) * t46 + (-t258 / 0.2e1 + t333 / 0.4e1) * t403 - (0.2e1 * t257 + t334) * t402 / 0.4e1 + (t322 + t332 / 0.4e1) * t250;
t273 = qJ(4) * t374;
t64 = 0.4e1 * pkin(5) * t475 + (t479 + t352 / 0.2e1 + t482) * t452 - t326 + (t475 + t487) * (0.4e1 * qJ(4) + 0.2e1 * t452);
t56 = -t400 / 0.2e1 - t397 / 0.2e1 - t306 + t340;
t54 = -t348 * t476 + t250 * t480 + (t352 / 0.4e1 + t478) * t483;
t24 = t310 - t318;
t22 = t311 - t320 + t488;
t19 = t329 * t476 + t454 + t518;
t17 = t297 - t299 + t490;
t14 = t330 * t476 + t491 * t482 + t454 / 0.2e1 + t518;
t13 = -t296 + t297 + t489;
t9 = t293 + t154 * t354 - t152 * t391 / 0.2e1 + t283 * t309 + t323 + t488 * t280;
t7 = t293 - t154 * t385 / 0.2e1 + t152 * t356 + t265 * t309 + t324 + t488 * t264;
t3 = t289 - t290 + t302;
t11 = [qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t18 + qJD(5) * t10 + qJD(6) * t5, t443 + (t302 + m(4) * (-t184 * t285 - t455 * t491) * pkin(2) + mrSges(7,3) * t422 + t153 * t385 + m(5) * (t271 * t491 - t503) + m(7) * (t330 * t264 - t505) - t378 + t401 + t371 + (-t250 * t452 + t345) * mrSges(4,3) + (pkin(7) * mrSges(3,2) - Ifges(3,6)) * t286 + t373 - t372 + Ifges(3,5) * t456 - t271 * t428 - t155 * t392 - mrSges(3,1) * t375 + m(6) * t520 + t529) * qJD(2) + t3 * qJD(3) + t14 * qJD(4) + t13 * qJD(5) + t7 * qJD(6), t3 * qJD(2) + t19 * qJD(4) + t17 * qJD(5) + t9 * qJD(6) + t430 + (t305 - t377 + t398 + (-Ifges(6,5) + (-mrSges(5,2) + mrSges(6,3)) * qJ(4)) * t250 + (pkin(3) * mrSges(5,2) - Ifges(6,6) + t417) * t249 + 0.2e1 * (t280 * t329 - t506) * t476 + 0.2e1 * t497 * t482 + mrSges(7,3) * t421 + t153 * t384 - t155 * t391 + 0.2e1 * t519 * t479 + t529) * qJD(3), qJD(2) * t14 + qJD(3) * t19 + qJD(5) * t54 + qJD(6) * t22 + t414, qJD(2) * t13 + qJD(3) * t17 + qJD(4) * t54 + qJD(6) * t24 + t415, t416 + t7 * qJD(2) + t9 * qJD(3) + t22 * qJD(4) + t24 * qJD(5) + (-mrSges(7,1) * t46 - mrSges(7,2) * t45 + t317) * qJD(6); qJD(3) * t4 + qJD(4) * t15 - qJD(5) * t12 - qJD(6) * t6 - t443, -qJD(3) * t35 + qJD(4) * t116 + qJD(6) * t62 (m(6) * (t288 * t452 + t273) + (t280 * t347 + t455 * t283) * t474 + m(5) * (-pkin(3) * t452 + t273) + t292) * qJD(3) + t64 * qJD(4) + t56 * qJD(6) + t328, qJD(3) * t64 - t325, -t413, t56 * qJD(3) + (-t264 * t336 + t332) * qJD(6) + t327; -qJD(2) * t4 + qJD(4) * t20 - qJD(5) * t16 - qJD(6) * t8 - t430, -qJD(4) * t63 - qJD(6) * t55 - t328, qJD(4) * t187 + qJD(6) * t65, t307, -t407 (-t280 * t336 + t332) * qJD(6) - t308; -qJD(2) * t15 - qJD(3) * t20 - qJD(5) * t53 - qJD(6) * t21 - t414, qJD(3) * t63 + t325, -t307, 0, -t383, -qJD(6) * t336 - t406; qJD(2) * t12 + qJD(3) * t16 + qJD(4) * t53 - qJD(6) * t23 - t415, t413, t407, t383, 0, -qJD(6) * t335 - t405; qJD(2) * t6 + qJD(3) * t8 + qJD(4) * t21 + qJD(5) * t23 - t416, qJD(3) * t55 - t327, t308, t406, t405, 0;];
Cq  = t11;
