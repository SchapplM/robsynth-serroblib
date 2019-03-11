% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:14:09
% EndTime: 2019-03-08 20:14:18
% DurationCPUTime: 4.62s
% Computational Cost: add. (6145->506), mult. (14893->701), div. (0->0), fcn. (13331->8), ass. (0->277)
t424 = m(6) / 0.4e1 + m(7) / 0.4e1;
t532 = 0.2e1 * t424;
t329 = sin(qJ(5));
t317 = t329 * mrSges(7,1);
t332 = cos(qJ(5));
t318 = t332 * mrSges(7,2);
t515 = t318 + t317;
t483 = t515 / 0.2e1;
t460 = Ifges(7,4) * t329;
t283 = Ifges(7,1) * t332 - t460;
t330 = sin(qJ(4));
t333 = cos(qJ(4));
t205 = t330 * Ifges(7,5) + t333 * t283;
t461 = Ifges(6,4) * t329;
t285 = Ifges(6,1) * t332 - t461;
t206 = t330 * Ifges(6,5) + t333 * t285;
t387 = (Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1) * t330;
t531 = t387 + t205 / 0.2e1 + t206 / 0.2e1;
t418 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t277 = t329 * mrSges(6,1) + t332 * mrSges(6,2);
t482 = -t277 / 0.2e1;
t417 = -mrSges(7,2) / 0.2e1 - mrSges(6,2) / 0.2e1;
t519 = t417 * t332;
t530 = -t329 * t418 - t482 + t483 + t519;
t501 = m(7) * pkin(5);
t425 = t501 / 0.2e1;
t366 = t425 + t418;
t529 = -t329 * t366 + t519;
t232 = t515 * t330;
t233 = t277 * t330;
t469 = pkin(9) * t333;
t472 = pkin(4) * t330;
t271 = qJ(3) - t469 + t472;
t249 = t332 * t271;
t335 = -pkin(2) - pkin(8);
t394 = -t329 * t335 + pkin(5);
t429 = t332 * t333;
t411 = qJ(6) * t429;
t102 = t330 * t394 + t249 - t411;
t432 = t330 * t335;
t160 = t329 * t271 + t332 * t432;
t434 = t329 * t333;
t117 = -qJ(6) * t434 + t160;
t471 = pkin(5) * t329;
t393 = -t335 + t471;
t250 = t393 * t330;
t362 = t102 * t329 - t117 * t332 - t250;
t261 = mrSges(6,1) * t330 - mrSges(6,3) * t429;
t485 = -t261 / 0.2e1;
t455 = t330 * mrSges(7,1);
t260 = -mrSges(7,3) * t429 + t455;
t486 = -t260 / 0.2e1;
t399 = t486 + t485;
t257 = -mrSges(6,2) * t330 - mrSges(6,3) * t434;
t453 = t330 * mrSges(7,2);
t256 = -mrSges(7,3) * t434 - t453;
t489 = t256 / 0.2e1;
t401 = t489 + t257 / 0.2e1;
t504 = -m(7) / 0.2e1;
t528 = t362 * t504 + t233 / 0.2e1 + t232 / 0.2e1 + t401 * t332 + t399 * t329;
t527 = 0.4e1 * t424;
t506 = -m(6) / 0.2e1;
t526 = mrSges(7,3) / 0.2e1;
t525 = m(7) * (t102 * t332 + t117 * t329);
t468 = mrSges(6,2) + mrSges(7,2);
t524 = mrSges(6,3) + mrSges(7,3);
t467 = Ifges(6,5) + Ifges(7,5);
t466 = Ifges(6,6) + Ifges(7,6);
t523 = Ifges(6,3) + Ifges(7,3);
t414 = mrSges(7,1) + t501;
t323 = t329 ^ 2;
t325 = t332 ^ 2;
t427 = t323 + t325;
t521 = t330 * (-0.1e1 + t427);
t328 = cos(pkin(6));
t327 = sin(pkin(6));
t334 = cos(qJ(2));
t438 = t327 * t334;
t228 = t328 * t333 - t330 * t438;
t331 = sin(qJ(2));
t435 = t329 * t331;
t119 = t228 * t332 + t327 * t435;
t436 = t329 * t330;
t254 = -t333 * mrSges(7,2) + mrSges(7,3) * t436;
t255 = -t333 * mrSges(6,2) + mrSges(6,3) * t436;
t402 = t254 / 0.2e1 + t255 / 0.2e1;
t520 = t402 * t119;
t518 = t466 * t332;
t517 = t257 + t256;
t516 = t261 + t260;
t319 = Ifges(7,5) * t332;
t320 = Ifges(6,5) * t332;
t514 = -t319 - t320;
t321 = Ifges(7,4) * t332;
t513 = -Ifges(7,2) * t329 + t321;
t282 = Ifges(7,1) * t329 + t321;
t322 = Ifges(6,4) * t332;
t512 = -Ifges(6,2) * t329 + t322;
t284 = Ifges(6,1) * t329 + t322;
t288 = pkin(4) * t333 + pkin(9) * t330;
t253 = t332 * t288;
t428 = t333 * t335;
t172 = -t329 * t428 + t253;
t173 = t329 * t288 + t332 * t428;
t369 = -t172 * t329 + t173 * t332;
t433 = t330 * t332;
t103 = qJ(6) * t433 + t333 * t394 + t253;
t120 = qJ(6) * t436 + t173;
t511 = -t103 * t329 + t120 * t332;
t159 = -t329 * t432 + t249;
t116 = t159 - t411;
t510 = t116 / 0.2e1 - t102 / 0.2e1;
t508 = 2 * qJD(4);
t507 = m(5) / 0.2e1;
t505 = m(6) / 0.2e1;
t503 = m(7) / 0.2e1;
t502 = -pkin(5) / 0.2e1;
t499 = -t103 / 0.2e1;
t227 = t328 * t330 + t333 * t438;
t495 = t227 / 0.2e1;
t494 = t228 / 0.2e1;
t452 = t332 * mrSges(6,1);
t457 = t329 * mrSges(6,2);
t382 = t452 - t457;
t231 = t382 * t333;
t493 = -t231 / 0.2e1;
t234 = t515 * t333;
t492 = t234 / 0.2e1;
t258 = t333 * mrSges(7,1) + mrSges(7,3) * t433;
t488 = -t258 / 0.2e1;
t487 = t258 / 0.2e1;
t465 = -qJ(6) - pkin(9);
t272 = t465 * t329;
t484 = t272 / 0.2e1;
t481 = -t329 / 0.2e1;
t480 = t329 / 0.2e1;
t479 = -t332 / 0.2e1;
t251 = t393 * t333;
t476 = m(7) * t251;
t470 = pkin(5) * t332;
t308 = -pkin(4) - t470;
t475 = m(7) * t308;
t474 = m(7) * t330;
t473 = m(7) * t333;
t430 = t332 * t119;
t431 = t331 * t332;
t118 = -t228 * t329 + t327 * t431;
t437 = t329 * t118;
t371 = t430 - t437;
t361 = t228 - t371;
t16 = (-t227 * t521 - t361 * t333) * t532;
t464 = t16 * qJD(4);
t456 = t329 * mrSges(7,2);
t454 = t330 * mrSges(5,2);
t451 = t332 * mrSges(7,1);
t450 = -t382 - mrSges(5,1);
t448 = t118 * t332;
t447 = t119 * t329;
t396 = -t325 / 0.2e1 - t323 / 0.2e1;
t385 = t396 * mrSges(6,3);
t341 = -t401 * t329 + (mrSges(7,3) * t396 + t385) * t333;
t306 = mrSges(7,2) * t434;
t230 = mrSges(7,1) * t429 - t306;
t404 = t493 - t230 / 0.2e1;
t384 = t404 * t333;
t326 = t333 ^ 2;
t419 = t326 * t502;
t14 = t384 - t417 * t329 + t341 * t330 + ((t419 + t502) * m(7) - t418 + (t510 * m(7) + t399) * t330) * t332;
t445 = t14 * qJD(2);
t174 = (-t330 * t435 + t332 * t334) * t327;
t442 = t174 * t329;
t175 = (t329 * t334 + t330 * t431) * t327;
t441 = t175 * t332;
t439 = t327 * t331;
t413 = t333 * t439;
t18 = (t118 * t174 + t119 * t175 - t227 * t413) * t527 + m(5) * (-t227 * t333 + t228 * t330 + t438) * t439;
t440 = t18 * qJD(1);
t426 = qJD(4) * t333;
t422 = t250 * t503;
t420 = -t473 / 0.2e1;
t416 = mrSges(6,3) / 0.2e1 + t526;
t412 = t335 * t439;
t410 = t227 * t480;
t408 = -t439 / 0.2e1;
t235 = t333 * t277;
t403 = -t234 / 0.2e1 - t235 / 0.2e1;
t259 = t333 * mrSges(6,1) + mrSges(6,3) * t433;
t400 = t488 - t259 / 0.2e1;
t273 = -t451 + t456;
t398 = t382 / 0.2e1 - t273 / 0.2e1;
t395 = m(7) * (-t102 + t116);
t392 = t427 * mrSges(7,3);
t391 = t332 * t425;
t390 = -mrSges(6,1) - t414;
t386 = qJD(4) * (t273 + t450);
t383 = t333 * mrSges(5,1) - t454;
t280 = Ifges(6,2) * t332 + t461;
t278 = Ifges(7,2) * t332 + t460;
t275 = t465 * t332;
t368 = t441 - t442;
t337 = (pkin(4) * t413 + pkin(9) * t368) * t505 + (t174 * t272 - t175 * t275 - t308 * t413) * t503 + t524 * (t441 / 0.2e1 - t442 / 0.2e1);
t353 = t172 * t118 + t173 * t119 - t228 * t428;
t370 = -t159 * t329 + t160 * t332;
t355 = -t370 + t432;
t356 = t103 * t118 + t119 * t120 + t228 * t251;
t1 = t400 * t118 - t520 + t403 * t228 + t353 * t506 + t356 * t504 + t398 * t413 + t337 + (t355 * t506 + t528) * t227;
t203 = t330 * Ifges(7,6) + t333 * t513;
t204 = t330 * Ifges(6,6) + t333 * t512;
t354 = -t204 / 0.2e1 - t203 / 0.2e1 + (-Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t330;
t5 = -t160 * t255 - t159 * t259 - t120 * t256 - t117 * t254 - qJ(3) * t383 - t235 * t432 - t173 * t257 - t172 * t261 - t102 * t258 - t103 * t260 + t250 * t234 + t251 * t232 - m(6) * (t159 * t172 + t160 * t173) - m(7) * (t102 * t103 + t117 * t120 - t250 * t251) + (-Ifges(5,4) * t330 + t354 * t329 + t531 * t332) * t330 + (-t335 * t233 + (t466 * t329 + Ifges(5,4) + t514) * t333 + (m(6) * t335 ^ 2 + Ifges(5,1) - Ifges(5,2) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t325 + ((Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t329 + (-Ifges(6,4) - Ifges(7,4)) * t332) * t329 - t523) * t330) * t333;
t378 = -t1 * qJD(1) - t5 * qJD(2);
t363 = t395 / 0.2e1 + t486;
t358 = t485 + t363;
t338 = (t333 * t391 - t404) * t227 + (t416 * t434 + t401) * t118 + (-t416 * t429 + t358) * t119;
t6 = -t174 * t366 - t175 * t417 + t338;
t236 = t333 * t278;
t237 = t333 * t280;
t238 = t333 * t282;
t239 = t333 * t284;
t8 = t116 * t256 + t159 * t257 - t160 * t261 + t251 * t230 + (-t260 + t395) * t117 + (-t335 * t231 + (t102 * mrSges(7,3) + t159 * mrSges(6,3) + t236 / 0.2e1 + t237 / 0.2e1 - t531) * t329 + (-t117 * mrSges(7,3) - t160 * mrSges(6,3) - t238 / 0.2e1 - t239 / 0.2e1 + (t234 + t476) * pkin(5) + t354) * t332) * t333;
t377 = t6 * qJD(1) + t8 * qJD(2);
t372 = t447 + t448;
t340 = m(5) * t408 + (t504 + t506) * t372;
t291 = t326 * t439;
t324 = t330 ^ 2;
t348 = (t324 * t439 + t291) * t507 + (t330 * t368 + t291) * t532;
t20 = t340 + t348;
t364 = t330 * mrSges(5,1) + t333 * mrSges(5,2) + mrSges(4,3);
t27 = t516 * t332 + t517 * t329 + (m(5) + m(4)) * qJ(3) + t525 + m(6) * (t159 * t332 + t160 * t329) + t364;
t376 = qJD(1) * t20 - qJD(2) * t27;
t42 = (t329 * t256 + t332 * t260 + t525) * t333;
t50 = (t408 + t448 / 0.2e1 + t447 / 0.2e1) * t473;
t375 = -qJD(1) * t50 - qJD(2) * t42;
t17 = t361 * t227 * t527;
t374 = t17 * qJD(1) + t16 * qJD(3);
t367 = t272 * t329 + t275 * t332;
t162 = t414 * t429 - t306;
t240 = -m(7) * t471 - t515;
t365 = qJD(2) * t162 - qJD(4) * t240;
t359 = m(7) * (t272 * t332 - t275 * t329);
t336 = (t370 * t505 + t528) * t333 + (t402 * t332 + t400 * t329 + (t251 + t511) * t503 + (t369 - 0.2e1 * t428) * t505 - t403) * t330;
t10 = -t359 / 0.2e1 + t336;
t67 = t333 * t521 * t527;
t357 = t16 * qJD(1) + t10 * qJD(2) + t67 * qJD(3);
t352 = pkin(9) * t385 + t335 * t482;
t351 = -t284 / 0.4e1 - t282 / 0.4e1 - t512 / 0.4e1 - t513 / 0.4e1 + mrSges(7,3) * t484;
t11 = t530 * t227;
t25 = -pkin(4) * t277 + t308 * t515 + (t282 / 0.2e1 + t513 / 0.2e1 + t284 / 0.2e1 + t512 / 0.2e1) * t332 + (t283 / 0.2e1 - t278 / 0.2e1 + t285 / 0.2e1 - t280 / 0.2e1 + (t273 + t475) * pkin(5)) * t329;
t339 = (t320 / 0.4e1 + t319 / 0.4e1) * t330 - t363 * t275 + pkin(4) * t493 + t251 * t483 + t256 * t484 + t308 * t230 / 0.2e1;
t342 = t285 / 0.4e1 + t283 / 0.4e1 - t280 / 0.4e1 - t278 / 0.4e1 + t275 * t526 + (t475 / 0.2e1 + t273 / 0.2e1) * pkin(5);
t343 = (t476 / 0.2e1 + t492) * pkin(5) - t203 / 0.4e1 - t204 / 0.4e1 - t238 / 0.4e1 - t239 / 0.4e1 - pkin(9) * t257 / 0.2e1;
t344 = pkin(9) * t485 - t237 / 0.4e1 - t236 / 0.4e1 + t206 / 0.4e1 + t205 / 0.4e1 + t510 * mrSges(7,3);
t347 = mrSges(7,1) * t499 + t120 * mrSges(7,2) / 0.2e1 - t172 * mrSges(6,1) / 0.2e1 + t173 * mrSges(6,2) / 0.2e1;
t4 = (m(7) * t499 + t488) * pkin(5) + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t352) * t333 + (t333 * t342 + t344 + t387) * t332 + ((-0.3e1 / 0.4e1 * Ifges(6,6) - 0.3e1 / 0.4e1 * Ifges(7,6)) * t330 + t351 * t333 + t343) * t329 + t339 + t347;
t44 = t530 * t333;
t350 = t11 * qJD(1) + t4 * qJD(2) - t44 * qJD(3) + t25 * qJD(4);
t164 = (-0.1e1 / 0.2e1 - t396) * t474;
t346 = m(7) * ((-t272 * t333 + t117) * t332 + (t275 * t333 - t102) * t329);
t31 = (-t453 / 0.2e1 - t256 / 0.2e1) * t332 + (-t455 / 0.2e1 + t260 / 0.2e1) * t329 - t422 - t346 / 0.2e1;
t49 = 0.2e1 * (-t437 / 0.4e1 + t430 / 0.4e1 - t228 / 0.4e1) * m(7);
t83 = -m(7) * t367 + t392;
t349 = qJD(1) * t49 - qJD(2) * t31 + qJD(3) * t164 + qJD(4) * t83;
t304 = qJ(3) * t438;
t264 = t326 * t412;
t163 = (t427 + 0.1e1) * t474 / 0.2e1;
t51 = t372 * t420 + t408 * t473;
t48 = m(7) * t494 + t371 * t503;
t45 = t420 * t471 + t529 * t333 - (t277 + t515) * t333 / 0.2e1;
t32 = t260 * t481 + t346 / 0.2e1 + t332 * t489 - t422 + (-t318 / 0.2e1 - t317 / 0.2e1) * t330;
t19 = m(4) * t439 - t340 + t348;
t13 = m(7) * t332 * t419 + t391 - t456 / 0.2e1 + t451 / 0.2e1 + t452 / 0.2e1 - t457 / 0.2e1 + t384 + (t332 * t358 + t341) * t330;
t12 = t277 * t495 + t410 * t501 + (t483 - t529) * t227;
t9 = t359 / 0.2e1 + t336;
t7 = t174 * t425 + t338 + (mrSges(6,1) + mrSges(7,1)) * t174 / 0.2e1 - t468 * t175 / 0.2e1;
t3 = t339 + ((-Ifges(6,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t330 + t343) * t329 + t344 * t332 + t103 * t425 + pkin(5) * t487 - t347 + t466 * t436 / 0.2e1 - t467 * t433 / 0.2e1 + (t351 * t329 + t342 * t332 + t352 + t523 / 0.2e1) * t333;
t2 = t337 + (t227 * t355 + t353) * t505 + (t227 * t362 + t356) * t503 + t228 * t492 + t235 * t494 + (-t232 - t233) * t495 + (t383 / 0.2e1 - t454 / 0.2e1 + (mrSges(5,1) / 0.2e1 + t398) * t333) * t439 + t516 * t410 + t517 * t227 * t479 + t520 + (t487 + t259 / 0.2e1) * t118;
t15 = [qJD(2) * t18 + qJD(4) * t17, t19 * qJD(3) + t2 * qJD(4) + t7 * qJD(5) + t51 * qJD(6) + t440 + (t174 * t260 + t174 * t261 + t175 * t256 + t175 * t257 + ((-mrSges(3,2) + t364) * t334 + (-mrSges(3,1) + mrSges(4,2) + (-t234 - t235) * t333 + (-t324 - t326) * mrSges(5,3)) * t331) * t327 + 0.2e1 * (t159 * t174 + t160 * t175 + t264) * t505 + 0.2e1 * (t102 * t174 + t117 * t175 - t251 * t413) * t503 + 0.2e1 * (t324 * t412 + t264 + t304) * t507 + m(4) * (-pkin(2) * t439 + t304)) * qJD(2), qJD(2) * t19 + t464, t2 * qJD(2) + t12 * qJD(5) + t48 * qJD(6) + t228 * t386 + (-t524 * t427 + mrSges(5,2)) * qJD(4) * t227 + ((-pkin(9) * t227 * t427 - pkin(4) * t228) * t505 + (t227 * t367 + t228 * t308) * t503) * t508 + t374, t7 * qJD(2) + t12 * qJD(4) + (-t118 * t468 + t119 * t390) * qJD(5), qJD(2) * t51 + qJD(4) * t48; -qJD(3) * t20 - qJD(4) * t1 + qJD(5) * t6 - qJD(6) * t50 - t440, qJD(3) * t27 - qJD(4) * t5 + qJD(5) * t8 - qJD(6) * t42, qJD(4) * t9 + qJD(5) * t13 - t376, t9 * qJD(3) + t3 * qJD(5) + t32 * qJD(6) + (-t335 * mrSges(5,2) + t329 * t467 - Ifges(5,6) + t518) * t426 + t378 + (m(7) * (t103 * t272 - t120 * t275 - t250 * t308) + pkin(4) * t233 + t272 * t258 - t250 * t273 - t275 * t254 - t308 * t232 + (m(6) * t369 + t332 * t255 - t329 * t259) * pkin(9) + t511 * mrSges(7,3) + t369 * mrSges(6,3) + (-Ifges(5,5) + (-m(6) * pkin(4) + t450) * t335 + (t285 + t283) * t481 + (t278 + t280) * t480 + (t282 + t284 + t513 + t512) * t479) * t330) * qJD(4), t13 * qJD(3) + t3 * qJD(4) + t377 + (-mrSges(6,1) * t160 - mrSges(6,2) * t159 - mrSges(7,2) * t116 + (-t518 + (mrSges(7,3) * pkin(5) - t467) * t329) * t333 - t414 * t117) * qJD(5), qJD(4) * t32 + t375; qJD(2) * t20 + t464, qJD(4) * t10 + qJD(5) * t14 + t376, t67 * qJD(4), t45 * qJD(5) + t163 * qJD(6) + t330 * t386 + (mrSges(6,3) * t427 - mrSges(5,2) + t392) * t426 + ((t427 * t469 - t472) * t505 + (t308 * t330 - t333 * t367) * t503) * t508 + t357, t445 + t45 * qJD(4) + (t329 * t468 + t332 * t390) * qJD(5) * t330, t163 * qJD(4); qJD(2) * t1 + qJD(5) * t11 + qJD(6) * t49 - t374, -qJD(3) * t10 + qJD(5) * t4 - qJD(6) * t31 - t378, -qJD(5) * t44 + qJD(6) * t164 - t357, qJD(5) * t25 + qJD(6) * t83, t350 + (-mrSges(7,2) * t272 - mrSges(7,3) * t470 - pkin(9) * t452 + (mrSges(6,2) * pkin(9) - t466) * t329 + t414 * t275 - t514) * qJD(5), t349; -qJD(2) * t6 - qJD(4) * t11, -qJD(3) * t14 - qJD(4) * t4 - qJD(6) * t162 - t377, qJD(4) * t44 - t445, qJD(6) * t240 - t350, 0, -t365; qJD(2) * t50 - qJD(4) * t49, qJD(4) * t31 + qJD(5) * t162 - t375, -t164 * qJD(4), -qJD(5) * t240 - t349, t365, 0;];
Cq  = t15;
