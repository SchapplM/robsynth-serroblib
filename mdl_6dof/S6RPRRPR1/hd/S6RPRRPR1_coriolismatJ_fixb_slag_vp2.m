% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:47
% EndTime: 2019-03-09 04:57:59
% DurationCPUTime: 7.72s
% Computational Cost: add. (22506->408), mult. (41975->550), div. (0->0), fcn. (48204->10), ass. (0->259)
t286 = sin(qJ(4));
t287 = sin(qJ(3));
t289 = cos(qJ(4));
t290 = cos(qJ(3));
t261 = -t286 * t287 + t289 * t290;
t365 = -cos(pkin(10)) * pkin(1) - pkin(2);
t263 = -pkin(3) * t290 + t365;
t231 = -t261 * pkin(4) + t263;
t528 = m(6) * t231;
t262 = -t286 * t290 - t289 * t287;
t283 = sin(pkin(11));
t417 = cos(pkin(11));
t227 = t417 * t261 + t262 * t283;
t285 = sin(qJ(6));
t288 = cos(qJ(6));
t314 = t283 * t261 - t262 * t417;
t135 = -pkin(5) * t227 - pkin(9) * t314 + t231;
t272 = sin(pkin(10)) * pkin(1) + pkin(7);
t452 = pkin(8) + t272;
t350 = t452 * t287;
t487 = t452 * t290;
t221 = -t286 * t350 + t289 * t487;
t302 = t261 * qJ(5) + t221;
t494 = -t286 * t487 - t289 * t350;
t501 = t262 * qJ(5) + t494;
t510 = t283 * t501 + t302 * t417;
t53 = t135 * t285 + t288 * t510;
t424 = t288 * t53;
t52 = t135 * t288 - t285 * t510;
t331 = -t285 * t52 + t424;
t511 = -t283 * t302 + t417 * t501;
t409 = t511 * t314;
t305 = t227 * t331 - t409;
t516 = t510 * t227;
t527 = t305 - t516;
t502 = t494 * mrSges(5,2);
t507 = t221 * mrSges(5,1);
t521 = t511 * mrSges(6,2);
t526 = -t502 / 0.2e1 - t507 / 0.2e1 - t521 / 0.2e1;
t264 = -mrSges(7,1) * t288 + mrSges(7,2) * t285;
t515 = t510 * t264;
t522 = t510 * mrSges(6,1);
t525 = t515 / 0.2e1 - t522 / 0.2e1;
t278 = Ifges(7,4) * t288;
t484 = -Ifges(7,2) * t285 + t278;
t113 = Ifges(7,6) * t314 + t227 * t484;
t443 = Ifges(7,4) * t285;
t269 = Ifges(7,1) * t288 - t443;
t115 = Ifges(7,5) * t314 + t227 * t269;
t334 = Ifges(7,5) * t285 + Ifges(7,6) * t288;
t456 = t288 / 0.2e1;
t458 = t285 / 0.2e1;
t468 = t314 / 0.2e1;
t441 = Ifges(7,2) * t288;
t266 = t441 + t443;
t459 = -t285 / 0.2e1;
t257 = t266 * t459;
t445 = Ifges(7,1) * t285;
t268 = t278 + t445;
t495 = t268 * t456 + t257;
t307 = Ifges(5,5) * t261 + Ifges(5,6) * t262 - Ifges(6,6) * t314 + t113 * t456 + t115 * t458 + t334 * t468 + (Ifges(6,5) + t495) * t227;
t524 = t307 - t502 - t507 + t515 - t521 - t522;
t425 = t288 * mrSges(7,2);
t428 = t285 * mrSges(7,1);
t265 = t425 + t428;
t159 = t265 * t314;
t386 = t285 * t227;
t160 = -mrSges(7,2) * t314 - mrSges(7,3) * t386;
t397 = t227 * t288;
t162 = mrSges(7,1) * t314 - mrSges(7,3) * t397;
t165 = mrSges(6,1) * t314 + t227 * mrSges(6,2);
t230 = -t262 * mrSges(5,1) + t261 * mrSges(5,2);
t500 = Ifges(7,5) * t227;
t116 = t269 * t314 - t500;
t385 = t288 * t116;
t499 = Ifges(7,6) * t227;
t114 = t314 * t484 - t499;
t390 = t285 * t114;
t467 = -t227 / 0.2e1;
t277 = Ifges(7,5) * t288;
t485 = -Ifges(7,6) * t285 + t277;
t496 = t265 * t227;
t523 = t510 * t159 + t53 * t160 + t52 * t162 + t231 * t165 + t263 * t230 - t496 * t511 + (t485 * t467 + t385 / 0.2e1 - t390 / 0.2e1 + Ifges(6,4) * t227) * t227 + (-Ifges(6,4) * t314 + t113 * t459 + t115 * t456 + t468 * t485 + (-Ifges(6,2) - Ifges(7,3) + Ifges(6,1)) * t227) * t314;
t455 = pkin(3) * t289;
t276 = pkin(4) + t455;
t391 = t283 * t286;
t243 = -pkin(3) * t391 + t276 * t417;
t238 = -pkin(5) - t243;
t520 = t238 * t510;
t364 = t417 * pkin(4);
t273 = -t364 - pkin(5);
t519 = t273 * t510;
t518 = t285 * t511;
t517 = t288 * t511;
t412 = t510 * t511;
t347 = t417 * t286;
t244 = pkin(3) * t347 + t283 * t276;
t514 = t511 * t244;
t513 = t283 * t511 - t417 * t510;
t498 = t227 * mrSges(6,3);
t492 = mrSges(6,3) * t314;
t491 = -mrSges(6,1) + t264;
t488 = t314 * t264;
t250 = (t289 * t417 - t391) * pkin(3);
t281 = t285 ^ 2;
t282 = t288 ^ 2;
t376 = t281 + t282;
t345 = t376 * t250;
t447 = mrSges(7,3) * t314;
t371 = t285 * t447;
t161 = mrSges(7,2) * t227 - t371;
t163 = -mrSges(7,1) * t227 - t288 * t447;
t453 = pkin(4) * t283;
t271 = pkin(9) + t453;
t298 = (t269 / 0.4e1 - t266 / 0.4e1 - t441 / 0.4e1) * t288 + (-t268 / 0.4e1 - t484 / 0.4e1 - t278 / 0.2e1 - t445 / 0.4e1) * t285;
t351 = t282 / 0.2e1 + t281 / 0.2e1;
t340 = mrSges(7,3) * t351;
t363 = -t511 * t265 / 0.2e1;
t460 = -t273 / 0.2e1;
t294 = (-t271 * t340 + t298) * t314 + t363 + t488 * t460;
t313 = 0.3e1 / 0.4e1 * t499 - t114 / 0.4e1;
t454 = pkin(4) * t262;
t141 = pkin(5) * t314 - pkin(9) * t227 - t454;
t66 = t141 * t285 + t517;
t469 = t66 / 0.2e1;
t65 = t141 * t288 - t518;
t470 = -t65 / 0.2e1;
t320 = mrSges(7,1) * t470 + mrSges(7,2) * t469;
t339 = t116 / 0.4e1 - t500 / 0.2e1;
t461 = -t271 / 0.2e1;
t439 = Ifges(7,3) * t314;
t466 = -t227 / 0.4e1;
t480 = t277 * t466 - t439 / 0.2e1;
t12 = (t161 * t461 + t313) * t285 + (t163 * t461 + t339) * t288 + t294 + t320 + t480;
t342 = t269 * t458 + t456 * t484 + t495;
t392 = t273 * t265;
t164 = t342 + t392;
t319 = -t425 / 0.2e1 - t428 / 0.2e1;
t311 = t319 * t227;
t88 = t496 / 0.2e1 + t311;
t379 = t88 * qJD(2);
t462 = -t250 / 0.2e1;
t465 = -t238 / 0.2e1;
t79 = (t465 + t460) * t265 + (mrSges(7,2) * t462 - t268 / 0.2e1 - t484 / 0.2e1) * t288 + (mrSges(7,1) * t462 - t269 / 0.2e1 + t266 / 0.2e1) * t285;
t483 = -t12 * qJD(1) + t79 * qJD(3) - t164 * qJD(4) + t379;
t239 = pkin(9) + t244;
t295 = (-t239 * t340 + t298) * t314 + t363 + t488 * t465;
t280 = t287 * pkin(3);
t139 = t141 + t280;
t59 = t139 * t288 - t518;
t60 = t139 * t285 + t517;
t321 = -t59 * mrSges(7,1) / 0.2e1 + t60 * mrSges(7,2) / 0.2e1;
t464 = -t239 / 0.2e1;
t10 = (t161 * t464 + t313) * t285 + (t163 * t464 + t339) * t288 + t295 + t321 + t480;
t395 = t238 * t265;
t148 = t342 + t395;
t482 = t10 * qJD(1) + t148 * qJD(3) - t379;
t383 = t288 * t161;
t387 = t285 * t163;
t317 = t383 / 0.2e1 - t387 / 0.2e1;
t384 = t288 * t160;
t388 = t285 * t162;
t481 = t384 / 0.2e1 - t388 / 0.2e1;
t329 = -t65 * t285 + t66 * t288;
t479 = t160 * t458 + t162 * t456;
t478 = (-mrSges(5,1) * t286 - mrSges(5,2) * t289) * pkin(3);
t477 = -Ifges(7,6) * t386 / 0.2e1 + Ifges(7,5) * t397 / 0.2e1 + t439 / 0.2e1 + t485 * t466 + t385 / 0.4e1 - t390 / 0.4e1;
t396 = t238 * t496;
t400 = t227 * t243;
t405 = t314 * t244;
t249 = (t283 * t289 + t347) * pkin(3);
t408 = t511 * t249;
t463 = t249 / 0.2e1;
t472 = m(7) / 0.2e1;
t473 = m(6) / 0.2e1;
t476 = (t514 - t408 + (-t243 + t250) * t510) * t473 + (t239 * t329 + t250 * t331 - t408 + t520) * t472 + t396 / 0.2e1 + t159 * t463 + (-t227 * t462 + t314 * t463 - t400 / 0.2e1 - t405 / 0.2e1) * mrSges(6,3) + t526;
t475 = 0.2e1 * qJD(3);
t471 = m(6) * pkin(4);
t457 = -t288 / 0.2e1;
t306 = m(7) * (t376 - 0.1e1) * t227 * t314;
t381 = t306 * qJD(2);
t89 = -t496 / 0.2e1 + t311;
t451 = t89 * qJD(6) + t381;
t450 = -t88 * qJD(6) - t381;
t444 = Ifges(5,4) * t262;
t427 = t285 * t59;
t426 = t285 * t60;
t423 = t288 * t59;
t422 = t288 * t60;
t297 = t159 * t468 + t317 * t227 + t314 * t481 + t467 * t496;
t330 = t422 - t427;
t7 = (t314 * t330 + t527) * t472 + t297;
t419 = t7 * qJD(1);
t8 = (t314 * t329 + t527) * t472 + t297;
t418 = t8 * qJD(1);
t344 = t376 * t271;
t322 = (t227 * t344 + t273 * t314) * t472 + (t227 * t283 - t314 * t417) * pkin(4) * t473;
t316 = t161 * t459 + t163 * t457;
t22 = t488 * t467 + (t351 * t447 - t316) * t314;
t416 = qJD(1) * t22;
t128 = t244 * t227 - t243 * t314;
t237 = t280 - t454;
t399 = t227 * t281;
t215 = mrSges(7,3) * t399;
t398 = t227 * t282;
t216 = mrSges(7,3) * t398;
t338 = t488 / 0.2e1 + t215 / 0.2e1 + t216 / 0.2e1;
t300 = t338 - t165 - t479;
t69 = t238 * t314 + (t398 + t399) * t239;
t17 = 0.2e1 * (t69 / 0.4e1 - t426 / 0.4e1 - t423 / 0.4e1) * m(7) + 0.2e1 * (t128 / 0.4e1 - t237 / 0.4e1) * m(6) + t300;
t407 = t17 * qJD(1);
t374 = -t471 / 0.2e1;
t301 = (t285 * t66 + t288 * t65) * t472 + t262 * t374;
t20 = t300 - t301 + t322;
t406 = t20 * qJD(1);
t403 = t227 * t249;
t25 = t311 - t317;
t394 = t25 * qJD(1);
t393 = t273 * t496;
t40 = t306 / 0.2e1;
t380 = t40 * qJD(1);
t375 = qJD(3) + qJD(4);
t373 = mrSges(7,3) * t427;
t372 = mrSges(7,3) * t422;
t370 = mrSges(7,3) * t470;
t369 = mrSges(7,3) * t469;
t368 = t271 * t388;
t367 = t271 * t384;
t343 = t453 * t492;
t337 = Ifges(5,4) * t261 + (-Ifges(5,1) + Ifges(5,2)) * t262;
t166 = -mrSges(6,1) * t227 + mrSges(6,2) * t314;
t3 = m(5) * t263 * t280 + m(7) * (t52 * t59 + t53 * t60 - t412) + (t365 * mrSges(4,2) + Ifges(4,4) * t290 + (Ifges(4,1) - Ifges(4,2)) * t287) * t290 + t60 * t161 + t59 * t163 + (-mrSges(5,1) * t280 + t337) * t261 + (-mrSges(5,2) * t280 - t444) * t262 + (mrSges(4,1) * t365 - Ifges(4,4) * t287) * t287 + (t166 + t528) * t237 + t523;
t333 = t3 * qJD(1) + t7 * qJD(2);
t4 = t337 * t261 - t454 * t528 + m(7) * (t52 * t65 + t53 * t66 - t412) + t66 * t161 + t65 * t163 + (-pkin(4) * t166 - t444) * t262 + t523;
t332 = t4 * qJD(1) + t8 * qJD(2);
t328 = t364 * t498;
t6 = -t511 * t488 + t53 * t163 + (t116 * t458 + t114 * t456 + mrSges(7,3) * t424 + t334 * t467 + (-t268 * t457 + t257) * t314) * t314 + (-t371 - t161) * t52;
t327 = -t6 * qJD(1) - t22 * qJD(2);
t13 = (t159 + t492) * t314 + (t383 - t387 + t498) * t227 + m(7) * t305 + m(6) * (-t409 + t516);
t326 = -qJD(1) * t13 - qJD(2) * t40;
t323 = t128 * t473 + t472 * t69;
t312 = t338 + t479;
t296 = (t250 * t314 + t128 - t403) * t473 + (t314 * t345 - t403 + t69) * t472;
t24 = t296 - t322;
t42 = t491 * t249 + t478 + (mrSges(7,3) * t376 - mrSges(6,2)) * t250 + m(6) * (-t243 * t249 + t244 * t250) + m(7) * (t238 * t249 + t239 * t345);
t291 = -m(7) * (t271 * t330 + t519) / 0.2e1 - t393 / 0.2e1 + t513 * t374 + t368 / 0.2e1 - t367 / 0.2e1 + t373 / 0.2e1 - t372 / 0.2e1 + t343 / 0.2e1 + t328 / 0.2e1 - t525 - t526;
t5 = t481 * t239 + t317 * t250 + t285 * t370 + t288 * t369 + t291 + t476 + t525;
t310 = t5 * qJD(1) + t24 * qJD(2) + t42 * qJD(3);
t304 = -t165 + t488 + t215 + t216 - t230;
t80 = t395 / 0.2e1 + t392 / 0.2e1 + t319 * t250 + t342;
t26 = t311 + t317;
t21 = t301 + t312 + t322;
t19 = t296 + t304 + t322;
t18 = (t423 + t426) * t472 + t237 * t473 + t312 + t323;
t11 = t271 * t316 + t294 - t320 + t477;
t9 = t239 * t316 + t295 - t321 + t477;
t2 = t307 + (t250 * t161 / 0.2e1 + t239 * t160 / 0.2e1 + t369) * t288 + (t162 * t464 + t163 * t462 + t370) * t285 + (t264 / 0.2e1 - mrSges(6,1) / 0.2e1) * t510 - t291 + t476;
t1 = qJD(3) * t7 + qJD(4) * t8 + qJD(5) * t40 - qJD(6) * t22;
t14 = [qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t13 - qJD(6) * t6, t1, t2 * qJD(4) + t18 * qJD(5) + t9 * qJD(6) + ((t239 * t330 + t520) * t472 + (-t243 * t510 + t514) * t473) * t475 + t333 + (Ifges(4,5) * t290 - Ifges(4,6) * t287 + t372 - t373 + t396 + (m(5) * (-t221 * t289 + t286 * t494) + (-t261 * t289 + t262 * t286) * mrSges(5,3)) * pkin(3) + (-mrSges(4,1) * t290 + mrSges(4,2) * t287) * t272 + (t384 - t388) * t239 + (-t400 - t405) * mrSges(6,3) + t524) * qJD(3), t2 * qJD(3) + (-t328 + t513 * t471 + m(7) * (t271 * t329 + t519) + t393 - t343 + t367 - t368 + t329 * mrSges(7,3) + t524) * qJD(4) + t21 * qJD(5) + t11 * qJD(6) + t332, qJD(3) * t18 + qJD(4) * t21 + qJD(6) * t26 - t326, t9 * qJD(3) + t11 * qJD(4) + t26 * qJD(5) + (-mrSges(7,1) * t53 - mrSges(7,2) * t52 - t314 * t334) * qJD(6) + t327; t1, t375 * t306, t419 + (-t287 * mrSges(4,1) - t290 * mrSges(4,2) + t304) * qJD(3) + t19 * qJD(4) + ((pkin(3) * t261 * t286 + t262 * t455) * m(5) / 0.2e1 + t323) * t475 + t451, t19 * qJD(3) + t418 + t451 + (t304 + 0.2e1 * t322) * qJD(4), t380, qJD(6) * t488 + t375 * t89 - t416; qJD(4) * t5 + qJD(5) * t17 + qJD(6) * t10 - t333, qJD(4) * t24 - t419 + t450, qJD(4) * t42 + qJD(6) * t148 (mrSges(7,3) * t345 + (m(7) * t344 + t283 * t471 - mrSges(6,2)) * t250 + (m(7) * t273 - t417 * t471 + t491) * t249 + t478) * qJD(4) + t80 * qJD(6) + t310, t407, t80 * qJD(4) + (t239 * t264 + t485) * qJD(6) + t482; -qJD(3) * t5 + qJD(5) * t20 + qJD(6) * t12 - t332, -qJD(3) * t24 - t418 + t450, -qJD(6) * t79 - t310, t164 * qJD(6), t406 (t264 * t271 + t485) * qJD(6) - t483; -qJD(3) * t17 - qJD(4) * t20 - qJD(6) * t25 + t326, -t380, -t407, -t406, 0, -qJD(6) * t265 - t394; -qJD(3) * t10 - qJD(4) * t12 + qJD(5) * t25 - t327, t375 * t88 + t416, qJD(4) * t79 - t482, t483, t394, 0;];
Cq  = t14;
