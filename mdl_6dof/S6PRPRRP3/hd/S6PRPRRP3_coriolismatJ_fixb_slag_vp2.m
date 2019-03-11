% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:51
% EndTime: 2019-03-08 20:05:04
% DurationCPUTime: 6.77s
% Computational Cost: add. (11323->510), mult. (26114->708), div. (0->0), fcn. (28756->10), ass. (0->261)
t471 = -m(7) / 0.2e1;
t290 = cos(qJ(5));
t444 = -t290 / 0.2e1;
t286 = sin(pkin(11));
t288 = cos(pkin(11));
t439 = sin(qJ(4));
t441 = cos(qJ(4));
t263 = -t286 * t441 - t288 * t439;
t289 = sin(qJ(5));
t424 = -qJ(6) - pkin(9);
t266 = t424 * t289;
t269 = t424 * t290;
t261 = t286 * t439 - t288 * t441;
t433 = pkin(5) * t261;
t396 = t263 * t290;
t276 = -pkin(3) * t288 - pkin(2);
t428 = pkin(9) * t263;
t185 = pkin(4) * t261 + t276 + t428;
t425 = pkin(8) + qJ(3);
t265 = t425 * t288;
t348 = t425 * t286;
t204 = t265 * t441 - t348 * t439;
t97 = t185 * t290 - t204 * t289;
t74 = qJ(6) * t396 + t97;
t59 = t74 + t433;
t397 = t263 * t289;
t98 = t185 * t289 + t204 * t290;
t75 = qJ(6) * t397 + t98;
t325 = t289 * t59 - t290 * t75;
t190 = -mrSges(7,2) * t261 + mrSges(7,3) * t397;
t194 = mrSges(7,1) * t261 + mrSges(7,3) * t396;
t443 = t290 / 0.2e1;
t447 = -t289 / 0.2e1;
t387 = t190 * t443 + t194 * t447;
t507 = ((t266 * t290 - t269 * t289) * t263 - t325) * t471 - t387;
t416 = t261 * mrSges(6,2);
t191 = mrSges(6,3) * t397 - t416;
t506 = t191 * t444 - t387;
t505 = Ifges(6,1) + Ifges(7,1);
t504 = Ifges(6,2) + Ifges(7,2);
t489 = mrSges(6,3) + mrSges(7,3);
t488 = Ifges(6,5) + Ifges(7,5);
t426 = Ifges(6,6) + Ifges(7,6);
t384 = t286 ^ 2 + t288 ^ 2;
t503 = t384 * mrSges(4,3);
t502 = t290 * t426;
t501 = -t190 - t191;
t417 = t261 * mrSges(6,1);
t195 = mrSges(6,3) * t396 + t417;
t500 = t194 + t195;
t418 = Ifges(7,4) * t290;
t329 = -Ifges(7,2) * t289 + t418;
t420 = Ifges(6,4) * t290;
t331 = -Ifges(6,2) * t289 + t420;
t499 = t331 + t329;
t419 = Ifges(7,4) * t289;
t333 = t290 * Ifges(7,1) - t419;
t421 = Ifges(6,4) * t289;
t335 = t290 * Ifges(6,1) - t421;
t498 = t335 + t333;
t287 = sin(pkin(6));
t291 = cos(qJ(2));
t394 = t287 * t291;
t221 = t261 * t394;
t440 = sin(qJ(2));
t365 = t287 * t440;
t176 = t221 * t289 + t290 * t365;
t177 = -t221 * t290 + t289 * t365;
t473 = -m(6) / 0.2e1;
t496 = -(t471 + t473) * (t176 * t290 + t177 * t289) + (m(4) + m(5)) * t365 / 0.2e1;
t493 = t261 / 0.2e1;
t284 = t289 ^ 2;
t285 = t290 ^ 2;
t383 = t284 + t285;
t490 = t383 / 0.2e1;
t203 = t265 * t439 + t348 * t441;
t139 = -pkin(5) * t397 + t203;
t438 = m(7) * t139;
t469 = m(7) * pkin(5);
t487 = -mrSges(7,1) - t469;
t268 = -mrSges(6,1) * t290 + mrSges(6,2) * t289;
t486 = t268 - mrSges(5,1);
t124 = t261 * Ifges(7,6) - t263 * t329;
t126 = t261 * Ifges(6,6) - t263 * t331;
t485 = t124 + t126;
t128 = t261 * Ifges(7,5) - t263 * t333;
t130 = t261 * Ifges(6,5) - t263 * t335;
t484 = t128 + t130;
t399 = t261 * t289;
t188 = mrSges(7,2) * t263 + mrSges(7,3) * t399;
t189 = mrSges(6,2) * t263 + mrSges(6,3) * t399;
t483 = t188 + t189;
t482 = -t266 * t289 - t269 * t290;
t279 = t289 * mrSges(7,1);
t280 = t290 * mrSges(7,2);
t385 = t280 + t279;
t410 = cos(pkin(6));
t244 = t286 * t410 + t288 * t365;
t303 = t286 * t365 - t288 * t410;
t481 = t288 * t244 + t286 * t303;
t480 = t290 * t504 + t419 + t421;
t479 = t289 * t505 + t418 + t420;
t454 = -t194 / 0.2e1;
t352 = -t195 / 0.2e1 + t454;
t354 = t191 / 0.2e1 + t190 / 0.2e1;
t181 = t385 * t263;
t422 = mrSges(6,2) * t290;
t337 = mrSges(6,1) * t289 + t422;
t182 = t337 * t263;
t476 = mrSges(5,3) * t263 + t181 + t182;
t373 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t475 = t498 * t493 + (-t373 + t488 / 0.2e1) * t263;
t474 = 2 * qJD(4);
t472 = m(6) / 0.2e1;
t470 = m(7) / 0.2e1;
t468 = mrSges(6,1) / 0.2e1;
t467 = -mrSges(7,2) / 0.2e1;
t466 = -mrSges(7,3) / 0.2e1;
t434 = pkin(4) * t263;
t198 = pkin(9) * t261 - t434;
t101 = t198 * t290 + t203 * t289;
t398 = t261 * t290;
t62 = -pkin(5) * t263 + qJ(6) * t398 + t101;
t465 = t62 / 0.2e1;
t464 = pkin(4) * mrSges(6,1);
t463 = pkin(4) * mrSges(6,2);
t154 = t244 * t439 + t303 * t441;
t461 = t154 / 0.2e1;
t155 = t244 * t441 - t303 * t439;
t460 = t155 / 0.2e1;
t459 = t181 / 0.2e1;
t458 = -t181 / 0.2e1;
t192 = -mrSges(7,1) * t263 + mrSges(7,3) * t398;
t455 = t192 / 0.2e1;
t452 = -t203 / 0.2e1;
t267 = -mrSges(7,1) * t290 + mrSges(7,2) * t289;
t449 = -t267 / 0.2e1;
t431 = pkin(5) * t290;
t277 = -pkin(4) - t431;
t448 = t277 / 0.2e1;
t445 = t289 / 0.2e1;
t436 = m(7) * t263;
t435 = m(7) * t277;
t432 = pkin(5) * t289;
t430 = pkin(9) * t191;
t429 = pkin(9) * t195;
t427 = Ifges(7,4) + Ifges(6,4);
t423 = -t59 + t74;
t415 = t261 * mrSges(5,3);
t413 = t266 * mrSges(7,3);
t412 = t269 * mrSges(7,3);
t411 = t59 * t290;
t248 = t261 * mrSges(5,2);
t350 = -t285 / 0.2e1 - t284 / 0.2e1;
t339 = t350 * mrSges(6,3);
t349 = pkin(9) * t383;
t293 = (mrSges(7,3) * t350 + t339) * t261 + (-t261 * t349 + t434) * t472 + (-t261 * t482 - t263 * t277) * t470;
t102 = t198 * t289 - t203 * t290;
t79 = qJ(6) * t399 + t102;
t298 = (t101 * t290 + t102 * t289) * t473 + (t289 * t79 + t290 * t62) * t471;
t351 = -t268 / 0.2e1 + t449;
t193 = -mrSges(6,1) * t263 + mrSges(6,3) * t398;
t353 = -t192 / 0.2e1 - t193 / 0.2e1;
t355 = -t188 / 0.2e1 - t189 / 0.2e1;
t14 = t248 + t353 * t290 + t355 * t289 + (mrSges(5,1) + t351) * t263 + t293 + t298;
t409 = qJD(2) * t14;
t408 = t139 * t289;
t220 = t263 * t394;
t407 = t154 * t220;
t406 = t154 * t263;
t405 = t176 * t289;
t404 = t177 * t290;
t120 = t155 * t290 - t289 * t394;
t390 = t290 * t120;
t119 = -t155 * t289 - t290 * t394;
t392 = t289 * t119;
t321 = -t390 + t392;
t381 = m(6) / 0.4e1 + m(7) / 0.4e1;
t19 = 0.4e1 * t381 * (t155 + t321) * t154;
t403 = t19 * qJD(1);
t366 = t287 ^ 2 * t440;
t20 = m(5) * (-t155 * t221 - t291 * t366 - t407) + m(4) * (t287 * t481 - t366) * t291 + (m(6) + m(7)) * (t119 * t176 + t120 * t177 - t407);
t402 = t20 * qJD(1);
t401 = t203 * t220;
t400 = t203 * t263;
t386 = (t280 / 0.2e1 + t279 / 0.2e1) * t261;
t111 = 0.2e1 * (t284 / 0.4e1 + t285 / 0.4e1 + 0.1e1 / 0.4e1) * t436;
t382 = t111 * qJD(2);
t380 = pkin(5) * t396;
t376 = t220 * t470;
t375 = mrSges(6,2) / 0.2e1 + mrSges(7,2) / 0.2e1;
t374 = t466 - mrSges(6,3) / 0.2e1;
t372 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t371 = -Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1;
t370 = t59 / 0.2e1 - t74 / 0.2e1;
t369 = m(7) * t423;
t250 = -m(7) * t432 - t385;
t364 = t154 * t445;
t363 = -t397 / 0.2e1;
t359 = t396 / 0.2e1;
t347 = t384 * qJ(3);
t345 = -mrSges(7,3) * pkin(5) + t488;
t344 = t469 / 0.2e1 + t468;
t243 = mrSges(7,1) * t396;
t340 = mrSges(7,2) * t397 - t243;
t140 = -pkin(5) * t399 + t204;
t338 = t140 * t470 - mrSges(7,1) * t399 / 0.2e1 + t398 * t467;
t179 = t385 * t261;
t180 = t337 * t261;
t292 = -(-mrSges(5,1) / 0.2e1 - t351) * t220 + (pkin(4) * t220 + (t404 - t405) * pkin(9)) * t472 + (t176 * t266 - t177 * t269 - t220 * t277) * t470 + t221 * mrSges(5,2) / 0.2e1 + t489 * (t404 / 0.2e1 - t405 / 0.2e1);
t307 = t101 * t119 + t102 * t120 + t155 * t203;
t308 = t119 * t62 + t120 * t79 + t139 * t155;
t313 = (-t263 * mrSges(5,1) - t248) * t394;
t324 = t289 * t97 - t290 * t98;
t314 = t204 + t324;
t315 = t140 + t325;
t2 = t313 / 0.2e1 + (t459 + t182 / 0.2e1) * t155 + t355 * t120 + t353 * t119 + t307 * t473 + t308 * t471 + t292 + (t179 / 0.2e1 + t180 / 0.2e1 + t354 * t290 + t352 * t289 + t314 * t473 + t315 * t471) * t154;
t304 = -t499 * t261 / 0.2e1 + (t372 - t426 / 0.2e1) * t263;
t3 = t101 * t195 + t102 * t191 - t139 * t179 - t140 * t181 - t203 * t180 - t204 * t182 + t75 * t188 + t98 * t189 + t79 * t190 + t59 * t192 + t97 * t193 + t62 * t194 - t276 * t248 + m(6) * (t101 * t97 + t102 * t98 + t203 * t204) + m(7) * (t139 * t140 + t59 * t62 + t75 * t79) + (Ifges(5,4) * t261 + (-t128 / 0.2e1 - t130 / 0.2e1 + t373 * t261) * t290 + (t124 / 0.2e1 + t126 / 0.2e1 - t372 * t261) * t289) * t261 + (-t276 * mrSges(5,1) - Ifges(5,4) * t263 + t475 * t290 + t304 * t289 + (-Ifges(6,3) - Ifges(7,3) + Ifges(5,1) - Ifges(5,2)) * t261) * t263;
t327 = -t2 * qJD(1) + t3 * qJD(2);
t294 = t444 * t479 + t445 * t480;
t6 = -t139 * t243 + t74 * t190 + t97 * t191 - t98 * t195 + (t369 - t194) * t75 + (t203 * t268 + mrSges(7,2) * t408 + t294 * t263 + (t181 - t438) * t431 - t325 * mrSges(7,3) - t324 * mrSges(6,3) + (t289 * t488 + t502) * t493 + t484 * t445 + t485 * t443) * t263;
t318 = mrSges(7,1) / 0.2e1 + t344;
t300 = t176 * t318 - t177 * t375;
t316 = -t369 / 0.2e1 + t194 / 0.2e1;
t7 = t243 * t461 + (-t289 * t375 + t290 * t344) * t406 + (t195 / 0.2e1 + t374 * t396 + t316) * t120 + (-t374 * t397 - t354) * t119 + t300;
t326 = -qJD(1) * t7 + qJD(2) * t6;
t11 = t503 + t476 * t263 + (t289 * t500 + t290 * t501 + t415) * t261 + m(7) * (-t139 * t263 + t261 * t325) + m(6) * (t261 * t324 - t400) + m(5) * (-t204 * t261 - t400) + m(4) * t347;
t295 = m(5) * (-t155 * t261 - t406) / 0.2e1 + m(4) * t481 / 0.2e1 + 0.2e1 * t381 * (t261 * t321 - t406);
t16 = t295 - t496;
t323 = qJD(1) * t16 + qJD(2) * t11;
t29 = (m(7) * (t289 * t75 + t411) + t289 * t190 + t290 * t194) * t263;
t309 = (t119 * t290 + t120 * t289) * t436;
t43 = -t376 - t309 / 0.2e1;
t322 = -qJD(1) * t43 + qJD(2) * t29;
t118 = m(7) * t380 - t340;
t317 = qJD(2) * t118 + qJD(4) * t250;
t311 = t263 * t268;
t12 = (t416 / 0.2e1 - t354) * t290 + (t417 / 0.2e1 + (t433 / 0.2e1 + t370) * m(7) - t352) * t289 + t386;
t310 = t12 * qJD(2);
t306 = -t139 * t385 / 0.2e1 - t266 * t190 / 0.2e1;
t160 = m(7) * t482 + mrSges(7,3) * t383;
t26 = t338 + t507;
t46 = 0.2e1 * (-t392 / 0.4e1 + t390 / 0.4e1 - t155 / 0.4e1) * m(7);
t305 = qJD(1) * t46 - qJD(2) * t26 + qJD(4) * t160;
t37 = -t277 * t385 - t267 * t432 + (-t427 * t290 + t463) * t290 + (t464 - pkin(5) * t435 + t427 * t289 + (t504 - t505) * t290) * t289;
t296 = (m(7) * t465 + t455) * pkin(5) + t101 * t468 - t102 * mrSges(6,2) / 0.2e1 + mrSges(7,1) * t465 + t79 * t467;
t4 = t243 * t448 - t316 * t269 + (pkin(9) * t339 + t371) * t263 + (t430 / 0.2e1 + mrSges(6,1) * t452 + t126 / 0.4e1 + t124 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(6,6) + 0.3e1 / 0.4e1 * Ifges(7,6)) * t261 + (t459 - t438 / 0.2e1) * pkin(5) + (t463 / 0.2e1 + t413 / 0.2e1 + t277 * t467 + (-Ifges(7,1) / 0.2e1 + Ifges(6,2) / 0.4e1 + Ifges(7,2) / 0.4e1 - Ifges(6,1) / 0.2e1) * t289) * t263) * t289 + (-t130 / 0.4e1 - t128 / 0.4e1 + mrSges(6,2) * t452 + t429 / 0.2e1 + (-0.3e1 / 0.4e1 * Ifges(6,5) - 0.3e1 / 0.4e1 * Ifges(7,5)) * t261 + t370 * mrSges(7,3) + (t412 / 0.2e1 - t464 / 0.2e1 + (Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.4e1 - Ifges(6,2) / 0.2e1 - Ifges(7,2) / 0.2e1) * t290 + (-0.3e1 / 0.2e1 * Ifges(6,4) - 0.3e1 / 0.2e1 * Ifges(7,4)) * t289 + (t435 / 0.2e1 + t267 / 0.2e1) * pkin(5)) * t263) * t290 + t296 + t306;
t301 = t4 * qJD(2) + t37 * qJD(4);
t110 = (-0.1e1 / 0.2e1 + t490) * t436;
t45 = m(7) * t460 + t321 * t471;
t44 = t309 / 0.2e1 - t376;
t27 = t338 - t507;
t18 = t351 * t263 + t293 - t298 + t483 * t445 + (t192 + t193) * t443;
t15 = t295 + t496;
t13 = t369 * t445 + t195 * t447 + (t422 / 0.2e1 + t344 * t289) * t261 + t386 - t506;
t10 = t364 * t469 + (t289 * t318 + t290 * t375) * t154 + (t337 + t385) * t461;
t8 = -t154 * t380 * t470 + t300 + (t311 + t340) * t461 + (t423 * t470 + t352) * t120 + t354 * t119 + t489 * (t119 * t363 + t120 * t359);
t5 = -pkin(4) * t311 / 0.2e1 + t296 - t306 + t203 * t337 / 0.2e1 + t74 * mrSges(7,3) * t443 + (-t289 * t372 + t290 * t373) * t261 + t371 * t263 + t363 * t413 + t429 * t444 + t430 * t447 + t340 * t448 + t380 * t449 - t269 * t454 + t432 * t458 + t411 * t466 + (-t423 * t269 + (-t277 * t396 + t408) * pkin(5)) * t470 + (-t426 * t289 + t488 * t290) * t261 / 0.4e1 - t485 * t289 / 0.4e1 + t484 * t290 / 0.4e1 - t498 * t396 / 0.4e1 + mrSges(6,3) * t428 * t490 + (-t412 + t480) * t359 + (t479 / 0.2e1 + t499 / 0.4e1) * t397;
t1 = t292 + t307 * t472 + t308 * t470 + t195 * t364 - t313 / 0.2e1 + t155 * t458 - t182 * t460 + t483 * t120 / 0.2e1 + (-t179 - t180) * t461 + (t455 + t193 / 0.2e1) * t119 + (t314 * t472 + t315 * t470 + t506) * t154;
t9 = [qJD(2) * t20 + qJD(4) * t19, t15 * qJD(3) + t1 * qJD(4) + t8 * qJD(5) + t44 * qJD(6) + t402 + (-m(6) * t401 + t221 * t415 + m(5) * (-t204 * t221 - t401) + m(4) * (-pkin(2) * t440 + t291 * t347) * t287 + (m(5) * t276 - mrSges(4,1) * t288 + mrSges(5,1) * t261 + mrSges(4,2) * t286 - mrSges(5,2) * t263 - mrSges(3,1)) * t365 - (-t476 + t438) * t220 + (m(6) * t98 + m(7) * t75 - t501) * t177 + (m(6) * t97 + m(7) * t59 + t500) * t176 + (-mrSges(3,2) + t503) * t394) * qJD(2), qJD(2) * t15, t403 + t1 * qJD(2) + t10 * qJD(5) + t45 * qJD(6) + ((-pkin(4) * t155 - t154 * t349) * t472 + (-t154 * t482 + t155 * t277) * t470) * t474 + ((t267 + t486) * t155 + (-t383 * t489 + mrSges(5,2)) * t154) * qJD(4), t8 * qJD(2) + t10 * qJD(4) + ((-mrSges(6,2) - mrSges(7,2)) * t119 + (-mrSges(6,1) + t487) * t120) * qJD(5), qJD(2) * t44 + qJD(4) * t45; qJD(3) * t16 - qJD(4) * t2 - qJD(5) * t7 - qJD(6) * t43 - t402, qJD(3) * t11 + qJD(4) * t3 + qJD(5) * t6 + qJD(6) * t29, qJD(4) * t18 + qJD(5) * t13 + qJD(6) * t110 + t323, t18 * qJD(3) + t5 * qJD(5) + t27 * qJD(6) + ((-pkin(4) * t204 + (-t101 * t289 + t102 * t290) * pkin(9)) * t472 + (t140 * t277 + t266 * t62 - t269 * t79) * t470) * t474 + t327 + (t203 * mrSges(5,2) + Ifges(5,6) * t263 + pkin(4) * t180 + t140 * t267 - t277 * t179 - t269 * t188 + t266 * t192 + (t102 * mrSges(6,3) + t79 * mrSges(7,3) + pkin(9) * t189 + t304) * t290 + (-mrSges(6,3) * t101 - mrSges(7,3) * t62 - pkin(9) * t193 - t475) * t289 + (-Ifges(5,5) + t294) * t261 + t486 * t204) * qJD(4), t13 * qJD(3) + t5 * qJD(4) + t326 + (-mrSges(6,1) * t98 - mrSges(6,2) * t97 - mrSges(7,2) * t74 + (t289 * t345 + t502) * t263 + t487 * t75) * qJD(5), qJD(3) * t110 + qJD(4) * t27 + t322; -qJD(2) * t16, -qJD(4) * t14 - qJD(5) * t12 + qJD(6) * t111 - t323, 0, -t409 (-t337 + t250) * qJD(5) - t310, t382; qJD(2) * t2 + qJD(6) * t46 - t403, qJD(3) * t14 - qJD(5) * t4 - qJD(6) * t26 - t327, t409, -qJD(5) * t37 + qJD(6) * t160, -t301 + (-mrSges(7,2) * t266 + (mrSges(6,2) * pkin(9) - t426) * t289 + (-mrSges(6,1) * pkin(9) + t345) * t290 - t487 * t269) * qJD(5), t305; qJD(2) * t7, qJD(3) * t12 + qJD(4) * t4 + qJD(6) * t118 - t326, t310, qJD(6) * t250 + t301, 0, t317; qJD(2) * t43 - qJD(4) * t46, -qJD(3) * t111 + qJD(4) * t26 - qJD(5) * t118 - t322, -t382, -qJD(5) * t250 - t305, -t317, 0;];
Cq  = t9;
