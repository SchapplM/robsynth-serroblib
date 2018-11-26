% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:28:28
% EndTime: 2018-11-23 18:28:43
% DurationCPUTime: 14.77s
% Computational Cost: add. (18468->712), mult. (44218->942), div. (0->0), fcn. (31825->8), ass. (0->323)
t516 = Ifges(6,1) + Ifges(7,1);
t514 = Ifges(7,4) + Ifges(6,5);
t348 = sin(qJ(3));
t349 = sin(qJ(2));
t352 = cos(qJ(3));
t353 = cos(qJ(2));
t314 = t348 * t353 + t352 * t349;
t301 = t314 * qJD(1);
t345 = qJD(2) + qJD(3);
t347 = sin(qJ(4));
t351 = cos(qJ(4));
t277 = t301 * t351 + t345 * t347;
t346 = sin(qJ(5));
t350 = cos(qJ(5));
t364 = t301 * t347 - t345 * t351;
t200 = t346 * t277 + t350 * t364;
t402 = qJD(1) * t353;
t403 = qJD(1) * t349;
t300 = -t348 * t403 + t352 * t402;
t296 = qJD(4) - t300;
t288 = qJD(5) + t296;
t357 = t350 * t277 - t346 * t364;
t433 = Ifges(6,4) * t357;
t102 = -Ifges(6,2) * t200 + Ifges(6,6) * t288 + t433;
t481 = -pkin(8) - pkin(7);
t330 = t481 * t353;
t318 = qJD(1) * t330;
t302 = t348 * t318;
t328 = t481 * t349;
t317 = qJD(1) * t328;
t307 = qJD(2) * pkin(2) + t317;
t263 = t352 * t307 + t302;
t234 = -t345 * pkin(3) - t263;
t189 = pkin(4) * t364 + t234;
t340 = -pkin(2) * t353 - pkin(1);
t326 = qJD(1) * t340;
t225 = -t300 * pkin(3) - t301 * pkin(9) + t326;
t407 = t352 * t318;
t264 = t348 * t307 - t407;
t235 = pkin(9) * t345 + t264;
t146 = t351 * t225 - t235 * t347;
t123 = -pkin(10) * t277 + t146;
t109 = pkin(4) * t296 + t123;
t147 = t347 * t225 + t351 * t235;
t124 = -pkin(10) * t364 + t147;
t409 = t350 * t124;
t37 = t109 * t346 + t409;
t35 = qJ(6) * t288 + t37;
t73 = t200 * pkin(5) - qJ(6) * t357 + t189;
t197 = Ifges(7,5) * t357;
t99 = Ifges(7,6) * t288 + Ifges(7,3) * t200 + t197;
t551 = t189 * mrSges(6,1) + t73 * mrSges(7,1) - t35 * mrSges(7,2) - t37 * mrSges(6,3) - t102 / 0.2e1 + t99 / 0.2e1;
t515 = -Ifges(6,4) + Ifges(7,5);
t513 = Ifges(7,2) + Ifges(6,3);
t512 = -Ifges(6,6) + Ifges(7,6);
t258 = pkin(3) * t301 - pkin(9) * t300;
t171 = t351 * t258 - t263 * t347;
t416 = t300 * t351;
t378 = t301 * pkin(4) - pkin(10) * t416;
t480 = -pkin(10) - pkin(9);
t388 = qJD(4) * t480;
t550 = t351 * t388 - t171 - t378;
t172 = t347 * t258 + t351 * t263;
t417 = t300 * t347;
t396 = pkin(10) * t417;
t549 = -t347 * t388 + t172 - t396;
t423 = t124 * t346;
t36 = t109 * t350 - t423;
t504 = qJD(6) - t36;
t34 = -pkin(5) * t288 + t504;
t548 = t189 * mrSges(6,2) + t34 * mrSges(7,2) - t36 * mrSges(6,3) - t73 * mrSges(7,3);
t460 = t288 / 0.2e1;
t469 = t357 / 0.2e1;
t472 = t200 / 0.2e1;
t473 = -t200 / 0.2e1;
t547 = -Ifges(6,2) * t473 + Ifges(7,3) * t472 + t460 * t512 + t469 * t515 + t551;
t461 = -t288 / 0.2e1;
t470 = -t357 / 0.2e1;
t546 = -Ifges(6,2) * t472 + Ifges(7,3) * t473 + t461 * t512 + t470 * t515 - t551;
t198 = Ifges(6,4) * t200;
t432 = Ifges(7,5) * t200;
t503 = t514 * t288 + t357 * t516 - t198 + t432;
t287 = pkin(4) * t417;
t400 = qJD(4) * t347;
t342 = pkin(4) * t400;
t545 = t342 - t287;
t389 = qJD(2) * t481;
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
t466 = t253 / 0.2e1;
t485 = t68 / 0.2e1;
t486 = -t68 / 0.2e1;
t487 = t67 / 0.2e1;
t541 = mrSges(6,1) * t105 + mrSges(7,1) * t14 - mrSges(7,2) * t2 - Ifges(6,4) * t67 / 0.2e1 - Ifges(6,6) * t253 / 0.2e1 + 0.2e1 * Ifges(7,3) * t485 - mrSges(6,3) * t7 + (-t486 + t485) * Ifges(6,2) + (t515 + Ifges(7,5)) * t487 + (t512 + Ifges(7,6)) * t466;
t540 = -Ifges(6,4) * t472 - Ifges(7,5) * t473 - t514 * t461 - t516 * t470 + t548;
t539 = t503 / 0.2e1;
t538 = t514 * t253 + t515 * t68 + t516 * t67;
t537 = Ifges(4,5) * t345;
t535 = t345 * Ifges(4,6);
t327 = t480 * t347;
t344 = t351 * pkin(10);
t329 = pkin(9) * t351 + t344;
t280 = t327 * t346 + t329 * t350;
t507 = -qJD(5) * t280 + t346 * t549 + t350 * t550;
t362 = t350 * t327 - t329 * t346;
t506 = qJD(5) * t362 + t346 * t550 - t350 * t549;
t410 = t346 * t351;
t313 = t347 * t350 + t410;
t223 = t313 * t300;
t411 = t346 * t347;
t311 = -t350 * t351 + t411;
t224 = t311 * t300;
t493 = qJD(4) + qJD(5);
t268 = t493 * t311;
t269 = t493 * t313;
t534 = -qJD(6) * t313 + t545 + (-t224 + t268) * qJ(6) + (-t223 + t269) * pkin(5);
t265 = t317 * t348 - t407;
t431 = pkin(2) * qJD(3);
t533 = t348 * t431 - t265;
t419 = t270 * t347;
t532 = t314 * t399 - t419;
t531 = -t347 * t43 + t351 * t42;
t529 = Ifges(5,5) * t277 - t364 * Ifges(5,6) + Ifges(5,3) * t296 + t200 * t512 + t288 * t513 + t357 * t514;
t117 = pkin(5) * t357 + qJ(6) * t200;
t520 = t364 / 0.2e1;
t511 = -t264 + t534;
t49 = t123 * t350 - t423;
t510 = pkin(4) * t397 + qJD(6) - t49;
t293 = t301 * qJ(6);
t509 = -t293 + t506;
t448 = pkin(5) * t301;
t508 = t448 - t507;
t376 = mrSges(5,1) * t347 + mrSges(5,2) * t351;
t505 = t234 * t376;
t108 = -mrSges(5,1) * t186 + mrSges(5,2) * t185;
t502 = m(5) * t178 + t108;
t501 = t533 + t534;
t262 = t312 * pkin(3) - t314 * pkin(9) + t340;
t281 = t328 * t348 - t330 * t352;
t190 = t351 * t262 - t281 * t347;
t413 = t314 * t351;
t136 = pkin(4) * t312 - pkin(10) * t413 + t190;
t274 = t351 * t281;
t191 = t347 * t262 + t274;
t414 = t314 * t347;
t156 = -pkin(10) * t414 + t191;
t500 = t346 * t136 + t350 * t156;
t499 = Ifges(5,5) * t185 + Ifges(5,6) * t186;
t498 = t533 + t545;
t497 = t352 * t328 + t330 * t348;
t496 = t513 * t253 + t512 * t68 + t514 * t67;
t495 = -t146 * t347 + t147 * t351;
t494 = qJD(4) * (-t146 * t351 - t147 * t347);
t8 = -qJD(5) * t37 + t27 * t350 - t29 * t346;
t184 = pkin(3) * t271 + pkin(9) * t270 + t393;
t321 = t349 * t389;
t322 = t353 * t389;
t205 = qJD(3) * t497 + t321 * t352 + t322 * t348;
t381 = t351 * t184 - t205 * t347;
t418 = t270 * t351;
t33 = pkin(10) * t418 + pkin(4) * t271 + (-t274 + (pkin(10) * t314 - t262) * t347) * qJD(4) + t381;
t69 = t347 * t184 + t351 * t205 + t262 * t399 - t281 * t400;
t40 = -pkin(10) * t532 + t69;
t12 = -qJD(5) * t500 + t33 * t350 - t346 * t40;
t492 = t314 * t493;
t490 = m(4) / 0.2e1;
t479 = pkin(1) * mrSges(3,1);
t478 = pkin(1) * mrSges(3,2);
t475 = t185 / 0.2e1;
t474 = t186 / 0.2e1;
t463 = -t277 / 0.2e1;
t462 = t277 / 0.2e1;
t459 = -t296 / 0.2e1;
t458 = -t300 / 0.2e1;
t456 = t301 / 0.2e1;
t453 = -t347 / 0.2e1;
t452 = t351 / 0.2e1;
t451 = m(4) * t326;
t449 = pkin(2) * t352;
t447 = t42 * mrSges(5,2);
t446 = t43 * mrSges(5,1);
t336 = pkin(2) * t348 + pkin(9);
t444 = -pkin(10) - t336;
t45 = mrSges(6,1) * t253 - mrSges(6,3) * t67;
t46 = -t253 * mrSges(7,1) + t67 * mrSges(7,2);
t443 = t46 - t45;
t44 = -mrSges(7,2) * t68 + mrSges(7,3) * t253;
t47 = -mrSges(6,2) * t253 - mrSges(6,3) * t68;
t442 = t47 + t44;
t441 = mrSges(4,3) * t300;
t440 = mrSges(4,3) * t301;
t439 = mrSges(6,3) * t200;
t438 = mrSges(6,3) * t357;
t437 = Ifges(3,4) * t349;
t436 = Ifges(5,4) * t277;
t435 = Ifges(5,4) * t347;
t434 = Ifges(5,4) * t351;
t430 = t301 * Ifges(4,4);
t427 = Ifges(3,5) * qJD(2);
t426 = Ifges(3,6) * qJD(2);
t425 = qJD(2) * mrSges(3,1);
t424 = qJD(2) * mrSges(3,2);
t420 = t178 * t497;
t230 = pkin(2) * t403 + t258;
t266 = t317 * t352 + t302;
t167 = t351 * t230 - t266 * t347;
t126 = t167 + t378;
t168 = t347 * t230 + t351 * t266;
t135 = -t396 + t168;
t59 = t346 * t126 + t350 * t135;
t159 = -mrSges(7,2) * t200 + mrSges(7,3) * t288;
t160 = -mrSges(6,2) * t288 - t439;
t406 = t159 + t160;
t161 = mrSges(6,1) * t288 - t438;
t162 = -mrSges(7,1) * t288 + mrSges(7,2) * t357;
t405 = t161 - t162;
t404 = -mrSges(4,1) * t345 + mrSges(5,1) * t364 + t277 * mrSges(5,2) + t440;
t392 = t352 * t431;
t390 = Ifges(5,3) * t253 + t499;
t339 = -pkin(4) * t351 - pkin(3);
t386 = t427 / 0.2e1;
t385 = -t426 / 0.2e1;
t174 = -t364 * Ifges(5,2) + Ifges(5,6) * t296 + t436;
t384 = t174 * t453;
t276 = Ifges(5,4) * t364;
t175 = t277 * Ifges(5,1) + t296 * Ifges(5,5) - t276;
t383 = t175 * t452;
t382 = qJD(4) * t444;
t226 = pkin(4) * t414 - t497;
t377 = mrSges(5,1) * t351 - mrSges(5,2) * t347;
t375 = Ifges(5,1) * t351 - t435;
t374 = Ifges(5,1) * t347 + t434;
t373 = -Ifges(5,2) * t347 + t434;
t372 = Ifges(5,2) * t351 + t435;
t371 = Ifges(5,5) * t351 - Ifges(5,6) * t347;
t370 = Ifges(5,5) * t347 + Ifges(5,6) * t351;
t58 = t126 * t350 - t135 * t346;
t129 = mrSges(5,1) * t253 - mrSges(5,3) * t185;
t130 = -mrSges(5,2) * t253 + mrSges(5,3) * t186;
t367 = -t347 * t129 + t351 * t130;
t78 = t136 * t350 - t156 * t346;
t309 = t444 * t347;
t310 = t336 * t351 + t344;
t363 = t350 * t309 - t310 * t346;
t255 = t309 * t346 + t310 * t350;
t11 = t136 * t397 - t156 * t398 + t346 * t33 + t350 * t40;
t259 = pkin(5) * t311 - qJ(6) * t313 + t339;
t358 = -t347 * t392 + t351 * t382;
t206 = qJD(3) * t281 + t321 * t348 - t352 * t322;
t4 = -pkin(5) * t253 - t8;
t356 = t8 * mrSges(6,1) - t4 * mrSges(7,1) - t7 * mrSges(6,2) + t2 * mrSges(7,3) + t496;
t128 = pkin(4) * t532 + t206;
t355 = m(5) * (t494 + t531);
t231 = t300 * Ifges(4,2) + t430 + t535;
t294 = Ifges(4,4) * t300;
t232 = Ifges(4,1) * t301 + t294 + t537;
t82 = t185 * Ifges(5,4) + t186 * Ifges(5,2) + t253 * Ifges(5,6);
t83 = Ifges(5,1) * t185 + Ifges(5,4) * t186 + Ifges(5,5) * t253;
t354 = t503 * (-t268 / 0.2e1 + t224 / 0.2e1) + (t294 + t232) * t458 + t547 * t269 - (Ifges(6,4) * t473 + Ifges(7,5) * t472 + t514 * t460 + t516 * t469 + t548) * t268 + t541 * t311 + t540 * t224 - t177 * mrSges(4,2) - t175 * t416 / 0.2e1 + t174 * t417 / 0.2e1 + t546 * t223 + (t538 / 0.2e1 + t105 * mrSges(6,2) + t4 * mrSges(7,2) - t8 * mrSges(6,3) - t14 * mrSges(7,3) + Ifges(6,4) * t486 + Ifges(7,5) * t485 + t514 * t466 + t516 * t487) * t313 + (t383 + t384 + t505) * qJD(4) + (-mrSges(4,1) - t377) * t178 + (t277 * t375 + t296 * t371) * qJD(4) / 0.2e1 - t373 * t359 / 0.2e1 - (Ifges(4,1) * t300 - t430 + t529) * t301 / 0.2e1 + (t146 * t416 + t147 * t417 + t531) * mrSges(5,3) - Ifges(4,5) * t252 - Ifges(4,6) * t253 + (-t505 - t326 * mrSges(4,2) - t537 / 0.2e1 + t371 * t459 + t375 * t463 + t373 * t520) * t300 + (-t36 * mrSges(6,1) + t34 * mrSges(7,1) + Ifges(6,6) * t472 + Ifges(7,6) * t473 + t37 * mrSges(6,2) - t35 * mrSges(7,3) + t147 * mrSges(5,2) - t146 * mrSges(5,1) - t326 * mrSges(4,1) + t535 / 0.2e1 - Ifges(4,2) * t458 + Ifges(5,3) * t459 + Ifges(5,5) * t463 + Ifges(5,6) * t520 + t514 * t470 + t513 * t461) * t301 + t263 * t441 + t82 * t452 + t231 * t456 + t370 * t466 + t372 * t474 + t374 * t475 + t347 * t83 / 0.2e1;
t341 = Ifges(3,4) * t402;
t337 = -pkin(4) * t350 - pkin(5);
t333 = pkin(4) * t346 + qJ(6);
t325 = mrSges(3,3) * t402 - t424;
t324 = -mrSges(3,3) * t403 + t425;
t323 = t339 - t449;
t299 = Ifges(3,1) * t403 + t341 + t427;
t298 = t426 + (t353 * Ifges(3,2) + t437) * qJD(1);
t283 = -mrSges(4,2) * t345 + t441;
t282 = t347 * t382 + t351 * t392;
t257 = -mrSges(4,1) * t300 + mrSges(4,2) * t301;
t241 = t311 * t314;
t240 = t313 * t314;
t237 = t259 - t449;
t218 = mrSges(5,1) * t296 - mrSges(5,3) * t277;
t217 = -t296 * mrSges(5,2) - mrSges(5,3) * t364;
t210 = t287 + t264;
t141 = qJD(5) * t255 + t282 * t346 - t350 * t358;
t140 = qJD(5) * t363 + t350 * t282 + t346 * t358;
t120 = pkin(5) * t240 + qJ(6) * t241 + t226;
t119 = mrSges(6,1) * t200 + mrSges(6,2) * t357;
t118 = mrSges(7,1) * t200 - mrSges(7,3) * t357;
t97 = pkin(4) * t277 + t117;
t96 = -t270 * t410 + (t413 * t493 - t419) * t350 - t411 * t492;
t95 = t311 * t270 - t313 * t492;
t72 = -pkin(5) * t312 - t78;
t71 = qJ(6) * t312 + t500;
t70 = -qJD(4) * t191 + t381;
t51 = -t58 - t448;
t50 = t293 + t59;
t48 = t123 * t346 + t409;
t25 = mrSges(6,1) * t68 + mrSges(6,2) * t67;
t24 = mrSges(7,1) * t68 - mrSges(7,3) * t67;
t19 = pkin(5) * t96 - qJ(6) * t95 + qJD(6) * t241 + t128;
t10 = -pkin(5) * t271 - t12;
t9 = qJ(6) * t271 + qJD(6) * t312 + t11;
t1 = [(t390 + t496 + t499) * t312 / 0.2e1 - t312 * t447 + t547 * t96 + t541 * t240 + (t14 * t241 + t2 * t312 + t271 * t35 - t73 * t95) * mrSges(7,3) + (-t105 * t241 + t189 * t95 - t271 * t37 - t312 * t7) * mrSges(6,2) + (-Ifges(7,5) * t241 + Ifges(7,6) * t312) * t485 + (-Ifges(6,4) * t241 + Ifges(6,6) * t312) * t486 + t8 * (mrSges(6,1) * t312 + mrSges(6,3) * t241) + t4 * (-mrSges(7,1) * t312 - mrSges(7,2) * t241) + (Ifges(7,5) * t95 + Ifges(7,6) * t271) * t472 + (Ifges(6,4) * t95 + Ifges(6,6) * t271) * t473 + (-pkin(7) * t325 - t298 / 0.2e1 + t385 + (-0.2e1 * t479 - 0.3e1 / 0.2e1 * t437 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t353) * qJD(1) + (t257 + qJD(1) * (mrSges(4,1) * t312 + mrSges(4,2) * t314) + 0.2e1 * t451) * pkin(2)) * t401 + m(7) * (t10 * t34 + t120 * t14 + t19 * t73 + t2 * t71 + t35 * t9 + t4 * t72) + (-Ifges(4,4) * t314 - t281 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t312 + t340 * mrSges(4,1)) * t253 + t404 * t206 + t9 * t159 + t11 * t160 + t12 * t161 + t10 * t162 - (mrSges(4,2) * t340 - mrSges(4,3) * t497 - Ifges(4,4) * t312) * t252 - t497 * t108 + (t83 * t452 - Ifges(4,1) * t252 + t375 * t475 + t373 * t474 + t82 * t453 + (mrSges(4,3) + t376) * t178 + (-t347 * t42 - t351 * t43) * mrSges(5,3) + (t175 * t453 - t351 * t174 / 0.2e1 + t234 * t377 + t374 * t463 + t370 * t459 + t372 * t520 - t495 * mrSges(5,3)) * qJD(4)) * t314 + t146 * (mrSges(5,1) * t271 + mrSges(5,3) * t418) + t147 * (-mrSges(5,2) * t271 + mrSges(5,3) * t419) + t296 * (-Ifges(5,5) * t418 + Ifges(5,6) * t419 + Ifges(5,3) * t271) / 0.2e1 - t364 * (-Ifges(5,4) * t418 + Ifges(5,2) * t419 + Ifges(5,6) * t271) / 0.2e1 + t234 * (-mrSges(5,1) * t419 - mrSges(5,2) * t418) + (-Ifges(5,1) * t418 + Ifges(5,4) * t419 + Ifges(5,5) * t271) * t462 + m(4) * (t177 * t281 + t205 * t264 - t206 * t263 - t420) + m(5) * (t146 * t70 + t147 * t69 + t190 * t43 + t191 * t42 + t206 * t234 - t420) + t128 * t119 + t19 * t118 + t120 * t24 + t78 * t45 + t71 * t44 + t72 * t46 + (t299 / 0.2e1 - pkin(7) * t324 + t386 + (-0.2e1 * t478 + 0.3e1 / 0.2e1 * Ifges(3,4) * t353) * qJD(1)) * t353 * qJD(2) + m(6) * (t105 * t226 + t11 * t37 + t12 * t36 + t128 * t189 + t500 * t7 + t78 * t8) + t500 * t47 + (-t177 * t312 + t263 * t270 - t264 * t271) * mrSges(4,3) + t300 * (-Ifges(4,4) * t270 - Ifges(4,2) * t271) / 0.2e1 + t326 * (mrSges(4,1) * t271 - mrSges(4,2) * t270) + t345 * (-Ifges(4,5) * t270 - Ifges(4,6) * t271) / 0.2e1 + (-Ifges(4,1) * t270 - Ifges(4,4) * t271) * t456 + (t271 * t513 + t514 * t95) * t460 + (-t241 * t514 + t312 * t513 + t314 * t371) * t466 + (-t241 * t516 + t312 * t514) * t487 + (t271 * t514 + t516 * t95) * t469 + t190 * t129 + t191 * t130 + t69 * t217 + t70 * t218 + t226 * t25 + t529 * t271 / 0.2e1 - t270 * t383 - t270 * t384 - t270 * t232 / 0.2e1 + t36 * (mrSges(6,1) * t271 - mrSges(6,3) * t95) + t34 * (-mrSges(7,1) * t271 + mrSges(7,2) * t95) - t271 * t231 / 0.2e1 + t205 * t283 - t538 * t241 / 0.2e1 + t95 * t539 + t312 * t446; t501 * t118 + t502 * (-pkin(3) - t449) + t442 * t255 + ((-t299 / 0.2e1 - t341 / 0.2e1 + t386 + qJD(1) * t478 + (t324 - t425) * pkin(7)) * t353 + (t385 + t298 / 0.2e1 + (t479 + t437 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t353) * qJD(1) + (t325 + t424) * pkin(7) + (-t257 - t451) * pkin(2)) * t349) * qJD(1) - m(4) * (-t263 * t265 + t264 * t266) - t404 * t265 - t405 * t141 + t406 * t140 - t50 * t159 - t59 * t160 - t58 * t161 - t51 * t162 + ((t252 * t352 - t253 * t348) * mrSges(4,3) + (t404 * t348 + (t217 * t351 - t218 * t347 + t283) * t352) * qJD(3)) * pkin(2) + 0.2e1 * ((t177 * t348 - t178 * t352) * t490 + (m(5) * (t234 * t348 + t495 * t352) / 0.2e1 + (-t263 * t348 + t264 * t352) * t490) * qJD(3)) * pkin(2) + t498 * t119 - m(5) * (t146 * t167 + t147 * t168 + t234 * t265) + ((-t217 * t347 - t218 * t351) * qJD(4) + t367 + t355) * t336 + (t14 * t237 + t2 * t255 - t363 * t4 + t501 * t73 + (t140 - t50) * t35 + (t141 - t51) * t34) * m(7) - t443 * t363 + (t105 * t323 + t363 * t8 + t255 * t7 + (t140 - t59) * t37 + (-t141 - t58) * t36 + t498 * t189) * m(6) + t354 - t168 * t217 - t167 * t218 + t237 * t24 - t266 * t283 + t323 * t25 + t264 * t440 + mrSges(5,3) * t494; t442 * t280 - t443 * t362 + (-t404 + t440) * t264 + ((-mrSges(5,3) * t146 - pkin(9) * t218) * t351 + (-mrSges(5,3) * t147 + pkin(4) * t119 - pkin(9) * t217) * t347) * qJD(4) + t367 * pkin(9) + t511 * t118 + t508 * t162 + t507 * t161 + t506 * t160 - m(5) * (t146 * t171 + t147 * t172 + t234 * t264) + t509 * t159 + t354 - t210 * t119 - t172 * t217 - t171 * t218 + t259 * t24 - t263 * t283 + t339 * t25 + pkin(9) * t355 - t502 * pkin(3) + (t14 * t259 + t2 * t280 + t34 * t508 + t35 * t509 - t362 * t4 + t511 * t73) * m(7) + (t105 * t339 + t362 * t8 + t280 * t7 + t506 * t37 + t507 * t36 + (-t210 + t342) * t189) * m(6); (m(7) * t34 * t398 - t277 * t119 + t346 * t47 + t350 * t45 + (t160 * t350 - t346 * t405) * qJD(5) + (0.2e1 * t189 * t463 + t346 * t7 + t350 * t8 - t36 * t398 + t37 * t397) * m(6)) * pkin(4) + (t2 * t333 + t337 * t4 - t34 * t48 + t35 * t510 - t73 * t97) * m(7) + t510 * t159 + t356 + t390 + t446 - t447 + (-t146 * t364 + t147 * t277) * mrSges(5,3) - t234 * (t277 * mrSges(5,1) - mrSges(5,2) * t364) - m(6) * (-t36 * t48 + t37 * t49) + t405 * t48 - t49 * t160 - t97 * t118 + t546 * t357 - t146 * t217 + t147 * t218 + t333 * t44 + t337 * t46 + (-Ifges(5,5) * t364 - Ifges(5,6) * t277) * t459 + t174 * t462 + (-Ifges(5,1) * t364 - t436) * t463 + (-Ifges(5,2) * t277 + t175 - t276) * t520 + (t539 + t540) * t200; (t405 + t438) * t37 + (-t406 - t439) * t36 + t356 + qJD(6) * t159 - t117 * t118 - pkin(5) * t46 + qJ(6) * t44 + (t200 * t34 + t35 * t357) * mrSges(7,2) - t73 * (mrSges(7,1) * t357 + mrSges(7,3) * t200) + (Ifges(7,3) * t357 - t432) * t473 + t102 * t469 - t189 * (mrSges(6,1) * t357 - mrSges(6,2) * t200) + (-t200 * t514 + t357 * t512) * t461 + (-pkin(5) * t4 + qJ(6) * t2 - t117 * t73 - t34 * t37 + t35 * t504) * m(7) + (-Ifges(6,2) * t357 - t198 + t503) * t472 + (-t200 * t516 + t197 - t433 + t99) * t470; t357 * t118 - t288 * t159 + 0.2e1 * (t4 / 0.2e1 + t73 * t469 + t35 * t461) * m(7) + t46;];
tauc  = t1(:);
