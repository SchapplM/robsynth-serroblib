% Calculate vector of inverse dynamics joint torques for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR14_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR14_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:24
% EndTime: 2024-09-27 18:42:34
% DurationCPUTime: 6.94s
% Computational Cost: add. (12619->586), mult. (22837->800), div. (0->0), fcn. (16528->26), ass. (0->343)
t364 = sin(qJ(2));
t368 = cos(qJ(3));
t369 = cos(qJ(2));
t360 = cos(pkin(5));
t431 = qJD(3) * t368;
t412 = t360 * t431;
t363 = sin(qJ(3));
t448 = t360 * t363;
t462 = pkin(1) * qJD(1);
t525 = pkin(2) * t412 - (-t364 * t448 + t368 * t369) * t462;
t359 = sin(pkin(5));
t485 = pkin(8) + pkin(9);
t420 = t359 * t485;
t398 = t363 * t420;
t524 = -qJD(3) * t398 + t525;
t447 = t360 * t368;
t378 = pkin(1) * (-t363 * t369 - t364 * t447);
t249 = qJD(1) * t378;
t327 = pkin(2) * t448;
t523 = (-t368 * t420 - t327) * qJD(3) - t249;
t328 = pkin(2) * t447;
t350 = t360 * pkin(3);
t222 = t328 + t350 - t398;
t449 = t359 * t368;
t278 = pkin(8) * t449 + t327;
t325 = pkin(9) * t449;
t240 = t325 + t278;
t362 = sin(qJ(4));
t367 = cos(qJ(4));
t429 = qJD(4) * t367;
t430 = qJD(4) * t362;
t496 = t222 * t429 - t240 * t430 + t523 * t362 + t524 * t367;
t144 = t362 * t222 + t367 * t240;
t495 = -qJD(4) * t144 - t524 * t362 + t523 * t367;
t486 = m(5) * pkin(3);
t522 = -mrSges(4,1) - t486;
t445 = t362 * t363;
t386 = t367 * t368 - t445;
t258 = t386 * t359;
t167 = (qJD(3) + qJD(4)) * t258;
t470 = pkin(10) * t167;
t521 = -t470 + t495;
t443 = t362 * t368;
t387 = t363 * t367 + t443;
t379 = t387 * qJD(4);
t168 = (-t387 * qJD(3) - t379) * t359;
t164 = t168 * pkin(10);
t520 = -t164 - t496;
t356 = qJD(1) + qJD(2);
t226 = t356 * t258;
t354 = qJDD(1) + qJDD(2);
t432 = qJD(3) * t363;
t241 = (t354 * t368 - t356 * t432) * t359;
t242 = (t354 * t363 + t356 * t431) * t359;
t115 = qJD(4) * t226 + t241 * t362 + t242 * t367;
t452 = t356 * t359;
t116 = t241 * t367 - t242 * t362 - t379 * t452;
t361 = sin(qJ(5));
t366 = cos(qJ(5));
t259 = t387 * t359;
t227 = t356 * t259;
t405 = t366 * t226 - t227 * t361;
t43 = qJD(5) * t405 + t115 * t366 + t116 * t361;
t519 = Ifges(6,5) * t43;
t139 = t226 * t361 + t227 * t366;
t44 = -t139 * qJD(5) - t115 * t361 + t116 * t366;
t518 = Ifges(6,6) * t44;
t517 = Ifges(4,5) * t242;
t516 = Ifges(5,5) * t115;
t515 = Ifges(4,6) * t241;
t514 = Ifges(5,6) * t116;
t302 = t360 * t354 + qJDD(3);
t513 = Ifges(4,3) * t302;
t296 = qJDD(4) + t302;
t512 = Ifges(5,3) * t296;
t288 = qJDD(5) + t296;
t511 = Ifges(6,3) * t288;
t131 = Ifges(6,4) * t405;
t304 = t360 * t356 + qJD(3);
t301 = qJD(4) + t304;
t291 = qJD(5) + t301;
t71 = Ifges(6,1) * t139 + Ifges(6,5) * t291 + t131;
t510 = t131 + t71;
t355 = t359 ^ 2;
t451 = t356 * t355;
t509 = -mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t423 = t364 * t462;
t282 = pkin(8) * t452 + t423;
t422 = t369 * t462;
t297 = pkin(2) * t356 + t422;
t418 = t297 * t448;
t178 = t282 * t368 + t418;
t419 = qJD(2) * t462;
t458 = pkin(1) * qJDD(1);
t285 = t364 * t458 + t369 * t419;
t471 = pkin(8) * t359;
t243 = t354 * t471 + t285;
t284 = -t364 * t419 + t369 * t458;
t257 = pkin(2) * t354 + t284;
t108 = -t178 * qJD(3) - t243 * t363 + t257 * t447;
t79 = pkin(3) * t302 - pkin(9) * t242 + t108;
t107 = t368 * t243 + t257 * t448 - t282 * t432 + t297 * t412;
t88 = pkin(9) * t241 + t107;
t266 = t297 * t447;
t395 = pkin(9) * t452 + t282;
t158 = -t395 * t363 + t266;
t140 = pkin(3) * t304 + t158;
t159 = t395 * t368 + t418;
t151 = t367 * t159;
t93 = t140 * t362 + t151;
t21 = -t93 * qJD(4) - t362 * t88 + t367 * t79;
t14 = pkin(4) * t296 - pkin(10) * t115 + t21;
t20 = t140 * t429 - t159 * t430 + t362 * t79 + t367 * t88;
t15 = pkin(10) * t116 + t20;
t469 = pkin(10) * t226;
t66 = t93 + t469;
t460 = t361 * t66;
t221 = t227 * pkin(10);
t149 = t362 * t159;
t92 = t367 * t140 - t149;
t65 = -t221 + t92;
t63 = pkin(4) * t301 + t65;
t28 = t366 * t63 - t460;
t4 = t28 * qJD(5) + t14 * t361 + t15 * t366;
t459 = t366 * t66;
t29 = t361 * t63 + t459;
t5 = -t29 * qJD(5) + t14 * t366 - t15 * t361;
t508 = t5 * mrSges(6,1) - t4 * mrSges(6,2);
t507 = t21 * mrSges(5,1) - t20 * mrSges(5,2);
t506 = t108 * mrSges(4,1) - t107 * mrSges(4,2);
t463 = Ifges(6,4) * t139;
t70 = Ifges(6,2) * t405 + Ifges(6,6) * t291 + t463;
t505 = t70 / 0.2e1;
t484 = -t405 / 0.2e1;
t349 = t360 * pkin(4);
t504 = m(6) * t386 * t349;
t143 = t367 * t222 - t240 * t362;
t407 = -pkin(10) * t259 + t349;
t112 = t143 + t407;
t251 = t258 * pkin(10);
t119 = t251 + t144;
t57 = t112 * t366 - t119 * t361;
t503 = t57 * qJD(5) + t521 * t361 - t520 * t366;
t58 = t112 * t361 + t119 * t366;
t502 = -t58 * qJD(5) + t520 * t361 + t521 * t366;
t352 = t368 * pkin(3);
t339 = t352 + pkin(2);
t501 = mrSges(6,2) * t405;
t500 = Ifges(6,1) * t405;
t499 = Ifges(6,5) * t405;
t474 = pkin(3) * t367;
t338 = pkin(4) + t474;
t427 = qJD(5) * t366;
t428 = qJD(5) * t361;
t444 = t362 * t366;
t100 = -t158 * t362 - t151;
t72 = t100 - t469;
t101 = t367 * t158 - t149;
t73 = -t221 + t101;
t498 = t361 * t73 - t366 * t72 - t338 * t428 + (-t362 * t427 + (-t361 * t367 - t444) * qJD(4)) * pkin(3);
t446 = t361 * t362;
t497 = -t361 * t72 - t366 * t73 + t338 * t427 + (-t362 * t428 + (t366 * t367 - t446) * qJD(4)) * pkin(3);
t494 = t297 * (mrSges(4,1) * t363 + mrSges(4,2) * t368);
t466 = Ifges(4,4) * t363;
t493 = t363 * (Ifges(4,1) * t368 - t466) * t451;
t340 = pkin(1) * t369 + pkin(2);
t293 = t340 * t447;
t303 = pkin(1) * t364 + t471;
t408 = -pkin(9) * t359 - t303;
t176 = t408 * t363 + t293 + t350;
t292 = t340 * t448;
t229 = t368 * t303 + t292;
t196 = t325 + t229;
t118 = t362 * t176 + t367 * t196;
t492 = -t278 * qJD(3) - t249;
t413 = t359 * t432;
t491 = -pkin(8) * t413 + t525;
t357 = qJ(3) + qJ(4);
t344 = pkin(5) + t357;
t322 = cos(t344) / 0.2e1;
t409 = pkin(5) - t357;
t332 = cos(t409);
t394 = sin(t409);
t425 = sin(t344) / 0.2e1;
t334 = qJ(5) + t344;
t315 = cos(t334) / 0.2e1;
t396 = -qJ(5) + t409;
t324 = cos(t396);
t385 = sin(t396);
t426 = sin(t334) / 0.2e1;
t438 = (t426 + t385 / 0.2e1) * mrSges(6,1) + (t315 - t324 / 0.2e1) * mrSges(6,2);
t490 = -(t425 + t394 / 0.2e1) * mrSges(5,1) - (t322 - t332 / 0.2e1) * mrSges(5,2) - t438;
t281 = t332 / 0.2e1 + t322;
t358 = qJ(1) + qJ(2);
t346 = sin(t358);
t345 = sin(t357);
t348 = cos(t358);
t455 = t345 * t348;
t209 = -t281 * t346 - t455;
t280 = t425 - t394 / 0.2e1;
t347 = cos(t357);
t388 = t280 * t346 - t347 * t348;
t270 = t324 / 0.2e1 + t315;
t351 = qJ(5) + t357;
t335 = sin(t351);
t193 = -t270 * t346 - t335 * t348;
t269 = t426 - t385 / 0.2e1;
t336 = cos(t351);
t390 = t269 * t346 - t336 * t348;
t441 = t193 * mrSges(6,1) + t390 * mrSges(6,2);
t489 = -t209 * mrSges(5,1) - t388 * mrSges(5,2) - t441;
t456 = t345 * t346;
t208 = -t281 * t348 + t456;
t389 = -t280 * t348 - t346 * t347;
t192 = -t270 * t348 + t335 * t346;
t391 = -t269 * t348 - t336 * t346;
t442 = -t192 * mrSges(6,1) + t391 * mrSges(6,2);
t488 = t208 * mrSges(5,1) - t389 * mrSges(5,2) - t442;
t487 = m(3) * pkin(1);
t483 = -t139 / 0.2e1;
t482 = t139 / 0.2e1;
t480 = t227 / 0.2e1;
t479 = -t291 / 0.2e1;
t477 = mrSges(6,3) * t28;
t476 = mrSges(6,3) * t29;
t475 = pkin(3) * t363;
t473 = pkin(4) * t227;
t472 = pkin(4) * t258;
t467 = mrSges(5,3) * t226;
t465 = Ifges(4,4) * t368;
t464 = Ifges(5,4) * t227;
t461 = t227 * mrSges(5,3);
t454 = t346 * t359;
t450 = t359 * t363;
t433 = qJD(2) * t369;
t436 = t368 * pkin(1) * t433 + t340 * t412;
t435 = t348 * pkin(2) + pkin(8) * t454;
t434 = qJD(2) * t364;
t424 = t511 + t518 + t519;
t421 = pkin(1) * t434;
t417 = t356 * t450;
t416 = t356 * t449;
t415 = t512 + t514 + t516;
t414 = t513 + t515 + t517;
t411 = t450 / 0.2e1;
t410 = t431 / 0.2e1;
t316 = pkin(3) * t413;
t148 = -pkin(4) * t168 + t316;
t117 = t367 * t176 - t196 * t362;
t337 = pkin(4) * t367 + pkin(3);
t215 = -t359 * (pkin(10) + t485) + (pkin(4) * t443 + t337 * t363) * t360;
t294 = pkin(4) * t347 + t339;
t406 = -t215 * t346 + t348 * t294;
t283 = pkin(3) * t448 - t420;
t404 = -t283 * t346 + t348 * t339;
t403 = pkin(3) * t417;
t402 = mrSges(4,3) * t417;
t401 = mrSges(4,3) * t416;
t400 = t359 * t423;
t399 = t360 * t421;
t290 = t339 * t359;
t267 = (-t340 - t352) * t359;
t157 = -pkin(3) * t241 - t257 * t359;
t102 = t251 + t118;
t98 = t117 + t407;
t47 = t102 * t366 + t361 * t98;
t46 = -t102 * t361 + t366 * t98;
t275 = (-mrSges(4,1) * t368 + mrSges(4,2) * t363) * t359;
t392 = -qJD(3) * t359 * t494 - t257 * t275;
t160 = t258 * t366 - t259 * t361;
t161 = t258 * t361 + t259 * t366;
t384 = t424 + t508;
t225 = (-t356 * t352 - t297) * t359;
t383 = -pkin(4) * t445 + t337 * t368;
t382 = -t346 * t363 + t348 * t447;
t254 = -t346 * t447 - t348 * t363;
t381 = (t368 * Ifges(4,2) + t466) * t359;
t177 = -t282 * t363 + t266;
t380 = (-t177 * t368 - t178 * t363) * mrSges(4,3);
t141 = (t408 * qJD(3) - t399) * t363 + t436;
t375 = qJD(2) * t378;
t142 = t375 + (t408 * t368 - t292) * qJD(3);
t48 = t367 * t141 + t362 * t142 + t176 * t429 - t196 * t430;
t377 = t304 * t359 * (Ifges(4,5) * t368 - Ifges(4,6) * t363);
t49 = -qJD(4) * t118 - t141 * t362 + t367 * t142;
t255 = -t346 * t448 + t348 * t368;
t374 = -t348 * mrSges(3,1) - t255 * mrSges(4,1) + t388 * mrSges(5,1) + t390 * mrSges(6,1) + t346 * mrSges(3,2) - t254 * mrSges(4,2) - t209 * mrSges(5,2) - t193 * mrSges(6,2) + t509 * t454;
t253 = -t346 * t368 - t348 * t448;
t373 = -t253 * mrSges(4,1) - t389 * mrSges(5,1) - t391 * mrSges(6,1) + t382 * mrSges(4,2) - t208 * mrSges(5,2) - t192 * mrSges(6,2) + (m(4) * pkin(2) + m(5) * t339 + m(6) * t294 + mrSges(3,1)) * t346 + (m(5) * t283 + m(6) * t215 + mrSges(3,2) + (-m(4) * pkin(8) + t509) * t359) * t348;
t123 = Ifges(5,2) * t226 + Ifges(5,6) * t301 + t464;
t220 = Ifges(5,4) * t226;
t124 = Ifges(5,1) * t227 + Ifges(5,5) * t301 + t220;
t147 = -pkin(4) * t226 + t225;
t372 = t405 * t477 - t225 * (mrSges(5,1) * t227 + mrSges(5,2) * t226) - t147 * t501 + t93 * t461 + t123 * t480 - t227 * (Ifges(5,1) * t226 - t464) / 0.2e1 - t301 * (Ifges(5,5) * t226 - Ifges(5,6) * t227) / 0.2e1 + t384 + t415 + t500 * t483 + t92 * t467 + t499 * t479 + t510 * t484 - (-Ifges(5,2) * t227 + t124 + t220) * t226 / 0.2e1 + (-t147 * mrSges(6,1) - Ifges(6,4) * t483 - Ifges(6,2) * t484 - Ifges(6,6) * t479 + t476 + t505) * t139 + t507;
t194 = t304 * Ifges(4,6) + t356 * t381;
t300 = Ifges(4,4) * t416;
t195 = Ifges(4,1) * t417 + t304 * Ifges(4,5) + t300;
t67 = t160 * qJD(5) + t167 * t366 + t168 * t361;
t68 = -t161 * qJD(5) - t167 * t361 + t168 * t366;
t74 = -pkin(4) * t116 + t157;
t371 = (-Ifges(4,2) * t363 + t465) * t410 * t451 + t67 * t71 / 0.2e1 + t147 * (-mrSges(6,1) * t68 + mrSges(6,2) * t67) + t167 * t124 / 0.2e1 + t168 * t123 / 0.2e1 + t225 * (-mrSges(5,1) * t168 + mrSges(5,2) * t167) + t226 * (Ifges(5,4) * t167 + Ifges(5,2) * t168) / 0.2e1 + t284 * mrSges(3,1) + t291 * (Ifges(6,5) * t67 + Ifges(6,6) * t68) / 0.2e1 + t301 * (Ifges(5,5) * t167 + Ifges(5,6) * t168) / 0.2e1 + Ifges(3,3) * t354 - t194 * t413 / 0.2e1 + (Ifges(4,4) * t242 + Ifges(4,2) * t241 + Ifges(4,6) * t302) * t449 / 0.2e1 + (t107 * t449 - t108 * t450) * mrSges(4,3) + (t515 / 0.2e1 + t512 / 0.2e1 + t511 / 0.2e1 + t514 / 0.2e1 + t516 / 0.2e1 + t518 / 0.2e1 + t519 / 0.2e1 + t517 / 0.2e1 + t513 / 0.2e1 + t506 + t507 + t508) * t360 + t68 * t505 + t405 * (Ifges(6,4) * t67 + Ifges(6,2) * t68) / 0.2e1 + t68 * t476 + (Ifges(5,1) * t167 + Ifges(5,4) * t168) * t480 + (Ifges(6,1) * t67 + Ifges(6,4) * t68) * t482 + (Ifges(4,1) * t242 + Ifges(4,4) * t241 + Ifges(4,5) * t302) * t411 + (-t157 * mrSges(5,1) + t20 * mrSges(5,3) + Ifges(5,4) * t115 + Ifges(5,2) * t116 + Ifges(5,6) * t296) * t258 + (t242 * (t363 * Ifges(4,1) + t465) / 0.2e1 + t302 * (Ifges(4,5) * t363 + Ifges(4,6) * t368) / 0.2e1 + t195 * t410) * t359 + (t157 * mrSges(5,2) - t21 * mrSges(5,3) + Ifges(5,1) * t115 + Ifges(5,4) * t116 + Ifges(5,5) * t296) * t259 - t67 * t477 + (t493 + t377) * qJD(3) / 0.2e1 + (t424 + t415 + t414) * t360 / 0.2e1 + (-t92 * t167 + t93 * t168) * mrSges(5,3) + (t74 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t43 + Ifges(6,4) * t44 + Ifges(6,5) * t288) * t161 + (-mrSges(6,1) * t74 + mrSges(6,3) * t4 + Ifges(6,4) * t43 + Ifges(6,2) * t44 + Ifges(6,6) * t288) * t160 + t241 * t381 / 0.2e1;
t370 = cos(qJ(1));
t365 = sin(qJ(1));
t353 = t370 * pkin(1);
t317 = t359 * t421;
t298 = -pkin(4) * t345 - t475;
t277 = -pkin(8) * t450 + t328;
t272 = pkin(3) * t444 + t338 * t361;
t271 = -pkin(3) * t446 + t338 * t366;
t263 = t317 + t316;
t239 = t356 * t275;
t238 = t383 * t360;
t237 = -mrSges(4,2) * t304 + t401;
t236 = mrSges(4,1) * t304 - t402;
t228 = -t303 * t363 + t293;
t197 = -t290 - t472;
t180 = mrSges(4,1) * t302 - mrSges(4,3) * t242;
t179 = -mrSges(4,2) * t302 + mrSges(4,3) * t241;
t175 = t267 - t472;
t174 = t403 + t473;
t173 = mrSges(5,1) * t301 - t461;
t172 = -mrSges(5,2) * t301 + t467;
t154 = -qJD(3) * t229 + t375;
t153 = (-qJD(3) * t303 - t399) * t363 + t436;
t152 = -mrSges(4,1) * t241 + mrSges(4,2) * t242;
t146 = -mrSges(5,1) * t226 + mrSges(5,2) * t227;
t129 = t148 + t317;
t121 = mrSges(6,1) * t291 - mrSges(6,3) * t139;
t120 = -mrSges(6,2) * t291 + mrSges(6,3) * t405;
t106 = -mrSges(5,2) * t296 + mrSges(5,3) * t116;
t105 = mrSges(5,1) * t296 - mrSges(5,3) * t115;
t81 = -mrSges(6,1) * t405 + mrSges(6,2) * t139;
t59 = -mrSges(5,1) * t116 + mrSges(5,2) * t115;
t40 = t49 - t470;
t39 = t164 + t48;
t38 = -mrSges(6,2) * t288 + mrSges(6,3) * t44;
t37 = mrSges(6,1) * t288 - mrSges(6,3) * t43;
t32 = t366 * t65 - t460;
t31 = -t361 * t65 - t459;
t13 = -mrSges(6,1) * t44 + mrSges(6,2) * t43;
t7 = -t47 * qJD(5) - t361 * t39 + t366 * t40;
t6 = t46 * qJD(5) + t361 * t40 + t366 * t39;
t1 = [t46 * t37 + t47 * t38 + t117 * t105 + t118 * t106 + t6 * t120 + t7 * t121 + t129 * t81 + t48 * t172 + t49 * t173 + t175 * t13 + t228 * t180 + t229 * t179 + t154 * t236 + t153 * t237 + t263 * t146 + t267 * t59 - t285 * mrSges(3,2) + m(4) * (t107 * t229 + t108 * t228 + t153 * t178 + t154 * t177 + (t257 * t340 - t297 * t421) * t355) + (qJD(3) * t380 - t152 * t340 + t239 * t421 + t392) * t359 + ((-t354 * t364 - t356 * t433) * mrSges(3,2) + (t354 * t369 - t356 * t434) * mrSges(3,1)) * pkin(1) + Ifges(2,3) * qJDD(1) + m(6) * (t129 * t147 + t175 * t74 + t28 * t7 + t29 * t6 + t4 * t47 + t46 * t5) + m(5) * (t117 * t21 + t118 * t20 + t157 * t267 + t225 * t263 + t48 * t93 + t49 * t92) + (t284 * t369 + t285 * t364) * t487 + t371 + ((mrSges(2,1) + (m(3) + m(4) + m(5) + m(6)) * pkin(1)) * t365 + t373 + mrSges(2,2) * t370) * g(1) + (t374 - m(4) * (t353 + t435) - m(6) * (t353 + t406) - m(5) * (t353 + t404) + t365 * mrSges(2,2) + (-mrSges(2,1) - t487) * t370) * g(2); t57 * t37 + t58 * t38 + t143 * t105 + t144 * t106 + t148 * t81 + t197 * t13 + t277 * t180 + t278 * t179 - t290 * t59 + (t356 * t422 - t285) * mrSges(3,2) + t491 * t237 + t492 * t236 + t502 * t121 + t503 * t120 + t374 * g(2) + t371 + t495 * t173 + t496 * t172 + t373 * g(1) + t356 * mrSges(3,1) * t423 + (-pkin(2) * t152 + (t146 * t475 + t380) * qJD(3) + (-t146 - t239 - t81) * t423 + t392) * t359 + (-t406 * g(2) + t197 * t74 + t4 * t58 + t5 * t57 + t503 * t29 + t502 * t28 + (t148 - t400) * t147) * m(6) + (-t404 * g(2) + t143 * t21 + t144 * t20 - t157 * t290 + t496 * t93 + t495 * t92 + (t316 - t400) * t225) * m(5) + (-t435 * g(2) + t107 * t278 + t108 * t277 + (pkin(2) * t257 + t297 * t423) * t355 + t491 * t178 + t492 * t177) * m(4); -t101 * t172 - t100 * t173 - t174 * t81 + t271 * t37 + t272 * t38 - t146 * t403 - m(5) * (t100 * t92 + t101 * t93 + t225 * t403) + t414 + t451 * t494 - (-Ifges(4,2) * t417 + t195 + t300) * t416 / 0.2e1 + t105 * t474 + (t20 * t362 + t21 * t367 + (-t362 * t92 + t367 * t93) * qJD(4)) * t486 + t372 + t497 * t120 + t498 * t121 + (-t147 * t174 + t271 * t5 + t272 * t4 + t498 * t28 + t497 * t29) * m(6) + (-t377 / 0.2e1 + t194 * t411 - t493 / 0.2e1) * t356 + (t255 * mrSges(4,2) - m(6) * (-t238 * t346 + t298 * t348) + t522 * t254 + t489) * g(1) + (-m(6) * t383 * t359 - t449 * t486 + t275 + t490) * g(3) + (t236 + t402) * t178 + (-t237 + t401) * t177 + (t362 * t106 + t172 * t429 - t173 * t430) * pkin(3) + t506 + (-t253 * mrSges(4,2) - m(6) * (t238 * t348 + t298 * t346) + t488 + t522 * t382) * g(2); -t81 * t473 - t32 * t120 - t31 * t121 - t92 * t172 + t93 * t173 + t372 - m(6) * (t147 * t473 + t28 * t31 + t29 * t32) + (-m(6) * t472 + t490) * g(3) + (-t348 * t504 + t488) * g(2) + (t346 * t504 + t489) * g(1) + (t120 * t427 - t121 * t428 + t361 * t38 + t366 * t37 + (g(1) * t455 + t361 * t4 + t366 * t5 + (-t28 * t361 + t29 * t366) * qJD(5) + g(2) * t456) * m(6)) * pkin(4); -t147 * (mrSges(6,1) * t139 + t501) + (-t463 + t500) * t483 + t70 * t482 + (-Ifges(6,6) * t139 + t499) * t479 - t28 * t120 + t29 * t121 - g(1) * t441 - g(2) * t442 - g(3) * t438 + (t139 * t29 + t28 * t405) * mrSges(6,3) + t384 + (-Ifges(6,2) * t139 + t510) * t484;];
tau = t1;
