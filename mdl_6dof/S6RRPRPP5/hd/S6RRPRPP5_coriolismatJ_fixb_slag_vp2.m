% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:34
% EndTime: 2019-03-09 10:03:47
% DurationCPUTime: 6.66s
% Computational Cost: add. (6728->585), mult. (13103->733), div. (0->0), fcn. (10130->4), ass. (0->274)
t311 = cos(qJ(2));
t293 = t311 * qJ(5);
t310 = cos(qJ(4));
t280 = t310 * t293;
t308 = sin(qJ(4));
t450 = pkin(4) + pkin(5);
t484 = t450 * t308;
t381 = t311 * t484;
t490 = t280 - t381;
t292 = t310 * qJ(5);
t489 = t292 - t484;
t294 = qJ(5) * t308;
t483 = t450 * t310;
t488 = t294 + t483;
t306 = t308 ^ 2;
t307 = t310 ^ 2;
t394 = t306 + t307;
t369 = -pkin(4) * t308 + t292;
t458 = m(7) / 0.2e1;
t460 = m(6) / 0.2e1;
t487 = t369 * t460 + t458 * t489;
t451 = pkin(3) + pkin(7);
t486 = mrSges(6,2) + mrSges(5,3);
t475 = Ifges(6,4) + Ifges(5,5);
t485 = Ifges(5,6) + Ifges(7,6);
t309 = sin(qJ(2));
t404 = t308 * t309;
t386 = mrSges(7,3) * t404;
t410 = t311 * mrSges(7,1);
t208 = -t386 - t410;
t481 = t208 / 0.2e1;
t480 = -t309 / 0.2e1;
t479 = t309 / 0.2e1;
t478 = -t311 / 0.2e1;
t437 = t311 / 0.2e1;
t313 = -pkin(2) - pkin(8);
t370 = -qJ(3) * t309 - pkin(1);
t193 = t313 * t311 + t370;
t248 = t451 * t309;
t396 = t308 * t193 - t310 * t248;
t403 = t308 * t311;
t334 = qJ(6) * t403 + t396;
t46 = -t450 * t309 + t334;
t399 = t310 * t311;
t279 = qJ(6) * t399;
t291 = t309 * qJ(5);
t84 = t193 * t310 + t248 * t308;
t71 = t291 + t84;
t53 = t279 + t71;
t477 = m(7) * (t46 * t308 + t53 * t310);
t476 = -mrSges(5,1) - mrSges(6,1);
t474 = -Ifges(4,6) - Ifges(3,4);
t401 = t309 * t310;
t215 = -t311 * mrSges(5,2) + mrSges(5,3) * t401;
t219 = mrSges(6,2) * t401 + t311 * mrSges(6,3);
t472 = t215 + t219;
t217 = -mrSges(5,2) * t309 - mrSges(5,3) * t399;
t295 = t309 * mrSges(6,3);
t385 = mrSges(6,2) * t399;
t218 = t295 - t385;
t395 = t217 + t218;
t297 = t310 * mrSges(6,3);
t354 = -t308 * mrSges(6,1) + t297;
t300 = Ifges(6,5) * t310;
t244 = -Ifges(6,1) * t308 + t300;
t301 = Ifges(7,4) * t310;
t242 = -Ifges(7,1) * t308 + t301;
t469 = -Ifges(6,3) * t308 - t300;
t468 = -Ifges(7,2) * t308 - t301;
t467 = Ifges(6,2) + Ifges(5,3) + Ifges(7,3);
t431 = mrSges(6,2) - mrSges(7,3);
t466 = -qJ(5) * t431 - Ifges(7,6);
t465 = t311 * t394;
t132 = -t309 * Ifges(7,5) + t311 * t242;
t134 = t309 * Ifges(6,4) + t311 * t244;
t423 = Ifges(5,4) * t310;
t352 = Ifges(5,1) * t308 + t423;
t136 = t309 * Ifges(5,5) - t311 * t352;
t363 = -Ifges(5,5) / 0.2e1 + Ifges(7,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t464 = -t363 * t309 + t132 / 0.2e1 + t134 / 0.2e1 + t136 / 0.2e1;
t463 = 0.2e1 * m(7);
t462 = 2 * qJD(4);
t461 = -m(6) / 0.2e1;
t459 = -m(7) / 0.2e1;
t457 = -mrSges(7,3) / 0.2e1;
t456 = -Ifges(6,6) / 0.2e1;
t368 = t309 * pkin(2) - qJ(3) * t311;
t202 = pkin(8) * t309 + t368;
t249 = t451 * t311;
t86 = -t308 * t202 + t249 * t310;
t47 = -qJ(6) * t404 - t450 * t311 - t86;
t455 = t47 / 0.2e1;
t76 = -pkin(4) * t311 - t86;
t454 = t76 / 0.2e1;
t453 = t84 / 0.2e1;
t452 = m(6) + m(7);
t105 = pkin(4) * t399 + t308 * t293 + t249;
t448 = t105 / 0.2e1;
t447 = -t490 / 0.2e1;
t413 = t310 * mrSges(7,1);
t415 = t308 * mrSges(7,2);
t353 = t413 + t415;
t182 = t353 * t311;
t446 = -t182 / 0.2e1;
t281 = mrSges(6,2) * t404;
t411 = t311 * mrSges(6,1);
t210 = t281 - t411;
t445 = t210 / 0.2e1;
t414 = t309 * mrSges(7,1);
t211 = mrSges(7,3) * t403 - t414;
t444 = t211 / 0.2e1;
t283 = mrSges(7,3) * t399;
t296 = t309 * mrSges(7,2);
t216 = t283 + t296;
t443 = -t216 / 0.2e1;
t442 = -t354 / 0.2e1;
t441 = t249 / 0.2e1;
t440 = t308 / 0.2e1;
t439 = -t310 / 0.2e1;
t438 = t310 / 0.2e1;
t226 = qJ(3) - t369;
t436 = m(6) * t226;
t435 = m(7) * t488;
t434 = m(7) * t310;
t432 = pkin(4) * t310;
t430 = t46 - t334;
t62 = t279 + t84;
t429 = t53 - t62;
t428 = t71 - t84;
t72 = -pkin(4) * t309 + t396;
t427 = t72 - t396;
t424 = Ifges(5,4) * t308;
t422 = Ifges(7,4) * t308;
t421 = Ifges(6,5) * t308;
t212 = mrSges(5,1) * t309 + mrSges(5,3) * t403;
t213 = -mrSges(6,1) * t309 - mrSges(6,2) * t403;
t227 = -pkin(2) * t311 + t370;
t371 = m(4) * t227 + t311 * mrSges(4,2) - mrSges(4,3) * t309;
t9 = -t212 * t404 + ((t211 + t213) * t308 + (t216 + t395) * t310 + t477 - m(6) * (-t72 * t308 - t71 * t310) - m(5) * (-t308 * t396 - t84 * t310) + t371) * t309;
t418 = qJD(1) * t9;
t367 = t294 + t451;
t104 = (-t367 - t432) * t309;
t228 = t310 * mrSges(6,1) + t308 * mrSges(6,3);
t178 = t228 * t309;
t179 = t353 * t309;
t230 = t310 * mrSges(5,1) - t308 * mrSges(5,2);
t180 = t230 * t309;
t181 = t228 * t311;
t209 = t311 * mrSges(5,1) - mrSges(5,3) * t404;
t214 = t311 * mrSges(7,2) - mrSges(7,3) * t401;
t328 = pkin(5) * t399 + t105;
t330 = t248 * mrSges(5,2) + Ifges(7,5) * t437 + t475 * t478 + (-t242 - t244 + t352) * t480;
t347 = -Ifges(6,3) * t310 + t421;
t348 = -Ifges(7,2) * t310 + t422;
t349 = Ifges(5,2) * t310 + t424;
t331 = -t248 * mrSges(5,1) + Ifges(6,6) * t437 + t349 * t480 + (t347 + t348) * t479 + t485 * t478;
t126 = t309 * Ifges(6,6) - t311 * t347;
t128 = -t309 * Ifges(7,6) - t311 * t348;
t130 = Ifges(5,6) * t309 - t311 * t349;
t359 = -t126 / 0.2e1 - t128 / 0.2e1 + t130 / 0.2e1;
t383 = Ifges(7,6) / 0.2e1 + t456;
t362 = Ifges(5,6) / 0.2e1 + t383;
t87 = t310 * t202 + t308 * t249;
t74 = t293 + t87;
t55 = -qJ(6) * t401 + t74;
t88 = (t367 + t483) * t309;
t3 = -t105 * t178 - t328 * t179 + t104 * t181 - t88 * t182 + t46 * t208 - t396 * t209 + t72 * t210 + t47 * t211 + t86 * t212 + t76 * t213 + t53 * t214 + t84 * t215 + t55 * t216 + t87 * t217 + t74 * t218 + t71 * t219 - t249 * t180 + t371 * t368 + m(5) * (-t248 * t249 - t396 * t86 + t84 * t87) + m(6) * (t104 * t105 + t71 * t74 + t72 * t76) + m(7) * (-t328 * t88 + t46 * t47 + t53 * t55) + (-pkin(1) * mrSges(3,1) - t227 * mrSges(4,2) + t474 * t309 + (t309 * t362 + t359) * t310 + t464 * t308) * t309 + (-pkin(1) * mrSges(3,2) - t227 * mrSges(4,3) + t331 * t310 + t330 * t308 + (-Ifges(3,2) + Ifges(3,1) + Ifges(4,2) - Ifges(4,3) + t467) * t309 + (t363 * t308 - t362 * t310 - t474) * t311) * t311;
t417 = t3 * qJD(1);
t412 = t310 * mrSges(5,2);
t183 = -pkin(4) * t403 + t280;
t184 = t354 * t311;
t185 = mrSges(7,1) * t403 - mrSges(7,2) * t399;
t355 = -t308 * mrSges(5,1) - t412;
t186 = t355 * t311;
t187 = t469 * t311;
t188 = t468 * t311;
t241 = -Ifges(5,2) * t308 + t423;
t189 = t311 * t241;
t243 = Ifges(7,1) * t310 + t422;
t190 = t311 * t243;
t245 = Ifges(6,1) * t310 + t421;
t191 = t311 * t245;
t247 = Ifges(5,1) * t310 - t424;
t192 = t311 * t247;
t284 = Ifges(5,6) * t403;
t4 = t490 * t182 + t183 * t181 + t105 * t184 - t328 * t185 + t62 * t211 - t334 * t216 + t249 * t186 + t284 * t479 + (-t212 + t213) * t84 - t395 * t396 + m(6) * (t105 * t183 - t396 * t71 + t72 * t84) + m(7) * (t328 * t490 - t334 * t53 + t46 * t62) + ((t71 * mrSges(6,2) + t84 * mrSges(5,3) - t53 * mrSges(7,3) + t190 / 0.2e1 + t191 / 0.2e1 + t192 / 0.2e1 + t383 * t309 + t359) * t308 + (-t72 * mrSges(6,2) - t396 * mrSges(5,3) + t46 * mrSges(7,3) + t187 / 0.2e1 + t188 / 0.2e1 + t189 / 0.2e1 - t464) * t310) * t311;
t409 = t4 * qJD(1);
t364 = mrSges(6,3) / 0.2e1 + mrSges(7,2) / 0.2e1 - mrSges(5,2) / 0.2e1;
t384 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t366 = mrSges(5,1) / 0.2e1 + t384;
t317 = (t366 * t308 - t364 * t310 - t487) * t309;
t322 = (t430 * t308 + t429 * t310) * t459 + (t427 * t308 + t428 * t310) * t461;
t373 = -t217 / 0.2e1 - t218 / 0.2e1;
t375 = t212 / 0.2e1 - t213 / 0.2e1;
t5 = (t443 + t373) * t310 + (-t211 / 0.2e1 + t375) * t308 + t317 + t322 + (-mrSges(6,2) / 0.2e1 + mrSges(7,3) / 0.2e1 - mrSges(5,3) / 0.2e1) * t465;
t408 = t5 * qJD(1);
t14 = (t216 + t218) * t309 + (t181 + t182) * t403 + m(7) * (t309 * t53 + t328 * t403) + m(6) * (t105 * t403 + t309 * t71);
t407 = qJD(1) * t14;
t400 = t310 * t216;
t405 = t308 * t211;
t22 = (t400 + t405 + t477) * t311;
t406 = qJD(1) * t22;
t397 = qJ(6) + t313;
t393 = qJD(2) * t310;
t392 = qJD(4) * t308;
t250 = t452 * t308;
t195 = t309 * t250;
t391 = t195 * qJD(1);
t372 = t307 / 0.2e1 + t306 / 0.2e1;
t200 = (0.1e1 / 0.2e1 + t372) * m(7);
t390 = t200 * qJD(2);
t389 = mrSges(7,1) - t476;
t388 = m(7) * t403;
t387 = t88 * t458;
t374 = -t214 / 0.2e1 - t219 / 0.2e1;
t298 = t310 * mrSges(7,2);
t361 = -t297 - t298 + t412;
t233 = -t308 * mrSges(7,1) + t298;
t360 = (t442 - t233 / 0.2e1) * t308;
t357 = t469 / 0.2e1 + t468 / 0.2e1 + t241 / 0.2e1;
t356 = t243 / 0.2e1 + t245 / 0.2e1 + t247 / 0.2e1;
t344 = t308 * t74 - t310 * t76;
t343 = t308 * t87 + t310 * t86;
t177 = qJ(3) - t489;
t235 = t294 + t432;
t10 = qJ(3) * t230 - t488 * t233 - t235 * t354 + (m(6) * t235 + t228) * t226 + (t353 + t435) * t177 + (t242 / 0.2e1 - t352 / 0.2e1 + t244 / 0.2e1 - t357) * t310 + (t349 / 0.2e1 - t348 / 0.2e1 - t347 / 0.2e1 - t356) * t308;
t220 = t397 * t308;
t221 = t397 * t310;
t299 = Ifges(6,6) * t310;
t315 = (t105 * t235 + t183 * t226) * t460 + (t177 * t490 + t430 * t220 + t429 * t221 + t328 * t488) * t458 + qJ(3) * t186 / 0.2e1 + t228 * t448 + t233 * t447 - t177 * t185 / 0.2e1 + t183 * t442 - t488 * t446 + t220 * t444 + t221 * t216 / 0.2e1 + t226 * t184 / 0.2e1 + t235 * t181 / 0.2e1 + t230 * t441 + t309 * t299 / 0.4e1 + t328 * t353 / 0.2e1;
t316 = (-pkin(4) * t76 + qJ(5) * t74) * t461 + (qJ(5) * t55 - t450 * t47) * t459 + pkin(4) * t445 + t450 * t481 + mrSges(7,1) * t455 - t55 * mrSges(7,2) / 0.2e1 - t74 * mrSges(6,3) / 0.2e1 + mrSges(6,1) * t454 - t86 * mrSges(5,1) / 0.2e1 + t87 * mrSges(5,2) / 0.2e1;
t318 = (-t62 / 0.2e1 + t53 / 0.2e1) * mrSges(7,3) + (-t71 / 0.2e1 + t453) * mrSges(6,2) + (t428 * t460 - t373) * t313 + t126 / 0.4e1 + t128 / 0.4e1 - t130 / 0.4e1 - t190 / 0.4e1 - t191 / 0.4e1 - t192 / 0.4e1;
t319 = (-t334 / 0.2e1 + t46 / 0.2e1) * mrSges(7,3) + (t396 / 0.2e1 - t72 / 0.2e1) * mrSges(6,2) + (t427 * t460 - t375) * t313 - t132 / 0.4e1 - t134 / 0.4e1 - t136 / 0.4e1 + t187 / 0.4e1 + t188 / 0.4e1 + t189 / 0.4e1;
t324 = t221 * t457 - t247 / 0.4e1 - t245 / 0.4e1 - t243 / 0.4e1 + t349 / 0.4e1 - t348 / 0.4e1 - t347 / 0.4e1;
t325 = t220 * t457 + t352 / 0.4e1 - t244 / 0.4e1 - t242 / 0.4e1 + t241 / 0.4e1 + t468 / 0.4e1 + t469 / 0.4e1;
t327 = t486 * t372 * t313;
t2 = t315 + ((-0.3e1 / 0.4e1 * Ifges(5,6) - 0.3e1 / 0.4e1 * Ifges(7,6) + Ifges(6,6) / 0.2e1) * t309 + t324 * t311 + t318) * t310 + ((-0.3e1 / 0.4e1 * Ifges(5,5) + 0.3e1 / 0.4e1 * Ifges(7,5) - 0.3e1 / 0.4e1 * Ifges(6,4)) * t309 + t325 * t311 + t319) * t308 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t327) * t311 + t374 * qJ(5) + t316;
t342 = t2 * qJD(1) + t10 * qJD(2);
t36 = mrSges(4,3) + (m(5) + m(4)) * qJ(3) + t436 + m(7) * t177 + t389 * t308 + t361;
t321 = -m(5) * t343 / 0.2e1 + t344 * t461 + (t55 * t308 - t47 * t310) * t459;
t333 = m(5) * t441 + m(6) * t448;
t8 = (t481 + t445 - t209 / 0.2e1 + t366 * t311) * t310 + (-t215 / 0.2e1 + t364 * t311 + t374) * t308 + ((-t220 * t310 + t221 * t308) * t309 + t328) * t458 + t321 + t333;
t341 = -qJD(1) * t8 - qJD(2) * t36;
t320 = (-t181 / 0.2e1 + t446) * t310 + (-t105 * t310 + (t226 * t311 + t309 * t313) * t308) * t460 + (t177 * t403 + t220 * t309 - t310 * t328) * t458;
t336 = m(6) * t454 + m(7) * t455;
t11 = t386 - t281 + (t360 + t384) * t311 + t320 - t336;
t38 = t177 * t434 + (-t233 - t354 + t436) * t310;
t340 = qJD(1) * t11 - qJD(2) * t38;
t326 = m(7) * ((t220 * t311 - t46) * t310 + (-t221 * t311 + t53) * t308);
t16 = (t414 / 0.2e1 + t444) * t310 + (t296 / 0.2e1 + t443) * t308 + t387 - t326 / 0.2e1;
t51 = m(7) * (t220 * t308 + t221 * t310) + t394 * mrSges(7,3);
t339 = -qJD(1) * t16 + qJD(2) * t51;
t37 = (-t381 / 0.4e1 + t280 / 0.4e1 + t490 / 0.4e1) * t463 - t185;
t54 = (t483 / 0.4e1 + t294 / 0.4e1 + t488 / 0.4e1) * t463 + t353;
t338 = qJD(1) * t37 + qJD(2) * t54;
t332 = 0.2e1 * t291 + t84;
t323 = -t295 - t296 + t332 * t461 + (t279 + t332) * t459;
t335 = m(6) * t453 + t62 * t458;
t20 = t323 + t335;
t251 = t452 * qJ(5) + mrSges(7,2) + mrSges(6,3);
t337 = qJD(1) * t20 - qJD(4) * t251;
t329 = mrSges(6,2) * pkin(4) - mrSges(7,3) * t450 + Ifges(7,5) - t475;
t201 = t394 * t458 + t459;
t194 = (qJD(1) * t403 - t393) * m(7);
t82 = m(7) * t220 + (m(6) * t313 - t431) * t308;
t77 = t488 * t458 - t435 / 0.2e1;
t45 = m(7) * t447 + t458 * t490;
t18 = t283 - t323 + t335 - t385;
t17 = t326 / 0.2e1 + t211 * t439 + t216 * t440 + t387 + (t415 / 0.2e1 + t413 / 0.2e1) * t309;
t13 = -t410 / 0.2e1 - t411 / 0.2e1 + t311 * t360 + t320 + t336;
t7 = t209 * t438 + (-t220 * t401 + t221 * t404 + t105) * t458 + (m(4) * pkin(7) + mrSges(4,1) + t364 * t308 + (pkin(5) * t458 + t366) * t310) * t311 - t321 + t333 + (t208 + t210) * t439 + (t214 + t472) * t440;
t6 = t405 / 0.2e1 + t400 / 0.2e1 + t213 * t440 - t308 * t212 / 0.2e1 + t317 - t322 + t395 * t438 + t457 * t465 + t486 * t394 * t437;
t1 = t315 + ((-Ifges(6,4) / 0.4e1 - Ifges(5,5) / 0.4e1 + Ifges(7,5) / 0.4e1) * t309 + t319) * t308 + ((-Ifges(5,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t309 + t318) * t310 - t316 + (t308 * t325 + t310 * t324 + t327) * t311 + (t214 + t219) * qJ(5) / 0.2e1 + t467 * t437 + (-Ifges(7,5) / 0.2e1 + t475 / 0.2e1) * t404 + (t456 + t485 / 0.2e1) * t401;
t12 = [qJD(2) * t3 - qJD(3) * t9 + qJD(4) * t4 + qJD(5) * t14 + qJD(6) * t22, t417 + t7 * qJD(3) + t1 * qJD(4) + t13 * qJD(5) + t17 * qJD(6) + (-t86 * mrSges(5,3) + t76 * mrSges(6,2) - t47 * mrSges(7,3) + (t209 - t210) * t313 - t330) * t393 + (-qJ(3) * t180 - t104 * t354 - t177 * t179 - t226 * t178 - t221 * t208 + t220 * t214 + t88 * t233 + (-t74 * mrSges(6,2) - t87 * mrSges(5,3) + t55 * mrSges(7,3) + t313 * t472 + t331) * t308 + 0.2e1 * (t104 * t226 + t313 * t344) * t460 + m(5) * (-qJ(3) * t248 + t313 * t343) + 0.2e1 * (-t177 * t88 + t220 * t55 - t221 * t47) * t458 + (-pkin(2) * mrSges(4,1) - Ifges(4,4) + Ifges(3,5) - t363 * t310 - t362 * t308 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2)) * pkin(7)) * t311 + (-qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6) + t357 * t310 + t356 * t308 + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(7)) * t309) * qJD(2), qJD(2) * t7 + qJD(4) * t6 - t418, t409 + t1 * qJD(2) + t6 * qJD(3) + t18 * qJD(5) + t45 * qJD(6) + ((-qJ(5) * t334 - t450 * t62) * t458 + (-pkin(4) * t84 - qJ(5) * t396) * t460) * t462 + (-t62 * mrSges(7,1) - t334 * mrSges(7,2) + t284 + ((-Ifges(6,6) - t466) * t308 + t329 * t310) * t311 + t476 * t84 - (-mrSges(5,2) + mrSges(6,3)) * t396) * qJD(4), qJD(2) * t13 + qJD(4) * t18 + t407, qJD(2) * t17 + qJD(4) * t45 + t406; qJD(3) * t8 + qJD(4) * t2 + qJD(5) * t11 - qJD(6) * t16 - t417, qJD(3) * t36 + qJD(4) * t10 - qJD(5) * t38 + qJD(6) * t51, qJD(6) * t201 - t341, t82 * qJD(5) + t77 * qJD(6) + t329 * t392 + t342 + (-t220 * mrSges(7,1) + t221 * mrSges(7,2) + t299 + m(7) * (qJ(5) * t221 - t220 * t450) + (-Ifges(5,6) + t466) * t310 + (m(6) * t369 + t354 + t355) * t313) * qJD(4), qJD(4) * t82 + t340, qJD(3) * t201 + qJD(4) * t77 + t339; -qJD(2) * t8 - qJD(4) * t5 + qJD(5) * t195 + t418, qJD(6) * t200 + t341, 0, -t361 * qJD(4) + t250 * qJD(5) - t389 * t392 + t462 * t487 - t408, qJD(4) * t250 + t391, t390; -qJD(2) * t2 + qJD(3) * t5 - qJD(5) * t20 + qJD(6) * t37 - t409, t54 * qJD(6) - t342, t408, t251 * qJD(5), -t337, t338; -qJD(2) * t11 - qJD(3) * t195 + qJD(4) * t20 + qJD(6) * t388 - t407, -qJD(6) * t434 - t340, -t391, t337, 0, t194; qJD(2) * t16 - qJD(4) * t37 - qJD(5) * t388 - t406, -qJD(3) * t200 - qJD(4) * t54 + qJD(5) * t434 - t339, -t390, -t338, -t194, 0;];
Cq  = t12;
