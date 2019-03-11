% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:19:58
% EndTime: 2019-03-09 05:20:10
% DurationCPUTime: 7.14s
% Computational Cost: add. (22909->487), mult. (40314->625), div. (0->0), fcn. (46762->8), ass. (0->286)
t284 = sin(qJ(6));
t281 = t284 ^ 2;
t287 = cos(qJ(6));
t282 = t287 ^ 2;
t396 = t281 + t282;
t541 = -t284 / 0.2e1;
t496 = t287 / 0.2e1;
t285 = sin(qJ(4));
t286 = sin(qJ(3));
t288 = cos(qJ(3));
t494 = cos(qJ(4));
t262 = -t285 * t288 - t286 * t494;
t271 = t286 * pkin(3) + qJ(2);
t233 = -pkin(4) * t262 + t271;
t540 = m(6) * t233;
t283 = sin(pkin(10));
t289 = -pkin(1) - pkin(7);
t491 = pkin(8) - t289;
t346 = t491 * t494;
t360 = t285 * t491;
t347 = t288 * t360;
t226 = -t286 * t346 - t347;
t250 = t262 * qJ(5);
t306 = t250 + t226;
t451 = cos(pkin(10));
t260 = t288 * t346;
t224 = t286 * t360 - t260;
t261 = t285 * t286 - t288 * t494;
t249 = t261 * qJ(5);
t522 = t249 + t224;
t122 = -t283 * t306 + t451 * t522;
t215 = t283 * t261 + t451 * t262;
t443 = t122 * t215;
t460 = t287 * mrSges(7,1);
t465 = t284 * mrSges(7,2);
t344 = t460 - t465;
t523 = -t451 * t261 + t262 * t283;
t532 = t344 * t523;
t539 = -t532 / 0.2e1;
t459 = t287 * mrSges(7,2);
t466 = t284 * mrSges(7,1);
t265 = t459 + t466;
t533 = t265 * t523;
t538 = -t533 / 0.2e1;
t537 = t215 * t344;
t355 = t396 * t215;
t529 = t283 * t522 + t306 * t451;
t447 = t529 * t523;
t467 = t282 * mrSges(7,3);
t468 = t281 * mrSges(7,3);
t536 = -t467 / 0.2e1 - t468 / 0.2e1;
t391 = t494 * pkin(3);
t275 = t391 + pkin(4);
t358 = t451 * t285;
t241 = pkin(3) * t358 + t283 * t275;
t422 = t241 * t523;
t415 = t283 * t285;
t240 = -pkin(3) * t415 + t275 * t451;
t423 = t240 * t215;
t399 = t423 + t422;
t535 = t396 * t523;
t484 = Ifges(7,4) * t284;
t343 = Ifges(7,1) * t287 - t484;
t278 = Ifges(7,4) * t287;
t481 = Ifges(7,2) * t284;
t521 = t278 - t481;
t266 = Ifges(7,2) * t287 + t484;
t267 = Ifges(7,1) * t284 + t278;
t531 = t266 * t541 + t267 * t496;
t534 = t284 / 0.2e1;
t349 = t343 * t534 + t521 * t496 + t531;
t483 = Ifges(7,5) * t215;
t480 = Ifges(7,6) * t215;
t473 = t215 * mrSges(6,3);
t474 = t523 * mrSges(6,3);
t280 = t288 * pkin(3);
t492 = t261 * pkin(4);
t236 = t280 - t492;
t464 = t284 * mrSges(7,3);
t149 = -mrSges(7,2) * t523 - t215 * t464;
t488 = mrSges(7,3) * t215;
t151 = mrSges(7,1) * t523 - t287 * t488;
t400 = t149 * t534 + t151 * t496;
t506 = m(7) / 0.2e1;
t508 = m(6) / 0.2e1;
t359 = t491 * t286;
t227 = t285 * t359 - t260;
t193 = t249 + t227;
t225 = t359 * t494 + t347;
t307 = -t250 + t225;
t124 = t193 * t451 + t283 * t307;
t130 = pkin(5) * t523 - pkin(9) * t215 - t492;
t129 = t130 + t280;
t59 = -t124 * t284 + t129 * t287;
t60 = t124 * t287 + t129 * t284;
t530 = t236 * t508 + (t284 * t60 + t287 * t59) * t506 + t400;
t528 = t282 / 0.2e1;
t527 = -mrSges(6,1) - t344;
t330 = (t528 + t281 / 0.2e1) * mrSges(7,3) * t523;
t524 = -pkin(3) * t261 * t285 + t262 * t391;
t429 = t523 * t284;
t150 = mrSges(7,2) * t215 - mrSges(7,3) * t429;
t428 = t523 * t287;
t152 = -mrSges(7,1) * t215 - mrSges(7,3) * t428;
t493 = pkin(4) * t283;
t269 = pkin(9) + t493;
t303 = t483 + (0.5e1 / 0.4e1 * t484 + t266 / 0.4e1 + (-Ifges(7,1) / 0.2e1 + Ifges(7,2) / 0.4e1) * t287) * t523;
t304 = -t480 + (t267 / 0.4e1 + t278 / 0.4e1 + (Ifges(7,1) / 0.4e1 - Ifges(7,2) / 0.2e1) * t284) * t523;
t501 = t523 / 0.2e1;
t382 = Ifges(7,3) * t501;
t441 = t122 * t265;
t326 = t382 + t441 / 0.2e1;
t375 = t451 * pkin(4);
t270 = -t375 - pkin(5);
t419 = t270 * t532;
t497 = t269 / 0.2e1;
t503 = -mrSges(7,2) / 0.2e1;
t504 = mrSges(7,1) / 0.2e1;
t63 = -t122 * t284 + t130 * t287;
t64 = t122 * t287 + t130 * t284;
t517 = t64 * t503 + t63 * t504;
t12 = -t419 / 0.2e1 + t269 * t330 + (t150 * t497 + t304) * t284 + (t152 * t497 + t303) * t287 + t326 + t517;
t417 = t270 * t265;
t156 = t349 + t417;
t327 = -t459 / 0.2e1 - t466 / 0.2e1;
t320 = t327 * t523;
t84 = t533 / 0.2e1 + t320;
t401 = t84 * qJD(2);
t237 = -pkin(5) - t240;
t248 = (t451 * t494 - t415) * pkin(3);
t319 = t327 * t248;
t74 = (-t237 / 0.2e1 - t270 / 0.2e1) * t265 + t319 - t349;
t520 = t12 * qJD(1) + t74 * qJD(3) - t156 * qJD(4) + t401;
t238 = pkin(9) + t241;
t426 = t237 * t532;
t500 = t238 / 0.2e1;
t516 = t60 * t503 + t59 * t504;
t10 = -t426 / 0.2e1 + t238 * t330 + (t150 * t500 + t304) * t284 + (t152 * t500 + t303) * t287 + t326 + t516;
t424 = t237 * t265;
t137 = t349 + t424;
t519 = -t10 * qJD(1) + t137 * qJD(3) - t401;
t408 = t287 * t150;
t411 = t284 * t152;
t325 = t408 / 0.2e1 - t411 / 0.2e1;
t409 = t287 * t149;
t412 = t284 * t151;
t518 = t409 / 0.2e1 - t412 / 0.2e1;
t515 = -t286 * mrSges(4,1) - t288 * mrSges(4,2);
t335 = -t284 * t63 + t287 * t64;
t514 = mrSges(7,3) * t535 + t537;
t413 = t284 * t150;
t369 = -t413 / 0.2e1;
t513 = t523 * t536 + t369;
t511 = (-t285 * mrSges(5,1) - t494 * mrSges(5,2)) * pkin(3);
t128 = -pkin(5) * t215 - pkin(9) * t523 + t233;
t57 = t128 * t287 - t284 * t529;
t58 = t128 * t284 + t287 * t529;
t337 = -t284 * t57 + t287 * t58;
t147 = t265 * t215;
t425 = t237 * t147;
t247 = (t283 * t494 + t358) * pkin(3);
t442 = t122 * t247;
t470 = t226 * mrSges(5,1);
t472 = t224 * mrSges(5,2);
t476 = t122 * mrSges(6,2);
t498 = -t248 / 0.2e1;
t499 = t247 / 0.2e1;
t510 = (t122 * t241 - t442 + (-t240 + t248) * t529) * t508 + (t237 * t529 + t238 * t335 + t248 * t337 - t442) * t506 - t476 / 0.2e1 - t472 / 0.2e1 - t470 / 0.2e1 + t425 / 0.2e1 + t533 * t499 + (-t215 * t498 + t523 * t499 - t422 / 0.2e1 - t423 / 0.2e1) * mrSges(6,3);
t509 = m(5) / 0.2e1;
t507 = -m(7) / 0.2e1;
t505 = m(6) * pkin(4);
t295 = (-t473 / 0.2e1 - t147 / 0.2e1 + t325) * t523 + (t474 / 0.2e1 - t518 + t538) * t215;
t312 = t337 * t523 + t443;
t454 = t60 * t287;
t456 = t59 * t284;
t336 = t454 - t456;
t397 = t225 + t226;
t398 = t224 - t227;
t121 = t193 * t283 - t307 * t451;
t440 = t121 * t523;
t6 = (-t215 * t336 + t312 - t440) * t506 + (-t124 * t215 - t440 + t443 + t447) * t508 + (-t261 * t397 + t262 * t398) * t509 + t295;
t8 = (-t215 * t335 + t312 - t447) * t506 + t295;
t502 = t6 * qJD(3) + t8 * qJD(4);
t37 = m(7) * (t215 - t355) * t523;
t403 = t37 * qJD(2);
t85 = t538 + t320;
t490 = t85 * qJD(6) + t403;
t489 = -t84 * qJD(6) - t403;
t485 = Ifges(5,4) * t261;
t482 = Ifges(7,5) * t287;
t479 = Ifges(7,6) * t284;
t478 = t529 * mrSges(6,1);
t477 = t121 * mrSges(6,1);
t475 = t124 * mrSges(6,2);
t471 = t225 * mrSges(5,1);
t469 = t227 * mrSges(5,2);
t455 = t6 * qJD(1);
t453 = t8 * qJD(1);
t9 = -t122 * t532 + t57 * t150 - t58 * t152 + ((-mrSges(7,3) * t58 - Ifges(7,4) * t428 + t480) * t287 + (t57 * mrSges(7,3) + t483 + (t484 + (-Ifges(7,1) + Ifges(7,2)) * t287) * t523) * t284) * t523;
t452 = t9 * qJD(1);
t444 = t122 * t523;
t14 = (t533 + t474) * t523 + (t408 - t411 + t473) * t215 + m(7) * (t215 * t337 - t444) + m(6) * (t215 * t529 - t444);
t450 = qJD(1) * t14;
t158 = -mrSges(6,1) * t215 + mrSges(6,2) * t523;
t357 = -t262 * mrSges(5,1) - t261 * mrSges(5,2);
t322 = -t158 - t357;
t308 = -t322 - t515;
t406 = t287 * t152;
t26 = t413 + t406 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(7) * (t284 * t58 + t287 * t57) + t540 + m(5) * t271 + t308;
t449 = qJD(1) * t26;
t446 = t529 * t344;
t445 = t122 * t121;
t439 = t121 * t344;
t157 = mrSges(6,1) * t523 + t215 * mrSges(6,2);
t376 = t539 + t396 * t488 / 0.2e1;
t297 = (t241 * t215 - t240 * t523) * t508 + (t237 * t523 + t238 * t355) * t506 + t376;
t18 = t297 - t157 - t530;
t438 = t18 * qJD(1);
t353 = t396 * t269;
t298 = (t215 * t353 + t270 * t523) * t506 + (t215 * t283 - t451 * t523) * t505 / 0.2e1;
t394 = -t505 / 0.2e1;
t305 = (t284 * t64 + t287 * t63) * t506 + t261 * t394 + t400;
t21 = t532 / 0.2e1 - t298 + t305 + t157 + t536 * t215;
t437 = t21 * qJD(1);
t436 = t523 ^ 2;
t433 = t247 * t523;
t364 = -t406 / 0.2e1;
t300 = -(t364 + t369 - t330) * t215 + t523 * t539;
t328 = -t465 / 0.2e1 + t460 / 0.2e1;
t22 = t300 - t328;
t427 = t22 * qJD(1);
t321 = t327 * t215;
t27 = t321 - t325;
t420 = t27 * qJD(1);
t418 = t270 * t147;
t404 = t287 * t269;
t301 = (-t215 * t355 - t436) * t506 + (-t215 ^ 2 - t436) * t508;
t329 = -m(6) / 0.2e1 + t396 * t507;
t43 = t301 + t329;
t402 = t43 * qJD(1);
t395 = qJD(3) + qJD(4);
t392 = mrSges(7,3) * t454;
t386 = t64 * mrSges(7,3) / 0.2e1;
t384 = t269 * t412;
t379 = -t464 / 0.2e1;
t378 = t464 / 0.2e1;
t377 = -t237 * t215 + t238 * t535;
t363 = -t404 / 0.2e1;
t354 = t396 * t248;
t352 = t474 * t493;
t345 = Ifges(5,4) * t262 + (-Ifges(5,1) + Ifges(5,2)) * t261;
t341 = -t479 + t482;
t340 = Ifges(7,5) * t284 + Ifges(7,6) * t287;
t100 = Ifges(7,5) * t523 + t215 * t343;
t99 = Ifges(7,6) * t523 + t215 * t521;
t294 = (Ifges(6,1) + Ifges(7,1) * t528 + (-t278 + t481 / 0.2e1) * t284) * t215 - Ifges(6,4) * t523 + t341 * t501 + t99 * t541 + t100 * t496;
t302 = -t529 * t474 + t58 * t149 + t57 * t151 + t233 * t157 + t271 * (-mrSges(5,1) * t261 + mrSges(5,2) * t262) - (t147 + t473) * t122;
t309 = (Ifges(6,2) + Ifges(7,3)) * t523 + (-Ifges(6,4) + t341) * t215;
t1 = t302 + (t121 * mrSges(6,3) + t294) * t523 + (-mrSges(5,3) * t398 + t345) * t262 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t286 + (-Ifges(4,1) + Ifges(4,2)) * t288) * t286 - (-t124 * mrSges(6,3) + t309) * t215 + (mrSges(5,3) * t397 - t485) * t261 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t288 + pkin(3) * t357) * t288 + t236 * t158 + t59 * t152 + t121 * t533 + t60 * t150 + m(5) * (t224 * t225 + t226 * t227 + t271 * t280) + m(6) * (t124 * t529 + t233 * t236 - t445) + m(7) * (t57 * t59 + t58 * t60 - t445);
t339 = t1 * qJD(1) + t6 * qJD(2);
t3 = t529 * t533 + t64 * t150 + t63 * t152 + (-pkin(4) * t158 - t485) * t261 - t492 * t540 + m(7) * (-t122 * t529 + t57 * t63 + t58 * t64) + t345 * t262 + (mrSges(6,3) * t529 + t294) * t523 - (-t122 * mrSges(6,3) + t309) * t215 + t302;
t338 = t3 * qJD(1) + t8 * qJD(2);
t333 = t375 * t473;
t293 = (-t215 * t354 + t377 - t433) * t506 + (-t215 * t248 + t399 - t433) * t508;
t313 = m(7) * (-t215 * t270 + t353 * t523);
t315 = (t215 * t451 + t283 * t523) * t505;
t299 = -t313 / 0.2e1 - t315 / 0.2e1;
t25 = t293 + t299;
t291 = (t121 * t270 + t269 * t336) * t507 + t477 / 0.2e1 + t439 / 0.2e1 + t475 / 0.2e1 - t471 / 0.2e1 + t469 / 0.2e1 - t418 / 0.2e1 + (-t121 * t451 + t124 * t283) * t394 + t384 / 0.2e1 + t149 * t363 + t59 * t378 - t392 / 0.2e1 + t352 / 0.2e1 + t333 / 0.2e1;
t4 = -t446 / 0.2e1 - t478 / 0.2e1 + t291 + t287 * t386 + t63 * t379 + t325 * t248 + t518 * t238 + t510;
t45 = t527 * t247 + t511 + (mrSges(7,3) * t396 - mrSges(6,2)) * t248 + m(6) * (-t240 * t247 + t241 * t248) + m(7) * (t237 * t247 + t238 * t354);
t318 = t4 * qJD(1) + t25 * qJD(2) + t45 * qJD(3);
t314 = t340 * t501 + Ifges(5,6) * t261 + Ifges(5,5) * t262 + t100 * t534 + t99 * t496 - Ifges(6,6) * t523 + (Ifges(6,5) + t531) * t215;
t311 = t322 + t514;
t296 = -t441 / 0.2e1 - t266 * t428 / 0.2e1 - t284 * (t521 * t523 - t480) / 0.4e1 + t382 + (0.2e1 * t343 * t523 - t483) * t287 / 0.4e1 + (t379 + t378) * t58 - (0.2e1 * t267 + t521) * t429 / 0.4e1 + (-t341 / 0.4e1 + t482 / 0.2e1 - t479 / 0.2e1) * t215;
t75 = t424 / 0.2e1 + t417 / 0.2e1 + t319 + t349;
t42 = t301 - t329;
t28 = t321 + t325;
t24 = t298 + t305 + t376;
t23 = t300 + t328;
t20 = t297 + t530;
t19 = t293 - t299 + t311;
t13 = t296 + t419 / 0.2e1 + t152 * t363 + t513 * t269 + t517;
t11 = t296 + t426 / 0.2e1 + (t364 + t513) * t238 + t516;
t2 = (t152 * t498 - t238 * t151 / 0.2e1 - t63 * mrSges(7,3) / 0.2e1) * t284 + (-t344 / 0.2e1 - mrSges(6,1) / 0.2e1) * t529 + (t248 * t150 / 0.2e1 + t149 * t500 + t386) * t287 - t291 + t314 + t510;
t5 = [qJD(2) * t26 + qJD(3) * t1 + qJD(4) * t3 + qJD(5) * t14 + qJD(6) * t9, qJD(5) * t42 + qJD(6) * t23 + t449 + t502 (m(5) * (t225 * t494 + t227 * t285) * pkin(3) - Ifges(4,6) * t288 - Ifges(4,5) * t286 - t439 + t425 + t471 - t469 - t477 - t475 + m(7) * t121 * t237 + t314 + m(6) * (-t121 * t240 + t124 * t241) + t392 - mrSges(7,3) * t456 + t515 * t289 + (m(7) * t336 + t409 - t412) * t238 - t399 * mrSges(6,3) - t524 * mrSges(5,3)) * qJD(3) + t2 * qJD(4) + t20 * qJD(5) + t11 * qJD(6) + t339, t2 * qJD(3) + (-t333 + (t122 * t283 - t451 * t529) * t505 + t418 - t446 - t472 - t470 - t476 - t478 + m(7) * (t269 * t335 + t270 * t529) + t314 + t149 * t404 - t384 - t352 + t335 * mrSges(7,3)) * qJD(4) + t24 * qJD(5) + t13 * qJD(6) + t338, qJD(2) * t42 + qJD(3) * t20 + qJD(4) * t24 + qJD(6) * t28 + t450, t452 + t23 * qJD(2) + t11 * qJD(3) + t13 * qJD(4) + t28 * qJD(5) + (-mrSges(7,1) * t58 - mrSges(7,2) * t57 - t340 * t523) * qJD(6); qJD(5) * t43 + qJD(6) * t22 - t449 + t502, t395 * t37, t455 + (-t308 + t514) * qJD(3) + t19 * qJD(4) + 0.2e1 * (t377 * t506 + t399 * t508 + t509 * t524) * qJD(3) + t490, t453 + t19 * qJD(3) + (t313 + t315 + t311) * qJD(4) + t490, t402, qJD(6) * t537 + t395 * t85 + t427; qJD(4) * t4 + qJD(5) * t18 - qJD(6) * t10 - t339, qJD(4) * t25 - t455 + t489, qJD(4) * t45 + qJD(6) * t137, t75 * qJD(6) + t318 + (t511 + (m(7) * t270 - t451 * t505 + t527) * t247 + (m(7) * t353 + t283 * t505 - mrSges(6,2) + t467 + t468) * t248) * qJD(4), t438, t75 * qJD(4) + (-t238 * t344 + t341) * qJD(6) + t519; -qJD(3) * t4 - qJD(5) * t21 - qJD(6) * t12 - t338, -qJD(3) * t25 - t453 + t489, -qJD(6) * t74 - t318, t156 * qJD(6), -t437 (-t269 * t344 + t341) * qJD(6) - t520; -qJD(2) * t43 - qJD(3) * t18 + qJD(4) * t21 - qJD(6) * t27 - t450, -t402, -t438, t437, 0, -qJD(6) * t265 - t420; -qJD(2) * t22 + qJD(3) * t10 + qJD(4) * t12 + qJD(5) * t27 - t452, t395 * t84 - t427, qJD(4) * t74 - t519, t520, t420, 0;];
Cq  = t5;
