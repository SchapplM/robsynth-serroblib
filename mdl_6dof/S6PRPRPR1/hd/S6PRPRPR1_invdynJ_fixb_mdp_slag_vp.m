% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:28:06
% EndTime: 2019-03-08 19:28:12
% DurationCPUTime: 4.48s
% Computational Cost: add. (2863->394), mult. (6728->548), div. (0->0), fcn. (5621->16), ass. (0->193)
t510 = qJ(5) + pkin(8);
t412 = sin(qJ(2));
t415 = cos(qJ(2));
t404 = sin(pkin(6));
t471 = qJD(2) * t404;
t454 = qJD(1) * t471;
t464 = qJDD(1) * t404;
t509 = t412 * t464 + t415 * t454;
t408 = cos(pkin(6));
t402 = sin(pkin(11));
t406 = cos(pkin(11));
t437 = t402 * t415 + t406 * t412;
t349 = t437 * t408;
t480 = t415 * t406;
t360 = t402 * t412 - t480;
t403 = sin(pkin(10));
t407 = cos(pkin(10));
t310 = t349 * t407 - t360 * t403;
t311 = t403 * t349 + t360 * t407;
t508 = g(1) * t311 - g(2) * t310;
t348 = t437 * t404;
t411 = sin(qJ(4));
t414 = cos(qJ(4));
t327 = -t348 * t411 + t408 * t414;
t328 = t348 * t414 + t408 * t411;
t401 = sin(pkin(12));
t405 = cos(pkin(12));
t282 = t327 * t401 + t328 * t405;
t487 = t404 * t412;
t347 = t402 * t487 - t404 * t480;
t413 = cos(qJ(6));
t340 = t347 * t413;
t410 = sin(qJ(6));
t507 = -t282 * t410 + t340;
t472 = qJD(1) * t404;
t456 = t412 * t472;
t370 = t406 * t456;
t455 = t415 * t472;
t342 = t402 * t455 + t370;
t469 = qJD(4) * t411;
t506 = pkin(4) * t469 - t342;
t380 = t408 * qJDD(1) + qJDD(3);
t368 = t414 * t380;
t377 = t415 * t464;
t346 = qJDD(2) * pkin(2) - t412 * t454 + t377;
t304 = t402 * t346 + t406 * t509;
t300 = qJDD(2) * pkin(8) + t304;
t382 = qJD(1) * t408 + qJD(3);
t422 = qJ(5) * qJDD(2) + qJD(2) * qJD(5) + qJD(4) * t382 + t300;
t365 = qJD(2) * pkin(2) + t455;
t332 = t402 * t365 + t370;
t446 = qJD(2) * t510 + t332;
t435 = t446 * qJD(4);
t267 = qJDD(4) * pkin(4) - t422 * t411 - t414 * t435 + t368;
t268 = (t380 - t435) * t411 + t422 * t414;
t259 = t267 * t405 - t268 * t401;
t257 = -qJDD(4) * pkin(5) - t259;
t470 = qJD(2) * t411;
t484 = t405 * t414;
t352 = qJD(2) * t484 - t401 * t470;
t351 = qJD(6) - t352;
t361 = t401 * t414 + t405 * t411;
t354 = t361 * qJD(2);
t386 = pkin(4) * t401 + pkin(9);
t398 = qJ(4) + pkin(12);
t393 = sin(t398);
t394 = cos(t398);
t488 = t404 * t407;
t490 = t403 * t404;
t505 = (pkin(4) * t470 + pkin(5) * t354 - pkin(9) * t352 + qJD(6) * t386) * t351 + g(1) * (t311 * t393 + t394 * t490) + g(2) * (-t310 * t393 - t394 * t488) + g(3) * (-t348 * t393 + t394 * t408) + t257;
t260 = t401 * t267 + t405 * t268;
t258 = qJDD(4) * pkin(9) + t260;
t301 = t414 * t382 - t446 * t411;
t296 = qJD(4) * pkin(4) + t301;
t302 = t382 * t411 + t446 * t414;
t498 = t302 * t401;
t273 = t296 * t405 - t498;
t271 = -qJD(4) * pkin(5) - t273;
t369 = t402 * t456;
t331 = t365 * t406 - t369;
t391 = pkin(4) * t414 + pkin(3);
t320 = -t391 * qJD(2) + qJD(5) - t331;
t278 = -pkin(5) * t352 - pkin(9) * t354 + t320;
t359 = t401 * t411 - t484;
t502 = pkin(2) * t406;
t436 = -t391 - t502;
t314 = pkin(5) * t359 - pkin(9) * t361 + t436;
t387 = pkin(2) * t402 + pkin(8);
t479 = qJ(5) + t387;
t358 = t479 * t414;
t450 = t479 * t411;
t318 = t405 * t358 - t401 * t450;
t353 = t361 * qJD(4);
t462 = qJDD(2) * t414;
t463 = qJDD(2) * t411;
t441 = -t401 * t463 + t405 * t462;
t321 = qJD(2) * t353 + qJDD(6) - t441;
t356 = t359 * qJD(4);
t428 = -g(3) * t348 + t508;
t448 = qJD(4) * t479;
t341 = qJD(5) * t414 - t411 * t448;
t345 = t406 * t455 - t369;
t425 = -qJD(5) * t411 - t414 * t448;
t475 = t405 * t341 + t359 * t345 + t401 * t425;
t504 = -(qJD(6) * t278 + t258) * t359 + t257 * t361 - t271 * t356 + (-qJD(6) * t314 - t475) * t351 - t318 * t321 + t428;
t429 = t360 * t408;
t309 = -t403 * t437 - t407 * t429;
t312 = t403 * t429 - t407 * t437;
t427 = g(1) * t312 + g(2) * t309 - g(3) * t347;
t503 = t314 * t321 - t427 * t394;
t501 = g(3) * t327;
t465 = qJD(2) * qJD(4);
t452 = t414 * t465;
t453 = t411 * t465;
t323 = t361 * qJDD(2) - t401 * t453 + t405 * t452;
t466 = t413 * qJD(4);
t457 = qJD(6) * t466 + t410 * qJDD(4) + t413 * t323;
t467 = qJD(6) * t410;
t279 = -t354 * t467 + t457;
t500 = t279 * t410;
t491 = t354 * t410;
t334 = -t466 + t491;
t496 = t334 * t351;
t495 = t334 * t354;
t336 = qJD(4) * t410 + t354 * t413;
t494 = t336 * t351;
t493 = t336 * t354;
t492 = t347 * t410;
t489 = t403 * t412;
t486 = t404 * t414;
t485 = t404 * t415;
t294 = t405 * t302;
t483 = t408 * t412;
t482 = t408 * t415;
t481 = t410 * t321;
t315 = t413 * t321;
t478 = qJDD(1) - g(3);
t477 = t279 * t359 + t336 * t353;
t274 = t401 * t296 + t294;
t476 = t341 * t401 - t361 * t345 - t405 * t425;
t474 = pkin(5) * t353 + pkin(9) * t356 + t506;
t399 = t411 ^ 2;
t473 = -t414 ^ 2 + t399;
t468 = qJD(6) * t361;
t460 = t361 * t481;
t459 = t361 * t315;
t458 = t407 * t482;
t449 = -t413 * qJDD(4) + t323 * t410;
t447 = t351 * t413;
t442 = g(1) * t403 - g(2) * t407;
t303 = t346 * t406 - t402 * t509;
t272 = qJD(4) * pkin(9) + t274;
t262 = t272 * t413 + t278 * t410;
t440 = t272 * t410 - t278 * t413;
t280 = t336 * qJD(6) + t449;
t439 = -t280 * t359 - t334 * t353;
t438 = t282 * t413 + t492;
t434 = t315 + (t352 * t410 - t467) * t351;
t433 = -t403 * t482 - t407 * t412;
t432 = t356 * t410 - t413 * t468;
t431 = t356 * t413 + t361 * t467;
t430 = -g(1) * t490 + g(2) * t488 - g(3) * t408;
t329 = -qJD(2) * pkin(3) - t331;
t424 = -qJD(2) * t329 - t300 - t508;
t276 = t301 * t405 - t498;
t423 = -t386 * t321 + (t271 + t276) * t351;
t389 = -pkin(3) - t502;
t421 = -qJDD(4) * t387 + (qJD(2) * t389 + t329 + t345) * qJD(4);
t420 = -g(1) * t433 - g(3) * t485;
t283 = pkin(4) * t453 - t391 * qJDD(2) + qJDD(5) - t303;
t416 = qJD(4) ^ 2;
t419 = -qJD(2) * t342 + t387 * t416 - t303 + t427 + (-pkin(3) + t389) * qJDD(2);
t417 = qJD(2) ^ 2;
t388 = -pkin(4) * t405 - pkin(5);
t374 = pkin(2) * t458;
t373 = qJDD(4) * t414 - t411 * t416;
t372 = qJDD(4) * t411 + t414 * t416;
t344 = t360 * t471;
t343 = qJD(2) * t348;
t325 = t348 * t394 + t393 * t408;
t322 = -qJD(4) * t354 + t441;
t317 = t358 * t401 + t405 * t450;
t290 = t327 * qJD(4) - t344 * t414;
t289 = -t328 * qJD(4) + t344 * t411;
t288 = -t311 * t394 + t393 * t490;
t286 = t310 * t394 - t393 * t488;
t281 = -t405 * t327 + t328 * t401;
t275 = t301 * t401 + t294;
t270 = t289 * t401 + t290 * t405;
t269 = -t405 * t289 + t290 * t401;
t265 = -pkin(5) * t322 - pkin(9) * t323 + t283;
t263 = t413 * t265;
t1 = [t478 * MDP(1) + (-t303 * t347 + t304 * t348 - t331 * t343 - t332 * t344 + t380 * t408 - g(3)) * MDP(5) + (qJD(4) * t289 + qJDD(4) * t327 - t347 * t462) * MDP(11) + (-qJD(4) * t290 - qJDD(4) * t328 + t347 * t463) * MDP(12) + (t269 * t354 + t270 * t352 + t281 * t323 + t282 * t322) * MDP(13) + (-t259 * t281 + t260 * t282 - t269 * t273 + t270 * t274 + t283 * t347 + t320 * t343 - g(3)) * MDP(14) + ((-t438 * qJD(6) - t270 * t410 + t343 * t413) * t351 + t507 * t321 + t269 * t334 + t281 * t280) * MDP(20) + (-(qJD(6) * t507 + t270 * t413 + t343 * t410) * t351 - t438 * t321 + t269 * t336 + t281 * t279) * MDP(21) + ((-t343 * t414 + t347 * t469) * MDP(11) + (qJD(4) * t347 * t414 + t343 * t411) * MDP(12)) * qJD(2) + ((qJDD(2) * t415 - t412 * t417) * MDP(3) + (-qJDD(2) * t412 - t415 * t417) * MDP(4)) * t404; qJDD(2) * MDP(2) + (t377 - g(2) * (t458 - t489) + t420) * MDP(3) + (-g(1) * (t403 * t483 - t407 * t415) - g(2) * (-t403 * t415 - t407 * t483) - t478 * t487) * MDP(4) + (-g(2) * t374 + t331 * t342 - t332 * t345 + (g(2) * t489 + t303 * t406 + t304 * t402 + t420) * pkin(2)) * MDP(5) + (qJDD(2) * t399 + 0.2e1 * t411 * t452) * MDP(6) + 0.2e1 * (t411 * t462 - t473 * t465) * MDP(7) + t372 * MDP(8) + t373 * MDP(9) + (t421 * t411 - t419 * t414) * MDP(11) + (t419 * t411 + t421 * t414) * MDP(12) + (-t259 * t361 - t260 * t359 + t273 * t356 - t274 * t353 + t317 * t323 + t318 * t322 + t475 * t352 + t476 * t354 + t428) * MDP(13) + (t260 * t318 - t259 * t317 + t283 * t436 - g(1) * (t433 * pkin(2) - t311 * t510 + t312 * t391) - g(2) * (-pkin(2) * t489 + t309 * t391 + t310 * t510 + t374) - g(3) * (pkin(2) * t485 - t347 * t391 + t348 * t510) + t506 * t320 + t475 * t274 - t476 * t273) * MDP(14) + (t279 * t361 * t413 - t431 * t336) * MDP(15) + (-(-t334 * t413 - t336 * t410) * t356 + (-t500 - t280 * t413 + (t334 * t410 - t336 * t413) * qJD(6)) * t361) * MDP(16) + (-t431 * t351 + t459 + t477) * MDP(17) + (t432 * t351 + t439 - t460) * MDP(18) + (t321 * t359 + t351 * t353) * MDP(19) + (-t440 * t353 + t263 * t359 + t317 * t280 + t476 * t334 + (t474 * t351 + (t271 * t361 - t272 * t359 - t318 * t351) * qJD(6) + t503) * t413 + t504 * t410) * MDP(20) + (-t262 * t353 + t317 * t279 + t476 * t336 + (-(-qJD(6) * t272 + t265) * t359 - t271 * t468 + (qJD(6) * t318 - t474) * t351 - t503) * t410 + t504 * t413) * MDP(21); (t430 + t380) * MDP(5) + t373 * MDP(11) - t372 * MDP(12) + (t322 * t361 + t323 * t359 - t352 * t356 + t353 * t354) * MDP(13) + (-t259 * t359 + t260 * t361 - t273 * t353 - t274 * t356 + t430) * MDP(14) + (-t439 - t460) * MDP(20) + (-t459 + t477) * MDP(21) + (t432 * MDP(20) + t431 * MDP(21)) * t351; MDP(8) * t463 + MDP(9) * t462 + qJDD(4) * MDP(10) + (t424 * t411 - t442 * t486 + t368 - t501) * MDP(11) + (g(3) * t328 + (t442 * t404 - t380) * t411 + t424 * t414) * MDP(12) + ((t274 - t275) * t354 + (t273 - t276) * t352 + (t322 * t401 - t323 * t405) * pkin(4)) * MDP(13) + (t273 * t275 - t274 * t276 + (t260 * t401 + t259 * t405 - t320 * t470 - g(1) * (t311 * t411 + t403 * t486) - g(2) * (-t310 * t411 - t407 * t486) - t501) * pkin(4)) * MDP(14) + (t336 * t447 + t500) * MDP(15) + ((t279 - t496) * t413 + (-t280 - t494) * t410) * MDP(16) + (t351 * t447 + t481 - t493) * MDP(17) + (t434 + t495) * MDP(18) - t351 * t354 * MDP(19) + (-t275 * t334 + t388 * t280 + t354 * t440 + t423 * t410 - t413 * t505) * MDP(20) + (t262 * t354 - t275 * t336 + t388 * t279 + t410 * t505 + t423 * t413) * MDP(21) + (-t411 * t414 * MDP(6) + t473 * MDP(7)) * t417; (-t352 ^ 2 - t354 ^ 2) * MDP(13) + (t273 * t354 - t274 * t352 + t283 + t427) * MDP(14) + (t434 - t495) * MDP(20) + (-t351 ^ 2 * t413 - t481 - t493) * MDP(21); t336 * t334 * MDP(15) + (-t334 ^ 2 + t336 ^ 2) * MDP(16) + (t457 + t496) * MDP(17) + (-t449 + t494) * MDP(18) + t321 * MDP(19) + (-t410 * t258 + t263 + t262 * t351 - t271 * t336 - g(1) * (-t288 * t410 - t312 * t413) - g(2) * (-t286 * t410 - t309 * t413) - g(3) * (-t325 * t410 + t340)) * MDP(20) + (-t413 * t258 - t410 * t265 - t440 * t351 + t271 * t334 - g(1) * (-t288 * t413 + t312 * t410) - g(2) * (-t286 * t413 + t309 * t410) - g(3) * (-t325 * t413 - t492)) * MDP(21) + (-MDP(17) * t491 - t336 * MDP(18) - t262 * MDP(20) + t440 * MDP(21)) * qJD(6);];
tau  = t1;
