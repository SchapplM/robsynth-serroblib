% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPPR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:00:08
% EndTime: 2019-03-09 03:00:14
% DurationCPUTime: 4.29s
% Computational Cost: add. (1571->428), mult. (2710->497), div. (0->0), fcn. (1383->6), ass. (0->198)
t355 = cos(qJ(3));
t446 = qJD(1) * t355;
t431 = qJD(1) * qJD(3);
t411 = t355 * t431;
t352 = sin(qJ(3));
t428 = qJDD(1) * t352;
t493 = t411 + t428;
t357 = -pkin(3) - pkin(4);
t358 = -pkin(1) - pkin(7);
t305 = t358 * qJD(1) + qJD(2);
t276 = (qJ(5) * qJD(1) + t305) * t355;
t492 = qJD(4) - t276;
t266 = qJD(3) * t357 + t492;
t458 = qJ(5) + t358;
t356 = cos(qJ(1));
t341 = g(2) * t356;
t353 = sin(qJ(1));
t342 = g(1) * t353;
t491 = t342 - t341;
t490 = g(1) * t356 + g(2) * t353;
t488 = qJDD(3) * t357;
t487 = t358 * qJDD(1);
t430 = qJD(1) * qJD(5);
t486 = -t493 * qJ(5) - t352 * t430;
t439 = qJD(3) * t355;
t293 = t305 * t439;
t304 = qJDD(2) + t487;
t295 = t352 * t304;
t344 = qJDD(3) * qJ(4);
t345 = qJD(3) * qJD(4);
t261 = t293 + t295 + t344 + t345;
t255 = -t261 + t486;
t252 = qJDD(3) * pkin(5) - t255;
t297 = t352 * t305;
t347 = qJD(3) * qJ(4);
t280 = t297 + t347;
t447 = qJD(1) * t352;
t326 = qJ(5) * t447;
t271 = -t326 - t280;
t267 = qJD(3) * pkin(5) - t271;
t333 = t355 * qJ(4);
t346 = -pkin(8) + t357;
t374 = pkin(5) * t355 + t346 * t352 - qJ(2);
t269 = t333 + t374;
t405 = qJD(3) * t458;
t272 = -qJD(5) * t355 + t352 * t405;
t329 = t355 * qJDD(1);
t412 = t352 * t431;
t288 = -qJDD(6) - t329 + t412;
t299 = t458 * t355;
t312 = qJD(6) + t446;
t308 = qJ(5) * t412;
t441 = qJD(3) * t352;
t292 = t305 * t441;
t296 = t355 * t304;
t410 = qJDD(4) + t292 - t296;
t362 = t308 + (-qJ(5) * qJDD(1) - t430) * t355 + t410;
t251 = t346 * qJDD(3) + t362;
t327 = qJ(4) * t446;
t432 = qJD(5) + t327;
t260 = t374 * qJD(1) + t432;
t399 = -qJD(6) * t260 - t251;
t485 = -(qJD(6) * t269 + t272) * t312 + t252 * t352 - t299 * t288 + (qJD(3) * t267 + t399) * t355;
t480 = pkin(3) * t352;
t408 = qJ(2) + t480;
t279 = t408 * qJD(1) - t327;
t351 = sin(qJ(6));
t354 = cos(qJ(6));
t484 = (MDP(27) * t354 - MDP(28) * t351) * t312 - t279 * MDP(17);
t475 = qJ(4) * t352;
t380 = t355 * t357 - t475;
t474 = qJDD(3) * pkin(3);
t262 = t410 - t474;
t401 = -t305 * t355 + qJD(4);
t478 = qJD(3) * pkin(3);
t278 = t401 - t478;
t483 = (t278 * t352 + t280 * t355) * qJD(3) + t261 * t352 - t262 * t355;
t466 = t312 * t351;
t377 = qJD(6) * t466 + t288 * t354;
t437 = qJD(6) * t354;
t378 = t288 * t351 - t312 * t437;
t435 = t351 * qJD(3);
t291 = t354 * t447 - t435;
t442 = qJD(3) * t291;
t440 = qJD(3) * t354;
t289 = t351 * t447 + t440;
t443 = qJD(3) * t289;
t482 = (-t378 + t443) * MDP(27) - (t377 - t442) * MDP(28);
t350 = qJ(4) + pkin(5);
t373 = t346 * t355 - t350 * t352;
t479 = g(3) * t355;
t481 = (t373 * qJD(1) + qJD(6) * t346) * t312 + t352 * t491 - t252 + t479;
t339 = g(3) * t352;
t477 = pkin(1) * qJDD(1);
t360 = qJD(1) ^ 2;
t476 = qJ(2) * t360;
t457 = t493 * t354;
t257 = -t289 * qJD(6) - qJDD(3) * t351 + t457;
t473 = t257 * t351;
t472 = t269 * t288;
t471 = t289 * t312;
t470 = t289 * t352;
t469 = t291 * t312;
t468 = t291 * t352;
t467 = t291 * t354;
t465 = t352 * t353;
t464 = t352 * t356;
t463 = t353 * t355;
t462 = t354 * t355;
t461 = t355 * t356;
t460 = t355 * t360;
t359 = qJD(3) ^ 2;
t459 = t358 * t359;
t456 = pkin(3) * t463 + qJ(4) * t465;
t455 = g(1) * t461 + g(2) * t463;
t331 = t355 * qJD(4);
t454 = -qJ(4) * t329 - qJD(1) * t331;
t453 = t356 * pkin(1) + t353 * qJ(2);
t348 = t352 ^ 2;
t349 = t355 ^ 2;
t450 = -t348 - t349;
t449 = t348 - t349;
t448 = t359 + t360;
t445 = qJD(3) * t271;
t444 = qJD(3) * t280;
t438 = qJD(6) * t352;
t420 = t357 * t352;
t392 = -qJ(2) + t420;
t270 = t392 * qJD(1) + t432;
t433 = -qJD(5) - t270;
t429 = qJDD(1) * qJ(2);
t427 = qJDD(3) * t352;
t426 = qJDD(3) * t355;
t425 = qJDD(3) * t358;
t424 = -MDP(14) + MDP(19);
t423 = MDP(15) - MDP(20);
t422 = MDP(16) + MDP(18);
t421 = 0.2e1 * qJD(1) * qJD(2);
t419 = t355 * t466;
t418 = t312 * t462;
t417 = t352 * t460;
t416 = t266 * t441 - t491;
t415 = t357 * MDP(21);
t414 = t312 * t435;
t413 = t312 * t440;
t409 = qJDD(5) - t454;
t407 = t333 - t480;
t391 = t333 + t420;
t287 = -qJ(2) + t391;
t403 = qJD(1) * t287 + t270;
t254 = t362 + t488;
t402 = -t254 - t445;
t363 = t373 * qJD(3) - qJD(2);
t248 = t363 * qJD(1) + t374 * qJDD(1) + t409;
t265 = t346 * qJD(3) + t492;
t400 = qJD(6) * t265 - t248;
t396 = qJDD(2) - t477;
t395 = pkin(3) * t465 + t356 * pkin(7) + t453;
t394 = g(1) * t463 - g(2) * t461 - t296 - t339;
t393 = -g(2) * t464 - t295 + t479;
t388 = pkin(3) * t355 + t475;
t385 = -qJDD(4) - t394;
t300 = qJ(2) - t407;
t384 = (qJD(1) * t300 + t279) * qJD(3);
t334 = t356 * qJ(2);
t382 = pkin(3) * t464 - qJ(4) * t461 + t334;
t381 = t421 + 0.2e1 * t429;
t379 = t292 - t385;
t366 = t380 * qJD(3) - qJD(2);
t253 = t366 * qJD(1) + t392 * qJDD(1) + t409;
t268 = t331 + t366;
t376 = qJD(1) * t268 + qJDD(1) * t287 + t253;
t372 = t388 * qJD(3) + qJD(2);
t371 = t381 - t459;
t370 = t379 - t474;
t256 = t372 * qJD(1) + t408 * qJDD(1) + t454;
t274 = -t331 + t372;
t367 = -qJD(1) * t274 - qJDD(1) * t300 - t256 + t459;
t365 = -g(1) * t465 + 0.2e1 * t344 + 0.2e1 * t345 - t393;
t275 = t297 + t326;
t364 = t346 * t288 + (-t267 + t275) * t312;
t328 = t349 * qJDD(1);
t325 = qJD(6) * t435;
t310 = t355 * t425;
t298 = t458 * t352;
t294 = t388 * qJD(1);
t284 = t351 * t353 - t354 * t461;
t283 = t351 * t461 + t353 * t354;
t282 = t351 * t356 + t353 * t462;
t281 = t351 * t463 - t354 * t356;
t277 = t380 * qJD(1);
t273 = qJD(5) * t352 + t355 * t405;
t259 = t331 + t363;
t258 = t351 * t428 + qJDD(3) * t354 - t325 + (t352 * t437 + t355 * t435) * qJD(1);
t250 = t260 * t351 + t265 * t354;
t249 = t260 * t354 - t265 * t351;
t247 = t354 * t248;
t1 = [(-t352 * t359 + t426) * MDP(9) + (-t355 * t359 - t427) * MDP(10) + (-0.2e1 * t352 * t411 + t328) * MDP(7) + (0.2e1 * qJ(2) * t411 + t310 + (t371 - t490) * t352) * MDP(12) + (t310 + t355 * t384 + (-t367 - t490) * t352) * MDP(14) + (-qJDD(3) * t299 + (t403 * t355 + t272) * qJD(3) + (t376 + t490) * t352) * MDP(19) + (t381 - t490) * MDP(5) + t490 * MDP(3) + t491 * MDP(2) + (qJDD(2) - t491 - 0.2e1 * t477) * MDP(4) + (t450 * t487 - t483 + t491) * MDP(15) + ((qJDD(1) * t298 - t255 + (-qJD(3) * t299 + t273) * qJD(1)) * t352 + (qJDD(1) * t299 + (qJD(3) * t298 - t272) * qJD(1) + t402) * t355 + t416) * MDP(20) + (-t249 * t441 - g(1) * t284 + g(2) * t282 + t247 * t355 + t298 * t258 + t273 * t289 + (t259 * t312 - t472 + (-t265 * t355 + t267 * t352 + t299 * t312) * qJD(6)) * t354 + t485 * t351) * MDP(27) + (t250 * t441 - g(1) * t283 - g(2) * t281 + t298 * t257 + t273 * t291 + (-(qJD(6) * t299 + t259) * t312 + t472 + t400 * t355 - t267 * t438) * t351 + t485 * t354) * MDP(28) + (t256 * t300 + t279 * t274 - g(1) * (t358 * t353 + t382) - g(2) * (-t353 * t333 + t395) + t483 * t358) * MDP(17) + (-t254 * t299 + t266 * t272 - t255 * t298 - t271 * t273 + t253 * t287 + t270 * t268 - g(1) * (pkin(4) * t464 + t382) - g(2) * (-qJ(5) * t356 + t395) + (-g(1) * t458 - g(2) * (pkin(4) * t352 - t333)) * t353) * MDP(21) + qJDD(1) * MDP(1) + (t257 * t352 * t354 + (-t351 * t438 + t354 * t439) * t291) * MDP(22) + (-t288 * t355 - t312 * t441) * MDP(26) + ((t257 + t413) * t355 + (-t377 - t442) * t352) * MDP(24) + ((-t258 - t414) * t355 + (t378 + t443) * t352) * MDP(25) + 0.2e1 * (-t352 * t329 + t449 * t431) * MDP(8) + (-t396 * pkin(1) - g(1) * (-pkin(1) * t353 + t334) - g(2) * t453 + (t421 + t429) * qJ(2)) * MDP(6) + (qJDD(3) * t298 + t376 * t355 + (-t403 * t352 + t273) * qJD(3) + t455) * MDP(18) + ((-0.2e1 * qJ(2) * t431 - t425) * t352 + t371 * t355 - t455) * MDP(13) + ((t384 + t425) * t352 + t367 * t355 + t455) * MDP(16) + ((-t289 * t354 - t291 * t351) * t439 + (-t473 - t258 * t354 + (t289 * t351 - t467) * qJD(6)) * t352) * MDP(23); qJDD(1) * MDP(4) - t360 * MDP(5) + (t396 - t491 - t476) * MDP(6) - t491 * MDP(17) + t416 * MDP(21) + t423 * (-qJDD(1) * t348 - t328) + (t270 * MDP(21) + t484) * qJD(1) + ((qJD(3) * t278 + t261) * MDP(17) - t255 * MDP(21) + (t258 - t414) * MDP(27) + (t257 - t413) * MDP(28)) * t352 + (MDP(12) - t424) * (-t448 * t352 + t426) + (-MDP(13) + t422) * (t448 * t355 + t427) + ((-t262 + t444) * MDP(17) + t402 * MDP(21) + t482) * t355; MDP(7) * t417 - t449 * t360 * MDP(8) + MDP(9) * t329 - MDP(10) * t428 + qJDD(3) * MDP(11) + (-qJ(2) * t460 - t394) * MDP(12) + ((t476 + t342) * t352 + t393) * MDP(13) + (0.2e1 * t474 + (-t279 * t355 - t294 * t352) * qJD(1) + t385) * MDP(14) + (-t388 * qJDD(1) + ((t280 - t347) * t355 + (-qJD(4) + t278 + t478) * t352) * qJD(1)) * MDP(15) + ((-t279 * t352 + t294 * t355) * qJD(1) + t365) * MDP(16) + (-t262 * pkin(3) - g(1) * t456 - g(3) * t407 + t261 * qJ(4) - t278 * t297 - t279 * t294 + t401 * t280 + t388 * t341) * MDP(17) + (-qJD(3) * t276 + t293 + (t270 * t352 - t277 * t355) * qJD(1) + t365 - t486) * MDP(18) + (-qJ(5) * t329 - qJD(3) * t275 + t308 + 0.2e1 * t488 + (-t277 * t352 + t433 * t355) * qJD(1) + t379) * MDP(19) + (-t380 * qJDD(1) + (t271 + t275 + t347) * t446) * MDP(20) + (t254 * t357 - t255 * qJ(4) - t266 * t275 - t270 * t277 - g(1) * (pkin(4) * t463 + t456) - g(3) * t391 - t492 * t271 - t380 * t341) * MDP(21) + (-t312 * t467 - t473) * MDP(22) + ((-t257 + t471) * t354 + (t258 + t469) * t351) * MDP(23) + ((-t418 + t468) * qJD(1) + t378) * MDP(24) + ((t419 - t470) * qJD(1) + t377) * MDP(25) + t312 * MDP(26) * t447 + (t249 * t447 + t350 * t258 + t289 * t492 + t364 * t351 - t354 * t481) * MDP(27) + (-t250 * t447 + t350 * t257 + t291 * t492 + t351 * t481 + t364 * t354) * MDP(28); (t370 - t444) * MDP(17) + (-qJDD(3) * pkin(4) + t308 + t370 + t445) * MDP(21) + t422 * (-t349 * t360 - t359) + t424 * (qJDD(3) - t417) + ((-MDP(21) * qJ(5) + t423) * qJDD(1) + (t433 * MDP(21) - t484) * qJD(1)) * t355 - t482; t329 * MDP(18) + (t409 + t490) * MDP(21) - t377 * MDP(27) + t378 * MDP(28) + t450 * MDP(20) * t360 + (-qJ(2) * MDP(21) + (MDP(19) + t415) * t352) * qJDD(1) + ((t266 * t355 + t271 * t352 - qJD(2)) * MDP(21) + (-t419 - t470) * MDP(27) + (-t418 - t468) * MDP(28) + ((-MDP(21) * qJ(4) - 0.2e1 * MDP(18)) * t352 + (0.2e1 * MDP(19) + t415) * t355) * qJD(3)) * qJD(1); t291 * t289 * MDP(22) + (-t289 ^ 2 + t291 ^ 2) * MDP(23) + (t457 + t471) * MDP(24) + (t325 + t469) * MDP(25) - t288 * MDP(26) + (-g(1) * t281 + g(2) * t283 + t250 * t312 - t267 * t291 + t247) * MDP(27) + (-g(1) * t282 - g(2) * t284 + t249 * t312 + t267 * t289) * MDP(28) + (-qJDD(3) * MDP(25) + (-t251 + t339) * MDP(28) + (-qJD(3) * MDP(24) - MDP(25) * t447 - MDP(27) * t265 - MDP(28) * t260) * qJD(6)) * t354 + ((-qJD(1) * t438 - qJDD(3)) * MDP(24) - t493 * MDP(25) + (t399 + t339) * MDP(27) + t400 * MDP(28)) * t351;];
tau  = t1;
