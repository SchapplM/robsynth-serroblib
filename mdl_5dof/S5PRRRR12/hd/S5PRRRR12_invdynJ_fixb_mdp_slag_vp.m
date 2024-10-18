% Calculate vector of inverse dynamics joint torques for
% S5PRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR12_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:23
% EndTime: 2024-09-28 18:08:25
% DurationCPUTime: 0.89s
% Computational Cost: add. (2471->323), mult. (4810->497), div. (0->0), fcn. (3914->26), ass. (0->190)
t369 = sin(pkin(5));
t489 = g(3) * t369;
t376 = sin(qJ(2));
t428 = qJDD(1) * t376;
t380 = cos(qJ(2));
t439 = qJD(2) * t380;
t395 = qJD(1) * t439 + t428;
t488 = t395 * t369;
t375 = sin(qJ(3));
t379 = cos(qJ(3));
t312 = (-t375 * t376 + t379 * t380) * t369;
t403 = t375 * t380 + t376 * t379;
t313 = t403 * t369;
t367 = sin(pkin(11));
t370 = cos(pkin(11));
t487 = (g(1) * t367 - g(2) * t370) * t369;
t442 = qJD(1) * t369;
t326 = qJD(2) * pkin(2) + t380 * t442;
t440 = qJD(1) * t376;
t419 = t369 * t440;
t297 = t379 * t326 - t375 * t419;
t374 = sin(qJ(4));
t436 = qJD(4) * t374;
t298 = t326 * t375 + t379 * t419;
t378 = cos(qJ(4));
t448 = t378 * t298;
t486 = pkin(3) * t436 - t297 * t374 - t448;
t435 = qJD(4) * t378;
t471 = t298 * t374;
t485 = -pkin(3) * t435 + t297 * t378 - t471;
t363 = qJD(2) + qJD(3);
t366 = qJ(2) + qJ(3);
t355 = pkin(5) + t366;
t345 = qJ(4) + t355;
t335 = cos(t345) / 0.2e1;
t356 = pkin(5) - t366;
t346 = -qJ(4) + t356;
t340 = cos(t346);
t318 = t340 / 0.2e1 + t335;
t334 = sin(t346) / 0.2e1;
t339 = sin(t345);
t360 = qJ(4) + t366;
t347 = sin(t360);
t484 = -g(3) * (t339 / 0.2e1 + t334) - g(2) * (t318 * t370 - t347 * t367) - g(1) * (-t318 * t367 - t347 * t370);
t276 = t312 * t378 - t313 * t374;
t368 = sin(pkin(6));
t371 = cos(pkin(6));
t372 = cos(pkin(5));
t269 = -t276 * t368 + t371 * t372;
t354 = qJD(4) + t363;
t463 = t354 * t371;
t329 = qJD(5) + t463;
t404 = t276 * t371 + t368 * t372;
t464 = t354 * t368;
t483 = t269 * t464 - t329 * t404;
t460 = t369 * t380;
t336 = qJDD(1) * t460;
t303 = qJDD(2) * pkin(2) - qJD(2) * t419 + t336;
t299 = t379 * t303;
t361 = qJDD(2) + qJDD(3);
t438 = qJD(3) * t326;
t418 = t375 * t438;
t437 = qJD(3) * t379;
t259 = -t418 + pkin(3) * t361 + t299 + (-t375 * t428 + (-t375 * t439 - t376 * t437) * qJD(1)) * t369;
t416 = qJD(3) * t440;
t409 = t369 * t416;
t325 = t375 * t409;
t470 = t303 * t375;
t260 = t470 - t325 + (t438 + t488) * t379;
t284 = pkin(3) * t363 + t297;
t411 = qJD(4) * t284 + t260;
t482 = -t374 * t259 - t378 * t411;
t481 = pkin(2) * t361;
t353 = qJDD(4) + t361;
t480 = pkin(3) * t353;
t479 = pkin(4) * t353;
t279 = t363 * t312;
t280 = t363 * t313;
t253 = qJD(4) * t276 + t279 * t378 - t280 * t374;
t475 = t253 * t329;
t266 = -t374 * t284 - t448;
t261 = pkin(10) * t464 - t266;
t373 = sin(qJ(5));
t474 = t261 * t373;
t377 = cos(qJ(5));
t473 = t261 * t377;
t472 = t266 * t354;
t350 = pkin(2) * t379 + pkin(3);
t451 = t374 * t375;
t408 = -pkin(2) * t451 + t378 * t350;
t314 = pkin(4) + t408;
t469 = t314 * t353;
t466 = t353 * t371;
t327 = qJDD(5) + t466;
t468 = t327 * t371;
t467 = t353 * t368;
t362 = t368 ^ 2;
t465 = t354 * t362;
t462 = t368 * t373;
t461 = t368 * t377;
t459 = t371 * t373;
t458 = t371 * t377;
t457 = t372 * t373;
t456 = t372 * t376;
t455 = t372 * t377;
t454 = t372 * t380;
t453 = t373 * t377;
t450 = t375 * t378;
t289 = t298 * t436;
t251 = pkin(10) * t467 - t289 - t482;
t449 = t377 * t251;
t447 = qJDD(1) - g(3);
t300 = t403 * t442;
t301 = qJD(1) * t312;
t446 = -t300 * t374 + t301 * t378 - t350 * t435 - (-t375 * t436 + (t378 * t379 - t451) * qJD(3)) * pkin(2);
t445 = t300 * t378 + t301 * t374 - t350 * t436 + (-t375 * t435 + (-t374 * t379 - t450) * qJD(3)) * pkin(2);
t444 = pkin(2) * t450 + t374 * t350;
t364 = t373 ^ 2;
t443 = -t377 ^ 2 + t364;
t441 = qJD(1) * t372;
t434 = qJD(5) * t354;
t433 = qJD(5) * t373;
t432 = qJD(5) * t377;
t431 = t327 * MDP(15);
t430 = qJD(5) - t329;
t429 = qJDD(1) * t372;
t427 = pkin(4) * t465;
t424 = t269 * t467;
t349 = pkin(3) * t378 + pkin(4);
t422 = t362 * t349 * t353;
t421 = t371 * t457;
t420 = t371 * t455;
t265 = t378 * t284 - t471;
t262 = pkin(4) * t354 + t265;
t258 = -t262 * t368 + t371 * t441;
t417 = qJD(5) * t258 * t368;
t415 = t368 * t429;
t414 = t445 * t354;
t413 = t327 + t466;
t412 = t329 + t463;
t410 = qJD(1) * t430;
t257 = t378 * t259;
t388 = qJD(4) * t266 - t374 * t260 + t257;
t252 = t388 + t479;
t406 = -t373 * t251 + t252 * t458 + t377 * t415;
t277 = t312 * t374 + t313 * t378;
t401 = t353 * MDP(8) + t371 * t431 + ((t377 * t413 - t412 * t433) * MDP(14) + (t373 * t413 + t412 * t432) * MDP(13)) * t368 + (0.2e1 * (t353 * t453 - t434 * t443) * MDP(12) + (0.2e1 * t354 * t373 * t432 + t353 * t364) * MDP(11)) * t362;
t400 = t262 * t458 - t474;
t399 = -t262 * t459 - t473;
t398 = t262 * t371 + t368 * t441;
t397 = pkin(4) * t458 - pkin(10) * t462;
t396 = pkin(4) * t459 + pkin(10) * t461;
t394 = t361 * MDP(5) + t401;
t393 = -t314 * t432 - t373 * t445;
t316 = t334 - t339 / 0.2e1;
t348 = cos(t360);
t392 = -g(1) * (-t316 * t367 - t348 * t370) - g(2) * (t316 * t370 - t348 * t367) - g(3) * (t335 - t340 / 0.2e1) + t289;
t338 = cos(t355) / 0.2e1;
t344 = cos(t356);
t324 = t344 / 0.2e1 + t338;
t337 = sin(t356) / 0.2e1;
t343 = sin(t355);
t357 = sin(t366);
t391 = -g(1) * (-t324 * t367 - t357 * t370) - g(2) * (t324 * t370 - t357 * t367) - g(3) * (t343 / 0.2e1 + t337) + t299;
t322 = t337 - t343 / 0.2e1;
t358 = cos(t366);
t390 = -g(1) * (-t322 * t367 - t358 * t370) - g(2) * (t322 * t370 - t358 * t367) - g(3) * (t338 - t344 / 0.2e1) + t325;
t304 = t367 * t421 - t370 * t377;
t306 = t367 * t377 + t370 * t421;
t308 = t367 * t459 - t370 * t455;
t310 = t367 * t455 + t370 * t459;
t389 = -g(1) * (t304 * t347 - t310 * t348) - g(2) * (-t306 * t347 - t308 * t348) - (-t347 * t459 + t348 * t377) * t489 + ((-t373 * t398 - t473) * qJD(5) + t406) * t371 + t373 * t417;
t249 = -t252 * t368 + t371 * t429;
t305 = -t367 * t420 - t370 * t373;
t307 = -t367 * t373 + t370 * t420;
t309 = -t367 * t458 - t370 * t457;
t311 = t367 * t457 - t370 * t458;
t387 = -g(1) * (-t305 * t347 + t311 * t348) - g(2) * (-t307 * t347 + t309 * t348) - (-t347 * t458 - t348 * t373) * t489 - (t449 + (t371 * t252 + t415) * t373 + (t377 * t398 - t474) * qJD(5)) * t371 + t249 * t462 + t377 * t417;
t386 = -t249 * t461 + t389;
t385 = -t488 + (-pkin(2) * t363 - t326) * qJD(3);
t384 = t392 + t482;
t383 = t388 + t484;
t382 = (-t349 * t433 - t377 * t486) * MDP(16) + (-t349 * t432 + t373 * t486) * MDP(17);
t381 = qJD(2) ^ 2;
t359 = t368 * pkin(10);
t331 = t368 * t455;
t330 = pkin(3) * t374 + t359;
t302 = t359 + t444;
t254 = -qJD(4) * t277 - t279 * t374 - t280 * t378;
t1 = [t447 * MDP(1) + (-t280 * t363 + t312 * t361) * MDP(6) + (-t279 * t363 - t313 * t361) * MDP(7) + (t254 * t354 + t276 * t353) * MDP(9) + (-t253 * t354 - t277 * t353) * MDP(10) + (-t373 * t475 + (-t277 * t373 + t331) * t327 + (-t424 + t276 * t468 + (t329 * t371 + t465) * t254) * t377 + (-t377 * t277 * t329 + t373 * t483) * qJD(5)) * MDP(16) + ((qJD(5) * t483 - t277 * t327 - t475) * t377 + (-(-qJD(5) * t277 + t254 * t371) * t329 - t404 * t327 - t254 * t465 + t424) * t373) * MDP(17) + ((qJDD(2) * t380 - t376 * t381) * MDP(3) + (-qJDD(2) * t376 - t380 * t381) * MDP(4)) * t369; qJDD(2) * MDP(2) + (t336 - g(1) * (-t367 * t454 - t370 * t376) - g(2) * (-t367 * t376 + t370 * t454) - g(3) * t460) * MDP(3) + (-g(1) * (t367 * t456 - t370 * t380) - g(2) * (-t367 * t380 - t370 * t456) - t447 * t376 * t369) * MDP(4) + (t300 * t363 + (-t409 + t481) * t379 + t385 * t375 + t391) * MDP(6) + (t301 * t363 + (-t303 - t481) * t375 + t385 * t379 + t390) * MDP(7) + (t353 * t408 + t383 + t414) * MDP(9) + (-t353 * t444 + t446 * t354 + t384) * MDP(10) + ((-t362 * t314 * t434 - t302 * t327 + (-qJD(5) * t314 * t371 + t446) * t329) * t373 + (t314 * t468 - t249 * t368 + (-qJD(5) * t302 + t371 * t445) * t329 + (t414 + t469) * t362) * t377 + t389) * MDP(16) + (-(t302 * t377 + t314 * t459) * t327 + (t302 * t433 + t371 * t393 + t377 * t446) * t329 + (t393 * t354 - t373 * t469) * t362 + t387) * MDP(17) + t394; (t297 * t363 - t326 * t437 + t390 - t470) * MDP(7) + (t298 * t363 + t391 - t418) * MDP(6) + t392 * MDP(10) + (-(t330 * t377 + t349 * t459) * t327 - t373 * t422 + t387) * MDP(17) + ((-t330 * t373 + t349 * t458) * t327 + t377 * t422 + t386) * MDP(16) + ((-t330 * t432 + t373 * t485) * MDP(16) + (t330 * t433 + t377 * t485) * MDP(17) + t382 * t371) * t329 + (MDP(10) * t485 - MDP(9) * t486 + t382 * t362) * t354 + t394 + (-t395 * MDP(6) * t375 + (-MDP(6) * t416 - MDP(7) * t395) * t379) * t369 + (t257 + t484) * MDP(9) + ((-qJD(4) * t298 + t480) * MDP(9) - t411 * MDP(10)) * t378 + (-t411 * MDP(9) + (-t259 - t480) * MDP(10)) * t374; (t383 - t472) * MDP(9) + (t265 * t354 + t384) * MDP(10) + (t397 * t327 - (-t265 * t373 + t266 * t458) * t329 + t386) * MDP(16) + (-t396 * t327 + (t265 * t377 + t266 * t459) * t329 + t387) * MDP(17) + ((-t329 * t396 - t373 * t427) * MDP(16) + (-t329 * t397 - t377 * t427) * MDP(17)) * qJD(5) + t401 + (t377 * MDP(16) - t373 * MDP(17)) * t362 * (-t472 + t479); t431 + (-t399 * t329 - g(1) * (t305 * t348 + t311 * t347) - g(2) * (t307 * t348 + t309 * t347) - g(3) * (t331 + (-t347 * t373 + t348 * t458) * t369) + t406) * MDP(16) + (-t449 - t252 * t459 + t400 * t329 - g(1) * (t304 * t348 + t310 * t347) - g(2) * (-t306 * t348 + t308 * t347) - (-t347 * t377 - t348 * t459) * t489) * MDP(17) + (MDP(16) * t399 - MDP(17) * t400) * qJD(5) + (-MDP(11) * t453 + MDP(12) * t443) * t362 * t354 ^ 2 + ((t353 * MDP(14) - MDP(16) * t487 - t372 * MDP(17) * t410 + (MDP(13) * t430 - t258 * MDP(17)) * t354) * t377 + (t353 * MDP(13) + MDP(17) * t487 + (-MDP(16) * t410 - MDP(17) * t447) * t372 + (-MDP(14) * t430 - t258 * MDP(16)) * t354) * t373) * t368;];
tau = t1;
