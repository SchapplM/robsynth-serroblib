% Calculate vector of inverse dynamics joint torques for
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPPRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S6PPPRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:45
% EndTime: 2019-03-08 18:40:48
% DurationCPUTime: 3.43s
% Computational Cost: add. (3223->353), mult. (8535->534), div. (0->0), fcn. (9030->18), ass. (0->186)
t350 = sin(qJ(4));
t347 = cos(pkin(6));
t326 = qJD(1) * t347 + qJD(2);
t337 = sin(pkin(14));
t338 = sin(pkin(13));
t342 = sin(pkin(6));
t343 = cos(pkin(14));
t344 = cos(pkin(13));
t346 = cos(pkin(7));
t442 = t344 * t346;
t381 = (-t337 * t338 + t343 * t442) * t342;
t341 = sin(pkin(7));
t445 = t341 * t343;
t292 = qJD(1) * t381 + t326 * t445;
t446 = t341 * t342;
t421 = t344 * t446;
t308 = -qJD(1) * t421 + t326 * t346 + qJD(3);
t340 = sin(pkin(8));
t345 = cos(pkin(8));
t399 = t292 * t345 + t308 * t340;
t480 = t399 * t350;
t451 = t338 * t343;
t380 = (t337 * t442 + t451) * t342;
t452 = t337 * t341;
t293 = qJD(1) * t380 + t326 * t452;
t353 = cos(qJ(4));
t378 = t350 * t293 - t399 * t353;
t268 = t293 * t353 + t480;
t323 = t347 * qJDD(1) + qJDD(2);
t291 = qJDD(1) * t380 + t323 * t452;
t290 = qJDD(1) * t381 + t323 * t445;
t307 = -qJDD(1) * t421 + t346 * t323 + qJDD(3);
t400 = t290 * t345 + t307 * t340;
t436 = qJD(4) * t353;
t372 = qJD(4) * t480 + t350 * t291 + t293 * t436 - t400 * t353;
t311 = t346 * t347 - t421;
t444 = t341 * t347;
t370 = t343 * t444 + t381;
t367 = t370 * t345;
t447 = t340 * t353;
t301 = t342 * t451 + (t342 * t442 + t444) * t337;
t461 = t301 * t350;
t271 = -t311 * t447 - t353 * t367 + t461;
t339 = sin(pkin(12));
t468 = cos(pkin(12));
t412 = t468 * t344;
t385 = -t339 * t338 + t347 * t412;
t414 = t342 * t468;
t365 = -t341 * t414 + t385 * t346;
t413 = t468 * t338;
t384 = t339 * t344 + t347 * t413;
t282 = t365 * t337 + t384 * t343;
t449 = t339 * t347;
t383 = -t344 * t449 - t413;
t368 = t339 * t446 + t383 * t346;
t382 = -t338 * t449 + t412;
t283 = t368 * t337 + t382 * t343;
t359 = t382 * t337 - t368 * t343;
t450 = t339 * t342;
t369 = t383 * t341 - t346 * t450;
t475 = t369 * t340 + t359 * t345;
t358 = t384 * t337 - t365 * t343;
t366 = t385 * t341 + t346 * t414;
t476 = t366 * t340 + t358 * t345;
t404 = g(1) * (t283 * t350 + t475 * t353) + g(2) * (t282 * t350 + t476 * t353);
t388 = g(3) * t271 + t404;
t479 = qJD(4) * t268 - t372 + t388;
t352 = cos(qJ(5));
t426 = qJD(4) * qJD(5);
t417 = t352 * t426;
t349 = sin(qJ(5));
t424 = qJDD(4) * t349;
t391 = -t417 - t424;
t478 = -qJD(5) * qJD(6) + t391;
t471 = -t291 * t353 - t400 * t350;
t252 = qJDD(4) * pkin(10) - qJD(4) * t378 - t471;
t266 = qJD(4) * pkin(10) + t268;
t278 = -t292 * t340 + t308 * t345;
t257 = t266 * t352 + t278 * t349;
t277 = -t290 * t340 + t307 * t345;
t465 = t277 * t352;
t239 = -qJDD(5) * pkin(5) + qJD(5) * t257 + t252 * t349 - t465;
t437 = qJD(4) * t352;
t327 = -qJD(6) + t437;
t360 = t311 * t340 + t367;
t272 = t301 * t353 + t360 * t350;
t361 = t311 * t345 - t370 * t340;
t258 = t272 * t349 - t361 * t352;
t261 = t282 * t353 - t476 * t350;
t263 = t283 * t353 - t475 * t350;
t275 = t358 * t340 - t366 * t345;
t276 = t359 * t340 - t369 * t345;
t389 = g(1) * (-t263 * t349 + t276 * t352) + g(2) * (-t261 * t349 + t275 * t352) - g(3) * t258;
t406 = pkin(5) * t349 - pkin(11) * t352;
t474 = (pkin(11) * qJD(6) + t406 * qJD(4)) * t327 - t239 - t389;
t423 = t352 * qJDD(4);
t390 = -t349 * t426 + t423;
t314 = qJDD(6) - t390;
t319 = t406 * qJD(5);
t321 = -pkin(5) * t352 - pkin(11) * t349 - pkin(4);
t473 = (t268 - t319) * t327 + t321 * t314;
t443 = t343 * t345;
t472 = (-t337 * t350 + t353 * t443) * t341 + t346 * t447;
t255 = qJD(5) * pkin(11) + t257;
t470 = (pkin(10) * t327 + t255) * qJD(6) + t388;
t469 = qJD(4) * pkin(4);
t466 = t277 * t349;
t348 = sin(qJ(6));
t351 = cos(qJ(6));
t407 = t348 * qJDD(5) - t478 * t351;
t429 = qJD(6) * t349;
t416 = qJD(4) * t429;
t297 = -t348 * t416 + t407;
t462 = t297 * t348;
t459 = t314 * t348;
t458 = t314 * t351;
t427 = t351 * qJD(5);
t439 = qJD(4) * t349;
t315 = t348 * t439 - t427;
t457 = t315 * t327;
t433 = qJD(5) * t348;
t317 = t351 * t439 + t433;
t456 = t317 * t327;
t455 = t317 * t351;
t453 = t327 * t352;
t448 = t340 * t350;
t335 = t349 ^ 2;
t440 = -t352 ^ 2 + t335;
t435 = qJD(5) * t315;
t434 = qJD(5) * t317;
t432 = qJD(5) * t349;
t431 = qJD(5) * t352;
t430 = qJD(6) * t348;
t428 = qJD(6) * t351;
t420 = t340 * t436;
t419 = t327 * t433;
t418 = t327 * t427;
t265 = t378 - t469;
t410 = -qJD(4) * t265 - t252;
t264 = t321 * qJD(4) + t378;
t403 = t266 * t349 - t278 * t352;
t409 = -qJDD(5) * pkin(11) + t403 * qJD(5) - qJD(6) * t264 - t252 * t352 - t466;
t405 = g(1) * (t263 * t352 + t276 * t349) + g(2) * (t261 * t352 + t275 * t349);
t259 = t272 * t352 + t361 * t349;
t247 = t259 * t351 + t271 * t348;
t246 = -t259 * t348 + t271 * t351;
t304 = t346 * t448 + (t337 * t353 + t350 * t443) * t341;
t310 = -t340 * t445 + t345 * t346;
t285 = t304 * t352 + t310 * t349;
t402 = t285 * t351 - t348 * t472;
t401 = -t285 * t348 - t351 * t472;
t284 = t304 * t349 - t310 * t352;
t313 = t345 * t349 + t352 * t448;
t312 = -t345 * t352 + t349 * t448;
t397 = -t327 * t428 + t459;
t396 = t327 * t430 + t458;
t395 = -MDP(12) * t352 + MDP(13) * t349 - MDP(5);
t394 = -g(1) * t450 + g(2) * t414 - g(3) * t347;
t387 = g(1) * t263 + g(2) * t261 + g(3) * t272;
t386 = -qJD(6) * t255 - t404;
t377 = t405 + t409;
t376 = -pkin(10) * qJDD(5) + (t265 - t378 - t469) * qJD(5);
t374 = qJD(6) * t321 * t327 - t387;
t254 = -qJD(5) * pkin(5) + t403;
t373 = -pkin(11) * t314 + (-t254 + t403) * t327;
t371 = -pkin(10) * t314 + qJD(5) * t254 + t327 * t378 - t409;
t354 = qJD(5) ^ 2;
t363 = 0.2e1 * qJDD(4) * pkin(4) - pkin(10) * t354 + t479;
t355 = qJD(4) ^ 2;
t331 = t351 * qJDD(5);
t306 = t313 * qJD(5) + t349 * t420;
t305 = -t312 * qJD(5) + t352 * t420;
t300 = t304 * qJD(4);
t299 = t472 * qJD(4);
t298 = t351 * t416 - t331 + (t424 + (qJD(6) + t437) * qJD(5)) * t348;
t274 = t285 * qJD(5) + t299 * t349;
t273 = -t284 * qJD(5) + t299 * t352;
t270 = t272 * qJD(4);
t269 = (t360 * t353 - t461) * qJD(4);
t245 = -t258 * qJD(5) + t269 * t352;
t244 = t259 * qJD(5) + t269 * t349;
t243 = qJD(4) * t319 + t321 * qJDD(4) + t372;
t242 = t351 * t243;
t241 = t255 * t351 + t264 * t348;
t240 = -t255 * t348 + t264 * t351;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t323 * t347 - g(3) + (t338 ^ 2 + t344 ^ 2) * t342 ^ 2 * qJDD(1)) * MDP(2) + (t290 * t370 + t291 * t301 + t307 * t311 - g(3)) * MDP(3) + (-qJD(5) * t244 - qJDD(5) * t258) * MDP(12) + (-qJD(5) * t245 - qJDD(5) * t259) * MDP(13) + (-(-t247 * qJD(6) - t245 * t348 + t270 * t351) * t327 + t246 * t314 + t244 * t315 + t258 * t298) * MDP(19) + ((t246 * qJD(6) + t245 * t351 + t270 * t348) * t327 - t247 * t314 + t244 * t317 + t258 * t297) * MDP(20) + (-t272 * MDP(6) + t395 * t271) * qJDD(4) + (-t270 * MDP(5) - t269 * MDP(6) + (-t270 * t352 + t271 * t432) * MDP(12) + (t270 * t349 + t271 * t431) * MDP(13)) * qJD(4); (t394 + t323) * MDP(2) + (t307 * t346 + (t290 * t343 + t291 * t337) * t341 + t394) * MDP(3) + (-qJD(5) * t274 - qJDD(5) * t284) * MDP(12) + (-qJD(5) * t273 - qJDD(5) * t285) * MDP(13) + (-(-t402 * qJD(6) - t273 * t348 + t300 * t351) * t327 + t401 * t314 + t274 * t315 + t284 * t298) * MDP(19) + ((t401 * qJD(6) + t273 * t351 + t300 * t348) * t327 - t402 * t314 + t274 * t317 + t284 * t297) * MDP(20) + (-MDP(6) * t304 - t395 * t472) * qJDD(4) + (-t300 * MDP(5) - t299 * MDP(6) + (-t300 * t352 - t432 * t472) * MDP(12) + (t300 * t349 - t431 * t472) * MDP(13)) * qJD(4); (g(1) * t369 + g(2) * t366 - g(3) * t311 + t307) * MDP(3) + (-qJD(5) * t306 - qJDD(5) * t312) * MDP(12) + (-qJD(5) * t305 - qJDD(5) * t313) * MDP(13) + (-(-t305 * t348 - t313 * t428) * t327 - t313 * t459 + t306 * t315 + t312 * t298) * MDP(19) + ((t305 * t351 - t313 * t430) * t327 - t313 * t458 + t306 * t317 + t312 * t297) * MDP(20) + ((-qJDD(4) * MDP(6) + (-MDP(19) * t351 + MDP(20) * t348) * t327 * qJD(4) + t395 * t355) * t350 + (t390 * MDP(12) + t391 * MDP(13) - t396 * MDP(19) + t397 * MDP(20) + qJDD(4) * MDP(5) - t355 * MDP(6)) * t353) * t340; qJDD(4) * MDP(4) + t479 * MDP(5) + (t387 + t471) * MDP(6) + (qJDD(4) * t335 + 0.2e1 * t349 * t417) * MDP(7) + 0.2e1 * (t349 * t423 - t440 * t426) * MDP(8) + (qJDD(5) * t349 + t352 * t354) * MDP(9) + (qJDD(5) * t352 - t349 * t354) * MDP(10) + (t376 * t349 + t363 * t352) * MDP(12) + (-t363 * t349 + t376 * t352) * MDP(13) + (t297 * t349 * t351 + (-t348 * t429 + t352 * t427) * t317) * MDP(14) + ((-t315 * t351 - t317 * t348) * t431 + (-t462 - t298 * t351 + (t315 * t348 - t455) * qJD(6)) * t349) * MDP(15) + ((-t297 - t418) * t352 + (t396 + t434) * t349) * MDP(16) + ((t298 + t419) * t352 + (-t397 - t435) * t349) * MDP(17) + (-t314 * t352 - t327 * t432) * MDP(18) + (t473 * t351 + t374 * t348 + (pkin(10) * t435 + t371 * t348 + t470 * t351 - t242) * t352 + (t254 * t428 + t240 * qJD(5) + t239 * t348 + t378 * t315 + (t298 - t419) * pkin(10)) * t349) * MDP(19) + (-t473 * t348 + t374 * t351 + (pkin(10) * t434 + t371 * t351 + (t243 - t470) * t348) * t352 + (-t254 * t430 - t241 * qJD(5) + t239 * t351 + t378 * t317 + (t297 - t418) * pkin(10)) * t349) * MDP(20); MDP(9) * t424 + MDP(10) * t423 + qJDD(5) * MDP(11) + (t410 * t349 - t389 + t465) * MDP(12) + (g(3) * t259 + t410 * t352 + t405 - t466) * MDP(13) + (-t327 * t455 + t462) * MDP(14) + ((t297 + t457) * t351 + (-t298 + t456) * t348) * MDP(15) + ((-t317 * t349 + t351 * t453) * qJD(4) + t397) * MDP(16) + ((t315 * t349 - t348 * t453) * qJD(4) + t396) * MDP(17) + t327 * MDP(18) * t439 + (-pkin(5) * t298 - t240 * t439 - t257 * t315 + t373 * t348 + t474 * t351) * MDP(19) + (-pkin(5) * t297 + t241 * t439 - t257 * t317 - t474 * t348 + t373 * t351) * MDP(20) + (-t349 * t352 * MDP(7) + t440 * MDP(8)) * t355; t315 * t317 * MDP(14) + (-t315 ^ 2 + t317 ^ 2) * MDP(15) + (t407 - t457) * MDP(16) + (t331 - t456) * MDP(17) + t314 * MDP(18) + (-g(3) * t246 - t241 * t327 - t254 * t317 + t242) * MDP(19) + (g(3) * t247 - t240 * t327 + t254 * t315) * MDP(20) + (-MDP(17) * t416 + t386 * MDP(19) + t377 * MDP(20)) * t351 + (-MDP(16) * t416 + t478 * MDP(17) + t377 * MDP(19) + (-t243 - t386) * MDP(20)) * t348;];
tau  = t1;
