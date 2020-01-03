% Calculate vector of inverse dynamics joint torques for
% S5RRPPR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:17
% EndTime: 2019-12-31 19:33:29
% DurationCPUTime: 6.59s
% Computational Cost: add. (3230->395), mult. (7624->534), div. (0->0), fcn. (5617->14), ass. (0->173)
t394 = cos(qJ(2));
t452 = cos(pkin(8));
t421 = t452 * t394;
t368 = qJD(1) * t421;
t387 = sin(pkin(8));
t391 = sin(qJ(2));
t435 = qJD(1) * t391;
t340 = t387 * t435 - t368;
t336 = qJD(5) + t340;
t356 = t387 * t394 + t391 * t452;
t343 = t356 * qJD(1);
t386 = sin(pkin(9));
t388 = cos(pkin(9));
t328 = t388 * qJD(2) - t343 * t386;
t393 = cos(qJ(5));
t327 = qJD(2) * t386 + t343 * t388;
t390 = sin(qJ(5));
t445 = t327 * t390;
t470 = t328 * t393 - t445;
t471 = t470 * t336;
t454 = qJ(3) + pkin(6);
t423 = qJD(2) * t454;
t338 = -qJD(3) * t391 - t394 * t423;
t363 = t454 * t391;
t308 = qJDD(2) * pkin(2) + qJD(1) * t338 - qJDD(1) * t363;
t337 = qJD(3) * t394 - t391 * t423;
t364 = t454 * t394;
t317 = qJD(1) * t337 + qJDD(1) * t364;
t267 = t308 * t452 - t387 * t317;
t266 = -qJDD(2) * pkin(3) + qJDD(4) - t267;
t383 = qJ(2) + pkin(8);
t378 = sin(t383);
t380 = cos(t383);
t392 = sin(qJ(1));
t395 = cos(qJ(1));
t418 = g(1) * t395 + g(2) * t392;
t402 = -g(3) * t380 + t378 * t418;
t466 = t266 - t402;
t409 = -t327 * t393 - t328 * t390;
t468 = t336 * t409;
t357 = t386 * t393 + t388 * t390;
t347 = t357 * qJD(5);
t437 = t357 * t340 + t347;
t465 = g(1) * t392 - g(2) * t395;
t467 = t380 * t465;
t342 = t356 * qJD(2);
t430 = qJDD(1) * t391;
t314 = qJD(1) * t342 - qJDD(1) * t421 + t387 * t430;
t309 = qJDD(5) + t314;
t355 = t386 * t390 - t393 * t388;
t438 = t336 * t355;
t464 = -t309 * t357 + t336 * t438;
t463 = pkin(2) * t394;
t462 = pkin(7) * t388;
t458 = g(3) * t378;
t456 = g(3) * t394;
t370 = pkin(2) * t387 + qJ(4);
t455 = pkin(7) + t370;
t453 = qJD(2) * pkin(2);
t451 = qJDD(1) * pkin(1);
t450 = t470 * t343;
t449 = t409 * t343;
t447 = t314 * t386;
t446 = t314 * t388;
t444 = t340 * t386;
t405 = -t387 * t391 + t421;
t345 = t405 * qJD(2);
t443 = t345 * t386;
t442 = t356 * t386;
t441 = t356 * t388;
t440 = t380 * t392;
t439 = t380 * t395;
t360 = qJD(1) * t364;
t348 = t387 * t360;
t431 = qJD(1) * qJD(2);
t425 = t391 * t431;
t315 = qJD(2) * t368 + qJDD(1) * t356 - t387 * t425;
t374 = pkin(1) + t463;
t428 = pkin(2) * t425 + qJDD(3);
t403 = -qJDD(1) * t374 + t428;
t252 = pkin(3) * t314 - qJ(4) * t315 - qJD(4) * t343 + t403;
t268 = t387 * t308 + t452 * t317;
t263 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t268;
t241 = t386 * t252 + t388 * t263;
t427 = t391 * t453;
t278 = pkin(3) * t342 - qJ(4) * t345 - qJD(4) * t356 + t427;
t298 = t337 * t452 + t387 * t338;
t254 = t386 * t278 + t388 * t298;
t362 = -qJD(1) * t374 + qJD(3);
t286 = pkin(3) * t340 - qJ(4) * t343 + t362;
t359 = qJD(1) * t363;
t353 = -t359 + t453;
t422 = t452 * t360;
t313 = t387 * t353 + t422;
t305 = qJD(2) * qJ(4) + t313;
t258 = t386 * t286 + t388 * t305;
t296 = pkin(2) * t435 + pkin(3) * t343 + qJ(4) * t340;
t319 = -t359 * t452 - t348;
t265 = t386 * t296 + t388 * t319;
t311 = -pkin(3) * t405 - qJ(4) * t356 - t374;
t324 = -t387 * t363 + t364 * t452;
t270 = t386 * t311 + t388 * t324;
t384 = t391 ^ 2;
t436 = -t394 ^ 2 + t384;
t434 = qJD(5) * t390;
t433 = qJD(5) * t393;
t312 = t353 * t452 - t348;
t299 = -qJD(2) * pkin(3) + qJD(4) - t312;
t432 = -qJD(4) + t299;
t429 = qJDD(1) * t394;
t293 = -t388 * qJDD(2) + t315 * t386;
t294 = qJDD(2) * t386 + t315 * t388;
t426 = -t390 * t293 + t393 * t294 + t328 * t433;
t240 = t388 * t252 - t263 * t386;
t236 = pkin(4) * t314 - pkin(7) * t294 + t240;
t239 = -pkin(7) * t293 + t241;
t420 = t393 * t236 - t239 * t390;
t253 = t388 * t278 - t298 * t386;
t257 = t388 * t286 - t305 * t386;
t419 = t393 * t293 + t390 * t294;
t264 = t388 * t296 - t319 * t386;
t269 = t388 * t311 - t324 * t386;
t297 = t337 * t387 - t452 * t338;
t318 = -t359 * t387 + t422;
t323 = t452 * t363 + t364 * t387;
t373 = -pkin(2) * t452 - pkin(3);
t416 = -t355 * t309 - t336 * t437;
t415 = pkin(3) * t380 + qJ(4) * t378;
t414 = t236 * t390 + t239 * t393;
t413 = -t240 * t388 - t241 * t386;
t245 = pkin(4) * t340 - pkin(7) * t327 + t257;
t249 = pkin(7) * t328 + t258;
t237 = t245 * t393 - t249 * t390;
t238 = t245 * t390 + t249 * t393;
t255 = -pkin(4) * t405 - pkin(7) * t441 + t269;
t259 = -pkin(7) * t442 + t270;
t412 = t255 * t393 - t259 * t390;
t411 = t255 * t390 + t259 * t393;
t410 = -t257 * t386 + t258 * t388;
t408 = -0.2e1 * pkin(1) * t431 - pkin(6) * qJDD(2);
t352 = t455 * t388;
t407 = pkin(4) * t343 + qJD(4) * t386 + qJD(5) * t352 + t340 * t462 + t264;
t351 = t455 * t386;
t406 = pkin(7) * t444 - qJD(4) * t388 + qJD(5) * t351 + t265;
t242 = -t327 * t434 + t426;
t396 = qJD(2) ^ 2;
t400 = -pkin(6) * t396 + 0.2e1 * t451 + t465;
t397 = qJD(1) ^ 2;
t399 = pkin(1) * t397 - pkin(6) * qJDD(1) + t418;
t398 = t266 * t356 + t299 * t345 - t418;
t243 = -qJD(5) * t409 + t419;
t382 = pkin(9) + qJ(5);
t379 = cos(t382);
t377 = sin(t382);
t366 = t395 * t374;
t361 = -t388 * pkin(4) + t373;
t339 = t340 ^ 2;
t332 = t377 * t392 + t379 * t439;
t331 = -t377 * t439 + t379 * t392;
t330 = t377 * t395 - t379 * t440;
t329 = t377 * t440 + t379 * t395;
t301 = t355 * t356;
t300 = t357 * t356;
t292 = pkin(4) * t442 + t323;
t279 = -pkin(4) * t444 + t318;
t272 = pkin(4) * t443 + t297;
t271 = -pkin(4) * t328 + t299;
t262 = t345 * t357 + t433 * t441 - t434 * t442;
t261 = -t345 * t355 - t347 * t356;
t247 = t293 * pkin(4) + t266;
t246 = -pkin(7) * t443 + t254;
t244 = pkin(4) * t342 - t345 * t462 + t253;
t1 = [qJDD(1) * MDP(1) + t465 * MDP(2) + t418 * MDP(3) + (qJDD(1) * t384 + 0.2e1 * t394 * t425) * MDP(4) + 0.2e1 * (t391 * t429 - t431 * t436) * MDP(5) + (qJDD(2) * t391 + t394 * t396) * MDP(6) + (qJDD(2) * t394 - t391 * t396) * MDP(7) + (t391 * t408 + t394 * t400) * MDP(9) + (-t391 * t400 + t394 * t408) * MDP(10) + (-t267 * t356 + t268 * t405 + t297 * t343 - t298 * t340 - t312 * t345 - t313 * t342 - t314 * t324 + t315 * t323 - t418) * MDP(11) + (t268 * t324 + t313 * t298 - t267 * t323 - t312 * t297 - t403 * t374 + t362 * t427 - g(1) * (-t374 * t392 + t395 * t454) - g(2) * (t392 * t454 + t366)) * MDP(12) + (-t240 * t405 + t253 * t340 + t257 * t342 + t269 * t314 + t323 * t293 - t297 * t328 + t386 * t398 + t388 * t467) * MDP(13) + (t241 * t405 - t254 * t340 - t258 * t342 - t270 * t314 + t323 * t294 + t297 * t327 - t386 * t467 + t388 * t398) * MDP(14) + (-t253 * t327 + t254 * t328 - t269 * t294 - t270 * t293 + t465 * t378 + t413 * t356 + (-t257 * t388 - t258 * t386) * t345) * MDP(15) + (-g(2) * t366 + t240 * t269 + t241 * t270 + t257 * t253 + t258 * t254 + t266 * t323 + t299 * t297 + (-g(1) * t454 - g(2) * t415) * t395 + (-g(1) * (-t374 - t415) - g(2) * t454) * t392) * MDP(16) + (-t242 * t301 - t261 * t409) * MDP(17) + (-t242 * t300 + t243 * t301 + t261 * t470 + t262 * t409) * MDP(18) + (-t242 * t405 + t261 * t336 - t301 * t309 - t342 * t409) * MDP(19) + (t243 * t405 - t262 * t336 - t300 * t309 + t342 * t470) * MDP(20) + (-t309 * t405 + t336 * t342) * MDP(21) + ((t244 * t393 - t246 * t390) * t336 + t412 * t309 - t420 * t405 + t237 * t342 - t272 * t470 + t292 * t243 + t247 * t300 + t271 * t262 - g(1) * t330 - g(2) * t332 + (t238 * t405 - t336 * t411) * qJD(5)) * MDP(22) + (-(t244 * t390 + t246 * t393) * t336 - t411 * t309 + t414 * t405 - t238 * t342 - t272 * t409 + t292 * t242 - t247 * t301 + t271 * t261 - g(1) * t329 - g(2) * t331 + (t237 * t405 - t336 * t412) * qJD(5)) * MDP(23); MDP(6) * t430 + MDP(7) * t429 + qJDD(2) * MDP(8) + (t391 * t399 - t456) * MDP(9) + (g(3) * t391 + t394 * t399) * MDP(10) + ((t313 - t318) * t343 + (-t312 + t319) * t340 + (-t314 * t387 - t315 * t452) * pkin(2)) * MDP(11) + (t312 * t318 - t313 * t319 + (t452 * t267 - t456 + t268 * t387 + (-qJD(1) * t362 + t418) * t391) * pkin(2)) * MDP(12) + (-t370 * t447 - t257 * t343 + t293 * t373 + t318 * t328 + (t386 * t432 - t264) * t340 - t466 * t388) * MDP(13) + (-t370 * t446 + t258 * t343 + t294 * t373 - t318 * t327 + (t388 * t432 + t265) * t340 + t466 * t386) * MDP(14) + (-t458 + t264 * t327 - t265 * t328 - t418 * t380 + (qJD(4) * t328 - t257 * t340 - t293 * t370 + t241) * t388 + (qJD(4) * t327 - t258 * t340 + t294 * t370 - t240) * t386) * MDP(15) + (t266 * t373 - t258 * t265 - t257 * t264 - t299 * t318 - g(3) * (t415 + t463) + (-t240 * t386 + t241 * t388) * t370 + t410 * qJD(4) + t418 * (pkin(2) * t391 + pkin(3) * t378 - qJ(4) * t380)) * MDP(16) + (t242 * t357 + t409 * t438) * MDP(17) + (-t242 * t355 - t243 * t357 + t409 * t437 - t438 * t470) * MDP(18) + (t449 - t464) * MDP(19) + (t416 - t450) * MDP(20) - t336 * t343 * MDP(21) + ((-t351 * t393 - t352 * t390) * t309 + t361 * t243 + t247 * t355 - t237 * t343 + t279 * t470 + (t390 * t406 - t393 * t407) * t336 + t437 * t271 + t402 * t379) * MDP(22) + (-(-t351 * t390 + t352 * t393) * t309 + t361 * t242 + t247 * t357 + t238 * t343 + t279 * t409 + (t390 * t407 + t393 * t406) * t336 - t438 * t271 - t402 * t377) * MDP(23) + (-MDP(4) * t391 * t394 + MDP(5) * t436) * t397; (-t343 ^ 2 - t339) * MDP(11) + (-pkin(2) * t429 + t312 * t343 + t428 - t451 - t465) * MDP(12) + (t328 * t343 + t446) * MDP(13) + (-t327 * t343 - t339 * t388 - t447) * MDP(14) + (-t293 * t386 - t294 * t388) * MDP(15) + (-t299 * t343 - t413 - t465) * MDP(16) + (t416 + t450) * MDP(22) + (t449 + t464) * MDP(23) + (t313 * MDP(12) + (t327 * t386 + t328 * t388) * MDP(15) + t410 * MDP(16) - MDP(13) * t444) * t340; (t327 * t340 + t293) * MDP(13) + (t328 * t340 + t294) * MDP(14) + (-t327 ^ 2 - t328 ^ 2) * MDP(15) + (t257 * t327 - t258 * t328 + t466) * MDP(16) + (t243 - t468) * MDP(22) + (t242 + t471) * MDP(23); t409 * t470 * MDP(17) + (t409 ^ 2 - t470 ^ 2) * MDP(18) + (t426 - t471) * MDP(19) + (-t419 - t468) * MDP(20) + t309 * MDP(21) + (-g(1) * t331 + g(2) * t329 + t238 * t336 + t271 * t409 + t377 * t458 + t420) * MDP(22) + (g(1) * t332 - g(2) * t330 + t237 * t336 - t271 * t470 + t379 * t458 - t414) * MDP(23) + (-MDP(19) * t445 + MDP(20) * t409 - MDP(22) * t238 - MDP(23) * t237) * qJD(5);];
tau = t1;
