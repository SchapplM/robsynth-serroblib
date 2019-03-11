% Calculate vector of inverse dynamics joint torques for
% S6PRPPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPPRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRPPRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:11
% EndTime: 2019-03-08 19:23:15
% DurationCPUTime: 3.00s
% Computational Cost: add. (1621->360), mult. (3329->495), div. (0->0), fcn. (2670->12), ass. (0->169)
t329 = sin(pkin(11));
t331 = sin(pkin(6));
t332 = cos(pkin(11));
t338 = sin(qJ(2));
t341 = cos(qJ(2));
t277 = (t329 * t338 + t332 * t341) * t331;
t330 = sin(pkin(10));
t333 = cos(pkin(10));
t334 = cos(pkin(6));
t426 = t334 * t341;
t284 = t330 * t338 - t333 * t426;
t427 = t334 * t338;
t285 = t330 * t341 + t333 * t427;
t286 = t330 * t426 + t333 * t338;
t287 = -t330 * t427 + t333 * t341;
t442 = g(1) * (-t286 * t332 + t287 * t329) + g(2) * (-t284 * t332 + t285 * t329);
t357 = g(3) * t277 + t442;
t418 = qJD(1) * t331;
t391 = t341 * t418;
t395 = t338 * t418;
t267 = t329 * t391 - t332 * t395;
t379 = qJD(3) * t329 - t267;
t448 = t379 * qJD(2) - t357;
t337 = sin(qJ(5));
t409 = qJD(6) * t337;
t447 = qJD(2) * t409 + qJDD(5);
t342 = -pkin(2) - pkin(3);
t431 = t331 * t338;
t394 = qJD(2) * t431;
t308 = qJD(1) * t394;
t403 = qJDD(1) * t331;
t311 = t341 * t403;
t387 = qJDD(3) + t308 - t311;
t264 = qJDD(2) * t342 + t387;
t310 = t338 * t403;
t402 = qJDD(2) * qJ(3);
t265 = t402 + t310 + (qJD(3) + t391) * qJD(2);
t230 = t264 * t332 - t329 * t265;
t228 = qJDD(2) * pkin(4) - t230;
t299 = -t329 * qJ(3) + t332 * t342;
t295 = pkin(4) - t299;
t300 = t332 * qJ(3) + t329 * t342;
t296 = -pkin(8) + t300;
t343 = qJD(5) ^ 2;
t446 = qJDD(2) * t295 - t296 * t343 + t228 + t448;
t231 = t329 * t264 + t332 * t265;
t229 = -qJDD(2) * pkin(8) + t231;
t370 = qJD(3) - t391;
t288 = qJD(2) * t342 + t370;
t298 = qJD(2) * qJ(3) + t395;
t258 = t329 * t288 + t332 * t298;
t254 = -qJD(2) * pkin(8) + t258;
t314 = -qJD(1) * t334 + qJD(4);
t340 = cos(qJ(5));
t243 = t254 * t340 + t314 * t337;
t312 = -qJDD(1) * t334 + qJDD(4);
t434 = t312 * t340;
t223 = -qJDD(5) * pkin(5) + qJD(5) * t243 + t229 * t337 - t434;
t317 = qJD(2) * t340 + qJD(6);
t248 = t284 * t329 + t285 * t332;
t252 = t286 * t329 + t287 * t332;
t429 = t331 * t341;
t397 = t329 * t429;
t278 = t332 * t431 - t397;
t259 = t278 * t337 + t334 * t340;
t430 = t331 * t340;
t359 = g(1) * (-t252 * t337 - t330 * t430) + g(2) * (-t248 * t337 + t333 * t430) - g(3) * t259;
t374 = -pkin(5) * t337 + pkin(9) * t340;
t445 = t317 * (pkin(9) * qJD(6) + t374 * qJD(2)) + t223 + t359;
t375 = pkin(5) * t340 + pkin(9) * t337;
t271 = t295 + t375;
t405 = qJD(2) * qJD(5);
t390 = t337 * t405;
t398 = t340 * qJDD(2);
t360 = -t390 + t398;
t291 = -qJDD(6) - t360;
t444 = (-qJD(5) * t374 - t379) * t317 + t271 * t291;
t237 = qJD(5) * pkin(9) + t243;
t436 = t296 * t317;
t443 = -qJD(6) * (t237 + t436) - t357;
t441 = qJDD(2) * pkin(2);
t339 = cos(qJ(6));
t389 = t340 * t405;
t401 = qJDD(2) * t337;
t361 = t389 + t401;
t336 = sin(qJ(6));
t404 = qJD(5) * qJD(6);
t396 = t447 * t336 + t339 * t404;
t255 = -t339 * t361 + t396;
t440 = t255 * t336;
t406 = t339 * qJD(5);
t417 = qJD(2) * t337;
t293 = t336 * t417 + t406;
t438 = t293 * t317;
t407 = t336 * qJD(5);
t294 = t339 * t417 - t407;
t437 = t294 * t317;
t435 = t312 * t337;
t433 = t317 * t339;
t432 = t331 * t337;
t428 = t332 * t337;
t425 = t336 * t340;
t424 = t339 * t340;
t423 = qJDD(1) - g(3);
t421 = pkin(2) * t429 + qJ(3) * t431;
t326 = t337 ^ 2;
t420 = -t340 ^ 2 + t326;
t344 = qJD(2) ^ 2;
t419 = t343 + t344;
t415 = qJD(5) * t293;
t414 = qJD(5) * t294;
t413 = qJD(5) * t296;
t412 = qJD(5) * t337;
t411 = qJD(5) * t340;
t410 = qJD(6) * t336;
t408 = qJD(6) * t339;
t400 = qJDD(5) * t337;
t399 = qJDD(5) * t340;
t393 = t317 * t407;
t392 = t317 * t406;
t385 = -t284 * pkin(2) + qJ(3) * t285;
t384 = -t286 * pkin(2) + qJ(3) * t287;
t257 = t288 * t332 - t329 * t298;
t253 = qJD(2) * pkin(4) - t257;
t380 = qJD(2) * t253 - t229;
t270 = qJD(1) * t277;
t378 = t332 * qJD(3) - t270;
t244 = qJD(2) * t375 + t253;
t369 = t254 * t337 - t314 * t340;
t377 = -qJDD(5) * pkin(9) + qJD(5) * t369 - qJD(6) * t244 - t229 * t340 - t435;
t376 = 0.2e1 * t389;
t373 = g(1) * (t252 * t340 - t330 * t432) + g(2) * (t248 * t340 + t333 * t432);
t372 = t255 + t392;
t318 = t336 * t404;
t256 = t336 * t401 + qJDD(5) * t339 - t318 + (t337 * t408 + t340 * t407) * qJD(2);
t371 = -t256 + t393;
t368 = t257 * t329 - t258 * t332;
t260 = t278 * t340 - t334 * t337;
t235 = t260 * t339 + t277 * t336;
t234 = -t260 * t336 + t277 * t339;
t367 = (-qJD(2) * pkin(2) + t370) * t338 + t298 * t341;
t365 = -g(1) * t286 - g(2) * t284 + g(3) * t429;
t364 = -t336 * t291 + t317 * t408;
t363 = t339 * t291 + t317 * t410;
t362 = t311 - t365;
t358 = -g(1) * t252 - g(2) * t248 - g(3) * t278;
t356 = -qJD(6) * t237 - t442;
t354 = -qJDD(3) + t362;
t353 = t373 + t377;
t352 = g(1) * t287 + g(2) * t285 + g(3) * t431 - t310;
t351 = -t364 - t415;
t350 = t363 - t414;
t349 = -qJD(6) * t271 * t317 - t358;
t236 = -qJD(5) * pkin(5) + t369;
t348 = pkin(9) * t291 + (t236 - t369) * t317;
t347 = -qJDD(5) * t296 + (-qJD(2) * t295 - t253 - t378) * qJD(5);
t345 = -t236 * qJD(5) + t296 * t291 - t317 * t378 + t377;
t304 = -t337 * t343 + t399;
t303 = -t340 * t343 - t400;
t269 = qJD(2) * t277;
t268 = qJD(2) * t397 - t332 * t394;
t266 = t387 - t441;
t233 = -qJD(5) * t259 + t269 * t340;
t232 = qJD(5) * t260 + t269 * t337;
t227 = pkin(5) * t360 + pkin(9) * t361 + t228;
t226 = t339 * t227;
t225 = t339 * t237 + t336 * t244;
t224 = -t336 * t237 + t339 * t244;
t1 = [t423 * MDP(1) + (qJDD(1) * t334 ^ 2 + t265 * t431 - t266 * t429 - g(3)) * MDP(7) + (-t230 * t277 + t231 * t278 - t257 * t268 + t258 * t269 - t312 * t334 - g(3)) * MDP(10) + (-qJD(5) * t232 - qJDD(5) * t259) * MDP(16) + (-qJD(5) * t233 - qJDD(5) * t260) * MDP(17) + ((-qJD(6) * t235 - t233 * t336 + t268 * t339) * t317 - t234 * t291 - t232 * t293 - t259 * t256) * MDP(23) + (-(qJD(6) * t234 + t233 * t339 + t268 * t336) * t317 + t235 * t291 - t232 * t294 + t259 * t255) * MDP(24) + (t278 * MDP(9) + (MDP(16) * t340 - MDP(17) * t337 + MDP(8)) * t277) * qJDD(2) + (t268 * MDP(8) + t269 * MDP(9) + (t268 * t340 - t277 * t412) * MDP(16) + (-t268 * t337 - t277 * t411) * MDP(17)) * qJD(2) + ((-MDP(4) + MDP(6)) * (qJDD(2) * t338 + t341 * t344) + (MDP(3) + MDP(5)) * (qJDD(2) * t341 - t338 * t344) + t367 * MDP(7) * qJD(2)) * t331; qJDD(2) * MDP(2) + t362 * MDP(3) + t352 * MDP(4) + (t354 + 0.2e1 * t441) * MDP(5) + (0.2e1 * qJD(2) * qJD(3) - t352 + 0.2e1 * t402) * MDP(6) + (-t266 * pkin(2) - g(1) * t384 - g(2) * t385 - g(3) * t421 + t265 * qJ(3) + t298 * qJD(3) - t367 * t418) * MDP(7) + (-qJDD(2) * t299 - t230 + t448) * MDP(8) + (qJD(2) * t378 + qJDD(2) * t300 + t231 + t358) * MDP(9) + (t231 * t300 + t230 * t299 - t258 * t270 + t257 * t267 - g(1) * (-pkin(3) * t286 + t384) - g(2) * (-pkin(3) * t284 + t385) - g(3) * (pkin(3) * t429 + t421) - t368 * qJD(3)) * MDP(10) + (qJDD(2) * t326 + t337 * t376) * MDP(11) + 0.2e1 * (t337 * t398 - t405 * t420) * MDP(12) + t303 * MDP(13) - t304 * MDP(14) + (t347 * t337 + t446 * t340) * MDP(16) + (-t446 * t337 + t347 * t340) * MDP(17) + (-t255 * t337 * t339 + (-t336 * t409 + t340 * t406) * t294) * MDP(18) + ((-t293 * t339 - t294 * t336) * t411 + (t440 - t256 * t339 + (t293 * t336 - t294 * t339) * qJD(6)) * t337) * MDP(19) + ((t255 - t392) * t340 + (t363 + t414) * t337) * MDP(20) + ((t256 + t393) * t340 + (t364 - t415) * t337) * MDP(21) + (-t291 * t340 - t317 * t412) * MDP(22) + (-t444 * t339 + t349 * t336 + (-t293 * t413 + t345 * t336 + t443 * t339 + t226) * t340 + (-t236 * t408 - t223 * t336 - t296 * t256 - t378 * t293 + (t336 * t436 - t224) * qJD(5)) * t337) * MDP(23) + (t444 * t336 + t349 * t339 + (-t294 * t413 + t345 * t339 + (-t227 - t443) * t336) * t340 + (t236 * t410 - t223 * t339 + t296 * t255 - t378 * t294 + (t296 * t433 + t225) * qJD(5)) * t337) * MDP(24); -qJDD(2) * MDP(5) - t344 * MDP(6) + (-qJD(2) * t298 + t308 - t354 - t441) * MDP(7) + (-qJDD(2) * t332 - t329 * t344) * MDP(8) + (qJDD(2) * t329 - t332 * t344) * MDP(9) + (qJD(2) * t368 + t230 * t332 + t231 * t329 + t365) * MDP(10) + ((0.2e1 * t390 - t398) * t332 + (-t340 * t419 - t400) * t329) * MDP(16) + ((t376 + t401) * t332 + (t337 * t419 - t399) * t329) * MDP(17) + (t363 * t332 + (t337 * t371 + t340 * t351) * t329 + (-(t329 * t339 - t332 * t425) * t317 + t293 * t428) * qJD(2)) * MDP(23) + (t364 * t332 + (t337 * t372 + t340 * t350) * t329 + ((t329 * t336 + t332 * t424) * t317 + t294 * t428) * qJD(2)) * MDP(24); (qJDD(4) - t423 * t334 + (g(1) * t330 - g(2) * t333) * t331) * MDP(10) + t304 * MDP(16) + t303 * MDP(17) + (-MDP(23) * t371 - MDP(24) * t372) * t340 + (MDP(23) * t351 + MDP(24) * t350) * t337; -MDP(13) * t401 - MDP(14) * t398 + qJDD(5) * MDP(15) + (t337 * t380 - t359 + t434) * MDP(16) + (g(3) * t260 + t380 * t340 + t373 - t435) * MDP(17) + (-t294 * t433 + t440) * MDP(18) + ((t255 + t438) * t339 + (t256 + t437) * t336) * MDP(19) + ((-t294 * t337 + t317 * t424) * qJD(2) + t364) * MDP(20) + ((t293 * t337 - t317 * t425) * qJD(2) - t363) * MDP(21) + t317 * MDP(22) * t417 + (pkin(5) * t256 + t224 * t417 + t243 * t293 + t348 * t336 - t445 * t339) * MDP(23) + (-pkin(5) * t255 - t225 * t417 + t243 * t294 + t445 * t336 + t348 * t339) * MDP(24) + (-MDP(11) * t337 * t340 + MDP(12) * t420) * t344; t294 * t293 * MDP(18) + (-t293 ^ 2 + t294 ^ 2) * MDP(19) + (t396 - t438) * MDP(20) + (-t318 - t437) * MDP(21) - t291 * MDP(22) + (-g(3) * t234 + t225 * t317 + t236 * t294 + t226) * MDP(23) + (g(3) * t235 + t224 * t317 - t236 * t293) * MDP(24) + (-t361 * MDP(20) + t447 * MDP(21) + t356 * MDP(23) + t353 * MDP(24)) * t339 + (t361 * MDP(21) + t353 * MDP(23) + (-t227 - t356) * MDP(24)) * t336;];
tau  = t1;
