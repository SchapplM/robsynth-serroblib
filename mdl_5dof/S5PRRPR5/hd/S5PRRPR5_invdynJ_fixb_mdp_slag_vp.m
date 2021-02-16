% Calculate vector of inverse dynamics joint torques for
% S5PRRPR5
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:53
% EndTime: 2021-01-15 16:05:03
% DurationCPUTime: 3.70s
% Computational Cost: add. (2057->376), mult. (4963->522), div. (0->0), fcn. (3948->14), ass. (0->175)
t332 = qJ(4) + pkin(7);
t328 = sin(pkin(10));
t337 = cos(qJ(3));
t334 = sin(qJ(3));
t432 = cos(pkin(10));
t380 = t432 * t334;
t301 = t328 * t337 + t380;
t294 = t301 * qJD(3);
t379 = t432 * t337;
t401 = qJDD(2) * t334;
t368 = -qJDD(2) * t379 + t328 * t401;
t259 = qJD(2) * t294 + t368;
t255 = qJDD(5) + t259;
t319 = pkin(3) * t337 + pkin(2);
t421 = t328 * t334;
t358 = t379 - t421;
t257 = -pkin(4) * t358 - pkin(8) * t301 - t319;
t325 = qJ(3) + pkin(10);
t321 = cos(t325);
t331 = cos(pkin(5));
t338 = cos(qJ(2));
t433 = cos(pkin(9));
t381 = t433 * t338;
t329 = sin(pkin(9));
t335 = sin(qJ(2));
t419 = t329 * t335;
t289 = -t331 * t381 + t419;
t382 = t433 * t335;
t418 = t329 * t338;
t291 = t331 * t418 + t382;
t370 = g(1) * t291 + g(2) * t289;
t330 = sin(pkin(5));
t415 = t330 * t338;
t352 = g(3) * t415 - t370;
t348 = t352 * t321;
t442 = t257 * t255 - t348;
t295 = t301 * qJD(2);
t399 = t331 * qJDD(1);
t310 = t337 * t399;
t403 = qJD(1) * qJD(2);
t280 = qJDD(2) * pkin(7) + (qJDD(1) * t335 + t338 * t403) * t330;
t410 = qJD(1) * t331;
t344 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t410 + t280;
t409 = qJD(1) * t335;
t392 = t330 * t409;
t375 = t332 * qJD(2) + t392;
t366 = t375 * qJD(3);
t233 = qJDD(3) * pkin(3) - t344 * t334 - t337 * t366 + t310;
t234 = (-t366 + t399) * t334 + t344 * t337;
t224 = t432 * t233 - t328 * t234;
t222 = -qJDD(3) * pkin(4) - t224;
t311 = qJD(2) * t379;
t408 = qJD(2) * t334;
t292 = t328 * t408 - t311;
t286 = qJD(5) + t292;
t438 = pkin(3) * t328;
t316 = pkin(8) + t438;
t288 = t331 * t419 - t381;
t290 = t331 * t382 + t418;
t320 = sin(t325);
t383 = t330 * t433;
t417 = t330 * t335;
t420 = t329 * t330;
t357 = g(1) * (t288 * t320 + t321 * t420) + g(2) * (-t290 * t320 - t321 * t383) + g(3) * (-t320 * t417 + t321 * t331);
t398 = pkin(3) * t408;
t441 = (pkin(4) * t295 + pkin(8) * t292 + qJD(5) * t316 + t398) * t286 + t222 + t357;
t225 = t328 * t233 + t432 * t234;
t223 = qJDD(3) * pkin(8) + t225;
t270 = -t375 * t334 + t337 * t410;
t434 = qJD(3) * pkin(3);
t268 = t270 + t434;
t271 = t334 * t410 + t375 * t337;
t422 = t328 * t271;
t239 = t432 * t268 - t422;
t235 = -qJD(3) * pkin(4) - t239;
t391 = qJD(1) * t415;
t285 = -t319 * qJD(2) + qJD(4) - t391;
t245 = pkin(4) * t292 - pkin(8) * t295 + t285;
t305 = t332 * t337;
t275 = t432 * t305 - t332 * t421;
t297 = t358 * qJD(3);
t371 = g(1) * t288 - g(2) * t290;
t353 = -g(3) * t417 + t371;
t384 = qJD(3) * t332;
t287 = qJD(4) * t337 - t334 * t384;
t355 = -qJD(4) * t334 - t337 * t384;
t412 = t432 * t287 + t328 * t355 - t358 * t391;
t440 = (qJD(5) * t245 + t223) * t358 + t222 * t301 + t235 * t297 + (-qJD(5) * t257 - t412) * t286 - t275 * t255 + t353;
t388 = t335 * t403;
t308 = t330 * t388;
t339 = qJD(3) ^ 2;
t385 = qJDD(1) * t415;
t439 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t339 + (-g(3) * t338 + t388) * t330 - t308 + t370 + t385;
t437 = pkin(3) * t334;
t298 = t331 * t337 - t334 * t417;
t436 = g(3) * t298;
t435 = qJD(2) * pkin(2);
t402 = qJD(2) * qJD(3);
t387 = t334 * t402;
t346 = t301 * qJDD(2) - t328 * t387;
t260 = qJD(3) * t311 + t346;
t333 = sin(qJ(5));
t336 = cos(qJ(5));
t404 = t336 * qJD(3);
t394 = qJD(5) * t404 + t333 * qJDD(3) + t336 * t260;
t405 = qJD(5) * t333;
t237 = -t295 * t405 + t394;
t430 = t237 * t333;
t429 = t255 * t333;
t423 = t295 * t333;
t276 = -t404 + t423;
t427 = t276 * t286;
t426 = t276 * t295;
t278 = qJD(3) * t333 + t295 * t336;
t425 = t278 * t286;
t424 = t278 * t295;
t416 = t330 * t337;
t252 = t336 * t255;
t414 = qJDD(1) - g(3);
t413 = t287 * t328 - t301 * t391 - t432 * t355;
t266 = t432 * t271;
t240 = t328 * t268 + t266;
t326 = t334 ^ 2;
t411 = -t337 ^ 2 + t326;
t407 = qJD(2) * t335;
t406 = qJD(5) * t301;
t400 = qJDD(2) * t337;
t397 = t334 * t434;
t396 = t333 * t415;
t395 = t336 * t415;
t393 = t432 * pkin(3);
t390 = t330 * t407;
t389 = qJD(2) * t415;
t386 = t337 * t402;
t378 = t330 * t414;
t377 = -t336 * qJDD(3) + t260 * t333;
t376 = t286 * t336;
t372 = t334 * t389;
t369 = pkin(4) * t294 - pkin(8) * t297 - t392 + t397;
t236 = qJD(3) * pkin(8) + t240;
t227 = t236 * t336 + t245 * t333;
t367 = t236 * t333 - t245 * t336;
t340 = qJD(2) ^ 2;
t365 = qJDD(2) * t338 - t335 * t340;
t364 = t252 + (-t292 * t333 - t405) * t286;
t363 = -g(1) * t329 + t433 * g(2);
t299 = t331 * t334 + t335 * t416;
t251 = t328 * t298 + t432 * t299;
t361 = -t251 * t333 - t395;
t360 = -t251 * t336 + t396;
t264 = -t288 * t321 + t320 * t420;
t359 = t297 * t336 - t301 * t405;
t356 = t365 * t330;
t354 = t299 * qJD(3);
t304 = -t391 - t435;
t350 = -qJD(2) * t304 - t280 - t371;
t349 = pkin(3) * t387 - t319 * qJDD(2) + qJDD(4) + t308;
t244 = t432 * t270 - t422;
t347 = -t316 * t255 + (t235 + t244) * t286;
t343 = -pkin(7) * qJDD(3) + (t304 + t391 - t435) * qJD(3);
t258 = t349 - t385;
t342 = -t354 - t372;
t317 = -t393 - pkin(4);
t283 = t320 * t331 + t321 * t417;
t274 = t305 * t328 + t332 * t380;
t269 = t298 * qJD(3) + t337 * t389;
t262 = t290 * t321 - t320 * t383;
t250 = -t432 * t298 + t299 * t328;
t243 = t432 * t269 + t328 * t342;
t242 = t270 * t328 + t266;
t241 = t269 * t328 - t432 * t342;
t238 = t278 * qJD(5) + t377;
t229 = pkin(4) * t259 - pkin(8) * t260 + t258;
t228 = t336 * t229;
t1 = [t414 * MDP(1) + MDP(3) * t356 + (-qJDD(2) * t335 - t338 * t340) * t330 * MDP(4) + (t298 * qJDD(3) + t337 * t356 + (-t354 - 0.2e1 * t372) * qJD(3)) * MDP(10) + (-qJD(3) * t269 - qJDD(3) * t299 + (-t365 * t334 - t338 * t386) * t330) * MDP(11) + (-qJD(3) * t241 - qJDD(3) * t250 + (-t259 * t338 + t292 * t407) * t330) * MDP(12) + (-qJD(3) * t243 - qJDD(3) * t251 + (-t260 * t338 + t295 * t407) * t330) * MDP(13) + (t241 * t295 - t243 * t292 + t250 * t260 - t251 * t259) * MDP(14) + (-t224 * t250 + t225 * t251 - t239 * t241 + t240 * t243 - g(3) + (-t258 * t338 + t285 * t407) * t330) * MDP(15) + ((t360 * qJD(5) - t243 * t333 + t336 * t390) * t286 + t361 * t255 + t241 * t276 + t250 * t238) * MDP(21) + (-(t361 * qJD(5) + t243 * t336 + t333 * t390) * t286 + t360 * t255 + t241 * t278 + t250 * t237) * MDP(22); qJDD(2) * MDP(2) + (t414 * t415 + t370) * MDP(3) + (-t335 * t378 - t371) * MDP(4) + (qJDD(2) * t326 + 0.2e1 * t334 * t386) * MDP(5) + 0.2e1 * (t334 * t400 - t411 * t402) * MDP(6) + (qJDD(3) * t334 + t337 * t339) * MDP(7) + (qJDD(3) * t337 - t334 * t339) * MDP(8) + (t343 * t334 + t439 * t337) * MDP(10) + (-t439 * t334 + t343 * t337) * MDP(11) + (-t292 * t392 - qJDD(3) * t274 - t258 * t358 - t259 * t319 + t285 * t294 - t348 + (t292 * t437 - t413) * qJD(3)) * MDP(12) + (-t295 * t392 - qJDD(3) * t275 + t258 * t301 - t260 * t319 + t285 * t297 + t352 * t320 + (t295 * t437 - t412) * qJD(3)) * MDP(13) + (-t224 * t301 + t225 * t358 - t239 * t297 - t240 * t294 - t259 * t275 + t260 * t274 - t412 * t292 + t413 * t295 + t353) * MDP(14) + (t225 * t275 - t224 * t274 - t258 * t319 + t285 * t397 - g(1) * (-t288 * t332 - t291 * t319) - g(2) * (-t289 * t319 + t290 * t332) + t412 * t240 - t413 * t239 + (-t285 * t409 - g(3) * (t319 * t338 + t332 * t335)) * t330) * MDP(15) + (t237 * t301 * t336 + t359 * t278) * MDP(16) + ((-t276 * t336 - t278 * t333) * t297 + (-t430 - t238 * t336 + (t276 * t333 - t278 * t336) * qJD(5)) * t301) * MDP(17) + (-t237 * t358 + t301 * t252 + t278 * t294 + t359 * t286) * MDP(18) + (-t301 * t429 + t238 * t358 - t276 * t294 + (-t297 * t333 - t336 * t406) * t286) * MDP(19) + (-t255 * t358 + t286 * t294) * MDP(20) + (-t367 * t294 - t228 * t358 + t274 * t238 + t413 * t276 + (t369 * t286 + (t235 * t301 + t236 * t358 - t275 * t286) * qJD(5) + t442) * t336 + t440 * t333) * MDP(21) + (-t227 * t294 + t274 * t237 + t413 * t278 + ((-qJD(5) * t236 + t229) * t358 - t235 * t406 + (qJD(5) * t275 - t369) * t286 - t442) * t333 + t440 * t336) * MDP(22); MDP(7) * t401 + MDP(8) * t400 + qJDD(3) * MDP(9) + (t350 * t334 + t363 * t416 + t310 - t436) * MDP(10) + (g(3) * t299 + (-t363 * t330 - t399) * t334 + t350 * t337) * MDP(11) + (t242 * qJD(3) - t285 * t295 + (t432 * qJDD(3) - t292 * t408) * pkin(3) - t357 + t224) * MDP(12) + (t244 * qJD(3) + t285 * t292 + g(1) * t264 + g(2) * t262 + g(3) * t283 + (-qJDD(3) * t328 - t295 * t408) * pkin(3) - t225) * MDP(13) + ((t240 - t242) * t295 + (-t239 + t244) * t292 + (-t259 * t328 - t432 * t260) * pkin(3)) * MDP(14) + (t224 * t393 + t225 * t438 + t239 * t242 - t240 * t244 - t285 * t398 + (-g(1) * (t288 * t334 + t329 * t416) - g(2) * (-t290 * t334 - t337 * t383) - t436) * pkin(3)) * MDP(15) + (t278 * t376 + t430) * MDP(16) + ((t237 - t427) * t336 + (-t238 - t425) * t333) * MDP(17) + (t286 * t376 - t424 + t429) * MDP(18) + (t364 + t426) * MDP(19) - t286 * t295 * MDP(20) + (t317 * t238 - t242 * t276 + t295 * t367 + t347 * t333 - t441 * t336) * MDP(21) + (t227 * t295 + t317 * t237 - t242 * t278 + t441 * t333 + t347 * t336) * MDP(22) + (-t334 * t337 * MDP(5) + t411 * MDP(6)) * t340; (0.2e1 * t295 * qJD(3) + t368) * MDP(12) + ((t311 - t292) * qJD(3) + t346) * MDP(13) + (-t292 ^ 2 - t295 ^ 2) * MDP(14) + (t239 * t295 + t240 * t292 - t338 * t378 + t349 - t370) * MDP(15) + (t364 - t426) * MDP(21) + (-t286 ^ 2 * t336 - t424 - t429) * MDP(22); t278 * t276 * MDP(16) + (-t276 ^ 2 + t278 ^ 2) * MDP(17) + (t394 + t427) * MDP(18) + (-t377 + t425) * MDP(19) + t255 * MDP(20) + (-t333 * t223 + t228 + t227 * t286 - t235 * t278 - g(1) * (-t264 * t333 + t291 * t336) - g(2) * (-t262 * t333 + t289 * t336) - g(3) * (-t283 * t333 - t395)) * MDP(21) + (-t336 * t223 - t333 * t229 - t367 * t286 + t235 * t276 - g(1) * (-t264 * t336 - t291 * t333) - g(2) * (-t262 * t336 - t289 * t333) - g(3) * (-t283 * t336 + t396)) * MDP(22) + (-MDP(18) * t423 - t278 * MDP(19) - t227 * MDP(21) + t367 * MDP(22)) * qJD(5);];
tau = t1;
