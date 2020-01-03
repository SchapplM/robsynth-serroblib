% Calculate vector of inverse dynamics joint torques for
% S5RPRPR15
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR15_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR15_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR15_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:32
% EndTime: 2019-12-31 18:37:38
% DurationCPUTime: 3.94s
% Computational Cost: add. (1794->376), mult. (3525->505), div. (0->0), fcn. (2296->10), ass. (0->166)
t326 = sin(pkin(8));
t327 = cos(pkin(8));
t328 = sin(qJ(5));
t331 = cos(qJ(5));
t426 = -t326 * t328 + t331 * t327;
t425 = t426 * qJD(5);
t329 = sin(qJ(3));
t395 = qJD(1) * t329;
t311 = qJD(5) + t395;
t332 = cos(qJ(3));
t394 = qJD(1) * t332;
t376 = t326 * t394;
t389 = t327 * qJD(3);
t289 = t376 - t389;
t375 = t327 * t394;
t393 = qJD(3) * t326;
t291 = t375 + t393;
t352 = t289 * t328 - t291 * t331;
t429 = t352 * t311;
t330 = sin(qJ(1));
t333 = cos(qJ(1));
t427 = g(1) * t330 - g(2) * t333;
t428 = t427 * t327;
t295 = t326 * t331 + t327 * t328;
t343 = t295 * qJD(5);
t422 = g(3) * t329;
t424 = t332 * t427 - t422;
t334 = -pkin(1) - pkin(6);
t421 = g(3) * t332;
t420 = pkin(7) + qJ(4);
t419 = (pkin(1) * qJDD(1));
t336 = qJD(1) ^ 2;
t418 = qJ(2) * t336;
t417 = qJ(4) * t332;
t276 = t331 * t289;
t412 = t291 * t328;
t242 = t276 + t412;
t416 = t242 * t311;
t307 = t334 * qJDD(1) + qJDD(2);
t308 = t334 * qJD(1) + qJD(2);
t392 = qJD(3) * t329;
t359 = -qJDD(3) * pkin(3) + t308 * t392 + qJDD(4);
t252 = -t307 * t332 + t359;
t415 = t252 * t326;
t414 = t252 * t327;
t413 = t252 * t332;
t387 = qJD(1) * qJD(3);
t371 = t332 * t387;
t383 = qJDD(1) * t329;
t344 = t371 + t383;
t293 = qJDD(5) + t344;
t411 = t293 * t426;
t410 = t293 * t295;
t409 = t308 * t332;
t407 = t326 * t332;
t406 = t327 * t332;
t298 = t329 * t308;
t405 = t329 * t330;
t404 = t329 * t333;
t403 = t329 * t334;
t361 = pkin(3) * t332 + qJ(4) * t329;
t274 = t361 * qJD(3) - qJD(4) * t332 + qJD(2);
t360 = pkin(3) * t329 - t417;
t299 = qJ(2) + t360;
t236 = t274 * qJD(1) + t299 * qJDD(1);
t247 = qJDD(3) * qJ(4) + t307 * t329 + (qJD(4) + t409) * qJD(3);
t225 = t326 * t236 + t327 * t247;
t283 = t299 * qJD(1);
t284 = qJD(3) * qJ(4) + t298;
t239 = t326 * t283 + t327 * t284;
t346 = t426 * t329;
t401 = qJD(1) * t346 + t425;
t342 = t295 * qJD(1);
t400 = t329 * t342 + t343;
t297 = t361 * qJD(1);
t251 = t326 * t297 + t308 * t406;
t391 = qJD(3) * t332;
t373 = t334 * t391;
t246 = t326 * t274 + t327 * t373;
t257 = t326 * t299 + t327 * t403;
t382 = qJDD(1) * t332;
t399 = t326 * qJDD(3) + t327 * t382;
t398 = t333 * pkin(1) + t330 * qJ(2);
t324 = t329 ^ 2;
t325 = t332 ^ 2;
t397 = t324 - t325;
t335 = qJD(3) ^ 2;
t396 = -t335 - t336;
t277 = -qJD(3) * pkin(3) + qJD(4) - t409;
t388 = -qJD(4) + t277;
t386 = qJDD(1) * qJ(2);
t385 = qJDD(1) * t326;
t384 = qJDD(1) * t327;
t381 = qJDD(3) * t329;
t380 = pkin(7) * t327 * t329;
t379 = 0.2e1 * qJD(1) * qJD(2);
t365 = -qJDD(3) * t327 + t326 * t382;
t372 = t329 * t387;
t258 = t326 * t372 - t365;
t259 = t327 * t372 - t399;
t378 = -qJD(5) * t276 + t328 * t258 - t331 * t259;
t377 = t326 * t395;
t374 = t326 * t392;
t370 = -t326 * t334 + pkin(4);
t369 = pkin(4) * t326 - t334;
t224 = t327 * t236 - t247 * t326;
t220 = t344 * pkin(4) + pkin(7) * t259 + t224;
t221 = pkin(7) * t258 + t225;
t367 = t331 * t220 - t221 * t328;
t366 = -t331 * t258 - t328 * t259;
t238 = t327 * t283 - t284 * t326;
t364 = g(1) * t333 + g(2) * t330;
t250 = t327 * t297 - t308 * t407;
t362 = qJDD(2) - t427;
t358 = t220 * t328 + t221 * t331;
t357 = -t224 * t326 + t225 * t327;
t226 = pkin(4) * t395 - pkin(7) * t291 + t238;
t228 = -pkin(7) * t289 + t239;
t217 = t226 * t331 - t228 * t328;
t218 = t226 * t328 + t228 * t331;
t356 = t238 * t327 + t239 * t326;
t355 = -t238 * t326 + t239 * t327;
t287 = t327 * t299;
t241 = -pkin(7) * t406 + t370 * t329 + t287;
t249 = -pkin(7) * t407 + t257;
t354 = t241 * t331 - t249 * t328;
t353 = t241 * t328 + t249 * t331;
t351 = -t307 + t427;
t350 = t427 * t326;
t304 = t420 * t327;
t349 = qJD(4) * t326 + qJD(5) * t304 + (pkin(4) * t332 + t380) * qJD(1) + t250;
t303 = t420 * t326;
t348 = pkin(7) * t377 - qJD(4) * t327 + qJD(5) * t303 + t251;
t347 = t311 * t426;
t345 = 0.2e1 * qJ(2) * t387 + qJDD(3) * t334;
t222 = -qJD(5) * t412 + t378;
t341 = t351 + t418;
t339 = -t329 * t427 - t421;
t338 = -t364 + t379 + 0.2e1 * t386;
t223 = -t352 * qJD(5) + t366;
t337 = -t334 * t335 + t338;
t323 = pkin(8) + qJ(5);
t320 = t333 * qJ(2);
t318 = qJDD(3) * t332;
t317 = cos(t323);
t316 = sin(t323);
t313 = -pkin(4) * t327 - pkin(3);
t288 = t369 * t332;
t278 = t369 * t392;
t273 = t426 * t332;
t272 = t295 * t332;
t268 = -t316 * t330 + t317 * t404;
t267 = t316 * t404 + t317 * t330;
t266 = t316 * t333 + t317 * t405;
t265 = -t316 * t405 + t317 * t333;
t264 = -pkin(4) * t377 + t298;
t261 = t327 * t274;
t256 = -t326 * t403 + t287;
t248 = pkin(4) * t289 + t277;
t245 = -t326 * t373 + t261;
t235 = pkin(7) * t374 + t246;
t234 = -t328 * t329 * t389 - t331 * t374 + t332 * t425;
t233 = -qJD(3) * t346 - t332 * t343;
t229 = t261 + (t370 * t332 + t380) * qJD(3);
t227 = -pkin(4) * t258 + t252;
t1 = [qJDD(1) * MDP(1) + t427 * MDP(2) + t364 * MDP(3) + (t362 - (2 * t419)) * MDP(4) + t338 * MDP(5) + (-(qJDD(2) - t419) * pkin(1) - g(1) * (-pkin(1) * t330 + t320) - g(2) * t398 + (t379 + t386) * qJ(2)) * MDP(6) + (qJDD(1) * t325 - 0.2e1 * t329 * t371) * MDP(7) + 0.2e1 * (-t329 * t382 + t397 * t387) * MDP(8) + (-t329 * t335 + t318) * MDP(9) + (-t332 * t335 - t381) * MDP(10) + (t337 * t329 + t345 * t332) * MDP(12) + (-t345 * t329 + t337 * t332) * MDP(13) + (t350 + (t415 + t334 * t258 + (qJD(1) * t256 + t238) * qJD(3)) * t332 + (t245 * qJD(1) + t256 * qJDD(1) + t224 - t364 * t327 + (-t277 * t326 + t289 * t334) * qJD(3)) * t329) * MDP(14) + (t428 + (t414 + t334 * t259 + (-qJD(1) * t257 - t239) * qJD(3)) * t332 + (-t246 * qJD(1) - t257 * qJDD(1) - t225 + t364 * t326 + (-t277 * t327 + t291 * t334) * qJD(3)) * t329) * MDP(15) + (-t245 * t291 - t246 * t289 + t256 * t259 + t257 * t258 + t356 * t392 + (-t224 * t327 - t225 * t326 + t364) * t332) * MDP(16) + (t225 * t257 + t239 * t246 + t224 * t256 + t238 * t245 - g(1) * (pkin(3) * t404 - t333 * t417 + t320) - g(2) * (pkin(6) * t333 + t398) + (t277 * t392 - t413) * t334 + (-g(1) * t334 - g(2) * t360) * t330) * MDP(17) + (t222 * t273 - t233 * t352) * MDP(18) + (-t222 * t272 - t223 * t273 - t233 * t242 + t234 * t352) * MDP(19) + (t222 * t329 + t233 * t311 + t273 * t293 - t352 * t391) * MDP(20) + (-t223 * t329 - t234 * t311 - t242 * t391 - t272 * t293) * MDP(21) + (t293 * t329 + t311 * t391) * MDP(22) + ((t229 * t331 - t235 * t328) * t311 + t354 * t293 + t367 * t329 + t217 * t391 - t278 * t242 + t288 * t223 + t227 * t272 + t248 * t234 - g(1) * t268 - g(2) * t266 + (-t218 * t329 - t353 * t311) * qJD(5)) * MDP(23) + (-(t229 * t328 + t235 * t331) * t311 - t353 * t293 - t358 * t329 - t218 * t391 + t278 * t352 + t288 * t222 + t227 * t273 + t248 * t233 + g(1) * t267 - g(2) * t265 + (-t217 * t329 - t354 * t311) * qJD(5)) * MDP(24); qJDD(1) * MDP(4) - t336 * MDP(5) + (t362 - t418 - t419) * MDP(6) + (t396 * t329 + t318) * MDP(12) + (t396 * t332 - t381) * MDP(13) + (-t324 * t385 + t258 * t332 + (-t327 * t336 + (t289 - 0.2e1 * t376) * qJD(3)) * t329) * MDP(14) + (-t324 * t384 + t259 * t332 + (t326 * t336 + (t291 - 0.2e1 * t375) * qJD(3)) * t329) * MDP(15) + ((qJD(1) * t291 + t258 * t329 - t289 * t391) * t327 + (qJD(1) * t289 - t259 * t329 + t291 * t391) * t326) * MDP(16) + (-t413 + t357 * t329 - t356 * qJD(1) + (t277 * t329 + t355 * t332) * qJD(3) - t427) * MDP(17) + (-qJD(1) * t347 + (-t295 * t311 * qJD(3) - t223) * t332 + (qJD(3) * t242 - t311 * t425 - t410) * t329) * MDP(23) + (t311 * t342 + (-qJD(3) * t347 - t222) * t332 + (-qJD(3) * t352 + t311 * t343 - t411) * t329) * MDP(24); MDP(9) * t382 - MDP(10) * t383 + qJDD(3) * MDP(11) + (-t341 * t332 + t422) * MDP(12) + (t341 * t329 + t421) * MDP(13) + (pkin(3) * t258 - t414 + (-t428 + (-qJ(4) * t393 - t238) * qJD(1)) * t332 + (-qJ(4) * t385 + g(3) * t327 - t289 * t308 + (t388 * t326 - t250) * qJD(1)) * t329) * MDP(14) + (pkin(3) * t259 + t415 + (t350 + (-qJ(4) * t389 + t239) * qJD(1)) * t332 + (-qJ(4) * t384 - g(3) * t326 - t291 * t308 + (t388 * t327 + t251) * qJD(1)) * t329) * MDP(15) + (t250 * t291 + t251 * t289 + (qJ(4) * t258 - qJD(4) * t289 - t238 * t395 + t225) * t327 + (-qJ(4) * t259 + qJD(4) * t291 - t239 * t395 - t224) * t326 + t339) * MDP(16) + (-t277 * t298 - t238 * t250 - t239 * t251 + t355 * qJD(4) + (-t252 - t424) * pkin(3) + (t339 + t357) * qJ(4)) * MDP(17) + (t222 * t295 - t352 * t401) * MDP(18) + (t222 * t426 - t223 * t295 - t401 * t242 + t352 * t400) * MDP(19) + (t401 * t311 + t352 * t394 + t410) * MDP(20) + (t242 * t394 - t400 * t311 + t411) * MDP(21) - t311 * MDP(22) * t394 + ((-t303 * t331 - t304 * t328) * t293 + t313 * t223 - t227 * t426 - t217 * t394 - t264 * t242 + (t348 * t328 - t349 * t331) * t311 + t400 * t248 - t424 * t317) * MDP(23) + (-(-t303 * t328 + t304 * t331) * t293 + t313 * t222 + t227 * t295 + t218 * t394 + t264 * t352 + (t349 * t328 + t348 * t331) * t311 + t401 * t248 + t424 * t316) * MDP(24) + (t332 * t329 * MDP(7) - t397 * MDP(8)) * t336; ((t291 - t393) * t395 + t365) * MDP(14) + ((-t289 - t389) * t395 + t399) * MDP(15) + (-t289 ^ 2 - t291 ^ 2) * MDP(16) + (t238 * t291 + t239 * t289 + t351 * t332 + t359 - t422) * MDP(17) + (t223 - t429) * MDP(23) + (t222 - t416) * MDP(24); -t352 * t242 * MDP(18) + (-t242 ^ 2 + t352 ^ 2) * MDP(19) + (t378 + t416) * MDP(20) + (-t366 - t429) * MDP(21) + t293 * MDP(22) + (-g(1) * t265 - g(2) * t267 + t218 * t311 + t248 * t352 + t316 * t421 + t367) * MDP(23) + (g(1) * t266 - g(2) * t268 + t217 * t311 + t242 * t248 + t317 * t421 - t358) * MDP(24) + (-MDP(20) * t412 + t352 * MDP(21) - t218 * MDP(23) - t217 * MDP(24)) * qJD(5);];
tau = t1;
