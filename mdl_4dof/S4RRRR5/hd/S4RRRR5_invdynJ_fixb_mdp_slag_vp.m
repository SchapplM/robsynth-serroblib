% Calculate vector of inverse dynamics joint torques for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:21
% EndTime: 2019-12-31 17:28:26
% DurationCPUTime: 3.61s
% Computational Cost: add. (1676->356), mult. (3793->492), div. (0->0), fcn. (2620->10), ass. (0->158)
t330 = sin(qJ(2));
t379 = qJDD(1) * t330;
t316 = pkin(5) * t379;
t334 = cos(qJ(2));
t380 = qJD(1) * qJD(2);
t367 = t334 * t380;
t268 = -qJDD(2) * pkin(2) + pkin(5) * t367 + t316;
t392 = qJD(1) * t334;
t311 = -qJD(3) + t392;
t331 = sin(qJ(1));
t335 = cos(qJ(1));
t358 = g(1) * t335 + g(2) * t331;
t415 = g(3) * t334;
t341 = t330 * t358 - t415;
t431 = qJD(3) * pkin(6) * t311 - t268 + t341;
t329 = sin(qJ(3));
t393 = qJD(1) * t330;
t375 = t329 * t393;
t333 = cos(qJ(3));
t381 = t333 * qJD(2);
t286 = t375 - t381;
t332 = cos(qJ(4));
t389 = qJD(2) * t329;
t288 = t333 * t393 + t389;
t328 = sin(qJ(4));
t409 = t288 * t328;
t238 = t332 * t286 + t409;
t306 = -qJD(4) + t311;
t430 = t238 * t306;
t352 = t286 * t328 - t332 * t288;
t429 = t306 * t352;
t385 = qJD(3) * t330;
t428 = -qJD(1) * t385 + qJDD(2);
t296 = -pkin(2) * t334 - pkin(6) * t330 - pkin(1);
t280 = t296 * qJD(1);
t318 = pkin(5) * t392;
t300 = qJD(2) * pkin(6) + t318;
t247 = t280 * t329 + t300 * t333;
t230 = -pkin(7) * t286 + t247;
t383 = qJD(4) * t328;
t228 = t230 * t383;
t299 = -qJD(2) * pkin(2) + pkin(5) * t393;
t252 = pkin(3) * t286 + t299;
t327 = qJ(3) + qJ(4);
t323 = sin(t327);
t324 = cos(t327);
t404 = t331 * t334;
t259 = t323 * t335 - t324 * t404;
t402 = t334 * t335;
t261 = t323 * t331 + t324 * t402;
t416 = g(3) * t330;
t427 = g(1) * t261 - g(2) * t259 + t238 * t252 + t324 * t416 + t228;
t258 = t323 * t404 + t324 * t335;
t260 = -t323 * t402 + t324 * t331;
t235 = qJD(3) * t381 + (t367 + t379) * t333 + t428 * t329;
t359 = pkin(2) * t330 - pkin(6) * t334;
t294 = t359 * qJD(2);
t250 = qJD(1) * t294 + qJDD(1) * t296;
t243 = t333 * t250;
t321 = t334 * qJDD(1);
t423 = -t330 * t380 + t321;
t267 = t423 * pkin(5) + qJDD(2) * pkin(6);
t283 = qJDD(3) - t423;
t216 = pkin(3) * t283 - pkin(7) * t235 - qJD(3) * t247 - t329 * t267 + t243;
t236 = (qJD(2) * (qJD(3) + t392) + t379) * t329 - t428 * t333;
t384 = qJD(3) * t333;
t386 = qJD(3) * t329;
t344 = t329 * t250 + t333 * t267 + t280 * t384 - t300 * t386;
t217 = -pkin(7) * t236 + t344;
t365 = t332 * t216 - t328 * t217;
t426 = -g(1) * t260 + g(2) * t258 + t252 * t352 + t323 * t416 + t365;
t279 = qJDD(4) + t283;
t425 = t279 * MDP(22) + (-t238 ^ 2 + t352 ^ 2) * MDP(19) - t238 * MDP(18) * t352;
t290 = t328 * t333 + t329 * t332;
t263 = t290 * t330;
t422 = -t329 * t385 + t334 * t381;
t421 = qJD(3) + qJD(4);
t364 = t235 * t328 + t332 * t236;
t219 = -qJD(4) * t352 + t364;
t419 = pkin(6) + pkin(7);
t418 = pkin(3) * t329;
t417 = pkin(5) * t329;
t246 = t333 * t280 - t300 * t329;
t229 = -pkin(7) * t288 + t246;
t227 = -pkin(3) * t311 + t229;
t414 = t227 * t332;
t413 = t230 * t332;
t412 = t235 * t329;
t411 = t286 * t311;
t410 = t288 * t311;
t408 = t288 * t333;
t407 = t329 * t330;
t406 = t329 * t334;
t405 = t330 * t333;
t403 = t333 * t334;
t289 = t328 * t329 - t332 * t333;
t345 = t289 * t334;
t401 = qJD(1) * t345 - t421 * t289;
t400 = (-t392 + t421) * t290;
t399 = t329 * t294 + t296 * t384;
t388 = qJD(2) * t330;
t398 = t333 * t294 + t388 * t417;
t291 = t359 * qJD(1);
t397 = pkin(5) * t375 + t333 * t291;
t312 = pkin(5) * t403;
t396 = t329 * t296 + t312;
t325 = t330 ^ 2;
t395 = -t334 ^ 2 + t325;
t391 = qJD(2) * t286;
t390 = qJD(2) * t288;
t387 = qJD(2) * t334;
t382 = qJD(4) * t332;
t377 = t332 * t235 - t328 * t236 - t286 * t382;
t376 = qJD(3) * t419;
t374 = t311 * t381;
t373 = t329 * t387;
t371 = t311 * t386;
t370 = t311 * t384;
t363 = -qJD(3) * t280 - t267;
t362 = qJD(4) * t227 + t217;
t360 = pkin(3) * t386 - t392 * t418 - t318;
t357 = g(1) * t331 - g(2) * t335;
t356 = t300 * t384 - t243;
t302 = t419 * t333;
t351 = pkin(3) * t330 - pkin(7) * t403;
t355 = qJD(1) * t351 + qJD(4) * t302 + t333 * t376 + t397;
t275 = t329 * t291;
t301 = t419 * t329;
t354 = -qJD(4) * t301 - t275 - (-t405 * pkin(5) - pkin(7) * t406) * qJD(1) - t329 * t376;
t353 = -pkin(6) * t283 + qJD(3) * t299;
t221 = t227 * t328 + t413;
t349 = -0.2e1 * pkin(1) * t380 - pkin(5) * qJDD(2);
t348 = t283 * t329 - t370;
t347 = t283 * t333 + t371;
t218 = -t288 * t383 + t377;
t337 = qJD(1) ^ 2;
t343 = pkin(1) * t337 + t358;
t342 = t384 * t330 + t373;
t336 = qJD(2) ^ 2;
t338 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t336 + t357;
t315 = -pkin(3) * t333 - pkin(2);
t295 = (pkin(5) + t418) * t330;
t285 = t333 * t296;
t274 = t329 * t331 + t333 * t402;
t273 = -t329 * t402 + t331 * t333;
t272 = t329 * t335 - t331 * t403;
t271 = t329 * t404 + t333 * t335;
t264 = t289 * t330;
t253 = pkin(3) * t342 + pkin(5) * t387;
t251 = -pkin(7) * t407 + t396;
t245 = -pkin(7) * t405 + t285 + (-pkin(3) - t417) * t334;
t226 = pkin(3) * t236 + t268;
t225 = -t383 * t407 + (t421 * t405 + t373) * t332 + t422 * t328;
t224 = -qJD(2) * t345 - t421 * t263;
t223 = -t342 * pkin(7) + (-t330 * t381 - t334 * t386) * pkin(5) + t399;
t222 = t351 * qJD(2) + (-t312 + (pkin(7) * t330 - t296) * t329) * qJD(3) + t398;
t220 = -t230 * t328 + t414;
t1 = [qJDD(1) * MDP(1) + t357 * MDP(2) + t358 * MDP(3) + (qJDD(1) * t325 + 0.2e1 * t330 * t367) * MDP(4) + 0.2e1 * (t321 * t330 - t380 * t395) * MDP(5) + (qJDD(2) * t330 + t334 * t336) * MDP(6) + (qJDD(2) * t334 - t330 * t336) * MDP(7) + (t330 * t349 + t334 * t338) * MDP(9) + (-t330 * t338 + t334 * t349) * MDP(10) + (t235 * t405 + t422 * t288) * MDP(11) + ((-t286 * t333 - t288 * t329) * t387 + (-t412 - t236 * t333 + (t286 * t329 - t408) * qJD(3)) * t330) * MDP(12) + ((-t235 - t374) * t334 + (t347 + t390) * t330) * MDP(13) + ((t311 * t389 + t236) * t334 + (-t348 - t391) * t330) * MDP(14) + (-t283 * t334 - t311 * t388) * MDP(15) + (-(-t296 * t386 + t398) * t311 + t285 * t283 - g(1) * t272 - g(2) * t274 + ((t370 + t391) * pkin(5) + (-pkin(5) * t283 + qJD(2) * t299 - t363) * t329 + t356) * t334 + (pkin(5) * t236 + t246 * qJD(2) + t268 * t329 + t299 * t384) * t330) * MDP(16) + (t399 * t311 - t396 * t283 - g(1) * t271 - g(2) * t273 + (t299 * t381 + (-t371 + t390) * pkin(5) + t344) * t334 + (-t299 * t386 - t247 * qJD(2) + t268 * t333 + (t235 - t374) * pkin(5)) * t330) * MDP(17) + (-t218 * t264 - t224 * t352) * MDP(18) + (-t218 * t263 + t219 * t264 - t224 * t238 + t225 * t352) * MDP(19) + (-t218 * t334 - t224 * t306 - t264 * t279 - t352 * t388) * MDP(20) + (t219 * t334 + t225 * t306 - t238 * t388 - t263 * t279) * MDP(21) + (-t279 * t334 - t306 * t388) * MDP(22) + (-(t222 * t332 - t223 * t328) * t306 + (t245 * t332 - t251 * t328) * t279 - t365 * t334 + t220 * t388 + t253 * t238 + t295 * t219 + t226 * t263 + t252 * t225 - g(1) * t259 - g(2) * t261 + (-(-t245 * t328 - t251 * t332) * t306 + t221 * t334) * qJD(4)) * MDP(23) + (-t221 * t388 - g(1) * t258 - g(2) * t260 + t295 * t218 + t252 * t224 - t226 * t264 - t228 * t334 - t253 * t352 + ((-qJD(4) * t251 + t222) * t306 - t245 * t279 + t216 * t334) * t328 + ((qJD(4) * t245 + t223) * t306 - t251 * t279 + t362 * t334) * t332) * MDP(24); MDP(6) * t379 + MDP(7) * t321 + qJDD(2) * MDP(8) + (t330 * t343 - t316 - t415) * MDP(9) + (t416 + (-pkin(5) * qJDD(1) + t343) * t334) * MDP(10) + (-t311 * t408 + t412) * MDP(11) + ((t235 + t411) * t333 + (-t236 + t410) * t329) * MDP(12) + ((-t288 * t330 + t311 * t403) * qJD(1) + t348) * MDP(13) + ((t286 * t330 - t311 * t406) * qJD(1) + t347) * MDP(14) + (-pkin(2) * t236 + t397 * t311 + t353 * t329 + (-t246 * t330 + (-pkin(5) * t286 - t299 * t329) * t334) * qJD(1) + t431 * t333) * MDP(16) + (-pkin(2) * t235 - t275 * t311 + t353 * t333 + (-t299 * t403 + t247 * t330 + (-t288 * t334 + t311 * t405) * pkin(5)) * qJD(1) - t431 * t329) * MDP(17) + (t218 * t290 - t352 * t401) * MDP(18) + (-t218 * t289 - t219 * t290 - t238 * t401 + t352 * t400) * MDP(19) + (t279 * t290 - t306 * t401) * MDP(20) + (-t279 * t289 + t306 * t400) * MDP(21) + ((-t301 * t332 - t302 * t328) * t279 + t315 * t219 + t226 * t289 + (t328 * t354 + t332 * t355) * t306 + t400 * t252 + t360 * t238 + t341 * t324) * MDP(23) + (-(-t301 * t328 + t302 * t332) * t279 + t315 * t218 + t226 * t290 + (-t328 * t355 + t332 * t354) * t306 + t401 * t252 - t360 * t352 - t341 * t323) * MDP(24) + (t311 * MDP(15) + MDP(20) * t352 + t238 * MDP(21) + t306 * MDP(22) - t220 * MDP(23) + t221 * MDP(24)) * t393 + (-t330 * t334 * MDP(4) + t395 * MDP(5)) * t337; t288 * t286 * MDP(11) + (-t286 ^ 2 + t288 ^ 2) * MDP(12) + (t235 - t411) * MDP(13) + (-t236 - t410) * MDP(14) + t283 * MDP(15) + (-g(1) * t273 + g(2) * t271 - t247 * t311 - t288 * t299 + (t363 + t416) * t329 - t356) * MDP(16) + (g(1) * t274 - g(2) * t272 + g(3) * t405 - t246 * t311 + t286 * t299 - t344) * MDP(17) + (t218 - t430) * MDP(20) + (-t219 + t429) * MDP(21) + ((-t229 * t328 - t413) * t306 - t221 * qJD(4) + (-t238 * t288 + t332 * t279 + t306 * t383) * pkin(3) + t426) * MDP(23) + ((t230 * t306 - t216) * t328 + (-t229 * t306 - t362) * t332 + (-t328 * t279 + t288 * t352 + t306 * t382) * pkin(3) + t427) * MDP(24) + t425; (t377 - t430) * MDP(20) + (-t364 + t429) * MDP(21) + (-t221 * t306 + t426) * MDP(23) + (-t328 * t216 - t332 * t217 - t220 * t306 + t427) * MDP(24) + (-MDP(20) * t409 + t352 * MDP(21) - t221 * MDP(23) - MDP(24) * t414) * qJD(4) + t425;];
tau = t1;
