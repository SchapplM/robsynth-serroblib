% Calculate vector of inverse dynamics joint torques for
% S5PRRPR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:50
% EndTime: 2021-01-15 15:53:00
% DurationCPUTime: 3.89s
% Computational Cost: add. (1675->332), mult. (3770->445), div. (0->0), fcn. (2834->14), ass. (0->148)
t336 = qJ(4) + pkin(6);
t332 = sin(pkin(9));
t334 = cos(pkin(9));
t341 = cos(qJ(3));
t389 = qJD(2) * t341;
t381 = t334 * t389;
t338 = sin(qJ(3));
t390 = qJD(2) * t338;
t293 = t332 * t390 - t381;
t340 = cos(qJ(5));
t283 = t340 * t293;
t303 = t332 * t341 + t334 * t338;
t296 = t303 * qJD(2);
t337 = sin(qJ(5));
t403 = t296 * t337;
t250 = -t283 - t403;
t328 = qJD(3) + qJD(5);
t404 = t250 * t328;
t366 = t293 * t337 - t340 * t296;
t405 = t366 * t328;
t342 = cos(qJ(2));
t333 = sin(pkin(8));
t335 = cos(pkin(8));
t373 = g(1) * t335 + g(2) * t333;
t358 = t373 * t342;
t339 = sin(qJ(2));
t411 = g(3) * t339;
t351 = t358 + t411;
t359 = t373 * t339;
t410 = g(3) * t342;
t350 = t359 - t410;
t387 = qJD(1) * qJD(2);
t301 = qJDD(2) * pkin(6) + qJDD(1) * t339 + t342 * t387;
t354 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + t301;
t392 = qJD(1) * t339;
t376 = t336 * qJD(2) + t392;
t364 = qJD(3) * t376;
t235 = qJDD(3) * pkin(3) - t338 * t354 - t341 * t364;
t236 = -t338 * t364 + t341 * t354;
t220 = t334 * t235 - t332 * t236;
t386 = qJD(2) * qJD(3);
t380 = t338 * t386;
t311 = t332 * t380;
t379 = t341 * t386;
t256 = t303 * qJDD(2) + t334 * t379 - t311;
t216 = qJDD(3) * pkin(4) - pkin(7) * t256 + t220;
t221 = t332 * t235 + t334 * t236;
t295 = t303 * qJD(3);
t384 = qJDD(2) * t341;
t314 = t334 * t384;
t385 = qJDD(2) * t338;
t255 = qJD(2) * t295 + t332 * t385 - t314;
t217 = -pkin(7) * t255 + t221;
t288 = t376 * t338;
t408 = qJD(3) * pkin(3);
t275 = -t288 + t408;
t289 = t376 * t341;
t400 = t334 * t289;
t238 = t332 * t275 + t400;
t415 = pkin(7) * t293;
t227 = t238 - t415;
t323 = pkin(3) * t341 + pkin(2);
t391 = qJD(1) * t342;
t299 = -qJD(2) * t323 + qJD(4) - t391;
t259 = pkin(4) * t293 + t299;
t329 = qJ(3) + pkin(9);
t326 = qJ(5) + t329;
t318 = sin(t326);
t319 = cos(t326);
t388 = qJD(5) * t337;
t399 = t335 * t342;
t401 = t333 * t342;
t422 = -t259 * t250 - g(1) * (-t318 * t333 - t319 * t399) - g(2) * (t318 * t335 - t319 * t401) - t337 * t216 - t340 * t217 + t227 * t388 + t319 * t411;
t421 = t259 * t366 - g(1) * (-t318 * t399 + t319 * t333) - g(2) * (-t318 * t401 - t319 * t335) + t340 * t216 - t337 * t217 + t318 * t411;
t327 = qJDD(3) + qJDD(5);
t420 = t327 * MDP(20) + t250 * MDP(16) * t366 + (-t250 ^ 2 + t366 ^ 2) * MDP(17);
t378 = qJD(3) * t336;
t290 = qJD(4) * t341 - t338 * t378;
t291 = -qJD(4) * t338 - t341 * t378;
t396 = -t290 * t332 + t334 * t291 + t303 * t391;
t402 = t332 * t338;
t302 = -t334 * t341 + t402;
t355 = t302 * t342;
t395 = qJD(1) * t355 + t334 * t290 + t332 * t291;
t419 = t338 * t408 - t392;
t298 = t302 * qJD(3);
t343 = qJD(3) ^ 2;
t374 = -qJDD(1) * t342 + t339 * t387;
t418 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t343 + t339 * (t373 + t387) - t374 - t410;
t377 = t340 * t255 + t256 * t337;
t219 = -qJD(5) * t366 + t377;
t417 = pkin(3) * t332;
t416 = pkin(3) * t338;
t414 = pkin(7) * t296;
t409 = qJD(2) * pkin(2);
t271 = t332 * t289;
t237 = t334 * t275 - t271;
t225 = qJD(3) * pkin(4) + t237 - t414;
t406 = t225 * t340;
t344 = qJD(2) ^ 2;
t398 = t341 * t344;
t397 = qJDD(1) - g(3);
t240 = -t334 * t288 - t271;
t309 = t336 * t338;
t310 = t336 * t341;
t262 = -t332 * t309 + t334 * t310;
t330 = t338 ^ 2;
t394 = -t341 ^ 2 + t330;
t382 = -qJD(5) * t283 - t337 * t255 + t340 * t256;
t239 = t288 * t332 - t400;
t261 = -t334 * t309 - t310 * t332;
t375 = pkin(4) * t295 + t419;
t372 = g(1) * t333 - g(2) * t335;
t246 = -pkin(7) * t302 + t262;
t371 = -pkin(7) * t298 + qJD(5) * t246 - t396;
t245 = -pkin(7) * t303 + t261;
t370 = pkin(7) * t295 - qJD(5) * t245 - t395;
t369 = -t225 * t337 - t227 * t340;
t286 = t303 * t339;
t287 = t302 * t339;
t368 = -t286 * t340 + t287 * t337;
t367 = -t286 * t337 - t287 * t340;
t257 = t340 * t302 + t303 * t337;
t258 = -t302 * t337 + t303 * t340;
t363 = qJDD(3) * t341 - t338 * t343;
t362 = qJDD(3) * t338 + t341 * t343;
t320 = pkin(3) * t334 + pkin(4);
t361 = t320 * t337 + t340 * t417;
t360 = t320 * t340 - t337 * t417;
t357 = t372 * t341;
t356 = pkin(3) * t380 + qJDD(4) + t374;
t218 = -t296 * t388 + t382;
t313 = -t391 - t409;
t347 = -pkin(6) * qJDD(3) + (t313 + t391 - t409) * qJD(3);
t260 = -qJDD(2) * t323 + t356;
t345 = -t313 * qJD(2) - t301 + t351;
t325 = cos(t329);
t324 = sin(t329);
t279 = pkin(4) * t302 - t323;
t265 = pkin(3) * t390 + pkin(4) * t296;
t242 = -qJD(2) * t355 - qJD(3) * t286;
t241 = -t296 * t342 + t298 * t339;
t229 = t240 - t414;
t228 = t239 + t415;
t226 = pkin(4) * t255 + t260;
t223 = qJD(5) * t258 + t340 * t295 - t298 * t337;
t222 = -qJD(5) * t257 - t295 * t337 - t298 * t340;
t1 = [t397 * MDP(1) + (qJD(3) * t241 - qJDD(3) * t286) * MDP(12) + (-qJD(3) * t242 + qJDD(3) * t287) * MDP(13) + (-t241 * t296 - t242 * t293 + t255 * t287 + t256 * t286) * MDP(14) + (-t220 * t286 - t221 * t287 + t237 * t241 + t238 * t242 - g(3)) * MDP(15) + ((-qJD(5) * t367 + t241 * t340 - t242 * t337) * t328 + t368 * t327) * MDP(21) + (-(qJD(5) * t368 + t241 * t337 + t242 * t340) * t328 - t367 * t327) * MDP(22) + (qJDD(2) * MDP(3) - t344 * MDP(4) + (-0.2e1 * t380 + t384) * MDP(10) + (-0.2e1 * t379 - t385) * MDP(11) - t255 * MDP(12) - t256 * MDP(13) - t260 * MDP(15) - t219 * MDP(21) - t218 * MDP(22)) * t342 + (-t344 * MDP(3) - qJDD(2) * MDP(4) + (-t362 - t398) * MDP(10) + (t338 * t344 - t363) * MDP(11) + (t293 * MDP(12) + t296 * MDP(13) + t299 * MDP(15) - MDP(21) * t250 - MDP(22) * t366) * qJD(2)) * t339; qJDD(2) * MDP(2) + (t342 * t397 + t359) * MDP(3) + (-t339 * t397 + t358) * MDP(4) + (qJDD(2) * t330 + 0.2e1 * t338 * t379) * MDP(5) + 0.2e1 * (t338 * t384 - t386 * t394) * MDP(6) + t362 * MDP(7) + t363 * MDP(8) + (t347 * t338 + t341 * t418) * MDP(10) + (-t338 * t418 + t347 * t341) * MDP(11) + (-t293 * t392 + qJDD(3) * t261 - t255 * t323 + t260 * t302 + t295 * t299 + t350 * t325 + (t293 * t416 + t396) * qJD(3)) * MDP(12) + (-t296 * t392 - qJDD(3) * t262 - t256 * t323 + t260 * t303 - t298 * t299 - t350 * t324 + (t296 * t416 - t395) * qJD(3)) * MDP(13) + (-t220 * t303 - t221 * t302 + t237 * t298 - t238 * t295 - t255 * t262 - t256 * t261 - t293 * t395 - t296 * t396 - t351) * MDP(14) + (t221 * t262 + t220 * t261 - t260 * t323 - g(3) * (t323 * t342 + t336 * t339) + t419 * t299 + t395 * t238 + t396 * t237 + t373 * (t323 * t339 - t336 * t342)) * MDP(15) + (t218 * t258 - t222 * t366) * MDP(16) + (-t218 * t257 - t219 * t258 + t222 * t250 + t223 * t366) * MDP(17) + (t222 * t328 + t258 * t327) * MDP(18) + (-t223 * t328 - t257 * t327) * MDP(19) + ((t245 * t340 - t246 * t337) * t327 + t279 * t219 + t226 * t257 + t259 * t223 + (t337 * t370 - t340 * t371) * t328 - t375 * t250 + t350 * t319) * MDP(21) + (-(t245 * t337 + t246 * t340) * t327 + t279 * t218 + t226 * t258 + t259 * t222 + (t337 * t371 + t340 * t370) * t328 - t375 * t366 - t350 * t318) * MDP(22); -t338 * MDP(5) * t398 + t394 * MDP(6) * t344 + MDP(7) * t385 + MDP(8) * t384 + qJDD(3) * MDP(9) + (t338 * t345 - t357) * MDP(10) + (t338 * t372 + t341 * t345) * MDP(11) + (-t239 * qJD(3) - t299 * t296 - g(1) * (-t324 * t399 + t325 * t333) - g(2) * (-t324 * t401 - t325 * t335) + t324 * t411 + (qJDD(3) * t334 - t293 * t390) * pkin(3) + t220) * MDP(12) + (t240 * qJD(3) + t299 * t293 - g(1) * (-t324 * t333 - t325 * t399) - g(2) * (t324 * t335 - t325 * t401) + t325 * t411 + (-qJDD(3) * t332 - t296 * t390) * pkin(3) - t221) * MDP(13) + ((t238 + t239) * t296 + (-t237 + t240) * t293 + (-t255 * t332 - t256 * t334) * pkin(3)) * MDP(14) + (-t237 * t239 - t238 * t240 + (t220 * t334 + t221 * t332 - t357 + (-t299 * qJD(2) + t351) * t338) * pkin(3)) * MDP(15) + (t218 - t404) * MDP(18) + (-t219 - t405) * MDP(19) + (t360 * t327 - (t228 * t340 - t229 * t337) * t328 + t265 * t250 + (-t328 * t361 + t369) * qJD(5) + t421) * MDP(21) + (-t361 * t327 + (t228 * t337 + t229 * t340) * t328 + t265 * t366 + (-t328 * t360 - t406) * qJD(5) + t422) * MDP(22) + t420; -t314 * MDP(12) - t311 * MDP(13) + (-t293 ^ 2 - t296 ^ 2) * MDP(14) + (t237 * t296 + t238 * t293 - t350 + t356) * MDP(15) + (t219 - t405) * MDP(21) + (t218 + t404) * MDP(22) + (MDP(12) * t402 + t303 * MDP(13) - t323 * MDP(15)) * qJDD(2) + ((t332 * t389 + t334 * t390 + t296) * MDP(12) + (-t293 + t381) * MDP(13)) * qJD(3); (t382 - t404) * MDP(18) + (-t377 - t405) * MDP(19) + (-t328 * t369 + t421) * MDP(21) + ((-t227 * t337 + t406) * t328 + t422) * MDP(22) + (-MDP(18) * t403 + t366 * MDP(19) + t369 * MDP(21) - MDP(22) * t406) * qJD(5) + t420;];
tau = t1;
