% Calculate vector of inverse dynamics joint torques for
% S5PRRRP5
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:50
% EndTime: 2021-01-15 16:33:58
% DurationCPUTime: 2.93s
% Computational Cost: add. (1807->316), mult. (3942->396), div. (0->0), fcn. (2776->10), ass. (0->155)
t329 = cos(qJ(3));
t328 = sin(qJ(2));
t393 = qJD(1) * t328;
t298 = qJD(2) * pkin(6) + t393;
t370 = pkin(7) * qJD(2) + t298;
t262 = t370 * t329;
t326 = sin(qJ(4));
t250 = t326 * t262;
t327 = sin(qJ(3));
t261 = t370 * t327;
t253 = qJD(3) * pkin(3) - t261;
t423 = cos(qJ(4));
t369 = t423 * t253 - t250;
t279 = t326 * t329 + t423 * t327;
t269 = t279 * qJD(2);
t408 = t269 * qJ(5);
t228 = -t408 + t369;
t372 = qJDD(2) * t423;
t383 = qJDD(2) * t327;
t359 = t326 * t383 - t329 * t372;
t320 = qJD(3) + qJD(4);
t436 = t320 * t279;
t236 = qJD(2) * t436 + t359;
t330 = cos(qJ(2));
t324 = sin(pkin(8));
t325 = cos(pkin(8));
t362 = g(1) * t325 + g(2) * t324;
t354 = t362 * t330;
t417 = g(3) * t328;
t435 = t354 + t417;
t317 = t329 * pkin(3);
t415 = pkin(2) + t317;
t377 = t423 * t329;
t402 = t326 * t327;
t352 = t377 - t402;
t264 = t352 * t328;
t386 = qJD(1) * qJD(2);
t275 = qJDD(2) * pkin(6) + qJDD(1) * t328 + t330 * t386;
t385 = qJD(2) * qJD(3);
t373 = t329 * t385;
t238 = -qJD(3) * t298 * t329 + qJDD(3) * pkin(3) - t327 * t275 + (-t373 - t383) * pkin(7);
t374 = t327 * t385;
t382 = qJDD(2) * t329;
t389 = qJD(3) * t327;
t239 = -t298 * t389 + t329 * t275 + (-t374 + t382) * pkin(7);
t434 = t423 * t238 - t326 * t239;
t424 = pkin(6) + pkin(7);
t378 = qJD(3) * t424;
t284 = t327 * t378;
t285 = t329 * t378;
t392 = qJD(1) * t330;
t290 = t424 * t327;
t291 = t424 * t329;
t395 = -t326 * t290 + t423 * t291;
t433 = -t395 * qJD(4) + t279 * t392 + t326 * t284 - t423 * t285;
t347 = t330 * t352;
t375 = t423 * qJD(4);
t388 = qJD(4) * t326;
t432 = qJD(1) * t347 + t423 * t284 + t326 * t285 + t290 * t375 + t291 * t388;
t318 = qJDD(3) + qJDD(4);
t422 = pkin(3) * t320;
t431 = -t326 * pkin(3) * t318 - t375 * t422;
t314 = t318 * pkin(4);
t360 = t320 * t402;
t363 = qJD(2) * t377;
t364 = t320 * t363 + t326 * t382 + t327 * t372;
t235 = qJD(2) * t360 - t364;
t412 = t235 * qJ(5);
t430 = t314 + t412;
t429 = t423 * qJD(3) + t375;
t323 = qJ(3) + qJ(4);
t315 = sin(t323);
t316 = cos(t323);
t404 = t325 * t330;
t405 = t324 * t330;
t428 = t315 * t417 - g(2) * (-t315 * t405 - t316 * t325) - g(1) * (-t315 * t404 + t316 * t324);
t308 = t328 * t386;
t331 = qJD(3) ^ 2;
t384 = qJDD(1) * t330;
t416 = g(3) * t330;
t427 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t331 + (t362 + t386) * t328 - t308 + t384 - t416;
t400 = qJDD(1) - g(3);
t426 = t362 * t328 + t400 * t330;
t425 = t269 ^ 2;
t414 = qJD(2) * pkin(2);
t411 = t236 * qJ(5);
t391 = qJD(2) * t327;
t267 = t326 * t391 - t363;
t410 = t267 * qJ(5);
t409 = t267 * t320;
t407 = t269 * t320;
t406 = t316 * t328;
t332 = qJD(2) ^ 2;
t401 = t329 * t332;
t399 = -qJ(5) * t436 + qJD(5) * t352 - t432;
t244 = -t429 * t329 + t360;
t398 = t244 * qJ(5) - t279 * qJD(5) + t433;
t227 = pkin(4) * t320 + t228;
t397 = -t228 + t227;
t396 = -t423 * t261 - t250;
t289 = pkin(4) * t316 + t317;
t321 = t327 ^ 2;
t394 = -t329 ^ 2 + t321;
t390 = qJD(2) * t328;
t276 = -qJD(2) * t415 - t392;
t371 = pkin(4) * t267 + qJD(5);
t243 = t276 + t371;
t387 = qJD(5) + t243;
t381 = pkin(3) * t391;
t380 = pkin(3) * t389;
t252 = t423 * t262;
t368 = t261 * t326 - t252;
t367 = -t423 * t290 - t291 * t326;
t242 = pkin(4) * t436 + t380;
t366 = t242 - t393;
t365 = t320 * t327;
t361 = g(1) * t324 - g(2) * t325;
t357 = -t316 * t416 + t362 * t406;
t356 = qJDD(3) * t329 - t327 * t331;
t355 = qJDD(3) * t327 + t329 * t331;
t353 = -t326 * t253 - t252;
t350 = pkin(3) * t374 - qJDD(2) * t415 + t308;
t266 = t267 ^ 2;
t348 = t269 * t267 * MDP(12) + (-t326 * qJD(2) * t365 + t364 + t409) * MDP(14) + (-t236 + t407) * MDP(15) + (-t266 + t425) * MDP(13) + t318 * MDP(16);
t299 = -t392 - t414;
t345 = -pkin(6) * qJDD(3) + (t299 + t392 - t414) * qJD(3);
t344 = pkin(4) * t236 + qJDD(5) + t350;
t343 = t353 * qJD(4) + t434;
t341 = t326 * t238 + t423 * t239 + t253 * t375 - t262 * t388;
t340 = t315 * t416 + (-qJD(1) * t269 - t362 * t315) * t328;
t339 = -t299 * qJD(2) - t275 + t435;
t337 = t343 + t428;
t336 = -g(1) * (-t315 * t324 - t316 * t404) - g(2) * (t315 * t325 - t316 * t405) + g(3) * t406 - t341;
t335 = -t276 * t269 + t337;
t334 = t276 * t267 + t336;
t333 = t387 * t267 + t336 + t411;
t319 = -qJ(5) - t424;
t312 = t423 * pkin(3) + pkin(4);
t288 = -pkin(3) * t327 - pkin(4) * t315;
t283 = pkin(2) + t289;
t263 = t279 * t328;
t256 = -pkin(4) * t352 - t415;
t248 = pkin(4) * t269 + t381;
t246 = t350 - t384;
t241 = qJ(5) * t352 + t395;
t240 = -qJ(5) * t279 + t367;
t233 = -t320 * t264 - t330 * t269;
t232 = qJD(2) * t347 - t328 * t436;
t231 = -t408 + t396;
t230 = t368 + t410;
t229 = -t353 - t410;
t224 = t344 - t384;
t219 = -t267 * qJD(5) + t341 - t411;
t218 = -t269 * qJD(5) + t343 + t430;
t1 = [t400 * MDP(1) + (-t232 * t267 - t233 * t269 - t235 * t263 - t236 * t264) * MDP(21) + (-t218 * t263 + t219 * t264 + t227 * t233 + t229 * t232 - g(3)) * MDP(22) + (MDP(17) + MDP(19)) * (t233 * t320 - t236 * t330 - t263 * t318 + t267 * t390) + (MDP(18) + MDP(20)) * (-t232 * t320 + t235 * t330 - t264 * t318 + t269 * t390) + (qJDD(2) * MDP(3) - t332 * MDP(4) + (-0.2e1 * t374 + t382) * MDP(10) + (-0.2e1 * t373 - t383) * MDP(11) - t224 * MDP(22)) * t330 + (-t332 * MDP(3) - qJDD(2) * MDP(4) + (-t355 - t401) * MDP(10) + (t327 * t332 - t356) * MDP(11) + qJD(2) * t243 * MDP(22)) * t328; qJDD(2) * MDP(2) + t426 * MDP(3) + (-t400 * t328 + t354) * MDP(4) + (qJDD(2) * t321 + 0.2e1 * t327 * t373) * MDP(5) + 0.2e1 * (t327 * t382 - t394 * t385) * MDP(6) + t355 * MDP(7) + t356 * MDP(8) + (t345 * t327 + t427 * t329) * MDP(10) + (-t427 * t327 + t345 * t329) * MDP(11) + (-t235 * t279 - t244 * t269) * MDP(12) + (-t235 * t352 - t236 * t279 + t244 * t267 - t269 * t436) * MDP(13) + (-t244 * t320 + t279 * t318) * MDP(14) + (t318 * t352 - t320 * t436) * MDP(15) + (-t415 * t236 + t276 * t436 - t246 * t352 + t367 * t318 + t357 + t433 * t320 + (t380 - t393) * t267) * MDP(17) + (t235 * t415 - t276 * t244 + t246 * t279 + t269 * t380 - t395 * t318 + t432 * t320 + t340) * MDP(18) + (-t224 * t352 + t236 * t256 + t240 * t318 + t243 * t436 + t366 * t267 + t398 * t320 + t357) * MDP(19) + (t224 * t279 - t235 * t256 - t241 * t318 + t242 * t269 - t243 * t244 - t399 * t320 + t340) * MDP(20) + (-t218 * t279 + t219 * t352 + t227 * t244 - t229 * t436 + t235 * t240 - t236 * t241 - t399 * t267 - t398 * t269 - t435) * MDP(21) + (t219 * t241 + t218 * t240 + t224 * t256 - g(3) * (t283 * t330 - t319 * t328) + t366 * t243 + t399 * t229 + t398 * t227 + t362 * (t283 * t328 + t319 * t330)) * MDP(22); -t327 * MDP(5) * t401 + t394 * MDP(6) * t332 + MDP(7) * t383 + MDP(8) * t382 + qJDD(3) * MDP(9) + (t339 * t327 - t361 * t329) * MDP(10) + (t361 * t327 + t339 * t329) * MDP(11) + (-t368 * t320 + (-t267 * t391 + t423 * t318 - t320 * t388) * pkin(3) + t335) * MDP(17) + (-t269 * t381 + t396 * t320 + t334 + t431) * MDP(18) + (-t230 * t320 - t248 * t267 + t312 * t318 - t387 * t269 + (-t252 + (-t253 - t422) * t326) * qJD(4) + t428 + t430 + t434) * MDP(19) + (t231 * t320 - t248 * t269 + t333 + t431) * MDP(20) + (t312 * t235 + (t229 + t230) * t269 + (-t227 + t231) * t267 + (-t236 * t326 + (-t423 * t267 + t269 * t326) * qJD(4)) * pkin(3)) * MDP(21) + (t218 * t312 - t229 * t231 - t227 * t230 - t243 * t248 - g(1) * (t288 * t404 + t289 * t324) - g(2) * (t288 * t405 - t289 * t325) - t288 * t417 + (t219 * t326 + (-t227 * t326 + t423 * t229) * qJD(4)) * pkin(3)) * MDP(22) + t348; (-t353 * t320 + t335) * MDP(17) + (t369 * t320 + t334) * MDP(18) + (t412 + t229 * t320 + 0.2e1 * t314 + (-t243 - t371) * t269 + t337) * MDP(19) + (-t425 * pkin(4) + t228 * t320 + t333) * MDP(20) + (pkin(4) * t235 - t397 * t267) * MDP(21) + (t397 * t229 + (-t243 * t269 + t218 + t428) * pkin(4)) * MDP(22) + t348; (t359 + t407) * MDP(19) + (t364 - t409) * MDP(20) + (-t266 - t425) * MDP(21) + (t227 * t269 + t229 * t267 + t344 - t426) * MDP(22) + (t429 * t327 * MDP(19) + (t320 * MDP(19) * t329 - MDP(20) * t365) * t326) * qJD(2);];
tau = t1;
