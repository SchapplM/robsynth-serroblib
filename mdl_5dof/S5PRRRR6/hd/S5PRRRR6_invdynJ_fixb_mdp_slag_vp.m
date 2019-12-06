% Calculate vector of inverse dynamics joint torques for
% S5PRRRR6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:22
% EndTime: 2019-12-05 17:10:26
% DurationCPUTime: 2.70s
% Computational Cost: add. (1437->256), mult. (2316->344), div. (0->0), fcn. (1716->14), ass. (0->143)
t331 = qJ(2) + qJ(3);
t321 = cos(t331);
t421 = g(3) * t321;
t319 = sin(t331);
t332 = sin(pkin(9));
t333 = cos(pkin(9));
t373 = g(1) * t333 + g(2) * t332;
t438 = t373 * t319;
t439 = -t438 + t421;
t339 = cos(qJ(4));
t397 = qJD(4) * t339;
t437 = qJD(5) * t339 + t397;
t341 = cos(qJ(2));
t317 = t341 * qJDD(1);
t337 = sin(qJ(2));
t393 = qJD(1) * qJD(2);
t268 = qJDD(2) * pkin(2) - t337 * t393 + t317;
t400 = qJD(1) * t341;
t307 = qJD(2) * pkin(2) + t400;
t336 = sin(qJ(3));
t340 = cos(qJ(3));
t401 = qJD(1) * t337;
t380 = qJD(3) * t401;
t399 = qJD(3) * t336;
t434 = qJDD(1) * t337 + t341 * t393;
t370 = t307 * t399 + (-t268 + t380) * t340 + t434 * t336;
t325 = qJDD(2) + qJDD(3);
t425 = pkin(3) * t325;
t236 = t370 - t425;
t436 = t236 + t421;
t338 = cos(qJ(5));
t433 = t437 * t338;
t405 = t338 * t339;
t334 = sin(qJ(5));
t335 = sin(qJ(4));
t409 = t334 * t335;
t274 = -t405 + t409;
t276 = t334 * t339 + t335 * t338;
t327 = qJD(2) + qJD(3);
t326 = qJD(4) + qJD(5);
t430 = t326 * t276;
t231 = t274 * t325 + t327 * t430;
t275 = t336 * t337 - t340 * t341;
t431 = t275 * t327;
t359 = t326 * t274;
t277 = t336 * t341 + t337 * t340;
t266 = t277 * qJD(1);
t375 = pkin(2) * t399 - t266;
t313 = pkin(2) * t336 + pkin(7);
t426 = pkin(2) * t340;
t314 = -pkin(3) - t426;
t342 = qJD(4) ^ 2;
t429 = t313 * t342 + t314 * t325 + t375 * t327;
t428 = (qJD(3) * t307 + t434) * t340 + t268 * t336;
t427 = -pkin(7) - pkin(8);
t424 = pkin(3) * t327;
t310 = g(3) * t319;
t420 = -pkin(8) - t313;
t261 = t307 * t336 + t340 * t401;
t254 = pkin(7) * t327 + t261;
t379 = pkin(8) * t327 + t254;
t242 = t379 * t339;
t419 = t242 * t338;
t418 = t261 * t327;
t324 = qJDD(4) + qJDD(5);
t416 = t274 * t324;
t415 = t276 * t324;
t414 = t321 * t332;
t413 = t321 * t333;
t412 = t325 * t335;
t411 = t325 * t339;
t410 = t327 * t335;
t406 = t335 * t339;
t404 = qJDD(1) - g(3);
t398 = qJD(4) * t335;
t389 = pkin(4) * t398;
t403 = t389 + t375;
t328 = t335 ^ 2;
t402 = -t339 ^ 2 + t328;
t396 = qJD(5) * t334;
t395 = qJD(5) * t338;
t390 = qJD(3) * t426;
t388 = t327 * t409;
t309 = t336 * t401;
t260 = t307 * t340 - t309;
t253 = -t260 - t424;
t387 = t253 * t397 + t436 * t335;
t386 = t253 * t398 + t339 * t438;
t315 = -pkin(4) * t339 - pkin(3);
t385 = qJD(4) * t427;
t384 = t327 * t397;
t377 = qJD(4) * t420;
t304 = t336 * t380;
t376 = g(1) * t413 + g(2) * t414 + t304 + t310;
t374 = -t261 + t389;
t372 = g(1) * t332 - g(2) * t333;
t241 = t379 * t335;
t240 = qJD(4) * pkin(4) - t241;
t369 = -t240 * t334 - t419;
t248 = t327 * t277;
t368 = -t248 * t327 - t275 * t325;
t269 = t420 * t335;
t322 = t339 * pkin(8);
t270 = t313 * t339 + t322;
t367 = t269 * t338 - t270 * t334;
t366 = t269 * t334 + t270 * t338;
t297 = t427 * t335;
t298 = pkin(7) * t339 + t322;
t365 = t297 * t338 - t298 * t334;
t364 = t297 * t334 + t298 * t338;
t361 = t327 * t398 - t411;
t357 = pkin(7) * t342 - t418 - t425;
t228 = t361 * pkin(4) + t236;
t243 = t315 * t327 - t260;
t330 = qJ(4) + qJ(5);
t320 = cos(t330);
t356 = t228 * t274 + t243 * t430 - t320 * t439;
t355 = t277 * t342 - t368;
t230 = t276 * t325 - t326 * t388 + t433 * t327;
t262 = -t327 * t405 + t388;
t264 = t276 * t327;
t354 = t264 * t262 * MDP(15) + (t262 * t326 + t230) * MDP(17) + (t264 * t326 - t231) * MDP(18) + (-t262 ^ 2 + t264 ^ 2) * MDP(16) + t324 * MDP(19);
t353 = -pkin(7) * qJDD(4) + (t260 - t424) * qJD(4);
t352 = 0.2e1 * t431 * qJD(4) - qJDD(4) * t277;
t351 = (-t230 * t274 - t231 * t276 + t262 * t359 - t264 * t430) * MDP(16) + (t230 * t276 - t264 * t359) * MDP(15) + (-t326 * t359 + t415) * MDP(17) + (-t326 * t430 - t416) * MDP(18) + 0.2e1 * (-t402 * t327 * qJD(4) + t325 * t406) * MDP(9) + (t325 * t328 + 0.2e1 * t335 * t384) * MDP(8) + (qJDD(4) * t335 + t339 * t342) * MDP(10) + (qJDD(4) * t339 - t335 * t342) * MDP(11) + t325 * MDP(5);
t350 = -t370 - t439;
t318 = sin(t330);
t348 = t228 * t276 - t243 * t359 + t439 * t318;
t267 = t340 * t400 - t309;
t347 = -qJDD(4) * t313 + (t314 * t327 + t267 - t390) * qJD(4);
t235 = pkin(7) * t325 - t304 + t428;
t346 = -t253 * t327 + t373 * t321 - t235 + t310;
t222 = -t254 * t397 + qJDD(4) * pkin(4) - t335 * t235 + (-t384 - t412) * pkin(8);
t345 = -g(1) * (-t318 * t332 - t320 * t413) - g(2) * (t318 * t333 - t320 * t414) + t243 * t262 + t242 * t396 + t320 * t310 + (-t242 * t326 - t222) * t334;
t223 = -t361 * pkin(8) + t339 * t235 - t254 * t398;
t344 = -g(1) * (-t318 * t413 + t320 * t332) - g(2) * (-t318 * t414 - t320 * t333) + t369 * qJD(5) + t338 * t222 - t334 * t223 - t243 * t264 + t318 * t310;
t343 = qJD(2) ^ 2;
t290 = t315 - t426;
t281 = t339 * t385;
t280 = t335 * t385;
t252 = -t335 * t390 + t339 * t377;
t251 = t335 * t377 + t339 * t390;
t1 = [t404 * MDP(1) + (qJDD(2) * t341 - t337 * t343) * MDP(3) + (-qJDD(2) * t337 - t341 * t343) * MDP(4) + t368 * MDP(6) + (-t277 * t325 + t327 * t431) * MDP(7) + (t352 * t335 - t355 * t339) * MDP(13) + (t355 * t335 + t352 * t339) * MDP(14) + (t275 * t231 + t248 * t262 + t431 * t430 + ((t334 * t398 + t335 * t396 - t433) * t326 - t415) * t277) * MDP(20) + (t275 * t230 + t248 * t264 - t431 * t359 + (-(-t437 * t334 - t335 * t395 - t338 * t398) * t326 + t416) * t277) * MDP(21); qJDD(2) * MDP(2) + (t267 * t327 + (-pkin(2) * t325 - t268) * t336 + ((-pkin(2) * t327 - t307) * qJD(3) - t434) * t340 + t376) * MDP(7) + (t347 * t335 + (-t436 - t429) * t339 + t386) * MDP(13) + (t347 * t339 + (-t438 + t429) * t335 + t387) * MDP(14) + (-(t367 * qJD(5) + t251 * t338 + t252 * t334) * t326 - t366 * t324 + t290 * t230 - t267 * t359 + t403 * t264 + t348) * MDP(21) + (-g(3) * t341 + t373 * t337 + t317) * MDP(3) + (t266 * t327 + (t325 * t340 - t327 * t399) * pkin(2) + t350) * MDP(6) + (-t404 * t337 + t373 * t341) * MDP(4) + t351 + ((-t366 * qJD(5) - t251 * t334 + t252 * t338) * t326 + t367 * t324 + t290 * t231 + t267 * t430 + t403 * t262 + t356) * MDP(20); (-(t365 * qJD(5) + t280 * t338 + t281 * t334) * t326 - t364 * t324 + t315 * t230 + t374 * t264 - t260 * t359 + t348) * MDP(21) + (t353 * t339 + (-t438 + t357) * t335 + t387) * MDP(14) + (t260 * t327 + t376 - t428) * MDP(7) + (t350 + t418) * MDP(6) + ((-t364 * qJD(5) - t280 * t334 + t281 * t338) * t326 + t365 * t324 + t315 * t231 + t374 * t262 + t260 * t430 + t356) * MDP(20) + (t353 * t335 + (-t357 - t436) * t339 + t386) * MDP(13) + t351; MDP(10) * t412 + MDP(11) * t411 + qJDD(4) * MDP(12) + (t346 * t335 - t372 * t339) * MDP(13) + (t372 * t335 + t346 * t339) * MDP(14) + (-(t241 * t334 - t419) * t326 + (-t262 * t410 + t338 * t324 - t326 * t396) * pkin(4) + t344) * MDP(20) + ((-qJD(5) * t240 - t241 * t326 - t223) * t338 + (-t264 * t410 - t334 * t324 - t326 * t395) * pkin(4) + t345) * MDP(21) + t354 + (-MDP(8) * t406 + t402 * MDP(9)) * t327 ^ 2; (-t369 * t326 + t344) * MDP(20) + ((-t223 + (-qJD(5) + t326) * t240) * t338 + t345) * MDP(21) + t354;];
tau = t1;
