% Calculate vector of inverse dynamics joint torques for
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:35
% EndTime: 2019-12-31 19:58:40
% DurationCPUTime: 3.18s
% Computational Cost: add. (2697->355), mult. (6277->464), div. (0->0), fcn. (4381->10), ass. (0->163)
t333 = sin(pkin(8));
t334 = cos(pkin(8));
t338 = sin(qJ(2));
t341 = cos(qJ(2));
t306 = t333 * t341 + t334 * t338;
t296 = t306 * qJD(1);
t377 = qJDD(1) * t341;
t378 = qJDD(1) * t338;
t360 = -t333 * t378 + t334 * t377;
t270 = -qJD(2) * t296 + t360;
t379 = qJD(1) * qJD(2);
t371 = t338 * t379;
t313 = t333 * t371;
t370 = t341 * t379;
t357 = t334 * t370 - t313;
t376 = pkin(2) * t371 + qJDD(3);
t411 = t341 * pkin(2);
t323 = pkin(1) + t411;
t426 = -pkin(7) * t306 - t323;
t230 = -t270 * pkin(3) - t357 * pkin(7) + t426 * qJDD(1) + t376;
t340 = cos(qJ(4));
t227 = t340 * t230;
t384 = qJD(1) * t338;
t400 = t334 * t341;
t294 = qJD(1) * t400 - t333 * t384;
t310 = -t323 * qJD(1) + qJD(3);
t243 = -pkin(3) * t294 - pkin(7) * t296 + t310;
t336 = -qJ(3) - pkin(6);
t311 = t336 * t338;
t308 = qJD(1) * t311;
t410 = qJD(2) * pkin(2);
t302 = t308 + t410;
t312 = t336 * t341;
t309 = qJD(1) * t312;
t401 = t334 * t309;
t269 = t333 * t302 - t401;
t259 = qJD(2) * pkin(7) + t269;
t337 = sin(qJ(4));
t229 = t243 * t337 + t259 * t340;
t367 = qJD(2) * t336;
t293 = -qJD(3) * t338 + t341 * t367;
t264 = qJDD(2) * pkin(2) + t293 * qJD(1) + qJDD(1) * t311;
t292 = qJD(3) * t341 + t338 * t367;
t271 = t292 * qJD(1) - qJDD(1) * t312;
t236 = t333 * t264 + t334 * t271;
t234 = qJDD(2) * pkin(7) + t236;
t346 = t306 * qJDD(1) + t357;
t380 = t340 * qJD(2);
t382 = qJD(4) * t337;
t238 = -qJD(4) * t380 - t337 * qJDD(2) + t296 * t382 - t340 * t346;
t295 = t306 * qJD(2);
t265 = qJD(1) * t295 + qJDD(4) - t360;
t281 = qJD(2) * t337 + t296 * t340;
t215 = pkin(4) * t265 + qJ(5) * t238 - t229 * qJD(4) - qJD(5) * t281 - t234 * t337 + t227;
t279 = t296 * t337 - t380;
t223 = -qJ(5) * t279 + t229;
t287 = qJD(4) - t294;
t431 = t287 * t223 + t215;
t239 = t281 * qJD(4) - t340 * qJDD(2) + t337 * t346;
t381 = qJD(4) * t340;
t353 = t337 * t230 + t340 * t234 + t243 * t381 - t259 * t382;
t216 = -qJ(5) * t239 - qJD(5) * t279 + t353;
t228 = t340 * t243 - t259 * t337;
t222 = -qJ(5) * t281 + t228;
t220 = pkin(4) * t287 + t222;
t430 = -t287 * t220 + t216;
t319 = pkin(2) * t333 + pkin(7);
t393 = qJ(5) + t319;
t429 = qJ(5) * t294 - qJD(4) * t393;
t339 = sin(qJ(1));
t342 = cos(qJ(1));
t425 = g(1) * t339 - g(2) * t342;
t363 = g(1) * t342 + g(2) * t339;
t330 = qJ(2) + pkin(8);
t325 = cos(t330);
t394 = t340 * t342;
t397 = t337 * t339;
t288 = t325 * t397 + t394;
t395 = t339 * t340;
t396 = t337 * t342;
t290 = -t325 * t396 + t395;
t424 = -g(1) * t290 + g(2) * t288;
t235 = t264 * t334 - t333 * t271;
t233 = -qJDD(2) * pkin(3) - t235;
t324 = sin(t330);
t349 = g(3) * t325 - t363 * t324;
t383 = qJD(4) * t287;
t423 = t319 * t383 + t233 + t349;
t422 = t281 ^ 2;
t421 = pkin(2) * t334;
t414 = g(3) * t324;
t412 = g(3) * t341;
t408 = t238 * t337;
t407 = t279 * t294;
t406 = t279 * t296;
t405 = t281 * t287;
t404 = t281 * t296;
t403 = t306 * t337;
t402 = t306 * t340;
t299 = t333 * t309;
t399 = t337 * t265;
t398 = t337 * t287;
t255 = t340 * t265;
t277 = t311 * t333 - t312 * t334;
t274 = t340 * t277;
t392 = -t222 + t220;
t391 = -t337 * t239 - t279 * t381;
t390 = t294 * t398 + t255;
t251 = pkin(2) * t384 + pkin(3) * t296 - pkin(7) * t294;
t273 = t308 * t334 + t299;
t389 = t337 * t251 + t340 * t273;
t305 = t333 * t338 - t400;
t267 = pkin(3) * t305 + t426;
t388 = t337 * t267 + t274;
t387 = qJD(5) * t340 + t429 * t337 - t389;
t247 = t340 * t251;
t386 = -pkin(4) * t296 - t247 + t429 * t340 + (-qJD(5) + t273) * t337;
t331 = t338 ^ 2;
t385 = -t341 ^ 2 + t331;
t375 = t338 * t410;
t250 = t292 * t334 + t293 * t333;
t298 = t305 * qJD(2);
t252 = pkin(3) * t295 + pkin(7) * t298 + t375;
t374 = t340 * t250 + t337 * t252 + t267 * t381;
t322 = pkin(4) * t340 + pkin(3);
t373 = t306 * t381;
t369 = pkin(4) * t337 - t336;
t249 = t292 * t333 - t334 * t293;
t268 = t302 * t334 + t299;
t272 = t308 * t333 - t401;
t276 = -t334 * t311 - t312 * t333;
t365 = t287 * t340;
t364 = -qJD(4) * t243 - t234;
t361 = -t259 * t381 + t227;
t335 = -qJ(5) - pkin(7);
t359 = t322 * t325 - t324 * t335;
t358 = qJ(5) * t298 - qJD(5) * t306;
t258 = -qJD(2) * pkin(3) - t268;
t356 = -0.2e1 * pkin(1) * t379 - pkin(6) * qJDD(2);
t355 = -t298 * t337 + t373;
t354 = -t298 * t340 - t306 * t382;
t352 = t287 * t258 - t319 * t265;
t350 = -t323 * qJDD(1) + t376;
t343 = qJD(2) ^ 2;
t348 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t343 + t425;
t344 = qJD(1) ^ 2;
t347 = pkin(1) * t344 - pkin(6) * qJDD(1) + t363;
t219 = pkin(4) * t239 + qJDD(5) + t233;
t320 = -pkin(3) - t421;
t314 = t342 * t323;
t304 = t393 * t340;
t303 = t393 * t337;
t291 = t325 * t394 + t397;
t289 = -t325 * t395 + t396;
t278 = t279 ^ 2;
t257 = t340 * t267;
t248 = t340 * t252;
t240 = pkin(4) * t279 + qJD(5) + t258;
t231 = -qJ(5) * t403 + t388;
t224 = pkin(4) * t305 - qJ(5) * t402 - t277 * t337 + t257;
t218 = -qJ(5) * t373 + (-qJD(4) * t277 + t358) * t337 + t374;
t217 = pkin(4) * t295 - t250 * t337 + t248 + t358 * t340 + (-t274 + (qJ(5) * t306 - t267) * t337) * qJD(4);
t1 = [qJDD(1) * MDP(1) + t425 * MDP(2) + t363 * MDP(3) + (qJDD(1) * t331 + 0.2e1 * t338 * t370) * MDP(4) + 0.2e1 * (t338 * t377 - t385 * t379) * MDP(5) + (qJDD(2) * t338 + t341 * t343) * MDP(6) + (qJDD(2) * t341 - t338 * t343) * MDP(7) + (t356 * t338 + t348 * t341) * MDP(9) + (-t348 * t338 + t356 * t341) * MDP(10) + (-t235 * t306 - t236 * t305 + t249 * t296 + t250 * t294 + t268 * t298 - t269 * t295 + t277 * t270 + t276 * t346 - t363) * MDP(11) + (t236 * t277 + t269 * t250 - t235 * t276 - t268 * t249 - t350 * t323 + t310 * t375 - g(1) * (-t323 * t339 - t336 * t342) - g(2) * (-t336 * t339 + t314)) * MDP(12) + (-t238 * t402 + t354 * t281) * MDP(13) + (-(-t279 * t340 - t281 * t337) * t298 + (t408 - t239 * t340 + (t279 * t337 - t281 * t340) * qJD(4)) * t306) * MDP(14) + (-t238 * t305 + t306 * t255 + t281 * t295 + t354 * t287) * MDP(15) + (-t239 * t305 - t279 * t295 - t355 * t287 - t306 * t399) * MDP(16) + (t265 * t305 + t287 * t295) * MDP(17) + ((-t277 * t381 + t248) * t287 + t257 * t265 + t361 * t305 + t228 * t295 + t249 * t279 + t276 * t239 + t258 * t373 - g(1) * t289 - g(2) * t291 + ((-qJD(4) * t267 - t250) * t287 - t277 * t265 + t364 * t305 + t233 * t306 - t258 * t298) * t337) * MDP(18) + (-(-t277 * t382 + t374) * t287 - t388 * t265 - t353 * t305 - t229 * t295 + t249 * t281 - t276 * t238 + t233 * t402 - g(1) * t288 - g(2) * t290 + t354 * t258) * MDP(19) + (-t217 * t281 - t218 * t279 + t224 * t238 - t231 * t239 + t425 * t324 - (-t220 * t340 - t223 * t337) * t298 + (-t215 * t340 - t216 * t337 + (t220 * t337 - t223 * t340) * qJD(4)) * t306) * MDP(20) + (t216 * t231 + t223 * t218 + t215 * t224 + t220 * t217 + t219 * (pkin(4) * t403 + t276) + t240 * (t355 * pkin(4) + t249) - g(2) * t314 + (-g(1) * t369 - g(2) * t359) * t342 + (-g(1) * (-t323 - t359) - g(2) * t369) * t339) * MDP(21); MDP(6) * t378 + MDP(7) * t377 + qJDD(2) * MDP(8) + (t347 * t338 - t412) * MDP(9) + (g(3) * t338 + t347 * t341) * MDP(10) + ((t269 - t272) * t296 + (-t273 + t268) * t294 + (t333 * t270 + (-t333 * t377 + t313 + (-t370 - t378) * t334) * t334) * pkin(2)) * MDP(11) + (t268 * t272 - t269 * t273 + (-t412 + t235 * t334 + t236 * t333 + (-qJD(1) * t310 + t363) * t338) * pkin(2)) * MDP(12) + (t281 * t365 - t408) * MDP(13) + ((-t238 + t407) * t340 - t281 * t398 + t391) * MDP(14) + (t287 * t365 + t399 - t404) * MDP(15) + (-t287 * t382 + t390 + t406) * MDP(16) - t287 * t296 * MDP(17) + (-t228 * t296 + t320 * t239 - t247 * t287 - t272 * t279 + (t273 * t287 + t352) * t337 - t423 * t340) * MDP(18) + (t229 * t296 - t320 * t238 - t272 * t281 + t389 * t287 + t423 * t337 + t352 * t340) * MDP(19) + (-t238 * t303 - t239 * t304 - t387 * t279 - t386 * t281 - t363 * t325 - t431 * t337 + t430 * t340 - t414) * MDP(20) + (t216 * t304 - t215 * t303 + t219 * (-t322 - t421) - g(3) * (t359 + t411) + (pkin(4) * t398 - t272) * t240 + t387 * t223 + t386 * t220 + t363 * (pkin(2) * t338 + t322 * t324 + t325 * t335)) * MDP(21) + (-t338 * t341 * MDP(4) + t385 * MDP(5)) * t344; (-t294 ^ 2 - t296 ^ 2) * MDP(11) + (t268 * t296 - t269 * t294 + t350 - t425) * MDP(12) + (t390 - t406) * MDP(18) - MDP(19) * t404 + t391 * MDP(20) + (-t240 * t296 - t425) * MDP(21) + (-MDP(18) * t383 - t265 * MDP(19) + MDP(20) * t405 + t430 * MDP(21)) * t337 + ((t238 + t407) * MDP(20) + t431 * MDP(21) - t287 ^ 2 * MDP(19)) * t340; t281 * t279 * MDP(13) + (-t278 + t422) * MDP(14) + (t279 * t287 - t238) * MDP(15) + (-t239 + t405) * MDP(16) + t265 * MDP(17) + (t229 * t287 - t258 * t281 + (t364 + t414) * t337 + t361 + t424) * MDP(18) + (g(1) * t291 - g(2) * t289 + t228 * t287 + t258 * t279 + t340 * t414 - t353) * MDP(19) + (pkin(4) * t238 - t392 * t279) * MDP(20) + (t392 * t223 + (-t240 * t281 + t337 * t414 + t215 + t424) * pkin(4)) * MDP(21); (-t278 - t422) * MDP(20) + (t220 * t281 + t223 * t279 + t219 + t349) * MDP(21);];
tau = t1;
