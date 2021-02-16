% Calculate Coriolis joint torque vector for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:27
% EndTime: 2021-01-16 02:21:36
% DurationCPUTime: 3.82s
% Computational Cost: add. (2351->361), mult. (6205->489), div. (0->0), fcn. (4647->10), ass. (0->168)
t328 = sin(pkin(11));
t333 = sin(qJ(3));
t396 = qJD(2) * t333;
t330 = cos(pkin(11));
t336 = cos(qJ(3));
t407 = t330 * t336;
t298 = -qJD(2) * t407 + t328 * t396;
t335 = cos(qJ(6));
t332 = sin(qJ(6));
t392 = qJD(3) * t332;
t279 = -t335 * t298 + t392;
t412 = t328 * t336;
t306 = t330 * t333 + t412;
t301 = t306 * qJD(2);
t440 = qJD(6) + t301;
t442 = t279 * t440;
t281 = qJD(3) * t335 + t298 * t332;
t370 = t440 * t281;
t334 = sin(qJ(2));
t329 = sin(pkin(6));
t398 = qJD(1) * t329;
t381 = t334 * t398;
t391 = qJD(3) * t333;
t441 = pkin(3) * t391 - t381;
t383 = MDP(12) - MDP(17);
t439 = MDP(13) - MDP(18);
t438 = MDP(14) + MDP(16);
t437 = t332 * t440;
t436 = (t333 ^ 2 - t336 ^ 2) * MDP(6);
t423 = qJ(4) + pkin(8);
t374 = qJD(3) * t423;
t294 = qJD(4) * t336 - t333 * t374;
t295 = -qJD(4) * t333 - t336 * t374;
t337 = cos(qJ(2));
t380 = t337 * t398;
t402 = t294 * t328 - t330 * t295 - t306 * t380;
t433 = -t328 * t333 + t407;
t401 = t294 * t330 + t295 * t328 - t433 * t380;
t390 = qJD(3) * t336;
t303 = -t328 * t391 + t330 * t390;
t434 = qJ(5) * t303 + qJD(5) * t306 - t441;
t309 = qJD(2) * pkin(8) + t381;
t369 = qJ(4) * qJD(2) + t309;
t331 = cos(pkin(6));
t397 = qJD(1) * t331;
t379 = t333 * t397;
t274 = t369 * t336 + t379;
t267 = t328 * t274;
t317 = t336 * t397;
t273 = -t369 * t333 + t317;
t244 = t273 * t330 - t267;
t387 = -qJD(5) + t244;
t271 = qJD(3) * pkin(3) + t273;
t408 = t330 * t274;
t240 = t328 * t271 + t408;
t237 = -qJD(3) * qJ(5) - t240;
t425 = pkin(5) * t298;
t231 = -t237 - t425;
t242 = t273 * t328 + t408;
t384 = qJD(2) * qJD(3);
t376 = t333 * t384;
t313 = t328 * t376;
t375 = t336 * t384;
t290 = t330 * t375 - t313;
t322 = -pkin(3) * t330 - pkin(4);
t318 = -pkin(9) + t322;
t432 = t318 * t290 + (t231 - t242 + t425) * t440;
t297 = t301 ^ 2;
t428 = pkin(4) + pkin(9);
t427 = pkin(3) * t333;
t300 = t306 * qJD(3);
t289 = qJD(2) * t300;
t426 = pkin(4) * t289;
t424 = pkin(5) * t301;
t422 = qJD(2) * pkin(2);
t360 = qJD(4) + t380;
t253 = (-t309 * t333 + t317) * qJD(3) + (-qJ(4) * t391 + t360 * t336) * qJD(2);
t254 = (-t309 * t336 - t379) * qJD(3) + (-qJ(4) * t390 - t360 * t333) * qJD(2);
t227 = t253 * t328 - t330 * t254;
t411 = t329 * t334;
t304 = t331 * t333 + t336 * t411;
t355 = t331 * t336 - t333 * t411;
t263 = t304 * t328 - t330 * t355;
t421 = t227 * t263;
t311 = t423 * t333;
t312 = t423 * t336;
t277 = t330 * t311 + t312 * t328;
t420 = t227 * t277;
t323 = -pkin(3) * t336 - pkin(2);
t291 = t323 * qJD(2) + qJD(4) - t380;
t342 = -qJ(5) * t301 + t291;
t252 = pkin(4) * t298 + t342;
t419 = t252 * t301;
t389 = qJD(6) * t332;
t388 = qJD(6) * t335;
t400 = t332 * t289 + t298 * t388;
t256 = -qJD(3) * t389 + t400;
t418 = t256 * t335;
t417 = t279 * t298;
t416 = t281 * t298;
t415 = t290 * t332;
t414 = t433 * t332;
t410 = t329 * t337;
t339 = qJD(2) ^ 2;
t409 = t329 * t339;
t338 = qJD(3) ^ 2;
t406 = t333 * t338;
t284 = t335 * t290;
t405 = t336 * t338;
t404 = pkin(5) * t303 + t402;
t403 = -pkin(5) * t300 + t401;
t228 = t330 * t253 + t328 * t254;
t395 = qJD(2) * t334;
t378 = t329 * t395;
t296 = pkin(3) * t376 + qJD(1) * t378;
t386 = t424 - t387;
t324 = pkin(3) * t396;
t382 = t334 * t409;
t377 = qJD(2) * t410;
t373 = qJ(5) * t298 + t324;
t225 = pkin(5) * t290 + t227;
t365 = -qJ(5) * t290 + t296;
t347 = -qJD(5) * t301 + t365;
t229 = t428 * t289 + t347;
t372 = t335 * t225 - t229 * t332;
t239 = t271 * t330 - t267;
t371 = t335 * t440;
t368 = MDP(25) * t440;
t367 = t336 * t377;
t366 = t333 * t377;
t364 = -t428 * t300 + t434;
t363 = -pkin(4) * t300 + t434;
t361 = qJD(5) - t239;
t226 = -qJD(3) * qJD(5) - t228;
t230 = -t428 * qJD(3) + t361 + t424;
t238 = t428 * t298 + t342;
t221 = t230 * t335 - t238 * t332;
t222 = t230 * t332 + t238 * t335;
t278 = -t311 * t328 + t312 * t330;
t358 = -qJ(5) * t306 + t323;
t357 = -t263 * t332 + t335 * t410;
t356 = t263 * t335 + t332 * t410;
t352 = t300 * t332 - t388 * t433;
t350 = -qJD(3) * t242 + t227;
t349 = qJD(2) * t422;
t223 = -pkin(5) * t289 - t226;
t348 = t223 + (-qJD(6) * t318 + t428 * t301 + t373) * t440;
t261 = pkin(5) * t306 + t277;
t346 = -t223 * t433 + t231 * t300 - t261 * t290;
t345 = t304 * qJD(3);
t344 = -0.2e1 * qJD(3) * t422;
t341 = -t345 - t366;
t340 = t227 * t306 + t277 * t290 - t278 * t289 - t401 * t298 + t402 * t301;
t320 = pkin(3) * t328 + qJ(5);
t292 = qJD(3) * t298;
t283 = t335 * t289;
t272 = t355 * qJD(3) + t367;
t265 = -pkin(4) * t433 + t358;
t264 = t330 * t304 + t328 * t355;
t262 = pkin(5) * t433 + t278;
t258 = pkin(4) * t301 + t373;
t257 = qJD(6) * t281 - t283;
t255 = -t428 * t433 + t358;
t243 = t330 * t272 + t328 * t341;
t241 = t272 * t328 - t330 * t341;
t236 = -qJD(3) * pkin(4) + t361;
t234 = t347 + t426;
t1 = [-MDP(3) * t382 - t337 * MDP(4) * t409 + (-t336 * t382 + (-t345 - 0.2e1 * t366) * qJD(3)) * MDP(10) + (t333 * t382 + (-t272 - t367) * qJD(3)) * MDP(11) + (t421 + t228 * t264 - t239 * t241 + t240 * t243 + (t291 * t395 - t296 * t337) * t329) * MDP(15) + (-t226 * t264 + t421 + t236 * t241 - t237 * t243 + (-t234 * t337 + t252 * t395) * t329) * MDP(19) + ((t357 * qJD(6) + t241 * t335 - t332 * t378) * t440 + t356 * t290 + t243 * t279 + t264 * t257) * MDP(25) + (-(t356 * qJD(6) + t241 * t332 + t335 * t378) * t440 + t357 * t290 + t243 * t281 + t264 * t256) * MDP(26) + t438 * (t241 * t301 - t243 * t298 + t263 * t290 - t264 * t289) - t439 * ((t290 * t337 - t301 * t395) * t329 + qJD(3) * t243) + t383 * ((-t289 * t337 + t298 * t395) * t329 - qJD(3) * t241); 0.2e1 * t333 * MDP(5) * t375 - 0.2e1 * t384 * t436 + MDP(7) * t405 - MDP(8) * t406 + (-pkin(8) * t405 + t333 * t344) * MDP(10) + (pkin(8) * t406 + t336 * t344) * MDP(11) + (-t298 * t381 + t289 * t323 + t291 * t300 - t296 * t433 + (t298 * t427 - t402) * qJD(3)) * MDP(12) + (-t301 * t381 + t290 * t323 + t291 * t303 + t296 * t306 + (t301 * t427 - t401) * qJD(3)) * MDP(13) + (t228 * t433 - t239 * t303 - t240 * t300 + t340) * MDP(14) + (t228 * t278 - t402 * t239 + t401 * t240 + t441 * t291 + t296 * t323 + t420) * MDP(15) + (-t226 * t433 + t236 * t303 + t237 * t300 + t340) * MDP(16) + (t402 * qJD(3) + t234 * t433 - t252 * t300 - t265 * t289 + t363 * t298) * MDP(17) + (t401 * qJD(3) - t234 * t306 - t252 * t303 - t265 * t290 + t363 * t301) * MDP(18) + (-t226 * t278 + t234 * t265 + t402 * t236 - t401 * t237 - t363 * t252 + t420) * MDP(19) + (-t256 * t414 + t352 * t281) * MDP(20) + ((-t279 * t332 + t281 * t335) * t300 - (t418 - t257 * t332 + (-t279 * t335 - t281 * t332) * qJD(6)) * t433) * MDP(21) + (t256 * t306 + t281 * t303 - t290 * t414 + t352 * t440) * MDP(22) + (-t433 * t284 - t257 * t306 - t279 * t303 + (t300 * t335 + t389 * t433) * t440) * MDP(23) + (t290 * t306 + t303 * t440) * MDP(24) + (-t255 * t415 + t372 * t306 + t221 * t303 + t262 * t257 - t346 * t335 + (t364 * t332 + t404 * t335) * t440 + t403 * t279 + ((-t255 * t335 - t261 * t332) * t440 - t222 * t306 - t231 * t414) * qJD(6)) * MDP(25) + (-t222 * t303 + t262 * t256 + t403 * t281 + (-t255 * t290 - (qJD(6) * t230 + t229) * t306 - t231 * qJD(6) * t433 + (-qJD(6) * t261 + t364) * t440) * t335 + (-(-qJD(6) * t238 + t225) * t306 + (qJD(6) * t255 - t404) * t440 + t346) * t332) * MDP(26); t339 * t436 + t336 * t349 * MDP(11) + (-t291 * t301 - t298 * t324 - t350) * MDP(12) + (qJD(3) * t244 + t291 * t298 - t301 * t324 - t228) * MDP(13) + ((t240 - t242) * t301 + (-t239 + t244) * t298 + (-t289 * t328 - t290 * t330) * pkin(3)) * MDP(14) + (t239 * t242 - t240 * t244 + (-t227 * t330 + t228 * t328 - t291 * t396) * pkin(3)) * MDP(15) + (-t289 * t320 + t290 * t322 + (-t237 - t242) * t301 + (t236 + t387) * t298) * MDP(16) + (t258 * t298 + t350 + t419) * MDP(17) + (-t252 * t298 + t258 * t301 + (0.2e1 * qJD(5) - t244) * qJD(3) + t228) * MDP(18) + (-t226 * t320 + t227 * t322 - t236 * t242 + t387 * t237 - t252 * t258) * MDP(19) + (-t332 * t370 + t418) * MDP(20) + ((-t257 - t370) * t335 + (-t256 + t442) * t332) * MDP(21) + (-t437 * t440 + t284 + t416) * MDP(22) + (-t371 * t440 - t415 - t417) * MDP(23) + t440 * t298 * MDP(24) + (t221 * t298 + t320 * t257 + t386 * t279 + t348 * t332 + t335 * t432) * MDP(25) + (-t222 * t298 + t320 * t256 + t386 * t281 - t332 * t432 + t348 * t335) * MDP(26) + (-t339 * t336 * MDP(5) + MDP(10) * t349) * t333; -t313 * MDP(13) + (t239 * t301 + t240 * t298 + t296) * MDP(15) + (t292 + t313) * MDP(18) + (t426 - t237 * t298 + (-qJD(5) - t236) * t301 + t365) * MDP(19) + (-t415 + t417) * MDP(25) + (-t284 + t416) * MDP(26) + (MDP(26) * t437 - t335 * t368) * t440 + (-t298 * MDP(13) + t383 * t301 + (t383 * t412 + (t383 * t333 + t439 * t336) * t330) * qJD(2)) * qJD(3) + t438 * (-t298 ^ 2 - t297); (t290 + t292) * MDP(16) - t301 * t298 * MDP(17) + (-t297 - t338) * MDP(18) + (qJD(3) * t237 + t227 + t419) * MDP(19) + (-qJD(3) * t279 + t284) * MDP(25) + (-qJD(3) * t281 - t415) * MDP(26) + (-MDP(26) * t371 - t332 * t368) * t440; t281 * t279 * MDP(20) + (-t279 ^ 2 + t281 ^ 2) * MDP(21) + (t400 + t442) * MDP(22) + (t283 + t370) * MDP(23) + t290 * MDP(24) + (t222 * t440 - t231 * t281 + t372) * MDP(25) + (t221 * t440 - t225 * t332 - t229 * t335 + t231 * t279) * MDP(26) + (-MDP(22) * t392 - t281 * MDP(23) - t222 * MDP(25) - t221 * MDP(26)) * qJD(6);];
tauc = t1;
