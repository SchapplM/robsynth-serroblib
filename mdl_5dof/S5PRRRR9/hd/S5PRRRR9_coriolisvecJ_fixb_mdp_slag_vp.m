% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:19
% EndTime: 2019-12-05 17:21:27
% DurationCPUTime: 3.58s
% Computational Cost: add. (1754->352), mult. (4495->508), div. (0->0), fcn. (3306->10), ass. (0->167)
t335 = sin(qJ(2));
t330 = sin(pkin(5));
t400 = qJD(1) * t330;
t381 = t335 * t400;
t309 = qJD(2) * pkin(7) + t381;
t331 = cos(pkin(5));
t334 = sin(qJ(3));
t338 = cos(qJ(3));
t399 = qJD(1) * t338;
t275 = -t334 * t309 + t331 * t399;
t267 = -qJD(3) * pkin(3) - t275;
t333 = sin(qJ(4));
t337 = cos(qJ(4));
t386 = t337 * qJD(3);
t397 = qJD(2) * t334;
t298 = t333 * t397 - t386;
t243 = pkin(4) * t298 + t267;
t336 = cos(qJ(5));
t395 = qJD(3) * t333;
t300 = t337 * t397 + t395;
t332 = sin(qJ(5));
t422 = t300 * t332;
t247 = t298 * t336 + t422;
t444 = t243 * t247;
t396 = qJD(2) * t338;
t322 = -qJD(4) + t396;
t318 = -qJD(5) + t322;
t443 = t247 * t318;
t350 = t298 * t332 - t300 * t336;
t442 = t318 * t350;
t385 = qJD(2) * qJD(3);
t369 = t338 * t385;
t392 = qJD(4) * t333;
t373 = t334 * t392;
t384 = qJD(3) * qJD(4);
t270 = -qJD(2) * t373 + (t369 + t384) * t337;
t418 = t331 * t334;
t320 = qJD(1) * t418;
t276 = t309 * t338 + t320;
t268 = qJD(3) * pkin(8) + t276;
t311 = -pkin(3) * t338 - pkin(8) * t334 - pkin(2);
t339 = cos(qJ(2));
t380 = t339 * t400;
t278 = qJD(2) * t311 - t380;
t417 = t333 * t278;
t240 = t268 * t337 + t417;
t398 = qJD(2) * t330;
t377 = t339 * t398;
t394 = qJD(3) * t334;
t244 = -t309 * t394 + (qJD(3) * t331 + t377) * t399;
t357 = pkin(3) * t334 - pkin(8) * t338;
t306 = t357 * qJD(3);
t274 = (t306 + t381) * qJD(2);
t363 = t244 * t333 - t274 * t337;
t342 = -qJD(4) * t240 - t363;
t370 = t334 * t385;
t221 = pkin(4) * t370 - pkin(9) * t270 + t342;
t391 = qJD(4) * t337;
t372 = t334 * t391;
t393 = qJD(3) * t338;
t376 = t333 * t393;
t344 = t372 + t376;
t271 = qJD(2) * t344 + t333 * t384;
t345 = t244 * t337 - t268 * t392 + t274 * t333 + t278 * t391;
t224 = -pkin(9) * t271 + t345;
t365 = t221 * t336 - t332 * t224;
t441 = t243 * t350 + t365;
t367 = MDP(23) * t394;
t440 = qJD(2) * t367 + (-t247 ^ 2 + t350 ^ 2) * MDP(20) - t247 * t350 * MDP(19);
t439 = MDP(5) * t334;
t302 = t332 * t337 + t333 * t336;
t281 = t302 * t334;
t328 = t334 ^ 2;
t438 = (-t338 ^ 2 + t328) * MDP(6);
t411 = t338 * t339;
t432 = pkin(7) * t333;
t437 = (-t333 * t411 + t335 * t337) * t400 - t337 * t306 - t394 * t432;
t436 = -(t333 * t335 + t337 * t411) * t400 + t333 * t306 + t311 * t391;
t435 = t338 * t386 - t373;
t434 = qJD(4) + qJD(5);
t362 = t270 * t332 + t271 * t336;
t226 = -qJD(5) * t350 + t362;
t433 = pkin(8) + pkin(9);
t431 = qJD(2) * pkin(2);
t239 = -t268 * t333 + t278 * t337;
t231 = -pkin(9) * t300 + t239;
t227 = -pkin(4) * t322 + t231;
t430 = t227 * t336;
t232 = -pkin(9) * t298 + t240;
t429 = t232 * t336;
t359 = t334 * t377;
t245 = qJD(1) * t359 + qJD(3) * t320 + t309 * t393;
t428 = t245 * t333;
t427 = t245 * t337;
t426 = t267 * t333;
t425 = t270 * t333;
t424 = t298 * t322;
t423 = t300 * t322;
t421 = t322 * t337;
t420 = t330 * t335;
t419 = t330 * t339;
t416 = t333 * t334;
t415 = t333 * t338;
t414 = t334 * t337;
t340 = qJD(3) ^ 2;
t413 = t334 * t340;
t412 = t337 * t338;
t410 = t338 * t340;
t323 = pkin(7) * t412;
t348 = pkin(4) * t334 - pkin(9) * t412;
t409 = -t348 * qJD(3) - (-t323 + (pkin(9) * t334 - t311) * t333) * qJD(4) + t437;
t408 = -t344 * pkin(9) + (-t334 * t386 - t338 * t392) * pkin(7) + t436;
t301 = t332 * t333 - t336 * t337;
t346 = t301 * t338;
t407 = qJD(2) * t346 - t301 * t434;
t406 = (-t396 + t434) * t302;
t303 = t357 * qJD(2);
t405 = t275 * t337 + t303 * t333;
t402 = t311 * t333 + t323;
t390 = qJD(5) * t332;
t389 = qJD(5) * t336;
t387 = t267 * qJD(4);
t383 = t270 * t336 - t271 * t332 - t298 * t389;
t382 = qJD(4) * t433;
t378 = t335 * t398;
t374 = t322 * t392;
t371 = t333 * t396;
t366 = MDP(16) * t394;
t228 = t232 * t390;
t364 = t332 * t221 - t228;
t361 = -t275 * t333 + t303 * t337;
t360 = qJD(5) * t227 + t224;
t358 = -t276 + (-t371 + t392) * pkin(4);
t310 = -t380 - t431;
t356 = -t310 - t380;
t315 = t433 * t337;
t355 = qJD(2) * t348 + qJD(5) * t315 + t337 * t382 + t361;
t314 = t433 * t333;
t354 = pkin(9) * t371 - qJD(5) * t314 - t333 * t382 - t405;
t223 = t227 * t332 + t429;
t297 = t337 * t311;
t251 = -pkin(9) * t414 + t297 + (-pkin(4) - t432) * t338;
t259 = -pkin(9) * t416 + t402;
t353 = t251 * t332 + t259 * t336;
t288 = t338 * t420 + t418;
t257 = -t288 * t333 - t337 * t419;
t347 = -t288 * t337 + t333 * t419;
t352 = t257 * t336 + t332 * t347;
t351 = t257 * t332 - t336 * t347;
t349 = qJD(2) * t328 - t322 * t338;
t287 = -t331 * t338 + t334 * t420;
t225 = -t300 * t390 + t383;
t343 = qJD(3) * (-t356 - t431);
t341 = qJD(2) ^ 2;
t326 = -pkin(4) * t337 - pkin(3);
t307 = (pkin(4) * t333 + pkin(7)) * t334;
t282 = t301 * t334;
t277 = pkin(4) * t344 + pkin(7) * t393;
t256 = qJD(3) * t288 + t359;
t255 = -qJD(3) * t287 + t338 * t377;
t237 = -t390 * t416 + (t414 * t434 + t376) * t336 + t435 * t332;
t236 = -qJD(3) * t346 - t281 * t434;
t235 = pkin(4) * t271 + t245;
t230 = qJD(4) * t257 + t255 * t337 + t333 * t378;
t229 = qJD(4) * t347 - t255 * t333 + t337 * t378;
t222 = -t232 * t332 + t430;
t1 = [(-t229 * t322 + t256 * t298 + t271 * t287) * MDP(17) + (t230 * t322 + t256 * t300 + t270 * t287) * MDP(18) + (t256 * t247 + t287 * t226 - (-qJD(5) * t351 + t229 * t336 - t230 * t332) * t318) * MDP(24) + (-t256 * t350 + t287 * t225 + (qJD(5) * t352 + t229 * t332 + t230 * t336) * t318) * MDP(25) + (-t256 * MDP(10) - t255 * MDP(11) + (MDP(17) * t257 + MDP(18) * t347 + MDP(24) * t352 - MDP(25) * t351) * t397) * qJD(3) + ((-MDP(10) * t334 - MDP(11) * t338) * t339 * t385 + (-t339 * MDP(4) + (-MDP(10) * t338 + MDP(11) * t334 - MDP(3)) * t335) * t341) * t330; 0.2e1 * t369 * t439 - 0.2e1 * t385 * t438 + MDP(7) * t410 - MDP(8) * t413 + (-pkin(7) * t410 + t334 * t343) * MDP(10) + (pkin(7) * t413 + t338 * t343) * MDP(11) + (t270 * t414 + t300 * t435) * MDP(12) + ((-t298 * t337 - t300 * t333) * t393 + (-t425 - t271 * t337 + (t298 * t333 - t300 * t337) * qJD(4)) * t334) * MDP(13) + (t322 * t373 - t270 * t338 + (t300 * t334 + t337 * t349) * qJD(3)) * MDP(14) + (t322 * t372 + t271 * t338 + (-t298 * t334 - t333 * t349) * qJD(3)) * MDP(15) + (-t322 - t396) * t366 + ((t311 * t392 + t437) * t322 + ((pkin(7) * t298 + t426) * qJD(3) + (t417 + (pkin(7) * t322 + t268) * t337) * qJD(4) + t363) * t338 + (-t298 * t380 + t337 * t387 + pkin(7) * t271 + t428 + ((-pkin(7) * t415 + t297) * qJD(2) + t239) * qJD(3)) * t334) * MDP(17) + (t436 * t322 + (t267 * t386 + (qJD(3) * t300 - t374) * pkin(7) + t345) * t338 + (-t300 * t380 - t333 * t387 + pkin(7) * t270 + t427 + (-pkin(7) * t421 - qJD(2) * t402 - t240) * qJD(3)) * t334) * MDP(18) + (-t225 * t282 - t236 * t350) * MDP(19) + (-t225 * t281 + t226 * t282 - t236 * t247 + t237 * t350) * MDP(20) + (-t225 * t338 - t236 * t318 + (-qJD(2) * t282 - t350) * t394) * MDP(21) + (t226 * t338 + t237 * t318 + (-qJD(2) * t281 - t247) * t394) * MDP(22) + (-t318 - t396) * t367 + (t277 * t247 + t307 * t226 + t235 * t281 + t243 * t237 - t365 * t338 + (t332 * t408 + t336 * t409) * t318 + (t223 * t338 + t318 * t353) * qJD(5) + (-t247 * t380 + ((t251 * t336 - t259 * t332) * qJD(2) + t222) * qJD(3)) * t334) * MDP(24) + (-t277 * t350 + t307 * t225 - t235 * t282 + t243 * t236 + (t360 * t336 + t364) * t338 + ((qJD(5) * t251 + t408) * t336 + (-qJD(5) * t259 - t409) * t332) * t318 + (t350 * t380 + (-qJD(2) * t353 - t223) * qJD(3)) * t334) * MDP(25); (qJD(3) * t276 - t245) * MDP(10) + t356 * t396 * MDP(11) + (-t300 * t421 + t425) * MDP(12) + ((t270 + t424) * t337 + (-t271 + t423) * t333) * MDP(13) + (-t322 * t391 + (t322 * t412 + (-t300 + t395) * t334) * qJD(2)) * MDP(14) + (t374 + (-t322 * t415 + (t298 + t386) * t334) * qJD(2)) * MDP(15) + (-pkin(3) * t271 - t427 + t361 * t322 - t276 * t298 + (pkin(8) * t421 + t426) * qJD(4) + (-t239 * t334 + (-pkin(8) * t394 - t267 * t338) * t333) * qJD(2)) * MDP(17) + (-pkin(3) * t270 + t428 - t405 * t322 - t276 * t300 + (-pkin(8) * t322 * t333 + t267 * t337) * qJD(4) + (-t267 * t412 + (-pkin(8) * t386 + t240) * t334) * qJD(2)) * MDP(18) + (t225 * t302 - t350 * t407) * MDP(19) + (-t225 * t301 - t226 * t302 - t247 * t407 + t350 * t406) * MDP(20) + (t326 * t226 + t235 * t301 + t406 * t243 + t358 * t247) * MDP(24) + (t326 * t225 + t235 * t302 + t407 * t243 - t350 * t358) * MDP(25) + (-t407 * MDP(21) + t406 * MDP(22) + (t332 * t354 + t336 * t355) * MDP(24) + (-t332 * t355 + t336 * t354) * MDP(25)) * t318 + (-t310 * MDP(10) + t322 * MDP(16) + (qJD(3) * t302 + t350) * MDP(21) + (-qJD(3) * t301 + t247) * MDP(22) + t318 * MDP(23) + ((-t314 * t336 - t315 * t332) * qJD(3) - t222) * MDP(24) + (-(-t314 * t332 + t315 * t336) * qJD(3) + t223) * MDP(25)) * t397 + (-t338 * t439 + t438) * t341; t300 * t298 * MDP(12) + (-t298 ^ 2 + t300 ^ 2) * MDP(13) + (t270 - t424) * MDP(14) + (-t271 - t423) * MDP(15) + qJD(2) * t366 + (-t240 * t322 - t267 * t300 + t342) * MDP(17) + (-t239 * t322 + t267 * t298 - t345) * MDP(18) + (t225 - t443) * MDP(21) + (-t226 + t442) * MDP(22) + ((-t231 * t332 - t429) * t318 - t223 * qJD(5) + (-t247 * t300 + t318 * t390 + t336 * t370) * pkin(4) + t441) * MDP(24) + (t444 + t228 + (t232 * t318 - t221) * t332 + (-t231 * t318 - t360) * t336 + (t300 * t350 + t318 * t389 - t332 * t370) * pkin(4)) * MDP(25) + t440; (t383 - t443) * MDP(21) + (-t362 + t442) * MDP(22) + (-t223 * t318 + t441) * MDP(24) + (-t222 * t318 - t336 * t224 - t364 + t444) * MDP(25) + (-MDP(21) * t422 + MDP(22) * t350 - MDP(24) * t223 - MDP(25) * t430) * qJD(5) + t440;];
tauc = t1;
