% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:47
% EndTime: 2019-12-31 20:26:56
% DurationCPUTime: 4.24s
% Computational Cost: add. (4117->399), mult. (12620->567), div. (0->0), fcn. (10113->10), ass. (0->179)
t361 = sin(pkin(5));
t362 = cos(pkin(10));
t369 = cos(qJ(2));
t436 = t362 * t369;
t413 = t361 * t436;
t347 = qJD(1) * t413;
t360 = sin(pkin(10));
t366 = sin(qJ(2));
t427 = qJD(1) * t361;
t411 = t366 * t427;
t335 = -t360 * t411 + t347;
t332 = qJD(4) - t335;
t385 = t360 * t369 + t362 * t366;
t338 = t385 * t427;
t368 = cos(qJ(4));
t363 = cos(pkin(5));
t426 = qJD(1) * t363;
t399 = qJD(2) + t426;
t348 = t368 * t399;
t365 = sin(qJ(4));
t305 = t338 * t365 - t348;
t304 = qJD(5) + t305;
t425 = qJD(2) * t361;
t410 = t366 * t425;
t461 = pkin(2) * t410;
t460 = MDP(5) * (t366 ^ 2 - t369 ^ 2);
t459 = MDP(6) * t369;
t452 = pkin(7) + qJ(3);
t408 = t452 * t366;
t397 = t361 * t408;
t453 = pkin(1) * t369;
t322 = (pkin(2) + t453) * t363 - t397;
t454 = pkin(1) * t366;
t415 = t363 * t454;
t438 = t361 * t369;
t333 = t438 * t452 + t415;
t290 = t360 * t322 + t362 * t333;
t281 = pkin(8) * t363 + t290;
t439 = t361 * t366;
t340 = t360 * t439 - t413;
t341 = t385 * t361;
t394 = (-pkin(2) * t369 - pkin(1)) * t361;
t296 = pkin(3) * t340 - pkin(8) * t341 + t394;
t458 = t368 * t281 + t365 * t296;
t357 = t361 ^ 2;
t381 = -pkin(7) * t438 - t415;
t457 = -t357 * t454 + t381 * t363;
t336 = qJD(2) * t341;
t330 = qJD(1) * t336;
t417 = t363 * t453;
t352 = qJD(1) * t417;
t308 = qJD(2) * pkin(2) + t352 + (t363 * pkin(2) - t397) * qJD(1);
t324 = t333 * qJD(1);
t437 = t362 * t324;
t267 = t360 * t308 + t437;
t260 = pkin(8) * t399 + t267;
t384 = qJD(1) * t394;
t342 = qJD(3) + t384;
t278 = -pkin(3) * t335 - pkin(8) * t338 + t342;
t246 = t260 * t368 + t278 * t365;
t350 = t352 * qJD(2);
t372 = (-qJD(2) * t408 + qJD(3) * t369) * t361;
t301 = qJD(1) * t372 + t350;
t310 = -qJD(2) * t333 - qJD(3) * t439;
t371 = qJD(1) * t310;
t257 = t362 * t301 + t360 * t371;
t395 = qJD(1) * t410;
t331 = qJD(2) * t347 - t360 * t395;
t349 = pkin(2) * t395;
t279 = pkin(3) * t330 - pkin(8) * t331 + t349;
t404 = t257 * t365 - t368 * t279;
t455 = -qJD(4) * t246 - t404;
t234 = -pkin(4) * t330 - t455;
t375 = -t368 * t338 - t365 * t399;
t456 = t304 * (-pkin(4) * t375 + pkin(9) * t304) + t234;
t270 = -qJD(4) * t375 + t365 * t331;
t370 = qJD(1) ^ 2;
t423 = qJD(4) * t365;
t269 = qJD(4) * t348 + t368 * t331 - t338 * t423;
t364 = sin(qJ(5));
t367 = cos(qJ(5));
t419 = qJD(5) * t367;
t412 = t367 * t269 + t364 * t330 + t332 * t419;
t420 = qJD(5) * t364;
t242 = t375 * t420 + t412;
t451 = t242 * t364;
t443 = t375 * t364;
t271 = -t367 * t332 - t443;
t450 = t271 * t304;
t273 = t332 * t364 - t367 * t375;
t449 = t273 * t304;
t448 = t273 * t335;
t401 = t304 * t367;
t447 = t305 * t332;
t446 = t305 * t338;
t445 = t375 * t332;
t444 = t375 * t338;
t442 = t332 * t365;
t441 = t335 * t368;
t440 = t357 * t370;
t312 = t360 * t324;
t435 = t364 * t270;
t433 = t367 * t270;
t323 = -qJD(1) * t397 + t352;
t282 = t323 * t360 + t437;
t432 = t282 - t332 * (pkin(4) * t365 - pkin(9) * t368);
t283 = t323 * t362 - t312;
t292 = pkin(2) * t411 + pkin(3) * t338 - pkin(8) * t335;
t430 = t368 * t283 + t365 * t292;
t429 = t368 * t330 + t335 * t442;
t424 = qJD(4) * t364;
t422 = qJD(4) * t367;
t421 = qJD(4) * t368;
t418 = qJD(1) * qJD(2);
t414 = t369 * t440;
t356 = -pkin(2) * t362 - pkin(3);
t409 = t361 * t363 * t370;
t407 = t357 * t418;
t380 = t368 * t257 - t260 * t423 + t278 * t421 + t365 * t279;
t233 = pkin(9) * t330 + t380;
t256 = t301 * t360 - t362 * t371;
t238 = pkin(4) * t270 - pkin(9) * t269 + t256;
t405 = -t233 * t364 + t367 * t238;
t403 = t269 * t364 - t367 * t330;
t353 = qJD(2) * t417;
t309 = t353 + t372;
t261 = t309 * t360 - t362 * t310;
t266 = t362 * t308 - t312;
t289 = t322 * t362 - t360 * t333;
t402 = t368 * t332;
t343 = -pkin(4) * t368 - pkin(9) * t365 + t356;
t400 = pkin(9) * t338 - qJD(5) * t343 + t430;
t396 = t369 * t407;
t392 = t233 * t367 + t238 * t364;
t241 = pkin(9) * t332 + t246;
t259 = -pkin(3) * t399 - t266;
t244 = t305 * pkin(4) + pkin(9) * t375 + t259;
t232 = t241 * t367 + t244 * t364;
t391 = t241 * t364 - t244 * t367;
t252 = pkin(9) * t340 + t458;
t280 = -pkin(3) * t363 - t289;
t315 = t341 * t365 - t363 * t368;
t316 = t341 * t368 + t363 * t365;
t253 = pkin(4) * t315 - pkin(9) * t316 + t280;
t390 = t252 * t367 + t253 * t364;
t389 = -t252 * t364 + t253 * t367;
t245 = -t260 * t365 + t278 * t368;
t262 = t309 * t362 + t310 * t360;
t337 = (-t360 * t366 + t436) * t425;
t293 = pkin(3) * t336 - pkin(8) * t337 + t461;
t388 = -t262 * t365 + t293 * t368;
t386 = -t281 * t365 + t296 * t368;
t295 = t316 * t367 + t340 * t364;
t294 = t316 * t364 - t340 * t367;
t383 = -t304 * t419 - t435;
t382 = t304 * t420 - t433;
t379 = t368 * t262 - t281 * t423 + t365 * t293 + t296 * t421;
t355 = pkin(2) * t360 + pkin(8);
t378 = t259 * t332 - t355 * t330;
t298 = t338 * t364 + t367 * t441;
t376 = -t420 * t365 + t367 * t421 - t298;
t374 = qJD(4) * t271 + t383;
t240 = -pkin(4) * t332 - t245;
t373 = -pkin(9) * t270 + (t240 + t245) * t304;
t297 = -t367 * t338 + t364 * t441;
t288 = -qJD(4) * t315 + t337 * t368;
t287 = qJD(4) * t316 + t337 * t365;
t264 = t273 * t423;
t251 = -pkin(4) * t340 - t386;
t250 = -qJD(5) * t294 + t288 * t367 + t336 * t364;
t249 = qJD(5) * t295 + t288 * t364 - t336 * t367;
t247 = -pkin(4) * t338 + t283 * t365 - t292 * t368;
t243 = t273 * qJD(5) + t403;
t239 = pkin(4) * t287 - pkin(9) * t288 + t261;
t236 = -pkin(4) * t336 + qJD(4) * t458 - t388;
t235 = pkin(9) * t336 + t379;
t230 = -qJD(5) * t232 + t405;
t229 = -qJD(5) * t391 + t392;
t1 = [0.2e1 * t366 * MDP(4) * t396 - 0.2e1 * t407 * t460 + (t381 * qJD(2) ^ 2 + 0.2e1 * t418 * t457) * MDP(9) + (-0.2e1 * pkin(1) * t396 - (-pkin(7) * t410 + t353) * t399 - (-pkin(7) * t395 + t350) * t363) * MDP(10) + (t256 * t341 - t257 * t340 + t261 * t338 + t262 * t335 - t266 * t337 - t267 * t336 - t289 * t331 - t290 * t330) * MDP(11) + (-t256 * t289 + t257 * t290 - t266 * t261 + t267 * t262 + (t342 + t384) * t461) * MDP(12) + (t269 * t316 - t288 * t375) * MDP(13) + (-t269 * t315 - t270 * t316 + t287 * t375 - t288 * t305) * MDP(14) + (t269 * t340 + t288 * t332 + t316 * t330 - t336 * t375) * MDP(15) + (-t270 * t340 - t287 * t332 - t305 * t336 - t315 * t330) * MDP(16) + (t330 * t340 + t332 * t336) * MDP(17) + (t388 * t332 + t386 * t330 - t404 * t340 + t245 * t336 + t261 * t305 + t280 * t270 + t256 * t315 + t259 * t287 + (-t246 * t340 - t332 * t458) * qJD(4)) * MDP(18) + (-t246 * t336 + t256 * t316 + t259 * t288 - t261 * t375 + t280 * t269 - t330 * t458 - t332 * t379 - t340 * t380) * MDP(19) + (t242 * t295 + t250 * t273) * MDP(20) + (-t242 * t294 - t243 * t295 - t249 * t273 - t250 * t271) * MDP(21) + (t242 * t315 + t250 * t304 + t270 * t295 + t273 * t287) * MDP(22) + (-t243 * t315 - t249 * t304 - t270 * t294 - t271 * t287) * MDP(23) + (t270 * t315 + t287 * t304) * MDP(24) + ((-qJD(5) * t390 - t235 * t364 + t239 * t367) * t304 + t389 * t270 + t230 * t315 - t391 * t287 + t236 * t271 + t251 * t243 + t234 * t294 + t240 * t249) * MDP(25) + (-(qJD(5) * t389 + t235 * t367 + t239 * t364) * t304 - t390 * t270 - t229 * t315 - t232 * t287 + t236 * t273 + t251 * t242 + t234 * t295 + t240 * t250) * MDP(26) + (-MDP(7) * t410 + t425 * t459) * (qJD(2) + 0.2e1 * t426); t440 * t460 - t409 * t459 + (pkin(1) * t414 + (-pkin(7) * t411 + t352) * t426) * MDP(10) + ((t267 - t282) * t338 + (t266 - t283) * t335 + (-t330 * t360 - t331 * t362) * pkin(2)) * MDP(11) + (t266 * t282 - t267 * t283 + (-t256 * t362 + t257 * t360 - t342 * t411) * pkin(2)) * MDP(12) + (t269 * t365 - t375 * t402) * MDP(13) + ((t269 - t447) * t368 + (-t270 + t445) * t365) * MDP(14) + (t365 * t330 + t332 * t402 + t444) * MDP(15) + (-t332 * t423 + t429 + t446) * MDP(16) - t332 * t338 * MDP(17) + (-t245 * t338 + t356 * t270 - t282 * t305 + (-t256 + (-qJD(4) * t355 - t292) * t332) * t368 + (t283 * t332 + t378) * t365) * MDP(18) + (t246 * t338 + t256 * t365 + t356 * t269 + t282 * t375 + (t355 * t423 + t430) * t332 + t378 * t368) * MDP(19) + (t242 * t365 * t367 + t273 * t376) * MDP(20) + (t271 * t298 + t273 * t297 + (-t271 * t367 - t273 * t364) * t421 + (-t451 - t243 * t367 + (t271 * t364 - t273 * t367) * qJD(5)) * t365) * MDP(21) + (-t242 * t368 + t264 + (t433 - t448) * t365 + t376 * t304) * MDP(22) + (t243 * t368 + (-t364 * t421 + t297) * t304 + (-t271 * t332 + t383) * t365) * MDP(23) + (-t270 * t368 + t304 * t442) * MDP(24) + (t343 * t433 - t240 * t297 - t247 * t271 + (t364 * t400 - t367 * t432) * t304 + (t240 * t424 + t355 * t374 - t230) * t368 + (t240 * t419 + t391 * t335 + t234 * t364 + t355 * t243 + (t304 * t355 * t364 - t391) * qJD(4)) * t365) * MDP(25) + (-t343 * t435 - t240 * t298 - t247 * t273 + (t364 * t432 + t367 * t400) * t304 + (t240 * t422 + t229 + (qJD(4) * t273 + t382) * t355) * t368 + (-t240 * t420 + t232 * t335 + t234 * t367 + t242 * t355 + (t355 * t401 - t232) * qJD(4)) * t365) * MDP(26) - t457 * MDP(9) * t370 + (-MDP(4) * t414 + MDP(7) * t409) * t366; (-t335 ^ 2 - t338 ^ 2) * MDP(11) + (t266 * t338 - t267 * t335 + t349) * MDP(12) + (t429 - t446) * MDP(18) + MDP(19) * t444 + t297 * t304 * MDP(25) + (t298 * t304 + t264) * MDP(26) + ((-t304 * t424 - t243) * MDP(25) + (-t304 * t422 - t242) * MDP(26) - t332 ^ 2 * MDP(19)) * t368 + (-qJD(4) * t332 * MDP(18) - t330 * MDP(19) + (-t271 * t335 + t374) * MDP(25) + (t382 - t448) * MDP(26)) * t365; -t305 ^ 2 * MDP(14) + (t269 + t447) * MDP(15) + (-t270 - t445) * MDP(16) + t330 * MDP(17) + (t246 * t332 + t455) * MDP(18) + (t245 * t332 + t259 * t305 - t380) * MDP(19) + (t273 * t401 + t451) * MDP(20) + ((t242 - t450) * t367 + (-t243 - t449) * t364) * MDP(21) + (t304 * t401 + t435) * MDP(22) + (-t304 ^ 2 * t364 + t433) * MDP(23) + (-pkin(4) * t243 - t246 * t271 + t373 * t364 - t367 * t456) * MDP(25) + (-pkin(4) * t242 - t246 * t273 + t364 * t456 + t373 * t367) * MDP(26) - (MDP(13) * t305 - MDP(14) * t375 - t259 * MDP(18) - MDP(22) * t273 + t271 * MDP(23) - t304 * MDP(24) + MDP(25) * t391 + t232 * MDP(26)) * t375; t273 * t271 * MDP(20) + (-t271 ^ 2 + t273 ^ 2) * MDP(21) + (t412 + t450) * MDP(22) + (-t403 + t449) * MDP(23) + t270 * MDP(24) + (t232 * t304 - t240 * t273 + t405) * MDP(25) + (t240 * t271 - t304 * t391 - t392) * MDP(26) + (MDP(22) * t443 - MDP(23) * t273 - MDP(25) * t232 + MDP(26) * t391) * qJD(5);];
tauc = t1;
