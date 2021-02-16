% Calculate Coriolis joint torque vector for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:55
% EndTime: 2021-01-15 23:11:05
% DurationCPUTime: 3.09s
% Computational Cost: add. (4022->336), mult. (10573->456), div. (0->0), fcn. (7720->8), ass. (0->165)
t407 = (qJD(1) * qJD(2));
t448 = -2 * t407;
t359 = sin(qJ(2));
t447 = MDP(4) * t359;
t442 = pkin(6) + pkin(7);
t335 = t442 * t359;
t329 = qJD(1) * t335;
t439 = qJD(2) * pkin(2);
t323 = -t329 + t439;
t406 = qJD(2) * t442;
t390 = qJD(1) * t406;
t324 = t359 * t390;
t361 = cos(qJ(3));
t446 = (qJD(3) * t323 - t324) * t361;
t362 = cos(qJ(2));
t445 = (t359 ^ 2 - t362 ^ 2) * MDP(5);
t358 = sin(qJ(3));
t411 = qJD(1) * t362;
t412 = qJD(1) * t359;
t318 = -t358 * t411 - t361 * t412;
t313 = t318 * qJ(4);
t336 = t442 * t362;
t331 = qJD(1) * t336;
t319 = t358 * t331;
t399 = t361 * t323 - t319;
t276 = t313 + t399;
t328 = t358 * t362 + t359 * t361;
t444 = qJD(1) * t328;
t353 = qJD(2) + qJD(3);
t299 = t353 * t328;
t294 = t299 * qJD(1);
t404 = t361 * t411;
t405 = t358 * t412;
t316 = -t404 + t405;
t325 = t362 * t390;
t410 = qJD(3) * t358;
t397 = -t358 * t325 - t331 * t410;
t243 = -qJ(4) * t294 - qJD(4) * t316 + t397 + t446;
t356 = sin(pkin(9));
t403 = t362 * t407;
t293 = qJD(3) * t404 - t353 * t405 + t361 * t403;
t420 = t361 * t331;
t383 = -t323 * t358 - t420;
t398 = t358 * t324 - t361 * t325;
t370 = t383 * qJD(3) + t398;
t367 = -qJ(4) * t293 + qJD(4) * t318 + t370;
t438 = cos(pkin(9));
t232 = t438 * t243 + t356 * t367;
t327 = t358 * t359 - t361 * t362;
t330 = t359 * t406;
t332 = t362 * t406;
t425 = t335 * t361;
t374 = -qJD(3) * t425 - t361 * t330 - t358 * t332 - t336 * t410;
t254 = -qJ(4) * t299 - qJD(4) * t327 + t374;
t298 = t353 * t327;
t381 = t335 * t358 - t336 * t361;
t369 = t381 * qJD(3) + t330 * t358 - t332 * t361;
t366 = qJ(4) * t298 - qJD(4) * t328 + t369;
t236 = t438 * t254 + t356 * t366;
t272 = pkin(3) * t353 + t276;
t437 = qJ(4) * t316;
t277 = -t383 - t437;
t424 = t356 * t277;
t248 = t438 * t272 - t424;
t246 = -t353 * pkin(4) - t248;
t290 = -t438 * t316 + t318 * t356;
t350 = -pkin(2) * t362 - pkin(1);
t334 = t350 * qJD(1);
t300 = pkin(3) * t316 + qJD(4) + t334;
t375 = -t356 * t316 - t438 * t318;
t253 = -pkin(4) * t290 - pkin(8) * t375 + t300;
t296 = t438 * t327 + t328 * t356;
t297 = -t356 * t327 + t438 * t328;
t302 = pkin(3) * t327 + t350;
t259 = pkin(4) * t296 - pkin(8) * t297 + t302;
t268 = -t438 * t298 - t356 * t299;
t285 = qJD(5) - t290;
t231 = t243 * t356 - t438 * t367;
t286 = -qJ(4) * t327 - t381;
t372 = -qJ(4) * t328 - t336 * t358 - t425;
t261 = t438 * t286 + t356 * t372;
t265 = t293 * t356 + t438 * t294;
t388 = t231 * t297 - t261 * t265;
t443 = t246 * t268 - (qJD(5) * t259 + t236) * t285 - (qJD(5) * t253 + t232) * t296 + t388;
t441 = pkin(3) * t318;
t440 = pkin(2) * qJD(3);
t357 = sin(qJ(5));
t409 = qJD(5) * t357;
t266 = t438 * t293 - t356 * t294;
t360 = cos(qJ(5));
t408 = qJD(5) * t360;
t417 = t360 * t266 + t353 * t408;
t240 = -t375 * t409 + t417;
t436 = t240 * t357;
t435 = t246 * t290;
t434 = t246 * t297;
t433 = t259 * t265;
t432 = t266 * t357;
t427 = t375 * t357;
t279 = -t360 * t353 + t427;
t431 = t279 * t285;
t430 = t279 * t375;
t281 = t353 * t357 + t360 * t375;
t429 = t281 * t285;
t428 = t281 * t375;
t426 = t334 * t318;
t423 = t356 * t358;
t422 = t357 * t265;
t363 = qJD(2) ^ 2;
t421 = t359 * t363;
t262 = t360 * t265;
t419 = t362 * t363;
t364 = qJD(1) ^ 2;
t418 = t362 * t364;
t273 = t438 * t277;
t249 = t356 * t272 + t273;
t414 = -t361 * t329 - t319;
t282 = t313 + t414;
t382 = t329 * t358 - t420;
t373 = t382 + t437;
t401 = t438 * t358;
t416 = -t282 * t356 + t438 * t373 + (t356 * t361 + t401) * t440;
t415 = t438 * t282 + t356 * t373 - (t438 * t361 - t423) * t440;
t349 = pkin(2) * t361 + pkin(3);
t312 = pkin(2) * t401 + t356 * t349;
t352 = t359 * t439;
t351 = pkin(2) * t412;
t402 = -pkin(2) * t353 - t323;
t283 = pkin(3) * t294 + qJD(2) * t351;
t292 = pkin(3) * t299 + t352;
t400 = pkin(1) * t448;
t396 = t360 * t285;
t258 = pkin(4) * t375 - pkin(8) * t290 - t441;
t306 = pkin(8) + t312;
t392 = qJD(5) * t306 + t258 + t351;
t346 = pkin(3) * t356 + pkin(8);
t391 = qJD(5) * t346 + t258;
t247 = pkin(8) * t353 + t249;
t234 = t247 * t360 + t253 * t357;
t389 = t231 * t357 + t234 * t375 + t246 * t408;
t387 = -t265 * t306 - t435;
t386 = -t265 * t346 - t435;
t385 = t247 * t357 - t253 * t360;
t384 = t248 * t290 + t249 * t375;
t380 = t262 + (t290 * t357 - t409) * t285;
t379 = -t231 * t360 + t246 * t409 + t375 * t385;
t378 = -t300 * t375 - t231;
t377 = t334 * t316 - t397;
t376 = t268 * t360 - t297 * t409;
t311 = -pkin(2) * t423 + t438 * t349;
t241 = t281 * qJD(5) + t432;
t368 = -t318 * t316 * MDP(11) - t285 * t375 * MDP(26) + ((t240 - t431) * t360 + (-t241 - t429) * t357) * MDP(23) + (t380 + t430) * MDP(25) + (t285 * t396 + t422 - t428) * MDP(24) + (t281 * t396 + t436) * MDP(22) + t293 * MDP(13) + (-t316 ^ 2 + t318 ^ 2) * MDP(12) + (t316 * MDP(13) + (-t318 - t444) * MDP(14)) * t353;
t365 = -t300 * t290 - t232;
t347 = -t438 * pkin(3) - pkin(4);
t305 = -pkin(4) - t311;
t301 = t351 - t441;
t267 = -t298 * t356 + t438 * t299;
t260 = t286 * t356 - t438 * t372;
t251 = t438 * t276 - t424;
t250 = t276 * t356 + t273;
t239 = pkin(4) * t267 - pkin(8) * t268 + t292;
t238 = pkin(4) * t265 - pkin(8) * t266 + t283;
t237 = t360 * t238;
t235 = t254 * t356 - t438 * t366;
t1 = [0.2e1 * t403 * t447 + t445 * t448 + MDP(6) * t419 - MDP(7) * t421 + (-pkin(6) * t419 + t359 * t400) * MDP(9) + (pkin(6) * t421 + t362 * t400) * MDP(10) + (t293 * t328 + t298 * t318) * MDP(11) + (-t293 * t327 - t294 * t328 + t298 * t316 + t299 * t318) * MDP(12) + (t350 * t294 + t334 * t299 + (qJD(1) * t327 + t316) * t352) * MDP(16) + (t350 * t293 - t334 * t298 + (-t318 + t444) * t352) * MDP(17) + (t265 * t302 + t267 * t300 + t283 * t296 - t290 * t292) * MDP(18) + (t266 * t302 + t268 * t300 + t283 * t297 + t292 * t375) * MDP(19) + (-t232 * t296 + t235 * t375 + t236 * t290 - t248 * t268 - t249 * t267 + t260 * t266 + t388) * MDP(20) + (t231 * t260 + t232 * t261 - t235 * t248 + t236 * t249 + t283 * t302 + t292 * t300) * MDP(21) + (t240 * t297 * t360 + t376 * t281) * MDP(22) + ((-t279 * t360 - t281 * t357) * t268 + (-t436 - t241 * t360 + (t279 * t357 - t281 * t360) * qJD(5)) * t297) * MDP(23) + (t240 * t296 + t297 * t262 + t267 * t281 + t376 * t285) * MDP(24) + (-t297 * t422 - t241 * t296 - t267 * t279 + (-t268 * t357 - t297 * t408) * t285) * MDP(25) + (t265 * t296 + t267 * t285) * MDP(26) + (-t385 * t267 + t235 * t279 + t237 * t296 + t260 * t241 + (t239 * t285 + t433 + (-t247 * t296 - t261 * t285 + t434) * qJD(5)) * t360 + t443 * t357) * MDP(27) + (-t234 * t267 + t235 * t281 + t260 * t240 + (-(-qJD(5) * t261 + t239) * t285 - t433 - (-qJD(5) * t247 + t238) * t296 - qJD(5) * t434) * t357 + t443 * t360) * MDP(28) + (-t298 * MDP(13) - t299 * MDP(14) + t369 * MDP(16) - t374 * MDP(17) - t235 * MDP(18) - t236 * MDP(19)) * t353; (-t265 * t312 - t266 * t311 - t290 * t415 + t375 * t416 + t384) * MDP(20) + (-t231 * t311 + t232 * t312 - t416 * t248 - t415 * t249 - t300 * t301) * MDP(21) + t368 + (t305 * t240 + t387 * t360 + t416 * t281 + (t392 * t357 + t415 * t360) * t285 + t389) * MDP(28) + (-t301 * t375 + t415 * t353 + t365) * MDP(19) + (-t316 * t351 + t426 - t382 * t353 + (t402 * t358 - t420) * qJD(3) + t398) * MDP(16) + (t305 * t241 + t387 * t357 + t416 * t279 + (t415 * t357 - t392 * t360) * t285 + t379) * MDP(27) + (t318 * t351 + t414 * t353 + (t402 * qJD(3) + t324) * t361 + t377) * MDP(17) + (t290 * t301 - t416 * t353 + t378) * MDP(18) + t364 * t445 - t418 * t447 + (t364 * t359 * MDP(9) + MDP(10) * t418) * pkin(1); t368 + (t347 * t241 - t250 * t279 + t386 * t357 + (t357 * t251 - t391 * t360) * t285 + t379) * MDP(27) + (t347 * t240 - t250 * t281 + t386 * t360 + (t360 * t251 + t391 * t357) * t285 + t389) * MDP(28) + (t251 * t353 + t375 * t441 + t365) * MDP(19) + (-t250 * t375 - t251 * t290 + (-t265 * t356 - t438 * t266) * pkin(3) + t384) * MDP(20) + (-t383 * t353 + t370 + t426) * MDP(16) + (t248 * t250 - t249 * t251 + (-t438 * t231 + t232 * t356 + t300 * t318) * pkin(3)) * MDP(21) + (t250 * t353 - t290 * t441 + t378) * MDP(18) + (t399 * t353 + t377 - t446) * MDP(17); (t353 * t375 + t265) * MDP(18) + (t290 * t353 + t266) * MDP(19) + (-t290 ^ 2 - t375 ^ 2) * MDP(20) + (t248 * t375 - t249 * t290 + t283) * MDP(21) + (t380 - t430) * MDP(27) + (-t285 ^ 2 * t360 - t422 - t428) * MDP(28); t281 * t279 * MDP(22) + (-t279 ^ 2 + t281 ^ 2) * MDP(23) + (t417 + t431) * MDP(24) + (t429 - t432) * MDP(25) + t265 * MDP(26) + (-t232 * t357 + t234 * t285 - t246 * t281 + t237) * MDP(27) + (-t232 * t360 - t238 * t357 + t246 * t279 - t285 * t385) * MDP(28) + (-MDP(24) * t427 - t281 * MDP(25) - t234 * MDP(27) + t385 * MDP(28)) * qJD(5);];
tauc = t1;
