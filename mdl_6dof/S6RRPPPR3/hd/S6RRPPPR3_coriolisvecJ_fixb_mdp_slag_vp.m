% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:44
% EndTime: 2019-03-09 08:15:51
% DurationCPUTime: 3.47s
% Computational Cost: add. (2063->407), mult. (4610->538), div. (0->0), fcn. (2645->6), ass. (0->182)
t354 = sin(pkin(9));
t362 = cos(qJ(2));
t421 = qJD(1) * t362;
t401 = t354 * t421;
t355 = cos(pkin(9));
t411 = t355 * qJD(2);
t309 = t401 - t411;
t360 = sin(qJ(2));
t407 = qJD(1) * qJD(2);
t399 = t360 * t407;
t448 = qJD(6) * t309 + t355 * t399;
t447 = -0.2e1 * t407;
t446 = t360 * MDP(4);
t352 = t360 ^ 2;
t353 = t362 ^ 2;
t445 = (t352 - t353) * MDP(5);
t339 = pkin(7) * t421;
t402 = qJ(4) * t421;
t318 = t339 - t402;
t351 = qJD(2) * qJ(3);
t306 = -t318 - t351;
t363 = -pkin(2) - pkin(3);
t403 = qJD(2) * t363;
t359 = sin(qJ(6));
t361 = cos(qJ(6));
t311 = t354 * t359 - t361 * t355;
t303 = t311 * qJD(6);
t415 = qJD(4) * t362;
t418 = qJD(2) * t360;
t376 = pkin(7) * t418 + t415;
t328 = qJ(4) * t399;
t350 = qJD(2) * qJD(3);
t426 = t328 + t350;
t278 = t376 * qJD(1) - t426;
t444 = pkin(7) - qJ(4);
t349 = -qJ(5) + t363;
t443 = pkin(8) - t349;
t358 = qJ(3) + pkin(4);
t442 = qJD(2) * pkin(2);
t422 = qJD(1) * t360;
t307 = -qJD(1) * pkin(1) - pkin(2) * t421 - qJ(3) * t422;
t284 = pkin(3) * t421 + qJD(4) - t307;
t386 = pkin(4) * t360 + qJ(5) * t362;
t268 = t386 * qJD(1) + t284;
t334 = qJ(4) * t422;
t405 = pkin(7) * t422;
t387 = qJD(3) + t405;
t378 = -t334 + t387;
t282 = t349 * qJD(2) + t378;
t239 = t355 * t268 - t282 * t354;
t441 = t239 * t362;
t240 = t354 * t268 + t355 * t282;
t440 = t240 * t362;
t345 = t362 * qJ(4);
t325 = pkin(7) * t362 - t345;
t439 = t278 * t325;
t365 = qJD(1) ^ 2;
t438 = t352 * t365;
t437 = t354 * t360;
t436 = t354 * t362;
t435 = t355 * t360;
t434 = t355 * t362;
t364 = qJD(2) ^ 2;
t433 = t360 * t364;
t432 = t362 * t364;
t431 = t362 * t365;
t372 = pkin(4) * t362 + t349 * t360;
t367 = t372 * qJD(2) + qJD(5) * t362;
t343 = t360 * qJD(3);
t398 = t362 * t407;
t425 = qJ(3) * t398 + qJD(1) * t343;
t252 = t367 * qJD(1) + t425;
t330 = pkin(7) * t398;
t416 = qJD(4) * t360;
t417 = qJD(2) * t362;
t287 = t330 + (-qJ(4) * t417 - t416) * qJD(1);
t279 = -qJD(2) * qJD(5) + t287;
t235 = t354 * t252 + t355 * t279;
t312 = t354 * t361 + t355 * t359;
t331 = qJD(6) + t422;
t430 = t303 * t331 - t312 * t398;
t424 = qJ(3) * t417 + t343;
t260 = t367 + t424;
t299 = t444 * t417 - t416;
t243 = t354 * t260 + t355 * t299;
t429 = t448 * t361;
t373 = t311 * t360;
t286 = qJD(1) * t373;
t428 = t303 + t286;
t285 = t312 * t422;
t304 = t312 * qJD(6);
t427 = -t304 - t285;
t336 = qJ(3) * t421;
t271 = t372 * qJD(1) + t336;
t251 = t354 * t271 + t355 * t318;
t322 = -t362 * pkin(2) - t360 * qJ(3) - pkin(1);
t308 = t362 * pkin(3) - t322;
t281 = t386 + t308;
t324 = t444 * t360;
t259 = t354 * t281 + t355 * t324;
t420 = qJD(2) * t311;
t419 = qJD(2) * t354;
t236 = pkin(8) * t309 + t240;
t414 = qJD(6) * t236;
t291 = qJD(2) * pkin(4) + qJD(5) - t306;
t412 = t291 * MDP(22);
t404 = pkin(5) * t354 - pkin(7);
t388 = t404 * t360;
t410 = -qJD(1) * t388 + qJD(3) - t334;
t316 = -t334 + t405;
t409 = qJD(3) + t316;
t408 = -qJD(4) - t284;
t406 = pkin(8) * t437;
t400 = t355 * t421;
t321 = t387 - t442;
t397 = -t321 - t442;
t396 = pkin(1) * t447;
t234 = t355 * t252 - t279 * t354;
t242 = t355 * t260 - t299 * t354;
t250 = t355 * t271 - t318 * t354;
t258 = t355 * t281 - t324 * t354;
t391 = t360 * t403;
t273 = qJD(1) * t391 + t425;
t283 = t391 + t424;
t395 = qJD(1) * t283 + t273;
t394 = qJD(1) * t308 + t284;
t393 = qJD(1) * t322 + t307;
t389 = t354 * t399;
t232 = -pkin(8) * t389 + t235;
t310 = t400 + t419;
t233 = pkin(5) * t422 + pkin(8) * t310 + t239;
t392 = -qJD(6) * t233 - t232;
t228 = t233 * t361 - t236 * t359;
t229 = t233 * t359 + t236 * t361;
t385 = t234 * t355 + t235 * t354;
t384 = -t234 * t354 + t235 * t355;
t383 = -t239 * t355 - t240 * t354;
t382 = t239 * t354 - t240 * t355;
t249 = pkin(5) * t360 + pkin(8) * t434 + t258;
t253 = pkin(8) * t436 + t259;
t381 = t249 * t361 - t253 * t359;
t380 = t249 * t359 + t253 * t361;
t267 = t309 * t359 - t310 * t361;
t379 = pkin(5) * t362 - pkin(8) * t435;
t288 = pkin(2) * t399 - t425;
t300 = pkin(2) * t418 - t424;
t377 = -pkin(7) * t364 - qJD(1) * t300 - t288;
t315 = t443 * t355;
t375 = -t379 * qJD(1) + qJD(5) * t354 + qJD(6) * t315 - t250;
t314 = t443 * t354;
t374 = -qJD(1) * t406 + qJD(5) * t355 - qJD(6) * t314 + t251;
t371 = t379 * qJD(2);
t370 = qJD(6) * t310 - t389;
t369 = qJD(2) * t388 - t415;
t368 = -t349 * t417 + (qJD(2) * t358 + qJD(5) - t291) * t360;
t320 = -pkin(7) * t399 + t350;
t323 = t339 + t351;
t366 = t320 * t362 + (t321 * t362 + (-t323 + t339) * t360) * qJD(2);
t245 = qJD(2) * t285 + t267 * qJD(6);
t346 = 0.2e1 * t350;
t337 = qJ(4) * t418;
t327 = pkin(5) * t355 + t358;
t317 = pkin(2) * t422 - t336;
t302 = -t404 * t362 - t345;
t301 = t361 * t309;
t298 = -t337 + t376;
t297 = t363 * t422 + t336;
t295 = t311 * t362;
t294 = t312 * t362;
t292 = t403 + t378;
t277 = t337 + t369;
t265 = -t310 * t359 - t301;
t264 = t369 * qJD(1) + t426;
t263 = -pkin(5) * t309 + t291;
t257 = -t362 * t303 - t312 * t418;
t256 = -qJD(2) * t373 + t362 * t304;
t244 = t370 * t359 + t429;
t238 = -qJD(2) * t406 + t243;
t237 = t242 + t371;
t231 = qJD(1) * t371 + t234;
t230 = t361 * t231;
t1 = [t445 * t447 + (t331 + t422) * MDP(27) * t417 + (t395 * t360 + (t394 * t362 - t298) * qJD(2)) * MDP(15) + (-t395 * t362 + (t394 * t360 + t299) * qJD(2)) * MDP(16) + (t278 * t362 - t287 * t360 + (-t292 * t362 - t306 * t360) * qJD(2) + (t298 * t362 - t299 * t360 + (-t324 * t362 + t325 * t360) * qJD(2)) * qJD(1)) * MDP(17) + (t234 * t258 + t235 * t259 + t239 * t242 + t240 * t243 - t291 * t298 - t439) * MDP(22) + (t273 * t308 + t283 * t284 + t287 * t324 + t292 * t299 + t298 * t306 - t439) * MDP(18) + (t278 * t434 + t298 * t310 + (-qJD(1) * t243 - t235) * t360 + (t291 * t435 - t440 + (-t259 * t362 + t325 * t435) * qJD(1)) * qJD(2)) * MDP(20) + (t278 * t436 + t298 * t309 + (qJD(1) * t242 + t234) * t360 + (t291 * t437 + t441 + (t258 * t362 + t325 * t437) * qJD(1)) * qJD(2)) * MDP(19) + (-pkin(7) * t432 + t360 * t396) * MDP(9) - MDP(7) * t433 + (pkin(7) * t433 + t362 * t396) * MDP(10) + (t377 * t362 + t393 * t418) * MDP(11) + (t244 * t360 + t256 * t331 + (qJD(1) * t295 + t267) * t417) * MDP(25) + (-t245 * t360 + t257 * t331 + (qJD(1) * t294 - t265) * t417) * MDP(26) + (-(t237 * t359 + t238 * t361) * t331 - (t231 * t359 + t232 * t361) * t360 + t277 * t267 + t302 * t244 + t264 * t295 + t263 * t256 + (-t228 * t360 - t381 * t331) * qJD(6) + (-t380 * qJD(1) - t229) * t417) * MDP(29) + ((t237 * t361 - t238 * t359) * t331 + (-t232 * t359 + t230) * t360 + t277 * t265 + t302 * t245 - t264 * t294 - t263 * t257 + (-t229 * t360 - t380 * t331) * qJD(6) + (t381 * qJD(1) + t228) * t417) * MDP(28) + (t377 * t360 - t393 * t417) * MDP(13) + (t242 * t310 + t243 * t309 + t385 * t362 + ((-t258 * t355 - t259 * t354) * qJD(1) + t383) * t418) * MDP(21) + (t244 * t295 + t256 * t267) * MDP(23) + (t244 * t294 - t245 * t295 - t256 * t265 + t257 * t267) * MDP(24) + t366 * MDP(12) + (t366 * pkin(7) + t288 * t322 + t300 * t307) * MDP(14) + 0.2e1 * t398 * t446 + MDP(6) * t432; -t431 * t446 + t365 * t445 + t346 * MDP(13) + (qJ(3) * t320 + qJD(3) * t323 - t307 * t317) * MDP(14) + (qJD(2) * t316 + t328 + t346) * MDP(15) + (-qJD(2) * t318 + t330) * MDP(16) + (-qJ(3) * t278 - t284 * t297 + t287 * t363 - t292 * t318 - t409 * t306) * MDP(18) + (-t278 * t355 - t409 * t309) * MDP(19) + (t278 * t354 - t409 * t310) * MDP(20) + (-t250 * t310 - t251 * t309 + (-qJD(5) * t309 + t239 * t422 - t235) * t355 + (qJD(5) * t310 + t240 * t422 + t234) * t354) * MDP(21) + (t382 * qJD(5) - t239 * t250 - t240 * t251 - t278 * t358 + t409 * t291 + t384 * t349) * MDP(22) + (-t244 * t312 + t428 * t267) * MDP(23) + (t244 * t311 + t245 * t312 - t428 * t265 - t427 * t267) * MDP(24) + (-t267 * t421 + t286 * t331 + t430) * MDP(25) + (-t427 * t331 + (t265 + t420) * t421) * MDP(26) - t331 * MDP(27) * t421 + (t327 * t245 - t264 * t311 + (t374 * t359 + t375 * t361) * t331 + t410 * t265 + t427 * t263 + ((t314 * t361 + t315 * t359) * qJD(2) - t228) * t421) * MDP(28) + (t327 * t244 - t264 * t312 + (-t375 * t359 + t374 * t361) * t331 + t410 * t267 + t428 * t263 + (-(t314 * t359 - t315 * t361) * qJD(2) + t229) * t421) * MDP(29) + (t365 * t360 * MDP(9) + MDP(10) * t431) * pkin(1) + ((-t307 * t360 + t317 * t362) * MDP(11) + ((t323 - t351) * t360 + (qJD(3) + t397) * t362) * MDP(12) + (t307 * t362 + t317 * t360) * MDP(13) + (t323 * t360 + t397 * t362) * pkin(7) * MDP(14) + (t408 * t362 + (-pkin(7) * qJD(2) - t297) * t360) * MDP(15) + ((-qJ(4) * qJD(2) + t297) * t362 + t408 * t360) * MDP(16) + (t292 - t409 - t403) * t362 * MDP(17) + (-t250 * t360 + t368 * t354 - t441) * MDP(19) + (t251 * t360 + t368 * t355 + t440) * MDP(20)) * qJD(1); t384 * MDP(22) + t430 * MDP(28) + (-MDP(19) * t355 + MDP(20) * t354) * t438 + (t286 * MDP(28) - t427 * MDP(29)) * t331 + (MDP(14) + MDP(18)) * t330 + (MDP(13) + MDP(15)) * (-t364 - t438) + (-t323 * MDP(14) + (t306 - t402) * MDP(18) + (t309 - t401) * MDP(19) + (t310 - t400) * MDP(20) - t412 - t265 * MDP(28) + (t311 * t421 - t267) * MDP(29)) * qJD(2) + ((-MDP(11) + MDP(16)) * t431 + (t307 * MDP(14) + t408 * MDP(18) + (-t309 * t354 - t310 * t355) * MDP(21) + t383 * MDP(22)) * qJD(1)) * t360; t425 * MDP(18) + t385 * MDP(22) + (-t353 * MDP(17) + (-MDP(19) * t354 - MDP(20) * t355 - MDP(17)) * t352) * t365 + (t427 * MDP(28) + t428 * MDP(29)) * t331 + ((t292 * MDP(18) + (t309 * t355 - t310 * t354) * MDP(21) - t382 * MDP(22) + (0.2e1 * MDP(16) + t363 * MDP(18) + (-t354 ^ 2 - t355 ^ 2) * MDP(21)) * qJD(2)) * t360 + (0.2e1 * qJD(2) * MDP(15) - t306 * MDP(18) + (-t309 + t411) * MDP(19) + (-t310 - t419) * MDP(20) + t412 + (t265 - t420) * MDP(28) + (-qJD(2) * t312 + t267) * MDP(29)) * t362) * qJD(1); (-t309 ^ 2 - t310 ^ 2) * MDP(21) + (-t239 * t310 - t240 * t309 - t278) * MDP(22) + (t267 * t331 + t245) * MDP(28) + (t301 * t331 + (-t389 + (qJD(6) + t331) * t310) * t359 + t429) * MDP(29) + ((-t310 + t419) * MDP(19) + (t309 + t411) * MDP(20)) * t422; -t265 ^ 2 * MDP(24) + (t265 * t331 + t429) * MDP(25) + MDP(27) * t398 + (t229 * t331 + t230) * MDP(28) + (t228 * t331 + t263 * t265) * MDP(29) + (MDP(23) * t265 + MDP(24) * t267 + MDP(26) * t331 - MDP(28) * t263) * t267 + (t370 * MDP(26) - MDP(28) * t414 + t392 * MDP(29)) * t361 + (t370 * MDP(25) - t448 * MDP(26) + t392 * MDP(28) + (-t231 + t414) * MDP(29)) * t359;];
tauc  = t1;
