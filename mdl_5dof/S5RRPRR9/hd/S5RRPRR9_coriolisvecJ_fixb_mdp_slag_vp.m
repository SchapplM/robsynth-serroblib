% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:22:06
% EndTime: 2019-12-31 20:22:13
% DurationCPUTime: 3.93s
% Computational Cost: add. (2951->352), mult. (7610->482), div. (0->0), fcn. (5685->8), ass. (0->167)
t358 = sin(pkin(9));
t362 = sin(qJ(2));
t399 = t362 * qJD(1);
t359 = cos(pkin(9));
t365 = cos(qJ(2));
t416 = t359 * t365;
t325 = qJD(1) * t416 - t358 * t399;
t437 = qJD(4) + qJD(5);
t447 = t325 - t437;
t337 = t358 * t365 + t359 * t362;
t327 = t337 * qJD(1);
t361 = sin(qJ(4));
t364 = cos(qJ(4));
t398 = t364 * qJD(2);
t305 = t327 * t361 - t398;
t363 = cos(qJ(5));
t307 = qJD(2) * t361 + t327 * t364;
t360 = sin(qJ(5));
t423 = t307 * t360;
t249 = t363 * t305 + t423;
t320 = qJD(4) - t325;
t316 = qJD(5) + t320;
t446 = t249 * t316;
t376 = t305 * t360 - t363 * t307;
t445 = t316 * t376;
t340 = t360 * t364 + t361 * t363;
t406 = t447 * t340;
t404 = qJD(4) * t361;
t420 = t325 * t361;
t444 = t404 - t420;
t394 = -pkin(2) * t365 - pkin(1);
t380 = t394 * qJD(1);
t343 = qJD(3) + t380;
t268 = -pkin(3) * t325 - pkin(7) * t327 + t343;
t434 = -qJ(3) - pkin(6);
t345 = t434 * t362;
t341 = qJD(1) * t345;
t433 = qJD(2) * pkin(2);
t333 = t341 + t433;
t346 = t434 * t365;
t342 = qJD(1) * t346;
t417 = t359 * t342;
t292 = t358 * t333 - t417;
t287 = qJD(2) * pkin(7) + t292;
t244 = t268 * t361 + t287 * t364;
t239 = -pkin(8) * t305 + t244;
t402 = qJD(5) * t360;
t236 = t239 * t402;
t330 = t358 * t342;
t291 = t333 * t359 + t330;
t286 = -qJD(2) * pkin(3) - t291;
t247 = pkin(4) * t305 + t286;
t443 = t247 * t249 + t236;
t397 = qJD(1) * qJD(2);
t390 = t365 * t397;
t391 = t362 * t397;
t318 = -t358 * t391 + t359 * t390;
t263 = qJD(4) * t398 + t364 * t318 - t327 * t404;
t326 = t337 * qJD(2);
t317 = qJD(1) * t326;
t349 = pkin(2) * t391;
t267 = pkin(3) * t317 - pkin(7) * t318 + t349;
t259 = t364 * t267;
t389 = qJD(2) * t434;
t323 = qJD(3) * t365 + t362 * t389;
t311 = t323 * qJD(1);
t324 = -t362 * qJD(3) + t365 * t389;
t369 = t324 * qJD(1);
t262 = t359 * t311 + t358 * t369;
t368 = -qJD(4) * t244 - t262 * t361 + t259;
t225 = pkin(4) * t317 - pkin(8) * t263 + t368;
t264 = qJD(4) * t307 + t318 * t361;
t403 = qJD(4) * t364;
t372 = t364 * t262 + t361 * t267 + t268 * t403 - t287 * t404;
t228 = -pkin(8) * t264 + t372;
t386 = t363 * t225 - t360 * t228;
t442 = t247 * t376 + t386;
t441 = t317 * MDP(24) + (-t249 ^ 2 + t376 ^ 2) * MDP(21) - t249 * t376 * MDP(20);
t440 = -0.2e1 * t397;
t439 = MDP(4) * t362;
t438 = MDP(5) * (t362 ^ 2 - t365 ^ 2);
t282 = t340 * t337;
t336 = t358 * t362 - t416;
t329 = t336 * qJD(2);
t412 = t364 * t329;
t373 = -t337 * t404 - t412;
t339 = t360 * t361 - t363 * t364;
t407 = t447 * t339;
t436 = -t316 * t407 - t317 * t340;
t385 = t263 * t360 + t363 * t264;
t232 = -qJD(5) * t376 + t385;
t351 = pkin(2) * t358 + pkin(7);
t435 = pkin(8) + t351;
t243 = t364 * t268 - t287 * t361;
t238 = -pkin(8) * t307 + t243;
t233 = pkin(4) * t320 + t238;
t432 = t233 * t363;
t431 = t239 * t363;
t430 = t249 * t327;
t429 = t376 * t327;
t428 = t263 * t361;
t427 = t305 * t320;
t426 = t305 * t327;
t425 = t307 * t320;
t424 = t307 * t327;
t419 = t337 * t361;
t418 = t337 * t364;
t415 = t361 * t317;
t414 = t361 * t329;
t366 = qJD(2) ^ 2;
t413 = t362 * t366;
t304 = t345 * t358 - t346 * t359;
t297 = t364 * t304;
t309 = t364 * t317;
t411 = t365 * t366;
t367 = qJD(1) ^ 2;
t410 = t365 * t367;
t279 = pkin(2) * t399 + pkin(3) * t327 - pkin(7) * t325;
t295 = t341 * t359 + t330;
t409 = t361 * t279 + t364 * t295;
t290 = pkin(3) * t336 - pkin(7) * t337 + t394;
t408 = t361 * t290 + t297;
t401 = qJD(5) * t363;
t396 = t362 * t433;
t395 = t363 * t263 - t360 * t264 - t305 * t401;
t352 = -pkin(2) * t359 - pkin(3);
t392 = t337 * t403;
t388 = qJD(4) * t435;
t387 = pkin(1) * t440;
t261 = t311 * t358 - t359 * t369;
t277 = t323 * t358 - t359 * t324;
t294 = t341 * t358 - t417;
t303 = -t359 * t345 - t346 * t358;
t384 = t320 * t364;
t383 = qJD(5) * t233 + t228;
t382 = t444 * pkin(4) - t294;
t381 = t406 * t316 - t339 * t317;
t272 = t364 * t279;
t335 = t435 * t364;
t379 = pkin(4) * t327 + qJD(5) * t335 - t295 * t361 + t272 + (-pkin(8) * t325 + t388) * t364;
t334 = t435 * t361;
t378 = -pkin(8) * t420 + qJD(5) * t334 + t361 * t388 + t409;
t227 = t233 * t360 + t431;
t377 = t261 * t337 - t304 * t317;
t375 = -t444 * t320 + t309;
t374 = t392 - t414;
t278 = t323 * t359 + t324 * t358;
t280 = pkin(3) * t326 + pkin(7) * t329 + t396;
t371 = t364 * t278 + t361 * t280 + t290 * t403 - t304 * t404;
t231 = -t307 * t402 + t395;
t370 = t286 * t320 - t351 * t317;
t344 = -pkin(4) * t364 + t352;
t293 = t317 * t336;
t285 = t364 * t290;
t283 = t339 * t337;
t276 = pkin(4) * t419 + t303;
t273 = t364 * t280;
t246 = pkin(4) * t374 + t277;
t245 = -pkin(8) * t419 + t408;
t241 = pkin(4) * t264 + t261;
t240 = pkin(4) * t336 - pkin(8) * t418 - t304 * t361 + t285;
t235 = -t402 * t419 + (t437 * t418 - t414) * t363 + t373 * t360;
t234 = -t437 * t282 + t339 * t329;
t230 = -pkin(8) * t374 + t371;
t229 = pkin(8) * t412 + pkin(4) * t326 - t278 * t361 + t273 + (-t297 + (pkin(8) * t337 - t290) * t361) * qJD(4);
t226 = -t239 * t360 + t432;
t1 = [0.2e1 * t390 * t439 + t438 * t440 + MDP(6) * t411 - MDP(7) * t413 + (-pkin(6) * t411 + t362 * t387) * MDP(9) + (pkin(6) * t413 + t365 * t387) * MDP(10) + (-t262 * t336 + t277 * t327 + t278 * t325 + t291 * t329 - t292 * t326 + t303 * t318 + t377) * MDP(11) + (t261 * t303 + t262 * t304 - t291 * t277 + t292 * t278 + (t343 + t380) * t396) * MDP(12) + (t263 * t418 + t307 * t373) * MDP(13) + (-(-t305 * t364 - t307 * t361) * t329 + (-t428 - t264 * t364 + (t305 * t361 - t307 * t364) * qJD(4)) * t337) * MDP(14) + (t263 * t336 + t307 * t326 + t309 * t337 + t320 * t373) * MDP(15) + (-t264 * t336 - t305 * t326 - t320 * t374 - t337 * t415) * MDP(16) + (t320 * t326 + t293) * MDP(17) + ((-t304 * t403 + t273) * t320 + t285 * t317 + (-t287 * t403 + t259) * t336 + t243 * t326 + t277 * t305 + t303 * t264 + t286 * t392 + ((-qJD(4) * t290 - t278) * t320 + (-qJD(4) * t268 - t262) * t336 - t286 * t329 + t377) * t361) * MDP(18) + (-t244 * t326 + t261 * t418 + t303 * t263 + t277 * t307 + t286 * t373 - t317 * t408 - t320 * t371 - t336 * t372) * MDP(19) + (-t231 * t283 - t234 * t376) * MDP(20) + (-t231 * t282 + t232 * t283 - t234 * t249 + t235 * t376) * MDP(21) + (t231 * t336 + t234 * t316 - t283 * t317 - t326 * t376) * MDP(22) + (-t232 * t336 - t235 * t316 - t249 * t326 - t282 * t317) * MDP(23) + (t316 * t326 + t293) * MDP(24) + ((t229 * t363 - t230 * t360) * t316 + (t240 * t363 - t245 * t360) * t317 + t386 * t336 + t226 * t326 + t246 * t249 + t276 * t232 + t241 * t282 + t247 * t235 + ((-t240 * t360 - t245 * t363) * t316 - t227 * t336) * qJD(5)) * MDP(25) + (-t227 * t326 + t276 * t231 + t247 * t234 + t236 * t336 - t241 * t283 - t246 * t376 + (-(-qJD(5) * t245 + t229) * t316 - t240 * t317 - t225 * t336) * t360 + (-(qJD(5) * t240 + t230) * t316 - t245 * t317 - t383 * t336) * t363) * MDP(26); -t410 * t439 + t367 * t438 + ((t291 - t295) * t325 + (-t317 * t358 - t318 * t359) * pkin(2)) * MDP(11) + (t291 * t294 - t292 * t295 + (-t261 * t359 + t262 * t358 - t343 * t399) * pkin(2)) * MDP(12) + (t307 * t384 + t428) * MDP(13) + ((t263 - t427) * t364 + (-t264 - t425) * t361) * MDP(14) + (t320 * t384 + t415 - t424) * MDP(15) + (t375 + t426) * MDP(16) + (-t261 * t364 + t352 * t264 - t294 * t305 + (-t351 * t403 - t272) * t320 + (t295 * t320 + t370) * t361) * MDP(18) + (t261 * t361 + t352 * t263 - t294 * t307 + (t351 * t404 + t409) * t320 + t370 * t364) * MDP(19) + (t231 * t340 - t376 * t407) * MDP(20) + (-t231 * t339 - t232 * t340 - t249 * t407 - t376 * t406) * MDP(21) + (t429 - t436) * MDP(22) + (t381 + t430) * MDP(23) + ((-t334 * t363 - t335 * t360) * t317 + t344 * t232 + t241 * t339 + (t360 * t378 - t363 * t379) * t316 + t382 * t249 - t406 * t247) * MDP(25) + (-(-t334 * t360 + t335 * t363) * t317 + t344 * t231 + t241 * t340 + (t360 * t379 + t363 * t378) * t316 - t382 * t376 + t407 * t247) * MDP(26) - ((-t292 + t294) * MDP(11) + t320 * MDP(17) + t243 * MDP(18) - t244 * MDP(19) + t316 * MDP(24) + t226 * MDP(25) - t227 * MDP(26)) * t327 + (MDP(9) * t362 * t367 + MDP(10) * t410) * pkin(1); (-t325 ^ 2 - t327 ^ 2) * MDP(11) + (t291 * t327 - t292 * t325 + t349) * MDP(12) + (t375 - t426) * MDP(18) + (-t320 ^ 2 * t364 - t415 - t424) * MDP(19) + (t381 - t430) * MDP(25) + (t429 + t436) * MDP(26); t307 * t305 * MDP(13) + (-t305 ^ 2 + t307 ^ 2) * MDP(14) + (t263 + t427) * MDP(15) + (-t264 + t425) * MDP(16) + t317 * MDP(17) + (t244 * t320 - t286 * t307 + t368) * MDP(18) + (t243 * t320 + t286 * t305 - t372) * MDP(19) + (t231 + t446) * MDP(22) + (-t232 - t445) * MDP(23) + (-(-t238 * t360 - t431) * t316 - t227 * qJD(5) + (-t249 * t307 - t316 * t402 + t363 * t317) * pkin(4) + t442) * MDP(25) + ((-t239 * t316 - t225) * t360 + (t238 * t316 - t383) * t363 + (t307 * t376 - t316 * t401 - t360 * t317) * pkin(4) + t443) * MDP(26) + t441; (t395 + t446) * MDP(22) + (-t385 - t445) * MDP(23) + (t227 * t316 + t442) * MDP(25) + (-t360 * t225 + t226 * t316 - t363 * t228 + t443) * MDP(26) + (-MDP(22) * t423 + MDP(23) * t376 - MDP(25) * t227 - MDP(26) * t432) * qJD(5) + t441;];
tauc = t1;
