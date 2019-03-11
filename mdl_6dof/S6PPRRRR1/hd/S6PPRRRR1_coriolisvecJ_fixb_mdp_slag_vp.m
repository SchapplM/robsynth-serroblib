% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:44
% EndTime: 2019-03-08 19:01:48
% DurationCPUTime: 3.20s
% Computational Cost: add. (3048->324), mult. (7944->482), div. (0->0), fcn. (6856->14), ass. (0->163)
t433 = pkin(9) + pkin(10);
t338 = cos(pkin(6));
t334 = sin(pkin(7));
t346 = cos(qJ(3));
t418 = t334 * t346;
t333 = sin(pkin(13));
t335 = sin(pkin(6));
t342 = sin(qJ(3));
t336 = cos(pkin(13));
t337 = cos(pkin(7));
t417 = t336 * t337;
t437 = (-t333 * t342 + t346 * t417) * t335;
t440 = t338 * t418 + t437;
t323 = qJD(1) * t338 + qJD(2);
t439 = qJD(1) * t437 + t323 * t418;
t341 = sin(qJ(4));
t345 = cos(qJ(4));
t438 = MDP(7) * (t341 ^ 2 - t345 ^ 2);
t408 = qJD(1) * t335;
t389 = t336 * t408;
t299 = t323 * t337 - t334 * t389;
t352 = (t333 * t346 + t342 * t417) * t335;
t419 = t334 * t342;
t394 = t323 * t419;
t280 = qJD(1) * t352 + t394;
t384 = qJD(3) * t433 + t280;
t257 = t299 * t345 - t341 * t384;
t271 = t439 * qJD(3);
t243 = qJD(4) * t257 + t345 * t271;
t256 = qJD(4) * pkin(4) + t257;
t344 = cos(qJ(5));
t436 = (qJD(5) * t256 + t243) * t344;
t376 = t337 * t389;
t390 = t333 * t408;
t403 = qJD(3) * t346;
t405 = qJD(3) * t342;
t272 = qJD(3) * t394 + t376 * t405 + t390 * t403;
t435 = qJD(3) * t280 - t272;
t402 = qJD(4) * t341;
t395 = pkin(4) * t402;
t375 = -t280 + t395;
t330 = qJD(4) + qJD(5);
t258 = t299 * t341 + t345 * t384;
t244 = -qJD(4) * t258 - t341 * t271;
t340 = sin(qJ(5));
t401 = qJD(5) * t340;
t382 = t244 * t340 - t258 * t401;
t223 = t382 + t436;
t427 = t258 * t344;
t235 = t256 * t340 + t427;
t383 = t243 * t340 - t244 * t344;
t224 = qJD(5) * t235 + t383;
t428 = t258 * t340;
t234 = t256 * t344 - t428;
t232 = -pkin(5) * t330 - t234;
t279 = -t342 * t390 + (t323 * t334 + t376) * t346;
t329 = -pkin(4) * t345 - pkin(3);
t270 = qJD(3) * t329 - t279;
t404 = qJD(3) * t345;
t322 = t344 * t404;
t406 = qJD(3) * t341;
t387 = t340 * t406;
t304 = -t322 + t387;
t306 = -t340 * t404 - t344 * t406;
t255 = pkin(5) * t304 + pkin(11) * t306 + t270;
t310 = t340 * t345 + t341 * t344;
t290 = t330 * t310;
t282 = t290 * qJD(3);
t309 = t340 * t341 - t344 * t345;
t286 = pkin(5) * t309 - pkin(11) * t310 + t329;
t289 = t330 * t309;
t315 = t433 * t341;
t316 = t433 * t345;
t298 = -t315 * t340 + t316 * t344;
t301 = qJD(6) + t304;
t391 = qJD(4) * t433;
t312 = t341 * t391;
t313 = t345 * t391;
t363 = -t315 * t344 - t316 * t340;
t414 = -qJD(5) * t363 - t279 * t309 + t312 * t344 + t313 * t340;
t434 = -(qJD(6) * t255 + t223) * t309 + t224 * t310 - t232 * t289 + (-qJD(6) * t286 + t414) * t301 - t298 * t282;
t432 = pkin(4) * t344;
t431 = qJD(3) * pkin(3);
t430 = t232 * t304;
t429 = t232 * t310;
t339 = sin(qJ(6));
t400 = qJD(6) * t339;
t397 = qJD(3) * qJD(4);
t386 = t345 * t397;
t281 = qJD(5) * t322 - t330 * t387 + t344 * t386;
t343 = cos(qJ(6));
t399 = qJD(6) * t343;
t411 = t281 * t343 + t330 * t399;
t260 = t306 * t400 + t411;
t426 = t260 * t339;
t425 = t281 * t339;
t424 = t286 * t282;
t420 = t306 * t339;
t293 = -t330 * t343 - t420;
t423 = t293 * t301;
t364 = t306 * t343 - t330 * t339;
t422 = t364 * t301;
t416 = t339 * t282;
t415 = t343 * t282;
t413 = pkin(5) * t290 + pkin(11) * t289 + t375;
t412 = qJD(5) * t298 - t279 * t310 - t312 * t340 + t313 * t344;
t409 = MDP(11) * t341;
t398 = qJD(6) * t346;
t396 = pkin(4) * t406;
t388 = t334 * t403;
t385 = -pkin(4) * t330 - t256;
t381 = t301 * t343;
t273 = -t279 - t431;
t380 = -qJD(3) * t273 - t271;
t285 = -pkin(5) * t306 + pkin(11) * t304;
t327 = pkin(4) * t340 + pkin(11);
t377 = qJD(6) * t327 + t285 + t396;
t269 = qJD(3) * t395 + t272;
t236 = t257 * t340 + t427;
t374 = pkin(4) * t401 - t236;
t237 = t257 * t344 - t428;
t373 = -qJD(5) * t432 + t237;
t233 = pkin(11) * t330 + t235;
t229 = t233 * t343 + t255 * t339;
t372 = t224 * t339 - t229 * t306 + t232 * t399;
t370 = -t282 * t327 + t430;
t369 = t233 * t339 - t255 * t343;
t288 = t338 * t419 + t352;
t300 = -t334 * t335 * t336 + t337 * t338;
t266 = -t288 * t341 + t300 * t345;
t267 = t288 * t345 + t300 * t341;
t247 = t266 * t340 + t267 * t344;
t368 = t247 * t343 - t339 * t440;
t367 = -t247 * t339 - t343 * t440;
t366 = t266 * t344 - t267 * t340;
t302 = t337 * t345 - t341 * t419;
t303 = t337 * t341 + t345 * t419;
t365 = t302 * t344 - t303 * t340;
t277 = t302 * t340 + t303 * t344;
t347 = qJD(4) ^ 2;
t362 = pkin(9) * t347 - t435;
t361 = qJD(4) * (t273 + t279 - t431);
t359 = -t224 * t343 + t232 * t400 - t306 * t369;
t358 = t270 * t306 - t383;
t357 = t270 * t304 - t382;
t356 = -t289 * t343 - t310 * t400;
t261 = -qJD(6) * t364 + t425;
t350 = ((t260 - t423) * t343 + (-t261 + t422) * t339) * MDP(21) + (-t364 * t381 + t426) * MDP(20) + (-t301 ^ 2 * t339 - t293 * t306 + t415) * MDP(23) + (t301 * t381 - t306 * t364 + t416) * MDP(22) + t281 * MDP(15) + (-t304 ^ 2 + t306 ^ 2) * MDP(14) + (-MDP(13) * t304 + MDP(24) * t301) * t306 + (t304 * MDP(15) + (-qJD(3) * t310 - t306) * MDP(16)) * t330;
t348 = qJD(3) ^ 2;
t328 = -pkin(5) - t432;
t292 = -qJD(4) * t303 - t341 * t388;
t291 = qJD(4) * t302 + t345 * t388;
t284 = t288 * qJD(3);
t283 = t440 * qJD(3);
t251 = qJD(4) * t266 + t283 * t345;
t250 = -qJD(4) * t267 - t283 * t341;
t249 = qJD(5) * t277 + t291 * t340 - t292 * t344;
t248 = qJD(5) * t365 + t291 * t344 + t292 * t340;
t245 = pkin(5) * t282 - pkin(11) * t281 + t269;
t242 = t343 * t245;
t226 = qJD(5) * t247 - t250 * t344 + t251 * t340;
t225 = qJD(5) * t366 + t250 * t340 + t251 * t344;
t1 = [(-t226 * t330 - t282 * t440 + t284 * t304) * MDP(18) + (-t225 * t330 - t281 * t440 - t284 * t306) * MDP(19) + ((-qJD(6) * t368 - t225 * t339 + t284 * t343) * t301 + t367 * t282 + t226 * t293 - t366 * t261) * MDP(25) + (-(qJD(6) * t367 + t225 * t343 + t284 * t339) * t301 - t368 * t282 - t226 * t364 - t366 * t260) * MDP(26) + (MDP(11) * t250 - MDP(12) * t251) * qJD(4) + (-t284 * MDP(4) - t283 * MDP(5) + (-t284 * t345 - t402 * t440) * MDP(11) + (-qJD(4) * t345 * t440 + t284 * t341) * MDP(12)) * qJD(3); ((-t248 * t339 - t277 * t399) * t301 - t277 * t416 + t249 * t293 - t365 * t261) * MDP(25) + (-(t248 * t343 - t277 * t400) * t301 - t277 * t415 - t249 * t364 - t365 * t260) * MDP(26) + (-MDP(18) * t249 - MDP(19) * t248) * t330 + (MDP(11) * t292 - MDP(12) * t291) * qJD(4) + ((-t282 * t346 + t304 * t405) * MDP(18) + (-t281 * t346 - t306 * t405) * MDP(19) + ((t339 * t398 + t343 * t405) * t301 - t346 * t415) * MDP(25) + (-(t339 * t405 - t343 * t398) * t301 + t346 * t416) * MDP(26) + (-MDP(12) * t345 - t409) * t346 * t397 + (-t346 * MDP(5) + (-MDP(11) * t345 + MDP(12) * t341 - MDP(4)) * t342) * t348) * t334; t435 * MDP(4) + (t279 - t439) * qJD(3) * MDP(5) - 0.2e1 * t397 * t438 + (t281 * t310 + t289 * t306) * MDP(13) + (-t281 * t309 - t282 * t310 + t289 * t304 + t290 * t306) * MDP(14) + (t269 * t309 + t270 * t290 + t282 * t329 + t304 * t375) * MDP(18) + (t269 * t310 - t270 * t289 + t281 * t329 - t306 * t375) * MDP(19) + (t260 * t310 * t343 - t356 * t364) * MDP(20) + (-(-t293 * t343 + t339 * t364) * t289 + (-t426 - t261 * t343 + (t293 * t339 + t343 * t364) * qJD(6)) * t310) * MDP(21) + (t260 * t309 - t290 * t364 + t301 * t356 + t310 * t415) * MDP(22) + (-t310 * t416 - t261 * t309 - t290 * t293 + (t289 * t339 - t310 * t399) * t301) * MDP(23) + (t282 * t309 + t290 * t301) * MDP(24) + (-t369 * t290 + t242 * t309 - t363 * t261 + t412 * t293 + (t424 + t413 * t301 + (-t233 * t309 - t298 * t301 + t429) * qJD(6)) * t343 + t434 * t339) * MDP(25) + (-t229 * t290 - t363 * t260 - t412 * t364 + (-t424 - (-qJD(6) * t233 + t245) * t309 - qJD(6) * t429 + (qJD(6) * t298 - t413) * t301) * t339 + t434 * t343) * MDP(26) + (-MDP(11) * t362 + MDP(12) * t361 + t347 * MDP(8)) * t345 + (MDP(11) * t361 + MDP(12) * t362 + 0.2e1 * MDP(6) * t386 - t347 * MDP(9)) * t341 + (-MDP(15) * t289 - MDP(16) * t290 - MDP(18) * t412 + MDP(19) * t414) * t330; (t306 * t396 + t237 * t330 + (qJD(5) * t385 - t243) * t344 + t357) * MDP(19) + (t328 * t260 + t370 * t343 - t374 * t364 + (t339 * t377 + t343 * t373) * t301 + t372) * MDP(26) + t348 * t438 + (-t304 * t396 + t236 * t330 + (t340 * t385 - t427) * qJD(5) + t358) * MDP(18) + t350 + (t328 * t261 + t370 * t339 + t374 * t293 + (t339 * t373 - t343 * t377) * t301 + t359) * MDP(25) + t380 * t409 + (-MDP(6) * t341 * t348 + MDP(12) * t380) * t345; (t358 + (-qJD(5) + t330) * t235) * MDP(18) + (t234 * t330 + t357 - t436) * MDP(19) + (-pkin(5) * t261 - (-t234 * t339 + t285 * t343) * t301 - t235 * t293 + t339 * t430 + (-t301 * t399 - t416) * pkin(11) + t359) * MDP(25) + (-pkin(5) * t260 + (t234 * t343 + t285 * t339) * t301 + t235 * t364 + t343 * t430 + (t301 * t400 - t415) * pkin(11) + t372) * MDP(26) + t350; -t364 * t293 * MDP(20) + (-t293 ^ 2 + t364 ^ 2) * MDP(21) + (t411 + t423) * MDP(22) + (-t422 - t425) * MDP(23) + t282 * MDP(24) + (-t223 * t339 + t229 * t301 + t232 * t364 + t242) * MDP(25) + (-t223 * t343 + t232 * t293 - t245 * t339 - t301 * t369) * MDP(26) + (MDP(22) * t420 + MDP(23) * t364 - MDP(25) * t229 + MDP(26) * t369) * qJD(6);];
tauc  = t1;
