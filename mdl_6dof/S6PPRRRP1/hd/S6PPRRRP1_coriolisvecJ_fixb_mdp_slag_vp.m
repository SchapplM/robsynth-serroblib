% Calculate Coriolis joint torque vector for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:23
% EndTime: 2021-01-16 00:50:32
% DurationCPUTime: 3.77s
% Computational Cost: add. (3211->367), mult. (8473->520), div. (0->0), fcn. (7103->12), ass. (0->165)
t339 = cos(pkin(6));
t326 = qJD(1) * t339 + qJD(2);
t334 = sin(pkin(12));
t336 = sin(pkin(6));
t342 = sin(qJ(3));
t345 = cos(qJ(3));
t337 = cos(pkin(12));
t338 = cos(pkin(7));
t413 = t337 * t338;
t349 = (t334 * t345 + t342 * t413) * t336;
t335 = sin(pkin(7));
t415 = t335 * t342;
t280 = qJD(1) * t349 + t326 * t415;
t341 = sin(qJ(4));
t344 = cos(qJ(4));
t363 = pkin(4) * t341 - pkin(10) * t344;
t318 = t363 * qJD(4);
t449 = t280 - t318;
t414 = t335 * t345;
t444 = (-t334 * t342 + t345 * t413) * t336;
t448 = t339 * t414 + t444;
t447 = t344 * MDP(12);
t446 = qJD(1) * t444 + t326 * t414;
t445 = MDP(6) * t341;
t332 = t341 ^ 2;
t443 = (-t344 ^ 2 + t332) * MDP(7);
t278 = qJD(3) * pkin(9) + t280;
t399 = qJD(1) * t336;
t379 = t337 * t399;
t296 = t326 * t338 - t335 * t379;
t262 = t344 * t278 + t341 * t296;
t259 = qJD(4) * pkin(10) + t262;
t364 = t338 * t379;
t380 = t334 * t399;
t279 = -t342 * t380 + t345 * (t326 * t335 + t364);
t321 = -pkin(4) * t344 - pkin(10) * t341 - pkin(3);
t270 = qJD(3) * t321 - t279;
t340 = sin(qJ(5));
t343 = cos(qJ(5));
t238 = -t259 * t340 + t343 * t270;
t393 = qJD(4) * t340;
t397 = qJD(3) * t341;
t314 = t343 * t397 + t393;
t235 = -qJ(6) * t314 + t238;
t395 = qJD(3) * t344;
t327 = -qJD(5) + t395;
t234 = -pkin(5) * t327 + t235;
t442 = t234 - t235;
t392 = qJD(4) * t341;
t412 = t340 * t344;
t434 = pkin(9) * t340;
t441 = -t279 * t412 + t449 * t343 - t392 * t434;
t389 = qJD(5) * t343;
t410 = t343 * t344;
t440 = -t279 * t410 + t321 * t389 - t449 * t340;
t439 = -t341 * t278 + t296 * t344;
t396 = qJD(3) * t342;
t378 = t335 * t396;
t394 = qJD(3) * t345;
t276 = t326 * t378 + t364 * t396 + t380 * t394;
t438 = qJD(3) * t280 - t276;
t390 = qJD(5) * t340;
t373 = t341 * t390;
t386 = qJD(3) * qJD(4);
t370 = t344 * t386;
t385 = qJD(4) * qJD(5);
t402 = (t370 + t385) * t343;
t293 = qJD(3) * t373 - t402;
t376 = t340 * t397;
t387 = t343 * qJD(4);
t312 = t376 - t387;
t437 = t312 * t327 - t293;
t372 = t341 * t389;
t391 = qJD(4) * t344;
t352 = t340 * t391 + t372;
t294 = qJD(3) * t352 + t340 * t385;
t436 = -t314 * t327 + t294;
t384 = MDP(18) + MDP(20);
t435 = pkin(5) * t312;
t433 = -qJ(6) - pkin(10);
t432 = qJD(3) * pkin(3);
t431 = qJ(6) * t341;
t424 = t270 * t340;
t239 = t259 * t343 + t424;
t236 = -qJ(6) * t312 + t239;
t430 = t236 * t327;
t275 = t446 * qJD(3);
t242 = t341 * t275 + t278 * t391 + t296 * t392;
t237 = pkin(5) * t294 + t242;
t429 = t237 * t340;
t428 = t237 * t343;
t427 = t242 * t340;
t426 = t242 * t343;
t258 = -qJD(4) * pkin(4) - t439;
t425 = t258 * t340;
t423 = t279 * t312;
t422 = t279 * t314;
t421 = t293 * t340;
t417 = t314 * t340;
t416 = t327 * t343;
t411 = t341 * t343;
t360 = pkin(5) * t341 - qJ(6) * t410;
t317 = t363 * qJD(3);
t366 = t343 * t317 - t340 * t439;
t369 = qJD(5) * t433;
t409 = qJD(3) * t360 + qJD(6) * t340 - t343 * t369 + t366;
t375 = t340 * t395;
t388 = qJD(6) * t343;
t405 = t340 * t317 + t343 * t439;
t408 = -qJ(6) * t375 - t340 * t369 - t388 + t405;
t328 = pkin(9) * t410;
t407 = t341 * t388 - t360 * qJD(4) - (-t328 + (-t321 + t431) * t340) * qJD(5) + t441;
t406 = -(-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t411 - (-qJD(6) * t341 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t344) * t340 - t440;
t401 = t340 * t321 + t328;
t383 = MDP(19) + MDP(21);
t377 = t335 * t394;
t374 = t327 * t390;
t371 = t341 * t386;
t368 = -qJD(6) - t435;
t241 = qJD(4) * t439 + t275 * t344;
t266 = qJD(3) * t318 + t276;
t367 = t241 * t340 - t343 * t266;
t365 = t343 * t241 - t259 * t390 + t340 * t266 + t270 * t389;
t285 = t339 * t415 + t349;
t299 = -t335 * t336 * t337 + t338 * t339;
t269 = t285 * t344 + t299 * t341;
t248 = t269 * t343 - t340 * t448;
t247 = -t269 * t340 - t343 * t448;
t268 = t285 * t341 - t299 * t344;
t361 = qJD(3) * t332 - t327 * t344;
t346 = qJD(4) ^ 2;
t359 = pkin(9) * t346 - t438;
t277 = -t279 - t432;
t358 = qJD(4) * (t277 + t279 - t432);
t301 = t338 * t341 + t344 * t415;
t288 = -t301 * t340 - t343 * t414;
t357 = -t301 * t343 + t340 * t414;
t300 = -t338 * t344 + t341 * t415;
t355 = qJ(6) * t293 - t367;
t354 = -qJ(6) * t294 + t365;
t347 = qJD(3) ^ 2;
t330 = -pkin(5) * t343 - pkin(4);
t323 = t433 * t343;
t322 = t433 * t340;
t319 = (pkin(5) * t340 + pkin(9)) * t341;
t311 = t343 * t321;
t309 = t312 ^ 2;
t295 = pkin(5) * t352 + pkin(9) * t391;
t290 = -t340 * t431 + t401;
t287 = qJD(4) * t301 + t341 * t377;
t286 = -qJD(4) * t300 + t344 * t377;
t283 = -qJ(6) * t411 + t311 + (-pkin(5) - t434) * t344;
t282 = t285 * qJD(3);
t281 = t448 * qJD(3);
t257 = qJD(5) * t357 - t286 * t340 + t343 * t378;
t256 = qJD(5) * t288 + t286 * t343 + t340 * t378;
t252 = pkin(5) * t375 + t262;
t249 = t258 - t368;
t246 = qJD(4) * t269 + t281 * t341;
t245 = -qJD(4) * t268 + t281 * t344;
t231 = -qJD(5) * t248 - t245 * t340 + t282 * t343;
t230 = qJD(5) * t247 + t245 * t343 + t282 * t340;
t229 = -qJD(6) * t312 + t354;
t228 = pkin(5) * t371 - qJD(5) * t239 - qJD(6) * t314 + t355;
t1 = [(-t230 * t312 - t231 * t314 + t247 * t293 - t248 * t294) * MDP(22) + (t228 * t247 + t229 * t248 + t230 * t236 + t231 * t234 + t237 * t268 + t246 * t249) * MDP(23) + t384 * (-t231 * t327 + t246 * t312 + t247 * t371 + t268 * t294) + t383 * (t230 * t327 + t246 * t314 - t248 * t371 - t268 * t293) + (-MDP(11) * t246 - MDP(12) * t245) * qJD(4) + (-t282 * MDP(4) - t281 * MDP(5) + (-t282 * t344 - t392 * t448) * MDP(11) + (t282 * t341 - t391 * t448) * MDP(12)) * qJD(3); (-t256 * t312 - t257 * t314 + t288 * t293 + t294 * t357) * MDP(22) + (t228 * t288 - t229 * t357 + t234 * t257 + t236 * t256 + t237 * t300 + t249 * t287) * MDP(23) + t384 * (-t257 * t327 + t287 * t312 + t288 * t371 + t294 * t300) + t383 * (t256 * t327 + t287 * t314 - t293 * t300 + t357 * t371) + (-MDP(11) * t287 - MDP(12) * t286) * qJD(4) + ((-t341 * MDP(11) - t447) * t345 * t386 + (-t345 * MDP(5) + (-t344 * MDP(11) + t341 * MDP(12) - MDP(4)) * t342) * t347) * t335; t438 * MDP(4) + (t279 - t446) * qJD(3) * MDP(5) + 0.2e1 * t370 * t445 - 0.2e1 * t386 * t443 + (t341 * t358 - t344 * t359) * MDP(11) + (t341 * t359 + t344 * t358) * MDP(12) + (-t293 * t411 + (t344 * t387 - t373) * t314) * MDP(13) + ((-t312 * t343 - t417) * t391 + (t421 - t294 * t343 + (t312 * t340 - t314 * t343) * qJD(5)) * t341) * MDP(14) + (t327 * t373 + t293 * t344 + (t314 * t341 + t343 * t361) * qJD(4)) * MDP(15) + (t327 * t372 + t294 * t344 + (-t312 * t341 - t340 * t361) * qJD(4)) * MDP(16) + (-t327 - t395) * MDP(17) * t392 + ((t321 * t390 + t441) * t327 + ((pkin(9) * t312 + t425) * qJD(4) + (t424 + (pkin(9) * t327 + t259) * t343) * qJD(5) + t367) * t344 + (t258 * t389 + pkin(9) * t294 + t427 - t423 + ((-pkin(9) * t412 + t311) * qJD(3) + t238) * qJD(4)) * t341) * MDP(18) + (t440 * t327 + (t258 * t387 + (qJD(4) * t314 - t374) * pkin(9) + t365) * t344 + (-t258 * t390 - pkin(9) * t293 + t426 - t422 + (-pkin(9) * t416 - qJD(3) * t401 - t239) * qJD(4)) * t341) * MDP(19) + (t294 * t319 + t295 * t312 + (t249 * t393 - t228) * t344 + t407 * t327 + (t249 * t389 + t429 - t423 + (qJD(3) * t283 + t234) * qJD(4)) * t341) * MDP(20) + (-t293 * t319 + t295 * t314 + (t249 * t387 + t229) * t344 - t406 * t327 + (-t249 * t390 + t428 - t422 + (-qJD(3) * t290 - t236) * qJD(4)) * t341) * MDP(21) + (t283 * t293 - t290 * t294 + t407 * t314 + t406 * t312 + (-t234 * t343 - t236 * t340) * t391 + (-t228 * t343 - t229 * t340 + (t234 * t340 - t236 * t343) * qJD(5)) * t341) * MDP(22) + (t228 * t283 + t229 * t290 + t237 * t319 + (-t279 * t341 + t295) * t249 - t406 * t236 - t407 * t234) * MDP(23) + (t344 * MDP(8) - t341 * MDP(9)) * t346; (qJD(4) * t262 - t277 * t397 - t242) * MDP(11) + (-qJD(3) * t277 - t275) * t447 + (-t314 * t416 - t421) * MDP(13) + (-t340 * t436 + t343 * t437) * MDP(14) + (-t327 * t389 + (t327 * t410 + (-t314 + t393) * t341) * qJD(3)) * MDP(15) + (t374 + (-t327 * t412 + (t312 + t387) * t341) * qJD(3)) * MDP(16) + t327 * MDP(17) * t397 + (-pkin(4) * t294 - t426 + t366 * t327 - t262 * t312 + (pkin(10) * t416 + t425) * qJD(5) + (-t238 * t341 + (-pkin(10) * t392 - t258 * t344) * t340) * qJD(3)) * MDP(18) + (pkin(4) * t293 + t427 - t405 * t327 - t262 * t314 + (-pkin(10) * t327 * t340 + t258 * t343) * qJD(5) + (-t258 * t410 + (-pkin(10) * t387 + t239) * t341) * qJD(3)) * MDP(19) + (-t428 - t252 * t312 + t294 * t330 + t409 * t327 + (t249 + t435) * t390 + (-t249 * t412 + (qJD(4) * t322 - t234) * t341) * qJD(3)) * MDP(20) + (t429 - t252 * t314 - t293 * t330 - t408 * t327 + (pkin(5) * t417 + t249 * t343) * qJD(5) + (-t249 * t410 + (qJD(4) * t323 + t236) * t341) * qJD(3)) * MDP(21) + (t293 * t322 + t294 * t323 + t409 * t314 + t408 * t312 + (t234 * t327 + t229) * t343 + (-t228 + t430) * t340) * MDP(22) + (t228 * t322 - t229 * t323 + t237 * t330 + (pkin(5) * t390 - t252) * t249 - t408 * t236 - t409 * t234) * MDP(23) + (-t344 * t445 + t443) * t347; -t309 * MDP(14) + t402 * MDP(15) + (-t239 * t327 - t367) * MDP(18) + (-t238 * t327 - t365) * MDP(19) + (t355 - t430) * MDP(20) + (-t235 * t327 - t354) * MDP(21) + pkin(5) * t293 * MDP(22) + (pkin(5) * t228 + t236 * t442) * MDP(23) + (-t327 * MDP(15) + t258 * MDP(19) + (qJD(6) + t249) * MDP(21) - t442 * MDP(22)) * t312 + (-MDP(16) * t412 + (0.2e1 * MDP(20) * pkin(5) + MDP(17)) * t341) * t386 + (-MDP(15) * t376 - MDP(16) * t314 - t239 * t384) * qJD(5) + (t312 * MDP(13) - t327 * MDP(16) - t258 * MDP(18) + (-t249 + t368) * MDP(20) - pkin(5) * t249 * MDP(23) + (-pkin(5) * MDP(21) + MDP(14)) * t314) * t314; t436 * MDP(20) + t437 * MDP(21) + (-t314 ^ 2 - t309) * MDP(22) + (t234 * t314 + t236 * t312 + t237) * MDP(23);];
tauc = t1;
