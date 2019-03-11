% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:44:04
% EndTime: 2019-03-08 21:44:10
% DurationCPUTime: 2.84s
% Computational Cost: add. (1956->346), mult. (4627->482), div. (0->0), fcn. (2994->8), ass. (0->166)
t439 = pkin(4) + pkin(8);
t345 = -pkin(3) - pkin(9);
t340 = sin(qJ(3));
t406 = qJD(2) * t340;
t341 = sin(qJ(2));
t337 = sin(pkin(6));
t410 = qJD(1) * t337;
t383 = t341 * t410;
t308 = qJD(2) * pkin(8) + t383;
t343 = cos(qJ(3));
t338 = cos(pkin(6));
t409 = qJD(1) * t338;
t276 = t340 * t308 - t343 * t409;
t448 = qJD(4) + t276;
t442 = pkin(4) * t406 + t448;
t256 = t345 * qJD(3) + t442;
t371 = -qJ(4) * t340 - pkin(2);
t294 = t345 * t343 + t371;
t344 = cos(qJ(2));
t382 = t344 * t410;
t267 = t294 * qJD(2) - t382;
t339 = sin(qJ(5));
t342 = cos(qJ(5));
t240 = t256 * t339 + t267 * t342;
t390 = qJD(2) * qJD(3);
t374 = t340 * t390;
t319 = t339 * t374;
t403 = qJD(3) * t339;
t405 = qJD(2) * t343;
t297 = t342 * t405 + t403;
t398 = qJD(5) * t297;
t272 = -t319 + t398;
t377 = t339 * t405;
t401 = qJD(3) * t342;
t299 = -t377 + t401;
t400 = qJD(3) * t343;
t407 = qJD(2) * t337;
t376 = qJD(1) * t407;
t365 = t344 * t376;
t404 = qJD(3) * t338;
t449 = qJD(1) * t404 + t365;
t254 = t308 * t400 + t340 * t449;
t373 = t343 * t390;
t247 = pkin(4) * t373 + t254;
t363 = pkin(9) * t340 - qJ(4) * t343;
t399 = qJD(4) * t340;
t349 = t363 * qJD(3) - t399;
t413 = pkin(3) * t374 + t341 * t376;
t255 = t349 * qJD(2) + t413;
t367 = t342 * t247 - t255 * t339;
t229 = pkin(5) * t373 + qJ(6) * t272 - t240 * qJD(5) - qJD(6) * t299 + t367;
t235 = -qJ(6) * t297 + t240;
t327 = qJD(5) + t406;
t451 = t327 * t235 + t229;
t396 = qJD(5) * t342;
t412 = qJD(5) * t377 + t342 * t374;
t273 = qJD(3) * t396 - t412;
t385 = -t339 * t247 - t342 * t255 - t256 * t396;
t397 = qJD(5) * t339;
t230 = -qJ(6) * t273 - qJD(6) * t297 - t267 * t397 - t385;
t435 = t267 * t339;
t239 = t342 * t256 - t435;
t234 = -qJ(6) * t299 + t239;
t233 = pkin(5) * t327 + t234;
t450 = -t327 * t233 + t230;
t335 = t340 ^ 2;
t336 = t343 ^ 2;
t447 = MDP(6) * (t335 - t336);
t432 = t299 * t327;
t305 = t439 * t400;
t423 = t342 * t344;
t444 = -(-t339 * t341 + t340 * t423) * t410 + t342 * t305;
t402 = qJD(3) * t340;
t329 = pkin(3) * t402;
t279 = t329 + t349;
t314 = t439 * t340;
t426 = t339 * t344;
t443 = (t340 * t426 + t341 * t342) * t410 - t342 * t279 + t294 * t397 - t339 * t305 - t314 * t396;
t414 = t342 * t294 + t339 * t314;
t277 = t343 * t308 + t340 * t409;
t334 = qJD(3) * qJ(4);
t270 = -t334 - t277;
t269 = -qJD(3) * pkin(3) + t448;
t354 = -qJ(4) * t400 - t399;
t260 = t354 * qJD(2) + t413;
t285 = t329 + t354;
t346 = qJD(3) ^ 2;
t441 = (-t285 + t383) * qJD(2) - pkin(8) * t346 - t260;
t440 = t299 ^ 2;
t438 = qJD(2) * pkin(2);
t333 = qJD(3) * qJD(4);
t384 = -t308 * t402 + t343 * t449;
t249 = -t333 - t384;
t243 = -pkin(4) * t374 - t249;
t437 = t243 * t339;
t436 = t243 * t342;
t434 = t272 * t342;
t312 = -pkin(3) * t343 + t371;
t408 = qJD(2) * t312;
t278 = -t382 + t408;
t433 = t278 * t341;
t431 = t299 * t343;
t430 = t327 * t342;
t429 = t327 * t345;
t428 = t337 * t341;
t427 = t339 * t340;
t425 = t340 * t346;
t424 = t342 * t343;
t422 = t343 * t346;
t421 = qJ(6) - t345;
t368 = qJ(6) * t343 - t294;
t394 = qJD(6) * t343;
t420 = pkin(5) * t400 + t368 * t396 + (-qJ(6) * t402 - qJD(5) * t314 - t279 + t394) * t339 + t444;
t395 = qJD(5) * t343;
t379 = t339 * t395;
t419 = -t342 * t394 + (t340 * t401 + t379) * qJ(6) - t443;
t418 = t233 - t234;
t266 = pkin(4) * t405 + t277;
t330 = pkin(3) * t406;
t283 = t363 * qJD(2) + t330;
t417 = t339 * t266 + t342 * t283;
t366 = t342 * t266 - t283 * t339;
t416 = -qJD(6) * t342 + t421 * t397 - (pkin(5) * t343 - qJ(6) * t427) * qJD(2) - t366;
t311 = t421 * t342;
t415 = -qJ(6) * t342 * t406 - qJD(5) * t311 - qJD(6) * t339 - t417;
t315 = t439 * t343;
t257 = t334 + t266;
t393 = t257 * qJD(5);
t389 = -MDP(10) + MDP(13);
t388 = MDP(11) - MDP(14);
t387 = t340 * t428;
t347 = qJD(2) ^ 2;
t386 = t340 * t343 * t347;
t381 = t341 * t407;
t380 = t344 * t407;
t378 = t342 * t395;
t372 = MDP(20) * t400;
t362 = -qJD(2) * t336 + t327 * t340;
t361 = t327 ^ 2;
t358 = qJD(3) * t276 + t384;
t357 = qJD(3) * t277 - t254;
t288 = -t338 * t343 + t387;
t356 = -t288 * t339 + t337 * t423;
t263 = t288 * t342 + t337 * t426;
t289 = t338 * t340 + t343 * t428;
t355 = t257 * t340 + t345 * t400;
t309 = -t382 - t438;
t351 = qJD(3) * (t309 + t382 - t438);
t350 = qJD(3) * (-t278 - t382 - t408);
t238 = pkin(5) * t273 + t243;
t348 = -t249 * t343 + t254 * t340 + (t269 * t343 + t270 * t340) * qJD(3);
t320 = t342 * t373;
t310 = t421 * t339;
t304 = t439 * t402;
t303 = -qJ(4) * t405 + t330;
t302 = t342 * t314;
t296 = t297 ^ 2;
t268 = t278 * t406;
t262 = qJD(3) * t289 + t340 * t380;
t261 = -qJD(3) * t387 + (t380 + t404) * t343;
t253 = -qJ(6) * t424 + t414;
t248 = pkin(5) * t340 + t368 * t339 + t302;
t246 = pkin(5) * t297 + qJD(6) + t257;
t237 = t263 * qJD(5) + t262 * t339 + t342 * t381;
t236 = t356 * qJD(5) + t262 * t342 - t339 * t381;
t1 = [(-t249 * t289 + t254 * t288 - t261 * t270 + t262 * t269) * MDP(15) + (t236 * t327 + t261 * t297 + t273 * t289) * MDP(21) + (-t237 * t327 + t261 * t299 - t272 * t289) * MDP(22) + (-t236 * t299 - t237 * t297 + t263 * t272 + t273 * t356) * MDP(23) + (t229 * t263 - t230 * t356 + t233 * t236 + t235 * t237 + t238 * t289 + t246 * t261) * MDP(24) + (t261 * t343 + t262 * t340) * MDP(12) * qJD(2) + (t389 * t262 - t388 * t261 + (-t289 * t340 * MDP(12) + (MDP(12) * t288 + MDP(21) * t263 + MDP(22) * t356) * t343) * qJD(2)) * qJD(3) + (-t260 * t344 * MDP(15) + (MDP(15) * t433 + (t389 * t340 - t388 * t343) * t344 * qJD(3)) * qJD(2) + (-t344 * MDP(4) + (t388 * t340 + t389 * t343 - MDP(3)) * t341) * t347) * t337; 0.2e1 * t340 * MDP(5) * t373 - 0.2e1 * t390 * t447 + MDP(7) * t422 - MDP(8) * t425 + (-pkin(8) * t422 + t340 * t351) * MDP(10) + (pkin(8) * t425 + t343 * t351) * MDP(11) + ((-t335 - t336) * t365 + t348) * MDP(12) + (t340 * t350 - t441 * t343) * MDP(13) + (t441 * t340 + t343 * t350) * MDP(14) + (t260 * t312 + t278 * t285 + (-t433 + (-t269 * t340 + t270 * t343) * t344) * t410 + t348 * pkin(8)) * MDP(15) + (t272 * t339 * t343 + (t339 * t402 - t378) * t299) * MDP(16) + ((-t297 * t339 + t299 * t342) * t402 + (t434 + t273 * t339 + (t297 * t342 + t299 * t339) * qJD(5)) * t343) * MDP(17) + (-t327 * t378 - t272 * t340 + (t362 * t339 + t431) * qJD(3)) * MDP(18) + (t327 * t379 - t273 * t340 + (-t297 * t343 + t362 * t342) * qJD(3)) * MDP(19) + (t327 + t406) * t372 + (t315 * t273 - t304 * t297 + (-t257 * t401 + t367) * t340 + (-t279 * t339 + t444) * t327 + (-t240 * t340 - t414 * t327) * qJD(5) + (-t297 * t382 - t339 * t393 + t436 + ((-t294 * t339 + t302) * qJD(2) + t239) * qJD(3)) * t343) * MDP(21) + (-t315 * t272 - t304 * t299 + ((qJD(3) * t257 + qJD(5) * t267) * t339 + t385) * t340 + t443 * t327 + (-t299 * t382 - t342 * t393 - t437 + (-t414 * qJD(2) - t240) * qJD(3)) * t343) * MDP(22) + (t248 * t272 - t253 * t273 - t420 * t299 - t419 * t297 + (-t233 * t339 + t235 * t342) * t402 + (t229 * t339 - t230 * t342 + (t233 * t342 + t235 * t339) * qJD(5)) * t343) * MDP(23) + (t230 * t253 + t229 * t248 + t238 * (pkin(5) * t424 + t315) + t419 * t235 + t420 * t233 + ((-pkin(5) * t397 - t382) * t343 + (-pkin(5) * t342 - t439) * t402) * t246) * MDP(24); -MDP(5) * t386 + t347 * t447 + (-t309 * t406 + t357) * MDP(10) + (-t309 * t405 - t358) * MDP(11) + (-t303 * t405 + t268 - t357) * MDP(13) + (0.2e1 * t333 + (t278 * t343 + t303 * t340) * qJD(2) + t358) * MDP(14) + (-pkin(3) * t254 - qJ(4) * t249 - t269 * t277 - t270 * t448 - t278 * t303) * MDP(15) + (-t339 * t432 - t434) * MDP(16) + ((-t273 - t432) * t342 + (t297 * t327 + t272) * t339) * MDP(17) + (-t327 * t397 + t320 + (-t327 * t427 - t431) * qJD(2)) * MDP(18) + (-t327 * t396 + (-t340 * t430 + (t297 - t403) * t343) * qJD(2)) * MDP(19) - t327 * MDP(20) * t405 + (qJ(4) * t273 + t437 - t366 * t327 + t442 * t297 + (t257 * t342 - t339 * t429) * qJD(5) + (-t239 * t343 + t355 * t342) * qJD(2)) * MDP(21) + (-qJ(4) * t272 + t436 + t417 * t327 + t442 * t299 + (-t257 * t339 - t342 * t429) * qJD(5) + (t240 * t343 - t355 * t339) * qJD(2)) * MDP(22) + (-t272 * t311 + t273 * t310 - t415 * t297 - t416 * t299 - t339 * t450 - t342 * t451) * MDP(23) + (-t230 * t310 - t229 * t311 + t238 * (pkin(5) * t339 + qJ(4)) + (pkin(5) * t430 + t442) * t246 + t415 * t235 + t416 * t233) * MDP(24); MDP(13) * t386 + (-t335 * t347 - t346) * MDP(14) + (t268 + t254) * MDP(15) + t320 * MDP(21) + (t270 * MDP(15) - t297 * MDP(21) - t299 * MDP(22) - t246 * MDP(24)) * qJD(3) + ((-t297 * t406 + t272 - t398) * MDP(23) + t451 * MDP(24) - MDP(22) * t361) * t342 + (-MDP(22) * t373 + (-t273 + t432) * MDP(23) + t450 * MDP(24) - MDP(21) * t361) * t339; (-t296 + t440) * MDP(17) + t319 * MDP(18) + (t412 + t432) * MDP(19) + qJD(2) * t372 + (t240 * t327 - t257 * t299 + t367) * MDP(21) + (t239 * t327 + t385) * MDP(22) + t418 * MDP(24) * t235 + (t272 * MDP(23) + (-t246 * t299 + t229) * MDP(24)) * pkin(5) + (t299 * MDP(16) + t327 * MDP(18) + t257 * MDP(22) - t418 * MDP(23)) * t297 + (-t297 * MDP(18) - MDP(19) * t401 - t240 * MDP(21) + MDP(22) * t435) * qJD(5); (-t296 - t440) * MDP(23) + (t233 * t299 + t235 * t297 + t238) * MDP(24);];
tauc  = t1;
