% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:32:23
% EndTime: 2019-03-08 21:32:30
% DurationCPUTime: 3.77s
% Computational Cost: add. (4193->394), mult. (10661->540), div. (0->0), fcn. (8044->10), ass. (0->164)
t357 = sin(qJ(3));
t360 = cos(qJ(3));
t441 = -qJ(4) - pkin(8);
t390 = qJD(3) * t441;
t322 = qJD(4) * t360 + t357 * t390;
t323 = -qJD(4) * t357 + t360 * t390;
t352 = sin(pkin(11));
t354 = cos(pkin(11));
t282 = t322 * t354 + t323 * t352;
t427 = t354 * t360;
t333 = t352 * t357 - t427;
t361 = cos(qJ(2));
t353 = sin(pkin(6));
t410 = qJD(1) * t353;
t397 = t361 * t410;
t306 = t333 * t397;
t417 = t282 + t306;
t407 = qJD(2) * t357;
t326 = qJD(2) * t427 - t352 * t407;
t321 = qJD(5) - t326;
t451 = MDP(5) * t357;
t450 = MDP(6) * (t357 ^ 2 - t360 ^ 2);
t334 = t352 * t360 + t354 * t357;
t418 = t322 * t352 - t354 * t323 - t334 * t397;
t327 = t334 * qJD(3);
t330 = t333 * qJD(3);
t406 = qJD(3) * t357;
t400 = pkin(3) * t406;
t284 = pkin(4) * t327 + pkin(9) * t330 + t400;
t399 = -pkin(3) * t360 - pkin(2);
t293 = pkin(4) * t333 - pkin(9) * t334 + t399;
t338 = t441 * t357;
t339 = t441 * t360;
t308 = t338 * t352 - t339 * t354;
t356 = sin(qJ(5));
t359 = cos(qJ(5));
t358 = sin(qJ(2));
t398 = t358 * t410;
t404 = qJD(5) * t359;
t405 = qJD(5) * t356;
t449 = -t293 * t404 + t308 * t405 - t417 * t359 + (-t284 + t398) * t356;
t415 = t356 * t293 + t359 * t308;
t336 = qJD(2) * pkin(8) + t398;
t381 = qJD(4) + t397;
t355 = cos(pkin(6));
t409 = qJD(1) * t355;
t396 = t357 * t409;
t448 = (-t336 * t360 - t396) * qJD(3) + (-qJD(3) * t360 * qJ(4) - t357 * t381) * qJD(2);
t447 = MDP(10) * t357 + MDP(11) * t360;
t446 = MDP(19) + MDP(21);
t445 = -MDP(20) + MDP(23);
t386 = qJ(4) * qJD(2) + t336;
t304 = t360 * t386 + t396;
t294 = t352 * t304;
t343 = t360 * t409;
t303 = -t357 * t386 + t343;
t298 = qJD(3) * pkin(3) + t303;
t257 = t298 * t354 - t294;
t254 = -qJD(3) * pkin(4) - t257;
t328 = t334 * qJD(2);
t403 = t359 * qJD(3);
t309 = t328 * t356 - t403;
t311 = qJD(3) * t356 + t328 * t359;
t242 = pkin(5) * t309 - qJ(6) * t311 + t254;
t319 = qJD(2) * t327;
t346 = pkin(3) * t352 + pkin(9);
t431 = t346 * t319;
t444 = t242 * t321 - t431;
t401 = qJD(2) * qJD(3);
t391 = t360 * t401;
t392 = t357 * t401;
t378 = -t352 * t392 + t354 * t391;
t277 = qJD(5) * t311 + t356 * t378;
t443 = t311 ^ 2;
t442 = pkin(5) * t319;
t440 = qJD(2) * pkin(2);
t439 = qJ(6) * t319;
t275 = (-t336 * t357 + t343) * qJD(3) + (-qJ(4) * t406 + t360 * t381) * qJD(2);
t245 = t275 * t352 - t354 * t448;
t276 = -qJD(5) * t403 + t328 * t405 - t359 * t378;
t230 = pkin(5) * t277 + qJ(6) * t276 - qJD(6) * t311 + t245;
t438 = t230 * t356;
t428 = t354 * t304;
t258 = t352 * t298 + t428;
t255 = qJD(3) * pkin(9) + t258;
t320 = qJD(2) * t399 + qJD(4) - t397;
t271 = -pkin(4) * t326 - pkin(9) * t328 + t320;
t240 = t255 * t359 + t271 * t356;
t437 = t240 * t321;
t436 = t276 * t356;
t435 = t309 * t311;
t434 = t309 * t326;
t387 = t311 * t321;
t433 = t321 * t356;
t432 = t334 * t359;
t430 = t353 * t358;
t429 = t353 * t361;
t314 = t356 * t319;
t362 = qJD(3) ^ 2;
t426 = t357 * t362;
t315 = t359 * t319;
t425 = t360 * t362;
t424 = qJ(6) * t327 + qJD(6) * t333 - t449;
t285 = -t306 * t356 - t359 * t398;
t423 = -pkin(5) * t327 + qJD(5) * t415 + t282 * t356 - t284 * t359 - t285;
t382 = pkin(5) * t356 - qJ(6) * t359;
t383 = pkin(5) * t359 + qJ(6) * t356;
t422 = -t382 * t330 + (qJD(5) * t383 - qJD(6) * t359) * t334 + t418;
t260 = t303 * t352 + t428;
t421 = qJD(6) * t356 - t321 * t382 + t260;
t262 = t303 * t354 - t294;
t283 = pkin(3) * t407 + pkin(4) * t328 - pkin(9) * t326;
t420 = t359 * t262 + t356 * t283;
t419 = -t356 * t277 - t309 * t404;
t416 = t326 * t433 + t315;
t414 = t321 * t404 + t314;
t408 = qJD(2) * t353;
t395 = t358 * t408;
t325 = pkin(3) * t392 + qJD(1) * t395;
t239 = -t255 * t356 + t271 * t359;
t402 = qJD(6) - t239;
t347 = -pkin(3) * t354 - pkin(4);
t394 = t361 * t408;
t393 = t346 * t405;
t246 = t354 * t275 + t352 * t448;
t268 = t319 * pkin(4) - pkin(9) * t378 + t325;
t372 = t359 * t246 - t255 * t405 + t356 * t268 + t271 * t404;
t228 = qJD(6) * t321 + t372 + t439;
t233 = -pkin(5) * t321 + t402;
t389 = -t233 * t326 + t228;
t385 = t356 * t246 + t255 * t404 - t359 * t268 + t271 * t405;
t229 = t385 - t442;
t234 = qJ(6) * t321 + t240;
t388 = t234 * t326 + t229;
t307 = -t354 * t338 - t339 * t352;
t380 = t233 * t359 - t234 * t356;
t379 = t245 * t334 - t308 * t319;
t331 = t355 * t357 + t360 * t430;
t377 = t355 * t360 - t357 * t430;
t289 = t354 * t331 + t352 * t377;
t273 = t289 * t356 + t359 * t429;
t274 = t289 * t359 - t356 * t429;
t376 = -t330 * t356 + t334 * t404;
t375 = t330 * t359 + t334 * t405;
t373 = t242 * t311 + t385;
t370 = t254 * t321 - t431;
t366 = -0.2e1 * qJD(3) * t440;
t363 = qJD(2) ^ 2;
t332 = t347 - t383;
t302 = -qJD(3) * t331 - t357 * t394;
t301 = qJD(3) * t377 + t360 * t394;
t288 = t331 * t352 - t354 * t377;
t270 = pkin(5) * t311 + qJ(6) * t309;
t263 = t334 * t382 + t307;
t261 = t301 * t354 + t302 * t352;
t259 = t301 * t352 - t354 * t302;
t252 = -pkin(5) * t333 - t293 * t359 + t308 * t356;
t251 = qJ(6) * t333 + t415;
t248 = t309 * t321 - t276;
t241 = -pkin(5) * t328 + t262 * t356 - t283 * t359;
t238 = qJ(6) * t328 + t420;
t237 = qJD(5) * t274 + t261 * t356 - t359 * t395;
t236 = -qJD(5) * t273 + t261 * t359 + t356 * t395;
t1 = [qJD(3) * t302 * MDP(10) - qJD(3) * t301 * MDP(11) + (t259 * t328 + t261 * t326 + t288 * t378 - t289 * t319) * MDP(12) + (t245 * t288 + t246 * t289 - t257 * t259 + t258 * t261) * MDP(13) + (-t236 * t309 + t237 * t311 - t273 * t276 - t274 * t277) * MDP(22) + (t228 * t274 + t229 * t273 + t230 * t288 + t233 * t237 + t234 * t236 + t242 * t259) * MDP(24) + (-t325 * t361 * MDP(13) + (t320 * t358 * MDP(13) - qJD(3) * t361 * t447) * qJD(2) + (-t361 * MDP(4) + (-MDP(10) * t360 + MDP(11) * t357 - MDP(3)) * t358) * t363) * t353 + t446 * (-t237 * t321 + t259 * t309 - t273 * t319 + t288 * t277) + t445 * (t236 * t321 - t259 * t311 + t274 * t319 + t276 * t288); 0.2e1 * t391 * t451 - 0.2e1 * t401 * t450 + MDP(7) * t425 - MDP(8) * t426 + (-pkin(8) * t425 + t357 * t366) * MDP(10) + (pkin(8) * t426 + t360 * t366) * MDP(11) + (-t246 * t333 + t257 * t330 - t258 * t327 + t307 * t378 + t326 * t417 + t328 * t418 + t379) * MDP(12) + (t246 * t308 + t245 * t307 + t325 * t399 + (-t398 + t400) * t320 + t417 * t258 - t418 * t257) * MDP(13) + (-t276 * t432 - t311 * t375) * MDP(14) + (-(-t309 * t359 - t311 * t356) * t330 + (t436 - t277 * t359 + (t309 * t356 - t311 * t359) * qJD(5)) * t334) * MDP(15) + (-t276 * t333 + t311 * t327 + t315 * t334 - t321 * t375) * MDP(16) + (-t277 * t333 - t309 * t327 - t314 * t334 - t321 * t376) * MDP(17) + (t319 * t333 + t321 * t327) * MDP(18) + (-t385 * t333 + t239 * t327 + t307 * t277 + t285 * t321 + t418 * t309 + ((-qJD(5) * t308 + t284) * t321 + t293 * t319 + t254 * qJD(5) * t334) * t359 + ((-qJD(5) * t293 - t282) * t321 - t254 * t330 + t379) * t356) * MDP(19) + (-t240 * t327 + t245 * t432 - t375 * t254 - t307 * t276 + t418 * t311 - t415 * t319 + t321 * t449 - t372 * t333) * MDP(20) + (-t229 * t333 - t233 * t327 + t242 * t376 - t252 * t319 + t263 * t277 + t309 * t422 - t321 * t423 + t334 * t438) * MDP(21) + (-t251 * t277 - t252 * t276 - t380 * t330 + t423 * t311 - t424 * t309 + (-t228 * t356 + t229 * t359 + (-t233 * t356 - t234 * t359) * qJD(5)) * t334) * MDP(22) + (t228 * t333 - t230 * t432 + t234 * t327 + t242 * t375 + t251 * t319 + t263 * t276 - t311 * t422 + t321 * t424) * MDP(23) + (t228 * t251 + t229 * t252 + t230 * t263 + t233 * t423 + t234 * t424 + t242 * t422) * MDP(24); ((t258 - t260) * t328 + (-t262 + t257) * t326 + (-t352 * t319 - t354 * t378) * pkin(3)) * MDP(12) + (t257 * t260 - t258 * t262 + (-t245 * t354 + t246 * t352 - t320 * t407) * pkin(3)) * MDP(13) + (t359 * t387 - t436) * MDP(14) + ((-t276 + t434) * t359 - t311 * t433 + t419) * MDP(15) + (-t321 * t326 * t359 - t311 * t328 + t414) * MDP(16) + (t309 * t328 - t321 * t405 + t416) * MDP(17) - t321 * t328 * MDP(18) + (-t239 * t328 - t260 * t309 + t347 * t277 + (-t245 + (-qJD(5) * t346 - t283) * t321) * t359 + (t262 * t321 + t370) * t356) * MDP(19) + (t240 * t328 + t245 * t356 - t260 * t311 - t347 * t276 + (t393 + t420) * t321 + t370 * t359) * MDP(20) + (-t230 * t359 + t233 * t328 + t277 * t332 + (-t346 * t404 + t241) * t321 - t421 * t309 + t444 * t356) * MDP(21) + (t238 * t309 - t241 * t311 + (-t277 * t346 + (t311 * t346 + t233) * qJD(5) + t389) * t359 + (-t276 * t346 + (t309 * t346 - t234) * qJD(5) + t388) * t356) * MDP(22) + (-t438 - t234 * t328 + t276 * t332 + (-t238 - t393) * t321 + t421 * t311 - t444 * t359) * MDP(23) + (t230 * t332 - t233 * t241 - t234 * t238 - t421 * t242 + (qJD(5) * t380 + t228 * t359 + t229 * t356) * t346) * MDP(24) + t447 * qJD(2) * t440 + (-t360 * t451 + t450) * t363; -t326 ^ 2 * MDP(12) + (-t258 * t326 + t325) * MDP(13) + t416 * MDP(19) + t419 * MDP(22) + t414 * MDP(23) + (-MDP(12) * t328 + t257 * MDP(13) - t242 * MDP(24) - t309 * t446 + t311 * t445) * t328 + (t319 * MDP(21) + (t276 + t434) * MDP(22) + (qJD(5) * t234 - t388) * MDP(24) + (-MDP(20) * t321 - t326 * MDP(23)) * t321) * t359 + (-t319 * MDP(20) + (qJD(5) * t233 + t389) * MDP(24) + MDP(22) * t387 + (-qJD(5) * MDP(19) - MDP(21) * t321) * t321) * t356; MDP(14) * t435 + (-t309 ^ 2 + t443) * MDP(15) + t248 * MDP(16) + (-t277 + t387) * MDP(17) + t319 * MDP(18) + (-t254 * t311 - t385 + t437) * MDP(19) + (t239 * t321 + t254 * t309 - t372) * MDP(20) + (-t270 * t309 - t373 + t437 + 0.2e1 * t442) * MDP(21) + (pkin(5) * t276 - qJ(6) * t277 + (t234 - t240) * t311 + (t233 - t402) * t309) * MDP(22) + (0.2e1 * t439 - t242 * t309 + t270 * t311 + (0.2e1 * qJD(6) - t239) * t321 + t372) * MDP(23) + (-pkin(5) * t229 + qJ(6) * t228 - t233 * t240 + t234 * t402 - t242 * t270) * MDP(24); (-qJD(3) * t328 + t435) * MDP(21) + t248 * MDP(22) + (-t321 ^ 2 - t443) * MDP(23) + (-t234 * t321 + t373 - t442) * MDP(24);];
tauc  = t1;
