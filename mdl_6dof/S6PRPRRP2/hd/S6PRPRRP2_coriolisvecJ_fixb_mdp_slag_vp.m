% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:16
% EndTime: 2019-03-08 20:03:21
% DurationCPUTime: 3.06s
% Computational Cost: add. (2612->357), mult. (6424->492), div. (0->0), fcn. (4737->10), ass. (0->153)
t316 = cos(pkin(11));
t320 = sin(qJ(2));
t315 = sin(pkin(6));
t378 = qJD(1) * t315;
t356 = t320 * t378;
t300 = t316 * t356;
t314 = sin(pkin(11));
t323 = cos(qJ(2));
t355 = t323 * t378;
t273 = t314 * t355 + t300;
t319 = sin(qJ(4));
t322 = cos(qJ(4));
t344 = pkin(4) * t319 - pkin(9) * t322;
t296 = t344 * qJD(4);
t416 = t273 - t296;
t415 = MDP(12) * t322;
t376 = qJD(2) * t322;
t306 = -qJD(5) + t376;
t414 = MDP(6) * t319;
t312 = t319 ^ 2;
t413 = MDP(7) * (-t322 ^ 2 + t312);
t278 = (t314 * t320 - t316 * t323) * t315;
t275 = qJD(2) * t278;
t269 = qJD(1) * t275;
t298 = qJD(2) * pkin(2) + t355;
t267 = t314 * t298 + t300;
t264 = qJD(2) * pkin(8) + t267;
t317 = cos(pkin(6));
t305 = qJD(1) * t317 + qJD(3);
t409 = -t319 * t264 + t305 * t322;
t225 = qJD(4) * t409 - t269 * t322;
t249 = t322 * t264 + t319 * t305;
t244 = qJD(4) * pkin(9) + t249;
t299 = t314 * t356;
t266 = t298 * t316 - t299;
t334 = -pkin(4) * t322 - pkin(9) * t319 - pkin(3);
t252 = qJD(2) * t334 - t266;
t279 = (t314 * t323 + t316 * t320) * t315;
t253 = (qJD(1) * t279 + t296) * qJD(2);
t318 = sin(qJ(5));
t321 = cos(qJ(5));
t370 = qJD(5) * t321;
t371 = qJD(5) * t318;
t345 = -t318 * t225 - t244 * t370 - t252 * t371 + t321 * t253;
t365 = qJD(2) * qJD(4);
t349 = t319 * t365;
t214 = -pkin(5) * t349 - t345;
t221 = t244 * t321 + t252 * t318;
t219 = -qJ(6) * t306 + t221;
t412 = t219 * t306 + t214;
t276 = t316 * t355 - t299;
t391 = t318 * t322;
t411 = -t276 * t391 + t416 * t321;
t408 = pkin(2) * t316;
t288 = t334 - t408;
t389 = t321 * t322;
t410 = -t276 * t389 + t288 * t370 - t416 * t318;
t308 = pkin(2) * t314 + pkin(8);
t381 = t318 * t288 + t308 * t389;
t363 = MDP(18) + MDP(20);
t362 = MDP(19) - MDP(22);
t373 = qJD(4) * t322;
t374 = qJD(4) * t319;
t226 = t264 * t373 - t319 * t269 + t305 * t374;
t348 = t322 * t365;
t352 = t319 * t371;
t364 = qJD(4) * qJD(5);
t270 = qJD(2) * t352 + (-t348 - t364) * t321;
t347 = t318 * t364;
t351 = t319 * t370;
t271 = t347 + (t318 * t373 + t351) * qJD(2);
t375 = qJD(4) * t318;
t377 = qJD(2) * t319;
t292 = t321 * t377 + t375;
t217 = pkin(5) * t271 + qJ(6) * t270 - qJD(6) * t292 + t226;
t407 = t217 * t318;
t406 = t217 * t321;
t404 = t226 * t318;
t403 = t226 * t321;
t243 = -qJD(4) * pkin(4) - t409;
t402 = t243 * t318;
t401 = t243 * t321;
t367 = t321 * qJD(4);
t290 = t318 * t377 - t367;
t400 = t276 * t290;
t399 = t276 * t292;
t398 = t288 * t321;
t397 = t290 * t306;
t394 = t306 * t321;
t393 = t308 * t318;
t392 = t308 * t321;
t390 = t319 * t321;
t388 = (-t308 * t371 - qJD(6)) * t322 + (qJ(6) - t392) * t374 + t410;
t346 = pkin(5) + t393;
t387 = qJD(5) * t381 - t346 * t374 + t411;
t341 = pkin(5) * t318 - qJ(6) * t321;
t386 = qJD(6) * t318 + t306 * t341 + t249;
t295 = t344 * qJD(2);
t385 = t318 * t295 + t321 * t409;
t354 = t322 * t367;
t384 = -t271 * t390 - t290 * t354;
t350 = t312 * t365;
t303 = t321 * t350;
t382 = -t306 * t354 + t303;
t379 = MDP(21) * t318;
t372 = qJD(5) * t290;
t369 = qJD(6) * t306;
t368 = t243 * qJD(5);
t220 = -t244 * t318 + t252 * t321;
t366 = qJD(6) - t220;
t361 = pkin(9) * t306 * t318;
t360 = pkin(9) * t394;
t359 = pkin(9) * t374;
t358 = pkin(9) * t367;
t357 = t321 * t225 + t252 * t370 + t318 * t253;
t353 = t306 * t370;
t342 = pkin(5) * t321 + qJ(6) * t318;
t328 = t244 * t371 - t357;
t213 = qJ(6) * t349 - t328 - t369;
t340 = t213 * t321 + t214 * t318;
t218 = pkin(5) * t306 + t366;
t339 = t218 * t321 - t219 * t318;
t338 = t218 * t318 + t219 * t321;
t337 = t295 * t321 - t318 * t409;
t262 = t279 * t322 + t317 * t319;
t236 = t262 * t321 + t278 * t318;
t235 = t262 * t318 - t278 * t321;
t261 = t279 * t319 - t317 * t322;
t332 = t308 + t341;
t274 = qJD(2) * t279;
t268 = qJD(1) * t274;
t324 = qJD(4) ^ 2;
t331 = qJD(2) * t273 - t308 * t324 - t268;
t263 = -qJD(2) * pkin(3) - t266;
t330 = qJD(4) * (qJD(2) * (-pkin(3) - t408) + t263 + t276);
t326 = -t220 * t306 + t328;
t325 = qJD(2) ^ 2;
t302 = -pkin(4) - t342;
t282 = t292 * t374;
t272 = t332 * t319;
t260 = pkin(5) * t292 + qJ(6) * t290;
t255 = t322 * t346 - t398;
t254 = -qJ(6) * t322 + t381;
t247 = -t270 - t397;
t240 = (qJD(5) * t342 - qJD(6) * t321) * t319 + t332 * t373;
t234 = -qJD(4) * t261 - t275 * t322;
t233 = qJD(4) * t262 - t275 * t319;
t230 = -pkin(5) * t377 - t337;
t229 = qJ(6) * t377 + t385;
t224 = pkin(5) * t290 - qJ(6) * t292 + t243;
t216 = -qJD(5) * t235 + t234 * t321 + t274 * t318;
t215 = qJD(5) * t236 + t234 * t318 - t274 * t321;
t1 = [(-t266 * t274 - t267 * t275 + t268 * t278 - t269 * t279) * MDP(5) + (t215 * t292 - t216 * t290 - t235 * t270 - t236 * t271) * MDP(21) + (t213 * t236 + t214 * t235 + t215 * t218 + t216 * t219 + t217 * t261 + t224 * t233) * MDP(23) + (-MDP(3) * t320 - MDP(4) * t323) * t325 * t315 + (-MDP(11) * t233 - MDP(12) * t234) * qJD(4) + ((-t322 * MDP(11) + t319 * MDP(12)) * t274 + (t278 * t415 + (MDP(11) * t278 - t235 * t363 - t236 * t362) * t319) * qJD(4)) * qJD(2) + t363 * (t215 * t306 + t233 * t290 + t261 * t271) + t362 * (t216 * t306 + t233 * t292 - t261 * t270); (t266 * t273 - t267 * t276 + (-t268 * t316 - t269 * t314) * pkin(2)) * MDP(5) + 0.2e1 * t348 * t414 - 0.2e1 * t365 * t413 + (t319 * t330 + t322 * t331) * MDP(11) + (-t319 * t331 + t322 * t330) * MDP(12) + (-t270 * t390 + (-t352 + t354) * t292) * MDP(13) + (-t292 * t351 + (-t292 * t373 + (t270 + t372) * t319) * t318 + t384) * MDP(14) + (t270 * t322 + t306 * t352 + t282 + t382) * MDP(15) + (t306 * t351 + t271 * t322 + (-t319 * t290 + (-qJD(2) * t312 + t306 * t322) * t318) * qJD(4)) * MDP(16) + (-t306 - t376) * MDP(17) * t374 + ((t288 * t371 + t411) * t306 + (t308 * t353 + (t290 * t308 + t402) * qJD(4) - t345) * t322 + (t321 * t368 + t404 + t308 * t271 - t400 + (-t306 * t393 + (-t308 * t391 + t398) * qJD(2) + t220) * qJD(4)) * t319) * MDP(18) + (t410 * t306 + ((-t306 * t308 - t244) * t371 + (t292 * t308 + t401) * qJD(4) + t357) * t322 + (-t318 * t368 + t403 - t308 * t270 - t399 + (-qJD(2) * t381 - t306 * t392 - t221) * qJD(4)) * t319) * MDP(19) + (t240 * t290 + t271 * t272 + (t224 * t375 + t214) * t322 + t387 * t306 + (t224 * t370 + t407 - t400 + (-qJD(2) * t255 - t218) * qJD(4)) * t319) * MDP(20) + (-t254 * t271 - t255 * t270 + t387 * t292 - t388 * t290 + t339 * t373 + (-qJD(5) * t338 - t213 * t318 + t214 * t321) * t319) * MDP(21) + (-t240 * t292 + t272 * t270 + (-t224 * t367 - t213) * t322 - t388 * t306 + (t224 * t371 - t406 + t399 + (qJD(2) * t254 + t219) * qJD(4)) * t319) * MDP(22) + (t213 * t254 + t214 * t255 + t217 * t272 + (-t319 * t276 + t240) * t224 + t388 * t219 + t387 * t218) * MDP(23) + (MDP(8) * t322 - MDP(9) * t319) * t324; (t282 - t303) * MDP(19) + t384 * MDP(21) + t382 * MDP(22) + (-t324 * MDP(12) - t217 * MDP(23) - t363 * t271 + t362 * t270 + (t292 * t379 + t338 * MDP(23) + (t321 * MDP(19) + t318 * t363) * t306) * qJD(4)) * t322 + (-t324 * MDP(11) - t270 * t379 - qJD(4) * t292 * MDP(22) + (qJD(4) * t224 + t340) * MDP(23) + ((t290 * t318 + t292 * t321) * MDP(21) + t339 * MDP(23) + (-t318 * t362 + t321 * t363) * t306) * qJD(5)) * t319 + t363 * (t290 * t374 - t318 * t350); (qJD(4) * t249 - t263 * t377 - t226) * MDP(11) + (-qJD(2) * t263 + t269) * t415 + (-t270 * t318 - t292 * t394) * MDP(13) + ((-t270 + t397) * t321 + (t292 * t306 - t271) * t318) * MDP(14) + (-t353 + (t306 * t389 + (-t292 + t375) * t319) * qJD(2)) * MDP(15) + (t306 * t371 + (-t306 * t391 + (t290 + t367) * t319) * qJD(2)) * MDP(16) + t306 * MDP(17) * t377 + (-pkin(4) * t271 - t403 + t337 * t306 - t249 * t290 + (t360 + t402) * qJD(5) + (-t220 * t319 + (-t243 * t322 - t359) * t318) * qJD(2)) * MDP(18) + (pkin(4) * t270 + t404 - t385 * t306 - t249 * t292 + (-t361 + t401) * qJD(5) + (-t243 * t389 + (t221 - t358) * t319) * qJD(2)) * MDP(19) + (-t406 - t230 * t306 + t271 * t302 - t386 * t290 + (t224 * t318 + t360) * qJD(5) + (t218 * t319 + (-t224 * t322 - t359) * t318) * qJD(2)) * MDP(20) + (t229 * t290 - t230 * t292 + (t213 - t306 * t218 + (qJD(5) * t292 - t271) * pkin(9)) * t321 + ((-t270 + t372) * pkin(9) + t412) * t318) * MDP(21) + (-t407 + t229 * t306 + t270 * t302 + t386 * t292 + (-t224 * t321 + t361) * qJD(5) + (t224 * t389 + (-t219 + t358) * t319) * qJD(2)) * MDP(22) + (t217 * t302 - t218 * t230 - t219 * t229 - t386 * t224 + (qJD(5) * t339 + t340) * pkin(9)) * MDP(23) + (-t322 * t414 + t413) * t325; t247 * MDP(15) - MDP(16) * t347 + t326 * MDP(19) + (pkin(5) * t270 - qJ(6) * t271) * MDP(21) + (-t326 - 0.2e1 * t369) * MDP(22) + (-pkin(5) * t214 + qJ(6) * t213 - t218 * t221 + t219 * t366 - t224 * t260) * MDP(23) + (-t306 * MDP(16) - t243 * MDP(18) - t224 * MDP(20) + (t219 - t221) * MDP(21) + t260 * MDP(22) + MDP(14) * t292) * t292 + (-MDP(16) * t351 + (-MDP(16) * t391 + (0.2e1 * pkin(5) * MDP(20) + 0.2e1 * qJ(6) * MDP(22) + MDP(17)) * t319) * qJD(4)) * qJD(2) + (t292 * MDP(13) + t243 * MDP(19) - t260 * MDP(20) + (t218 - t366) * MDP(21) - t224 * MDP(22) - MDP(14) * t290) * t290 + t363 * (-t221 * t306 + t345); (t290 * t292 - t349) * MDP(20) + t247 * MDP(21) + (-t292 ^ 2 - t306 ^ 2) * MDP(22) + (t224 * t292 + t412) * MDP(23);];
tauc  = t1;
