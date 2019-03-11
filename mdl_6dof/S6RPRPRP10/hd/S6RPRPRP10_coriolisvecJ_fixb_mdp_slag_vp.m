% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:38
% EndTime: 2019-03-09 03:32:44
% DurationCPUTime: 3.03s
% Computational Cost: add. (2271->376), mult. (4339->481), div. (0->0), fcn. (2205->4), ass. (0->160)
t364 = MDP(23) + MDP(25);
t317 = -pkin(1) - pkin(7);
t292 = qJD(1) * t317 + qJD(2);
t306 = qJD(3) * qJD(4);
t315 = cos(qJ(3));
t376 = qJD(3) * t315;
t262 = -t292 * t376 - t306;
t281 = t315 * t292;
t409 = qJD(3) * pkin(3);
t345 = -qJD(4) + t409;
t264 = -t281 - t345;
t313 = sin(qJ(3));
t280 = t313 * t292;
t308 = qJD(3) * qJ(4);
t266 = -t280 - t308;
t399 = t266 * t315;
t417 = qJD(3) * (t313 * (-t264 + t281) + t399) + t262 * t313;
t380 = qJD(1) * t315;
t298 = qJD(5) + t380;
t416 = t298 ^ 2;
t310 = t313 ^ 2;
t311 = t315 ^ 2;
t415 = MDP(8) * (t310 - t311);
t414 = MDP(6) * qJ(2) + MDP(5);
t413 = t313 * pkin(3) + qJ(2);
t407 = qJ(4) * t315;
t336 = pkin(8) * t313 - t407;
t268 = t336 + t413;
t410 = pkin(4) - t317;
t285 = t410 * t315;
t312 = sin(qJ(5));
t314 = cos(qJ(5));
t385 = t314 * t268 + t312 * t285;
t316 = -pkin(3) - pkin(8);
t375 = qJD(3) * t316;
t412 = qJD(4) - t281;
t363 = MDP(24) - MDP(27);
t381 = qJD(1) * t313;
t346 = qJD(5) * t381;
t366 = qJD(1) * qJD(3);
t347 = t315 * t366;
t365 = qJD(3) * qJD(5);
t243 = -t314 * t346 + (-t347 + t365) * t312;
t290 = t314 * t347;
t354 = t312 * t381;
t377 = qJD(3) * t314;
t274 = t354 + t377;
t373 = qJD(5) * t274;
t244 = -t290 + t373;
t246 = -pkin(4) * t347 - t262;
t216 = pkin(5) * t244 + qJ(6) * t243 - qJD(6) * t274 + t246;
t406 = t216 * t312;
t405 = t216 * t314;
t404 = t243 * t314;
t403 = t244 * t312;
t402 = t246 * t312;
t401 = t246 * t314;
t379 = qJD(3) * t312;
t272 = -t314 * t381 + t379;
t398 = t272 * t298;
t397 = t272 * t312;
t396 = t272 * t313;
t395 = t272 * t314;
t394 = t274 * t312;
t393 = t274 * t314;
t392 = t298 * t315;
t391 = t298 * t316;
t390 = t314 * t315;
t319 = qJD(1) ^ 2;
t389 = t315 * t319;
t318 = qJD(3) ^ 2;
t388 = t317 * t318;
t338 = pkin(5) * t314 + qJ(6) * t312;
t327 = -pkin(4) - t338;
t322 = t327 * t315;
t387 = qJD(1) * t322 - qJD(5) * t338 + qJD(6) * t314 - t412;
t278 = pkin(3) * t380 + qJ(4) * t381;
t258 = pkin(8) * t380 + t278;
t259 = -pkin(4) * t381 + t280;
t386 = t314 * t258 + t312 * t259;
t384 = pkin(3) * t381 + qJD(1) * qJ(2);
t382 = t318 + t319;
t378 = qJD(3) * t313;
t374 = qJD(4) * t315;
t372 = qJD(5) * t312;
t371 = qJD(5) * t314;
t370 = qJD(6) * t298;
t251 = t259 + t308;
t224 = pkin(5) * t272 - qJ(6) * t274 + t251;
t369 = t224 * MDP(28);
t368 = pkin(4) * t380 + t412;
t240 = t368 + t375;
t249 = qJD(1) * t336 + t384;
t222 = t240 * t314 - t249 * t312;
t367 = qJD(6) - t222;
t307 = qJD(1) * qJD(2);
t362 = 0.2e1 * qJD(1);
t360 = t312 * t391;
t359 = t314 * t391;
t358 = t313 * t389;
t335 = (qJD(3) * pkin(8) - qJD(4)) * t315;
t348 = t313 * t366;
t355 = pkin(3) * t347 + qJ(4) * t348 + t307;
t233 = qJD(1) * t335 + t355;
t276 = t292 * t378;
t253 = -pkin(4) * t348 + t276;
t357 = -t314 * t233 - t240 * t371 - t312 * t253;
t353 = t313 * t375;
t352 = t314 * t375;
t351 = t298 * t372;
t350 = t313 * t371;
t349 = pkin(3) * t376 + qJ(4) * t378 + qJD(2);
t265 = -qJ(4) * t380 + t384;
t283 = -t407 + t413;
t343 = qJD(1) * t283 + t265;
t342 = t312 * t233 + t240 * t372 + t249 * t371 - t314 * t253;
t341 = qJ(6) * t348;
t324 = t249 * t372 + t357;
t214 = -t324 - t341 + t370;
t219 = -pkin(5) * t298 + t367;
t340 = -t219 * t380 - t214;
t337 = pkin(5) * t312 - qJ(6) * t314;
t296 = pkin(5) * t348;
t215 = t296 + t342;
t334 = -t214 * t312 + t215 * t314;
t223 = t240 * t312 + t249 * t314;
t220 = qJ(6) * t298 + t223;
t333 = t219 * t314 - t220 * t312;
t332 = t219 * t312 + t220 * t314;
t331 = -t258 * t312 + t259 * t314;
t330 = -t268 * t312 + t285 * t314;
t329 = -qJD(1) * t310 + t392;
t328 = t274 * t298;
t242 = -qJD(1) * t374 + t355;
t256 = t349 - t374;
t326 = -qJD(1) * t256 - t242 + t388;
t325 = t223 * t298 - t342;
t245 = t335 + t349;
t270 = t410 * t378;
t323 = t314 * t245 - t268 * t372 - t312 * t270 + t285 * t371;
t321 = t222 * t298 + t324;
t320 = (-t394 + t395) * MDP(26) - t332 * MDP(28) + (t312 * t364 + t314 * t363) * t298;
t301 = t313 * t317;
t295 = t317 * t376;
t288 = t312 * t348;
t284 = -pkin(4) * t313 + t301;
t282 = qJ(4) + t337;
t271 = -pkin(4) * t376 + t295;
t254 = t265 * t380;
t241 = t313 * t327 + t301;
t234 = pkin(5) * t274 + qJ(6) * t272;
t229 = -pkin(5) * t315 - t330;
t228 = qJ(6) * t315 + t385;
t227 = -t243 + t398;
t226 = pkin(5) * t381 - t331;
t225 = -qJ(6) * t381 + t386;
t221 = t295 + (qJD(5) * t337 - qJD(6) * t312) * t313 + qJD(3) * t322;
t218 = pkin(5) * t378 + qJD(5) * t385 + t245 * t312 + t270 * t314;
t217 = -qJ(6) * t378 + qJD(6) * t315 + t323;
t1 = [-0.2e1 * t313 * MDP(7) * t347 + 0.2e1 * t366 * t415 + (-t313 * t388 + (qJ(2) * t376 + qJD(2) * t313) * t362) * MDP(12) + (-t315 * t388 + (-qJ(2) * t378 + qJD(2) * t315) * t362) * MDP(13) + t417 * MDP(14) + (t313 * t326 - t343 * t376) * MDP(15) + (t315 * t326 + t343 * t378) * MDP(16) + (t242 * t283 + t256 * t265 - t317 * t417) * MDP(17) + (-t243 * t312 * t313 + (t312 * t376 + t350) * t274) * MDP(18) + ((t393 - t397) * t376 + (-t404 - t403 + (-t394 - t395) * qJD(5)) * t313) * MDP(19) + (t298 * t350 - t243 * t315 + (-t274 * t313 + t312 * t329) * qJD(3)) * MDP(20) + (-t313 * t351 - t244 * t315 + (t314 * t329 + t396) * qJD(3)) * MDP(21) + (-t298 - t380) * MDP(22) * t378 + (-t342 * t315 + t271 * t272 + t284 * t244 + (-qJD(5) * t285 - t245) * t298 * t312 + ((-qJD(5) * t268 - t270) * t298 - t251 * t376) * t314 + (t251 * t372 - t401 + (-qJD(1) * t330 - t222) * qJD(3)) * t313) * MDP(23) + (-t323 * t298 + t271 * t274 - t284 * t243 + ((qJD(3) * t251 + qJD(5) * t249) * t312 + t357) * t315 + (t251 * t371 + t402 + (qJD(1) * t385 + t223) * qJD(3)) * t313) * MDP(24) + (-t218 * t298 + t221 * t272 + t241 * t244 + (-t224 * t377 - t215) * t315 + (t224 * t372 - t405 + (qJD(1) * t229 + t219) * qJD(3)) * t313) * MDP(25) + (-t217 * t272 + t218 * t274 - t228 * t244 - t229 * t243 + t332 * t376 + (qJD(5) * t333 + t214 * t314 + t215 * t312) * t313) * MDP(26) + (t217 * t298 - t221 * t274 + t241 * t243 + (-t224 * t379 + t214) * t315 + (-t224 * t371 - t406 + (-qJD(1) * t228 - t220) * qJD(3)) * t313) * MDP(27) + (t214 * t228 + t215 * t229 + t216 * t241 + t217 * t220 + t218 * t219 + t221 * t224) * MDP(28) + 0.2e1 * t414 * t307 + (-MDP(10) * t315 - MDP(9) * t313) * t318; -t414 * t319 + (-t265 * MDP(17) + t320) * qJD(1) + t364 * (t272 * t376 + (t244 + t290) * t313) + ((-MDP(13) + MDP(16)) * t382 + (t403 - t404) * MDP(26) + t334 * MDP(28) + t320 * qJD(5) + ((-t266 - t280) * MDP(17) + t369 - t363 * (-t274 + t354)) * qJD(3)) * t315 + ((-MDP(12) + MDP(15)) * t382 - t262 * MDP(17) + t216 * MDP(28) - t363 * t243 + (t264 * MDP(17) + (-t393 - t397) * MDP(26) - t333 * MDP(28) + (-t312 * t363 + t314 * t364) * t298) * qJD(3)) * t313; MDP(7) * t358 - t319 * t415 + ((-t266 - t308) * t315 + (t264 + t345) * t313) * qJD(1) * MDP(14) + (t278 * t381 + t254) * MDP(15) + (0.2e1 * t306 + (-t265 * t313 + t278 * t315) * qJD(1)) * MDP(16) + (-qJ(4) * t262 - qJD(4) * t266 - t265 * t278 + (t399 + (-t264 - t409) * t313) * t292) * MDP(17) + (-t312 * t328 - t404) * MDP(18) + ((-t244 - t328) * t314 + (t243 + t398) * t312) * MDP(19) + (-t351 + (-t312 * t392 + (t274 - t377) * t313) * qJD(1)) * MDP(20) + (-t298 * t371 + t288 + (-t298 * t390 - t396) * qJD(1)) * MDP(21) + t298 * MDP(22) * t381 + (qJ(4) * t244 + t402 - t331 * t298 + t368 * t272 + (t251 * t314 - t360) * qJD(5) + (t251 * t390 + (t222 - t352) * t313) * qJD(1)) * MDP(23) + (-qJ(4) * t243 + t401 + t386 * t298 + t368 * t274 + (-t251 * t312 - t359) * qJD(5) + (-t223 * t313 + (-t251 * t315 + t353) * t312) * qJD(1)) * MDP(24) + (t406 + t226 * t298 + t244 * t282 - t387 * t272 + (t224 * t314 - t360) * qJD(5) + (t224 * t390 + (-t219 - t352) * t313) * qJD(1)) * MDP(25) + (t225 * t272 - t226 * t274 + (-t220 * t380 + t243 * t316 + t215 + (-t272 * t316 - t220) * qJD(5)) * t314 + (-t244 * t316 + (t274 * t316 - t219) * qJD(5) + t340) * t312) * MDP(26) + (-t405 - t225 * t298 + t243 * t282 + t387 * t274 + (t224 * t312 + t359) * qJD(5) + (t220 * t313 + (t224 * t315 - t353) * t312) * qJD(1)) * MDP(27) + (t216 * t282 - t219 * t226 - t220 * t225 - t387 * t224 + (qJD(5) * t332 - t334) * t316) * MDP(28) + (MDP(13) * t313 * t319 - MDP(12) * t389) * qJ(2); -MDP(15) * t358 + (-t311 * t319 - t318) * MDP(16) + (t254 + t276) * MDP(17) + t363 * t288 + (t266 * MDP(17) - t272 * t364 - t274 * t363 - t369) * qJD(3) + ((-qJD(5) * t272 + t243) * MDP(26) + (qJD(5) * t220 - t215) * MDP(28) + ((-MDP(26) * t272 + MDP(28) * t220) * t315 - t364 * t378) * qJD(1) - t363 * t416) * t314 + ((t274 * t380 - t244 + t373) * MDP(26) + (qJD(5) * t219 - t340) * MDP(28) - t364 * t416) * t312; t227 * MDP(20) + (-t312 * t346 - t314 * t365 + t290) * MDP(21) - MDP(22) * t348 + t325 * MDP(23) + t321 * MDP(24) + (-0.2e1 * t296 + t325) * MDP(25) + (pkin(5) * t243 - qJ(6) * t244) * MDP(26) + (-t321 - 0.2e1 * t341 + 0.2e1 * t370) * MDP(27) + (-pkin(5) * t215 + qJ(6) * t214 - t219 * t223 + t220 * t367 - t224 * t234) * MDP(28) + (t298 * MDP(21) - t251 * MDP(23) - t224 * MDP(25) + (t220 - t223) * MDP(26) + t234 * MDP(27) + MDP(19) * t274) * t274 + (t274 * MDP(18) + t251 * MDP(24) - t234 * MDP(25) + (t219 - t367) * MDP(26) - t224 * MDP(27) - MDP(19) * t272) * t272; (t272 * t274 + t348) * MDP(25) + t227 * MDP(26) + (-t274 ^ 2 - t416) * MDP(27) + (-t220 * t298 + t224 * t274 + t215) * MDP(28);];
tauc  = t1;
