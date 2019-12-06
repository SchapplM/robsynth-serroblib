% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR3
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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:43:02
% EndTime: 2019-12-05 18:43:07
% DurationCPUTime: 1.82s
% Computational Cost: add. (2044->235), mult. (3576->322), div. (0->0), fcn. (2480->8), ass. (0->139)
t351 = sin(pkin(9));
t352 = cos(pkin(9));
t354 = sin(qJ(3));
t357 = cos(qJ(3));
t324 = -t351 * t354 + t352 * t357;
t348 = qJD(1) + qJD(2);
t306 = t324 * t348;
t356 = cos(qJ(5));
t292 = t356 * t306;
t325 = t351 * t357 + t352 * t354;
t317 = t325 * qJD(3);
t298 = t348 * t317;
t391 = qJD(3) * t357;
t378 = t348 * t391;
t392 = qJD(3) * t354;
t379 = t348 * t392;
t299 = -t351 * t379 + t352 * t378;
t307 = t325 * t348;
t353 = sin(qJ(5));
t389 = qJD(5) * t353;
t221 = qJD(5) * t292 - t353 * t298 + t356 * t299 - t307 * t389;
t366 = t306 * t353 + t356 * t307;
t222 = t366 * qJD(5) + t356 * t298 + t299 * t353;
t253 = -t307 * t353 + t292;
t347 = qJD(3) + qJD(5);
t405 = t253 * t347;
t406 = t366 * t347;
t425 = (-t222 + t406) * MDP(19) + (-t253 ^ 2 + t366 ^ 2) * MDP(17) - t253 * MDP(16) * t366 + (t221 - t405) * MDP(18);
t424 = qJ(4) + pkin(7);
t355 = sin(qJ(2));
t408 = pkin(1) * qJD(1);
t384 = t355 * t408;
t376 = t424 * t348 + t384;
t294 = t376 * t354;
t285 = qJD(3) * pkin(3) - t294;
t295 = t376 * t357;
t403 = t352 * t295;
t245 = t351 * t285 + t403;
t413 = pkin(8) * t306;
t231 = t245 + t413;
t380 = -pkin(3) * t357 - pkin(2);
t358 = cos(qJ(2));
t383 = t358 * t408;
t305 = t380 * t348 + qJD(4) - t383;
t260 = -pkin(4) * t306 + t305;
t423 = t231 * t389 - t260 * t253;
t420 = MDP(7) * t354;
t419 = (t354 ^ 2 - t357 ^ 2) * MDP(8);
t344 = t357 * qJD(4);
t377 = qJD(3) * t424;
t314 = -t354 * t377 + t344;
t315 = -qJD(4) * t354 - t357 * t377;
t397 = -t314 * t351 + t352 * t315 + t325 * t383;
t396 = t352 * t314 + t351 * t315 - t324 * t383;
t418 = t354 * MDP(12) + t357 * MDP(13);
t417 = qJD(5) - t347;
t407 = pkin(1) * qJD(2);
t381 = qJD(1) * t407;
t372 = t358 * t381;
t362 = qJD(4) * t348 + t372;
t364 = qJD(3) * t376;
t258 = -t354 * t364 + t362 * t357;
t259 = -t362 * t354 - t357 * t364;
t223 = -t258 * t351 + t352 * t259;
t215 = -pkin(8) * t299 + t223;
t224 = t352 * t258 + t351 * t259;
t216 = -pkin(8) * t298 + t224;
t416 = t356 * t215 - t353 * t216 - t260 * t366;
t415 = pkin(1) * t358;
t414 = pkin(3) * t351;
t412 = pkin(8) * t307;
t318 = t324 * qJD(3);
t411 = pkin(8) * t318;
t410 = pkin(8) * t325;
t404 = t348 * t354;
t279 = t351 * t295;
t359 = qJD(3) ^ 2;
t402 = t354 * t359;
t401 = t357 * t359;
t340 = pkin(1) * t355 + pkin(7);
t400 = -qJ(4) - t340;
t365 = t324 * t356 - t325 * t353;
t236 = t365 * qJD(5) - t317 * t353 + t318 * t356;
t337 = t355 * t381;
t316 = pkin(3) * t379 + t337;
t265 = pkin(4) * t298 + t316;
t269 = t324 * t353 + t325 * t356;
t399 = t260 * t236 + t265 * t269;
t237 = t269 * qJD(5) + t356 * t317 + t318 * t353;
t398 = t260 * t237 - t265 * t365;
t373 = qJD(3) * t400;
t382 = t358 * t407;
t274 = t354 * t373 + t357 * t382 + t344;
t275 = (-qJD(4) - t382) * t354 + t357 * t373;
t241 = t352 * t274 + t351 * t275;
t247 = -t352 * t294 - t279;
t322 = t400 * t354;
t345 = t357 * qJ(4);
t323 = t340 * t357 + t345;
t267 = t351 * t322 + t352 * t323;
t331 = -pkin(2) * t348 - t383;
t395 = t331 * t391 + t354 * t337;
t334 = t424 * t354;
t335 = pkin(7) * t357 + t345;
t284 = -t351 * t334 + t352 * t335;
t390 = qJD(3) * t358;
t387 = t357 * MDP(12);
t385 = -qJD(2) + t348;
t342 = pkin(3) * t392;
t291 = pkin(4) * t317 + t342;
t240 = -t274 * t351 + t352 * t275;
t244 = t352 * t285 - t279;
t246 = t294 * t351 - t403;
t266 = t352 * t322 - t323 * t351;
t283 = -t352 * t334 - t335 * t351;
t371 = -t223 * t325 + t224 * t324 - t244 * t318 - t245 * t317;
t370 = t291 - t384;
t321 = t324 * pkin(8);
t369 = qJD(5) * (t321 + t284) + t411 - t397;
t313 = t317 * pkin(8);
t368 = -qJD(5) * (t283 - t410) + t313 - t396;
t228 = qJD(3) * pkin(4) + t244 - t412;
t367 = -t353 * t228 - t356 * t231;
t301 = -pkin(4) * t324 + t380;
t361 = -MDP(10) * t402 + (t221 * t365 - t222 * t269 + t236 * t253 - t237 * t366) * MDP(17) + (t221 * t269 + t236 * t366) * MDP(16) - 0.2e1 * t348 * qJD(3) * t419 + 0.2e1 * t378 * t420 + MDP(9) * t401 + (t236 * MDP(18) - t237 * MDP(19)) * t347;
t343 = t355 * t407;
t341 = -pkin(2) - t415;
t338 = pkin(3) * t352 + pkin(4);
t319 = t331 * t392;
t290 = t301 - t415;
t278 = t291 + t343;
t273 = pkin(3) * t404 + pkin(4) * t307;
t249 = t321 + t267;
t248 = t266 - t410;
t233 = t247 - t412;
t232 = t246 - t413;
t226 = -t313 + t241;
t225 = t240 - t411;
t1 = [(-t240 * t307 + t241 * t306 - t266 * t299 - t267 * t298 + t371) * MDP(14) - t337 * MDP(5) + (t278 * t366 + t290 * t221 - (t225 * t353 + t226 * t356 + (t248 * t356 - t249 * t353) * qJD(5)) * t347 + t399) * MDP(22) + (-t278 * t253 + t290 * t222 + (t225 * t356 - t226 * t353 + (-t248 * t353 - t249 * t356) * qJD(5)) * t347 + t398) * MDP(21) + (-t340 * t401 + t341 * t379 + t319) * MDP(12) + (t340 * t402 + t341 * t378 + t395) * MDP(13) + (((-qJD(1) - t348) * MDP(6) - t418 * qJD(3)) * t358 + (-qJD(1) * t387 + (t354 * MDP(13) - MDP(5) - t387) * t348) * t355) * t407 + t361 + (t224 * t267 + t245 * t241 + t223 * t266 + t244 * t240 + t316 * (t380 - t415) + t305 * (t343 + t342)) * MDP(15); (-t283 * t299 - t284 * t298 + t396 * t306 - t397 * t307 + t371) * MDP(14) + (-pkin(2) * t378 + pkin(7) * t402 + (-t355 * t404 + t357 * t390) * t408 + t395) * MDP(13) + (-pkin(2) * t379 - pkin(7) * t401 + t319 + (t385 * t357 * t355 + t354 * t390) * t408) * MDP(12) + (t224 * t284 + t223 * t283 + t316 * t380 + (t342 - t384) * t305 + t396 * t245 + t397 * t244) * MDP(15) + (t301 * t221 + (t369 * t353 + t368 * t356) * t347 + t370 * t366 + t399) * MDP(22) + (t301 * t222 + (t368 * t353 - t369 * t356) * t347 - t370 * t253 + t398) * MDP(21) + t385 * MDP(6) * t383 + t361 + (t348 * t384 - t337) * MDP(5); ((t245 + t246) * t307 + (t244 - t247) * t306 + (-t298 * t351 - t299 * t352) * pkin(3)) * MDP(14) + (-t244 * t246 - t245 * t247 + (t223 * t352 + t224 * t351 - t305 * t404) * pkin(3)) * MDP(15) + (t273 * t253 - (t232 * t356 - t233 * t353) * t347 + ((-t338 * t353 - t356 * t414) * t347 + t367) * qJD(5) + t416) * MDP(21) + (-t356 * t216 - t353 * t215 - t273 * t366 + (t232 * t353 + t233 * t356) * t347 + (-(t338 * t356 - t353 * t414) * t347 - t356 * t228) * qJD(5) + t423) * MDP(22) + t418 * (-t331 * t348 - t372) + (-t357 * t420 + t419) * t348 ^ 2 + t425; (-t306 ^ 2 - t307 ^ 2) * MDP(14) + (t244 * t307 - t245 * t306 + t316) * MDP(15) + (t222 + t406) * MDP(21) + (t221 + t405) * MDP(22); (t417 * t367 + t416) * MDP(21) + ((-t231 * t347 - t215) * t353 + (-t417 * t228 - t216) * t356 + t423) * MDP(22) + t425;];
tauc = t1;
