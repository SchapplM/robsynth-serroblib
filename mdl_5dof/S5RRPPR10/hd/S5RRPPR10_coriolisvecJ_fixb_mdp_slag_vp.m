% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:47
% EndTime: 2019-12-31 19:44:53
% DurationCPUTime: 2.96s
% Computational Cost: add. (1492->361), mult. (3813->490), div. (0->0), fcn. (2439->6), ass. (0->167)
t370 = (qJD(1) * qJD(2));
t412 = -2 * t370;
t324 = sin(qJ(2));
t411 = MDP(4) * t324;
t326 = cos(qJ(2));
t320 = t326 ^ 2;
t410 = MDP(5) * (t324 ^ 2 - t320);
t374 = qJD(4) * t326;
t380 = qJD(2) * t324;
t408 = qJ(4) * t380 - t374;
t321 = sin(pkin(8));
t383 = qJD(1) * t324;
t363 = t321 * t383;
t322 = cos(pkin(8));
t371 = t322 * qJD(2);
t281 = t363 - t371;
t361 = t322 * t383;
t381 = qJD(2) * t321;
t283 = t361 + t381;
t323 = sin(qJ(5));
t325 = cos(qJ(5));
t238 = t281 * t323 + t283 * t325;
t358 = t326 * t370;
t304 = t321 * t358;
t349 = t322 * t358;
t332 = -t325 * t304 + t323 * t349;
t214 = t238 * qJD(5) + t332;
t279 = t283 ^ 2;
t407 = pkin(3) + pkin(4);
t406 = -pkin(7) + qJ(3);
t405 = qJD(2) * pkin(2);
t404 = qJ(3) * t324;
t403 = qJ(4) * t322;
t296 = -pkin(2) * t326 - pkin(1) - t404;
t275 = t296 * qJD(1);
t382 = qJD(1) * t326;
t315 = pkin(6) * t382;
t302 = qJD(2) * qJ(3) + t315;
t240 = t275 * t322 - t321 * t302;
t231 = pkin(3) * t382 + qJD(4) - t240;
t402 = t231 * t324;
t241 = t321 * t275 + t322 * t302;
t233 = -qJ(4) * t382 + t241;
t401 = t233 * t324;
t397 = t283 * t323;
t236 = -t325 * t281 + t397;
t311 = qJD(5) + t382;
t400 = t236 * t311;
t399 = t238 * t311;
t345 = pkin(2) * t324 - qJ(3) * t326;
t269 = t345 * qJD(2) - qJD(3) * t324;
t398 = t269 * t322;
t289 = t345 * qJD(1);
t396 = t289 * t322;
t395 = t321 * t324;
t394 = t321 * t325;
t393 = t321 * t326;
t392 = t322 * t324;
t391 = t322 * t326;
t327 = qJD(2) ^ 2;
t390 = t324 * t327;
t389 = t326 * t327;
t328 = qJD(1) ^ 2;
t388 = t326 * t328;
t285 = t321 * t323 + t322 * t325;
t335 = t285 * t326;
t387 = -qJD(1) * t335 - t285 * qJD(5);
t362 = t321 * t382;
t367 = t323 * t391;
t372 = qJD(5) * t325;
t373 = qJD(5) * t323;
t386 = -qJD(1) * t367 + t321 * t372 - t322 * t373 + t325 * t362;
t260 = t269 * qJD(1);
t314 = pkin(6) * t383;
t293 = (qJD(3) - t314) * qJD(2);
t229 = t321 * t260 + t322 * t293;
t310 = pkin(6) * t358;
t385 = -pkin(3) * t304 - t310;
t309 = pkin(6) * t391;
t257 = t321 * t296 + t309;
t379 = qJD(2) * t326;
t378 = qJD(3) * t283;
t377 = qJD(3) * t322;
t376 = qJD(4) * t283;
t375 = qJD(4) * t321;
t369 = pkin(7) * t391;
t308 = pkin(6) * t393;
t368 = pkin(6) * t380;
t359 = t324 * t370;
t366 = qJ(4) * t359 + t229;
t365 = t281 * t372 + t323 * t304 + t325 * t349;
t364 = -pkin(6) * t321 - pkin(3);
t360 = qJD(4) * t392;
t357 = MDP(23) * t383;
t356 = qJ(4) * t321 + pkin(2);
t353 = -qJD(3) + t405;
t347 = -t314 + t353;
t355 = t347 - t405;
t354 = pkin(1) * t412;
t228 = t260 * t322 - t321 * t293;
t210 = (-t407 * t324 - t369) * t370 - t228;
t211 = (pkin(7) * t381 - qJD(4)) * t382 + t366;
t352 = t325 * t210 - t211 * t323;
t256 = t296 * t322 - t308;
t339 = -t407 * t321 + t403;
t351 = -t339 * t382 + t315 + t375;
t344 = pkin(3) * t321 - t403;
t255 = t344 * t382 + t315;
t350 = t255 + t375;
t348 = t364 * t324;
t250 = -qJ(4) * t326 + t257;
t343 = t210 * t323 + t211 * t325;
t212 = pkin(4) * t382 - pkin(7) * t283 + t231;
t217 = pkin(7) * t281 + t233;
t207 = t212 * t325 - t217 * t323;
t208 = t212 * t323 + t217 * t325;
t318 = t326 * pkin(3);
t234 = pkin(4) * t326 + t308 + t318 + (-pkin(7) * t324 - t296) * t322;
t239 = pkin(7) * t395 + t250;
t342 = t234 * t325 - t239 * t323;
t341 = t234 * t323 + t239 * t325;
t286 = -t322 * t323 + t394;
t273 = t321 * t289;
t253 = -pkin(6) * t361 + t273;
t261 = t321 * t269;
t247 = -t322 * t368 + t261;
t340 = pkin(6) + t344;
t300 = t406 * t322;
t329 = -t369 + (-pkin(4) + t364) * t324;
t338 = -t329 * qJD(1) + qJD(3) * t321 - qJD(5) * t300 + t396;
t299 = t406 * t321;
t312 = qJ(4) * t383;
t336 = -pkin(6) * t392 + pkin(7) * t393;
t337 = t336 * qJD(1) - qJD(5) * t299 + t273 + t312 - t377;
t265 = t285 * t324;
t334 = t283 * t373 - t365;
t333 = -pkin(6) + t339;
t331 = qJ(4) * t283 + t347;
t330 = -qJ(4) * t349 - t385;
t303 = qJD(3) * t362;
t294 = -pkin(3) * t322 - t356;
t274 = t407 * t322 + t356;
t267 = t281 * t382;
t266 = t281 * t377;
t264 = t323 * t392 - t324 * t394;
t262 = t340 * t324;
t252 = pkin(6) * t363 + t396;
t251 = -t256 + t318;
t249 = t333 * t324;
t246 = t321 * t368 + t398;
t245 = qJD(1) * t348 - t396;
t244 = t253 + t312;
t243 = t340 * t379 - t360;
t235 = qJD(2) * t348 - t398;
t230 = t333 * t379 + t360;
t227 = t247 + t408;
t226 = pkin(3) * t281 - t331;
t225 = t286 * t324 * qJD(5) + qJD(2) * t335;
t224 = qJD(2) * t367 + qJD(5) * t265 - t379 * t394;
t223 = t330 - t376;
t221 = -pkin(3) * t359 - t228;
t220 = t336 * qJD(2) + t261 + t408;
t219 = t329 * qJD(2) - t398;
t218 = t376 + (-pkin(4) * t321 + t403) * t358 + t385;
t216 = -qJD(1) * t374 + t366;
t215 = -t407 * t281 + t331;
t1 = [0.2e1 * t358 * t411 + t410 * t412 + MDP(6) * t389 - MDP(7) * t390 + (-pkin(6) * t389 + t324 * t354) * MDP(9) + (pkin(6) * t390 + t326 * t354) * MDP(10) + ((-qJD(1) * t246 - t228) * t326 + ((pkin(6) * t281 - t321 * t347) * t326 + (t240 + (t256 + 0.2e1 * t308) * qJD(1)) * t324) * qJD(2)) * MDP(11) + ((qJD(1) * t247 + t229) * t326 + ((pkin(6) * t283 - t322 * t347) * t326 + (-t241 + (-t257 + 0.2e1 * t309) * qJD(1)) * t324) * qJD(2)) * MDP(12) + (-t246 * t283 - t247 * t281 + (-t228 * t322 - t229 * t321) * t324 + (-t240 * t322 - t241 * t321 + (-t256 * t322 - t257 * t321) * qJD(1)) * t379) * MDP(13) + (t228 * t256 + t229 * t257 + t240 * t246 + t241 * t247 + (-t347 + t314) * pkin(6) * t379) * MDP(14) + (t223 * t395 + t243 * t281 + (qJD(1) * t235 + t221) * t326 + (t226 * t393 - t402 + (-t251 * t324 + t262 * t393) * qJD(1)) * qJD(2)) * MDP(15) + (-t227 * t281 + t235 * t283 + (-t216 * t321 + t221 * t322) * t324 + (t231 * t322 - t233 * t321 + (-t250 * t321 + t251 * t322) * qJD(1)) * t379) * MDP(16) + (-t223 * t392 - t243 * t283 + (-qJD(1) * t227 - t216) * t326 + (-t226 * t391 + t401 + (t250 * t324 - t262 * t391) * qJD(1)) * qJD(2)) * MDP(17) + (t216 * t250 + t221 * t251 + t223 * t262 + t226 * t243 + t227 * t233 + t231 * t235) * MDP(18) + (t225 * t238 - t265 * t334) * MDP(19) + (-t214 * t265 - t224 * t238 - t225 * t236 + t264 * t334) * MDP(20) + (-t334 * t326 + t225 * t311 + (-qJD(1) * t265 - t238) * t380) * MDP(21) + (-t214 * t326 - t224 * t311 + (qJD(1) * t264 + t236) * t380) * MDP(22) + (-t311 - t382) * MDP(23) * t380 + ((t219 * t325 - t220 * t323) * t311 + t352 * t326 + t230 * t236 + t249 * t214 + t218 * t264 + t215 * t224 + (-t208 * t326 - t341 * t311) * qJD(5) + (-t342 * qJD(1) - t207) * t380) * MDP(24) + (-(t219 * t323 + t220 * t325) * t311 - t343 * t326 + t230 * t238 - t249 * t334 + t218 * t265 + t215 * t225 + (-t207 * t326 - t342 * t311) * qJD(5) + (t341 * qJD(1) + t208) * t380) * MDP(25); -t388 * t411 + t328 * t410 + (t303 + ((-qJ(3) * t381 - t240) * t324 + (t252 + t355 * t321 + (-t281 - t371) * pkin(6)) * t326) * qJD(1)) * MDP(11) + ((-qJ(3) * t371 + t241) * t324 + (-t253 + (-t283 + t381) * pkin(6) + (t347 - t353) * t322) * t326) * qJD(1) * MDP(12) + (t252 * t283 + t253 * t281 - t266 + (t240 * t382 + t229) * t322 + (t241 * t382 - t228 + t378) * t321) * MDP(13) + (-t240 * t252 - t241 * t253 + (-t240 * t321 + t241 * t322) * qJD(3) + (-t228 * t321 + t229 * t322) * qJ(3) + t355 * t315) * MDP(14) + (-t223 * t322 + t303 - t350 * t281 + (t402 - t245 * t326 + (-t226 * t326 + (t294 * t326 - t404) * qJD(2)) * t321) * qJD(1)) * MDP(15) + (t244 * t281 - t245 * t283 - t266 + (-t231 * t382 + t216) * t322 + (t233 * t382 + t221 + t378) * t321) * MDP(16) + (-t223 * t321 + t350 * t283 + (-t401 + t244 * t326 + (qJ(3) * t380 + (-qJD(2) * t294 - qJD(3) + t226) * t326) * t322) * qJD(1)) * MDP(17) + (qJ(3) * t216 * t322 + t223 * t294 - t226 * t255 - t231 * t245 + (-t244 + t377) * t233 + (qJ(3) * t221 + qJD(3) * t231 - qJD(4) * t226) * t321) * MDP(18) + (t387 * t238 - t286 * t334) * MDP(19) + (-t214 * t286 - t387 * t236 - t386 * t238 + t285 * t334) * MDP(20) + (t387 * t311 + (-qJD(2) * t286 + t238) * t383) * MDP(21) + (-t386 * t311 + (qJD(2) * t285 - t236) * t383) * MDP(22) + t311 * t357 + (t274 * t214 + t218 * t285 + (t337 * t323 + t338 * t325) * t311 + t351 * t236 + t386 * t215 + (-(t299 * t325 - t300 * t323) * qJD(2) + t207) * t383) * MDP(24) + (-t274 * t334 + t218 * t286 + (-t338 * t323 + t337 * t325) * t311 + t351 * t238 + t387 * t215 + ((t299 * t323 + t300 * t325) * qJD(2) - t208) * t383) * MDP(25) + (t328 * t324 * MDP(9) + MDP(10) * t388) * pkin(1); (t240 * t283 + t241 * t281 + t310) * MDP(14) + (t233 * t281 + (-qJD(4) - t231) * t283 + t330) * MDP(18) + (-t214 - t399) * MDP(24) + (t334 + t400) * MDP(25) + (MDP(11) + MDP(15)) * (-t283 * t382 + t304) + (MDP(12) - MDP(17)) * (t267 + t349) + (MDP(13) + MDP(16)) * (-t281 ^ 2 - t279); t283 * t281 * MDP(15) - t267 * MDP(16) + (-t320 * t328 - t279) * MDP(17) + (t226 * t283 - t228) * MDP(18) + (-t236 * t283 - t311 * t373) * MDP(24) + (-t238 * t283 - t311 * t372) * MDP(25) + ((t233 * MDP(18) + (-MDP(24) * t323 - MDP(25) * t325) * t311) * t326 + (MDP(16) * t391 + (-MDP(18) * pkin(3) - MDP(24) * t325 + MDP(25) * t323 - MDP(15)) * t324) * qJD(2)) * qJD(1); t238 * t236 * MDP(19) + (-t236 ^ 2 + t238 ^ 2) * MDP(20) + (t365 + t400) * MDP(21) + (-t332 + t399) * MDP(22) - qJD(2) * t357 + (t208 * t311 - t215 * t238 + t352) * MDP(24) + (t207 * t311 + t215 * t236 - t343) * MDP(25) + (-MDP(21) * t397 - t238 * MDP(22) - t208 * MDP(24) - t207 * MDP(25)) * qJD(5);];
tauc = t1;
