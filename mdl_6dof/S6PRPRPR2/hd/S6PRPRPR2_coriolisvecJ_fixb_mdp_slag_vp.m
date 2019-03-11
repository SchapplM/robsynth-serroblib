% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:33:00
% EndTime: 2019-03-08 19:33:07
% DurationCPUTime: 3.36s
% Computational Cost: add. (2207->343), mult. (5761->507), div. (0->0), fcn. (4523->12), ass. (0->158)
t342 = cos(pkin(11));
t346 = sin(qJ(2));
t340 = sin(pkin(6));
t399 = qJD(1) * t340;
t388 = t346 * t399;
t324 = t342 * t388;
t339 = sin(pkin(11));
t349 = cos(qJ(2));
t387 = t349 * t399;
t286 = t339 * t387 + t324;
t345 = sin(qJ(4));
t348 = cos(qJ(4));
t375 = pkin(4) * t345 - qJ(5) * t348;
t303 = qJD(4) * t375 - qJD(5) * t345;
t422 = t286 - t303;
t336 = t345 ^ 2;
t421 = MDP(7) * (-t348 ^ 2 + t336);
t295 = (t339 * t346 - t342 * t349) * t340;
t323 = t339 * t388;
t289 = t342 * t387 - t323;
t338 = sin(pkin(12));
t332 = pkin(2) * t339 + pkin(8);
t395 = qJD(4) * t345;
t385 = t332 * t395;
t311 = t338 * t385;
t341 = cos(pkin(12));
t411 = t338 * t348;
t407 = -t289 * t411 + t422 * t341 - t311;
t322 = qJD(2) * pkin(2) + t387;
t282 = t339 * t322 + t324;
t279 = qJD(2) * pkin(8) + t282;
t343 = cos(pkin(6));
t329 = qJD(1) * t343 + qJD(3);
t420 = -t345 * t279 + t329 * t348;
t409 = t341 * t348;
t419 = -t289 * t409 - t422 * t338;
t418 = pkin(2) * t342;
t417 = pkin(9) + qJ(5);
t288 = qJD(2) * t295;
t284 = qJD(1) * t288;
t394 = qJD(4) * t348;
t232 = t279 * t394 - t345 * t284 + t329 * t395;
t416 = t232 * t338;
t415 = t232 * t341;
t344 = sin(qJ(6));
t347 = cos(qJ(6));
t317 = t338 * t347 + t341 * t344;
t412 = t338 * t344;
t316 = -t347 * t341 + t412;
t357 = t316 * t348;
t392 = qJD(6) * t345;
t261 = -qJD(4) * t357 - t317 * t392;
t397 = qJD(2) * t348;
t330 = -qJD(6) + t397;
t414 = t261 * t330;
t410 = t341 * t345;
t230 = -t284 * t348 + (qJD(5) + t420) * qJD(4);
t296 = (t339 * t349 + t342 * t346) * t340;
t254 = (qJD(1) * t296 + t303) * qJD(2);
t221 = t341 * t230 + t338 * t254;
t364 = pkin(5) * t345 - pkin(9) * t409;
t356 = t364 * qJD(4);
t408 = -t356 + t407;
t406 = (-pkin(9) * t411 - t332 * t410) * qJD(4) + t419;
t378 = t341 * t385;
t405 = t378 - t419;
t258 = t348 * t279 + t345 * t329;
t253 = qJD(4) * qJ(5) + t258;
t281 = t322 * t342 - t323;
t363 = -pkin(4) * t348 - qJ(5) * t345 - pkin(3);
t263 = qJD(2) * t363 - t281;
t225 = t341 * t253 + t338 * t263;
t358 = t317 * t348;
t354 = qJD(4) * t358;
t391 = qJD(6) * t347;
t262 = t391 * t410 - t392 * t412 + t354;
t297 = t317 * t345;
t389 = qJD(2) * qJD(4);
t384 = t345 * t389;
t404 = t262 * t330 - t297 * t384;
t319 = t375 * qJD(2);
t235 = t338 * t319 + t341 * t420;
t403 = qJD(2) * t357 - t316 * qJD(6);
t402 = -qJD(2) * t358 + t317 * qJD(6);
t310 = t363 - t418;
t275 = t338 * t310 + t332 * t409;
t390 = t341 * qJD(4);
t398 = qJD(2) * t345;
t312 = t338 * t398 - t390;
t383 = t348 * t389;
t377 = t341 * t383;
t401 = -t312 * t391 + t347 * t377;
t396 = qJD(4) * t338;
t223 = -pkin(9) * t312 + t225;
t393 = qJD(6) * t223;
t386 = t338 * t397;
t382 = MDP(21) * t395;
t381 = pkin(5) * t338 + t332;
t380 = -qJD(4) * pkin(4) + qJD(5);
t220 = -t230 * t338 + t341 * t254;
t224 = -t253 * t338 + t341 * t263;
t234 = t341 * t319 - t338 * t420;
t376 = t338 * t383;
t219 = -pkin(9) * t376 + t221;
t314 = t341 * t398 + t396;
t222 = -pkin(5) * t397 - pkin(9) * t314 + t224;
t379 = -qJD(6) * t222 - t219;
t374 = -t220 * t338 + t221 * t341;
t215 = t222 * t347 - t223 * t344;
t216 = t222 * t344 + t223 * t347;
t373 = -t224 * t338 + t225 * t341;
t277 = t296 * t348 + t343 * t345;
t241 = -t277 * t338 + t295 * t341;
t242 = t277 * t341 + t295 * t338;
t372 = t241 * t347 - t242 * t344;
t371 = t241 * t344 + t242 * t347;
t300 = t341 * t310;
t260 = -pkin(9) * t410 + t300 + (-t332 * t338 - pkin(5)) * t348;
t265 = -pkin(9) * t338 * t345 + t275;
t370 = t260 * t347 - t265 * t344;
t369 = t260 * t344 + t265 * t347;
t276 = t296 * t345 - t343 * t348;
t367 = t312 * t344 - t314 * t347;
t365 = t338 * MDP(13) + t341 * MDP(14);
t287 = qJD(2) * t296;
t283 = qJD(1) * t287;
t350 = qJD(4) ^ 2;
t362 = qJD(2) * t286 - t332 * t350 - t283;
t278 = -qJD(2) * pkin(3) - t281;
t361 = qJD(4) * (qJD(2) * (-pkin(3) - t418) + t278 + t289);
t328 = t417 * t341;
t360 = qJD(2) * t364 + qJD(5) * t338 + qJD(6) * t328 + t234;
t327 = t417 * t338;
t359 = pkin(9) * t386 + qJD(5) * t341 - qJD(6) * t327 - t235;
t251 = -t420 + t380;
t355 = -qJD(6) * t314 - t376;
t353 = -qJ(5) * t395 + (-t251 + t380) * t348;
t304 = t347 * t312;
t268 = t314 * t344 + t304;
t352 = t312 * MDP(13) + t314 * MDP(14) + t251 * MDP(16) + t268 * MDP(22);
t238 = qJD(2) * t354 - qJD(6) * t367;
t351 = qJD(2) ^ 2;
t334 = -pkin(5) * t341 - pkin(4);
t302 = t381 * t345;
t298 = t316 * t345;
t294 = t381 * t394;
t274 = -t332 * t411 + t300;
t264 = t367 * t395;
t245 = pkin(5) * t386 + t258;
t240 = -qJD(4) * t276 - t288 * t348;
t237 = t344 * t355 + t401;
t236 = pkin(5) * t312 + t251;
t228 = pkin(5) * t376 + t232;
t227 = t240 * t341 + t287 * t338;
t226 = -t240 * t338 + t287 * t341;
t218 = qJD(2) * t356 + t220;
t217 = t347 * t218;
t1 = [(-t281 * t287 - t282 * t288 + t283 * t295 - t284 * t296) * MDP(5) - t240 * qJD(4) * MDP(12) + (-t226 * t314 - t227 * t312) * MDP(15) + (t220 * t241 + t221 * t242 + t224 * t226 + t225 * t227 + t232 * t276) * MDP(16) + (-(-qJD(6) * t371 + t226 * t347 - t227 * t344) * t330 + t276 * t238) * MDP(22) + ((qJD(6) * t372 + t226 * t344 + t227 * t347) * t330 + t276 * t237) * MDP(23) + (-MDP(3) * t346 - MDP(4) * t349) * t351 * t340 + (-MDP(11) * qJD(4) - MDP(23) * t367 + t352) * (qJD(4) * t277 - t288 * t345) + (t287 * t345 * MDP(12) + (-MDP(11) * t287 - MDP(13) * t226 + MDP(14) * t227) * t348 + ((t295 * MDP(12) + (-t241 * t341 - t242 * t338) * MDP(15) + t365 * t276) * t348 + (t295 * MDP(11) + t241 * MDP(13) - t242 * MDP(14) + MDP(22) * t372 - MDP(23) * t371) * t345) * qJD(4)) * qJD(2); (t281 * t286 - t282 * t289 + (-t283 * t342 - t284 * t339) * pkin(2)) * MDP(5) - 0.2e1 * t389 * t421 + (t407 * t314 + t405 * t312 + (-t224 * t341 - t225 * t338 + (-t274 * t341 - t275 * t338) * qJD(2)) * t394) * MDP(15) + (t332 * t394 * t251 + t220 * t274 + t221 * t275 - t407 * t224 - t405 * t225) * MDP(16) + (-t237 * t298 - t261 * t367) * MDP(17) + (-t237 * t297 + t238 * t298 - t261 * t268 + t262 * t367) * MDP(18) + (-t298 * t384 - t264 - t414) * MDP(19) + (-t268 * t395 + t404) * MDP(20) - t397 * t382 + (t228 * t297 + t236 * t262 + t302 * t238 + t294 * t268) * MDP(22) + (-t228 * t298 + t236 * t261 + t302 * t237 - t294 * t367) * MDP(23) + (-t382 + (qJD(6) * t369 + t344 * t406 + t347 * t408) * MDP(22) + (qJD(6) * t370 - t344 * t408 + t347 * t406) * MDP(23)) * t330 + (t350 * MDP(8) + t362 * MDP(11) + t361 * MDP(12) + (-t220 + (t251 * t338 + t312 * t332) * qJD(4) + (t311 + t407) * qJD(2)) * MDP(13) + (t221 + (t251 * t341 + t314 * t332) * qJD(4) + (t378 - t405) * qJD(2)) * MDP(14) - t237 * MDP(19) + t238 * MDP(20) + (qJD(6) * t216 + t219 * t344 - t217) * MDP(22) + (qJD(6) * t215 + t218 * t344 + t219 * t347) * MDP(23)) * t348 + (0.2e1 * MDP(6) * t383 - t350 * MDP(9) + t361 * MDP(11) - t362 * MDP(12) + (t416 - t289 * t312 + (qJD(2) * t274 + t224) * qJD(4)) * MDP(13) + (t415 - t289 * t314 + (-qJD(2) * t275 - t225) * qJD(4)) * MDP(14) + (-t220 * t341 - t221 * t338) * MDP(15) + (t232 * t332 - t251 * t289) * MDP(16) + (-t289 * t268 + (qJD(2) * t370 + t215) * qJD(4)) * MDP(22) + (t289 * t367 + (-qJD(2) * t369 - t216) * qJD(4)) * MDP(23)) * t345; t404 * MDP(22) + (-t264 + t414) * MDP(23) + (-t350 * MDP(12) - t232 * MDP(16) - t238 * MDP(22) - t237 * MDP(23)) * t348 + (-t350 * MDP(11) + MDP(16) * t374) * t345 + (-t365 * t336 * qJD(2) + (MDP(23) * qJD(2) * t298 + t352) * t345 + ((-t312 * t341 + t314 * t338) * MDP(15) + t373 * MDP(16)) * t348) * qJD(4); (qJD(4) * t258 - t278 * t398 - t232) * MDP(11) + (-qJD(2) * t278 + t284) * t348 * MDP(12) + (-t415 - t258 * t312 + (-t224 * t345 + t234 * t348 + t338 * t353) * qJD(2)) * MDP(13) + (t416 - t258 * t314 + (t225 * t345 - t235 * t348 + t341 * t353) * qJD(2)) * MDP(14) + (t234 * t314 + t235 * t312 + (-qJD(5) * t312 + t224 * t397 + t221) * t341 + (qJD(5) * t314 + t225 * t397 - t220) * t338) * MDP(15) + (-pkin(4) * t232 + qJ(5) * t374 + qJD(5) * t373 - t224 * t234 - t225 * t235 - t251 * t258) * MDP(16) + (t237 * t317 - t367 * t403) * MDP(17) + (-t237 * t316 - t238 * t317 - t268 * t403 + t367 * t402) * MDP(18) + (-t403 * t330 + (qJD(4) * t317 + t367) * t398) * MDP(19) + (t402 * t330 + (-qJD(4) * t316 + t268) * t398) * MDP(20) + t330 * MDP(21) * t398 + (t228 * t316 + t334 * t238 - t245 * t268 + (t344 * t359 + t347 * t360) * t330 + t402 * t236 + ((-t327 * t347 - t328 * t344) * qJD(4) - t215) * t398) * MDP(22) + (t228 * t317 + t334 * t237 + t245 * t367 + (-t344 * t360 + t347 * t359) * t330 + t403 * t236 + (-(-t327 * t344 + t328 * t347) * qJD(4) + t216) * t398) * MDP(23) + (-t345 * t348 * MDP(6) + t421) * t351; (-t312 ^ 2 - t314 ^ 2) * MDP(15) + (t224 * t314 + t225 * t312 + t232) * MDP(16) + (t330 * t367 + t238) * MDP(22) + (t304 * t330 + (-t376 + (-qJD(6) + t330) * t314) * t344 + t401) * MDP(23) + ((-t314 + t396) * MDP(13) + (t312 + t390) * MDP(14)) * t397; -t268 ^ 2 * MDP(18) + (-t268 * t330 + t401) * MDP(19) + qJD(2) * t382 + (-t216 * t330 + t217) * MDP(22) + (-t215 * t330 + t236 * t268) * MDP(23) - (MDP(17) * t268 - MDP(18) * t367 - MDP(20) * t330 - t236 * MDP(22)) * t367 + (MDP(20) * t355 - MDP(22) * t393 + MDP(23) * t379) * t347 + (t355 * MDP(19) + (qJD(6) * t312 - t377) * MDP(20) + t379 * MDP(22) + (-t218 + t393) * MDP(23)) * t344;];
tauc  = t1;
