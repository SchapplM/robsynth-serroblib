% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:24
% EndTime: 2019-03-09 02:06:30
% DurationCPUTime: 2.80s
% Computational Cost: add. (2589->348), mult. (4962->477), div. (0->0), fcn. (2800->6), ass. (0->144)
t304 = cos(qJ(4));
t302 = sin(qJ(4));
t303 = cos(qJ(5));
t356 = qJD(5) * t303;
t333 = t302 * t356;
t301 = sin(qJ(5));
t352 = t301 * qJD(4);
t402 = t304 * t352 + t333;
t365 = qJD(1) * t302;
t364 = qJD(1) * t304;
t286 = qJD(5) + t364;
t401 = qJ(2) * MDP(6) + MDP(5);
t349 = qJD(1) * qJD(4);
t329 = t302 * t349;
t285 = pkin(5) * t329;
t299 = cos(pkin(9));
t350 = qJD(1) * qJD(2);
t331 = t299 * t350;
t305 = -pkin(1) - pkin(2);
t284 = t305 * qJD(1) + qJD(2);
t298 = sin(pkin(9));
t366 = qJ(2) * qJD(1);
t263 = t298 * t284 + t299 * t366;
t253 = -qJD(1) * pkin(7) + t263;
t398 = qJD(3) * t304 - t302 * t253;
t226 = qJD(4) * t398 + t304 * t331;
t295 = t302 * qJD(3);
t238 = t304 * t253 + t295;
t233 = qJD(4) * pkin(8) + t238;
t262 = t284 * t299 - t298 * t366;
t252 = qJD(1) * pkin(3) - t262;
t324 = pkin(4) * t304 + pkin(8) * t302;
t234 = t324 * qJD(1) + t252;
t323 = -pkin(4) * t302 + pkin(8) * t304;
t267 = t298 * qJD(2) + t323 * qJD(4);
t254 = t267 * qJD(1);
t357 = qJD(5) * t301;
t325 = -t301 * t226 - t233 * t356 - t234 * t357 + t303 * t254;
t208 = t285 - t325;
t215 = t233 * t303 + t234 * t301;
t213 = qJ(6) * t286 + t215;
t400 = -t213 * t286 + t208;
t380 = t301 * t304;
t268 = t298 * t380 + t299 * t303;
t361 = qJD(4) * t302;
t336 = t298 * t361;
t378 = t303 * t304;
t374 = -t268 * qJD(5) - t303 * t336 - (t298 * t301 + t299 * t378) * qJD(1);
t269 = t298 * t378 - t299 * t301;
t399 = t269 * qJD(5) - t301 * t336 - (-t298 * t303 + t299 * t380) * qJD(1);
t397 = MDP(10) * t302;
t296 = t302 ^ 2;
t396 = (-t304 ^ 2 + t296) * MDP(11);
t348 = MDP(22) + MDP(24);
t347 = MDP(23) - MDP(26);
t359 = qJD(4) * t304;
t227 = qJD(4) * t295 + t253 * t359 + t302 * t331;
t328 = t304 * t349;
t360 = qJD(4) * t303;
t271 = t301 * t365 + t360;
t358 = qJD(5) * t271;
t249 = -t303 * t328 + t358;
t289 = qJD(5) * t352;
t250 = t402 * qJD(1) - t289;
t272 = t303 * t365 - t352;
t209 = -pkin(5) * t250 - qJ(6) * t249 + qJD(6) * t272 + t227;
t395 = -t209 * MDP(27) - t347 * t249;
t232 = -qJD(4) * pkin(4) - t398;
t217 = -pkin(5) * t271 + qJ(6) * t272 + t232;
t394 = t217 * MDP(27) - t347 * t272;
t392 = t209 * t301;
t391 = t209 * t303;
t389 = t227 * t301;
t388 = t227 * t303;
t387 = t232 * t301;
t386 = t232 * t303;
t385 = t249 * t301;
t327 = -t298 * qJ(2) + t299 * t305;
t273 = pkin(3) - t327;
t255 = t273 + t324;
t384 = t255 * t303;
t383 = t271 * t286;
t382 = t286 * t301;
t381 = t286 * t303;
t379 = t302 * t303;
t321 = pkin(5) * t301 - qJ(6) * t303;
t377 = qJD(6) * t301 - t286 * t321 + t238;
t275 = t323 * qJD(1);
t376 = t301 * t275 + t303 * t398;
t370 = t299 * qJ(2) + t298 * t305;
t274 = -pkin(7) + t370;
t372 = t301 * t255 + t274 * t378;
t306 = qJD(4) ^ 2;
t307 = qJD(1) ^ 2;
t368 = t306 + t307;
t367 = MDP(25) * t303;
t363 = qJD(2) * t299;
t355 = qJD(6) * t286;
t214 = -t233 * t301 + t234 * t303;
t351 = qJD(6) - t214;
t346 = pkin(8) * t382;
t345 = pkin(8) * t381;
t344 = pkin(8) * t360;
t343 = 0.2e1 * t350;
t342 = t286 * t380;
t341 = t286 * t378;
t340 = -t303 * t226 - t234 * t356 - t301 * t254;
t339 = -t402 * t272 + t302 * t385;
t337 = t304 * t363;
t338 = t255 * t356 + t301 * t267 + t303 * t337;
t334 = t302 * t357;
t332 = t286 * t356;
t330 = t296 * t349;
t326 = t298 * t343;
t322 = -pkin(5) * t303 - qJ(6) * t301;
t312 = t233 * t357 + t340;
t207 = -qJ(6) * t329 - t312 + t355;
t320 = t207 * t303 + t208 * t301;
t212 = -pkin(5) * t286 + t351;
t319 = t212 * t303 - t213 * t301;
t318 = t212 * t301 + t213 * t303;
t317 = t275 * t303 - t301 * t398;
t315 = t274 - t321;
t314 = -t274 * t306 + t326;
t313 = t215 * t286 + t325;
t311 = t303 * t359 - t334;
t310 = qJD(4) * (-qJD(1) * t273 - t252 - t363);
t309 = t394 * qJD(4);
t308 = t214 * t286 + t312;
t282 = t303 * t330;
t281 = t301 * t330;
t279 = t301 * pkin(8) * t329;
t277 = -pkin(4) + t322;
t240 = -pkin(5) * t272 - qJ(6) * t271;
t235 = t315 * t302;
t225 = t249 - t383;
t221 = -t384 + (t274 * t301 - pkin(5)) * t304;
t220 = qJ(6) * t304 + t372;
t219 = pkin(5) * t365 - t317;
t218 = -qJ(6) * t365 + t376;
t216 = t315 * t359 + (t322 * qJD(5) + qJD(6) * t303 + t363) * t302;
t211 = pkin(5) * t361 + (qJD(5) * t274 * t304 - t267) * t303 + (qJD(5) * t255 - t274 * t361 + t337) * t301;
t210 = (-t274 * t357 + qJD(6)) * t304 + (-t274 * t303 - qJ(6)) * t361 + t338;
t1 = [MDP(7) * t326 + 0.2e1 * MDP(8) * t331 + (-t262 * t298 + t263 * t299 + (-t298 * t327 + t299 * t370) * qJD(1)) * qJD(2) * MDP(9) + 0.2e1 * t328 * t397 - 0.2e1 * t349 * t396 + (t302 * t310 + t314 * t304) * MDP(15) + (-t314 * t302 + t304 * t310) * MDP(16) + (-t249 * t379 + t311 * t272) * MDP(17) + (-t250 * t379 - t311 * t271 + t339) * MDP(18) + (t286 * t334 + t249 * t304 + t282 + (t272 * t302 - t341) * qJD(4)) * MDP(19) + (t302 * t332 + t250 * t304 - t281 + (-t271 * t302 + t342) * qJD(4)) * MDP(20) + (-t286 - t364) * MDP(21) * t361 + ((-t255 * t357 + t267 * t303) * t286 + ((-t274 * t356 - t301 * t363) * t286 + (-t271 * t274 - t387) * qJD(4) + t325) * t304 + (-t271 * t363 - t232 * t356 - t389 - t274 * t250 + (t274 * t382 - (-t274 * t380 + t384) * qJD(1) - t214) * qJD(4)) * t302) * MDP(22) + (-t338 * t286 + ((t274 * t286 + t233) * t357 + (-t272 * t274 - t386) * qJD(4) + t340) * t304 + (-t272 * t363 + t232 * t357 - t388 + t274 * t249 + (t372 * qJD(1) + t274 * t381 + t215) * qJD(4)) * t302) * MDP(23) + (-t211 * t286 - t216 * t271 - t235 * t250 + (-t217 * t352 - t208) * t304 + (-t217 * t356 - t392 + (qJD(1) * t221 + t212) * qJD(4)) * t302) * MDP(24) + (t210 * t271 - t211 * t272 + t220 * t250 + t221 * t249 - t319 * t359 + (t318 * qJD(5) + t207 * t301 - t208 * t303) * t302) * MDP(25) + (t210 * t286 + t216 * t272 - t235 * t249 + (t217 * t360 + t207) * t304 + (-t217 * t357 + t391 + (-qJD(1) * t220 - t213) * qJD(4)) * t302) * MDP(26) + (t207 * t220 + t208 * t221 + t209 * t235 + t210 * t213 + t211 * t212 + t216 * t217) * MDP(27) + t401 * t343 + (-t304 * MDP(12) + t302 * MDP(13)) * t306; (t249 * t268 + t250 * t269 + t374 * t271 - t272 * t399) * MDP(25) + (t207 * t269 + t208 * t268 + t399 * t212 + t374 * t213) * MDP(27) - t401 * t307 + (-t307 * MDP(8) + (0.2e1 * MDP(16) * t359 - MDP(9) * t263) * qJD(1)) * t299 + (-t347 * t374 - t348 * t399) * t286 + ((t348 * t268 + t347 * t269) * qJD(4) + (0.2e1 * qJD(4) * MDP(15) + t348 * t271 - t394) * t299) * t365 + (qJD(1) * t262 * MDP(9) - t307 * MDP(7) + (-t368 * MDP(15) + t309) * t304 + (t368 * MDP(16) - t395) * t302 + t348 * (-t250 * t302 - t271 * t359)) * t298; t339 * MDP(25) + t347 * t282 + (-t306 * MDP(16) + t348 * t250 + (t271 * t367 + t318 * MDP(27) + (-t348 * t301 - t347 * t303) * t286) * qJD(4) + t395) * t304 + (-t306 * MDP(15) + t250 * t367 + t320 * MDP(27) + t309 + (-t271 * t301 * MDP(25) + t319 * MDP(27) + (t347 * t301 - t348 * t303) * t286) * qJD(5)) * t302 + t348 * (-t271 * t361 + t281); (qJD(4) * t238 + t252 * t365 - t227) * MDP(15) + (t252 - t363) * t364 * MDP(16) + (-t272 * t381 + t385) * MDP(17) + ((t249 + t383) * t303 + (t272 * t286 + t250) * t301) * MDP(18) + (t332 + (t341 + (-t272 - t352) * t302) * qJD(1)) * MDP(19) + (-t286 * t357 + (-t342 + (t271 - t360) * t302) * qJD(1)) * MDP(20) + t286 * MDP(21) * t365 + (t279 + pkin(4) * t250 - t388 - t317 * t286 + t238 * t271 + (-t345 + t387) * qJD(5) + (t214 * t302 + t232 * t380) * qJD(1)) * MDP(22) + (-pkin(4) * t249 + t389 + t376 * t286 + t238 * t272 + (t346 + t386) * qJD(5) + (t232 * t378 + (-t215 + t344) * t302) * qJD(1)) * MDP(23) + (-t391 + t219 * t286 - t277 * t250 + t279 + t377 * t271 + (t217 * t301 - t345) * qJD(5) + (-t212 * t302 + t217 * t380) * qJD(1)) * MDP(24) + (-t218 * t271 + t219 * t272 + (t207 + t286 * t212 + (-qJD(5) * t272 + t250) * pkin(8)) * t303 + ((t249 - t358) * pkin(8) + t400) * t301) * MDP(25) + (-t392 - t218 * t286 - t277 * t249 - t377 * t272 + (-t217 * t303 - t346) * qJD(5) + (-t217 * t378 + (t213 - t344) * t302) * qJD(1)) * MDP(26) + (t209 * t277 - t212 * t219 - t213 * t218 - t377 * t217 + (t319 * qJD(5) + t320) * pkin(8)) * MDP(27) + (-t304 * t397 + t396) * t307; t225 * MDP(19) - t289 * MDP(20) + t313 * MDP(22) + t308 * MDP(23) + (-0.2e1 * t285 + t313) * MDP(24) + (-pkin(5) * t249 + qJ(6) * t250) * MDP(25) + (-t308 + 0.2e1 * t355) * MDP(26) + (-pkin(5) * t208 + qJ(6) * t207 - t212 * t215 + t351 * t213 - t217 * t240) * MDP(27) + (-t286 * MDP(20) + t232 * MDP(22) + t217 * MDP(24) + (-t213 + t215) * MDP(25) - t240 * MDP(26) + MDP(18) * t272) * t272 + (MDP(20) * t333 + (MDP(20) * t380 + (-0.2e1 * qJ(6) * MDP(26) - MDP(21)) * t302) * qJD(4)) * qJD(1) + (t272 * MDP(17) - t232 * MDP(23) + t240 * MDP(24) + (-t212 + t351) * MDP(25) + t217 * MDP(26) - MDP(18) * t271) * t271; (t271 * t272 + t329) * MDP(24) + t225 * MDP(25) + (-t272 ^ 2 - t286 ^ 2) * MDP(26) + (-t217 * t272 + t400) * MDP(27);];
tauc  = t1;
