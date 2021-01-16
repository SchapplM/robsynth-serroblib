% Calculate Coriolis joint torque vector for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:51:58
% EndTime: 2021-01-16 01:52:06
% DurationCPUTime: 2.72s
% Computational Cost: add. (2204->359), mult. (4876->489), div. (0->0), fcn. (3221->8), ass. (0->163)
t298 = cos(qJ(4));
t291 = t298 ^ 2;
t295 = sin(qJ(4));
t405 = (t295 ^ 2 - t291) * MDP(9);
t300 = -pkin(2) - pkin(8);
t299 = cos(qJ(2));
t292 = sin(pkin(6));
t358 = qJD(1) * t292;
t339 = t299 * t358;
t313 = qJD(3) - t339;
t260 = t300 * qJD(2) + t313;
t293 = cos(pkin(6));
t376 = t293 * t298;
t285 = qJD(1) * t376;
t241 = t295 * t260 + t285;
t233 = qJD(4) * pkin(9) + t241;
t277 = pkin(4) * t295 - pkin(9) * t298 + qJ(3);
t296 = sin(qJ(2));
t340 = t296 * t358;
t251 = t277 * qJD(2) + t340;
t294 = sin(qJ(5));
t297 = cos(qJ(5));
t219 = -t233 * t294 + t297 * t251;
t353 = qJD(2) * t298;
t334 = t297 * t353;
t352 = qJD(4) * t294;
t273 = t334 + t352;
t213 = -qJ(6) * t273 + t219;
t354 = qJD(2) * t295;
t286 = qJD(5) + t354;
t210 = pkin(5) * t286 + t213;
t404 = t210 - t213;
t345 = t297 * qJD(4);
t287 = qJD(5) * t345;
t347 = qJD(5) * t298;
t331 = t294 * t347;
t333 = t295 * t345;
t306 = t331 + t333;
t244 = t306 * qJD(2) - t287;
t335 = t294 * t353;
t271 = t335 - t345;
t403 = -t271 * t286 - t244;
t344 = qJD(2) * qJD(4);
t329 = t295 * t344;
t282 = t294 * t329;
t245 = qJD(5) * t273 - t282;
t402 = t273 * t286 + t245;
t315 = pkin(4) * t298 + pkin(9) * t295;
t267 = t315 * qJD(4) + qJD(3);
t374 = t294 * t296;
t401 = -(-t295 * t374 + t297 * t299) * t358 + t297 * t267;
t348 = qJD(5) * t297;
t369 = t296 * t297;
t400 = (t294 * t299 + t295 * t369) * t358 - t298 * t300 * t345 - t294 * t267 - t277 * t348;
t370 = t295 * t300;
t361 = t294 * t277 + t297 * t370;
t357 = qJD(1) * t295;
t240 = t260 * t298 - t293 * t357;
t343 = MDP(20) + MDP(22);
t399 = pkin(5) * t271;
t398 = -qJ(6) - pkin(9);
t397 = qJD(2) * pkin(2);
t387 = t251 * t294;
t220 = t233 * t297 + t387;
t214 = -qJ(6) * t271 + t220;
t396 = t214 * t286;
t355 = qJD(2) * t292;
t338 = t296 * t355;
t316 = t298 * t338;
t351 = qJD(4) * t295;
t362 = -qJD(4) * t285 - t260 * t351;
t225 = -qJD(1) * t316 - t362;
t217 = pkin(5) * t245 + t225;
t395 = t217 * t294;
t394 = t217 * t297;
t393 = t225 * t294;
t392 = t225 * t297;
t232 = -qJD(4) * pkin(4) - t240;
t391 = t232 * t294;
t390 = t232 * t297;
t389 = t244 * t294;
t388 = t245 * t297;
t384 = t271 * t294;
t383 = t271 * t297;
t381 = t273 * t294;
t380 = t273 * t297;
t379 = t286 * t294;
t378 = t286 * t297;
t377 = t292 * t299;
t375 = t294 * t295;
t373 = t294 * t298;
t372 = t294 * t300;
t371 = t295 * t297;
t368 = t297 * t298;
t326 = pkin(5) - t372;
t346 = qJD(6) * t297;
t349 = qJD(5) * t294;
t367 = qJ(6) * t333 - t361 * qJD(5) + (qJ(6) * t349 + t326 * qJD(4) - t346) * t298 + t401;
t330 = t297 * t347;
t366 = -qJ(6) * t330 + (-qJD(6) * t298 + (qJ(6) * qJD(4) - qJD(5) * t300) * t295) * t294 - t400;
t275 = t315 * qJD(2);
t322 = -t240 * t294 + t297 * t275;
t325 = qJD(5) * t398;
t365 = (pkin(5) * t298 + qJ(6) * t371) * qJD(2) + t322 + qJD(6) * t294 - t297 * t325;
t336 = t294 * t354;
t363 = t297 * t240 + t294 * t275;
t364 = qJ(6) * t336 - t294 * t325 - t346 + t363;
t301 = qJD(4) ^ 2;
t302 = qJD(2) ^ 2;
t359 = -t301 - t302;
t356 = qJD(2) * qJ(3);
t350 = qJD(4) * t298;
t342 = MDP(21) + MDP(23);
t337 = t299 * t355;
t332 = t286 * t349;
t328 = t298 * t344;
t327 = MDP(19) * t350;
t324 = -qJD(6) - t399;
t224 = t260 * t350 + (-qJD(4) * t293 + t338) * t357;
t242 = (t267 + t339) * qJD(2);
t323 = -t224 * t294 + t297 * t242;
t321 = t286 + t354;
t320 = qJD(5) * t295 + qJD(2);
t319 = pkin(5) * t328;
t318 = -t297 * t224 + t233 * t349 - t294 * t242 - t251 * t348;
t317 = t298 * t340;
t276 = t340 + t356;
t314 = -t276 + t340;
t312 = t210 * t297 + t214 * t294;
t311 = t210 * t294 - t214 * t297;
t310 = qJD(2) * t291 - t286 * t295;
t309 = -pkin(9) * t350 + t232 * t295;
t259 = -t295 * t377 + t376;
t237 = -t259 * t294 + t292 * t369;
t238 = t259 * t297 + t292 * t374;
t258 = t293 * t295 + t298 * t377;
t308 = qJ(6) * t244 + t323;
t307 = -qJ(6) * t245 - t318;
t305 = t314 - t356;
t268 = (qJD(3) + t339) * qJD(2);
t304 = t313 * qJD(2) - t300 * t301 + t268;
t303 = (t380 + t384) * MDP(24) - t312 * MDP(25);
t288 = -pkin(5) * t297 - pkin(4);
t280 = t398 * t297;
t279 = t398 * t294;
t270 = t313 - t397;
t269 = (pkin(5) * t294 - t300) * t298;
t266 = t271 ^ 2;
t265 = t297 * t277;
t250 = t273 * t317;
t249 = t271 * t317;
t246 = t300 * t351 + (-t294 * t351 + t330) * pkin(5);
t236 = t259 * qJD(4) - t316;
t235 = -t258 * qJD(4) + t295 * t338;
t234 = -qJ(6) * t373 + t361;
t228 = -qJ(6) * t368 + t326 * t295 + t265;
t226 = -pkin(5) * t336 + t241;
t222 = t232 - t324;
t216 = t237 * qJD(5) + t235 * t297 + t294 * t337;
t215 = -t238 * qJD(5) - t235 * t294 + t297 * t337;
t207 = -qJD(6) * t271 + t307;
t206 = -t220 * qJD(5) - qJD(6) * t273 + t308 + t319;
t1 = [(-t215 * t273 - t216 * t271 + t237 * t244 - t238 * t245) * MDP(24) + (t206 * t237 + t207 * t238 + t210 * t215 + t214 * t216 + t217 * t258 + t222 * t236) * MDP(25) + t342 * (-t216 * t286 + t236 * t273 - t238 * t328 - t244 * t258) + t343 * (t215 * t286 + t236 * t271 + t237 * t328 + t245 * t258) + (-MDP(13) * t236 - MDP(14) * t235) * qJD(4) + ((qJD(2) * t276 * MDP(7) + (MDP(13) * t295 + MDP(14) * t298 - MDP(4) + MDP(6)) * t302) * t299 + (t268 * MDP(7) + (-MDP(3) + MDP(5)) * t302 + ((MDP(13) * t298 - MDP(14) * t295) * qJD(4) + (t270 - t339) * MDP(7)) * qJD(2)) * t296) * t292; 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (qJ(3) * t268 + qJD(3) * t276 + (-t276 * t299 + (-t270 - t397) * t296) * t358) * MDP(7) + 0.2e1 * t344 * t405 - t301 * t298 * MDP(11) - t305 * t350 * MDP(13) + (t304 * t298 + t305 * t351) * MDP(14) + (-t244 * t368 - t306 * t273) * MDP(15) + ((t381 + t383) * t351 + (t389 - t388 + (-t380 + t384) * qJD(5)) * t298) * MDP(16) + (-t286 * t331 + (t273 * t298 + t310 * t297) * qJD(4)) * MDP(17) + (-t286 * t330 + (-t271 * t298 - t310 * t294) * qJD(4)) * MDP(18) + t321 * t327 + (t249 + (-t277 * t349 + t401) * t286 + (t232 * t348 + t393 - t300 * t245 + (-t286 * t372 + (-t294 * t370 + t265) * qJD(2) + t219) * qJD(4)) * t298) * MDP(20) + (t250 + t400 * t286 + (-t232 * t349 + t392 + t300 * t244 + (-t361 * qJD(2) - t220) * qJD(4)) * t298) * MDP(21) + (t245 * t269 + t246 * t271 + t249 + t367 * t286 + (t222 * t348 + t395 + (qJD(2) * t228 + t210) * qJD(4)) * t298) * MDP(22) + (-t244 * t269 + t246 * t273 + t250 - t366 * t286 + (-t222 * t349 + t394 + (-qJD(2) * t234 - t214) * qJD(4)) * t298) * MDP(23) + (t228 * t244 - t234 * t245 - t367 * t273 - t366 * t271 + t312 * t351 + (t311 * qJD(5) - t206 * t297 - t207 * t294) * t298) * MDP(24) + (t206 * t228 + t207 * t234 + t217 * t269 + (t246 + t317) * t222 + t366 * t214 + t367 * t210) * MDP(25) + (-0.2e1 * MDP(8) * t328 - t301 * MDP(10) + t304 * MDP(13) - t244 * MDP(17) - t245 * MDP(18) + ((t271 * t300 - t391) * qJD(4) + (-t387 + (-t286 * t300 - t233) * t297) * qJD(5) + t323) * MDP(20) + (t300 * t332 + (t273 * t300 - t390) * qJD(4) + t318) * MDP(21) + (-t222 * t352 + t206) * MDP(22) + (-t222 * t345 - t207) * MDP(23)) * t295; -t302 * MDP(6) + t343 * (-t245 * t298 - t320 * t378 + (t271 * t295 - t321 * t373) * qJD(4)) + t342 * (t244 * t298 + t320 * t379 + (-t286 * t368 + (t273 - t334) * t295) * qJD(4)) + (t314 * MDP(7) + t303) * qJD(2) + (t359 * MDP(14) - t217 * MDP(25) + ((t381 - t383) * MDP(24) - t311 * MDP(25)) * qJD(4)) * t298 + (t359 * MDP(13) + (-t388 - t389) * MDP(24) + (qJD(4) * t222 - t206 * t294 + t207 * t297) * MDP(25) + t303 * qJD(5)) * t295; (qJD(4) * t241 + t314 * t353 + t362) * MDP(13) - t314 * t354 * MDP(14) + (t273 * t378 - t389) * MDP(15) + (-t294 * t402 + t297 * t403) * MDP(16) + (t286 * t348 + (t286 * t371 + (-t273 + t352) * t298) * qJD(2)) * MDP(17) + (-t332 + (-t286 * t375 + (t271 + t345) * t298) * qJD(2)) * MDP(18) - t286 * MDP(19) * t353 + (-pkin(4) * t245 - t392 - t322 * t286 - t241 * t271 + (-pkin(9) * t378 + t391) * qJD(5) + (-t219 * t298 + t309 * t294) * qJD(2)) * MDP(20) + (pkin(4) * t244 + t393 + t363 * t286 - t241 * t273 + (pkin(9) * t379 + t390) * qJD(5) + (t220 * t298 + t309 * t297) * qJD(2)) * MDP(21) + (-t394 - t226 * t271 + t245 * t288 - t365 * t286 + (t222 + t399) * t349 + (t222 * t375 + (qJD(4) * t279 - t210) * t298) * qJD(2)) * MDP(22) + (t395 - t226 * t273 - t244 * t288 + t364 * t286 + (pkin(5) * t381 + t222 * t297) * qJD(5) + (t222 * t371 + (qJD(4) * t280 + t214) * t298) * qJD(2)) * MDP(23) + (t244 * t279 + t245 * t280 + t365 * t273 + t364 * t271 + (-t286 * t210 + t207) * t297 + (-t206 - t396) * t294) * MDP(24) + (t206 * t279 - t207 * t280 + t217 * t288 + (pkin(5) * t349 - t226) * t222 - t364 * t214 - t365 * t210) * MDP(25) + (t298 * t295 * MDP(8) - t405) * t302; -t266 * MDP(16) + (-t297 * t329 + t287) * MDP(17) + t282 * MDP(18) + qJD(2) * t327 + (t220 * t286 + t323) * MDP(20) + (t219 * t286 + t318) * MDP(21) + (t308 + 0.2e1 * t319 + t396) * MDP(22) + (t213 * t286 - t307) * MDP(23) + pkin(5) * t244 * MDP(24) + (pkin(5) * t206 + t214 * t404) * MDP(25) + (t286 * MDP(17) + t232 * MDP(21) + (qJD(6) + t222) * MDP(23) - t404 * MDP(24)) * t271 + (-MDP(17) * t335 - t273 * MDP(18) - t220 * t343) * qJD(5) + (t271 * MDP(15) + t286 * MDP(18) - t232 * MDP(20) + (-t222 + t324) * MDP(22) - pkin(5) * t222 * MDP(25) + (-MDP(23) * pkin(5) + MDP(16)) * t273) * t273; t402 * MDP(22) + t403 * MDP(23) + (-t273 ^ 2 - t266) * MDP(24) + (t210 * t273 + t214 * t271 + t217) * MDP(25);];
tauc = t1;
