% Calculate vector of inverse dynamics joint torques for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:27
% EndTime: 2019-12-05 16:46:30
% DurationCPUTime: 1.88s
% Computational Cost: add. (1419->255), mult. (2333->302), div. (0->0), fcn. (1562->10), ass. (0->144)
t285 = qJ(2) + qJ(3);
t278 = sin(t285);
t286 = sin(pkin(8));
t287 = cos(pkin(8));
t321 = g(1) * t287 + g(2) * t286;
t391 = t321 * t278;
t396 = 2 * qJD(4);
t291 = cos(qJ(4));
t395 = pkin(4) * t291 + pkin(3);
t272 = g(3) * t278;
t279 = cos(t285);
t394 = -t321 * t279 - t272;
t282 = qJD(2) + qJD(3);
t288 = sin(qJ(4));
t312 = qJ(5) * t288 + t395;
t393 = t282 * t312;
t281 = qJDD(2) + qJDD(3);
t392 = t312 * t281;
t290 = sin(qJ(2));
t293 = cos(qJ(2));
t340 = qJD(1) * qJD(2);
t305 = -qJDD(1) * t290 - t293 * t340;
t289 = sin(qJ(3));
t292 = cos(qJ(3));
t243 = t289 * t290 - t292 * t293;
t387 = t243 * t282;
t349 = qJD(1) * t290;
t329 = qJD(3) * t349;
t260 = t289 * t329;
t277 = t293 * qJDD(1);
t240 = qJDD(2) * pkin(2) - t290 * t340 + t277;
t348 = qJD(1) * t293;
t265 = qJD(2) * pkin(2) + t348;
t384 = -t292 * (qJD(3) * t265 - t305) - t240 * t289;
t208 = pkin(7) * t281 - t260 - t384;
t206 = t291 * t208;
t338 = qJDD(4) * qJ(5);
t230 = t265 * t289 + t292 * t349;
t223 = pkin(7) * t282 + t230;
t372 = t223 * t288;
t203 = t338 + t206 + (qJD(5) - t372) * qJD(4);
t205 = t288 * t208;
t342 = qJD(4) * t291;
t373 = qJDD(4) * pkin(4);
t385 = t223 * t342 - t373;
t204 = qJDD(5) + t205 + t385;
t386 = t203 * t291 + t204 * t288;
t383 = pkin(2) * t292;
t382 = pkin(3) * t281;
t381 = pkin(3) * t282;
t294 = qJD(4) ^ 2;
t379 = pkin(7) * t294;
t376 = g(3) * t279;
t375 = qJD(4) * pkin(4);
t374 = pkin(7) * qJDD(4);
t371 = t223 * t291;
t267 = t289 * t349;
t229 = t265 * t292 - t267;
t370 = t229 * t282;
t369 = t230 * t282;
t273 = pkin(2) * t289 + pkin(7);
t367 = t273 * t294;
t364 = t279 * t288;
t363 = t281 * t288;
t362 = t282 * t288;
t361 = t282 * t291;
t360 = t286 * t288;
t359 = t286 * t291;
t358 = t287 * t288;
t357 = t287 * t291;
t356 = t288 * t291;
t355 = qJDD(1) - g(3);
t341 = qJD(5) * t288;
t343 = qJD(4) * t288;
t236 = pkin(4) * t343 - qJ(5) * t342 - t341;
t347 = qJD(3) * t289;
t336 = pkin(2) * t347;
t225 = t236 + t336;
t244 = t289 * t293 + t290 * t292;
t238 = t244 * qJD(1);
t354 = t225 - t238;
t353 = t236 - t230;
t352 = t391 * t291;
t283 = t288 ^ 2;
t284 = t291 ^ 2;
t351 = t283 - t284;
t350 = t283 + t284;
t346 = qJD(3) * t292;
t345 = qJD(4) * qJ(5);
t344 = qJD(4) * t282;
t337 = qJDD(4) * t273;
t335 = pkin(2) * t346;
t280 = t282 ^ 2;
t334 = t280 * t356;
t333 = -g(3) * t364 + t391 * t288;
t331 = t282 * t347;
t317 = -t265 * t347 + (t240 - t329) * t292 + t305 * t289;
t209 = -t317 - t382;
t328 = -t209 - t376;
t327 = t350 * t281;
t326 = qJD(5) + t372;
t325 = t229 * t343 + t230 * t361 + t352;
t239 = t292 * t348 - t267;
t324 = t238 * t361 + t239 * t343 + t352;
t323 = t260 - t394;
t322 = t379 - t382;
t320 = pkin(4) * t288 - qJ(5) * t291;
t222 = -t229 - t381;
t319 = t209 * t288 + t222 * t342 - t333;
t318 = 0.2e1 * (t281 * t356 - t351 * t344) * MDP(9) + (t281 * t283 + 0.2e1 * t342 * t362) * MDP(8) + (qJDD(4) * t288 + t291 * t294) * MDP(10) + (qJDD(4) * t291 - t288 * t294) * MDP(11) + t281 * MDP(5);
t213 = t326 - t375;
t214 = t345 + t371;
t316 = t213 * t288 + t214 * t291;
t217 = t282 * t244;
t315 = -t217 * t282 - t243 * t281;
t274 = -pkin(3) - t383;
t314 = t274 * t281 + t367;
t313 = t350 * MDP(16) - MDP(7);
t233 = t279 * t359 - t358;
t235 = t279 * t357 + t360;
t311 = g(1) * t235 + g(2) * t233 - t206;
t201 = (qJD(4) * t320 - t341) * t282 - t392 - t317;
t310 = -t201 - t379 + t392;
t241 = -t312 - t383;
t309 = t241 * t282 - t335;
t308 = t274 * t282 - t335;
t307 = -t241 * t281 - t201 - t367;
t232 = t279 * t360 + t357;
t234 = t279 * t358 - t359;
t306 = g(1) * t234 + g(2) * t232 + t288 * t272 - t205;
t304 = t244 * t294 - t315;
t303 = -qJDD(5) + t306;
t302 = -qJDD(4) * t244 + t387 * t396;
t301 = -g(3) * t293 + t290 * t321;
t300 = t317 - t376 + t391;
t298 = t213 * t342 - t214 * t343 + t386 + t394;
t297 = (t213 * t291 - t214 * t288) * qJD(4) + t386;
t296 = -g(3) * (qJ(5) * t364 + t395 * t279) + t312 * t391 + t394 * pkin(7);
t295 = qJD(2) ^ 2;
t237 = t320 * t282;
t219 = t222 * t343;
t211 = -t229 - t393;
t210 = t211 * t343;
t1 = [t355 * MDP(1) + (qJDD(2) * t293 - t290 * t295) * MDP(3) + (-qJDD(2) * t290 - t293 * t295) * MDP(4) + t315 * MDP(6) + (t201 * t243 + t211 * t217 - g(3)) * MDP(18) + (MDP(13) + MDP(15)) * (t288 * t302 - t291 * t304) + (MDP(14) - MDP(17)) * (t288 * t304 + t291 * t302) - (MDP(18) * t316 + t282 * t313) * t387 + (MDP(18) * t297 + t281 * t313) * t244; qJDD(2) * MDP(2) + (t277 + t301) * MDP(3) + (-t355 * t290 + t321 * t293) * MDP(4) + (t238 * t282 + (t281 * t292 - t331) * pkin(2) + t300) * MDP(6) + (t239 * t282 + (-pkin(2) * t281 - t240) * t289 + ((-pkin(2) * t282 - t265) * qJD(3) + t305) * t292 + t323) * MDP(7) + (t219 + (qJD(4) * t308 - t337) * t288 + (-pkin(2) * t331 - t314 + t328) * t291 + t324) * MDP(13) + ((-t337 + (t239 + t308) * qJD(4)) * t291 + ((-t238 + t336) * t282 + t314) * t288 + t319) * MDP(14) + (t210 + (qJD(4) * t309 - t337) * t288 + (-t225 * t282 + t307 - t376) * t291 + t324) * MDP(15) + (t273 * t327 + t298 + (-t239 + t335) * t282 * t350) * MDP(16) + ((t337 + (-t211 - t239 - t309) * qJD(4)) * t291 + (-t354 * t282 + t307) * t288 + t333) * MDP(17) + (t201 * t241 - t316 * t239 + t354 * t211 + (t316 * t346 + t301) * pkin(2) + t297 * t273 + t296) * MDP(18) + t318; (t300 + t369) * MDP(6) + (t323 + t370 + t384) * MDP(7) + (t219 + (-pkin(3) * t344 - t374) * t288 + (-t322 + t328) * t291 + t325) * MDP(13) + ((-t374 + (t229 - t381) * qJD(4)) * t291 + (t322 - t369) * t288 + t319) * MDP(14) + (t210 + (-t312 * t344 - t374) * t288 + (-t236 * t282 + t310 - t376) * t291 + t325) * MDP(15) + (pkin(7) * t327 - t350 * t370 + t298) * MDP(16) + ((t374 + (-t211 - t229 + t393) * qJD(4)) * t291 + (-t353 * t282 + t310) * t288 + t333) * MDP(17) + (t297 * pkin(7) - t201 * t312 + t353 * t211 - t316 * t229 + t296) * MDP(18) + t318; -MDP(8) * t334 + t351 * MDP(9) * t280 + MDP(10) * t363 + t291 * t281 * MDP(11) + qJDD(4) * MDP(12) + (-t222 * t362 + t306) * MDP(13) + ((-t222 * t282 + t272) * t291 + t311) * MDP(14) + (0.2e1 * t373 + (-t211 * t288 + t237 * t291) * t282 + t303) * MDP(15) + (-t320 * t281 + ((t214 - t345) * t288 + (qJD(5) - t213 - t375) * t291) * t282) * MDP(16) + (-t291 * t272 + 0.2e1 * t338 + qJD(5) * t396 + (t211 * t291 + t237 * t288) * t282 - t311) * MDP(17) + (t203 * qJ(5) - t204 * pkin(4) - t211 * t237 - t213 * t371 - g(1) * (-pkin(4) * t234 + qJ(5) * t235) - g(2) * (-pkin(4) * t232 + qJ(5) * t233) + t320 * t272 + t326 * t214) * MDP(18); (-qJDD(4) - t334) * MDP(15) + MDP(16) * t363 + (-t280 * t283 - t294) * MDP(17) + (-qJD(4) * t214 + t211 * t362 - t303 + t385) * MDP(18);];
tau = t1;
