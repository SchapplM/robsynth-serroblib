% Calculate vector of inverse dynamics joint torques for
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RRRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:31
% EndTime: 2019-12-31 21:49:35
% DurationCPUTime: 1.98s
% Computational Cost: add. (2213->306), mult. (3136->338), div. (0->0), fcn. (1662->12), ass. (0->163)
t300 = cos(qJ(2));
t401 = pkin(1) * t300;
t280 = qJDD(1) * t401;
t289 = qJDD(1) + qJDD(2);
t296 = sin(qJ(2));
t394 = pkin(1) * qJD(1);
t351 = t296 * t394;
t237 = pkin(2) * t289 - qJD(2) * t351 + t280;
t290 = qJD(1) + qJD(2);
t364 = qJD(1) * t300;
t350 = pkin(1) * t364;
t247 = pkin(2) * t290 + t350;
t295 = sin(qJ(3));
t334 = qJD(3) * t351;
t257 = t295 * t334;
t299 = cos(qJ(3));
t343 = qJD(2) * t364;
t355 = qJDD(1) * t296;
t310 = (t343 + t355) * pkin(1);
t405 = -t295 * t237 - (qJD(3) * t247 + t310) * t299 + t257;
t284 = qJDD(3) + t289;
t210 = pkin(8) * t284 - t405;
t298 = cos(qJ(4));
t208 = t298 * t210;
t354 = qJDD(4) * qJ(5);
t232 = t247 * t295 + t299 * t351;
t285 = qJD(3) + t290;
t226 = pkin(8) * t285 + t232;
t294 = sin(qJ(4));
t388 = t226 * t294;
t205 = t354 + t208 + (qJD(5) - t388) * qJD(4);
t207 = t294 * t210;
t357 = qJD(4) * t298;
t391 = qJDD(4) * pkin(4);
t402 = t226 * t357 - t391;
t206 = qJDD(5) + t207 + t402;
t403 = t205 * t298 + t206 * t294;
t293 = qJ(1) + qJ(2);
t288 = qJ(3) + t293;
t275 = sin(t288);
t276 = cos(t288);
t369 = g(1) * t276 + g(2) * t275;
t400 = pkin(2) * t299;
t399 = pkin(3) * t284;
t398 = pkin(3) * t285;
t397 = pkin(4) * t298;
t302 = qJD(4) ^ 2;
t396 = pkin(8) * t302;
t267 = g(1) * t275;
t286 = sin(t293);
t272 = g(1) * t286;
t395 = g(2) * t276;
t393 = qJD(4) * pkin(4);
t392 = pkin(8) * qJDD(4);
t279 = pkin(2) + t401;
t361 = qJD(3) * t299;
t362 = qJD(3) * t295;
t374 = t295 * t296;
t218 = t279 * t361 + (-t296 * t362 + (t299 * t300 - t374) * qJD(2)) * pkin(1);
t390 = t218 * t285;
t373 = t296 * t299;
t321 = t295 * t300 + t373;
t219 = t279 * t362 + (qJD(2) * t321 + t296 * t361) * pkin(1);
t389 = t219 * t285;
t387 = t226 * t298;
t261 = t295 * t351;
t231 = t247 * t299 - t261;
t386 = t231 * t285;
t385 = t232 * t285;
t368 = pkin(1) * t373 + t295 * t279;
t239 = pkin(8) + t368;
t384 = t239 * t302;
t328 = qJ(5) * t294 + t397;
t248 = -pkin(3) - t328;
t383 = t248 * t284;
t382 = t248 * t285;
t381 = t275 * t294;
t380 = t276 * t294;
t277 = pkin(2) * t295 + pkin(8);
t379 = t277 * t302;
t378 = t284 * t294;
t377 = t285 * t294;
t376 = t285 * t298;
t375 = t294 * t298;
t356 = qJD(5) * t294;
t358 = qJD(4) * t294;
t242 = pkin(4) * t358 - qJ(5) * t357 - t356;
t349 = pkin(2) * t362;
t233 = t242 + t349;
t240 = t321 * t394;
t372 = t233 - t240;
t371 = t242 - t232;
t370 = g(1) * t381 - g(2) * t380;
t287 = cos(t293);
t367 = g(1) * t287 + g(2) * t286;
t291 = t294 ^ 2;
t292 = t298 ^ 2;
t366 = t291 - t292;
t365 = t291 + t292;
t360 = qJD(4) * qJ(5);
t359 = qJD(4) * t285;
t353 = qJDD(4) * t239;
t352 = qJDD(4) * t277;
t269 = pkin(1) * t374;
t348 = pkin(2) * t361;
t283 = t285 ^ 2;
t347 = t283 * t375;
t254 = t298 * t267;
t346 = t231 * t358 + t232 * t376 + t254;
t241 = t299 * t350 - t261;
t345 = t240 * t376 + t241 * t358 + t254;
t344 = t285 * t362;
t325 = -t295 * pkin(1) * t343 - qJDD(1) * t269 - t247 * t362 + (t237 - t334) * t299;
t211 = -t325 - t399;
t342 = -t211 - t395;
t341 = t365 * t284;
t338 = t279 * t299 - t269;
t228 = t248 - t338;
t340 = t228 * t285 - t218;
t337 = qJD(5) + t388;
t336 = qJD(1) * (-qJD(2) + t290);
t335 = qJD(2) * (-qJD(1) - t290);
t225 = -t231 - t398;
t333 = t211 * t294 + t225 * t357 - t370;
t332 = t275 * pkin(8) + qJ(5) * t380 + (pkin(3) + t397) * t276;
t331 = -g(2) * t287 + t272 + t280;
t330 = t396 - t399;
t327 = pkin(4) * t294 - qJ(5) * t298;
t326 = 0.2e1 * (t284 * t375 - t359 * t366) * MDP(11) + (t284 * t291 + 0.2e1 * t357 * t377) * MDP(10) + (qJDD(4) * t294 + t298 * t302) * MDP(12) + (qJDD(4) * t298 - t294 * t302) * MDP(13) + t284 * MDP(7);
t324 = pkin(2) * t287 + t332;
t216 = t337 - t393;
t217 = t360 + t387;
t323 = t216 * t294 + t217 * t298;
t278 = -pkin(3) - t400;
t322 = t278 * t284 + t379;
t320 = g(1) * t380 + g(2) * t381 - g(3) * t298 - t207;
t202 = (qJD(4) * t327 - t356) * t285 + t383 - t325;
t319 = -t202 - t383 - t396;
t318 = t289 * MDP(4) + t326;
t244 = t248 - t400;
t317 = t244 * t285 - t348;
t316 = t278 * t285 - t348;
t315 = -t244 * t284 - t202 - t379;
t314 = -qJDD(5) + t320;
t238 = -pkin(3) - t338;
t312 = t238 * t284 + t384 + t389;
t311 = t248 * t267;
t309 = t267 + t325 - t395;
t308 = -t353 + (t238 * t285 - t218) * qJD(4);
t213 = t219 + t242;
t307 = -t213 * t285 - t228 * t284 - t202 - t384;
t306 = t216 * t357 - t217 * t358 - t369 + t403;
t263 = t276 * pkin(8);
t305 = -g(1) * t263 - t311;
t304 = (t216 * t298 - t217 * t294) * qJD(4) + t403;
t303 = t369 + t405;
t301 = cos(qJ(1));
t297 = sin(qJ(1));
t236 = t327 * t285;
t222 = t225 * t358;
t214 = -t231 + t382;
t212 = t214 * t358;
t1 = [t318 + qJDD(1) * MDP(1) + (((-qJDD(1) - t289) * t296 + t300 * t335) * pkin(1) + t367) * MDP(6) + ((t289 * t300 + t296 * t335) * pkin(1) + t331) * MDP(5) + (t202 * t228 + t214 * t213 - g(1) * (-pkin(1) * t297 - pkin(2) * t286 + t263) - g(2) * (pkin(1) * t301 + t324) - t311 + t323 * t218 + t304 * t239) * MDP(20) + (t222 + t254 + t308 * t294 + (-t312 + t342) * t298) * MDP(15) + (t294 * t312 + t298 * t308 + t333) * MDP(16) + (t239 * t341 + t365 * t390 + t306) * MDP(18) + (t212 + t254 + (qJD(4) * t340 - t353) * t294 + (t307 - t395) * t298) * MDP(17) + ((t353 + (-t214 - t340) * qJD(4)) * t298 + t307 * t294 + t370) * MDP(19) + (t284 * t338 + t309 - t389) * MDP(8) + (-t284 * t368 + t303 - t390) * MDP(9) + (g(1) * t297 - g(2) * t301) * MDP(2) + (g(1) * t301 + g(2) * t297) * MDP(3); (pkin(1) * t296 * t336 + t331) * MDP(5) + ((t300 * t336 - t355) * pkin(1) + t367) * MDP(6) + (t240 * t285 + (t284 * t299 - t344) * pkin(2) + t309) * MDP(8) + (t241 * t285 + t257 + (-pkin(2) * t284 - t237) * t295 + ((-pkin(2) * t285 - t247) * qJD(3) - t310) * t299 + t369) * MDP(9) + (t222 + (qJD(4) * t316 - t352) * t294 + (-pkin(2) * t344 - t322 + t342) * t298 + t345) * MDP(15) + ((-t352 + (t241 + t316) * qJD(4)) * t298 + ((-t240 + t349) * t285 + t322) * t294 + t333) * MDP(16) + (t212 + (qJD(4) * t317 - t352) * t294 + (-t233 * t285 + t315 - t395) * t298 + t345) * MDP(17) + (t277 * t341 + t306 + (-t241 + t348) * t285 * t365) * MDP(18) + ((t352 + (-t214 - t241 - t317) * qJD(4)) * t298 + (-t285 * t372 + t315) * t294 + t370) * MDP(19) + (t202 * t244 - g(2) * t324 - t323 * t241 + t372 * t214 + (t323 * t361 + t272) * pkin(2) + t304 * t277 + t305) * MDP(20) + t318; (t309 + t385) * MDP(8) + (t303 + t386) * MDP(9) + (t222 + (-pkin(3) * t359 - t392) * t294 + (-t330 + t342) * t298 + t346) * MDP(15) + ((-t392 + (t231 - t398) * qJD(4)) * t298 + (t330 - t385) * t294 + t333) * MDP(16) + (t212 + (t248 * t359 - t392) * t294 + (-t242 * t285 + t319 - t395) * t298 + t346) * MDP(17) + (pkin(8) * t341 - t365 * t386 + t306) * MDP(18) + ((t392 + (-t214 - t231 - t382) * qJD(4)) * t298 + (-t285 * t371 + t319) * t294 + t370) * MDP(19) + (pkin(8) * t304 - g(2) * t332 + t202 * t248 + t214 * t371 - t231 * t323 + t305) * MDP(20) + t326; -MDP(10) * t347 + t366 * t283 * MDP(11) + MDP(12) * t378 + t298 * t284 * MDP(13) + qJDD(4) * MDP(14) + (-t225 * t377 + t320) * MDP(15) + (g(3) * t294 - t208 + (-t225 * t285 + t369) * t298) * MDP(16) + (0.2e1 * t391 + (-t214 * t294 + t236 * t298) * t285 + t314) * MDP(17) + (-t327 * t284 + ((t217 - t360) * t294 + (qJD(5) - t216 - t393) * t298) * t285) * MDP(18) + (0.2e1 * t354 + 0.2e1 * qJD(4) * qJD(5) + t208 + (t236 * t285 - g(3)) * t294 + (t214 * t285 - t369) * t298) * MDP(19) + (-t206 * pkin(4) - g(3) * t328 + t205 * qJ(5) - t214 * t236 - t216 * t387 + t217 * t337 + t369 * t327) * MDP(20); (-qJDD(4) - t347) * MDP(17) + MDP(18) * t378 + (-t283 * t291 - t302) * MDP(19) + (-qJD(4) * t217 + t214 * t377 - t314 + t402) * MDP(20);];
tau = t1;
