% Calculate vector of inverse dynamics joint torques for
% S5PRRRP6
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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:26
% EndTime: 2019-12-05 16:52:31
% DurationCPUTime: 2.72s
% Computational Cost: add. (1814->306), mult. (3776->389), div. (0->0), fcn. (2617->10), ass. (0->149)
t303 = sin(qJ(3));
t393 = pkin(7) + pkin(6);
t265 = t393 * t303;
t305 = cos(qJ(3));
t266 = t393 * t305;
t302 = sin(qJ(4));
t392 = cos(qJ(4));
t235 = -t302 * t265 + t392 * t266;
t299 = qJ(3) + qJ(4);
t293 = sin(t299);
t295 = qJDD(3) + qJDD(4);
t296 = qJD(3) + qJD(4);
t304 = sin(qJ(2));
t300 = sin(pkin(8));
t301 = cos(pkin(8));
t391 = g(1) * t301;
t339 = g(2) * t300 + t391;
t331 = t339 * t304;
t355 = qJD(3) * t393;
t261 = t303 * t355;
t262 = t305 * t355;
t306 = cos(qJ(2));
t353 = t392 * t305;
t375 = t302 * t303;
t326 = t353 - t375;
t322 = t306 * t326;
t327 = -t392 * t265 - t302 * t266;
t371 = -qJD(1) * t322 + t327 * qJD(4) - t392 * t261 - t302 * t262;
t388 = g(3) * t306;
t403 = t235 * t295 + (t331 - t388) * t293 + t371 * t296;
t348 = qJDD(2) * t392;
t359 = qJDD(2) * t303;
t260 = t302 * t305 + t392 * t303;
t402 = t296 * t260;
t215 = qJD(2) * t402 + t302 * t359 - t305 * t348;
t330 = t339 * t306;
t389 = g(3) * t304;
t401 = t330 + t389;
t251 = t326 * t304;
t290 = t295 * qJ(5);
t291 = t296 * qJD(5);
t399 = t290 + t291;
t398 = MDP(17) + MDP(19);
t357 = MDP(18) - MDP(21);
t292 = t295 * pkin(4);
t397 = qJDD(5) - t292;
t308 = qJD(3) ^ 2;
t361 = qJD(1) * qJD(2);
t345 = -qJDD(1) * t306 + t304 * t361;
t395 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t308 + (t339 + t361) * t304 - t345 - t388;
t254 = t260 * qJD(2);
t394 = t254 ^ 2;
t387 = qJD(2) * pkin(2);
t367 = qJD(1) * t304;
t272 = qJD(2) * pkin(6) + t367;
t347 = pkin(7) * qJD(2) + t272;
t248 = t347 * t303;
t240 = qJD(3) * pkin(3) - t248;
t249 = t347 * t305;
t376 = t302 * t249;
t221 = t392 * t240 - t376;
t385 = t221 * t296;
t354 = t392 * t249;
t222 = t302 * t240 + t354;
t384 = t222 * t296;
t341 = qJD(2) * t353;
t365 = qJD(2) * t303;
t352 = t302 * t365;
t252 = -t341 + t352;
t382 = t252 * t254;
t381 = t293 * t304;
t294 = cos(t299);
t380 = t294 * t304;
t379 = t300 * t294;
t378 = t300 * t306;
t377 = t301 * t306;
t374 = t303 * t306;
t309 = qJD(2) ^ 2;
t373 = t305 * t309;
t372 = qJDD(1) - g(3);
t323 = t306 * t260;
t370 = -qJD(1) * t323 + t235 * qJD(4) - t302 * t261 + t392 * t262;
t225 = -t392 * t248 - t376;
t351 = qJD(4) * t392;
t369 = pkin(3) * t351 + qJD(5) - t225;
t297 = t303 ^ 2;
t368 = -t305 ^ 2 + t297;
t366 = qJD(1) * t306;
t364 = qJD(3) * t303;
t363 = qJD(4) * t302;
t362 = qJD(5) - t221;
t360 = qJD(2) * qJD(3);
t358 = qJDD(2) * t305;
t356 = pkin(3) * t364;
t289 = pkin(3) * t305 + pkin(2);
t350 = t303 * t360;
t349 = t305 * t360;
t337 = t296 * t375;
t229 = -qJD(3) * t353 - t305 * t351 + t337;
t202 = pkin(4) * t402 + qJ(5) * t229 - qJD(5) * t260 + t356;
t346 = -t202 + t367;
t257 = qJDD(2) * pkin(6) + qJDD(1) * t304 + t306 * t361;
t220 = -qJD(3) * t272 * t305 + qJDD(3) * pkin(3) - t303 * t257 + (-t349 - t359) * pkin(7);
t223 = -t272 * t364 + t305 * t257 + (-t350 + t358) * pkin(7);
t344 = t302 * t220 + t392 * t223 + t240 * t351 - t249 * t363;
t343 = -t392 * t220 + t302 * t223 + t240 * t363 + t249 * t351;
t342 = -t296 * t341 - t302 * t358 - t303 * t348;
t224 = -t302 * t248 + t354;
t340 = -pkin(3) * t363 + t224;
t338 = g(1) * t300 - g(2) * t301;
t227 = pkin(4) * t254 + qJ(5) * t252;
t334 = qJDD(3) * t305 - t303 * t308;
t333 = qJDD(3) * t303 + t305 * t308;
t332 = pkin(4) * t294 + qJ(5) * t293 + t289;
t329 = t356 - t367;
t205 = -t342 + (t252 - t352) * t296;
t324 = MDP(12) * t382 + t205 * MDP(14) + (t254 * t296 - t215) * MDP(15) + (-t252 ^ 2 + t394) * MDP(13) + t295 * MDP(16);
t258 = -t289 * qJD(2) - t366;
t245 = -t301 * t293 + t294 * t378;
t247 = t293 * t300 + t294 * t377;
t320 = g(1) * t247 + g(2) * t245 + g(3) * t380 - t344;
t244 = t293 * t378 + t294 * t301;
t246 = t293 * t377 - t379;
t319 = g(1) * t246 + g(2) * t244 + g(3) * t381 - t343;
t231 = pkin(3) * t350 - t289 * qJDD(2) + t345;
t273 = -t366 - t387;
t318 = -pkin(6) * qJDD(3) + (t273 + t366 - t387) * qJD(3);
t317 = g(2) * t304 * t379 - t294 * t388 + t295 * t327 - t370 * t296 + t380 * t391;
t217 = pkin(4) * t252 - qJ(5) * t254 + t258;
t316 = -t217 * t252 - t320;
t315 = t258 * t252 + t320;
t314 = -t258 * t254 + t319;
t313 = -t273 * qJD(2) - t257 + t401;
t312 = -t217 * t254 + t319 - t397;
t311 = -g(1) * (-t246 * pkin(4) + qJ(5) * t247) - g(2) * (-t244 * pkin(4) + qJ(5) * t245) - g(3) * (-pkin(4) * t381 + qJ(5) * t380);
t288 = -t392 * pkin(3) - pkin(4);
t284 = pkin(3) * t302 + qJ(5);
t250 = t260 * t304;
t228 = -pkin(4) * t326 - qJ(5) * t260 - t289;
t226 = pkin(3) * t365 + t227;
t214 = qJD(2) * t337 + t342;
t211 = t296 * qJ(5) + t222;
t210 = -t296 * pkin(4) + t362;
t207 = qJD(2) * t323 + t251 * t296;
t206 = qJD(2) * t322 - t304 * t402;
t201 = pkin(4) * t215 + qJ(5) * t214 - qJD(5) * t254 + t231;
t200 = t343 + t397;
t199 = t344 + t399;
t1 = [t372 * MDP(1) + (-t206 * t252 + t207 * t254 - t214 * t250 - t215 * t251) * MDP(20) + (t199 * t251 + t200 * t250 + t206 * t211 + t207 * t210 - g(3)) * MDP(22) + (qJDD(2) * MDP(3) - t309 * MDP(4) + (-0.2e1 * t350 + t358) * MDP(10) + (-0.2e1 * t349 - t359) * MDP(11) - t201 * MDP(22) - t398 * t215 + t357 * t214) * t306 + (-t309 * MDP(3) - qJDD(2) * MDP(4) + (-t333 - t373) * MDP(10) + (t303 * t309 - t334) * MDP(11) + (t217 * MDP(22) + t357 * t254) * qJD(2)) * t304 - t357 * (t206 * t296 + t251 * t295) + t398 * (qJD(2) * t304 * t252 - t207 * t296 - t250 * t295); qJDD(2) * MDP(2) + (t372 * t306 + t331) * MDP(3) + (-t372 * t304 + t330) * MDP(4) + (qJDD(2) * t297 + 0.2e1 * t303 * t349) * MDP(5) + 0.2e1 * (t303 * t358 - t368 * t360) * MDP(6) + t333 * MDP(7) + t334 * MDP(8) + (t318 * t303 + t305 * t395) * MDP(10) + (-t303 * t395 + t318 * t305) * MDP(11) + (-t214 * t260 - t229 * t254) * MDP(12) + (-t214 * t326 - t215 * t260 + t229 * t252 - t254 * t402) * MDP(13) + (-t229 * t296 + t260 * t295) * MDP(14) + (t295 * t326 - t296 * t402) * MDP(15) + (-t215 * t289 - t231 * t326 + t329 * t252 + t258 * t402 + t317) * MDP(17) + (t214 * t289 - t229 * t258 + t231 * t260 + t329 * t254 - t403) * MDP(18) + (-t201 * t326 + t215 * t228 + t217 * t402 - t346 * t252 + t317) * MDP(19) + (t199 * t326 + t200 * t260 - t210 * t229 - t211 * t402 + t214 * t327 - t215 * t235 - t371 * t252 + t370 * t254 - t401) * MDP(20) + (-t201 * t260 + t214 * t228 + t217 * t229 + t346 * t254 + t403) * MDP(21) + (t199 * t235 - t200 * t327 + t201 * t228 + t217 * t202 + t371 * t211 + t370 * t210 + (-g(3) * t332 - t339 * t393) * t306 + (-g(3) * t393 - t217 * qJD(1) + t339 * t332) * t304) * MDP(22); -t303 * MDP(5) * t373 + t368 * MDP(6) * t309 + MDP(7) * t359 + MDP(8) * t358 + qJDD(3) * MDP(9) + (t313 * t303 - t338 * t305) * MDP(10) + (t338 * t303 + t313 * t305) * MDP(11) + (t224 * t296 + (-t252 * t365 + t392 * t295 - t296 * t363) * pkin(3) + t314) * MDP(17) + (t225 * t296 + (-t254 * t365 - t295 * t302 - t296 * t351) * pkin(3) + t315) * MDP(18) + (-t226 * t252 - t288 * t295 + t340 * t296 + t312) * MDP(19) + (-t214 * t288 - t215 * t284 + (t211 - t340) * t254 + (t210 - t369) * t252) * MDP(20) + (t226 * t254 + t284 * t295 + t369 * t296 + t316 + t399) * MDP(21) + (t199 * t284 + t200 * t288 - t217 * t226 - t210 * t224 + t369 * t211 + (t210 * t363 - g(1) * (t300 * t305 - t301 * t374) - g(2) * (-t300 * t374 - t301 * t305) + t303 * t389) * pkin(3) + t311) * MDP(22) + t324; (t314 + t384) * MDP(17) + (t315 + t385) * MDP(18) + (-t227 * t252 + t292 + t312 + t384) * MDP(19) + (pkin(4) * t214 - qJ(5) * t215 + (t211 - t222) * t254 + (t210 - t362) * t252) * MDP(20) + (t227 * t254 + 0.2e1 * t290 + 0.2e1 * t291 + t316 - t385) * MDP(21) + (-t200 * pkin(4) + t199 * qJ(5) - t210 * t222 + t362 * t211 - t217 * t227 + t311) * MDP(22) + t324; (-t295 + t382) * MDP(19) + t205 * MDP(20) + (-t296 ^ 2 - t394) * MDP(21) + (-t211 * t296 - t312) * MDP(22);];
tau = t1;
