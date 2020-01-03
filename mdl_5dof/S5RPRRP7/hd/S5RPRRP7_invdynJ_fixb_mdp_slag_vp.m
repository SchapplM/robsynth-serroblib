% Calculate vector of inverse dynamics joint torques for
% S5RPRRP7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:44
% EndTime: 2019-12-31 18:45:49
% DurationCPUTime: 3.54s
% Computational Cost: add. (2133->366), mult. (4253->466), div. (0->0), fcn. (2602->10), ass. (0->149)
t303 = sin(pkin(8));
t289 = pkin(1) * t303 + pkin(6);
t278 = t289 * qJDD(1);
t410 = -qJD(2) * qJD(3) - t278;
t300 = qJ(1) + pkin(8);
t293 = sin(t300);
t294 = cos(t300);
t337 = g(1) * t294 + g(2) * t293;
t280 = t289 * qJD(1);
t306 = sin(qJ(3));
t309 = cos(qJ(3));
t251 = qJD(2) * t309 - t306 * t280;
t409 = qJD(3) * t251;
t367 = qJD(1) * t309;
t408 = qJD(4) - t367;
t359 = qJD(4) * t306;
t407 = qJD(1) * t359 - qJDD(3);
t305 = sin(qJ(4));
t308 = cos(qJ(4));
t352 = qJDD(1) * t306;
t230 = (qJD(3) * (qJD(4) + t367) + t352) * t305 + t407 * t308;
t406 = t306 * t337;
t329 = pkin(3) * t309 + pkin(7) * t306 + pkin(2);
t304 = cos(pkin(8));
t401 = pkin(1) * t304;
t259 = -t329 - t401;
t379 = t308 * t309;
t372 = t305 * t259 + t289 * t379;
t241 = -qJD(3) * pkin(3) - t251;
t356 = t308 * qJD(3);
t368 = qJD(1) * t306;
t266 = t305 * t368 - t356;
t364 = qJD(3) * t305;
t268 = t308 * t368 + t364;
t217 = pkin(4) * t266 - qJ(5) * t268 + t241;
t297 = t309 * qJDD(1);
t354 = qJD(1) * qJD(3);
t262 = t306 * t354 + qJDD(4) - t297;
t399 = pkin(7) * t262;
t405 = -t217 * t408 + t399;
t351 = MDP(17) + MDP(19);
t404 = t308 * MDP(18) + t305 * t351;
t402 = t268 ^ 2;
t400 = pkin(4) * t262;
t396 = g(3) * t306;
t395 = g(3) * t309;
t394 = pkin(7) * qJD(4);
t393 = qJ(5) * t262;
t252 = t306 * qJD(2) + t309 * t280;
t242 = qJD(3) * pkin(7) + t252;
t243 = t259 * qJD(1);
t219 = t242 * t308 + t243 * t305;
t214 = qJ(5) * t408 + t219;
t392 = t214 * t408;
t391 = t219 * t408;
t390 = t266 * t268;
t389 = t266 * t408;
t388 = t268 * t408;
t387 = t268 * t306;
t386 = t268 * t308;
t340 = pkin(3) * t306 - pkin(7) * t309;
t272 = t340 * qJD(3);
t385 = t272 * t308;
t384 = t289 * t305;
t383 = t289 * t308;
t382 = t305 * t309;
t381 = t306 * t308;
t380 = t308 * t262;
t378 = qJDD(2) - g(3);
t346 = t309 * t356;
t377 = -t230 * t381 - t266 * t346;
t271 = t340 * qJD(1);
t376 = t308 * t251 + t305 * t271;
t375 = t306 * t380 + t346 * t408;
t358 = qJD(4) * t308;
t374 = t259 * t358 + t305 * t272;
t333 = pkin(4) * t305 - qJ(5) * t308;
t373 = -qJD(5) * t305 + t408 * t333 - t252;
t371 = t337 * t381;
t301 = t306 ^ 2;
t370 = -t309 ^ 2 + t301;
t369 = MDP(20) * t305;
t290 = -pkin(2) - t401;
t281 = qJD(1) * t290;
t365 = qJD(3) * t266;
t363 = qJD(3) * t306;
t362 = qJD(3) * t309;
t361 = qJD(4) * t266;
t360 = qJD(4) * t305;
t218 = -t242 * t305 + t243 * t308;
t355 = qJD(5) - t218;
t350 = MDP(18) - MDP(21);
t225 = qJDD(3) * pkin(7) + qJDD(2) * t306 + t278 * t309 + t409;
t231 = qJD(1) * t272 + qJDD(1) * t259;
t349 = t308 * t225 + t305 * t231 + t243 * t358;
t348 = -t280 * t362 + t410 * t306;
t347 = t305 * t359;
t345 = t309 * t354;
t343 = pkin(4) + t384;
t341 = t305 * t225 - t308 * t231 + t242 * t358 + t243 * t360;
t245 = t293 * t382 + t294 * t308;
t247 = -t293 * t308 + t294 * t382;
t339 = -g(1) * t245 + g(2) * t247;
t246 = t293 * t379 - t294 * t305;
t248 = t293 * t305 + t294 * t379;
t338 = g(1) * t246 - g(2) * t248;
t336 = g(1) * t293 - g(2) * t294;
t307 = sin(qJ(1));
t310 = cos(qJ(1));
t335 = g(1) * t307 - g(2) * t310;
t334 = pkin(4) * t308 + qJ(5) * t305;
t322 = -t242 * t360 + t349;
t208 = qJD(5) * t408 + t322 + t393;
t209 = qJDD(5) + t341 - t400;
t332 = t208 * t308 + t209 * t305;
t213 = -pkin(4) * t408 + t355;
t331 = t213 * t308 - t214 * t305;
t330 = t213 * t305 + t214 * t308;
t328 = pkin(3) + t334;
t327 = t394 * t408 + t395;
t325 = t289 + t333;
t324 = -t262 * t305 - t358 * t408;
t226 = -qJDD(3) * pkin(3) - qJDD(2) * t309 - t348;
t229 = -qJD(4) * t356 + (-t345 - t352) * t308 + t407 * t305;
t210 = pkin(4) * t230 + qJ(5) * t229 - qJD(5) * t268 + t226;
t323 = -t210 - t327;
t321 = -qJD(1) * t281 + t337;
t320 = t241 * t408 - t399;
t319 = 0.2e1 * qJD(3) * t281 - qJDD(3) * t289;
t318 = -t309 * t337 - t396;
t317 = g(1) * t247 + g(2) * t245 + t305 * t396 - t341;
t311 = qJD(3) ^ 2;
t316 = -0.2e1 * qJDD(1) * t290 - t289 * t311 + t336;
t315 = t217 * t268 + qJDD(5) - t317;
t314 = -g(1) * t248 - g(2) * t246 - g(3) * t381 + t322;
t277 = qJDD(3) * t309 - t306 * t311;
t276 = qJDD(3) * t306 + t309 * t311;
t254 = t268 * t363;
t239 = t325 * t306;
t235 = pkin(4) * t268 + qJ(5) * t266;
t234 = -t259 * t308 + t309 * t343;
t233 = -qJ(5) * t309 + t372;
t223 = -pkin(4) * t368 + t251 * t305 - t271 * t308;
t222 = qJ(5) * t368 + t376;
t216 = (qJD(4) * t334 - qJD(5) * t308) * t306 + t325 * t362;
t215 = -t229 + t389;
t212 = t372 * qJD(4) - t343 * t363 - t385;
t211 = (-t289 * t360 - qJD(5)) * t309 + (qJ(5) - t383) * t363 + t374;
t1 = [qJDD(1) * MDP(1) + t335 * MDP(2) + (g(1) * t310 + g(2) * t307) * MDP(3) + (t335 + (t303 ^ 2 + t304 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t301 + 0.2e1 * t306 * t345) * MDP(5) + 0.2e1 * (t297 * t306 - t354 * t370) * MDP(6) + t276 * MDP(7) + t277 * MDP(8) + (t306 * t319 + t309 * t316) * MDP(10) + (-t306 * t316 + t309 * t319) * MDP(11) + (-t229 * t381 + (t346 - t347) * t268) * MDP(12) + (-t358 * t387 + (-t268 * t362 + (t229 + t361) * t306) * t305 + t377) * MDP(13) + (t229 * t309 - t347 * t408 + t254 + t375) * MDP(14) + ((-t364 * t408 + t230) * t309 + (t324 - t365) * t306) * MDP(15) + (-t262 * t309 + t363 * t408) * MDP(16) + ((-t259 * t360 + t385) * t408 + t259 * t380 + (t241 * t364 + (t324 + t365) * t289 + t341) * t309 + (t241 * t358 + t226 * t305 + t289 * t230 + (t384 * t408 + t218) * qJD(3)) * t306 + t338) * MDP(17) + (-t374 * t408 - t372 * t262 + ((t289 * t408 - t242) * t360 + (t241 * t308 + t268 * t289) * qJD(3) + t349) * t309 + (-t241 * t360 + t226 * t308 - t289 * t229 + (t383 * t408 - t219) * qJD(3)) * t306 + t339) * MDP(18) + (-t212 * t408 + t216 * t266 + t230 * t239 - t234 * t262 + (t217 * t364 + t209) * t309 + (-qJD(3) * t213 + t210 * t305 + t217 * t358) * t306 + t338) * MDP(19) + (-t211 * t266 + t212 * t268 - t229 * t234 - t230 * t233 + t331 * t362 + (-qJD(4) * t330 - t208 * t305 + t209 * t308 + t336) * t306) * MDP(20) + (t211 * t408 - t216 * t268 + t229 * t239 + t233 * t262 + (-t217 * t356 - t208) * t309 + (qJD(3) * t214 - t210 * t308 + t217 * t360) * t306 - t339) * MDP(21) + (t208 * t233 + t214 * t211 + t210 * t239 + t217 * t216 + t209 * t234 + t213 * t212 - g(1) * (-pkin(1) * t307 - pkin(4) * t246 - qJ(5) * t245) - g(2) * (pkin(1) * t310 + pkin(4) * t248 + qJ(5) * t247) + (-g(1) * pkin(6) - g(2) * t329) * t294 + (-g(2) * pkin(6) + g(1) * t329) * t293) * MDP(22); t378 * MDP(4) + t277 * MDP(10) - t276 * MDP(11) + t254 * MDP(18) + t377 * MDP(20) + t375 * MDP(21) - g(3) * MDP(22) + t351 * t266 * t363 + (-t210 * MDP(22) - t351 * t230 + t350 * t229 + (t330 * MDP(22) + t268 * t369 - t404 * t408) * qJD(3)) * t309 + (-t229 * t369 - qJD(3) * t268 * MDP(21) + (qJD(3) * t217 + t332) * MDP(22) - t404 * t262 + ((t266 * t305 + t386) * MDP(20) + t331 * MDP(22) - (-t305 * t350 + t308 * t351) * t408) * qJD(4)) * t306; MDP(7) * t352 + MDP(8) * t297 + qJDD(3) * MDP(9) + (qJD(3) * t252 + t306 * t321 + t309 * t378 + t348) * MDP(10) + (t409 + (qJD(3) * t280 - t378) * t306 + (t321 + t410) * t309) * MDP(11) + (-t229 * t305 + t386 * t408) * MDP(12) + ((-t229 - t389) * t308 + (-t230 - t388) * t305) * MDP(13) + ((-t379 * t408 - t387) * qJD(1) - t324) * MDP(14) + (-t408 * t360 + t380 + (t266 * t306 + t382 * t408) * qJD(1)) * MDP(15) - t408 * MDP(16) * t368 + (-t218 * t368 - pkin(3) * t230 - t252 * t266 + (-t395 - t226 - (t271 + t394) * t408) * t308 + (t251 * t408 + t320) * t305 + t371) * MDP(17) + (pkin(3) * t229 + t376 * t408 + t219 * t368 - t252 * t268 + t320 * t308 + (t226 + t327 - t406) * t305) * MDP(18) + (t213 * t368 + t223 * t408 - t230 * t328 + t373 * t266 - t405 * t305 + t323 * t308 + t371) * MDP(19) + (t222 * t266 - t223 * t268 + (t208 + t408 * t213 + (qJD(4) * t268 - t230) * pkin(7)) * t308 + (t209 - t392 + (-t229 + t361) * pkin(7)) * t305 + t318) * MDP(20) + (-t214 * t368 - t222 * t408 - t229 * t328 - t373 * t268 + t405 * t308 + (t323 + t406) * t305) * MDP(21) + (-t213 * t223 - t214 * t222 + t373 * t217 + (qJD(4) * t331 + t318 + t332) * pkin(7) + (-t210 - t395 + t406) * t328) * MDP(22) + (-t306 * t309 * MDP(5) + t370 * MDP(6)) * qJD(1) ^ 2; MDP(12) * t390 + (-t266 ^ 2 + t402) * MDP(13) + t215 * MDP(14) + (-t230 + t388) * MDP(15) + t262 * MDP(16) + (-t241 * t268 + t317 + t391) * MDP(17) + (t218 * t408 + t241 * t266 - t314) * MDP(18) + (-t235 * t266 - t315 + t391 + 0.2e1 * t400) * MDP(19) + (pkin(4) * t229 - qJ(5) * t230 + (t214 - t219) * t268 + (t213 - t355) * t266) * MDP(20) + (0.2e1 * t393 - t217 * t266 + t235 * t268 - (-0.2e1 * qJD(5) + t218) * t408 + t314) * MDP(21) + (t208 * qJ(5) - t209 * pkin(4) - t217 * t235 - t213 * t219 - g(1) * (-pkin(4) * t247 + qJ(5) * t248) - g(2) * (-pkin(4) * t245 + qJ(5) * t246) + t333 * t396 + t355 * t214) * MDP(22); (-t262 + t390) * MDP(19) + t215 * MDP(20) + (-t408 ^ 2 - t402) * MDP(21) + (t315 - t392 - t400) * MDP(22);];
tau = t1;
