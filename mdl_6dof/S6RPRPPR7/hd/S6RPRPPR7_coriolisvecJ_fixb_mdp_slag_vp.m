% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:34
% EndTime: 2019-03-09 02:57:39
% DurationCPUTime: 2.66s
% Computational Cost: add. (2001->310), mult. (4288->403), div. (0->0), fcn. (2803->6), ass. (0->141)
t381 = MDP(14) + MDP(16);
t294 = sin(qJ(3));
t367 = sin(pkin(9));
t321 = t367 * t294;
t273 = qJD(1) * t321;
t296 = cos(qJ(3));
t368 = cos(pkin(9));
t322 = t368 * t296;
t316 = qJD(1) * t322;
t262 = t316 - t273;
t256 = qJD(6) + t262;
t380 = t256 ^ 2;
t379 = MDP(8) * (t294 ^ 2 - t296 ^ 2);
t378 = qJ(2) * MDP(6) + MDP(5);
t297 = -pkin(1) - pkin(7);
t272 = qJD(1) * t297 + qJD(2);
t340 = qJD(1) * t294;
t254 = -qJ(4) * t340 + t272 * t294;
t243 = t367 * t254;
t339 = qJD(1) * t296;
t255 = -qJ(4) * t339 + t296 * t272;
t225 = t255 * t368 - t243;
t329 = -qJD(5) + t225;
t333 = qJD(4) * t296;
t336 = qJD(3) * t294;
t233 = -t272 * t336 + (qJ(4) * t336 - t333) * qJD(1);
t334 = qJD(4) * t294;
t335 = qJD(3) * t296;
t234 = t272 * t335 + (-qJ(4) * t335 - t334) * qJD(1);
t207 = t367 * t233 + t368 * t234;
t246 = qJD(3) * pkin(3) + t255;
t217 = t246 * t368 - t243;
t323 = t368 * t254;
t218 = t367 * t246 + t323;
t319 = qJD(3) * t367;
t320 = qJD(3) * t368;
t260 = t294 * t319 - t296 * t320;
t261 = -t294 * t320 - t296 * t319;
t302 = t294 * t368 + t296 * t367;
t376 = -t207 * t302 - t217 * t261 + t218 * t260;
t202 = -qJD(3) * qJD(5) - t207;
t311 = qJD(5) - t217;
t212 = -qJD(3) * pkin(4) + t311;
t214 = -qJD(3) * qJ(5) - t218;
t375 = t202 * t302 + t212 * t261 - t214 * t260;
t258 = t302 * qJD(1);
t370 = t258 * pkin(5);
t201 = -t214 - t370;
t224 = t255 * t367 + t323;
t327 = qJD(1) * qJD(3);
t250 = t302 * t327;
t283 = -pkin(3) * t368 - pkin(4);
t277 = -pkin(8) + t283;
t374 = -t277 * t250 + (t201 - t224 + t370) * t256;
t257 = t262 ^ 2;
t372 = pkin(4) + pkin(8);
t271 = qJD(3) * t316;
t249 = qJD(3) * t273 - t271;
t371 = pkin(4) * t249;
t369 = t262 * pkin(5);
t206 = -t368 * t233 + t234 * t367;
t345 = qJ(4) - t297;
t269 = t345 * t294;
t270 = t345 * t296;
t231 = -t269 * t367 + t368 * t270;
t364 = t206 * t231;
t265 = -t321 + t322;
t363 = t206 * t265;
t346 = t294 * pkin(3) + qJ(2);
t315 = -qJ(5) * t265 + t346;
t213 = t302 * t372 + t315;
t360 = t213 * t250;
t293 = sin(qJ(6));
t332 = qJD(6) * t293;
t295 = cos(qJ(6));
t331 = qJD(6) * t295;
t344 = -t293 * t249 + t258 * t331;
t215 = -qJD(3) * t332 + t344;
t358 = t215 * t295;
t337 = qJD(3) * t293;
t235 = -t295 * t258 + t337;
t355 = t235 * t256;
t354 = t235 * t258;
t237 = qJD(3) * t295 + t258 * t293;
t353 = t237 * t256;
t352 = t237 * t258;
t351 = t250 * t265;
t350 = t250 * t293;
t349 = t302 * t293;
t242 = t295 * t250;
t299 = qJD(1) ^ 2;
t348 = t296 * t299;
t298 = qJD(3) ^ 2;
t347 = t297 * t298;
t289 = qJD(1) * qJD(2);
t324 = t296 * t327;
t343 = pkin(3) * t324 + t289;
t330 = pkin(3) * t335 + qJD(2);
t328 = t369 - t329;
t326 = 0.2e1 * qJD(1);
t268 = pkin(3) * t340 + qJD(1) * qJ(2) + qJD(4);
t318 = pkin(3) * t339 + qJ(5) * t258;
t317 = t256 * t293;
t314 = qJ(5) * t250 + t343;
t251 = t336 * t345 - t333;
t252 = -qJD(3) * t270 - t334;
t220 = -t368 * t251 + t252 * t367;
t199 = -qJD(3) * t372 + t311 + t369;
t312 = -qJ(5) * t262 + t268;
t203 = t258 * t372 + t312;
t192 = t199 * t295 - t203 * t293;
t193 = t199 * t293 + t203 * t295;
t310 = -t256 * t317 - t242;
t309 = -t260 * t293 + t302 * t331;
t195 = pkin(5) * t249 - t202;
t308 = t195 + (-qJD(6) * t277 + t262 * t372 + t318) * t256;
t307 = -qJD(5) * t262 + t314;
t306 = -qJ(5) * t261 - qJD(5) * t265 + t330;
t219 = pkin(4) * t258 + t312;
t305 = t219 * t262 + t206;
t222 = t265 * pkin(5) + t231;
t304 = t195 * t302 - t201 * t260 + t222 * t250;
t303 = -t380 * t295 + t350;
t221 = t251 * t367 + t252 * t368;
t232 = -t269 * t368 - t270 * t367;
t300 = t220 * t262 - t221 * t258 - t231 * t250 + t232 * t249 + t363;
t279 = pkin(3) * t367 + qJ(5);
t241 = t295 * t249;
t227 = pkin(4) * t302 + t315;
t226 = pkin(4) * t262 + t318;
t223 = -pkin(5) * t302 + t232;
t216 = t237 * qJD(6) + t241;
t211 = -pkin(4) * t260 + t306;
t205 = t260 * pkin(5) + t221;
t204 = t261 * pkin(5) + t220;
t200 = t307 - t371;
t198 = -t260 * t372 + t306;
t197 = -t250 * pkin(5) + t206;
t196 = t295 * t197;
t194 = -t249 * t372 + t307;
t1 = [0.2e1 * t327 * t379 - t298 * t296 * MDP(10) + qJ(2) * t335 * t326 * MDP(12) + (-t296 * t347 + (-qJ(2) * t336 + qJD(2) * t296) * t326) * MDP(13) + (t300 + t376) * MDP(14) + (t207 * t232 - t217 * t220 + t218 * t221 + t268 * t330 + t343 * t346 + t364) * MDP(15) + (t300 + t375) * MDP(16) + (qJD(3) * t220 - t200 * t302 - t211 * t258 + t219 * t260 + t227 * t249) * MDP(17) + (qJD(3) * t221 - t200 * t265 - t211 * t262 - t219 * t261 + t227 * t250) * MDP(18) + (t200 * t227 - t202 * t232 + t211 * t219 + t212 * t220 - t214 * t221 + t364) * MDP(19) + (t215 * t349 + t237 * t309) * MDP(20) + ((t235 * t293 - t237 * t295) * t260 - (-t358 + t216 * t293 + (t235 * t295 + t237 * t293) * qJD(6)) * t302) * MDP(21) + (t215 * t265 + t237 * t261 - t250 * t349 + t256 * t309) * MDP(22) + (-t302 * t242 - t216 * t265 - t235 * t261 + (-t260 * t295 - t302 * t332) * t256) * MDP(23) + (t256 * t261 - t351) * MDP(24) + (t192 * t261 + t196 * t265 + t205 * t235 + t223 * t216 + (-t194 * t265 - t198 * t256 + t360) * t293 + (t204 * t256 - t304) * t295 + ((-t213 * t295 - t222 * t293) * t256 - t193 * t265 + t201 * t349) * qJD(6)) * MDP(25) + (-t193 * t261 + t205 * t237 + t223 * t215 + (-(qJD(6) * t222 + t198) * t256 + t360 - (qJD(6) * t199 + t194) * t265 + t201 * qJD(6) * t302) * t295 + (-(-qJD(6) * t213 + t204) * t256 - (-qJD(6) * t203 + t197) * t265 + t304) * t293) * MDP(26) + 0.2e1 * t378 * t289 + (-0.2e1 * MDP(7) * t324 - t298 * MDP(9) + (qJD(2) * t326 - t347) * MDP(12)) * t294; (-qJD(1) * t268 - t363 - t376) * MDP(15) + (qJD(1) * t258 - qJD(3) * t261) * MDP(17) + (qJD(1) * t262 - qJD(3) * t260) * MDP(18) + (-qJD(1) * t219 - t363 - t375) * MDP(19) + (t216 * t302 - t235 * t260 + t265 * t242) * MDP(25) + (t215 * t302 - t237 * t260 - t265 * t350) * MDP(26) - t378 * t299 + ((qJD(1) * t293 - t261 * t295 + t265 * t332) * MDP(25) + (qJD(1) * t295 + t261 * t293 + t265 * t331) * MDP(26)) * t256 + (MDP(12) * t294 + MDP(13) * t296) * (-t298 - t299) + t381 * (t249 * t302 + t258 * t260 - t261 * t262 + t351); t294 * MDP(7) * t348 - t299 * t379 + ((t218 - t224) * t262 + (-t217 + t225) * t258 + (t249 * t367 + t250 * t368) * pkin(3)) * MDP(14) + (t217 * t224 - t218 * t225 + (-t206 * t368 + t207 * t367 - t268 * t339) * pkin(3)) * MDP(15) + (t249 * t279 - t250 * t283 + (-t214 - t224) * t262 + (t212 + t329) * t258) * MDP(16) + (-t224 * qJD(3) + t226 * t258 + t305) * MDP(17) + (-t219 * t258 + t226 * t262 + (0.2e1 * qJD(5) - t225) * qJD(3) + t207) * MDP(18) + (-t202 * t279 + t206 * t283 - t212 * t224 + t214 * t329 - t219 * t226) * MDP(19) + (-t237 * t317 + t358) * MDP(20) + ((-t216 - t353) * t295 + (-t215 + t355) * t293) * MDP(21) + (t310 + t352) * MDP(22) + (t303 - t354) * MDP(23) + t256 * t258 * MDP(24) + (t192 * t258 + t279 * t216 + t328 * t235 + t308 * t293 + t374 * t295) * MDP(25) + (-t193 * t258 + t279 * t215 + t328 * t237 - t374 * t293 + t308 * t295) * MDP(26) + (MDP(13) * t294 * t299 - MDP(12) * t348) * qJ(2); (t217 * t262 + t218 * t258 + t343) * MDP(15) + (-t271 + (t273 - t262) * qJD(3)) * MDP(17) + (qJD(3) * t258 + t250) * MDP(18) + (-t371 - t214 * t258 + (-qJD(5) - t212) * t262 + t314) * MDP(19) + (t303 + t354) * MDP(25) + (t380 * t293 + t242 + t352) * MDP(26) + t381 * (-t258 ^ 2 - t257); -t262 * t258 * MDP(17) + (-t257 - t298) * MDP(18) + (t214 * qJD(3) + t305) * MDP(19) + (-qJD(3) * t235 + t310) * MDP(25) + (-qJD(3) * t237 + t303) * MDP(26); t237 * t235 * MDP(20) + (-t235 ^ 2 + t237 ^ 2) * MDP(21) + (t344 + t355) * MDP(22) + (-t241 + t353) * MDP(23) - t250 * MDP(24) + (t193 * t256 - t194 * t293 - t201 * t237 + t196) * MDP(25) + (t192 * t256 - t194 * t295 - t197 * t293 + t201 * t235) * MDP(26) + (-MDP(22) * t337 - MDP(23) * t237 - MDP(25) * t193 - MDP(26) * t192) * qJD(6);];
tauc  = t1;
