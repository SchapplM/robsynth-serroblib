% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:11:03
% EndTime: 2019-03-08 21:11:09
% DurationCPUTime: 3.08s
% Computational Cost: add. (1250->348), mult. (3052->453), div. (0->0), fcn. (1918->8), ass. (0->154)
t289 = sin(qJ(2));
t283 = sin(pkin(6));
t358 = qJD(1) * t283;
t330 = t289 * t358;
t245 = qJD(2) * pkin(8) + t330;
t288 = sin(qJ(3));
t235 = t288 * t245;
t291 = cos(qJ(3));
t284 = cos(pkin(6));
t357 = qJD(1) * t284;
t258 = t291 * t357;
t365 = t235 - t258;
t388 = qJD(4) + t365;
t219 = -qJD(3) * pkin(3) + t388;
t293 = -pkin(3) - pkin(4);
t331 = qJD(3) * t293;
t355 = qJD(2) * t288;
t217 = -qJ(5) * t355 + t365;
t343 = -qJD(4) - t217;
t208 = t331 - t343;
t287 = sin(qJ(6));
t290 = cos(qJ(6));
t387 = t287 * MDP(25) + t290 * MDP(26);
t380 = pkin(8) - qJ(5);
t281 = t288 ^ 2;
t282 = t291 ^ 2;
t386 = MDP(6) * (t281 - t282);
t249 = t380 * t288;
t262 = qJD(6) + t355;
t279 = -pkin(9) + t293;
t301 = pkin(5) * t291 + t279 * t288;
t299 = t301 * qJD(3);
t271 = t288 * qJD(4);
t351 = qJD(3) * t291;
t361 = qJ(4) * t351 + t271;
t385 = (-qJD(6) * t249 + t299 + t330 + t361) * t262;
t372 = t283 * t289;
t234 = t284 * t288 + t291 * t372;
t292 = cos(qJ(2));
t356 = qJD(2) * t283;
t329 = t292 * t356;
t312 = t288 * t329;
t216 = qJD(3) * t234 + t312;
t371 = t283 * t292;
t384 = qJD(6) * t371 + t216;
t224 = t291 * t245 + t288 * t357;
t354 = qJD(2) * t291;
t218 = -qJ(5) * t354 + t224;
t280 = qJD(3) * qJ(4);
t211 = -t218 - t280;
t220 = t280 + t224;
t383 = (-t365 + t235) * qJD(3);
t382 = MDP(13) - MDP(18);
t324 = qJD(1) * t356;
t310 = t292 * t324;
t352 = qJD(3) * t288;
t323 = qJD(1) * t352;
t206 = t245 * t351 + t284 * t323 + t288 * t310;
t350 = qJD(5) * t288;
t202 = (-qJ(5) * t351 - t350) * qJD(2) + t206;
t257 = t292 * t358;
t285 = qJD(2) * pkin(2);
t246 = -t257 - t285;
t226 = -pkin(3) * t354 - qJ(4) * t355 + t246;
t213 = pkin(4) * t354 + qJD(5) - t226;
t309 = pkin(5) * t288 + pkin(9) * t291;
t203 = qJD(2) * t309 + t213;
t209 = qJD(3) * pkin(5) - t211;
t247 = -t291 * pkin(3) - t288 * qJ(4) - pkin(2);
t237 = t291 * pkin(4) - t247;
t225 = t309 + t237;
t230 = t351 * t380 - t350;
t381 = -(qJD(6) * t225 + t230) * t262 + (qJD(3) * t209 - qJD(6) * t203 + t257 * t262 - t202) * t288;
t278 = qJD(3) * qJD(4);
t364 = qJD(3) * t258 + t291 * t310;
t204 = -t245 * t352 + t278 + t364;
t341 = qJD(2) * qJD(3);
t322 = t288 * t341;
t260 = qJ(5) * t322;
t349 = qJD(5) * t291;
t197 = qJD(2) * t349 - t204 - t260;
t379 = t197 * t287;
t378 = t197 * t290;
t344 = t290 * qJD(3);
t347 = qJD(6) * t287;
t326 = t291 * t347;
t363 = qJD(2) * t326 + t290 * t322;
t221 = -qJD(6) * t344 + t363;
t377 = t221 * t287;
t238 = t287 * t354 - t344;
t376 = t238 * t262;
t375 = t238 * t291;
t353 = qJD(3) * t287;
t239 = t290 * t354 + t353;
t374 = t239 * t262;
t373 = t239 * t291;
t295 = qJD(2) ^ 2;
t370 = t283 * t295;
t369 = t287 * t262;
t294 = qJD(3) ^ 2;
t368 = t288 * t294;
t367 = t290 * t262;
t366 = t291 * t294;
t321 = t291 * t341;
t362 = qJ(4) * t321 + qJD(2) * t271;
t359 = t281 + t282;
t205 = qJD(3) * t279 - t343;
t348 = qJD(6) * t205;
t346 = qJD(6) * t290;
t345 = t209 * qJD(6);
t342 = -qJD(5) - t213;
t340 = MDP(12) - MDP(17);
t338 = pkin(3) * t352;
t337 = t288 * t369;
t336 = t288 * t367;
t335 = t288 * t372;
t334 = t289 * t370;
t333 = t288 * t291 * t295;
t332 = 0.2e1 * t278 + t364;
t327 = t262 * t347;
t325 = t262 * t346;
t320 = MDP(24) * t354;
t212 = (t330 + t338) * qJD(2) - t362;
t319 = -pkin(8) * t294 - t212;
t318 = t246 - t285;
t315 = qJD(2) * t247 + t226;
t313 = t288 * t331;
t195 = t203 * t290 - t205 * t287;
t196 = t203 * t287 + t205 * t290;
t308 = -t208 * t288 + t211 * t291;
t306 = qJD(2) * t282 - t262 * t288;
t305 = t290 * MDP(25) - t287 * MDP(26);
t303 = qJD(3) * t224 - t206;
t298 = t226 * MDP(15) - t262 * t305;
t297 = t204 * t291 + t206 * t288 + (t219 * t291 - t220 * t288) * qJD(3);
t286 = qJ(4) + pkin(5);
t266 = qJ(4) * t354;
t253 = t287 * t322;
t250 = t380 * t291;
t243 = t289 * t288 * t324;
t241 = t323 * t371;
t240 = pkin(3) * t355 - t266;
t233 = -t284 * t291 + t335;
t231 = t338 - t361;
t229 = t380 * t352 + t349;
t228 = t293 * t355 + t266;
t227 = t313 + t361;
t222 = qJD(6) * t239 - t253;
t215 = -qJD(3) * t335 + (qJD(3) * t284 + t329) * t291;
t214 = qJD(2) * t301 + t266;
t207 = (t313 - t330) * qJD(2) + t362;
t201 = (t299 - t330) * qJD(2) + t362;
t198 = t290 * t201;
t1 = [(t204 * t234 + t206 * t233 - t212 * t371 + t215 * t220 + t216 * t219) * MDP(15) + (-t197 * t234 + t202 * t233 + t207 * t371 + t208 * t216 - t211 * t215) * MDP(19) + ((-t233 * t346 - t384 * t287) * t262 - t215 * t238 - t234 * t222) * MDP(25) + (-(-t233 * t347 + t384 * t290) * t262 - t215 * t239 + t234 * t221) * MDP(26) + (-MDP(10) * t216 + MDP(16) * t215) * qJD(3) + (-MDP(4) * t292 + (-MDP(10) * t291 - MDP(16) * t288 - MDP(3)) * t289) * t370 + (-MDP(11) + MDP(14)) * (-t288 * t334 + (t291 * t329 + t215) * qJD(3)) + t340 * (-t291 * t334 + (-t216 - t312) * qJD(3)) + ((-MDP(19) * t213 + t298) * t372 + ((-MDP(10) * t371 - t234 * t382) * t288 + ((MDP(16) + t305) * t371 + (-t387 + t382) * t233) * t291) * qJD(3) + t382 * (t215 * t291 + t216 * t288)) * qJD(2); 0.2e1 * t288 * MDP(5) * t321 - 0.2e1 * t341 * t386 + MDP(7) * t366 - MDP(8) * t368 + (-pkin(8) * t366 + t318 * t352 + t241) * MDP(10) + (pkin(8) * t368 + (t318 + t257) * t351) * MDP(11) + (t241 + t315 * t352 + ((-t231 + t330) * qJD(2) + t319) * t291) * MDP(12) + (-t310 * t359 + t297) * MDP(13) + (t243 + (-qJD(2) * t231 + t319) * t288 + (-t315 - t257) * t351) * MDP(14) + (t212 * t247 + t226 * t231 + (-t226 * t289 + (-t219 * t288 - t220 * t291) * t292) * t358 + t297 * pkin(8)) * MDP(15) + (t243 + (qJD(2) * t227 + t207) * t288 + (-t229 + (qJD(2) * t237 + t213 - t257) * t291) * qJD(3)) * MDP(16) + (-t207 * t291 - t241 + (t213 * t288 + t230) * qJD(3) + (t237 * t352 + (-t227 - t330) * t291) * qJD(2)) * MDP(17) + (t197 * t291 - t202 * t288 + (-t208 * t291 - t211 * t288) * qJD(3) + (t229 * t291 - t230 * t288 + (-t249 * t291 + t250 * t288) * qJD(3) + t359 * t257) * qJD(2)) * MDP(18) + (-t197 * t250 + t202 * t249 + t207 * t237 + t208 * t230 + t211 * t229 + t213 * t227 + (t213 * t289 + t292 * t308) * t358) * MDP(19) + (-t221 * t290 * t291 + (-t288 * t344 - t326) * t239) * MDP(20) + ((t238 * t290 + t239 * t287) * t352 + (t377 - t222 * t290 + (t238 * t287 - t239 * t290) * qJD(6)) * t291) * MDP(21) + (t262 * t326 + t221 * t288 + (-t290 * t306 - t373) * qJD(3)) * MDP(22) + (t291 * t325 + t222 * t288 + (t287 * t306 + t375) * qJD(3)) * MDP(23) + (t262 + t355) * MDP(24) * t351 + (t198 * t288 - t250 * t222 + t229 * t238 + (-t288 * t348 + t385) * t290 + t381 * t287 + (t238 * t257 - t290 * t345 + t379 + ((t225 * t290 - t249 * t287) * qJD(2) + t195) * qJD(3)) * t291) * MDP(25) + (t250 * t221 + t229 * t239 + (-(t201 - t348) * t288 - t385) * t287 + t381 * t290 + (t239 * t257 + t287 * t345 + t378 + (-(t225 * t287 + t249 * t290) * qJD(2) - t196) * qJD(3)) * t291) * MDP(26); -MDP(5) * t333 + t295 * t386 + (-t246 * t355 + t303) * MDP(10) + (-t246 * t354 - t364 + t383) * MDP(11) + t303 * MDP(12) + (t332 - t383) * MDP(14) + (-pkin(3) * t206 + qJ(4) * t204 - t219 * t224 + t220 * t388 - t226 * t240) * MDP(15) + (t260 + (t217 - t235) * qJD(3) + t332) * MDP(16) + (-qJD(3) * t218 + t206) * MDP(17) + (-qJ(4) * t197 + t202 * t293 - t208 * t218 + t211 * t343 - t213 * t228) * MDP(19) + (t239 * t367 - t377) * MDP(20) + ((-t221 - t376) * t290 + (-t222 - t374) * t287) * MDP(21) - t325 * MDP(22) + t327 * MDP(23) - t262 * t320 + (-t286 * t222 - t378 - (t214 * t290 - t218 * t287) * t262 + t343 * t238 + (-t209 * t287 - t279 * t367) * qJD(6)) * MDP(25) + (t286 * t221 + t379 + (t214 * t287 + t218 * t290) * t262 + t343 * t239 + (-t209 * t290 + t279 * t369) * qJD(6)) * MDP(26) + (-t336 * MDP(22) + t337 * MDP(23) + t387 * (-t209 * t288 - t279 * t351) + (-t226 * MDP(12) + t240 * MDP(14) - t228 * MDP(16) + t342 * MDP(17)) * t288 + (t240 * MDP(12) + t226 * MDP(14) + t342 * MDP(16) + (-qJ(5) * qJD(3) + t228) * MDP(17) + (t239 - t353) * MDP(22) + (-t238 - t344) * MDP(23) - t195 * MDP(25) + t196 * MDP(26)) * t291) * qJD(2); (-qJD(3) * t220 + t206) * MDP(15) + (qJD(3) * t211 + t206) * MDP(19) + (qJD(3) * t238 - t325) * MDP(25) + (qJD(3) * t239 + t327) * MDP(26) - t340 * t333 + (MDP(14) + MDP(16)) * (-t281 * t295 - t294) + ((-MDP(19) * qJ(5) - t387) * t351 + (MDP(19) * t342 + t298) * t288) * qJD(2); t362 * MDP(19) - t359 * MDP(18) * t295 - t387 * t262 * qJD(6) + ((-t308 - t330) * MDP(19) + (-t337 - t375) * MDP(25) + (-t336 - t373) * MDP(26) + ((0.2e1 * MDP(16) + t305) * t291 + (MDP(19) * t293 + 0.2e1 * MDP(17)) * t288) * qJD(3)) * qJD(2); t239 * t238 * MDP(20) + (-t238 ^ 2 + t239 ^ 2) * MDP(21) + (t363 - t376) * MDP(22) + (-t253 - t374) * MDP(23) + qJD(3) * t320 + (t196 * t262 - t202 * t287 + t209 * t239 + t198) * MDP(25) + (t195 * t262 - t201 * t287 - t202 * t290 - t209 * t238) * MDP(26) + (-MDP(22) * t344 + MDP(23) * t239 - MDP(25) * t196 - MDP(26) * t195) * qJD(6);];
tauc  = t1;
