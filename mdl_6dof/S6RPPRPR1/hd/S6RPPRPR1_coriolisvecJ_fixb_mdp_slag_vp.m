% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:08
% EndTime: 2019-03-09 01:40:15
% DurationCPUTime: 3.95s
% Computational Cost: add. (3007->323), mult. (7406->451), div. (0->0), fcn. (5640->10), ass. (0->135)
t336 = cos(pkin(10));
t393 = cos(qJ(4));
t360 = qJD(1) * t393;
t322 = t336 * t360;
t333 = sin(pkin(10));
t339 = sin(qJ(4));
t378 = t339 * t333;
t361 = qJD(1) * t378;
t301 = -t322 + t361;
t298 = qJD(6) + t301;
t323 = sin(pkin(9)) * pkin(1) + qJ(3);
t318 = t323 * qJD(1);
t329 = t336 * qJD(2);
t389 = pkin(7) * qJD(1);
t289 = t329 + (-t318 - t389) * t333;
t294 = t333 * qJD(2) + t336 * t318;
t290 = t336 * t389 + t294;
t253 = t339 * t289 + t393 * t290;
t399 = qJD(4) * t253;
t359 = qJD(1) * (t333 ^ 2 + t336 ^ 2);
t398 = MDP(7) * t359;
t377 = t339 * t336;
t314 = t393 * t333 + t377;
t303 = t314 * qJD(1);
t332 = sin(pkin(11));
t335 = cos(pkin(11));
t285 = -t335 * qJD(4) + t303 * t332;
t340 = cos(qJ(6));
t397 = t340 * t285;
t338 = sin(qJ(6));
t313 = t332 * t340 + t335 * t338;
t371 = t298 * t313;
t287 = qJD(4) * t332 + t303 * t335;
t396 = t285 * t338 - t287 * t340;
t391 = pkin(7) + t323;
t308 = t391 * t333;
t309 = t391 * t336;
t395 = -t393 * t308 - t339 * t309;
t348 = -t393 * t289 + t339 * t290;
t307 = t314 * qJD(4);
t296 = qJD(1) * t307;
t376 = t340 * t335;
t379 = t332 * t338;
t311 = -t376 + t379;
t372 = t298 * t311;
t394 = -t296 * t313 + t372 * t298;
t300 = t301 ^ 2;
t392 = pkin(8) * t335;
t390 = pkin(8) + qJ(5);
t249 = t287 * t338 + t397;
t388 = t249 * t303;
t387 = t396 * t303;
t321 = qJD(4) * t322;
t295 = -qJD(4) * t361 + t321;
t386 = t295 * t332;
t385 = t295 * t335;
t383 = t301 * t332;
t347 = t393 * t336 - t378;
t306 = t347 * qJD(4);
t382 = t306 * t332;
t381 = t314 * t332;
t380 = t314 * t335;
t349 = -qJD(6) * t287 - t386;
t367 = qJD(6) * t340;
t373 = -t285 * t367 + t295 * t376;
t222 = t349 * t338 + t373;
t375 = -t222 * t347 - t307 * t396;
t368 = qJD(6) * t314;
t231 = t313 * t306 + t367 * t380 - t368 * t379;
t271 = t313 * t314;
t374 = -t231 * t298 - t271 * t296;
t343 = t347 * qJD(3);
t235 = qJD(1) * t343 + (qJD(5) - t348) * qJD(4);
t246 = pkin(4) * t296 - qJ(5) * t295 - qJD(5) * t303;
t211 = t335 * t235 + t332 * t246;
t245 = qJD(4) * qJ(5) + t253;
t315 = -cos(pkin(9)) * pkin(1) - pkin(3) * t336 - pkin(2);
t299 = t315 * qJD(1) + qJD(3);
t261 = pkin(4) * t301 - qJ(5) * t303 + t299;
t218 = t335 * t245 + t332 * t261;
t275 = pkin(4) * t303 + qJ(5) * t301;
t228 = t332 * t275 - t335 * t348;
t254 = t395 * qJD(4) + t343;
t262 = pkin(4) * t307 - qJ(5) * t306 - qJD(5) * t314;
t221 = t335 * t254 + t332 * t262;
t268 = -pkin(4) * t347 - qJ(5) * t314 + t315;
t274 = -t339 * t308 + t393 * t309;
t233 = t332 * t268 + t335 * t274;
t212 = -pkin(8) * t285 + t218;
t369 = qJD(6) * t212;
t365 = MDP(14) * qJD(4);
t364 = qJD(1) * qJD(3);
t210 = -t235 * t332 + t335 * t246;
t217 = -t245 * t332 + t335 * t261;
t227 = t335 * t275 + t332 * t348;
t220 = -t254 * t332 + t335 * t262;
t232 = t335 * t268 - t274 * t332;
t207 = -pkin(8) * t386 + t211;
t208 = pkin(5) * t301 - pkin(8) * t287 + t217;
t358 = -qJD(6) * t208 - t207;
t236 = t333 * qJD(3) * t360 + t364 * t377 + t399;
t357 = -t311 * t296 - t371 * t298;
t203 = t208 * t340 - t212 * t338;
t204 = t208 * t338 + t212 * t340;
t356 = -t217 * t332 + t218 * t335;
t219 = -pkin(5) * t347 - pkin(8) * t380 + t232;
t225 = -pkin(8) * t381 + t233;
t355 = t219 * t340 - t225 * t338;
t354 = t219 * t338 + t225 * t340;
t223 = -qJD(6) * t396 + t313 * t295;
t353 = t223 * t347 - t249 * t307;
t230 = -t311 * t306 - t313 * t368;
t272 = t311 * t314;
t352 = -t230 * t298 + t272 * t296;
t351 = -t285 * t335 + t287 * t332;
t350 = (-t318 * t333 + t329) * t333 - t294 * t336;
t317 = t390 * t335;
t346 = pkin(5) * t303 + qJD(5) * t332 + qJD(6) * t317 + t301 * t392 + t227;
t316 = t390 * t332;
t345 = pkin(8) * t383 - qJD(5) * t335 + qJD(6) * t316 + t228;
t243 = -qJD(4) * pkin(4) + qJD(5) + t348;
t344 = t236 * t314 + t243 * t306 - t295 * t395;
t342 = -pkin(4) * t295 - qJ(5) * t296 + (-qJD(5) + t243) * t301;
t255 = t314 * qJD(3) + t274 * qJD(4);
t326 = -pkin(5) * t335 - pkin(4);
t258 = pkin(5) * t381 - t395;
t238 = pkin(5) * t382 + t255;
t237 = -pkin(5) * t383 + t253;
t229 = t285 * pkin(5) + t243;
t226 = pkin(5) * t386 + t236;
t214 = -pkin(8) * t382 + t221;
t209 = pkin(5) * t307 - t306 * t392 + t220;
t206 = pkin(5) * t296 - pkin(8) * t385 + t210;
t205 = t340 * t206;
t1 = [(t295 * t314 + t303 * t306) * MDP(9) + (t295 * t347 - t296 * t314 - t301 * t306 - t303 * t307) * MDP(10) + (t296 * t315 + t299 * t307) * MDP(14) + (t295 * t315 + t299 * t306) * MDP(15) + (-t210 * t347 + t217 * t307 + t220 * t301 + t232 * t296 + t255 * t285 + t344 * t332) * MDP(16) + (t211 * t347 - t218 * t307 - t221 * t301 - t233 * t296 + t255 * t287 + t344 * t335) * MDP(17) + (-t220 * t287 - t221 * t285 + (-t210 * t314 - t217 * t306 - t232 * t295) * t335 + (-t211 * t314 - t218 * t306 - t233 * t295) * t332) * MDP(18) + (t210 * t232 + t211 * t233 + t217 * t220 + t218 * t221 - t236 * t395 + t243 * t255) * MDP(19) + (-t222 * t272 - t230 * t396) * MDP(20) + (-t222 * t271 + t223 * t272 - t230 * t249 + t231 * t396) * MDP(21) + (-t352 + t375) * MDP(22) + (t353 + t374) * MDP(23) + (-t296 * t347 + t298 * t307) * MDP(24) + ((t209 * t340 - t214 * t338) * t298 + t355 * t296 - (-t207 * t338 + t205) * t347 + t203 * t307 + t238 * t249 + t258 * t223 + t226 * t271 + t229 * t231 + (t204 * t347 - t354 * t298) * qJD(6)) * MDP(25) + (-(t209 * t338 + t214 * t340) * t298 - t354 * t296 + (t206 * t338 + t207 * t340) * t347 - t204 * t307 - t238 * t396 + t258 * t222 - t226 * t272 + t229 * t230 + (t203 * t347 - t355 * t298) * qJD(6)) * MDP(26) + (0.2e1 * t398 + (t323 * t359 - t350) * MDP(8)) * qJD(3) + (t306 * MDP(11) - t307 * MDP(12) - t255 * MDP(14) - t254 * MDP(15)) * qJD(4); -t307 * t365 + (t285 * t307 - t296 * t381 - t347 * t386) * MDP(16) + (t287 * t307 - t296 * t380 - t347 * t385) * MDP(17) + (-t210 * t381 + t211 * t380 - t236 * t347 + t243 * t307) * MDP(19) + (-t353 + t374) * MDP(25) + (t352 + t375) * MDP(26) + (-qJD(4) * MDP(15) + t351 * MDP(18) + t356 * MDP(19) + (-t332 * MDP(16) - t335 * MDP(17)) * t301) * t306; 0.2e1 * t303 * t365 + (t321 + (-t301 - t361) * qJD(4)) * MDP(15) + (-t285 * t303 + t296 * t335 - t300 * t332) * MDP(16) + (-t287 * t303 - t296 * t332 - t300 * t335) * MDP(17) + (t351 * t301 + (-t332 ^ 2 - t335 ^ 2) * t295) * MDP(18) + (t210 * t335 + t211 * t332 - t243 * t303 + t356 * t301) * MDP(19) + (t357 - t388) * MDP(25) + (t387 + t394) * MDP(26) + (t350 * MDP(8) - t398) * qJD(1); -t300 * MDP(10) + (t321 + (t301 - t361) * qJD(4)) * MDP(11) + (-t236 + t399) * MDP(14) + (t299 * t301 - t347 * t364) * MDP(15) + (-t227 * t301 - t236 * t335 - t253 * t285 + t342 * t332) * MDP(16) + (t228 * t301 + t236 * t332 - t253 * t287 + t342 * t335) * MDP(17) + (t227 * t287 + t228 * t285 + (-qJD(5) * t285 - t217 * t301 + t211) * t335 + (qJD(5) * t287 - t218 * t301 - t210) * t332) * MDP(18) + (-pkin(4) * t236 - t217 * t227 - t218 * t228 - t243 * t253 + t356 * qJD(5) + (-t210 * t332 + t211 * t335) * qJ(5)) * MDP(19) + (t222 * t313 + t372 * t396) * MDP(20) + (-t222 * t311 - t223 * t313 + t372 * t249 + t371 * t396) * MDP(21) + (t387 - t394) * MDP(22) + (t357 + t388) * MDP(23) + ((-t316 * t340 - t317 * t338) * t296 + t326 * t223 + t226 * t311 - t237 * t249 + (t345 * t338 - t346 * t340) * t298 + t371 * t229) * MDP(25) + (-(-t316 * t338 + t317 * t340) * t296 + t326 * t222 + t226 * t313 + t237 * t396 + (t346 * t338 + t345 * t340) * t298 - t372 * t229) * MDP(26) + (MDP(10) * t303 - t299 * MDP(14) - t217 * MDP(16) + t218 * MDP(17) - t298 * MDP(24) - t203 * MDP(25) + t204 * MDP(26) + t301 * MDP(9)) * t303; (t287 * t301 + t386) * MDP(16) + (-t285 * t301 + t385) * MDP(17) + (-t285 ^ 2 - t287 ^ 2) * MDP(18) + (t217 * t287 + t218 * t285 + t236) * MDP(19) + (-t396 * t298 + t223) * MDP(25) + (-t298 * t397 + (-t287 * t298 + t349) * t338 + t373) * MDP(26); -t249 ^ 2 * MDP(21) + (t249 * t298 + t373) * MDP(22) + t296 * MDP(24) + (t204 * t298 + t205) * MDP(25) + (t203 * t298 + t229 * t249) * MDP(26) - (MDP(20) * t249 - MDP(21) * t396 + MDP(23) * t298 - t229 * MDP(25)) * t396 + (t349 * MDP(23) - MDP(25) * t369 + t358 * MDP(26)) * t340 + (t349 * MDP(22) + (qJD(6) * t285 - t385) * MDP(23) + t358 * MDP(25) + (-t206 + t369) * MDP(26)) * t338;];
tauc  = t1;
