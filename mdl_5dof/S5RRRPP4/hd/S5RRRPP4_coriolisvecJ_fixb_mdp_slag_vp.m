% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRPP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:50
% EndTime: 2019-12-31 20:55:54
% DurationCPUTime: 1.87s
% Computational Cost: add. (2921->262), mult. (7671->342), div. (0->0), fcn. (5297->6), ass. (0->139)
t308 = cos(qJ(3));
t309 = cos(qJ(2));
t351 = qJD(1) * t309;
t343 = t308 * t351;
t306 = sin(qJ(3));
t307 = sin(qJ(2));
t352 = qJD(1) * t307;
t344 = t306 * t352;
t262 = -t343 + t344;
t264 = -t306 * t351 - t308 * t352;
t304 = sin(pkin(8));
t305 = cos(pkin(8));
t237 = -t305 * t262 + t264 * t304;
t301 = qJD(2) + qJD(3);
t376 = t237 * t301;
t326 = -t262 * t304 - t305 * t264;
t375 = t326 ^ 2;
t347 = qJD(1) * qJD(2);
t374 = -0.2e1 * t347;
t373 = MDP(4) * t307;
t372 = MDP(5) * (t307 ^ 2 - t309 ^ 2);
t369 = pkin(6) + pkin(7);
t283 = t369 * t307;
t277 = qJD(1) * t283;
t367 = qJD(2) * pkin(2);
t271 = -t277 + t367;
t345 = qJD(2) * t369;
t331 = qJD(1) * t345;
t272 = t307 * t331;
t371 = (qJD(3) * t271 - t272) * t308;
t259 = t264 * qJ(4);
t284 = t369 * t309;
t279 = qJD(1) * t284;
t265 = t306 * t279;
t354 = -t308 * t277 - t265;
t231 = t259 + t354;
t269 = t308 * t279;
t333 = t277 * t306 - t269;
t366 = qJ(4) * t262;
t322 = t333 + t366;
t346 = pkin(2) * t304 * t306;
t349 = qJD(3) * t308;
t355 = -qJD(3) * t346 - t304 * t322 + (pkin(2) * t349 - t231) * t305;
t336 = t308 * t271 - t265;
t228 = t259 + t336;
t276 = t306 * t309 + t307 * t308;
t370 = qJD(1) * t276;
t246 = t301 * t276;
t241 = t246 * qJD(1);
t273 = t309 * t331;
t350 = qJD(3) * t306;
t334 = -t306 * t273 - t279 * t350;
t200 = -qJ(4) * t241 - qJD(4) * t262 + t334 + t371;
t342 = t309 * t347;
t240 = qJD(3) * t343 - t301 * t344 + t308 * t342;
t325 = -t271 * t306 - t269;
t335 = t306 * t272 - t308 * t273;
t315 = qJD(3) * t325 + t335;
t313 = -qJ(4) * t240 + qJD(4) * t264 + t315;
t189 = t200 * t304 - t305 * t313;
t297 = -pkin(2) * t309 - pkin(1);
t282 = t297 * qJD(1);
t247 = pkin(3) * t262 + qJD(4) + t282;
t208 = -pkin(4) * t237 - qJ(5) * t326 + t247;
t321 = -t208 * t326 - t189;
t368 = pkin(3) * t264;
t275 = t306 * t307 - t308 * t309;
t324 = t283 * t306 - t284 * t308;
t233 = -qJ(4) * t275 - t324;
t317 = -qJ(4) * t276 - t283 * t308 - t284 * t306;
t215 = t233 * t304 - t305 * t317;
t365 = t189 * t215;
t229 = -t325 - t366;
t225 = t305 * t229;
t205 = t228 * t304 + t225;
t364 = t205 * t326;
t363 = t229 * t304;
t362 = t282 * t264;
t361 = t305 * t306;
t310 = qJD(2) ^ 2;
t360 = t307 * t310;
t359 = t309 * t310;
t311 = qJD(1) ^ 2;
t358 = t309 * t311;
t190 = t305 * t200 + t304 * t313;
t357 = -qJD(5) - t355;
t224 = pkin(3) * t301 + t228;
t204 = t304 * t224 + t225;
t356 = t231 * t304 - t305 * t322 - (t304 * t308 + t361) * qJD(3) * pkin(2);
t296 = pkin(2) * t308 + pkin(3);
t258 = pkin(2) * t361 + t304 * t296;
t206 = t228 * t305 - t363;
t348 = qJD(5) - t206;
t299 = t307 * t367;
t298 = pkin(2) * t352;
t341 = -pkin(2) * t301 - t271;
t340 = pkin(3) * t241 + qJD(2) * t298;
t339 = pkin(3) * t246 + t299;
t338 = t356 * t326;
t337 = pkin(1) * t374;
t217 = t240 * t304 + t305 * t241;
t330 = t208 * t237 + t190;
t203 = t224 * t305 - t363;
t201 = -pkin(4) * t301 + qJD(5) - t203;
t202 = qJ(5) * t301 + t204;
t329 = -t201 * t237 + t202 * t326;
t328 = t203 * t237 + t204 * t326;
t218 = t240 * t305 - t241 * t304;
t323 = pkin(3) * t275 + t297;
t257 = t296 * t305 - t346;
t320 = t282 * t262 - t334;
t319 = -t264 * t262 * MDP(11) + t240 * MDP(13) + (-t262 ^ 2 + t264 ^ 2) * MDP(12) + (t262 * MDP(13) + (-t264 - t370) * MDP(14)) * t301;
t278 = t307 * t345;
t280 = t309 * t345;
t318 = -t308 * t278 - t306 * t280 - t283 * t349 - t284 * t350;
t213 = pkin(4) * t326 - qJ(5) * t237 - t368;
t209 = -qJ(4) * t246 - qJD(4) * t275 + t318;
t245 = t301 * t275;
t314 = qJD(3) * t324 + t278 * t306 - t308 * t280;
t312 = qJ(4) * t245 - qJD(4) * t276 + t314;
t193 = t209 * t304 - t305 * t312;
t194 = t305 * t209 + t304 * t312;
t216 = t305 * t233 + t304 * t317;
t244 = -t275 * t304 + t276 * t305;
t316 = t189 * t244 + t193 * t326 + t194 * t237 + t215 * t218 - t216 * t217;
t192 = pkin(4) * t217 - qJ(5) * t218 - qJD(5) * t326 + t340;
t300 = t301 * qJD(5);
t294 = -pkin(3) * t305 - pkin(4);
t293 = pkin(3) * t304 + qJ(5);
t252 = -pkin(4) - t257;
t251 = qJ(5) + t258;
t243 = t305 * t275 + t276 * t304;
t220 = -t245 * t305 - t246 * t304;
t219 = -t245 * t304 + t305 * t246;
t214 = pkin(4) * t243 - qJ(5) * t244 + t323;
t212 = t213 + t298;
t195 = pkin(4) * t219 - qJ(5) * t220 - qJD(5) * t244 + t339;
t188 = t300 + t190;
t1 = [0.2e1 * t342 * t373 + t372 * t374 + MDP(6) * t359 - MDP(7) * t360 + (-pkin(6) * t359 + t307 * t337) * MDP(9) + (pkin(6) * t360 + t309 * t337) * MDP(10) + (t240 * t276 + t245 * t264) * MDP(11) + (-t240 * t275 - t241 * t276 + t245 * t262 + t246 * t264) * MDP(12) + (t297 * t241 + t282 * t246 + (qJD(1) * t275 + t262) * t299) * MDP(16) + (t297 * t240 - t282 * t245 + (-t264 + t370) * t299) * MDP(17) + (-t190 * t243 - t203 * t220 - t204 * t219 + t316) * MDP(18) + (t190 * t216 - t203 * t193 + t204 * t194 + t247 * t339 + t340 * t323 + t365) * MDP(19) + (t192 * t243 - t195 * t237 + t208 * t219 + t214 * t217) * MDP(20) + (-t188 * t243 + t201 * t220 - t202 * t219 + t316) * MDP(21) + (-t192 * t244 - t195 * t326 - t208 * t220 - t214 * t218) * MDP(22) + (t188 * t216 + t192 * t214 + t193 * t201 + t194 * t202 + t195 * t208 + t365) * MDP(23) + (-t245 * MDP(13) - t246 * MDP(14) + t314 * MDP(16) - t318 * MDP(17) - t193 * MDP(20) + t194 * MDP(22)) * t301; -t358 * t373 + t311 * t372 + (-t262 * t298 + t362 - t333 * t301 + (t306 * t341 - t269) * qJD(3) + t335) * MDP(16) + (t264 * t298 + t354 * t301 + (qJD(3) * t341 + t272) * t308 + t320) * MDP(17) + (-t217 * t258 - t218 * t257 + t237 * t355 + t328 - t338) * MDP(18) + (t190 * t258 - t189 * t257 - t247 * (t298 - t368) + t355 * t204 + t356 * t203) * MDP(19) + (t212 * t237 + t301 * t356 + t321) * MDP(20) + (-t217 * t251 + t218 * t252 - t237 * t357 + t329 - t338) * MDP(21) + (t212 * t326 - t301 * t357 + t300 + t330) * MDP(22) + (t188 * t251 + t189 * t252 - t201 * t356 - t202 * t357 - t208 * t212) * MDP(23) + t319 + (MDP(9) * t307 * t311 + MDP(10) * t358) * pkin(1); (-t301 * t325 + t315 + t362) * MDP(16) + (t301 * t336 + t320 - t371) * MDP(17) + (-t364 - t206 * t237 + (-t217 * t304 - t218 * t305) * pkin(3) + t328) * MDP(18) + (t203 * t205 - t204 * t206 + (-t189 * t305 + t190 * t304 + t247 * t264) * pkin(3)) * MDP(19) + (t205 * t301 + t213 * t237 + t321) * MDP(20) + (-t217 * t293 + t218 * t294 + t237 * t348 + t329 - t364) * MDP(21) + (-t206 * t301 + t213 * t326 + 0.2e1 * t300 + t330) * MDP(22) + (t188 * t293 + t189 * t294 - t201 * t205 + t202 * t348 - t208 * t213) * MDP(23) + t319; (t203 * t326 - t204 * t237 + t340) * MDP(19) + (t301 * t326 + t217) * MDP(20) + (-t218 - t376) * MDP(22) + (-t201 * t326 - t202 * t237 + t192) * MDP(23) + (MDP(18) + MDP(21)) * (-t237 ^ 2 - t375); -t326 * t237 * MDP(20) + (t218 - t376) * MDP(21) + (-t301 ^ 2 - t375) * MDP(22) + (-t202 * t301 - t321) * MDP(23);];
tauc = t1;
