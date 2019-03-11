% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:34
% EndTime: 2019-03-08 19:58:39
% DurationCPUTime: 2.50s
% Computational Cost: add. (1993->304), mult. (4984->441), div. (0->0), fcn. (3686->10), ass. (0->145)
t292 = cos(pkin(11));
t296 = sin(qJ(2));
t291 = sin(pkin(6));
t344 = qJD(1) * t291;
t327 = t296 * t344;
t277 = t292 * t327;
t290 = sin(pkin(11));
t299 = cos(qJ(2));
t326 = t299 * t344;
t246 = t290 * t326 + t277;
t295 = sin(qJ(4));
t298 = cos(qJ(4));
t313 = pkin(4) * t295 - pkin(9) * t298;
t273 = t313 * qJD(4);
t381 = -t246 + t273;
t380 = MDP(12) * t298;
t379 = MDP(6) * t295;
t251 = (t290 * t296 - t292 * t299) * t291;
t288 = t295 ^ 2;
t378 = (-t298 ^ 2 + t288) * MDP(7);
t276 = t290 * t327;
t249 = t292 * t326 - t276;
t307 = -pkin(4) * t298 - pkin(9) * t295 - pkin(3);
t373 = pkin(2) * t292;
t264 = t307 - t373;
t294 = sin(qJ(5));
t297 = cos(qJ(5));
t337 = qJD(5) * t297;
t357 = t297 * t298;
t377 = -t249 * t357 + t264 * t337 + t381 * t294;
t275 = qJD(2) * pkin(2) + t326;
t241 = t290 * t275 + t277;
t238 = qJD(2) * pkin(8) + t241;
t293 = cos(pkin(6));
t282 = qJD(1) * t293 + qJD(3);
t376 = -t295 * t238 + t282 * t298;
t340 = qJD(4) * t295;
t359 = t294 * t298;
t284 = pkin(2) * t290 + pkin(8);
t361 = t284 * t294;
t375 = t249 * t359 + t381 * t297 + t340 * t361;
t341 = qJD(4) * t294;
t343 = qJD(2) * t295;
t269 = t297 * t343 + t341;
t374 = t269 ^ 2;
t372 = pkin(5) * t295;
t371 = -qJ(6) - pkin(9);
t370 = qJ(6) * t295;
t248 = qJD(2) * t251;
t243 = qJD(1) * t248;
t339 = qJD(4) * t298;
t209 = t238 * t339 - t295 * t243 + t282 * t340;
t369 = t209 * t294;
t368 = t209 * t297;
t219 = -qJD(4) * pkin(4) - t376;
t367 = t219 * t294;
t366 = t219 * t297;
t342 = qJD(2) * t298;
t283 = -qJD(5) + t342;
t365 = t269 * t283;
t363 = t283 * t294;
t362 = t283 * t297;
t240 = t275 * t292 - t276;
t228 = t307 * qJD(2) - t240;
t360 = t294 * t228;
t358 = t295 * t297;
t223 = t298 * t238 + t295 * t282;
t220 = qJD(4) * pkin(9) + t223;
t202 = -t220 * t294 + t297 * t228;
t199 = -qJ(6) * t269 + t202;
t198 = -pkin(5) * t283 + t199;
t356 = t198 - t199;
t271 = t284 * t357;
t306 = -qJ(6) * t357 + t372;
t336 = qJD(6) * t297;
t355 = -t295 * t336 + t306 * qJD(4) + (-t271 + (-t264 + t370) * t294) * qJD(5) + t375;
t354 = (-qJ(6) * qJD(5) - qJD(4) * t284) * t358 + (-qJD(6) * t295 + (-qJ(6) * qJD(4) - qJD(5) * t284) * t298) * t294 + t377;
t272 = t313 * qJD(2);
t353 = t294 * t272 + t297 * t376;
t329 = qJD(4) * qJD(5);
t319 = t294 * t329;
t322 = t295 * t337;
t325 = t294 * t339;
t245 = t319 + (t322 + t325) * qJD(2);
t332 = t297 * qJD(4);
t267 = t294 * t343 - t332;
t324 = t298 * t332;
t352 = -t245 * t358 - t267 * t324;
t317 = qJD(5) * t371;
t350 = t336 - t353 + (qJ(6) * t342 + t317) * t294;
t314 = t297 * t272 - t294 * t376;
t349 = -t306 * qJD(2) - qJD(6) * t294 + t297 * t317 - t314;
t347 = t294 * t264 + t271;
t331 = qJD(2) * qJD(4);
t320 = t298 * t331;
t346 = (t320 + t329) * t297;
t338 = qJD(5) * t294;
t335 = t219 * qJD(5);
t323 = t295 * t338;
t244 = qJD(2) * t323 - t346;
t334 = t244 * MDP(20);
t333 = t294 * MDP(18);
t330 = qJD(4) * MDP(17);
t208 = qJD(4) * t376 - t243 * t298;
t252 = (t290 * t299 + t292 * t296) * t291;
t229 = (qJD(1) * t252 + t273) * qJD(2);
t328 = t297 * t208 + t228 * t337 + t294 * t229;
t321 = t288 * t331;
t318 = pkin(5) * t294 + t284;
t316 = t294 * t208 - t297 * t229;
t315 = t283 * t284 + t220;
t203 = t220 * t297 + t360;
t200 = -qJ(6) * t267 + t203;
t312 = -t198 * t297 - t200 * t294;
t311 = t198 * t294 - t200 * t297;
t236 = t252 * t298 + t293 * t295;
t215 = t236 * t297 + t251 * t294;
t214 = -t236 * t294 + t251 * t297;
t235 = t252 * t295 - t293 * t298;
t308 = qJD(2) * t288 - t283 * t298;
t201 = pkin(5) * t245 + t209;
t247 = qJD(2) * t252;
t242 = qJD(1) * t247;
t300 = qJD(4) ^ 2;
t305 = qJD(2) * t246 - t284 * t300 - t242;
t237 = -qJD(2) * pkin(3) - t240;
t304 = qJD(4) * (qJD(2) * (-pkin(3) - t373) + t237 + t249);
t303 = t220 * t338 - t328;
t302 = -t203 * qJD(5) - t316;
t301 = qJD(2) ^ 2;
t280 = t371 * t297;
t279 = t371 * t294;
t266 = t267 ^ 2;
t256 = t269 * t340;
t255 = t297 * t264;
t230 = -t294 * t370 + t347;
t227 = -qJ(6) * t358 + t255 + (-pkin(5) - t361) * t298;
t213 = -t235 * qJD(4) - t248 * t298;
t212 = t236 * qJD(4) - t248 * t295;
t211 = pkin(5) * t267 + qJD(6) + t219;
t197 = t214 * qJD(5) + t213 * t297 + t247 * t294;
t196 = -t215 * qJD(5) - t213 * t294 + t247 * t297;
t195 = -qJ(6) * t245 - qJD(6) * t267 - t303;
t194 = qJ(6) * t244 - qJD(6) * t269 + t331 * t372 + t302;
t1 = [(-t240 * t247 - t241 * t248 + t242 * t251 - t243 * t252) * MDP(5) + (-t196 * t283 + t212 * t267 + t235 * t245) * MDP(18) + (t197 * t283 + t212 * t269 - t235 * t244) * MDP(19) + (-t196 * t269 - t197 * t267 + t214 * t244 - t215 * t245) * MDP(20) + (t194 * t214 + t195 * t215 + t196 * t198 + t197 * t200 + t201 * t235 + t211 * t212) * MDP(21) + (-MDP(3) * t296 - MDP(4) * t299) * t301 * t291 + (-MDP(11) * t212 - MDP(12) * t213) * qJD(4) + ((-MDP(11) * t298 + MDP(12) * t295) * t247 + (t251 * t380 + (MDP(11) * t251 + MDP(18) * t214 - MDP(19) * t215) * t295) * qJD(4)) * qJD(2); (t240 * t246 - t241 * t249 + (-t242 * t292 - t243 * t290) * pkin(2)) * MDP(5) + 0.2e1 * t320 * t379 - 0.2e1 * t331 * t378 + (t295 * t304 + t305 * t298) * MDP(11) + (-t305 * t295 + t298 * t304) * MDP(12) + (-t244 * t358 + (-t323 + t324) * t269) * MDP(13) + (-t269 * t322 + (-t269 * t339 + (qJD(5) * t267 + t244) * t295) * t294 + t352) * MDP(14) + (t244 * t298 + t283 * t323 + t308 * t332 + t256) * MDP(15) + (t283 * t322 + t245 * t298 + (-t267 * t295 - t308 * t294) * qJD(4)) * MDP(16) + (-t283 - t342) * t295 * t330 + ((t264 * t338 - t375) * t283 + ((t267 * t284 + t367) * qJD(4) + (t315 * t297 + t360) * qJD(5) + t316) * t298 + (t297 * t335 + t369 + t284 * t245 - t249 * t267 + ((-t284 * t359 + t255) * qJD(2) + t202) * qJD(4)) * t295) * MDP(18) + (t377 * t283 + (-t315 * t338 + (t269 * t284 + t366) * qJD(4) + t328) * t298 + (-t294 * t335 + t368 - t284 * t244 - t249 * t269 + (-t347 * qJD(2) - t284 * t362 - t203) * qJD(4)) * t295) * MDP(19) + (t227 * t244 - t230 * t245 - t355 * t269 - t354 * t267 + t312 * t339 + (t311 * qJD(5) - t194 * t297 - t195 * t294) * t295) * MDP(20) + (t194 * t227 + t195 * t230 + t354 * t200 + t355 * t198 + t211 * t318 * t339 + (t201 * t318 + (pkin(5) * t337 - t249) * t211) * t295) * MDP(21) + (t298 * MDP(8) - t295 * MDP(9)) * t300; -t321 * t333 + (-t297 * t321 + t256) * MDP(19) + t352 * MDP(20) + (-t300 * MDP(12) - t245 * MDP(18) + t244 * MDP(19) - t201 * MDP(21) + (t269 * t294 * MDP(20) - t311 * MDP(21) + (t297 * MDP(19) + t333) * t283) * qJD(4)) * t298 + (-t300 * MDP(11) + qJD(4) * t267 * MDP(18) - t294 * t334 + (qJD(4) * t211 - t194 * t294 + t195 * t297) * MDP(21) + ((t267 * t294 + t269 * t297) * MDP(20) + t312 * MDP(21) + (t297 * MDP(18) - t294 * MDP(19)) * t283) * qJD(5)) * t295; (qJD(4) * t223 - t237 * t343 - t209) * MDP(11) + (-qJD(2) * t237 + t243) * t380 + (-t244 * t294 - t269 * t362) * MDP(13) + ((t267 * t283 - t244) * t297 + (-t245 + t365) * t294) * MDP(14) + (-t283 * t337 + (t283 * t357 + (-t269 + t341) * t295) * qJD(2)) * MDP(15) + (t283 * t338 + (-t283 * t359 + (t267 + t332) * t295) * qJD(2)) * MDP(16) + t283 * MDP(17) * t343 + (-pkin(4) * t245 - t368 + t314 * t283 - t223 * t267 + (pkin(9) * t362 + t367) * qJD(5) + (-t202 * t295 + (-pkin(9) * t340 - t219 * t298) * t294) * qJD(2)) * MDP(18) + (pkin(4) * t244 + t369 - t353 * t283 - t223 * t269 + (-pkin(9) * t363 + t366) * qJD(5) + (-t219 * t357 + (-pkin(9) * t332 + t203) * t295) * qJD(2)) * MDP(19) + (t244 * t279 + t245 * t280 - t349 * t269 - t350 * t267 + (t283 * t198 + t195) * t297 + (t283 * t200 - t194) * t294) * MDP(20) + (-t195 * t280 + t194 * t279 + t201 * (-pkin(5) * t297 - pkin(4)) + (-pkin(5) * t363 - t223) * t211 + t350 * t200 + t349 * t198) * MDP(21) + (-t298 * t379 + t378) * t301; (-t266 + t374) * MDP(14) + t346 * MDP(15) + (-t319 - t365) * MDP(16) + (-t203 * t283 - t219 * t269 + t302) * MDP(18) + (-t202 * t283 + t303) * MDP(19) + t356 * MDP(21) * t200 + (t334 + (-t211 * t269 + t194) * MDP(21)) * pkin(5) + (t269 * MDP(13) - t283 * MDP(15) + t219 * MDP(19) - t356 * MDP(20)) * t267 + (-MDP(16) * t325 + (t330 + (-t294 * MDP(15) - t297 * MDP(16)) * qJD(5)) * t295) * qJD(2); (-t266 - t374) * MDP(20) + (t198 * t269 + t200 * t267 + t201) * MDP(21);];
tauc  = t1;
