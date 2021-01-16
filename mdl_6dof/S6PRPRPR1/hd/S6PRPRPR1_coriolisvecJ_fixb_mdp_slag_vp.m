% Calculate Coriolis joint torque vector for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:06:27
% EndTime: 2021-01-16 01:06:36
% DurationCPUTime: 3.20s
% Computational Cost: add. (2258->304), mult. (5851->434), div. (0->0), fcn. (4641->12), ass. (0->143)
t291 = cos(pkin(11));
t295 = sin(qJ(2));
t290 = sin(pkin(6));
t348 = qJD(1) * t290;
t334 = t295 * t348;
t272 = t291 * t334;
t289 = sin(pkin(11));
t298 = cos(qJ(2));
t333 = t298 * t348;
t248 = t289 * t333 + t272;
t240 = qJD(2) * t248;
t252 = (t289 * t295 - t291 * t298) * t290;
t251 = qJD(2) * t252;
t294 = sin(qJ(4));
t297 = cos(qJ(4));
t377 = (t294 ^ 2 - t297 ^ 2) * MDP(7);
t369 = qJD(4) * pkin(4);
t376 = t294 * t369 - t248;
t292 = cos(pkin(6));
t277 = qJD(1) * t292 + qJD(3);
t270 = qJD(2) * pkin(2) + t333;
t238 = t289 * t270 + t272;
t236 = qJD(2) * pkin(8) + t238;
t326 = qJ(5) * qJD(2) + t236;
t222 = t297 * t277 - t326 * t294;
t356 = t294 * t277;
t223 = t326 * t297 + t356;
t241 = qJD(1) * t251;
t322 = qJD(2) * qJD(5) - t241;
t375 = -t223 * qJD(4) - t322 * t294;
t288 = sin(pkin(12));
t368 = cos(pkin(12));
t267 = t288 * t297 + t368 * t294;
t261 = t267 * qJD(2);
t203 = t222 * qJD(4) + t322 * t297;
t192 = t203 * t288 - t368 * t375;
t330 = t368 * t297;
t276 = qJD(2) * t330;
t347 = qJD(2) * t294;
t258 = t288 * t347 - t276;
t257 = qJD(6) + t258;
t280 = pkin(4) * t288 + pkin(9);
t339 = pkin(4) * t347;
t374 = (pkin(5) * t261 + pkin(9) * t258 + qJD(6) * t280 + t339) * t257 + t192;
t331 = t368 * t203;
t193 = t375 * t288 + t331;
t217 = t222 + t369;
t359 = t288 * t223;
t198 = t368 * t217 - t359;
t196 = -qJD(4) * pkin(5) - t198;
t271 = t289 * t334;
t237 = t270 * t291 - t271;
t335 = -pkin(4) * t297 - pkin(3);
t231 = t335 * qJD(2) + qJD(5) - t237;
t206 = pkin(5) * t258 - pkin(9) * t261 + t231;
t371 = pkin(2) * t291;
t274 = t335 - t371;
t306 = -t288 * t294 + t330;
t227 = -pkin(5) * t306 - pkin(9) * t267 + t274;
t263 = t306 * qJD(4);
t281 = pkin(2) * t289 + pkin(8);
t355 = qJ(5) + t281;
t265 = t355 * t297;
t329 = t355 * t294;
t229 = t368 * t265 - t288 * t329;
t260 = t267 * qJD(4);
t254 = qJD(2) * t260;
t320 = t192 * t267 - t229 * t254;
t328 = qJD(4) * t355;
t245 = qJD(5) * t297 - t294 * t328;
t250 = t291 * t333 - t271;
t305 = -qJD(5) * t294 - t297 * t328;
t353 = -t368 * t245 + t306 * t250 - t288 * t305;
t373 = (qJD(6) * t206 + t193) * t306 + t196 * t263 + (-qJD(6) * t227 + t353) * t257 + t320;
t372 = t294 * MDP(11) + t297 * MDP(12);
t370 = pkin(4) * t294;
t293 = sin(qJ(6));
t343 = qJD(6) * t293;
t340 = qJD(2) * qJD(4);
t332 = t294 * t340;
t275 = t288 * t332;
t255 = qJD(4) * t276 - t275;
t296 = cos(qJ(6));
t341 = t296 * qJD(4);
t350 = qJD(6) * t341 + t296 * t255;
t220 = -t261 * t343 + t350;
t367 = t220 * t293;
t366 = t227 * t254;
t360 = t261 * t293;
t242 = -t341 + t360;
t365 = t242 * t257;
t364 = t242 * t261;
t244 = qJD(4) * t293 + t261 * t296;
t363 = t244 * t257;
t362 = t244 * t261;
t361 = t255 * t293;
t357 = t293 * t254;
t246 = t296 * t254;
t354 = -t220 * t306 + t244 * t260;
t215 = t368 * t223;
t199 = t288 * t217 + t215;
t352 = t245 * t288 - t267 * t250 - t368 * t305;
t351 = pkin(5) * t260 - pkin(9) * t263 + t376;
t346 = qJD(2) * t297;
t345 = qJD(4) * t297;
t344 = qJD(6) * t267;
t337 = t267 * t357;
t336 = t267 * t246;
t233 = pkin(4) * t332 + t240;
t327 = t296 * t257;
t197 = qJD(4) * pkin(9) + t199;
t191 = t197 * t296 + t206 * t293;
t319 = t197 * t293 - t206 * t296;
t253 = (t289 * t298 + t291 * t295) * t290;
t234 = t253 * t297 + t292 * t294;
t315 = -t253 * t294 + t292 * t297;
t208 = t368 * t234 + t288 * t315;
t318 = t208 * t296 + t252 * t293;
t317 = -t208 * t293 + t252 * t296;
t221 = t244 * qJD(6) + t361;
t316 = t221 * t306 - t242 * t260;
t311 = t246 + (-t258 * t293 - t343) * t257;
t310 = -t263 * t293 - t296 * t344;
t309 = -t263 * t296 + t267 * t343;
t299 = qJD(4) ^ 2;
t308 = t281 * t299;
t235 = -qJD(2) * pkin(3) - t237;
t307 = qJD(4) * (qJD(2) * (-pkin(3) - t371) + t235 + t250);
t201 = t368 * t222 - t359;
t304 = -t280 * t254 + (t196 + t201) * t257;
t303 = -t234 * qJD(4) + t294 * t251;
t300 = qJD(2) ^ 2;
t282 = -t368 * pkin(4) - pkin(5);
t249 = qJD(2) * t253;
t228 = t265 * t288 + t368 * t329;
t211 = t315 * qJD(4) - t251 * t297;
t207 = t234 * t288 - t368 * t315;
t205 = pkin(5) * t254 - pkin(9) * t255 + t233;
t204 = t296 * t205;
t200 = t222 * t288 + t215;
t195 = t368 * t211 + t288 * t303;
t194 = t211 * t288 - t368 * t303;
t1 = [(-t237 * t249 - t238 * t251 + t240 * t252 - t241 * t253) * MDP(5) + (-t249 * t346 + (-t253 * t345 + (-qJD(4) * t292 + 0.2e1 * t251) * t294) * qJD(4)) * MDP(11) + (-qJD(4) * t211 + (t249 * t294 + t252 * t345) * qJD(2)) * MDP(12) + (-qJD(4) * t194 + t249 * t258 + t252 * t254) * MDP(13) + (-qJD(4) * t195 + t249 * t261 + t252 * t255) * MDP(14) + (t194 * t261 - t195 * t258 + t207 * t255 - t208 * t254) * MDP(15) + (t192 * t207 + t193 * t208 - t194 * t198 + t195 * t199 + t231 * t249 + t233 * t252) * MDP(16) + ((-t318 * qJD(6) - t195 * t293 + t249 * t296) * t257 + t317 * t254 + t194 * t242 + t207 * t221) * MDP(22) + (-(t317 * qJD(6) + t195 * t296 + t249 * t293) * t257 - t318 * t254 + t194 * t244 + t207 * t220) * MDP(23) + (-t295 * MDP(3) - t298 * MDP(4)) * t290 * t300; (t237 * t248 - t238 * t250 + (-t240 * t291 - t241 * t289) * pkin(2)) * MDP(5) - 0.2e1 * t340 * t377 + (t231 * t260 - t233 * t306 - t248 * t258 + t254 * t274 + (t258 * t370 - t352) * qJD(4)) * MDP(13) + (t231 * t263 + t233 * t267 - t248 * t261 + t255 * t274 + (t261 * t370 + t353) * qJD(4)) * MDP(14) + (t193 * t306 - t198 * t263 - t199 * t260 + t228 * t255 + t353 * t258 + t352 * t261 + t320) * MDP(15) + (t192 * t228 + t193 * t229 - t352 * t198 - t353 * t199 + t376 * t231 + t233 * t274) * MDP(16) + (t220 * t267 * t296 - t309 * t244) * MDP(17) + ((-t242 * t296 - t244 * t293) * t263 + (-t367 - t221 * t296 + (t242 * t293 - t244 * t296) * qJD(6)) * t267) * MDP(18) + (-t309 * t257 + t336 + t354) * MDP(19) + (t310 * t257 + t316 - t337) * MDP(20) + (-t254 * t306 + t257 * t260) * MDP(21) + (-t319 * t260 - t204 * t306 + t228 * t221 + t352 * t242 + (t366 + t351 * t257 + (t196 * t267 + t197 * t306 - t229 * t257) * qJD(6)) * t296 + t373 * t293) * MDP(22) + (-t191 * t260 + t228 * t220 + t352 * t244 + (-t366 + (-qJD(6) * t197 + t205) * t306 - t196 * t344 + (qJD(6) * t229 - t351) * t257) * t293 + t373 * t296) * MDP(23) + (t307 * MDP(11) + t308 * MDP(12) - t299 * MDP(9)) * t294 + (-t308 * MDP(11) + t307 * MDP(12) + 0.2e1 * MDP(6) * t332 + t299 * MDP(8)) * t297; (-t254 * t267 - t255 * t306 - t258 * t263 + t260 * t261) * MDP(15) + (-t192 * t306 + t193 * t267 - t198 * t260 + t199 * t263) * MDP(16) + (-t316 - t337) * MDP(22) + (-t336 + t354) * MDP(23) - t372 * t299 + (t310 * MDP(22) + t309 * MDP(23)) * t257 + (-t260 * MDP(13) - t263 * MDP(14)) * qJD(4); (qJD(4) * t200 - t231 * t261 - t258 * t339 - t192) * MDP(13) + (-t331 + t231 * t258 + (-qJD(2) * pkin(4) * t261 + t288 * t322) * t294 + (-t288 * (-qJ(5) * t346 - t236 * t297 - t356) + t201) * qJD(4)) * MDP(14) + ((t199 - t200) * t261 + (-t198 + t201) * t258 + (-t254 * t288 - t368 * t255) * pkin(4)) * MDP(15) + (t198 * t200 - t199 * t201 + (-t368 * t192 + t193 * t288 - t231 * t347) * pkin(4)) * MDP(16) + (t244 * t327 + t367) * MDP(17) + ((t220 - t365) * t296 + (-t221 - t363) * t293) * MDP(18) + (t257 * t327 + t357 - t362) * MDP(19) + (t311 + t364) * MDP(20) - t257 * t261 * MDP(21) + (-t200 * t242 + t282 * t221 + t261 * t319 + t304 * t293 - t374 * t296) * MDP(22) + (t191 * t261 - t200 * t244 + t282 * t220 + t374 * t293 + t304 * t296) * MDP(23) + t372 * (-qJD(2) * t235 + t241) + (-t294 * t297 * MDP(6) + t377) * t300; 0.2e1 * t261 * qJD(4) * MDP(13) + (-t275 + (t276 - t258) * qJD(4)) * MDP(14) + (-t258 ^ 2 - t261 ^ 2) * MDP(15) + (t198 * t261 + t199 * t258 + t233) * MDP(16) + (t311 - t364) * MDP(22) + (-t257 ^ 2 * t296 - t357 - t362) * MDP(23); t244 * t242 * MDP(17) + (-t242 ^ 2 + t244 ^ 2) * MDP(18) + (t350 + t365) * MDP(19) + (-t361 + t363) * MDP(20) + t254 * MDP(21) + (t191 * t257 - t193 * t293 - t196 * t244 + t204) * MDP(22) + (-t193 * t296 + t196 * t242 - t205 * t293 - t257 * t319) * MDP(23) + (-MDP(19) * t360 - t244 * MDP(20) - t191 * MDP(22) + t319 * MDP(23)) * qJD(6);];
tauc = t1;
