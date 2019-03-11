% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:39
% EndTime: 2019-03-09 08:35:46
% DurationCPUTime: 3.09s
% Computational Cost: add. (1950->361), mult. (4119->475), div. (0->0), fcn. (2106->4), ass. (0->165)
t298 = sin(qJ(2));
t357 = qJD(1) * t298;
t278 = pkin(7) * t357;
t254 = -qJ(4) * t357 + t278;
t334 = qJD(3) + t254;
t301 = -pkin(2) - pkin(3);
t338 = qJD(2) * t301;
t238 = t338 + t334;
t344 = qJD(1) * qJD(2);
t392 = -0.2e1 * t344;
t292 = t298 ^ 2;
t300 = cos(qJ(2));
t293 = t300 ^ 2;
t391 = (t292 - t293) * MDP(5);
t356 = qJD(1) * t300;
t279 = pkin(7) * t356;
t256 = -qJ(4) * t356 + t279;
t291 = qJD(2) * qJ(3);
t246 = -t256 - t291;
t270 = qJD(5) + t357;
t390 = t270 * MDP(21);
t333 = t298 * t344;
t267 = qJ(4) * t333;
t289 = qJD(2) * qJD(3);
t354 = qJD(2) * t298;
t341 = pkin(7) * t354;
t313 = qJD(4) * t300 + t341;
t227 = qJD(1) * t313 - t267 - t289;
t299 = cos(qJ(5));
t297 = sin(qJ(5));
t355 = qJD(2) * t297;
t252 = t299 * t356 + t355;
t389 = t252 ^ 2;
t388 = pkin(5) * t297;
t387 = pkin(7) - qJ(4);
t296 = qJ(3) + pkin(4);
t386 = qJD(2) * pkin(2);
t385 = qJ(6) * t300;
t384 = t227 * t297;
t383 = t227 * t299;
t342 = qJD(2) * qJD(5);
t272 = t299 * t342;
t349 = qJD(5) * t300;
t336 = t297 * t349;
t347 = t299 * qJD(2);
t337 = t298 * t347;
t309 = t336 + t337;
t229 = qJD(1) * t309 - t272;
t382 = t229 * t297;
t265 = t297 * t333;
t230 = qJD(5) * t252 - t265;
t381 = t230 * t299;
t251 = t297 * t356 - t347;
t380 = t251 * t297;
t379 = t251 * t299;
t378 = t252 * t270;
t377 = t252 * t297;
t376 = t252 * t299;
t375 = t270 * t297;
t374 = t270 * t298;
t373 = t270 * t299;
t372 = t298 * t299;
t302 = qJD(2) ^ 2;
t371 = t298 * t302;
t263 = t387 * t298;
t253 = t299 * t263;
t370 = t299 * t300;
t369 = t300 * t302;
t303 = qJD(1) ^ 2;
t368 = t300 * t303;
t290 = -pkin(8) + t301;
t367 = qJ(6) - t290;
t248 = -qJD(1) * pkin(1) - pkin(2) * t356 - qJ(3) * t357;
t234 = pkin(3) * t356 + qJD(4) - t248;
t322 = pkin(4) * t298 + pkin(8) * t300;
t219 = qJD(1) * t322 + t234;
t232 = qJD(2) * t290 + t334;
t203 = t299 * t219 - t232 * t297;
t201 = qJ(6) * t252 + t203;
t200 = pkin(5) * t270 + t201;
t366 = t200 - t201;
t327 = qJD(5) * t367;
t275 = qJ(3) * t356;
t311 = pkin(4) * t300 + t290 * t298;
t222 = qJD(1) * t311 + t275;
t363 = t297 * t222 + t299 * t256;
t365 = -qJD(6) * t299 - t363 + (qJ(6) * t357 + t327) * t297;
t315 = pkin(5) * t300 - qJ(6) * t372;
t329 = t299 * t222 - t256 * t297;
t364 = -qJD(1) * t315 + qJD(6) * t297 + t299 * t327 - t329;
t261 = -t300 * pkin(2) - t298 * qJ(3) - pkin(1);
t249 = t300 * pkin(3) - t261;
t231 = t322 + t249;
t362 = t297 * t231 + t253;
t282 = t298 * qJD(3);
t332 = t300 * t344;
t360 = qJ(3) * t332 + qJD(1) * t282;
t353 = qJD(2) * t300;
t359 = qJ(3) * t353 + t282;
t352 = qJD(4) * t298;
t351 = qJD(5) * t297;
t350 = qJD(5) * t299;
t348 = qJD(6) * t300;
t345 = -qJD(4) - t234;
t343 = qJD(2) * MDP(23);
t340 = t298 * t368;
t308 = t311 * qJD(2);
t218 = t308 + t359;
t244 = t353 * t387 - t352;
t339 = t297 * t218 + t231 * t350 + t299 * t244;
t335 = t299 * t349;
t331 = pkin(1) * t392;
t330 = qJD(3) - t386;
t328 = t299 * t231 - t263 * t297;
t323 = t298 * t338;
t223 = qJD(1) * t323 + t360;
t233 = t323 + t359;
t326 = qJD(1) * t233 + t223;
t325 = qJD(1) * t249 + t234;
t324 = qJD(1) * t261 + t248;
t204 = t219 * t297 + t232 * t299;
t212 = qJD(1) * t308 + t360;
t211 = t299 * t212;
t269 = pkin(7) * t332;
t235 = t269 + (-qJ(4) * t353 - t352) * qJD(1);
t307 = -qJD(5) * t204 - t235 * t297 + t211;
t196 = pkin(5) * t332 - qJ(6) * t229 + qJD(6) * t252 + t307;
t310 = -t297 * t212 - t219 * t350 + t232 * t351 - t299 * t235;
t197 = qJ(6) * t230 + qJD(6) * t251 - t310;
t321 = t196 * t299 + t197 * t297;
t202 = qJ(6) * t251 + t204;
t320 = -t200 * t299 - t202 * t297;
t319 = -t200 * t297 + t202 * t299;
t318 = qJD(1) * t293 - t374;
t241 = qJD(2) * pkin(4) - t246;
t217 = -pkin(5) * t251 + qJD(6) + t241;
t317 = -t246 * MDP(18) + t217 * MDP(27);
t316 = -t297 * MDP(24) - t299 * MDP(25);
t236 = pkin(2) * t333 - t360;
t245 = pkin(2) * t354 - t359;
t314 = -pkin(7) * t302 - qJD(1) * t245 - t236;
t312 = -t241 * t298 - t290 * t353;
t208 = -pkin(5) * t230 - t227;
t259 = -pkin(7) * t333 + t289;
t260 = t278 + t330;
t262 = t279 + t291;
t306 = t259 * t300 + (t260 * t300 + (-t262 + t279) * t298) * qJD(2);
t305 = (-t377 + t379) * MDP(26) + t319 * MDP(27) + t316 * t270;
t304 = (-t376 - t380) * MDP(26) + t320 * MDP(27) + (-t299 * MDP(24) + t297 * MDP(25)) * t270;
t286 = 0.2e1 * t289;
t285 = t300 * qJ(4);
t276 = qJ(4) * t354;
t264 = pkin(7) * t300 - t285;
t258 = t367 * t299;
t257 = t367 * t297;
t255 = pkin(2) * t357 - t275;
t250 = t251 ^ 2;
t243 = -t276 + t313;
t242 = t301 * t357 + t275;
t215 = t299 * t218;
t209 = t297 * t385 + t362;
t207 = pkin(5) * t298 + qJ(6) * t370 + t328;
t199 = qJ(6) * t335 + (-qJ(6) * t354 - qJD(5) * t263 + t348) * t297 + t339;
t198 = t299 * t348 - t244 * t297 + t215 + t315 * qJD(2) + (-t253 + (-t231 - t385) * t297) * qJD(5);
t1 = [0.2e1 * t298 * MDP(4) * t332 + t391 * t392 + MDP(6) * t369 - MDP(7) * t371 + (-pkin(7) * t369 + t298 * t331) * MDP(9) + (pkin(7) * t371 + t300 * t331) * MDP(10) + (t300 * t314 + t324 * t354) * MDP(11) + t306 * MDP(12) + (t298 * t314 - t324 * t353) * MDP(13) + (pkin(7) * t306 + t236 * t261 + t245 * t248) * MDP(14) + (t326 * t298 + (t300 * t325 - t243) * qJD(2)) * MDP(15) + (-t326 * t300 + (t298 * t325 + t244) * qJD(2)) * MDP(16) + (t227 * t300 - t235 * t298 + (-t238 * t300 - t246 * t298) * qJD(2) + (t243 * t300 - t244 * t298 + (-t263 * t300 + t264 * t298) * qJD(2)) * qJD(1)) * MDP(17) + (t223 * t249 - t227 * t264 + t233 * t234 + t235 * t263 + t238 * t244 + t243 * t246) * MDP(18) + (-t229 * t370 - t252 * t309) * MDP(19) + ((t377 + t379) * t354 + (t382 - t381 + (-t376 + t380) * qJD(5)) * t300) * MDP(20) + (t270 * t336 + t229 * t298 + (-t252 * t300 - t299 * t318) * qJD(2)) * MDP(21) + (t270 * t335 + t230 * t298 + (t251 * t300 + t297 * t318) * qJD(2)) * MDP(22) + (t270 + t357) * t300 * t343 + ((-t263 * t350 + t215) * t270 + (-t232 * t350 + t211) * t298 + t243 * t251 - t264 * t230 + ((-qJD(5) * t231 - t244) * t270 + (qJD(2) * t241 - qJD(5) * t219 - t235) * t298) * t297 + (-t241 * t350 + t384 + (qJD(1) * t328 + t203) * qJD(2)) * t300) * MDP(24) + (-(-t263 * t351 + t339) * t270 + t243 * t252 + t264 * t229 + (t241 * t347 + t310) * t298 + (t241 * t351 + t383 + (-qJD(1) * t362 - t204) * qJD(2)) * t300) * MDP(25) + (t198 * t252 + t199 * t251 - t207 * t229 + t209 * t230 + t320 * t354 + (qJD(5) * t319 + t321) * t300) * MDP(26) + (t197 * t209 + t202 * t199 + t196 * t207 + t200 * t198 - t208 * t285 + t217 * (t354 * t388 + t276 - t341) + (t208 * (pkin(7) - t388) + t217 * (-pkin(5) * t350 - qJD(4))) * t300) * MDP(27); -MDP(4) * t340 + t303 * t391 + t286 * MDP(13) + (qJ(3) * t259 + qJD(3) * t262 - t248 * t255) * MDP(14) + (qJD(2) * t254 + t267 + t286) * MDP(15) + (-qJD(2) * t256 + t269) * MDP(16) + (-qJ(3) * t227 - t234 * t242 + t235 * t301 - t238 * t256 - t246 * t334) * MDP(18) + (t252 * t373 - t382) * MDP(19) + ((-t251 * t270 - t229) * t299 + (-t230 - t378) * t297) * MDP(20) - t350 * t390 + t270 * t351 * MDP(22) - t270 * MDP(23) * t356 + (-t296 * t230 - t383 - t329 * t270 - t334 * t251 + (-t241 * t297 - t290 * t373) * qJD(5)) * MDP(24) + (t296 * t229 + t384 + t363 * t270 - t334 * t252 + (-t241 * t299 + t290 * t375) * qJD(5)) * MDP(25) + (-t229 * t257 - t230 * t258 + t364 * t252 + t365 * t251 + (t200 * t270 - t197) * t299 + (t202 * t270 + t196) * t297) * MDP(26) + (-t197 * t258 + t196 * t257 + t208 * (pkin(5) * t299 + t296) + (-pkin(5) * t375 + t334) * t217 + t365 * t202 + t364 * t200) * MDP(27) + (MDP(9) * t298 * t303 + MDP(10) * t368) * pkin(1) + ((-t248 * t298 + t255 * t300) * MDP(11) + ((t262 - t291) * t298 + (-t260 + t330) * t300) * MDP(12) + (t248 * t300 + t255 * t298) * MDP(13) + (t262 * t298 + (-t260 - t386) * t300) * pkin(7) * MDP(14) + (t345 * t300 + (-pkin(7) * qJD(2) - t242) * t298) * MDP(15) + ((-qJ(4) * qJD(2) + t242) * t300 + t345 * t298) * MDP(16) + (-t270 * t372 + (t252 - t355) * t300) * MDP(21) + (t297 * t374 + (-t251 - t347) * t300) * MDP(22) + (-t203 * t300 + t297 * t312) * MDP(24) + (t204 * t300 + t299 * t312) * MDP(25)) * qJD(1); (t381 + t382) * MDP(26) + (-t196 * t297 + t197 * t299) * MDP(27) + (-MDP(11) + MDP(16)) * t340 + (MDP(14) + MDP(18)) * t269 + (MDP(13) + MDP(15)) * (-t292 * t303 - t302) + (-t262 * MDP(14) + t251 * MDP(24) + t252 * MDP(25) - t317) * qJD(2) + t304 * qJD(5) + ((-MDP(18) * qJ(4) + t316) * t353 + (t248 * MDP(14) + MDP(18) * t345 + t304) * t298) * qJD(1); t360 * MDP(18) + (-t229 * t299 + t230 * t297) * MDP(26) + t321 * MDP(27) + (-t292 - t293) * MDP(17) * t303 + t305 * qJD(5) + ((0.2e1 * qJD(2) * MDP(15) + (-t251 + t347) * MDP(24) + (-t252 - t355) * MDP(25) + t317) * t300 + (t238 * MDP(18) + (MDP(18) * t301 + 0.2e1 * MDP(16)) * qJD(2) + t305) * t298) * qJD(1); (-t250 + t389) * MDP(20) - t272 * MDP(21) + (t297 * t342 - t265 - t378) * MDP(22) + (t204 * t270 + t241 * t252 + t307) * MDP(24) + (t203 * t270 + t310) * MDP(25) + t366 * MDP(27) * t202 + (-t229 * MDP(26) + (t217 * t252 + t196) * MDP(27)) * pkin(5) + (t252 * MDP(19) - t241 * MDP(25) + MDP(26) * t366 - t390) * t251 + (MDP(21) * t337 + (t343 + (t297 * MDP(21) + t299 * MDP(22)) * qJD(5)) * t300) * qJD(1); (-t250 - t389) * MDP(26) + (-t200 * t252 - t202 * t251 + t208) * MDP(27);];
tauc  = t1;
