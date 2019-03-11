% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:31
% EndTime: 2019-03-08 18:54:36
% DurationCPUTime: 2.73s
% Computational Cost: add. (2466->312), mult. (6543->462), div. (0->0), fcn. (5476->12), ass. (0->154)
t309 = cos(pkin(6));
t297 = qJD(1) * t309 + qJD(2);
t304 = sin(pkin(12));
t306 = sin(pkin(6));
t312 = sin(qJ(3));
t315 = cos(qJ(3));
t307 = cos(pkin(12));
t308 = cos(pkin(7));
t384 = t307 * t308;
t320 = (t304 * t315 + t312 * t384) * t306;
t305 = sin(pkin(7));
t386 = t305 * t312;
t253 = qJD(1) * t320 + t297 * t386;
t311 = sin(qJ(4));
t314 = cos(qJ(4));
t333 = pkin(4) * t311 - pkin(10) * t314;
t290 = t333 * qJD(4);
t410 = -t253 + t290;
t385 = t305 * t315;
t405 = (-t304 * t312 + t315 * t384) * t306;
t409 = t309 * t385 + t405;
t368 = MDP(12) * t314;
t408 = qJD(1) * t405 + t297 * t385;
t407 = MDP(6) * t311;
t302 = t311 ^ 2;
t406 = MDP(7) * (-t314 ^ 2 + t302);
t367 = qJD(1) * t306;
t346 = t307 * t367;
t334 = t308 * t346;
t347 = t304 * t367;
t252 = -t312 * t347 + t315 * (t297 * t305 + t334);
t292 = -pkin(4) * t314 - pkin(10) * t311 - pkin(3);
t310 = sin(qJ(5));
t313 = cos(qJ(5));
t357 = qJD(5) * t313;
t380 = t313 * t314;
t404 = -t252 * t380 + t292 * t357 + t310 * t410;
t251 = qJD(3) * pkin(9) + t253;
t268 = t297 * t308 - t305 * t346;
t403 = -t251 * t311 + t268 * t314;
t364 = qJD(3) * t312;
t345 = t305 * t364;
t362 = qJD(3) * t315;
t249 = t297 * t345 + t334 * t364 + t347 * t362;
t402 = qJD(3) * t253 - t249;
t360 = qJD(4) * t311;
t382 = t310 * t314;
t398 = pkin(9) * t310;
t401 = t252 * t382 + t313 * t410 + t360 * t398;
t361 = qJD(4) * t310;
t365 = qJD(3) * t311;
t286 = t313 * t365 + t361;
t400 = t286 ^ 2;
t399 = pkin(5) * t311;
t397 = -qJ(6) - pkin(10);
t396 = qJD(3) * pkin(3);
t395 = qJ(6) * t311;
t248 = t408 * qJD(3);
t359 = qJD(4) * t314;
t217 = t248 * t311 + t251 * t359 + t268 * t360;
t394 = t217 * t310;
t393 = t217 * t313;
t231 = -qJD(4) * pkin(4) - t403;
t392 = t231 * t310;
t358 = qJD(5) * t310;
t341 = t311 * t358;
t353 = qJD(3) * qJD(4);
t339 = t314 * t353;
t351 = qJD(4) * qJD(5);
t371 = (t339 + t351) * t313;
t266 = qJD(3) * t341 - t371;
t391 = t266 * t310;
t363 = qJD(3) * t314;
t298 = -qJD(5) + t363;
t389 = t286 * t298;
t388 = t298 * t310;
t387 = t298 * t313;
t243 = qJD(3) * t292 - t252;
t383 = t310 * t243;
t381 = t311 * t313;
t235 = t251 * t314 + t268 * t311;
t232 = qJD(4) * pkin(10) + t235;
t213 = -t232 * t310 + t243 * t313;
t210 = -qJ(6) * t286 + t213;
t209 = -pkin(5) * t298 + t210;
t379 = t209 - t210;
t299 = pkin(9) * t380;
t330 = -qJ(6) * t380 + t399;
t356 = qJD(6) * t313;
t378 = -t311 * t356 + t330 * qJD(4) + (-t299 + (-t292 + t395) * t310) * qJD(5) + t401;
t289 = t333 * qJD(3);
t377 = t289 * t310 + t313 * t403;
t376 = (-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t381 + (-qJD(6) * t311 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t314) * t310 + t404;
t337 = qJD(5) * t397;
t375 = t356 - t377 + (qJ(6) * t363 + t337) * t310;
t335 = t289 * t313 - t310 * t403;
t374 = -qJD(3) * t330 - qJD(6) * t310 + t313 * t337 - t335;
t370 = t292 * t310 + t299;
t355 = t231 * qJD(5);
t354 = t313 * qJD(4);
t352 = qJD(4) * MDP(17);
t348 = pkin(5) * t310 + pkin(9);
t344 = t305 * t362;
t343 = t310 * t359;
t342 = t298 * t358;
t340 = t311 * t357;
t338 = t310 * t351;
t216 = qJD(4) * t403 + t248 * t314;
t239 = qJD(3) * t290 + t249;
t336 = t310 * t216 - t239 * t313;
t214 = t232 * t313 + t383;
t258 = t309 * t386 + t320;
t271 = -t305 * t306 * t307 + t308 * t309;
t242 = t258 * t314 + t271 * t311;
t223 = t242 * t313 - t310 * t409;
t222 = -t242 * t310 - t313 * t409;
t241 = t258 * t311 - t271 * t314;
t331 = qJD(3) * t302 - t298 * t314;
t267 = t338 + (t340 + t343) * qJD(3);
t212 = pkin(5) * t267 + t217;
t316 = qJD(4) ^ 2;
t329 = pkin(9) * t316 - t402;
t250 = -t252 - t396;
t328 = qJD(4) * (t250 + t252 - t396);
t273 = t308 * t311 + t314 * t386;
t261 = -t273 * t310 - t313 * t385;
t327 = -t273 * t313 + t310 * t385;
t272 = -t308 * t314 + t311 * t386;
t325 = -MDP(11) * t314 + MDP(12) * t311 - MDP(4);
t324 = -t216 * t313 + t232 * t358 - t239 * t310 - t243 * t357;
t318 = -qJD(5) * t214 - t336;
t317 = qJD(3) ^ 2;
t294 = t397 * t313;
t293 = t397 * t310;
t284 = t310 * t365 - t354;
t283 = t313 * t292;
t281 = t284 ^ 2;
t263 = -t310 * t395 + t370;
t260 = qJD(4) * t273 + t311 * t344;
t259 = -qJD(4) * t272 + t314 * t344;
t256 = -qJ(6) * t381 + t283 + (-pkin(5) - t398) * t314;
t255 = t258 * qJD(3);
t254 = t409 * qJD(3);
t230 = qJD(5) * t327 - t259 * t310 + t313 * t345;
t229 = qJD(5) * t261 + t259 * t313 + t310 * t345;
t224 = pkin(5) * t284 + qJD(6) + t231;
t221 = -qJD(4) * t241 + t254 * t314;
t220 = qJD(4) * t242 + t254 * t311;
t211 = -qJ(6) * t284 + t214;
t208 = qJD(5) * t222 + t221 * t313 + t255 * t310;
t207 = -qJD(5) * t223 - t221 * t310 + t255 * t313;
t206 = -qJ(6) * t267 - qJD(6) * t284 - t324;
t205 = qJ(6) * t266 - qJD(6) * t286 + t353 * t399 + t318;
t1 = [(-t207 * t298 + t220 * t284 + t241 * t267) * MDP(18) + (t208 * t298 + t220 * t286 - t241 * t266) * MDP(19) + (-t207 * t286 - t208 * t284 + t222 * t266 - t223 * t267) * MDP(20) + (t205 * t222 + t206 * t223 + t207 * t209 + t208 * t211 + t212 * t241 + t220 * t224) * MDP(21) + (-MDP(11) * t220 - MDP(12) * t221) * qJD(4) + (-t254 * MDP(5) + t325 * t255 + (-t409 * t368 + (-MDP(11) * t409 + MDP(18) * t222 - MDP(19) * t223) * t311) * qJD(4)) * qJD(3); (-t230 * t298 + t260 * t284 + t267 * t272) * MDP(18) + (t229 * t298 + t260 * t286 - t266 * t272) * MDP(19) + (-t229 * t284 - t230 * t286 + t261 * t266 + t267 * t327) * MDP(20) + (t205 * t261 - t206 * t327 + t209 * t230 + t211 * t229 + t212 * t272 + t224 * t260) * MDP(21) + (-t260 * MDP(11) - t259 * MDP(12) + (MDP(18) * t261 + MDP(19) * t327) * t365) * qJD(4) + ((-MDP(11) * t311 - t368) * t315 * t353 + (-MDP(5) * t315 + t312 * t325) * t317) * t305; t402 * MDP(4) + (t252 - t408) * qJD(3) * MDP(5) + 0.2e1 * t339 * t407 - 0.2e1 * t353 * t406 + (t311 * t328 - t314 * t329) * MDP(11) + (t311 * t329 + t314 * t328) * MDP(12) + (-t266 * t381 + (t314 * t354 - t341) * t286) * MDP(13) + ((-t284 * t313 - t286 * t310) * t359 + (t391 - t267 * t313 + (t284 * t310 - t286 * t313) * qJD(5)) * t311) * MDP(14) + (t298 * t341 + t266 * t314 + (t286 * t311 + t313 * t331) * qJD(4)) * MDP(15) + (t298 * t340 + t267 * t314 + (-t284 * t311 - t310 * t331) * qJD(4)) * MDP(16) + (-t298 - t363) * t311 * t352 + ((t292 * t358 - t401) * t298 + ((pkin(9) * t284 + t392) * qJD(4) + (t383 + (pkin(9) * t298 + t232) * t313) * qJD(5) + t336) * t314 + (t313 * t355 + pkin(9) * t267 + t394 - t252 * t284 + ((-pkin(9) * t382 + t283) * qJD(3) + t213) * qJD(4)) * t311) * MDP(18) + (t404 * t298 + (t231 * t354 + (qJD(4) * t286 - t342) * pkin(9) - t324) * t314 + (-t310 * t355 - pkin(9) * t266 + t393 - t252 * t286 + (-pkin(9) * t387 - qJD(3) * t370 - t214) * qJD(4)) * t311) * MDP(19) + (t256 * t266 - t263 * t267 - t378 * t286 - t376 * t284 + (-t209 * t313 - t211 * t310) * t359 + (-t205 * t313 - t206 * t310 + (t209 * t310 - t211 * t313) * qJD(5)) * t311) * MDP(20) + (t205 * t256 + t206 * t263 + t376 * t211 + t378 * t209 + t224 * t348 * t359 + (t212 * t348 + (pkin(5) * t357 - t252) * t224) * t311) * MDP(21) + (MDP(8) * t314 - MDP(9) * t311) * t316; (qJD(4) * t235 - t250 * t365 - t217) * MDP(11) + (-qJD(3) * t250 - t248) * t368 + (-t286 * t387 - t391) * MDP(13) + ((t284 * t298 - t266) * t313 + (-t267 + t389) * t310) * MDP(14) + (-t298 * t357 + (t298 * t380 + (-t286 + t361) * t311) * qJD(3)) * MDP(15) + (t342 + (-t298 * t382 + (t284 + t354) * t311) * qJD(3)) * MDP(16) + t298 * MDP(17) * t365 + (-pkin(4) * t267 - t393 + t335 * t298 - t235 * t284 + (pkin(10) * t387 + t392) * qJD(5) + (-t213 * t311 + (-pkin(10) * t360 - t231 * t314) * t310) * qJD(3)) * MDP(18) + (pkin(4) * t266 + t394 - t377 * t298 - t235 * t286 + (-pkin(10) * t388 + t231 * t313) * qJD(5) + (-t231 * t380 + (-pkin(10) * t354 + t214) * t311) * qJD(3)) * MDP(19) + (t266 * t293 + t267 * t294 - t374 * t286 - t375 * t284 + (t209 * t298 + t206) * t313 + (t211 * t298 - t205) * t310) * MDP(20) + (-t206 * t294 + t205 * t293 + t212 * (-pkin(5) * t313 - pkin(4)) + (-pkin(5) * t388 - t235) * t224 + t375 * t211 + t374 * t209) * MDP(21) + (-t314 * t407 + t406) * t317; (-t281 + t400) * MDP(14) + t371 * MDP(15) + (-t338 - t389) * MDP(16) + (-t214 * t298 - t231 * t286 + t318) * MDP(18) + (-t213 * t298 + t324) * MDP(19) + t379 * MDP(21) * t211 + (t266 * MDP(20) + (-t224 * t286 + t205) * MDP(21)) * pkin(5) + (MDP(13) * t286 - MDP(15) * t298 + MDP(19) * t231 - MDP(20) * t379) * t284 + (-MDP(16) * t343 + (t352 + (-MDP(15) * t310 - MDP(16) * t313) * qJD(5)) * t311) * qJD(3); (-t281 - t400) * MDP(20) + (t209 * t286 + t211 * t284 + t212) * MDP(21);];
tauc  = t1;
