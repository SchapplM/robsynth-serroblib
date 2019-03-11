% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:51:07
% EndTime: 2019-03-08 18:51:11
% DurationCPUTime: 2.62s
% Computational Cost: add. (1616->301), mult. (4348->444), div. (0->0), fcn. (3601->12), ass. (0->155)
t291 = sin(pkin(12));
t293 = sin(pkin(6));
t299 = sin(qJ(3));
t302 = cos(qJ(3));
t294 = cos(pkin(12));
t295 = cos(pkin(7));
t381 = t294 * t295;
t311 = (t291 * t302 + t299 * t381) * t293;
t296 = cos(pkin(6));
t283 = qJD(1) * t296 + qJD(2);
t292 = sin(pkin(7));
t383 = t292 * t299;
t349 = t283 * t383;
t241 = qJD(1) * t311 + t349;
t239 = qJD(3) * pkin(9) + t241;
t375 = qJD(1) * t293;
t342 = t294 * t375;
t256 = t283 * t295 - t292 * t342;
t298 = sin(qJ(4));
t301 = cos(qJ(4));
t378 = -t298 * t239 + t301 * t256;
t401 = -qJD(5) + t378;
t219 = -qJD(4) * pkin(4) - t401;
t382 = t292 * t302;
t398 = (-t291 * t299 + t302 * t381) * t293;
t402 = t296 * t382 + t398;
t400 = qJD(1) * t398 + t283 * t382;
t289 = t298 ^ 2;
t290 = t301 ^ 2;
t399 = (t289 - t290) * MDP(7);
t328 = t295 * t342;
t343 = t291 * t375;
t369 = qJD(3) * t302;
t371 = qJD(3) * t299;
t237 = qJD(3) * t349 + t328 * t371 + t343 * t369;
t397 = qJD(3) * t241 - t237;
t367 = qJD(4) * t298;
t396 = -pkin(4) * t367 + t241;
t222 = t301 * t239 + t298 * t256;
t220 = -qJD(4) * qJ(5) - t222;
t372 = qJD(3) * t298;
t350 = pkin(5) * t372;
t355 = t350 - t401;
t240 = -t299 * t343 + (t283 * t292 + t328) * t302;
t365 = qJD(5) * t298;
t366 = qJD(4) * t301;
t316 = -qJ(5) * t366 - t365;
t353 = qJD(3) * qJD(4);
t337 = t298 * t353;
t329 = pkin(4) * t337 + t237;
t224 = t316 * qJD(3) + t329;
t379 = -t316 + t396;
t304 = qJD(4) ^ 2;
t393 = pkin(9) * t304;
t395 = t379 * qJD(3) - t224 - t393;
t303 = -pkin(4) - pkin(10);
t394 = pkin(5) + pkin(9);
t392 = qJD(3) * pkin(3);
t236 = t400 * qJD(3);
t344 = -t301 * t236 + t239 * t367 - t256 * t366;
t208 = (qJD(5) - t350) * qJD(4) - t344;
t297 = sin(qJ(6));
t391 = t208 * t297;
t300 = cos(qJ(6));
t390 = t208 * t300;
t357 = t300 * qJD(3);
t359 = t297 * qJD(4);
t265 = t301 * t357 + t359;
t277 = t297 * t337;
t253 = -qJD(6) * t265 + t277;
t389 = t253 * t300;
t285 = qJD(6) + t372;
t388 = t265 * t285;
t360 = t297 * qJD(3);
t338 = t301 * t360;
t356 = t300 * qJD(4);
t267 = -t338 + t356;
t387 = t267 * t285;
t386 = t267 * t301;
t385 = t285 * t298;
t384 = t285 * t303;
t327 = pkin(10) * t298 - qJ(5) * t301;
t309 = t327 * qJD(4) - t365;
t380 = -t309 + t396;
t335 = t300 * t353;
t377 = qJD(6) * t338 + t298 * t335;
t333 = -qJ(5) * t298 - pkin(3);
t273 = -pkin(4) * t301 + t333;
t373 = qJD(3) * t273;
t370 = qJD(3) * t301;
t364 = qJD(6) * t297;
t363 = qJD(6) * t300;
t362 = qJD(6) * t302;
t286 = pkin(5) * t370;
t216 = -t220 + t286;
t361 = t216 * qJD(6);
t358 = t300 * MDP(23);
t352 = MDP(11) - MDP(14);
t351 = MDP(12) - MDP(15);
t276 = t394 * t301;
t346 = t298 * t383;
t305 = qJD(3) ^ 2;
t345 = t298 * t301 * t305;
t212 = t298 * t236 + t239 * t366 + t256 * t367;
t341 = t292 * t369;
t340 = t285 * t364;
t339 = t301 * t363;
t336 = t301 * t353;
t334 = MDP(21) * t370;
t210 = pkin(5) * t336 + t212;
t223 = t309 * qJD(3) + t329;
t330 = t300 * t210 - t223 * t297;
t215 = t303 * qJD(4) + t355;
t264 = t303 * t301 + t333;
t225 = t264 * qJD(3) - t240;
t206 = t215 * t300 - t297 * t225;
t207 = t215 * t297 + t225 * t300;
t246 = t296 * t383 + t311;
t259 = -t292 * t293 * t294 + t295 * t296;
t228 = t246 * t301 + t259 * t298;
t275 = t394 * t298;
t326 = t264 * t300 + t275 * t297;
t325 = -qJD(3) * t290 + t385;
t324 = t285 * t297;
t322 = qJD(4) * t222 - t212;
t321 = t393 - t397;
t238 = -t240 - t392;
t320 = qJD(4) * (t238 + t240 - t392);
t261 = t295 * t298 + t301 * t383;
t318 = t216 * t298 + t303 * t366;
t229 = -t240 + t373;
t317 = qJD(4) * (-t229 - t240 - t373);
t315 = t300 * MDP(22) - t297 * MDP(23) + MDP(13);
t310 = t297 * MDP(22) - t351 + t358;
t308 = t351 * t298 - t352 * t301 - MDP(4);
t211 = -qJD(4) * qJD(5) + t344;
t306 = -t211 * t301 + t212 * t298 + (t219 * t301 + t220 * t298) * qJD(4);
t288 = pkin(4) * t372;
t279 = t301 * t335;
t271 = qJD(4) * t276;
t270 = t394 * t367;
t268 = -qJ(5) * t370 + t288;
t260 = -t295 * t301 + t346;
t257 = t327 * qJD(3) + t288;
t254 = qJD(6) * t356 - t377;
t248 = qJD(4) * t261 + t298 * t341;
t247 = qJD(4) * t346 - t295 * t366 - t301 * t341;
t244 = t246 * qJD(3);
t243 = t402 * qJD(3);
t227 = t246 * t298 - t259 * t301;
t226 = t229 * t372;
t218 = t286 + t222;
t214 = -t246 * t367 + (qJD(4) * t259 + t243) * t301;
t213 = qJD(4) * t228 + t243 * t298;
t1 = [(-t211 * t228 + t212 * t227 + t213 * t219 - t214 * t220 - t224 * t402 + t229 * t244) * MDP(16) + ((t213 * t300 - t244 * t297 + (-t227 * t297 + t300 * t402) * qJD(6)) * t285 + t214 * t265 + t228 * t254) * MDP(22) + (-(t213 * t297 + t244 * t300 + (t227 * t300 + t297 * t402) * qJD(6)) * t285 + t214 * t267 + t228 * t253) * MDP(23) + (-t352 * t213 - t351 * t214) * qJD(4) + (-t243 * MDP(5) + (t213 * t298 + t214 * t301) * MDP(13) + t308 * t244 + ((-t228 * MDP(13) - t352 * t402) * t298 + (t315 * t227 + t310 * t402) * t301) * qJD(4)) * qJD(3); (-t211 * t261 + t212 * t260 + t219 * t248 + t220 * t247) * MDP(16) + ((t248 * t300 - t260 * t364) * t285 - t247 * t265 + t261 * t254) * MDP(22) + (-(t248 * t297 + t260 * t363) * t285 - t247 * t267 + t261 * t253) * MDP(23) + (-t247 * t301 + t248 * t298) * MDP(13) * qJD(3) + ((-t224 * t302 + t229 * t371) * MDP(16) + ((-t299 * t360 + t300 * t362) * MDP(22) - (t297 * t362 + t299 * t357) * MDP(23)) * t285 + (-t302 * MDP(5) + t308 * t299) * t305) * t292 + (-t352 * t248 + t351 * t247 + ((-t261 * MDP(13) - t352 * t382) * t298 + (t315 * t260 + t310 * t382) * t301) * qJD(3)) * qJD(4); t397 * MDP(4) + (t240 - t400) * qJD(3) * MDP(5) + 0.2e1 * t298 * MDP(6) * t336 - 0.2e1 * t353 * t399 + (t298 * t320 - t321 * t301) * MDP(11) + (t321 * t298 + t301 * t320) * MDP(12) + ((-t289 - t290) * t240 * qJD(3) + t306) * MDP(13) + (t298 * t317 - t301 * t395) * MDP(14) + (t298 * t395 + t301 * t317) * MDP(15) + (t224 * t273 + (-t219 * t298 + t220 * t301) * t240 - t379 * t229 + t306 * pkin(9)) * MDP(16) + (-t253 * t297 * t301 + (t298 * t359 - t339) * t267) * MDP(17) + ((-t265 * t297 + t267 * t300) * t367 + (-t389 + t254 * t297 + (t265 * t300 + t267 * t297) * qJD(6)) * t301) * MDP(18) + (-t285 * t339 + t253 * t298 + (t325 * t297 + t386) * qJD(4)) * MDP(19) + (t301 * t340 - t254 * t298 + (-t265 * t301 + t325 * t300) * qJD(4)) * MDP(20) + (t285 + t372) * MDP(21) * t366 + (t276 * t254 - t270 * t265 + (-t216 * t356 + t330) * t298 + ((-t240 * t298 + t271) * t300 + t380 * t297) * t285 + (-t207 * t298 - t326 * t285) * qJD(6) + (-t297 * t361 + t390 - t240 * t265 + ((-t264 * t297 + t275 * t300) * qJD(3) + t206) * qJD(4)) * t301) * MDP(22) + (t276 * t253 - t270 * t267 + (-(qJD(6) * t215 + t223) * t298 + (-qJD(6) * t275 + t380) * t285) * t300 + (-(-qJD(6) * t264 + t271) * t285 + (qJD(4) * t216 + qJD(6) * t225 + t240 * t285 - t210) * t298) * t297 + (-t300 * t361 - t391 - t240 * t267 + (-t326 * qJD(3) - t207) * qJD(4)) * t301) * MDP(23) + (t301 * MDP(8) - t298 * MDP(9)) * t304; -MDP(6) * t345 + t305 * t399 + (-t238 * t372 + t322) * MDP(11) + (qJD(4) * t378 - t238 * t370 + t344) * MDP(12) + (-t268 * t370 + t226 - t322) * MDP(14) + ((0.2e1 * qJD(5) - t378) * qJD(4) + (t229 * t301 + t268 * t298) * qJD(3) - t344) * MDP(15) + (-pkin(4) * t212 - qJ(5) * t211 - t219 * t222 + t220 * t401 - t229 * t268) * MDP(16) + (-t267 * t324 + t389) * MDP(17) + ((-t254 - t387) * t300 + (-t253 + t388) * t297) * MDP(18) + (-t340 + t279 + (-t297 * t385 - t386) * qJD(3)) * MDP(19) + (-t285 * t363 + (-t300 * t385 + (t265 - t359) * t301) * qJD(3)) * MDP(20) - t285 * t334 + (qJ(5) * t254 + t391 - (t218 * t300 - t257 * t297) * t285 + t355 * t265 + (t216 * t300 - t297 * t384) * qJD(6) + (-t206 * t301 + t318 * t300) * qJD(3)) * MDP(22) + (qJ(5) * t253 + t390 + (t218 * t297 + t257 * t300) * t285 + t355 * t267 + (-t216 * t297 - t300 * t384) * qJD(6) + (t207 * t301 - t318 * t297) * qJD(3)) * MDP(23); MDP(14) * t345 + (-t289 * t305 - t304) * MDP(15) + (t226 + t212) * MDP(16) + t279 * MDP(22) + (t220 * MDP(16) - t265 * MDP(22) + (-t267 - t338) * MDP(23)) * qJD(4) + (-MDP(22) * t324 - t285 * t358) * t285; t267 * t265 * MDP(17) + (-t265 ^ 2 + t267 ^ 2) * MDP(18) + (t277 + t388) * MDP(19) + (t377 + t387) * MDP(20) + qJD(4) * t334 + (t207 * t285 - t216 * t267 + t330) * MDP(22) + (t206 * t285 - t210 * t297 + t216 * t265 - t223 * t300) * MDP(23) + (-t265 * MDP(19) - MDP(20) * t356 - t207 * MDP(22) - t206 * MDP(23)) * qJD(6);];
tauc  = t1;
