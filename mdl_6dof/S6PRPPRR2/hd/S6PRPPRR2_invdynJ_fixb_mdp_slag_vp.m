% Calculate vector of inverse dynamics joint torques for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6PRPPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:08
% EndTime: 2019-03-08 19:20:12
% DurationCPUTime: 2.88s
% Computational Cost: add. (1395->340), mult. (3123->464), div. (0->0), fcn. (2624->12), ass. (0->155)
t312 = sin(qJ(2));
t315 = cos(qJ(2));
t306 = sin(pkin(6));
t375 = qJD(1) * t306;
t350 = qJD(2) * t375;
t363 = qJDD(1) * t306;
t411 = t312 * t363 + t315 * t350;
t304 = sin(pkin(11));
t307 = cos(pkin(11));
t382 = t315 * t307;
t391 = t306 * t312;
t410 = -t304 * t391 + t306 * t382;
t288 = t315 * t363;
t259 = qJDD(2) * pkin(2) - t312 * t350 + t288;
t230 = t259 * t307 - t411 * t304;
t334 = qJDD(4) - t230;
t405 = -pkin(3) - pkin(8);
t219 = t405 * qJDD(2) + t334;
t354 = t315 * t375;
t275 = qJD(2) * pkin(2) + t354;
t355 = t312 * t375;
t279 = t304 * t355;
t249 = t275 * t307 - t279;
t340 = qJD(4) - t249;
t243 = t405 * qJD(2) + t340;
t309 = cos(pkin(6));
t291 = qJD(1) * t309 + qJD(3);
t311 = sin(qJ(5));
t314 = cos(qJ(5));
t229 = t243 * t311 + t291 * t314;
t289 = t309 * qJDD(1) + qJDD(3);
t396 = t289 * t311;
t210 = -qJDD(5) * pkin(5) + qJD(5) * t229 - t219 * t314 + t396;
t293 = qJD(2) * t311 + qJD(6);
t305 = sin(pkin(10));
t308 = cos(pkin(10));
t267 = t304 * t312 - t382;
t330 = t267 * t309;
t337 = t304 * t315 + t307 * t312;
t235 = -t305 * t337 - t308 * t330;
t238 = t305 * t330 - t308 * t337;
t338 = -t309 * t311 - t314 * t410;
t392 = t306 * t311;
t329 = g(1) * (-t238 * t314 - t305 * t392) + g(2) * (-t235 * t314 + t308 * t392) + g(3) * t338;
t342 = pkin(5) * t314 + pkin(9) * t311;
t409 = (pkin(9) * qJD(6) + t342 * qJD(2)) * t293 + t210 + t329;
t258 = t307 * t354 - t279;
t336 = pkin(5) * t311 - pkin(9) * t314 + qJ(4);
t404 = pkin(2) * t304;
t263 = t336 + t404;
t364 = qJD(2) * qJD(5);
t348 = t314 * t364;
t362 = qJDD(2) * t311;
t265 = qJDD(6) + t348 + t362;
t266 = t342 * qJD(5) + qJD(4);
t408 = (t258 - t266) * t293 - t263 * t265;
t386 = t309 * t315;
t387 = t309 * t312;
t378 = -t304 * t386 - t307 * t387;
t239 = -t308 * t267 + t305 * t378;
t234 = t305 * t267 + t308 * t378;
t250 = t275 * t304 + t307 * t355;
t245 = qJD(2) * qJ(4) + t250;
t261 = t337 * t306;
t255 = qJD(1) * t261;
t296 = -pkin(2) * t307 - pkin(3);
t292 = -pkin(8) + t296;
t294 = qJ(4) + t404;
t407 = qJDD(5) * t292 + (qJD(2) * t294 + t245 - t255) * qJD(5);
t227 = qJD(5) * pkin(9) + t229;
t328 = -g(1) * t239 + g(2) * t234 - g(3) * t261;
t406 = -(t292 * t293 + t227) * qJD(6) + t328;
t403 = t289 * t309 - g(3);
t402 = qJDD(2) * pkin(3);
t313 = cos(qJ(6));
t366 = t313 * qJD(5);
t351 = t311 * t366;
t310 = sin(qJ(6));
t368 = qJD(6) * t310;
t326 = -t314 * t368 - t351;
t361 = qJDD(2) * t314;
t356 = qJD(6) * t366 + t310 * qJDD(5) + t313 * t361;
t241 = t326 * qJD(2) + t356;
t401 = t241 * t310;
t371 = qJD(5) * t310;
t374 = qJD(2) * t314;
t271 = t313 * t374 + t371;
t349 = t311 * t364;
t335 = -t313 * qJDD(5) + (-t349 + t361) * t310;
t242 = t271 * qJD(6) + t335;
t400 = t242 * t311;
t353 = t310 * t374;
t269 = t353 - t366;
t398 = t269 * t293;
t397 = t271 * t293;
t395 = t289 * t314;
t393 = t305 * t312;
t390 = t306 * t314;
t389 = t306 * t315;
t385 = t310 * t293;
t384 = t313 * t293;
t383 = t313 * t314;
t381 = qJDD(1) - g(3);
t369 = qJD(5) * t314;
t380 = t241 * t311 + t271 * t369;
t303 = t314 ^ 2;
t377 = t311 ^ 2 - t303;
t316 = qJD(5) ^ 2;
t317 = qJD(2) ^ 2;
t376 = -t316 - t317;
t373 = qJD(5) * t269;
t372 = qJD(5) * t292;
t370 = qJD(5) * t311;
t367 = qJD(6) * t313;
t365 = qJD(4) - t258;
t357 = t308 * t386;
t231 = t304 * t259 + t411 * t307;
t352 = t293 * t371;
t344 = qJD(2) * t245 - t219;
t233 = t336 * qJD(2) + t250;
t212 = t227 * t313 + t233 * t310;
t339 = t227 * t310 - t233 * t313;
t228 = t243 * t314 - t291 * t311;
t247 = t309 * t314 - t311 * t410;
t218 = t247 * t313 + t261 * t310;
t217 = -t247 * t310 + t261 * t313;
t333 = -t305 * t386 - t308 * t312;
t332 = t310 * t265 + t293 * t367;
t331 = -t313 * t265 + t293 * t368;
t327 = g(1) * t238 + g(2) * t235 + g(3) * t410;
t324 = -qJD(6) * t263 * t293 - t327;
t323 = (-t332 + t373) * MDP(21);
t226 = -qJD(5) * pkin(5) - t228;
t322 = -pkin(9) * t265 + (t226 + t228) * t293;
t321 = -g(1) * t333 - g(3) * t389;
t209 = qJDD(5) * pkin(9) + qJD(5) * t228 + t219 * t311 + t395;
t320 = -qJD(5) * t226 - qJD(6) * t233 + t255 * t293 - t265 * t292 - t209;
t319 = t327 + t334;
t220 = qJDD(2) * qJ(4) + qJD(4) * qJD(2) + t231;
t318 = t365 * qJD(2) + qJDD(2) * t294 - t292 * t316 + t220 + t328;
t283 = pkin(2) * t357;
t282 = qJDD(5) * t314 - t311 * t316;
t281 = -qJDD(5) * t311 - t314 * t316;
t264 = t311 * t352;
t257 = t410 * qJD(2);
t256 = qJD(2) * t261;
t244 = -qJD(2) * pkin(3) + t340;
t225 = t334 - t402;
t224 = t235 * t311 + t308 * t390;
t222 = -t238 * t311 + t305 * t390;
t216 = t247 * qJD(5) - t256 * t314;
t215 = t338 * qJD(5) + t256 * t311;
t214 = t266 * qJD(2) + t336 * qJDD(2) + t231;
t213 = t313 * t214;
t1 = [t381 * MDP(1) + (t230 * t410 + t231 * t261 - t249 * t256 + t250 * t257 + t403) * MDP(5) + (t220 * t261 - t225 * t410 + t244 * t256 + t245 * t257 + t403) * MDP(8) + (-qJD(5) * t216 + qJDD(5) * t338) * MDP(14) + (-qJD(5) * t215 - qJDD(5) * t247) * MDP(15) + ((-t218 * qJD(6) - t215 * t310 + t257 * t313) * t293 + t217 * t265 + t216 * t269 - t338 * t242) * MDP(21) + (-(t217 * qJD(6) + t215 * t313 + t257 * t310) * t293 - t218 * t265 + t216 * t271 - t338 * t241) * MDP(22) + (-t410 * MDP(6) + (MDP(14) * t311 + MDP(15) * t314 + MDP(7)) * t261) * qJDD(2) + (t256 * MDP(6) + t257 * MDP(7) + (t257 * t311 + t261 * t369) * MDP(14) + (t257 * t314 - t261 * t370) * MDP(15)) * qJD(2) + ((qJDD(2) * t315 - t312 * t317) * MDP(3) + (-qJDD(2) * t312 - t315 * t317) * MDP(4)) * t306; qJDD(2) * MDP(2) + (t288 - g(2) * (t357 - t393) + t321) * MDP(3) + (-g(1) * (t305 * t387 - t308 * t315) - g(2) * (-t305 * t315 - t308 * t387) - t381 * t391) * MDP(4) + (-g(2) * t283 + t249 * t255 - t250 * t258 + (g(2) * t393 + t230 * t307 + t231 * t304 + t321) * pkin(2)) * MDP(5) + (-qJD(2) * t255 + (-pkin(3) + t296) * qJDD(2) + t319) * MDP(6) + ((qJ(4) + t294) * qJDD(2) + (0.2e1 * qJD(4) - t258) * qJD(2) + t328 + t231) * MDP(7) + (t220 * t294 + t225 * t296 - t244 * t255 - g(1) * (t333 * pkin(2) + pkin(3) * t238 + qJ(4) * t239) - g(2) * (-pkin(2) * t393 + pkin(3) * t235 - qJ(4) * t234 + t283) - g(3) * (pkin(2) * t389 + pkin(3) * t410 + qJ(4) * t261) + t365 * t245) * MDP(8) + (qJDD(2) * t303 - 0.2e1 * t311 * t348) * MDP(9) + 0.2e1 * (-t311 * t361 + t377 * t364) * MDP(10) + t282 * MDP(11) + t281 * MDP(12) + (t318 * t311 + t407 * t314) * MDP(14) + (-t407 * t311 + t318 * t314) * MDP(15) + (t241 * t383 + t326 * t271) * MDP(16) + ((t269 * t313 + t271 * t310) * t370 + (-t401 - t242 * t313 + (t269 * t310 - t271 * t313) * qJD(6)) * t314) * MDP(17) + (t265 * t383 + t326 * t293 + t380) * MDP(18) + (-t400 + t264 + (-t332 - t373) * t314) * MDP(19) + (t265 * t311 + t293 * t369) * MDP(20) + (-t408 * t313 + t324 * t310 + (t269 * t372 + t320 * t310 + t406 * t313 + t213) * t311 + (t226 * t367 + t210 * t310 - t292 * t242 + t255 * t269 + (-t292 * t385 - t339) * qJD(5)) * t314) * MDP(21) + (t408 * t310 + t324 * t313 + (t271 * t372 + t320 * t313 + (-t214 - t406) * t310) * t311 + (-t226 * t368 + t210 * t313 - t292 * t241 + t255 * t271 + (-t292 * t384 - t212) * qJD(5)) * t314) * MDP(22); t281 * MDP(14) - t282 * MDP(15) + (t264 + t400) * MDP(21) + (t293 * t351 + t380) * MDP(22) + (MDP(5) + MDP(8)) * (-g(3) * t309 + (-g(1) * t305 + g(2) * t308) * t306 + t289) + (t331 * MDP(22) + t323) * t314; qJDD(2) * MDP(6) - t317 * MDP(7) + (t319 - t402) * MDP(8) + (-t245 * MDP(8) + (-MDP(21) * t313 + MDP(22) * t310) * t293) * qJD(2) + (qJDD(5) * MDP(14) + t376 * MDP(15) + (-t242 - t352) * MDP(21) + (-t293 * t366 - t241) * MDP(22)) * t314 + (t376 * MDP(14) - qJDD(5) * MDP(15) + t323 + (qJD(5) * t271 + t331) * MDP(22)) * t311; MDP(11) * t361 - MDP(12) * t362 + qJDD(5) * MDP(13) + (-t344 * t314 - t329 - t396) * MDP(14) + (g(1) * t222 - g(2) * t224 + g(3) * t247 + t344 * t311 - t395) * MDP(15) + (t271 * t384 + t401) * MDP(16) + ((t241 - t398) * t313 + (-t242 - t397) * t310) * MDP(17) + ((-t271 * t314 + t311 * t384) * qJD(2) + t332) * MDP(18) + ((t269 * t314 - t311 * t385) * qJD(2) - t331) * MDP(19) - t293 * MDP(20) * t374 + (-pkin(5) * t242 - t229 * t269 + t322 * t310 - t409 * t313 + t339 * t374) * MDP(21) + (-pkin(5) * t241 + t212 * t374 - t229 * t271 + t409 * t310 + t322 * t313) * MDP(22) + (t314 * t311 * MDP(9) - t377 * MDP(10)) * t317; t271 * t269 * MDP(16) + (-t269 ^ 2 + t271 ^ 2) * MDP(17) + (-t313 * t349 + t356 + t398) * MDP(18) + (-t335 + t397) * MDP(19) + t265 * MDP(20) + (-t310 * t209 + t213 + t212 * t293 - t226 * t271 - g(1) * (-t222 * t310 + t239 * t313) - g(2) * (t224 * t310 - t234 * t313) - g(3) * t217) * MDP(21) + (-t313 * t209 - t310 * t214 - t339 * t293 + t226 * t269 - g(1) * (-t222 * t313 - t239 * t310) - g(2) * (t224 * t313 + t234 * t310) + g(3) * t218) * MDP(22) + (-MDP(18) * t353 - t271 * MDP(19) - t212 * MDP(21) + t339 * MDP(22)) * qJD(6);];
tau  = t1;
