% Calculate vector of inverse dynamics joint torques for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:30
% EndTime: 2019-12-05 16:23:36
% DurationCPUTime: 2.57s
% Computational Cost: add. (1479->276), mult. (3348->375), div. (0->0), fcn. (2524->12), ass. (0->139)
t386 = qJ(4) + pkin(6);
t312 = sin(pkin(9));
t314 = cos(pkin(9));
t318 = sin(qJ(3));
t321 = cos(qJ(3));
t343 = t312 * t318 - t314 * t321;
t280 = t343 * qJD(2);
t320 = cos(qJ(5));
t270 = t320 * t280;
t289 = t312 * t321 + t314 * t318;
t282 = t289 * qJD(2);
t317 = sin(qJ(5));
t380 = t282 * t317;
t237 = -t270 - t380;
t309 = qJD(3) + qJD(5);
t381 = t237 * t309;
t322 = cos(qJ(2));
t313 = sin(pkin(8));
t315 = cos(pkin(8));
t353 = g(1) * t315 + g(2) * t313;
t335 = t353 * t322;
t319 = sin(qJ(2));
t388 = g(3) * t319;
t329 = t335 + t388;
t387 = g(3) * t322;
t398 = t353 * t319;
t403 = t398 - t387;
t374 = qJDD(1) - g(3);
t402 = t374 * t322 + t398;
t366 = qJD(1) * qJD(2);
t287 = qJDD(2) * pkin(6) + qJDD(1) * t319 + t322 * t366;
t333 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + t287;
t369 = qJD(1) * t319;
t355 = t386 * qJD(2) + t369;
t342 = qJD(3) * t355;
t222 = qJDD(3) * pkin(3) - t333 * t318 - t321 * t342;
t223 = -t318 * t342 + t333 * t321;
t207 = t314 * t222 - t223 * t312;
t365 = qJD(2) * qJD(3);
t358 = t321 * t365;
t359 = t318 * t365;
t243 = t289 * qJDD(2) - t312 * t359 + t314 * t358;
t203 = qJDD(3) * pkin(4) - pkin(7) * t243 + t207;
t208 = t312 * t222 + t314 * t223;
t281 = t289 * qJD(3);
t242 = -qJD(2) * t281 - t343 * qJDD(2);
t204 = pkin(7) * t242 + t208;
t275 = t355 * t318;
t384 = qJD(3) * pkin(3);
t262 = -t275 + t384;
t276 = t355 * t321;
t378 = t314 * t276;
t225 = t312 * t262 + t378;
t392 = pkin(7) * t280;
t214 = t225 - t392;
t306 = pkin(3) * t321 + pkin(2);
t368 = qJD(1) * t322;
t285 = -t306 * qJD(2) + qJD(4) - t368;
t246 = pkin(4) * t280 + t285;
t307 = qJ(3) + pkin(9) + qJ(5);
t301 = sin(t307);
t302 = cos(t307);
t367 = qJD(5) * t317;
t377 = t315 * t322;
t379 = t313 * t322;
t401 = -t246 * t237 - g(1) * (-t301 * t313 - t302 * t377) - g(2) * (t301 * t315 - t302 * t379) - t317 * t203 - t320 * t204 + t214 * t367 + t302 * t388;
t308 = qJDD(3) + qJDD(5);
t346 = -t280 * t317 + t320 * t282;
t400 = t308 * MDP(18) + (-t237 ^ 2 + t346 ^ 2) * MDP(15) - t237 * MDP(14) * t346;
t382 = t346 * t309;
t357 = qJD(3) * t386;
t277 = qJD(4) * t321 - t318 * t357;
t278 = -qJD(4) * t318 - t321 * t357;
t373 = -t277 * t312 + t314 * t278 + t289 * t368;
t372 = t314 * t277 + t312 * t278 + t343 * t368;
t396 = t318 * t384 - t369;
t284 = t343 * qJD(3);
t395 = -t246 * t346 - g(1) * (-t301 * t377 + t302 * t313) - g(2) * (-t301 * t379 - t302 * t315) + t320 * t203 - t317 * t204 + t301 * t388;
t304 = t319 * t366;
t323 = qJD(3) ^ 2;
t364 = qJDD(1) * t322;
t394 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t323 + (t353 + t366) * t319 - t304 + t364 - t387;
t356 = -t320 * t242 + t243 * t317;
t206 = t346 * qJD(5) + t356;
t393 = pkin(3) * t312;
t391 = pkin(7) * t282;
t385 = qJD(2) * pkin(2);
t258 = t312 * t276;
t224 = t314 * t262 - t258;
t212 = qJD(3) * pkin(4) + t224 - t391;
t376 = t320 * t212;
t324 = qJD(2) ^ 2;
t375 = t321 * t324;
t227 = -t314 * t275 - t258;
t295 = t386 * t318;
t296 = t386 * t321;
t249 = -t312 * t295 + t314 * t296;
t310 = t318 ^ 2;
t371 = -t321 ^ 2 + t310;
t363 = qJDD(2) * t318;
t362 = qJDD(2) * t321;
t360 = -qJD(5) * t270 + t317 * t242 + t320 * t243;
t226 = t275 * t312 - t378;
t248 = -t314 * t295 - t296 * t312;
t354 = pkin(4) * t281 + t396;
t352 = g(1) * t313 - g(2) * t315;
t233 = -pkin(7) * t343 + t249;
t351 = -pkin(7) * t284 + qJD(5) * t233 - t373;
t232 = -pkin(7) * t289 + t248;
t350 = pkin(7) * t281 - qJD(5) * t232 - t372;
t349 = -t317 * t212 - t320 * t214;
t273 = t289 * t319;
t274 = t343 * t319;
t348 = -t273 * t320 + t274 * t317;
t347 = -t273 * t317 - t274 * t320;
t345 = -t289 * t317 - t320 * t343;
t245 = t289 * t320 - t317 * t343;
t341 = qJDD(3) * t321 - t318 * t323;
t340 = qJDD(3) * t318 + t321 * t323;
t303 = pkin(3) * t314 + pkin(4);
t339 = t303 * t317 + t320 * t393;
t338 = t303 * t320 - t317 * t393;
t334 = t352 * t321;
t205 = -t282 * t367 + t360;
t330 = pkin(3) * t359 - t306 * qJDD(2) + qJDD(4) + t304;
t299 = -t368 - t385;
t327 = -pkin(6) * qJDD(3) + (t299 + t368 - t385) * qJD(3);
t247 = t330 - t364;
t325 = -t299 * qJD(2) - t287 + t329;
t266 = pkin(4) * t343 - t306;
t252 = pkin(3) * qJD(2) * t318 + pkin(4) * t282;
t229 = -t322 * t280 - t319 * t281;
t228 = -t322 * t282 + t319 * t284;
t216 = t227 - t391;
t215 = t226 + t392;
t213 = -pkin(4) * t242 + t247;
t210 = t245 * qJD(5) + t320 * t281 - t284 * t317;
t209 = t345 * qJD(5) - t281 * t317 - t284 * t320;
t1 = [t374 * MDP(1) + (-t228 * t282 - t229 * t280 - t242 * t274 + t243 * t273) * MDP(12) + (-t207 * t273 - t208 * t274 + t224 * t228 + t225 * t229 - g(3)) * MDP(13) + ((-t347 * qJD(5) + t228 * t320 - t229 * t317) * t309 + t348 * t308) * MDP(19) + (-(t348 * qJD(5) + t228 * t317 + t229 * t320) * t309 - t347 * t308) * MDP(20) + (qJDD(2) * MDP(3) - t324 * MDP(4) + (-0.2e1 * t359 + t362) * MDP(10) + (-0.2e1 * t358 - t363) * MDP(11) - t247 * MDP(13) - t206 * MDP(19) - t205 * MDP(20)) * t322 + (-t324 * MDP(3) - qJDD(2) * MDP(4) + (-t340 - t375) * MDP(10) + (t318 * t324 - t341) * MDP(11) + (t285 * MDP(13) - MDP(19) * t237 + MDP(20) * t346) * qJD(2)) * t319; qJDD(2) * MDP(2) + t402 * MDP(3) + (-t374 * t319 + t335) * MDP(4) + (qJDD(2) * t310 + 0.2e1 * t318 * t358) * MDP(5) + 0.2e1 * (t318 * t362 - t371 * t365) * MDP(6) + t340 * MDP(7) + t341 * MDP(8) + (t327 * t318 + t394 * t321) * MDP(10) + (-t394 * t318 + t327 * t321) * MDP(11) + (-t207 * t289 - t208 * t343 + t224 * t284 - t225 * t281 + t242 * t249 - t243 * t248 - t372 * t280 - t373 * t282 - t329) * MDP(12) + (t208 * t249 + t207 * t248 - t247 * t306 - g(3) * (t306 * t322 + t319 * t386) + t396 * t285 + t372 * t225 + t373 * t224 + t353 * (t306 * t319 - t322 * t386)) * MDP(13) + (t205 * t245 + t209 * t346) * MDP(14) + (t205 * t345 - t206 * t245 + t209 * t237 - t210 * t346) * MDP(15) + (t209 * t309 + t245 * t308) * MDP(16) + (-t210 * t309 + t308 * t345) * MDP(17) + ((t232 * t320 - t233 * t317) * t308 + t266 * t206 - t213 * t345 + t246 * t210 + (t350 * t317 - t351 * t320) * t309 - t354 * t237 + t403 * t302) * MDP(19) + (-(t232 * t317 + t233 * t320) * t308 + t266 * t205 + t213 * t245 + t246 * t209 + (t351 * t317 + t350 * t320) * t309 + t354 * t346 - t403 * t301) * MDP(20); -t318 * MDP(5) * t375 + t371 * MDP(6) * t324 + MDP(7) * t363 + MDP(8) * t362 + qJDD(3) * MDP(9) + (t325 * t318 - t334) * MDP(10) + (t352 * t318 + t325 * t321) * MDP(11) + ((t225 + t226) * t282 - (t224 - t227) * t280 + (t242 * t312 - t243 * t314) * pkin(3)) * MDP(12) + (-t224 * t226 - t225 * t227 + (t207 * t314 + t208 * t312 - t334 + (-t285 * qJD(2) + t329) * t318) * pkin(3)) * MDP(13) + (t205 - t381) * MDP(16) + (-t206 + t382) * MDP(17) + (t338 * t308 - (t215 * t320 - t216 * t317) * t309 + t252 * t237 + (-t339 * t309 + t349) * qJD(5) + t395) * MDP(19) + (-t339 * t308 + (t215 * t317 + t216 * t320) * t309 - t252 * t346 + (-t338 * t309 - t376) * qJD(5) + t401) * MDP(20) + t400; (-t280 ^ 2 - t282 ^ 2) * MDP(12) + (t224 * t282 + t225 * t280 + t330 - t402) * MDP(13) + (t206 + t382) * MDP(19) + (t205 + t381) * MDP(20); (t360 - t381) * MDP(16) + (-t356 + t382) * MDP(17) + (-t349 * t309 + t395) * MDP(19) + ((-t214 * t317 + t376) * t309 + t401) * MDP(20) + (-MDP(16) * t380 - t346 * MDP(17) + t349 * MDP(19) - MDP(20) * t376) * qJD(5) + t400;];
tau = t1;
