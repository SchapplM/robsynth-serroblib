% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:54
% EndTime: 2019-03-09 02:08:59
% DurationCPUTime: 3.24s
% Computational Cost: add. (1899->370), mult. (3319->457), div. (0->0), fcn. (1837->6), ass. (0->161)
t303 = pkin(1) + qJ(3);
t305 = sin(qJ(4));
t308 = cos(qJ(4));
t340 = pkin(4) * t305 - pkin(8) * t308;
t268 = t340 + t303;
t239 = t268 * qJD(1) - qJD(2);
t281 = qJ(2) * qJD(1) + qJD(3);
t274 = -pkin(7) * qJD(1) + t281;
t267 = t305 * t274;
t245 = qJD(4) * pkin(8) + t267;
t304 = sin(qJ(5));
t307 = cos(qJ(5));
t224 = t239 * t304 + t245 * t307;
t341 = pkin(4) * t308 + pkin(8) * t305;
t258 = t341 * qJD(4) + qJD(3);
t294 = qJDD(1) * qJ(3);
t300 = qJDD(1) * pkin(1);
t359 = qJDD(2) - t300;
t350 = -t294 + t359;
t228 = t258 * qJD(1) + t340 * qJDD(1) - t350;
t226 = t307 * t228;
t366 = t307 * qJD(4);
t368 = qJD(5) * t304;
t318 = t305 * t366 + t308 * t368;
t361 = qJDD(1) * t308;
t229 = t318 * qJD(1) - qJD(5) * t366 - t304 * qJDD(4) - t307 * t361;
t295 = qJD(1) * qJD(2);
t296 = qJ(2) * qJDD(1);
t349 = qJDD(3) + t295 + t296;
t264 = -pkin(7) * qJDD(1) + t349;
t369 = qJD(4) * t308;
t234 = qJDD(4) * pkin(8) + t264 * t305 + t274 * t369;
t365 = qJD(1) * qJD(4);
t352 = t308 * t365;
t362 = qJDD(1) * t305;
t257 = qJDD(5) + t352 + t362;
t371 = qJD(4) * t304;
t375 = qJD(1) * t308;
t262 = t307 * t375 + t371;
t214 = pkin(5) * t257 + qJ(6) * t229 - t224 * qJD(5) - qJD(6) * t262 - t234 * t304 + t226;
t376 = qJD(1) * t305;
t354 = t304 * t376;
t230 = -qJD(4) * t354 + t262 * qJD(5) - t307 * qJDD(4) + t304 * t361;
t260 = t304 * t375 - t366;
t367 = qJD(5) * t307;
t357 = -t304 * t228 - t307 * t234 - t239 * t367;
t320 = -t245 * t368 - t357;
t215 = -qJ(6) * t230 - qJD(6) * t260 + t320;
t309 = cos(qJ(1));
t292 = g(2) * t309;
t306 = sin(qJ(1));
t379 = g(1) * t306 - t292;
t422 = -t214 * t307 - t215 * t304 - t379;
t302 = -pkin(7) + qJ(2);
t353 = t302 * t369;
t421 = qJD(2) * t305 + t353;
t386 = t307 * t309;
t393 = t304 * t306;
t247 = t305 * t393 - t386;
t388 = t306 * t307;
t389 = t305 * t309;
t249 = -t304 * t389 - t388;
t420 = -g(1) * t249 + g(2) * t247;
t419 = -g(1) * t309 - g(2) * t306;
t418 = qJD(1) * t303;
t370 = qJD(4) * t305;
t346 = -qJDD(4) * pkin(4) + t274 * t370;
t396 = t264 * t308;
t233 = t346 - t396;
t278 = qJD(5) + t376;
t409 = g(3) * t305;
t417 = qJD(5) * pkin(8) * t278 - t419 * t308 + t233 - t409;
t275 = -qJD(2) + t418;
t416 = (qJD(2) + t275 + t418) * qJD(4) + qJDD(4) * t302;
t223 = t307 * t239 - t245 * t304;
t220 = -qJ(6) * t262 + t223;
t219 = pkin(5) * t278 + t220;
t221 = -qJ(6) * t260 + t224;
t333 = t219 * t304 - t221 * t307;
t398 = t262 * t304;
t400 = t260 * t307;
t313 = (-t398 + t400) * MDP(24) + t333 * MDP(25) + (t304 * MDP(22) + t307 * MDP(23)) * t278;
t415 = t262 ^ 2;
t289 = 0.2e1 * t295;
t414 = pkin(5) * t304;
t408 = g(3) * t308;
t407 = qJ(6) + pkin(8);
t406 = t229 * t304;
t405 = t230 * t307;
t404 = t257 * t304;
t403 = t257 * t307;
t402 = t260 * t278;
t401 = t260 * t304;
t399 = t262 * t278;
t397 = t262 * t307;
t395 = t274 * t308;
t394 = t278 * t304;
t392 = t304 * t308;
t391 = t304 * t309;
t390 = t305 * t307;
t387 = t307 * t308;
t385 = -t220 + t219;
t347 = qJD(5) * t407;
t266 = t341 * qJD(1);
t382 = t304 * t266 + t274 * t387;
t384 = -qJ(6) * t354 + qJD(6) * t307 - t304 * t347 - t382;
t252 = t307 * t266;
t383 = -qJD(6) * t304 - t307 * t347 + t274 * t392 - t252 - (pkin(5) * t308 + qJ(6) * t390) * qJD(1);
t381 = t304 * t268 + t302 * t390;
t380 = t309 * pkin(1) + t306 * qJ(2);
t299 = t308 ^ 2;
t378 = t305 ^ 2 - t299;
t310 = qJD(4) ^ 2;
t311 = qJD(1) ^ 2;
t377 = -t310 - t311;
t373 = qJD(4) * t260;
t372 = qJD(4) * t262;
t364 = qJD(3) * qJD(1);
t363 = qJDD(1) * t303;
t358 = qJ(6) * t387;
t355 = t309 * qJ(3) + t380;
t351 = qJDD(2) - t379;
t348 = -t302 + t414;
t345 = t278 * t302 + t245;
t344 = -qJD(5) * t239 - t234;
t343 = t304 * t258 + t268 * t367 + t421 * t307;
t342 = -t300 + t351;
t246 = -qJD(4) * pkin(4) - t395;
t335 = -t214 * t304 + t215 * t307;
t334 = t219 * t307 + t221 * t304;
t280 = pkin(5) * t307 + pkin(4);
t330 = -t280 * t305 + t308 * t407;
t329 = -t294 + t342;
t328 = -MDP(22) * t307 + MDP(23) * t304;
t326 = -t264 - t419;
t324 = t278 * t367 + t404;
t323 = -t278 * t368 + t403;
t322 = t289 + 0.2e1 * t296 + t419;
t319 = pkin(5) * t230 + qJDD(6) + t346;
t317 = -pkin(8) * t257 + t278 * t246;
t316 = qJD(1) * t275 + t326;
t315 = -qJD(6) * t308 + (qJ(6) * qJD(4) - qJD(5) * t302) * t305;
t265 = -t350 + t364;
t314 = -t302 * t310 + t265 + t363 + t364 + t379;
t312 = (t397 + t401) * MDP(24) - t334 * MDP(25) + t328 * t278;
t288 = t309 * qJ(2);
t284 = qJDD(4) * t308;
t271 = t407 * t307;
t270 = t407 * t304;
t256 = t260 ^ 2;
t254 = t307 * t268;
t250 = t305 * t386 - t393;
t248 = -t305 * t388 - t391;
t243 = t307 * t258;
t236 = -qJ(6) * t392 + t381;
t235 = pkin(5) * t260 + qJD(6) + t246;
t232 = -t358 + t254 + (-t302 * t304 + pkin(5)) * t305;
t218 = t319 - t396;
t217 = -qJD(5) * t358 + t315 * t304 + t343;
t216 = pkin(5) * t369 + t243 + t315 * t307 + ((qJ(6) * t308 - t268) * qJD(5) - t421) * t304;
t1 = [(-g(1) * t248 - g(2) * t250 + t243 * t278 + t254 * t257 + (t302 * t373 - t345 * t367 + t226) * t305 + (-qJD(2) * t260 + qJD(4) * t223 - t230 * t302 + t246 * t367) * t308 + ((-qJD(5) * t268 - t353) * t278 + t233 * t308 + (-qJD(2) * t278 - qJD(4) * t246 - t257 * t302 + t344) * t305) * t304) * MDP(22) + (-t343 * t278 - t381 * t257 - g(1) * t247 - g(2) * t249 + (t345 * t368 + (-t246 * t307 + t262 * t302) * qJD(4) + t357) * t305 + (-qJD(2) * t262 - qJD(4) * t224 + t229 * t302 + t233 * t307 - t246 * t368) * t308) * MDP(23) + (-t216 * t262 - t217 * t260 + t229 * t232 - t230 * t236 + t334 * t370 + (t333 * qJD(5) + t422) * t308) * MDP(24) + (t215 * t236 + t221 * t217 + t214 * t232 + t219 * t216 - g(1) * (-pkin(5) * t391 - pkin(7) * t309 + t288) - g(2) * (t280 * t389 + t355) - t235 * t348 * t370 + (t218 * t348 + t235 * (pkin(5) * t367 - qJD(2)) + t407 * t292) * t308 + (-g(1) * (t330 - t303) - g(2) * (-pkin(7) - t414)) * t306) * MDP(25) + qJDD(1) * MDP(1) + (-0.2e1 * t300 + t351) * MDP(4) + t322 * MDP(5) + (-t359 * pkin(1) - g(1) * (-pkin(1) * t306 + t288) - g(2) * t380 + (t289 + t296) * qJ(2)) * MDP(6) + (qJDD(3) + t322) * MDP(7) + (-t329 + t363 + 0.2e1 * t364) * MDP(8) + (t265 * t303 + t275 * qJD(3) + t349 * qJ(2) + t281 * qJD(2) - g(1) * (-t303 * t306 + t288) - g(2) * t355) * MDP(9) + (qJDD(1) * t299 - 0.2e1 * t305 * t352) * MDP(10) + 0.2e1 * (-t305 * t361 + t378 * t365) * MDP(11) + (-t305 * t310 + t284) * MDP(12) + (-qJDD(4) * t305 - t308 * t310) * MDP(13) + (t314 * t305 + t416 * t308) * MDP(15) + (-t416 * t305 + t314 * t308) * MDP(16) + (-t229 * t387 - t318 * t262) * MDP(17) + ((t398 + t400) * t370 + (t406 - t405 + (-t397 + t401) * qJD(5)) * t308) * MDP(18) + ((-t278 * t366 - t229) * t305 + (t323 + t372) * t308) * MDP(19) + ((t278 * t371 - t230) * t305 + (-t324 - t373) * t308) * MDP(20) + (t257 * t305 + t278 * t369) * MDP(21) + t379 * MDP(2) - t419 * MDP(3); t342 * MDP(6) + t329 * MDP(9) + (-t229 * t307 + t230 * t304) * MDP(24) + t422 * MDP(25) + (-MDP(6) * qJ(2) - MDP(5) - MDP(7)) * t311 + t328 * t257 + t313 * qJD(5) + (-t305 * MDP(15) - t308 * MDP(16) + MDP(4) - MDP(8)) * qJDD(1) + ((-qJD(3) - t281) * MDP(9) + (-0.2e1 * MDP(15) * qJD(4) + t260 * MDP(22) + t262 * MDP(23) + t235 * MDP(25)) * t308 + (0.2e1 * qJD(4) * MDP(16) + t313) * t305) * qJD(1); qJDD(1) * MDP(7) - t311 * MDP(8) + (t419 + t349) * MDP(9) + t284 * MDP(15) + t419 * MDP(25) + (-t275 * MDP(9) + t312) * qJD(1) + (t377 * MDP(16) - t230 * MDP(22) + t229 * MDP(23) - t218 * MDP(25) - t313 * qJD(4)) * t308 + (t377 * MDP(15) - qJDD(4) * MDP(16) + (t373 - t404) * MDP(22) + (t372 - t403) * MDP(23) + (-t405 - t406) * MDP(24) + (qJD(4) * t235 + t335) * MDP(25) + t312 * qJD(5)) * t305; MDP(12) * t361 - MDP(13) * t362 + qJDD(4) * MDP(14) + (-t316 * t308 + t409) * MDP(15) + (t316 * t305 + t408) * MDP(16) + (t278 * t397 - t406) * MDP(17) + ((-t229 - t402) * t307 + (-t230 - t399) * t304) * MDP(18) + ((-t262 * t308 + t278 * t390) * qJD(1) + t324) * MDP(19) + ((t260 * t308 - t305 * t394) * qJD(1) + t323) * MDP(20) - t278 * MDP(21) * t375 + (-t223 * t375 - t260 * t267 - pkin(4) * t230 - t252 * t278 + (t278 * t395 + t317) * t304 - t417 * t307) * MDP(22) + (pkin(4) * t229 + t224 * t375 - t262 * t267 + t382 * t278 + t417 * t304 + t317 * t307) * MDP(23) + (-t408 - t229 * t270 - t230 * t271 - t383 * t262 - t384 * t260 - t334 * qJD(5) + (-t334 * qJD(1) + t419) * t305 + t335) * MDP(24) + (t215 * t271 - t214 * t270 - t218 * t280 - g(3) * t330 + (pkin(5) * t394 - t267) * t235 + t384 * t221 + t383 * t219 + t419 * (t280 * t308 + t305 * t407)) * MDP(25) + (t308 * t305 * MDP(10) - t378 * MDP(11)) * t311; t262 * t260 * MDP(17) + (-t256 + t415) * MDP(18) + (-t229 + t402) * MDP(19) + (-t230 + t399) * MDP(20) + t257 * MDP(21) + (-t245 * t367 + t224 * t278 - t246 * t262 + t226 + (t344 + t408) * t304 + t420) * MDP(22) + (g(1) * t250 - g(2) * t248 + g(3) * t387 + t223 * t278 + t246 * t260 - t320) * MDP(23) + (pkin(5) * t229 - t385 * t260) * MDP(24) + (t385 * t221 + (g(3) * t392 - t235 * t262 + t214 + t420) * pkin(5)) * MDP(25); (-t256 - t415) * MDP(24) + (t219 * t262 + t221 * t260 + t326 * t308 + t319 - t409) * MDP(25);];
tau  = t1;
