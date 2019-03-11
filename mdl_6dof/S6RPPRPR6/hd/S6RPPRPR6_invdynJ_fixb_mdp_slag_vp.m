% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:43
% EndTime: 2019-03-09 01:51:46
% DurationCPUTime: 2.84s
% Computational Cost: add. (1197->367), mult. (2007->426), div. (0->0), fcn. (1035->6), ass. (0->169)
t309 = cos(qJ(4));
t372 = qJD(1) * qJD(4);
t352 = t309 * t372;
t306 = sin(qJ(4));
t369 = qJDD(1) * t306;
t425 = t352 + t369;
t280 = qJ(2) * qJD(1) + qJD(3);
t264 = -pkin(7) * qJD(1) + t280;
t234 = (pkin(5) * qJD(1) - t264) * t309;
t374 = qJD(5) + t234;
t377 = qJD(6) * t306;
t424 = qJD(1) * t377 + qJDD(4);
t307 = sin(qJ(1));
t310 = cos(qJ(1));
t394 = g(1) * t307 - g(2) * t310;
t393 = g(1) * t310 + g(2) * t307;
t289 = t306 * pkin(4);
t304 = pkin(1) + qJ(3);
t412 = qJ(5) * t309;
t335 = t412 - t304;
t423 = t335 - t289;
t282 = t309 * qJDD(1);
t353 = t306 * t372;
t422 = -t353 + t282;
t247 = -qJDD(6) - t422;
t308 = cos(qJ(6));
t237 = t308 * t247;
t387 = qJD(1) * t309;
t271 = qJD(6) + t387;
t305 = sin(qJ(6));
t402 = t271 * t305;
t320 = -qJD(6) * t402 - t237;
t421 = qJD(1) * t304;
t257 = t306 * t264;
t388 = qJD(1) * t306;
t233 = -pkin(5) * t388 + t257;
t386 = qJD(4) * qJ(5);
t228 = t233 + t386;
t311 = -pkin(4) - pkin(8);
t420 = -t311 * t247 + (t228 - t233) * t271;
t297 = qJD(1) * qJD(2);
t298 = qJ(2) * qJDD(1);
t349 = qJDD(3) + t297 + t298;
t252 = -pkin(7) * qJDD(1) + t349;
t240 = t306 * t252;
t368 = qJDD(4) * qJ(5);
t341 = -t240 - t368;
t403 = t264 * t309;
t222 = (-qJD(5) - t403) * qJD(4) + t341;
t340 = -t252 * t309 + qJDD(5);
t382 = qJD(4) * t306;
t329 = t264 * t382 + t340;
t411 = qJDD(4) * pkin(4);
t223 = t329 - t411;
t343 = -qJD(5) + t403;
t413 = qJD(4) * pkin(4);
t236 = -t343 - t413;
t239 = -t257 - t386;
t419 = (t236 * t306 - t239 * t309) * qJD(4) - t222 * t306 - t223 * t309;
t265 = -qJD(2) + t421;
t303 = -pkin(7) + qJ(2);
t367 = qJDD(4) * t303;
t418 = (qJD(2) + t265 + t421) * qJD(4) + t367;
t375 = pkin(4) * t388 - qJD(2);
t229 = -t335 * qJD(1) + t375;
t417 = (-qJD(1) * t423 + qJD(2) + t229) * qJD(4) + t367;
t287 = 0.2e1 * t297;
t416 = g(3) * t306;
t415 = g(3) * t309;
t414 = pkin(5) - t303;
t218 = -pkin(5) * t369 + (qJD(5) - t234) * qJD(4) - t341;
t410 = t218 * t306;
t339 = t425 * t305 + t424 * t308;
t371 = qJD(4) * qJD(6);
t220 = -t305 * t371 + t339;
t409 = t220 * t308;
t408 = t247 * t305;
t383 = qJD(4) * t305;
t249 = -t308 * t388 + t383;
t407 = t249 * t271;
t406 = t249 * t306;
t381 = qJD(4) * t308;
t251 = t305 * t388 + t381;
t405 = t251 * t271;
t404 = t251 * t306;
t401 = t305 * t306;
t400 = t307 * t309;
t399 = t308 * t309;
t398 = t309 * t310;
t397 = t425 * t308;
t396 = -g(1) * t398 - g(2) * t400;
t254 = pkin(4) * t387 + qJ(5) * t388;
t395 = t310 * pkin(1) + t307 * qJ(2);
t300 = t306 ^ 2;
t301 = t309 ^ 2;
t392 = -t300 - t301;
t312 = qJD(4) ^ 2;
t313 = qJD(1) ^ 2;
t391 = t312 + t313;
t390 = MDP(27) * t308;
t389 = qJD(1) * t265;
t385 = qJD(4) * t249;
t384 = qJD(4) * t251;
t380 = qJD(4) * t309;
t379 = qJD(5) * t309;
t319 = pkin(8) * t306 - t335;
t224 = t319 * qJD(1) + t375;
t378 = qJD(6) * t224;
t376 = qJD(6) * t308;
t373 = qJ(5) * qJDD(1);
t296 = qJD(3) * qJD(1);
t370 = qJDD(1) * t304;
t366 = qJDD(4) * t306;
t365 = qJDD(4) * t309;
t364 = MDP(15) - MDP(18);
t363 = MDP(16) - MDP(19);
t302 = qJDD(1) * pkin(1);
t362 = t302 - qJDD(2);
t361 = t309 * t402;
t360 = t271 * t399;
t359 = t306 * t309 * t313;
t358 = t310 * qJ(3) + t395;
t357 = t271 * t383;
t356 = t271 * t381;
t354 = pkin(4) * t380 + qJ(5) * t382 + qJD(3);
t350 = qJDD(2) - t394;
t348 = qJD(4) * t414;
t346 = qJD(1) * t254 - g(3);
t345 = qJD(4) * pkin(8) - qJD(5);
t344 = -qJD(4) * t239 - t223;
t295 = qJDD(1) * qJ(3);
t336 = -t295 - t296 - t362;
t317 = -t425 * pkin(4) - qJ(5) * t353 + t336;
t215 = pkin(8) * t369 + (t345 * qJD(1) - t373) * t309 - t317;
t225 = t311 * qJD(4) + t374;
t342 = qJD(6) * t225 + t215;
t338 = 0.2e1 * t298 + t287 - t393;
t337 = -t302 + t350;
t334 = t396 + t416;
t332 = pkin(4) * t309 + qJ(5) * t306;
t328 = -t378 - t416;
t214 = t224 * t308 + t225 * t305;
t325 = -t236 * t309 - t239 * t306;
t323 = -t295 + t337;
t321 = t271 * t376 - t408;
t318 = -t303 * t312 + t394;
t219 = (-qJD(1) * qJD(5) - t373) * t309 - t317;
t232 = t354 - t379;
t316 = -qJD(1) * t232 + qJDD(1) * t423 - t219 - t318;
t315 = -t336 + t318 + t370 + t296;
t314 = -t415 + t218 - t393 * t306 + (pkin(8) * t387 - qJD(6) * t311 + t254) * t271;
t286 = t310 * qJ(2);
t283 = t301 * t313;
t260 = t414 * t309;
t259 = t414 * t306;
t244 = t305 * t400 - t308 * t310;
t243 = t305 * t310 + t307 * t399;
t242 = t305 * t398 + t307 * t308;
t241 = t305 * t307 - t308 * t398;
t238 = t289 + t319;
t231 = qJD(2) * t306 - t309 * t348;
t230 = -qJD(2) * t309 - t306 * t348;
t227 = t345 * t309 + t354;
t226 = t229 * t387;
t221 = t251 * qJD(6) + qJDD(4) * t305 - t397;
t217 = t422 * pkin(5) + t311 * qJDD(4) + t329;
t216 = t308 * t217;
t213 = -t224 * t305 + t225 * t308;
t1 = [((-t227 * t305 + t230 * t308) * t271 - (-t238 * t305 + t260 * t308) * t247 + (-t215 * t305 + t216) * t309 + t231 * t249 - t259 * t221 - t308 * t410 - g(1) * t244 + g(2) * t242 + (-t213 * t306 - t228 * t399) * qJD(4) + ((-t238 * t308 - t260 * t305) * t271 - t214 * t309 + t228 * t401) * qJD(6)) * MDP(26) + (t214 * t382 - g(1) * t243 - g(2) * t241 - t259 * t220 + t231 * t251 + (-(qJD(6) * t260 + t227) * t271 + t238 * t247 - t342 * t309 + t228 * t377) * t308 + (-(-qJD(6) * t238 + t230) * t271 + t260 * t247 + t410 + (qJD(4) * t228 - t217 + t378) * t309) * t305) * MDP(27) + (-t219 * t423 + t229 * t232 - g(1) * (-pkin(7) * t310 + t286) - g(2) * (-qJ(5) * t398 + t310 * t289 + t358) + t325 * qJD(2) + (g(2) * pkin(7) - g(1) * t423) * t307 + t419 * t303) * MDP(20) + (-t336 * t304 + t265 * qJD(3) + t349 * qJ(2) + t280 * qJD(2) - g(1) * (-t304 * t307 + t286) - g(2) * t358) * MDP(9) + qJDD(1) * MDP(1) + (-t247 * t309 - t271 * t382) * MDP(25) + ((t220 + t357) * t309 + (t321 - t384) * t306) * MDP(23) + ((-t221 + t356) * t309 + (t320 + t385) * t306) * MDP(24) + (0.2e1 * t296 - t323 + t370) * MDP(8) + 0.2e1 * (-t306 * t282 + (t300 - t301) * t372) * MDP(11) + (-t306 * t312 + t365) * MDP(12) + (-t309 * t312 - t366) * MDP(13) + t393 * MDP(3) + t394 * MDP(2) + (t417 * t306 + t316 * t309) * MDP(19) + (t316 * t306 - t417 * t309) * MDP(18) + (t315 * t306 + t418 * t309) * MDP(15) + (-t418 * t306 + t315 * t309) * MDP(16) + (t393 + (qJDD(1) * t303 + t297) * t392 - t419) * MDP(17) + (t362 * pkin(1) - g(1) * (-pkin(1) * t307 + t286) - g(2) * t395 + (t287 + t298) * qJ(2)) * MDP(6) + (qJDD(1) * t301 - 0.2e1 * t306 * t352) * MDP(10) + (-0.2e1 * t302 + t350) * MDP(4) + (qJDD(3) + t338) * MDP(7) + t338 * MDP(5) + ((-t249 * t305 + t251 * t308) * t380 + (t409 - t221 * t305 + (-t249 * t308 - t251 * t305) * qJD(6)) * t306) * MDP(22) + (t220 * t401 + (t305 * t380 + t306 * t376) * t251) * MDP(21); t337 * MDP(6) + (-t296 + t323) * MDP(9) + t283 * MDP(17) + (t317 - t394) * MDP(20) + t321 * MDP(26) + t320 * MDP(27) + (MDP(17) * t300 - qJ(2) * MDP(6) - MDP(5) - MDP(7)) * t313 - t364 * (0.2e1 * t352 + t369) + t363 * (-t282 + 0.2e1 * t353) + (MDP(20) * t412 + MDP(4) - MDP(8)) * qJDD(1) + (-t280 * MDP(9) + (-t325 + t379) * MDP(20) + (t360 - t406) * MDP(26) + (-t361 - t404) * MDP(27)) * qJD(1); -t313 * MDP(8) + (t349 - t393) * MDP(9) - t393 * MDP(20) + t364 * (-t391 * t306 + t365) - t363 * (t391 * t309 + t366) + (-t229 * MDP(20) - t265 * MDP(9) + (MDP(26) * t305 + t390) * t271) * qJD(1) + ((qJD(4) * t236 - t222) * MDP(20) + (t221 + t356) * MDP(26) + (t220 - t357) * MDP(27)) * t306 + (t392 * MDP(17) + MDP(7)) * qJDD(1) + (t344 * MDP(20) + (-t320 + t385) * MDP(26) + (t321 + t384) * MDP(27)) * t309; MDP(10) * t359 + (-t300 * t313 + t283) * MDP(11) + MDP(12) * t282 - MDP(13) * t369 + qJDD(4) * MDP(14) + ((t252 - t389) * t309 + t334) * MDP(15) + (t415 - t240 + (t393 + t389) * t306) * MDP(16) + (-t332 * qJDD(1) + ((-t239 - t386) * t309 + (-qJD(5) + t236 + t413) * t306) * qJD(1)) * MDP(17) + (t346 * t306 + t226 + t340 - t396 - 0.2e1 * t411) * MDP(18) + (0.2e1 * t368 + 0.2e1 * qJD(4) * qJD(5) + t240 + t346 * t309 + (-qJD(1) * t229 - t393) * t306) * MDP(19) + (-t222 * qJ(5) - t223 * pkin(4) - t229 * t254 - t236 * t257 - g(3) * (-t289 + t412) + t343 * t239 - t393 * t332) * MDP(20) + (-t251 * t402 + t409) * MDP(21) + ((-t221 - t405) * t308 + (-t220 + t407) * t305) * MDP(22) + ((-t361 + t404) * qJD(1) + t320) * MDP(23) + ((-t360 - t406) * qJD(1) - t321) * MDP(24) + t271 * MDP(25) * t388 + (qJ(5) * t221 + t213 * t388 + t374 * t249 + t314 * t305 + t420 * t308) * MDP(26) + (qJ(5) * t220 - t214 * t388 + t374 * t251 - t420 * t305 + t314 * t308) * MDP(27); MDP(17) * t282 + (qJDD(4) - t359) * MDP(18) + (-t283 - t312) * MDP(19) + (t226 - t334 - t344) * MDP(20) + (-t237 - t385) * MDP(26) + (-t384 + t408) * MDP(27) + (-MDP(26) * t402 - t271 * t390) * t271; t251 * t249 * MDP(21) + (-t249 ^ 2 + t251 ^ 2) * MDP(22) + (t339 + t407) * MDP(23) + (t397 + t405) * MDP(24) - t247 * MDP(25) + (-g(1) * t241 + g(2) * t243 + t214 * t271 - t228 * t251 + t216) * MDP(26) + (-g(1) * t242 - g(2) * t244 + t213 * t271 + t228 * t249) * MDP(27) + (-MDP(24) * t371 + t328 * MDP(26) - t342 * MDP(27)) * t308 + (-MDP(23) * t371 - t424 * MDP(24) - t342 * MDP(26) + (-t217 - t328) * MDP(27)) * t305;];
tau  = t1;
