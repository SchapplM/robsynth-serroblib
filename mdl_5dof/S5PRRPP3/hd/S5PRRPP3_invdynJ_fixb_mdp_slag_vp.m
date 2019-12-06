% Calculate vector of inverse dynamics joint torques for
% S5PRRPP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:45
% EndTime: 2019-12-05 16:13:51
% DurationCPUTime: 3.71s
% Computational Cost: add. (1735->390), mult. (3792->502), div. (0->0), fcn. (2577->8), ass. (0->158)
t333 = sin(pkin(7));
t335 = cos(pkin(7));
t364 = g(1) * t335 + g(2) * t333;
t430 = MDP(12) + MDP(16);
t436 = MDP(14) + MDP(17);
t337 = sin(qJ(2));
t435 = t337 * t364;
t339 = cos(qJ(2));
t391 = qJD(1) * qJD(2);
t291 = qJDD(2) * pkin(6) + qJDD(1) * t337 + t339 * t391;
t336 = sin(qJ(3));
t275 = t336 * t291;
t402 = qJD(1) * t337;
t311 = qJD(2) * pkin(6) + t402;
t338 = cos(qJ(3));
t396 = qJD(3) * t338;
t242 = -qJDD(3) * pkin(3) + t311 * t396 + qJDD(4) + t275;
t414 = t336 * t339;
t280 = t333 * t414 + t335 * t338;
t284 = -t333 * t338 + t335 * t414;
t366 = g(1) * t284 + g(2) * t280;
t415 = t336 * t337;
t348 = g(3) * t415 + t366;
t343 = t242 - t348;
t363 = pkin(3) * t336 - qJ(4) * t338;
t278 = qJD(3) * t363 - qJD(4) * t336;
t332 = sin(pkin(8));
t334 = cos(pkin(8));
t411 = t338 * t339;
t287 = t332 * t337 + t334 * t411;
t434 = -t287 * qJD(1) + t332 * t278;
t384 = t332 * t411;
t433 = qJD(1) * t384 + (t278 - t402) * t334;
t326 = t334 * qJDD(3);
t390 = qJD(2) * qJD(3);
t371 = t338 * t390;
t389 = qJDD(2) * t336;
t353 = t371 + t389;
t257 = t332 * t353 - t326;
t258 = qJDD(3) * t332 + t334 * t353;
t431 = pkin(4) * t257 - qJ(5) * t258;
t386 = MDP(13) - MDP(18);
t398 = qJD(3) * t332;
t400 = qJD(2) * t336;
t297 = t334 * t400 + t398;
t293 = t297 ^ 2;
t428 = pkin(4) * t334;
t427 = pkin(6) * t338;
t426 = pkin(6) * t339;
t423 = g(3) * t339;
t422 = qJD(2) * pkin(2);
t421 = qJ(4) * t334;
t300 = t363 * qJD(2);
t417 = t300 * t334;
t356 = pkin(3) * t338 + qJ(4) * t336 + pkin(2);
t416 = t356 * t334;
t301 = t336 * t311;
t413 = t337 * t334;
t412 = t337 * t338;
t302 = t338 * t311;
t410 = qJDD(1) - g(3);
t324 = t337 * t391;
t369 = -qJDD(1) * t339 + t324;
t229 = qJD(2) * t278 - qJDD(2) * t356 + t369;
t239 = qJDD(3) * qJ(4) + t291 * t338 + (qJD(4) - t301) * qJD(3);
t225 = t332 * t229 + t334 * t239;
t393 = qJD(5) * t338;
t397 = qJD(3) * t336;
t409 = -t393 + (-pkin(6) * t334 + qJ(5)) * t397 + t434;
t380 = pkin(6) * t332 + pkin(4);
t408 = -t380 * t397 - t433;
t385 = pkin(6) * t397;
t407 = t332 * t385 + t433;
t406 = -t334 * t385 + t434;
t401 = qJD(1) * t339;
t268 = -qJD(2) * t356 - t401;
t289 = qJD(3) * qJ(4) + t302;
t232 = t332 * t268 + t334 * t289;
t405 = t364 * t415;
t264 = -t332 * t356 + t334 * t427;
t330 = t336 ^ 2;
t331 = t338 ^ 2;
t404 = t330 - t331;
t340 = qJD(3) ^ 2;
t341 = qJD(2) ^ 2;
t403 = t340 + t341;
t399 = qJD(2) * t338;
t395 = qJD(4) * t334;
t394 = qJD(5) * t297;
t392 = t334 * qJD(3);
t388 = qJDD(2) * t338;
t387 = qJDD(2) * t339;
t383 = t336 * t413;
t373 = qJ(4) * t388;
t375 = t332 * t399;
t381 = g(3) * t383 + qJD(4) * t375 + t332 * t373;
t379 = qJ(4) * t397;
t295 = t332 * t400 - t392;
t378 = t295 * t401;
t377 = t297 * t401;
t376 = t336 * t401;
t374 = t337 * t397;
t372 = t336 * t390;
t224 = t334 * t229 - t332 * t239;
t368 = t339 * pkin(2) + pkin(3) * t411 + t337 * pkin(6) + qJ(4) * t414;
t367 = -pkin(3) * t415 + qJ(4) * t412;
t281 = t333 * t411 - t335 * t336;
t285 = t333 * t336 + t335 * t411;
t365 = g(1) * t285 + g(2) * t281;
t279 = -qJD(3) * pkin(3) + qJD(4) + t301;
t362 = pkin(4) * t332 - qJ(5) * t334;
t312 = -t401 - t422;
t360 = g(3) * t337 - qJD(2) * t312;
t231 = t268 * t334 - t289 * t332;
t359 = qJ(4) * t258 + qJD(4) * t297;
t355 = pkin(6) + t362;
t354 = pkin(4) * t388 + qJDD(5) - t224;
t283 = -t332 * t339 + t334 * t412;
t282 = t332 * t412 + t334 * t339;
t352 = t372 - t388;
t259 = t282 * t333;
t261 = t282 * t335;
t286 = t384 - t413;
t351 = -g(1) * t261 - g(2) * t259 + g(3) * t286;
t260 = t283 * t333;
t262 = t283 * t335;
t350 = g(1) * t262 + g(2) * t260 - g(3) * t287;
t347 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t340 - t369 - t423;
t230 = pkin(4) * t295 - qJ(5) * t297 + t279;
t346 = t279 * MDP(15) + t230 * MDP(19) + t297 * t386;
t345 = -g(3) * t412 - t257 * t421 - t295 * t395 - t365;
t344 = -pkin(6) * qJDD(3) + (t312 + t401 - t422) * qJD(3);
t342 = t356 * t435;
t318 = t335 * t426;
t315 = t333 * t426;
t303 = -qJ(5) * t332 - pkin(3) - t428;
t288 = t332 * t300;
t277 = t284 * pkin(3);
t276 = t280 * pkin(3);
t274 = t295 * t399;
t269 = t355 * t336;
t263 = -t332 * t427 - t416;
t255 = t338 * t380 + t416;
t254 = -qJ(5) * t338 + t264;
t251 = qJD(2) * t287 - t334 * t374;
t250 = -qJD(2) * t413 - t332 * t374 + t339 * t375;
t247 = t362 * t399 + t302;
t246 = -t301 * t334 + t288;
t245 = t301 * t332 + t417;
t244 = -qJD(5) * t334 * t336 + t355 * t396;
t241 = -t417 + (-pkin(4) * qJD(2) - t311 * t332) * t336;
t240 = t288 + (qJ(5) * qJD(2) - t311 * t334) * t336;
t228 = -qJ(5) * t399 + t232;
t227 = pkin(4) * t399 + qJD(5) - t231;
t223 = t242 - t394 + t431;
t222 = -pkin(4) * t372 + t354;
t221 = qJ(5) * t352 - qJD(2) * t393 + t225;
t1 = [t410 * MDP(1) + (-t337 * t341 + t387) * MDP(3) + (-qJDD(2) * t337 - t339 * t341) * MDP(4) + (-t224 * t282 + t225 * t283 - t231 * t250 + t232 * t251 - g(3)) * MDP(15) + (t221 * t283 + t222 * t282 + t227 * t250 + t228 * t251 - g(3)) * MDP(19) + ((qJDD(2) * MDP(10) - 0.2e1 * MDP(11) * t390) * t339 + (-MDP(10) * t403 - qJDD(3) * MDP(11) + qJD(3) * t346) * t337 + t430 * (qJD(2) * t250 + qJDD(2) * t282) + t386 * (qJD(2) * t251 + qJDD(2) * t283)) * t338 + (-MDP(11) * t387 + (-qJDD(3) * MDP(10) + MDP(11) * t403 + t242 * MDP(15) + t223 * MDP(19) + t258 * t386) * t337 + (t346 * t339 + (-0.2e1 * t339 * MDP(10) - t282 * t430 - t386 * t283) * qJD(3)) * qJD(2)) * t336 + t430 * (t257 * t415 + (t337 * t396 + t339 * t400) * t295) + t436 * (t250 * t297 - t251 * t295 - t283 * t257 + t258 * t282); qJDD(2) * MDP(2) + (t339 * t410 + t435) * MDP(3) + (-t337 * t410 + t339 * t364) * MDP(4) + (qJDD(2) * t330 + 0.2e1 * t336 * t371) * MDP(5) + 0.2e1 * (t336 * t388 - t390 * t404) * MDP(6) + (qJDD(3) * t336 + t338 * t340) * MDP(7) + (qJDD(3) * t338 - t336 * t340) * MDP(8) + (t344 * t336 + ((t364 + t391) * t337 + t347) * t338) * MDP(10) + (t344 * t338 + (-t347 - t324) * t336 - t405) * MDP(11) + ((-t378 + pkin(6) * t257 + t242 * t332 + (qJD(2) * t263 + t231) * qJD(3)) * t336 + (-qJDD(2) * t263 - t224 + (pkin(6) * t295 + t279 * t332) * qJD(3) - t407 * qJD(2)) * t338 + t350) * MDP(12) + ((-t377 + pkin(6) * t258 + t242 * t334 + (-qJD(2) * t264 - t232) * qJD(3)) * t336 + (qJDD(2) * t264 + t225 + (pkin(6) * t297 + t279 * t334) * qJD(3) + t406 * qJD(2)) * t338 + t351) * MDP(13) + (-t257 * t264 - t258 * t263 - t407 * t297 - t406 * t295 + (-t231 * t334 - t232 * t332) * t396 + (-t224 * t334 - t225 * t332 - t423) * t336 + t405) * MDP(14) + (t225 * t264 + t224 * t263 - t279 * t376 - g(1) * t318 - g(2) * t315 - g(3) * t368 + t406 * t232 + t407 * t231 + (t242 * t336 + t279 * t396) * pkin(6) + t342) * MDP(15) + (t244 * t295 + t257 * t269 + (-t378 + t223 * t332 + (-qJD(2) * t255 - t227) * qJD(3)) * t336 + (qJD(2) * t408 + qJDD(2) * t255 + t230 * t398 + t222) * t338 + t350) * MDP(16) + (-t254 * t257 + t255 * t258 + t408 * t297 - t409 * t295 + (t227 * t334 - t228 * t332) * t396 + (-t221 * t332 + t222 * t334 - t423) * t336 + t405) * MDP(17) + (-t244 * t297 - t258 * t269 + (t377 - t223 * t334 + (qJD(2) * t254 + t228) * qJD(3)) * t336 + (-qJD(2) * t409 - qJDD(2) * t254 - t230 * t392 - t221) * t338 - t351) * MDP(18) + (t221 * t254 + t223 * t269 + t222 * t255 - g(1) * (-pkin(4) * t262 - qJ(5) * t261 + t318) - g(2) * (-pkin(4) * t260 - qJ(5) * t259 + t315) - g(3) * (pkin(4) * t287 + qJ(5) * t286 + t368) + (t244 - t376) * t230 + t409 * t228 + t408 * t227 + t342) * MDP(19); MDP(7) * t389 + MDP(8) * t388 + qJDD(3) * MDP(9) + (t336 * t360 - t275 + t366) * MDP(10) + ((-t291 + t360) * t338 + t365) * MDP(11) + (-t295 * t302 - pkin(3) * t257 + (-t242 + t366) * t334 + (-t231 * t336 + t245 * t338 + (-t279 * t338 - t379) * t332) * qJD(2) + t381) * MDP(12) + (-pkin(3) * t258 + (qJDD(2) * t421 - t297 * t311) * t338 + t343 * t332 + (t232 * t336 - t246 * t338 + (-t379 + (qJD(4) - t279) * t338) * t334) * qJD(2)) * MDP(13) + (t245 * t297 + t246 * t295 + (t231 * t399 + t225) * t334 + (t232 * t399 - t224 + t359) * t332 + t345) * MDP(14) + (-t242 * pkin(3) - t232 * t246 - t231 * t245 - t279 * t302 + g(1) * t277 + g(2) * t276 - g(3) * t367 + (-t231 * t332 + t232 * t334) * qJD(4) + (-t224 * t332 + t225 * t334 - t365) * qJ(4)) * MDP(15) + (t257 * t303 + (-qJD(5) * t332 - t247) * t295 + (-t223 + t366) * t334 + (t227 * t336 - t241 * t338 + (-t230 * t338 - t379) * t332) * qJD(2) + t381) * MDP(16) + (t240 * t295 - t241 * t297 + (-t227 * t399 + t221) * t334 + (t228 * t399 + t222 + t359) * t332 + t345) * MDP(17) + (-t334 * t373 + t247 * t297 - t258 * t303 + (-t223 + t348 + t394) * t332 + (-t228 * t336 + t240 * t338 + (t379 + (-qJD(4) + t230) * t338) * t334) * qJD(2)) * MDP(18) + (t221 * t421 + t223 * t303 - t230 * t247 - t227 * t241 - g(1) * (qJ(4) * t285 - t284 * t428 - t277) - g(2) * (qJ(4) * t281 - t280 * t428 - t276) - g(3) * (-pkin(4) * t383 + t367) + (-t240 + t395) * t228 + (t222 * qJ(4) + qJ(5) * t348 + t227 * qJD(4) - t230 * qJD(5)) * t332) * MDP(19) + (-MDP(5) * t336 * t338 + MDP(6) * t404) * t341; (t231 * t297 + t232 * t295 + t343) * MDP(15) + (t228 * t295 + (-qJD(5) - t227) * t297 + t343 + t431) * MDP(19) + t386 * (t274 + t258) + t436 * (-t295 ^ 2 - t293) + (t399 * (-t297 + t398) + t332 * t389 - t326) * t430; (t295 * t297 - t352) * MDP(16) + (-t274 + t258) * MDP(17) + (-t331 * t341 - t293) * MDP(18) + (t230 * t297 - g(1) * (t285 * t332 - t335 * t413) - g(2) * (t281 * t332 - t333 * t413) - g(3) * t282 + (-pkin(4) * t397 + t228 * t338) * qJD(2) + t354) * MDP(19);];
tau = t1;
