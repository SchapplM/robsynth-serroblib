% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:32
% EndTime: 2019-12-05 18:28:39
% DurationCPUTime: 3.54s
% Computational Cost: add. (2909->257), mult. (7828->353), div. (0->0), fcn. (6071->8), ass. (0->146)
t366 = sin(pkin(9));
t367 = cos(pkin(9));
t370 = sin(qJ(2));
t373 = cos(qJ(2));
t348 = t366 * t373 + t367 * t370;
t338 = t348 * qJD(2);
t328 = qJD(1) * t338;
t404 = qJD(1) * qJD(2);
t401 = t373 * t404;
t402 = t370 * t404;
t329 = -t366 * t402 + t367 * t401;
t347 = -t366 * t370 + t367 * t373;
t337 = t347 * qJD(1);
t339 = t348 * qJD(1);
t369 = sin(qJ(4));
t372 = cos(qJ(4));
t409 = qJD(4) * t372;
t410 = qJD(4) * t369;
t248 = -t369 * t328 + t372 * t329 + t337 * t409 - t339 * t410;
t385 = t337 * t369 + t372 * t339;
t249 = qJD(4) * t385 + t372 * t328 + t329 * t369;
t296 = t372 * t337 - t339 * t369;
t363 = qJD(2) + qJD(4);
t420 = t296 * t363;
t421 = t385 * t363;
t371 = cos(qJ(5));
t287 = t296 * t371;
t368 = sin(qJ(5));
t408 = qJD(5) * t368;
t231 = qJD(5) * t287 + t371 * t248 - t368 * t249 - t385 * t408;
t251 = t296 * t368 + t371 * t385;
t255 = -t368 * t385 + t287;
t376 = -qJD(5) * t251 - t248 * t368 - t371 * t249;
t362 = qJD(5) + t363;
t422 = t255 * t362;
t423 = t251 * t362;
t441 = (t231 - t422) * MDP(22) + (t376 + t423) * MDP(23) + (t251 ^ 2 - t255 ^ 2) * MDP(21) - t255 * t251 * MDP(20);
t451 = t441 + (-t249 + t421) * MDP(16) + (-t296 ^ 2 + t385 ^ 2) * MDP(14) - t296 * t385 * MDP(13) + (t248 - t420) * MDP(15);
t425 = -qJ(3) - pkin(6);
t356 = t425 * t373;
t353 = qJD(1) * t356;
t342 = t366 * t353;
t355 = t425 * t370;
t352 = qJD(1) * t355;
t424 = qJD(2) * pkin(2);
t346 = t352 + t424;
t298 = t367 * t346 + t342;
t426 = pkin(7) * t339;
t274 = qJD(2) * pkin(3) + t298 - t426;
t419 = t367 * t353;
t299 = t366 * t346 - t419;
t427 = pkin(7) * t337;
t279 = t299 + t427;
t388 = -t274 * t369 - t279 * t372;
t446 = pkin(8) * t296;
t238 = -t388 + t446;
t403 = -pkin(2) * t373 - pkin(1);
t390 = t403 * qJD(1);
t354 = qJD(3) + t390;
t304 = -pkin(3) * t337 + t354;
t265 = -pkin(4) * t296 + t304;
t450 = t238 * t408 - t265 * t255;
t399 = qJD(2) * t425;
t334 = qJD(3) * t373 + t370 * t399;
t316 = t334 * qJD(1);
t335 = -qJD(3) * t370 + t373 * t399;
t317 = t335 * qJD(1);
t280 = -t316 * t366 + t367 * t317;
t263 = -pkin(7) * t329 + t280;
t281 = t367 * t316 + t366 * t317;
t264 = -pkin(7) * t328 + t281;
t380 = -t372 * (qJD(4) * t274 + t264) - t369 * t263 + t279 * t410;
t227 = -pkin(8) * t249 - t380;
t379 = qJD(4) * t388 + t372 * t263 - t369 * t264;
t228 = -pkin(8) * t248 + t379;
t444 = -t368 * t227 + t371 * t228 - t265 * t251;
t445 = (-t238 * t362 - t228) * t368 + t450;
t443 = -t304 * t296 + t380;
t439 = -0.2e1 * t404;
t428 = pkin(4) * t385;
t438 = pkin(8) * t385;
t437 = MDP(4) * t370;
t436 = MDP(5) * (t370 ^ 2 - t373 ^ 2);
t302 = -t352 * t366 + t419;
t282 = t302 - t427;
t303 = t367 * t352 + t342;
t283 = t303 - t426;
t359 = pkin(2) * t367 + pkin(3);
t429 = pkin(2) * t366;
t383 = t359 * t372 - t369 * t429;
t433 = -t383 * qJD(4) + t369 * t282 + t372 * t283;
t333 = t359 * t369 + t372 * t429;
t432 = -t333 * qJD(4) - t372 * t282 + t283 * t369;
t431 = qJD(5) - t362;
t430 = -t304 * t385 + t379;
t374 = qJD(2) ^ 2;
t418 = t370 * t374;
t417 = t371 * t238;
t416 = t373 * t374;
t375 = qJD(1) ^ 2;
t415 = t373 * t375;
t414 = t432 + t446;
t413 = t433 - t438;
t290 = t367 * t334 + t366 * t335;
t307 = t366 * t355 - t367 * t356;
t405 = t370 * qJD(1);
t361 = t370 * t424;
t358 = pkin(2) * t402;
t305 = pkin(3) * t328 + t358;
t312 = pkin(3) * t338 + t361;
t311 = pkin(2) * t405 + pkin(3) * t339;
t394 = t372 * t274 - t279 * t369;
t237 = t394 - t438;
t235 = pkin(4) * t363 + t237;
t400 = -pkin(4) * t362 - t235;
t398 = pkin(1) * t439;
t289 = -t334 * t366 + t367 * t335;
t306 = t367 * t355 + t356 * t366;
t389 = -t368 * t235 - t417;
t291 = -pkin(7) * t348 + t306;
t292 = pkin(7) * t347 + t307;
t387 = -t291 * t369 - t292 * t372;
t301 = t347 * t369 + t348 * t372;
t384 = t372 * t347 - t348 * t369;
t260 = t301 * t368 - t371 * t384;
t261 = t301 * t371 + t368 * t384;
t319 = -pkin(3) * t347 + t403;
t341 = t347 * qJD(2);
t270 = -pkin(7) * t341 + t289;
t271 = -pkin(7) * t338 + t290;
t381 = t369 * t270 + t372 * t271 + t291 * t409 - t292 * t410;
t377 = qJD(4) * t387 + t372 * t270 - t271 * t369;
t332 = pkin(4) + t383;
t275 = -pkin(4) * t384 + t319;
t266 = t311 + t428;
t259 = qJD(4) * t301 + t372 * t338 + t341 * t369;
t258 = qJD(4) * t384 - t338 * t369 + t341 * t372;
t244 = pkin(4) * t259 + t312;
t243 = pkin(4) * t249 + t305;
t242 = pkin(8) * t384 - t387;
t241 = -pkin(8) * t301 + t291 * t372 - t292 * t369;
t234 = qJD(5) * t261 + t258 * t368 + t371 * t259;
t233 = -qJD(5) * t260 + t258 * t371 - t259 * t368;
t230 = -pkin(8) * t258 + t377;
t229 = -pkin(8) * t259 + t381;
t1 = [0.2e1 * t401 * t437 + t436 * t439 + MDP(6) * t416 - MDP(7) * t418 + (-pkin(6) * t416 + t370 * t398) * MDP(9) + (pkin(6) * t418 + t373 * t398) * MDP(10) + (-t280 * t348 + t281 * t347 - t289 * t339 + t290 * t337 - t298 * t341 - t299 * t338 - t306 * t329 - t307 * t328) * MDP(11) + (t280 * t306 + t281 * t307 + t298 * t289 + t299 * t290 + (t354 + t390) * t361) * MDP(12) + (t248 * t301 + t258 * t385) * MDP(13) + (t248 * t384 - t249 * t301 + t258 * t296 - t259 * t385) * MDP(14) + (t319 * t249 + t304 * t259 - t296 * t312 - t305 * t384) * MDP(18) + (t319 * t248 + t304 * t258 + t305 * t301 + t312 * t385) * MDP(19) + (t231 * t261 + t233 * t251) * MDP(20) + (-t231 * t260 + t233 * t255 - t234 * t251 + t261 * t376) * MDP(21) + (t265 * t234 + t243 * t260 - t244 * t255 - t275 * t376) * MDP(25) + (t275 * t231 + t265 * t233 + t243 * t261 + t244 * t251) * MDP(26) + (t258 * MDP(15) - t259 * MDP(16) + t377 * MDP(18) - t381 * MDP(19)) * t363 + (t233 * MDP(22) - t234 * MDP(23) + (-t229 * t368 + t230 * t371 + (-t241 * t368 - t242 * t371) * qJD(5)) * MDP(25) + (-t229 * t371 - t230 * t368 - (t241 * t371 - t242 * t368) * qJD(5)) * MDP(26)) * t362; -t415 * t437 + t375 * t436 + ((t299 + t302) * t339 + (t298 - t303) * t337 + (-t328 * t366 - t329 * t367) * pkin(2)) * MDP(11) + (-t298 * t302 - t299 * t303 + (t280 * t367 + t281 * t366 - t354 * t405) * pkin(2)) * MDP(12) + (t296 * t311 + t363 * t432 + t430) * MDP(18) + (-t311 * t385 + t363 * t433 + t443) * MDP(19) + (t266 * t255 + (t368 * t413 + t371 * t414) * t362 + ((-t332 * t368 - t333 * t371) * t362 + t389) * qJD(5) + t444) * MDP(25) + (-t266 * t251 + (-t228 + (qJD(5) * t333 - t414) * t362) * t368 + (-qJD(5) * t235 - t227 + (-qJD(5) * t332 + t413) * t362) * t371 + t450) * MDP(26) + (t375 * t370 * MDP(9) + MDP(10) * t415) * pkin(1) + t451; (-t337 ^ 2 - t339 ^ 2) * MDP(11) + (t298 * t339 - t299 * t337 + t358) * MDP(12) + (t249 + t421) * MDP(18) + (t248 + t420) * MDP(19) + (-t376 + t423) * MDP(25) + (t231 + t422) * MDP(26); (-t363 * t388 + t430) * MDP(18) + (t363 * t394 + t443) * MDP(19) + (t255 * t428 - (-t237 * t368 - t417) * t362 + (t368 * t400 - t417) * qJD(5) + t444) * MDP(25) + (-t251 * t428 + (qJD(5) * t400 + t237 * t362 - t227) * t371 + t445) * MDP(26) + t451; (t389 * t431 + t444) * MDP(25) + ((-t235 * t431 - t227) * t371 + t445) * MDP(26) + t441;];
tauc = t1;
