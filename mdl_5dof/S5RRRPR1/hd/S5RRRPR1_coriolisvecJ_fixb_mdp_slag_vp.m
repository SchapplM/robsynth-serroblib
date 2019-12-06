% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:58
% EndTime: 2019-12-05 18:39:04
% DurationCPUTime: 3.05s
% Computational Cost: add. (3322->282), mult. (8953->387), div. (0->0), fcn. (6610->8), ass. (0->158)
t373 = qJD(2) + qJD(3);
t382 = cos(qJ(3));
t383 = cos(qJ(2));
t417 = qJD(1) * qJD(2);
t413 = t383 * t417;
t421 = qJD(1) * t383;
t414 = t382 * t421;
t379 = sin(qJ(3));
t380 = sin(qJ(2));
t422 = qJD(1) * t380;
t415 = t379 * t422;
t321 = qJD(3) * t414 - t373 * t415 + t382 * t413;
t351 = t379 * t383 + t380 * t382;
t326 = t373 * t351;
t322 = t326 * qJD(1);
t376 = sin(pkin(9));
t377 = cos(pkin(9));
t280 = -t321 * t376 - t322 * t377;
t281 = t321 * t377 - t322 * t376;
t337 = -t414 + t415;
t339 = -t379 * t421 - t382 * t422;
t317 = t337 * t376 + t339 * t377;
t378 = sin(qJ(5));
t381 = cos(qJ(5));
t419 = qJD(5) * t378;
t396 = -t337 * t377 + t339 * t376;
t431 = t381 * t396;
t237 = qJD(5) * t431 + t378 * t280 + t381 * t281 + t317 * t419;
t454 = -t381 * t317 + t378 * t396;
t238 = qJD(5) * t454 - t381 * t280 + t281 * t378;
t275 = t317 * t378 + t431;
t372 = qJD(5) + t373;
t438 = t275 * t372;
t439 = t454 * t372;
t457 = (-t238 + t439) * MDP(23) + (-t275 ^ 2 + t454 ^ 2) * MDP(21) - t275 * t454 * MDP(20) + (t237 - t438) * MDP(22);
t335 = t339 * qJ(4);
t446 = pkin(6) + pkin(7);
t359 = t446 * t383;
t354 = qJD(1) * t359;
t340 = t379 * t354;
t358 = t446 * t380;
t352 = qJD(1) * t358;
t441 = qJD(2) * pkin(2);
t346 = -t352 + t441;
t405 = t382 * t346 - t340;
t303 = t335 + t405;
t294 = pkin(3) * t373 + t303;
t344 = t382 * t354;
t395 = -t346 * t379 - t344;
t440 = qJ(4) * t337;
t304 = -t395 - t440;
t434 = t377 * t304;
t259 = t376 * t294 + t434;
t443 = pkin(8) * t396;
t245 = t259 + t443;
t369 = -pkin(2) * t383 - pkin(1);
t357 = t369 * qJD(1);
t327 = pkin(3) * t337 + qJD(4) + t357;
t287 = -pkin(4) * t396 + t327;
t407 = t245 * t419 - t287 * t275;
t416 = qJD(2) * t446;
t400 = qJD(1) * t416;
t348 = t383 * t400;
t420 = qJD(3) * t379;
t403 = -t379 * t348 - t354 * t420;
t347 = t380 * t400;
t449 = (qJD(3) * t346 - t347) * t382;
t253 = -qJ(4) * t322 - qJD(4) * t337 + t403 + t449;
t404 = t379 * t347 - t382 * t348;
t388 = qJD(3) * t395 + t404;
t254 = -qJ(4) * t321 + qJD(4) * t339 + t388;
t235 = -t253 * t376 + t377 * t254;
t230 = -pkin(8) * t281 + t235;
t236 = t377 * t253 + t376 * t254;
t231 = pkin(8) * t280 + t236;
t392 = t381 * t230 - t378 * t231 - t287 * t454;
t314 = t317 * pkin(8);
t453 = -0.2e1 * t417;
t452 = MDP(4) * t380;
t451 = MDP(5) * (t380 ^ 2 - t383 ^ 2);
t402 = t352 * t379 - t344;
t306 = t402 + t440;
t424 = -t382 * t352 - t340;
t307 = t335 + t424;
t433 = t377 * t379;
t442 = pkin(2) * qJD(3);
t426 = -t377 * t306 + t307 * t376 + (-t376 * t382 - t433) * t442;
t435 = t376 * t379;
t425 = -t376 * t306 - t377 * t307 + (t377 * t382 - t435) * t442;
t448 = qJD(1) * t351;
t447 = qJD(5) - t372;
t445 = pkin(3) * t339;
t444 = pkin(3) * t376;
t437 = t357 * t339;
t436 = t358 * t382;
t295 = t376 * t304;
t384 = qJD(2) ^ 2;
t432 = t380 * t384;
t430 = t383 * t384;
t385 = qJD(1) ^ 2;
t429 = t383 * t385;
t428 = t443 + t426;
t427 = t314 - t425;
t350 = t379 * t380 - t382 * t383;
t353 = t380 * t416;
t355 = t383 * t416;
t390 = -qJD(3) * t436 - t382 * t353 - t379 * t355 - t359 * t420;
t266 = -qJ(4) * t326 - qJD(4) * t350 + t390;
t325 = t373 * t350;
t394 = t358 * t379 - t359 * t382;
t387 = qJD(3) * t394 + t353 * t379 - t382 * t355;
t267 = qJ(4) * t325 - qJD(4) * t351 + t387;
t242 = t377 * t266 + t376 * t267;
t261 = t377 * t303 - t295;
t315 = -qJ(4) * t351 - t359 * t379 - t436;
t316 = -qJ(4) * t350 - t394;
t271 = t376 * t315 + t377 * t316;
t371 = t380 * t441;
t370 = pkin(2) * t422;
t412 = -pkin(2) * t373 - t346;
t411 = pkin(3) * t322 + qJD(2) * t370;
t410 = pkin(3) * t326 + t371;
t409 = pkin(1) * t453;
t241 = -t266 * t376 + t377 * t267;
t258 = t377 * t294 - t295;
t260 = -t303 * t376 - t434;
t270 = t377 * t315 - t316 * t376;
t368 = pkin(2) * t382 + pkin(3);
t333 = -pkin(2) * t435 + t377 * t368;
t291 = -pkin(4) * t317 - t445;
t243 = pkin(4) * t373 + t258 + t314;
t399 = -t378 * t243 - t381 * t245;
t398 = t258 * t396 - t259 * t317;
t323 = -t350 * t377 - t351 * t376;
t324 = -t350 * t376 + t351 * t377;
t397 = t323 * t381 - t324 * t378;
t283 = t323 * t378 + t324 * t381;
t393 = pkin(3) * t350 + t369;
t391 = t357 * t337 - t403;
t386 = -t339 * t337 * MDP(11) + t321 * MDP(13) + (-t337 ^ 2 + t339 ^ 2) * MDP(12) + (t337 * MDP(13) + (-t339 - t448) * MDP(14)) * t373 + t457;
t366 = pkin(3) * t377 + pkin(4);
t334 = pkin(2) * t433 + t368 * t376;
t328 = pkin(4) + t333;
t299 = -pkin(4) * t323 + t393;
t288 = t291 + t370;
t286 = -t325 * t377 - t326 * t376;
t285 = t325 * t376 - t326 * t377;
t265 = -pkin(4) * t285 + t410;
t257 = pkin(8) * t323 + t271;
t256 = -pkin(8) * t324 + t270;
t255 = -pkin(4) * t280 + t411;
t247 = t314 + t261;
t246 = t260 - t443;
t240 = qJD(5) * t283 - t381 * t285 + t286 * t378;
t239 = qJD(5) * t397 + t285 * t378 + t286 * t381;
t234 = pkin(8) * t285 + t242;
t233 = -pkin(8) * t286 + t241;
t1 = [0.2e1 * t413 * t452 + t451 * t453 + MDP(6) * t430 - MDP(7) * t432 + (-pkin(6) * t430 + t380 * t409) * MDP(9) + (pkin(6) * t432 + t383 * t409) * MDP(10) + (t321 * t351 + t325 * t339) * MDP(11) + (-t321 * t350 - t322 * t351 + t325 * t337 + t339 * t326) * MDP(12) + (t369 * t322 + t357 * t326 + (qJD(1) * t350 + t337) * t371) * MDP(16) + (t369 * t321 - t357 * t325 + (-t339 + t448) * t371) * MDP(17) + (-t235 * t324 + t236 * t323 + t241 * t317 + t242 * t396 - t258 * t286 + t259 * t285 - t270 * t281 + t271 * t280) * MDP(18) + (t235 * t270 + t236 * t271 + t258 * t241 + t259 * t242 + t327 * t410 + t393 * t411) * MDP(19) + (t237 * t283 + t239 * t454) * MDP(20) + (t237 * t397 - t238 * t283 + t239 * t275 - t240 * t454) * MDP(21) + (t299 * t238 + t287 * t240 - t255 * t397 - t265 * t275) * MDP(25) + (t299 * t237 + t287 * t239 + t255 * t283 + t265 * t454) * MDP(26) + (-t325 * MDP(13) - t326 * MDP(14) + t387 * MDP(16) - t390 * MDP(17)) * t373 + (t239 * MDP(22) - t240 * MDP(23) + (t233 * t381 - t234 * t378 + (-t256 * t378 - t257 * t381) * qJD(5)) * MDP(25) + (-t233 * t378 - t234 * t381 - (t256 * t381 - t257 * t378) * qJD(5)) * MDP(26)) * t372; (t288 * t275 + (t427 * t378 + t428 * t381) * t372 + ((-t328 * t378 - t334 * t381) * t372 + t399) * qJD(5) + t392) * MDP(25) + (-t288 * t454 + (-t230 + (qJD(5) * t334 - t428) * t372) * t378 + (-qJD(5) * t243 - t231 + (-qJD(5) * t328 + t427) * t372) * t381 + t407) * MDP(26) + (t236 * t334 + t235 * t333 - t327 * (t370 - t445) + t425 * t259 + t426 * t258) * MDP(19) - t429 * t452 + t386 + (t339 * t370 + t424 * t373 + (qJD(3) * t412 + t347) * t382 + t391) * MDP(17) + (-t337 * t370 + t437 - t402 * t373 + (t379 * t412 - t344) * qJD(3) + t404) * MDP(16) + (t280 * t334 - t281 * t333 + t317 * t426 + t396 * t425 + t398) * MDP(18) + t385 * t451 + (MDP(9) * t380 * t385 + MDP(10) * t429) * pkin(1); (-t373 * t395 + t388 + t437) * MDP(16) + (t373 * t405 + t391 - t449) * MDP(17) + (-t260 * t317 - t261 * t396 + (t280 * t376 - t281 * t377) * pkin(3) + t398) * MDP(18) + (-t258 * t260 - t259 * t261 + (t235 * t377 + t236 * t376 + t327 * t339) * pkin(3)) * MDP(19) + (t291 * t275 - (t246 * t381 - t247 * t378) * t372 + ((-t366 * t378 - t381 * t444) * t372 + t399) * qJD(5) + t392) * MDP(25) + (-t381 * t231 - t378 * t230 - t291 * t454 + (t246 * t378 + t247 * t381) * t372 + (-(t366 * t381 - t378 * t444) * t372 - t381 * t243) * qJD(5) + t407) * MDP(26) + t386; (-t317 ^ 2 - t396 ^ 2) * MDP(18) + (-t258 * t317 - t259 * t396 + t411) * MDP(19) + (t238 + t439) * MDP(25) + (t237 + t438) * MDP(26); (t447 * t399 + t392) * MDP(25) + ((-t245 * t372 - t230) * t378 + (-t447 * t243 - t231) * t381 + t407) * MDP(26) + t457;];
tauc = t1;
