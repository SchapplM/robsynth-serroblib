% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:47
% EndTime: 2019-03-09 02:42:53
% DurationCPUTime: 4.28s
% Computational Cost: add. (2763->405), mult. (5812->499), div. (0->0), fcn. (3998->14), ass. (0->184)
t362 = sin(pkin(10));
t368 = sin(qJ(3));
t433 = qJD(1) * t368;
t364 = cos(pkin(10));
t371 = cos(qJ(3));
t439 = t364 * t371;
t315 = -qJD(1) * t439 + t362 * t433;
t367 = sin(qJ(6));
t370 = cos(qJ(6));
t297 = qJD(3) * t367 - t370 * t315;
t324 = t362 * t371 + t364 * t368;
t318 = t324 * qJD(1);
t473 = qJD(6) + t318;
t474 = t297 * t473;
t412 = t370 * t473;
t299 = qJD(3) * t370 + t315 * t367;
t411 = t473 * t299;
t447 = t473 * t367;
t363 = sin(pkin(9));
t342 = pkin(1) * t363 + pkin(7);
t438 = qJ(4) + t342;
t472 = MDP(12) + MDP(14);
t358 = qJ(1) + pkin(9);
t349 = sin(t358);
t351 = cos(t358);
t405 = g(1) * t349 - g(2) * t351;
t365 = cos(pkin(9));
t344 = -pkin(1) * t365 - pkin(2);
t355 = t371 * pkin(3);
t470 = t344 - t355;
t424 = qJD(1) * qJD(3);
t416 = t368 * t424;
t332 = t362 * t416;
t415 = t371 * t424;
t290 = qJDD(1) * t324 + t364 * t415 - t332;
t469 = -qJ(5) * t290 - qJD(5) * t318;
t357 = qJ(3) + pkin(10);
t348 = sin(t357);
t350 = cos(t357);
t468 = pkin(4) * t350 + qJ(5) * t348;
t406 = g(1) * t351 + g(2) * t349;
t410 = t438 * qJD(1);
t301 = qJD(2) * t368 + t371 * t410;
t292 = t362 * t301;
t300 = t371 * qJD(2) - t368 * t410;
t271 = t300 * t364 - t292;
t427 = -qJD(5) + t271;
t295 = qJD(3) * pkin(3) + t300;
t440 = t364 * t301;
t267 = t362 * t295 + t440;
t261 = -qJD(3) * qJ(5) - t267;
t461 = pkin(5) * t315;
t250 = -t261 - t461;
t270 = t300 * t362 + t440;
t288 = qJDD(6) + t290;
t343 = -pkin(3) * t364 - pkin(4);
t336 = -pkin(8) + t343;
t467 = t336 * t288 + (t250 - t270 + t461) * t473;
t317 = t324 * qJD(3);
t441 = t362 * t368;
t323 = -t439 + t441;
t428 = qJD(6) * t370;
t393 = t317 * t367 + t323 * t428;
t446 = t323 * t367;
t466 = -t288 * t446 - t393 * t473;
t313 = t318 ^ 2;
t464 = pkin(4) + pkin(8);
t421 = qJDD(1) * t371;
t333 = t364 * t421;
t422 = qJDD(1) * t368;
t289 = qJD(1) * t317 + t362 * t422 - t333;
t463 = pkin(4) * t289;
t460 = pkin(5) * t318;
t339 = g(3) * t350;
t456 = g(3) * t371;
t453 = qJDD(3) * pkin(4);
t417 = t370 * qJDD(3) + t367 * t289 + t315 * t428;
t423 = qJD(3) * qJD(6);
t254 = -t367 * t423 + t417;
t452 = t254 * t370;
t391 = -qJ(5) * t324 + t470;
t272 = t323 * t464 + t391;
t451 = t272 * t288;
t450 = t288 * t367;
t449 = t297 * t315;
t448 = t299 * t315;
t445 = t349 * t367;
t444 = t349 * t370;
t443 = t351 * t367;
t442 = t351 * t370;
t282 = t370 * t288;
t437 = qJDD(2) - g(3);
t432 = qJD(3) * t368;
t320 = qJD(3) * t439 - t362 * t432;
t436 = t254 * t324 + t299 * t320;
t353 = t371 * qJDD(2);
t328 = t342 * qJDD(1);
t384 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t328;
t395 = t410 * qJD(3);
t264 = qJDD(3) * pkin(3) - t368 * t384 - t371 * t395 + t353;
t268 = (qJDD(2) - t395) * t368 + t384 * t371;
t246 = t362 * t264 + t364 * t268;
t345 = t355 + pkin(2);
t372 = cos(qJ(1));
t435 = t372 * pkin(1) + t351 * t345;
t360 = t368 ^ 2;
t434 = -t371 ^ 2 + t360;
t331 = qJD(1) * t344;
t314 = qJD(1) * t470 + qJD(4);
t380 = -qJ(5) * t318 + t314;
t258 = t315 * t464 + t380;
t430 = qJD(6) * t258;
t429 = qJD(6) * t323;
t426 = t460 - t427;
t420 = pkin(3) * t416 + qJDD(4);
t347 = pkin(3) * t432;
t418 = qJDD(3) * qJ(5) + t246;
t414 = pkin(3) * t433 + qJ(5) * t315;
t245 = t264 * t364 - t362 * t268;
t266 = t295 * t364 - t292;
t413 = qJD(3) * t438;
t302 = qJD(4) * t371 - t368 * t413;
t303 = -qJD(4) * t368 - t371 * t413;
t274 = t302 * t362 - t364 * t303;
t321 = t438 * t368;
t322 = t438 * t371;
t285 = t364 * t321 + t322 * t362;
t383 = qJDD(1) * t470 + t420;
t377 = t383 + t469;
t242 = t289 * t464 + t377;
t403 = qJD(5) - t266;
t248 = -qJD(3) * t464 + t403 + t460;
t409 = qJD(6) * t248 + t242;
t408 = MDP(23) * t473;
t369 = sin(qJ(1));
t404 = g(1) * t369 - g(2) * t372;
t366 = -qJ(4) - pkin(7);
t402 = -pkin(1) * t369 - t351 * t366;
t400 = qJDD(5) - t245;
t399 = -t430 + t339;
t241 = t248 * t367 + t258 * t370;
t284 = t370 * t289;
t255 = qJD(6) * t299 + qJDD(3) * t367 - t284;
t397 = -t255 * t324 - t297 * t320;
t275 = t302 * t364 + t303 * t362;
t286 = -t321 * t362 + t322 * t364;
t394 = t420 - t405;
t243 = -qJD(3) * qJD(5) - t418;
t392 = -qJ(5) * t320 - qJD(5) * t324 + t347;
t390 = t323 * t282 + t317 * t412 - t429 * t447;
t239 = -pkin(5) * t289 - t243;
t277 = pkin(5) * t324 + t285;
t388 = t239 * t323 + t250 * t317 - t277 * t288;
t387 = -qJD(1) * t331 - t328 + t406;
t386 = 0.2e1 * qJD(3) * t331 - qJDD(3) * t342;
t385 = -g(3) * t348 - t350 * t406;
t373 = qJD(3) ^ 2;
t381 = -0.2e1 * qJDD(1) * t344 - t342 * t373 + t405;
t379 = t274 * t318 - t275 * t315 + t285 * t290 - t286 * t289 - t406;
t276 = pkin(4) * t315 + t380;
t378 = t276 * t318 - t348 * t406 + t339 + t400;
t376 = t239 + (-qJD(6) * t336 + t318 * t464 + t414) * t473 + t385;
t338 = pkin(3) * t362 + qJ(5);
t327 = qJDD(3) * t371 - t368 * t373;
t326 = qJDD(3) * t368 + t371 * t373;
t311 = qJD(3) * t315;
t307 = -t348 * t445 + t442;
t306 = t348 * t444 + t443;
t305 = t348 * t443 + t444;
t304 = t348 * t442 - t445;
t281 = pkin(4) * t323 + t391;
t279 = pkin(4) * t318 + t414;
t278 = -pkin(5) * t323 + t286;
t269 = pkin(4) * t317 + t392;
t260 = -qJD(3) * pkin(4) + t403;
t257 = -pkin(5) * t317 + t275;
t256 = pkin(5) * t320 + t274;
t251 = t317 * t464 + t392;
t247 = t377 + t463;
t244 = t400 - t453;
t240 = t248 * t370 - t258 * t367;
t238 = pkin(5) * t290 - qJDD(3) * t464 + t400;
t237 = t370 * t238;
t1 = [qJDD(1) * MDP(1) + t404 * MDP(2) + (g(1) * t372 + g(2) * t369) * MDP(3) + (t404 + (t363 ^ 2 + t365 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t360 + 0.2e1 * t368 * t415) * MDP(5) + 0.2e1 * (t368 * t421 - t424 * t434) * MDP(6) + t326 * MDP(7) + t327 * MDP(8) + (t368 * t386 + t371 * t381) * MDP(10) + (-t368 * t381 + t371 * t386) * MDP(11) + (-t245 * t324 - t246 * t323 - t266 * t320 - t267 * t317 + t379) * MDP(12) + (t246 * t286 + t267 * t275 - t245 * t285 - t266 * t274 + t383 * t470 + t314 * t347 - g(1) * (-t345 * t349 + t402) - g(2) * (-t349 * t366 + t435)) * MDP(13) + (t243 * t323 + t244 * t324 + t260 * t320 + t261 * t317 + t379) * MDP(14) + (qJD(3) * t274 + qJDD(3) * t285 - t247 * t323 - t269 * t315 - t276 * t317 - t281 * t289 - t350 * t405) * MDP(15) + (qJD(3) * t275 + qJDD(3) * t286 - t247 * t324 - t269 * t318 - t276 * t320 - t281 * t290 + t348 * t405) * MDP(16) + (t247 * t281 + t276 * t269 - t243 * t286 - t261 * t275 + t244 * t285 + t260 * t274 - g(1) * t402 - g(2) * (t351 * t468 + t435) + (-g(1) * (-t345 - t468) + g(2) * t366) * t349) * MDP(17) + (t254 * t446 + t299 * t393) * MDP(18) + ((-t297 * t367 + t299 * t370) * t317 + (t452 - t255 * t367 + (-t297 * t370 - t299 * t367) * qJD(6)) * t323) * MDP(19) + (t436 - t466) * MDP(20) + (t390 + t397) * MDP(21) + (t288 * t324 + t320 * t473) * MDP(22) + (-g(1) * t307 - g(2) * t305 + t237 * t324 + t240 * t320 + t278 * t255 + t257 * t297 + (-t242 * t324 - t251 * t473 - t451) * t367 + (t256 * t473 - t388) * t370 + ((-t272 * t370 - t277 * t367) * t473 - t241 * t324 + t250 * t446) * qJD(6)) * MDP(23) + (g(1) * t306 - g(2) * t304 - t241 * t320 + t278 * t254 + t257 * t299 + (-(qJD(6) * t277 + t251) * t473 - t451 - t409 * t324 + t250 * t429) * t370 + (-(-qJD(6) * t272 + t256) * t473 - (t238 - t430) * t324 + t388) * t367) * MDP(24); t437 * MDP(4) + t327 * MDP(10) - t326 * MDP(11) + (-t245 * t323 + t246 * t324 - t266 * t317 + t267 * t320 - g(3)) * MDP(13) + (qJD(3) * t317 + qJDD(3) * t323) * MDP(15) + (qJD(3) * t320 + qJDD(3) * t324) * MDP(16) + (-t243 * t324 + t244 * t323 + t260 * t317 - t261 * t320 - g(3)) * MDP(17) + (t390 - t397) * MDP(23) + (t436 + t466) * MDP(24) + t472 * (-t289 * t324 + t290 * t323 - t315 * t320 + t317 * t318); MDP(7) * t422 + MDP(8) * t421 + qJDD(3) * MDP(9) + (t368 * t387 + t353 - t456) * MDP(10) + (-t368 * t437 + t371 * t387) * MDP(11) + ((t267 - t270) * t318 + (-t266 + t271) * t315 + (-t289 * t362 - t290 * t364) * pkin(3)) * MDP(12) + (t266 * t270 - t267 * t271 + (-t456 + t245 * t364 + t246 * t362 + (-qJD(1) * t314 + t406) * t368) * pkin(3)) * MDP(13) + (-t289 * t338 + t290 * t343 + (-t261 - t270) * t318 + (t260 + t427) * t315) * MDP(14) + (-qJD(3) * t270 + t279 * t315 + (-pkin(4) + t343) * qJDD(3) + t378) * MDP(15) + (qJDD(3) * t338 - t276 * t315 + t279 * t318 + (0.2e1 * qJD(5) - t271) * qJD(3) + t385 + t418) * MDP(16) + (-t243 * t338 + t244 * t343 - t276 * t279 - t260 * t270 - g(3) * (t355 + t468) + t427 * t261 + t406 * (pkin(3) * t368 + pkin(4) * t348 - qJ(5) * t350)) * MDP(17) + (-t367 * t411 + t452) * MDP(18) + ((-t255 - t411) * t370 + (-t254 + t474) * t367) * MDP(19) + (-t447 * t473 + t282 + t448) * MDP(20) + (-t412 * t473 - t449 - t450) * MDP(21) + t473 * t315 * MDP(22) + (t240 * t315 + t338 * t255 + t426 * t297 + t376 * t367 + t370 * t467) * MDP(23) + (-t241 * t315 + t338 * t254 + t426 * t299 - t367 * t467 + t376 * t370) * MDP(24) + (-MDP(5) * t368 * t371 + MDP(6) * t434) * qJD(1) ^ 2; (t266 * t318 + t267 * t315 + t394) * MDP(13) + t333 * MDP(15) + (t311 + t332) * MDP(16) + (-t260 * t318 - t261 * t315 + t394 + t463 + t469) * MDP(17) + (t449 - t450) * MDP(23) + (-t282 + t448) * MDP(24) + (MDP(24) * t447 - t370 * t408) * t473 + (-t318 * MDP(15) + (-MDP(15) * t324 - MDP(16) * t439) * qJD(1)) * qJD(3) + (-MDP(15) * t441 - t324 * MDP(16) + (MDP(13) + MDP(17)) * t470) * qJDD(1) + t472 * (-t315 ^ 2 - t313); (t311 + t290) * MDP(14) + (-t315 * t318 + qJDD(3)) * MDP(15) + (-t313 - t373) * MDP(16) + (qJD(3) * t261 + t378 - t453) * MDP(17) + (-qJD(3) * t297 + t282) * MDP(23) + (-qJD(3) * t299 - t450) * MDP(24) + (-MDP(24) * t412 - t367 * t408) * t473; t299 * t297 * MDP(18) + (-t297 ^ 2 + t299 ^ 2) * MDP(19) + (t417 + t474) * MDP(20) + (t284 + t411) * MDP(21) + t288 * MDP(22) + (-g(1) * t304 - g(2) * t306 + t241 * t473 - t250 * t299 + t237) * MDP(23) + (g(1) * t305 - g(2) * t307 + t240 * t473 + t250 * t297) * MDP(24) + (-MDP(21) * t423 + MDP(23) * t399 - MDP(24) * t409) * t370 + (-MDP(20) * t423 + (-qJD(6) * t315 - qJDD(3)) * MDP(21) - t409 * MDP(23) + (-t238 - t399) * MDP(24)) * t367;];
tau  = t1;
