% Calculate vector of inverse dynamics joint torques for
% S5PRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:49
% EndTime: 2019-12-05 16:38:00
% DurationCPUTime: 5.46s
% Computational Cost: add. (2155->445), mult. (5231->643), div. (0->0), fcn. (4157->12), ass. (0->182)
t368 = sin(qJ(3));
t429 = qJD(2) * qJD(3);
t413 = t368 * t429;
t371 = cos(qJ(3));
t426 = qJDD(2) * t371;
t385 = t413 - t426;
t397 = pkin(3) * t368 - qJ(4) * t371;
t317 = qJD(3) * t397 - qJD(4) * t368;
t365 = cos(pkin(10));
t372 = cos(qJ(2));
t364 = sin(pkin(5));
t441 = qJD(1) * t364;
t420 = t372 * t441;
t369 = sin(qJ(2));
t421 = t369 * t441;
t363 = sin(pkin(10));
t459 = t363 * t371;
t474 = t420 * t459 + (t317 - t421) * t365;
t453 = t365 * t372;
t308 = (t363 * t369 + t371 * t453) * t364;
t473 = qJD(1) * t308 - t363 * t317;
t337 = qJD(2) * pkin(7) + t421;
t366 = cos(pkin(5));
t452 = t366 * t371;
t472 = qJD(1) * t452 - t368 * t337;
t471 = -qJDD(3) * pkin(3) + qJDD(4);
t373 = qJD(3) ^ 2;
t430 = qJD(1) * qJD(2);
t414 = t369 * t430;
t457 = t364 * t372;
t396 = -qJDD(1) * t457 + t364 * t414;
t465 = sin(pkin(9));
t406 = t465 * t369;
t466 = cos(pkin(9));
t407 = t466 * t372;
t318 = -t366 * t407 + t406;
t468 = g(2) * t318;
t405 = t465 * t372;
t408 = t466 * t369;
t320 = t366 * t405 + t408;
t469 = g(1) * t320;
t402 = t468 + t469;
t470 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t373 + t364 * (-g(3) * t372 + t414) - t396 + t402;
t467 = qJD(2) * pkin(2);
t436 = qJD(3) * t363;
t438 = qJD(2) * t368;
t334 = t365 * t438 + t436;
t370 = cos(qJ(5));
t461 = t334 * t370;
t392 = pkin(3) * t371 + qJ(4) * t368 + pkin(2);
t460 = t392 * t365;
t458 = t364 * t369;
t367 = sin(qJ(5));
t456 = t365 * t367;
t455 = t365 * t370;
t454 = t365 * t371;
t359 = t365 * qJDD(3);
t412 = t371 * t429;
t427 = qJDD(2) * t368;
t386 = t412 + t427;
t306 = t363 * t386 - t359;
t300 = qJDD(5) + t306;
t451 = t367 * t300;
t450 = t367 * t368;
t449 = t368 * t370;
t448 = t370 * t300;
t447 = t370 * t371;
t446 = qJDD(1) - g(3);
t312 = qJDD(2) * pkin(7) + (qJDD(1) * t369 + t372 * t430) * t364;
t428 = qJDD(1) * t366;
t411 = t368 * t428;
t251 = t411 + qJDD(3) * qJ(4) + t312 * t371 + (qJD(4) + t472) * qJD(3);
t261 = qJD(2) * t317 - qJDD(2) * t392 + t396;
t245 = t365 * t251 + t363 * t261;
t422 = pkin(7) * t363 + pkin(4);
t435 = qJD(3) * t368;
t445 = -t422 * t435 - t474;
t440 = qJD(1) * t368;
t352 = t366 * t440;
t304 = t371 * t337 + t352;
t291 = qJD(3) * qJ(4) + t304;
t305 = -qJD(2) * t392 - t420;
t257 = t365 * t291 + t363 * t305;
t424 = pkin(7) * t435;
t444 = t363 * t424 + t474;
t443 = -t365 * t424 - t473;
t336 = t397 * qJD(2);
t269 = t363 * t336 + t365 * t472;
t310 = pkin(7) * t454 - t363 * t392;
t361 = t368 ^ 2;
t442 = -t371 ^ 2 + t361;
t439 = qJD(2) * t364;
t437 = qJD(2) * t371;
t434 = qJD(3) * t371;
t433 = qJD(5) * t367;
t432 = qJD(5) * t370;
t360 = t365 * qJD(3);
t332 = t363 * t438 - t360;
t431 = -qJD(5) - t332;
t425 = qJDD(3) * t363;
t423 = t368 * t457;
t418 = t369 * t439;
t417 = t372 * t439;
t416 = t367 * t437;
t415 = qJ(4) * t426;
t410 = t364 * t466;
t409 = t364 * t465;
t241 = pkin(8) * t385 + t245;
t382 = qJD(3) * t352 + t368 * t312 + t337 * t434 - t371 * t428;
t254 = t382 + t471;
t307 = t365 * t386 + t425;
t247 = pkin(4) * t306 - pkin(8) * t307 + t254;
t404 = -t367 * t241 + t370 * t247;
t403 = pkin(4) * t363 - pkin(8) * t365;
t319 = t366 * t408 + t405;
t321 = -t366 * t406 + t407;
t401 = g(1) * t321 + g(2) * t319;
t400 = -qJD(5) * t416 + t307 * t367 - t370 * t413;
t399 = qJ(4) * qJD(5) * t365 + t403 * t437 + t304;
t393 = pkin(7) + t403;
t314 = t393 * t368;
t398 = -qJD(5) * t314 - (-pkin(7) * t365 + pkin(8)) * t435 + t473;
t395 = t370 * t241 + t367 * t247;
t244 = -t251 * t363 + t261 * t365;
t253 = -pkin(8) * t437 + t257;
t288 = -qJD(3) * pkin(3) + qJD(4) - t472;
t255 = pkin(4) * t332 - pkin(8) * t334 + t288;
t243 = t253 * t370 + t255 * t367;
t394 = t253 * t367 - t255 * t370;
t324 = t366 * t368 + t371 * t458;
t281 = t324 * t365 - t363 * t457;
t323 = t368 * t458 - t452;
t263 = t281 * t370 + t323 * t367;
t262 = -t281 * t367 + t323 * t370;
t256 = -t291 * t363 + t305 * t365;
t268 = t336 * t365 - t363 * t472;
t391 = t370 * t307 + t385 * t367;
t389 = t365 * t447 + t450;
t326 = t365 * t450 + t447;
t295 = t334 * t367 + t370 * t437;
t339 = -pkin(4) * t365 - pkin(8) * t363 - pkin(3);
t387 = pkin(8) * t438 - qJD(4) * t365 - qJD(5) * t339 + t269;
t282 = t319 * t368 + t371 * t410;
t284 = t321 * t368 - t371 * t409;
t384 = g(1) * t284 + g(2) * t282 + g(3) * t323;
t283 = t319 * t371 - t368 * t410;
t285 = t321 * t371 + t368 * t409;
t383 = g(1) * t285 + g(2) * t283 + g(3) * t324;
t381 = g(3) * t457 - t402;
t380 = -t254 + t384;
t379 = -qJ(4) * t435 + (qJD(4) - t288) * t371;
t294 = -pkin(8) * t371 + t310;
t378 = qJD(5) * t294 + t368 * t420 - t393 * t434;
t338 = -t420 - t467;
t376 = -pkin(7) * qJDD(3) + (t338 + t420 - t467) * qJD(3);
t375 = -t382 + t384;
t374 = qJD(2) ^ 2;
t327 = t365 * t449 - t367 * t371;
t316 = t389 * qJD(2);
t315 = t365 * t416 - t370 * t438;
t309 = -pkin(7) * t459 - t460;
t297 = -t416 + t461;
t293 = t371 * t422 + t460;
t287 = qJD(3) * t324 + t368 * t417;
t286 = -qJD(3) * t323 + t371 * t417;
t280 = t324 * t363 + t364 * t453;
t274 = -t370 * t435 - t371 * t433 + (t367 * t434 + t368 * t432) * t365;
t273 = qJD(3) * t389 - qJD(5) * t326;
t271 = -t320 * t454 + t321 * t363;
t270 = -t318 * t454 + t319 * t363;
t267 = t286 * t365 + t363 * t418;
t266 = t286 * t363 - t365 * t418;
t264 = -pkin(4) * t438 - t268;
t260 = t285 * t365 + t320 * t363;
t259 = t283 * t365 + t318 * t363;
t252 = pkin(4) * t437 - t256;
t249 = (qJD(5) * t334 + t426) * t370 + t400;
t248 = -qJD(5) * t295 + t391;
t240 = -pkin(4) * t385 - t244;
t239 = -t243 * qJD(5) + t404;
t238 = -t394 * qJD(5) + t395;
t1 = [t446 * MDP(1) + (-qJD(3) * t287 - qJDD(3) * t323) * MDP(10) + (-qJD(3) * t286 - qJDD(3) * t324) * MDP(11) + (t280 * t426 + t287 * t332 + t306 * t323) * MDP(12) + (t281 * t426 + t287 * t334 + t307 * t323) * MDP(13) + (t266 * t334 - t267 * t332 + t280 * t307 - t281 * t306) * MDP(14) + (-t244 * t280 + t245 * t281 + t254 * t323 - t256 * t266 + t257 * t267 + t287 * t288 - g(3)) * MDP(15) + (-(-qJD(5) * t263 - t267 * t367 + t287 * t370) * t431 + t262 * t300 + t266 * t295 + t280 * t249) * MDP(21) + ((qJD(5) * t262 + t267 * t370 + t287 * t367) * t431 - t263 * t300 + t266 * t297 + t280 * t248) * MDP(22) + ((t266 * t371 - t280 * t435) * MDP(12) + (t267 * t371 - t281 * t435) * MDP(13)) * qJD(2) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t371 + MDP(11) * t368 - MDP(3)) * t374) * t369 + (-MDP(10) * t385 - MDP(11) * t386 + qJDD(2) * MDP(3) - t374 * MDP(4)) * t372) * t364; qJDD(2) * MDP(2) + (t446 * t457 + t402) * MDP(3) + (-t446 * t458 + t401) * MDP(4) + (qJDD(2) * t361 + 0.2e1 * t368 * t412) * MDP(5) + 0.2e1 * (t368 * t426 - t429 * t442) * MDP(6) + (qJDD(3) * t368 + t371 * t373) * MDP(7) + (qJDD(3) * t371 - t368 * t373) * MDP(8) + (t376 * t368 + t470 * t371) * MDP(10) + (-t470 * t368 + t376 * t371) * MDP(11) + (-g(1) * t271 - g(2) * t270 - g(3) * t308 + (-t332 * t420 + pkin(7) * t306 + t254 * t363 + (qJD(2) * t309 + t256) * qJD(3)) * t368 + (-qJDD(2) * t309 - t244 + (pkin(7) * t332 + t288 * t363) * qJD(3) - t444 * qJD(2)) * t371) * MDP(12) + ((-g(3) * t458 - t401) * t365 + (-t334 * t420 + pkin(7) * t307 + t254 * t365 + (-qJD(2) * t310 - t257) * qJD(3)) * t368 + (t310 * qJDD(2) + t245 + (pkin(7) * t334 + t288 * t365) * qJD(3) + t443 * qJD(2) + t381 * t363) * t371) * MDP(13) + (-t306 * t310 - t307 * t309 - t444 * t334 - t443 * t332 + (-t256 * t365 - t257 * t363) * t434 + (-t244 * t365 - t245 * t363 - t381) * t368) * MDP(14) + (t244 * t309 + t245 * t310 + t443 * t257 + t444 * t256 + t392 * t469 + t392 * t468 + (t254 * t368 + t288 * t434 - t401) * pkin(7) + (-g(3) * pkin(7) * t369 + (-g(3) * t392 - t288 * t440) * t372) * t364) * MDP(15) + (t248 * t327 + t273 * t297) * MDP(16) + (-t248 * t326 - t249 * t327 - t273 * t295 - t274 * t297) * MDP(17) + (-t273 * t431 + t300 * t327 + (t248 * t368 + t297 * t434) * t363) * MDP(18) + (t274 * t431 - t300 * t326 + (-t249 * t368 - t295 * t434) * t363) * MDP(19) + (t300 * t368 - t431 * t434) * t363 * MDP(20) + ((-t294 * t367 + t314 * t370) * t300 + t293 * t249 + t240 * t326 + t252 * t274 - g(1) * (t271 * t370 - t320 * t450) - g(2) * (t270 * t370 - t318 * t450) - g(3) * (t308 * t370 + t367 * t423) + (t239 * t368 - t394 * t434) * t363 - (t367 * t398 - t370 * t378) * t431 + t445 * t295) * MDP(21) + (-(t294 * t370 + t314 * t367) * t300 + t293 * t248 + t240 * t327 + t252 * t273 - g(1) * (-t271 * t367 - t320 * t449) - g(2) * (-t270 * t367 - t318 * t449) - g(3) * (-t308 * t367 + t370 * t423) + (-t238 * t368 - t243 * t434) * t363 - (t367 * t378 + t370 * t398) * t431 + t445 * t297) * MDP(22); MDP(7) * t427 + MDP(8) * t426 + qJDD(3) * MDP(9) + (qJD(3) * t304 - t338 * t438 + t375) * MDP(10) + (-t411 + (-qJD(2) * t338 - t312) * t371 + t383) * MDP(11) + (t363 * t415 - pkin(3) * t306 - t304 * t332 + t380 * t365 + (-t256 * t368 + t268 * t371 + t363 * t379) * qJD(2)) * MDP(12) + (t365 * t415 - pkin(3) * t307 - t304 * t334 - t380 * t363 + (t257 * t368 - t269 * t371 + t365 * t379) * qJD(2)) * MDP(13) + (t268 * t334 + t269 * t332 + (-qJ(4) * t306 - qJD(4) * t332 + t256 * t437 + t245) * t365 + (qJ(4) * t307 + qJD(4) * t334 + t257 * t437 - t244) * t363 - t383) * MDP(14) + (-t256 * t268 - t257 * t269 - t288 * t304 + (-t256 * t363 + t257 * t365) * qJD(4) + t380 * pkin(3) + (-t244 * t363 + t245 * t365 - t383) * qJ(4)) * MDP(15) + (t248 * t363 * t370 + (-t363 * t433 - t316) * t297) * MDP(16) + (t295 * t316 + t297 * t315 + (-t248 * t367 - t249 * t370 + (t295 * t367 - t297 * t370) * qJD(5)) * t363) * MDP(17) + (-t248 * t365 + t316 * t431 + (-t297 * t437 + t431 * t433 + t448) * t363) * MDP(18) + (t249 * t365 - t315 * t431 + (t295 * t437 + t431 * t432 - t451) * t363) * MDP(19) + (t363 * t431 * t437 - t300 * t365) * MDP(20) + ((-qJ(4) * t456 + t339 * t370) * t300 - t239 * t365 - t264 * t295 - t252 * t315 - g(1) * (-t284 * t455 + t285 * t367) - g(2) * (-t282 * t455 + t283 * t367) - g(3) * (-t323 * t455 + t324 * t367) - (t367 * t387 - t370 * t399) * t431 + (qJ(4) * t249 + qJD(4) * t295 + t240 * t367 + t252 * t432 + t394 * t437) * t363) * MDP(21) + (-(qJ(4) * t455 + t339 * t367) * t300 + t238 * t365 - t264 * t297 - t252 * t316 - g(1) * (t284 * t456 + t285 * t370) - g(2) * (t282 * t456 + t283 * t370) - g(3) * (t323 * t456 + t324 * t370) - (t367 * t399 + t370 * t387) * t431 + (qJ(4) * t248 + qJD(4) * t297 + t240 * t370 + t243 * t437 - t252 * t433) * t363) * MDP(22) + (-MDP(5) * t368 * t371 + MDP(6) * t442) * t374; (t363 * t427 - t359) * MDP(12) + (t365 * t427 + t425) * MDP(13) + (-t332 ^ 2 - t334 ^ 2) * MDP(14) + (t256 * t334 + t257 * t332 - t375 + t471) * MDP(15) + (-t295 * t334 + t448) * MDP(21) + (-t297 * t334 - t451) * MDP(22) + ((-t334 + t436) * MDP(12) + (t332 + t360) * MDP(13)) * t437 - (MDP(21) * t367 + MDP(22) * t370) * t431 ^ 2; t297 * t295 * MDP(16) + (-t295 ^ 2 + t297 ^ 2) * MDP(17) + (-t295 * t431 + t391) * MDP(18) + (-t297 * t431 - t370 * t426 - t400) * MDP(19) + t300 * MDP(20) + (-t243 * t431 - t252 * t297 - g(1) * (-t260 * t367 + t284 * t370) - g(2) * (-t259 * t367 + t282 * t370) - g(3) * t262 + t404) * MDP(21) + (t394 * t431 + t252 * t295 - g(1) * (-t260 * t370 - t284 * t367) - g(2) * (-t259 * t370 - t282 * t367) + g(3) * t263 - t395) * MDP(22) + (-MDP(18) * t295 - MDP(19) * t461 - MDP(21) * t243 + MDP(22) * t394) * qJD(5);];
tau = t1;
