% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:36
% EndTime: 2019-03-09 04:48:44
% DurationCPUTime: 4.32s
% Computational Cost: add. (4424->416), mult. (9454->563), div. (0->0), fcn. (5837->6), ass. (0->177)
t398 = sin(qJ(4));
t400 = cos(qJ(4));
t446 = t400 * qJD(3);
t401 = cos(qJ(3));
t456 = qJD(1) * t401;
t364 = t398 * t456 - t446;
t455 = qJD(3) * t398;
t366 = t400 * t456 + t455;
t396 = sin(pkin(9));
t397 = cos(pkin(9));
t317 = t397 * t364 + t366 * t396;
t399 = sin(qJ(3));
t457 = qJD(1) * t399;
t387 = qJD(4) + t457;
t499 = t317 * t387;
t362 = t396 * t400 + t397 * t398;
t348 = t362 * qJD(1);
t479 = t397 * t400;
t480 = t396 * t398;
t361 = -t479 + t480;
t451 = qJD(4) * t399;
t453 = qJD(3) * t401;
t493 = -t361 * t453 - t362 * t451 - t348;
t347 = t362 * qJD(4);
t463 = t399 * t348 + t347;
t427 = t398 * t457;
t428 = qJD(1) * t479;
t450 = qJD(4) * t400;
t452 = qJD(4) * t398;
t492 = t396 * t452 - t397 * t450;
t462 = -t396 * t427 + t399 * t428 - t492;
t442 = 0.2e1 * qJD(1);
t498 = qJD(2) * t442;
t491 = MDP(21) + MDP(24);
t415 = -t364 * t396 + t397 * t366;
t497 = t415 ^ 2;
t395 = t401 ^ 2;
t496 = MDP(8) * (t399 ^ 2 - t395);
t402 = -pkin(1) - pkin(7);
t382 = qJD(1) * t402 + qJD(2);
t482 = t382 * t401;
t355 = -qJD(3) * pkin(3) - t482;
t322 = pkin(4) * t364 + qJD(5) + t355;
t271 = pkin(5) * t317 - qJ(6) * t415 + t322;
t495 = t271 * t415;
t494 = MDP(6) * qJ(2) + MDP(5);
t433 = t401 * t446;
t435 = t398 * t453;
t465 = qJD(1) * t480 - t396 * t433 - t397 * t435 + t399 * t492 - t428;
t418 = pkin(3) * t401 + pkin(8) * t399;
t368 = t418 * qJD(1);
t353 = t400 * t368;
t475 = t399 * t400;
t477 = t398 * t401;
t301 = -t382 * t477 + t353 + (pkin(4) * t401 + qJ(5) * t475) * qJD(1);
t473 = t400 * t401;
t461 = t398 * t368 + t382 * t473;
t310 = qJ(5) * t427 + t461;
t490 = -qJ(5) - pkin(8);
t423 = qJD(4) * t490;
t448 = qJD(5) * t400;
t341 = t398 * t423 + t448;
t406 = -qJD(5) * t398 + t400 * t423;
t464 = (-t301 + t406) * t397 + (t310 - t341) * t396;
t370 = pkin(3) * t399 - pkin(8) * t401 + qJ(2);
t474 = t399 * t402;
t460 = t398 * t370 + t400 * t474;
t351 = t370 * qJD(1);
t369 = t399 * t382;
t354 = qJD(3) * pkin(8) + t369;
t313 = t351 * t398 + t354 * t400;
t360 = qJD(3) * t418 + qJD(2);
t340 = t360 * qJD(1);
t405 = -t313 * qJD(4) + t400 * t340;
t299 = -qJ(5) * t364 + t313;
t295 = t397 * t299;
t312 = t400 * t351 - t354 * t398;
t298 = -qJ(5) * t366 + t312;
t269 = t298 * t396 + t295;
t488 = t269 * t415;
t487 = t299 * t396;
t443 = qJD(3) * qJD(4);
t392 = t400 * t443;
t449 = qJD(4) * t401;
t431 = t398 * t449;
t434 = t399 * t446;
t407 = -t431 - t434;
t327 = qJD(1) * t407 + t392;
t486 = t327 * t398;
t485 = t355 * t398;
t484 = t364 * t387;
t483 = t366 * t387;
t481 = t387 * t400;
t478 = t398 * t387;
t476 = t398 * t402;
t472 = t401 * t402;
t404 = qJD(1) ^ 2;
t471 = t401 * t404;
t403 = qJD(3) ^ 2;
t470 = t402 * t403;
t263 = -qJ(5) * t327 - qJD(5) * t366 + (pkin(4) * qJD(1) - t382 * t398) * t453 + t405;
t379 = qJD(3) * t427;
t328 = t366 * qJD(4) - t379;
t409 = t398 * t340 + t351 * t450 - t354 * t452 + t382 * t433;
t268 = -qJ(5) * t328 - qJD(5) * t364 + t409;
t469 = -t397 * t263 + t396 * t268;
t252 = t396 * t263 + t397 * t268;
t275 = t396 * t301 + t397 * t310;
t272 = qJ(6) * t456 + t275;
t308 = t397 * t341 + t396 * t406;
t468 = -t272 + t308;
t467 = -pkin(5) * t456 + t464;
t343 = t400 * t360;
t424 = pkin(4) - t476;
t279 = qJ(5) * t434 + t343 - t460 * qJD(4) + (qJ(5) * t452 + qJD(3) * t424 - t448) * t401;
t430 = t400 * t449;
t438 = t398 * t360 + t370 * t450 + t402 * t433;
t282 = -qJ(5) * t430 + (-qJD(5) * t401 + (qJ(5) * qJD(3) - qJD(4) * t402) * t399) * t398 + t438;
t257 = t396 * t279 + t397 * t282;
t466 = qJD(6) * t362 + t369 + t462 * qJ(6) - t463 * pkin(5) + (-t427 - t452) * pkin(4);
t291 = pkin(4) * t387 + t298;
t267 = t396 * t291 + t295;
t359 = t400 * t370;
t315 = -qJ(5) * t473 + t399 * t424 + t359;
t323 = -qJ(5) * t477 + t460;
t288 = t396 * t315 + t397 * t323;
t458 = -t403 - t404;
t454 = qJD(3) * t399;
t447 = qJD(6) * t387;
t270 = t298 * t397 - t487;
t445 = qJD(6) - t270;
t444 = qJD(1) * qJD(3);
t440 = t398 * t474;
t425 = t401 * t444;
t439 = qJ(6) * t425 + t252;
t437 = -pkin(4) * t400 - pkin(3);
t436 = t398 * t454;
t426 = t490 * t398;
t309 = pkin(4) * t328 + t382 * t454;
t293 = t327 * t396 + t397 * t328;
t422 = pkin(4) * t477 - t472;
t421 = t364 + t446;
t420 = -t366 + t455;
t294 = t327 * t397 - t328 * t396;
t377 = t490 * t400;
t325 = -t377 * t396 - t397 * t426;
t326 = -t397 * t377 + t396 * t426;
t417 = -t326 * t293 + t294 * t325 - t308 * t317;
t256 = t279 * t397 - t282 * t396;
t266 = t291 * t397 - t487;
t287 = t315 * t397 - t323 * t396;
t414 = qJD(1) * t395 - t387 * t399;
t413 = -MDP(19) * t398 - MDP(20) * t400;
t412 = -pkin(8) * t453 + t355 * t399;
t410 = t402 * t454 + (t430 - t436) * pkin(4);
t250 = -pkin(5) * t425 + t469;
t335 = t362 * t399;
t253 = pkin(5) * t293 - qJ(6) * t294 - qJD(6) * t415 + t309;
t391 = -pkin(4) * t397 - pkin(5);
t389 = pkin(4) * t396 + qJ(6);
t338 = -t396 * t477 + t397 * t473;
t337 = t361 * t399;
t336 = t362 * t401;
t314 = pkin(5) * t361 - qJ(6) * t362 + t437;
t306 = t347 * t401 - t396 * t436 + t397 * t434;
t304 = qJD(3) * t335 + t361 * t449;
t297 = pkin(5) * t336 - qJ(6) * t338 + t422;
t283 = -pkin(5) * t399 - t287;
t281 = qJ(6) * t399 + t288;
t276 = pkin(4) * t366 + pkin(5) * t415 + qJ(6) * t317;
t260 = qJ(6) * t387 + t267;
t259 = -pkin(5) * t387 + qJD(6) - t266;
t258 = -pkin(5) * t304 + qJ(6) * t306 - qJD(6) * t338 + t410;
t255 = -pkin(5) * t453 - t256;
t254 = qJ(6) * t453 + qJD(6) * t399 + t257;
t249 = t439 + t447;
t1 = [0.2e1 * t444 * t496 - t403 * t401 * MDP(10) + (-t401 * t470 + (-qJ(2) * t454 + qJD(2) * t401) * t442) * MDP(13) + (t327 * t473 + t407 * t366) * MDP(14) + ((t364 * t400 + t366 * t398) * t454 + (-t486 - t328 * t400 + (t364 * t398 - t366 * t400) * qJD(4)) * t401) * MDP(15) + (-t387 * t431 + (t366 * t401 + t400 * t414) * qJD(3)) * MDP(16) + (-t387 * t430 + (-t364 * t401 - t398 * t414) * qJD(3)) * MDP(17) + (-t328 * t472 + t343 * t387 + (t355 * t473 - t387 * t460) * qJD(4) + (-t387 * t476 + (t359 - t440) * qJD(1) + t312) * t453) * MDP(19) + (-(-qJD(4) * t440 + t438) * t387 + (-t402 * t327 - t355 * t452) * t401 + (-qJD(1) * t460 - t313) * t453) * MDP(20) + (-t252 * t336 - t256 * t415 - t257 * t317 + t266 * t306 + t267 * t304 - t287 * t294 - t288 * t293 + t338 * t469) * MDP(21) + (t252 * t288 + t266 * t256 + t267 * t257 - t287 * t469 + t309 * t422 + t322 * t410) * MDP(22) + (t253 * t336 - t255 * t387 + t258 * t317 - t271 * t304 + t293 * t297) * MDP(23) + (-t249 * t336 + t250 * t338 - t254 * t317 + t255 * t415 - t259 * t306 + t260 * t304 - t281 * t293 + t283 * t294) * MDP(24) + (-t253 * t338 + t254 * t387 - t258 * t415 + t271 * t306 - t294 * t297) * MDP(25) + (t249 * t281 + t250 * t283 + t253 * t297 + t254 * t260 + t255 * t259 + t258 * t271) * MDP(26) + (qJ(2) * t442 * MDP(12) + (t387 + t457) * MDP(18) + (-qJD(1) * t283 - t259) * MDP(23) + (qJD(1) * t281 + t260) * MDP(25)) * t453 + t494 * t498 + (-0.2e1 * MDP(7) * t425 - t403 * MDP(9) + (-t470 + t498) * MDP(12) + t327 * MDP(16) - t328 * MDP(17) + ((t364 * t402 - t485) * qJD(3) + t405) * MDP(19) + (-t409 + (t402 * t366 + (-t355 + t482) * t400) * qJD(3)) * MDP(20) - t250 * MDP(23) + t249 * MDP(25)) * t399; (-t252 * t337 + t266 * t465 + t267 * t493 + t335 * t469) * MDP(22) + (-t249 * t337 + t250 * t335 - t259 * t465 + t260 * t493) * MDP(26) - t494 * t404 + (t458 * MDP(12) + (t364 * MDP(19) + t366 * MDP(20) + MDP(22) * t322 + MDP(23) * t317 - MDP(25) * t415 + MDP(26) * t271) * qJD(3)) * t399 + (t465 * MDP(23) + t493 * MDP(25) + (-MDP(19) * t400 + MDP(20) * t398) * (qJD(1) + t451)) * t387 + (t458 * MDP(13) - t328 * MDP(19) - t327 * MDP(20) - t309 * MDP(22) - t293 * MDP(23) + t294 * MDP(25) - t253 * MDP(26) + (t413 * t387 + (-t335 * MDP(23) - t337 * MDP(25) + t399 * t413) * qJD(1)) * qJD(3)) * t401 + (t337 * t293 + t294 * t335 - t317 * t493 - t415 * t465) * t491; t399 * MDP(7) * t471 - t404 * t496 + (t366 * t481 + t486) * MDP(14) + ((t327 - t484) * t400 + (-t328 - t483) * t398) * MDP(15) + (t387 * t450 + (t387 * t475 + t401 * t420) * qJD(1)) * MDP(16) + (-t387 * t452 + (-t399 * t478 + t401 * t421) * qJD(1)) * MDP(17) - t387 * MDP(18) * t456 + (-pkin(3) * t328 - t353 * t387 + (t387 * t477 - t399 * t421) * t382 + (-pkin(8) * t481 + t485) * qJD(4) + (-t312 * t401 + t398 * t412) * qJD(1)) * MDP(19) + (-pkin(3) * t327 + t461 * t387 + t420 * t369 + (pkin(8) * t478 + t355 * t400) * qJD(4) + (t313 * t401 + t400 * t412) * qJD(1)) * MDP(20) + (-t252 * t361 - t266 * t462 - t267 * t463 + t275 * t317 + t362 * t469 - t415 * t464 + t417) * MDP(21) + (t252 * t326 + t469 * t325 + t309 * t437 + (pkin(4) * t478 - t369) * t322 + (t308 - t275) * t267 + t464 * t266) * MDP(22) + (t253 * t361 + t293 * t314 + t467 * t387 - t466 * t317 + t463 * t271 + (-qJD(3) * t325 + t259) * t456) * MDP(23) + (-t249 * t361 + t250 * t362 + t259 * t462 - t260 * t463 + t272 * t317 - t415 * t467 + t417) * MDP(24) + (-t253 * t362 - t294 * t314 + t468 * t387 + t466 * t415 - t462 * t271 + (qJD(3) * t326 - t260) * t456) * MDP(25) + (t249 * t326 + t250 * t325 + t253 * t314 - t259 * t467 + t260 * t468 - t271 * t466) * MDP(26) + (MDP(13) * t399 * t404 - MDP(12) * t471) * qJ(2); t366 * t364 * MDP(14) + (-t364 ^ 2 + t366 ^ 2) * MDP(15) + (t392 + t484) * MDP(16) + (-t398 * t443 + t379 + t483) * MDP(17) + (t313 * t387 - t355 * t366 - t382 * t435 + t405) * MDP(19) + (t312 * t387 + t355 * t364 - t409) * MDP(20) + (t267 * t415 - t488) * MDP(21) + (t266 * t269 - t267 * t270) * MDP(22) + (t269 * t387 - t469 - t495) * MDP(23) + (t260 * t415 - t293 * t389 + t294 * t391 - t488) * MDP(24) + (-t270 * t387 + t276 * t415 + t439 + 0.2e1 * t447) * MDP(25) + (t249 * t389 + t250 * t391 - t259 * t269 + t260 * t445 - t271 * t276) * MDP(26) + (-MDP(16) * t434 + ((-t398 * MDP(16) - t400 * MDP(17)) * qJD(4) + (MDP(18) + (pkin(5) - t391) * MDP(23) + t389 * MDP(25)) * qJD(3)) * t401) * qJD(1) + ((-t293 * t396 - t294 * t397) * MDP(21) + (t252 * t396 - t322 * t366 - t397 * t469) * MDP(22)) * pkin(4) + ((-t266 + t270) * MDP(21) - t276 * MDP(23) + (t259 - t445) * MDP(24) - t271 * MDP(25)) * t317; (t266 * t415 + t267 * t317 + t309) * MDP(22) + (t387 * t415 + t293) * MDP(23) + (-t294 + t499) * MDP(25) + (-t259 * t415 + t260 * t317 + t253) * MDP(26) + t491 * (-t317 ^ 2 - t497); (t317 * t415 - t425) * MDP(23) + (t294 + t499) * MDP(24) + (-t387 ^ 2 - t497) * MDP(25) + (-t260 * t387 + t250 + t495) * MDP(26);];
tauc  = t1;
