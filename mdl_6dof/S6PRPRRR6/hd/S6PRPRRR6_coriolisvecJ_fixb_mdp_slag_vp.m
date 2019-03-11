% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:57
% EndTime: 2019-03-08 20:49:04
% DurationCPUTime: 3.94s
% Computational Cost: add. (2140->395), mult. (4928->559), div. (0->0), fcn. (3525->10), ass. (0->184)
t369 = -pkin(2) - pkin(8);
t368 = cos(qJ(2));
t359 = sin(pkin(6));
t438 = qJD(1) * t359;
t414 = t368 * t438;
t389 = qJD(3) - t414;
t317 = qJD(2) * t369 + t389;
t360 = cos(pkin(6));
t367 = cos(qJ(4));
t363 = sin(qJ(4));
t437 = qJD(1) * t363;
t294 = t317 * t367 - t360 * t437;
t281 = -qJD(4) * pkin(4) - t294;
t362 = sin(qJ(5));
t366 = cos(qJ(5));
t423 = t366 * qJD(4);
t433 = qJD(2) * t367;
t327 = t362 * t433 - t423;
t268 = pkin(5) * t327 + t281;
t365 = cos(qJ(6));
t410 = t366 * t433;
t431 = qJD(4) * t362;
t329 = t410 + t431;
t361 = sin(qJ(6));
t463 = t329 * t361;
t274 = t365 * t327 + t463;
t487 = t268 * t274;
t434 = qJD(2) * t363;
t352 = qJD(5) + t434;
t348 = qJD(6) + t352;
t486 = t274 * t348;
t385 = t327 * t361 - t365 * t329;
t485 = t348 * t385;
t330 = t361 * t362 - t365 * t366;
t477 = qJD(5) + qJD(6);
t484 = t477 * t330;
t331 = t361 * t366 + t362 * t365;
t372 = t477 * t331;
t483 = MDP(13) * t363 + MDP(14) * t367;
t426 = qJD(5) * t367;
t407 = t362 * t426;
t377 = -t363 * t423 - t407;
t298 = qJD(2) * t377 + qJD(5) * t423;
t458 = t360 * t367;
t350 = qJD(1) * t458;
t295 = t363 * t317 + t350;
t282 = qJD(4) * pkin(9) + t295;
t336 = pkin(4) * t363 - pkin(9) * t367 + qJ(3);
t364 = sin(qJ(2));
t415 = t364 * t438;
t303 = qJD(2) * t336 + t415;
t467 = t303 * t362;
t262 = t282 * t366 + t467;
t435 = qJD(2) * t359;
t413 = t364 * t435;
t429 = qJD(4) * t367;
t269 = t317 * t429 + (-qJD(4) * t360 + t413) * t437;
t393 = pkin(4) * t367 + pkin(9) * t363;
t323 = qJD(4) * t393 + qJD(3);
t296 = (t323 + t414) * qJD(2);
t401 = -t269 * t362 + t366 * t296;
t373 = -qJD(5) * t262 + t401;
t422 = qJD(2) * qJD(4);
t405 = t367 * t422;
t246 = pkin(5) * t405 - pkin(10) * t298 + t373;
t430 = qJD(4) * t363;
t409 = t362 * t430;
t299 = -qJD(2) * t409 + t329 * qJD(5);
t427 = qJD(5) * t366;
t419 = -t366 * t269 - t362 * t296 - t303 * t427;
t428 = qJD(5) * t362;
t380 = -t282 * t428 - t419;
t249 = -pkin(10) * t299 + t380;
t403 = t365 * t246 - t361 * t249;
t482 = t268 * t385 + t403;
t481 = MDP(26) * t405 + (-t274 ^ 2 + t385 ^ 2) * MDP(23) - t274 * MDP(22) * t385;
t358 = t367 ^ 2;
t480 = MDP(9) * (t363 ^ 2 - t358);
t457 = t362 * t364;
t479 = -(-t363 * t457 + t366 * t368) * t438 + t366 * t323;
t452 = t364 * t366;
t478 = (t362 * t368 + t363 * t452) * t438 - t367 * t369 * t423 - t362 * t323 - t336 * t427;
t399 = t298 * t361 + t365 * t299;
t251 = -qJD(6) * t385 + t399;
t476 = pkin(9) + pkin(10);
t475 = qJD(2) * pkin(2);
t261 = -t282 * t362 + t366 * t303;
t254 = -pkin(10) * t329 + t261;
t252 = pkin(5) * t352 + t254;
t474 = t252 * t365;
t255 = -pkin(10) * t327 + t262;
t473 = t255 * t365;
t395 = t367 * t413;
t445 = -qJD(4) * t350 - t317 * t430;
t270 = -qJD(1) * t395 - t445;
t472 = t270 * t362;
t471 = t270 * t366;
t470 = t281 * t362;
t469 = t281 * t366;
t468 = t298 * t362;
t465 = t327 * t352;
t464 = t329 * t352;
t462 = t352 * t362;
t461 = t352 * t363;
t460 = t352 * t366;
t459 = t359 * t368;
t456 = t362 * t367;
t455 = t362 * t369;
t454 = t363 * t366;
t453 = t363 * t369;
t451 = t366 * t367;
t347 = t366 * t453;
t404 = pkin(5) - t455;
t421 = pkin(10) * t454;
t450 = (-t347 + (pkin(10) * t367 - t336) * t362) * qJD(5) + (t367 * t404 + t421) * qJD(4) + t479;
t406 = t366 * t426;
t376 = -t406 + t409;
t420 = t362 * t453;
t449 = -pkin(10) * t376 + qJD(5) * t420 + t478;
t332 = t393 * qJD(2);
t448 = t366 * t294 + t362 * t332;
t381 = t363 * t330;
t447 = -qJD(2) * t381 - t484;
t379 = t331 * qJD(2);
t446 = t363 * t379 + t372;
t444 = t362 * t336 + t347;
t436 = qJD(2) * qJ(3);
t432 = qJD(4) * t348;
t425 = qJD(6) * t361;
t424 = qJD(6) * t365;
t418 = t365 * t298 - t361 * t299 - t327 * t424;
t416 = qJD(5) * t476;
t412 = t368 * t435;
t411 = t362 * t434;
t253 = t255 * t425;
t402 = t361 * t246 - t253;
t400 = t352 * t369 + t282;
t398 = -t294 * t362 + t366 * t332;
t397 = qJD(6) * t252 + t249;
t396 = qJD(5) * t363 + qJD(2);
t346 = t363 * t405;
t394 = -t295 + (t411 + t428) * pkin(5);
t335 = t415 + t436;
t392 = -t335 + t415;
t343 = t476 * t366;
t391 = qJD(6) * t343 + (pkin(5) * t367 + t421) * qJD(2) + t398 + t366 * t416;
t342 = t476 * t362;
t390 = pkin(10) * t411 + qJD(6) * t342 + t362 * t416 + t448;
t248 = t252 * t361 + t473;
t322 = t366 * t336;
t272 = -pkin(10) * t451 + t363 * t404 + t322;
t285 = -pkin(10) * t456 + t444;
t388 = t272 * t361 + t285 * t365;
t316 = -t363 * t459 + t458;
t288 = -t316 * t362 + t359 * t452;
t289 = t316 * t366 + t359 * t457;
t387 = t288 * t365 - t289 * t361;
t386 = t288 * t361 + t289 * t365;
t384 = qJD(2) * t358 - t461;
t383 = -pkin(9) * t429 + t281 * t363;
t315 = t360 * t363 + t367 * t459;
t382 = t392 * qJD(2);
t250 = -t329 * t425 + t418;
t378 = t330 * qJD(2);
t375 = t392 - t436;
t324 = (qJD(3) + t414) * qJD(2);
t370 = qJD(4) ^ 2;
t374 = qJD(2) * t389 - t369 * t370 + t324;
t371 = qJD(2) ^ 2;
t355 = -pkin(5) * t366 - pkin(4);
t326 = t389 - t475;
t325 = (pkin(5) * t362 - t369) * t367;
t309 = t330 * t367;
t308 = t331 * t367;
t300 = -pkin(5) * t376 + t369 * t430;
t287 = qJD(4) * t316 - t395;
t286 = -qJD(4) * t315 + t363 * t413;
t264 = -t425 * t456 + (t477 * t451 - t409) * t365 + t377 * t361;
t263 = qJD(4) * t381 - t367 * t372;
t260 = pkin(5) * t299 + t270;
t257 = qJD(5) * t288 + t286 * t366 + t362 * t412;
t256 = -qJD(5) * t289 - t286 * t362 + t366 * t412;
t247 = -t255 * t361 + t474;
t1 = [(t256 * t352 + t287 * t327 + t299 * t315) * MDP(20) + (-t257 * t352 + t287 * t329 + t298 * t315) * MDP(21) + ((-qJD(6) * t386 + t256 * t365 - t257 * t361) * t348 + t287 * t274 + t315 * t251) * MDP(27) + (-(qJD(6) * t387 + t256 * t361 + t257 * t365) * t348 - t287 * t385 + t315 * t250) * MDP(28) + (-t287 * MDP(13) - t286 * MDP(14) + (t288 * MDP(20) - t289 * MDP(21) + MDP(27) * t387 - MDP(28) * t386) * t433) * qJD(4) + ((qJD(2) * t335 * MDP(7) + (-MDP(4) + MDP(6) + t483) * t371) * t368 + (t324 * MDP(7) + (-MDP(3) + MDP(5)) * t371 + ((MDP(13) * t367 - MDP(14) * t363) * qJD(4) + (t326 - t414) * MDP(7)) * qJD(2)) * t364) * t359; 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (qJ(3) * t324 + qJD(3) * t335 + (-t335 * t368 + (-t326 - t475) * t364) * t438) * MDP(7) - 0.2e1 * MDP(8) * t346 + 0.2e1 * t422 * t480 + (t363 * t374 - t375 * t429) * MDP(13) + (t367 * t374 + t375 * t430) * MDP(14) + (t298 * t451 + t329 * t377) * MDP(15) + ((t327 * t366 + t329 * t362) * t430 + (-t468 - t299 * t366 + (t327 * t362 - t329 * t366) * qJD(5)) * t367) * MDP(16) + (-t352 * t407 + t298 * t363 + (t329 * t367 + t366 * t384) * qJD(4)) * MDP(17) + (-t352 * t406 - t299 * t363 + (-t327 * t367 - t362 * t384) * qJD(4)) * MDP(18) + (t352 * t429 + t346) * MDP(19) + ((-t336 * t428 + t479) * t352 + ((t327 * t369 - t470) * qJD(4) + (-t366 * t400 - t467) * qJD(5) + t401) * t363 + (t327 * t415 + t281 * t427 + t472 - t369 * t299 + (-t352 * t455 + (t322 - t420) * qJD(2) + t261) * qJD(4)) * t367) * MDP(20) + (t478 * t352 + (t400 * t428 + (t329 * t369 - t469) * qJD(4) + t419) * t363 + (t329 * t415 - t281 * t428 + t471 - t369 * t298 + (-qJD(2) * t444 - t262) * qJD(4)) * t367) * MDP(21) + (-t250 * t309 - t263 * t385) * MDP(22) + (-t250 * t308 + t251 * t309 - t263 * t274 + t264 * t385) * MDP(23) + (t250 * t363 + t263 * t348 + (-qJD(2) * t309 - t385) * t429) * MDP(24) + (-t251 * t363 - t264 * t348 + (-qJD(2) * t308 - t274) * t429) * MDP(25) + (t348 * t429 + t346) * MDP(26) + (t403 * t363 + t300 * t274 + t325 * t251 + t260 * t308 + t268 * t264 + (t361 * t449 + t365 * t450) * t348 + (-t248 * t363 - t348 * t388) * qJD(6) + (t274 * t415 + ((t272 * t365 - t285 * t361) * qJD(2) + t247) * qJD(4)) * t367) * MDP(27) + (-(t397 * t365 + t402) * t363 - t300 * t385 + t325 * t250 - t260 * t309 + t268 * t263 + ((-qJD(6) * t272 + t449) * t365 + (qJD(6) * t285 - t450) * t361) * t348 + (-t385 * t415 + (-qJD(2) * t388 - t248) * qJD(4)) * t367) * MDP(28) + (-MDP(10) * t363 - MDP(11) * t367) * t370; -t371 * MDP(6) + MDP(7) * t382 + (-t299 * t367 - t396 * t460 + (t327 * t363 + (-t352 - t434) * t456) * qJD(4)) * MDP(20) + (-t298 * t367 + t396 * t462 + (-t352 * t451 + (t329 - t410) * t363) * qJD(4)) * MDP(21) + (t348 * t378 + (-t331 * t432 - t251) * t367 + ((-t331 * t433 + t274) * qJD(4) + t348 * t484) * t363) * MDP(27) + (t348 * t379 + (t330 * t432 - t250) * t367 + (t372 * t348 + (t367 * t378 - t385) * qJD(4)) * t363) * MDP(28) + t483 * (-t370 - t371); (qJD(4) * t295 + t367 * t382 + t445) * MDP(13) - t392 * t434 * MDP(14) + (t329 * t460 + t468) * MDP(15) + ((t298 - t465) * t366 + (-t299 - t464) * t362) * MDP(16) + (t352 * t427 + (t352 * t454 + (-t329 + t431) * t367) * qJD(2)) * MDP(17) + (-t352 * t428 + (-t362 * t461 + (t327 + t423) * t367) * qJD(2)) * MDP(18) + (-pkin(4) * t299 - t471 - t398 * t352 - t295 * t327 + (-pkin(9) * t460 + t470) * qJD(5) + (-t261 * t367 + t362 * t383) * qJD(2)) * MDP(20) + (-pkin(4) * t298 + t472 + t448 * t352 - t295 * t329 + (pkin(9) * t462 + t469) * qJD(5) + (t262 * t367 + t366 * t383) * qJD(2)) * MDP(21) + (t250 * t331 - t385 * t447) * MDP(22) + (-t250 * t330 - t251 * t331 - t274 * t447 + t385 * t446) * MDP(23) + (t355 * t251 + t260 * t330 + t446 * t268 + t394 * t274) * MDP(27) + (t355 * t250 + t260 * t331 + t447 * t268 - t385 * t394) * MDP(28) + (t447 * MDP(24) - t446 * MDP(25) + (t361 * t390 - t365 * t391) * MDP(27) + (t361 * t391 + t365 * t390) * MDP(28)) * t348 + (-t352 * MDP(19) + (qJD(4) * t331 + t385) * MDP(24) + (-qJD(4) * t330 + t274) * MDP(25) - t348 * MDP(26) + ((-t342 * t365 - t343 * t361) * qJD(4) - t247) * MDP(27) + (-(-t342 * t361 + t343 * t365) * qJD(4) + t248) * MDP(28)) * t433 + (MDP(8) * t363 * t367 - t480) * t371; t329 * t327 * MDP(15) + (-t327 ^ 2 + t329 ^ 2) * MDP(16) + (t298 + t465) * MDP(17) + (-t299 + t464) * MDP(18) + MDP(19) * t405 + (t262 * t352 - t281 * t329 + t373) * MDP(20) + (t261 * t352 + t281 * t327 - t380) * MDP(21) + (t250 + t486) * MDP(24) + (-t251 - t485) * MDP(25) + (-(-t254 * t361 - t473) * t348 - t248 * qJD(6) + (-t274 * t329 - t348 * t425 + t365 * t405) * pkin(5) + t482) * MDP(27) + (t487 + t253 + (-t255 * t348 - t246) * t361 + (t254 * t348 - t397) * t365 + (t329 * t385 - t348 * t424 - t361 * t405) * pkin(5)) * MDP(28) + t481; (t418 + t486) * MDP(24) + (-t399 - t485) * MDP(25) + (t248 * t348 + t482) * MDP(27) + (t247 * t348 - t365 * t249 - t402 + t487) * MDP(28) + (-MDP(24) * t463 + MDP(25) * t385 - MDP(27) * t248 - MDP(28) * t474) * qJD(6) + t481;];
tauc  = t1;
