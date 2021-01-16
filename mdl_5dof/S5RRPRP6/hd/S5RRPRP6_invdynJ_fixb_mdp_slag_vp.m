% Calculate vector of inverse dynamics joint torques for
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:30:11
% EndTime: 2021-01-15 20:30:25
% DurationCPUTime: 5.20s
% Computational Cost: add. (3683->448), mult. (8513->543), div. (0->0), fcn. (5985->10), ass. (0->193)
t388 = sin(pkin(8));
t395 = cos(qJ(2));
t392 = sin(qJ(2));
t483 = cos(pkin(8));
t434 = t483 * t392;
t355 = t388 * t395 + t434;
t345 = t355 * qJD(1);
t391 = sin(qJ(4));
t394 = cos(qJ(4));
t329 = qJD(2) * t391 + t345 * t394;
t507 = t329 * t391;
t433 = t483 * t395;
t370 = qJD(1) * t433;
t452 = qJD(1) * t392;
t342 = t388 * t452 - t370;
t336 = qJD(4) + t342;
t506 = t336 * t507;
t385 = qJ(2) + pkin(8);
t379 = sin(t385);
t393 = sin(qJ(1));
t396 = cos(qJ(1));
t426 = g(1) * t396 + g(2) * t393;
t417 = t426 * t379;
t344 = t355 * qJD(2);
t447 = qJDD(1) * t392;
t423 = -qJDD(1) * t433 + t388 * t447;
t318 = qJD(1) * t344 + t423;
t449 = qJD(1) * qJD(2);
t439 = t392 * t449;
t405 = qJDD(1) * t355 - t388 * t439;
t402 = qJD(2) * t370 + t405;
t445 = pkin(2) * t439 + qJDD(3);
t446 = qJDD(1) * t395;
t480 = qJDD(1) * pkin(1);
t273 = -pkin(2) * t446 + t318 * pkin(3) - pkin(7) * t402 + t445 - t480;
t270 = t394 * t273;
t485 = t395 * pkin(2);
t378 = pkin(1) + t485;
t361 = -qJD(1) * t378 + qJD(3);
t288 = pkin(3) * t342 - pkin(7) * t345 + t361;
t390 = -qJ(3) - pkin(6);
t440 = t390 * t392;
t357 = qJD(1) * t440;
t484 = qJD(2) * pkin(2);
t351 = t357 + t484;
t363 = t390 * t395;
t358 = qJD(1) * t363;
t435 = t483 * t358;
t317 = t388 * t351 - t435;
t307 = qJD(2) * pkin(7) + t317;
t272 = t288 * t391 + t307 * t394;
t436 = qJD(2) * t390;
t409 = -qJD(3) * t392 + t395 * t436;
t312 = qJDD(2) * pkin(2) + qJD(1) * t409 + qJDD(1) * t440;
t341 = qJD(3) * t395 + t392 * t436;
t319 = qJD(1) * t341 - qJDD(1) * t363;
t280 = t388 * t312 + t483 * t319;
t278 = qJDD(2) * pkin(7) + t280;
t404 = -qJD(4) * t272 - t278 * t391 + t270;
t448 = qJD(4) * qJD(2);
t451 = qJD(4) * t391;
t412 = t391 * qJDD(2) - t345 * t451 + (t402 + t448) * t394;
t482 = qJ(5) * t412;
t313 = qJDD(4) + t318;
t497 = pkin(4) * t313;
t256 = -qJD(5) * t329 + t404 - t482 + t497;
t327 = -t394 * qJD(2) + t345 * t391;
t266 = -qJ(5) * t327 + t272;
t479 = t266 * t336;
t505 = t256 + t479;
t495 = g(1) * t393;
t504 = -g(2) * t396 + t495;
t503 = t329 * qJD(4);
t380 = cos(t385);
t462 = t394 * t396;
t465 = t391 * t393;
t337 = t380 * t465 + t462;
t463 = t393 * t394;
t464 = t391 * t396;
t339 = -t380 * t464 + t463;
t489 = g(3) * t391;
t502 = -g(1) * t339 + g(2) * t337 + t379 * t489;
t491 = g(3) * t379;
t501 = t426 * t380 + t491;
t500 = t329 ^ 2;
t499 = pkin(2) * t388;
t498 = pkin(2) * t392;
t490 = g(3) * t380;
t488 = g(3) * t395;
t487 = t327 * pkin(4);
t486 = t394 * pkin(4);
t400 = -t394 * qJDD(2) + t391 * t402;
t283 = t400 + t503;
t481 = qJ(5) * t283;
t478 = t412 * t391;
t477 = t327 * t336;
t476 = t327 * t342;
t475 = t327 * t345;
t474 = t329 * t336;
t473 = t329 * t345;
t471 = t342 * t391;
t470 = t355 * t391;
t469 = t355 * t394;
t377 = pkin(3) + t486;
t468 = t377 * t380;
t467 = t379 * t396;
t348 = t388 * t358;
t466 = t391 * t313;
t303 = t394 * t313;
t325 = -t363 * t483 + t388 * t440;
t322 = t394 * t325;
t373 = pkin(7) + t499;
t461 = qJ(5) + t373;
t271 = t394 * t288 - t307 * t391;
t265 = -qJ(5) * t329 + t271;
t263 = pkin(4) * t336 + t265;
t460 = -t265 + t263;
t450 = qJD(4) * t394;
t459 = -t391 * t283 - t327 * t450;
t279 = t483 * t312 - t388 * t319;
t298 = pkin(2) * t452 + pkin(3) * t345 + pkin(7) * t342;
t321 = t357 * t483 + t348;
t458 = t391 * t298 + t394 * t321;
t413 = -t388 * t392 + t433;
t315 = -pkin(3) * t413 - pkin(7) * t355 - t378;
t457 = t391 * t315 + t322;
t431 = qJD(4) * t461;
t456 = -qJ(5) * t471 + qJD(5) * t394 - t391 * t431 - t458;
t292 = t394 * t298;
t455 = -pkin(4) * t345 - t292 + (-qJ(5) * t342 - t431) * t394 + (-qJD(5) + t321) * t391;
t454 = (g(1) * t462 + g(2) * t463) * t379;
t386 = t392 ^ 2;
t453 = -t395 ^ 2 + t386;
t444 = t392 * t484;
t297 = t341 * t483 + t388 * t409;
t347 = t413 * qJD(2);
t299 = pkin(3) * t344 - pkin(7) * t347 + t444;
t443 = t394 * t297 + t391 * t299 + t315 * t450;
t442 = t483 * pkin(2);
t441 = t355 * t450;
t277 = -qJDD(2) * pkin(3) - t279;
t260 = pkin(4) * t283 + qJDD(5) + t277;
t437 = -t260 - t490;
t432 = qJD(5) + t487;
t296 = t341 * t388 - t483 * t409;
t320 = t357 * t388 - t435;
t324 = -t363 * t388 - t390 * t434;
t430 = t336 * t394;
t429 = t391 * t273 + t394 * t278 + t288 * t450 - t307 * t451;
t374 = -t442 - pkin(3);
t428 = -g(1) * t337 - g(2) * t339;
t338 = -t380 * t463 + t464;
t340 = t380 * t462 + t465;
t427 = -g(1) * t338 - g(2) * t340;
t424 = qJD(4) * t373 * t336 + t277;
t316 = t351 * t483 + t348;
t389 = -qJ(5) - pkin(7);
t422 = t379 * t389 - t468;
t257 = -qJD(5) * t327 + t429 - t481;
t421 = -t263 * t336 + t257;
t420 = -qJ(5) * t347 - qJD(5) * t355;
t418 = t303 + (-t451 - t471) * t336;
t416 = -0.2e1 * pkin(1) * t449 - pkin(6) * qJDD(2);
t415 = t347 * t391 + t441;
t414 = t347 * t394 - t355 * t451;
t306 = -qJD(2) * pkin(3) - t316;
t411 = t306 * t336 - t373 * t313;
t335 = -qJDD(1) * t378 + t445;
t397 = qJD(2) ^ 2;
t408 = -pkin(6) * t397 + 0.2e1 * t480 + t504;
t398 = qJD(1) ^ 2;
t407 = pkin(1) * t398 - pkin(6) * qJDD(1) + t426;
t406 = g(1) * t340 - g(2) * t338 + t394 * t491 - t429;
t401 = t404 + t502;
t399 = t345 * t450 + t391 * t448 + t400;
t376 = t390 * t396;
t368 = t396 * t378;
t367 = t380 * t489;
t362 = t374 - t486;
t353 = t461 * t394;
t352 = t461 * t391;
t326 = t327 ^ 2;
t305 = t394 * t315;
t294 = pkin(4) * t470 + t324;
t293 = t394 * t299;
t285 = -pkin(4) * t471 + t320;
t284 = t306 + t432;
t276 = pkin(4) * t415 + t296;
t274 = -qJ(5) * t470 + t457;
t267 = -pkin(4) * t413 - qJ(5) * t469 - t325 * t391 + t305;
t259 = -qJ(5) * t441 + (-qJD(4) * t325 + t420) * t391 + t443;
t258 = pkin(4) * t344 - t297 * t391 + t293 + t420 * t394 + (-t322 + (qJ(5) * t355 - t315) * t391) * qJD(4);
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t386 + 0.2e1 * t395 * t439) * MDP(4) + 0.2e1 * (t392 * t446 - t449 * t453) * MDP(5) + (qJDD(2) * t392 + t395 * t397) * MDP(6) + (qJDD(2) * t395 - t392 * t397) * MDP(7) + (t392 * t416 + t395 * t408) * MDP(9) + (-t408 * t392 + t416 * t395) * MDP(10) + (-qJDD(2) * t324 - t318 * t378 - t335 * t413 + t344 * t361 + t504 * t380 + (t342 * t498 - t296) * qJD(2)) * MDP(11) + (g(2) * t467 - t297 * qJD(2) - t325 * qJDD(2) + t335 * t355 + t345 * t444 + t361 * t347 - t378 * t402 - t379 * t495) * MDP(12) + (-t279 * t355 + t280 * t413 + t296 * t345 - t297 * t342 - t316 * t347 - t317 * t344 - t325 * t318 + t324 * t402 - t426) * MDP(13) + (t280 * t325 + t317 * t297 - t279 * t324 - t316 * t296 - t335 * t378 + t361 * t444 - g(1) * (-t378 * t393 - t376) - g(2) * (-t390 * t393 + t368)) * MDP(14) + (t329 * t414 + t412 * t469) * MDP(15) + ((-t327 * t394 - t507) * t347 + (-t478 - t283 * t394 + (t327 * t391 - t329 * t394) * qJD(4)) * t355) * MDP(16) + (t303 * t355 + t329 * t344 + t336 * t414 - t412 * t413) * MDP(17) + (t283 * t413 - t327 * t344 - t336 * t415 - t355 * t466) * MDP(18) + (-t313 * t413 + t336 * t344) * MDP(19) + ((-t325 * t450 + t293) * t336 + t305 * t313 - (-t307 * t450 + t270) * t413 + t271 * t344 + t296 * t327 + t324 * t283 + t306 * t441 + ((-qJD(4) * t315 - t297) * t336 - t325 * t313 - (-qJD(4) * t288 - t278) * t413 + t277 * t355 + t306 * t347) * t391 + t427) * MDP(20) + (-(-t325 * t451 + t443) * t336 - t457 * t313 + t429 * t413 - t272 * t344 + t296 * t329 + t324 * t412 + t277 * t469 + t414 * t306 + t428) * MDP(21) + (-t256 * t413 + t258 * t336 + t260 * t470 + t263 * t344 + t267 * t313 + t276 * t327 + t283 * t294 + t284 * t415 + t427) * MDP(22) + (t257 * t413 - t259 * t336 + t260 * t469 - t266 * t344 - t274 * t313 + t276 * t329 + t284 * t414 + t294 * t412 + t428) * MDP(23) + (-t258 * t329 - t259 * t327 - t267 * t412 - t274 * t283 + t504 * t379 + (-t263 * t394 - t266 * t391) * t347 + (-t256 * t394 - t257 * t391 + (t263 * t391 - t266 * t394) * qJD(4)) * t355) * MDP(24) + (t257 * t274 + t266 * t259 + t256 * t267 + t263 * t258 + t260 * t294 + t284 * t276 - g(1) * (pkin(4) * t464 - t376) - g(2) * (-t389 * t467 + t396 * t468 + t368) + (-g(1) * (-t378 + t422) - g(2) * (pkin(4) * t391 - t390)) * t393) * MDP(25) + t504 * MDP(2) + t426 * MDP(3); MDP(6) * t447 + MDP(7) * t446 + qJDD(2) * MDP(8) + (t392 * t407 - t488) * MDP(9) + (g(3) * t392 + t395 * t407) * MDP(10) + (-t490 + t320 * qJD(2) - t361 * t345 + t417 + (qJDD(2) * t483 - t342 * t452) * pkin(2) + t279) * MDP(11) + (qJD(2) * t321 + t342 * t361 + (-qJDD(2) * t388 - t345 * t452) * pkin(2) - t280 + t501) * MDP(12) + (-t318 * t499 - t402 * t442 - (-t317 + t320) * t345 + (-t316 + t321) * t342) * MDP(13) + (t316 * t320 - t317 * t321 + (t483 * t279 - t488 + t280 * t388 + (-qJD(1) * t361 + t426) * t392) * pkin(2)) * MDP(14) + (t329 * t430 + t478) * MDP(15) + ((t412 - t476) * t394 - t506 + t459) * MDP(16) + (t336 * t430 + t466 - t473) * MDP(17) + (t418 + t475) * MDP(18) - t336 * t345 * MDP(19) + (-t271 * t345 + t374 * t283 - t292 * t336 - t320 * t327 + (-t424 - t490) * t394 + (t321 * t336 + t411) * t391 + t454) * MDP(20) + (t374 * t412 + t458 * t336 + t272 * t345 - t320 * t329 + t367 + t411 * t394 + (-t417 + t424) * t391) * MDP(21) + (-t263 * t345 + t283 * t362 - t285 * t327 - t313 * t352 + t437 * t394 + t455 * t336 + (t284 * t342 + (t284 + t487) * qJD(4)) * t391 + t454) * MDP(22) + (t266 * t345 + t412 * t362 - t285 * t329 - t313 * t353 + t367 - t456 * t336 + t284 * t430 + (pkin(4) * t503 + t260 - t417) * t391) * MDP(23) + (-t283 * t353 - t456 * t327 - t455 * t329 + t352 * t412 - t391 * t505 + t421 * t394 - t501) * MDP(24) + (t257 * t353 - t256 * t352 + t260 * t362 - g(3) * (-t422 + t485) + (pkin(4) * t451 - t285) * t284 + t456 * t266 + t455 * t263 + t426 * (t377 * t379 + t380 * t389 + t498)) * MDP(25) + (-MDP(4) * t392 * t395 + MDP(5) * t453) * t398; (0.2e1 * qJD(2) * t345 + t423) * MDP(11) + ((t370 - t342) * qJD(2) + t405) * MDP(12) + (-t342 ^ 2 - t345 ^ 2) * MDP(13) + (t316 * t345 + t317 * t342 + t335 - t504) * MDP(14) + ((-t412 - t476) * t394 + t506 + t459) * MDP(24) + (-t284 * t345 + t421 * t391 + t394 * t505 - t504) * MDP(25) + (MDP(21) + MDP(23)) * (-t336 ^ 2 * t394 - t466 - t473) + (MDP(20) + MDP(22)) * (t418 - t475); t329 * t327 * MDP(15) + (-t326 + t500) * MDP(16) + (t412 + t477) * MDP(17) + (-t399 + t474) * MDP(18) + t313 * MDP(19) + (t272 * t336 - t306 * t329 + t401) * MDP(20) + (t271 * t336 + t306 * t327 + t406) * MDP(21) + (0.2e1 * t497 - t482 + t479 + (-t284 - t432) * t329 + t401) * MDP(22) + (-pkin(4) * t500 + t481 + t265 * t336 + (qJD(5) + t284) * t327 + t406) * MDP(23) + (-pkin(4) * t412 - t327 * t460) * MDP(24) + (t460 * t266 + (-t284 * t329 + t256 + t502) * pkin(4)) * MDP(25); (t399 + t474) * MDP(22) + (t412 - t477) * MDP(23) + (-t326 - t500) * MDP(24) + (t263 * t329 + t266 * t327 - t417 - t437) * MDP(25);];
tau = t1;
