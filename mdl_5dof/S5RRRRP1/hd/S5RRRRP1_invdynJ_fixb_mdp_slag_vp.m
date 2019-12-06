% Calculate vector of inverse dynamics joint torques for
% S5RRRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:09
% EndTime: 2019-12-05 18:46:15
% DurationCPUTime: 3.54s
% Computational Cost: add. (3607->346), mult. (8570->446), div. (0->0), fcn. (6172->12), ass. (0->171)
t420 = qJ(2) + qJ(3);
t414 = qJ(4) + t420;
t402 = sin(t414);
t403 = cos(t414);
t424 = sin(qJ(1));
t428 = cos(qJ(1));
t465 = g(1) * t428 + g(2) * t424;
t534 = -g(3) * t403 + t465 * t402;
t426 = cos(qJ(3));
t422 = sin(qJ(3));
t423 = sin(qJ(2));
t496 = qJD(1) * t423;
t482 = t422 * t496;
t427 = cos(qJ(2));
t495 = qJD(1) * t427;
t356 = -t426 * t495 + t482;
t358 = -t422 * t495 - t426 * t496;
t421 = sin(qJ(4));
t425 = cos(qJ(4));
t325 = -t425 * t356 + t358 * t421;
t321 = t325 ^ 2;
t416 = qJDD(2) + qJDD(3);
t409 = qJDD(4) + t416;
t417 = qJD(2) + qJD(3);
t410 = qJD(4) + t417;
t487 = qJDD(1) * t427;
t490 = qJD(1) * qJD(2);
t480 = t427 * t490;
t488 = qJDD(1) * t423;
t456 = -t480 - t488;
t489 = qJD(1) * qJD(3);
t530 = -t427 * t489 + t456;
t310 = -t417 * t482 + t422 * t487 - t426 * t530;
t368 = t422 * t427 + t423 * t426;
t452 = t368 * qJD(2);
t440 = -t368 * qJD(3) - t452;
t435 = t440 * qJD(1);
t459 = t422 * t423 - t426 * t427;
t431 = -t459 * qJDD(1) + t435;
t461 = t356 * t421 + t425 * t358;
t445 = t461 * qJD(4) - t310 * t421 + t425 * t431;
t491 = qJD(4) * t425;
t492 = qJD(4) * t421;
t453 = t425 * t310 - t356 * t491 + t358 * t492 + t421 * t431;
t521 = t461 ^ 2;
t533 = t409 * MDP(22) + t325 * MDP(18) * t461 + (-t325 * t410 + t453) * MDP(20) + (-t410 * t461 + t445) * MDP(21) + (-t321 + t521) * MDP(19);
t415 = t427 * pkin(2);
t527 = -pkin(1) - t415;
t511 = qJ(5) * t325;
t318 = t461 * qJ(5);
t343 = t459 * pkin(3) + t527;
t481 = t423 * t490;
t398 = pkin(2) * t481;
t294 = -pkin(3) * t435 + t343 * qJDD(1) + t398;
t531 = -pkin(4) * t445 + qJDD(5) + t294;
t380 = t527 * qJD(1);
t529 = qJDD(1) * t527;
t340 = pkin(3) * t356 + t380;
t520 = pkin(6) + pkin(7);
t338 = qJDD(2) * pkin(2) + t520 * t456;
t455 = -t481 + t487;
t339 = t520 * t455;
t382 = t520 * t427;
t374 = qJD(1) * t382;
t363 = t426 * t374;
t381 = t520 * t423;
t372 = qJD(1) * t381;
t512 = qJD(2) * pkin(2);
t365 = -t372 + t512;
t460 = -t422 * t365 - t363;
t444 = t460 * qJD(3) + t426 * t338 - t422 * t339;
t272 = pkin(3) * t416 - pkin(8) * t310 + t444;
t494 = qJD(3) * t422;
t354 = t374 * t494;
t470 = -qJD(3) * t365 - t339;
t276 = -t354 + (pkin(8) * t530 + t338) * t422 + ((-t423 * t489 + t455) * pkin(8) - t470) * t426;
t351 = t358 * pkin(8);
t359 = t422 * t374;
t472 = t426 * t365 - t359;
t308 = t351 + t472;
t297 = pkin(3) * t417 + t308;
t518 = pkin(8) * t356;
t309 = -t460 - t518;
t522 = (qJD(4) * t297 + t276) * t425 + t421 * t272 - t309 * t492;
t439 = g(3) * t402 - t340 * t325 + t403 * t465 - t522;
t301 = t425 * t309;
t462 = -t421 * t297 - t301;
t446 = t462 * qJD(4) + t425 * t272 - t421 * t276;
t433 = t340 * t461 + t446 + t534;
t370 = t426 * t381;
t524 = -pkin(8) * t368 - t422 * t382;
t319 = -t370 + t524;
t500 = -t422 * t381 + t426 * t382;
t320 = -t459 * pkin(8) + t500;
t504 = t421 * t319 + t425 * t320;
t471 = t372 * t422 - t363;
t312 = t471 + t518;
t501 = -t426 * t372 - t359;
t313 = t351 + t501;
t405 = t426 * pkin(2) + pkin(3);
t510 = t421 * t422;
t526 = -t405 * t491 - (-t422 * t492 + (t425 * t426 - t510) * qJD(3)) * pkin(2) + t421 * t312 + t425 * t313;
t508 = t422 * t425;
t525 = -t405 * t492 + (-t422 * t491 + (-t421 * t426 - t508) * qJD(3)) * pkin(2) - t425 * t312 + t313 * t421;
t519 = pkin(4) * t461;
t299 = t421 * t309;
t477 = t425 * t297 - t299;
t274 = t477 + t318;
t270 = pkin(4) * t410 + t274;
t507 = t270 - t274;
t506 = t425 * t308 - t299;
t503 = -t318 - t526;
t502 = t511 + t525;
t412 = cos(t420);
t499 = pkin(3) * t412 + pkin(4) * t403;
t418 = t423 ^ 2;
t498 = -t427 ^ 2 + t418;
t493 = qJD(3) * t426;
t408 = t423 * t512;
t483 = qJD(2) * t520;
t373 = t423 * t483;
t375 = t427 * t483;
t486 = -t426 * t373 - t422 * t375 - t381 * t493;
t485 = t415 + t499;
t341 = pkin(2) * t496 - pkin(3) * t358;
t476 = -t308 * t421 - t301;
t473 = t425 * t319 - t320 * t421;
t352 = t398 + t529;
t468 = t352 + t529;
t467 = -pkin(2) * t510 + t425 * t405;
t411 = sin(t420);
t466 = -pkin(3) * t411 - pkin(4) * t402;
t464 = g(1) * t424 - g(2) * t428;
t275 = -t462 + t511;
t463 = t270 * t325 - t275 * t461;
t458 = -0.2e1 * pkin(1) * t490 - pkin(6) * qJDD(2);
t457 = t425 * t459;
t290 = -pkin(8) * t452 + t524 * qJD(3) + t486;
t337 = t417 * t459;
t364 = t426 * t375;
t291 = pkin(8) * t337 - t500 * qJD(3) + t373 * t422 - t364;
t454 = t425 * t290 + t421 * t291 + t319 * t491 - t320 * t492;
t429 = qJD(2) ^ 2;
t449 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t429 + t464;
t430 = qJD(1) ^ 2;
t448 = pkin(1) * t430 - pkin(6) * qJDD(1) + t465;
t336 = t425 * t368 - t421 * t459;
t447 = 0.2e1 * t380 * t417;
t443 = -t504 * qJD(4) - t290 * t421 + t425 * t291;
t438 = g(3) * t411 - t422 * t338 + t380 * t356 + t412 * t465 + t470 * t426 + t354;
t436 = -t358 * t356 * MDP(11) + (t356 * t417 + t310) * MDP(13) + (-t358 * t417 + t431) * MDP(14) + (-t356 ^ 2 + t358 ^ 2) * MDP(12) + t416 * MDP(15) + t533;
t328 = -t440 * pkin(3) + t408;
t432 = -g(3) * t412 + t380 * t358 + t411 * t465 + t444;
t413 = -qJ(5) - pkin(8) - t520;
t404 = pkin(3) * t425 + pkin(4);
t353 = pkin(2) * t508 + t405 * t421;
t349 = pkin(4) + t467;
t344 = pkin(1) + t485;
t335 = t368 * t421 + t457;
t295 = -pkin(4) * t325 + qJD(5) + t340;
t286 = t336 * qJD(4) - t421 * t337 - t425 * t440;
t285 = qJD(4) * t457 + t425 * t337 + t368 * t492 - t421 * t440;
t284 = -qJ(5) * t335 + t504;
t283 = -qJ(5) * t336 + t473;
t278 = t318 + t506;
t277 = t476 - t511;
t267 = qJ(5) * t285 - qJD(5) * t336 + t443;
t266 = -qJ(5) * t286 - qJD(5) * t335 + t454;
t265 = qJ(5) * t445 + qJD(5) * t325 + t522;
t264 = pkin(4) * t409 - qJ(5) * t453 + qJD(5) * t461 + t446;
t1 = [qJDD(1) * MDP(1) + t464 * MDP(2) + t465 * MDP(3) + (qJDD(1) * t418 + 0.2e1 * t423 * t480) * MDP(4) + 0.2e1 * (t423 * t487 - t498 * t490) * MDP(5) + (qJDD(2) * t423 + t427 * t429) * MDP(6) + (qJDD(2) * t427 - t423 * t429) * MDP(7) + (t458 * t423 + t449 * t427) * MDP(9) + (-t449 * t423 + t458 * t427) * MDP(10) + (t310 * t368 + t337 * t358) * MDP(11) + (-t310 * t459 + t337 * t356 - t358 * t440 + t431 * t368) * MDP(12) + (-t337 * t417 + t368 * t416) * MDP(13) + (-t459 * t416 + t440 * t417) * MDP(14) + (t356 * t408 - t364 * t417 - t370 * t416 + t464 * t412 + (-qJD(3) * t382 * t417 + t447 * t423 - t468 * t427) * t426 + ((qJD(3) * t381 + t373) * t417 - t382 * t416 + t468 * t423 + t447 * t427) * t422) * MDP(16) + (-t358 * t408 + t527 * t310 + t352 * t368 - t380 * t337 - (-t382 * t494 + t486) * t417 - t500 * t416 - t464 * t411) * MDP(17) + (t285 * t461 + t336 * t453) * MDP(18) + (-t285 * t325 + t286 * t461 - t335 * t453 + t336 * t445) * MDP(19) + (-t285 * t410 + t336 * t409) * MDP(20) + (-t286 * t410 - t335 * t409) * MDP(21) + (t340 * t286 + t294 * t335 - t325 * t328 - t343 * t445 + t464 * t403 + t473 * t409 + t443 * t410) * MDP(23) + (-t340 * t285 + t294 * t336 - t328 * t461 + t343 * t453 - t464 * t402 - t504 * t409 - t454 * t410) * MDP(24) + (-t264 * t336 - t265 * t335 + t266 * t325 + t267 * t461 + t270 * t285 - t275 * t286 - t283 * t453 + t284 * t445 - t465) * MDP(25) + (t265 * t284 + t275 * t266 + t264 * t283 + t270 * t267 + t531 * (t335 * pkin(4) + t343) + t295 * (t286 * pkin(4) + t328) - g(1) * (-t344 * t424 - t413 * t428) - g(2) * (t344 * t428 - t413 * t424)) * MDP(26); (t265 * t353 + t264 * t349 - t295 * (t341 - t519) - g(3) * t485 - t465 * (-pkin(2) * t423 + t466) + t503 * t275 + t502 * t270) * MDP(26) + (t501 * t417 + (t358 * t496 - t422 * t416 - t417 * t493) * pkin(2) + t438) * MDP(17) + (t325 * t503 - t349 * t453 + t353 * t445 + t461 * t502 + t463) * MDP(25) + (-g(3) * t427 + t448 * t423) * MDP(9) + (g(3) * t423 + t448 * t427) * MDP(10) + (t341 * t461 - t353 * t409 + t526 * t410 + t439) * MDP(24) + qJDD(2) * MDP(8) + MDP(7) * t487 + MDP(6) * t488 + t436 + (-t471 * t417 + (-t356 * t496 + t426 * t416 - t417 * t494) * pkin(2) + t432) * MDP(16) + (t325 * t341 + t467 * t409 + t525 * t410 + t433) * MDP(23) + (-t423 * t427 * MDP(4) + t498 * MDP(5)) * t430; (t264 * t404 - t275 * t278 - t270 * t277 + t295 * t519 - g(3) * t499 - t465 * t466 + (t265 * t421 + t295 * t358 + (-t270 * t421 + t275 * t425) * qJD(4)) * pkin(3)) * MDP(26) + (t506 * t410 + (-t358 * t461 - t409 * t421 - t410 * t491) * pkin(3) + t439) * MDP(24) + (-t277 * t461 - t278 * t325 - t453 * t404 + (t445 * t421 + (t325 * t425 - t421 * t461) * qJD(4)) * pkin(3) + t463) * MDP(25) + (t472 * t417 + t438) * MDP(17) + t436 + (-t460 * t417 + t432) * MDP(16) + (-t476 * t410 + (-t325 * t358 + t409 * t425 - t410 * t492) * pkin(3) + t433) * MDP(23); (-t462 * t410 + t433) * MDP(23) + (t477 * t410 + t439) * MDP(24) + (-pkin(4) * t453 + t325 * t507) * MDP(25) + (t507 * t275 + (t295 * t461 + t264 + t534) * pkin(4)) * MDP(26) + t533; (-t321 - t521) * MDP(25) + (-t270 * t461 - t275 * t325 - t464 + t531) * MDP(26);];
tau = t1;
