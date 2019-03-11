% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRRPP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:49
% EndTime: 2019-03-08 22:52:59
% DurationCPUTime: 4.84s
% Computational Cost: add. (3233->463), mult. (7868->592), div. (0->0), fcn. (5340->8), ass. (0->191)
t552 = MDP(20) - MDP(25);
t411 = sin(qJ(3));
t414 = cos(qJ(3));
t444 = pkin(3) * t411 - pkin(9) * t414;
t377 = t444 * qJD(3);
t382 = -pkin(3) * t414 - pkin(9) * t411 - pkin(2);
t410 = sin(qJ(4));
t412 = sin(qJ(2));
t413 = cos(qJ(4));
t479 = qJD(4) * t413;
t408 = sin(pkin(6));
t488 = qJD(1) * t408;
t415 = cos(qJ(2));
t504 = t414 * t415;
t551 = -(t410 * t412 + t413 * t504) * t488 + t410 * t377 + t382 * t479;
t484 = qJD(2) * t414;
t396 = -qJD(4) + t484;
t469 = qJD(2) * qJD(3);
t455 = t411 * t469;
t550 = qJ(5) * t455 - t396 * qJD(5);
t505 = t413 * t414;
t399 = pkin(8) * t505;
t480 = qJD(4) * t410;
t549 = qJD(4) * t399 + t382 * t480 - (t410 * t504 - t412 * t413) * t488;
t533 = MDP(19) + MDP(23);
t548 = -MDP(17) - t533;
t547 = MDP(21) + MDP(24);
t546 = MDP(5) * t411;
t406 = t411 ^ 2;
t545 = MDP(6) * (-t414 ^ 2 + t406);
t462 = t412 * t488;
t379 = qJD(2) * pkin(8) + t462;
t409 = cos(pkin(6));
t486 = qJD(2) * t408;
t459 = t415 * t486;
t482 = qJD(3) * t411;
t487 = qJD(1) * t414;
t316 = -t379 * t482 + (qJD(3) * t409 + t459) * t487;
t511 = t409 * t411;
t392 = qJD(1) * t511;
t351 = t414 * t379 + t392;
t338 = qJD(3) * pkin(9) + t351;
t349 = (t377 + t462) * qJD(2);
t461 = t415 * t488;
t352 = qJD(2) * t382 - t461;
t450 = -t410 * t316 - t338 * t479 + t413 * t349 - t352 * t480;
t281 = -pkin(4) * t455 - t450;
t302 = t413 * t338 + t410 * t352;
t389 = t396 * qJ(5);
t296 = -t389 + t302;
t544 = t296 * t396 + t281;
t468 = qJD(3) * qJD(4);
t453 = t410 * t468;
t456 = t411 * t479;
t481 = qJD(3) * t414;
t346 = qJD(2) * (t410 * t481 + t456) + t453;
t473 = t413 * qJD(3);
t485 = qJD(2) * t411;
t372 = t410 * t485 - t473;
t542 = t346 * qJ(6) + t372 * qJD(6);
t454 = t414 * t469;
t457 = t411 * t480;
t345 = qJD(2) * t457 + (-t454 - t468) * t413;
t515 = t372 * t396;
t541 = -t345 + t515;
t483 = qJD(3) * t410;
t374 = t413 * t485 + t483;
t540 = t374 * t396 - t346;
t539 = -t377 * t413 + t549;
t538 = -qJD(5) * t410 - t351;
t350 = -t411 * t379 + t409 * t487;
t537 = 0.2e1 * t550;
t534 = qJ(5) * t482 - qJD(5) * t414 + t551;
t301 = -t410 * t338 + t413 * t352;
t471 = qJD(5) - t301;
t526 = qJ(5) * t410;
t531 = pkin(4) + pkin(5);
t532 = -t413 * t531 - t526;
t370 = t374 ^ 2;
t530 = pkin(9) - qJ(6);
t529 = qJD(2) * pkin(2);
t528 = qJ(5) * t346;
t527 = qJ(5) * t372;
t525 = qJ(5) * t413;
t524 = qJ(6) * t411;
t448 = t411 * t459;
t317 = qJD(1) * t448 + qJD(3) * t392 + t379 * t481;
t283 = t346 * pkin(4) + t345 * qJ(5) - t374 * qJD(5) + t317;
t523 = t283 * t410;
t522 = t283 * t413;
t294 = qJ(6) * t372 + t302;
t289 = t294 - t389;
t521 = t289 * t396;
t519 = t317 * t410;
t518 = t317 * t413;
t445 = qJD(3) * pkin(3) + t350;
t517 = t445 * t410;
t516 = t345 * t410;
t513 = t408 * t412;
t512 = t408 * t415;
t510 = t410 * t396;
t509 = t410 * t414;
t508 = t411 * t413;
t417 = qJD(3) ^ 2;
t507 = t411 * t417;
t506 = t413 * t396;
t503 = t414 * t417;
t463 = -pkin(8) * t410 - pkin(4);
t476 = qJD(6) * t413;
t502 = (-qJ(6) * t481 - t377) * t413 + (qJ(6) * t480 - t476 + (-pkin(5) + t463) * qJD(3)) * t411 + t549;
t501 = (-pkin(8) * qJD(3) + qJ(6) * qJD(4)) * t508 + (qJD(6) * t411 + (-pkin(8) * qJD(4) + qJ(6) * qJD(3)) * t414) * t410 + t534;
t385 = t530 * t413;
t376 = t444 * qJD(2);
t452 = -t410 * t350 + t376 * t413;
t500 = (-qJ(6) * t505 - t411 * t531) * qJD(2) - t452 - qJD(4) * t385 + qJD(6) * t410;
t499 = (-t411 * t473 - t414 * t480) * pkin(8) + t534;
t493 = t413 * t350 + t410 * t376;
t305 = qJ(5) * t485 + t493;
t498 = qJ(6) * t410 * t484 + t480 * t530 + t305 + t476;
t497 = t463 * t482 + t539;
t435 = -t410 * t531 + t525;
t496 = t396 * t435 + t538;
t440 = pkin(4) * t410 - t525;
t494 = t396 * t440 - t538;
t490 = t410 * t382 + t399;
t477 = qJD(5) * t413;
t475 = t445 * qJD(4);
t474 = t411 * MDP(16);
t293 = qJ(6) * t374 + t301;
t472 = qJD(5) - t293;
t426 = qJ(5) * t374 + t445;
t292 = -t372 * t531 + qJD(6) + t426;
t470 = qJD(6) + t292;
t467 = pkin(9) * t510;
t466 = pkin(9) * t506;
t465 = pkin(9) * t482;
t464 = pkin(9) * t473;
t460 = t412 * t486;
t458 = t396 * t480;
t398 = pkin(8) * t509;
t451 = t382 * t413 - t398;
t449 = t411 * t461;
t447 = t372 * t461;
t446 = t374 * t461;
t343 = -qJ(5) * t414 + t490;
t380 = -t461 - t529;
t442 = -t380 - t461;
t441 = pkin(4) * t413 + t526;
t295 = pkin(4) * t396 + t471;
t439 = t295 * t413 - t296 * t410;
t438 = qJD(2) * t406 - t396 * t414;
t437 = pkin(8) + t440;
t360 = t414 * t513 + t511;
t327 = t360 * t410 + t413 * t512;
t328 = t360 * t413 - t410 * t512;
t359 = -t409 * t414 + t411 * t513;
t282 = -pkin(5) * t346 - t283;
t433 = -t282 * t410 - t292 * t479;
t432 = t282 * t413 - t292 * t480;
t430 = qJ(6) * t345 - t450;
t428 = -t413 * t316 + t338 * t480 - t410 * t349 - t352 * t479;
t427 = -pkin(8) + t435;
t280 = -t428 + t550;
t423 = qJD(3) * (-t442 - t529);
t422 = -t301 * t396 + t428;
t421 = t345 + t515;
t419 = -t455 * t531 + t430;
t418 = qJD(2) ^ 2;
t405 = t414 * pkin(4);
t384 = t530 * t410;
t381 = -pkin(3) - t441;
t369 = pkin(3) - t532;
t354 = t437 * t411;
t344 = t405 - t451;
t342 = t427 * t411;
t326 = qJD(3) * t360 + t448;
t325 = -qJD(3) * t359 + t414 * t459;
t322 = pkin(4) * t374 + t527;
t320 = t410 * t524 + t343;
t318 = pkin(5) * t414 + t398 + t405 + (-t382 - t524) * t413;
t310 = -t374 * t531 - t527;
t308 = (qJD(4) * t441 - t477) * t411 + t437 * t481;
t307 = -pkin(4) * t485 - t452;
t303 = pkin(4) * t372 - t426;
t297 = (qJD(4) * t532 + t477) * t411 + t427 * t481;
t291 = -qJD(4) * t327 + t325 * t413 + t410 * t460;
t290 = qJD(4) * t328 + t325 * t410 - t413 * t460;
t284 = t396 * t531 + t472;
t279 = -qJD(6) * t374 + t419;
t278 = t280 + t542;
t1 = [(t280 * t328 + t281 * t327 + t283 * t359 + t290 * t295 + t291 * t296 + t303 * t326) * MDP(22) + (t278 * t328 + t279 * t327 - t282 * t359 + t284 * t290 + t289 * t291 - t292 * t326) * MDP(26) + (t548 * t327 * t485 - t326 * MDP(10) - t325 * MDP(11)) * qJD(3) + ((-MDP(10) * t411 - MDP(11) * t414) * t415 * t469 + (-t415 * MDP(4) + (-MDP(10) * t414 + MDP(11) * t411 - MDP(3)) * t412) * t418) * t408 + t552 * (t290 * t374 - t291 * t372 - t327 * t345 - t328 * t346) + (MDP(18) - t547) * (t291 * t396 + t326 * t374 - t328 * t455 - t345 * t359) - t548 * (t290 * t396 + t326 * t372 + t359 * t346); 0.2e1 * t454 * t546 - 0.2e1 * t469 * t545 + MDP(7) * t503 - MDP(8) * t507 + (-pkin(8) * t503 + t411 * t423) * MDP(10) + (pkin(8) * t507 + t414 * t423) * MDP(11) + (-t345 * t508 + (t414 * t473 - t457) * t374) * MDP(12) + ((-t372 * t413 - t374 * t410) * t481 + (t516 - t346 * t413 + (t372 * t410 - t374 * t413) * qJD(4)) * t411) * MDP(13) + (t396 * t457 + t345 * t414 + (t374 * t411 + t413 * t438) * qJD(3)) * MDP(14) + (t396 * t456 + t346 * t414 + (-t372 * t411 - t410 * t438) * qJD(3)) * MDP(15) + (-t396 - t484) * qJD(3) * t474 + (t539 * t396 + ((pkin(8) * t372 - t517) * qJD(3) - t450) * t414 + (-t447 - t413 * t475 + pkin(8) * t346 + t519 + (-pkin(8) * t510 + qJD(2) * t451 + t301) * qJD(3)) * t411) * MDP(17) + (t551 * t396 + (-t445 * t473 + (qJD(3) * t374 - t458) * pkin(8) - t428) * t414 + (-t446 + t410 * t475 - pkin(8) * t345 + t518 + (-pkin(8) * t506 - qJD(2) * t490 - t302) * qJD(3)) * t411) * MDP(18) + (t308 * t372 + t346 * t354 + (t303 * t483 + t281) * t414 + t497 * t396 + (-t447 + t303 * t479 + t523 + (-qJD(2) * t344 - t295) * qJD(3)) * t411) * MDP(19) + (-t343 * t346 - t344 * t345 + t497 * t374 - t499 * t372 + t439 * t481 + (-t280 * t410 + t281 * t413 + (-t295 * t410 - t296 * t413) * qJD(4)) * t411) * MDP(20) + (-t308 * t374 + t345 * t354 + (-t303 * t473 - t280) * t414 - t499 * t396 + (t446 + t303 * t480 - t522 + (qJD(2) * t343 + t296) * qJD(3)) * t411) * MDP(21) + (t280 * t343 + t281 * t344 + t283 * t354 + (t308 - t449) * t303 + t499 * t296 + t497 * t295) * MDP(22) + (-t297 * t372 - t342 * t346 + (-t292 * t483 + t279) * t414 + t502 * t396 + (-t447 + (-qJD(2) * t318 - t284) * qJD(3) + t433) * t411) * MDP(23) + (t297 * t374 - t342 * t345 + (t292 * t473 - t278) * t414 - t501 * t396 + (t446 + (qJD(2) * t320 + t289) * qJD(3) + t432) * t411) * MDP(24) + (t318 * t345 + t320 * t346 - t502 * t374 + t501 * t372 + (-t284 * t413 + t289 * t410) * t481 + (t278 * t410 - t279 * t413 + (t284 * t410 + t289 * t413) * qJD(4)) * t411) * MDP(25) + (t278 * t320 + t279 * t318 + t282 * t342 + (t297 + t449) * t292 + t501 * t289 + t502 * t284) * MDP(26); (qJD(3) * t351 - t380 * t485 - t317) * MDP(10) + t442 * t484 * MDP(11) + (-t374 * t506 - t516) * MDP(12) + (t410 * t540 + t413 * t541) * MDP(13) + (-t396 * t479 + (t396 * t505 + (-t374 + t483) * t411) * qJD(2)) * MDP(14) + (t458 + (-t396 * t509 + (t372 + t473) * t411) * qJD(2)) * MDP(15) + t396 * qJD(2) * t474 + (-pkin(3) * t346 - t518 + t452 * t396 - t351 * t372 + (t466 - t517) * qJD(4) + (-t301 * t411 + (t414 * t445 - t465) * t410) * qJD(2)) * MDP(17) + (pkin(3) * t345 + t519 - t493 * t396 - t351 * t374 + (-t413 * t445 - t467) * qJD(4) + (t445 * t505 + (t302 - t464) * t411) * qJD(2)) * MDP(18) + (-t522 - t307 * t396 + t346 * t381 - t494 * t372 + (t303 * t410 + t466) * qJD(4) + (t295 * t411 + (-t303 * t414 - t465) * t410) * qJD(2)) * MDP(19) + (t305 * t372 - t307 * t374 + (t280 - t396 * t295 + (qJD(4) * t374 - t346) * pkin(9)) * t413 + ((qJD(4) * t372 - t345) * pkin(9) + t544) * t410) * MDP(20) + (-t523 + t305 * t396 + t345 * t381 + t494 * t374 + (-t303 * t413 + t467) * qJD(4) + (t303 * t505 + (-t296 + t464) * t411) * qJD(2)) * MDP(21) + (t283 * t381 - t295 * t307 - t296 * t305 - t494 * t303 + (qJD(4) * t439 + t280 * t413 + t281 * t410) * pkin(9)) * MDP(22) + (-t346 * t369 - t500 * t396 + t496 * t372 + (t292 * t509 + (-qJD(3) * t384 + t284) * t411) * qJD(2) + t432) * MDP(23) + (-t345 * t369 + t498 * t396 - t496 * t374 + (-t292 * t505 + (qJD(3) * t385 - t289) * t411) * qJD(2) - t433) * MDP(24) + (t345 * t384 + t346 * t385 + t500 * t374 - t498 * t372 + (t284 * t396 - t278) * t413 + (-t279 - t521) * t410) * MDP(25) + (t278 * t385 + t279 * t384 + t282 * t369 - t284 * t500 - t289 * t498 - t292 * t496) * MDP(26) + (-t414 * t546 + t545) * t418; -t421 * MDP(14) - MDP(15) * t453 + t422 * MDP(18) + (pkin(4) * t345 - t528) * MDP(20) + (-t422 + t537) * MDP(21) + (-pkin(4) * t281 + qJ(5) * t280 - t295 * t302 + t296 * t471 - t303 * t322) * MDP(22) + (-t294 * t396 - t430) * MDP(23) + (t293 * t396 - t428 + t537 + t542) * MDP(24) + (-t345 * t531 + t528) * MDP(25) + (qJ(5) * t278 - t279 * t531 - t284 * t294 + t289 * t472 - t292 * t310) * MDP(26) + (-t396 * MDP(15) + t445 * MDP(17) - t303 * MDP(19) + (t296 - t302) * MDP(20) + t322 * MDP(21) + t470 * MDP(23) - t310 * MDP(24) + (-t289 + t294) * MDP(25) + MDP(13) * t374) * t374 + (-MDP(15) * t456 + (-MDP(15) * t509 + (0.2e1 * pkin(4) * MDP(19) + 0.2e1 * t531 * MDP(23) + MDP(16)) * t411) * qJD(3)) * qJD(2) + (t374 * MDP(12) - t445 * MDP(18) - t322 * MDP(19) + (t295 - t471) * MDP(20) - t303 * MDP(21) + t310 * MDP(23) + t292 * MDP(24) + (-t284 + t472) * MDP(25) - MDP(13) * t372) * t372 + (MDP(17) + MDP(19)) * (-t302 * t396 + t450); (t303 * t374 + t544) * MDP(22) + (-t470 * t374 + t419 + t521) * MDP(26) + t533 * (t372 * t374 - t455) + t547 * (-t396 ^ 2 - t370) - t552 * t421; t540 * MDP(23) + t541 * MDP(24) + (-t372 ^ 2 - t370) * MDP(25) + (t284 * t374 - t289 * t372 + t282) * MDP(26);];
tauc  = t1;
