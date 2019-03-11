% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRPP3
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
%   see S6PRRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:37
% EndTime: 2019-03-08 22:58:46
% DurationCPUTime: 4.66s
% Computational Cost: add. (3249->474), mult. (7873->602), div. (0->0), fcn. (5333->8), ass. (0->192)
t415 = sin(qJ(3));
t418 = cos(qJ(3));
t443 = pkin(3) * t415 - pkin(9) * t418;
t375 = t443 * qJD(3);
t381 = -pkin(3) * t418 - pkin(9) * t415 - pkin(2);
t414 = sin(qJ(4));
t416 = sin(qJ(2));
t417 = cos(qJ(4));
t484 = qJD(4) * t417;
t411 = sin(pkin(6));
t493 = qJD(1) * t411;
t419 = cos(qJ(2));
t508 = t418 * t419;
t556 = -(t414 * t416 + t417 * t508) * t493 + t414 * t375 + t381 * t484;
t457 = t415 * t484;
t486 = qJD(3) * t418;
t555 = t414 * t486 + t457;
t485 = qJD(4) * t414;
t458 = t415 * t485;
t482 = t417 * qJD(3);
t554 = -t418 * t482 + t458;
t541 = MDP(21) + MDP(24);
t553 = MDP(18) - t541;
t466 = t416 * t493;
t377 = qJD(2) * pkin(8) + t466;
t412 = cos(pkin(6));
t514 = t412 * t415;
t394 = qJD(1) * t514;
t351 = t418 * t377 + t394;
t340 = qJD(3) * pkin(9) + t351;
t465 = t419 * t493;
t352 = qJD(2) * t381 - t465;
t301 = t340 * t414 - t417 * t352;
t488 = qJD(3) * t414;
t490 = qJD(2) * t415;
t370 = t417 * t490 + t488;
t433 = pkin(5) * t370 + t301;
t480 = qJD(5) + t433;
t552 = -MDP(20) + MDP(25);
t551 = MDP(5) * t418;
t409 = t415 ^ 2;
t550 = MDP(6) * (-t418 ^ 2 + t409);
t489 = qJD(2) * t418;
t397 = -qJD(4) + t489;
t386 = qJD(5) * t397;
t478 = qJD(2) * qJD(3);
t402 = t415 * t478;
t396 = qJ(5) * t402;
t477 = qJD(3) * qJD(4);
t456 = t414 * t477;
t346 = qJD(2) * t555 + t456;
t491 = qJD(2) * t411;
t462 = t419 * t491;
t487 = qJD(3) * t415;
t492 = qJD(1) * t418;
t319 = -t377 * t487 + (qJD(3) * t412 + t462) * t492;
t349 = (t375 + t466) * qJD(2);
t450 = t417 * t319 - t340 * t485 + t414 * t349 + t352 * t484;
t429 = -pkin(5) * t346 + t450;
t280 = -t386 + t396 + t429;
t536 = pkin(4) + qJ(6);
t285 = t397 * t536 + t480;
t549 = -t285 * t397 + t280;
t449 = -t414 * t319 - t340 * t484 + t417 * t349 - t352 * t485;
t283 = -pkin(4) * t402 - t449;
t302 = t417 * t340 + t414 * t352;
t297 = qJ(5) * t397 - t302;
t548 = -t297 * t397 + t283;
t509 = t417 * t418;
t401 = pkin(8) * t509;
t547 = -(t414 * t508 - t416 * t417) * t493 + qJD(4) * t401 - t375 * t417 + t381 * t485;
t546 = qJD(5) * t418 - t556;
t350 = -t415 * t377 + t412 * t492;
t544 = qJD(5) * t414 + t351 + (t414 * t489 - t485) * pkin(4);
t543 = MDP(17) - MDP(20);
t542 = MDP(19) + MDP(23);
t395 = t397 ^ 2;
t538 = pkin(5) + pkin(9);
t368 = t414 * t490 - t482;
t537 = pkin(5) * t368;
t535 = qJD(2) * pkin(2);
t534 = qJ(5) * t346;
t533 = qJ(5) * t368;
t532 = qJ(5) * t417;
t446 = t415 * t462;
t320 = qJD(1) * t446 + qJD(3) * t394 + t377 * t486;
t345 = qJD(2) * t554 - t417 * t477;
t425 = qJ(5) * t345 - qJD(5) * t370 + t320;
t284 = pkin(4) * t346 + t425;
t531 = t284 * t414;
t530 = t284 * t417;
t525 = t320 * t414;
t524 = t320 * t417;
t339 = -qJD(3) * pkin(3) - t350;
t523 = t339 * t414;
t522 = t345 * t414;
t521 = t368 * t370;
t520 = t368 * t397;
t519 = t370 * t397;
t518 = t397 * t414;
t517 = t397 * t417;
t516 = t411 * t416;
t515 = t411 * t419;
t513 = t414 * t415;
t512 = t414 * t418;
t511 = t415 * t417;
t420 = qJD(3) ^ 2;
t510 = t415 * t420;
t507 = t418 * t420;
t467 = -pkin(8) * t414 - pkin(4);
t474 = pkin(5) * t509;
t506 = -pkin(5) * t458 + qJD(6) * t418 + (t474 + (-qJ(6) + t467) * t415) * qJD(3) + t547;
t400 = pkin(8) * t512;
t475 = pkin(5) * t512;
t505 = (-pkin(5) * t511 - t400) * qJD(4) + (-t475 + (-pkin(8) * t417 + qJ(5)) * t415) * qJD(3) - t546;
t383 = t538 * t417;
t372 = t443 * qJD(2);
t454 = -t414 * t350 + t372 * t417;
t504 = -(-t415 * t536 + t474) * qJD(2) + t454 + qJD(4) * t383;
t503 = -qJ(5) * t487 + (t482 * t415 + t418 * t485) * pkin(8) + t546;
t498 = t417 * t350 + t414 * t372;
t502 = (qJ(5) * t415 - t475) * qJD(2) + t498 + t538 * t485;
t501 = t467 * t487 + t547;
t436 = qJ(6) * t414 - t532;
t428 = t436 * t418;
t500 = qJD(2) * t428 - qJD(4) * t436 + qJD(6) * t417 + t544;
t499 = qJ(5) * t484 - t489 * t532 + t544;
t496 = pkin(4) * t513 + t415 * pkin(8);
t495 = t414 * t381 + t401;
t483 = t415 * MDP(16);
t481 = -qJD(5) - t301;
t294 = t302 - t537;
t479 = -qJD(6) - t294;
t473 = pkin(9) * t518;
t472 = pkin(9) * t517;
t471 = pkin(9) * t487;
t470 = pkin(9) * t482;
t469 = t417 * t515;
t463 = t416 * t491;
t459 = t397 * t485;
t455 = -qJ(5) * t414 - pkin(3);
t452 = t381 * t417 - t400;
t448 = pkin(4) * t555 + pkin(8) * t486 + qJ(5) * t458;
t447 = t415 * t465;
t445 = t368 * t465;
t444 = t370 * t465;
t343 = qJ(5) * t418 - t495;
t378 = -t465 - t535;
t439 = -t378 - t465;
t437 = -t386 + t450;
t296 = pkin(4) * t397 - t481;
t435 = t296 * t417 + t297 * t414;
t434 = qJD(2) * t409 - t397 * t418;
t358 = t418 * t516 + t514;
t330 = t358 * t417 - t414 * t515;
t357 = -t412 * t418 + t415 * t516;
t281 = qJD(6) * t368 + t346 * t536 + t425;
t423 = -qJ(5) * t370 + t339;
t295 = t368 * t536 + t423;
t432 = t281 * t414 + t295 * t484;
t431 = -t281 * t417 + t295 * t485;
t430 = pkin(5) * t345 + t449;
t424 = qJD(3) * (-t439 - t535);
t310 = -t345 - t520;
t422 = -t402 * t536 - t430;
t421 = qJD(2) ^ 2;
t408 = t418 * pkin(4);
t393 = 0.2e1 * t396;
t382 = t538 * t414;
t379 = -pkin(4) * t417 + t455;
t361 = -t417 * t536 + t455;
t354 = -qJ(5) * t511 + t496;
t344 = t408 - t452;
t338 = t415 * t436 + t496;
t329 = t358 * t414 + t469;
t328 = qJD(3) * t358 + t446;
t327 = -qJD(3) * t357 + t418 * t462;
t325 = pkin(4) * t370 + t533;
t322 = -pkin(5) * t513 - t343;
t321 = qJ(6) * t418 + t400 + t408 + (pkin(5) * t415 - t381) * t417;
t311 = t370 * t536 + t533;
t309 = (-qJ(5) * t486 - qJD(5) * t415) * t417 + t448;
t308 = -pkin(4) * t490 - t454;
t307 = -qJ(5) * t490 - t498;
t303 = pkin(4) * t368 + t423;
t292 = qJD(3) * t428 + (qJD(6) * t414 + (qJ(6) * qJD(4) - qJD(5)) * t417) * t415 + t448;
t290 = -qJD(4) * t469 + t327 * t417 - t358 * t485 + t414 * t463;
t289 = qJD(4) * t330 + t327 * t414 - t417 * t463;
t287 = qJD(6) - t297 - t537;
t282 = -t396 - t437;
t279 = qJD(6) * t397 + t422;
t1 = [(-t282 * t330 + t283 * t329 + t284 * t357 + t289 * t296 - t290 * t297 + t303 * t328) * MDP(22) + (t279 * t329 + t280 * t330 + t281 * t357 + t285 * t289 + t287 * t290 + t295 * t328) * MDP(26) + (-t328 * MDP(10) - t327 * MDP(11) + (-t553 * t330 + (-MDP(25) - t543) * t329) * t490) * qJD(3) + ((-MDP(10) * t415 - MDP(11) * t418) * t419 * t478 + (-t419 * MDP(4) + (-MDP(10) * t418 + MDP(11) * t415 - MDP(3)) * t416) * t421) * t411 + t542 * (t289 * t370 - t290 * t368 - t329 * t345 - t330 * t346) + (MDP(17) + t552) * (t289 * t397 + t328 * t368 + t357 * t346) + t553 * (t290 * t397 + t328 * t370 - t357 * t345); 0.2e1 * t402 * t551 - 0.2e1 * t478 * t550 + MDP(7) * t507 - MDP(8) * t510 + (-pkin(8) * t507 + t415 * t424) * MDP(10) + (pkin(8) * t510 + t418 * t424) * MDP(11) + (-t345 * t511 - t370 * t554) * MDP(12) + ((-t368 * t417 - t370 * t414) * t486 + (t522 - t346 * t417 + (t368 * t414 - t370 * t417) * qJD(4)) * t415) * MDP(13) + (t397 * t458 + t345 * t418 + (t370 * t415 + t417 * t434) * qJD(3)) * MDP(14) + (t397 * t457 + t346 * t418 + (-t368 * t415 - t414 * t434) * qJD(3)) * MDP(15) + (-t397 - t489) * qJD(3) * t483 + (t547 * t397 + ((pkin(8) * t368 + t523) * qJD(3) - t449) * t418 + (-t445 + t339 * t484 + pkin(8) * t346 + t525 + (-pkin(8) * t518 + qJD(2) * t452 - t301) * qJD(3)) * t415) * MDP(17) + (t556 * t397 + (t339 * t482 + (qJD(3) * t370 - t459) * pkin(8) + t450) * t418 + (-t444 - t339 * t485 - pkin(8) * t345 + t524 + (-pkin(8) * t517 - qJD(2) * t495 - t302) * qJD(3)) * t415) * MDP(18) + (t343 * t346 - t344 * t345 + t501 * t370 + t503 * t368 + t435 * t486 + (t282 * t414 + t283 * t417 + (-t296 * t414 + t297 * t417) * qJD(4)) * t415) * MDP(19) + (-t309 * t368 - t346 * t354 + (-t303 * t488 - t283) * t418 - t501 * t397 + (t445 - t303 * t484 - t531 + (qJD(2) * t344 + t296) * qJD(3)) * t415) * MDP(20) + (-t309 * t370 + t345 * t354 + (-t303 * t482 + t282) * t418 + t503 * t397 + (t444 + t303 * t485 - t530 + (-qJD(2) * t343 - t297) * qJD(3)) * t415) * MDP(21) + (t282 * t343 + t283 * t344 + t284 * t354 + (t309 - t447) * t303 + t503 * t297 + t501 * t296) * MDP(22) + (-t321 * t345 - t322 * t346 + t506 * t370 - t505 * t368 + (t285 * t417 - t287 * t414) * t486 + (t279 * t417 - t280 * t414 + (-t285 * t414 - t287 * t417) * qJD(4)) * t415) * MDP(23) + (-t292 * t370 + t338 * t345 + (-t295 * t482 - t280) * t418 - t505 * t397 + (t444 + (qJD(2) * t322 + t287) * qJD(3) + t431) * t415) * MDP(24) + (t292 * t368 + t338 * t346 + (t295 * t488 + t279) * t418 + t506 * t397 + (-t445 + (-qJD(2) * t321 - t285) * qJD(3) + t432) * t415) * MDP(25) + (t279 * t321 + t280 * t322 + t281 * t338 + (t292 - t447) * t295 + t505 * t287 + t506 * t285) * MDP(26); (qJD(3) * t351 - t378 * t490 - t320) * MDP(10) + t439 * t489 * MDP(11) + (-t370 * t517 - t522) * MDP(12) + ((-t345 + t520) * t417 + (-t346 + t519) * t414) * MDP(13) + (-t397 * t484 + (t397 * t509 + (-t370 + t488) * t415) * qJD(2)) * MDP(14) + (t459 + (-t397 * t512 + (t368 + t482) * t415) * qJD(2)) * MDP(15) + t397 * qJD(2) * t483 + (-pkin(3) * t346 - t524 + t454 * t397 - t351 * t368 + (t472 + t523) * qJD(4) + (t301 * t415 + (-t339 * t418 - t471) * t414) * qJD(2)) * MDP(17) + (pkin(3) * t345 + t525 - t498 * t397 - t351 * t370 + (t339 * t417 - t473) * qJD(4) + (-t339 * t509 + (t302 - t470) * t415) * qJD(2)) * MDP(18) + (-t307 * t368 - t308 * t370 + (-t282 - t397 * t296 + (qJD(4) * t370 - t346) * pkin(9)) * t417 + ((qJD(4) * t368 - t345) * pkin(9) + t548) * t414) * MDP(19) + (t530 + t308 * t397 - t346 * t379 + t499 * t368 + (-t303 * t414 - t472) * qJD(4) + (-t296 * t415 + (t303 * t418 + t471) * t414) * qJD(2)) * MDP(20) + (-t531 - t307 * t397 + t345 * t379 + t499 * t370 + (-t303 * t417 + t473) * qJD(4) + (t303 * t509 + (t297 + t470) * t415) * qJD(2)) * MDP(21) + (t284 * t379 - t296 * t308 - t297 * t307 - t499 * t303 + (qJD(4) * t435 - t282 * t417 + t283 * t414) * pkin(9)) * MDP(22) + (-t345 * t382 - t346 * t383 + t504 * t370 + t502 * t368 + t549 * t417 + (t287 * t397 + t279) * t414) * MDP(23) + (t345 * t361 + t502 * t397 + t500 * t370 + (t295 * t509 + (qJD(3) * t383 - t287) * t415) * qJD(2) - t432) * MDP(24) + (t346 * t361 + t504 * t397 - t500 * t368 + (-t295 * t512 + (-qJD(3) * t382 + t285) * t415) * qJD(2) + t431) * MDP(25) + (t279 * t382 + t280 * t383 + t281 * t361 + t285 * t504 - t287 * t502 - t295 * t500) * MDP(26) + (-t415 * t551 + t550) * t421; t310 * MDP(14) - MDP(15) * t456 - t450 * MDP(18) + (pkin(4) * t345 - t534) * MDP(19) + (t393 + t437) * MDP(21) + (-pkin(4) * t283 - qJ(5) * t282 - t296 * t302 + t297 * t481 - t303 * t325) * MDP(22) + (t345 * t536 - t534) * MDP(23) + (-0.2e1 * t386 + t393 + t429) * MDP(24) + t430 * MDP(25) + (qJ(5) * t280 - t279 * t536 + t285 * t479 + t287 * t480 - t295 * t311) * MDP(26) + (t301 * MDP(18) + t481 * MDP(21) - t433 * MDP(24) + (-0.2e1 * qJD(6) - t294) * MDP(25) - t543 * t302) * t397 + (-t397 * MDP(15) - t339 * MDP(17) + (-t297 - t302) * MDP(19) + t303 * MDP(20) + t325 * MDP(21) + (t287 + t479) * MDP(23) + t311 * MDP(24) - t295 * MDP(25) + MDP(13) * t370) * t370 + (-MDP(15) * t457 + (-MDP(15) * t512 + (-0.2e1 * pkin(4) * MDP(20) + 0.2e1 * t536 * MDP(25) + MDP(16)) * t415) * qJD(3)) * qJD(2) + (t370 * MDP(12) + t339 * MDP(18) + (t296 + t481) * MDP(19) + t325 * MDP(20) - t303 * MDP(21) + (t285 - t480) * MDP(23) - t295 * MDP(24) - t311 * MDP(25) - MDP(13) * t368) * t368 + t543 * t449; (t303 * t370 + t548) * MDP(22) + (t295 * t370 + (qJD(6) + t287) * t397 + t422) * MDP(26) + t552 * (-t402 + t521) + t542 * t310 + t541 * (-t370 ^ 2 - t395); (-t346 - t519) * MDP(23) + (t402 + t521) * MDP(24) + (-t368 ^ 2 - t395) * MDP(25) + (-t295 * t368 + t549) * MDP(26);];
tauc  = t1;
