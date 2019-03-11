% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RRPRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:48:00
% EndTime: 2019-03-09 09:48:10
% DurationCPUTime: 5.14s
% Computational Cost: add. (7865->450), mult. (19858->590), div. (0->0), fcn. (14409->8), ass. (0->192)
t450 = sin(pkin(9));
t452 = cos(pkin(9));
t454 = sin(qJ(2));
t456 = cos(qJ(2));
t428 = t450 * t456 + t452 * t454;
t449 = sin(pkin(10));
t451 = cos(pkin(10));
t453 = sin(qJ(4));
t455 = cos(qJ(4));
t532 = -t449 * t453 + t451 * t455;
t367 = t532 * t428;
t417 = t428 * qJD(1);
t491 = t455 * qJD(2);
t387 = t417 * t453 - t491;
t389 = qJD(2) * t453 + t417 * t455;
t334 = t451 * t387 + t389 * t449;
t492 = t454 * qJD(1);
t513 = t452 * t456;
t414 = qJD(1) * t513 - t450 * t492;
t409 = qJD(4) - t414;
t540 = t334 * t409;
t427 = t449 * t455 + t451 * t453;
t413 = t427 * qJD(4);
t502 = t427 * t414 - t413;
t494 = qJD(4) * t455;
t495 = qJD(4) * t453;
t501 = t532 * t414 + t449 * t495 - t451 * t494;
t426 = t450 * t454 - t513;
t419 = t426 * qJD(2);
t482 = t428 * t494;
t539 = -t419 * t453 + t482;
t467 = -t387 * t449 + t451 * t389;
t538 = t467 ^ 2;
t489 = qJD(1) * qJD(2);
t537 = -0.2e1 * t489;
t536 = MDP(4) * t454;
t535 = MDP(5) * (t454 ^ 2 - t456 ^ 2);
t530 = -qJ(3) - pkin(7);
t434 = t530 * t456;
t431 = qJD(1) * t434;
t420 = t450 * t431;
t433 = t530 * t454;
t430 = qJD(1) * t433;
t529 = qJD(2) * pkin(2);
t423 = t430 + t529;
t377 = t423 * t452 + t420;
t370 = -qJD(2) * pkin(3) - t377;
t332 = pkin(4) * t387 + qJD(5) + t370;
t303 = pkin(5) * t334 - qJ(6) * t467 + t332;
t534 = t303 * t467;
t474 = t409 * t455;
t362 = pkin(2) * t492 + pkin(3) * t417 - pkin(8) * t414;
t358 = t455 * t362;
t380 = t430 * t452 + t420;
t314 = -qJ(5) * t414 * t455 + pkin(4) * t417 - t380 * t453 + t358;
t499 = t453 * t362 + t455 * t380;
t520 = t414 * t453;
t323 = -qJ(5) * t520 + t499;
t442 = pkin(2) * t450 + pkin(8);
t508 = qJ(5) + t442;
t475 = qJD(4) * t508;
t394 = qJD(5) * t455 - t453 * t475;
t459 = -qJD(5) * t453 - t455 * t475;
t503 = (-t314 + t459) * t451 + (t323 - t394) * t449;
t514 = t452 * t431;
t379 = t430 * t450 - t514;
t533 = -t379 + (t495 - t520) * pkin(4);
t531 = MDP(20) + MDP(23);
t483 = -pkin(2) * t456 - pkin(1);
t470 = t483 * qJD(1);
t432 = qJD(3) + t470;
t352 = -pkin(3) * t414 - pkin(8) * t417 + t432;
t378 = t450 * t423 - t514;
t371 = qJD(2) * pkin(8) + t378;
t326 = t352 * t453 + t371 * t455;
t316 = -qJ(5) * t387 + t326;
t311 = t451 * t316;
t325 = t455 * t352 - t371 * t453;
t315 = -qJ(5) * t389 + t325;
t291 = t315 * t449 + t311;
t528 = t291 * t467;
t527 = t316 * t449;
t480 = t456 * t489;
t481 = t454 * t489;
t405 = -t450 * t481 + t452 * t480;
t497 = qJD(4) * t491 + t455 * t405;
t348 = -t417 * t495 + t497;
t526 = t348 * t453;
t525 = t387 * t409;
t524 = t387 * t417;
t523 = t389 * t409;
t522 = t389 * t417;
t521 = t405 * t453;
t518 = t428 * t453;
t517 = t428 * t455;
t415 = t428 * qJD(2);
t404 = qJD(1) * t415;
t512 = t453 * t404;
t457 = qJD(2) ^ 2;
t511 = t454 * t457;
t385 = t433 * t450 - t434 * t452;
t381 = t455 * t385;
t395 = t455 * t404;
t510 = t456 * t457;
t458 = qJD(1) ^ 2;
t509 = t456 * t458;
t438 = pkin(2) * t481;
t351 = pkin(3) * t404 - pkin(8) * t405 + t438;
t344 = t455 * t351;
t479 = qJD(2) * t530;
t411 = qJD(3) * t456 + t454 * t479;
t397 = t411 * qJD(1);
t412 = -qJD(3) * t454 + t456 * t479;
t398 = t412 * qJD(1);
t347 = t397 * t452 + t398 * t450;
t476 = -t347 * t453 + t344;
t281 = pkin(4) * t404 - qJ(5) * t348 - t326 * qJD(4) - qJD(5) * t389 + t476;
t349 = qJD(4) * t389 + t521;
t485 = t455 * t347 + t453 * t351 + t352 * t494;
t462 = -t371 * t495 + t485;
t286 = -qJ(5) * t349 - qJD(5) * t387 + t462;
t507 = -t451 * t281 + t449 * t286;
t272 = t449 * t281 + t451 * t286;
t296 = t449 * t314 + t451 * t323;
t289 = qJ(6) * t417 + t296;
t342 = t451 * t394 + t449 * t459;
t506 = -t289 + t342;
t505 = -pkin(5) * t417 + t503;
t488 = t454 * t529;
t363 = pkin(3) * t415 + pkin(8) * t419 + t488;
t359 = t455 * t363;
t361 = t411 * t452 + t412 * t450;
t376 = pkin(3) * t426 - pkin(8) * t428 + t483;
t466 = qJ(5) * t419 - qJD(5) * t428;
t294 = pkin(4) * t415 - t361 * t453 + t359 + t466 * t455 + (-t381 + (qJ(5) * t428 - t376) * t453) * qJD(4);
t484 = t455 * t361 + t453 * t363 + t376 * t494;
t298 = -qJ(5) * t482 + (-qJD(4) * t385 + t466) * t453 + t484;
t276 = t449 * t294 + t451 * t298;
t504 = pkin(5) * t502 - qJ(6) * t501 + qJD(6) * t427 - t533;
t308 = pkin(4) * t409 + t315;
t288 = t449 * t308 + t311;
t369 = t455 * t376;
t322 = pkin(4) * t426 - qJ(5) * t517 - t385 * t453 + t369;
t498 = t453 * t376 + t381;
t328 = -qJ(5) * t518 + t498;
t302 = t449 * t322 + t451 * t328;
t500 = t409 * t520 + t395;
t346 = t450 * t397 - t452 * t398;
t493 = qJD(6) * t409;
t292 = t315 * t451 - t527;
t490 = qJD(6) - t292;
t486 = t404 * qJ(6) + t272;
t444 = -pkin(2) * t452 - pkin(3);
t478 = t508 * t453;
t477 = pkin(1) * t537;
t318 = t348 * t449 + t451 * t349;
t360 = t411 * t450 - t452 * t412;
t384 = -t452 * t433 - t434 * t450;
t324 = pkin(4) * t349 + t346;
t270 = -pkin(5) * t404 + t507;
t319 = t348 * t451 - t349 * t449;
t424 = t508 * t455;
t372 = t424 * t449 + t451 * t478;
t373 = t451 * t424 - t449 * t478;
t473 = -t373 * t318 + t319 * t372 - t342 * t334;
t471 = pkin(4) * t518 + t384;
t275 = t294 * t451 - t298 * t449;
t287 = t308 * t451 - t527;
t301 = t322 * t451 - t328 * t449;
t468 = t346 * t428 - t385 * t404;
t465 = -pkin(4) * t455 + t444;
t464 = pkin(4) * t539 + t360;
t463 = -t419 * t455 - t428 * t495;
t461 = t409 * t370 - t442 * t404;
t277 = pkin(5) * t318 - qJ(6) * t319 - qJD(6) * t467 + t324;
t443 = -pkin(4) * t451 - pkin(5);
t439 = pkin(4) * t449 + qJ(6);
t366 = t427 * t428;
t364 = -pkin(5) * t532 - qJ(6) * t427 + t465;
t330 = t428 * t413 + t419 * t532;
t329 = -qJD(4) * t367 + t427 * t419;
t309 = pkin(5) * t366 - qJ(6) * t367 + t471;
t305 = pkin(4) * t389 + pkin(5) * t467 + qJ(6) * t334;
t300 = -pkin(5) * t426 - t301;
t299 = qJ(6) * t426 + t302;
t283 = qJ(6) * t409 + t288;
t282 = -pkin(5) * t409 + qJD(6) - t287;
t280 = -pkin(5) * t329 + qJ(6) * t330 - qJD(6) * t367 + t464;
t274 = -pkin(5) * t415 - t275;
t273 = qJ(6) * t415 + qJD(6) * t426 + t276;
t269 = t486 + t493;
t1 = [0.2e1 * t480 * t536 + t535 * t537 + MDP(6) * t510 - MDP(7) * t511 + (-pkin(7) * t510 + t454 * t477) * MDP(9) + (pkin(7) * t511 + t456 * t477) * MDP(10) + (-t347 * t426 + t360 * t417 + t361 * t414 + t377 * t419 - t378 * t415 + t384 * t405 + t468) * MDP(11) + (t346 * t384 + t347 * t385 - t377 * t360 + t378 * t361 + (t432 + t470) * t488) * MDP(12) + (t348 * t517 + t463 * t389) * MDP(13) + (-(-t387 * t455 - t389 * t453) * t419 + (-t526 - t349 * t455 + (t387 * t453 - t389 * t455) * qJD(4)) * t428) * MDP(14) + (t348 * t426 + t389 * t415 + t428 * t395 + t463 * t409) * MDP(15) + (-t349 * t426 - t387 * t415 - t409 * t539 - t428 * t512) * MDP(16) + (t404 * t426 + t409 * t415) * MDP(17) + ((-t385 * t494 + t359) * t409 + t369 * t404 + (-t371 * t494 + t344) * t426 + t325 * t415 + t360 * t387 + t384 * t349 + t370 * t482 + ((-qJD(4) * t376 - t361) * t409 + (-qJD(4) * t352 - t347) * t426 - t370 * t419 + t468) * t453) * MDP(18) + (-(-t385 * t495 + t484) * t409 - t498 * t404 - t462 * t426 - t326 * t415 + t360 * t389 + t384 * t348 + t346 * t517 + t463 * t370) * MDP(19) + (-t272 * t366 - t275 * t467 - t276 * t334 + t287 * t330 + t288 * t329 - t301 * t319 - t302 * t318 + t367 * t507) * MDP(20) + (t272 * t302 + t287 * t275 + t288 * t276 - t301 * t507 + t324 * t471 + t332 * t464) * MDP(21) + (-t270 * t426 - t274 * t409 + t277 * t366 + t280 * t334 - t282 * t415 - t300 * t404 - t303 * t329 + t309 * t318) * MDP(22) + (-t269 * t366 + t270 * t367 - t273 * t334 + t274 * t467 - t282 * t330 + t283 * t329 - t299 * t318 + t300 * t319) * MDP(23) + (t269 * t426 + t273 * t409 - t277 * t367 - t280 * t467 + t283 * t415 + t299 * t404 + t303 * t330 - t309 * t319) * MDP(24) + (t269 * t299 + t270 * t300 + t273 * t283 + t274 * t282 + t277 * t309 + t280 * t303) * MDP(25); -t509 * t536 + t458 * t535 + ((t378 - t379) * t417 + (t377 - t380) * t414 + (-t404 * t450 - t405 * t452) * pkin(2)) * MDP(11) + (t377 * t379 - t378 * t380 + (-t346 * t452 + t347 * t450 - t432 * t492) * pkin(2)) * MDP(12) + (t389 * t474 + t526) * MDP(13) + ((t348 - t525) * t455 + (-t349 - t523) * t453) * MDP(14) + (t409 * t474 + t512 - t522) * MDP(15) + (-t409 * t495 + t500 + t524) * MDP(16) - t409 * t417 * MDP(17) + (-t325 * t417 - t346 * t455 + t444 * t349 - t379 * t387 + (-t442 * t494 - t358) * t409 + (t380 * t409 + t461) * t453) * MDP(18) + (t326 * t417 + t346 * t453 + t444 * t348 - t379 * t389 + (t442 * t495 + t499) * t409 + t461 * t455) * MDP(19) + (t272 * t532 + t501 * t287 + t502 * t288 + t296 * t334 + t427 * t507 - t467 * t503 + t473) * MDP(20) + (t272 * t373 + t507 * t372 + t324 * t465 + t533 * t332 + (t342 - t296) * t288 + t503 * t287) * MDP(21) + (-t277 * t532 + t282 * t417 - t502 * t303 + t318 * t364 - t504 * t334 - t372 * t404 + t505 * t409) * MDP(22) + (t269 * t532 + t270 * t427 - t501 * t282 + t502 * t283 + t289 * t334 - t467 * t505 + t473) * MDP(23) + (-t277 * t427 - t283 * t417 + t501 * t303 - t319 * t364 + t373 * t404 + t506 * t409 + t467 * t504) * MDP(24) + (t269 * t373 + t270 * t372 + t277 * t364 - t505 * t282 + t506 * t283 - t504 * t303) * MDP(25) + (t458 * t454 * MDP(9) + MDP(10) * t509) * pkin(1); (-t414 ^ 2 - t417 ^ 2) * MDP(11) + (t377 * t417 - t378 * t414 + t438) * MDP(12) + (t500 - t524) * MDP(18) + (-t512 - t522) * MDP(19) + (t272 * t427 + t502 * t287 - t501 * t288 - t332 * t417 - t507 * t532) * MDP(21) + (-t334 * t417 + t404 * t532) * MDP(22) + t427 * t404 * MDP(24) + (t269 * t427 - t270 * t532 - t502 * t282 - t501 * t283 - t303 * t417) * MDP(25) + (t417 * MDP(24) - t502 * t531) * t467 + (-MDP(18) * t495 - MDP(19) * t474 + t502 * MDP(22) - t501 * MDP(24)) * t409 + t531 * (-t427 * t318 - t319 * t532 + t334 * t501); t389 * t387 * MDP(13) + (-t387 ^ 2 + t389 ^ 2) * MDP(14) + (t497 + t525) * MDP(15) + (-t521 + t523) * MDP(16) + (t326 * t409 - t370 * t389 + t476) * MDP(18) + (t325 * t409 + t370 * t387 - t485) * MDP(19) + (t288 * t467 - t528) * MDP(20) + (t287 * t291 - t288 * t292) * MDP(21) + (t291 * t409 - t507 - t534) * MDP(22) + (t283 * t467 - t318 * t439 + t319 * t443 - t528) * MDP(23) + (-t292 * t409 + t305 * t467 + t486 + 0.2e1 * t493) * MDP(24) + (t269 * t439 + t270 * t443 - t282 * t291 + t490 * t283 - t303 * t305) * MDP(25) + (MDP(17) + (pkin(5) - t443) * MDP(22) + t439 * MDP(24)) * t404 + ((-MDP(16) * t417 - MDP(18) * t371) * t455 + (-t417 * MDP(15) - MDP(16) * qJD(2) - MDP(18) * t352 + MDP(19) * t371) * t453) * qJD(4) + ((-t318 * t449 - t319 * t451) * MDP(20) + (t272 * t449 - t332 * t389 - t451 * t507) * MDP(21)) * pkin(4) + ((-t287 + t292) * MDP(20) - t305 * MDP(22) + (t282 - t490) * MDP(23) - t303 * MDP(24)) * t334; (t287 * t467 + t288 * t334 + t324) * MDP(21) + (t409 * t467 + t318) * MDP(22) + (-t319 + t540) * MDP(24) + (-t282 * t467 + t283 * t334 + t277) * MDP(25) + t531 * (-t334 ^ 2 - t538); (-qJD(2) * t417 + t334 * t467) * MDP(22) + (t319 + t540) * MDP(23) + (-t409 ^ 2 - t538) * MDP(24) + (-t283 * t409 + t270 + t534) * MDP(25);];
tauc  = t1;
