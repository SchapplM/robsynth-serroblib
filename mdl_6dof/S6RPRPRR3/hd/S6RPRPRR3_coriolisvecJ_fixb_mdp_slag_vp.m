% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:36
% EndTime: 2019-03-09 03:42:47
% DurationCPUTime: 6.39s
% Computational Cost: add. (4189->436), mult. (10106->602), div. (0->0), fcn. (7441->10), ass. (0->188)
t479 = cos(qJ(3));
t534 = qJD(1) * t479;
t460 = -qJD(5) + t534;
t456 = -qJD(6) + t460;
t470 = sin(pkin(11));
t472 = cos(pkin(11));
t522 = t472 * qJD(3);
t476 = sin(qJ(3));
t535 = qJD(1) * t476;
t432 = t470 * t535 - t522;
t532 = qJD(3) * t470;
t434 = t472 * t535 + t532;
t475 = sin(qJ(5));
t478 = cos(qJ(5));
t380 = t478 * t432 + t434 * t475;
t477 = cos(qJ(6));
t379 = t432 * t475 - t434 * t478;
t474 = sin(qJ(6));
t554 = t379 * t474;
t573 = -t477 * t380 + t554;
t572 = t456 * t573;
t500 = pkin(3) * t476 - qJ(4) * t479;
t444 = t500 * qJD(1);
t462 = sin(pkin(10)) * pkin(1) + pkin(7);
t451 = t462 * qJD(1);
t562 = qJD(2) * t479 - t476 * t451;
t372 = t470 * t444 + t472 * t562;
t518 = t470 * t534;
t354 = -pkin(8) * t518 + t372;
t574 = qJD(4) * t472 - t354;
t571 = t379 * t460;
t570 = t380 * t460;
t495 = t379 * t477 + t380 * t474;
t569 = t456 * t495;
t547 = t478 * t472;
t553 = t470 * t475;
t439 = -t547 + t553;
t489 = t439 * t479;
t539 = qJD(1) * t489 - t439 * qJD(5);
t440 = t470 * t478 + t472 * t475;
t490 = t440 * t479;
t538 = -qJD(1) * t490 + t440 * qJD(5);
t467 = t476 * qJD(2);
t418 = t479 * t451 + t467;
t402 = qJD(3) * qJ(4) + t418;
t463 = -cos(pkin(10)) * pkin(1) - pkin(2);
t428 = -pkin(3) * t479 - qJ(4) * t476 + t463;
t405 = t428 * qJD(1);
t344 = -t402 * t470 + t472 * t405;
t324 = -pkin(4) * t534 - pkin(8) * t434 + t344;
t345 = t472 * t402 + t470 * t405;
t327 = -pkin(8) * t432 + t345;
t303 = t324 * t475 + t327 * t478;
t300 = -pkin(9) * t380 + t303;
t525 = qJD(6) * t474;
t298 = t300 * t525;
t511 = -qJD(3) * pkin(3) + qJD(4);
t399 = -t562 + t511;
t373 = pkin(4) * t432 + t399;
t323 = pkin(5) * t380 + t373;
t567 = -t323 * t573 + t298;
t521 = qJD(1) * qJD(3);
t515 = t479 * t521;
t503 = t470 * t515;
t526 = qJD(5) * t478;
t335 = -t432 * t526 + t515 * t547 + (-qJD(5) * t434 - t503) * t475;
t396 = (qJD(4) + t562) * qJD(3);
t421 = t500 * qJD(3) - qJD(4) * t476;
t406 = t421 * qJD(1);
t339 = -t396 * t470 + t472 * t406;
t549 = t472 * t479;
t491 = pkin(4) * t476 - pkin(8) * t549;
t486 = t491 * qJD(3);
t325 = qJD(1) * t486 + t339;
t340 = t472 * t396 + t470 * t406;
t328 = -pkin(8) * t503 + t340;
t509 = t478 * t325 - t328 * t475;
t482 = -t303 * qJD(5) + t509;
t516 = t476 * t521;
t291 = pkin(5) * t516 - pkin(9) * t335 + t482;
t485 = qJD(3) * t490;
t560 = t379 * qJD(5);
t336 = qJD(1) * t485 - t560;
t528 = qJD(5) * t475;
t488 = t324 * t526 + t475 * t325 - t327 * t528 + t478 * t328;
t292 = -pkin(9) * t336 + t488;
t510 = t477 * t291 - t474 * t292;
t566 = t323 * t495 + t510;
t531 = qJD(3) * t476;
t514 = MDP(27) * t531;
t565 = qJD(1) * t514 + (t495 ^ 2 - t573 ^ 2) * MDP(24) + t573 * t495 * MDP(23);
t564 = MDP(5) * t476;
t468 = t476 ^ 2;
t563 = (-t479 ^ 2 + t468) * MDP(6);
t416 = t472 * t428;
t550 = t472 * t476;
t362 = -pkin(8) * t550 + t416 + (-t462 * t470 - pkin(4)) * t479;
t443 = t462 * t549;
t387 = t470 * t428 + t443;
t552 = t470 * t476;
t370 = -pkin(8) * t552 + t387;
t540 = t475 * t362 + t478 * t370;
t559 = pkin(8) + qJ(4);
t448 = t559 * t470;
t449 = t559 * t472;
t537 = -t475 * t448 + t478 * t449;
t371 = t472 * t444 - t470 * t562;
t341 = t491 * qJD(1) + t371;
t492 = qJD(4) * t470 + qJD(5) * t449;
t561 = -t448 * t526 + t574 * t478 + (-t341 - t492) * t475;
t508 = t474 * t335 + t477 * t336;
t296 = -t495 * qJD(6) + t508;
t302 = t478 * t324 - t327 * t475;
t299 = pkin(9) * t379 + t302;
t297 = -pkin(5) * t460 + t299;
t558 = t297 * t477;
t557 = t300 * t477;
t413 = t440 * t476;
t414 = t439 * t476;
t359 = t477 * t413 - t414 * t474;
t527 = qJD(5) * t476;
t366 = -qJD(3) * t489 - t440 * t527;
t367 = t526 * t550 - t527 * t553 + t485;
t304 = -t359 * qJD(6) + t366 * t477 - t367 * t474;
t556 = t304 * t456;
t555 = t366 * t460;
t551 = t470 * t479;
t453 = t476 * t462;
t480 = qJD(3) ^ 2;
t548 = t476 * t480;
t546 = t479 * t480;
t360 = -t413 * t474 - t414 * t477;
t305 = t360 * qJD(6) + t366 * t474 + t477 * t367;
t545 = t305 * t456 - t359 * t516;
t384 = t477 * t439 + t440 * t474;
t544 = -t384 * qJD(6) - t474 * t538 + t477 * t539;
t385 = -t439 * t474 + t440 * t477;
t543 = t385 * qJD(6) + t474 * t539 + t477 * t538;
t541 = t367 * t460 - t413 * t516;
t530 = qJD(3) * t479;
t407 = qJD(3) * t467 + t451 * t530;
t517 = t462 * t531;
t377 = t472 * t421 + t470 * t517;
t411 = (pkin(4) * t470 + t462) * t530;
t420 = pkin(4) * t552 + t453;
t452 = qJD(1) * t463;
t524 = qJD(6) * t477;
t520 = t462 * t551;
t519 = t477 * t335 - t474 * t336 - t380 * t524;
t388 = pkin(4) * t503 + t407;
t390 = pkin(4) * t518 + t418;
t464 = -pkin(4) * t472 - pkin(3);
t513 = MDP(20) * t535;
t512 = pkin(5) * t538 - t390;
t352 = t486 + t377;
t409 = t470 * t421;
t361 = t409 + (-pkin(8) * t551 - t472 * t453) * qJD(3);
t507 = t478 * t352 - t361 * t475;
t506 = t478 * t362 - t370 * t475;
t505 = -t478 * t448 - t449 * t475;
t504 = qJD(6) * t297 + t292;
t338 = t478 * t341;
t364 = -pkin(9) * t439 + t537;
t502 = pkin(5) * t535 + pkin(9) * t539 + t440 * qJD(4) + t537 * qJD(5) + qJD(6) * t364 - t354 * t475 + t338;
t363 = -pkin(9) * t440 + t505;
t501 = -pkin(9) * t538 + qJD(6) * t363 + t561;
t289 = t297 * t474 + t557;
t308 = -pkin(5) * t479 + pkin(9) * t414 + t506;
t309 = -pkin(9) * t413 + t540;
t499 = t308 * t474 + t309 * t477;
t498 = -t339 * t470 + t340 * t472;
t497 = -t344 * t470 + t345 * t472;
t493 = 0.2e1 * qJD(3) * t452;
t487 = t475 * t352 + t478 * t361 + t362 * t526 - t370 * t528;
t295 = t379 * t525 + t519;
t483 = -qJ(4) * t531 + (-t399 + t511) * t479;
t408 = pkin(5) * t439 + t464;
t386 = t416 - t520;
t378 = -t472 * t517 + t409;
t376 = pkin(5) * t413 + t420;
t368 = t379 * t531;
t332 = pkin(5) * t367 + t411;
t311 = t495 * t531;
t310 = pkin(5) * t336 + t388;
t294 = -pkin(9) * t367 + t487;
t293 = pkin(5) * t531 - pkin(9) * t366 - qJD(5) * t540 + t507;
t288 = -t300 * t474 + t558;
t1 = [-0.2e1 * t521 * t563 + (-t462 * t546 + t476 * t493) * MDP(10) + (t462 * t548 + t479 * t493) * MDP(11) + (t407 * t552 + (-qJD(1) * t377 - t339) * t479 + ((t399 * t470 + t432 * t462) * t479 + (t344 + (t386 + t520) * qJD(1)) * t476) * qJD(3)) * MDP(12) + (t407 * t550 + (qJD(1) * t378 + t340) * t479 + ((t399 * t472 + t434 * t462) * t479 + (-t345 + (-t387 + t443) * qJD(1)) * t476) * qJD(3)) * MDP(13) + (-t377 * t434 - t378 * t432 + (-t339 * t472 - t340 * t470) * t476 + (-t344 * t472 - t345 * t470 + (-t386 * t472 - t387 * t470) * qJD(1)) * t530) * MDP(14) + (t339 * t386 + t340 * t387 + t344 * t377 + t345 * t378 + (t399 * t530 + t407 * t476) * t462) * MDP(15) + (-t335 * t414 - t366 * t379) * MDP(16) + (-t335 * t413 + t336 * t414 - t366 * t380 + t367 * t379) * MDP(17) + (-t335 * t479 - t414 * t516 - t368 - t555) * MDP(18) + (t336 * t479 - t380 * t531 + t541) * MDP(19) + (-t460 - t534) * MDP(20) * t531 + (-t507 * t460 - t509 * t479 + t411 * t380 + t420 * t336 + t388 * t413 + t373 * t367 + (t303 * t479 + t460 * t540) * qJD(5) + (t506 * qJD(1) + t302) * t531) * MDP(21) + (t487 * t460 + t488 * t479 - t411 * t379 + t420 * t335 - t388 * t414 + t373 * t366 + (-t540 * qJD(1) - t303) * t531) * MDP(22) + (t295 * t360 - t304 * t495) * MDP(23) + (-t295 * t359 - t296 * t360 + t304 * t573 + t305 * t495) * MDP(24) + (-t295 * t479 + t360 * t516 - t311 - t556) * MDP(25) + (t296 * t479 + t531 * t573 + t545) * MDP(26) + (-t456 - t534) * t514 + (-(t293 * t477 - t294 * t474) * t456 - t510 * t479 - t332 * t573 + t376 * t296 + t310 * t359 + t323 * t305 + (t289 * t479 + t499 * t456) * qJD(6) + ((t308 * t477 - t309 * t474) * qJD(1) + t288) * t531) * MDP(28) + (t376 * t295 - t298 * t479 + t323 * t304 + t310 * t360 - t332 * t495 + ((-qJD(6) * t309 + t293) * t456 + t291 * t479) * t474 + ((qJD(6) * t308 + t294) * t456 + t504 * t479) * t477 + (-t499 * qJD(1) - t289) * t531) * MDP(29) + MDP(7) * t546 - MDP(8) * t548 + 0.2e1 * t515 * t564; t541 * MDP(21) + (-t368 + t555) * MDP(22) + t545 * MDP(28) + (-t311 + t556) * MDP(29) + (-t480 * MDP(11) - t407 * MDP(15) - t336 * MDP(21) - t335 * MDP(22) - t296 * MDP(28) - t295 * MDP(29)) * t479 + (-t480 * MDP(10) + t498 * MDP(15)) * t476 + ((-t470 * MDP(12) - t472 * MDP(13)) * t468 * qJD(1) + ((-t432 * t472 + t434 * t470) * MDP(14) + t497 * MDP(15)) * t479 + (t432 * MDP(12) + t434 * MDP(13) + t399 * MDP(15) + t380 * MDP(21) - t573 * MDP(28) + (t414 * MDP(22) - t360 * MDP(29)) * qJD(1)) * t476) * qJD(3); (qJD(3) * t418 - t452 * t535 - t407) * MDP(10) - t452 * t534 * MDP(11) + (-t407 * t472 - t418 * t432 + (-t344 * t476 + t371 * t479 + t483 * t470) * qJD(1)) * MDP(12) + (t407 * t470 - t418 * t434 + (t345 * t476 - t372 * t479 + t483 * t472) * qJD(1)) * MDP(13) + (t371 * t434 + t372 * t432 + (-qJD(4) * t432 + t344 * t534 + t340) * t472 + (qJD(4) * t434 + t345 * t534 - t339) * t470) * MDP(14) + (-pkin(3) * t407 + t498 * qJ(4) + t497 * qJD(4) - t344 * t371 - t345 * t372 - t399 * t418) * MDP(15) + (t335 * t440 - t379 * t539) * MDP(16) + (-t335 * t439 - t336 * t440 + t379 * t538 - t539 * t380) * MDP(17) + (-t539 * t460 + (qJD(3) * t440 + t379) * t535) * MDP(18) + (t538 * t460 + (-qJD(3) * t439 + t380) * t535) * MDP(19) + t460 * t513 + (t464 * t336 - t390 * t380 + t388 * t439 + (t338 + t492 * t478 + (-qJD(5) * t448 + t574) * t475) * t460 + t538 * t373 + (t505 * qJD(3) - t302) * t535) * MDP(21) + (t464 * t335 + t390 * t379 + t388 * t440 + t561 * t460 + t539 * t373 + (-t537 * qJD(3) + t303) * t535) * MDP(22) + (t295 * t385 - t495 * t544) * MDP(23) + (-t295 * t384 - t296 * t385 + t495 * t543 + t544 * t573) * MDP(24) + (-t544 * t456 + (qJD(3) * t385 + t495) * t535) * MDP(25) + (t543 * t456 + (-qJD(3) * t384 - t573) * t535) * MDP(26) + t456 * MDP(27) * t535 + (t408 * t296 + t310 * t384 + (t501 * t474 + t502 * t477) * t456 + t543 * t323 - t512 * t573 + ((t363 * t477 - t364 * t474) * qJD(3) - t288) * t535) * MDP(28) + (t408 * t295 + t310 * t385 + (-t502 * t474 + t501 * t477) * t456 + t544 * t323 - t512 * t495 + (-(t363 * t474 + t364 * t477) * qJD(3) + t289) * t535) * MDP(29) + (-t479 * t564 + t563) * qJD(1) ^ 2; (-t432 ^ 2 - t434 ^ 2) * MDP(14) + (t344 * t434 + t345 * t432 + t407) * MDP(15) + (t336 + t571) * MDP(21) + (t335 + t570) * MDP(22) + (t296 + t569) * MDP(28) + (t295 - t572) * MDP(29) + ((-t434 + t532) * MDP(12) + (t432 + t522) * MDP(13)) * t534; -t379 * t380 * MDP(16) + (t379 ^ 2 - t380 ^ 2) * MDP(17) + (t335 - t570) * MDP(18) + (-t440 * t515 + t560 + t571) * MDP(19) + qJD(3) * t513 + (-t303 * t460 + t373 * t379 + t482) * MDP(21) + (-t302 * t460 + t373 * t380 - t488) * MDP(22) + (t295 + t572) * MDP(25) + (-t296 + t569) * MDP(26) + ((-t299 * t474 - t557) * t456 - t289 * qJD(6) + (-t379 * t573 + t456 * t525 + t477 * t516) * pkin(5) + t566) * MDP(28) + ((t300 * t456 - t291) * t474 + (-t299 * t456 - t504) * t477 + (-t379 * t495 + t456 * t524 - t474 * t516) * pkin(5) + t567) * MDP(29) + t565; (t519 + t572) * MDP(25) + (-t508 + t569) * MDP(26) + (-t289 * t456 + t566) * MDP(28) + (-t288 * t456 - t474 * t291 - t477 * t292 + t567) * MDP(29) + (MDP(25) * t554 + t495 * MDP(26) - t289 * MDP(28) - MDP(29) * t558) * qJD(6) + t565;];
tauc  = t1;
