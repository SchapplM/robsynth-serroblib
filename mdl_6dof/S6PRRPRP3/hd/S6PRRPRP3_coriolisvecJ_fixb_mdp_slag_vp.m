% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:30
% EndTime: 2019-03-08 21:38:40
% DurationCPUTime: 5.26s
% Computational Cost: add. (4258->441), mult. (10766->611), div. (0->0), fcn. (8026->10), ass. (0->182)
t441 = sin(qJ(3));
t444 = cos(qJ(3));
t467 = pkin(3) * t441 - qJ(4) * t444;
t388 = qJD(3) * t467 - qJD(4) * t441;
t436 = sin(pkin(11));
t503 = qJD(3) * t441;
t493 = pkin(8) * t503;
t421 = t436 * t493;
t438 = cos(pkin(11));
t442 = sin(qJ(2));
t437 = sin(pkin(6));
t510 = qJD(1) * t437;
t445 = cos(qJ(2));
t525 = t444 * t445;
t516 = t438 * t388 - (-t436 * t525 + t438 * t442) * t510 + t421;
t550 = t436 * t388 - (t436 * t442 + t438 * t525) * t510;
t528 = t438 * t444;
t461 = pkin(4) * t441 - pkin(9) * t528;
t454 = t461 * qJD(3);
t549 = t454 + t516;
t529 = t438 * t441;
t532 = t436 * t444;
t548 = (-pkin(8) * t529 - pkin(9) * t532) * qJD(3) + t550;
t507 = qJD(2) * t441;
t402 = qJD(3) * t438 - t436 * t507;
t403 = qJD(3) * t436 + t438 * t507;
t440 = sin(qJ(5));
t443 = cos(qJ(5));
t348 = -t443 * t402 + t403 * t440;
t506 = qJD(2) * t444;
t428 = -qJD(5) + t506;
t547 = t348 * t428;
t526 = t443 * t438;
t405 = t436 * t440 - t526;
t500 = qJD(5) * t443;
t501 = qJD(5) * t440;
t539 = -t436 * t501 + t438 * t500;
t513 = t405 * t506 + t539;
t406 = t436 * t443 + t438 * t440;
t391 = t406 * qJD(5);
t457 = t406 * t444;
t512 = -qJD(2) * t457 + t391;
t546 = MDP(5) * t441;
t545 = MDP(6) * (t441 ^ 2 - t444 ^ 2);
t413 = -pkin(3) * t444 - qJ(4) * t441 - pkin(2);
t401 = t438 * t413;
t355 = -pkin(9) * t529 + t401 + (-pkin(8) * t436 - pkin(4)) * t444;
t380 = pkin(8) * t528 + t436 * t413;
t533 = t436 * t441;
t366 = -pkin(9) * t533 + t380;
t517 = t440 * t355 + t443 * t366;
t544 = t517 * qJD(5) + t548 * t440 - t549 * t443;
t543 = t355 * t500 - t366 * t501 + t549 * t440 + t548 * t443;
t407 = t467 * qJD(2);
t491 = t442 * t510;
t411 = qJD(2) * pkin(8) + t491;
t439 = cos(pkin(6));
t509 = qJD(1) * t439;
t540 = -t441 * t411 + t444 * t509;
t335 = t438 * t407 - t436 * t540;
t318 = qJD(2) * t461 + t335;
t336 = t436 * t407 + t438 * t540;
t488 = t436 * t506;
t325 = -pkin(9) * t488 + t336;
t538 = pkin(9) + qJ(4);
t416 = t538 * t436;
t417 = t538 * t438;
t462 = -t416 * t443 - t417 * t440;
t542 = -qJD(4) * t405 + qJD(5) * t462 - t440 * t318 - t443 * t325;
t368 = -t416 * t440 + t417 * t443;
t541 = qJD(4) * t406 + qJD(5) * t368 + t318 * t443 - t325 * t440;
t495 = MDP(21) + MDP(23);
t494 = MDP(22) - MDP(25);
t496 = qJD(2) * qJD(3);
t481 = t444 * t496;
t471 = t436 * t481;
t502 = qJD(3) * t444;
t473 = t502 * t526;
t514 = qJD(2) * t473 + t402 * t500;
t316 = t440 * (qJD(5) * t403 + t471) - t514;
t537 = qJD(2) * pkin(2);
t424 = t441 * t509;
t377 = t444 * t411 + t424;
t373 = qJD(3) * qJ(4) + t377;
t490 = t445 * t510;
t378 = qJD(2) * t413 - t490;
t320 = t438 * t373 + t436 * t378;
t311 = pkin(9) * t402 + t320;
t536 = t311 * t440;
t508 = qJD(2) * t437;
t486 = t445 * t508;
t472 = qJD(1) * t486;
t344 = qJD(3) * t424 + t411 * t502 + t441 * t472;
t535 = t344 * t436;
t534 = t344 * t438;
t531 = t437 * t442;
t530 = t437 * t445;
t446 = qJD(3) ^ 2;
t527 = t441 * t446;
t524 = t444 * t446;
t523 = qJ(6) * t503 - qJD(6) * t444 + t543;
t522 = -pkin(5) * t503 + t544;
t521 = qJ(6) * t507 - t542;
t520 = pkin(5) * t507 + t541;
t357 = pkin(4) * t488 + t377;
t519 = -t512 * pkin(5) + t513 * qJ(6) + qJD(6) * t406 + t357;
t341 = t444 * t472 + (qJD(4) + t540) * qJD(3);
t358 = (t388 + t491) * qJD(2);
t303 = t438 * t341 + t436 * t358;
t478 = t438 * t493;
t515 = -t478 + t550;
t485 = t436 * t502;
t397 = pkin(4) * t485 + pkin(8) * t502;
t408 = pkin(4) * t533 + t441 * pkin(8);
t505 = qJD(3) * t462;
t504 = qJD(3) * t368;
t499 = qJD(6) * t428;
t498 = t441 * MDP(20);
t319 = -t373 * t436 + t438 * t378;
t304 = -pkin(4) * t506 - pkin(9) * t403 + t319;
t287 = t304 * t443 - t536;
t497 = qJD(6) - t287;
t302 = -t341 * t436 + t438 * t358;
t296 = qJD(2) * t454 + t302;
t299 = -pkin(9) * t471 + t303;
t492 = t440 * t296 + t443 * t299 + t304 * t500;
t431 = -pkin(4) * t438 - pkin(3);
t487 = t442 * t508;
t482 = t441 * t496;
t480 = -qJD(3) * pkin(3) + qJD(4);
t479 = t495 * t440;
t477 = t443 * t296 - t440 * t299 - t304 * t501 - t311 * t500;
t327 = pkin(4) * t471 + t344;
t476 = t441 * t490;
t475 = t348 * t490;
t351 = t402 * t440 + t403 * t443;
t474 = t351 * t490;
t412 = -t490 - t537;
t468 = -t412 - t490;
t288 = t304 * t440 + t311 * t443;
t464 = t355 * t443 - t366 * t440;
t394 = t439 * t441 + t444 * t531;
t362 = -t394 * t436 - t438 * t530;
t363 = t394 * t438 - t436 * t530;
t463 = t362 * t443 - t363 * t440;
t313 = t362 * t440 + t363 * t443;
t460 = -t287 * t428 - t492;
t393 = -t439 * t444 + t441 * t531;
t456 = -t311 * t501 + t492;
t452 = qJD(3) * t457;
t369 = -t540 + t480;
t279 = -pkin(5) * t482 - t477;
t451 = qJD(3) * (-t468 - t537);
t450 = -t440 * t494 + t443 * t495;
t449 = -qJ(4) * t503 + (-t369 + t480) * t444;
t337 = -pkin(4) * t402 + t369;
t317 = qJD(2) * t452 + qJD(5) * t351;
t280 = pkin(5) * t317 + qJ(6) * t316 - qJD(6) * t351 + t327;
t447 = qJD(2) ^ 2;
t386 = t405 * t441;
t385 = t406 * t441;
t379 = -pkin(8) * t532 + t401;
t365 = qJD(3) * t394 + t441 * t486;
t364 = -qJD(3) * t393 + t444 * t486;
t346 = pkin(5) * t405 - qJ(6) * t406 + t431;
t339 = t539 * t441 + t452;
t338 = t391 * t441 + t440 * t485 - t473;
t334 = t364 * t438 + t436 * t487;
t333 = -t364 * t436 + t438 * t487;
t324 = pkin(5) * t385 + qJ(6) * t386 + t408;
t310 = pkin(5) * t351 + qJ(6) * t348;
t309 = pkin(5) * t444 - t464;
t308 = -qJ(6) * t444 + t517;
t293 = -t316 - t547;
t292 = pkin(5) * t348 - qJ(6) * t351 + t337;
t289 = pkin(5) * t339 + qJ(6) * t338 + qJD(6) * t386 + t397;
t286 = qJD(5) * t313 - t333 * t443 + t334 * t440;
t285 = qJD(5) * t463 + t333 * t440 + t334 * t443;
t284 = -qJ(6) * t428 + t288;
t282 = pkin(5) * t428 + t497;
t278 = qJ(6) * t482 + t456 - t499;
t1 = [-t364 * qJD(3) * MDP(11) + (-t333 * t403 + t334 * t402) * MDP(14) + (t302 * t362 + t303 * t363 + t319 * t333 + t320 * t334 + t344 * t393) * MDP(15) + (-t285 * t348 + t286 * t351 - t313 * t317 + t316 * t463) * MDP(24) + (t278 * t313 - t279 * t463 + t280 * t393 + t282 * t286 + t284 * t285) * MDP(26) + (-MDP(4) * t445 + (-MDP(10) * t444 + MDP(11) * t441 - MDP(3)) * t442) * t447 * t437 + (-qJD(3) * MDP(10) - t402 * MDP(12) + t403 * MDP(13) + t369 * MDP(15) + t292 * MDP(26) + t351 * t494) * t365 + ((-MDP(12) * t333 + MDP(13) * t334) * t444 + ((-MDP(11) * t530 + (-t362 * t438 - t363 * t436) * MDP(14) + (t436 * MDP(12) + t438 * MDP(13)) * t393) * t444 + (-MDP(10) * t530 + t362 * MDP(12) - t363 * MDP(13) - t313 * t494 + t463 * t495) * t441) * qJD(3)) * qJD(2) + t495 * (t286 * t428 + t393 * t317 + t365 * t348) + t494 * (t285 * t428 - t316 * t393); 0.2e1 * t481 * t546 - 0.2e1 * t496 * t545 + MDP(7) * t524 - MDP(8) * t527 + (-pkin(8) * t524 + t441 * t451) * MDP(10) + (pkin(8) * t527 + t444 * t451) * MDP(11) + ((t402 * t490 + t535 + (qJD(2) * t379 + t319) * qJD(3)) * t441 + (-t302 + (-pkin(8) * t402 + t369 * t436) * qJD(3) + (t421 - t516) * qJD(2)) * t444) * MDP(12) + ((-t403 * t490 + t534 + (-qJD(2) * t380 - t320) * qJD(3)) * t441 + (t303 + (pkin(8) * t403 + t369 * t438) * qJD(3) + (t478 + t515) * qJD(2)) * t444) * MDP(13) + ((-t302 * t438 - t303 * t436) * t441 - t516 * t403 + t515 * t402 + (-t319 * t438 - t320 * t436 + (-t379 * t438 - t380 * t436) * qJD(2)) * t502) * MDP(14) + (-t369 * t476 + t302 * t379 + t303 * t380 + t515 * t320 + t516 * t319 + (t344 * t441 + t369 * t502) * pkin(8)) * MDP(15) + (t316 * t386 - t338 * t351) * MDP(16) + (t316 * t385 + t317 * t386 + t338 * t348 - t339 * t351) * MDP(17) + (t316 * t444 + t338 * t428 + (-qJD(2) * t386 + t351) * t503) * MDP(18) + (t317 * t444 + t339 * t428 + (-qJD(2) * t385 - t348) * t503) * MDP(19) + (-t428 - t506) * qJD(3) * t498 + (-t477 * t444 + t397 * t348 + t408 * t317 + t327 * t385 + t337 * t339 + t544 * t428 + (-t475 + (qJD(2) * t464 + t287) * qJD(3)) * t441) * MDP(21) + (t456 * t444 + t397 * t351 - t408 * t316 - t327 * t386 - t337 * t338 + t543 * t428 + (-t474 + (-qJD(2) * t517 - t288) * qJD(3)) * t441) * MDP(22) + (t279 * t444 + t280 * t385 + t289 * t348 + t292 * t339 + t317 * t324 + t522 * t428 + (-t475 + (-qJD(2) * t309 - t282) * qJD(3)) * t441) * MDP(23) + (-t278 * t385 - t279 * t386 - t282 * t338 - t284 * t339 - t308 * t317 - t309 * t316 - t348 * t523 + t351 * t522) * MDP(24) + (-t278 * t444 + t280 * t386 - t289 * t351 + t292 * t338 + t316 * t324 - t523 * t428 + (t474 + (qJD(2) * t308 + t284) * qJD(3)) * t441) * MDP(25) + (t278 * t308 + t279 * t309 + t280 * t324 + (t289 - t476) * t292 + t523 * t284 + t522 * t282) * MDP(26); (qJD(3) * t377 - t412 * t507 - t344) * MDP(10) + t468 * t506 * MDP(11) + (-t534 + t377 * t402 + (-t319 * t441 + t335 * t444 + t436 * t449) * qJD(2)) * MDP(12) + (t535 - t377 * t403 + (t320 * t441 - t336 * t444 + t438 * t449) * qJD(2)) * MDP(13) + (t335 * t403 - t336 * t402 + (qJD(4) * t402 + t319 * t506 + t303) * t438 + (qJD(4) * t403 + t320 * t506 - t302) * t436) * MDP(14) + (-pkin(3) * t344 - t319 * t335 - t320 * t336 - t369 * t377 + (-t319 * t436 + t320 * t438) * qJD(4) + (-t302 * t436 + t303 * t438) * qJ(4)) * MDP(15) + (-t316 * t406 + t351 * t513) * MDP(16) + (t316 * t405 - t317 * t406 - t348 * t513 - t351 * t512) * MDP(17) + (-t513 * t428 + (qJD(3) * t406 - t351) * t507) * MDP(18) + (t512 * t428 + (-qJD(3) * t405 + t348) * t507) * MDP(19) + t428 * qJD(2) * t498 + (t431 * t317 + t327 * t405 - t357 * t348 + t541 * t428 + t512 * t337 + (-t287 + t505) * t507) * MDP(21) + (-t431 * t316 + t327 * t406 - t357 * t351 + t542 * t428 + t513 * t337 + (t288 - t504) * t507) * MDP(22) + (t280 * t405 + t317 * t346 + t520 * t428 - t519 * t348 + t512 * t292 + (t282 + t505) * t507) * MDP(23) + (-t278 * t405 + t279 * t406 + t282 * t513 - t284 * t512 + t316 * t462 - t317 * t368 + t348 * t521 + t351 * t520) * MDP(24) + (-t280 * t406 + t316 * t346 + t521 * t428 + t519 * t351 - t513 * t292 + (-t284 + t504) * t507) * MDP(25) + (t278 * t368 - t279 * t462 + t280 * t346 + t282 * t520 - t284 * t521 - t292 * t519) * MDP(26) + (-t444 * t546 + t545) * t447; (-t402 ^ 2 - t403 ^ 2) * MDP(14) + (t319 * t403 - t320 * t402 + t344) * MDP(15) - t348 ^ 2 * MDP(24) + (t284 * t348 + t280) * MDP(26) + (-t351 * MDP(24) - t282 * MDP(26) - t428 * t495) * t351 + (t402 * t479 + t403 * t450) * qJD(5) + (-t403 * MDP(12) - t402 * MDP(13) + ((MDP(13) + t479) * t438 + (MDP(12) + t450) * t436) * qJD(3)) * t506 - t494 * (-t514 - t547); t293 * MDP(18) + t460 * MDP(22) + (pkin(5) * t316 - qJ(6) * t317) * MDP(24) + (-t460 - 0.2e1 * t499) * MDP(25) + (-pkin(5) * t279 + qJ(6) * t278 - t282 * t288 + t284 * t497 - t292 * t310) * MDP(26) + (-MDP(19) * t351 + t494 * t536) * qJD(5) + (-t428 * MDP(19) - t337 * MDP(21) - t292 * MDP(23) + (t284 - t288) * MDP(24) + t310 * MDP(25) + MDP(17) * t351) * t351 + (-MDP(19) * t457 + (0.2e1 * pkin(5) * MDP(23) + 0.2e1 * qJ(6) * MDP(25) + MDP(20)) * t441) * t496 + (t351 * MDP(16) + t337 * MDP(22) - t310 * MDP(23) + (t282 - t497) * MDP(24) - t292 * MDP(25) - MDP(17) * t348) * t348 + t495 * (-t288 * t428 + t477); (t348 * t351 - t482) * MDP(23) + t293 * MDP(24) + (-t351 ^ 2 - t428 ^ 2) * MDP(25) + (t284 * t428 + t292 * t351 + t279) * MDP(26);];
tauc  = t1;
