% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:18
% EndTime: 2019-03-08 21:54:26
% DurationCPUTime: 4.73s
% Computational Cost: add. (4360->363), mult. (11386->517), div. (0->0), fcn. (9244->12), ass. (0->179)
t445 = cos(qJ(6));
t494 = qJD(6) * t445;
t437 = sin(pkin(12));
t439 = cos(pkin(12));
t443 = sin(qJ(3));
t447 = cos(qJ(3));
t416 = -t437 * t443 + t439 * t447;
t409 = t416 * qJD(2);
t446 = cos(qJ(5));
t397 = t446 * t409;
t417 = t437 * t447 + t439 * t443;
t411 = t417 * qJD(2);
t442 = sin(qJ(5));
t364 = -t411 * t442 + t397;
t553 = t364 * t445;
t557 = t494 - t553;
t532 = -qJ(4) - pkin(8);
t480 = qJD(3) * t532;
t405 = qJD(4) * t447 + t443 * t480;
t406 = -qJD(4) * t443 + t447 * t480;
t448 = cos(qJ(2));
t438 = sin(pkin(6));
t501 = qJD(1) * t438;
t486 = t448 * t501;
t508 = -t405 * t437 + t439 * t406 + t417 * t486;
t507 = t439 * t405 + t437 * t406 - t416 * t486;
t462 = t409 * t442 + t446 * t411;
t441 = sin(qJ(6));
t495 = qJD(6) * t441;
t410 = t417 * qJD(3);
t400 = qJD(2) * t410;
t490 = qJD(2) * qJD(3);
t481 = t447 * t490;
t482 = t443 * t490;
t401 = -t437 * t482 + t439 * t481;
t496 = qJD(5) * t442;
t324 = qJD(5) * t397 - t442 * t400 + t446 * t401 - t411 * t496;
t434 = qJD(3) + qJD(5);
t509 = t445 * t324 + t434 * t494;
t305 = -t462 * t495 + t509;
t304 = t305 * t445;
t352 = t434 * t441 + t445 * t462;
t528 = t324 * t441;
t306 = qJD(6) * t352 + t528;
t520 = t462 * t441;
t350 = -t445 * t434 + t520;
t556 = -t441 * t306 - t557 * t350 + t304;
t303 = t305 * t441;
t325 = qJD(5) * t462 + t446 * t400 + t401 * t442;
t491 = -qJD(6) + t364;
t318 = t441 * t325;
t510 = -t491 * t494 + t318;
t522 = t364 * t434;
t524 = t462 * t434;
t526 = t352 * t462;
t555 = (-t325 + t524) * MDP(17) - t364 ^ 2 * MDP(15) + (-MDP(14) * t364 + MDP(15) * t462 + MDP(25) * t491) * t462 + (t324 - t522) * MDP(16) + (t557 * t352 + t303) * MDP(21) + (t491 * t553 + t510 - t526) * MDP(23);
t554 = t364 * t441;
t413 = t416 * qJD(3);
t552 = pkin(9) * t413 - t508;
t551 = -pkin(9) * t410 + t507;
t444 = sin(qJ(2));
t487 = t444 * t501;
t497 = qJD(3) * t443;
t550 = pkin(3) * t497 - t487;
t330 = pkin(5) * t462 - pkin(10) * t364;
t421 = qJD(2) * pkin(8) + t487;
t440 = cos(pkin(6));
t500 = qJD(1) * t440;
t428 = t447 * t500;
t467 = qJD(4) + t486;
t348 = (-t421 * t443 + t428) * qJD(3) + (-qJ(4) * t497 + t447 * t467) * qJD(2);
t485 = t443 * t500;
t349 = (-t421 * t447 - t485) * qJD(3) + (-qJ(4) * qJD(3) * t447 - t443 * t467) * qJD(2);
t311 = -t348 * t437 + t439 * t349;
t308 = -pkin(9) * t401 + t311;
t312 = t439 * t348 + t437 * t349;
t309 = -pkin(9) * t400 + t312;
t474 = qJ(4) * qJD(2) + t421;
t380 = t447 * t474 + t485;
t372 = t437 * t380;
t379 = -t443 * t474 + t428;
t376 = qJD(3) * pkin(3) + t379;
t333 = t439 * t376 - t372;
t533 = pkin(9) * t411;
t321 = qJD(3) * pkin(4) + t333 - t533;
t518 = t439 * t380;
t334 = t437 * t376 + t518;
t534 = pkin(9) * t409;
t326 = t334 + t534;
t283 = (qJD(5) * t321 + t309) * t446 + t308 * t442 - t326 * t496;
t488 = -pkin(3) * t447 - pkin(2);
t402 = qJD(2) * t488 + qJD(4) - t486;
t367 = -pkin(4) * t409 + t402;
t548 = -t364 * t367 - t283;
t544 = MDP(5) * t443;
t543 = MDP(6) * (t443 ^ 2 - t447 ^ 2);
t527 = t350 * t462;
t469 = pkin(4) * t410 + t550;
t320 = t445 * t325;
t541 = -t491 * t495 - t320;
t540 = MDP(10) * t443 + MDP(11) * t447;
t298 = t321 * t442 + t326 * t446;
t284 = qJD(5) * t298 - t446 * t308 + t309 * t442;
t297 = t321 * t446 - t326 * t442;
t295 = -pkin(5) * t434 - t297;
t296 = pkin(10) * t434 + t298;
t310 = -pkin(5) * t364 - pkin(10) * t462 + t367;
t465 = t296 * t441 - t310 * t445;
t539 = -t284 * t445 + t295 * t495 + t462 * t465;
t286 = t296 * t445 + t310 * t441;
t538 = t284 * t441 + t286 * t462 + t295 * t494;
t537 = -t367 * t462 - t284;
t371 = t416 * t442 + t417 * t446;
t391 = -pkin(4) * t416 + t488;
t461 = t446 * t416 - t417 * t442;
t314 = -pkin(5) * t461 - pkin(10) * t371 + t391;
t423 = t532 * t443;
t424 = t532 * t447;
t383 = t439 * t423 + t424 * t437;
t356 = -pkin(9) * t417 + t383;
t384 = t437 * t423 - t439 * t424;
t357 = pkin(9) * t416 + t384;
t317 = t356 * t442 + t357 * t446;
t331 = qJD(5) * t461 - t410 * t442 + t413 * t446;
t464 = t356 * t446 - t357 * t442;
t514 = -qJD(5) * t464 + t442 * t552 - t446 * t551;
t536 = (qJD(6) * t310 + t283) * t461 + t284 * t371 + t295 * t331 - (-qJD(6) * t314 + t514) * t491 - t317 * t325;
t535 = pkin(3) * t437;
t531 = qJD(2) * pkin(2);
t530 = t295 * t371;
t529 = t314 * t325;
t525 = t352 * t441;
t519 = t438 * t444;
t449 = qJD(3) ^ 2;
t517 = t443 * t449;
t516 = t447 * t449;
t513 = qJD(5) * t317 + t442 * t551 + t446 * t552;
t336 = -t379 * t437 - t518;
t327 = t336 - t534;
t338 = t439 * t379 - t372;
t328 = t338 - t533;
t431 = pkin(3) * t439 + pkin(4);
t458 = t431 * t446 - t442 * t535;
t512 = -t458 * qJD(5) + t327 * t442 + t328 * t446;
t459 = t431 * t442 + t446 * t535;
t506 = t459 * qJD(5) + t327 * t446 - t328 * t442;
t407 = pkin(3) * t482 + qJD(2) * t487;
t499 = qJD(2) * t443;
t498 = qJD(2) * t444;
t493 = qJD(6) * t448;
t484 = qJD(2) * t438 * t448;
t385 = pkin(3) * t499 + pkin(4) * t411;
t475 = t491 * t441;
t404 = pkin(10) + t459;
t471 = qJD(6) * t404 + t330 + t385;
t366 = pkin(4) * t400 + t407;
t332 = qJD(5) * t371 + t446 * t410 + t413 * t442;
t470 = pkin(5) * t332 - pkin(10) * t331 + t469;
t466 = -t295 * t364 - t325 * t404;
t414 = t440 * t447 - t443 * t519;
t415 = t440 * t443 + t447 * t519;
t359 = t414 * t439 - t415 * t437;
t360 = t414 * t437 + t415 * t439;
t463 = t359 * t446 - t360 * t442;
t323 = t359 * t442 + t360 * t446;
t460 = -t491 * t554 - t541;
t457 = t331 * t445 - t371 * t495;
t455 = -0.2e1 * qJD(3) * t531;
t450 = qJD(2) ^ 2;
t403 = -pkin(5) - t458;
t378 = -qJD(3) * t415 - t443 * t484;
t377 = qJD(3) * t414 + t447 * t484;
t337 = t377 * t439 + t378 * t437;
t335 = -t377 * t437 + t378 * t439;
t292 = pkin(5) * t325 - pkin(10) * t324 + t366;
t291 = t445 * t292;
t288 = qJD(5) * t323 - t335 * t446 + t337 * t442;
t287 = qJD(5) * t463 + t335 * t442 + t337 * t446;
t1 = [(-t335 * t411 + t337 * t409 - t359 * t401 - t360 * t400) * MDP(12) + (t311 * t359 + t312 * t360 + t333 * t335 + t334 * t337) * MDP(13) + (-(-t287 * t441 - t323 * t494) * t491 - t323 * t318 + t288 * t350 - t463 * t306) * MDP(26) + ((t287 * t445 - t323 * t495) * t491 - t323 * t320 + t288 * t352 - t463 * t305) * MDP(27) + (-MDP(19) * t288 - MDP(20) * t287) * t434 + (MDP(10) * t378 - MDP(11) * t377) * qJD(3) + ((t402 * t498 - t407 * t448) * MDP(13) + (-t325 * t448 - t364 * t498) * MDP(19) + (-t324 * t448 + t462 * t498) * MDP(20) + (-(t441 * t493 + t445 * t498) * t491 - t448 * t320) * MDP(26) + ((t441 * t498 - t445 * t493) * t491 + t448 * t318) * MDP(27) - t540 * t448 * t490 + (-t448 * MDP(4) + (-MDP(10) * t447 + MDP(11) * t443 - MDP(3)) * t444) * t450) * t438; 0.2e1 * t481 * t544 - 0.2e1 * t490 * t543 + MDP(7) * t516 - MDP(8) * t517 + (-pkin(8) * t516 + t443 * t455) * MDP(10) + (pkin(8) * t517 + t447 * t455) * MDP(11) + (-t311 * t417 + t312 * t416 - t333 * t413 - t334 * t410 - t383 * t401 - t384 * t400 + t409 * t507 - t411 * t508) * MDP(12) + (t311 * t383 + t312 * t384 + t508 * t333 + t507 * t334 + t402 * t550 + t407 * t488) * MDP(13) + (t324 * t371 + t331 * t462) * MDP(14) + (t324 * t461 - t325 * t371 + t331 * t364 - t332 * t462) * MDP(15) + (t325 * t391 + t332 * t367 - t364 * t469 - t366 * t461) * MDP(19) + (t324 * t391 + t331 * t367 + t366 * t371 + t462 * t469) * MDP(20) + (t304 * t371 + t352 * t457) * MDP(21) + ((-t350 * t445 - t525) * t331 + (-t303 - t306 * t445 + (t350 * t441 - t352 * t445) * qJD(6)) * t371) * MDP(22) + (-t305 * t461 + t320 * t371 + t332 * t352 - t457 * t491) * MDP(23) + (-t371 * t318 + t306 * t461 - t332 * t350 - (-t331 * t441 - t371 * t494) * t491) * MDP(24) + (-t325 * t461 - t332 * t491) * MDP(25) + (-t465 * t332 - t291 * t461 - t464 * t306 + t513 * t350 + (t529 - t470 * t491 + (t296 * t461 + t317 * t491 + t530) * qJD(6)) * t445 + t536 * t441) * MDP(26) + (-t286 * t332 - t464 * t305 + t513 * t352 + (-t529 + (-qJD(6) * t296 + t292) * t461 - qJD(6) * t530 - (qJD(6) * t317 - t470) * t491) * t441 + t536 * t445) * MDP(27) + (MDP(16) * t331 - MDP(17) * t332 - MDP(19) * t513 + MDP(20) * t514) * t434; ((t334 + t336) * t411 + (t333 - t338) * t409 + (-t400 * t437 - t401 * t439) * pkin(3)) * MDP(12) + (-t333 * t336 - t334 * t338 + (t311 * t439 + t312 * t437 - t402 * t499) * pkin(3)) * MDP(13) + (t364 * t385 - t434 * t506 + t537) * MDP(19) + (-t385 * t462 + t434 * t512 + t548) * MDP(20) + (t491 * t525 + t556) * MDP(22) + (t460 + t527) * MDP(24) + (t403 * t306 + t466 * t441 + t506 * t350 - (t441 * t512 - t445 * t471) * t491 + t539) * MDP(26) + (t403 * t305 + t466 * t445 + t506 * t352 - (t441 * t471 + t445 * t512) * t491 + t538) * MDP(27) + t540 * qJD(2) * t531 + (-t447 * t544 + t543) * t450 + t555; (-t409 ^ 2 - t411 ^ 2) * MDP(12) + (t333 * t411 - t334 * t409 + t407) * MDP(13) + (t325 + t524) * MDP(19) + (t324 + t522) * MDP(20) + (t460 - t527) * MDP(26) + (-t445 * t491 ^ 2 - t318 - t526) * MDP(27); (t298 * t434 + t537) * MDP(19) + (t297 * t434 + t548) * MDP(20) + (t352 * t475 + t556) * MDP(22) + (-t475 * t491 + t320 + t527) * MDP(24) + (-pkin(5) * t306 + (-t297 * t441 + t330 * t445) * t491 - t298 * t350 - t295 * t554 - t510 * pkin(10) + t539) * MDP(26) + (-pkin(5) * t305 - (t297 * t445 + t330 * t441) * t491 - t298 * t352 - t295 * t553 + t541 * pkin(10) + t538) * MDP(27) + t555; t352 * t350 * MDP(21) + (-t350 ^ 2 + t352 ^ 2) * MDP(22) + (-t350 * t491 + t509) * MDP(23) + (-t352 * t491 - t528) * MDP(24) + t325 * MDP(25) + (-t283 * t441 - t286 * t491 - t295 * t352 + t291) * MDP(26) + (-t283 * t445 - t292 * t441 + t295 * t350 + t465 * t491) * MDP(27) + (-MDP(23) * t520 - MDP(24) * t352 - MDP(26) * t286 + MDP(27) * t465) * qJD(6);];
tauc  = t1;
