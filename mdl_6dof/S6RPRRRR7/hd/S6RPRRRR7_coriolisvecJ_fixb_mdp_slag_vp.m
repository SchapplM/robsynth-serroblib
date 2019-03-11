% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:18:04
% EndTime: 2019-03-09 07:18:14
% DurationCPUTime: 5.05s
% Computational Cost: add. (5836->378), mult. (12779->509), div. (0->0), fcn. (9342->8), ass. (0->177)
t432 = cos(qJ(6));
t491 = qJD(6) * t432;
t430 = sin(qJ(4));
t434 = cos(qJ(4));
t435 = cos(qJ(3));
t498 = qJD(1) * t435;
t431 = sin(qJ(3));
t499 = qJD(1) * t431;
t391 = t430 * t499 - t434 * t498;
t392 = -t430 * t498 - t434 * t499;
t429 = sin(qJ(5));
t433 = cos(qJ(5));
t354 = t391 * t429 + t433 * t392;
t545 = t354 * t432;
t550 = t491 - t545;
t396 = t430 * t435 + t431 * t434;
t423 = qJD(3) + qJD(4);
t365 = t423 * t396;
t495 = qJD(4) * t430;
t497 = qJD(3) * t431;
t512 = t434 * t435;
t366 = t423 * t512 - t430 * t497 - t431 * t495;
t452 = t430 * t431 - t512;
t454 = t433 * t396 - t429 * t452;
t319 = qJD(5) * t454 + t433 * t365 + t366 * t429;
t421 = qJD(5) + t423;
t549 = t319 * t421;
t535 = qJD(1) * t452;
t440 = t423 * t535;
t436 = -pkin(1) - pkin(7);
t405 = qJD(1) * t436 + qJD(2);
t481 = pkin(8) * qJD(1) - t405;
t376 = t481 * t497;
t382 = -pkin(8) * t499 + t405 * t431;
t473 = t430 * t376 - t382 * t495;
t383 = -pkin(8) * t498 + t435 * t405;
t375 = qJD(3) * pkin(3) + t383;
t496 = qJD(3) * t435;
t377 = t481 * t496;
t536 = t434 * (qJD(4) * t375 - t377);
t304 = pkin(9) * t440 + t473 + t536;
t386 = t391 * pkin(9);
t372 = t430 * t382;
t475 = t434 * t375 - t372;
t330 = t386 + t475;
t328 = pkin(4) * t423 + t330;
t358 = t365 * qJD(1);
t373 = t434 * t382;
t456 = -t375 * t430 - t373;
t474 = t434 * t376 + t430 * t377;
t443 = qJD(4) * t456 + t474;
t305 = pkin(9) * t358 + t443;
t532 = pkin(9) * t392;
t331 = -t456 + t532;
t494 = qJD(5) * t429;
t478 = t305 * t429 - t331 * t494;
t277 = (qJD(5) * t328 + t304) * t433 + t478;
t402 = pkin(3) * t499 + qJD(1) * qJ(2);
t367 = -pkin(4) * t392 + t402;
t521 = t354 * t367;
t548 = -t277 - t521;
t524 = t331 * t433;
t302 = t328 * t429 + t524;
t479 = t304 * t429 - t433 * t305;
t278 = qJD(5) * t302 + t479;
t455 = t433 * t391 - t392 * t429;
t523 = t455 * t367;
t547 = t523 - t278;
t428 = sin(qJ(6));
t492 = qJD(6) * t428;
t493 = qJD(5) * t433;
t313 = -t433 * t358 + t391 * t494 + t392 * t493 + t429 * t440;
t505 = t432 * t313 + t421 * t491;
t292 = t455 * t492 + t505;
t290 = t292 * t428;
t291 = t292 * t432;
t340 = t421 * t428 - t432 * t455;
t527 = t313 * t428;
t293 = t340 * qJD(6) + t527;
t444 = qJD(5) * t455 + t358 * t429 + t433 * t440;
t312 = t432 * t444;
t519 = t455 * t428;
t338 = -t432 * t421 - t519;
t489 = -qJD(6) + t354;
t310 = t428 * t444;
t506 = -t489 * t491 - t310;
t544 = t489 * t428;
t546 = -t354 ^ 2 * MDP(22) + (-t354 * t421 + t313) * MDP(23) + t444 * MDP(24) + (t354 * MDP(21) + MDP(22) * t455 - t421 * MDP(24) - MDP(32) * t489) * t455 + (t550 * t340 + t290) * MDP(28) + (t340 * t455 + t489 * t545 + t506) * MDP(30) + (-t338 * t455 - t489 * t544 - t312) * MDP(31) + (-t428 * t293 - t550 * t338 + t340 * t544 + t291) * MDP(29);
t525 = t331 * t429;
t301 = t328 * t433 - t525;
t299 = -pkin(5) * t421 - t301;
t529 = t299 * t354;
t543 = qJ(2) * MDP(6) + MDP(5);
t300 = pkin(10) * t421 + t302;
t315 = -pkin(5) * t354 + pkin(10) * t455 + t367;
t287 = t300 * t432 + t315 * t428;
t461 = t278 * t428 - t287 * t455 + t299 * t491;
t458 = t300 * t428 - t315 * t432;
t451 = -t278 * t432 + t299 * t492 - t458 * t455;
t327 = -pkin(5) * t455 - pkin(10) * t354;
t539 = MDP(7) * t431;
t538 = MDP(8) * (t431 ^ 2 - t435 ^ 2);
t361 = t396 * t429 + t433 * t452;
t441 = -qJD(5) * t361 - t365 * t429 + t433 * t366;
t537 = t421 * t441;
t531 = pkin(8) - t436;
t394 = t531 * t497;
t401 = t531 * t435;
t395 = qJD(3) * t401;
t400 = t531 * t431;
t516 = t401 * t434;
t448 = -qJD(4) * t516 + t430 * t394 - t434 * t395 + t400 * t495;
t322 = -pkin(9) * t366 + t448;
t453 = t400 * t434 + t401 * t430;
t442 = qJD(4) * t453 + t434 * t394 + t395 * t430;
t323 = pkin(9) * t365 + t442;
t344 = pkin(9) * t452 + t400 * t430 - t516;
t345 = -pkin(9) * t396 - t453;
t457 = t344 * t433 - t345 * t429;
t282 = qJD(5) * t457 + t322 * t433 + t323 * t429;
t325 = t344 * t429 + t345 * t433;
t413 = t431 * pkin(3) + qJ(2);
t378 = pkin(4) * t396 + t413;
t326 = pkin(5) * t454 + pkin(10) * t361 + t378;
t534 = -t278 * t361 - t299 * t319 + t325 * t444 + (qJD(6) * t326 + t282) * t489 - t454 * (qJD(6) * t315 + t277);
t533 = pkin(4) * t391;
t528 = t299 * t361;
t526 = t326 * t444;
t518 = t365 * t423;
t517 = t366 * t423;
t515 = t402 * t391;
t514 = t429 * t430;
t513 = t430 * t433;
t438 = qJD(1) ^ 2;
t511 = t435 * t438;
t437 = qJD(3) ^ 2;
t510 = t436 * t437;
t472 = -t383 * t430 - t373;
t332 = t472 - t532;
t503 = t434 * t383 - t372;
t333 = t386 + t503;
t416 = pkin(3) * t434 + pkin(4);
t507 = t332 * t429 + t333 * t433 - t416 * t493 - (-t430 * t494 + (t433 * t434 - t514) * qJD(4)) * pkin(3);
t504 = t332 * t433 - t333 * t429 + t416 * t494 + (t430 * t493 + (t429 * t434 + t513) * qJD(4)) * pkin(3);
t424 = qJD(1) * qJD(2);
t488 = qJD(1) * qJD(3);
t484 = t435 * t488;
t399 = pkin(3) * t484 + t424;
t406 = pkin(3) * t496 + qJD(2);
t487 = 0.2e1 * qJD(1);
t419 = pkin(3) * t498;
t483 = -pkin(3) * t423 - t375;
t482 = -pkin(4) * t421 - t328;
t321 = t327 - t533;
t385 = pkin(3) * t513 + t416 * t429 + pkin(10);
t465 = qJD(6) * t385 + t321 + t419;
t414 = pkin(4) * t429 + pkin(10);
t464 = qJD(6) * t414 + t321;
t306 = t330 * t429 + t524;
t463 = pkin(4) * t494 - t306;
t307 = t330 * t433 - t525;
t462 = -pkin(4) * t493 + t307;
t347 = pkin(4) * t366 + t406;
t460 = t385 * t444 - t529;
t459 = t414 * t444 - t529;
t450 = -t402 * t392 - t473;
t449 = -t319 * t432 + t361 * t492;
t342 = -pkin(4) * t440 + t399;
t439 = t391 * t392 * MDP(14) + (-t392 * t423 - t358) * MDP(16) + (-t391 * t423 + t440) * MDP(17) + (t391 ^ 2 - t392 ^ 2) * MDP(15) + t546;
t415 = -pkin(4) * t433 - pkin(5);
t384 = pkin(3) * t514 - t416 * t433 - pkin(5);
t370 = t419 - t533;
t288 = pkin(5) * t441 + pkin(10) * t319 + t347;
t285 = -pkin(5) * t444 - t313 * pkin(10) + t342;
t284 = t432 * t285;
t283 = qJD(5) * t325 + t322 * t429 - t323 * t433;
t1 = [(-t441 * t489 - t444 * t454) * MDP(32) + (-t283 * t421 + t342 * t454 - t347 * t354 + t367 * t441 - t378 * t444) * MDP(26) + (-t313 * t454 - t319 * t354 - t361 * t444 + t441 * t455) * MDP(22) + (-(-t338 * t432 - t340 * t428) * t319 - (-t290 - t293 * t432 + (t338 * t428 - t340 * t432) * qJD(6)) * t361) * MDP(29) + (-t313 * t361 + t319 * t455) * MDP(21) + (-t282 * t421 + t313 * t378 - t319 * t367 - t342 * t361 - t347 * t455) * MDP(27) + (-t361 * t310 - t293 * t454 - t441 * t338 - (t319 * t428 + t361 * t491) * t489) * MDP(31) + (t292 * t454 + t312 * t361 + t340 * t441 - t449 * t489) * MDP(30) + (t358 * t396 - t365 * t392 + t391 * t366 - t440 * t452) * MDP(15) + (t358 * t452 + t365 * t391) * MDP(14) + (-t413 * t358 - t402 * t365 - t406 * t391 - t399 * t452 - t423 * t448) * MDP(20) + (-t435 * MDP(10) - t431 * MDP(9)) * t437 + (t283 * t338 + t284 * t454 - t458 * t441 - t457 * t293 + (-t288 * t489 - t526 + (-t300 * t454 + t325 * t489 - t528) * qJD(6)) * t432 + t534 * t428) * MDP(33) + (t283 * t340 - t287 * t441 - t457 * t292 + ((-qJD(6) * t325 + t288) * t489 + t526 - (-qJD(6) * t300 + t285) * t454 + qJD(6) * t528) * t428 + t534 * t432) * MDP(34) + (-t291 * t361 + t340 * t449) * MDP(28) + (t402 * t366 - t406 * t392 + t399 * t396 + (-t413 * t535 + t442) * t423) * MDP(19) - MDP(23) * t549 + 0.2e1 * t543 * t424 - MDP(17) * t517 + (-t435 * t510 + (-qJ(2) * t497 + qJD(2) * t435) * t487) * MDP(13) + (-t431 * t510 + (qJ(2) * t496 + qJD(2) * t431) * t487) * MDP(12) + 0.2e1 * t488 * t538 - 0.2e1 * t484 * t539 - MDP(24) * t537 - MDP(16) * t518; (qJD(1) * t392 - t518) * MDP(19) + (qJD(1) * t391 - t517) * MDP(20) + (qJD(1) * t354 - t549) * MDP(26) + (qJD(1) * t455 - t537) * MDP(27) + (t293 * t361 + t310 * t454 + t319 * t338) * MDP(33) + (t292 * t361 + t312 * t454 + t319 * t340) * MDP(34) - t543 * t438 - ((-qJD(1) * t432 - t428 * t441 - t454 * t491) * MDP(33) + (qJD(1) * t428 - t432 * t441 + t454 * t492) * MDP(34)) * t489 + (MDP(12) * t431 + MDP(13) * t435) * (-t437 - t438); (t391 * t419 + t503 * t423 + (qJD(4) * t483 + t377) * t434 + t450) * MDP(20) + (t392 * t419 + t515 - t472 * t423 + (t430 * t483 - t373) * qJD(4) + t474) * MDP(19) + (t370 * t455 + t421 * t507 + t548) * MDP(27) + (t354 * t370 - t421 * t504 + t547) * MDP(26) + t439 + (t384 * t293 + t460 * t428 + t504 * t338 - (t428 * t507 - t432 * t465) * t489 + t451) * MDP(33) + (t384 * t292 + t460 * t432 + t504 * t340 - (t428 * t465 + t432 * t507) * t489 + t461) * MDP(34) - t438 * t538 + t511 * t539 + (MDP(13) * t431 * t438 - MDP(12) * t511) * qJ(2); (-t455 * t533 + t307 * t421 - t521 + (qJD(5) * t482 - t304) * t433 - t478) * MDP(27) + (-t354 * t533 + t306 * t421 + t523 + (t429 * t482 - t524) * qJD(5) - t479) * MDP(26) + t439 + (t415 * t293 + t459 * t428 + t463 * t338 - (t428 * t462 - t432 * t464) * t489 + t451) * MDP(33) + (t415 * t292 + t459 * t432 + t463 * t340 - (t428 * t464 + t432 * t462) * t489 + t461) * MDP(34) + (-t423 * t456 + t443 + t515) * MDP(19) + (t423 * t475 + t450 - t536) * MDP(20); (t302 * t421 + t547) * MDP(26) + (t301 * t421 + t548) * MDP(27) + (-pkin(5) * t293 + (-t301 * t428 + t327 * t432) * t489 - t302 * t338 - t428 * t529 - t506 * pkin(10) + t451) * MDP(33) + (-pkin(5) * t292 - (t301 * t432 + t327 * t428) * t489 - t302 * t340 - t299 * t545 + (-t489 * t492 + t312) * pkin(10) + t461) * MDP(34) + t546; t340 * t338 * MDP(28) + (-t338 ^ 2 + t340 ^ 2) * MDP(29) + (-t338 * t489 + t505) * MDP(30) + (-t340 * t489 - t527) * MDP(31) - t444 * MDP(32) + (-t277 * t428 - t287 * t489 - t299 * t340 + t284) * MDP(33) + (-t277 * t432 - t285 * t428 + t299 * t338 + t458 * t489) * MDP(34) + (MDP(30) * t519 - MDP(31) * t340 - MDP(33) * t287 + MDP(34) * t458) * qJD(6);];
tauc  = t1;
