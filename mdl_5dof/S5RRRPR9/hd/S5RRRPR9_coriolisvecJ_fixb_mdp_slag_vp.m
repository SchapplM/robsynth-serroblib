% Calculate Coriolis joint torque vector for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:24:16
% EndTime: 2021-01-15 23:24:39
% DurationCPUTime: 7.10s
% Computational Cost: add. (4120->441), mult. (10450->612), div. (0->0), fcn. (7310->8), ass. (0->190)
t480 = sin(qJ(3));
t481 = sin(qJ(2));
t533 = qJD(1) * t481;
t518 = t480 * t533;
t483 = cos(qJ(3));
t523 = t483 * qJD(2);
t441 = t518 - t523;
t531 = qJD(2) * t480;
t443 = t483 * t533 + t531;
t477 = sin(pkin(9));
t478 = cos(pkin(9));
t386 = t478 * t441 + t443 * t477;
t482 = cos(qJ(5));
t551 = t482 * t386;
t479 = sin(qJ(5));
t496 = -t441 * t477 + t478 * t443;
t565 = t496 * t479;
t335 = t551 + t565;
t484 = cos(qJ(2));
t532 = qJD(1) * t484;
t463 = -qJD(3) + t532;
t455 = -qJD(5) + t463;
t566 = t335 * t455;
t522 = qJD(1) * qJD(2);
t512 = t484 * t522;
t527 = qJD(3) * t481;
t515 = t480 * t527;
t521 = qJD(2) * qJD(3);
t401 = -qJD(1) * t515 + (t512 + t521) * t483;
t471 = pkin(6) * t532;
t453 = qJD(2) * pkin(7) + t471;
t448 = -pkin(2) * t484 - pkin(7) * t481 - pkin(1);
t430 = t448 * qJD(1);
t557 = t480 * t430;
t394 = t453 * t483 + t557;
t501 = pkin(2) * t481 - pkin(7) * t484;
t445 = t501 * qJD(2);
t431 = qJD(1) * t445;
t513 = t481 * t522;
t503 = pkin(6) * t513;
t540 = -t483 * t431 - t480 * t503;
t488 = -qJD(3) * t394 - t540;
t314 = pkin(3) * t513 - qJ(4) * t401 - qJD(4) * t443 + t488;
t526 = qJD(3) * t483;
t514 = t481 * t526;
t529 = qJD(2) * t484;
t517 = t480 * t529;
t576 = t514 + t517;
t402 = qJD(1) * t576 + t480 * t521;
t528 = qJD(3) * t480;
t494 = t430 * t526 + t480 * t431 - t453 * t528;
t487 = -t483 * t503 + t494;
t319 = -qJ(4) * t402 - qJD(4) * t441 + t487;
t297 = t478 * t314 - t319 * t477;
t352 = t401 * t478 - t402 * t477;
t295 = pkin(4) * t513 - pkin(8) * t352 + t297;
t298 = t477 * t314 + t478 * t319;
t351 = t401 * t477 + t478 * t402;
t296 = -pkin(8) * t351 + t298;
t393 = t483 * t430 - t453 * t480;
t361 = -qJ(4) * t443 + t393;
t353 = -pkin(3) * t463 + t361;
t362 = -qJ(4) * t441 + t394;
t558 = t478 * t362;
t316 = t477 * t353 + t558;
t578 = pkin(8) * t386;
t307 = t316 - t578;
t524 = qJD(5) * t479;
t306 = t307 * t524;
t452 = -qJD(2) * pkin(2) + pkin(6) * t533;
t399 = pkin(3) * t441 + qJD(4) + t452;
t343 = pkin(4) * t386 + t399;
t584 = -t479 * t295 - t482 * t296 + t335 * t343 + t306;
t497 = t386 * t479 - t482 * t496;
t583 = MDP(26) * t513 + (-t335 ^ 2 + t497 ^ 2) * MDP(23) - t335 * MDP(22) * t497;
t567 = t455 * t497;
t507 = t482 * t295 - t479 * t296;
t581 = t343 * t497 + t507;
t550 = t483 * t484;
t493 = pkin(3) * t481 - qJ(4) * t550;
t568 = qJ(4) + pkin(7);
t509 = qJD(3) * t568;
t444 = t501 * qJD(1);
t537 = pkin(6) * t518 + t483 * t444;
t580 = qJD(1) * t493 + qJD(4) * t480 + t483 * t509 + t537;
t426 = t480 * t444;
t525 = qJD(4) * t483;
t554 = t481 * t483;
t555 = t480 * t484;
t579 = t426 + (-t554 * pkin(6) - qJ(4) * t555) * qJD(1) + t480 * t509 - t525;
t435 = t477 * t483 + t478 * t480;
t424 = t435 * qJD(3);
t542 = t435 * t532 - t424;
t434 = t477 * t480 - t478 * t483;
t577 = t463 * t434;
t575 = -0.2e1 * t522;
t574 = pkin(8) * t496;
t573 = MDP(4) * t481;
t475 = t481 ^ 2;
t572 = MDP(5) * (-t484 ^ 2 + t475);
t547 = -t579 * t477 + t580 * t478;
t546 = t580 * t477 + t579 * t478;
t502 = -t471 + (-t480 * t532 + t528) * pkin(3);
t506 = t482 * t351 + t352 * t479;
t302 = -qJD(5) * t497 + t506;
t570 = pkin(3) * t477;
t569 = pkin(6) * t480;
t564 = t401 * t480;
t563 = t441 * t463;
t562 = t443 * t463;
t561 = t452 * t480;
t560 = t452 * t483;
t559 = t463 * t483;
t357 = t477 * t362;
t556 = t480 * t481;
t485 = qJD(2) ^ 2;
t553 = t481 * t485;
t315 = t478 * t353 - t357;
t305 = -pkin(4) * t463 + t315 - t574;
t552 = t482 * t305;
t549 = t484 * t485;
t486 = qJD(1) ^ 2;
t548 = t484 * t486;
t382 = t482 * t434 + t435 * t479;
t545 = -qJD(5) * t382 + t479 * t542 + t482 * t577;
t383 = -t434 * t479 + t435 * t482;
t544 = qJD(5) * t383 + t479 * t577 - t482 * t542;
t465 = pkin(6) * t550;
t530 = qJD(2) * t481;
t538 = t483 * t445 + t530 * t569;
t330 = -t481 * t525 + t493 * qJD(2) + (-t465 + (qJ(4) * t481 - t448) * t480) * qJD(3) + t538;
t539 = t480 * t445 + t448 * t526;
t339 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t554 + (-qJD(4) * t481 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t484) * t480 + t539;
t304 = t477 * t330 + t478 * t339;
t322 = t478 * t361 - t357;
t437 = t483 * t448;
t390 = -qJ(4) * t554 + t437 + (-pkin(3) - t569) * t484;
t536 = t480 * t448 + t465;
t395 = -qJ(4) * t556 + t536;
t341 = t477 * t390 + t478 * t395;
t543 = -pkin(4) * t542 + t502;
t450 = t568 * t480;
t451 = t568 * t483;
t398 = -t477 * t450 + t478 * t451;
t446 = pkin(3) * t556 + t481 * pkin(6);
t519 = -qJD(5) * t551 - t479 * t351 + t482 * t352;
t404 = pkin(3) * t576 + pkin(6) * t529;
t469 = -pkin(3) * t483 - pkin(2);
t516 = t484 * t523;
t510 = MDP(15) * t530;
t384 = pkin(3) * t402 + pkin(6) * t512;
t508 = pkin(1) * t575;
t303 = t478 * t330 - t339 * t477;
t321 = -t361 * t477 - t558;
t340 = t478 * t390 - t395 * t477;
t397 = -t478 * t450 - t451 * t477;
t505 = t441 + t523;
t504 = -t443 + t531;
t373 = -pkin(8) * t434 + t398;
t500 = pkin(4) * t533 + pkin(8) * t577 + qJD(5) * t373 + t547;
t372 = -pkin(8) * t435 + t397;
t499 = pkin(8) * t542 + qJD(5) * t372 - t546;
t293 = t479 * t305 + t482 * t307;
t415 = t434 * t481;
t323 = -pkin(4) * t484 + pkin(8) * t415 + t340;
t414 = t435 * t481;
t325 = -pkin(8) * t414 + t341;
t498 = t323 * t479 + t325 * t482;
t364 = t482 * t414 - t415 * t479;
t365 = -t414 * t479 - t415 * t482;
t495 = qJD(1) * t475 - t463 * t484;
t467 = pkin(3) * t478 + pkin(4);
t492 = t467 * t479 + t482 * t570;
t491 = t467 * t482 - t479 * t570;
t301 = -t496 * t524 + t519;
t408 = pkin(4) * t434 + t469;
t391 = pkin(4) * t414 + t446;
t367 = t424 * t481 + t477 * t517 - t478 * t516;
t366 = t434 * t527 - t435 * t529;
t354 = pkin(3) * t443 + pkin(4) * t496;
t342 = -pkin(4) * t366 + t404;
t324 = pkin(4) * t351 + t384;
t311 = t322 - t574;
t310 = t321 + t578;
t309 = qJD(5) * t365 - t482 * t366 - t367 * t479;
t308 = -qJD(5) * t364 + t366 * t479 - t367 * t482;
t300 = pkin(8) * t366 + t304;
t299 = pkin(4) * t530 + pkin(8) * t367 + t303;
t292 = -t307 * t479 + t552;
t1 = [t572 * t575 + (-pkin(6) * t549 + t481 * t508) * MDP(9) + (pkin(6) * t553 + t484 * t508) * MDP(10) + (t401 * t554 + (-t515 + t516) * t443) * MDP(11) + ((-t441 * t483 - t443 * t480) * t529 + (-t564 - t402 * t483 + (t441 * t480 - t443 * t483) * qJD(3)) * t481) * MDP(12) + (t463 * t515 - t401 * t484 + (t443 * t481 + t483 * t495) * qJD(2)) * MDP(13) + (t463 * t514 + t402 * t484 + (-t441 * t481 - t480 * t495) * qJD(2)) * MDP(14) + (-t463 - t532) * t510 + (-(-t448 * t528 + t538) * t463 + (t452 * t526 + pkin(6) * t402 + (t437 * qJD(1) + t393) * qJD(2)) * t481 + ((pkin(6) * t441 + t561) * qJD(2) + (t557 + (pkin(6) * t463 + t453) * t483) * qJD(3) + t540) * t484) * MDP(16) + ((-pkin(6) * t484 * t528 + t539) * t463 + t494 * t484 + (pkin(6) * t401 - t452 * t528) * t481 + ((pkin(6) * t443 + t560) * t484 + (-pkin(6) * t559 - qJD(1) * t536 - t394) * t481) * qJD(2)) * MDP(17) + (-t297 * t484 - t303 * t463 + t351 * t446 - t366 * t399 + t384 * t414 + t386 * t404 + (qJD(1) * t340 + t315) * t530) * MDP(18) + (t298 * t484 + t304 * t463 + t352 * t446 - t367 * t399 - t384 * t415 + t496 * t404 + (-qJD(1) * t341 - t316) * t530) * MDP(19) + (t297 * t415 - t298 * t414 - t303 * t496 - t304 * t386 + t315 * t367 + t316 * t366 - t340 * t352 - t341 * t351) * MDP(20) + (t297 * t340 + t298 * t341 + t303 * t315 + t304 * t316 + t384 * t446 + t399 * t404) * MDP(21) + (t301 * t365 - t308 * t497) * MDP(22) + (-t301 * t364 - t302 * t365 - t308 * t335 + t309 * t497) * MDP(23) + (-t301 * t484 - t308 * t455 + (qJD(1) * t365 - t497) * t530) * MDP(24) + (t302 * t484 + t309 * t455 + (-qJD(1) * t364 - t335) * t530) * MDP(25) + (-t455 - t532) * MDP(26) * t530 + (-(t299 * t482 - t300 * t479) * t455 - t507 * t484 + t342 * t335 + t391 * t302 + t324 * t364 + t343 * t309 + (t293 * t484 + t455 * t498) * qJD(5) + ((t323 * t482 - t325 * t479) * qJD(1) + t292) * t530) * MDP(27) + (t391 * t301 - t306 * t484 + t343 * t308 + t324 * t365 - t342 * t497 + ((-qJD(5) * t325 + t299) * t455 + t295 * t484) * t479 + ((qJD(5) * t323 + t300) * t455 + (qJD(5) * t305 + t296) * t484) * t482 + (-qJD(1) * t498 - t293) * t530) * MDP(28) + MDP(6) * t549 - MDP(7) * t553 + 0.2e1 * t512 * t573; -t548 * t573 + t486 * t572 + (-t443 * t559 + t564) * MDP(11) + ((t401 + t563) * t483 + (-t402 + t562) * t480) * MDP(12) + (-t463 * t526 + (t463 * t550 + t481 * t504) * qJD(1)) * MDP(13) + (t463 * t528 + (-t463 * t555 + t481 * t505) * qJD(1)) * MDP(14) + (-pkin(2) * t402 + t537 * t463 + (pkin(7) * t559 + t561) * qJD(3) + ((-pkin(7) * t531 - t393) * t481 + (-pkin(6) * t505 - t561) * t484) * qJD(1)) * MDP(16) + (-pkin(2) * t401 - t426 * t463 + (-pkin(7) * t463 * t480 + t560) * qJD(3) + (-t452 * t550 + (-pkin(7) * t523 + t394) * t481 + (t463 * t554 + t504 * t484) * pkin(6)) * qJD(1)) * MDP(17) + (t351 * t469 + t384 * t434 + t502 * t386 - t542 * t399 + t547 * t463) * MDP(18) + (t352 * t469 + t384 * t435 + t399 * t577 - t546 * t463 + t496 * t502) * MDP(19) + (-t297 * t435 - t298 * t434 - t315 * t577 + t316 * t542 - t351 * t398 - t352 * t397 + t386 * t546 + t496 * t547) * MDP(20) + (t297 * t397 + t298 * t398 - t315 * t547 - t316 * t546 + t384 * t469 + t399 * t502) * MDP(21) + (t301 * t383 - t497 * t545) * MDP(22) + (-t301 * t382 - t302 * t383 - t335 * t545 + t497 * t544) * MDP(23) + (t408 * t302 + t324 * t382 + t543 * t335 + t544 * t343) * MDP(27) + (t408 * t301 + t324 * t383 + t545 * t343 - t497 * t543) * MDP(28) + (-t545 * MDP(24) + t544 * MDP(25) + (t479 * t499 + t482 * t500) * MDP(27) + (-t479 * t500 + t482 * t499) * MDP(28)) * t455 + (t463 * MDP(15) + (qJD(2) * t397 - t315) * MDP(18) + (-qJD(2) * t398 + t316) * MDP(19) + (qJD(2) * t383 + t497) * MDP(24) + (-qJD(2) * t382 + t335) * MDP(25) + t455 * MDP(26) + ((t372 * t482 - t373 * t479) * qJD(2) - t292) * MDP(27) + (-(t372 * t479 + t373 * t482) * qJD(2) + t293) * MDP(28)) * t533 + (MDP(9) * t481 * t486 + MDP(10) * t548) * pkin(1); t443 * t441 * MDP(11) + (-t441 ^ 2 + t443 ^ 2) * MDP(12) + (t401 - t563) * MDP(13) + (-t402 - t562) * MDP(14) + qJD(1) * t510 + (-t394 * t463 - t443 * t452 + t488) * MDP(16) + (-t393 * t463 + t441 * t452 - t487) * MDP(17) + (t321 * t463 - t496 * t399 + (-t386 * t443 + t478 * t513) * pkin(3) + t297) * MDP(18) + (-t322 * t463 + t386 * t399 + (-t443 * t496 - t477 * t513) * pkin(3) - t298) * MDP(19) + ((-t351 * t477 - t352 * t478) * pkin(3) + (t316 + t321) * t496 + (-t315 + t322) * t386) * MDP(20) + (-t315 * t321 - t316 * t322 + (t297 * t478 + t298 * t477 - t399 * t443) * pkin(3)) * MDP(21) + (t301 - t566) * MDP(24) + (-t302 + t567) * MDP(25) + (t491 * t513 + (t310 * t482 - t311 * t479) * t455 - t354 * t335 + (t455 * t492 - t293) * qJD(5) + t581) * MDP(27) + (-t492 * t513 - (t310 * t479 + t311 * t482) * t455 + t354 * t497 + (t455 * t491 - t552) * qJD(5) + t584) * MDP(28) + t583; (-t463 * t496 + t351) * MDP(18) + (t386 * t463 + t352) * MDP(19) + (-t386 ^ 2 - t496 ^ 2) * MDP(20) + (t315 * t496 + t316 * t386 + t384) * MDP(21) + (t302 + t567) * MDP(27) + (t301 + t566) * MDP(28); (t519 - t566) * MDP(24) + (-t506 + t567) * MDP(25) + (-t293 * t455 + t581) * MDP(27) + (-t292 * t455 + t584) * MDP(28) + (-MDP(24) * t565 + t497 * MDP(25) - MDP(27) * t293 - MDP(28) * t552) * qJD(5) + t583;];
tauc = t1;
