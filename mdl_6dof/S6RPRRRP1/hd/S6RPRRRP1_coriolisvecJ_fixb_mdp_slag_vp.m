% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:15
% EndTime: 2019-03-09 05:57:23
% DurationCPUTime: 4.15s
% Computational Cost: add. (5692->401), mult. (12870->521), div. (0->0), fcn. (8731->8), ass. (0->175)
t406 = cos(qJ(5));
t468 = qJD(5) * t406;
t404 = sin(qJ(4));
t407 = cos(qJ(3));
t509 = cos(qJ(4));
t456 = qJD(1) * t509;
t405 = sin(qJ(3));
t473 = qJD(1) * t405;
t374 = t404 * t473 - t407 * t456;
t490 = t374 * t406;
t521 = t468 + t490;
t403 = sin(qJ(5));
t469 = qJD(5) * t403;
t491 = t374 * t403;
t520 = -t469 - t491;
t382 = t404 * t407 + t405 * t509;
t464 = qJD(3) + qJD(4);
t352 = t464 * t382;
t519 = t352 * qJD(1);
t392 = sin(pkin(10)) * pkin(1) + pkin(7);
t506 = pkin(8) + t392;
t373 = qJD(5) + t374;
t486 = t406 * t519;
t518 = pkin(9) * (t373 * t469 - t486);
t517 = MDP(6) * (t405 ^ 2 - t407 ^ 2);
t516 = t407 * MDP(5);
t452 = t506 * qJD(1);
t364 = t407 * qJD(2) - t405 * t452;
t355 = t364 * qJD(3);
t365 = qJD(2) * t405 + t407 * t452;
t356 = t365 * qJD(3);
t504 = qJD(3) * pkin(3);
t359 = t364 + t504;
t455 = qJD(4) * t509;
t470 = qJD(4) * t404;
t283 = t355 * t509 - t404 * t356 + t359 * t455 - t365 * t470;
t430 = -t404 * t405 + t407 * t509;
t419 = t430 * qJD(4);
t302 = t519 * pkin(4) + (-pkin(9) * t419 + (pkin(3) * t405 - pkin(9) * t430) * qJD(3)) * qJD(1);
t358 = t509 * t365;
t321 = t359 * t404 + t358;
t316 = pkin(9) * t464 + t321;
t472 = qJD(1) * t407;
t376 = -t404 * t472 - t405 * t456;
t393 = -cos(pkin(10)) * pkin(1) - pkin(2);
t383 = -pkin(3) * t407 + t393;
t377 = t383 * qJD(1);
t334 = pkin(4) * t374 + pkin(9) * t376 + t377;
t424 = t283 * t406 + t302 * t403 - t316 * t469 + t334 * t468;
t503 = qJ(6) * t519;
t264 = qJD(6) * t373 + t424 + t503;
t447 = t283 * t403 - t302 * t406 + t316 * t468 + t334 * t469;
t508 = pkin(5) * t519;
t266 = t447 - t508;
t440 = t264 * t406 + t266 * t403;
t351 = qJD(3) * t430 + t419;
t411 = t351 * qJD(1);
t449 = t406 * t464;
t310 = -qJD(5) * t449 - t376 * t469 - t406 * t411;
t422 = t406 * t376 - t403 * t464;
t481 = t310 * t430 - t352 * t422;
t311 = -qJD(5) * t422 + t403 * t411;
t360 = -t376 * t403 - t449;
t515 = -t311 * t430 + t352 * t360;
t327 = t364 * t404 + t358;
t446 = pkin(3) * t470 - t327;
t341 = -pkin(4) * t430 - pkin(9) * t382 + t383;
t379 = t506 * t405;
t380 = t506 * t407;
t344 = -t379 * t404 + t380 * t509;
t476 = t341 * t403 + t344 * t406;
t514 = pkin(5) * t520 + qJ(6) * t521 + qJD(6) * t403;
t513 = -t379 * t509 - t380 * t404;
t512 = MDP(24) + MDP(26);
t511 = t422 ^ 2;
t510 = t373 ^ 2;
t507 = pkin(5) * t376;
t505 = pkin(3) * qJD(4);
t284 = t355 * t404 + t356 * t509 + t359 * t470 + t365 * t455;
t268 = pkin(5) * t311 + qJ(6) * t310 + qJD(6) * t422 + t284;
t502 = t268 * t403;
t289 = t316 * t406 + t334 * t403;
t501 = t289 * t373;
t500 = t310 * t403;
t395 = pkin(3) * t404 + pkin(9);
t499 = t519 * t395;
t498 = t351 * t403;
t497 = t351 * t406;
t496 = t360 * t373;
t495 = t360 * t403;
t494 = t422 * t360;
t493 = t422 * t373;
t492 = t422 * t406;
t451 = t373 * t406;
t489 = t382 * t406;
t488 = t403 * t519;
t357 = t404 * t365;
t408 = qJD(3) ^ 2;
t487 = t405 * t408;
t485 = t407 * t408;
t484 = t321 + t514;
t483 = t514 - t446;
t482 = -t311 * t489 - t360 * t497;
t320 = t359 * t509 - t357;
t346 = -pkin(4) * t376 + pkin(9) * t374;
t479 = t320 * t406 + t346 * t403;
t328 = t364 * t509 - t357;
t338 = pkin(3) * t473 + t346;
t478 = t328 * t406 + t338 * t403;
t477 = t351 * t451 + t382 * t486;
t474 = MDP(27) * t403;
t386 = qJD(1) * t393;
t288 = -t316 * t403 + t334 * t406;
t467 = qJD(6) - t288;
t466 = qJD(1) * qJD(3);
t463 = t509 * pkin(3);
t461 = t405 * t504;
t460 = t403 * t509;
t459 = t406 * t509;
t457 = t382 * t469;
t454 = t405 * t466;
t453 = qJD(3) * t506;
t450 = t352 * t464;
t448 = pkin(3) * t455;
t315 = -pkin(4) * t464 - t320;
t444 = t284 * t403 - t289 * t376 + t315 * t468;
t443 = pkin(5) * t406 + qJ(6) * t403;
t442 = pkin(5) * t403 - qJ(6) * t406;
t277 = -pkin(5) * t373 + t467;
t278 = qJ(6) * t373 + t289;
t439 = t277 * t406 - t278 * t403;
t438 = t277 * t403 + t278 * t406;
t437 = t315 * t374 - t499;
t436 = -t320 * t403 + t346 * t406;
t435 = -t492 + t495;
t434 = 0.2e1 * qJD(3) * t386;
t384 = -pkin(4) - t443;
t290 = pkin(5) * t360 + qJ(6) * t422 + t315;
t433 = -t268 * t406 - t277 * t376 + t290 * t469;
t432 = t278 * t376 - t290 * t490 - t502;
t431 = -t284 * t406 + t288 * t376 + t315 * t469;
t428 = t382 * t468 + t498;
t427 = t457 - t497;
t426 = -t290 * t422 + t447;
t425 = t377 * t376 - t284;
t370 = t405 * t453;
t371 = t407 * t453;
t305 = qJD(4) * t513 - t370 * t509 - t371 * t404;
t309 = pkin(4) * t352 - pkin(9) * t351 + t461;
t423 = t305 * t406 + t309 * t403 + t341 * t468 - t344 * t469;
t421 = -MDP(25) * t406 - t403 * t512;
t420 = (-t373 * t468 - t488) * pkin(9);
t418 = -t395 * t469 + t406 * t448;
t417 = t277 * t521 + t278 * t520 + t440;
t416 = qJD(5) * t439 + t440;
t415 = qJD(5) * t435 - t311 * t406 - t500;
t414 = ((-t310 - t496) * t406 + (-t311 + t493) * t403) * MDP(20) + (-t422 * t451 - t500) * MDP(19) + (-t360 * t376 - t403 * t510 + t486) * MDP(22) + (t373 * t451 - t376 * t422 + t488) * MDP(21) + (t374 * t464 + t411) * MDP(14) + (-t376 * t464 - t519) * MDP(15) + (-t374 ^ 2 + t376 ^ 2) * MDP(13) + (-MDP(12) * t374 + MDP(23) * t373) * t376;
t413 = t377 * t374 - t283;
t306 = qJD(4) * t344 - t370 * t404 + t371 * t509;
t396 = -t463 - pkin(4);
t378 = -t463 + t384;
t369 = t376 * qJ(6);
t326 = -pkin(5) * t422 + qJ(6) * t360;
t308 = t382 * t442 - t513;
t296 = pkin(5) * t430 - t341 * t406 + t344 * t403;
t295 = -qJ(6) * t430 + t476;
t293 = -t310 + t496;
t292 = -t436 + t507;
t291 = -t369 + t479;
t287 = t328 * t403 - t338 * t406 + t507;
t286 = -t369 + t478;
t271 = t442 * t351 + (qJD(5) * t443 - qJD(6) * t406) * t382 + t306;
t270 = -pkin(5) * t352 + qJD(5) * t476 + t305 * t403 - t309 * t406;
t269 = qJ(6) * t352 - qJD(6) * t430 + t423;
t1 = [0.2e1 * t454 * t516 - 0.2e1 * t466 * t517 + MDP(7) * t485 - MDP(8) * t487 + (-t392 * t485 + t405 * t434) * MDP(10) + (t392 * t487 + t407 * t434) * MDP(11) + (-t376 * t351 + t382 * t411) * MDP(12) + (-t351 * t374 + t376 * t352 - t382 * t519 + t411 * t430) * MDP(13) + t351 * t464 * MDP(14) - MDP(15) * t450 + (t383 * t519 + t377 * t352 - t306 * t464 + (-qJD(1) * t430 + t374) * t461) * MDP(17) + (pkin(3) * t382 * t454 - t305 * t464 + t377 * t351 - t376 * t461 + t383 * t411) * MDP(18) + (-t310 * t489 + t422 * t427) * MDP(19) + (t422 * t498 + (t500 + (t492 + t495) * qJD(5)) * t382 + t482) * MDP(20) + (-t373 * t457 + t477 + t481) * MDP(21) + (-t373 * t428 - t382 * t488 - t515) * MDP(22) + (t352 * t373 - t430 * t519) * MDP(23) + (t447 * t430 + t288 * t352 + t306 * t360 - t513 * t311 + ((-qJD(5) * t344 + t309) * t373 + t341 * t519 + t315 * qJD(5) * t382) * t406 + ((-qJD(5) * t341 - t305) * t373 - t344 * t519 + t284 * t382 + t315 * t351) * t403) * MDP(24) + (t284 * t489 - t289 * t352 - t306 * t422 + t310 * t513 - t315 * t427 - t373 * t423 + t424 * t430 - t476 * t519) * MDP(25) + (t266 * t430 - t270 * t373 + t271 * t360 - t277 * t352 + t290 * t428 - t296 * t519 + t308 * t311 + t382 * t502) * MDP(26) + (-t269 * t360 - t270 * t422 - t295 * t311 - t296 * t310 + t439 * t351 + (-qJD(5) * t438 - t264 * t403 + t266 * t406) * t382) * MDP(27) + (-t264 * t430 - t268 * t489 + t269 * t373 + t271 * t422 + t278 * t352 + t290 * t427 + t295 * t519 + t308 * t310) * MDP(28) + (t264 * t295 + t266 * t296 + t268 * t308 + t269 * t278 + t270 * t277 + t271 * t290) * MDP(29); -MDP(17) * t450 + t481 * MDP(25) + t482 * MDP(27) + (t477 - t481) * MDP(28) + (-t268 * t430 + t290 * t352) * MDP(29) + (-MDP(10) * t405 - MDP(11) * t407) * t408 + (-MDP(18) * t464 + MDP(29) * t438 + t373 * t421 - t422 * t474) * t351 + (-t310 * t474 + t440 * MDP(29) + t421 * t519 + (t435 * MDP(27) + t439 * MDP(29) + (-t512 * t406 + (MDP(25) - MDP(28)) * t403) * t373) * qJD(5)) * t382 + t512 * t515; (t328 * t464 + (t376 * t473 - t455 * t464) * pkin(3) + t413) * MDP(18) + (t396 * t311 + t437 * t403 + t446 * t360 + ((-qJD(5) * t395 - t338) * t406 + (-t448 + t328) * t403) * t373 + t431) * MDP(24) + (-t396 * t310 + t437 * t406 - t446 * t422 + (-t418 + t478) * t373 + t444) * MDP(25) + (t327 * t464 + (-t374 * t473 - t464 * t470) * pkin(3) + t425) * MDP(17) + (t268 * t378 - t277 * t287 - t278 * t286 - t483 * t290 + (t277 * t460 + t278 * t459) * t505 + t416 * t395) * MDP(29) + (t378 * t311 + (t290 * t374 - t499) * t403 - t483 * t360 + (-t395 * t468 - t403 * t448 + t287) * t373 + t433) * MDP(26) + (t286 * t360 + t287 * t422 + (-t360 * t459 - t422 * t460) * t505 + t415 * t395 + t417) * MDP(27) + (t378 * t310 + (-qJD(5) * t290 + t499) * t406 - t483 * t422 + (-t286 + t418) * t373 + t432) * MDP(28) + t414 + (-t405 * t516 + t517) * qJD(1) ^ 2 + (-MDP(10) * t473 - MDP(11) * t472) * t386; (t320 * t464 + t413) * MDP(18) + (pkin(4) * t310 + t315 * t490 + t321 * t422 + t373 * t479 + t444 + t518) * MDP(25) + (t321 * t464 + t425) * MDP(17) + (pkin(9) * t415 + t291 * t360 + t292 * t422 + t417) * MDP(27) + (-t290 * t468 - t291 * t373 + t310 * t384 - t422 * t484 + t432 - t518) * MDP(28) + (t290 * t491 + t292 * t373 + t311 * t384 - t360 * t484 + t420 + t433) * MDP(26) + t414 + (-pkin(4) * t311 + t315 * t491 - t321 * t360 - t373 * t436 + t420 + t431) * MDP(24) + (pkin(9) * t416 + t268 * t384 - t277 * t292 - t278 * t291 - t290 * t484) * MDP(29); -MDP(19) * t494 + (-t360 ^ 2 + t511) * MDP(20) + t293 * MDP(21) + (-t311 - t493) * MDP(22) + t519 * MDP(23) + (t315 * t422 - t447 + t501) * MDP(24) + (t288 * t373 + t315 * t360 - t424) * MDP(25) + (-t326 * t360 - t426 + t501 + 0.2e1 * t508) * MDP(26) + (pkin(5) * t310 - qJ(6) * t311 - (t278 - t289) * t422 + (t277 - t467) * t360) * MDP(27) + (0.2e1 * t503 - t290 * t360 - t326 * t422 + (0.2e1 * qJD(6) - t288) * t373 + t424) * MDP(28) + (-pkin(5) * t266 + qJ(6) * t264 - t277 * t289 + t278 * t467 - t290 * t326) * MDP(29); t293 * MDP(27) + (-t510 - t511) * MDP(28) + (-t278 * t373 + t426 - t508) * MDP(29) + (-t494 - t519) * MDP(26);];
tauc  = t1;
