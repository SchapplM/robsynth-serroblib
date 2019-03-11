% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:40
% EndTime: 2019-03-09 03:59:48
% DurationCPUTime: 4.63s
% Computational Cost: add. (3762->403), mult. (8422->546), div. (0->0), fcn. (6130->8), ass. (0->182)
t405 = sin(qJ(3));
t408 = cos(qJ(3));
t495 = sin(pkin(10));
t496 = cos(pkin(10));
t371 = -t495 * t405 + t496 * t408;
t366 = t371 * qJD(1);
t404 = sin(qJ(5));
t407 = cos(qJ(5));
t451 = t407 * qJD(3);
t339 = t366 * t404 - t451;
t416 = t405 * t496 + t408 * t495;
t500 = t416 * qJD(1);
t508 = qJD(5) + t500;
t513 = t339 * t508;
t406 = cos(qJ(6));
t341 = qJD(3) * t404 + t366 * t407;
t403 = sin(qJ(6));
t484 = t341 * t403;
t285 = t406 * t339 + t484;
t353 = qJD(6) + t508;
t512 = t285 * t353;
t424 = t339 * t403 - t406 * t341;
t511 = t353 * t424;
t455 = qJD(5) * t404;
t479 = t500 * t404;
t510 = t455 + t479;
t354 = t366 * qJD(3);
t433 = t508 * t407;
t509 = -t354 * t404 - t433 * t508;
t409 = -pkin(1) - pkin(7);
t381 = qJD(1) * t409 + qJD(2);
t460 = qJD(1) * t408;
t359 = -qJ(4) * t460 + t408 * t381;
t349 = qJD(3) * pkin(3) + t359;
t461 = qJD(1) * t405;
t358 = -qJ(4) * t461 + t381 * t405;
t440 = t496 * t358;
t306 = t495 * t349 + t440;
t297 = qJD(3) * pkin(8) + t306;
t376 = pkin(3) * t461 + qJD(1) * qJ(2) + qJD(4);
t307 = pkin(4) * t500 - pkin(8) * t366 + t376;
t276 = t297 * t407 + t307 * t404;
t269 = -pkin(9) * t339 + t276;
t453 = qJD(6) * t403;
t266 = t269 * t453;
t346 = t495 * t358;
t305 = t349 * t496 - t346;
t296 = -qJD(3) * pkin(4) - t305;
t279 = t339 * pkin(5) + t296;
t507 = t279 * t285 + t266;
t355 = qJD(3) * t500;
t298 = qJD(5) * t451 - t407 * t355 - t366 * t455;
t457 = qJD(4) * t405;
t458 = qJD(3) * t408;
t338 = t381 * t458 + (-qJ(4) * t458 - t457) * qJD(1);
t456 = qJD(4) * t408;
t459 = qJD(3) * t405;
t412 = -t381 * t459 + (qJ(4) * t459 - t456) * qJD(1);
t283 = t338 * t496 + t412 * t495;
t399 = qJD(1) * qJD(2);
t448 = qJD(1) * qJD(3);
t442 = t408 * t448;
t465 = pkin(3) * t442 + t399;
t294 = pkin(4) * t354 + pkin(8) * t355 + t465;
t290 = t407 * t294;
t413 = -qJD(5) * t276 - t283 * t404 + t290;
t257 = pkin(5) * t354 - pkin(9) * t298 + t413;
t480 = t355 * t404;
t299 = t341 * qJD(5) - t480;
t454 = qJD(5) * t407;
t419 = t407 * t283 + t404 * t294 - t297 * t455 + t307 * t454;
t260 = -pkin(9) * t299 + t419;
t435 = t406 * t257 - t403 * t260;
t506 = t279 * t424 + t435;
t505 = (-t285 ^ 2 + t424 ^ 2) * MDP(24) + t354 * MDP(27) - t285 * MDP(23) * t424;
t504 = MDP(8) * (t405 ^ 2 - t408 ^ 2);
t474 = t403 * t407;
t374 = t404 * t406 + t474;
t321 = t374 * t371;
t499 = qJD(5) + qJD(6);
t332 = t499 * t374;
t373 = t403 * t404 - t406 * t407;
t502 = -t332 * t353 - t373 * t354;
t501 = MDP(12) * qJ(2);
t331 = t499 * t373;
t467 = -t373 * t500 - t331;
t482 = t354 * t374;
t498 = -t353 * t467 - t482;
t434 = t298 * t403 + t406 * t299;
t264 = -qJD(6) * t424 + t434;
t390 = pkin(3) * t495 + pkin(8);
t497 = pkin(9) + t390;
t411 = qJD(1) ^ 2;
t494 = qJ(2) * t411;
t275 = -t297 * t404 + t407 * t307;
t268 = -pkin(9) * t341 + t275;
t265 = pkin(5) * t508 + t268;
t493 = t265 * t406;
t492 = t269 * t406;
t282 = t338 * t495 - t496 * t412;
t491 = t282 * t371;
t490 = t285 * t366;
t489 = t424 * t366;
t488 = t298 * t404;
t470 = qJ(4) - t409;
t377 = t470 * t405;
t378 = t470 * t408;
t334 = -t377 * t496 - t378 * t495;
t487 = t334 * t354;
t486 = t339 * t366;
t485 = t341 * t366;
t436 = qJD(3) * t495;
t437 = qJD(3) * t496;
t364 = t405 * t436 - t408 * t437;
t483 = t353 * t364;
t365 = -t405 * t437 - t408 * t436;
t478 = t365 * t404;
t477 = t365 * t407;
t476 = t371 * t404;
t475 = t371 * t407;
t328 = t416 * t354;
t327 = t407 * t334;
t344 = t407 * t354;
t473 = t408 * t411;
t410 = qJD(3) ^ 2;
t472 = t409 * t410;
t471 = t405 * pkin(3) + qJ(2);
t317 = t359 * t496 - t346;
t319 = pkin(3) * t460 + pkin(4) * t366 + pkin(8) * t500;
t469 = t407 * t317 + t404 * t319;
t326 = pkin(4) * t416 - pkin(8) * t371 + t471;
t468 = t404 * t326 + t327;
t314 = t374 * t500;
t466 = t332 + t314;
t452 = qJD(6) * t406;
t450 = pkin(3) * t458 + qJD(2);
t447 = 0.2e1 * qJD(1);
t445 = t406 * t298 - t403 * t299 - t339 * t452;
t444 = t371 * t455;
t443 = t371 * t454;
t441 = qJD(5) * t497;
t432 = qJD(6) * t265 + t260;
t431 = -qJD(5) * t416 - qJD(1);
t391 = -pkin(3) * t496 - pkin(4);
t316 = t359 * t495 + t440;
t430 = pkin(5) * t510 - t316;
t427 = -t314 * t353 + t502;
t313 = t407 * t319;
t370 = t497 * t407;
t426 = pkin(5) * t366 + qJD(6) * t370 - t317 * t404 + t313 + (pkin(9) * t500 + t441) * t407;
t369 = t497 * t404;
t425 = pkin(9) * t479 + qJD(6) * t369 + t404 * t441 + t469;
t356 = t459 * t470 - t456;
t357 = -qJD(3) * t378 - t457;
t310 = -t496 * t356 + t357 * t495;
t333 = -t377 * t495 + t496 * t378;
t259 = t265 * t403 + t492;
t423 = -t508 * t510 + t344;
t422 = t443 + t478;
t421 = -t444 + t477;
t420 = t353 * t374;
t311 = t356 * t495 + t357 * t496;
t318 = -pkin(4) * t364 - pkin(8) * t365 + t450;
t418 = t407 * t311 + t404 * t318 + t326 * t454 - t334 * t455;
t263 = -t341 * t453 + t445;
t417 = t296 * t508 - t390 * t354;
t414 = t283 * t416 + t305 * t365 - t306 * t364 - t491;
t379 = -t407 * pkin(5) + t391;
t324 = t407 * t326;
t322 = t373 * t371;
t309 = t407 * t318;
t303 = pkin(5) * t476 + t333;
t278 = pkin(5) * t422 + t310;
t277 = -pkin(9) * t476 + t468;
t274 = pkin(5) * t416 - pkin(9) * t475 - t334 * t404 + t324;
t272 = t299 * pkin(5) + t282;
t271 = t365 * t474 - t403 * t444 - t453 * t476 + (t475 * t499 + t478) * t406;
t270 = -t321 * t499 - t373 * t365;
t262 = -pkin(9) * t422 + t418;
t261 = -pkin(9) * t477 - pkin(5) * t364 - t311 * t404 + t309 + (-t327 + (pkin(9) * t371 - t326) * t404) * qJD(5);
t258 = -t269 * t403 + t493;
t1 = [-t410 * t408 * MDP(10) + 0.2e1 * t448 * t504 + t458 * t447 * t501 + (-t408 * t472 + (-qJ(2) * t459 + qJD(2) * t408) * t447) * MDP(13) + (t310 * t366 - t311 * t500 - t333 * t355 - t414 - t487) * MDP(14) + (t282 * t333 + t283 * t334 - t305 * t310 + t306 * t311 + t376 * t450 + t465 * t471) * MDP(15) + (t298 * t475 + t341 * t421) * MDP(16) + ((-t339 * t407 - t341 * t404) * t365 + (-t488 - t299 * t407 + (t339 * t404 - t341 * t407) * qJD(5)) * t371) * MDP(17) + (t298 * t416 - t341 * t364 + t344 * t371 + t421 * t508) * MDP(18) + (-t299 * t416 + t339 * t364 - t354 * t476 - t422 * t508) * MDP(19) + (-t364 * t508 + t328) * MDP(20) + ((-t334 * t454 + t309) * t508 + t324 * t354 + (-t297 * t454 + t290) * t416 - t275 * t364 + t310 * t339 + t333 * t299 + t296 * t443 + ((-qJD(5) * t326 - t311) * t508 - t487 + (-qJD(5) * t307 - t283) * t416 + t491 + t296 * t365) * t404) * MDP(21) + (t276 * t364 + t282 * t475 + t296 * t421 + t333 * t298 + t310 * t341 - t354 * t468 - t416 * t419 - t418 * t508) * MDP(22) + (-t263 * t322 - t270 * t424) * MDP(23) + (-t263 * t321 + t264 * t322 - t270 * t285 + t271 * t424) * MDP(24) + (t263 * t416 + t270 * t353 - t322 * t354 + t364 * t424) * MDP(25) + (-t264 * t416 - t271 * t353 + t285 * t364 - t321 * t354) * MDP(26) + (t328 - t483) * MDP(27) + ((t261 * t406 - t262 * t403) * t353 + (t274 * t406 - t277 * t403) * t354 + t435 * t416 - t258 * t364 + t278 * t285 + t303 * t264 + t272 * t321 + t279 * t271 + ((-t274 * t403 - t277 * t406) * t353 - t259 * t416) * qJD(6)) * MDP(28) + (t259 * t364 + t303 * t263 + t266 * t416 + t279 * t270 - t272 * t322 - t278 * t424 + (-(-qJD(6) * t277 + t261) * t353 - t274 * t354 - t257 * t416) * t403 + (-(qJD(6) * t274 + t262) * t353 - t277 * t354 - t432 * t416) * t406) * MDP(29) + 0.2e1 * (MDP(6) * qJ(2) + MDP(5)) * t399 + (-t410 * MDP(9) + (qJD(2) * t447 - t472) * MDP(12) - 0.2e1 * MDP(7) * t442) * t405; -t411 * MDP(5) - MDP(6) * t494 + (t355 * t371 + t364 * t500 - t365 * t366 - t328) * MDP(14) + (-qJD(1) * t376 + t414) * MDP(15) + (-t404 * t328 - t299 * t371 - t339 * t365 + (t364 * t404 + t407 * t431) * t508) * MDP(21) + (-t416 * t344 - t298 * t371 - t341 * t365 + (t364 * t407 - t404 * t431) * t508) * MDP(22) + (-t371 * t264 - t365 * t285 + t364 * t420 + t373 * t353 * qJD(1) - (-t331 * t353 + t482) * t416) * MDP(28) + (qJD(1) * t420 - t371 * t263 + t365 * t424 - t373 * t483 - t416 * t502) * MDP(29) + (t405 * MDP(12) + t408 * MDP(13)) * (-t410 - t411); -t411 * t504 - t473 * t501 + (-(t305 - t317) * t500 + (-t354 * t495 + t355 * t496) * pkin(3)) * MDP(14) + (t305 * t316 - t306 * t317 + (-t282 * t496 + t283 * t495 - t376 * t460) * pkin(3)) * MDP(15) + (t341 * t433 + t488) * MDP(16) + ((t298 - t513) * t407 + (-t341 * t508 - t299) * t404) * MDP(17) + (-t485 - t509) * MDP(18) + (t423 + t486) * MDP(19) + (-t282 * t407 + t391 * t299 - t316 * t339 + (-t390 * t454 - t313) * t508 + (t317 * t508 + t417) * t404) * MDP(21) + (t282 * t404 + t391 * t298 - t316 * t341 + (t390 * t455 + t469) * t508 + t417 * t407) * MDP(22) + (t263 * t374 - t424 * t467) * MDP(23) + (-t263 * t373 - t264 * t374 - t285 * t467 + t424 * t466) * MDP(24) + (t489 - t498) * MDP(25) + (t427 + t490) * MDP(26) + ((-t369 * t406 - t370 * t403) * t354 + t379 * t264 + t272 * t373 + (t403 * t425 - t406 * t426) * t353 + t430 * t285 + t466 * t279) * MDP(28) + (-(-t369 * t403 + t370 * t406) * t354 + t379 * t263 + t272 * t374 + (t403 * t426 + t406 * t425) * t353 - t430 * t424 + t467 * t279) * MDP(29) + (MDP(13) * t494 + MDP(7) * t473) * t405 + ((t306 - t316) * MDP(14) - t508 * MDP(20) - t275 * MDP(21) + t276 * MDP(22) - t353 * MDP(27) - t258 * MDP(28) + t259 * MDP(29)) * t366; (-t366 ^ 2 - t500 ^ 2) * MDP(14) + (t305 * t366 + t306 * t500 + t465) * MDP(15) + (t423 - t486) * MDP(21) + (-t485 + t509) * MDP(22) + (t427 - t490) * MDP(28) + (t489 + t498) * MDP(29); t341 * t339 * MDP(16) + (-t339 ^ 2 + t341 ^ 2) * MDP(17) + (t298 + t513) * MDP(18) + (t480 + (-qJD(5) + t508) * t341) * MDP(19) + t354 * MDP(20) + (t276 * t508 - t296 * t341 + t413) * MDP(21) + (t275 * t508 + t296 * t339 - t419) * MDP(22) + (t263 + t512) * MDP(25) + (-t264 - t511) * MDP(26) + (-(-t268 * t403 - t492) * t353 - t259 * qJD(6) + (-t285 * t341 - t353 * t453 + t354 * t406) * pkin(5) + t506) * MDP(28) + ((-t269 * t353 - t257) * t403 + (t268 * t353 - t432) * t406 + (t341 * t424 - t353 * t452 - t354 * t403) * pkin(5) + t507) * MDP(29) + t505; (t445 + t512) * MDP(25) + (-t434 - t511) * MDP(26) + (t259 * t353 + t506) * MDP(28) + (-t403 * t257 + t258 * t353 - t406 * t260 + t507) * MDP(29) + (-MDP(25) * t484 + MDP(26) * t424 - MDP(28) * t259 - MDP(29) * t493) * qJD(6) + t505;];
tauc  = t1;
