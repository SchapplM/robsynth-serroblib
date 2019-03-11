% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:24
% EndTime: 2019-03-09 09:15:34
% DurationCPUTime: 5.88s
% Computational Cost: add. (4106->408), mult. (9889->550), div. (0->0), fcn. (7248->8), ass. (0->184)
t449 = sin(pkin(10));
t450 = cos(pkin(10));
t457 = cos(qJ(2));
t503 = qJD(1) * t457;
t454 = sin(qJ(2));
t504 = qJD(1) * t454;
t387 = t449 * t504 + t450 * t503;
t390 = -t449 * t503 + t450 * t504;
t453 = sin(qJ(5));
t456 = cos(qJ(5));
t341 = t387 * t453 - t390 * t456;
t452 = sin(qJ(6));
t497 = qJD(6) * t452;
t444 = qJD(2) - qJD(5);
t455 = cos(qJ(6));
t492 = qJD(1) * qJD(2);
t489 = t454 * t492;
t422 = t450 * t489;
t488 = t457 * t492;
t380 = t449 * t488 - t422;
t404 = t449 * t454 + t450 * t457;
t396 = t404 * qJD(2);
t381 = qJD(1) * t396;
t498 = qJD(5) * t456;
t499 = qJD(5) * t453;
t465 = t453 * t380 - t456 * t381 + t387 * t498 + t390 * t499;
t496 = qJD(6) * t455;
t512 = -t444 * t496 - t455 * t465;
t290 = t341 * t497 + t512;
t326 = -t341 * t455 - t444 * t452;
t522 = t465 * t452;
t291 = t326 * qJD(6) - t522;
t306 = -qJD(5) * t341 + t456 * t380 + t381 * t453;
t304 = t455 * t306;
t517 = t452 * t306;
t533 = t456 * t387 + t390 * t453;
t544 = qJD(6) + t533;
t520 = t326 * t544;
t519 = t341 * t452;
t324 = t455 * t444 - t519;
t521 = t324 * t544;
t525 = t290 * t452;
t545 = t444 * t533;
t547 = t341 * t444;
t548 = t326 * t341;
t549 = t324 * t341;
t554 = t544 * t455;
t555 = t544 ^ 2;
t563 = -((t291 + t520) * t452 - (t290 - t521) * t455) * MDP(27) + (t326 * t554 + t525) * MDP(26) + (t544 * t554 + t517 + t548) * MDP(28) - (t452 * t555 - t304 + t549) * MDP(29) + (t341 ^ 2 - t533 ^ 2) * MDP(20) - (t465 + t545) * MDP(21) + (-t306 + t547) * MDP(22) + (-MDP(19) * t533 + MDP(30) * t544) * t341;
t551 = pkin(5) * t341;
t436 = pkin(7) * t504;
t411 = qJ(4) * t504 - t436;
t458 = -pkin(2) - pkin(3);
t490 = t458 * qJD(2);
t376 = qJD(3) + t490 - t411;
t437 = pkin(7) * t503;
t413 = -qJ(4) * t503 + t437;
t446 = qJD(2) * qJ(3);
t400 = t413 + t446;
t332 = t450 * t376 - t400 * t449;
t528 = pkin(8) * t390;
t318 = -qJD(2) * pkin(4) + t332 - t528;
t333 = t449 * t376 + t450 * t400;
t529 = pkin(8) * t387;
t319 = t333 - t529;
t295 = t318 * t453 + t319 * t456;
t293 = -pkin(9) * t444 + t295;
t401 = -qJD(1) * pkin(1) - pkin(2) * t503 - qJ(3) * t504;
t366 = pkin(3) * t503 + qJD(4) - t401;
t340 = pkin(4) * t387 + t366;
t296 = pkin(5) * t533 + pkin(9) * t341 + t340;
t284 = t293 * t455 + t296 * t452;
t550 = t284 * t341;
t473 = t293 * t452 - t296 * t455;
t546 = t341 * t473;
t502 = qJD(2) * t454;
t527 = pkin(7) - qJ(4);
t383 = -qJD(4) * t457 - t502 * t527;
t445 = qJD(2) * qJD(3);
t362 = qJD(1) * t383 + t445;
t429 = pkin(7) * t488;
t500 = qJD(4) * t454;
t501 = qJD(2) * t457;
t370 = t429 + (-qJ(4) * t501 - t500) * qJD(1);
t484 = t362 * t449 - t450 * t370;
t313 = -pkin(8) * t381 - t484;
t323 = t450 * t362 + t449 * t370;
t314 = -pkin(8) * t380 + t323;
t281 = t453 * t313 + t456 * t314 + t318 * t498 - t319 * t499;
t542 = t340 * t533 - t281;
t282 = -t313 * t456 + t453 * t314 + t318 * t499 + t319 * t498;
t541 = t340 * t341 - t282;
t535 = -0.2e1 * t492;
t447 = t454 ^ 2;
t534 = MDP(5) * (-t457 ^ 2 + t447);
t434 = qJ(3) * t503;
t382 = t458 * t504 + t434;
t348 = -pkin(4) * t390 + t382;
t414 = -qJ(3) * t449 + t450 * t458;
t410 = -pkin(4) + t414;
t415 = qJ(3) * t450 + t449 * t458;
t469 = t410 * t453 + t415 * t456;
t352 = -pkin(9) + t469;
t532 = t544 * (-pkin(9) * t533 + qJD(6) * t352 + t348 + t551) - t282;
t531 = t544 * (t544 * pkin(9) - t551) + t282;
t421 = t527 * t457;
t384 = qJD(2) * t421 - t500;
t334 = -t383 * t449 + t450 * t384;
t320 = -pkin(8) * t396 + t334;
t335 = t450 * t383 + t449 * t384;
t395 = t449 * t501 - t450 * t502;
t321 = -pkin(8) * t395 + t335;
t420 = t527 * t454;
t358 = t450 * t420 - t421 * t449;
t405 = -t449 * t457 + t450 * t454;
t336 = -pkin(8) * t405 + t358;
t359 = t449 * t420 + t450 * t421;
t337 = -pkin(8) * t404 + t359;
t471 = t336 * t456 - t337 * t453;
t287 = qJD(5) * t471 + t320 * t453 + t321 * t456;
t294 = t318 * t456 - t319 * t453;
t292 = pkin(5) * t444 - t294;
t349 = t456 * t404 + t405 * t453;
t350 = -t404 * t453 + t405 * t456;
t418 = -t457 * pkin(2) - t454 * qJ(3) - pkin(1);
t402 = t457 * pkin(3) - t418;
t355 = pkin(4) * t404 + t402;
t300 = pkin(5) * t349 - pkin(9) * t350 + t355;
t302 = t336 * t453 + t337 * t456;
t311 = -qJD(5) * t349 - t395 * t453 + t396 * t456;
t530 = -(qJD(6) * t300 + t287) * t544 - t302 * t306 - (qJD(6) * t296 + t281) * t349 + t282 * t350 + t292 * t311;
t526 = qJD(2) * pkin(2);
t524 = t292 * t350;
t523 = t300 * t306;
t459 = qJD(2) ^ 2;
t516 = t454 * t459;
t515 = t457 * t459;
t460 = qJD(1) ^ 2;
t514 = t457 * t460;
t353 = -t411 * t449 + t450 * t413;
t327 = t353 - t529;
t354 = t450 * t411 + t449 * t413;
t328 = t354 + t528;
t403 = t449 * t453 - t450 * t456;
t470 = t410 * t456 - t415 * t453;
t513 = qJD(3) * t403 - qJD(5) * t470 + t327 * t453 + t328 * t456;
t406 = t449 * t456 + t450 * t453;
t511 = qJD(3) * t406 + qJD(5) * t469 + t327 * t456 - t328 * t453;
t510 = t444 * t403;
t509 = t444 * t406;
t440 = t454 * qJD(3);
t508 = qJ(3) * t488 + qJD(1) * t440;
t507 = qJ(3) * t501 + t440;
t491 = t454 * t514;
t486 = pkin(1) * t535;
t485 = qJD(3) - t526;
t480 = qJD(1) * t418 + t401;
t479 = qJD(3) * t449 + t353;
t478 = qJD(3) * t450 - t354;
t474 = t454 * t490;
t472 = t332 * t449 - t333 * t450;
t468 = qJD(6) * t406 + t504;
t371 = pkin(2) * t489 - t508;
t385 = pkin(2) * t502 - t507;
t467 = -pkin(7) * t459 - qJD(1) * t385 - t371;
t466 = t311 * t455 - t350 * t497;
t365 = t474 + t507;
t357 = qJD(1) * t474 + t508;
t339 = pkin(4) * t395 + t365;
t463 = -pkin(9) * t306 + (t292 + t294) * t544;
t331 = pkin(4) * t380 + t357;
t462 = -t352 * t306 + (-t292 + t513) * t544;
t416 = -pkin(7) * t489 + t445;
t417 = t436 + t485;
t419 = t437 + t446;
t461 = t416 * t457 + (t417 * t457 + (-t419 + t437) * t454) * qJD(2);
t412 = pkin(2) * t504 - t434;
t351 = pkin(5) - t470;
t312 = qJD(5) * t350 + t456 * t395 + t396 * t453;
t289 = pkin(5) * t312 - pkin(9) * t311 + t339;
t288 = qJD(5) * t302 - t320 * t456 + t321 * t453;
t286 = pkin(5) * t306 + pkin(9) * t465 + t331;
t285 = t455 * t286;
t1 = [((-t324 * t455 - t326 * t452) * t311 + (-t525 - t291 * t455 + (t324 * t452 - t326 * t455) * qJD(6)) * t350) * MDP(27) + (-pkin(7) * t515 + t454 * t486) * MDP(9) - MDP(7) * t516 + (pkin(7) * t516 + t457 * t486) * MDP(10) + (-MDP(21) * t311 + MDP(22) * t312 + MDP(24) * t288 + MDP(25) * t287) * t444 + (-t350 * t517 - t291 * t349 - t312 * t324 + (-t311 * t452 - t350 * t496) * t544) * MDP(29) + (t290 * t349 + t304 * t350 + t312 * t326 + t466 * t544) * MDP(28) + (-t284 * t312 + t288 * t326 - t471 * t290 + (-(-qJD(6) * t302 + t289) * t544 - t523 - (-qJD(6) * t293 + t286) * t349 - qJD(6) * t524) * t452 + t530 * t455) * MDP(32) + (-t473 * t312 + t285 * t349 + t288 * t324 - t471 * t291 + (t289 * t544 + t523 + (-t293 * t349 - t302 * t544 + t524) * qJD(6)) * t455 + t530 * t452) * MDP(31) + (t306 * t349 + t312 * t544) * MDP(30) + (t457 * t467 + t480 * t502) * MDP(11) + (t454 * t467 - t480 * t501) * MDP(13) + (t290 * t350 * t455 + t326 * t466) * MDP(26) + (pkin(7) * t461 + t371 * t418 + t385 * t401) * MDP(14) + t461 * MDP(12) + (t323 * t359 + t332 * t334 + t333 * t335 + t357 * t402 - t358 * t484 + t365 * t366) * MDP(18) + (-t323 * t404 - t332 * t396 - t333 * t395 - t334 * t390 - t335 * t387 - t358 * t381 - t359 * t380 + t405 * t484) * MDP(17) + (t306 * t355 + t312 * t340 + t331 * t349 + t339 * t533) * MDP(24) + (t311 * t340 + t331 * t350 - t339 * t341 - t355 * t465) * MDP(25) + (-t311 * t341 - t350 * t465) * MDP(19) + (-t306 * t350 - t311 * t533 + t312 * t341 + t349 * t465) * MDP(20) + 0.2e1 * t454 * MDP(4) * t488 + t534 * t535 + MDP(6) * t515 + (-qJD(2) * t334 + t357 * t404 + t365 * t387 + t366 * t395 + t380 * t402) * MDP(15) + (qJD(2) * t335 + t357 * t405 + t365 * t390 + t366 * t396 + t381 * t402) * MDP(16); (qJD(2) * t478 - t366 * t387 - t382 * t390 + t323) * MDP(16) - MDP(4) * t491 + (qJD(2) * t479 + t366 * t390 - t382 * t387 + t484) * MDP(15) + (-t380 * t415 - t381 * t414 + (-t333 + t479) * t390 + (t332 - t478) * t387) * MDP(17) + 0.2e1 * t445 * MDP(13) + (qJ(3) * t416 + qJD(3) * t419 - t401 * t412) * MDP(14) + (t351 * t290 + t511 * t326 + t532 * t452 + t462 * t455 + t550) * MDP(32) + (t351 * t291 + t511 * t324 + t462 * t452 - t532 * t455 + t546) * MDP(31) + (-t348 * t533 + t444 * t511 - t541) * MDP(24) + (t341 * t348 - t444 * t513 - t542) * MDP(25) + (-qJD(3) * t472 + t323 * t415 - t332 * t353 - t333 * t354 - t366 * t382 - t414 * t484) * MDP(18) + ((t419 * t454 + (-t417 - t526) * t457) * pkin(7) * MDP(14) + (-t401 * t454 + t412 * t457) * MDP(11) + (t401 * t457 + t412 * t454) * MDP(13) + ((t419 - t446) * t454 + (-t417 + t485) * t457) * MDP(12)) * qJD(1) + t460 * t534 + (MDP(9) * t454 * t460 + MDP(10) * t514) * pkin(1) - t563; -MDP(11) * t491 + (-t447 * t460 - t459) * MDP(13) + (-qJD(2) * t419 + t401 * t504 + t429) * MDP(14) + (-t387 * t504 - t449 * t459) * MDP(15) + (-t390 * t504 - t450 * t459) * MDP(16) + (-t380 * t449 - t381 * t450 + (t387 * t450 - t390 * t449) * qJD(2)) * MDP(17) + (qJD(2) * t472 + t323 * t449 - t366 * t504 - t450 * t484) * MDP(18) + (-t444 * t509 - t504 * t533) * MDP(24) + (t341 * t504 + t444 * t510) * MDP(25) + (-t406 * t517 + t403 * t291 - t509 * t324 + (-t452 * t510 - t455 * t468) * t544) * MDP(31) + (-t406 * t304 + t403 * t290 - t509 * t326 + (t452 * t468 - t455 * t510) * t544) * MDP(32); -t422 * MDP(15) + (-t387 ^ 2 - t390 ^ 2) * MDP(17) + (t332 * t390 + t333 * t387 + t508) * MDP(18) + (t306 + t547) * MDP(24) + (-t465 + t545) * MDP(25) + (t304 + t549) * MDP(31) + (-t517 + t548) * MDP(32) + (-t390 * MDP(15) + t387 * MDP(16) + ((MDP(15) * t449 + MDP(16) * t450) * t457 + (t449 * MDP(16) + MDP(18) * t458) * t454) * qJD(1)) * qJD(2) - (t452 * MDP(31) + t455 * MDP(32)) * t555; (-t295 * t444 + t541) * MDP(24) + (-t294 * t444 + t542) * MDP(25) + (-pkin(5) * t291 - t295 * t324 + t463 * t452 - t531 * t455 - t546) * MDP(31) + (-pkin(5) * t290 - t295 * t326 + t531 * t452 + t463 * t455 - t550) * MDP(32) + t563; t326 * t324 * MDP(26) + (-t324 ^ 2 + t326 ^ 2) * MDP(27) + (t512 + t521) * MDP(28) + (t520 + t522) * MDP(29) + t306 * MDP(30) + (-t281 * t452 + t284 * t544 - t292 * t326 + t285) * MDP(31) + (-t281 * t455 - t286 * t452 + t292 * t324 - t473 * t544) * MDP(32) + (MDP(28) * t519 - MDP(29) * t326 - MDP(31) * t284 + MDP(32) * t473) * qJD(6);];
tauc  = t1;
