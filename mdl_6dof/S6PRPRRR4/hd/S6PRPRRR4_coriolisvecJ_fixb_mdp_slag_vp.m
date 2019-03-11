% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:39:04
% EndTime: 2019-03-08 20:39:15
% DurationCPUTime: 5.19s
% Computational Cost: add. (3674->406), mult. (9532->551), div. (0->0), fcn. (7819->12), ass. (0->189)
t423 = cos(pkin(12));
t519 = cos(qJ(4));
t469 = t519 * t423;
t412 = qJD(2) * t469;
t421 = sin(pkin(12));
t427 = sin(qJ(4));
t494 = t427 * t421;
t466 = qJD(2) * t494;
t388 = -t412 + t466;
t522 = qJD(5) + qJD(6);
t539 = t388 + t522;
t441 = t469 - t494;
t518 = pkin(8) + qJ(3);
t404 = t518 * t421;
t405 = t518 * t423;
t523 = -t519 * t404 - t427 * t405;
t328 = t441 * qJD(3) + t523 * qJD(4);
t422 = sin(pkin(6));
t431 = cos(qJ(2));
t497 = t422 * t431;
t434 = t441 * t497;
t368 = qJD(1) * t434;
t538 = -t328 + t368;
t391 = t441 * qJD(4);
t397 = t519 * t421 + t427 * t423;
t392 = t397 * qJD(4);
t428 = sin(qJ(2));
t482 = qJD(1) * t422;
t468 = t428 * t482;
t537 = pkin(4) * t392 - pkin(9) * t391 - t468;
t403 = qJD(2) * qJ(3) + t468;
t424 = cos(pkin(6));
t481 = qJD(1) * t424;
t411 = t423 * t481;
t517 = pkin(8) * qJD(2);
t360 = t411 + (-t403 - t517) * t421;
t373 = t423 * t403 + t421 * t481;
t361 = t423 * t517 + t373;
t307 = t427 * t360 + t519 * t361;
t536 = qJD(4) * t307;
t524 = t519 * t360 - t427 * t361;
t302 = -qJD(4) * pkin(4) - t524;
t390 = qJD(2) * t397;
t426 = sin(qJ(5));
t430 = cos(qJ(5));
t476 = t430 * qJD(4);
t369 = t390 * t426 - t476;
t297 = t369 * pkin(5) + t302;
t429 = cos(qJ(6));
t371 = qJD(4) * t426 + t390 * t430;
t425 = sin(qJ(6));
t506 = t371 * t425;
t311 = t429 * t369 + t506;
t535 = t297 * t311;
t385 = qJD(5) + t388;
t381 = qJD(6) + t385;
t534 = t311 * t381;
t447 = t369 * t425 - t429 * t371;
t533 = t381 * t447;
t496 = t425 * t430;
t400 = t426 * t429 + t496;
t485 = t539 * t400;
t480 = qJD(5) * t426;
t503 = t388 * t426;
t532 = t480 + t503;
t409 = qJD(4) * t412;
t382 = -qJD(4) * t466 + t409;
t324 = qJD(5) * t476 + t430 * t382 - t390 * t480;
t383 = qJD(2) * t392;
t303 = qJD(4) * pkin(9) + t307;
t414 = -pkin(3) * t423 - pkin(2);
t467 = t431 * t482;
t451 = qJD(3) - t467;
t384 = t414 * qJD(2) + t451;
t315 = pkin(4) * t388 - pkin(9) * t390 + t384;
t289 = t303 * t430 + t315 * t426;
t395 = (qJD(3) + t467) * qJD(2);
t529 = t441 * t395;
t294 = qJD(4) * t524 + t529;
t498 = t422 * t428;
t465 = qJD(2) * t498;
t408 = qJD(1) * t465;
t323 = pkin(4) * t383 - pkin(9) * t382 + t408;
t317 = t430 * t323;
t433 = -qJD(5) * t289 - t294 * t426 + t317;
t272 = pkin(5) * t383 - pkin(10) * t324 + t433;
t325 = qJD(5) * t371 + t382 * t426;
t479 = qJD(5) * t430;
t437 = t430 * t294 - t303 * t480 + t315 * t479 + t426 * t323;
t273 = -pkin(10) * t325 + t437;
t461 = t429 * t272 - t425 * t273;
t531 = t297 * t447 + t461;
t530 = t383 * MDP(27) + (-t311 ^ 2 + t447 ^ 2) * MDP(24) - t311 * MDP(23) * t447;
t342 = t400 * t397;
t446 = (-t403 * t421 + t411) * t421 - t373 * t423;
t528 = t446 * t431;
t484 = t421 ^ 2 + t423 ^ 2;
t527 = t484 * MDP(7);
t366 = -t427 * t404 + t519 * t405;
t435 = t397 * t497;
t488 = -qJD(1) * t435 + qJD(3) * t397 + qJD(4) * t366;
t526 = t368 * t426 + t537 * t430;
t349 = -pkin(4) * t441 - pkin(9) * t397 + t414;
t525 = -t349 * t479 + t366 * t480 - t537 * t426 + t538 * t430;
t399 = t425 * t426 - t429 * t430;
t486 = t539 * t399;
t521 = t486 * t381 - t383 * t400;
t459 = t324 * t425 + t429 * t325;
t280 = -qJD(6) * t447 + t459;
t520 = pkin(9) + pkin(10);
t516 = qJD(2) * pkin(2);
t288 = -t303 * t426 + t430 * t315;
t283 = -pkin(10) * t371 + t288;
t278 = pkin(5) * t385 + t283;
t515 = t278 * t429;
t284 = -pkin(10) * t369 + t289;
t514 = t284 * t429;
t513 = t311 * t390;
t512 = t447 * t390;
t511 = t324 * t426;
t510 = t369 * t385;
t509 = t369 * t390;
t508 = t371 * t385;
t507 = t371 * t390;
t502 = t391 * t426;
t501 = t391 * t430;
t500 = t397 * t426;
t499 = t397 * t430;
t495 = t426 * t383;
t357 = t430 * t366;
t375 = t430 * t383;
t492 = -pkin(10) * t501 + pkin(5) * t392 - t328 * t426 + (-t357 + (pkin(10) * t397 - t349) * t426) * qJD(5) + t526;
t463 = t397 * t479;
t440 = t463 + t502;
t491 = pkin(10) * t440 + t525;
t490 = pkin(5) * t440 + t488;
t347 = pkin(4) * t390 + pkin(9) * t388;
t489 = t426 * t347 + t430 * t524;
t487 = t426 * t349 + t357;
t478 = qJD(6) * t425;
t477 = qJD(6) * t429;
t474 = t429 * t324 - t425 * t325 - t369 * t477;
t473 = qJD(5) * t520;
t464 = t397 * t480;
t462 = t484 * t395;
t281 = t284 * t478;
t460 = t425 * t272 - t281;
t458 = t385 * t430;
t456 = qJD(6) * t278 + t273;
t295 = t397 * t395 + t536;
t455 = t532 * pkin(5) - t307;
t454 = -t485 * t381 - t399 * t383;
t340 = t430 * t347;
t407 = t520 * t430;
t453 = pkin(5) * t390 + qJD(6) * t407 - t524 * t426 + t340 + (pkin(10) * t388 + t473) * t430;
t406 = t520 * t426;
t452 = pkin(10) * t503 + qJD(6) * t406 + t426 * t473 + t489;
t275 = t278 * t425 + t514;
t345 = t430 * t349;
t296 = -pkin(5) * t441 - pkin(10) * t499 - t366 * t426 + t345;
t298 = -pkin(10) * t500 + t487;
t450 = t296 * t425 + t298 * t429;
t386 = -t421 * t498 + t423 * t424;
t387 = t421 * t424 + t423 * t498;
t336 = t427 * t386 + t519 * t387;
t318 = -t336 * t426 - t430 * t497;
t444 = -t336 * t430 + t426 * t497;
t449 = t318 * t429 + t425 * t444;
t448 = t318 * t425 - t429 * t444;
t445 = -t532 * t385 + t375;
t442 = t519 * t386 - t427 * t387;
t439 = -t464 + t501;
t438 = -pkin(9) * t383 + t385 * t302;
t279 = -t371 * t478 + t474;
t432 = qJD(2) ^ 2;
t417 = -pkin(5) * t430 - pkin(4);
t398 = t451 - t516;
t355 = t383 * t441;
t343 = t399 * t397;
t330 = pkin(5) * t500 - t523;
t305 = qJD(2) * t435 + qJD(4) * t336;
t304 = qJD(2) * t434 + qJD(4) * t442;
t292 = t391 * t496 - t425 * t464 - t478 * t500 + (t522 * t499 + t502) * t429;
t291 = -t522 * t342 - t399 * t391;
t287 = qJD(5) * t444 - t304 * t426 + t430 * t465;
t286 = qJD(5) * t318 + t304 * t430 + t426 * t465;
t282 = pkin(5) * t325 + t295;
t274 = -t284 * t425 + t515;
t1 = [(t287 * t385 + t305 * t369 + t318 * t383 - t325 * t442) * MDP(21) + (-t286 * t385 + t305 * t371 - t324 * t442 + t383 * t444) * MDP(22) + ((-qJD(6) * t448 - t286 * t425 + t287 * t429) * t381 + t449 * t383 + t305 * t311 - t442 * t280) * MDP(28) + (-(qJD(6) * t449 + t286 * t429 + t287 * t425) * t381 - t448 * t383 - t305 * t447 - t442 * t279) * MDP(29) + (-t386 * t421 + t387 * t423) * MDP(8) * t395 + (-MDP(14) * t305 - MDP(15) * t304) * qJD(4) + ((-MDP(14) * t383 - MDP(15) * t382) * t431 + (-MDP(8) * t528 + (t388 * MDP(14) + t390 * MDP(15) + (t398 - t467) * MDP(8)) * t428) * qJD(2) + ((-MDP(4) + t527) * t431 + (-MDP(5) * t423 + MDP(6) * t421 - MDP(3)) * t428) * t432) * t422; (t451 * qJD(2) * t484 + t462) * MDP(7) + (-t446 * qJD(3) + qJ(3) * t462 + (t528 + (-t398 - t516) * t428) * t482) * MDP(8) + (t382 * t397 + t390 * t391) * MDP(9) + (t382 * t441 - t383 * t397 - t388 * t391 - t390 * t392) * MDP(10) + (t383 * t414 + t384 * t392 + (-qJD(2) * t441 - t388) * t468) * MDP(14) + (t382 * t414 + t384 * t391) * MDP(15) + (t324 * t499 + t371 * t439) * MDP(16) + ((-t369 * t430 - t371 * t426) * t391 + (-t511 - t325 * t430 + (t369 * t426 - t371 * t430) * qJD(5)) * t397) * MDP(17) + (-t324 * t441 + t371 * t392 + t397 * t375 + t385 * t439) * MDP(18) + (t325 * t441 - t369 * t392 - t385 * t440 - t397 * t495) * MDP(19) + (t385 * t392 - t355) * MDP(20) + (t345 * t383 - (-t303 * t479 + t317) * t441 + t288 * t392 - t523 * t325 + t302 * t463 + (-t366 * t479 + t526) * t385 + t488 * t369 + ((-qJD(5) * t349 - t328) * t385 - t366 * t383 - (-qJD(5) * t315 - t294) * t441 + t295 * t397 + t302 * t391) * t426) * MDP(21) + (-t289 * t392 + t295 * t499 + t439 * t302 - t324 * t523 + t488 * t371 - t487 * t383 + t525 * t385 + t437 * t441) * MDP(22) + (-t279 * t343 - t291 * t447) * MDP(23) + (-t279 * t342 + t280 * t343 - t291 * t311 + t292 * t447) * MDP(24) + (-t279 * t441 + t291 * t381 - t343 * t383 - t392 * t447) * MDP(25) + (t280 * t441 - t292 * t381 - t311 * t392 - t342 * t383) * MDP(26) + (t381 * t392 - t355) * MDP(27) + ((t296 * t429 - t298 * t425) * t383 - t461 * t441 + t274 * t392 + t330 * t280 + t282 * t342 + t297 * t292 + (t491 * t425 + t492 * t429) * t381 + t490 * t311 + (t275 * t441 - t381 * t450) * qJD(6)) * MDP(28) + (-t450 * t383 + (t456 * t429 + t460) * t441 - t275 * t392 + t330 * t279 - t282 * t343 + t297 * t291 + ((-qJD(6) * t296 + t491) * t429 + (qJD(6) * t298 - t492) * t425) * t381 - t490 * t447) * MDP(29) + (t391 * MDP(11) - t392 * MDP(12) - t488 * MDP(14) + t538 * MDP(15)) * qJD(4); -t432 * t527 + (qJD(2) * t446 + t408) * MDP(8) + 0.2e1 * t390 * qJD(4) * MDP(14) + (t409 + (-t388 - t466) * qJD(4)) * MDP(15) + (t445 - t509) * MDP(21) + (-t385 ^ 2 * t430 - t495 - t507) * MDP(22) + (t454 - t513) * MDP(28) + (t512 + t521) * MDP(29); -t388 ^ 2 * MDP(10) + (t409 + (t388 - t466) * qJD(4)) * MDP(11) + (-t295 + t536) * MDP(14) + (t384 * t388 - t529) * MDP(15) + (t371 * t458 + t511) * MDP(16) + ((t324 - t510) * t430 + (-t325 - t508) * t426) * MDP(17) + (t385 * t458 + t495 - t507) * MDP(18) + (t445 + t509) * MDP(19) + (-pkin(4) * t325 - t295 * t430 - t307 * t369 + (-pkin(9) * t479 - t340) * t385 + (t385 * t524 + t438) * t426) * MDP(21) + (-pkin(4) * t324 + t295 * t426 - t307 * t371 + (pkin(9) * t480 + t489) * t385 + t438 * t430) * MDP(22) + (t279 * t400 + t447 * t486) * MDP(23) + (-t279 * t399 - t280 * t400 + t486 * t311 + t447 * t485) * MDP(24) + (t512 - t521) * MDP(25) + (t454 + t513) * MDP(26) + ((-t406 * t429 - t407 * t425) * t383 + t417 * t280 + t282 * t399 + (t425 * t452 - t429 * t453) * t381 + t455 * t311 + t485 * t297) * MDP(28) + (-(-t406 * t425 + t407 * t429) * t383 + t417 * t279 + t282 * t400 + (t425 * t453 + t429 * t452) * t381 - t455 * t447 - t486 * t297) * MDP(29) + (MDP(10) * t390 - t384 * MDP(14) - t385 * MDP(20) - t288 * MDP(21) + t289 * MDP(22) - t381 * MDP(27) - t274 * MDP(28) + t275 * MDP(29) + t388 * MDP(9)) * t390; t371 * t369 * MDP(16) + (-t369 ^ 2 + t371 ^ 2) * MDP(17) + (t324 + t510) * MDP(18) + (-t325 + t508) * MDP(19) + t383 * MDP(20) + (t289 * t385 - t302 * t371 + t433) * MDP(21) + (t288 * t385 + t302 * t369 - t437) * MDP(22) + (t279 + t534) * MDP(25) + (-t280 - t533) * MDP(26) + (-(-t283 * t425 - t514) * t381 - t275 * qJD(6) + (-t311 * t371 - t381 * t478 + t383 * t429) * pkin(5) + t531) * MDP(28) + (t535 + t281 + (-t284 * t381 - t272) * t425 + (t283 * t381 - t456) * t429 + (t371 * t447 - t381 * t477 - t383 * t425) * pkin(5)) * MDP(29) + t530; (t474 + t534) * MDP(25) + (-t459 - t533) * MDP(26) + (t275 * t381 + t531) * MDP(28) + (-t429 * t273 + t274 * t381 - t460 + t535) * MDP(29) + (-MDP(25) * t506 + MDP(26) * t447 - MDP(28) * t275 - MDP(29) * t515) * qJD(6) + t530;];
tauc  = t1;
