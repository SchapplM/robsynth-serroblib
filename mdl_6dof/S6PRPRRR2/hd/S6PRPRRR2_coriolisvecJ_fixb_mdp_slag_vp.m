% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRR2
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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:30:03
% EndTime: 2019-03-08 20:30:11
% DurationCPUTime: 4.21s
% Computational Cost: add. (2469->386), mult. (6199->548), div. (0->0), fcn. (4788->12), ass. (0->186)
t385 = cos(pkin(12));
t390 = sin(qJ(2));
t384 = sin(pkin(6));
t457 = qJD(1) * t384;
t439 = t390 * t457;
t364 = t385 * t439;
t383 = sin(pkin(12));
t394 = cos(qJ(2));
t438 = t394 * t457;
t326 = t383 * t438 + t364;
t389 = sin(qJ(4));
t393 = cos(qJ(4));
t415 = pkin(4) * t389 - pkin(9) * t393;
t360 = t415 * qJD(4);
t504 = t326 - t360;
t362 = qJD(2) * pkin(2) + t438;
t320 = t383 * t362 + t364;
t313 = qJD(2) * pkin(8) + t320;
t386 = cos(pkin(6));
t373 = qJD(1) * t386 + qJD(3);
t494 = -t389 * t313 + t373 * t393;
t287 = -qJD(4) * pkin(4) - t494;
t388 = sin(qJ(5));
t392 = cos(qJ(5));
t445 = t392 * qJD(4);
t456 = qJD(2) * t389;
t351 = t388 * t456 - t445;
t276 = pkin(5) * t351 + t287;
t391 = cos(qJ(6));
t454 = qJD(4) * t388;
t353 = t392 * t456 + t454;
t387 = sin(qJ(6));
t477 = t353 * t387;
t302 = t391 * t351 + t477;
t503 = t276 * t302;
t455 = qJD(2) * t393;
t374 = -qJD(5) + t455;
t371 = -qJD(6) + t374;
t502 = t302 * t371;
t408 = t351 * t387 - t391 * t353;
t501 = t371 * t408;
t458 = MDP(12) * t393;
t444 = qJD(2) * qJD(4);
t431 = t393 * t444;
t451 = qJD(5) * t388;
t435 = t389 * t451;
t443 = qJD(4) * qJD(5);
t323 = -qJD(2) * t435 + (t431 + t443) * t392;
t291 = t393 * t313 + t389 * t373;
t288 = qJD(4) * pkin(9) + t291;
t363 = t383 * t439;
t319 = t362 * t385 - t363;
t405 = -pkin(4) * t393 - pkin(9) * t389 - pkin(3);
t296 = qJD(2) * t405 - t319;
t473 = t388 * t296;
t265 = t288 * t392 + t473;
t332 = (t383 * t390 - t385 * t394) * t384;
t328 = qJD(2) * t332;
t322 = qJD(1) * t328;
t272 = qJD(4) * t494 - t322 * t393;
t333 = (t383 * t394 + t385 * t390) * t384;
t297 = (qJD(1) * t333 + t360) * qJD(2);
t424 = t388 * t272 - t392 * t297;
t397 = -qJD(5) * t265 - t424;
t432 = t389 * t444;
t253 = pkin(5) * t432 - pkin(10) * t323 + t397;
t450 = qJD(5) * t392;
t434 = t389 * t450;
t452 = qJD(4) * t393;
t437 = t388 * t452;
t398 = t434 + t437;
t324 = qJD(2) * t398 + t388 * t443;
t442 = t392 * t272 + t296 * t450 + t388 * t297;
t400 = -t288 * t451 + t442;
t254 = -pkin(10) * t324 + t400;
t427 = t391 * t253 - t387 * t254;
t500 = t276 * t408 + t427;
t429 = MDP(24) * t456;
t499 = qJD(4) * t429 + (-t302 ^ 2 + t408 ^ 2) * MDP(21) - t302 * t408 * MDP(20);
t498 = MDP(6) * t389;
t381 = t389 ^ 2;
t497 = MDP(7) * (-t393 ^ 2 + t381);
t355 = t387 * t392 + t388 * t391;
t335 = t355 * t389;
t329 = t385 * t438 - t363;
t453 = qJD(4) * t389;
t471 = t388 * t393;
t375 = pkin(2) * t383 + pkin(8);
t474 = t375 * t388;
t496 = -t329 * t471 + t504 * t392 - t453 * t474;
t490 = pkin(2) * t385;
t349 = t405 - t490;
t469 = t392 * t393;
t495 = -t329 * t469 + t349 * t450 - t504 * t388;
t493 = t393 * t445 - t435;
t492 = qJD(5) + qJD(6);
t421 = t323 * t387 + t391 * t324;
t267 = -qJD(6) * t408 + t421;
t491 = pkin(9) + pkin(10);
t264 = -t288 * t388 + t392 * t296;
t261 = -pkin(10) * t353 + t264;
t259 = -pkin(5) * t374 + t261;
t489 = t259 * t391;
t262 = -pkin(10) * t351 + t265;
t488 = t262 * t391;
t487 = t267 * t393;
t273 = t313 * t452 - t389 * t322 + t373 * t453;
t486 = t273 * t388;
t485 = t273 * t392;
t354 = t387 * t388 - t391 * t392;
t401 = t354 * t393;
t277 = -qJD(4) * t401 - t335 * t492;
t484 = t277 * t371;
t483 = t287 * t388;
t482 = t287 * t392;
t481 = t323 * t388;
t480 = t324 * t393;
t479 = t351 * t374;
t478 = t353 * t374;
t475 = t374 * t392;
t472 = t388 * t389;
t470 = t389 * t392;
t356 = t375 * t469;
t404 = pkin(5) * t389 - pkin(10) * t469;
t468 = -t404 * qJD(4) - (-t356 + (pkin(10) * t389 - t349) * t388) * qJD(5) + t496;
t467 = (-t389 * t445 - t393 * t451) * t375 - t398 * pkin(10) + t495;
t449 = qJD(6) * t387;
t278 = -t449 * t472 + (t470 * t492 + t437) * t391 + t493 * t387;
t466 = t278 * t371 - t335 * t432;
t357 = t415 * qJD(2);
t465 = t388 * t357 + t392 * t494;
t464 = qJD(2) * t401 - t354 * t492;
t463 = (-t455 + t492) * t355;
t460 = t388 * t349 + t356;
t448 = qJD(6) * t391;
t447 = t287 * qJD(5);
t441 = t391 * t323 - t387 * t324 - t351 * t448;
t440 = qJD(5) * t491;
t433 = t388 * t455;
t428 = MDP(17) * t453;
t260 = t262 * t449;
t426 = t387 * t253 - t260;
t266 = -t353 * t449 + t441;
t425 = -t266 * t393 - t408 * t453;
t423 = t374 * t375 + t288;
t422 = t392 * t357 - t388 * t494;
t420 = -t323 * t393 + t353 * t453;
t419 = qJD(6) * t259 + t254;
t418 = t374 * t435;
t417 = t374 * t434;
t416 = -t291 + (-t433 + t451) * pkin(5);
t369 = t491 * t392;
t414 = qJD(2) * t404 + qJD(6) * t369 + t392 * t440 + t422;
t368 = t491 * t388;
t413 = pkin(10) * t433 - qJD(6) * t368 - t388 * t440 - t465;
t256 = t259 * t387 + t488;
t311 = t333 * t393 + t386 * t389;
t281 = -t311 * t388 + t332 * t392;
t282 = t311 * t392 + t332 * t388;
t412 = t281 * t391 - t282 * t387;
t411 = t281 * t387 + t282 * t391;
t338 = t392 * t349;
t295 = -pkin(10) * t470 + t338 + (-pkin(5) - t474) * t393;
t299 = -pkin(10) * t472 + t460;
t410 = t295 * t387 + t299 * t391;
t310 = t333 * t389 - t386 * t393;
t406 = qJD(2) * t381 - t374 * t393;
t327 = qJD(2) * t333;
t321 = qJD(1) * t327;
t395 = qJD(4) ^ 2;
t403 = qJD(2) * t326 - t375 * t395 - t321;
t312 = -qJD(2) * pkin(3) - t319;
t402 = qJD(4) * (qJD(2) * (-pkin(3) - t490) + t312 + t329);
t399 = t406 * t388;
t396 = qJD(2) ^ 2;
t379 = -pkin(5) * t392 - pkin(4);
t341 = (pkin(5) * t388 + t375) * t389;
t336 = t354 * t389;
t315 = pkin(5) * t398 + t375 * t452;
t280 = -qJD(4) * t310 - t328 * t393;
t279 = qJD(4) * t311 - t328 * t389;
t263 = pkin(5) * t324 + t273;
t258 = qJD(5) * t281 + t280 * t392 + t327 * t388;
t257 = -qJD(5) * t282 - t280 * t388 + t327 * t392;
t255 = -t262 * t387 + t489;
t1 = [(-t319 * t327 - t320 * t328 + t321 * t332 - t322 * t333) * MDP(5) + (-t257 * t374 + t279 * t351 + t310 * t324) * MDP(18) + (t258 * t374 + t279 * t353 + t310 * t323) * MDP(19) + (-(-qJD(6) * t411 + t257 * t391 - t258 * t387) * t371 + t279 * t302 + t310 * t267) * MDP(25) + ((qJD(6) * t412 + t257 * t387 + t258 * t391) * t371 - t279 * t408 + t310 * t266) * MDP(26) + (-MDP(3) * t390 - MDP(4) * t394) * t396 * t384 + (-MDP(11) * t279 - MDP(12) * t280) * qJD(4) + ((-MDP(11) * t393 + MDP(12) * t389) * t327 + (t332 * t458 + (t332 * MDP(11) + t281 * MDP(18) - t282 * MDP(19) + MDP(25) * t412 - MDP(26) * t411) * t389) * qJD(4)) * qJD(2); (t319 * t326 - t320 * t329 + (-t321 * t385 - t322 * t383) * pkin(2)) * MDP(5) + 0.2e1 * t431 * t498 - 0.2e1 * t444 * t497 + (t389 * t402 + t393 * t403) * MDP(11) + (-t389 * t403 + t393 * t402) * MDP(12) + (t323 * t470 + t353 * t493) * MDP(13) + ((-t351 * t392 - t353 * t388) * t452 + (-t481 - t324 * t392 + (t351 * t388 - t353 * t392) * qJD(5)) * t389) * MDP(14) + (t406 * t445 + t418 + t420) * MDP(15) + (t417 + t480 + (-t351 * t389 - t399) * qJD(4)) * MDP(16) + (-t374 - t455) * t428 + ((t349 * t451 + t496) * t374 + ((t351 * t375 + t483) * qJD(4) + (t392 * t423 + t473) * qJD(5) + t424) * t393 + (t392 * t447 + t486 + t375 * t324 - t329 * t351 + ((-t375 * t471 + t338) * qJD(2) + t264) * qJD(4)) * t389) * MDP(18) + (t495 * t374 + (-t423 * t451 + (t353 * t375 + t482) * qJD(4) + t442) * t393 + (-t388 * t447 + t485 + t375 * t323 - t329 * t353 + (-qJD(2) * t460 - t375 * t475 - t265) * qJD(4)) * t389) * MDP(19) + (-t266 * t336 - t277 * t408) * MDP(20) + (-t266 * t335 + t267 * t336 - t277 * t302 + t278 * t408) * MDP(21) + (-t336 * t432 + t425 - t484) * MDP(22) + (-t302 * t453 + t466 + t487) * MDP(23) + (-t371 - t455) * MDP(24) * t453 + (-t427 * t393 + t315 * t302 + t341 * t267 + t263 * t335 + t276 * t278 + (t387 * t467 + t391 * t468) * t371 + (t256 * t393 + t371 * t410) * qJD(6) + (-t329 * t302 + ((t295 * t391 - t299 * t387) * qJD(2) + t255) * qJD(4)) * t389) * MDP(25) + ((t419 * t391 + t426) * t393 - t315 * t408 + t341 * t266 - t263 * t336 + t276 * t277 + ((qJD(6) * t295 + t467) * t391 + (-qJD(6) * t299 - t468) * t387) * t371 + (t329 * t408 + (-qJD(2) * t410 - t256) * qJD(4)) * t389) * MDP(26) + (MDP(8) * t393 - MDP(9) * t389) * t395; (t417 - t480) * MDP(18) + (-t418 + t420) * MDP(19) + (t466 - t487) * MDP(25) + (t425 + t484) * MDP(26) + (-MDP(11) * t389 - t458) * t395 + (-t406 * MDP(19) * t392 - MDP(18) * t399 + (MDP(26) * qJD(2) * t336 + t351 * MDP(18) + t302 * MDP(25)) * t389) * qJD(4); (qJD(4) * t291 - t312 * t456 - t273) * MDP(11) + (-qJD(2) * t312 + t322) * t458 + (-t353 * t475 + t481) * MDP(13) + ((t323 + t479) * t392 + (-t324 + t478) * t388) * MDP(14) + (-t374 * t450 + (t374 * t469 + (-t353 + t454) * t389) * qJD(2)) * MDP(15) + (t374 * t451 + (-t374 * t471 + (t351 + t445) * t389) * qJD(2)) * MDP(16) + t374 * MDP(17) * t456 + (-pkin(4) * t324 - t485 + t422 * t374 - t291 * t351 + (pkin(9) * t475 + t483) * qJD(5) + (-t264 * t389 + (-pkin(9) * t453 - t287 * t393) * t388) * qJD(2)) * MDP(18) + (-pkin(4) * t323 + t486 - t465 * t374 - t291 * t353 + (-pkin(9) * t374 * t388 + t482) * qJD(5) + (-t287 * t469 + (-pkin(9) * t445 + t265) * t389) * qJD(2)) * MDP(19) + (t266 * t355 - t408 * t464) * MDP(20) + (-t266 * t354 - t267 * t355 - t302 * t464 + t408 * t463) * MDP(21) + (-t464 * t371 + (qJD(4) * t355 + t408) * t456) * MDP(22) + (t463 * t371 + (-qJD(4) * t354 + t302) * t456) * MDP(23) + t371 * t429 + (t263 * t354 + t379 * t267 + (t387 * t413 + t391 * t414) * t371 + t416 * t302 + t463 * t276 + ((-t368 * t391 - t369 * t387) * qJD(4) - t255) * t456) * MDP(25) + (t263 * t355 + t379 * t266 + (-t387 * t414 + t391 * t413) * t371 - t416 * t408 + t464 * t276 + (-(-t368 * t387 + t369 * t391) * qJD(4) + t256) * t456) * MDP(26) + (-t393 * t498 + t497) * t396; t353 * t351 * MDP(13) + (-t351 ^ 2 + t353 ^ 2) * MDP(14) + (t323 - t479) * MDP(15) + (-t324 - t478) * MDP(16) + qJD(2) * t428 + (-t265 * t374 - t287 * t353 + t397) * MDP(18) + (-t264 * t374 + t287 * t351 - t400) * MDP(19) + (t266 - t502) * MDP(22) + (-t267 + t501) * MDP(23) + ((-t261 * t387 - t488) * t371 - t256 * qJD(6) + (-t302 * t353 + t371 * t449 + t391 * t432) * pkin(5) + t500) * MDP(25) + (t503 + t260 + (t262 * t371 - t253) * t387 + (-t261 * t371 - t419) * t391 + (t353 * t408 + t371 * t448 - t387 * t432) * pkin(5)) * MDP(26) + t499; (t441 - t502) * MDP(22) + (-t421 + t501) * MDP(23) + (-t256 * t371 + t500) * MDP(25) + (-t391 * t254 - t255 * t371 - t426 + t503) * MDP(26) + (-MDP(22) * t477 + MDP(23) * t408 - MDP(25) * t256 - MDP(26) * t489) * qJD(6) + t499;];
tauc  = t1;
