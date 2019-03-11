% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:29
% EndTime: 2019-03-09 04:33:37
% DurationCPUTime: 4.46s
% Computational Cost: add. (3033->446), mult. (6767->557), div. (0->0), fcn. (3950->6), ass. (0->181)
t511 = MDP(20) - MDP(25);
t388 = cos(qJ(3));
t456 = qJD(1) * t388;
t369 = -qJD(4) + t456;
t386 = sin(qJ(3));
t441 = qJD(1) * qJD(3);
t429 = t386 * t441;
t510 = qJ(5) * t429 - t369 * qJD(5);
t509 = t369 * MDP(15);
t508 = MDP(21) + MDP(24);
t507 = MDP(5) * t386;
t381 = t386 ^ 2;
t506 = MDP(6) * (-t388 ^ 2 + t381);
t371 = sin(pkin(9)) * pkin(1) + pkin(7);
t357 = t371 * qJD(1);
t377 = t386 * qJD(2);
t324 = t388 * t357 + t377;
t316 = qJD(3) * pkin(8) + t324;
t500 = t388 * qJD(2) - t386 * t357;
t317 = t500 * qJD(3);
t372 = -cos(pkin(9)) * pkin(1) - pkin(2);
t338 = -pkin(3) * t388 - pkin(8) * t386 + t372;
t319 = t338 * qJD(1);
t419 = pkin(3) * t386 - pkin(8) * t388;
t355 = t419 * qJD(3);
t337 = qJD(1) * t355;
t385 = sin(qJ(4));
t387 = cos(qJ(4));
t450 = qJD(4) * t387;
t451 = qJD(4) * t385;
t422 = t316 * t450 + t385 * t317 + t319 * t451 - t387 * t337;
t264 = -pkin(4) * t429 + t422;
t283 = t387 * t316 + t385 * t319;
t365 = t369 * qJ(5);
t277 = -t365 + t283;
t505 = t277 * t369 + t264;
t440 = qJD(3) * qJD(4);
t427 = t385 * t440;
t431 = t386 * t450;
t453 = qJD(3) * t388;
t434 = t385 * t453;
t308 = qJD(1) * (t431 + t434) + t427;
t445 = t387 * qJD(3);
t457 = qJD(1) * t386;
t348 = t385 * t457 - t445;
t504 = t308 * qJ(6) + t348 * qJD(6);
t428 = t388 * t441;
t432 = t386 * t451;
t307 = qJD(1) * t432 + (-t428 - t440) * t387;
t479 = t348 * t369;
t503 = -t307 + t479;
t455 = qJD(3) * t385;
t350 = t387 * t457 + t455;
t502 = t350 * t369 - t308;
t501 = -qJD(5) * t385 - t324;
t499 = 0.2e1 * t510;
t282 = -t385 * t316 + t387 * t319;
t443 = qJD(5) - t282;
t497 = MDP(19) + MDP(23);
t491 = qJ(5) * t385;
t495 = pkin(4) + pkin(5);
t496 = -t387 * t495 - t491;
t345 = t350 ^ 2;
t494 = pkin(8) - qJ(6);
t493 = qJ(5) * t308;
t492 = qJ(5) * t348;
t490 = qJ(5) * t387;
t489 = qJ(6) * t386;
t318 = qJD(3) * t377 + t357 * t453;
t267 = t308 * pkin(4) + t307 * qJ(5) - t350 * qJD(5) + t318;
t488 = t267 * t385;
t487 = t267 * t387;
t275 = qJ(6) * t348 + t283;
t269 = t275 - t365;
t486 = t269 * t369;
t484 = t307 * t385;
t420 = qJD(3) * pkin(3) + t500;
t483 = t420 * t385;
t482 = t420 * t387;
t481 = t318 * t385;
t480 = t318 * t387;
t477 = t369 * t387;
t476 = t371 * t369;
t475 = t371 * t385;
t474 = t385 * t388;
t473 = t386 * t387;
t390 = qJD(3) ^ 2;
t472 = t386 * t390;
t471 = t387 * t388;
t470 = t388 * t390;
t361 = t494 * t387;
t354 = t419 * qJD(1);
t424 = t354 * t387 - t385 * t500;
t469 = (-qJ(6) * t471 - t386 * t495) * qJD(1) - t424 - qJD(4) * t361 + qJD(6) * t385;
t464 = t385 * t354 + t387 * t500;
t287 = qJ(5) * t457 + t464;
t447 = qJD(6) * t387;
t468 = qJ(6) * t385 * t456 + t451 * t494 + t287 + t447;
t405 = -t385 * t495 + t490;
t467 = t369 * t405 + t501;
t416 = pkin(4) * t385 - t490;
t466 = t369 * t416 - t501;
t433 = t388 * t445;
t465 = -t308 * t473 - t348 * t433;
t353 = t371 * t471;
t463 = qJD(4) * t353 + t338 * t451;
t462 = t338 * t450 + t385 * t355;
t430 = t381 * t441;
t363 = t387 * t430;
t461 = -t369 * t433 + t363;
t460 = t385 * t338 + t353;
t458 = MDP(20) * t385;
t358 = qJD(1) * t372;
t454 = qJD(3) * t386;
t452 = qJD(4) * t348;
t448 = qJD(5) * t387;
t446 = t386 * MDP(16);
t274 = qJ(6) * t350 + t282;
t444 = qJD(5) - t274;
t397 = qJ(5) * t350 + t420;
t273 = -t348 * t495 + qJD(6) + t397;
t442 = qJD(6) + t273;
t439 = pkin(8) * t369 * t385;
t438 = pkin(8) * t477;
t437 = pkin(8) * t454;
t436 = pkin(8) * t445;
t435 = t387 * t317 + t319 * t450 + t385 * t337;
t426 = -pkin(4) - t475;
t425 = MDP(18) - t508;
t352 = t371 * t474;
t423 = t338 * t387 - t352;
t343 = t369 * t431;
t421 = t385 * t430;
t292 = -qJ(5) * t388 + t460;
t418 = t355 * t387 - t463;
t417 = pkin(4) * t387 + t491;
t399 = t316 * t451 - t435;
t263 = -t399 + t510;
t415 = t263 * t387 + t264 * t385;
t268 = t369 * t495 + t444;
t414 = t268 * t387 - t269 * t385;
t413 = t268 * t385 + t269 * t387;
t276 = pkin(4) * t369 + t443;
t412 = t276 * t387 - t277 * t385;
t411 = t276 * t385 + t277 * t387;
t408 = 0.2e1 * qJD(3) * t358;
t407 = qJ(5) * t454 - qJD(5) * t388 + t462;
t406 = t371 + t416;
t262 = -pkin(5) * t308 - t267;
t404 = -t262 * t385 - t273 * t450;
t403 = t262 * t387 - t273 * t451;
t401 = qJ(6) * t307 + t422;
t398 = t369 * t434 + t343 - t421;
t396 = -t371 + t405;
t394 = -t282 * t369 + t399;
t393 = t307 + t479;
t392 = -t429 * t495 + t401;
t380 = t388 * pkin(4);
t360 = t494 * t385;
t356 = -pkin(3) - t417;
t339 = pkin(3) - t496;
t329 = t348 * t454;
t328 = t350 * t454;
t309 = t406 * t386;
t296 = pkin(4) * t350 + t492;
t294 = t396 * t386;
t293 = t380 - t423;
t290 = t385 * t489 + t292;
t289 = -pkin(4) * t457 - t424;
t286 = -t350 * t495 - t492;
t284 = pkin(5) * t388 + t352 + t380 + (-t338 - t489) * t387;
t281 = pkin(4) * t348 - t397;
t279 = (qJD(4) * t417 - t448) * t386 + t406 * t453;
t272 = t426 * t454 - t418;
t271 = (t496 * qJD(4) + t448) * t386 + t396 * t453;
t270 = (-t386 * t445 - t451 * t388) * t371 + t407;
t266 = (qJ(6) * qJD(4) - qJD(3) * t371) * t473 + (qJD(6) * t386 + (qJ(6) * qJD(3) - qJD(4) * t371) * t388) * t385 + t407;
t265 = (-qJ(6) * t453 - t355) * t387 + (qJ(6) * t451 - t447 + (-pkin(5) + t426) * qJD(3)) * t386 + t463;
t261 = -qJD(6) * t350 + t392;
t260 = t263 + t504;
t1 = [0.2e1 * t428 * t507 - 0.2e1 * t441 * t506 + MDP(7) * t470 - MDP(8) * t472 + (-t371 * t470 + t386 * t408) * MDP(10) + (t371 * t472 + t388 * t408) * MDP(11) + (-t307 * t473 + (-t432 + t433) * t350) * MDP(12) + (-t350 * t431 + (-t350 * t453 + (t307 + t452) * t386) * t385 + t465) * MDP(13) + (t307 * t388 + t369 * t432 + t328 + t461) * MDP(14) + (t343 + t308 * t388 + (-t348 * t386 + (-qJD(1) * t381 + t369 * t388) * t385) * qJD(3)) * MDP(15) + (-t369 - t456) * qJD(3) * t446 + (-t418 * t369 + ((t348 * t371 - t483) * qJD(3) + t422) * t388 + (-t420 * t450 + t371 * t308 + t481 + (qJD(1) * t423 - t369 * t475 + t282) * qJD(3)) * t386) * MDP(17) + (t462 * t369 + ((-t316 - t476) * t451 + (t350 * t371 - t482) * qJD(3) + t435) * t388 + (t420 * t451 - t371 * t307 + t480 + (-qJD(1) * t460 - t387 * t476 - t283) * qJD(3)) * t386) * MDP(18) + (t272 * t369 + t279 * t348 + t308 * t309 + (t281 * t455 + t264) * t388 + (t281 * t450 + t488 + (-qJD(1) * t293 - t276) * qJD(3)) * t386) * MDP(19) + (-t270 * t348 + t272 * t350 - t292 * t308 - t293 * t307 + t412 * t453 + (-qJD(4) * t411 - t263 * t385 + t264 * t387) * t386) * MDP(20) + (-t270 * t369 - t279 * t350 + t307 * t309 + (-t281 * t445 - t263) * t388 + (t281 * t451 - t487 + (qJD(1) * t292 + t277) * qJD(3)) * t386) * MDP(21) + (t263 * t292 + t264 * t293 + t267 * t309 + t270 * t277 + t272 * t276 + t279 * t281) * MDP(22) + (t265 * t369 - t271 * t348 - t294 * t308 + (-t273 * t455 + t261) * t388 + ((-qJD(1) * t284 - t268) * qJD(3) + t404) * t386) * MDP(23) + (-t266 * t369 + t271 * t350 - t294 * t307 + (t273 * t445 - t260) * t388 + ((qJD(1) * t290 + t269) * qJD(3) + t403) * t386) * MDP(24) + (-t265 * t350 + t266 * t348 + t284 * t307 + t290 * t308 - t414 * t453 + (qJD(4) * t413 + t260 * t385 - t261 * t387) * t386) * MDP(25) + (t260 * t290 + t261 * t284 + t262 * t294 + t265 * t268 + t266 * t269 + t271 * t273) * MDP(26); (t329 - t421) * MDP(17) + (t328 - t363) * MDP(18) + (t329 + t398) * MDP(19) + t465 * MDP(20) + t461 * MDP(21) + t398 * MDP(23) + (-t328 + t461) * MDP(24) + (-t390 * MDP(11) - t267 * MDP(22) + t262 * MDP(26) + (-MDP(17) - t497) * t308 + t425 * t307 + (t350 * t458 + t411 * MDP(22) + (t348 * t387 - t350 * t385) * MDP(25) + t413 * MDP(26) + (t385 * MDP(17) + t387 * MDP(18)) * t369) * qJD(3)) * t388 + (-t390 * MDP(10) - t307 * t458 + t415 * MDP(22) + (t308 * t387 + t484) * MDP(25) + (t260 * t387 + t261 * t385) * MDP(26) + (-t350 * MDP(21) + t281 * MDP(22) + t348 * MDP(23) - t273 * MDP(26)) * qJD(3) + (t412 * MDP(22) + t414 * MDP(26) + (t387 * MDP(17) - t385 * t425) * t369 + t511 * (t348 * t385 + t350 * t387)) * qJD(4)) * t386; (qJD(3) * t324 - t358 * t457 - t318) * MDP(10) - t358 * t456 * MDP(11) + (-t350 * t477 - t484) * MDP(12) + (t502 * t385 + t503 * t387) * MDP(13) - t369 * t450 * MDP(14) + t451 * t509 + (-pkin(3) * t308 - t480 + t424 * t369 - t324 * t348 + (t438 - t483) * qJD(4)) * MDP(17) + (pkin(3) * t307 + t481 - t464 * t369 - t324 * t350 + (-t439 - t482) * qJD(4)) * MDP(18) + (-t487 - t289 * t369 + t308 * t356 - t466 * t348 + (t281 * t385 + t438) * qJD(4)) * MDP(19) + (t287 * t348 - t289 * t350 + (t263 - t369 * t276 + (qJD(4) * t350 - t308) * pkin(8)) * t387 + ((-t307 + t452) * pkin(8) + t505) * t385) * MDP(20) + (-t488 + t287 * t369 + t307 * t356 + t466 * t350 + (-t281 * t387 + t439) * qJD(4)) * MDP(21) + (t267 * t356 - t276 * t289 - t277 * t287 - t466 * t281 + (qJD(4) * t412 + t415) * pkin(8)) * MDP(22) + (-t308 * t339 + t467 * t348 - t469 * t369 + t403) * MDP(23) + (-t307 * t339 - t467 * t350 + t468 * t369 - t404) * MDP(24) + (t307 * t360 + t308 * t361 + t469 * t350 - t468 * t348 + (t268 * t369 - t260) * t387 + (-t261 - t486) * t385) * MDP(25) + (t260 * t361 + t261 * t360 + t262 * t339 - t268 * t469 - t269 * t468 - t467 * t273) * MDP(26) + ((t369 * t471 + (-t350 + t455) * t386) * MDP(14) + (-t369 * t474 + (t348 + t445) * t386) * MDP(15) + t369 * t446 + (-t282 * t386 + (t388 * t420 - t437) * t385) * MDP(17) + (t420 * t471 + (t283 - t436) * t386) * MDP(18) + (t276 * t386 + (-t281 * t388 - t437) * t385) * MDP(19) + (t281 * t471 + (-t277 + t436) * t386) * MDP(21) + (t273 * t474 + (-qJD(3) * t360 + t268) * t386) * MDP(23) + (-t273 * t471 + (qJD(3) * t361 - t269) * t386) * MDP(24) + (-t388 * t507 + t506) * qJD(1)) * qJD(1); -t393 * MDP(14) - MDP(15) * t427 + t394 * MDP(18) + (pkin(4) * t307 - t493) * MDP(20) + (-t394 + t499) * MDP(21) + (-pkin(4) * t264 + qJ(5) * t263 - t276 * t283 + t277 * t443 - t281 * t296) * MDP(22) + (-t275 * t369 - t401) * MDP(23) + (t274 * t369 - t399 + t499 + t504) * MDP(24) + (-t307 * t495 + t493) * MDP(25) + (qJ(5) * t260 - t261 * t495 - t268 * t275 + t269 * t444 - t273 * t286) * MDP(26) + (-t509 + t420 * MDP(17) - t281 * MDP(19) + (t277 - t283) * MDP(20) + t296 * MDP(21) + t442 * MDP(23) - t286 * MDP(24) + (-t269 + t275) * MDP(25) + MDP(13) * t350) * t350 + (-MDP(15) * t431 + (-MDP(15) * t474 + (0.2e1 * pkin(4) * MDP(19) + 0.2e1 * t495 * MDP(23) + MDP(16)) * t386) * qJD(3)) * qJD(1) + (t350 * MDP(12) - t420 * MDP(18) - t296 * MDP(19) + (t276 - t443) * MDP(20) - t281 * MDP(21) + t286 * MDP(23) + t273 * MDP(24) + (-t268 + t444) * MDP(25) - MDP(13) * t348) * t348 + (MDP(17) + MDP(19)) * (-t283 * t369 - t422); (t281 * t350 + t505) * MDP(22) + (-t442 * t350 + t392 + t486) * MDP(26) + t497 * (t348 * t350 - t429) + t508 * (-t369 ^ 2 - t345) - t511 * t393; t502 * MDP(23) + t503 * MDP(24) + (-t348 ^ 2 - t345) * MDP(25) + (t268 * t350 - t269 * t348 + t262) * MDP(26);];
tauc  = t1;
