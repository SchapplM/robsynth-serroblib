% Calculate Coriolis joint torque vector for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:18
% EndTime: 2021-01-16 02:06:32
% DurationCPUTime: 6.10s
% Computational Cost: add. (3773->419), mult. (10046->604), div. (0->0), fcn. (7853->12), ass. (0->185)
t412 = sin(qJ(2));
t407 = sin(pkin(6));
t463 = qJD(1) * t407;
t451 = t412 * t463;
t411 = sin(qJ(3));
t459 = qJD(3) * t411;
t504 = pkin(3) * t459 - t451;
t414 = cos(qJ(3));
t493 = cos(pkin(11));
t443 = t493 * t414;
t394 = qJD(2) * t443;
t406 = sin(pkin(11));
t461 = qJD(2) * t411;
t371 = t406 * t461 - t394;
t367 = qJD(6) + t371;
t384 = t406 * t414 + t493 * t411;
t373 = t384 * qJD(3);
t424 = -t406 * t411 + t443;
t376 = t424 * qJD(3);
t503 = pkin(4) * t373 - qJ(5) * t376 - qJD(5) * t384 + t504;
t495 = qJ(4) + pkin(8);
t445 = qJD(3) * t495;
t368 = qJD(4) * t414 - t411 * t445;
t423 = -qJD(4) * t411 - t414 * t445;
t415 = cos(qJ(2));
t450 = t415 * t463;
t468 = t493 * t368 + t406 * t423 - t424 * t450;
t374 = t384 * qJD(2);
t405 = sin(pkin(12));
t408 = cos(pkin(12));
t359 = t408 * qJD(3) - t374 * t405;
t413 = cos(qJ(6));
t502 = t413 * t359;
t410 = sin(qJ(6));
t385 = t405 * t413 + t408 * t410;
t465 = t367 * t385;
t358 = qJD(3) * t405 + t374 * t408;
t308 = t358 * t413 + t359 * t410;
t501 = MDP(5) * t414;
t500 = (t411 ^ 2 - t414 ^ 2) * MDP(6);
t472 = -t468 * t405 + t503 * t408;
t471 = t503 * t405 + t468 * t408;
t469 = t368 * t406 - t384 * t450 - t493 * t423;
t364 = qJD(2) * t373;
t476 = t413 * t408;
t481 = t405 * t410;
t383 = -t476 + t481;
t466 = t367 * t383;
t499 = -t364 * t385 + t466 * t367;
t370 = t371 ^ 2;
t498 = pkin(3) * t411;
t497 = pkin(9) * t408;
t397 = pkin(3) * t406 + qJ(5);
t496 = pkin(9) + t397;
t494 = qJD(2) * pkin(2);
t387 = qJD(2) * pkin(8) + t451;
t409 = cos(pkin(6));
t462 = qJD(1) * t409;
t395 = t414 * t462;
t436 = qJD(4) + t450;
t318 = (-t387 * t411 + t395) * qJD(3) + (-qJ(4) * t459 + t436 * t414) * qJD(2);
t449 = t411 * t462;
t319 = (-t387 * t414 - t449) * qJD(3) + (-qJ(4) * qJD(3) * t414 - t436 * t411) * qJD(2);
t276 = t318 * t406 - t493 * t319;
t480 = t407 * t412;
t379 = t409 * t411 + t414 * t480;
t428 = t409 * t414 - t411 * t480;
t333 = t379 * t406 - t493 * t428;
t492 = t276 * t333;
t390 = t495 * t411;
t391 = t495 * t414;
t354 = t493 * t390 + t391 * t406;
t491 = t276 * t354;
t306 = t358 * t410 - t502;
t490 = t306 * t374;
t489 = t308 * t374;
t455 = qJD(2) * qJD(3);
t446 = t411 * t455;
t392 = t406 * t446;
t365 = qJD(3) * t394 - t392;
t487 = t365 * t405;
t486 = t365 * t408;
t485 = t371 * t405;
t484 = t376 * t405;
t483 = t384 * t405;
t482 = t384 * t408;
t442 = qJ(4) * qJD(2) + t387;
t350 = t442 * t414 + t449;
t341 = t406 * t350;
t479 = t407 * t415;
t417 = qJD(2) ^ 2;
t478 = t407 * t417;
t416 = qJD(3) ^ 2;
t477 = t411 * t416;
t475 = t414 * t416;
t474 = pkin(5) * t373 - t376 * t497 + t472;
t473 = pkin(9) * t484 - t471;
t277 = t493 * t318 + t406 * t319;
t273 = qJD(3) * qJD(5) + t277;
t460 = qJD(2) * t412;
t448 = t407 * t460;
t369 = pkin(3) * t446 + qJD(1) * t448;
t290 = pkin(4) * t364 - qJ(5) * t365 - qJD(5) * t374 + t369;
t260 = t408 * t273 + t405 * t290;
t349 = -t442 * t411 + t395;
t345 = qJD(3) * pkin(3) + t349;
t444 = t493 * t350;
t298 = t406 * t345 + t444;
t293 = qJD(3) * qJ(5) + t298;
t401 = -pkin(3) * t414 - pkin(2);
t366 = t401 * qJD(2) + qJD(4) - t450;
t315 = pkin(4) * t371 - qJ(5) * t374 + t366;
t267 = t408 * t293 + t405 * t315;
t302 = t493 * t349 - t341;
t454 = pkin(3) * t461;
t328 = pkin(4) * t374 + qJ(5) * t371 + t454;
t272 = t408 * t302 + t405 * t328;
t470 = pkin(5) * t484 + t469;
t337 = -pkin(4) * t424 - qJ(5) * t384 + t401;
t355 = -t406 * t390 + t493 * t391;
t295 = t405 * t337 + t408 * t355;
t456 = qJD(6) * t413;
t467 = t359 * t456 + t365 * t476;
t262 = pkin(9) * t359 + t267;
t458 = qJD(6) * t262;
t457 = qJD(6) * t384;
t452 = t412 * t478;
t447 = qJD(2) * t479;
t259 = -t273 * t405 + t408 * t290;
t266 = -t293 * t405 + t408 * t315;
t271 = -t302 * t405 + t408 * t328;
t294 = t408 * t337 - t355 * t405;
t300 = t349 * t406 + t444;
t258 = -pkin(9) * t487 + t260;
t261 = pkin(5) * t371 - pkin(9) * t358 + t266;
t441 = -qJD(6) * t261 - t258;
t440 = t414 * t447;
t439 = t411 * t447;
t400 = -t493 * pkin(3) - pkin(4);
t438 = -t383 * t364 - t465 * t367;
t297 = t493 * t345 - t341;
t254 = t261 * t413 - t262 * t410;
t255 = t261 * t410 + t262 * t413;
t435 = -t266 * t405 + t267 * t408;
t434 = t276 * t384 + t354 * t365;
t280 = -pkin(5) * t424 - pkin(9) * t482 + t294;
t282 = -pkin(9) * t483 + t295;
t433 = t280 * t413 - t282 * t410;
t432 = t280 * t410 + t282 * t413;
t334 = t493 * t379 + t406 * t428;
t316 = -t334 * t405 - t408 * t479;
t317 = t334 * t408 - t405 * t479;
t431 = t316 * t413 - t317 * t410;
t430 = t316 * t410 + t317 * t413;
t429 = -qJD(6) * t358 - t487;
t427 = qJD(2) * t494;
t381 = t496 * t408;
t426 = pkin(5) * t374 + qJD(5) * t405 + qJD(6) * t381 + t371 * t497 + t271;
t380 = t496 * t405;
t425 = pkin(9) * t485 - qJD(5) * t408 + qJD(6) * t380 + t272;
t292 = -qJD(3) * pkin(4) + qJD(5) - t297;
t422 = t292 * t376 + t434;
t421 = t379 * qJD(3);
t420 = -0.2e1 * qJD(3) * t494;
t419 = -t364 * t397 + t365 * t400 + (-qJD(5) + t292) * t371;
t418 = -t421 - t439;
t275 = t308 * qJD(6) + t385 * t365;
t389 = -t408 * pkin(5) + t400;
t348 = t428 * qJD(3) + t440;
t332 = t383 * t384;
t331 = t385 * t384;
t325 = pkin(5) * t483 + t354;
t301 = t493 * t348 + t406 * t418;
t299 = t348 * t406 - t493 * t418;
t287 = t301 * t408 + t405 * t448;
t286 = -t301 * t405 + t408 * t448;
t285 = -pkin(5) * t485 + t300;
t284 = t385 * t376 + t456 * t482 - t457 * t481;
t283 = -t383 * t376 - t385 * t457;
t281 = -pkin(5) * t359 + t292;
t274 = t429 * t410 + t467;
t268 = pkin(5) * t487 + t276;
t257 = pkin(5) * t364 - pkin(9) * t486 + t259;
t256 = t413 * t257;
t1 = [-MDP(3) * t452 - t415 * MDP(4) * t478 + (-t414 * t452 + (-t421 - 0.2e1 * t439) * qJD(3)) * MDP(10) + (t411 * t452 + (-t348 - t440) * qJD(3)) * MDP(11) + (-qJD(3) * t299 + (-t364 * t415 + t371 * t460) * t407) * MDP(12) + (-qJD(3) * t301 + (-t365 * t415 + t374 * t460) * t407) * MDP(13) + (t299 * t374 - t301 * t371 + t333 * t365 - t334 * t364) * MDP(14) + (t492 + t277 * t334 - t297 * t299 + t298 * t301 + (t366 * t460 - t369 * t415) * t407) * MDP(15) + (t286 * t371 - t299 * t359 + t316 * t364 + t333 * t487) * MDP(16) + (-t287 * t371 + t299 * t358 - t317 * t364 + t333 * t486) * MDP(17) + (-t286 * t358 + t287 * t359 + (-t316 * t408 - t317 * t405) * t365) * MDP(18) + (t259 * t316 + t260 * t317 + t266 * t286 + t267 * t287 + t292 * t299 + t492) * MDP(19) + ((-t430 * qJD(6) + t286 * t413 - t287 * t410) * t367 + t431 * t364 + t299 * t306 + t333 * t275) * MDP(25) + (-(t431 * qJD(6) + t286 * t410 + t287 * t413) * t367 - t430 * t364 + t299 * t308 + t333 * t274) * MDP(26); 0.2e1 * t446 * t501 - 0.2e1 * t455 * t500 + MDP(7) * t475 - MDP(8) * t477 + (-pkin(8) * t475 + t411 * t420) * MDP(10) + (pkin(8) * t477 + t414 * t420) * MDP(11) + (-t371 * t451 + t364 * t401 + t366 * t373 - t369 * t424 + (t371 * t498 - t469) * qJD(3)) * MDP(12) + (-t374 * t451 + t365 * t401 + t366 * t376 + t369 * t384 + (t374 * t498 - t468) * qJD(3)) * MDP(13) + (t277 * t424 - t297 * t376 - t298 * t373 - t355 * t364 - t468 * t371 + t469 * t374 + t434) * MDP(14) + (t277 * t355 - t469 * t297 + t468 * t298 + t366 * t504 + t369 * t401 + t491) * MDP(15) + (-t259 * t424 + t266 * t373 + t294 * t364 - t359 * t469 + t472 * t371 + t422 * t405) * MDP(16) + (t260 * t424 - t267 * t373 - t295 * t364 + t469 * t358 - t471 * t371 + t422 * t408) * MDP(17) + (-t472 * t358 + t471 * t359 + (-t259 * t384 - t266 * t376 - t294 * t365) * t408 + (-t260 * t384 - t267 * t376 - t295 * t365) * t405) * MDP(18) + (t259 * t294 + t260 * t295 + t472 * t266 + t471 * t267 + t469 * t292 + t491) * MDP(19) + (-t274 * t332 + t283 * t308) * MDP(20) + (-t274 * t331 + t275 * t332 - t283 * t306 - t284 * t308) * MDP(21) + (-t274 * t424 + t283 * t367 + t308 * t373 - t332 * t364) * MDP(22) + (t275 * t424 - t284 * t367 - t306 * t373 - t331 * t364) * MDP(23) + (-t364 * t424 + t367 * t373) * MDP(24) + (t433 * t364 - (-t258 * t410 + t256) * t424 + t254 * t373 + t325 * t275 + t268 * t331 + t281 * t284 + (t473 * t410 + t474 * t413) * t367 + t470 * t306 + (t255 * t424 - t432 * t367) * qJD(6)) * MDP(25) + (-t432 * t364 + (t257 * t410 + t258 * t413) * t424 - t255 * t373 + t325 * t274 - t268 * t332 + t281 * t283 + (-t474 * t410 + t473 * t413) * t367 + t470 * t308 + (t254 * t424 - t433 * t367) * qJD(6)) * MDP(26); t417 * t500 + t414 * t427 * MDP(11) + (qJD(3) * t300 - t366 * t374 - t371 * t454 - t276) * MDP(12) + (qJD(3) * t302 + t366 * t371 - t374 * t454 - t277) * MDP(13) + ((t298 - t300) * t374 + (-t297 + t302) * t371 + (-t364 * t406 - t493 * t365) * pkin(3)) * MDP(14) + (t297 * t300 - t298 * t302 + (-t493 * t276 + t277 * t406 - t366 * t461) * pkin(3)) * MDP(15) + (-t266 * t374 - t271 * t371 - t276 * t408 + t300 * t359 + t419 * t405) * MDP(16) + (t267 * t374 + t272 * t371 + t276 * t405 - t300 * t358 + t419 * t408) * MDP(17) + (t271 * t358 - t272 * t359 + (qJD(5) * t359 - t266 * t371 + t260) * t408 + (qJD(5) * t358 - t267 * t371 - t259) * t405) * MDP(18) + (-t266 * t271 - t267 * t272 + t276 * t400 - t292 * t300 + (-t259 * t405 + t260 * t408) * t397 + t435 * qJD(5)) * MDP(19) + (t274 * t385 - t466 * t308) * MDP(20) + (-t274 * t383 - t275 * t385 + t466 * t306 - t465 * t308) * MDP(21) + (-t489 - t499) * MDP(22) + (t438 + t490) * MDP(23) - t367 * t374 * MDP(24) + ((-t380 * t413 - t381 * t410) * t364 + t389 * t275 + t268 * t383 - t254 * t374 - t285 * t306 + (t425 * t410 - t426 * t413) * t367 + t465 * t281) * MDP(25) + (-(-t380 * t410 + t381 * t413) * t364 + t389 * t274 + t268 * t385 + t255 * t374 - t285 * t308 + (t426 * t410 + t425 * t413) * t367 - t466 * t281) * MDP(26) + (MDP(10) * t427 - t417 * t501) * t411; 0.2e1 * t374 * qJD(3) * MDP(12) + (-t392 + (t394 - t371) * qJD(3)) * MDP(13) + (-t374 ^ 2 - t370) * MDP(14) + (t297 * t374 + t298 * t371 + t369) * MDP(15) + (t359 * t374 + t364 * t408 - t370 * t405) * MDP(16) + (-t358 * t374 - t364 * t405 - t370 * t408) * MDP(17) + ((t358 * t405 + t359 * t408) * t371 + (-t405 ^ 2 - t408 ^ 2) * t365) * MDP(18) + (t259 * t408 + t260 * t405 - t292 * t374 + t435 * t371) * MDP(19) + (t438 - t490) * MDP(25) + (-t489 + t499) * MDP(26); (t358 * t371 + t487) * MDP(16) + (t359 * t371 + t486) * MDP(17) + (-t358 ^ 2 - t359 ^ 2) * MDP(18) + (t266 * t358 - t267 * t359 + t276) * MDP(19) + (t308 * t367 + t275) * MDP(25) + (t367 * t502 + (-t358 * t367 + t429) * t410 + t467) * MDP(26); -t306 ^ 2 * MDP(21) + (t306 * t367 + t467) * MDP(22) + t364 * MDP(24) + (t255 * t367 + t256) * MDP(25) + (t254 * t367 + t281 * t306) * MDP(26) + (MDP(20) * t306 + t308 * MDP(21) + MDP(23) * t367 - t281 * MDP(25)) * t308 + (t429 * MDP(23) - MDP(25) * t458 + t441 * MDP(26)) * t413 + (t429 * MDP(22) + (-qJD(6) * t359 - t486) * MDP(23) + t441 * MDP(25) + (-t257 + t458) * MDP(26)) * t410;];
tauc = t1;
