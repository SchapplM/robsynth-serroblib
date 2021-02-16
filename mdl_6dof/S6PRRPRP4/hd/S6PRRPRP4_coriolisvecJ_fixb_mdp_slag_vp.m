% Calculate Coriolis joint torque vector for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:12:58
% EndTime: 2021-01-16 03:13:09
% DurationCPUTime: 4.11s
% Computational Cost: add. (2483->395), mult. (5835->528), div. (0->0), fcn. (3797->8), ass. (0->185)
t483 = pkin(4) + pkin(8);
t377 = -pkin(3) - pkin(9);
t372 = sin(qJ(3));
t448 = qJD(2) * t372;
t373 = sin(qJ(2));
t369 = sin(pkin(6));
t452 = qJD(1) * t369;
t423 = t373 * t452;
t338 = qJD(2) * pkin(8) + t423;
t375 = cos(qJ(3));
t370 = cos(pkin(6));
t451 = qJD(1) * t370;
t305 = t372 * t338 - t375 * t451;
t493 = qJD(4) + t305;
t435 = pkin(4) * t448 + t493;
t283 = t377 * qJD(3) + t435;
t410 = -qJ(4) * t372 - pkin(2);
t324 = t377 * t375 + t410;
t376 = cos(qJ(2));
t422 = t376 * t452;
t295 = t324 * qJD(2) - t422;
t371 = sin(qJ(5));
t374 = cos(qJ(5));
t264 = t283 * t371 + t295 * t374;
t447 = qJD(2) * t375;
t416 = t371 * t447;
t443 = qJD(3) * t374;
t329 = -t416 + t443;
t433 = qJD(2) * qJD(3);
t413 = t372 * t433;
t349 = t371 * t413;
t445 = qJD(3) * t371;
t327 = t374 * t447 + t445;
t440 = qJD(5) * t327;
t301 = -t349 + t440;
t442 = qJD(3) * t375;
t449 = qJD(2) * t369;
t415 = qJD(1) * t449;
t399 = t376 * t415;
t446 = qJD(3) * t370;
t494 = qJD(1) * t446 + t399;
t280 = t338 * t442 + t494 * t372;
t412 = t375 * t433;
t273 = pkin(4) * t412 + t280;
t397 = pkin(9) * t372 - qJ(4) * t375;
t441 = qJD(4) * t372;
t381 = t397 * qJD(3) - t441;
t456 = pkin(3) * t413 + t373 * t415;
t281 = t381 * qJD(2) + t456;
t405 = t374 * t273 - t281 * t371;
t390 = qJ(6) * t301 + t405;
t403 = pkin(5) * t412;
t251 = -t264 * qJD(5) - qJD(6) * t329 + t390 + t403;
t259 = -qJ(6) * t327 + t264;
t358 = qJD(5) + t448;
t481 = t259 * t358;
t496 = t251 + t481;
t438 = qJD(5) * t374;
t455 = qJD(5) * t416 + t374 * t413;
t302 = qJD(3) * t438 - t455;
t439 = qJD(5) * t371;
t402 = -t371 * t273 - t374 * t281 - t283 * t438 + t295 * t439;
t385 = -qJ(6) * t302 - t402;
t252 = -qJD(6) * t327 + t385;
t263 = t374 * t283 - t295 * t371;
t258 = -qJ(6) * t329 + t263;
t257 = pkin(5) * t358 + t258;
t495 = -t257 * t358 + t252;
t367 = t372 ^ 2;
t368 = t375 ^ 2;
t492 = MDP(6) * (t367 - t368);
t475 = t329 * t358;
t490 = t257 - t258;
t489 = t327 * t358 + t301;
t488 = t302 + t475;
t335 = t483 * t442;
t466 = t374 * t376;
t487 = -(-t371 * t373 + t372 * t466) * t452 + t374 * t335;
t444 = qJD(3) * t372;
t361 = pkin(3) * t444;
t308 = t361 + t381;
t344 = t483 * t372;
t470 = t371 * t376;
t486 = (t372 * t470 + t373 * t374) * t452 - t374 * t308 + t324 * t439 - t371 * t335 - t344 * t438;
t457 = t374 * t324 + t371 * t344;
t306 = t375 * t338 + t372 * t451;
t366 = qJD(3) * qJ(4);
t299 = -t366 - t306;
t298 = -qJD(3) * pkin(3) + t493;
t430 = MDP(21) + MDP(23);
t386 = -qJ(4) * t442 - t441;
t288 = t386 * qJD(2) + t456;
t314 = t361 + t386;
t378 = qJD(3) ^ 2;
t484 = (-t314 + t423) * qJD(2) - pkin(8) * t378 - t288;
t482 = qJD(2) * pkin(2);
t365 = qJD(3) * qJD(4);
t425 = -t338 * t444 + t494 * t375;
t275 = -t365 - t425;
t269 = -pkin(4) * t413 - t275;
t480 = t269 * t371;
t479 = t269 * t374;
t478 = t301 * t374;
t342 = -pkin(3) * t375 + t410;
t450 = qJD(2) * t342;
t307 = -t422 + t450;
t477 = t307 * t373;
t474 = t329 * t375;
t473 = t358 * t377;
t472 = t369 * t373;
t471 = t371 * t372;
t469 = t372 * t374;
t468 = t372 * t378;
t467 = t374 * t375;
t465 = t375 * t378;
t464 = qJ(6) - t377;
t406 = qJ(6) * t375 - t324;
t436 = qJD(6) * t375;
t463 = pkin(5) * t442 + t406 * t438 + (-qJ(6) * t444 - qJD(5) * t344 - t308 + t436) * t371 + t487;
t437 = qJD(5) * t375;
t419 = t371 * t437;
t462 = -t374 * t436 + (t372 * t443 + t419) * qJ(6) - t486;
t294 = pkin(4) * t447 + t306;
t362 = pkin(3) * t448;
t312 = t397 * qJD(2) + t362;
t404 = t374 * t294 - t312 * t371;
t461 = (pkin(5) * t375 - qJ(6) * t471) * qJD(2) + t404 + qJD(6) * t374 - t464 * t439;
t341 = t464 * t374;
t458 = t371 * t294 + t374 * t312;
t460 = qJ(6) * t374 * t448 + qJD(5) * t341 + qJD(6) * t371 + t458;
t424 = -pkin(5) * t374 - pkin(4);
t459 = pkin(5) * t438 - t424 * t448 + t493;
t345 = t483 * t375;
t285 = t366 + t294;
t408 = -pkin(5) * t327 - qJD(6);
t272 = t285 - t408;
t453 = MDP(26) * t272;
t434 = MDP(12) * qJD(2);
t432 = -MDP(10) + MDP(13);
t431 = MDP(11) - MDP(14);
t429 = MDP(22) + MDP(24);
t428 = t358 * t469;
t427 = t372 * t472;
t379 = qJD(2) ^ 2;
t426 = t372 * t375 * t379;
t421 = t373 * t449;
t420 = t376 * t449;
t418 = t358 * t438;
t417 = t374 * t437;
t411 = MDP(20) * t447;
t401 = t327 * t422;
t400 = t329 * t422;
t396 = -qJD(2) * t368 + t358 * t372;
t395 = t358 * t371;
t393 = qJD(3) * t305 + t425;
t392 = qJD(3) * t306 - t280;
t318 = -t370 * t375 + t427;
t391 = -t318 * t371 + t369 * t466;
t291 = t318 * t374 + t369 * t470;
t319 = t370 * t372 + t375 * t472;
t389 = t285 * t372 + t377 * t442;
t262 = pkin(5) * t302 + t269;
t388 = -t262 * t371 - t272 * t438;
t387 = t262 * t374 - t272 * t439;
t339 = -t422 - t482;
t383 = qJD(3) * (t339 + t422 - t482);
t382 = qJD(3) * (-t307 - t422 - t450);
t380 = -t275 * t375 + t280 * t372 + (t298 * t375 + t299 * t372) * qJD(3);
t359 = pkin(5) * t371 + qJ(4);
t350 = t374 * t412;
t340 = t464 * t371;
t334 = t483 * t444;
t333 = -qJ(4) * t447 + t362;
t332 = t374 * t344;
t326 = t327 ^ 2;
t317 = pkin(5) * t467 + t345;
t297 = t307 * t448;
t296 = -pkin(5) * t419 + (-pkin(8) + t424) * t444;
t290 = t319 * qJD(3) + t372 * t420;
t289 = -qJD(3) * t427 + (t420 + t446) * t375;
t279 = -qJ(6) * t467 + t457;
t274 = pkin(5) * t372 + t406 * t371 + t332;
t261 = t291 * qJD(5) + t290 * t371 + t374 * t421;
t260 = t391 * qJD(5) + t290 * t374 - t371 * t421;
t1 = [(-t275 * t319 + t280 * t318 - t289 * t299 + t290 * t298) * MDP(15) + (-t260 * t329 - t261 * t327 + t291 * t301 + t302 * t391) * MDP(25) + (t251 * t291 - t252 * t391 + t257 * t260 + t259 * t261 + t262 * t319 + t272 * t289) * MDP(26) + t429 * (-t261 * t358 + t289 * t329 - t301 * t319 + t391 * t412) + t430 * (t260 * t358 + t289 * t327 + t291 * t412 + t302 * t319) + (t289 * t375 + t290 * t372) * t434 + (t432 * t290 - t431 * t289 + (t318 * t375 - t319 * t372) * t434) * qJD(3) + (-t288 * t376 * MDP(15) + (MDP(15) * t477 + (t432 * t372 - t431 * t375) * t376 * qJD(3)) * qJD(2) + (-t376 * MDP(4) + (t431 * t372 + t432 * t375 - MDP(3)) * t373) * t379) * t369; 0.2e1 * t372 * MDP(5) * t412 - 0.2e1 * t433 * t492 + MDP(7) * t465 - MDP(8) * t468 + (-pkin(8) * t465 + t372 * t383) * MDP(10) + (pkin(8) * t468 + t375 * t383) * MDP(11) + ((-t367 - t368) * t399 + t380) * MDP(12) + (t372 * t382 - t484 * t375) * MDP(13) + (t484 * t372 + t375 * t382) * MDP(14) + (t288 * t342 + t307 * t314 + (-t477 + (-t298 * t372 + t299 * t375) * t376) * t452 + t380 * pkin(8)) * MDP(15) + (t301 * t371 * t375 + (t371 * t444 - t417) * t329) * MDP(16) + ((-t327 * t371 + t329 * t374) * t444 + (t478 + t302 * t371 + (t327 * t374 + t329 * t371) * qJD(5)) * t375) * MDP(17) + (-t358 * t417 - t301 * t372 + (t396 * t371 + t474) * qJD(3)) * MDP(18) + (t358 * t419 - t302 * t372 + (-t327 * t375 + t396 * t374) * qJD(3)) * MDP(19) + (t358 + t448) * MDP(20) * t442 + (t345 * t302 - t334 * t327 + (-t285 * t443 + t405) * t372 + (-t308 * t371 + t487) * t358 + (-t264 * t372 - t457 * t358) * qJD(5) + (-t401 - t285 * t439 + t479 + ((-t324 * t371 + t332) * qJD(2) + t263) * qJD(3)) * t375) * MDP(21) + (-t345 * t301 - t334 * t329 + (t285 * t445 + t402) * t372 + t486 * t358 + (-t400 - t285 * t438 - t480 + (-t457 * qJD(2) - t264) * qJD(3)) * t375) * MDP(22) + (t296 * t327 + t302 * t317 + (-t272 * t443 + t251) * t372 + t463 * t358 + (-t401 + (qJD(2) * t274 + t257) * qJD(3) + t387) * t375) * MDP(23) + (t296 * t329 - t301 * t317 + (t272 * t445 - t252) * t372 - t462 * t358 + (-t400 + (-qJD(2) * t279 - t259) * qJD(3) + t388) * t375) * MDP(24) + (t274 * t301 - t279 * t302 - t463 * t329 - t462 * t327 + (-t257 * t371 + t259 * t374) * t444 + (t251 * t371 - t252 * t374 + (t257 * t374 + t259 * t371) * qJD(5)) * t375) * MDP(25) + (t251 * t274 + t252 * t279 + t262 * t317 + (-t375 * t422 + t296) * t272 + t462 * t259 + t463 * t257) * MDP(26); -MDP(5) * t426 + t379 * t492 + (-t339 * t448 + t392) * MDP(10) + (-t339 * t447 - t393) * MDP(11) + (-t333 * t447 + t297 - t392) * MDP(13) + (0.2e1 * t365 + (t307 * t375 + t333 * t372) * qJD(2) + t393) * MDP(14) + (-pkin(3) * t280 - qJ(4) * t275 - t298 * t306 - t299 * t493 - t307 * t333) * MDP(15) + (-t329 * t395 - t478) * MDP(16) + (t489 * t371 - t488 * t374) * MDP(17) + (-t358 * t439 + t350 + (-t358 * t471 - t474) * qJD(2)) * MDP(18) + (-t418 + (-t428 + (t327 - t445) * t375) * qJD(2)) * MDP(19) - t358 * t411 + (qJ(4) * t302 + t480 - t404 * t358 + t435 * t327 + (t285 * t374 - t371 * t473) * qJD(5) + (-t263 * t375 + t374 * t389) * qJD(2)) * MDP(21) + (-qJ(4) * t301 + t479 + t458 * t358 + t435 * t329 + (-t285 * t371 - t374 * t473) * qJD(5) + (t264 * t375 - t371 * t389) * qJD(2)) * MDP(22) + (t302 * t359 - t461 * t358 + t459 * t327 + (t272 * t469 + (-qJD(3) * t341 - t257) * t375) * qJD(2) - t388) * MDP(23) + (-t301 * t359 + t460 * t358 + t459 * t329 + (-t272 * t471 + (qJD(3) * t340 + t259) * t375) * qJD(2) + t387) * MDP(24) + (-t301 * t341 + t302 * t340 + t460 * t327 + t461 * t329 - t495 * t371 - t496 * t374) * MDP(25) + (-t251 * t341 - t252 * t340 - t257 * t461 - t460 * t259 + t262 * t359 + t272 * t459) * MDP(26); MDP(13) * t426 + (-t367 * t379 - t378) * MDP(14) + (qJD(3) * t299 + t280 + t297) * MDP(15) - qJD(3) * t453 + t429 * (-t418 - qJD(3) * t329 + (-t371 * t442 - t428) * qJD(2)) + t430 * (-qJD(3) * t327 - t358 * t395 + t350) + ((-t327 * t448 + t301 - t440) * MDP(25) + t496 * MDP(26)) * t374 + ((-t302 + t475) * MDP(25) + t495 * MDP(26)) * t371; -t326 * MDP(17) + t349 * MDP(18) + t455 * MDP(19) + qJD(3) * t411 + (t264 * t358 + t405) * MDP(21) + (t263 * t358 + t402) * MDP(22) + (t390 + 0.2e1 * t403 + t481) * MDP(23) + (t258 * t358 - t385) * MDP(24) + t301 * pkin(5) * MDP(25) + (pkin(5) * t251 + t490 * t259) * MDP(26) + (t358 * MDP(18) + t285 * MDP(22) + (qJD(6) + t272) * MDP(24) - t490 * MDP(25)) * t327 + (-t327 * MDP(18) - MDP(19) * t443 - t430 * t264) * qJD(5) + (t327 * MDP(16) + t358 * MDP(19) - t285 * MDP(21) + (-t272 + t408) * MDP(23) - pkin(5) * t453 + (-MDP(24) * pkin(5) + MDP(17)) * t329) * t329; t488 * MDP(23) - t489 * MDP(24) + (-t329 ^ 2 - t326) * MDP(25) + (t257 * t329 + t259 * t327 + t262) * MDP(26);];
tauc = t1;
