% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:35
% EndTime: 2019-03-08 22:14:44
% DurationCPUTime: 4.20s
% Computational Cost: add. (2368->364), mult. (5653->499), div. (0->0), fcn. (4080->10), ass. (0->174)
t385 = cos(qJ(5));
t381 = sin(qJ(5));
t386 = cos(qJ(3));
t455 = qJD(2) * t386;
t435 = t381 * t455;
t382 = sin(qJ(3));
t458 = qJD(2) * t382;
t328 = t385 * t458 - t435;
t380 = sin(qJ(6));
t449 = qJD(6) * t380;
t332 = t381 * t382 + t385 * t386;
t401 = t332 * qJD(5);
t444 = qJD(2) * qJD(3);
t431 = t386 * t444;
t432 = t382 * t444;
t465 = t381 * t432 + t385 * t431;
t293 = -qJD(2) * t401 + t465;
t372 = qJD(3) - qJD(5);
t384 = cos(qJ(6));
t448 = qJD(6) * t384;
t468 = t384 * t293 - t372 * t448;
t268 = -t328 * t449 + t468;
t311 = t328 * t384 - t372 * t380;
t480 = t293 * t380;
t269 = qJD(6) * t311 + t480;
t450 = qJD(5) * t385;
t452 = qJD(3) * t386;
t398 = t381 * t452 + t382 * t450;
t294 = qJD(2) * t398 - qJD(5) * t435 - t385 * t432;
t477 = t328 * t380;
t309 = t384 * t372 + t477;
t473 = t384 * t294;
t475 = t380 * t294;
t482 = t268 * t380;
t493 = t332 * qJD(2);
t503 = qJD(6) + t493;
t495 = t503 * t311;
t496 = t503 ^ 2;
t510 = t309 * t503;
t511 = -((t269 + t495) * t380 + (-t268 + t510) * t384) * MDP(24) + (t384 * t495 + t482) * MDP(23) + (-t311 * t328 + t384 * t496 + t475) * MDP(25) - (-t309 * t328 + t380 * t496 - t473) * MDP(26) - (t328 * t372 + t294) * MDP(19) - (-t328 ^ 2 + t493 ^ 2) * MDP(17) + (MDP(16) * t493 - MDP(27) * t503) * t328;
t383 = sin(qJ(2));
t377 = sin(pkin(6));
t460 = qJD(1) * t377;
t437 = t383 * t460;
t342 = qJD(2) * pkin(8) + t437;
t329 = t382 * t342;
t378 = cos(pkin(6));
t459 = qJD(1) * t378;
t356 = t386 * t459;
t314 = -t329 + t356;
t504 = qJD(4) - t314;
t505 = -pkin(9) * t458 + t504;
t375 = t382 ^ 2;
t376 = t386 ^ 2;
t498 = MDP(6) * (t375 - t376);
t388 = -pkin(3) - pkin(4);
t438 = t388 * qJD(3);
t292 = t438 + t505;
t315 = t386 * t342 + t382 * t459;
t303 = -pkin(9) * t455 + t315;
t374 = qJD(3) * qJ(4);
t296 = t303 + t374;
t266 = t292 * t385 - t296 * t381;
t264 = pkin(5) * t372 - t266;
t497 = t264 * t503;
t416 = t382 * t438;
t368 = t382 * qJD(4);
t462 = qJ(4) * t452 + t368;
t412 = t416 + t462 + t437;
t308 = t374 + t315;
t305 = -qJD(3) * pkin(3) + t504;
t373 = qJD(3) * qJD(4);
t387 = cos(qJ(2));
t355 = t387 * t460;
t414 = qJD(2) * t355;
t466 = qJD(3) * t356 + t386 * t414;
t439 = t373 + t466;
t453 = qJD(3) * t382;
t279 = (pkin(9) * qJD(2) - t342) * t453 + t439;
t433 = qJD(1) * t453;
t289 = t342 * t452 + t378 * t433 + t382 * t414;
t281 = -pkin(9) * t431 + t289;
t451 = qJD(5) * t381;
t256 = t381 * t279 - t281 * t385 + t292 * t451 + t296 * t450;
t295 = pkin(5) * t328 + pkin(10) * t493;
t364 = qJ(4) * t455;
t320 = t388 * t458 + t364;
t411 = qJ(4) * t385 + t381 * t388;
t336 = -pkin(10) + t411;
t492 = t503 * (qJD(6) * t336 - t295 + t320) - t256;
t491 = t256 + (pkin(10) * qJD(6) + t295) * t503;
t255 = t385 * t279 + t381 * t281 + t292 * t450 - t296 * t451;
t379 = qJD(2) * pkin(2);
t343 = -t355 - t379;
t316 = -pkin(3) * t455 - qJ(4) * t458 + t343;
t304 = pkin(4) * t455 - t316;
t272 = pkin(5) * t493 - pkin(10) * t328 + t304;
t344 = -t386 * pkin(3) - t382 * qJ(4) - pkin(2);
t331 = t386 * pkin(4) - t344;
t404 = t381 * t386 - t382 * t385;
t282 = pkin(5) * t332 + pkin(10) * t404 + t331;
t298 = qJD(3) * t332 - t401;
t486 = pkin(8) - pkin(9);
t345 = t486 * t382;
t346 = t486 * t386;
t313 = t345 * t381 + t346 * t385;
t337 = t486 * t453;
t338 = qJD(3) * t346;
t405 = t345 * t385 - t346 * t381;
t470 = qJD(5) * t405 - t332 * t355 - t337 * t385 + t338 * t381;
t489 = -(qJD(6) * t272 + t255) * t332 - t256 * t404 + t264 * t298 + (-qJD(6) * t282 - t470) * t503 - t313 * t294;
t487 = -(t314 + t329) * qJD(3) + t466;
t267 = t292 * t381 + t296 * t385;
t265 = -pkin(10) * t372 + t267;
t409 = t265 * t380 - t272 * t384;
t485 = t409 * t328;
t258 = t265 * t384 + t272 * t380;
t484 = t258 * t328;
t483 = t264 * t404;
t481 = t282 * t294;
t479 = t493 * t372;
t476 = t377 * t386;
t389 = qJD(3) ^ 2;
t474 = t382 * t389;
t472 = t386 * t389;
t410 = -qJ(4) * t381 + t385 * t388;
t471 = -qJD(5) * t410 + t303 * t381 - t505 * t385;
t469 = qJD(5) * t313 - t337 * t381 - t338 * t385 - t404 * t355;
t467 = qJD(5) * t411 + t303 * t385 + t505 * t381;
t463 = -qJ(4) * t431 - qJD(2) * t368;
t457 = qJD(2) * t383;
t456 = qJD(2) * t384;
t454 = qJD(2) * t387;
t447 = qJD(6) * t387;
t443 = -MDP(10) - MDP(12);
t442 = pkin(3) * t453;
t441 = t377 * t382 * t383;
t390 = qJD(2) ^ 2;
t440 = t382 * t386 * t390;
t436 = t382 * t454;
t430 = t343 - t379;
t422 = qJD(2) * t344 + t316;
t419 = t372 ^ 2;
t418 = t372 * t385;
t417 = qJD(3) * t443;
t415 = t454 * t476;
t299 = -t385 * t453 - t386 * t451 + t398;
t413 = pkin(5) * t299 - pkin(10) * t298 + t412;
t324 = -t378 * t386 + t441;
t325 = t378 * t382 + t383 * t476;
t406 = t324 * t385 - t325 * t381;
t287 = t324 * t381 + t325 * t385;
t403 = qJD(3) * t315 - t289;
t402 = t298 * t384 + t404 * t449;
t397 = t304 * t328 + t256;
t396 = -pkin(10) * t294 + t266 * t503 + t497;
t395 = -t304 * t493 + t255;
t297 = (t437 + t442) * qJD(2) + t463;
t321 = t442 - t462;
t394 = -pkin(8) * t389 - t297 + (-t321 + t437) * qJD(2);
t393 = -t336 * t294 + t471 * t503 - t497;
t283 = -t342 * t453 + t439;
t392 = t283 * t386 + t289 * t382 + (t305 * t386 - t308 * t382) * qJD(3);
t290 = (t416 - t437) * qJD(2) - t463;
t339 = t377 * t387 * t433;
t335 = pkin(5) - t410;
t334 = pkin(3) * t458 - t364;
t301 = qJD(3) * t325 + t377 * t436;
t300 = -qJD(3) * t324 + t415;
t262 = qJD(5) * t406 + t300 * t385 + t301 * t381;
t261 = qJD(5) * t287 + t300 * t381 - t301 * t385;
t260 = pkin(5) * t294 - pkin(10) * t293 + t290;
t259 = t384 * t260;
t1 = [(t283 * t325 + t289 * t324 + t300 * t308 + t301 * t305) * MDP(15) + ((-t262 * t380 - t287 * t448) * t503 - t287 * t475 + t261 * t309 - t406 * t269) * MDP(28) + (-(t262 * t384 - t287 * t449) * t503 - t287 * t473 + t261 * t311 - t406 * t268) * MDP(29) + (MDP(21) * t261 + MDP(22) * t262) * t372 + (-MDP(11) + MDP(14)) * (-t390 * t441 + (t300 + t415) * qJD(3)) + t301 * t417 + (t300 * t386 + t301 * t382 + (t324 * t386 - t325 * t382) * qJD(3)) * MDP(13) * qJD(2) + ((-t297 * t387 + t316 * t457) * MDP(15) + (t294 * t387 - t457 * t493) * MDP(21) + (t293 * t387 - t328 * t457) * MDP(22) + ((-t380 * t447 - t383 * t456) * t503 + t387 * t473) * MDP(28) + (-(-t380 * t457 + t384 * t447) * t503 - t387 * t475) * MDP(29) + t417 * t436 + (-t387 * MDP(4) + (t386 * t443 - MDP(3)) * t383) * t390) * t377; 0.2e1 * t382 * MDP(5) * t431 - 0.2e1 * t444 * t498 + MDP(7) * t472 - MDP(8) * t474 + (-pkin(8) * t472 + t430 * t453 + t339) * MDP(10) + (pkin(8) * t474 + (t430 + t355) * t452) * MDP(11) + (t386 * t394 + t422 * t453 + t339) * MDP(12) + ((-t375 - t376) * t414 + t392) * MDP(13) + ((-t422 - t355) * t452 + t394 * t382) * MDP(14) + (t297 * t344 + t316 * t321 + (-t316 * t383 + (-t305 * t382 - t308 * t386) * t387) * t460 + t392 * pkin(8)) * MDP(15) + (-t293 * t404 + t298 * t328) * MDP(16) + (-t293 * t332 + t294 * t404 - t298 * t493 - t299 * t328) * MDP(17) + (t290 * t332 + t294 * t331 + t299 * t304 + t412 * t493) * MDP(21) + (-t290 * t404 + t293 * t331 + t298 * t304 + t328 * t412) * MDP(22) + (-t268 * t384 * t404 + t311 * t402) * MDP(23) + ((-t309 * t384 - t311 * t380) * t298 - (-t482 - t269 * t384 + (t309 * t380 - t311 * t384) * qJD(6)) * t404) * MDP(24) + (t268 * t332 + t299 * t311 + t402 * t503 - t404 * t473) * MDP(25) + (t404 * t475 - t269 * t332 - t299 * t309 + (-t298 * t380 + t404 * t448) * t503) * MDP(26) + (t294 * t332 + t299 * t503) * MDP(27) + (-t409 * t299 + t259 * t332 - t405 * t269 + t469 * t309 + (t481 + t413 * t503 + (-t265 * t332 - t313 * t503 - t483) * qJD(6)) * t384 + t489 * t380) * MDP(28) + (-t258 * t299 - t405 * t268 + t469 * t311 + (-t481 - (-qJD(6) * t265 + t260) * t332 + qJD(6) * t483 + (qJD(6) * t313 - t413) * t503) * t380 + t489 * t384) * MDP(29) + (-MDP(18) * t298 + MDP(19) * t299 + MDP(21) * t469 + MDP(22) * t470) * t372; -MDP(5) * t440 + t390 * t498 + (-t343 * t458 + t403) * MDP(10) + (-t343 * t455 - t487) * MDP(11) + ((-t316 * t382 + t334 * t386) * qJD(2) + t403) * MDP(12) + (0.2e1 * t373 + (t316 * t386 + t334 * t382) * qJD(2) + t487) * MDP(14) + (-pkin(3) * t289 + qJ(4) * t283 - t305 * t315 + t504 * t308 - t316 * t334) * MDP(15) + (qJD(5) * t493 - t465 + t479) * MDP(18) + (-t320 * t493 + t372 * t467 + t397) * MDP(21) + (-t320 * t328 - t372 * t471 + t395) * MDP(22) + (t335 * t269 + t467 * t309 + t393 * t380 - t492 * t384 - t485) * MDP(28) + (t335 * t268 + t467 * t311 + t492 * t380 + t393 * t384 - t484) * MDP(29) - t511; -MDP(12) * t440 + (-t375 * t390 - t389) * MDP(14) + (-qJD(3) * t308 + t316 * t458 + t289) * MDP(15) + (-t381 * t419 - t458 * t493) * MDP(21) + (-t328 * t458 - t385 * t419) * MDP(22) + (-t385 * t269 + (t380 * t418 - t382 * t456) * t503 + (-t309 * t372 - t448 * t503 - t475) * t381) * MDP(28) + (-t385 * t268 + (t380 * t458 + t384 * t418) * t503 + (-t311 * t372 + t449 * t503 - t473) * t381) * MDP(29); (t293 - t479) * MDP(18) + (-t267 * t372 - t397) * MDP(21) + (-t266 * t372 - t395) * MDP(22) + (-pkin(5) * t269 - t267 * t309 + t396 * t380 - t491 * t384 + t485) * MDP(28) + (-pkin(5) * t268 - t267 * t311 + t491 * t380 + t396 * t384 + t484) * MDP(29) + t511; t311 * t309 * MDP(23) + (-t309 ^ 2 + t311 ^ 2) * MDP(24) + (t468 + t510) * MDP(25) + (-t480 + t495) * MDP(26) + t294 * MDP(27) + (-t255 * t380 + t258 * t503 - t264 * t311 + t259) * MDP(28) + (-t255 * t384 - t260 * t380 + t264 * t309 - t409 * t503) * MDP(29) + (-MDP(25) * t477 - MDP(26) * t311 - MDP(28) * t258 + MDP(29) * t409) * qJD(6);];
tauc  = t1;
