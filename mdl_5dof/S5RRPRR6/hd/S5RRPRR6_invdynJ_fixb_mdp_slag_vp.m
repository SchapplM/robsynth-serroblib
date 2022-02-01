% Calculate vector of inverse dynamics joint torques for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:32
% EndTime: 2022-01-20 11:17:37
% DurationCPUTime: 3.36s
% Computational Cost: add. (2504->335), mult. (3820->460), div. (0->0), fcn. (2513->14), ass. (0->190)
t376 = qJDD(1) + qJDD(2);
t379 = qJD(1) + qJD(2);
t384 = sin(pkin(9));
t390 = cos(qJ(5));
t391 = cos(qJ(4));
t453 = qJD(4) + qJD(5);
t424 = t453 * t391;
t386 = sin(qJ(5));
t479 = t386 * t391;
t387 = sin(qJ(4));
t491 = t376 * t387;
t459 = qJD(4) * t387;
t436 = t384 * t459;
t480 = t386 * t387;
t443 = t384 * t480;
t510 = -qJD(5) * t443 - t386 * t436;
t266 = t510 * t379 + (t376 * t479 + (t379 * t424 + t491) * t390) * t384;
t388 = sin(qJ(2));
t498 = pkin(1) * qJD(1);
t450 = t388 * t498;
t340 = qJ(3) * t379 + t450;
t377 = t384 ^ 2;
t385 = cos(pkin(9));
t487 = t379 * t385;
t348 = -qJD(4) + t487;
t455 = qJD(4) + t348;
t511 = t340 * (t377 * t379 + t385 * t455);
t392 = cos(qJ(2));
t463 = qJD(1) * t392;
t420 = -pkin(1) * t463 + qJD(3);
t488 = t379 * t384;
t343 = -pkin(3) * t385 - pkin(7) * t384 - pkin(2);
t458 = qJD(4) * t391;
t460 = qJD(3) * t385;
t481 = t385 * t392;
t507 = -(t387 * t388 + t391 * t481) * t498 + t343 * t458 + t391 * t460;
t502 = pkin(1) * t392;
t327 = t343 - t502;
t363 = pkin(1) * t388 + qJ(3);
t482 = t385 * t391;
t474 = t387 * t327 + t363 * t482;
t506 = qJ(3) * t482 + t387 * t343;
t383 = qJ(1) + qJ(2);
t372 = sin(t383);
t366 = g(1) * t372;
t462 = qJD(2) * t388;
t448 = pkin(1) * t462;
t470 = -qJD(1) * t448 + qJDD(1) * t502;
t505 = t366 + t470;
t503 = qJD(4) * t506 + (-t387 * t481 + t388 * t391) * t498 + t387 * t460;
t501 = pkin(2) * t376;
t374 = cos(t383);
t500 = g(2) * t374;
t499 = g(3) * t384;
t298 = t343 * t379 + t420;
t446 = t340 * t482;
t408 = -t298 * t387 - t446;
t485 = t384 * t387;
t452 = pkin(8) * t485;
t279 = -t379 * t452 - t408;
t497 = t279 * t390;
t342 = -qJD(5) + t348;
t496 = t342 * t385;
t495 = t363 * t387;
t494 = t372 * t385;
t493 = t374 * t385;
t492 = t376 * t385;
t490 = t376 * t391;
t454 = qJDD(1) * t388;
t461 = qJD(2) * t392;
t306 = qJ(3) * t376 + qJD(3) * t379 + (qJD(1) * t461 + t454) * pkin(1);
t299 = t377 * t306;
t486 = t379 * t387;
t484 = t384 * t391;
t483 = t385 * t387;
t478 = t387 * t391;
t477 = t390 * t391;
t472 = t374 * pkin(2) + t372 * qJ(3);
t471 = g(1) * t374 + g(2) * t372;
t378 = t385 ^ 2;
t469 = t377 + t378;
t381 = t391 ^ 2;
t468 = t387 ^ 2 - t381;
t467 = MDP(10) * t377;
t466 = MDP(11) * t377;
t465 = MDP(12) * t384;
t464 = MDP(13) * t384;
t457 = qJD(5) * t386;
t347 = -qJDD(4) + t492;
t341 = -qJDD(5) + t347;
t334 = t341 * MDP(21);
t456 = t347 * MDP(14);
t451 = pkin(8) * t484;
t447 = qJ(3) * t483;
t444 = t379 * t484;
t442 = t384 * t477;
t433 = qJDD(3) - t470;
t290 = t343 * t376 + t433;
t441 = t387 * t290 + t298 * t458 + t306 * t482;
t356 = t461 * pkin(1) + qJD(3);
t440 = t327 * t458 + t356 * t482 + t387 * t448;
t439 = t379 * t462;
t437 = t379 * t458;
t435 = t385 * t459;
t434 = t340 * t459;
t432 = -pkin(2) * t372 + t374 * qJ(3);
t312 = t433 - t501;
t431 = -t312 - t500;
t430 = t306 * t469;
t429 = t356 * t469;
t428 = t469 * t376;
t427 = t347 + t492;
t426 = t379 * t455;
t402 = -t385 * t434 + t441;
t407 = t437 + t491;
t260 = -pkin(8) * t384 * t407 + t402;
t293 = t391 * t298;
t278 = -pkin(8) * t444 - t340 * t483 + t293;
t268 = -pkin(4) * t348 + t278;
t425 = qJD(5) * t268 + t260;
t423 = t378 * t306 + t299 - t471;
t287 = t391 * t290;
t419 = -t306 * t483 + t287;
t331 = t391 * t343;
t291 = -t451 + t331 + (-qJ(3) * t387 - pkin(4)) * t385;
t418 = qJD(5) * t291 + (-t447 - t451) * qJD(4) + t507;
t297 = -t452 + t506;
t352 = pkin(8) * t436;
t417 = qJD(5) * t297 - t352 + t503;
t416 = -t268 * t386 - t497;
t324 = t391 * t327;
t281 = -t451 + t324 + (-pkin(4) - t495) * t385;
t289 = -t452 + t474;
t415 = t281 * t390 - t289 * t386;
t414 = t281 * t386 + t289 * t390;
t413 = t387 * t390 + t479;
t412 = -t477 + t480;
t411 = qJD(4) * (t348 + t487);
t353 = t384 * pkin(4) * t458;
t410 = -t420 * t384 - t353;
t316 = t372 * t483 + t374 * t391;
t318 = t372 * t391 - t374 * t483;
t406 = -g(1) * t316 - g(2) * t318 + t391 * t299 + t402 * t385;
t317 = -t372 * t482 + t374 * t387;
t319 = t372 * t387 + t374 * t482;
t405 = t377 * t340 * t458 - g(1) * t317 - g(2) * t319 + t387 * t299;
t404 = t379 * t450 - t500;
t259 = -t376 * t451 - pkin(4) * t347 + (-t446 + (pkin(8) * t488 - t298) * t387) * qJD(4) + t419;
t273 = t279 * t457;
t280 = (pkin(4) * t407 + t306) * t384;
t397 = t453 * t413;
t282 = t397 * t384;
t301 = (pkin(4) * t486 + t340) * t384;
t382 = qJ(4) + qJ(5);
t371 = sin(t382);
t373 = cos(t382);
t307 = t371 * t494 + t373 * t374;
t309 = -t371 * t493 + t372 * t373;
t322 = t412 * t384;
t403 = -g(1) * t307 - g(2) * t309 + (t386 * t259 + t425 * t390 - t273) * t385 - t280 * t322 - t301 * t282;
t265 = t376 * t442 + (-t376 * t480 - t379 * t397) * t384;
t303 = (-t442 + t443) * t379;
t304 = t413 * t488;
t401 = -t303 * t304 * MDP(17) + (-t304 * t342 + t265) * MDP(19) + (t303 * t342 - t266) * MDP(20) + (t303 ^ 2 - t304 ^ 2) * MDP(18) - t334;
t375 = t379 ^ 2;
t400 = -t348 ^ 2 - t375 * t377;
t253 = qJD(5) * t416 + t390 * t259 - t386 * t260;
t283 = t384 * t390 * t424 + t510;
t308 = t371 * t374 - t373 * t494;
t310 = t371 * t372 + t373 * t493;
t321 = t413 * t384;
t399 = -g(1) * t308 - g(2) * t310 - t253 * t385 + t280 * t321 + t301 * t283;
t398 = t420 * t469;
t396 = t273 + t373 * t499 + g(1) * t310 + (t279 * t342 - t259) * t386 - g(2) * t308 + t301 * t304;
t395 = (-t265 * t321 + t266 * t322 + t282 * t304 + t283 * t303) * MDP(18) + (-t265 * t385 + t282 * t342 + t322 * t341) * MDP(19) + (t266 * t385 + t283 * t342 + t321 * t341) * MDP(20) + (-t265 * t322 + t282 * t303) * MDP(17) + (t387 * t411 - t391 * t427) * t465 + (t387 * t427 + t391 * t411) * t464 + 0.2e1 * (qJD(4) * t379 * t468 - t376 * t478) * t466 + (t376 * t381 - 0.2e1 * t387 * t437) * t467 + t376 * MDP(4) + (t334 + t456) * t385;
t394 = -g(1) * t309 + g(2) * t307 + t301 * t303 + t371 * t499 + t253;
t393 = cos(qJ(1));
t389 = sin(qJ(1));
t368 = -pkin(2) - t502;
t359 = pkin(4) * t485;
t355 = t391 * t448;
t351 = g(1) * t494;
t336 = qJ(3) * t384 + t359;
t335 = -pkin(2) * t379 + t420;
t325 = t363 * t384 + t359;
t311 = t356 * t384 + t353;
t275 = -t474 * qJD(4) - t356 * t483 + t352 + t355;
t274 = (-t363 * t483 - t451) * qJD(4) + t440;
t262 = qJD(4) * t408 + t419;
t1 = [((-t363 * t435 + t440) * t348 + t474 * t347 + ((t356 * t379 + t363 * t376) * t391 + (-t363 * t379 - t340) * t459) * t377 + t406) * MDP(16) + (-(-t327 * t459 + t355) * t348 - t324 * t347 + (-(-t356 * t387 - t363 * t458) * t348 + t347 * t495 - t262) * t385 + (t356 * t486 + t363 * t407) * t377 + t405) * MDP(15) + (g(1) * t389 - g(2) * t393) * MDP(2) + (g(1) * t393 + g(2) * t389) * MDP(3) + (-(-qJD(5) * t414 - t274 * t386 + t275 * t390) * t342 - t415 * t341 + t311 * t304 + t325 * t266 + t399) * MDP(22) + (t363 * t428 + t379 * t429 + t423) * MDP(8) + (t312 * t368 + t335 * t448 - g(1) * (-pkin(1) * t389 + t432) - g(2) * (pkin(1) * t393 + t472) + t340 * t429 + t363 * t430) * MDP(9) + (-t500 + (t376 * t392 - t439) * pkin(1) + t505) * MDP(5) + (((-qJDD(1) - t376) * t388 + (-qJD(1) - t379) * t461) * pkin(1) + t471) * MDP(6) + (t351 + (-pkin(1) * t439 - t368 * t376 + t431) * t385) * MDP(7) + qJDD(1) * MDP(1) + t395 + ((qJD(5) * t415 + t274 * t390 + t275 * t386) * t342 + t414 * t341 - t311 * t303 + t325 * t265 + t403) * MDP(23); ((t291 * t386 + t297 * t390) * t341 + t336 * t265 + (-t386 * t417 + t390 * t418) * t342 + t410 * t303 + t403) * MDP(23) + (-(t331 - t447) * t347 - t262 * t385 + t503 * t348 + (qJ(3) * t491 + (qJ(3) * t458 + t387 * t420) * t379) * t377 + t405) * MDP(15) + (t404 + t505) * MDP(5) + (t506 * t347 + (-qJ(3) * t435 + t507) * t348 + (qJ(3) * t490 - t434 + (-qJ(3) * t459 + t391 * t420) * t379) * t377 + t406) * MDP(16) + (-t312 * pkin(2) - g(1) * t432 - g(2) * t472 + qJ(3) * t430 - t335 * t450 + t340 * t398) * MDP(9) + (qJ(3) * t428 + t398 * t379 + t423) * MDP(8) + (-(t291 * t390 - t297 * t386) * t341 + t336 * t266 + (t386 * t418 + t390 * t417) * t342 - t410 * t304 + t399) * MDP(22) + (t351 + (-t312 + t404 + t501) * t385) * MDP(7) + ((-t454 + (-qJD(2) + t379) * t463) * pkin(1) + t471) * MDP(6) + t395; -MDP(7) * t492 - t469 * t375 * MDP(8) + (-t340 * t379 * t469 - t366 - t431) * MDP(9) + (-t347 * t391 + t387 * t400) * MDP(15) + (t347 * t387 + t391 * t400) * MDP(16) + (t397 * t342 + t412 * t341 + (-t384 * t304 - t413 * t496) * t379) * MDP(22) + (t303 * t488 + t413 * t341 + (-t453 * t342 + t496 * t379) * t412) * MDP(23); (-t387 * t426 + t490) * t465 + (-t391 * t426 - t491) * t464 - t456 + (-g(1) * t318 + g(2) * t316 + t287 - t391 * t511 + (-t298 * t455 - t385 * t306 + t499) * t387) * MDP(15) + (g(1) * t319 - g(2) * t317 + g(3) * t484 - t293 * t348 + t387 * t511 - t441) * MDP(16) + ((-t278 * t386 - t497) * t342 + (-t304 * t444 - t341 * t390 + t342 * t457) * pkin(4) + t394) * MDP(22) + ((-t278 * t342 - t425) * t390 + (qJD(5) * t342 * t390 + t303 * t444 + t341 * t386) * pkin(4) + t396) * MDP(23) + t401 + (-t468 * t466 + t467 * t478) * t375; (t342 * t416 + t394) * MDP(22) + ((-t260 + (-qJD(5) - t342) * t268) * t390 + t396) * MDP(23) + t401;];
tau = t1;
