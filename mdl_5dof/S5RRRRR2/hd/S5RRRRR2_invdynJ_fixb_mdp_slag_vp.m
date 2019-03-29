% Calculate vector of inverse dynamics joint torques for
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:39
% EndTime: 2019-03-29 15:26:42
% DurationCPUTime: 2.85s
% Computational Cost: add. (2536->336), mult. (4631->484), div. (0->0), fcn. (3452->14), ass. (0->162)
t382 = cos(qJ(2));
t369 = qJD(1) + qJD(2);
t416 = qJD(1) * (-qJD(2) + t369);
t377 = sin(qJ(2));
t430 = qJDD(1) * t377;
t373 = qJ(1) + qJ(2);
t362 = sin(t373);
t364 = cos(t373);
t444 = g(1) * t364 + g(2) * t362;
t387 = (t382 * t416 - t430) * pkin(1) + t444;
t380 = cos(qJ(4));
t381 = cos(qJ(3));
t449 = t380 * t381;
t375 = sin(qJ(4));
t376 = sin(qJ(3));
t453 = t375 * t376;
t319 = -t449 + t453;
t367 = qJDD(1) + qJDD(2);
t320 = t375 * t381 + t376 * t380;
t368 = qJD(3) + qJD(4);
t485 = t368 * t320;
t273 = t319 * t367 + t369 * t485;
t372 = qJ(3) + qJ(4);
t361 = sin(t372);
t463 = t361 * t364;
t464 = t361 * t362;
t484 = g(1) * t463 + g(2) * t464;
t355 = g(1) * t362;
t429 = qJDD(1) * t382;
t359 = pkin(1) * t429;
t481 = t355 + t359;
t439 = qJD(2) * t382;
t395 = qJD(1) * t439 + t430;
t479 = pkin(1) * t381;
t483 = t395 * t479;
t480 = pkin(1) * t377;
t478 = pkin(2) * t381;
t477 = g(2) * t364;
t353 = g(3) * t361;
t363 = cos(t372);
t476 = g(3) * t363;
t475 = pkin(1) * qJD(1);
t314 = t320 * t369;
t343 = t369 * t449;
t436 = qJD(3) * t381;
t420 = t369 * t436;
t424 = t369 * t453;
t272 = qJD(4) * t343 + t320 * t367 - t368 * t424 + t380 * t420;
t366 = qJDD(3) + qJDD(4);
t374 = sin(qJ(5));
t379 = cos(qJ(5));
t432 = qJD(5) * t379;
t422 = t379 * t272 + t374 * t366 + t368 * t432;
t433 = qJD(5) * t374;
t259 = -t314 * t433 + t422;
t474 = t259 * t374;
t271 = qJDD(5) + t273;
t473 = t271 * t374;
t472 = t271 * t379;
t471 = t271 * t381;
t293 = t314 * t374 - t379 * t368;
t312 = -t343 + t424;
t304 = qJD(5) + t312;
t470 = t293 * t304;
t295 = t314 * t379 + t368 * t374;
t469 = t295 * t304;
t296 = t368 * t319;
t468 = t296 * t374;
t467 = t296 * t379;
t466 = t320 * t374;
t465 = t320 * t379;
t462 = t362 * t363;
t461 = t363 * t364;
t460 = t363 * t374;
t459 = t363 * t379;
t458 = t364 * t374;
t457 = t364 * t379;
t456 = t367 * t381;
t455 = t367 * t382;
t454 = t369 * t376;
t450 = t376 * t381;
t288 = qJDD(3) * pkin(2) + (-t376 * t430 + (-t376 * t439 - t377 * t436) * qJD(1)) * pkin(1);
t440 = qJD(1) * t377;
t426 = pkin(1) * t440;
t414 = t376 * t426;
t407 = qJD(3) * t414;
t448 = -t380 * t288 - t375 * t407;
t425 = qJD(2) * t480;
t412 = qJD(1) * t425;
t446 = (t412 + t477) * t376;
t445 = t481 * t381;
t370 = t376 ^ 2;
t443 = -t381 ^ 2 + t370;
t442 = MDP(14) * t312;
t441 = MDP(25) * t304;
t438 = qJD(3) * t369;
t437 = qJD(3) * t376;
t326 = qJD(3) * pkin(2) - t414;
t435 = qJD(4) * t326;
t434 = qJD(4) * t380;
t431 = qJD(5) * t381;
t428 = qJDD(3) * t381;
t427 = -qJDD(1) - t367;
t419 = t375 * t435;
t263 = t419 + (t375 * t430 + (t375 * t439 + t377 * t434) * qJD(1)) * t479 + t448;
t413 = t381 * t426;
t298 = -t380 * t326 + t375 * t413;
t290 = t298 * t432;
t423 = g(3) * t460 + t263 * t374 + t290;
t289 = t298 * t433;
t421 = t484 * t379 + t289;
t406 = qJD(4) * t413;
t409 = -t380 * t407 + (t288 - t406) * t375;
t262 = (t435 + t483) * t380 + t409;
t418 = -t262 + t353;
t417 = t379 * t304;
t415 = (-qJD(1) - t369) * qJD(2);
t408 = -0.2e1 * t438;
t302 = t320 * t426;
t405 = t298 * t312 - t302 * t304;
t328 = t380 * t413;
t299 = t326 * t375 + t328;
t316 = -t369 * t478 - t382 * t475;
t280 = t299 * t379 + t316 * t374;
t279 = -t299 * t374 + t316 * t379;
t311 = t319 * t480;
t344 = -pkin(1) * t382 - t478;
t404 = -t311 * t379 + t344 * t374;
t403 = t311 * t374 + t344 * t379;
t402 = qJD(5) * t375 + t454;
t401 = -t320 * t433 - t467;
t292 = t412 - t359 + (t369 * t437 - t456) * pkin(2);
t400 = -g(1) * t464 + g(2) * t463 + t292 * t320 - t316 * t296;
t399 = g(1) * t462 - g(2) * t461 + t292 * t319 + t316 * t485;
t398 = -t448 - t476 + t484;
t397 = t320 * t382;
t396 = t319 * t382;
t351 = t379 * t366;
t260 = t295 * qJD(5) + t272 * t374 - t351;
t394 = ((t259 - t470) * t379 + (-t260 - t469) * t374) * MDP(22) + (t295 * t417 + t474) * MDP(21) + (-t304 ^ 2 * t374 + t293 * t314 + t472) * MDP(24) + (-t295 * t314 + t304 * t417 + t473) * MDP(23) + (t312 * t368 + t272) * MDP(16) + (t314 * t368 - t273) * MDP(17) + (-t312 ^ 2 + t314 ^ 2) * MDP(15) + t366 * MDP(18);
t393 = t416 * t480 - t477;
t305 = t362 * t460 + t457;
t307 = t362 * t379 - t363 * t458;
t392 = -g(1) * t305 - g(2) * t307 - (t279 * qJD(5) + t262 * t379 + t292 * t374) * t319 + t263 * t465 - t280 * t485 - t298 * t467;
t287 = t379 * t292;
t306 = -t362 * t459 + t458;
t308 = t362 * t374 + t363 * t457;
t390 = -g(1) * t306 - g(2) * t308 + (-t280 * qJD(5) - t262 * t374 + t287) * t319 + t263 * t466 + t279 * t485 - t298 * t468 + t320 * t290;
t389 = g(1) * t461 + g(2) * t462 + t312 * t316 + t353 - t409;
t384 = qJD(3) ^ 2;
t388 = (-(-t293 * t379 - t295 * t374) * t296 + (-t474 - t260 * t379 + (t293 * t374 - t295 * t379) * qJD(5)) * t320) * MDP(22) + (t259 * t319 + t271 * t465 + t295 * t485 + t401 * t304) * MDP(23) + (-t271 * t466 - t260 * t319 - t293 * t485 + (-t320 * t432 + t468) * t304) * MDP(24) + (-t272 * t319 - t273 * t320 + t296 * t312 - t314 * t485) * MDP(15) + (t259 * t465 + t401 * t295) * MDP(21) + (t271 * t319 + t304 * t485) * MDP(25) + (t272 * t320 - t296 * t314) * MDP(14) + (-t296 * t368 + t320 * t366) * MDP(16) + (-t319 * t366 - t368 * t485) * MDP(17) + 0.2e1 * (t367 * t450 - t443 * t438) * MDP(8) + (t367 * t370 + 0.2e1 * t376 * t420) * MDP(7) + (-t376 * t384 + t428) * MDP(10) + (qJDD(3) * t376 + t381 * t384) * MDP(9) + t367 * MDP(4);
t385 = (-pkin(2) * t368 - t326) * qJD(4) - t483;
t383 = cos(qJ(1));
t378 = sin(qJ(1));
t327 = pkin(2) * t437 + t425;
t310 = t320 * t480;
t303 = t396 * t475;
t301 = t397 * t475;
t300 = t375 * t414 - t328;
t270 = (qJD(2) * t397 - t296 * t377) * pkin(1);
t269 = (-qJD(2) * t396 - t377 * t485) * pkin(1);
t1 = [(-(t269 * t379 + t327 * t374) * t304 - t404 * t271 + t270 * t295 + t310 * t259 + (-t298 * t466 - t403 * t304) * qJD(5) + t392) * MDP(27) + (-t376 * t355 + ((-t428 + (qJD(2) * t369 + t384) * t376) * t377 + (t427 * t376 + t381 * t408) * t382) * pkin(1) + t446) * MDP(13) + (-t381 * t477 + ((t455 + (-t384 + t415) * t377) * t381 + (-qJDD(3) * t377 + t382 * t408) * t376) * pkin(1) + t445) * MDP(12) + (-t477 + (t377 * t415 + t455) * pkin(1) + t481) * MDP(5) + ((t427 * t377 + t382 * t415) * pkin(1) + t444) * MDP(6) + ((-t404 * qJD(5) - t269 * t374 + t327 * t379) * t304 + t403 * t271 + t270 * t293 + t310 * t260 + t390) * MDP(26) + qJDD(1) * MDP(1) + t388 + (-t270 * t368 + t273 * t344 - t310 * t366 + t312 * t327 + t399) * MDP(19) + (-t269 * t368 + t272 * t344 + t311 * t366 + t314 * t327 + t400) * MDP(20) + (g(1) * t383 + g(2) * t378) * MDP(3) + (g(1) * t378 - g(2) * t383) * MDP(2); (-t320 * t289 + (-t303 * t379 + t374 * t426) * t304 - t301 * t295 + (t374 * t471 + (-t374 * t437 + t379 * t431) * t304) * pkin(2) + t392) * MDP(27) + (-t312 * t426 + t301 * t368 + (-t273 * t381 + t312 * t437) * pkin(2) + t399) * MDP(19) + (-t314 * t426 - t303 * t368 + (-t272 * t381 + t314 * t437) * pkin(2) + t400) * MDP(20) + ((-t355 + (-t369 * t440 - t429) * pkin(1)) * t376 + t446) * MDP(13) + (t381 * t393 + t445) * MDP(12) + t387 * MDP(6) + (t393 + t481) * MDP(5) + t388 + (-(t303 * t374 + t379 * t426) * t304 - t301 * t293 + (-t379 * t471 + (t374 * t431 + t379 * t437) * t304) * pkin(2) + t390) * MDP(26); (-t302 * t368 + (-t314 * t454 - t366 * t375) * pkin(2) + t385 * t380 + t389) * MDP(20) + t394 + (t280 * t314 + t300 * t295 + t405 * t379 - t444 * t374 * t361 + (-t380 * t259 + (qJD(4) * t295 - t472) * t375 + (t402 * t374 - t379 * t434) * t304) * pkin(2) + t423) * MDP(27) + qJDD(3) * MDP(11) + MDP(10) * t456 + t376 * t367 * MDP(9) + t314 * t442 - t314 * t441 + (-t380 * t406 - t300 * t368 - t314 * t316 + (-t312 * t454 + t366 * t380) * pkin(2) + t385 * t375 + t398) * MDP(19) + (-g(3) * t381 + t376 * t387) * MDP(12) + (g(3) * t376 + t381 * t387) * MDP(13) + (-t279 * t314 + t300 * t293 + (-t263 - t476) * t379 + t405 * t374 + (-t380 * t260 + (qJD(4) * t293 - t473) * t375 + (-t374 * t434 - t402 * t379) * t304) * pkin(2) + t421) * MDP(26) + (-MDP(7) * t450 + t443 * MDP(8)) * t369 ^ 2; (t299 * t368 + t398 - t419) * MDP(19) + (-t326 * t434 + t389) * MDP(20) + (-g(3) * t459 - t263 * t379 - t293 * t299 + t421) * MDP(26) + (-t295 * t299 - t374 * t484 + t423) * MDP(27) - (MDP(19) * t316 + MDP(26) * t279 - MDP(27) * t280 + t441 - t442) * t314 + (-t368 * MDP(20) + (MDP(26) * t374 + MDP(27) * t379) * (-t304 + t312)) * t298 + (-t395 * MDP(19) * t375 + (-qJD(4) * MDP(19) * t440 - t395 * MDP(20)) * t380) * t479 + t394; t295 * t293 * MDP(21) + (-t293 ^ 2 + t295 ^ 2) * MDP(22) + (t422 + t470) * MDP(23) + (t351 + t469) * MDP(24) + t271 * MDP(25) + (-g(1) * t307 + g(2) * t305 + t280 * t304 - t295 * t298 + t287) * MDP(26) + (g(1) * t308 - g(2) * t306 + t279 * t304 + t293 * t298) * MDP(27) + (t418 * MDP(27) + (-MDP(24) * t314 - MDP(26) * t299 - MDP(27) * t316) * qJD(5)) * t379 + (-qJD(5) * t314 * MDP(23) + (-qJD(5) * t368 - t272) * MDP(24) + (-qJD(5) * t316 + t418) * MDP(26) + (qJD(5) * t299 - t292) * MDP(27)) * t374;];
tau  = t1;
