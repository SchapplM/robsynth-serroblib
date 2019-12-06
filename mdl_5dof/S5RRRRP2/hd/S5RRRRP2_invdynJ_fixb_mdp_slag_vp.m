% Calculate vector of inverse dynamics joint torques for
% S5RRRRP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:48:21
% EndTime: 2019-12-05 18:48:27
% DurationCPUTime: 2.45s
% Computational Cost: add. (2551->298), mult. (3847->360), div. (0->0), fcn. (2445->12), ass. (0->166)
t405 = cos(qJ(3));
t396 = qJD(1) + qJD(2);
t403 = sin(qJ(2));
t499 = pkin(1) * t403;
t461 = qJD(1) * t499;
t346 = pkin(7) * t396 + t461;
t447 = pkin(8) * t396 + t346;
t307 = t447 * t405;
t401 = sin(qJ(4));
t300 = t401 * t307;
t402 = sin(qJ(3));
t306 = t447 * t402;
t303 = qJD(3) * pkin(3) - t306;
t500 = cos(qJ(4));
t445 = t500 * t303 - t300;
t336 = t401 * t405 + t500 * t402;
t316 = t336 * t396;
t489 = t316 * qJ(5);
t509 = t489 - t445;
t399 = qJ(3) + qJ(4);
t385 = sin(t399);
t387 = cos(t399);
t400 = qJ(1) + qJ(2);
t386 = sin(t400);
t388 = cos(t400);
t472 = -g(2) * t386 + g(3) * t388;
t508 = -g(1) * t387 + t472 * t385;
t507 = g(2) * t388 + g(3) * t386;
t449 = qJD(4) * t500;
t453 = t500 * t405;
t506 = -qJD(3) * t453 - t405 * t449;
t483 = t401 * t402;
t425 = t453 - t483;
t390 = t405 * pkin(3);
t494 = pkin(2) + t390;
t501 = -pkin(7) - pkin(8);
t455 = qJD(3) * t501;
t343 = t402 * t455;
t344 = t405 * t455;
t406 = cos(qJ(2));
t468 = qJD(1) * t406;
t460 = pkin(1) * t468;
t360 = t501 * t402;
t389 = t405 * pkin(8);
t361 = pkin(7) * t405 + t389;
t473 = t401 * t360 + t500 * t361;
t505 = -t473 * qJD(4) + t336 * t460 - t401 * t343 + t500 * t344;
t463 = qJD(4) * t401;
t504 = -t500 * t343 - t401 * t344 - t360 * t449 + t361 * t463 + t425 * t460;
t377 = pkin(7) + t499;
t493 = -pkin(8) - t377;
t332 = t493 * t402;
t333 = t377 * t405 + t389;
t474 = t401 * t332 + t500 * t333;
t465 = qJD(3) * t402;
t383 = pkin(3) * t465;
t503 = t383 - t461;
t395 = qJD(3) + qJD(4);
t502 = t316 ^ 2;
t498 = pkin(1) * t406;
t393 = qJDD(1) + qJDD(2);
t497 = pkin(2) * t393;
t496 = pkin(2) * t396;
t492 = qJ(5) * t336;
t458 = t396 * t483;
t314 = -t396 * t453 + t458;
t318 = -t396 * t494 - t460;
t287 = pkin(4) * t314 + qJD(5) + t318;
t491 = t287 * t316;
t490 = t314 * qJ(5);
t488 = t386 * t387;
t487 = t387 * t388;
t486 = t393 * t402;
t485 = t393 * t405;
t484 = t396 * t402;
t481 = t402 * t405;
t294 = t395 * t336;
t477 = -t294 * qJ(5) + qJD(5) * t425;
t480 = t477 - t504;
t430 = t395 * t483;
t293 = t430 + t506;
t428 = t293 * qJ(5) - t336 * qJD(5);
t479 = t428 + t505;
t268 = pkin(4) * t395 - t509;
t478 = t268 + t509;
t476 = -t500 * t306 - t300;
t467 = qJD(2) * t403;
t384 = pkin(1) * t467;
t470 = -qJD(1) * t384 + qJDD(1) * t498;
t319 = -t470 - t497;
t347 = -t460 - t496;
t464 = qJD(3) * t405;
t475 = t319 * t402 + t347 * t464;
t471 = pkin(4) * t387 + t390;
t397 = t402 ^ 2;
t469 = -t405 ^ 2 + t397;
t466 = qJD(2) * t406;
t462 = qJDD(1) * t403;
t459 = pkin(1) * t466;
t457 = t347 * t465 + t507 * t405;
t379 = -pkin(2) - t498;
t302 = t500 * t307;
t452 = t396 * t467;
t451 = t396 * t465;
t450 = t396 * t464;
t448 = pkin(4) * t294 + t383;
t446 = qJD(3) * t493;
t444 = t306 * t401 - t302;
t443 = t500 * t332 - t333 * t401;
t442 = t500 * t360 - t361 * t401;
t342 = pkin(2) + t471;
t394 = -qJ(5) + t501;
t441 = -t342 * t388 + t386 * t394;
t440 = t396 * t461;
t290 = pkin(3) * t451 - t393 * t494 - t470;
t439 = g(2) * t487 + g(3) * t488 - t290 * t425 + t318 * t294;
t438 = -t336 * t393 + t506 * t396;
t437 = -t470 - t507;
t436 = -pkin(4) * t425 - t494;
t431 = t425 * t393;
t429 = -t342 * t386 - t388 * t394;
t273 = t396 * t430 + t438;
t392 = qJDD(3) + qJDD(4);
t320 = pkin(7) * t393 + (qJD(1) * t466 + t462) * pkin(1);
t277 = -t346 * t464 + qJDD(3) * pkin(3) - t320 * t402 + (-t450 - t486) * pkin(8);
t278 = -t346 * t465 + t320 * t405 + (-t451 + t485) * pkin(8);
t426 = -t401 * t303 - t302;
t414 = qJD(4) * t426 + t500 * t277 - t401 * t278;
t255 = t392 * pkin(4) + t273 * qJ(5) - t316 * qJD(5) + t414;
t274 = t294 * t396 - t431;
t411 = t401 * t277 + t500 * t278 + t303 * t449 - t307 * t463;
t256 = -t274 * qJ(5) - t314 * qJD(5) + t411;
t270 = -t426 - t490;
t427 = -t255 * t336 + t256 * t425 + t268 * t293 - t270 * t294 - t472;
t304 = t402 * t446 + t405 * t459;
t305 = -t402 * t459 + t405 * t446;
t424 = t500 * t304 + t401 * t305 + t332 * t449 - t333 * t463;
t313 = t314 ^ 2;
t422 = t316 * t314 * MDP(14) + (-t438 + (t314 - t458) * t395) * MDP(16) + t431 * MDP(17) + (-t313 + t502) * MDP(15) + t392 * MDP(18);
t421 = -t347 * t396 - t320 + t472;
t420 = t290 * t336 - t318 * t293 - t385 * t507;
t408 = qJD(3) ^ 2;
t419 = (-t273 * t425 - t274 * t336 + t293 * t314 - t294 * t316) * MDP(15) + (-t273 * t336 - t293 * t316) * MDP(14) + (-t293 * t395 + t336 * t392) * MDP(16) + (-t294 * t395 + t392 * t425) * MDP(17) + 0.2e1 * (-t469 * t396 * qJD(3) + t393 * t481) * MDP(8) + (t393 * t397 + 0.2e1 * t402 * t450) * MDP(7) + (qJDD(3) * t405 - t402 * t408) * MDP(10) + (qJDD(3) * t402 + t405 * t408) * MDP(9) + t393 * MDP(4);
t418 = -pkin(7) * t408 + t440 + t497;
t417 = -pkin(1) * t452 - t377 * t408 - t379 * t393;
t416 = -pkin(7) * qJDD(3) + (t460 - t496) * qJD(3);
t415 = -qJDD(3) * t377 + (t379 * t396 - t459) * qJD(3);
t265 = pkin(4) * t274 + qJDD(5) + t290;
t413 = -t474 * qJD(4) - t401 * t304 + t500 * t305;
t410 = g(1) * t385 - g(2) * t488 + g(3) * t487 + t318 * t314 - t411;
t409 = -t318 * t316 + t414 + t508;
t407 = cos(qJ(1));
t404 = sin(qJ(1));
t378 = t500 * pkin(3) + pkin(4);
t355 = t379 - t390;
t345 = t384 + t383;
t331 = t425 * qJ(5);
t289 = t331 + t473;
t288 = t442 - t492;
t282 = t331 + t474;
t281 = t443 - t492;
t272 = -t489 + t476;
t271 = t444 + t490;
t259 = t413 + t428;
t258 = t424 + t477;
t1 = [(g(2) * t407 + g(3) * t404) * MDP(2) + (-g(2) * t404 + g(3) * t407) * MDP(3) + (t256 * t282 + t270 * t258 + t255 * t281 + t268 * t259 + t265 * (t436 - t498) + t287 * (t384 + t448) - g(2) * (-pkin(1) * t407 + t441) - g(3) * (-pkin(1) * t404 + t429)) * MDP(22) + qJDD(1) * MDP(1) + (t415 * t405 + (-t417 - t507) * t402 + t475) * MDP(13) + (t415 * t402 + (-t319 + t417) * t405 + t457) * MDP(12) + t419 + (((-qJDD(1) - t393) * t403 + (-qJD(1) - t396) * t466) * pkin(1) + t472) * MDP(6) + (t355 * t274 + t345 * t314 + t392 * t443 + t395 * t413 + t439) * MDP(19) + (-t355 * t273 + t345 * t316 - t474 * t392 - t424 * t395 + t420) * MDP(20) + (-t258 * t314 - t259 * t316 + t273 * t281 - t274 * t282 + t427) * MDP(21) + ((t393 * t406 - t452) * pkin(1) - t437) * MDP(5); (-t437 + t440) * MDP(5) + (t416 * t405 + (-t418 - t507) * t402 + t475) * MDP(13) + (t416 * t402 + (-t319 + t418) * t405 + t457) * MDP(12) + t419 + (-t274 * t494 + t503 * t314 + t392 * t442 + t505 * t395 + t439) * MDP(19) + ((-t462 + (-qJD(2) + t396) * t468) * pkin(1) + t472) * MDP(6) + (t256 * t289 + t255 * t288 + t265 * t436 - g(2) * t441 - g(3) * t429 + (t448 - t461) * t287 + t480 * t270 + t479 * t268) * MDP(22) + (t273 * t494 + t503 * t316 - t473 * t392 + t504 * t395 + t420) * MDP(20) + (t273 * t288 - t274 * t289 - t480 * t314 - t479 * t316 + t427) * MDP(21); MDP(9) * t486 + MDP(10) * t485 + qJDD(3) * MDP(11) + (-g(1) * t405 + t402 * t421) * MDP(12) + (g(1) * t402 + t405 * t421) * MDP(13) + (-t444 * t395 + (-t314 * t484 + t500 * t392 - t395 * t463) * pkin(3) + t409) * MDP(19) + (t476 * t395 + (-t316 * t484 - t401 * t392 - t395 * t449) * pkin(3) + t410) * MDP(20) + (t378 * t273 + (t270 + t271) * t316 + (-t268 + t272) * t314 + (-t274 * t401 + (-t500 * t314 + t316 * t401) * qJD(4)) * pkin(3)) * MDP(21) + (t255 * t378 - t270 * t272 - t268 * t271 - pkin(4) * t491 - g(1) * t471 - t472 * (-pkin(3) * t402 - pkin(4) * t385) + (-t287 * t484 + t256 * t401 + (-t268 * t401 + t500 * t270) * qJD(4)) * pkin(3)) * MDP(22) + t422 + (-MDP(7) * t481 + t469 * MDP(8)) * t396 ^ 2; (-t395 * t426 + t409) * MDP(19) + (t395 * t445 + t410) * MDP(20) + (pkin(4) * t273 - t478 * t314) * MDP(21) + (t478 * t270 + (t255 - t491 + t508) * pkin(4)) * MDP(22) + t422; (-t313 - t502) * MDP(21) + (t268 * t316 + t270 * t314 + t265 - t507) * MDP(22);];
tau = t1;
