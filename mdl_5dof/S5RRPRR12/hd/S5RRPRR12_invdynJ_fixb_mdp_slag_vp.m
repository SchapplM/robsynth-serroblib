% Calculate vector of inverse dynamics joint torques for
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR12_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:47
% EndTime: 2019-12-31 20:30:53
% DurationCPUTime: 4.28s
% Computational Cost: add. (2250->384), mult. (4791->488), div. (0->0), fcn. (3259->8), ass. (0->182)
t385 = sin(qJ(4));
t389 = cos(qJ(4));
t390 = cos(qJ(2));
t463 = qJD(1) * t390;
t386 = sin(qJ(2));
t464 = qJD(1) * t386;
t319 = -t385 * t463 + t389 * t464;
t453 = qJD(1) * qJD(2);
t445 = t386 * t453;
t451 = qJDD(1) * t390;
t511 = t445 - t451;
t444 = t390 * t453;
t452 = qJDD(1) * t386;
t512 = t444 + t452;
t407 = t385 * t511 + t389 * t512;
t323 = t385 * t386 + t389 * t390;
t410 = t323 * qJD(4);
t274 = -qJD(1) * t410 + t407;
t375 = qJDD(2) - qJDD(4);
t376 = qJD(2) - qJD(4);
t384 = sin(qJ(5));
t388 = cos(qJ(5));
t455 = qJD(5) * t388;
t449 = t388 * t274 - t384 * t375 - t376 * t455;
t456 = qJD(5) * t384;
t265 = -t319 * t456 + t449;
t297 = t319 * t388 - t376 * t384;
t350 = t388 * t375;
t266 = t297 * qJD(5) + t274 * t384 + t350;
t459 = qJD(4) * t389;
t460 = qJD(4) * t385;
t461 = qJD(2) * t390;
t517 = t385 * t461 + t386 * t459 - t390 * t460;
t275 = t517 * qJD(1) + t323 * qJDD(1) - t389 * t445;
t295 = t319 * t384 + t388 * t376;
t272 = qJDD(5) + t275;
t472 = t388 * t272;
t478 = t384 * t272;
t482 = t265 * t384;
t499 = t323 * qJD(1);
t509 = qJD(5) + t499;
t502 = t509 * t297;
t504 = t509 ^ 2;
t518 = t295 * t509;
t519 = -((t266 + t502) * t384 + (-t265 + t518) * t388) * MDP(23) + (t388 * t502 + t482) * MDP(22) + (-t297 * t319 + t388 * t504 + t478) * MDP(24) - (-t295 * t319 + t384 * t504 - t472) * MDP(25) - (t319 * t376 + t275) * MDP(18) - (-t319 ^ 2 + t499 ^ 2) * MDP(16) - t375 * MDP(19) + (MDP(15) * t499 - MDP(26) * t509) * t319;
t366 = pkin(6) * t464;
t510 = -pkin(7) * t464 + qJD(3) + t366;
t466 = t390 * pkin(2) + t386 * qJ(3);
t333 = -pkin(1) - t466;
t392 = -pkin(2) - pkin(3);
t448 = t392 * qJD(2);
t306 = t448 + t510;
t367 = pkin(6) * t463;
t330 = -pkin(7) * t463 + t367;
t379 = qJD(2) * qJ(3);
t320 = t330 + t379;
t280 = t306 * t389 - t320 * t385;
t277 = pkin(4) * t376 - t280;
t503 = t509 * t277;
t467 = t389 * qJ(3) + t385 * t392;
t391 = cos(qJ(1));
t487 = g(2) * t391;
t387 = sin(qJ(1));
t490 = g(1) * t387;
t500 = -t487 + t490;
t488 = g(2) * t387;
t489 = g(1) * t391;
t428 = t488 + t489;
t286 = pkin(4) * t319 + pkin(8) * t499;
t359 = qJ(3) * t463;
t310 = t392 * t464 + t359;
t327 = -pkin(8) + t467;
t351 = pkin(6) * t444;
t361 = pkin(6) * t452;
t443 = qJDD(3) + t351 + t361;
t287 = -pkin(7) * t512 + t392 * qJDD(2) + t443;
t362 = pkin(6) * t451;
t377 = qJDD(2) * qJ(3);
t378 = qJD(2) * qJD(3);
t304 = -pkin(6) * t445 + t362 + t377 + t378;
t288 = pkin(7) * t511 + t304;
t415 = -t389 * t287 + t385 * t288 + t306 * t460 + t320 * t459;
t262 = pkin(4) * t375 + t415;
t473 = t387 * t390;
t476 = t386 * t389;
t311 = t385 * t473 - t387 * t476;
t475 = t386 * t391;
t477 = t385 * t390;
t313 = -t389 * t475 + t391 * t477;
t408 = g(1) * t313 + g(2) * t311 + g(3) * t323;
t404 = -t262 + t408;
t497 = (qJD(5) * t327 - t286 + t310) * t509 + t404;
t496 = (pkin(8) * qJD(5) + t286) * t509 - t404;
t321 = -qJD(1) * pkin(1) - pkin(2) * t463 - qJ(3) * t464;
t486 = pkin(6) * qJDD(2);
t494 = (qJD(1) * t333 + t321) * qJD(2) - t486;
t462 = qJD(2) * t386;
t492 = pkin(6) - pkin(7);
t329 = t492 * t462;
t336 = t492 * t390;
t331 = qJD(2) * t336;
t335 = t492 * t386;
t418 = t335 * t389 - t336 * t385;
t268 = t418 * qJD(4) - t329 * t389 + t331 * t385;
t322 = t390 * pkin(3) - t333;
t324 = t476 - t477;
t279 = pkin(4) * t323 - pkin(8) * t324 + t322;
t292 = t323 * qJD(2) - t410;
t299 = t335 * t385 + t336 * t389;
t303 = pkin(3) * t463 - t321;
t270 = pkin(4) * t499 - pkin(8) * t319 + t303;
t411 = t385 * t287 + t389 * t288 + t306 * t459 - t320 * t460;
t435 = -pkin(8) * t375 + qJD(5) * t270 + t411;
t493 = t262 * t324 - t299 * t272 + t277 * t292 - (qJD(5) * t279 + t268) * t509 - t435 * t323 + t489;
t312 = t323 * t387;
t491 = g(1) * t312;
t382 = qJDD(1) * pkin(1);
t485 = qJDD(2) * pkin(2);
t281 = t306 * t385 + t320 * t389;
t278 = -pkin(8) * t376 + t281;
t263 = t270 * t388 - t278 * t384;
t484 = t263 * t319;
t264 = t270 * t384 + t278 * t388;
t483 = t264 * t319;
t481 = t277 * t324;
t480 = t499 * t376;
t394 = qJD(1) ^ 2;
t474 = t386 * t394;
t422 = -qJ(3) * t385 + t389 * t392;
t471 = t422 * qJD(4) - t330 * t385 + t389 * t510;
t470 = t467 * qJD(4) + t330 * t389 + t385 * t510;
t370 = t386 * qJD(3);
t468 = qJ(3) * t461 + t370;
t380 = t386 ^ 2;
t381 = t390 ^ 2;
t465 = t380 - t381;
t458 = qJD(5) * t278;
t457 = qJD(5) * t319;
t450 = t390 * t474;
t441 = -qJD(2) * pkin(2) + qJD(3);
t432 = t376 ^ 2;
t431 = t389 * t376;
t430 = t386 * t448;
t393 = qJD(2) ^ 2;
t429 = pkin(6) * t393 + t487;
t427 = g(2) * t312 + g(3) * t324;
t426 = pkin(2) * t386 - qJ(3) * t390;
t425 = t279 * t272 + t491;
t424 = pkin(2) * t451 + qJ(3) * t512 + qJD(1) * t370 + t382;
t423 = -t458 - t487;
t332 = t366 + t441;
t334 = t367 + t379;
t419 = t332 * t390 - t334 * t386;
t417 = g(1) * t475 - g(3) * t390 + t386 * t488 - t361;
t414 = -0.2e1 * pkin(1) * t453 - t486;
t413 = t292 * t388 - t324 * t456;
t412 = -qJDD(3) + t417;
t302 = t430 + t468;
t406 = -t429 + 0.2e1 * t382;
t403 = t427 - t435;
t402 = -pkin(8) * t272 + t280 * t509 + t503;
t399 = -t327 * t272 - t471 * t509 - t503;
t285 = pkin(2) * t445 - t424;
t315 = pkin(2) * t462 - t468;
t398 = -qJD(1) * t315 - qJDD(1) * t333 - t285 - t429;
t276 = pkin(3) * t451 + qJD(1) * t430 + t424;
t397 = t303 * t319 - t408 + t415;
t314 = t323 * t391;
t396 = -g(1) * t314 - t303 * t499 + t411 - t427;
t309 = t443 - t485;
t395 = t419 * qJD(2) + t304 * t390 + t309 * t386 - t428;
t353 = g(1) * t473;
t326 = pkin(4) - t422;
t325 = pkin(2) * t464 - t359;
t294 = t314 * t388 - t384 * t387;
t293 = -t314 * t384 - t387 * t388;
t291 = -t389 * t462 + t517;
t269 = t299 * qJD(4) - t329 * t385 - t331 * t389;
t267 = pkin(4) * t291 - pkin(8) * t292 + t302;
t260 = pkin(4) * t275 - pkin(8) * t274 + t276;
t259 = t388 * t260;
t1 = [((-t295 * t388 - t297 * t384) * t292 + (-t482 - t266 * t388 + (t295 * t384 - t297 * t388) * qJD(5)) * t324) * MDP(23) + 0.2e1 * (t386 * t451 - t465 * t453) * MDP(5) + (-g(2) * t314 + t269 * t376 + t275 * t322 + t276 * t323 + t291 * t303 + t302 * t499 - t375 * t418 + t491) * MDP(20) + (-t274 * t323 - t275 * t324 - t291 * t319 - t292 * t499) * MDP(16) + (t395 * pkin(6) + t321 * t315 + (t285 - t500) * t333) * MDP(14) + (qJDD(1) * t380 + 0.2e1 * t386 * t444) * MDP(4) + (qJDD(2) * t386 + t390 * t393) * MDP(6) + (qJDD(2) * t390 - t386 * t393) * MDP(7) + (-t324 * t478 - t266 * t323 - t291 * t295 + (-t292 * t384 - t324 * t455) * t509) * MDP(25) + (t265 * t323 + t291 * t297 + t324 * t472 + t413 * t509) * MDP(24) + (-g(2) * t293 - t264 * t291 - t418 * t265 + t269 * t297 + (-(-qJD(5) * t299 + t267) * t509 - (t260 - t458) * t323 - qJD(5) * t481 - t425) * t384 + t493 * t388) * MDP(28) + (-g(2) * t294 + t259 * t323 + t263 * t291 - t418 * t266 + t269 * t295 + (t267 * t509 + (-t278 * t323 - t299 * t509 + t481) * qJD(5) + t425) * t388 + t493 * t384) * MDP(27) + (t272 * t323 + t291 * t509) * MDP(26) + qJDD(1) * MDP(1) + ((t380 + t381) * qJDD(1) * pkin(6) + t395) * MDP(12) + t500 * MDP(2) + (-g(1) * t311 + g(2) * t313 + t268 * t376 + t274 * t322 + t276 * t324 + t292 * t303 + t299 * t375 + t302 * t319) * MDP(21) + (t291 * t376 + t323 * t375) * MDP(18) + (-t292 * t376 - t324 * t375) * MDP(17) + (t274 * t324 + t292 * t319) * MDP(15) + t428 * MDP(3) + (t414 * t390 + (-t406 - t490) * t386) * MDP(10) + (t265 * t324 * t388 + t413 * t297) * MDP(22) + (t414 * t386 + t406 * t390 + t353) * MDP(9) + (t494 * t386 + t398 * t390 + t353) * MDP(11) + (-t494 * t390 + (t398 + t490) * t386) * MDP(13); (0.2e1 * t485 + (-t321 * t386 + t325 * t390) * qJD(1) + t412) * MDP(11) + (pkin(1) * t474 + t417) * MDP(9) + (-t310 * t319 + t467 * t375 + t471 * t376 + t396) * MDP(21) + (-t310 * t499 - t422 * t375 + t470 * t376 + t397) * MDP(20) + (qJD(4) * t499 - t407 + t480) * MDP(17) + (-t426 * qJDD(1) + ((t334 - t379) * t386 + (-t332 + t441) * t390) * qJD(1)) * MDP(12) - MDP(4) * t450 + (t326 * t265 + t470 * t297 + t497 * t384 + t399 * t388 - t483) * MDP(28) + (t326 * t266 + t470 * t295 + t399 * t384 - t497 * t388 + t484) * MDP(27) + qJDD(2) * MDP(8) + (-t419 * qJD(1) * pkin(6) - t309 * pkin(2) - g(3) * t466 + t304 * qJ(3) + t334 * qJD(3) - t321 * t325 + t428 * t426) * MDP(14) + (t362 + 0.2e1 * t377 + 0.2e1 * t378 + (qJD(1) * t325 - g(3)) * t386 + (qJD(1) * t321 - t428) * t390) * MDP(13) + (g(3) * t386 - t362 + (pkin(1) * t394 + t428) * t390) * MDP(10) + MDP(7) * t451 + t465 * MDP(5) * t394 + MDP(6) * t452 - t519; (-qJDD(2) - t450) * MDP(11) + MDP(12) * t452 + (-t380 * t394 - t393) * MDP(13) + (-qJD(2) * t334 + t321 * t464 + t351 - t412 - t485) * MDP(14) + (-t375 * t389 - t385 * t432 - t464 * t499) * MDP(20) + (-t319 * t464 + t375 * t385 - t389 * t432) * MDP(21) + (-t389 * t266 + (t384 * t431 - t388 * t464) * t509 + (-t376 * t295 - t455 * t509 - t478) * t385) * MDP(27) + (-t389 * t265 + (t384 * t464 + t388 * t431) * t509 + (-t376 * t297 + t456 * t509 - t472) * t385) * MDP(28); (t274 - t480) * MDP(17) + (-t281 * t376 - t397) * MDP(20) + (-t280 * t376 - t396) * MDP(21) + (-pkin(4) * t266 - t281 * t295 + t402 * t384 - t496 * t388 - t484) * MDP(27) + (-pkin(4) * t265 - t281 * t297 + t496 * t384 + t402 * t388 + t483) * MDP(28) + t519; t297 * t295 * MDP(22) + (-t295 ^ 2 + t297 ^ 2) * MDP(23) + (t449 + t518) * MDP(24) + (-t350 + t502) * MDP(25) + t272 * MDP(26) + (-g(1) * t293 + t264 * t509 - t277 * t297 + t259) * MDP(27) + (g(1) * t294 + t263 * t509 + t277 * t295) * MDP(28) + (-MDP(25) * t457 + t423 * MDP(27) + t403 * MDP(28)) * t388 + (-MDP(24) * t457 + (qJD(5) * t376 - t274) * MDP(25) + t403 * MDP(27) + (-t260 - t423) * MDP(28)) * t384;];
tau = t1;
