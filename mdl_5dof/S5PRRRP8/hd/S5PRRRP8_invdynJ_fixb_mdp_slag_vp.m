% Calculate vector of inverse dynamics joint torques for
% S5PRRRP8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:47
% EndTime: 2019-12-05 17:00:55
% DurationCPUTime: 4.80s
% Computational Cost: add. (2337->410), mult. (5392->543), div. (0->0), fcn. (4119->10), ass. (0->174)
t367 = cos(qJ(3));
t424 = qJD(2) * qJD(3);
t411 = t367 * t424;
t364 = sin(qJ(3));
t422 = qJDD(2) * t364;
t492 = -t411 - t422;
t437 = qJD(2) * t367;
t491 = qJD(4) - t437;
t430 = qJD(4) * t364;
t490 = qJD(2) * t430 - qJDD(3);
t363 = sin(qJ(4));
t366 = cos(qJ(4));
t281 = t363 * (qJD(3) * (qJD(4) + t437) + t422) + t490 * t366;
t401 = pkin(3) * t364 - pkin(8) * t367;
t334 = t401 * qJD(3);
t395 = pkin(3) * t367 + pkin(8) * t364 + pkin(2);
t431 = qJD(4) * t363;
t361 = sin(pkin(5));
t441 = qJD(1) * t361;
t365 = sin(qJ(2));
t368 = cos(qJ(2));
t450 = t367 * t368;
t485 = t363 * t450 - t365 * t366;
t489 = -t334 * t366 - t395 * t431 - t485 * t441;
t306 = (t363 * t365 + t366 * t450) * t361;
t429 = qJD(4) * t366;
t488 = -qJD(1) * t306 + t363 * t334 - t395 * t429;
t335 = qJD(2) * pkin(7) + t365 * t441;
t362 = cos(pkin(5));
t454 = t362 * t367;
t487 = qJD(1) * t454 - t364 * t335;
t356 = t367 * qJDD(2);
t486 = -t364 * t424 + t356;
t470 = cos(pkin(9));
t407 = t470 * t365;
t360 = sin(pkin(9));
t458 = t360 * t368;
t316 = t362 * t407 + t458;
t406 = t470 * t368;
t459 = t360 * t365;
t318 = -t362 * t459 + t406;
t457 = t361 * t365;
t320 = t364 * t457 - t454;
t408 = t361 * t470;
t456 = t361 * t367;
t484 = g(3) * t320 - g(2) * (-t316 * t364 - t367 * t408) - g(1) * (-t318 * t364 + t360 * t456);
t295 = -qJD(3) * pkin(3) - t487;
t427 = t366 * qJD(3);
t438 = qJD(2) * t364;
t329 = t363 * t438 - t427;
t434 = qJD(3) * t363;
t331 = t366 * t438 + t434;
t261 = pkin(4) * t329 - qJ(5) * t331 + t295;
t327 = qJDD(4) - t486;
t478 = pkin(8) * t327;
t483 = -t261 * t491 + t478;
t369 = qJD(3) ^ 2;
t425 = qJD(1) * qJD(2);
t413 = t365 * t425;
t455 = t361 * t368;
t397 = -qJDD(1) * t455 + t361 * t413;
t315 = -t362 * t406 + t459;
t474 = g(2) * t315;
t317 = t362 * t458 + t407;
t476 = g(1) * t317;
t400 = t474 + t476;
t482 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t369 + t361 * (-g(3) * t368 + t413) - t397 + t400;
t480 = t331 ^ 2;
t479 = pkin(4) * t327;
t472 = pkin(8) * qJD(4);
t471 = qJD(2) * pkin(2);
t469 = qJ(5) * t327;
t440 = qJD(1) * t364;
t347 = t362 * t440;
t303 = t367 * t335 + t347;
t296 = qJD(3) * pkin(8) + t303;
t420 = t368 * t441;
t304 = -qJD(2) * t395 - t420;
t260 = t296 * t366 + t304 * t363;
t252 = qJ(5) * t491 + t260;
t467 = t252 * t491;
t466 = t260 * t491;
t280 = -qJD(4) * t427 + t490 * t363 + t492 * t366;
t465 = t280 * t363;
t464 = t329 * t331;
t463 = t329 * t491;
t462 = t331 * t491;
t461 = t331 * t366;
t460 = t395 * t366;
t453 = t363 * t367;
t451 = t366 * t367;
t449 = qJDD(1) - g(3);
t428 = qJD(4) * t367;
t433 = qJD(3) * t364;
t448 = qJ(5) * t433 - qJD(5) * t367 + (-t363 * t428 - t364 * t427) * pkin(7) + t488;
t447 = -pkin(4) * t433 + (-t363 * t433 + t366 * t428) * pkin(7) + t489;
t333 = t401 * qJD(2);
t446 = t363 * t333 + t366 * t487;
t398 = pkin(4) * t363 - qJ(5) * t366;
t445 = -qJD(5) * t363 + t491 * t398 - t303;
t443 = pkin(7) * t451 - t363 * t395;
t358 = t364 ^ 2;
t442 = -t367 ^ 2 + t358;
t439 = qJD(2) * t361;
t436 = qJD(3) * t329;
t435 = qJD(3) * t331;
t432 = qJD(3) * t367;
t259 = -t296 * t363 + t304 * t366;
t426 = qJD(5) - t259;
t423 = qJDD(1) * t362;
t418 = t365 * t439;
t417 = t368 * t439;
t416 = t491 * t434;
t415 = t491 * t427;
t414 = t491 * t431;
t409 = t364 * t423;
t308 = qJDD(2) * pkin(7) + (qJDD(1) * t365 + t368 * t425) * t361;
t255 = qJDD(3) * pkin(8) + qJD(3) * t487 + t308 * t367 + t409;
t273 = qJD(2) * t334 - qJDD(2) * t395 + t397;
t404 = t363 * t255 - t366 * t273 + t296 * t429 + t304 * t431;
t403 = t329 * t420;
t402 = t331 * t420;
t399 = pkin(4) * t366 + qJ(5) * t363;
t251 = -pkin(4) * t491 + t426;
t396 = t251 * t366 - t252 * t363;
t394 = pkin(3) + t399;
t393 = pkin(7) + t398;
t321 = t362 * t364 + t365 * t456;
t290 = t321 * t363 + t366 * t455;
t291 = t321 * t366 - t363 * t455;
t389 = -t327 * t363 - t429 * t491;
t388 = t327 * t366 - t414;
t386 = t366 * t255 + t363 * t273 - t296 * t431 + t304 * t429;
t274 = -t315 * t453 - t316 * t366;
t276 = -t317 * t453 - t318 * t366;
t305 = t485 * t361;
t385 = g(1) * t276 + g(2) * t274 + g(3) * t305;
t275 = -t315 * t451 + t316 * t363;
t277 = -t317 * t451 + t318 * t363;
t384 = -g(1) * t277 - g(2) * t275 - g(3) * t306;
t285 = t316 * t367 - t364 * t408;
t287 = t360 * t361 * t364 + t318 * t367;
t382 = g(1) * t287 + g(2) * t285 + g(3) * t321;
t381 = qJD(3) * t347 + t364 * t308 + t335 * t432 - t367 * t423;
t379 = t295 * t491 - t478;
t256 = -qJDD(3) * pkin(3) + t381;
t377 = -t472 * t491 + t484;
t263 = t285 * t363 - t315 * t366;
t265 = t287 * t363 - t317 * t366;
t375 = g(1) * t265 + g(2) * t263 + g(3) * t290 - t404;
t248 = pkin(4) * t281 + qJ(5) * t280 - qJD(5) * t331 + t256;
t374 = -t248 + t377;
t336 = -t420 - t471;
t373 = -pkin(7) * qJDD(3) + (t336 + t420 - t471) * qJD(3);
t264 = t285 * t366 + t315 * t363;
t266 = t287 * t366 + t317 * t363;
t372 = -g(1) * t266 - g(2) * t264 - g(3) * t291 + t386;
t371 = t261 * t331 + qJDD(5) - t375;
t370 = qJD(2) ^ 2;
t313 = t393 * t364;
t298 = t460 + (pkin(7) * t363 + pkin(4)) * t367;
t297 = -qJ(5) * t367 + t443;
t289 = qJD(3) * t321 + t364 * t417;
t288 = -qJD(3) * t320 + t367 * t417;
t283 = pkin(4) * t331 + qJ(5) * t329;
t272 = (qJD(4) * t399 - qJD(5) * t366) * t364 + t393 * t432;
t271 = -pkin(4) * t438 - t333 * t366 + t363 * t487;
t270 = qJ(5) * t438 + t446;
t258 = -t280 + t463;
t250 = -qJD(4) * t290 + t288 * t366 + t363 * t418;
t249 = qJD(4) * t291 + t288 * t363 - t366 * t418;
t247 = qJDD(5) + t404 - t479;
t246 = qJD(5) * t491 + t386 + t469;
t1 = [t449 * MDP(1) + (-qJD(3) * t289 - qJDD(3) * t320) * MDP(10) + (-qJD(3) * t288 - qJDD(3) * t321) * MDP(11) + (t249 * t331 - t250 * t329 - t280 * t290 - t281 * t291) * MDP(20) + (t246 * t291 + t247 * t290 + t248 * t320 + t249 * t251 + t250 * t252 + t261 * t289 - g(3)) * MDP(22) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t367 + MDP(11) * t364 - MDP(3)) * t370) * t365 + (t486 * MDP(10) + t492 * MDP(11) + qJDD(2) * MDP(3) - t370 * MDP(4)) * t368) * t361 + (MDP(17) + MDP(19)) * (-t249 * t491 + t320 * t281 + t289 * t329 - t290 * t327) + (MDP(18) - MDP(21)) * (-t250 * t491 - t280 * t320 + t289 * t331 - t291 * t327); qJDD(2) * MDP(2) + (t449 * t455 + t400) * MDP(3) + (g(1) * t318 + g(2) * t316 - t449 * t457) * MDP(4) + (qJDD(2) * t358 + 0.2e1 * t364 * t411) * MDP(5) + 0.2e1 * (t356 * t364 - t424 * t442) * MDP(6) + (qJDD(3) * t364 + t367 * t369) * MDP(7) + (qJDD(3) * t367 - t364 * t369) * MDP(8) + (t373 * t364 + t367 * t482) * MDP(10) + (-t364 * t482 + t373 * t367) * MDP(11) + (-t280 * t364 * t366 + (-t363 * t430 + t427 * t367) * t331) * MDP(12) + ((-t329 * t366 - t331 * t363) * t432 + (t465 - t281 * t366 + (t329 * t363 - t461) * qJD(4)) * t364) * MDP(13) + ((t280 + t415) * t367 + (t388 + t435) * t364) * MDP(14) + ((t281 - t416) * t367 + (t389 - t436) * t364) * MDP(15) + (-t327 * t367 + t433 * t491) * MDP(16) + (-t327 * t460 - t489 * t491 + (t295 * t434 + (t389 + t436) * pkin(7) + t404) * t367 + (-t403 + t295 * t429 + t259 * qJD(3) + t256 * t363 + (t281 + t416) * pkin(7)) * t364 + t384) * MDP(17) + (-t443 * t327 - t488 * t491 + (t295 * t427 + (t414 + t435) * pkin(7) + t386) * t367 + (-t402 - t295 * t431 - t260 * qJD(3) + t256 * t366 + (-t280 + t415) * pkin(7)) * t364 + t385) * MDP(18) + (t272 * t329 + t281 * t313 - t298 * t327 + (t261 * t434 + t247) * t367 - t447 * t491 + (-qJD(3) * t251 + t248 * t363 + t261 * t429 - t403) * t364 + t384) * MDP(19) + (-t280 * t298 - t281 * t297 + t447 * t331 - t448 * t329 + t396 * t432 + (-g(3) * t455 - t246 * t363 + t247 * t366 + (-t251 * t363 - t252 * t366) * qJD(4) + t400) * t364) * MDP(20) + (-t272 * t331 + t280 * t313 + t297 * t327 + (-t261 * t427 - t246) * t367 + t448 * t491 + (qJD(3) * t252 - t248 * t366 + t261 * t431 + t402) * t364 - t385) * MDP(21) + (t246 * t297 + t248 * t313 + t261 * t272 + t247 * t298 - g(1) * (pkin(4) * t277 + pkin(7) * t318 + qJ(5) * t276) - g(2) * (pkin(4) * t275 + pkin(7) * t316 + qJ(5) * t274) - g(3) * (pkin(4) * t306 + qJ(5) * t305) + t395 * t476 + t395 * t474 + t448 * t252 + t447 * t251 + (-g(3) * pkin(7) * t365 + (-g(3) * t395 - t261 * t440) * t368) * t361) * MDP(22); MDP(7) * t422 + MDP(8) * t356 + qJDD(3) * MDP(9) + (qJD(3) * t303 - t336 * t438 - t381 + t484) * MDP(10) + (-t409 + (-qJD(2) * t336 - t308) * t367 + t382) * MDP(11) + (t461 * t491 - t465) * MDP(12) + ((-t280 - t463) * t366 + (-t281 - t462) * t363) * MDP(13) + ((-t331 * t364 - t451 * t491) * qJD(2) - t389) * MDP(14) + ((t329 * t364 + t453 * t491) * qJD(2) + t388) * MDP(15) - t491 * MDP(16) * t438 + (-t259 * t438 - pkin(3) * t281 - t303 * t329 + (t487 * t491 + t379) * t363 + (-t256 - (t333 + t472) * t491 + t484) * t366) * MDP(17) + (pkin(3) * t280 + t446 * t491 + t260 * t438 - t303 * t331 + t379 * t366 + (t256 - t377) * t363) * MDP(18) + (t251 * t438 + t271 * t491 - t281 * t394 + t445 * t329 - t363 * t483 + t374 * t366) * MDP(19) + (t270 * t329 - t271 * t331 + (t246 + t491 * t251 + (qJD(4) * t331 - t281) * pkin(8)) * t366 + (t247 - t467 + (qJD(4) * t329 - t280) * pkin(8)) * t363 - t382) * MDP(20) + (-t252 * t438 - t270 * t491 - t280 * t394 - t445 * t331 + t374 * t363 + t366 * t483) * MDP(21) + (-t251 * t271 - t252 * t270 + t445 * t261 + (qJD(4) * t396 + t246 * t366 + t247 * t363 - t382) * pkin(8) + (-t248 + t484) * t394) * MDP(22) + (-MDP(5) * t364 * t367 + MDP(6) * t442) * t370; MDP(12) * t464 + (-t329 ^ 2 + t480) * MDP(13) + t258 * MDP(14) + (-t281 + t462) * MDP(15) + t327 * MDP(16) + (-t295 * t331 + t375 + t466) * MDP(17) + (t259 * t491 + t295 * t329 - t372) * MDP(18) + (-t283 * t329 - t371 + t466 + 0.2e1 * t479) * MDP(19) + (pkin(4) * t280 - qJ(5) * t281 + (t252 - t260) * t331 + (t251 - t426) * t329) * MDP(20) + (0.2e1 * t469 - t261 * t329 + t283 * t331 - (-0.2e1 * qJD(5) + t259) * t491 + t372) * MDP(21) + (t246 * qJ(5) - t247 * pkin(4) - t261 * t283 - t251 * t260 - g(1) * (-pkin(4) * t265 + qJ(5) * t266) - g(2) * (-pkin(4) * t263 + qJ(5) * t264) - g(3) * (-pkin(4) * t290 + qJ(5) * t291) + t426 * t252) * MDP(22); (-t327 + t464) * MDP(19) + t258 * MDP(20) + (-t491 ^ 2 - t480) * MDP(21) + (t371 - t467 - t479) * MDP(22);];
tau = t1;
