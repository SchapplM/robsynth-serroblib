% Calculate vector of inverse dynamics joint torques for
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:47:45
% EndTime: 2021-01-15 19:48:05
% DurationCPUTime: 7.20s
% Computational Cost: add. (3358->428), mult. (7908->567), div. (0->0), fcn. (5805->14), ass. (0->181)
t403 = cos(qJ(2));
t463 = cos(pkin(8));
t431 = t463 * t403;
t376 = qJD(1) * t431;
t396 = sin(pkin(8));
t400 = sin(qJ(2));
t446 = qJD(1) * t400;
t348 = t396 * t446 - t376;
t345 = qJD(5) + t348;
t432 = t463 * t400;
t364 = t396 * t403 + t432;
t351 = t364 * qJD(1);
t395 = sin(pkin(9));
t397 = cos(pkin(9));
t336 = t397 * qJD(2) - t351 * t395;
t402 = cos(qJ(5));
t335 = qJD(2) * t395 + t351 * t397;
t399 = sin(qJ(5));
t457 = t335 * t399;
t479 = t336 * t402 - t457;
t481 = t345 * t479;
t398 = -qJ(3) - pkin(6);
t434 = qJD(2) * t398;
t413 = -qJD(3) * t400 + t403 * t434;
t437 = t398 * t400;
t314 = qJDD(2) * pkin(2) + qJD(1) * t413 + qJDD(1) * t437;
t346 = qJD(3) * t403 + t400 * t434;
t373 = t398 * t403;
t323 = qJD(1) * t346 - qJDD(1) * t373;
t272 = t463 * t314 - t396 * t323;
t271 = -qJDD(2) * pkin(3) + qJDD(4) - t272;
t392 = qJ(2) + pkin(8);
t387 = sin(t392);
t401 = sin(qJ(1));
t404 = cos(qJ(1));
t428 = g(1) * t404 + g(2) * t401;
t389 = cos(t392);
t467 = g(3) * t389;
t412 = t387 * t428 - t467;
t480 = t271 - t412;
t419 = -t335 * t402 - t336 * t399;
t477 = t345 * t419;
t365 = t395 * t402 + t397 * t399;
t355 = t365 * qJD(5);
t448 = t365 * t348 + t355;
t475 = g(1) * t401 - g(2) * t404;
t476 = t475 * t387;
t418 = t475 * t389;
t350 = t364 * qJD(2);
t441 = qJDD(1) * t400;
t425 = -qJDD(1) * t431 + t396 * t441;
t320 = qJD(1) * t350 + t425;
t315 = qJDD(5) + t320;
t363 = t395 * t399 - t402 * t397;
t449 = t345 * t363;
t474 = -t315 * t365 + t345 * t449;
t468 = g(3) * t387;
t473 = -t428 * t389 - t468;
t347 = t348 ^ 2;
t472 = pkin(2) * t400;
t471 = pkin(2) * t403;
t470 = pkin(7) * t397;
t466 = g(3) * t403;
t378 = pkin(2) * t396 + qJ(4);
t465 = pkin(7) + t378;
t464 = qJD(2) * pkin(2);
t462 = t479 * t351;
t461 = t419 * t351;
t459 = t320 * t395;
t458 = t320 * t397;
t456 = t348 * t395;
t414 = -t396 * t400 + t431;
t353 = t414 * qJD(2);
t455 = t353 * t395;
t454 = t364 * t395;
t453 = t364 * t397;
t452 = t389 * t401;
t451 = t389 * t404;
t368 = qJD(1) * t373;
t356 = t396 * t368;
t450 = t398 * t401;
t442 = qJD(1) * qJD(2);
t436 = t400 * t442;
t407 = qJDD(1) * t364 - t396 * t436;
t321 = qJD(2) * t376 + t407;
t383 = pkin(1) + t471;
t344 = pkin(2) * t436 - qJDD(1) * t383 + qJDD(3);
t257 = pkin(3) * t320 - qJ(4) * t321 - qJD(4) * t351 + t344;
t273 = t396 * t314 + t463 * t323;
t268 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t273;
t246 = t395 * t257 + t397 * t268;
t439 = t400 * t464;
t283 = pkin(3) * t350 - qJ(4) * t353 - qJD(4) * t364 + t439;
t303 = t346 * t463 + t396 * t413;
t259 = t395 * t283 + t397 * t303;
t371 = -qJD(1) * t383 + qJD(3);
t291 = pkin(3) * t348 - qJ(4) * t351 + t371;
t367 = qJD(1) * t437;
t361 = t367 + t464;
t433 = t463 * t368;
t319 = t396 * t361 - t433;
t311 = qJD(2) * qJ(4) + t319;
t263 = t395 * t291 + t397 * t311;
t301 = pkin(2) * t446 + pkin(3) * t351 + qJ(4) * t348;
t325 = t367 * t463 + t356;
t270 = t395 * t301 + t397 * t325;
t317 = -pkin(3) * t414 - qJ(4) * t364 - t383;
t332 = -t373 * t463 + t396 * t437;
t275 = t395 * t317 + t397 * t332;
t393 = t400 ^ 2;
t447 = -t403 ^ 2 + t393;
t445 = qJD(5) * t399;
t444 = qJD(5) * t402;
t318 = t361 * t463 + t356;
t304 = -qJD(2) * pkin(3) + qJD(4) - t318;
t443 = -qJD(4) + t304;
t440 = qJDD(1) * t403;
t298 = -t397 * qJDD(2) + t321 * t395;
t299 = qJDD(2) * t395 + t321 * t397;
t438 = -t399 * t298 + t402 * t299 + t336 * t444;
t245 = t397 * t257 - t268 * t395;
t241 = pkin(4) * t320 - pkin(7) * t299 + t245;
t244 = -pkin(7) * t298 + t246;
t430 = t402 * t241 - t244 * t399;
t258 = t397 * t283 - t303 * t395;
t262 = t397 * t291 - t311 * t395;
t429 = t402 * t298 + t399 * t299;
t269 = t397 * t301 - t325 * t395;
t274 = t397 * t317 - t332 * t395;
t302 = t346 * t396 - t463 * t413;
t324 = t367 * t396 - t433;
t331 = -t373 * t396 - t398 * t432;
t381 = -pkin(2) * t463 - pkin(3);
t426 = -t363 * t315 - t345 * t448;
t424 = t241 * t399 + t244 * t402;
t423 = -t245 * t397 - t246 * t395;
t250 = pkin(4) * t348 - pkin(7) * t335 + t262;
t254 = pkin(7) * t336 + t263;
t242 = t250 * t402 - t254 * t399;
t243 = t250 * t399 + t254 * t402;
t260 = -pkin(4) * t414 - pkin(7) * t453 + t274;
t264 = -pkin(7) * t454 + t275;
t422 = t260 * t402 - t264 * t399;
t421 = t260 * t399 + t264 * t402;
t420 = -t262 * t395 + t263 * t397;
t417 = -0.2e1 * pkin(1) * t442 - pkin(6) * qJDD(2);
t360 = t465 * t397;
t416 = pkin(4) * t351 + qJD(4) * t395 + qJD(5) * t360 + t348 * t470 + t269;
t359 = t465 * t395;
t415 = pkin(7) * t456 - qJD(4) * t397 + qJD(5) * t359 + t270;
t247 = -t335 * t445 + t438;
t405 = qJD(2) ^ 2;
t410 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t405 + t475;
t406 = qJD(1) ^ 2;
t409 = pkin(1) * t406 - pkin(6) * qJDD(1) + t428;
t408 = t271 * t364 + t304 * t353 - t428;
t248 = -qJD(5) * t419 + t429;
t391 = pkin(9) + qJ(5);
t388 = cos(t391);
t386 = sin(t391);
t382 = t398 * t404;
t372 = -t396 * pkin(3) + qJ(4) * t463;
t370 = -t397 * pkin(4) + t381;
t369 = pkin(3) * t463 + qJ(4) * t396 + pkin(2);
t340 = t386 * t401 + t388 * t451;
t339 = -t386 * t451 + t388 * t401;
t338 = t386 * t404 - t388 * t452;
t337 = t386 * t452 + t388 * t404;
t328 = t369 * t403 + t372 * t400 + pkin(1);
t306 = t363 * t364;
t305 = t365 * t364;
t297 = pkin(4) * t454 + t331;
t284 = -pkin(4) * t456 + t324;
t277 = pkin(4) * t455 + t302;
t276 = -pkin(4) * t336 + t304;
t267 = t353 * t365 + t444 * t453 - t445 * t454;
t266 = -t353 * t363 - t355 * t364;
t252 = pkin(4) * t298 + t271;
t251 = -pkin(7) * t455 + t259;
t249 = pkin(4) * t350 - t353 * t470 + t258;
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t393 + 0.2e1 * t403 * t436) * MDP(4) + 0.2e1 * (t400 * t440 - t442 * t447) * MDP(5) + (qJDD(2) * t400 + t403 * t405) * MDP(6) + (qJDD(2) * t403 - t400 * t405) * MDP(7) + (t400 * t417 + t403 * t410) * MDP(9) + (-t400 * t410 + t403 * t417) * MDP(10) + (-qJDD(2) * t331 - t320 * t383 - t344 * t414 + t350 * t371 + t418 + (t348 * t472 - t302) * qJD(2)) * MDP(11) + (-qJDD(2) * t332 - t321 * t383 + t344 * t364 + t353 * t371 - t476 + (t351 * t472 - t303) * qJD(2)) * MDP(12) + (-t272 * t364 + t273 * t414 + t302 * t351 - t303 * t348 - t318 * t353 - t319 * t350 - t320 * t332 + t321 * t331 - t428) * MDP(13) + (t273 * t332 + t319 * t303 - t272 * t331 - t318 * t302 - t344 * t383 + t371 * t439 - g(1) * (-t383 * t401 - t382) - g(2) * (t383 * t404 - t450)) * MDP(14) + (-t245 * t414 + t258 * t348 + t262 * t350 + t274 * t320 + t331 * t298 - t302 * t336 + t395 * t408 + t397 * t418) * MDP(15) + (t246 * t414 - t259 * t348 - t263 * t350 - t275 * t320 + t331 * t299 + t302 * t335 - t395 * t418 + t397 * t408) * MDP(16) + (-t258 * t335 + t259 * t336 - t274 * t299 - t275 * t298 + t476 + t423 * t364 + (-t262 * t397 - t263 * t395) * t353) * MDP(17) + (t246 * t275 + t263 * t259 + t245 * t274 + t262 * t258 + t271 * t331 + t304 * t302 - g(1) * (-t328 * t401 - t382) - g(2) * (t328 * t404 - t450)) * MDP(18) + (-t247 * t306 - t266 * t419) * MDP(19) + (-t247 * t305 + t248 * t306 + t266 * t479 + t267 * t419) * MDP(20) + (-t247 * t414 + t266 * t345 - t306 * t315 - t350 * t419) * MDP(21) + (t248 * t414 - t267 * t345 - t305 * t315 + t350 * t479) * MDP(22) + (-t315 * t414 + t345 * t350) * MDP(23) + ((t249 * t402 - t251 * t399) * t345 + t422 * t315 - t430 * t414 + t242 * t350 - t277 * t479 + t297 * t248 + t252 * t305 + t276 * t267 - g(1) * t338 - g(2) * t340 + (t243 * t414 - t345 * t421) * qJD(5)) * MDP(24) + (-(t249 * t399 + t251 * t402) * t345 - t421 * t315 + t424 * t414 - t243 * t350 - t277 * t419 + t297 * t247 - t252 * t306 + t276 * t266 - g(1) * t337 - g(2) * t339 + (t242 * t414 - t345 * t422) * qJD(5)) * MDP(25) + t475 * MDP(2) + t428 * MDP(3); MDP(6) * t441 + MDP(7) * t440 + qJDD(2) * MDP(8) + (t400 * t409 - t466) * MDP(9) + (g(3) * t400 + t403 * t409) * MDP(10) + (t324 * qJD(2) - t371 * t351 + (qJDD(2) * t463 - t348 * t446) * pkin(2) + t412 + t272) * MDP(11) + (qJD(2) * t325 + t348 * t371 + (-qJDD(2) * t396 - t351 * t446) * pkin(2) - t273 - t473) * MDP(12) + ((t319 - t324) * t351 + (-t318 + t325) * t348 + (-t320 * t396 - t321 * t463) * pkin(2)) * MDP(13) + (t318 * t324 - t319 * t325 + (t463 * t272 - t466 + t273 * t396 + (-qJD(1) * t371 + t428) * t400) * pkin(2)) * MDP(14) + (-t378 * t459 - t262 * t351 + t298 * t381 + t324 * t336 + (t395 * t443 - t269) * t348 - t480 * t397) * MDP(15) + (-t378 * t458 + t263 * t351 + t299 * t381 - t324 * t335 + (t397 * t443 + t270) * t348 + t480 * t395) * MDP(16) + (t269 * t335 - t270 * t336 + (qJD(4) * t336 - t262 * t348 - t298 * t378 + t246) * t397 + (qJD(4) * t335 - t263 * t348 + t299 * t378 - t245) * t395 + t473) * MDP(17) + (t271 * t381 - t263 * t270 - t262 * t269 - t304 * t324 - g(3) * (pkin(3) * t389 + qJ(4) * t387 + t471) + (-t245 * t395 + t246 * t397) * t378 - t428 * (-t369 * t400 + t372 * t403) + t420 * qJD(4)) * MDP(18) + (t247 * t365 + t419 * t449) * MDP(19) + (-t247 * t363 - t248 * t365 + t419 * t448 - t449 * t479) * MDP(20) + (t461 - t474) * MDP(21) + (t426 - t462) * MDP(22) - t345 * t351 * MDP(23) + ((-t359 * t402 - t360 * t399) * t315 + t370 * t248 + t252 * t363 - t242 * t351 + t284 * t479 + (t399 * t415 - t402 * t416) * t345 + t448 * t276 + t412 * t388) * MDP(24) + (-(-t359 * t399 + t360 * t402) * t315 + t370 * t247 + t252 * t365 + t243 * t351 + t284 * t419 + (t399 * t416 + t402 * t415) * t345 - t449 * t276 - t412 * t386) * MDP(25) + (-MDP(4) * t400 * t403 + MDP(5) * t447) * t406; (0.2e1 * t351 * qJD(2) + t425) * MDP(11) + ((t376 - t348) * qJD(2) + t407) * MDP(12) + (-t351 ^ 2 - t347) * MDP(13) + (t318 * t351 + t319 * t348 + t344 - t475) * MDP(14) + (t336 * t351 - t347 * t395 + t458) * MDP(15) + (-t335 * t351 - t347 * t397 - t459) * MDP(16) + (-t298 * t395 - t299 * t397 + (t335 * t395 + t336 * t397) * t348) * MDP(17) + (-t304 * t351 + t348 * t420 - t423 - t475) * MDP(18) + (t426 + t462) * MDP(24) + (t461 + t474) * MDP(25); (t335 * t348 + t298) * MDP(15) + (t336 * t348 + t299) * MDP(16) + (-t335 ^ 2 - t336 ^ 2) * MDP(17) + (t262 * t335 - t263 * t336 - t364 * t428 + t271 + t467) * MDP(18) + (t248 - t477) * MDP(24) + (t247 + t481) * MDP(25); t419 * t479 * MDP(19) + (t419 ^ 2 - t479 ^ 2) * MDP(20) + (t438 - t481) * MDP(21) + (-t429 - t477) * MDP(22) + t315 * MDP(23) + (-g(1) * t339 + g(2) * t337 + t243 * t345 + t276 * t419 + t386 * t468 + t430) * MDP(24) + (g(1) * t340 - g(2) * t338 + t242 * t345 - t276 * t479 + t388 * t468 - t424) * MDP(25) + (-MDP(21) * t457 + MDP(22) * t419 - MDP(24) * t243 - MDP(25) * t242) * qJD(5);];
tau = t1;
