% Calculate vector of inverse dynamics joint torques for
% S5RRPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:55:04
% EndTime: 2019-12-31 19:55:09
% DurationCPUTime: 3.42s
% Computational Cost: add. (3296->338), mult. (7660->410), div. (0->0), fcn. (5543->12), ass. (0->163)
t389 = sin(pkin(8));
t390 = cos(pkin(8));
t393 = sin(qJ(2));
t395 = cos(qJ(2));
t343 = t389 * t395 + t390 * t393;
t443 = qJD(1) * qJD(2);
t435 = t395 * t443;
t436 = t393 * t443;
t304 = qJDD(1) * t343 - t389 * t436 + t390 * t435;
t425 = t389 * t393 - t390 * t395;
t335 = t425 * qJD(1);
t336 = t343 * qJD(1);
t392 = sin(qJ(4));
t418 = t343 * qJD(2);
t410 = qJD(1) * t418;
t416 = t425 * qJDD(1);
t399 = -t416 - t410;
t471 = cos(qJ(4));
t437 = qJD(4) * t471;
t446 = qJD(4) * t392;
t419 = t471 * t304 - t335 * t437 - t336 * t446 + t392 * t399;
t294 = -t471 * t335 - t336 * t392;
t385 = qJD(2) + qJD(4);
t458 = t294 * t385;
t249 = t419 - t458;
t383 = qJDD(2) + qJDD(4);
t420 = t392 * t335 - t336 * t471;
t403 = qJD(4) * t420 - t392 * t304 + t471 * t399;
t459 = t294 ^ 2;
t460 = t420 * t385;
t484 = t420 ^ 2;
t486 = t294 * t420;
t487 = t249 * MDP(15) + MDP(13) * t486 + t383 * MDP(17) + (t403 - t460) * MDP(16) + (-t459 + t484) * MDP(14);
t381 = t395 * pkin(2);
t374 = t381 + pkin(1);
t352 = -qJD(1) * t374 + qJD(3);
t310 = pkin(3) * t335 + t352;
t262 = -pkin(4) * t294 + qJ(5) * t420 + t310;
t483 = t262 * t294;
t482 = t310 * t294;
t481 = t310 * t420;
t386 = qJ(2) + pkin(8);
t380 = qJ(4) + t386;
t370 = sin(t380);
t371 = cos(t380);
t450 = t371 * pkin(4) + t370 * qJ(5);
t480 = pkin(3) * cos(t386) + t381 + t450;
t379 = t383 * pkin(4);
t474 = qJDD(5) - t379;
t479 = -t262 * t420 + t474;
t269 = -pkin(4) * t420 - qJ(5) * t294;
t464 = qJ(3) + pkin(6);
t360 = t464 * t393;
t348 = qJD(1) * t360;
t361 = t464 * t395;
t349 = qJD(1) * t361;
t453 = t390 * t349;
t308 = t348 * t389 - t453;
t469 = pkin(7) * t335;
t283 = t308 + t469;
t338 = t389 * t349;
t309 = -t390 * t348 - t338;
t468 = pkin(7) * t336;
t284 = t309 - t468;
t372 = pkin(2) * t390 + pkin(3);
t470 = pkin(2) * t389;
t439 = t392 * t470;
t477 = qJD(4) * t439 + t392 * t283 + t284 * t471 - t372 * t437;
t449 = t392 * t372 + t471 * t470;
t375 = t383 * qJ(5);
t378 = t385 * qJD(5);
t476 = t375 + t378;
t317 = pkin(3) * t425 - t374;
t394 = sin(qJ(1));
t396 = cos(qJ(1));
t475 = g(1) * t394 - g(2) * t396;
t432 = qJD(2) * t464;
t332 = qJD(3) * t395 - t393 * t432;
t333 = -qJD(3) * t393 - t395 * t432;
t285 = -t389 * t332 + t390 * t333;
t417 = t425 * qJD(2);
t276 = pkin(7) * t417 + t285;
t286 = t390 * t332 + t389 * t333;
t277 = -pkin(7) * t418 + t286;
t311 = -t390 * t360 - t361 * t389;
t287 = -pkin(7) * t343 + t311;
t312 = -t389 * t360 + t390 * t361;
t288 = -pkin(7) * t425 + t312;
t421 = t287 * t471 - t392 * t288;
t246 = qJD(4) * t421 + t392 * t276 + t277 * t471;
t268 = t392 * t287 + t288 * t471;
t473 = t246 * t385 + t268 * t383 + t370 * t475;
t466 = g(3) * t395;
t465 = t393 * pkin(2);
t463 = qJD(2) * pkin(2);
t342 = -t348 + t463;
t302 = t390 * t342 - t338;
t280 = qJD(2) * pkin(3) + t302 - t468;
t303 = t389 * t342 + t453;
t282 = t303 - t469;
t260 = t392 * t280 + t282 * t471;
t462 = t260 * t385;
t457 = t370 * t394;
t456 = t370 * t396;
t455 = t371 * t394;
t454 = t371 * t396;
t301 = qJDD(2) * pkin(2) + qJD(1) * t333 - qJDD(1) * t360;
t307 = qJD(1) * t332 + qJDD(1) * t361;
t271 = t389 * t301 + t390 * t307;
t452 = qJD(5) - t477;
t451 = qJD(4) * t449 + t283 * t471 - t392 * t284;
t387 = t393 ^ 2;
t447 = -t395 ^ 2 + t387;
t259 = t280 * t471 - t392 * t282;
t444 = qJD(5) - t259;
t442 = qJDD(1) * t393;
t441 = qJDD(1) * t395;
t440 = pkin(2) * t436 + qJDD(3);
t377 = t393 * t463;
t313 = pkin(3) * t336 + qJD(1) * t465;
t434 = -pkin(4) * t370 - pkin(3) * sin(t386) - t465;
t270 = t390 * t301 - t307 * t389;
t258 = qJDD(2) * pkin(3) - pkin(7) * t304 + t270;
t261 = pkin(7) * t399 + t271;
t430 = t392 * t258 + t471 * t261 + t280 * t437 - t282 * t446;
t429 = -t471 * t258 + t392 * t261 + t280 * t446 + t282 * t437;
t428 = g(1) * t396 + g(2) * t394;
t423 = pkin(1) + t480;
t422 = -0.2e1 * pkin(1) * t443 - pkin(6) * qJDD(2);
t415 = t372 * t471 - t439;
t414 = g(1) * t454 + g(2) * t455 + g(3) * t370 - t430;
t412 = -qJDD(1) * t374 + t440;
t411 = t471 * t425;
t409 = g(1) * t456 + g(2) * t457 - g(3) * t371 - t429;
t314 = pkin(3) * t418 + t377;
t397 = qJD(2) ^ 2;
t407 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t397 + t475;
t398 = qJD(1) ^ 2;
t406 = pkin(1) * t398 - pkin(6) * qJDD(1) + t428;
t247 = qJD(4) * t268 - t276 * t471 + t392 * t277;
t405 = g(1) * t455 - g(2) * t454 - t247 * t385 + t383 * t421;
t404 = t259 * t385 + t414;
t306 = t343 * t471 - t392 * t425;
t402 = -t385 * t451 + t409;
t401 = -t409 + t479;
t281 = pkin(3) * t410 + qJDD(1) * t317 + t440;
t245 = -pkin(4) * t403 - qJ(5) * t419 + qJD(5) * t420 + t281;
t384 = -pkin(7) - t464;
t354 = qJ(5) * t454;
t353 = qJ(5) * t455;
t330 = -pkin(4) - t415;
t329 = qJ(5) + t449;
t305 = t343 * t392 + t411;
t273 = qJD(4) * t306 - t392 * t417 + t418 * t471;
t272 = t343 * t446 + t385 * t411 + t392 * t418;
t266 = t305 * pkin(4) - t306 * qJ(5) + t317;
t265 = t269 + t313;
t251 = t385 * qJ(5) + t260;
t250 = -t385 * pkin(4) + t444;
t248 = t273 * pkin(4) + t272 * qJ(5) - t306 * qJD(5) + t314;
t244 = t429 + t474;
t243 = t430 + t476;
t1 = [qJDD(1) * MDP(1) + t475 * MDP(2) + t428 * MDP(3) + (qJDD(1) * t387 + 0.2e1 * t393 * t435) * MDP(4) + 0.2e1 * (t393 * t441 - t443 * t447) * MDP(5) + (qJDD(2) * t393 + t395 * t397) * MDP(6) + (qJDD(2) * t395 - t393 * t397) * MDP(7) + (t393 * t422 + t395 * t407) * MDP(9) + (-t393 * t407 + t395 * t422) * MDP(10) + (-t286 * t335 - t271 * t425 - t285 * t336 - t311 * t304 - t270 * t343 - t312 * t416 + (t302 * t425 - t303 * t343 - t312 * t336) * qJD(2) - t428) * MDP(11) + (t271 * t312 + t303 * t286 + t270 * t311 + t302 * t285 - t412 * t374 + t352 * t377 - g(1) * (-t374 * t394 + t396 * t464) - g(2) * (t374 * t396 + t394 * t464)) * MDP(12) + (t272 * t420 + t306 * t419) * MDP(13) + (-t272 * t294 + t273 * t420 - t305 * t419 + t306 * t403) * MDP(14) + (-t272 * t385 + t306 * t383) * MDP(15) + (-t273 * t385 - t305 * t383) * MDP(16) + (t273 * t310 + t281 * t305 - t294 * t314 - t317 * t403 + t405) * MDP(18) + (-t272 * t310 + t281 * t306 - t314 * t420 + t317 * t419 - t473) * MDP(19) + (t245 * t305 - t248 * t294 + t262 * t273 - t266 * t403 + t405) * MDP(20) + (-t243 * t305 + t244 * t306 + t246 * t294 - t247 * t420 - t250 * t272 - t251 * t273 + t268 * t403 - t419 * t421 - t428) * MDP(21) + (-t245 * t306 + t248 * t420 + t262 * t272 - t266 * t419 + t473) * MDP(22) + (t243 * t268 - t244 * t421 + t245 * t266 + t251 * t246 + t250 * t247 + t262 * t248 + (g(1) * t384 - g(2) * t423) * t396 + (g(1) * t423 + g(2) * t384) * t394) * MDP(23); MDP(6) * t442 + MDP(7) * t441 + qJDD(2) * MDP(8) + (t393 * t406 - t466) * MDP(9) + (g(3) * t393 + t395 * t406) * MDP(10) + ((t303 + t308) * t336 - (-t309 + t302) * t335 + (-t390 * t304 + ((-t435 - t442) * t389 + (-t436 + t441) * t390) * t389) * pkin(2)) * MDP(11) + (-t302 * t308 - t303 * t309 + (-t466 + t270 * t390 + t271 * t389 + (-qJD(1) * t352 + t428) * t393) * pkin(2)) * MDP(12) + (t294 * t313 + t383 * t415 + t402 + t481) * MDP(18) + (t313 * t420 - t449 * t383 + t385 * t477 + t414 - t482) * MDP(19) + (t265 * t294 - t330 * t383 + t402 - t479) * MDP(20) + (t403 * t329 + t330 * t419 + (-t251 - t451) * t420 + (-t250 + t452) * t294) * MDP(21) + (-t265 * t420 + t329 * t383 + t385 * t452 - t414 + t476 + t483) * MDP(22) + (t243 * t329 + t244 * t330 - t262 * t265 - g(1) * (t396 * t434 + t354) - g(2) * (t394 * t434 + t353) - g(3) * t480 + t452 * t251 + t451 * t250) * MDP(23) + (-t393 * t395 * MDP(4) + t447 * MDP(5)) * t398 + t487; (-t335 ^ 2 - t336 ^ 2) * MDP(11) + (t302 * t336 + t303 * t335 + t412 - t475) * MDP(12) + (-t459 - t484) * MDP(21) + (t250 * t420 - t251 * t294 + t245 - t475) * MDP(23) + (MDP(19) - MDP(22)) * (t419 + t458) + (MDP(18) + MDP(20)) * (-t403 - t460); (t409 + t462 + t481) * MDP(18) + (t404 - t482) * MDP(19) + (t269 * t294 + t379 - t401 + t462) * MDP(20) + (-pkin(4) * t419 + qJ(5) * t403 - (t251 - t260) * t420 - (t250 - t444) * t294) * MDP(21) + (-t269 * t420 + 0.2e1 * t375 + 0.2e1 * t378 - t404 + t483) * MDP(22) + (t243 * qJ(5) - t244 * pkin(4) - t262 * t269 - t250 * t260 - g(1) * (-pkin(4) * t456 + t354) - g(2) * (-pkin(4) * t457 + t353) - g(3) * t450 + t444 * t251) * MDP(23) + t487; (-t383 + t486) * MDP(20) + t249 * MDP(21) + (-t385 ^ 2 - t484) * MDP(22) + (-t251 * t385 + t401) * MDP(23);];
tau = t1;
