% Calculate vector of inverse dynamics joint torques for
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:54
% EndTime: 2019-12-31 20:55:59
% DurationCPUTime: 3.72s
% Computational Cost: add. (3761->358), mult. (8859->439), div. (0->0), fcn. (6167->12), ass. (0->172)
t402 = qJ(2) + qJ(3);
t393 = sin(t402);
t394 = cos(t402);
t407 = sin(qJ(1));
t410 = cos(qJ(1));
t442 = g(1) * t410 + g(2) * t407;
t501 = -g(3) * t394 + t393 * t442;
t408 = cos(qJ(3));
t409 = cos(qJ(2));
t469 = qJD(1) * t409;
t456 = t408 * t469;
t405 = sin(qJ(3));
t406 = sin(qJ(2));
t470 = qJD(1) * t406;
t457 = t405 * t470;
t329 = -t456 + t457;
t331 = -t405 * t469 - t408 * t470;
t403 = sin(pkin(8));
t404 = cos(pkin(8));
t303 = -t404 * t329 + t331 * t403;
t399 = qJD(2) + qJD(3);
t500 = t303 * t399;
t465 = qJD(1) * qJD(2);
t454 = t409 * t465;
t464 = qJDD(1) * t406;
t498 = t454 + t464;
t433 = -t329 * t403 - t404 * t331;
t497 = t433 ^ 2;
t493 = pkin(6) + pkin(7);
t395 = t409 * pkin(2);
t484 = pkin(1) + t395;
t325 = t331 * qJ(4);
t359 = t493 * t409;
t349 = qJD(1) * t359;
t332 = t405 * t349;
t358 = t493 * t406;
t347 = qJD(1) * t358;
t474 = -t408 * t347 - t332;
t298 = t325 + t474;
t336 = t408 * t349;
t447 = t347 * t405 - t336;
t482 = qJ(4) * t329;
t430 = t447 + t482;
t462 = pkin(2) * t403 * t405;
t467 = qJD(3) * t408;
t475 = -qJD(3) * t462 - t403 * t430 + (t467 * pkin(2) - t298) * t404;
t483 = qJD(2) * pkin(2);
t339 = -t347 + t483;
t448 = t408 * t339 - t332;
t293 = t325 + t448;
t473 = -t405 * t358 + t408 * t359;
t392 = pkin(8) + t402;
t378 = sin(t392);
t379 = cos(t392);
t438 = t379 * pkin(4) + t378 * qJ(5);
t496 = g(1) * t407 - g(2) * t410;
t439 = t484 * qJDD(1);
t342 = t405 * t409 + t406 * t408;
t311 = t399 * t342;
t463 = qJDD(1) * t409;
t437 = t405 * t464 - t408 * t463;
t296 = qJD(1) * t311 + t437;
t495 = pkin(3) * t296 + qJDD(4);
t295 = qJD(3) * t456 - t399 * t457 + t405 * t463 + t498 * t408;
t397 = qJDD(2) + qJDD(3);
t312 = qJDD(2) * pkin(2) - t493 * t498;
t455 = t406 * t465;
t314 = t493 * (-t455 + t463);
t432 = -t339 * t405 - t336;
t421 = qJD(3) * t432 + t408 * t312 - t405 * t314;
t261 = pkin(3) * t397 - qJ(4) * t295 + qJD(4) * t331 + t421;
t468 = qJD(3) * t405;
t494 = (qJD(3) * t339 + t314) * t408 + t405 * t312 - t349 * t468;
t264 = -qJ(4) * t296 - qJD(4) * t329 + t494;
t251 = t404 * t261 - t403 * t264;
t250 = -t397 * pkin(4) + qJDD(5) - t251;
t357 = t484 * qJD(1);
t313 = pkin(3) * t329 + qJD(4) - t357;
t274 = -pkin(4) * t303 - qJ(5) * t433 + t313;
t419 = -g(3) * t379 - t274 * t433 + t442 * t378 - t250;
t491 = pkin(3) * t331;
t490 = pkin(3) * t393;
t489 = pkin(4) * t378;
t481 = qJ(5) * t379;
t294 = -t432 - t482;
t289 = t404 * t294;
t269 = t293 * t403 + t289;
t480 = t269 * t433;
t479 = t294 * t403;
t478 = t404 * t405;
t252 = t403 * t261 + t404 * t264;
t288 = pkin(3) * t399 + t293;
t268 = t403 * t288 + t289;
t477 = qJD(5) + t475;
t476 = -t298 * t403 + t404 * t430 + (t403 * t408 + t478) * qJD(3) * pkin(2);
t385 = pkin(2) * t408 + pkin(3);
t324 = pkin(2) * t478 + t403 * t385;
t382 = pkin(3) * t394;
t472 = t382 + t395;
t400 = t406 ^ 2;
t471 = -t409 ^ 2 + t400;
t270 = t293 * t404 - t479;
t466 = qJD(5) - t270;
t389 = t406 * t483;
t461 = t397 * qJ(5) + t252;
t460 = t382 + t438;
t458 = qJD(2) * t493;
t453 = pkin(3) * t311 + t389;
t351 = -pkin(2) * t406 - t490;
t452 = t351 - t489;
t450 = t476 * t433;
t271 = t295 * t403 + t404 * t296;
t446 = -t408 * t358 - t359 * t405;
t341 = t405 * t406 - t408 * t409;
t444 = pkin(3) * t341 - t484;
t443 = -t489 - t490;
t376 = pkin(2) * t455;
t440 = t376 + t495;
t267 = t288 * t404 - t479;
t265 = -pkin(4) * t399 + qJD(5) - t267;
t266 = qJ(5) * t399 + t268;
t436 = -t265 * t303 + t266 * t433;
t435 = t267 * t303 + t268 * t433;
t272 = t295 * t404 - t296 * t403;
t323 = t385 * t404 - t462;
t431 = -0.2e1 * pkin(1) * t465 - pkin(6) * qJDD(2);
t429 = -qJ(4) * t342 + t446;
t326 = t376 - t439;
t348 = t406 * t458;
t350 = t409 * t458;
t426 = -t408 * t348 - t405 * t350 - t358 * t467 - t359 * t468;
t279 = pkin(4) * t433 - qJ(5) * t303 - t491;
t425 = -t331 * t329 * MDP(11) + (t329 * t399 + t295) * MDP(13) + (-t437 + (-qJD(1) * t342 - t331) * t399) * MDP(14) + (-t329 ^ 2 + t331 ^ 2) * MDP(12) + t397 * MDP(15);
t411 = qJD(2) ^ 2;
t423 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t411 + t496;
t412 = qJD(1) ^ 2;
t422 = pkin(1) * t412 - pkin(6) * qJDD(1) + t442;
t420 = -t473 * qJD(3) + t348 * t405 - t408 * t350;
t275 = -qJ(4) * t311 - qJD(4) * t341 + t426;
t310 = t399 * t341;
t414 = qJ(4) * t310 - qJD(4) * t342 + t420;
t255 = t275 * t403 - t404 * t414;
t256 = t404 * t275 + t403 * t414;
t300 = -qJ(4) * t341 + t473;
t281 = t300 * t403 - t404 * t429;
t282 = t404 * t300 + t403 * t429;
t418 = t255 * t433 + t256 * t303 - t282 * t271 + t272 * t281 - t442;
t417 = pkin(4) * t271 - qJ(5) * t272 - qJD(5) * t433 + t440;
t416 = -g(3) * t378 + t274 * t303 - t379 * t442 + t461;
t415 = g(3) * t393 - t357 * t329 + t442 * t394 - t494;
t413 = -t357 * t331 + t421 + t501;
t398 = -qJ(4) - t493;
t390 = t399 * qJD(5);
t388 = pkin(2) * t470;
t380 = -pkin(3) * t404 - pkin(4);
t377 = pkin(3) * t403 + qJ(5);
t353 = t410 * t481;
t352 = t407 * t481;
t346 = pkin(1) + t472;
t338 = t410 * t346;
t319 = -pkin(4) - t323;
t318 = qJ(5) + t324;
t308 = -t341 * t403 + t342 * t404;
t307 = t404 * t341 + t342 * t403;
t286 = -t310 * t404 - t311 * t403;
t285 = -t310 * t403 + t404 * t311;
t280 = pkin(4) * t307 - qJ(5) * t308 + t444;
t278 = t279 + t388;
t257 = pkin(4) * t285 - qJ(5) * t286 - qJD(5) * t308 + t453;
t253 = -t439 + t417;
t249 = t390 + t461;
t1 = [qJDD(1) * MDP(1) + t496 * MDP(2) + t442 * MDP(3) + (qJDD(1) * t400 + 0.2e1 * t406 * t454) * MDP(4) + 0.2e1 * (t406 * t463 - t465 * t471) * MDP(5) + (qJDD(2) * t406 + t409 * t411) * MDP(6) + (qJDD(2) * t409 - t406 * t411) * MDP(7) + (t406 * t431 + t409 * t423) * MDP(9) + (-t406 * t423 + t409 * t431) * MDP(10) + (t295 * t342 + t310 * t331) * MDP(11) + (-t295 * t341 - t296 * t342 + t310 * t329 + t311 * t331) * MDP(12) + (-t310 * t399 + t342 * t397) * MDP(13) + (-t311 * t399 - t341 * t397) * MDP(14) + (-t296 * t484 - t357 * t311 + t326 * t341 + t329 * t389 + t394 * t496 + t397 * t446 + t399 * t420) * MDP(16) + (-t295 * t484 + t357 * t310 + t326 * t342 - t331 * t389 - t393 * t496 - t397 * t473 - t399 * t426) * MDP(17) + (-t251 * t308 - t252 * t307 - t267 * t286 - t268 * t285 + t418) * MDP(18) + (t252 * t282 + t268 * t256 - t251 * t281 - t267 * t255 + (t326 + t495) * t444 + t313 * t453 - g(1) * (-t346 * t407 - t398 * t410) - g(2) * (-t398 * t407 + t338)) * MDP(19) + (t253 * t307 - t255 * t399 - t257 * t303 + t271 * t280 + t274 * t285 - t281 * t397 + t379 * t496) * MDP(20) + (-t249 * t307 + t250 * t308 + t265 * t286 - t266 * t285 + t418) * MDP(21) + (-t253 * t308 + t256 * t399 - t257 * t433 - t272 * t280 - t274 * t286 + t282 * t397 + t378 * t496) * MDP(22) + (-g(2) * t338 + t249 * t282 + t250 * t281 + t253 * t280 + t265 * t255 + t266 * t256 + t274 * t257 + (g(1) * t398 - g(2) * t438) * t410 + (-g(1) * (-t346 - t438) + g(2) * t398) * t407) * MDP(23); (t252 * t324 + t251 * t323 - t313 * (t388 - t491) - g(3) * t472 - t442 * t351 + t475 * t268 - t476 * t267) * MDP(19) + (-t271 * t318 + t272 * t319 + t303 * t477 + t436 + t450) * MDP(21) + (-t447 * t399 + (-t329 * t470 + t397 * t408 - t399 * t468) * pkin(2) + t413) * MDP(16) + (t249 * t318 + t250 * t319 - t274 * t278 - g(1) * (t410 * t452 + t353) - g(2) * (t407 * t452 + t352) - g(3) * (t395 + t460) + t477 * t266 + t476 * t265) * MDP(23) + (-g(3) * t409 + t406 * t422) * MDP(9) + (g(3) * t406 + t409 * t422) * MDP(10) + qJDD(2) * MDP(8) + (t474 * t399 + (t331 * t470 - t397 * t405 - t399 * t467) * pkin(2) + t415) * MDP(17) + (t278 * t433 + t318 * t397 + t399 * t477 + t390 + t416) * MDP(22) + (t278 * t303 - t319 * t397 - t399 * t476 + t419) * MDP(20) + (-t271 * t324 - t272 * t323 + t303 * t475 + t435 + t450) * MDP(18) + MDP(6) * t464 + MDP(7) * t463 + t425 + (-MDP(4) * t406 * t409 + MDP(5) * t471) * t412; (-t399 * t432 + t413) * MDP(16) + (t399 * t448 + t415) * MDP(17) + (-t480 - t270 * t303 + (-t271 * t403 - t272 * t404) * pkin(3) + t435) * MDP(18) + (t267 * t269 - t268 * t270 + (t251 * t404 + t252 * t403 + t313 * t331 + t501) * pkin(3)) * MDP(19) + (t269 * t399 + t279 * t303 - t380 * t397 + t419) * MDP(20) + (-t271 * t377 + t272 * t380 + t303 * t466 + t436 - t480) * MDP(21) + (-t270 * t399 + t279 * t433 + t377 * t397 + 0.2e1 * t390 + t416) * MDP(22) + (t249 * t377 + t250 * t380 - t274 * t279 - t265 * t269 - g(1) * (t410 * t443 + t353) - g(2) * (t407 * t443 + t352) - g(3) * t460 + t466 * t266) * MDP(23) + t425; (t267 * t433 - t268 * t303 + t440 - t496) * MDP(19) + (t399 * t433 + t271) * MDP(20) + (-t272 - t500) * MDP(22) + (-t265 * t433 - t266 * t303 + t417 - t496) * MDP(23) - (MDP(19) + MDP(23)) * t439 + (MDP(18) + MDP(21)) * (-t303 ^ 2 - t497); (-t303 * t433 - t397) * MDP(20) + (t272 - t500) * MDP(21) + (-t399 ^ 2 - t497) * MDP(22) + (-t266 * t399 - t419) * MDP(23);];
tau = t1;
