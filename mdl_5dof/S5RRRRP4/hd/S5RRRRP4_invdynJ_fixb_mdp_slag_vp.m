% Calculate vector of inverse dynamics joint torques for
% S5RRRRP4
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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:21
% EndTime: 2019-12-31 21:51:25
% DurationCPUTime: 2.69s
% Computational Cost: add. (3264->348), mult. (4805->409), div. (0->0), fcn. (3006->12), ass. (0->180)
t408 = qJ(1) + qJ(2);
t398 = cos(t408);
t411 = sin(qJ(2));
t477 = qJD(2) * t411;
t392 = pkin(1) * t477;
t414 = cos(qJ(2));
t508 = pkin(1) * t414;
t480 = -qJD(1) * t392 + qJDD(1) * t508;
t521 = -g(2) * t398 + t480;
t402 = qJDD(1) + qJDD(2);
t507 = pkin(2) * t402;
t520 = -t507 - t521;
t413 = cos(qJ(3));
t510 = cos(qJ(4));
t457 = qJD(4) * t510;
t461 = t510 * t413;
t519 = -qJD(3) * t461 - t413 * t457;
t409 = sin(qJ(4));
t410 = sin(qJ(3));
t490 = t409 * t410;
t436 = t461 - t490;
t339 = t409 * t413 + t510 * t410;
t404 = qJD(1) + qJD(2);
t329 = t339 * t404;
t403 = qJD(3) + qJD(4);
t441 = t403 * t490;
t450 = -t339 * t402 + t519 * t404;
t278 = t404 * t441 + t450;
t313 = t403 * t339;
t443 = t436 * t402;
t279 = t313 * t404 - t443;
t505 = pkin(3) * t413;
t386 = pkin(2) + t505;
t475 = qJD(3) * t410;
t459 = t404 * t475;
t305 = pkin(3) * t459 - t386 * t402 - t480;
t264 = pkin(4) * t279 + qJ(5) * t278 - qJD(5) * t329 + t305;
t467 = t404 * t490;
t327 = -t404 * t461 + t467;
t478 = qJD(1) * t414;
t469 = pkin(1) * t478;
t331 = -t386 * t404 - t469;
t285 = pkin(4) * t327 - qJ(5) * t329 + t331;
t518 = -t264 * t436 + t285 * t313;
t312 = t441 + t519;
t517 = -t264 * t339 + t285 * t312;
t516 = -t305 * t436 + t331 * t313;
t515 = t305 * t339 - t331 * t312;
t407 = qJ(3) + qJ(4);
t395 = sin(t407);
t397 = cos(t407);
t482 = t397 * pkin(4) + t395 * qJ(5);
t396 = sin(t408);
t481 = g(1) * t398 + g(2) * t396;
t401 = qJDD(3) + qJDD(4);
t389 = t401 * qJ(5);
t393 = t403 * qJD(5);
t514 = t389 + t393;
t394 = t401 * pkin(4);
t513 = qJDD(5) - t394;
t416 = -pkin(8) - pkin(7);
t366 = t416 * t410;
t399 = t413 * pkin(8);
t367 = pkin(7) * t413 + t399;
t318 = t409 * t366 + t510 * t367;
t498 = t395 * t398;
t499 = t395 * t396;
t484 = g(1) * t499 - g(2) * t498;
t464 = qJD(3) * t416;
t342 = t410 * t464;
t343 = t413 * t464;
t437 = t510 * t366 - t409 * t367;
t487 = t437 * qJD(4) + t510 * t342 + t409 * t343 - t436 * t469;
t512 = t318 * t401 + t487 * t403 + t484;
t511 = t329 ^ 2;
t509 = pkin(1) * t411;
t506 = pkin(2) * t404;
t380 = g(1) * t396;
t384 = pkin(7) + t509;
t503 = -pkin(8) - t384;
t470 = qJD(1) * t509;
t345 = pkin(7) * t404 + t470;
t456 = pkin(8) * t404 + t345;
t321 = t456 * t410;
t316 = qJD(3) * pkin(3) - t321;
t322 = t456 * t413;
t491 = t409 * t322;
t291 = t510 * t316 - t491;
t502 = t291 * t403;
t463 = t510 * t322;
t292 = t409 * t316 + t463;
t501 = t292 * t403;
t500 = t327 * t329;
t497 = t396 * t397;
t496 = t397 * t398;
t495 = t398 * t416;
t494 = t402 * t410;
t493 = t402 * t413;
t492 = t404 * t410;
t488 = t410 * t413;
t486 = t318 * qJD(4) - t339 * t469 + t409 * t342 - t510 * t343;
t346 = -t469 - t506;
t485 = t346 * t475 + t413 * t380;
t297 = -t510 * t321 - t491;
t483 = pkin(3) * t457 + qJD(5) - t297;
t405 = t410 ^ 2;
t479 = -t413 ^ 2 + t405;
t476 = qJD(2) * t414;
t474 = qJD(3) * t413;
t473 = qJD(4) * t409;
t472 = qJD(5) - t291;
t471 = qJDD(1) * t411;
t468 = pkin(1) * t476;
t391 = pkin(3) * t475;
t466 = t346 * t474 + t520 * t410;
t465 = pkin(4) * t496 + qJ(5) * t498 + t398 * t386;
t460 = t404 * t477;
t458 = t404 * t474;
t454 = qJD(3) * t503;
t453 = t404 * t470;
t333 = pkin(7) * t402 + (qJD(1) * t476 + t471) * pkin(1);
t288 = -t345 * t474 + qJDD(3) * pkin(3) - t333 * t410 + (-t458 - t494) * pkin(8);
t293 = -t345 * t475 + t333 * t413 + (-t459 + t493) * pkin(8);
t452 = t409 * t288 + t510 * t293 + t316 * t457 - t322 * t473;
t451 = -t510 * t288 + t409 * t293 + t316 * t473 + t322 * t457;
t447 = g(1) * t497 - g(2) * t496;
t277 = pkin(4) * t313 + qJ(5) * t312 - qJD(5) * t339 + t391;
t446 = -t277 + t470;
t296 = -t409 * t321 + t463;
t445 = pkin(3) * t473 - t296;
t444 = -pkin(3) * t410 - pkin(4) * t395;
t303 = pkin(4) * t329 + qJ(5) * t327;
t440 = t380 + t521;
t263 = t452 + t514;
t265 = t451 + t513;
t280 = -t403 * pkin(4) + t472;
t284 = t403 * qJ(5) + t292;
t439 = t263 * t436 + t265 * t339 - t280 * t312 - t284 * t313 - t481;
t336 = t503 * t410;
t337 = t384 * t413 + t399;
t438 = t510 * t336 - t409 * t337;
t307 = t409 * t336 + t510 * t337;
t435 = t391 - t470;
t434 = g(1) * t496 + g(2) * t497 + g(3) * t395 - t452;
t319 = t410 * t454 + t413 * t468;
t320 = -t410 * t468 + t413 * t454;
t267 = t438 * qJD(4) + t510 * t319 + t409 * t320;
t433 = t267 * t403 + t307 * t401 + t484;
t276 = -t450 + (t327 - t467) * t403;
t432 = MDP(14) * t500 + t276 * MDP(16) + t443 * MDP(17) + (-t327 ^ 2 + t511) * MDP(15) + t401 * MDP(18);
t308 = -pkin(4) * t436 - qJ(5) * t339 - t386;
t431 = -t346 * t404 - t333 + t481;
t430 = g(1) * t498 + g(2) * t499 - g(3) * t397 - t451;
t417 = qJD(3) ^ 2;
t429 = (-t278 * t436 - t279 * t339 + t312 * t327 - t313 * t329) * MDP(15) + (-t278 * t339 - t312 * t329) * MDP(14) + (-t312 * t403 + t339 * t401) * MDP(16) + (-t313 * t403 + t401 * t436) * MDP(17) + 0.2e1 * (-t479 * t404 * qJD(3) + t402 * t488) * MDP(8) + (t402 * t405 + 0.2e1 * t410 * t458) * MDP(7) + (qJDD(3) * t413 - t410 * t417) * MDP(10) + (qJDD(3) * t410 + t413 * t417) * MDP(9) + t402 * MDP(4);
t428 = pkin(7) * t417 - t453 - t507;
t268 = t307 * qJD(4) + t409 * t319 - t510 * t320;
t427 = -t268 * t403 + t401 * t438 + t447;
t387 = -pkin(2) - t508;
t426 = pkin(1) * t460 + t384 * t417 + t387 * t402;
t425 = -t285 * t327 - t434;
t424 = t331 * t327 + t434;
t423 = -pkin(7) * qJDD(3) + (t469 - t506) * qJD(3);
t422 = -t331 * t329 + t430;
t421 = t401 * t437 - t486 * t403 + t447;
t420 = -qJDD(3) * t384 + (t387 * t404 - t468) * qJD(3);
t419 = (-g(1) * (-t386 - t482) + g(2) * t416) * t396;
t418 = t285 * t329 - t430 + t513;
t415 = cos(qJ(1));
t412 = sin(qJ(1));
t385 = -t510 * pkin(3) - pkin(4);
t376 = pkin(3) * t409 + qJ(5);
t361 = -t386 - t508;
t349 = qJ(5) * t496;
t347 = qJ(5) * t497;
t344 = t392 + t391;
t304 = t308 - t508;
t298 = pkin(3) * t492 + t303;
t273 = t277 + t392;
t1 = [t429 + (t420 * t410 + (-t426 - t520) * t413 + t485) * MDP(12) + (t420 * t413 + (t426 - t380) * t410 + t466) * MDP(13) + (t279 * t361 + t327 * t344 + t427 + t516) * MDP(19) + (t273 * t327 + t279 * t304 + t427 + t518) * MDP(21) + (t263 * t307 + t284 * t267 + t264 * t304 + t285 * t273 - t265 * t438 + t280 * t268 - g(1) * (-pkin(1) * t412 - t495) - g(2) * (pkin(1) * t415 + t465) + t419) * MDP(24) + (-t278 * t361 + t329 * t344 - t433 + t515) * MDP(20) + (-t267 * t327 + t268 * t329 + t278 * t438 - t279 * t307 + t439) * MDP(22) + (-t273 * t329 + t278 * t304 + t433 + t517) * MDP(23) + qJDD(1) * MDP(1) + (g(1) * t415 + g(2) * t412) * MDP(3) + (((-qJDD(1) - t402) * t411 + (-qJD(1) - t404) * t476) * pkin(1) + t481) * MDP(6) + ((t402 * t414 - t460) * pkin(1) + t440) * MDP(5) + (g(1) * t412 - g(2) * t415) * MDP(2); t429 + (t423 * t410 + (-t428 - t520) * t413 + t485) * MDP(12) + (t423 * t413 + (t428 - t380) * t410 + t466) * MDP(13) + (-t279 * t386 + t435 * t327 + t421 + t516) * MDP(19) + (t279 * t308 - t446 * t327 + t421 + t518) * MDP(21) + (g(1) * t495 - g(2) * t465 + t263 * t318 + t264 * t308 - t265 * t437 + t486 * t280 + t487 * t284 - t446 * t285 + t419) * MDP(24) + (t278 * t386 + t435 * t329 - t512 + t515) * MDP(20) + (t278 * t437 - t279 * t318 - t487 * t327 + t486 * t329 + t439) * MDP(22) + (t278 * t308 + t446 * t329 + t512 + t517) * MDP(23) + ((-t471 + (-qJD(2) + t404) * t478) * pkin(1) + t481) * MDP(6) + (t440 + t453) * MDP(5); MDP(9) * t494 + MDP(10) * t493 + qJDD(3) * MDP(11) + (-g(3) * t413 + t431 * t410) * MDP(12) + (g(3) * t410 + t431 * t413) * MDP(13) + (t296 * t403 + (-t327 * t492 + t510 * t401 - t403 * t473) * pkin(3) + t422) * MDP(19) + (t297 * t403 + (-t329 * t492 - t401 * t409 - t403 * t457) * pkin(3) + t424) * MDP(20) + (-t298 * t327 - t385 * t401 - t445 * t403 - t418) * MDP(21) + (-t278 * t385 - t279 * t376 + (t284 + t445) * t329 + (t280 - t483) * t327) * MDP(22) + (t298 * t329 + t376 * t401 + t483 * t403 + t425 + t514) * MDP(23) + (t263 * t376 + t265 * t385 - t285 * t298 - g(1) * (t444 * t398 + t349) - g(2) * (t444 * t396 + t347) - g(3) * (t482 + t505) + t483 * t284 + t445 * t280) * MDP(24) + t432 + (-MDP(7) * t488 + t479 * MDP(8)) * t404 ^ 2; (t422 + t501) * MDP(19) + (t424 + t502) * MDP(20) + (-t303 * t327 + t394 - t418 + t501) * MDP(21) + (pkin(4) * t278 - qJ(5) * t279 + (t284 - t292) * t329 + (t280 - t472) * t327) * MDP(22) + (t303 * t329 + 0.2e1 * t389 + 0.2e1 * t393 + t425 - t502) * MDP(23) + (t263 * qJ(5) - t265 * pkin(4) - t285 * t303 - t280 * t292 - g(1) * (-pkin(4) * t498 + t349) - g(2) * (-pkin(4) * t499 + t347) - g(3) * t482 + t472 * t284) * MDP(24) + t432; (-t401 + t500) * MDP(21) + t276 * MDP(22) + (-t403 ^ 2 - t511) * MDP(23) + (-t284 * t403 + t418) * MDP(24);];
tau = t1;
