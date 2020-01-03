% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:34
% EndTime: 2019-12-31 22:12:47
% DurationCPUTime: 6.88s
% Computational Cost: add. (4330->442), mult. (11561->613), div. (0->0), fcn. (8665->8), ass. (0->178)
t410 = sin(pkin(5));
t417 = cos(qJ(2));
t476 = qJD(1) * t417;
t402 = t410 * t476;
t434 = t402 - qJD(3);
t414 = sin(qJ(2));
t411 = cos(pkin(5));
t477 = qJD(1) * t411;
t465 = pkin(1) * t477;
t400 = t414 * t465;
t374 = pkin(7) * t402 + t400;
t413 = sin(qJ(3));
t416 = cos(qJ(3));
t520 = t374 + t434 * (pkin(3) * t413 - pkin(9) * t416);
t448 = qJD(2) + t477;
t478 = qJD(1) * t410;
t462 = t414 * t478;
t421 = -t413 * t462 + t416 * t448;
t467 = qJD(1) * qJD(2);
t456 = t410 * t467;
t440 = t417 * t456;
t431 = t416 * t440;
t419 = t421 * qJD(3) + t431;
t519 = -qJD(4) * t434 + t419;
t509 = -qJ(5) - pkin(9);
t518 = qJ(5) * t421 + qJD(4) * t509;
t407 = t410 ^ 2;
t517 = -0.2e1 * t407 * t467;
t516 = MDP(5) * (t414 ^ 2 - t417 ^ 2);
t497 = t410 * t414;
t403 = pkin(7) * t497;
t511 = pkin(1) * t417;
t368 = t403 + (-pkin(2) - t511) * t411;
t379 = -t411 * t416 + t413 * t497;
t380 = t411 * t413 + t416 * t497;
t310 = pkin(3) * t379 - pkin(9) * t380 + t368;
t496 = t410 * t417;
t466 = pkin(7) * t496;
t512 = pkin(1) * t414;
t369 = t466 + (pkin(8) + t512) * t411;
t370 = (-pkin(2) * t417 - pkin(8) * t414 - pkin(1)) * t410;
t486 = t416 * t369 + t413 * t370;
t312 = -pkin(9) * t496 + t486;
t412 = sin(qJ(4));
t415 = cos(qJ(4));
t488 = t412 * t310 + t415 * t312;
t474 = qJD(3) * t413;
t464 = pkin(8) * t474;
t515 = t412 * t464 - t415 * t520;
t371 = -pkin(7) * t462 + t417 * t465;
t428 = (pkin(2) * t414 - pkin(8) * t417) * t410;
t372 = qJD(1) * t428;
t485 = t416 * t371 + t413 * t372;
t314 = pkin(9) * t462 + t485;
t395 = -pkin(3) * t416 - pkin(9) * t413 - pkin(2);
t470 = qJD(4) * t415;
t514 = t415 * t314 - t395 * t470 + t412 * t520;
t358 = t448 * t413 + t416 * t462;
t325 = t415 * t358 - t412 * t434;
t513 = t325 ^ 2;
t418 = qJD(1) ^ 2;
t510 = pkin(4) * t412;
t508 = pkin(8) * qJD(3);
t506 = qJ(5) * t413;
t441 = t414 * t456;
t341 = pkin(8) * t448 + t374;
t351 = qJD(1) * t370;
t373 = qJD(2) * t428;
t365 = qJD(1) * t373;
t375 = (t411 * t511 - t403) * qJD(2);
t366 = qJD(1) * t375;
t472 = qJD(3) * t416;
t445 = -t341 * t472 - t351 * t474 + t416 * t365 - t413 * t366;
t280 = -pkin(3) * t441 - t445;
t505 = t280 * t412;
t504 = t280 * t415;
t471 = qJD(4) * t412;
t288 = t358 * t471 - t412 * t441 - t519 * t415;
t503 = t288 * t412;
t323 = t358 * t412 + t415 * t434;
t352 = qJD(4) - t421;
t502 = t323 * t352;
t501 = t325 * t352;
t432 = t413 * t440;
t327 = t358 * qJD(3) + t432;
t500 = t327 * t412;
t499 = t327 * t415;
t498 = t407 * t418;
t495 = t413 * t415;
t494 = t415 * t416;
t493 = t415 * t417;
t340 = -pkin(2) * t448 - t371;
t298 = -pkin(3) * t421 - t358 * pkin(9) + t340;
t307 = t416 * t341 + t413 * t351;
t301 = -pkin(9) * t434 + t307;
t272 = t415 * t298 - t301 * t412;
t268 = -qJ(5) * t325 + t272;
t266 = pkin(4) * t352 + t268;
t492 = t266 - t268;
t346 = (t412 * t414 + t416 * t493) * t478;
t405 = pkin(8) * t494;
t444 = t413 * t402;
t469 = qJD(5) * t415;
t491 = -pkin(4) * t444 + qJ(5) * t346 + t314 * t412 - t413 * t469 + (pkin(4) * t413 - qJ(5) * t494) * qJD(3) + (-t405 + (-t395 + t506) * t412) * qJD(4) + t515;
t443 = t416 * t402;
t345 = t412 * t443 - t415 * t462;
t490 = qJ(5) * t345 + (-qJ(5) * qJD(4) - t508) * t495 + (-qJD(5) * t413 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t416) * t412 - t514;
t306 = -t413 * t341 + t416 * t351;
t317 = pkin(3) * t358 - pkin(9) * t421;
t489 = t415 * t306 + t412 * t317;
t484 = t518 * t412 + t469 - t489;
t316 = t415 * t317;
t483 = -pkin(4) * t358 - t316 + t518 * t415 + (-qJD(5) + t306) * t412;
t480 = t412 * t395 + t405;
t475 = qJD(2) * t414;
t473 = qJD(3) * t415;
t468 = t340 * qJD(3);
t463 = t412 * t496;
t461 = t410 * t475;
t460 = qJD(2) * t496;
t459 = t352 * t471;
t454 = t415 * t310 - t312 * t412;
t453 = -t413 * t369 + t370 * t416;
t452 = t417 * t434;
t451 = t352 * t415;
t450 = qJD(3) * t434;
t446 = MDP(4) * t407 * t414 * t417;
t439 = MDP(15) * t462;
t438 = pkin(1) * t517;
t311 = pkin(3) * t496 - t453;
t436 = -t412 * t472 + t345;
t435 = t415 * t472 - t346;
t433 = pkin(7) * t440;
t273 = t298 * t412 + t301 * t415;
t359 = t413 * t371;
t430 = pkin(3) * t462 - t359;
t429 = -t369 * t472 - t370 * t474 + t373 * t416 - t413 * t375;
t334 = t380 * t412 + t410 * t493;
t427 = -t352 * t470 - t500;
t300 = pkin(3) * t434 - t306;
t426 = -pkin(9) * t327 + t300 * t352;
t425 = t341 * t474 - t351 * t472 - t413 * t365 - t416 * t366;
t424 = -t369 * t474 + t370 * t472 + t413 * t373 + t416 * t375;
t279 = pkin(9) * t441 - t425;
t287 = t327 * pkin(3) - pkin(9) * t419 + qJD(2) * t400 + t433;
t264 = t415 * t279 + t412 * t287 + t298 * t470 - t301 * t471;
t283 = pkin(9) * t461 + t424;
t332 = qJD(3) * t380 + t413 * t460;
t333 = -qJD(3) * t379 + t416 * t460;
t376 = (t411 * t512 + t466) * qJD(2);
t292 = pkin(3) * t332 - pkin(9) * t333 + t376;
t423 = t415 * t283 + t412 * t292 + t310 * t470 - t312 * t471;
t422 = pkin(1) * (-t411 * t467 + t498);
t284 = -pkin(3) * t461 - t429;
t265 = -qJD(4) * t273 - t279 * t412 + t415 * t287;
t420 = -t488 * qJD(4) - t283 * t412 + t415 * t292;
t289 = t358 * t470 + t519 * t412 - t415 * t441;
t267 = pkin(4) * t289 + t280;
t397 = t509 * t415;
t396 = t509 * t412;
t386 = t415 * t395;
t367 = qJD(1) * t376;
t338 = -t412 * t506 + t480;
t335 = t380 * t415 - t463;
t328 = -qJ(5) * t495 + t386 + (-pkin(8) * t412 - pkin(4)) * t416;
t322 = t323 ^ 2;
t313 = -t372 * t416 - t430;
t297 = -qJD(4) * t334 + t333 * t415 + t412 * t461;
t296 = -qJD(4) * t463 + t333 * t412 + t380 * t470 - t415 * t461;
t281 = t323 * pkin(4) + qJD(5) + t300;
t274 = -qJ(5) * t334 + t488;
t270 = pkin(4) * t379 - qJ(5) * t335 + t454;
t269 = -qJ(5) * t323 + t273;
t263 = -qJ(5) * t296 - qJD(5) * t334 + t423;
t262 = pkin(4) * t332 - qJ(5) * t297 - qJD(5) * t335 + t420;
t261 = -qJ(5) * t289 - qJD(5) * t323 + t264;
t260 = pkin(4) * t327 + qJ(5) * t288 - qJD(5) * t325 + t265;
t1 = [0.2e1 * t446 * t467 + t516 * t517 + (-t367 * t411 - t376 * t448 + t414 * t438) * MDP(9) + (-t366 * t411 - t375 * t448 + t417 * t438) * MDP(10) + (t358 * t333 + t380 * t419) * MDP(11) + (-t380 * t327 - t358 * t332 + t333 * t421 - t379 * t419) * MDP(12) + (-t333 * t434 + t358 * t461 + t380 * t441 - t419 * t496) * MDP(13) + (t332 * t434 + (t327 * t417 + (-qJD(1) * t379 + t421) * t475) * t410) * MDP(14) + (-t407 * t476 - t410 * t434) * MDP(15) * t475 + (-t429 * t434 - t376 * t421 + t368 * t327 + t367 * t379 + t340 * t332 + (-t445 * t417 + (qJD(1) * t453 + t306) * t475) * t410) * MDP(16) + (-t307 * t461 + t340 * t333 + t376 * t358 + t367 * t380 + t368 * t419 + t424 * t434 - t425 * t496 - t441 * t486) * MDP(17) + (-t288 * t335 + t297 * t325) * MDP(18) + (t288 * t334 - t289 * t335 - t296 * t325 - t297 * t323) * MDP(19) + (-t288 * t379 + t297 * t352 + t325 * t332 + t327 * t335) * MDP(20) + (-t289 * t379 - t296 * t352 - t323 * t332 - t327 * t334) * MDP(21) + (t327 * t379 + t332 * t352) * MDP(22) + (t265 * t379 + t272 * t332 + t280 * t334 + t284 * t323 + t311 * t289 + t300 * t296 + t327 * t454 + t352 * t420) * MDP(23) + (-t264 * t379 - t273 * t332 + t280 * t335 + t284 * t325 - t311 * t288 + t300 * t297 - t327 * t488 - t352 * t423) * MDP(24) + (-t260 * t335 - t261 * t334 - t262 * t325 - t263 * t323 - t266 * t297 - t269 * t296 + t270 * t288 - t274 * t289) * MDP(25) + (t261 * t274 + t269 * t263 + t260 * t270 + t266 * t262 + t267 * (pkin(4) * t334 + t311) + t281 * (pkin(4) * t296 + t284)) * MDP(26) + (MDP(6) * t460 - MDP(7) * t461) * (qJD(2) + 0.2e1 * t477); t498 * t516 + (t374 * t448 + t414 * t422 - t433) * MDP(9) + (pkin(7) * t441 + t371 * t448 + t417 * t422) * MDP(10) + (-qJD(3) * t413 ^ 2 * t462 + ((qJD(3) * t448 + t440) * t413 - t434 * t358) * t416) * MDP(11) + (-t413 * t327 + t419 * t416 + (t444 - t474) * t358 - (t443 - t472) * t421) * MDP(12) + (-t416 * t450 + (t416 * t452 + (qJD(2) * t413 - t358) * t414) * t478) * MDP(13) + (t413 * t450 + (-t413 * t452 + (qJD(2) * t416 - t421) * t414) * t478) * MDP(14) + t434 * t439 + (-pkin(2) * t327 + t413 * t468 - t359 * t434 + t374 * t421 + (pkin(8) * t450 + t372 * t434 - t367) * t416 + (-t306 * t414 + (-pkin(8) * t475 - t340 * t417) * t413) * t478) * MDP(16) + (-pkin(2) * t419 + t307 * t462 - t340 * t443 - t374 * t358 + t367 * t413 + (-t464 - t485) * t434 + (-pkin(8) * t441 + t468) * t416) * MDP(17) + (-t288 * t495 + (-t413 * t471 + t435) * t325) * MDP(18) + (t323 * t346 + t325 * t345 + (-t323 * t415 - t325 * t412) * t472 + (t503 - t289 * t415 + (t323 * t412 - t325 * t415) * qJD(4)) * t413) * MDP(19) + (t288 * t416 + t435 * t352 + (-t325 * t434 - t459 + t499) * t413) * MDP(20) + (t289 * t416 + t436 * t352 + (t323 * t434 + t427) * t413) * MDP(21) + (-t352 * t413 * t434 - t327 * t416) * MDP(22) + (-t300 * t345 - t313 * t323 + t386 * t327 + ((-qJD(4) * t395 + t314) * t412 + t515) * t352 + (t300 * t412 * qJD(3) - t265 + (qJD(3) * t323 + t427) * pkin(8)) * t416 + (pkin(8) * t289 - t272 * t434 + t300 * t470 + t505) * t413) * MDP(23) + (-t480 * t327 - t313 * t325 - t300 * t346 + t514 * t352 + (t300 * t473 + t264 + (qJD(3) * t325 + t459) * pkin(8)) * t416 + (-t300 * t471 + t504 + t434 * t273 + (t352 * t473 - t288) * pkin(8)) * t413) * MDP(24) + (t266 * t346 + t269 * t345 + t288 * t328 - t289 * t338 - t491 * t325 - t490 * t323 + (-t266 * t415 - t269 * t412) * t472 + (-t260 * t415 - t261 * t412 + (t266 * t412 - t269 * t415) * qJD(4)) * t413) * MDP(25) + (t260 * t328 + t261 * t338 + t267 * (pkin(8) + t510) * t413 + ((t372 + t508) * t416 + (t413 * t470 - t436) * pkin(4) + t430) * t281 + t490 * t269 + t491 * t266) * MDP(26) + (-t446 + (-t417 * MDP(6) + t414 * MDP(7)) * t410 * t411) * t418; -t421 ^ 2 * MDP(12) + (t421 * t402 + t431) * MDP(13) - t432 * MDP(14) + qJD(2) * t439 + (-t307 * t434 + t445) * MDP(16) + (-t306 * t434 - t340 * t421 + t425) * MDP(17) + (t325 * t451 - t503) * MDP(18) + ((-t288 - t502) * t415 + (-t289 - t501) * t412) * MDP(19) + (t352 * t451 + t500) * MDP(20) + (-t352 ^ 2 * t412 + t499) * MDP(21) + (-pkin(3) * t289 - t504 - t307 * t323 + (-pkin(9) * t470 - t316) * t352 + (t306 * t352 + t426) * t412) * MDP(23) + (pkin(3) * t288 + t505 - t307 * t325 + (pkin(9) * t471 + t489) * t352 + t426 * t415) * MDP(24) + (t288 * t396 + t289 * t397 - t483 * t325 - t484 * t323 + (-t266 * t352 + t261) * t415 + (-t269 * t352 - t260) * t412) * MDP(25) + (-t261 * t397 + t260 * t396 + t267 * (-pkin(4) * t415 - pkin(3)) + (t352 * t510 - t307) * t281 + t484 * t269 + t483 * t266) * MDP(26) + (-MDP(11) * t421 + MDP(12) * t358 - t402 * MDP(14) - t340 * MDP(16) - t325 * MDP(20) + t323 * MDP(21) - t352 * MDP(22) - t272 * MDP(23) + t273 * MDP(24)) * t358; t325 * t323 * MDP(18) + (-t322 + t513) * MDP(19) + (-t288 + t502) * MDP(20) + (-t289 + t501) * MDP(21) + t327 * MDP(22) + (t273 * t352 - t300 * t325 + t265) * MDP(23) + (t272 * t352 + t300 * t323 - t264) * MDP(24) + (pkin(4) * t288 - t323 * t492) * MDP(25) + (t492 * t269 + (-t281 * t325 + t260) * pkin(4)) * MDP(26); (-t322 - t513) * MDP(25) + (t266 * t325 + t269 * t323 + t267) * MDP(26);];
tauc = t1;
