% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRP11
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:47
% EndTime: 2019-12-31 22:19:02
% DurationCPUTime: 8.38s
% Computational Cost: add. (5815->505), mult. (15401->680), div. (0->0), fcn. (11481->8), ass. (0->187)
t420 = sin(pkin(5));
t427 = cos(qJ(2));
t492 = qJD(1) * t427;
t412 = t420 * t492;
t448 = t412 - qJD(3);
t424 = sin(qJ(2));
t421 = cos(pkin(5));
t493 = qJD(1) * t421;
t479 = pkin(1) * t493;
t411 = t424 * t479;
t389 = pkin(7) * t412 + t411;
t423 = sin(qJ(3));
t426 = cos(qJ(3));
t541 = t389 + t448 * (pkin(3) * t423 - pkin(9) * t426);
t489 = qJD(3) * t423;
t540 = pkin(8) * t489;
t464 = qJD(2) + t493;
t494 = qJD(1) * t420;
t476 = t424 * t494;
t430 = -t423 * t476 + t426 * t464;
t481 = qJD(1) * qJD(2);
t470 = t420 * t481;
t456 = t427 * t470;
t441 = t426 * t456;
t429 = qJD(3) * t430 + t441;
t538 = -qJD(4) * t448 + t429;
t537 = t423 * t412 - t489;
t366 = qJD(4) - t430;
t417 = t420 ^ 2;
t536 = -0.2e1 * t417 * t481;
t535 = MDP(5) * (t424 ^ 2 - t427 ^ 2);
t353 = pkin(8) * t464 + t389;
t384 = (-pkin(2) * t427 - pkin(8) * t424 - pkin(1)) * t420;
t365 = qJD(1) * t384;
t438 = (pkin(2) * t424 - pkin(8) * t427) * t420;
t388 = qJD(2) * t438;
t379 = qJD(1) * t388;
t511 = t420 * t424;
t413 = pkin(7) * t511;
t526 = pkin(1) * t427;
t390 = (t421 * t526 - t413) * qJD(2);
t380 = qJD(1) * t390;
t487 = qJD(3) * t426;
t434 = t353 * t489 - t365 * t487 - t423 * t379 - t426 * t380;
t457 = t424 * t470;
t291 = pkin(9) * t457 - t434;
t372 = t423 * t464 + t426 * t476;
t442 = t423 * t456;
t341 = t372 * qJD(3) + t442;
t447 = pkin(7) * t456;
t303 = t341 * pkin(3) - pkin(9) * t429 + qJD(2) * t411 + t447;
t386 = -pkin(7) * t476 + t427 * t479;
t352 = -pkin(2) * t464 - t386;
t314 = -pkin(3) * t430 - t372 * pkin(9) + t352;
t322 = t426 * t353 + t423 * t365;
t316 = -pkin(9) * t448 + t322;
t422 = sin(qJ(4));
t425 = cos(qJ(4));
t485 = qJD(4) * t425;
t486 = qJD(4) * t422;
t274 = -t422 * t291 + t425 * t303 - t314 * t486 - t316 * t485;
t525 = pkin(4) * t341;
t272 = -t274 - t525;
t284 = t314 * t422 + t316 * t425;
t280 = qJ(5) * t366 + t284;
t534 = -t280 * t366 + t272;
t382 = t413 + (-pkin(2) - t526) * t421;
t393 = -t421 * t426 + t423 * t511;
t394 = t421 * t423 + t426 * t511;
t325 = pkin(3) * t393 - pkin(9) * t394 + t382;
t510 = t420 * t427;
t480 = pkin(7) * t510;
t527 = pkin(1) * t424;
t383 = t480 + (pkin(8) + t527) * t421;
t499 = t426 * t383 + t423 * t384;
t327 = -pkin(9) * t510 + t499;
t533 = t422 * t325 + t425 * t327;
t387 = qJD(1) * t438;
t498 = t426 * t386 + t423 * t387;
t330 = pkin(9) * t476 + t498;
t408 = -pkin(3) * t426 - pkin(9) * t423 - pkin(2);
t532 = t425 * t330 - t408 * t485 + t541 * t422;
t321 = -t423 * t353 + t426 * t365;
t315 = pkin(3) * t448 - t321;
t337 = t372 * t422 + t425 * t448;
t339 = t425 * t372 - t422 * t448;
t282 = t337 * pkin(4) - t339 * qJ(5) + t315;
t523 = pkin(9) * t341;
t531 = t282 * t366 - t523;
t433 = -t383 * t489 + t384 * t487 + t423 * t388 + t426 * t390;
t491 = qJD(2) * t424;
t475 = t420 * t491;
t299 = pkin(9) * t475 + t433;
t474 = qJD(2) * t510;
t345 = qJD(3) * t394 + t423 * t474;
t346 = -qJD(3) * t393 + t426 * t474;
t391 = (t421 * t527 + t480) * qJD(2);
t308 = pkin(3) * t345 - pkin(9) * t346 + t391;
t530 = -qJD(4) * t533 - t299 * t422 + t308 * t425;
t529 = t339 ^ 2;
t528 = t366 ^ 2;
t428 = qJD(1) ^ 2;
t524 = pkin(8) * t426;
t522 = qJ(5) * t341;
t461 = -t353 * t487 - t365 * t489 + t426 * t379 - t423 * t380;
t292 = -pkin(3) * t457 - t461;
t304 = t372 * t486 - t422 * t457 - t425 * t538;
t305 = t372 * t485 + t422 * t538 - t425 * t457;
t276 = pkin(4) * t305 + qJ(5) * t304 - qJD(5) * t339 + t292;
t521 = t276 * t422;
t520 = t276 * t425;
t518 = t282 * t339;
t517 = t292 * t422;
t516 = t304 * t422;
t515 = t337 * t366;
t514 = t339 * t337;
t513 = t339 * t366;
t512 = t417 * t428;
t509 = t422 * t341;
t508 = t425 * t341;
t507 = t425 * t427;
t484 = qJD(4) * t426;
t488 = qJD(3) * t425;
t506 = qJD(5) * t426 - (-t422 * t484 - t423 * t488) * pkin(8) + t532 + t537 * qJ(5);
t505 = -t408 * t486 + (t330 + t540) * t422 + (-pkin(8) * t484 - t541) * t425 - t537 * pkin(4);
t449 = pkin(4) * t422 - qJ(5) * t425;
t504 = qJD(5) * t422 - t366 * t449 + t322;
t373 = t423 * t386;
t329 = -pkin(3) * t476 - t387 * t426 + t373;
t459 = t426 * t412;
t359 = t422 * t459 - t425 * t476;
t360 = (t422 * t424 + t426 * t507) * t494;
t440 = pkin(8) + t449;
t450 = pkin(4) * t425 + qJ(5) * t422;
t503 = pkin(4) * t359 - qJ(5) * t360 + t329 - (qJD(4) * t450 - qJD(5) * t425) * t423 - t440 * t487;
t333 = pkin(3) * t372 - pkin(9) * t430;
t502 = t425 * t321 + t422 * t333;
t496 = t422 * t408 + t425 * t524;
t490 = qJD(3) * t422;
t483 = t352 * qJD(3);
t283 = t314 * t425 - t316 * t422;
t482 = qJD(5) - t283;
t478 = pkin(9) * t486;
t477 = t422 * t510;
t473 = t366 * t486;
t469 = -t423 * t383 + t384 * t426;
t468 = t427 * t448;
t467 = t366 * t425;
t466 = qJD(3) * t448;
t462 = MDP(4) * t417 * t424 * t427;
t455 = MDP(15) * t476;
t454 = pkin(1) * t536;
t326 = pkin(3) * t510 - t469;
t452 = t422 * t487 - t359;
t451 = t425 * t487 - t360;
t279 = -pkin(4) * t366 + t482;
t446 = t279 * t425 - t280 * t422;
t444 = t325 * t425 - t327 * t422;
t439 = -t383 * t487 - t384 * t489 + t388 * t426 - t423 * t390;
t347 = t394 * t422 + t420 * t507;
t437 = -t366 * t485 - t509;
t436 = t315 * t366 - t523;
t435 = t284 * t366 + t274;
t273 = t425 * t291 + t422 * t303 + t314 * t485 - t316 * t486;
t432 = t425 * t299 + t422 * t308 + t325 * t485 - t327 * t486;
t431 = pkin(1) * (-t421 * t481 + t512);
t300 = -pkin(3) * t475 - t439;
t404 = -pkin(3) - t450;
t385 = t440 * t423;
t381 = qJD(1) * t391;
t355 = -t408 * t425 + (pkin(8) * t422 + pkin(4)) * t426;
t354 = -qJ(5) * t426 + t496;
t348 = t394 * t425 - t477;
t312 = -qJD(4) * t347 + t346 * t425 + t422 * t475;
t311 = -qJD(4) * t477 + t346 * t422 + t394 * t485 - t425 * t475;
t307 = pkin(4) * t339 + qJ(5) * t337;
t293 = pkin(4) * t347 - qJ(5) * t348 + t326;
t288 = -pkin(4) * t393 - t444;
t287 = qJ(5) * t393 + t533;
t286 = -pkin(4) * t372 + t321 * t422 - t333 * t425;
t285 = qJ(5) * t372 + t502;
t281 = -t304 + t515;
t278 = pkin(4) * t311 - qJ(5) * t312 - qJD(5) * t348 + t300;
t277 = -pkin(4) * t345 - t530;
t275 = qJ(5) * t345 + qJD(5) * t393 + t432;
t271 = qJD(5) * t366 + t273 + t522;
t1 = [t535 * t536 + (-t381 * t421 - t391 * t464 + t424 * t454) * MDP(9) + (-t380 * t421 - t390 * t464 + t427 * t454) * MDP(10) + (t372 * t346 + t394 * t429) * MDP(11) + (-t394 * t341 - t372 * t345 + t346 * t430 - t393 * t429) * MDP(12) + (-t346 * t448 + t372 * t475 + t394 * t457 - t429 * t510) * MDP(13) + (t345 * t448 + (t341 * t427 + (-qJD(1) * t393 + t430) * t491) * t420) * MDP(14) + (-t417 * t492 - t420 * t448) * MDP(15) * t491 + (-t439 * t448 - t391 * t430 + t382 * t341 + t381 * t393 + t352 * t345 + (-t461 * t427 + (qJD(1) * t469 + t321) * t491) * t420) * MDP(16) + (-t322 * t475 + t352 * t346 + t391 * t372 + t381 * t394 + t382 * t429 + t433 * t448 - t434 * t510 - t457 * t499) * MDP(17) + (-t304 * t348 + t312 * t339) * MDP(18) + (t304 * t347 - t305 * t348 - t311 * t339 - t312 * t337) * MDP(19) + (-t304 * t393 + t312 * t366 + t339 * t345 + t341 * t348) * MDP(20) + (-t305 * t393 - t311 * t366 - t337 * t345 - t341 * t347) * MDP(21) + (t341 * t393 + t345 * t366) * MDP(22) + (t274 * t393 + t283 * t345 + t292 * t347 + t300 * t337 + t326 * t305 + t315 * t311 + t444 * t341 + t366 * t530) * MDP(23) + (-t273 * t393 - t284 * t345 + t292 * t348 + t300 * t339 - t326 * t304 + t315 * t312 - t341 * t533 - t366 * t432) * MDP(24) + (-t272 * t393 + t276 * t347 - t277 * t366 + t278 * t337 - t279 * t345 + t282 * t311 - t288 * t341 + t293 * t305) * MDP(25) + (-t271 * t347 + t272 * t348 - t275 * t337 + t277 * t339 + t279 * t312 - t280 * t311 - t287 * t305 - t288 * t304) * MDP(26) + (t271 * t393 + t275 * t366 - t276 * t348 - t278 * t339 + t280 * t345 - t282 * t312 + t287 * t341 + t293 * t304) * MDP(27) + (t271 * t287 + t272 * t288 + t275 * t280 + t276 * t293 + t277 * t279 + t278 * t282) * MDP(28) + 0.2e1 * t462 * t481 + (MDP(6) * t474 - MDP(7) * t475) * (qJD(2) + 0.2e1 * t493); (t408 * t508 - t315 * t359 - t329 * t337 + (-t541 * t425 + (-qJD(4) * t408 + t330) * t422) * t366 + (t315 * t490 - t274 + (qJD(3) * t337 + t437) * pkin(8)) * t426 + (t315 * t485 + t517 - t448 * t283 + (t366 * t490 + t305) * pkin(8)) * t423) * MDP(23) + (-t496 * t341 - t329 * t339 - t315 * t360 + t532 * t366 + (t315 * t488 + t273 + (qJD(3) * t339 + t473) * pkin(8)) * t426 + (-t315 * t486 + t292 * t425 + t448 * t284 + (t366 * t488 - t304) * pkin(8)) * t423) * MDP(24) + (t272 * t426 + t305 * t385 - t341 * t355 + t505 * t366 - t503 * t337 + t452 * t282 + (t279 * t448 + t282 * t485 + t521) * t423) * MDP(25) + (-t279 * t360 + t280 * t359 - t304 * t355 - t305 * t354 - t505 * t339 + t506 * t337 + t446 * t487 + (-t271 * t422 + t272 * t425 + (-t279 * t422 - t280 * t425) * qJD(4)) * t423) * MDP(26) + (-t271 * t426 + t304 * t385 + t341 * t354 - t506 * t366 + t503 * t339 - t451 * t282 + (-t280 * t448 + t282 * t486 - t520) * t423) * MDP(27) + (t271 * t354 + t272 * t355 + t276 * t385 - t279 * t505 - t280 * t506 - t282 * t503) * MDP(28) + t512 * t535 + (t389 * t464 + t424 * t431 - t447) * MDP(9) + (pkin(7) * t457 + t386 * t464 + t427 * t431) * MDP(10) + (-qJD(3) * t423 ^ 2 * t476 + ((qJD(3) * t464 + t456) * t423 - t448 * t372) * t426) * MDP(11) + (-t423 * t341 + t426 * t429 + t537 * t372 - (t459 - t487) * t430) * MDP(12) + (-t426 * t466 + (t426 * t468 + (qJD(2) * t423 - t372) * t424) * t494) * MDP(13) + (t423 * t466 + (-t423 * t468 + (qJD(2) * t426 - t430) * t424) * t494) * MDP(14) + (-pkin(2) * t341 + t423 * t483 - t373 * t448 + t389 * t430 + (pkin(8) * t466 + t387 * t448 - t381) * t426 + (-t321 * t424 + (-pkin(8) * t491 - t352 * t427) * t423) * t494) * MDP(16) + (-pkin(2) * t429 + t322 * t476 - t352 * t459 - t389 * t372 + t381 * t423 + t426 * t483 - t457 * t524 + (-t498 - t540) * t448) * MDP(17) + (-t304 * t423 * t425 + (-t423 * t486 + t451) * t339) * MDP(18) + (t337 * t360 + t339 * t359 + (-t337 * t425 - t339 * t422) * t487 + (t516 - t305 * t425 + (t337 * t422 - t339 * t425) * qJD(4)) * t423) * MDP(19) + (t304 * t426 + t451 * t366 + (-t339 * t448 - t473 + t508) * t423) * MDP(20) + (t305 * t426 - t452 * t366 + (t337 * t448 + t437) * t423) * MDP(21) + (-t366 * t423 * t448 - t341 * t426) * MDP(22) + t448 * t455 + (-t462 + (-MDP(6) * t427 + MDP(7) * t424) * t420 * t421) * t428; -t430 ^ 2 * MDP(12) + (t430 * t412 + t441) * MDP(13) - t442 * MDP(14) + qJD(2) * t455 + (-t322 * t448 + t461) * MDP(16) + (-t321 * t448 - t352 * t430 + t434) * MDP(17) + (t339 * t467 - t516) * MDP(18) + ((-t304 - t515) * t425 + (-t305 - t513) * t422) * MDP(19) + (t366 * t467 + t509) * MDP(20) + (-t422 * t528 + t508) * MDP(21) + (-pkin(3) * t305 - t322 * t337 + (-t292 + (-pkin(9) * qJD(4) - t333) * t366) * t425 + (t321 * t366 + t436) * t422) * MDP(23) + (pkin(3) * t304 + t517 - t322 * t339 + (t478 + t502) * t366 + t436 * t425) * MDP(24) + (-t520 + t305 * t404 + (-pkin(9) * t485 + t286) * t366 - t504 * t337 + t531 * t422) * MDP(25) + (t285 * t337 - t286 * t339 + (t271 + t366 * t279 + (qJD(4) * t339 - t305) * pkin(9)) * t425 + ((qJD(4) * t337 - t304) * pkin(9) + t534) * t422) * MDP(26) + (-t521 + t304 * t404 + (-t285 - t478) * t366 + t504 * t339 - t531 * t425) * MDP(27) + (t276 * t404 - t279 * t286 - t280 * t285 - t504 * t282 + (qJD(4) * t446 + t271 * t425 + t272 * t422) * pkin(9)) * MDP(28) + (-MDP(11) * t430 + MDP(12) * t372 - t412 * MDP(14) - t352 * MDP(16) - t339 * MDP(20) + t337 * MDP(21) - t366 * MDP(22) - t283 * MDP(23) + t284 * MDP(24) + t279 * MDP(25) - t280 * MDP(27)) * t372; MDP(18) * t514 + (-t337 ^ 2 + t529) * MDP(19) + t281 * MDP(20) + (-t305 + t513) * MDP(21) + t341 * MDP(22) + (-t315 * t339 + t435) * MDP(23) + (t283 * t366 + t315 * t337 - t273) * MDP(24) + (-t307 * t337 + t435 - t518 + 0.2e1 * t525) * MDP(25) + (pkin(4) * t304 - qJ(5) * t305 + (t280 - t284) * t339 + (t279 - t482) * t337) * MDP(26) + (0.2e1 * t522 - t282 * t337 + t307 * t339 + (0.2e1 * qJD(5) - t283) * t366 + t273) * MDP(27) + (-pkin(4) * t272 + qJ(5) * t271 - t279 * t284 + t280 * t482 - t282 * t307) * MDP(28); (-t341 + t514) * MDP(25) + t281 * MDP(26) + (-t528 - t529) * MDP(27) + (t518 + t534) * MDP(28);];
tauc = t1;
