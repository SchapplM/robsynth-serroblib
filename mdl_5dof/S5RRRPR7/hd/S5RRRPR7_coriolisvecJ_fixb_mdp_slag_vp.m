% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:25
% EndTime: 2019-12-31 21:17:34
% DurationCPUTime: 4.36s
% Computational Cost: add. (3812->356), mult. (9507->500), div. (0->0), fcn. (6836->8), ass. (0->166)
t413 = sin(qJ(3));
t416 = cos(qJ(2));
t485 = cos(qJ(3));
t450 = qJD(1) * t485;
t414 = sin(qJ(2));
t462 = qJD(1) * t414;
t494 = -t413 * t462 + t416 * t450;
t363 = qJD(5) - t494;
t471 = t413 * t416;
t370 = -qJD(1) * t471 - t414 * t450;
t407 = qJD(2) + qJD(3);
t410 = sin(pkin(9));
t411 = cos(pkin(9));
t353 = t370 * t410 + t411 * t407;
t415 = cos(qJ(5));
t493 = t353 * t415;
t412 = sin(qJ(5));
t379 = t410 * t415 + t411 * t412;
t492 = t363 * t379;
t469 = t415 * t411;
t472 = t410 * t412;
t378 = -t469 + t472;
t465 = t363 * t378;
t351 = t370 * t411 - t407 * t410;
t491 = t351 * t415 - t353 * t412;
t457 = qJD(1) * qJD(2);
t490 = -0.2e1 * t457;
t489 = (t414 ^ 2 - t416 ^ 2) * MDP(5);
t336 = t494 * t407;
t381 = t485 * t414 + t471;
t348 = t407 * t381;
t337 = t348 * qJD(1);
t448 = t414 * t457;
t281 = pkin(2) * t448 + pkin(3) * t337 - qJ(4) * t336 + qJD(4) * t370;
t486 = -pkin(7) - pkin(6);
t391 = t486 * t414;
t383 = qJD(1) * t391;
t483 = qJD(2) * pkin(2);
t375 = t383 + t483;
t453 = qJD(2) * t486;
t442 = qJD(1) * t453;
t376 = t414 * t442;
t377 = t416 * t442;
t392 = t486 * t416;
t385 = qJD(1) * t392;
t449 = qJD(3) * t485;
t461 = qJD(3) * t413;
t421 = t375 * t449 + t485 * t376 + t413 * t377 + t385 * t461;
t297 = t407 * qJD(4) + t421;
t261 = t411 * t281 - t297 * t410;
t262 = t410 * t281 + t411 * t297;
t438 = -t261 * t410 + t262 * t411;
t338 = -pkin(3) * t370 - qJ(4) * t494;
t327 = pkin(2) * t462 + t338;
t371 = t413 * t385;
t345 = t485 * t383 + t371;
t295 = t411 * t327 - t345 * t410;
t395 = pkin(2) * t449 + qJD(4);
t445 = t395 * t410 + t295;
t296 = t410 * t327 + t411 * t345;
t444 = t395 * t411 - t296;
t372 = t485 * t385;
t344 = t413 * t383 - t372;
t439 = pkin(2) * t461 - t344;
t488 = t485 * t391 + t413 * t392;
t487 = qJD(1) * t381;
t484 = t411 * pkin(4);
t406 = t411 * pkin(8);
t481 = t336 * t410;
t480 = t336 * t411;
t432 = -t413 * t414 + t485 * t416;
t347 = t407 * t432;
t479 = t347 * t410;
t478 = t494 * t410;
t477 = t494 * t411;
t476 = t381 * t410;
t475 = t381 * t411;
t417 = qJD(2) ^ 2;
t470 = t414 * t417;
t468 = t416 * t417;
t418 = qJD(1) ^ 2;
t467 = t416 * t418;
t455 = t414 * t483;
t293 = pkin(3) * t348 - qJ(4) * t347 - qJD(4) * t381 + t455;
t384 = t414 * t453;
t386 = t416 * t453;
t309 = qJD(3) * t488 + t485 * t384 + t413 * t386;
t265 = t410 * t293 + t411 * t309;
t405 = -pkin(2) * t416 - pkin(1);
t390 = t405 * qJD(1);
t322 = -pkin(3) * t494 + qJ(4) * t370 + t390;
t342 = t413 * t375 - t372;
t328 = qJ(4) * t407 + t342;
t287 = t410 * t322 + t411 * t328;
t458 = qJD(5) * t415;
t464 = t336 * t469 + t353 * t458;
t341 = t485 * t375 + t371;
t299 = t410 * t338 + t411 * t341;
t340 = -pkin(3) * t432 - qJ(4) * t381 + t405;
t355 = t413 * t391 - t485 * t392;
t304 = t410 * t340 + t411 * t355;
t271 = pkin(8) * t353 + t287;
t460 = qJD(5) * t271;
t459 = qJD(5) * t381;
t456 = pkin(8) * t478;
t447 = pkin(1) * t490;
t301 = t375 * t461 + t413 * t376 - t485 * t377 - t385 * t449;
t446 = -t287 * t370 + t301 * t410;
t264 = t411 * t293 - t309 * t410;
t286 = t411 * t322 - t328 * t410;
t298 = t411 * t338 - t341 * t410;
t303 = t411 * t340 - t355 * t410;
t256 = -pkin(8) * t481 + t262;
t268 = -pkin(4) * t494 + pkin(8) * t351 + t286;
t443 = -qJD(5) * t268 - t256;
t404 = -t485 * pkin(2) - pkin(3);
t441 = -t370 * pkin(4) - pkin(8) * t477;
t357 = pkin(4) * t478;
t440 = -t357 + t439;
t254 = t268 * t415 - t271 * t412;
t255 = t268 * t412 + t271 * t415;
t283 = -pkin(4) * t432 - pkin(8) * t475 + t303;
t292 = -pkin(8) * t476 + t304;
t437 = t283 * t415 - t292 * t412;
t436 = t283 * t412 + t292 * t415;
t435 = t286 * t370 - t301 * t411;
t434 = qJD(5) * t351 - t481;
t433 = t286 * t477 + t287 * t478 + t438;
t401 = pkin(2) * t413 + qJ(4);
t374 = t401 * t411 + t406;
t431 = qJD(5) * t374 + t441 + t445;
t373 = (-pkin(8) - t401) * t410;
t430 = -qJD(5) * t373 - t444 - t456;
t389 = qJ(4) * t411 + t406;
t429 = qJD(4) * t410 + qJD(5) * t389 + t298 + t441;
t388 = (-pkin(8) - qJ(4)) * t410;
t428 = -qJD(4) * t411 - qJD(5) * t388 + t299 - t456;
t427 = t370 * t390 - t301;
t276 = pkin(4) * t481 + t301;
t326 = -t407 * pkin(3) + qJD(4) - t341;
t302 = -pkin(4) * t353 + t326;
t426 = t254 * t370 + t276 * t378 + t492 * t302;
t425 = -t255 * t370 + t276 * t379 - t465 * t302;
t424 = t301 * t381 + t326 * t347 - t336 * t488;
t423 = -pkin(3) * t336 - qJ(4) * t337 - (-qJD(4) + t326) * t494;
t422 = t336 * t404 - t337 * t401 - (t326 - t395) * t494;
t266 = t434 * t412 + t464;
t267 = -qJD(5) * t491 + t379 * t336;
t306 = -t351 * t412 - t493;
t420 = (-t266 * t378 - t267 * t379 + t465 * t306 + t491 * t492) * MDP(23) + (t266 * t379 + t465 * t491) * MDP(22) + (t337 * t379 - t465 * t363 - t370 * t491) * MDP(24) + (-t306 * t370 - t337 * t378 - t363 * t492) * MDP(25) + t336 * MDP(13) + (t370 ^ 2 - t494 ^ 2) * MDP(12) + (MDP(11) * t494 + t363 * MDP(26)) * t370 + (-t494 * MDP(13) + (-t370 - t487) * MDP(14)) * t407;
t419 = -t390 * t494 - t421;
t310 = t355 * qJD(3) + t413 * t384 - t485 * t386;
t402 = -pkin(3) - t484;
t387 = t404 - t484;
t332 = t378 * t381;
t331 = t379 * t381;
t324 = pkin(4) * t476 - t488;
t313 = t357 + t342;
t285 = pkin(4) * t479 + t310;
t274 = t379 * t347 + t458 * t475 - t459 * t472;
t273 = -t378 * t347 - t379 * t459;
t263 = -pkin(8) * t479 + t265;
t259 = pkin(4) * t348 - t347 * t406 + t264;
t253 = pkin(4) * t337 - pkin(8) * t480 + t261;
t252 = t415 * t253;
t1 = [0.2e1 * t416 * MDP(4) * t448 + t489 * t490 + MDP(6) * t468 - MDP(7) * t470 + (-pkin(6) * t468 + t414 * t447) * MDP(9) + (pkin(6) * t470 + t416 * t447) * MDP(10) + (t336 * t381 - t347 * t370) * MDP(11) + (t336 * t432 - t337 * t381 + t347 * t494 + t348 * t370) * MDP(12) + (t337 * t405 + t348 * t390 + (-qJD(1) * t432 - t494) * t455) * MDP(16) + (t336 * t405 + t347 * t390 + (-t370 + t487) * t455) * MDP(17) + (-t261 * t432 - t264 * t494 + t286 * t348 + t303 * t337 - t310 * t353 + t424 * t410) * MDP(18) + (t262 * t432 + t265 * t494 - t287 * t348 - t304 * t337 - t310 * t351 + t424 * t411) * MDP(19) + (t264 * t351 + t265 * t353 + (-t261 * t381 - t286 * t347 - t303 * t336) * t411 + (-t262 * t381 - t287 * t347 - t304 * t336) * t410) * MDP(20) + (t261 * t303 + t262 * t304 + t264 * t286 + t265 * t287 - t301 * t488 + t310 * t326) * MDP(21) + (-t266 * t332 - t273 * t491) * MDP(22) + (-t266 * t331 + t267 * t332 - t273 * t306 + t274 * t491) * MDP(23) + (-t266 * t432 + t273 * t363 - t332 * t337 - t348 * t491) * MDP(24) + (t267 * t432 - t274 * t363 - t306 * t348 - t331 * t337) * MDP(25) + (-t337 * t432 + t348 * t363) * MDP(26) + ((t259 * t415 - t263 * t412) * t363 + t437 * t337 - (-t256 * t412 + t252) * t432 + t254 * t348 + t285 * t306 + t324 * t267 + t276 * t331 + t302 * t274 + (t255 * t432 - t436 * t363) * qJD(5)) * MDP(27) + (-(t259 * t412 + t263 * t415) * t363 - t436 * t337 + (t253 * t412 + t256 * t415) * t432 - t255 * t348 - t285 * t491 + t324 * t266 - t276 * t332 + t302 * t273 + (t254 * t432 - t437 * t363) * qJD(5)) * MDP(28) + (t347 * MDP(13) - t348 * MDP(14) - t310 * MDP(16) - t309 * MDP(17)) * t407; t420 + (t344 * t407 + (-t407 * t461 + t462 * t494) * pkin(2) + t427) * MDP(16) - t414 * MDP(4) * t467 + ((t373 * t415 - t374 * t412) * t337 + t387 * t267 + (t430 * t412 - t431 * t415) * t363 + t440 * t306 + t426) * MDP(27) + (-(t373 * t412 + t374 * t415) * t337 + t387 * t266 + (t431 * t412 + t430 * t415) * t363 - t440 * t491 + t425) * MDP(28) + (t345 * t407 + (t370 * t462 - t407 * t449) * pkin(2) + t419) * MDP(17) + (t295 * t494 - t353 * t439 + t422 * t410 + t435) * MDP(18) + (-t351 * t445 + t353 * t444 + t433) * MDP(20) + (-t296 * t494 - t351 * t439 + t422 * t411 + t446) * MDP(19) + (-t445 * t286 + t444 * t287 + t301 * t404 + t439 * t326 + t438 * t401) * MDP(21) + t418 * t489 + (t418 * t414 * MDP(9) + MDP(10) * t467) * pkin(1); t420 + (-pkin(3) * t301 - t286 * t298 - t287 * t299 - t326 * t342 + (-t286 * t410 + t287 * t411) * qJD(4) + t438 * qJ(4)) * MDP(21) + (-(t388 * t412 + t389 * t415) * t337 + t402 * t266 + t313 * t491 + (t429 * t412 + t428 * t415) * t363 + t425) * MDP(28) + (t341 * t407 + t419) * MDP(17) + (t342 * t407 + t427) * MDP(16) + (-t298 * t351 - t299 * t353 + (-t351 * t410 + t353 * t411) * qJD(4) + t433) * MDP(20) + ((t388 * t415 - t389 * t412) * t337 + t402 * t267 - t313 * t306 + (t428 * t412 - t429 * t415) * t363 + t426) * MDP(27) + (t298 * t494 + t342 * t353 + t423 * t410 + t435) * MDP(18) + (-t299 * t494 + t342 * t351 + t423 * t411 + t446) * MDP(19); (t351 * t494 + t481) * MDP(18) + (-t353 * t494 + t480) * MDP(19) + (-t351 ^ 2 - t353 ^ 2) * MDP(20) + (-t286 * t351 - t287 * t353 + t301) * MDP(21) + (-t491 * t363 + t267) * MDP(27) + (t363 * t493 + (t351 * t363 + t434) * t412 + t464) * MDP(28); -t306 ^ 2 * MDP(23) + (t306 * t363 + t464) * MDP(24) + t337 * MDP(26) + (t255 * t363 + t252) * MDP(27) + (t254 * t363 + t302 * t306) * MDP(28) - (MDP(22) * t306 - MDP(23) * t491 + t363 * MDP(25) - t302 * MDP(27)) * t491 + (t434 * MDP(25) - MDP(27) * t460 + t443 * MDP(28)) * t415 + (t434 * MDP(24) + (-qJD(5) * t353 - t480) * MDP(25) + t443 * MDP(27) + (-t253 + t460) * MDP(28)) * t412;];
tauc = t1;
