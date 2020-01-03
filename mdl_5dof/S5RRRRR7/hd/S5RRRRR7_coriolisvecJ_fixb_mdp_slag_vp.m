% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:35
% EndTime: 2019-12-31 22:22:43
% DurationCPUTime: 4.19s
% Computational Cost: add. (4609->343), mult. (11765->464), div. (0->0), fcn. (8834->8), ass. (0->171)
t397 = cos(qJ(5));
t452 = qJD(5) * t397;
t395 = sin(qJ(3));
t396 = sin(qJ(2));
t458 = qJD(1) * t396;
t445 = t395 * t458;
t399 = cos(qJ(3));
t400 = cos(qJ(2));
t457 = qJD(1) * t400;
t446 = t399 * t457;
t353 = -t445 + t446;
t354 = -t395 * t457 - t399 * t458;
t394 = sin(qJ(4));
t398 = cos(qJ(4));
t328 = t398 * t353 + t354 * t394;
t505 = t328 * t397;
t508 = t452 - t505;
t366 = t395 * t400 + t396 * t399;
t390 = qJD(2) + qJD(3);
t339 = t390 * t366;
t334 = t339 * qJD(1);
t490 = pkin(6) + pkin(7);
t447 = qJD(2) * t490;
t424 = qJD(1) * t447;
t363 = t400 * t424;
t374 = t490 * t400;
t369 = qJD(1) * t374;
t456 = qJD(3) * t395;
t434 = -t395 * t363 - t369 * t456;
t373 = t490 * t396;
t367 = qJD(1) * t373;
t487 = qJD(2) * pkin(2);
t361 = -t367 + t487;
t362 = t396 * t424;
t497 = (qJD(3) * t361 - t362) * t399;
t286 = -pkin(8) * t334 + t434 + t497;
t350 = t354 * pkin(8);
t355 = t395 * t369;
t436 = t399 * t361 - t355;
t309 = t350 + t436;
t305 = pkin(3) * t390 + t309;
t449 = qJD(1) * qJD(2);
t444 = t400 * t449;
t333 = qJD(3) * t446 - t390 * t445 + t399 * t444;
t359 = t399 * t369;
t415 = -t361 * t395 - t359;
t435 = t395 * t362 - t399 * t363;
t406 = t415 * qJD(3) + t435;
t287 = -pkin(8) * t333 + t406;
t488 = pkin(8) * t353;
t310 = -t415 + t488;
t455 = qJD(4) * t394;
t438 = t287 * t394 - t310 * t455;
t254 = (qJD(4) * t305 + t286) * t398 + t438;
t386 = -pkin(2) * t400 - pkin(1);
t372 = t386 * qJD(1);
t340 = -pkin(3) * t353 + t372;
t478 = t328 * t340;
t507 = -t254 - t478;
t416 = t353 * t394 - t398 * t354;
t393 = sin(qJ(5));
t453 = qJD(5) * t393;
t454 = qJD(4) * t398;
t281 = t398 * t333 - t394 * t334 + t353 * t454 + t354 * t455;
t389 = qJD(4) + t390;
t464 = t397 * t281 + t389 * t452;
t269 = -t416 * t453 + t464;
t267 = t269 * t393;
t268 = t269 * t397;
t313 = t389 * t393 + t397 * t416;
t484 = t281 * t393;
t270 = t313 * qJD(5) + t484;
t282 = t416 * qJD(4) + t333 * t394 + t398 * t334;
t278 = t397 * t282;
t476 = t416 * t393;
t311 = -t397 * t389 + t476;
t450 = -qJD(5) + t328;
t276 = t393 * t282;
t465 = -t450 * t452 + t276;
t504 = t450 * t393;
t506 = -t282 * MDP(21) - t328 ^ 2 * MDP(19) + (-t328 * t389 + t281) * MDP(20) + (t450 * t505 + t465) * MDP(27) + (-t328 * MDP(18) + t416 * MDP(19) + t389 * MDP(21) - t313 * MDP(27) + MDP(29) * t450) * t416 + (t508 * t313 + t267) * MDP(25) + (t311 * t416 - t450 * t504 + t278) * MDP(28) + (-t393 * t270 - t508 * t311 + t313 * t504 + t268) * MDP(26);
t482 = t310 * t394;
t283 = t305 * t398 - t482;
t279 = -pkin(4) * t389 - t283;
t486 = t279 * t328;
t481 = t310 * t398;
t284 = t305 * t394 + t481;
t439 = t286 * t394 - t398 * t287;
t255 = t284 * qJD(4) + t439;
t480 = t416 * t340;
t503 = -t255 - t480;
t302 = pkin(4) * t416 - pkin(9) * t328;
t500 = -0.2e1 * t449;
t498 = MDP(4) * t396;
t496 = (t396 ^ 2 - t400 ^ 2) * MDP(5);
t495 = qJD(1) * t366;
t280 = pkin(9) * t389 + t284;
t290 = -pkin(4) * t328 - pkin(9) * t416 + t340;
t418 = t280 * t393 - t290 * t397;
t413 = -t255 * t397 + t279 * t453 + t416 * t418;
t265 = t280 * t397 + t290 * t393;
t421 = t255 * t393 + t265 * t416 + t279 * t452;
t368 = t396 * t447;
t370 = t400 * t447;
t474 = t373 * t399;
t410 = -qJD(3) * t474 - t399 * t368 - t395 * t370 - t374 * t456;
t297 = -pkin(8) * t339 + t410;
t365 = t395 * t396 - t399 * t400;
t338 = t390 * t365;
t414 = t373 * t395 - t374 * t399;
t405 = t414 * qJD(3) + t368 * t395 - t399 * t370;
t298 = pkin(8) * t338 + t405;
t319 = -pkin(8) * t366 - t374 * t395 - t474;
t320 = -pkin(8) * t365 - t414;
t417 = t319 * t398 - t320 * t394;
t261 = t417 * qJD(4) + t297 * t398 + t298 * t394;
t336 = t398 * t365 + t366 * t394;
t294 = -t336 * qJD(4) - t338 * t398 - t339 * t394;
t337 = -t365 * t394 + t366 * t398;
t343 = pkin(3) * t365 + t386;
t299 = pkin(4) * t336 - pkin(9) * t337 + t343;
t301 = t319 * t394 + t320 * t398;
t492 = t255 * t337 + t279 * t294 - t301 * t282 + (qJD(5) * t299 + t261) * t450 - (qJD(5) * t290 + t254) * t336;
t489 = pkin(3) * t354;
t485 = t279 * t337;
t483 = t299 * t282;
t475 = t372 * t354;
t473 = t394 * t395;
t472 = t395 * t398;
t401 = qJD(2) ^ 2;
t471 = t396 * t401;
t469 = t400 * t401;
t402 = qJD(1) ^ 2;
t468 = t400 * t402;
t433 = t367 * t395 - t359;
t315 = t433 - t488;
t461 = -t399 * t367 - t355;
t316 = t350 + t461;
t385 = pkin(2) * t399 + pkin(3);
t463 = t315 * t394 + t316 * t398 - t385 * t454 - (-t395 * t455 + (t398 * t399 - t473) * qJD(3)) * pkin(2);
t462 = t315 * t398 - t316 * t394 + t385 * t455 + (t395 * t454 + (t394 * t399 + t472) * qJD(3)) * pkin(2);
t388 = t396 * t487;
t387 = pkin(2) * t458;
t443 = -pkin(2) * t390 - t361;
t442 = -pkin(3) * t389 - t305;
t317 = pkin(3) * t334 + qJD(2) * t387;
t330 = pkin(3) * t339 + t388;
t441 = pkin(1) * t500;
t296 = t302 - t489;
t349 = pkin(2) * t472 + t385 * t394 + pkin(9);
t426 = qJD(5) * t349 + t296 + t387;
t383 = pkin(3) * t394 + pkin(9);
t425 = qJD(5) * t383 + t296;
t288 = t309 * t394 + t481;
t423 = pkin(3) * t455 - t288;
t289 = t309 * t398 - t482;
t422 = -pkin(3) * t454 + t289;
t420 = -t282 * t349 - t486;
t419 = -t282 * t383 - t486;
t412 = -t372 * t353 - t434;
t411 = t294 * t397 - t337 * t453;
t403 = t354 * t353 * MDP(11) + t333 * MDP(13) + (-t353 ^ 2 + t354 ^ 2) * MDP(12) + (-t353 * MDP(13) + (-t354 - t495) * MDP(14)) * t390 + t506;
t384 = -pkin(3) * t398 - pkin(4);
t348 = pkin(2) * t473 - t385 * t398 - pkin(4);
t341 = t387 - t489;
t295 = t337 * qJD(4) - t338 * t394 + t398 * t339;
t263 = pkin(4) * t295 - pkin(9) * t294 + t330;
t262 = t301 * qJD(4) + t297 * t394 - t298 * t398;
t260 = pkin(4) * t282 - pkin(9) * t281 + t317;
t259 = t397 * t260;
t1 = [MDP(6) * t469 + ((-t311 * t397 - t313 * t393) * t294 + (-t267 - t270 * t397 + (t311 * t393 - t313 * t397) * qJD(5)) * t337) * MDP(26) + (t262 * t313 - t265 * t295 - t417 * t269 + ((-qJD(5) * t301 + t263) * t450 - t483 - (-qJD(5) * t280 + t260) * t336 - qJD(5) * t485) * t393 + t492 * t397) * MDP(31) + (t259 * t336 + t262 * t311 - t418 * t295 - t417 * t270 + (-t263 * t450 + t483 + (-t280 * t336 + t301 * t450 + t485) * qJD(5)) * t397 + t492 * t393) * MDP(30) + (t337 * t268 + t411 * t313) * MDP(25) + (-t337 * t276 - t270 * t336 - t295 * t311 - (-t294 * t393 - t337 * t452) * t450) * MDP(28) + (t269 * t336 + t337 * t278 + t295 * t313 - t411 * t450) * MDP(27) - MDP(7) * t471 + (pkin(6) * t471 + t400 * t441) * MDP(10) + (-pkin(6) * t469 + t396 * t441) * MDP(9) + (t386 * t333 - t372 * t338 + (-t354 + t495) * t388) * MDP(17) + (t386 * t334 + t372 * t339 + (qJD(1) * t365 - t353) * t388) * MDP(16) + 0.2e1 * t444 * t498 + t496 * t500 + (t282 * t336 - t295 * t450) * MDP(29) + (t281 * t337 + t294 * t416) * MDP(18) + (-t281 * t336 - t282 * t337 + t294 * t328 - t295 * t416) * MDP(19) + (t333 * t366 + t338 * t354) * MDP(11) + (-t333 * t365 - t334 * t366 - t338 * t353 + t339 * t354) * MDP(12) + (t282 * t343 + t295 * t340 + t317 * t336 - t328 * t330) * MDP(23) + (t281 * t343 + t294 * t340 + t317 * t337 + t330 * t416) * MDP(24) + (-t338 * MDP(13) - t339 * MDP(14) + t405 * MDP(16) - t410 * MDP(17)) * t390 + (t294 * MDP(20) - t295 * MDP(21) - t262 * MDP(23) - t261 * MDP(24)) * t389; (t348 * t270 + t420 * t393 + t462 * t311 - (t463 * t393 - t426 * t397) * t450 + t413) * MDP(30) + (t348 * t269 + t420 * t397 + t462 * t313 - (t426 * t393 + t463 * t397) * t450 + t421) * MDP(31) + (t328 * t341 - t462 * t389 + t503) * MDP(23) + (-t341 * t416 + t463 * t389 + t507) * MDP(24) + (t354 * t387 + t461 * t390 + (t443 * qJD(3) + t362) * t399 + t412) * MDP(17) + (t353 * t387 + t475 - t433 * t390 + (t443 * t395 - t359) * qJD(3) + t435) * MDP(16) + t403 - t468 * t498 + t402 * t496 + (t402 * t396 * MDP(9) + MDP(10) * t468) * pkin(1); (t384 * t270 + t419 * t393 + t423 * t311 - (t422 * t393 - t425 * t397) * t450 + t413) * MDP(30) + (t384 * t269 + t419 * t397 + t423 * t313 - (t425 * t393 + t422 * t397) * t450 + t421) * MDP(31) + (-t328 * t489 + t288 * t389 - t480 + (t442 * t394 - t481) * qJD(4) - t439) * MDP(23) + (-t415 * t390 + t406 + t475) * MDP(16) + (t436 * t390 + t412 - t497) * MDP(17) + (t416 * t489 + t289 * t389 - t478 + (t442 * qJD(4) - t286) * t398 - t438) * MDP(24) + t403; (t284 * t389 + t503) * MDP(23) + (t283 * t389 + t507) * MDP(24) + (-pkin(4) * t270 + (-t283 * t393 + t302 * t397) * t450 - t284 * t311 - t393 * t486 - t465 * pkin(9) + t413) * MDP(30) + (-pkin(4) * t269 - (t283 * t397 + t302 * t393) * t450 - t284 * t313 - t279 * t505 + (-t450 * t453 - t278) * pkin(9) + t421) * MDP(31) + t506; t313 * t311 * MDP(25) + (-t311 ^ 2 + t313 ^ 2) * MDP(26) + (-t311 * t450 + t464) * MDP(27) + (-t313 * t450 - t484) * MDP(28) + t282 * MDP(29) + (-t254 * t393 - t265 * t450 - t279 * t313 + t259) * MDP(30) + (-t254 * t397 - t260 * t393 + t279 * t311 + t418 * t450) * MDP(31) + (-MDP(27) * t476 - t313 * MDP(28) - t265 * MDP(30) + t418 * MDP(31)) * qJD(5);];
tauc = t1;
