% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:12
% EndTime: 2019-03-09 08:48:19
% DurationCPUTime: 5.08s
% Computational Cost: add. (4058->379), mult. (10225->494), div. (0->0), fcn. (7693->8), ass. (0->170)
t440 = sin(qJ(6));
t483 = qJD(6) * t440;
t437 = sin(pkin(10));
t438 = cos(pkin(10));
t442 = sin(qJ(2));
t445 = cos(qJ(2));
t411 = t437 * t445 + t438 * t442;
t397 = t411 * qJD(2);
t380 = qJD(1) * t397;
t478 = qJD(1) * qJD(2);
t475 = t442 * t478;
t420 = t437 * t475;
t474 = t445 * t478;
t381 = t438 * t474 - t420;
t487 = qJD(1) * t445;
t476 = t438 * t487;
t488 = qJD(1) * t442;
t395 = t437 * t488 - t476;
t398 = t411 * qJD(1);
t441 = sin(qJ(5));
t444 = cos(qJ(5));
t484 = qJD(5) * t444;
t485 = qJD(5) * t441;
t300 = t441 * t380 + t444 * t381 + t395 * t484 - t398 * t485;
t432 = qJD(2) - qJD(5);
t443 = cos(qJ(6));
t482 = qJD(6) * t443;
t492 = t443 * t300 - t432 * t482;
t525 = t395 * t441 + t444 * t398;
t286 = -t483 * t525 + t492;
t330 = -t432 * t440 + t443 * t525;
t514 = t300 * t440;
t287 = qJD(6) * t330 + t514;
t301 = qJD(5) * t525 - t444 * t380 + t381 * t441;
t497 = t443 * t301;
t527 = -t444 * t395 + t398 * t441;
t535 = qJD(6) + t527;
t470 = t483 * t535 - t497;
t461 = t440 * t527 * t535 + t470;
t506 = t527 * t432;
t507 = t525 * t432;
t509 = t330 * t525;
t510 = t330 * t535;
t505 = t525 * t440;
t327 = t443 * t432 + t505;
t511 = t327 * t525;
t512 = t327 * t535;
t517 = t286 * t440;
t499 = t440 * t301;
t545 = t535 * t443;
t541 = t535 * t545 + t499;
t553 = -((t287 + t510) * t440 - (t286 - t512) * t443) * MDP(25) + (t330 * t545 + t517) * MDP(24) + (-t509 + t541) * MDP(26) + (t300 - t506) * MDP(19) - (t301 + t507) * MDP(20) - (t461 - t511) * MDP(27) + (t525 ^ 2 - t527 ^ 2) * MDP(18) + (MDP(17) * t527 - MDP(28) * t535) * t525;
t518 = -qJ(3) - pkin(7);
t418 = t518 * t442;
t415 = qJD(1) * t418;
t419 = t518 * t445;
t416 = qJD(1) * t419;
t502 = t437 * t416;
t355 = t438 * t415 + t502;
t480 = qJD(4) - t355;
t539 = pkin(5) * t525;
t417 = -qJD(1) * pkin(1) - pkin(2) * t487 + qJD(3);
t334 = t395 * pkin(3) - t398 * qJ(4) + t417;
t315 = -pkin(4) * t395 - t334;
t285 = pkin(5) * t527 - pkin(9) * t525 + t315;
t407 = qJD(2) * pkin(2) + t415;
t350 = t407 * t438 + t502;
t460 = qJD(4) - t350;
t520 = pkin(8) * t398;
t317 = -t520 + (-pkin(3) - pkin(4)) * qJD(2) + t460;
t501 = t438 * t416;
t351 = t437 * t407 - t501;
t348 = qJD(2) * qJ(4) + t351;
t521 = pkin(8) * t395;
t323 = t348 + t521;
t292 = t317 * t441 + t323 * t444;
t290 = -pkin(9) * t432 + t292;
t276 = t285 * t443 - t290 * t440;
t538 = t276 * t525;
t277 = t285 * t440 + t290 * t443;
t537 = t277 * t525;
t536 = -t520 + t480;
t473 = qJD(2) * t518;
t390 = qJD(3) * t445 + t442 * t473;
t368 = t390 * qJD(1);
t391 = -qJD(3) * t442 + t445 * t473;
t369 = t391 * qJD(1);
t329 = t438 * t368 + t437 * t369;
t433 = qJD(2) * qJD(4);
t324 = t433 + t329;
t309 = pkin(8) * t380 + t324;
t326 = t368 * t437 - t438 * t369;
t311 = -pkin(8) * t381 + t326;
t278 = t444 * t309 + t441 * t311 + t317 * t484 - t323 * t485;
t533 = t315 * t527 - t278;
t279 = t441 * t309 - t311 * t444 + t317 * t485 + t323 * t484;
t532 = -t315 * t525 - t279;
t530 = -0.2e1 * t478;
t529 = MDP(4) * t442;
t528 = MDP(5) * (t442 ^ 2 - t445 ^ 2);
t336 = pkin(2) * t488 + t398 * pkin(3) + t395 * qJ(4);
t318 = -pkin(4) * t398 - t336;
t429 = -pkin(2) * t438 - pkin(3);
t424 = -pkin(4) + t429;
t427 = pkin(2) * t437 + qJ(4);
t455 = t424 * t441 + t427 * t444;
t367 = -pkin(9) + t455;
t524 = t535 * (-pkin(9) * t527 + qJD(6) * t367 + t318 - t539) - t279;
t523 = t535 * (pkin(9) * t535 + t539) + t279;
t337 = t390 * t437 - t438 * t391;
t486 = qJD(2) * t442;
t500 = t438 * t445;
t400 = qJD(2) * t500 - t437 * t486;
t319 = -pkin(8) * t400 + t337;
t338 = t438 * t390 + t437 * t391;
t320 = pkin(8) * t397 + t338;
t358 = -t438 * t418 - t419 * t437;
t339 = -pkin(8) * t411 + t358;
t359 = t437 * t418 - t438 * t419;
t410 = t437 * t442 - t500;
t340 = pkin(8) * t410 + t359;
t459 = t339 * t444 - t340 * t441;
t283 = qJD(5) * t459 + t319 * t441 + t320 * t444;
t291 = t317 * t444 - t323 * t441;
t289 = pkin(5) * t432 - t291;
t519 = -t445 * pkin(2) - pkin(1);
t349 = t410 * pkin(3) - t411 * qJ(4) + t519;
t325 = -pkin(4) * t410 - t349;
t353 = t410 * t441 + t411 * t444;
t457 = t444 * t410 - t411 * t441;
t293 = -pkin(5) * t457 - pkin(9) * t353 + t325;
t297 = t339 * t441 + t340 * t444;
t306 = qJD(5) * t457 + t397 * t441 + t400 * t444;
t522 = t279 * t353 + t289 * t306 - t297 * t301 - (qJD(6) * t293 + t283) * t535 + (qJD(6) * t285 + t278) * t457;
t393 = t398 ^ 2;
t516 = t289 * t353;
t515 = t293 * t301;
t513 = t326 * t358;
t446 = qJD(2) ^ 2;
t498 = t442 * t446;
t495 = t445 * t446;
t447 = qJD(1) ^ 2;
t494 = t445 * t447;
t354 = t415 * t437 - t501;
t331 = t354 + t521;
t456 = t424 * t444 - t427 * t441;
t493 = -qJD(5) * t456 + t331 * t441 - t444 * t536;
t491 = qJD(5) * t455 + t331 * t444 + t441 * t536;
t477 = pkin(2) * t486;
t471 = pkin(1) * t530;
t463 = t432 ^ 2;
t462 = t444 * t432;
t425 = pkin(2) * t475;
t312 = t380 * pkin(3) - t381 * qJ(4) - t398 * qJD(4) + t425;
t454 = t334 * t398 + t326;
t453 = t306 * t443 - t353 * t483;
t299 = -pkin(4) * t380 - t312;
t322 = t397 * pkin(3) - t400 * qJ(4) - t411 * qJD(4) + t477;
t308 = -pkin(4) * t397 - t322;
t451 = -pkin(9) * t301 + (t289 + t291) * t535;
t450 = t326 * t411 + t337 * t398 - t338 * t395 + t358 * t381 - t359 * t380;
t448 = -t367 * t301 + (-t289 + t493) * t535;
t366 = pkin(5) - t456;
t342 = -qJD(2) * pkin(3) + t460;
t307 = qJD(5) * t353 - t444 * t397 + t400 * t441;
t284 = qJD(5) * t297 - t319 * t444 + t320 * t441;
t282 = pkin(5) * t307 - pkin(9) * t306 + t308;
t281 = pkin(5) * t301 - pkin(9) * t300 + t299;
t280 = t443 * t281;
t1 = [(t299 * t353 + t300 * t325 + t306 * t315 + t308 * t525) * MDP(23) + (-t299 * t457 + t301 * t325 + t307 * t315 + t308 * t527) * MDP(22) + (t286 * t353 * t443 + t330 * t453) * MDP(24) + (-t301 * t457 + t307 * t535) * MDP(28) + (t300 * t457 - t301 * t353 - t306 * t527 - t307 * t525) * MDP(18) + (t300 * t353 + t306 * t525) * MDP(17) + 0.2e1 * t474 * t529 + (t312 * t349 + t322 * t334 + t324 * t359 + t337 * t342 + t338 * t348 + t513) * MDP(16) + (-t329 * t410 - t350 * t400 - t351 * t397 + t450) * MDP(11) + (-t324 * t410 + t342 * t400 - t348 * t397 + t450) * MDP(14) + t528 * t530 + (-pkin(7) * t495 + t442 * t471) * MDP(9) + (-t286 * t457 + t307 * t330 + t353 * t497 + t453 * t535) * MDP(26) - MDP(7) * t498 + (pkin(7) * t498 + t445 * t471) * MDP(10) + (-t353 * t499 + t287 * t457 - t307 * t327 + (-t306 * t440 - t353 * t482) * t535) * MDP(27) + (-qJD(2) * t337 + t312 * t410 + t322 * t395 + t334 * t397 + t349 * t380) * MDP(13) + (-t277 * t307 + t284 * t330 - t459 * t286 + (-(-qJD(6) * t297 + t282) * t535 - t515 + (-qJD(6) * t290 + t281) * t457 - qJD(6) * t516) * t440 + t522 * t443) * MDP(30) + (t276 * t307 - t280 * t457 + t284 * t327 - t459 * t287 + (t282 * t535 + t515 + (t290 * t457 - t297 * t535 + t516) * qJD(6)) * t443 + t522 * t440) * MDP(29) + ((-t327 * t443 - t330 * t440) * t306 + (-t517 - t287 * t443 + (t327 * t440 - t330 * t443) * qJD(6)) * t353) * MDP(25) + (t513 + t329 * t359 - t350 * t337 + t351 * t338 + (qJD(1) * t519 + t417) * t477) * MDP(12) + MDP(6) * t495 + (qJD(2) * t338 - t312 * t411 - t322 * t398 - t334 * t400 - t349 * t381) * MDP(15) + (-MDP(19) * t306 + MDP(20) * t307 + MDP(22) * t284 + MDP(23) * t283) * t432; (qJD(2) * t354 - t336 * t395 - t454) * MDP(13) + t447 * t528 + (MDP(9) * t442 * t447 + MDP(10) * t494) * pkin(1) + (t324 * t427 + t326 * t429 - t334 * t336 - t342 * t354 + t348 * t480) * MDP(16) + (t366 * t286 + t491 * t330 + t440 * t524 + t448 * t443 - t537) * MDP(30) + (t366 * t287 + t491 * t327 + t448 * t440 - t443 * t524 + t538) * MDP(29) + (-t318 * t527 + t432 * t491 - t532) * MDP(22) + (-t318 * t525 - t432 * t493 - t533) * MDP(23) + (t350 * t354 - t351 * t355 + (-t326 * t438 + t329 * t437 - t417 * t488) * pkin(2)) * MDP(12) + ((t351 - t354) * t398 + (-t350 + t355) * t395 + (-t380 * t437 - t381 * t438) * pkin(2)) * MDP(11) + (-t380 * t427 + t381 * t429 + (t348 - t354) * t398 + (t342 - t480) * t395) * MDP(14) + (-qJD(2) * t355 - t334 * t395 + t336 * t398 + t329 + 0.2e1 * t433) * MDP(15) - t494 * t529 - t553; (t350 * t398 + t351 * t395 + t425) * MDP(12) + t420 * MDP(15) + (-t342 * t398 + t348 * t395 + t312) * MDP(16) + (-t301 + t507) * MDP(22) + (-t300 - t506) * MDP(23) + (t461 + t511) * MDP(29) + (t509 + t541) * MDP(30) + ((t437 * t487 + t438 * t488 + t398) * MDP(13) + (t395 - t476) * MDP(15)) * qJD(2) + (MDP(11) + MDP(14)) * (-t395 ^ 2 - t393); t398 * t395 * MDP(13) + (-t420 + (t395 + t476) * qJD(2)) * MDP(14) + (-t393 - t446) * MDP(15) + (-qJD(2) * t348 + t454) * MDP(16) + (-t398 * t527 - t441 * t463) * MDP(22) + (-t398 * t525 - t444 * t463) * MDP(23) + (-t444 * t287 + (-t398 * t443 + t440 * t462) * t535 + (-t327 * t432 - t482 * t535 - t499) * t441) * MDP(29) + (-t444 * t286 + (t398 * t440 + t443 * t462) * t535 + (-t330 * t432 + t470) * t441) * MDP(30); (-t292 * t432 + t532) * MDP(22) + (-t291 * t432 + t533) * MDP(23) + (-pkin(5) * t287 - t292 * t327 + t451 * t440 - t443 * t523 - t538) * MDP(29) + (-pkin(5) * t286 - t292 * t330 + t440 * t523 + t451 * t443 + t537) * MDP(30) + t553; t330 * t327 * MDP(24) + (-t327 ^ 2 + t330 ^ 2) * MDP(25) + (t492 + t512) * MDP(26) + (t510 - t514) * MDP(27) + t301 * MDP(28) + (t277 * t535 - t278 * t440 - t289 * t330 + t280) * MDP(29) + (t276 * t535 - t278 * t443 - t281 * t440 + t289 * t327) * MDP(30) + (-MDP(26) * t505 - MDP(27) * t330 - MDP(29) * t277 - MDP(30) * t276) * qJD(6);];
tauc  = t1;
