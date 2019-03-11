% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:30:06
% EndTime: 2019-03-09 04:30:11
% DurationCPUTime: 3.96s
% Computational Cost: add. (4457->411), mult. (10372->545), div. (0->0), fcn. (6629->8), ass. (0->176)
t416 = sin(qJ(4));
t418 = cos(qJ(4));
t461 = t418 * qJD(3);
t417 = sin(qJ(3));
t471 = qJD(1) * t417;
t375 = t416 * t471 - t461;
t468 = qJD(3) * t416;
t377 = t418 * t471 + t468;
t412 = sin(pkin(10));
t414 = cos(pkin(10));
t326 = t414 * t375 + t377 * t412;
t419 = cos(qJ(3));
t470 = qJD(1) * t419;
t398 = -qJD(4) + t470;
t515 = t326 * t398;
t371 = t412 * t418 + t414 * t416;
t360 = t371 * qJD(4);
t479 = t371 * t470 - t360;
t465 = qJD(4) * t416;
t446 = t412 * t465;
t451 = t416 * t470;
t464 = qJD(4) * t418;
t493 = t414 * t418;
t478 = -t412 * t451 - t414 * t464 + t470 * t493 + t446;
t447 = t417 * t464;
t466 = qJD(3) * t419;
t450 = t416 * t466;
t509 = t447 + t450;
t431 = -t375 * t412 + t414 * t377;
t514 = t431 ^ 2;
t513 = MDP(5) * t417;
t410 = t417 ^ 2;
t512 = MDP(6) * (-t419 ^ 2 + t410);
t403 = sin(pkin(9)) * pkin(1) + pkin(7);
t387 = t403 * qJD(1);
t510 = qJD(2) * t419 - t417 * t387;
t339 = -qJD(3) * pkin(3) - t510;
t323 = pkin(4) * t375 + qJD(5) + t339;
t281 = pkin(5) * t326 - qJ(6) * t431 + t323;
t511 = t281 * t431;
t487 = t418 * t419;
t427 = pkin(4) * t417 - qJ(5) * t487;
t434 = pkin(3) * t417 - pkin(8) * t419;
t379 = t434 * qJD(1);
t437 = t418 * t379 - t416 * t510;
t309 = qJD(1) * t427 + t437;
t480 = t416 * t379 + t418 * t510;
t314 = -qJ(5) * t451 + t480;
t506 = -qJ(5) - pkin(8);
t441 = qJD(4) * t506;
t463 = qJD(5) * t418;
t357 = t416 * t441 + t463;
t423 = -qJD(5) * t416 + t418 * t441;
t481 = (-t309 + t423) * t414 + (t314 - t357) * t412;
t408 = t417 * qJD(2);
t351 = t419 * t387 + t408;
t508 = -t351 + (-t451 + t465) * pkin(4);
t507 = MDP(19) + MDP(22);
t340 = qJD(3) * pkin(8) + t351;
t405 = -cos(pkin(9)) * pkin(1) - pkin(2);
t367 = -pkin(3) * t419 - pkin(8) * t417 + t405;
t345 = t367 * qJD(1);
t492 = t416 * t345;
t312 = t340 * t418 + t492;
t303 = -qJ(5) * t375 + t312;
t300 = t414 * t303;
t311 = -t340 * t416 + t418 * t345;
t302 = -qJ(5) * t377 + t311;
t275 = t302 * t412 + t300;
t505 = t275 * t431;
t504 = t303 * t412;
t448 = t417 * t465;
t459 = qJD(1) * qJD(3);
t443 = t419 * t459;
t458 = qJD(3) * qJD(4);
t473 = (t443 + t458) * t418;
t335 = -qJD(1) * t448 + t473;
t503 = t335 * t416;
t442 = t416 * t458;
t336 = qJD(1) * t509 + t442;
t502 = t336 * t419;
t501 = t339 * t416;
t500 = t339 * t418;
t344 = qJD(3) * t408 + t387 * t466;
t499 = t344 * t416;
t498 = t344 * t418;
t497 = t375 * t398;
t496 = t377 * t398;
t495 = t398 * t418;
t494 = t403 * t416;
t491 = t416 * t417;
t490 = t416 * t419;
t489 = t417 * t418;
t420 = qJD(3) ^ 2;
t488 = t417 * t420;
t486 = t419 * t420;
t343 = t510 * qJD(3);
t380 = t434 * qJD(3);
t366 = qJD(1) * t380;
t438 = t416 * t343 - t418 * t366;
t422 = -qJD(4) * t312 - t438;
t444 = t417 * t459;
t272 = pkin(4) * t444 - qJ(5) * t335 - qJD(5) * t377 + t422;
t455 = t418 * t343 + t345 * t464 + t416 * t366;
t426 = -t340 * t465 + t455;
t278 = -qJ(5) * t336 - qJD(5) * t375 + t426;
t485 = -t414 * t272 + t412 * t278;
t260 = t412 * t272 + t414 * t278;
t283 = t412 * t309 + t414 * t314;
t279 = qJ(6) * t471 + t283;
t320 = t414 * t357 + t412 * t423;
t484 = t279 - t320;
t483 = pkin(5) * t471 - t481;
t378 = t403 * t487;
t467 = qJD(3) * t417;
t476 = t418 * t380 + t467 * t494;
t286 = -t417 * t463 + t427 * qJD(3) + (-t378 + (qJ(5) * t417 - t367) * t416) * qJD(4) + t476;
t477 = t367 * t464 + t416 * t380;
t291 = (-qJ(5) * qJD(4) - qJD(3) * t403) * t489 + (-qJD(5) * t417 + (-qJ(5) * qJD(3) - qJD(4) * t403) * t419) * t416 + t477;
t265 = t412 * t286 + t414 * t291;
t482 = t479 * pkin(5) - t478 * qJ(6) + qJD(6) * t371 - t508;
t299 = -pkin(4) * t398 + t302;
t271 = t412 * t299 + t300;
t355 = t418 * t367;
t318 = -qJ(5) * t489 + t355 + (-pkin(4) - t494) * t419;
t475 = t416 * t367 + t378;
t322 = -qJ(5) * t491 + t475;
t290 = t412 * t318 + t414 * t322;
t474 = pkin(4) * t491 + t417 * t403;
t388 = qJD(1) * t405;
t462 = qJD(6) * t398;
t276 = t302 * t414 - t504;
t460 = qJD(6) - t276;
t456 = qJ(6) * t444 + t260;
t453 = t509 * pkin(4) + t403 * t466;
t452 = -pkin(4) * t418 - pkin(3);
t449 = t419 * t461;
t445 = t506 * t416;
t306 = t335 * t412 + t414 * t336;
t440 = -t335 * t419 + t377 * t467;
t439 = t398 * t403 + t340;
t436 = t398 * t448;
t435 = t398 * t447;
t313 = pkin(4) * t336 + t344;
t307 = t335 * t414 - t336 * t412;
t389 = t506 * t418;
t331 = -t389 * t412 - t414 * t445;
t332 = -t414 * t389 + t412 * t445;
t433 = -t332 * t306 + t307 * t331 - t320 * t326;
t264 = t286 * t414 - t291 * t412;
t270 = t299 * t414 - t504;
t289 = t318 * t414 - t322 * t412;
t430 = qJD(1) * t410 - t398 * t419;
t428 = 0.2e1 * qJD(3) * t388;
t258 = -pkin(5) * t444 + t485;
t425 = t430 * t416;
t261 = pkin(5) * t306 - qJ(6) * t307 - qJD(6) * t431 + t313;
t404 = -pkin(4) * t414 - pkin(5);
t401 = pkin(4) * t412 + qJ(6);
t370 = t412 * t416 - t493;
t349 = -t412 * t491 + t414 * t489;
t348 = t371 * t417;
t324 = pkin(5) * t370 - qJ(6) * t371 + t452;
t317 = t360 * t417 + t412 * t450 - t414 * t449;
t316 = -t412 * t449 - t414 * t509 + t417 * t446;
t304 = pkin(5) * t348 - qJ(6) * t349 + t474;
t292 = pkin(4) * t377 + pkin(5) * t431 + qJ(6) * t326;
t287 = pkin(5) * t419 - t289;
t285 = -qJ(6) * t419 + t290;
t277 = -pkin(5) * t316 + qJ(6) * t317 - qJD(6) * t349 + t453;
t267 = -qJ(6) * t398 + t271;
t266 = pkin(5) * t398 + qJD(6) - t270;
t263 = -pkin(5) * t467 - t264;
t262 = qJ(6) * t467 - qJD(6) * t419 + t265;
t257 = t456 - t462;
t1 = [0.2e1 * t443 * t513 - 0.2e1 * t459 * t512 + MDP(7) * t486 - MDP(8) * t488 + (-t403 * t486 + t417 * t428) * MDP(10) + (t403 * t488 + t419 * t428) * MDP(11) + (t335 * t489 + (-t448 + t449) * t377) * MDP(12) + ((-t375 * t418 - t377 * t416) * t466 + (-t503 - t336 * t418 + (t375 * t416 - t377 * t418) * qJD(4)) * t417) * MDP(13) + (t430 * t461 + t436 + t440) * MDP(14) + (t435 + t502 + (-t375 * t417 - t425) * qJD(3)) * MDP(15) + (-t398 - t470) * MDP(16) * t467 + (-(-t367 * t465 + t476) * t398 + ((t375 * t403 + t501) * qJD(3) + (t418 * t439 + t492) * qJD(4) + t438) * t419 + (t339 * t464 + t403 * t336 + t499 + ((-t403 * t490 + t355) * qJD(1) + t311) * qJD(3)) * t417) * MDP(17) + (t477 * t398 + (-t439 * t465 + (t377 * t403 + t500) * qJD(3) + t455) * t419 + (-t339 * t465 + t403 * t335 + t498 + (-qJD(1) * t475 - t403 * t495 - t312) * qJD(3)) * t417) * MDP(18) + (-t260 * t348 - t264 * t431 - t265 * t326 + t270 * t317 + t271 * t316 - t289 * t307 - t290 * t306 + t349 * t485) * MDP(19) + (t260 * t290 + t270 * t264 + t271 * t265 - t289 * t485 + t313 * t474 + t323 * t453) * MDP(20) + (t258 * t419 + t261 * t348 + t263 * t398 + t277 * t326 - t281 * t316 + t304 * t306 + (-qJD(1) * t287 - t266) * t467) * MDP(21) + (-t257 * t348 + t258 * t349 - t262 * t326 + t263 * t431 - t266 * t317 + t267 * t316 - t285 * t306 + t287 * t307) * MDP(22) + (-t257 * t419 - t261 * t349 - t262 * t398 - t277 * t431 + t281 * t317 - t304 * t307 + (qJD(1) * t285 + t267) * t467) * MDP(23) + (t257 * t285 + t258 * t287 + t261 * t304 + t262 * t267 + t263 * t266 + t277 * t281) * MDP(24); (t435 - t502) * MDP(17) + (-t436 + t440) * MDP(18) + (t260 * t349 + t270 * t316 - t271 * t317 - t313 * t419 + t348 * t485) * MDP(20) + (-t306 * t419 - t316 * t398) * MDP(21) + (t307 * t419 + t317 * t398) * MDP(23) + (t257 * t349 + t258 * t348 - t261 * t419 - t266 * t316 - t267 * t317) * MDP(24) + (-MDP(10) * t417 - MDP(11) * t419) * t420 + (-t430 * MDP(18) * t418 - MDP(17) * t425 + (t375 * MDP(17) + t323 * MDP(20) + (-qJD(1) * t348 + t326) * MDP(21) + (qJD(1) * t349 - t431) * MDP(23) + t281 * MDP(24)) * t417) * qJD(3) + t507 * (-t349 * t306 + t307 * t348 - t316 * t431 + t317 * t326); (qJD(3) * t351 - t388 * t471 - t344) * MDP(10) - t388 * t470 * MDP(11) + (-t377 * t495 + t503) * MDP(12) + ((t335 + t497) * t418 + (-t336 + t496) * t416) * MDP(13) + (-t398 * t464 + (t398 * t487 + (-t377 + t468) * t417) * qJD(1)) * MDP(14) + (t398 * t465 + (-t398 * t490 + (t375 + t461) * t417) * qJD(1)) * MDP(15) + t398 * MDP(16) * t471 + (-pkin(3) * t336 - t498 + t437 * t398 - t351 * t375 + (pkin(8) * t495 + t501) * qJD(4) + (-t311 * t417 + (-pkin(8) * t467 - t339 * t419) * t416) * qJD(1)) * MDP(17) + (-pkin(3) * t335 + t499 - t480 * t398 - t351 * t377 + (-pkin(8) * t398 * t416 + t500) * qJD(4) + (-t339 * t487 + (-pkin(8) * t461 + t312) * t417) * qJD(1)) * MDP(18) + (-t260 * t370 + t270 * t478 + t271 * t479 + t283 * t326 + t371 * t485 - t431 * t481 + t433) * MDP(19) + (t260 * t332 + t485 * t331 + t313 * t452 + t508 * t323 + (t320 - t283) * t271 + t481 * t270) * MDP(20) + (t261 * t370 + t306 * t324 + t483 * t398 - t482 * t326 - t479 * t281 + (-qJD(3) * t331 + t266) * t471) * MDP(21) + (-t257 * t370 + t258 * t371 - t266 * t478 + t267 * t479 + t279 * t326 + t431 * t483 + t433) * MDP(22) + (-t261 * t371 - t307 * t324 + t484 * t398 + t482 * t431 + t478 * t281 + (qJD(3) * t332 - t267) * t471) * MDP(23) + (t257 * t332 + t258 * t331 + t261 * t324 + t266 * t483 - t267 * t484 - t281 * t482) * MDP(24) + (-t419 * t513 + t512) * qJD(1) ^ 2; t377 * t375 * MDP(12) + (-t375 ^ 2 + t377 ^ 2) * MDP(13) + (t473 - t497) * MDP(14) + (-t442 - t496) * MDP(15) + (-t312 * t398 - t339 * t377 + t422) * MDP(17) + (-t311 * t398 + t339 * t375 - t426) * MDP(18) + (t271 * t431 - t505) * MDP(19) + (t270 * t275 - t271 * t276) * MDP(20) + (-t275 * t398 - t485 - t511) * MDP(21) + (t267 * t431 - t306 * t401 + t307 * t404 - t505) * MDP(22) + (t276 * t398 + t292 * t431 + t456 - 0.2e1 * t462) * MDP(23) + (t257 * t401 + t258 * t404 - t266 * t275 + t267 * t460 - t281 * t292) * MDP(24) + (-MDP(15) * t450 + ((-MDP(14) * t416 - MDP(15) * t418) * qJD(4) + (MDP(16) + (pkin(5) - t404) * MDP(21) + t401 * MDP(23)) * qJD(3)) * t417) * qJD(1) + ((-t306 * t412 - t307 * t414) * MDP(19) + (t260 * t412 - t323 * t377 - t414 * t485) * MDP(20)) * pkin(4) + ((-t270 + t276) * MDP(19) - t292 * MDP(21) + (t266 - t460) * MDP(22) - t281 * MDP(23)) * t326; (t270 * t431 + t271 * t326 + t313) * MDP(20) + (-t398 * t431 + t306) * MDP(21) + (-t307 - t515) * MDP(23) + (-t266 * t431 + t267 * t326 + t261) * MDP(24) + t507 * (-t326 ^ 2 - t514); (t326 * t431 - t444) * MDP(21) + (t307 - t515) * MDP(22) + (-t398 ^ 2 - t514) * MDP(23) + (t267 * t398 + t258 + t511) * MDP(24);];
tauc  = t1;
