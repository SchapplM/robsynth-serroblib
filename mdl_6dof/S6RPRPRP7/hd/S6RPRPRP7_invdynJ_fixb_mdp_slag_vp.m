% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:54
% EndTime: 2019-03-09 03:23:01
% DurationCPUTime: 4.25s
% Computational Cost: add. (3889->426), mult. (7663->535), div. (0->0), fcn. (5116->10), ass. (0->190)
t421 = -pkin(1) - pkin(7);
t371 = t421 * qJDD(1) + qJDD(2);
t419 = cos(qJ(3));
t363 = t419 * t371;
t372 = t421 * qJD(1) + qJD(2);
t416 = sin(qJ(3));
t469 = qJDD(1) * t419;
t471 = qJD(1) * qJD(4);
t472 = qJD(1) * qJD(3);
t480 = qJD(3) * t416;
t302 = -t419 * t471 - t372 * t480 + qJDD(3) * pkin(3) + t363 + (t416 * t472 - t469) * qJ(4);
t479 = qJD(3) * t419;
t311 = (-qJ(4) * qJD(1) + t372) * t479 + (-qJ(4) * qJDD(1) + t371 - t471) * t416;
t516 = sin(pkin(9));
t517 = cos(pkin(9));
t285 = t516 * t302 + t517 * t311;
t281 = qJDD(3) * pkin(8) + t285;
t481 = qJD(1) * t419;
t347 = -qJ(4) * t481 + t419 * t372;
t343 = qJD(3) * pkin(3) + t347;
t482 = qJD(1) * t416;
t346 = -qJ(4) * t482 + t372 * t416;
t460 = t517 * t346;
t305 = t516 * t343 + t460;
t297 = qJD(3) * pkin(8) + t305;
t428 = t517 * t416 + t516 * t419;
t353 = t428 * qJD(1);
t361 = -t516 * t416 + t517 * t419;
t356 = t361 * qJD(1);
t365 = pkin(3) * t482 + qJD(1) * qJ(2) + qJD(4);
t306 = pkin(4) * t353 - pkin(8) * t356 + t365;
t415 = sin(qJ(5));
t418 = cos(qJ(5));
t283 = t297 * t418 + t306 * t415;
t457 = qJD(3) * t517;
t536 = qJD(1) * t457 + qJDD(1) * t516;
t456 = qJD(3) * t516;
t537 = qJD(1) * t456 - qJDD(1) * t517;
t326 = t537 * t416 - t536 * t419;
t427 = t536 * t416 + t537 * t419;
t408 = qJDD(1) * qJ(2);
t409 = qJD(1) * qJD(2);
t462 = t419 * t472;
t470 = qJDD(1) * t416;
t440 = qJDD(4) + t408 + t409 + (t462 + t470) * pkin(3);
t288 = -t326 * pkin(4) + t427 * pkin(8) + t440;
t287 = t418 * t288;
t475 = t418 * qJD(3);
t477 = qJD(5) * t415;
t292 = -qJD(5) * t475 - t415 * qJDD(3) + t356 * t477 + t418 * t427;
t324 = -qJDD(5) + t326;
t334 = qJD(3) * t415 + t356 * t418;
t268 = -pkin(5) * t324 + qJ(6) * t292 - t283 * qJD(5) - qJD(6) * t334 - t281 * t415 + t287;
t332 = t356 * t415 - t475;
t276 = -qJ(6) * t332 + t283;
t535 = qJD(5) + t353;
t539 = t276 * t535 + t268;
t425 = -t418 * qJDD(3) - t415 * t427;
t293 = t334 * qJD(5) + t425;
t476 = qJD(5) * t418;
t433 = t418 * t281 + t415 * t288 - t297 * t477 + t306 * t476;
t269 = -qJ(6) * t293 - qJD(6) * t332 + t433;
t282 = -t297 * t415 + t418 * t306;
t275 = -qJ(6) * t334 + t282;
t273 = pkin(5) * t535 + t275;
t538 = -t273 * t535 + t269;
t423 = qJD(1) ^ 2;
t417 = sin(qJ(1));
t420 = cos(qJ(1));
t531 = g(1) * t417 - g(2) * t420;
t431 = -qJ(2) * t423 - t531;
t403 = t416 * pkin(3);
t414 = -qJ(4) - pkin(7);
t532 = t420 * t403 + t417 * t414;
t407 = qJ(3) + pkin(9);
t395 = sin(t407);
t497 = t418 * t420;
t500 = t415 * t417;
t349 = -t395 * t500 + t497;
t498 = t417 * t418;
t499 = t415 * t420;
t351 = t395 * t499 + t498;
t530 = -g(1) * t349 - g(2) * t351;
t284 = t517 * t302 - t516 * t311;
t280 = -qJDD(3) * pkin(4) - t284;
t386 = t516 * pkin(3) + pkin(8);
t396 = cos(t407);
t426 = -g(3) * t395 + t396 * t531;
t478 = qJD(5) * t535;
t529 = t386 * t478 + t280 + t426;
t449 = g(1) * t420 + g(2) * t417;
t467 = 0.2e1 * t409;
t528 = 0.2e1 * t408 + t467 - t449;
t444 = t273 * t418 + t276 * t415;
t505 = t334 * t418;
t509 = t332 * t415;
t527 = -MDP(23) * (t505 + t509) + MDP(24) * t444 + t535 * (t418 * MDP(21) - t415 * MDP(22));
t526 = t334 ^ 2;
t524 = pkin(5) * t415;
t520 = g(3) * t396;
t519 = g(3) * t416;
t518 = t418 * pkin(5);
t515 = pkin(1) * qJDD(1);
t513 = t292 * t415;
t512 = t293 * t418;
t511 = t332 * t353;
t510 = t332 * t356;
t508 = t332 * t418;
t507 = t334 * t356;
t506 = t334 * t415;
t504 = t353 * t415;
t503 = t361 * t415;
t502 = t361 * t418;
t501 = t415 * t324;
t319 = t418 * t324;
t495 = qJ(4) - t421;
t366 = t495 * t416;
t367 = t495 * t419;
t329 = -t517 * t366 - t516 * t367;
t327 = t418 * t329;
t496 = qJ(2) + t403;
t494 = qJ(6) + t386;
t493 = -t275 + t273;
t492 = -t415 * t293 - t332 * t476;
t340 = t516 * t346;
t315 = t517 * t347 - t340;
t317 = pkin(3) * t481 + pkin(4) * t356 + pkin(8) * t353;
t491 = t418 * t315 + t415 * t317;
t490 = -t504 * t535 - t319;
t325 = pkin(4) * t428 - pkin(8) * t361 + t496;
t489 = t415 * t325 + t327;
t455 = qJD(5) * t494;
t488 = -qJ(6) * t504 + qJD(6) * t418 - t415 * t455 - t491;
t313 = t418 * t317;
t487 = -pkin(5) * t356 - t313 + (-qJ(6) * t353 - t455) * t418 + (-qJD(6) + t315) * t415;
t486 = t420 * pkin(1) + t417 * qJ(2);
t412 = t419 ^ 2;
t484 = t416 ^ 2 - t412;
t422 = qJD(3) ^ 2;
t483 = -t422 - t423;
t474 = pkin(3) * t479 + qJD(2);
t468 = qJDD(3) * t416;
t344 = -qJD(4) * t419 + t495 * t480;
t345 = -qJD(3) * t367 - qJD(4) * t416;
t310 = t516 * t344 + t517 * t345;
t354 = t416 * t456 - t419 * t457;
t355 = -t416 * t457 - t419 * t456;
t316 = -pkin(4) * t354 - pkin(8) * t355 + t474;
t466 = t418 * t310 + t415 * t316 + t325 * t476;
t465 = t417 * t403 + t486;
t464 = t361 * t476;
t402 = t420 * qJ(2);
t461 = -pkin(1) * t417 + t402;
t452 = t535 * t418;
t451 = -qJD(5) * t306 - t281;
t450 = qJDD(2) - t515;
t387 = -t517 * pkin(3) - pkin(4);
t445 = -t297 * t476 + t287;
t304 = t517 * t343 - t340;
t309 = -t517 * t344 + t516 * t345;
t314 = t516 * t347 + t460;
t328 = -t516 * t366 + t517 * t367;
t443 = t273 * t415 - t276 * t418;
t391 = pkin(4) + t518;
t413 = -qJ(6) - pkin(8);
t441 = t391 * t395 + t396 * t413;
t439 = -qJ(6) * t355 - qJD(6) * t361;
t437 = t415 * MDP(21) + t418 * MDP(22);
t436 = t355 * t415 + t464;
t435 = t355 * t418 - t361 * t477;
t434 = 0.2e1 * qJ(2) * t472 + qJDD(3) * t421;
t296 = -qJD(3) * pkin(4) - t304;
t432 = t296 * t535 + t386 * t324;
t430 = t284 * t361 + t304 * t355 - t531;
t272 = t293 * pkin(5) + qJDD(6) + t280;
t424 = -t421 * t422 + t528;
t399 = qJDD(3) * t419;
t359 = t494 * t418;
t358 = t494 * t415;
t352 = t395 * t497 - t500;
t350 = t395 * t498 + t499;
t331 = t332 ^ 2;
t321 = t418 * t325;
t308 = t418 * t316;
t290 = t332 * pkin(5) + qJD(6) + t296;
t289 = -qJ(6) * t503 + t489;
t279 = pkin(5) * t428 - qJ(6) * t502 - t329 * t415 + t321;
t271 = -qJ(6) * t464 + (-qJD(5) * t329 + t439) * t415 + t466;
t270 = -pkin(5) * t354 - t310 * t415 + t308 + t439 * t418 + (-t327 + (qJ(6) * t361 - t325) * t415) * qJD(5);
t1 = [qJDD(1) * MDP(1) + t531 * MDP(2) + t449 * MDP(3) + (qJDD(2) - t531 - 0.2e1 * t515) * MDP(4) + t528 * MDP(5) + (-t450 * pkin(1) - g(1) * t461 - g(2) * t486 + (t467 + t408) * qJ(2)) * MDP(6) + (qJDD(1) * t412 - 0.2e1 * t416 * t462) * MDP(7) + 0.2e1 * (-t416 * t469 + t472 * t484) * MDP(8) + (-t416 * t422 + t399) * MDP(9) + (-t419 * t422 - t468) * MDP(10) + (t416 * t424 + t419 * t434) * MDP(12) + (-t416 * t434 + t419 * t424) * MDP(13) + (-t285 * t428 + t305 * t354 + t309 * t356 - t310 * t353 + t329 * t326 - t328 * t427 - t430) * MDP(14) + (t285 * t329 + t305 * t310 - t284 * t328 - t304 * t309 + t440 * t496 + t365 * t474 - g(1) * (t461 + t532) - g(2) * (-t414 * t420 + t465)) * MDP(15) + (-t292 * t502 + t334 * t435) * MDP(16) + ((-t506 - t508) * t355 + (t513 - t512 + (-t505 + t509) * qJD(5)) * t361) * MDP(17) + (-t292 * t428 - t319 * t361 - t334 * t354 + t435 * t535) * MDP(18) + (-t293 * t428 + t332 * t354 + t361 * t501 - t436 * t535) * MDP(19) + (-t324 * t428 - t354 * t535) * MDP(20) + ((-t329 * t476 + t308) * t535 - t321 * t324 + t445 * t428 - t282 * t354 + t309 * t332 + t328 * t293 + t296 * t464 - g(1) * t352 - g(2) * t350 + ((-qJD(5) * t325 - t310) * t535 + t329 * t324 + t451 * t428 + t280 * t361 + t296 * t355) * t415) * MDP(21) + (-(-t329 * t477 + t466) * t535 + t489 * t324 - t433 * t428 + t283 * t354 + t309 * t334 - t328 * t292 + t280 * t502 + g(1) * t351 - g(2) * t349 + t435 * t296) * MDP(22) + (-t270 * t334 - t271 * t332 + t279 * t292 - t289 * t293 + t449 * t396 - t444 * t355 + (qJD(5) * t443 - t268 * t418 - t269 * t415) * t361) * MDP(23) + (t269 * t289 + t276 * t271 + t268 * t279 + t273 * t270 + t272 * (pkin(5) * t503 + t328) + t290 * (pkin(5) * t436 + t309) - g(1) * (t402 + t441 * t420 + (-pkin(1) - t524) * t417 + t532) - g(2) * ((-t414 + t524) * t420 + t441 * t417 + t465)) * MDP(24); qJDD(1) * MDP(4) - t423 * MDP(5) + (t450 + t431) * MDP(6) + (t416 * t483 + t399) * MDP(12) + (t419 * t483 - t468) * MDP(13) + (-t355 * t356 + t361 * t427) * MDP(14) + t430 * MDP(15) + (-t293 * t361 - t332 * t355) * MDP(21) + (t292 * t361 - t334 * t355) * MDP(22) + (-t272 * t361 - t290 * t355 - t531) * MDP(24) + (MDP(14) * t353 - t305 * MDP(15) + (-t506 + t508) * MDP(23) + t443 * MDP(24) + t437 * t535) * t354 + (-t365 * MDP(15) - t527) * qJD(1) - (-t326 * MDP(14) - t285 * MDP(15) + (t512 + t513) * MDP(23) + (t268 * t415 - t269 * t418) * MDP(24) - t437 * t324 + t527 * qJD(5)) * t428; MDP(9) * t469 - MDP(10) * t470 + qJDD(3) * MDP(11) + (t419 * t431 + t363 + t519) * MDP(12) + (g(3) * t419 + (-t371 - t431) * t416) * MDP(13) + ((t305 - t314) * t356 - (-t315 + t304) * t353 + (t326 * t516 + t427 * t517) * pkin(3)) * MDP(14) + (t304 * t314 - t305 * t315 + (t516 * t285 + t517 * t284 + t519 + (-qJD(1) * t365 - t531) * t419) * pkin(3)) * MDP(15) + (t334 * t452 - t513) * MDP(16) + ((-t292 - t511) * t418 - t535 * t506 + t492) * MDP(17) + (t452 * t535 - t501 - t507) * MDP(18) + (-t477 * t535 + t490 + t510) * MDP(19) - t535 * t356 * MDP(20) + (-t282 * t356 + t387 * t293 - t313 * t535 - t314 * t332 + (t315 * t535 + t432) * t415 - t529 * t418) * MDP(21) + (t283 * t356 - t387 * t292 - t314 * t334 + t529 * t415 + t432 * t418 + t491 * t535) * MDP(22) + (-t292 * t358 - t293 * t359 - t488 * t332 - t487 * t334 - t531 * t395 - t539 * t415 + t538 * t418 - t520) * MDP(23) + (t269 * t359 - t268 * t358 + t272 * (t387 - t518) - g(3) * (-t403 - t441) + (t524 * t535 - t314) * t290 + t488 * t276 + t487 * t273 - t531 * (pkin(3) * t419 + t391 * t396 - t395 * t413)) * MDP(24) + (t419 * t416 * MDP(7) - t484 * MDP(8)) * t423; (-t353 ^ 2 - t356 ^ 2) * MDP(14) + (t304 * t356 + t305 * t353 + t440 - t449) * MDP(15) + (t490 - t510) * MDP(21) - MDP(22) * t507 + t492 * MDP(23) + (-t290 * t356 - t449) * MDP(24) + (MDP(23) * t334 * t535 - MDP(21) * t478 + t324 * MDP(22) + t538 * MDP(24)) * t415 + ((t292 - t511) * MDP(23) + t539 * MDP(24) - t535 ^ 2 * MDP(22)) * t418; t334 * t332 * MDP(16) + (-t331 + t526) * MDP(17) + (t332 * t535 - t292) * MDP(18) + (-t425 + (-qJD(5) + t535) * t334) * MDP(19) - t324 * MDP(20) + (t283 * t535 - t296 * t334 + (t451 + t520) * t415 + t445 + t530) * MDP(21) + (g(1) * t350 - g(2) * t352 + t282 * t535 + t296 * t332 + t418 * t520 - t433) * MDP(22) + (pkin(5) * t292 - t332 * t493) * MDP(23) + (t493 * t276 + (-t290 * t334 + t415 * t520 + t268 + t530) * pkin(5)) * MDP(24); (-t331 - t526) * MDP(23) + (t273 * t334 + t276 * t332 + t272 + t426) * MDP(24);];
tau  = t1;
