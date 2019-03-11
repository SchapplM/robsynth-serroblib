% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:39:08
% EndTime: 2019-03-09 03:39:16
% DurationCPUTime: 4.36s
% Computational Cost: add. (3799->378), mult. (9163->511), div. (0->0), fcn. (6778->10), ass. (0->179)
t394 = sin(pkin(11));
t396 = cos(pkin(11));
t403 = cos(qJ(3));
t451 = qJD(1) * t403;
t400 = sin(qJ(3));
t452 = qJD(1) * t400;
t364 = -t394 * t452 + t396 * t451;
t488 = qJD(5) + qJD(6);
t498 = t364 - t488;
t374 = t394 * t403 + t396 * t400;
t366 = t374 * qJD(1);
t399 = sin(qJ(5));
t402 = cos(qJ(5));
t445 = t402 * qJD(3);
t344 = t366 * t399 - t445;
t401 = cos(qJ(6));
t346 = qJD(3) * t399 + t366 * t402;
t398 = sin(qJ(6));
t475 = t346 * t398;
t288 = t401 * t344 + t475;
t360 = qJD(5) - t364;
t356 = qJD(6) + t360;
t497 = t288 * t356;
t418 = t344 * t398 - t401 * t346;
t496 = t356 * t418;
t385 = sin(pkin(10)) * pkin(1) + pkin(7);
t462 = qJ(4) + t385;
t466 = t398 * t402;
t376 = t399 * t401 + t466;
t455 = t498 * t376;
t450 = qJD(5) * t399;
t472 = t364 * t399;
t495 = t450 - t472;
t428 = t462 * qJD(1);
t347 = t403 * qJD(2) - t428 * t400;
t485 = qJD(3) * pkin(3);
t342 = t347 + t485;
t348 = t400 * qJD(2) + t428 * t403;
t467 = t396 * t348;
t285 = t394 * t342 + t467;
t282 = qJD(3) * pkin(8) + t285;
t387 = -cos(pkin(10)) * pkin(1) - pkin(2);
t416 = -pkin(3) * t403 + t387;
t408 = t416 * qJD(1);
t363 = qJD(4) + t408;
t303 = -pkin(4) * t364 - pkin(8) * t366 + t363;
t270 = t282 * t402 + t303 * t399;
t262 = -pkin(9) * t344 + t270;
t448 = qJD(6) * t398;
t260 = t262 * t448;
t339 = t394 * t348;
t284 = t342 * t396 - t339;
t281 = -qJD(3) * pkin(4) - t284;
t274 = pkin(5) * t344 + t281;
t494 = t274 * t288 + t260;
t444 = qJD(1) * qJD(3);
t435 = t403 * t444;
t436 = t400 * t444;
t358 = -t394 * t436 + t396 * t435;
t306 = qJD(5) * t445 + t402 * t358 - t366 * t450;
t365 = t374 * qJD(3);
t357 = qJD(1) * t365;
t443 = qJD(1) * qJD(4);
t330 = qJD(3) * t347 + t403 * t443;
t489 = -t348 * qJD(3) - t400 * t443;
t277 = t396 * t330 + t394 * t489;
t382 = pkin(3) * t436;
t309 = pkin(4) * t357 - pkin(8) * t358 + t382;
t305 = t402 * t309;
t406 = -t270 * qJD(5) - t277 * t399 + t305;
t250 = pkin(5) * t357 - pkin(9) * t306 + t406;
t307 = qJD(5) * t346 + t358 * t399;
t449 = qJD(5) * t402;
t411 = t402 * t277 - t282 * t450 + t303 * t449 + t399 * t309;
t253 = -pkin(9) * t307 + t411;
t432 = t401 * t250 - t398 * t253;
t493 = t274 * t418 + t432;
t492 = t357 * MDP(25) + (-t288 ^ 2 + t418 ^ 2) * MDP(22) - t288 * MDP(21) * t418;
t491 = MDP(5) * t400;
t324 = t376 * t374;
t490 = (t400 ^ 2 - t403 ^ 2) * MDP(6);
t375 = t398 * t399 - t401 * t402;
t456 = t498 * t375;
t487 = -t456 * t356 - t357 * t376;
t431 = t306 * t398 + t401 * t307;
t259 = -t418 * qJD(6) + t431;
t384 = pkin(3) * t394 + pkin(8);
t486 = pkin(9) + t384;
t269 = -t282 * t399 + t402 * t303;
t261 = -pkin(9) * t346 + t269;
t257 = pkin(5) * t360 + t261;
t484 = t257 * t401;
t483 = t262 * t401;
t482 = t288 * t366;
t481 = t418 * t366;
t480 = t306 * t399;
t479 = t344 * t360;
t478 = t344 * t366;
t477 = t346 * t360;
t476 = t346 * t366;
t373 = t394 * t400 - t396 * t403;
t368 = t373 * qJD(3);
t471 = t368 * t399;
t470 = t368 * t402;
t469 = t374 * t399;
t468 = t374 * t402;
t465 = t399 * t357;
t404 = qJD(3) ^ 2;
t464 = t400 * t404;
t369 = t462 * t400;
t370 = t462 * t403;
t328 = -t369 * t394 + t370 * t396;
t322 = t402 * t328;
t351 = t402 * t357;
t463 = t403 * t404;
t447 = qJD(6) * t401;
t439 = t401 * t306 - t398 * t307 - t344 * t447;
t258 = -t346 * t448 + t439;
t461 = t258 * t373 - t365 * t418;
t438 = t374 * t450;
t267 = -t368 * t466 - t398 * t438 - t448 * t469 + (t468 * t488 - t471) * t401;
t460 = -t267 * t356 - t324 * t357;
t459 = t306 * t373 + t346 * t365;
t293 = t347 * t396 - t339;
t319 = pkin(3) * t452 + pkin(4) * t366 - pkin(8) * t364;
t458 = t402 * t293 + t399 * t319;
t323 = pkin(4) * t373 - pkin(8) * t374 + t416;
t457 = t399 * t323 + t322;
t379 = qJD(1) * t387;
t442 = t400 * t485;
t441 = t374 * t465;
t440 = t374 * t351;
t386 = -pkin(3) * t396 - pkin(4);
t437 = t374 * t449;
t433 = qJD(5) * t486;
t276 = t330 * t394 - t396 * t489;
t292 = t347 * t394 + t467;
t430 = qJD(3) * t462;
t349 = qJD(4) * t403 - t400 * t430;
t350 = -qJD(4) * t400 - t403 * t430;
t301 = t349 * t394 - t396 * t350;
t327 = t396 * t369 + t370 * t394;
t429 = t360 * t402;
t427 = qJD(6) * t257 + t253;
t426 = pkin(5) * t495 - t292;
t425 = t356 * t455 - t375 * t357;
t312 = t402 * t319;
t372 = t486 * t402;
t424 = pkin(5) * t366 + qJD(6) * t372 - t293 * t399 + t312 + (-pkin(9) * t364 + t433) * t402;
t371 = t486 * t399;
t423 = -pkin(9) * t472 + qJD(6) * t371 + t399 * t433 + t458;
t252 = t257 * t398 + t483;
t422 = -t259 * t373 - t288 * t365;
t266 = -t324 * t488 + t375 * t368;
t325 = t375 * t374;
t421 = -t266 * t356 + t325 * t357;
t420 = t276 * t374 - t328 * t357;
t419 = -t307 * t373 - t344 * t365;
t415 = 0.2e1 * qJD(3) * t379;
t414 = -t360 * t495 + t351;
t413 = t437 - t471;
t412 = t438 + t470;
t302 = t349 * t396 + t350 * t394;
t320 = pkin(4) * t365 + pkin(8) * t368 + t442;
t410 = t402 * t302 + t399 * t320 + t323 * t449 - t328 * t450;
t409 = t360 * t281 - t384 * t357;
t377 = -pkin(5) * t402 + t386;
t331 = t357 * t373;
t318 = t402 * t323;
t313 = t402 * t320;
t297 = pkin(5) * t469 + t327;
t273 = t413 * pkin(5) + t301;
t272 = -pkin(9) * t469 + t457;
t271 = pkin(5) * t373 - pkin(9) * t468 - t328 * t399 + t318;
t265 = pkin(5) * t307 + t276;
t255 = -t413 * pkin(9) + t410;
t254 = pkin(9) * t470 + pkin(5) * t365 - t302 * t399 + t313 + (-t322 + (pkin(9) * t374 - t323) * t399) * qJD(5);
t251 = -t262 * t398 + t484;
t1 = [0.2e1 * t435 * t491 - 0.2e1 * t444 * t490 + MDP(7) * t463 - MDP(8) * t464 + (-t385 * t463 + t400 * t415) * MDP(10) + (t385 * t464 + t403 * t415) * MDP(11) + (-t277 * t373 + t284 * t368 - t285 * t365 + t301 * t366 + t302 * t364 + t327 * t358 + t420) * MDP(12) + (t276 * t327 + t277 * t328 - t284 * t301 + t285 * t302 + (t363 + t408) * t442) * MDP(13) + (t306 * t468 - t412 * t346) * MDP(14) + (-(-t344 * t402 - t346 * t399) * t368 + (-t480 - t307 * t402 + (t344 * t399 - t346 * t402) * qJD(5)) * t374) * MDP(15) + (-t412 * t360 + t440 + t459) * MDP(16) + (-t413 * t360 + t419 - t441) * MDP(17) + (t360 * t365 + t331) * MDP(18) + ((-t328 * t449 + t313) * t360 + t318 * t357 + (-t282 * t449 + t305) * t373 + t269 * t365 + t301 * t344 + t327 * t307 + t281 * t437 + ((-qJD(5) * t323 - t302) * t360 + (-qJD(5) * t303 - t277) * t373 - t281 * t368 + t420) * t399) * MDP(19) + (-t270 * t365 + t276 * t468 - t412 * t281 + t301 * t346 + t327 * t306 - t457 * t357 - t410 * t360 - t411 * t373) * MDP(20) + (-t258 * t325 - t266 * t418) * MDP(21) + (-t258 * t324 + t259 * t325 - t266 * t288 + t267 * t418) * MDP(22) + (-t421 + t461) * MDP(23) + (t422 + t460) * MDP(24) + (t356 * t365 + t331) * MDP(25) + ((t254 * t401 - t255 * t398) * t356 + (t271 * t401 - t272 * t398) * t357 + t432 * t373 + t251 * t365 + t273 * t288 + t297 * t259 + t265 * t324 + t274 * t267 + ((-t271 * t398 - t272 * t401) * t356 - t252 * t373) * qJD(6)) * MDP(26) + (-t252 * t365 + t297 * t258 + t260 * t373 - t265 * t325 + t274 * t266 - t273 * t418 + (-(-qJD(6) * t272 + t254) * t356 - t271 * t357 - t250 * t373) * t398 + (-(qJD(6) * t271 + t255) * t356 - t272 * t357 - t427 * t373) * t401) * MDP(27); (-t357 * t374 + t358 * t373 - t364 * t368 + t365 * t366) * MDP(12) + (t276 * t373 + t277 * t374 - t284 * t365 - t285 * t368) * MDP(13) + (-t419 - t441) * MDP(19) + (-t440 + t459) * MDP(20) + (-t422 + t460) * MDP(26) + (t421 + t461) * MDP(27) + (-MDP(10) * t400 - MDP(11) * t403) * t404 + (-t413 * MDP(19) + t412 * MDP(20)) * t360; ((t284 - t293) * t364 + (-t357 * t394 - t358 * t396) * pkin(3)) * MDP(12) + (t284 * t292 - t285 * t293 + (-t276 * t396 + t277 * t394 - t363 * t452) * pkin(3)) * MDP(13) + (t346 * t429 + t480) * MDP(14) + ((t306 - t479) * t402 + (-t307 - t477) * t399) * MDP(15) + (t360 * t429 + t465 - t476) * MDP(16) + (t414 + t478) * MDP(17) + (-t276 * t402 - t292 * t344 + t386 * t307 + (-t384 * t449 - t312) * t360 + (t293 * t360 + t409) * t399) * MDP(19) + (t276 * t399 - t292 * t346 + t386 * t306 + (t384 * t450 + t458) * t360 + t409 * t402) * MDP(20) + (t258 * t376 - t418 * t456) * MDP(21) + (-t258 * t375 - t259 * t376 - t456 * t288 - t418 * t455) * MDP(22) + (t481 - t487) * MDP(23) + (t425 + t482) * MDP(24) + ((-t371 * t401 - t372 * t398) * t357 + t377 * t259 + t265 * t375 + (t423 * t398 - t424 * t401) * t356 + t426 * t288 - t455 * t274) * MDP(26) + (-(-t371 * t398 + t372 * t401) * t357 + t377 * t258 + t265 * t376 + (t424 * t398 + t423 * t401) * t356 - t426 * t418 + t456 * t274) * MDP(27) + (-t403 * t491 + t490) * qJD(1) ^ 2 + (-MDP(10) * t452 - MDP(11) * t451) * t379 - ((-t285 + t292) * MDP(12) + t360 * MDP(18) + t269 * MDP(19) - t270 * MDP(20) + t356 * MDP(25) + t251 * MDP(26) - t252 * MDP(27)) * t366; (-t364 ^ 2 - t366 ^ 2) * MDP(12) + (t284 * t366 - t285 * t364 + t382) * MDP(13) + (t414 - t478) * MDP(19) + (-t360 ^ 2 * t402 - t465 - t476) * MDP(20) + (t425 - t482) * MDP(26) + (t481 + t487) * MDP(27); t346 * t344 * MDP(14) + (-t344 ^ 2 + t346 ^ 2) * MDP(15) + (t306 + t479) * MDP(16) + (-t307 + t477) * MDP(17) + t357 * MDP(18) + (t270 * t360 - t281 * t346 + t406) * MDP(19) + (t269 * t360 + t281 * t344 - t411) * MDP(20) + (t258 + t497) * MDP(23) + (-t259 - t496) * MDP(24) + (-(-t261 * t398 - t483) * t356 - t252 * qJD(6) + (-t288 * t346 - t356 * t448 + t357 * t401) * pkin(5) + t493) * MDP(26) + ((-t262 * t356 - t250) * t398 + (t261 * t356 - t427) * t401 + (t346 * t418 - t356 * t447 - t357 * t398) * pkin(5) + t494) * MDP(27) + t492; (t439 + t497) * MDP(23) + (-t431 - t496) * MDP(24) + (t252 * t356 + t493) * MDP(26) + (-t398 * t250 + t251 * t356 - t401 * t253 + t494) * MDP(27) + (-MDP(23) * t475 + t418 * MDP(24) - t252 * MDP(26) - MDP(27) * t484) * qJD(6) + t492;];
tauc  = t1;
