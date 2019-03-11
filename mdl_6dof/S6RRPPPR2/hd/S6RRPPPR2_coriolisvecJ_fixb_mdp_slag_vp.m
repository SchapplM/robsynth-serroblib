% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:16
% EndTime: 2019-03-09 08:12:25
% DurationCPUTime: 5.86s
% Computational Cost: add. (3612->407), mult. (9218->538), div. (0->0), fcn. (6674->8), ass. (0->180)
t411 = sin(pkin(9));
t413 = cos(pkin(9));
t415 = sin(qJ(2));
t417 = cos(qJ(2));
t386 = t411 * t417 + t413 * t415;
t375 = t386 * qJD(1);
t365 = qJD(6) + t375;
t453 = qJD(1) * qJD(2);
t448 = t415 * t453;
t395 = t411 * t448;
t447 = t417 * t453;
t362 = t413 * t447 - t395;
t410 = sin(pkin(10));
t412 = cos(pkin(10));
t414 = sin(qJ(6));
t416 = cos(qJ(6));
t385 = t410 * t416 + t412 * t414;
t465 = t365 * t385;
t497 = -t410 * t414 + t412 * t416;
t502 = t497 * t362 - t365 * t465;
t470 = t413 * t417;
t384 = t411 * t415 - t470;
t328 = t497 * t384;
t460 = qJD(6) * t416;
t461 = qJD(6) * t414;
t464 = t497 * t375 - t410 * t461 + t412 * t460;
t440 = -t385 * t362 - t365 * t464;
t458 = t415 * qJD(1);
t372 = -qJD(1) * t470 + t411 * t458;
t344 = qJD(2) * t410 - t412 * t372;
t346 = qJD(2) * t412 + t372 * t410;
t432 = t344 * t414 - t346 * t416;
t501 = t365 * t432;
t500 = -0.2e1 * t453;
t499 = MDP(4) * t415;
t498 = MDP(5) * (t415 ^ 2 - t417 ^ 2);
t485 = -qJ(3) - pkin(7);
t394 = t485 * t417;
t390 = qJD(1) * t394;
t380 = t411 * t390;
t393 = t485 * t415;
t389 = qJD(1) * t393;
t339 = t389 * t413 + t380;
t495 = qJD(4) - t339;
t383 = qJD(2) * pkin(2) + t389;
t471 = t413 * t390;
t333 = t411 * t383 - t471;
t330 = -qJD(2) * qJ(4) - t333;
t490 = pkin(4) * t372;
t306 = qJD(5) - t330 - t490;
t403 = -pkin(2) * t413 - pkin(3);
t399 = -qJ(5) + t403;
t374 = t386 * qJD(2);
t361 = qJD(1) * t374;
t401 = pkin(2) * t411 + qJ(4);
t479 = t361 * t401;
t494 = -t362 * t399 + t375 * (qJD(5) - t306) + t479;
t491 = t375 ^ 2;
t493 = -t362 * t410 - t412 * t491;
t489 = pkin(4) * t375;
t488 = pkin(8) * t412;
t487 = pkin(3) + qJ(5);
t486 = -pkin(8) + t399;
t480 = t346 * t414;
t303 = t416 * t344 + t480;
t484 = t303 * t365;
t483 = t303 * t372;
t482 = t432 * t372;
t446 = qJD(2) * t485;
t366 = qJD(3) * t417 + t415 * t446;
t356 = t366 * qJD(1);
t367 = -qJD(3) * t415 + t417 * t446;
t357 = t367 * qJD(1);
t313 = t356 * t411 - t413 * t357;
t342 = -t413 * t393 - t394 * t411;
t481 = t313 * t342;
t478 = t361 * t412;
t476 = t362 * t412;
t418 = qJD(2) ^ 2;
t469 = t415 * t418;
t468 = t417 * t418;
t419 = qJD(1) ^ 2;
t467 = t417 * t419;
t400 = pkin(2) * t448;
t444 = -qJ(4) * t362 + t400;
t427 = -qJD(4) * t375 + t444;
t275 = qJD(5) * t372 + t361 * t487 + t427;
t287 = pkin(4) * t362 - qJD(2) * qJD(5) + t313;
t258 = t412 * t275 + t410 * t287;
t462 = qJD(2) * t415;
t377 = qJD(2) * t470 - t411 * t462;
t405 = pkin(2) * t462;
t426 = -qJ(4) * t377 - qJD(4) * t386 + t405;
t280 = qJD(5) * t384 + t374 * t487 + t426;
t323 = t366 * t411 - t413 * t367;
t301 = pkin(4) * t377 + t323;
t263 = t412 * t280 + t410 * t301;
t450 = -pkin(2) * t417 - pkin(1);
t439 = t450 * qJD(1);
t392 = qJD(3) + t439;
t421 = -qJ(4) * t375 + t392;
t294 = t372 * t487 + t421;
t332 = t383 * t413 + t380;
t438 = qJD(4) - t332;
t298 = -qJD(2) * t487 + t438 + t489;
t268 = t412 * t294 + t410 * t298;
t443 = pkin(2) * t458 + qJ(4) * t372;
t299 = t375 * t487 + t443;
t338 = t389 * t411 - t471;
t315 = t338 - t490;
t271 = t412 * t299 + t410 * t315;
t429 = -qJ(4) * t386 + t450;
t312 = t384 * t487 + t429;
t325 = pkin(4) * t386 + t342;
t277 = t412 * t312 + t410 * t325;
t314 = t413 * t356 + t411 * t357;
t459 = t375 * MDP(14);
t449 = -pkin(5) * t412 - pkin(4);
t456 = -t375 * t449 + t495;
t455 = t489 + t495;
t452 = -t344 * t460 + t361 * t385;
t406 = qJD(2) * qJD(4);
t451 = t406 + t314;
t445 = pkin(1) * t500;
t286 = t412 * t287;
t254 = pkin(5) * t362 + t286 + (-pkin(8) * t361 - t275) * t410;
t255 = pkin(8) * t478 + t258;
t442 = t416 * t254 - t255 * t414;
t267 = -t294 * t410 + t412 * t298;
t441 = t497 * t361;
t437 = t254 * t414 + t255 * t416;
t257 = -t275 * t410 + t286;
t436 = t257 * t412 + t258 * t410;
t260 = pkin(5) * t375 - pkin(8) * t346 + t267;
t261 = -pkin(8) * t344 + t268;
t251 = t260 * t416 - t261 * t414;
t252 = t260 * t414 + t261 * t416;
t321 = t412 * t325;
t266 = pkin(5) * t386 + t321 + (-pkin(8) * t384 - t312) * t410;
t269 = t384 * t488 + t277;
t435 = t266 * t416 - t269 * t414;
t434 = t266 * t414 + t269 * t416;
t433 = -t267 * t412 - t268 * t410;
t324 = t366 * t413 + t367 * t411;
t343 = t393 * t411 - t394 * t413;
t288 = -pkin(4) * t361 + t451;
t317 = pkin(3) * t372 + t421;
t428 = t317 * t375 + t313;
t311 = t412 * t315;
t370 = t486 * t410;
t425 = qJD(5) * t412 + qJD(6) * t370 - pkin(5) * t372 + t311 + (-pkin(8) * t375 - t299) * t410;
t371 = t486 * t412;
t424 = qJD(5) * t410 - qJD(6) * t371 + t375 * t488 + t271;
t329 = t385 * t384;
t272 = -t346 * t461 + t452;
t326 = -pkin(4) * t384 + t343;
t422 = t288 * t384 + t306 * t374 + t326 * t361;
t273 = -qJD(6) * t432 - t441;
t420 = t313 * t386 + t323 * t375 - t324 * t372 + t342 * t362 - t343 * t361;
t391 = pkin(5) * t410 + t401;
t364 = qJD(2) * t372;
t331 = pkin(3) * t384 + t429;
t327 = -qJD(2) * pkin(3) + t438;
t322 = pkin(3) * t375 + t443;
t307 = pkin(3) * t374 + t426;
t302 = -pkin(4) * t374 + t324;
t300 = t384 * t449 + t343;
t297 = t412 * t301;
t291 = pkin(3) * t361 + t427;
t284 = t374 * t449 + t324;
t283 = pkin(5) * t344 + t306;
t282 = qJD(6) * t329 - t374 * t497;
t281 = qJD(6) * t328 + t374 * t385;
t278 = t361 * t449 + t451;
t276 = -t312 * t410 + t321;
t270 = -t299 * t410 + t311;
t262 = -t280 * t410 + t297;
t259 = t374 * t488 + t263;
t256 = pkin(5) * t377 + t297 + (-pkin(8) * t374 - t280) * t410;
t1 = [0.2e1 * t447 * t499 + t498 * t500 + MDP(6) * t468 - MDP(7) * t469 + (-pkin(7) * t468 + t415 * t445) * MDP(9) + (pkin(7) * t469 + t417 * t445) * MDP(10) + (-t314 * t384 - t332 * t377 - t333 * t374 + t420) * MDP(11) + (t481 + t314 * t343 - t332 * t323 + t333 * t324 + (t392 + t439) * t405) * MDP(12) + (t327 * t377 + t330 * t374 - t384 * t451 + t420) * MDP(13) + (qJD(2) * t323 - t291 * t384 - t307 * t372 - t317 * t374 - t331 * t361) * MDP(14) + (qJD(2) * t324 - t291 * t386 - t307 * t375 - t317 * t377 - t331 * t362) * MDP(15) + (t291 * t331 + t307 * t317 + t323 * t327 - t324 * t330 + t343 * t451 + t481) * MDP(16) + (t257 * t386 + t262 * t375 + t267 * t377 + t276 * t362 + t302 * t344 - t412 * t422) * MDP(17) + (-t258 * t386 - t263 * t375 - t268 * t377 - t277 * t362 + t302 * t346 + t410 * t422) * MDP(18) + (-t262 * t346 - t263 * t344 + (t258 * t384 + t268 * t374 + t277 * t361) * t412 + (-t257 * t384 - t267 * t374 - t276 * t361) * t410) * MDP(19) + (t257 * t276 + t258 * t277 + t262 * t267 + t263 * t268 + t288 * t326 + t302 * t306) * MDP(20) + (t272 * t329 - t281 * t432) * MDP(21) + (t272 * t328 - t273 * t329 - t281 * t303 + t282 * t432) * MDP(22) + (t272 * t386 + t281 * t365 + t329 * t362 - t377 * t432) * MDP(23) + (-t273 * t386 - t282 * t365 - t303 * t377 + t328 * t362) * MDP(24) + (t362 * t386 + t365 * t377) * MDP(25) + ((t256 * t416 - t259 * t414) * t365 + t435 * t362 + t442 * t386 + t251 * t377 + t284 * t303 + t300 * t273 - t278 * t328 + t283 * t282 + (-t252 * t386 - t365 * t434) * qJD(6)) * MDP(26) + (-(t256 * t414 + t259 * t416) * t365 - t434 * t362 - t437 * t386 - t252 * t377 - t284 * t432 + t300 * t272 + t278 * t329 + t283 * t281 + (-t251 * t386 - t365 * t435) * qJD(6)) * MDP(27); -t467 * t499 + t419 * t498 + ((t333 - t338) * t375 + (-t332 + t339) * t372 + (-t361 * t411 - t362 * t413) * pkin(2)) * MDP(11) + (t332 * t338 - t333 * t339 + (-t313 * t413 + t314 * t411 - t392 * t458) * pkin(2)) * MDP(12) + (-t479 + t362 * t403 + (-t330 - t338) * t375 + (t327 - t495) * t372) * MDP(13) + (-qJD(2) * t338 + t322 * t372 + t428) * MDP(14) + (-qJD(2) * t339 - t317 * t372 + t322 * t375 + t314 + 0.2e1 * t406) * MDP(15) + (t313 * t403 - t317 * t322 - t327 * t338 - t330 * t495 + t401 * t451) * MDP(16) + (t267 * t372 - t270 * t375 + t288 * t410 + t455 * t344 - t412 * t494) * MDP(17) + (-t268 * t372 + t271 * t375 + t288 * t412 + t455 * t346 + t410 * t494) * MDP(18) + (t270 * t346 + t271 * t344 + (qJD(5) * t346 - t268 * t375 - t257) * t412 + (qJD(5) * t344 + t267 * t375 - t258) * t410) * MDP(19) + (qJD(5) * t433 - t267 * t270 - t268 * t271 + t288 * t401 + t306 * t455 + t399 * t436) * MDP(20) + (t272 * t497 + t432 * t465) * MDP(21) + (-t272 * t385 - t273 * t497 + t303 * t465 + t432 * t464) * MDP(22) + (-t482 + t502) * MDP(23) + (t440 - t483) * MDP(24) + t365 * t372 * MDP(25) + ((-t370 * t414 + t371 * t416) * t362 + t391 * t273 + t278 * t385 + t251 * t372 + (t414 * t424 - t416 * t425) * t365 + t456 * t303 + t464 * t283) * MDP(26) + (-(t370 * t416 + t371 * t414) * t362 + t391 * t272 + t278 * t497 - t252 * t372 + (t414 * t425 + t416 * t424) * t365 - t456 * t432 - t465 * t283) * MDP(27) + (MDP(9) * t415 * t419 + MDP(10) * t467) * pkin(1); (t333 * t372 + t400) * MDP(12) + (t364 + t395) * MDP(15) + (-t330 * t372 + t444) * MDP(16) + (t344 * t372 + t493) * MDP(17) + (t346 * t372 - t476) * MDP(18) + (-t257 * t410 + t258 * t412 + t306 * t372) * MDP(20) + (t440 + t483) * MDP(26) + (-t482 - t502) * MDP(27) + (pkin(3) * MDP(16) + (t410 ^ 2 + t412 ^ 2) * MDP(19)) * t361 + (t332 * MDP(12) + (-qJD(4) - t327) * MDP(16) + (t344 * t410 + t346 * t412) * MDP(19) + t433 * MDP(20) + t410 * MDP(18) * t375) * t375 + (-t459 + (-MDP(14) * t386 - MDP(15) * t470) * qJD(1)) * qJD(2) + (MDP(11) + MDP(13)) * (-t372 ^ 2 - t491); (t362 + t364) * MDP(13) - t372 * t459 + (-t491 - t418) * MDP(15) + (qJD(2) * t330 + t428) * MDP(16) + (-qJD(2) * t344 - t410 * t491 + t476) * MDP(17) + (-qJD(2) * t346 + t493) * MDP(18) + (-t344 * t412 + t346 * t410) * t375 * MDP(19) + (-qJD(2) * t306 + (-t267 * t410 + t268 * t412) * t375 + t436) * MDP(20) + (-qJD(2) * t303 + t502) * MDP(26) + (qJD(2) * t432 + t440) * MDP(27); (t346 * t375 - t478) * MDP(17) + (-t344 * t375 + t361 * t410) * MDP(18) + (-t344 ^ 2 - t346 ^ 2) * MDP(19) + (t267 * t346 + t268 * t344 + t288) * MDP(20) + (t273 - t501) * MDP(26) + (t272 - t484) * MDP(27); -t432 * t303 * MDP(21) + (-t303 ^ 2 + t432 ^ 2) * MDP(22) + (t452 + t484) * MDP(23) + (t441 - t501) * MDP(24) + t362 * MDP(25) + (t252 * t365 + t283 * t432 + t442) * MDP(26) + (t251 * t365 + t283 * t303 - t437) * MDP(27) + (-MDP(23) * t480 + MDP(24) * t432 - MDP(26) * t252 - MDP(27) * t251) * qJD(6);];
tauc  = t1;
