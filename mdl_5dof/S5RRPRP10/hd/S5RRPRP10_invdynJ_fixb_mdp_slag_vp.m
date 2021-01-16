% Calculate vector of inverse dynamics joint torques for
% S5RRPRP10
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:55
% EndTime: 2021-01-15 21:03:07
% DurationCPUTime: 3.98s
% Computational Cost: add. (2214->419), mult. (4531->508), div. (0->0), fcn. (2592->6), ass. (0->196)
t485 = pkin(3) + pkin(6);
t365 = sin(qJ(4));
t369 = cos(qJ(2));
t445 = qJD(1) * t369;
t416 = t365 * t445;
t368 = cos(qJ(4));
t440 = qJD(2) * t368;
t315 = -t416 + t440;
t371 = -pkin(2) - pkin(7);
t366 = sin(qJ(2));
t476 = qJ(3) * t366;
t309 = t371 * t369 - pkin(1) - t476;
t284 = t309 * qJD(1);
t446 = qJD(1) * t366;
t350 = pkin(6) * t446;
t490 = qJD(3) + t350;
t431 = pkin(3) * t446 + t490;
t288 = t371 * qJD(2) + t431;
t265 = t284 * t368 + t288 * t365;
t430 = qJD(1) * qJD(2);
t415 = t366 * t430;
t338 = pkin(2) * t415;
t475 = qJ(3) * t369;
t396 = pkin(7) * t366 - t475;
t438 = qJD(3) * t366;
t381 = t396 * qJD(2) - t438;
t263 = qJD(1) * t381 + qJDD(1) * t309 + t338;
t414 = t369 * t430;
t427 = qJDD(1) * t366;
t386 = t414 + t427;
t337 = pkin(6) * t414;
t347 = pkin(6) * t427;
t413 = qJDD(3) + t337 + t347;
t274 = t386 * pkin(3) + t371 * qJDD(2) + t413;
t409 = -t263 * t365 + t368 * t274;
t378 = -t265 * qJD(4) + t409;
t426 = qJDD(1) * t369;
t442 = qJD(2) * t365;
t313 = t368 * t445 + t442;
t437 = qJD(4) * t313;
t450 = t368 * qJDD(2) + t365 * t415;
t271 = t365 * t426 + t437 - t450;
t474 = qJ(5) * t271;
t312 = qJDD(4) + t386;
t483 = pkin(4) * t312;
t251 = -qJD(5) * t315 + t378 + t474 + t483;
t258 = -qJ(5) * t313 + t265;
t340 = qJD(4) + t446;
t471 = t258 * t340;
t497 = t251 + t471;
t435 = qJD(4) * t368;
t436 = qJD(4) * t365;
t407 = t368 * t263 + t365 * t274 - t284 * t436 + t288 * t435;
t395 = -qJD(4) * t416 + qJDD(2) * t365 - t368 * t415;
t429 = qJD(2) * qJD(4);
t272 = (t426 + t429) * t368 + t395;
t473 = qJ(5) * t272;
t252 = -qJD(5) * t313 + t407 - t473;
t264 = -t284 * t365 + t368 * t288;
t257 = -qJ(5) * t315 + t264;
t256 = pkin(4) * t340 + t257;
t496 = -t340 * t256 + t252;
t468 = t315 * t340;
t495 = -t272 + t468;
t367 = sin(qJ(1));
t370 = cos(qJ(1));
t400 = g(1) * t370 + g(2) * t367;
t493 = t400 * t369;
t492 = qJ(5) - t371;
t330 = t485 * t366;
t451 = t368 * t309 + t365 * t330;
t297 = t368 * t312;
t491 = -t340 * t436 + t297;
t489 = pkin(4) * t272 + qJDD(5);
t458 = t368 * t370;
t301 = -t365 * t367 + t366 * t458;
t460 = t367 * t368;
t303 = t365 * t370 + t366 * t460;
t459 = t368 * t369;
t488 = -g(1) * t301 - g(2) * t303 + g(3) * t459;
t351 = pkin(6) * t445;
t322 = pkin(3) * t445 + t351;
t362 = qJD(2) * qJ(3);
t299 = t362 + t322;
t487 = t340 * t299 + t371 * t312;
t486 = t315 ^ 2;
t481 = g(1) * t367;
t479 = g(2) * t370;
t478 = g(3) * t366;
t477 = pkin(6) * qJDD(2);
t472 = qJDD(2) * pkin(2);
t470 = t271 * t368;
t469 = t313 * t340;
t467 = t365 * t312;
t466 = t365 * t366;
t465 = t365 * t369;
t464 = t366 * t367;
t463 = t366 * t368;
t462 = t366 * t370;
t373 = qJD(1) ^ 2;
t461 = t366 * t373;
t455 = -t257 + t256;
t354 = pkin(2) * t446;
t293 = t396 * qJD(1) + t354;
t454 = t368 * t293 + t365 * t322;
t307 = t368 * t322;
t453 = -qJD(5) * t368 + t492 * t436 + t293 * t365 - t307 - (pkin(4) * t369 - qJ(5) * t466) * qJD(1);
t325 = t492 * t368;
t452 = -qJ(5) * t368 * t446 - qJD(4) * t325 - qJD(5) * t365 - t454;
t422 = -pkin(4) * t368 - pkin(3);
t449 = pkin(4) * t435 - t422 * t446 + t490;
t331 = t485 * t369;
t363 = t366 ^ 2;
t364 = t369 ^ 2;
t448 = t363 - t364;
t398 = pkin(2) * t369 + t476;
t327 = pkin(1) + t398;
t447 = qJD(1) * t327;
t444 = qJD(2) * t313;
t443 = qJD(2) * t315;
t441 = qJD(2) * t366;
t439 = qJD(2) * t369;
t434 = qJD(4) * t369;
t433 = qJD(4) * t371;
t432 = qJD(5) * t369;
t428 = qJDD(1) * t327;
t425 = t369 * t461;
t348 = pkin(6) * t426;
t360 = qJDD(2) * qJ(3);
t361 = qJD(2) * qJD(3);
t424 = t348 + t360 + t361;
t423 = -g(1) * t462 - g(2) * t464 + g(3) * t369;
t421 = qJD(2) * t485;
t412 = -pkin(4) * t313 - qJD(5);
t276 = -t412 + t299;
t420 = t276 * t436;
t419 = t276 * t435;
t417 = t365 * t434;
t411 = -qJD(2) * pkin(2) + qJD(3);
t410 = qJ(5) * t369 - t309;
t346 = pkin(4) * t365 + qJ(3);
t408 = t346 * t366 + t369 * t492;
t345 = pkin(6) - t422;
t406 = -t347 - t423;
t405 = pkin(3) * t426 + t424;
t404 = qJD(1) * t421;
t372 = qJD(2) ^ 2;
t403 = pkin(6) * t372 + t479;
t402 = g(1) * t303 - g(2) * t301;
t302 = t365 * t462 + t460;
t304 = -t365 * t464 + t458;
t401 = -g(1) * t304 - g(2) * t302;
t399 = -t479 + t481;
t397 = -pkin(2) * t366 + t475;
t326 = t350 + t411;
t328 = -t351 - t362;
t394 = t326 * t369 + t328 * t366;
t393 = t340 * t365;
t390 = -0.2e1 * pkin(1) * t430 - t477;
t389 = -t446 * t447 + qJDD(3) - t406;
t388 = -t340 * t435 - t467;
t387 = -qJ(3) * t439 - t438;
t353 = pkin(2) * t441;
t281 = t353 + t381;
t323 = t485 * t439;
t385 = t368 * t281 - t309 * t436 + t365 * t323 + t330 * t435;
t384 = 0.2e1 * qJDD(1) * pkin(1) - t403;
t383 = 0.2e1 * t447 * qJD(2) + t477;
t382 = -t478 - t493;
t275 = -t366 * t404 + t405;
t380 = t275 + t382;
t255 = t275 + t489;
t379 = t255 + t382;
t273 = t387 * qJD(1) + t338 - t428;
t296 = t353 + t387;
t377 = qJD(1) * t296 + t273 + t403 - t428;
t376 = g(1) * t302 - g(2) * t304 - g(3) * t465 - t407;
t375 = t378 + t488;
t285 = pkin(6) * t415 - t424;
t291 = t413 - t472;
t374 = t394 * qJD(2) - t285 * t369 + t291 * t366 - t400;
t343 = t369 * t481;
t324 = t492 * t365;
t321 = t366 * t421;
t319 = -qJ(3) * t445 + t354;
t318 = t368 * t330;
t311 = t313 ^ 2;
t308 = t368 * t323;
t298 = pkin(4) * t459 + t331;
t292 = pkin(1) + t408;
t277 = -pkin(4) * t417 - t345 * t441;
t268 = -qJ(5) * t459 + t451;
t267 = pkin(4) * t366 + t410 * t365 + t318;
t254 = -t368 * t432 + (t366 * t440 + t417) * qJ(5) + t385;
t253 = pkin(4) * t439 + t308 + t410 * t435 + (-qJ(5) * t441 - qJD(4) * t330 - t281 + t432) * t365;
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t363 + 0.2e1 * t366 * t414) * MDP(4) + 0.2e1 * (t366 * t426 - t448 * t430) * MDP(5) + (qJDD(2) * t366 + t369 * t372) * MDP(6) + (qJDD(2) * t369 - t366 * t372) * MDP(7) + (t390 * t366 + t384 * t369 + t343) * MDP(9) + (t390 * t369 + (-t384 - t481) * t366) * MDP(10) + ((t363 + t364) * qJDD(1) * pkin(6) + t374) * MDP(11) + (t366 * t383 + t369 * t377 - t343) * MDP(12) + (t383 * t369 + (-t377 + t481) * t366) * MDP(13) + (-t447 * t296 + (-t273 + t399) * t327 + t374 * pkin(6)) * MDP(14) + (t271 * t465 + (t365 * t441 - t368 * t434) * t315) * MDP(15) + ((-t313 * t365 + t315 * t368) * t441 + (t470 + t272 * t365 + (t313 * t368 + t315 * t365) * qJD(4)) * t369) * MDP(16) + ((t340 * t442 - t271) * t366 + (t388 + t443) * t369) * MDP(17) + ((t340 * t440 - t272) * t366 + (-t444 - t491) * t369) * MDP(18) + (t312 * t366 + t340 * t439) * MDP(19) + ((-t281 * t365 + t308) * t340 + (-t309 * t365 + t318) * t312 + t409 * t366 - t321 * t313 + t331 * t272 + t275 * t459 + (t264 * t369 - t299 * t463) * qJD(2) + (-t265 * t366 - t299 * t465 - t451 * t340) * qJD(4) + t401) * MDP(20) + (-t385 * t340 - t451 * t312 - t321 * t315 - t331 * t271 + (t299 * t442 - t407) * t366 + (-qJD(2) * t265 - t275 * t365 - t299 * t435) * t369 + t402) * MDP(21) + (t253 * t340 + t267 * t312 + t272 * t298 + t277 * t313 + (-t276 * t440 + t251) * t366 + (qJD(2) * t256 + t255 * t368 - t420) * t369 + t401) * MDP(22) + (-t254 * t340 - t268 * t312 - t271 * t298 + t277 * t315 + (t276 * t442 - t252) * t366 + (-qJD(2) * t258 - t255 * t365 - t419) * t369 + t402) * MDP(23) + (-t253 * t315 - t254 * t313 + t267 * t271 - t268 * t272 + t343 + (-t256 * t365 + t258 * t368) * t441 + (-t479 + t251 * t365 - t252 * t368 + (t256 * t368 + t258 * t365) * qJD(4)) * t369) * MDP(24) + (t252 * t268 + t258 * t254 + t251 * t267 + t256 * t253 + t255 * t298 + t276 * t277 - g(1) * (-t292 * t367 + t345 * t370) - g(2) * (t292 * t370 + t345 * t367)) * MDP(25) + t399 * MDP(2) + t400 * MDP(3); -MDP(4) * t425 + t448 * MDP(5) * t373 + MDP(6) * t427 + MDP(7) * t426 + qJDD(2) * MDP(8) + (pkin(1) * t461 + t406) * MDP(9) + (t478 - t348 + (pkin(1) * t373 + t400) * t369) * MDP(10) + (t397 * qJDD(1) + ((-t328 - t362) * t366 + (-t326 + t411) * t369) * qJD(1)) * MDP(11) + (-t319 * t445 + t389 - 0.2e1 * t472) * MDP(12) + (t348 + 0.2e1 * t360 + 0.2e1 * t361 + (qJD(1) * t319 - g(3)) * t366 + (-qJD(1) * t447 - t400) * t369) * MDP(13) + (-t394 * qJD(1) * pkin(6) - t291 * pkin(2) - g(3) * t398 - t285 * qJ(3) - t328 * qJD(3) + t319 * t447 - t400 * t397) * MDP(14) + (-t315 * t393 - t470) * MDP(15) + ((-t272 - t468) * t368 + (t271 + t469) * t365) * MDP(16) + ((-t315 * t369 - t340 * t466) * qJD(1) + t491) * MDP(17) + ((t313 * t369 - t340 * t463) * qJD(1) + t388) * MDP(18) - t340 * MDP(19) * t445 + (-t264 * t445 + qJ(3) * t272 - t307 * t340 + t431 * t313 + t487 * t368 + ((t293 - t433) * t340 + t380) * t365) * MDP(20) + (-qJ(3) * t271 + t454 * t340 + t265 * t445 + t431 * t315 - t487 * t365 + (-t340 * t433 + t380) * t368) * MDP(21) + (t419 + t272 * t346 - t312 * t325 + t453 * t340 + t449 * t313 + (-t256 * t369 + t276 * t463) * qJD(1) + t379 * t365) * MDP(22) + (-t420 - t271 * t346 + t312 * t324 - t452 * t340 + t449 * t315 + (t258 * t369 - t276 * t466) * qJD(1) + t379 * t368) * MDP(23) + (-t271 * t325 + t272 * t324 - t452 * t313 - t453 * t315 - t496 * t365 - t497 * t368 - t423) * MDP(24) + (-t252 * t324 - t251 * t325 + t255 * t346 - g(3) * t408 - t400 * (t346 * t369 - t366 * t492) + t449 * t276 + t452 * t258 + t453 * t256) * MDP(25); MDP(11) * t427 + (qJDD(2) + t425) * MDP(12) + (-t363 * t373 - t372) * MDP(13) + (qJD(2) * t328 + t337 + t389 - t472) * MDP(14) + (-qJD(2) * t276 + t423) * MDP(25) + (MDP(21) + MDP(23)) * (-t340 ^ 2 * t368 - t443 - t467) + (MDP(20) + MDP(22)) * (-t340 * t393 + t297 - t444) + ((-t313 * t446 + t271 - t437) * MDP(24) + t497 * MDP(25)) * t368 + (t495 * MDP(24) + t496 * MDP(25)) * t365; t315 * t313 * MDP(15) + (-t311 + t486) * MDP(16) + (t469 - t271) * MDP(17) + t495 * MDP(18) + t312 * MDP(19) + (t265 * t340 - t299 * t315 + t375) * MDP(20) + (t264 * t340 + t299 * t313 + t376) * MDP(21) + (0.2e1 * t483 + t474 + t471 + (-t276 + t412) * t315 + t375) * MDP(22) + (-pkin(4) * t486 + t473 + t257 * t340 + (qJD(5) + t276) * t313 + t376) * MDP(23) + (pkin(4) * t271 - t455 * t313) * MDP(24) + (t455 * t258 + (-t276 * t315 + t251 + t488) * pkin(4)) * MDP(25); (t368 * t429 + t395 + t468) * MDP(22) + (-t365 * t429 + t450 - t469) * MDP(23) + (-t311 - t486) * MDP(24) + (qJDD(1) * t368 * MDP(22) + (-qJD(1) * t435 - qJDD(1) * t365) * MDP(23)) * t369 + (t256 * t315 + t258 * t313 + t405 - t493 + (-g(3) - t404) * t366 + t489) * MDP(25);];
tau = t1;
