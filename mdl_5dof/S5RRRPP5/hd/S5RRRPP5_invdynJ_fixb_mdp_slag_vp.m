% Calculate vector of inverse dynamics joint torques for
% S5RRRPP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:35
% EndTime: 2019-12-31 20:58:41
% DurationCPUTime: 3.74s
% Computational Cost: add. (2485->382), mult. (5394->424), div. (0->0), fcn. (3513->8), ass. (0->172)
t400 = sin(qJ(3));
t401 = sin(qJ(2));
t403 = cos(qJ(2));
t491 = cos(qJ(3));
t334 = t400 * t403 + t401 * t491;
t323 = t334 * qJD(1);
t389 = t403 * pkin(2);
t504 = t389 + pkin(1);
t343 = t504 * qJD(1);
t506 = -qJ(4) * t323 - t343;
t394 = qJDD(2) + qJDD(3);
t395 = qJD(2) + qJD(3);
t500 = t394 * qJ(4) + t395 * qJD(4);
t457 = qJDD(1) * t403;
t459 = qJD(1) * qJD(2);
t505 = t401 * t459 - t457;
t492 = pkin(7) + pkin(6);
t300 = t395 * t334;
t443 = qJDD(1) * t491;
t458 = qJDD(1) * t401;
t429 = t400 * t458 - t403 * t443;
t280 = qJD(1) * t300 + t429;
t451 = t491 * t403;
t437 = qJD(1) * t451;
t463 = qJD(1) * t401;
t321 = t400 * t463 - t437;
t503 = t280 * qJ(5) + t321 * qJD(5);
t344 = t492 * t401;
t345 = t492 * t403;
t306 = -t400 * t344 + t491 * t345;
t501 = 0.2e1 * t500;
t398 = qJ(2) + qJ(3);
t387 = sin(t398);
t388 = cos(t398);
t466 = t388 * pkin(3) + t387 * qJ(4);
t402 = sin(qJ(1));
t404 = cos(qJ(1));
t444 = g(1) * t402 - g(2) * t404;
t392 = g(1) * t404;
t465 = g(2) * t402 + t392;
t447 = t491 * qJD(3);
t499 = -t491 * qJD(2) - t447;
t337 = qJD(1) * t344;
t487 = qJD(2) * pkin(2);
t331 = -t337 + t487;
t339 = qJD(1) * t345;
t471 = t400 * t339;
t293 = t491 * t331 - t471;
t498 = qJD(4) - t293;
t385 = t394 * pkin(3);
t497 = qJDD(4) - t385;
t474 = t387 * t404;
t475 = t387 * t402;
t496 = g(1) * t474 + g(2) * t475 - g(3) * t388;
t360 = t447 * pkin(2) + qJD(4);
t372 = pkin(2) * t400 + qJ(4);
t495 = t360 * t395 + t372 * t394 + t500;
t494 = t321 ^ 2;
t320 = t323 ^ 2;
t493 = pkin(3) + pkin(4);
t490 = pkin(2) * t401;
t489 = pkin(4) * t323;
t373 = t388 * pkin(4);
t486 = qJ(4) * t280;
t399 = qJDD(1) * pkin(1);
t484 = t280 * t372;
t483 = t293 * t395;
t330 = t491 * t339;
t294 = t400 * t331 + t330;
t482 = t294 * t395;
t297 = -t400 * t337 + t330;
t481 = t297 * t395;
t298 = -t491 * t337 - t471;
t480 = t298 * t395;
t479 = t321 * qJ(5);
t478 = t321 * t323;
t477 = t321 * t395;
t473 = t388 * t402;
t472 = t388 * t404;
t470 = t400 * t401;
t469 = qJ(5) - t492;
t291 = t323 * pkin(3) + t321 * qJ(4);
t315 = t323 * qJ(5);
t284 = t315 + t298;
t468 = t360 - t284;
t467 = t360 - t298;
t396 = t401 ^ 2;
t464 = -t403 ^ 2 + t396;
t462 = qJD(3) * t400;
t277 = t315 + t293;
t461 = qJD(4) - t277;
t456 = t401 * t487;
t455 = pkin(2) * t462;
t317 = t505 * pkin(2) - t399;
t454 = t389 + t466;
t452 = qJD(2) * t492;
t450 = t395 * t462;
t445 = t403 * t459;
t442 = t395 * t401;
t303 = qJDD(2) * pkin(2) + t492 * (-t445 - t458);
t304 = t492 * t505;
t441 = t400 * t303 - t491 * t304 + t331 * t447 - t339 * t462;
t440 = -t491 * t303 - t400 * t304 + t331 * t462 + t339 * t447;
t439 = t395 * t437 + t400 * t457 + t401 * t443;
t379 = -t491 * pkin(2) - pkin(3);
t438 = g(2) * (pkin(3) * t472 + qJ(4) * t474 + t404 * t504);
t436 = g(1) * t475 - g(2) * t474;
t435 = g(1) * t473 - g(2) * t472;
t434 = qJ(4) * t334 + t504;
t283 = t297 + t479;
t433 = -t283 + t455;
t432 = -t297 + t455;
t431 = -pkin(3) * t387 - t490;
t278 = t294 + t479;
t430 = t395 * t470;
t287 = t463 * pkin(2) + t291;
t427 = -t504 - t466;
t279 = qJD(1) * t430 - t439;
t256 = t280 * pkin(3) + t279 * qJ(4) - t323 * qJD(4) + t317;
t258 = t441 + t500;
t426 = -0.2e1 * pkin(1) * t459 - pkin(6) * qJDD(2);
t305 = t344 * t491 + t400 * t345;
t259 = t440 + t497;
t338 = t401 * t452;
t340 = t403 * t452;
t268 = -t491 * t338 - t400 * t340 - t344 * t447 - t345 * t462;
t423 = -g(1) * t472 - g(2) * t473 - g(3) * t387 + t441;
t422 = -t440 + t496;
t265 = -qJD(1) * t400 * t442 + t439 + t477;
t421 = MDP(11) * t478 + t265 * MDP(13) - t429 * MDP(14) + (t320 - t494) * MDP(12) + t394 * MDP(15);
t299 = t403 * t499 + t430;
t420 = -qJ(4) * t299 + qJD(4) * t334 - t456;
t254 = -t394 * pkin(4) + t279 * qJ(5) - t323 * qJD(5) + t259;
t407 = qJD(2) ^ 2;
t419 = -pkin(6) * t407 + 0.2e1 * t399 + t444;
t408 = qJD(1) ^ 2;
t418 = pkin(1) * t408 - pkin(6) * qJDD(1) + t465;
t417 = t268 * t395 + t306 * t394 + t436;
t269 = qJD(3) * t306 - t400 * t338 + t491 * t340;
t416 = -t269 * t395 - t305 * t394 + t435;
t253 = -pkin(4) * t280 + qJDD(5) - t256;
t285 = pkin(3) * t321 + t506;
t415 = -t285 * t321 + t423;
t414 = -t343 * t321 - t423;
t413 = t343 * t323 + t422;
t412 = t465 * t493 * t387;
t267 = -t321 * t493 + qJD(5) - t506;
t411 = t267 * t321 + t423 + t503;
t410 = -t285 * t323 + t422 - t497;
t409 = -t267 * t323 + t254 - t496;
t386 = t395 * qJ(4);
t368 = -pkin(4) + t379;
t349 = qJ(4) * t472;
t347 = qJ(4) * t473;
t346 = pkin(2) * t450;
t333 = -t451 + t470;
t292 = pkin(3) * t333 - t434;
t290 = t386 + t294;
t289 = qJ(5) * t333 + t306;
t288 = -t334 * qJ(5) + t305;
t286 = -pkin(3) * t395 + t498;
t282 = -t333 * t493 + t434;
t276 = -t291 - t489;
t271 = t278 + t386;
t270 = -t287 - t489;
t266 = -t395 * t493 - t315 + t498;
t262 = pkin(3) * t300 - t420;
t261 = t299 * qJ(5) - t334 * qJD(5) + t269;
t260 = qJ(5) * t300 + qJD(5) * t333 + t268;
t257 = -t300 * t493 + t420;
t255 = t258 + t503;
t1 = [qJDD(1) * MDP(1) + t444 * MDP(2) + t465 * MDP(3) + (qJDD(1) * t396 + 0.2e1 * t401 * t445) * MDP(4) + 0.2e1 * (t401 * t457 - t459 * t464) * MDP(5) + (qJDD(2) * t401 + t403 * t407) * MDP(6) + (qJDD(2) * t403 - t401 * t407) * MDP(7) + (t401 * t426 + t403 * t419) * MDP(9) + (-t401 * t419 + t403 * t426) * MDP(10) + (-t279 * t334 - t299 * t323) * MDP(11) + (t279 * t333 - t280 * t334 + t299 * t321 - t300 * t323) * MDP(12) + (-t299 * t395 + t334 * t394) * MDP(13) + (-t300 * t395 - t333 * t394) * MDP(14) + (-t280 * t504 - t300 * t343 + t317 * t333 + t321 * t456 + t416) * MDP(16) + (t279 * t504 + t299 * t343 + t317 * t334 + t323 * t456 - t417) * MDP(17) + (t256 * t333 + t262 * t321 + t280 * t292 + t285 * t300 + t416) * MDP(18) + (-t258 * t333 + t259 * t334 - t268 * t321 + t269 * t323 - t279 * t305 - t280 * t306 - t286 * t299 - t290 * t300 - t465) * MDP(19) + (-t256 * t334 - t262 * t323 + t279 * t292 + t285 * t299 + t417) * MDP(20) + (t258 * t306 + t290 * t268 + t256 * t292 + t285 * t262 + t259 * t305 + t286 * t269 - t492 * t392 - t438 + (-g(1) * t427 - g(2) * t492) * t402) * MDP(21) + (-t253 * t333 - t257 * t321 - t261 * t395 - t267 * t300 - t280 * t282 - t288 * t394 + t435) * MDP(22) + (t253 * t334 + t257 * t323 + t260 * t395 - t267 * t299 - t279 * t282 + t289 * t394 + t436) * MDP(23) + (-t254 * t334 + t255 * t333 + t260 * t321 - t261 * t323 + t266 * t299 + t271 * t300 + t279 * t288 + t280 * t289 + t465) * MDP(24) + (t255 * t289 + t271 * t260 + t254 * t288 + t266 * t261 + t253 * t282 + t267 * t257 - t438 + (g(1) * t469 - g(2) * t373) * t404 + (-g(1) * (t427 - t373) + g(2) * t469) * t402) * MDP(25); (t258 * t372 + t259 * t379 - t285 * t287 - g(1) * (t404 * t431 + t349) - g(2) * (t402 * t431 + t347) - g(3) * t454 + t467 * t290 + t432 * t286) * MDP(21) + (t279 * t368 + t484 + (-t271 - t433) * t323 + (-t266 + t468) * t321) * MDP(24) + MDP(6) * t458 + MDP(7) * t457 + (-g(3) * t403 + t401 * t418) * MDP(9) + (g(3) * t401 + t403 * t418) * MDP(10) + (-t279 * t379 - t484 + (t290 + t432) * t323 + (t286 - t467) * t321) * MDP(19) + (t255 * t372 + t254 * t368 - t267 * t270 - g(1) * (-t404 * t490 + t349) - g(2) * (-t402 * t490 + t347) - g(3) * (t373 + t454) + t412 + t468 * t271 + t433 * t266) * MDP(25) + (t481 + (-t321 * t463 + t394 * t491 - t450) * pkin(2) + t413) * MDP(16) + (t480 + (-t323 * t463 - t394 * t400 - t395 * t447) * pkin(2) + t414) * MDP(17) + (t287 * t323 + t415 - t480 + t495) * MDP(20) + (t270 * t321 + t283 * t395 - t368 * t394 - t346 - t409) * MDP(22) + qJDD(2) * MDP(8) + t421 + (-t287 * t321 - t379 * t394 - t346 + t410 + t481) * MDP(18) + (-t270 * t323 - t284 * t395 + t411 + t495) * MDP(23) + (-t401 * t403 * MDP(4) + MDP(5) * t464) * t408; (t413 + t482) * MDP(16) + (t414 + t483) * MDP(17) + (-t291 * t321 + t385 + t410 + t482) * MDP(18) + (pkin(3) * t279 - t486 + (t290 - t294) * t323 + (t286 - t498) * t321) * MDP(19) + (t291 * t323 + t415 - t483 + t501) * MDP(20) + (t258 * qJ(4) - t259 * pkin(3) - t285 * t291 - t286 * t294 - g(1) * (-pkin(3) * t474 + t349) - g(2) * (-pkin(3) * t475 + t347) - g(3) * t466 + t498 * t290) * MDP(21) + (t276 * t321 + t278 * t395 + t394 * t493 - t409) * MDP(22) + (-t276 * t323 - t277 * t395 + t411 + t501) * MDP(23) + (t486 - t279 * t493 + (-t271 + t278) * t323 + (-t266 + t461) * t321) * MDP(24) + (t255 * qJ(4) - t254 * t493 - t266 * t278 - t267 * t276 - g(1) * t349 - g(2) * t347 - g(3) * (t373 + t466) + t412 + t461 * t271) * MDP(25) + t421; t265 * MDP(19) + (-t290 * t395 - t410) * MDP(21) + (t279 - t477) * MDP(24) + (-t271 * t395 + t409) * MDP(25) + (MDP(18) + MDP(22)) * (-t394 + t478) + (MDP(20) + MDP(23)) * (-t395 ^ 2 - t320); (-t323 * t395 - t429) * MDP(22) + (t439 - t477) * MDP(23) + (-t320 - t494) * MDP(24) + (t266 * t323 - t271 * t321 + t253 + t444) * MDP(25) + (t499 * t401 * MDP(22) + (-MDP(22) * t395 * t403 - MDP(23) * t442) * t400) * qJD(1);];
tau = t1;
