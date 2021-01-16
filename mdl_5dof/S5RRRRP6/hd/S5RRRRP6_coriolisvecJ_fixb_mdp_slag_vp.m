% Calculate Coriolis joint torque vector for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:10:33
% EndTime: 2021-01-16 00:10:43
% DurationCPUTime: 3.90s
% Computational Cost: add. (4113->371), mult. (9973->485), div. (0->0), fcn. (6831->6), ass. (0->163)
t381 = cos(qJ(2));
t466 = pkin(7) + pkin(6);
t362 = t466 * t381;
t354 = qJD(1) * t362;
t378 = sin(qJ(3));
t342 = t378 * t354;
t379 = sin(qJ(2));
t361 = t466 * t379;
t352 = qJD(1) * t361;
t465 = cos(qJ(3));
t315 = -t352 * t465 - t342;
t421 = qJD(3) * t465;
t474 = -pkin(2) * t421 + t315;
t377 = sin(qJ(4));
t434 = qJD(4) * t377;
t422 = qJD(1) * t465;
t436 = qJD(1) * t379;
t339 = t378 * t436 - t381 * t422;
t456 = t339 * t377;
t476 = t434 + t456;
t451 = t378 * t381;
t341 = -qJD(1) * t451 - t379 * t422;
t380 = cos(qJ(4));
t399 = -t378 * t379 + t381 * t465;
t393 = t399 * qJD(3);
t318 = qJD(2) * t399 + t393;
t385 = t318 * qJD(1);
t431 = qJD(2) + qJD(3);
t413 = qJD(4) * t431;
t401 = t341 * t434 + (t385 + t413) * t380;
t322 = -t341 * t377 - t380 * t431;
t336 = qJD(4) + t339;
t458 = t322 * t336;
t475 = t401 - t458;
t432 = qJD(1) * qJD(2);
t473 = -0.2e1 * t432;
t472 = MDP(5) * (t379 ^ 2 - t381 ^ 2);
t343 = t465 * t354;
t314 = -t378 * t352 + t343;
t435 = qJD(3) * t378;
t409 = pkin(2) * t435 - t314;
t471 = t476 * pkin(4);
t309 = -pkin(3) * t341 + pkin(8) * t339;
t298 = pkin(2) * t436 + t309;
t470 = t377 * t298 + t380 * t474;
t469 = -t465 * t361 - t378 * t362;
t468 = qJ(5) * t456 - t380 * qJD(5);
t384 = t377 * t385;
t394 = t380 * t341 - t377 * t431;
t280 = -qJD(4) * t394 + t384;
t467 = t394 ^ 2;
t350 = t379 * t465 + t451;
t319 = t431 * t350;
t308 = t319 * qJD(1);
t464 = pkin(4) * t308;
t463 = t380 * pkin(4);
t462 = -qJ(5) - pkin(8);
t461 = qJD(2) * pkin(2);
t460 = t401 * t377;
t459 = t308 * t380;
t457 = t394 * t336;
t455 = t339 * t380;
t454 = t350 * t377;
t453 = t350 * t380;
t452 = t377 * t308;
t382 = qJD(2) ^ 2;
t450 = t379 * t382;
t374 = t380 * qJ(5);
t327 = -t378 * t361 + t362 * t465;
t320 = t380 * t327;
t449 = t381 * t382;
t383 = qJD(1) ^ 2;
t448 = t381 * t383;
t369 = pkin(2) * t378 + pkin(8);
t447 = -qJ(5) - t369;
t372 = -pkin(2) * t381 - pkin(1);
t360 = t372 * qJD(1);
t296 = pkin(3) * t339 + pkin(8) * t341 + t360;
t344 = -t352 + t461;
t313 = t378 * t344 + t343;
t300 = pkin(8) * t431 + t313;
t263 = t380 * t296 - t300 * t377;
t255 = qJ(5) * t394 + t263;
t254 = pkin(4) * t336 + t255;
t446 = t254 - t255;
t295 = t380 * t298;
t407 = -t341 * pkin(4) + t339 * t374;
t415 = qJD(4) * t447;
t445 = -t380 * t415 + t295 + t407 + (qJD(5) - t474) * t377;
t312 = t344 * t465 - t342;
t416 = t380 * t309 - t312 * t377;
t419 = qJD(4) * t462;
t444 = t377 * qJD(5) - t380 * t419 + t407 + t416;
t443 = -t377 * t415 + t468 + t470;
t440 = t377 * t309 + t380 * t312;
t442 = -t377 * t419 + t440 + t468;
t441 = t471 + t409;
t311 = -pkin(3) * t399 - pkin(8) * t350 + t372;
t438 = t377 * t311 + t320;
t433 = qJD(4) * t380;
t428 = t379 * t461;
t277 = pkin(3) * t319 - pkin(8) * t318 + t428;
t425 = qJD(2) * t466;
t353 = t379 * t425;
t355 = t381 * t425;
t282 = t469 * qJD(3) - t465 * t353 - t378 * t355;
t426 = t377 * t277 + t380 * t282 + t311 * t433;
t423 = t350 * t433;
t299 = -pkin(3) * t431 - t312;
t292 = t299 * t433;
t420 = t379 * t432;
t418 = t322 * pkin(4) + qJD(5);
t417 = pkin(1) * t473;
t414 = t336 * t380;
t268 = t308 * pkin(3) + (-pkin(8) * t393 + (t379 * pkin(2) - pkin(8) * t399) * qJD(2)) * qJD(1);
t410 = qJD(1) * t425;
t345 = t379 * t410;
t346 = t381 * t410;
t273 = t344 * t421 - t345 * t465 - t378 * t346 - t354 * t435;
t411 = t377 * t268 + t380 * t273 + t296 * t433 - t300 * t434;
t274 = t344 * t435 - t378 * t345 + t465 * t346 + t354 * t421;
t370 = -pkin(2) * t465 - pkin(3);
t408 = -t313 + t471;
t264 = t296 * t377 + t300 * t380;
t406 = -t264 * t341 + t274 * t377 + t292;
t256 = -qJ(5) * t322 + t264;
t404 = -t254 * t380 - t256 * t377;
t403 = t299 * t339 - t308 * t369;
t402 = -qJ(5) * t318 - qJD(5) * t350;
t400 = t263 * t341 - t274 * t380 + t299 * t434;
t398 = t318 * t377 + t423;
t397 = t318 * t380 - t350 * t434;
t250 = pkin(4) * t280 + t274;
t396 = qJ(5) * t280 - t411;
t395 = t360 * t341 - t274;
t278 = t299 + t418;
t392 = t250 * t377 - t256 * t341 + (t433 + t455) * t278;
t391 = -t250 * t380 + t254 * t341 + t278 * t476;
t266 = t380 * t268;
t390 = -qJD(4) * t264 - t273 * t377 + t266;
t388 = -qJ(5) * t401 + t390;
t241 = qJD(5) * t394 + t388 + t464;
t243 = -qJD(5) * t322 - t396;
t389 = qJD(4) * t404 - t241 * t377 + t243 * t380 - t254 * t455 - t256 * t456;
t387 = (t475 * t380 + (-t280 + t457) * t377) * MDP(19) + (-t394 * t414 + t460) * MDP(18) + (-t336 ^ 2 * t377 - t322 * t341 + t459) * MDP(21) + (t336 * t414 - t341 * t394 + t452) * MDP(20) + t385 * MDP(13) + (-t339 ^ 2 + t341 ^ 2) * MDP(12) + (-MDP(11) * t339 + t336 * MDP(22)) * t341 + (t339 * MDP(13) + (-qJD(1) * t350 - t341) * MDP(14)) * t431;
t386 = t360 * t339 - t273;
t283 = qJD(3) * t327 - t378 * t353 + t355 * t465;
t371 = -pkin(3) - t463;
t359 = pkin(8) * t380 + t374;
t358 = t462 * t377;
t357 = t370 - t463;
t348 = t369 * t380 + t374;
t347 = t447 * t377;
t321 = t322 ^ 2;
t306 = t380 * t311;
t297 = pkin(4) * t454 - t469;
t276 = t380 * t277;
t267 = -qJ(5) * t454 + t438;
t261 = -pkin(4) * t399 - t327 * t377 - t350 * t374 + t306;
t259 = pkin(4) * t398 + t283;
t246 = -qJ(5) * t423 + (-qJD(4) * t327 + t402) * t377 + t426;
t244 = pkin(4) * t319 - t282 * t377 + t276 + t402 * t380 + (-t320 + (qJ(5) * t350 - t311) * t377) * qJD(4);
t1 = [0.2e1 * t381 * MDP(4) * t420 + t472 * t473 + MDP(6) * t449 - MDP(7) * t450 + (-pkin(6) * t449 + t379 * t417) * MDP(9) + (pkin(6) * t450 + t381 * t417) * MDP(10) + (-t341 * t318 + t350 * t385) * MDP(11) + (-t350 * t308 - t318 * t339 + t341 * t319 + t385 * t399) * MDP(12) + (t372 * t308 + t360 * t319 + (-qJD(1) * t399 + t339) * t428) * MDP(16) + (pkin(2) * t350 * t420 + t360 * t318 - t341 * t428 + t372 * t385) * MDP(17) + (-t394 * t397 + t401 * t453) * MDP(18) + ((-t322 * t380 + t377 * t394) * t318 + (-t460 - t280 * t380 + (t322 * t377 + t380 * t394) * qJD(4)) * t350) * MDP(19) + (t308 * t453 - t319 * t394 + t336 * t397 - t399 * t401) * MDP(20) + (t280 * t399 - t319 * t322 - t336 * t398 - t350 * t452) * MDP(21) + (-t308 * t399 + t319 * t336) * MDP(22) + ((-t327 * t433 + t276) * t336 + t306 * t308 - (-t300 * t433 + t266) * t399 + t263 * t319 + t283 * t322 - t469 * t280 + t350 * t292 + ((-qJD(4) * t311 - t282) * t336 - t327 * t308 - (-qJD(4) * t296 - t273) * t399 + t274 * t350 + t299 * t318) * t377) * MDP(23) + (-(-t327 * t434 + t426) * t336 - t438 * t308 + t411 * t399 - t264 * t319 - t283 * t394 - t469 * t401 + t274 * t453 + t397 * t299) * MDP(24) + (-t241 * t399 + t244 * t336 + t250 * t454 + t254 * t319 + t259 * t322 + t261 * t308 + t278 * t398 + t280 * t297) * MDP(25) + (t243 * t399 - t246 * t336 + t250 * t453 - t256 * t319 - t259 * t394 - t267 * t308 + t278 * t397 + t297 * t401) * MDP(26) + (t244 * t394 - t246 * t322 - t261 * t401 - t267 * t280 + t404 * t318 + (-t241 * t380 - t243 * t377 + (t254 * t377 - t256 * t380) * qJD(4)) * t350) * MDP(27) + (t241 * t261 + t243 * t267 + t244 * t254 + t246 * t256 + t250 * t297 + t259 * t278) * MDP(28) + (t318 * MDP(13) - t319 * MDP(14) - MDP(16) * t283 - MDP(17) * t282) * t431; t387 + (t315 * t431 + (t341 * t436 - t421 * t431) * pkin(2) + t386) * MDP(17) + (-t280 * t348 + t322 * t443 - t347 * t401 - t394 * t445 + t389) * MDP(27) + (-t308 * t348 + t336 * t443 + t357 * t401 - t394 * t441 + t392) * MDP(26) + (t314 * t431 + (-t339 * t436 - t431 * t435) * pkin(2) + t395) * MDP(16) + (t370 * t280 + t403 * t377 + t409 * t322 + (-t369 * t433 + t474 * t377 - t295) * t336 + t400) * MDP(23) + (t370 * t401 + t403 * t380 - t409 * t394 + (t369 * t434 + t470) * t336 + t406) * MDP(24) + (t280 * t357 + t308 * t347 + t322 * t441 - t336 * t445 + t391) * MDP(25) + (t241 * t347 + t243 * t348 + t250 * t357 - t254 * t445 - t256 * t443 + t278 * t441) * MDP(28) + t383 * t472 - t379 * MDP(4) * t448 + (MDP(9) * t379 * t383 + MDP(10) * t448) * pkin(1); t387 + (t312 * t431 + t386) * MDP(17) + (t280 * t371 + t308 * t358 + t322 * t408 - t336 * t444 + t391) * MDP(25) + (-t308 * t359 + t336 * t442 + t371 * t401 - t394 * t408 + t392) * MDP(26) + (-pkin(3) * t280 - t416 * t336 - t313 * t322 + t299 * t456 + (-t433 * t336 - t452) * pkin(8) + t400) * MDP(23) + (-pkin(3) * t401 + t440 * t336 + t313 * t394 + t299 * t455 + (t434 * t336 - t459) * pkin(8) + t406) * MDP(24) + (-t280 * t359 + t322 * t442 - t358 * t401 - t394 * t444 + t389) * MDP(27) + (t241 * t358 + t243 * t359 + t250 * t371 - t254 * t444 - t256 * t442 + t278 * t408) * MDP(28) + (t313 * t431 + t395) * MDP(16); -t394 * t322 * MDP(18) + (-t321 + t467) * MDP(19) + (t401 + t458) * MDP(20) + (-t280 - t457) * MDP(21) + t308 * MDP(22) + (t264 * t336 + t299 * t394 + t390) * MDP(23) + (t263 * t336 + t299 * t322 - t411) * MDP(24) + (0.2e1 * t464 + t256 * t336 - (-t278 - t418) * t394 + t388) * MDP(25) + (-pkin(4) * t467 + t255 * t336 + (qJD(5) + t278) * t322 + t396) * MDP(26) + (-pkin(4) * t401 - t322 * t446) * MDP(27) + (t446 * t256 + (t278 * t394 + t241) * pkin(4)) * MDP(28); t475 * MDP(26) + (-t321 - t467) * MDP(27) + (-t254 * t394 + t256 * t322 + t250) * MDP(28) + (-t341 * t433 + t377 * t413 + t384 - t457) * MDP(25);];
tauc = t1;
