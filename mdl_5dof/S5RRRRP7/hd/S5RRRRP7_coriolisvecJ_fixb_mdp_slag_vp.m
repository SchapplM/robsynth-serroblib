% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRP7
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
%   see S5RRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:52
% EndTime: 2019-12-31 21:57:59
% DurationCPUTime: 3.68s
% Computational Cost: add. (4151->372), mult. (10103->487), div. (0->0), fcn. (6866->6), ass. (0->162)
t362 = cos(qJ(4));
t416 = qJD(4) * t362;
t360 = sin(qJ(3));
t363 = cos(qJ(2));
t450 = cos(qJ(3));
t405 = qJD(1) * t450;
t361 = sin(qJ(2));
t419 = qJD(1) * t361;
t327 = t360 * t419 - t363 * t405;
t433 = t327 * t362;
t462 = t416 + t433;
t359 = sin(qJ(4));
t417 = qJD(4) * t359;
t434 = t327 * t359;
t461 = -t417 - t434;
t430 = t360 * t363;
t337 = t450 * t361 + t430;
t413 = qJD(2) + qJD(3);
t308 = t413 * t337;
t460 = t308 * qJD(1);
t414 = qJD(1) * qJD(2);
t459 = -0.2e1 * t414;
t324 = qJD(4) + t327;
t428 = t362 * t460;
t458 = (t324 * t417 - t428) * pkin(8);
t457 = MDP(5) * (t361 ^ 2 - t363 ^ 2);
t385 = -t360 * t361 + t450 * t363;
t375 = t385 * qJD(3);
t262 = t460 * pkin(3) + (-pkin(8) * t375 + (t361 * pkin(2) - t385 * pkin(8)) * qJD(2)) * qJD(1);
t451 = -pkin(7) - pkin(6);
t345 = t451 * t361;
t339 = qJD(1) * t345;
t446 = qJD(2) * pkin(2);
t332 = t339 + t446;
t409 = qJD(2) * t451;
t397 = qJD(1) * t409;
t333 = t361 * t397;
t334 = t363 * t397;
t346 = t451 * t363;
t341 = qJD(1) * t346;
t404 = t450 * qJD(3);
t418 = qJD(3) * t360;
t269 = t332 * t404 + t450 * t333 + t360 * t334 + t341 * t418;
t329 = -qJD(1) * t430 - t361 * t405;
t355 = -pkin(2) * t363 - pkin(1);
t344 = t355 * qJD(1);
t290 = pkin(3) * t327 + pkin(8) * t329 + t344;
t331 = t450 * t341;
t304 = t360 * t332 - t331;
t293 = t413 * pkin(8) + t304;
t379 = t359 * t262 + t362 * t269 + t290 * t416 - t293 * t417;
t445 = qJ(5) * t460;
t235 = qJD(5) * t324 + t379 + t445;
t398 = -t362 * t262 + t359 * t269 + t290 * t417 + t293 * t416;
t449 = pkin(4) * t460;
t238 = t398 - t449;
t456 = t235 * t362 + t238 * t359;
t305 = t360 * t339 - t331;
t396 = pkin(2) * t418 - t305;
t302 = -pkin(3) * t385 - pkin(8) * t337 + t355;
t315 = t360 * t345 - t450 * t346;
t421 = t359 * t302 + t362 * t315;
t455 = t461 * pkin(4) + t462 * qJ(5) + qJD(5) * t359;
t454 = t450 * t345 + t360 * t346;
t307 = t385 * qJD(2) + t375;
t367 = t307 * qJD(1);
t377 = t362 * t329 - t359 * t413;
t275 = -t377 * qJD(4) + t359 * t367;
t453 = t377 ^ 2;
t452 = t324 ^ 2;
t448 = pkin(4) * t329;
t447 = pkin(2) * qJD(3);
t270 = t332 * t418 + t360 * t333 - t450 * t334 - t341 * t404;
t400 = t362 * t413;
t274 = -qJD(4) * t400 - t329 * t417 - t362 * t367;
t240 = pkin(4) * t275 + qJ(5) * t274 + qJD(5) * t377 + t270;
t444 = t240 * t359;
t330 = t360 * t341;
t303 = t450 * t332 + t330;
t292 = -t413 * pkin(3) - t303;
t310 = -t329 * t359 - t400;
t255 = t310 * pkin(4) + qJ(5) * t377 + t292;
t443 = t255 * t377;
t442 = t274 * t359;
t441 = t275 * t362;
t353 = pkin(2) * t360 + pkin(8);
t440 = t460 * t353;
t439 = t310 * t324;
t438 = t310 * t359;
t437 = t377 * t310;
t436 = t377 * t324;
t435 = t377 * t362;
t432 = t337 * t362;
t431 = t359 * t460;
t364 = qJD(2) ^ 2;
t429 = t361 * t364;
t427 = t363 * t364;
t365 = qJD(1) ^ 2;
t426 = t363 * t365;
t425 = t304 + t455;
t424 = t455 - t396;
t300 = -pkin(3) * t329 + pkin(8) * t327;
t423 = t359 * t300 + t362 * t303;
t291 = pkin(2) * t419 + t300;
t306 = t450 * t339 + t330;
t422 = t359 * t291 + t362 * t306;
t258 = t290 * t362 - t293 * t359;
t415 = qJD(5) - t258;
t412 = t450 * pkin(2);
t410 = t361 * t446;
t408 = t359 * t450;
t407 = t362 * t450;
t403 = t361 * t414;
t402 = pkin(1) * t459;
t401 = t324 * t362;
t399 = pkin(2) * t404;
t259 = t290 * t359 + t293 * t362;
t394 = -t259 * t329 + t270 * t359 + t292 * t416;
t393 = t362 * pkin(4) + t359 * qJ(5);
t392 = pkin(4) * t359 - qJ(5) * t362;
t250 = -pkin(4) * t324 + t415;
t251 = qJ(5) * t324 + t259;
t391 = t250 * t362 - t251 * t359;
t390 = t292 * t327 - t440;
t389 = t300 * t362 - t303 * t359;
t343 = -pkin(3) - t393;
t388 = t251 * t329 - t255 * t433 - t444;
t387 = -t240 * t362 - t250 * t329 + t255 * t417;
t386 = t258 * t329 - t270 * t362 + t292 * t417;
t383 = t307 * t359 + t337 * t416;
t382 = -t307 * t362 + t337 * t417;
t381 = t259 * t324 - t398;
t380 = t344 * t329 - t270;
t273 = pkin(3) * t308 - pkin(8) * t307 + t410;
t340 = t361 * t409;
t342 = t363 * t409;
t278 = t454 * qJD(3) + t450 * t340 + t360 * t342;
t378 = t359 * t273 + t362 * t278 + t302 * t416 - t315 * t417;
t376 = (-t324 * t416 - t431) * pkin(8);
t374 = -t353 * t417 + t362 * t399;
t373 = t462 * t250 + t461 * t251 + t456;
t372 = t391 * qJD(4) + t456;
t371 = -t442 - t441 + (-t435 + t438) * qJD(4);
t370 = ((-t274 - t439) * t362 + (-t275 + t436) * t359) * MDP(19) + (-t377 * t401 - t442) * MDP(18) + (-t310 * t329 - t452 * t359 + t428) * MDP(21) + (t324 * t401 - t329 * t377 + t431) * MDP(20) + (t327 * t413 + t367) * MDP(13) + (-t329 * t413 - t460) * MDP(14) + (-t327 ^ 2 + t329 ^ 2) * MDP(12) + (-MDP(11) * t327 + t324 * MDP(22)) * t329;
t369 = t344 * t327 - t269;
t279 = t315 * qJD(3) + t360 * t340 - t450 * t342;
t354 = -t412 - pkin(3);
t335 = -t412 + t343;
t322 = t329 * qJ(5);
t282 = -pkin(4) * t377 + qJ(5) * t310;
t277 = t392 * t337 - t454;
t264 = pkin(4) * t385 - t302 * t362 + t315 * t359;
t263 = -qJ(5) * t385 + t421;
t257 = -t389 + t448;
t256 = -t322 + t423;
t254 = -t291 * t362 + t306 * t359 + t448;
t253 = -t322 + t422;
t249 = -t274 + t439;
t242 = t392 * t307 + (t393 * qJD(4) - qJD(5) * t362) * t337 + t279;
t241 = -pkin(4) * t308 + t421 * qJD(4) - t273 * t362 + t278 * t359;
t239 = qJ(5) * t308 - qJD(5) * t385 + t378;
t1 = [0.2e1 * t363 * MDP(4) * t403 + t457 * t459 + MDP(6) * t427 - MDP(7) * t429 + (-pkin(6) * t427 + t361 * t402) * MDP(9) + (pkin(6) * t429 + t363 * t402) * MDP(10) + (-t329 * t307 + t337 * t367) * MDP(11) + (-t307 * t327 + t329 * t308 - t337 * t460 + t367 * t385) * MDP(12) + (t355 * t460 + t344 * t308 + (-qJD(1) * t385 + t327) * t410) * MDP(16) + (pkin(2) * t337 * t403 + t344 * t307 - t329 * t410 + t355 * t367) * MDP(17) + (-t274 * t432 + t377 * t382) * MDP(18) + ((-t310 * t362 + t359 * t377) * t307 + (t442 - t441 + (t435 + t438) * qJD(4)) * t337) * MDP(19) + (t274 * t385 - t308 * t377 - t382 * t324 + t337 * t428) * MDP(20) + (t275 * t385 - t308 * t310 - t383 * t324 - t337 * t431) * MDP(21) + (t308 * t324 - t385 * t460) * MDP(22) + (t398 * t385 + t258 * t308 + t279 * t310 - t454 * t275 + ((-qJD(4) * t315 + t273) * t324 + t302 * t460 + t292 * qJD(4) * t337) * t362 + ((-qJD(4) * t302 - t278) * t324 - t315 * t460 + t270 * t337 + t292 * t307) * t359) * MDP(23) + (-t259 * t308 + t270 * t432 + t274 * t454 - t279 * t377 - t382 * t292 - t378 * t324 + t379 * t385 - t421 * t460) * MDP(24) + (t238 * t385 - t241 * t324 + t242 * t310 - t250 * t308 + t383 * t255 - t264 * t460 + t275 * t277 + t337 * t444) * MDP(25) + (-t239 * t310 - t241 * t377 - t263 * t275 - t264 * t274 + t391 * t307 + (-t235 * t359 + t238 * t362 + (-t250 * t359 - t251 * t362) * qJD(4)) * t337) * MDP(26) + (-t235 * t385 + t239 * t324 - t240 * t432 + t242 * t377 + t251 * t308 + t382 * t255 + t263 * t460 + t274 * t277) * MDP(27) + (t235 * t263 + t238 * t264 + t239 * t251 + t240 * t277 + t241 * t250 + t242 * t255) * MDP(28) + (t307 * MDP(13) - t308 * MDP(14) - t279 * MDP(16) - t278 * MDP(17)) * t413; (t240 * t335 - t250 * t254 - t251 * t253 - t424 * t255 + (t250 * t408 + t251 * t407) * t447 + t372 * t353) * MDP(28) + (t335 * t275 + (t255 * t327 - t440) * t359 - t424 * t310 + (-t353 * t416 - t359 * t399 + t254) * t324 + t387) * MDP(25) + (t305 * t413 + (-t327 * t419 - t413 * t418) * pkin(2) + t380) * MDP(16) + t365 * t457 - t361 * MDP(4) * t426 + (t253 * t310 + t254 * t377 + (-t310 * t407 - t377 * t408) * t447 + t371 * t353 + t373) * MDP(26) + (t335 * t274 + (-qJD(4) * t255 + t440) * t362 - t424 * t377 + (-t253 + t374) * t324 + t388) * MDP(27) + (t306 * t413 + (t329 * t419 - t413 * t404) * pkin(2) + t369) * MDP(17) + (t354 * t275 + t390 * t359 + t396 * t310 + ((-qJD(4) * t353 - t291) * t362 + (-t399 + t306) * t359) * t324 + t386) * MDP(23) + (-t354 * t274 + t390 * t362 - t396 * t377 + (-t374 + t422) * t324 + t394) * MDP(24) + t370 + (t365 * t361 * MDP(9) + MDP(10) * t426) * pkin(1); (-pkin(3) * t275 + t292 * t434 - t304 * t310 - t389 * t324 + t376 + t386) * MDP(23) + (pkin(3) * t274 + t292 * t433 + t304 * t377 + t423 * t324 + t394 + t458) * MDP(24) + (t304 * t413 + t380) * MDP(16) + (t371 * pkin(8) + t256 * t310 + t257 * t377 + t373) * MDP(26) + (-t255 * t416 - t256 * t324 + t274 * t343 - t377 * t425 + t388 - t458) * MDP(27) + (t255 * t434 + t257 * t324 + t275 * t343 - t425 * t310 + t376 + t387) * MDP(25) + (t372 * pkin(8) + t240 * t343 - t250 * t257 - t251 * t256 - t425 * t255) * MDP(28) + (t303 * t413 + t369) * MDP(17) + t370; -MDP(18) * t437 + (-t310 ^ 2 + t453) * MDP(19) + t249 * MDP(20) + (-t275 - t436) * MDP(21) + t460 * MDP(22) + (t292 * t377 + t381) * MDP(23) + (t258 * t324 + t292 * t310 - t379) * MDP(24) + (-t282 * t310 + t381 + t443 + 0.2e1 * t449) * MDP(25) + (pkin(4) * t274 - qJ(5) * t275 - (t251 - t259) * t377 + (t250 - t415) * t310) * MDP(26) + (0.2e1 * t445 - t255 * t310 - t282 * t377 + (0.2e1 * qJD(5) - t258) * t324 + t379) * MDP(27) + (-pkin(4) * t238 + qJ(5) * t235 - t250 * t259 + t415 * t251 - t255 * t282) * MDP(28); t249 * MDP(26) + (-t452 - t453) * MDP(27) + (-t251 * t324 + t238 - t443) * MDP(28) + (-t437 - t460) * MDP(25);];
tauc = t1;
