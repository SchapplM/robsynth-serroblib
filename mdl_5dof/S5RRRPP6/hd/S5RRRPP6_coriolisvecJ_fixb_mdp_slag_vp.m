% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:08
% EndTime: 2019-12-31 21:02:14
% DurationCPUTime: 3.99s
% Computational Cost: add. (3153->367), mult. (7854->502), div. (0->0), fcn. (5020->6), ass. (0->166)
t374 = sin(qJ(3));
t417 = qJD(3) * t374;
t377 = cos(qJ(2));
t421 = qJD(1) * t377;
t466 = -t374 * t421 + t417;
t375 = sin(qJ(2));
t422 = qJD(1) * t375;
t403 = t374 * t422;
t376 = cos(qJ(3));
t411 = t376 * qJD(2);
t335 = t403 - t411;
t420 = qJD(2) * t374;
t337 = t376 * t422 + t420;
t372 = sin(pkin(8));
t373 = cos(pkin(8));
t293 = t373 * t335 + t337 * t372;
t357 = -qJD(3) + t421;
t465 = t293 * t357;
t331 = t372 * t376 + t373 * t374;
t321 = t331 * qJD(3);
t433 = t331 * t421 - t321;
t415 = qJD(3) * t376;
t447 = t373 * t376;
t432 = t466 * t372 - t373 * t415 + t421 * t447;
t398 = t375 * t415;
t418 = qJD(2) * t377;
t401 = t374 * t418;
t464 = t398 + t401;
t385 = -t335 * t372 + t373 * t337;
t463 = t385 ^ 2;
t409 = qJD(1) * qJD(2);
t462 = -0.2e1 * t409;
t345 = -qJD(2) * pkin(2) + pkin(6) * t422;
t306 = pkin(3) * t335 + qJD(4) + t345;
t255 = pkin(4) * t293 - qJ(5) * t385 + t306;
t461 = t255 * t385;
t460 = t375 * MDP(4);
t370 = t375 ^ 2;
t459 = (-t377 ^ 2 + t370) * MDP(5);
t441 = t376 * t377;
t382 = pkin(3) * t375 - qJ(4) * t441;
t388 = pkin(2) * t375 - pkin(7) * t377;
t338 = t388 * qJD(1);
t427 = pkin(6) * t403 + t376 * t338;
t287 = t382 * qJD(1) + t427;
t323 = t374 * t338;
t443 = t375 * t376;
t444 = t374 * t377;
t297 = t323 + (-pkin(6) * t443 - qJ(4) * t444) * qJD(1);
t456 = -qJ(4) - pkin(7);
t393 = qJD(3) * t456;
t414 = qJD(4) * t376;
t318 = t374 * t393 + t414;
t381 = -qJD(4) * t374 + t376 * t393;
t434 = (-t287 + t381) * t373 + (t297 - t318) * t372;
t366 = pkin(6) * t421;
t458 = t466 * pkin(3) - t366;
t457 = pkin(6) * t374;
t346 = qJD(2) * pkin(7) + t366;
t342 = -pkin(2) * t377 - pkin(7) * t375 - pkin(1);
t327 = t342 * qJD(1);
t446 = t374 * t327;
t301 = t346 * t376 + t446;
t281 = -qJ(4) * t335 + t301;
t278 = t373 * t281;
t300 = t376 * t327 - t346 * t374;
t280 = -qJ(4) * t337 + t300;
t253 = t280 * t372 + t278;
t455 = t253 * t385;
t454 = t281 * t372;
t416 = qJD(3) * t375;
t399 = t374 * t416;
t395 = t377 * t409;
t408 = qJD(2) * qJD(3);
t428 = (t395 + t408) * t376;
t307 = -qJD(1) * t399 + t428;
t453 = t307 * t374;
t452 = t335 * t357;
t451 = t337 * t357;
t450 = t345 * t374;
t449 = t345 * t376;
t448 = t357 * t376;
t445 = t374 * t375;
t378 = qJD(2) ^ 2;
t442 = t375 * t378;
t440 = t377 * t378;
t379 = qJD(1) ^ 2;
t439 = t377 * t379;
t339 = t388 * qJD(2);
t328 = qJD(1) * t339;
t396 = t375 * t409;
t389 = pkin(6) * t396;
t431 = -t376 * t328 - t374 * t389;
t380 = -t301 * qJD(3) - t431;
t247 = pkin(3) * t396 - qJ(4) * t307 - qJD(4) * t337 + t380;
t394 = t374 * t408;
t308 = qJD(1) * t464 + t394;
t383 = -t327 * t415 - t374 * t328 + t346 * t417;
t252 = -qJ(4) * t308 - qJD(4) * t335 - t376 * t389 - t383;
t438 = -t373 * t247 + t372 * t252;
t236 = t372 * t247 + t373 * t252;
t262 = t372 * t287 + t373 * t297;
t257 = qJ(5) * t422 + t262;
t285 = t373 * t318 + t372 * t381;
t437 = t257 - t285;
t436 = pkin(4) * t422 - t434;
t360 = pkin(6) * t441;
t419 = qJD(2) * t375;
t429 = t376 * t339 + t419 * t457;
t263 = -t375 * t414 + t382 * qJD(2) + (-t360 + (qJ(4) * t375 - t342) * t374) * qJD(3) + t429;
t430 = t374 * t339 + t342 * t415;
t268 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t443 + (-qJD(4) * t375 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t377) * t374 + t430;
t241 = t372 * t263 + t373 * t268;
t435 = pkin(4) * t433 - qJ(5) * t432 + qJD(5) * t331 - t458;
t276 = -pkin(3) * t357 + t280;
t249 = t372 * t276 + t278;
t333 = t376 * t342;
t298 = -qJ(4) * t443 + t333 + (-pkin(3) - t457) * t377;
t425 = t374 * t342 + t360;
t302 = -qJ(4) * t445 + t425;
t270 = t372 * t298 + t373 * t302;
t424 = pkin(3) * t445 + t375 * pkin(6);
t413 = qJD(5) * t357;
t412 = t375 * MDP(15);
t254 = t280 * t373 - t454;
t410 = qJD(5) - t254;
t406 = qJ(5) * t396 + t236;
t405 = pkin(3) * t464 + pkin(6) * t418;
t404 = -pkin(3) * t376 - pkin(2);
t400 = t377 * t411;
t397 = t456 * t374;
t291 = pkin(3) * t308 + pkin(6) * t395;
t392 = pkin(1) * t462;
t274 = t307 * t372 + t373 * t308;
t391 = t335 + t411;
t390 = -t337 + t420;
t275 = t307 * t373 - t308 * t372;
t344 = t456 * t376;
t304 = -t344 * t372 - t373 * t397;
t305 = -t373 * t344 + t372 * t397;
t387 = -t305 * t274 + t275 * t304 - t285 * t293;
t240 = t263 * t373 - t268 * t372;
t248 = t276 * t373 - t454;
t269 = t298 * t373 - t302 * t372;
t330 = t372 * t374 - t447;
t384 = qJD(1) * t370 - t357 * t377;
t234 = -pkin(4) * t396 + t438;
t237 = pkin(4) * t274 - qJ(5) * t275 - qJD(5) * t385 + t291;
t363 = -pkin(3) * t373 - pkin(4);
t361 = pkin(3) * t372 + qJ(5);
t314 = -t372 * t445 + t373 * t443;
t313 = t331 * t375;
t289 = pkin(4) * t330 - qJ(5) * t331 + t404;
t283 = t375 * t321 + t372 * t401 - t373 * t400;
t282 = t330 * t416 - t331 * t418;
t277 = pkin(4) * t313 - qJ(5) * t314 + t424;
t267 = pkin(4) * t377 - t269;
t266 = -qJ(5) * t377 + t270;
t256 = pkin(3) * t337 + pkin(4) * t385 + qJ(5) * t293;
t244 = -qJ(5) * t357 + t249;
t243 = pkin(4) * t357 + qJD(5) - t248;
t242 = -pkin(4) * t282 + qJ(5) * t283 - qJD(5) * t314 + t405;
t239 = -pkin(4) * t419 - t240;
t238 = qJ(5) * t419 - qJD(5) * t377 + t241;
t233 = t406 - t413;
t1 = [0.2e1 * t395 * t460 + t459 * t462 + MDP(6) * t440 - MDP(7) * t442 + (-pkin(6) * t440 + t375 * t392) * MDP(9) + (pkin(6) * t442 + t377 * t392) * MDP(10) + (t307 * t443 + (-t399 + t400) * t337) * MDP(11) + ((-t335 * t376 - t337 * t374) * t418 + (-t453 - t308 * t376 + (t335 * t374 - t337 * t376) * qJD(3)) * t375) * MDP(12) + (t357 * t399 - t307 * t377 + (t337 * t375 + t384 * t376) * qJD(2)) * MDP(13) + (t357 * t398 + t308 * t377 + (-t335 * t375 - t384 * t374) * qJD(2)) * MDP(14) + (-t357 - t421) * qJD(2) * t412 + (-(-t342 * t417 + t429) * t357 + (t345 * t415 + pkin(6) * t308 + (qJD(1) * t333 + t300) * qJD(2)) * t375 + ((pkin(6) * t335 + t450) * qJD(2) + (t446 + (pkin(6) * t357 + t346) * t376) * qJD(3) + t431) * t377) * MDP(16) + ((-pkin(6) * t377 * t417 + t430) * t357 - t383 * t377 + (pkin(6) * t307 - t345 * t417) * t375 + ((pkin(6) * t337 + t449) * t377 + (-pkin(6) * t448 - t425 * qJD(1) - t301) * t375) * qJD(2)) * MDP(17) + (-t236 * t313 - t240 * t385 - t241 * t293 + t248 * t283 + t249 * t282 - t269 * t275 - t270 * t274 + t314 * t438) * MDP(18) + (t236 * t270 + t248 * t240 + t249 * t241 - t269 * t438 + t291 * t424 + t306 * t405) * MDP(19) + (t234 * t377 + t237 * t313 + t239 * t357 + t242 * t293 - t255 * t282 + t274 * t277 + (-qJD(1) * t267 - t243) * t419) * MDP(20) + (-t233 * t313 + t234 * t314 - t238 * t293 + t239 * t385 - t243 * t283 + t244 * t282 - t266 * t274 + t267 * t275) * MDP(21) + (-t233 * t377 - t237 * t314 - t238 * t357 - t242 * t385 + t255 * t283 - t275 * t277 + (qJD(1) * t266 + t244) * t419) * MDP(22) + (t233 * t266 + t234 * t267 + t237 * t277 + t238 * t244 + t239 * t243 + t242 * t255) * MDP(23); -t439 * t460 + t379 * t459 + (-t337 * t448 + t453) * MDP(11) + ((t307 + t452) * t376 + (-t308 + t451) * t374) * MDP(12) + (-t357 * t415 + (t357 * t441 + t390 * t375) * qJD(1)) * MDP(13) + (t357 * t417 + (-t357 * t444 + t391 * t375) * qJD(1)) * MDP(14) + t357 * qJD(1) * t412 + (-pkin(2) * t308 + t427 * t357 + (pkin(7) * t448 + t450) * qJD(3) + ((-pkin(7) * t420 - t300) * t375 + (-t391 * pkin(6) - t450) * t377) * qJD(1)) * MDP(16) + (-pkin(2) * t307 - t323 * t357 + (-pkin(7) * t357 * t374 + t449) * qJD(3) + (-t345 * t441 + (-pkin(7) * t411 + t301) * t375 + (t357 * t443 + t390 * t377) * pkin(6)) * qJD(1)) * MDP(17) + (-t236 * t330 + t432 * t248 + t433 * t249 + t262 * t293 + t331 * t438 - t385 * t434 + t387) * MDP(18) + (t236 * t305 + t438 * t304 + t291 * t404 + t458 * t306 + (t285 - t262) * t249 + t434 * t248) * MDP(19) + (t237 * t330 + t274 * t289 + t436 * t357 - t435 * t293 - t433 * t255 + (-qJD(2) * t304 + t243) * t422) * MDP(20) + (-t233 * t330 + t234 * t331 - t432 * t243 + t433 * t244 + t257 * t293 + t385 * t436 + t387) * MDP(21) + (-t237 * t331 - t275 * t289 + t437 * t357 + t435 * t385 + t432 * t255 + (qJD(2) * t305 - t244) * t422) * MDP(22) + (t233 * t305 + t234 * t304 + t237 * t289 + t436 * t243 - t437 * t244 - t435 * t255) * MDP(23) + (t379 * t375 * MDP(9) + MDP(10) * t439) * pkin(1); t337 * t335 * MDP(11) + (-t335 ^ 2 + t337 ^ 2) * MDP(12) + (t428 - t452) * MDP(13) + (-t394 - t451) * MDP(14) + (-t301 * t357 - t337 * t345 + t380) * MDP(16) + (-t300 * t357 + t335 * t345 + t383) * MDP(17) + (t249 * t385 - t455) * MDP(18) + (t248 * t253 - t249 * t254) * MDP(19) + (-t253 * t357 - t438 - t461) * MDP(20) + (t244 * t385 - t274 * t361 + t275 * t363 - t455) * MDP(21) + (t254 * t357 + t256 * t385 + t406 - 0.2e1 * t413) * MDP(22) + (t233 * t361 + t234 * t363 - t243 * t253 + t410 * t244 - t255 * t256) * MDP(23) + (-MDP(14) * t401 + ((-t374 * MDP(13) - t376 * MDP(14)) * qJD(3) + (MDP(15) + t376 * pkin(6) * MDP(17) + (pkin(4) - t363) * MDP(20) + t361 * MDP(22)) * qJD(2)) * t375) * qJD(1) + ((-t274 * t372 - t275 * t373) * MDP(18) + (t236 * t372 - t306 * t337 - t373 * t438) * MDP(19)) * pkin(3) + ((-t248 + t254) * MDP(18) - t256 * MDP(20) + (t243 - t410) * MDP(21) - t255 * MDP(22)) * t293; (t248 * t385 + t249 * t293 + t291) * MDP(19) + (-t357 * t385 + t274) * MDP(20) + (-t275 - t465) * MDP(22) + (-t243 * t385 + t244 * t293 + t237) * MDP(23) + (MDP(18) + MDP(21)) * (-t293 ^ 2 - t463); (t293 * t385 - t396) * MDP(20) + (t275 - t465) * MDP(21) + (-t357 ^ 2 - t463) * MDP(22) + (t244 * t357 + t234 + t461) * MDP(23);];
tauc = t1;
