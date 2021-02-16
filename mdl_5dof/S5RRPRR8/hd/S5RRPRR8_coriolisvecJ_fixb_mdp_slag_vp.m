% Calculate Coriolis joint torque vector for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:20
% EndTime: 2021-01-15 21:35:32
% DurationCPUTime: 3.81s
% Computational Cost: add. (3267->319), mult. (8568->437), div. (0->0), fcn. (6557->8), ass. (0->156)
t379 = cos(qJ(5));
t415 = qJD(5) * t379;
t374 = sin(pkin(9));
t375 = cos(pkin(9));
t381 = cos(qJ(2));
t418 = qJD(1) * t381;
t410 = t375 * t418;
t378 = sin(qJ(2));
t419 = qJD(1) * t378;
t341 = -t374 * t419 + t410;
t380 = cos(qJ(4));
t329 = t380 * t341;
t377 = sin(qJ(4));
t352 = t374 * t381 + t375 * t378;
t421 = qJD(1) * t352;
t300 = -t377 * t421 + t329;
t465 = t300 * t379;
t469 = t415 - t465;
t392 = t341 * t377 + t380 * t421;
t376 = sin(qJ(5));
t416 = qJD(5) * t376;
t342 = t352 * qJD(2);
t332 = qJD(1) * t342;
t412 = qJD(1) * qJD(2);
t408 = t378 * t412;
t361 = t374 * t408;
t407 = t381 * t412;
t333 = t375 * t407 - t361;
t417 = qJD(4) * t377;
t268 = qJD(4) * t329 - t377 * t332 + t380 * t333 - t417 * t421;
t371 = qJD(2) + qJD(4);
t425 = t379 * t268 + t371 * t415;
t250 = -t392 * t416 + t425;
t249 = t250 * t379;
t288 = t371 * t376 + t379 * t392;
t442 = t268 * t376;
t251 = t288 * qJD(5) + t442;
t434 = t392 * t376;
t285 = -t379 * t371 + t434;
t468 = -t376 * t251 - t469 * t285 + t249;
t248 = t250 * t376;
t269 = qJD(4) * t392 + t380 * t332 + t333 * t377;
t413 = -qJD(5) + t300;
t265 = t376 * t269;
t426 = -t413 * t415 + t265;
t436 = t300 * t371;
t438 = t392 * t371;
t440 = t288 * t392;
t467 = (-t269 + t438) * MDP(18) - t300 ^ 2 * MDP(16) + (-t300 * MDP(15) + MDP(16) * t392 + MDP(26) * t413) * t392 + (t268 - t436) * MDP(17) + (t469 * t288 + t248) * MDP(22) + (t413 * t465 + t426 - t440) * MDP(24);
t466 = t300 * t376;
t271 = pkin(4) * t392 - pkin(8) * t300;
t446 = -qJ(3) - pkin(6);
t406 = qJD(2) * t446;
t338 = qJD(3) * t381 + t378 * t406;
t320 = t338 * qJD(1);
t339 = -qJD(3) * t378 + t381 * t406;
t321 = t339 * qJD(1);
t284 = -t320 * t374 + t375 * t321;
t275 = -pkin(7) * t333 + t284;
t287 = t375 * t320 + t374 * t321;
t276 = -pkin(7) * t332 + t287;
t360 = t446 * t381;
t357 = qJD(1) * t360;
t346 = t374 * t357;
t359 = t446 * t378;
t356 = qJD(1) * t359;
t445 = qJD(2) * pkin(2);
t350 = t356 + t445;
t302 = t375 * t350 + t346;
t447 = pkin(7) * t421;
t280 = qJD(2) * pkin(3) + t302 - t447;
t433 = t375 * t357;
t303 = t374 * t350 - t433;
t448 = pkin(7) * t341;
t283 = t303 + t448;
t238 = t380 * (qJD(4) * t280 + t276) + t275 * t377 - t283 * t417;
t368 = -pkin(2) * t381 - pkin(1);
t420 = qJD(1) * t368;
t358 = qJD(3) + t420;
t308 = -pkin(3) * t341 + t358;
t463 = -t300 * t308 - t238;
t460 = -0.2e1 * t412;
t458 = MDP(4) * t378;
t457 = MDP(5) * (t378 ^ 2 - t381 ^ 2);
t441 = t285 * t392;
t267 = t379 * t269;
t455 = -t413 * t416 - t267;
t257 = t280 * t377 + t283 * t380;
t239 = qJD(4) * t257 - t380 * t275 + t276 * t377;
t256 = t280 * t380 - t283 * t377;
t254 = -pkin(4) * t371 - t256;
t255 = pkin(8) * t371 + t257;
t258 = -pkin(4) * t300 - pkin(8) * t392 + t308;
t394 = t255 * t376 - t258 * t379;
t454 = -t239 * t379 + t254 * t416 + t392 * t394;
t241 = t255 * t379 + t258 * t376;
t453 = t239 * t376 + t241 * t392 + t254 * t415;
t452 = -t308 * t392 - t239;
t292 = -t338 * t374 + t375 * t339;
t351 = t374 * t378 - t375 * t381;
t345 = t351 * qJD(2);
t277 = pkin(7) * t345 + t292;
t293 = t375 * t338 + t374 * t339;
t278 = -pkin(7) * t342 + t293;
t310 = t375 * t359 + t360 * t374;
t294 = -pkin(7) * t352 + t310;
t311 = t374 * t359 - t375 * t360;
t295 = -pkin(7) * t351 + t311;
t393 = t294 * t380 - t295 * t377;
t242 = qJD(4) * t393 + t277 * t377 + t278 * t380;
t304 = t380 * t351 + t352 * t377;
t305 = -t351 * t377 + t352 * t380;
t323 = pkin(3) * t351 + t368;
t262 = pkin(4) * t304 - pkin(8) * t305 + t323;
t264 = t294 * t377 + t295 * t380;
t272 = -qJD(4) * t304 - t342 * t377 - t345 * t380;
t451 = t239 * t305 + t254 * t272 - t264 * t269 + (qJD(5) * t262 + t242) * t413 - (qJD(5) * t258 + t238) * t304;
t450 = pkin(2) * t374;
t449 = pkin(2) * t378;
t444 = t254 * t305;
t443 = t262 * t269;
t439 = t288 * t376;
t382 = qJD(2) ^ 2;
t432 = t378 * t382;
t431 = t381 * t382;
t383 = qJD(1) ^ 2;
t430 = t381 * t383;
t306 = -t356 * t374 + t433;
t289 = t306 - t448;
t307 = t375 * t356 + t346;
t290 = t307 - t447;
t367 = pkin(2) * t375 + pkin(3);
t389 = t367 * t380 - t377 * t450;
t427 = -t389 * qJD(4) + t289 * t377 + t290 * t380;
t390 = t367 * t377 + t380 * t450;
t424 = t390 * qJD(4) + t289 * t380 - t290 * t377;
t370 = t378 * t445;
t369 = pkin(2) * t419;
t364 = pkin(2) * t408;
t309 = pkin(3) * t332 + t364;
t316 = pkin(3) * t342 + t370;
t315 = pkin(3) * t421 + t369;
t405 = pkin(1) * t460;
t400 = t413 * t376;
t337 = pkin(8) + t390;
t396 = qJD(5) * t337 + t271 + t315;
t395 = -t254 * t300 - t269 * t337;
t391 = -t413 * t466 - t455;
t388 = t272 * t379 - t305 * t416;
t336 = -pkin(4) - t389;
t273 = qJD(4) * t305 + t380 * t342 - t345 * t377;
t246 = pkin(4) * t273 - pkin(8) * t272 + t316;
t245 = pkin(4) * t269 - pkin(8) * t268 + t309;
t244 = t379 * t245;
t243 = qJD(4) * t264 - t277 * t380 + t278 * t377;
t1 = [0.2e1 * t407 * t458 + t457 * t460 + MDP(6) * t431 - MDP(7) * t432 + (-pkin(6) * t431 + t378 * t405) * MDP(9) + (pkin(6) * t432 + t381 * t405) * MDP(10) + (t332 * t368 + t342 * t358 + (t292 + (qJD(1) * t351 - t341) * t449) * qJD(2)) * MDP(11) + (t333 * t368 - t345 * t358 + (0.2e1 * t421 * t449 - t293) * qJD(2)) * MDP(12) + (-t284 * t352 - t287 * t351 - t292 * t421 + t293 * t341 + t302 * t345 - t303 * t342 - t310 * t333 - t311 * t332) * MDP(13) + (t284 * t310 + t287 * t311 + t292 * t302 + t293 * t303 + (t358 + t420) * t370) * MDP(14) + (t268 * t305 + t272 * t392) * MDP(15) + (-t268 * t304 - t269 * t305 + t272 * t300 - t273 * t392) * MDP(16) + (t269 * t323 + t273 * t308 - t300 * t316 + t304 * t309) * MDP(20) + (t268 * t323 + t272 * t308 + t305 * t309 + t316 * t392) * MDP(21) + (t249 * t305 + t288 * t388) * MDP(22) + ((-t285 * t379 - t439) * t272 + (-t248 - t251 * t379 + (t285 * t376 - t288 * t379) * qJD(5)) * t305) * MDP(23) + (t250 * t304 + t267 * t305 + t273 * t288 - t388 * t413) * MDP(24) + (-t305 * t265 - t251 * t304 - t273 * t285 - (-t272 * t376 - t305 * t415) * t413) * MDP(25) + (t269 * t304 - t273 * t413) * MDP(26) + (-t394 * t273 + t243 * t285 + t244 * t304 - t393 * t251 + (-t246 * t413 + t443 + (-t255 * t304 + t264 * t413 + t444) * qJD(5)) * t379 + t451 * t376) * MDP(27) + (-t241 * t273 + t243 * t288 - t393 * t250 + ((-qJD(5) * t264 + t246) * t413 - t443 - (-qJD(5) * t255 + t245) * t304 - qJD(5) * t444) * t376 + t451 * t379) * MDP(28) + (t272 * MDP(17) - t273 * MDP(18) - t243 * MDP(20) - t242 * MDP(21)) * t371; t383 * t457 - t430 * t458 + (t336 * t251 + t395 * t376 + t424 * t285 - (t376 * t427 - t379 * t396) * t413 + t454) * MDP(27) + (t336 * t250 + t395 * t379 + t424 * t288 - (t376 * t396 + t379 * t427) * t413 + t453) * MDP(28) + (t300 * t315 - t371 * t424 + t452) * MDP(20) + (-qJD(2) * t306 + t341 * t369 - t358 * t421 + t284) * MDP(11) + (qJD(2) * t307 - t341 * t358 - t369 * t421 - t287) * MDP(12) + ((t303 + t306) * t421 + (t302 - t307) * t341 + (-t332 * t374 - t333 * t375) * pkin(2)) * MDP(13) + (MDP(9) * t378 * t383 + MDP(10) * t430) * pkin(1) + (-t315 * t392 + t371 * t427 + t463) * MDP(21) + (-t302 * t306 - t303 * t307 + (t284 * t375 + t287 * t374 - t358 * t419) * pkin(2)) * MDP(14) + (t413 * t439 + t468) * MDP(23) + (t391 + t441) * MDP(25) + t467; -t361 * MDP(12) + (-t341 ^ 2 - t421 ^ 2) * MDP(13) + (t302 * t421 - t303 * t341 + t364) * MDP(14) + (t269 + t438) * MDP(20) + (t268 + t436) * MDP(21) + (t391 - t441) * MDP(27) + (-t379 * t413 ^ 2 - t265 - t440) * MDP(28) + ((t374 * t418 + t375 * t419 + t421) * MDP(11) + (t341 + t410) * MDP(12)) * qJD(2); (t257 * t371 + t452) * MDP(20) + (t256 * t371 + t463) * MDP(21) + (t288 * t400 + t468) * MDP(23) + (-t400 * t413 + t267 + t441) * MDP(25) + (-pkin(4) * t251 + (-t256 * t376 + t271 * t379) * t413 - t257 * t285 - t254 * t466 - t426 * pkin(8) + t454) * MDP(27) + (-pkin(4) * t250 - (t256 * t379 + t271 * t376) * t413 - t257 * t288 - t254 * t465 + t455 * pkin(8) + t453) * MDP(28) + t467; t288 * t285 * MDP(22) + (-t285 ^ 2 + t288 ^ 2) * MDP(23) + (-t285 * t413 + t425) * MDP(24) + (-t288 * t413 - t442) * MDP(25) + t269 * MDP(26) + (-t238 * t376 - t241 * t413 - t254 * t288 + t244) * MDP(27) + (-t238 * t379 - t245 * t376 + t254 * t285 + t394 * t413) * MDP(28) + (-MDP(24) * t434 - MDP(25) * t288 - MDP(27) * t241 + MDP(28) * t394) * qJD(5);];
tauc = t1;
