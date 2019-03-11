% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:45
% EndTime: 2019-03-09 02:54:51
% DurationCPUTime: 4.07s
% Computational Cost: add. (3351->366), mult. (7461->517), div. (0->0), fcn. (5257->8), ass. (0->165)
t379 = sin(qJ(3));
t381 = cos(qJ(3));
t454 = sin(pkin(9));
t455 = cos(pkin(9));
t345 = -t454 * t379 + t455 * t381;
t338 = t345 * qJD(1);
t376 = sin(pkin(10));
t377 = cos(pkin(10));
t314 = -t377 * qJD(3) + t338 * t376;
t380 = cos(qJ(6));
t464 = t380 * t314;
t316 = qJD(3) * t376 + t338 * t377;
t378 = sin(qJ(6));
t463 = t314 * t378 - t316 * t380;
t462 = MDP(8) * (t379 ^ 2 - t381 ^ 2);
t326 = t338 * qJD(3);
t389 = t379 * t455 + t381 * t454;
t459 = t389 * qJD(1);
t332 = qJD(6) + t459;
t348 = t376 * t380 + t377 * t378;
t341 = t348 * qJD(6);
t437 = t380 * t377;
t438 = t376 * t378;
t347 = -t437 + t438;
t461 = -t347 * t326 - t341 * t332;
t460 = MDP(12) * qJ(2);
t340 = t347 * qJD(6);
t431 = -t347 * t459 - t340;
t446 = t326 * t348;
t458 = -t332 * t431 - t446;
t382 = -pkin(1) - pkin(7);
t457 = pkin(8) * t377;
t361 = pkin(3) * t454 + qJ(5);
t456 = pkin(8) + t361;
t384 = qJD(1) ^ 2;
t453 = qJ(2) * t384;
t355 = qJD(1) * t382 + qJD(2);
t421 = qJD(4) * t379;
t422 = qJD(3) * t381;
t313 = t355 * t422 + (-qJ(4) * t422 - t421) * qJD(1);
t420 = qJD(4) * t381;
t423 = qJD(3) * t379;
t385 = -t355 * t423 + (qJ(4) * t423 - t420) * qJD(1);
t272 = t313 * t454 - t455 * t385;
t433 = qJ(4) - t382;
t351 = t433 * t379;
t352 = t433 * t381;
t310 = -t351 * t454 + t455 * t352;
t452 = t272 * t310;
t451 = t272 * t345;
t275 = t316 * t378 + t464;
t450 = t275 * t338;
t449 = t463 * t338;
t327 = qJD(3) * t459;
t448 = t310 * t327;
t447 = t326 * t389;
t445 = t327 * t376;
t444 = t327 * t377;
t407 = qJD(3) * t454;
t408 = qJD(3) * t455;
t336 = t379 * t407 - t381 * t408;
t443 = t332 * t336;
t442 = t459 * t376;
t337 = -t379 * t408 - t381 * t407;
t441 = t337 * t376;
t440 = t345 * t376;
t439 = t345 * t377;
t436 = t381 * t384;
t383 = qJD(3) ^ 2;
t435 = t382 * t383;
t434 = t379 * pkin(3) + qJ(2);
t372 = qJD(1) * qJD(2);
t415 = qJD(1) * qJD(3);
t412 = t381 * t415;
t429 = pkin(3) * t412 + t372;
t265 = pkin(4) * t326 + qJ(5) * t327 - qJD(5) * t338 + t429;
t273 = t455 * t313 + t454 * t385;
t268 = qJD(3) * qJD(5) + t273;
t242 = t376 * t265 + t377 * t268;
t424 = qJD(1) * t381;
t331 = -qJ(4) * t424 + t381 * t355;
t324 = qJD(3) * pkin(3) + t331;
t425 = qJD(1) * t379;
t330 = -qJ(4) * t425 + t355 * t379;
t411 = t455 * t330;
t287 = t454 * t324 + t411;
t280 = qJD(3) * qJ(5) + t287;
t350 = pkin(3) * t425 + qJD(1) * qJ(2) + qJD(4);
t288 = pkin(4) * t459 - qJ(5) * t338 + t350;
t253 = t377 * t280 + t376 * t288;
t416 = pkin(3) * t422 + qJD(2);
t274 = -pkin(4) * t336 - qJ(5) * t337 - qJD(5) * t345 + t416;
t328 = t423 * t433 - t420;
t329 = -qJD(3) * t352 - t421;
t294 = t328 * t454 + t329 * t455;
t248 = t376 * t274 + t377 * t294;
t321 = t454 * t330;
t296 = t331 * t455 - t321;
t297 = pkin(3) * t424 + pkin(4) * t338 + qJ(5) * t459;
t255 = t377 * t296 + t376 * t297;
t302 = pkin(4) * t389 - qJ(5) * t345 + t434;
t311 = -t351 * t455 - t352 * t454;
t262 = t376 * t302 + t377 * t311;
t417 = qJD(6) * t380;
t432 = -t314 * t417 - t327 * t437;
t291 = t348 * t459;
t430 = t341 + t291;
t426 = qJD(1) * t459;
t245 = -pkin(8) * t314 + t253;
t419 = qJD(6) * t245;
t418 = qJD(6) * t345;
t414 = 0.2e1 * qJD(1);
t241 = t377 * t265 - t268 * t376;
t247 = t377 * t274 - t294 * t376;
t252 = -t280 * t376 + t377 * t288;
t254 = -t296 * t376 + t377 * t297;
t261 = t377 * t302 - t311 * t376;
t238 = pkin(8) * t445 + t242;
t240 = pkin(5) * t459 - pkin(8) * t316 + t252;
t406 = -qJD(6) * t240 - t238;
t366 = -pkin(3) * t455 - pkin(4);
t403 = -t291 * t332 + t461;
t286 = t324 * t455 - t321;
t293 = -t455 * t328 + t329 * t454;
t295 = t331 * t454 + t411;
t234 = t240 * t380 - t245 * t378;
t235 = t240 * t378 + t245 * t380;
t402 = t241 * t389 - t252 * t336;
t401 = -t242 * t389 + t253 * t336;
t251 = pkin(5) * t389 - pkin(8) * t439 + t261;
t256 = -pkin(8) * t440 + t262;
t400 = t251 * t380 - t256 * t378;
t399 = t251 * t378 + t256 * t380;
t398 = -t252 * t376 + t253 * t377;
t279 = -qJD(3) * pkin(4) + qJD(5) - t286;
t397 = -t279 * t337 - t451;
t396 = t314 * t377 - t316 * t376;
t395 = -qJD(6) * t316 + t445;
t343 = t456 * t377;
t394 = pkin(5) * t338 + qJD(5) * t376 + qJD(6) * t343 + t457 * t459 + t254;
t342 = t456 * t376;
t393 = pkin(8) * t442 - qJD(5) * t377 + qJD(6) * t342 + t255;
t392 = t332 * t348;
t391 = -t397 - t448;
t390 = t327 * t345 + t336 * t459 - t447;
t388 = -t326 * t361 - t327 * t366 + (-qJD(5) + t279) * t459;
t386 = t273 * t389 + t286 * t337 - t287 * t336 - t451;
t250 = -qJD(6) * t463 - t327 * t348;
t353 = -t377 * pkin(5) + t366;
t333 = t459 ^ 2;
t299 = t347 * t345;
t298 = t348 * t345;
t284 = pkin(5) * t440 + t310;
t269 = -pkin(5) * t442 + t295;
t267 = pkin(5) * t441 + t293;
t260 = t314 * pkin(5) + t279;
t259 = -pkin(5) * t445 + t272;
t258 = t337 * t348 + t417 * t439 - t418 * t438;
t257 = -t337 * t347 - t348 * t418;
t249 = t395 * t378 + t432;
t243 = -pkin(8) * t441 + t248;
t239 = -pkin(5) * t336 - t337 * t457 + t247;
t237 = pkin(5) * t326 + pkin(8) * t444 + t241;
t236 = t380 * t237;
t1 = [0.2e1 * t415 * t462 - t383 * t381 * MDP(10) + t422 * t414 * t460 + (-t381 * t435 + (-qJ(2) * t423 + qJD(2) * t381) * t414) * MDP(13) + (t293 * t338 - t294 * t459 - t311 * t326 - t386 - t448) * MDP(14) + (t273 * t311 - t286 * t293 + t287 * t294 + t350 * t416 + t429 * t434 + t452) * MDP(15) + (t247 * t459 + t261 * t326 + t293 * t314 + t376 * t391 + t402) * MDP(16) + (-t248 * t459 - t262 * t326 + t293 * t316 + t377 * t391 + t401) * MDP(17) + (-t247 * t316 - t248 * t314 + (-t241 * t345 - t252 * t337 + t261 * t327) * t377 + (-t242 * t345 - t253 * t337 + t262 * t327) * t376) * MDP(18) + (t241 * t261 + t242 * t262 + t247 * t252 + t248 * t253 + t279 * t293 + t452) * MDP(19) + (-t249 * t299 - t257 * t463) * MDP(20) + (-t249 * t298 + t250 * t299 - t257 * t275 + t258 * t463) * MDP(21) + (t249 * t389 + t257 * t332 - t299 * t326 + t336 * t463) * MDP(22) + (-t250 * t389 - t258 * t332 + t275 * t336 - t298 * t326) * MDP(23) + (-t443 + t447) * MDP(24) + ((t239 * t380 - t243 * t378) * t332 + t400 * t326 + (-t238 * t378 + t236) * t389 - t234 * t336 + t267 * t275 + t284 * t250 + t259 * t298 + t260 * t258 + (-t235 * t389 - t332 * t399) * qJD(6)) * MDP(25) + (-(t239 * t378 + t243 * t380) * t332 - t399 * t326 - (t237 * t378 + t238 * t380) * t389 + t235 * t336 - t267 * t463 + t284 * t249 - t259 * t299 + t260 * t257 + (-t234 * t389 - t332 * t400) * qJD(6)) * MDP(26) + 0.2e1 * (MDP(6) * qJ(2) + MDP(5)) * t372 + (-0.2e1 * MDP(7) * t412 - t383 * MDP(9) + (qJD(2) * t414 - t435) * MDP(12)) * t379; -t384 * MDP(5) - MDP(6) * t453 + (-t337 * t338 + t390) * MDP(14) + (-qJD(1) * t350 + t386) * MDP(15) + (-t314 * t337 + t376 * t390 - t377 * t426) * MDP(16) + (-t316 * t337 + t376 * t426 + t377 * t390) * MDP(17) + (t396 * t336 + (t314 * t376 + t316 * t377) * qJD(1)) * MDP(18) + ((-qJD(1) * t252 - t401) * t377 + (-qJD(1) * t253 - t402) * t376 + t397) * MDP(19) + (-t345 * t250 - t337 * t275 + t336 * t392 + t347 * t332 * qJD(1) - (-t332 * t340 + t446) * t389) * MDP(25) + (qJD(1) * t392 - t345 * t249 + t337 * t463 - t347 * t443 - t389 * t461) * MDP(26) + (MDP(12) * t379 + MDP(13) * t381) * (-t383 - t384); -t384 * t462 - t436 * t460 + ((t287 - t295) * t338 + (-t286 + t296) * t459 + (-t326 * t454 + t327 * t455) * pkin(3)) * MDP(14) + (t286 * t295 - t287 * t296 + (-t272 * t455 + t273 * t454 - t350 * t424) * pkin(3)) * MDP(15) + (-t252 * t338 - t254 * t459 - t272 * t377 - t295 * t314 + t376 * t388) * MDP(16) + (t253 * t338 + t255 * t459 + t272 * t376 - t295 * t316 + t377 * t388) * MDP(17) + (t254 * t316 + t255 * t314 + (-qJD(5) * t314 - t252 * t459 + t242) * t377 + (qJD(5) * t316 - t253 * t459 - t241) * t376) * MDP(18) + (-t252 * t254 - t253 * t255 + t272 * t366 - t279 * t295 + (-t241 * t376 + t242 * t377) * t361 + t398 * qJD(5)) * MDP(19) + (t249 * t348 - t431 * t463) * MDP(20) + (-t249 * t347 - t250 * t348 - t275 * t431 + t430 * t463) * MDP(21) + (t449 - t458) * MDP(22) + (t403 + t450) * MDP(23) - t332 * t338 * MDP(24) + ((-t342 * t380 - t343 * t378) * t326 + t353 * t250 + t259 * t347 - t234 * t338 - t269 * t275 + (t378 * t393 - t380 * t394) * t332 + t430 * t260) * MDP(25) + (-(-t342 * t378 + t343 * t380) * t326 + t353 * t249 + t259 * t348 + t235 * t338 + t269 * t463 + (t378 * t394 + t380 * t393) * t332 + t431 * t260) * MDP(26) + (MDP(13) * t453 + MDP(7) * t436) * t379; (-t338 ^ 2 - t333) * MDP(14) + (t286 * t338 + t429) * MDP(15) + (-t314 * t338 + t326 * t377) * MDP(16) + (-t316 * t338 - t326 * t376 - t333 * t377) * MDP(17) + (t241 * t377 + t242 * t376 - t279 * t338) * MDP(19) + (t403 - t450) * MDP(25) + (t449 + t458) * MDP(26) + (t376 ^ 2 + t377 ^ 2) * MDP(18) * t327 + (t287 * MDP(15) - MDP(16) * t442 - MDP(18) * t396 + t398 * MDP(19)) * t459; (t316 * t459 - t445) * MDP(16) + (-t314 * t459 - t444) * MDP(17) + (-t314 ^ 2 - t316 ^ 2) * MDP(18) + (t252 * t316 + t253 * t314 + t272) * MDP(19) + (-t332 * t463 + t250) * MDP(25) + (-t332 * t464 + (-t316 * t332 + t395) * t378 + t432) * MDP(26); -t275 ^ 2 * MDP(21) + (t275 * t332 + t432) * MDP(22) + t326 * MDP(24) + (t235 * t332 + t236) * MDP(25) + (t234 * t332 + t260 * t275) * MDP(26) - (MDP(20) * t275 - MDP(21) * t463 + t332 * MDP(23) - t260 * MDP(25)) * t463 + (MDP(23) * t395 - MDP(25) * t419 + MDP(26) * t406) * t380 + (t395 * MDP(22) + (qJD(6) * t314 + t444) * MDP(23) + t406 * MDP(25) + (-t237 + t419) * MDP(26)) * t378;];
tauc  = t1;
