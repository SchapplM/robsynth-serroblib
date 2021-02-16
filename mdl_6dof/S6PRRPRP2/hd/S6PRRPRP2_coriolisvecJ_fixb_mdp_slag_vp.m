% Calculate Coriolis joint torque vector for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:57
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:55:01
% EndTime: 2021-01-16 02:55:12
% DurationCPUTime: 4.66s
% Computational Cost: add. (4321->429), mult. (11029->571), div. (0->0), fcn. (8308->10), ass. (0->175)
t364 = sin(qJ(3));
t367 = cos(qJ(3));
t459 = qJ(4) + pkin(8);
t408 = qJD(3) * t459;
t329 = qJD(4) * t367 - t364 * t408;
t360 = sin(pkin(11));
t378 = -qJD(4) * t364 - t367 * t408;
t457 = cos(pkin(11));
t290 = t457 * t329 + t360 * t378;
t405 = t457 * t367;
t385 = -t360 * t364 + t405;
t368 = cos(qJ(2));
t361 = sin(pkin(6));
t429 = qJD(1) * t361;
t415 = t368 * t429;
t313 = t385 * t415;
t432 = t290 - t313;
t349 = qJD(2) * t405;
t427 = qJD(2) * t364;
t332 = t360 * t427 - t349;
t328 = qJD(5) + t332;
t468 = MDP(6) * (t364 ^ 2 - t367 ^ 2);
t341 = t360 * t367 + t457 * t364;
t433 = t329 * t360 - t341 * t415 - t457 * t378;
t334 = t341 * qJD(3);
t337 = t385 * qJD(3);
t425 = qJD(3) * t364;
t418 = pkin(3) * t425;
t292 = pkin(4) * t334 - pkin(9) * t337 + t418;
t356 = -pkin(3) * t367 - pkin(2);
t301 = -pkin(4) * t385 - pkin(9) * t341 + t356;
t345 = t459 * t364;
t346 = t459 * t367;
t315 = -t360 * t345 + t457 * t346;
t363 = sin(qJ(5));
t366 = cos(qJ(5));
t365 = sin(qJ(2));
t416 = t365 * t429;
t423 = qJD(5) * t366;
t424 = qJD(5) * t363;
t467 = -t301 * t423 + t315 * t424 - t432 * t366 + (-t292 + t416) * t363;
t431 = t363 * t301 + t366 * t315;
t343 = qJD(2) * pkin(8) + t416;
t362 = cos(pkin(6));
t428 = qJD(1) * t362;
t414 = t364 * t428;
t389 = -t343 * t367 - t414;
t395 = qJD(4) + t415;
t411 = qJD(3) * t367 * qJ(4);
t466 = t389 * qJD(3) + (-t395 * t364 - t411) * qJD(2);
t335 = t341 * qJD(2);
t465 = qJD(3) * t335;
t402 = qJ(4) * qJD(2) + t343;
t311 = t402 * t367 + t414;
t302 = t360 * t311;
t350 = t367 * t428;
t310 = -t402 * t364 + t350;
t306 = qJD(3) * pkin(3) + t310;
t265 = t457 * t306 - t302;
t262 = -qJD(3) * pkin(4) - t265;
t422 = t366 * qJD(3);
t316 = t335 * t363 - t422;
t318 = qJD(3) * t363 + t335 * t366;
t250 = t316 * pkin(5) - t318 * qJ(6) + t262;
t326 = qJD(2) * t334;
t353 = pkin(3) * t360 + pkin(9);
t445 = t353 * t326;
t464 = t328 * t250 - t445;
t420 = qJD(2) * qJD(3);
t409 = t364 * t420;
t347 = t360 * t409;
t379 = qJD(3) * t349 - t347;
t285 = t318 * qJD(5) + t363 * t379;
t463 = t318 ^ 2;
t462 = t328 ^ 2;
t461 = pkin(3) * t364;
t460 = pkin(5) * t326;
t458 = qJD(2) * pkin(2);
t456 = qJ(6) * t326;
t283 = (-t343 * t364 + t350) * qJD(3) + (-qJ(4) * t425 + t395 * t367) * qJD(2);
t253 = t283 * t360 - t457 * t466;
t284 = -qJD(5) * t422 + t335 * t424 - t366 * t379;
t238 = pkin(5) * t285 + qJ(6) * t284 - qJD(6) * t318 + t253;
t455 = t238 * t363;
t406 = t457 * t311;
t266 = t360 * t306 + t406;
t263 = qJD(3) * pkin(9) + t266;
t327 = t356 * qJD(2) + qJD(4) - t415;
t279 = pkin(4) * t332 - pkin(9) * t335 + t327;
t248 = t263 * t366 + t279 * t363;
t242 = qJ(6) * t328 + t248;
t454 = t242 * t328;
t453 = t248 * t328;
t452 = t284 * t363;
t451 = t316 * t332;
t450 = t316 * t335;
t449 = t318 * t316;
t404 = t318 * t328;
t448 = t318 * t335;
t447 = t328 * t363;
t446 = t341 * t366;
t444 = t361 * t365;
t443 = t361 * t368;
t370 = qJD(2) ^ 2;
t442 = t361 * t370;
t321 = t363 * t326;
t369 = qJD(3) ^ 2;
t441 = t364 * t369;
t322 = t366 * t326;
t440 = t367 * t369;
t439 = qJ(6) * t334 - qJD(6) * t385 - t467;
t293 = t313 * t363 - t366 * t416;
t438 = -pkin(5) * t334 + t431 * qJD(5) + t290 * t363 - t292 * t366 - t293;
t396 = pkin(5) * t363 - qJ(6) * t366;
t397 = t366 * pkin(5) + t363 * qJ(6);
t437 = t396 * t337 + (t397 * qJD(5) - qJD(6) * t366) * t341 + t433;
t268 = t310 * t360 + t406;
t436 = qJD(6) * t363 - t328 * t396 + t268;
t270 = t457 * t310 - t302;
t419 = pkin(3) * t427;
t291 = pkin(4) * t335 + pkin(9) * t332 + t419;
t435 = t366 * t270 + t363 * t291;
t434 = -t363 * t285 - t316 * t423;
t426 = qJD(2) * t365;
t413 = t361 * t426;
t331 = pkin(3) * t409 + qJD(1) * t413;
t247 = -t263 * t363 + t279 * t366;
t421 = qJD(6) - t247;
t417 = t365 * t442;
t412 = qJD(2) * t443;
t410 = t353 * t424;
t407 = t457 * t283;
t314 = t457 * t345 + t346 * t360;
t254 = t466 * t360 + t407;
t276 = t326 * pkin(4) - t379 * pkin(9) + t331;
t401 = t363 * t254 + t263 * t423 - t366 * t276 + t279 * t424;
t400 = t367 * t412;
t399 = t364 * t412;
t354 = -t457 * pkin(3) - pkin(4);
t241 = -pkin(5) * t328 + t421;
t394 = t241 * t366 - t242 * t363;
t393 = t253 * t341 - t315 * t326;
t392 = t321 + (t332 * t366 + t423) * t328;
t391 = -t328 * t424 - t332 * t447 + t322;
t338 = t362 * t364 + t367 * t444;
t390 = t362 * t367 - t364 * t444;
t297 = t457 * t338 + t360 * t390;
t281 = t297 * t363 + t366 * t443;
t282 = t297 * t366 - t363 * t443;
t388 = t337 * t363 + t341 * t423;
t387 = -t337 * t366 + t341 * t424;
t386 = qJD(2) * t458;
t384 = t250 * t318 + t401;
t383 = t366 * t254 - t263 * t424 + t363 * t276 + t279 * t423;
t381 = t328 * t262 - t445;
t377 = t338 * qJD(3);
t375 = -0.2e1 * qJD(3) * t458;
t371 = -t377 - t399;
t339 = -t397 + t354;
t309 = t390 * qJD(3) + t400;
t296 = t338 * t360 - t457 * t390;
t278 = pkin(5) * t318 + qJ(6) * t316;
t271 = t396 * t341 + t314;
t269 = t457 * t309 + t360 * t371;
t267 = t309 * t360 - t457 * t371;
t260 = pkin(5) * t385 - t301 * t366 + t315 * t363;
t259 = -qJ(6) * t385 + t431;
t256 = t316 * t328 - t284;
t249 = -pkin(5) * t335 + t270 * t363 - t291 * t366;
t246 = qJ(6) * t335 + t435;
t245 = t282 * qJD(5) + t269 * t363 - t366 * t413;
t244 = -t281 * qJD(5) + t269 * t366 + t363 * t413;
t237 = t401 - t460;
t236 = qJD(6) * t328 + t383 + t456;
t1 = [-MDP(3) * t417 - t368 * MDP(4) * t442 + (-t367 * t417 + (-t377 - 0.2e1 * t399) * qJD(3)) * MDP(10) + (t364 * t417 + (-t309 - t400) * qJD(3)) * MDP(11) + (-qJD(3) * t267 + (-t326 * t368 + t332 * t426) * t361) * MDP(12) + (-t269 * qJD(3) + (t335 * t426 - t368 * t379) * t361) * MDP(13) + (t267 * t335 - t269 * t332 + t296 * t379 - t297 * t326) * MDP(14) + (t253 * t296 + t254 * t297 - t265 * t267 + t266 * t269 + (t327 * t426 - t331 * t368) * t361) * MDP(15) + (-t244 * t316 + t245 * t318 - t281 * t284 - t282 * t285) * MDP(24) + (t236 * t282 + t237 * t281 + t238 * t296 + t241 * t245 + t242 * t244 + t250 * t267) * MDP(26) + (MDP(21) + MDP(23)) * (-t245 * t328 + t267 * t316 - t281 * t326 + t296 * t285) + (-MDP(22) + MDP(25)) * (t244 * t328 - t267 * t318 + t282 * t326 + t284 * t296); 0.2e1 * t367 * MDP(5) * t409 - 0.2e1 * t420 * t468 + MDP(7) * t440 - MDP(8) * t441 + (-pkin(8) * t440 + t364 * t375) * MDP(10) + (pkin(8) * t441 + t367 * t375) * MDP(11) + (-t332 * t416 + t326 * t356 + t327 * t334 - t331 * t385 + (t332 * t461 - t433) * qJD(3)) * MDP(12) + (-t335 * t416 + t327 * t337 + t331 * t341 - t356 * t347 + (t335 * t461 + t356 * t349 - t432) * qJD(3)) * MDP(13) + (t254 * t385 - t265 * t337 - t266 * t334 + t314 * t379 - t432 * t332 + t433 * t335 + t393) * MDP(14) + (t253 * t314 + t254 * t315 + t331 * t356 + (-t416 + t418) * t327 + t432 * t266 - t433 * t265) * MDP(15) + (-t284 * t446 - t387 * t318) * MDP(16) + ((-t316 * t366 - t318 * t363) * t337 + (t452 - t285 * t366 + (t316 * t363 - t318 * t366) * qJD(5)) * t341) * MDP(17) + (t284 * t385 + t318 * t334 + t341 * t322 - t387 * t328) * MDP(18) + (t285 * t385 - t316 * t334 - t341 * t321 - t388 * t328) * MDP(19) + (-t326 * t385 + t328 * t334) * MDP(20) + (t401 * t385 + t247 * t334 + t314 * t285 + t293 * t328 + t433 * t316 + ((-qJD(5) * t315 + t292) * t328 + t301 * t326 + t262 * qJD(5) * t341) * t366 + ((-qJD(5) * t301 - t290) * t328 + t262 * t337 + t393) * t363) * MDP(21) + (-t248 * t334 + t253 * t446 - t387 * t262 - t314 * t284 + t433 * t318 - t431 * t326 + t467 * t328 + t383 * t385) * MDP(22) + (t237 * t385 - t241 * t334 + t388 * t250 - t260 * t326 + t271 * t285 + t437 * t316 - t438 * t328 + t341 * t455) * MDP(23) + (-t259 * t285 - t260 * t284 + t394 * t337 + t438 * t318 - t439 * t316 + (-t236 * t363 + t237 * t366 + (-t241 * t363 - t242 * t366) * qJD(5)) * t341) * MDP(24) + (-t236 * t385 - t238 * t446 + t242 * t334 + t387 * t250 + t259 * t326 + t271 * t284 - t437 * t318 + t439 * t328) * MDP(25) + (t236 * t259 + t237 * t260 + t238 * t271 + t438 * t241 + t439 * t242 + t437 * t250) * MDP(26); t370 * t468 + t364 * MDP(10) * t386 + (qJD(3) * t268 - t327 * t335 - t332 * t419 - t253) * MDP(12) + (-t407 + t327 * t332 + (-t360 * t389 + t270) * qJD(3) + (t360 * t411 + (-pkin(3) * t335 + t360 * t395) * t364) * qJD(2)) * MDP(13) + ((-t268 + t266) * t335 + (-t265 + t270) * t332 + (-t360 * t326 - t457 * t379) * pkin(3)) * MDP(14) + (t265 * t268 - t266 * t270 + (-t457 * t253 + t254 * t360 - t327 * t427) * pkin(3)) * MDP(15) + (t366 * t404 - t452) * MDP(16) + ((-t284 - t451) * t366 - t318 * t447 + t434) * MDP(17) + (t392 - t448) * MDP(18) + (t391 + t450) * MDP(19) - t328 * t335 * MDP(20) + (-t247 * t335 - t268 * t316 + t354 * t285 + (-t253 + (-qJD(5) * t353 - t291) * t328) * t366 + (t270 * t328 + t381) * t363) * MDP(21) + (t248 * t335 + t253 * t363 - t268 * t318 - t354 * t284 + (t410 + t435) * t328 + t381 * t366) * MDP(22) + (-t238 * t366 + t241 * t335 + t285 * t339 + (-t353 * t423 + t249) * t328 - t436 * t316 + t464 * t363) * MDP(23) + (t246 * t316 - t249 * t318 + (t241 * t332 - t285 * t353 + t236 + (t318 * t353 + t241) * qJD(5)) * t366 + (-t242 * t332 - t284 * t353 + t237 + (t316 * t353 - t242) * qJD(5)) * t363) * MDP(24) + (-t455 - t242 * t335 + t284 * t339 + (-t246 - t410) * t328 + t436 * t318 - t464 * t366) * MDP(25) + (t238 * t339 - t241 * t249 - t242 * t246 - t436 * t250 + (t394 * qJD(5) + t236 * t366 + t237 * t363) * t353) * MDP(26) + (-t364 * t370 * MDP(5) + t386 * MDP(11)) * t367; 0.2e1 * MDP(12) * t465 + (-t347 + (t349 - t332) * qJD(3)) * MDP(13) + (-t332 ^ 2 - t335 ^ 2) * MDP(14) + (t265 * t335 + t266 * t332 + t331) * MDP(15) + (t391 - t450) * MDP(21) + (-t462 * t366 - t321 - t448) * MDP(22) + (-t328 * t447 + t322 - t450) * MDP(23) + ((t284 - t451) * t366 + t363 * t404 + t434) * MDP(24) + (t392 + t448) * MDP(25) + (-t250 * t335 + (-t237 + t454) * t366 + (t328 * t241 + t236) * t363) * MDP(26); MDP(16) * t449 + (-t316 ^ 2 + t463) * MDP(17) + t256 * MDP(18) + (-t285 + t404) * MDP(19) + t326 * MDP(20) + (-t262 * t318 - t401 + t453) * MDP(21) + (t247 * t328 + t262 * t316 - t383) * MDP(22) + (-t278 * t316 - t384 + t453 + 0.2e1 * t460) * MDP(23) + (pkin(5) * t284 - qJ(6) * t285 + (t242 - t248) * t318 + (t241 - t421) * t316) * MDP(24) + (0.2e1 * t456 - t250 * t316 + t278 * t318 + (0.2e1 * qJD(6) - t247) * t328 + t383) * MDP(25) + (-pkin(5) * t237 + qJ(6) * t236 - t241 * t248 + t421 * t242 - t250 * t278) * MDP(26); (t449 - t465) * MDP(23) + t256 * MDP(24) + (-t462 - t463) * MDP(25) + (t384 - t454 - t460) * MDP(26);];
tauc = t1;
