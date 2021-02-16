% Calculate Coriolis joint torque vector for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:49:09
% EndTime: 2021-01-15 22:49:21
% DurationCPUTime: 4.05s
% Computational Cost: add. (3720->310), mult. (10029->417), div. (0->0), fcn. (7348->8), ass. (0->161)
t389 = qJD(2) + qJD(3);
t398 = cos(qJ(3));
t399 = cos(qJ(2));
t432 = qJD(1) * qJD(2);
t428 = t399 * t432;
t436 = qJD(1) * t399;
t429 = t398 * t436;
t395 = sin(qJ(3));
t396 = sin(qJ(2));
t437 = qJD(1) * t396;
t430 = t395 * t437;
t332 = qJD(3) * t429 - t389 * t430 + t398 * t428;
t368 = t395 * t399 + t396 * t398;
t338 = t389 * t368;
t333 = t338 * qJD(1);
t392 = sin(pkin(9));
t393 = cos(pkin(9));
t288 = t332 * t392 + t393 * t333;
t289 = t332 * t393 - t333 * t392;
t394 = sin(qJ(5));
t397 = cos(qJ(5));
t354 = -t429 + t430;
t356 = -t395 * t436 - t398 * t437;
t412 = -t354 * t392 - t393 * t356;
t434 = qJD(5) * t394;
t329 = -t393 * t354 + t356 * t392;
t446 = t329 * t397;
t245 = qJD(5) * t446 - t394 * t288 + t397 * t289 - t412 * t434;
t280 = t329 * t394 + t397 * t412;
t283 = -t394 * t412 + t446;
t403 = -qJD(5) * t280 - t397 * t288 - t289 * t394;
t388 = qJD(5) + t389;
t453 = t283 * t388;
t454 = t280 * t388;
t472 = (t403 + t454) * MDP(25) + (t280 ^ 2 - t283 ^ 2) * MDP(23) - t283 * t280 * MDP(22) + (t245 - t453) * MDP(24);
t351 = t356 * qJ(4);
t461 = pkin(6) + pkin(7);
t376 = t461 * t399;
t371 = qJD(1) * t376;
t357 = t395 * t371;
t375 = t461 * t396;
t369 = qJD(1) * t375;
t456 = qJD(2) * pkin(2);
t363 = -t369 + t456;
t422 = t398 * t363 - t357;
t311 = t351 + t422;
t302 = pkin(3) * t389 + t311;
t361 = t398 * t371;
t411 = -t363 * t395 - t361;
t455 = qJ(4) * t354;
t312 = -t411 - t455;
t449 = t393 * t312;
t267 = t392 * t302 + t449;
t458 = pkin(8) * t329;
t253 = t267 + t458;
t385 = -pkin(2) * t399 - pkin(1);
t374 = t385 * qJD(1);
t339 = pkin(3) * t354 + qJD(4) + t374;
t295 = -pkin(4) * t329 + t339;
t424 = t253 * t434 - t295 * t283;
t431 = qJD(2) * t461;
t417 = qJD(1) * t431;
t365 = t399 * t417;
t435 = qJD(3) * t395;
t420 = -t395 * t365 - t371 * t435;
t364 = t396 * t417;
t465 = (qJD(3) * t363 - t364) * t398;
t261 = -qJ(4) * t333 - qJD(4) * t354 + t420 + t465;
t421 = t395 * t364 - t398 * t365;
t405 = t411 * qJD(3) + t421;
t262 = -qJ(4) * t332 + qJD(4) * t356 + t405;
t243 = -t261 * t392 + t393 * t262;
t238 = -pkin(8) * t289 + t243;
t244 = t393 * t261 + t392 * t262;
t239 = -pkin(8) * t288 + t244;
t409 = t397 * t238 - t394 * t239 - t295 * t280;
t469 = -0.2e1 * t432;
t323 = t412 * pkin(8);
t467 = t396 * MDP(4);
t464 = (t396 ^ 2 - t399 ^ 2) * MDP(5);
t419 = t369 * t395 - t361;
t314 = t419 + t455;
t439 = -t398 * t369 - t357;
t315 = t351 + t439;
t448 = t393 * t395;
t457 = pkin(2) * qJD(3);
t441 = t393 * t314 - t315 * t392 - (-t392 * t398 - t448) * t457;
t450 = t392 * t395;
t440 = t392 * t314 + t393 * t315 - (t393 * t398 - t450) * t457;
t463 = qJD(1) * t368;
t462 = qJD(5) - t388;
t460 = pkin(3) * t356;
t459 = pkin(3) * t392;
t452 = t374 * t356;
t451 = t375 * t398;
t303 = t392 * t312;
t400 = qJD(2) ^ 2;
t447 = t396 * t400;
t445 = t399 * t400;
t401 = qJD(1) ^ 2;
t444 = t399 * t401;
t443 = t458 - t441;
t442 = -t323 + t440;
t367 = t395 * t396 - t398 * t399;
t370 = t396 * t431;
t372 = t399 * t431;
t406 = -qJD(3) * t451 - t398 * t370 - t395 * t372 - t376 * t435;
t274 = -qJ(4) * t338 - qJD(4) * t367 + t406;
t337 = t389 * t367;
t410 = t375 * t395 - t376 * t398;
t404 = t410 * qJD(3) + t370 * t395 - t398 * t372;
t275 = qJ(4) * t337 - qJD(4) * t368 + t404;
t250 = t393 * t274 + t392 * t275;
t269 = t393 * t311 - t303;
t324 = -qJ(4) * t368 - t376 * t395 - t451;
t325 = -qJ(4) * t367 - t410;
t279 = t392 * t324 + t393 * t325;
t387 = t396 * t456;
t386 = pkin(2) * t437;
t427 = -pkin(2) * t389 - t363;
t316 = pkin(3) * t333 + qJD(2) * t386;
t331 = pkin(3) * t338 + t387;
t426 = pkin(1) * t469;
t249 = -t274 * t392 + t393 * t275;
t266 = t393 * t302 - t303;
t268 = -t311 * t392 - t449;
t278 = t393 * t324 - t325 * t392;
t384 = pkin(2) * t398 + pkin(3);
t349 = -pkin(2) * t450 + t393 * t384;
t299 = pkin(4) * t412 - t460;
t416 = -t329 * t339 - t244;
t251 = pkin(4) * t389 + t266 - t323;
t415 = -t394 * t251 - t397 * t253;
t414 = t266 * t329 + t267 * t412;
t335 = t393 * t367 + t368 * t392;
t336 = -t367 * t392 + t368 * t393;
t290 = t335 * t397 + t336 * t394;
t291 = -t335 * t394 + t336 * t397;
t341 = pkin(3) * t367 + t385;
t408 = -t339 * t412 + t243;
t407 = t374 * t354 - t420;
t402 = -t356 * t354 * MDP(11) + t332 * MDP(13) + (-t354 ^ 2 + t356 ^ 2) * MDP(12) + (t354 * MDP(13) + (-t356 - t463) * MDP(14)) * t389 + t472;
t383 = pkin(3) * t393 + pkin(4);
t350 = pkin(2) * t448 + t384 * t392;
t344 = pkin(4) + t349;
t340 = t386 - t460;
t307 = pkin(4) * t335 + t341;
t296 = t299 + t386;
t294 = -t337 * t393 - t338 * t392;
t293 = -t337 * t392 + t393 * t338;
t273 = pkin(4) * t293 + t331;
t265 = -pkin(8) * t335 + t279;
t264 = -pkin(8) * t336 + t278;
t263 = pkin(4) * t288 + t316;
t255 = -t323 + t269;
t254 = t268 - t458;
t248 = t291 * qJD(5) + t397 * t293 + t294 * t394;
t247 = -t290 * qJD(5) - t293 * t394 + t294 * t397;
t242 = -pkin(8) * t293 + t250;
t241 = -pkin(8) * t294 + t249;
t1 = [0.2e1 * t428 * t467 + t464 * t469 + MDP(6) * t445 - MDP(7) * t447 + (-pkin(6) * t445 + t396 * t426) * MDP(9) + (pkin(6) * t447 + t399 * t426) * MDP(10) + (t332 * t368 + t337 * t356) * MDP(11) + (-t332 * t367 - t333 * t368 + t337 * t354 + t338 * t356) * MDP(12) + (t385 * t333 + t374 * t338 + (qJD(1) * t367 + t354) * t387) * MDP(16) + (t385 * t332 - t374 * t337 + (-t356 + t463) * t387) * MDP(17) + (t288 * t341 + t293 * t339 + t316 * t335 - t329 * t331) * MDP(18) + (t289 * t341 + t294 * t339 + t316 * t336 + t331 * t412) * MDP(19) + (-t243 * t336 - t244 * t335 - t249 * t412 + t250 * t329 - t266 * t294 - t267 * t293 - t278 * t289 - t279 * t288) * MDP(20) + (t243 * t278 + t244 * t279 + t249 * t266 + t250 * t267 + t316 * t341 + t331 * t339) * MDP(21) + (t245 * t291 + t247 * t280) * MDP(22) + (-t245 * t290 + t247 * t283 - t248 * t280 + t291 * t403) * MDP(23) + (t295 * t248 + t263 * t290 - t273 * t283 - t307 * t403) * MDP(27) + (t307 * t245 + t295 * t247 + t263 * t291 + t273 * t280) * MDP(28) + (-t337 * MDP(13) - t338 * MDP(14) + t404 * MDP(16) - t406 * MDP(17) + t249 * MDP(18) - t250 * MDP(19)) * t389 + (t247 * MDP(24) - t248 * MDP(25) + (t241 * t397 - t242 * t394 + (-t264 * t394 - t265 * t397) * qJD(5)) * MDP(27) + (-t241 * t394 - t242 * t397 - (t264 * t397 - t265 * t394) * qJD(5)) * MDP(28)) * t388; (t296 * t283 + (t442 * t394 + t443 * t397) * t388 + ((-t344 * t394 - t350 * t397) * t388 + t415) * qJD(5) + t409) * MDP(27) + (-t296 * t280 + (-t238 + (qJD(5) * t350 - t443) * t388) * t394 + (-qJD(5) * t251 - t239 + (-qJD(5) * t344 + t442) * t388) * t397 + t424) * MDP(28) + (t356 * t386 + t439 * t389 + (t427 * qJD(3) + t364) * t398 + t407) * MDP(17) + (-t354 * t386 + t452 - t419 * t389 + (t427 * t395 - t361) * qJD(3) + t421) * MDP(16) + (-t288 * t350 - t289 * t349 - t329 * t440 + t412 * t441 + t414) * MDP(20) + (t243 * t349 + t244 * t350 - t441 * t266 - t440 * t267 - t339 * t340) * MDP(21) + (t329 * t340 - t441 * t389 + t408) * MDP(18) + (-t340 * t412 + t440 * t389 + t416) * MDP(19) - t444 * t467 + t402 + t401 * t464 + (t401 * t396 * MDP(9) + MDP(10) * t444) * pkin(1); (-t397 * t239 - t394 * t238 - t299 * t280 + (t254 * t394 + t255 * t397) * t388 + (-(t383 * t397 - t394 * t459) * t388 - t397 * t251) * qJD(5) + t424) * MDP(28) + (t299 * t283 - (t254 * t397 - t255 * t394) * t388 + ((-t383 * t394 - t397 * t459) * t388 + t415) * qJD(5) + t409) * MDP(27) + (-t411 * t389 + t405 + t452) * MDP(16) + (t422 * t389 + t407 - t465) * MDP(17) + (-t268 * t389 - t329 * t460 + t408) * MDP(18) + (t269 * t389 + t412 * t460 + t416) * MDP(19) + (t268 * t412 - t269 * t329 + (-t288 * t392 - t289 * t393) * pkin(3) + t414) * MDP(20) + (-t266 * t268 - t267 * t269 + (t243 * t393 + t244 * t392 + t339 * t356) * pkin(3)) * MDP(21) + t402; (t389 * t412 + t288) * MDP(18) + (t329 * t389 + t289) * MDP(19) + (-t329 ^ 2 - t412 ^ 2) * MDP(20) + (t266 * t412 - t267 * t329 + t316) * MDP(21) + (-t403 + t454) * MDP(27) + (t245 + t453) * MDP(28); (t462 * t415 + t409) * MDP(27) + ((-t253 * t388 - t238) * t394 + (-t462 * t251 - t239) * t397 + t424) * MDP(28) + t472;];
tauc = t1;
