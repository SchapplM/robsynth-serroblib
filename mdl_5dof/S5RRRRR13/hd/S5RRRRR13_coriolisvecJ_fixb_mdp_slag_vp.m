% Calculate Coriolis joint torque vector for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR13_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR13_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:17
% EndTime: 2024-09-27 17:32:20
% DurationCPUTime: 0.94s
% Computational Cost: add. (2939->252), mult. (5463->352), div. (0->0), fcn. (3466->10), ass. (0->153)
t327 = sin(qJ(5));
t331 = cos(qJ(5));
t322 = qJD(1) + qJD(2);
t318 = qJD(3) + t322;
t325 = sin(pkin(5));
t332 = cos(qJ(4));
t420 = t325 * t332;
t393 = t318 * t420;
t328 = sin(qJ(4));
t424 = t318 * t325;
t394 = t328 * t424;
t441 = -t327 * t394 + t331 * t393;
t329 = sin(qJ(3));
t334 = cos(qJ(2));
t330 = sin(qJ(2));
t333 = cos(qJ(3));
t416 = t330 * t333;
t350 = -t329 * t334 - t416;
t429 = pkin(1) * qJD(1);
t295 = t350 * t429;
t410 = qJD(3) * t329;
t401 = pkin(2) * t410;
t364 = t295 + t401;
t321 = t325 ^ 2;
t421 = t321 * t332;
t440 = MDP(10) * t328 * t421 - (t328 ^ 2 - t332 ^ 2) * MDP(11) * t321;
t434 = qJD(4) + qJD(5);
t402 = t334 * t429;
t306 = t322 * pkin(2) + t402;
t403 = t330 * t429;
t286 = -t306 * t329 - t333 * t403;
t407 = qJD(4) * t328;
t400 = pkin(4) * t407;
t439 = t325 * (t286 + t400);
t289 = (t327 * t328 - t331 * t332) * t325;
t352 = t327 * t332 + t328 * t331;
t290 = t352 * t325;
t371 = t329 * t403;
t285 = t333 * t306 - t371;
t326 = cos(pkin(5));
t406 = qJD(4) * t332;
t387 = t326 * t406;
t419 = t326 * t328;
t437 = -pkin(3) * t387 + t285 * t332 + t286 * t419;
t417 = t329 * t330;
t349 = t333 * t334 - t417;
t296 = t349 * t429;
t418 = t326 * t332;
t436 = t295 * t418 - t296 * t328 - (-t328 * t333 - t329 * t418) * qJD(3) * pkin(2);
t315 = pkin(2) * t333 + pkin(3);
t409 = qJD(3) * t333;
t435 = t295 * t419 - t315 * t387 + (-pkin(2) * t409 + t296) * t332;
t273 = pkin(9) * t424 - t286;
t365 = pkin(10) * t424 + t273;
t432 = pkin(3) * t318;
t275 = t285 + t432;
t398 = t275 * t419;
t239 = t365 * t332 + t398;
t369 = qJD(2) * t402;
t413 = (qJD(2) + qJD(3)) * t371;
t257 = (qJD(3) * t306 + t369) * t333 - t413;
t433 = t330 * MDP(5) + t334 * MDP(6);
t431 = pkin(4) * t332;
t430 = pkin(10) * t325;
t336 = (t350 * qJD(2) - t330 * t409) * pkin(1);
t335 = qJD(1) * t336;
t258 = -t306 * t410 + t335;
t386 = t275 * t406;
t391 = t332 * t257 + t258 * t419 + t326 * t386;
t343 = t273 * t407 - t391;
t428 = t343 * t326;
t427 = t239 * t331;
t256 = (-t318 * t431 - t275) * t325;
t280 = -t327 * t393 - t331 * t394;
t426 = t256 * t280;
t425 = t273 * t332;
t423 = t318 * t326;
t422 = t321 * t328;
t253 = t258 * t418;
t380 = -t328 * t257 + t253;
t415 = ((-t398 - t425) * qJD(4) + t380) * t326 + t258 * t421;
t307 = qJD(4) + t423;
t405 = qJD(4) - t307;
t404 = pkin(3) * t419;
t399 = t325 * (-pkin(9) - pkin(10));
t316 = pkin(1) * t334 + pkin(2);
t292 = -pkin(1) * t417 + t316 * t333 + pkin(3);
t397 = t292 * t419;
t396 = t315 * t419;
t395 = t318 * t421;
t242 = (t318 * t400 - t258) * t325;
t264 = t434 * t290;
t270 = t275 * t418;
t348 = t365 * t328;
t238 = t270 - t348;
t233 = pkin(4) * t307 + t238;
t353 = -t327 * t233 - t427;
t222 = -qJD(4) * t348 + t391;
t223 = -qJD(4) * t239 + t380;
t381 = -t327 * t222 + t331 * t223;
t338 = t353 * qJD(5) + t381;
t392 = t242 * t289 + t256 * t264 + t338 * t326;
t271 = t316 * t409 + (t349 * qJD(2) - t330 * t410) * pkin(1);
t272 = -t316 * t410 + t336;
t390 = t332 * t271 + t272 * t419 + t292 * t387;
t389 = MDP(12) * t420;
t388 = t325 * t407;
t385 = -t275 - t432;
t305 = qJD(5) + t307;
t384 = -pkin(4) * t305 - t233;
t319 = t325 * pkin(9);
t288 = pkin(1) * t416 + t316 * t329 + t319;
t383 = -t288 - t430;
t310 = pkin(2) * t329 + t319;
t382 = -t310 - t430;
t379 = -t272 * t318 - t258;
t378 = t286 * t318 - t258;
t377 = -t271 * t328 + t272 * t418;
t376 = -t292 * t318 - t275;
t311 = pkin(4) * t388;
t375 = t364 * t325 + t311;
t373 = pkin(4) * t394;
t372 = t326 * t401;
t370 = t328 * t399;
t234 = qJD(5) * t327 * t239;
t263 = t434 * t289;
t360 = -(t327 * t223 - t234 + (qJD(5) * t233 + t222) * t331) * t326 - t256 * t263 + t242 * t290;
t246 = -t285 * t328 + t286 * t418;
t314 = pkin(10) * t420;
t345 = -pkin(9) * t420 - t404;
t359 = qJD(5) * (t314 - t345) + t246 - (t332 * t399 - t404) * qJD(4);
t320 = t326 * pkin(4);
t358 = -qJD(5) * (pkin(3) * t418 + t320 + t370) - qJD(4) * t370 + t437;
t357 = -qJD(5) * (t315 * t418 + t382 * t328 + t320) - (t382 * qJD(4) - t372) * t328 + t435;
t346 = -t310 * t332 - t396;
t356 = qJD(5) * (t314 - t346) - (t382 * t332 - t396) * qJD(4) + t436;
t355 = t383 * t328;
t354 = (-pkin(2) * t318 - t306) * qJD(3);
t347 = -t288 * t332 - t397;
t243 = t441 * t434;
t344 = t280 * t441 * MDP(17) + (-t305 * t441 + t243) * MDP(19) + (-t352 * t424 * t434 - t280 * t305) * MDP(20) + (t280 ^ 2 - t441 ^ 2) * MDP(18);
t244 = t318 * t264;
t341 = (-t243 * t289 - t244 * t290 - t263 * t441 + t264 * t280) * MDP(18) + (t243 * t290 + t263 * t280) * MDP(17) + (t243 * t326 - t263 * t305) * MDP(19) + (-t244 * t326 - t264 * t305) * MDP(20) + 0.2e1 * t440 * qJD(4) * t318 + (-MDP(13) * t388 + qJD(4) * t389) * (t307 + t423);
t339 = t234 - t256 * t441 + (-t239 * t305 - t223) * t327;
t304 = (-pkin(3) - t431) * t325;
t293 = (-t315 - t431) * t325;
t277 = (-t292 - t431) * t325;
t255 = -t272 * t325 + t311;
t250 = t314 - t347;
t245 = t292 * t418 + t320 + t355;
t228 = (t383 * t332 - t397) * qJD(4) + t377;
t227 = qJD(4) * t355 + t390;
t1 = [-t379 * MDP(8) + (-t271 * t318 - t257) * MDP(9) + ((-t227 * t327 + t228 * t331 + (-t245 * t327 - t250 * t331) * qJD(5)) * t305 - t255 * t441 + t277 * t244 + t392) * MDP(22) + (-(t227 * t331 + t228 * t327 + (t245 * t331 - t250 * t327) * qJD(5)) * t305 - t255 * t280 + t277 * t243 + t360) * MDP(23) + (-(-t288 * t407 + t390) * t307 + t428 + (t379 * t328 + t376 * t406) * t321) * MDP(16) + (t377 * t307 + t272 * t395 + (t347 * t307 + t376 * t422) * qJD(4) + t415) * MDP(15) + t341 + t433 * pkin(1) * qJD(2) * (-qJD(1) - t322); (t428 + ((qJD(4) * t310 + t372) * t328 + t435) * t307 + (-t386 - t258 * t328 + (-t315 * t406 + t364 * t328) * t318) * t321) * MDP(16) + (t296 * t318 + (t354 - t369) * t333 + t413) * MDP(9) + (t293 * t243 + (t356 * t327 + t357 * t331) * t305 - t375 * t280 + t360) * MDP(23) + (t293 * t244 + (t357 * t327 - t356 * t331) * t305 - t375 * t441 + t392) * MDP(22) + (-t295 * t318 + t329 * t354 + t335) * MDP(8) + ((t346 * qJD(4) - t436) * t307 + (-t275 * t407 + (-t315 * t407 - t364 * t332) * t318) * t321 + t415) * MDP(15) + t341 + t433 * (-qJD(2) + t322) * t429; -t378 * MDP(8) + (t285 * t318 - t257) * MDP(9) + (-t286 * t395 - t246 * t307 + (t345 * t307 + t385 * t422) * qJD(4) + t415) * MDP(15) + (t428 + (pkin(9) * t388 + t437) * t307 + (t378 * t328 + t385 * t406) * t321) * MDP(16) + (t304 * t244 + (t358 * t327 - t359 * t331) * t305 - t441 * t439 + t392) * MDP(22) + (t304 * t243 + (t359 * t327 + t358 * t331) * t305 - t280 * t439 + t360) * MDP(23) + t341; t405 * t318 * t389 - t405 * MDP(13) * t394 + (t253 - t405 * t425 + (-t257 + (t318 * t321 - t405 * t326) * t275) * t328) * MDP(15) + ((-t273 * t328 + t270) * t307 + t275 * t395 + t343) * MDP(16) + (-(-t238 * t327 - t427) * t305 + t441 * t373 + t426 + (t384 * t327 - t427) * qJD(5) + t381) * MDP(22) + (t280 * t373 + (t384 * qJD(5) + t238 * t305 - t222) * t331 + t339) * MDP(23) + t344 - t440 * t318 ^ 2; (-t353 * t305 + t338 + t426) * MDP(22) + ((-t222 + (-qJD(5) + t305) * t233) * t331 + t339) * MDP(23) + t344;];
tauc = t1;
