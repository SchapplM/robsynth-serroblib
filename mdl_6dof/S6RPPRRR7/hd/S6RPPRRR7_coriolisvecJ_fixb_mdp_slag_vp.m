% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:55
% EndTime: 2019-03-09 02:34:03
% DurationCPUTime: 3.83s
% Computational Cost: add. (3821->325), mult. (8776->432), div. (0->0), fcn. (6816->8), ass. (0->148)
t353 = cos(qJ(6));
t399 = qJD(6) * t353;
t347 = sin(pkin(10));
t348 = cos(pkin(10));
t355 = cos(qJ(4));
t405 = qJD(1) * t355;
t352 = sin(qJ(4));
t406 = qJD(1) * t352;
t317 = -t347 * t405 - t348 * t406;
t389 = t347 * t406;
t392 = t348 * t405;
t318 = -t389 + t392;
t351 = sin(qJ(5));
t354 = cos(qJ(5));
t290 = -t354 * t317 + t318 * t351;
t451 = t290 * t353;
t456 = t399 + t451;
t322 = -t347 * t352 + t348 * t355;
t368 = t317 * t351 + t354 * t318;
t350 = sin(qJ(6));
t400 = qJD(6) * t350;
t321 = t347 * t355 + t348 * t352;
t319 = t321 * qJD(4);
t311 = qJD(1) * t319;
t327 = qJD(4) * t389;
t403 = qJD(4) * t355;
t391 = t348 * t403;
t312 = qJD(1) * t391 - t327;
t401 = qJD(5) * t354;
t402 = qJD(5) * t351;
t264 = -t354 * t311 - t351 * t312 + t317 * t401 - t318 * t402;
t344 = qJD(4) + qJD(5);
t410 = t353 * t264 + t344 * t399;
t247 = -t368 * t400 + t410;
t244 = t247 * t353;
t284 = t344 * t350 + t353 * t368;
t430 = t264 * t350;
t248 = qJD(6) * t284 + t430;
t420 = t368 * t350;
t282 = -t353 * t344 + t420;
t455 = -t350 * t248 - t456 * t282 + t244;
t243 = t247 * t350;
t265 = qJD(5) * t368 - t311 * t351 + t354 * t312;
t260 = t350 * t265;
t448 = -qJD(6) - t290;
t411 = -t399 * t448 + t260;
t414 = t353 * t448;
t422 = t290 * t344;
t424 = t368 * t344;
t426 = t284 * t368;
t454 = (-t265 + t424) * MDP(21) - t290 ^ 2 * MDP(19) + (-t290 * t414 + t411 - t426) * MDP(27) + (t290 * MDP(18) + MDP(19) * t368 + MDP(29) * t448) * t368 + (t264 + t422) * MDP(20) + (t456 * t284 + t243) * MDP(25);
t349 = -pkin(1) - qJ(3);
t440 = qJD(1) * t349;
t328 = qJD(2) + t440;
t386 = -pkin(7) * qJD(1) + t328;
t313 = t386 * t347;
t314 = t386 * t348;
t442 = -t313 * t352 + t355 * t314;
t278 = -pkin(8) * t318 + t442;
t277 = qJD(4) * pkin(4) + t278;
t369 = -t313 * t355 - t314 * t352;
t279 = pkin(8) * t317 - t369;
t429 = t279 * t351;
t251 = t277 * t354 - t429;
t249 = -pkin(5) * t344 - t251;
t453 = t249 * t290;
t346 = qJD(1) * qJ(2);
t340 = qJD(3) + t346;
t341 = t347 * pkin(3);
t325 = qJD(1) * t341 + t340;
t298 = -pkin(4) * t317 + t325;
t452 = t290 * t298;
t379 = t448 * t350;
t450 = t322 * qJD(3);
t409 = t347 ^ 2 + t348 ^ 2;
t449 = qJD(3) * t409;
t269 = pkin(5) * t368 + pkin(9) * t290;
t447 = qJ(2) * MDP(6) + t347 * MDP(7) + t348 * MDP(8) + MDP(5);
t427 = t282 * t368;
t262 = t353 * t265;
t441 = -t400 * t448 - t262;
t428 = t279 * t354;
t252 = t277 * t351 + t428;
t362 = t321 * qJD(3);
t266 = -pkin(8) * t312 - qJD(1) * t362 + qJD(4) * t442;
t364 = t450 * qJD(1);
t267 = pkin(8) * t311 + qJD(4) * t369 - t364;
t384 = t266 * t351 - t354 * t267;
t234 = qJD(5) * t252 + t384;
t250 = pkin(9) * t344 + t252;
t255 = pkin(5) * t290 - pkin(9) * t368 + t298;
t371 = t250 * t350 - t255 * t353;
t439 = -t234 * t353 + t249 * t400 + t368 * t371;
t236 = t250 * t353 + t255 * t350;
t438 = t234 * t350 + t236 * t368 + t249 * t399;
t437 = -t368 * t298 - t384;
t383 = t267 * t351 - t279 * t402;
t233 = (qJD(5) * t277 + t266) * t354 + t383;
t404 = qJD(4) * t352;
t320 = -t347 * t404 + t391;
t367 = t354 * t321 + t322 * t351;
t272 = qJD(5) * t367 + t354 * t319 + t320 * t351;
t434 = -pkin(7) + t349;
t323 = t434 * t347;
t324 = t434 * t348;
t359 = -t323 * t404 + t324 * t403 - t362;
t275 = -pkin(8) * t320 + t359;
t366 = -t323 * t355 - t324 * t352;
t358 = qJD(4) * t366 - t450;
t276 = pkin(8) * t319 + t358;
t285 = -pkin(8) * t322 - t323 * t352 + t324 * t355;
t286 = -pkin(8) * t321 - t366;
t370 = t285 * t354 - t286 * t351;
t237 = qJD(5) * t370 + t275 * t354 + t276 * t351;
t258 = t285 * t351 + t286 * t354;
t294 = t321 * t351 - t322 * t354;
t336 = qJ(2) + t341;
t304 = pkin(4) * t321 + t336;
t259 = pkin(5) * t367 + pkin(9) * t294 + t304;
t436 = -t234 * t294 - t249 * t272 - t258 * t265 + (qJD(6) * t259 + t237) * t448 - t367 * (qJD(6) * t255 + t233);
t435 = pkin(4) * t318;
t432 = t249 * t294;
t431 = t259 * t265;
t425 = t284 * t350;
t407 = qJD(1) * t321;
t395 = MDP(13) * qJD(4);
t345 = qJD(1) * qJD(2);
t387 = -pkin(4) * t344 - t277;
t299 = pkin(4) * t312 + t345;
t306 = pkin(4) * t320 + qJD(2);
t380 = qJD(1) * t409;
t338 = pkin(4) * t351 + pkin(9);
t375 = qJD(6) * t338 + t269 + t435;
t253 = t278 * t351 + t428;
t374 = pkin(4) * t402 - t253;
t254 = t278 * t354 - t429;
t373 = -pkin(4) * t401 + t254;
t372 = -t265 * t338 + t453;
t365 = t290 * t379 - t441;
t363 = -t272 * t353 + t294 * t400;
t357 = -qJD(5) * t294 - t319 * t351 + t354 * t320;
t356 = qJD(1) ^ 2;
t339 = -pkin(4) * t354 - pkin(5);
t241 = pkin(5) * t357 + pkin(9) * t272 + t306;
t240 = pkin(5) * t265 - pkin(9) * t264 + t299;
t239 = t353 * t240;
t238 = qJD(5) * t258 + t275 * t351 - t276 * t354;
t1 = [((t340 + t346) * qJD(2) + (-t328 - t440) * t449) * MDP(10) + (-t236 * t357 + t238 * t284 - t370 * t247 + ((-qJD(6) * t258 + t241) * t448 - t431 - (-qJD(6) * t250 + t240) * t367 + qJD(6) * t432) * t350 + t436 * t353) * MDP(31) + (-t371 * t357 + t238 * t282 + t239 * t367 - t370 * t248 + (-t241 * t448 + t431 + (-t250 * t367 + t258 * t448 - t432) * qJD(6)) * t353 + t436 * t350) * MDP(30) - t319 * t395 + (t294 * t260 - t248 * t367 - t357 * t282 - (t272 * t350 + t294 * t399) * t448) * MDP(28) + (-t264 * t294 - t272 * t368) * MDP(18) + (t264 * t304 - t272 * t298 - t294 * t299 + t306 * t368) * MDP(24) + (-t264 * t367 + t265 * t294 + t272 * t290 - t357 * t368) * MDP(19) + (-t244 * t294 + t284 * t363) * MDP(25) + (-(-t282 * t353 - t425) * t272 - (-t243 - t248 * t353 + (t282 * t350 - t284 * t353) * qJD(6)) * t294) * MDP(26) + (-t336 * t311 - t325 * t319 - t359 * qJD(4) + (qJD(1) * t322 + t318) * qJD(2)) * MDP(17) + (-t311 * t322 - t318 * t319) * MDP(11) + (t311 * t321 - t312 * t322 - t317 * t319 - t318 * t320) * MDP(12) + (t265 * t367 - t357 * t448) * MDP(29) + (t247 * t367 - t262 * t294 + t284 * t357 - t363 * t448) * MDP(27) + (t265 * t304 + t290 * t306 + t298 * t357 + t299 * t367) * MDP(23) - t320 * qJD(4) * MDP(14) + 0.2e1 * qJD(3) * MDP(9) * t380 + (t336 * t312 + t325 * t320 + t358 * qJD(4) + (-t317 + t407) * qJD(2)) * MDP(16) + (-MDP(20) * t272 - MDP(21) * t357 - MDP(23) * t238 - MDP(24) * t237) * t344 + 0.2e1 * t447 * t345; (t248 * t294 - t260 * t367 + t272 * t282) * MDP(30) + (t247 * t294 - t262 * t367 + t272 * t284) * MDP(31) + (-MDP(23) * t272 - MDP(24) * t357) * t344 + (-MDP(16) * t319 - MDP(17) * t320) * qJD(4) + ((-t340 - t449) * MDP(10) + t317 * MDP(16) - t318 * MDP(17) - t290 * MDP(23) - t368 * MDP(24)) * qJD(1) - t447 * t356 - ((-qJD(1) * t353 - t350 * t357 - t367 * t399) * MDP(30) + (qJD(1) * t350 - t353 * t357 + t367 * t400) * MDP(31)) * t448; (t328 * t380 + t345) * MDP(10) - t327 * MDP(16) + (t265 + t424) * MDP(23) + (t264 - t422) * MDP(24) + (t365 - t427) * MDP(30) + (-t414 * t448 - t260 - t426) * MDP(31) - t409 * MDP(9) * t356 + ((t318 + t392) * MDP(16) + 0.2e1 * t317 * MDP(17)) * qJD(4); (-t290 * t435 + t253 * t344 + (t351 * t387 - t428) * qJD(5) + t437) * MDP(23) + (t339 * t248 + t372 * t350 + t374 * t282 - (t350 * t373 - t353 * t375) * t448 + t439) * MDP(30) + (t339 * t247 + t372 * t353 + t374 * t284 - (t350 * t375 + t353 * t373) * t448 + t438) * MDP(31) + (-t368 * t435 + t254 * t344 + t452 + (qJD(5) * t387 - t266) * t354 - t383) * MDP(24) + (-t325 * t318 - t364) * MDP(16) + (t425 * t448 + t455) * MDP(26) - t318 * t317 * MDP(11) + (t327 + (t318 - t392) * qJD(4)) * MDP(14) + (-t317 - t407) * t395 + (qJD(3) * t407 - t325 * t317) * MDP(17) + (-t317 ^ 2 + t318 ^ 2) * MDP(12) + (t365 + t427) * MDP(28) + t454; ((-qJD(5) + t344) * t252 + t437) * MDP(23) + (t251 * t344 - t233 + t452) * MDP(24) + (t284 * t379 + t455) * MDP(26) + (-t379 * t448 + t262 + t427) * MDP(28) + (-pkin(5) * t248 + (-t251 * t350 + t269 * t353) * t448 - t252 * t282 + t350 * t453 - t411 * pkin(9) + t439) * MDP(30) + (-pkin(5) * t247 - (t251 * t353 + t269 * t350) * t448 - t252 * t284 + t249 * t451 + t441 * pkin(9) + t438) * MDP(31) + t454; t284 * t282 * MDP(25) + (-t282 ^ 2 + t284 ^ 2) * MDP(26) + (-t282 * t448 + t410) * MDP(27) + (-t284 * t448 - t430) * MDP(28) + t265 * MDP(29) + (-t233 * t350 - t236 * t448 - t249 * t284 + t239) * MDP(30) + (-t233 * t353 - t240 * t350 + t249 * t282 + t371 * t448) * MDP(31) + (-MDP(27) * t420 - MDP(28) * t284 - MDP(30) * t236 + MDP(31) * t371) * qJD(6);];
tauc  = t1;
