% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRP1
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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:07
% EndTime: 2019-12-05 18:46:11
% DurationCPUTime: 2.37s
% Computational Cost: add. (2893->270), mult. (7557->361), div. (0->0), fcn. (5391->6), ass. (0->152)
t336 = sin(qJ(3));
t339 = cos(qJ(3));
t340 = cos(qJ(2));
t393 = qJD(1) * t340;
t337 = sin(qJ(2));
t394 = qJD(1) * t337;
t301 = -t336 * t393 - t339 * t394;
t295 = t301 * pkin(8);
t420 = pkin(6) + pkin(7);
t319 = t420 * t340;
t315 = qJD(1) * t319;
t302 = t336 * t315;
t318 = t420 * t337;
t313 = qJD(1) * t318;
t413 = qJD(2) * pkin(2);
t308 = -t313 + t413;
t373 = t339 * t308 - t302;
t255 = t295 + t373;
t332 = qJD(2) + qJD(3);
t247 = pkin(3) * t332 + t255;
t306 = t339 * t315;
t361 = -t308 * t336 - t306;
t384 = t339 * t393;
t385 = t336 * t394;
t299 = -t384 + t385;
t416 = pkin(8) * t299;
t256 = -t361 - t416;
t338 = cos(qJ(4));
t251 = t338 * t256;
t335 = sin(qJ(4));
t364 = -t335 * t247 - t251;
t312 = t336 * t340 + t337 * t339;
t352 = t312 * qJD(2);
t345 = -t312 * qJD(3) - t352;
t344 = t345 * qJD(1);
t386 = qJD(2) * t420;
t366 = qJD(1) * t386;
t310 = t340 * t366;
t392 = qJD(3) * t336;
t371 = -t336 * t310 - t315 * t392;
t309 = t337 * t366;
t428 = (qJD(3) * t308 - t309) * t339;
t234 = pkin(8) * t344 + t371 + t428;
t388 = qJD(2) * qJD(1);
t383 = t340 * t388;
t278 = qJD(3) * t384 - t332 * t385 + t339 * t383;
t372 = t336 * t309 - t339 * t310;
t348 = t361 * qJD(3) + t372;
t235 = -pkin(8) * t278 + t348;
t379 = -t335 * t234 + t338 * t235;
t350 = t364 * qJD(4) + t379;
t328 = -t340 * pkin(2) - pkin(1);
t317 = t328 * qJD(1);
t285 = pkin(3) * t299 + t317;
t362 = t299 * t335 + t338 * t301;
t411 = t285 * t362;
t434 = t350 + t411;
t391 = qJD(4) * t335;
t378 = t335 * t235 - t256 * t391;
t351 = -(qJD(4) * t247 + t234) * t338 - t378;
t272 = -t338 * t299 + t301 * t335;
t410 = t285 * t272;
t433 = t351 - t410;
t268 = t272 ^ 2;
t331 = qJD(4) + t332;
t349 = t362 * qJD(4) - t278 * t335 + t338 * t344;
t390 = qJD(4) * t338;
t353 = t338 * t278 - t299 * t390 + t301 * t391 + t335 * t344;
t421 = t362 ^ 2;
t432 = t272 * t362 * MDP(18) + (-t272 * t331 + t353) * MDP(20) + (-t331 * t362 + t349) * MDP(21) + (-t268 + t421) * MDP(19);
t412 = qJ(5) * t272;
t262 = t362 * qJ(5);
t430 = -0.2e1 * t388;
t429 = t337 * MDP(4);
t427 = (t337 ^ 2 - t340 ^ 2) * MDP(5);
t329 = pkin(2) * t394;
t260 = -pkin(3) * t344 + qJD(2) * t329;
t426 = -pkin(4) * t349 + t260;
t370 = t313 * t336 - t306;
t258 = t370 + t416;
t396 = -t339 * t313 - t302;
t259 = t295 + t396;
t327 = pkin(2) * t339 + pkin(3);
t407 = t335 * t336;
t425 = -t327 * t390 - (-t336 * t391 + (t338 * t339 - t407) * qJD(3)) * pkin(2) + t335 * t258 + t338 * t259;
t405 = t336 * t338;
t424 = -t327 * t391 + (-t336 * t390 + (-t335 * t339 - t405) * qJD(3)) * pkin(2) - t338 * t258 + t259 * t335;
t423 = -pkin(8) * t312 - t336 * t319;
t314 = t337 * t386;
t316 = t340 * t386;
t360 = t336 * t318 - t339 * t319;
t422 = t360 * qJD(3) + t314 * t336 - t339 * t316;
t419 = pkin(2) * t337;
t418 = pkin(3) * t301;
t417 = pkin(4) * t362;
t409 = t317 * t301;
t408 = t318 * t339;
t249 = t335 * t256;
t341 = qJD(2) ^ 2;
t404 = t337 * t341;
t403 = t340 * t341;
t342 = qJD(1) ^ 2;
t402 = t340 * t342;
t377 = t338 * t247 - t249;
t223 = t377 + t262;
t221 = pkin(4) * t331 + t223;
t401 = t221 - t223;
t400 = t338 * t255 - t249;
t398 = -t262 - t425;
t397 = t412 + t424;
t330 = t337 * t413;
t387 = -qJD(3) * t408 - t339 * t314 - t336 * t316;
t382 = -pkin(2) * t332 - t308;
t381 = -pkin(3) * t331 - t247;
t286 = t329 - t418;
t380 = pkin(1) * t430;
t376 = -t255 * t335 - t251;
t224 = -t364 + t412;
t365 = t221 * t272 - t224 * t362;
t263 = -t408 + t423;
t359 = t336 * t337 - t339 * t340;
t264 = -t359 * pkin(8) - t360;
t363 = -t263 * t335 - t264 * t338;
t358 = t317 * t299 - t371;
t357 = t317 * t312;
t356 = t328 * t312;
t355 = t338 * t359;
t241 = -pkin(8) * t352 + qJD(3) * t423 + t387;
t284 = t332 * t359;
t242 = pkin(8) * t284 + t422;
t354 = t338 * t241 + t335 * t242 + t263 * t390 - t264 * t391;
t283 = t338 * t312 - t335 * t359;
t288 = t359 * pkin(3) + t328;
t347 = t363 * qJD(4) - t241 * t335 + t338 * t242;
t346 = -t301 * t299 * MDP(11) + (t299 * t332 + t278) * MDP(13) + (-t301 * t332 + t344) * MDP(14) + (-t299 ^ 2 + t301 ^ 2) * MDP(12) + t432;
t275 = -t345 * pkin(3) + t330;
t326 = pkin(3) * t338 + pkin(4);
t296 = pkin(2) * t405 + t327 * t335;
t293 = -pkin(2) * t407 + t327 * t338 + pkin(4);
t282 = t312 * t335 + t355;
t243 = -pkin(4) * t272 + qJD(5) + t285;
t237 = t283 * qJD(4) - t335 * t284 - t338 * t345;
t236 = qJD(4) * t355 + t338 * t284 + t312 * t391 - t335 * t345;
t230 = -qJ(5) * t282 - t363;
t229 = -qJ(5) * t283 + t263 * t338 - t264 * t335;
t226 = t262 + t400;
t225 = t376 - t412;
t218 = qJ(5) * t236 - qJD(5) * t283 + t347;
t217 = -qJ(5) * t237 - qJD(5) * t282 + t354;
t216 = -qJ(5) * t353 + qJD(5) * t362 + t350;
t215 = qJ(5) * t349 + qJD(5) * t272 - t351;
t1 = [0.2e1 * t383 * t429 + t427 * t430 + MDP(6) * t403 - MDP(7) * t404 + (-pkin(6) * t403 + t337 * t380) * MDP(9) + (pkin(6) * t404 + t340 * t380) * MDP(10) + (t278 * t312 + t284 * t301) * MDP(11) + (-t278 * t359 + t284 * t299 - t301 * t345 + t312 * t344) * MDP(12) + (t357 * qJD(3) + (t299 * t419 + t357) * qJD(2) + (qJD(3) * t356 + (t359 * t419 + t356) * qJD(2)) * qJD(1)) * MDP(16) + (t328 * t278 - t317 * t284 + (qJD(1) * t312 - t301) * t330) * MDP(17) + (t236 * t362 + t283 * t353) * MDP(18) + (-t236 * t272 + t237 * t362 - t282 * t353 + t283 * t349) * MDP(19) + (t285 * t237 + t260 * t282 - t272 * t275 - t288 * t349) * MDP(23) + (-t285 * t236 + t260 * t283 - t275 * t362 + t288 * t353) * MDP(24) + (-t215 * t282 - t216 * t283 + t217 * t272 + t218 * t362 + t221 * t236 - t224 * t237 - t229 * t353 + t230 * t349) * MDP(25) + (t215 * t230 + t224 * t217 + t216 * t229 + t221 * t218 + t426 * (t282 * pkin(4) + t288) + t243 * (t237 * pkin(4) + t275)) * MDP(26) + (-t284 * MDP(13) + t345 * MDP(14) + t422 * MDP(16) + (t319 * t392 - t387) * MDP(17)) * t332 + (-t236 * MDP(20) - t237 * MDP(21) + t347 * MDP(23) - t354 * MDP(24)) * t331; t346 + (t272 * t286 + t331 * t424 + t434) * MDP(23) + (t286 * t362 + t331 * t425 + t433) * MDP(24) + (t215 * t296 + t216 * t293 - t243 * (t286 - t417) + t398 * t224 + t397 * t221) * MDP(26) - t402 * t429 + (t301 * t329 + t396 * t332 + (t382 * qJD(3) + t309) * t339 + t358) * MDP(17) + (-t299 * t329 + t409 - t370 * t332 + (t382 * t336 - t306) * qJD(3) + t372) * MDP(16) + (t272 * t398 - t293 * t353 + t296 * t349 + t362 * t397 + t365) * MDP(25) + t342 * t427 + (t342 * t337 * MDP(9) + MDP(10) * t402) * pkin(1); (-t361 * t332 + t348 + t409) * MDP(16) + (t373 * t332 + t358 - t428) * MDP(17) + (-t272 * t418 + t411 - t376 * t331 + (t381 * t335 - t251) * qJD(4) + t379) * MDP(23) + (-t362 * t418 - t410 + t400 * t331 + (t381 * qJD(4) - t234) * t338 - t378) * MDP(24) + (-t225 * t362 - t226 * t272 - t353 * t326 + (t349 * t335 + (t272 * t338 - t335 * t362) * qJD(4)) * pkin(3) + t365) * MDP(25) + (t243 * t417 + t216 * t326 - t221 * t225 - t224 * t226 + (t215 * t335 + t243 * t301 + (-t221 * t335 + t224 * t338) * qJD(4)) * pkin(3)) * MDP(26) + t346; (-t364 * t331 + t434) * MDP(23) + (t377 * t331 + t433) * MDP(24) + (-pkin(4) * t353 + t272 * t401) * MDP(25) + (t401 * t224 + (t243 * t362 + t216) * pkin(4)) * MDP(26) + t432; (-t268 - t421) * MDP(25) + (-t221 * t362 - t224 * t272 + t426) * MDP(26);];
tauc = t1;
