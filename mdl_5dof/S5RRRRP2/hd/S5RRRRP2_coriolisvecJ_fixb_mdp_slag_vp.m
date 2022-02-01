% Calculate Coriolis joint torque vector for
% S5RRRRP2
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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:31
% EndTime: 2022-01-20 11:49:35
% DurationCPUTime: 2.05s
% Computational Cost: add. (2499->258), mult. (4214->330), div. (0->0), fcn. (2668->6), ass. (0->150)
t338 = cos(qJ(3));
t331 = qJD(1) + qJD(2);
t336 = sin(qJ(2));
t421 = pkin(1) * qJD(1);
t382 = t336 * t421;
t434 = pkin(7) + pkin(8);
t373 = t434 * t331 + t382;
t276 = t373 * t338;
t334 = sin(qJ(4));
t268 = t334 * t276;
t335 = sin(qJ(3));
t275 = t373 * t335;
t271 = qJD(3) * pkin(3) - t275;
t337 = cos(qJ(4));
t368 = t337 * t271 - t268;
t409 = t335 * t337;
t302 = t334 * t338 + t409;
t285 = t302 * t331;
t414 = t285 * qJ(5);
t229 = -t414 + t368;
t330 = qJD(3) + qJD(4);
t259 = t330 * t302;
t246 = t259 * t331;
t359 = qJD(3) * t373;
t339 = cos(qJ(2));
t420 = pkin(1) * qJD(2);
t379 = qJD(1) * t420;
t363 = t339 * t379;
t255 = -t335 * t359 + t338 * t363;
t256 = -t335 * t363 - t338 * t359;
t390 = qJD(4) * t334;
t349 = -(qJD(4) * t271 + t255) * t337 - t334 * t256 + t276 * t390;
t433 = -qJ(5) * t246 - t349;
t389 = qJD(4) * t337;
t392 = qJD(3) * t338;
t432 = t337 * t392 + t338 * t389;
t364 = t330 * MDP(21);
t431 = MDP(7) * t335;
t430 = MDP(8) * (t335 ^ 2 - t338 ^ 2);
t378 = qJD(3) * t434;
t303 = t335 * t378;
t304 = t338 * t378;
t313 = t434 * t335;
t328 = t338 * pkin(8);
t314 = pkin(7) * t338 + t328;
t356 = t313 * t334 - t314 * t337;
t381 = t339 * t421;
t429 = qJD(4) * t356 + t302 * t381 + t334 * t303 - t337 * t304;
t410 = t334 * t335;
t301 = -t337 * t338 + t410;
t428 = -t301 * t381 + t337 * t303 + t334 * t304 + t313 * t389 + t314 * t390;
t427 = t335 * MDP(12) + t338 * MDP(13);
t426 = t285 ^ 2;
t424 = pkin(1) * t339;
t423 = pkin(3) * t330;
t322 = pkin(1) * t336 + pkin(7);
t422 = -pkin(8) - t322;
t360 = t330 * t410;
t395 = t432 * t331;
t245 = t331 * t360 - t395;
t419 = qJ(5) * t245;
t283 = t301 * t331;
t417 = qJ(5) * t283;
t416 = qJ(5) * t302;
t415 = t283 * t330;
t325 = -pkin(3) * t338 - pkin(2);
t287 = t325 * t331 - t381;
t412 = t287 * t285;
t411 = t331 * t335;
t340 = qJD(3) ^ 2;
t408 = t335 * t340;
t270 = t337 * t276;
t406 = t338 * t340;
t399 = -t259 * qJ(5) - t301 * qJD(5);
t405 = t399 - t428;
t258 = t360 - t432;
t355 = qJ(5) * t258 - qJD(5) * t302;
t404 = t355 + t429;
t228 = pkin(4) * t330 + t229;
t403 = t228 - t229;
t319 = t336 * t379;
t393 = qJD(3) * t335;
t377 = t331 * t393;
t289 = pkin(3) * t377 + t319;
t240 = pkin(4) * t246 + t289;
t371 = pkin(4) * t283 + qJD(5);
t248 = t287 + t371;
t402 = t240 * t301 + t248 * t259;
t401 = t240 * t302 - t248 * t258;
t400 = t287 * t259 + t289 * t301;
t398 = -t287 * t258 + t289 * t302;
t397 = -t337 * t275 - t268;
t307 = -pkin(2) * t331 - t381;
t396 = t307 * t392 + t335 * t319;
t391 = qJD(3) * t339;
t387 = t338 * MDP(12);
t385 = -qJD(2) + t331;
t384 = qJD(5) + t248;
t383 = pkin(3) * t411;
t380 = t339 * t420;
t326 = pkin(3) * t393;
t376 = t331 * t392;
t252 = pkin(4) * t259 + t326;
t372 = qJD(3) * t422;
t370 = -t334 * t255 + t337 * t256;
t367 = t275 * t334 - t270;
t365 = t330 * t335;
t219 = -qJD(5) * t283 + t433;
t358 = -t271 * t334 - t270;
t347 = qJD(4) * t358 + t370;
t343 = t347 + t419;
t220 = -qJD(5) * t285 + t343;
t230 = -t358 - t417;
t362 = -t219 * t301 - t220 * t302 + t228 * t258 - t230 * t259;
t361 = t252 - t382;
t297 = t422 * t335;
t298 = t322 * t338 + t328;
t357 = -t297 * t334 - t298 * t337;
t282 = pkin(4) * t301 + t325;
t281 = t283 ^ 2;
t354 = t285 * t283 * MDP(14) + (-t331 * t334 * t365 + t395 + t415) * MDP(16) + (-t281 + t426) * MDP(15);
t272 = t335 * t372 + t338 * t380;
t273 = -t335 * t380 + t338 * t372;
t353 = t337 * t272 + t334 * t273 + t297 * t389 - t298 * t390;
t351 = t326 - t382;
t348 = -MDP(10) * t408 + (t245 * t301 - t246 * t302 + t258 * t283 - t259 * t285) * MDP(15) + (-t245 * t302 - t258 * t285) * MDP(14) - 0.2e1 * t331 * qJD(3) * t430 + 0.2e1 * t376 * t431 + MDP(9) * t406 + (-t258 * MDP(16) - t259 * MDP(17)) * t330;
t346 = qJD(4) * t357 - t272 * t334 + t337 * t273;
t344 = t287 * t283 + t349;
t342 = (-t270 + (-t271 - t423) * t334) * qJD(4) + t370;
t341 = t283 * t384 - t433;
t327 = t336 * t420;
t324 = -pkin(2) - t424;
t323 = pkin(3) * t337 + pkin(4);
t315 = t389 * t423;
t310 = t325 - t424;
t305 = t327 + t326;
t296 = t301 * qJ(5);
t290 = t307 * t393;
t274 = t282 - t424;
t262 = pkin(4) * t285 + t383;
t250 = -t296 - t356;
t249 = -t313 * t337 - t314 * t334 - t416;
t247 = t252 + t327;
t242 = -t296 - t357;
t241 = t297 * t337 - t298 * t334 - t416;
t234 = -t414 + t397;
t233 = t367 + t417;
t222 = t346 + t355;
t221 = t353 + t399;
t1 = [(-t310 * t245 + t305 * t285 - t330 * t353 + t398) * MDP(20) + (-t221 * t330 - t245 * t274 + t247 * t285 + t401) * MDP(22) + (-t221 * t283 - t222 * t285 + t241 * t245 - t242 * t246 + t362) * MDP(23) + (t310 * t246 + t305 * t283 + t330 * t346 + t400) * MDP(19) - t319 * MDP(5) + (-t322 * t406 + t324 * t377 + t290) * MDP(12) + (t322 * t408 + t324 * t376 + t396) * MDP(13) + (t222 * t330 + t246 * t274 + t247 * t283 + t402) * MDP(21) + (t219 * t242 + t220 * t241 + t221 * t230 + t222 * t228 + t240 * t274 + t247 * t248) * MDP(24) + t348 + (((-qJD(1) - t331) * MDP(6) - t427 * qJD(3)) * t339 + (-qJD(1) * t387 + (t335 * MDP(13) - MDP(5) - t387) * t331) * t336) * t420; (t219 * t250 + t220 * t249 + t228 * t404 + t230 * t405 + t240 * t282 + t248 * t361) * MDP(24) + (-pkin(2) * t377 - pkin(7) * t406 + t290 + (t336 * t338 * t385 + t335 * t391) * t421) * MDP(12) + (-pkin(2) * t376 + pkin(7) * t408 + (-t336 * t411 + t338 * t391) * t421 + t396) * MDP(13) + (t246 * t282 + t283 * t361 + t330 * t404 + t402) * MDP(21) + t385 * MDP(6) * t381 + (-t325 * t245 + t351 * t285 + t428 * t330 + t398) * MDP(20) + (-t245 * t282 + t285 * t361 - t330 * t405 + t401) * MDP(22) + (t245 * t249 - t246 * t250 - t283 * t405 - t285 * t404 + t362) * MDP(23) + (t325 * t246 + t351 * t283 + t429 * t330 + t400) * MDP(19) + (t331 * t382 - t319) * MDP(5) + t348; (-t283 * t383 - t330 * t367 + t342 - t412) * MDP(19) + (-t285 * t383 + t330 * t397 - t315 + t344) * MDP(20) + (-t233 * t330 - t262 * t283 - t285 * t384 + t342 + t419) * MDP(21) + (t234 * t330 - t262 * t285 - t315 + t341) * MDP(22) + (t245 * t323 + (t230 + t233) * t285 + (-t228 + t234) * t283 + (-t246 * t334 + (-t283 * t337 + t285 * t334) * qJD(4)) * pkin(3)) * MDP(23) + (t220 * t323 - t228 * t233 - t230 * t234 - t248 * t262 + (t219 * t334 + (-t228 * t334 + t230 * t337) * qJD(4)) * pkin(3)) * MDP(24) + t354 + t427 * (-t307 * t331 - t363) + (-t338 * t431 + t430) * t331 ^ 2; (-t330 * t358 + t347 - t412) * MDP(19) + (t330 * t368 + t344) * MDP(20) + (t230 * t330 + (-t248 - t371) * t285 + t343) * MDP(21) + (-pkin(4) * t426 + t229 * t330 + t341) * MDP(22) + (pkin(4) * t245 - t283 * t403) * MDP(23) + (t403 * t230 + (-t248 * t285 + t220) * pkin(4)) * MDP(24) + t354; t285 * t364 + (t395 - t415) * MDP(22) + (-t281 - t426) * MDP(23) + (t228 * t285 + t230 * t283 + t240) * MDP(24) + (t364 * t409 + (-MDP(22) * t365 + t338 * t364) * t334) * t331;];
tauc = t1;
