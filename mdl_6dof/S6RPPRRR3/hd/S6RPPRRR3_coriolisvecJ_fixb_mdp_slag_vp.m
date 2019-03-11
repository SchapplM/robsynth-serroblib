% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:52
% EndTime: 2019-03-09 02:23:58
% DurationCPUTime: 3.10s
% Computational Cost: add. (2023->346), mult. (4333->489), div. (0->0), fcn. (2817->8), ass. (0->162)
t335 = sin(qJ(5));
t338 = cos(qJ(5));
t386 = t338 * qJD(4);
t339 = cos(qJ(4));
t398 = qJD(1) * t339;
t301 = t335 * t398 - t386;
t337 = cos(qJ(6));
t377 = t338 * t398;
t396 = qJD(4) * t335;
t303 = t377 + t396;
t334 = sin(qJ(6));
t419 = t303 * t334;
t248 = t337 * t301 + t419;
t336 = sin(qJ(4));
t399 = qJD(1) * t336;
t321 = qJD(5) + t399;
t318 = qJD(6) + t321;
t445 = t248 * t318;
t304 = t334 * t335 - t337 * t338;
t436 = qJD(5) + qJD(6);
t444 = t304 * t436;
t354 = t301 * t334 - t337 * t303;
t443 = t318 * t354;
t305 = t334 * t338 + t335 * t337;
t342 = t436 * t305;
t395 = qJD(4) * t336;
t374 = t335 * t395;
t391 = qJD(5) * t339;
t344 = -t338 * t391 + t374;
t320 = -cos(pkin(10)) * pkin(1) - pkin(2) - pkin(7);
t298 = qJD(1) * t320 + qJD(3);
t329 = t339 * qJD(2);
t272 = t336 * t298 + t329;
t264 = qJD(4) * pkin(8) + t272;
t322 = sin(pkin(10)) * pkin(1) + qJ(3);
t295 = pkin(4) * t336 - pkin(8) * t339 + t322;
t275 = t295 * qJD(1);
t423 = t275 * t335;
t240 = t264 * t338 + t423;
t235 = -pkin(9) * t301 + t240;
t390 = qJD(6) * t334;
t232 = t235 * t390;
t437 = -t336 * qJD(2) + t298 * t339;
t263 = -qJD(4) * pkin(4) - t437;
t244 = pkin(5) * t301 + t263;
t442 = t244 * t248 + t232;
t376 = t335 * t391;
t345 = -t336 * t386 - t376;
t268 = qJD(1) * t345 + qJD(5) * t386;
t265 = t437 * qJD(4);
t358 = pkin(4) * t339 + pkin(8) * t336;
t299 = qJD(4) * t358 + qJD(3);
t280 = t299 * qJD(1);
t364 = -t265 * t335 + t338 * t280;
t343 = -qJD(5) * t240 + t364;
t385 = qJD(1) * qJD(4);
t371 = t339 * t385;
t224 = pkin(5) * t371 - pkin(9) * t268 + t343;
t373 = t335 * t399;
t269 = -qJD(4) * t373 + qJD(5) * t303;
t392 = qJD(5) * t338;
t382 = -t338 * t265 - t275 * t392 - t335 * t280;
t393 = qJD(5) * t335;
t348 = -t264 * t393 - t382;
t225 = -pkin(9) * t269 + t348;
t367 = t337 * t224 - t334 * t225;
t441 = t244 * t354 + t367;
t369 = MDP(26) * t398;
t440 = qJD(4) * t369 + (-t248 ^ 2 + t354 ^ 2) * MDP(23) - t248 * t354 * MDP(22);
t331 = t339 ^ 2;
t439 = MDP(9) * (t336 ^ 2 - t331);
t309 = qJD(1) * t322;
t438 = t309 * MDP(7);
t366 = t268 * t334 + t337 * t269;
t228 = -qJD(6) * t354 + t366;
t435 = 0.2e1 * qJD(3);
t434 = pkin(8) + pkin(9);
t433 = t228 * t336;
t239 = -t264 * t335 + t338 * t275;
t234 = -pkin(9) * t303 + t239;
t231 = pkin(5) * t321 + t234;
t432 = t231 * t337;
t431 = t235 * t337;
t349 = t304 * t336;
t236 = qJD(4) * t349 - t339 * t342;
t430 = t236 * t318;
t429 = t263 * t335;
t428 = t263 * t338;
t266 = qJD(4) * t329 + t298 * t395;
t427 = t266 * t335;
t426 = t266 * t338;
t425 = t268 * t335;
t424 = t269 * t336;
t421 = t301 * t321;
t420 = t303 * t321;
t418 = t321 * t335;
t417 = t321 * t338;
t416 = t335 * t336;
t415 = t335 * t339;
t414 = t336 * t338;
t340 = qJD(4) ^ 2;
t413 = t336 * t340;
t412 = t338 * t339;
t411 = t339 * t268;
t410 = t339 * t340;
t389 = qJD(6) * t337;
t381 = t337 * t268 - t334 * t269 - t301 * t389;
t227 = -t303 * t390 + t381;
t394 = qJD(4) * t339;
t409 = t227 * t336 - t354 * t394;
t237 = -t390 * t415 + (t412 * t436 - t374) * t337 + t345 * t334;
t283 = t305 * t339;
t408 = -t237 * t318 - t283 * t371;
t407 = -qJD(1) * t349 - t444;
t347 = t305 * qJD(1);
t406 = t336 * t347 + t342;
t405 = t268 * t336 + t303 * t394;
t306 = t358 * qJD(1);
t404 = t335 * t306 + t338 * t437;
t300 = t320 * t414;
t403 = t335 * t295 + t300;
t400 = qJD(1) * t331;
t397 = qJD(4) * t318;
t387 = t263 * qJD(5);
t384 = pkin(9) * t414;
t383 = t320 * t416;
t380 = t339 * t320 * t386 + t295 * t392 + t335 * t299;
t379 = qJD(5) * t434;
t378 = t335 * t400;
t370 = MDP(19) * t398;
t368 = -t320 * t335 + pkin(5);
t365 = t320 * t321 + t264;
t363 = t338 * t306 - t335 * t437;
t362 = qJD(6) * t231 + t225;
t361 = qJD(5) * t336 + qJD(1);
t360 = t321 * t376;
t317 = t336 * t371;
t359 = -t272 + (t373 + t393) * pkin(5);
t315 = t434 * t338;
t357 = qJD(6) * t315 + (pkin(5) * t339 + t384) * qJD(1) + t363 + t338 * t379;
t314 = t434 * t335;
t356 = pkin(9) * t373 + qJD(6) * t314 + t335 * t379 + t404;
t222 = t231 * t334 + t431;
t282 = t338 * t295;
t243 = -pkin(9) * t412 + t336 * t368 + t282;
t245 = -pkin(9) * t415 + t403;
t355 = t243 * t334 + t245 * t337;
t353 = -t321 * t336 + t400;
t351 = t344 * t321;
t350 = -pkin(8) * t394 + t263 * t336;
t346 = t304 * qJD(1);
t341 = qJD(1) ^ 2;
t326 = -pkin(5) * t338 - pkin(4);
t290 = t338 * t299;
t285 = (pkin(5) * t335 - t320) * t339;
t284 = t304 * t339;
t255 = -pkin(5) * t344 + t320 * t395;
t241 = pkin(5) * t269 + t266;
t230 = pkin(9) * t344 - qJD(5) * t383 + t380;
t229 = t290 + (-t300 + (pkin(9) * t339 - t295) * t335) * qJD(5) + (t339 * t368 + t384) * qJD(4);
t221 = -t235 * t334 + t432;
t1 = [qJD(1) * MDP(6) * t435 + t435 * t438 - 0.2e1 * MDP(8) * t317 + 0.2e1 * t385 * t439 - MDP(10) * t413 - MDP(11) * t410 + (t309 * t394 - t320 * t413 + (t322 * t394 + t336 * t435) * qJD(1)) * MDP(13) + (-t309 * t395 - t320 * t410 + (-t322 * t395 + t339 * t435) * qJD(1)) * MDP(14) + (t303 * t345 + t338 * t411) * MDP(15) + ((t301 * t338 + t303 * t335) * t395 + (-t425 - t269 * t338 + (t301 * t335 - t303 * t338) * qJD(5)) * t339) * MDP(16) + (t353 * t386 - t360 + t405) * MDP(17) + (-t424 + (-t301 * t339 - t378) * qJD(4) + t351) * MDP(18) + (t321 * t394 + t317) * MDP(19) + ((-t295 * t393 + t290) * t321 + ((t301 * t320 - t429) * qJD(4) + (-t338 * t365 - t423) * qJD(5) + t364) * t336 + (t338 * t387 + t427 - t320 * t269 + (-t320 * t418 + (t282 - t383) * qJD(1) + t239) * qJD(4)) * t339) * MDP(20) + (-t380 * t321 + (t365 * t393 + (t303 * t320 - t428) * qJD(4) + t382) * t336 + (-t335 * t387 + t426 - t320 * t268 + (-qJD(1) * t403 - t240) * qJD(4)) * t339) * MDP(21) + (-t227 * t284 - t236 * t354) * MDP(22) + (-t227 * t283 + t228 * t284 - t236 * t248 + t237 * t354) * MDP(23) + (-t284 * t371 + t409 + t430) * MDP(24) + (-t248 * t394 + t408 - t433) * MDP(25) + (t318 * t394 + t317) * MDP(26) + ((t229 * t337 - t230 * t334) * t318 + t367 * t336 + t255 * t248 + t285 * t228 + t241 * t283 + t244 * t237 + (-t222 * t336 - t318 * t355) * qJD(6) + ((t243 * t337 - t245 * t334) * qJD(1) + t221) * t394) * MDP(27) + (t285 * t227 + t232 * t336 + t244 * t236 - t241 * t284 - t255 * t354 + (-(-qJD(6) * t245 + t229) * t318 - t224 * t336) * t334 + (-(qJD(6) * t243 + t230) * t318 - t362 * t336) * t337 + (-qJD(1) * t355 - t222) * t394) * MDP(28); (t351 + t424) * MDP(20) + (t360 + t405) * MDP(21) + (t408 + t433) * MDP(27) + (t409 - t430) * MDP(28) + (-MDP(13) * t339 + MDP(14) * t336) * t340 + (-MDP(20) * t378 + (MDP(28) * qJD(1) * t284 + t301 * MDP(20) + t248 * MDP(27)) * t339 - t353 * MDP(21) * t338) * qJD(4); -t341 * MDP(6) - qJD(1) * t438 + (-t339 * t269 - t361 * t417 + (t301 * t336 + (-t321 - t399) * t415) * qJD(4)) * MDP(20) + (-t411 + t361 * t418 + (-t321 * t412 + (t303 - t377) * t336) * qJD(4)) * MDP(21) + (t318 * t346 + (-t305 * t397 - t228) * t339 + ((-t305 * t398 + t248) * qJD(4) + t318 * t444) * t336) * MDP(27) + (t318 * t347 + (t304 * t397 - t227) * t339 + (t342 * t318 + (t339 * t346 - t354) * qJD(4)) * t336) * MDP(28) + (MDP(13) * t336 + MDP(14) * t339) * (-t340 - t341); (qJD(4) * t272 - t309 * t398 - t266) * MDP(13) + t309 * t399 * MDP(14) + (t303 * t417 + t425) * MDP(15) + ((t268 - t421) * t338 + (-t269 - t420) * t335) * MDP(16) + (t321 * t392 + (t321 * t414 + (-t303 + t396) * t339) * qJD(1)) * MDP(17) + (-t321 * t393 + (-t321 * t416 + (t301 + t386) * t339) * qJD(1)) * MDP(18) - t321 * t370 + (-pkin(4) * t269 - t426 - t363 * t321 - t272 * t301 + (-pkin(8) * t417 + t429) * qJD(5) + (-t239 * t339 + t335 * t350) * qJD(1)) * MDP(20) + (-pkin(4) * t268 + t427 + t404 * t321 - t272 * t303 + (pkin(8) * t418 + t428) * qJD(5) + (t240 * t339 + t338 * t350) * qJD(1)) * MDP(21) + (t227 * t305 - t354 * t407) * MDP(22) + (-t227 * t304 - t228 * t305 - t248 * t407 + t354 * t406) * MDP(23) + (t407 * t318 + (qJD(4) * t305 + t354) * t398) * MDP(24) + (-t406 * t318 + (-qJD(4) * t304 + t248) * t398) * MDP(25) - t318 * t369 + (t326 * t228 + t241 * t304 + (t334 * t356 - t337 * t357) * t318 + t359 * t248 + t406 * t244 + ((-t314 * t337 - t315 * t334) * qJD(4) - t221) * t398) * MDP(27) + (t326 * t227 + t241 * t305 + (t334 * t357 + t337 * t356) * t318 - t359 * t354 + t407 * t244 + (-(-t314 * t334 + t315 * t337) * qJD(4) + t222) * t398) * MDP(28) + (t339 * t336 * MDP(8) - t439) * t341; t303 * t301 * MDP(15) + (-t301 ^ 2 + t303 ^ 2) * MDP(16) + (t268 + t421) * MDP(17) + (-t269 + t420) * MDP(18) + qJD(4) * t370 + (t240 * t321 - t263 * t303 + t343) * MDP(20) + (t239 * t321 + t263 * t301 - t348) * MDP(21) + (t227 + t445) * MDP(24) + (-t228 - t443) * MDP(25) + (-(-t234 * t334 - t431) * t318 - t222 * qJD(6) + (-t248 * t303 - t318 * t390 + t337 * t371) * pkin(5) + t441) * MDP(27) + ((-t235 * t318 - t224) * t334 + (t234 * t318 - t362) * t337 + (t303 * t354 - t318 * t389 - t334 * t371) * pkin(5) + t442) * MDP(28) + t440; (t381 + t445) * MDP(24) + (-t366 - t443) * MDP(25) + (t222 * t318 + t441) * MDP(27) + (t221 * t318 - t334 * t224 - t337 * t225 + t442) * MDP(28) + (-MDP(24) * t419 + MDP(25) * t354 - MDP(27) * t222 - MDP(28) * t432) * qJD(6) + t440;];
tauc  = t1;
