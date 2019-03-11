% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:45
% EndTime: 2019-03-09 02:18:52
% DurationCPUTime: 3.42s
% Computational Cost: add. (3795->301), mult. (9338->407), div. (0->0), fcn. (7408->10), ass. (0->141)
t327 = cos(qJ(6));
t372 = qJD(6) * t327;
t320 = sin(pkin(11));
t326 = sin(qJ(4));
t322 = cos(pkin(11));
t329 = cos(qJ(4));
t386 = t322 * t329;
t340 = t320 * t326 - t386;
t297 = t340 * qJD(1);
t387 = t320 * t329;
t305 = t322 * t326 + t387;
t298 = t305 * qJD(1);
t325 = sin(qJ(5));
t328 = cos(qJ(5));
t272 = t328 * t297 + t298 * t325;
t417 = t272 * t327;
t422 = t372 + t417;
t343 = -t297 * t325 + t328 * t298;
t324 = sin(qJ(6));
t373 = qJD(6) * t324;
t376 = qJD(4) * t329;
t378 = qJD(1) * t322;
t308 = t376 * t378;
t377 = qJD(1) * t326;
t364 = t320 * t377;
t293 = -qJD(4) * t364 + t308;
t300 = t305 * qJD(4);
t294 = qJD(1) * t300;
t374 = qJD(5) * t328;
t375 = qJD(5) * t325;
t250 = t328 * t293 - t325 * t294 - t297 * t374 - t298 * t375;
t319 = qJD(4) + qJD(5);
t381 = t327 * t250 + t319 * t372;
t236 = -t343 * t373 + t381;
t235 = t236 * t327;
t268 = t319 * t324 + t327 * t343;
t399 = t250 * t324;
t237 = t268 * qJD(6) + t399;
t389 = t343 * t324;
t266 = -t327 * t319 + t389;
t421 = -t324 * t237 - t422 * t266 + t235;
t234 = t236 * t324;
t251 = t343 * qJD(5) + t293 * t325 + t328 * t294;
t247 = t324 * t251;
t415 = -qJD(6) - t272;
t382 = -t372 * t415 + t247;
t391 = t272 * t319;
t393 = t343 * t319;
t395 = t268 * t343;
t420 = (-t251 + t393) * MDP(19) - t272 ^ 2 * MDP(17) + (t272 * MDP(16) + MDP(17) * t343 + MDP(27) * t415) * t343 + (t250 + t391) * MDP(18) + (t422 * t268 + t234) * MDP(23) + (-t415 * t417 + t382 - t395) * MDP(25);
t312 = sin(pkin(10)) * pkin(1) + qJ(3);
t307 = t312 * qJD(1);
t316 = t322 * qJD(2);
t281 = t316 + (-pkin(7) * qJD(1) - t307) * t320;
t289 = t320 * qJD(2) + t322 * t307;
t282 = pkin(7) * t378 + t289;
t410 = t329 * t281 - t282 * t326;
t256 = -pkin(8) * t298 + t410;
t255 = qJD(4) * pkin(4) + t256;
t345 = -t281 * t326 - t282 * t329;
t257 = -pkin(8) * t297 - t345;
t398 = t257 * t325;
t229 = t255 * t328 - t398;
t227 = -pkin(5) * t319 - t229;
t419 = t227 * t272;
t306 = -cos(pkin(10)) * pkin(1) - pkin(3) * t322 - pkin(2);
t295 = t306 * qJD(1) + qJD(3);
t276 = pkin(4) * t297 + t295;
t418 = t272 * t276;
t356 = t415 * t324;
t253 = pkin(5) * t343 + pkin(9) * t272;
t396 = t266 * t343;
t249 = t327 * t251;
t409 = -t373 * t415 - t249;
t408 = qJD(3) * t297;
t397 = t257 * t328;
t230 = t255 * t325 + t397;
t244 = -pkin(8) * t294 + t410 * qJD(4) - t408;
t336 = t305 * qJD(3);
t334 = qJD(1) * t336;
t245 = -pkin(8) * t293 + t345 * qJD(4) - t334;
t360 = t244 * t325 - t328 * t245;
t216 = t230 * qJD(5) + t360;
t228 = pkin(9) * t319 + t230;
t239 = pkin(5) * t272 - pkin(9) * t343 + t276;
t348 = t228 * t324 - t239 * t327;
t407 = -t216 * t327 + t227 * t373 + t343 * t348;
t218 = t228 * t327 + t239 * t324;
t406 = t216 * t324 + t218 * t343 + t227 * t372;
t405 = -t343 * t276 - t360;
t359 = t245 * t325 - t257 * t375;
t215 = (qJD(5) * t255 + t244) * t328 + t359;
t402 = pkin(7) + t312;
t301 = t402 * t320;
t302 = t402 * t322;
t333 = -t301 * t376 + qJD(3) * t386 + (-qJD(3) * t320 - qJD(4) * t302) * t326;
t260 = -pkin(8) * t300 + t333;
t299 = t340 * qJD(4);
t342 = t301 * t326 - t302 * t329;
t331 = t342 * qJD(4) - t336;
t261 = pkin(8) * t299 + t331;
t264 = -pkin(8) * t305 - t301 * t329 - t302 * t326;
t265 = -pkin(8) * t340 - t342;
t346 = t264 * t328 - t265 * t325;
t219 = t346 * qJD(5) + t260 * t328 + t261 * t325;
t241 = t264 * t325 + t265 * t328;
t278 = t305 * t328 - t325 * t340;
t280 = pkin(4) * t340 + t306;
t341 = -t305 * t325 - t328 * t340;
t246 = -pkin(5) * t341 - pkin(9) * t278 + t280;
t258 = t341 * qJD(5) - t299 * t328 - t300 * t325;
t404 = t216 * t278 + t227 * t258 - t241 * t251 + (qJD(6) * t246 + t219) * t415 + (qJD(6) * t239 + t215) * t341;
t403 = pkin(4) * t298;
t401 = t227 * t278;
t400 = t246 * t251;
t394 = t268 * t324;
t259 = t278 * qJD(5) - t299 * t325 + t328 * t300;
t384 = -t236 * t341 + t268 * t259;
t380 = t320 ^ 2 + t322 ^ 2;
t367 = t278 * t247;
t366 = t278 * t249;
t362 = -pkin(4) * t319 - t255;
t357 = qJD(1) * t380;
t313 = pkin(4) * t325 + pkin(9);
t352 = qJD(6) * t313 + t253 + t403;
t231 = t256 * t325 + t397;
t351 = pkin(4) * t375 - t231;
t232 = t256 * t328 - t398;
t350 = -pkin(4) * t374 + t232;
t349 = -t251 * t313 + t419;
t347 = t237 * t341 - t259 * t266;
t344 = (-t307 * t320 + t316) * t320 - t289 * t322;
t339 = t272 * t356 - t409;
t338 = -t258 * t324 - t278 * t372;
t337 = -t258 * t327 + t278 * t373;
t314 = -pkin(4) * t328 - pkin(5);
t224 = pkin(4) * t300 + pkin(5) * t259 - pkin(9) * t258;
t222 = pkin(4) * t294 + pkin(5) * t251 - pkin(9) * t250;
t221 = t327 * t222;
t220 = t241 * qJD(5) + t260 * t325 - t261 * t328;
t1 = [(t293 * t305 - t298 * t299) * MDP(9) + (-t293 * t340 - t294 * t305 + t297 * t299 - t298 * t300) * MDP(10) + (t306 * t294 + t295 * t300) * MDP(14) + (t306 * t293 - t295 * t299) * MDP(15) + (t250 * t278 + t258 * t343) * MDP(16) + (t250 * t341 - t251 * t278 - t258 * t272 - t259 * t343) * MDP(17) + (t251 * t280 + t259 * t276 + (t272 * t300 - t294 * t341) * pkin(4)) * MDP(21) + (t250 * t280 + t258 * t276 + (t278 * t294 + t300 * t343) * pkin(4)) * MDP(22) + (t278 * t235 - t337 * t268) * MDP(23) + ((-t266 * t327 - t394) * t258 + (-t234 - t237 * t327 + (t266 * t324 - t268 * t327) * qJD(6)) * t278) * MDP(24) + (t337 * t415 + t366 + t384) * MDP(25) + (-t338 * t415 + t347 - t367) * MDP(26) + (-t251 * t341 - t259 * t415) * MDP(27) + (-t348 * t259 + t220 * t266 - t221 * t341 - t346 * t237 + (-t224 * t415 + t400 + (t228 * t341 + t241 * t415 + t401) * qJD(6)) * t327 + t404 * t324) * MDP(28) + (-t218 * t259 + t220 * t268 - t346 * t236 + ((-qJD(6) * t241 + t224) * t415 - t400 + (-qJD(6) * t228 + t222) * t341 - qJD(6) * t401) * t324 + t404 * t327) * MDP(29) + (t258 * MDP(18) - t259 * MDP(19) - t220 * MDP(21) - t219 * MDP(22)) * t319 + (0.2e1 * MDP(7) * t357 + (t312 * t357 - t344) * MDP(8)) * qJD(3) + (-t299 * MDP(11) - t300 * MDP(12) + MDP(14) * t331 - MDP(15) * t333) * qJD(4); (-t347 - t367) * MDP(28) + (-t366 + t384) * MDP(29) + (-MDP(21) * t259 - MDP(22) * t258) * t319 - (t338 * MDP(28) + t337 * MDP(29)) * t415 + (-MDP(14) * t300 + MDP(15) * t299) * qJD(4); t308 * MDP(15) + (t251 + t393) * MDP(21) + (t250 - t391) * MDP(22) + (t339 - t396) * MDP(28) + (-t327 * t415 ^ 2 - t247 - t395) * MDP(29) - t380 * MDP(7) * qJD(1) ^ 2 + t344 * MDP(8) * qJD(1) + ((qJD(1) * t387 + t322 * t377 + t298) * MDP(14) + (-t297 - t364) * MDP(15)) * qJD(4); t298 * t297 * MDP(9) + (-t297 ^ 2 + t298 ^ 2) * MDP(10) + (t308 + (t297 - t364) * qJD(4)) * MDP(11) + (-t295 * t298 - t334) * MDP(14) + (t295 * t297 + t408) * MDP(15) + (-t272 * t403 + t231 * t319 + (t362 * t325 - t397) * qJD(5) + t405) * MDP(21) + (-t343 * t403 + t232 * t319 + t418 + (t362 * qJD(5) - t244) * t328 - t359) * MDP(22) + (t394 * t415 + t421) * MDP(24) + (t339 + t396) * MDP(26) + (t314 * t237 + t349 * t324 + t351 * t266 - (t350 * t324 - t352 * t327) * t415 + t407) * MDP(28) + (t314 * t236 + t349 * t327 + t351 * t268 - (t352 * t324 + t350 * t327) * t415 + t406) * MDP(29) + t420; ((-qJD(5) + t319) * t230 + t405) * MDP(21) + (t229 * t319 - t215 + t418) * MDP(22) + (t268 * t356 + t421) * MDP(24) + (-t356 * t415 + t249 + t396) * MDP(26) + (-pkin(5) * t237 + (-t229 * t324 + t253 * t327) * t415 - t230 * t266 + t324 * t419 - t382 * pkin(9) + t407) * MDP(28) + (-pkin(5) * t236 - (t229 * t327 + t253 * t324) * t415 - t230 * t268 + t227 * t417 + t409 * pkin(9) + t406) * MDP(29) + t420; t268 * t266 * MDP(23) + (-t266 ^ 2 + t268 ^ 2) * MDP(24) + (-t266 * t415 + t381) * MDP(25) + (-t268 * t415 - t399) * MDP(26) + t251 * MDP(27) + (-t215 * t324 - t218 * t415 - t227 * t268 + t221) * MDP(28) + (-t215 * t327 - t222 * t324 + t227 * t266 + t348 * t415) * MDP(29) + (-MDP(25) * t389 - t268 * MDP(26) - t218 * MDP(28) + t348 * MDP(29)) * qJD(6);];
tauc  = t1;
