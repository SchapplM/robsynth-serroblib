% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:49
% EndTime: 2019-03-09 01:53:58
% DurationCPUTime: 4.15s
% Computational Cost: add. (3057->352), mult. (6980->487), div. (0->0), fcn. (5164->8), ass. (0->149)
t335 = sin(pkin(10));
t337 = cos(pkin(10));
t336 = sin(pkin(9));
t341 = sin(qJ(4));
t379 = qJD(1) * t341;
t368 = t336 * t379;
t338 = cos(pkin(9));
t343 = cos(qJ(4));
t378 = qJD(1) * t343;
t370 = t338 * t378;
t413 = t368 - t370;
t282 = -t337 * qJD(4) - t335 * t413;
t342 = cos(qJ(6));
t414 = t342 * t282;
t340 = sin(qJ(6));
t310 = t335 * t342 + t337 * t340;
t303 = t310 * qJD(6);
t309 = t336 * t343 + t338 * t341;
t349 = qJD(1) * t309;
t382 = t310 * t349 + t303;
t412 = MDP(6) * qJ(2) + t336 * MDP(7) + t338 * MDP(8) + MDP(5);
t284 = qJD(4) * t335 - t337 * t413;
t411 = t282 * t340 - t284 * t342;
t297 = qJD(6) + t349;
t352 = t297 * t310;
t339 = -pkin(1) - qJ(3);
t409 = qJD(1) * t339;
t318 = qJD(2) + t409;
t367 = -pkin(7) * qJD(1) + t318;
t295 = t367 * t336;
t296 = t367 * t338;
t410 = -t341 * t295 + t296 * t343;
t404 = -pkin(7) + t339;
t311 = t404 * t336;
t312 = t404 * t338;
t278 = t341 * t311 - t312 * t343;
t381 = t336 ^ 2 + t338 ^ 2;
t408 = qJD(3) * t381;
t385 = t342 * t337;
t388 = t335 * t340;
t307 = -t385 + t388;
t407 = t307 * qJD(6);
t383 = -t307 * t349 - t407;
t317 = qJD(4) * t368;
t374 = qJD(4) * t343;
t369 = t338 * t374;
t294 = qJD(1) * t369 - t317;
t395 = t294 * t310;
t406 = -t297 * t383 - t395;
t405 = pkin(8) * t337;
t330 = t336 * pkin(3);
t403 = pkin(8) + qJ(5);
t246 = t284 * t340 + t414;
t401 = t246 * t413;
t400 = t411 * t413;
t293 = qJD(4) * t349;
t308 = t336 * t341 - t343 * t338;
t399 = t293 * t308;
t398 = t293 * t335;
t397 = t293 * t337;
t396 = t294 * t309;
t393 = t349 * t335;
t375 = qJD(4) * t341;
t304 = -t336 * t374 - t338 * t375;
t392 = t304 * t335;
t277 = t307 * t294;
t391 = t308 * t335;
t390 = t308 * t337;
t323 = qJ(2) + t330;
t350 = t309 * qJD(3);
t240 = -qJD(1) * t350 + (qJD(5) + t410) * qJD(4);
t333 = qJD(1) * qJD(2);
t245 = pkin(4) * t294 + qJ(5) * t293 + qJD(5) * t413 + t333;
t217 = t337 * t240 + t335 * t245;
t305 = -t336 * t375 + t369;
t251 = pkin(4) * t305 - qJ(5) * t304 + qJD(5) * t308 + qJD(2);
t257 = -qJD(4) * t278 - t350;
t224 = t335 * t251 + t337 * t257;
t266 = t343 * t295 + t341 * t296;
t261 = qJD(4) * qJ(5) + t266;
t334 = qJD(1) * qJ(2);
t328 = qJD(3) + t334;
t313 = qJD(1) * t330 + t328;
t262 = pkin(4) * t349 + qJ(5) * t413 + t313;
t228 = t337 * t261 + t335 * t262;
t274 = -pkin(4) * t413 + qJ(5) * t349;
t232 = t335 * t274 + t337 * t410;
t273 = pkin(4) * t309 + qJ(5) * t308 + t323;
t279 = t311 * t343 + t312 * t341;
t237 = t335 * t273 + t337 * t279;
t372 = qJD(6) * t342;
t384 = -t282 * t372 - t293 * t385;
t380 = qJD(1) * t349;
t377 = qJD(4) * t304;
t376 = qJD(4) * t305;
t219 = -pkin(8) * t282 + t228;
t373 = qJD(6) * t219;
t216 = -t240 * t335 + t337 * t245;
t223 = t337 * t251 - t257 * t335;
t227 = -t261 * t335 + t337 * t262;
t231 = t337 * t274 - t335 * t410;
t236 = t337 * t273 - t279 * t335;
t366 = qJD(1) * t381;
t213 = pkin(8) * t398 + t217;
t215 = pkin(5) * t349 - pkin(8) * t284 + t227;
t365 = -qJD(6) * t215 - t213;
t364 = -t297 * t382 - t277;
t209 = t215 * t342 - t219 * t340;
t210 = t215 * t340 + t219 * t342;
t363 = t216 * t309 + t227 * t305;
t362 = -t217 * t309 - t228 * t305;
t226 = pkin(5) * t309 + pkin(8) * t390 + t236;
t229 = pkin(8) * t391 + t237;
t361 = t226 * t342 - t229 * t340;
t360 = t226 * t340 + t229 * t342;
t359 = -t227 * t335 + t228 * t337;
t256 = -qJD(4) * pkin(4) + qJD(5) - t410;
t346 = qJD(3) * t413 - t295 * t374 - t296 * t375;
t358 = -t256 * t304 - t308 * t346;
t357 = -t282 * t337 + t284 * t335;
t355 = -qJD(6) * t284 + t398;
t315 = t403 * t337;
t354 = -pkin(5) * t413 + qJD(5) * t335 + qJD(6) * t315 + t349 * t405 + t231;
t314 = t403 * t335;
t353 = pkin(8) * t393 - qJD(5) * t337 + qJD(6) * t314 + t232;
t351 = t297 * t307;
t348 = -t278 * t293 - t358;
t347 = -t305 * t349 - t396 - t399;
t345 = pkin(4) * t293 - qJ(5) * t294 + (-qJD(5) + t256) * t349;
t222 = -qJD(6) * t411 - t293 * t310;
t258 = -qJD(3) * t308 + qJD(4) * t279;
t344 = qJD(1) ^ 2;
t326 = -pkin(5) * t337 - pkin(4);
t298 = t349 ^ 2;
t272 = t307 * t308;
t271 = t310 * t308;
t259 = -pkin(5) * t391 + t278;
t244 = -pkin(5) * t393 + t266;
t239 = pkin(5) * t392 + t258;
t235 = pkin(5) * t282 + t256;
t234 = qJD(6) * t308 * t388 + t304 * t310 - t372 * t390;
t233 = t303 * t308 - t304 * t307;
t230 = -pkin(5) * t398 - t346;
t221 = t355 * t340 + t384;
t218 = -pkin(8) * t392 + t224;
t214 = pkin(5) * t305 - t304 * t405 + t223;
t212 = pkin(5) * t294 + pkin(8) * t397 + t216;
t211 = t342 * t212;
t1 = [0.2e1 * qJD(3) * MDP(9) * t366 + ((t328 + t334) * qJD(2) + (-t318 - t409) * t408) * MDP(10) + (-t304 * t413 + t399) * MDP(11) + (t293 * t309 + t294 * t308 - t304 * t349 + t305 * t413) * MDP(12) + MDP(13) * t377 - MDP(14) * t376 + (0.2e1 * t349 * qJD(2) - qJD(4) * t258 + t294 * t323 + t305 * t313) * MDP(16) + (-qJD(4) * t257 - t293 * t323 + t304 * t313 + (-qJD(1) * t308 - t413) * qJD(2)) * MDP(17) + (t223 * t349 + t236 * t294 + t258 * t282 + t335 * t348 + t363) * MDP(18) + (-t224 * t349 - t237 * t294 + t258 * t284 + t337 * t348 + t362) * MDP(19) + (-t223 * t284 - t224 * t282 + (t216 * t308 - t227 * t304 + t236 * t293) * t337 + (t217 * t308 - t228 * t304 + t237 * t293) * t335) * MDP(20) + (t216 * t236 + t217 * t237 + t223 * t227 + t224 * t228 + t256 * t258 - t278 * t346) * MDP(21) + (t221 * t272 - t233 * t411) * MDP(22) + (t221 * t271 - t222 * t272 - t233 * t246 + t234 * t411) * MDP(23) + (t221 * t309 + t233 * t297 + t272 * t294 - t305 * t411) * MDP(24) + (-t222 * t309 - t234 * t297 - t246 * t305 + t271 * t294) * MDP(25) + (t297 * t305 + t396) * MDP(26) + ((t214 * t342 - t218 * t340) * t297 + t361 * t294 + (-t213 * t340 + t211) * t309 + t209 * t305 + t239 * t246 + t259 * t222 - t230 * t271 + t235 * t234 + (-t210 * t309 - t297 * t360) * qJD(6)) * MDP(27) + (-(t214 * t340 + t218 * t342) * t297 - t360 * t294 - (t212 * t340 + t213 * t342) * t309 - t210 * t305 - t239 * t411 + t259 * t221 + t230 * t272 + t235 * t233 + (-t209 * t309 - t297 * t361) * qJD(6)) * MDP(28) + 0.2e1 * t412 * t333; (-t328 - t408) * qJD(1) * MDP(10) + (t377 - t380) * MDP(16) + (qJD(1) * t413 - t376) * MDP(17) + (-t282 * t304 + t335 * t347 - t337 * t380) * MDP(18) + (-t284 * t304 + t335 * t380 + t337 * t347) * MDP(19) + (t357 * t305 + (t282 * t335 + t284 * t337) * qJD(1)) * MDP(20) + ((-qJD(1) * t227 - t362) * t337 + (-qJD(1) * t228 - t363) * t335 + t358) * MDP(21) + (t308 * t222 - t304 * t246 - t305 * t352 + qJD(1) * t351 + (t297 * t407 - t395) * t309) * MDP(27) + (t308 * t221 + t304 * t411 + t305 * t351 + qJD(1) * t352 + (qJD(6) * t352 + t277) * t309) * MDP(28) - t412 * t344; (t318 * t366 + t333) * MDP(10) - t317 * MDP(16) + (t282 * t413 + t294 * t337) * MDP(18) + (t284 * t413 - t294 * t335 - t298 * t337) * MDP(19) + (t216 * t337 + t217 * t335 + t256 * t413) * MDP(21) + (t364 + t401) * MDP(27) + (-t400 + t406) * MDP(28) - t381 * MDP(9) * t344 - (-t335 ^ 2 - t337 ^ 2) * MDP(20) * t293 + (-MDP(18) * t393 + MDP(20) * t357 + t359 * MDP(21)) * t349 + ((-t413 + t370) * MDP(16) + (-t336 * t378 - t338 * t379 - t349) * MDP(17)) * qJD(4); (t413 ^ 2 - t298) * MDP(12) + (t317 + (-t413 - t370) * qJD(4)) * MDP(14) + (qJD(4) * t266 + t313 * t413 + t346) * MDP(16) + (t227 * t413 - t266 * t282 + t335 * t345 + t337 * t346) * MDP(18) + (-t228 * t413 - t266 * t284 - t335 * t346 + t337 * t345) * MDP(19) + (t231 * t284 + t232 * t282 + (-qJD(5) * t282 + t217) * t337 + (qJD(5) * t284 - t216) * t335) * MDP(20) + (pkin(4) * t346 - t227 * t231 - t228 * t232 - t256 * t266 + t359 * qJD(5) + (-t216 * t335 + t217 * t337) * qJ(5)) * MDP(21) + (t221 * t310 - t383 * t411) * MDP(22) + (-t221 * t307 - t222 * t310 - t246 * t383 + t382 * t411) * MDP(23) + (-t400 - t406) * MDP(24) + (t364 - t401) * MDP(25) + t297 * t413 * MDP(26) + ((-t314 * t342 - t315 * t340) * t294 + t326 * t222 + t230 * t307 + t209 * t413 - t244 * t246 + (t340 * t353 - t342 * t354) * t297 + t382 * t235) * MDP(27) + (-(-t314 * t340 + t315 * t342) * t294 + t326 * t221 + t230 * t310 - t210 * t413 + t244 * t411 + (t340 * t354 + t342 * t353) * t297 + t383 * t235) * MDP(28) + (-t413 * MDP(11) + (qJD(3) + t313) * MDP(17) - t231 * MDP(18) + t232 * MDP(19) + (-t227 * t337 - t228 * t335) * MDP(20)) * t349; (t284 * t349 - t398) * MDP(18) + (-t282 * t349 - t397) * MDP(19) + (-t282 ^ 2 - t284 ^ 2) * MDP(20) + (t227 * t284 + t228 * t282 - t346) * MDP(21) + (-t297 * t411 + t222) * MDP(27) + (-t297 * t414 + (-t284 * t297 + t355) * t340 + t384) * MDP(28); -t246 ^ 2 * MDP(23) + (t246 * t297 + t384) * MDP(24) + t294 * MDP(26) + (t210 * t297 + t211) * MDP(27) + (t209 * t297 + t235 * t246) * MDP(28) - (MDP(22) * t246 - MDP(23) * t411 + MDP(25) * t297 - MDP(27) * t235) * t411 + (MDP(25) * t355 - MDP(27) * t373 + MDP(28) * t365) * t342 + (t355 * MDP(24) + (qJD(6) * t282 + t397) * MDP(25) + t365 * MDP(27) + (-t212 + t373) * MDP(28)) * t340;];
tauc  = t1;
