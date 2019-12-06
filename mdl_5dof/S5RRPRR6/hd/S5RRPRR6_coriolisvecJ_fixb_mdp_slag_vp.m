% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:36:15
% EndTime: 2019-12-05 18:36:20
% DurationCPUTime: 1.79s
% Computational Cost: add. (1619->219), mult. (2929->327), div. (0->0), fcn. (1839->8), ass. (0->138)
t297 = sin(qJ(5));
t295 = sin(pkin(9));
t298 = sin(qJ(4));
t368 = qJD(4) * t298;
t348 = t295 * t368;
t380 = t297 * t298;
t353 = t295 * t380;
t408 = -qJD(5) * t353 - t297 * t348;
t301 = cos(qJ(4));
t407 = t298 * MDP(13) + t301 * MDP(14);
t290 = t295 ^ 2;
t389 = t290 * t301;
t406 = -MDP(11) * t298 * t389 + (t298 ^ 2 - t301 ^ 2) * MDP(12) * t290;
t296 = cos(pkin(9));
t299 = sin(qJ(2));
t302 = cos(qJ(2));
t405 = t302 * MDP(6) + t299 * (t296 * MDP(7) + MDP(5));
t394 = pkin(1) * qJD(1);
t327 = -t302 * t394 + qJD(3);
t292 = qJD(1) + qJD(2);
t388 = t292 * t295;
t273 = -t296 * pkin(3) - t295 * pkin(7) - pkin(2);
t367 = qJD(4) * t301;
t370 = qJD(3) * t296;
t381 = t296 * t302;
t402 = -(t298 * t299 + t301 * t381) * t394 + t273 * t367 + t301 * t370;
t365 = qJD(4) + qJD(5);
t382 = t296 * t301;
t358 = qJ(3) * t382;
t399 = (t273 * t298 + t358) * qJD(4) + (-t298 * t381 + t299 * t301) * t394 + t298 * t370;
t397 = pkin(1) * t299;
t396 = pkin(1) * t302;
t395 = pkin(8) * t295;
t393 = pkin(1) * qJD(2);
t271 = qJ(3) * t292 + t299 * t394;
t386 = t292 * t298;
t243 = (pkin(4) * t386 + t271) * t295;
t300 = cos(qJ(5));
t378 = t300 * t301;
t352 = t295 * t378;
t333 = t292 * t352;
t244 = t292 * t353 - t333;
t392 = t243 * t244;
t387 = t292 * t296;
t280 = -qJD(4) + t387;
t272 = -qJD(5) + t280;
t391 = t272 * t296;
t359 = qJD(1) * t393;
t267 = t292 * qJD(3) + t302 * t359;
t258 = t290 * t267;
t390 = t290 * t292;
t385 = t292 * t301;
t384 = t295 * t301;
t383 = t296 * t298;
t241 = t273 * t292 + t327;
t357 = t271 * t382;
t311 = -t241 * t298 - t357;
t364 = pkin(8) * t388;
t223 = -t298 * t364 - t311;
t379 = t300 * t223;
t346 = t271 * t368;
t334 = t299 * t359;
t350 = t241 * t367 + t267 * t382 + t298 * t334;
t377 = t267 * t389 + (-t296 * t346 + t350) * t296;
t291 = t296 ^ 2;
t376 = t291 * t267 + t258;
t374 = t408 * t292;
t373 = t290 + t291;
t366 = qJD(4) + t280;
t363 = pkin(8) * t384;
t361 = t299 * t393;
t360 = pkin(4) * t367;
t288 = qJ(3) + t397;
t356 = t288 * t382;
t355 = t290 * t385;
t354 = t292 * t384;
t307 = -pkin(8) * t354 - t271 * t383;
t209 = t307 * qJD(4) + t350;
t326 = -t267 * t383 + t301 * t334;
t210 = (-t357 + (-t241 + t364) * t298) * qJD(4) + t326;
t238 = t301 * t241;
t222 = t238 + t307;
t216 = -pkin(4) * t280 + t222;
t219 = qJD(5) * t297 * t223;
t318 = t297 * t301 + t298 * t300;
t303 = t365 * t318;
t229 = t303 * t295;
t239 = (t292 * t360 + t267) * t295;
t317 = t378 - t380;
t256 = t317 * t295;
t351 = (t297 * t210 - t219 + (qJD(5) * t216 + t209) * t300) * t296 - t243 * t229 + t239 * t256;
t262 = t273 - t396;
t285 = t302 * t393 + qJD(3);
t349 = t262 * t367 + t285 * t382 + t298 * t361;
t347 = t296 * t368;
t344 = pkin(4) * t272 - t216;
t343 = t262 - t395;
t342 = t273 - t395;
t341 = t373 * t302;
t340 = t373 * t267;
t339 = t373 * t285;
t338 = -t297 * t209 + t300 * t210;
t337 = t373 * qJD(3);
t336 = t388 * t397;
t335 = pkin(4) * t354;
t325 = -t285 * t383 + t301 * t361;
t319 = -t216 * t297 - t379;
t203 = t319 * qJD(5) + t338;
t313 = t365 * t352;
t230 = t313 + t408;
t255 = t318 * t295;
t324 = -t203 * t296 + t243 * t230 + t239 * t255;
t214 = t311 * qJD(4) + t326;
t323 = t290 * t271 * t367 - t214 * t296 + t298 * t258;
t321 = qJD(5) * (t342 * t301 + (-qJ(3) * t298 - pkin(4)) * t296) + (-qJ(3) * t383 - t363) * qJD(4) + t402;
t281 = pkin(8) * t348;
t320 = qJD(5) * (t342 * t298 + t358) - t281 + t399;
t282 = t295 * t360;
t316 = -t327 * t295 - t282;
t310 = -t262 * t298 - t356;
t224 = t292 * t229;
t245 = t318 * t388;
t308 = -t244 * t245 * MDP(18) + (-t245 * t272 - t224) * MDP(20) + (t244 * t272 - t365 * t333 - t374) * MDP(21) + (t244 ^ 2 - t245 ^ 2) * MDP(19);
t225 = t292 * t313 + t374;
t306 = (t224 * t255 - t225 * t256 + t229 * t245 + t230 * t244) * MDP(19) + (-t224 * t256 + t229 * t244) * MDP(18) + (t224 * t296 + t229 * t272) * MDP(20) + (t225 * t296 + t230 * t272) * MDP(21) + (0.2e1 * t406 * t292 + t407 * t295 * (t280 + t387)) * qJD(4);
t304 = t219 + (t223 * t272 - t210) * t297 + t243 * t245;
t289 = t292 ^ 2;
t287 = t298 * t295 * pkin(4);
t277 = t295 * t334;
t269 = qJ(3) * t295 + t287;
t268 = -t292 * pkin(2) + t327;
t257 = t288 * t295 + t287;
t247 = t285 * t295 + t282;
t233 = t343 * t298 + t356;
t228 = t343 * t301 + (-t288 * t298 - pkin(4)) * t296;
t221 = t310 * qJD(4) + t281 + t325;
t220 = (-t288 * t383 - t363) * qJD(4) + t349;
t1 = [(-(-t220 * t297 + t221 * t300 + (-t228 * t297 - t233 * t300) * qJD(5)) * t272 + t247 * t245 + t257 * t225 + t324) * MDP(23) + ((t220 * t300 + t221 * t297 + (t228 * t300 - t233 * t297) * qJD(5)) * t272 - t247 * t244 - t257 * t224 + t351) * MDP(24) + (t292 * t339 + t376) * MDP(9) + (t290 * t285 * t386 - t325 * t280 + (-t310 * t280 + t288 * t355) * qJD(4) + t323) * MDP(16) + (t271 * t339 + t288 * t340 + (t268 + (-pkin(2) - t396) * qJD(1)) * t361) * MDP(10) + (qJD(2) * t336 + t277) * MDP(8) + ((-t288 * t347 + t349) * t280 + (t285 * t385 + (-t288 * t292 - t271) * t368) * t290 + t377) * MDP(17) + t306 + t405 * (-qJD(1) - t292) * t393; (t269 * t225 + (t321 * t297 + t320 * t300) * t272 - t316 * t245 + t324) * MDP(23) + (t271 * t337 + qJ(3) * t340 + ((-pkin(2) * qJD(2) - t268) * t299 - t271 * t341) * t394) * MDP(10) + ((-t341 * t394 + t337) * t292 + t376) * MDP(9) + (-t269 * t224 + (-t320 * t297 + t321 * t300) * t272 + t316 * t244 + t351) * MDP(24) + (-qJD(1) * t336 + t277) * MDP(8) + (t399 * t280 + (qJ(3) * t367 + t327 * t298) * t390 + t323) * MDP(16) + t306 + ((-qJ(3) * t347 + t402) * t280 + (-t346 + (-qJ(3) * t368 + t327 * t301) * t292) * t290 + t377) * MDP(17) + t405 * (-qJD(2) + t292) * t394; -t373 * t289 * MDP(9) + (-t373 * t292 * t271 + t334) * MDP(10) + (t303 * t272 + (-t295 * t245 - t318 * t391) * t292) * MDP(23) + (t244 * t388 + (t365 * t272 - t391 * t292) * t317) * MDP(24) + (t298 * MDP(16) + t301 * MDP(17)) * (-t280 ^ 2 - t289 * t290); (-t271 * t355 + t311 * t280 + t214) * MDP(16) + (-t238 * t280 + (t366 * t296 + t390) * t298 * t271 - t350) * MDP(17) + ((-t222 * t297 - t379) * t272 - t245 * t335 + t392 + (t344 * t297 - t379) * qJD(5) + t338) * MDP(23) + (t244 * t335 + (t344 * qJD(5) - t222 * t272 - t209) * t300 + t304) * MDP(24) + t308 - t407 * t366 * t388 - t406 * t289; (t319 * t272 + t203 + t392) * MDP(23) + ((-t209 + (-qJD(5) - t272) * t216) * t300 + t304) * MDP(24) + t308;];
tauc = t1;
