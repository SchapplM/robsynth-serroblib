% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:56
% EndTime: 2019-12-05 16:33:04
% DurationCPUTime: 2.70s
% Computational Cost: add. (1435->308), mult. (3843->465), div. (0->0), fcn. (2866->10), ass. (0->142)
t301 = sin(qJ(3));
t379 = MDP(5) * t301;
t304 = cos(qJ(3));
t378 = (t301 ^ 2 - t304 ^ 2) * MDP(6);
t324 = pkin(3) * t301 - qJ(4) * t304;
t256 = t324 * qJD(3) - qJD(4) * t301;
t296 = sin(pkin(10));
t348 = qJD(3) * t301;
t341 = pkin(7) * t348;
t287 = t296 * t341;
t298 = cos(pkin(10));
t302 = sin(qJ(2));
t297 = sin(pkin(5));
t354 = qJD(1) * t297;
t305 = cos(qJ(2));
t364 = t304 * t305;
t360 = t298 * t256 - (-t296 * t364 + t298 * t302) * t354 + t287;
t377 = t296 * t256 - (t296 * t302 + t298 * t364) * t354;
t339 = t302 * t354;
t280 = qJD(2) * pkin(7) + t339;
t299 = cos(pkin(5));
t353 = qJD(1) * t299;
t376 = -t301 * t280 + t304 * t353;
t375 = pkin(8) + qJ(4);
t374 = qJD(2) * pkin(2);
t288 = t301 * t353;
t352 = qJD(2) * t297;
t334 = t305 * t352;
t328 = qJD(1) * t334;
t347 = qJD(3) * t304;
t220 = qJD(3) * t288 + t280 * t347 + t301 * t328;
t373 = t220 * t296;
t372 = t220 * t298;
t300 = sin(qJ(5));
t371 = t296 * t300;
t370 = t296 * t304;
t369 = t297 * t302;
t368 = t297 * t305;
t367 = t298 * t301;
t366 = t298 * t304;
t306 = qJD(3) ^ 2;
t365 = t301 * t306;
t363 = t304 * t306;
t218 = t304 * t328 + (qJD(4) + t376) * qJD(3);
t228 = (t256 + t339) * qJD(2);
t199 = t298 * t218 + t296 * t228;
t318 = pkin(4) * t301 - pkin(8) * t366;
t312 = t318 * qJD(3);
t362 = t312 + t360;
t361 = (-pkin(7) * t367 - pkin(8) * t370) * qJD(3) + t377;
t245 = t304 * t280 + t288;
t241 = qJD(3) * qJ(4) + t245;
t282 = -pkin(3) * t304 - qJ(4) * t301 - pkin(2);
t338 = t305 * t354;
t246 = t282 * qJD(2) - t338;
t206 = t298 * t241 + t296 * t246;
t329 = t298 * t341;
t359 = -t329 + t377;
t276 = t324 * qJD(2);
t214 = t296 * t276 + t298 * t376;
t303 = cos(qJ(5));
t274 = -t303 * t298 + t371;
t313 = t274 * t304;
t358 = qJD(2) * t313 - t274 * qJD(5);
t275 = t296 * t303 + t298 * t300;
t314 = t275 * t304;
t357 = -qJD(2) * t314 + t275 * qJD(5);
t343 = t298 * qJD(3);
t351 = qJD(2) * t301;
t270 = t296 * t351 - t343;
t342 = qJD(2) * qJD(3);
t333 = t304 * t342;
t327 = t298 * t333;
t344 = qJD(5) * t303;
t356 = -t270 * t344 + t303 * t327;
t248 = pkin(7) * t366 + t296 * t282;
t350 = qJD(2) * t304;
t349 = qJD(3) * t296;
t201 = -pkin(8) * t270 + t206;
t346 = qJD(5) * t201;
t345 = qJD(5) * t301;
t340 = pkin(4) * t296 + pkin(7);
t336 = t296 * t350;
t335 = t302 * t352;
t332 = MDP(20) * t348;
t331 = -qJD(3) * pkin(3) + qJD(4);
t198 = -t218 * t296 + t298 * t228;
t205 = -t241 * t296 + t298 * t246;
t213 = t298 * t276 - t296 * t376;
t326 = t296 * t333;
t197 = -pkin(8) * t326 + t199;
t272 = t298 * t351 + t349;
t200 = -pkin(4) * t350 - pkin(8) * t272 + t205;
t330 = -qJD(5) * t200 - t197;
t281 = -t338 - t374;
t325 = -t281 - t338;
t193 = t200 * t303 - t201 * t300;
t194 = t200 * t300 + t201 * t303;
t269 = t298 * t282;
t226 = -pkin(8) * t367 + t269 + (-pkin(7) * t296 - pkin(4)) * t304;
t236 = -pkin(8) * t296 * t301 + t248;
t323 = t226 * t303 - t236 * t300;
t322 = t226 * t300 + t236 * t303;
t262 = t299 * t301 + t304 * t369;
t232 = -t262 * t296 - t298 * t368;
t233 = t262 * t298 - t296 * t368;
t321 = t232 * t303 - t233 * t300;
t320 = t232 * t300 + t233 * t303;
t319 = t270 * t300 - t272 * t303;
t261 = -t299 * t304 + t301 * t369;
t285 = t375 * t298;
t316 = t318 * qJD(2) + qJD(4) * t296 + qJD(5) * t285 + t213;
t284 = t375 * t296;
t315 = pkin(8) * t336 + qJD(4) * t298 - qJD(5) * t284 - t214;
t311 = -qJD(5) * t272 - t326;
t310 = qJD(3) * t314;
t237 = -t376 + t331;
t309 = qJD(3) * (-t325 - t374);
t308 = -qJ(4) * t348 + (-t237 + t331) * t304;
t203 = qJD(2) * t310 - t319 * qJD(5);
t307 = qJD(2) ^ 2;
t292 = -pkin(4) * t298 - pkin(3);
t290 = -qJD(5) + t350;
t277 = t340 * t301;
t265 = t340 * t347;
t257 = t303 * t270;
t254 = t274 * t301;
t253 = t275 * t301;
t247 = -pkin(7) * t370 + t269;
t234 = -t261 * qJD(3) + t304 * t334;
t227 = pkin(4) * t336 + t245;
t221 = t272 * t300 + t257;
t217 = t344 * t367 - t345 * t371 + t310;
t216 = -qJD(3) * t313 - t275 * t345;
t215 = pkin(4) * t270 + t237;
t212 = t234 * t298 + t296 * t335;
t211 = -t234 * t296 + t298 * t335;
t208 = pkin(4) * t326 + t220;
t202 = t311 * t300 + t356;
t196 = qJD(2) * t312 + t198;
t195 = t303 * t196;
t1 = [-t234 * qJD(3) * MDP(11) + (-t211 * t272 - t212 * t270) * MDP(14) + (t198 * t232 + t199 * t233 + t205 * t211 + t206 * t212 + t220 * t261) * MDP(15) + (-(-t320 * qJD(5) + t211 * t303 - t212 * t300) * t290 + t261 * t203) * MDP(21) + (t261 * t202 + (t321 * qJD(5) + t211 * t300 + t212 * t303) * t290) * MDP(22) + (-MDP(4) * t305 + (-MDP(10) * t304 + MDP(11) * t301 - MDP(3)) * t302) * t307 * t297 + (-MDP(10) * qJD(3) + MDP(12) * t270 + MDP(13) * t272 + MDP(15) * t237 + MDP(21) * t221 - MDP(22) * t319) * (t262 * qJD(3) + t301 * t334) + ((-MDP(12) * t211 + MDP(13) * t212) * t304 + ((-MDP(11) * t368 + (-t232 * t298 - t233 * t296) * MDP(14) + (t296 * MDP(12) + t298 * MDP(13)) * t261) * t304 + (-MDP(10) * t368 + t232 * MDP(12) - t233 * MDP(13) + t321 * MDP(21) - t320 * MDP(22)) * t301) * qJD(3)) * qJD(2); 0.2e1 * t333 * t379 - 0.2e1 * t342 * t378 + MDP(7) * t363 - MDP(8) * t365 + (-pkin(7) * t363 + t301 * t309) * MDP(10) + (pkin(7) * t365 + t304 * t309) * MDP(11) + ((-t270 * t338 + t373 + (qJD(2) * t247 + t205) * qJD(3)) * t301 + (-t198 + (pkin(7) * t270 + t237 * t296) * qJD(3) + (t287 - t360) * qJD(2)) * t304) * MDP(12) + ((-t272 * t338 + t372 + (-qJD(2) * t248 - t206) * qJD(3)) * t301 + (t199 + (pkin(7) * t272 + t237 * t298) * qJD(3) + (t329 + t359) * qJD(2)) * t304) * MDP(13) + ((-t198 * t298 - t199 * t296) * t301 - t360 * t272 - t359 * t270 + (-t205 * t298 - t206 * t296 + (-t247 * t298 - t248 * t296) * qJD(2)) * t347) * MDP(14) + (-t237 * t301 * t338 + t198 * t247 + t199 * t248 + t359 * t206 + t360 * t205 + (t220 * t301 + t237 * t347) * pkin(7)) * MDP(15) + (-t202 * t254 - t216 * t319) * MDP(16) + (-t202 * t253 + t203 * t254 - t216 * t221 + t217 * t319) * MDP(17) + (-t202 * t304 - t216 * t290 + (-qJD(2) * t254 - t319) * t348) * MDP(18) + (t203 * t304 + t217 * t290 + (-qJD(2) * t253 - t221) * t348) * MDP(19) + (-t290 - t350) * t332 + (-(-t197 * t300 + t195) * t304 + t265 * t221 + t277 * t203 + t208 * t253 + t215 * t217 + (t361 * t300 - t362 * t303) * t290 + (t194 * t304 + t322 * t290) * qJD(5) + (-t221 * t338 + (t323 * qJD(2) + t193) * qJD(3)) * t301) * MDP(21) + (-t265 * t319 + t277 * t202 - t208 * t254 + t215 * t216 + (t196 * t300 + t197 * t303) * t304 + (t362 * t300 + t361 * t303) * t290 + (t193 * t304 + t323 * t290) * qJD(5) + (t319 * t338 + (-t322 * qJD(2) - t194) * qJD(3)) * t301) * MDP(22); (qJD(3) * t245 - t281 * t351 - t220) * MDP(10) + t325 * t350 * MDP(11) + (-t372 - t245 * t270 + (-t205 * t301 + t213 * t304 + t308 * t296) * qJD(2)) * MDP(12) + (t373 - t245 * t272 + (t206 * t301 - t214 * t304 + t308 * t298) * qJD(2)) * MDP(13) + (t213 * t272 + t214 * t270 + (-qJD(4) * t270 + t205 * t350 + t199) * t298 + (qJD(4) * t272 + t206 * t350 - t198) * t296) * MDP(14) + (-pkin(3) * t220 - t205 * t213 - t206 * t214 - t237 * t245 + (-t205 * t296 + t206 * t298) * qJD(4) + (-t198 * t296 + t199 * t298) * qJ(4)) * MDP(15) + (t202 * t275 - t319 * t358) * MDP(16) + (-t202 * t274 - t203 * t275 - t358 * t221 + t319 * t357) * MDP(17) + (-t358 * t290 + (qJD(3) * t275 + t319) * t351) * MDP(18) + (t357 * t290 + (-qJD(3) * t274 + t221) * t351) * MDP(19) + t290 * MDP(20) * t351 + (t292 * t203 + t208 * t274 - t227 * t221 + (t315 * t300 + t316 * t303) * t290 + t357 * t215 + ((-t284 * t303 - t285 * t300) * qJD(3) - t193) * t351) * MDP(21) + (t292 * t202 + t208 * t275 + t227 * t319 + (-t316 * t300 + t315 * t303) * t290 + t358 * t215 + (-(-t284 * t300 + t285 * t303) * qJD(3) + t194) * t351) * MDP(22) + (-t304 * t379 + t378) * t307; (-t270 ^ 2 - t272 ^ 2) * MDP(14) + (t205 * t272 + t206 * t270 + t220) * MDP(15) + (t319 * t290 + t203) * MDP(21) + (t257 * t290 + (-t326 + (-qJD(5) + t290) * t272) * t300 + t356) * MDP(22) + ((-t272 + t349) * MDP(12) + (t270 + t343) * MDP(13)) * t350; -t221 ^ 2 * MDP(17) + (-t221 * t290 + t356) * MDP(18) + qJD(2) * t332 + (-t194 * t290 + t195) * MDP(21) + (-t193 * t290 + t215 * t221) * MDP(22) - (MDP(16) * t221 - MDP(17) * t319 - MDP(19) * t290 - MDP(21) * t215) * t319 + (t311 * MDP(19) - MDP(21) * t346 + t330 * MDP(22)) * t303 + (t311 * MDP(18) + (qJD(5) * t270 - t327) * MDP(19) + t330 * MDP(21) + (-t196 + t346) * MDP(22)) * t300;];
tauc = t1;
