% Calculate vector of inverse dynamics joint torques for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:51
% EndTime: 2019-12-05 17:47:56
% DurationCPUTime: 2.16s
% Computational Cost: add. (1473->261), mult. (2887->333), div. (0->0), fcn. (1982->10), ass. (0->131)
t325 = qJD(1) ^ 2;
t319 = sin(qJ(1));
t322 = cos(qJ(1));
t364 = g(1) * t319 - g(2) * t322;
t330 = -qJ(2) * t325 - t364;
t314 = sin(pkin(8));
t315 = cos(pkin(8));
t318 = sin(qJ(3));
t321 = cos(qJ(3));
t336 = t314 * t321 + t315 * t318;
t264 = t336 * qJD(1);
t320 = cos(qJ(5));
t254 = t320 * t264;
t358 = qJD(1) * t321;
t359 = qJD(1) * t318;
t267 = -t314 * t359 + t315 * t358;
t317 = sin(qJ(5));
t371 = t267 * t317;
t226 = t254 + t371;
t308 = qJD(3) + qJD(5);
t373 = t226 * t308;
t323 = -pkin(1) - pkin(6);
t281 = qJDD(1) * t323 + qJDD(2);
t274 = t321 * t281;
t282 = qJD(1) * t323 + qJD(2);
t350 = qJDD(1) * t321;
t352 = qJD(1) * qJD(4);
t353 = qJD(1) * qJD(3);
t357 = qJD(3) * t318;
t215 = -t321 * t352 - t282 * t357 + qJDD(3) * pkin(3) + t274 + (t318 * t353 - t350) * qJ(4);
t356 = qJD(3) * t321;
t220 = (-qJ(4) * qJD(1) + t282) * t356 + (-qJ(4) * qJDD(1) + t281 - t352) * t318;
t200 = t315 * t215 - t220 * t314;
t351 = qJDD(1) * t318;
t234 = t314 * t351 - t315 * t350 + t336 * t353;
t196 = qJDD(3) * pkin(4) + pkin(7) * t234 + t200;
t201 = t314 * t215 + t315 * t220;
t335 = t314 * t318 - t315 * t321;
t233 = -qJDD(1) * t336 + t335 * t353;
t197 = pkin(7) * t233 + t201;
t261 = -qJ(4) * t358 + t321 * t282;
t251 = qJD(3) * pkin(3) + t261;
t260 = -qJ(4) * t359 + t282 * t318;
t369 = t315 * t260;
t217 = t314 * t251 + t369;
t380 = pkin(7) * t264;
t206 = t217 - t380;
t276 = pkin(3) * t359 + qJD(1) * qJ(2) + qJD(4);
t239 = pkin(4) * t264 + t276;
t300 = qJ(3) + pkin(8) + qJ(5);
t291 = sin(t300);
t292 = cos(t300);
t355 = qJD(5) * t317;
t388 = g(3) * t292 - t317 * t196 - t320 * t197 + t206 * t355 + t239 * t226 + t364 * t291;
t307 = qJDD(3) + qJDD(5);
t337 = -t264 * t317 + t320 * t267;
t387 = t307 * MDP(20) + t226 * MDP(16) * t337 + (-t226 ^ 2 + t337 ^ 2) * MDP(17);
t374 = t337 * t308;
t236 = -t317 * t336 - t320 * t335;
t235 = -t317 * t335 + t320 * t336;
t265 = t314 * t357 - t315 * t356;
t266 = t336 * qJD(3);
t202 = -qJD(5) * t235 + t265 * t317 - t320 * t266;
t385 = t202 * t308 + t236 * t307;
t309 = qJDD(1) * qJ(2);
t343 = g(1) * t322 + g(2) * t319;
t310 = qJD(1) * qJD(2);
t348 = 0.2e1 * t310;
t384 = 0.2e1 * t309 + t348 - t343;
t383 = g(3) * t291 + t320 * t196 - t317 * t197 - t239 * t337 - t364 * t292;
t345 = t320 * t233 + t234 * t317;
t199 = qJD(5) * t337 - t345;
t381 = pkin(3) * t314;
t379 = pkin(7) * t267;
t378 = g(3) * t318;
t303 = t318 * pkin(3);
t377 = pkin(1) * qJDD(1);
t247 = t314 * t260;
t216 = t315 * t251 - t247;
t205 = qJD(3) * pkin(4) + t216 - t379;
t368 = t320 * t205;
t367 = qJ(2) + t303;
t366 = qJ(4) - t323;
t258 = -qJD(4) * t321 + t357 * t366;
t278 = t366 * t321;
t259 = -qJD(3) * t278 - qJD(4) * t318;
t219 = t314 * t258 + t315 * t259;
t224 = t315 * t261 - t247;
t277 = t366 * t318;
t238 = -t315 * t277 - t314 * t278;
t365 = t322 * pkin(1) + t319 * qJ(2);
t313 = t321 ^ 2;
t363 = t318 ^ 2 - t313;
t324 = qJD(3) ^ 2;
t362 = -t324 - t325;
t360 = qJD(1) * t276;
t354 = pkin(3) * t356 + qJD(2);
t349 = qJDD(3) * t318;
t347 = -qJD(5) * t254 + t317 * t233 - t320 * t234;
t346 = t321 * t353;
t218 = t315 * t258 - t259 * t314;
t223 = -t261 * t314 - t369;
t237 = t277 * t314 - t315 * t278;
t344 = qJDD(2) - t377;
t203 = t236 * qJD(5) - t320 * t265 - t266 * t317;
t341 = -t203 * t308 - t235 * t307;
t340 = -t317 * t205 - t320 * t206;
t221 = pkin(7) * t335 + t237;
t222 = -pkin(7) * t336 + t238;
t339 = t221 * t320 - t222 * t317;
t338 = t221 * t317 + t222 * t320;
t334 = qJDD(4) + t309 + t310 + (t346 + t351) * pkin(3);
t293 = pkin(3) * t315 + pkin(4);
t333 = t293 * t317 + t320 * t381;
t332 = t293 * t320 - t317 * t381;
t331 = 0.2e1 * qJ(2) * t353 + qJDD(3) * t323;
t198 = -t267 * t355 + t347;
t327 = -t200 * t335 + t201 * t336 - t216 * t266 - t217 * t265 - t364;
t326 = -t323 * t324 + t384;
t316 = -qJ(4) - pkin(6);
t302 = t322 * qJ(2);
t299 = qJDD(3) * t321;
t252 = pkin(4) * t336 + t367;
t243 = pkin(3) * t358 + pkin(4) * t267;
t240 = -pkin(4) * t265 + t354;
t211 = -pkin(4) * t233 + t334;
t210 = t224 - t379;
t209 = t223 + t380;
t208 = pkin(7) * t265 + t219;
t207 = pkin(7) * t266 + t218;
t1 = [qJDD(1) * MDP(1) + t364 * MDP(2) + t343 * MDP(3) + (qJDD(2) - t364 - 0.2e1 * t377) * MDP(4) + t384 * MDP(5) + (-t344 * pkin(1) - g(1) * (-pkin(1) * t319 + t302) - g(2) * t365 + (t348 + t309) * qJ(2)) * MDP(6) + (qJDD(1) * t313 - 0.2e1 * t318 * t346) * MDP(7) + 0.2e1 * (-t318 * t350 + t353 * t363) * MDP(8) + (-t318 * t324 + t299) * MDP(9) + (-t321 * t324 - t349) * MDP(10) + (t318 * t326 + t321 * t331) * MDP(12) + (-t318 * t331 + t321 * t326) * MDP(13) + (-t218 * t267 - t219 * t264 + t233 * t238 + t234 * t237 - t327) * MDP(14) + (t201 * t238 + t217 * t219 + t200 * t237 + t216 * t218 + t334 * t367 + t276 * t354 - g(1) * (t322 * t303 + t302 + (-pkin(1) + t316) * t319) - g(2) * (t303 * t319 - t316 * t322 + t365)) * MDP(15) + (t198 * t236 + t202 * t337) * MDP(16) + (-t198 * t235 - t199 * t236 - t202 * t226 - t203 * t337) * MDP(17) + t385 * MDP(18) + t341 * MDP(19) + (t240 * t226 + t252 * t199 + t211 * t235 + t239 * t203 + (-qJD(5) * t338 + t207 * t320 - t208 * t317) * t308 + t339 * t307 - t343 * t291) * MDP(21) + (t240 * t337 + t252 * t198 + t211 * t236 + t239 * t202 - (qJD(5) * t339 + t207 * t317 + t208 * t320) * t308 - t338 * t307 - t343 * t292) * MDP(22); qJDD(1) * MDP(4) - t325 * MDP(5) + (t344 + t330) * MDP(6) + (t318 * t362 + t299) * MDP(12) + (t321 * t362 - t349) * MDP(13) + (t233 * t336 - t234 * t335 + t264 * t265 + t266 * t267) * MDP(14) + (t327 - t360) * MDP(15) + (-qJD(1) * t226 + t385) * MDP(21) + (-qJD(1) * t337 + t341) * MDP(22); MDP(9) * t350 - MDP(10) * t351 + qJDD(3) * MDP(11) + (t321 * t330 + t274 + t378) * MDP(12) + (g(3) * t321 + (-t281 - t330) * t318) * MDP(13) + ((t217 + t223) * t267 - (t216 - t224) * t264 + (t233 * t314 + t234 * t315) * pkin(3)) * MDP(14) + (-t216 * t223 - t217 * t224 + (t378 + t200 * t315 + t201 * t314 + (-t364 - t360) * t321) * pkin(3)) * MDP(15) + (t198 + t373) * MDP(18) + (-t199 + t374) * MDP(19) + (t332 * t307 - t243 * t226 - (t209 * t320 - t210 * t317) * t308 + (-t308 * t333 + t340) * qJD(5) + t383) * MDP(21) + (-t333 * t307 - t243 * t337 + (t209 * t317 + t210 * t320) * t308 + (-t308 * t332 - t368) * qJD(5) + t388) * MDP(22) + (t321 * t318 * MDP(7) - MDP(8) * t363) * t325 + t387; (-t264 ^ 2 - t267 ^ 2) * MDP(14) + (t216 * t267 + t217 * t264 + t334 - t343) * MDP(15) + (t199 + t374) * MDP(21) + (t198 - t373) * MDP(22); (t347 + t373) * MDP(18) + (t345 + t374) * MDP(19) + (-t308 * t340 + t383) * MDP(21) + ((-t206 * t317 + t368) * t308 + t388) * MDP(22) + (-MDP(18) * t371 - t337 * MDP(19) + t340 * MDP(21) - MDP(22) * t368) * qJD(5) + t387;];
tau = t1;
