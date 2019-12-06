% Calculate vector of inverse dynamics joint torques for
% S5PRPRR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:58:02
% EndTime: 2019-12-05 15:58:06
% DurationCPUTime: 2.81s
% Computational Cost: add. (1576->319), mult. (3775->435), div. (0->0), fcn. (3162->14), ass. (0->153)
t294 = sin(pkin(10));
t297 = cos(pkin(10));
t366 = t297 * MDP(5);
t393 = -t294 * MDP(6) + t366;
t300 = sin(qJ(4));
t303 = cos(qJ(4));
t266 = t294 * t303 + t297 * t300;
t263 = t266 * qJD(4);
t349 = qJDD(2) * t303;
t281 = t297 * t349;
t350 = qJDD(2) * t300;
t329 = -t294 * t350 + t281;
t229 = qJD(2) * t263 - t329;
t224 = qJDD(5) + t229;
t265 = t294 * t300 - t303 * t297;
t284 = -pkin(3) * t297 - pkin(2);
t227 = pkin(4) * t265 - pkin(8) * t266 + t284;
t293 = pkin(10) + qJ(4);
t287 = cos(t293);
t298 = cos(pkin(5));
t304 = cos(qJ(2));
t385 = cos(pkin(9));
t338 = t385 * t304;
t295 = sin(pkin(9));
t301 = sin(qJ(2));
t370 = t295 * t301;
t255 = -t298 * t338 + t370;
t339 = t385 * t301;
t369 = t295 * t304;
t257 = t298 * t369 + t339;
t333 = g(1) * t257 + g(2) * t255;
t296 = sin(pkin(5));
t367 = t296 * t304;
t312 = g(3) * t367 - t333;
t308 = t312 * t287;
t392 = t227 * t224 - t308;
t357 = qJD(2) * t303;
t282 = t297 * t357;
t358 = qJD(2) * t300;
t343 = t294 * t358;
t259 = -t282 + t343;
t252 = qJD(5) + t259;
t362 = t294 ^ 2 + t297 ^ 2;
t391 = t362 * MDP(7);
t302 = cos(qJ(5));
t219 = t302 * t224;
t299 = sin(qJ(5));
t356 = qJD(5) * t299;
t390 = t252 * t356 - t219;
t261 = t266 * qJD(2);
t344 = qJD(1) * t367;
t330 = qJD(3) - t344;
t360 = qJD(1) * t301;
t345 = t296 * t360;
t268 = qJD(2) * qJ(3) + t345;
t361 = qJD(1) * t298;
t278 = t297 * t361;
t234 = t278 + (-pkin(7) * qJD(2) - t268) * t294;
t247 = t297 * t268 + t294 * t361;
t359 = qJD(2) * t297;
t235 = pkin(7) * t359 + t247;
t211 = t234 * t300 + t235 * t303;
t351 = qJDD(2) * qJ(3);
t353 = qJDD(1) * t296;
t245 = t301 * t353 + t351 + (qJD(3) + t344) * qJD(2);
t352 = qJDD(1) * t298;
t276 = t297 * t352;
t384 = pkin(7) * qJDD(2);
t215 = t276 + (-t245 - t384) * t294;
t222 = t297 * t245 + t294 * t352;
t216 = t297 * t384 + t222;
t327 = -t215 * t303 + t216 * t300;
t199 = -qJDD(4) * pkin(4) + qJD(4) * t211 + t327;
t256 = t298 * t339 + t369;
t258 = -t298 * t370 + t338;
t286 = sin(t293);
t340 = t296 * t385;
t368 = t296 * t301;
t371 = t295 * t296;
t316 = g(1) * (-t258 * t286 + t287 * t371) + g(2) * (-t256 * t286 - t287 * t340) + g(3) * (-t286 * t368 + t287 * t298);
t389 = (pkin(4) * t261 + t252 * pkin(8)) * t252 + t199 + t316;
t210 = t234 * t303 - t235 * t300;
t326 = t215 * t300 + t216 * t303;
t198 = qJDD(4) * pkin(8) + qJD(4) * t210 + t326;
t206 = -qJD(4) * pkin(4) - t210;
t251 = t284 * qJD(2) + t330;
t212 = pkin(4) * t259 - pkin(8) * t261 + t251;
t386 = pkin(7) + qJ(3);
t269 = t386 * t294;
t270 = t386 * t297;
t239 = -t269 * t300 + t270 * t303;
t262 = t265 * qJD(4);
t332 = g(1) * t258 + g(2) * t256;
t348 = g(3) * t368;
t313 = t265 * t367;
t323 = -t269 * t303 - t270 * t300;
t364 = -qJD(1) * t313 + t265 * qJD(3) - t323 * qJD(4);
t388 = -(qJD(5) * t212 + t198) * t265 + t199 * t266 - t206 * t262 + (-qJD(5) * t227 + t364) * t252 - t239 * t224 - t348 - t332;
t383 = qJDD(2) * pkin(2);
t346 = qJD(4) * t282 + t294 * t349 + t297 * t350;
t228 = -qJD(4) * t343 + t346;
t354 = t302 * qJD(4);
t347 = qJD(5) * t354 + t299 * qJDD(4) + t302 * t228;
t204 = -t261 * t356 + t347;
t382 = t204 * t299;
t381 = t206 * t266;
t380 = t224 * t299;
t373 = t261 * t299;
t242 = -t354 + t373;
t378 = t242 * t252;
t377 = t242 * t261;
t244 = qJD(4) * t299 + t261 * t302;
t376 = t244 * t252;
t375 = t244 * t261;
t374 = (-t268 * t294 + t278) * t294;
t365 = qJDD(1) - g(3);
t314 = t266 * t367;
t363 = -qJD(1) * t314 + t266 * qJD(3) + t239 * qJD(4);
t355 = qJD(5) * t302;
t341 = qJD(2) * t360;
t336 = -t302 * qJDD(4) + t228 * t299;
t335 = t252 * t302;
t331 = pkin(4) * t263 + pkin(8) * t262 - t345;
t207 = qJD(4) * pkin(8) + t211;
t201 = t207 * t302 + t212 * t299;
t328 = t207 * t299 - t212 * t302;
t325 = -t247 * t297 + t374;
t253 = -t294 * t368 + t297 * t298;
t254 = t294 * t298 + t297 * t368;
t324 = t253 * t303 - t254 * t300;
t218 = t253 * t300 + t254 * t303;
t318 = t296 * t341 - t304 * t353 + qJDD(3);
t248 = t318 - t383;
t322 = -t248 + t333;
t321 = -t299 * t259 * t252 - t390;
t319 = MDP(3) + t393;
t317 = -t262 * t302 - t266 * t356;
t315 = -qJD(2) * t374 + t247 * t359;
t309 = -pkin(8) * t224 + (t206 + t210) * t252;
t221 = -t245 * t294 + t276;
t307 = -t221 * t294 + t222 * t297 - t332;
t236 = t284 * qJDD(2) + t318;
t305 = qJD(2) ^ 2;
t267 = -qJD(2) * pkin(2) + t330;
t250 = t286 * t298 + t287 * t368;
t233 = t258 * t287 + t286 * t371;
t231 = t256 * t287 - t286 * t340;
t209 = qJD(2) * t314 + t218 * qJD(4);
t208 = -qJD(2) * t313 + t324 * qJD(4);
t205 = qJD(5) * t244 + t336;
t203 = pkin(4) * t229 - pkin(8) * t228 + t236;
t202 = t302 * t203;
t1 = [t365 * MDP(1) + (t221 * t253 + t222 * t254 - g(3)) * MDP(8) + (-qJD(4) * t209 + qJDD(4) * t324) * MDP(14) + (-qJD(4) * t208 - qJDD(4) * t218) * MDP(15) + ((-t208 * t299 - t218 * t355) * t252 - t218 * t380 + t209 * t242 - t324 * t205) * MDP(21) + (-(t208 * t302 - t218 * t356) * t252 - t218 * t219 + t209 * t244 - t324 * t204) * MDP(22) + (-t253 * t294 + t254 * t297) * MDP(7) * qJDD(2) + ((-qJDD(2) * MDP(4) - t319 * t305 + (MDP(14) * t259 + MDP(15) * t261 + MDP(8) * t267 + (MDP(21) * t302 - MDP(22) * t299) * t252) * qJD(2)) * t301 + ((-t248 + t315) * MDP(8) - t229 * MDP(14) - t228 * MDP(15) + t390 * MDP(21) + (t252 * t355 + t380) * MDP(22) + (-MDP(4) + t391) * t305 + t319 * qJDD(2)) * t304) * t296; qJDD(2) * MDP(2) + (t365 * t367 + t333) * MDP(3) + (-t365 * t368 + t332) * MDP(4) + (-t348 + t307 + (qJD(2) * t330 + t351) * t362) * MDP(7) + (-t325 * qJD(3) + t322 * pkin(2) + t307 * qJ(3) + (-g(3) * (pkin(2) * t304 + qJ(3) * t301) + (-t267 * t301 + t325 * t304) * qJD(1)) * t296) * MDP(8) + (t228 * t266 - t261 * t262) * MDP(9) + (-t228 * t265 - t229 * t266 + t259 * t262 - t261 * t263) * MDP(10) + (-qJD(4) * t262 + qJDD(4) * t266) * MDP(11) + (-qJD(4) * t263 - qJDD(4) * t265) * MDP(12) + (-t363 * qJD(4) + qJDD(4) * t323 + t229 * t284 + t236 * t265 + t251 * t263 - t259 * t345 - t308) * MDP(14) + (t364 * qJD(4) - qJDD(4) * t239 + t228 * t284 + t236 * t266 - t251 * t262 - t261 * t345 + t312 * t286) * MDP(15) + (t204 * t266 * t302 + t317 * t244) * MDP(16) + (-(-t242 * t302 - t244 * t299) * t262 + (-t382 - t205 * t302 + (t242 * t299 - t244 * t302) * qJD(5)) * t266) * MDP(17) + (t204 * t265 + t266 * t219 + t244 * t263 + t317 * t252) * MDP(18) + (-t266 * t380 - t205 * t265 - t242 * t263 + (t262 * t299 - t266 * t355) * t252) * MDP(19) + (t224 * t265 + t252 * t263) * MDP(20) + (-t328 * t263 + t202 * t265 - t323 * t205 + t363 * t242 + (t331 * t252 + (-t207 * t265 - t239 * t252 + t381) * qJD(5) + t392) * t302 + t388 * t299) * MDP(21) + (-t201 * t263 - t323 * t204 + t363 * t244 + (-(-qJD(5) * t207 + t203) * t265 - qJD(5) * t381 + (qJD(5) * t239 - t331) * t252 - t392) * t299 + t388 * t302) * MDP(22) + t393 * ((-g(3) * t304 + t341) * t296 + t322 + t383); (t312 - t315 + t318) * MDP(8) - t281 * MDP(14) + t346 * MDP(15) + (t321 - t377) * MDP(21) + (-t252 ^ 2 * t302 - t375 - t380) * MDP(22) - t305 * t391 + (-t366 - pkin(2) * MDP(8) + (MDP(14) * t300 + MDP(6)) * t294) * qJDD(2) + ((t294 * t357 + t297 * t358 + t261) * MDP(14) + (-t259 - t343) * MDP(15)) * qJD(4); -t259 ^ 2 * MDP(10) + ((t259 - t343) * qJD(4) + t346) * MDP(11) + t329 * MDP(12) + qJDD(4) * MDP(13) + (-t316 - t327) * MDP(14) + (g(1) * t233 + g(2) * t231 + g(3) * t250 + t251 * t259 - t326) * MDP(15) + (t244 * t335 + t382) * MDP(16) + ((t204 - t378) * t302 + (-t205 - t376) * t299) * MDP(17) + (t252 * t335 - t375 + t380) * MDP(18) + (t321 + t377) * MDP(19) + (-pkin(4) * t205 - t211 * t242 + t309 * t299 - t302 * t389) * MDP(21) + (-pkin(4) * t204 - t211 * t244 + t299 * t389 + t309 * t302) * MDP(22) + (MDP(10) * t261 - t251 * MDP(14) - t252 * MDP(20) + MDP(21) * t328 + t201 * MDP(22) + t259 * MDP(9)) * t261; t244 * t242 * MDP(16) + (-t242 ^ 2 + t244 ^ 2) * MDP(17) + (t347 + t378) * MDP(18) + (-t336 + t376) * MDP(19) + t224 * MDP(20) + (-t299 * t198 + t202 + t201 * t252 - t206 * t244 - g(1) * (-t233 * t299 + t257 * t302) - g(2) * (-t231 * t299 + t255 * t302) - g(3) * (-t250 * t299 - t302 * t367)) * MDP(21) + (-t302 * t198 - t299 * t203 - t328 * t252 + t206 * t242 - g(1) * (-t233 * t302 - t257 * t299) - g(2) * (-t231 * t302 - t255 * t299) - g(3) * (-t250 * t302 + t299 * t367)) * MDP(22) + (-MDP(18) * t373 - t244 * MDP(19) - t201 * MDP(21) + t328 * MDP(22)) * qJD(5);];
tau = t1;
