% Calculate vector of inverse dynamics joint torques for
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:47
% EndTime: 2019-12-31 19:49:50
% DurationCPUTime: 1.43s
% Computational Cost: add. (1411->250), mult. (2234->289), div. (0->0), fcn. (1233->12), ass. (0->138)
t286 = cos(qJ(2));
t364 = pkin(1) * t286;
t265 = qJDD(1) * t364;
t275 = qJDD(1) + qJDD(2);
t283 = sin(qJ(2));
t354 = pkin(1) * qJD(2);
t323 = qJD(1) * t354;
t219 = pkin(2) * t275 - t283 * t323 + t265;
t280 = sin(pkin(8));
t281 = cos(pkin(8));
t329 = qJDD(1) * t283;
t371 = pkin(1) * t329 + t286 * t323;
t201 = t280 * t219 + t281 * t371;
t372 = pkin(7) * t275 + qJD(3) * qJD(4) + t201;
t285 = cos(qJ(4));
t276 = qJD(1) + qJD(2);
t355 = pkin(1) * qJD(1);
t324 = t286 * t355;
t232 = pkin(2) * t276 + t324;
t325 = t283 * t355;
t248 = t281 * t325;
t215 = t232 * t280 + t248;
t211 = pkin(7) * t276 + t215;
t282 = sin(qJ(4));
t352 = t211 * t282;
t203 = qJD(3) * t285 - t352;
t370 = qJD(5) - t203;
t277 = t282 ^ 2;
t278 = t285 ^ 2;
t336 = t277 + t278;
t369 = t276 * t336;
t320 = t282 * qJDD(3) + t285 * t372;
t328 = qJDD(4) * qJ(5);
t188 = t328 + (qJD(5) - t352) * qJD(4) + t320;
t332 = qJD(4) * t285;
t309 = -t285 * qJDD(3) + t211 * t332 + t282 * t372;
t353 = (qJDD(4) * pkin(4));
t365 = qJDD(5) - t353;
t189 = t309 + t365;
t368 = t188 * t285 + t189 * t282;
t199 = -qJD(4) * pkin(4) + t370;
t204 = qJD(3) * t282 + t211 * t285;
t202 = qJD(4) * qJ(5) + t204;
t279 = qJ(1) + qJ(2);
t271 = sin(t279);
t272 = cos(t279);
t367 = g(1) * t271 - g(2) * t272;
t270 = pkin(8) + t279;
t254 = sin(t270);
t255 = cos(t270);
t366 = g(1) * t255 + g(2) * t254;
t363 = pkin(2) * t271;
t362 = pkin(2) * t281;
t361 = pkin(4) * t285;
t360 = g(1) * t254;
t357 = g(2) * t255;
t264 = pkin(2) + t364;
t343 = t281 * t283;
t339 = pkin(1) * t343 + t264 * t280;
t221 = pkin(7) + t339;
t288 = qJD(4) ^ 2;
t351 = t221 * t288;
t222 = t280 * t324 + t248;
t350 = t222 * t276;
t349 = t254 * t282;
t348 = t255 * t282;
t256 = pkin(2) * t280 + pkin(7);
t347 = t256 * t288;
t346 = t275 * t282;
t345 = t276 * t282;
t344 = t280 * t283;
t342 = t282 * t285;
t331 = qJD(5) * t282;
t333 = qJD(4) * t282;
t226 = pkin(4) * t333 - qJ(5) * t332 - t331;
t341 = t226 - t222;
t340 = g(1) * t349 - g(2) * t348;
t338 = g(1) * t272 + g(2) * t271;
t337 = t277 - t278;
t334 = qJD(4) * t276;
t327 = qJDD(4) * t221;
t326 = qJDD(4) * t256;
t274 = t276 ^ 2;
t321 = t274 * t342;
t247 = t280 * t325;
t224 = t281 * t324 - t247;
t240 = t285 * t360;
t319 = t224 * t333 + t285 * t350 + t240;
t200 = t219 * t281 - t280 * t371;
t195 = -pkin(3) * t275 - t200;
t318 = -t195 - t357;
t317 = t336 * t275;
t315 = t203 + t352;
t304 = qJ(5) * t282 + t361;
t298 = -pkin(3) - t304;
t313 = -pkin(1) * t344 + t264 * t281;
t212 = t298 - t313;
t225 = (t281 * t286 - t344) * t354;
t314 = t212 * t276 - t225;
t214 = t232 * t281 - t247;
t312 = qJD(1) * (-qJD(2) + t276);
t311 = qJD(2) * (-qJD(1) - t276);
t210 = -pkin(3) * t276 - t214;
t308 = t195 * t282 + t210 * t332 - t340;
t307 = t265 + t367;
t284 = sin(qJ(1));
t306 = -pkin(1) * t284 - t363;
t303 = pkin(4) * t282 - qJ(5) * t285;
t234 = qJDD(4) * t282 + t285 * t288;
t235 = qJDD(4) * t285 - t282 * t288;
t302 = 0.2e1 * (t275 * t342 - t334 * t337) * MDP(9) + (t275 * t277 + 0.2e1 * t332 * t345) * MDP(8) + t234 * MDP(10) + t235 * MDP(11) + t275 * MDP(4);
t263 = pkin(2) * t272;
t301 = t254 * pkin(7) + qJ(5) * t348 + t263 + (pkin(3) + t361) * t255;
t300 = t199 * t282 + t202 * t285;
t257 = -pkin(3) - t362;
t299 = t257 * t275 + t347;
t190 = (qJD(4) * t303 - t331) * t276 + t298 * t275 - t200;
t229 = t298 - t362;
t297 = -t229 * t275 - t190 - t347;
t223 = (t280 * t286 + t343) * t354;
t220 = -pkin(3) - t313;
t296 = t220 * t275 + t223 * t276 + t351;
t295 = t298 * t360;
t294 = g(1) * t348 + g(2) * t349 - g(3) * t285 - t309;
t293 = -t327 + (t220 * t276 - t225) * qJD(4);
t205 = t223 + t226;
t292 = -t205 * t276 - t212 * t275 - t190 - t351;
t291 = qJD(4) * t204 + t294;
t290 = t199 * t332 - t202 * t333 - t366 + t368;
t289 = (t199 * t285 - t202 * t282) * qJD(4) + t368;
t287 = cos(qJ(1));
t273 = t287 * pkin(1);
t250 = t255 * pkin(7);
t227 = t303 * t276;
t206 = t210 * t333;
t198 = t276 * t298 - t214;
t194 = t198 * t333;
t1 = [qJDD(1) * MDP(1) + (g(1) * t284 - g(2) * t287) * MDP(2) + (g(1) * t287 + g(2) * t284) * MDP(3) + ((t275 * t286 + t283 * t311) * pkin(1) + t307) * MDP(5) + (((-qJDD(1) - t275) * t283 + t286 * t311) * pkin(1) + t338) * MDP(6) + (t201 * t339 + t215 * t225 + t200 * t313 - t214 * t223 - g(1) * t306 - g(2) * (t263 + t273)) * MDP(7) + (t206 + t240 + t293 * t282 + (-t296 + t318) * t285) * MDP(13) + (t282 * t296 + t285 * t293 + t308) * MDP(14) + (t194 + t240 + (qJD(4) * t314 - t327) * t282 + (t292 - t357) * t285) * MDP(15) + (t221 * t317 + t225 * t369 + t290) * MDP(16) + ((t327 + (-t198 - t314) * qJD(4)) * t285 + t292 * t282 + t340) * MDP(17) + (t190 * t212 + t198 * t205 - g(1) * (t250 + t306) - g(2) * (t273 + t301) - t295 + t300 * t225 + t289 * t221) * MDP(18) + t302; (pkin(1) * t283 * t312 + t307) * MDP(5) + ((t286 * t312 - t329) * pkin(1) + t338) * MDP(6) + (t214 * t222 - t215 * t224 + (t200 * t281 + t201 * t280 + t367) * pkin(2)) * MDP(7) + (t206 + (t257 * t334 - t326) * t282 + (-t299 + t318) * t285 + t319) * MDP(13) + ((-t326 + (t257 * t276 + t224) * qJD(4)) * t285 + (t299 - t350) * t282 + t308) * MDP(14) + (t194 + (t229 * t334 - t326) * t282 + (-t226 * t276 + t297 - t357) * t285 + t319) * MDP(15) + (-t224 * t369 + t256 * t317 + t290) * MDP(16) + ((t326 + (-t229 * t276 - t198 - t224) * qJD(4)) * t285 + (-t276 * t341 + t297) * t282 + t340) * MDP(17) + (t190 * t229 - g(1) * (t250 - t363) - g(2) * t301 - t295 - t300 * t224 + t341 * t198 + t289 * t256) * MDP(18) + t302; (qJDD(3) - g(3)) * MDP(7) + (qJD(4) * t300 + t188 * t282 - t189 * t285 - g(3)) * MDP(18) + (MDP(13) + MDP(15)) * t235 + (-MDP(14) + MDP(17)) * t234; -MDP(8) * t321 + t337 * MDP(9) * t274 + MDP(10) * t346 + qJDD(4) * MDP(12) + (-t210 * t345 + t291) * MDP(13) + (g(3) * t282 + t315 * qJD(4) + (-t210 * t276 + t366) * t285 - t320) * MDP(14) + ((2 * t353) - qJDD(5) + (-t198 * t282 + t227 * t285) * t276 + t291) * MDP(15) + (0.2e1 * t328 + (t227 * t276 - g(3)) * t282 + (t198 * t276 - t366) * t285 + (0.2e1 * qJD(5) - t315) * qJD(4) + t320) * MDP(17) + (-t189 * pkin(4) - g(3) * t304 + t188 * qJ(5) - t198 * t227 - t199 * t204 + t202 * t370 + t366 * t303) * MDP(18) + (MDP(11) * t285 - MDP(16) * t303) * t275; (-qJDD(4) - t321) * MDP(15) + MDP(16) * t346 + (-t274 * t277 - t288) * MDP(17) + (-qJD(4) * t202 + t198 * t345 - t294 + t365) * MDP(18);];
tau = t1;
