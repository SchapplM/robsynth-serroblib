% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PPRR3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:30
% EndTime: 2019-12-31 16:17:32
% DurationCPUTime: 1.93s
% Computational Cost: add. (3039->256), mult. (5744->287), div. (0->0), fcn. (3848->6), ass. (0->149)
t326 = cos(qJ(3));
t328 = qJD(3) ^ 2;
t324 = sin(qJ(3));
t349 = t324 * qJDD(3);
t293 = t326 * t328 + t349;
t318 = g(3) - qJDD(1);
t277 = pkin(4) * t293 + t326 * t318;
t320 = sin(pkin(6));
t321 = cos(pkin(6));
t348 = t326 * qJDD(3);
t294 = t324 * t328 - t348;
t338 = pkin(4) * t294 + t324 * t318;
t372 = t321 * t293 + t320 * t294;
t380 = -qJ(1) * t372 + t321 * t277 + t320 * t338;
t297 = t320 * g(1) - t321 * g(2);
t291 = -qJDD(2) + t297;
t298 = t321 * g(1) + t320 * g(2);
t256 = t326 * t291 - t324 * t298;
t257 = -t324 * t291 - t326 * t298;
t231 = t326 * t256 - t324 * t257;
t337 = t324 * t256 + t326 * t257;
t379 = t320 * t231 - t321 * t337;
t378 = t321 * t231 + t320 * t337;
t254 = -t320 * t293 + t321 * t294;
t377 = qJ(1) * t254 + t320 * t277 - t321 * t338;
t370 = pkin(1) + pkin(2);
t251 = -t328 * pkin(3) + qJDD(3) * pkin(5) + t257;
t323 = sin(qJ(4));
t325 = cos(qJ(4));
t236 = t323 * t251 - t325 * t318;
t237 = t325 * t251 + t323 * t318;
t217 = t325 * t236 - t323 * t237;
t369 = pkin(3) * t217;
t218 = t323 * t236 + t325 * t237;
t250 = -qJDD(3) * pkin(3) - t328 * pkin(5) + t256;
t212 = t326 * t218 + t324 * t250;
t368 = pkin(4) * t212;
t367 = pkin(4) * t231;
t366 = pkin(4) * t337;
t365 = pkin(5) * t326;
t363 = qJ(2) * t318;
t316 = t323 ^ 2;
t362 = t316 * t328;
t359 = t320 * t318;
t309 = t321 * t318;
t358 = t323 * t250;
t307 = t323 * t328 * t325;
t299 = qJDD(4) + t307;
t357 = t323 * t299;
t300 = qJDD(4) - t307;
t356 = t323 * t300;
t355 = t325 * t250;
t354 = t325 * t299;
t353 = t325 * t300;
t317 = t325 ^ 2;
t352 = t316 + t317;
t351 = qJD(3) * qJD(4);
t350 = t323 * qJDD(3);
t312 = t325 * qJDD(3);
t347 = t323 * t351;
t311 = t325 * t351;
t211 = t324 * t218 - t326 * t250;
t346 = -pkin(4) * t211 - t324 * t369;
t283 = t320 * t298;
t345 = t321 * t291 - t283;
t344 = t321 * t297 - t283;
t284 = t321 * t298;
t343 = -t320 * t291 - t284;
t342 = -t320 * t297 - t284;
t341 = t324 * t307;
t340 = t326 * t307;
t339 = pkin(3) * t250 - pkin(5) * t218;
t336 = -pkin(3) * t326 - pkin(5) * t324 - pkin(2);
t315 = t317 * t328;
t327 = qJD(4) ^ 2;
t304 = -t315 - t327;
t264 = t323 * t304 + t354;
t221 = -pkin(3) * t264 + t236;
t225 = -pkin(5) * t264 + t358;
t268 = t325 * t304 - t357;
t290 = t312 - 0.2e1 * t347;
t238 = t324 * t268 + t326 * t290;
t335 = -pkin(4) * t238 - t324 * t221 + t326 * t225;
t302 = -t327 - t362;
t266 = t325 * t302 - t356;
t222 = -pkin(3) * t266 + t237;
t226 = -pkin(5) * t266 + t355;
t270 = -t323 * t302 - t353;
t287 = 0.2e1 * t311 + t350;
t239 = t324 * t270 - t326 * t287;
t334 = -pkin(4) * t239 - t324 * t222 + t326 * t226;
t333 = pkin(3) * t290 + pkin(5) * t268 - t355;
t332 = pkin(3) * t287 - pkin(5) * t270 - t358;
t240 = t326 * t268 - t324 * t290;
t331 = -pkin(4) * t240 - t326 * t221 - t324 * t225;
t241 = t326 * t270 + t324 * t287;
t330 = -pkin(4) * t241 - t326 * t222 - t324 * t226;
t292 = t352 * qJDD(3);
t295 = t315 + t362;
t329 = pkin(3) * t295 + pkin(5) * t292 + t218;
t303 = t315 - t327;
t301 = t327 - t362;
t296 = -t315 + t362;
t289 = t312 - t347;
t288 = t311 + t350;
t286 = t352 * t351;
t274 = t324 * qJDD(4) + t326 * t286;
t273 = t325 * t288 - t316 * t351;
t272 = -t326 * qJDD(4) + t324 * t286;
t271 = -t323 * t289 - t317 * t351;
t269 = -t323 * t301 + t354;
t267 = t325 * t303 - t356;
t265 = t325 * t301 + t357;
t263 = t323 * t303 + t353;
t262 = (t288 + t311) * t323;
t261 = (-t289 + t347) * t325;
t260 = pkin(1) * t291 - qJ(2) * t298;
t259 = t326 * t292 - t324 * t295;
t258 = t324 * t292 + t326 * t295;
t253 = -t323 * t287 + t325 * t290;
t252 = t325 * t287 + t323 * t290;
t249 = t326 * t273 - t341;
t248 = t326 * t271 + t341;
t247 = t324 * t273 + t340;
t246 = t324 * t271 - t340;
t245 = t326 * t269 + t323 * t349;
t244 = t326 * t267 + t324 * t312;
t243 = t324 * t269 - t323 * t348;
t242 = t324 * t267 - t325 * t348;
t235 = t326 * t253 + t324 * t296;
t234 = t324 * t253 - t326 * t296;
t228 = qJ(2) * t294 + t370 * t293 + t257;
t227 = -qJ(2) * t293 + t370 * t294 + t256;
t220 = t363 + t367;
t219 = t370 * t318 - t366;
t214 = -pkin(4) * t258 + t326 * t217;
t213 = pkin(4) * t259 + t324 * t217;
t210 = qJ(2) * t337 + t231 * t370;
t209 = qJ(2) * t259 - t370 * t258 - t329;
t208 = qJ(2) * t241 - t370 * t239 + t332;
t207 = qJ(2) * t240 - t370 * t238 - t333;
t206 = qJ(2) * t266 + t334;
t205 = qJ(2) * t264 + t335;
t204 = t370 * t266 + t330;
t203 = t370 * t264 + t331;
t202 = -(qJ(2) - t365) * t217 + t346;
t201 = -t368 - (pkin(1) - t336) * t217;
t200 = qJ(2) * t212 - t370 * t211 + t339;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t359, -t309, -t344, -qJ(1) * t344, 0, 0, 0, 0, 0, 0, -t359, -t345, t309, -qJ(1) * t345 + (-t320 * pkin(1) + t321 * qJ(2)) * t318, 0, 0, -t254, 0, -t372, 0, -t377, t380, t378, -qJ(1) * t378 - t320 * t219 + t321 * t220, t320 * t247 + t321 * t249, t320 * t234 + t321 * t235, t320 * t243 + t321 * t245, t320 * t246 + t321 * t248, t320 * t242 + t321 * t244, t320 * t272 + t321 * t274, t321 * t205 - t320 * t203 - qJ(1) * (-t321 * t238 + t320 * t240), t321 * t206 - t320 * t204 - qJ(1) * (-t321 * t239 + t320 * t241), t321 * t214 + t320 * t213 - qJ(1) * (-t321 * t258 + t320 * t259), t321 * t202 - t320 * t201 - qJ(1) * (-t321 * t211 + t320 * t212); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t309, -t359, t342, qJ(1) * t342, 0, 0, 0, 0, 0, 0, t309, t343, t359, qJ(1) * t343 + (t321 * pkin(1) + t320 * qJ(2)) * t318, 0, 0, -t372, 0, t254, 0, t380, t377, t379, -qJ(1) * t379 + t321 * t219 + t320 * t220, -t321 * t247 + t320 * t249, -t321 * t234 + t320 * t235, -t321 * t243 + t320 * t245, -t321 * t246 + t320 * t248, -t321 * t242 + t320 * t244, -t321 * t272 + t320 * t274, t320 * t205 + t321 * t203 + qJ(1) * (t320 * t238 + t321 * t240), t320 * t206 + t321 * t204 + qJ(1) * (t320 * t239 + t321 * t241), t320 * t214 - t321 * t213 + qJ(1) * (t320 * t258 + t321 * t259), t320 * t202 + t321 * t201 + qJ(1) * (t320 * t211 + t321 * t212); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t297, t298, 0, 0, 0, 0, 0, 0, 0, 0, t291, 0, -t298, t260, 0, 0, 0, 0, 0, -qJDD(3), t227, t228, 0, t210, -t262, -t252, -t265, t261, -t263, 0, t207, t208, t209, t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t318, -t297, 0, 0, 0, 0, 0, 0, 0, 0, -t291, t318, t363, 0, 0, -t294, 0, -t293, 0, t338, t277, t231, t220, t249, t235, t245, t248, t244, t274, t205, t206, t214, t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t318, 0, -t298, 0, 0, 0, 0, 0, 0, 0, t318, -t298, 0, pkin(1) * t318, 0, 0, -t293, 0, t294, 0, t277, -t338, -t337, t219, -t247, -t234, -t243, -t246, -t242, -t272, t203, t204, -t213, t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, t298, 0, 0, 0, 0, 0, 0, 0, 0, t291, 0, -t298, t260, 0, 0, 0, 0, 0, -qJDD(3), t227, t228, 0, t210, -t262, -t252, -t265, t261, -t263, 0, t207, t208, t209, t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t291, t318, 0, 0, 0, -t294, 0, -t293, 0, t338, t277, t231, t367, t249, t235, t245, t248, t244, t274, t335, t334, t214, t217 * t365 + t346; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, 0, -t298, 0, 0, 0, 0, 0, 0, -qJDD(3), pkin(2) * t294 + t256, pkin(2) * t293 + t257, 0, pkin(2) * t231, -t262, -t252, -t265, t261, -t263, 0, -pkin(2) * t238 - t333, -pkin(2) * t239 + t332, -pkin(2) * t258 - t329, -pkin(2) * t211 + t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t318, t298, 0, 0, 0, 0, t293, 0, -t294, 0, -t277, t338, t337, -pkin(2) * t318 + t366, t247, t234, t243, t246, t242, t272, -pkin(2) * t264 - t331, -pkin(2) * t266 - t330, t213, -t217 * t336 + t368; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t328, 0, 0, t318, t256, 0, t273, t253, t269, t271, t267, t286, t225, t226, t217, pkin(5) * t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t328, 0, qJDD(3), 0, -t318, 0, t257, 0, t307, -t296, -t350, -t307, -t312, -qJDD(4), t221, t222, 0, t369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t256, -t257, 0, 0, t262, t252, t265, -t261, t263, 0, t333, -t332, t329, -t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, t290, t299, -t311, t303, t311, 0, t250, t236, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t347, t287, t301, t289, t300, -t347, -t250, 0, t237, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t307, t296, t350, t307, t312, qJDD(4), -t236, -t237, 0, 0;];
m_new_reg = t1;