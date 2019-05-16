% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RPPR2
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RPPR2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:11:20
% EndTime: 2019-05-04 19:11:22
% DurationCPUTime: 1.61s
% Computational Cost: add. (6083->178), mult. (9275->170), div. (0->0), fcn. (3832->6), ass. (0->96)
t280 = -qJD(1) + qJD(4);
t278 = t280 ^ 2;
t279 = qJDD(1) - qJDD(4);
t285 = sin(qJ(4));
t287 = cos(qJ(4));
t260 = t278 * t287 - t279 * t285;
t262 = t278 * t285 + t279 * t287;
t283 = sin(pkin(6));
t284 = cos(pkin(6));
t232 = t260 * t284 - t262 * t283;
t282 = g(3) + qJDD(3);
t245 = pkin(5) * t260 + t282 * t287;
t334 = pkin(5) * t262 + t282 * t285;
t219 = qJ(3) * t232 + t245 * t284 - t283 * t334;
t286 = sin(qJ(1));
t288 = cos(qJ(1));
t236 = t260 * t283 + t262 * t284;
t331 = t232 * t286 - t236 * t288;
t339 = qJ(3) * t236 + t245 * t283 + t284 * t334;
t349 = -pkin(4) * t331 + t286 * t219 - t288 * t339;
t340 = t232 * t288 + t236 * t286;
t348 = -pkin(4) * t340 + t288 * t219 + t286 * t339;
t306 = (qJD(2) * qJD(1));
t276 = 2 * t306;
t289 = qJD(1) ^ 2;
t272 = t288 * g(1) + t286 * g(2);
t281 = qJDD(1) * qJ(2);
t294 = -t272 + t281;
t319 = pkin(1) + pkin(2);
t248 = -t319 * t289 + t276 + t294;
t271 = t286 * g(1) - t288 * g(2);
t300 = -qJDD(2) + t271;
t292 = -t289 * qJ(2) - t300;
t253 = -t319 * qJDD(1) + t292;
t229 = t283 * t248 - t284 * t253;
t227 = -qJDD(1) * pkin(3) - t229;
t230 = t284 * t248 + t283 * t253;
t228 = -pkin(3) * t289 + t230;
t210 = -t287 * t227 + t285 * t228;
t211 = t285 * t227 + t287 * t228;
t304 = t210 * t285 + t287 * t211;
t205 = t210 * t287 - t211 * t285;
t311 = t205 * t284;
t197 = -t283 * t304 + t311;
t312 = t205 * t283;
t336 = t284 * t304 + t312;
t346 = t197 * t286 - t288 * t336;
t345 = t197 * t288 + t286 * t336;
t267 = -t283 * qJDD(1) + t284 * t289;
t251 = qJ(3) * t267 + t282 * t284;
t268 = qJDD(1) * t284 + t283 * t289;
t299 = qJ(3) * t268 + t282 * t283;
t321 = t288 * t267 + t268 * t286;
t335 = -pkin(4) * t321 + t288 * t251 + t286 * t299;
t214 = t229 * t284 - t230 * t283;
t298 = t229 * t283 + t230 * t284;
t333 = t214 * t286 - t288 * t298;
t332 = t214 * t288 + t286 * t298;
t238 = -t267 * t286 + t288 * t268;
t330 = pkin(4) * t238 + t286 * t251 - t288 * t299;
t318 = pkin(1) * t282;
t317 = pkin(3) * t205;
t315 = qJ(2) * t282;
t314 = qJ(3) * t214;
t313 = qJDD(1) * pkin(1);
t275 = pkin(2) * t282;
t305 = qJ(3) * t298 - t275;
t254 = pkin(1) * t289 - t294 - (2 * t306);
t255 = -t292 + t313;
t303 = -t288 * t254 - t255 * t286;
t302 = -t271 * t286 - t288 * t272;
t301 = -pkin(3) * t260 - t211;
t269 = qJDD(1) * t286 + t288 * t289;
t257 = -pkin(4) * t269 + g(3) * t288;
t270 = qJDD(1) * t288 - t286 * t289;
t256 = pkin(4) * t270 + g(3) * t286;
t297 = t254 * t286 - t255 * t288;
t295 = t271 * t288 - t272 * t286;
t293 = -pkin(3) * t262 - t210;
t202 = -pkin(3) * t282 + pkin(5) * t304;
t291 = pkin(5) * t311 + qJ(3) * t197 - t202 * t283;
t290 = pkin(5) * t312 + qJ(3) * t336 + t202 * t284 - t275;
t265 = t300 + 0.2e1 * t313;
t258 = -t272 + t276 + 0.2e1 * t281;
t231 = pkin(1) * t255 - qJ(2) * t254;
t221 = -qJ(2) * t267 + t319 * t268 + t229;
t220 = qJ(2) * t268 + t319 * t267 + t230;
t209 = t314 + t315;
t207 = -t305 + t318;
t201 = -qJ(2) * t232 + t236 * t319 - t293;
t200 = qJ(2) * t236 + t232 * t319 - t301;
t199 = qJ(2) * t298 + t214 * t319;
t194 = t291 + t315;
t193 = -t290 + t318;
t192 = qJ(2) * t336 + t197 * t319 + t317;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t270, 0, -t269, 0, -t256, -t257, -t295, -pkin(4) * t295, 0, t270, 0, 0, t269, 0, -t256, t297, t257, pkin(4) * t297 + (-pkin(1) * t286 + qJ(2) * t288) * g(3), 0, 0, -t238, 0, -t321, 0, -t330, t335, t332, -pkin(4) * t332 - t286 * t207 + t288 * t209, 0, 0, t331, 0, -t340, 0, -t349, t348, t345, -pkin(4) * t345 - t286 * t193 + t288 * t194; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t269, 0, t270, 0, t257, -t256, t302, pkin(4) * t302, 0, t269, 0, 0, -t270, 0, t257, t303, t256, pkin(4) * t303 + (pkin(1) * t288 + qJ(2) * t286) * g(3), 0, 0, -t321, 0, t238, 0, t335, t330, t333, -pkin(4) * t333 + t288 * t207 + t286 * t209, 0, 0, -t340, 0, -t331, 0, t348, t349, t346, -pkin(4) * t346 + t288 * t193 + t286 * t194; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t271, t272, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t265, 0, t258, t231, 0, 0, 0, 0, 0, qJDD(1), t221, t220, 0, t199, 0, 0, 0, 0, 0, t279, t201, t200, 0, t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t289, 0, 0, -g(3), -t271, 0, 0, qJDD(1), 0, 0, t289, 0, 0, -t255, g(3), qJ(2) * g(3), 0, 0, -t268, 0, -t267, 0, t299, t251, t214, t209, 0, 0, -t236, 0, -t232, 0, t339, t219, t197, t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, 0, qJDD(1), 0, g(3), 0, -t272, 0, 0, t289, 0, 0, -qJDD(1), 0, g(3), -t254, 0, pkin(1) * g(3), 0, 0, -t267, 0, t268, 0, t251, -t299, -t298, t207, 0, 0, -t232, 0, t236, 0, t219, -t339, -t336, t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t271, t272, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t265, 0, t258, t231, 0, 0, 0, 0, 0, qJDD(1), t221, t220, 0, t199, 0, 0, 0, 0, 0, t279, t201, t200, 0, t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t289, 0, 0, -t255, g(3), 0, 0, 0, -t268, 0, -t267, 0, t299, t251, t214, t314, 0, 0, -t236, 0, -t232, 0, t339, t219, t197, t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t255, 0, -t254, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(2) * t268 + t229, pkin(2) * t267 + t230, 0, pkin(2) * t214, 0, 0, 0, 0, 0, t279, pkin(2) * t236 - t293, pkin(2) * t232 - t301, 0, pkin(2) * t197 + t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t289, 0, 0, qJDD(1), 0, -g(3), t254, 0, 0, 0, 0, t267, 0, -t268, 0, -t251, t299, t298, t305, 0, 0, t232, 0, -t236, 0, -t219, t339, t336, t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t289, 0, 0, t282, t229, 0, 0, 0, -t262, 0, -t260, 0, t334, t245, t205, pkin(5) * t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, 0, -qJDD(1), 0, -t282, 0, t230, 0, 0, 0, t260, 0, -t262, 0, -t245, t334, t304, t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), -t229, -t230, 0, 0, 0, 0, 0, 0, 0, -t279, t293, t301, 0, -t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t279, 0, -t278, 0, 0, t282, t210, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, 0, -t279, 0, -t282, 0, t211, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t279, -t210, -t211, 0, 0;];
m_new_reg  = t1;
