% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PRPR2
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
%   pkin=[a2,a3,a4,d2,d4,theta3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:00
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PRPR2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:59:57
% EndTime: 2019-05-04 18:59:58
% DurationCPUTime: 1.66s
% Computational Cost: add. (5738->158), mult. (7775->154), div. (0->0), fcn. (5022->6), ass. (0->95)
t297 = qJD(2) + qJD(4);
t295 = t297 ^ 2;
t296 = qJDD(2) + qJDD(4);
t302 = sin(qJ(4));
t304 = cos(qJ(4));
t273 = t304 * t295 + t302 * t296;
t276 = t302 * t295 - t304 * t296;
t300 = sin(pkin(6));
t301 = cos(pkin(6));
t237 = t301 * t273 - t300 * t276;
t298 = g(3) - qJDD(3);
t259 = pkin(5) * t273 - t304 * t298;
t354 = pkin(5) * t276 - t302 * t298;
t211 = qJ(3) * t237 + t301 * t259 - t300 * t354;
t241 = t300 * t273 + t301 * t276;
t303 = sin(qJ(2));
t305 = cos(qJ(2));
t217 = t305 * t237 - t303 * t241;
t360 = qJ(3) * t241 + t300 * t259 + t301 * t354;
t197 = pkin(4) * t217 + t305 * t211 - t303 * t360;
t221 = t303 * t237 + t305 * t241;
t366 = pkin(4) * t221 + t303 * t211 + t305 * t360;
t299 = g(2) - qJDD(1);
t284 = t303 * g(1) - t305 * t299;
t278 = qJDD(2) * pkin(2) + t284;
t285 = t305 * g(1) + t303 * t299;
t306 = qJD(2) ^ 2;
t279 = -t306 * pkin(2) - t285;
t243 = -t301 * t278 + t300 * t279;
t233 = qJDD(2) * pkin(3) - t243;
t244 = t300 * t278 + t301 * t279;
t234 = -t306 * pkin(3) + t244;
t214 = -t304 * t233 + t302 * t234;
t215 = t302 * t233 + t304 * t234;
t328 = t302 * t214 + t304 * t215;
t202 = t304 * t214 - t302 * t215;
t338 = t301 * t202;
t191 = -t300 * t328 + t338;
t339 = t300 * t202;
t356 = t301 * t328 + t339;
t363 = t303 * t191 + t305 * t356;
t186 = t305 * t191 - t303 * t356;
t327 = t300 * t243 + t301 * t244;
t225 = t301 * t243 - t300 * t244;
t332 = t305 * t225;
t206 = -t303 * t327 + t332;
t336 = t303 * t225;
t357 = t305 * t327 + t336;
t286 = t300 * qJDD(2) + t301 * t306;
t263 = qJ(3) * t286 - t301 * t298;
t287 = t301 * qJDD(2) - t300 * t306;
t318 = -qJ(3) * t287 - t300 * t298;
t324 = -t303 * t286 + t305 * t287;
t355 = -pkin(4) * t324 + t303 * t263 + t305 * t318;
t250 = t305 * t286 + t303 * t287;
t230 = pkin(4) * t250 + t305 * t263 - t303 * t318;
t247 = t305 * t284 - t303 * t285;
t343 = pkin(1) * t247;
t342 = pkin(2) * t298;
t341 = pkin(4) * t247;
t340 = qJ(1) * g(3);
t199 = pkin(3) * t202;
t331 = -pkin(2) * t191 - t199;
t330 = -pkin(2) * t286 - t244;
t216 = pkin(2) * t225;
t329 = pkin(1) * t206 + t216;
t325 = -t303 * t284 - t305 * t285;
t322 = pkin(1) * t186 - t331;
t289 = t305 * qJDD(2) - t303 * t306;
t321 = -pkin(4) * t289 - t303 * g(3);
t320 = -pkin(3) * t276 - t214;
t319 = pkin(2) * t287 - t243;
t317 = pkin(1) * t250 - t330;
t288 = t303 * qJDD(2) + t305 * t306;
t316 = pkin(1) * t288 - t285;
t315 = -pkin(2) * t241 + t320;
t198 = pkin(3) * t298 + pkin(5) * t328;
t181 = pkin(5) * t339 + qJ(3) * t356 + t301 * t198 + t342;
t183 = pkin(5) * t338 + qJ(3) * t191 - t300 * t198;
t314 = pkin(4) * t186 - t303 * t181 + t305 * t183;
t313 = -pkin(3) * t273 - t215;
t312 = -pkin(1) * t324 - t319;
t311 = -pkin(1) * t289 - t284;
t310 = -pkin(2) * t237 + t313;
t309 = pkin(1) * t221 - t315;
t213 = qJ(3) * t327 + t342;
t308 = pkin(4) * t206 + qJ(3) * t332 - t303 * t213;
t307 = pkin(1) * t217 - t310;
t294 = pkin(1) * t298;
t292 = qJ(1) * t298;
t270 = -pkin(4) * t288 + t305 * g(3);
t245 = pkin(1) * g(3) + pkin(4) * t325;
t193 = pkin(4) * t357 + qJ(3) * t336 + t305 * t213 + t294;
t180 = pkin(4) * t363 + t305 * t181 + t303 * t183 + t294;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, 0, -t299, -g(3), -t340, 0, 0, t289, 0, -t288, 0, t321, -t270, -t247, -t340 - t341, 0, 0, t324, 0, -t250, 0, t355, t230, t206, -t292 + t308, 0, 0, -t221, 0, -t217, 0, t366, t197, t186, -t292 + t314; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t288, 0, t289, 0, t270, t321, t325, t245, 0, 0, t250, 0, t324, 0, -t230, t355, t357, t193, 0, 0, t217, 0, -t221, 0, -t197, t366, t363, t180; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, -t299, 0, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, qJDD(2), qJ(1) * t288 - t311, qJ(1) * t289 - t316, 0, -qJ(1) * t325 + t343, 0, 0, 0, 0, 0, qJDD(2), qJ(1) * t250 - t312, qJ(1) * t324 - t317, 0, -qJ(1) * t357 - t329, 0, 0, 0, 0, 0, t296, qJ(1) * t217 - t309, -qJ(1) * t221 - t307, 0, -qJ(1) * t363 - t322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, -g(3), 0, 0, 0, t289, 0, -t288, 0, t321, -t270, -t247, -t341, 0, 0, t324, 0, -t250, 0, t355, t230, t206, t308, 0, 0, -t221, 0, -t217, 0, t366, t197, t186, t314; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, 0, -g(1), 0, 0, 0, 0, 0, 0, -qJDD(2), t311, t316, 0, -t343, 0, 0, 0, 0, 0, -qJDD(2), t312, t317, 0, t329, 0, 0, 0, 0, 0, -t296, t309, t307, 0, t322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t288, 0, t289, 0, t270, t321, t325, t245, 0, 0, t250, 0, t324, 0, -t230, t355, t357, t193, 0, 0, t217, 0, -t221, 0, -t197, t366, t363, t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t306, 0, 0, -g(3), -t284, 0, 0, 0, t287, 0, -t286, 0, t318, t263, t225, qJ(3) * t225, 0, 0, -t241, 0, -t237, 0, t360, t211, t191, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t306, 0, qJDD(2), 0, g(3), 0, -t285, 0, 0, 0, t286, 0, t287, 0, -t263, t318, t327, t213, 0, 0, t237, 0, -t241, 0, -t211, t360, t356, t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t284, t285, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t319, t330, 0, -t216, 0, 0, 0, 0, 0, t296, t315, t310, 0, t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t306, 0, 0, -t298, t243, 0, 0, 0, -t276, 0, -t273, 0, t354, t259, t202, pkin(5) * t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t306, 0, qJDD(2), 0, t298, 0, t244, 0, 0, 0, t273, 0, -t276, 0, -t259, t354, t328, t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t243, -t244, 0, 0, 0, 0, 0, 0, 0, t296, t320, t313, 0, -t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296, 0, -t295, 0, 0, -t298, t214, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, 0, t296, 0, t298, 0, t215, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296, -t214, -t215, 0, 0;];
m_new_reg  = t1;
