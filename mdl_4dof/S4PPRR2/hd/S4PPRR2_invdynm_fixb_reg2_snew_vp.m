% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PPRR2
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
%   pkin=[a2,a3,a4,d3,d4,theta2]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PPRR2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:50:58
% EndTime: 2019-05-04 18:50:59
% DurationCPUTime: 1.50s
% Computational Cost: add. (4682->150), mult. (6145->141), div. (0->0), fcn. (4978->6), ass. (0->88)
t277 = qJD(3) + qJD(4);
t275 = t277 ^ 2;
t276 = qJDD(3) + qJDD(4);
t281 = sin(qJ(4));
t283 = cos(qJ(4));
t256 = t283 * t275 + t281 * t276;
t259 = t281 * t275 - t283 * t276;
t282 = sin(qJ(3));
t284 = cos(qJ(3));
t220 = t284 * t256 - t282 * t259;
t278 = g(3) - qJDD(2);
t244 = pkin(5) * t256 - t283 * t278;
t333 = pkin(5) * t259 - t281 * t278;
t196 = pkin(4) * t220 + t284 * t244 - t282 * t333;
t224 = t282 * t256 + t284 * t259;
t279 = sin(pkin(6));
t280 = cos(pkin(6));
t197 = t280 * t220 - t279 * t224;
t339 = pkin(4) * t224 + t282 * t244 + t284 * t333;
t182 = qJ(2) * t197 + t280 * t196 - t279 * t339;
t201 = t279 * t220 + t280 * t224;
t345 = qJ(2) * t201 + t279 * t196 + t280 * t339;
t310 = g(2) - qJDD(1);
t262 = t279 * g(1) - t280 * t310;
t263 = t280 * g(1) + t279 * t310;
t301 = t284 * t262 + t282 * t263;
t227 = qJDD(3) * pkin(3) + t301;
t285 = qJD(3) ^ 2;
t308 = -t282 * t262 + t284 * t263;
t228 = -t285 * pkin(3) - t308;
t209 = -t283 * t227 + t281 * t228;
t210 = t281 * t227 + t283 * t228;
t305 = t281 * t209 + t283 * t210;
t187 = t283 * t209 - t281 * t210;
t311 = t284 * t187;
t176 = -t282 * t305 + t311;
t313 = t282 * t187;
t336 = t284 * t305 + t313;
t342 = t279 * t176 + t280 * t336;
t171 = t280 * t176 - t279 * t336;
t304 = -t282 * t301 - t284 * t308;
t214 = t282 * t308 - t284 * t301;
t316 = t280 * t214;
t191 = -t279 * t304 + t316;
t319 = t279 * t214;
t335 = t280 * t304 + t319;
t266 = t282 * qJDD(3) + t284 * t285;
t248 = pkin(4) * t266 - t284 * t278;
t267 = t284 * qJDD(3) - t282 * t285;
t295 = -pkin(4) * t267 - t282 * t278;
t300 = -t279 * t266 + t280 * t267;
t334 = -qJ(2) * t300 + t279 * t248 + t280 * t295;
t235 = t280 * t266 + t279 * t267;
t208 = qJ(2) * t235 + t280 * t248 - t279 * t295;
t230 = t280 * t262 - t279 * t263;
t322 = pkin(1) * t230;
t321 = pkin(2) * t278;
t320 = qJ(2) * t230;
t268 = t279 * t278;
t315 = t280 * t278;
t184 = pkin(3) * t187;
t309 = -pkin(2) * t176 - t184;
t307 = -pkin(2) * t266 + t308;
t211 = pkin(2) * t214;
t306 = pkin(1) * t191 + t211;
t302 = -t279 * t262 - t280 * t263;
t298 = pkin(1) * t171 - t309;
t297 = -pkin(3) * t259 - t209;
t296 = pkin(2) * t267 + t301;
t294 = pkin(1) * t235 - t307;
t293 = -pkin(2) * t224 + t297;
t183 = pkin(3) * t278 + pkin(5) * t305;
t166 = pkin(4) * t336 + pkin(5) * t313 + t284 * t183 + t321;
t168 = pkin(4) * t176 + pkin(5) * t311 - t282 * t183;
t292 = qJ(2) * t171 - t279 * t166 + t280 * t168;
t291 = -pkin(3) * t256 - t210;
t290 = -pkin(1) * t300 - t296;
t289 = -pkin(2) * t220 + t291;
t288 = pkin(1) * t201 - t293;
t204 = pkin(4) * t304 + t321;
t287 = pkin(4) * t316 + qJ(2) * t191 - t279 * t204;
t286 = pkin(1) * t197 - t289;
t274 = pkin(1) * t278;
t272 = qJ(1) * t278;
t216 = qJ(2) * t302 + t274;
t178 = pkin(4) * t319 + qJ(2) * t335 + t280 * t204 + t274;
t165 = qJ(2) * t342 + t280 * t166 + t279 * t168 + t274;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, 0, -t310, -g(3), -qJ(1) * g(3), 0, 0, 0, 0, 0, 0, -t268, -t315, -t230, -t272 - t320, 0, 0, t300, 0, -t235, 0, t334, t208, t191, -t272 + t287, 0, 0, -t201, 0, -t197, 0, t345, t182, t171, -t272 + t292; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t315, -t268, t302, t216, 0, 0, t235, 0, t300, 0, -t208, t334, t335, t178, 0, 0, t197, 0, -t201, 0, -t182, t345, t342, t165; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, -t310, 0, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, 0, t262, t263, 0, -qJ(1) * t302 + t322, 0, 0, 0, 0, 0, qJDD(3), qJ(1) * t235 - t290, qJ(1) * t300 - t294, 0, -qJ(1) * t335 - t306, 0, 0, 0, 0, 0, t276, qJ(1) * t197 - t288, -qJ(1) * t201 - t286, 0, -qJ(1) * t342 - t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t310, -g(3), 0, 0, 0, 0, 0, 0, 0, -t268, -t315, -t230, -t320, 0, 0, t300, 0, -t235, 0, t334, t208, t191, t287, 0, 0, -t201, 0, -t197, 0, t345, t182, t171, t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, -t262, -t263, 0, -t322, 0, 0, 0, 0, 0, -qJDD(3), t290, t294, 0, t306, 0, 0, 0, 0, 0, -t276, t288, t286, 0, t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t315, -t268, t302, t216, 0, 0, t235, 0, t300, 0, -t208, t334, t335, t178, 0, 0, t197, 0, -t201, 0, -t182, t345, t342, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, -t262, 0, 0, 0, t267, 0, -t266, 0, t295, t248, t214, pkin(4) * t214, 0, 0, -t224, 0, -t220, 0, t339, t196, t176, t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, 0, -t263, 0, 0, 0, t266, 0, t267, 0, -t248, t295, t304, t204, 0, 0, t220, 0, -t224, 0, -t196, t339, t336, t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262, t263, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t296, t307, 0, -t211, 0, 0, 0, 0, 0, t276, t293, t289, 0, t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t285, 0, 0, -t278, -t301, 0, 0, 0, -t259, 0, -t256, 0, t333, t244, t187, pkin(5) * t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, 0, qJDD(3), 0, t278, 0, -t308, 0, 0, 0, t256, 0, -t259, 0, -t244, t333, t305, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t301, t308, 0, 0, 0, 0, 0, 0, 0, t276, t297, t291, 0, -t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276, 0, -t275, 0, 0, -t278, t209, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, 0, t276, 0, t278, 0, t210, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276, -t209, -t210, 0, 0;];
m_new_reg  = t1;
