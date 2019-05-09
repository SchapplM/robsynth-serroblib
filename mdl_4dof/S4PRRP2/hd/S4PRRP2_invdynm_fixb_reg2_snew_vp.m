% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PRRP2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PRRP2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:03:31
% EndTime: 2019-05-04 19:03:32
% DurationCPUTime: 1.22s
% Computational Cost: add. (2642->151), mult. (3228->123), div. (0->0), fcn. (2026->4), ass. (0->87)
t242 = qJD(2) + qJD(3);
t240 = t242 ^ 2;
t243 = g(3) - qJDD(4);
t234 = -qJ(4) * t240 + t243;
t247 = cos(qJ(3));
t241 = qJDD(2) + qJDD(3);
t245 = sin(qJ(3));
t271 = t241 * t245;
t224 = t240 * t247 + t271;
t292 = pkin(5) * t224;
t194 = -qJ(4) * t271 + t234 * t247 - t292;
t246 = sin(qJ(2));
t248 = cos(qJ(2));
t270 = t241 * t247;
t227 = t240 * t245 - t270;
t291 = pkin(5) * t227;
t294 = -qJ(4) * t270 - t234 * t245 + t291;
t202 = t224 * t246 + t227 * t248;
t300 = pkin(4) * t202;
t303 = -t194 * t246 + t248 * t294 + t300;
t198 = t224 * t248 - t227 * t246;
t301 = pkin(4) * t198;
t164 = t194 * t248 + t246 * t294 - t301;
t213 = g(3) * t247 - t292;
t293 = -g(3) * t245 + t291;
t302 = -t213 * t246 + t248 * t293 + t300;
t174 = t213 * t248 + t246 * t293 - t301;
t244 = g(2) - qJDD(1);
t231 = g(1) * t248 + t244 * t246;
t250 = qJD(2) ^ 2;
t229 = -pkin(2) * t250 - t231;
t230 = g(1) * t246 - t244 * t248;
t253 = qJDD(2) * pkin(2) + t230;
t204 = t229 * t245 - t247 * t253;
t205 = t229 * t247 + t245 * t253;
t266 = t204 * t245 + t205 * t247;
t182 = t204 * t247 - t205 * t245;
t274 = t182 * t248;
t161 = -t246 * t266 + t274;
t275 = t182 * t246;
t295 = t248 * t266 + t275;
t189 = -pkin(3) * t240 + t205;
t237 = t241 * pkin(3);
t188 = t204 - t237;
t272 = t188 * t247;
t172 = -t189 * t245 + t272;
t273 = t188 * t245;
t267 = t189 * t247 + t273;
t290 = t172 * t246 + t248 * t267;
t157 = t172 * t248 - t246 * t267;
t208 = t230 * t248 - t231 * t246;
t285 = pkin(1) * t208;
t282 = pkin(4) * t208;
t277 = qJ(1) * g(3);
t276 = qJ(4) * t241;
t187 = pkin(3) * t188;
t269 = -pkin(2) * t172 - t187;
t179 = pkin(2) * t182;
t268 = pkin(1) * t161 + t179;
t264 = -t230 * t246 - t231 * t248;
t263 = pkin(1) * t157 - t269;
t233 = qJDD(2) * t248 - t246 * t250;
t262 = -pkin(4) * t233 - g(3) * t246;
t223 = pkin(2) * t227;
t261 = t204 + t223;
t236 = 0.2e1 * t237;
t260 = -t204 + t236;
t232 = qJDD(2) * t246 + t248 * t250;
t259 = pkin(1) * t232 - t231;
t258 = pkin(1) * t202 + t261;
t185 = pkin(3) * t243 + qJ(4) * t189;
t152 = pkin(2) * t243 + pkin(5) * t267 + qJ(4) * t273 + t185 * t247;
t154 = pkin(5) * t172 + qJ(4) * t272 - t185 * t245;
t257 = pkin(4) * t157 - t152 * t246 + t154 * t248;
t254 = -pkin(1) * t233 - t230;
t178 = pkin(2) * g(3) + pkin(5) * t266;
t252 = pkin(4) * t161 + pkin(5) * t274 - t178 * t246;
t251 = qJ(1) * t198 - t258;
t184 = -pkin(2) * t224 - t205;
t169 = pkin(1) * t198 - t184;
t249 = pkin(1) * g(3);
t219 = -pkin(4) * t232 + g(3) * t248;
t206 = pkin(4) * t264 + t249;
t163 = -qJ(1) * t202 - t169;
t151 = pkin(4) * t295 + pkin(5) * t275 + t178 * t248 + t249;
t150 = pkin(1) * t243 + pkin(4) * t290 + t152 * t248 + t154 * t246;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, 0, -t244, -g(3), -t277, 0, 0, t233, 0, -t232, 0, t262, -t219, -t208, -t277 - t282, 0, 0, -t202, 0, -t198, 0, t302, -t174, t161, t252 - t277, 0, 0, -t202, 0, -t198, 0, t303, -t164, t157, -qJ(1) * t243 + t257; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t232, 0, t233, 0, t219, t262, t264, t206, 0, 0, t198, 0, -t202, 0, t174, t302, t295, t151, 0, 0, t198, 0, -t202, 0, t164, t303, t290, t150; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, -t244, 0, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, qJDD(2), qJ(1) * t232 - t254, qJ(1) * t233 - t259, 0, -qJ(1) * t264 + t285, 0, 0, 0, 0, 0, t241, t251, t163, 0, -qJ(1) * t295 - t268, 0, 0, 0, 0, 0, t241, t236 + t251, t163, 0, -qJ(1) * t290 - t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, -g(3), 0, 0, 0, t233, 0, -t232, 0, t262, -t219, -t208, -t282, 0, 0, -t202, 0, -t198, 0, t302, -t174, t161, t252, 0, 0, -t202, 0, -t198, 0, t303, -t164, t157, t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, 0, -g(1), 0, 0, 0, 0, 0, 0, -qJDD(2), t254, t259, 0, -t285, 0, 0, 0, 0, 0, -t241, t258, t169, 0, t268, 0, 0, 0, 0, 0, -t241, -0.2e1 * t237 + t258, t169, 0, t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t232, 0, t233, 0, t219, t262, t264, t206, 0, 0, t198, 0, -t202, 0, t174, t302, t295, t151, 0, 0, t198, 0, -t202, 0, t164, t303, t290, t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t250, 0, 0, -g(3), -t230, 0, 0, 0, -t227, 0, -t224, 0, t293, -t213, t182, pkin(5) * t182, 0, 0, -t227, 0, -t224, 0, t294, -t194, t172, t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250, 0, qJDD(2), 0, g(3), 0, -t231, 0, 0, 0, t224, 0, -t227, 0, t213, t293, t266, t178, 0, 0, t224, 0, -t227, 0, t194, t294, t267, t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t230, t231, 0, 0, 0, 0, 0, 0, 0, t241, -t261, t184, 0, -t179, 0, 0, 0, 0, 0, t241, -t223 + t260, t184, 0, t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, 0, -t240, 0, 0, -g(3), t204, 0, 0, 0, t241, 0, -t240, 0, -t276, -t234, t188, qJ(4) * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, 0, t241, 0, g(3), 0, t205, 0, 0, 0, t240, 0, t241, 0, t234, -t276, t189, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, -t204, -t205, 0, 0, 0, 0, 0, 0, 0, t241, t260, -t205, 0, -t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, 0, -t240, 0, 0, -t243, t188, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, 0, t241, 0, t243, 0, t189, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, -t188, -t189, 0, 0;];
m_new_reg  = t1;
