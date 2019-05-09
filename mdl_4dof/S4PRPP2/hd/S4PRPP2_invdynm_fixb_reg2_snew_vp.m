% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PRPP2
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
%   pkin=[a2,a3,a4,d2,theta3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PRPP2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:54:38
% EndTime: 2019-05-04 18:54:39
% DurationCPUTime: 1.16s
% Computational Cost: add. (2208->144), mult. (3184->120), div. (0->0), fcn. (1898->4), ass. (0->78)
t228 = sin(pkin(5));
t229 = cos(pkin(5));
t232 = qJD(2) ^ 2;
t210 = t228 * qJDD(2) + t229 * t232;
t211 = t229 * qJDD(2) - t228 * t232;
t230 = sin(qJ(2));
t231 = cos(qJ(2));
t250 = t231 * t210 + t230 * t211;
t226 = g(2) - qJDD(1);
t208 = t230 * g(1) - t231 * t226;
t199 = qJDD(2) * pkin(2) + t208;
t209 = t231 * g(1) + t230 * t226;
t200 = -t232 * pkin(2) - t209;
t172 = t228 * t199 + t229 * t200;
t273 = pkin(2) * t210;
t256 = -t273 - t172;
t240 = pkin(1) * t250 - t256;
t249 = -t230 * t210 + t231 * t211;
t287 = qJ(1) * t249 - t240;
t225 = g(3) - qJDD(3);
t187 = qJ(3) * t211 + t228 * t225;
t254 = -qJ(3) * t210 + t229 * t225;
t286 = -pkin(4) * t250 - t230 * t187 + t231 * t254;
t171 = -t229 * t199 + t228 * t200;
t252 = t228 * t171 + t229 * t172;
t157 = t229 * t171 - t228 * t172;
t265 = t231 * t157;
t143 = -t230 * t252 + t265;
t269 = t230 * t157;
t284 = t231 * t252 + t269;
t160 = pkin(4) * t249 + t231 * t187 + t230 * t254;
t260 = (qJD(4) * qJD(2));
t221 = 2 * t260;
t257 = t221 + t172;
t259 = qJDD(2) * qJ(4);
t167 = -t232 * pkin(3) + t257 + t259;
t227 = qJDD(2) * pkin(3);
t243 = qJDD(4) + t171;
t169 = -t232 * qJ(4) - t227 + t243;
t150 = t228 * t167 - t229 * t169;
t253 = t229 * t167 + t228 * t169;
t281 = -t230 * t150 + t231 * t253;
t138 = t231 * t150 + t230 * t253;
t176 = t231 * t208 - t230 * t209;
t274 = pkin(1) * t176;
t272 = pkin(4) * t176;
t271 = qJ(1) * g(3);
t261 = -pkin(3) * t169 + qJ(4) * t167;
t258 = pkin(2) * t150 + t261;
t154 = pkin(2) * t157;
t255 = pkin(1) * t143 + t154;
t251 = -t230 * t208 - t231 * t209;
t220 = 0.2e1 * t259;
t246 = t220 + t257;
t213 = t231 * qJDD(2) - t230 * t232;
t245 = -pkin(4) * t213 - t230 * g(3);
t206 = pkin(2) * t211;
t244 = t171 - t206;
t242 = -pkin(1) * t138 - t258;
t212 = t230 * qJDD(2) + t231 * t232;
t239 = pkin(1) * t212 - t209;
t238 = -pkin(1) * t249 + t244;
t145 = qJ(3) * t253 + (pkin(3) * t229 + qJ(4) * t228 + pkin(2)) * t225;
t147 = -qJ(3) * t150 + (-pkin(3) * t228 + qJ(4) * t229) * t225;
t237 = -pkin(4) * t138 - t230 * t145 + t231 * t147;
t223 = 0.2e1 * t227;
t236 = t223 - t243;
t235 = -pkin(1) * t213 - t208;
t234 = qJ(1) * t250 - t238;
t153 = pkin(2) * t225 + qJ(3) * t252;
t233 = pkin(4) * t143 + qJ(3) * t265 - t230 * t153;
t219 = pkin(1) * t225;
t217 = qJ(1) * t225;
t197 = -pkin(4) * t212 + t231 * g(3);
t173 = pkin(1) * g(3) + pkin(4) * t251;
t136 = pkin(4) * t284 + qJ(3) * t269 + t231 * t153 + t219;
t135 = pkin(4) * t281 + t231 * t145 + t230 * t147 + t219;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, 0, -t226, -g(3), -t271, 0, 0, t213, 0, -t212, 0, t245, -t197, -t176, -t271 - t272, 0, 0, t249, 0, -t250, 0, -t160, -t286, t143, -t217 + t233, 0, t249, 0, 0, t250, 0, -t160, -t138, t286, -t217 + t237; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t212, 0, t213, 0, t197, t245, t251, t173, 0, 0, t250, 0, t249, 0, t286, -t160, t284, t136, 0, t250, 0, 0, -t249, 0, t286, t281, t160, t135; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, -t226, 0, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, qJDD(2), qJ(1) * t212 - t235, qJ(1) * t213 - t239, 0, -qJ(1) * t251 + t274, 0, 0, 0, 0, 0, qJDD(2), t234, t287, 0, -qJ(1) * t284 - t255, 0, 0, 0, qJDD(2), 0, 0, -qJDD(4) + t223 + t234, 0, t220 + t221 - t287, -qJ(1) * t281 - t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t226, -g(3), 0, 0, 0, t213, 0, -t212, 0, t245, -t197, -t176, -t272, 0, 0, t249, 0, -t250, 0, -t160, -t286, t143, t233, 0, t249, 0, 0, t250, 0, -t160, -t138, t286, t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, 0, -g(1), 0, 0, 0, 0, 0, 0, -qJDD(2), t235, t239, 0, -t274, 0, 0, 0, 0, 0, -qJDD(2), t238, t240, 0, t255, 0, 0, 0, -qJDD(2), 0, 0, qJDD(4) - 0.2e1 * t227 + t238, 0, -t240 - 0.2e1 * t259 - (2 * t260), t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t212, 0, t213, 0, t197, t245, t251, t173, 0, 0, t250, 0, t249, 0, t286, -t160, t284, t136, 0, t250, 0, 0, -t249, 0, t286, t281, t160, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t232, 0, 0, -g(3), -t208, 0, 0, 0, t211, 0, -t210, 0, -t187, -t254, t157, qJ(3) * t157, 0, t211, 0, 0, t210, 0, -t187, -t150, t254, t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, 0, qJDD(2), 0, g(3), 0, -t209, 0, 0, 0, t210, 0, t211, 0, t254, -t187, t252, t153, 0, t210, 0, 0, -t211, 0, t254, t253, t187, t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t208, t209, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t244, t256, 0, -t154, 0, 0, 0, qJDD(2), 0, 0, t206 + t236, 0, t273 + t246, t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t232, 0, 0, -t225, t171, 0, 0, qJDD(2), 0, 0, t232, 0, 0, t169, t225, qJ(4) * t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, 0, qJDD(2), 0, t225, 0, t172, 0, 0, t232, 0, 0, -qJDD(2), 0, t225, t167, 0, pkin(3) * t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t171, -t172, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t236, 0, t246, t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t232, 0, 0, t169, t225, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, -t169, 0, t167, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t232, 0, 0, qJDD(2), 0, -t225, -t167, 0, 0;];
m_new_reg  = t1;
