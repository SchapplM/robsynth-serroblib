% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PPRR1
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
% Datum: 2019-05-04 18:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PPRR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:47:49
% EndTime: 2019-05-04 18:47:51
% DurationCPUTime: 1.32s
% Computational Cost: add. (3227->160), mult. (4729->149), div. (0->0), fcn. (3800->6), ass. (0->84)
t229 = qJD(3) + qJD(4);
t227 = t229 ^ 2;
t228 = qJDD(3) + qJDD(4);
t233 = sin(qJ(4));
t235 = cos(qJ(4));
t207 = t227 * t235 + t228 * t233;
t210 = t227 * t233 - t228 * t235;
t234 = sin(qJ(3));
t236 = cos(qJ(3));
t183 = t207 * t236 - t210 * t234;
t230 = g(3) - qJDD(1);
t201 = pkin(5) * t207 + t230 * t235;
t278 = pkin(5) * t210 + t230 * t233;
t169 = pkin(4) * t183 + t201 * t236 - t234 * t278;
t231 = sin(pkin(6));
t232 = cos(pkin(6));
t187 = t207 * t234 + t210 * t236;
t275 = t183 * t231 - t187 * t232;
t283 = pkin(4) * t187 + t201 * t234 + t236 * t278;
t293 = -qJ(1) * t275 + t231 * t169 - t232 * t283;
t284 = t183 * t232 + t187 * t231;
t292 = -qJ(1) * t284 + t232 * t169 + t231 * t283;
t221 = g(1) * t231 - t232 * g(2);
t218 = -qJDD(2) + t221;
t222 = g(1) * t232 + g(2) * t231;
t196 = t236 * t218 - t222 * t234;
t192 = qJDD(3) * pkin(3) - t196;
t197 = -t234 * t218 - t236 * t222;
t237 = qJD(3) ^ 2;
t193 = -pkin(3) * t237 + t197;
t175 = -t235 * t192 + t193 * t233;
t176 = t233 * t192 + t235 * t193;
t249 = t175 * t233 + t235 * t176;
t164 = t175 * t235 - t176 * t233;
t256 = t164 * t236;
t156 = -t234 * t249 + t256;
t257 = t164 * t234;
t280 = t236 * t249 + t257;
t290 = t156 * t231 - t232 * t280;
t289 = t156 * t232 + t231 * t280;
t219 = qJDD(3) * t234 + t236 * t237;
t205 = pkin(4) * t219 + t230 * t236;
t220 = -t236 * qJDD(3) + t234 * t237;
t243 = pkin(4) * t220 + t230 * t234;
t265 = t232 * t219 + t220 * t231;
t279 = -qJ(1) * t265 + t232 * t205 + t231 * t243;
t181 = t196 * t236 - t197 * t234;
t242 = t196 * t234 + t197 * t236;
t277 = t181 * t231 - t232 * t242;
t276 = t181 * t232 + t231 * t242;
t194 = -t219 * t231 + t232 * t220;
t274 = qJ(1) * t194 + t231 * t205 - t232 * t243;
t263 = pkin(1) + pkin(2);
t262 = pkin(1) * t230;
t261 = pkin(3) * t164;
t260 = pkin(4) * t181;
t258 = qJ(2) * t230;
t251 = t231 * t230;
t223 = t232 * t230;
t225 = pkin(2) * t230;
t250 = pkin(4) * t242 - t225;
t215 = t231 * t222;
t248 = t218 * t232 - t215;
t247 = t221 * t232 - t215;
t216 = t232 * t222;
t246 = -t218 * t231 - t216;
t245 = -t221 * t231 - t216;
t244 = -pkin(3) * t207 - t176;
t240 = -pkin(3) * t210 - t175;
t161 = -pkin(3) * t230 + pkin(5) * t249;
t239 = pkin(4) * t156 + pkin(5) * t256 - t161 * t234;
t238 = pkin(4) * t280 + pkin(5) * t257 + t161 * t236 - t225;
t198 = pkin(1) * t218 - qJ(2) * t222;
t178 = -qJ(2) * t219 + t263 * t220 + t196;
t177 = qJ(2) * t220 + t263 * t219 + t197;
t174 = t258 + t260;
t172 = -t250 + t262;
t160 = qJ(2) * t242 + t181 * t263;
t159 = -qJ(2) * t183 + t187 * t263 - t240;
t158 = qJ(2) * t187 + t183 * t263 - t244;
t153 = t239 + t258;
t152 = -t238 + t262;
t151 = qJ(2) * t280 + t156 * t263 + t261;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t251, -t223, -t247, -qJ(1) * t247, 0, 0, 0, 0, 0, 0, -t251, -t248, t223, -qJ(1) * t248 + (-pkin(1) * t231 + qJ(2) * t232) * t230, 0, 0, -t194, 0, -t265, 0, -t274, t279, t276, -qJ(1) * t276 - t231 * t172 + t232 * t174, 0, 0, t275, 0, -t284, 0, -t293, t292, t289, -qJ(1) * t289 - t231 * t152 + t232 * t153; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t223, -t251, t245, qJ(1) * t245, 0, 0, 0, 0, 0, 0, t223, t246, t251, qJ(1) * t246 + (pkin(1) * t232 + qJ(2) * t231) * t230, 0, 0, -t265, 0, t194, 0, t279, t274, t277, -qJ(1) * t277 + t232 * t172 + t231 * t174, 0, 0, -t284, 0, -t275, 0, t292, t293, t290, -qJ(1) * t290 + t232 * t152 + t231 * t153; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t221, t222, 0, 0, 0, 0, 0, 0, 0, 0, t218, 0, -t222, t198, 0, 0, 0, 0, 0, -qJDD(3), t178, t177, 0, t160, 0, 0, 0, 0, 0, -t228, t159, t158, 0, t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, -t221, 0, 0, 0, 0, 0, 0, 0, 0, -t218, t230, t258, 0, 0, -t220, 0, -t219, 0, t243, t205, t181, t174, 0, 0, -t187, 0, -t183, 0, t283, t169, t156, t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, 0, -t222, 0, 0, 0, 0, 0, 0, 0, t230, -t222, 0, t262, 0, 0, -t219, 0, t220, 0, t205, -t243, -t242, t172, 0, 0, -t183, 0, t187, 0, t169, -t283, -t280, t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t222, 0, 0, 0, 0, 0, 0, 0, 0, t218, 0, -t222, t198, 0, 0, 0, 0, 0, -qJDD(3), t178, t177, 0, t160, 0, 0, 0, 0, 0, -t228, t159, t158, 0, t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t218, t230, 0, 0, 0, -t220, 0, -t219, 0, t243, t205, t181, t260, 0, 0, -t187, 0, -t183, 0, t283, t169, t156, t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, 0, -t222, 0, 0, 0, 0, 0, 0, -qJDD(3), pkin(2) * t220 + t196, pkin(2) * t219 + t197, 0, pkin(2) * t181, 0, 0, 0, 0, 0, -t228, pkin(2) * t187 - t240, pkin(2) * t183 - t244, 0, pkin(2) * t156 + t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, t222, 0, 0, 0, 0, t219, 0, -t220, 0, -t205, t243, t242, t250, 0, 0, t183, 0, -t187, 0, -t169, t283, t280, t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t237, 0, 0, t230, t196, 0, 0, 0, -t210, 0, -t207, 0, t278, t201, t164, pkin(5) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, 0, qJDD(3), 0, -t230, 0, t197, 0, 0, 0, t207, 0, -t210, 0, -t201, t278, t249, t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t196, -t197, 0, 0, 0, 0, 0, 0, 0, t228, t240, t244, 0, -t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, 0, -t227, 0, 0, t230, t175, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, 0, t228, 0, -t230, 0, t176, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, -t175, -t176, 0, 0;];
m_new_reg  = t1;
