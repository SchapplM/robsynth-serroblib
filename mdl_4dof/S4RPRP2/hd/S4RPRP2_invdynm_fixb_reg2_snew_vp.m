% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RPRP2
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
%   pkin=[a2,a3,a4,d1,d3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RPRP2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:15:24
% EndTime: 2019-05-04 19:15:26
% DurationCPUTime: 1.21s
% Computational Cost: add. (2691->163), mult. (3788->139), div. (0->0), fcn. (1516->4), ass. (0->90)
t237 = sin(qJ(1));
t239 = cos(qJ(1));
t232 = -qJD(1) + qJD(3);
t230 = t232 ^ 2;
t234 = g(3) + qJDD(4);
t223 = qJ(4) * t230 + t234;
t236 = sin(qJ(3));
t231 = qJDD(1) - qJDD(3);
t238 = cos(qJ(3));
t261 = t238 * t231;
t250 = t236 * t230 + t261;
t283 = pkin(5) * t250;
t290 = qJ(4) * t261 + t236 * t223 + t283;
t263 = t236 * t231;
t214 = t238 * t230 - t263;
t291 = pkin(5) * t214;
t293 = -qJ(4) * t263 + t238 * t223 + t291;
t284 = t239 * t214 + t237 * t250;
t302 = pkin(4) * t284;
t304 = t237 * t290 + t239 * t293 - t302;
t289 = t236 * g(3) + t283;
t292 = t238 * g(3) + t291;
t303 = t237 * t289 + t239 * t292 - t302;
t195 = t237 * t214 - t239 * t250;
t273 = pkin(4) * t195;
t299 = t237 * t292 - t239 * t289 - t273;
t298 = t237 * t293 - t239 * t290 - t273;
t260 = (qJD(2) * qJD(1));
t228 = 2 * t260;
t241 = qJD(1) ^ 2;
t222 = t239 * g(1) + t237 * g(2);
t233 = qJDD(1) * qJ(2);
t249 = -t222 + t233;
t276 = pkin(1) + pkin(2);
t199 = -t276 * t241 + t228 + t249;
t221 = t237 * g(1) - t239 * g(2);
t255 = -qJDD(2) + t221;
t245 = -t241 * qJ(2) - t255;
t242 = -t276 * qJDD(1) + t245;
t184 = t238 * t199 + t236 * t242;
t182 = -t230 * pkin(3) + t184;
t183 = t236 * t199 - t238 * t242;
t226 = t231 * pkin(3);
t181 = t183 + t226;
t262 = t238 * t181;
t168 = -t236 * t182 + t262;
t264 = t236 * t181;
t254 = t238 * t182 + t264;
t288 = t237 * t168 - t239 * t254;
t287 = t239 * t168 + t237 * t254;
t174 = t238 * t183 - t236 * t184;
t253 = t236 * t183 + t238 * t184;
t286 = t237 * t174 - t239 * t253;
t285 = t239 * t174 + t237 * t253;
t275 = pkin(3) * t181;
t272 = pkin(5) * t174;
t266 = qJ(4) * t231;
t265 = qJDD(1) * pkin(1);
t204 = t241 * pkin(1) - t249 - (2 * t260);
t206 = -t245 + t265;
t259 = -t239 * t204 - t237 * t206;
t257 = -t237 * t221 - t239 * t222;
t256 = -pkin(2) * g(3) + pkin(5) * t253;
t219 = t237 * qJDD(1) + t239 * t241;
t208 = -pkin(4) * t219 + t239 * g(3);
t220 = t239 * qJDD(1) - t237 * t241;
t207 = pkin(4) * t220 + t237 * g(3);
t252 = t237 * t204 - t239 * t206;
t251 = t239 * t221 - t237 * t222;
t248 = pkin(2) * t250 + t183;
t178 = pkin(1) * t250 - qJ(2) * t214 + t248;
t180 = -pkin(3) * t234 + qJ(4) * t182;
t244 = pkin(5) * t168 + qJ(4) * t262 - t236 * t180;
t243 = -pkin(5) * t254 - qJ(4) * t264 - t238 * t180;
t240 = pkin(1) * g(3);
t235 = qJ(2) * g(3);
t225 = 0.2e1 * t226;
t217 = t255 + 0.2e1 * t265;
t212 = -t222 + t228 + 0.2e1 * t233;
t185 = pkin(1) * t206 - qJ(2) * t204;
t179 = pkin(2) * t214 + t184;
t177 = qJ(2) * t250 + t214 * t276 + t184;
t176 = t178 + t225;
t171 = t235 + t272;
t170 = t240 - t256;
t165 = qJ(2) * t234 + t244;
t164 = t276 * t234 + t243;
t163 = qJ(2) * t253 + t174 * t276;
t162 = qJ(2) * t254 + t168 * t276 + t275;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t220, 0, -t219, 0, -t207, -t208, -t251, -pkin(4) * t251, 0, t220, 0, 0, t219, 0, -t207, t252, t208, pkin(4) * t252 + (-t237 * pkin(1) + t239 * qJ(2)) * g(3), 0, 0, t195, 0, -t284, 0, -t299, t303, t285, -pkin(4) * t285 - t237 * t170 + t239 * t171, 0, 0, t195, 0, -t284, 0, -t298, t304, t287, -pkin(4) * t287 - t237 * t164 + t239 * t165; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t219, 0, t220, 0, t208, -t207, t257, pkin(4) * t257, 0, t219, 0, 0, -t220, 0, t208, t259, t207, pkin(4) * t259 + (t239 * pkin(1) + t237 * qJ(2)) * g(3), 0, 0, -t284, 0, -t195, 0, t303, t299, t286, -pkin(4) * t286 + t239 * t170 + t237 * t171, 0, 0, -t284, 0, -t195, 0, t304, t298, t288, -pkin(4) * t288 + t239 * t164 + t237 * t165; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t221, t222, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t217, 0, t212, t185, 0, 0, 0, 0, 0, t231, t178, t177, 0, t163, 0, 0, 0, 0, 0, t231, t176, t177, 0, t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t241, 0, 0, -g(3), -t221, 0, 0, qJDD(1), 0, 0, t241, 0, 0, -t206, g(3), t235, 0, 0, -t250, 0, -t214, 0, t289, t292, t174, t171, 0, 0, -t250, 0, -t214, 0, t290, t293, t168, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, 0, qJDD(1), 0, g(3), 0, -t222, 0, 0, t241, 0, 0, -qJDD(1), 0, g(3), -t204, 0, t240, 0, 0, -t214, 0, t250, 0, t292, -t289, -t253, t170, 0, 0, -t214, 0, t250, 0, t293, -t290, -t254, t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t221, t222, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t217, 0, t212, t185, 0, 0, 0, 0, 0, t231, t178, t177, 0, t163, 0, 0, 0, 0, 0, t231, t176, t177, 0, t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t241, 0, 0, -t206, g(3), 0, 0, 0, -t250, 0, -t214, 0, t289, t292, t174, t272, 0, 0, -t250, 0, -t214, 0, t290, t293, t168, t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t206, 0, -t204, 0, 0, 0, 0, 0, 0, t231, t248, t179, 0, pkin(2) * t174, 0, 0, 0, 0, 0, t231, t225 + t248, t179, 0, pkin(2) * t168 + t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, 0, 0, qJDD(1), 0, -g(3), t204, 0, 0, 0, 0, t214, 0, -t250, 0, -t292, t289, t253, t256, 0, 0, t214, 0, -t250, 0, -t293, t290, t254, -pkin(2) * t234 - t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, 0, -t230, 0, 0, g(3), t183, 0, 0, 0, -t231, 0, -t230, 0, t266, t223, t181, qJ(4) * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, 0, -t231, 0, -g(3), 0, t184, 0, 0, 0, t230, 0, -t231, 0, -t223, t266, t182, t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, -t183, -t184, 0, 0, 0, 0, 0, 0, 0, -t231, -t183 - 0.2e1 * t226, -t184, 0, -t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, 0, -t230, 0, 0, t234, t181, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, 0, -t231, 0, -t234, 0, t182, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, -t181, -t182, 0, 0;];
m_new_reg  = t1;
