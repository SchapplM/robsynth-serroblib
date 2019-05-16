% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PPRP2
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
%   pkin=[a2,a3,a4,d3,theta2]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PPRP2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:44:35
% EndTime: 2019-05-04 18:44:37
% DurationCPUTime: 1.06s
% Computational Cost: add. (1824->136), mult. (2562->107), div. (0->0), fcn. (1854->4), ass. (0->71)
t212 = sin(qJ(3));
t213 = cos(qJ(3));
t214 = qJD(3) ^ 2;
t193 = qJDD(3) * t212 + t213 * t214;
t194 = qJDD(3) * t213 - t212 * t214;
t210 = sin(pkin(5));
t211 = cos(pkin(5));
t229 = t211 * t193 + t194 * t210;
t243 = g(2) - qJDD(1);
t188 = t210 * g(1) - t211 * t243;
t189 = g(1) * t211 + t210 * t243;
t241 = -t212 * t188 + t213 * t189;
t255 = pkin(2) * t193;
t236 = -t255 + t241;
t220 = pkin(1) * t229 - t236;
t228 = -t193 * t210 + t211 * t194;
t269 = qJ(1) * t228 - t220;
t208 = g(3) - qJDD(2);
t175 = pkin(4) * t194 + t208 * t212;
t234 = -pkin(4) * t193 + t213 * t208;
t268 = -qJ(2) * t229 - t175 * t210 + t211 * t234;
t230 = t213 * t188 + t189 * t212;
t232 = -t212 * t230 - t213 * t241;
t151 = t212 * t241 - t213 * t230;
t251 = t151 * t211;
t131 = -t210 * t232 + t251;
t252 = t151 * t210;
t266 = t211 * t232 + t252;
t143 = qJ(2) * t228 + t175 * t211 + t210 * t234;
t240 = (qJD(4) * qJD(3));
t204 = 2 * t240;
t237 = t204 - t241;
t239 = qJDD(3) * qJ(4);
t155 = -pkin(3) * t214 + t237 + t239;
t209 = qJDD(3) * pkin(3);
t223 = -qJDD(4) + t230;
t157 = -t214 * qJ(4) - t209 - t223;
t138 = t155 * t212 - t157 * t213;
t233 = t213 * t155 + t157 * t212;
t263 = -t138 * t210 + t211 * t233;
t126 = t138 * t211 + t210 * t233;
t161 = t188 * t211 - t189 * t210;
t256 = pkin(1) * t161;
t253 = qJ(2) * t161;
t195 = t210 * t208;
t244 = t211 * t208;
t242 = -pkin(3) * t157 + qJ(4) * t155;
t238 = pkin(2) * t138 + t242;
t148 = pkin(2) * t151;
t235 = pkin(1) * t131 + t148;
t231 = -t188 * t210 - t211 * t189;
t203 = 0.2e1 * t239;
t225 = t203 + t237;
t191 = pkin(2) * t194;
t224 = t191 + t230;
t222 = -pkin(1) * t126 - t238;
t219 = pkin(1) * t228 + t224;
t206 = 0.2e1 * t209;
t218 = t206 + t223;
t133 = pkin(4) * t233 + (pkin(3) * t213 + qJ(4) * t212 + pkin(2)) * t208;
t135 = -pkin(4) * t138 + (-pkin(3) * t212 + qJ(4) * t213) * t208;
t217 = -qJ(2) * t126 - t133 * t210 + t211 * t135;
t216 = qJ(1) * t229 + t219;
t141 = pkin(2) * t208 + pkin(4) * t232;
t215 = pkin(4) * t251 + qJ(2) * t131 - t141 * t210;
t202 = pkin(1) * t208;
t200 = qJ(1) * t208;
t158 = qJ(2) * t231 + t202;
t124 = pkin(4) * t252 + qJ(2) * t266 + t141 * t211 + t202;
t123 = qJ(2) * t263 + t133 * t211 + t135 * t210 + t202;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, 0, -t243, -g(3), -qJ(1) * g(3), 0, 0, 0, 0, 0, 0, -t195, -t244, -t161, -t200 - t253, 0, 0, t228, 0, -t229, 0, -t143, -t268, t131, -t200 + t215, 0, t228, 0, 0, t229, 0, -t143, -t126, t268, -t200 + t217; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t244, -t195, t231, t158, 0, 0, t229, 0, t228, 0, t268, -t143, t266, t124, 0, t229, 0, 0, -t228, 0, t268, t263, t143, t123; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, -t243, 0, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, 0, t188, t189, 0, -qJ(1) * t231 + t256, 0, 0, 0, 0, 0, qJDD(3), t216, t269, 0, -qJ(1) * t266 - t235, 0, 0, 0, qJDD(3), 0, 0, -qJDD(4) + t206 + t216, 0, t203 + t204 - t269, -qJ(1) * t263 - t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, -g(3), 0, 0, 0, 0, 0, 0, 0, -t195, -t244, -t161, -t253, 0, 0, t228, 0, -t229, 0, -t143, -t268, t131, t215, 0, t228, 0, 0, t229, 0, -t143, -t126, t268, t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, -t188, -t189, 0, -t256, 0, 0, 0, 0, 0, -qJDD(3), -t219, t220, 0, t235, 0, 0, 0, -qJDD(3), 0, 0, qJDD(4) - 0.2e1 * t209 - t219, 0, -t220 - 0.2e1 * t239 - (2 * t240), t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t244, -t195, t231, t158, 0, 0, t229, 0, t228, 0, t268, -t143, t266, t124, 0, t229, 0, 0, -t228, 0, t268, t263, t143, t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, -t188, 0, 0, 0, t194, 0, -t193, 0, -t175, -t234, t151, pkin(4) * t151, 0, t194, 0, 0, t193, 0, -t175, -t138, t234, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, 0, -t189, 0, 0, 0, t193, 0, t194, 0, t234, -t175, t232, t141, 0, t193, 0, 0, -t194, 0, t234, t233, t175, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, t189, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t224, t236, 0, -t148, 0, 0, 0, qJDD(3), 0, 0, t191 + t218, 0, t255 + t225, t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t214, 0, 0, -t208, -t230, 0, 0, qJDD(3), 0, 0, t214, 0, 0, t157, t208, qJ(4) * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, 0, qJDD(3), 0, t208, 0, -t241, 0, 0, t214, 0, 0, -qJDD(3), 0, t208, t155, 0, pkin(3) * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t230, t241, 0, 0, 0, 0, 0, qJDD(3), 0, 0, t218, 0, t225, t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, 0, t214, 0, 0, t157, t208, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, 0, -t157, 0, t155, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214, 0, 0, qJDD(3), 0, -t208, -t155, 0, 0;];
m_new_reg  = t1;
