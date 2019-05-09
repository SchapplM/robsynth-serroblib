% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PRPP3
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PRPP3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP3_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_invdynm_fixb_reg2_snew_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:55:54
% EndTime: 2019-05-04 18:55:55
% DurationCPUTime: 0.80s
% Computational Cost: add. (889->137), mult. (1216->91), div. (0->0), fcn. (506->2), ass. (0->74)
t136 = (qJDD(2) * qJ(3));
t167 = (qJD(3) * qJD(2));
t183 = 2 * t136 + 2 * t167;
t182 = pkin(2) + pkin(3);
t139 = g(2) - qJDD(1);
t140 = sin(qJ(2));
t141 = cos(qJ(2));
t120 = t140 * g(1) - t141 * t139;
t121 = t141 * g(1) + t140 * t139;
t97 = t141 * t120 - t140 * t121;
t181 = pkin(1) * t97;
t180 = pkin(4) * t97;
t143 = qJD(2) ^ 2;
t166 = t140 * qJDD(2);
t122 = t141 * t143 + t166;
t179 = pkin(1) * t122;
t152 = -qJDD(3) + t120;
t145 = t143 * qJ(3) + t152;
t109 = t182 * qJDD(2) + t145;
t178 = pkin(3) * t109;
t177 = pkin(4) * t122;
t165 = t141 * qJDD(2);
t123 = -t140 * t143 + t165;
t176 = pkin(4) * t123;
t175 = qJ(1) * g(3);
t174 = t140 * g(3);
t173 = t141 * g(3);
t153 = t121 - (2 * t167);
t150 = -t136 + t153;
t106 = t182 * t143 + t150;
t172 = qJ(4) * t106;
t171 = qJ(4) * t109;
t170 = qJDD(2) * pkin(2);
t169 = qJDD(2) * pkin(3);
t138 = g(3) + qJDD(4);
t110 = t143 * pkin(2) + t150;
t112 = t145 + t170;
t168 = pkin(2) * t112 - qJ(3) * t110;
t88 = t140 * t110 - t141 * t112;
t164 = pkin(4) * t88 + qJ(3) * t173;
t163 = -t140 * t106 + t141 * t109;
t162 = -t141 * t110 - t140 * t112;
t161 = -t140 * t120 - t141 * t121;
t160 = pkin(1) * t88 - t168;
t159 = pkin(2) * t109 - qJ(3) * t106 + t178;
t126 = qJ(4) * t143 + t138;
t158 = qJ(4) * t166 + t141 * t126 - t177;
t157 = t173 - t177;
t114 = t174 + t176;
t82 = t141 * t106 + t140 * t109;
t91 = t182 * t138 + t172;
t95 = qJ(3) * t138 + t171;
t156 = -pkin(4) * t163 - t140 * t91 + t141 * t95;
t155 = -t121 + t179;
t119 = pkin(1) * t123;
t154 = t119 + t120;
t151 = -pkin(1) * t163 - t159;
t134 = 0.2e1 * t170;
t149 = t134 + t152;
t148 = qJ(1) * t123 - t155;
t147 = -t119 - 0.2e1 * t170 - t152;
t146 = qJ(1) * t122 + t154;
t144 = -qJDD(3) + t134 + t146;
t142 = pkin(1) * g(3);
t137 = qJ(4) * qJDD(2);
t133 = 0.2e1 * t169;
t113 = -t121 + t183;
t101 = qJ(4) * t165 - t140 * t126 - t176;
t94 = -(2 * t136) + t153 - t179;
t92 = pkin(4) * t161 + t142;
t90 = -t148 + t183;
t81 = pkin(4) * t162 + t142 + (pkin(2) * t141 + qJ(3) * t140) * g(3);
t80 = pkin(1) * t138 - pkin(4) * t82 + t140 * t95 + t141 * t91;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, 0, -t139, -g(3), -t175, 0, 0, t123, 0, -t122, 0, -t114, -t157, -t97, -t175 - t180, 0, t123, 0, 0, t122, 0, -t114, t88, t157, (-pkin(2) * t140 - qJ(1)) * g(3) + t164, 0, 0, -t123, 0, -t122, 0, t101, t158, t163, -qJ(1) * t138 + t156; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t122, 0, t123, 0, t157, -t114, t161, t92, 0, t122, 0, 0, -t123, 0, t157, t162, t114, t81, 0, 0, -t122, 0, t123, 0, t158, -t101, t82, t80; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, -t139, 0, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, qJDD(2), t146, t148, 0, -qJ(1) * t161 + t181, 0, 0, 0, qJDD(2), 0, 0, t144, 0, t90, -qJ(1) * t162 - t160, 0, 0, 0, 0, 0, qJDD(2), t133 + t144, t90, 0, qJ(1) * t82 - t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -g(3), 0, 0, 0, t123, 0, -t122, 0, -t114, -t157, -t97, -t180, 0, t123, 0, 0, t122, 0, -t114, t88, t157, -pkin(2) * t174 + t164, 0, 0, -t123, 0, -t122, 0, t101, t158, t163, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, -g(1), 0, 0, 0, 0, 0, 0, -qJDD(2), -t154, t155, 0, -t181, 0, 0, 0, -qJDD(2), 0, 0, t147, 0, t94, t160, 0, 0, 0, 0, 0, -qJDD(2), t147 - 0.2e1 * t169, t94, 0, t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t122, 0, t123, 0, t157, -t114, t161, t92, 0, t122, 0, 0, -t123, 0, t157, t162, t114, t81, 0, 0, -t122, 0, t123, 0, t158, -t101, t82, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t143, 0, 0, -g(3), -t120, 0, 0, qJDD(2), 0, 0, t143, 0, 0, -t112, g(3), qJ(3) * g(3), 0, 0, -qJDD(2), 0, -t143, 0, t137, t126, t109, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, 0, qJDD(2), 0, g(3), 0, -t121, 0, 0, t143, 0, 0, -qJDD(2), 0, g(3), -t110, 0, pkin(2) * g(3), 0, 0, -t143, 0, qJDD(2), 0, t126, -t137, t106, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t120, t121, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t149, 0, t113, t168, 0, 0, 0, 0, 0, qJDD(2), t133 + t149, t113, 0, t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t143, 0, 0, -t112, g(3), 0, 0, 0, -qJDD(2), 0, -t143, 0, t137, t126, t109, t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t112, 0, -t110, 0, 0, 0, 0, 0, 0, qJDD(2), t112 + t133, -t110, 0, t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, 0, 0, qJDD(2), 0, -g(3), t110, 0, 0, 0, 0, t143, 0, -qJDD(2), 0, -t126, t137, -t106, -pkin(3) * t138 - t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t143, 0, 0, t138, t109, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, 0, -qJDD(2), 0, -t138, 0, -t106, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), -t109, t106, 0, 0;];
m_new_reg  = t1;
