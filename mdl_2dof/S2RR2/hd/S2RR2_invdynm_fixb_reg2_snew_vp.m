% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S2RR2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
%
% Output:
% m_new_reg [(3*3)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S2RR2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_invdynm_fixb_reg2_snew_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_invdynm_fixb_reg2_snew_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR2_invdynm_fixb_reg2_snew_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_invdynm_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:19:16
% EndTime: 2019-05-04 18:19:16
% DurationCPUTime: 0.28s
% Computational Cost: add. (298->83), mult. (741->104), div. (0->0), fcn. (480->4), ass. (0->70)
t151 = sin(qJ(1));
t153 = cos(qJ(1));
t137 = t153 * g(1) - t151 * g(3);
t170 = qJDD(1) * pkin(1);
t125 = -t137 + t170;
t150 = sin(qJ(2));
t152 = cos(qJ(2));
t121 = t152 * g(2) + t150 * t125;
t122 = -t150 * g(2) + t152 * t125;
t105 = t121 * t152 - t122 * t150;
t173 = pkin(1) * t105;
t172 = t151 * g(2);
t171 = t153 * g(2);
t148 = t150 ^ 2;
t155 = qJD(1) ^ 2;
t169 = t148 * t155;
t136 = t151 * g(1) + t153 * g(3);
t126 = t155 * pkin(1) + t136;
t168 = t150 * t126;
t142 = t152 * t155 * t150;
t134 = qJDD(2) + t142;
t167 = t150 * t134;
t135 = qJDD(2) - t142;
t166 = t150 * t135;
t165 = t152 * t126;
t164 = t152 * t135;
t149 = t152 ^ 2;
t163 = t148 + t149;
t162 = qJD(1) * qJD(2);
t161 = t151 * qJDD(1);
t145 = t152 * qJDD(1);
t160 = t153 * qJDD(1);
t159 = t150 * t162;
t143 = t152 * t162;
t158 = t150 * t121 + t152 * t122;
t157 = t151 * t142;
t156 = t153 * t142;
t154 = qJD(2) ^ 2;
t146 = t149 * t155;
t144 = t150 * qJDD(1);
t141 = -t146 - t154;
t140 = t146 - t154;
t139 = -t154 - t169;
t138 = t154 - t169;
t133 = -t146 + t169;
t132 = -t153 * t155 - t161;
t131 = -t151 * t155 + t160;
t130 = t145 - 0.2e1 * t159;
t129 = t145 - t159;
t128 = t144 + t143;
t127 = t144 + 0.2e1 * t143;
t124 = t152 * t134;
t123 = t163 * t162;
t120 = t152 * t128 - t148 * t162;
t119 = -t150 * t129 - t149 * t162;
t117 = -t150 * t138 + t124;
t116 = t152 * t140 - t166;
t115 = t152 * t138 + t167;
t114 = t150 * t140 + t164;
t113 = (t128 + t143) * t150;
t112 = (t129 - t159) * t152;
t111 = -t150 * t127 + t152 * t130;
t110 = t127 * t152 + t130 * t150;
t109 = -t165 - pkin(1) * (t139 * t152 - t166);
t108 = t165 + pkin(1) * (t141 * t152 - t167);
t107 = -t168 + pkin(1) * (-t139 * t150 - t164);
t106 = -t168 - pkin(1) * (t141 * t150 + t124);
t103 = pkin(1) * t158;
t102 = t163 * t170 + t158;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t131, 0, t132, 0, -t172, -t171, -t136 * t153 + t137 * t151, 0, t120 * t153 - t157, t111 * t153 + t133 * t151, t117 * t153 + t150 * t161, t119 * t153 + t157, t116 * t153 + t151 * t145, qJDD(2) * t151 + t123 * t153, t106 * t153 - t121 * t151, t109 * t153 - t122 * t151, t153 * t105, t153 * t173; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, qJDD(1), t136, t137, 0, 0, t113, t110, t115, t112, t114, 0, t108, t107, t102, t103; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, t132, 0, -t131, 0, -t171, t172, t136 * t151 + t137 * t153, 0, -t120 * t151 - t156, -t111 * t151 + t133 * t153, -t117 * t151 + t150 * t160, -t119 * t151 + t156, -t116 * t151 + t152 * t160, qJDD(2) * t153 - t123 * t151, -t106 * t151 - t121 * t153, -t109 * t151 - t122 * t153, -t151 * t105, -t151 * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t155, 0, 0, -g(2), -t136, 0, t120, t111, t117, t119, t116, t123, t106, t109, t105, t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, 0, qJDD(1), 0, g(2), 0, -t137, 0, t142, -t133, -t144, -t142, -t145, -qJDD(2), t121, t122, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t136, t137, 0, 0, t113, t110, t115, t112, t114, 0, t108, t107, t102, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t130, t134, -t143, t140, t143, 0, -t126, t121, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t127, t138, t129, t135, -t159, t126, 0, t122, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, t133, t144, t142, t145, qJDD(2), -t121, -t122, 0, 0;];
m_new_reg  = t1;
