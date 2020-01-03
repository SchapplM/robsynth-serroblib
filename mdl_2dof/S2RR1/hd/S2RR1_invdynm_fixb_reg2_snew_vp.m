% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S2RR1
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
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S2RR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynm_fixb_reg2_snew_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynm_fixb_reg2_snew_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynm_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:11
% EndTime: 2020-01-03 11:19:11
% DurationCPUTime: 0.29s
% Computational Cost: add. (298->86), mult. (741->105), div. (0->0), fcn. (480->4), ass. (0->70)
t148 = sin(qJ(1));
t150 = cos(qJ(1));
t135 = t148 * g(1) + t150 * g(3);
t171 = qJDD(1) * pkin(1);
t124 = -t135 - t171;
t147 = sin(qJ(2));
t149 = cos(qJ(2));
t121 = -t149 * g(2) + t147 * t124;
t122 = t147 * g(2) + t149 * t124;
t153 = t149 * t121 - t147 * t122;
t175 = pkin(1) * t153;
t106 = t147 * t121 + t149 * t122;
t174 = pkin(1) * t106;
t173 = t148 * g(2);
t172 = t150 * g(2);
t145 = t147 ^ 2;
t152 = qJD(1) ^ 2;
t170 = t145 * t152;
t146 = t149 ^ 2;
t169 = t146 * t152;
t136 = t150 * g(1) - t148 * g(3);
t125 = t152 * pkin(1) + t136;
t168 = t147 * t125;
t141 = t149 * t152 * t147;
t133 = qJDD(2) + t141;
t167 = t147 * t133;
t134 = qJDD(2) - t141;
t166 = t147 * t134;
t165 = t149 * t125;
t164 = t149 * t133;
t163 = t149 * t134;
t162 = -t145 - t146;
t161 = qJD(1) * qJD(2);
t160 = t147 * qJDD(1);
t159 = t148 * qJDD(1);
t158 = t149 * qJDD(1);
t157 = t150 * qJDD(1);
t142 = t147 * t161;
t156 = t149 * t161;
t155 = t148 * t141;
t154 = t150 * t141;
t151 = qJD(2) ^ 2;
t140 = -t151 - t169;
t139 = -t151 + t169;
t138 = -t151 - t170;
t137 = t151 - t170;
t132 = (t145 - t146) * t152;
t131 = -t148 * t152 + t157;
t130 = t150 * t152 + t159;
t129 = t142 - t158;
t128 = 0.2e1 * t142 - t158;
t127 = -0.2e1 * t156 - t160;
t126 = -t156 - t160;
t123 = t162 * t161;
t120 = t149 * t126 + t145 * t161;
t119 = -t147 * t129 + t146 * t161;
t118 = -t147 * t137 + t164;
t117 = t149 * t139 - t166;
t116 = -t149 * t137 - t167;
t115 = -t147 * t139 - t163;
t114 = (-t126 + t156) * t147;
t113 = (-t129 - t142) * t149;
t112 = -t147 * t127 + t149 * t128;
t111 = -t149 * t127 - t147 * t128;
t110 = -t165 + pkin(1) * (-t149 * t138 + t166);
t109 = -t165 - pkin(1) * (t149 * t140 - t167);
t108 = t168 - pkin(1) * (-t147 * t138 - t163);
t107 = -t168 + pkin(1) * (-t147 * t140 - t164);
t103 = -t162 * t171 - t106;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t130, 0, t131, 0, t172, -t173, -t150 * t135 + t148 * t136, 0, t148 * t120 - t154, t148 * t112 + t150 * t132, t148 * t118 - t147 * t157, t148 * t119 + t154, t148 * t117 - t149 * t157, t150 * qJDD(2) + t148 * t123, t148 * t107 - t150 * t121, t148 * t110 - t150 * t122, t148 * t153, t148 * t175; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, qJDD(1), -t136, t135, 0, 0, t114, t111, t116, t113, t115, 0, t109, t108, t103, -t174; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, t131, 0, -t130, 0, -t173, -t172, t148 * t135 + t150 * t136, 0, t150 * t120 + t155, t150 * t112 - t148 * t132, t150 * t118 + t147 * t159, t150 * t119 - t155, t150 * t117 + t148 * t158, -t148 * qJDD(2) + t150 * t123, t150 * t107 + t148 * t121, t150 * t110 + t148 * t122, t150 * t153, t150 * t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t152, 0, 0, -g(2), t136, 0, t120, t112, t118, t119, t117, t123, t107, t110, t153, t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, 0, qJDD(1), 0, g(2), 0, -t135, 0, -t141, t132, -t160, t141, -t158, qJDD(2), -t121, -t122, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t136, t135, 0, 0, t114, t111, t116, t113, t115, 0, t109, t108, t103, -t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, t128, t133, t156, t139, -t156, 0, -t125, t121, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, t127, t137, t129, t134, t142, t125, 0, t122, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t132, -t160, t141, -t158, qJDD(2), -t121, -t122, 0, 0;];
m_new_reg = t1;
