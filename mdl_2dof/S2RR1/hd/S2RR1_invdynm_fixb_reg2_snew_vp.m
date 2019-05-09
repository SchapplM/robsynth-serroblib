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
% Datum: 2019-05-04 18:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
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
% StartTime: 2019-05-04 18:17:53
% EndTime: 2019-05-04 18:17:53
% DurationCPUTime: 0.30s
% Computational Cost: add. (300->88), mult. (741->105), div. (0->0), fcn. (480->4), ass. (0->70)
t149 = sin(qJ(1));
t151 = cos(qJ(1));
t134 = t149 * g(1) + t151 * g(3);
t172 = qJDD(1) * pkin(1);
t123 = t134 - t172;
t148 = sin(qJ(2));
t150 = cos(qJ(2));
t120 = -t150 * g(2) + t148 * t123;
t121 = t148 * g(2) + t150 * t123;
t154 = t150 * t120 - t148 * t121;
t175 = pkin(1) * t154;
t105 = t148 * t120 + t150 * t121;
t174 = pkin(1) * t105;
t173 = t151 * g(2);
t145 = t148 ^ 2;
t153 = qJD(1) ^ 2;
t171 = t145 * t153;
t146 = t150 ^ 2;
t170 = t146 * t153;
t135 = t151 * g(1) - t149 * g(3);
t124 = -t153 * pkin(1) + t135;
t169 = t148 * t124;
t140 = t150 * t153 * t148;
t132 = qJDD(2) + t140;
t168 = t148 * t132;
t133 = qJDD(2) - t140;
t167 = t148 * t133;
t166 = t150 * t124;
t165 = t150 * t132;
t164 = t150 * t133;
t163 = -t145 - t146;
t162 = qJD(1) * qJD(2);
t161 = t148 * qJDD(1);
t160 = t149 * qJDD(1);
t159 = t150 * qJDD(1);
t158 = t151 * qJDD(1);
t141 = t148 * t162;
t157 = t150 * t162;
t156 = t149 * t140;
t155 = t151 * t140;
t152 = qJD(2) ^ 2;
t143 = t149 * g(2);
t139 = -t152 - t170;
t138 = -t152 + t170;
t137 = -t152 - t171;
t136 = t152 - t171;
t131 = (t145 - t146) * t153;
t130 = t149 * t153 - t158;
t129 = t151 * t153 + t160;
t128 = t141 - t159;
t127 = 0.2e1 * t141 - t159;
t126 = -0.2e1 * t157 - t161;
t125 = -t157 - t161;
t122 = t163 * t162;
t119 = t150 * t125 + t145 * t162;
t118 = -t148 * t128 + t146 * t162;
t117 = -t148 * t136 + t165;
t116 = t150 * t138 - t167;
t115 = -t150 * t136 - t168;
t114 = -t148 * t138 - t164;
t113 = (-t125 + t157) * t148;
t112 = (-t128 - t141) * t150;
t111 = -t148 * t126 + t150 * t127;
t110 = -t150 * t126 - t148 * t127;
t109 = t166 - pkin(1) * (t150 * t139 - t168);
t108 = t166 + pkin(1) * (-t150 * t137 + t167);
t107 = t169 + pkin(1) * (-t148 * t139 - t165);
t106 = -t169 - pkin(1) * (-t148 * t137 - t164);
t102 = -t163 * t172 - t105;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, -t129, 0, t130, 0, -t173, t143, -t151 * t134 + t149 * t135, 0, -t149 * t119 + t155, -t149 * t111 - t151 * t131, -t149 * t117 + t148 * t158, -t149 * t118 - t155, -t149 * t116 + t150 * t158, -t151 * qJDD(2) - t149 * t122, -t149 * t107 + t151 * t120, -t149 * t108 + t151 * t121, -t149 * t154, -t149 * t175; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, qJDD(1), t135, -t134, 0, 0, t113, t110, t115, t112, t114, 0, t109, t106, t102, -t174; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, t130, 0, t129, 0, t143, t173, t149 * t134 + t151 * t135, 0, -t151 * t119 - t156, -t151 * t111 + t149 * t131, -t151 * t117 - t148 * t160, -t151 * t118 + t156, -t151 * t116 - t149 * t159, t149 * qJDD(2) - t151 * t122, -t151 * t107 - t149 * t120, -t151 * t108 - t149 * t121, -t151 * t154, -t151 * t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t153, 0, 0, -g(2), -t135, 0, t119, t111, t117, t118, t116, t122, t107, t108, t154, t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, 0, qJDD(1), 0, g(2), 0, t134, 0, -t140, t131, -t161, t140, -t159, qJDD(2), -t120, -t121, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t135, -t134, 0, 0, t113, t110, t115, t112, t114, 0, t109, t106, t102, -t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t127, t132, t157, t138, -t157, 0, t124, t120, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t126, t136, t128, t133, t141, -t124, 0, t121, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, t131, -t161, t140, -t159, qJDD(2), -t120, -t121, 0, 0;];
m_new_reg  = t1;
