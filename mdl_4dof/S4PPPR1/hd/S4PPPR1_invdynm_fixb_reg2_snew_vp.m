% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PPPR1
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
%   pkin=[a2,a3,a4,d4,theta1]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PPPR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:35:37
% EndTime: 2019-05-04 18:35:39
% DurationCPUTime: 0.82s
% Computational Cost: add. (980->121), mult. (1329->96), div. (0->0), fcn. (1092->4), ass. (0->63)
t139 = sin(qJ(4));
t140 = cos(qJ(4));
t171 = qJD(4) ^ 2;
t123 = t140 * qJDD(4) - t139 * t171;
t124 = t139 * qJDD(4) + t140 * t171;
t137 = sin(pkin(5));
t138 = cos(pkin(5));
t104 = t138 * t123 - t137 * t124;
t136 = g(3) - qJDD(1);
t152 = -pkin(4) * t124 + t140 * t136;
t180 = pkin(2) * t124 - t152;
t173 = -pkin(4) * t123 - t139 * t136;
t97 = -pkin(2) * t123 + t173;
t182 = -qJ(1) * t104 + t137 * t180 + t138 * t97;
t103 = t137 * t123 + t138 * t124;
t181 = -qJ(1) * t103 + t137 * t97 - t138 * t180;
t125 = t137 * g(1) - t138 * g(2);
t121 = -qJDD(2) + t125;
t126 = t138 * g(1) + t137 * g(2);
t122 = -qJDD(3) + t126;
t101 = -t139 * t121 + t140 * t122;
t102 = -t140 * t121 - t139 * t122;
t143 = t139 * t101 + t140 * t102;
t92 = t140 * t101 - t139 * t102;
t179 = t137 * t92 + t138 * t143;
t178 = t137 * t143 - t138 * t92;
t108 = t138 * t121;
t112 = t137 * t126;
t174 = t112 - t108;
t133 = qJ(3) * t136;
t172 = (pkin(2) + pkin(4)) * t92 - t133;
t169 = pkin(2) * t143;
t168 = pkin(1) * t136;
t167 = pkin(2) * t121;
t163 = pkin(1) + qJ(3);
t161 = qJ(2) * t136;
t160 = qJ(3) * t121;
t159 = t137 * t121;
t127 = t137 * t136;
t156 = t138 * t136;
t88 = pkin(4) * t143;
t151 = pkin(3) * t136 + t88;
t150 = pkin(2) * t122 - t133;
t149 = -t137 * t122 + t108;
t148 = t138 * t125 - t112;
t113 = t138 * t126;
t147 = -t113 - t159;
t146 = -t137 * t125 - t113;
t144 = pkin(3) * t123 - t101;
t142 = t138 * t122 + t159;
t141 = -pkin(3) * t124 - t102;
t114 = pkin(1) * t121;
t107 = t161 - t167;
t106 = -t150 + t168;
t105 = -qJ(2) * t126 + t114;
t94 = -qJ(2) * t122 + t114 + t160;
t89 = pkin(3) * t92;
t87 = qJ(2) * t123 + t163 * t124 + t144;
t86 = -qJ(2) * t124 + t163 * t123 + t141;
t85 = t168 - t172;
t84 = t169 + t88 + (pkin(3) + qJ(2)) * t136;
t83 = -qJ(2) * t92 - t143 * t163 - t89;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t127, -t156, -t148, -qJ(1) * t148, 0, 0, 0, 0, 0, 0, t174, t127, t156, qJ(1) * t174 + (-t137 * pkin(1) + t138 * qJ(2)) * t136, 0, 0, 0, 0, 0, 0, t156, t149, -t127, -qJ(1) * t149 - t137 * t106 + t138 * t107, 0, 0, t103, 0, t104, 0, t181, t182, t179, qJ(1) * t179 - t137 * t85 + t138 * t84; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t156, -t127, t146, qJ(1) * t146, 0, 0, 0, 0, 0, 0, t147, -t156, t127, qJ(1) * t147 + (t138 * pkin(1) + t137 * qJ(2)) * t136, 0, 0, 0, 0, 0, 0, t127, t142, t156, -qJ(1) * t142 + t138 * t106 + t137 * t107, 0, 0, -t104, 0, t103, 0, -t182, t181, t178, qJ(1) * t178 + t137 * t84 + t138 * t85; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t125, t126, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t126, t105, 0, 0, 0, 0, 0, 0, -t122, 0, t121, t94, 0, 0, 0, 0, 0, qJDD(4), t87, t86, 0, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, -t125, 0, 0, 0, 0, 0, 0, 0, -t121, 0, t136, t161, 0, 0, 0, 0, 0, 0, t136, t121, 0, t107, 0, 0, t124, 0, t123, 0, -t180, t97, t143, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, -t126, 0, 0, 0, 0, 0, 0, 0, -t126, -t136, 0, t168, 0, 0, 0, 0, 0, 0, 0, t122, t136, t106, 0, 0, -t123, 0, t124, 0, -t97, -t180, -t92, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t126, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t126, t105, 0, 0, 0, 0, 0, 0, -t122, 0, t121, t94, 0, 0, 0, 0, 0, qJDD(4), t87, t86, 0, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t126, 0, 0, 0, 0, 0, 0, 0, -t122, 0, t121, t160, 0, 0, 0, 0, 0, qJDD(4), qJ(3) * t124 + t144, qJ(3) * t123 + t141, 0, -qJ(3) * t143 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, -t136, 0, 0, 0, 0, 0, 0, 0, -t136, -t121, 0, t167, 0, 0, -t124, 0, -t123, 0, t180, -t97, -t143, -t151 - t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, t136, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t136, t150, 0, 0, t123, 0, -t124, 0, t97, t180, t92, t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t136, 0, 0, 0, t123, 0, -t124, 0, t173, -t152, t92, pkin(4) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, 0, -t121, 0, 0, 0, 0, 0, 0, -qJDD(4), -t144, -t141, 0, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t121, 0, 0, 0, 0, t124, 0, t123, 0, t152, t173, t143, t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), 0, -t171, 0, 0, -t136, t101, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, 0, qJDD(4), 0, t136, 0, t102, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -t101, -t102, 0, 0;];
m_new_reg  = t1;
