% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PPRP3
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
%   pkin=[a2,a3,a4,d3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PPRP3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP3_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_invdynm_fixb_reg2_snew_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:46:16
% EndTime: 2019-05-04 18:46:17
% DurationCPUTime: 0.78s
% Computational Cost: add. (813->128), mult. (961->73), div. (0->0), fcn. (572->2), ass. (0->67)
t126 = cos(qJ(3));
t127 = qJD(3) ^ 2;
t125 = sin(qJ(3));
t147 = t125 * qJDD(3);
t106 = t126 * t127 + t147;
t101 = pkin(4) * t106;
t120 = g(3) + qJDD(4);
t112 = qJ(4) * t127 + t120;
t145 = qJ(4) * t147 + t126 * t112 + t101;
t160 = pkin(1) * t106;
t168 = t160 + t145;
t146 = t126 * qJDD(3);
t107 = -t125 * t127 + t146;
t104 = pkin(1) * t107;
t100 = pkin(4) * t107;
t165 = -qJ(4) * t146 + t125 * t112 - t100;
t78 = t104 - t165;
t166 = t125 * g(3) - t100;
t87 = t104 - t166;
t95 = t126 * g(3) + t101;
t167 = t160 + t95;
t164 = pkin(2) * g(3);
t163 = -pkin(1) - pkin(4);
t123 = qJDD(3) * pkin(3);
t121 = g(2) - qJDD(1);
t122 = g(1) - qJDD(2);
t139 = t125 * t121 - t126 * t122;
t94 = -t123 - t139;
t153 = t126 * t94;
t98 = -t126 * t121 - t125 * t122;
t96 = -t127 * pkin(3) + t98;
t73 = t125 * t96 - t153;
t162 = pkin(2) * t73;
t83 = t125 * t98 + t126 * t139;
t161 = pkin(2) * t83;
t159 = pkin(1) * t121;
t158 = pkin(2) * t106;
t157 = qJ(1) * g(3);
t156 = pkin(2) + qJ(1);
t136 = -t125 * t139 + t126 * t98;
t155 = qJ(2) * t136;
t154 = t125 * t94;
t152 = qJ(2) * t121;
t148 = qJ(4) * qJDD(3);
t81 = pkin(4) * t136;
t143 = pkin(1) * t136 + t81;
t137 = t126 * t96 + t154;
t93 = pkin(3) * t94;
t140 = qJ(2) * t137 + t93;
t138 = -pkin(2) * t107 - t139;
t135 = qJ(2) * t106 - t138;
t92 = -pkin(3) * t120 + qJ(4) * t96;
t134 = qJ(4) * t153 - t125 * t92;
t132 = -pkin(2) * t120 + pkin(4) * t137 + qJ(4) * t154 + t126 * t92;
t131 = -qJ(2) * t107 + t98;
t130 = -qJ(1) * t107 - t135;
t129 = -pkin(1) * t137 - t132;
t124 = qJ(2) * g(3);
t119 = -0.2e1 * t123;
t118 = 0.2e1 * t123;
t111 = pkin(1) * t122 + t124;
t91 = t98 + t158;
t80 = -t131 - t158;
t70 = t106 * t156 + t131;
t69 = t163 * t83 + t124;
t68 = qJ(2) * t120 + t163 * t73 + t134;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t121, 0, -g(3), -t157, 0, 0, 0, 0, 0, 0, -g(3), t121, 0, -t157 - t159, 0, 0, t106, 0, t107, 0, -t167, -t87, t136, -g(3) * t156 + t143, 0, 0, t106, 0, t107, 0, -t168, -t78, t137, -qJ(1) * t120 - t129; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(1), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, g(3), t111, 0, 0, t107, 0, -t106, 0, -t87, t167, -t83, t69, 0, 0, t107, 0, -t106, 0, -t78, t168, -t73, t68; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, 0, t122, 0, -t121, qJ(1) * t122 - t152, 0, 0, 0, 0, 0, -qJDD(3), t130, t70, 0, -t156 * t83 + t155, 0, 0, 0, 0, 0, -qJDD(3), t119 + t130, t70, 0, -t156 * t73 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -g(1), 0, 0, 0, 0, 0, 0, 0, -t122, 0, t121, t152, 0, 0, 0, 0, 0, qJDD(3), t135, t80, 0, -t155 + t161, 0, 0, 0, 0, 0, qJDD(3), t118 + t135, t80, 0, -t140 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, g(3), 0, 0, 0, 0, 0, 0, 0, g(3), -t121, 0, t159, 0, 0, -t106, 0, -t107, 0, t167, t87, -t136, -t143 + t164, 0, 0, -t106, 0, -t107, 0, t168, t78, -t137, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, g(3), t111, 0, 0, t107, 0, -t106, 0, -t87, t167, -t83, t69, 0, 0, t107, 0, -t106, 0, -t78, t168, -t73, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, g(3), 0, 0, 0, t107, 0, -t106, 0, t166, t95, -t83, -pkin(4) * t83, 0, 0, t107, 0, -t106, 0, t165, t145, -t73, -pkin(4) * t73 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, 0, -t121, 0, 0, 0, 0, 0, 0, -qJDD(3), t138, t91, 0, -t161, 0, 0, 0, 0, 0, -qJDD(3), t119 + t138, t91, 0, t93 - t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), t121, 0, 0, 0, 0, t106, 0, t107, 0, -t95, t166, t136, t81 - t164, 0, 0, t106, 0, t107, 0, -t145, t165, t137, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t127, 0, 0, g(3), -t139, 0, 0, 0, qJDD(3), 0, -t127, 0, -t148, t112, t94, qJ(4) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, 0, qJDD(3), 0, -g(3), 0, t98, 0, 0, 0, t127, 0, qJDD(3), 0, -t112, -t148, t96, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t139, -t98, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t118 + t139, -t98, 0, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t127, 0, 0, t120, t94, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, 0, qJDD(3), 0, -t120, 0, t96, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t94, -t96, 0, 0;];
m_new_reg  = t1;
