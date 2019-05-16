% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S3RPP1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
%
% Output:
% m_new_reg [(3*4)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S3RPP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_invdynm_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_invdynm_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_invdynm_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:27:38
% EndTime: 2019-05-04 18:27:39
% DurationCPUTime: 0.42s
% Computational Cost: add. (419->92), mult. (671->70), div. (0->0), fcn. (242->2), ass. (0->54)
t139 = 2 * qJD(1);
t109 = sin(qJ(1));
t110 = cos(qJ(1));
t97 = t110 * g(1) + t109 * g(2);
t118 = (qJD(2) * t139) - t97;
t96 = t109 * g(1) - g(2) * t110;
t122 = -qJDD(2) + t96;
t115 = (qJD(3) * t139) + t122;
t112 = qJD(1) ^ 2;
t132 = t112 * qJ(2);
t113 = t115 + t132;
t134 = pkin(1) + qJ(3);
t80 = qJDD(1) * t134 + t113;
t138 = pkin(2) * t80;
t129 = t110 * qJDD(1);
t94 = -t109 * t112 + t129;
t137 = pkin(3) * t94;
t130 = t109 * qJDD(1);
t95 = t110 * t112 + t130;
t136 = pkin(3) * t95;
t135 = t112 * pkin(1);
t133 = qJDD(1) * pkin(1);
t106 = qJDD(1) * qJ(2);
t114 = -t106 - t118;
t81 = t112 * t134 - qJDD(3) + t114;
t128 = -t109 * t80 - t110 * t81;
t84 = t114 + t135;
t86 = t122 + t132 + t133;
t127 = -t109 * t86 - t110 * t84;
t126 = -t109 * t96 - t110 * t97;
t125 = pkin(2) * t81 - qJ(3) * g(3);
t124 = g(3) * t109 + t137;
t123 = g(3) * t110 - t136;
t121 = t109 * t81 - t110 * t80;
t120 = t109 * t84 - t110 * t86;
t119 = t109 * t97 - t110 * t96;
t100 = pkin(2) * t112 - g(3);
t117 = -pkin(2) * t130 - t100 * t110 - t136;
t116 = qJDD(3) + t118;
t111 = pkin(1) * g(3);
t108 = qJ(2) * g(3);
t107 = pkin(2) * qJDD(1);
t102 = 0.2e1 * t106;
t101 = 0.2e1 * qJDD(1) * qJ(3);
t92 = -t122 - 0.2e1 * t133;
t91 = t102 + t118;
t85 = t102 + t116;
t83 = t101 + t115 + 0.2e1 * t133;
t79 = -pkin(2) * t129 + t100 * t109 - t137;
t77 = t108 - t138;
t76 = t111 - t125;
t75 = pkin(1) * t86 - qJ(2) * t84;
t74 = -qJ(2) * t81 + t134 * t80;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t94, 0, -t95, 0, -t124, -t123, t119, pkin(3) * t119, 0, -t94, t95, 0, 0, 0, t120, t124, t123, pkin(3) * t120 + (-pkin(1) * t109 + qJ(2) * t110) * g(3), 0, t95, t94, 0, 0, 0, t121, t117, t79, pkin(3) * t121 - t109 * t76 + t110 * t77; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t95, 0, t94, 0, t123, -t124, t126, pkin(3) * t126, 0, -t95, -t94, 0, 0, 0, t127, -t123, t124, pkin(3) * t127 + (pkin(1) * t110 + qJ(2) * t109) * g(3), 0, -t94, t95, 0, 0, 0, t128, -t79, t117, pkin(3) * t128 + t109 * t77 + t110 * t76; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t96, t97, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t92, t91, t75, qJDD(1), 0, 0, 0, 0, 0, 0, t85, t83, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t112, 0, 0, -g(3), -t96, 0, 0, -qJDD(1), t112, 0, 0, 0, -t86, 0, g(3), t108, 0, t112, qJDD(1), 0, 0, 0, -t80, -t100, -t107, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, 0, qJDD(1), 0, g(3), 0, -t97, 0, 0, -t112, -qJDD(1), 0, 0, 0, -t84, -g(3), 0, t111, 0, -qJDD(1), t112, 0, 0, 0, -t81, t107, -t100, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t96, t97, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t92, t91, t75, qJDD(1), 0, 0, 0, 0, 0, 0, t85, t83, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t86, -t84, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t106 + t116 - t135, t101 + t113 + t133, qJ(3) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t112, 0, 0, 0, t86, 0, -g(3), 0, 0, -t112, -qJDD(1), 0, 0, 0, t80, t100, t107, t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, qJDD(1), 0, 0, 0, t84, g(3), 0, 0, 0, qJDD(1), -t112, 0, 0, 0, t81, -t107, t100, t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t81, t80, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t112, 0, 0, 0, t81, 0, -g(3), 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, qJDD(1), 0, 0, 0, -t80, g(3), 0, 0;];
m_new_reg  = t1;
