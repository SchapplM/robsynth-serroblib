% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:56
% EndTime: 2019-12-31 19:52:58
% DurationCPUTime: 0.65s
% Computational Cost: add. (1718->132), mult. (2260->137), div. (0->0), fcn. (1078->6), ass. (0->98)
t86 = qJD(1) + qJD(2);
t84 = t86 ^ 2;
t90 = sin(qJ(4));
t93 = cos(qJ(4));
t71 = t90 * t84 * t93;
t66 = qJDD(4) + t71;
t137 = t90 * t66;
t89 = t93 ^ 2;
t140 = t89 * t84;
t96 = qJD(4) ^ 2;
t70 = t96 + t140;
t38 = t93 * t70 + t137;
t125 = qJD(4) * t86;
t123 = t90 * t125;
t85 = qJDD(1) + qJDD(2);
t133 = t93 * t85;
t56 = -0.2e1 * t123 + t133;
t91 = sin(qJ(2));
t94 = cos(qJ(2));
t151 = pkin(1) * (t94 * t38 + t91 * t56);
t149 = -pkin(7) - pkin(2);
t150 = -qJ(3) * t56 + t149 * t38;
t127 = t93 * qJ(5);
t148 = pkin(4) * t90 - t127;
t146 = pkin(1) * (t91 * t84 - t94 * t85);
t145 = pkin(1) * (t94 * t84 + t91 * t85);
t143 = t85 * pkin(2);
t92 = sin(qJ(1));
t95 = cos(qJ(1));
t121 = t92 * g(1) - t95 * g(2);
t64 = qJDD(1) * pkin(1) + t121;
t115 = t95 * g(1) + t92 * g(2);
t65 = -qJD(1) ^ 2 * pkin(1) - t115;
t30 = t94 * t64 - t91 * t65;
t118 = qJDD(3) - t143 - t30;
t26 = -t84 * qJ(3) + t118;
t23 = -t85 * pkin(7) + t26;
t15 = t90 * g(3) + t93 * t23;
t16 = t93 * g(3) - t90 * t23;
t6 = t93 * t15 - t90 * t16;
t142 = t148 * t84;
t88 = t90 ^ 2;
t141 = t88 * t84;
t138 = t90 * t56;
t136 = t90 * t85;
t67 = qJDD(4) - t71;
t134 = t93 * t67;
t128 = t85 * qJ(3);
t31 = t91 * t64 + t94 * t65;
t109 = t84 * pkin(2) - t128 - t31;
t126 = (qJD(3) * t86);
t76 = 2 * t126;
t24 = -t109 + t76;
t132 = -pkin(2) * t26 + qJ(3) * t24;
t130 = t88 + t89;
t63 = t130 * t84;
t131 = t63 - t96;
t124 = qJD(5) * t93;
t122 = t93 * t125;
t120 = 0.2e1 * t128 + t76 + t31;
t107 = t84 * pkin(7) + t109;
t20 = -t107 + t76;
t119 = qJ(3) * t20 + t149 * t6;
t69 = -t96 - t141;
t37 = t90 * t69 + t134;
t53 = 0.2e1 * t122 + t136;
t117 = qJ(3) * t53 + t149 * t37;
t61 = t130 * t85;
t116 = -qJ(3) * t63 - t149 * t61;
t114 = pkin(4) * t93 + qJ(5) * t90;
t27 = t93 * t53 + t138;
t41 = -t93 * (-t96 + t141) + t137;
t111 = -t143 + t118;
t110 = t93 * t20 - t150;
t54 = -t122 - t136;
t55 = -t123 + t133;
t101 = -t54 * pkin(4) - t55 * qJ(5) - t107;
t97 = 0.2e1 * t86 * t124 - t101 - (2 * t126);
t108 = t93 * (-pkin(4) * t122 + (t56 - t123) * qJ(5) + t97) - pkin(4) * t138 + t150;
t106 = t90 * t20 + t117;
t105 = -t93 * t142 - qJDD(5) + t15;
t100 = qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t90 * t142 - t16;
t99 = -qJDD(4) * pkin(4) - t105;
t104 = -t90 * (t131 * pkin(4) + t100) + t93 * (t131 * qJ(5) + t99) + t116;
t103 = t116 - t6;
t12 = -t96 * pkin(4) + t100;
t13 = t96 * qJ(5) - t99;
t2 = t90 * t12 + t93 * t13;
t9 = t76 + (t114 * qJD(4) - 0.2e1 * t124) * t86 + t101;
t102 = t149 * t2 + (qJ(3) + t148) * t9;
t98 = -t53 * t127 - t90 * (-qJ(5) * t123 + (-t53 - t122) * pkin(4) + t97) + t117;
t62 = (-t88 + t89) * t84;
t42 = t134 - t90 * (t96 - t140);
t33 = (t55 - t123) * t93;
t32 = (-t54 + t122) * t90;
t29 = pkin(1) * (t94 * t61 - t91 * t63);
t19 = pkin(1) * (-t94 * t37 + t91 * t53);
t1 = [0, 0, 0, 0, 0, qJDD(1), t121, t115, 0, 0, 0, 0, 0, 0, 0, t85, t30 - t146, -t31 - t145, 0, pkin(1) * (t94 * t30 + t91 * t31), t85, 0, 0, 0, 0, 0, 0, t111 + t146, t120 + t145, pkin(1) * (t91 * t24 - t94 * t26) + t132, t33, -t27, t42, t32, -t41, 0, t106 + t19, t110 + t151, t103 + t29, pkin(1) * (t91 * t20 - t94 * t6) + t119, t33, t42, t27, 0, t41, t32, t19 + t98, t104 + t29, t108 - t151, pkin(1) * (-t94 * t2 + t91 * t9) + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t30, -t31, 0, 0, t85, 0, 0, 0, 0, 0, 0, t111, t120, t132, t33, -t27, t42, t32, -t41, 0, t106, t110, t103, t119, t33, t42, t27, 0, t41, t32, t98, t104, t108, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t84, t26, 0, 0, 0, 0, 0, 0, t37, -t38, -t61, t6, 0, 0, 0, 0, 0, 0, t37, -t61, t38, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t62, t133, -t71, -t136, qJDD(4), t15, t16, 0, 0, t71, t133, -t62, qJDD(4), t136, -t71, (t69 + t96) * qJ(5) + (qJDD(4) + t67) * pkin(4) + t105, -t114 * t85, qJ(5) * t66 + (t70 - t96) * pkin(4) + t100, pkin(4) * t13 + qJ(5) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t133, -t70, -t13;];
tauJ_reg = t1;
