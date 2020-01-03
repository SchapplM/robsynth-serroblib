% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRR2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:20
% EndTime: 2019-12-31 17:23:22
% DurationCPUTime: 0.79s
% Computational Cost: add. (921->166), mult. (1409->236), div. (0->0), fcn. (934->12), ass. (0->128)
t92 = sin(qJ(3));
t95 = cos(qJ(4));
t145 = t95 * t92;
t91 = sin(qJ(4));
t96 = cos(qJ(3));
t40 = t91 * t96 + t145;
t86 = qJD(1) + qJD(2);
t30 = t40 * t86;
t90 = qJ(1) + qJ(2);
t80 = cos(t90);
t161 = g(2) * t80;
t138 = pkin(1) * qJD(1);
t93 = sin(qJ(2));
t125 = t93 * t138;
t97 = cos(qJ(2));
t158 = t97 * pkin(1);
t140 = -qJD(2) * t125 + qJDD(1) * t158;
t84 = qJDD(1) + qJDD(2);
t160 = t84 * pkin(2);
t32 = -t140 - t160;
t163 = t32 + t161;
t78 = sin(t90);
t141 = g(1) * t80 + g(2) * t78;
t85 = qJD(3) + qJD(4);
t162 = -pkin(6) - pkin(7);
t70 = g(1) * t78;
t159 = t86 * pkin(2);
t73 = t93 * pkin(1) + pkin(6);
t157 = -pkin(7) - t73;
t144 = t95 * t96;
t129 = t86 * t144;
t149 = t91 * t92;
t130 = t86 * t149;
t28 = -t129 + t130;
t156 = t30 * t28;
t89 = qJ(3) + qJ(4);
t77 = sin(t89);
t155 = t77 * t78;
t154 = t77 * t80;
t79 = cos(t89);
t153 = t78 * t79;
t152 = t79 * t80;
t151 = t85 * t97;
t150 = t86 * t92;
t147 = t92 * t84;
t44 = t86 * pkin(6) + t125;
t121 = pkin(7) * t86 + t44;
t26 = t121 * t96;
t146 = t95 * t26;
t143 = t96 * t84;
t134 = qJD(3) * t92;
t137 = qJD(1) * t97;
t124 = pkin(1) * t137;
t45 = -t124 - t159;
t142 = t45 * t134 + t96 * t70;
t87 = t92 ^ 2;
t139 = -t96 ^ 2 + t87;
t136 = qJD(2) * t93;
t135 = qJD(2) * t97;
t133 = qJD(3) * t96;
t132 = qJD(4) * t91;
t131 = qJDD(1) * t93;
t128 = t45 * t133 + t163 * t92;
t127 = pkin(1) * t135;
t126 = pkin(3) * t134;
t123 = t86 * t136;
t122 = t86 * t133;
t75 = -t96 * pkin(3) - pkin(2);
t119 = qJD(3) * t162;
t118 = qJD(3) * t157;
t117 = t86 * t125;
t116 = -t95 * t143 + t91 * t147;
t25 = t121 * t92;
t22 = qJD(3) * pkin(3) - t25;
t115 = -t91 * t22 - t146;
t37 = t157 * t92;
t81 = t96 * pkin(7);
t38 = t96 * t73 + t81;
t114 = t95 * t37 - t91 * t38;
t113 = t91 * t37 + t95 * t38;
t63 = t162 * t92;
t64 = t96 * pkin(6) + t81;
t112 = t95 * t63 - t91 * t64;
t111 = t91 * t63 + t95 * t64;
t39 = -t144 + t149;
t110 = t140 + t70 - t161;
t107 = t86 * t134 - t143;
t17 = t107 * pkin(3) + t32;
t18 = t85 * t39;
t31 = t75 * t86 - t124;
t109 = -g(1) * t155 + g(2) * t154 + t17 * t40 - t31 * t18;
t19 = t85 * t40;
t108 = g(1) * t153 - g(2) * t152 + t17 * t39 + t31 * t19;
t33 = t84 * pkin(6) + (qJD(1) * t135 + t131) * pkin(1);
t106 = -t45 * t86 + t141 - t33;
t99 = qJD(3) ^ 2;
t105 = pkin(6) * t99 - t117 - t160;
t74 = -pkin(2) - t158;
t104 = pkin(1) * t123 + t73 * t99 + t74 * t84;
t5 = qJD(4) * t129 + t95 * t122 - t85 * t130 + t91 * t143 + t84 * t145;
t103 = -pkin(6) * qJDD(3) + (t124 - t159) * qJD(3);
t102 = -qJDD(3) * t73 + (t74 * t86 - t127) * qJD(3);
t9 = -t44 * t133 + qJDD(3) * pkin(3) - t92 * t33 + (-t122 - t147) * pkin(7);
t101 = t31 * t28 + t26 * t132 + g(2) * t153 + g(1) * t152 + g(3) * t77 + (-t26 * t85 - t9) * t91;
t10 = -t107 * pkin(7) - t44 * t134 + t96 * t33;
t100 = g(1) * t154 + g(2) * t155 - g(3) * t79 + t115 * qJD(4) - t91 * t10 - t31 * t30 + t95 * t9;
t98 = cos(qJ(1));
t94 = sin(qJ(1));
t83 = qJDD(3) + qJDD(4);
t82 = t86 ^ 2;
t56 = t75 - t158;
t55 = qJDD(3) * t96 - t99 * t92;
t54 = qJDD(3) * t92 + t99 * t96;
t43 = pkin(1) * t136 + t126;
t42 = t96 * t119;
t41 = t92 * t119;
t34 = 0.2e1 * t92 * t122 + t87 * t84;
t24 = t96 * t118 - t92 * t127;
t23 = t92 * t118 + t96 * t127;
t20 = -0.2e1 * t139 * t86 * qJD(3) + 0.2e1 * t92 * t143;
t12 = -t19 * t85 - t39 * t83;
t11 = -t18 * t85 + t40 * t83;
t8 = -t28 ^ 2 + t30 ^ 2;
t6 = t19 * t86 + t116;
t3 = t28 * t85 + t5;
t2 = -t30 * t18 + t5 * t40;
t1 = t18 * t28 - t30 * t19 - t5 * t39 - t40 * t6;
t4 = [qJDD(1), g(1) * t94 - g(2) * t98, g(1) * t98 + g(2) * t94, t84, (t84 * t97 - t123) * pkin(1) + t110, ((-qJDD(1) - t84) * t93 + (-qJD(1) - t86) * t135) * pkin(1) + t141, t34, t20, t54, t55, 0, t102 * t92 + (-t104 - t163) * t96 + t142, t102 * t96 + (t104 - t70) * t92 + t128, t2, t1, t11, t12, 0, t43 * t28 + t56 * t6 + (-t113 * qJD(4) - t91 * t23 + t95 * t24) * t85 + t114 * t83 + t108, t43 * t30 + t56 * t5 - (t114 * qJD(4) + t95 * t23 + t91 * t24) * t85 - t113 * t83 + t109; 0, 0, 0, t84, t110 + t117, (-t131 + (-qJD(2) + t86) * t137) * pkin(1) + t141, t34, t20, t54, t55, 0, t103 * t92 + (-t105 - t163) * t96 + t142, t103 * t96 + (t105 - t70) * t92 + t128, t2, t1, t11, t12, 0, t28 * t126 + t75 * t6 + (-t111 * qJD(4) - t91 * t41 + t95 * t42) * t85 + t112 * t83 + (t40 * t151 - t93 * t28) * t138 + t108, t30 * t126 + t75 * t5 - (t112 * qJD(4) + t95 * t41 + t91 * t42) * t85 - t111 * t83 + (-t39 * t151 - t93 * t30) * t138 + t109; 0, 0, 0, 0, 0, 0, -t92 * t82 * t96, t139 * t82, t147, t143, qJDD(3), -g(3) * t96 + t106 * t92, g(3) * t92 + t106 * t96, t156, t8, t3, -t116, t83, -(t91 * t25 - t146) * t85 + (-t85 * t132 - t28 * t150 + t95 * t83) * pkin(3) + t100, (-qJD(4) * t22 - t25 * t85 - t10) * t95 + (-qJD(4) * t95 * t85 - t30 * t150 - t91 * t83) * pkin(3) + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t8, t3, -t116, t83, -t115 * t85 + t100, (-t10 + (-qJD(4) + t85) * t22) * t95 + t101;];
tau_reg = t4;
