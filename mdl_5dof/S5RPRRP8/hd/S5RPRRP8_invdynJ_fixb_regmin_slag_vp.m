% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:32
% EndTime: 2019-12-31 18:47:34
% DurationCPUTime: 0.98s
% Computational Cost: add. (1436->192), mult. (2071->214), div. (0->0), fcn. (1102->6), ass. (0->119)
t147 = sin(qJ(1));
t148 = cos(qJ(1));
t67 = sin(qJ(3));
t69 = cos(qJ(3));
t32 = -t147 * t67 - t148 * t69;
t33 = -t147 * t69 + t148 * t67;
t94 = g(1) * t32 + g(2) * t33;
t110 = qJD(1) - qJD(3);
t167 = t110 ^ 2;
t62 = qJDD(1) - qJDD(3);
t66 = sin(qJ(4));
t68 = cos(qJ(4));
t92 = -t68 * pkin(4) - t66 * qJ(5);
t88 = pkin(3) - t92;
t166 = t88 * t62;
t164 = t110 * t88;
t117 = qJ(2) * qJD(1);
t70 = -pkin(1) - pkin(2);
t41 = t70 * qJD(1) + qJD(2);
t24 = -t67 * t117 + t69 * t41;
t10 = -t24 + t164;
t165 = t110 * t10;
t163 = qJD(3) * t117 - t70 * qJDD(1) - qJDD(2);
t162 = qJD(4) * t110;
t112 = qJD(1) * qJD(2);
t113 = qJ(2) * qJDD(1);
t161 = qJD(3) * t41 + t112 + t113;
t160 = t69 * t167;
t116 = qJD(4) * qJ(5);
t25 = t69 * t117 + t67 * t41;
t18 = -pkin(7) * t110 + t25;
t136 = t68 * t18;
t13 = t116 + t136;
t120 = qJDD(4) * pkin(4);
t122 = qJD(4) * t68;
t83 = t161 * t69 - t163 * t67;
t8 = -t62 * pkin(7) + t83;
t6 = t66 * t8;
t3 = t18 * t122 + qJDD(5) - t120 + t6;
t159 = -t13 * qJD(4) + t3;
t130 = t69 * qJ(2) + t67 * t70;
t100 = -qJD(4) * pkin(4) + qJD(5);
t137 = t66 * t18;
t12 = t100 + t137;
t111 = qJDD(4) * qJ(5);
t7 = t68 * t8;
t2 = t111 + t7 + (qJD(5) - t137) * qJD(4);
t151 = t2 * t68;
t77 = (t12 * t68 - t13 * t66) * qJD(4) + t3 * t66 + t151;
t73 = t77 + t94;
t103 = -t67 * qJ(2) + t69 * t70;
t89 = t12 * t66 + t13 * t68;
t95 = g(1) * t33 - g(2) * t32;
t158 = g(3) * t68 - t94 * t66;
t35 = -pkin(7) + t130;
t115 = qJDD(4) * t35;
t20 = -t103 + t88;
t22 = t69 * qJD(2) + t103 * qJD(3);
t157 = (-t110 * t20 - t10 - t22) * qJD(4) - t115;
t91 = pkin(4) * t66 - qJ(5) * t68;
t26 = t91 * qJD(4) - t66 * qJD(5);
t156 = -qJDD(4) * t67 + 0.2e1 * t69 * t162;
t150 = t62 * pkin(3);
t149 = t110 * pkin(3);
t144 = t22 * t110;
t23 = t67 * qJD(2) + t130 * qJD(3);
t143 = t23 * t110;
t142 = t24 * t110;
t141 = t25 * t110;
t140 = t26 * t110;
t138 = t110 * t66;
t135 = t68 * t62;
t134 = t69 * t62;
t132 = t66 * t24 * qJD(4) - t68 * t141;
t131 = t26 - t25;
t129 = t148 * pkin(1) + t147 * qJ(2);
t128 = g(1) * t147 - g(2) * t148;
t64 = t66 ^ 2;
t65 = t68 ^ 2;
t127 = t64 - t65;
t126 = t64 + t65;
t125 = pkin(1) * qJDD(1);
t124 = pkin(7) * qJDD(4);
t109 = t66 * t167 * t68;
t108 = 0.2e1 * t112;
t17 = -t24 + t149;
t107 = t17 + t149;
t106 = t126 * t62;
t104 = t10 + t164;
t101 = -t6 + t158;
t98 = qJDD(2) - t125;
t96 = -t147 * pkin(1) + t148 * qJ(2);
t87 = t161 * t67 + t163 * t69;
t71 = qJD(4) ^ 2;
t86 = pkin(7) * t71 - t95;
t85 = -t35 * t71 - t95;
t84 = g(1) * t148 + g(2) * t147;
t34 = pkin(3) - t103;
t82 = -t115 + (-t110 * t34 - t17 - t22) * qJD(4);
t9 = t87 + t150;
t81 = -t86 - t9 - t150;
t80 = -t134 + (-t71 - t167) * t67;
t1 = -t140 + t87 + t166;
t79 = -t1 - t86 - t166;
t78 = t87 - t95;
t11 = t23 - t26;
t76 = t11 * t110 + t20 * t62 + t1 + t85;
t75 = -t34 * t62 - t143 - t85 - t9;
t74 = t83 + t94;
t72 = qJD(1) ^ 2;
t43 = t66 * t62;
t39 = qJDD(4) * t68 - t71 * t66;
t38 = qJDD(4) * t66 + t71 * t68;
t29 = t91 * t110;
t21 = -0.2e1 * t122 * t138 - t64 * t62;
t15 = t127 * t162 - t66 * t135;
t5 = t156 * t66 + t80 * t68;
t4 = -t156 * t68 + t80 * t66;
t14 = [qJDD(1), t128, t84, -qJDD(2) + 0.2e1 * t125 + t128, t108 - t84 + 0.2e1 * t113, -t98 * pkin(1) - g(1) * t96 - g(2) * t129 + (t108 + t113) * qJ(2), t62, -t103 * t62 + t143 + t78, t130 * t62 + t144 + t74, -t21, -0.2e1 * t15, -t38, -t39, 0, t82 * t66 - t75 * t68, t75 * t66 + t82 * t68, t157 * t66 + t76 * t68, -t35 * t106 - t126 * t144 - t73, -t157 * t68 + t76 * t66, t1 * t20 + t10 * t11 - g(1) * (-t147 * pkin(2) + t32 * pkin(7) + t88 * t33 + t96) - g(2) * (t148 * pkin(2) + t33 * pkin(7) - t88 * t32 + t129) + t89 * t22 + (t12 * t122 + t159 * t66 + t151) * t35; 0, 0, 0, -qJDD(1), -t72, -t72 * qJ(2) - t128 + t98, 0, -t167 * t67 - t134, t67 * t62 - t160, 0, 0, 0, 0, 0, t5, -t4, t5, -t67 * t106 + t126 * t160, t4, (-t110 * t89 - t1) * t69 + (t77 - t165) * t67 - t128; 0, 0, 0, 0, 0, 0, -t62, -t78 - t141, -t74 - t142, t21, 0.2e1 * t15, t38, t39, 0, (t107 * qJD(4) - t124) * t66 + t81 * t68 + t132, (-t124 + (t107 + t24) * qJD(4)) * t68 + (-t81 + t141) * t66, (t104 * qJD(4) - t124) * t66 + (t79 + t140) * t68 + t132, -pkin(7) * t106 + t126 * t142 + t73, (t124 + (-t104 - t24) * qJD(4)) * t68 + (t110 * t131 + t79) * t66, t73 * pkin(7) + t131 * t10 - t89 * t24 + (-t1 + t95) * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, t127 * t167, -t43, -t135, qJDD(4), t17 * t138 + t101, -g(3) * t66 - t7 + (t110 * t17 - t94) * t68, 0.2e1 * t120 - qJDD(5) - (-t10 * t66 - t29 * t68) * t110 + t101, t91 * t62 - ((t13 - t116) * t66 + (t100 - t12) * t68) * t110, 0.2e1 * t111 + 0.2e1 * qJD(4) * qJD(5) + t7 + (t110 * t29 + g(3)) * t66 + (t94 - t165) * t68, t2 * qJ(5) - t3 * pkin(4) + t10 * t29 - t12 * t136 - g(3) * t92 + (qJD(5) + t137) * t13 - t94 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t109, -t43, -t167 * t64 - t71, -t10 * t138 - t158 + t159;];
tau_reg = t14;
