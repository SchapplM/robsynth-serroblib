% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR11
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:59
% EndTime: 2019-12-31 18:06:01
% DurationCPUTime: 1.63s
% Computational Cost: add. (1551->275), mult. (2710->346), div. (0->0), fcn. (1451->6), ass. (0->148)
t68 = cos(qJ(4));
t132 = qJD(4) * t68;
t56 = (qJDD(1) * qJ(2));
t57 = (qJD(1) * qJD(2));
t109 = qJDD(3) + t56 + t57;
t29 = -pkin(6) * qJDD(1) + t109;
t39 = qJ(2) * qJD(1) + qJD(3);
t34 = -pkin(6) * qJD(1) + t39;
t65 = sin(qJ(4));
t12 = qJDD(4) * pkin(7) + t34 * t132 + t29 * t65;
t63 = pkin(1) + qJ(3);
t168 = pkin(7) * t68;
t98 = pkin(4) * t65 - t168;
t32 = t98 + t63;
t17 = qJD(1) * t32 - qJD(2);
t155 = t34 * t65;
t18 = qJD(4) * pkin(7) + t155;
t64 = sin(qJ(5));
t67 = cos(qJ(5));
t5 = t17 * t67 - t18 * t64;
t61 = qJDD(1) * pkin(1);
t119 = t61 - qJDD(2);
t55 = qJDD(1) * qJ(3);
t108 = t55 + t119;
t99 = pkin(4) * t68 + pkin(7) * t65;
t25 = t99 * qJD(4) + qJD(3);
t8 = t25 * qJD(1) + t98 * qJDD(1) + t108;
t1 = qJD(5) * t5 + t67 * t12 + t64 * t8;
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t141 = g(1) * t66 - g(2) * t69;
t6 = t17 * t64 + t18 * t67;
t7 = t67 * t8;
t2 = -qJD(5) * t6 - t64 * t12 + t7;
t93 = t5 * t64 - t6 * t67;
t176 = qJD(5) * t93 - t1 * t64 - t2 * t67 - t141;
t118 = qJD(1) * qJD(4);
t105 = t65 * t118;
t120 = t68 * qJDD(1);
t134 = qJD(4) * t64;
t137 = qJD(1) * t68;
t28 = t67 * t137 + t134;
t131 = qJD(5) * t28;
t10 = -t67 * qJDD(4) + (-t105 + t120) * t64 + t131;
t37 = qJD(1) * t65 + qJD(5);
t158 = t28 * t37;
t175 = t158 - t10;
t173 = g(1) * t69 + g(2) * t66;
t170 = qJD(1) * t63;
t35 = -qJD(2) + t170;
t92 = -t35 * qJD(1) - t173;
t125 = t67 * qJD(4);
t26 = t64 * t137 - t125;
t86 = t26 * t37;
t127 = qJD(5) * t68;
t80 = t65 * t125 + t64 * t127;
t9 = qJD(1) * t80 - qJD(5) * t125 - t64 * qJDD(4) - t67 * t120;
t174 = -t9 + t86;
t59 = t65 ^ 2;
t60 = t68 ^ 2;
t139 = t59 + t60;
t172 = t139 * t34;
t133 = qJD(4) * t65;
t11 = -qJDD(4) * pkin(4) + t34 * t133 - t29 * t68;
t167 = g(3) * t65;
t75 = -t173 * t68 - t11 + t167;
t171 = -pkin(7) * qJD(5) * t37 + t75;
t62 = -pkin(6) + qJ(2);
t169 = qJD(4) * (qJD(2) + t35 + t170) + qJDD(4) * t62;
t49 = 2 * t57;
t166 = g(3) * t68;
t165 = t5 * t37;
t164 = t6 * t37;
t163 = t69 * pkin(6);
t162 = t9 * t64;
t161 = t10 * t67;
t160 = t26 * t68;
t159 = t28 * t26;
t157 = t28 * t67;
t156 = t28 * t68;
t154 = t34 * t68;
t71 = qJD(1) ^ 2;
t153 = t59 * t71;
t104 = t68 * t118;
t121 = t65 * qJDD(1);
t24 = qJDD(5) + t104 + t121;
t152 = t64 * t24;
t151 = t64 * t65;
t150 = t64 * t66;
t149 = t65 * t67;
t148 = t65 * t69;
t147 = t66 * t67;
t146 = t67 * t24;
t145 = t67 * t68;
t144 = t67 * t69;
t143 = t68 * t10;
t142 = t69 * pkin(1) + t66 * qJ(2);
t70 = qJD(4) ^ 2;
t138 = -t70 - t71;
t136 = qJD(4) * t26;
t135 = qJD(4) * t28;
t130 = qJD(5) * t64;
t129 = qJD(5) * t65;
t128 = qJD(5) * t67;
t124 = qJDD(1) * t63;
t122 = qJDD(4) * t65;
t117 = qJD(3) * qJD(1);
t116 = t37 * t151;
t115 = t37 * t149;
t114 = t68 * t71 * t65;
t112 = t69 * qJ(3) + t142;
t110 = qJDD(2) - t141;
t106 = t139 * t29;
t103 = (2 * t56) + t49 - t173;
t102 = qJD(1) + t129;
t101 = -t61 + t110;
t100 = t65 * t104;
t95 = -t62 * t129 + t25;
t94 = t5 * t67 + t6 * t64;
t47 = t69 * qJ(2);
t91 = -t63 * t66 + t47;
t89 = -t55 + t101;
t30 = t108 + t117;
t88 = t35 * qJD(3) + t30 * t63;
t84 = -qJD(5) * t17 - t12 + t166;
t83 = t37 * t128 + t152;
t82 = t37 * t130 - t146;
t79 = -t29 - t92;
t19 = -qJD(4) * pkin(4) - t154;
t78 = -pkin(7) * t24 + t37 * t19;
t77 = -t173 * t65 - t166;
t76 = qJD(2) * t65 + qJD(5) * t32 + t62 * t132;
t74 = -t94 * qJD(5) + t1 * t67 - t2 * t64;
t72 = -t62 * t70 + t117 + t124 + t141 + t30;
t44 = t60 * t71;
t42 = qJDD(4) * t68;
t31 = t99 * qJD(1);
t23 = t65 * t144 - t150;
t22 = -t64 * t148 - t147;
t21 = -t65 * t147 - t64 * t69;
t20 = t65 * t150 - t144;
t16 = t62 * t149 + t32 * t64;
t15 = -t62 * t151 + t32 * t67;
t14 = t34 * t145 + t31 * t64;
t13 = -t64 * t154 + t31 * t67;
t4 = -t76 * t64 + t95 * t67;
t3 = t95 * t64 + t76 * t67;
t27 = [0, 0, 0, 0, 0, qJDD(1), t141, t173, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t61 + t110, t103, t119 * pkin(1) - g(1) * (-pkin(1) * t66 + t47) - g(2) * t142 + (t56 + t49) * qJ(2), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) + t103, 0.2e1 * t117 - t89 + t124, -g(1) * t91 - g(2) * t112 + t109 * qJ(2) + t39 * qJD(2) + t88, qJDD(1) * t60 - 0.2e1 * t100, -0.2e1 * t65 * t120 + 0.2e1 * (t59 - t60) * t118, -t65 * t70 + t42, qJDD(1) * t59 + 0.2e1 * t100, -t68 * t70 - t122, 0, t169 * t68 + t72 * t65, -t169 * t65 + t72 * t68, t173 + t139 * (-qJDD(1) * t62 - t29 - t57), -g(1) * (t91 - t163) - g(2) * (-t66 * pkin(6) + t112) + t62 * t106 + qJD(2) * t172 + t88, -t145 * t9 - t28 * t80, (t26 * t67 + t28 * t64) * t133 + (-t161 + t162 + (t26 * t64 - t157) * qJD(5)) * t68, (-t125 * t37 - t9) * t65 + (-t82 + t135) * t68, t64 * t143 + (t67 * t127 - t133 * t64) * t26, (t134 * t37 - t10) * t65 + (-t83 - t136) * t68, t132 * t37 + t24 * t65, -g(1) * t21 - g(2) * t23 + t15 * t24 + t37 * t4 + (t2 + (-t19 * t64 + t26 * t62) * qJD(4)) * t65 + (-qJD(2) * t26 + qJD(4) * t5 - t10 * t62 + t11 * t64 + t128 * t19) * t68, -g(1) * t20 - g(2) * t22 - t16 * t24 - t3 * t37 + (-t1 + (-t19 * t67 + t28 * t62) * qJD(4)) * t65 + (-qJD(2) * t28 - qJD(4) * t6 + t11 * t67 - t130 * t19 + t62 * t9) * t68, -t10 * t16 + t94 * t133 + t15 * t9 + t176 * t68 - t26 * t3 - t28 * t4, t1 * t16 + t6 * t3 + t2 * t15 + t5 * t4 - t11 * t68 * t62 - g(1) * (t47 - t163) - g(2) * (pkin(4) * t148 - t69 * t168 + t112) + (-qJD(2) * t68 + t133 * t62) * t19 + (g(2) * pkin(6) + g(1) * t32) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t71, -qJ(2) * t71 + t101, 0, 0, 0, 0, 0, 0, 0, -t71, -qJDD(1), (-qJD(3) - t39) * qJD(1) + t89, 0, 0, 0, 0, 0, 0, -0.2e1 * t104 - t121, 0.2e1 * t105 - t120, t44 + t153, (-qJD(3) - t172) * qJD(1) + t89, 0, 0, 0, 0, 0, 0, (t116 + t160) * qJD(1) + t82, (t115 + t156) * qJD(1) + t83, t174 * t67 - t175 * t64, (t19 * t68 + t65 * t93) * qJD(1) + t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t71, t92 + t109, 0, 0, 0, 0, 0, 0, t138 * t65 + t42, t138 * t68 - t122, -t139 * qJDD(1), t106 + t92, 0, 0, 0, 0, 0, 0, -t143 + (t136 - t152) * t65 + (-t102 * t67 - t132 * t64) * t37, t68 * t9 + (t135 - t146) * t65 + (t102 * t64 - t125 * t68) * t37, (-t10 * t65 + t102 * t28 - t132 * t26) * t67 + (t102 * t26 + t132 * t28 - t9 * t65) * t64, -t94 * qJD(1) + (-qJD(4) * t93 - t11) * t68 + (qJD(4) * t19 + t74) * t65 - t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t44 - t153, t120, -t114, -t121, qJDD(4), -t68 * t79 + t167, t65 * t79 + t166, 0, 0, t37 * t157 - t162, (-t9 - t86) * t67 + (-t10 - t158) * t64, (t115 - t156) * qJD(1) + t83, t64 * t86 - t161, (-t116 + t160) * qJD(1) - t82, -t37 * t137, -pkin(4) * t10 - t13 * t37 - t5 * t137 - t26 * t155 + t171 * t67 + t78 * t64, pkin(4) * t9 + t6 * t137 + t14 * t37 - t28 * t155 - t171 * t64 + t78 * t67, t13 * t28 + t14 * t26 + (t1 - t165 + (-t10 + t131) * pkin(7)) * t67 + (-t2 - t164 + (qJD(5) * t26 - t9) * pkin(7)) * t64 + t77, -t19 * t155 - t5 * t13 - t6 * t14 + t75 * pkin(4) + (t74 + t77) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, -t26 ^ 2 + t28 ^ 2, t174, -t159, t175, t24, -g(1) * t22 + g(2) * t20 - t128 * t18 - t19 * t28 + t64 * t84 + t164 + t7, g(1) * t23 - g(2) * t21 + t19 * t26 + t165 + (qJD(5) * t18 - t8) * t64 + t84 * t67, 0, 0;];
tau_reg = t27;
