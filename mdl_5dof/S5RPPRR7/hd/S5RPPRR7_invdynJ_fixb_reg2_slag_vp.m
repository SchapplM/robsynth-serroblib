% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:58
% EndTime: 2019-12-31 18:00:00
% DurationCPUTime: 1.47s
% Computational Cost: add. (1910->269), mult. (3413->348), div. (0->0), fcn. (1955->10), ass. (0->144)
t69 = qJ(1) + pkin(8);
t62 = sin(t69);
t63 = cos(t69);
t179 = g(1) * t62 - g(2) * t63;
t74 = sin(pkin(8));
t54 = pkin(1) * t74 + qJ(3);
t47 = qJD(1) * t54;
t103 = -t47 * qJD(1) - t179;
t80 = cos(qJ(4));
t148 = qJD(4) * t80;
t135 = t80 * qJDD(1);
t79 = cos(qJ(5));
t140 = t79 * qJD(4);
t76 = sin(qJ(5));
t143 = qJD(5) * t80;
t77 = sin(qJ(4));
t94 = t140 * t77 + t143 * t76;
t17 = qJD(1) * t94 - qJD(5) * t140 - t76 * qJDD(4) - t79 * t135;
t151 = qJD(1) * t80;
t43 = qJD(4) * t76 + t151 * t79;
t157 = t43 * t148 - t17 * t77;
t51 = qJD(1) * t77 + qJD(5);
t75 = cos(pkin(8));
t59 = -pkin(1) * t75 - pkin(2);
t50 = -pkin(6) + t59;
t37 = qJD(1) * t50 + qJD(3);
t24 = qJD(2) * t80 + t37 * t77;
t142 = t24 * qJD(4);
t35 = qJDD(1) * t50 + qJDD(3);
t32 = t80 * t35;
t10 = -t77 * qJDD(2) - t142 + t32;
t8 = -qJDD(4) * pkin(4) - t10;
t90 = g(3) * t77 - t179 * t80 - t8;
t178 = -pkin(7) * qJD(5) * t51 + t90;
t149 = qJD(4) * t77;
t125 = t76 * t149;
t146 = qJD(5) * t43;
t18 = -qJD(1) * t125 - t79 * qJDD(4) + t135 * t76 + t146;
t152 = pkin(1) * qJDD(1);
t109 = pkin(4) * t77 - pkin(7) * t80;
t177 = qJDD(1) * t59;
t134 = qJD(1) * qJD(4);
t117 = t80 * t134;
t136 = t77 * qJDD(1);
t39 = qJDD(5) + t117 + t136;
t160 = t79 * t39;
t176 = t80 * t160 - t51 * t94;
t22 = qJD(4) * pkin(7) + t24;
t34 = t109 + t54;
t25 = t34 * qJD(1);
t5 = -t22 * t76 + t25 * t79;
t6 = t22 * t79 + t25 * t76;
t104 = t5 * t76 - t6 * t79;
t175 = qJD(4) * t104 + t8;
t23 = -qJD(2) * t77 + t37 * t80;
t129 = -t80 * qJDD(2) - t37 * t148 - t77 * t35;
t133 = qJD(2) * qJD(4);
t9 = -t133 * t77 - t129;
t87 = -(t23 * t77 - t24 * t80) * qJD(4) + t10 * t80 + t9 * t77;
t174 = 0.2e1 * qJD(4) * t47 + qJDD(4) * t50;
t173 = -pkin(2) - pkin(6);
t170 = g(3) * t80;
t169 = t5 * t51;
t168 = t6 * t51;
t41 = t151 * t76 - t140;
t167 = t41 * t51;
t166 = t43 * t41;
t165 = t43 * t51;
t164 = t76 * t39;
t163 = t76 * t77;
t162 = t77 * t18;
t161 = t77 * t79;
t159 = t80 * t17;
t158 = t80 * t18;
t60 = t74 * t152;
t155 = qJDD(1) * qJ(3) + t60;
t71 = t77 ^ 2;
t72 = t80 ^ 2;
t154 = t71 - t72;
t82 = qJD(4) ^ 2;
t83 = qJD(1) ^ 2;
t153 = -t82 - t83;
t150 = qJD(4) * t41;
t147 = qJD(5) * t41;
t145 = qJD(5) * t76;
t144 = qJD(5) * t79;
t73 = qJDD(2) - g(3);
t138 = qJDD(4) * t77;
t137 = qJDD(4) * t80;
t132 = qJD(3) * qJD(1);
t130 = t80 * t83 * t77;
t81 = cos(qJ(1));
t127 = t81 * pkin(1) + t63 * pkin(2) + t62 * qJ(3);
t124 = t76 * t148;
t123 = t41 * t149;
t122 = t43 * t149;
t121 = t80 * t140;
t120 = t43 * t143;
t78 = sin(qJ(1));
t118 = -t78 * pkin(1) + t63 * qJ(3);
t115 = -t17 + t147;
t114 = (-t71 - t72) * qJDD(1);
t112 = qJD(5) * t77 + qJD(1);
t111 = t77 * t117;
t110 = pkin(4) * t80 + pkin(7) * t77;
t108 = g(1) * t63 + g(2) * t62;
t106 = g(1) * t78 - g(2) * t81;
t105 = t5 * t79 + t6 * t76;
t100 = g(2) * (t63 * pkin(6) + t127);
t38 = t132 + t155;
t99 = t47 * qJD(3) + t38 * t54;
t20 = t50 * t161 + t34 * t76;
t19 = -t50 * t163 + t34 * t79;
t7 = qJDD(4) * pkin(7) + t9;
t97 = -qJD(5) * t25 + t170 - t7;
t96 = t144 * t51 + t164;
t95 = qJDD(3) + t177;
t93 = qJDD(1) * t54 - t108;
t40 = qJD(4) * t110 + qJD(3);
t21 = -qJD(4) * pkin(4) - t23;
t92 = -pkin(7) * t39 + t21 * t51;
t91 = -t179 * t77 - t170;
t89 = t133 - t103;
t16 = qJD(1) * t40 + qJDD(1) * t109 + t155;
t1 = qJD(5) * t5 + t76 * t16 + t79 * t7;
t12 = t79 * t16;
t2 = -qJD(5) * t6 - t76 * t7 + t12;
t88 = -qJD(5) * t105 + t1 * t79 - t2 * t76;
t86 = -t50 * t82 + t132 + t38 + t93;
t85 = qJD(4) * t21 + t88;
t46 = -t77 * t82 + t137;
t45 = -t80 * t82 - t138;
t44 = t110 * qJD(1);
t36 = t51 * t125;
t29 = t63 * t161 - t62 * t76;
t28 = t63 * t163 + t62 * t79;
t27 = t62 * t161 + t63 * t76;
t26 = -t62 * t163 + t63 * t79;
t14 = t23 * t79 + t44 * t76;
t13 = -t23 * t76 + t44 * t79;
t11 = t79 * t158;
t4 = -qJD(5) * t20 - t124 * t50 + t79 * t40;
t3 = qJD(5) * t19 + t121 * t50 + t76 * t40;
t15 = [0, 0, 0, 0, 0, qJDD(1), t106, g(1) * t81 + g(2) * t78, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t152 * t75 + t179, t108 - 0.2e1 * t60, 0, (t106 + (t74 ^ 2 + t75 ^ 2) * t152) * pkin(1), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - t179 + 0.2e1 * t177, t93 + 0.2e1 * t132 + t155, t95 * t59 - g(1) * (-t62 * pkin(2) + t118) - g(2) * t127 + t99, qJDD(1) * t72 - 0.2e1 * t111, 0.2e1 * t154 * t134 - 0.2e1 * t77 * t135, t46, qJDD(1) * t71 + 0.2e1 * t111, t45, 0, t174 * t80 + t86 * t77, -t174 * t77 + t86 * t80, t114 * t50 + t179 - t87, -g(1) * (t173 * t62 + t118) - t100 + t87 * t50 + t99, -t79 * t159 - t43 * t94, -t11 + (-t120 + t123) * t79 + (t122 + (t17 + t147) * t80) * t76, t157 + t176, t76 * t158 + (t143 * t79 - t125) * t41, -t162 + t36 + (-t96 - t150) * t80, t148 * t51 + t39 * t77, -g(1) * t29 - g(2) * t27 + t19 * t39 + t4 * t51 + (t2 + (-t21 * t76 + t41 * t50) * qJD(4)) * t77 + (qJD(4) * t5 + t144 * t21 - t18 * t50 + t76 * t8) * t80, g(1) * t28 - g(2) * t26 - t20 * t39 - t3 * t51 + (-t1 + (-t21 * t79 + t43 * t50) * qJD(4)) * t77 + (-qJD(4) * t6 - t145 * t21 + t17 * t50 + t79 * t8) * t80, t17 * t19 - t18 * t20 - t3 * t41 - t4 * t43 + t105 * t149 + (qJD(5) * t104 - t1 * t76 - t2 * t79 + t108) * t80, t1 * t20 + t6 * t3 + t2 * t19 + t5 * t4 - g(1) * (t109 * t63 + t118) - t100 + (-g(1) * t173 - g(2) * t109) * t62 + (t149 * t21 - t8 * t80) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, 0, 0, 0, 0, t45, -t46, 0, -t10 * t77 + t9 * t80 - g(3) + (-t23 * t80 - t24 * t77) * qJD(4), 0, 0, 0, 0, 0, 0, t162 + t36 + (-t96 + t150) * t80, t157 - t176, -t11 + (t120 + t123) * t79 + (t115 * t80 - t122) * t76, t175 * t77 + t85 * t80 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t83, t103 + t95, 0, 0, 0, 0, 0, 0, t153 * t77 + t137, t153 * t80 - t138, t114, t103 + t87, 0, 0, 0, 0, 0, 0, -t158 + (t150 - t164) * t77 + (-t112 * t79 - t124) * t51, t159 + (qJD(4) * t43 - t160) * t77 + (t112 * t76 - t121) * t51, (t112 * t43 - t148 * t41 - t162) * t79 + (t112 * t41 + t157) * t76, -t105 * qJD(1) - t175 * t80 + t85 * t77 - t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t154 * t83, t135, -t130, -t136, qJDD(4), t142 + t32 + (-qJD(4) * t37 - t73) * t77 - t89 * t80, t23 * qJD(4) + t77 * t89 + t129 + t170, 0, 0, t79 * t165 - t17 * t76, (-t17 - t167) * t79 + (-t18 - t165) * t76, (t51 * t161 - t43 * t80) * qJD(1) + t96, t76 * t167 - t18 * t79, -t51 * t145 + t160 + (-t51 * t163 + t41 * t80) * qJD(1), -t51 * t151, -pkin(4) * t18 - t13 * t51 - t5 * t151 + t178 * t79 - t24 * t41 + t92 * t76, pkin(4) * t17 + t14 * t51 + t6 * t151 - t178 * t76 - t24 * t43 + t92 * t79, t13 * t43 + t14 * t41 + (t1 - t169 + (-t18 + t146) * pkin(7)) * t79 + (pkin(7) * t115 - t168 - t2) * t76 + t91, -t5 * t13 - t6 * t14 - t21 * t24 + t90 * pkin(4) + (t88 + t91) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, -t41 ^ 2 + t43 ^ 2, t167 - t17, -t166, t165 - t18, t39, -g(1) * t26 - g(2) * t28 - t144 * t22 - t21 * t43 + t76 * t97 + t12 + t168, g(1) * t27 - g(2) * t29 + t21 * t41 + t169 + (qJD(5) * t22 - t16) * t76 + t97 * t79, 0, 0;];
tau_reg = t15;
