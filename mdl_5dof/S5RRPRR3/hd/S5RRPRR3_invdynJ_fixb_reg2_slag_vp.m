% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:40
% EndTime: 2020-01-03 12:00:41
% DurationCPUTime: 1.08s
% Computational Cost: add. (2923->215), mult. (4964->264), div. (0->0), fcn. (2944->16), ass. (0->150)
t100 = qJ(1) + qJ(2);
t90 = pkin(9) + t100;
t83 = qJ(4) + t90;
t74 = sin(t83);
t75 = cos(t83);
t137 = -g(2) * t75 - g(3) * t74;
t95 = qJDD(1) + qJDD(2);
t87 = qJDD(4) + t95;
t180 = t87 * pkin(4);
t104 = sin(qJ(4));
t108 = cos(qJ(4));
t101 = sin(pkin(9));
t102 = cos(pkin(9));
t105 = sin(qJ(2));
t178 = pkin(1) * t105;
t149 = qJD(1) * t178;
t109 = cos(qJ(2));
t152 = qJD(1) * t109;
t96 = qJD(1) + qJD(2);
t61 = pkin(1) * t152 + pkin(2) * t96;
t35 = -t101 * t149 + t102 * t61;
t33 = pkin(3) * t96 + t35;
t36 = t101 * t61 + t102 * t149;
t18 = t104 * t33 + t108 * t36;
t150 = qJDD(1) * t105;
t122 = pkin(1) * (qJD(2) * t152 + t150);
t184 = pkin(2) * t95;
t177 = pkin(1) * t109;
t85 = qJDD(1) * t177;
t46 = -qJD(2) * t149 + t184 + t85;
t29 = -t101 * t122 + t102 * t46;
t21 = t95 * pkin(3) + t29;
t161 = t101 * t46;
t186 = t102 * t122;
t30 = t161 + t186;
t8 = -t18 * qJD(4) - t104 * t30 + t108 * t21;
t6 = -t180 - t8;
t131 = t137 - t6;
t154 = t102 * t105;
t128 = pkin(1) * (-t101 * t109 - t154);
t50 = qJD(1) * t128;
t155 = t101 * t105;
t127 = pkin(1) * (t102 * t109 - t155);
t52 = qJD(1) * t127;
t176 = pkin(2) * t101;
t175 = pkin(2) * t102;
t79 = pkin(3) + t175;
t55 = -t104 * t176 + t108 * t79;
t171 = t55 * qJD(4) - t104 * t50 - t108 * t52;
t88 = qJD(4) + t96;
t187 = t171 * t88;
t111 = qJD(5) ^ 2;
t57 = t104 * t79 + t108 * t176;
t170 = t57 * qJD(4) - t104 * t52 + t108 * t50;
t146 = t170 * t88;
t48 = -pkin(4) - t55;
t49 = pkin(8) + t57;
t185 = t111 * t49 + t48 * t87 + t146;
t17 = -t104 * t36 + t108 * t33;
t103 = sin(qJ(5));
t107 = cos(qJ(5));
t16 = pkin(8) * t88 + t18;
t11 = qJD(3) * t107 - t103 * t16;
t159 = t107 * t16;
t12 = qJD(3) * t103 + t159;
t157 = qJD(5) * t11;
t142 = -t17 * qJD(4) - t104 * t21 - t108 * t30;
t5 = pkin(8) * t87 - t142;
t2 = t103 * qJDD(3) + t107 * t5 + t157;
t156 = qJD(5) * t12;
t89 = t107 * qJDD(3);
t3 = -t103 * t5 - t156 + t89;
t117 = -t3 * t103 + t2 * t107 + (-t103 * t12 - t107 * t11) * qJD(5);
t69 = g(3) * t75;
t168 = -g(2) * t74 + t69;
t183 = pkin(4) * t88;
t84 = pkin(2) + t177;
t54 = -pkin(1) * t155 + t102 * t84;
t47 = pkin(3) + t54;
t56 = pkin(1) * t154 + t101 * t84;
t27 = -t104 * t56 + t108 * t47;
t51 = qJD(2) * t128;
t53 = qJD(2) * t127;
t9 = t27 * qJD(4) + t104 * t51 + t108 * t53;
t179 = t88 * t9;
t28 = t104 * t47 + t108 * t56;
t10 = t28 * qJD(4) + t104 * t53 - t108 * t51;
t174 = t10 * t88;
t173 = t17 * t88;
t172 = t18 * t88;
t169 = t75 * pkin(4) + t74 * pkin(8);
t77 = sin(t90);
t72 = pkin(3) * t77;
t91 = sin(t100);
t81 = pkin(2) * t91;
t167 = t72 + t81;
t78 = cos(t90);
t73 = pkin(3) * t78;
t92 = cos(t100);
t82 = pkin(2) * t92;
t166 = t73 + t82;
t106 = sin(qJ(1));
t93 = t106 * pkin(1);
t165 = t81 + t93;
t110 = cos(qJ(1));
t94 = t110 * pkin(1);
t164 = t82 + t94;
t97 = t103 ^ 2;
t98 = t107 ^ 2;
t163 = t97 - t98;
t162 = t97 + t98;
t153 = t103 * t107;
t99 = qJDD(3) - g(1);
t151 = qJD(5) * t107;
t86 = t88 ^ 2;
t148 = t86 * t153;
t147 = g(2) * t91 - g(3) * t92;
t145 = t162 * t87;
t15 = -t17 - t183;
t144 = -t131 * t103 + t15 * t151;
t141 = t166 + t169;
t140 = qJD(1) * (-qJD(2) + t96);
t139 = qJD(2) * (-qJD(1) - t96);
t138 = t103 * t88 * t151;
t136 = -g(2) * t92 - g(3) * t91;
t135 = -g(2) * t110 - g(3) * t106;
t67 = t74 * pkin(4);
t133 = -pkin(8) * t75 + t167 + t67;
t132 = t103 * t11 - t107 * t12;
t130 = t142 - t168;
t129 = t136 + t85;
t125 = pkin(8) * t111 - t172 - t180;
t23 = -pkin(4) - t27;
t24 = pkin(8) + t28;
t124 = t111 * t24 + t23 * t87 + t174;
t123 = -pkin(8) * qJDD(5) + (t17 - t183) * qJD(5);
t121 = -qJDD(5) * t24 + (t23 * t88 - t9) * qJD(5);
t120 = -qJDD(5) * t49 + (t48 * t88 - t171) * qJD(5);
t118 = -qJD(3) * qJD(5) - t15 * t88 - t168 - t5;
t116 = t168 + t117;
t115 = g(2) * t77 - g(3) * t78 - t186;
t114 = -g(2) * t78 - g(3) * t77 + t29;
t113 = t137 + t8;
t63 = qJDD(5) * t107 - t103 * t111;
t62 = qJDD(5) * t103 + t107 * t111;
t39 = t87 * t98 - 0.2e1 * t138;
t38 = t87 * t97 + 0.2e1 * t138;
t32 = -0.2e1 * t163 * t88 * qJD(5) + 0.2e1 * t87 * t153;
t13 = t15 * qJD(5) * t103;
t1 = [0, 0, 0, 0, 0, qJDD(1), t135, g(2) * t106 - g(3) * t110, 0, 0, 0, 0, 0, 0, 0, t95, (t105 * t139 + t109 * t95) * pkin(1) + t129, ((-qJDD(1) - t95) * t105 + t109 * t139) * pkin(1) + t147, 0, (t135 + (t105 ^ 2 + t109 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t95, t51 * t96 + t54 * t95 + t114, -t53 * t96 - t56 * t95 + t115 - t161, 0, -g(2) * t164 - g(3) * t165 + t29 * t54 + t30 * t56 + t35 * t51 + t36 * t53, 0, 0, 0, 0, 0, t87, t27 * t87 + t113 - t174, -t28 * t87 + t130 - t179, 0, -t142 * t28 + t18 * t9 + t8 * t27 - t17 * t10 - g(2) * (t73 + t164) - g(3) * (t72 + t165), t38, t32, t62, t39, t63, 0, t13 + t121 * t103 + (-t124 + t131) * t107, t103 * t124 + t107 * t121 + t144, t145 * t24 + t162 * t179 + t116, t6 * t23 + t15 * t10 - g(2) * (t94 + t141) - g(3) * (t133 + t93) - t132 * t9 + t117 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t140 * t178 + t129, (t109 * t140 - t150) * pkin(1) + t147, 0, 0, 0, 0, 0, 0, 0, t95, t95 * t175 - t50 * t96 + t114, t52 * t96 + (-t46 - t184) * t101 + t115, 0, -t35 * t50 - t36 * t52 + (t101 * t30 + t102 * t29 + t136) * pkin(2), 0, 0, 0, 0, 0, t87, t55 * t87 + t113 - t146, -t57 * t87 + t130 - t187, 0, -g(2) * t166 - g(3) * t167 - t142 * t57 - t17 * t170 + t171 * t18 + t8 * t55, t38, t32, t62, t39, t63, 0, t13 + t120 * t103 + (t131 - t185) * t107, t185 * t103 + t120 * t107 + t144, t49 * t145 + t162 * t187 + t116, t6 * t48 - g(2) * t141 - g(3) * t133 + t170 * t15 + ((t2 - t157) * t49 + t171 * t12) * t107 + ((-t3 - t156) * t49 - t171 * t11) * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, 0, 0, 0, 0, 0, t63, -t62, 0, -qJD(5) * t132 + t2 * t103 + t3 * t107 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t113 + t172, t130 + t173, 0, 0, t38, t32, t62, t39, t63, 0, t13 + t123 * t103 + (-t125 + t131) * t107, t103 * t125 + t107 * t123 + t144, pkin(8) * t145 - t162 * t173 + t116, -t6 * pkin(4) - t15 * t18 - g(2) * t169 - g(3) * t67 + t132 * t17 + (t117 + t69) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t163 * t86, t103 * t87, t148, t107 * t87, qJDD(5), -g(1) * t107 + t89 + (t12 - t159) * qJD(5) + t118 * t103, t157 + (qJD(5) * t16 - t99) * t103 + t118 * t107, 0, 0;];
tau_reg = t1;
