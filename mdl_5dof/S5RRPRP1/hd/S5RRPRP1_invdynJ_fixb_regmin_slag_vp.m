% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:59
% EndTime: 2022-01-20 10:20:04
% DurationCPUTime: 1.11s
% Computational Cost: add. (1382->237), mult. (2201->277), div. (0->0), fcn. (1262->12), ass. (0->158)
t185 = qJDD(3) - g(3);
t103 = cos(qJ(2));
t159 = pkin(1) * qJD(2);
t135 = qJD(1) * t159;
t100 = sin(qJ(2));
t144 = qJDD(1) * t100;
t184 = pkin(1) * t144 + t103 * t135;
t95 = qJ(1) + qJ(2);
t85 = sin(t95);
t86 = cos(t95);
t183 = g(1) * t85 - g(2) * t86;
t102 = cos(qJ(4));
t92 = qJD(1) + qJD(2);
t158 = qJ(5) * t92;
t160 = pkin(1) * qJD(1);
t138 = t103 * t160;
t49 = t92 * pkin(2) + t138;
t139 = t100 * t160;
t97 = cos(pkin(8));
t66 = t97 * t139;
t96 = sin(pkin(8));
t29 = t96 * t49 + t66;
t24 = t92 * pkin(7) + t29;
t131 = t24 + t158;
t117 = t131 * t102;
t182 = pkin(2) * t85;
t99 = sin(qJ(4));
t93 = t99 ^ 2;
t181 = pkin(4) * t93;
t82 = pkin(8) + t95;
t70 = sin(t82);
t180 = g(1) * t70;
t71 = cos(t82);
t179 = g(1) * t71;
t178 = g(2) * t70;
t177 = g(2) * t71;
t91 = qJDD(1) + qJDD(2);
t175 = t91 * pkin(3);
t174 = t97 * pkin(2);
t84 = t102 * qJD(3);
t12 = -t131 * t99 + t84;
t157 = qJD(4) * pkin(4);
t9 = t12 + t157;
t173 = -t12 + t9;
t172 = g(3) * t102;
t101 = sin(qJ(1));
t171 = t101 * pkin(1);
t170 = t102 * pkin(4);
t169 = t103 * pkin(1);
t65 = t96 * t139;
t28 = t97 * t49 - t65;
t23 = -t92 * pkin(3) - t28;
t168 = t23 * t92;
t78 = pkin(3) + t170;
t50 = -t78 - t174;
t167 = t50 * t91;
t166 = t71 * t99;
t165 = t92 * t99;
t69 = t99 * t91;
t155 = t100 * t97;
t79 = pkin(2) + t169;
t164 = pkin(1) * t155 + t96 * t79;
t163 = g(1) * t86 + g(2) * t85;
t94 = t102 ^ 2;
t162 = -t93 - t94;
t161 = t93 - t94;
t156 = t100 * t96;
t154 = t102 * t70;
t153 = t102 * t91;
t152 = t102 * t92;
t151 = t102 * t99;
t40 = pkin(7) + t164;
t150 = -qJ(5) - t40;
t72 = t96 * pkin(2) + pkin(7);
t149 = -qJ(5) - t72;
t148 = qJD(4) * t92;
t147 = qJDD(4) * pkin(4);
t146 = t99 * qJD(4);
t145 = t102 * qJD(4);
t18 = -t78 * t92 + qJD(5) - t28;
t80 = qJDD(1) * t169;
t37 = t91 * pkin(2) - t100 * t135 + t80;
t16 = -t184 * t96 + t97 * t37;
t137 = t92 * t146;
t62 = pkin(4) * t137;
t111 = qJDD(5) - t16 + t62;
t5 = -t78 * t91 + t111;
t55 = g(2) * t166;
t143 = t18 * t145 + t5 * t99 + t55;
t10 = -t16 - t175;
t142 = t10 * t99 + t23 * t145 + t55;
t41 = t96 * t138 + t66;
t43 = t97 * t138 - t65;
t56 = g(1) * t154;
t141 = t43 * t146 + t41 * t152 + t56;
t17 = t184 * t97 + t96 * t37;
t140 = pkin(4) * t146;
t136 = -t5 - t177;
t133 = -t10 - t177;
t132 = -pkin(1) * t156 + t97 * t79;
t130 = qJD(4) * t150;
t129 = qJD(4) * t149;
t128 = 0.2e1 * t92 * t145;
t127 = qJD(1) * (-qJD(2) + t92);
t126 = qJD(2) * (-qJD(1) - t92);
t11 = t91 * pkin(7) + t17;
t125 = -qJD(4) * qJD(3) - t11;
t39 = -pkin(3) - t132;
t124 = t80 + t183;
t77 = pkin(2) * t86;
t98 = -qJ(5) - pkin(7);
t122 = -t70 * t98 + t71 * t78 + t77;
t121 = -t178 - t179;
t120 = -t41 * t92 - t180;
t42 = (t103 * t96 + t155) * t159;
t30 = t42 + t140;
t34 = t39 - t170;
t119 = t30 * t92 + t34 * t91;
t105 = qJD(4) ^ 2;
t73 = -pkin(3) - t174;
t118 = -t105 * t72 - t73 * t91;
t81 = t102 * qJDD(3);
t116 = g(1) * t166 + t99 * t178 - t172 + t81;
t113 = -qJ(5) * t91 + t125;
t107 = qJD(5) * t92 - t113;
t21 = t24 * t146;
t3 = -t21 + (-qJ(5) * t148 + qJDD(3)) * t99 + t107 * t102;
t115 = t3 * t102 + t121;
t114 = -qJDD(4) * t72 + t73 * t148;
t112 = g(2) * t154 + t102 * t179 - t185 * t99 + t21;
t110 = -t70 * t78 - t71 * t98 - t182;
t109 = t105 * t40 + t39 * t91 + t42 * t92;
t44 = (t103 * t97 - t156) * t159;
t108 = -qJDD(4) * t40 + (t39 * t92 - t44) * qJD(4);
t106 = (-qJD(5) - t18) * t92 + t113;
t104 = cos(qJ(1));
t90 = t92 ^ 2;
t88 = t104 * pkin(1);
t87 = t102 * qJ(5);
t83 = t102 * qJD(5);
t52 = qJDD(4) * t102 - t105 * t99;
t51 = qJDD(4) * t99 + t105 * t102;
t47 = t102 * t72 + t87;
t46 = t149 * t99;
t38 = t99 * t128 + t93 * t91;
t36 = -t99 * qJD(5) + t102 * t129;
t35 = t99 * t129 + t83;
t32 = t43 * t145;
t26 = t102 * t40 + t87;
t25 = t150 * t99;
t22 = -0.2e1 * t161 * t148 + 0.2e1 * t91 * t151;
t19 = t23 * t146;
t14 = t18 * t146;
t13 = t99 * qJD(3) + t117;
t7 = (-qJD(5) - t44) * t99 + t102 * t130;
t6 = t102 * t44 + t99 * t130 + t83;
t2 = -qJD(4) * t117 - t107 * t99 + t147 + t81;
t1 = [qJDD(1), g(1) * t101 - g(2) * t104, g(1) * t104 + g(2) * t101, t91, (t100 * t126 + t103 * t91) * pkin(1) + t124, ((-qJDD(1) - t91) * t100 + t103 * t126) * pkin(1) + t163, t17 * t164 + t29 * t44 + t16 * t132 - t28 * t42 - g(1) * (-t171 - t182) - g(2) * (t77 + t88), t38, t22, t51, t52, 0, t19 + t56 + t108 * t99 + (-t109 + t133) * t102, t108 * t102 + (t109 - t180) * t99 + t142, t25 * qJDD(4) + t14 + t56 + (t34 * t165 + t7) * qJD(4) + (-t119 + t136) * t102, -t26 * qJDD(4) + (t34 * t152 - t6) * qJD(4) + (t119 - t180) * t99 + t143, (t26 * t91 + t6 * t92 + (-t25 * t92 - t9) * qJD(4)) * t102 + (-t25 * t91 - t7 * t92 - t2 + (-t26 * t92 - t13) * qJD(4)) * t99 + t115, t3 * t26 + t13 * t6 + t2 * t25 + t9 * t7 + t5 * t34 + t18 * t30 - g(1) * (t110 - t171) - g(2) * (t122 + t88); 0, 0, 0, t91, t100 * pkin(1) * t127 + t124, (t103 * t127 - t144) * pkin(1) + t163, t28 * t41 - t29 * t43 + (t16 * t97 + t17 * t96 + t183) * pkin(2), t38, t22, t51, t52, 0, t19 + t114 * t99 + (t118 + t133) * t102 + t141, t32 + t114 * t102 + (-t118 + t120) * t99 + t142, t46 * qJDD(4) + t14 + (t50 * t165 + t36) * qJD(4) + (t136 - t62 - t167) * t102 + t141, -t47 * qJDD(4) + t32 + (t120 + t167) * t99 + (-t35 + (t102 * t50 + t181) * t92) * qJD(4) + t143, (-qJD(4) * t9 + t47 * t91) * t102 + (-t13 * qJD(4) - t46 * t91 - t2) * t99 + (t102 * t35 - t36 * t99 + t162 * t43 + (-t102 * t46 - t47 * t99) * qJD(4)) * t92 + t115, t3 * t47 + t2 * t46 + t5 * t50 - g(1) * t110 - g(2) * t122 + (t99 * t43 + t36) * t9 + (-t41 + t140) * t18 + (-t102 * t43 + t35) * t13; 0, 0, 0, 0, 0, 0, t185, 0, 0, 0, 0, 0, t52, -t51, t52, -t51, 0, t2 * t102 + t3 * t99 - g(3) + (t13 * t102 - t9 * t99) * qJD(4); 0, 0, 0, 0, 0, 0, 0, -t90 * t151, t161 * t90, t69, t153, qJDD(4), (-t11 - t168) * t99 + t116, (-t99 * t24 + t84) * qJD(4) + (t125 - t168) * t102 + t112, 0.2e1 * t147 + (t13 - t117) * qJD(4) + (t90 * t170 + t106) * t99 + t116, -t90 * t181 + (t99 * t158 + t12) * qJD(4) + t106 * t102 + t112, -pkin(4) * t69 + (-t157 + t173) * t152, t173 * t13 + (-t172 + t2 + (-t18 * t92 - t121) * t99) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t137 - t153, t69 + t128, t162 * t90, t9 * t165 - t175 - t180 + t177 + (-pkin(4) * t91 - t13 * t92) * t102 + t111;];
tau_reg = t1;
