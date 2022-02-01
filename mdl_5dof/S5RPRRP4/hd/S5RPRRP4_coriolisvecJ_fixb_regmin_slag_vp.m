% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:32:52
% EndTime: 2022-01-23 09:32:55
% DurationCPUTime: 1.17s
% Computational Cost: add. (2029->201), mult. (5435->285), div. (0->0), fcn. (3782->6), ass. (0->144)
t141 = qJD(3) + qJD(4);
t95 = sin(pkin(8));
t182 = t141 * t95;
t91 = t95 ^ 2;
t181 = 0.2e1 * t91;
t100 = cos(qJ(3));
t97 = sin(qJ(4));
t98 = sin(qJ(3));
t99 = cos(qJ(4));
t70 = t97 * t100 + t99 * t98;
t60 = t70 * t95;
t96 = cos(pkin(8));
t92 = t96 ^ 2;
t180 = t181 + t92;
t179 = t100 * t98;
t158 = qJ(2) * t98;
t174 = pkin(7) * t95;
t109 = -t100 * t174 - t96 * t158;
t73 = -t96 * pkin(2) - t95 * pkin(6) - pkin(1);
t63 = t73 * qJD(1) + qJD(2);
t59 = t100 * t63;
t42 = t109 * qJD(1) + t59;
t147 = t96 * qJD(1);
t84 = -qJD(3) + t147;
t28 = -t84 * pkin(3) + t42;
t153 = qJ(2) * t100;
t138 = t96 * t153;
t152 = qJD(1) * t95;
t139 = t98 * t152;
t167 = t98 * t63;
t43 = -pkin(7) * t139 + qJD(1) * t138 + t167;
t35 = t97 * t43;
t130 = t99 * t28 - t35;
t135 = t100 * t152;
t121 = t99 * t135;
t122 = t97 * t139;
t56 = t121 - t122;
t50 = t56 * qJ(5);
t6 = t130 - t50;
t134 = -t73 + t174;
t178 = t134 * t98 - t138;
t177 = t56 ^ 2;
t79 = -qJD(4) + t84;
t5 = -t79 * pkin(4) + t6;
t176 = t5 - t6;
t175 = pkin(3) * t79;
t110 = qJD(1) * t70;
t53 = t95 * t110;
t173 = t53 * t79;
t172 = t56 * t53;
t171 = t56 * t79;
t64 = pkin(3) * t139 + qJ(2) * t152;
t170 = t64 * t56;
t169 = t95 * t98;
t168 = t97 * t98;
t37 = t99 * t43;
t166 = t99 * t42 - t35;
t154 = t99 * t100;
t69 = t154 - t168;
t165 = (t141 - t147) * t69;
t48 = t141 * t70;
t164 = t96 * t110 - t48;
t142 = qJD(1) * qJD(2);
t133 = t96 * t142;
t143 = qJD(3) * t100;
t163 = t100 * t133 + t63 * t143;
t144 = qJD(2) * t100;
t162 = t73 * t143 + t96 * t144;
t161 = t141 * t122;
t123 = pkin(3) * t135;
t62 = qJD(3) * t123 + t95 * t142;
t66 = (pkin(3) * t143 + qJD(2)) * t95;
t68 = pkin(3) * t169 + t95 * qJ(2);
t160 = t91 + t92;
t159 = -t100 ^ 2 + t98 ^ 2;
t33 = t141 * t60;
t26 = qJD(1) * t33;
t157 = t26 * qJ(5);
t156 = t53 * qJ(5);
t101 = qJD(1) ^ 2;
t155 = t91 * t101;
t151 = qJD(2) * t98;
t150 = qJD(3) * t98;
t149 = qJD(4) * t97;
t148 = qJD(4) * t99;
t146 = qJD(3) + t84;
t126 = -t53 * pkin(4) - qJD(5);
t40 = -t126 + t64;
t145 = qJD(5) + t40;
t137 = t96 * t151;
t136 = qJ(2) * t150;
t132 = qJD(1) * qJD(3) * t91;
t108 = t109 * qJD(3);
t23 = qJD(1) * t108 + t163;
t24 = -t63 * t150 + (-t137 + (pkin(7) * t169 - t138) * qJD(3)) * qJD(1);
t131 = -t97 * t23 + t99 * t24;
t38 = t108 + t162;
t39 = t178 * qJD(3) - t137;
t129 = -t97 * t38 + t99 * t39;
t128 = -t97 * t42 - t37;
t127 = t160 * t101;
t125 = -t28 * t148 + t43 * t149 - t99 * t23 - t97 * t24;
t124 = qJD(1) * t146;
t116 = t154 * t182;
t27 = qJD(1) * t116 - t161;
t16 = t27 * pkin(4) + t62;
t120 = t95 * t124;
t119 = -t97 * t28 - t37;
t44 = (-pkin(3) - t158) * t96 - t134 * t100;
t118 = t178 * t99 - t97 * t44;
t117 = 0.2e1 * t160 * t142;
t115 = t64 * t53 + t125;
t114 = t27 * qJ(5) + t125;
t113 = qJD(3) * t95 * (t84 + t147);
t112 = t44 * t148 + t149 * t178 + t99 * t38 + t97 * t39;
t107 = -t84 ^ 2 - t155;
t106 = t145 * t53 + t114;
t105 = t119 * qJD(4) + t131;
t104 = t105 + t157;
t103 = (-t37 + (-t28 + t175) * t97) * qJD(4) + t131;
t102 = t48 * t152;
t88 = t99 * pkin(3) + pkin(4);
t67 = t148 * t175;
t61 = t69 * t95;
t51 = t53 ^ 2;
t46 = t56 * pkin(4) + t123;
t45 = t60 * pkin(4) + t68;
t34 = -t168 * t182 + t116;
t18 = t34 * pkin(4) + t66;
t17 = -t51 + t177;
t15 = -t141 * t121 + t161 - t171;
t14 = -t102 - t173;
t13 = -t60 * qJ(5) - t118;
t12 = -t96 * pkin(4) - t61 * qJ(5) + t178 * t97 + t99 * t44;
t11 = -t53 * t152 - t164 * t79;
t10 = -t56 * t152 + t165 * t79;
t9 = -t50 + t166;
t8 = t128 + t156;
t7 = -t119 - t156;
t4 = t33 * qJ(5) + t118 * qJD(4) - t61 * qJD(5) + t129;
t3 = -t34 * qJ(5) - t60 * qJD(5) + t112;
t2 = -t56 * qJD(5) + t104;
t1 = -t53 * qJD(5) - t114;
t19 = [0, 0, 0, 0, t117, qJ(2) * t117, -0.2e1 * t132 * t179, 0.2e1 * t159 * t132, t98 * t113, t100 * t113, 0, t84 * t137 + (-(-t98 * t73 - t138) * t84 + t96 * t167) * qJD(3) + t180 * qJD(1) * (qJ(2) * t143 + t151), (-t96 * t136 + t162) * t84 + t163 * t96 + (-t180 * t136 + t144 * t181) * qJD(1), -t26 * t61 - t56 * t33, t26 * t60 - t61 * t27 + t33 * t53 - t56 * t34, t26 * t96 + t33 * t79, t27 * t96 + t34 * t79, 0, -t129 * t79 - t131 * t96 + t66 * t53 + t68 * t27 + t62 * t60 + t64 * t34 + (-t118 * t79 - t119 * t96) * qJD(4), t112 * t79 - t125 * t96 - t68 * t26 - t64 * t33 + t66 * t56 + t62 * t61, t16 * t60 + t18 * t53 - t2 * t96 + t45 * t27 + t40 * t34 - t4 * t79, t1 * t96 + t16 * t61 + t18 * t56 - t45 * t26 + t3 * t79 - t40 * t33, -t1 * t60 + t12 * t26 - t13 * t27 - t2 * t61 - t3 * t53 + t5 * t33 - t7 * t34 - t4 * t56, t1 * t13 + t2 * t12 + t16 * t45 + t40 * t18 + t7 * t3 + t5 * t4; 0, 0, 0, 0, -t127, -qJ(2) * t127, 0, 0, 0, 0, 0, t107 * t98, t107 * t100, 0, 0, 0, 0, 0, t11, t10, t11, t10, -t164 * t56 - t165 * t53 + t69 * t26 - t70 * t27, t1 * t70 - t40 * t152 + t164 * t5 + t165 * t7 + t2 * t69; 0, 0, 0, 0, 0, 0, t155 * t179, -t159 * t155, -t98 * t120, -t100 * t120, 0, (-t146 * t63 - t133) * t98 + (-t96 * t124 - t155) * t153, -t59 * t84 + (t146 * t147 + t155) * t158 - t163, t172, t17, t14, t15, 0, -t53 * t123 + t128 * t79 + t103 - t170, -t56 * t123 - t166 * t79 + t115 + t67, -t145 * t56 - t46 * t53 + t8 * t79 + t103 + t157, -t46 * t56 - t9 * t79 + t106 + t67, t88 * t26 + (t7 + t8) * t56 + (-t5 + t9) * t53 + (-t27 * t97 + (-t53 * t99 + t56 * t97) * qJD(4)) * pkin(3), t2 * t88 - t40 * t46 - t5 * t8 - t7 * t9 + (t1 * t97 + (-t5 * t97 + t7 * t99) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t17, t14, t15, 0, t119 * t79 + t105 - t170, -t130 * t79 + t115, -t7 * t79 + (t126 - t40) * t56 + t104, -t177 * pkin(4) - t6 * t79 + t106, t26 * pkin(4) - t176 * t53, t176 * t7 + (-t40 * t56 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 - t171, -t102 + t173, -t51 - t177, t5 * t56 + t7 * t53 + t16;];
tauc_reg = t19;
