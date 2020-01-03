% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR13_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:27
% EndTime: 2019-12-31 19:15:32
% DurationCPUTime: 1.66s
% Computational Cost: add. (1562->261), mult. (3536->400), div. (0->0), fcn. (2338->6), ass. (0->144)
t98 = cos(qJ(4));
t145 = t98 * qJD(3);
t99 = cos(qJ(3));
t157 = qJD(1) * t99;
t95 = sin(qJ(4));
t64 = t95 * t157 - t145;
t136 = t98 * t157;
t155 = qJD(3) * t95;
t66 = t136 + t155;
t94 = sin(qJ(5));
t97 = cos(qJ(5));
t116 = t64 * t94 - t97 * t66;
t23 = t97 * t64 + t66 * t94;
t187 = t116 * t23;
t180 = qJD(4) + qJD(5);
t67 = t94 * t95 - t97 * t98;
t186 = t180 * t67;
t68 = t94 * t98 + t95 * t97;
t103 = t180 * t68;
t185 = t116 ^ 2 - t23 ^ 2;
t148 = qJD(5) * t97;
t149 = qJD(5) * t94;
t150 = qJD(4) * t99;
t132 = t95 * t150;
t96 = sin(qJ(3));
t109 = -t96 * t145 - t132;
t34 = qJD(1) * t109 + qJD(4) * t145;
t154 = qJD(3) * t96;
t135 = t95 * t154;
t35 = -qJD(1) * t135 + qJD(4) * t66;
t6 = -t64 * t148 - t66 * t149 + t97 * t34 - t94 * t35;
t158 = qJD(1) * t96;
t86 = qJD(4) + t158;
t84 = qJD(5) + t86;
t184 = t23 * t84 + t6;
t73 = pkin(3) * t96 - pkin(7) * t99 + qJ(2);
t53 = t73 * qJD(1);
t100 = -pkin(1) - pkin(6);
t85 = t100 * qJD(1) + qJD(2);
t72 = t96 * t85;
t56 = qJD(3) * pkin(7) + t72;
t20 = t53 * t95 + t56 * t98;
t15 = -pkin(8) * t64 + t20;
t13 = t15 * t149;
t173 = t85 * t99;
t57 = -qJD(3) * pkin(3) - t173;
t30 = pkin(4) * t64 + t57;
t183 = t23 * t30 + t13;
t104 = t116 * qJD(5) - t94 * t34 - t97 * t35;
t182 = -t116 * t84 + t104;
t120 = pkin(3) * t99 + pkin(7) * t96;
t62 = t120 * qJD(3) + qJD(2);
t44 = t62 * qJD(1);
t38 = t98 * t44;
t106 = -t20 * qJD(4) + t38;
t153 = qJD(3) * t99;
t4 = -pkin(8) * t34 + (pkin(4) * qJD(1) - t85 * t95) * t153 + t106;
t133 = t99 * t145;
t151 = qJD(4) * t98;
t152 = qJD(4) * t95;
t112 = t85 * t133 + t53 * t151 - t56 * t152 + t95 * t44;
t5 = -pkin(8) * t35 + t112;
t129 = t97 * t4 - t94 * t5;
t19 = t98 * t53 - t56 * t95;
t14 = -pkin(8) * t66 + t19;
t12 = pkin(4) * t86 + t14;
t178 = t15 * t97;
t2 = t12 * t94 + t178;
t181 = -t2 * qJD(5) + t30 * t116 + t129;
t139 = 0.2e1 * qJD(1);
t159 = t100 * t96;
t83 = t98 * t159;
t162 = t95 * t73 + t83;
t179 = pkin(7) + pkin(8);
t177 = t34 * t95;
t176 = t57 * t95;
t175 = t64 * t86;
t174 = t66 * t86;
t172 = t86 * t95;
t171 = t86 * t96;
t170 = t86 * t98;
t169 = t95 * t99;
t168 = t96 * t98;
t167 = t98 * t99;
t166 = t99 * t35;
t113 = t67 * t96;
t165 = -qJD(1) * t113 - t186;
t110 = t68 * qJD(1);
t164 = t96 * t110 + t103;
t69 = t120 * qJD(1);
t163 = t85 * t167 + t95 * t69;
t93 = t99 ^ 2;
t161 = t96 ^ 2 - t93;
t160 = t100 * t95;
t156 = qJD(3) * t84;
t101 = qJD(3) ^ 2;
t147 = t100 * t101;
t102 = qJD(1) ^ 2;
t146 = t102 * qJ(2);
t144 = -t101 - t102;
t143 = qJ(2) * qJD(3);
t142 = qJD(1) * qJD(3);
t141 = pkin(8) * t168;
t140 = t85 * t169;
t138 = t95 * t159;
t137 = t95 * t158;
t131 = t98 * t150;
t130 = qJD(2) * t139;
t128 = qJD(4) * t179;
t87 = t99 * t142;
t127 = pkin(4) - t160;
t125 = qJD(5) * t12 + t5;
t124 = t64 + t145;
t123 = -t66 + t155;
t122 = qJD(4) * t96 + qJD(1);
t121 = -t72 + (t137 + t152) * pkin(4);
t55 = t98 * t69;
t80 = t179 * t98;
t119 = qJD(5) * t80 - t140 + t55 + (pkin(4) * t99 + t141) * qJD(1) + t98 * t128;
t79 = t179 * t95;
t118 = pkin(8) * t137 + qJD(5) * t79 + t95 * t128 + t163;
t61 = t98 * t73;
t21 = -pkin(8) * t167 + t127 * t96 + t61;
t29 = -pkin(8) * t169 + t162;
t117 = t21 * t94 + t29 * t97;
t115 = qJD(1) * t93 - t171;
t114 = -pkin(7) * t153 + t57 * t96;
t111 = qJD(1) * t67;
t108 = -t131 + t135;
t107 = -qJD(4) * t138 + t100 * t133 + t73 * t151 + t95 * t62;
t90 = -pkin(4) * t98 - pkin(3);
t82 = t96 * t87;
t63 = (pkin(4) * t95 - t100) * t99;
t49 = t98 * t62;
t46 = t67 * t99;
t45 = t68 * t99;
t36 = -pkin(4) * t108 + t100 * t154;
t17 = pkin(4) * t35 + t85 * t154;
t11 = -t149 * t169 + (t180 * t167 - t135) * t97 + t109 * t94;
t10 = qJD(3) * t113 - t103 * t99;
t9 = pkin(8) * t108 + t107;
t8 = t49 + (-t83 + (pkin(8) * t99 - t73) * t95) * qJD(4) + (t127 * t99 + t141) * qJD(3);
t1 = t12 * t97 - t15 * t94;
t3 = [0, 0, 0, 0, t130, qJ(2) * t130, -0.2e1 * t82, 0.2e1 * t161 * t142, -t101 * t96, -t101 * t99, 0, -t96 * t147 + (qJD(2) * t96 + t99 * t143) * t139, -t99 * t147 + (qJD(2) * t99 - t96 * t143) * t139, t109 * t66 + t34 * t167, (t64 * t98 + t66 * t95) * t154 + (-t177 - t35 * t98 + (t64 * t95 - t66 * t98) * qJD(4)) * t99, -t86 * t132 + t34 * t96 + (t115 * t98 + t66 * t99) * qJD(3), -t86 * t131 - t35 * t96 + (-t115 * t95 - t64 * t99) * qJD(3), t86 * t153 + t82, -t100 * t166 + t38 * t96 + t49 * t86 + (-t162 * t86 + t57 * t167 - t20 * t96) * qJD(4) + ((t100 * t64 - t176) * t96 + (-t86 * t160 + (t61 - t138) * qJD(1) + t19) * t99) * qJD(3), -t107 * t86 - t112 * t96 + (-t100 * t34 - t57 * t152) * t99 + ((-t162 * qJD(1) - t20) * t99 + (t100 * t66 + (-t57 + t173) * t98) * t96) * qJD(3), -t10 * t116 - t46 * t6, -t10 * t23 - t104 * t46 + t11 * t116 - t45 * t6, t10 * t84 + t6 * t96 + (-qJD(1) * t46 - t116) * t153, -t11 * t84 + t104 * t96 + (-qJD(1) * t45 - t23) * t153, t84 * t153 + t82, (t97 * t8 - t94 * t9) * t84 + t129 * t96 + t36 * t23 - t63 * t104 + t17 * t45 + t30 * t11 + (-t117 * t84 - t2 * t96) * qJD(5) + ((t21 * t97 - t29 * t94) * qJD(1) + t1) * t153, t30 * t10 + t13 * t96 - t17 * t46 - t36 * t116 + t63 * t6 + (-(-qJD(5) * t29 + t8) * t84 - t4 * t96) * t94 + (-(qJD(5) * t21 + t9) * t84 - t125 * t96) * t97 + (-qJD(1) * t117 - t2) * t153; 0, 0, 0, 0, -t102, -t146, 0, 0, 0, 0, 0, t144 * t96, t144 * t99, 0, 0, 0, 0, 0, -t166 - t122 * t170 + (t64 * t96 + (-t86 - t158) * t169) * qJD(3), -t34 * t99 + t122 * t172 + (-t86 * t167 + (t66 - t136) * t96) * qJD(3), 0, 0, 0, 0, 0, t84 * t111 + (-t68 * t156 + t104) * t99 + ((-t68 * t157 + t23) * qJD(3) + t84 * t186) * t96, t84 * t110 + (t67 * t156 - t6) * t99 + (t103 * t84 + (t111 * t99 - t116) * qJD(3)) * t96; 0, 0, 0, 0, 0, 0, t99 * t102 * t96, -t161 * t102, 0, 0, 0, -t99 * t146, t96 * t146, t66 * t170 + t177, (t34 - t175) * t98 + (-t35 - t174) * t95, t86 * t151 + (t123 * t99 + t86 * t168) * qJD(1), -t86 * t152 + (t124 * t99 - t95 * t171) * qJD(1), -t86 * t157, -pkin(3) * t35 - t55 * t86 + (-t124 * t96 + t86 * t169) * t85 + (-pkin(7) * t170 + t176) * qJD(4) + (t114 * t95 - t19 * t99) * qJD(1), -pkin(3) * t34 + t163 * t86 + t123 * t72 + (pkin(7) * t172 + t57 * t98) * qJD(4) + (t114 * t98 + t20 * t99) * qJD(1), -t116 * t165 + t6 * t68, t104 * t68 + t116 * t164 - t165 * t23 - t6 * t67, t165 * t84 + (qJD(3) * t68 + t116) * t157, -t164 * t84 + (-qJD(3) * t67 + t23) * t157, -t84 * t157, t17 * t67 - t90 * t104 + (t118 * t94 - t119 * t97) * t84 + t164 * t30 + t121 * t23 + ((-t79 * t97 - t80 * t94) * qJD(3) - t1) * t157, t17 * t68 + t90 * t6 + (t118 * t97 + t119 * t94) * t84 + t165 * t30 - t121 * t116 + (-(-t79 * t94 + t80 * t97) * qJD(3) + t2) * t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t64, -t64 ^ 2 + t66 ^ 2, t34 + t175, t174 - t35, t87, -qJD(3) * t140 + t20 * t86 - t57 * t66 + t106, t19 * t86 + t57 * t64 - t112, -t187, t185, t184, t182, t87, -(-t14 * t94 - t178) * t84 + (-t149 * t84 - t23 * t66 + t87 * t97) * pkin(4) + t181, (-t15 * t84 - t4) * t94 + (t14 * t84 - t125) * t97 + (t116 * t66 - t148 * t84 - t87 * t94) * pkin(4) + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187, t185, t184, t182, t87, t2 * t84 + t181, t1 * t84 - t125 * t97 - t94 * t4 + t183;];
tauc_reg = t3;
