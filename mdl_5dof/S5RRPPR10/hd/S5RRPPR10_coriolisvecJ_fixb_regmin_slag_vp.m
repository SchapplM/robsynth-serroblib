% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:40
% EndTime: 2019-12-31 19:44:47
% DurationCPUTime: 2.13s
% Computational Cost: add. (1438->304), mult. (3750->426), div. (0->0), fcn. (2439->6), ass. (0->166)
t120 = cos(qJ(2));
t165 = t120 * qJD(1);
t105 = qJD(5) + t165;
t206 = qJD(5) - t105;
t118 = sin(qJ(2));
t115 = sin(pkin(8));
t116 = cos(pkin(8));
t117 = sin(qJ(5));
t119 = cos(qJ(5));
t80 = t119 * t115 - t117 * t116;
t58 = t80 * t118;
t163 = qJD(1) * qJD(2);
t205 = -0.2e1 * t163;
t174 = qJD(1) * t118;
t154 = t116 * t174;
t168 = t115 * qJD(2);
t77 = t154 + t168;
t73 = t77 ^ 2;
t155 = t115 * t174;
t166 = t116 * qJD(2);
t75 = t155 - t166;
t204 = -t75 ^ 2 - t73;
t164 = t120 * qJD(4);
t173 = qJD(2) * t118;
t203 = qJ(4) * t173 - t164;
t152 = t120 * t163;
t145 = t116 * t152;
t32 = t117 * t75 + t119 * t77;
t98 = t115 * t152;
t8 = qJD(5) * t32 + t117 * t145 - t119 * t98;
t202 = pkin(3) + pkin(4);
t201 = -pkin(7) + qJ(3);
t79 = t117 * t115 + t119 * t116;
t129 = t79 * t120;
t200 = -qJD(1) * t129 - t79 * qJD(5);
t156 = t115 * t165;
t183 = t117 * t120;
t159 = t116 * t183;
t169 = qJD(5) * t119;
t170 = qJD(5) * t117;
t199 = -qJD(1) * t159 + t115 * t169 - t116 * t170 + t119 * t156;
t139 = pkin(2) * t118 - qJ(3) * t120;
t63 = qJD(2) * t139 - t118 * qJD(3);
t54 = t63 * qJD(1);
t108 = pkin(6) * t174;
t87 = (qJD(3) - t108) * qJD(2);
t23 = t115 * t54 + t116 * t87;
t182 = t118 * qJ(3);
t90 = -t120 * pkin(2) - pkin(1) - t182;
t69 = t90 * qJD(1);
t109 = pkin(6) * t165;
t96 = qJD(2) * qJ(3) + t109;
t35 = t115 * t69 + t116 * t96;
t198 = qJD(2) * pkin(2);
t197 = t116 * t63;
t83 = t139 * qJD(1);
t196 = t116 * t83;
t34 = -t115 * t96 + t116 * t69;
t25 = pkin(3) * t165 + qJD(4) - t34;
t195 = t118 * t25;
t27 = -qJ(4) * t165 + t35;
t194 = t118 * t27;
t30 = t117 * t77 - t119 * t75;
t193 = t30 * t105;
t192 = t32 * t105;
t185 = t116 * t120;
t103 = pkin(6) * t185;
t51 = t115 * t90 + t103;
t104 = pkin(6) * t152;
t191 = -pkin(3) * t98 - t104;
t190 = qJ(4) * t116;
t189 = qJD(3) * t77;
t188 = t115 * t118;
t187 = t115 * t120;
t186 = t116 * t118;
t180 = t119 * t120;
t122 = qJD(1) ^ 2;
t179 = t120 * t122;
t121 = qJD(2) ^ 2;
t178 = t121 * t118;
t177 = t121 * t120;
t176 = t77 * qJD(4);
t114 = t120 ^ 2;
t175 = t118 ^ 2 - t114;
t172 = qJD(2) * t120;
t171 = qJD(3) * t116;
t167 = t115 * qJD(4);
t162 = pkin(7) * t185;
t102 = pkin(6) * t187;
t151 = t118 * t163;
t161 = qJ(4) * t151 + t23;
t160 = pkin(6) * t173;
t158 = -pkin(6) * t115 - pkin(3);
t22 = -t115 * t87 + t116 * t54;
t4 = (-t202 * t118 - t162) * t163 - t22;
t5 = (pkin(7) * t168 - qJD(4)) * t165 + t161;
t157 = -t117 * t5 + t119 * t4;
t153 = qJD(4) * t186;
t150 = t115 * qJ(4) + pkin(2);
t147 = -qJD(3) + t198;
t141 = -t108 + t147;
t149 = t141 - t198;
t50 = t116 * t90 - t102;
t148 = pkin(1) * t205;
t134 = -t202 * t115 + t190;
t146 = -t134 * t165 + t109 + t167;
t144 = -t77 * t165 + t98;
t44 = -t120 * qJ(4) + t51;
t143 = t118 * t158;
t142 = t117 * t4 + t119 * t5;
t11 = t75 * pkin(7) + t27;
t6 = pkin(4) * t165 - t77 * pkin(7) + t25;
t2 = t119 * t11 + t117 * t6;
t140 = t117 * t11 - t119 * t6;
t138 = pkin(3) * t115 - t190;
t112 = t120 * pkin(3);
t28 = t120 * pkin(4) + t102 + t112 + (-pkin(7) * t118 - t90) * t116;
t33 = pkin(7) * t188 + t44;
t137 = t117 * t33 - t119 * t28;
t136 = t117 * t28 + t119 * t33;
t67 = t115 * t83;
t47 = -pkin(6) * t154 + t67;
t55 = t115 * t63;
t41 = -t116 * t160 + t55;
t135 = pkin(6) + t138;
t124 = -t162 + (-pkin(4) + t158) * t118;
t94 = t201 * t116;
t133 = -t124 * qJD(1) + qJD(3) * t115 - qJD(5) * t94 + t196;
t106 = qJ(4) * t174;
t130 = -pkin(6) * t186 + pkin(7) * t187;
t93 = t201 * t115;
t132 = qJD(1) * t130 - qJD(5) * t93 + t106 - t171 + t67;
t131 = -t117 * t98 - t119 * t145 - t75 * t169 + t77 * t170;
t59 = t79 * t118;
t127 = -pkin(6) + t134;
t126 = t77 * qJ(4) + t141;
t125 = -qJ(4) * t145 - t191;
t97 = qJD(3) * t156;
t88 = -t116 * pkin(3) - t150;
t68 = t202 * t116 + t150;
t61 = t75 * t165;
t60 = t75 * t171;
t56 = t135 * t118;
t49 = t138 * t165 + t109;
t46 = pkin(6) * t155 + t196;
t45 = t112 - t50;
t43 = t127 * t118;
t42 = t61 + t145;
t40 = t115 * t160 + t197;
t39 = qJD(1) * t143 - t196;
t38 = t106 + t47;
t37 = t135 * t172 - t153;
t29 = qJD(2) * t143 - t197;
t24 = t127 * t172 + t153;
t21 = t41 + t203;
t20 = t75 * pkin(3) - t126;
t19 = qJD(2) * t129 + qJD(5) * t58;
t18 = qJD(2) * t159 + qJD(5) * t59 - t168 * t180;
t17 = t125 - t176;
t15 = -pkin(3) * t151 - t22;
t14 = qJD(2) * t130 + t203 + t55;
t13 = t124 * qJD(2) - t197;
t12 = t176 + (-pkin(4) * t115 + t190) * t152 + t191;
t10 = -qJD(1) * t164 + t161;
t9 = -t202 * t75 + t126;
t1 = [0, 0, 0, 0.2e1 * t120 * t151, t175 * t205, t177, -t178, 0, -pkin(6) * t177 + t118 * t148, pkin(6) * t178 + t120 * t148, (-qJD(1) * t40 - t22) * t120 + ((pkin(6) * t75 - t115 * t141) * t120 + (t34 + (t50 + 0.2e1 * t102) * qJD(1)) * t118) * qJD(2), (qJD(1) * t41 + t23) * t120 + ((pkin(6) * t77 - t116 * t141) * t120 + (-t35 + (-t51 + 0.2e1 * t103) * qJD(1)) * t118) * qJD(2), -t40 * t77 - t41 * t75 + (-t115 * t23 - t116 * t22) * t118 + (-t115 * t35 - t116 * t34 + (-t115 * t51 - t116 * t50) * qJD(1)) * t172, t22 * t50 + t23 * t51 + t34 * t40 + t35 * t41 + (-t141 + t108) * pkin(6) * t172, t17 * t188 + t37 * t75 + (qJD(1) * t29 + t15) * t120 + (t20 * t187 - t195 + (-t118 * t45 + t187 * t56) * qJD(1)) * qJD(2), -t21 * t75 + t29 * t77 + (-t10 * t115 + t116 * t15) * t118 + (-t115 * t27 + t116 * t25 + (-t115 * t44 + t116 * t45) * qJD(1)) * t172, -t17 * t186 - t37 * t77 + (-qJD(1) * t21 - t10) * t120 + (-t20 * t185 + t194 + (t118 * t44 - t185 * t56) * qJD(1)) * qJD(2), t10 * t44 + t15 * t45 + t17 * t56 + t20 * t37 + t27 * t21 + t25 * t29, -t131 * t59 + t32 * t19, -t131 * t58 - t32 * t18 - t19 * t30 - t59 * t8, t19 * t105 - t131 * t120 + (-qJD(1) * t59 - t32) * t173, -t18 * t105 - t8 * t120 + (-qJD(1) * t58 + t30) * t173, (-t105 - t165) * t173, (-t117 * t14 + t119 * t13) * t105 + t157 * t120 + t24 * t30 + t43 * t8 - t12 * t58 + t9 * t18 + (-t105 * t136 - t120 * t2) * qJD(5) + (qJD(1) * t137 + t140) * t173, -(t117 * t13 + t119 * t14) * t105 - t142 * t120 + t24 * t32 - t43 * t131 + t12 * t59 + t9 * t19 + (t105 * t137 + t120 * t140) * qJD(5) + (qJD(1) * t136 + t2) * t173; 0, 0, 0, -t118 * t179, t175 * t122, 0, 0, 0, t122 * pkin(1) * t118, pkin(1) * t179, t97 + ((-qJ(3) * t168 - t34) * t118 + (t46 + t149 * t115 + (-t75 - t166) * pkin(6)) * t120) * qJD(1), ((-qJ(3) * t166 + t35) * t118 + (-t47 + (-t77 + t168) * pkin(6) + (-t147 + t141) * t116) * t120) * qJD(1), t46 * t77 + t47 * t75 - t60 + (t34 * t165 + t23) * t116 + (t35 * t165 + t189 - t22) * t115, -t34 * t46 - t35 * t47 + (-t115 * t34 + t116 * t35) * qJD(3) + (-t22 * t115 + t23 * t116) * qJ(3) + t149 * t109, -t17 * t116 + t97 + (-t49 - t167) * t75 + (t195 - t120 * t39 + (-t120 * t20 + (t120 * t88 - t182) * qJD(2)) * t115) * qJD(1), t38 * t75 - t39 * t77 - t60 + (-t165 * t25 + t10) * t116 + (t165 * t27 + t15 + t189) * t115, t49 * t77 + (-t17 + t176) * t115 + (-t194 + t120 * t38 + (qJ(3) * t173 + (-qJD(2) * t88 - qJD(3) + t20) * t120) * t116) * qJD(1), t17 * t88 - t20 * t49 - t25 * t39 - t27 * t38 + (qJ(3) * t10 + qJD(3) * t27) * t116 + (qJ(3) * t15 + qJD(3) * t25 - qJD(4) * t20) * t115, -t131 * t80 + t200 * t32, t131 * t79 - t199 * t32 - t200 * t30 - t80 * t8, t200 * t105 + (-qJD(2) * t80 + t32) * t174, -t199 * t105 + (qJD(2) * t79 - t30) * t174, t105 * t174, t12 * t79 + t68 * t8 + t199 * t9 + t146 * t30 + (t117 * t132 + t119 * t133) * t105 + (-(-t117 * t94 + t119 * t93) * qJD(2) - t140) * t174, t12 * t80 - t68 * t131 + t200 * t9 + t146 * t32 + (-t117 * t133 + t119 * t132) * t105 + ((t117 * t93 + t119 * t94) * qJD(2) - t2) * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t42, t204, t34 * t77 + t35 * t75 + t104, t144, t204, -t42, t27 * t75 + (-qJD(4) - t25) * t77 + t125, 0, 0, 0, 0, 0, -t8 - t192, t131 + t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t75 - t151, -t61 + t145, -t114 * t122 - t73, t20 * t77 + (-pkin(3) * t173 + t120 * t27) * qJD(1) - t22, 0, 0, 0, 0, 0, -t105 * t170 - t77 * t30 + (-t105 * t183 - t119 * t173) * qJD(1), -t105 * t169 - t77 * t32 + (-t105 * t180 + t117 * t173) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t30, -t30 ^ 2 + t32 ^ 2, -t131 + t193, -t8 + t192, -t151, -t206 * t2 - t9 * t32 + t157, t206 * t140 + t9 * t30 - t142;];
tauc_reg = t1;
