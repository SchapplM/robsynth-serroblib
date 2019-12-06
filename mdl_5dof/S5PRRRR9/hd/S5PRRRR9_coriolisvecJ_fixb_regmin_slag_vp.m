% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:10
% EndTime: 2019-12-05 17:21:19
% DurationCPUTime: 2.26s
% Computational Cost: add. (1701->292), mult. (4428->443), div. (0->0), fcn. (3306->10), ass. (0->162)
t113 = sin(qJ(5));
t117 = cos(qJ(5));
t114 = sin(qJ(4));
t118 = cos(qJ(4));
t165 = t118 * qJD(3);
t115 = sin(qJ(3));
t174 = qJD(2) * t115;
t79 = t114 * t174 - t165;
t172 = qJD(3) * t114;
t81 = t118 * t174 + t172;
t133 = t113 * t79 - t117 * t81;
t28 = t113 * t81 + t117 * t79;
t219 = t133 * t28;
t169 = qJD(4) * t114;
t151 = t115 * t169;
t119 = cos(qJ(3));
t153 = t119 * t165;
t218 = t153 - t151;
t217 = t133 ^ 2 - t28 ^ 2;
t112 = cos(pkin(5));
t176 = qJD(1) * t119;
t116 = sin(qJ(2));
t111 = sin(pkin(5));
t177 = qJD(1) * t111;
t149 = t116 * t177;
t90 = qJD(2) * pkin(7) + t149;
t56 = t112 * t176 - t115 * t90;
t48 = -qJD(3) * pkin(3) - t56;
t24 = pkin(4) * t79 + t48;
t120 = cos(qJ(2));
t159 = t120 * t177;
t92 = -pkin(3) * t119 - pkin(8) * t115 - pkin(2);
t59 = qJD(2) * t92 - t159;
t197 = t114 * t59;
t186 = t112 * t115;
t101 = qJD(1) * t186;
t57 = t119 * t90 + t101;
t49 = qJD(3) * pkin(8) + t57;
t21 = t118 * t49 + t197;
t13 = -pkin(9) * t79 + t21;
t167 = qJD(5) * t113;
t9 = t13 * t167;
t216 = t24 * t28 + t9;
t166 = qJD(5) * t117;
t163 = qJD(3) * qJD(4);
t51 = qJD(2) * t218 + t118 * t163;
t168 = qJD(4) * t118;
t150 = t115 * t168;
t170 = qJD(3) * t119;
t154 = t114 * t170;
t127 = t150 + t154;
t52 = qJD(2) * t127 + t114 * t163;
t6 = -t113 * t52 + t117 * t51 - t166 * t79 - t167 * t81;
t173 = qJD(2) * t119;
t102 = -qJD(4) + t173;
t99 = -qJD(5) + t102;
t215 = -t28 * t99 + t6;
t123 = qJD(5) * t133 - t113 * t51 - t117 * t52;
t214 = t133 * t99 + t123;
t164 = qJD(2) * qJD(3);
t104 = t115 * t164;
t175 = qJD(2) * t111;
t156 = t120 * t175;
t171 = qJD(3) * t115;
t25 = -t90 * t171 + (qJD(3) * t112 + t156) * t176;
t138 = pkin(3) * t115 - pkin(8) * t119;
t87 = t138 * qJD(3);
t55 = (t87 + t149) * qJD(2);
t145 = t114 * t25 - t118 * t55;
t124 = -qJD(4) * t21 - t145;
t2 = pkin(4) * t104 - pkin(9) * t51 + t124;
t129 = t114 * t55 + t118 * t25 + t168 * t59 - t169 * t49;
t5 = -pkin(9) * t52 + t129;
t161 = -t113 * t5 + t117 * t2;
t196 = t117 * t13;
t20 = -t114 * t49 + t118 * t59;
t12 = -pkin(9) * t81 + t20;
t8 = -pkin(4) * t102 + t12;
t4 = t113 * t8 + t196;
t213 = -qJD(5) * t4 + t24 * t133 + t161;
t83 = t113 * t118 + t114 * t117;
t62 = t83 * t115;
t181 = t119 * t120;
t208 = pkin(7) * t114;
t212 = (-t114 * t181 + t116 * t118) * t177 - t118 * t87 - t171 * t208;
t211 = -(t114 * t116 + t118 * t181) * t177 + t114 * t87 + t92 * t168;
t210 = qJD(4) + qJD(5);
t209 = pkin(8) + pkin(9);
t182 = t118 * t119;
t103 = pkin(7) * t182;
t131 = pkin(4) * t115 - pkin(9) * t182;
t207 = -t131 * qJD(3) - (-t103 + (pkin(9) * t115 - t92) * t114) * qJD(4) + t212;
t206 = -t127 * pkin(9) + (-t115 * t165 - t119 * t169) * pkin(7) + t211;
t82 = t113 * t114 - t117 * t118;
t128 = t82 * t119;
t205 = qJD(2) * t128 - t210 * t82;
t204 = (-t173 + t210) * t83;
t84 = t138 * qJD(2);
t203 = t114 * t84 + t118 * t56;
t201 = qJD(2) * pkin(2);
t200 = t102 * t79;
t199 = t102 * t81;
t198 = t114 * t48;
t143 = t115 * t156;
t26 = qJD(1) * t143 + qJD(3) * t101 + t170 * t90;
t195 = t26 * t114;
t194 = t26 * t118;
t193 = t51 * t114;
t191 = t114 * t92 + t103;
t190 = t102 * t118;
t189 = t111 * t116;
t188 = t111 * t120;
t122 = qJD(2) ^ 2;
t187 = t111 * t122;
t185 = t114 * t115;
t184 = t114 * t119;
t183 = t115 * t118;
t121 = qJD(3) ^ 2;
t180 = t121 * t115;
t179 = t121 * t119;
t109 = t115 ^ 2;
t178 = -t119 ^ 2 + t109;
t162 = t116 * t187;
t160 = qJD(4) * t209;
t157 = t116 * t175;
t155 = t114 * t173;
t152 = t102 * t169;
t148 = qJD(5) * t8 + t5;
t144 = -t114 * t56 + t118 * t84;
t142 = t119 * t156;
t141 = -t57 + (-t155 + t169) * pkin(4);
t96 = t209 * t118;
t140 = qJD(2) * t131 + qJD(5) * t96 + t118 * t160 + t144;
t95 = t209 * t114;
t139 = pkin(9) * t155 - qJD(5) * t95 - t114 * t160 - t203;
t91 = -t159 - t201;
t137 = -t91 - t159;
t78 = t118 * t92;
t32 = -pkin(9) * t183 + t78 + (-pkin(4) - t208) * t119;
t40 = -pkin(9) * t185 + t191;
t136 = t113 * t32 + t117 * t40;
t69 = t119 * t189 + t186;
t130 = t114 * t188 - t118 * t69;
t38 = -t114 * t69 - t118 * t188;
t135 = t113 * t130 + t117 * t38;
t134 = t113 * t38 - t117 * t130;
t132 = qJD(2) * t109 - t102 * t119;
t68 = -t112 * t119 + t115 * t189;
t126 = qJD(3) * (-t137 - t201);
t107 = -pkin(4) * t118 - pkin(3);
t88 = (pkin(4) * t114 + pkin(7)) * t115;
t63 = t82 * t115;
t58 = pkin(4) * t127 + pkin(7) * t170;
t37 = qJD(3) * t69 + t143;
t36 = -qJD(3) * t68 + t142;
t18 = -t167 * t185 + (t183 * t210 + t154) * t117 + t218 * t113;
t17 = -qJD(3) * t128 - t210 * t62;
t16 = pkin(4) * t52 + t26;
t11 = qJD(4) * t38 + t114 * t157 + t118 * t36;
t10 = qJD(4) * t130 - t114 * t36 + t118 * t157;
t3 = -t113 * t13 + t117 * t8;
t1 = [0, 0, -t162, -t120 * t187, 0, 0, 0, 0, 0, -t119 * t162 + (-t37 - t143) * qJD(3), t115 * t162 + (-t36 - t142) * qJD(3), 0, 0, 0, 0, 0, -t10 * t102 + t104 * t38 + t37 * t79 + t52 * t68, t102 * t11 + t104 * t130 + t37 * t81 + t51 * t68, 0, 0, 0, 0, 0, t37 * t28 - t68 * t123 - (-qJD(5) * t134 + t117 * t10 - t113 * t11) * t99 + t135 * t104, -t37 * t133 + t68 * t6 + (qJD(5) * t135 + t113 * t10 + t117 * t11) * t99 - t134 * t104; 0, 0, 0, 0, 0.2e1 * t119 * t104, -0.2e1 * t178 * t164, t179, -t180, 0, -pkin(7) * t179 + t115 * t126, pkin(7) * t180 + t119 * t126, t81 * t153 + (t118 * t51 - t169 * t81) * t115, (-t114 * t81 - t118 * t79) * t170 + (-t193 - t118 * t52 + (t114 * t79 - t118 * t81) * qJD(4)) * t115, t102 * t151 - t119 * t51 + (t115 * t81 + t118 * t132) * qJD(3), t102 * t150 + t119 * t52 + (-t114 * t132 - t115 * t79) * qJD(3), (-t102 - t173) * t171, (t169 * t92 + t212) * t102 + ((pkin(7) * t79 + t198) * qJD(3) + (t197 + (pkin(7) * t102 + t49) * t118) * qJD(4) + t145) * t119 + (-t79 * t159 + t48 * t168 + pkin(7) * t52 + t195 + ((-pkin(7) * t184 + t78) * qJD(2) + t20) * qJD(3)) * t115, t211 * t102 + (t48 * t165 + (qJD(3) * t81 - t152) * pkin(7) + t129) * t119 + (-t81 * t159 - t48 * t169 + pkin(7) * t51 + t194 + (-pkin(7) * t190 - qJD(2) * t191 - t21) * qJD(3)) * t115, -t133 * t17 - t6 * t63, -t123 * t63 + t133 * t18 - t17 * t28 - t6 * t62, -t119 * t6 - t17 * t99 + (-qJD(2) * t63 - t133) * t171, -t119 * t123 + t18 * t99 + (-qJD(2) * t62 - t28) * t171, (-t99 - t173) * t171, t58 * t28 - t88 * t123 + t16 * t62 + t24 * t18 - t161 * t119 + (t113 * t206 + t117 * t207) * t99 + (t119 * t4 + t136 * t99) * qJD(5) + (-t28 * t159 + ((-t113 * t40 + t117 * t32) * qJD(2) + t3) * qJD(3)) * t115, -t9 * t119 - t16 * t63 + t24 * t17 - t58 * t133 + t88 * t6 + (t2 * t119 + (-qJD(5) * t40 - t207) * t99) * t113 + (t148 * t119 + (qJD(5) * t32 + t206) * t99) * t117 + (t133 * t159 + (-qJD(2) * t136 - t4) * qJD(3)) * t115; 0, 0, 0, 0, -t115 * t122 * t119, t178 * t122, 0, 0, 0, qJD(3) * t57 - t174 * t91 - t26, t137 * t173, -t190 * t81 + t193, (t51 + t200) * t118 + (-t52 + t199) * t114, -t102 * t168 + (t102 * t182 + (-t81 + t172) * t115) * qJD(2), t152 + (-t102 * t184 + (t79 + t165) * t115) * qJD(2), t102 * t174, -pkin(3) * t52 - t194 + t144 * t102 - t57 * t79 + (pkin(8) * t190 + t198) * qJD(4) + (-t20 * t115 + (-pkin(8) * t171 - t119 * t48) * t114) * qJD(2), -pkin(3) * t51 + t195 - t203 * t102 - t57 * t81 + (-pkin(8) * t102 * t114 + t118 * t48) * qJD(4) + (-t48 * t182 + (-pkin(8) * t165 + t21) * t115) * qJD(2), -t133 * t205 + t6 * t83, t123 * t83 + t133 * t204 - t205 * t28 - t6 * t82, -t205 * t99 + (qJD(3) * t83 + t133) * t174, t204 * t99 + (-qJD(3) * t82 + t28) * t174, t99 * t174, -t107 * t123 + t16 * t82 + (t113 * t139 + t117 * t140) * t99 + t141 * t28 + t204 * t24 + ((-t113 * t96 - t117 * t95) * qJD(3) - t3) * t174, t107 * t6 + t16 * t83 + (-t113 * t140 + t117 * t139) * t99 - t141 * t133 + t205 * t24 + (-(-t113 * t95 + t117 * t96) * qJD(3) + t4) * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t79, -t79 ^ 2 + t81 ^ 2, t51 - t200, -t199 - t52, t104, -t102 * t21 - t48 * t81 + t124, -t102 * t20 + t48 * t79 - t129, -t219, t217, t215, t214, t104, (-t113 * t12 - t196) * t99 + (t104 * t117 + t167 * t99 - t28 * t81) * pkin(4) + t213, (t13 * t99 - t2) * t113 + (-t12 * t99 - t148) * t117 + (-t104 * t113 + t133 * t81 + t166 * t99) * pkin(4) + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, t217, t215, t214, t104, -t4 * t99 + t213, -t113 * t2 - t117 * t148 - t3 * t99 + t216;];
tauc_reg = t1;
