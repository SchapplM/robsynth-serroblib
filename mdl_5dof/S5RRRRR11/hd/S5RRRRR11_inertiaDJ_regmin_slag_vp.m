% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:42
% EndTime: 2019-12-31 22:43:51
% DurationCPUTime: 2.79s
% Computational Cost: add. (3329->313), mult. (9355->589), div. (0->0), fcn. (8982->10), ass. (0->165)
t128 = sin(qJ(3));
t214 = -0.4e1 * t128;
t126 = sin(qJ(5));
t127 = sin(qJ(4));
t130 = cos(qJ(5));
t131 = cos(qJ(4));
t100 = t126 * t131 + t130 * t127;
t87 = t100 * t128;
t132 = cos(qJ(3));
t182 = qJD(3) * t132;
t154 = t131 * t182;
t180 = qJD(4) * t127;
t213 = -t128 * t180 + t154;
t122 = t131 ^ 2;
t187 = t127 ^ 2 - t122;
t151 = t187 * qJD(4);
t212 = qJD(4) + qJD(5);
t146 = -t132 * pkin(3) - t128 * pkin(9);
t104 = -pkin(2) + t146;
t178 = qJD(4) * t132;
t162 = t127 * t178;
t183 = qJD(3) * t131;
t136 = t128 * t183 + t162;
t145 = pkin(3) * t128 - pkin(9) * t132;
t139 = t145 * qJD(3);
t179 = qJD(4) * t131;
t54 = pkin(8) * t136 - t104 * t179 - t127 * t139;
t124 = sin(pkin(5));
t211 = 0.2e1 * t124;
t210 = pkin(9) + pkin(10);
t125 = cos(pkin(5));
t129 = sin(qJ(2));
t196 = t124 * t129;
t89 = -t125 * t132 + t128 * t196;
t209 = t89 * pkin(4);
t208 = pkin(1) * t129;
t207 = pkin(4) * t130;
t206 = pkin(8) * t124;
t205 = pkin(8) * t127;
t133 = cos(qJ(2));
t81 = pkin(7) * t196 + (-pkin(1) * t133 - pkin(2)) * t125;
t90 = t125 * t128 + t132 * t196;
t50 = t89 * pkin(3) - t90 * pkin(9) + t81;
t195 = t124 * t133;
t172 = pkin(7) * t195;
t82 = t172 + (pkin(8) + t208) * t125;
t83 = (-pkin(2) * t133 - pkin(8) * t129 - pkin(1)) * t124;
t204 = t128 * t83 + t132 * t82;
t52 = -pkin(9) * t195 + t204;
t21 = t127 * t50 + t131 * t52;
t64 = t90 * t127 + t131 * t195;
t19 = -t64 * pkin(10) + t21;
t202 = t130 * t19;
t192 = t127 * t128;
t188 = t131 * t132;
t114 = pkin(8) * t188;
t80 = t127 * t104 + t114;
t66 = -pkin(10) * t192 + t80;
t201 = t130 * t66;
t185 = qJD(2) * t129;
t160 = t124 * t185;
t184 = qJD(2) * t133;
t159 = t124 * t184;
t197 = qJD(3) * t89;
t63 = t132 * t159 - t197;
t32 = -qJD(4) * t64 + t127 * t160 + t63 * t131;
t200 = t32 * t127;
t199 = t32 * t131;
t118 = qJD(3) * t128;
t158 = t127 * t118;
t198 = pkin(8) * t158 + t131 * t139;
t107 = t210 * t127;
t194 = t126 * t107;
t191 = t128 * t131;
t189 = t130 * t131;
t121 = t128 ^ 2;
t186 = -t132 ^ 2 + t121;
t181 = qJD(3) * t133;
t177 = qJD(5) * t126;
t176 = qJD(5) * t130;
t175 = -0.2e1 * pkin(2) * qJD(3);
t174 = -0.2e1 * pkin(3) * qJD(4);
t173 = t132 * t205;
t171 = pkin(4) * t180;
t170 = pkin(8) * t182;
t169 = pkin(4) * t177;
t168 = pkin(4) * t176;
t167 = t127 * t195;
t62 = qJD(3) * t90 + t128 * t159;
t86 = (t125 * t208 + t172) * qJD(2);
t134 = t62 * pkin(3) - t63 * pkin(9) + t86;
t84 = (pkin(2) * t129 - pkin(8) * t133) * t124 * qJD(2);
t85 = -t125 * pkin(1) * t184 + pkin(7) * t160;
t25 = t82 * t118 - t128 * t84 + t132 * t85 - t83 * t182;
t23 = pkin(9) * t160 - t25;
t11 = -t21 * qJD(4) - t127 * t23 + t131 * t134;
t5 = t62 * pkin(4) - t32 * pkin(10) + t11;
t10 = -t127 * t134 - t131 * t23 - t50 * t179 + t52 * t180;
t31 = -qJD(4) * t167 + t63 * t127 - t131 * t160 + t90 * t179;
t7 = -t31 * pkin(10) - t10;
t166 = -t126 * t7 + t130 * t5;
t119 = t124 ^ 2;
t164 = t119 * t184;
t161 = t131 * t178;
t157 = t127 * t182;
t156 = t127 * t179;
t155 = t128 * t182;
t35 = (pkin(4) * t128 - pkin(10) * t188) * qJD(3) + (-t114 + (pkin(10) * t128 - t104) * t127) * qJD(4) + t198;
t135 = t128 * t179 + t157;
t41 = -pkin(10) * t135 - t54;
t153 = -t126 * t41 + t130 * t35;
t20 = -t127 * t52 + t131 * t50;
t152 = -t128 * t82 + t132 * t83;
t150 = t186 * qJD(3);
t149 = 0.2e1 * t155;
t148 = t129 * t164;
t147 = t127 * t154;
t51 = pkin(3) * t195 - t152;
t65 = t90 * t131 - t167;
t18 = -t65 * pkin(10) + t20 + t209;
t9 = t126 * t18 + t202;
t98 = t131 * t104;
t59 = -pkin(10) * t191 + t98 + (-pkin(4) - t205) * t132;
t37 = t126 * t59 + t201;
t39 = t126 * t65 + t130 * t64;
t40 = -t126 * t64 + t130 * t65;
t144 = -t127 * t65 - t131 * t64;
t26 = -t83 * t118 + t128 * t85 + t132 * t84 - t82 * t182;
t108 = t210 * t131;
t68 = t130 * t108 - t194;
t99 = t126 * t127 - t189;
t1 = -t126 * t5 - t130 * t7 - t18 * t176 + t19 * t177;
t24 = -pkin(3) * t160 - t26;
t143 = t24 * t127 + t51 * t179;
t142 = -t24 * t131 + t51 * t180;
t141 = t127 * t62 + t89 * t179;
t140 = -t131 * t62 + t89 * t180;
t14 = -t126 * t35 - t130 * t41 - t59 * t176 + t66 * t177;
t138 = t128 * t181 + t132 * t185;
t137 = t128 * t185 - t132 * t181;
t117 = -t131 * pkin(4) - pkin(3);
t111 = -0.2e1 * t155;
t101 = (pkin(4) * t127 + pkin(8)) * t128;
t88 = t99 * t128;
t79 = t98 - t173;
t75 = pkin(4) * t135 + t170;
t67 = -t130 * t107 - t126 * t108;
t61 = t212 * t100;
t60 = t212 * t99;
t56 = 0.2e1 * t89 * t62;
t55 = -t80 * qJD(4) + t198;
t53 = t89 * t118 - t62 * t132;
t46 = -t68 * qJD(5) + (-t210 * t189 + t194) * qJD(4);
t45 = t100 * qJD(4) * t210 + t107 * t176 + t108 * t177;
t43 = -t177 * t192 + (t212 * t191 + t157) * t130 + t213 * t126;
t42 = -t99 * t182 - t212 * t87;
t36 = -t126 * t66 + t130 * t59;
t27 = t64 * pkin(4) + t51;
t16 = t31 * pkin(4) + t24;
t15 = -qJD(5) * t37 + t153;
t13 = qJD(5) * t40 + t126 * t32 + t130 * t31;
t12 = -qJD(5) * t39 - t126 * t31 + t130 * t32;
t8 = -t126 * t19 + t130 * t18;
t2 = -t9 * qJD(5) + t166;
t3 = [0, 0, 0, 0.2e1 * t148, 0.2e1 * (-t129 ^ 2 + t133 ^ 2) * t119 * qJD(2), 0.2e1 * t125 * t159, -0.2e1 * t125 * t160, 0, -0.2e1 * t119 * pkin(1) * t185 - 0.2e1 * t86 * t125, -0.2e1 * pkin(1) * t164 + 0.2e1 * t85 * t125, 0.2e1 * t90 * t63, -0.2e1 * t90 * t62 - 0.2e1 * t63 * t89, (-t133 * t63 + t90 * t185) * t211, (t133 * t62 - t89 * t185) * t211, -0.2e1 * t148, 0.2e1 * t81 * t62 + 0.2e1 * t86 * t89 + 0.2e1 * (-t26 * t133 + t152 * t185) * t124, 0.2e1 * t81 * t63 + 0.2e1 * t86 * t90 + 0.2e1 * (-t25 * t133 - t204 * t185) * t124, 0.2e1 * t65 * t32, -0.2e1 * t65 * t31 - 0.2e1 * t32 * t64, 0.2e1 * t32 * t89 + 0.2e1 * t65 * t62, -0.2e1 * t31 * t89 - 0.2e1 * t64 * t62, t56, 0.2e1 * t11 * t89 + 0.2e1 * t20 * t62 + 0.2e1 * t24 * t64 + 0.2e1 * t51 * t31, 0.2e1 * t10 * t89 - 0.2e1 * t21 * t62 + 0.2e1 * t24 * t65 + 0.2e1 * t51 * t32, 0.2e1 * t40 * t12, -0.2e1 * t12 * t39 - 0.2e1 * t40 * t13, 0.2e1 * t12 * t89 + 0.2e1 * t40 * t62, -0.2e1 * t13 * t89 - 0.2e1 * t39 * t62, t56, 0.2e1 * t27 * t13 + 0.2e1 * t16 * t39 + 0.2e1 * t2 * t89 + 0.2e1 * t8 * t62, 0.2e1 * t1 * t89 + 0.2e1 * t27 * t12 + 0.2e1 * t16 * t40 - 0.2e1 * t9 * t62; 0, 0, 0, 0, 0, t159, -t160, 0, -t86, t85, t63 * t128 + t90 * t182, -t128 * t62 + t63 * t132 + (-t128 * t90 - t132 * t89) * qJD(3), t137 * t124, t138 * t124, 0, -pkin(2) * t62 + t81 * t118 - t86 * t132 - t137 * t206, -pkin(2) * t63 + t86 * t128 - t138 * t206 + t81 * t182, t65 * t154 + (-t65 * t180 + t199) * t128, t144 * t182 + (-t200 - t131 * t31 + (t127 * t64 - t131 * t65) * qJD(4)) * t128, (t89 * t183 - t32) * t132 + (qJD(3) * t65 - t140) * t128, (-t127 * t197 + t31) * t132 + (-qJD(3) * t64 - t141) * t128, t53, t55 * t89 + t79 * t62 + (-t11 + (pkin(8) * t64 + t127 * t51) * qJD(3)) * t132 + (pkin(8) * t31 + qJD(3) * t20 + t143) * t128, t54 * t89 - t80 * t62 + (-t10 + (pkin(8) * t65 + t131 * t51) * qJD(3)) * t132 + (pkin(8) * t32 - qJD(3) * t21 - t142) * t128, -t12 * t88 + t40 * t42, -t12 * t87 + t88 * t13 - t42 * t39 - t40 * t43, t40 * t118 - t12 * t132 + t42 * t89 - t88 * t62, -t39 * t118 + t13 * t132 - t43 * t89 - t87 * t62, t53, t101 * t13 + t8 * t118 - t2 * t132 + t15 * t89 + t16 * t87 + t27 * t43 + t36 * t62 + t75 * t39, -t1 * t132 + t101 * t12 - t9 * t118 + t14 * t89 - t16 * t88 + t27 * t42 - t37 * t62 + t75 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, -0.2e1 * t150, 0, 0, 0, t128 * t175, t132 * t175, -0.2e1 * t121 * t156 + 0.2e1 * t122 * t155, 0.2e1 * t121 * t151 + t147 * t214, 0.2e1 * t128 * t162 + 0.2e1 * t186 * t183, -0.2e1 * t127 * t150 + 0.2e1 * t128 * t161, t111, 0.2e1 * t79 * t118 - 0.2e1 * t55 * t132 + 0.2e1 * (t121 * t179 + t127 * t149) * pkin(8), -0.2e1 * t80 * t118 - 0.2e1 * t54 * t132 + 0.2e1 * (-t121 * t180 + t131 * t149) * pkin(8), -0.2e1 * t88 * t42, -0.2e1 * t42 * t87 + 0.2e1 * t88 * t43, -0.2e1 * t88 * t118 - 0.2e1 * t42 * t132, -0.2e1 * t87 * t118 + 0.2e1 * t43 * t132, t111, 0.2e1 * t101 * t43 + 0.2e1 * t36 * t118 - 0.2e1 * t15 * t132 + 0.2e1 * t75 * t87, 0.2e1 * t101 * t42 - 0.2e1 * t37 * t118 - 0.2e1 * t14 * t132 - 0.2e1 * t75 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t62, t160, t26, t25, t65 * t179 + t200, qJD(4) * t144 - t127 * t31 + t199, t141, -t140, 0, -pkin(3) * t31 - pkin(9) * t141 + t142, -pkin(3) * t32 + pkin(9) * t140 + t143, t12 * t100 - t40 * t60, -t100 * t13 - t12 * t99 + t60 * t39 - t40 * t61, t100 * t62 - t60 * t89, -t61 * t89 - t99 * t62, 0, t117 * t13 + t16 * t99 + t171 * t39 + t27 * t61 + t46 * t89 + t67 * t62, t16 * t100 + t117 * t12 + t171 * t40 - t27 * t60 + t45 * t89 - t68 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, -t118, 0, -t170, pkin(8) * t118, -t128 * t151 + t147, t156 * t214 - t187 * t182, t158 - t161, t136, 0, (pkin(9) * t188 + (-pkin(3) * t131 + t205) * t128) * qJD(4) + (t127 * t146 - t114) * qJD(3), (pkin(8) * t191 + t127 * t145) * qJD(4) + (t131 * t146 + t173) * qJD(3), t42 * t100 + t88 * t60, -t100 * t43 - t42 * t99 + t60 * t87 + t88 * t61, t100 * t118 + t60 * t132, -t99 * t118 + t61 * t132, 0, t101 * t61 + t117 * t43 + t67 * t118 - t46 * t132 + t171 * t87 + t75 * t99, t75 * t100 - t101 * t60 + t117 * t42 - t68 * t118 - t45 * t132 - t171 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t156, -0.2e1 * t151, 0, 0, 0, t127 * t174, t131 * t174, -0.2e1 * t100 * t60, -0.2e1 * t100 * t61 + 0.2e1 * t60 * t99, 0, 0, 0, 0.2e1 * t117 * t61 + 0.2e1 * t171 * t99, 0.2e1 * t100 * t171 - 0.2e1 * t117 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, t62, t11, t10, 0, 0, t12, -t13, t62, t62 * t207 + (-t202 + (-t18 - t209) * t126) * qJD(5) + t166, (-t126 * t62 - t176 * t89) * pkin(4) + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, -t135, t118, t55, t54, 0, 0, t42, -t43, t118, t118 * t207 + (-t201 + (t132 * pkin(4) - t59) * t126) * qJD(5) + t153, (-t126 * t118 + t132 * t176) * pkin(4) + t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t180, 0, -pkin(9) * t179, pkin(9) * t180, 0, 0, -t60, -t61, 0, t46, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t169, -0.2e1 * t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t13, t62, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t43, t118, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t61, 0, t46, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, -t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
