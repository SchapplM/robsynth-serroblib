% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:51:50
% EndTime: 2019-03-09 20:51:57
% DurationCPUTime: 2.87s
% Computational Cost: add. (3780->311), mult. (8350->479), div. (0->0), fcn. (7390->6), ass. (0->170)
t209 = qJD(2) + qJD(3);
t123 = sin(qJ(3));
t124 = sin(qJ(2));
t126 = cos(qJ(3));
t127 = cos(qJ(2));
t93 = t123 * t127 + t126 * t124;
t52 = t209 * t93;
t92 = t123 * t124 - t126 * t127;
t211 = t52 * qJ(5) + t92 * qJD(5);
t125 = cos(qJ(4));
t115 = qJD(4) * t125;
t122 = sin(qJ(4));
t213 = qJ(6) * t115 + t122 * qJD(6);
t212 = 0.2e1 * t211;
t116 = t122 * qJ(5);
t210 = t125 * pkin(4) + t116;
t120 = t122 ^ 2;
t121 = t125 ^ 2;
t183 = t120 - t121;
t155 = t183 * qJD(4);
t177 = t125 * qJD(6);
t180 = qJD(4) * t122;
t166 = t93 * t180;
t51 = t209 * t92;
t32 = t125 * t51 + t166;
t179 = t124 * qJD(2);
t171 = pkin(2) * t179;
t23 = t52 * pkin(3) + t51 * pkin(9) + t171;
t205 = -pkin(8) - pkin(7);
t100 = t205 * t124;
t101 = t205 * t127;
t161 = qJD(2) * t205;
t152 = t127 * t161;
t153 = t124 * t161;
t181 = qJD(3) * t126;
t182 = qJD(3) * t123;
t29 = -t100 * t181 - t101 * t182 - t123 * t152 - t126 * t153;
t111 = -t127 * pkin(2) - pkin(1);
t39 = t92 * pkin(3) - t93 * pkin(9) + t111;
t59 = t123 * t100 - t126 * t101;
t10 = -t59 * t115 + t122 * t29 + t125 * t23 - t39 * t180;
t43 = t52 * pkin(4);
t7 = -t43 - t10;
t208 = -t32 * qJ(6) + t93 * t177 - t7;
t128 = -pkin(4) - pkin(5);
t164 = qJD(4) * t93 * qJ(5);
t207 = -t128 * t51 - t164;
t30 = t59 * qJD(3) + t123 * t153 - t126 * t152;
t192 = t51 * qJ(5);
t206 = -t192 + (qJD(4) * t128 + qJD(5)) * t93;
t129 = 0.2e1 * qJD(5);
t204 = pkin(9) * t52;
t203 = pkin(9) * t92;
t186 = t125 * qJ(5);
t58 = -t126 * t100 - t123 * t101;
t19 = (t128 * t122 + t186) * t93 - t58;
t18 = t19 * t115;
t8 = t207 * t122 + t206 * t125 - t30;
t202 = t8 * t122 + t18;
t201 = t93 * t51;
t200 = pkin(9) - qJ(6);
t47 = t58 * t115;
t199 = t30 * t122 + t47;
t170 = pkin(2) * t182;
t75 = pkin(4) * t180 - qJ(5) * t115 - t122 * qJD(5);
t61 = -pkin(5) * t180 - t75;
t55 = t61 - t170;
t117 = t125 * pkin(5);
t110 = -t126 * pkin(2) - pkin(3);
t83 = t110 - t210;
t68 = t117 - t83;
t60 = t68 * t115;
t198 = t55 * t122 + t60;
t197 = t122 * t39 + t125 * t59;
t174 = pkin(3) + t210;
t84 = t117 + t174;
t66 = t84 * t115;
t196 = t61 * t122 + t66;
t62 = t75 + t170;
t195 = -t62 - t75;
t109 = t123 * pkin(2) + pkin(9);
t194 = t109 * t52;
t193 = t109 * t92;
t191 = t110 * t115 + t122 * t170;
t190 = qJD(4) * t19;
t145 = t122 * pkin(4) - t186;
t28 = t145 * t93 + t58;
t189 = qJD(4) * t28;
t188 = t122 * qJ(6);
t187 = t122 * t125;
t185 = -qJ(6) + t109;
t178 = t125 * qJD(5);
t176 = t127 * qJD(2);
t175 = -0.2e1 * pkin(1) * qJD(2);
t16 = t92 * qJ(5) + t197;
t173 = pkin(3) * t180;
t172 = pkin(3) * t115;
t169 = pkin(2) * t181;
t168 = pkin(9) * t180;
t167 = pkin(9) * t115;
t165 = t93 * t115;
t24 = t28 * t180;
t46 = t58 * t180;
t163 = t68 * t180;
t162 = t84 * t180;
t160 = t122 * t115;
t53 = t122 * t59;
t159 = t125 * t39 - t53;
t158 = -0.4e1 * t93 * t187;
t154 = qJ(6) * t180 - t177;
t63 = t109 * t180 - t125 * t169;
t49 = t154 - t63;
t86 = t185 * t122;
t157 = qJD(4) * t86 + t49;
t74 = t154 - t168;
t98 = t200 * t122;
t156 = qJD(4) * t98 + t74;
t151 = -t19 * t51 + t8 * t93;
t133 = t192 + (pkin(4) * qJD(4) - qJD(5)) * t93;
t140 = -t51 * pkin(4) + t164;
t11 = t140 * t122 + t133 * t125 + t30;
t150 = t11 * t93 - t28 * t51;
t149 = t30 * t93 - t58 * t51;
t148 = -t51 * t92 + t93 * t52;
t147 = -t51 * t68 + t55 * t93;
t146 = -t51 * t84 + t61 * t93;
t144 = -t110 * t93 + t193;
t17 = -t92 * pkin(4) - t159;
t143 = t122 * t17 + t125 * t16;
t142 = -t122 * t16 + t125 * t17;
t141 = t110 * t180 - t125 * t170;
t72 = (t120 + t121) * t169;
t139 = -t125 * t52 + t92 * t180;
t9 = -t39 * t115 - t122 * t23 + t125 * t29 + t59 * t180;
t137 = t174 * t51 + t75 * t93 - t204;
t136 = -t11 + (-t174 * t93 - t203) * qJD(4);
t135 = -t11 + (t83 * t93 - t193) * qJD(4);
t64 = t109 * t115 + t122 * t169;
t132 = -t51 * t188 + t213 * t93 - t9;
t73 = -t210 * qJD(4) + t178;
t131 = -t92 * t169 - t51 * t83 + t62 * t93 - t194;
t4 = -t9 + t211;
t1 = t142 * qJD(4) + t7 * t122 + t4 * t125;
t130 = -t194 - t110 * t51 + (t123 * t93 - t126 * t92) * qJD(3) * pkin(2);
t119 = qJ(5) * t129;
t104 = 0.2e1 * t160;
t99 = t200 * t125;
t91 = -0.2e1 * t155;
t90 = t93 ^ 2;
t88 = t99 * t180;
t87 = t185 * t125;
t81 = t174 * t180;
t76 = t167 - t213;
t69 = -t178 + (-t125 * t128 + t116) * qJD(4);
t67 = t87 * t180;
t65 = t83 * t180;
t57 = t61 * t125;
t50 = t64 - t213;
t45 = t55 * t125;
t34 = t122 * t51 - t165;
t33 = t92 * t115 + t122 * t52;
t22 = -t93 * t155 - t51 * t187;
t15 = qJD(4) * t158 + t183 * t51;
t14 = t93 * t188 + t16;
t13 = t14 * t180;
t12 = t53 + t128 * t92 + (-qJ(6) * t93 - t39) * t125;
t6 = t8 * t125;
t3 = t132 + t211;
t2 = -t52 * pkin(5) - t208;
t5 = [0, 0, 0, 0.2e1 * t124 * t176, 0.2e1 * (-t124 ^ 2 + t127 ^ 2) * qJD(2), 0, 0, 0, t124 * t175, t127 * t175, -0.2e1 * t201, -0.2e1 * t148, 0, 0, 0, 0.2e1 * t111 * t52 + 0.2e1 * t92 * t171, -0.2e1 * t111 * t51 + 0.2e1 * t93 * t171, -0.2e1 * t121 * t201 - 0.2e1 * t90 * t160, 0.2e1 * t90 * t155 - t158 * t51, 0.2e1 * t125 * t148 - 0.2e1 * t166 * t92, -0.2e1 * t148 * t122 - 0.2e1 * t165 * t92, 0.2e1 * t92 * t52, 0.2e1 * t10 * t92 + 0.2e1 * t149 * t122 + 0.2e1 * t159 * t52 + 0.2e1 * t47 * t93, 0.2e1 * t149 * t125 - 0.2e1 * t197 * t52 - 0.2e1 * t93 * t46 + 0.2e1 * t9 * t92, 0.2e1 * t122 * t150 + 0.2e1 * t165 * t28 - 0.2e1 * t17 * t52 - 0.2e1 * t7 * t92, -0.2e1 * t142 * t51 + 0.2e1 * (-qJD(4) * t143 - t122 * t4 + t125 * t7) * t93, -0.2e1 * t125 * t150 + 0.2e1 * t16 * t52 + 0.2e1 * t24 * t93 + 0.2e1 * t4 * t92, 0.2e1 * t28 * t11 + 0.2e1 * t16 * t4 + 0.2e1 * t17 * t7, -0.2e1 * t12 * t52 - 0.2e1 * t122 * t151 - 0.2e1 * t18 * t93 - 0.2e1 * t2 * t92, 0.2e1 * t125 * t151 + 0.2e1 * t14 * t52 - 0.2e1 * t166 * t19 + 0.2e1 * t3 * t92, -0.2e1 * (-t12 * t125 + t122 * t14) * t51 + 0.2e1 * (t122 * t3 - t125 * t2 + (t12 * t122 + t125 * t14) * qJD(4)) * t93, 0.2e1 * t12 * t2 + 0.2e1 * t14 * t3 + 0.2e1 * t19 * t8; 0, 0, 0, 0, 0, t176, -t179, 0, -pkin(7) * t176, pkin(7) * t179, 0, 0, -t51, -t52, 0, -t30, t29, t22, t15, t33, -t139, 0, t46 + (-qJD(4) * t144 - t30) * t125 + t130 * t122, t125 * t130 + t144 * t180 + t199, t122 * t131 + t125 * t135 + t24, t1, t135 * t122 + (-t131 - t189) * t125, t1 * t109 + t11 * t83 + t143 * t169 + t28 * t62, -t93 * t60 - t50 * t92 - t86 * t52 + t6 + (-t147 - t190) * t122, t125 * t147 - t163 * t93 + t49 * t92 + t87 * t52 + t202, t13 + (t157 * t93 - t51 * t87 - t2) * t122 + (-t50 * t93 + t51 * t86 - t3 + (t87 * t93 - t12) * qJD(4)) * t125, t12 * t50 + t14 * t49 + t19 * t55 + t2 * t86 + t3 * t87 + t8 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t170, -0.2e1 * t169, t104, t91, 0, 0, 0, 0.2e1 * t141, 0.2e1 * t191, -0.2e1 * t62 * t125 + 0.2e1 * t65, 0.2e1 * t72, -0.2e1 * t115 * t83 - 0.2e1 * t62 * t122, 0.2e1 * t109 * t72 + 0.2e1 * t83 * t62, 0.2e1 * t45 - 0.2e1 * t163, 0.2e1 * t198, -0.2e1 * t50 * t122 - 0.2e1 * t125 * t157 + 0.2e1 * t67, 0.2e1 * t87 * t49 + 0.2e1 * t86 * t50 + 0.2e1 * t68 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t52, 0, -t30, t29, t22, t15, t33, -t139, 0, t46 + (pkin(3) * t51 - t204) * t122 + (-t30 + (-pkin(3) * t93 - t203) * qJD(4)) * t125, pkin(3) * t32 + pkin(9) * t139 + t199, t122 * t137 + t125 * t136 + t24, t1, t136 * t122 + (-t137 - t189) * t125, pkin(9) * t1 - t11 * t174 + t28 * t75, -t93 * t66 - t98 * t52 - t76 * t92 + t6 + (-t146 - t190) * t122, t125 * t146 - t162 * t93 + t99 * t52 + t74 * t92 + t202, t13 + (t156 * t93 - t51 * t99 - t2) * t122 + (t51 * t98 - t76 * t93 - t3 + (t93 * t99 - t12) * qJD(4)) * t125, t12 * t76 + t14 * t74 + t19 * t61 + t2 * t98 + t3 * t99 + t8 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, -t169, t104, t91, 0, 0, 0, t141 - t173, -t172 + t191, t195 * t125 + t65 - t81, t72, t195 * t122 + (-t83 + t174) * t115, pkin(9) * t72 - t174 * t62 + t83 * t75, t45 + t57 + (-t68 - t84) * t180, t196 + t198, t67 + t88 + (-t50 - t76) * t122 + (-t49 - t74 + (-t86 - t98) * qJD(4)) * t125, t49 * t99 + t50 * t98 + t55 * t84 + t68 * t61 + t87 * t74 + t86 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t91, 0, 0, 0, -0.2e1 * t173, -0.2e1 * t172, -0.2e1 * t75 * t125 - 0.2e1 * t81, 0, 0.2e1 * t115 * t174 - 0.2e1 * t75 * t122, -0.2e1 * t174 * t75, 0.2e1 * t57 - 0.2e1 * t162, 0.2e1 * t196, -0.2e1 * t76 * t122 - 0.2e1 * t125 * t156 + 0.2e1 * t88, 0.2e1 * t84 * t61 + 0.2e1 * t99 * t74 + 0.2e1 * t98 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t34, t52, t10, t9, -t7 + t43, t122 * t133 - t125 * t140, -t9 + t212, -t7 * pkin(4) + t4 * qJ(5) + t16 * qJD(5) (pkin(5) - t128) * t52 + t208, t132 + t212, t206 * t122 - t207 * t125, t3 * qJ(5) + t14 * qJD(5) + t2 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t180, 0, -t64, t63, -t64, t73, -t63, t109 * t73 - t145 * t169, -t50, t49, t69, t49 * qJ(5) + t87 * qJD(5) + t50 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t180, 0, -t167, t168, -t167, t73, -t168, t73 * pkin(9), -t76, t74, t69, t74 * qJ(5) + t99 * qJD(5) + t76 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t119, 0, t129, 0, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t32, 0, t7, -t52, 0, t32, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, 0, t64, 0, 0, -t115, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, 0, t167, 0, 0, -t115, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t115, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t115, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
