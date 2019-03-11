% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:23:12
% EndTime: 2019-03-09 22:23:20
% DurationCPUTime: 3.15s
% Computational Cost: add. (6937->292), mult. (17035->533), div. (0->0), fcn. (15846->10), ass. (0->169)
t152 = sin(qJ(2));
t151 = sin(qJ(3));
t156 = cos(qJ(2));
t194 = t156 * qJD(2);
t183 = t151 * t194;
t155 = cos(qJ(3));
t198 = qJD(3) * t155;
t224 = t152 * t198 + t183;
t223 = -0.4e1 * t152;
t150 = sin(qJ(4));
t154 = cos(qJ(4));
t115 = t150 * t155 + t154 * t151;
t101 = t115 * t152;
t219 = pkin(8) + pkin(9);
t127 = t219 * t151;
t128 = t219 * t155;
t203 = -t150 * t127 + t154 * t128;
t169 = -t156 * pkin(2) - t152 * pkin(8);
t123 = -pkin(1) + t169;
t206 = t155 * t156;
t134 = pkin(7) * t206;
t202 = t151 * t123 + t134;
t180 = t155 * t194;
t199 = qJD(3) * t151;
t222 = -t152 * t199 + t180;
t144 = t152 ^ 2;
t172 = (-t156 ^ 2 + t144) * qJD(2);
t145 = t155 ^ 2;
t201 = t151 ^ 2 - t145;
t173 = t201 * qJD(3);
t221 = qJD(3) + qJD(4);
t168 = pkin(2) * t152 - pkin(8) * t156;
t121 = t168 * qJD(2);
t141 = t152 * qJD(2);
t197 = qJD(3) * t156;
t187 = t151 * t197;
t159 = t155 * t141 + t187;
t70 = pkin(7) * t159 - t151 * t121 - t123 * t198;
t149 = sin(qJ(6));
t153 = cos(qJ(6));
t147 = sin(pkin(11));
t148 = cos(pkin(11));
t114 = t150 * t151 - t154 * t155;
t102 = t114 * t152;
t113 = t155 * t123;
t207 = t152 * t155;
t217 = pkin(7) * t151;
t87 = -pkin(9) * t207 + t113 + (-pkin(3) - t217) * t156;
t208 = t151 * t152;
t92 = -pkin(9) * t208 + t202;
t176 = -t150 * t92 + t154 * t87;
t46 = -t156 * pkin(4) + t102 * qJ(5) + t176;
t89 = t154 * t92;
t215 = t150 * t87 + t89;
t52 = -t101 * qJ(5) + t215;
t29 = -t147 * t52 + t148 * t46;
t74 = -t147 * t101 - t148 * t102;
t20 = -t156 * pkin(5) - t74 * pkin(10) + t29;
t30 = t147 * t46 + t148 * t52;
t73 = -t148 * t101 + t147 * t102;
t21 = t73 * pkin(10) + t30;
t213 = t149 * t21;
t65 = -t221 * t101 - t114 * t194;
t196 = qJD(4) * t150;
t66 = -t196 * t208 + (t221 * t207 + t183) * t154 + t222 * t150;
t33 = -t147 * t66 + t148 * t65;
t184 = t151 * t141;
t204 = pkin(7) * t184 + t155 * t121;
t60 = (pkin(3) * t152 - pkin(9) * t206) * qJD(2) + (-t134 + (pkin(9) * t152 - t123) * t151) * qJD(3) + t204;
t64 = -pkin(9) * t224 - t70;
t177 = -t150 * t64 + t154 * t60;
t28 = -t215 * qJD(4) + t177;
t14 = pkin(4) * t141 - t65 * qJ(5) + t102 * qJD(5) + t28;
t195 = qJD(4) * t154;
t27 = -t150 * t60 - t154 * t64 - t87 * t195 + t92 * t196;
t16 = -t66 * qJ(5) - t101 * qJD(5) - t27;
t8 = t148 * t14 - t147 * t16;
t6 = pkin(5) * t141 - t33 * pkin(10) + t8;
t32 = -t147 * t65 - t148 * t66;
t9 = t147 * t14 + t148 * t16;
t7 = t32 * pkin(10) + t9;
t220 = (-qJD(6) * t20 - t7) * t153 + qJD(6) * t213 - t149 * t6;
t218 = pkin(4) * t147;
t216 = t154 * pkin(3);
t189 = qJD(3) * t219;
t119 = t151 * t189;
t120 = t155 * t189;
t67 = t154 * t119 + t150 * t120 + t127 * t195 + t128 * t196;
t91 = t221 * t115;
t40 = -t91 * qJ(5) - t114 * qJD(5) - t67;
t68 = -t203 * qJD(4) + t150 * t119 - t154 * t120;
t90 = t221 * t114;
t41 = t90 * qJ(5) - t115 * qJD(5) + t68;
t25 = t147 * t41 + t148 * t40;
t174 = -t154 * t127 - t150 * t128;
t78 = -t115 * qJ(5) + t174;
t79 = -t114 * qJ(5) + t203;
t49 = t147 * t78 + t148 * t79;
t214 = pkin(3) * qJD(4);
t135 = t148 * pkin(4) + pkin(5);
t137 = pkin(4) + t216;
t210 = t147 * t150;
t103 = -pkin(3) * t210 + t148 * t137;
t98 = pkin(5) + t103;
t212 = -t135 - t98;
t209 = t148 * t150;
t122 = pkin(3) * t208 + t152 * pkin(7);
t193 = -0.2e1 * pkin(1) * qJD(2);
t192 = -0.2e1 * pkin(2) * qJD(3);
t140 = pkin(3) * t199;
t191 = pkin(3) * t196;
t190 = pkin(3) * t195;
t139 = pkin(7) * t194;
t93 = t224 * pkin(3) + t139;
t138 = -t155 * pkin(3) - pkin(2);
t185 = t155 * t197;
t182 = t151 * t198;
t181 = t152 * t194;
t80 = t91 * pkin(4) + t140;
t104 = pkin(3) * t209 + t147 * t137;
t178 = t104 + t218;
t24 = -t147 * t40 + t148 * t41;
t48 = -t147 * t79 + t148 * t78;
t100 = (t148 * t154 - t210) * t214;
t99 = (-t147 * t154 - t209) * t214;
t175 = -t149 * t100 + t153 * t99;
t171 = 0.2e1 * t181;
t170 = t151 * t180;
t88 = t101 * pkin(4) + t122;
t167 = t149 * t20 + t153 * t21;
t86 = -t147 * t114 + t148 * t115;
t34 = -t86 * pkin(10) + t48;
t85 = -t148 * t114 - t147 * t115;
t35 = t85 * pkin(10) + t49;
t166 = t149 * t35 - t153 * t34;
t165 = t149 * t34 + t153 * t35;
t42 = t149 * t74 - t153 * t73;
t43 = t149 * t73 + t153 * t74;
t55 = t149 * t86 - t153 * t85;
t56 = t149 * t85 + t153 * t86;
t164 = -t153 * t100 - t149 * t99;
t163 = t153 * t104 + t149 * t98;
t162 = t149 * t104 - t153 * t98;
t53 = t66 * pkin(4) + t93;
t97 = t114 * pkin(4) + t138;
t161 = t149 * t135 + t153 * t218;
t160 = -t153 * t135 + t149 * t218;
t2 = -qJD(6) * t167 - t149 * t7 + t153 * t6;
t132 = -0.2e1 * t181;
t96 = t161 * qJD(6);
t95 = t160 * qJD(6);
t71 = -t202 * qJD(3) + t204;
t69 = -t85 * pkin(5) + t97;
t62 = -t147 * t91 - t148 * t90;
t61 = t147 * t90 - t148 * t91;
t54 = -t73 * pkin(5) + t88;
t51 = -qJD(6) * t163 + t175;
t50 = qJD(6) * t162 + t164;
t39 = -t61 * pkin(5) + t80;
t26 = -t32 * pkin(5) + t53;
t23 = qJD(6) * t56 + t149 * t62 - t153 * t61;
t22 = -qJD(6) * t55 + t149 * t61 + t153 * t62;
t18 = t61 * pkin(10) + t25;
t17 = -t62 * pkin(10) + t24;
t11 = qJD(6) * t43 + t149 * t33 - t153 * t32;
t10 = -qJD(6) * t42 + t149 * t32 + t153 * t33;
t4 = -qJD(6) * t165 - t149 * t18 + t153 * t17;
t3 = qJD(6) * t166 - t149 * t17 - t153 * t18;
t1 = [0, 0, 0, t171, -0.2e1 * t172, 0, 0, 0, t152 * t193, t156 * t193, -0.2e1 * t144 * t182 + 0.2e1 * t145 * t181, 0.2e1 * t144 * t173 + t170 * t223, 0.2e1 * t152 * t187 + 0.2e1 * t155 * t172, -0.2e1 * t151 * t172 + 0.2e1 * t152 * t185, t132, 0.2e1 * t113 * t141 - 0.2e1 * t71 * t156 + 0.2e1 * (t144 * t198 + t151 * t181) * pkin(7), -0.2e1 * t70 * t156 - 0.2e1 * t202 * t141 + 0.2e1 * (-t144 * t199 + t155 * t171) * pkin(7), -0.2e1 * t102 * t65, -0.2e1 * t65 * t101 + 0.2e1 * t102 * t66, -0.2e1 * t102 * t141 - 0.2e1 * t65 * t156, -0.2e1 * t101 * t141 + 0.2e1 * t66 * t156, t132, 0.2e1 * t93 * t101 + 0.2e1 * t122 * t66 + 0.2e1 * t176 * t141 - 0.2e1 * t28 * t156, -0.2e1 * t93 * t102 + 0.2e1 * t122 * t65 - 0.2e1 * t215 * t141 - 0.2e1 * t27 * t156, -0.2e1 * t29 * t33 + 0.2e1 * t30 * t32 + 0.2e1 * t9 * t73 - 0.2e1 * t8 * t74, 0.2e1 * t29 * t8 + 0.2e1 * t30 * t9 + 0.2e1 * t88 * t53, 0.2e1 * t43 * t10, -0.2e1 * t10 * t42 - 0.2e1 * t43 * t11, -0.2e1 * t10 * t156 + 0.2e1 * t43 * t141, 0.2e1 * t11 * t156 - 0.2e1 * t42 * t141, t132, -0.2e1 * t2 * t156 + 0.2e1 * (t153 * t20 - t213) * t141 + 0.2e1 * t26 * t42 + 0.2e1 * t54 * t11, 0.2e1 * t54 * t10 - 0.2e1 * t167 * t141 - 0.2e1 * t156 * t220 + 0.2e1 * t26 * t43; 0, 0, 0, 0, 0, t194, -t141, 0, -t139, pkin(7) * t141, -t152 * t173 + t170, t182 * t223 - t201 * t194, t184 - t185, t159, 0 (pkin(8) * t206 + (-pkin(2) * t155 + t217) * t152) * qJD(3) + (t151 * t169 - t134) * qJD(2) (pkin(7) * t207 + t151 * t168) * qJD(3) + (t155 * t169 + t156 * t217) * qJD(2), t102 * t90 + t65 * t115, t90 * t101 + t102 * t91 - t65 * t114 - t115 * t66, t115 * t141 + t90 * t156, -t114 * t141 + t91 * t156, 0, t101 * t140 + t93 * t114 + t122 * t91 + t138 * t66 + t174 * t141 - t68 * t156, -t102 * t140 + t93 * t115 - t122 * t90 + t138 * t65 - t203 * t141 - t67 * t156, -t24 * t74 + t25 * t73 - t29 * t62 + t30 * t61 + t49 * t32 - t48 * t33 - t8 * t86 + t9 * t85, t29 * t24 + t30 * t25 + t8 * t48 + t9 * t49 + t53 * t97 + t88 * t80, t10 * t56 + t43 * t22, -t10 * t55 - t56 * t11 - t22 * t42 - t43 * t23, t56 * t141 - t22 * t156, -t55 * t141 + t23 * t156, 0, t69 * t11 - t166 * t141 - t4 * t156 + t54 * t23 + t26 * t55 + t39 * t42, t69 * t10 - t165 * t141 - t3 * t156 + t54 * t22 + t26 * t56 + t39 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t182, -0.2e1 * t173, 0, 0, 0, t151 * t192, t155 * t192, -0.2e1 * t115 * t90, 0.2e1 * t90 * t114 - 0.2e1 * t115 * t91, 0, 0, 0, 0.2e1 * t114 * t140 + 0.2e1 * t138 * t91, 0.2e1 * t115 * t140 - 0.2e1 * t138 * t90, -0.2e1 * t24 * t86 + 0.2e1 * t25 * t85 - 0.2e1 * t48 * t62 + 0.2e1 * t49 * t61, 0.2e1 * t48 * t24 + 0.2e1 * t49 * t25 + 0.2e1 * t97 * t80, 0.2e1 * t56 * t22, -0.2e1 * t22 * t55 - 0.2e1 * t56 * t23, 0, 0, 0, 0.2e1 * t69 * t23 + 0.2e1 * t39 * t55, 0.2e1 * t69 * t22 + 0.2e1 * t39 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, -t224, t141, t71, t70, 0, 0, t65, -t66, t141, t141 * t216 + (-t89 + (t156 * pkin(3) - t87) * t150) * qJD(4) + t177 (-t150 * t141 + t156 * t195) * pkin(3) + t27, t100 * t73 - t103 * t33 + t104 * t32 - t99 * t74, t30 * t100 + t8 * t103 + t9 * t104 + t29 * t99, 0, 0, t10, -t11, t141, -t162 * t141 - t51 * t156 + t2, -t163 * t141 - t50 * t156 + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, -t199, 0, -pkin(8) * t198, pkin(8) * t199, 0, 0, -t90, -t91, 0, t68, t67, t100 * t85 - t103 * t62 + t104 * t61 - t99 * t86, t49 * t100 + t24 * t103 + t25 * t104 + t48 * t99, 0, 0, t22, -t23, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t191, -0.2e1 * t190, 0, 0.2e1 * t104 * t100 + 0.2e1 * t103 * t99, 0, 0, 0, 0, 0, 0.2e1 * t51, 0.2e1 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t66, t141, t28, t27 (t147 * t32 - t148 * t33) * pkin(4) (t147 * t9 + t148 * t8) * pkin(4), 0, 0, t10, -t11, t141, -t160 * t141 + t96 * t156 + t2, -t161 * t141 - t95 * t156 + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t91, 0, t68, t67 (t147 * t61 - t148 * t62) * pkin(4) (t147 * t25 + t148 * t24) * pkin(4), 0, 0, t22, -t23, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, -t190, 0 (t100 * t147 + t148 * t99) * pkin(4), 0, 0, 0, 0, 0 (t212 * t149 - t178 * t153) * qJD(6) + t175 (t178 * t149 + t212 * t153) * qJD(6) + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t96, 0.2e1 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, 0, 0, 0, 0, t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t141, t2, t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t23, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
