% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:09:57
% EndTime: 2019-03-09 10:10:07
% DurationCPUTime: 3.83s
% Computational Cost: add. (9193->270), mult. (19845->481), div. (0->0), fcn. (20562->10), ass. (0->169)
t138 = sin(pkin(11));
t136 = t138 ^ 2;
t139 = cos(pkin(11));
t137 = t139 ^ 2;
t198 = t136 + t137;
t234 = t198 * qJD(5);
t235 = 0.2e1 * t234;
t141 = sin(qJ(4));
t207 = sin(pkin(10));
t187 = t207 * pkin(2);
t127 = t141 * t187;
t208 = cos(pkin(10));
t130 = t208 * pkin(2) + pkin(3);
t230 = cos(qJ(4));
t184 = qJD(4) * t230;
t97 = qJD(4) * t127 - t130 * t184;
t95 = qJD(5) - t97;
t233 = t198 * t95;
t142 = sin(qJ(2));
t143 = cos(qJ(2));
t112 = t208 * t142 + t207 * t143;
t161 = t207 * t142 - t208 * t143;
t159 = t230 * t161;
t174 = qJD(2) * t207;
t175 = qJD(2) * t208;
t165 = -t142 * t174 + t143 * t175;
t197 = qJD(4) * t141;
t200 = t142 * t175 + t143 * t174;
t54 = qJD(4) * t159 + t112 * t197 + t141 * t200 - t230 * t165;
t79 = t230 * t112 - t141 * t161;
t55 = t79 * qJD(4) + t141 * t165 + t230 * t200;
t78 = t112 * t141 + t159;
t171 = t54 * t78 - t55 * t79;
t232 = 0.2e1 * t171;
t106 = t141 * t130 + t230 * t187;
t140 = sin(qJ(6));
t229 = cos(qJ(6));
t188 = t229 * t139;
t231 = -t138 * t140 + t188;
t228 = t139 * pkin(5);
t223 = -qJ(3) - pkin(7);
t179 = qJD(2) * t223;
t162 = -t142 * qJD(3) + t143 * t179;
t163 = qJD(3) * t143 + t142 * t179;
t68 = t208 * t162 - t207 * t163;
t149 = -pkin(8) * t165 + t68;
t69 = t207 * t162 + t208 * t163;
t152 = -t200 * pkin(8) + t69;
t119 = t223 * t142;
t120 = t223 * t143;
t89 = t208 * t119 + t207 * t120;
t70 = -t112 * pkin(8) + t89;
t90 = t207 * t119 - t208 * t120;
t71 = -pkin(8) * t161 + t90;
t45 = t141 * t70 + t230 * t71;
t26 = qJD(4) * t45 + t141 * t152 - t230 * t149;
t44 = t141 * t71 - t230 * t70;
t227 = t44 * t26;
t98 = t106 * qJD(4);
t226 = t44 * t98;
t225 = t79 * t98;
t102 = qJ(5) + t106;
t224 = -pkin(9) - t102;
t222 = qJ(5) + pkin(9);
t189 = t229 * t138;
t201 = t140 * t139;
t113 = t189 + t201;
t109 = t113 * qJD(6);
t214 = t138 * t54;
t19 = -pkin(5) * t214 + t26;
t213 = t138 * t79;
t36 = pkin(5) * t213 + t44;
t221 = t36 * t109 - t19 * t231;
t183 = qJD(6) * t229;
t196 = qJD(6) * t140;
t186 = t138 * t196;
t108 = -t139 * t183 + t186;
t220 = -t36 * t108 + t19 * t113;
t23 = -t54 * t189 - t79 * t186 + (-t140 * t54 + t79 * t183) * t139;
t48 = t113 * t79;
t219 = t108 * t48 - t113 * t23;
t132 = -t143 * pkin(2) - pkin(1);
t96 = pkin(3) * t161 + t132;
t156 = t78 * pkin(4) + t96;
t154 = -t79 * qJ(5) + t156;
t28 = t138 * t154 + t139 * t45;
t105 = t230 * t130 - t127;
t103 = -pkin(4) - t105;
t94 = t103 - t228;
t218 = t94 * t109 - t231 * t98;
t217 = -t94 * t108 + t98 * t113;
t215 = t138 * t45;
t50 = t138 * t55;
t212 = t139 * t26;
t211 = t139 * t54;
t210 = t98 * t138;
t209 = t98 * t139;
t206 = t231 * t109;
t205 = t113 * t108;
t131 = -pkin(4) - t228;
t204 = t131 * t108;
t203 = t131 * t109;
t195 = t142 * qJD(2);
t194 = t143 * qJD(2);
t40 = 0.2e1 * t78 * t55;
t193 = -0.2e1 * t79 * t54;
t192 = -0.2e1 * pkin(1) * qJD(2);
t191 = pkin(2) * t195;
t46 = t138 * t211;
t153 = t78 * pkin(5) - t215 + (-t222 * t79 + t156) * t139;
t150 = t229 * t153;
t24 = -pkin(9) * t213 + t28;
t11 = -t140 * t24 + t150;
t151 = t140 * t153;
t12 = t229 * t24 + t151;
t148 = -t141 * t149 - t230 * t152 - t70 * t184 + t71 * t197;
t147 = t138 * t148;
t93 = t200 * pkin(3) + t191;
t158 = t55 * pkin(4) - t79 * qJD(5) + t93;
t155 = t222 * t54 + t158;
t144 = t55 * pkin(5) + t139 * t155 + t147;
t146 = t139 * t148;
t145 = t138 * t155 - t146;
t3 = -qJD(6) * t150 - t140 * t144 - t229 * t145 + t24 * t196;
t4 = -qJD(6) * t151 - t140 * t145 + t229 * t144 - t24 * t183;
t190 = t11 * t108 - t12 * t109 - t4 * t113 - t231 * t3;
t185 = t142 * t194;
t182 = t222 * t138;
t181 = t224 * t138;
t180 = t229 * qJD(5);
t166 = t229 * t181;
t135 = t139 * pkin(9);
t88 = t102 * t139 + t135;
t38 = -qJD(6) * t166 - t95 * t188 + (qJD(6) * t88 + t138 * t95) * t140;
t39 = -t88 * t183 - t95 * t201 + (-t224 * t196 - t229 * t95) * t138;
t58 = -t140 * t88 + t166;
t59 = t140 * t181 + t229 * t88;
t177 = t58 * t108 - t59 * t109 - t39 * t113 - t231 * t38;
t118 = t139 * qJ(5) + t135;
t167 = t229 * t182;
t63 = qJD(6) * t167 - t139 * t180 + (t138 * qJD(5) + qJD(6) * t118) * t140;
t64 = -t118 * t183 - qJD(5) * t201 + (t222 * t196 - t180) * t138;
t86 = -t140 * t118 - t167;
t87 = t229 * t118 - t140 * t182;
t176 = t86 * t108 - t87 * t109 - t64 * t113 - t231 * t63;
t172 = t26 * t79 - t44 * t54;
t170 = t108 * t78 - t113 * t55;
t22 = t79 * t109 + t231 * t54;
t49 = t231 * t79;
t169 = -t109 * t49 - t22 * t231;
t157 = t54 * qJ(5) + t158;
t13 = t139 * t157 + t147;
t14 = t138 * t157 - t146;
t5 = -t13 * t138 + t139 * t14;
t27 = t139 * t154 - t215;
t168 = -t138 * t27 + t139 * t28;
t164 = pkin(4) * t54 - qJ(5) * t55 - qJD(5) * t78;
t160 = -t102 * t55 - t103 * t54 - t78 * t95 + t225;
t82 = -0.2e1 * t205;
t81 = -0.2e1 * t206;
t56 = -0.2e1 * t108 * t231 - 0.2e1 * t113 * t109;
t51 = t139 * t55;
t33 = -t109 * t78 + t231 * t55;
t29 = (t136 - t137) * t54;
t25 = t26 * t138;
t16 = t109 * t48 - t23 * t231;
t15 = -t108 * t49 - t113 * t22;
t6 = t169 + t219;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t185, 0.2e1 * (-t142 ^ 2 + t143 ^ 2) * qJD(2), 0, -0.2e1 * t185, 0, 0, t142 * t192, t143 * t192, 0, 0, 0.2e1 * t112 * t165, -0.2e1 * t112 * t200 - 0.2e1 * t165 * t161, 0, 0.2e1 * t161 * t200, 0, 0, 0.2e1 * t132 * t200 + 0.2e1 * t161 * t191, 0.2e1 * t112 * t191 + 0.2e1 * t132 * t165, -0.2e1 * t68 * t112 - 0.2e1 * t69 * t161 - 0.2e1 * t89 * t165 - 0.2e1 * t90 * t200, 0.2e1 * t132 * t191 + 0.2e1 * t68 * t89 + 0.2e1 * t69 * t90, t193, t232, 0, t40, 0, 0, 0.2e1 * t55 * t96 + 0.2e1 * t78 * t93, -0.2e1 * t54 * t96 + 0.2e1 * t79 * t93, 0.2e1 * t148 * t78 - 0.2e1 * t45 * t55 + 0.2e1 * t172, -0.2e1 * t148 * t45 + 0.2e1 * t96 * t93 + 0.2e1 * t227, t137 * t193, 0.4e1 * t79 * t46, -0.2e1 * t171 * t139, t136 * t193, t138 * t232, t40, 0.2e1 * t13 * t78 + 0.2e1 * t138 * t172 + 0.2e1 * t27 * t55, 0.2e1 * t139 * t172 - 0.2e1 * t14 * t78 - 0.2e1 * t28 * t55, 0.2e1 * (-t13 * t79 + t27 * t54) * t139 + 0.2e1 * (-t14 * t79 + t28 * t54) * t138, 0.2e1 * t13 * t27 + 0.2e1 * t14 * t28 + 0.2e1 * t227, -0.2e1 * t49 * t22, 0.2e1 * t22 * t48 - 0.2e1 * t23 * t49, -0.2e1 * t22 * t78 + 0.2e1 * t49 * t55, 0.2e1 * t48 * t23, -0.2e1 * t23 * t78 - 0.2e1 * t48 * t55, t40, 0.2e1 * t11 * t55 + 0.2e1 * t19 * t48 + 0.2e1 * t23 * t36 + 0.2e1 * t4 * t78, -0.2e1 * t12 * t55 + 0.2e1 * t19 * t49 - 0.2e1 * t22 * t36 + 0.2e1 * t3 * t78, 0.2e1 * t11 * t22 - 0.2e1 * t12 * t23 + 0.2e1 * t3 * t48 - 0.2e1 * t4 * t49, 0.2e1 * t11 * t4 - 0.2e1 * t12 * t3 + 0.2e1 * t19 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, 0, -t195, 0, -pkin(7) * t194, pkin(7) * t195, 0, 0, 0, 0, t165, 0, -t200, 0, t68, -t69 (-t208 * t165 - t207 * t200) * pkin(2) (t207 * t69 + t208 * t68) * pkin(2), 0, 0, -t54, 0, -t55, 0, -t26, t148, t105 * t54 - t106 * t55 + t78 * t97 + t225, -t26 * t105 - t106 * t148 - t45 * t97 + t226, -t46, t29, t50, t46, t51, 0, t138 * t160 - t212, t139 * t160 + t25, t5, t102 * t5 + t103 * t26 + t168 * t95 + t226, t15, t6, -t170, t16, t33, 0, t23 * t94 + t39 * t78 + t48 * t98 + t55 * t58 + t221, -t22 * t94 + t38 * t78 + t49 * t98 - t55 * t59 + t220, t22 * t58 - t23 * t59 + t38 * t48 - t39 * t49 + t190, t11 * t39 - t12 * t38 + t19 * t94 - t3 * t59 + t36 * t98 + t4 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t98, 0.2e1 * t97, 0, -0.2e1 * t105 * t98 - 0.2e1 * t106 * t97, 0, 0, 0, 0, 0, 0, -0.2e1 * t209, 0.2e1 * t210, 0.2e1 * t233, 0.2e1 * t102 * t233 + 0.2e1 * t103 * t98, t82, t56, 0, t81, 0, 0, 0.2e1 * t218, 0.2e1 * t217, 0.2e1 * t177, -0.2e1 * t38 * t59 + 0.2e1 * t39 * t58 + 0.2e1 * t94 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, t165, 0, t191, 0, 0, 0, 0, 0, 0, t55, -t54, 0, t93, 0, 0, 0, 0, 0, 0, t51, -t50, t198 * t54, t13 * t139 + t138 * t14, 0, 0, 0, 0, 0, 0, t33, t170, -t169 + t219, -t108 * t12 - t109 * t11 - t113 * t3 + t231 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108 * t59 - t109 * t58 - t113 * t38 + t231 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t205 - 0.2e1 * t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, 0, -t55, 0, -t26, t148, 0, 0, -t46, t29, t50, t46, t51, 0, t138 * t164 - t212, t139 * t164 + t25, t5, -pkin(4) * t26 + qJ(5) * t5 + qJD(5) * t168, t15, t6, -t170, t16, t33, 0, t131 * t23 + t55 * t86 + t64 * t78 + t221, -t131 * t22 - t55 * t87 + t63 * t78 + t220, t22 * t86 - t23 * t87 + t48 * t63 - t49 * t64 + t190, t11 * t64 - t12 * t63 + t131 * t19 - t3 * t87 + t4 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, t97, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t210, t234 + t233, -pkin(4) * t98 + qJ(5) * t233 + t102 * t234, t82, t56, 0, t81, 0, 0, t203 + t218, -t204 + t217, t176 + t177, t131 * t98 - t38 * t87 + t39 * t86 + t58 * t64 - t59 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108 * t87 - t109 * t86 - t113 * t63 + t231 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, qJ(5) * t235, t82, t56, 0, t81, 0, 0, 0.2e1 * t203, -0.2e1 * t204, 0.2e1 * t176, -0.2e1 * t63 * t87 + 0.2e1 * t64 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214, -t211, 0, t26, 0, 0, 0, 0, 0, 0, t23, -t22, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, 0, 0, 0, 0, 0, t109, -t108, 0, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t108, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, -t23, t55, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, 0, -t109, 0, t39, t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, t108, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, 0, -t109, 0, t64, t63, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;