% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:16:17
% EndTime: 2019-03-09 22:16:25
% DurationCPUTime: 2.68s
% Computational Cost: add. (3662->292), mult. (8324->468), div. (0->0), fcn. (7678->8), ass. (0->177)
t154 = sin(qJ(4));
t148 = t154 * qJ(5);
t157 = cos(qJ(4));
t185 = t157 * pkin(4) + t148;
t153 = sin(qJ(6));
t246 = cos(qJ(6));
t118 = t153 * t154 + t246 * t157;
t120 = -t153 * t157 + t246 * t154;
t151 = t154 ^ 2;
t152 = t157 ^ 2;
t221 = t151 - t152;
t192 = t221 * qJD(4);
t254 = qJD(2) + qJD(3);
t196 = qJD(6) * t246;
t197 = qJD(4) * t246;
t147 = qJD(4) * t157;
t199 = t153 * t147;
t217 = qJD(6) * t153;
t253 = t157 * t217 - t199 + (-t196 + t197) * t154;
t187 = t157 * t197;
t218 = qJD(4) * t154;
t72 = t118 * qJD(6) - t153 * t218 - t187;
t155 = sin(qJ(3));
t158 = cos(qJ(3));
t159 = cos(qJ(2));
t250 = -pkin(8) - pkin(7);
t203 = qJD(2) * t250;
t189 = t159 * t203;
t156 = sin(qJ(2));
t190 = t156 * t203;
t127 = t250 * t156;
t129 = t250 * t159;
t83 = t155 * t127 - t158 * t129;
t46 = t83 * qJD(3) + t155 * t190 - t158 * t189;
t215 = t157 * qJD(5);
t252 = t185 * qJD(4) - t215;
t119 = t155 * t156 - t158 * t159;
t121 = t155 * t159 + t158 * t156;
t160 = -pkin(4) - pkin(5);
t142 = -t159 * pkin(2) - pkin(1);
t67 = t119 * pkin(3) - t121 * pkin(9) + t142;
t76 = t154 * t83;
t21 = t76 + (-pkin(10) * t121 - t67) * t157 + t160 * t119;
t227 = t121 * t154;
t236 = t154 * t67 + t157 * t83;
t31 = t119 * qJ(5) + t236;
t24 = pkin(10) * t227 + t31;
t177 = t153 * t21 + t246 * t24;
t74 = t254 * t119;
t232 = t157 * t74;
t173 = t121 * t218 + t232;
t216 = t156 * qJD(2);
t208 = pkin(2) * t216;
t75 = t254 * t121;
t37 = t75 * pkin(3) + t74 * pkin(9) + t208;
t219 = qJD(3) * t158;
t220 = qJD(3) * t155;
t45 = -t127 * t219 - t129 * t220 - t155 * t189 - t158 * t190;
t194 = t83 * t147 - t154 * t45 - t157 * t37 + t67 * t218;
t6 = t173 * pkin(10) + t160 * t75 + t194;
t233 = t154 * t74;
t174 = t121 * t147 - t233;
t108 = t119 * qJD(5);
t15 = -t67 * t147 - t154 * t37 + t157 * t45 + t83 * t218;
t68 = t75 * qJ(5);
t9 = t108 - t15 + t68;
t8 = pkin(10) * t174 + t9;
t2 = -t177 * qJD(6) - t153 * t8 + t246 * t6;
t251 = 0.2e1 * qJD(5);
t249 = pkin(9) - pkin(10);
t248 = pkin(9) * t75;
t247 = t75 * pkin(4);
t245 = pkin(9) * t119;
t140 = t155 * pkin(2) + pkin(9);
t244 = -pkin(10) + t140;
t223 = t157 * qJ(5);
t176 = t160 * t154 + t223;
t14 = -t176 * t74 + (t215 + (t160 * t157 - t148) * qJD(4)) * t121 - t46;
t82 = -t158 * t127 - t155 * t129;
t33 = t121 * t176 - t82;
t243 = t14 * t118 - t253 * t33;
t242 = t14 * t120 - t33 * t72;
t241 = t82 * t147 + t46 * t154;
t207 = pkin(2) * t220;
t145 = t154 * qJD(5);
t100 = pkin(4) * t218 - qJ(5) * t147 - t145;
t84 = -pkin(5) * t218 - t100;
t81 = t84 - t207;
t141 = -t158 * pkin(2) - pkin(3);
t109 = t141 - t185;
t149 = t157 * pkin(5);
t93 = -t109 + t149;
t240 = t81 * t118 - t253 * t93;
t239 = t81 * t120 - t93 * t72;
t211 = pkin(3) + t185;
t110 = t149 + t211;
t238 = -t110 * t253 + t84 * t118;
t237 = -t110 * t72 + t84 * t120;
t235 = t121 * t74;
t234 = t140 * t75;
t231 = t157 * t75;
t85 = t100 + t207;
t230 = -t100 - t85;
t184 = t154 * pkin(4) - t223;
t42 = t184 * t121 + t82;
t229 = qJD(4) * t42;
t228 = t119 * t140;
t224 = t154 * t157;
t222 = t141 * t147 + t154 * t207;
t214 = t159 * qJD(2);
t213 = -0.2e1 * pkin(1) * qJD(2);
t212 = 0.2e1 * t119 * t75;
t210 = pkin(3) * t218;
t209 = pkin(3) * t147;
t206 = pkin(2) * t219;
t205 = pkin(9) * t218;
t204 = pkin(9) * t147;
t38 = t42 * t218;
t69 = t82 * t218;
t198 = t154 * t147;
t113 = t244 * t157;
t195 = t157 * t67 - t76;
t193 = -0.4e1 * t121 * t224;
t191 = t154 * t206;
t143 = pkin(10) * t218;
t186 = t143 - t205;
t90 = t140 * t218 - t157 * t206;
t32 = -t119 * pkin(4) - t195;
t183 = t154 * t32 + t157 * t31;
t182 = -t154 * t31 + t157 * t32;
t181 = -t121 * t141 + t228;
t98 = (t151 + t152) * t206;
t180 = t141 * t218 - t157 * t207;
t178 = t143 - t90;
t1 = -t153 * t6 - t21 * t196 + t24 * t217 - t246 * t8;
t175 = t119 * t218 - t231;
t172 = t246 * qJ(5) + t153 * t160;
t112 = t244 * t154;
t171 = t153 * t112 + t246 * t113;
t126 = t249 * t154;
t128 = t249 * t157;
t170 = t153 * t126 + t246 * t128;
t168 = t100 * t121 + t211 * t74 - t248;
t17 = t252 * t121 - t184 * t74 + t46;
t166 = -t17 + (-t121 * t211 - t245) * qJD(4);
t165 = -t17 + (t109 * t121 - t228) * qJD(4);
t163 = qJD(4) * t113 + t191;
t162 = -t109 * t74 - t119 * t206 + t121 * t85 - t234;
t12 = t194 - t247;
t3 = t182 * qJD(4) + t12 * t154 + t9 * t157;
t161 = -t234 - t141 * t74 + (-t119 * t158 + t121 * t155) * qJD(3) * pkin(2);
t132 = 0.2e1 * t198;
t117 = -0.2e1 * t192;
t116 = t121 ^ 2;
t107 = t211 * t218;
t96 = t153 * qJD(5) + qJD(6) * t172;
t95 = qJ(5) * t217 - t246 * qJD(5) - t160 * t196;
t92 = t109 * t218;
t91 = t140 * t147 + t191;
t60 = t118 * t121;
t59 = t120 * t121;
t54 = -0.2e1 * t120 * t72;
t49 = t119 * t147 + t154 * t75;
t44 = t170 * qJD(6) + t153 * t186 - t249 * t187;
t43 = -t126 * t196 + t128 * t217 - t246 * t186 - t249 * t199;
t36 = -t121 * t192 - t74 * t224;
t30 = t119 * t72 - t75 * t120;
t29 = t118 * t75 - t119 * t253;
t28 = qJD(4) * t193 + t221 * t74;
t27 = 0.2e1 * t118 * t72 + 0.2e1 * t120 * t253;
t26 = t171 * qJD(6) + t153 * t178 - t246 * t163;
t25 = -t112 * t196 + t113 * t217 - t153 * t163 - t246 * t178;
t19 = t120 * t74 + t72 * t121;
t18 = t118 * t74 + t253 * t121;
t13 = -t18 * t120 - t60 * t72;
t4 = t118 * t18 - t19 * t120 + t253 * t60 - t59 * t72;
t5 = [0, 0, 0, 0.2e1 * t156 * t214, 0.2e1 * (-t156 ^ 2 + t159 ^ 2) * qJD(2), 0, 0, 0, t156 * t213, t159 * t213, -0.2e1 * t235, 0.2e1 * t74 * t119 - 0.2e1 * t121 * t75, 0, 0, 0, 0.2e1 * t119 * t208 + 0.2e1 * t142 * t75, 0.2e1 * t121 * t208 - 0.2e1 * t142 * t74, -0.2e1 * t116 * t198 - 0.2e1 * t152 * t235, 0.2e1 * t116 * t192 - t74 * t193, -0.2e1 * t119 * t173 + 0.2e1 * t121 * t231, -0.2e1 * t119 * t174 - 0.2e1 * t75 * t227, t212, -0.2e1 * t119 * t194 + 0.2e1 * t241 * t121 + 0.2e1 * t195 * t75 - 0.2e1 * t82 * t233, 0.2e1 * t15 * t119 - 0.2e1 * t236 * t75 - 0.2e1 * t82 * t232 + 0.2e1 * (t46 * t157 - t69) * t121, -0.2e1 * t42 * t233 - 0.2e1 * t12 * t119 - 0.2e1 * t32 * t75 + 0.2e1 * (t42 * t147 + t17 * t154) * t121, -0.2e1 * t182 * t74 + 0.2e1 * (-qJD(4) * t183 + t12 * t157 - t154 * t9) * t121, 0.2e1 * t42 * t232 + 0.2e1 * t9 * t119 + 0.2e1 * t31 * t75 + 0.2e1 * (-t17 * t157 + t38) * t121, 0.2e1 * t32 * t12 + 0.2e1 * t42 * t17 + 0.2e1 * t31 * t9, -0.2e1 * t60 * t18, -0.2e1 * t59 * t18 - 0.2e1 * t19 * t60, 0.2e1 * t119 * t18 - 0.2e1 * t75 * t60, 0.2e1 * t19 * t119 - 0.2e1 * t59 * t75, t212, -0.2e1 * t2 * t119 - 0.2e1 * (-t153 * t24 + t246 * t21) * t75 - 0.2e1 * t14 * t59 + 0.2e1 * t33 * t19, -0.2e1 * t1 * t119 + 0.2e1 * t14 * t60 + 0.2e1 * t177 * t75 - 0.2e1 * t33 * t18; 0, 0, 0, 0, 0, t214, -t216, 0, -pkin(7) * t214, pkin(7) * t216, 0, 0, -t74, -t75, 0, -t46, t45, t36, t28, t49, -t175, 0, t69 + (-qJD(4) * t181 - t46) * t157 + t161 * t154, t157 * t161 + t181 * t218 + t241, t154 * t162 + t157 * t165 + t38, t3, t165 * t154 + (-t162 - t229) * t157, t17 * t109 + t140 * t3 + t183 * t206 + t42 * t85, t13, t4, t30, t29, 0, t26 * t119 - (t246 * t112 - t153 * t113) * t75 - t81 * t59 + t93 * t19 + t243, -t25 * t119 + t171 * t75 - t93 * t18 + t81 * t60 + t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t207, -0.2e1 * t206, t132, t117, 0, 0, 0, 0.2e1 * t180, 0.2e1 * t222, -0.2e1 * t85 * t157 + 0.2e1 * t92, 0.2e1 * t98, -0.2e1 * t109 * t147 - 0.2e1 * t85 * t154, 0.2e1 * t109 * t85 + 0.2e1 * t140 * t98, t54, t27, 0, 0, 0, 0.2e1 * t240, 0.2e1 * t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t75, 0, -t46, t45, t36, t28, t49, -t175, 0, t69 + (pkin(3) * t74 - t248) * t154 + (-t46 + (-pkin(3) * t121 - t245) * qJD(4)) * t157, pkin(3) * t173 + pkin(9) * t175 + t241, t154 * t168 + t157 * t166 + t38, t3, t166 * t154 + (-t168 - t229) * t157, pkin(9) * t3 + t42 * t100 - t17 * t211, t13, t4, t30, t29, 0, t44 * t119 - (t246 * t126 - t153 * t128) * t75 - t84 * t59 + t110 * t19 + t243, -t110 * t18 - t43 * t119 + t170 * t75 + t84 * t60 + t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, -t206, t132, t117, 0, 0, 0, t180 - t210, -t209 + t222, t230 * t157 - t107 + t92, t98, t230 * t154 + (-t109 + t211) * t147, pkin(9) * t98 + t109 * t100 - t211 * t85, t54, t27, 0, 0, 0, t238 + t240, t237 + t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t117, 0, 0, 0, -0.2e1 * t210, -0.2e1 * t209, -0.2e1 * t100 * t157 - 0.2e1 * t107, 0, -0.2e1 * t100 * t154 + 0.2e1 * t147 * t211, -0.2e1 * t211 * t100, t54, t27, 0, 0, 0, 0.2e1 * t238, 0.2e1 * t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, -t174, t75, -t194, t15, -t194 + 0.2e1 * t247, t185 * t74 + (qJD(4) * t184 - t145) * t121, 0.2e1 * t108 - t15 + 0.2e1 * t68, -t12 * pkin(4) + t9 * qJ(5) + t31 * qJD(5), 0, 0, t18, t19, t75, t96 * t119 - (-t153 * qJ(5) + t246 * t160) * t75 - t2, -t95 * t119 + t172 * t75 - t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, -t218, 0, -t91, t90, -t91, -t252, -t90, -t140 * t252 - t184 * t206, 0, 0, t72, -t253, 0, t26, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, -t218, 0, -t204, t205, -t204, -t252, -t205, -t252 * pkin(9), 0, 0, t72, -t253, 0, t44, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, qJ(5) * t251, 0, 0, 0, 0, 0, 0.2e1 * t96, -0.2e1 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t173, 0, t12, 0, 0, 0, 0, 0, t119 * t217 - t246 * t75, t119 * t196 + t153 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, t91, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, t204, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, -t75, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t253, 0, -t26, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t253, 0, -t44, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, -t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
