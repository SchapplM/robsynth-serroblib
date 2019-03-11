% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:55
% EndTime: 2019-03-09 03:20:00
% DurationCPUTime: 3.40s
% Computational Cost: add. (5607->371), mult. (14262->456), div. (0->0), fcn. (10448->6), ass. (0->202)
t156 = cos(pkin(9));
t249 = cos(qJ(3));
t200 = t249 * t156;
t185 = qJD(1) * t200;
t155 = sin(pkin(9));
t158 = sin(qJ(3));
t218 = t155 * t158;
t199 = qJD(1) * t218;
t122 = -t185 + t199;
t157 = sin(qJ(5));
t159 = cos(qJ(5));
t100 = qJD(3) * t157 - t122 * t159;
t131 = t155 * t249 + t156 * t158;
t257 = t131 * qJD(1);
t264 = qJD(5) + t257;
t267 = t100 * t264;
t127 = t131 * qJD(3);
t164 = qJD(1) * t127;
t206 = qJD(5) * t159;
t207 = qJD(5) * t157;
t79 = qJD(3) * t207 - t122 * t206 - t157 * t164;
t270 = -t79 + t267;
t269 = t79 + t267;
t102 = qJD(3) * t159 + t122 * t157;
t266 = t102 * t264;
t80 = qJD(5) * t102 - t159 * t164;
t41 = -t80 + t266;
t150 = -t156 * pkin(2) - pkin(1);
t134 = qJD(1) * t150 + qJD(2);
t166 = -qJ(4) * t257 + t134;
t250 = pkin(3) + pkin(8);
t57 = t122 * t250 + t166;
t246 = pkin(7) + qJ(2);
t137 = t246 * t155;
t132 = qJD(1) * t137;
t138 = t246 * t156;
t133 = qJD(1) * t138;
t94 = t132 * t249 + t158 * t133;
t260 = qJD(4) + t94;
t213 = pkin(4) * t257 + t260;
t64 = -qJD(3) * t250 + t213;
t24 = t157 * t64 + t159 * t57;
t17 = -qJ(6) * t100 + t24;
t268 = t264 * t17;
t143 = qJD(3) * t185;
t208 = qJD(3) * t158;
t198 = t155 * t208;
t109 = qJD(1) * t198 - t143;
t165 = t257 * qJD(3);
t62 = pkin(3) * t165 + t109 * qJ(4) - qJD(4) * t257;
t44 = pkin(8) * t164 + t62;
t197 = qJD(2) * t249;
t184 = qJD(1) * t197;
t202 = qJD(1) * qJD(2);
t195 = t158 * t202;
t196 = qJD(3) * t249;
t69 = -t132 * t208 + t133 * t196 + t155 * t184 + t156 * t195;
t49 = -pkin(4) * t109 + t69;
t7 = -qJD(5) * t24 - t157 * t44 + t159 * t49;
t254 = t24 * t264 + t7;
t194 = -t157 * t49 - t159 * t44 - t206 * t64 + t207 * t57;
t23 = -t157 * t57 + t159 * t64;
t181 = -t23 * t264 - t194;
t16 = -qJ(6) * t102 + t23;
t15 = pkin(5) * t264 + t16;
t175 = qJ(6) * t80 + t194;
t2 = -qJD(6) * t100 - t175;
t255 = t15 * t264 - t2;
t162 = qJ(6) * t79 + t7;
t248 = pkin(5) * t109;
t1 = -qJD(6) * t102 + t162 - t248;
t256 = t1 + t268;
t265 = -t157 * t255 + t159 * t256;
t105 = t159 * t109;
t190 = t157 * t264;
t174 = -t190 * t264 - t105;
t126 = -t156 * t196 + t198;
t96 = t109 * t131;
t262 = -t126 * t257 - t96;
t189 = t159 * t264;
t120 = t257 ^ 2;
t251 = t122 ^ 2;
t259 = -t251 - t120;
t258 = -t251 + t120;
t154 = qJD(3) * qJ(4);
t95 = -t132 * t158 + t133 * t249;
t78 = -pkin(4) * t122 + t95;
t70 = t154 + t78;
t253 = t109 * t250 + t264 * t70;
t83 = t158 * (qJD(2) * t155 + qJD(3) * t138) + t137 * t196 - t156 * t197;
t252 = t102 ^ 2;
t97 = t137 * t249 + t138 * t158;
t247 = t69 * t97;
t245 = t15 - t16;
t226 = t122 * qJ(4);
t73 = t250 * t257 + t226;
t34 = t157 * t78 + t159 * t73;
t130 = -t200 + t218;
t176 = -t131 * qJ(4) + t150;
t76 = t130 * t250 + t176;
t87 = pkin(4) * t131 + t97;
t39 = t157 * t87 + t159 * t76;
t82 = pkin(3) * t122 + t166;
t244 = t257 * t82;
t153 = qJD(3) * qJD(4);
t186 = t132 * t196 + t133 * t208 + t155 * t195 - t156 * t184;
t65 = -t153 + t186;
t45 = -pkin(4) * t164 - t65;
t25 = t80 * pkin(5) + t45;
t243 = t157 * t25;
t242 = t157 * t45;
t241 = t159 * t25;
t240 = t159 * t45;
t239 = t159 * t79;
t238 = t80 * t157;
t204 = t159 * qJD(6);
t215 = qJ(6) + t250;
t234 = qJ(6) * t257;
t72 = t159 * t78;
t237 = t215 * t207 - t204 + t122 * pkin(5) - t72 - (-t73 - t234) * t157;
t136 = t215 * t159;
t236 = -qJD(5) * t136 - t157 * qJD(6) - t159 * t234 - t34;
t201 = -pkin(5) * t159 - pkin(4);
t235 = pkin(5) * t206 - t201 * t257 + t260;
t233 = qJD(3) * t83;
t231 = t102 * t100;
t229 = t102 * t122;
t227 = t264 * t122;
t225 = t122 * t100;
t224 = t122 * t257;
t223 = t122 * t127;
t221 = t127 * t157;
t220 = t127 * t159;
t219 = t130 * t159;
t217 = t157 * t109;
t98 = -t137 * t158 + t138 * t249;
t84 = qJD(2) * t131 + qJD(3) * t98;
t216 = t84 * qJD(3);
t211 = t155 ^ 2 + t156 ^ 2;
t210 = qJD(3) * t126;
t209 = qJD(3) * t127;
t205 = qJD(5) * t250;
t193 = pkin(5) * t100 + qJD(6);
t192 = -qJ(6) * t130 - t76;
t191 = t211 * qJD(1) ^ 2;
t180 = -t157 * t23 + t159 * t24;
t179 = qJ(4) * t126 - qJD(4) * t131;
t178 = 0.2e1 * t211 * t202;
t52 = t127 * t250 + t179;
t59 = -t126 * pkin(4) + t84;
t11 = t157 * t59 + t159 * t52 + t206 * t87 - t207 * t76;
t173 = t130 * t206 + t221;
t172 = t130 * t207 - t220;
t171 = qJD(3) * t95 - t69;
t170 = -qJD(3) * t94 + t186;
t169 = t109 * t130 + t126 * t122 - t127 * t257;
t168 = -t189 * t264 + t217;
t163 = -t97 * t109 + t83 * t122 + t69 * t131 + t257 * t84;
t58 = -t127 * pkin(4) - t83;
t149 = pkin(5) * t157 + qJ(4);
t135 = t215 * t157;
t112 = qJD(3) * t122;
t99 = t100 ^ 2;
t93 = pkin(3) * t130 + t176;
t92 = pkin(3) * t257 + t226;
t91 = -t154 - t95;
t90 = t109 - t112;
t89 = -qJD(3) * pkin(3) + t260;
t88 = -pkin(4) * t130 + t98;
t86 = t159 * t87;
t75 = pkin(3) * t127 + t179;
t74 = t159 * t80;
t63 = t130 * t201 + t98;
t56 = t159 * t59;
t54 = -t126 * t264 - t96;
t48 = -t99 + t252;
t43 = t193 + t70;
t38 = -t157 * t76 + t86;
t37 = pkin(5) * t172 + t58;
t36 = -qJD(3) * t102 + t168;
t35 = -qJD(3) * t100 + t174;
t33 = -t157 * t73 + t72;
t32 = qJ(6) * t219 + t39;
t31 = t168 - t225;
t30 = t168 + t225;
t29 = t229 - t174;
t28 = t174 + t229;
t27 = t159 * t267 + t238;
t26 = -t102 * t190 - t239;
t21 = t131 * pkin(5) + t157 * t192 + t86;
t19 = t100 * t172 - t219 * t80;
t18 = -t79 * t157 * t130 + t102 * t173;
t14 = -t102 * t126 - t130 * t217 - t131 * t79 + t173 * t264;
t13 = t100 * t126 - t105 * t130 - t80 * t131 - t172 * t264;
t12 = -qJD(5) * t39 - t157 * t52 + t56;
t10 = t157 * t41 - t159 * t270;
t9 = -t102 * t189 + t157 * t269 - t74;
t8 = t157 * t270 + t159 * t266 - t74;
t5 = -qJ(6) * t172 + t130 * t204 + t11;
t4 = (-t100 * t157 + t102 * t159) * t127 + (-t238 - t239 + (-t100 * t159 - t102 * t157) * qJD(5)) * t130;
t3 = -pkin(5) * t126 + t56 + t192 * t206 + (-qJ(6) * t127 - qJD(5) * t87 - qJD(6) * t130 - t52) * t157;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, qJ(2) * t178, t262, -t131 * t165 + t169, -t210, t130 * t165 + t223, -t209, 0, t134 * t127 + t150 * t164 - t216, -t109 * t150 - t126 * t134 + t233, -t94 * t126 - t95 * t127 + t130 * t186 - t165 * t98 + t163, -t186 * t98 - t83 * t95 + t84 * t94 + t247, 0, t210, t209, t262, -t131 * t164 + t169, t130 * t164 + t223, -t89 * t126 + t91 * t127 + t65 * t130 - t164 * t98 + t163, -t75 * t122 - t82 * t127 - t62 * t130 - t165 * t93 + t216, t109 * t93 + t126 * t82 - t131 * t62 - t257 * t75 - t233, t62 * t93 - t65 * t98 + t75 * t82 + t83 * t91 + t84 * t89 + t247, t18, t4, t14, t19, t13, t54, -t70 * t220 + t100 * t58 - t109 * t38 + t264 * t12 - t126 * t23 + t131 * t7 + t80 * t88 + (t207 * t70 - t240) * t130, t70 * t221 + t102 * t58 + t109 * t39 - t11 * t264 + t126 * t24 + t131 * t194 - t79 * t88 + (t206 * t70 + t242) * t130, -t100 * t11 - t102 * t12 + t38 * t79 - t39 * t80 + t180 * t127 + (-t157 * t7 - t159 * t194 + (-t157 * t24 - t159 * t23) * qJD(5)) * t130, t11 * t24 + t12 * t23 - t194 * t39 + t38 * t7 + t45 * t88 + t58 * t70, t18, t4, t14, t19, t13, t54, -t43 * t220 + t1 * t131 + t100 * t37 - t109 * t21 + t264 * t3 - t126 * t15 + t63 * t80 + (t207 * t43 - t241) * t130, t43 * t221 + t102 * t37 + t109 * t32 - t264 * t5 + t126 * t17 - t131 * t2 - t63 * t79 + (t206 * t43 + t243) * t130, -t100 * t5 - t102 * t3 + t21 * t79 - t32 * t80 + (-t15 * t157 + t159 * t17) * t127 + (-t1 * t157 + t159 * t2 + (-t15 * t159 - t157 * t17) * qJD(5)) * t130, t1 * t21 + t15 * t3 + t17 * t5 + t2 * t32 + t25 * t63 + t37 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, -qJ(2) * t191, 0, 0, 0, 0, 0, 0, 0.2e1 * t165, t143 + (-t122 - t199) * qJD(3), t259, t122 * t95 - t257 * t94, 0, 0, 0, 0, 0, 0, t259, -0.2e1 * t165, t109 + t112, -t91 * t122 - t257 * t89 + t62, 0, 0, 0, 0, 0, 0, t30, t29, t8, t122 * t70 - t157 * t254 + t159 * t181, 0, 0, 0, 0, 0, 0, t30, t29, t8, t122 * t43 - t157 * t256 - t159 * t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, t258, t143 + (t122 - t199) * qJD(3), -t224, 0, 0, -t134 * t257 + t171, t122 * t134 + t170, 0, 0, 0, t90, 0, t224, t258, -t224, pkin(3) * t109 - qJ(4) * t164 - (t91 + t95) * t257 + (t89 - t260) * t122, t122 * t92 - t171 + t244, -t122 * t82 + t257 * t92 + 0.2e1 * t153 - t170, -pkin(3) * t69 - qJ(4) * t65 - t260 * t91 - t82 * t92 - t89 * t95, t26, t9, t28, t27, t31, t227, qJ(4) * t80 + t122 * t23 + t242 + (t157 * t205 - t33) * t264 + t213 * t100 + t253 * t159, -qJ(4) * t79 - t122 * t24 + t240 + (t159 * t205 + t34) * t264 + t213 * t102 - t253 * t157, t100 * t34 + t102 * t33 + (-t257 * t24 - t250 * t79 - t7 + (t100 * t250 - t24) * qJD(5)) * t159 + (t257 * t23 + t250 * t80 + t194 + (-t102 * t250 + t23) * qJD(5)) * t157, qJ(4) * t45 - t23 * t33 - t24 * t34 + t213 * t70 - (qJD(5) * t180 - t157 * t194 + t159 * t7) * t250, t26, t9, t28, t27, t31, t227, t100 * t235 + t109 * t136 + t122 * t15 + t149 * t80 + t189 * t43 + t237 * t264 + t243, t102 * t235 - t109 * t135 - t122 * t17 - t149 * t79 - t190 * t43 - t236 * t264 + t241, -t100 * t236 - t102 * t237 + t135 * t80 - t136 * t79 - t265, -t1 * t136 - t135 * t2 + t149 * t25 + t15 * t237 + t17 * t236 + t235 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t224, -qJD(3) ^ 2 - t120, qJD(3) * t91 + t244 + t69, 0, 0, 0, 0, 0, 0, t35, t36, t10, -qJD(3) * t70 + t157 * t181 + t159 * t254, 0, 0, 0, 0, 0, 0, t35, t36, t10, -qJD(3) * t43 + t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, t48, t270, -t231, t41, -t109, -t102 * t70 + t254, t100 * t70 - t181, 0, 0, t231, t48, t270, -t231, t41, -t109, -0.2e1 * t248 + t268 + (-t193 - t43) * t102 + t162, -pkin(5) * t252 + t264 * t16 + (qJD(6) + t43) * t100 + t175, pkin(5) * t79 - t100 * t245, t245 * t17 + (-t102 * t43 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 + t266, -t269, -t99 - t252, t17 * t100 + t15 * t102 + t25;];
tauc_reg  = t6;
