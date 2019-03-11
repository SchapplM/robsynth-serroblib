% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRRPR2
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
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:59:24
% EndTime: 2019-03-09 21:59:36
% DurationCPUTime: 5.74s
% Computational Cost: add. (13623->381), mult. (29577->624), div. (0->0), fcn. (29541->10), ass. (0->208)
t162 = sin(pkin(11));
t160 = t162 ^ 2;
t163 = cos(pkin(11));
t161 = t163 ^ 2;
t252 = t160 + t161;
t294 = t252 * qJD(5);
t295 = 0.2e1 * t294;
t285 = cos(qJ(3));
t245 = t285 * pkin(2);
t154 = t245 + pkin(3);
t165 = sin(qJ(4));
t227 = t285 * qJD(3);
t216 = pkin(2) * t227;
t284 = cos(qJ(4));
t232 = qJD(4) * t284;
t282 = sin(qJ(3));
t233 = qJD(3) * t282;
t90 = (qJD(4) * t282 + t233) * pkin(2) * t165 - t154 * t232 - t284 * t216;
t89 = qJD(5) - t90;
t293 = t252 * t89;
t249 = qJD(4) * t165;
t156 = pkin(3) * t249;
t211 = t284 * t282;
t91 = ((t165 * t285 + t211) * qJD(3) + qJD(4) * t211) * pkin(2) + t154 * t249;
t292 = -t91 - t156;
t214 = pkin(3) * t232;
t143 = t214 + qJD(5);
t291 = t252 * t143;
t166 = sin(qJ(2));
t167 = cos(qJ(2));
t197 = -t166 * t285 - t167 * t282;
t191 = t197 * qJD(3);
t178 = qJD(2) * t197 + t191;
t238 = t282 * t166;
t241 = t285 * t167;
t198 = t238 - t241;
t187 = t284 * t198;
t254 = qJD(2) * t238 + t166 * t233;
t288 = -(qJD(2) * t285 + t227) * t167 + t254;
t43 = qJD(4) * t187 - t165 * t178 - t197 * t249 + t284 * t288;
t192 = t165 * t198;
t44 = -qJD(4) * t192 - t165 * t288 - t178 * t284 - t197 * t232;
t95 = -t165 * t197 + t187;
t96 = -t197 * t284 - t192;
t209 = t43 * t95 - t44 * t96;
t290 = 0.2e1 * t209;
t286 = -pkin(8) - pkin(7);
t164 = sin(qJ(6));
t283 = cos(qJ(6));
t239 = t283 * t163;
t289 = -t162 * t164 + t239;
t280 = pkin(3) * t165;
t279 = t163 * pkin(5);
t278 = t166 * pkin(2);
t213 = t286 * t282;
t131 = t166 * t213;
t139 = t286 * t167;
t184 = -qJD(3) * t131 + t139 * t227;
t205 = qJD(2) * t213;
t212 = t285 * t286;
t206 = qJD(2) * t212;
t172 = pkin(9) * t288 - t166 * t205 + t167 * t206 + t184;
t132 = t166 * t212;
t64 = -qJD(3) * t132 - t139 * t233 - t166 * t206 - t167 * t205;
t47 = pkin(9) * t178 - t64;
t105 = t139 * t282 + t132;
t80 = pkin(9) * t197 + t105;
t106 = -t139 * t285 + t131;
t81 = -pkin(9) * t198 + t106;
t50 = t165 * t80 + t284 * t81;
t23 = qJD(4) * t50 + t165 * t47 - t172 * t284;
t49 = t165 * t81 - t284 * t80;
t277 = t49 * t23;
t276 = t49 * t91;
t275 = t91 * t96;
t118 = pkin(2) * t211 + t154 * t165;
t114 = qJ(5) + t118;
t274 = -pkin(10) - t114;
t151 = qJ(5) + t280;
t273 = -pkin(10) - t151;
t272 = qJ(5) + pkin(10);
t240 = t283 * t162;
t259 = t164 * t163;
t129 = t240 + t259;
t120 = t129 * qJD(6);
t266 = t162 * t43;
t19 = -pkin(5) * t266 + t23;
t264 = t162 * t96;
t38 = pkin(5) * t264 + t49;
t271 = t120 * t38 - t19 * t289;
t231 = qJD(6) * t283;
t248 = qJD(6) * t164;
t236 = t162 * t248;
t119 = -t163 * t231 + t236;
t270 = -t119 * t38 + t129 * t19;
t155 = -t167 * pkin(2) - pkin(1);
t110 = pkin(3) * t198 + t155;
t181 = t95 * pkin(4) + t110;
t180 = -t96 * qJ(5) + t181;
t28 = t162 * t180 + t163 * t50;
t243 = t282 * pkin(2);
t117 = t154 * t284 - t165 * t243;
t115 = -pkin(4) - t117;
t109 = t115 - t279;
t269 = t109 * t120 - t289 * t91;
t268 = -t109 * t119 + t129 * t91;
t265 = t162 * t50;
t263 = t163 * t43;
t262 = t23 * t163;
t261 = t91 * t163;
t244 = t284 * pkin(3);
t153 = -t244 - pkin(4);
t137 = t153 - t279;
t258 = t120 * t137 - t156 * t289;
t257 = -t119 * t137 + t129 * t156;
t251 = qJD(2) * t166;
t250 = qJD(2) * t167;
t37 = 0.2e1 * t95 * t44;
t247 = -0.2e1 * t96 * t43;
t246 = -0.2e1 * pkin(1) * qJD(2);
t242 = pkin(2) * t251;
t39 = t162 * t263;
t235 = t166 * t250;
t176 = t95 * pkin(5) - t265 + (-t272 * t96 + t181) * t163;
t174 = t283 * t176;
t25 = -pkin(10) * t264 + t28;
t13 = -t164 * t25 + t174;
t175 = t164 * t176;
t14 = t25 * t283 + t175;
t190 = t232 * t80 - t249 * t81 + t284 * t47;
t173 = t165 * (-(t167 * t227 - t254) * pkin(9) + t184) + t190;
t177 = t165 * (-t131 + (-pkin(9) + t286) * t241);
t186 = -pkin(3) * t197 + t278;
t170 = t163 * t173 + (t162 * t186 + t163 * t177) * qJD(2);
t188 = pkin(3) * t191;
t208 = t44 * pkin(4) - t96 * qJD(5);
t168 = (t272 * t43 - t188 + t208) * t162 + t170;
t171 = -t165 * t172 - t190;
t199 = t43 * qJ(5) + t208;
t169 = t162 * t171 + t163 * (-pkin(3) * t178 + t199 + t242) + pkin(10) * t263 + t44 * pkin(5);
t3 = -qJD(6) * t174 - t164 * t169 - t168 * t283 + t248 * t25;
t4 = -qJD(6) * t175 - t164 * t168 + t169 * t283 - t231 * t25;
t234 = t119 * t13 - t120 * t14 - t129 * t4 - t289 * t3;
t230 = t272 * t162;
t229 = t274 * t162;
t228 = t273 * t162;
t226 = t283 * qJD(5);
t159 = t163 * pkin(10);
t104 = t114 * t163 + t159;
t202 = t283 * t229;
t33 = -t89 * t239 - qJD(6) * t202 + (qJD(6) * t104 + t162 * t89) * t164;
t34 = -t104 * t231 - t89 * t259 + (-t248 * t274 - t283 * t89) * t162;
t60 = -t104 * t164 + t202;
t61 = t104 * t283 + t164 * t229;
t224 = t119 * t60 - t120 * t61 - t129 * t34 - t289 * t33;
t124 = t151 * t163 + t159;
t201 = t283 * t228;
t57 = -qJD(6) * t201 - t143 * t239 + (qJD(6) * t124 + t143 * t162) * t164;
t58 = -t124 * t231 - t143 * t259 + (-t143 * t283 - t248 * t273) * t162;
t87 = -t164 * t124 + t201;
t88 = t124 * t283 + t164 * t228;
t223 = t119 * t87 - t120 * t88 - t129 * t58 - t289 * t57;
t138 = qJ(5) * t163 + t159;
t203 = t283 * t230;
t100 = -t164 * t138 - t203;
t101 = t138 * t283 - t164 * t230;
t70 = qJD(6) * t203 - t163 * t226 + (qJD(5) * t162 + qJD(6) * t138) * t164;
t71 = -t138 * t231 - qJD(5) * t259 + (t248 * t272 - t226) * t162;
t222 = t100 * t119 - t101 * t120 - t129 * t71 - t289 * t70;
t220 = t252 * t151;
t218 = t49 * t156;
t217 = t163 * t156;
t215 = pkin(2) * t233;
t210 = t23 * t96 - t49 * t43;
t179 = -(-t166 * t227 - t167 * t233) * pkin(3) + t199;
t10 = t162 * t179 + t170;
t9 = -t162 * t173 + t163 * t179 + (-t162 * t177 + t163 * t186) * qJD(2);
t5 = t10 * t163 - t162 * t9;
t27 = t163 * t180 - t265;
t207 = -t162 * t27 + t163 * t28;
t200 = pkin(4) * t43 - qJ(5) * t44 - qJD(5) * t95;
t193 = t155 * t197;
t189 = -t114 * t44 - t115 * t43 - t89 * t95 + t275;
t185 = -t143 * t95 - t151 * t44 - t153 * t43 + t156 * t96;
t152 = -pkin(4) - t279;
t141 = t162 * t156;
t108 = t152 * t120;
t107 = t152 * t119;
t99 = -0.2e1 * t129 * t119;
t98 = -0.2e1 * t289 * t120;
t86 = qJD(2) * t186 - t188;
t82 = t91 * t162;
t65 = -t106 * qJD(3) + (t167 * t212 - t131) * qJD(2);
t59 = -0.2e1 * t119 * t289 - 0.2e1 * t120 * t129;
t54 = t289 * t96;
t53 = t129 * t96;
t42 = t163 * t44;
t41 = t162 * t44;
t30 = -t120 * t95 + t289 * t44;
t29 = -t119 * t95 + t129 * t44;
t26 = (t160 - t161) * t43;
t22 = t23 * t162;
t21 = -t43 * t240 - t96 * t236 + (-t164 * t43 + t231 * t96) * t163;
t20 = t96 * t120 + t289 * t43;
t16 = t120 * t53 - t21 * t289;
t15 = -t119 * t54 - t129 * t20;
t6 = t119 * t53 - t120 * t54 - t129 * t21 - t20 * t289;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t235, 0.2e1 * (-t166 ^ 2 + t167 ^ 2) * qJD(2), 0, -0.2e1 * t235, 0, 0, t166 * t246, t167 * t246, 0, 0, 0.2e1 * t197 * t288, -0.2e1 * t178 * t197 + 0.2e1 * t198 * t288, 0, -0.2e1 * t198 * t178, 0, 0, -0.2e1 * qJD(3) * t193 + 0.2e1 * (t198 * t278 - t193) * qJD(2), -0.2e1 * t155 * t288 - 0.2e1 * t197 * t242, 0.2e1 * t105 * t288 + 0.2e1 * t106 * t178 + 0.2e1 * t197 * t65 + 0.2e1 * t198 * t64, 0.2e1 * t105 * t65 - 0.2e1 * t106 * t64 + 0.2e1 * t155 * t242, t247, t290, 0, t37, 0, 0, 0.2e1 * t110 * t44 + 0.2e1 * t86 * t95, -0.2e1 * t110 * t43 + 0.2e1 * t86 * t96, 0.2e1 * t171 * t95 - 0.2e1 * t50 * t44 + 0.2e1 * t210, 0.2e1 * t110 * t86 + 0.2e1 * t277 + 0.2e1 * (t173 + t165 * (-t131 + (-pkin(9) * t285 + t212) * t167) * qJD(2)) * t50, t161 * t247, 0.4e1 * t96 * t39, -0.2e1 * t209 * t163, t160 * t247, t162 * t290, t37, 0.2e1 * t162 * t210 + 0.2e1 * t27 * t44 + 0.2e1 * t9 * t95, -0.2e1 * t10 * t95 + 0.2e1 * t163 * t210 - 0.2e1 * t28 * t44, 0.2e1 * (t27 * t43 - t9 * t96) * t163 + 0.2e1 * (-t10 * t96 + t28 * t43) * t162, 0.2e1 * t10 * t28 + 0.2e1 * t27 * t9 + 0.2e1 * t277, -0.2e1 * t54 * t20, 0.2e1 * t20 * t53 - 0.2e1 * t21 * t54, -0.2e1 * t20 * t95 + 0.2e1 * t44 * t54, 0.2e1 * t53 * t21, -0.2e1 * t21 * t95 - 0.2e1 * t44 * t53, t37, 0.2e1 * t13 * t44 + 0.2e1 * t19 * t53 + 0.2e1 * t21 * t38 + 0.2e1 * t4 * t95, -0.2e1 * t14 * t44 + 0.2e1 * t19 * t54 - 0.2e1 * t20 * t38 + 0.2e1 * t3 * t95, 0.2e1 * t13 * t20 - 0.2e1 * t14 * t21 + 0.2e1 * t3 * t53 - 0.2e1 * t4 * t54, 0.2e1 * t13 * t4 - 0.2e1 * t14 * t3 + 0.2e1 * t19 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250, 0, -t251, 0, -pkin(7) * t250, pkin(7) * t251, 0, 0, 0, 0, -t288, 0, t178, 0, t65, t64, t178 * t243 - t197 * t215 - t198 * t216 + t245 * t288 (-t64 * t282 + t285 * t65 + (-t105 * t282 + t106 * t285) * qJD(3)) * pkin(2), 0, 0, -t43, 0, -t44, 0, -t23, t171, t117 * t43 - t118 * t44 + t90 * t95 + t275, -t23 * t117 - t118 * t171 - t50 * t90 + t276, -t39, t26, t41, t39, t42, 0, t162 * t189 - t262, t163 * t189 + t22, t5, t114 * t5 + t115 * t23 + t207 * t89 + t276, t15, t6, t29, t16, t30, 0, t109 * t21 + t34 * t95 + t44 * t60 + t53 * t91 + t271, -t109 * t20 + t33 * t95 - t44 * t61 + t54 * t91 + t270, t20 * t60 - t21 * t61 + t33 * t53 - t34 * t54 + t234, t109 * t19 + t13 * t34 - t14 * t33 - t3 * t61 + t38 * t91 + t4 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t215, -0.2e1 * t216, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t91, 0.2e1 * t90, 0, -0.2e1 * t117 * t91 - 0.2e1 * t118 * t90, 0, 0, 0, 0, 0, 0, -0.2e1 * t261, 0.2e1 * t82, 0.2e1 * t293, 0.2e1 * t114 * t293 + 0.2e1 * t115 * t91, t99, t59, 0, t98, 0, 0, 0.2e1 * t269, 0.2e1 * t268, 0.2e1 * t224, 0.2e1 * t109 * t91 - 0.2e1 * t33 * t61 + 0.2e1 * t34 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t288, 0, t178, 0, t65, t64, 0, 0, 0, 0, -t43, 0, -t44, 0, -t23, t171 (t284 * t43 - t165 * t44 + (t165 * t96 - t284 * t95) * qJD(4)) * pkin(3), -t171 * t280 + t214 * t50 - t23 * t244 + t218, -t39, t26, t41, t39, t42, 0, t162 * t185 - t262, t163 * t185 + t22, t5, t143 * t207 + t151 * t5 + t153 * t23 + t218, t15, t6, t29, t16, t30, 0, t137 * t21 + t156 * t53 + t44 * t87 + t58 * t95 + t271, -t137 * t20 + t156 * t54 - t44 * t88 + t57 * t95 + t270, t20 * t87 - t21 * t88 + t53 * t57 - t54 * t58 + t234, t13 * t58 + t137 * t19 - t14 * t57 + t156 * t38 - t3 * t88 + t4 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, -t216, 0, 0, 0, 0, 0, 0, 0, 0, t292, -t214 + t90, 0 (-t284 * t91 - t165 * t90 + (-t117 * t165 + t118 * t284) * qJD(4)) * pkin(3), 0, 0, 0, 0, 0, 0, t292 * t163, t141 + t82, t291 + t293, t114 * t291 + t115 * t156 + t153 * t91 + t220 * t89, t99, t59, 0, t98, 0, 0, t258 + t269, t257 + t268, t223 + t224, t109 * t156 + t137 * t91 - t33 * t88 + t34 * t87 - t57 * t61 + t58 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t156, -0.2e1 * t214, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t217, 0.2e1 * t141, 0.2e1 * t291, 0.2e1 * t143 * t220 + 0.2e1 * t153 * t156, t99, t59, 0, t98, 0, 0, 0.2e1 * t258, 0.2e1 * t257, 0.2e1 * t223, 0.2e1 * t137 * t156 - 0.2e1 * t57 * t88 + 0.2e1 * t58 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, -t44, 0, -t23, t171, 0, 0, -t39, t26, t41, t39, t42, 0, t162 * t200 - t262, t163 * t200 + t22, t5, -pkin(4) * t23 + qJ(5) * t5 + qJD(5) * t207, t15, t6, t29, t16, t30, 0, t100 * t44 + t152 * t21 + t71 * t95 + t271, -t101 * t44 - t152 * t20 + t70 * t95 + t270, t100 * t20 - t101 * t21 + t53 * t70 - t54 * t71 + t234, t100 * t4 - t101 * t3 + t13 * t71 - t14 * t70 + t152 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t90, 0, 0, 0, 0, 0, 0, 0, 0, -t261, t82, t294 + t293, -pkin(4) * t91 + qJ(5) * t293 + t114 * t294, t99, t59, 0, t98, 0, 0, t108 + t269, -t107 + t268, t222 + t224, t100 * t34 - t101 * t33 + t152 * t91 + t60 * t71 - t61 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t214, 0, 0, 0, 0, 0, 0, 0, 0, -t217, t141, t294 + t291, -pkin(4) * t156 + qJ(5) * t291 + t151 * t294, t99, t59, 0, t98, 0, 0, t108 + t258, -t107 + t257, t222 + t223, t100 * t58 - t101 * t57 + t152 * t156 - t70 * t88 + t71 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, qJ(5) * t295, t99, t59, 0, t98, 0, 0, 0.2e1 * t108, -0.2e1 * t107, 0.2e1 * t222, 0.2e1 * t100 * t71 - 0.2e1 * t101 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t266, -t263, 0, t23, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, 0, 0, 0, 0, 0, 0, t120, -t119, 0, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, 0, 0, 0, 0, 0, 0, t120, -t119, 0, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t119, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, -t21, t44, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, 0, -t120, 0, t34, t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, 0, -t120, 0, t58, t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, 0, -t120, 0, t71, t70, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
