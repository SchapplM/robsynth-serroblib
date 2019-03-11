% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:10
% EndTime: 2019-03-09 18:23:25
% DurationCPUTime: 5.71s
% Computational Cost: add. (7640->353), mult. (15879->565), div. (0->0), fcn. (15119->8), ass. (0->194)
t144 = sin(qJ(2));
t146 = cos(qJ(2));
t268 = cos(qJ(3));
t233 = t268 * t146;
t266 = sin(qJ(3));
t111 = t266 * t144 - t233;
t290 = -0.4e1 * t111;
t271 = pkin(3) + pkin(9);
t231 = t266 * t146;
t178 = -t268 * t144 - t231;
t173 = t178 * qJD(2);
t276 = t178 * qJD(3);
t161 = t276 + t173;
t160 = t111 * t161;
t289 = -0.2e1 * t160;
t270 = -pkin(8) - pkin(7);
t208 = t270 * t266;
t191 = qJD(2) * t208;
t209 = t270 * t268;
t192 = qJD(2) * t209;
t89 = -t144 * t209 - t270 * t231;
t55 = t89 * qJD(3) - t144 * t192 - t146 * t191;
t195 = t270 * t233;
t56 = -t146 * t192 - qJD(3) * t195 + (qJD(3) * t208 + t191) * t144;
t218 = t266 * qJD(3);
t221 = t268 * qJD(3);
t82 = -qJD(2) * t233 - t146 * t221 + (t266 * qJD(2) + t218) * t144;
t90 = t144 * t208 - t195;
t288 = 0.2e1 * t55 * t111 + 0.2e1 * t161 * t90 - 0.2e1 * t178 * t56 - 0.2e1 * t89 * t82;
t287 = 0.2e1 * t82 * t111 - 0.2e1 * t161 * t178;
t267 = cos(qJ(6));
t219 = t267 * qJD(6);
t286 = t267 * qJD(5) + t219;
t136 = pkin(2) * t218;
t143 = sin(qJ(5));
t145 = cos(qJ(5));
t265 = sin(qJ(6));
t229 = t265 * t145;
t110 = t267 * t143 + t229;
t112 = -t265 * t143 + t267 * t145;
t274 = qJD(5) + qJD(6);
t80 = t274 * t110;
t222 = qJD(6) * t265;
t277 = t265 * qJD(5) + t222;
t81 = -t277 * t143 + t145 * t286;
t29 = ((t267 * t110 - t265 * t112) * qJD(6) + t265 * t81 - t267 * t80) * pkin(5);
t174 = -pkin(4) * t178 + t89;
t193 = t82 * qJ(4) + qJD(4) * t178;
t243 = t144 * qJD(2);
t237 = pkin(2) * t243;
t285 = -qJD(5) * t174 + t271 * t161 - t193 - t237;
t134 = -pkin(2) * t146 - pkin(1);
t188 = qJ(4) * t178 + t134;
t170 = -t271 * t111 - t188;
t284 = -t82 * pkin(4) + qJD(5) * t170 + t56;
t141 = t143 ^ 2;
t142 = t145 ^ 2;
t247 = t141 - t142;
t275 = t247 * qJD(5);
t254 = t112 * t80;
t255 = t110 * t81;
t18 = -0.2e1 * t254 + 0.2e1 * t255;
t213 = pkin(2) * t221;
t122 = t213 + qJD(4);
t116 = t122 * t145;
t238 = t266 * pkin(2);
t128 = t238 + qJ(4);
t245 = qJD(5) * t143;
t283 = -t128 * t245 + t116;
t282 = t141 + t142;
t158 = t145 * t161;
t281 = -t111 * t245 - t158;
t159 = t143 * t161;
t244 = qJD(5) * t145;
t280 = t111 * t244 - t159;
t239 = t268 * pkin(2);
t133 = -t239 - pkin(3);
t199 = pkin(9) - t133;
t183 = qJD(5) * t199;
t197 = t143 * t136;
t279 = t145 * t183 - t197;
t196 = t145 * t136;
t278 = t143 * t183 + t196;
t273 = t111 * t286 - t161 * t265;
t250 = t143 * t111;
t33 = t143 * t170 + t145 * t174;
t153 = -pkin(5) * t178 - pkin(10) * t250 + t33;
t152 = t267 * t153;
t252 = t111 * t145;
t34 = t143 * t174 - t145 * t170;
t26 = pkin(10) * t252 + t34;
t14 = -t265 * t26 + t152;
t151 = t265 * t153;
t15 = t267 * t26 + t151;
t11 = t285 * t143 + t284 * t145;
t269 = t82 * pkin(5);
t148 = -t280 * pkin(10) + t11 - t269;
t10 = -t284 * t143 + t285 * t145;
t150 = t281 * pkin(10) - t10;
t3 = -qJD(6) * t152 - t265 * t148 - t267 * t150 + t26 * t222;
t4 = -qJD(6) * t151 + t267 * t148 - t265 * t150 - t26 * t219;
t223 = t3 * t110 - t4 * t112 + t14 * t80 - t15 * t81;
t189 = pkin(10) + t199;
t105 = t189 * t143;
t164 = t189 * t245 + t196;
t182 = t145 * t189;
t165 = -qJD(5) * t182 + t197;
t169 = t267 * t182;
t31 = qJD(6) * t169 - t105 * t222 - t265 * t164 - t267 * t165;
t168 = t265 * t182;
t32 = qJD(6) * t168 + t105 * t219 + t267 * t164 - t265 * t165;
t72 = t265 * t105 - t169;
t73 = -t267 * t105 - t168;
t216 = t31 * t110 - t32 * t112 + t72 * t80 - t73 * t81;
t241 = pkin(10) + t271;
t118 = t241 * t143;
t214 = t241 * t145;
t187 = t267 * t214;
t200 = t241 * t245;
t51 = -t118 * t222 + t274 * t187 - t265 * t200;
t186 = t265 * t214;
t52 = t118 * t219 + t274 * t186 + t267 * t200;
t85 = t265 * t118 - t187;
t86 = -t267 * t118 - t186;
t215 = t51 * t110 - t52 * t112 + t85 * t80 - t86 * t81;
t98 = t282 * t136;
t147 = 0.2e1 * qJD(4);
t264 = t144 * pkin(2);
t38 = pkin(4) * t161 - t55;
t23 = -t281 * pkin(5) + t38;
t57 = (-pkin(5) * t145 - pkin(4)) * t111 + t90;
t263 = t23 * t110 + t57 * t81;
t262 = t23 * t112 - t57 * t80;
t69 = -t111 * pkin(4) + t90;
t261 = t38 * t143 + t69 * t244;
t135 = pkin(5) * t244;
t107 = t135 + t122;
t139 = t143 * pkin(5);
t119 = t139 + t128;
t260 = t107 * t110 + t119 * t81;
t259 = t107 * t112 - t119 * t80;
t123 = qJD(4) + t135;
t129 = qJ(4) + t139;
t258 = t123 * t110 + t129 * t81;
t257 = t123 * t112 - t129 * t80;
t251 = t128 * t122;
t117 = t128 * t244;
t249 = t122 * t143 + t117;
t137 = qJD(4) * t143;
t248 = qJ(4) * t244 + t137;
t138 = qJD(4) * t145;
t242 = t146 * qJD(2);
t240 = -0.2e1 * pkin(1) * qJD(2);
t65 = 0.2e1 * t178 * t82;
t236 = t143 * t271;
t235 = t145 * t271;
t234 = qJD(5) * t271;
t225 = t144 * t242;
t224 = t143 * t244;
t212 = pkin(5) * t219;
t211 = pkin(5) * t222;
t207 = t178 * t234;
t108 = t111 ^ 2;
t205 = t108 * t224;
t202 = t38 * t145 - t69 * t245;
t201 = -t55 * t90 + t56 * t89;
t156 = t161 * t267;
t22 = t274 * t111 * t229 + t273 * t143 + t145 * t156;
t70 = t112 * t111;
t17 = t110 * t22 - t70 * t81;
t21 = t143 * t156 - t273 * t145 + t277 * t250;
t71 = t110 * t111;
t16 = -t112 * t21 - t71 * t80;
t42 = t110 * t82 + t178 * t81;
t190 = t122 * qJ(4) + t128 * qJD(4);
t185 = -t143 * t82 - t178 * t244;
t59 = -t145 * t82 + t178 * t245;
t184 = t82 * t199;
t157 = t143 * t158;
t140 = qJ(4) * t147;
t121 = -0.2e1 * t224;
t120 = 0.2e1 * t224;
t109 = 0.2e1 * t275;
t79 = pkin(3) * t111 + t188;
t64 = -0.2e1 * t254;
t63 = 0.2e1 * t255;
t48 = -t111 * t275 - t157;
t45 = -pkin(3) * t276 + (-t178 * pkin(3) + t264) * qJD(2) + t193;
t41 = -t112 * t82 + t178 * t80;
t40 = t224 * t290 + (qJD(2) + qJD(3)) * t247 * t178;
t39 = 0.2e1 * t110 * t80 - 0.2e1 * t112 * t81;
t6 = t110 * t21 - t112 * t22 - t70 * t80 - t71 * t81;
t5 = -t10 * t143 - t33 * t245 + (t34 * qJD(5) + t11) * t145;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t225, 0.2e1 * (-t144 ^ 2 + t146 ^ 2) * qJD(2), 0, -0.2e1 * t225, 0, 0, t144 * t240, t146 * t240, 0, 0, t65, t287, 0, t289, 0, 0, -0.2e1 * t134 * t276 + 0.2e1 * (t111 * t264 - t134 * t178) * qJD(2), -0.2e1 * t134 * t82 - 0.2e1 * t178 * t237, t288, 0.2e1 * t134 * t237 + 0.2e1 * t201, 0, 0, 0, t65, t287, t289, t288, -0.2e1 * t45 * t111 + 0.2e1 * t161 * t79, 0.2e1 * t178 * t45 + 0.2e1 * t79 * t82, 0.2e1 * t45 * t79 + 0.2e1 * t201, -0.2e1 * t141 * t160 + 0.2e1 * t205, -0.2e1 * t108 * t275 + t157 * t290, 0.2e1 * t111 * t185 + 0.2e1 * t159 * t178, -0.2e1 * t142 * t160 - 0.2e1 * t205, 0.2e1 * t111 * t59 + 0.2e1 * t158 * t178, t65, -0.2e1 * t11 * t178 - 0.2e1 * t38 * t252 - 0.2e1 * t281 * t69 - 0.2e1 * t33 * t82, -0.2e1 * t10 * t178 + 0.2e1 * t250 * t38 + 0.2e1 * t280 * t69 + 0.2e1 * t34 * t82, -0.2e1 * t10 * t252 - 0.2e1 * t11 * t250 - 0.2e1 * t280 * t33 + 0.2e1 * t281 * t34, -0.2e1 * t10 * t34 + 0.2e1 * t11 * t33 + 0.2e1 * t38 * t69, -0.2e1 * t71 * t21, -0.2e1 * t21 * t70 - 0.2e1 * t22 * t71, 0.2e1 * t178 * t21 - 0.2e1 * t71 * t82, -0.2e1 * t70 * t22, 0.2e1 * t178 * t22 - 0.2e1 * t70 * t82, t65, -0.2e1 * t14 * t82 - 0.2e1 * t178 * t4 + 0.2e1 * t22 * t57 - 0.2e1 * t23 * t70, 0.2e1 * t15 * t82 - 0.2e1 * t178 * t3 - 0.2e1 * t21 * t57 + 0.2e1 * t23 * t71, 0.2e1 * t14 * t21 - 0.2e1 * t15 * t22 - 0.2e1 * t3 * t70 - 0.2e1 * t4 * t71, 0.2e1 * t14 * t4 - 0.2e1 * t15 * t3 + 0.2e1 * t23 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, 0, -t243, 0, -pkin(7) * t242, pkin(7) * t243, 0, 0, 0, 0, -t82, 0, t161, 0, -t56, t55, -t111 * t213 - t136 * t178 + t161 * t238 + t82 * t239 (-t266 * t55 - t268 * t56 + (t266 * t89 + t268 * t90) * qJD(3)) * pkin(2), 0, t82, -t161, 0, 0, 0, -t122 * t111 + t128 * t173 - t133 * t82 + (t128 - t238) * t276, t56, -t55, t90 * t122 - t55 * t128 + t56 * t133 + t89 * t136, t48, t40, t59, -t48, -t185, 0, -t283 * t111 + t128 * t158 + t145 * t184 - t178 * t278 + t261, t111 * t117 + t122 * t250 - t128 * t159 - t143 * t184 - t178 * t279 + t202, -t5, t69 * t122 + t38 * t128 + (-t11 * t199 + t136 * t33 - t183 * t34) * t145 + (t10 * t199 + t136 * t34 + t183 * t33) * t143, t16, t6, t41, t17, t42, 0, -t107 * t70 + t119 * t22 - t178 * t32 - t72 * t82 + t263, t107 * t71 - t119 * t21 - t178 * t31 + t73 * t82 + t262, t21 * t72 - t22 * t73 - t31 * t70 - t32 * t71 + t223, t107 * t57 + t119 * t23 + t14 * t32 - t15 * t31 - t3 * t73 + t4 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t136, -0.2e1 * t213, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t136, 0.2e1 * t122, 0.2e1 * t133 * t136 + 0.2e1 * t251, t121, t109, 0, t120, 0, 0, 0.2e1 * t249, 0.2e1 * t283, -0.2e1 * t98, -0.2e1 * t282 * qJD(3) * t199 * t238 + 0.2e1 * t251, t64, t39, 0, t63, 0, 0, 0.2e1 * t260, 0.2e1 * t259, 0.2e1 * t216, 0.2e1 * t107 * t119 - 0.2e1 * t31 * t73 + 0.2e1 * t32 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, 0, t161, 0, -t56, t55, 0, 0, 0, t82, -t161, 0, 0, 0, pkin(3) * t82 + qJ(4) * t161 - qJD(4) * t111, t56, -t55, -pkin(3) * t56 - qJ(4) * t55 + qJD(4) * t90, t48, t40, t59, -t48, -t185, 0, -t281 * qJ(4) - t111 * t138 - t143 * t207 + t235 * t82 + t261, t280 * qJ(4) + t111 * t137 - t145 * t207 - t236 * t82 + t202, -t5, t10 * t236 - t11 * t235 + t38 * qJ(4) + t69 * qJD(4) + (-t235 * t34 + t236 * t33) * qJD(5), t16, t6, t41, t17, t42, 0, -t123 * t70 + t129 * t22 - t178 * t52 - t82 * t85 + t263, t123 * t71 - t129 * t21 - t178 * t51 + t82 * t86 + t262, t21 * t85 - t22 * t86 - t51 * t70 - t52 * t71 + t223, t123 * t57 + t129 * t23 + t14 * t52 - t15 * t51 - t3 * t86 + t4 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, -t213, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t147 + t213, -pkin(3) * t136 + t190, t121, t109, 0, t120, 0, 0, t248 + t249, t116 + t138 + (-qJ(4) - t128) * t245, -t98, -t271 * t98 + t190, t64, t39, 0, t63, 0, 0, t258 + t260, t257 + t259, t215 + t216, t107 * t129 + t119 * t123 - t31 * t86 + t32 * t85 - t51 * t73 + t52 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t140, t121, t109, 0, t120, 0, 0, 0.2e1 * t248, -0.2e1 * qJ(4) * t245 + 0.2e1 * t138, 0, t140, t64, t39, 0, t63, 0, 0, 0.2e1 * t258, 0.2e1 * t257, 0.2e1 * t215, 0.2e1 * t123 * t129 - 0.2e1 * t51 * t86 + 0.2e1 * t52 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, 0, 0, t56, 0, 0, 0, 0, 0, 0, t59, -t185, 0, t5, 0, 0, 0, 0, 0, 0, t41, t42, -t16 - t17, -t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t280, 0, t281, -t82, t11, t10, 0, 0, 0, 0, -t21, 0, -t22, -t82, t178 * t211 - t267 * t269 + t4 (t178 * t219 + t265 * t82) * pkin(5) + t3 (t267 * t21 - t265 * t22 + (t265 * t71 + t267 * t70) * qJD(6)) * pkin(5) (t267 * t4 - t265 * t3 + (-t265 * t14 + t267 * t15) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, 0, -t244, 0, t278, t279, 0, 0, 0, 0, -t80, 0, -t81, 0, t32, t31, -t29 (t267 * t32 - t265 * t31 + (-t265 * t72 + t267 * t73) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, 0, -t244, 0, t143 * t234, t145 * t234, 0, 0, 0, 0, -t80, 0, -t81, 0, t52, t51, -t29 (t267 * t52 - t265 * t51 + (-t265 * t85 + t267 * t86) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, -t244, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t81, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t211, -0.2e1 * t212, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, -t22, -t82, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, 0, -t81, 0, t32, t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, 0, -t81, 0, t52, t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t81, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, -t212, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
