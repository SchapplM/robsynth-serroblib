% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP11_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_inertiaDJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:41:13
% EndTime: 2019-03-09 06:41:37
% DurationCPUTime: 9.51s
% Computational Cost: add. (15428->500), mult. (45156->874), div. (0->0), fcn. (47933->12), ass. (0->261)
t141 = sin(pkin(6));
t139 = sin(pkin(12));
t142 = cos(pkin(12));
t287 = cos(pkin(7));
t305 = sin(qJ(3));
t216 = t287 * t305;
t306 = cos(qJ(3));
t187 = t306 * t139 + t142 * t216;
t140 = sin(pkin(7));
t288 = cos(pkin(6));
t231 = t140 * t288;
t202 = t305 * t231;
t314 = t187 * t141 + t202;
t325 = t314 * qJD(3);
t230 = t141 * t287;
t316 = t142 * t230 + t231;
t180 = t316 * t306;
t283 = t139 * t141;
t225 = t305 * t283;
t83 = t225 - t180;
t332 = 0.2e1 * t325 * t83;
t252 = pkin(1) * t288;
t189 = t288 * pkin(2) + t142 * t252;
t178 = (-t287 * pkin(9) - qJ(2)) * t283 + t189;
t253 = -pkin(2) * t142 - pkin(1);
t284 = t139 * t140;
t191 = (-pkin(9) * t284 + t253) * t141;
t331 = t140 * t191 + t178 * t287;
t144 = sin(qJ(4));
t133 = qJD(4) * t144;
t146 = cos(qJ(4));
t330 = -t83 * t133 + t146 * t325;
t282 = t141 * t142;
t194 = t140 * t282 - t288 * t287;
t160 = -t144 * t194 + t146 * t314;
t171 = pkin(9) * t316 + qJ(2) * t282 + t139 * t252;
t315 = t305 * t171 - t306 * t331;
t60 = pkin(3) * t194 + t315;
t73 = t144 * t314 + t146 * t194;
t151 = t73 * pkin(4) - pkin(11) * t160 + t60;
t62 = t306 * t171 + t305 * t331;
t154 = -pkin(10) * t194 + t62;
t273 = qJD(2) * t141;
t223 = t273 * t284;
t115 = qJD(3) * t225;
t312 = -qJD(3) * t180 + t115;
t327 = -pkin(3) * t325 - pkin(10) * t312 + qJD(4) * t154 - t223;
t201 = qJD(2) * t139 * t230;
t248 = t142 * t273;
t152 = qJD(3) * t315 + t305 * t201 - t306 * t248;
t308 = t83 * pkin(3);
t156 = -pkin(10) * t314 - t140 * t178 + t287 * t191 + t308;
t328 = -qJD(4) * t156 + t152;
t22 = t327 * t144 + t328 * t146;
t329 = -pkin(11) * t325 - qJD(5) * t151 + t22;
t250 = t140 * t305;
t102 = t144 * t287 + t146 * t250;
t235 = qJD(3) * t306;
t218 = t146 * t235;
t229 = t146 * t287;
t249 = t144 * t305;
t101 = t140 * t249 - t229;
t265 = t101 * qJD(4);
t175 = t140 * t218 - t265;
t220 = t140 * t235;
t272 = qJD(4) * t102;
t85 = t144 * t220 + t272;
t197 = t85 * t144 + t146 * t265;
t326 = -t102 * t133 + t175 * t146 + t197;
t323 = -0.4e1 * t144;
t322 = 0.2e1 * t194;
t321 = qJD(3) * t62;
t143 = sin(qJ(5));
t135 = t143 ^ 2;
t145 = cos(qJ(5));
t137 = t145 ^ 2;
t275 = t135 - t137;
t320 = qJD(5) * t160 - t325;
t269 = qJD(4) * t146;
t319 = t144 * t325 + t83 * t269;
t190 = t194 * qJD(4);
t317 = -t312 - t190;
t169 = qJD(4) * t314;
t63 = t144 * t317 + t146 * t169;
t51 = t73 * t133 - t146 * t63;
t313 = t144 * qJD(6) + (pkin(10) * qJD(5) + qJ(6) * qJD(4)) * t146;
t138 = t146 ^ 2;
t311 = 0.2e1 * qJD(5);
t310 = t63 * pkin(5);
t65 = t143 * t160 - t83 * t145;
t309 = t65 * pkin(5);
t307 = t83 * pkin(4);
t304 = pkin(5) * t145;
t303 = pkin(10) * t143;
t23 = t328 * t144 - t327 * t146;
t19 = -pkin(4) * t325 - t23;
t79 = t144 * t169;
t159 = t146 * t317 - t79;
t268 = qJD(5) * t143;
t44 = t143 * t159 + t145 * t320 + t83 * t268;
t7 = t44 * pkin(5) + t19;
t302 = t143 * t7;
t301 = t145 * t7;
t300 = -qJ(6) - pkin(11);
t38 = t144 * t156 + t146 * t154;
t36 = pkin(11) * t83 + t38;
t17 = t143 * t151 + t145 * t36;
t299 = t101 * t85;
t132 = qJD(5) * t145;
t45 = t83 * t132 - t143 * t320 + t145 * t159;
t298 = t143 * t45;
t297 = t143 * t65;
t66 = t83 * t143 + t145 * t160;
t296 = t143 * t66;
t295 = t145 * t44;
t294 = t145 * t45;
t293 = t145 * t66;
t243 = t144 * t132;
t105 = t143 * t269 + t243;
t259 = pkin(10) * t269;
t95 = pkin(5) * t105 + t259;
t290 = t95 * t143;
t289 = t95 * t145;
t119 = t300 * t143;
t286 = t119 * t144;
t120 = t300 * t145;
t285 = t120 * t144;
t281 = t143 * t144;
t280 = t143 * t146;
t279 = t144 * t145;
t278 = t145 * t146;
t215 = -pkin(4) * t146 - t144 * pkin(11);
t203 = -pkin(3) + t215;
t193 = qJD(5) * t203;
t214 = pkin(4) * t144 - pkin(11) * t146;
t196 = qJD(4) * t214;
t277 = -t143 * t196 - t145 * t193;
t238 = t143 * t133;
t276 = pkin(10) * t238 + t145 * t196;
t128 = pkin(10) * t278;
t97 = t143 * t203 + t128;
t136 = t144 ^ 2;
t274 = t136 - t138;
t271 = qJD(4) * t143;
t270 = qJD(4) * t145;
t267 = qJD(5) * t144;
t266 = qJD(5) * t146;
t52 = 0.2e1 * t73 * t63;
t263 = -0.2e1 * pkin(3) * qJD(4);
t262 = pkin(10) * t280;
t261 = -0.2e1 * t268;
t260 = pkin(5) * t133;
t258 = pkin(5) * t268;
t58 = t144 * t154;
t26 = -t146 * t156 - t307 + t309 + t58;
t255 = t26 * t268;
t251 = t140 * t306;
t247 = t145 * t269;
t246 = t101 * t268;
t245 = t143 * t267;
t244 = t143 * t266;
t242 = t145 * t266;
t114 = (pkin(5) * t143 + pkin(10)) * t144;
t239 = t114 * t268;
t237 = t143 * t132;
t236 = t144 * t269;
t234 = qJD(3) * t305;
t233 = qJD(4) * t306;
t130 = -pkin(4) - t304;
t232 = -t130 + t304;
t16 = -t143 * t36 + t145 * t151;
t228 = t274 * qJD(4);
t227 = 0.2e1 * t236;
t54 = t306 * t201 + t305 * t248 + t321;
t226 = t54 * t251;
t224 = t143 * t251;
t222 = t143 * t247;
t221 = t136 * t237;
t219 = t140 * t234;
t13 = -qJ(6) * t65 + t17;
t8 = pkin(5) * t73 - qJ(6) * t66 + t16;
t213 = -t13 * t143 - t145 * t8;
t212 = pkin(5) * t135 + t130 * t145;
t211 = -t143 * t17 - t145 * t16;
t210 = -t145 * t65 - t296;
t192 = t143 * t102 + t145 * t251;
t87 = t145 * t102 - t224;
t209 = -t143 * t87 + t145 * t192;
t111 = t145 * t203;
t96 = t111 - t262;
t208 = -t143 * t97 - t145 * t96;
t207 = -t23 * t144 - t22 * t146;
t206 = -0.2e1 * t288 * t273;
t185 = qJ(2) * t284 + t287 * t253;
t186 = t140 * t189;
t37 = -t58 + (-pkin(10) * t202 - t186 + t308 + (-pkin(10) * t187 + t185) * t141) * t146;
t35 = -t37 - t307;
t200 = t35 * t132 + t143 * t19;
t199 = -t145 * t19 + t35 * t268;
t50 = t73 * t132 + t143 * t63;
t198 = -t145 * t63 + t73 * t268;
t74 = t101 * t132 + t143 * t85;
t75 = -t145 * t85 + t246;
t149 = t63 * pkin(4) - pkin(11) * t159 + t54;
t3 = -t143 * t149 + t329 * t145 + t36 * t268;
t195 = 0.2e1 * (t139 ^ 2 + t142 ^ 2) * t141 ^ 2 * qJD(2);
t103 = -t245 + t247;
t104 = t144 * t270 + t244;
t4 = -t36 * t132 + t329 * t143 + t145 * t149;
t182 = t211 * qJD(5) - t143 * t4 - t145 * t3;
t68 = t192 * qJD(5) - t143 * t219 - t145 * t175;
t69 = qJD(5) * t224 - t102 * t132 - t143 * t175 + t145 * t219;
t46 = t209 * qJD(5) - t143 * t69 - t145 * t68;
t76 = pkin(10) * t104 + t277;
t77 = -t97 * qJD(5) + t276;
t181 = t208 * qJD(5) - t143 * t77 - t145 * t76;
t179 = t180 * t83;
t161 = qJ(6) * t245 - t143 * t193 - t145 * t313 + t276;
t147 = -t45 * qJ(6) - t66 * qJD(6) + t4;
t126 = -0.2e1 * t236;
t125 = -0.2e1 * t237;
t124 = 0.2e1 * t237;
t112 = -0.2e1 * t275 * qJD(5);
t106 = t238 - t242;
t100 = -qJD(6) * t143 + t300 * t132;
t99 = -t145 * qJD(6) - t300 * t268;
t94 = 0.2e1 * t137 * t236 - 0.2e1 * t221;
t93 = 0.2e1 * t135 * t236 + 0.2e1 * t221;
t92 = t267 * t275 - t222;
t91 = t237 * t323 - t275 * t269;
t90 = -0.2e1 * t143 * t228 + 0.2e1 * t144 * t242;
t89 = 0.2e1 * t144 * t244 + 0.2e1 * t274 * t270;
t88 = -qJ(6) * t281 + t97;
t84 = t275 * t136 * t311 + t222 * t323;
t82 = -qJ(6) * t279 + t111 + (-pkin(5) - t303) * t146;
t71 = t141 * t185 - t186;
t70 = (pkin(10) * qJD(4) + qJ(6) * qJD(5)) * t279 + t313 * t143 + t277;
t67 = t161 + t260;
t48 = (t143 * t265 - t69) * t146 + (-qJD(4) * t192 + t74) * t144;
t47 = (t145 * t265 - t68) * t146 + (-qJD(4) * t87 - t75) * t144;
t39 = -0.2e1 * t192 * t69 - 0.2e1 * t68 * t87 + 0.2e1 * t299;
t33 = t209 * t269 + (t143 * t68 - t145 * t69 + (-t143 * t192 - t145 * t87) * qJD(5)) * t144;
t32 = 0.2e1 * t66 * t45;
t31 = 0.2e1 * t65 * t44;
t28 = t66 * t132 + t298;
t27 = t65 * t268 - t295;
t25 = t66 * t247 + (-t66 * t268 + t294) * t144;
t24 = t65 * t243 + (t144 * t44 + t65 * t269) * t143;
t21 = 0.2e1 * t45 * t73 + 0.2e1 * t63 * t66;
t20 = -0.2e1 * t44 * t73 - 0.2e1 * t63 * t65;
t15 = (-t73 * t271 + t44) * t146 + (-qJD(4) * t65 - t50) * t144;
t14 = (t73 * t270 - t45) * t146 + (qJD(4) * t66 - t198) * t144;
t12 = t101 * t44 - t192 * t63 + t65 * t85 + t69 * t73;
t11 = t101 * t45 - t63 * t87 + t66 * t85 + t68 * t73;
t10 = -0.2e1 * t44 * t66 - 0.2e1 * t45 * t65;
t9 = t210 * qJD(5) - t143 * t44 + t294;
t6 = t192 * t45 - t44 * t87 + t65 * t68 - t66 * t69;
t5 = t210 * t269 + (-t298 - t295 + (-t293 + t297) * qJD(5)) * t144;
t2 = t44 * qJ(6) + t65 * qJD(6) + t3;
t1 = t147 + t310;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139 * t206, t142 * t206, t195, qJ(2) * t195, -0.2e1 * t314 * t312, 0.2e1 * t115 * t83 + 0.2e1 * (-t314 ^ 2 - t179) * qJD(3), t312 * t322, t332, t325 * t322, 0, 0.2e1 * t194 * t54 + 0.2e1 * t83 * t223 + 0.2e1 * t325 * t71, -0.2e1 * t152 * t194 + 0.2e1 * t223 * t314 - 0.2e1 * t312 * t71, 0.2e1 * t152 * t83 - 0.2e1 * t312 * t315 + 0.2e1 * t314 * t54 - 0.2e1 * t325 * t62, 0.2e1 * (t62 * (-t139 * t216 + t306 * t142) + t71 * t284) * t273 - 0.2e1 * (-t54 + t321) * t315, 0.2e1 * t160 * t159, -0.2e1 * t159 * t73 - 0.2e1 * t160 * t63, 0.2e1 * (-t79 + (-t115 - t190) * t146) * t83 + 0.2e1 * (t146 * t179 + t160 * t314) * qJD(3), t52, -0.2e1 * t325 * t73 - 0.2e1 * t63 * t83, t332, 0.2e1 * t23 * t83 + 0.2e1 * t325 * t37 + 0.2e1 * t54 * t73 + 0.2e1 * t60 * t63, 0.2e1 * t159 * t60 + 0.2e1 * t160 * t54 + 0.2e1 * t22 * t83 - 0.2e1 * t325 * t38, -0.2e1 * t159 * t37 - 0.2e1 * t160 * t23 + 0.2e1 * t22 * t73 - 0.2e1 * t38 * t63, -0.2e1 * t22 * t38 + 0.2e1 * t23 * t37 + 0.2e1 * t54 * t60, t32, t10, t21, t31, t20, t52, 0.2e1 * t16 * t63 + 0.2e1 * t19 * t65 + 0.2e1 * t35 * t44 + 0.2e1 * t4 * t73, -0.2e1 * t17 * t63 + 0.2e1 * t19 * t66 + 0.2e1 * t3 * t73 + 0.2e1 * t35 * t45, -0.2e1 * t16 * t45 - 0.2e1 * t17 * t44 + 0.2e1 * t3 * t65 - 0.2e1 * t4 * t66, 0.2e1 * t16 * t4 - 0.2e1 * t17 * t3 + 0.2e1 * t19 * t35, t32, t10, t21, t31, t20, t52, 0.2e1 * t1 * t73 + 0.2e1 * t26 * t44 + 0.2e1 * t63 * t8 + 0.2e1 * t65 * t7, -0.2e1 * t13 * t63 + 0.2e1 * t2 * t73 + 0.2e1 * t26 * t45 + 0.2e1 * t66 * t7, -0.2e1 * t1 * t66 - 0.2e1 * t13 * t44 + 0.2e1 * t2 * t65 - 0.2e1 * t45 * t8, 0.2e1 * t1 * t8 - 0.2e1 * t13 * t2 + 0.2e1 * t26 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t194 * t250 + t287 * t314) * qJD(3), t194 * t220 - t287 * t312, -t83 * t220 + t251 * t312, t140 * t201 - t152 * t250 + t219 * t315 + t62 * t220 - t226, 0, 0, 0, 0, 0, 0, -t101 * t325 + t73 * t219 - t63 * t251 - t85 * t83, -t102 * t325 - t159 * t251 + t160 * t219 - t175 * t83, t101 * t159 - t102 * t63 + t160 * t85 - t175 * t73, -t23 * t101 - t22 * t102 + t175 * t38 + t219 * t60 - t37 * t85 - t226, 0, 0, 0, 0, 0, 0, t12, t11, t6, t101 * t19 + t16 * t69 - t17 * t68 - t192 * t4 - t3 * t87 + t35 * t85, 0, 0, 0, 0, 0, 0, t12, t11, t6, -t1 * t192 + t101 * t7 - t13 * t68 - t2 * t87 + t26 * t85 + t69 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t229 * t272 + 0.2e1 * t299 + 0.2e1 * (t102 * (-qJD(4) * t249 + t218) - t306 * t219) * t140, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t312, 0, -t325, 0, -t54, t152, 0, 0, t138 * t169 + (-t79 + (-0.2e1 * t190 - t312) * t146) * t144, -t160 * t133 - t144 * t63 + t146 * t159 - t73 * t269, t319, t51, t330, 0, -pkin(3) * t63 - pkin(10) * t319 + t60 * t133 - t54 * t146, -pkin(3) * t159 - pkin(10) * t330 + t54 * t144 + t60 * t269, -t38 * t133 + t160 * t259 - t37 * t269 + t207 + (t144 * t159 + t51) * pkin(10), -t54 * pkin(3) + ((-t38 * t144 - t37 * t146) * qJD(4) + t207) * pkin(10), t25, t5, t14, t24, t15, t51, t96 * t63 + t77 * t73 + (-t4 + (pkin(10) * t65 + t143 * t35) * qJD(4)) * t146 + (pkin(10) * t44 + qJD(4) * t16 + t200) * t144, -t97 * t63 + t76 * t73 + (-t3 + (pkin(10) * t66 + t145 * t35) * qJD(4)) * t146 + (pkin(10) * t45 - qJD(4) * t17 - t199) * t144, -t97 * t44 - t96 * t45 + t76 * t65 - t77 * t66 + t211 * t269 + (t143 * t3 - t145 * t4 + (t143 * t16 - t145 * t17) * qJD(5)) * t144, t16 * t77 - t17 * t76 - t3 * t97 + t4 * t96 + (t144 * t19 + t269 * t35) * pkin(10), t25, t5, t14, t24, t15, t51, t114 * t44 + t82 * t63 + t95 * t65 + t67 * t73 + (t26 * t271 - t1) * t146 + (qJD(4) * t8 + t132 * t26 + t302) * t144, t114 * t45 - t88 * t63 + t95 * t66 + t70 * t73 + (t26 * t270 - t2) * t146 + (-qJD(4) * t13 - t255 + t301) * t144, -t88 * t44 - t82 * t45 + t70 * t65 - t67 * t66 + t213 * t269 + (-t1 * t145 + t143 * t2 + (-t13 * t145 + t143 * t8) * qJD(5)) * t144, t1 * t82 + t114 * t7 - t13 * t70 - t2 * t88 + t26 * t95 + t67 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, -t220, 0, 0, 0, 0, 0, 0, 0, 0 (-t144 * t233 - t146 * t234) * t140 (t144 * t234 - t146 * t233) * t140, t326, -pkin(3) * t219 + t326 * pkin(10), 0, 0, 0, 0, 0, 0, t48, t47, t33, pkin(10) * t197 - t192 * t77 - t68 * t97 + t69 * t96 - t87 * t76, 0, 0, 0, 0, 0, 0, t48, t47, t33, t101 * t95 + t114 * t85 - t192 * t67 - t68 * t88 + t69 * t82 - t70 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, -0.2e1 * t228, 0, t126, 0, 0, t144 * t263, t146 * t263, 0, 0, t94, t84, t89, t93, t90, t126, 0.2e1 * t96 * t133 - 0.2e1 * t146 * t77 + 0.2e1 * (t132 * t136 + t143 * t227) * pkin(10), -0.2e1 * t97 * t133 - 0.2e1 * t146 * t76 + 0.2e1 * (-t136 * t268 + t145 * t227) * pkin(10), 0.2e1 * t208 * t269 + 0.2e1 * (t143 * t76 - t145 * t77 + (t143 * t96 - t145 * t97) * qJD(5)) * t144, 0.2e1 * pkin(10) ^ 2 * t236 - 0.2e1 * t97 * t76 + 0.2e1 * t96 * t77, t94, t84, t89, t93, t90, t126, 0.2e1 * (t114 * t271 - t67) * t146 + 0.2e1 * (qJD(4) * t82 + t114 * t132 + t290) * t144, 0.2e1 * (t114 * t270 - t70) * t146 + 0.2e1 * (-qJD(4) * t88 - t239 + t289) * t144, 0.2e1 * (-t143 * t88 - t145 * t82) * t269 + 0.2e1 * (t143 * t70 - t145 * t67 + (t143 * t82 - t145 * t88) * qJD(5)) * t144, 0.2e1 * t114 * t95 + 0.2e1 * t67 * t82 - 0.2e1 * t70 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, 0, -t63, t325, t23, t22, 0, 0, t28, t9, t50, t27, -t198, 0, -pkin(4) * t44 - pkin(11) * t50 + t199, -pkin(4) * t45 + pkin(11) * t198 + t200 (t298 - t295 + (t293 + t297) * qJD(5)) * pkin(11) + t182, -pkin(4) * t19 + pkin(11) * t182, t28, t9, t50, t27, -t198, 0, t100 * t73 + t119 * t63 + t130 * t44 - t301 + (t26 + t309) * t268, t120 * t63 + t130 * t45 + t302 + t73 * t99 + (pkin(5) * t296 + t145 * t26) * qJD(5), qJD(5) * t213 - t1 * t143 - t100 * t66 - t119 * t45 + t120 * t44 - t145 * t2 + t65 * t99, pkin(5) * t255 + t1 * t119 + t100 * t8 + t120 * t2 - t13 * t99 + t130 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t175, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, t46, -pkin(4) * t85 + pkin(11) * t46, 0, 0, 0, 0, 0, 0, t75, t74, t46, pkin(5) * t246 - t100 * t192 + t119 * t69 + t120 * t68 + t130 * t85 - t87 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t269, 0, -t133, 0, -t259, pkin(10) * t133, 0, 0, -t92, t91, t106, t92, t104, 0 (pkin(11) * t278 + (-pkin(4) * t145 + t303) * t144) * qJD(5) + (t143 * t215 - t128) * qJD(4) (pkin(10) * t279 + t143 * t214) * qJD(5) + (t145 * t215 + t262) * qJD(4), t181, -pkin(4) * t259 + pkin(11) * t181, -t92, t91, t106, t92, t104, 0, -t100 * t146 - t289 + (t130 * t280 + t286) * qJD(4) + (t114 * t143 + t144 * t212) * qJD(5), t290 - t146 * t99 + (t130 * t278 + t285) * qJD(4) + (t114 * t145 + t232 * t281) * qJD(5) (-t119 * t269 - t100 * t144 - t70 + (-t82 + t285) * qJD(5)) * t145 + (t120 * t269 + t144 * t99 - t67 + (-t88 + t286) * qJD(5)) * t143, pkin(5) * t239 + t100 * t82 + t119 * t67 + t120 * t70 + t130 * t95 - t88 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t112, 0, t125, 0, 0, pkin(4) * t261, -0.2e1 * pkin(4) * t132, 0, 0, t124, t112, 0, t125, 0, 0, t232 * t261, t212 * t311, -0.2e1 * t100 * t143 - 0.2e1 * t145 * t99 + 0.2e1 * (-t119 * t145 + t120 * t143) * qJD(5), 0.2e1 * t100 * t119 + 0.2e1 * t120 * t99 + 0.2e1 * t130 * t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, -t44, t63, t4, t3, 0, 0, 0, 0, t45, 0, -t44, t63, t147 + 0.2e1 * t310, t2, -t45 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t68, 0, 0, 0, 0, 0, 0, 0, 0, t69, t68, 0, t69 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, -t105, t133, t77, t76, 0, 0, 0, 0, t103, 0, -t105, t133, t161 + 0.2e1 * t260, t70, -t103 * pkin(5), t67 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, 0, -t268, 0, -pkin(11) * t132, pkin(11) * t268, 0, 0, 0, 0, t132, 0, -t268, 0, t100, t99, -pkin(5) * t132, t100 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t103, 0, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, t132, 0, t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t18;