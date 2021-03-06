% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRP13_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:01:33
% EndTime: 2019-03-09 13:01:49
% DurationCPUTime: 6.27s
% Computational Cost: add. (5959->457), mult. (15304->783), div. (0->0), fcn. (13930->8), ass. (0->253)
t140 = sin(pkin(6));
t146 = cos(qJ(2));
t246 = qJD(2) * t146;
t125 = t140 * t246;
t142 = sin(qJ(4));
t145 = cos(qJ(4));
t259 = t140 * t146;
t268 = cos(pkin(6));
t181 = -t142 * t259 + t268 * t145;
t143 = sin(qJ(2));
t224 = pkin(1) * t268;
t201 = t143 * t224;
t288 = pkin(3) + pkin(8);
t231 = t140 * t288;
t167 = t146 * t231 + t201;
t206 = t268 * qJ(3);
t69 = t206 + t167;
t92 = t268 * t142 + t145 * t259;
t152 = t92 * pkin(4) - t181 * pkin(10) + t69;
t267 = qJ(3) * t143;
t289 = pkin(2) + pkin(9);
t169 = t140 * (-t289 * t146 - pkin(1) - t267);
t295 = qJD(2) * t167 - qJD(4) * t169;
t200 = t146 * t224;
t168 = -t268 * pkin(2) - t200;
t159 = -t268 * pkin(9) + t168;
t155 = t143 * t231 + t159;
t245 = qJD(3) * t143;
t266 = qJ(3) * t146;
t296 = -qJD(4) * t155 - t140 * (-t245 + (t289 * t143 - t266) * qJD(2));
t24 = -t142 * t295 + t145 * t296;
t297 = -pkin(10) * t125 - qJD(5) * t152 + t24;
t141 = sin(qJ(5));
t239 = qJD(5) * t145;
t218 = t141 * t239;
t144 = cos(qJ(5));
t242 = qJD(4) * t144;
t172 = t142 * t242 + t218;
t247 = qJD(2) * t143;
t214 = t140 * t247;
t205 = qJD(4) * t268;
t222 = qJD(4) * t259;
t252 = -t142 * t205 - t145 * t222;
t173 = t142 * t214 + t252;
t294 = t173 * t145;
t136 = t141 ^ 2;
t138 = t144 ^ 2;
t251 = t136 - t138;
t260 = t140 * t143;
t235 = pkin(8) * t260;
t170 = -t200 + t235;
t86 = t170 * qJD(2);
t293 = t143 ^ 2;
t292 = 0.2e1 * t140;
t291 = 0.2e1 * qJD(3);
t290 = 0.2e1 * qJD(5);
t64 = -t142 * t222 + (t205 - t214) * t145;
t287 = t64 * pkin(5);
t174 = t181 * t141;
t255 = t143 * t144;
t65 = -t140 * t255 + t174;
t286 = t65 * pkin(5);
t285 = pkin(5) * t144;
t284 = pkin(10) * t142;
t283 = t142 * pkin(4);
t25 = t142 * t296 + t145 * t295;
t21 = -pkin(4) * t125 - t25;
t66 = t141 * t260 + t181 * t144;
t37 = qJD(5) * t66 - t144 * t125 + t173 * t141;
t8 = t37 * pkin(5) + t21;
t282 = t8 * t141;
t281 = t8 * t144;
t280 = -qJ(6) - pkin(10);
t43 = t142 * t155 + t145 * t169;
t40 = pkin(10) * t260 + t43;
t23 = t141 * t152 + t144 * t40;
t279 = t141 * t65;
t278 = t141 * t66;
t277 = t142 * t69;
t101 = pkin(8) * t259 + t201;
t87 = t101 * qJD(2);
t276 = t143 * t87;
t275 = t144 * t65;
t274 = t144 * t66;
t273 = t37 * t144;
t132 = qJD(5) * t144;
t38 = t252 * t144 - qJD(5) * t174 + (t143 * t132 + (t141 * t146 + t142 * t255) * qJD(2)) * t140;
t272 = t38 * t141;
t271 = t64 * t142;
t241 = qJD(4) * t289;
t220 = t142 * t241;
t243 = qJD(4) * t142;
t213 = t141 * t243;
t216 = t144 * t239;
t99 = -t213 + t216;
t73 = t99 * pkin(5) - t220;
t270 = t73 * t141;
t269 = t73 * t144;
t265 = qJD(4) * t65;
t264 = qJD(4) * t66;
t114 = t280 * t141;
t263 = t114 * t145;
t115 = t280 * t144;
t262 = t115 * t145;
t130 = -pkin(4) - t285;
t261 = t130 * t144;
t258 = t141 * t145;
t257 = t141 * t289;
t256 = t142 * t289;
t254 = t144 * t145;
t253 = t145 * t289;
t121 = t144 * t256;
t193 = -pkin(10) * t145 + t283;
t178 = qJ(3) + t193;
t78 = t141 * t178 - t121;
t250 = t136 + t138;
t137 = t142 ^ 2;
t139 = t145 ^ 2;
t249 = t137 - t139;
t248 = t137 + t139;
t244 = qJD(4) * t141;
t133 = qJD(4) * t145;
t240 = qJD(5) * t141;
t238 = qJD(5) * t289;
t237 = t141 * qJD(6);
t236 = t144 * qJD(6);
t51 = 0.2e1 * t92 * t64;
t234 = -0.2e1 * t240;
t233 = pkin(5) * t240;
t104 = t144 * t178;
t194 = pkin(4) * t145 + t284;
t171 = t194 * qJD(4) + qJD(3);
t211 = t144 * t133;
t232 = -qJD(5) * t104 - t141 * t171 + t211 * t289;
t230 = t92 * t244;
t229 = t92 * t242;
t61 = t142 * t169;
t39 = t61 - t145 * t159 + (-t145 * t288 - pkin(4)) * t260;
t30 = t39 + t286;
t228 = t30 * t240;
t227 = t289 * t260;
t226 = t141 * t256;
t225 = qJ(6) * t254;
t135 = t140 ^ 2;
t223 = t135 * t246;
t219 = t145 * t241;
t217 = t141 * t238;
t105 = (pkin(5) * t141 + t289) * t145;
t215 = t105 * t240;
t212 = t141 * t132;
t210 = t142 * t133;
t209 = t87 * t268;
t208 = -t130 + t285;
t207 = pkin(5) + t257;
t22 = -t141 * t40 + t144 * t152;
t204 = t250 * t145;
t203 = t268 * qJD(3);
t202 = t249 * qJD(4);
t124 = 0.2e1 * t210;
t199 = t289 * t125;
t198 = t144 * t213;
t197 = t139 * t212;
t196 = t143 * t223;
t195 = qJD(2) * t140 * t268;
t192 = -pkin(2) * t146 - t267;
t191 = pkin(5) * t136 + t261;
t14 = pkin(5) * t92 - qJ(6) * t66 + t22;
t15 = -qJ(6) * t65 + t23;
t190 = t14 * t144 + t141 * t15;
t189 = t14 * t141 - t144 * t15;
t188 = t141 * t23 + t144 * t22;
t187 = t141 * t22 - t144 * t23;
t55 = t207 * t142 + t104 - t225;
t63 = -qJ(6) * t258 + t78;
t186 = t141 * t63 + t144 * t55;
t185 = t141 * t55 - t144 * t63;
t184 = t275 + t278;
t77 = t104 + t226;
t183 = t141 * t78 + t144 * t77;
t182 = t141 * t77 - t144 * t78;
t180 = t143 * t195;
t179 = t146 * t195;
t48 = t92 * t133 + t271;
t177 = t39 * t132 + t21 * t141;
t176 = -t21 * t144 + t39 * t240;
t47 = t92 * t132 + t141 * t64;
t175 = -t144 * t64 + t92 * t240;
t150 = t203 - t252 * pkin(10) + t64 * pkin(4) + (t200 + (-t284 - t288) * t260) * qJD(2);
t3 = -t141 * t150 + t144 * t297 + t40 * t240;
t98 = t142 * t132 + t141 * t133;
t4 = -t40 * t132 + t141 * t297 + t144 * t150;
t164 = -t188 * qJD(5) - t4 * t141 - t3 * t144;
t148 = -t38 * qJ(6) - t66 * qJD(6) + t4;
t1 = t148 + t287;
t2 = t37 * qJ(6) + t65 * qJD(6) + t3;
t163 = -t190 * qJD(5) - t1 * t141 - t2 * t144;
t42 = t145 * t155 - t61;
t162 = -t24 * t142 + t25 * t145 + (-t42 * t142 + t43 * t145) * qJD(4);
t161 = t272 - t273 + (t274 + t279) * qJD(5);
t49 = -t142 * t217 + t232;
t156 = -t78 * qJD(5) + t144 * t171;
t50 = t141 * t219 + t156;
t160 = -t183 * qJD(5) - t50 * t141 - t49 * t144;
t88 = -t280 * t240 - t236;
t89 = t280 * t132 - t237;
t158 = -t89 * t141 - t88 * t144 + (-t114 * t144 + t115 * t141) * qJD(5);
t153 = qJ(6) * t172 - t145 * t236 + t156;
t134 = qJ(3) * t291;
t123 = -0.2e1 * t212;
t122 = 0.2e1 * t212;
t108 = -0.2e1 * t196;
t107 = 0.2e1 * t196;
t106 = -0.2e1 * t251 * qJD(5);
t97 = t248 * t132;
t95 = t142 * t240 - t211;
t94 = t248 * t240;
t93 = qJD(4) * t204;
t85 = 0.2e1 * (t146 ^ 2 - t293) * t135 * qJD(2);
t84 = t168 + t235;
t83 = (-pkin(1) + t192) * t140;
t82 = -t206 - t101;
t80 = (-t143 * t243 + t145 * t246) * t140;
t79 = (-t143 * t133 - t142 * t246) * t140;
t76 = -t203 + t86;
t75 = -0.2e1 * t138 * t210 - 0.2e1 * t197;
t74 = -0.2e1 * t136 * t210 + 0.2e1 * t197;
t72 = t239 * t251 + t198;
t71 = -0.4e1 * t145 * t212 + t251 * t243;
t70 = (-t245 + (pkin(2) * t143 - t266) * qJD(2)) * t140;
t68 = 0.2e1 * t141 * t202 - 0.2e1 * t142 * t216;
t67 = -0.2e1 * t142 * t218 - 0.2e1 * t249 * t242;
t60 = t251 * t139 * t290 + 0.4e1 * t145 * t198;
t59 = (-0.1e1 + t250) * t124;
t57 = t203 + (-t288 * t260 + t200) * qJD(2);
t34 = t145 * t237 - qJ(6) * t213 + (t225 - t226) * qJD(5) + t232;
t33 = t207 * t133 + t153;
t29 = 0.2e1 * t66 * t38;
t28 = 0.2e1 * t65 * t37;
t27 = t66 * t132 + t272;
t26 = t65 * t240 - t273;
t19 = 0.2e1 * t38 * t92 + 0.2e1 * t64 * t66;
t18 = -0.2e1 * t37 * t92 - 0.2e1 * t64 * t65;
t17 = -t66 * t218 + (t145 * t38 - t66 * t243) * t144;
t16 = t65 * t216 + (t145 * t37 - t65 * t243) * t141;
t13 = (-t37 + t230) * t142 + (-t47 - t265) * t145;
t12 = (-t37 - t230) * t145 + (-t47 + t265) * t142;
t11 = (t38 - t229) * t142 + (-t175 + t264) * t145;
t10 = (-t38 - t229) * t145 + (t175 + t264) * t142;
t9 = -0.2e1 * t37 * t66 - 0.2e1 * t38 * t65;
t7 = -t184 * qJD(5) - t141 * t37 + t38 * t144;
t6 = t184 * t243 + (-t272 - t273 + (-t274 + t279) * qJD(5)) * t145;
t5 = (-t275 + t278) * t133 + t161 * t142;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t85, 0.2e1 * t179, t108, -0.2e1 * t180, 0, -0.2e1 * t135 * pkin(1) * t247 - 0.2e1 * t209, -0.2e1 * pkin(1) * t223 + 0.2e1 * t86 * t268 (t276 - t146 * t86 + (-t101 * t143 + t146 * t170) * qJD(2)) * t292, -0.2e1 * t101 * t86 + 0.2e1 * t170 * t87, 0, -0.2e1 * t179, 0.2e1 * t180, t107, t85, t108 (t276 - t146 * t76 + (t143 * t82 + t146 * t84) * qJD(2)) * t292, 0.2e1 * t209 + 0.2e1 * (t146 * t70 - t83 * t247) * t140, -0.2e1 * t76 * t268 + 0.2e1 * (-t143 * t70 - t83 * t246) * t140, 0.2e1 * t70 * t83 + 0.2e1 * t76 * t82 + 0.2e1 * t84 * t87, 0.2e1 * t181 * t173, -0.2e1 * t173 * t92 - 0.2e1 * t181 * t64 (t252 * t143 + (t142 * t293 * t140 + t146 * t181) * qJD(2)) * t292, t51 (-t143 * t64 - t246 * t92) * t292, t107, 0.2e1 * t57 * t92 + 0.2e1 * t69 * t64 + 0.2e1 * (t143 * t25 + t246 * t42) * t140, 0.2e1 * t57 * t181 + 0.2e1 * t69 * t252 + 0.2e1 * (t24 * t143 + (t143 * t277 - t146 * t43) * qJD(2)) * t140, -0.2e1 * t173 * t42 - 0.2e1 * t181 * t25 + 0.2e1 * t24 * t92 - 0.2e1 * t43 * t64, -0.2e1 * t24 * t43 + 0.2e1 * t25 * t42 + 0.2e1 * t57 * t69, t29, t9, t19, t28, t18, t51, 0.2e1 * t21 * t65 + 0.2e1 * t22 * t64 + 0.2e1 * t37 * t39 + 0.2e1 * t4 * t92, 0.2e1 * t21 * t66 - 0.2e1 * t23 * t64 + 0.2e1 * t3 * t92 + 0.2e1 * t38 * t39, -0.2e1 * t22 * t38 - 0.2e1 * t23 * t37 + 0.2e1 * t3 * t65 - 0.2e1 * t4 * t66, 0.2e1 * t21 * t39 + 0.2e1 * t22 * t4 - 0.2e1 * t23 * t3, t29, t9, t19, t28, t18, t51, 0.2e1 * t1 * t92 + 0.2e1 * t14 * t64 + 0.2e1 * t30 * t37 + 0.2e1 * t65 * t8, -0.2e1 * t15 * t64 + 0.2e1 * t2 * t92 + 0.2e1 * t30 * t38 + 0.2e1 * t66 * t8, -0.2e1 * t1 * t66 - 0.2e1 * t14 * t38 - 0.2e1 * t15 * t37 + 0.2e1 * t2 * t65, 0.2e1 * t1 * t14 - 0.2e1 * t15 * t2 + 0.2e1 * t30 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, 0, -t214, 0, -t87, t86, 0, 0, 0, -t125, t214, 0, 0, 0 (t192 * qJD(2) + qJD(3) * t146) * t140, t87, 0.2e1 * t203 - t86, -pkin(2) * t87 - qJ(3) * t76 - qJD(3) * t82, -t181 * t243 + t294, -t145 * t64 - t173 * t142 + (t142 * t92 - t145 * t181) * qJD(4), t80, t48, t79, 0, -t145 * t199 + qJ(3) * t64 + qJD(3) * t92 + t57 * t142 + (t142 * t227 + t145 * t69) * qJD(4), t142 * t199 + qJD(3) * t181 + qJ(3) * t173 + t57 * t145 + (t145 * t227 - t277) * qJD(4) (t289 * t64 + t24) * t142 + (t173 * t289 - t25) * t145 + ((t289 * t92 - t43) * t145 + (-t181 * t289 + t42) * t142) * qJD(4), t57 * qJ(3) + t69 * qJD(3) - t162 * t289, t17, t6, t11, t16, t13, t48, t50 * t92 + t77 * t64 + (t4 + (-t141 * t39 - t289 * t65) * qJD(4)) * t142 + (qJD(4) * t22 + t289 * t37 + t177) * t145, t49 * t92 - t78 * t64 + (t3 + (-t144 * t39 - t289 * t66) * qJD(4)) * t142 + (-qJD(4) * t23 + t289 * t38 - t176) * t145, -t37 * t78 - t38 * t77 + t49 * t65 - t50 * t66 + t188 * t243 + (qJD(5) * t187 + t141 * t3 - t144 * t4) * t145, t22 * t50 - t23 * t49 - t3 * t78 + t4 * t77 - (-t21 * t145 + t243 * t39) * t289, t17, t6, t11, t16, t13, t48, t105 * t37 + t33 * t92 + t55 * t64 + t65 * t73 + (-t244 * t30 + t1) * t142 + (qJD(4) * t14 + t132 * t30 + t282) * t145, t105 * t38 + t34 * t92 - t63 * t64 + t66 * t73 + (-t242 * t30 + t2) * t142 + (-qJD(4) * t15 - t228 + t281) * t145, -t33 * t66 + t34 * t65 - t37 * t63 - t38 * t55 + t190 * t243 + (qJD(5) * t189 - t1 * t144 + t141 * t2) * t145, t1 * t55 + t105 * t8 + t14 * t33 - t15 * t34 - t2 * t63 + t30 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, t134, -0.2e1 * t210, 0.2e1 * t202, 0, t124, 0, 0, 0.2e1 * qJ(3) * t133 + 0.2e1 * qJD(3) * t142, -0.2e1 * qJ(3) * t243 + 0.2e1 * qJD(3) * t145, 0, t134, t75, t60, t67, t74, t68, t124, 0.2e1 * t139 * t144 * t238 + 0.2e1 * t50 * t142 + 0.2e1 * (t77 - 0.2e1 * t226) * t133, -0.2e1 * t139 * t217 + 0.2e1 * t49 * t142 + 0.2e1 * (-t78 - 0.2e1 * t121) * t133, 0.2e1 * t183 * t243 + 0.2e1 * (qJD(5) * t182 + t141 * t49 - t144 * t50) * t145, -0.2e1 * t210 * t289 ^ 2 - 0.2e1 * t78 * t49 + 0.2e1 * t77 * t50, t75, t60, t67, t74, t68, t124, 0.2e1 * (-t105 * t244 + t33) * t142 + 0.2e1 * (qJD(4) * t55 + t105 * t132 + t270) * t145, 0.2e1 * (-t105 * t242 + t34) * t142 + 0.2e1 * (-qJD(4) * t63 - t215 + t269) * t145, 0.2e1 * t186 * t243 + 0.2e1 * (qJD(5) * t185 + t141 * t34 - t144 * t33) * t145, 0.2e1 * t105 * t73 + 0.2e1 * t33 * t55 - 0.2e1 * t34 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, 0, 0, t87, 0, 0, 0, 0, 0, 0, t80, t79, -t271 - t294 + (t142 * t181 - t92 * t145) * qJD(4), t162, 0, 0, 0, 0, 0, 0, t12, t10, t5 (-qJD(4) * t187 - t21) * t145 + (qJD(4) * t39 + t164) * t142, 0, 0, 0, 0, 0, 0, t12, t10, t5 (-qJD(4) * t189 - t8) * t145 + (qJD(4) * t30 + t163) * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t94, 0, -t182 * t133 + (t160 + 0.2e1 * t219) * t142, 0, 0, 0, 0, 0, 0, -t97, t94, 0 (-qJD(4) * t185 - t73) * t145 + (qJD(4) * t105 - qJD(5) * t186 - t33 * t141 - t34 * t144) * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, 0, -t64, t125, t25, t24, 0, 0, t27, t7, t47, t26, -t175, 0, -pkin(4) * t37 - pkin(10) * t47 + t176, -pkin(4) * t38 + pkin(10) * t175 + t177, pkin(10) * t161 + t164, -t21 * pkin(4) + pkin(10) * t164, t27, t7, t47, t26, -t175, 0, t114 * t64 + t130 * t37 - t281 + t89 * t92 + (t30 + t286) * t240, t115 * t64 + t130 * t38 + t282 + t88 * t92 + (pkin(5) * t278 + t144 * t30) * qJD(5), -t114 * t38 + t115 * t37 + t65 * t88 - t66 * t89 + t163, pkin(5) * t228 + t1 * t114 + t115 * t2 + t130 * t8 + t14 * t89 - t15 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, 0, -t133, 0, t220, t219, 0, 0, -t72, t71, t98, t72, -t95, 0 (t141 * t253 - t144 * t194) * qJD(5) + (t141 * t193 + t121) * qJD(4) (t141 * t194 + t144 * t253) * qJD(5) + (-pkin(10) * t254 + (pkin(4) * t144 - t257) * t142) * qJD(4), t160, pkin(4) * t220 + pkin(10) * t160, -t72, t71, t98, t72, -t95, 0, t89 * t142 - t269 + (-t130 * t141 * t142 + t263) * qJD(4) + (t105 * t141 + t145 * t191) * qJD(5), t270 + t88 * t142 + (-t142 * t261 + t262) * qJD(4) + (t105 * t144 + t208 * t258) * qJD(5) (t114 * t243 - t145 * t89 - t34 + (-t55 + t262) * qJD(5)) * t144 + (-t115 * t243 + t145 * t88 - t33 + (-t63 + t263) * qJD(5)) * t141, pkin(5) * t215 + t114 * t33 + t115 * t34 + t130 * t73 + t55 * t89 - t63 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, -t133, 0, 0, 0, 0, 0, 0, 0, 0, -t172, -t99, t93 (pkin(10) * t204 - t283) * qJD(4), 0, 0, 0, 0, 0, 0, -t172, -t99, t93 (-t233 + (-t114 * t141 - t115 * t144) * qJD(4)) * t145 + (qJD(4) * t130 + t158) * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t106, 0, t123, 0, 0, pkin(4) * t234, -0.2e1 * pkin(4) * t132, 0, 0, t122, t106, 0, t123, 0, 0, t208 * t234, t191 * t290, 0.2e1 * t158, 0.2e1 * t114 * t89 + 0.2e1 * t115 * t88 + 0.2e1 * t130 * t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, -t37, t64, t4, t3, 0, 0, 0, 0, t38, 0, -t37, t64, t148 + 0.2e1 * t287, t2, -t38 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, 0, -t99, t133, t50, t49, 0, 0, 0, 0, -t172, 0, -t99, t133 (0.2e1 * pkin(5) + t257) * t133 + t153, t34, t172 * pkin(5), t33 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, t95, 0, 0, 0, 0, 0, 0, 0, 0, -t98, t95, 0, -t98 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, 0, -t240, 0, -pkin(10) * t132, pkin(10) * t240, 0, 0, 0, 0, t132, 0, -t240, 0, t89, t88, -pkin(5) * t132, t89 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t38, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, -t172, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, t132, 0, t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t20;
