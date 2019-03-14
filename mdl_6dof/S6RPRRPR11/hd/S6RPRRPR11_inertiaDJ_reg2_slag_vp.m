% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR11_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_inertiaDJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:44:33
% EndTime: 2019-03-09 05:45:03
% DurationCPUTime: 12.09s
% Computational Cost: add. (19096->567), mult. (55743->1020), div. (0->0), fcn. (59754->14), ass. (0->239)
t124 = cos(qJ(4));
t118 = sin(pkin(12));
t120 = sin(pkin(6));
t122 = cos(pkin(12));
t283 = cos(pkin(7));
t293 = sin(qJ(3));
t224 = t283 * t293;
t214 = t122 * t224;
t119 = sin(pkin(7));
t284 = cos(pkin(6));
t238 = t284 * t119;
t215 = t293 * t238;
t295 = cos(qJ(3));
t298 = t120 * (t295 * t118 + t214) + t215;
t304 = t298 * qJD(3);
t306 = t124 * t304;
t258 = pkin(1) * t284;
t202 = pkin(2) * t284 + t122 * t258;
t278 = t118 * t120;
t178 = (-pkin(9) * t283 - qJ(2)) * t278 + t202;
t260 = -pkin(2) * t122 - pkin(1);
t279 = t118 * t119;
t206 = (-pkin(9) * t279 + t260) * t120;
t305 = t119 * t206 + t178 * t283;
t240 = t120 * t283;
t276 = t120 * t122;
t170 = t118 * t258 + qJ(2) * t276 + (t122 * t240 + t238) * pkin(9);
t43 = t295 * t170 + t305 * t293;
t303 = qJD(3) * t43;
t121 = cos(pkin(13));
t117 = sin(pkin(13));
t225 = t283 * t295;
t184 = t225 * t276 + t295 * t238;
t180 = t124 * t184;
t149 = (t117 * t298 + t121 * t180) * qJD(3);
t208 = t119 * t276 - t283 * t284;
t205 = t208 * qJD(4);
t123 = sin(qJ(4));
t167 = qJD(4) * t298;
t68 = t123 * t167;
t236 = t293 * t278;
t95 = qJD(3) * t236;
t145 = (-t68 + (-t95 - t205) * t124) * t121 + t149;
t207 = t236 - t184;
t204 = qJD(4) * t207;
t302 = t123 * t304 + t124 * t204;
t292 = sin(qJ(6));
t242 = qJD(6) * t292;
t294 = cos(qJ(6));
t243 = qJD(6) * t294;
t84 = t117 * t242 - t121 * t243;
t256 = t119 * t293;
t86 = t123 * t256 - t124 * t283;
t91 = t292 * t117 - t294 * t121;
t271 = qJD(4) * t123;
t299 = -qJD(3) * t184 + t95;
t47 = t124 * t167 + (-t299 - t205) * t123;
t55 = t123 * t298 + t124 * t208;
t29 = -t124 * t47 + t55 * t271;
t116 = t124 ^ 2;
t301 = (t123 ^ 2 - t116) * qJD(4);
t300 = t293 * t170 - t305 * t295;
t297 = -0.2e1 * t84;
t255 = t294 * t117;
t92 = t292 * t121 + t255;
t85 = t92 * qJD(6);
t296 = 0.2e1 * t85;
t291 = pkin(4) * t123;
t290 = pkin(4) * t124;
t289 = pkin(10) * t123;
t257 = t119 * t295;
t230 = qJD(3) * t257;
t87 = t123 * t283 + t124 * t256;
t75 = qJD(4) * t87 + t123 * t230;
t61 = t86 * t75;
t288 = qJ(5) + pkin(11);
t140 = -pkin(10) * t208 + t43;
t142 = pkin(3) * t207 - pkin(10) * t298 - t119 * t178 + t206 * t283;
t23 = t123 * t142 + t124 * t140;
t133 = qJ(5) * t207 + t23;
t154 = -t123 * t208 + t124 * t298;
t41 = pkin(3) * t208 + t300;
t136 = t55 * pkin(4) - qJ(5) * t154 + t41;
t13 = t117 * t136 + t121 * t133;
t287 = t121 * t47;
t285 = t75 * t123;
t274 = t121 * t124;
t107 = pkin(10) * t274;
t222 = -t123 * qJ(5) - t290;
t212 = -pkin(3) + t222;
t80 = t117 * t212 + t107;
t282 = qJ(5) * t117;
t281 = t117 * t123;
t280 = t117 * t124;
t277 = t119 * qJ(2);
t275 = t121 * t123;
t272 = qJD(2) * t120;
t270 = qJD(4) * t124;
t269 = qJD(5) * t117;
t268 = qJD(5) * t121;
t267 = qJD(5) * t124;
t266 = t123 * qJD(5);
t32 = 0.2e1 * t55 * t47;
t265 = -0.2e1 * pkin(3) * qJD(4);
t264 = pkin(10) * t280;
t263 = pkin(10) * t275;
t262 = pkin(10) * t271;
t111 = pkin(10) * t270;
t259 = pkin(10) * t117 + pkin(5);
t252 = t122 * t272;
t251 = t117 * t270;
t250 = t117 * t266;
t249 = t121 * t270;
t248 = t121 * t266;
t247 = t123 * t270;
t246 = qJD(3) * t293;
t245 = qJD(4) * t295;
t244 = qJD(5) * t292;
t241 = t294 * qJD(5);
t237 = 0.2e1 * t247;
t234 = t272 * t279;
t233 = t117 * t249;
t229 = t119 * t246;
t112 = t117 ^ 2;
t114 = t121 ^ 2;
t228 = 0.2e1 * (t112 + t114) * qJD(5);
t227 = -0.2e1 * t301;
t226 = t288 * t292;
t139 = qJD(4) * t140;
t141 = qJD(4) * t142;
t134 = -t123 * t139 + t124 * t141;
t203 = t207 * qJD(5);
t131 = t95 * t289 + t134 + t203;
t135 = -t124 * t300 + t123 * (pkin(3) * t298 - pkin(10) * t184) + t298 * qJ(5);
t138 = -qJ(5) * t180 + t43;
t151 = t47 * pkin(4) - qJD(5) * t154;
t193 = -t124 * t205 - t68;
t185 = -t95 * t124 + t193;
t147 = -qJ(5) * t185 + t151;
t196 = -t118 * t224 + t295 * t122;
t179 = t123 * t279 + t124 * t196;
t197 = t118 * t225 + t293 * t122;
t7 = -t117 * t131 + t121 * t147 + (-t117 * t179 + t121 * t197) * t272 + (-t117 * t135 + t121 * t138) * qJD(3);
t8 = t121 * t131 + t117 * t147 + (t117 * t197 + t121 * t179) * t272 + (t117 * t138 + t121 * t135) * qJD(3);
t223 = -t7 * t117 + t8 * t121;
t221 = -qJ(5) * t124 + t291;
t74 = qJD(4) * t86 - t124 * t230;
t62 = t117 * t74 + t121 * t229;
t63 = t117 * t229 - t121 * t74;
t220 = -t117 * t62 + t121 * t63;
t70 = -t248 + (pkin(10) * t281 + t121 * t221) * qJD(4);
t71 = -t250 + (t117 * t221 - t263) * qJD(4);
t219 = -t117 * t70 + t121 * t71;
t213 = qJD(2) * t118 * t240;
t137 = qJD(3) * t300 + t293 * t213 - t295 * t252;
t148 = pkin(3) * t304 + pkin(10) * t299 + t234;
t17 = -t123 * t148 + t124 * t137 - t134;
t18 = (-t139 + t148) * t124 + (t137 - t141) * t123;
t218 = -t18 * t123 - t17 * t124;
t217 = -0.2e1 * t284 * t272;
t216 = t288 * t255;
t211 = t283 * t260;
t210 = -t124 * t288 + t291;
t209 = 0.2e1 * (t118 ^ 2 + t122 ^ 2) * t120 ^ 2 * qJD(2);
t72 = -t117 * t87 - t121 * t257;
t73 = -t117 * t257 + t121 * t87;
t52 = t292 * t72 + t294 * t73;
t199 = t293 * pkin(3) - t295 * pkin(10) + t277;
t198 = t154 * t117 - t121 * t207;
t194 = t119 * t202;
t190 = qJD(4) * t198;
t189 = t285 - t124 * t74 + (-t123 * t87 + t124 * t86) * qJD(4);
t188 = t294 * t198;
t187 = t292 * t198;
t186 = t124 * t190;
t183 = -pkin(10) * t214 + t211;
t181 = -t259 * t124 + (-t123 * t288 - pkin(3) - t290) * t121;
t175 = t294 * t181;
t174 = t292 * t181;
t173 = t250 - (t117 * t210 - t263) * qJD(4);
t168 = -t248 + (t121 * t210 + t123 * t259) * qJD(4);
t164 = -0.2e1 * t304;
t153 = -t124 * t299 + t193;
t158 = t153 * t117 - t121 * t304;
t157 = t117 * t158;
t156 = t158 * t121;
t40 = t123 * t140;
t155 = t124 * (-pkin(3) * t184 - pkin(10) * t215 - t194) - t40;
t152 = qJD(3) * t154;
t144 = t145 * t117;
t143 = t145 * t121;
t14 = -pkin(4) * t304 - t18;
t37 = t295 * t213 + t293 * t252 + t303;
t132 = -qJ(5) * t153 + t151 + t37;
t12 = -t117 * t133 + t121 * t136;
t48 = t117 * t207 + t121 * t154;
t130 = t55 * pkin(5) - t48 * pkin(11) + t12;
t129 = t294 * t130;
t128 = t292 * t130;
t127 = qJ(5) * t304 - t17 + t203;
t126 = -pkin(11) * t158 + t117 * t132 + t121 * t127;
t125 = t47 * pkin(5) - t145 * pkin(11) - t117 * t127 + t121 * t132;
t110 = -pkin(5) * t121 - pkin(4);
t103 = -0.2e1 * t247;
t99 = t288 * t121;
t94 = (pkin(5) * t117 + pkin(10)) * t123;
t88 = pkin(5) * t251 + t111;
t83 = t91 * t123;
t82 = t92 * t123;
t79 = t121 * t212 - t264;
t78 = -t117 * t226 + t294 * t99;
t77 = -t292 * t99 - t216;
t76 = -pkin(11) * t281 + t80;
t65 = -t123 * t84 + t270 * t92;
t64 = t123 * t85 - t294 * t249 + t292 * t251;
t60 = -t99 * t243 - t121 * t244 + (qJD(6) * t226 - t241) * t117;
t59 = qJD(6) * t216 + t117 * t244 - t121 * t241 + t242 * t99;
t53 = -t194 + (t118 * t277 + t211) * t120;
t51 = -t292 * t73 + t294 * t72;
t50 = t294 * t76 + t174;
t49 = -t292 * t76 + t175;
t31 = -qJD(6) * t174 + t294 * t168 + t292 * t173 - t76 * t243;
t30 = -qJD(6) * t175 - t292 * t168 + t294 * t173 + t76 * t242;
t28 = -t52 * qJD(6) - t292 * t63 + t294 * t62;
t27 = t73 * t242 - t72 * t243 - t292 * t62 - t294 * t63;
t26 = t294 * t48 - t187;
t25 = t292 * t48 + t188;
t22 = t124 * (t118 * t199 + t183) * t120 + t155;
t21 = t184 * pkin(4) + (-t124 * t183 + (-t293 * pkin(4) - t124 * t199) * t118) * t120 - t155;
t19 = -pkin(4) * t207 + pkin(5) * t198 - t124 * t142 + t40;
t16 = -qJD(6) * t187 + t292 * t145 + t294 * t158 + t48 * t243;
t15 = qJD(6) * t188 - t294 * t145 + t292 * t158 + t48 * t242;
t11 = -pkin(11) * t198 + t13;
t10 = pkin(5) * t158 + t14;
t4 = t294 * t11 + t128;
t3 = -t292 * t11 + t129;
t2 = -qJD(6) * t128 - t11 * t243 + t294 * t125 - t292 * t126;
t1 = -qJD(6) * t129 + t11 * t242 - t292 * t125 - t294 * t126;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118 * t217, t122 * t217, t209, qJ(2) * t209, -0.2e1 * t298 * t299, 0.2e1 * t207 * t299 - 0.2e1 * t298 * t304, 0.2e1 * t299 * t208, -t207 * t164, -t208 * t164, 0, 0.2e1 * t207 * t234 + 0.2e1 * t208 * t37 + 0.2e1 * t304 * t53, -0.2e1 * t137 * t208 + 0.2e1 * t234 * t298 - 0.2e1 * t299 * t53, 0.2e1 * t137 * t207 + 0.2e1 * t298 * t37 - 0.2e1 * t299 * t300 - 0.2e1 * t304 * t43, 0.2e1 * (t196 * t43 + t279 * t53) * t272 - 0.2e1 * (-t37 + t303) * t300, 0.2e1 * t154 * t153, -0.2e1 * t153 * t55 - 0.2e1 * t154 * t47, 0.2e1 * t152 * t298 + 0.2e1 * t153 * t207, t32, -0.2e1 * t207 * t47 - 0.2e1 * t304 * t55, 0.2e1 * t207 * t304, 0.2e1 * t18 * t207 + 0.2e1 * t22 * t304 + 0.2e1 * t37 * t55 + 0.2e1 * t41 * t47, 0.2e1 * t153 * t41 + 0.2e1 * t154 * t37 + 0.2e1 * t17 * t207 - 0.2e1 * t23 * t304, -0.2e1 * t153 * t22 - 0.2e1 * t154 * t18 + 0.2e1 * t17 * t55 - 0.2e1 * t23 * t47, -0.2e1 * t17 * t23 + 0.2e1 * t18 * t22 + 0.2e1 * t37 * t41, 0.2e1 * t48 * t145, -0.2e1 * t145 * t198 - 0.2e1 * t158 * t48, 0.2e1 * t48 * t47 + 0.2e1 * (t121 * t185 + t149) * t55, 0.2e1 * t198 * t158, -0.2e1 * t158 * t55 - 0.2e1 * t198 * t47, t32, 0.2e1 * t12 * t47 + 0.2e1 * t14 * t198 + 0.2e1 * t158 * t21 + 0.2e1 * t7 * t55, -0.2e1 * t13 * t47 + 0.2e1 * t14 * t48 + 0.2e1 * t145 * t21 - 0.2e1 * t8 * t55, -0.2e1 * t12 * t145 - 0.2e1 * t13 * t158 - 0.2e1 * t198 * t8 - 0.2e1 * t7 * t48, 0.2e1 * t12 * t7 + 0.2e1 * t13 * t8 + 0.2e1 * t14 * t21, -0.2e1 * t26 * t15, 0.2e1 * t15 * t25 - 0.2e1 * t16 * t26, -0.2e1 * t15 * t55 + 0.2e1 * t26 * t47, 0.2e1 * t25 * t16, -0.2e1 * t16 * t55 - 0.2e1 * t25 * t47, t32, 0.2e1 * t10 * t25 + 0.2e1 * t16 * t19 + 0.2e1 * t2 * t55 + 0.2e1 * t3 * t47, 0.2e1 * t1 * t55 + 0.2e1 * t10 * t26 - 0.2e1 * t15 * t19 - 0.2e1 * t4 * t47, 0.2e1 * t1 * t25 + 0.2e1 * t15 * t3 - 0.2e1 * t16 * t4 - 0.2e1 * t2 * t26, -0.2e1 * t1 * t4 + 0.2e1 * t10 * t19 + 0.2e1 * t2 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t208 * t256 + t283 * t298) * qJD(3), t208 * t230 - t283 * t299, -t207 * t230 + t257 * t299, t119 * t213 - t137 * t256 + t229 * t300 + t230 * t43 - t257 * t37, 0, 0, 0, 0, 0, 0, -t207 * t75 + t229 * t55 - t257 * t47 - t304 * t86, t152 * t256 - t153 * t257 + t207 * t74 - t304 * t87, t153 * t86 + t154 * t75 - t87 * t47 + t74 * t55, -t17 * t87 - t18 * t86 - t22 * t75 - t23 * t74 + (t41 * t246 - t37 * t295) * t119, 0, 0, 0, 0, 0, 0, t158 * t86 + t198 * t75 + t72 * t47 + t62 * t55, t145 * t86 - t73 * t47 + t75 * t48 - t63 * t55, -t145 * t72 - t158 * t73 - t198 * t63 - t62 * t48, t12 * t62 + t13 * t63 + t14 * t86 + t21 * t75 + t7 * t72 + t73 * t8, 0, 0, 0, 0, 0, 0, t16 * t86 + t25 * t75 + t28 * t55 + t47 * t51, -t15 * t86 + t26 * t75 + t27 * t55 - t47 * t52, t15 * t51 - t16 * t52 + t25 * t27 - t26 * t28, -t1 * t52 + t10 * t86 + t19 * t75 + t2 * t51 - t27 * t4 + t28 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t119 ^ 2 * t295 * t246 - 0.2e1 * t87 * t74 + 0.2e1 * t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t62 * t72 + 0.2e1 * t63 * t73 + 0.2e1 * t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t27 * t52 + 0.2e1 * t28 * t51 + 0.2e1 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, 0, -t304, 0, -t37, t137, 0, 0, t116 * t167 + (-t68 + (-0.2e1 * t205 - t299) * t124) * t123, -t123 * t47 + t124 * t153 - t154 * t271 - t270 * t55, t302, t29, -t123 * t204 + t306, 0, -pkin(3) * t47 - pkin(10) * t302 - t37 * t124 + t271 * t41, -pkin(3) * t153 - pkin(10) * t306 + t37 * t123 + t207 * t262 + t270 * t41, pkin(10) * t29 + t111 * t154 + t153 * t289 - t22 * t270 - t23 * t271 + t218, -t37 * pkin(3) + ((-t23 * t123 - t22 * t124) * qJD(4) + t218) * pkin(10), t123 * t143 + t249 * t48, -t121 * t186 - t251 * t48 + (-t144 - t156) * t123, -t124 * t145 + t249 * t55 + t271 * t48 + t275 * t47, t117 * t186 + t123 * t157, -t123 * t190 + t124 * t158 - t251 * t55 - t281 * t47, t29, t111 * t198 + t12 * t271 - t7 * t124 + t14 * t281 + t158 * t289 + t21 * t251 + t79 * t47 + t70 * t55, t111 * t48 + t8 * t124 - t13 * t271 + t14 * t275 + t145 * t289 + t21 * t249 - t80 * t47 - t71 * t55, -t12 * t249 - t13 * t251 - t145 * t79 - t158 * t80 - t198 * t71 - t275 * t7 - t281 * t8 - t70 * t48, t12 * t70 + t13 * t71 + t7 * t79 + t8 * t80 + (t14 * t123 + t21 * t270) * pkin(10), t15 * t83 - t26 * t64, t15 * t82 + t16 * t83 + t25 * t64 - t26 * t65, t124 * t15 + t26 * t271 - t83 * t47 - t64 * t55, t16 * t82 + t25 * t65, t124 * t16 - t25 * t271 - t82 * t47 - t65 * t55, t29, t10 * t82 - t124 * t2 + t94 * t16 + t19 * t65 + t88 * t25 + t271 * t3 + t31 * t55 + t49 * t47, -t1 * t124 - t10 * t83 - t94 * t15 - t19 * t64 + t88 * t26 - t271 * t4 + t30 * t55 - t50 * t47, t1 * t82 + t15 * t49 - t16 * t50 + t2 * t83 + t25 * t30 - t26 * t31 + t3 * t64 - t4 * t65, -t1 * t50 + t10 * t94 + t19 * t88 + t2 * t49 + t3 * t31 - t30 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t229, -t230, 0, 0, 0, 0, 0, 0, 0, 0 (-t123 * t245 - t124 * t246) * t119 (t123 * t246 - t124 * t245) * t119, t189, -pkin(3) * t229 + pkin(10) * t189, 0, 0, 0, 0, 0, 0, t75 * t281 - t124 * t62 + (t123 * t72 + t280 * t86) * qJD(4), t75 * t275 + t124 * t63 + (-t123 * t73 + t274 * t86) * qJD(4) (-t117 * t63 - t121 * t62) * t123 + (-t117 * t73 - t121 * t72) * t270, t62 * t79 + t63 * t80 + t72 * t70 + t73 * t71 + (t270 * t86 + t285) * pkin(10), 0, 0, 0, 0, 0, 0, -t124 * t28 + t271 * t51 + t86 * t65 + t75 * t82, -t124 * t27 - t271 * t52 - t86 * t64 - t75 * t83, t27 * t82 + t28 * t83 + t51 * t64 - t52 * t65, -t27 * t50 + t28 * t49 - t30 * t52 + t31 * t51 + t75 * t94 + t86 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, t227, 0, t103, 0, 0, t123 * t265, t124 * t265, 0, 0, t114 * t237, -0.4e1 * t123 * t233, 0.2e1 * t121 * t301, t112 * t237, t117 * t227, t103, -0.2e1 * t124 * t70 + 0.2e1 * (t79 + 0.2e1 * t264) * t271, 0.2e1 * t124 * t71 + 0.2e1 * (-t80 + 0.2e1 * t107) * t271, 0.2e1 * (-t117 * t71 - t121 * t70) * t123 + 0.2e1 * (-t117 * t80 - t121 * t79) * t270, 0.2e1 * pkin(10) ^ 2 * t247 + 0.2e1 * t79 * t70 + 0.2e1 * t80 * t71, 0.2e1 * t83 * t64, 0.2e1 * t64 * t82 + 0.2e1 * t65 * t83, 0.2e1 * t124 * t64 - 0.2e1 * t271 * t83, 0.2e1 * t82 * t65, 0.2e1 * t124 * t65 - 0.2e1 * t271 * t82, t103, -0.2e1 * t124 * t31 + 0.2e1 * t271 * t49 + 0.2e1 * t94 * t65 + 0.2e1 * t88 * t82, -0.2e1 * t124 * t30 - 0.2e1 * t271 * t50 - 0.2e1 * t94 * t64 - 0.2e1 * t88 * t83, 0.2e1 * t30 * t82 + 0.2e1 * t31 * t83 + 0.2e1 * t49 * t64 - 0.2e1 * t50 * t65, -0.2e1 * t30 * t50 + 0.2e1 * t31 * t49 + 0.2e1 * t88 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, 0, -t47, t304, t18, t17, 0, 0, t144, -t157 + t143, t117 * t47, -t156, t287, 0, -pkin(4) * t158 - t14 * t121 - t269 * t55 - t282 * t47, -pkin(4) * t145 - qJ(5) * t287 + t14 * t117 - t268 * t55, -qJ(5) * t156 + t145 * t282 - t198 * t268 + t269 * t48 + t223, -pkin(4) * t14 + (-t117 * t12 + t121 * t13) * qJD(5) + t223 * qJ(5), -t15 * t92 - t26 * t84, t15 * t91 - t16 * t92 + t25 * t84 - t26 * t85, t47 * t92 - t55 * t84, t16 * t91 + t25 * t85, -t47 * t91 - t55 * t85, 0, t10 * t91 + t110 * t16 + t19 * t85 + t47 * t77 + t55 * t60, t10 * t92 - t110 * t15 - t19 * t84 - t47 * t78 + t55 * t59, t1 * t91 + t15 * t77 - t16 * t78 - t2 * t92 + t25 * t59 - t26 * t60 + t3 * t84 - t4 * t85, -t1 * t78 + t10 * t110 + t2 * t77 + t3 * t60 - t4 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t74, 0, 0, 0, 0, 0, 0, 0, 0, -t75 * t121, t75 * t117, t220, -pkin(4) * t75 + (-t117 * t72 + t121 * t73) * qJD(5) + t220 * qJ(5), 0, 0, 0, 0, 0, 0, t75 * t91 + t85 * t86, t75 * t92 - t84 * t86, t27 * t91 - t28 * t92 + t51 * t84 - t52 * t85, t110 * t75 - t27 * t78 + t28 * t77 + t51 * t60 - t52 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270, 0, -t271, 0, -t111, t262, 0, 0, t233 (-t112 + t114) * t270, t117 * t271, -t233, t121 * t271, 0, t117 * t267 + (t117 * t222 - t107) * qJD(4), t121 * t267 + (t121 * t222 + t264) * qJD(4), t219, -pkin(4) * t111 + (-t117 * t79 + t80 * t121) * qJD(5) + t219 * qJ(5), -t64 * t92 + t83 * t84, t64 * t91 - t65 * t92 + t82 * t84 + t83 * t85, t124 * t84 + t271 * t92, t65 * t91 + t82 * t85, t124 * t85 - t271 * t91, 0, t110 * t65 - t124 * t60 + t271 * t77 + t94 * t85 + t88 * t91, -t110 * t64 - t124 * t59 - t271 * t78 - t94 * t84 + t88 * t92, t30 * t91 - t31 * t92 + t49 * t84 - t50 * t85 + t59 * t82 + t60 * t83 + t64 * t77 - t65 * t78, t110 * t88 - t30 * t78 + t31 * t77 + t49 * t60 - t50 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, qJ(5) * t228, t92 * t297, 0.2e1 * t84 * t91 - 0.2e1 * t85 * t92, 0, t91 * t296, 0, 0, t110 * t296, t110 * t297, 0.2e1 * t59 * t91 - 0.2e1 * t60 * t92 + 0.2e1 * t77 * t84 - 0.2e1 * t78 * t85, -0.2e1 * t59 * t78 + 0.2e1 * t60 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t145, 0, t14, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t249, 0, t111, 0, 0, 0, 0, 0, 0, t65, -t64, 0, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t84, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, -t16, t47, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, 0, -t65, t271, t31, t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, 0, -t85, 0, t60, t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;