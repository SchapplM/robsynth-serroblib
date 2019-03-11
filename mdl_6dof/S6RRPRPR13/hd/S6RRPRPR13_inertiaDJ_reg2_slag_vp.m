% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPR13_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:29:34
% EndTime: 2019-03-09 11:29:57
% DurationCPUTime: 9.09s
% Computational Cost: add. (7529->519), mult. (19219->930), div. (0->0), fcn. (18075->10), ass. (0->250)
t118 = sin(pkin(11));
t123 = cos(qJ(4));
t287 = pkin(2) + pkin(9);
t300 = t123 * t287;
t301 = t118 * t300;
t120 = cos(pkin(11));
t241 = t120 * t300;
t284 = sin(qJ(6));
t285 = cos(qJ(6));
t168 = -t284 * t118 + t285 * t120;
t119 = sin(pkin(6));
t122 = sin(qJ(2));
t124 = cos(qJ(2));
t270 = cos(pkin(6));
t240 = pkin(1) * t270;
t212 = t124 * t240;
t165 = -t270 * pkin(2) - t212;
t156 = t270 * pkin(9) - t165;
t286 = pkin(3) + pkin(8);
t247 = t119 * t286;
t145 = t122 * t247 - t156;
t254 = qJD(3) * t122;
t268 = qJ(3) * t124;
t299 = -qJD(4) * t145 - t119 * (-t254 + (t287 * t122 - t268) * qJD(2));
t213 = t122 * t240;
t215 = t124 * t247;
t160 = t213 + t215;
t199 = -t287 * t124 - pkin(1);
t269 = qJ(3) * t122;
t166 = t119 * (t199 - t269);
t298 = qJD(2) * t160 - qJD(4) * t166;
t228 = qJD(6) * t284;
t229 = qJD(6) * t285;
t79 = t118 * t228 - t120 * t229;
t297 = t123 * t79;
t255 = qJD(2) * t124;
t104 = t119 * t255;
t121 = sin(qJ(4));
t256 = qJD(2) * t122;
t233 = t119 * t256;
t222 = qJD(4) * t270;
t263 = t119 * t124;
t235 = qJD(4) * t263;
t271 = -t121 * t222 - t123 * t235;
t173 = t121 * t233 + t271;
t153 = t118 * t104 + t120 * t173;
t296 = t153 * t118;
t74 = t168 * t121;
t295 = t173 * t123;
t207 = t120 * t104;
t51 = t173 * t118;
t179 = -t51 + t207;
t294 = t179 * t123;
t75 = t168 * t123;
t264 = t119 * t122;
t249 = pkin(8) * t264;
t167 = -t212 + t249;
t76 = t167 * qJD(2);
t293 = (t121 ^ 2 - t123 ^ 2) * qJD(4);
t292 = -0.2e1 * t79;
t239 = t285 * t118;
t88 = t284 * t120 + t239;
t80 = t88 * qJD(6);
t291 = 0.2e1 * t80;
t115 = t120 ^ 2;
t290 = t122 ^ 2;
t289 = 0.2e1 * t119;
t288 = 0.2e1 * qJD(3);
t283 = pkin(4) * t121;
t282 = pkin(4) * t123;
t281 = t120 * pkin(5);
t280 = pkin(10) + qJ(5);
t152 = t121 * t156;
t174 = -t123 * qJ(3) + t121 * t286;
t169 = qJ(5) + t174;
t163 = t169 * t122;
t178 = t123 * t199;
t132 = -t152 + (t178 + t163) * t119;
t224 = t270 * qJ(3);
t162 = t224 + t213;
t186 = -t121 * t263 + t270 * t123;
t81 = t270 * t121 + t123 * t263;
t141 = t81 * pkin(4) - qJ(5) * t186 + t162;
t135 = t215 + t141;
t18 = t118 * t135 + t120 * t132;
t279 = t121 * t80;
t248 = pkin(8) * t263;
t83 = t213 + t248;
t77 = t83 * qJD(2);
t278 = t122 * t77;
t277 = t123 * t80;
t20 = t299 * t121 + t298 * t123;
t217 = pkin(4) * t104;
t16 = -t20 - t217;
t276 = t16 * t120;
t275 = t16 * t123;
t56 = t118 * t264 + t120 * t186;
t274 = t56 * t118;
t57 = -t121 * t235 + (t222 - t233) * t123;
t273 = t57 * t121;
t60 = t224 + t160;
t272 = t60 * t121;
t260 = t121 * t287;
t102 = t120 * t260;
t190 = -qJ(5) * t123 + t283;
t177 = qJ(3) + t190;
t63 = t118 * t177 - t102;
t267 = qJ(5) * t121;
t266 = t118 * t121;
t265 = t118 * t123;
t262 = t120 * t121;
t261 = t120 * t123;
t113 = t118 ^ 2;
t258 = t113 + t115;
t253 = qJD(3) * t123;
t111 = qJD(4) * t121;
t252 = qJD(4) * t123;
t251 = qJD(4) * t287;
t250 = qJD(5) * t121;
t36 = 0.2e1 * t81 * t57;
t245 = t123 * t286;
t244 = t88 * t252;
t243 = t118 * t260;
t242 = t287 * t264;
t114 = t119 ^ 2;
t236 = t114 * t255;
t108 = t121 * t251;
t234 = t118 * t111;
t232 = t120 * t111;
t231 = t121 * t252;
t230 = qJD(5) * t284;
t227 = t77 * t270;
t226 = t285 * qJD(5);
t225 = t118 * t287 + pkin(5);
t223 = qJD(2) * t270;
t221 = t258 * t123;
t220 = t270 * qJD(3);
t219 = t258 * qJD(5);
t218 = -0.2e1 * t231;
t103 = 0.2e1 * t231;
t216 = -qJD(5) * t123 + qJD(3);
t214 = -t267 - t286;
t209 = t287 * t104;
t208 = t122 * t236;
t206 = t118 * t232;
t205 = t285 * t111;
t204 = t284 * t111;
t201 = 0.2e1 * t219;
t200 = 0.2e1 * t293;
t198 = -t245 - pkin(4);
t197 = t280 * t284;
t70 = t186 * t118;
t196 = t120 * t264 - t70;
t195 = t122 * t223;
t194 = t124 * t223;
t142 = t57 * pkin(4) - t271 * qJ(5) - t186 * qJD(5) + t220;
t134 = pkin(1) * t194 + t142;
t151 = t123 * t156;
t136 = pkin(1) * t121 * t195 - qJD(4) * t151;
t137 = -t199 * t111 + (-t253 + qJD(5) + (qJ(3) * t121 + t245) * qJD(4)) * t122;
t164 = t124 * t169;
t7 = -t118 * t136 + t120 * t134 + (-t118 * t137 + (-t118 * t164 + (t120 * t214 - t301) * t122) * qJD(2)) * t119;
t131 = t118 * t134 + t120 * t136;
t133 = t120 * t137;
t155 = (t118 * t214 + t241) * t122;
t8 = (t133 + (t120 * t164 + t155) * qJD(2)) * t119 + t131;
t193 = -t7 * t118 + t8 * t120;
t192 = -pkin(2) * t124 - t269;
t191 = t267 + t282;
t17 = t118 * t152 + t120 * t141 + (pkin(1) * t265 + (t120 * t286 + t287 * t265) * t124 - t118 * t163) * t119;
t189 = -t118 * t17 + t120 * t18;
t183 = t120 * t216;
t48 = t183 + (t120 * t191 + t301) * qJD(4);
t184 = t118 * t216;
t49 = t184 + (t118 * t191 - t241) * qJD(4);
t188 = -t118 * t48 + t120 * t49;
t62 = t120 * t177 + t243;
t187 = -t118 * t62 + t120 * t63;
t185 = -qJ(5) * t57 - qJD(5) * t81;
t182 = t280 * t239;
t181 = t119 * t195;
t180 = t119 * t194;
t176 = t120 * t196;
t31 = t81 * t252 + t273;
t175 = t280 * t121 + t282;
t172 = t179 * t120;
t171 = t179 * t121;
t170 = t285 * t196;
t161 = t176 + t274;
t27 = t284 * t196 + t285 * t56;
t19 = -t298 * t121 + t299 * t123;
t55 = t121 * t166;
t29 = t123 * t145 - t55;
t30 = -t152 + (t122 * t174 + t178) * t119;
t157 = -t19 * t121 + t20 * t123 + (-t29 * t121 + t30 * t123) * qJD(4);
t150 = t225 * t121 + (-t280 * t123 + qJ(3) + t283) * t120;
t149 = t123 * t153;
t146 = t55 + t151;
t144 = t285 * t150;
t143 = t284 * t150;
t139 = t184 + (t118 * t175 - t241) * qJD(4);
t138 = t183 + (t175 * t120 + t225 * t123) * qJD(4);
t130 = t81 * pkin(5) - t56 * pkin(10) - t118 * t132 + t120 * t135;
t129 = t285 * t130;
t128 = t284 * t130;
t127 = -t51 * pkin(10) + (t133 + ((pkin(10) + t169) * t120 * t124 + t155) * qJD(2)) * t119 + t131;
t126 = -t118 * (qJ(5) * t104 + qJD(5) * t264 - t19) + t120 * ((t214 * t264 + t212) * qJD(2) + t142) - t153 * pkin(10) + t57 * pkin(5);
t112 = qJ(3) * t288;
t110 = -pkin(4) - t281;
t99 = t280 * t120;
t90 = -0.2e1 * t208;
t89 = 0.2e1 * t208;
t86 = (pkin(5) * t118 + t287) * t123;
t78 = -pkin(5) * t234 - t108;
t73 = t88 * t123;
t72 = t88 * t121;
t71 = 0.2e1 * (t124 ^ 2 - t290) * t114 * qJD(2);
t69 = t165 + t249;
t68 = (-pkin(1) + t192) * t119;
t67 = -t162 - t248;
t66 = (-t122 * t111 + t123 * t255) * t119;
t65 = (-t121 * t255 - t122 * t252) * t119;
t64 = -t220 + t76;
t61 = (-t254 + (pkin(2) * t122 - t268) * qJD(2)) * t119;
t59 = -t118 * t197 + t285 * t99;
t58 = -t284 * t99 - t182;
t54 = -pkin(10) * t265 + t63;
t53 = t220 + (-t286 * t264 + t212) * qJD(2);
t43 = -t118 * t205 - t120 * t204 - t297;
t42 = -qJD(6) * t74 - t244;
t41 = -t118 * t204 + t120 * t205 + t277;
t40 = -qJD(4) * t75 + t279;
t38 = -t99 * t229 - t120 * t230 + (qJD(6) * t197 - t226) * t118;
t37 = qJD(6) * t182 + t118 * t230 - t120 * t226 + t228 * t99;
t26 = t284 * t56 - t170;
t25 = t198 * t264 + t146;
t24 = t285 * t54 + t143;
t23 = -t284 * t54 + t144;
t21 = t70 * pkin(5) + (t198 - t281) * t264 + t146;
t15 = t27 * qJD(6) + t284 * t153 - t285 * t179;
t14 = -qJD(6) * t170 - t285 * t153 - t284 * t179 + t56 * t228;
t13 = -qJD(6) * t143 + t285 * t138 - t284 * t139 - t54 * t229;
t12 = -qJD(6) * t144 - t284 * t138 - t285 * t139 + t54 * t228;
t11 = pkin(10) * t196 + t18;
t10 = -pkin(5) * t179 + t16;
t4 = t285 * t11 + t128;
t3 = -t284 * t11 + t129;
t2 = -qJD(6) * t128 - t11 * t229 + t285 * t126 - t284 * t127;
t1 = -qJD(6) * t129 + t11 * t228 - t284 * t126 - t285 * t127;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t71, 0.2e1 * t180, t90, -0.2e1 * t181, 0, -0.2e1 * pkin(1) * t114 * t256 - 0.2e1 * t227, -0.2e1 * pkin(1) * t236 + 0.2e1 * t76 * t270 (t278 - t124 * t76 + (-t122 * t83 + t124 * t167) * qJD(2)) * t289, 0.2e1 * t167 * t77 - 0.2e1 * t76 * t83, 0, -0.2e1 * t180, 0.2e1 * t181, t89, t71, t90 (t278 - t124 * t64 + (t122 * t67 + t124 * t69) * qJD(2)) * t289, 0.2e1 * t227 + 0.2e1 * (t124 * t61 - t256 * t68) * t119, -0.2e1 * t64 * t270 + 0.2e1 * (-t122 * t61 - t255 * t68) * t119, 0.2e1 * t61 * t68 + 0.2e1 * t64 * t67 + 0.2e1 * t69 * t77, 0.2e1 * t186 * t173, -0.2e1 * t173 * t81 - 0.2e1 * t186 * t57 (t271 * t122 + (t121 * t290 * t119 + t124 * t186) * qJD(2)) * t289, t36 (-t122 * t57 - t255 * t81) * t289, t89, 0.2e1 * t53 * t81 + 0.2e1 * t57 * t60 + 0.2e1 * (t122 * t20 + t255 * t29) * t119, 0.2e1 * t53 * t186 + 0.2e1 * t60 * t271 + 0.2e1 * (t19 * t122 + (t122 * t272 - t124 * t30) * qJD(2)) * t119, -0.2e1 * t173 * t29 - 0.2e1 * t186 * t20 + 0.2e1 * t19 * t81 - 0.2e1 * t30 * t57, -0.2e1 * t19 * t30 + 0.2e1 * t20 * t29 + 0.2e1 * t53 * t60, 0.2e1 * t56 * t153, 0.2e1 * t153 * t196 + 0.2e1 * t179 * t56, 0.2e1 * t153 * t81 + 0.2e1 * t56 * t57, 0.2e1 * t196 * t179, 0.2e1 * t179 * t81 + 0.2e1 * t196 * t57, t36, -0.2e1 * t16 * t196 + 0.2e1 * t17 * t57 - 0.2e1 * t179 * t25 + 0.2e1 * t7 * t81, 0.2e1 * t153 * t25 + 0.2e1 * t16 * t56 - 0.2e1 * t18 * t57 - 0.2e1 * t8 * t81, -0.2e1 * t153 * t17 + 0.2e1 * t179 * t18 + 0.2e1 * t196 * t8 - 0.2e1 * t7 * t56, 0.2e1 * t16 * t25 + 0.2e1 * t17 * t7 + 0.2e1 * t18 * t8, -0.2e1 * t27 * t14, 0.2e1 * t14 * t26 - 0.2e1 * t15 * t27, -0.2e1 * t14 * t81 + 0.2e1 * t27 * t57, 0.2e1 * t26 * t15, -0.2e1 * t15 * t81 - 0.2e1 * t26 * t57, t36, 0.2e1 * t10 * t26 + 0.2e1 * t15 * t21 + 0.2e1 * t2 * t81 + 0.2e1 * t3 * t57, 0.2e1 * t1 * t81 + 0.2e1 * t10 * t27 - 0.2e1 * t14 * t21 - 0.2e1 * t4 * t57, 0.2e1 * t1 * t26 + 0.2e1 * t14 * t3 - 0.2e1 * t15 * t4 - 0.2e1 * t2 * t27, -0.2e1 * t1 * t4 + 0.2e1 * t10 * t21 + 0.2e1 * t2 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, -t233, 0, -t77, t76, 0, 0, 0, -t104, t233, 0, 0, 0 (qJD(2) * t192 + qJD(3) * t124) * t119, t77, 0.2e1 * t220 - t76, -pkin(2) * t77 - qJ(3) * t64 - qJD(3) * t67, -t111 * t186 + t295, -t123 * t57 - t173 * t121 + (t121 * t81 - t123 * t186) * qJD(4), t66, t31, t65, 0, -t123 * t209 + qJ(3) * t57 + qJD(3) * t81 + t53 * t121 + (t121 * t242 + t123 * t60) * qJD(4), t121 * t209 + qJD(3) * t186 + qJ(3) * t173 + t53 * t123 + (t123 * t242 - t272) * qJD(4) (t287 * t57 + t19) * t121 + (t173 * t287 - t20) * t123 + ((t287 * t81 - t30) * t123 + (-t186 * t287 + t29) * t121) * qJD(4), t53 * qJ(3) + t60 * qJD(3) - t157 * t287 (-t111 * t56 + t149) * t120 (t172 - t296) * t123 + (-t176 + t274) * t111, t57 * t261 + t153 * t121 + (t123 * t56 - t262 * t81) * qJD(4) (t111 * t196 - t294) * t118, -t57 * t265 + t171 + (t123 * t196 + t266 * t81) * qJD(4), t31, t7 * t121 + t48 * t81 + t62 * t57 + (t16 * t118 - t179 * t287) * t123 + (t17 * t123 + (-t25 * t118 + t196 * t287) * t121) * qJD(4), -t8 * t121 - t49 * t81 - t63 * t57 + (t153 * t287 + t276) * t123 + (-t18 * t123 + (-t120 * t25 - t287 * t56) * t121) * qJD(4), t49 * t196 + t63 * t179 - t48 * t56 - t62 * t153 + (-t8 * t118 - t7 * t120) * t123 + (t118 * t18 + t120 * t17) * t111, t17 * t48 + t18 * t49 + t7 * t62 + t8 * t63 - (t111 * t25 - t275) * t287, -t14 * t75 - t27 * t41, t14 * t73 - t15 * t75 + t26 * t41 - t27 * t43, -t121 * t14 + t252 * t27 - t41 * t81 + t57 * t75, t15 * t73 + t26 * t43, -t121 * t15 - t252 * t26 - t43 * t81 - t57 * t73, t31, t10 * t73 + t121 * t2 + t13 * t81 + t15 * t86 + t21 * t43 + t23 * t57 + t252 * t3 + t26 * t78, t1 * t121 + t10 * t75 + t12 * t81 - t14 * t86 - t21 * t41 - t24 * t57 - t252 * t4 + t27 * t78, t1 * t73 + t12 * t26 - t13 * t27 + t14 * t23 - t15 * t24 - t2 * t75 + t3 * t41 - t4 * t43, -t1 * t24 + t10 * t86 - t12 * t4 + t13 * t3 + t2 * t23 + t21 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, t112, t218, t200, 0, t103, 0, 0, 0.2e1 * qJ(3) * t252 + 0.2e1 * qJD(3) * t121, -0.2e1 * qJ(3) * t111 + 0.2e1 * t253, 0, t112, t115 * t218, 0.4e1 * t123 * t206, -0.2e1 * t120 * t293, t113 * t218, t118 * t200, t103, 0.2e1 * t48 * t121 + 0.2e1 * (t62 - 0.2e1 * t243) * t252, -0.2e1 * t49 * t121 + 0.2e1 * (-t63 - 0.2e1 * t102) * t252, 0.2e1 * (-t118 * t49 - t120 * t48) * t123 + 0.2e1 * (t118 * t63 + t120 * t62) * t111, -0.2e1 * t231 * t287 ^ 2 + 0.2e1 * t62 * t48 + 0.2e1 * t63 * t49, -0.2e1 * t75 * t41, 0.2e1 * t41 * t73 - 0.2e1 * t43 * t75, -0.2e1 * t121 * t41 + 0.2e1 * t252 * t75, 0.2e1 * t73 * t43, -0.2e1 * t121 * t43 - 0.2e1 * t252 * t73, t103, 0.2e1 * t121 * t13 + 0.2e1 * t23 * t252 + 0.2e1 * t43 * t86 + 0.2e1 * t73 * t78, 0.2e1 * t12 * t121 - 0.2e1 * t24 * t252 - 0.2e1 * t41 * t86 + 0.2e1 * t75 * t78, 0.2e1 * t12 * t73 - 0.2e1 * t13 * t75 + 0.2e1 * t23 * t41 - 0.2e1 * t24 * t43, -0.2e1 * t12 * t24 + 0.2e1 * t13 * t23 + 0.2e1 * t78 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0, t77, 0, 0, 0, 0, 0, 0, t66, t65, -t273 - t295 + (t121 * t186 - t81 * t123) * qJD(4), t157, 0, 0, 0, 0, 0, 0, -t57 * t266 + t294 + (-t121 * t196 - t265 * t81) * qJD(4), -t57 * t262 - t149 + (t121 * t56 - t261 * t81) * qJD(4), t120 * t171 + t121 * t296 + t161 * t252, -t275 + t193 * t121 + (t121 * t25 + t123 * t189) * qJD(4), 0, 0, 0, 0, 0, 0, t111 * t26 - t123 * t15 + t42 * t81 - t57 * t72, t111 * t27 + t123 * t14 + t40 * t81 - t57 * t74, -t14 * t72 - t15 * t74 + t26 * t40 - t27 * t42, -t1 * t74 - t10 * t123 + t111 * t21 - t2 * t72 + t3 * t42 - t4 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188 * t121 + (t187 + 0.2e1 * t260) * t252, 0, 0, 0, 0, 0, 0, t121 * t42 - t123 * t43 + (t121 * t73 - t123 * t72) * qJD(4), t121 * t40 + t123 * t41 + (t121 * t75 - t123 * t74) * qJD(4), t40 * t73 - t41 * t72 - t42 * t75 - t43 * t74, t111 * t86 - t12 * t74 - t123 * t78 - t13 * t72 + t23 * t42 - t24 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t258) * t103, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t40 * t74 - 0.2e1 * t42 * t72 - 0.2e1 * t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, 0, -t57, t104, t20, t19, 0, 0, t296, t173 * t115 + (-t51 + 0.2e1 * t207) * t118, t118 * t57, t172, t120 * t57, 0, pkin(4) * t179 + t118 * t185 - t276 (t16 - t217) * t118 + (-pkin(4) * t173 + t185) * t120, t161 * qJD(5) + (t172 + t296) * qJ(5) + t193, -pkin(4) * t16 + qJ(5) * t193 + qJD(5) * t189, -t14 * t88 - t27 * t79, -t14 * t168 - t15 * t88 + t26 * t79 - t27 * t80, t57 * t88 - t79 * t81, -t15 * t168 + t26 * t80, t168 * t57 - t80 * t81, 0, -t10 * t168 + t110 * t15 + t21 * t80 + t38 * t81 + t57 * t58, t10 * t88 - t110 * t14 - t21 * t79 + t37 * t81 - t57 * t59, -t1 * t168 + t14 * t58 - t15 * t59 - t2 * t88 + t26 * t37 - t27 * t38 + t3 * t79 - t4 * t80, -t1 * t59 + t10 * t110 + t2 * t58 + t3 * t38 - t37 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, 0, -t252, 0, t108, t123 * t251, 0, 0, -t206 (t113 - t115) * t111, t118 * t252, t206, t120 * t252, 0, -t118 * t250 + (t118 * t190 + t102) * qJD(4), -t120 * t250 + (t120 * t190 - t243) * qJD(4), t188, pkin(4) * t108 + qJ(5) * t188 + qJD(5) * t187, -t41 * t88 - t75 * t79, -t168 * t41 - t43 * t88 + t73 * t79 - t75 * t80, -t121 * t79 + t244, -t168 * t43 + t73 * t80, t168 * t252 - t279, 0, t110 * t43 + t121 * t38 - t168 * t78 + t252 * t58 + t80 * t86, -t110 * t41 + t121 * t37 - t252 * t59 + t78 * t88 - t79 * t86, -t12 * t168 - t13 * t88 + t23 * t79 - t24 * t80 + t37 * t73 - t38 * t75 + t41 * t58 - t43 * t59, t110 * t78 - t12 * t59 + t13 * t58 + t23 * t38 - t24 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t252, 0, 0, 0, 0, 0, 0, 0, 0, -t232, t234, qJD(4) * t221, t121 * t219 + (qJ(5) * t221 - t283) * qJD(4), 0, 0, 0, 0, 0, 0, -t111 * t168 - t277, t111 * t88 + t297, -t168 * t40 - t42 * t88 - t72 * t79 - t74 * t80, t110 * t111 - t37 * t74 - t38 * t72 - t40 * t59 + t42 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, qJ(5) * t201, t88 * t292, -0.2e1 * t168 * t79 - 0.2e1 * t80 * t88, 0, -t168 * t291, 0, 0, t110 * t291, t110 * t292, -0.2e1 * t168 * t37 - 0.2e1 * t38 * t88 + 0.2e1 * t58 * t79 - 0.2e1 * t59 * t80, -0.2e1 * t37 * t59 + 0.2e1 * t38 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, t153, 0, t16, 0, 0, 0, 0, 0, 0, t15, -t14, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, -t232, 0, -t108, 0, 0, 0, 0, 0, 0, t43, -t41, 0, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, -t15, t57, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, -t43, t252, t13, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, 0, -t80, 0, t38, t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
