% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_inertiaDJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:00:01
% EndTime: 2019-03-09 01:00:22
% DurationCPUTime: 8.19s
% Computational Cost: add. (11798->525), mult. (32717->931), div. (0->0), fcn. (33740->14), ass. (0->228)
t148 = sin(qJ(3));
t152 = cos(qJ(3));
t260 = cos(pkin(7));
t219 = t148 * t260;
t142 = sin(pkin(7));
t258 = t142 * t152;
t105 = pkin(2) * t219 + pkin(9) * t258;
t171 = pkin(10) * t260 + t105;
t251 = qJD(3) * t142;
t295 = qJD(4) * t171 - (pkin(3) * t148 - pkin(10) * t152) * t251;
t218 = t152 * t260;
t259 = t142 * t148;
t180 = -pkin(2) * t218 + pkin(9) * t259;
t174 = t180 * qJD(3);
t191 = (-pkin(3) * t152 - pkin(10) * t148 - pkin(2)) * t142;
t294 = -qJD(4) * t191 + t174;
t225 = t152 * t251;
t293 = qJD(4) * t260 + t225;
t145 = sin(qJ(6));
t140 = t145 ^ 2;
t150 = cos(qJ(6));
t141 = t150 ^ 2;
t253 = t140 - t141;
t290 = t253 * qJD(6);
t147 = sin(qJ(4));
t151 = cos(qJ(4));
t102 = t147 * t260 + t151 * t259;
t146 = sin(qJ(5));
t226 = qJD(4) * t259;
t233 = t147 * t293 + t151 * t226;
t246 = qJD(5) * t146;
t284 = cos(qJ(5));
t86 = t147 * t226 - t151 * t293;
t175 = t102 * t246 + t146 * t233 + t284 * t86;
t199 = t147 * t259 - t151 * t260;
t184 = t284 * t199;
t165 = qJD(5) * t184 + t175;
t74 = t102 * t284 - t146 * t199;
t45 = qJD(5) * t74 - t146 * t86 + t233 * t284;
t98 = qJD(3) * t105;
t70 = pkin(4) * t233 + t98;
t156 = t45 * pkin(5) + pkin(12) * t165 + t70;
t244 = qJD(6) * t145;
t250 = qJD(3) * t148;
t129 = t142 * t250;
t215 = pkin(4) * t129;
t47 = t147 * t294 - t151 * t295;
t155 = t86 * pkin(11) + t215 + t47;
t46 = t147 * t295 + t151 * t294;
t163 = pkin(11) * t233 + t46;
t71 = -t147 * t171 + t151 * t191;
t166 = -t102 * pkin(11) + t71;
t164 = pkin(4) * t258 - t166;
t72 = t147 * t191 + t151 * t171;
t178 = -pkin(11) * t199 + t72;
t36 = -t146 * t178 - t164 * t284;
t14 = -qJD(5) * t36 - t146 * t155 + t284 * t163;
t73 = t146 * t102 + t184;
t96 = -pkin(3) * t260 + t180;
t76 = pkin(4) * t199 + t96;
t162 = t73 * pkin(5) - t74 * pkin(12) + t76;
t292 = -pkin(12) * t129 - qJD(6) * t162 + t14;
t53 = t284 * t178;
t278 = -t146 * t164 + t53;
t35 = -pkin(12) * t258 + t278;
t3 = -t145 * t156 + t150 * t292 + t35 * t244;
t2 = t3 * t150;
t19 = -t145 * t35 + t150 * t162;
t20 = t145 * t162 + t150 * t35;
t203 = t145 * t20 + t150 * t19;
t138 = qJD(6) * t150;
t4 = -t35 * t138 + t145 * t292 + t150 * t156;
t172 = -qJD(6) * t203 - t4 * t145 - t2;
t159 = t146 * t163 + t155 * t284;
t15 = -qJD(5) * t278 + t159;
t13 = -pkin(5) * t129 - t15;
t34 = pkin(5) * t258 - t36;
t220 = -t13 * t150 + t244 * t34;
t149 = sin(qJ(2));
t153 = cos(qJ(2));
t217 = t153 * t260;
t291 = -t148 * t149 + t152 * t217;
t289 = qJD(4) + qJD(5);
t288 = 0.2e1 * t142;
t287 = 0.2e1 * t152;
t286 = 0.2e1 * qJD(3);
t285 = -pkin(11) - pkin(10);
t283 = pkin(10) * t142;
t143 = sin(pkin(6));
t252 = qJD(2) * t143;
t224 = t149 * t252;
t212 = t142 * t224;
t144 = cos(pkin(6));
t62 = t144 * t225 + (t291 * qJD(3) + (-t149 * t219 + t152 * t153) * qJD(2)) * t143;
t179 = -t147 * t62 + t151 * t212;
t101 = -t142 * t143 * t153 + t144 * t260;
t186 = t148 * t217 + t149 * t152;
t81 = t143 * t186 + t144 * t259;
t65 = t101 * t147 + t151 * t81;
t167 = -qJD(4) * t65 + t179;
t64 = t101 * t151 - t147 * t81;
t40 = qJD(4) * t64 + t147 * t212 + t151 * t62;
t43 = t146 * t64 + t284 * t65;
t18 = qJD(5) * t43 + t146 * t40 - t167 * t284;
t42 = t146 * t65 - t284 * t64;
t282 = t18 * t42;
t232 = qJD(4) * t285;
t111 = t147 * t232;
t228 = t284 * t151;
t209 = qJD(4) * t228;
t123 = t285 * t151;
t257 = t146 * t147;
t88 = -t123 * t284 + t257 * t285;
t60 = qJD(5) * t88 + t146 * t111 - t209 * t285;
t229 = t284 * t147;
t87 = -t146 * t123 - t229 * t285;
t281 = t60 * t87;
t61 = t144 * t129 + (t186 * qJD(3) + (t148 * t153 + t149 * t218) * qJD(2)) * t143;
t80 = -t143 * t291 - t144 * t258;
t48 = t80 * t61;
t280 = t80 * t98;
t279 = t13 * t145 + t138 * t34;
t277 = t138 * t87 + t145 * t60;
t276 = pkin(4) * qJD(5);
t256 = t146 * t151;
t109 = t229 + t256;
t221 = qJD(5) * t284;
t84 = -t151 * t221 + t257 * t289 - t209;
t275 = t109 * t84;
t66 = t145 * t74 + t150 * t258;
t26 = qJD(6) * t66 - t129 * t145 + t150 * t165;
t273 = t145 * t26;
t234 = t145 * t258;
t27 = -qJD(6) * t234 - t129 * t150 + t138 * t74 - t145 * t165;
t272 = t145 * t27;
t271 = t145 * t66;
t270 = t145 * t84;
t85 = t289 * t109;
t269 = t145 * t85;
t268 = t146 * t42;
t267 = t146 * t87;
t266 = t147 * t86;
t265 = t150 * t26;
t264 = t150 * t27;
t67 = t150 * t74 - t234;
t263 = t150 * t67;
t262 = t150 * t84;
t261 = t150 * t85;
t136 = -pkin(4) * t284 - pkin(5);
t237 = pkin(4) * t246;
t254 = t136 * t138 + t145 * t237;
t249 = qJD(4) * t147;
t248 = qJD(4) * t151;
t247 = qJD(4) * t152;
t243 = 0.2e1 * t73 * t45;
t242 = -0.2e1 * pkin(3) * qJD(4);
t108 = -t228 + t257;
t241 = 0.2e1 * t108 * t85;
t240 = pkin(4) * t249;
t239 = pkin(5) * t244;
t238 = pkin(5) * t138;
t236 = t145 * t262;
t235 = t80 * t249;
t82 = t87 * t244;
t137 = -pkin(4) * t151 - pkin(3);
t231 = t145 * t284;
t230 = t150 * t284;
t139 = t142 ^ 2;
t227 = t139 * t250;
t223 = t145 * t138;
t222 = t147 * t248;
t214 = pkin(4) * t221;
t213 = t139 * t224;
t106 = t109 ^ 2;
t211 = t106 * t223;
t210 = t152 * t227;
t207 = t233 * t151;
t206 = t260 * t251;
t205 = t18 * t87 + t42 * t60;
t204 = t108 * t45 + t73 * t85;
t28 = -t145 * t43 + t150 * t80;
t29 = t145 * t80 + t150 * t43;
t202 = t145 * t29 + t150 * t28;
t185 = -pkin(5) * t108 + pkin(12) * t109 - t137;
t177 = t150 * t185;
t54 = -t145 * t88 - t177;
t55 = -t145 * t185 + t150 * t88;
t201 = t145 * t55 + t150 * t54;
t200 = t145 * t67 + t150 * t66;
t135 = pkin(4) * t146 + pkin(12);
t198 = t108 * t135 - t109 * t136;
t197 = t136 * t244 - t150 * t237;
t8 = t138 * t42 + t145 * t18;
t9 = -t150 * t18 + t244 * t42;
t31 = t138 * t73 + t145 * t45;
t196 = -t150 * t45 + t244 * t73;
t195 = t108 * t244 - t261;
t194 = -t109 * t138 + t270;
t193 = t109 * t244 + t262;
t192 = t199 * t147;
t189 = (t140 + t141) * t284;
t188 = t147 * t247 + t151 * t250;
t187 = t147 * t250 - t151 * t247;
t183 = pkin(5) * t85 + pkin(12) * t84 + t240;
t182 = (-t108 * t284 + t109 * t146) * qJD(5);
t17 = qJD(5) * t42 - t146 * t167 - t284 * t40;
t5 = -t138 * t80 - t145 * t61 + t150 * t17 + t244 * t43;
t6 = -qJD(6) * t29 + t145 * t17 + t150 * t61;
t1 = -qJD(6) * t202 - t145 * t6 - t150 * t5;
t59 = qJD(5) * t87 - t111 * t284 - t232 * t256;
t24 = qJD(6) * t177 - t145 * t183 + t150 * t59 + t244 * t88;
t25 = -qJD(6) * t55 + t145 * t59 + t150 * t183;
t7 = -qJD(6) * t201 - t145 * t25 - t150 * t24;
t170 = -t47 * t147 - t46 * t151 + (-t72 * t147 - t71 * t151) * qJD(4);
t168 = pkin(4) * t182 - t135 * t85 - t136 * t84;
t158 = -t147 * t179 + t40 * t151 - t248 * t64;
t128 = -0.2e1 * t223;
t127 = 0.2e1 * t223;
t113 = -0.2e1 * t210;
t107 = -0.2e1 * t290;
t100 = t189 * t276;
t69 = t108 * t138 + t269;
t52 = t109 * t290 + t236;
t49 = -0.4e1 * t109 * t223 + t253 * t84;
t23 = t244 * t66 - t264;
t22 = t138 * t67 - t273;
t10 = -qJD(6) * t200 - t265 - t272;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t101 * t212 + 0.2e1 * t62 * t81 + 0.2e1 * t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t167 * t64 + 0.2e1 * t65 * t40 + 0.2e1 * t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17 * t43 + 0.2e1 * t282 + 0.2e1 * t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t28 * t6 - 0.2e1 * t29 * t5 + 0.2e1 * t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224, -t153 * t252, 0, 0, 0, 0, 0, 0, 0, 0, t101 * t129 - t152 * t213 - t260 * t61, t101 * t225 + t148 * t213 - t260 * t62 (t148 * t61 + t152 * t62 + (-t148 * t81 + t152 * t80) * qJD(3)) * t142, -pkin(2) * t213 + t62 * t105 - t174 * t81 + t180 * t61 + t280, 0, 0, 0, 0, 0, 0, t129 * t64 - t167 * t258 + t199 * t61 + t233 * t80, t102 * t61 - t80 * t86 + (t152 * t40 - t250 * t65) * t142, -t102 * t167 - t199 * t40 - t233 * t65 + t64 * t86, t167 * t71 + t40 * t72 - t65 * t46 + t64 * t47 + t61 * t96 + t280, 0, 0, 0, 0, 0, 0, t45 * t80 + t61 * t73 + (t152 * t18 - t250 * t42) * t142, t61 * t74 - t80 * t165 + (-t152 * t17 - t250 * t43) * t142, -t165 * t42 + t17 * t73 + t18 * t74 - t43 * t45, -t14 * t43 - t15 * t42 - t17 * t278 - t18 * t36 + t61 * t76 + t70 * t80, 0, 0, 0, 0, 0, 0, t18 * t66 + t27 * t42 + t28 * t45 + t6 * t73, t18 * t67 - t26 * t42 - t29 * t45 + t5 * t73, t26 * t28 - t27 * t29 + t5 * t66 - t6 * t67, t13 * t42 + t18 * t34 + t19 * t6 - t20 * t5 + t28 * t4 - t29 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t210 (-t148 ^ 2 + t152 ^ 2) * t139 * t286, t206 * t287, t113, -0.2e1 * t148 * t206, 0, -0.2e1 * pkin(2) * t227 - 0.2e1 * t260 * t98 (-t139 * pkin(2) * t152 + t180 * t260) * t286, -0.2e1 * t105 * t129 - 0.2e1 * t174 * t258 + 0.2e1 * t180 * t225 + 0.2e1 * t259 * t98, -0.2e1 * t105 * t174 + 0.2e1 * t180 * t98, -0.2e1 * t102 * t86, -0.2e1 * t102 * t233 + 0.2e1 * t199 * t86 (t102 * t250 + t152 * t86) * t288, 0.2e1 * t199 * t233 (t152 * t233 - t199 * t250) * t288, t113, 0.2e1 * t98 * t199 + 0.2e1 * t96 * t233 + 0.2e1 * (-t152 * t47 + t250 * t71) * t142, 0.2e1 * t102 * t98 - 0.2e1 * t86 * t96 + 0.2e1 * (-t152 * t46 - t250 * t72) * t142, -0.2e1 * t102 * t47 + 0.2e1 * t199 * t46 - 0.2e1 * t233 * t72 + 0.2e1 * t71 * t86, -0.2e1 * t46 * t72 + 0.2e1 * t47 * t71 + 0.2e1 * t96 * t98, -0.2e1 * t74 * t165, 0.2e1 * t165 * t73 - 0.2e1 * t45 * t74 (t152 * t165 + t250 * t74) * t288, t243 (t152 * t45 - t250 * t73) * t288, t113, 0.2e1 * t45 * t76 + 0.2e1 * t70 * t73 + 0.2e1 * (-t15 * t152 + t250 * t36) * t142, 0.2e1 * t70 * t74 - 0.2e1 * t76 * t165 + 0.2e1 * (-t14 * t152 - t250 * t278) * t142, 0.2e1 * t14 * t73 - 0.2e1 * t15 * t74 + 0.2e1 * t165 * t36 - 0.2e1 * t278 * t45, -0.2e1 * t14 * t278 + 0.2e1 * t15 * t36 + 0.2e1 * t70 * t76, -0.2e1 * t67 * t26, 0.2e1 * t26 * t66 - 0.2e1 * t27 * t67, -0.2e1 * t26 * t73 + 0.2e1 * t45 * t67, 0.2e1 * t66 * t27, -0.2e1 * t27 * t73 - 0.2e1 * t45 * t66, t243, 0.2e1 * t13 * t66 + 0.2e1 * t19 * t45 + 0.2e1 * t27 * t34 + 0.2e1 * t4 * t73, 0.2e1 * t13 * t67 - 0.2e1 * t20 * t45 - 0.2e1 * t26 * t34 + 0.2e1 * t3 * t73, 0.2e1 * t19 * t26 - 0.2e1 * t20 * t27 + 0.2e1 * t3 * t66 - 0.2e1 * t4 * t67, 0.2e1 * t13 * t34 + 0.2e1 * t19 * t4 - 0.2e1 * t20 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t62, 0, 0, 0, 0, 0, 0, 0, 0, -t151 * t61 + t235, t147 * t61 + t248 * t80, t158, -t61 * pkin(3) + pkin(10) * t158, 0, 0, 0, 0, 0, 0, t108 * t61 + t80 * t85, t109 * t61 - t80 * t84, t108 * t17 + t109 * t18 - t42 * t84 - t43 * t85, pkin(4) * t235 + t137 * t61 - t17 * t88 - t43 * t59 + t205, 0, 0, 0, 0, 0, 0, t108 * t6 + t109 * t8 - t270 * t42 + t28 * t85, t108 * t5 - t109 * t9 - t262 * t42 - t29 * t85, t202 * t84 + (t145 * t5 - t150 * t6 + (t145 * t28 - t150 * t29) * qJD(6)) * t109, -t24 * t29 + t25 * t28 - t5 * t55 + t54 * t6 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225, 0, -t129, 0, -t98, t174, 0, 0, t102 * t248 - t266, -t147 * t233 - t86 * t151 + (-t102 * t147 - t151 * t199) * qJD(4), t187 * t142, qJD(4) * t192 - t207, t188 * t142, 0, -pkin(3) * t233 - t98 * t151 - t187 * t283 + t249 * t96, pkin(3) * t86 + t147 * t98 - t188 * t283 + t248 * t96 (-t266 - t207 + (t151 * t102 + t192) * qJD(4)) * pkin(10) + t170, -pkin(3) * t98 + pkin(10) * t170, -t109 * t165 - t74 * t84, t108 * t165 - t109 * t45 + t73 * t84 - t74 * t85 (t109 * t250 + t152 * t84) * t142, t204 (-t108 * t250 + t152 * t85) * t142, 0, t73 * t240 + t108 * t70 + t137 * t45 + t76 * t85 + (t152 * t60 - t250 * t87) * t142, t74 * t240 - t137 * t165 + t70 * t109 - t76 * t84 + (-t152 * t59 - t250 * t88) * t142, t108 * t14 - t109 * t15 - t165 * t87 - t278 * t85 + t36 * t84 - t45 * t88 + t59 * t73 + t60 * t74, t137 * t70 - t14 * t88 - t15 * t87 + t240 * t76 - t278 * t59 - t36 * t60, -t67 * t262 + (-t244 * t67 - t265) * t109, t200 * t84 + (t273 - t264 + (-t263 + t271) * qJD(6)) * t109, -t108 * t26 - t109 * t196 - t262 * t73 + t67 * t85, -t66 * t270 + (t138 * t66 + t272) * t109, -t108 * t27 - t109 * t31 + t270 * t73 - t66 * t85, t204, t108 * t4 + t109 * t279 + t19 * t85 + t25 * t73 + t27 * t87 - t270 * t34 + t45 * t54 + t60 * t66, t108 * t3 - t109 * t220 - t20 * t85 + t24 * t73 - t26 * t87 - t262 * t34 - t45 * t55 + t60 * t67, t24 * t66 - t25 * t67 + t26 * t54 - t27 * t55 + t203 * t84 + (t145 * t3 - t150 * t4 + (t145 * t19 - t150 * t20) * qJD(6)) * t109, t13 * t87 + t19 * t25 - t20 * t24 - t3 * t55 + t34 * t60 + t4 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t222, 0.2e1 * (-t147 ^ 2 + t151 ^ 2) * qJD(4), 0, -0.2e1 * t222, 0, 0, t147 * t242, t151 * t242, 0, 0, -0.2e1 * t275, 0.2e1 * t108 * t84 - 0.2e1 * t109 * t85, 0, t241, 0, 0, 0.2e1 * t108 * t240 + 0.2e1 * t137 * t85, 0.2e1 * t109 * t240 - 0.2e1 * t137 * t84, 0.2e1 * t108 * t59 + 0.2e1 * t109 * t60 - 0.2e1 * t84 * t87 - 0.2e1 * t85 * t88, 0.2e1 * t137 * t240 - 0.2e1 * t59 * t88 + 0.2e1 * t281, -0.2e1 * t141 * t275 - 0.2e1 * t211, 0.2e1 * t106 * t290 + 0.4e1 * t109 * t236, -0.2e1 * t108 * t193 + 0.2e1 * t109 * t261, -0.2e1 * t140 * t275 + 0.2e1 * t211, 0.2e1 * t108 * t194 - 0.2e1 * t109 * t269, t241, 0.2e1 * t108 * t25 + 0.2e1 * t109 * t277 - 0.2e1 * t270 * t87 + 0.2e1 * t54 * t85, -0.2e1 * t87 * t262 + 0.2e1 * t108 * t24 - 0.2e1 * t55 * t85 + 0.2e1 * (t60 * t150 - t82) * t109, 0.2e1 * t201 * t84 + 0.2e1 * (t145 * t24 - t150 * t25 + (t145 * t54 - t150 * t55) * qJD(6)) * t109, -0.2e1 * t24 * t55 + 0.2e1 * t25 * t54 + 0.2e1 * t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, -t40, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0 (-t284 * t18 - t146 * t17 + (t284 * t43 + t268) * qJD(5)) * pkin(4), 0, 0, 0, 0, 0, 0, t9, t8, t1, t18 * t136 + (t230 * t29 - t231 * t28 + t268) * t276 + t1 * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, 0, -t233, t129, t47, t46, 0, 0, 0, 0, -t165, 0, -t45, t129 (-t146 * t166 - t53) * qJD(5) + (t246 * t287 + t250 * t284) * t142 * pkin(4) + t159, -t146 * t215 + t214 * t258 + t14 (-t146 * t45 + t284 * t175 + (t146 * t74 + (-t73 + t184) * t284) * qJD(5)) * pkin(4) (t284 * t15 - t14 * t146 + (-t146 * t36 + t278 * t284) * qJD(5)) * pkin(4), t22, t10, t31, t23, -t196, 0, t136 * t27 - t31 * t135 + (t146 * t66 - t231 * t73) * t276 + t220, -t136 * t26 + t196 * t135 + (t146 * t67 - t230 * t73) * t276 + t279, -t2 + (-t66 * t214 - t135 * t27 + (t135 * t67 - t19) * qJD(6)) * t150 + (t67 * t214 - t135 * t26 - t4 + (t135 * t66 - t20) * qJD(6)) * t145, t13 * t136 + (t146 * t34 - t19 * t231 + t20 * t230) * t276 + t172 * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, 0, -t249, 0, -pkin(10) * t248, pkin(10) * t249, 0, 0, 0, 0, -t84, 0, -t85, 0, -t60, t59 (-t146 * t85 + t284 * t84 + t182) * pkin(4) (-t284 * t60 - t146 * t59 + (t284 * t88 + t267) * qJD(5)) * pkin(4), -t52, t49, t69, t52, -t195, 0, t82 + (-qJD(6) * t198 - t60) * t150 + t168 * t145, t150 * t168 + t198 * t244 + t277, t7, t60 * t136 + (t230 * t55 - t231 * t54 + t267) * t276 + t7 * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t237, -0.2e1 * t214, 0, 0, t127, t107, 0, t128, 0, 0, 0.2e1 * t197, 0.2e1 * t254, 0.2e1 * t100, 0.2e1 * (t135 * t189 + t136 * t146) * t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0, 0, 0, 0, 0, 0, 0, 0, t9, t8, t1, -pkin(5) * t18 + pkin(12) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, 0, -t45, t129, t15, t14, 0, 0, t22, t10, t31, t23, -t196, 0, -pkin(5) * t27 - pkin(12) * t31 + t220, pkin(5) * t26 + pkin(12) * t196 + t279 (-t273 - t264 + (t263 + t271) * qJD(6)) * pkin(12) + t172, -pkin(5) * t13 + pkin(12) * t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, 0, -t85, 0, -t60, t59, 0, 0, -t52, t49, t69, t52, -t195, 0, t82 + (pkin(5) * t84 - pkin(12) * t85) * t145 + (-t60 + (-pkin(5) * t109 - pkin(12) * t108) * qJD(6)) * t150, pkin(5) * t193 + pkin(12) * t195 + t277, t7, -pkin(5) * t60 + pkin(12) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, -t214, 0, 0, t127, t107, 0, t128, 0, 0, t197 - t239, -t238 + t254, t100 (-pkin(5) * t146 + pkin(12) * t189) * t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t107, 0, t128, 0, 0, -0.2e1 * t239, -0.2e1 * t238, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, -t27, t45, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, 0, t194, t85, t25, t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, 0, -t244, 0, -t135 * t138 - t145 * t214, t135 * t244 - t150 * t214, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, 0, -t244, 0, -pkin(12) * t138, pkin(12) * t244, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t11;
