% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR11_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:42
% EndTime: 2019-12-31 22:44:02
% DurationCPUTime: 7.19s
% Computational Cost: add. (7702->472), mult. (20561->874), div. (0->0), fcn. (19463->10), ass. (0->218)
t111 = sin(pkin(5));
t113 = sin(qJ(3));
t257 = cos(pkin(5));
t200 = t257 * t113;
t114 = sin(qJ(2));
t116 = cos(qJ(3));
t252 = t114 * t116;
t157 = t111 * t252 + t200;
t256 = t111 * t114;
t75 = t113 * t256 - t116 * t257;
t271 = t75 * pkin(3);
t117 = cos(qJ(2));
t224 = pkin(1) * t257;
t189 = t117 * t224;
t153 = -pkin(2) * t257 - t189;
t234 = pkin(7) * t256;
t70 = t153 + t234;
t131 = -pkin(9) * t157 + t271 + t70;
t190 = t114 * t224;
t152 = pkin(8) * t257 + t190;
t255 = t111 * t117;
t233 = pkin(7) * t255;
t142 = t152 + t233;
t139 = qJD(3) * t142;
t177 = pkin(2) * t114 - pkin(8) * t117;
t247 = qJD(2) * t111;
t162 = t177 * t247;
t154 = -t189 + t234;
t149 = t154 * qJD(2);
t227 = -pkin(8) * t114 - pkin(1);
t169 = -pkin(2) * t117 + t227;
t163 = t169 * t111;
t282 = -qJD(3) * t163 + t149;
t21 = t282 * t116 + (t139 - t162) * t113;
t246 = qJD(2) * t114;
t211 = t111 * t246;
t285 = -pkin(9) * t211 - qJD(4) * t131 + t21;
t148 = t116 * t152;
t126 = t148 + (t113 * t227 + (-pkin(2) * t113 + pkin(7) * t116 - pkin(9)) * t117) * t111;
t199 = qJD(2) * t257;
t178 = t114 * t199;
t245 = qJD(2) * t117;
t210 = t111 * t245;
t198 = t257 * qJD(3);
t220 = qJD(3) * t256;
t88 = t113 * t220;
t276 = -t116 * (t198 + t210) + t88;
t150 = t157 * qJD(3);
t55 = t113 * t210 + t150;
t284 = -pkin(1) * t178 - t55 * pkin(3) - pkin(7) * t210 - pkin(9) * t276 + qJD(4) * t126;
t112 = sin(qJ(4));
t115 = cos(qJ(4));
t151 = qJD(4) * t157;
t95 = t115 * t255;
t134 = -qJD(4) * t95 - t115 * t276 + (-t151 + t211) * t112;
t132 = t134 * t112;
t240 = qJD(4) * t115;
t56 = -t112 * t255 + t115 * t157;
t283 = t56 * t240 + t132;
t261 = t112 * t157 + t95;
t201 = t261 * qJD(4);
t241 = qJD(4) * t112;
t214 = t117 * t241;
t262 = t112 * t276 - t115 * t151;
t275 = -t111 * (t115 * t246 + t214) - t262;
t281 = t112 * t201 - t115 * t275;
t270 = cos(qJ(5));
t204 = t270 * qJD(5);
t280 = qJD(4) * t270 + t204;
t243 = qJD(3) * t116;
t158 = t112 * t243 + t113 * t240;
t217 = t115 * t243;
t279 = t113 * t241 - t217;
t278 = -t115 * t55 + t75 * t241;
t105 = qJD(3) * t113;
t43 = t75 * t105 - t116 * t55;
t107 = t112 ^ 2;
t109 = t115 ^ 2;
t249 = t107 - t109;
t197 = qJD(4) * t249;
t277 = qJD(4) + qJD(5);
t239 = qJD(4) * t116;
t215 = t112 * t239;
t244 = qJD(3) * t115;
t159 = t113 * t244 + t215;
t267 = pkin(3) * t116;
t176 = -pkin(9) * t113 - t267;
t168 = -pkin(2) + t176;
t164 = t115 * t168;
t175 = pkin(3) * t113 - pkin(9) * t116;
t165 = t175 * t112;
t44 = pkin(8) * t159 - qJD(3) * t165 - qJD(4) * t164;
t273 = -pkin(10) - pkin(9);
t274 = t113 * t273 - pkin(2) - t267;
t110 = t116 ^ 2;
t272 = t55 * pkin(4);
t269 = sin(qJ(5));
t268 = pkin(3) * t115;
t266 = pkin(8) * t111;
t265 = pkin(8) * t112;
t264 = pkin(8) * t113;
t263 = pkin(8) * t116;
t18 = t112 * t131 + t115 * t126;
t251 = t115 * t116;
t101 = pkin(8) * t251;
t67 = t112 * t168 + t101;
t254 = t112 * t113;
t253 = t113 * t115;
t250 = t116 * t117;
t108 = t113 ^ 2;
t248 = t108 - t110;
t242 = qJD(3) * t117;
t48 = 0.2e1 * t75 * t55;
t238 = -0.2e1 * pkin(2) * qJD(3);
t237 = -0.2e1 * pkin(3) * qJD(4);
t236 = -t113 * t282 + t116 * t139;
t235 = t112 * t263;
t232 = pkin(4) * t241;
t231 = pkin(8) * t243;
t104 = -pkin(4) * t115 - pkin(3);
t225 = -pkin(4) - t265;
t223 = t270 * t115;
t222 = t269 * t112;
t106 = t111 ^ 2;
t221 = t106 * t245;
t218 = t113 * t242;
t212 = t115 * t239;
t209 = t112 * t240;
t208 = t113 * t243;
t207 = qJD(3) * t270;
t206 = qJD(3) * t269;
t203 = t269 * qJD(5);
t202 = t261 * qJD(3);
t196 = t248 * qJD(3);
t195 = 0.2e1 * t208;
t194 = pkin(4) * t204;
t193 = pkin(4) * t203;
t192 = t273 * t270;
t191 = t273 * t269;
t188 = t113 * t222;
t187 = t115 * t208;
t186 = t108 * t209;
t185 = t114 * t221;
t184 = t116 * t207;
t183 = t116 * t206;
t182 = t270 * t261;
t181 = t269 * t261;
t180 = t116 * t202;
t8 = t112 * t284 + t115 * t285;
t9 = t112 * t285 - t115 * t284;
t174 = -t9 * t112 - t8 * t115;
t66 = t164 - t235;
t173 = -t112 * t67 - t115 * t66;
t22 = t116 * t162 - t236;
t172 = -t22 * t113 - t21 * t116;
t171 = t112 * t192;
t170 = t112 * t191;
t20 = (-t114 * pkin(3) - t116 * t177) * t247 + t236;
t46 = -t113 * t142 + t116 * t163;
t42 = pkin(3) * t255 - t46;
t167 = t20 * t112 + t240 * t42;
t166 = t112 * t55 + t240 * t75;
t118 = -pkin(10) * t134 + t272 + t9;
t121 = t75 * pkin(4) - t56 * pkin(10) - t112 * t126 + t115 * t131;
t120 = t270 * t121;
t123 = -pkin(10) * t275 - t8;
t16 = -pkin(10) * t261 + t18;
t1 = -qJD(5) * t120 - t118 * t269 - t123 * t270 + t16 * t203;
t124 = (-t112 * t274 - t101) * qJD(4) + (t273 * t251 + (-t225 + t268) * t113) * qJD(3);
t133 = -pkin(10) * t158 - t44;
t138 = t115 * t274 + t116 * t225;
t137 = t270 * t138;
t57 = -pkin(10) * t254 + t67;
t12 = -qJD(5) * t137 - t124 * t269 - t133 * t270 + t203 * t57;
t85 = t112 * t270 + t115 * t269;
t161 = t113 * t246 - t116 * t242;
t80 = t190 + t233;
t93 = t273 * t115;
t59 = -t270 * t93 + t170;
t147 = pkin(7) * t250 + t113 * t169;
t45 = -t67 * qJD(4) + (pkin(8) * t254 + t115 * t175) * qJD(3);
t141 = qJD(4) * t173 - t112 * t45 - t115 * t44;
t136 = t269 * t138;
t135 = t112 * t275 + t115 * t201;
t54 = t277 * t85;
t31 = t270 * t57 + t136;
t130 = t124 * t270 - t133 * t269;
t128 = t115 * t134 - t241 * t56;
t119 = t269 * t121;
t2 = -qJD(5) * t119 + t118 * t270 - t123 * t269 - t16 * t204;
t97 = -0.2e1 * t208;
t87 = -0.2e1 * t185;
t86 = (pkin(4) * t112 + pkin(8)) * t113;
t84 = t222 - t223;
t74 = t113 * t223 - t188;
t73 = t85 * t113;
t72 = qJD(2) * t80;
t64 = pkin(4) * t158 + t231;
t60 = -t112 * t217 + t113 * t197;
t58 = t269 * t93 + t171;
t53 = (qJD(4) * t269 + t203) * t112 - t280 * t115;
t47 = t111 * t147 + t148;
t39 = -t59 * qJD(5) + (t115 * t192 - t170) * qJD(4);
t38 = -t171 * t277 - t191 * t240 - t203 * t93;
t36 = t112 * t184 + (t113 * t280 + t183) * t115 - t277 * t188;
t35 = t112 * t183 + t113 * t54 - t115 * t184;
t34 = t270 * t56 - t181;
t33 = t269 * t56 + t182;
t30 = -t269 * t57 + t137;
t23 = pkin(4) * t261 + t42;
t17 = -t112 * t148 + t115 * (-pkin(9) * t200 + t153 + t271) + (-t112 * (-t117 * pkin(9) + t147) + t115 * (pkin(7) * t114 - pkin(9) * t252)) * t111;
t14 = -t262 * pkin(4) + (-pkin(4) * t214 + (pkin(8) * t250 + (-pkin(2) * t116 + t104) * t114) * qJD(2)) * t111 + t236;
t13 = -qJD(5) * t31 + t130;
t11 = -qJD(5) * t181 + t134 * t269 + t204 * t56 + t270 * t275;
t10 = qJD(5) * t182 - t134 * t270 + t203 * t56 + t269 * t275;
t7 = t16 * t270 + t119;
t6 = -t16 * t269 + t120;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t185, 0.2e1 * (-t114 ^ 2 + t117 ^ 2) * t106 * qJD(2), 0.2e1 * t199 * t255, t87, -0.2e1 * t111 * t178, 0, -0.2e1 * pkin(1) * t106 * t246 - 0.2e1 * t257 * t72, -0.2e1 * pkin(1) * t221 + 0.2e1 * t149 * t257, -0.2e1 * t149 * t255 + 0.2e1 * t154 * t210 - 0.2e1 * t211 * t80 + 0.2e1 * t256 * t72, -0.2e1 * t149 * t80 + 0.2e1 * t154 * t72, -0.2e1 * t157 * t276, -0.2e1 * t157 * t55 + 0.2e1 * t276 * t75, 0.2e1 * t157 * t211 + 0.2e1 * t255 * t276, t48, 0.2e1 * (t117 * t55 - t246 * t75) * t111, t87, 0.2e1 * t70 * t55 + 0.2e1 * t72 * t75 + 0.2e1 * (-t117 * t22 + t246 * t46) * t111, 0.2e1 * t157 * t72 - 0.2e1 * t21 * t255 - 0.2e1 * t211 * t47 - 0.2e1 * t276 * t70, -0.2e1 * t157 * t22 + 0.2e1 * t21 * t75 + 0.2e1 * t276 * t46 - 0.2e1 * t47 * t55, -0.2e1 * t21 * t47 + 0.2e1 * t22 * t46 + 0.2e1 * t70 * t72, 0.2e1 * t56 * t134, -0.2e1 * t134 * t261 - 0.2e1 * t275 * t56, 0.2e1 * t134 * t75 + 0.2e1 * t56 * t55, 0.2e1 * t261 * t275, -0.2e1 * t261 * t55 - 0.2e1 * t275 * t75, t48, 0.2e1 * t17 * t55 + 0.2e1 * t20 * t261 + 0.2e1 * t275 * t42 + 0.2e1 * t9 * t75, 0.2e1 * t134 * t42 - 0.2e1 * t18 * t55 + 0.2e1 * t20 * t56 + 0.2e1 * t8 * t75, -0.2e1 * t134 * t17 - 0.2e1 * t18 * t275 + 0.2e1 * t261 * t8 - 0.2e1 * t9 * t56, 0.2e1 * t17 * t9 - 0.2e1 * t18 * t8 + 0.2e1 * t20 * t42, -0.2e1 * t34 * t10, 0.2e1 * t10 * t33 - 0.2e1 * t11 * t34, -0.2e1 * t10 * t75 + 0.2e1 * t34 * t55, 0.2e1 * t33 * t11, -0.2e1 * t11 * t75 - 0.2e1 * t33 * t55, t48, 0.2e1 * t11 * t23 + 0.2e1 * t14 * t33 + 0.2e1 * t2 * t75 + 0.2e1 * t55 * t6, 0.2e1 * t1 * t75 - 0.2e1 * t10 * t23 + 0.2e1 * t14 * t34 - 0.2e1 * t55 * t7, 0.2e1 * t1 * t33 + 0.2e1 * t10 * t6 - 0.2e1 * t11 * t7 - 0.2e1 * t2 * t34, -0.2e1 * t1 * t7 + 0.2e1 * t14 * t23 + 0.2e1 * t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, 0, -t211, 0, -t72, t149, 0, 0, t110 * t220 + (-t88 + (0.2e1 * t198 + t210) * t116) * t113, -t116 * t276 - t75 * t243 + (-t150 - t55) * t113, t161 * t111, t43, (t116 * t246 + t218) * t111, 0, -pkin(2) * t55 + t105 * t70 - t72 * t116 - t161 * t266, pkin(2) * t276 + t72 * t113 - t211 * t263 - t218 * t266 + t243 * t70, pkin(8) * t43 - t47 * t105 + t157 * t231 - t46 * t243 - t264 * t276 + t172, -pkin(2) * t72 + ((-t113 * t47 - t116 * t46) * qJD(3) + t172) * pkin(8), t113 * t128 + t217 * t56, -t115 * t180 - t158 * t56 + (-t132 + t281) * t113, t56 * t105 - t116 * t134 + t55 * t253 - t279 * t75, t112 * t180 + t113 * t135, (-t112 * qJD(3) * t75 + t275) * t116 + (-t202 - t166) * t113, t43, t45 * t75 + t66 * t55 + (-t9 + (pkin(8) * t261 + t42 * t112) * qJD(3)) * t116 + (pkin(8) * t275 + t17 * qJD(3) + t167) * t113, -t18 * t105 - t8 * t116 + t134 * t264 + t20 * t253 + t231 * t56 - t279 * t42 + t44 * t75 - t67 * t55, -t66 * t134 - t158 * t18 + t17 * t279 - t9 * t253 + t8 * t254 + t44 * t261 - t275 * t67 - t45 * t56, t17 * t45 - t18 * t44 + t66 * t9 - t67 * t8 + (t20 * t113 + t243 * t42) * pkin(8), -t10 * t74 - t34 * t35, t10 * t73 - t11 * t74 + t33 * t35 - t34 * t36, t10 * t116 + t105 * t34 - t35 * t75 + t55 * t74, t11 * t73 + t33 * t36, -t105 * t33 + t11 * t116 - t36 * t75 - t55 * t73, t43, t105 * t6 + t11 * t86 - t116 * t2 + t13 * t75 + t14 * t73 + t23 * t36 + t30 * t55 + t33 * t64, -t1 * t116 - t10 * t86 - t105 * t7 + t12 * t75 + t14 * t74 - t23 * t35 - t31 * t55 + t34 * t64, t1 * t73 + t10 * t30 - t11 * t31 + t12 * t33 - t13 * t34 - t2 * t74 + t35 * t6 - t36 * t7, -t1 * t31 - t12 * t7 + t13 * t6 + t14 * t86 + t2 * t30 + t23 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, -0.2e1 * t196, 0, t97, 0, 0, t113 * t238, t116 * t238, 0, 0, 0.2e1 * t109 * t208 - 0.2e1 * t186, 0.2e1 * t108 * t197 - 0.4e1 * t112 * t187, 0.2e1 * t113 * t215 + 0.2e1 * t244 * t248, 0.2e1 * t107 * t208 + 0.2e1 * t186, -0.2e1 * t112 * t196 + 0.2e1 * t113 * t212, t97, 0.2e1 * t66 * t105 - 0.2e1 * t116 * t45 + 0.2e1 * (t108 * t240 + t112 * t195) * pkin(8), -0.2e1 * t67 * t105 - 0.2e1 * t116 * t44 + 0.2e1 * (-t108 * t241 + 0.2e1 * t187) * pkin(8), 0.2e1 * t173 * t243 + 0.2e1 * (t112 * t44 - t115 * t45 + (t112 * t66 - t115 * t67) * qJD(4)) * t113, 0.2e1 * pkin(8) ^ 2 * t208 - 0.2e1 * t44 * t67 + 0.2e1 * t45 * t66, -0.2e1 * t74 * t35, 0.2e1 * t35 * t73 - 0.2e1 * t36 * t74, 0.2e1 * t105 * t74 + 0.2e1 * t116 * t35, 0.2e1 * t73 * t36, -0.2e1 * t105 * t73 + 0.2e1 * t116 * t36, t97, 0.2e1 * t105 * t30 - 0.2e1 * t116 * t13 + 0.2e1 * t36 * t86 + 0.2e1 * t64 * t73, -0.2e1 * t105 * t31 - 0.2e1 * t116 * t12 - 0.2e1 * t35 * t86 + 0.2e1 * t64 * t74, 0.2e1 * t12 * t73 - 0.2e1 * t13 * t74 + 0.2e1 * t30 * t35 - 0.2e1 * t31 * t36, -0.2e1 * t12 * t31 + 0.2e1 * t13 * t30 + 0.2e1 * t64 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t276, 0, -t55, t211, t22, t21, 0, 0, t283, t128 - t135, t166, t281, -t278, 0, -pkin(3) * t275 - pkin(9) * t166 - t20 * t115 + t241 * t42, -pkin(3) * t134 + pkin(9) * t278 + t167, -t17 * t240 - t18 * t241 + t174 + (t281 + t283) * pkin(9), -pkin(3) * t20 + ((-t18 * t112 - t17 * t115) * qJD(4) + t174) * pkin(9), -t10 * t85 - t34 * t53, t10 * t84 - t11 * t85 + t33 * t53 - t34 * t54, -t53 * t75 + t55 * t85, t11 * t84 + t33 * t54, -t54 * t75 - t55 * t84, 0, t104 * t11 + t14 * t84 + t23 * t54 + t232 * t33 + t39 * t75 + t55 * t58, -t10 * t104 + t14 * t85 - t23 * t53 + t232 * t34 + t38 * t75 - t55 * t59, t1 * t84 + t10 * t58 - t11 * t59 - t2 * t85 + t33 * t38 - t34 * t39 + t53 * t6 - t54 * t7, -t1 * t59 + t104 * t14 + t2 * t58 + t23 * t232 - t38 * t7 + t39 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, 0, -t105, 0, -t231, pkin(8) * t105, 0, 0, -t60, -0.4e1 * t113 * t209 - t243 * t249, t105 * t112 - t212, t60, t159, 0, (pkin(9) * t251 + (t265 - t268) * t113) * qJD(4) + (t112 * t176 - t101) * qJD(3), (pkin(8) * t253 + t165) * qJD(4) + (t115 * t176 + t235) * qJD(3), t141, -pkin(3) * t231 + pkin(9) * t141, -t35 * t85 - t53 * t74, t35 * t84 - t36 * t85 + t53 * t73 - t54 * t74, t105 * t85 + t116 * t53, t36 * t84 + t54 * t73, -t105 * t84 + t116 * t54, 0, t104 * t36 + t105 * t58 - t116 * t39 + t232 * t73 + t54 * t86 + t64 * t84, -t104 * t35 - t105 * t59 - t116 * t38 + t232 * t74 - t53 * t86 + t64 * t85, t12 * t84 - t13 * t85 + t30 * t53 - t31 * t54 + t35 * t58 - t36 * t59 + t38 * t73 - t39 * t74, t104 * t64 - t12 * t59 + t13 * t58 + t232 * t86 + t30 * t39 - t31 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t209, -0.2e1 * t197, 0, -0.2e1 * t209, 0, 0, t112 * t237, t115 * t237, 0, 0, -0.2e1 * t85 * t53, 0.2e1 * t53 * t84 - 0.2e1 * t54 * t85, 0, 0.2e1 * t84 * t54, 0, 0, 0.2e1 * t104 * t54 + 0.2e1 * t232 * t84, -0.2e1 * t104 * t53 + 0.2e1 * t232 * t85, 0.2e1 * t38 * t84 - 0.2e1 * t39 * t85 + 0.2e1 * t53 * t58 - 0.2e1 * t54 * t59, 0.2e1 * t104 * t232 - 0.2e1 * t38 * t59 + 0.2e1 * t39 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, 0, -t275, t55, t9, t8, 0, 0, 0, 0, -t10, 0, -t11, t55, -t193 * t75 + t270 * t272 + t2, (-t204 * t75 - t269 * t55) * pkin(4) + t1, (t270 * t10 - t269 * t11 + (t269 * t34 - t270 * t33) * qJD(5)) * pkin(4), (t270 * t2 - t269 * t1 + (-t269 * t6 + t270 * t7) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t279, 0, -t158, t105, t45, t44, 0, 0, 0, 0, -t35, 0, -t36, t105, pkin(4) * t113 * t207 - qJD(5) * t136 + t116 * t193 - t204 * t57 + t130, (-t113 * t206 + t116 * t204) * pkin(4) + t12, (t270 * t35 - t269 * t36 + (t269 * t74 - t270 * t73) * qJD(5)) * pkin(4), (t270 * t13 - t269 * t12 + (-t269 * t30 + t270 * t31) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, 0, -t241, 0, -pkin(9) * t240, pkin(9) * t241, 0, 0, 0, 0, -t53, 0, -t54, 0, t39, t38, (t270 * t53 - t269 * t54 + (t269 * t85 - t270 * t84) * qJD(5)) * pkin(4), (t270 * t39 - t269 * t38 + (-t269 * t58 + t270 * t59) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t193, -0.2e1 * t194, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, -t11, t55, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, -t36, t105, t13, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, 0, -t54, 0, t39, t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, -t194, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
