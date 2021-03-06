% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR12_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:55:12
% EndTime: 2019-12-31 22:55:54
% DurationCPUTime: 16.51s
% Computational Cost: add. (13152->572), mult. (40133->1072), div. (0->0), fcn. (39747->12), ass. (0->264)
t102 = sin(pkin(5));
t105 = sin(qJ(3));
t110 = cos(qJ(2));
t291 = cos(pkin(6));
t231 = t110 * t291;
t106 = sin(qJ(2));
t109 = cos(qJ(3));
t283 = t106 * t109;
t175 = t105 * t231 + t283;
t161 = t175 * qJD(3);
t101 = sin(pkin(6));
t292 = cos(pkin(5));
t234 = t101 * t292;
t215 = t105 * t234;
t197 = qJD(3) * t215;
t232 = t109 * t291;
t174 = t105 * t110 + t106 * t232;
t319 = t174 * qJD(2);
t124 = t197 - (-t319 - t161) * t102;
t213 = t102 * t231;
t214 = t109 * t234;
t284 = t105 * t106;
t251 = t102 * t284;
t57 = -t109 * t213 - t214 + t251;
t335 = 0.2e1 * t57 * t124;
t108 = cos(qJ(4));
t104 = sin(qJ(4));
t273 = qJD(4) * t104;
t333 = -t124 * t108 + t57 * t273;
t193 = t102 * (-t291 * pkin(9) - pkin(8));
t248 = pkin(1) * t292;
t221 = t110 * t248;
t142 = t106 * t193 + t221;
t247 = t292 * pkin(2);
t250 = -pkin(2) * t110 - pkin(1);
t289 = t101 * t106;
t122 = t291 * (t247 + t142) + t101 * (-pkin(9) * t289 + t250) * t102;
t163 = t213 + t234;
t222 = t106 * t248;
t285 = t102 * t110;
t171 = -pkin(8) * t285 - t222;
t134 = pkin(9) * t163 - t171;
t30 = -t105 * t134 + t109 * t122;
t162 = t175 * t102;
t287 = t101 * t109;
t330 = t102 * t319;
t331 = -t287 * t330 + ((-t214 + t57) * t105 - t109 * t162) * qJD(3) * t101;
t182 = t101 * t285 - t292 * t291;
t26 = pkin(3) * t182 - t30;
t58 = t215 + t162;
t46 = t104 * t58 + t108 * t182;
t47 = -t104 * t182 + t58 * t108;
t115 = t46 * pkin(4) - t47 * pkin(11) + t26;
t176 = t182 * pkin(10);
t31 = t105 * t122 + t109 * t134;
t117 = -t176 + t31;
t141 = qJD(2) * (t110 * t193 - t222);
t280 = qJD(2) * t102;
t168 = (-pkin(9) * t101 * t110 + pkin(2) * t106) * t280;
t286 = t102 * t106;
t216 = t291 * t286;
t198 = qJD(2) * t216;
t278 = qJD(2) * t110;
t238 = t102 * t278;
t44 = -t105 * t198 - qJD(3) * t251 + (qJD(3) * t163 + t238) * t109;
t323 = -pkin(3) * t124 + t44 * pkin(10) + qJD(4) * t117 + t101 * t141 - t291 * t168;
t217 = t57 * pkin(3) - t58 * pkin(10);
t167 = t221 + t247;
t41 = -t101 * t167 + (pkin(8) * t289 + t291 * t250) * t102;
t126 = t217 + t41;
t317 = -qJD(2) * t142 - t122 * qJD(3);
t318 = qJD(3) * t134 - t101 * t168 - t291 * t141;
t20 = t318 * t105 + t317 * t109;
t279 = qJD(2) * t106;
t239 = t102 * t279;
t218 = t101 * t239;
t324 = -pkin(10) * t218 - qJD(4) * t126 + t20;
t5 = t323 * t104 + t324 * t108;
t325 = -pkin(11) * t124 - qJD(5) * t115 + t5;
t290 = t101 * t105;
t169 = -pkin(2) * t232 + pkin(9) * t290;
t68 = -t291 * pkin(3) + t169;
t254 = t104 * t290;
t73 = -t108 * t291 + t254;
t288 = t101 * t108;
t74 = t104 * t291 + t105 * t288;
t137 = t73 * pkin(4) - t74 * pkin(11) + t68;
t233 = t105 * t291;
t220 = pkin(2) * t233;
t166 = t291 * pkin(10) + t220;
t155 = qJD(4) * t166;
t69 = t169 * qJD(3);
t138 = -t104 * t155 - t108 * t69;
t210 = pkin(3) * t105 - pkin(10) * t109;
t249 = -pkin(10) * t105 - pkin(2);
t200 = pkin(3) * t109 - t249;
t311 = pkin(9) * t109;
t151 = t104 * t311 + t108 * t200;
t320 = t151 * qJD(4);
t322 = -(-t320 + (t105 * pkin(11) + t104 * t210) * qJD(3)) * t101 - t138 - qJD(5) * t137;
t275 = qJD(3) * t109;
t240 = t101 * t275;
t281 = t74 * qJD(4);
t59 = t104 * t240 + t281;
t191 = -t108 * t59 + t73 * t273;
t103 = sin(qJ(5));
t97 = t103 ^ 2;
t107 = cos(qJ(5));
t99 = t107 ^ 2;
t307 = t97 - t99;
t235 = qJD(5) * t307;
t181 = t104 * t200;
t321 = (t108 * t311 - t181) * qJD(4);
t228 = t291 * qJD(4);
t172 = t228 + t240;
t89 = qJD(4) * t254;
t316 = -t108 * t172 + t89;
t95 = t101 ^ 2;
t315 = t109 ^ 2;
t314 = 0.2e1 * t101;
t313 = pkin(2) * t101;
t312 = pkin(9) * t104;
t310 = pkin(10) * t101;
t309 = pkin(10) * t104;
t308 = pkin(10) * t108;
t17 = t104 * t126 + t108 * t117;
t259 = pkin(9) * t287;
t53 = -t101 * t181 + t108 * (t166 + t259);
t23 = t104 * t218 - t58 * t273 + (-qJD(4) * t182 + t44) * t108;
t269 = qJD(5) * t103;
t13 = -t47 * t269 + (qJD(5) * t57 + t23) * t107 + t124 * t103;
t306 = t103 * t13;
t33 = t103 * t47 - t57 * t107;
t305 = t103 * t33;
t276 = qJD(3) * t105;
t241 = t101 * t276;
t60 = t103 * t74 + t107 * t287;
t38 = qJD(5) * t60 - t103 * t241 + t107 * t316;
t304 = t103 * t38;
t303 = t103 * t60;
t302 = t104 * t23;
t34 = t103 * t57 + t107 * t47;
t12 = qJD(5) * t34 + t103 * t23 - t124 * t107;
t301 = t107 * t12;
t300 = t107 * t13;
t299 = t107 * t34;
t298 = t107 * t38;
t252 = t103 * t287;
t268 = qJD(5) * t107;
t39 = -qJD(5) * t252 - t103 * t316 - t107 * t241 + t74 * t268;
t297 = t107 * t39;
t61 = t107 * t74 - t252;
t296 = t107 * t61;
t22 = qJD(4) * t47 + t44 * t104 - t108 * t218;
t295 = t108 * t22;
t98 = t104 ^ 2;
t293 = t108 ^ 2 - t98;
t282 = t107 * t108;
t274 = qJD(4) * t103;
t272 = qJD(4) * t107;
t271 = qJD(4) * t108;
t270 = qJD(4) * t109;
t267 = qJD(5) * t108;
t266 = 0.2e1 * t46 * t22;
t265 = 0.2e1 * t73 * t59;
t264 = -0.2e1 * pkin(3) * qJD(4);
t263 = -0.2e1 * pkin(4) * qJD(5);
t262 = t103 * t308;
t261 = pkin(8) * t286;
t260 = pkin(10) * t282;
t258 = pkin(10) * t271;
t96 = t102 ^ 2;
t257 = t96 * t278;
t256 = t95 * t275;
t253 = t102 * t283;
t246 = t104 * t270;
t245 = t107 * t271;
t244 = t103 * t267;
t243 = t104 * t268;
t242 = t107 * t267;
t237 = t103 * t268;
t236 = t104 * t271;
t230 = qJD(3) * t291;
t229 = t293 * qJD(4);
t227 = 0.2e1 * t236;
t226 = t95 * t239;
t225 = t98 * t237;
t224 = t105 * t256;
t223 = t106 * t257;
t219 = t107 * t236;
t212 = t292 * t280;
t211 = t105 * t230;
t209 = -pkin(4) * t108 - pkin(11) * t104;
t208 = pkin(4) * t104 - pkin(11) * t108;
t15 = pkin(11) * t57 + t17;
t7 = -t103 * t15 + t107 * t115;
t8 = t103 * t115 + t107 * t15;
t207 = -t103 * t8 - t107 * t7;
t206 = t22 * t73 + t46 * t59;
t49 = -pkin(11) * t287 + t53;
t27 = -t103 * t49 + t107 * t137;
t28 = t103 * t137 + t107 * t49;
t205 = -t103 * t28 - t107 * t27;
t204 = -t103 * t34 - t107 * t33;
t203 = -t103 * t61 - t107 * t60;
t199 = pkin(3) - t209;
t180 = t107 * t199;
t65 = -t180 - t262;
t66 = -t103 * t199 + t260;
t202 = -t103 * t66 - t107 * t65;
t183 = qJD(3) * t210;
t35 = (-t104 * t183 + t320) * t101 - t138;
t139 = t104 * t69 - t108 * t155;
t36 = (t108 * t183 - t321) * t101 + t139;
t201 = -t36 * t104 - t35 * t108;
t196 = -t102 * pkin(1) - pkin(2) * t285;
t145 = t167 - t261;
t16 = -t104 * (-t109 * t171 + t145 * t233 + t196 * t290 - t176) + t108 * (-t101 * t145 + t291 * t196 + t217) - (t214 + (t109 * t231 + (-t291 ^ 2 - t95) * t284) * t102) * t312;
t14 = -t57 * pkin(4) - t16;
t6 = t324 * t104 - t323 * t108;
t4 = -pkin(4) * t124 - t6;
t195 = t103 * t4 + t14 * t268;
t194 = -t107 * t4 + t14 * t269;
t192 = t46 * t273 - t295;
t190 = t103 * t22 + t46 * t268;
t189 = -t107 * t22 + t46 * t269;
t32 = (t321 + (-t105 * pkin(4) - t108 * t210) * qJD(3)) * t101 - t139;
t157 = t104 * t166;
t48 = t157 + (-t108 * t249 + (pkin(3) * t108 + pkin(4) + t312) * t109) * t101;
t188 = t103 * t32 + t48 * t268;
t187 = -t107 * t32 + t48 * t269;
t186 = t103 * t59 + t73 * t268;
t185 = -t107 * t59 + t73 * t269;
t184 = t208 * t103;
t178 = t104 * t276 - t108 * t270;
t177 = t104 * t272 + t244;
t173 = t101 * t182;
t78 = t220 + t259;
t170 = -t221 + t261;
t164 = qJD(3) * t173;
t21 = t317 * t105 - t318 * t109;
t19 = -pkin(3) * t218 - t21;
t112 = t22 * pkin(4) - t23 * pkin(11) + t19;
t1 = -t103 * t112 + t325 * t107 + t15 * t269;
t2 = t325 * t103 + t107 * t112 - t15 * t268;
t147 = t207 * qJD(5) - t1 * t107 - t103 * t2;
t146 = -t104 * t6 - t108 * t5 + (-t104 * t17 - t108 * t16) * qJD(4);
t127 = pkin(2) * t211 + t59 * pkin(4) + pkin(9) * t240 + pkin(11) * t316;
t10 = -t103 * t127 + t322 * t107 + t49 * t269;
t11 = t322 * t103 + t107 * t127 - t49 * t268;
t144 = t205 * qJD(5) - t10 * t107 - t103 * t11;
t50 = pkin(10) * t177 - qJD(4) * t184 + qJD(5) * t180;
t51 = -t66 * qJD(5) + (t103 * t309 + t107 * t208) * qJD(4);
t143 = t202 * qJD(5) - t103 * t51 - t107 * t50;
t120 = t124 * t104 + t57 * t271;
t93 = -0.2e1 * t236;
t84 = -0.2e1 * t224;
t72 = t171 * qJD(2);
t71 = t170 * qJD(2);
t70 = t78 * qJD(3);
t62 = -t103 * t245 + t104 * t235;
t52 = -t101 * t151 - t157;
t45 = (pkin(2) * t216 - t101 * t171) * qJD(2);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t223, 0.2e1 * (-t106 ^ 2 + t110 ^ 2) * t96 * qJD(2), 0.2e1 * t110 * t212, -0.2e1 * t223, -0.2e1 * t106 * t212, 0, -0.2e1 * pkin(1) * t279 * t96 + 0.2e1 * t292 * t72, -0.2e1 * pkin(1) * t257 + 0.2e1 * t292 * t71, 0.2e1 * (-t106 * t72 - t110 * t71 + (t106 * t171 + t110 * t170) * qJD(2)) * t102, -0.2e1 * t170 * t72 + 0.2e1 * t171 * t71, 0.2e1 * t58 * t44, -0.2e1 * t124 * t58 - 0.2e1 * t44 * t57, -0.2e1 * t182 * t44 + 0.2e1 * t218 * t58, t335, 0.2e1 * t124 * t182 - 0.2e1 * t218 * t57, -0.2e1 * t173 * t239, -0.2e1 * t21 * t182 + 0.2e1 * t45 * t57 + 0.2e1 * t41 * t197 + 0.2e1 * (t41 * t161 + (t174 * t41 + t289 * t30) * qJD(2)) * t102, -0.2e1 * t182 * t20 - 0.2e1 * t218 * t31 + 0.2e1 * t41 * t44 + 0.2e1 * t45 * t58, -0.2e1 * t124 * t31 + 0.2e1 * t20 * t57 - 0.2e1 * t21 * t58 - 0.2e1 * t30 * t44, -0.2e1 * t20 * t31 + 0.2e1 * t21 * t30 + 0.2e1 * t41 * t45, 0.2e1 * t47 * t23, -0.2e1 * t22 * t47 - 0.2e1 * t23 * t46, 0.2e1 * t124 * t47 + 0.2e1 * t23 * t57, t266, -0.2e1 * t124 * t46 - 0.2e1 * t22 * t57, t335, 0.2e1 * t124 * t16 + 0.2e1 * t19 * t46 + 0.2e1 * t26 * t22 + 0.2e1 * t6 * t57, -0.2e1 * t124 * t17 + 0.2e1 * t19 * t47 + 0.2e1 * t26 * t23 + 0.2e1 * t5 * t57, -0.2e1 * t16 * t23 - 0.2e1 * t17 * t22 + 0.2e1 * t46 * t5 - 0.2e1 * t47 * t6, 0.2e1 * t16 * t6 - 0.2e1 * t17 * t5 + 0.2e1 * t19 * t26, 0.2e1 * t34 * t13, -0.2e1 * t12 * t34 - 0.2e1 * t13 * t33, 0.2e1 * t13 * t46 + 0.2e1 * t22 * t34, 0.2e1 * t33 * t12, -0.2e1 * t12 * t46 - 0.2e1 * t22 * t33, t266, 0.2e1 * t12 * t14 + 0.2e1 * t2 * t46 + 0.2e1 * t22 * t7 + 0.2e1 * t33 * t4, 0.2e1 * t1 * t46 + 0.2e1 * t13 * t14 - 0.2e1 * t22 * t8 + 0.2e1 * t34 * t4, 0.2e1 * t1 * t33 - 0.2e1 * t12 * t8 - 0.2e1 * t13 * t7 - 0.2e1 * t2 * t34, -0.2e1 * t1 * t8 + 0.2e1 * t14 * t4 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, 0, -t239, 0, t72, t71, 0, 0, (t105 * t44 + t275 * t58) * t101, ((-qJD(3) * t57 + t44) * t109 + (-t330 + (-t105 * t163 - t253 - t58) * qJD(3)) * t105) * t101, t105 * t226 - t109 * t164 + t291 * t44, t331, t105 * t164 + t109 * t226 - t124 * t291, t101 * t198, -t124 * t313 - t169 * t218 + t182 * t70 + t21 * t291 + t241 * t41 - t287 * t45, -t69 * t182 + t20 * t291 + (-pkin(2) * t44 + t45 * t105 - t239 * t78 + t275 * t41) * t101, -t124 * t78 + t169 * t44 - t20 * t287 - t21 * t290 - t240 * t30 - t241 * t31 + t69 * t57 + t70 * t58, -t169 * t21 - t20 * t78 - t30 * t70 - t31 * t69 - t45 * t313, t23 * t74 - t316 * t47, -t74 * t22 - t23 * t73 + t316 * t46 - t47 * t59, t124 * t74 - t23 * t287 + t241 * t47 - t316 * t57, t206, t22 * t287 - t59 * t57 - t73 * t330 + (-t73 * t253 + (-t46 * t101 - t163 * t73) * t105) * qJD(3), t331, t124 * t52 + t16 * t241 + t19 * t73 + t68 * t22 + t26 * t59 - t287 * t6 + t36 * t57 + t70 * t46, -t124 * t53 - t17 * t241 + t19 * t74 + t68 * t23 - t26 * t316 - t287 * t5 + t35 * t57 + t70 * t47, t16 * t316 - t17 * t59 - t53 * t22 - t52 * t23 + t35 * t46 - t36 * t47 + t5 * t73 - t6 * t74, t16 * t36 - t17 * t35 + t19 * t68 + t26 * t70 - t5 * t53 + t52 * t6, t13 * t61 - t34 * t38, -t12 * t61 - t13 * t60 + t33 * t38 - t34 * t39, t13 * t73 + t22 * t61 + t34 * t59 - t38 * t46, t12 * t60 + t33 * t39, -t12 * t73 - t22 * t60 - t33 * t59 - t39 * t46, t206, t11 * t46 + t12 * t48 + t14 * t39 + t2 * t73 + t22 * t27 + t32 * t33 + t4 * t60 + t59 * t7, t1 * t73 + t10 * t46 + t13 * t48 - t14 * t38 - t22 * t28 + t32 * t34 + t4 * t61 - t59 * t8, t1 * t60 + t10 * t33 - t11 * t34 - t12 * t28 - t13 * t27 - t2 * t61 + t38 * t7 - t39 * t8, -t1 * t28 - t10 * t8 + t11 * t7 + t14 * t32 + t2 * t27 + t4 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t224, 0.2e1 * (-t105 ^ 2 + t315) * t95 * qJD(3), 0.2e1 * t230 * t287, t84, -0.2e1 * t101 * t211, 0, -0.2e1 * pkin(2) * t276 * t95 - 0.2e1 * t291 * t70, -0.2e1 * pkin(2) * t256 + 0.2e1 * t291 * t69, (t105 * t70 - t109 * t69 + (-t105 * t78 + t109 * t169) * qJD(3)) * t314, 0.2e1 * t169 * t70 - 0.2e1 * t69 * t78, -0.2e1 * t74 * t316, 0.2e1 * t316 * t73 - 0.2e1 * t74 * t59, (-(t108 * t228 - t89) * t109 + (t105 * t74 - t315 * t288) * qJD(3)) * t314, t265, (t109 * t59 - t276 * t73) * t314, t84, 0.2e1 * t59 * t68 + 0.2e1 * t70 * t73 + 0.2e1 * (-t109 * t36 + t276 * t52) * t101, -0.2e1 * t241 * t53 - 0.2e1 * t287 * t35 - 0.2e1 * t316 * t68 + 0.2e1 * t70 * t74, 0.2e1 * t316 * t52 + 0.2e1 * t35 * t73 - 0.2e1 * t36 * t74 - 0.2e1 * t53 * t59, -0.2e1 * t35 * t53 + 0.2e1 * t36 * t52 + 0.2e1 * t68 * t70, -0.2e1 * t61 * t38, 0.2e1 * t38 * t60 - 0.2e1 * t39 * t61, -0.2e1 * t38 * t73 + 0.2e1 * t59 * t61, 0.2e1 * t60 * t39, -0.2e1 * t39 * t73 - 0.2e1 * t59 * t60, t265, 0.2e1 * t11 * t73 + 0.2e1 * t27 * t59 + 0.2e1 * t32 * t60 + 0.2e1 * t39 * t48, 0.2e1 * t10 * t73 - 0.2e1 * t28 * t59 + 0.2e1 * t32 * t61 - 0.2e1 * t38 * t48, 0.2e1 * t10 * t60 - 0.2e1 * t11 * t61 + 0.2e1 * t27 * t38 - 0.2e1 * t28 * t39, -0.2e1 * t10 * t28 + 0.2e1 * t11 * t27 + 0.2e1 * t32 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t124, t218, t21, t20, 0, 0, t271 * t47 + t302, -t104 * t22 + t108 * t23 + (-t104 * t47 - t108 * t46) * qJD(4), t120, t192, -t333, 0, -pkin(3) * t22 - pkin(10) * t120 - t19 * t108 + t26 * t273, -pkin(3) * t23 + t333 * pkin(10) + t19 * t104 + t26 * t271, (t302 - t295 + (t104 * t46 + t108 * t47) * qJD(4)) * pkin(10) + t146, -pkin(3) * t19 + pkin(10) * t146, t34 * t245 + (-t269 * t34 + t300) * t104, t204 * t271 + (-t306 - t301 + (-t299 + t305) * qJD(5)) * t104, (t272 * t46 - t13) * t108 + (qJD(4) * t34 - t189) * t104, t33 * t243 + (t104 * t12 + t271 * t33) * t103, (-t274 * t46 + t12) * t108 + (-qJD(4) * t33 - t190) * t104, t192, t22 * t65 + t46 * t51 + (-t2 + (pkin(10) * t33 + t103 * t14) * qJD(4)) * t108 + (pkin(10) * t12 + qJD(4) * t7 + t195) * t104, -t22 * t66 + t46 * t50 + (-t1 + (pkin(10) * t34 + t107 * t14) * qJD(4)) * t108 + (pkin(10) * t13 - qJD(4) * t8 - t194) * t104, -t12 * t66 - t13 * t65 + t33 * t50 - t34 * t51 + t207 * t271 + (t1 * t103 - t107 * t2 + (t103 * t7 - t107 * t8) * qJD(5)) * t104, -t1 * t66 + t2 * t65 - t50 * t8 + t51 * t7 + (t104 * t4 + t14 * t271) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, 0, -t241, 0, -t70, t69, 0, 0, -t89 * t104 + (t104 * t172 + t281) * t108, -t104 * t59 - t108 * t316 - t271 * t73 - t273 * t74, t178 * t101, t191, (t108 * t276 + t246) * t101, 0, -pkin(3) * t59 - t108 * t70 - t178 * t310 + t273 * t68, pkin(3) * t316 + t70 * t104 - t241 * t308 - t246 * t310 + t271 * t68, t191 * pkin(10) + t258 * t74 - t271 * t52 - t273 * t53 - t309 * t316 + t201, -pkin(3) * t70 + ((-t53 * t104 - t52 * t108) * qJD(4) + t201) * pkin(10), t61 * t245 + (-t269 * t61 - t298) * t104, t203 * t271 + (t304 - t297 + (-t296 + t303) * qJD(5)) * t104, (t272 * t73 + t38) * t108 + (qJD(4) * t61 - t185) * t104, t60 * t243 + (t104 * t39 + t271 * t60) * t103, (-t274 * t73 + t39) * t108 + (-qJD(4) * t60 - t186) * t104, t191, t51 * t73 + t59 * t65 + (-t11 + (pkin(10) * t60 + t103 * t48) * qJD(4)) * t108 + (pkin(10) * t39 + qJD(4) * t27 + t188) * t104, t50 * t73 - t59 * t66 + (-t10 + (pkin(10) * t61 + t107 * t48) * qJD(4)) * t108 + (-pkin(10) * t38 - qJD(4) * t28 - t187) * t104, t38 * t65 - t39 * t66 + t50 * t60 - t51 * t61 + t205 * t271 + (t10 * t103 - t107 * t11 + (t103 * t27 - t107 * t28) * qJD(5)) * t104, -t10 * t66 + t11 * t65 + t27 * t51 - t28 * t50 + (t104 * t32 + t271 * t48) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, 0.2e1 * t229, 0, t93, 0, 0, t104 * t264, t108 * t264, 0, 0, 0.2e1 * t236 * t99 - 0.2e1 * t225, -0.4e1 * t103 * t219 + 0.2e1 * t235 * t98, 0.2e1 * t104 * t244 - 0.2e1 * t272 * t293, 0.2e1 * t236 * t97 + 0.2e1 * t225, 0.2e1 * t103 * t229 + 0.2e1 * t104 * t242, t93, 0.2e1 * t65 * t273 - 0.2e1 * t108 * t51 + 0.2e1 * (t103 * t227 + t268 * t98) * pkin(10), -0.2e1 * t66 * t273 - 0.2e1 * t108 * t50 + 0.2e1 * (-t269 * t98 + 0.2e1 * t219) * pkin(10), 0.2e1 * t202 * t271 + 0.2e1 * (t103 * t50 - t107 * t51 + (t103 * t65 - t107 * t66) * qJD(5)) * t104, 0.2e1 * pkin(10) ^ 2 * t236 - 0.2e1 * t50 * t66 + 0.2e1 * t51 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, -t22, t124, t6, t5, 0, 0, t268 * t34 + t306, qJD(5) * t204 - t103 * t12 + t300, t190, t269 * t33 - t301, -t189, 0, -pkin(4) * t12 - pkin(11) * t190 + t194, -pkin(4) * t13 + pkin(11) * t189 + t195, (t306 - t301 + (t299 + t305) * qJD(5)) * pkin(11) + t147, -pkin(4) * t4 + pkin(11) * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t316, 0, -t59, t241, t36, t35, 0, 0, t268 * t61 - t304, qJD(5) * t203 - t103 * t39 - t298, t186, t269 * t60 - t297, -t185, 0, -pkin(4) * t39 - pkin(11) * t186 + t187, pkin(4) * t38 + pkin(11) * t185 + t188, (-t304 - t297 + (t296 + t303) * qJD(5)) * pkin(11) + t144, -pkin(4) * t32 + pkin(11) * t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, 0, -t273, 0, -t258, pkin(10) * t273, 0, 0, -t62, -0.4e1 * t104 * t237 - t271 * t307, t103 * t273 - t242, t62, t177, 0, (pkin(11) * t282 + (-pkin(4) * t107 + pkin(10) * t103) * t104) * qJD(5) + (t103 * t209 - t260) * qJD(4), (t107 * t309 + t184) * qJD(5) + (t107 * t209 + t262) * qJD(4), t143, -pkin(4) * t258 + pkin(11) * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t237, -0.2e1 * t235, 0, -0.2e1 * t237, 0, 0, t103 * t263, t107 * t263, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, -t12, t22, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, -t39, t59, t11, t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104 * t269 + t245, 0, -t103 * t271 - t243, t273, t51, t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, 0, -t269, 0, -pkin(11) * t268, pkin(11) * t269, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
