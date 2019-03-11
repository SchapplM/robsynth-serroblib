% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPR13_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:04:04
% EndTime: 2019-03-10 00:04:21
% DurationCPUTime: 5.63s
% Computational Cost: add. (5754->462), mult. (15545->805), div. (0->0), fcn. (14675->10), ass. (0->218)
t176 = cos(qJ(4));
t163 = qJD(4) * t176;
t172 = sin(qJ(6));
t173 = sin(qJ(4));
t288 = cos(qJ(6));
t229 = qJD(4) * t288;
t303 = -t163 * t172 + t173 * t229;
t174 = sin(qJ(3));
t164 = qJD(3) * t174;
t177 = cos(qJ(3));
t257 = qJD(4) * t177;
t240 = t176 * t257;
t302 = t164 * t173 - t240;
t301 = -0.4e1 * t174;
t171 = sin(pkin(6));
t175 = sin(qJ(2));
t273 = t171 * t175;
t275 = cos(pkin(6));
t112 = t174 * t273 - t177 * t275;
t259 = qJD(4) * t173;
t113 = t174 * t275 + t177 * t273;
t178 = cos(qJ(2));
t264 = qJD(2) * t178;
t238 = t171 * t264;
t74 = qJD(3) * t113 + t174 * t238;
t194 = t112 * t259 - t176 * t74;
t300 = pkin(10) * t194;
t272 = t171 * t178;
t75 = -qJD(3) * t112 + t177 * t238;
t299 = -qJD(4) * t272 + t75;
t298 = -t172 * t259 - t176 * t229;
t169 = t176 ^ 2;
t267 = t173 ^ 2 - t169;
t225 = t267 * qJD(4);
t228 = qJD(6) * t288;
t256 = qJD(6) * t172;
t297 = -t173 * t228 + t176 * t256 + t303;
t125 = t172 * t173 + t176 * t288;
t72 = qJD(6) * t125 + t298;
t241 = t173 * t257;
t262 = qJD(3) * t176;
t189 = t174 * t262 + t241;
t285 = pkin(10) * t177;
t213 = pkin(3) * t174 - t285;
t127 = t213 * qJD(3);
t283 = t174 * pkin(10);
t214 = -pkin(3) * t177 - t283;
t134 = -pkin(2) + t214;
t268 = t127 * t173 + t134 * t163;
t65 = pkin(9) * t189 - t268;
t246 = pkin(1) * t275;
t98 = -pkin(2) * t275 + pkin(8) * t273 - t178 * t246;
t60 = pkin(3) * t112 - pkin(10) * t113 + t98;
t100 = (-pkin(2) * t178 - pkin(9) * t175 - pkin(1)) * t171;
t186 = pkin(8) * t272 + t175 * t246;
t99 = pkin(9) * t275 + t186;
t281 = t100 * t174 + t177 * t99;
t62 = -pkin(10) * t272 + t281;
t25 = -t173 * t62 + t176 * t60;
t291 = pkin(4) + pkin(5);
t77 = t113 * t176 - t173 * t272;
t20 = -pkin(11) * t77 - t112 * t291 - t25;
t26 = t173 * t60 + t176 * t62;
t23 = qJ(5) * t112 + t26;
t76 = t113 * t173 + t176 * t272;
t21 = pkin(11) * t76 + t23;
t10 = t172 * t20 + t21 * t288;
t265 = qJD(2) * t175;
t239 = t171 * t265;
t105 = (pkin(2) * t175 - pkin(9) * t178) * t171 * qJD(2);
t226 = qJD(2) * t275;
t215 = t178 * t226;
t106 = -pkin(1) * t215 + pkin(8) * t239;
t261 = qJD(3) * t177;
t32 = -t100 * t261 - t105 * t174 + t106 * t177 + t164 * t99;
t30 = pkin(10) * t239 - t32;
t107 = t186 * qJD(2);
t36 = t74 * pkin(3) - t75 * pkin(10) + t107;
t227 = t163 * t62 + t173 * t30 - t176 * t36 + t259 * t60;
t44 = -t113 * t259 + t173 * t239 + t176 * t299;
t4 = -t44 * pkin(11) - t291 * t74 + t227;
t43 = t113 * t163 + t173 * t299 - t176 * t239;
t104 = t112 * qJD(5);
t12 = -t163 * t60 - t173 * t36 - t176 * t30 + t259 * t62;
t71 = t74 * qJ(5);
t8 = t104 - t12 + t71;
t6 = pkin(11) * t43 + t8;
t2 = -qJD(6) * t10 - t172 * t6 + t288 * t4;
t286 = pkin(9) * t173;
t158 = t177 * t286;
t165 = t177 * pkin(4);
t284 = pkin(11) * t174;
t68 = t177 * pkin(5) + t158 + t165 + (-t134 - t284) * t176;
t269 = t176 * t177;
t159 = pkin(9) * t269;
t96 = t134 * t173 + t159;
t81 = -qJ(5) * t177 + t96;
t69 = t173 * t284 + t81;
t41 = t172 * t68 + t288 * t69;
t222 = -pkin(9) * t302 - t176 * t127 + t134 * t259;
t258 = qJD(4) * t174;
t242 = t173 * t258;
t42 = pkin(11) * t242 + (-pkin(11) * t269 - t174 * t291) * qJD(3) + t222;
t161 = qJ(5) * t164;
t270 = t174 * t176;
t45 = t161 + (-pkin(9) * qJD(3) + pkin(11) * qJD(4)) * t270 + (-qJD(5) + (-pkin(9) * qJD(4) + pkin(11) * qJD(3)) * t173) * t177 + t268;
t16 = -qJD(6) * t41 - t172 * t45 + t288 * t42;
t271 = t173 * qJ(5);
t296 = -t176 * t291 - t271;
t274 = qJ(5) * t176;
t211 = pkin(4) * t173 - t274;
t204 = pkin(9) + t211;
t101 = t204 * t174;
t255 = t173 * qJD(5);
t111 = qJD(4) * t211 - t255;
t212 = pkin(4) * t176 + t271;
t130 = -pkin(3) - t212;
t295 = qJD(3) * (-t130 * t177 + t283) - qJD(4) * t101 - t111 * t174;
t254 = t176 * qJD(5);
t294 = qJD(4) * t212 - t254;
t293 = 0.2e1 * t171;
t292 = 0.2e1 * qJD(5);
t290 = pkin(10) - pkin(11);
t289 = t74 * pkin(4);
t287 = pkin(9) * t171;
t282 = t100 * t177 - t174 * t99;
t280 = t173 * t76;
t279 = t176 * t43;
t278 = t176 * t77;
t277 = t44 * t173;
t276 = t44 * t176;
t168 = t174 ^ 2;
t266 = -t177 ^ 2 + t168;
t263 = qJD(3) * t173;
t260 = qJD(3) * t178;
t253 = t177 * qJD(5);
t252 = -0.2e1 * pkin(2) * qJD(3);
t251 = -0.2e1 * pkin(3) * qJD(4);
t250 = 0.2e1 * t112 * t74;
t61 = pkin(3) * t272 - t282;
t249 = pkin(4) * t164;
t248 = pkin(10) * t259;
t247 = pkin(10) * t163;
t245 = t288 * t173;
t166 = t171 ^ 2;
t244 = t166 * t264;
t234 = t173 * t261;
t233 = t173 * t163;
t232 = t174 * t261;
t231 = t176 * t261;
t230 = qJD(3) * t288;
t33 = -t100 * t164 + t105 * t177 + t106 * t174 - t261 * t99;
t95 = t134 * t176 - t158;
t224 = t266 * qJD(3);
t223 = 0.2e1 * t232;
t221 = t175 * t244;
t220 = t173 * t231;
t219 = t177 * t230;
t24 = -pkin(4) * t112 - t25;
t209 = -t173 * t23 + t176 * t24;
t208 = -t173 * t77 - t176 * t76;
t82 = t165 - t95;
t207 = -t173 * t81 + t176 * t82;
t206 = t77 * qJ(5) - t61;
t1 = -t172 * t4 - t20 * t228 + t21 * t256 - t288 * t6;
t46 = t172 * t77 - t288 * t76;
t47 = t172 * t76 + t288 * t77;
t188 = pkin(3) * t239 + t33;
t180 = t44 * qJ(5) + t77 * qJD(5) + t188;
t14 = pkin(4) * t43 - t180;
t27 = pkin(4) * t76 - t206;
t201 = -t14 * t173 - t163 * t27;
t200 = -t14 * t176 + t259 * t27;
t199 = t163 * t61 - t173 * t188;
t198 = t176 * t188 + t259 * t61;
t197 = -t173 * t291 + t274;
t196 = t112 * t164 - t177 * t74;
t195 = t112 * t163 + t173 * t74;
t15 = -t172 * t42 - t228 * t68 + t256 * t69 - t288 * t45;
t129 = qJ(5) * t288 - t172 * t291;
t142 = t290 * t173;
t143 = t290 * t176;
t80 = t142 * t172 + t143 * t288;
t193 = -pkin(9) + t197;
t192 = t195 * pkin(10);
t191 = t174 * t260 + t177 * t265;
t190 = t174 * t265 - t177 * t260;
t11 = t227 - t289;
t182 = qJD(4) * t209 + t11 * t173 + t8 * t176;
t57 = t161 - t65 - t253;
t63 = t222 - t249;
t181 = qJD(4) * t207 + t63 * t173 + t57 * t176;
t150 = -0.2e1 * t232;
t149 = pkin(10) * t240;
t128 = -qJ(5) * t172 - t288 * t291;
t126 = -t172 * t176 + t245;
t122 = pkin(3) - t296;
t114 = t231 - t242;
t109 = t125 * t174;
t108 = t172 * t270 - t174 * t245;
t103 = qJD(5) * t172 + qJD(6) * t129;
t102 = qJ(5) * t256 - qJD(5) * t288 + t228 * t291;
t90 = qJD(4) * t197 + t255;
t79 = t142 * t288 - t143 * t172;
t78 = t193 * t174;
t64 = t174 * t294 + t204 * t261;
t54 = t80 * qJD(6) + t290 * t298;
t53 = -t142 * t228 + t143 * t256 + t290 * t303;
t50 = t172 * t231 - t173 * t219 + t174 * t72;
t49 = -t172 * t234 + t174 * t297 - t176 * t219;
t48 = (qJD(4) * t296 + t254) * t174 + t193 * t261;
t40 = -t172 * t69 + t288 * t68;
t22 = -t291 * t76 + t206;
t18 = -qJD(6) * t46 + t43 * t172 + t288 * t44;
t17 = qJD(6) * t47 + t172 * t44 - t288 * t43;
t9 = -t172 * t21 + t20 * t288;
t7 = -t291 * t43 + t180;
t3 = [0, 0, 0, 0.2e1 * t221, 0.2e1 * (-t175 ^ 2 + t178 ^ 2) * t166 * qJD(2), t215 * t293, -0.2e1 * t226 * t273, 0, -0.2e1 * pkin(1) * t166 * t265 - 0.2e1 * t107 * t275, -0.2e1 * pkin(1) * t244 + 0.2e1 * t106 * t275, 0.2e1 * t113 * t75, -0.2e1 * t112 * t75 - 0.2e1 * t113 * t74 (t113 * t265 - t178 * t75) * t293 (-t112 * t265 + t178 * t74) * t293, -0.2e1 * t221, 0.2e1 * t107 * t112 + 0.2e1 * t98 * t74 + 0.2e1 * (-t33 * t178 + t265 * t282) * t171, 0.2e1 * t107 * t113 + 0.2e1 * t98 * t75 + 0.2e1 * (-t32 * t178 - t265 * t281) * t171, 0.2e1 * t77 * t44, -0.2e1 * t43 * t77 - 0.2e1 * t44 * t76, 0.2e1 * t112 * t44 + 0.2e1 * t74 * t77, -0.2e1 * t112 * t43 - 0.2e1 * t74 * t76, t250, -0.2e1 * t112 * t227 - 0.2e1 * t188 * t76 + 0.2e1 * t25 * t74 + 0.2e1 * t43 * t61, 0.2e1 * t112 * t12 - 0.2e1 * t188 * t77 - 0.2e1 * t26 * t74 + 0.2e1 * t44 * t61, -0.2e1 * t11 * t112 + 0.2e1 * t14 * t76 - 0.2e1 * t24 * t74 + 0.2e1 * t27 * t43, 0.2e1 * t11 * t77 - 0.2e1 * t23 * t43 + 0.2e1 * t24 * t44 - 0.2e1 * t76 * t8, 0.2e1 * t112 * t8 - 0.2e1 * t14 * t77 + 0.2e1 * t23 * t74 - 0.2e1 * t27 * t44, 0.2e1 * t11 * t24 + 0.2e1 * t14 * t27 + 0.2e1 * t23 * t8, 0.2e1 * t47 * t18, -0.2e1 * t17 * t47 - 0.2e1 * t18 * t46, -0.2e1 * t112 * t18 - 0.2e1 * t47 * t74, 0.2e1 * t112 * t17 + 0.2e1 * t46 * t74, t250, -0.2e1 * t112 * t2 + 0.2e1 * t17 * t22 + 0.2e1 * t46 * t7 - 0.2e1 * t74 * t9, -0.2e1 * t1 * t112 + 0.2e1 * t10 * t74 + 0.2e1 * t18 * t22 + 0.2e1 * t47 * t7; 0, 0, 0, 0, 0, t238, -t239, 0, -t107, t106, t113 * t261 + t174 * t75, -t174 * t74 + t75 * t177 + (-t112 * t177 - t113 * t174) * qJD(3), t190 * t171, t191 * t171, 0, -pkin(2) * t74 - t107 * t177 + t164 * t98 - t190 * t287, -pkin(2) * t75 + t107 * t174 - t191 * t287 + t261 * t98, t77 * t231 + (-t259 * t77 + t276) * t174, t208 * t261 + (-t277 - t279 + (-t278 + t280) * qJD(4)) * t174 (t112 * t262 - t44) * t177 + (qJD(3) * t77 - t194) * t174 (-t112 * t263 + t43) * t177 + (-qJD(3) * t76 - t195) * t174, t196, -t222 * t112 + t95 * t74 + (t227 + (pkin(9) * t76 + t173 * t61) * qJD(3)) * t177 + (pkin(9) * t43 + qJD(3) * t25 + t199) * t174, t65 * t112 - t96 * t74 + (-t12 + (pkin(9) * t77 + t176 * t61) * qJD(3)) * t177 + (pkin(9) * t44 - qJD(3) * t26 - t198) * t174, t101 * t43 - t63 * t112 + t64 * t76 - t82 * t74 + (t263 * t27 + t11) * t177 + (-qJD(3) * t24 - t201) * t174, -t81 * t43 + t82 * t44 - t57 * t76 + t63 * t77 + t209 * t261 + (t11 * t176 - t173 * t8 + (-t173 * t24 - t176 * t23) * qJD(4)) * t174, -t101 * t44 + t57 * t112 - t64 * t77 + t81 * t74 + (-t262 * t27 - t8) * t177 + (qJD(3) * t23 + t200) * t174, t101 * t14 + t11 * t82 + t23 * t57 + t24 * t63 + t27 * t64 + t8 * t81, t109 * t18 - t47 * t49, -t108 * t18 - t109 * t17 + t46 * t49 - t47 * t50, -t109 * t74 + t112 * t49 - t164 * t47 + t177 * t18, t108 * t74 + t112 * t50 + t164 * t46 - t17 * t177, t196, t108 * t7 - t112 * t16 - t164 * t9 + t17 * t78 + t177 * t2 + t22 * t50 - t40 * t74 + t46 * t48, t1 * t177 + t10 * t164 + t109 * t7 - t112 * t15 + t18 * t78 - t22 * t49 + t41 * t74 + t47 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, -0.2e1 * t224, 0, 0, 0, t174 * t252, t177 * t252, -0.2e1 * t168 * t233 + 0.2e1 * t169 * t232, 0.2e1 * t168 * t225 + t220 * t301, 0.2e1 * t174 * t241 + 0.2e1 * t262 * t266, -0.2e1 * t173 * t224 + 0.2e1 * t174 * t240, t150, 0.2e1 * t95 * t164 + 0.2e1 * t222 * t177 + 0.2e1 * (t163 * t168 + t173 * t223) * pkin(9), -0.2e1 * t96 * t164 - 0.2e1 * t65 * t177 + 0.2e1 * (-t168 * t259 + t176 * t223) * pkin(9), 0.2e1 * (t101 * t263 + t63) * t177 + 0.2e1 * (-qJD(3) * t82 + t101 * t163 + t64 * t173) * t174, 0.2e1 * t207 * t261 + 0.2e1 * (-t173 * t57 + t176 * t63 + (-t173 * t82 - t176 * t81) * qJD(4)) * t174, 0.2e1 * (-t101 * t262 - t57) * t177 + 0.2e1 * (qJD(3) * t81 + t101 * t259 - t64 * t176) * t174, 0.2e1 * t101 * t64 + 0.2e1 * t57 * t81 + 0.2e1 * t63 * t82, -0.2e1 * t109 * t49, 0.2e1 * t108 * t49 - 0.2e1 * t109 * t50, -0.2e1 * t109 * t164 - 0.2e1 * t177 * t49, 0.2e1 * t108 * t164 - 0.2e1 * t177 * t50, t150, 0.2e1 * t108 * t48 + 0.2e1 * t16 * t177 - 0.2e1 * t164 * t40 + 0.2e1 * t50 * t78, 0.2e1 * t109 * t48 + 0.2e1 * t15 * t177 + 0.2e1 * t164 * t41 - 0.2e1 * t49 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t74, t239, t33, t32, t163 * t77 + t277, qJD(4) * t208 - t173 * t43 + t276, t195, -t194, 0, -pkin(3) * t43 - t192 + t198, -pkin(3) * t44 + t199 + t300, t111 * t76 + t130 * t43 - t192 + t200 (t277 - t279 + (t278 + t280) * qJD(4)) * pkin(10) + t182, -t111 * t77 - t130 * t44 + t201 - t300, pkin(10) * t182 + t27 * t111 + t14 * t130, t126 * t18 - t47 * t72, -t125 * t18 - t126 * t17 + t297 * t47 + t46 * t72, t112 * t72 - t126 * t74, -t112 * t297 + t125 * t74, 0, t112 * t54 + t122 * t17 + t125 * t7 - t22 * t297 + t46 * t90 - t74 * t79, -t112 * t53 + t122 * t18 + t126 * t7 - t22 * t72 + t47 * t90 + t74 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, -t164, 0, -pkin(9) * t261, pkin(9) * t164, -t174 * t225 + t220, t233 * t301 - t261 * t267, t302, t189, 0, t149 + (-pkin(3) * t176 + t286) * t258 + (t173 * t214 - t159) * qJD(3) (pkin(9) * t270 + t173 * t213) * qJD(4) + (t176 * t214 + t158) * qJD(3), t149 + (t130 * t258 - t64) * t176 - t295 * t173, t181 (-t64 + (t130 * t174 + t285) * qJD(4)) * t173 + t295 * t176, pkin(10) * t181 + t101 * t111 + t64 * t130, -t109 * t72 - t126 * t49, t108 * t72 + t109 * t297 + t125 * t49 - t126 * t50, -t126 * t164 - t177 * t72, t125 * t164 + t177 * t297, 0, t108 * t90 + t122 * t50 + t125 * t48 - t164 * t79 - t177 * t54 - t297 * t78, t109 * t90 - t122 * t49 + t126 * t48 + t164 * t80 + t177 * t53 - t72 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t233, -0.2e1 * t225, 0, 0, 0, t173 * t251, t176 * t251, -0.2e1 * t111 * t176 + 0.2e1 * t130 * t259, 0, -0.2e1 * t111 * t173 - 0.2e1 * t130 * t163, 0.2e1 * t130 * t111, -0.2e1 * t126 * t72, 0.2e1 * t125 * t72 + 0.2e1 * t126 * t297, 0, 0, 0, -0.2e1 * t122 * t297 + 0.2e1 * t125 * t90, -0.2e1 * t122 * t72 + 0.2e1 * t126 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, t74, -t227, t12, -t227 + 0.2e1 * t289, -pkin(4) * t44 - qJ(5) * t43 - qJD(5) * t76, 0.2e1 * t104 - t12 + 0.2e1 * t71, -pkin(4) * t11 + qJ(5) * t8 + qJD(5) * t23, 0, 0, -t18, t17, t74, t103 * t112 - t128 * t74 - t2, -t102 * t112 + t129 * t74 - t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t163 * t174 - t234, t164, -t222, t65, -t222 + 0.2e1 * t249 (-pkin(4) * t261 - qJ(5) * t258) * t176 + (-qJ(5) * t261 + (pkin(4) * qJD(4) - qJD(5)) * t174) * t173, 0.2e1 * t161 - t65 - 0.2e1 * t253, -pkin(4) * t63 + qJ(5) * t57 + qJD(5) * t81, 0, 0, t49, t50, t164, -t103 * t177 - t128 * t164 - t16, t102 * t177 + t129 * t164 - t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, -t259, 0, -t247, t248, -t247, -t294, -t248, -t294 * pkin(10), 0, 0, t72, -t297, 0, t54, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292, qJ(5) * t292, 0, 0, 0, 0, 0, 0.2e1 * t103, -0.2e1 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t44, 0, t11, 0, 0, 0, 0, 0, t112 * t256 - t288 * t74, t112 * t228 + t172 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t114, 0, t63, 0, 0, 0, 0, 0, -t174 * t230 - t177 * t256, t164 * t172 - t177 * t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, 0, t247, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, -t74, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t50, -t164, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t297, 0, -t54, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256, -t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
