% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:09
% EndTime: 2019-03-08 21:54:20
% DurationCPUTime: 4.62s
% Computational Cost: add. (5500->404), mult. (13203->565), div. (0->0), fcn. (11023->16), ass. (0->233)
t195 = cos(qJ(6));
t264 = qJD(6) * t195;
t185 = sin(pkin(12));
t188 = cos(pkin(12));
t193 = sin(qJ(3));
t197 = cos(qJ(3));
t151 = -t185 * t193 + t188 * t197;
t144 = t151 * qJD(2);
t196 = cos(qJ(5));
t131 = t196 * t144;
t152 = t185 * t197 + t188 * t193;
t146 = t152 * qJD(2);
t192 = sin(qJ(5));
t87 = -t146 * t192 + t131;
t334 = t195 * t87;
t341 = t264 - t334;
t182 = qJD(3) + qJD(5);
t294 = t182 * t87;
t145 = t152 * qJD(3);
t100 = -qJD(2) * t145 + t151 * qJDD(2);
t262 = qJD(2) * qJD(3);
t249 = t197 * t262;
t250 = t193 * t262;
t101 = t152 * qJDD(2) - t185 * t250 + t188 * t249;
t266 = qJD(5) * t192;
t36 = qJD(5) * t131 + t192 * t100 + t196 * t101 - t146 * t266;
t340 = t36 - t294;
t303 = qJ(4) + pkin(8);
t243 = qJD(3) * t303;
t137 = t197 * qJD(4) - t193 * t243;
t138 = -t193 * qJD(4) - t197 * t243;
t187 = sin(pkin(6));
t198 = cos(qJ(2));
t275 = t187 * t198;
t254 = qJD(1) * t275;
t285 = -t137 * t185 + t188 * t138 + t152 * t254;
t284 = t188 * t137 + t185 * t138 - t151 * t254;
t273 = -qJD(6) + t87;
t339 = t273 + qJD(6);
t191 = sin(qJ(6));
t181 = qJDD(3) + qJDD(5);
t226 = t144 * t192 + t196 * t146;
t265 = qJD(6) * t191;
t24 = t191 * t181 + t182 * t264 + t195 * t36 - t226 * t265;
t21 = t24 * t195;
t73 = t182 * t191 + t195 * t226;
t25 = qJD(6) * t73 - t195 * t181 + t191 * t36;
t71 = -t195 * t182 + t191 * t226;
t338 = -t191 * t25 - t341 * t71 + t21;
t20 = t24 * t191;
t337 = t341 * t73 + t20;
t37 = t226 * qJD(5) - t196 * t100 + t192 * t101;
t34 = qJDD(6) + t37;
t31 = t191 * t34;
t299 = -t264 * t273 + t31;
t304 = t73 * t226;
t336 = t273 * t334 + t299 - t304;
t309 = pkin(9) * t146;
t194 = sin(qJ(2));
t268 = qJD(1) * t194;
t255 = t187 * t268;
t233 = t303 * qJD(2) + t255;
t189 = cos(pkin(6));
t269 = qJD(1) * t189;
t112 = t193 * t269 + t233 * t197;
t104 = t185 * t112;
t111 = -t233 * t193 + t197 * t269;
t297 = qJD(3) * pkin(3);
t108 = t111 + t297;
t59 = t188 * t108 - t104;
t43 = qJD(3) * pkin(4) - t309 + t59;
t310 = pkin(9) * t144;
t274 = t188 * t112;
t60 = t185 * t108 + t274;
t46 = t60 + t310;
t22 = -t192 * t46 + t196 * t43;
t17 = -pkin(5) * t182 - t22;
t335 = t17 * t87;
t283 = cos(pkin(11));
t241 = t283 * t194;
t186 = sin(pkin(11));
t278 = t186 * t198;
t141 = t189 * t241 + t278;
t240 = t283 * t198;
t279 = t186 * t194;
t143 = -t189 * t279 + t240;
t180 = qJ(3) + pkin(12) + qJ(5);
t174 = sin(t180);
t175 = cos(t180);
t242 = t187 * t283;
t277 = t187 * t194;
t280 = t186 * t187;
t215 = -g(3) * (-t174 * t277 + t175 * t189) - g(2) * (-t141 * t174 - t175 * t242) - g(1) * (-t143 * t174 + t175 * t280);
t261 = t189 * qJDD(1);
t167 = t197 * t261;
t263 = qJD(1) * qJD(2);
t123 = qJDD(2) * pkin(8) + (qJDD(1) * t194 + t198 * t263) * t187;
t206 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t269 + t123;
t223 = t233 * qJD(3);
t54 = qJDD(3) * pkin(3) - t206 * t193 - t197 * t223 + t167;
t55 = (-t223 + t261) * t193 + t206 * t197;
t29 = -t185 * t55 + t188 * t54;
t13 = qJDD(3) * pkin(4) - pkin(9) * t101 + t29;
t30 = t185 * t54 + t188 * t55;
t14 = pkin(9) * t100 + t30;
t23 = t192 * t43 + t196 * t46;
t318 = -t23 * qJD(5) + t196 * t13 - t192 * t14;
t3 = -pkin(5) * t181 - t318;
t214 = t215 - t3;
t333 = t226 * t87;
t148 = t151 * qJD(3);
t332 = pkin(9) * t148 - t285;
t331 = pkin(9) * t145 - t284;
t330 = t191 * t273;
t295 = t182 * t226;
t329 = -t37 + t295;
t140 = -t189 * t240 + t279;
t142 = t189 * t278 + t241;
t231 = g(1) * t142 + g(2) * t140;
t213 = g(3) * t275 - t231;
t208 = t213 * t175;
t103 = t151 * t192 + t152 * t196;
t177 = pkin(3) * t197 + pkin(2);
t125 = -pkin(4) * t151 - t177;
t225 = t196 * t151 - t152 * t192;
t39 = -pkin(5) * t225 - pkin(10) * t103 + t125;
t327 = t39 * t34 - t208;
t326 = t226 ^ 2 - t87 ^ 2;
t52 = pkin(5) * t226 - pkin(10) * t87;
t121 = t174 * t189 + t175 * t277;
t315 = (qJD(5) * t43 + t14) * t196 + t192 * t13 - t46 * t266;
t93 = t141 * t175 - t174 * t242;
t95 = t143 * t175 + t174 * t280;
t134 = -t177 * qJD(2) + qJD(4) - t254;
t96 = -pkin(4) * t144 + t134;
t325 = g(1) * t95 + g(2) * t93 + g(3) * t121 - t96 * t87 - t315;
t305 = t71 * t226;
t322 = t273 * t226;
t32 = t195 * t34;
t219 = -t265 * t273 - t32;
t179 = t193 * t297;
t228 = pkin(4) * t145 + t179 - t255;
t311 = pkin(3) * t188;
t176 = pkin(4) + t311;
t312 = pkin(3) * t185;
t271 = t192 * t176 + t196 * t312;
t18 = pkin(10) * t182 + t23;
t35 = -pkin(5) * t87 - pkin(10) * t226 + t96;
t227 = t18 * t191 - t195 * t35;
t321 = t17 * t265 + t226 * t227;
t5 = t18 * t195 + t191 * t35;
t320 = t17 * t264 - t214 * t191 + t5 * t226;
t319 = -t226 * t96 + t215 + t318;
t2 = pkin(10) * t181 + t315;
t230 = g(1) * t143 + g(2) * t141;
t212 = -g(3) * t277 - t230;
t159 = t303 * t193;
t160 = t303 * t197;
t116 = -t188 * t159 - t160 * t185;
t77 = -pkin(9) * t152 + t116;
t117 = -t185 * t159 + t188 * t160;
t78 = pkin(9) * t151 + t117;
t41 = t192 * t78 - t196 * t77;
t302 = t41 * qJD(5) + t332 * t192 + t331 * t196;
t42 = t192 * t77 + t196 * t78;
t56 = t225 * qJD(5) - t145 * t192 + t148 * t196;
t317 = (qJD(6) * t35 + t2) * t225 + t17 * t56 + t3 * t103 - (-qJD(6) * t39 + t302) * t273 - t42 * t34 + t212;
t251 = t194 * t263;
t165 = t187 * t251;
t199 = qJD(3) ^ 2;
t247 = qJDD(1) * t275;
t316 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t199 + (-g(3) * t198 + t251) * t187 - t165 + t231 + t247;
t149 = t189 * t197 - t193 * t277;
t307 = g(3) * t149;
t301 = t42 * qJD(5) - t331 * t192 + t332 * t196;
t298 = qJD(2) * pkin(2);
t296 = t103 * t17;
t293 = t191 * t73;
t291 = t195 * t73;
t290 = t195 * t273;
t221 = t176 * t196 - t192 * t312;
t62 = -t111 * t185 - t274;
t49 = t62 - t310;
t64 = t188 * t111 - t104;
t50 = t64 - t309;
t287 = -t221 * qJD(5) + t192 * t49 + t196 * t50;
t286 = t271 * qJD(5) - t192 * t50 + t196 * t49;
t276 = t187 * t197;
t272 = qJDD(1) - g(3);
t183 = t193 ^ 2;
t270 = -t197 ^ 2 + t183;
t267 = qJD(2) * t194;
t260 = t197 * qJDD(2);
t178 = t193 * qJD(2) * pkin(3);
t257 = t191 * t275;
t256 = t195 * t275;
t253 = t187 * t267;
t252 = qJD(2) * t275;
t248 = t198 * t262;
t118 = pkin(4) * t146 + t178;
t236 = t187 * t272;
t136 = pkin(10) + t271;
t234 = qJD(6) * t136 + t118 + t52;
t57 = t103 * qJD(5) + t196 * t145 + t148 * t192;
t229 = pkin(5) * t57 - pkin(10) * t56 + t228;
t150 = t189 * t193 + t194 * t276;
t82 = t149 * t188 - t150 * t185;
t83 = t149 * t185 + t150 * t188;
t44 = t192 * t83 - t196 * t82;
t45 = t192 * t82 + t196 * t83;
t224 = -t330 * t87 - t219;
t200 = qJD(2) ^ 2;
t222 = qJDD(2) * t198 - t194 * t200;
t220 = -g(1) * t186 + t283 * g(2);
t218 = -t191 * t45 - t256;
t217 = -t195 * t45 + t257;
t158 = -t254 - t298;
t211 = -qJD(2) * t158 - t123 + t230;
t210 = -t136 * t34 - t273 * t287 - t335;
t209 = pkin(3) * t250 - t177 * qJDD(2) + qJDD(4) + t165;
t204 = -pkin(8) * qJDD(3) + (t158 + t254 - t298) * qJD(3);
t99 = t209 - t247;
t58 = -pkin(4) * t100 + t99;
t135 = -pkin(5) - t221;
t110 = -qJD(3) * t150 - t193 * t252;
t109 = qJD(3) * t149 + t197 * t252;
t63 = t109 * t188 + t110 * t185;
t61 = -t109 * t185 + t110 * t188;
t9 = t45 * qJD(5) + t192 * t63 - t196 * t61;
t8 = -t44 * qJD(5) + t192 * t61 + t196 * t63;
t7 = pkin(5) * t37 - pkin(10) * t36 + t58;
t6 = t195 * t7;
t1 = [t272, 0, t222 * t187 (-qJDD(2) * t194 - t198 * t200) * t187, 0, 0, 0, 0, 0, qJD(3) * t110 + qJDD(3) * t149 + (-t193 * t248 + t222 * t197) * t187, -qJD(3) * t109 - qJDD(3) * t150 + (-t193 * t222 - t197 * t248) * t187, t100 * t83 - t101 * t82 + t144 * t63 - t146 * t61, t29 * t82 + t30 * t83 + t59 * t61 + t60 * t63 - g(3) + (t134 * t267 - t198 * t99) * t187, 0, 0, 0, 0, 0, -t181 * t44 - t182 * t9 + (-t198 * t37 - t267 * t87) * t187, -t181 * t45 - t182 * t8 + (-t198 * t36 + t226 * t267) * t187, 0, 0, 0, 0, 0 -(qJD(6) * t217 - t191 * t8 + t195 * t253) * t273 + t218 * t34 + t9 * t71 + t44 * t25 (qJD(6) * t218 + t191 * t253 + t195 * t8) * t273 + t217 * t34 + t9 * t73 + t44 * t24; 0, qJDD(2), t272 * t275 + t231, -t194 * t236 + t230, qJDD(2) * t183 + 0.2e1 * t193 * t249, 0.2e1 * t193 * t260 - 0.2e1 * t270 * t262, qJDD(3) * t193 + t197 * t199, qJDD(3) * t197 - t193 * t199, 0, t204 * t193 + t316 * t197, -t316 * t193 + t204 * t197, t100 * t117 - t101 * t116 + t284 * t144 - t145 * t60 - t285 * t146 - t148 * t59 + t151 * t30 - t152 * t29 + t212, t30 * t117 + t29 * t116 - t99 * t177 + t134 * t179 - g(1) * (-t142 * t177 + t143 * t303) - g(2) * (-t140 * t177 + t141 * t303) + t284 * t60 + t285 * t59 + (-t134 * t268 - g(3) * (t177 * t198 + t194 * t303)) * t187, t103 * t36 + t226 * t56, -t103 * t37 + t225 * t36 - t226 * t57 + t56 * t87, t103 * t181 + t182 * t56, t181 * t225 - t182 * t57, 0, t125 * t37 - t181 * t41 - t301 * t182 - t225 * t58 - t228 * t87 + t57 * t96 - t208, t103 * t58 + t125 * t36 + t213 * t174 - t181 * t42 + t302 * t182 + t226 * t228 + t56 * t96, t56 * t291 + (-t265 * t73 + t21) * t103 (-t195 * t71 - t293) * t56 + (-t20 - t195 * t25 + (t191 * t71 - t291) * qJD(6)) * t103, -t103 * t219 - t225 * t24 - t56 * t290 + t57 * t73, -t299 * t103 + t225 * t25 + t330 * t56 - t57 * t71, -t225 * t34 - t273 * t57, -t6 * t225 + t41 * t25 - t227 * t57 + t301 * t71 + (-t229 * t273 + (t18 * t225 + t273 * t42 + t296) * qJD(6) + t327) * t195 + t317 * t191, t41 * t24 - t5 * t57 + t301 * t73 + ((-qJD(6) * t18 + t7) * t225 - qJD(6) * t296 - (qJD(6) * t42 - t229) * t273 - t327) * t191 + t317 * t195; 0, 0, 0, 0, -t193 * t200 * t197, t270 * t200, t193 * qJDD(2), t260, qJDD(3), t193 * t211 + t220 * t276 + t167 - t307, g(3) * t150 + (-t187 * t220 - t261) * t193 + t211 * t197 (t60 + t62) * t146 + (t59 - t64) * t144 + (t100 * t185 - t101 * t188) * pkin(3), -t134 * t178 + t29 * t311 + t30 * t312 - t59 * t62 - t60 * t64 + (-g(1) * (-t143 * t193 + t186 * t276) - g(2) * (-t141 * t193 - t197 * t242) - t307) * pkin(3), -t333, t326, t340, t329, t181, t118 * t87 + t221 * t181 - t286 * t182 + t319, -t118 * t226 - t271 * t181 + t287 * t182 + t325, t337, t273 * t293 + t338, t336, t224 + t305, t322, t135 * t25 + t286 * t71 + t210 * t191 + (t234 * t273 + t214) * t195 + t321, t135 * t24 + t210 * t195 - t234 * t330 + t286 * t73 + t320; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144 ^ 2 - t146 ^ 2, -t144 * t60 + t146 * t59 - t198 * t236 + t209 - t231, 0, 0, 0, 0, 0, t37 + t295, t36 + t294, 0, 0, 0, 0, 0, t224 - t305, -t273 * t290 - t304 - t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t333, t326, t340, t329, t181, t182 * t23 + t319, t182 * t22 + t325, t337, t330 * t73 + t338, t336, -t273 * t330 + t305 + t32, t322, -pkin(5) * t25 - t23 * t71 + (-pkin(10) * t34 - t22 * t273 - t335) * t191 + (-(-pkin(10) * qJD(6) - t52) * t273 + t214) * t195 + t321, -pkin(5) * t24 - (t191 * t52 + t195 * t22) * t273 - t23 * t73 - t17 * t334 + t219 * pkin(10) + t320; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, -t273 * t71 + t24, -t273 * t73 - t25, t34, -t191 * t2 + t6 - t17 * t73 - g(1) * (t142 * t195 - t191 * t95) - g(2) * (t140 * t195 - t191 * t93) - g(3) * (-t121 * t191 - t256) - t339 * t5, -t195 * t2 - t191 * t7 + t17 * t71 - g(1) * (-t142 * t191 - t195 * t95) - g(2) * (-t140 * t191 - t195 * t93) - g(3) * (-t121 * t195 + t257) + t339 * t227;];
tau_reg  = t1;
