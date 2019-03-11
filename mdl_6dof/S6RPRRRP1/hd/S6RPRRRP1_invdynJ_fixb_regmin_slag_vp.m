% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:10
% EndTime: 2019-03-09 05:57:21
% DurationCPUTime: 4.62s
% Computational Cost: add. (7146->441), mult. (14776->552), div. (0->0), fcn. (10279->14), ass. (0->245)
t184 = sin(qJ(4));
t188 = cos(qJ(3));
t328 = cos(qJ(4));
t259 = t328 * t188;
t185 = sin(qJ(3));
t273 = t185 * qJDD(1);
t284 = t184 * t188;
t130 = t328 * t185 + t284;
t271 = qJD(3) + qJD(4);
t94 = t271 * t130;
t77 = qJD(1) * t94 - qJDD(1) * t259 + t184 * t273;
t74 = qJDD(5) + t77;
t177 = qJ(1) + pkin(10);
t170 = sin(t177);
t171 = cos(t177);
t241 = g(1) * t171 + g(2) * t170;
t187 = cos(qJ(5));
t276 = qJD(5) * t187;
t256 = qJD(1) * t328;
t279 = qJD(1) * t185;
t118 = t184 * t279 - t188 * t256;
t294 = t118 * t187;
t349 = t276 + t294;
t183 = sin(qJ(5));
t277 = qJD(5) * t183;
t348 = t118 * t183 + t277;
t181 = sin(pkin(10));
t160 = pkin(1) * t181 + pkin(7);
t321 = pkin(8) + t160;
t180 = qJ(3) + qJ(4);
t174 = sin(t180);
t289 = t171 * t174;
t291 = t170 * t174;
t344 = g(1) * t289 + g(2) * t291;
t182 = cos(pkin(10));
t161 = -t182 * pkin(1) - pkin(2);
t176 = t188 * pkin(3);
t335 = t161 - t176;
t343 = -pkin(9) * t130 + t335;
t175 = cos(t180);
t342 = t241 * t175;
t117 = qJD(5) + t118;
t250 = t321 * qJD(1);
t105 = t185 * qJD(2) + t188 * t250;
t255 = qJD(4) * t328;
t278 = qJD(4) * t184;
t172 = t188 * qJDD(2);
t141 = t160 * qJDD(1);
t249 = pkin(8) * qJDD(1) + t141;
t68 = qJDD(3) * pkin(3) - t105 * qJD(3) - t249 * t185 + t172;
t104 = t188 * qJD(2) - t250 * t185;
t75 = t104 * qJD(3) + t185 * qJDD(2) + t249 * t188;
t314 = qJD(3) * pkin(3);
t99 = t104 + t314;
t207 = -t105 * t278 + t184 * t68 + t99 * t255 + t328 * t75;
t270 = qJDD(3) + qJDD(4);
t13 = pkin(9) * t270 + t207;
t274 = qJD(1) * qJD(3);
t254 = t185 * t274;
t158 = pkin(3) * t254;
t220 = -t184 * t185 + t259;
t93 = t271 * t220;
t198 = t93 * qJD(1);
t26 = t77 * pkin(4) - pkin(9) * t198 + t343 * qJDD(1) + t158;
t98 = t328 * t105;
t63 = t184 * t99 + t98;
t57 = pkin(9) * t271 + t63;
t120 = -qJD(1) * t284 - t185 * t256;
t121 = t335 * qJD(1);
t79 = t118 * pkin(4) + t120 * pkin(9) + t121;
t223 = t187 * t13 + t183 * t26 + t79 * t276 - t277 * t57;
t304 = t74 * qJ(6);
t2 = qJD(6) * t117 + t223 + t304;
t252 = t183 * t13 - t187 * t26 + t57 * t276 + t79 * t277;
t330 = pkin(5) * t74;
t4 = qJDD(6) + t252 - t330;
t341 = t4 * t183 + t2 * t187;
t307 = t187 * t74;
t338 = (t117 * t277 - t307) * pkin(9);
t247 = t187 * t271;
t100 = -t183 * t120 - t247;
t195 = t130 * qJDD(1) + t198;
t218 = t187 * t120 - t183 * t271;
t43 = -qJD(5) * t218 + t195 * t183 - t187 * t270;
t337 = -t94 * t100 + t220 * t43;
t69 = t184 * t104 + t98;
t244 = pkin(3) * t278 - t69;
t336 = pkin(5) * t348 - qJ(6) * t349 - qJD(6) * t183;
t334 = t175 * pkin(4) + t174 * pkin(9);
t333 = t130 * t270 + t271 * t93;
t332 = t218 ^ 2;
t331 = t117 ^ 2;
t329 = pkin(9) * t74;
t327 = pkin(5) * t120;
t164 = g(3) * t174;
t323 = g(3) * t175;
t251 = t105 * t255 + t184 * t75 + t99 * t278 - t328 * t68;
t14 = -pkin(4) * t270 + t251;
t42 = -qJD(5) * t247 - t120 * t277 - t183 * t270 - t187 * t195;
t5 = t43 * pkin(5) + t42 * qJ(6) + qJD(6) * t218 + t14;
t322 = t5 * t183;
t292 = t130 * t187;
t306 = t187 * t93;
t320 = -t100 * t306 - t43 * t292;
t97 = t184 * t105;
t62 = t328 * t99 - t97;
t90 = -pkin(4) * t120 + pkin(9) * t118;
t319 = t183 * t90 + t187 * t62;
t70 = t328 * t104 - t97;
t84 = pkin(3) * t279 + t90;
t318 = t183 * t84 + t187 * t70;
t317 = -t218 * t94 + t220 * t42;
t86 = -pkin(4) * t220 + t343;
t123 = t321 * t185;
t124 = t321 * t188;
t89 = -t184 * t123 + t328 * t124;
t316 = t183 * t86 + t187 * t89;
t315 = pkin(9) * qJD(5);
t31 = t183 * t79 + t187 * t57;
t313 = t117 * t31;
t56 = -pkin(4) * t271 - t62;
t32 = t100 * pkin(5) + qJ(6) * t218 + t56;
t312 = t118 * t32;
t311 = t118 * t56;
t166 = pkin(3) * t184 + pkin(9);
t310 = t166 * t74;
t309 = t183 * t74;
t308 = t183 * t93;
t305 = t42 * t183;
t303 = t336 + t244;
t302 = -t63 + t336;
t301 = t100 * t117;
t300 = t100 * t183;
t299 = t218 * t100;
t298 = t218 * t117;
t297 = t218 * t187;
t296 = t117 * t120;
t248 = t117 * t187;
t293 = t120 * t118;
t287 = t175 * t183;
t286 = t175 * t187;
t285 = t183 * qJ(6);
t30 = -t183 * t57 + t187 * t79;
t283 = qJD(6) - t30;
t282 = qJDD(2) - g(3);
t281 = t344 * t187;
t178 = t185 ^ 2;
t280 = -t188 ^ 2 + t178;
t144 = qJD(1) * t161;
t272 = t188 * qJDD(1);
t269 = t328 * pkin(3);
t268 = t185 * t314;
t266 = t218 * t308;
t52 = t74 * t292;
t27 = t32 * t277;
t265 = t32 * t276;
t53 = t56 * t277;
t264 = t164 + t342;
t263 = g(3) * t287 - t344 * t183;
t262 = -t5 - t323;
t261 = t183 * t328;
t260 = t187 * t328;
t258 = t130 * t277;
t257 = -t14 - t323;
t253 = qJD(3) * t321;
t246 = pkin(3) * t255;
t245 = pkin(5) * t286 + t175 * t285 + t334;
t108 = t170 * t287 + t171 * t187;
t110 = -t170 * t187 + t171 * t287;
t243 = -g(1) * t108 + g(2) * t110;
t109 = t170 * t286 - t171 * t183;
t111 = t170 * t183 + t171 * t286;
t242 = g(1) * t109 - g(2) * t111;
t240 = g(1) * t170 - g(2) * t171;
t186 = sin(qJ(1));
t189 = cos(qJ(1));
t239 = g(1) * t186 - g(2) * t189;
t237 = t187 * pkin(5) + t285;
t236 = pkin(5) * t183 - qJ(6) * t187;
t234 = -t310 + t311;
t20 = -pkin(5) * t117 + t283;
t21 = qJ(6) * t117 + t31;
t233 = -t183 * t21 + t187 * t20;
t232 = t183 * t20 + t187 * t21;
t229 = -t20 * t120 + t27 + t281;
t228 = t30 * t120 + t281 + t53;
t227 = pkin(4) + t237;
t226 = t176 + pkin(2) + t334;
t224 = t258 - t306;
t114 = t185 * t253;
t115 = t188 * t253;
t221 = -t328 * t123 - t184 * t124;
t45 = t221 * qJD(4) - t328 * t114 - t184 * t115;
t51 = pkin(4) * t94 - pkin(9) * t93 + t268;
t222 = t183 * t51 + t187 * t45 + t86 * t276 - t277 * t89;
t219 = -t117 * t258 + t93 * t248 + t52;
t217 = -t31 * t120 + t14 * t183 + t56 * t276 + t263;
t215 = -qJD(1) * t144 - t141 + t241;
t214 = t21 * t120 - t32 * t294 - t263 - t322;
t213 = 0.2e1 * qJD(3) * t144 - qJDD(3) * t160;
t212 = g(1) * t110 + g(2) * t108 + t183 * t164 - t252;
t211 = -t166 * t277 + t187 * t246;
t210 = -t305 + (-t297 + t300) * qJD(5);
t209 = t121 * t120 - t251 - t323 + t344;
t191 = qJD(3) ^ 2;
t208 = -0.2e1 * qJDD(1) * t161 - t160 * t191 + t240;
t206 = (-t130 * t276 - t308) * t117 - t130 * t309;
t205 = qJD(5) * t233 + t341;
t204 = -t187 * t43 + t210;
t203 = -t218 * t32 + qJDD(6) - t212;
t202 = t20 * t349 - t21 * t348 - t264 + t341;
t201 = t206 - t337;
t200 = -g(1) * t111 - g(2) * t109 - t187 * t164 + t223;
t46 = t89 * qJD(4) - t184 * t114 + t328 * t115;
t199 = t121 * t118 - t207 + t264;
t197 = t241 * t174 * t227 - t342 * pkin(9);
t192 = qJD(1) ^ 2;
t190 = -pkin(8) - pkin(7);
t167 = -t269 - pkin(4);
t140 = qJDD(3) * t188 - t185 * t191;
t139 = qJDD(3) * t185 + t188 * t191;
t122 = -t269 - t227;
t113 = t120 * qJ(6);
t106 = qJDD(1) * t335 + t158;
t80 = -t118 ^ 2 + t120 ^ 2;
t76 = t220 * t270 - t271 * t94;
t66 = -pkin(5) * t218 + qJ(6) * t100;
t50 = t130 * t236 - t221;
t49 = -t120 * t271 - t77;
t48 = t118 * t271 + t195;
t40 = pkin(5) * t220 + t183 * t89 - t187 * t86;
t39 = -qJ(6) * t220 + t316;
t34 = t183 * t62 - t187 * t90 + t327;
t33 = -t113 + t319;
t29 = t183 * t70 - t187 * t84 + t327;
t28 = -t113 + t318;
t23 = -t42 + t301;
t18 = t117 * t248 - t120 * t218 + t309;
t17 = -t100 * t120 - t331 * t183 + t307;
t15 = -t218 * t248 - t305;
t9 = t236 * t93 + (qJD(5) * t237 - qJD(6) * t187) * t130 + t46;
t8 = -t94 * pkin(5) + t316 * qJD(5) + t183 * t45 - t187 * t51;
t7 = qJ(6) * t94 - qJD(6) * t220 + t222;
t6 = (-t42 - t301) * t187 + (-t43 + t298) * t183;
t1 = [qJDD(1), t239, g(1) * t189 + g(2) * t186 (t239 + (t181 ^ 2 + t182 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t178 + 0.2e1 * t188 * t254, 0.2e1 * t185 * t272 - 0.2e1 * t274 * t280, t139, t140, 0, t185 * t213 + t188 * t208, -t185 * t208 + t188 * t213, -t120 * t93 + t195 * t130, -t93 * t118 + t120 * t94 - t130 * t77 + t195 * t220, t333, t76, 0, -t106 * t220 + t118 * t268 + t121 * t94 + t175 * t240 + t221 * t270 - t271 * t46 + t335 * t77, -g(1) * t291 + g(2) * t289 + t106 * t130 - t120 * t268 + t121 * t93 + t195 * t335 - t270 * t89 - t271 * t45, t218 * t224 - t42 * t292, t266 + (t305 + (t297 + t300) * qJD(5)) * t130 + t320, t219 + t317, t206 + t337, t117 * t94 - t220 * t74, t252 * t220 + t30 * t94 + t46 * t100 - t221 * t43 + ((-qJD(5) * t89 + t51) * t117 + t86 * t74 + t56 * qJD(5) * t130) * t187 + ((-qJD(5) * t86 - t45) * t117 - t89 * t74 + t14 * t130 + t56 * t93) * t183 + t242, -t222 * t117 - t316 * t74 + t223 * t220 - t31 * t94 - t46 * t218 + t221 * t42 + t56 * t306 + (t14 * t187 - t53) * t130 + t243, t32 * t308 + t9 * t100 - t8 * t117 + t4 * t220 - t20 * t94 - t40 * t74 + t50 * t43 + (t265 + t322) * t130 + t242, -t7 * t100 - t8 * t218 - t39 * t43 - t40 * t42 + t233 * t93 + t240 * t174 + (-qJD(5) * t232 - t183 * t2 + t187 * t4) * t130, -t32 * t306 + t9 * t218 + t7 * t117 - t2 * t220 + t21 * t94 + t39 * t74 + t50 * t42 + (-t5 * t187 + t27) * t130 - t243, t2 * t39 + t21 * t7 + t5 * t50 + t32 * t9 + t4 * t40 + t20 * t8 - g(1) * (-t186 * pkin(1) - t109 * pkin(5) - t108 * qJ(6)) - g(2) * (t189 * pkin(1) + t111 * pkin(5) + t110 * qJ(6)) + (g(1) * t190 - g(2) * t226) * t171 + (g(1) * t226 + g(2) * t190) * t170; 0, 0, 0, t282, 0, 0, 0, 0, 0, t140, -t139, 0, 0, 0, 0, 0, t76, -t333, 0, 0, 0, 0, 0, t201, t117 * t224 + t317 - t52, t201, t130 * t210 - t266 + t320, t219 - t317, t130 * t205 - t220 * t5 + t232 * t93 + t32 * t94 - g(3); 0, 0, 0, 0, -t185 * t192 * t188, t280 * t192, t273, t272, qJDD(3), -g(3) * t188 + t185 * t215 + t172, -t282 * t185 + t215 * t188, -t293, t80, t48, t49, t270, t69 * t271 + (-t118 * t279 + t328 * t270 - t271 * t278) * pkin(3) + t209, t70 * t271 + (t120 * t279 - t184 * t270 - t255 * t271) * pkin(3) + t199, t15, t6, t18, t17, t296, t167 * t43 + t257 * t187 + t234 * t183 + t244 * t100 + ((-qJD(5) * t166 - t84) * t187 + (-t246 + t70) * t183) * t117 + t228, -t167 * t42 + t234 * t187 - t244 * t218 + (-t211 + t318) * t117 + t217, t122 * t43 + t262 * t187 + (-t310 + t312) * t183 + t303 * t100 + (-t166 * t276 - t183 * t246 + t29) * t117 + t229, t28 * t100 + t29 * t218 + (-t100 * t260 - t218 * t261) * qJD(4) * pkin(3) + t204 * t166 + t202, t122 * t42 + (-qJD(5) * t32 + t310) * t187 + t303 * t218 + (-t28 + t211) * t117 + t214, t5 * t122 - t21 * t28 - t20 * t29 - g(3) * (t176 + t245) + t303 * t32 + (t241 * t185 + (t20 * t261 + t21 * t260) * qJD(4)) * pkin(3) + t205 * t166 + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t293, t80, t48, t49, t270, t271 * t63 + t209, t271 * t62 + t199, t15, t6, t18, t17, t296, -pkin(4) * t43 - t63 * t100 + (t117 * t62 + t311 - t329) * t183 + ((-t90 - t315) * t117 + t257) * t187 + t228, pkin(4) * t42 + t319 * t117 + t218 * t63 + t56 * t294 + t217 + t338, t34 * t117 - t227 * t43 + (t312 - t329) * t183 + t302 * t100 + (-t117 * t315 + t262) * t187 + t229, pkin(9) * t204 + t33 * t100 + t218 * t34 + t202, -t33 * t117 + t218 * t302 - t227 * t42 + t214 - t265 - t338, pkin(9) * t205 - g(3) * t245 - t20 * t34 - t21 * t33 - t227 * t5 + t302 * t32 + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, -t100 ^ 2 + t332, t23, -t43 - t298, t74, t218 * t56 + t212 + t313, t100 * t56 + t117 * t30 - t200, -t100 * t66 - t203 + t313 + 0.2e1 * t330, pkin(5) * t42 - t43 * qJ(6) - (t21 - t31) * t218 + (t20 - t283) * t100, 0.2e1 * t304 - t32 * t100 - t66 * t218 + (0.2e1 * qJD(6) - t30) * t117 + t200, t2 * qJ(6) - t4 * pkin(5) - t32 * t66 - t20 * t31 - g(1) * (-pkin(5) * t110 + qJ(6) * t111) - g(2) * (-pkin(5) * t108 + qJ(6) * t109) + t283 * t21 + t236 * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299 - t74, t23, -t331 - t332, -t117 * t21 + t203 - t330;];
tau_reg  = t1;
