% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPPR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:12:36
% EndTime: 2019-05-05 03:12:52
% DurationCPUTime: 6.53s
% Computational Cost: add. (12547->400), mult. (26713->540), div. (0->0), fcn. (18696->12), ass. (0->248)
t204 = sin(pkin(11));
t206 = cos(pkin(11));
t210 = sin(qJ(3));
t258 = qJD(2) * qJD(3);
t192 = t210 * t258;
t213 = cos(qJ(3));
t256 = t213 * qJDD(2);
t177 = -t192 + t256;
t263 = qJD(2) * t210;
t169 = -qJD(3) * t206 + t204 * t263;
t171 = qJD(3) * t204 + t206 * t263;
t282 = t171 * t169;
t229 = t177 - t282;
t271 = t206 * t229;
t167 = t171 ^ 2;
t201 = t213 ^ 2;
t215 = qJD(2) ^ 2;
t197 = t201 * t215;
t304 = -t167 - t197;
t103 = -t204 * t304 + t271;
t205 = sin(pkin(6));
t207 = cos(pkin(6));
t211 = sin(qJ(2));
t214 = cos(qJ(2));
t249 = t213 * t258;
t257 = t210 * qJDD(2);
t176 = t249 + t257;
t150 = qJDD(3) * t204 + t176 * t206;
t260 = t213 * qJD(2);
t251 = t169 * t260;
t222 = t150 + t251;
t71 = t103 * t213 + t210 * t222;
t276 = t204 * t229;
t91 = -t206 * t304 - t276;
t359 = (t211 * t71 + t214 * t91) * t205 + t207 * (t103 * t210 - t213 * t222);
t358 = pkin(2) * t91 + pkin(8) * t71;
t354 = pkin(3) * t91;
t297 = t169 ^ 2;
t114 = -t297 - t167;
t158 = t171 * t260;
t245 = -qJDD(3) * t206 + t176 * t204;
t121 = -t245 - t158;
t223 = -t150 + t251;
t329 = t121 * t206 - t204 * t223;
t342 = t114 * t210 + t213 * t329;
t353 = pkin(8) * t342;
t352 = qJ(4) * t91;
t350 = qJ(4) * t103;
t330 = t121 * t204 + t206 * t223;
t348 = t207 * (-t114 * t213 + t210 * t329) + (t211 * t342 - t214 * t330) * t205;
t316 = t158 - t245;
t347 = t316 * pkin(4);
t127 = t177 + t282;
t270 = t206 * t127;
t300 = -t297 - t197;
t312 = t204 * t300 - t270;
t275 = t204 * t127;
t311 = t206 * t300 + t275;
t325 = -t210 * t316 + t213 * t311;
t345 = -pkin(2) * t312 + pkin(8) * t325;
t344 = -pkin(3) * t114 + qJ(4) * t329;
t154 = t297 - t197;
t343 = -t213 * t121 + t210 * (t154 * t206 + t276);
t340 = t207 * (t210 * t311 + t213 * t316) + (t211 * t325 - t214 * t312) * t205;
t337 = pkin(3) * t312;
t336 = qJ(4) * t312;
t155 = -t167 + t197;
t332 = t155 * t206 - t275;
t331 = pkin(3) * t316 + qJ(4) * t311;
t328 = t154 * t204 - t271;
t273 = t206 * t316;
t278 = t204 * t222;
t326 = t210 * (-t273 + t278) + t213 * (t167 - t297);
t323 = t210 * (-t155 * t204 - t270) + t213 * t223;
t209 = sin(qJ(6));
t212 = cos(qJ(6));
t135 = -t169 * t212 + t171 * t209;
t137 = t169 * t209 + t171 * t212;
t100 = t137 * t135;
t172 = qJDD(6) + t177;
t305 = -t100 + t172;
t320 = t209 * t305;
t319 = t212 * t305;
t283 = sin(pkin(10));
t284 = cos(pkin(10));
t226 = g(1) * t283 - g(2) * t284;
t264 = -g(3) + qJDD(1);
t218 = -t205 * t226 + t207 * t264;
t217 = t210 * t218;
t181 = -g(1) * t284 - g(2) * t283;
t224 = t207 * t226;
t314 = t205 * t264 + t224;
t131 = t214 * t181 + t211 * t314;
t109 = -t215 * pkin(2) + qJDD(2) * pkin(8) + t131;
t238 = -pkin(3) * t213 - qJ(4) * t210;
t174 = t238 * qJD(2);
t243 = qJD(2) * t174 + t109;
t296 = qJD(3) ^ 2;
t76 = -pkin(3) * t296 + qJDD(3) * qJ(4) + t213 * t243 + t217;
t239 = t211 * t181 - t214 * t314;
t108 = -qJDD(2) * pkin(2) - pkin(8) * t215 + t239;
t237 = t176 + t249;
t85 = -t237 * qJ(4) + (-t177 + t192) * pkin(3) + t108;
t248 = t204 * t76 - t206 * t85;
t228 = pkin(4) * t177 - qJ(5) * t197 + qJDD(5) + t248;
t219 = pkin(9) * t223 + t228;
t140 = pkin(4) * t169 - qJ(5) * t171;
t295 = 2 * qJD(4);
t259 = t295 + t140;
t315 = t219 + (pkin(5) * t169 + t259) * t171 + t177 * pkin(5);
t189 = qJD(6) + t260;
t117 = t189 * t135;
t87 = -qJD(6) * t135 + t150 * t212 + t209 * t245;
t309 = -t117 + t87;
t306 = t243 * t210;
t292 = t204 * t85 + t206 * t76;
t301 = -qJ(5) * t177 - 0.2e1 * qJD(5) * t260 - t140 * t169 + t292;
t279 = t204 * t316;
t299 = t206 * t222 + t279;
t146 = t213 * t218;
t241 = -qJDD(3) * pkin(3) - qJ(4) * t296 + qJDD(4) - t146;
t227 = t150 * qJ(5) - t241 + t347;
t269 = t210 * t109;
t298 = -(qJ(5) * t169 * t213 - t174 * t210) * qJD(2) - t227 + t269;
t133 = t135 ^ 2;
t134 = t137 ^ 2;
t187 = t189 ^ 2;
t294 = pkin(4) + pkin(5);
t75 = t241 + t306;
t291 = t204 * t75;
t290 = t206 * t75;
t234 = pkin(5) * t260 - pkin(9) * t171;
t261 = qJD(5) * t171;
t47 = -0.2e1 * t261 + t298;
t37 = pkin(5) * t245 + pkin(9) * t297 - t171 * t234 + t47;
t289 = t209 * t37;
t89 = t100 + t172;
t288 = t209 * t89;
t262 = qJD(4) * t169;
t161 = -0.2e1 * t262;
t235 = t161 + t301;
t35 = -pkin(4) * t197 + t235;
t30 = -pkin(5) * t297 + pkin(9) * t245 - t234 * t260 + t35;
t287 = t212 * t30;
t286 = t212 * t37;
t285 = t212 * t89;
t281 = t189 * t209;
t280 = t189 * t212;
t188 = t210 * t215 * t213;
t183 = qJDD(3) + t188;
t267 = t210 * t183;
t182 = -t188 + qJDD(3);
t265 = t213 * t182;
t255 = t171 * t295;
t46 = t161 + t292;
t254 = t213 * t100;
t253 = t213 * t282;
t252 = pkin(4) * t206 + pkin(3);
t45 = t248 + t255;
t26 = t204 * t45 + t206 * t46;
t11 = t209 * t30 - t212 * t315;
t97 = -t146 + t269;
t98 = t213 * t109 + t217;
t50 = t210 * t97 + t213 * t98;
t246 = t209 * t150 - t212 * t245;
t242 = t204 * t251;
t240 = t210 * (t150 * t206 + t158 * t204) - t253;
t12 = t209 * t315 + t287;
t8 = -t11 * t212 + t12 * t209;
t9 = t209 * t11 + t212 * t12;
t25 = t204 * t46 - t206 * t45;
t236 = -pkin(2) + t238;
t233 = -t206 * t245 - t242;
t151 = t206 * t158;
t232 = t151 + t242;
t225 = (-qJD(6) + t189) * t137 - t246;
t164 = t213 * t177;
t220 = t164 + t210 * (t169 * t206 - t171 * t204) * t260;
t36 = t259 * t171 + t228;
t216 = t210 * (t204 * t245 - t206 * t251) + t253;
t200 = t210 ^ 2;
t196 = t200 * t215;
t186 = -t197 - t296;
t185 = -t196 - t296;
t180 = t196 + t197;
t179 = (t200 + t201) * qJDD(2);
t178 = -0.2e1 * t192 + t256;
t175 = 0.2e1 * t249 + t257;
t160 = 0.2e1 * t261;
t145 = -t185 * t210 - t265;
t144 = t186 * t213 - t267;
t112 = -t134 + t187;
t111 = t133 - t187;
t107 = t150 * t204 - t151;
t106 = -t134 - t187;
t99 = t134 - t133;
t95 = -t187 - t133;
t86 = -qJD(6) * t137 - t246;
t80 = (-t135 * t212 + t137 * t209) * t189;
t79 = (t135 * t209 + t137 * t212) * t189;
t73 = -t133 - t134;
t68 = t117 + t87;
t63 = (qJD(6) + t189) * t137 + t246;
t62 = t111 * t212 - t288;
t61 = -t112 * t209 + t319;
t60 = -t111 * t209 - t285;
t59 = -t112 * t212 - t320;
t58 = -t137 * t281 + t212 * t87;
t57 = -t137 * t280 - t209 * t87;
t56 = t135 * t280 - t209 * t86;
t55 = -t135 * t281 - t212 * t86;
t52 = -t106 * t209 - t285;
t51 = t106 * t212 - t288;
t49 = t212 * t95 - t320;
t48 = t209 * t95 + t319;
t43 = t160 - t298 + t347;
t42 = t160 - t306 + (t222 + t251) * qJ(5) + t227;
t41 = t209 * t68 + t212 * t225;
t40 = -t209 * t309 - t212 * t63;
t39 = t209 * t225 - t212 * t68;
t38 = t209 * t63 - t212 * t309;
t34 = t204 * t51 + t206 * t52;
t33 = t204 * t52 - t206 * t51;
t32 = -qJ(5) * t114 + t36;
t31 = (-t114 - t197) * pkin(4) + t235;
t29 = t204 * t48 + t206 * t49;
t28 = t204 * t49 - t206 * t48;
t24 = -t210 * t309 + t213 * t34;
t23 = -t210 * t63 + t213 * t29;
t22 = t210 * t75 + t213 * t26;
t21 = t204 * t39 + t206 * t41;
t20 = t204 * t41 - t206 * t39;
t19 = -pkin(9) * t51 + qJ(5) * t309 - t286;
t18 = t204 * t36 + t206 * t35;
t17 = t204 * t35 - t206 * t36;
t16 = -pkin(9) * t48 + qJ(5) * t63 - t289;
t15 = t21 * t213 - t210 * t73;
t14 = -pkin(9) * t52 + t294 * t309 + t289;
t13 = -pkin(9) * t49 + t294 * t63 - t286;
t10 = t18 * t213 + t210 * t47;
t7 = -pkin(9) * t8 - qJ(5) * t37;
t6 = -pkin(9) * t39 + qJ(5) * t73 - t8;
t5 = -pkin(9) * t41 + t294 * t73 - t9;
t4 = -pkin(9) * t9 - t294 * t37;
t3 = t204 * t8 + t206 * t9;
t2 = t204 * t9 - t206 * t8;
t1 = t210 * t37 + t213 * t3;
t27 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t264, 0, 0, 0, 0, 0, 0, (qJDD(2) * t214 - t211 * t215) * t205, (-qJDD(2) * t211 - t214 * t215) * t205, 0, t207 ^ 2 * t264 + (t211 * t131 - t214 * t239 - t224) * t205, 0, 0, 0, 0, 0, 0, t207 * (t183 * t213 + t186 * t210) + (t144 * t211 + t178 * t214) * t205, t207 * (-t182 * t210 + t185 * t213) + (t145 * t211 - t175 * t214) * t205, (t179 * t211 + t180 * t214) * t205, t207 * (t210 * t98 - t213 * t97) + (-t108 * t214 + t211 * t50) * t205, 0, 0, 0, 0, 0, 0, t340, t359, t348, t207 * (t210 * t26 - t213 * t75) + (t211 * t22 - t214 * t25) * t205, 0, 0, 0, 0, 0, 0, t340, t348, -t359, t207 * (t18 * t210 - t213 * t47) + (t10 * t211 - t17 * t214) * t205, 0, 0, 0, 0, 0, 0, t207 * (t210 * t29 + t213 * t63) + (t211 * t23 - t214 * t28) * t205, t207 * (t210 * t34 + t213 * t309) + (t211 * t24 - t214 * t33) * t205, t207 * (t21 * t210 + t213 * t73) + (t15 * t211 - t20 * t214) * t205, t207 * (t210 * t3 - t213 * t37) + (t1 * t211 - t2 * t214) * t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t239, -t131, 0, 0, t237 * t210, t175 * t213 + t178 * t210, t267 + t213 * (-t196 + t296), -t210 * t249 + t164, t210 * (t197 - t296) + t265, 0, pkin(2) * t178 + pkin(8) * t144 - t108 * t213, -pkin(2) * t175 + pkin(8) * t145 + t108 * t210, pkin(2) * t180 + pkin(8) * t179 + t50, -pkin(2) * t108 + pkin(8) * t50, t240, -t326, t323, t216, t343, t220, t210 * (t291 - t336) + t213 * (t45 - t337) + t345, t210 * (t290 + t352) + t213 * (t46 + t354) + t358, -t210 * t25 + t236 * t330 + t353, pkin(8) * t22 + t236 * t25, t240, t323, t326, t220, -t343, t216, t210 * (qJ(5) * t273 - t204 * t43 - t336) + t213 * (pkin(4) * t127 - qJ(5) * t300 - t337 + t36) + t345, t210 * (-qJ(4) * t330 - t204 * t31 + t206 * t32) + t213 * (-pkin(3) * t330 - pkin(4) * t223 - qJ(5) * t121) - pkin(2) * t330 + t353, t210 * (-pkin(4) * t278 + t206 * t42 - t352) + t213 * (-t354 + qJ(5) * t229 + 0.2e1 * t262 + (t304 + t197) * pkin(4) - t301) - t358, t210 * (-qJ(4) * t17 + (pkin(4) * t204 - qJ(5) * t206) * t47) + t213 * (-pkin(3) * t17 + pkin(4) * t36 - qJ(5) * t35) - pkin(2) * t17 + pkin(8) * t10, t210 * (-t204 * t57 + t206 * t58) + t254, t210 * (-t204 * t38 + t206 * t40) + t213 * t99, t210 * (-t204 * t59 + t206 * t61) + t213 * t68, t210 * (-t204 * t55 + t206 * t56) - t254, t210 * (-t204 * t60 + t206 * t62) + t213 * t225, t210 * (-t204 * t79 + t206 * t80) + t213 * t172, t210 * (-qJ(4) * t28 - t13 * t204 + t16 * t206) + t213 * (-pkin(3) * t28 - qJ(5) * t49 + t294 * t48 - t11) - pkin(2) * t28 + pkin(8) * t23, t210 * (-qJ(4) * t33 - t14 * t204 + t19 * t206) + t213 * (-pkin(3) * t33 + pkin(4) * t51 - qJ(5) * t52 - t287 - t209 * (t171 * t140 + t219 + t255) + (-t127 * t209 + t51) * pkin(5)) - pkin(2) * t33 + pkin(8) * t24, t210 * (-qJ(4) * t20 - t204 * t5 + t206 * t6) + t213 * (-pkin(3) * t20 - qJ(5) * t41 + t294 * t39) - pkin(2) * t20 + pkin(8) * t15, t210 * (-qJ(4) * t2 - t204 * t4 + t206 * t7) + t213 * (-pkin(3) * t2 - qJ(5) * t9 + t294 * t8) - pkin(2) * t2 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, t196 - t197, t257, t188, t256, qJDD(3), -t97, -t98, 0, 0, t107, t299, t332, t233, t328, t232, -t290 + t331, -pkin(3) * t222 + t291 + t350, t26 + t344, -pkin(3) * t75 + qJ(4) * t26, t107, t332, -t299, t232, -t328, t233, qJ(5) * t279 + t206 * t43 + t331, t204 * t32 + t206 * t31 + t344, t204 * t42 + t222 * t252 - t350, qJ(4) * t18 + (-qJ(5) * t204 - t252) * t47, t204 * t58 + t206 * t57, t204 * t40 + t206 * t38, t204 * t61 + t206 * t59, t204 * t56 + t206 * t55, t204 * t62 + t206 * t60, t204 * t80 + t206 * t79, pkin(3) * t63 + qJ(4) * t29 + t13 * t206 + t16 * t204, pkin(3) * t309 + qJ(4) * t34 + t14 * t206 + t19 * t204, pkin(3) * t73 + qJ(4) * t21 + t204 * t6 + t206 * t5, -pkin(3) * t37 + qJ(4) * t3 + t204 * t7 + t206 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t316, t222, t114, t75, 0, 0, 0, 0, 0, 0, -t316, t114, -t222, t47, 0, 0, 0, 0, 0, 0, -t63, -t309, -t73, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, -t223, t304, t36, 0, 0, 0, 0, 0, 0, t48, t51, t39, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t99, t68, -t100, t225, t172, -t11, -t12, 0, 0;];
tauJ_reg  = t27;
