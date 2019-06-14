% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPR10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 23:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPR10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:56:10
% EndTime: 2019-05-05 23:56:21
% DurationCPUTime: 4.84s
% Computational Cost: add. (13203->370), mult. (25608->456), div. (0->0), fcn. (16353->8), ass. (0->254)
t209 = cos(qJ(4));
t210 = cos(qJ(3));
t255 = qJD(1) * qJD(3);
t246 = t210 * t255;
t206 = sin(qJ(3));
t253 = t206 * qJDD(1);
t183 = -t246 - t253;
t178 = qJDD(4) - t183;
t205 = sin(qJ(4));
t258 = qJD(1) * t210;
t179 = -t209 * qJD(3) + t205 * t258;
t181 = t205 * qJD(3) + t209 * t258;
t281 = t181 * t179;
t128 = -t281 - t178;
t270 = t205 * t128;
t177 = t181 ^ 2;
t195 = t206 * qJD(1) + qJD(4);
t296 = t195 ^ 2;
t307 = -t177 - t296;
t90 = t209 * t307 + t270;
t347 = pkin(3) * t90;
t346 = pkin(8) * t90;
t263 = t209 * t128;
t92 = -t205 * t307 + t263;
t345 = pkin(8) * t92;
t344 = qJ(2) * t90;
t343 = t206 * t92;
t197 = t210 * qJDD(1);
t247 = t206 * t255;
t184 = t197 - t247;
t241 = t209 * qJDD(3) - t205 * t184;
t132 = t181 * qJD(4) - t241;
t166 = t195 * t181;
t105 = t132 - t166;
t297 = t179 ^ 2;
t158 = t297 - t296;
t342 = t206 * t105 + t210 * (-t209 * t158 - t270);
t341 = -t205 * t158 + t263;
t228 = -t205 * qJDD(3) - t209 * t184;
t133 = -t179 * qJD(4) - t228;
t282 = t179 * t195;
t309 = t133 - t282;
t273 = t205 * t309;
t306 = t177 - t297;
t310 = t132 + t166;
t339 = -t206 * t306 + t210 * (t209 * t310 + t273);
t293 = pkin(7) + pkin(1);
t304 = -t281 + t178;
t119 = t209 * t304;
t303 = -t296 - t297;
t317 = t205 * t303 + t119;
t269 = t205 * t304;
t316 = t209 * t303 - t269;
t329 = t206 * t316 - t210 * t310;
t338 = qJ(2) * t317 - t293 * t329;
t337 = pkin(3) * t317;
t336 = pkin(8) * t316;
t335 = pkin(8) * t317;
t159 = -t177 + t296;
t331 = t209 * t159 + t269;
t308 = t133 + t282;
t328 = t210 * (-t205 * t159 + t119) + t206 * t308;
t305 = t177 + t297;
t327 = pkin(3) * t305;
t204 = sin(qJ(6));
t172 = -qJDD(6) + t178;
t208 = cos(qJ(6));
t138 = -t208 * t179 + t204 * t181;
t140 = t204 * t179 + t208 * t181;
t95 = t140 * t138;
t313 = -t172 - t95;
t326 = t204 * t313;
t322 = t208 * t313;
t320 = t210 * t305;
t318 = t309 * qJ(5);
t315 = -t205 * t310 + t209 * t309;
t192 = -qJD(6) + t195;
t120 = t138 * t192;
t78 = -t138 * qJD(6) + t204 * t132 + t208 * t133;
t314 = t120 + t78;
t213 = qJD(1) ^ 2;
t312 = t293 * t213;
t207 = sin(qJ(1));
t211 = cos(qJ(1));
t244 = t207 * g(1) - t211 * g(2);
t234 = qJDD(2) - t244;
t261 = t213 * qJ(2);
t220 = t234 - t261;
t162 = -t293 * qJDD(1) + t220;
t142 = t206 * g(3) + t210 * t162;
t212 = qJD(3) ^ 2;
t238 = pkin(3) * t206 - pkin(8) * t210;
t224 = t213 * t238;
t114 = qJDD(3) * pkin(3) + t212 * pkin(8) - t210 * t224 + t142;
t311 = -pkin(4) * t166 + t114;
t141 = t179 * pkin(4) - t181 * qJ(5);
t254 = qJD(2) * qJD(1);
t199 = 0.2e1 * t254;
t201 = qJDD(1) * qJ(2);
t236 = t211 * g(1) + t207 * g(2);
t226 = -t201 + t236;
t222 = t199 - t226;
t230 = -t184 + t247;
t231 = -t183 + t246;
t100 = t231 * pkin(3) + t230 * pkin(8) + t222 - t312;
t143 = t210 * g(3) - t206 * t162;
t115 = -t212 * pkin(3) + qJDD(3) * pkin(8) - t206 * t224 - t143;
t72 = t205 * t100 + t209 * t115;
t240 = -t178 * qJ(5) + t179 * t141 - t72;
t302 = -(t307 + t296) * pkin(4) - qJ(5) * t128 - t240;
t294 = pkin(4) + pkin(5);
t71 = -t209 * t100 + t205 * t115;
t44 = -t178 * pkin(4) - qJ(5) * t296 + t181 * t141 + qJDD(5) + t71;
t32 = -t304 * pkin(5) - t308 * pkin(9) + t44;
t155 = -t195 * pkin(5) - t181 * pkin(9);
t257 = qJD(5) * t195;
t189 = 0.2e1 * t257;
t233 = t189 - t240;
t43 = -pkin(4) * t296 + t233;
t35 = -t297 * pkin(5) + t132 * pkin(9) + t195 * t155 + t43;
t17 = t204 * t35 - t208 * t32;
t18 = t204 * t32 + t208 * t35;
t8 = -t208 * t17 + t204 * t18;
t9 = t204 * t17 + t208 * t18;
t301 = qJ(5) * t9 - t294 * t8;
t137 = t140 ^ 2;
t190 = t192 ^ 2;
t113 = -t137 - t190;
t82 = t172 - t95;
t287 = t204 * t82;
t62 = t208 * t113 + t287;
t285 = t208 * t82;
t63 = -t204 * t113 + t285;
t300 = qJ(5) * t63 - t294 * t62 + t18;
t136 = t138 ^ 2;
t89 = -t190 - t136;
t48 = t204 * t89 + t322;
t49 = t208 * t89 - t326;
t299 = qJ(5) * t49 - t294 * t48 + t17;
t243 = -t208 * t132 + t204 * t133;
t55 = (qJD(6) + t192) * t140 + t243;
t58 = -t120 + t78;
t27 = -t204 * t55 - t208 * t58;
t29 = t204 * t58 - t208 * t55;
t298 = qJ(5) * t29 - t294 * t27;
t295 = 0.2e1 * qJD(5);
t292 = pkin(4) * t209;
t290 = t132 * pkin(4);
t215 = t311 + t318;
t36 = -t290 - t132 * pkin(5) - t297 * pkin(9) + (t295 + t155) * t181 + t215;
t288 = t204 * t36;
t286 = t208 * t36;
t284 = qJ(5) * t209;
t283 = qJDD(1) * pkin(1);
t280 = t192 * t204;
t279 = t192 * t208;
t278 = t195 * t205;
t277 = t195 * t209;
t202 = t206 ^ 2;
t276 = t202 * t213;
t203 = t210 ^ 2;
t275 = t203 * t213;
t272 = t205 * t308;
t271 = t205 * t114;
t249 = t206 * t213 * t210;
t266 = t206 * (qJDD(3) + t249);
t264 = t209 * t114;
t262 = t210 * (qJDD(3) - t249);
t259 = t202 + t203;
t252 = t181 * t295;
t251 = t206 * t95;
t250 = t206 * t281;
t248 = t179 * t277;
t245 = -qJ(5) * t205 - pkin(3);
t38 = t205 * t71 + t209 * t72;
t156 = t181 * t278;
t239 = t210 * (t209 * t133 - t156) + t250;
t237 = -pkin(4) * t44 + qJ(5) * t43;
t235 = -t209 * t132 + t179 * t278;
t232 = -pkin(4) * t308 - qJ(5) * t105;
t229 = t205 * t72 - t209 * t71;
t99 = t210 * t142 - t206 * t143;
t227 = qJ(2) + t238;
t221 = (-t179 * t205 - t181 * t209) * t195;
t218 = t210 * (t156 - t248) + t206 * t178;
t217 = t210 * (t205 * t132 + t248) - t250;
t216 = pkin(4) * t304 + qJ(5) * t303 - t44;
t214 = t215 + t252;
t186 = t259 * qJDD(1);
t185 = t197 - 0.2e1 * t247;
t182 = 0.2e1 * t246 + t253;
t167 = -t220 + t283;
t157 = t226 - 0.2e1 * t254 + t312;
t153 = -t266 + t210 * (-t212 - t275);
t152 = t206 * (-t212 - t276) + t262;
t117 = -t137 + t190;
t116 = t136 - t190;
t111 = (qJD(4) + t195) * t179 + t228;
t106 = (-qJD(4) + t195) * t181 + t241;
t102 = t209 * t308;
t101 = t205 * t133 + t181 * t277;
t94 = t137 - t136;
t81 = (t138 * t208 - t140 * t204) * t192;
t80 = (-t138 * t204 - t140 * t208) * t192;
t79 = -t136 - t137;
t77 = -t140 * qJD(6) - t243;
t76 = t209 * t106 + t272;
t75 = -t209 * t105 + t272;
t73 = -t205 * t105 - t102;
t69 = t208 * t116 + t287;
t68 = -t204 * t117 + t322;
t67 = -t204 * t116 + t285;
t66 = -t208 * t117 - t326;
t64 = t210 * t111 + t343;
t60 = t210 * t309 - t343;
t54 = (qJD(6) - t192) * t140 + t243;
t53 = t140 * t280 + t208 * t78;
t52 = t140 * t279 - t204 * t78;
t51 = -t138 * t279 - t204 * t77;
t50 = t138 * t280 - t208 * t77;
t47 = t206 * t76 + t320;
t46 = t206 * t75 + t320;
t45 = t214 - t290;
t42 = qJ(5) * t305 + t44;
t41 = (t305 - t296) * pkin(4) + t233;
t40 = (-t310 - t132) * pkin(4) + t214;
t39 = t252 - t290 + t311 + 0.2e1 * t318;
t34 = t205 * t62 + t209 * t63;
t33 = t205 * t63 - t209 * t62;
t30 = t210 * t114 + t206 * t38;
t28 = -t204 * t314 - t208 * t54;
t26 = t204 * t54 - t208 * t314;
t25 = t205 * t48 + t209 * t49;
t24 = t205 * t49 - t209 * t48;
t23 = t205 * t44 + t209 * t43;
t22 = t205 * t43 - t209 * t44;
t21 = t206 * t34 + t210 * t314;
t20 = t206 * t25 + t210 * t54;
t19 = -pkin(9) * t62 + qJ(5) * t314 + t286;
t16 = t206 * t23 + t210 * t45;
t15 = -pkin(9) * t48 + qJ(5) * t54 + t288;
t14 = t205 * t27 + t209 * t29;
t13 = t205 * t29 - t209 * t27;
t12 = -pkin(9) * t63 + t294 * t314 - t288;
t11 = t206 * t14 + t210 * t79;
t10 = -pkin(9) * t49 + t294 * t54 + t286;
t7 = -pkin(9) * t8 + qJ(5) * t36;
t6 = -pkin(9) * t27 + qJ(5) * t79 - t8;
t5 = -pkin(9) * t29 + t294 * t79 - t9;
t4 = -pkin(9) * t9 + t294 * t36;
t3 = t205 * t8 + t209 * t9;
t2 = t205 * t9 - t209 * t8;
t1 = t206 * t3 + t210 * t36;
t31 = [0, 0, 0, 0, 0, qJDD(1), t244, t236, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t234 - 0.2e1 * t283, t199 + 0.2e1 * t201 - t236, pkin(1) * t167 + qJ(2) * (-t213 * pkin(1) + t222), -t230 * t210, -t210 * t182 - t206 * t185, t262 - t206 * (t212 - t275), t231 * t206, t210 * (-t212 + t276) - t266, 0, qJ(2) * t182 - t293 * t152 - t206 * t157, qJ(2) * t185 - t293 * t153 - t210 * t157, t293 * t186 - t259 * t261 - t99, -qJ(2) * t157 - t293 * t99, t239, -t339, t328, t217, -t342, t218, t210 * (-t271 - t335) - t206 * (t71 - t337) + t338, t210 * (-t264 - t346) - t206 * (t72 - t347) + t344 - t293 * t64, -t210 * t229 + t227 * (t205 * t106 - t102) - t293 * t47, t227 * t229 - t293 * t30, t239, t328, t339, t218, t342, t217, t210 * (-t205 * t40 - t284 * t310 - t335) - t206 * (-t216 - t337) + t338, t210 * (-pkin(8) * t73 - t205 * t41 + t209 * t42) - t206 * (-pkin(3) * t73 - t232) + qJ(2) * t73 - t293 * t46, t210 * (-pkin(4) * t273 + t209 * t39 + t346) - t206 * (-0.2e1 * t257 - t302 + t347) - t344 - t293 * t60, t210 * (-pkin(8) * t22 + (-pkin(4) * t205 + t284) * t45) - t206 * (-pkin(3) * t22 - t237) + qJ(2) * t22 - t293 * t16, t210 * (-t205 * t52 + t209 * t53) - t251, t210 * (-t205 * t26 + t209 * t28) - t206 * t94, t210 * (-t205 * t66 + t209 * t68) - t206 * t58, t210 * (-t205 * t50 + t209 * t51) + t251, t210 * (-t205 * t67 + t209 * t69) + t206 * t55, t210 * (-t205 * t80 + t209 * t81) + t206 * t172, t210 * (-pkin(8) * t24 - t205 * t10 + t209 * t15) - t206 * (-pkin(3) * t24 - t299) + qJ(2) * t24 - t293 * t20, t210 * (-pkin(8) * t33 - t205 * t12 + t209 * t19) - t206 * (-pkin(3) * t33 - t300) + qJ(2) * t33 - t293 * t21, t210 * (-pkin(8) * t13 - t205 * t5 + t209 * t6) - t206 * (-pkin(3) * t13 - t298) + qJ(2) * t13 - t293 * t11, t210 * (-pkin(8) * t2 - t205 * t4 + t209 * t7) - t206 * (-pkin(3) * t2 - t301) + qJ(2) * t2 - t293 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t213, -t167, 0, 0, 0, 0, 0, 0, t152, t153, -t186, t99, 0, 0, 0, 0, 0, 0, t329, t64, t47, t30, 0, 0, 0, 0, 0, 0, t329, t46, t60, t16, 0, 0, 0, 0, 0, 0, t20, t21, t11, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, (-t202 + t203) * t213, t197, -t249, -t253, qJDD(3), t142, t143, 0, 0, t101, t315, t331, t235, -t341, t221, -pkin(3) * t310 + t264 + t336, pkin(3) * t111 - t271 + t345, pkin(8) * t76 + t327 + t38, pkin(3) * t114 + pkin(8) * t38, t101, t331, -t315, t221, t341, t235, t209 * t40 + t245 * t310 + t336, pkin(8) * t75 + t205 * t42 + t209 * t41 + t327, -t345 + t205 * t39 + (pkin(3) + t292) * t309, pkin(8) * t23 + (-t245 + t292) * t45, t205 * t53 + t209 * t52, t205 * t28 + t209 * t26, t205 * t68 + t209 * t66, t205 * t51 + t209 * t50, t205 * t69 + t209 * t67, t205 * t81 + t209 * t80, pkin(3) * t54 + pkin(8) * t25 + t209 * t10 + t205 * t15, pkin(3) * t314 + pkin(8) * t34 + t209 * t12 + t205 * t19, pkin(3) * t79 + pkin(8) * t14 + t205 * t6 + t209 * t5, pkin(3) * t36 + pkin(8) * t3 + t205 * t7 + t209 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t281, t306, t308, -t281, -t105, t178, -t71, -t72, 0, 0, t281, t308, -t306, t178, t105, -t281, t216, t232, t189 + t302, t237, -t95, -t94, -t58, t95, t55, t172, t299, t300, t298, t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t304, t308, t307, t44, 0, 0, 0, 0, 0, 0, t48, t62, t27, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t94, t58, -t95, -t55, -t172, -t17, -t18, 0, 0;];
tauJ_reg  = t31;
