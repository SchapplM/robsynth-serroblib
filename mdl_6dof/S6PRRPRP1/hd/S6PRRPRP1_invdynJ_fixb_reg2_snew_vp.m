% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:40:32
% EndTime: 2019-05-05 03:40:44
% DurationCPUTime: 3.97s
% Computational Cost: add. (16256->393), mult. (35529->546), div. (0->0), fcn. (26031->12), ass. (0->246)
t217 = sin(pkin(11));
t219 = cos(pkin(11));
t227 = cos(qJ(3));
t269 = qJD(2) * t227;
t224 = sin(qJ(3));
t270 = qJD(2) * t224;
t186 = t217 * t269 + t219 * t270;
t223 = sin(qJ(5));
t226 = cos(qJ(5));
t170 = -t226 * qJD(3) + t223 * t186;
t172 = t223 * qJD(3) + t226 * t186;
t140 = t172 * t170;
t261 = qJD(2) * qJD(3);
t250 = t227 * t261;
t260 = t224 * qJDD(2);
t191 = t250 + t260;
t209 = t227 * qJDD(2);
t251 = t224 * t261;
t242 = t209 - t251;
t246 = t217 * t191 - t219 * t242;
t160 = qJDD(5) + t246;
t308 = -t140 + t160;
t316 = pkin(5) * t308;
t163 = t219 * t191 + t217 * t242;
t133 = -t170 * qJD(5) + t223 * qJDD(3) + t226 * t163;
t184 = t217 * t270 - t219 * t269;
t181 = qJD(5) + t184;
t150 = t181 * t170;
t107 = t133 + t150;
t315 = qJ(6) * t107;
t161 = t186 * t184;
t307 = qJDD(3) - t161;
t314 = t217 * t307;
t313 = t219 * t307;
t278 = t223 * t308;
t275 = t226 * t308;
t218 = sin(pkin(6));
t220 = cos(pkin(6));
t286 = sin(pkin(10));
t287 = cos(pkin(10));
t237 = t286 * g(1) - t287 * g(2);
t236 = t220 * t237;
t273 = -g(3) + qJDD(1);
t312 = t218 * t273 + t236;
t169 = t172 ^ 2;
t180 = t181 ^ 2;
t136 = -t169 - t180;
t168 = t170 ^ 2;
t247 = -t226 * qJDD(3) + t223 * t163;
t132 = -t172 * qJD(5) - t247;
t145 = t181 * pkin(5) - t172 * qJ(6);
t154 = t184 * pkin(4) - t186 * pkin(9);
t229 = qJD(3) ^ 2;
t195 = -t287 * g(1) - t286 * g(2);
t225 = sin(qJ(2));
t228 = cos(qJ(2));
t153 = t228 * t195 + t225 * t312;
t230 = qJD(2) ^ 2;
t232 = -t230 * pkin(2) + qJDD(2) * pkin(8) + t153;
t233 = -t218 * t237 + t220 * t273;
t127 = t224 * t233 + t227 * t232;
t196 = qJD(3) * pkin(3) - qJ(4) * t270;
t214 = t227 ^ 2;
t212 = t214 * t230;
t114 = -pkin(3) * t212 + t242 * qJ(4) - qJD(3) * t196 + t127;
t126 = t224 * t232 - t227 * t233;
t205 = t224 * t230 * t227;
t197 = qJDD(3) + t205;
t231 = -t126 + (-t191 + t250) * qJ(4) + t197 * pkin(3);
t70 = -0.2e1 * qJD(4) * t184 + t219 * t114 + t217 * t231;
t58 = -t229 * pkin(4) + qJDD(3) * pkin(9) - t184 * t154 + t70;
t243 = t225 * t195 - t228 * t312;
t146 = -qJDD(2) * pkin(2) - t230 * pkin(8) + t243;
t120 = -t242 * pkin(3) - qJ(4) * t212 + t196 * t270 + qJDD(4) + t146;
t263 = t186 * qJD(3);
t141 = t246 + t263;
t264 = t184 * qJD(3);
t244 = -t163 + t264;
t81 = pkin(4) * t141 + t244 * pkin(9) + t120;
t43 = t223 * t81 + t226 * t58;
t239 = t132 * qJ(6) - 0.2e1 * qJD(6) * t170 - t181 * t145 + t43;
t311 = -t239 + (t136 + t168) * pkin(5);
t309 = t133 - t150;
t104 = (qJD(5) - t181) * t172 + t247;
t182 = t184 ^ 2;
t183 = t186 ^ 2;
t128 = -t180 - t168;
t85 = t223 * t128 + t275;
t306 = pkin(4) * t85;
t118 = t140 + t160;
t279 = t223 * t118;
t90 = t226 * t136 - t279;
t305 = pkin(4) * t90;
t265 = qJD(6) * t172;
t165 = -0.2e1 * t265;
t42 = t223 * t58 - t226 * t81;
t238 = -t315 - t42 + t316;
t26 = t165 + t238;
t304 = pkin(5) * t26;
t73 = -t104 * t223 - t226 * t107;
t303 = pkin(9) * t73;
t302 = pkin(9) * t85;
t301 = pkin(9) * t90;
t300 = pkin(4) * t217;
t299 = pkin(5) * t107;
t125 = -t168 - t169;
t75 = -t104 * t226 + t223 * t107;
t54 = -t219 * t125 + t217 * t75;
t55 = t217 * t125 + t219 * t75;
t32 = -t224 * t54 + t227 * t55;
t298 = -pkin(2) * t73 + pkin(8) * t32;
t103 = (qJD(5) + t181) * t172 + t247;
t86 = t226 * t128 - t278;
t61 = -t219 * t103 + t217 * t86;
t62 = t217 * t103 + t219 * t86;
t35 = -t224 * t61 + t227 * t62;
t297 = -pkin(2) * t85 + pkin(8) * t35;
t276 = t226 * t118;
t91 = -t223 * t136 - t276;
t66 = t217 * t91 - t219 * t309;
t67 = t217 * t309 + t219 * t91;
t37 = -t224 * t66 + t227 * t67;
t296 = -pkin(2) * t90 + pkin(8) * t37;
t295 = qJ(4) * t54;
t294 = qJ(4) * t61;
t293 = qJ(4) * t66;
t292 = t223 * t26;
t248 = t217 * t114 - t219 * t231;
t240 = -qJDD(3) * pkin(4) - t229 * pkin(9) + t248;
t245 = (0.2e1 * qJD(4) + t154) * t186;
t57 = t245 + t240;
t291 = t223 * t57;
t267 = qJD(4) * t186;
t69 = t248 + 0.2e1 * t267;
t44 = t217 * t70 - t219 * t69;
t290 = t224 * t44;
t289 = t226 * t26;
t288 = t226 * t57;
t285 = t181 * t223;
t284 = t181 * t226;
t283 = t217 * t120;
t157 = qJDD(3) + t161;
t282 = t217 * t157;
t281 = t219 * t120;
t280 = t219 * t157;
t277 = t224 * t197;
t198 = qJDD(3) - t205;
t274 = t227 * t198;
t259 = pkin(3) * t54 - pkin(4) * t125 + pkin(9) * t75;
t258 = pkin(3) * t61 - pkin(4) * t103 + pkin(9) * t86;
t257 = pkin(3) * t66 - pkin(4) * t309 + pkin(9) * t91;
t256 = t217 * t140;
t255 = t219 * t140;
t254 = -pkin(4) * t219 - pkin(3);
t253 = -pkin(3) * t85 + qJ(4) * t62;
t252 = -pkin(3) * t90 + qJ(4) * t67;
t45 = t217 * t69 + t219 * t70;
t18 = t223 * t42 + t226 * t43;
t87 = t224 * t126 + t227 * t127;
t192 = t209 - 0.2e1 * t251;
t12 = t217 * t18 - t219 * t57;
t13 = t219 * t18 + t217 * t57;
t3 = -t224 * t12 + t227 * t13;
t17 = t223 * t43 - t226 * t42;
t142 = -t246 + t263;
t235 = t238 + t316;
t234 = -t132 * pkin(5) - t168 * qJ(6) + t172 * t145 + qJDD(6) + t240;
t46 = t245 + t234;
t213 = t224 ^ 2;
t210 = t213 * t230;
t204 = -t212 - t229;
t203 = -t210 - t229;
t194 = t210 + t212;
t193 = (t213 + t214) * qJDD(2);
t190 = 0.2e1 * t250 + t260;
t179 = -0.2e1 * t267;
t177 = -t183 - t229;
t176 = -t183 + t229;
t175 = t182 - t229;
t174 = -t224 * t203 - t274;
t173 = t227 * t204 - t277;
t166 = 0.2e1 * t265;
t155 = -t229 - t182;
t148 = -t169 + t180;
t147 = t168 - t180;
t144 = t163 + t264;
t139 = -t182 - t183;
t137 = t169 - t168;
t135 = -t217 * t177 - t280;
t134 = t219 * t177 - t282;
t122 = t219 * t155 - t314;
t121 = t217 * t155 + t313;
t116 = (-t170 * t226 + t172 * t223) * t181;
t115 = (-t170 * t223 - t172 * t226) * t181;
t113 = t219 * t142 + t217 * t144;
t112 = t217 * t142 - t219 * t144;
t100 = t226 * t133 - t172 * t285;
t99 = t223 * t133 + t172 * t284;
t98 = -t223 * t132 + t170 * t284;
t97 = t226 * t132 + t170 * t285;
t96 = -t224 * t134 + t227 * t135;
t95 = t226 * t147 - t279;
t94 = -t223 * t148 + t275;
t93 = t223 * t147 + t276;
t92 = t226 * t148 + t278;
t82 = -t224 * t121 + t227 * t122;
t78 = -pkin(5) * t309 - qJ(6) * t118;
t77 = -t224 * t112 + t227 * t113;
t76 = -t226 * t103 - t223 * t309;
t74 = -t223 * t103 + t226 * t309;
t63 = t224 * (t219 * t116 + t217 * t160) + t227 * (t217 * t116 - t219 * t160);
t52 = qJ(4) * t55;
t51 = -pkin(4) * t73 + t299;
t50 = t224 * (t219 * t100 + t256) + t227 * (t217 * t100 - t255);
t49 = t224 * (t219 * t98 - t256) + t227 * (t217 * t98 + t255);
t48 = t288 - t301;
t47 = t291 - t302;
t40 = t224 * (-t217 * t104 + t219 * t95) + t227 * (t219 * t104 + t217 * t95);
t39 = t224 * (t217 * t107 + t219 * t94) + t227 * (-t219 * t107 + t217 * t94);
t38 = -qJ(6) * t136 + t46;
t33 = t224 * (t217 * t137 + t219 * t76) + t227 * (-t219 * t137 + t217 * t76);
t30 = t43 - t305;
t29 = -pkin(5) * t103 + qJ(6) * t128 - t186 * t154 + t179 - t234;
t28 = t42 - t306;
t27 = -t168 * pkin(5) + t239;
t25 = t166 - t238 + t315;
t24 = -qJ(6) * t104 + (-t125 - t168) * pkin(5) + t239;
t23 = -t223 * t78 + t226 * t38 - t301;
t22 = -t305 - t311;
t21 = -qJ(6) * t275 - t223 * t29 - t302;
t20 = t166 - t235 - t306;
t19 = t227 * t45 - t290;
t16 = -pkin(5) * t46 + qJ(6) * t27;
t15 = -t17 - t303;
t14 = t220 * (t224 * t67 + t227 * t66) + (t225 * t37 - t228 * t90) * t218;
t11 = t220 * (t224 * t62 + t227 * t61) + (t225 * t35 - t228 * t85) * t218;
t10 = t226 * t27 - t292;
t9 = t223 * t27 + t289;
t8 = t220 * (t224 * t55 + t227 * t54) + (t225 * t32 - t228 * t73) * t218;
t7 = t219 * t10 + t217 * t46;
t6 = t217 * t10 - t219 * t46;
t5 = -t223 * t24 + t226 * t25 - t303;
t4 = -pkin(4) * t9 - t304;
t2 = -pkin(9) * t9 - qJ(6) * t289 - t223 * t16;
t1 = -t224 * t6 + t227 * t7;
t31 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t273, 0, 0, 0, 0, 0, 0, (qJDD(2) * t228 - t225 * t230) * t218, (-qJDD(2) * t225 - t228 * t230) * t218, 0, t220 ^ 2 * t273 + (t225 * t153 - t228 * t243 - t236) * t218, 0, 0, 0, 0, 0, 0, t220 * (t227 * t197 + t224 * t204) + (t225 * t173 + t228 * t192) * t218, t220 * (-t224 * t198 + t227 * t203) + (t225 * t174 - t228 * t190) * t218, (t193 * t225 + t194 * t228) * t218, t220 * (-t227 * t126 + t224 * t127) + (-t228 * t146 + t225 * t87) * t218, 0, 0, 0, 0, 0, 0, t220 * (t227 * t121 + t224 * t122) + (-t228 * t141 + t225 * t82) * t218, t220 * (t227 * t134 + t224 * t135) + (t225 * t96 + t228 * t244) * t218, t220 * (t227 * t112 + t224 * t113) + (-t228 * t139 + t225 * t77) * t218, t220 * (t224 * t45 + t227 * t44) + (-t228 * t120 + t225 * t19) * t218, 0, 0, 0, 0, 0, 0, t11, t14, t8, t220 * (t227 * t12 + t224 * t13) + (-t228 * t17 + t225 * t3) * t218, 0, 0, 0, 0, 0, 0, t11, t14, t8, t220 * (t224 * t7 + t227 * t6) + (t225 * t1 - t228 * t9) * t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t243, -t153, 0, 0, (t191 + t250) * t224, t227 * t190 + t224 * t192, t277 + t227 * (-t210 + t229), t192 * t227, t224 * (t212 - t229) + t274, 0, pkin(2) * t192 + pkin(8) * t173 - t227 * t146, -pkin(2) * t190 + pkin(8) * t174 + t224 * t146, pkin(2) * t194 + pkin(8) * t193 + t87, -pkin(2) * t146 + pkin(8) * t87, t224 * (t219 * t163 - t217 * t263) + t227 * (t217 * t163 + t219 * t263), t224 * (-t219 * t141 + t217 * t244) + t227 * (-t217 * t141 - t219 * t244), t224 * (-t217 * t176 + t313) + t227 * (t219 * t176 + t314), t224 * (t217 * t246 + t219 * t264) + t227 * (t217 * t264 - t219 * t246), t224 * (t219 * t175 - t282) + t227 * (t217 * t175 + t280), (t224 * (-t184 * t219 + t186 * t217) + t227 * (-t184 * t217 - t186 * t219)) * qJD(3), t224 * (-qJ(4) * t121 + t283) + t227 * (-pkin(3) * t141 + qJ(4) * t122 - t281) - pkin(2) * t141 + pkin(8) * t82, t224 * (-qJ(4) * t134 + t281) + t227 * (pkin(3) * t244 + qJ(4) * t135 + t283) + pkin(2) * t244 + pkin(8) * t96, t224 * (-qJ(4) * t112 - t44) + t227 * (-pkin(3) * t139 + qJ(4) * t113 + t45) - pkin(2) * t139 + pkin(8) * t77, -qJ(4) * t290 + t227 * (-pkin(3) * t120 + qJ(4) * t45) - pkin(2) * t120 + pkin(8) * t19, t50, t33, t39, t49, t40, t63, t224 * (-t217 * t28 + t219 * t47 - t294) + t227 * (t217 * t47 + t219 * t28 + t253) + t297, t224 * (-t217 * t30 + t219 * t48 - t293) + t227 * (t217 * t48 + t219 * t30 + t252) + t296, t224 * (t219 * t15 + t73 * t300 - t295) + t227 * (t217 * t15 + t254 * t73 + t52) + t298, (t224 * (-pkin(9) * t219 + t300) + t227 * (-pkin(9) * t217 + t254) - pkin(2)) * t17 + (pkin(8) + qJ(4)) * t3, t50, t33, t39, t49, t40, t63, t224 * (-t217 * t20 + t219 * t21 - t294) + t227 * (t219 * t20 + t217 * t21 + t253) + t297, t224 * (-t217 * t22 + t219 * t23 - t293) + t227 * (t217 * t23 + t219 * t22 + t252) + t296, t224 * (-t217 * t51 + t219 * t5 - t295) + t227 * (-pkin(3) * t73 + t217 * t5 + t219 * t51 + t52) + t298, t224 * (-qJ(4) * t6 + t219 * t2 - t217 * t4) + t227 * (-pkin(3) * t9 + qJ(4) * t7 + t217 * t2 + t219 * t4) - pkin(2) * t9 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, t210 - t212, t260, t205, t209, qJDD(3), -t126, -t127, 0, 0, t161, t183 - t182, t144, -t161, t142, qJDD(3), pkin(3) * t121 + t179 - t248, pkin(3) * t134 - t70, pkin(3) * t112, pkin(3) * t44, t99, t74, t92, t97, t93, t115, t258 - t288, t257 + t291, t18 + t259, pkin(3) * t12 - pkin(4) * t57 + pkin(9) * t18, t99, t74, t92, t97, t93, t115, -qJ(6) * t278 + t226 * t29 + t258, t223 * t38 + t226 * t78 + t257, t223 * t25 + t226 * t24 + t259, pkin(3) * t6 - pkin(4) * t46 + pkin(9) * t10 - qJ(6) * t292 + t226 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, -t244, t139, t120, 0, 0, 0, 0, 0, 0, t85, t90, t73, t17, 0, 0, 0, 0, 0, 0, t85, t90, t73, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, t137, t107, -t140, -t104, t160, -t42, -t43, 0, 0, t140, t137, t107, -t140, -t104, t160, t165 + t235, t311, -t299, t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t309, t125, t46;];
tauJ_reg  = t31;
