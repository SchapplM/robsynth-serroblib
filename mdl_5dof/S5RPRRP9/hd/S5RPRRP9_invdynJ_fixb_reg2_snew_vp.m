% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:13
% EndTime: 2019-12-31 18:49:22
% DurationCPUTime: 5.22s
% Computational Cost: add. (10563->342), mult. (26163->446), div. (0->0), fcn. (19663->8), ass. (0->217)
t202 = sin(pkin(8));
t203 = cos(pkin(8));
t205 = sin(qJ(3));
t208 = cos(qJ(3));
t207 = cos(qJ(4));
t198 = qJDD(3) + qJDD(4);
t181 = (-t202 * t205 + t203 * t208) * qJD(1);
t223 = t202 * t208 + t203 * t205;
t183 = t223 * qJD(1);
t204 = sin(qJ(4));
t154 = -t207 * t181 + t183 * t204;
t156 = t204 * t181 + t183 * t207;
t263 = t156 * t154;
t121 = -t263 - t198;
t257 = t204 * t121;
t153 = t156 ^ 2;
t201 = qJD(3) + qJD(4);
t283 = t201 ^ 2;
t291 = -t153 - t283;
t85 = t207 * t291 + t257;
t253 = t207 * t121;
t87 = -t204 * t291 + t253;
t41 = t205 * t85 - t208 * t87;
t52 = t205 * t87 + t208 * t85;
t333 = qJ(2) * (t202 * t52 + t203 * t41);
t332 = pkin(6) * t41;
t331 = pkin(6) * t52;
t284 = t154 ^ 2;
t142 = t284 - t283;
t92 = -t204 * t142 + t253;
t96 = -t207 * t142 - t257;
t329 = t202 * (t205 * t92 - t208 * t96) - t203 * (t205 * t96 + t208 * t92);
t328 = pkin(3) * t85;
t327 = pkin(7) * t85;
t326 = pkin(7) * t87;
t242 = t203 * qJDD(1);
t243 = t202 * qJDD(1);
t134 = -t205 * t243 + t208 * t242;
t246 = t183 * qJD(3);
t165 = t134 - t246;
t180 = t223 * qJDD(1);
t247 = t181 * qJD(3);
t167 = t180 + t247;
t224 = t204 * t165 + t207 * t167;
t103 = -t154 * qJD(4) + t224;
t260 = t201 * t154;
t293 = t103 - t260;
t148 = t201 * t156;
t232 = -t207 * t165 + t167 * t204;
t220 = qJD(4) * t156 + t232;
t77 = t148 + t220;
t45 = -t204 * t77 + t207 * t293;
t272 = t204 * t293;
t47 = t207 * t77 + t272;
t323 = t202 * (t205 * t45 + t208 * t47) + t203 * (t205 * t47 - t208 * t45);
t289 = -t263 + t198;
t256 = t204 * t289;
t288 = -t283 - t284;
t294 = t207 * t288 - t256;
t111 = t207 * t289;
t295 = t204 * t288 + t111;
t301 = t205 * t294 + t208 * t295;
t320 = pkin(6) * t301;
t302 = -t205 * t295 + t208 * t294;
t319 = pkin(6) * t302;
t316 = qJ(2) * (-t202 * t301 + t203 * t302);
t143 = -t153 + t283;
t303 = t207 * t143 + t256;
t304 = -t204 * t143 + t111;
t315 = t202 * (-t205 * t303 + t208 * t304) + t203 * (t205 * t304 + t208 * t303);
t100 = -t284 - t153;
t314 = pkin(1) * t100;
t313 = pkin(2) * t100;
t312 = pkin(3) * t100;
t311 = pkin(3) * t295;
t310 = pkin(7) * t294;
t309 = pkin(7) * t295;
t261 = t183 * t181;
t296 = qJDD(3) + t261;
t307 = t205 * t296;
t305 = t208 * t296;
t210 = qJD(1) ^ 2;
t206 = sin(qJ(1));
t281 = cos(qJ(1));
t222 = t281 * g(1) + t206 * g(2);
t297 = -t210 * pkin(1) + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t222;
t298 = t293 * qJ(5);
t292 = t103 + t260;
t123 = t153 - t284;
t277 = t203 * pkin(2);
t279 = g(3) * t203;
t290 = -t279 + (-pkin(6) * qJDD(1) + t210 * t277 - t297) * t202;
t236 = g(1) * t206 - t281 * g(2);
t226 = -qJDD(2) + t236;
t250 = t210 * qJ(2);
t264 = qJDD(1) * pkin(1);
t177 = t226 + t250 + t264;
t199 = t202 ^ 2;
t200 = t203 ^ 2;
t249 = t210 * t200;
t287 = qJ(2) * t249 + t199 * t250 - t177 - t264;
t219 = (-t154 * t204 - t156 * t207) * t201;
t259 = t201 * t204;
t141 = t156 * t259;
t258 = t201 * t207;
t238 = t154 * t258;
t227 = t141 - t238;
t286 = t202 * (-t205 * t219 + t208 * t227) + t203 * (t205 * t227 + t208 * t219);
t221 = t204 * t220 + t238;
t229 = t154 * t259 - t207 * t220;
t285 = t202 * (-t205 * t229 + t208 * t221) + t203 * (t205 * t221 + t208 * t229);
t178 = t181 ^ 2;
t179 = t183 ^ 2;
t282 = 2 * qJD(5);
t280 = pkin(4) * t207;
t278 = t220 * pkin(4);
t122 = pkin(4) * t154 - qJ(5) * t156;
t230 = -g(3) * t202 + t297 * t203;
t149 = -pkin(2) * t249 + pkin(6) * t242 + t230;
t105 = t149 * t205 - t208 * t290;
t215 = (-t167 + t247) * pkin(7) + t296 * pkin(3) - t105;
t106 = t208 * t149 + t290 * t205;
t225 = qJD(3) * pkin(3) - pkin(7) * t183;
t66 = -t178 * pkin(3) + t165 * pkin(7) - qJD(3) * t225 + t106;
t38 = t204 * t215 + t207 * t66;
t228 = t198 * qJ(5) - t154 * t122 + t201 * t282 + t38;
t30 = -pkin(4) * t283 + t228;
t37 = t204 * t66 - t207 * t215;
t32 = -t198 * pkin(4) - qJ(5) * t283 + t122 * t156 + qJDD(5) + t37;
t276 = -pkin(4) * t32 + qJ(5) * t30;
t78 = -t148 + t220;
t275 = -pkin(4) * t292 - qJ(5) * t78;
t57 = -t208 * t105 + t205 * t106;
t274 = t202 * t57;
t271 = t204 * t292;
t248 = t199 + t200;
t160 = (pkin(1) + t277) * qJDD(1) + (t248 * pkin(6) + qJ(2)) * t210 + t226;
t97 = t165 * pkin(3) + t178 * pkin(7) - t183 * t225 + t160;
t270 = t204 * t97;
t18 = t204 * t38 - t207 * t37;
t269 = t205 * t18;
t267 = t207 * t97;
t266 = t208 * t18;
t265 = qJ(5) * t207;
t262 = t160 * t205;
t162 = qJDD(3) - t261;
t254 = t205 * t162;
t252 = t208 * t160;
t251 = t208 * t162;
t245 = qJD(4) + t201;
t237 = -qJ(5) * t204 - pkin(3);
t19 = t204 * t37 + t207 * t38;
t58 = t105 * t205 + t208 * t106;
t231 = t202 * (t297 * t202 + t279) + t203 * t230;
t217 = -pkin(4) * t291 - qJ(5) * t121 + t30;
t216 = pkin(4) * t289 + qJ(5) * t288 - t32;
t214 = -pkin(4) * t148 + t156 * t282 + t97;
t213 = t214 + t298;
t209 = qJD(3) ^ 2;
t195 = t200 * qJDD(1);
t194 = t199 * qJDD(1);
t185 = t248 * t210;
t173 = -t179 - t209;
t172 = -t179 + t209;
t171 = t178 - t209;
t166 = t180 + 0.2e1 * t247;
t164 = -t134 + 0.2e1 * t246;
t159 = -t209 - t178;
t136 = -t178 - t179;
t133 = -t205 * t173 - t251;
t132 = t208 * t173 - t254;
t131 = t208 * t134 + t205 * t180;
t130 = t205 * t134 - t208 * t180;
t129 = t208 * t159 - t307;
t128 = t205 * t159 + t305;
t80 = -t245 * t154 + t224;
t79 = (-qJD(4) + t201) * t156 - t232;
t76 = t245 * t156 + t232;
t73 = t207 * t292;
t72 = t207 * t103 - t141;
t71 = t204 * t103 + t156 * t258;
t56 = -t267 - t327;
t51 = -t270 - t309;
t50 = t207 * t79 + t271;
t48 = -t207 * t78 + t271;
t46 = t204 * t79 - t73;
t44 = -t204 * t78 - t73;
t35 = -pkin(3) * t80 - t270 + t326;
t34 = -pkin(3) * t77 + t267 + t310;
t33 = t213 - t278;
t28 = -qJ(5) * t100 + t32;
t27 = (-t100 - t283) * pkin(4) + t228;
t26 = (-t220 - t76) * pkin(4) + t213;
t25 = t214 - t278 + 0.2e1 * t298;
t24 = -t205 * t46 + t208 * t50;
t23 = -t205 * t44 + t208 * t48;
t22 = t205 * t50 + t208 * t46;
t21 = t205 * t48 + t208 * t44;
t20 = t202 * (-t205 * t71 + t208 * t72) + t203 * (t205 * t72 + t208 * t71);
t17 = pkin(3) * t97 + pkin(7) * t19;
t16 = -t204 * t26 - t76 * t265 - t309;
t15 = -pkin(4) * t272 + t207 * t25 + t327;
t14 = t204 * t32 + t207 * t30;
t13 = t204 * t30 - t207 * t32;
t12 = -pkin(7) * t46 - t18;
t11 = t207 * t26 + t237 * t76 + t310;
t10 = -t326 + t204 * t25 + (pkin(3) + t280) * t293;
t9 = pkin(7) * t50 + t19 - t312;
t8 = t208 * t19 - t269;
t7 = t205 * t19 + t266;
t6 = -pkin(7) * t44 - t204 * t27 + t207 * t28;
t5 = pkin(7) * t48 + t204 * t28 + t207 * t27 - t312;
t4 = -pkin(7) * t13 + (-pkin(4) * t204 + t265) * t33;
t3 = -t205 * t13 + t208 * t14;
t2 = t208 * t13 + t205 * t14;
t1 = pkin(7) * t14 + (-t237 + t280) * t33;
t29 = [0, 0, 0, 0, 0, qJDD(1), t236, t222, 0, 0, t194, 0.2e1 * t202 * t242, 0, t195, 0, 0, -t287 * t203, t287 * t202, pkin(1) * t185 + qJ(2) * (t195 + t194) + t231, pkin(1) * t177 + qJ(2) * t231, t202 * (t208 * t167 - t205 * t246) + t203 * (t205 * t167 + t208 * t246), t202 * (-t208 * t164 - t205 * t166) + t203 * (-t205 * t164 + t208 * t166), t202 * (-t205 * t172 + t305) + t203 * (t208 * t172 + t307), t202 * (-t205 * t165 - t208 * t247) + t203 * (t208 * t165 - t205 * t247), t202 * (t208 * t171 - t254) + t203 * (t205 * t171 + t251), (t202 * (t181 * t208 + t183 * t205) + t203 * (t181 * t205 - t183 * t208)) * qJD(3), t202 * (-pkin(6) * t128 - t262) + t203 * (-pkin(2) * t164 + pkin(6) * t129 + t252) - pkin(1) * t164 + qJ(2) * (-t128 * t202 + t129 * t203), t202 * (-pkin(6) * t132 - t252) + t203 * (-pkin(2) * t166 + pkin(6) * t133 - t262) - pkin(1) * t166 + qJ(2) * (-t132 * t202 + t133 * t203), t202 * (-pkin(6) * t130 - t57) + t203 * (-pkin(2) * t136 + pkin(6) * t131 + t58) - pkin(1) * t136 + qJ(2) * (-t130 * t202 + t131 * t203), -pkin(6) * t274 + t203 * (pkin(2) * t160 + pkin(6) * t58) + pkin(1) * t160 + qJ(2) * (t203 * t58 - t274), t20, -t323, t315, t285, t329, t286, t202 * (-t205 * t34 + t208 * t51 - t320) + t203 * (-pkin(2) * t77 + t205 * t51 + t208 * t34 + t319) - pkin(1) * t77 + t316, t202 * (-t205 * t35 + t208 * t56 - t331) + t203 * (-pkin(2) * t80 + t205 * t56 + t208 * t35 - t332) - pkin(1) * t80 - t333, t202 * (-pkin(6) * t22 + t208 * t12 - t205 * t9) + t203 * (pkin(6) * t24 + t205 * t12 + t208 * t9 - t313) - t314 + qJ(2) * (-t202 * t22 + t203 * t24), t202 * (-pkin(6) * t7 - pkin(7) * t266 - t205 * t17) + t203 * (pkin(2) * t97 + pkin(6) * t8 - pkin(7) * t269 + t208 * t17) + pkin(1) * t97 + qJ(2) * (-t202 * t7 + t203 * t8), t20, t315, t323, t286, -t329, t285, t202 * (-t205 * t11 + t208 * t16 - t320) + t203 * (-pkin(2) * t76 + t208 * t11 + t205 * t16 + t319) - pkin(1) * t76 + t316, t202 * (-pkin(6) * t21 - t205 * t5 + t208 * t6) + t203 * (pkin(6) * t23 + t205 * t6 + t208 * t5 - t313) - t314 + qJ(2) * (-t202 * t21 + t203 * t23), t202 * (-t205 * t10 + t208 * t15 + t331) + t203 * (pkin(2) * t293 + t208 * t10 + t205 * t15 + t332) + pkin(1) * t293 + t333, t202 * (-pkin(6) * t2 - t205 * t1 + t208 * t4) + t203 * (pkin(2) * t33 + pkin(6) * t3 + t208 * t1 + t205 * t4) + pkin(1) * t33 + qJ(2) * (-t202 * t2 + t203 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t242, t243, -t185, -t177, 0, 0, 0, 0, 0, 0, t164, t166, t136, -t160, 0, 0, 0, 0, 0, 0, t77, t80, t100, -t97, 0, 0, 0, 0, 0, 0, t76, t100, -t293, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261, t179 - t178, t180, t261, t134, qJDD(3), -t105, -t106, 0, 0, t263, t123, t292, -t263, -t78, t198, -t37 + t311, -t38 + t328, pkin(3) * t46, pkin(3) * t18, t263, t292, -t123, t198, t78, -t263, t216 + t311, pkin(3) * t44 + t275, t217 - t328, pkin(3) * t13 + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, t123, t292, -t263, -t78, t198, -t37, -t38, 0, 0, t263, t292, -t123, t198, t78, -t263, t216, t275, t217, t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t289, t292, t291, t32;];
tauJ_reg = t29;
