% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:29
% EndTime: 2019-03-09 04:33:37
% DurationCPUTime: 4.46s
% Computational Cost: add. (4116->495), mult. (8155->579), div. (0->0), fcn. (5111->10), ass. (0->235)
t180 = cos(qJ(4));
t177 = sin(qJ(4));
t272 = qJD(3) * t177;
t178 = sin(qJ(3));
t274 = qJD(1) * t178;
t106 = t180 * t274 + t272;
t181 = cos(qJ(3));
t273 = qJD(1) * t181;
t343 = qJD(4) - t273;
t292 = t106 * t343;
t259 = t178 * qJDD(1);
t267 = qJD(4) * t178;
t342 = qJD(1) * t267 - qJDD(3);
t43 = t177 * (qJD(3) * (qJD(4) + t273) + t259) + t342 * t180;
t345 = -t43 - t292;
t162 = t181 * qJDD(1);
t261 = qJD(1) * qJD(3);
t100 = t178 * t261 + qJDD(4) - t162;
t344 = t100 * qJ(5) + qJD(5) * t343;
t171 = qJ(1) + pkin(9);
t158 = cos(t171);
t288 = t158 * t178;
t157 = sin(t171);
t319 = g(2) * t157;
t333 = g(1) * t288 + t178 * t319;
t262 = t180 * qJD(3);
t104 = t177 * t274 - t262;
t266 = qJD(4) * t180;
t247 = t178 * t266;
t270 = qJD(3) * t181;
t271 = qJD(3) * t178;
t286 = t177 * t178;
t338 = t104 * t271 - (t177 * t270 + t247) * t343 - t100 * t286 - t181 * t43;
t296 = qJ(5) * t177;
t326 = pkin(4) + pkin(5);
t211 = -t326 * t180 - t296;
t93 = pkin(3) - t211;
t341 = t43 * qJ(6) + t104 * qJD(6);
t340 = 0.2e1 * t344;
t99 = t106 ^ 2;
t339 = -t343 ^ 2 - t99;
t175 = sin(pkin(9));
t151 = pkin(1) * t175 + pkin(7);
t121 = t151 * qJD(1);
t163 = t178 * qJD(2);
t78 = t181 * t121 + t163;
t337 = qJD(5) * t177 + t78;
t294 = t104 * t343;
t244 = t181 * t261;
t42 = -qJD(4) * t262 + (-t244 - t259) * t180 + t342 * t177;
t335 = -t42 - t294;
t62 = qJD(3) * pkin(8) + t78;
t166 = t178 * pkin(8);
t168 = t181 * pkin(3);
t252 = -pkin(2) - t168;
t216 = t252 - t166;
t176 = cos(pkin(9));
t325 = pkin(1) * t176;
t201 = t216 - t325;
t67 = t201 * qJD(1);
t24 = -t177 * t62 + t180 * t67;
t280 = qJD(5) - t24;
t334 = t181 * qJD(2) - t178 * t121;
t332 = g(1) * t158 + t319;
t89 = t100 * pkin(4);
t331 = t89 - qJDD(5);
t227 = qJD(3) * pkin(3) + t334;
t203 = qJ(5) * t106 + t227;
t23 = pkin(4) * t104 - t203;
t322 = pkin(8) * t100;
t330 = -t23 * t343 + t322;
t14 = -t326 * t104 + qJD(6) + t203;
t268 = qJD(4) * t177;
t119 = t151 * qJDD(1);
t260 = qJDD(2) * t178;
t34 = qJDD(3) * pkin(8) + qJD(3) * t334 + t119 * t181 + t260;
t226 = pkin(3) * t178 - pkin(8) * t181;
t111 = t226 * qJD(3);
t44 = qJD(1) * t111 + qJDD(1) * t201;
t238 = t177 * t34 - t180 * t44 + t62 * t266 + t67 * t268;
t285 = t177 * t181;
t69 = t157 * t285 + t158 * t180;
t71 = -t157 * t180 + t158 * t285;
t196 = g(1) * t71 + g(2) * t69 + g(3) * t286 - t238;
t194 = t196 + t331;
t305 = qJ(6) * t42;
t329 = (qJD(6) + t14) * t106 + t194 - t305;
t327 = t104 ^ 2;
t324 = pkin(4) * t177;
t323 = pkin(5) * t100;
t321 = g(1) * t157;
t318 = g(3) * t181;
t317 = pkin(8) - qJ(6);
t245 = t181 * t262;
t284 = t178 * t180;
t316 = -t104 * t245 - t43 * t284;
t25 = t177 * t67 + t180 * t62;
t257 = t326 * t177;
t295 = qJ(5) * t180;
t210 = -t257 + t295;
t315 = t210 * t343 + t337;
t110 = t226 * qJD(1);
t314 = t177 * t110 + t180 * t334;
t283 = t180 * t100;
t74 = t178 * t283;
t313 = t245 * t343 + t74;
t282 = t180 * t181;
t109 = t151 * t282;
t152 = -pkin(2) - t325;
t276 = t166 + t168;
t92 = t152 - t276;
t312 = qJD(4) * t109 + t92 * t268;
t311 = t177 * t111 + t92 * t266;
t263 = qJD(6) * t180;
t30 = qJ(5) * t274 + t314;
t310 = -qJ(6) * t177 * t273 - t317 * t268 - t263 - t30;
t221 = -t295 + t324;
t309 = t221 * t343 - t337;
t125 = t317 * t180;
t59 = t177 * t334;
t236 = -t110 * t180 + t59;
t256 = t326 * t178;
t308 = qJD(4) * t125 - qJD(6) * t177 - (-qJ(6) * t282 - t256) * qJD(1) - t236;
t307 = pkin(8) * qJD(4);
t306 = qJ(5) * t43;
t137 = t343 * qJ(5);
t16 = qJ(6) * t104 + t25;
t10 = t137 + t16;
t304 = t10 * t343;
t18 = t137 + t25;
t303 = t343 * t18;
t302 = t343 * t25;
t301 = t177 * t42;
t300 = t181 * t42;
t298 = t177 * t92 + t109;
t297 = qJ(5) * t104;
t293 = t106 * t104;
t291 = t343 * t180;
t290 = t151 * t177;
t289 = t157 * t181;
t287 = t158 * t181;
t15 = qJ(6) * t106 + t24;
t281 = qJD(5) - t15;
t278 = t333 * t177;
t277 = t333 * t180;
t172 = t178 ^ 2;
t275 = -t181 ^ 2 + t172;
t122 = qJD(1) * t152;
t269 = qJD(4) * t104;
t264 = qJD(5) * t180;
t258 = t177 * t44 + t180 * t34 + t67 * t266;
t255 = t14 * t268;
t254 = t14 * t266;
t253 = g(1) * t287 + g(2) * t289 + g(3) * t178;
t234 = -qJD(3) * t163 + t181 * qJDD(2) - t178 * t119 - t121 * t270;
t35 = -qJDD(3) * pkin(3) - t234;
t6 = t43 * pkin(4) + t42 * qJ(5) - t106 * qJD(5) + t35;
t3 = -pkin(5) * t43 + qJDD(6) - t6;
t251 = t3 - t318;
t250 = t106 * t270;
t249 = t343 * t268;
t248 = t177 * t267;
t242 = -pkin(4) - t290;
t70 = t157 * t282 - t158 * t177;
t241 = -t69 * pkin(4) + qJ(5) * t70;
t72 = t157 * t177 + t158 * t282;
t240 = -t71 * pkin(4) + qJ(5) * t72;
t239 = t106 * t271 + t300;
t108 = t151 * t285;
t237 = t180 * t92 - t108;
t235 = -t42 + t269;
t232 = pkin(4) * t282 + qJ(5) * t285 + t276;
t231 = t106 * t247;
t230 = g(1) * t69 - g(2) * t71;
t229 = g(1) * t70 - g(2) * t72;
t228 = t111 * t180 - t312;
t46 = -qJ(5) * t181 + t298;
t179 = sin(qJ(1));
t182 = cos(qJ(1));
t225 = g(1) * t179 - g(2) * t182;
t9 = -t326 * t343 + t281;
t224 = t10 * t180 + t177 * t9;
t223 = t10 * t177 - t180 * t9;
t222 = pkin(4) * t180 + t296;
t17 = -pkin(4) * t343 + t280;
t220 = t17 * t180 - t177 * t18;
t219 = t17 * t177 + t18 * t180;
t215 = qJ(5) * t271 - qJD(5) * t181 + t311;
t5 = t238 - t331;
t214 = pkin(3) + t222;
t213 = t307 * t343 + t318;
t209 = -t268 * t62 + t258;
t208 = t177 * t100 + t266 * t343;
t207 = -t213 - t6;
t206 = -t248 * t343 + t313;
t205 = -pkin(1) * t179 - t70 * pkin(4) + t158 * pkin(7) - qJ(5) * t69;
t202 = t245 - t248;
t4 = t209 + t344;
t200 = -t227 * t343 - t322;
t199 = t182 * pkin(1) + t158 * pkin(2) + pkin(3) * t287 + t72 * pkin(4) + t157 * pkin(7) + pkin(8) * t288 + qJ(5) * t71;
t184 = qJD(3) ^ 2;
t198 = 0.2e1 * qJDD(1) * t152 + t151 * t184 - t321;
t197 = -t100 + t293;
t195 = 0.2e1 * qJD(3) * t122 - qJDD(3) * t151;
t193 = t42 - t294;
t192 = qJD(4) * t220 + t5 * t177 + t4 * t180;
t190 = t106 * t23 - t194;
t189 = g(1) * t72 + g(2) * t70 + g(3) * t284 - t209;
t188 = t24 * t343 + t189;
t185 = qJD(1) ^ 2;
t167 = t181 * pkin(4);
t142 = qJ(5) * t284;
t132 = g(2) * t288;
t128 = pkin(8) * t287;
t126 = pkin(8) * t289;
t124 = t317 * t177;
t118 = qJDD(3) * t181 - t178 * t184;
t117 = qJDD(3) * t178 + t181 * t184;
t55 = -t142 + (t151 + t324) * t178;
t49 = pkin(4) * t106 + t297;
t48 = t142 + (-t151 - t257) * t178;
t47 = t167 - t237;
t36 = qJ(6) * t286 + t46;
t32 = -pkin(4) * t274 + t236;
t29 = -t326 * t106 - t297;
t26 = pkin(5) * t181 + t108 + t167 + (-qJ(6) * t178 - t92) * t180;
t21 = (qJD(4) * t222 - t264) * t178 + (t151 + t221) * t270;
t13 = t242 * t271 - t228;
t12 = (qJD(4) * t211 + t264) * t178 + (-t151 + t210) * t270;
t11 = (-t178 * t262 - t181 * t268) * t151 + t215;
t8 = (qJ(6) * qJD(4) - qJD(3) * t151) * t284 + (qJD(6) * t178 + (qJ(6) * qJD(3) - qJD(4) * t151) * t181) * t177 + t215;
t7 = (-qJ(6) * t270 - t111) * t180 + (qJ(6) * t268 - t263 + (-pkin(5) + t242) * qJD(3)) * t178 + t312;
t2 = t4 + t341;
t1 = -qJD(6) * t106 + t305 - t323 + t5;
t19 = [qJDD(1), t225, g(1) * t182 + g(2) * t179 (t225 + (t175 ^ 2 + t176 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t172 + 0.2e1 * t178 * t244, 0.2e1 * t162 * t178 - 0.2e1 * t261 * t275, t117, t118, 0, t195 * t178 + (-g(2) * t158 - t198) * t181, t178 * t198 + t181 * t195 + t132, t106 * t202 - t284 * t42, -t231 + (-t250 + (t42 + t269) * t178) * t177 + t316, t206 + t239 (-t272 * t343 + t43) * t181 + (-qJD(3) * t104 - t208) * t178, -t100 * t181 + t271 * t343, t228 * t343 + t237 * t100 + ((t104 * t151 - t177 * t227) * qJD(3) + t238) * t181 + (-t227 * t266 + t151 * t43 + t35 * t177 + (t290 * t343 + t24) * qJD(3)) * t178 + t229, -t311 * t343 - t298 * t100 + ((t151 * t343 - t62) * t268 + (t106 * t151 - t180 * t227) * qJD(3) + t258) * t181 + (t227 * t268 - t151 * t42 + t35 * t180 + (t151 * t291 - t25) * qJD(3)) * t178 - t230, -t100 * t47 + t104 * t21 - t13 * t343 + t43 * t55 + (t23 * t272 + t5) * t181 + (-qJD(3) * t17 + t177 * t6 + t23 * t266) * t178 + t229, -t104 * t11 + t106 * t13 - t42 * t47 - t43 * t46 - t132 + t220 * t270 + (-qJD(4) * t219 - t177 * t4 + t180 * t5 + t321) * t178, t100 * t46 - t106 * t21 + t11 * t343 + t42 * t55 + (-t23 * t262 - t4) * t181 + (qJD(3) * t18 - t180 * t6 + t23 * t268) * t178 + t230, -g(1) * t205 - g(2) * t199 + t18 * t11 + t17 * t13 + t23 * t21 - t216 * t321 + t4 * t46 + t5 * t47 + t6 * t55, -t100 * t26 - t104 * t12 - t343 * t7 - t43 * t48 + (-t14 * t272 + t1) * t181 + (-qJD(3) * t9 - t177 * t3 - t254) * t178 + t229, t100 * t36 + t106 * t12 + t343 * t8 - t42 * t48 + (t14 * t262 - t2) * t181 + (qJD(3) * t10 + t180 * t3 - t255) * t178 + t230, t104 * t8 - t106 * t7 + t26 * t42 + t36 * t43 + t132 + t223 * t270 + (qJD(4) * t224 - t1 * t180 + t177 * t2 - t321) * t178, t2 * t36 + t10 * t8 + t1 * t26 + t9 * t7 + t3 * t48 + t14 * t12 - g(1) * (-pkin(5) * t70 + t205) - g(2) * (pkin(5) * t72 - qJ(6) * t288 + t199) - (-t317 * t178 + t252) * t321; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, t118, -t117, 0, 0, 0, 0, 0, t338, -t202 * t343 + t239 - t74, t338, t231 + (t178 * t235 + t250) * t177 + t316, -t300 + (-qJD(3) * t106 - t249) * t178 + t313, -g(3) + (qJD(3) * t219 - t6) * t181 + (qJD(3) * t23 + t192) * t178, t338, t206 - t239 (t104 * t180 - t106 * t177) * t270 + (t301 + t180 * t43 + (-t104 * t177 - t106 * t180) * qJD(4)) * t178, -g(3) + (qJD(3) * t224 + t3) * t181 + (-qJD(3) * t14 - qJD(4) * t223 + t1 * t177 + t180 * t2) * t178; 0, 0, 0, 0, -t178 * t185 * t181, t275 * t185, t259, t162, qJDD(3), qJD(3) * t78 - t122 * t274 + t234 - t318 + t333, -t260 + (-qJD(1) * t122 - t119) * t181 + t253, t106 * t291 - t301, t177 * t345 + t335 * t180 (-t106 * t178 - t282 * t343) * qJD(1) + t208, -t249 + t283 + (t104 * t178 + t285 * t343) * qJD(1), -t343 * t274, -t24 * t274 - pkin(3) * t43 - t78 * t104 + t59 * t343 + (-t318 - t35 - (t110 + t307) * t343) * t180 + t200 * t177 + t277, pkin(3) * t42 + t314 * t343 + t25 * t274 - t78 * t106 + t200 * t180 + (t213 + t35) * t177 - t278, t309 * t104 + t17 * t274 - t177 * t330 + t207 * t180 - t214 * t43 + t32 * t343 + t277, t104 * t30 - t106 * t32 + (t4 + t343 * t17 + (qJD(4) * t106 - t43) * pkin(8)) * t180 + (pkin(8) * t235 - t303 + t5) * t177 - t253, -t309 * t106 + t207 * t177 - t18 * t274 + t180 * t330 - t214 * t42 - t30 * t343 + t278, t192 * pkin(8) - g(1) * t128 - g(2) * t126 - g(3) * t232 - t17 * t32 - t18 * t30 + t309 * t23 + (t178 * t332 - t6) * t214, -t255 - t100 * t124 - t43 * t93 + t251 * t180 - t308 * t343 - t315 * t104 + (t14 * t285 + t178 * t9) * qJD(1) + t277, t254 + t100 * t125 - t42 * t93 + t251 * t177 + t310 * t343 + t315 * t106 + (-t10 * t178 - t14 * t282) * qJD(1) + t278, t124 * t42 + t125 * t43 - t308 * t106 + t310 * t104 + (-t343 * t9 - t2) * t180 + (-t1 + t304) * t177 + t253, t2 * t125 + t1 * t124 + t3 * t93 - g(1) * (-qJ(6) * t287 + t128) - g(2) * (-qJ(6) * t289 + t126) - g(3) * (pkin(5) * t282 + t232) + t308 * t9 + t315 * t14 + t310 * t10 + (g(3) * qJ(6) + t332 * t93) * t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, t99 - t327, -t193, -t43 + t292, t100, t106 * t227 + t196 + t302, -t104 * t227 + t188, -t104 * t49 - t190 + t302 + t89, pkin(4) * t42 - t306 + (t18 - t25) * t106 + (t17 - t280) * t104, -t104 * t23 + t106 * t49 - t188 + t340, t4 * qJ(5) - t5 * pkin(4) - t23 * t49 - t17 * t25 - g(1) * t240 - g(2) * t241 - g(3) * (-pkin(4) * t286 + t142) + t280 * t18, t104 * t29 + t343 * t16 + (pkin(5) + t326) * t100 + t329, t104 * t14 - t106 * t29 - t15 * t343 - t189 + t340 + t341, t306 - t326 * t42 + (-t10 + t16) * t106 + (-t9 + t281) * t104, t2 * qJ(5) - t1 * t326 - t9 * t16 - t14 * t29 - g(1) * (-pkin(5) * t71 + t240) - g(2) * (-pkin(5) * t69 + t241) - g(3) * (-t177 * t256 + t142) + t281 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, -t193, t339, t190 - t303, t197, t339, t193, -t304 - t323 - t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t345, t335, -t99 - t327, -t10 * t104 + t106 * t9 + t251 + t333;];
tau_reg  = t19;
