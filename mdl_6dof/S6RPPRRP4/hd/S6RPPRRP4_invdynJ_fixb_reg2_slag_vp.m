% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:30
% EndTime: 2019-03-09 02:06:36
% DurationCPUTime: 4.20s
% Computational Cost: add. (5131->486), mult. (9195->580), div. (0->0), fcn. (5587->8), ass. (0->244)
t169 = sin(qJ(4));
t165 = sin(pkin(9));
t166 = cos(pkin(9));
t328 = sin(qJ(1));
t329 = cos(qJ(1));
t110 = -t328 * t165 - t329 * t166;
t111 = t329 * t165 - t328 * t166;
t228 = g(1) * t110 + g(2) * t111;
t171 = cos(qJ(4));
t322 = g(3) * t171;
t193 = -t169 * t228 + t322;
t168 = sin(qJ(5));
t275 = t168 * qJD(4);
t170 = cos(qJ(5));
t286 = qJD(1) * t170;
t113 = t169 * t286 - t275;
t280 = qJD(4) * t171;
t272 = qJD(1) * qJD(2);
t139 = t166 * t272;
t172 = -pkin(1) - pkin(2);
t129 = t172 * qJDD(1) + qJDD(2);
t268 = qJDD(1) * t166;
t294 = qJ(2) * t268 + t165 * t129;
t73 = t139 + t294;
t343 = qJDD(1) * pkin(7) - qJD(3) * qJD(4) - t73;
t130 = t172 * qJD(1) + qJD(2);
t288 = qJ(2) * qJD(1);
t90 = t165 * t130 + t166 * t288;
t78 = -qJD(1) * pkin(7) + t90;
t22 = t171 * qJDD(3) + t343 * t169 - t78 * t280;
t18 = -qJDD(4) * pkin(4) - t22;
t281 = qJD(4) * t170;
t287 = qJD(1) * t169;
t112 = t168 * t287 + t281;
t279 = qJD(5) * t112;
t271 = qJD(1) * qJD(4);
t243 = t171 * t271;
t265 = t169 * qJDD(1);
t340 = t243 + t265;
t48 = t168 * qJDD(4) - t340 * t170 + t279;
t277 = qJD(5) * t170;
t246 = t169 * t277;
t250 = t171 * t275;
t199 = t246 + t250;
t49 = qJD(1) * t199 - qJD(5) * t275 + t170 * qJDD(4) + t168 * t265;
t5 = -pkin(5) * t49 - qJ(6) * t48 + qJD(6) * t113 + t18;
t345 = -t5 + t193;
t285 = qJD(1) * t171;
t134 = qJD(5) + t285;
t303 = t113 * t134;
t304 = t112 * t134;
t344 = (t49 + t303) * t168 + (t48 + t304) * t170;
t298 = t168 * t171;
t104 = t165 * t298 + t166 * t170;
t249 = t169 * t281;
t296 = t170 * t171;
t318 = -qJD(5) * t104 - t165 * t249 - (t165 * t168 + t166 * t296) * qJD(1);
t254 = t166 * t287;
t339 = -t165 * t280 + t254;
t338 = -t18 + t193;
t337 = t49 - t303;
t153 = t171 * qJDD(1);
t244 = t169 * t271;
t335 = -t244 + t153;
t75 = t169 * t78;
t63 = qJD(3) * t171 - t75;
t57 = -qJD(4) * pkin(4) - t63;
t24 = -pkin(5) * t112 + qJ(6) * t113 + t57;
t109 = -qJDD(5) - t335;
t326 = pkin(8) * t109;
t334 = t134 * t24 + t326;
t278 = qJD(5) * t168;
t206 = t170 * t109 + t134 * t278;
t333 = qJD(1) * (-t112 * t169 + t134 * t298) + t206;
t258 = -t169 * qJDD(3) + t343 * t171;
t282 = qJD(4) * t169;
t21 = -t282 * t78 - t258;
t274 = t169 * qJD(3);
t64 = t171 * t78 + t274;
t183 = -(t169 * t64 + t171 * t63) * qJD(4) - t22 * t169 + t21 * t171;
t230 = -pkin(4) * t169 + pkin(8) * t171;
t100 = t165 * qJD(2) + qJD(4) * t230;
t120 = t166 * qJ(2) + t165 * t172;
t116 = -pkin(7) + t120;
t284 = qJD(2) * t166;
t253 = t171 * t284;
t276 = qJD(5) * t171;
t119 = -t165 * qJ(2) + t166 * t172;
t115 = pkin(3) - t119;
t321 = t169 * pkin(8);
t231 = pkin(4) * t171 + t321;
t79 = t115 + t231;
t10 = -t168 * (qJD(5) * t79 - t116 * t282 + t253) - t170 * (t116 * t276 - t100);
t252 = t134 * t281;
t283 = qJD(4) * t113;
t332 = t169 * (t206 - t283) - t171 * (t48 + t252);
t257 = t165 * t296;
t301 = t166 * t168;
t105 = t257 - t301;
t302 = t165 * t169;
t331 = t105 * t109 + t113 * t339 - t134 * t318 + t302 * t48;
t330 = t113 ^ 2;
t327 = pkin(5) * t109;
t325 = pkin(8) * t113;
t159 = g(3) * t169;
t320 = t5 * t170;
t248 = t170 * t280;
t297 = t169 * t170;
t41 = t49 * t297;
t319 = t112 * t248 + t41;
t317 = qJD(5) * t257 + t165 * t286 - t166 * t278 - t275 * t302 - t285 * t301;
t39 = t116 * t296 + t168 * t79;
t223 = pkin(5) * t168 - qJ(6) * t170;
t316 = qJD(5) * t223 - t168 * qJD(6) - t274 - (-qJD(1) * t223 + t78) * t171;
t58 = qJD(4) * pkin(8) + t64;
t89 = t166 * t130 - t165 * t288;
t77 = qJD(1) * pkin(3) - t89;
t59 = qJD(1) * t231 + t77;
t20 = t168 * t59 + t170 * t58;
t13 = qJ(6) * t134 + t20;
t315 = t13 * t134;
t314 = t134 * t20;
t313 = t170 * t79;
t312 = t18 * t170;
t311 = t49 * t170;
t118 = t230 * qJD(1);
t34 = t168 * t118 + t170 * t63;
t310 = pkin(1) * qJDD(1);
t309 = qJD(1) * t77;
t308 = t109 * qJ(6);
t307 = t110 * t171;
t306 = t111 * t171;
t305 = t112 * t113;
t174 = qJD(1) ^ 2;
t300 = t166 * t174;
t299 = t168 * t169;
t19 = -t168 * t58 + t170 * t59;
t295 = qJD(6) - t19;
t293 = t329 * pkin(1) + t328 * qJ(2);
t292 = g(1) * t328 - g(2) * t329;
t162 = t169 ^ 2;
t163 = t171 ^ 2;
t291 = t162 - t163;
t290 = t162 + t163;
t173 = qJD(4) ^ 2;
t289 = t173 + t174;
t273 = qJ(2) * qJDD(1);
t269 = qJDD(1) * t165;
t267 = qJDD(4) * t169;
t266 = qJDD(4) * t171;
t84 = t113 * t246;
t264 = t113 * t250 - t48 * t299 + t84;
t263 = t168 * t100 + t170 * t253 + t79 * t277;
t262 = g(3) * t296 - t228 * t297;
t261 = pkin(8) * qJD(5) * t134;
t260 = pkin(8) * t277;
t259 = 0.2e1 * t272;
t256 = t169 * t174 * t171;
t255 = t329 * pkin(2) + t293;
t247 = t169 * t278;
t245 = t112 ^ 2 - t330;
t137 = t165 * t272;
t242 = t112 * t282 - t49 * t171;
t17 = qJDD(4) * pkin(8) + t21;
t239 = -qJ(2) * t269 + t166 * t129;
t72 = -t137 + t239;
t66 = qJDD(1) * pkin(3) - t72;
t31 = t335 * pkin(4) + t340 * pkin(8) + t66;
t241 = t168 * t17 - t170 * t31 + t58 * t277 + t59 * t278;
t240 = t48 - t279;
t238 = qJDD(1) * t290;
t237 = qJDD(2) - t310;
t236 = t169 * t243;
t50 = -t110 * t170 + t111 * t298;
t54 = -t110 * t298 - t111 * t170;
t235 = g(1) * t50 + g(2) * t54;
t51 = t110 * t168 + t111 * t296;
t55 = -t110 * t296 + t111 * t168;
t234 = -g(1) * t51 - g(2) * t55;
t233 = (g(1) * t307 + g(2) * t306) * pkin(8);
t232 = -t328 * pkin(1) + t329 * qJ(2);
t229 = -g(1) * t111 + g(2) * t110;
t225 = t240 * pkin(8);
t224 = -pkin(5) * t170 - qJ(6) * t168;
t12 = -pkin(5) * t134 + t295;
t222 = t12 * t170 - t13 * t168;
t221 = t12 * t168 + t13 * t170;
t220 = t165 * t89 - t166 * t90;
t219 = t168 * t20 + t170 * t19;
t218 = -t168 * t19 + t170 * t20;
t215 = t169 * t63 - t171 * t64;
t33 = t118 * t170 - t168 * t63;
t212 = pkin(4) - t224;
t211 = t116 - t223;
t210 = t169 * t18 + t280 * t57;
t3 = t168 * t31 + t170 * t17 + t59 * t277 - t278 * t58;
t207 = -t168 * t109 + t134 * t277;
t205 = -t110 * pkin(3) + pkin(7) * t111 + t255;
t204 = -t112 * t247 + t319;
t203 = g(1) * t329 + g(2) * t328;
t202 = -t228 + t309;
t200 = -t247 + t248;
t197 = -t328 * pkin(2) + t232;
t196 = -pkin(4) * t307 - t110 * t321 + t205;
t195 = t134 * t57 + t326;
t194 = pkin(8) * t311 + t171 * t228 + t159;
t191 = -t168 * t304 + t311;
t189 = t104 * t48 + t105 * t49 + t318 * t112 - t113 * t317;
t188 = g(1) * t54 - g(2) * t50 - g(3) * t299 - t241;
t187 = t111 * pkin(3) + t110 * pkin(7) + t197;
t186 = -qJDD(4) * t116 + (-qJD(1) * t115 - t284 - t77) * qJD(4);
t1 = qJD(6) * t134 + t3 - t308;
t2 = qJDD(6) + t241 + t327;
t185 = qJD(5) * t222 + t1 * t170 + t2 * t168;
t184 = -qJD(5) * t219 + t168 * t241 + t3 * t170;
t182 = t104 * t109 + t339 * t112 - t134 * t317 - t49 * t302;
t181 = pkin(4) * t306 + t111 * t321 + t187;
t180 = t112 * t199 + t299 * t49;
t179 = t109 * t299 - t134 * t199;
t178 = -g(1) * t55 + g(2) * t51 + g(3) * t297 + t3;
t177 = qJDD(1) * t115 - t116 * t173 + t137 + t229 + t66;
t176 = -t113 * t24 + qJDD(6) - t188;
t175 = t179 - t242;
t125 = -t169 * t173 + t266;
t124 = -t171 * t173 - t267;
t122 = t134 * t287;
t70 = -t109 * t171 - t134 * t282;
t65 = -pkin(5) * t113 - qJ(6) * t112;
t62 = t211 * t169;
t38 = -t116 * t298 + t313;
t36 = -t313 + (t116 * t168 - pkin(5)) * t171;
t35 = qJ(6) * t171 + t39;
t30 = pkin(5) * t287 - t33;
t29 = -qJ(6) * t287 + t34;
t28 = t48 - t304;
t25 = (-t113 * t169 + t134 * t296) * qJD(1) + t207;
t23 = t211 * t280 + (qJD(5) * t224 + qJD(6) * t170 + t284) * t169;
t14 = t48 * t168 - t170 * t303;
t11 = t113 * t200 - t297 * t48;
t9 = (-t168 * t276 - t249) * t116 + t263;
t8 = pkin(5) * t282 - t10;
t7 = (t48 - t252) * t171 + (t206 + t283) * t169;
t6 = (-t116 * t278 + qJD(6)) * t171 + (-t116 * t170 - qJ(6)) * t282 + t263;
t4 = [0, 0, 0, 0, 0, qJDD(1), t292, t203, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + t292 + 0.2e1 * t310, 0, -t203 + t259 + 0.2e1 * t273, -t237 * pkin(1) - g(1) * t232 - g(2) * t293 + (t259 + t273) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), -qJDD(1) * t119 + 0.2e1 * t137 + t229 - t239, qJDD(1) * t120 + 0.2e1 * t139 + t228 + t294, 0, -g(1) * t197 - g(2) * t255 - qJD(2) * t220 + t72 * t119 + t73 * t120, qJDD(1) * t162 + 0.2e1 * t236, 0.2e1 * t153 * t169 - 0.2e1 * t271 * t291, t124, qJDD(1) * t163 - 0.2e1 * t236, -t125, 0, t169 * t186 + t171 * t177, -t169 * t177 + t171 * t186, -t116 * t238 - t139 * t290 - t183 - t228, t66 * t115 - g(1) * t187 - g(2) * t205 + (t77 * t165 - t166 * t215) * qJD(2) + t183 * t116, t11, -t204 - t264, t7, t180 (t134 * t275 + t49) * t171 + (-qJD(4) * t112 + t207) * t169, t70, t10 * t134 - t38 * t109 + (-t241 + (-t112 * t116 - t168 * t57) * qJD(4)) * t171 + (-qJD(4) * t19 - t112 * t284 - t116 * t49 - t18 * t168 - t277 * t57) * t169 + t234, t39 * t109 - t9 * t134 + (-t3 + (-t113 * t116 - t170 * t57) * qJD(4)) * t171 + (qJD(4) * t20 - t113 * t284 + t116 * t48 + t278 * t57 - t312) * t169 + t235, t10 * t113 + t9 * t112 - t38 * t48 + t39 * t49 + t219 * t280 + (qJD(5) * t218 + t168 * t3 - t170 * t241 + t229) * t169, t169 * t284 * t57 - g(1) * t181 - g(2) * t196 + t19 * t10 + t116 * t210 + t20 * t9 - t241 * t38 + t3 * t39, t11, t7, t112 * t200 + t264 + t41, t70, t179 + t242, t180, t36 * t109 - t23 * t112 - t8 * t134 - t62 * t49 + (-t24 * t275 - t2) * t171 + (qJD(4) * t12 - t5 * t168 - t24 * t277) * t169 + t234, t6 * t112 - t8 * t113 + t35 * t49 + t36 * t48 - t222 * t280 + (qJD(5) * t221 + t1 * t168 - t170 * t2 + t229) * t169, -t35 * t109 + t23 * t113 + t6 * t134 - t62 * t48 + (t24 * t281 + t1) * t171 + (-qJD(4) * t13 - t24 * t278 + t320) * t169 - t235, t1 * t35 + t13 * t6 + t5 * t62 + t24 * t23 + t2 * t36 + t12 * t8 - g(1) * (t51 * pkin(5) + t50 * qJ(6) + t181) - g(2) * (pkin(5) * t55 + qJ(6) * t54 + t196); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t174, -qJ(2) * t174 + t237 - t292, 0, 0, 0, 0, 0, 0, -t165 * t174 - t268, t269 - t300, 0, qJD(1) * t220 + t73 * t165 + t72 * t166 - t292, 0, 0, 0, 0, 0, 0 (0.2e1 * t244 - t153) * t166 + (-t171 * t289 - t267) * t165 (0.2e1 * t243 + t265) * t166 + (t169 * t289 - t266) * t165, -t165 * t238 + t290 * t300 (qJD(1) * t215 - t66) * t166 + (t183 - t309) * t165 - t292, 0, 0, 0, 0, 0, 0, t182, t331, t189, t104 * t241 + t3 * t105 + t165 * t210 - t19 * t317 + t20 * t318 - t254 * t57 - t292, 0, 0, 0, 0, 0, 0, t182, t189, -t331, -t24 * t254 + t1 * t105 + t2 * t104 + (t169 * t5 + t24 * t280) * t165 + t318 * t13 + t317 * t12 - t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + g(3), 0, 0, 0, 0, 0, 0, t125, t124, 0, -qJD(4) * t215 + t21 * t169 + t22 * t171 + g(3), 0, 0, 0, 0, 0, 0, t175, t332, -t84 + (-t113 * t280 + t169 * t240) * t168 + t319, g(3) + (qJD(4) * t218 - t18) * t171 + (qJD(4) * t57 + t184) * t169, 0, 0, 0, 0, 0, 0, t175, t204 - t264, -t332, g(3) + (qJD(4) * t221 - t5) * t171 + (qJD(4) * t24 + t185) * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256, t291 * t174, -t265, t256, -t153, qJDD(4), t64 * qJD(4) + t169 * t202 + t22 + t322, -t159 + (t63 + t75) * qJD(4) + t202 * t171 + t258, 0, 0, t14, t344, t25, t191, -t333, t122, t19 * t287 + pkin(4) * t49 + t64 * t112 - t312 + (-t33 - t260) * t134 + t195 * t168 + t262, -t20 * t287 - pkin(4) * t48 + t64 * t113 + t34 * t134 + t195 * t170 + (t261 - t338) * t168, -t34 * t112 - t33 * t113 + (-t19 * t285 + t3 + (-t19 - t325) * qJD(5)) * t170 + (t225 + t241 - t314) * t168 + t194, -t19 * t33 - t20 * t34 - t57 * t64 + t338 * pkin(4) + (t184 + t159) * pkin(8) + t233, t14, t25, -t344, t122, t333, t191, -t12 * t287 + t212 * t49 - t320 + (t30 - t260) * t134 - t316 * t112 + t334 * t168 + t262, -t29 * t112 + t30 * t113 + (t12 * t285 + t1 + (t12 - t325) * qJD(5)) * t170 + (t2 + t225 - t315) * t168 + t194, t13 * t287 + t212 * t48 - t29 * t134 + t316 * t113 - t334 * t170 + (-t261 + t345) * t168, -t12 * t30 - t13 * t29 + t316 * t24 + (t185 + t159) * pkin(8) + t233 + t345 * t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t305, -t245, t28, -t305, t337, -t109, t113 * t57 + t188 + t314, -t112 * t57 + t134 * t19 - t178, 0, 0, t305, t28, t245, -t109, -t337, -t305, t112 * t65 - t176 + t314 - 0.2e1 * t327, -pkin(5) * t48 + t49 * qJ(6) + (-t13 + t20) * t113 + (-t12 + t295) * t112, -0.2e1 * t308 + t24 * t112 - t65 * t113 + (0.2e1 * qJD(6) - t19) * t134 + t178, t1 * qJ(6) - t2 * pkin(5) - t24 * t65 - t12 * t20 - g(1) * (-pkin(5) * t54 + qJ(6) * t55) - g(2) * (pkin(5) * t50 - qJ(6) * t51) - t223 * t159 + t295 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 + t305, t28, -t134 ^ 2 - t330, t176 - t315 + t327;];
tau_reg  = t4;
