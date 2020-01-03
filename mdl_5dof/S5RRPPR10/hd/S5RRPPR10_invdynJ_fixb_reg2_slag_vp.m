% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:46
% EndTime: 2019-12-31 19:44:55
% DurationCPUTime: 4.79s
% Computational Cost: add. (3718->494), mult. (8617->634), div. (0->0), fcn. (5808->8), ass. (0->250)
t174 = sin(qJ(5));
t177 = cos(qJ(5));
t172 = sin(pkin(8));
t175 = sin(qJ(2));
t277 = qJD(1) * t175;
t259 = t172 * t277;
t173 = cos(pkin(8));
t268 = t173 * qJD(2);
t110 = t259 - t268;
t257 = t173 * t277;
t275 = qJD(2) * t172;
t112 = t257 + t275;
t214 = t110 * t174 + t112 * t177;
t159 = t173 * qJDD(2);
t178 = cos(qJ(2));
t266 = qJD(1) * qJD(2);
t251 = t178 * t266;
t265 = t175 * qJDD(1);
t199 = t251 + t265;
t77 = t172 * t199 - t159;
t78 = t172 * qJDD(2) + t173 * t199;
t14 = qJD(5) * t214 + t174 * t78 - t177 * t77;
t276 = qJD(1) * t178;
t144 = qJD(5) + t276;
t301 = t214 * t144;
t342 = -t14 + t301;
t169 = g(3) * t178;
t154 = pkin(6) * t265;
t206 = qJDD(2) * pkin(2) - pkin(6) * t251 - qJDD(3) - t154;
t253 = t206 - t169;
t176 = sin(qJ(1));
t179 = cos(qJ(1));
t232 = g(1) * t179 + g(2) * t176;
t327 = t232 * t175;
t184 = -t327 - t253;
t250 = t172 * t265;
t170 = t175 ^ 2;
t171 = t178 ^ 2;
t279 = t170 - t171;
t324 = qJD(1) * t279;
t341 = qJD(2) * (t110 * t175 + t172 * t324) - (t77 + t250) * t178;
t297 = t110 * t173;
t329 = (t112 * t172 + t297) * t178;
t339 = qJD(1) * t329 - t172 * t77 + t78 * t173;
t303 = t173 * t77;
t304 = t172 * t78;
t338 = qJD(2) * t329 + t175 * (t303 + t304);
t337 = -2 * pkin(1);
t336 = t214 ^ 2;
t215 = -t110 * t177 + t112 * t174;
t335 = t215 ^ 2;
t162 = t175 * qJ(3);
t165 = t178 * pkin(2);
t262 = -pkin(1) - t165;
t210 = t262 - t162;
t105 = t210 * qJD(1);
t157 = pkin(6) * t276;
t129 = qJD(2) * qJ(3) + t157;
t55 = t105 * t173 - t129 * t172;
t40 = pkin(3) * t276 + qJD(4) - t55;
t23 = pkin(4) * t276 - pkin(7) * t112 + t40;
t56 = t105 * t172 + t129 * t173;
t45 = -qJ(4) * t276 + t56;
t25 = pkin(7) * t110 + t45;
t220 = t174 * t25 - t177 * t23;
t161 = t178 * qJDD(1);
t226 = pkin(2) * t175 - qJ(3) * t178;
t95 = qJD(2) * t226 - t175 * qJD(3);
t47 = qJD(1) * t95 + qJDD(1) * t210;
t252 = t175 * t266;
t198 = -t252 + t161;
t84 = pkin(6) * t198 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t21 = -t172 * t84 + t173 * t47;
t207 = pkin(3) * t161 + qJDD(4) - t21;
t17 = -pkin(3) * t252 + t207;
t8 = pkin(4) * t198 - t78 * pkin(7) + t17;
t22 = t172 * t47 + t173 * t84;
t15 = qJ(4) * t252 + (-qJ(4) * qJDD(1) - qJD(1) * qJD(4)) * t178 + t22;
t9 = t77 * pkin(7) + t15;
t1 = -qJD(5) * t220 + t174 * t8 + t177 * t9;
t334 = t144 * t220 + t1;
t6 = t174 * t23 + t177 * t25;
t2 = -qJD(5) * t6 - t174 * t9 + t177 * t8;
t333 = t6 * t144 + t2;
t320 = pkin(3) + pkin(4);
t332 = t320 * t77;
t270 = qJD(5) * t177;
t271 = qJD(5) * t174;
t13 = -t110 * t270 + t112 * t271 - t174 * t77 - t177 * t78;
t305 = t144 * t215;
t331 = t13 - t305;
t108 = t112 ^ 2;
t326 = -t110 ^ 2 - t108;
t274 = qJD(2) * t175;
t325 = qJ(4) * t274 - qJD(4) * t178;
t247 = t172 * qJ(4) + pkin(2);
t104 = t173 * t320 + t247;
t248 = t173 * t161;
t181 = qJD(1) ^ 2;
t295 = t171 * t181;
t322 = -t172 * t295 + (t110 + t268) * t277 - t248;
t269 = t112 * qJD(4);
t300 = t78 * qJ(4);
t186 = t206 + t269 + t300;
t10 = t186 - t332;
t321 = t10 + t327;
t319 = pkin(3) * t77;
t318 = pkin(1) * t181;
t317 = pkin(6) * t112;
t316 = g(1) * t176;
t312 = t214 * t215;
t310 = pkin(7) - qJ(3);
t115 = t172 * t174 + t173 * t177;
t261 = -pkin(6) * t172 - pkin(3);
t289 = t173 * t178;
t188 = -pkin(7) * t289 + (-pkin(4) + t261) * t175;
t119 = t226 * qJD(1);
t291 = t173 * t119;
t32 = qJD(1) * t188 - t291;
t103 = t172 * t119;
t152 = qJ(4) * t277;
t290 = t173 * t175;
t292 = t172 * t178;
t203 = -pkin(6) * t290 + pkin(7) * t292;
t42 = qJD(1) * t203 + t103 + t152;
t125 = t310 * t172;
t126 = t310 * t173;
t64 = -t125 * t177 + t126 * t174;
t309 = qJD(3) * t115 + qJD(5) * t64 - t174 * t32 - t177 * t42;
t293 = t172 * t177;
t116 = -t173 * t174 + t293;
t65 = -t125 * t174 - t126 * t177;
t308 = qJD(3) * t116 - qJD(5) * t65 + t174 * t42 - t177 * t32;
t200 = t115 * t178;
t97 = t115 * qJD(5);
t307 = -qJD(1) * t200 - t97;
t256 = t173 * t276;
t258 = t172 * t276;
t306 = t172 * t270 - t173 * t271 - t174 * t256 + t177 * t258;
t302 = t173 * t95;
t299 = pkin(6) * qJDD(1);
t298 = qJ(4) * t173;
t294 = t172 * t175;
t288 = t175 * t176;
t287 = t175 * t179;
t286 = t176 * t173;
t285 = t176 * t178;
t284 = t178 * t179;
t283 = t179 * t172;
t281 = t165 + t162;
t122 = -pkin(1) - t281;
t80 = pkin(6) * t289 + t122 * t172;
t282 = g(1) * t288 - g(2) * t287;
t280 = pkin(1) * t179 + pkin(6) * t176;
t278 = t170 + t171;
t273 = qJD(2) * t178;
t272 = qJD(4) * t172;
t264 = pkin(6) * t274;
t263 = g(1) * t284 + g(2) * t285 + g(3) * t175;
t260 = qJ(3) * t274;
t255 = qJD(4) * t290;
t91 = t110 * t276;
t16 = -t186 + t319;
t254 = -t16 - t169;
t249 = t172 * t161;
t204 = -t172 * t320 + t298;
t245 = -t204 * t276 + t157 + t272;
t137 = pkin(6) * t292;
t79 = t122 * t173 - t137;
t241 = g(1) * t173 * t287 + g(2) * t175 * t286 + qJ(3) * t249 + qJD(3) * t258;
t240 = pkin(3) * t289 + qJ(4) * t292 + t281;
t239 = pkin(2) * t284 + qJ(3) * t287 + t280;
t238 = qJ(3) * t248;
t237 = t175 * t251;
t100 = t173 * t285 - t283;
t166 = t179 * pkin(6);
t99 = t172 * t285 + t173 * t179;
t236 = -pkin(3) * t100 - t99 * qJ(4) + t166;
t101 = t178 * t283 - t286;
t235 = -g(1) * t99 + g(2) * t101;
t234 = t261 * t175;
t102 = t172 * t176 + t173 * t284;
t233 = g(1) * t100 - g(2) * t102;
t231 = -g(2) * t179 + t316;
t229 = qJD(2) * pkin(2) - pkin(6) * t277 - qJD(3);
t67 = -qJ(4) * t178 + t80;
t228 = -qJ(3) * t303 - qJD(3) * t297 - t263;
t225 = pkin(3) * t172 - t298;
t224 = pkin(6) * t110 - t172 * t229;
t164 = t178 * pkin(3);
t46 = t178 * pkin(4) + t137 + t164 + (-pkin(7) * t175 - t122) * t173;
t54 = pkin(7) * t294 + t67;
t18 = -t174 * t54 + t177 * t46;
t19 = t174 * t46 + t177 * t54;
t219 = t100 * t177 + t174 * t99;
t218 = t100 * t174 - t177 * t99;
t217 = qJ(3) * t78 + qJD(3) * t112;
t83 = t172 * t95;
t62 = -t173 * t264 + t83;
t213 = qJD(1) * (-t112 + t275);
t212 = t144 ^ 2;
t180 = qJD(2) ^ 2;
t211 = qJDD(2) * t178 - t175 * t180;
t73 = -pkin(6) * t257 + t103;
t209 = pkin(3) * t173 + t247;
t208 = pkin(6) + t225;
t197 = -pkin(6) + t204;
t195 = -t172 * t91 - t303;
t194 = t112 * qJ(4) + t229;
t192 = pkin(3) * t102 + qJ(4) * t101 + t239;
t191 = (t110 * t273 + t175 * t77) * t172;
t190 = -g(1) * t101 - g(2) * t99 - g(3) * t294;
t189 = t210 * t316;
t185 = t178 * t213 - t159 + t250;
t183 = t184 - t300;
t142 = qJ(3) * t284;
t138 = qJ(3) * t285;
t136 = t175 * t181 * t178;
t114 = -qJDD(5) - t198;
t106 = qJDD(1) * t171 - 0.2e1 * t237;
t88 = t115 * t175;
t87 = t174 * t290 - t175 * t293;
t85 = t208 * t175;
t76 = t225 * t276 + t157;
t72 = pkin(6) * t259 + t291;
t71 = t164 - t79;
t66 = t197 * t175;
t61 = t172 * t264 + t302;
t60 = qJD(1) * t234 - t291;
t59 = t152 + t73;
t58 = t208 * t273 - t255;
t48 = qJD(2) * t234 - t302;
t44 = t101 * t174 + t102 * t177;
t43 = t101 * t177 - t102 * t174;
t39 = t197 * t273 + t255;
t38 = t62 + t325;
t37 = t91 + t78;
t36 = pkin(3) * t110 - t194;
t35 = qJD(5) * t116 * t175 + qJD(2) * t200;
t34 = t174 * t178 * t268 + t175 * t97 - t273 * t293;
t33 = t173 * t295 + t175 * t213 - t249;
t29 = -t112 * t256 + t304;
t28 = qJD(2) * t203 + t325 + t83;
t27 = qJD(2) * t188 - t302;
t26 = (t112 * t273 + t175 * t78) * t173;
t24 = -t110 * t320 + t194;
t20 = (-t173 * t265 - t78) * t178 + (t112 * t175 + t173 * t324) * qJD(2);
t4 = -qJD(5) * t19 - t174 * t28 + t177 * t27;
t3 = qJD(5) * t18 + t174 * t27 + t177 * t28;
t5 = [0, 0, 0, 0, 0, qJDD(1), t231, t232, 0, 0, qJDD(1) * t170 + 0.2e1 * t237, -0.2e1 * qJD(2) * t324 + 0.2e1 * t161 * t175, qJDD(2) * t175 + t178 * t180, t106, t211, 0, (-pkin(6) * qJDD(2) + t266 * t337) * t175 + (0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t180 + t231) * t178, -pkin(6) * t211 + t199 * t337 - t282, 0.2e1 * t278 * t299 - t232, -g(1) * (-t176 * pkin(1) + t166) - g(2) * t280 + (pkin(6) ^ 2 * t278 + (pkin(1) ^ 2)) * qJDD(1), t26, -t338, t20, t191, -t341, t106, (pkin(6) * t77 - t206 * t172 + (qJD(1) * t79 + t55) * qJD(2)) * t175 + (-qJD(1) * t61 + qJD(2) * t224 - qJDD(1) * t79 - t21) * t178 + t233, (pkin(6) * t78 - t206 * t173 + (-qJD(1) * t80 - t56) * qJD(2)) * t175 + (qJD(1) * t62 + qJDD(1) * t80 + t22 + (-t173 * t229 + t317) * qJD(2)) * t178 + t235, -t62 * t110 - t61 * t112 - t80 * t77 - t79 * t78 + (-t172 * t22 - t173 * t21) * t175 + (-t172 * t56 - t173 * t55) * t273 + t282, t22 * t80 + t56 * t62 + t21 * t79 + t55 * t61 - g(1) * t166 - g(2) * t239 - t189 + (-t175 * t206 - t229 * t273) * pkin(6), t26, t20, t338, t106, t341, t191, t110 * t58 + t77 * t85 + (t16 * t172 + (-qJD(1) * t71 - t40) * qJD(2)) * t175 + (qJD(1) * t48 + qJDD(1) * t71 + t275 * t36 + t17) * t178 + t233, -t38 * t110 + t48 * t112 - t67 * t77 + t71 * t78 + (-t15 * t172 + t17 * t173) * t175 + (-t172 * t45 + t173 * t40) * t273 + t282, -t112 * t58 - t78 * t85 + (-t16 * t173 + (qJD(1) * t67 + t45) * qJD(2)) * t175 + (-qJD(1) * t38 - qJDD(1) * t67 - t268 * t36 - t15) * t178 - t235, -g(1) * t236 - g(2) * t192 + t15 * t67 + t16 * t85 + t17 * t71 + t36 * t58 + t45 * t38 + t40 * t48 - t189, -t13 * t88 + t214 * t35, t13 * t87 - t14 * t88 - t214 * t34 - t215 * t35, -t114 * t88 - t13 * t178 + t144 * t35 - t214 * t274, t14 * t87 + t215 * t34, t114 * t87 - t14 * t178 - t144 * t34 + t215 * t274, -t114 * t178 - t144 * t274, g(1) * t219 - g(2) * t44 + t10 * t87 - t114 * t18 + t14 * t66 + t144 * t4 + t178 * t2 + t215 * t39 + t220 * t274 + t24 * t34, -g(1) * t218 - g(2) * t43 - t1 * t178 + t10 * t88 + t114 * t19 - t13 * t66 - t144 * t3 + t214 * t39 + t24 * t35 + t274 * t6, -t1 * t87 + t13 * t18 - t14 * t19 - t2 * t88 - t214 * t4 - t215 * t3 + t220 * t35 - t34 * t6 - t282, t1 * t19 + t6 * t3 + t2 * t18 - t220 * t4 + t10 * t66 + t24 * t39 - g(1) * (-t100 * pkin(4) + t236) - g(2) * (pkin(4) * t102 - pkin(7) * t287 + t192) - (t175 * t310 + t262) * t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, t279 * t181, t265, t136, t161, qJDD(2), -t154 - t169 + (t232 + t318) * t175, (-t299 + t318) * t178 + t263, 0, 0, t29, t339, t33, t195, t322, t136, -pkin(2) * t77 + t253 * t173 + ((-qJ(3) * t275 - t55) * t175 + (-t224 + t72) * t178) * qJD(1) + t241, t238 - pkin(2) * t78 + t184 * t172 + ((-qJ(3) * t268 + t56) * t175 + (-t317 - t73 + (qJD(3) + t229) * t173) * t178) * qJD(1), t110 * t73 + t112 * t72 + (t276 * t55 + t22) * t173 + (t276 * t56 - t21 + t217) * t172 + t228, t206 * pkin(2) - t56 * t73 - t55 * t72 + t229 * t157 - g(1) * (-pkin(2) * t287 + t142) - g(2) * (-pkin(2) * t288 + t138) - g(3) * t281 + (-t172 * t55 + t173 * t56) * qJD(3) + (-t172 * t21 + t173 * t22) * qJ(3), t29, t33, -t339, t136, -t322, t195, -t209 * t77 + t254 * t173 + (-t76 - t272) * t110 + (t175 * t40 - t178 * t60 + (-t178 * t36 - t260) * t172) * qJD(1) + t241, t110 * t59 - t112 * t60 + (-t276 * t40 + t15) * t173 + (t276 * t45 + t17 + t217) * t172 + t228, -t238 + t112 * t76 + t209 * t78 + (t254 + t269 + t327) * t172 + (-t175 * t45 + t178 * t59 + (t260 + (-qJD(3) + t36) * t178) * t173) * qJD(1), -t45 * t59 - t36 * t76 - t40 * t60 - g(1) * t142 - g(2) * t138 - g(3) * t240 + (qJ(3) * t15 + qJD(3) * t45) * t173 + (qJ(3) * t17 + qJD(3) * t40 - qJD(4) * t36) * t172 + (-t16 + t327) * t209, -t116 * t13 + t214 * t307, t13 * t115 - t116 * t14 - t214 * t306 - t215 * t307, -t114 * t116 + t144 * t307 + t214 * t277, t14 * t115 + t215 * t306, t114 * t115 - t144 * t306 - t215 * t277, t144 * t277, -g(3) * t200 + t104 * t14 - t64 * t114 + t115 * t321 + t308 * t144 + t215 * t245 - t220 * t277 + t306 * t24, -t6 * t277 - t104 * t13 + t65 * t114 - t309 * t144 + t307 * t24 + t245 * t214 + (-t169 + t321) * t116, -t1 * t115 - t116 * t2 + t13 * t64 - t14 * t65 - t214 * t308 - t215 * t309 + t220 * t307 - t306 * t6 + t263, t1 * t65 + t2 * t64 + t10 * t104 - g(1) * (-pkin(7) * t284 + t142) - g(2) * (-pkin(7) * t285 + t138) - g(3) * (pkin(4) * t289 + t240) + t309 * t6 - t308 * t220 + t245 * t24 + (g(3) * pkin(7) + t104 * t232) * t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t37, t326, t56 * t110 + t55 * t112 + t184, 0, 0, 0, 0, 0, 0, t185, t326, -t37, t319 + t110 * t45 + (-qJD(4) - t40) * t112 + t183, 0, 0, 0, 0, 0, 0, -t14 - t301, t13 + t305, t335 + t336, t214 * t220 - t215 * t6 + t183 - t269 + t332; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110 * t112 + t198, -t91 + t78, -t108 - t295, t36 * t112 + (-pkin(3) * t274 + t178 * t45) * qJD(1) + t190 + t207, 0, 0, 0, 0, 0, 0, -t112 * t215 - t177 * t114 - t174 * t212, -t112 * t214 + t174 * t114 - t177 * t212, t174 * t342 + t331 * t177, -t24 * t112 + t174 * t334 + t177 * t333 + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, -t335 + t336, -t331, -t312, t342, -t114, -g(1) * t43 + g(2) * t218 + g(3) * t87 - t214 * t24 + t333, g(1) * t44 + g(2) * t219 + g(3) * t88 + t215 * t24 - t334, 0, 0;];
tau_reg = t5;
