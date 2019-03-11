% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:24
% EndTime: 2019-03-09 02:01:33
% DurationCPUTime: 4.65s
% Computational Cost: add. (6951->486), mult. (14906->562), div. (0->0), fcn. (10861->14), ass. (0->243)
t180 = cos(pkin(10));
t316 = cos(qJ(4));
t250 = t316 * t180;
t152 = qJD(1) * t250;
t178 = sin(pkin(10));
t184 = sin(qJ(4));
t269 = t184 * t178;
t123 = qJD(1) * t269 - t152;
t118 = qJD(5) + t123;
t136 = t316 * t178 + t184 * t180;
t183 = sin(qJ(5));
t262 = qJD(5) * t183;
t248 = t136 * t262;
t245 = qJD(4) * t316;
t263 = qJD(4) * t184;
t249 = t178 * t263;
t127 = -t180 * t245 + t249;
t186 = cos(qJ(5));
t280 = t127 * t186;
t211 = t248 + t280;
t278 = t136 * t186;
t128 = t136 * qJD(4);
t257 = t178 * qJDD(1);
t228 = -qJDD(1) * t250 + t184 * t257;
t85 = qJD(1) * t128 + t228;
t82 = qJDD(5) + t85;
t209 = -t118 * t211 + t82 * t278;
t125 = t136 * qJD(1);
t105 = qJD(4) * t183 + t125 * t186;
t213 = t250 - t269;
t202 = -qJD(1) * t249 + qJDD(1) * t136;
t325 = qJD(4) * t152 + t202;
t332 = qJD(4) * qJD(5) + t325;
t47 = -t183 * qJDD(4) + t125 * t262 - t332 * t186;
t303 = t105 * t128 + t213 * t47;
t335 = t209 - t303;
t181 = cos(pkin(9));
t308 = t181 * pkin(1);
t160 = -pkin(2) - t308;
t258 = qJDD(1) * t160;
t138 = qJDD(3) + t258;
t177 = qJ(1) + pkin(9);
t167 = sin(t177);
t169 = cos(t177);
t244 = -g(1) * t167 + g(2) * t169;
t334 = -t138 - t244;
t331 = t125 * t128 - t213 * t325;
t179 = sin(pkin(9));
t154 = pkin(1) * t179 + qJ(3);
t146 = t154 * qJD(1);
t115 = t178 * qJD(2) + t180 * t146;
t300 = pkin(7) * qJD(1);
t102 = t180 * t300 + t115;
t165 = t180 * qJD(2);
t101 = t165 + (-t146 - t300) * t178;
t133 = qJD(1) * qJD(3) + qJDD(1) * t154;
t163 = t180 * qJDD(2);
t96 = t163 + (-pkin(7) * qJDD(1) - t133) * t178;
t108 = t178 * qJDD(2) + t180 * t133;
t256 = t180 * qJDD(1);
t97 = pkin(7) * t256 + t108;
t254 = -t101 * t245 - t184 * t96 - t316 * t97;
t21 = -t102 * t263 - t254;
t19 = qJDD(4) * pkin(8) + t21;
t261 = qJD(5) * t186;
t292 = pkin(1) * qJDD(1);
t251 = t181 * t292;
t36 = -qJDD(1) * pkin(2) - pkin(3) * t256 + t85 * pkin(4) - pkin(8) * t325 + qJDD(3) - t251;
t55 = t184 * t101 + t316 * t102;
t52 = qJD(4) * pkin(8) + t55;
t159 = t180 * pkin(3) + pkin(2);
t144 = -t159 - t308;
t120 = qJD(1) * t144 + qJD(3);
t62 = pkin(4) * t123 - pkin(8) * t125 + t120;
t242 = t183 * t19 - t186 * t36 + t52 * t261 + t62 * t262;
t318 = pkin(5) * t82;
t2 = qJDD(6) + t242 - t318;
t24 = t183 * t62 + t186 * t52;
t16 = qJ(6) * t118 + t24;
t298 = t118 * t16;
t330 = -t2 + t298;
t296 = t118 * t24;
t329 = -t242 + t296;
t176 = pkin(10) + qJ(4);
t166 = sin(t176);
t235 = g(1) * t169 + g(2) * t167;
t328 = t166 * t235;
t307 = pkin(7) + t154;
t129 = t307 * t178;
t130 = t307 * t180;
t327 = -t316 * t129 - t184 * t130;
t168 = cos(t176);
t265 = t168 * pkin(4) + t166 * pkin(8);
t247 = t136 * t261;
t281 = t127 * t183;
t212 = t247 - t281;
t279 = t136 * t183;
t197 = -t118 * t212 - t82 * t279;
t103 = -t186 * qJD(4) + t125 * t183;
t48 = -t186 * qJDD(4) + t125 * t261 + t332 * t183;
t306 = t128 * t103 - t213 * t48;
t326 = t197 - t306;
t324 = t125 * qJD(4);
t100 = t184 * t102;
t54 = t316 * t101 - t100;
t51 = -qJD(4) * pkin(4) - t54;
t27 = t103 * pkin(5) - t105 * qJ(6) + t51;
t317 = pkin(8) * t82;
t323 = t118 * t27 - t317;
t310 = g(3) * t166;
t322 = t235 * t168 + t310;
t72 = -pkin(4) * t213 - pkin(8) * t136 + t144;
t78 = -t184 * t129 + t316 * t130;
t302 = t183 * t72 + t186 * t78;
t58 = t213 * qJD(3) + qJD(4) * t327;
t84 = pkin(4) * t128 + pkin(8) * t127;
t10 = -qJD(5) * t302 - t183 * t58 + t186 * t84;
t321 = t105 ^ 2;
t320 = t118 ^ 2;
t319 = t125 ^ 2;
t185 = sin(qJ(1));
t315 = pkin(1) * t185;
t314 = pkin(8) * t105;
t309 = g(3) * t168;
t305 = -t103 * t261 - t183 * t48;
t83 = pkin(4) * t125 + pkin(8) * t123;
t30 = t183 * t83 + t186 * t54;
t37 = t48 * t278;
t304 = t103 * t280 - t37;
t301 = t127 * t123 - t136 * t85;
t299 = qJ(6) * t82;
t23 = -t183 * t52 + t186 * t62;
t297 = t118 * t23;
t295 = t183 * t47;
t73 = t183 * t82;
t294 = t186 * t48;
t74 = t186 * t82;
t229 = pkin(5) * t183 - qJ(6) * t186;
t293 = -qJD(6) * t183 + t118 * t229 - t55;
t291 = t103 * t118;
t290 = t103 * t125;
t289 = t103 * t186;
t288 = t105 * t103;
t287 = t105 * t125;
t286 = t105 * t183;
t285 = t118 * t125;
t240 = t118 * t183;
t284 = t118 * t186;
t283 = t125 * t123;
t277 = t166 * t167;
t276 = t166 * t169;
t275 = t167 * t183;
t274 = t167 * t186;
t273 = t168 * t169;
t272 = t169 * t183;
t271 = t169 * t186;
t268 = qJD(6) - t23;
t267 = (g(1) * t271 + g(2) * t274) * t166;
t187 = cos(qJ(1));
t173 = t187 * pkin(1);
t266 = t169 * t159 + t173;
t174 = t178 ^ 2;
t175 = t180 ^ 2;
t264 = t174 + t175;
t66 = t105 * t281;
t255 = t105 * t247 - t47 * t279 - t66;
t253 = pkin(8) * qJD(5) * t118;
t246 = t103 ^ 2 - t321;
t98 = t105 * t262;
t243 = -t186 * t47 - t98;
t22 = -t101 * t263 - t102 * t245 - t184 * t97 + t316 * t96;
t239 = pkin(4) * t273 + pkin(8) * t276 + t266;
t238 = -g(1) * t277 + g(2) * t276;
t110 = t168 * t275 + t271;
t112 = t168 * t272 - t274;
t237 = g(1) * t110 - g(2) * t112;
t111 = t168 * t274 - t272;
t113 = t168 * t271 + t275;
t236 = g(1) * t111 - g(2) * t113;
t233 = g(1) * t185 - g(2) * t187;
t232 = (qJD(5) * t103 - t47) * pkin(8);
t182 = -pkin(7) - qJ(3);
t231 = -t169 * t182 - t315;
t230 = pkin(5) * t186 + qJ(6) * t183;
t14 = -pkin(5) * t118 + t268;
t227 = t14 * t186 - t16 * t183;
t226 = -t14 * t183 - t16 * t186;
t225 = t183 * t24 + t186 * t23;
t224 = t183 * t23 - t186 * t24;
t29 = -t183 * t54 + t186 * t83;
t39 = -t183 * t78 + t186 * t72;
t107 = -t133 * t178 + t163;
t221 = -t107 * t178 + t108 * t180;
t220 = (-t146 * t178 + t165) * t178 - t115 * t180;
t218 = -t118 * t262 - t123 * t240 + t74;
t217 = t118 * t261 + t123 * t284 + t73;
t216 = pkin(4) + t230;
t215 = t253 + t309;
t20 = -qJDD(4) * pkin(4) - t22;
t3 = t183 * t36 + t186 * t19 + t62 * t261 - t262 * t52;
t9 = t183 * t84 + t186 * t58 + t72 * t261 - t262 * t78;
t214 = t118 * t51 - t317;
t210 = t103 * t248 + t304;
t208 = t20 + t215;
t207 = t103 * t240 - t294;
t206 = -t218 - t290;
t205 = -t258 + t334;
t117 = qJDD(1) * t144 + qJDD(3);
t204 = -t309 + t328;
t203 = g(1) * t112 + g(2) * t110 + t183 * t310 - t242;
t201 = -pkin(8) * t294 - t322;
t200 = (-g(1) * (-t159 - t265) + g(2) * t182) * t167;
t1 = qJD(6) * t118 + t299 + t3;
t199 = qJD(5) * t227 + t1 * t186 + t183 * t2;
t198 = -qJD(5) * t225 + t183 * t242 + t3 * t186;
t196 = t103 * t212 + t48 * t279;
t195 = t105 * t27 + qJDD(6) - t203;
t194 = t197 + t306;
t193 = -g(1) * t113 - g(2) * t111 - t186 * t310 + t3;
t59 = qJD(3) * t136 + qJD(4) * t78;
t190 = t105 * t118 - t48;
t141 = pkin(8) * t273;
t139 = t167 * t168 * pkin(8);
t122 = t123 ^ 2;
t87 = -qJD(4) * t128 + qJDD(4) * t213;
t86 = -qJD(4) * t127 + qJDD(4) * t136;
t60 = pkin(5) * t105 + qJ(6) * t103;
t46 = t136 * t229 - t327;
t41 = t118 * t128 - t213 * t82;
t32 = pkin(5) * t213 - t39;
t31 = -qJ(6) * t213 + t302;
t28 = -t47 + t291;
t26 = -pkin(5) * t125 - t29;
t25 = qJ(6) * t125 + t30;
t15 = t217 - t287;
t13 = t105 * t284 - t295;
t12 = -t229 * t127 + (qJD(5) * t230 - qJD(6) * t186) * t136 + t59;
t11 = -t105 * t211 - t47 * t278;
t8 = -pkin(5) * t128 - t10;
t7 = qJ(6) * t128 - qJD(6) * t213 + t9;
t6 = t209 + t303;
t5 = pkin(5) * t48 + qJ(6) * t47 - qJD(6) * t105 + t20;
t4 = [0, 0, 0, 0, 0, qJDD(1), t233, g(1) * t187 + g(2) * t185, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t244 + 0.2e1 * t251, -0.2e1 * t179 * t292 + t235, 0 (t233 + (t179 ^ 2 + t181 ^ 2) * t292) * pkin(1), t174 * qJDD(1), 0.2e1 * t178 * t256, 0, t175 * qJDD(1), 0, 0, t205 * t180, -t205 * t178, t133 * t264 + t221 - t235, t138 * t160 - g(1) * (-pkin(2) * t167 + qJ(3) * t169 - t315) - g(2) * (pkin(2) * t169 + qJ(3) * t167 + t173) + t221 * t154 - t220 * qJD(3), -t125 * t127 + t136 * t325, t301 - t331, t86, t123 * t128 - t213 * t85, t87, 0, -qJD(4) * t59 + qJDD(4) * t327 - t117 * t213 + t120 * t128 + t144 * t85 - t168 * t244, -t58 * qJD(4) - t78 * qJDD(4) + t117 * t136 - t120 * t127 + t144 * t325 + t238, -t58 * t123 + t59 * t125 + t54 * t127 - t55 * t128 - t22 * t136 + t21 * t213 - t325 * t327 - t78 * t85 - t235, t21 * t78 + t55 * t58 + t22 * t327 - t54 * t59 + t117 * t144 - g(1) * (-t159 * t167 + t231) - g(2) * (-t167 * t182 + t266) t11, t210 - t255, t6, t196, t326, t41, -t51 * t281 + t10 * t118 + t103 * t59 + t128 * t23 + t213 * t242 + t39 * t82 - t48 * t327 + (t183 * t20 + t261 * t51) * t136 + t236, -t51 * t280 + t105 * t59 - t118 * t9 - t128 * t24 + t213 * t3 - t302 * t82 + t47 * t327 + (t186 * t20 - t262 * t51) * t136 - t237, -t10 * t105 - t103 * t9 + t39 * t47 - t302 * t48 + t225 * t127 + (qJD(5) * t224 - t183 * t3 + t186 * t242) * t136 - t238, -g(1) * t231 - g(2) * t239 + t23 * t10 - t20 * t327 + t24 * t9 - t242 * t39 + t3 * t302 + t51 * t59 + t200, t11, t6, -t103 * t211 + t255 + t37, t41, -t326, t196, -t27 * t281 + t103 * t12 - t118 * t8 - t128 * t14 + t213 * t2 - t32 * t82 + t46 * t48 + (t183 * t5 + t261 * t27) * t136 + t236, -t103 * t7 + t105 * t8 - t31 * t48 - t32 * t47 - t227 * t127 + (qJD(5) * t226 - t1 * t183 + t186 * t2) * t136 - t238, t27 * t280 - t1 * t213 - t105 * t12 + t118 * t7 + t128 * t16 + t31 * t82 + t46 * t47 + (-t186 * t5 + t262 * t27) * t136 + t237, t1 * t31 + t16 * t7 + t5 * t46 + t27 * t12 + t2 * t32 + t14 * t8 - g(1) * (-pkin(5) * t111 - qJ(6) * t110 + t231) - g(2) * (pkin(5) * t113 + qJ(6) * t112 + t239) + t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t107 * t180 + t108 * t178 - g(3), 0, 0, 0, 0, 0, 0, t87, -t86, t301 + t331, -t127 * t55 - t128 * t54 + t136 * t21 + t213 * t22 - g(3), 0, 0, 0, 0, 0, 0, t194, -t335, -t66 + (-t295 + (t103 * t183 + t105 * t186) * qJD(5)) * t136 + t304, t127 * t224 + t128 * t51 + t136 * t198 - t20 * t213 - g(3), 0, 0, 0, 0, 0, 0, t194, t210 + t255, t335, t127 * t226 + t128 * t27 + t136 * t199 - t213 * t5 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256, t257, -t264 * qJD(1) ^ 2, qJD(1) * t220 - t334, 0, 0, 0, 0, 0, 0, t228 + 0.2e1 * t324 (t152 - t123) * qJD(4) + t202, -t122 - t319, t123 * t55 + t125 * t54 + t117 + t244, 0, 0, 0, 0, 0, 0, t218 - t290, -t186 * t320 - t287 - t73 (-t103 * t123 + t47) * t186 + t105 * t240 + t305, -t125 * t51 + t329 * t186 + (t3 - t297) * t183 + t244, 0, 0, 0, 0, 0, 0, -t183 * t320 - t290 + t74 (t286 - t289) * t123 - t243 + t305, t217 + t287, -t125 * t27 + t330 * t186 + (t118 * t14 + t1) * t183 + t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, -t122 + t319 (t152 + t123) * qJD(4) + t202, -t283, -t228, qJDD(4), qJD(4) * t55 - t120 * t125 + t204 + t22, t120 * t123 + (t54 + t100) * qJD(4) + t254 + t322, 0, 0, t13 (-t286 - t289) * t123 + t243 + t305, t15, t207, -t206, -t285, -pkin(4) * t48 - t103 * t55 - t118 * t29 - t125 * t23 + t183 * t214 - t186 * t208 + t267, pkin(4) * t47 - t105 * t55 + t118 * t30 + t125 * t24 + t214 * t186 + (t208 - t328) * t183, t103 * t30 + t105 * t29 + (-t123 * t23 + t3 + (-t23 + t314) * qJD(5)) * t186 + (t232 - t329) * t183 + t201, -t20 * pkin(4) - t24 * t30 - t23 * t29 - t51 * t55 - g(1) * (-pkin(4) * t276 + t141) - g(2) * (-pkin(4) * t277 + t139) - g(3) * t265 + t198 * pkin(8), t13, t15, t98 + (t105 * t123 + t48) * t183 + (t47 + t291) * t186, -t285, t206, t207, t118 * t26 + t125 * t14 - t216 * t48 + t293 * t103 + (-t215 - t5) * t186 + t323 * t183 + t267, t103 * t25 - t105 * t26 + (t123 * t14 + t1 + (t14 + t314) * qJD(5)) * t186 + (t232 - t330) * t183 + t201, -t118 * t25 - t125 * t16 - t216 * t47 - t293 * t105 - t323 * t186 + (t204 - t5 - t253) * t183, -t16 * t25 - t14 * t26 - g(1) * t141 - g(2) * t139 - g(3) * (t168 * t230 + t265) + t293 * t27 + t199 * pkin(8) + (-t5 + t328) * t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, -t246, t28, -t288, t190, t82, -t105 * t51 + t203 + t296, t103 * t51 - t193 + t297, 0, 0, t288, t28, t246, t82, -t190, -t288, -t103 * t60 - t195 + t296 + 0.2e1 * t318, pkin(5) * t47 - qJ(6) * t48 + (t16 - t24) * t105 + (t14 - t268) * t103, 0.2e1 * t299 - t103 * t27 + t105 * t60 + (0.2e1 * qJD(6) - t23) * t118 + t193, t1 * qJ(6) - t2 * pkin(5) - t27 * t60 - t14 * t24 - g(1) * (-pkin(5) * t112 + qJ(6) * t113) - g(2) * (-pkin(5) * t110 + qJ(6) * t111) + t229 * t310 + t268 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(5) - t228 + t288 - t324, t28, -t320 - t321, t195 - t298 - t318;];
tau_reg  = t4;
