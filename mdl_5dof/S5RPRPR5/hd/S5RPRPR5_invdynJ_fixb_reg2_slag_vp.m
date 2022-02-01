% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:26:00
% EndTime: 2022-01-23 09:26:10
% DurationCPUTime: 4.57s
% Computational Cost: add. (5356->435), mult. (13293->574), div. (0->0), fcn. (9793->14), ass. (0->248)
t200 = sin(pkin(8));
t199 = sin(pkin(9));
t201 = cos(pkin(9));
t204 = sin(qJ(3));
t207 = cos(qJ(3));
t143 = t199 * t207 + t201 * t204;
t221 = qJD(1) * t143;
t106 = t200 * t221;
t278 = qJD(1) * t200;
t256 = t204 * t278;
t233 = t199 * t256;
t276 = qJD(1) * t207;
t255 = t200 * t276;
t110 = t201 * t255 - t233;
t203 = sin(qJ(5));
t206 = cos(qJ(5));
t270 = qJD(5) * t206;
t271 = qJD(5) * t203;
t263 = t200 * qJDD(1);
t245 = t207 * t263;
t246 = t204 * t263;
t267 = qJD(1) * qJD(3);
t248 = t207 * t267;
t334 = t200 * t248 + t246;
t65 = qJD(3) * t233 - t199 * t245 - t334 * t201;
t266 = qJDD(1) * t204;
t328 = qJD(3) * t143;
t66 = (qJD(1) * t328 + t199 * t266) * t200 - t201 * t245;
t16 = t106 * t270 + t110 * t271 - t203 * t65 + t206 * t66;
t202 = cos(pkin(8));
t277 = qJD(1) * t202;
t169 = -qJD(3) + t277;
t162 = -qJD(5) + t169;
t51 = t206 * t106 + t110 * t203;
t310 = t162 * t51;
t339 = -t16 - t310;
t315 = t51 ^ 2;
t227 = -t106 * t203 + t206 * t110;
t316 = t227 ^ 2;
t338 = -t315 + t316;
t314 = t51 * t227;
t337 = t202 * t221 - t328;
t226 = t199 * t204 - t201 * t207;
t133 = t226 * qJD(3);
t336 = t226 * t277 - t133;
t17 = qJD(5) * t227 - t203 * t66 - t206 * t65;
t308 = t227 * t162;
t335 = -t17 - t308;
t152 = pkin(2) * t202 + pkin(6) * t200 + pkin(1);
t302 = qJ(2) * t207;
t168 = t202 * t302;
t105 = -t204 * t152 + t168;
t333 = qJD(3) * t105;
t321 = pkin(7) * t110;
t130 = -qJD(1) * t152 + qJD(2);
t117 = t207 * t130;
t294 = t200 * t207;
t258 = qJ(4) * t294;
t303 = qJ(2) * t204;
t260 = t202 * t303;
t219 = -t258 - t260;
t74 = qJD(1) * t219 + t117;
t64 = -t169 * pkin(3) + t74;
t257 = qJ(2) * t277;
t89 = t130 * t204 + t207 * t257;
t75 = -qJ(4) * t256 + t89;
t69 = t199 * t75;
t35 = t201 * t64 - t69;
t20 = -pkin(4) * t169 - t321 + t35;
t322 = pkin(7) * t106;
t309 = t201 * t75;
t36 = t199 * t64 + t309;
t22 = t36 - t322;
t262 = t202 * qJDD(1);
t167 = -qJDD(3) + t262;
t275 = qJD(2) * t202;
t252 = t204 * t275;
t272 = qJD(4) * t200;
t218 = -t207 * t272 - t252;
t127 = -qJDD(1) * t152 + qJDD(2);
t116 = t207 * t127;
t274 = qJD(3) * t130;
t228 = -t204 * t274 + t116;
t295 = t200 * t204;
t259 = qJ(4) * t295;
t25 = -t167 * pkin(3) + t219 * qJDD(1) + ((-t168 + t259) * qJD(3) + t218) * qJD(1) + t228;
t214 = qJD(3) * t219 - t204 * t272;
t268 = qJD(1) * qJD(2);
t249 = t207 * t268;
t273 = qJD(3) * t207;
t236 = qJDD(1) * t168 + t204 * t127 + t130 * t273 + t202 * t249;
t31 = -qJ(4) * t246 + qJD(1) * t214 + t236;
t12 = -t199 * t31 + t201 * t25;
t6 = -pkin(4) * t167 + pkin(7) * t66 + t12;
t13 = t199 * t25 + t201 * t31;
t7 = pkin(7) * t65 + t13;
t1 = (qJD(5) * t20 + t7) * t206 + t203 * t6 - t22 * t271;
t196 = qJ(3) + pkin(9);
t186 = qJ(5) + t196;
t177 = sin(t186);
t178 = cos(t186);
t208 = cos(qJ(1));
t205 = sin(qJ(1));
t292 = t202 * t205;
t100 = t177 * t208 - t178 * t292;
t291 = t202 * t208;
t102 = t177 * t205 + t178 * t291;
t317 = g(3) * t200;
t131 = pkin(3) * t256 + qJ(2) * t278 + qJD(4);
t73 = pkin(4) * t106 + t131;
t332 = g(1) * t102 - g(2) * t100 + t178 * t317 + t73 * t51 - t1;
t261 = 0.2e1 * qJ(2);
t330 = qJ(2) * qJDD(1) + t268;
t287 = t207 * t208;
t290 = t204 * t205;
t134 = t202 * t290 + t287;
t288 = t205 * t207;
t289 = t204 * t208;
t136 = -t202 * t289 + t288;
t329 = -g(1) * t136 + g(2) * t134;
t101 = -t177 * t291 + t178 * t205;
t9 = t20 * t203 + t206 * t22;
t2 = -qJD(5) * t9 - t203 * t7 + t206 * t6;
t99 = t177 * t292 + t178 * t208;
t327 = -g(1) * t101 + g(2) * t99 + t177 * t317 - t227 * t73 + t2;
t326 = t110 ^ 2;
t325 = t65 * pkin(4);
t324 = pkin(3) * t199;
t323 = pkin(3) * t204;
t319 = g(1) * t205;
t192 = g(2) * t208;
t190 = t207 * pkin(3);
t313 = qJ(4) + pkin(6);
t81 = -t143 * t203 - t206 * t226;
t312 = qJD(5) * t81 + t337 * t203 + t336 * t206;
t82 = t143 * t206 - t203 * t226;
t311 = -qJD(5) * t82 - t336 * t203 + t337 * t206;
t282 = -t152 * t273 + t207 * t275;
t62 = t214 + t282;
t63 = (-t168 + (qJ(4) * t200 + t152) * t204) * qJD(3) + t218;
t29 = t199 * t63 + t201 * t62;
t39 = t201 * t74 - t69;
t141 = t207 * t152;
t79 = -t258 - t141 + (-pkin(3) - t303) * t202;
t90 = -t259 + t105;
t44 = t199 * t79 + t201 * t90;
t88 = -t204 * t257 + t117;
t307 = t88 * t169;
t306 = t89 * t169;
t180 = pkin(3) * t201 + pkin(4);
t128 = t180 * t206 - t203 * t324;
t38 = -t199 * t74 - t309;
t26 = t38 + t322;
t27 = t39 - t321;
t305 = t128 * qJD(5) - t203 * t26 - t206 * t27;
t129 = t180 * t203 + t206 * t324;
t304 = -t129 * qJD(5) + t203 * t27 - t206 * t26;
t301 = qJDD(1) * pkin(1);
t300 = t106 * t169;
t299 = t110 * t106;
t298 = t110 * t169;
t297 = t167 * t202;
t194 = t200 ^ 2;
t209 = qJD(1) ^ 2;
t296 = t194 * t209;
t293 = t200 * t208;
t183 = sin(t196);
t147 = pkin(4) * t183 + t323;
t286 = t208 * t147;
t210 = qJ(2) ^ 2;
t239 = t268 * t261;
t264 = t194 * qJDD(1);
t283 = t194 * t239 + t210 * t264;
t184 = cos(t196);
t148 = pkin(4) * t184 + t190;
t139 = (pkin(3) * t273 + qJD(2)) * t200;
t144 = pkin(3) * t295 + t200 * qJ(2);
t187 = t205 * qJ(2);
t281 = t208 * pkin(1) + t187;
t195 = t202 ^ 2;
t280 = t194 + t195;
t197 = t204 ^ 2;
t198 = t207 ^ 2;
t279 = t197 - t198;
t265 = qJDD(1) * t207;
t254 = t106 * t278;
t253 = t110 * t278;
t250 = t204 * t268;
t244 = -t192 + t319;
t28 = -t199 * t62 + t201 * t63;
t43 = -t199 * t90 + t201 * t79;
t241 = t280 * t209;
t240 = qJD(1) * (-qJD(3) - t169);
t181 = qJDD(2) - t301;
t238 = t167 + t262;
t237 = t207 * t204 * t296;
t234 = qJD(3) * t260;
t231 = t204 * t248;
t171 = t200 * t319;
t230 = -g(2) * t293 + t171;
t229 = g(1) * t208 + g(2) * t205;
t126 = t226 * t200;
t34 = -pkin(4) * t202 + pkin(7) * t126 + t43;
t125 = t143 * t200;
t37 = -pkin(7) * t125 + t44;
t14 = -t203 * t37 + t206 * t34;
t15 = t203 * t34 + t206 * t37;
t68 = -t125 * t203 - t126 * t206;
t80 = t334 * pkin(3) + qJ(2) * t263 + t200 * t268 + qJDD(4);
t225 = qJD(3) * (t169 + t277);
t224 = t181 - t301 + t192;
t223 = t202 * t267 + t296;
t222 = t248 + t266;
t217 = -t169 ^ 2 - t296;
t215 = g(3) * t202 - t200 * t229 + t80;
t193 = -pkin(7) - t313;
t188 = t208 * qJ(2);
t179 = qJ(2) + t323;
t157 = -qJDD(5) + t167;
t145 = pkin(2) + t148;
t137 = t202 * t287 + t290;
t135 = -t202 * t288 + t289;
t124 = t313 * t200 + pkin(1) + (pkin(2) + t190) * t202;
t123 = t183 * t205 + t184 * t291;
t122 = -t183 * t291 + t184 * t205;
t121 = t183 * t208 - t184 * t292;
t120 = t183 * t292 + t184 * t208;
t113 = t200 * t328;
t109 = t200 * t133;
t104 = -t141 - t260;
t98 = t106 ^ 2;
t87 = -t252 - t333;
t86 = -t234 + t282;
t85 = pkin(3) * t255 + pkin(4) * t110;
t83 = pkin(4) * t125 + t144;
t76 = -pkin(4) * t109 + t139;
t67 = t206 * t125 - t126 * t203;
t46 = (-qJ(2) * t222 - t250) * t202 + t228;
t45 = -qJD(1) * t234 + t236;
t40 = t80 - t325;
t33 = qJD(5) * t68 - t206 * t109 - t203 * t113;
t32 = -t203 * t109 + t206 * t113 + t125 * t270 - t126 * t271;
t19 = pkin(7) * t109 + t29;
t18 = pkin(7) * t113 + t28;
t8 = t20 * t206 - t203 * t22;
t4 = -qJD(5) * t15 + t206 * t18 - t203 * t19;
t3 = qJD(5) * t14 + t203 * t18 + t206 * t19;
t5 = [0, 0, 0, 0, 0, qJDD(1), t244, t229, 0, 0, t264, 0.2e1 * t200 * t262, 0, t195 * qJDD(1), 0, 0, (-t224 + t319) * t202, t200 * t224 - t171, 0.2e1 * t330 * t280 - t229, -t181 * pkin(1) - g(1) * (-pkin(1) * t205 + t188) - g(2) * t281 + (qJDD(1) * t210 + t239) * t195 + t283, (qJDD(1) * t198 - 0.2e1 * t231) * t194, 0.2e1 * (-t204 * t265 + t279 * t267) * t194, (t204 * t225 - t207 * t238) * t200, (qJDD(1) * t197 + 0.2e1 * t231) * t194, (t204 * t238 + t207 * t225) * t200, t297, -g(1) * t135 - g(2) * t137 - t104 * t167 - t87 * t169 - t46 * t202 + (t222 * t261 + 0.2e1 * t250) * t194, -g(1) * t134 - g(2) * t136 + t105 * t167 + t86 * t169 + t45 * t202 + (0.2e1 * t249 + (-t204 * t267 + t265) * t261) * t194, t171 + (-t192 + (-qJD(3) * t89 - qJDD(1) * t104 - t46 + (-t87 - t333) * qJD(1)) * t207 + (qJD(3) * t88 - qJDD(1) * t105 - t45 + (qJD(3) * t104 - t86) * qJD(1)) * t204) * t200, t45 * t105 + t89 * t86 + t46 * t104 + t88 * t87 - g(1) * (-t152 * t205 + t188) - g(2) * (t152 * t208 + t187) + t283, -t110 * t113 + t126 * t66, t106 * t113 + t109 * t110 + t125 * t66 - t126 * t65, t113 * t169 + t126 * t167 + t202 * t66, -t106 * t109 - t125 * t65, -t109 * t169 + t125 * t167 - t202 * t65, t297, -g(1) * t121 - g(2) * t123 + t106 * t139 - t109 * t131 - t12 * t202 + t125 * t80 - t144 * t65 - t167 * t43 - t169 * t28, -g(1) * t120 - g(2) * t122 + t110 * t139 - t113 * t131 - t126 * t80 + t13 * t202 - t144 * t66 + t167 * t44 + t169 * t29, -t106 * t29 + t109 * t36 - t110 * t28 + t113 * t35 + t12 * t126 - t125 * t13 + t43 * t66 + t44 * t65 + t230, t13 * t44 + t36 * t29 + t12 * t43 + t35 * t28 + t80 * t144 + t131 * t139 - g(1) * (-t124 * t205 + t179 * t208) - g(2) * (t124 * t208 + t179 * t205), -t16 * t68 - t227 * t32, t16 * t67 - t17 * t68 - t227 * t33 + t32 * t51, -t157 * t68 + t16 * t202 + t162 * t32, t17 * t67 + t33 * t51, t157 * t67 + t162 * t33 + t17 * t202, t157 * t202, -g(1) * t100 - g(2) * t102 - t14 * t157 - t162 * t4 + t17 * t83 - t2 * t202 + t33 * t73 + t40 * t67 + t51 * t76, -g(1) * t99 - g(2) * t101 + t1 * t202 + t15 * t157 - t16 * t83 + t162 * t3 + t227 * t76 - t32 * t73 + t40 * t68, -t1 * t67 + t14 * t16 - t15 * t17 - t2 * t68 - t227 * t4 - t3 * t51 + t32 * t8 - t33 * t9 + t230, t1 * t15 + t9 * t3 + t2 * t14 + t8 * t4 + t40 * t83 + t73 * t76 - g(1) * (t188 + t286) - g(2) * (t145 * t291 - t193 * t293 + t281) + (-g(1) * (-t145 * t202 + t193 * t200 - pkin(1)) - g(2) * t147) * t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, t263, -t241, -qJ(2) * t241 + t181 - t244, 0, 0, 0, 0, 0, 0, -t207 * t167 + t204 * t217, t204 * t167 + t207 * t217, (-t197 - t198) * t263, -qJ(2) * t296 + (t46 - t306) * t207 + (t45 + t307) * t204 - t244, 0, 0, 0, 0, 0, 0, t167 * t226 - t169 * t337 - t254, t143 * t167 + t169 * t336 - t253, -t106 * t336 - t110 * t337 + t143 * t65 - t226 * t66, -t12 * t226 + t13 * t143 - t131 * t278 + t336 * t36 + t337 * t35 - t244, 0, 0, 0, 0, 0, 0, -t81 * t157 - t162 * t311 - t278 * t51, t82 * t157 + t162 * t312 - t227 * t278, t81 * t16 - t82 * t17 - t227 * t311 - t312 * t51, t1 * t82 + t2 * t81 - t278 * t73 + t311 * t8 + t312 * t9 - t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, -t279 * t296, (t204 * t240 + t265) * t200, -t237, (t207 * t240 - t266) * t200, -t167, -t306 + t116 - t223 * t302 + (-t330 * t202 - t274 + t317) * t204 + t329, g(1) * t137 - g(2) * t135 + g(3) * t294 + t223 * t303 - t236 - t307, 0, 0, t299, -t98 + t326, -t66 - t300, -t299, t65 - t298, -t167, t183 * t317 - g(1) * t122 + g(2) * t120 - t131 * t110 + t38 * t169 + (-t167 * t201 - t207 * t254) * pkin(3) + t12, t184 * t317 + g(1) * t123 - g(2) * t121 + t131 * t106 - t39 * t169 + (t167 * t199 - t207 * t253) * pkin(3) - t13, (t36 + t38) * t110 + (-t35 + t39) * t106 + (t199 * t65 + t201 * t66) * pkin(3), -t35 * t38 - t36 * t39 + (t13 * t199 + t12 * t201 + (g(3) * t204 - t131 * t276) * t200 + t329) * pkin(3), t314, t338, t339, -t314, t335, -t157, -t128 * t157 - t162 * t304 - t51 * t85 + t327, t129 * t157 + t162 * t305 - t227 * t85 + t332, t128 * t16 - t129 * t17 + (-t304 + t9) * t227 + (-t305 - t8) * t51, t1 * t129 + t2 * t128 - t73 * t85 - g(1) * (t148 * t205 - t202 * t286) - g(2) * (-t147 * t292 - t148 * t208) + t147 * t317 + t305 * t9 + t304 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 - t298, -t66 + t300, -t98 - t326, t36 * t106 + t35 * t110 + t215, 0, 0, 0, 0, 0, 0, t17 - t308, -t16 + t310, -t315 - t316, t227 * t8 + t51 * t9 + t215 - t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, t338, t339, -t314, t335, -t157, -t9 * t162 + t327, -t8 * t162 + t332, 0, 0;];
tau_reg = t5;
