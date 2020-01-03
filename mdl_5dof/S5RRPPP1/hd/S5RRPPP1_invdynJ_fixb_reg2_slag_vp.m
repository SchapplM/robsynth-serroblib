% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPP1
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:36
% EndTime: 2019-12-31 19:24:45
% DurationCPUTime: 4.61s
% Computational Cost: add. (4267->495), mult. (12196->641), div. (0->0), fcn. (9194->8), ass. (0->236)
t189 = sin(qJ(2));
t191 = cos(qJ(2));
t272 = t191 * qJDD(1);
t275 = qJD(1) * qJD(2);
t338 = t189 * t275 - t272;
t186 = sin(pkin(5));
t302 = t186 * t189;
t170 = qJ(3) * t302;
t287 = t191 * pkin(2) + t170;
t337 = -pkin(1) - t287;
t185 = sin(pkin(8));
t187 = cos(pkin(8));
t188 = cos(pkin(5));
t297 = t188 * t189;
t126 = t185 * t191 + t187 * t297;
t113 = t126 * qJD(2);
t295 = t188 * t191;
t266 = t187 * t295;
t273 = t189 * qJDD(1);
t274 = qJDD(2) * t186;
t227 = -qJDD(1) * t266 + t185 * t273 - t187 * t274;
t58 = qJD(1) * t113 + t227;
t282 = qJD(1) * t191;
t261 = t188 * t282;
t283 = qJD(1) * t189;
t263 = t185 * t283;
t280 = qJD(2) * t187;
t82 = -t186 * t280 - t187 * t261 + t263;
t336 = t58 * qJ(5) + t82 * qJD(5);
t262 = t186 * t283;
t277 = qJD(3) * t187;
t300 = t186 * t191;
t269 = qJ(3) * t300;
t322 = pkin(2) * t189;
t228 = -t269 + t322;
t135 = t228 * qJD(1);
t256 = qJ(3) * t188 + pkin(7);
t240 = qJD(1) * t256;
t136 = t189 * t240;
t137 = t191 * t240;
t304 = t185 * t188;
t305 = t185 * t186;
t46 = t135 * t305 - t187 * t136 - t137 * t304;
t310 = qJ(4) * t262 - qJD(4) * t188 - t186 * t277 + t46;
t141 = -t188 * qJD(2) + t186 * t282;
t125 = t185 * t295 + t187 * t189;
t255 = t188 * t263;
t258 = t191 * t275;
t59 = -qJD(2) * t255 + qJDD(1) * t125 + t185 * t274 + t187 * t258;
t281 = qJD(2) * t186;
t143 = t261 + t281;
t331 = t186 ^ 2 + t188 ^ 2;
t66 = t143 * t187 - t331 * t263;
t335 = -t141 * t66 + t59;
t138 = t141 ^ 2;
t83 = qJD(1) * t125 + t185 * t281;
t334 = -t83 ^ 2 - t138;
t145 = t256 * t189;
t333 = t187 * (-t145 * t188 + t186 * t337);
t316 = pkin(3) + qJ(5);
t95 = qJDD(2) * t188 + t338 * t186;
t332 = pkin(4) * t59 - t316 * t95;
t192 = cos(qJ(1));
t290 = t192 * t188;
t190 = sin(qJ(1));
t294 = t189 * t190;
t132 = t186 * t294 - t290;
t293 = t189 * t192;
t296 = t188 * t190;
t133 = t186 * t293 + t296;
t330 = -g(1) * t133 - g(2) * t132 + g(3) * t300;
t279 = qJD(2) * t189;
t276 = qJD(3) * t189;
t96 = qJD(2) * t228 - t186 * t276;
t239 = qJD(2) * t256;
t97 = qJD(3) * t295 - t189 * t239;
t98 = -t188 * t276 - t191 * t239;
t36 = t187 * t97 + t98 * t304 + t96 * t305;
t23 = -t186 * (qJ(4) * t279 - qJD(4) * t191) - t36;
t102 = pkin(7) * t282 + qJ(3) * t143;
t119 = qJD(2) * pkin(2) - t136;
t120 = t337 * qJD(1);
t39 = -t185 * t102 + (t119 * t188 + t120 * t186) * t187;
t52 = t143 * qJD(3) - t338 * pkin(7) + (-t188 * t338 + t274) * qJ(3);
t63 = qJDD(2) * pkin(2) + qJD(1) * t98 - qJDD(1) * t145;
t64 = qJD(1) * t96 + qJDD(1) * t337;
t9 = -t185 * t52 + (t186 * t64 + t188 * t63) * t187;
t268 = t185 * t297;
t114 = -qJD(2) * t268 + t191 * t280;
t8 = t186 * (-t191 * t59 + t83 * t279) - t114 * t141 + t125 * t95;
t124 = t185 * t189 - t266;
t7 = t186 * (-t191 * t58 + t82 * t279) - t141 * t113 + t95 * t124;
t111 = t126 * qJD(1);
t329 = t186 * (t187 * t95 + t82 * t283) - t141 * t111 - t188 * t58;
t112 = t187 * t282 - t255;
t328 = t186 * (-t185 * t95 + t83 * t283) - t112 * t141 - t59 * t188;
t327 = -t111 * t83 - t82 * t112 + (t185 * t58 - t187 * t59) * t186;
t323 = t95 * pkin(3);
t320 = g(1) * t190;
t318 = g(3) * t189;
t317 = t82 * t83;
t65 = t331 * t187 * t283 + t185 * t143;
t315 = t141 * t65;
t313 = t187 * t96;
t67 = t188 * t135 + t186 * t137;
t221 = -t112 * qJ(4) + t67;
t312 = -t316 * t111 + (-qJD(4) * t185 - qJD(5) * t187) * t186 - t221;
t116 = t185 * t136;
t298 = t187 * t188;
t245 = t137 * t298 - t116;
t260 = t316 * t189;
t278 = qJD(3) * t185;
t306 = t135 * t187;
t311 = -pkin(4) * t112 - qJD(5) * t188 - t245 + (qJD(1) * t260 + t278 + t306) * t186;
t309 = pkin(4) * t111 - t310;
t308 = pkin(7) * qJDD(1);
t307 = qJD(4) * t83;
t303 = t186 * t187;
t301 = t186 * t190;
t299 = t186 * t192;
t292 = t190 * t191;
t291 = t191 * t192;
t146 = t256 * t191;
t130 = t185 * t146;
t289 = pkin(3) * t300 + t130;
t129 = pkin(2) * t304 + qJ(3) * t303;
t181 = t192 * pkin(7);
t288 = qJ(3) * t290 + t181;
t286 = t192 * pkin(1) + t190 * pkin(7);
t183 = t189 ^ 2;
t184 = t191 ^ 2;
t285 = t183 - t184;
t284 = t183 + t184;
t10 = t187 * t52 + t63 * t304 + t64 * t305;
t40 = t187 * t102 + t119 * t304 + t120 * t305;
t271 = pkin(2) * t294;
t270 = pkin(2) * t293;
t194 = qJD(1) ^ 2;
t267 = t189 * t194 * t191;
t51 = -t145 * t304 + t187 * t146 + t305 * t337;
t265 = -pkin(2) * t187 - pkin(3);
t257 = -qJ(4) * t185 - pkin(2);
t53 = -t186 * t98 + t188 * t96;
t68 = t186 * t145 + t188 * t337;
t254 = t189 * t258;
t72 = t185 * t292 + (t188 * t294 + t299) * t187;
t74 = t126 * t192 - t187 * t301;
t253 = g(1) * t72 - g(2) * t74;
t73 = -t185 * t299 + t187 * t292 - t190 * t268;
t75 = t185 * t301 + t187 * t291 - t192 * t268;
t252 = g(1) * t73 - g(2) * t75;
t85 = t185 * t97;
t251 = -t98 * t298 + t85;
t5 = -t95 * qJ(4) + t141 * qJD(4) - t10;
t103 = t185 * t294 - t190 * t266;
t104 = t125 * t190;
t158 = t190 * t269;
t250 = -t104 * pkin(3) - qJ(4) * t103 + t158;
t105 = t185 * t293 - t192 * t266;
t106 = t125 * t192;
t160 = t192 * t269;
t249 = -t106 * pkin(3) - qJ(4) * t105 + t160;
t248 = g(1) * t132 - g(2) * t133;
t247 = g(1) * t192 + g(2) * t190;
t246 = -g(2) * t192 + t320;
t244 = -t65 * t83 + t66 * t82;
t107 = -qJ(4) * t188 - t129;
t28 = -t186 * t63 + t188 * t64 + qJDD(3);
t62 = -t186 * t119 + t188 * t120 + qJD(3);
t20 = t113 * t82 + t124 * t58;
t21 = t114 * t83 + t125 * t59;
t238 = pkin(2) * t291 + qJ(3) * t296 + t192 * t170 + t286;
t32 = qJ(4) * t141 - t40;
t231 = -t73 * pkin(3) - qJ(4) * t72 + t288;
t230 = pkin(4) * t300 - t322;
t127 = t187 * t191 - t268;
t229 = t127 * pkin(3) + qJ(4) * t126 + t287;
t24 = -t111 * t82 - t58 * t303;
t25 = -t112 * t83 + t59 * t305;
t222 = -0.2e1 * pkin(1) * t275 - pkin(7) * qJDD(2);
t220 = -t125 * qJ(4) + t68;
t44 = qJ(4) * t300 - t51;
t216 = g(1) * t105 + g(2) * t103 - g(3) * t126;
t215 = g(1) * t106 + g(2) * t104 - g(3) * t127;
t214 = -t83 * qJ(4) + t62;
t3 = -pkin(4) * t58 + qJDD(5) - t5;
t213 = t75 * pkin(3) + qJ(4) * t74 + t238;
t212 = -qJ(4) * t114 - qJD(4) * t125 + t53;
t211 = t337 * t320;
t210 = -t191 * t247 - t318;
t209 = qJDD(4) - t9;
t208 = qJD(4) - t39;
t193 = qJD(2) ^ 2;
t207 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t193 + t246;
t206 = pkin(1) * t194 + t247 - t308;
t204 = t113 * t83 + t114 * t82 + t124 * t59 + t125 * t58;
t203 = t95 - t317;
t4 = pkin(3) * t58 - qJ(4) * t59 + t28 - t307;
t201 = -t126 * t275 - t227;
t200 = -g(1) * t74 - g(2) * t72 - g(3) * t124 + t209;
t199 = t4 + t330;
t198 = t58 - t315;
t196 = -t141 * t82 + t59;
t167 = qJ(3) * t305;
t128 = pkin(2) * t298 - t167;
t109 = (-pkin(3) * t187 + t257) * t186;
t108 = t265 * t188 + t167;
t80 = (-t316 * t187 + t257) * t186;
t79 = pkin(4) * t303 - t107;
t69 = pkin(4) * t305 + t167 + (-qJ(5) + t265) * t188;
t61 = t141 * t262 + t188 * t95;
t56 = (-t141 * t279 - t191 * t95) * t186;
t50 = -t130 + t333;
t47 = t289 - t333;
t45 = t116 + (t135 * t186 - t137 * t188) * t187;
t43 = pkin(3) * t124 + t220;
t42 = (-pkin(3) * t283 - t306) * t186 + t245;
t38 = pkin(3) * t111 + t221;
t37 = -pkin(4) * t124 - t44;
t35 = -t85 + (t186 * t96 + t188 * t98) * t187;
t34 = t316 * t124 + t220;
t33 = t145 * t298 + pkin(4) * t125 + (qJ(5) * t191 - t187 * t337) * t186 + t289;
t31 = t141 * pkin(3) + t208;
t29 = (-pkin(3) * t279 - t313) * t186 + t251;
t22 = pkin(3) * t82 + t214;
t19 = pkin(3) * t113 + t212;
t18 = -pkin(4) * t82 + qJD(5) - t32;
t17 = -pkin(4) * t113 - t23;
t16 = pkin(4) * t114 + (-qJD(2) * t260 + qJD(5) * t191 - t313) * t186 + t251;
t15 = t316 * t82 + t214;
t14 = pkin(4) * t83 + t316 * t141 + t208;
t11 = qJD(5) * t124 + t316 * t113 + t212;
t6 = t209 - t323;
t2 = qJD(5) * t141 + t209 + t332;
t1 = t4 + t336;
t12 = [0, 0, 0, 0, 0, qJDD(1), t246, t247, 0, 0, qJDD(1) * t183 + 0.2e1 * t254, 0.2e1 * t189 * t272 - 0.2e1 * t285 * t275, qJDD(2) * t189 + t191 * t193, qJDD(1) * t184 - 0.2e1 * t254, qJDD(2) * t191 - t189 * t193, 0, t189 * t222 + t191 * t207, -t189 * t207 + t191 * t222, 0.2e1 * t284 * t308 - t247, -g(1) * (-pkin(1) * t190 + t181) - g(2) * t286 + (t284 * pkin(7) ^ 2 + pkin(1) ^ 2) * qJDD(1), t21, -t204, t8, t20, -t7, t56, t113 * t62 + t124 * t28 - t141 * t35 + t50 * t95 + t53 * t82 + t58 * t68 + (-t191 * t9 + t39 * t279) * t186 + t252, t114 * t62 + t125 * t28 + t141 * t36 - t51 * t95 + t53 * t83 + t59 * t68 + (t10 * t191 - t40 * t279) * t186 - t253, -t10 * t124 - t113 * t40 - t114 * t39 - t125 * t9 - t35 * t83 - t36 * t82 - t50 * t59 - t51 * t58 + t248, -g(1) * t288 - g(2) * t238 + t10 * t51 + t28 * t68 + t39 * t35 + t40 * t36 + t9 * t50 + t62 * t53 - t211, t56, -t8, t7, t21, -t204, t20, t113 * t32 + t114 * t31 + t124 * t5 + t125 * t6 + t23 * t82 + t29 * t83 + t44 * t58 + t47 * t59 + t248, -t113 * t22 - t124 * t4 - t141 * t29 - t19 * t82 - t43 * t58 + t47 * t95 + (-t191 * t6 + t31 * t279) * t186 - t252, -t114 * t22 - t125 * t4 + t141 * t23 - t19 * t83 - t43 * t59 - t44 * t95 + (t191 * t5 - t32 * t279) * t186 + t253, -g(1) * t231 - g(2) * t213 + t22 * t19 + t32 * t23 + t31 * t29 + t4 * t43 + t5 * t44 + t6 * t47 - t211, t56, t7, t8, t20, t204, t21, -t113 * t18 + t114 * t14 - t124 * t3 + t125 * t2 + t16 * t83 - t17 * t82 + t33 * t59 - t37 * t58 + t248, -t1 * t125 - t11 * t83 - t114 * t15 - t141 * t17 - t34 * t59 + t37 * t95 + (t18 * t279 - t191 * t3) * t186 + t253, t1 * t124 + t11 * t82 + t113 * t15 + t141 * t16 - t33 * t95 + t34 * t58 + (-t14 * t279 + t191 * t2) * t186 + t252, t1 * t34 + t15 * t11 + t2 * t33 + t14 * t16 + t3 * t37 + t18 * t17 - g(1) * (-pkin(4) * t132 - qJ(5) * t73 + t231) - g(2) * (pkin(4) * t133 + qJ(5) * t75 + t213) - t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t267, t285 * t194, t273, t267, t272, qJDD(2), -g(3) * t191 + t189 * t206, t191 * t206 + t318, 0, 0, t25, -t327, -t328, t24, t329, t61, -t111 * t62 + t128 * t95 + t141 * t45 + t188 * t9 - t67 * t82 + (-pkin(2) * t58 + t141 * t278 - t187 * t28 - t39 * t283) * t186 + t215, -t10 * t188 - t112 * t62 - t129 * t95 - t141 * t46 - t67 * t83 + (-pkin(2) * t59 + t141 * t277 + t185 * t28 + t40 * t283) * t186 - t216, t111 * t40 + t112 * t39 - t128 * t59 - t129 * t58 + t45 * t83 + t46 * t82 + (t10 * t187 - t185 * t9 + (t185 * t83 - t187 * t82) * qJD(3) + t210) * t186, t10 * t129 + t9 * t128 - t40 * t46 - t39 * t45 - t62 * t67 - g(1) * (t160 - t270) - g(2) * (t158 - t271) - g(3) * t287 + (-t28 * pkin(2) + (-t185 * t39 + t187 * t40) * qJD(3)) * t186, t61, t328, -t329, t25, -t327, t24, t107 * t58 + t108 * t59 - t111 * t32 - t112 * t31 - t42 * t83 + t310 * t82 + (-t187 * t5 + (qJD(3) * t83 + t6) * t185 + t210) * t186, t108 * t95 - t109 * t58 + t111 * t22 + t141 * t42 + t188 * t6 + t38 * t82 + (-t31 * t283 + t187 * t4 + (-qJD(3) * t141 + qJD(4) * t82) * t185) * t186 - t215, -t107 * t95 - t109 * t59 + t112 * t22 - t188 * t5 + t38 * t83 + t310 * t141 + (t32 * t283 + (-t4 + t307) * t185) * t186 + t216, t4 * t109 + t5 * t107 + t6 * t108 - t22 * t38 - t31 * t42 - g(1) * (t249 - t270) - g(2) * (t250 - t271) - g(3) * t229 + t310 * t32 + (qJD(3) * t31 - qJD(4) * t22) * t305, t61, -t329, -t328, t24, t327, t25, t111 * t18 - t112 * t14 - t58 * t79 + t59 * t69 + t311 * t83 - t309 * t82 + (t185 * t2 + t187 * t3 + t210) * t186, t112 * t15 + t188 * t3 - t59 * t80 + t79 * t95 - t312 * t83 + (-t1 * t185 - t18 * t283) * t186 - t309 * t141 + t216, -t111 * t15 - t188 * t2 + t58 * t80 - t69 * t95 + t312 * t82 + (-t1 * t187 + t14 * t283) * t186 + t311 * t141 + t215, t1 * t80 + t2 * t69 + t3 * t79 - g(1) * (-qJ(5) * t106 + t192 * t230 + t249) - g(2) * (-qJ(5) * t104 + t190 * t230 + t250) - g(3) * (pkin(4) * t302 + qJ(5) * t127 + t229) + t309 * t18 + t312 * t15 + t311 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t335, t244, t39 * t65 - t40 * t66 + t28 + t330, 0, 0, 0, 0, 0, 0, t244, t201 + t315, -t335, -t31 * t65 + t32 * t66 + t199, 0, 0, 0, 0, 0, 0, t244, -t335, t198, -t14 * t65 - t18 * t66 + t199 + t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, t203, t334, -t141 * t32 + t22 * t83 + t200 - t323, 0, 0, 0, 0, 0, 0, t196, t334, -t203, t15 * t83 + (qJD(5) + t18) * t141 + t200 + t332; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141 * t83 + t201, t95 + t317, -t82 ^ 2 - t138, -g(1) * t75 - g(2) * t73 - g(3) * t125 - t14 * t141 - t15 * t82 + t3;];
tau_reg = t12;
