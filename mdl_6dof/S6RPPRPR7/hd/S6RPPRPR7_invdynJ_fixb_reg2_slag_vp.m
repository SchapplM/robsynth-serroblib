% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:55
% EndTime: 2019-03-09 01:54:02
% DurationCPUTime: 5.55s
% Computational Cost: add. (8367->486), mult. (17147->603), div. (0->0), fcn. (12656->14), ass. (0->241)
t198 = cos(pkin(9));
t205 = cos(qJ(4));
t281 = t205 * t198;
t258 = qJD(1) * t281;
t203 = sin(qJ(4));
t196 = sin(pkin(9));
t271 = qJD(1) * t196;
t259 = t203 * t271;
t137 = t258 - t259;
t195 = sin(pkin(10));
t197 = cos(pkin(10));
t113 = qJD(4) * t195 + t137 * t197;
t202 = sin(qJ(6));
t322 = cos(qJ(6));
t127 = t195 * t137;
t337 = qJD(4) * t197 - t127;
t222 = t322 * t337;
t56 = -t113 * t202 + t222;
t345 = t56 ^ 2;
t147 = t196 * t205 + t198 * t203;
t134 = t147 * qJD(1);
t131 = qJD(6) + t134;
t344 = t56 * t131;
t55 = t113 * t322 + t202 * t337;
t343 = t55 ^ 2;
t191 = pkin(9) + qJ(4);
t181 = cos(t191);
t179 = sin(t191);
t319 = g(3) * t179;
t204 = sin(qJ(1));
t206 = cos(qJ(1));
t334 = g(1) * t204 - g(2) * t206;
t215 = t181 * t334 - t319;
t201 = -pkin(1) - qJ(3);
t326 = -qJD(1) * qJD(3) + qJDD(1) * t201;
t150 = qJDD(2) + t326;
t253 = -pkin(7) * qJDD(1) + t150;
t118 = t253 * t196;
t119 = t253 * t198;
t162 = qJD(1) * t201 + qJD(2);
t255 = -pkin(7) * qJD(1) + t162;
t128 = t255 * t196;
t129 = t255 * t198;
t269 = qJD(4) * t205;
t270 = qJD(4) * t203;
t250 = t118 * t203 - t119 * t205 + t128 * t269 + t129 * t270;
t331 = -qJDD(4) * pkin(4) + qJDD(5);
t34 = t250 + t331;
t342 = t34 + t215;
t318 = g(3) * t181;
t221 = t179 * t334 + t318;
t261 = t322 * t197;
t223 = -t195 * t202 + t261;
t268 = qJD(6) * t202;
t332 = qJD(6) * t261 - t195 * t268;
t299 = -t134 * t223 - t332;
t148 = t195 * t322 + t197 * t202;
t140 = t148 * qJD(6);
t298 = t134 * t148 + t140;
t141 = t147 * qJD(4);
t149 = -t196 * t203 + t281;
t125 = t203 * t128;
t80 = t205 * t129 - t125;
t65 = -qJD(4) * pkin(4) + qJD(5) - t80;
t213 = -t65 * t141 + t34 * t149 + t334;
t339 = t137 * t337;
t188 = t196 ^ 2;
t189 = t198 ^ 2;
t272 = t188 + t189;
t338 = t162 * t272;
t90 = t223 * t147;
t246 = g(1) * t206 + g(2) * t204;
t224 = t246 * t179;
t316 = -pkin(7) + t201;
t151 = t316 * t196;
t152 = t316 * t198;
t336 = -t151 * t203 + t152 * t205;
t182 = t196 * pkin(3);
t200 = -pkin(7) - qJ(3);
t335 = t182 * t206 + t200 * t204;
t192 = qJDD(1) * qJ(2);
t193 = qJD(1) * qJD(2);
t333 = t192 + t193;
t156 = qJDD(3) + t333;
t219 = -t246 + t156;
t264 = t198 * qJDD(1);
t265 = t196 * qJDD(1);
t243 = -t203 * t265 + t205 * t264;
t98 = qJD(1) * t141 - t243;
t82 = -qJDD(4) * t197 - t195 * t98;
t83 = qJDD(4) * t195 - t197 * t98;
t17 = -qJD(6) * t222 + t113 * t268 + t202 * t82 - t322 * t83;
t330 = -t17 * t223 - t298 * t55;
t216 = -qJD(4) * t259 + qJDD(1) * t147;
t99 = qJD(4) * t258 + t216;
t96 = qJDD(6) + t99;
t329 = -t131 * t299 + t148 * t96;
t132 = t134 ^ 2;
t307 = t195 * t99;
t328 = -t132 * t197 - t307;
t327 = t141 * t337 + t149 * t82;
t81 = t128 * t205 + t129 * t203;
t70 = qJD(4) * qJ(5) + t81;
t177 = qJD(1) * qJ(2) + qJD(3);
t153 = pkin(3) * t271 + t177;
t75 = pkin(4) * t134 - qJ(5) * t137 + t153;
t35 = -t195 * t70 + t197 * t75;
t21 = pkin(5) * t134 - pkin(8) * t113 + t35;
t36 = t195 * t75 + t197 * t70;
t23 = pkin(8) * t337 + t36;
t225 = t202 * t23 - t21 * t322;
t263 = t118 * t205 + t119 * t203 + t129 * t269;
t32 = qJDD(4) * qJ(5) + (qJD(5) - t125) * qJD(4) + t263;
t174 = pkin(3) * t265;
t144 = t174 + t156;
t33 = pkin(4) * t99 + qJ(5) * t98 - qJD(5) * t137 + t144;
t12 = -t195 * t32 + t197 * t33;
t6 = pkin(5) * t99 - pkin(8) * t83 + t12;
t13 = t195 * t33 + t197 * t32;
t9 = -pkin(8) * t82 + t13;
t1 = -qJD(6) * t225 + t202 * t6 + t322 * t9;
t325 = t137 ^ 2;
t324 = 0.2e1 * t193;
t323 = t82 * pkin(5);
t321 = pkin(5) * t195;
t320 = pkin(8) * t197;
t317 = t55 * t56;
t315 = pkin(8) + qJ(5);
t142 = -t196 * t270 + t198 * t269;
t59 = pkin(4) * t142 + qJ(5) * t141 - qJD(5) * t149 + qJD(2);
t66 = -qJD(3) * t147 + qJD(4) * t336;
t27 = t195 * t59 + t197 * t66;
t95 = pkin(4) * t137 + qJ(5) * t134;
t41 = t195 * t95 + t197 * t80;
t154 = t315 * t195;
t155 = t315 * t197;
t106 = -t154 * t322 - t155 * t202;
t291 = t134 * t197;
t40 = -t195 * t80 + t197 * t95;
t24 = pkin(5) * t137 + pkin(8) * t291 + t40;
t292 = t134 * t195;
t30 = pkin(8) * t292 + t41;
t314 = qJD(5) * t223 + qJD(6) * t106 - t202 * t24 - t30 * t322;
t107 = -t154 * t202 + t155 * t322;
t313 = -qJD(5) * t148 - qJD(6) * t107 + t202 * t30 - t24 * t322;
t93 = t197 * t99;
t312 = -t132 * t195 + t93;
t311 = t137 * t56;
t306 = t55 * t137;
t305 = t82 * t197;
t304 = t83 * t195;
t303 = t83 * t197;
t302 = t99 * t147;
t103 = t151 * t205 + t152 * t203;
t170 = qJ(2) + t182;
t94 = pkin(4) * t147 - qJ(5) * t149 + t170;
t48 = t103 * t197 + t195 * t94;
t301 = qJD(1) * t223 + qJD(6) * t90 + t142 * t148;
t88 = t148 * t147;
t300 = qJD(1) * t148 + qJD(6) * t88 - t142 * t223;
t297 = pkin(1) * qJDD(1);
t295 = t113 * t137;
t294 = t113 * t195;
t293 = t134 * t137;
t290 = t141 * t195;
t289 = t149 * t195;
t288 = t181 * t204;
t287 = t181 * t206;
t190 = pkin(10) + qJ(6);
t178 = sin(t190);
t283 = t204 * t178;
t180 = cos(t190);
t282 = t204 * t180;
t280 = t206 * t178;
t279 = t206 * t180;
t278 = t81 * qJD(4);
t277 = -qJD(5) + t65;
t276 = -qJD(4) * t141 + qJDD(4) * t149;
t275 = g(1) * t287 + g(2) * t288;
t274 = pkin(1) * t206 + qJ(2) * t204;
t262 = t182 * t204 + t274;
t184 = t206 * qJ(2);
t257 = -t204 * pkin(1) + t184;
t26 = -t195 * t66 + t197 * t59;
t256 = t202 * t83 + t322 * t82;
t47 = -t103 * t195 + t197 * t94;
t254 = t272 * t150;
t251 = qJDD(2) - t297;
t18 = qJD(6) * t55 + t256;
t249 = -t148 * t18 - t299 * t56;
t247 = -t131 * t298 + t223 * t96;
t244 = pkin(4) * t179 - qJ(5) * t181;
t242 = -t12 * t197 - t13 * t195;
t241 = t12 * t147 + t35 * t142;
t240 = -t13 * t147 - t36 * t142;
t238 = -t195 * t35 + t197 * t36;
t237 = -t141 * t113 + t149 * t83;
t236 = t113 * t142 + t83 * t147;
t235 = -t134 * t141 + t149 * t99;
t234 = t134 * t142 + t302;
t233 = -t137 * t141 - t149 * t98;
t173 = pkin(5) * t197 + pkin(4);
t231 = t173 * t179 - t181 * t315;
t230 = t197 * t337;
t229 = t257 + t335;
t228 = -qJD(4) * t142 - qJDD(4) * t147;
t227 = -t206 * t200 + t262;
t31 = pkin(5) * t147 - t149 * t320 + t47;
t39 = -pkin(8) * t289 + t48;
t14 = -t202 * t39 + t31 * t322;
t8 = t202 * t21 + t23 * t322;
t15 = t202 * t31 + t322 * t39;
t218 = t142 * t337 - t82 * t147;
t217 = t230 - t294;
t211 = -g(1) * t288 + g(2) * t287 - t250 + t319;
t2 = -qJD(6) * t8 - t202 * t9 + t322 * t6;
t210 = t219 + t333;
t37 = -t128 * t270 + t263;
t209 = -t141 * t80 + t142 * t81 + t147 * t37 - t149 * t250 - t334;
t67 = qJD(3) * t149 + qJD(4) * t103;
t208 = -t211 + t331;
t207 = qJD(1) ^ 2;
t124 = t179 * t279 - t283;
t123 = t179 * t280 + t282;
t122 = t179 * t282 + t280;
t121 = -t179 * t283 + t279;
t91 = t223 * t149;
t89 = t148 * t149;
t71 = t195 * t82;
t68 = pkin(5) * t289 - t336;
t51 = -pkin(5) * t292 + t81;
t50 = -pkin(5) * t290 + t67;
t46 = -pkin(5) * t337 + t65;
t45 = -t141 * t148 + t149 * t332;
t43 = t140 * t149 + t141 * t223;
t22 = pkin(8) * t290 + t27;
t20 = pkin(5) * t142 + t141 * t320 + t26;
t19 = t34 + t323;
t4 = -qJD(6) * t15 + t20 * t322 - t202 * t22;
t3 = qJD(6) * t14 + t202 * t20 + t22 * t322;
t5 = [0, 0, 0, 0, 0, qJDD(1), t334, t246, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t334 - 0.2e1 * t297, 0.2e1 * t192 + t324 - t246, -t251 * pkin(1) - g(1) * t257 - g(2) * t274 + (t192 + t324) * qJ(2), t189 * qJDD(1), -0.2e1 * t196 * t264, 0, t188 * qJDD(1), 0, 0, t210 * t196, t210 * t198, t334 + t272 * (-t150 - t326) t156 * qJ(2) + t177 * qJD(2) - g(1) * (t201 * t204 + t184) - g(2) * (qJ(3) * t206 + t274) + t201 * t254 - qJD(3) * t338, t233, -t137 * t142 + t147 * t98 - t235, t276, t234, t228, 0, qJD(2) * t134 - t67 * qJD(4) + qJDD(4) * t336 + t153 * t142 + t144 * t147 + t170 * t99 - t224, qJD(2) * t137 - qJD(4) * t66 - qJDD(4) * t103 - t141 * t153 + t144 * t149 - t170 * t98 - t275, -t103 * t99 - t134 * t66 + t137 * t67 + t336 * t98 - t209, -g(1) * t229 - g(2) * t227 + qJD(2) * t153 + t103 * t37 + t144 * t170 - t250 * t336 + t66 * t81 - t67 * t80, t237 * t197 (-t304 - t305) * t149 - t217 * t141, t197 * t235 + t236, t327 * t195, -t195 * t235 + t218, t234, t26 * t134 + t195 * t213 - t197 * t224 - t336 * t82 - t337 * t67 + t47 * t99 + t241, t67 * t113 - t27 * t134 + t195 * t224 + t197 * t213 - t336 * t83 - t48 * t99 + t240, t27 * t337 - t48 * t82 - t26 * t113 - t47 * t83 + t242 * t149 + (t195 * t36 + t197 * t35) * t141 + t275, t13 * t48 + t36 * t27 + t12 * t47 + t35 * t26 - t34 * t336 + t65 * t67 - g(1) * (t206 * t244 + t229) - g(2) * (t204 * t244 + t227) -t17 * t91 - t43 * t55, t17 * t89 - t18 * t91 - t43 * t56 - t45 * t55, -t131 * t43 + t142 * t55 - t147 * t17 + t91 * t96, t18 * t89 - t45 * t56, -t131 * t45 + t142 * t56 - t147 * t18 - t89 * t96, t131 * t142 + t147 * t96, -g(1) * t124 - g(2) * t122 + t131 * t4 + t14 * t96 - t142 * t225 + t147 * t2 + t18 * t68 + t19 * t89 + t45 * t46 - t50 * t56, g(1) * t123 - g(2) * t121 - t1 * t147 - t131 * t3 - t142 * t8 - t15 * t96 - t17 * t68 + t19 * t91 - t43 * t46 + t50 * t55, -t1 * t89 + t14 * t17 - t15 * t18 - t2 * t91 - t225 * t43 + t3 * t56 - t4 * t55 - t45 * t8 + t275, t1 * t15 + t8 * t3 + t2 * t14 - t225 * t4 + t19 * t68 + t46 * t50 - g(1) * (t184 + t335) - g(2) * t262 + (-g(1) * t231 - g(2) * (-t200 + t321)) * t206 + (-g(1) * (-pkin(1) - t321) - g(2) * t231) * t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t207, -qJ(2) * t207 + t251 - t334, 0, 0, 0, 0, 0, 0, -t207 * t196, -t207 * t198, -t272 * qJDD(1), -t177 * qJD(1) + t254 - t334, 0, 0, 0, 0, 0, 0, -qJD(1) * t134 + t276, -qJD(1) * t137 + t228, -t233 - t234, -qJD(1) * t153 + t209, 0, 0, 0, 0, 0, 0, -t195 * t302 + (-qJD(1) * t197 - t142 * t195) * t134 - t327, -t197 * t302 + (qJD(1) * t195 - t142 * t197) * t134 - t237 (qJD(1) * t113 + t218) * t197 + (-qJD(1) * t337 + t236) * t195 (-qJD(1) * t35 - t240) * t197 + (-qJD(1) * t36 - t241) * t195 - t213, 0, 0, 0, 0, 0, 0, -t131 * t301 - t141 * t56 - t149 * t18 - t88 * t96, t131 * t300 + t141 * t55 + t149 * t17 - t90 * t96, -t88 * t17 - t90 * t18 - t300 * t56 + t301 * t55, t1 * t90 + t46 * t141 - t19 * t149 - t2 * t88 + t225 * t301 - t300 * t8 - t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, t264, -t272 * t207, qJD(1) * t338 + t219, 0, 0, 0, 0, 0, 0 (t137 + t258) * qJD(4) + t216, -0.2e1 * t134 * qJD(4) + t243, -t132 - t325, t134 * t81 + t137 * t80 + t174 + t219, 0, 0, 0, 0, 0, 0, t312 + t339, -t295 + t328, -t303 - t71 + (t230 + t294) * t134, t134 * t238 - t65 * t137 - t242 - t246, 0, 0, 0, 0, 0, 0, t247 + t311, -t306 - t329, t249 - t330, t1 * t148 - t46 * t137 + t2 * t223 + t225 * t298 - t299 * t8 - t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, -t132 + t325, t243, -t293 (t137 - t258) * qJD(4) - t216, qJDD(4), -t137 * t153 + t211 + t278, t153 * t134 + (t80 + t125) * qJD(4) + t221 - t263, 0, 0, t113 * t291 + t304, t134 * t217 + t303 - t71, -t295 - t328, -t292 * t337 - t305, t312 - t339, -t293, -qJ(5) * t307 - pkin(4) * t82 - t81 * t127 - t35 * t137 + (t195 * t277 - t40) * t134 + (-t342 + t278) * t197, -qJ(5) * t93 - pkin(4) * t83 - t81 * t113 + t36 * t137 + (t197 * t277 + t41) * t134 + t342 * t195, t40 * t113 + t41 * t127 + (-qJ(5) * t82 - qJD(5) * t127 - t35 * t134 + t13 + (qJD(5) * t197 - t41) * qJD(4)) * t197 + (qJ(5) * t83 + qJD(5) * t113 - t134 * t36 - t12) * t195 - t221, -t35 * t40 - t36 * t41 - t65 * t81 + t238 * qJD(5) - t342 * pkin(4) + (-t12 * t195 + t13 * t197 - t221) * qJ(5), -t17 * t148 - t299 * t55, t249 + t330, -t306 + t329, -t18 * t223 - t298 * t56, t247 - t311, -t131 * t137, t106 * t96 + t131 * t313 + t137 * t225 - t173 * t18 - t180 * t215 - t19 * t223 + t298 * t46 + t51 * t56, -t107 * t96 - t131 * t314 + t8 * t137 + t19 * t148 + t173 * t17 + t178 * t215 - t299 * t46 - t51 * t55, t1 * t223 + t106 * t17 - t107 * t18 - t2 * t148 - t225 * t299 - t298 * t8 - t313 * t55 + t314 * t56 - t221, g(3) * t231 + t1 * t107 + t2 * t106 - t19 * t173 - t313 * t225 + t314 * t8 - t46 * t51 - t334 * (t173 * t181 + t179 * t315); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113 * t134 + t82, t134 * t337 + t83, -t113 ^ 2 - t337 ^ 2, t113 * t35 - t337 * t36 + t208, 0, 0, 0, 0, 0, 0, t55 * t131 + t18, -t17 + t344, -t343 - t345, -t225 * t55 - t56 * t8 + t208 + t323; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t317, t343 - t345, -t17 - t344, t317, -t256 + (-qJD(6) + t131) * t55, t96, -g(1) * t121 - g(2) * t123 + t8 * t131 + t178 * t318 - t46 * t55 + t2, g(1) * t122 - g(2) * t124 - t131 * t225 + t180 * t318 - t46 * t56 - t1, 0, 0;];
tau_reg  = t5;
