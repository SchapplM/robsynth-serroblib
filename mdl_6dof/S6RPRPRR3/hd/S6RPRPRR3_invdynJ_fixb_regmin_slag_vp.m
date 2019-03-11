% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:30
% EndTime: 2019-03-09 03:42:42
% DurationCPUTime: 5.26s
% Computational Cost: add. (5388->460), mult. (11700->630), div. (0->0), fcn. (8729->18), ass. (0->232)
t223 = cos(qJ(3));
t219 = sin(qJ(3));
t210 = qJ(1) + pkin(10);
t203 = sin(t210);
t205 = cos(t210);
t263 = g(1) * t205 + g(2) * t203;
t248 = t263 * t219;
t236 = -g(3) * t223 + t248;
t214 = sin(pkin(10));
t196 = pkin(1) * t214 + pkin(7);
t184 = t196 * qJD(1);
t297 = qJD(3) * t223;
t182 = t196 * qJDD(1);
t334 = -qJD(2) * qJD(3) - t182;
t285 = -t184 * t297 + t334 * t219;
t240 = -qJDD(3) * pkin(3) + qJDD(4) - t285;
t79 = -qJDD(2) * t223 + t240;
t341 = t79 - t236;
t215 = cos(pkin(11));
t213 = sin(pkin(11));
t301 = qJD(1) * t223;
t283 = t213 * t301;
t143 = qJD(2) * t223 - t219 * t184;
t260 = pkin(3) * t219 - qJ(4) * t223;
t173 = t260 * qJD(1);
t92 = t215 * t143 + t213 * t173;
t69 = -pkin(8) * t283 + t92;
t340 = qJD(4) * t215 - t69;
t192 = -qJD(5) + t301;
t189 = -qJD(6) + t192;
t290 = t215 * qJD(3);
t302 = qJD(1) * t219;
t159 = t213 * t302 - t290;
t299 = qJD(3) * t213;
t161 = t215 * t302 + t299;
t218 = sin(qJ(5));
t222 = cos(qJ(5));
t100 = t222 * t159 + t161 * t218;
t217 = sin(qJ(6));
t221 = cos(qJ(6));
t99 = t159 * t218 - t161 * t222;
t39 = t221 * t100 - t217 * t99;
t339 = t189 * t39;
t338 = t192 * t99;
t254 = t100 * t217 + t221 * t99;
t337 = t254 * t39;
t336 = t100 * t192;
t335 = t189 * t254;
t312 = t213 * t218;
t167 = -t222 * t215 + t312;
t242 = t167 * t223;
t306 = qJD(1) * t242 - t167 * qJD(5);
t168 = t213 * t222 + t215 * t218;
t243 = t168 * t223;
t305 = -qJD(1) * t243 + t168 * qJD(5);
t333 = t254 ^ 2 - t39 ^ 2;
t291 = qJD(6) * t221;
t292 = qJD(6) * t217;
t200 = t215 * qJDD(3);
t289 = qJD(1) * qJD(3);
t277 = t223 * t289;
t286 = t219 * qJDD(1);
t241 = t277 + t286;
t125 = t213 * t241 - t200;
t287 = qJDD(3) * t213;
t126 = t215 * t241 + t287;
t293 = qJD(5) * t222;
t295 = qJD(5) * t218;
t34 = -t218 * t125 + t222 * t126 - t159 * t293 - t161 * t295;
t35 = -qJD(5) * t99 + t222 * t125 + t126 * t218;
t6 = -t100 * t291 - t217 * t35 + t221 * t34 + t292 * t99;
t332 = t6 - t339;
t144 = t219 * qJD(2) + t223 * t184;
t130 = qJD(3) * qJ(4) + t144;
t251 = pkin(3) * t223 + qJ(4) * t219 + pkin(2);
t216 = cos(pkin(10));
t325 = pkin(1) * t216;
t155 = -t251 - t325;
t133 = t155 * qJD(1);
t59 = -t130 * t213 + t215 * t133;
t47 = -pkin(4) * t301 - pkin(8) * t161 + t59;
t60 = t215 * t130 + t213 * t133;
t49 = -pkin(8) * t159 + t60;
t17 = t218 * t47 + t222 * t49;
t13 = -pkin(9) * t100 + t17;
t11 = t13 * t292;
t209 = pkin(11) + qJ(5);
t207 = qJ(6) + t209;
t194 = sin(t207);
t195 = cos(t207);
t314 = t203 * t223;
t111 = t194 * t205 - t195 * t314;
t313 = t205 * t223;
t113 = t194 * t203 + t195 * t313;
t324 = g(3) * t219;
t127 = -qJD(3) * pkin(3) + qJD(4) - t143;
t93 = pkin(4) * t159 + t127;
t46 = pkin(5) * t100 + t93;
t331 = g(1) * t113 - g(2) * t111 + t195 * t324 + t39 * t46 + t11;
t110 = t194 * t314 + t195 * t205;
t112 = -t194 * t313 + t195 * t203;
t206 = t223 * qJDD(1);
t326 = t219 * t289 - t206;
t166 = qJDD(5) + t326;
t71 = qJDD(3) * qJ(4) + qJDD(2) * t219 + t182 * t223 + (qJD(4) + t143) * qJD(3);
t147 = qJD(3) * t260 - qJD(4) * t219;
t84 = qJD(1) * t147 + qJDD(1) * t155;
t29 = -t213 * t71 + t215 * t84;
t23 = t326 * pkin(4) - pkin(8) * t126 + t29;
t30 = t213 * t84 + t215 * t71;
t26 = -pkin(8) * t125 + t30;
t273 = -t218 * t26 + t222 * t23;
t231 = -qJD(5) * t17 + t273;
t2 = pkin(5) * t166 - pkin(9) * t34 + t231;
t246 = t218 * t23 + t222 * t26 + t47 * t293 - t295 * t49;
t3 = -pkin(9) * t35 + t246;
t284 = t221 * t2 - t217 * t3;
t16 = -t218 * t49 + t222 * t47;
t12 = pkin(9) * t99 + t16;
t10 = -pkin(5) * t192 + t12;
t316 = t13 * t221;
t5 = t217 * t10 + t316;
t330 = -g(1) * t112 + g(2) * t110 - qJD(6) * t5 + t194 * t324 + t46 * t254 + t284;
t230 = qJD(6) * t254 - t217 * t34 - t221 * t35;
t329 = t230 + t335;
t307 = qJDD(2) - g(3);
t328 = t307 * t223;
t322 = pkin(8) + qJ(4);
t180 = t322 * t213;
t181 = t322 * t215;
t304 = -t218 * t180 + t222 * t181;
t253 = qJD(4) * t213 + qJD(5) * t181;
t308 = t215 * t223;
t252 = pkin(4) * t219 - pkin(8) * t308;
t91 = -t143 * t213 + t215 * t173;
t56 = qJD(1) * t252 + t91;
t327 = -t180 * t293 + t340 * t222 + (-t253 - t56) * t218;
t154 = qJDD(6) + t166;
t139 = t168 * t219;
t140 = t167 * t219;
t78 = -t139 * t217 - t140 * t221;
t294 = qJD(5) * t219;
t86 = -qJD(3) * t242 - t168 * t294;
t309 = t215 * t219;
t87 = qJD(3) * t243 + t293 * t309 - t294 * t312;
t19 = qJD(6) * t78 + t217 * t86 + t221 * t87;
t77 = t221 * t139 - t140 * t217;
t321 = -t77 * t154 + t19 * t189;
t104 = t221 * t167 + t168 * t217;
t320 = -qJD(6) * t104 - t305 * t217 + t306 * t221;
t105 = -t167 * t217 + t168 * t221;
t319 = qJD(6) * t105 + t306 * t217 + t305 * t221;
t142 = t215 * t155;
t81 = -pkin(8) * t309 + t142 + (-t196 * t213 - pkin(4)) * t223;
t108 = t213 * t155 + t196 * t308;
t311 = t213 * t219;
t90 = -pkin(8) * t311 + t108;
t317 = t218 * t81 + t222 * t90;
t315 = -t139 * t166 + t87 * t192;
t310 = t213 * t223;
t186 = t219 * t196;
t298 = qJD(3) * t219;
t282 = t196 * t298;
t97 = t215 * t147 + t213 * t282;
t137 = (pkin(4) * t213 + t196) * t297;
t146 = pkin(4) * t311 + t186;
t211 = t219 ^ 2;
t303 = -t223 ^ 2 + t211;
t197 = -pkin(2) - t325;
t185 = qJD(1) * t197;
t109 = pkin(4) * t283 + t144;
t198 = -pkin(4) * t215 - pkin(3);
t281 = t305 * pkin(5) - t109;
t280 = -t223 * t6 - t254 * t298;
t279 = qJ(4) * t206;
t276 = t213 * t286;
t275 = t215 * t286;
t274 = qJD(6) * t10 + t3;
t67 = qJD(3) * t252 + t97;
t135 = t213 * t147;
t80 = t135 + (-pkin(8) * t310 - t215 * t186) * qJD(3);
t271 = -t218 * t80 + t222 * t67;
t270 = -t218 * t90 + t222 * t81;
t269 = -t223 * t34 - t298 * t99;
t267 = -t222 * t180 - t181 * t218;
t266 = qJD(1) * t303;
t55 = t222 * t56;
t83 = -pkin(9) * t167 + t304;
t265 = pkin(5) * t302 + t306 * pkin(9) + t168 * qJD(4) + t304 * qJD(5) + qJD(6) * t83 - t218 * t69 + t55;
t82 = -pkin(9) * t168 + t267;
t264 = -t305 * pkin(9) + qJD(6) * t82 + t327;
t262 = g(1) * t203 - g(2) * t205;
t220 = sin(qJ(1));
t224 = cos(qJ(1));
t261 = g(1) * t220 - g(2) * t224;
t18 = -qJD(6) * t77 - t217 * t87 + t221 * t86;
t259 = -t154 * t78 + t18 * t189;
t258 = -t29 * t213 + t30 * t215;
t257 = -t213 * t59 + t215 * t60;
t255 = t140 * t166 + t192 * t86;
t250 = t261 * pkin(1);
t249 = -t223 * t230 - t298 * t39;
t247 = -t100 * t298 + t223 * t35;
t245 = t218 * t67 + t222 * t80 + t81 * t293 - t295 * t90;
t239 = -qJD(1) * t185 + t263;
t238 = -qJ(4) * t298 + (qJD(4) - t127) * t223;
t237 = 0.2e1 * qJD(3) * t185 - qJDD(3) * t196;
t234 = -t223 * t263 - t324;
t225 = qJD(3) ^ 2;
t229 = -0.2e1 * qJDD(1) * t197 - t196 * t225 + t262;
t48 = pkin(4) * t125 + t79;
t226 = qJD(1) ^ 2;
t204 = cos(t209);
t202 = sin(t209);
t176 = qJDD(3) * t223 - t219 * t225;
t175 = qJDD(3) * t219 + t223 * t225;
t134 = pkin(5) * t167 + t198;
t124 = t202 * t203 + t204 * t313;
t123 = -t202 * t313 + t203 * t204;
t122 = t202 * t205 - t204 * t314;
t121 = t202 * t314 + t204 * t205;
t107 = -t196 * t310 + t142;
t98 = -t215 * t282 + t135;
t96 = pkin(5) * t139 + t146;
t50 = pkin(5) * t87 + t137;
t28 = -pkin(9) * t139 + t317;
t27 = -pkin(5) * t223 + pkin(9) * t140 + t270;
t14 = pkin(5) * t35 + t48;
t9 = -pkin(9) * t87 + t245;
t8 = pkin(5) * t298 - pkin(9) * t86 - qJD(5) * t317 + t271;
t4 = t221 * t10 - t13 * t217;
t1 = [qJDD(1), t261, g(1) * t224 + g(2) * t220 (t214 ^ 2 + t216 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t250, qJDD(1) * t211 + 0.2e1 * t219 * t277, -0.2e1 * qJD(3) * t266 + 0.2e1 * t206 * t219, t175, t176, 0, t219 * t237 + t223 * t229, -t219 * t229 + t223 * t237, -t263 * t213 + (t196 * t125 + t79 * t213 + (qJD(1) * t107 + t59) * qJD(3)) * t219 + (-t97 * qJD(1) - t107 * qJDD(1) - t29 + t262 * t215 + (t127 * t213 + t159 * t196) * qJD(3)) * t223, -t263 * t215 + (t126 * t196 + t215 * t79 + (-qJD(1) * t108 - t60) * qJD(3)) * t219 + (qJD(1) * t98 + qJDD(1) * t108 + t30 - t262 * t213 + (t127 * t215 + t161 * t196) * qJD(3)) * t223, -t107 * t126 - t108 * t125 - t159 * t98 - t161 * t97 + (-t213 * t60 - t215 * t59) * t297 + (-t213 * t30 - t215 * t29 + t262) * t219, t29 * t107 + t30 * t108 + t59 * t97 + t60 * t98 + t250 + (-g(1) * pkin(7) - g(2) * t251) * t205 + (-g(2) * pkin(7) + g(1) * t251) * t203 + (t127 * t297 + t219 * t79) * t196, -t140 * t34 - t86 * t99, -t100 * t86 - t139 * t34 + t140 * t35 + t87 * t99, -t255 + t269, t247 + t315, -t166 * t223 - t192 * t298, -t271 * t192 + t270 * t166 - t273 * t223 + t16 * t298 + t137 * t100 + t146 * t35 + t48 * t139 + t93 * t87 - g(1) * t122 - g(2) * t124 + (t17 * t223 + t192 * t317) * qJD(5), -g(1) * t121 - g(2) * t123 - t137 * t99 - t48 * t140 + t146 * t34 - t317 * t166 - t17 * t298 + t245 * t192 + t246 * t223 + t93 * t86, -t18 * t254 + t6 * t78, -t18 * t39 + t19 * t254 + t230 * t78 - t6 * t77, -t259 + t280, t249 + t321, -t154 * t223 - t189 * t298 -(-t217 * t9 + t221 * t8) * t189 + (-t217 * t28 + t221 * t27) * t154 - t284 * t223 + t4 * t298 + t50 * t39 - t96 * t230 + t14 * t77 + t46 * t19 - g(1) * t111 - g(2) * t113 + (-(-t217 * t27 - t221 * t28) * t189 + t5 * t223) * qJD(6), -t5 * t298 - g(1) * t110 - g(2) * t112 - t11 * t223 + t14 * t78 + t46 * t18 - t50 * t254 + t96 * t6 + ((-qJD(6) * t28 + t8) * t189 - t27 * t154 + t2 * t223) * t217 + ((qJD(6) * t27 + t9) * t189 - t28 * t154 + t274 * t223) * t221; 0, 0, 0, t307, 0, 0, 0, 0, 0, t176, -t175 (-t125 + t276) * t223 + (t159 * t219 - t213 * t266) * qJD(3) (-t126 + t275) * t223 + (t161 * t219 - t215 * t266) * qJD(3) (-t125 * t215 + t126 * t213) * t219 + (-t159 * t215 + t161 * t213) * t297, -t223 * t79 - g(3) + t258 * t219 + (t127 * t219 + t223 * t257) * qJD(3), 0, 0, 0, 0, 0, -t247 + t315, t255 + t269, 0, 0, 0, 0, 0, -t249 + t321, t259 + t280; 0, 0, 0, 0, -t219 * t226 * t223, t303 * t226, t286, t206, qJDD(3), qJD(3) * t144 + t239 * t219 + t285 + t328, qJD(3) * t143 + (qJD(3) * t184 - t307) * t219 + (t239 + t334) * t223, t213 * t279 - pkin(3) * t125 - t144 * t159 - t341 * t215 + (t213 * t238 - t219 * t59 + t223 * t91) * qJD(1), t215 * t279 - pkin(3) * t126 - t144 * t161 + t341 * t213 + (t215 * t238 + t219 * t60 - t223 * t92) * qJD(1), t159 * t92 + t161 * t91 + (-qJ(4) * t125 - qJD(4) * t159 + t301 * t59 + t30) * t215 + (qJ(4) * t126 + qJD(4) * t161 + t301 * t60 - t29) * t213 + t234, -t127 * t144 - t59 * t91 - t60 * t92 + t257 * qJD(4) - t341 * pkin(3) + (t234 + t258) * qJ(4), t168 * t34 - t306 * t99, -t306 * t100 - t167 * t34 - t168 * t35 + t305 * t99, t166 * t168 - t306 * t192 + t302 * t99, t100 * t302 - t166 * t167 + t305 * t192, t192 * t302, t267 * t166 + t198 * t35 + t48 * t167 - t16 * t302 - t109 * t100 + t305 * t93 + (t55 + t253 * t222 + (-qJD(5) * t180 + t340) * t218) * t192 + t236 * t204, t109 * t99 - t304 * t166 + t48 * t168 + t17 * t302 + t327 * t192 + t198 * t34 - t202 * t236 + t306 * t93, t105 * t6 - t254 * t320, -t104 * t6 + t105 * t230 + t254 * t319 - t320 * t39, t105 * t154 - t320 * t189 + t254 * t302, -t104 * t154 + t319 * t189 + t39 * t302, t189 * t302 (-t217 * t83 + t221 * t82) * t154 - t134 * t230 + t14 * t104 - t4 * t302 + t319 * t46 + t281 * t39 + (t217 * t264 + t221 * t265) * t189 + t236 * t195 -(t217 * t82 + t221 * t83) * t154 + t134 * t6 + t14 * t105 + t5 * t302 + t320 * t46 - t281 * t254 + (-t217 * t265 + t221 * t264) * t189 - t236 * t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276 - t200 + (-t161 + t299) * t301, t275 + t287 + (t159 + t290) * t301, -t159 ^ 2 - t161 ^ 2, t159 * t60 + t161 * t59 + t240 - t248 - t328, 0, 0, 0, 0, 0, t35 + t338, t34 + t336, 0, 0, 0, 0, 0, -t230 + t335, t6 + t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99 * t100, -t100 ^ 2 + t99 ^ 2, t34 - t336, -t35 + t338, t166, -g(1) * t123 + g(2) * t121 - t17 * t192 + t202 * t324 + t93 * t99 + t231, g(1) * t124 - g(2) * t122 + t100 * t93 - t16 * t192 + t204 * t324 - t246, -t337, t333, t332, t329, t154 (-t12 * t217 - t316) * t189 + (t221 * t154 + t189 * t292 + t39 * t99) * pkin(5) + t330 (t13 * t189 - t2) * t217 + (-t12 * t189 - t274) * t221 + (-t217 * t154 + t189 * t291 - t254 * t99) * pkin(5) + t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t337, t333, t332, t329, t154, -t189 * t5 + t330, -t189 * t4 - t217 * t2 - t221 * t274 + t331;];
tau_reg  = t1;
