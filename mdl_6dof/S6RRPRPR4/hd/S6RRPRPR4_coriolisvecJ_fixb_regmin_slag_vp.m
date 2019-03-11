% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:17
% EndTime: 2019-03-09 10:26:30
% DurationCPUTime: 4.78s
% Computational Cost: add. (10668->422), mult. (32202->594), div. (0->0), fcn. (26382->12), ass. (0->234)
t214 = cos(pkin(11));
t221 = cos(qJ(2));
t212 = sin(pkin(6));
t287 = qJD(1) * t212;
t272 = t221 * t287;
t193 = t214 * t272;
t211 = sin(pkin(11));
t218 = sin(qJ(2));
t273 = t218 * t287;
t177 = -t211 * t273 + t193;
t173 = qJD(4) - t177;
t217 = sin(qJ(4));
t220 = cos(qJ(4));
t204 = pkin(2) * t211 + pkin(9);
t294 = qJ(5) + t204;
t260 = qJD(4) * t294;
t215 = cos(pkin(6));
t321 = pkin(1) * t218;
t278 = t215 * t321;
t299 = t212 * t221;
t317 = pkin(8) + qJ(3);
t174 = t317 * t299 + t278;
t163 = t174 * qJD(1);
t151 = t211 * t163;
t320 = pkin(1) * t221;
t277 = t215 * t320;
t198 = qJD(1) * t277;
t270 = t317 * t218;
t256 = t212 * t270;
t162 = -qJD(1) * t256 + t198;
t111 = t162 * t214 - t151;
t243 = t211 * t221 + t214 * t218;
t180 = t243 * t287;
t126 = pkin(2) * t273 + pkin(3) * t180 - pkin(9) * t177;
t291 = t220 * t111 + t217 * t126;
t303 = t177 * t217;
t335 = qJ(5) * t303 + t220 * qJD(5) - t217 * t260 - t291;
t115 = t220 * t126;
t334 = -pkin(4) * t180 - t115 + (qJ(5) * t177 - t260) * t220 + (-qJD(5) + t111) * t217;
t284 = qJD(4) * t217;
t333 = t284 - t303;
t286 = qJD(1) * t215;
t200 = qJD(2) + t286;
t143 = t180 * t217 - t220 * t200;
t210 = sin(pkin(12));
t213 = cos(pkin(12));
t244 = -t180 * t220 - t200 * t217;
t261 = -t213 * t143 + t210 * t244;
t326 = qJD(6) - t261;
t216 = sin(qJ(6));
t219 = cos(qJ(6));
t245 = -t143 * t210 - t213 * t244;
t73 = -t219 * t173 + t216 * t245;
t332 = t326 * t73;
t189 = t210 * t217 - t213 * t220;
t123 = t189 * t177;
t187 = t189 * qJD(4);
t331 = t123 - t187;
t190 = t210 * t220 + t213 * t217;
t289 = t173 * t190;
t265 = t219 * t326;
t285 = qJD(2) * t212;
t271 = t218 * t285;
t254 = qJD(1) * t271;
t172 = qJD(2) * t193 - t211 * t254;
t283 = qJD(4) * t220;
t97 = t220 * t172 - t180 * t284 + t200 * t283;
t98 = -qJD(4) * t244 + t217 * t172;
t67 = t210 * t97 + t213 * t98;
t312 = t216 * t67;
t330 = -t265 * t326 - t312;
t329 = pkin(2) * t271;
t309 = t335 * t210 - t334 * t213;
t308 = t334 * t210 + t335 * t213;
t161 = (pkin(2) + t320) * t215 - t256;
t121 = t211 * t161 + t214 * t174;
t109 = pkin(9) * t215 + t121;
t297 = t214 * t221;
t300 = t212 * t218;
t182 = t211 * t300 - t212 * t297;
t183 = t243 * t212;
t253 = (-pkin(2) * t221 - pkin(1)) * t212;
t128 = t182 * pkin(3) - t183 * pkin(9) + t253;
t292 = t220 * t109 + t217 * t128;
t298 = t214 * t163;
t110 = t162 * t211 + t298;
t327 = t333 * pkin(4) - t110;
t203 = pkin(4) * t210 + pkin(10);
t178 = qJD(2) * t183;
t171 = qJD(1) * t178;
t195 = pkin(2) * t254;
t107 = pkin(3) * t171 - pkin(9) * t172 + t195;
t196 = qJD(2) * t198;
t226 = (-qJD(2) * t270 + qJD(3) * t221) * t212;
t137 = qJD(1) * t226 + t196;
t148 = -t174 * qJD(2) - qJD(3) * t300;
t138 = t148 * qJD(1);
t80 = t137 * t214 + t138 * t211;
t268 = t220 * t107 - t217 * t80;
t240 = qJD(1) * t253;
t185 = qJD(3) + t240;
t106 = -t177 * pkin(3) - t180 * pkin(9) + t185;
t146 = t200 * pkin(2) + t162;
t95 = t211 * t146 + t298;
t90 = pkin(9) * t200 + t95;
t64 = t106 * t217 + t220 * t90;
t225 = -qJD(4) * t64 + t268;
t16 = t171 * pkin(4) - t97 * qJ(5) + qJD(5) * t244 + t225;
t237 = t106 * t283 + t217 * t107 + t220 * t80 - t284 * t90;
t20 = -qJ(5) * t98 - qJD(5) * t143 + t237;
t5 = t16 * t213 - t20 * t210;
t3 = -pkin(5) * t171 - t5;
t325 = (-pkin(4) * t244 + pkin(5) * t245 - pkin(10) * t261 + qJD(6) * t203) * t326 + t3;
t282 = qJD(6) * t216;
t82 = -t123 * t219 + t180 * t216;
t231 = t187 * t219 + t190 * t282 + t82;
t302 = t190 * t219;
t324 = -t231 * t326 + t67 * t302;
t206 = -pkin(2) * t214 - pkin(3);
t242 = -pkin(4) * t220 + t206;
t136 = pkin(5) * t189 - pkin(10) * t190 + t242;
t188 = t294 * t220;
t264 = t294 * t217;
t142 = t213 * t188 - t210 * t264;
t323 = (-t289 * pkin(5) + t331 * pkin(10) + qJD(6) * t142 - t327) * t326 - t136 * t67;
t68 = -t210 * t98 + t213 * t97;
t266 = -t219 * t171 + t216 * t68;
t75 = t173 * t216 + t219 * t245;
t40 = qJD(6) * t75 + t266;
t322 = -t189 * t40 - t289 * t73;
t79 = t137 * t211 - t214 * t138;
t59 = pkin(4) * t98 + t79;
t15 = pkin(5) * t67 - pkin(10) * t68 + t59;
t63 = t220 * t106 - t217 * t90;
t51 = qJ(5) * t244 + t63;
t46 = pkin(4) * t173 + t51;
t52 = -qJ(5) * t143 + t64;
t49 = t213 * t52;
t22 = t210 * t46 + t49;
t18 = pkin(10) * t173 + t22;
t94 = t146 * t214 - t151;
t89 = -pkin(3) * t200 - t94;
t72 = pkin(4) * t143 + qJD(5) + t89;
t41 = -pkin(5) * t261 - pkin(10) * t245 + t72;
t249 = t18 * t216 - t219 * t41;
t6 = t210 * t16 + t213 * t20;
t4 = pkin(10) * t171 + t6;
t1 = -qJD(6) * t249 + t216 * t15 + t219 * t4;
t319 = t73 * t245;
t318 = t75 * t245;
t154 = t183 * t217 - t215 * t220;
t179 = (-t211 * t218 + t297) * t285;
t119 = -qJD(4) * t154 + t179 * t220;
t155 = t183 * t220 + t215 * t217;
t127 = pkin(3) * t178 - pkin(9) * t179 + t329;
t199 = qJD(2) * t277;
t147 = t199 + t226;
t92 = t147 * t214 + t148 * t211;
t267 = t220 * t127 - t217 * t92;
t27 = t178 * pkin(4) - t119 * qJ(5) - qJD(4) * t292 - t155 * qJD(5) + t267;
t118 = qJD(4) * t155 + t179 * t217;
t234 = -t109 * t284 + t217 * t127 + t128 * t283 + t220 * t92;
t33 = -qJ(5) * t118 - qJD(5) * t154 + t234;
t10 = t210 * t27 + t213 * t33;
t262 = -t109 * t217 + t220 * t128;
t55 = pkin(4) * t182 - qJ(5) * t155 + t262;
t61 = -qJ(5) * t154 + t292;
t35 = t210 * t55 + t213 * t61;
t315 = t142 * t67;
t313 = t210 * t52;
t281 = qJD(6) * t219;
t39 = t216 * t171 + t173 * t281 + t219 * t68 - t245 * t282;
t311 = t39 * t216;
t310 = pkin(5) * t180 + t309;
t307 = t143 * t173;
t306 = t143 * t180;
t305 = t244 * t173;
t304 = t244 * t180;
t207 = t212 ^ 2;
t222 = qJD(1) ^ 2;
t301 = t207 * t222;
t296 = t217 * t171;
t288 = t218 ^ 2 - t221 ^ 2;
t280 = qJD(2) - t200;
t279 = t207 * t321;
t274 = t221 * t301;
t269 = qJD(1) * qJD(2) * t207;
t91 = t147 * t211 - t214 * t148;
t120 = t161 * t214 - t211 * t174;
t259 = t173 * t220;
t258 = t200 + t286;
t257 = t39 * t189 + t289 * t75;
t255 = t221 * t269;
t12 = t18 * t219 + t216 * t41;
t9 = -t210 * t33 + t213 * t27;
t21 = t213 * t46 - t313;
t34 = -t210 * t61 + t213 * t55;
t30 = pkin(10) * t182 + t35;
t100 = -t154 * t210 + t155 * t213;
t108 = -pkin(3) * t215 - t120;
t228 = pkin(4) * t154 + t108;
t99 = t213 * t154 + t155 * t210;
t44 = pkin(5) * t99 - pkin(10) * t100 + t228;
t248 = t216 * t44 + t219 * t30;
t247 = -t216 * t30 + t219 * t44;
t77 = t100 * t219 + t182 * t216;
t76 = t100 * t216 - t182 * t219;
t241 = t219 * t67 + (t216 * t261 - t282) * t326;
t239 = t220 * t171 - t333 * t173;
t238 = pkin(4) * t118 + t91;
t236 = -pkin(8) * t299 - t278;
t235 = -pkin(8) * t254 + t196;
t233 = -t204 * t171 + t173 * t89;
t81 = -t123 * t216 - t219 * t180;
t232 = -t187 * t216 + t190 * t281 - t81;
t229 = t236 * t200;
t17 = -pkin(5) * t173 - t21;
t25 = t213 * t51 - t313;
t227 = -t203 * t67 + (t17 + t25) * t326;
t2 = -qJD(6) * t12 + t219 * t15 - t216 * t4;
t224 = -t315 + t3 * t190 + (pkin(10) * t180 - qJD(6) * t136 - t308) * t326;
t223 = -t190 * t312 - t232 * t326;
t205 = -pkin(4) * t213 - pkin(5);
t141 = t210 * t188 + t213 * t264;
t71 = -t118 * t210 + t119 * t213;
t70 = t213 * t118 + t119 * t210;
t43 = qJD(6) * t77 - t178 * t219 + t216 * t71;
t42 = -qJD(6) * t76 + t178 * t216 + t219 * t71;
t29 = -pkin(5) * t182 - t34;
t24 = t210 * t51 + t49;
t23 = pkin(5) * t70 - pkin(10) * t71 + t238;
t8 = pkin(10) * t178 + t10;
t7 = -pkin(5) * t178 - t9;
t11 = [0, 0, 0, 0.2e1 * t218 * t255, -0.2e1 * t288 * t269, t258 * t221 * t285, -t258 * t271, 0 (t229 + (t215 * t236 - 0.2e1 * t279) * qJD(1)) * qJD(2), -0.2e1 * pkin(1) * t255 - (-pkin(8) * t271 + t199) * t200 - t235 * t215, -t120 * t172 - t121 * t171 + t177 * t92 - t178 * t95 - t179 * t94 + t180 * t91 - t182 * t80 + t183 * t79, -t79 * t120 + t80 * t121 - t94 * t91 + t95 * t92 + (t185 + t240) * t329, -t119 * t244 + t155 * t97, t118 * t244 - t119 * t143 - t154 * t97 - t155 * t98, t119 * t173 + t155 * t171 - t178 * t244 + t182 * t97, -t118 * t173 - t143 * t178 - t154 * t171 - t182 * t98, t171 * t182 + t173 * t178, t267 * t173 + t262 * t171 + t268 * t182 + t63 * t178 + t91 * t143 + t108 * t98 + t79 * t154 + t89 * t118 + (-t173 * t292 - t182 * t64) * qJD(4), t108 * t97 + t89 * t119 + t79 * t155 - t171 * t292 - t173 * t234 - t64 * t178 - t182 * t237 - t244 * t91, t10 * t261 - t100 * t5 - t21 * t71 - t22 * t70 - t245 * t9 - t34 * t68 - t35 * t67 - t6 * t99, t22 * t10 + t21 * t9 + t228 * t59 + t238 * t72 + t5 * t34 + t6 * t35, t39 * t77 + t42 * t75, -t39 * t76 - t40 * t77 - t42 * t73 - t43 * t75, t326 * t42 + t39 * t99 + t67 * t77 + t70 * t75, -t326 * t43 - t40 * t99 - t67 * t76 - t70 * t73, t326 * t70 + t67 * t99 (-qJD(6) * t248 - t216 * t8 + t219 * t23) * t326 + t247 * t67 + t2 * t99 - t249 * t70 + t7 * t73 + t29 * t40 + t3 * t76 + t17 * t43 -(qJD(6) * t247 + t216 * t23 + t219 * t8) * t326 - t248 * t67 - t1 * t99 - t12 * t70 + t7 * t75 + t29 * t39 + t3 * t77 + t17 * t42; 0, 0, 0, -t218 * t274, t288 * t301, t280 * t272, -t280 * t273, 0, t222 * t279 + (qJD(2) * t236 - t229) * qJD(1), pkin(1) * t274 + (-pkin(8) * t273 + t198) * t200 - t235 (-t110 + t95) * t180 + (-t111 + t94) * t177 + (-t171 * t211 - t172 * t214) * pkin(2), t94 * t110 - t95 * t111 + (-t185 * t273 + t211 * t80 - t214 * t79) * pkin(2), t97 * t217 - t244 * t259 (t97 - t307) * t220 + (-t98 + t305) * t217, t173 * t259 + t296 + t304, t239 + t306, -t173 * t180, -t110 * t143 - t63 * t180 + t206 * t98 - t79 * t220 + (-t204 * t283 - t115) * t173 + (t111 * t173 + t233) * t217, t110 * t244 + t64 * t180 + t206 * t97 + t79 * t217 + (t204 * t284 + t291) * t173 + t233 * t220, t141 * t68 - t6 * t189 - t5 * t190 - t21 * t331 - t22 * t289 + t245 * t309 + t261 * t308 - t315, -t5 * t141 + t6 * t142 - t309 * t21 + t308 * t22 + t59 * t242 + t327 * t72, -t231 * t75 + t302 * t39, t82 * t73 + t75 * t81 - (-t216 * t75 - t219 * t73) * t187 + (-t311 - t219 * t40 + (t216 * t73 - t219 * t75) * qJD(6)) * t190, t257 + t324, t223 + t322, t67 * t189 + t289 * t326, t141 * t40 + t232 * t17 + t2 * t189 + t224 * t216 - t219 * t323 - t249 * t289 + t310 * t73, -t1 * t189 - t289 * t12 + t141 * t39 - t231 * t17 + t216 * t323 + t224 * t219 + t310 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177 ^ 2 - t180 ^ 2, -t177 * t95 + t180 * t94 + t195, 0, 0, 0, 0, 0, t239 - t306, -t173 ^ 2 * t220 - t296 + t304, t189 * t68 - t190 * t67 + t245 * t289 + t261 * t331, -t72 * t180 - t5 * t189 + t6 * t190 - t21 * t289 + t22 * t331, 0, 0, 0, 0, 0, t223 - t322, t257 - t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244 * t143, -t143 ^ 2 + t244 ^ 2, t97 + t307, -t98 - t305, t171, t64 * t173 + t244 * t89 + t225, t143 * t89 + t173 * t63 - t237 (-t210 * t67 - t213 * t68) * pkin(4) + (t21 - t25) * t261 + (t22 - t24) * t245, t21 * t24 - t22 * t25 + (t210 * t6 + t213 * t5 + t244 * t72) * pkin(4), t265 * t75 + t311 (t39 - t332) * t219 + (-t326 * t75 - t40) * t216, -t318 - t330, t241 + t319, -t326 * t245, t205 * t40 + t227 * t216 - t219 * t325 - t24 * t73 + t245 * t249, t12 * t245 + t205 * t39 + t216 * t325 + t227 * t219 - t24 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245 ^ 2 - t261 ^ 2, t21 * t245 - t22 * t261 + t59, 0, 0, 0, 0, 0, t241 - t319, -t318 + t330; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 * t73, -t73 ^ 2 + t75 ^ 2, t39 + t332, -t266 + (-qJD(6) + t326) * t75, t67, t12 * t326 - t17 * t75 + t2, t17 * t73 - t249 * t326 - t1;];
tauc_reg  = t11;
