% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:42
% EndTime: 2019-03-09 08:15:51
% DurationCPUTime: 4.03s
% Computational Cost: add. (2834->484), mult. (5740->574), div. (0->0), fcn. (3450->10), ass. (0->251)
t181 = cos(qJ(2));
t178 = sin(qJ(2));
t307 = g(3) * t178;
t179 = sin(qJ(1));
t160 = g(2) * t179;
t182 = cos(qJ(1));
t316 = g(1) * t182 + t160;
t191 = t181 * t316 + t307;
t314 = qJDD(2) * pkin(4) + qJDD(5);
t274 = t181 * qJD(4);
t277 = qJD(2) * t178;
t210 = pkin(7) * t277 + t274;
t272 = qJD(1) * qJD(2);
t252 = t178 * t272;
t115 = qJ(4) * t252;
t269 = t181 * qJDD(1);
t137 = pkin(7) * t269;
t165 = qJDD(2) * qJ(3);
t166 = qJD(2) * qJD(3);
t261 = t137 + t165 + t166;
t237 = -t115 - t261;
t39 = qJ(4) * t269 + qJD(1) * t210 + t237;
t31 = -t39 + t314;
t320 = t31 - t191;
t280 = qJD(1) * t178;
t121 = qJD(6) + t280;
t172 = sin(pkin(9));
t173 = cos(pkin(9));
t177 = sin(qJ(6));
t180 = cos(qJ(6));
t97 = t172 * t177 - t180 * t173;
t205 = t97 * t178;
t88 = t97 * qJD(6);
t303 = qJD(1) * t205 + t88;
t251 = t181 * t272;
t270 = t178 * qJDD(1);
t204 = t251 + t270;
t96 = qJDD(6) + t204;
t98 = t172 * t180 + t173 * t177;
t233 = t121 * t303 - t98 * t96;
t279 = qJD(1) * t181;
t255 = t172 * t279;
t275 = t173 * qJD(2);
t94 = t255 - t275;
t254 = t173 * t279;
t278 = qJD(2) * t172;
t95 = t254 + t278;
t244 = t177 * t95 + t180 * t94;
t319 = t121 * t244;
t312 = pkin(2) + pkin(3);
t318 = t178 * t312;
t256 = t312 * qJD(2);
t253 = t312 * qJDD(2);
t140 = pkin(7) * t279;
t104 = -qJ(4) * t279 + t140;
t167 = qJD(2) * qJ(3);
t91 = -t104 - t167;
t153 = t178 * pkin(4);
t224 = t181 * qJ(5) + t153;
t161 = g(1) * t179;
t308 = g(2) * t182;
t317 = t161 - t308;
t315 = -qJD(6) + t121;
t47 = t177 * t94 - t180 * t95;
t238 = t173 * qJDD(2) - t172 * t269;
t63 = t172 * t252 + t238;
t286 = -t172 * qJDD(2) - t173 * t269;
t64 = t173 * t252 + t286;
t9 = qJD(6) * t47 + t177 * t64 + t180 * t63;
t89 = t98 * qJD(6);
t302 = -t98 * t280 - t89;
t212 = -t302 * t121 + t97 * t96;
t148 = t178 * qJ(3);
t156 = t181 * pkin(2);
t284 = t156 + t148;
t106 = -pkin(1) - t284;
t298 = pkin(7) * qJDD(2);
t92 = -qJD(1) * pkin(1) - pkin(2) * t279 - qJ(3) * t280;
t313 = (qJD(1) * t106 + t92) * qJD(2) - t298;
t311 = pkin(8) * t178;
t310 = pkin(8) * t181;
t309 = g(2) * qJ(4);
t159 = g(3) * t181;
t155 = t181 * pkin(3);
t305 = pkin(7) - qJ(4);
t268 = qJ(5) + t312;
t304 = pkin(8) + t268;
t176 = qJ(3) + pkin(4);
t203 = pkin(4) * t181 - t268 * t178;
t188 = qJD(2) * t203 + t181 * qJD(5);
t147 = t178 * qJD(3);
t171 = qJDD(1) * pkin(1);
t220 = pkin(2) * t269 + t204 * qJ(3) + qJD(1) * t147 + t171;
t199 = pkin(3) * t269 + qJDD(4) + t220;
t13 = qJD(1) * t188 + qJDD(1) * t224 + t199;
t120 = pkin(7) * t251;
t136 = pkin(7) * t270;
t248 = qJDD(3) + t120 + t136;
t271 = qJD(1) * qJD(4);
t186 = -qJ(4) * t204 - t178 * t271 + t248;
t29 = -qJD(2) * qJD(5) - t268 * qJDD(2) + t186;
t7 = t172 * t13 + t173 * t29;
t66 = pkin(3) * t279 + qJD(4) - t92;
t48 = qJD(1) * t224 + t66;
t131 = qJ(4) * t280;
t266 = pkin(7) * t280;
t230 = qJD(3) + t266;
t213 = -t131 + t230;
t62 = -t268 * qJD(2) + t213;
t19 = t172 * t48 + t173 * t62;
t276 = qJD(2) * t181;
t285 = qJ(3) * t276 + t147;
t40 = t188 + t285;
t84 = -t178 * qJD(4) + t305 * t276;
t22 = t172 * t40 + t173 * t84;
t133 = qJ(3) * t279;
t51 = qJD(1) * t203 + t133;
t28 = t173 * t104 + t172 * t51;
t109 = t305 * t178;
t260 = t155 + t284;
t218 = t260 + t224;
t61 = pkin(1) + t218;
t38 = t173 * t109 + t172 * t61;
t301 = qJD(2) * pkin(2);
t299 = t47 * t121;
t297 = qJ(3) * t181;
t296 = qJDD(2) * pkin(2);
t168 = t178 ^ 2;
t185 = qJD(1) ^ 2;
t295 = t168 * t185;
t294 = t178 * t179;
t293 = t178 * t182;
t292 = t178 * t185;
t291 = t179 * t181;
t290 = t181 * t182;
t289 = t91 * qJD(2);
t258 = pkin(5) * t172 - pkin(7);
t232 = t178 * t258;
t288 = -qJD(1) * t232 + qJD(3) - t131;
t287 = -qJD(4) - t66;
t169 = t181 ^ 2;
t282 = t168 - t169;
t281 = t168 + t169;
t102 = -t131 + t266;
t273 = qJD(3) + t102;
t267 = t172 * t311;
t265 = t181 * t292;
t264 = t244 * t279;
t263 = t47 * t279;
t262 = t137 + 0.2e1 * t165 + 0.2e1 * t166;
t259 = -g(1) * t293 - g(2) * t294 + t159;
t6 = t173 * t13 - t172 * t29;
t2 = pkin(5) * t204 - t64 * pkin(8) + t6;
t5 = -pkin(8) * t63 + t7;
t257 = -t177 * t5 + t180 * t2;
t250 = t172 * t270;
t249 = t173 * t270;
t247 = -pkin(1) - t148;
t21 = -t172 * t84 + t173 * t40;
t18 = -t172 * t62 + t173 * t48;
t93 = pkin(1) + t260;
t243 = qJD(1) * t93 + t66;
t27 = -t172 * t104 + t173 * t51;
t37 = -t109 * t172 + t173 * t61;
t242 = t287 * t178;
t240 = g(1) * t268;
t239 = 0.2e1 * t251;
t236 = t182 * pkin(1) + pkin(2) * t290 + t179 * pkin(7) + qJ(3) * t293;
t235 = t136 + t259;
t234 = t178 * t256;
t184 = qJD(2) ^ 2;
t231 = pkin(7) * t184 + t308;
t227 = -t6 * t172 + t7 * t173;
t226 = t177 * t2 + t180 * t5;
t157 = t182 * pkin(7);
t225 = g(1) * (-t182 * qJ(4) + t157);
t12 = pkin(5) * t280 + pkin(8) * t95 + t18;
t14 = pkin(8) * t94 + t19;
t3 = t12 * t180 - t14 * t177;
t4 = t12 * t177 + t14 * t180;
t223 = t172 * t18 - t173 * t19;
t24 = pkin(5) * t178 + t173 * t310 + t37;
t30 = t172 * t310 + t38;
t222 = -t177 * t30 + t180 * t24;
t221 = t177 * t24 + t180 * t30;
t219 = pkin(3) * t290 + t236;
t105 = t230 - t301;
t107 = t140 + t167;
t217 = t105 * t181 - t107 * t178;
t216 = -qJDD(3) - t235;
t215 = pkin(5) * t181 - t173 * t311;
t214 = t247 - t156;
t211 = -0.2e1 * pkin(1) * t272 - t298;
t8 = t244 * qJD(6) - t177 * t63 + t180 * t64;
t209 = t120 - t216;
t100 = t304 * t172;
t208 = -qJD(1) * t267 + qJD(5) * t173 - qJD(6) * t100 + t28;
t101 = t304 * t173;
t207 = -qJD(1) * t215 + qJD(5) * t172 + qJD(6) * t101 - t27;
t206 = t178 * (-t172 * t19 - t173 * t18);
t202 = t7 * t172 + t6 * t173 - t308;
t201 = -t181 * t31 + t316;
t200 = -t231 + 0.2e1 * t171;
t198 = t172 * t295 - t249;
t197 = -t173 * t295 - t250;
t196 = -qJ(4) * qJDD(1) - t316;
t72 = qJD(2) * pkin(4) + qJD(5) - t91;
t195 = t268 * t276 + (qJD(5) - t72) * t178;
t26 = -qJD(1) * t234 + t199;
t65 = -t234 + t285;
t193 = -qJD(1) * t65 - qJDD(1) * t93 - t26 + t308;
t190 = -qJ(4) * t270 + t209;
t43 = pkin(2) * t252 - t220;
t85 = pkin(2) * t277 - t285;
t189 = -qJD(1) * t85 - qJDD(1) * t106 - t231 - t43;
t71 = -pkin(7) * t252 + t261;
t80 = t248 - t296;
t187 = qJD(2) * t217 + t80 * t178 + t71 * t181;
t164 = pkin(9) + qJ(6);
t150 = t181 * qJ(4);
t146 = cos(t164);
t145 = sin(t164);
t134 = qJ(4) * t277;
t126 = g(1) * t291;
t125 = g(1) * t294;
t119 = qJ(3) * t290;
t117 = qJ(3) * t291;
t114 = pkin(5) * t173 + t176;
t113 = -t184 - t295;
t110 = pkin(7) * t181 - t150;
t108 = qJDD(2) + t265;
t103 = pkin(2) * t280 - t133;
t87 = -t181 * t258 - t150;
t83 = -t134 + t210;
t82 = -t312 * t280 + t133;
t79 = t97 * t181;
t78 = t98 * t181;
t77 = -t256 + t213;
t76 = -t145 * t179 + t146 * t293;
t75 = -t145 * t293 - t146 * t179;
t74 = -t145 * t182 - t146 * t294;
t73 = t145 * t294 - t146 * t182;
t59 = qJD(2) * t232 + t134 - t274;
t44 = -pkin(5) * t94 + t72;
t36 = -t253 + t186;
t35 = -t181 * t88 - t98 * t277;
t34 = -qJD(2) * t205 + t181 * t89;
t17 = t63 * pkin(5) + t31;
t16 = -qJD(2) * t267 + t22;
t15 = qJD(2) * t215 + t21;
t1 = [qJDD(1), t317, t316, qJDD(1) * t168 + t178 * t239, 0.2e1 * t178 * t269 - 0.2e1 * t282 * t272, qJDD(2) * t178 + t181 * t184, qJDD(2) * t181 - t178 * t184, 0, t178 * t211 + t181 * t200 + t126, -t178 * t200 + t181 * t211 - t125, t313 * t178 + t189 * t181 + t126, pkin(7) * qJDD(1) * t281 + t187 - t316, t189 * t178 - t313 * t181 + t125, pkin(7) * t187 - g(1) * t157 - g(2) * t236 + t43 * t106 - t214 * t161 + t92 * t85, t110 * qJDD(2) + t125 + (t181 * t243 - t83) * qJD(2) - t193 * t178, t109 * qJDD(2) - t126 + (t178 * t243 + t84) * qJD(2) + t193 * t181 (-qJD(2) * t77 - qJDD(1) * t110 + t39 + (-qJD(2) * t109 + t83) * qJD(1)) * t181 + (-t289 - qJDD(1) * t109 - t36 + (qJD(2) * t110 - t84) * qJD(1)) * t178 + t316, t36 * t109 + t77 * t84 - t39 * t110 + t91 * t83 + t26 * t93 + t66 * t65 - t225 - g(2) * t219 + (-g(1) * (t214 - t155) + t309) * t179, t110 * t63 + t83 * t94 + (qJD(1) * t37 + t18) * t276 + t201 * t172 + (t21 * qJD(1) + t37 * qJDD(1) + t173 * t317 + t278 * t72 + t6) * t178, t110 * t64 + t83 * t95 + (-qJD(1) * t38 - t19) * t276 + t201 * t173 + (-t22 * qJD(1) - t38 * qJDD(1) - t172 * t317 + t275 * t72 - t7) * t178, qJD(2) * t206 + t181 * t202 + t21 * t95 + t22 * t94 - t37 * t64 - t38 * t63 + t126, t7 * t38 + t19 * t22 + t6 * t37 + t18 * t21 + t31 * t110 - t72 * t83 - t225 - g(2) * (pkin(4) * t293 + qJ(5) * t290 + t219) + (-g(1) * (t247 - t153) + t309 + t181 * t240) * t179, t34 * t47 + t79 * t8, t244 * t34 + t35 * t47 + t78 * t8 - t79 * t9, t121 * t34 + t178 * t8 + t276 * t47 + t79 * t96, t121 * t35 - t178 * t9 + t244 * t276 + t78 * t96, t121 * t276 + t178 * t96 (t180 * t15 - t177 * t16) * t121 + t222 * t96 + t257 * t178 + t3 * t276 - t59 * t244 + t87 * t9 - t17 * t78 - t44 * t35 - g(1) * t74 - g(2) * t76 + (-t121 * t221 - t178 * t4) * qJD(6) -(t177 * t15 + t180 * t16) * t121 - t221 * t96 - t226 * t178 - t4 * t276 + t59 * t47 + t87 * t8 + t17 * t79 + t44 * t34 - g(1) * t73 - g(2) * t75 + (-t121 * t222 - t178 * t3) * qJD(6); 0, 0, 0, -t265, t282 * t185, t270, t269, qJDD(2), pkin(1) * t292 - t235, t307 - t137 + (pkin(1) * t185 + t316) * t181, 0.2e1 * t296 + (t103 * t181 - t178 * t92) * qJD(1) + t216 (-pkin(2) * t178 + t297) * qJDD(1) + ((t107 - t167) * t178 + (qJD(3) - t105 - t301) * t181) * qJD(1) (qJD(1) * t103 - g(3)) * t178 + (qJD(1) * t92 - t316) * t181 + t262, t71 * qJ(3) + t107 * qJD(3) - t80 * pkin(2) - t92 * t103 - g(1) * (-pkin(2) * t293 + t119) - g(2) * (-pkin(2) * t294 + t117) - g(3) * t284 - t217 * qJD(1) * pkin(7), t102 * qJD(2) + t115 + (-g(3) + (-pkin(7) * qJD(2) - t82) * qJD(1)) * t178 + (qJD(1) * t287 + t196) * t181 + t262, -t104 * qJD(2) - 0.2e1 * t253 + ((-qJ(4) * qJD(2) + t82) * t181 + t242) * qJD(1) + t190 (-t297 + t318) * qJDD(1) + (-t273 + t77 + t256) * t279, -g(1) * t119 - g(2) * t117 - g(3) * t260 - t39 * qJ(3) - t77 * t104 - t273 * t91 - t312 * t36 + t316 * t318 - t66 * t82, t268 * t250 + t176 * t63 - t273 * t94 + t320 * t173 + (t172 * t195 - t178 * t27 - t18 * t181) * qJD(1), t268 * t249 + t176 * t64 - t273 * t95 - t320 * t172 + (t173 * t195 + t178 * t28 + t19 * t181) * qJD(1), -t27 * t95 - t28 * t94 + (-qJD(5) * t94 + t18 * t280 + t268 * t63 - t7) * t173 + (qJD(5) * t95 + t19 * t280 - t268 * t64 + t6) * t172 - t259, t31 * t176 - t19 * t28 - t18 * t27 - g(1) * (pkin(4) * t290 + t119) - g(2) * (pkin(4) * t291 + t117) - g(3) * t218 + t273 * t72 - t227 * t268 + t223 * qJD(5) + (t268 * t160 + t182 * t240) * t178, t303 * t47 - t8 * t98, t244 * t303 - t302 * t47 + t8 * t97 + t98 * t9, t233 - t263, t212 - t264, -t121 * t279 (t100 * t180 + t101 * t177) * t96 + t114 * t9 - t17 * t97 - t3 * t279 - t288 * t244 + t302 * t44 + (t177 * t208 + t180 * t207) * t121 - t191 * t146 -(t100 * t177 - t101 * t180) * t96 + t114 * t8 - t17 * t98 + t4 * t279 + t288 * t47 + t303 * t44 + (-t177 * t207 + t180 * t208) * t121 + t191 * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, t270, t113, -qJD(2) * t107 + t280 * t92 + t209 - t296, t113, t108, -t270, t289 - t253 + (-qJ(4) * t276 + t242) * qJD(1) + t190 (t94 - t255) * qJD(2) + t197 (t95 - t254) * qJD(2) + t198, t172 * t64 - t173 * t63 + (-t172 * t94 - t173 * t95) * t280, qJD(1) * t206 - t72 * qJD(2) + t227 + t259, 0, 0, 0, 0, 0, qJD(2) * t244 + t233, -qJD(2) * t47 + t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t239 + t270, 0.2e1 * t252 - t269, -t281 * t185 (-t181 * t91 + (t77 - t256) * t178) * qJD(1) + t199 + t317 (-t94 + t275) * t279 - t198 (-t95 - t278) * t279 + t197, -t172 * t63 - t173 * t64 + (-t172 * t95 + t173 * t94) * t280, t161 + (-t178 * t223 + t181 * t72) * qJD(1) + t202, 0, 0, 0, 0, 0, -t212 - t264, t233 + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t95 + t278) * t280 + t238 (t94 + t275) * t280 + t286, -t94 ^ 2 - t95 ^ 2, -t18 * t95 - t19 * t94 + (-pkin(7) * t272 - g(3)) * t178 + (t196 - t271) * t181 - t237 + t314, 0, 0, 0, 0, 0, t9 + t299, t8 + t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 * t244, -t244 ^ 2 + t47 ^ 2, t8 - t319, t299 - t9, t96, -g(1) * t75 + g(2) * t73 - t145 * t159 + t315 * t4 - t44 * t47 + t257, g(1) * t76 - g(2) * t74 - t146 * t159 - t244 * t44 + t315 * t3 - t226;];
tau_reg  = t1;
