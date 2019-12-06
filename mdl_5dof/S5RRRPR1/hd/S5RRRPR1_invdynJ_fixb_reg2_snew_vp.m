% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:55
% EndTime: 2019-12-05 18:39:07
% DurationCPUTime: 6.03s
% Computational Cost: add. (37204->448), mult. (86309->644), div. (0->0), fcn. (63074->10), ass. (0->278)
t335 = -2 * qJD(4);
t252 = sin(pkin(9));
t255 = sin(qJ(3));
t259 = cos(qJ(3));
t260 = cos(qJ(2));
t256 = sin(qJ(2));
t290 = qJD(1) * t256;
t223 = -t259 * t260 * qJD(1) + t255 * t290;
t225 = (t255 * t260 + t256 * t259) * qJD(1);
t253 = cos(pkin(9));
t204 = t253 * t223 + t225 * t252;
t206 = -t223 * t252 + t225 * t253;
t176 = t206 * t204;
t248 = qJDD(2) + qJDD(3);
t326 = -t176 + t248;
t334 = t252 * t326;
t333 = t253 * t326;
t254 = sin(qJ(5));
t258 = cos(qJ(5));
t166 = t258 * t204 + t206 * t254;
t168 = -t204 * t254 + t206 * t258;
t129 = t168 * t166;
t240 = qJDD(5) + t248;
t328 = -t129 + t240;
t332 = t254 * t328;
t212 = t225 * t223;
t325 = -t212 + t248;
t331 = t255 * t325;
t330 = t258 * t328;
t329 = t259 * t325;
t242 = t256 * qJDD(1);
t285 = qJD(1) * qJD(2);
t284 = t260 * t285;
t231 = t242 + t284;
t243 = t260 * qJDD(1);
t283 = t256 * t285;
t232 = t243 - t283;
t278 = t231 * t255 - t259 * t232;
t189 = -qJD(3) * t225 - t278;
t271 = t231 * t259 + t232 * t255;
t190 = -qJD(3) * t223 + t271;
t156 = t189 * t252 + t190 * t253;
t249 = qJD(2) + qJD(3);
t196 = t249 * t204;
t327 = t156 - t196;
t140 = t156 + t196;
t219 = t249 * t223;
t182 = t190 + t219;
t262 = qJD(1) ^ 2;
t293 = t256 * t262;
t257 = sin(qJ(1));
t322 = cos(qJ(1));
t269 = t322 * g(1) + t257 * g(2);
t314 = qJDD(1) * pkin(6);
t227 = -t262 * pkin(1) - t269 + t314;
t302 = t227 * t256;
t186 = qJDD(2) * pkin(2) - pkin(7) * t231 - t302 + (pkin(2) * t293 + pkin(7) * t285 - g(3)) * t260;
t216 = -t256 * g(3) + t260 * t227;
t251 = t260 ^ 2;
t245 = t251 * t262;
t268 = qJD(2) * pkin(2) - pkin(7) * t290;
t187 = -pkin(2) * t245 + t232 * pkin(7) - qJD(2) * t268 + t216;
t152 = -t259 * t186 + t187 * t255;
t116 = t325 * pkin(3) - t182 * qJ(4) - t152;
t153 = t255 * t186 + t259 * t187;
t221 = t223 ^ 2;
t274 = pkin(3) * t249 - qJ(4) * t225;
t118 = -t221 * pkin(3) + t189 * qJ(4) - t249 * t274 + t153;
t275 = t253 * t116 - t118 * t252 + t206 * t335;
t73 = t252 * t116 + t253 * t118 + t204 * t335;
t164 = t166 ^ 2;
t165 = t168 ^ 2;
t201 = t204 ^ 2;
t202 = t206 ^ 2;
t222 = t225 ^ 2;
t241 = qJD(5) + t249;
t239 = t241 ^ 2;
t324 = t249 ^ 2;
t58 = t326 * pkin(4) - t140 * pkin(8) + t275;
t155 = t253 * t189 - t190 * t252;
t276 = pkin(4) * t249 - pkin(8) * t206;
t59 = -t201 * pkin(4) + t155 * pkin(8) - t249 * t276 + t73;
t32 = t254 * t59 - t258 * t58;
t33 = t254 * t58 + t258 * t59;
t17 = t254 * t33 - t258 * t32;
t18 = t254 * t32 + t258 * t33;
t319 = t17 * t253;
t8 = t18 * t252 + t319;
t323 = pkin(3) * t8 + pkin(4) * t17;
t280 = -t258 * t155 + t156 * t254;
t264 = (-qJD(5) + t241) * t168 - t280;
t162 = t241 * t166;
t272 = t155 * t254 + t156 * t258;
t96 = -qJD(5) * t166 + t272;
t83 = t162 + t96;
t50 = t254 * t264 - t258 * t83;
t52 = t254 * t83 + t258 * t264;
t26 = t252 * t52 + t253 * t50;
t321 = pkin(3) * t26 + pkin(4) * t50;
t320 = t17 * t252;
t282 = g(1) * t257 - t322 * g(2);
t267 = qJDD(1) * pkin(1) + t282;
t192 = pkin(2) * t232 - t268 * t290 + (pkin(7) * t251 + pkin(6)) * t262 + t267;
t130 = pkin(3) * t189 + qJ(4) * t221 - t225 * t274 - qJDD(4) + t192;
t85 = pkin(4) * t155 + pkin(8) * t201 - t206 * t276 + t130;
t318 = t254 * t85;
t42 = t252 * t73 + t253 * t275;
t317 = t255 * t42;
t316 = t258 * t85;
t315 = t259 * t42;
t112 = -t152 * t259 + t153 * t255;
t313 = t112 * t256;
t123 = t129 + t240;
t312 = t123 * t254;
t311 = t123 * t258;
t310 = t130 * t252;
t309 = t130 * t253;
t171 = t176 + t248;
t308 = t171 * t252;
t307 = t171 * t253;
t306 = t192 * t255;
t305 = t192 * t259;
t209 = t212 + t248;
t304 = t209 * t255;
t303 = t209 * t259;
t301 = t241 * t254;
t300 = t241 * t258;
t299 = t249 * t206;
t298 = t249 * t252;
t297 = t249 * t253;
t296 = t249 * t255;
t295 = t249 * t259;
t238 = t260 * t293;
t294 = t256 * (qJDD(2) + t238);
t292 = t260 * (qJDD(2) - t238);
t287 = qJD(3) + t249;
t286 = qJD(5) + t241;
t43 = -t252 * t275 + t253 * t73;
t113 = t152 * t255 + t259 * t153;
t215 = g(3) * t260 + t302;
t279 = t256 * t215 + t260 * t216;
t193 = -t202 - t324;
t144 = t193 * t253 - t308;
t277 = pkin(3) * t144 - t73;
t121 = -t239 - t164;
t89 = t121 * t254 + t330;
t90 = t121 * t258 - t332;
t55 = t252 * t90 + t253 * t89;
t273 = pkin(3) * t55 + pkin(4) * t89 - t32;
t169 = -t324 - t201;
t126 = t169 * t252 + t333;
t270 = pkin(3) * t126 + t275;
t266 = t155 + t299;
t157 = -t165 - t239;
t105 = t157 * t258 - t312;
t106 = -t157 * t254 - t311;
t64 = t105 * t253 + t106 * t252;
t265 = pkin(3) * t64 + pkin(4) * t105 - t33;
t263 = (-qJD(3) + t249) * t225 - t278;
t261 = qJD(2) ^ 2;
t250 = t256 ^ 2;
t244 = t250 * t262;
t233 = t243 - 0.2e1 * t283;
t230 = t242 + 0.2e1 * t284;
t226 = pkin(6) * t262 + t267;
t218 = -t222 + t324;
t217 = t221 - t324;
t214 = -t222 - t324;
t211 = t222 - t221;
t207 = -t324 - t221;
t195 = -t202 + t324;
t194 = t201 - t324;
t191 = -t221 - t222;
t184 = -t214 * t255 - t303;
t183 = t214 * t259 - t304;
t181 = t190 - t219;
t180 = -t287 * t223 + t271;
t177 = t287 * t225 + t278;
t175 = t207 * t259 - t331;
t174 = t207 * t255 + t329;
t173 = t202 - t201;
t161 = -t165 + t239;
t160 = t164 - t239;
t159 = (-t204 * t253 + t206 * t252) * t249;
t158 = (-t204 * t252 - t206 * t253) * t249;
t151 = -t201 - t202;
t149 = t194 * t253 - t308;
t148 = -t195 * t252 + t333;
t147 = t194 * t252 + t307;
t146 = t195 * t253 + t334;
t145 = -t193 * t252 - t307;
t142 = t182 * t255 + t259 * t263;
t141 = -t182 * t259 + t255 * t263;
t135 = -t155 + t299;
t134 = t156 * t253 - t206 * t298;
t133 = t156 * t252 + t206 * t297;
t132 = -t155 * t252 + t204 * t297;
t131 = t155 * t253 + t204 * t298;
t128 = t165 - t164;
t127 = t169 * t253 - t334;
t120 = (-t166 * t258 + t168 * t254) * t241;
t119 = (-t166 * t254 - t168 * t258) * t241;
t111 = -t164 - t165;
t110 = t160 * t258 - t312;
t109 = -t161 * t254 + t330;
t108 = t160 * t254 + t311;
t107 = t161 * t258 + t332;
t103 = -t144 * t255 + t145 * t259;
t102 = t144 * t259 + t145 * t255;
t101 = -qJ(4) * t144 - t309;
t100 = t140 * t252 + t253 * t266;
t99 = -t135 * t253 - t252 * t327;
t98 = -t140 * t253 + t252 * t266;
t97 = -t135 * t252 + t253 * t327;
t95 = -qJD(5) * t168 - t280;
t94 = pkin(3) * t98;
t93 = -qJ(4) * t126 - t310;
t92 = -t126 * t255 + t127 * t259;
t91 = t126 * t259 + t127 * t255;
t87 = -t119 * t252 + t120 * t253;
t86 = t119 * t253 + t120 * t252;
t84 = -t286 * t166 + t272;
t82 = -t162 + t96;
t79 = t286 * t168 + t280;
t78 = -t168 * t301 + t258 * t96;
t77 = t168 * t300 + t254 * t96;
t76 = t166 * t300 - t254 * t95;
t75 = t166 * t301 + t258 * t95;
t74 = -pkin(3) * t327 + qJ(4) * t145 - t310;
t71 = -pkin(3) * t135 + qJ(4) * t127 + t309;
t69 = -t108 * t252 + t110 * t253;
t68 = -t107 * t252 + t109 * t253;
t67 = t108 * t253 + t110 * t252;
t66 = t107 * t253 + t109 * t252;
t65 = -t105 * t252 + t106 * t253;
t62 = t100 * t259 - t255 * t98;
t61 = t100 * t255 + t259 * t98;
t60 = -pkin(8) * t105 - t316;
t56 = -t252 * t89 + t253 * t90;
t53 = -pkin(8) * t89 - t318;
t51 = -t254 * t82 - t258 * t79;
t49 = -t254 * t79 + t258 * t82;
t47 = -t252 * t77 + t253 * t78;
t46 = -t252 * t75 + t253 * t76;
t45 = t252 * t78 + t253 * t77;
t44 = t252 * t76 + t253 * t75;
t41 = pkin(3) * t42;
t40 = -pkin(4) * t84 + pkin(8) * t106 - t318;
t39 = pkin(3) * t130 + qJ(4) * t43;
t38 = -t255 * t64 + t259 * t65;
t37 = t255 * t65 + t259 * t64;
t36 = -pkin(4) * t79 + pkin(8) * t90 + t316;
t35 = -qJ(4) * t98 - t42;
t34 = -pkin(3) * t151 + qJ(4) * t100 + t43;
t30 = -t255 * t55 + t259 * t56;
t29 = t255 * t56 + t259 * t55;
t28 = -t252 * t50 + t253 * t52;
t27 = -t252 * t49 + t253 * t51;
t25 = t252 * t51 + t253 * t49;
t23 = t259 * t43 - t317;
t22 = t255 * t43 + t315;
t21 = -qJ(4) * t64 - t252 * t40 + t253 * t60;
t20 = -qJ(4) * t55 - t252 * t36 + t253 * t53;
t19 = -pkin(3) * t84 + qJ(4) * t65 + t252 * t60 + t253 * t40;
t15 = -pkin(3) * t79 + qJ(4) * t56 + t252 * t53 + t253 * t36;
t14 = -t255 * t26 + t259 * t28;
t13 = t255 * t28 + t259 * t26;
t12 = pkin(4) * t85 + pkin(8) * t18;
t11 = -pkin(8) * t50 - t17;
t10 = -pkin(4) * t111 + pkin(8) * t52 + t18;
t9 = t18 * t253 - t320;
t6 = -qJ(4) * t26 - t10 * t252 + t11 * t253;
t5 = -pkin(3) * t111 + qJ(4) * t28 + t10 * t253 + t11 * t252;
t4 = -t255 * t8 + t259 * t9;
t3 = t255 * t9 + t259 * t8;
t2 = -pkin(8) * t319 - qJ(4) * t8 - t12 * t252;
t1 = pkin(3) * t85 - pkin(8) * t320 + qJ(4) * t9 + t12 * t253;
t7 = [0, 0, 0, 0, 0, qJDD(1), t282, t269, 0, 0, (t231 + t284) * t256, t230 * t260 + t233 * t256, t294 + t260 * (-t244 + t261), (t232 - t283) * t260, t256 * (t245 - t261) + t292, 0, t260 * t226 + pkin(1) * t233 + pkin(6) * (t260 * (-t245 - t261) - t294), -t256 * t226 - pkin(1) * t230 + pkin(6) * (-t292 - t256 * (-t244 - t261)), pkin(1) * (t244 + t245) + (t250 + t251) * t314 + t279, pkin(1) * t226 + pkin(6) * t279, t256 * (t190 * t259 - t225 * t296) + t260 * (t190 * t255 + t225 * t295), t256 * (-t177 * t259 - t181 * t255) + t260 * (-t177 * t255 + t181 * t259), t256 * (-t218 * t255 + t329) + t260 * (t218 * t259 + t331), t256 * (-t189 * t255 + t223 * t295) + t260 * (t189 * t259 + t223 * t296), t256 * (t217 * t259 - t304) + t260 * (t217 * t255 + t303), (t256 * (-t223 * t259 + t225 * t255) + t260 * (-t223 * t255 - t225 * t259)) * t249, t256 * (-pkin(7) * t174 - t306) + t260 * (-pkin(2) * t177 + pkin(7) * t175 + t305) - pkin(1) * t177 + pkin(6) * (-t174 * t256 + t175 * t260), t256 * (-pkin(7) * t183 - t305) + t260 * (-pkin(2) * t180 + pkin(7) * t184 - t306) - pkin(1) * t180 + pkin(6) * (-t183 * t256 + t184 * t260), t256 * (-pkin(7) * t141 - t112) + t260 * (-pkin(2) * t191 + pkin(7) * t142 + t113) - pkin(1) * t191 + pkin(6) * (-t141 * t256 + t142 * t260), -pkin(7) * t313 + t260 * (pkin(2) * t192 + pkin(7) * t113) + pkin(1) * t192 + pkin(6) * (t113 * t260 - t313), t256 * (-t133 * t255 + t134 * t259) + t260 * (t133 * t259 + t134 * t255), t256 * (-t255 * t97 + t259 * t99) + t260 * (t255 * t99 + t259 * t97), t256 * (-t146 * t255 + t148 * t259) + t260 * (t146 * t259 + t148 * t255), t256 * (-t131 * t255 + t132 * t259) + t260 * (t131 * t259 + t132 * t255), t256 * (-t147 * t255 + t149 * t259) + t260 * (t147 * t259 + t149 * t255), t256 * (-t158 * t255 + t159 * t259) + t260 * (t158 * t259 + t159 * t255), t256 * (-pkin(7) * t91 - t255 * t71 + t259 * t93) + t260 * (-pkin(2) * t135 + pkin(7) * t92 + t255 * t93 + t259 * t71) - pkin(1) * t135 + pkin(6) * (-t256 * t91 + t260 * t92), t256 * (-pkin(7) * t102 + t101 * t259 - t255 * t74) + t260 * (-pkin(2) * t327 + pkin(7) * t103 + t101 * t255 + t259 * t74) - pkin(1) * t327 + pkin(6) * (-t102 * t256 + t103 * t260), t256 * (-pkin(7) * t61 - t255 * t34 + t259 * t35) + t260 * (-pkin(2) * t151 + pkin(7) * t62 + t255 * t35 + t259 * t34) - pkin(1) * t151 + pkin(6) * (-t256 * t61 + t260 * t62), t256 * (-pkin(7) * t22 - qJ(4) * t315 - t255 * t39) + t260 * (pkin(2) * t130 + pkin(7) * t23 - qJ(4) * t317 + t259 * t39) + pkin(1) * t130 + pkin(6) * (-t22 * t256 + t23 * t260), t256 * (-t255 * t45 + t259 * t47) + t260 * (t255 * t47 + t259 * t45), t256 * (-t25 * t255 + t259 * t27) + t260 * (t25 * t259 + t255 * t27), t256 * (-t255 * t66 + t259 * t68) + t260 * (t255 * t68 + t259 * t66), t256 * (-t255 * t44 + t259 * t46) + t260 * (t255 * t46 + t259 * t44), t256 * (-t255 * t67 + t259 * t69) + t260 * (t255 * t69 + t259 * t67), t256 * (-t255 * t86 + t259 * t87) + t260 * (t255 * t87 + t259 * t86), t256 * (-pkin(7) * t29 - t15 * t255 + t20 * t259) + t260 * (-pkin(2) * t79 + pkin(7) * t30 + t15 * t259 + t20 * t255) - pkin(1) * t79 + pkin(6) * (-t256 * t29 + t260 * t30), t256 * (-pkin(7) * t37 - t19 * t255 + t21 * t259) + t260 * (-pkin(2) * t84 + pkin(7) * t38 + t19 * t259 + t21 * t255) - pkin(1) * t84 + pkin(6) * (-t256 * t37 + t260 * t38), t256 * (-pkin(7) * t13 - t255 * t5 + t259 * t6) + t260 * (-pkin(2) * t111 + pkin(7) * t14 + t255 * t6 + t259 * t5) - pkin(1) * t111 + pkin(6) * (-t13 * t256 + t14 * t260), t256 * (-pkin(7) * t3 - t1 * t255 + t2 * t259) + t260 * (pkin(2) * t85 + pkin(7) * t4 + t1 * t259 + t2 * t255) + pkin(1) * t85 + pkin(6) * (-t256 * t3 + t260 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, t244 - t245, t242, t238, t243, qJDD(2), -t215, -t216, 0, 0, t212, t211, t182, -t212, t263, t248, pkin(2) * t174 - t152, pkin(2) * t183 - t153, pkin(2) * t141, pkin(2) * t112, t176, t173, t140, -t176, t266, t248, pkin(2) * t91 + t270, pkin(2) * t102 + t277, pkin(2) * t61 + t94, pkin(2) * t22 + t41, t129, t128, t83, -t129, t264, t240, pkin(2) * t29 + t273, pkin(2) * t37 + t265, pkin(2) * t13 + t321, pkin(2) * t3 + t323; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t211, t182, -t212, t263, t248, -t152, -t153, 0, 0, t176, t173, t140, -t176, t266, t248, t270, t277, t94, t41, t129, t128, t83, -t129, t264, t240, t273, t265, t321, t323; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t327, t151, -t130, 0, 0, 0, 0, 0, 0, t79, t84, t111, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t128, t83, -t129, t264, t240, -t32, -t33, 0, 0;];
tauJ_reg = t7;
