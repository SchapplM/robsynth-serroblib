% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PPRRRR2
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:28
% EndTime: 2019-03-08 19:05:38
% DurationCPUTime: 4.89s
% Computational Cost: add. (4094->450), mult. (10244->685), div. (0->0), fcn. (9480->18), ass. (0->232)
t177 = sin(qJ(4));
t181 = cos(qJ(4));
t230 = pkin(4) * t177 - pkin(10) * t181;
t138 = t230 * qJD(4);
t303 = cos(pkin(6));
t153 = qJD(1) * t303 + qJD(2);
t172 = sin(pkin(6));
t170 = sin(pkin(13));
t178 = sin(qJ(3));
t182 = cos(qJ(3));
t173 = cos(pkin(13));
t174 = cos(pkin(7));
t290 = t173 * t174;
t211 = t170 * t182 + t178 * t290;
t199 = t211 * t172;
t171 = sin(pkin(7));
t292 = t171 * t178;
t75 = qJD(1) * t199 + t153 * t292;
t333 = t138 - t75;
t176 = sin(qJ(5));
t180 = cos(qJ(5));
t265 = t180 * qJD(4);
t279 = qJD(3) * t177;
t127 = t176 * t279 - t265;
t274 = qJD(4) * t176;
t129 = t180 * t279 + t274;
t175 = sin(qJ(6));
t179 = cos(qJ(6));
t219 = t127 * t175 - t179 * t129;
t81 = t179 * t127 + t129 * t175;
t332 = t219 * t81;
t151 = t303 * qJDD(1) + qJDD(2);
t280 = qJD(1) * t172;
t255 = t173 * t280;
t235 = t174 * t255;
t263 = qJDD(1) * t172;
t243 = t173 * t263;
t278 = qJD(3) * t178;
t254 = t171 * t278;
t256 = t170 * t280;
t293 = t170 * t178;
t260 = t172 * t293;
t276 = qJD(3) * t182;
t187 = -t182 * (t151 * t171 + t174 * t243) + qJDD(1) * t260 + t153 * t254 + t235 * t278 + t256 * t276;
t302 = cos(pkin(12));
t227 = t303 * t302;
t301 = sin(pkin(12));
t110 = t170 * t227 + t173 * t301;
t193 = t170 * t301 - t173 * t227;
t239 = t172 * t302;
t321 = t171 * t239 + t193 * t174;
t58 = t110 * t178 + t182 * t321;
t226 = t303 * t301;
t111 = -t170 * t226 + t173 * t302;
t194 = t170 * t302 + t173 * t226;
t238 = t172 * t301;
t320 = -t171 * t238 + t194 * t174;
t60 = t111 * t178 + t182 * t320;
t240 = t171 * t303;
t233 = t182 * t240;
t259 = t182 * t290;
t86 = -t172 * t259 - t233 + t260;
t207 = g(1) * t60 + g(2) * t58 + g(3) * t86;
t331 = t75 * qJD(3) - t187 + t207;
t269 = qJD(5) * t177;
t330 = -qJD(3) * t269 + qJDD(4);
t329 = t219 ^ 2 - t81 ^ 2;
t277 = qJD(3) * t181;
t154 = -qJD(5) + t277;
t149 = -qJD(6) + t154;
t266 = qJD(6) * t179;
t267 = qJD(6) * t175;
t264 = qJD(3) * qJD(4);
t246 = t181 * t264;
t262 = t177 * qJDD(3);
t76 = qJD(5) * t265 + (t246 + t262) * t180 + t330 * t176;
t77 = t176 * (qJD(4) * (qJD(5) + t277) + t262) - t330 * t180;
t22 = -t127 * t266 - t129 * t267 - t175 * t77 + t179 * t76;
t328 = -t149 * t81 + t22;
t169 = qJ(5) + qJ(6);
t164 = sin(t169);
t165 = cos(t169);
t103 = t153 * t174 - t171 * t255;
t70 = qJD(3) * pkin(9) + t75;
t52 = t177 * t103 + t181 * t70;
t50 = qJD(4) * pkin(10) + t52;
t142 = -pkin(4) * t181 - pkin(10) * t177 - pkin(3);
t74 = -t178 * t256 + t182 * (t153 * t171 + t235);
t65 = qJD(3) * t142 - t74;
t21 = t176 * t65 + t180 * t50;
t19 = -pkin(11) * t127 + t21;
t17 = t19 * t267;
t59 = t110 * t182 - t178 * t321;
t88 = t171 * t193 - t174 * t239;
t34 = t177 * t88 + t181 * t59;
t61 = t111 * t182 - t178 * t320;
t89 = t171 * t194 + t174 * t238;
t36 = t177 * t89 + t181 * t61;
t322 = t103 * t181 - t177 * t70;
t49 = -qJD(4) * pkin(4) - t322;
t39 = pkin(5) * t127 + t49;
t200 = -t172 * t173 * t171 + t174 * t303;
t87 = t178 * t240 + t199;
t64 = t177 * t200 + t87 * t181;
t327 = t39 * t81 - g(1) * (-t164 * t60 - t165 * t36) - g(2) * (-t164 * t58 - t165 * t34) - g(3) * (-t164 * t86 - t165 * t64) + t17;
t162 = t181 * qJDD(3);
t124 = t177 * t264 + qJDD(5) - t162;
t101 = t151 * t174 - t171 * t243;
t287 = t177 * t101;
t212 = t259 - t293;
t43 = qJDD(3) * pkin(9) + (t151 * t178 + t153 * t276) * t171 + (qJD(1) * qJD(3) * t212 + qJDD(1) * t211) * t172;
t10 = qJDD(4) * pkin(10) + qJD(4) * t322 + t181 * t43 + t287;
t26 = qJD(3) * t138 + qJDD(3) * t142 + t187;
t25 = t180 * t26;
t192 = -qJD(5) * t21 - t176 * t10 + t25;
t2 = t124 * pkin(5) - t76 * pkin(11) + t192;
t268 = qJD(5) * t180;
t261 = t180 * t10 + t176 * t26 + t65 * t268;
t270 = qJD(5) * t176;
t214 = t270 * t50 - t261;
t3 = -pkin(11) * t77 - t214;
t258 = -t175 * t3 + t179 * t2;
t20 = -t176 * t50 + t180 * t65;
t18 = -pkin(11) * t129 + t20;
t16 = -pkin(5) * t154 + t18;
t308 = t179 * t19;
t5 = t16 * t175 + t308;
t326 = t39 * t219 - g(1) * (-t164 * t36 + t165 * t60) - g(2) * (-t164 * t34 + t165 * t58) - g(3) * (-t164 * t64 + t165 * t86) - qJD(6) * t5 + t258;
t191 = qJD(6) * t219 - t175 * t76 - t179 * t77;
t325 = t149 * t219 + t191;
t273 = qJD(4) * t177;
t288 = t176 * t181;
t313 = pkin(9) * t176;
t324 = -t333 * t180 - t273 * t313 - t288 * t74;
t285 = t180 * t181;
t323 = t142 * t268 + t333 * t176 - t285 * t74;
t131 = t175 * t180 + t176 * t179;
t107 = t131 * t177;
t319 = -t176 * t269 + t181 * t265;
t318 = qJD(5) + qJD(6);
t242 = qJD(6) * t16 + t3;
t317 = t175 * t2 + t179 * t242;
t315 = (pkin(9) * t154 + t50) * qJD(5) + t207;
t314 = pkin(10) + pkin(11);
t155 = pkin(9) * t285;
t216 = pkin(5) * t177 - pkin(11) * t285;
t311 = -t216 * qJD(4) - (-t155 + (pkin(11) * t177 - t142) * t176) * qJD(5) + t324;
t272 = qJD(4) * t181;
t248 = t176 * t272;
t202 = t177 * t268 + t248;
t310 = -t202 * pkin(11) + (-t177 * t265 - t181 * t270) * pkin(9) + t323;
t309 = qJD(3) * pkin(3);
t307 = t76 * t176;
t135 = t230 * qJD(3);
t306 = t176 * t135 + t180 * t322;
t130 = t175 * t176 - t179 * t180;
t205 = t130 * t181;
t305 = qJD(3) * t205 - t130 * t318;
t304 = (-t277 + t318) * t131;
t298 = t127 * t154;
t297 = t129 * t154;
t296 = t129 * t180;
t295 = t164 * t181;
t294 = t165 * t181;
t291 = t171 * t182;
t289 = t176 * t177;
t286 = t177 * t180;
t282 = t176 * t142 + t155;
t167 = t177 ^ 2;
t281 = -t181 ^ 2 + t167;
t275 = qJD(4) * t127;
t271 = qJD(5) * t154;
t257 = qJD(5) * t314;
t253 = t171 * t276;
t252 = t176 * t277;
t251 = t154 * t265;
t245 = t182 * t264;
t234 = -t52 + (-t252 + t270) * pkin(5);
t146 = t314 * t176;
t229 = pkin(11) * t252 - qJD(6) * t146 - t176 * t257 - t306;
t118 = t180 * t135;
t147 = t314 * t180;
t228 = qJD(3) * t216 + qJD(6) * t147 - t176 * t322 + t180 * t257 + t118;
t37 = -t176 * t64 + t180 * t86;
t38 = t176 * t86 + t180 * t64;
t225 = -t175 * t38 + t179 * t37;
t224 = t175 * t37 + t179 * t38;
t126 = t180 * t142;
t85 = -pkin(11) * t286 + t126 + (-pkin(5) - t313) * t181;
t96 = -pkin(11) * t289 + t282;
t223 = t175 * t85 + t179 * t96;
t115 = t174 * t177 + t181 * t292;
t213 = -t115 * t180 + t176 * t291;
t94 = -t115 * t176 - t180 * t291;
t222 = t175 * t213 + t179 * t94;
t221 = t175 * t94 - t179 * t213;
t218 = -t101 * t181 + t103 * t273 + t177 * t43 + t70 * t272;
t184 = qJD(3) ^ 2;
t217 = qJDD(3) * t182 - t178 * t184;
t114 = -t174 * t181 + t177 * t292;
t210 = t176 * t124 - t154 * t268;
t209 = t180 * t124 + t154 * t270;
t63 = t87 * t177 - t181 * t200;
t208 = g(1) * (-t177 * t61 + t181 * t89) + g(2) * (-t177 * t59 + t181 * t88) - g(3) * t63;
t206 = g(1) * t61 + g(2) * t59 + g(3) * t87;
t11 = -qJDD(4) * pkin(4) + t218;
t198 = -pkin(10) * t124 - t154 * t49;
t69 = -t74 - t309;
t197 = -pkin(9) * qJDD(4) + (t69 + t74 - t309) * qJD(4);
t190 = pkin(10) * t271 - t11 - t208;
t183 = qJD(4) ^ 2;
t186 = 0.2e1 * qJDD(3) * pkin(3) - pkin(9) * t183 + t331;
t158 = -pkin(5) * t180 - pkin(4);
t139 = (pkin(5) * t176 + pkin(9)) * t177;
t121 = qJDD(6) + t124;
t108 = t130 * t177;
t102 = pkin(5) * t202 + pkin(9) * t272;
t93 = qJD(4) * t115 + t177 * t253;
t92 = -qJD(4) * t114 + t181 * t253;
t79 = t87 * qJD(3);
t78 = (t172 * t212 + t233) * qJD(3);
t56 = -t267 * t289 + (t286 * t318 + t248) * t179 + t319 * t175;
t55 = -qJD(4) * t205 - t107 * t318;
t48 = qJD(5) * t213 - t176 * t92 + t180 * t254;
t47 = qJD(5) * t94 + t176 * t254 + t180 * t92;
t32 = -qJD(4) * t63 + t78 * t181;
t31 = qJD(4) * t64 + t78 * t177;
t8 = pkin(5) * t77 + t11;
t7 = qJD(5) * t37 + t79 * t176 + t32 * t180;
t6 = -qJD(5) * t38 - t32 * t176 + t79 * t180;
t4 = t16 * t179 - t175 * t19;
t1 = [qJDD(1) - g(3), t151 * t303 - g(3) + (t170 ^ 2 + t173 ^ 2) * t172 ^ 2 * qJDD(1), 0, -qJD(3) * t79 - qJDD(3) * t86, -qJD(3) * t78 - qJDD(3) * t87, 0, 0, 0, 0, 0, -t86 * t162 - t31 * qJD(4) - t63 * qJDD(4) + (-t181 * t79 + t273 * t86) * qJD(3), t86 * t262 - t32 * qJD(4) - t64 * qJDD(4) + (t177 * t79 + t272 * t86) * qJD(3), 0, 0, 0, 0, 0, t124 * t37 + t127 * t31 - t154 * t6 + t63 * t77, -t124 * t38 + t129 * t31 + t154 * t7 + t63 * t76, 0, 0, 0, 0, 0 -(-qJD(6) * t224 - t175 * t7 + t179 * t6) * t149 + t225 * t121 + t31 * t81 - t63 * t191 (qJD(6) * t225 + t175 * t6 + t179 * t7) * t149 - t224 * t121 - t31 * t219 + t63 * t22; 0, -g(3) * t303 + (-g(1) * t301 + g(2) * t302) * t172 + t151, 0, t217 * t171 (-qJDD(3) * t178 - t182 * t184) * t171, 0, 0, 0, 0, 0, -t93 * qJD(4) - t114 * qJDD(4) + (-t177 * t245 + t181 * t217) * t171, -t92 * qJD(4) - t115 * qJDD(4) + (-t177 * t217 - t181 * t245) * t171, 0, 0, 0, 0, 0, t114 * t77 + t124 * t94 + t127 * t93 - t154 * t48, t114 * t76 + t124 * t213 + t129 * t93 + t154 * t47, 0, 0, 0, 0, 0 -(-qJD(6) * t221 - t175 * t47 + t179 * t48) * t149 + t222 * t121 + t93 * t81 - t114 * t191 (qJD(6) * t222 + t175 * t48 + t179 * t47) * t149 - t221 * t121 - t93 * t219 + t114 * t22; 0, 0, qJDD(3), t331, -t151 * t292 - t211 * t263 + (-t153 * t291 - t212 * t280 + t74) * qJD(3) + t206, qJDD(3) * t167 + 0.2e1 * t177 * t246, 0.2e1 * t162 * t177 - 0.2e1 * t264 * t281, qJDD(4) * t177 + t181 * t183, qJDD(4) * t181 - t177 * t183, 0, t177 * t197 + t181 * t186, -t177 * t186 + t181 * t197, t129 * t319 + t76 * t286 (-t127 * t180 - t129 * t176) * t272 + (-t307 - t180 * t77 + (t127 * t176 - t296) * qJD(5)) * t177 (-t76 - t251) * t181 + (qJD(4) * t129 + t209) * t177 (t154 * t274 + t77) * t181 + (-t210 - t275) * t177, -t124 * t181 - t154 * t273, t126 * t124 + t324 * t154 + (t142 * t271 - t206) * t176 + (pkin(9) * t275 - t25 + (-pkin(9) * t124 + qJD(4) * t49 + qJD(5) * t65 + t10) * t176 + t315 * t180) * t181 + (pkin(9) * t77 + qJD(4) * t20 + t11 * t176 - t127 * t74 + t268 * t49) * t177, -t282 * t124 + t323 * t154 - t206 * t180 + ((pkin(9) * t129 + t180 * t49) * qJD(4) - t315 * t176 + t261) * t181 + (-t49 * t270 - t21 * qJD(4) + t11 * t180 - t74 * t129 + (t76 - t251) * pkin(9)) * t177, -t108 * t22 - t219 * t55, -t107 * t22 - t108 * t191 + t219 * t56 - t55 * t81, -t108 * t121 - t149 * t55 - t181 * t22 - t219 * t273, -t107 * t121 + t149 * t56 - t181 * t191 - t273 * t81, -t121 * t181 - t149 * t273 (-t175 * t96 + t179 * t85) * t121 - t258 * t181 + t102 * t81 - t139 * t191 + t8 * t107 + t39 * t56 - g(1) * (t164 * t61 - t294 * t60) - g(2) * (t164 * t59 - t294 * t58) - g(3) * (t164 * t87 - t294 * t86) + (qJD(4) * t4 - t74 * t81) * t177 + (t310 * t175 + t311 * t179) * t149 + (t149 * t223 + t181 * t5) * qJD(6), -t223 * t121 + (-t17 + t317) * t181 - t102 * t219 + t139 * t22 - t8 * t108 + t39 * t55 - g(1) * (t165 * t61 + t295 * t60) - g(2) * (t165 * t59 + t295 * t58) - g(3) * (t165 * t87 + t295 * t86) + (-qJD(4) * t5 + t219 * t74) * t177 + ((qJD(6) * t85 + t310) * t179 + (-qJD(6) * t96 - t311) * t175) * t149; 0, 0, 0, 0, 0, -t181 * t184 * t177, t281 * t184, t262, t162, qJDD(4), qJD(4) * t52 - t279 * t69 - t208 - t218, g(1) * t36 + g(2) * t34 + g(3) * t64 - t287 + (-qJD(3) * t69 - t43) * t181, -t154 * t296 + t307 (t76 + t298) * t180 + (-t77 + t297) * t176 (-t129 * t177 + t154 * t285) * qJD(3) + t210 (t127 * t177 - t154 * t288) * qJD(3) + t209, t154 * t279, -t20 * t279 - pkin(4) * t77 + t118 * t154 - t52 * t127 + (-t154 * t322 + t198) * t176 + t190 * t180, -pkin(4) * t76 - t52 * t129 - t154 * t306 - t190 * t176 + t198 * t180 + t21 * t279, t22 * t131 - t219 * t305, -t22 * t130 + t131 * t191 + t219 * t304 - t305 * t81, t131 * t121 - t149 * t305 + t219 * t279, -t130 * t121 + t149 * t304 + t279 * t81, t149 * t279 (-t146 * t179 - t147 * t175) * t121 - t158 * t191 + t8 * t130 - t4 * t279 + t234 * t81 + t304 * t39 + (t175 * t229 + t179 * t228) * t149 - t208 * t165 -(-t146 * t175 + t147 * t179) * t121 + t158 * t22 + t8 * t131 + t5 * t279 - t234 * t219 + t305 * t39 + (-t175 * t228 + t179 * t229) * t149 + t208 * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129 * t127, -t127 ^ 2 + t129 ^ 2, t76 - t298, -t297 - t77, t124, -t21 * t154 - t49 * t129 - g(1) * (-t176 * t36 + t180 * t60) - g(2) * (-t176 * t34 + t180 * t58) - g(3) * t37 + t192, -t20 * t154 + t49 * t127 - g(1) * (-t176 * t60 - t180 * t36) - g(2) * (-t176 * t58 - t180 * t34) + g(3) * t38 + t214, -t332, t329, t328, t325, t121 (-t175 * t18 - t308) * t149 + (t121 * t179 - t129 * t81 + t149 * t267) * pkin(5) + t326 (t149 * t19 - t2) * t175 + (-t149 * t18 - t242) * t179 + (-t121 * t175 + t129 * t219 + t149 * t266) * pkin(5) + t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t332, t329, t328, t325, t121, -t5 * t149 + t326, -t4 * t149 - t317 + t327;];
tau_reg  = t1;
