% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:08
% EndTime: 2019-03-08 21:21:19
% DurationCPUTime: 4.59s
% Computational Cost: add. (3050->474), mult. (6941->657), div. (0->0), fcn. (5284->14), ass. (0->241)
t179 = sin(pkin(6));
t185 = sin(qJ(2));
t188 = cos(qJ(2));
t275 = qJD(1) * qJD(2);
t253 = t188 * t275;
t181 = cos(pkin(6));
t281 = qJD(3) * t181;
t313 = qJDD(2) * pkin(8);
t344 = qJD(1) * t281 + t313 + (qJDD(1) * t185 + t253) * t179;
t290 = qJD(1) * t179;
t263 = t185 * t290;
t134 = qJD(2) * pkin(8) + t263;
t184 = sin(qJ(3));
t187 = cos(qJ(3));
t289 = qJD(1) * t181;
t293 = -t134 * t184 + t187 * t289;
t340 = qJD(4) - t293;
t286 = qJD(2) * t184;
t160 = qJD(6) + t286;
t343 = t160 - qJD(6);
t330 = pkin(4) + pkin(8);
t248 = -qJ(4) * t184 - pkin(2);
t177 = sin(pkin(11));
t180 = cos(pkin(11));
t183 = sin(qJ(6));
t186 = cos(qJ(6));
t125 = t177 * t186 + t180 * t183;
t109 = t125 * qJD(6);
t205 = t125 * t184;
t321 = -qJD(2) * t205 - t109;
t284 = qJD(2) * t187;
t257 = t180 * t284;
t283 = qJD(3) * t177;
t120 = t257 + t283;
t255 = t177 * t284;
t282 = qJD(3) * t180;
t122 = -t255 + t282;
t222 = t120 * t183 - t122 * t186;
t342 = t160 * t222;
t68 = -qJD(3) * pkin(3) + t340;
t274 = qJD(2) * qJD(3);
t251 = t184 * t274;
t270 = t187 * qJDD(2);
t341 = -t251 + t270;
t302 = t179 * t187;
t114 = t181 * t184 + t185 * t302;
t190 = qJD(2) ^ 2;
t217 = qJDD(2) * t188 - t185 * t190;
t250 = t187 * t274;
t301 = t179 * t188;
t259 = qJD(2) * t301;
t303 = t179 * t185;
t267 = t184 * t303;
t61 = -qJD(3) * t267 + (t259 + t281) * t187;
t339 = qJD(3) * t61 + qJDD(3) * t114 + t179 * (t184 * t217 + t188 * t250);
t113 = -t181 * t187 + t267;
t62 = qJD(3) * t114 + t184 * t259;
t338 = -qJD(3) * t62 - qJDD(3) * t113 + t179 * (t187 * t217 - t188 * t251);
t279 = qJD(3) * t187;
t131 = t330 * t279;
t299 = t184 * t188;
t280 = qJD(3) * t184;
t165 = pkin(3) * t280;
t316 = qJ(4) * t187;
t223 = qJ(5) * t184 - t316;
t278 = qJD(4) * t184;
t195 = qJD(3) * t223 - qJD(5) * t187 - t278;
t66 = t165 + t195;
t325 = t180 * t131 - t177 * t66 - (-t177 * t185 + t180 * t299) * t290;
t324 = t177 * t131 + t180 * t66 - (t177 * t299 + t180 * t185) * t290;
t174 = qJD(3) * qJ(4);
t80 = t134 * t187 + t184 * t289;
t70 = -t174 - t80;
t296 = pkin(4) * t286 + t340;
t288 = qJD(1) * t188;
t262 = t179 * t288;
t136 = -pkin(3) * t187 + t248;
t287 = qJD(2) * t136;
t81 = -t262 + t287;
t337 = t286 * t81 + qJDD(4);
t336 = t134 * t279 + t184 * t344;
t254 = t185 * t275;
t228 = -qJDD(1) * t301 + t179 * t254;
t319 = cos(pkin(10));
t242 = t319 * t188;
t178 = sin(pkin(10));
t305 = t178 * t185;
t105 = -t181 * t242 + t305;
t243 = t319 * t185;
t304 = t178 * t188;
t107 = t181 * t304 + t243;
t234 = g(1) * t107 + g(2) * t105;
t189 = qJD(3) ^ 2;
t329 = pkin(8) * t189;
t335 = 0.2e1 * qJDD(2) * pkin(2) + t179 * (-g(3) * t188 + t254) - t228 + t234 - t329;
t182 = -pkin(3) - qJ(5);
t65 = pkin(4) * t284 + t80;
t46 = qJD(5) + t174 + t65;
t333 = -t182 * t279 + (qJD(5) - t46) * t184;
t209 = -qJ(4) * t279 - t278;
t102 = t165 + t209;
t201 = g(3) * t301 - t234;
t272 = qJDD(2) * t136;
t214 = pkin(3) * t251 + t228;
t29 = qJD(2) * t209 + t214 + t272;
t331 = qJD(2) * (-t102 + t263) - t201 - t272 - t29 - t329;
t82 = qJDD(3) * t177 + t180 * t341;
t238 = qJDD(3) * t180 - t177 * t270;
t83 = t177 * t251 + t238;
t11 = -qJD(6) * t222 + t183 * t83 + t186 * t82;
t328 = -pkin(9) + t182;
t273 = qJDD(1) * t181;
t203 = -t187 * t273 + t336;
t200 = qJDD(4) + t203;
t271 = t184 * qJDD(2);
t204 = t250 + t271;
t15 = pkin(4) * t204 - qJD(3) * qJD(5) + qJDD(3) * t182 + t200;
t117 = t182 * t187 + t248;
t20 = qJD(2) * t195 + qJDD(2) * t117 + t214;
t7 = t15 * t177 + t180 * t20;
t215 = -pkin(9) * t177 * t184 + pkin(5) * t187;
t327 = qJD(3) * t215 + t325;
t326 = -pkin(9) * t180 * t280 - t324;
t42 = qJD(3) * t182 + t296;
t63 = qJD(2) * t117 - t262;
t19 = t177 * t42 + t180 * t63;
t166 = pkin(3) * t286;
t97 = qJD(2) * t223 + t166;
t28 = t177 * t65 + t180 * t97;
t323 = qJD(2) * pkin(2);
t43 = t120 * t186 + t122 * t183;
t322 = t160 * t43;
t258 = t180 * t286;
t276 = qJD(6) * t186;
t277 = qJD(6) * t183;
t306 = t177 * t183;
t320 = -t177 * t277 + t180 * t276 + t186 * t258 - t286 * t306;
t145 = t330 * t184;
t51 = t117 * t180 + t145 * t177;
t318 = pkin(8) * qJDD(3);
t315 = qJD(3) * t80;
t312 = qJDD(3) * pkin(3);
t311 = t105 * t187;
t310 = t107 * t187;
t171 = pkin(11) + qJ(6);
t167 = sin(t171);
t309 = t167 * t184;
t168 = cos(t171);
t308 = t168 * t184;
t175 = t184 ^ 2;
t307 = t175 * t190;
t300 = t180 * t187;
t298 = t187 * t188;
t265 = -pkin(5) * t180 - pkin(4);
t297 = -t265 * t286 + t340;
t294 = qJDD(1) - g(3);
t146 = t330 * t187;
t176 = t187 ^ 2;
t292 = t175 - t176;
t291 = t175 + t176;
t285 = qJD(2) * t185;
t269 = -pkin(3) * t311 + t105 * t248;
t268 = -pkin(3) * t310 + t107 * t248;
t266 = t184 * t190 * t187;
t6 = t15 * t180 - t177 * t20;
t2 = pkin(5) * t204 - pkin(9) * t83 + t6;
t5 = -pkin(9) * t82 + t7;
t264 = -t183 * t5 + t186 * t2;
t261 = t187 * t288;
t260 = t179 * t285;
t249 = t182 * t271;
t106 = t181 * t243 + t304;
t244 = t179 * t319;
t57 = t106 * t184 + t187 * t244;
t58 = t106 * t187 - t184 * t244;
t247 = -pkin(3) * t57 + qJ(4) * t58;
t108 = -t181 * t305 + t242;
t59 = t108 * t184 - t178 * t302;
t60 = t178 * t179 * t184 + t108 * t187;
t246 = -pkin(3) * t59 + qJ(4) * t60;
t18 = -t177 * t63 + t180 * t42;
t27 = -t177 * t97 + t180 * t65;
t240 = -pkin(3) * t113 + qJ(4) * t114;
t239 = -t134 * t280 + t184 * t273 + t187 * t344;
t124 = qJDD(6) + t204;
t221 = -t180 * t186 + t306;
t237 = -t124 * t221 + t160 * t321;
t233 = g(1) * t108 + g(2) * t106;
t232 = t177 * t7 + t180 * t6;
t231 = t183 * t2 + t186 * t5;
t12 = -pkin(9) * t120 + t19;
t9 = pkin(5) * t286 - pkin(9) * t122 + t18;
t4 = t12 * t186 + t183 * t9;
t229 = t12 * t183 - t186 * t9;
t128 = t180 * t145;
t35 = pkin(5) * t184 + t128 + (pkin(9) * t187 - t117) * t177;
t40 = -pkin(9) * t300 + t51;
t227 = -t183 * t40 + t186 * t35;
t226 = t183 * t35 + t186 * t40;
t55 = t113 * t180 + t177 * t301;
t56 = t113 * t177 - t180 * t301;
t225 = -t183 * t56 + t186 * t55;
t224 = t183 * t55 + t186 * t56;
t220 = g(3) * (pkin(2) * t301 + pkin(8) * t303 + (pkin(3) * t298 + qJ(4) * t299) * t179);
t172 = qJDD(3) * qJ(4);
t173 = qJD(3) * qJD(4);
t22 = -t172 - t173 - t239;
t132 = t328 * t177;
t212 = qJD(2) * t215 + qJD(5) * t180 + qJD(6) * t132 + t27;
t133 = t328 * t180;
t211 = pkin(9) * t258 + qJD(5) * t177 - qJD(6) * t133 + t28;
t10 = -t120 * t276 - t122 * t277 - t183 * t82 + t186 * t83;
t210 = (-t177 * t18 + t180 * t19) * t184;
t208 = -t124 * t125 - t160 * t320;
t207 = g(1) * t59 + g(2) * t57 + g(3) * t113;
t206 = -g(1) * t60 - g(2) * t58 - g(3) * t114;
t94 = t221 * t187;
t16 = pkin(4) * t341 + qJDD(5) - t22;
t202 = t16 + t206;
t199 = g(3) * t303 - t16 * t187 + t233;
t135 = -t262 - t323;
t197 = -t318 + (t135 + t262 - t323) * qJD(3);
t196 = t318 + (-t262 - t81 - t287) * qJD(3);
t194 = -qJD(3) * t293 + t206 + t239;
t193 = -t203 + t207;
t24 = t200 - t312;
t191 = t24 * t184 - t22 * t187 + (t184 * t70 + t187 * t68) * qJD(3) - t233;
t161 = pkin(5) * t177 + qJ(4);
t130 = t330 * t280;
t129 = -qJ(4) * t284 + t166;
t104 = pkin(5) * t300 + t146;
t95 = t125 * t187;
t91 = (-pkin(8) + t265) * t280;
t50 = -t117 * t177 + t128;
t37 = t109 * t187 - t221 * t280;
t36 = qJD(3) * t205 + qJD(6) * t94;
t34 = t177 * t62 + t180 * t260;
t33 = -t177 * t260 + t180 * t62;
t30 = pkin(5) * t120 + t46;
t8 = pkin(5) * t82 + t16;
t1 = [t294, 0, t217 * t179 (-qJDD(2) * t185 - t188 * t190) * t179, 0, 0, 0, 0, 0, t338, -t339 (t113 * t184 + t114 * t187) * qJDD(2) + (t184 * t62 + t187 * t61 + (t113 * t187 - t114 * t184) * qJD(3)) * qJD(2), -t338, t339, t113 * t24 - t114 * t22 - t61 * t70 + t62 * t68 - g(3) + (-t188 * t29 + t285 * t81) * t179, t55 * t271 + t114 * t82 + t120 * t61 + (t184 * t33 + t279 * t55) * qJD(2), -t56 * t271 + t114 * t83 + t122 * t61 + (-t184 * t34 - t279 * t56) * qJD(2), -t120 * t34 - t122 * t33 - t55 * t83 - t56 * t82, t114 * t16 + t18 * t33 + t19 * t34 + t46 * t61 + t55 * t6 + t56 * t7 - g(3), 0, 0, 0, 0, 0 (-qJD(6) * t224 - t183 * t34 + t186 * t33) * t160 + t225 * t124 + t61 * t43 + t114 * t11 -(qJD(6) * t225 + t183 * t33 + t186 * t34) * t160 - t224 * t124 - t61 * t222 + t114 * t10; 0, qJDD(2), t294 * t301 + t234, -t294 * t303 + t233, qJDD(2) * t175 + 0.2e1 * t184 * t250, 0.2e1 * t184 * t270 - 0.2e1 * t274 * t292, qJDD(3) * t184 + t187 * t189, qJDD(3) * t187 - t184 * t189, 0, t184 * t197 + t187 * t335, -t184 * t335 + t187 * t197, t291 * t313 + (-g(3) * t185 - t253 * t291) * t179 + t191, t184 * t196 - t187 * t331, t184 * t331 + t187 * t196, t29 * t136 + t81 * t102 - g(1) * t268 - g(2) * t269 - t220 + (-t185 * t81 + (-t184 * t68 + t187 * t70) * t188) * t290 + t191 * pkin(8), -t130 * t120 + t146 * t82 + (-t120 * t262 + (qJD(2) * t50 + t18) * qJD(3)) * t187 - t199 * t180 + (qJD(2) * t325 + qJDD(2) * t50 - t177 * t201 - t282 * t46 + t6) * t184, -t130 * t122 + t146 * t83 + (-t122 * t262 + (-qJD(2) * t51 - t19) * qJD(3)) * t187 + t199 * t177 + (-qJD(2) * t324 - t51 * qJDD(2) - t180 * t201 + t283 * t46 - t7) * t184, -t50 * t83 - t51 * t82 - t325 * t122 - t324 * t120 + qJD(3) * t210 + (t177 * t6 - t180 * t7 - t201) * t187, t7 * t51 + t6 * t50 + t16 * t146 - t46 * t130 - g(1) * (-qJ(5) * t310 + t108 * t330 + t268) - g(2) * (-qJ(5) * t311 + t106 * t330 + t269) - t220 + t324 * t19 + t325 * t18 + (-t46 * t261 - g(3) * (pkin(4) * t185 + qJ(5) * t298)) * t179, -t10 * t95 - t222 * t36, t10 * t94 + t11 * t95 - t222 * t37 - t36 * t43, t10 * t184 - t124 * t95 + t160 * t36 - t222 * t279, -t11 * t184 + t124 * t94 + t160 * t37 - t279 * t43, t124 * t184 + t160 * t279, t227 * t124 + t264 * t184 - t229 * t279 + t91 * t43 + t104 * t11 - t8 * t94 - t30 * t37 - g(1) * (-t107 * t309 + t108 * t168) - g(2) * (-t105 * t309 + t106 * t168) + (t183 * t326 + t186 * t327) * t160 + (-t160 * t226 - t184 * t4) * qJD(6) + (-t43 * t261 - g(3) * (t167 * t299 + t168 * t185)) * t179, -t226 * t124 - t231 * t184 - t4 * t279 - t91 * t222 + t104 * t10 - t8 * t95 + t30 * t36 - g(1) * (-t107 * t308 - t108 * t167) - g(2) * (-t105 * t308 - t106 * t167) + (-t183 * t327 + t186 * t326) * t160 + (-t160 * t227 + t184 * t229) * qJD(6) + (t222 * t261 - g(3) * (-t167 * t185 + t168 * t299)) * t179; 0, 0, 0, 0, -t266, t292 * t190, t271, t270, qJDD(3), -t135 * t286 + t193 + t315, -t135 * t284 - t194 (-pkin(3) * t184 + t316) * qJDD(2), -0.2e1 * t312 - t315 + (-qJD(2) * t129 - t273) * t187 - t207 + t336 + t337, 0.2e1 * t172 + 0.2e1 * t173 + (t129 * t184 + t187 * t81) * qJD(2) + t194, -t24 * pkin(3) - g(1) * t246 - g(2) * t247 - g(3) * t240 - t22 * qJ(4) - t81 * t129 - t340 * t70 - t68 * t80, t180 * t249 + qJ(4) * t82 + t296 * t120 + t202 * t177 + (-t18 * t187 - t180 * t333 - t184 * t27) * qJD(2), -t177 * t249 + qJ(4) * t83 + t296 * t122 + t202 * t180 + (t177 * t333 + t184 * t28 + t19 * t187) * qJD(2), t120 * t28 + t122 * t27 + (qJD(5) * t122 - t182 * t83 - t19 * t286 - t6) * t180 + (qJD(5) * t120 + t18 * t286 - t182 * t82 - t7) * t177 + t207, t16 * qJ(4) - t19 * t28 - t18 * t27 - g(1) * (-qJ(5) * t59 + t246) - g(2) * (-qJ(5) * t57 + t247) - g(3) * (-qJ(5) * t113 + t240) + t296 * t46 + t232 * t182 + (-t177 * t19 - t18 * t180) * qJD(5), -t10 * t221 - t222 * t321, -t10 * t125 + t11 * t221 + t222 * t320 - t321 * t43, t222 * t284 + t237, t284 * t43 + t208, -t160 * t284 (-t132 * t183 + t133 * t186) * t124 + t161 * t11 + t8 * t125 + t229 * t284 + t297 * t43 + t320 * t30 + (t183 * t211 - t186 * t212) * t160 + t206 * t167 -(t132 * t186 + t133 * t183) * t124 + t161 * t10 - t8 * t221 + t4 * t284 - t297 * t222 + t321 * t30 + (t183 * t212 + t186 * t211) * t160 + t206 * t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, qJDD(3) + t266, -t189 - t307, qJD(3) * t70 - t193 - t312 + t337, t180 * t271 - t177 * t307 + (-t120 + t257) * qJD(3), -t177 * t271 - t180 * t307 + (-t122 - t255) * qJD(3), -t177 * t82 - t180 * t83 + (-t120 * t180 + t122 * t177) * t286, qJD(2) * t210 - qJD(3) * t46 - t207 + t232, 0, 0, 0, 0, 0, -qJD(3) * t43 + t237, qJD(3) * t222 + t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122 * t286 + t82 (-t120 + t283) * t286 + t238, -t120 ^ 2 - t122 ^ 2, t120 * t19 + t122 * t18 + t202, 0, 0, 0, 0, 0, t11 - t342, t10 - t322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222 * t43, t222 ^ 2 - t43 ^ 2, t10 + t322, -t11 - t342, t124, t30 * t222 - g(1) * (-t107 * t167 + t168 * t59) - g(2) * (-t105 * t167 + t168 * t57) - g(3) * (t113 * t168 + t167 * t301) + t264 + t343 * t4, t30 * t43 - g(1) * (-t107 * t168 - t167 * t59) - g(2) * (-t105 * t168 - t167 * t57) - g(3) * (-t113 * t167 + t168 * t301) - t231 - t343 * t229;];
tau_reg  = t1;
