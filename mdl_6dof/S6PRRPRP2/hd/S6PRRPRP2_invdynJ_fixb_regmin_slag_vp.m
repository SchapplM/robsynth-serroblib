% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:57
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:54:48
% EndTime: 2021-01-16 02:55:07
% DurationCPUTime: 5.42s
% Computational Cost: add. (5744->513), mult. (13399->673), div. (0->0), fcn. (10635->14), ass. (0->253)
t191 = cos(qJ(3));
t307 = cos(pkin(11));
t248 = t307 * t191;
t182 = sin(pkin(11));
t188 = sin(qJ(3));
t297 = t182 * t188;
t219 = t248 - t297;
t192 = cos(qJ(2));
t184 = sin(pkin(6));
t281 = qJD(1) * t184;
t262 = t192 * t281;
t105 = t219 * t262;
t186 = qJ(4) + pkin(8);
t254 = qJD(3) * t186;
t130 = qJD(4) * t191 - t188 * t254;
t213 = -qJD(4) * t188 - t191 * t254;
t69 = t130 * t307 + t182 * t213;
t339 = t69 - t105;
t163 = qJD(2) * t248;
t279 = qJD(2) * t188;
t136 = t182 * t279 - t163;
t129 = qJD(5) + t136;
t185 = cos(pkin(6));
t189 = sin(qJ(2));
t293 = t184 * t189;
t142 = t185 * t191 - t188 * t293;
t187 = sin(qJ(5));
t190 = cos(qJ(5));
t237 = pkin(5) * t187 - qJ(6) * t190;
t263 = t189 * t281;
t244 = qJD(2) * t186 + t263;
t280 = qJD(1) * t185;
t102 = -t188 * t244 + t191 * t280;
t103 = t188 * t280 + t191 * t244;
t250 = t307 * t103;
t47 = t102 * t182 + t250;
t311 = -qJD(6) * t187 + t129 * t237 - t47;
t154 = t186 * t191;
t107 = t154 * t307 - t186 * t297;
t275 = qJD(5) * t190;
t276 = qJD(5) * t187;
t249 = t307 * t188;
t146 = t182 * t191 + t249;
t138 = t146 * qJD(3);
t141 = t219 * qJD(3);
t318 = qJD(3) * pkin(3);
t267 = t188 * t318;
t71 = pkin(4) * t138 - pkin(9) * t141 + t267;
t173 = t191 * pkin(3) + pkin(2);
t84 = -pkin(4) * t219 - pkin(9) * t146 - t173;
t340 = t107 * t276 - t84 * t275 - t339 * t190 + (t263 - t71) * t187;
t309 = t130 * t182 - t146 * t262 - t213 * t307;
t179 = qJ(3) + pkin(11);
t174 = sin(t179);
t308 = cos(pkin(10));
t251 = t308 * t192;
t183 = sin(pkin(10));
t295 = t183 * t189;
t133 = -t185 * t251 + t295;
t252 = t308 * t189;
t294 = t183 * t192;
t135 = t185 * t294 + t252;
t239 = g(1) * t135 + g(2) * t133;
t291 = t184 * t192;
t211 = -g(3) * t291 + t239;
t338 = t211 * t174;
t160 = t187 * t291;
t292 = t184 * t191;
t143 = t185 * t188 + t189 * t292;
t76 = t142 * t182 + t143 * t307;
t60 = t190 * t76 - t160;
t273 = qJD(2) * qJD(3);
t258 = t188 * t273;
t206 = qJDD(2) * t146 - t182 * t258;
t337 = qJD(3) * t163 + t206;
t336 = t263 - t267;
t139 = t146 * qJD(2);
t335 = t139 * qJD(3);
t334 = pkin(3) * t258 + qJDD(4);
t108 = -qJD(3) * t190 + t139 * t187;
t110 = qJD(3) * t187 + t139 * t190;
t94 = t182 * t103;
t98 = t102 + t318;
t44 = t307 * t98 - t94;
t39 = -qJD(3) * pkin(4) - t44;
t21 = pkin(5) * t108 - qJ(6) * t110 + t39;
t326 = pkin(3) * t182;
t170 = pkin(9) + t326;
t270 = t188 * qJDD(2);
t235 = -qJDD(2) * t248 + t182 * t270;
t86 = qJD(2) * t138 + t235;
t82 = qJDD(5) + t86;
t315 = t170 * t82;
t333 = t129 * t21 - t315;
t256 = qJDD(1) * t291;
t274 = qJD(1) * qJD(2);
t259 = t189 * t274;
t236 = t184 * t259 - t256;
t306 = qJDD(2) * pkin(2);
t117 = t236 - t306;
t193 = qJD(3) ^ 2;
t332 = -pkin(8) * t193 + t184 * (-g(3) * t192 + t259) - t117 + t239 + t306;
t331 = t110 ^ 2;
t330 = t129 ^ 2;
t329 = pkin(5) * t82;
t328 = qJ(6) * t138 - qJD(6) * t219 - t340;
t321 = t107 * t190 + t187 * t84;
t72 = t105 * t187 - t190 * t263;
t327 = -pkin(5) * t138 + qJD(5) * t321 + t187 * t69 - t190 * t71 - t72;
t325 = pkin(3) * t188;
t271 = t185 * qJDD(1);
t162 = t191 * t271;
t118 = qJDD(2) * pkin(8) + (qJDD(1) * t189 + t192 * t274) * t184;
t203 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t280 + t118;
t233 = t244 * qJD(3);
t35 = qJDD(3) * pkin(3) - t188 * t203 - t191 * t233 + t162;
t36 = (-t233 + t271) * t188 + t203 * t191;
t10 = -t182 * t36 + t307 * t35;
t11 = t182 * t35 + t307 * t36;
t49 = t102 * t307 - t94;
t268 = pkin(3) * t279;
t70 = pkin(4) * t139 + pkin(9) * t136 + t268;
t322 = t187 * t70 + t190 * t49;
t45 = t182 * t98 + t250;
t320 = qJ(6) * t82;
t319 = qJD(2) * pkin(2);
t40 = qJD(3) * pkin(9) + t45;
t128 = -qJD(2) * t173 + qJD(4) - t262;
t58 = pkin(4) * t136 - pkin(9) * t139 + t128;
t19 = t187 * t58 + t190 * t40;
t13 = qJ(6) * t129 + t19;
t317 = t129 * t13;
t316 = t129 * t19;
t272 = qJD(3) * qJD(5);
t41 = -t187 * qJDD(3) + t139 * t276 + (-t337 - t272) * t190;
t314 = t187 * t41;
t77 = t187 * t82;
t78 = t190 * t82;
t195 = -qJDD(3) * t190 + t187 * t337;
t42 = qJD(5) * t110 + t195;
t312 = -t108 * t275 - t187 * t42;
t238 = pkin(5) * t190 + qJ(6) * t187;
t310 = t237 * t141 + (qJD(5) * t238 - qJD(6) * t190) * t146 + t309;
t305 = t108 * t136;
t304 = t108 * t139;
t303 = t110 * t108;
t246 = t110 * t129;
t302 = t110 * t139;
t301 = t129 * t187;
t300 = t146 * t190;
t175 = cos(t179);
t299 = t175 * t187;
t298 = t175 * t190;
t296 = t183 * t184;
t289 = t187 * t141;
t288 = t190 * t141;
t287 = t190 * t192;
t18 = -t187 * t40 + t190 * t58;
t286 = qJD(6) - t18;
t285 = qJDD(1) - g(3);
t132 = t185 * t295 - t251;
t284 = -t132 * t186 - t135 * t173;
t283 = t173 * t291 + t186 * t293;
t180 = t188 ^ 2;
t282 = -t191 ^ 2 + t180;
t278 = qJD(2) * t189;
t277 = qJD(5) * t170;
t269 = t191 * qJDD(2);
t265 = t184 * t287;
t264 = t307 * pkin(3);
t261 = t184 * t278;
t260 = qJD(2) * t291;
t257 = t191 * t273;
t26 = -pkin(3) * t269 + t86 * pkin(4) - pkin(9) * t337 + t117 + t334;
t9 = qJDD(3) * pkin(9) + t11;
t255 = t187 * t9 - t190 * t26 + t275 * t40 + t276 * t58;
t253 = t184 * t308;
t134 = t185 * t252 + t294;
t247 = -t133 * t173 + t134 * t186;
t87 = t132 * t175 - t174 * t296;
t106 = t154 * t182 + t186 * t249;
t243 = t188 * t260;
t242 = pkin(3) * t183 * t292 + t132 * t325;
t8 = -qJDD(3) * pkin(4) - t10;
t241 = pkin(4) * t175 + pkin(9) * t174;
t240 = g(1) * t132 - g(2) * t134;
t234 = t142 * pkin(3);
t194 = qJD(2) ^ 2;
t232 = qJDD(2) * t192 - t189 * t194;
t231 = -t129 * t276 - t136 * t301 + t78;
t230 = t77 + (t136 * t190 + t275) * t129;
t229 = pkin(4) + t238;
t226 = -g(1) * t183 + g(2) * t308;
t59 = t187 * t76 + t265;
t225 = t187 * t26 + t190 * t9 + t275 * t58 - t276 * t40;
t223 = -t146 * t276 + t288;
t101 = qJD(3) * t142 + t191 * t260;
t212 = t143 * qJD(3);
t199 = -t212 - t243;
t48 = t101 * t307 + t182 * t199;
t16 = qJD(5) * t60 + t187 * t48 - t190 * t261;
t46 = t101 * t182 - t199 * t307;
t75 = -t142 * t307 + t143 * t182;
t221 = t108 * t46 - t129 * t16 + t42 * t75 - t59 * t82;
t220 = t129 * t39 - t315;
t112 = t160 * t175 - t190 * t293;
t64 = -t133 * t299 - t134 * t190;
t66 = t132 * t190 - t135 * t299;
t218 = g(1) * t66 + g(2) * t64 + g(3) * t112;
t113 = (t175 * t287 + t187 * t189) * t184;
t65 = -t133 * t298 + t134 * t187;
t67 = -t132 * t187 - t135 * t298;
t217 = -g(1) * t67 - g(2) * t65 - g(3) * t113;
t120 = -t174 * t293 + t175 * t185;
t88 = -t134 * t174 - t175 * t253;
t90 = t132 * t174 + t175 * t296;
t216 = -g(1) * t90 - g(2) * t88 - g(3) * t120;
t121 = t174 * t185 + t175 * t293;
t89 = t134 * t175 - t174 * t253;
t215 = g(1) * t87 - g(2) * t89 - g(3) * t121;
t214 = t232 * t184;
t152 = -t262 - t319;
t210 = -qJD(2) * t152 - t118 - t240;
t209 = (-t134 * t188 - t191 * t253) * pkin(3);
t15 = -qJD(5) * t59 + t187 * t261 + t190 * t48;
t208 = t110 * t46 - t129 * t15 - t41 * t75 - t60 * t82;
t53 = -t133 * t190 + t187 * t89;
t55 = -t135 * t190 - t187 * t87;
t92 = t121 * t187 + t265;
t207 = g(1) * t55 + g(2) * t53 + g(3) * t92 - t255;
t205 = -t129 * t277 + t216;
t3 = pkin(5) * t42 + qJ(6) * t41 - qJD(6) * t110 + t8;
t202 = t205 - t3;
t201 = -pkin(8) * qJDD(3) + (t152 + t262 - t319) * qJD(3);
t54 = t133 * t187 + t190 * t89;
t56 = t135 * t187 - t190 * t87;
t93 = t121 * t190 - t160;
t200 = -g(1) * t56 - g(2) * t54 - g(3) * t93 + t225;
t85 = -qJDD(2) * t173 + t236 + t334;
t198 = t110 * t21 + qJDD(6) - t207;
t171 = -t264 - pkin(4);
t144 = -t264 - t229;
t57 = pkin(5) * t110 + qJ(6) * t108;
t50 = t146 * t237 + t106;
t30 = pkin(5) * t219 + t107 * t187 - t190 * t84;
t29 = -qJ(6) * t219 + t321;
t22 = t108 * t129 - t41;
t20 = -pkin(5) * t139 + t187 * t49 - t190 * t70;
t17 = qJ(6) * t139 + t322;
t12 = -pkin(5) * t129 + t286;
t2 = qJDD(6) + t255 - t329;
t1 = qJD(6) * t129 + t225 + t320;
t4 = [t285, 0, t214, (-qJDD(2) * t189 - t192 * t194) * t184, 0, 0, 0, 0, 0, t142 * qJDD(3) + t191 * t214 + (-t212 - 0.2e1 * t243) * qJD(3), -qJD(3) * t101 - qJDD(3) * t143 + (-t188 * t232 - t192 * t257) * t184, -qJD(3) * t46 - qJDD(3) * t75 + (t136 * t278 - t192 * t86) * t184, -t48 * qJD(3) - t76 * qJDD(3) + (t139 * t278 - t192 * t337) * t184, -t48 * t136 + t46 * t139 + t337 * t75 - t76 * t86, -t10 * t75 + t11 * t76 - t44 * t46 + t45 * t48 - g(3) + (t128 * t278 - t192 * t85) * t184, 0, 0, 0, 0, 0, t221, t208, t221, -t108 * t15 + t110 * t16 - t41 * t59 - t42 * t60, -t208, t1 * t60 + t12 * t16 + t13 * t15 + t2 * t59 + t21 * t46 + t3 * t75 - g(3); 0, qJDD(2), t211 + t256, -t285 * t293 - t240, qJDD(2) * t180 + 0.2e1 * t188 * t257, 0.2e1 * t188 * t269 - 0.2e1 * t273 * t282, qJDD(3) * t188 + t191 * t193, qJDD(3) * t191 - t188 * t193, 0, t188 * t201 + t191 * t332, -t188 * t332 + t191 * t201, -t136 * t263 - qJDD(3) * t106 + t128 * t138 - t219 * t85 - t173 * t86 + t211 * t175 + (t136 * t325 - t309) * qJD(3), -qJD(3) * t339 - t107 * qJDD(3) + t128 * t141 - t139 * t336 + t85 * t146 - t173 * t337 - t338, -g(3) * t293 - t10 * t146 + t106 * t337 - t107 * t86 + t11 * t219 - t136 * t339 - t45 * t138 + t139 * t309 - t44 * t141 + t240, -g(1) * t284 - g(2) * t247 - g(3) * t283 - t10 * t106 + t11 * t107 - t128 * t336 - t85 * t173 - t309 * t44 + t339 * t45, t110 * t223 - t300 * t41, (-t108 * t190 - t110 * t187) * t141 + (t314 - t190 * t42 + (t108 * t187 - t110 * t190) * qJD(5)) * t146, t110 * t138 + t129 * t223 + t219 * t41 + t300 * t82, -t146 * t77 - t108 * t138 + t219 * t42 + (-t146 * t275 - t289) * t129, t129 * t138 - t219 * t82, t255 * t219 + t18 * t138 + t106 * t42 + t72 * t129 + t309 * t108 + ((-qJD(5) * t107 + t71) * t129 + t84 * t82 + t39 * qJD(5) * t146) * t190 + ((-qJD(5) * t84 - t69) * t129 - t107 * t82 + t8 * t146 + t39 * t141) * t187 + t217, -t321 * t82 + t225 * t219 - t19 * t138 - t106 * t41 + t39 * t288 + (t8 * t190 - t276 * t39) * t146 + t340 * t129 + t309 * t110 + t218, t21 * t289 - t12 * t138 + t219 * t2 - t30 * t82 + t42 * t50 + (t187 * t3 + t21 * t275) * t146 - t327 * t129 + t310 * t108 + t217, -t29 * t42 - t30 * t41 + (t12 * t190 - t13 * t187) * t141 + t327 * t110 - t328 * t108 + t338 + (-t1 * t187 + t190 * t2 + (-t12 * t187 - t13 * t190) * qJD(5)) * t146, -t21 * t288 - t1 * t219 + t13 * t138 + t29 * t82 + t41 * t50 + (-t190 * t3 + t21 * t276) * t146 + t328 * t129 - t310 * t110 - t218, t1 * t29 + t3 * t50 + t2 * t30 - g(1) * (pkin(5) * t67 + qJ(6) * t66 - t135 * t241 + t284) - g(2) * (pkin(5) * t65 + qJ(6) * t64 - t133 * t241 + t247) - g(3) * (pkin(5) * t113 + qJ(6) * t112 + t241 * t291 + t283) + t310 * t21 + t328 * t13 + t327 * t12; 0, 0, 0, 0, -t188 * t194 * t191, t282 * t194, t270, t269, qJDD(3), -g(3) * t142 + t188 * t210 + t226 * t292 + t162, g(3) * t143 + (-t184 * t226 - t271) * t188 + t210 * t191, t47 * qJD(3) - t128 * t139 + (qJDD(3) * t307 - t136 * t279) * pkin(3) + t216 + t10, qJD(3) * t49 + t128 * t136 + (-qJDD(3) * t182 - t139 * t279) * pkin(3) - t215 - t11, -t337 * t264 - t86 * t326 - (-t45 + t47) * t139 + (-t44 + t49) * t136, -g(1) * t242 - g(2) * t209 - g(3) * t234 + t10 * t264 + t11 * t326 - t128 * t268 + t44 * t47 - t45 * t49, t190 * t246 - t314, (-t41 - t305) * t190 - t110 * t301 + t312, t230 - t302, t231 + t304, -t129 * t139, -t47 * t108 - t18 * t139 + t171 * t42 + (t49 * t129 + t220) * t187 + (-t8 + (-t70 - t277) * t129 + t216) * t190, -t171 * t41 + t322 * t129 + t19 * t139 - t47 * t110 + t220 * t190 + (-t205 + t8) * t187, t311 * t108 + t12 * t139 + t129 * t20 + t144 * t42 + t187 * t333 + t202 * t190, t108 * t17 - t110 * t20 + (t12 * t136 - t170 * t42 + t1 + (t110 * t170 + t12) * qJD(5)) * t190 + (-t13 * t136 - t170 * t41 + t2 + (t108 * t170 - t13) * qJD(5)) * t187 + t215, -t311 * t110 - t129 * t17 - t13 * t139 + t144 * t41 + t202 * t187 - t190 * t333, t3 * t144 - t13 * t17 - t12 * t20 - g(1) * (-pkin(9) * t87 + t229 * t90 + t242) - g(2) * (t89 * pkin(9) + t229 * t88 + t209) - g(3) * (pkin(9) * t121 + t120 * t229 + t234) + t311 * t21 + (t1 * t190 + t12 * t275 - t13 * t276 + t187 * t2) * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235 + 0.2e1 * t335, (t163 - t136) * qJD(3) + t206, -t136 ^ 2 - t139 ^ 2, t136 * t45 + t139 * t44 - t211 + t85, 0, 0, 0, 0, 0, t231 - t304, -t190 * t330 - t302 - t77, -t129 * t301 - t304 + t78, (t41 - t305) * t190 + t187 * t246 + t312, t230 + t302, -t139 * t21 + (-t2 + t317) * t190 + (t12 * t129 + t1) * t187 - t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t303, -t108 ^ 2 + t331, t22, -t139 * t275 - t187 * t272 - t195 + t246, t82, -t110 * t39 + t207 + t316, t108 * t39 + t129 * t18 - t200, -t108 * t57 - t198 + t316 + 0.2e1 * t329, pkin(5) * t41 - qJ(6) * t42 + (t13 - t19) * t110 + (t12 - t286) * t108, 0.2e1 * t320 - t108 * t21 + t110 * t57 + (0.2e1 * qJD(6) - t18) * t129 + t200, t1 * qJ(6) - t2 * pkin(5) - t21 * t57 - t12 * t19 - g(1) * (-pkin(5) * t55 + qJ(6) * t56) - g(2) * (-pkin(5) * t53 + qJ(6) * t54) - g(3) * (-pkin(5) * t92 + qJ(6) * t93) + t286 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(5) - t235 + t303 - t335, t22, -t330 - t331, t198 - t317 - t329;];
tau_reg = t4;
