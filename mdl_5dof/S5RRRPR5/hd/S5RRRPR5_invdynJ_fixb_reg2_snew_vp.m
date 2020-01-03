% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPR5
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:20
% EndTime: 2019-12-31 21:14:33
% DurationCPUTime: 5.84s
% Computational Cost: add. (28166->446), mult. (63816->641), div. (0->0), fcn. (46619->10), ass. (0->279)
t267 = cos(qJ(2));
t250 = t267 * qJDD(1);
t263 = sin(qJ(2));
t301 = qJD(1) * qJD(2);
t295 = t263 * t301;
t238 = t250 - t295;
t258 = t267 ^ 2;
t269 = qJD(1) ^ 2;
t264 = sin(qJ(1));
t335 = cos(qJ(1));
t293 = t264 * g(1) - t335 * g(2);
t275 = qJDD(1) * pkin(1) + t293;
t307 = qJD(1) * t263;
t276 = qJD(2) * pkin(2) - pkin(7) * t307;
t196 = t238 * pkin(2) - t276 * t307 + t275 + (pkin(7) * t258 + pkin(6)) * t269;
t346 = 2 * qJD(4);
t259 = sin(pkin(9));
t262 = sin(qJ(3));
t266 = cos(qJ(3));
t230 = t266 * t267 * qJD(1) - t262 * t307;
t231 = (t267 * t262 + t263 * t266) * qJD(1);
t260 = cos(pkin(9));
t209 = -t260 * t230 + t259 * t231;
t211 = t259 * t230 + t260 * t231;
t175 = t211 * t209;
t255 = qJDD(2) + qJDD(3);
t338 = -t175 + t255;
t345 = t259 * t338;
t344 = t260 * t338;
t261 = sin(qJ(5));
t249 = t263 * qJDD(1);
t294 = t267 * t301;
t237 = t249 + t294;
t285 = t262 * t237 - t266 * t238;
t193 = -t231 * qJD(3) - t285;
t279 = t266 * t237 + t262 * t238;
t194 = t230 * qJD(3) + t279;
t287 = -t260 * t193 + t259 * t194;
t157 = qJDD(5) + t287;
t256 = qJD(2) + qJD(3);
t265 = cos(qJ(5));
t190 = t261 * t211 - t265 * t256;
t192 = t265 * t211 + t261 * t256;
t161 = t192 * t190;
t339 = t157 - t161;
t343 = t261 * t339;
t217 = t230 * t231;
t337 = t217 + t255;
t342 = t262 * t337;
t341 = t265 * t339;
t340 = t266 * t337;
t228 = t230 ^ 2;
t281 = t256 * pkin(3) - t231 * qJ(4);
t129 = t193 * pkin(3) + t228 * qJ(4) - t231 * t281 - qJDD(4) + t196;
t328 = t256 * t211;
t136 = t287 + t328;
t171 = t209 * pkin(4) - t211 * pkin(8);
t336 = t256 ^ 2;
t313 = t263 * t269;
t277 = t335 * g(1) + t264 * g(2);
t331 = qJDD(1) * pkin(6);
t233 = -t269 * pkin(1) - t277 + t331;
t315 = t263 * t233;
t185 = qJDD(2) * pkin(2) - t237 * pkin(7) - t315 + (pkin(2) * t313 + pkin(7) * t301 - g(3)) * t267;
t221 = -t263 * g(3) + t267 * t233;
t252 = t258 * t269;
t186 = -pkin(2) * t252 + t238 * pkin(7) - qJD(2) * t276 + t221;
t156 = t262 * t185 + t266 * t186;
t118 = -t228 * pkin(3) + t193 * qJ(4) - t256 * t281 + t156;
t155 = -t266 * t185 + t262 * t186;
t224 = t256 * t230;
t181 = t194 - t224;
t271 = pkin(3) * t337 - qJ(4) * t181 - t155;
t70 = -0.2e1 * qJD(4) * t209 + t260 * t118 + t259 * t271;
t52 = -pkin(4) * t336 + t255 * pkin(8) - t209 * t171 + t70;
t159 = t259 * t193 + t260 * t194;
t200 = t256 * t209;
t140 = t159 - t200;
t72 = t136 * pkin(4) - pkin(8) * t140 - t129;
t33 = t261 * t52 - t265 * t72;
t34 = t261 * t72 + t265 * t52;
t17 = t261 * t33 + t265 * t34;
t290 = t259 * t118 - t260 * t271;
t69 = t211 * t346 + t290;
t205 = qJD(5) + t209;
t288 = t261 * t159 - t265 * t255;
t107 = (qJD(5) - t205) * t192 + t288;
t188 = t190 ^ 2;
t189 = t192 ^ 2;
t204 = t205 ^ 2;
t206 = t209 ^ 2;
t207 = t211 ^ 2;
t229 = t231 ^ 2;
t334 = pkin(4) * t259;
t51 = -t255 * pkin(4) - t336 * pkin(8) + (t346 + t171) * t211 + t290;
t48 = t261 * t51;
t36 = t259 * t70 - t260 * t69;
t333 = t262 * t36;
t49 = t265 * t51;
t332 = t266 * t36;
t330 = t205 * t261;
t329 = t205 * t265;
t327 = t256 * t259;
t326 = t256 * t260;
t325 = t256 * t262;
t324 = t256 * t266;
t323 = t259 * t129;
t169 = t175 + t255;
t322 = t259 * t169;
t321 = t260 * t129;
t320 = t260 * t169;
t120 = t157 + t161;
t319 = t261 * t120;
t318 = t262 * t196;
t214 = -t217 + t255;
t317 = t262 * t214;
t113 = -t266 * t155 + t262 * t156;
t316 = t263 * t113;
t244 = t267 * t313;
t314 = t263 * (qJDD(2) + t244);
t312 = t265 * t120;
t311 = t266 * t196;
t310 = t266 * t214;
t309 = t267 * (qJDD(2) - t244);
t304 = qJD(3) + t256;
t302 = qJD(5) + t205;
t12 = t259 * t17 - t260 * t51;
t300 = pkin(3) * t12 - pkin(4) * t51 + pkin(8) * t17;
t299 = t259 * t161;
t298 = t260 * t161;
t297 = -pkin(4) * t260 - pkin(3);
t37 = t259 * t69 + t260 * t70;
t280 = -t265 * t159 - t261 * t255;
t112 = t302 * t190 + t280;
t152 = -t189 - t204;
t85 = -t261 * t152 - t312;
t56 = t260 * t112 + t259 * t85;
t292 = pkin(3) * t56 + pkin(4) * t112 + pkin(8) * t85 + t48;
t108 = -t302 * t192 - t288;
t144 = -t204 - t188;
t80 = t265 * t144 - t343;
t46 = t260 * t108 + t259 * t80;
t291 = pkin(3) * t46 + pkin(4) * t108 + pkin(8) * t80 - t49;
t114 = t262 * t155 + t266 * t156;
t220 = t267 * g(3) + t315;
t286 = t263 * t220 + t267 * t221;
t197 = -t207 - t336;
t146 = t260 * t197 - t322;
t284 = pkin(3) * t146 - t70;
t135 = t188 + t189;
t128 = -t190 * qJD(5) - t280;
t166 = t205 * t190;
t111 = t128 + t166;
t66 = -t107 * t265 + t261 * t111;
t41 = t260 * t135 + t259 * t66;
t283 = pkin(3) * t41 + pkin(4) * t135 + pkin(8) * t66 + t17;
t16 = t261 * t34 - t265 * t33;
t167 = -t336 - t206;
t125 = t259 * t167 + t344;
t278 = pkin(3) * t125 - t69;
t274 = -t287 + t328;
t273 = (-qJD(3) + t256) * t231 - t285;
t268 = qJD(2) ^ 2;
t257 = t263 ^ 2;
t251 = t257 * t269;
t239 = t250 - 0.2e1 * t295;
t236 = t249 + 0.2e1 * t294;
t232 = t269 * pkin(6) + t275;
t223 = -t229 + t336;
t222 = t228 - t336;
t219 = -t229 - t336;
t216 = t229 - t228;
t212 = -t336 - t228;
t199 = -t207 + t336;
t198 = t206 - t336;
t195 = -t228 - t229;
t183 = -t262 * t219 - t310;
t182 = t266 * t219 - t317;
t180 = t194 + t224;
t179 = t304 * t230 + t279;
t176 = t304 * t231 + t285;
t174 = t266 * t212 - t342;
t173 = t262 * t212 + t340;
t172 = t207 - t206;
t165 = -t189 + t204;
t164 = t188 - t204;
t163 = (-t209 * t260 + t211 * t259) * t256;
t162 = (-t209 * t259 - t211 * t260) * t256;
t160 = t189 - t188;
t154 = -t206 - t207;
t151 = t260 * t198 - t322;
t150 = -t259 * t199 + t344;
t149 = t259 * t198 + t320;
t148 = t260 * t199 + t345;
t147 = -t259 * t197 - t320;
t143 = t262 * t181 + t266 * t273;
t142 = -t266 * t181 + t262 * t273;
t141 = t159 + t200;
t133 = t260 * t159 - t211 * t327;
t132 = t259 * t159 + t211 * t326;
t131 = t209 * t326 + t259 * t287;
t130 = t209 * t327 - t260 * t287;
t127 = -t192 * qJD(5) - t288;
t126 = t260 * t167 - t345;
t123 = (-t190 * t265 + t192 * t261) * t205;
t122 = (-t190 * t261 - t192 * t265) * t205;
t110 = t128 - t166;
t104 = t265 * t128 - t192 * t330;
t103 = t261 * t128 + t192 * t329;
t102 = -t261 * t127 + t190 * t329;
t101 = t265 * t127 + t190 * t330;
t100 = -t262 * t146 + t266 * t147;
t99 = t266 * t146 + t262 * t147;
t98 = -qJ(4) * t146 - t321;
t97 = t260 * t123 + t259 * t157;
t96 = t259 * t123 - t260 * t157;
t95 = t259 * t141 + t260 * t274;
t94 = -t260 * t136 - t259 * t140;
t93 = -t260 * t141 + t259 * t274;
t92 = -t259 * t136 + t260 * t140;
t91 = pkin(3) * t93;
t90 = t265 * t164 - t319;
t89 = -t261 * t165 + t341;
t88 = t261 * t164 + t312;
t87 = t265 * t165 + t343;
t86 = -qJ(4) * t125 - t323;
t84 = t265 * t152 - t319;
t82 = -t262 * t125 + t266 * t126;
t81 = t266 * t125 + t262 * t126;
t79 = t261 * t144 + t341;
t77 = t260 * t104 + t299;
t76 = t260 * t102 - t299;
t75 = t259 * t104 - t298;
t74 = t259 * t102 + t298;
t73 = -pkin(3) * t140 + qJ(4) * t147 - t323;
t68 = -pkin(3) * t136 + qJ(4) * t126 + t321;
t65 = t265 * t108 - t261 * t110;
t64 = -t107 * t261 - t265 * t111;
t63 = t261 * t108 + t265 * t110;
t61 = -t259 * t107 + t260 * t90;
t60 = t259 * t111 + t260 * t89;
t59 = t260 * t107 + t259 * t90;
t58 = -t260 * t111 + t259 * t89;
t57 = -t259 * t112 + t260 * t85;
t54 = -t262 * t93 + t266 * t95;
t53 = t262 * t95 + t266 * t93;
t47 = -t259 * t108 + t260 * t80;
t44 = t259 * t160 + t260 * t65;
t43 = -t260 * t160 + t259 * t65;
t42 = -t259 * t135 + t260 * t66;
t39 = -pkin(8) * t84 + t49;
t38 = -pkin(8) * t79 + t48;
t35 = pkin(3) * t36;
t30 = pkin(3) * t129 + qJ(4) * t37;
t29 = -t262 * t56 + t266 * t57;
t28 = t262 * t57 + t266 * t56;
t27 = -t262 * t46 + t266 * t47;
t26 = t262 * t47 + t266 * t46;
t25 = -qJ(4) * t93 - t36;
t24 = -t262 * t41 + t266 * t42;
t23 = t262 * t42 + t266 * t41;
t22 = -pkin(3) * t154 + qJ(4) * t95 + t37;
t21 = -pkin(4) * t84 + t34;
t20 = -pkin(4) * t79 + t33;
t19 = t266 * t37 - t333;
t18 = t262 * t37 + t332;
t14 = -pkin(8) * t64 - t16;
t13 = t260 * t17 + t259 * t51;
t10 = -qJ(4) * t56 - t259 * t21 + t260 * t39;
t9 = -qJ(4) * t46 - t259 * t20 + t260 * t38;
t8 = -pkin(3) * t84 + qJ(4) * t57 + t260 * t21 + t259 * t39;
t7 = -pkin(3) * t79 + qJ(4) * t47 + t260 * t20 + t259 * t38;
t6 = -qJ(4) * t41 + t260 * t14 + t64 * t334;
t5 = qJ(4) * t42 + t259 * t14 + t297 * t64;
t4 = -t262 * t12 + t266 * t13;
t3 = t266 * t12 + t262 * t13;
t2 = -qJ(4) * t12 + (-pkin(8) * t260 + t334) * t16;
t1 = qJ(4) * t13 + (-pkin(8) * t259 + t297) * t16;
t11 = [0, 0, 0, 0, 0, qJDD(1), t293, t277, 0, 0, (t237 + t294) * t263, t267 * t236 + t263 * t239, t314 + t267 * (-t251 + t268), (t238 - t295) * t267, t263 * (t252 - t268) + t309, 0, t267 * t232 + pkin(1) * t239 + pkin(6) * (t267 * (-t252 - t268) - t314), -t263 * t232 - pkin(1) * t236 + pkin(6) * (-t309 - t263 * (-t251 - t268)), pkin(1) * (t251 + t252) + (t257 + t258) * t331 + t286, pkin(1) * t232 + pkin(6) * t286, t263 * (t266 * t194 - t231 * t325) + t267 * (t262 * t194 + t231 * t324), t263 * (-t266 * t176 - t262 * t180) + t267 * (-t262 * t176 + t266 * t180), t263 * (-t262 * t223 + t340) + t267 * (t266 * t223 + t342), t263 * (-t262 * t193 - t230 * t324) + t267 * (t266 * t193 - t230 * t325), t263 * (t266 * t222 - t317) + t267 * (t262 * t222 + t310), (t263 * (t230 * t266 + t231 * t262) + t267 * (t230 * t262 - t231 * t266)) * t256, t263 * (-pkin(7) * t173 - t318) + t267 * (-pkin(2) * t176 + pkin(7) * t174 + t311) - pkin(1) * t176 + pkin(6) * (-t263 * t173 + t267 * t174), t263 * (-pkin(7) * t182 - t311) + t267 * (-pkin(2) * t179 + pkin(7) * t183 - t318) - pkin(1) * t179 + pkin(6) * (-t263 * t182 + t267 * t183), t263 * (-pkin(7) * t142 - t113) + t267 * (-pkin(2) * t195 + pkin(7) * t143 + t114) - pkin(1) * t195 + pkin(6) * (-t263 * t142 + t267 * t143), -pkin(7) * t316 + t267 * (pkin(2) * t196 + pkin(7) * t114) + pkin(1) * t196 + pkin(6) * (t267 * t114 - t316), t263 * (-t262 * t132 + t266 * t133) + t267 * (t266 * t132 + t262 * t133), t263 * (-t262 * t92 + t266 * t94) + t267 * (t262 * t94 + t266 * t92), t263 * (-t262 * t148 + t266 * t150) + t267 * (t266 * t148 + t262 * t150), t263 * (-t262 * t130 + t266 * t131) + t267 * (t266 * t130 + t262 * t131), t263 * (-t262 * t149 + t266 * t151) + t267 * (t266 * t149 + t262 * t151), t263 * (-t262 * t162 + t266 * t163) + t267 * (t266 * t162 + t262 * t163), t263 * (-pkin(7) * t81 - t262 * t68 + t266 * t86) + t267 * (-pkin(2) * t136 + pkin(7) * t82 + t262 * t86 + t266 * t68) - pkin(1) * t136 + pkin(6) * (-t263 * t81 + t267 * t82), t263 * (-pkin(7) * t99 - t262 * t73 + t266 * t98) + t267 * (-pkin(2) * t140 + pkin(7) * t100 + t262 * t98 + t266 * t73) - pkin(1) * t140 + pkin(6) * (t267 * t100 - t263 * t99), t263 * (-pkin(7) * t53 - t262 * t22 + t266 * t25) + t267 * (-pkin(2) * t154 + pkin(7) * t54 + t266 * t22 + t262 * t25) - pkin(1) * t154 + pkin(6) * (-t263 * t53 + t267 * t54), t263 * (-pkin(7) * t18 - qJ(4) * t332 - t262 * t30) + t267 * (pkin(2) * t129 + pkin(7) * t19 - qJ(4) * t333 + t266 * t30) + pkin(1) * t129 + pkin(6) * (-t263 * t18 + t267 * t19), t263 * (-t262 * t75 + t266 * t77) + t267 * (t262 * t77 + t266 * t75), t263 * (-t262 * t43 + t266 * t44) + t267 * (t262 * t44 + t266 * t43), t263 * (-t262 * t58 + t266 * t60) + t267 * (t262 * t60 + t266 * t58), t263 * (-t262 * t74 + t266 * t76) + t267 * (t262 * t76 + t266 * t74), t263 * (-t262 * t59 + t266 * t61) + t267 * (t262 * t61 + t266 * t59), t263 * (-t262 * t96 + t266 * t97) + t267 * (t262 * t97 + t266 * t96), t263 * (-pkin(7) * t26 - t262 * t7 + t266 * t9) + t267 * (-pkin(2) * t79 + pkin(7) * t27 + t262 * t9 + t266 * t7) - pkin(1) * t79 + pkin(6) * (-t263 * t26 + t267 * t27), t263 * (-pkin(7) * t28 + t266 * t10 - t262 * t8) + t267 * (-pkin(2) * t84 + pkin(7) * t29 + t262 * t10 + t266 * t8) - pkin(1) * t84 + pkin(6) * (-t263 * t28 + t267 * t29), t263 * (-pkin(7) * t23 - t262 * t5 + t266 * t6) + t267 * (-pkin(2) * t64 + pkin(7) * t24 + t262 * t6 + t266 * t5) - pkin(1) * t64 + pkin(6) * (-t263 * t23 + t267 * t24), t263 * (-pkin(7) * t3 - t262 * t1 + t266 * t2) + t267 * (-pkin(2) * t16 + pkin(7) * t4 + t266 * t1 + t262 * t2) - pkin(1) * t16 + pkin(6) * (-t263 * t3 + t267 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, t251 - t252, t249, t244, t250, qJDD(2), -t220, -t221, 0, 0, -t217, t216, t181, t217, t273, t255, pkin(2) * t173 - t155, pkin(2) * t182 - t156, pkin(2) * t142, pkin(2) * t113, t175, t172, t141, -t175, t274, t255, pkin(2) * t81 + t278, pkin(2) * t99 + t284, pkin(2) * t53 + t91, pkin(2) * t18 + t35, t103, t63, t87, t101, t88, t122, pkin(2) * t26 + t291, pkin(2) * t28 + t292, pkin(2) * t23 + t283, pkin(2) * t3 + t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, t216, t181, t217, t273, t255, -t155, -t156, 0, 0, t175, t172, t141, -t175, t274, t255, t278, t284, t91, t35, t103, t63, t87, t101, t88, t122, t291, t292, t283, t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t140, t154, -t129, 0, 0, 0, 0, 0, 0, t79, t84, t64, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, t160, t111, -t161, -t107, t157, -t33, -t34, 0, 0;];
tauJ_reg = t11;
