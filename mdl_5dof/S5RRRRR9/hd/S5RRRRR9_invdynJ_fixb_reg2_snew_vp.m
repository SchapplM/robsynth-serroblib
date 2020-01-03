% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:56
% EndTime: 2019-12-31 22:30:10
% DurationCPUTime: 6.01s
% Computational Cost: add. (38737->459), mult. (79167->631), div. (0->0), fcn. (57035->10), ass. (0->285)
t251 = sin(qJ(5));
t253 = sin(qJ(3));
t258 = cos(qJ(3));
t254 = sin(qJ(2));
t296 = qJD(1) * t254;
t221 = -qJD(2) * t258 + t253 * t296;
t222 = qJD(2) * t253 + t258 * t296;
t252 = sin(qJ(4));
t257 = cos(qJ(4));
t201 = t221 * t257 + t222 * t252;
t203 = -t221 * t252 + t222 * t257;
t256 = cos(qJ(5));
t166 = t201 * t256 + t203 * t251;
t168 = -t201 * t251 + t203 * t256;
t132 = t168 * t166;
t291 = qJD(1) * qJD(2);
t240 = t254 * t291;
t259 = cos(qJ(2));
t290 = t259 * qJDD(1);
t226 = -t240 + t290;
t220 = -qJDD(3) + t226;
t217 = -qJDD(4) + t220;
t213 = -qJDD(5) + t217;
t331 = -t132 - t213;
t338 = t251 * t331;
t170 = t203 * t201;
t330 = -t170 - t217;
t337 = t252 * t330;
t336 = t256 * t331;
t335 = t257 * t330;
t261 = qJD(1) ^ 2;
t238 = qJD(1) * t259 - qJD(3);
t234 = -qJD(4) + t238;
t229 = -qJD(5) + t234;
t154 = t166 * t229;
t241 = t254 * qJDD(1);
t284 = t259 * t291;
t225 = t241 + t284;
t269 = -t253 * qJDD(2) - t258 * t225;
t196 = -qJD(3) * t221 - t269;
t279 = -qJDD(2) * t258 + t253 * t225;
t266 = -qJD(3) * t222 - t279;
t144 = -qJD(4) * t201 + t196 * t257 + t252 * t266;
t281 = t252 * t196 - t257 * t266;
t265 = qJD(4) * t203 + t281;
t93 = -qJD(5) * t166 + t144 * t256 - t251 * t265;
t334 = t154 + t93;
t315 = t222 * t221;
t263 = -t220 - t315;
t333 = t253 * t263;
t332 = t258 * t263;
t187 = t201 * t234;
t129 = -t187 + t144;
t329 = t187 + t144;
t212 = t221 * t238;
t178 = t196 - t212;
t282 = t251 * t144 + t256 * t265;
t75 = (qJD(5) + t229) * t168 + t282;
t126 = (qJD(4) + t234) * t203 + t281;
t174 = (qJD(3) + t238) * t222 + t279;
t164 = t166 ^ 2;
t165 = t168 ^ 2;
t199 = t201 ^ 2;
t200 = t203 ^ 2;
t328 = t221 ^ 2;
t219 = t222 ^ 2;
t228 = t229 ^ 2;
t233 = t234 ^ 2;
t236 = t238 ^ 2;
t327 = qJD(2) ^ 2;
t255 = sin(qJ(1));
t260 = cos(qJ(1));
t283 = t255 * g(1) - g(2) * t260;
t215 = qJDD(1) * pkin(1) + t261 * pkin(6) + t283;
t272 = -t226 + t240;
t273 = t225 + t284;
t173 = pkin(2) * t272 - pkin(7) * t273 - t215;
t274 = g(1) * t260 + g(2) * t255;
t316 = qJDD(1) * pkin(6);
t216 = -pkin(1) * t261 - t274 + t316;
t275 = -pkin(2) * t259 - pkin(7) * t254;
t278 = t261 * t275 + t216;
t324 = t254 * g(3);
t183 = -pkin(2) * t327 + qJDD(2) * pkin(7) + t259 * t278 - t324;
t146 = -t173 * t258 + t253 * t183;
t110 = t263 * pkin(3) - pkin(8) * t178 - t146;
t147 = t173 * t253 + t183 * t258;
t209 = -pkin(3) * t238 - pkin(8) * t222;
t112 = -pkin(3) * t328 + pkin(8) * t266 + t209 * t238 + t147;
t68 = -t110 * t257 + t252 * t112;
t69 = t110 * t252 + t112 * t257;
t39 = t252 * t69 - t257 * t68;
t326 = pkin(3) * t39;
t89 = -t126 * t252 - t129 * t257;
t325 = pkin(3) * t89;
t323 = t259 * g(3);
t182 = -qJDD(2) * pkin(2) - pkin(7) * t327 + t254 * t278 + t323;
t135 = -pkin(3) * t266 - pkin(8) * t328 + t222 * t209 + t182;
t184 = -pkin(4) * t234 - pkin(9) * t203;
t80 = pkin(4) * t265 - pkin(9) * t199 + t203 * t184 + t135;
t322 = t251 * t80;
t54 = pkin(4) * t330 - pkin(9) * t129 - t68;
t58 = -pkin(4) * t199 - pkin(9) * t265 + t184 * t234 + t69;
t30 = t251 * t58 - t256 * t54;
t31 = t251 * t54 + t256 * t58;
t17 = t251 * t31 - t256 * t30;
t321 = t252 * t17;
t320 = t253 * t39;
t319 = t256 * t80;
t318 = t257 * t17;
t317 = t258 * t39;
t314 = t229 * t251;
t313 = t229 * t256;
t312 = t234 * t252;
t311 = t234 * t257;
t310 = t238 * t253;
t309 = t238 * t258;
t117 = -t132 + t213;
t308 = t251 * t117;
t307 = t252 * t135;
t156 = -t170 + t217;
t306 = t252 * t156;
t305 = t253 * t182;
t190 = t220 - t315;
t304 = t253 * t190;
t237 = t259 * t261 * t254;
t303 = t254 * (qJDD(2) + t237);
t302 = t256 * t117;
t301 = t257 * t135;
t300 = t257 * t156;
t299 = t258 * t182;
t298 = t258 * t190;
t297 = t259 * (-t237 + qJDD(2));
t295 = qJD(3) - t238;
t16 = pkin(4) * t17;
t18 = t251 * t30 + t256 * t31;
t7 = t18 * t252 + t318;
t289 = pkin(3) * t7 + t16;
t288 = t259 * t132;
t287 = t259 * t170;
t286 = t259 * t315;
t78 = -t154 + t93;
t47 = -t251 * t75 - t256 * t78;
t49 = t251 * t78 - t256 * t75;
t24 = t252 * t49 + t257 * t47;
t45 = pkin(4) * t47;
t285 = pkin(3) * t24 + t45;
t40 = t252 * t68 + t257 * t69;
t106 = t146 * t253 + t147 * t258;
t207 = t216 * t254 + t323;
t208 = t216 * t259 - t324;
t280 = t254 * t207 + t208 * t259;
t124 = -t228 - t164;
t86 = t124 * t251 + t336;
t277 = pkin(4) * t86 - t30;
t180 = -t200 - t233;
t133 = t180 * t257 + t306;
t276 = pkin(3) * t133 - t69;
t271 = t146 * t258 - t147 * t253;
t270 = -pkin(1) + t275;
t149 = -t165 - t228;
t99 = t149 * t256 + t308;
t268 = pkin(4) * t99 - t31;
t161 = -t233 - t199;
t115 = t161 * t252 + t335;
t267 = pkin(3) * t115 - t68;
t87 = t124 * t256 - t338;
t51 = t252 * t87 + t257 * t86;
t264 = pkin(3) * t51 + t277;
t100 = -t149 * t251 + t302;
t59 = t100 * t252 + t257 * t99;
t262 = pkin(3) * t59 + t268;
t248 = t259 ^ 2;
t247 = t254 ^ 2;
t245 = t248 * t261;
t243 = t247 * t261;
t227 = -0.2e1 * t240 + t290;
t224 = t241 + 0.2e1 * t284;
t211 = -t219 + t236;
t210 = -t236 + t328;
t205 = t219 - t328;
t204 = -t219 - t236;
t197 = -t236 - t328;
t189 = t219 + t328;
t186 = -t200 + t233;
t185 = t199 - t233;
t179 = t221 * t295 + t269;
t177 = t196 + t212;
t175 = -t222 * t295 - t279;
t169 = t200 - t199;
t163 = -t204 * t253 + t298;
t162 = t204 * t258 + t304;
t160 = t197 * t258 - t333;
t159 = t197 * t253 + t332;
t153 = -t165 + t228;
t152 = t164 - t228;
t151 = (t201 * t257 - t203 * t252) * t234;
t150 = (t201 * t252 + t203 * t257) * t234;
t148 = -t199 - t200;
t141 = -t174 * t258 + t178 * t253;
t139 = t185 * t257 + t306;
t138 = -t186 * t252 + t335;
t137 = t185 * t252 - t300;
t136 = t186 * t257 + t337;
t134 = -t180 * t252 + t300;
t131 = t164 - t165;
t125 = (qJD(4) - t234) * t203 + t281;
t123 = t144 * t257 + t203 * t312;
t122 = t144 * t252 - t203 * t311;
t121 = -t201 * t311 + t252 * t265;
t120 = -t201 * t312 - t257 * t265;
t116 = t161 * t257 - t337;
t114 = (t166 * t256 - t168 * t251) * t229;
t113 = (t166 * t251 + t168 * t256) * t229;
t107 = -t164 - t165;
t104 = t152 * t256 + t308;
t103 = -t153 * t251 + t336;
t102 = t152 * t251 - t302;
t101 = t153 * t256 + t338;
t97 = -pkin(8) * t133 + t301;
t96 = -t133 * t253 + t134 * t258;
t95 = t133 * t258 + t134 * t253;
t94 = -pkin(8) * t115 + t307;
t92 = -qJD(5) * t168 - t282;
t91 = -t126 * t257 + t129 * t252;
t90 = -t125 * t257 - t252 * t329;
t88 = -t125 * t252 + t257 * t329;
t84 = -t115 * t253 + t116 * t258;
t83 = t115 * t258 + t116 * t253;
t82 = -t113 * t252 + t114 * t257;
t81 = t113 * t257 + t114 * t252;
t74 = (qJD(5) - t229) * t168 + t282;
t73 = t168 * t314 + t256 * t93;
t72 = -t168 * t313 + t251 * t93;
t71 = -t166 * t313 - t251 * t92;
t70 = -t166 * t314 + t256 * t92;
t66 = -pkin(3) * t329 + pkin(8) * t134 + t307;
t65 = -pkin(3) * t125 + pkin(8) * t116 - t301;
t64 = -t102 * t252 + t104 * t257;
t63 = -t101 * t252 + t103 * t257;
t62 = t102 * t257 + t104 * t252;
t61 = t101 * t257 + t103 * t252;
t60 = t100 * t257 - t252 * t99;
t57 = -pkin(9) * t99 + t319;
t56 = -t253 * t89 + t258 * t91;
t55 = t253 * t91 + t258 * t89;
t52 = -t252 * t86 + t257 * t87;
t50 = -pkin(9) * t86 + t322;
t48 = -t251 * t334 - t256 * t74;
t46 = -t251 * t74 + t256 * t334;
t44 = -t252 * t72 + t257 * t73;
t43 = -t252 * t70 + t257 * t71;
t42 = t252 * t73 + t257 * t72;
t41 = t252 * t71 + t257 * t70;
t38 = -pkin(3) * t135 + pkin(8) * t40;
t37 = -pkin(4) * t334 + pkin(9) * t100 + t322;
t36 = -t253 * t59 + t258 * t60;
t35 = t253 * t60 + t258 * t59;
t34 = -pkin(4) * t74 + pkin(9) * t87 - t319;
t33 = -pkin(8) * t89 - t39;
t32 = -pkin(3) * t148 + pkin(8) * t91 + t40;
t28 = -t253 * t51 + t258 * t52;
t27 = t253 * t52 + t258 * t51;
t26 = -t252 * t47 + t257 * t49;
t25 = -t252 * t46 + t257 * t48;
t23 = t252 * t48 + t257 * t46;
t22 = t258 * t40 - t320;
t21 = t253 * t40 + t317;
t20 = -pkin(8) * t59 - t252 * t37 + t257 * t57;
t19 = -pkin(8) * t51 - t252 * t34 + t257 * t50;
t15 = -pkin(3) * t334 + pkin(8) * t60 + t252 * t57 + t257 * t37;
t14 = -pkin(3) * t74 + pkin(8) * t52 + t252 * t50 + t257 * t34;
t13 = -pkin(4) * t80 + pkin(9) * t18;
t12 = -t24 * t253 + t258 * t26;
t11 = t24 * t258 + t253 * t26;
t10 = -pkin(9) * t47 - t17;
t9 = -pkin(4) * t107 + pkin(9) * t49 + t18;
t8 = t18 * t257 - t321;
t6 = -pkin(8) * t24 + t10 * t257 - t252 * t9;
t5 = -pkin(3) * t107 + pkin(8) * t26 + t10 * t252 + t257 * t9;
t4 = -t253 * t7 + t258 * t8;
t3 = t253 * t8 + t258 * t7;
t2 = -pkin(8) * t7 - pkin(9) * t318 - t13 * t252;
t1 = -pkin(3) * t80 + pkin(8) * t8 - pkin(9) * t321 + t13 * t257;
t29 = [0, 0, 0, 0, 0, qJDD(1), t283, t274, 0, 0, t273 * t254, t224 * t259 + t227 * t254, t303 + t259 * (-t243 + t327), -t272 * t259, t254 * (t245 - t327) + t297, 0, t259 * t215 + pkin(1) * t227 + pkin(6) * (t259 * (-t245 - t327) - t303), -t254 * t215 - pkin(1) * t224 + pkin(6) * (-t297 - t254 * (-t243 - t327)), pkin(1) * (t243 + t245) + (t247 + t248) * t316 + t280, pkin(1) * t215 + pkin(6) * t280, t254 * (t196 * t258 + t222 * t310) - t286, t254 * (t175 * t258 - t177 * t253) - t259 * t205, t254 * (-t211 * t253 + t332) - t259 * t178, t254 * (-t221 * t309 - t253 * t266) + t286, t254 * (t210 * t258 + t304) + t259 * t174, t259 * t220 + t254 * (t221 * t258 - t222 * t253) * t238, t254 * (-pkin(7) * t159 + t305) + t259 * (-pkin(2) * t159 + t146) - pkin(1) * t159 + pkin(6) * (t160 * t259 - t175 * t254), t254 * (-pkin(7) * t162 + t299) + t259 * (-pkin(2) * t162 + t147) - pkin(1) * t162 + pkin(6) * (t163 * t259 - t179 * t254), t254 * t271 + pkin(6) * (t141 * t259 - t189 * t254) + t270 * (-t174 * t253 - t178 * t258), pkin(6) * (t106 * t259 + t182 * t254) - t270 * t271, t254 * (-t122 * t253 + t123 * t258) - t287, t254 * (-t253 * t88 + t258 * t90) - t259 * t169, t254 * (-t136 * t253 + t138 * t258) - t259 * t129, t254 * (-t120 * t253 + t121 * t258) + t287, t254 * (-t137 * t253 + t139 * t258) + t259 * t126, t254 * (-t150 * t253 + t151 * t258) + t259 * t217, t254 * (-pkin(7) * t83 - t253 * t65 + t258 * t94) + t259 * (-pkin(2) * t83 - t267) - pkin(1) * t83 + pkin(6) * (t125 * t254 + t259 * t84), t254 * (-pkin(7) * t95 - t253 * t66 + t258 * t97) + t259 * (-pkin(2) * t95 - t276) - pkin(1) * t95 + pkin(6) * (t254 * t329 + t259 * t96), t254 * (-pkin(7) * t55 - t253 * t32 + t258 * t33) + t259 * (-pkin(2) * t55 - t325) - pkin(1) * t55 + pkin(6) * (t148 * t254 + t259 * t56), t254 * (-pkin(7) * t21 - pkin(8) * t317 - t253 * t38) + t259 * (-pkin(2) * t21 - t326) - pkin(1) * t21 + pkin(6) * (t135 * t254 + t22 * t259), t254 * (-t253 * t42 + t258 * t44) - t288, t254 * (-t23 * t253 + t25 * t258) + t259 * t131, t254 * (-t253 * t61 + t258 * t63) - t259 * t78, t254 * (-t253 * t41 + t258 * t43) + t288, t254 * (-t253 * t62 + t258 * t64) + t259 * t75, t254 * (-t253 * t81 + t258 * t82) + t259 * t213, t254 * (-pkin(7) * t27 - t14 * t253 + t19 * t258) + t259 * (-pkin(2) * t27 - t264) - pkin(1) * t27 + pkin(6) * (t254 * t74 + t259 * t28), t254 * (-pkin(7) * t35 - t15 * t253 + t20 * t258) + t259 * (-pkin(2) * t35 - t262) - pkin(1) * t35 + pkin(6) * (t254 * t334 + t259 * t36), t254 * (-pkin(7) * t11 - t253 * t5 + t258 * t6) + t259 * (-pkin(2) * t11 - t285) - pkin(1) * t11 + pkin(6) * (t107 * t254 + t12 * t259), t254 * (-pkin(7) * t3 - t1 * t253 + t2 * t258) + t259 * (-pkin(2) * t3 - t289) - pkin(1) * t3 + pkin(6) * (t254 * t80 + t259 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, t243 - t245, t241, t237, t290, qJDD(2), -t207, -t208, 0, 0, t196 * t253 - t222 * t309, t175 * t253 + t177 * t258, t211 * t258 + t333, -t221 * t310 + t258 * t266, t210 * t253 - t298, (t221 * t253 + t222 * t258) * t238, pkin(2) * t175 + pkin(7) * t160 - t299, pkin(2) * t179 + pkin(7) * t163 + t305, pkin(2) * t189 + pkin(7) * t141 + t106, -pkin(2) * t182 + pkin(7) * t106, t122 * t258 + t123 * t253, t253 * t90 + t258 * t88, t136 * t258 + t138 * t253, t120 * t258 + t121 * t253, t137 * t258 + t139 * t253, t150 * t258 + t151 * t253, -pkin(2) * t125 + pkin(7) * t84 + t253 * t94 + t258 * t65, -pkin(2) * t329 + pkin(7) * t96 + t253 * t97 + t258 * t66, -pkin(2) * t148 + pkin(7) * t56 + t253 * t33 + t258 * t32, -pkin(2) * t135 + pkin(7) * t22 - pkin(8) * t320 + t258 * t38, t253 * t44 + t258 * t42, t23 * t258 + t25 * t253, t253 * t63 + t258 * t61, t253 * t43 + t258 * t41, t253 * t64 + t258 * t62, t253 * t82 + t258 * t81, -pkin(2) * t74 + pkin(7) * t28 + t14 * t258 + t19 * t253, -pkin(2) * t334 + pkin(7) * t36 + t15 * t258 + t20 * t253, -pkin(2) * t107 + pkin(7) * t12 + t253 * t6 + t258 * t5, -pkin(2) * t80 + pkin(7) * t4 + t1 * t258 + t2 * t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t315, t205, t178, -t315, -t174, -t220, -t146, -t147, 0, 0, t170, t169, t129, -t170, -t126, -t217, t267, t276, t325, t326, t132, -t131, t78, -t132, -t75, -t213, t264, t262, t285, t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t169, t129, -t170, -t126, -t217, -t68, -t69, 0, 0, t132, -t131, t78, -t132, -t75, -t213, t277, t268, t45, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, -t131, t78, -t132, -t75, -t213, -t30, -t31, 0, 0;];
tauJ_reg = t29;
