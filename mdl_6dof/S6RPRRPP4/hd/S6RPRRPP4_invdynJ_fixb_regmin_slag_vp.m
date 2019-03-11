% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:40:50
% EndTime: 2019-03-09 04:41:01
% DurationCPUTime: 4.81s
% Computational Cost: add. (9261->485), mult. (21909->603), div. (0->0), fcn. (16718->14), ass. (0->232)
t214 = cos(pkin(9));
t221 = cos(qJ(3));
t289 = t221 * t214;
t186 = qJD(1) * t289;
t212 = sin(pkin(9));
t218 = sin(qJ(3));
t299 = t212 * t218;
t269 = qJD(1) * t299;
t152 = t186 - t269;
t141 = qJD(4) - t152;
t209 = pkin(9) + qJ(3);
t199 = sin(t209);
t201 = cos(t209);
t219 = sin(qJ(1));
t222 = cos(qJ(1));
t258 = g(1) * t222 + g(2) * t219;
t232 = -g(3) * t201 + t258 * t199;
t280 = qJD(1) * qJD(2);
t329 = pkin(7) + qJ(2);
t338 = t329 * qJDD(1) + t280;
t135 = t338 * t212;
t136 = t338 * t214;
t172 = t329 * t212;
t166 = qJD(1) * t172;
t173 = t329 * t214;
t167 = qJD(1) * t173;
t284 = qJD(3) * t221;
t285 = qJD(3) * t218;
t244 = -t221 * t135 - t218 * t136 + t166 * t285 - t167 * t284;
t67 = -qJDD(3) * pkin(3) - t244;
t356 = -qJD(4) * pkin(8) * t141 + t232 - t67;
t164 = t212 * t221 + t214 * t218;
t211 = sin(pkin(10));
t213 = cos(pkin(10));
t217 = sin(qJ(4));
t220 = cos(qJ(4));
t343 = -t211 * t217 + t213 * t220;
t92 = t343 * t164;
t283 = qJD(4) * t217;
t310 = t152 * t217;
t355 = t283 - t310;
t281 = t220 * qJD(3);
t342 = t164 * qJD(1);
t122 = t217 * t342 - t281;
t124 = qJD(3) * t217 + t220 * t342;
t75 = t213 * t122 + t124 * t211;
t354 = t141 * t75;
t163 = t211 * t220 + t213 * t217;
t151 = t163 * qJD(4);
t318 = -t163 * t152 + t151;
t282 = qJD(4) * t220;
t317 = -t343 * t152 - t211 * t283 + t213 * t282;
t197 = pkin(4) * t220 + pkin(3);
t169 = t201 * t197;
t195 = pkin(2) * t214 + pkin(1);
t353 = -t195 - t169;
t315 = qJDD(1) * pkin(1);
t345 = g(1) * t219 - g(2) * t222;
t246 = -qJDD(2) + t315 + t345;
t162 = -t289 + t299;
t155 = t162 * qJD(3);
t268 = t164 * t282;
t352 = -t155 * t217 + t268;
t351 = qJD(3) * t342;
t251 = -t122 * t211 + t213 * t124;
t350 = t251 ^ 2;
t328 = qJ(5) + pkin(8);
t265 = qJD(4) * t328;
t147 = t220 * qJD(5) - t217 * t265;
t235 = -t217 * qJD(5) - t220 * t265;
t347 = -t166 * t221 - t218 * t167;
t107 = pkin(3) * t342 - pkin(8) * t152;
t96 = t220 * t107;
t51 = -qJ(5) * t152 * t220 + pkin(4) * t342 - t217 * t347 + t96;
t320 = t217 * t107 + t220 * t347;
t56 = -qJ(5) * t310 + t320;
t326 = (t235 - t51) * t213 + (-t147 + t56) * t211;
t118 = t172 * t221 + t218 * t173;
t114 = -t218 * t166 + t221 * t167;
t344 = t355 * pkin(4) - t114;
t341 = qJ(2) * qJDD(1);
t278 = t214 * qJDD(1);
t279 = t212 * qJDD(1);
t270 = qJD(3) * t186 + t218 * t278 + t221 * t279;
t111 = -qJD(3) * t269 + t270;
t71 = qJD(4) * t281 + t217 * qJDD(3) + t220 * t111 - t283 * t342;
t72 = qJD(4) * t124 - t220 * qJDD(3) + t217 * t111;
t39 = t211 * t71 + t213 * t72;
t40 = -t211 * t72 + t213 * t71;
t340 = pkin(5) * t39 - qJ(6) * t40 - qJD(6) * t251;
t210 = qJ(4) + pkin(10);
t202 = cos(t210);
t301 = t202 * t222;
t200 = sin(t210);
t302 = t200 * t219;
t128 = t201 * t302 + t301;
t288 = t222 * t200;
t292 = t219 * t202;
t130 = t201 * t288 - t292;
t156 = t164 * qJD(3);
t252 = t218 * t279 - t221 * t278;
t112 = qJD(1) * t156 + t252;
t104 = qJDD(4) + t112;
t106 = qJD(3) * pkin(8) + t114;
t170 = -qJD(1) * t195 + qJD(2);
t83 = -t152 * pkin(3) - pkin(8) * t342 + t170;
t59 = t106 * t220 + t217 * t83;
t168 = -qJDD(1) * t195 + qJDD(2);
t68 = t112 * pkin(3) - t111 * pkin(8) + t168;
t65 = t220 * t68;
t250 = -t218 * t135 + t221 * t136;
t66 = qJDD(3) * pkin(8) + qJD(3) * t347 + t250;
t12 = t104 * pkin(4) - t71 * qJ(5) - qJD(4) * t59 - t124 * qJD(5) - t217 * t66 + t65;
t238 = -t106 * t283 + t217 * t68 + t220 * t66 + t83 * t282;
t15 = -qJ(5) * t72 - qJD(5) * t122 + t238;
t3 = t213 * t12 - t211 * t15;
t272 = -qJDD(6) + t3;
t331 = g(3) * t199;
t105 = -qJD(3) * pkin(3) - t347;
t73 = pkin(4) * t122 + qJD(5) + t105;
t35 = pkin(5) * t75 - qJ(6) * t251 + t73;
t339 = g(1) * t130 + g(2) * t128 + t200 * t331 - t35 * t251 + t272;
t233 = t201 * t258 + t331;
t337 = t141 ^ 2;
t335 = pkin(5) * t104;
t50 = -qJ(5) * t122 + t59;
t46 = t213 * t50;
t58 = -t106 * t217 + t220 * t83;
t49 = -qJ(5) * t124 + t58;
t21 = t211 * t49 + t46;
t330 = t21 * t251;
t4 = t211 * t12 + t213 * t15;
t110 = pkin(3) * t162 - pkin(8) * t164 - t195;
t119 = -t172 * t218 + t173 * t221;
t115 = t220 * t119;
t247 = qJ(5) * t155 - qJD(5) * t164;
t84 = -t162 * qJD(2) - t118 * qJD(3);
t108 = pkin(3) * t156 + pkin(8) * t155;
t97 = t220 * t108;
t26 = t156 * pkin(4) - t217 * t84 + t97 + t247 * t220 + (-t115 + (qJ(5) * t164 - t110) * t217) * qJD(4);
t277 = t217 * t108 + t110 * t282 + t220 * t84;
t30 = -qJ(5) * t268 + (-qJD(4) * t119 + t247) * t217 + t277;
t9 = t211 * t26 + t213 * t30;
t44 = pkin(4) * t141 + t49;
t20 = t211 * t44 + t46;
t29 = t211 * t51 + t213 * t56;
t307 = t164 * t220;
t99 = t220 * t110;
t53 = pkin(4) * t162 - qJ(5) * t307 - t119 * t217 + t99;
t308 = t164 * t217;
t319 = t217 * t110 + t115;
t60 = -qJ(5) * t308 + t319;
t34 = t211 * t53 + t213 * t60;
t327 = t318 * pkin(5) - t317 * qJ(6) - qJD(6) * t163 + t344;
t325 = pkin(5) * t342 - t326;
t24 = qJ(6) * t342 + t29;
t89 = t213 * t147 + t211 * t235;
t324 = t89 - t24;
t322 = t211 * t50;
t321 = t71 * t217;
t314 = t122 * t141;
t313 = t122 * t342;
t312 = t124 * t141;
t311 = t124 * t342;
t304 = t199 * t219;
t303 = t199 * t222;
t295 = t217 * t104;
t294 = t217 * t219;
t293 = t217 * t222;
t291 = t219 * t220;
t93 = t220 * t104;
t290 = t220 * t222;
t22 = t213 * t49 - t322;
t287 = qJD(6) - t22;
t286 = t212 ^ 2 + t214 ^ 2;
t276 = t104 * qJ(6) + t4;
t273 = t201 * t293;
t267 = t328 * t217;
t264 = -qJD(4) * t83 - t66;
t263 = t286 * qJD(1) ^ 2;
t261 = t141 * t220;
t85 = t164 * qJD(2) - t172 * t285 + t173 * t284;
t260 = 0.2e1 * t286;
t259 = -g(1) * t304 + g(2) * t303;
t256 = -t106 * t282 + t65;
t255 = -t75 ^ 2 - t350;
t254 = pkin(4) * t308 + t118;
t253 = pkin(5) * t202 + qJ(6) * t200;
t8 = -t211 * t30 + t213 * t26;
t19 = t213 * t44 - t322;
t33 = -t211 * t60 + t213 * t53;
t245 = -t355 * t141 + t93;
t242 = t352 * pkin(4) + t85;
t142 = t201 * t294 + t290;
t239 = -t155 * t220 - t164 * t283;
t237 = -pkin(8) * t104 + t105 * t141;
t234 = t246 + t315;
t231 = pkin(4) * t294 + t329 * t219 - t353 * t222 + t328 * t303;
t230 = pkin(4) * t293 + t353 * t219 + t222 * t329 - t328 * t304;
t229 = -t163 * t39 + t251 * t318 - t317 * t75 - t343 * t40;
t227 = t260 * t280 - t258;
t38 = t72 * pkin(4) + qJDD(5) + t67;
t174 = t328 * t220;
t120 = t174 * t211 + t213 * t267;
t121 = t213 * t174 - t211 * t267;
t226 = t120 * t40 - t121 * t39 - t89 * t75 - t233;
t225 = t38 - t232;
t194 = -pkin(4) * t213 - pkin(5);
t190 = pkin(4) * t211 + qJ(6);
t188 = pkin(4) * t291;
t145 = t201 * t290 + t294;
t144 = -t273 + t291;
t143 = -t201 * t291 + t293;
t131 = t201 * t301 + t302;
t129 = t201 * t292 - t288;
t109 = -pkin(5) * t343 - qJ(6) * t163 - t197;
t91 = t163 * t164;
t63 = t151 * t164 + t343 * t155;
t62 = -qJD(4) * t92 + t155 * t163;
t45 = pkin(5) * t91 - qJ(6) * t92 + t254;
t42 = pkin(4) * t124 + pkin(5) * t251 + qJ(6) * t75;
t32 = -pkin(5) * t162 - t33;
t31 = qJ(6) * t162 + t34;
t18 = qJ(6) * t141 + t20;
t17 = -pkin(5) * t141 + qJD(6) - t19;
t16 = -pkin(5) * t62 + qJ(6) * t63 - qJD(6) * t92 + t242;
t7 = -pkin(5) * t156 - t8;
t6 = qJ(6) * t156 + qJD(6) * t162 + t9;
t5 = t38 + t340;
t2 = -t272 - t335;
t1 = qJD(6) * t141 + t276;
t10 = [qJDD(1), t345, t258, t234 * t214, -t234 * t212, t260 * t341 + t227, pkin(1) * t246 + (t286 * t341 + t227) * qJ(2), t111 * t164 - t155 * t342, -t111 * t162 - t112 * t164 - t152 * t155 - t156 * t342, -qJD(3) * t155 + qJDD(3) * t164, -qJD(3) * t156 - qJDD(3) * t162, 0, -t85 * qJD(3) - t118 * qJDD(3) - t195 * t112 + t170 * t156 + t168 * t162 + t201 * t345, -qJD(3) * t84 - qJDD(3) * t119 - t111 * t195 - t155 * t170 + t164 * t168 + t259, t124 * t239 + t307 * t71 -(-t122 * t220 - t124 * t217) * t155 + (-t321 - t220 * t72 + (t122 * t217 - t124 * t220) * qJD(4)) * t164, t124 * t156 + t141 * t239 + t71 * t162 + t164 * t93, -t122 * t156 - t352 * t141 - t72 * t162 - t164 * t295, t104 * t162 + t141 * t156 (-t119 * t282 + t97) * t141 + t99 * t104 + t256 * t162 + t58 * t156 + t85 * t122 + t118 * t72 + t105 * t268 - g(1) * t143 - g(2) * t145 + ((-qJD(4) * t110 - t84) * t141 - t119 * t104 + t264 * t162 + t67 * t164 - t105 * t155) * t217 -(-t119 * t283 + t277) * t141 - t319 * t104 - t238 * t162 - t59 * t156 + t85 * t124 + t118 * t71 + t67 * t307 - g(1) * t142 - g(2) * t144 + t239 * t105, t19 * t63 + t20 * t62 - t251 * t8 - t3 * t92 - t33 * t40 - t34 * t39 - t4 * t91 - t75 * t9 - t259, -g(1) * t230 - g(2) * t231 + t19 * t8 + t20 * t9 + t242 * t73 + t254 * t38 + t3 * t33 + t4 * t34, g(1) * t129 - g(2) * t131 - t104 * t32 - t141 * t7 - t156 * t17 + t16 * t75 - t162 * t2 - t35 * t62 + t39 * t45 + t5 * t91, -t1 * t91 - t17 * t63 + t18 * t62 + t2 * t92 + t251 * t7 - t31 * t39 + t32 * t40 - t6 * t75 - t259, g(1) * t128 - g(2) * t130 + t1 * t162 + t104 * t31 + t141 * t6 + t156 * t18 - t16 * t251 + t35 * t63 - t40 * t45 - t5 * t92, t1 * t31 + t18 * t6 + t5 * t45 + t35 * t16 + t2 * t32 + t17 * t7 - g(1) * (-t129 * pkin(5) - t128 * qJ(6) + t230) - g(2) * (pkin(5) * t131 + qJ(6) * t130 + t231); 0, 0, 0, -t278, t279, -t263, -qJ(2) * t263 - t246, 0, 0, 0, 0, 0, t252 + 0.2e1 * t351 (t152 - t269) * qJD(3) + t270, 0, 0, 0, 0, 0, t245 - t313, -t337 * t220 - t295 - t311, t229, t4 * t163 - t19 * t318 + t20 * t317 + t3 * t343 - t342 * t73 - t345, t104 * t343 - t141 * t318 - t342 * t75, t229, t163 * t104 + t141 * t317 + t251 * t342, t1 * t163 + t17 * t318 + t18 * t317 - t2 * t343 - t342 * t35 - t345; 0, 0, 0, 0, 0, 0, 0, -t342 * t152, -t152 ^ 2 + t342 ^ 2 (-t152 - t269) * qJD(3) + t270, -t252, qJDD(3), t114 * qJD(3) - t170 * t342 + t232 + t244, -t170 * t152 + t233 - t250, t124 * t261 + t321 (t71 - t314) * t220 + (-t72 - t312) * t217, t141 * t261 + t295 - t311, t245 + t313, -t141 * t342, -pkin(3) * t72 - t114 * t122 - t96 * t141 - t58 * t342 + (t141 * t347 + t237) * t217 + t356 * t220, -pkin(3) * t71 - t114 * t124 + t320 * t141 - t217 * t356 + t237 * t220 + t59 * t342, -t163 * t3 - t19 * t317 - t20 * t318 - t251 * t326 + t29 * t75 + t343 * t4 + t226, t4 * t121 - t3 * t120 - t38 * t197 - g(3) * (t199 * t328 + t169) + t344 * t73 + (t89 - t29) * t20 + t326 * t19 + t258 * (t197 * t199 - t201 * t328) -t120 * t104 + t109 * t39 - t141 * t325 + t17 * t342 + t202 * t232 + t318 * t35 + t327 * t75 - t343 * t5, t1 * t343 + t163 * t2 + t17 * t317 - t18 * t318 + t24 * t75 + t251 * t325 + t226, t104 * t121 - t109 * t40 + t141 * t324 - t163 * t5 - t18 * t342 + t200 * t232 - t251 * t327 - t317 * t35, -g(3) * t169 + t1 * t121 + t5 * t109 + t2 * t120 + t327 * t35 + t324 * t18 + t325 * t17 + (-g(3) * t253 - t258 * t328) * t201 + (-g(3) * t328 + t258 * (t197 + t253)) * t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124 * t122, -t122 ^ 2 + t124 ^ 2, t71 + t314, t312 - t72, t104, -g(1) * t144 + g(2) * t142 - t105 * t124 + t59 * t141 + (t264 + t331) * t217 + t256, g(1) * t145 - g(2) * t143 + t105 * t122 + t141 * t58 + t220 * t331 - t238, t20 * t251 - t330 + (-t211 * t39 - t213 * t40) * pkin(4) + (-t19 + t22) * t75, -g(1) * t188 + t19 * t21 - t20 * t22 + (g(2) * t290 - t73 * t124 + t4 * t211 + t3 * t213 + t217 * t233) * pkin(4), t21 * t141 - t42 * t75 + (pkin(5) - t194) * t104 + t339, t18 * t251 - t190 * t39 + t194 * t40 - t330 + (t17 - t287) * t75, -t202 * t331 - g(1) * t131 - g(2) * t129 + t190 * t104 - t35 * t75 + t42 * t251 + (0.2e1 * qJD(6) - t22) * t141 + t276, t1 * t190 + t2 * t194 - t35 * t42 - t17 * t21 - g(1) * (-pkin(4) * t273 - pkin(5) * t130 + qJ(6) * t131 + t188) - g(2) * (-pkin(4) * t142 - t128 * pkin(5) + t129 * qJ(6)) + t287 * t18 - (-pkin(4) * t217 - pkin(5) * t200 + qJ(6) * t202) * t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t255, t19 * t251 + t20 * t75 + t225, t141 * t251 + t39, t255, -t40 + t354, -t17 * t251 + t18 * t75 + t225 + t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251 * t75 - qJDD(4) - t252 - t351, t40 + t354, -t337 - t350, -t141 * t18 - t335 - t339;];
tau_reg  = t10;
