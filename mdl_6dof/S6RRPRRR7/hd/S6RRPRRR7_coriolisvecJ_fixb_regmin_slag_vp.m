% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:59:46
% EndTime: 2019-03-09 13:59:58
% DurationCPUTime: 4.47s
% Computational Cost: add. (5369->436), mult. (12128->588), div. (0->0), fcn. (8676->8), ass. (0->239)
t192 = sin(qJ(4));
t193 = sin(qJ(2));
t196 = cos(qJ(4));
t197 = cos(qJ(2));
t135 = t192 * t193 + t196 * t197;
t325 = t135 * qJD(1);
t342 = qJD(5) + t325;
t272 = qJD(1) * t197;
t273 = qJD(1) * t193;
t129 = -t192 * t272 + t196 * t273;
t184 = qJD(2) - qJD(4);
t191 = sin(qJ(5));
t195 = cos(qJ(5));
t94 = t129 * t191 + t195 * t184;
t352 = t342 * t94;
t266 = qJD(4) * t196;
t267 = qJD(4) * t192;
t269 = qJD(2) * t197;
t351 = t192 * t269 + t193 * t266 - t197 * t267;
t175 = pkin(7) * t273;
t142 = pkin(8) * t273 - t175;
t350 = qJD(3) - t142;
t198 = -pkin(2) - pkin(3);
t251 = t198 * qJD(2);
t113 = t251 + t350;
t271 = qJD(2) * t193;
t321 = pkin(7) - pkin(8);
t143 = t321 * t271;
t185 = qJD(2) * qJD(3);
t119 = -qJD(1) * t143 + t185;
t176 = pkin(7) * t272;
t144 = -pkin(8) * t272 + t176;
t186 = qJD(2) * qJ(3);
t130 = t144 + t186;
t259 = qJD(1) * qJD(2);
t244 = t197 * t259;
t167 = pkin(7) * t244;
t132 = -pkin(8) * t244 + t167;
t34 = t113 * t267 + t192 * t119 + t130 * t266 - t132 * t196;
t210 = t135 * qJD(4);
t245 = t193 * t259;
t278 = t192 * t245 + t196 * t244;
t79 = -qJD(1) * t210 + t278;
t304 = t191 * t79;
t96 = t129 * t195 - t184 * t191;
t40 = qJD(5) * t96 + t304;
t11 = pkin(5) * t40 + t34;
t194 = cos(qJ(6));
t190 = sin(qJ(6));
t283 = t190 * t195;
t137 = t191 * t194 + t283;
t284 = t190 * t191;
t134 = -t194 * t195 + t284;
t264 = qJD(5) * t195;
t327 = -qJD(6) * t195 - t264;
t208 = t327 * t194;
t258 = qJD(5) + qJD(6);
t313 = t134 * t325 + t258 * t284 + t208;
t72 = t113 * t196 - t192 * t130;
t59 = pkin(4) * t184 - t72;
t35 = pkin(5) * t94 + t59;
t131 = -qJD(1) * pkin(1) - pkin(2) * t272 - qJ(3) * t273;
t107 = pkin(3) * t272 - t131;
t54 = pkin(4) * t325 - pkin(9) * t129 + t107;
t73 = t192 * t113 + t196 * t130;
t60 = -pkin(9) * t184 + t73;
t21 = -t191 * t60 + t195 * t54;
t15 = -pkin(10) * t96 + t21;
t10 = pkin(5) * t342 + t15;
t22 = t191 * t54 + t195 * t60;
t16 = -pkin(10) * t94 + t22;
t305 = t16 * t194;
t6 = t10 * t190 + t305;
t349 = -t11 * t137 - t6 * t129 + t313 * t35;
t311 = (t258 + t325) * t137;
t5 = t10 * t194 - t16 * t190;
t348 = -t11 * t134 + t5 * t129 - t311 * t35;
t118 = qJD(6) + t342;
t219 = t190 * t94 - t194 * t96;
t80 = t351 * qJD(1) - t196 * t245;
t306 = t137 * t80;
t347 = t118 * t313 - t129 * t219 - t306;
t307 = t134 * t80;
t44 = t190 * t96 + t194 * t94;
t346 = t118 * t311 - t129 * t44 + t307;
t262 = qJD(6) * t194;
t263 = qJD(6) * t190;
t265 = qJD(5) * t191;
t39 = -t129 * t265 - t184 * t264 + t195 * t79;
t8 = -t190 * t40 + t194 * t39 - t94 * t262 - t263 * t96;
t345 = t8 * t137 + t219 * t313;
t202 = qJD(6) * t219 - t190 * t39 - t194 * t40;
t344 = -t134 * t8 + t137 * t202 + t219 * t311 + t313 * t44;
t343 = t219 * t44;
t326 = -t192 * qJ(3) + t196 * t198;
t114 = t196 * qJD(3) + t326 * qJD(4);
t84 = t142 * t196 + t144 * t192;
t296 = t114 - t84;
t341 = t219 ^ 2 - t44 ^ 2;
t340 = t118 * t44 + t8;
t14 = t16 * t263;
t339 = t35 * t44 + t14;
t234 = t193 * t251;
t179 = t193 * qJD(3);
t276 = qJ(3) * t244 + qJD(1) * t179;
t93 = qJD(1) * t234 + t276;
t27 = pkin(4) * t80 - pkin(9) * t79 + t93;
t25 = t195 * t27;
t33 = t113 * t266 + t196 * t119 - t130 * t267 + t192 * t132;
t203 = -qJD(5) * t22 - t191 * t33 + t25;
t2 = pkin(5) * t80 - pkin(10) * t39 + t203;
t214 = t191 * t27 + t195 * t33 + t54 * t264 - t265 * t60;
t3 = -pkin(10) * t40 + t214;
t253 = -t190 * t3 + t194 * t2;
t338 = -qJD(6) * t6 + t35 * t219 + t253;
t297 = t39 * t191;
t330 = t195 * t342;
t337 = t330 * t96 + t297;
t336 = -t118 * t219 + t202;
t303 = t191 * t80;
t335 = -t129 * t96 + t330 * t342 + t303;
t334 = (-t39 + t352) * t195 + (t342 * t96 + t40) * t191;
t333 = -0.2e1 * t259;
t332 = t342 * t59;
t218 = qJ(3) * t196 + t192 * t198;
t295 = qJD(4) * t218 + t144 * t196 + t350 * t192;
t331 = t129 * t184 + t80;
t136 = -t192 * t197 + t193 * t196;
t77 = t137 * t136;
t153 = t321 * t193;
t155 = t321 * t197;
t97 = -t153 * t196 + t192 * t155;
t329 = -t190 * t265 - t191 * t263;
t290 = t325 * t191;
t328 = (t265 + t290) * pkin(5);
t301 = t195 * t80;
t322 = t191 * t342 ^ 2 - t129 * t94 - t301;
t320 = pkin(9) + pkin(10);
t319 = pkin(5) * t195;
t141 = -pkin(9) + t218;
t317 = pkin(10) - t141;
t81 = pkin(4) * t129 + pkin(9) * t325;
t316 = t191 * t81 + t195 * t72;
t171 = qJ(3) * t272;
t121 = t198 * t273 + t171;
t55 = t121 - t81;
t315 = t191 * t55 + t195 * t84;
t150 = -t197 * pkin(2) - t193 * qJ(3) - pkin(1);
t133 = t197 * pkin(3) - t150;
t70 = pkin(4) * t135 - pkin(9) * t136 + t133;
t98 = t153 * t192 + t155 * t196;
t92 = t195 * t98;
t310 = t191 * t70 + t92;
t309 = -t328 + t295;
t308 = qJD(2) * pkin(2);
t89 = qJD(2) * t135 - t210;
t302 = t191 * t89;
t300 = t195 * t89;
t299 = t195 * t96;
t298 = t34 * t195;
t294 = t118 * t129;
t293 = t342 * t129;
t291 = t325 * t184;
t289 = t129 * t325;
t287 = t136 * t191;
t286 = t136 * t195;
t200 = qJD(1) ^ 2;
t281 = t197 * t200;
t199 = qJD(2) ^ 2;
t280 = t199 * t193;
t279 = t199 * t197;
t275 = qJ(3) * t269 + t179;
t187 = t193 ^ 2;
t274 = -t197 ^ 2 + t187;
t270 = qJD(2) * t196;
t268 = qJD(4) * t118;
t256 = pkin(10) * t290;
t254 = t193 * t281;
t252 = qJD(5) * t320;
t250 = t136 * t264;
t246 = -t129 ^ 2 + t325 ^ 2;
t243 = qJD(6) * t10 + t3;
t241 = qJD(5) * t317;
t240 = pkin(1) * t333;
t239 = qJD(3) - t308;
t236 = qJD(1) * t150 + t131;
t235 = t184 ^ 2;
t233 = -t73 + t328;
t140 = pkin(4) - t326;
t108 = t317 * t191;
t232 = -qJD(6) * t108 - t195 * t114 - t191 * t241 - t256 + t315;
t109 = t317 * t195;
t217 = -pkin(10) * t195 * t325 - pkin(5) * t129;
t53 = t195 * t55;
t231 = -qJD(6) * t109 + t296 * t191 - t195 * t241 + t217 + t53;
t152 = t320 * t191;
t229 = qJD(6) * t152 + t191 * t252 + t256 + t316;
t154 = t320 * t195;
t76 = t195 * t81;
t228 = qJD(6) * t154 - t191 * t72 + t195 * t252 - t217 + t76;
t223 = t22 * t129 + t34 * t191;
t221 = -t21 * t129 - t298;
t106 = pkin(2) * t245 - t276;
t122 = pkin(2) * t271 - t275;
t216 = -pkin(7) * t199 - qJD(1) * t122 - t106;
t215 = t250 + t302;
t102 = t234 + t275;
t88 = -t193 * t270 + t351;
t32 = pkin(4) * t88 - pkin(9) * t89 + t102;
t145 = qJD(2) * t155;
t49 = -t97 * qJD(4) - t143 * t196 + t145 * t192;
t213 = t191 * t32 + t195 * t49 + t70 * t264 - t265 * t98;
t212 = -pkin(9) * t80 + t332;
t211 = -t141 * t80 - t332;
t206 = t107 * t129 + t34;
t205 = -t107 * t325 + t33;
t50 = qJD(4) * t98 - t143 * t192 - t145 * t196;
t148 = -pkin(7) * t245 + t185;
t149 = t175 + t239;
t151 = t176 + t186;
t201 = t148 * t197 + (t149 * t197 + (-t151 + t176) * t193) * qJD(2);
t174 = -pkin(4) - t319;
t139 = pkin(2) * t273 - t171;
t128 = t191 * t273 + t195 * t270;
t125 = -t191 * t270 + t195 * t273;
t124 = t140 + t319;
t78 = t134 * t136;
t71 = pkin(5) * t287 + t97;
t62 = t195 * t70;
t57 = t80 * t135;
t30 = t195 * t32;
t28 = -pkin(10) * t287 + t310;
t23 = pkin(5) * t135 - pkin(10) * t286 - t191 * t98 + t62;
t19 = pkin(5) * t215 + t50;
t13 = t89 * t283 + (t258 * t286 + t302) * t194 + t329 * t136;
t12 = -t134 * t89 - t258 * t77;
t7 = -pkin(10) * t215 + t213;
t4 = -pkin(10) * t300 + pkin(5) * t88 - t191 * t49 + t30 + (-t92 + (pkin(10) * t136 - t70) * t191) * qJD(5);
t1 = [0, 0, 0, 0.2e1 * t193 * t244, t274 * t333, t279, -t280, 0, -pkin(7) * t279 + t193 * t240, pkin(7) * t280 + t197 * t240, t197 * t216 + t236 * t271, t201, t193 * t216 - t236 * t269, pkin(7) * t201 + t106 * t150 + t122 * t131, t129 * t89 + t136 * t79, -t129 * t88 - t135 * t79 - t136 * t80 - t325 * t89, -t89 * t184, t88 * t184, 0, t102 * t325 + t107 * t88 + t133 * t80 + t135 * t93 + t184 * t50, t102 * t129 + t107 * t89 + t133 * t79 + t136 * t93 + t184 * t49, t89 * t299 + (t195 * t39 - t265 * t96) * t136 (-t191 * t96 - t195 * t94) * t89 + (-t297 - t195 * t40 + (t191 * t94 - t299) * qJD(5)) * t136, t80 * t286 + t135 * t39 + t88 * t96 + (-t136 * t265 + t300) * t342, -t135 * t40 - t215 * t342 - t287 * t80 - t88 * t94, t342 * t88 + t57 (-t264 * t98 + t30) * t342 + t62 * t80 + (-t264 * t60 + t25) * t135 + t21 * t88 + t50 * t94 + t97 * t40 + t59 * t250 + ((-qJD(5) * t70 - t49) * t342 - t98 * t80 + (-qJD(5) * t54 - t33) * t135 + t34 * t136 + t59 * t89) * t191, -t213 * t342 - t310 * t80 - t214 * t135 - t22 * t88 + t50 * t96 + t97 * t39 + t59 * t300 + (-t265 * t59 + t298) * t136, -t12 * t219 - t78 * t8, -t12 * t44 + t13 * t219 - t202 * t78 - t77 * t8, t118 * t12 + t135 * t8 - t219 * t88 - t78 * t80, -t118 * t13 + t135 * t202 - t44 * t88 - t77 * t80, t118 * t88 + t57 (-t190 * t7 + t194 * t4) * t118 + (-t190 * t28 + t194 * t23) * t80 + t253 * t135 + t5 * t88 + t19 * t44 - t71 * t202 + t11 * t77 + t35 * t13 + ((-t190 * t23 - t194 * t28) * t118 - t6 * t135) * qJD(6), -t11 * t78 + t35 * t12 + t14 * t135 - t19 * t219 - t6 * t88 + t71 * t8 + (-(-qJD(6) * t28 + t4) * t118 - t23 * t80 - t2 * t135) * t190 + (-(qJD(6) * t23 + t7) * t118 - t28 * t80 - t243 * t135) * t194; 0, 0, 0, -t254, t274 * t200, 0, 0, 0, t200 * pkin(1) * t193, pkin(1) * t281 (-t131 * t193 + t139 * t197) * qJD(1) ((t151 - t186) * t193 + (-t149 + t239) * t197) * qJD(1), 0.2e1 * t185 + (t131 * t197 + t139 * t193) * qJD(1), qJ(3) * t148 + qJD(3) * t151 - t131 * t139 + (t151 * t193 + (-t149 - t308) * t197) * qJD(1) * pkin(7), -t289, t246, qJD(4) * t325 - t278 + t291, t331, 0, -t121 * t325 + t184 * t295 + t206, -t121 * t129 + t184 * t296 + t205, -t337, t334, -t335, t322, t293, t140 * t40 + t295 * t94 + (-t141 * t264 - t53) * t342 + (-t296 * t342 + t211) * t191 - t221, t140 * t39 + t295 * t96 + (t141 * t265 + t315) * t342 + (-t114 * t342 + t211) * t195 - t223, -t345, -t344, t347, t346, t294 (t108 * t194 + t109 * t190) * t80 - t124 * t202 + t309 * t44 + (t190 * t232 - t194 * t231) * t118 + t348 -(t108 * t190 - t109 * t194) * t80 + t124 * t8 - t309 * t219 + (t190 * t231 + t194 * t232) * t118 + t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, 0, -t187 * t200 - t199, -qJD(2) * t151 + t131 * t273 + t167, 0, 0, 0, 0, 0, -t192 * t235 - t273 * t325, -t129 * t273 - t196 * t235, 0, 0, 0, 0, 0, -t196 * t40 + (-t191 * t266 - t125) * t342 + (-t184 * t94 - t264 * t342 - t303) * t192, -t196 * t39 + (-t195 * t266 + t128) * t342 + (-t184 * t96 + t265 * t342 - t301) * t192, 0, 0, 0, 0, 0 -(t125 * t194 - t128 * t190) * t118 + (-t137 * t268 + t202) * t196 + ((t208 - t329) * t118 - t306 - t184 * t44) * t192 (t125 * t190 + t128 * t194) * t118 + (t134 * t268 - t8) * t196 + (-(t190 * t327 - t191 * t262 - t194 * t265) * t118 + t307 + t184 * t219) * t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, -t246, t79 - t291, -t331, 0, -t184 * t73 - t206, -t184 * t72 - t205, t337, -t334, t335, -t322, -t293, -pkin(4) * t40 - t73 * t94 + (-pkin(9) * t264 - t76) * t342 + (t342 * t72 + t212) * t191 + t221, -pkin(4) * t39 - t73 * t96 + (pkin(9) * t265 + t316) * t342 + t212 * t195 + t223, t345, t344, -t347, -t346, -t294 (-t152 * t194 - t154 * t190) * t80 - t174 * t202 + t233 * t44 + (t190 * t229 - t194 * t228) * t118 - t348 -(-t152 * t190 + t154 * t194) * t80 + t174 * t8 - t233 * t219 + (t190 * t228 + t194 * t229) * t118 - t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t94, -t94 ^ 2 + t96 ^ 2, t39 + t352, -t304 + (-qJD(5) + t342) * t96, t80, t22 * t342 - t59 * t96 + t203, t21 * t342 + t59 * t94 - t214, -t343, t341, t340, t336, t80 -(-t15 * t190 - t305) * t118 + (-t118 * t263 + t194 * t80 - t44 * t96) * pkin(5) + t338 (-t118 * t16 - t2) * t190 + (t118 * t15 - t243) * t194 + (-t118 * t262 - t190 * t80 + t219 * t96) * pkin(5) + t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t343, t341, t340, t336, t80, t118 * t6 + t338, t118 * t5 - t190 * t2 - t194 * t243 + t339;];
tauc_reg  = t1;
