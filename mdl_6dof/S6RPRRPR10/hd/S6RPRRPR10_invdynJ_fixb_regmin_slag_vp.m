% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% tau_reg [6x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:13
% EndTime: 2019-03-09 05:37:24
% DurationCPUTime: 5.08s
% Computational Cost: add. (3783->492), mult. (7307->634), div. (0->0), fcn. (4759->8), ass. (0->241)
t162 = sin(qJ(3));
t307 = g(3) * t162;
t166 = cos(qJ(3));
t167 = cos(qJ(1));
t163 = sin(qJ(1));
t308 = g(1) * t163;
t328 = g(2) * t167 - t308;
t323 = t328 * t166;
t336 = t307 + t323;
t268 = qJD(1) * t162;
t140 = qJD(4) + t268;
t161 = sin(qJ(4));
t316 = pkin(4) + pkin(5);
t240 = t316 * t161;
t165 = cos(qJ(4));
t290 = qJ(5) * t165;
t192 = -t290 + t240;
t326 = qJD(4) - qJD(6);
t253 = t165 * qJD(3);
t267 = qJD(1) * t166;
t107 = t161 * t267 - t253;
t265 = qJD(3) * t161;
t109 = t165 * t267 + t265;
t160 = sin(qJ(6));
t164 = cos(qJ(6));
t201 = -t164 * t107 + t109 * t160;
t56 = t107 * t160 + t109 * t164;
t335 = -t201 ^ 2 + t56 ^ 2;
t260 = qJD(4) * t162;
t223 = qJD(1) + t260;
t266 = qJD(3) * t109;
t251 = qJD(1) * qJD(3);
t232 = t166 * t251;
t248 = t162 * qJDD(1);
t105 = qJDD(4) + t232 + t248;
t289 = t105 * t165;
t258 = qJD(4) * t166;
t234 = t161 * t258;
t183 = t162 * t253 + t234;
t247 = t166 * qJDD(1);
t44 = qJD(1) * t183 - qJD(4) * t253 - t161 * qJDD(3) - t165 * t247;
t334 = (t161 * t223 - t166 * t253) * t140 + (t266 - t289) * t162 + t166 * t44;
t276 = t167 * t165;
t279 = t163 * t161;
t90 = t162 * t279 - t276;
t278 = t163 * t165;
t91 = t161 * t167 + t162 * t278;
t207 = t160 * t91 - t164 * t90;
t281 = t162 * t167;
t92 = t161 * t281 + t278;
t93 = t162 * t276 - t279;
t229 = t160 * t93 - t92 * t164;
t169 = -pkin(1) - pkin(7);
t135 = t169 * qJD(1) + qJD(2);
t97 = -qJD(3) * pkin(3) - t135 * t166;
t182 = qJ(5) * t109 - t97;
t23 = -t316 * t107 + t182;
t259 = qJD(4) * t165;
t261 = qJD(4) * t161;
t220 = pkin(3) * t166 + pkin(8) * t162;
t106 = qJD(3) * t220 + qJD(2);
t309 = pkin(8) * t166;
t219 = pkin(3) * t162 - t309;
t119 = qJ(2) + t219;
t51 = qJD(1) * t106 + qJDD(1) * t119;
t130 = t169 * qJDD(1) + qJDD(2);
t263 = qJD(3) * t166;
t64 = qJDD(3) * pkin(8) + t130 * t162 + t135 * t263;
t89 = t119 * qJD(1);
t118 = t162 * t135;
t96 = qJD(3) * pkin(8) + t118;
t228 = -t161 * t64 + t165 * t51 - t96 * t259 - t89 * t261;
t213 = qJDD(5) - t228;
t3 = pkin(9) * t44 - t316 * t105 + t213;
t264 = qJD(3) * t162;
t237 = t161 * t264;
t262 = qJD(4) * t109;
t45 = -qJD(1) * t237 - t165 * qJDD(3) + t161 * t247 + t262;
t127 = t140 * qJD(5);
t245 = -t161 * t51 - t165 * t64 - t89 * t259;
t191 = -t261 * t96 - t245;
t95 = t105 * qJ(5);
t9 = t127 + t191 + t95;
t5 = pkin(9) * t45 + t9;
t238 = -t160 * t5 + t164 * t3;
t277 = t165 * t166;
t283 = t161 * t166;
t82 = t160 * t277 - t164 * t283;
t333 = -g(1) * t207 + g(2) * t229 - g(3) * t82 + t23 * t56 - t238;
t330 = t56 * t201;
t41 = -t161 * t96 + t165 * t89;
t274 = qJD(5) - t41;
t329 = qJD(5) * t161 + t118;
t171 = qJD(1) ^ 2;
t184 = -qJ(2) * t171 + t328;
t252 = qJD(6) - t140;
t327 = t252 - qJD(6);
t254 = qJD(6) * t164;
t255 = qJD(6) * t160;
t12 = t107 * t254 - t109 * t255 + t160 * t45 - t164 * t44;
t325 = t201 * t252 + t12;
t310 = pkin(8) * t105;
t33 = pkin(4) * t107 - t182;
t324 = t140 * t33 - t310;
t291 = qJ(5) * t161;
t320 = -t316 * t165 - t291;
t301 = pkin(8) * qJD(4);
t242 = t140 * t301;
t319 = -t242 + t336;
t206 = t160 * t92 + t164 * t93;
t40 = t160 * t90 + t164 * t91;
t112 = t160 * t161 + t164 * t165;
t83 = t112 * t166;
t318 = -g(1) * t40 + g(2) * t206 - g(3) * t83 - t23 * t201;
t317 = t109 ^ 2;
t315 = pkin(8) - pkin(9);
t312 = pkin(4) * t105;
t311 = pkin(4) * t161;
t306 = g(3) * t166;
t186 = qJD(1) * t112;
t59 = t326 * t112;
t305 = t162 * t186 + t59;
t60 = t160 * t259 + t161 * t254 - t164 * t261 - t165 * t255;
t200 = t160 * t165 - t161 * t164;
t80 = t200 * t162;
t304 = qJD(1) * t80 + t60;
t303 = -t140 * t192 + t329;
t42 = t161 * t89 + t165 * t96;
t211 = -t290 + t311;
t302 = t140 * t211 - t329;
t129 = t140 * qJ(5);
t31 = t129 + t42;
t300 = t140 * t31;
t299 = t140 * t42;
t27 = pkin(9) * t107 + t42;
t20 = t129 + t27;
t298 = t160 * t20;
t296 = t161 * t44;
t114 = t220 * qJD(1);
t294 = t161 * t114 + t135 * t277;
t293 = pkin(1) * qJDD(1);
t292 = qJ(5) * t107;
t288 = t107 * t140;
t287 = t109 * t107;
t286 = t109 * t140;
t285 = t109 * t165;
t284 = t161 * t162;
t282 = t162 * t165;
t280 = t162 * t169;
t275 = -pkin(9) * t109 + t274;
t273 = g(2) * t166 * t276 + g(3) * t282;
t272 = t161 * t119 + t165 * t280;
t271 = t167 * pkin(1) + t163 * qJ(2);
t159 = t166 ^ 2;
t270 = t162 ^ 2 - t159;
t170 = qJD(3) ^ 2;
t269 = -t170 - t171;
t256 = qJD(5) * t165;
t250 = qJDD(1) * qJ(2);
t249 = qJDD(3) * t162;
t17 = -t316 * t140 + t275;
t246 = t160 * t3 + t164 * t5 + t17 * t254;
t244 = t166 * t308;
t235 = t169 * t263;
t243 = t161 * t106 + t119 * t259 + t165 * t235;
t241 = 0.2e1 * qJD(1) * qJD(2);
t48 = qJ(5) * t267 + t294;
t126 = t315 * t165;
t66 = t162 * qJ(5) + t272;
t233 = t169 * t260;
t230 = -t160 * t44 - t164 * t45;
t102 = t135 * t283;
t227 = -t114 * t165 + t102;
t132 = t161 * t280;
t226 = t119 * t165 - t132;
t225 = t252 ^ 2;
t224 = qJDD(2) - t293;
t222 = g(1) * t92 + g(2) * t90;
t221 = -g(1) * t93 - g(2) * t91;
t218 = g(1) * t167 + g(2) * t163;
t117 = t166 * t130;
t216 = qJDD(3) * pkin(3) - t135 * t264 + t117;
t125 = t315 * t161;
t215 = pkin(9) * t161 * t268 + qJD(6) * t125 - t315 * t261 - t48;
t185 = pkin(9) * t282 - t316 * t166;
t214 = qJD(1) * t185 - t326 * t126 + t227;
t212 = pkin(4) * t165 + t291;
t8 = t160 * t17 + t164 * t20;
t38 = t132 + (-pkin(9) * t166 - t119) * t165 - t316 * t162;
t50 = pkin(9) * t283 + t66;
t209 = -t160 * t50 + t164 * t38;
t208 = t160 * t38 + t164 * t50;
t30 = -pkin(4) * t140 + t274;
t205 = t161 * t31 - t165 * t30;
t204 = t161 * t30 + t165 * t31;
t203 = qJ(5) * t164 - t160 * t316;
t202 = -qJ(5) * t160 - t164 * t316;
t197 = -t119 * t261 - t161 * t235 + (t106 - t233) * t165;
t196 = pkin(3) + t212;
t193 = t20 * t255 - t246;
t190 = t105 * t161 + t140 * t259;
t189 = -t140 * t261 + t289;
t188 = t252 * t200;
t187 = 0.2e1 * qJ(2) * t251 + qJDD(3) * t169;
t181 = t140 * t97 - t310;
t180 = g(1) * t90 - g(2) * t92 + g(3) * t283 + t228;
t179 = t162 * t328 - t306;
t19 = qJ(5) * t263 + t162 * qJD(5) - t161 * t233 + t243;
t178 = -qJ(5) * t44 + qJD(5) * t109 + t216;
t177 = -t218 + t241 + 0.2e1 * t250;
t13 = qJD(6) * t56 + t230;
t10 = t213 - t312;
t176 = -qJD(4) * t205 + t10 * t161 + t9 * t165;
t175 = t109 * t33 + qJDD(5) - t180;
t174 = -t169 * t170 + t177;
t173 = g(1) * t91 - g(2) * t93 + g(3) * t277 + t140 * t41 - t191;
t172 = -t105 * t284 - t166 * t45 + t107 * t264 + (-t161 * t263 - t165 * t223) * t140;
t155 = t167 * qJ(2);
t150 = qJDD(3) * t166;
t139 = qJ(5) * t277;
t101 = pkin(3) - t320;
t100 = -qJDD(6) + t105;
t81 = t112 * t162;
t71 = -t139 + (-t169 + t311) * t166;
t67 = -pkin(4) * t162 - t226;
t65 = t139 + (t169 - t240) * t166;
t61 = pkin(4) * t109 + t292;
t49 = -pkin(4) * t267 + t227;
t34 = -t316 * t109 - t292;
t29 = (qJD(4) * t212 - t256) * t166 + (t169 - t211) * t264;
t25 = t288 - t44;
t24 = -pkin(4) * t263 - t197;
t22 = t326 * t166 * t200 - t112 * t264;
t21 = qJD(3) * t80 + t166 * t59;
t18 = (t320 * qJD(4) + t256) * t166 + (-t169 + t192) * t264;
t15 = (t165 * t258 - t237) * pkin(9) + t19;
t14 = pkin(9) * t234 + qJD(3) * t185 - t197;
t11 = pkin(4) * t45 - t178;
t7 = t164 * t17 - t298;
t6 = -t316 * t45 + t178;
t1 = [qJDD(1), -t328, t218, qJDD(2) + t328 - 0.2e1 * t293, t177, -t224 * pkin(1) - g(1) * (-pkin(1) * t163 + t155) - g(2) * t271 + (t241 + t250) * qJ(2), qJDD(1) * t159 - 0.2e1 * t162 * t232, -0.2e1 * t162 * t247 + 0.2e1 * t251 * t270, -t162 * t170 + t150, -t166 * t170 - t249, 0, t162 * t174 + t166 * t187, -t162 * t187 + t166 * t174, -t109 * t183 - t277 * t44 (t107 * t165 + t109 * t161) * t264 + (t296 - t165 * t45 + (t107 * t161 - t285) * qJD(4)) * t166 (-t140 * t253 - t44) * t162 + (t189 + t266) * t166 (t140 * t265 - t45) * t162 + (-qJD(3) * t107 - t190) * t166, t105 * t162 + t140 * t263, t197 * t140 + t226 * t105 + ((t107 * t169 - t161 * t97) * qJD(3) + t228) * t162 + (qJD(3) * t41 - t161 * t216 - t169 * t45 + t259 * t97) * t166 + t221, -t243 * t140 - t272 * t105 + ((t140 * t169 + t96) * t261 + (t109 * t169 - t165 * t97) * qJD(3) + t245) * t162 + (-qJD(3) * t42 - t165 * t216 + t169 * t44 - t261 * t97) * t166 + t222, -t105 * t67 + t107 * t29 - t140 * t24 + t45 * t71 + (-t265 * t33 - t10) * t162 + (-qJD(3) * t30 + t11 * t161 + t259 * t33) * t166 + t221, -t107 * t19 + t109 * t24 - t44 * t67 - t45 * t66 + t205 * t264 + (-qJD(4) * t204 + t10 * t165 - t161 * t9 + t218) * t166, t105 * t66 - t109 * t29 + t140 * t19 + t44 * t71 + (t253 * t33 + t9) * t162 + (qJD(3) * t31 - t11 * t165 + t261 * t33) * t166 - t222, t9 * t66 + t31 * t19 + t11 * t71 + t33 * t29 + t10 * t67 + t30 * t24 - g(1) * (pkin(3) * t281 + pkin(4) * t93 + qJ(5) * t92 - t167 * t309 + t155) - g(2) * (pkin(4) * t91 + pkin(7) * t167 + qJ(5) * t90 + t271) + (-g(1) * t169 - g(2) * t219) * t163, t12 * t83 + t22 * t56, -t12 * t82 - t13 * t83 - t201 * t22 + t21 * t56, -t100 * t83 - t12 * t162 + t22 * t252 - t263 * t56, t100 * t82 + t13 * t162 + t201 * t263 + t21 * t252, t100 * t162 - t252 * t263 (t14 * t164 - t15 * t160) * t252 - t209 * t100 - t238 * t162 - t7 * t263 + t18 * t201 + t65 * t13 + t6 * t82 - t23 * t21 - g(1) * t206 - g(2) * t40 + (t162 * t8 - t208 * t252) * qJD(6) -(qJD(6) * t209 + t14 * t160 + t15 * t164) * t252 + t208 * t100 - t193 * t162 + t8 * t263 + t18 * t56 + t65 * t12 + t6 * t83 + t23 * t22 + g(1) * t229 + g(2) * t207; 0, 0, 0, qJDD(1), -t171, t224 + t184, 0, 0, 0, 0, 0, t162 * t269 + t150, t166 * t269 - t249, 0, 0, 0, 0, 0, t172, t334, t172 (-t107 * t263 + t109 * t223 - t162 * t45) * t165 + (t107 * t223 + t109 * t263 - t162 * t44) * t161, -t334, -t205 * qJD(1) + (qJD(3) * t204 - t11) * t166 + (qJD(3) * t33 + t176) * t162 + t328, 0, 0, 0, 0, 0, t80 * t100 + t252 * t186 + (-qJD(3) * t188 + t13) * t166 + (-qJD(3) * t201 + t252 * t59) * t162, t81 * t100 - qJD(1) * t188 + (-qJD(3) * t112 * t252 + t12) * t166 + (-qJD(3) * t56 - t252 * t60) * t162; 0, 0, 0, 0, 0, 0, t166 * t171 * t162, -t270 * t171, t247, -t248, qJDD(3), t166 * t184 + t117 + t307, t306 + (-t130 - t184) * t162, t140 * t285 - t296 (-t44 - t288) * t165 + (-t45 - t286) * t161 (-t109 * t166 + t140 * t282) * qJD(1) + t190 (t107 * t166 - t140 * t284) * qJD(1) + t189, -t140 * t267, -t41 * t267 - t107 * t118 - pkin(3) * t45 + t102 * t140 + (-t244 + t216 + (-t114 - t301) * t140) * t165 + t181 * t161 + t273, pkin(3) * t44 + t294 * t140 + t42 * t267 - t109 * t118 + t181 * t165 + (-t319 - t216) * t161, t30 * t267 - t196 * t45 + t140 * t49 + t302 * t107 + (-t11 - t242 - t244) * t165 + t324 * t161 + t273, t107 * t48 - t109 * t49 + (t9 + t140 * t30 + (-t45 + t262) * pkin(8)) * t165 + (t10 - t300 + (qJD(4) * t107 - t44) * pkin(8)) * t161 + t179, -t31 * t267 - t196 * t44 - t140 * t48 - t302 * t109 - t324 * t165 + (-t11 + t319) * t161, -t30 * t49 - t31 * t48 + t302 * t33 + (t176 + t179) * pkin(8) + (-t11 + t336) * t196, -t12 * t200 + t305 * t56, -t112 * t12 + t13 * t200 - t201 * t305 - t304 * t56, t100 * t200 + t252 * t305 + t267 * t56, t100 * t112 - t201 * t267 - t252 * t304, t252 * t267 -(t125 * t164 - t126 * t160) * t100 + t101 * t13 + g(3) * t81 + t303 * t201 + t304 * t23 - (t160 * t215 + t164 * t214) * t252 + t7 * t267 + (t6 + t323) * t112 (t125 * t160 + t126 * t164) * t100 + t101 * t12 - t6 * t200 - g(3) * t80 + t303 * t56 + t305 * t23 - (-t160 * t214 + t164 * t215) * t252 + (-t8 * qJD(1) - t200 * t328) * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, -t107 ^ 2 + t317, t25, t286 - t45, t105, -t109 * t97 + t180 + t299, t107 * t97 + t173, -t107 * t61 - t175 + t299 + 0.2e1 * t312, pkin(4) * t44 - qJ(5) * t45 + (t31 - t42) * t109 + (t30 - t274) * t107, -t107 * t33 + t109 * t61 + 0.2e1 * t127 - t173 + 0.2e1 * t95, t9 * qJ(5) - t10 * pkin(4) - t33 * t61 - t30 * t42 - g(1) * (-pkin(4) * t90 + qJ(5) * t91) - g(2) * (pkin(4) * t92 - qJ(5) * t93) - g(3) * (-pkin(4) * t283 + t139) + t274 * t31, -t330, -t335, -t325, -t252 * t56 + t13, t100, -t202 * t100 - t34 * t201 - (t160 * t275 + t164 * t27) * t252 + (-t203 * t252 + t8) * qJD(6) + t333, t203 * t100 - t34 * t56 - (-t160 * t27 + t164 * t275) * t252 + (-t202 * t252 - t298) * qJD(6) + t246 + t318; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105 + t287, t25, -t140 ^ 2 - t317, t175 - t300 - t312, 0, 0, 0, 0, 0, -t164 * t100 - t109 * t201 - t160 * t225, t160 * t100 - t109 * t56 - t164 * t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t330, t335, t325, t327 * t56 - t230, -t100, t327 * t8 - t333, t252 * t7 + t193 - t318;];
tau_reg  = t1;
