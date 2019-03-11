% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:55:54
% EndTime: 2019-03-08 20:56:10
% DurationCPUTime: 5.70s
% Computational Cost: add. (8624->446), mult. (27987->695), div. (0->0), fcn. (25004->16), ass. (0->234)
t173 = sin(pkin(14));
t176 = sin(pkin(6));
t177 = cos(pkin(14));
t187 = cos(qJ(2));
t183 = sin(qJ(2));
t303 = cos(pkin(7));
t250 = t183 * t303;
t191 = t176 * (-t173 * t187 - t177 * t250);
t133 = qJD(1) * t191;
t174 = sin(pkin(8));
t178 = cos(pkin(8));
t175 = sin(pkin(7));
t282 = qJD(1) * t176;
t264 = t183 * t282;
t241 = t175 * t264;
t280 = qJD(3) * t175;
t299 = t173 * t174;
t213 = t133 * t174 - t178 * t241 + t280 * t299;
t186 = cos(qJ(4));
t248 = t303 * t174;
t232 = t186 * t248;
t289 = t178 * t186;
t267 = t177 * t289;
t329 = t175 * t267 + t232;
t192 = t176 * (-t173 * t250 + t177 * t187);
t137 = qJD(1) * t192;
t182 = sin(qJ(4));
t196 = t133 * t178 + t174 * t241;
t290 = t178 * t182;
t198 = t175 * (-t173 * t290 + t177 * t186);
t294 = t175 * t178;
t193 = (t177 * t294 + t248) * pkin(10);
t255 = t173 * t303;
t295 = t175 * t177;
t284 = pkin(2) * t255 + qJ(3) * t295;
t122 = t193 + t284;
t253 = t177 * t303;
t168 = pkin(2) * t253;
t266 = t303 * pkin(3);
t127 = t266 + t168 + (-pkin(10) * t178 - qJ(3)) * t175 * t173;
t214 = -pkin(3) * t177 - pkin(10) * t299;
t140 = (-pkin(2) + t214) * t175;
t216 = t127 * t178 + t140 * t174;
t317 = -t182 * t122 + t216 * t186;
t328 = qJD(3) * t198 + t317 * qJD(4) - t137 * t186 - t182 * t196;
t298 = t173 * t182;
t117 = (t232 + (t267 - t298) * t175) * qJD(4);
t126 = (t173 * t186 + t177 * t290) * t175 + t182 * t248;
t118 = t126 * qJD(4);
t327 = -pkin(4) * t118 + pkin(11) * t117 - t213;
t281 = qJD(2) * t175;
t153 = qJ(3) * t281 + t264;
t263 = t187 * t282;
t158 = qJD(2) * pkin(2) + t263;
t179 = cos(pkin(6));
t293 = t175 * t179;
t265 = qJD(1) * t293;
t95 = -t173 * t153 + t158 * t253 + t177 * t265;
t83 = (-t173 * pkin(10) * t294 + t266) * qJD(2) + t95;
t251 = t179 * t303;
t271 = qJD(1) * t251 + qJD(3);
t97 = (t214 * qJD(2) - t158) * t175 + t271;
t224 = t174 * t97 + t178 * t83;
t96 = t177 * t153 + t158 * t255 + t173 * t265;
t82 = qJD(2) * t193 + t96;
t322 = t182 * t82 - t224 * t186;
t116 = t126 * qJD(2);
t185 = cos(qJ(5));
t261 = t174 * t281;
t159 = t177 * t261;
t247 = qJD(2) * t303;
t209 = t178 * t247 - t159;
t200 = -qJD(4) - t209;
t139 = t185 * t200;
t181 = sin(qJ(5));
t85 = t116 * t181 + t139;
t84 = qJD(6) + t85;
t138 = qJD(2) * t198;
t278 = qJD(4) * t186;
t326 = -t174 * t278 + t138;
t262 = t173 * t281;
t325 = qJD(2) * t329 - t182 * t262;
t113 = qJD(5) - t325;
t189 = t122 * t186 + t216 * t182;
t197 = (t173 * t289 + t177 * t182) * t175;
t310 = qJD(3) * t197 + qJD(4) * t189 - t137 * t182 + t186 * t196;
t274 = qJD(5) * t185;
t276 = qJD(5) * t181;
t125 = t175 * t298 - t329;
t88 = -t127 * t174 + t178 * t140;
t52 = pkin(4) * t125 - pkin(11) * t126 + t88;
t252 = t178 * t303;
t146 = t174 * t295 - t252;
t57 = -pkin(11) * t146 + t189;
t324 = t327 * t181 - t328 * t185 - t52 * t274 + t57 * t276;
t323 = -t328 * t181 - t327 * t185;
t171 = t175 ^ 2;
t321 = t171 * (t173 ^ 2 + t177 ^ 2);
t109 = qJD(2) * t118;
t297 = t174 * t182;
t81 = t186 * t82;
t37 = t83 * t290 + t97 * t297 + t81;
t33 = -pkin(11) * t200 + t37;
t58 = -t174 * t83 + t178 * t97;
t38 = -pkin(4) * t325 - pkin(11) * t116 + t58;
t14 = t181 * t38 + t185 * t33;
t148 = (t263 + t280) * qJD(2);
t292 = t176 * t183;
t260 = qJD(2) * t292;
t235 = qJD(1) * t260;
t201 = t303 * t235;
t119 = -t173 * t148 - t177 * t201;
t120 = t177 * t148 - t173 * t201;
t215 = t175 * t235;
t205 = t174 * t215;
t316 = (t178 * t119 + t205) * t182 + t186 * t120;
t22 = -t322 * qJD(4) + t316;
t108 = t325 * qJD(4);
t92 = -t119 * t174 + t178 * t215;
t45 = pkin(4) * t109 - pkin(11) * t108 + t92;
t256 = t181 * t22 - t185 * t45;
t4 = -pkin(5) * t109 + t14 * qJD(5) + t256;
t87 = t185 * t116 - t181 * t200;
t319 = (pkin(5) * t87 + t84 * pkin(12)) * t84 + t4;
t318 = t224 * t182 + t81;
t61 = qJD(5) * t87 + t181 * t108;
t180 = sin(qJ(6));
t184 = cos(qJ(6));
t60 = -qJD(5) * t139 + t185 * t108 - t116 * t276;
t64 = t113 * t180 + t184 * t87;
t28 = t64 * qJD(6) - t184 * t109 + t180 * t60;
t249 = t187 * t303;
t124 = t177 * t292 + (t176 * t249 + t293) * t173;
t123 = t177 * t293 + (-t173 * t183 + t177 * t249) * t176;
t147 = -t176 * t187 * t175 + t251;
t217 = t123 * t178 + t147 * t174;
t315 = -t124 * t182 + t217 * t186;
t228 = -t119 * t289 + t182 * t120 - t186 * t205;
t23 = t318 * qJD(4) + t228;
t10 = pkin(5) * t61 - pkin(12) * t60 + t23;
t12 = pkin(12) * t113 + t14;
t32 = pkin(4) * t200 + t322;
t17 = t85 * pkin(5) - t87 * pkin(12) + t32;
t226 = t12 * t180 - t17 * t184;
t207 = t181 * t45 + t185 * t22 + t38 * t274 - t33 * t276;
t3 = pkin(12) * t109 + t207;
t1 = -t226 * qJD(6) + t10 * t180 + t184 * t3;
t219 = t181 * t52 + t185 * t57;
t314 = -pkin(5) * t118 + t219 * qJD(5) - t323;
t62 = -t184 * t113 + t180 * t87;
t313 = t62 * t84;
t312 = t64 * t84;
t78 = pkin(4) * t116 - pkin(11) * t325;
t311 = t181 * t78 - t185 * t322;
t309 = t113 * t85;
t272 = qJD(6) * t184;
t273 = qJD(6) * t180;
t27 = t180 * t109 + t113 * t272 + t184 * t60 - t87 * t273;
t308 = t180 * t27;
t307 = t180 * t61;
t306 = t184 * t61;
t305 = t87 * t113;
t304 = -t37 + t113 * (pkin(5) * t181 - pkin(12) * t185);
t302 = t325 * t185;
t301 = t119 * t173;
t296 = t174 * t186;
t188 = qJD(2) ^ 2;
t291 = t176 * t188;
t149 = -t178 * t185 + t181 * t297;
t240 = t173 * t261;
t286 = t149 * qJD(5) + t181 * t240 + t326 * t185;
t150 = t178 * t181 + t185 * t297;
t285 = t150 * qJD(5) - t326 * t181 + t185 * t240;
t279 = qJD(4) * t182;
t277 = qJD(5) * t180;
t275 = qJD(5) * t184;
t270 = t84 * t277;
t269 = t84 * t275;
t268 = t183 * t291;
t254 = t175 * t303;
t246 = t184 * t84;
t163 = -pkin(5) * t185 - pkin(12) * t181 - pkin(4);
t245 = pkin(12) * t116 - qJD(6) * t163 + t311;
t244 = t113 * t185;
t243 = t171 * t268;
t238 = t175 * t260;
t56 = pkin(4) * t146 - t317;
t90 = t126 * t181 + t185 * t146;
t91 = t126 * t185 - t146 * t181;
t29 = pkin(5) * t90 - pkin(12) * t91 + t56;
t236 = -pkin(12) * t118 - qJD(6) * t29 + t324;
t25 = pkin(12) * t125 + t219;
t65 = -t90 * qJD(5) + t117 * t185;
t66 = t91 * qJD(5) + t117 * t181;
t234 = -pkin(5) * t66 + pkin(12) * t65 + qJD(6) * t25 - t310;
t231 = t188 * t254;
t230 = qJD(3) * t254;
t6 = t12 * t184 + t17 * t180;
t225 = -t173 * t95 + t177 * t96;
t68 = t124 * t186 + t217 * t182;
t89 = -t123 * t174 + t147 * t178;
t42 = t181 * t89 + t185 * t68;
t223 = -t180 * t315 + t184 * t42;
t222 = -t180 * t42 - t184 * t315;
t13 = -t181 * t33 + t185 * t38;
t220 = -t181 * t57 + t185 * t52;
t41 = t181 * t68 - t185 * t89;
t218 = t184 * t125 - t180 * t91;
t70 = t125 * t180 + t184 * t91;
t212 = qJD(6) * t296 + t286;
t211 = -t84 * t272 - t307;
t210 = -t84 * t273 + t306;
t202 = -pkin(11) * t109 + t113 * t32;
t134 = qJD(2) * t197;
t199 = -qJD(6) * t150 + t174 * t279 - t134;
t135 = qJD(2) * t191;
t195 = t135 * t178 + t174 * t238;
t11 = -pkin(5) * t113 - t13;
t194 = -pkin(12) * t61 + (t11 + t13) * t84;
t2 = -t6 * qJD(6) + t184 * t10 - t180 * t3;
t136 = qJD(2) * t192;
t132 = -t158 * t175 + t271;
t102 = -t135 * t174 + t178 * t238;
t74 = t116 * t180 + t184 * t302;
t73 = -t184 * t116 + t180 * t302;
t40 = t315 * qJD(4) + t136 * t186 + t195 * t182;
t39 = qJD(4) * t68 + t136 * t182 - t186 * t195;
t31 = t70 * qJD(6) - t184 * t118 + t180 * t65;
t30 = t218 * qJD(6) + t118 * t180 + t184 * t65;
t24 = -pkin(5) * t125 - t220;
t19 = -pkin(5) * t116 - t181 * t322 - t185 * t78;
t16 = -t41 * qJD(5) + t102 * t181 + t185 * t40;
t15 = t42 * qJD(5) - t102 * t185 + t181 * t40;
t5 = [0, 0, -t268, -t187 * t291, t135 * t247 - t177 * t243, -t136 * t247 + t173 * t243 (-t135 * t173 + t136 * t177) * t281, t119 * t123 + t120 * t124 + t135 * t95 + t136 * t96 + (qJD(1) * t147 + t132) * t238, 0, 0, 0, 0, 0, -t102 * t325 + t89 * t109 + t200 * t39, t102 * t116 + t89 * t108 + t200 * t40, 0, 0, 0, 0, 0, -t109 * t41 - t113 * t15 - t315 * t61 + t39 * t85, -t109 * t42 - t113 * t16 - t315 * t60 + t39 * t87, 0, 0, 0, 0, 0 (-qJD(6) * t223 - t16 * t180 + t184 * t39) * t84 + t222 * t61 + t15 * t62 + t41 * t28 -(qJD(6) * t222 + t16 * t184 + t180 * t39) * t84 - t223 * t61 + t15 * t64 + t41 * t27; 0, 0, 0, 0, t119 * t303 + (-t303 * t133 - t173 * t230) * qJD(2), -t120 * t303 + (t303 * t137 - t177 * t230) * qJD(2) (t120 * t177 - t301) * t175 + ((t133 * t173 - t137 * t177) * t175 + qJD(3) * t321) * qJD(2), t120 * t284 + t119 * t168 - t171 * pkin(2) * t235 - t96 * t137 - t95 * t133 + (-qJ(3) * t301 + t225 * qJD(3) - t132 * t264) * t175, t108 * t126 + t116 * t117, -t108 * t125 - t109 * t126 - t116 * t118 + t117 * t325, -t146 * t108 - t117 * t200, t146 * t109 + t118 * t200, 0, t88 * t109 + t58 * t118 + t92 * t125 + t23 * t146 + t310 * t200 - t213 * t325, t88 * t108 + t213 * t116 + t58 * t117 + t92 * t126 + t22 * t146 + t328 * t200, t60 * t91 + t65 * t87, -t60 * t90 - t61 * t91 - t65 * t85 - t66 * t87, t109 * t91 + t113 * t65 + t118 * t87 + t125 * t60, -t109 * t90 - t113 * t66 - t118 * t85 - t125 * t61, t109 * t125 + t113 * t118, t220 * t109 - t256 * t125 + t13 * t118 + t56 * t61 + t23 * t90 + t32 * t66 + t310 * t85 + t323 * t113 + (-t219 * t113 - t14 * t125) * qJD(5), -t219 * t109 + t324 * t113 - t14 * t118 - t207 * t125 + t23 * t91 + t310 * t87 + t32 * t65 + t56 * t60, t27 * t70 + t30 * t64, t218 * t27 - t28 * t70 - t30 * t62 - t31 * t64, t27 * t90 + t30 * t84 + t61 * t70 + t64 * t66, t218 * t61 - t28 * t90 - t31 * t84 - t62 * t66, t61 * t90 + t66 * t84 (-t180 * t25 + t184 * t29) * t61 + t2 * t90 - t226 * t66 + t24 * t28 - t4 * t218 + t11 * t31 + (t180 * t236 - t184 * t234) * t84 + t314 * t62 -(t180 * t29 + t184 * t25) * t61 - t1 * t90 - t6 * t66 + t24 * t27 + t4 * t70 + t11 * t30 + (t180 * t234 + t184 * t236) * t84 + t314 * t64; 0, 0, 0, 0, t173 * t231, t177 * t231, -t188 * t321 (-t225 + t264) * t281, 0, 0, 0, 0, 0, t178 * t109 - t134 * t200 + (t200 * t279 + t262 * t325) * t174, t178 * t108 - t138 * t200 + (-t116 * t262 + t200 * t278) * t174, 0, 0, 0, 0, 0, -t109 * t149 - t134 * t85 + (-t186 * t61 + t85 * t279) * t174 - t285 * t113, -t109 * t150 - t134 * t87 + (-t186 * t60 + t87 * t279) * t174 + t286 * t113, 0, 0, 0, 0, 0 (-t150 * t180 - t184 * t296) * t61 + t149 * t28 + (t180 * t212 + t184 * t199) * t84 + t285 * t62 -(t150 * t184 - t180 * t296) * t61 + t149 * t27 + (-t180 * t199 + t184 * t212) * t84 + t285 * t64; 0, 0, 0, 0, 0, 0, 0, 0, -t325 * t116, t116 ^ 2 - t325 ^ 2, t200 * t325 + t108 (qJD(4) - t159) * t116 + (t116 * t252 - t118) * qJD(2), 0, -t58 * t116 + t37 * t209 + (-t318 + t37) * qJD(4) - t228, -t322 * t209 - t325 * t58 - t316, t181 * t60 + t87 * t244 (t60 - t309) * t185 + (-t61 - t305) * t181, t181 * t109 + t113 * t244 - t116 * t87, -t113 ^ 2 * t181 + t185 * t109 + t116 * t85, -t113 * t116, -pkin(4) * t61 - t13 * t116 - t37 * t85 + (-t23 + (-pkin(11) * qJD(5) - t78) * t113) * t185 + (-t113 * t322 + t202) * t181, -pkin(4) * t60 + t14 * t116 + t23 * t181 - t37 * t87 + (pkin(11) * t276 + t311) * t113 + t202 * t185, t181 * t184 * t27 + (-t181 * t273 + t184 * t274 - t74) * t64, t62 * t74 + t64 * t73 + (-t180 * t64 - t184 * t62) * t274 + (-t308 - t184 * t28 + (t180 * t62 - t184 * t64) * qJD(6)) * t181, -t74 * t84 + (-t27 + t269) * t185 + (t113 * t64 + t210) * t181, t73 * t84 + (t28 - t270) * t185 + (-t113 * t62 + t211) * t181, t113 * t181 * t84 - t185 * t61, t163 * t306 - t11 * t73 - t19 * t62 + (t180 * t245 + t184 * t304) * t84 + (t11 * t277 - t2 + (qJD(5) * t62 + t211) * pkin(11)) * t185 + (t11 * t272 + t4 * t180 - t113 * t226 + (t28 + t270) * pkin(11)) * t181, -t163 * t307 - t11 * t74 - t19 * t64 + (-t180 * t304 + t184 * t245) * t84 + (t11 * t275 + t1 + (qJD(5) * t64 - t210) * pkin(11)) * t185 + (-t11 * t273 + t4 * t184 - t113 * t6 + (t27 + t269) * pkin(11)) * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87 * t85, -t85 ^ 2 + t87 ^ 2, t60 + t309, t305 - t61, t109, -t32 * t87 - t256 + (-qJD(5) + t113) * t14, t113 * t13 + t32 * t85 - t207, t64 * t246 + t308 (t27 - t313) * t184 + (-t28 - t312) * t180, t246 * t84 - t64 * t87 + t307, -t180 * t84 ^ 2 + t62 * t87 + t306, -t84 * t87, -pkin(5) * t28 - t14 * t62 + t194 * t180 - t184 * t319 + t226 * t87, -pkin(5) * t27 - t14 * t64 + t180 * t319 + t194 * t184 + t6 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t62, -t62 ^ 2 + t64 ^ 2, t27 + t313, -t28 + t312, t61, -t11 * t64 + t6 * t84 + t2, t11 * t62 - t226 * t84 - t1;];
tauc_reg  = t5;
