% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:17:57
% EndTime: 2019-03-09 00:18:06
% DurationCPUTime: 4.35s
% Computational Cost: add. (5679->438), mult. (13912->583), div. (0->0), fcn. (10221->10), ass. (0->219)
t176 = sin(qJ(3));
t179 = cos(qJ(3));
t207 = pkin(3) * t176 - pkin(9) * t179;
t132 = t207 * qJD(3);
t175 = sin(qJ(4));
t177 = sin(qJ(2));
t178 = cos(qJ(4));
t252 = qJD(3) * t176;
t172 = sin(pkin(6));
t258 = qJD(1) * t172;
t180 = cos(qJ(2));
t266 = t179 * t180;
t298 = pkin(8) * t175;
t322 = t178 * t132 + t252 * t298 - (-t175 * t266 + t177 * t178) * t258;
t138 = -pkin(3) * t179 - pkin(9) * t176 - pkin(2);
t247 = qJD(4) * t178;
t321 = -(t175 * t177 + t178 * t266) * t258 + t175 * t132 + t138 * t247;
t268 = t178 * t179;
t161 = pkin(8) * t268;
t205 = pkin(4) * t176 - pkin(10) * t268;
t190 = t205 * qJD(3);
t320 = t190 + (-t161 + (pkin(10) * t176 - t138) * t175) * qJD(4) + t322;
t249 = qJD(4) * t175;
t251 = qJD(3) * t178;
t250 = qJD(3) * t179;
t227 = t175 * t250;
t229 = t176 * t247;
t314 = t227 + t229;
t319 = -t314 * pkin(10) + (-t176 * t251 - t179 * t249) * pkin(8) + t321;
t254 = qJD(2) * t179;
t232 = t175 * t254;
t302 = -pkin(10) - pkin(9);
t238 = qJD(4) * t302;
t237 = t177 * t258;
t135 = qJD(2) * pkin(8) + t237;
t173 = cos(pkin(6));
t257 = qJD(1) * t179;
t102 = -t176 * t135 + t173 * t257;
t129 = t207 * qJD(2);
t280 = t178 * t102 + t175 * t129;
t318 = pkin(10) * t232 + t175 * t238 - t280;
t219 = -t102 * t175 + t178 * t129;
t317 = -qJD(2) * t205 + t178 * t238 - t219;
t255 = qJD(2) * t176;
t228 = t175 * t255;
t125 = -t228 + t251;
t253 = qJD(3) * t175;
t126 = t178 * t255 + t253;
t174 = sin(qJ(5));
t299 = cos(qJ(5));
t194 = t174 * t125 + t299 * t126;
t288 = qJD(3) * pkin(3);
t94 = -t102 - t288;
t65 = -pkin(4) * t125 + t94;
t72 = -t299 * t125 + t126 * t174;
t24 = pkin(5) * t72 - qJ(6) * t194 + t65;
t316 = t24 * t72;
t315 = t65 * t72;
t294 = t194 * t72;
t272 = t174 * t175;
t192 = t299 * t178 - t272;
t306 = qJD(4) + qJD(5);
t225 = t299 * qJD(5);
t307 = t299 * qJD(4) + t225;
t282 = t178 * t307 - t192 * t254 - t272 * t306;
t128 = t174 * t178 + t299 * t175;
t81 = t306 * t128;
t281 = -t128 * t254 + t81;
t303 = t194 ^ 2;
t313 = -t72 ^ 2 + t303;
t159 = -qJD(4) + t254;
t147 = -qJD(5) + t159;
t248 = qJD(4) * t176;
t230 = t175 * t248;
t189 = t178 * t250 - t230;
t243 = qJD(3) * qJD(4);
t223 = t178 * t243;
t185 = qJD(2) * t189 + t223;
t224 = qJD(2) * t248;
t239 = qJD(2) * t227 + t175 * t243 + t178 * t224;
t246 = qJD(5) * t174;
t26 = -t125 * t225 + t126 * t246 + t174 * t239 - t299 * t185;
t21 = -t147 * t72 - t26;
t43 = pkin(5) * t194 + qJ(6) * t72;
t295 = t24 * t194;
t311 = t65 * t194;
t124 = t178 * t138;
t269 = t176 * t178;
t77 = -pkin(10) * t269 + t124 + (-pkin(4) - t298) * t179;
t260 = t175 * t138 + t161;
t270 = t175 * t176;
t86 = -pkin(10) * t270 + t260;
t310 = t174 * t320 + t77 * t225 - t246 * t86 + t319 * t299;
t142 = t302 * t175;
t143 = t302 * t178;
t193 = t299 * t142 + t174 * t143;
t309 = t193 * qJD(5) + t174 * t317 + t299 * t318;
t98 = t174 * t142 - t299 * t143;
t308 = t98 * qJD(5) + t174 * t318 - t299 * t317;
t236 = t180 * t258;
t105 = qJD(2) * t138 - t236;
t271 = t175 * t105;
t273 = t173 * t176;
t154 = qJD(1) * t273;
t103 = t179 * t135 + t154;
t95 = qJD(3) * pkin(9) + t103;
t58 = t178 * t95 + t271;
t42 = pkin(10) * t125 + t58;
t209 = -t103 + (-t232 + t249) * pkin(4);
t27 = qJD(5) * t194 + t174 * t185 + t299 * t239;
t305 = -t147 * t194 - t27;
t198 = t174 * t77 + t299 * t86;
t304 = qJD(5) * t198 + t319 * t174 - t299 * t320;
t301 = -pkin(5) * t252 + t304;
t300 = qJ(6) * t252 - qJD(6) * t179 + t310;
t297 = pkin(9) * t159;
t293 = qJ(6) * t255 - t309;
t292 = pkin(5) * t255 + t308;
t291 = t281 * pkin(5) - t282 * qJ(6) - qJD(6) * t128 + t209;
t289 = qJD(2) * pkin(2);
t287 = t174 * t42;
t286 = t175 * t94;
t211 = t176 * t236;
t67 = qJD(2) * t211 + qJD(3) * t154 + t135 * t250;
t285 = t67 * t175;
t284 = t67 * t178;
t57 = t178 * t105 - t175 * t95;
t41 = -pkin(10) * t126 + t57;
t20 = t299 * t41 - t287;
t283 = pkin(4) * t225 + qJD(6) - t20;
t279 = qJD(3) * t193;
t278 = qJD(3) * t98;
t277 = t126 * t159;
t276 = t172 * t177;
t275 = t172 * t180;
t182 = qJD(2) ^ 2;
t274 = t172 * t182;
t267 = t179 * t159;
t181 = qJD(3) ^ 2;
t265 = t181 * t176;
t264 = t181 * t179;
t32 = -pkin(4) * t159 + t41;
t14 = t299 * t32 - t287;
t263 = qJD(6) - t14;
t133 = pkin(4) * t270 + t176 * pkin(8);
t170 = t176 ^ 2;
t259 = -t179 ^ 2 + t170;
t256 = qJD(2) * t172;
t245 = t170 * qJD(2);
t244 = qJD(2) * qJD(3);
t241 = t299 * t42;
t240 = t177 * t274;
t104 = t314 * pkin(4) + pkin(8) * t250;
t167 = -pkin(4) * t178 - pkin(3);
t234 = t177 * t256;
t233 = t180 * t256;
t231 = t159 * t249;
t164 = t176 * t244;
t101 = (t132 + t237) * qJD(2);
t66 = -t135 * t252 + (qJD(3) * t173 + t233) * t257;
t222 = -t178 * t101 + t175 * t66;
t13 = qJD(2) * t190 - qJD(4) * t42 - t222;
t196 = t175 * t101 + t105 * t247 + t178 * t66 - t249 * t95;
t18 = -pkin(10) * t239 + t196;
t221 = -t174 * t13 - t299 * t18 - t32 * t225 + t42 * t246;
t220 = -t299 * t13 + t174 * t18 + t42 * t225 + t32 * t246;
t218 = t159 + t254;
t217 = -t125 + t251;
t216 = qJD(4) + t254;
t215 = pkin(5) * t164;
t214 = t194 * t236;
t213 = t176 * t233;
t212 = t179 * t233;
t210 = t299 * t250;
t19 = t174 * t41 + t241;
t208 = pkin(4) * t246 - t19;
t136 = -t236 - t289;
t206 = -t136 - t236;
t137 = t147 * qJD(6);
t155 = qJ(6) * t164;
t1 = t155 - t137 - t221;
t204 = -t14 * t147 + t221;
t15 = t174 * t32 + t241;
t203 = -t147 * t15 - t220;
t201 = -t174 * t86 + t299 * t77;
t115 = t179 * t276 + t273;
t197 = -t115 * t178 + t175 * t275;
t84 = -t115 * t175 - t178 * t275;
t199 = t174 * t197 + t299 * t84;
t47 = t174 * t84 - t197 * t299;
t114 = -t173 * t179 + t176 * t276;
t191 = t216 * t253;
t49 = pkin(4) * t239 + t67;
t2 = -t215 + t220;
t82 = -qJD(3) * t114 + t212;
t36 = qJD(4) * t197 - t175 * t82 + t178 * t234;
t37 = qJD(4) * t84 + t175 * t234 + t178 * t82;
t5 = t47 * qJD(5) + t174 * t37 - t299 * t36;
t83 = qJD(3) * t115 + t213;
t188 = t114 * t27 + t147 * t5 + t164 * t199 + t83 * t72;
t187 = qJD(3) * (-t206 - t289);
t4 = t199 * qJD(5) + t174 * t36 + t299 * t37;
t186 = t114 * t26 - t147 * t4 + t164 * t47 - t194 * t83;
t166 = -t299 * pkin(4) - pkin(5);
t162 = pkin(4) * t174 + qJ(6);
t109 = t192 * t176;
t108 = t128 * t176;
t70 = -pkin(5) * t192 - qJ(6) * t128 + t167;
t62 = pkin(5) * t108 - qJ(6) * t109 + t133;
t52 = t175 * t210 - t174 * t230 - t246 * t270 + (t174 * t250 + t176 * t307) * t178;
t51 = t174 * t227 + t176 * t81 - t178 * t210;
t40 = t179 * pkin(5) - t201;
t39 = -qJ(6) * t179 + t198;
t29 = pkin(4) * t126 + t43;
t10 = -t147 * qJ(6) + t15;
t9 = pkin(5) * t52 + qJ(6) * t51 - qJD(6) * t109 + t104;
t8 = t147 * pkin(5) + t263;
t3 = t27 * pkin(5) + t26 * qJ(6) - qJD(6) * t194 + t49;
t6 = [0, 0, -t240, -t180 * t274, 0, 0, 0, 0, 0, -t179 * t240 + (-t83 - t213) * qJD(3), t176 * t240 + (-t82 - t212) * qJD(3), 0, 0, 0, 0, 0, t114 * t239 - t83 * t125 - t36 * t159 + t164 * t84, t114 * t223 + t83 * t126 + t37 * t159 + (t114 * t189 + t197 * t252) * qJD(2), 0, 0, 0, 0, 0, t188, -t186, t188, t194 * t5 + t199 * t26 - t27 * t47 - t4 * t72, t186, t1 * t47 + t10 * t4 + t114 * t3 - t199 * t2 + t24 * t83 + t5 * t8; 0, 0, 0, 0, 0.2e1 * t179 * t164, -0.2e1 * t259 * t244, t264, -t265, 0, -pkin(8) * t264 + t176 * t187, pkin(8) * t265 + t179 * t187, t126 * t189 + t185 * t269 (t178 * t125 - t126 * t175) * t250 + ((-t125 + t228) * t249 + (-t126 * qJD(4) - t191 - t239) * t178) * t176, t218 * t230 + (t126 * t176 + (t245 + (-t159 - t216) * t179) * t178) * qJD(3), t159 * t229 + t239 * t179 + (t125 * t176 + (-t245 + t267) * t175) * qJD(3), -t218 * t252 (t138 * t249 - t322) * t159 + ((-pkin(8) * t125 + t286) * qJD(3) + (t271 + (pkin(8) * t159 + t95) * t178) * qJD(4) + t222) * t179 + (pkin(8) * t239 + t285 + t94 * t247 + t125 * t236 + ((-t179 * t298 + t124) * qJD(2) + t57) * qJD(3)) * t176, t321 * t159 + (t94 * t251 + (qJD(3) * t126 - t231) * pkin(8) + t196) * t179 + (-t126 * t236 + t284 + (-pkin(8) * t255 - t94) * t249 + (-t260 * qJD(2) - t58 + (-t159 + t216) * pkin(8) * t178) * qJD(3)) * t176, -t109 * t26 - t194 * t51, t108 * t26 - t109 * t27 - t194 * t52 + t51 * t72, t147 * t51 + t179 * t26 + (qJD(2) * t109 + t194) * t252, t147 * t52 + t179 * t27 + (-qJD(2) * t108 - t72) * t252 (-t147 - t254) * t252, t49 * t108 + t133 * t27 + t14 * t252 + t201 * t164 + t220 * t179 + t65 * t52 + (t104 - t211) * t72 + t304 * t147, -t221 * t179 + t104 * t194 - t133 * t26 + t49 * t109 - t65 * t51 + t310 * t147 + (-t214 + (-qJD(2) * t198 - t15) * qJD(3)) * t176, t108 * t3 + t179 * t2 + t24 * t52 + t27 * t62 + t72 * t9 + t301 * t147 + (-t72 * t236 + (-qJD(2) * t40 - t8) * qJD(3)) * t176, -t1 * t108 - t10 * t52 + t109 * t2 + t194 * t301 - t26 * t40 - t27 * t39 - t300 * t72 - t51 * t8, -t1 * t179 - t109 * t3 + t24 * t51 + t26 * t62 - t194 * t9 - t300 * t147 + (t214 + (qJD(2) * t39 + t10) * qJD(3)) * t176, t1 * t39 + t2 * t40 + t3 * t62 + t301 * t8 + (t9 - t211) * t24 + t300 * t10; 0, 0, 0, 0, -t176 * t182 * t179, t259 * t182, 0, 0, 0, qJD(3) * t103 - t136 * t255 - t67, t206 * t254, -t175 ^ 2 * t224 + (t191 - t277) * t178 (-t239 + t277) * t175 + ((t125 + t251) * qJD(4) + (t179 * t217 - t230) * qJD(2)) * t178, -t159 * t247 + (t178 * t267 + (-t126 + t253) * t176) * qJD(2), t231 + (-t175 * t267 + t176 * t217) * qJD(2), t159 * t255, -pkin(3) * t239 - t284 + t219 * t159 + t103 * t125 + (t178 * t297 + t286) * qJD(4) + (-t57 * t176 + (-pkin(9) * t252 - t179 * t94) * t175) * qJD(2), t285 - t280 * t159 - t103 * t126 + (-t175 * t297 + (t94 - t288) * t178) * qJD(4) + ((-t94 - t288) * t268 + (pkin(3) * t249 - pkin(9) * t251 + t58) * t176) * qJD(2), -t128 * t26 + t194 * t282, -t128 * t27 - t192 * t26 - t194 * t281 - t282 * t72, -t282 * t147 + (qJD(3) * t128 - t194) * t255, t281 * t147 + (qJD(3) * t192 + t72) * t255, t147 * t255, -t49 * t192 + t167 * t27 + t209 * t72 + t281 * t65 + t308 * t147 + (-t14 + t279) * t255, t49 * t128 - t167 * t26 + t209 * t194 + t282 * t65 + t309 * t147 + (t15 - t278) * t255, -t192 * t3 + t27 * t70 + t291 * t72 + t281 * t24 + t292 * t147 + (t8 + t279) * t255, t1 * t192 - t281 * t10 + t128 * t2 + t193 * t26 + t194 * t292 - t27 * t98 + t282 * t8 + t293 * t72, -t128 * t3 + t26 * t70 - t291 * t194 - t282 * t24 + t293 * t147 + (-t10 + t278) * t255, t1 * t98 - t10 * t293 - t193 * t2 + t291 * t24 + t292 * t8 + t3 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 * t125, -t125 ^ 2 + t126 ^ 2, t125 * t159 + t185, -t239 - t277, t164, -t126 * t94 - t222 + (-qJD(4) - t159) * t58, -t125 * t94 - t159 * t57 - t196, t294, t313, t21, t305, t164, -t19 * t147 - t311 + (-t126 * t72 + t147 * t246 + t299 * t164) * pkin(4) - t220, -t20 * t147 + t315 + (-t126 * t194 + t147 * t225 - t164 * t174) * pkin(4) + t221, -t295 - t29 * t72 + t208 * t147 + (pkin(5) - t166) * t164 - t220, -t162 * t27 - t166 * t26 + (t10 + t208) * t194 + (-t283 + t8) * t72, -t147 * t283 + t162 * t164 + t194 * t29 + t1 - t316, t1 * t162 + t10 * t283 + t166 * t2 + t208 * t8 - t24 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t294, t313, t21, t305, t164, t203 - t311, t204 + t315, -t43 * t72 + t203 + 0.2e1 * t215 - t295, pkin(5) * t26 - qJ(6) * t27 + (t10 - t15) * t194 + (t8 - t263) * t72, t194 * t43 - 0.2e1 * t137 + 0.2e1 * t155 - t204 - t316, -pkin(5) * t2 + qJ(6) * t1 + t10 * t263 - t15 * t8 - t24 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164 + t294, t21, -t147 ^ 2 - t303, t10 * t147 + t2 + t295;];
tauc_reg  = t6;
