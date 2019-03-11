% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:04:00
% EndTime: 2019-03-09 22:04:11
% DurationCPUTime: 4.17s
% Computational Cost: add. (7813->384), mult. (19828->495), div. (0->0), fcn. (14831->8), ass. (0->219)
t182 = sin(qJ(2));
t184 = cos(qJ(2));
t297 = sin(qJ(3));
t249 = qJD(1) * t297;
t299 = cos(qJ(3));
t250 = qJD(1) * t299;
t139 = -t182 * t249 + t184 * t250;
t140 = -t182 * t250 - t184 * t249;
t181 = sin(qJ(4));
t298 = cos(qJ(4));
t113 = t298 * t139 + t140 * t181;
t180 = sin(qJ(6));
t177 = qJD(2) + qJD(3);
t209 = t297 * t182 - t299 * t184;
t194 = t177 * t209;
t190 = t298 * t194;
t151 = t299 * t182 + t297 * t184;
t196 = t177 * t151;
t193 = t196 * qJD(1);
t247 = qJD(4) * t298;
t260 = qJD(4) * t181;
t215 = -qJD(1) * t190 + t139 * t247 + t140 * t260 - t181 * t193;
t216 = t181 * t139 - t298 * t140;
t305 = qJD(6) + t216;
t183 = cos(qJ(6));
t315 = t305 * t183;
t208 = -t180 * t215 - t305 * t315;
t175 = qJD(4) + t177;
t92 = t113 * t183 + t175 * t180;
t14 = t92 * t113 + t208;
t316 = t305 * t180;
t46 = t183 * t215;
t218 = -t305 * t316 + t46;
t94 = -t113 * t180 + t175 * t183;
t13 = -t94 * t113 + t218;
t274 = qJ(5) * t113;
t78 = pkin(4) * t216 - t274;
t323 = qJD(6) - t305;
t258 = qJD(6) * t183;
t259 = qJD(6) * t180;
t191 = t298 * t196;
t192 = t194 * qJD(1);
t189 = qJD(1) * t191 - t181 * t192;
t53 = qJD(4) * t216 + t189;
t25 = -t113 * t258 - t175 * t259 + t180 * t53;
t24 = t25 * t183;
t16 = -t316 * t94 + t24;
t47 = t183 * t53;
t26 = qJD(6) * t94 - t47;
t321 = t305 * t92;
t1 = -t315 * t94 + (-t25 + t321) * t180 - t183 * t26;
t322 = t113 * pkin(5);
t136 = t140 * pkin(9);
t301 = pkin(7) + pkin(8);
t159 = t301 * t184;
t154 = qJD(1) * t159;
t141 = t297 * t154;
t158 = t301 * t182;
t152 = qJD(1) * t158;
t289 = qJD(2) * pkin(2);
t147 = -t152 + t289;
t226 = t299 * t147 - t141;
t90 = t136 + t226;
t83 = t177 * pkin(3) + t90;
t145 = t299 * t154;
t211 = t297 * t147 + t145;
t291 = t139 * pkin(9);
t91 = t211 + t291;
t89 = t298 * t91;
t55 = t181 * t83 + t89;
t48 = -qJ(5) * t175 - t55;
t31 = -t48 + t322;
t243 = t305 * t31;
t252 = qJD(2) * t301;
t232 = qJD(1) * t252;
t148 = t182 * t232;
t149 = t184 * t232;
t246 = t297 * qJD(3);
t248 = qJD(3) * t299;
t201 = t147 * t248 - t299 * t148 - t297 * t149 - t154 * t246;
t60 = -pkin(9) * t193 + t201;
t229 = t297 * t148 - t299 * t149;
t61 = pkin(9) * t192 - t147 * t246 - t154 * t248 + t229;
t12 = t181 * t60 + t91 * t247 + t83 * t260 - t298 * t61;
t172 = -t184 * pkin(2) - pkin(1);
t157 = qJD(1) * t172;
t124 = -t139 * pkin(3) + t157;
t202 = -qJ(5) * t216 + t124;
t64 = -pkin(4) * t113 + t202;
t288 = t216 * t64;
t222 = t12 + t288;
t171 = t299 * pkin(2) + pkin(3);
t227 = t297 * t152 - t145;
t96 = t227 - t291;
t264 = -t299 * t152 - t141;
t97 = t136 + t264;
t318 = t181 * t96 + t298 * t97 - t171 * t247 + ((t297 * qJD(4) + t246) * t181 - t298 * t248) * pkin(2);
t317 = t216 * t113;
t276 = t305 * t113;
t88 = t181 * t91;
t54 = -t298 * t83 + t88;
t265 = qJD(5) + t54;
t303 = t216 ^ 2;
t32 = -t113 ^ 2 + t303;
t293 = t216 * pkin(5);
t266 = t293 + t265;
t302 = pkin(4) + pkin(10);
t27 = -t302 * t175 + t266;
t38 = -t113 * t302 + t202;
t18 = t180 * t27 + t183 * t38;
t244 = -t181 * t61 - t83 * t247 + t91 * t260 - t298 * t60;
t11 = -t175 * qJD(5) + t244;
t4 = -pkin(5) * t53 - t11;
t251 = t18 * t113 + t4 * t183;
t213 = t64 * t113 - t11;
t219 = -t124 * t113 + t244;
t314 = -t113 * t175 + t215;
t224 = t180 * t38 - t183 * t27;
t313 = t224 * t113 + t4 * t180 + t31 * t258;
t255 = qJD(1) * qJD(2);
t311 = -0.2e1 * t255;
t309 = t216 * t31;
t281 = -qJD(5) + t318;
t233 = t298 * t297;
t279 = -t181 * t97 + t298 * t96 + t171 * t260 + (qJD(4) * t233 + (t299 * t181 + t233) * qJD(3)) * pkin(2);
t63 = t298 * t90 - t88;
t277 = pkin(3) * t247 + qJD(5) - t63;
t62 = t181 * t90 + t89;
t230 = pkin(3) * t260 - t62;
t308 = t124 * t216;
t307 = t279 * t175;
t30 = -t139 * t260 + t140 * t247 + t175 * t216 - t189;
t295 = pkin(3) * t140;
t203 = t298 * t209;
t122 = t181 * t151 + t203;
t205 = t181 * t209;
t123 = t298 * t151 - t205;
t127 = pkin(3) * t209 + t172;
t200 = -t123 * qJ(5) + t127;
t43 = t302 * t122 + t200;
t290 = t43 * t215;
t103 = -t151 * pkin(9) - t299 * t158 - t297 * t159;
t210 = t297 * t158 - t299 * t159;
t104 = -pkin(9) * t209 - t210;
t153 = t182 * t252;
t155 = t184 * t252;
t206 = -t299 * t153 - t297 * t155 - t158 * t248 - t159 * t246;
t73 = -pkin(9) * t196 + t206;
t228 = t297 * t153 - t299 * t155;
t74 = pkin(9) * t194 + t158 * t246 - t159 * t248 + t228;
t20 = -t103 * t247 + t104 * t260 - t181 * t74 - t298 * t73;
t287 = t175 * t20;
t217 = t181 * t103 + t298 * t104;
t21 = t217 * qJD(4) + t181 * t73 - t298 * t74;
t286 = t175 * t21;
t284 = t180 * t94;
t283 = t302 * t215;
t282 = t293 - t281;
t280 = -t322 + t279;
t278 = t293 + t277;
t271 = t122 * t180;
t270 = t140 * t139;
t187 = qJD(1) ^ 2;
t269 = t184 * t187;
t186 = qJD(2) ^ 2;
t268 = t186 * t182;
t267 = t186 * t184;
t262 = t182 ^ 2 - t184 ^ 2;
t261 = qJD(1) * t182;
t254 = t297 * pkin(2);
t174 = t182 * t289;
t173 = pkin(2) * t261;
t245 = t182 * t255;
t242 = pkin(1) * t311;
t105 = t216 * pkin(10);
t134 = -t298 * t171 + t181 * t254 - pkin(4);
t130 = -pkin(10) + t134;
t70 = -t295 + t78;
t67 = t173 + t70;
t241 = -qJD(6) * t130 + t105 + t67;
t170 = -t298 * pkin(3) - pkin(4);
t168 = -pkin(10) + t170;
t240 = -qJD(6) * t168 + t105 + t70;
t239 = t305 * t302 - t274;
t165 = pkin(2) * t245;
t231 = -t322 + t230;
t76 = -t298 * t103 + t181 * t104;
t44 = -pkin(4) * t175 + t265;
t225 = -t44 * t113 - t216 * t48;
t223 = t183 * t309 + t313;
t221 = t175 * t55 - t12;
t220 = -t12 - t308;
t49 = t123 * pkin(5) + t76;
t69 = -qJD(4) * t205 + t151 * t247 - t181 * t194 + t191;
t214 = t4 * t122 - t215 * t49 + t31 * t69;
t212 = t157 * t140 + t229;
t199 = -t157 * t139 - t201;
t116 = pkin(3) * t196 + t174;
t98 = pkin(3) * t193 + t165;
t68 = qJD(4) * t203 + t151 * t260 + t181 * t196 + t190;
t19 = t69 * pkin(4) + t68 * qJ(5) - t123 * qJD(5) + t116;
t15 = t53 * pkin(4) - qJ(5) * t215 - qJD(5) * t216 + t98;
t169 = pkin(3) * t181 + qJ(5);
t133 = pkin(2) * t233 + t181 * t171 + qJ(5);
t125 = t173 - t295;
t95 = -t139 ^ 2 + t140 ^ 2;
t87 = -t140 * t177 - t193;
t86 = -t139 * t177 - t192;
t75 = t122 * pkin(4) + t200;
t50 = -t122 * pkin(5) + t217;
t35 = t215 * t123;
t34 = t55 + t322;
t10 = -t68 * pkin(5) + t21;
t9 = -pkin(5) * t69 - t20;
t8 = t69 * pkin(10) + t19;
t7 = t53 * pkin(10) + t15;
t6 = pkin(5) * t215 + t12;
t5 = t183 * t6;
t2 = [0, 0, 0, 0.2e1 * t184 * t245, t262 * t311, t267, -t268, 0, -pkin(7) * t267 + t182 * t242, pkin(7) * t268 + t184 * t242, t140 * t194 - t151 * t192 (qJD(1) * (-t151 ^ 2 + t209 ^ 2) - t139 * t209 + t140 * t151) * t177, -t194 * t177, -t196 * t177, 0, -t139 * t174 + t209 * t165 + (qJD(3) * t210 + 0.2e1 * t157 * t151 + t228) * t177, -t140 * t174 + t151 * t165 - t157 * t194 - t172 * t192 - t177 * t206, -t216 * t68 + t35, -t113 * t68 - t122 * t215 - t123 * t53 - t216 * t69, -t68 * t175, -t69 * t175, 0, -t113 * t116 + t122 * t98 + t124 * t69 + t127 * t53 - t286, t116 * t216 + t123 * t98 - t124 * t68 + t127 * t215 + t287, t11 * t122 - t113 * t20 + t12 * t123 + t21 * t216 + t215 * t76 - t217 * t53 - t44 * t68 + t48 * t69, t113 * t19 - t122 * t15 - t53 * t75 - t64 * t69 + t286, -t123 * t15 - t19 * t216 - t215 * t75 + t64 * t68 - t287, -t11 * t217 + t12 * t76 + t15 * t75 + t19 * t64 + t20 * t48 + t21 * t44, t69 * t284 + (t25 * t180 + t258 * t94) * t122 (-t180 * t92 + t183 * t94) * t69 + (-t180 * t26 + t24 + (-t183 * t92 - t284) * qJD(6)) * t122, t215 * t271 + t25 * t123 - t94 * t68 + (t122 * t258 + t180 * t69) * t305, t122 * t46 - t26 * t123 + t92 * t68 + (-t122 * t259 + t183 * t69) * t305, -t305 * t68 + t35, t5 * t123 + t224 * t68 + t50 * t26 + t9 * t92 + (-t7 * t123 - t305 * t8 - t290) * t180 + (t10 * t305 - t214) * t183 + ((-t180 * t49 - t183 * t43) * t305 - t18 * t123 + t31 * t271) * qJD(6), t18 * t68 + t50 * t25 + t9 * t94 + (-(qJD(6) * t49 + t8) * t305 - t290 - (qJD(6) * t27 + t7) * t123 + t31 * qJD(6) * t122) * t183 + (-(-qJD(6) * t43 + t10) * t305 - (-qJD(6) * t38 + t6) * t123 + t214) * t180; 0, 0, 0, -t182 * t269, t262 * t187, 0, 0, 0, t187 * pkin(1) * t182, pkin(1) * t269, t270, t95, t86, t87, 0, t139 * t173 - t227 * t177 + (-t177 * t254 - t211) * qJD(3) + t212, t264 * t177 + (t140 * t261 - t177 * t248) * pkin(2) + t199, -t317, t32, t314, t30, 0, t113 * t125 + t220 - t307, -t125 * t216 + t175 * t318 + t219, -t113 * t281 - t133 * t53 + t134 * t215 + t216 * t279 + t225, -t113 * t67 + t222 + t307, -t281 * t175 + t216 * t67 + t213, -t11 * t133 + t12 * t134 + t279 * t44 + t281 * t48 - t64 * t67, t16, t1, t13, t14, -t276, t130 * t46 + t133 * t26 + t282 * t92 + (t241 * t180 + t280 * t183) * t305 + t223, t133 * t25 + t282 * t94 + t241 * t315 + (-t130 * t215 - t280 * t305 - t243) * t180 + t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270, t95, t86, t87, 0, t212 + (-qJD(3) + t177) * t211, t177 * t226 + t199, -t317, t32, t314, t30, 0, t62 * t175 + (-t113 * t140 - t175 * t260) * pkin(3) + t220, t63 * t175 + (t140 * t216 - t175 * t247) * pkin(3) + t219, t113 * t277 - t169 * t53 + t170 * t215 + t216 * t230 + t225, -t113 * t70 + t175 * t230 + t222, t277 * t175 + t216 * t70 + t213, -t11 * t169 + t12 * t170 + t230 * t44 - t277 * t48 - t64 * t70, t16, t1, t13, t14, -t276, t168 * t46 + t169 * t26 + t278 * t92 + (t180 * t240 + t183 * t231) * t305 + t223, t169 * t25 + t278 * t94 + t240 * t315 + (-t168 * t215 - t231 * t305 - t243) * t180 + t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t317, t32, t314, t30, 0, t221 - t308, -t175 * t54 + t219, -pkin(4) * t215 - qJ(5) * t53 + (-t48 - t55) * t216 - (t44 - t265) * t113, -t113 * t78 - t221 + t288, t265 * t175 + t216 * t78 + t213, -t12 * pkin(4) - t11 * qJ(5) - t265 * t48 - t44 * t55 - t64 * t78, t16, t1, t13, t14, -t276, qJ(5) * t26 + t266 * t92 + (-t283 + t309) * t183 + (t180 * t239 - t183 * t34) * t305 + t313, qJ(5) * t25 + t266 * t94 + t239 * t315 + (t305 * t34 - t243 + t283) * t180 + t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, t317, -t175 ^ 2 - t303, t175 * t48 + t222, 0, 0, 0, 0, 0, -t175 * t92 + t218, -t175 * t94 + t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94 * t92, -t92 ^ 2 + t94 ^ 2, t25 + t321, -t323 * t94 + t47, t215, -t18 * t323 - t180 * t7 - t31 * t94 + t5, -t180 * t6 - t183 * t7 + t224 * t323 + t31 * t92;];
tauc_reg  = t2;
