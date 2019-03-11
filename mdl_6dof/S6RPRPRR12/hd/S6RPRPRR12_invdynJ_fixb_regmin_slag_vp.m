% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% tau_reg [6x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:20:10
% EndTime: 2019-03-09 04:20:15
% DurationCPUTime: 3.78s
% Computational Cost: add. (2826->428), mult. (5332->548), div. (0->0), fcn. (3395->10), ass. (0->237)
t174 = -pkin(1) - pkin(7);
t295 = pkin(4) - t174;
t165 = sin(qJ(6));
t166 = sin(qJ(5));
t169 = cos(qJ(6));
t170 = cos(qJ(5));
t104 = t165 * t170 + t166 * t169;
t171 = cos(qJ(3));
t257 = qJD(1) * t171;
t303 = qJD(5) + qJD(6);
t293 = (-t257 - t303) * t104;
t268 = t169 * t170;
t273 = t165 * t166;
t198 = -t268 + t273;
t147 = t171 * qJDD(1);
t167 = sin(qJ(3));
t244 = qJD(1) * qJD(3);
t227 = t167 * t244;
t308 = -t227 + t147;
t98 = -qJDD(5) - t308;
t90 = -qJDD(6) + t98;
t319 = t198 * t90;
t253 = qJD(3) * t170;
t258 = qJD(1) * t167;
t101 = t166 * t258 + t253;
t230 = t170 * t258;
t255 = qJD(3) * t166;
t99 = -t230 + t255;
t200 = t169 * t101 - t165 * t99;
t38 = t101 * t165 + t169 * t99;
t318 = t200 * t38;
t226 = t171 * t244;
t240 = t167 * qJDD(1);
t317 = t226 + t240;
t172 = cos(qJ(1));
t154 = g(2) * t172;
t168 = sin(qJ(1));
t155 = g(1) * t168;
t309 = t155 - t154;
t263 = pkin(3) * t258 + qJD(1) * qJ(2);
t81 = -qJ(4) * t257 + t263;
t316 = -t81 * qJD(1) - t309;
t315 = t200 ^ 2 - t38 ^ 2;
t136 = qJD(5) + t257;
t125 = qJD(6) + t136;
t247 = qJD(6) * t169;
t248 = qJD(6) * t165;
t251 = qJD(5) * t166;
t34 = -qJD(3) * t251 + qJD(5) * t230 + t170 * qJDD(3) + t166 * t317;
t35 = qJD(5) * t101 + t166 * qJDD(3) - t170 * t317;
t6 = -t101 * t248 - t165 * t35 + t169 * t34 - t99 * t247;
t314 = t125 * t38 + t6;
t173 = -pkin(3) - pkin(8);
t126 = t174 * qJD(1) + qJD(2);
t109 = t171 * t126;
t304 = qJD(4) - t109;
t264 = pkin(4) * t257 + t304;
t52 = t173 * qJD(3) + t264;
t280 = qJ(4) * t171;
t207 = pkin(8) * t167 - t280;
t54 = qJD(1) * t207 + t263;
t20 = t166 * t52 + t170 * t54;
t14 = -pkin(9) * t99 + t20;
t12 = t14 * t248;
t164 = qJ(5) + qJ(6);
t148 = sin(t164);
t297 = g(3) * t167;
t160 = qJD(3) * qJ(4);
t108 = t167 * t126;
t71 = -pkin(4) * t258 + t108;
t56 = t160 + t71;
t36 = pkin(5) * t99 + t56;
t149 = cos(t164);
t266 = t171 * t172;
t63 = -t148 * t266 - t149 * t168;
t269 = t168 * t171;
t65 = -t148 * t269 + t149 * t172;
t313 = g(1) * t65 - g(2) * t63 + t148 * t297 + t36 * t38 + t12;
t77 = t170 * t98;
t311 = -t136 * t251 - t77;
t277 = qJD(3) * t99;
t310 = t277 + t77;
t70 = t104 * t167;
t157 = qJDD(1) * qJ(2);
t159 = qJD(1) * qJD(2);
t238 = 0.2e1 * t159;
t307 = t238 + 0.2e1 * t157;
t250 = qJD(5) * t170;
t306 = qJD(6) * t170 + t250;
t152 = t167 * pkin(3);
t305 = -t280 + t152;
t302 = t136 * t56 - t173 * t98;
t201 = pkin(3) * t317 + qJ(4) * t227 + t157 + t159;
t219 = qJD(3) * pkin(8) - qJD(4);
t245 = qJ(4) * qJDD(1);
t21 = pkin(8) * t240 + (qJD(1) * t219 - t245) * t171 + t201;
t254 = qJD(3) * t167;
t102 = t126 * t254;
t123 = t174 * qJDD(1) + qJDD(2);
t202 = -t171 * t123 + qJDD(4) + t102;
t25 = pkin(4) * t308 + t173 * qJDD(3) + t202;
t223 = -t166 * t21 + t170 * t25;
t182 = -qJD(5) * t20 + t223;
t2 = -t98 * pkin(5) - t34 * pkin(9) + t182;
t239 = -t166 * t25 - t170 * t21 - t52 * t250;
t3 = -pkin(9) * t35 - t251 * t54 - t239;
t235 = -t165 * t3 + t169 * t2;
t19 = -t166 * t54 + t170 * t52;
t13 = -pkin(9) * t101 + t19;
t11 = pkin(5) * t136 + t13;
t287 = t14 * t169;
t5 = t11 * t165 + t287;
t62 = t148 * t168 - t149 * t266;
t64 = -t148 * t172 - t149 * t269;
t301 = -g(1) * t64 + g(2) * t62 - qJD(6) * t5 - t149 * t297 - t36 * t200 + t235;
t7 = qJD(6) * t200 + t165 * t34 + t169 * t35;
t300 = t125 * t200 - t7;
t220 = -qJD(3) * pkin(3) + qJD(4);
t78 = -t109 + t220;
t82 = -t108 - t160;
t204 = t167 * t78 - t171 * t82;
t107 = t167 * t123;
t156 = qJDD(3) * qJ(4);
t158 = qJD(3) * qJD(4);
t252 = qJD(3) * t171;
t42 = -t126 * t252 - t107 - t156 - t158;
t275 = qJDD(3) * pkin(3);
t43 = t202 - t275;
t180 = qJD(3) * t204 - t42 * t167 - t43 * t171;
t110 = qJ(2) + t305;
t241 = qJDD(3) * t174;
t299 = (qJD(1) * t110 + t81) * qJD(3) + t241;
t296 = g(3) * t171;
t294 = pkin(9) - t173;
t187 = -t165 * t251 - t166 * t248;
t234 = t170 * t257;
t292 = t169 * t234 - t257 * t273 + t268 * t303 + t187;
t106 = pkin(3) * t257 + qJ(4) * t258;
t68 = pkin(8) * t257 + t106;
t291 = t166 * t71 + t170 * t68;
t91 = qJ(2) + t152 + t207;
t114 = t295 * t171;
t92 = t166 * t114;
t290 = t170 * t91 + t92;
t289 = t104 * t90;
t288 = t136 * t99;
t286 = t166 * t98;
t284 = t34 * t170;
t236 = -pkin(5) * t170 - pkin(4);
t283 = pkin(5) * t250 - t236 * t257 + t304;
t282 = pkin(1) * qJDD(1);
t176 = qJD(1) ^ 2;
t281 = qJ(2) * t176;
t279 = qJD(3) * t38;
t278 = qJD(3) * t200;
t276 = qJD(5) * t54;
t274 = t101 * t136;
t272 = t166 * t167;
t271 = t166 * t171;
t270 = t167 * t170;
t267 = t170 * t171;
t262 = t172 * pkin(1) + t168 * qJ(2);
t162 = t167 ^ 2;
t163 = t171 ^ 2;
t260 = t162 - t163;
t175 = qJD(3) ^ 2;
t259 = t175 + t176;
t256 = qJD(3) * t101;
t249 = qJD(5) * t173;
t243 = qJDD(3) * t167;
t242 = qJDD(3) * t171;
t237 = t171 * t176 * t167;
t233 = t166 * t252;
t232 = t170 * t252;
t229 = pkin(3) * t252 + qJ(4) * t254 + qJD(2);
t228 = -pkin(9) * t167 - t91;
t113 = t294 * t170;
t225 = qJD(6) * t11 + t3;
t53 = t171 * t219 + t229;
t94 = t295 * t254;
t222 = -t166 * t53 - t170 * t94;
t221 = qJD(1) * t106 - g(3);
t216 = (-t162 - t163) * qJDD(1);
t215 = qJDD(2) - t282;
t214 = qJD(5) * t171 + qJD(1);
t213 = t125 * t293 + t319;
t95 = t295 * t252;
t111 = t294 * t166;
t196 = -pkin(5) * t167 - pkin(9) * t271;
t58 = t170 * t71;
t212 = qJD(1) * t196 - qJD(6) * t111 - t166 * t68 - t294 * t251 + t58;
t211 = pkin(9) * t234 + t113 * t303 + t291;
t210 = g(1) * t172 + g(2) * t168;
t208 = pkin(3) * t171 + qJ(4) * t167;
t206 = -t281 + t154;
t199 = t256 - t286;
t197 = t136 * t166;
t195 = -t136 * t250 + t286;
t193 = t114 * t250 - t166 * t94 + t170 * t53 - t251 * t91;
t192 = -t292 * t125 + t289;
t191 = t125 * t198;
t190 = 0.2e1 * qJ(2) * t244 + t241;
t189 = qJD(3) * t104;
t188 = -t174 * t175 - t210;
t186 = -t167 * t251 + t232;
t137 = g(1) * t269;
t185 = qJDD(4) + t137 + t81 * t257 + (-t123 - t154) * t171;
t184 = -t167 * t309 - t296;
t26 = -pkin(4) * t317 - t42;
t183 = t26 + t184;
t179 = t188 + t307;
t28 = (-qJD(1) * qJD(4) - t245) * t171 + t201;
t66 = -qJD(4) * t171 + t229;
t178 = -qJD(1) * t66 - qJDD(1) * t110 - t188 - t28;
t151 = t172 * qJ(2);
t142 = t167 * t174;
t138 = pkin(5) * t166 + qJ(4);
t112 = -pkin(4) * t167 + t142;
t93 = t170 * t114;
t89 = -t167 * t259 + t242;
t88 = t171 * t259 + t243;
t87 = -t166 * t269 + t170 * t172;
t86 = -t166 * t172 - t168 * t267;
t85 = -t166 * t266 - t168 * t170;
t84 = t166 * t168 - t170 * t266;
t73 = t167 * t236 + t142;
t69 = t165 * t272 - t167 * t268;
t46 = -pkin(5) * t186 - t95;
t32 = pkin(9) * t270 + t290;
t27 = t171 * pkin(5) + t166 * t228 + t93;
t16 = t165 * t233 - t169 * t232 + t303 * t70;
t15 = -t167 * t198 * t303 + t171 * t189;
t10 = pkin(5) * t35 + t26;
t9 = pkin(9) * t186 + t193;
t8 = t196 * qJD(3) + (t170 * t228 - t92) * qJD(5) + t222;
t4 = t11 * t169 - t14 * t165;
t1 = [qJDD(1), t309, t210, qJDD(2) - t309 - 0.2e1 * t282, -t210 + t307, -t215 * pkin(1) - g(1) * (-t168 * pkin(1) + t151) - g(2) * t262 + (t238 + t157) * qJ(2), qJDD(1) * t163 - 0.2e1 * t167 * t226, -0.2e1 * t147 * t167 + 0.2e1 * t244 * t260, -t167 * t175 + t242, -t171 * t175 - t243, 0, t167 * t179 + t171 * t190, -t167 * t190 + t171 * t179, t174 * t216 - t180 + t309, t178 * t167 - t171 * t299, t167 * t299 + t178 * t171, t28 * t110 + t81 * t66 - g(1) * (-qJ(4) * t266 + t172 * t152 + t151) - g(2) * (t172 * pkin(7) + t262) + (-g(1) * t174 - g(2) * t305) * t168 + t180 * t174, t34 * t272 + (t167 * t250 + t233) * t101 (t101 * t170 - t166 * t99) * t252 + (-t166 * t35 + t284 + (-t101 * t166 - t170 * t99) * qJD(5)) * t167 (t136 * t255 + t34) * t171 + (-t195 - t256) * t167 (t136 * t253 - t35) * t171 + (t277 + t311) * t167, -t136 * t254 - t171 * t98, t222 * t136 - (-t166 * t91 + t93) * t98 + t223 * t171 - t95 * t99 + t112 * t35 - t26 * t270 - g(1) * t85 - g(2) * t87 + (-t167 * t19 - t56 * t267) * qJD(3) + (-t290 * t136 - t20 * t171 + t56 * t272) * qJD(5), -t193 * t136 + t290 * t98 - t95 * t101 + t112 * t34 - g(1) * t84 - g(2) * t86 + ((qJD(3) * t56 + t276) * t166 + t239) * t171 + (qJD(3) * t20 + t26 * t166 + t250 * t56) * t167, t15 * t200 + t6 * t70, -t15 * t38 - t16 * t200 - t6 * t69 - t7 * t70, t125 * t15 + t171 * t6 - t200 * t254 - t70 * t90, -t125 * t16 - t171 * t7 + t254 * t38 + t69 * t90, -t125 * t254 - t171 * t90 (-t165 * t9 + t169 * t8) * t125 - (-t165 * t32 + t169 * t27) * t90 + t235 * t171 - t4 * t254 + t46 * t38 + t73 * t7 + t10 * t69 + t36 * t16 - g(1) * t63 - g(2) * t65 + ((-t165 * t27 - t169 * t32) * t125 - t5 * t171) * qJD(6), t5 * t254 - g(1) * t62 - g(2) * t64 + t10 * t70 + t12 * t171 + t36 * t15 + t46 * t200 + t73 * t6 + (-(-qJD(6) * t32 + t8) * t125 + t27 * t90 - t2 * t171) * t165 + (-(qJD(6) * t27 + t9) * t125 + t32 * t90 - t225 * t171) * t169; 0, 0, 0, qJDD(1), -t176, t215 - t309 - t281, 0, 0, 0, 0, 0, t89, -t88, t216, -t89, t88, t180 + t316, 0, 0, 0, 0, 0, t167 * t35 + t310 * t171 + (t166 * t214 + t167 * t253) * t136, t167 * t34 + t199 * t171 + (-t166 * t254 + t170 * t214) * t136, 0, 0, 0, 0, 0, t104 * t125 * qJD(1) + (-qJD(3) * t191 + t7) * t167 + ((t165 * t306 + t166 * t247 + t169 * t251) * t125 - t319 + t279) * t171, -qJD(1) * t191 + (-t125 * t189 + t6) * t167 + (-(-t169 * t306 - t187) * t125 - t289 + t278) * t171; 0, 0, 0, 0, 0, 0, t237, -t260 * t176, t147, -t240, qJDD(3), t297 - t137 + (t123 + t206) * t171, t296 - t107 + (-t206 + t155) * t167, -t208 * qJDD(1) + ((-t82 - t160) * t171 + (-t220 + t78) * t167) * qJD(1), t167 * t221 + t185 - 0.2e1 * t275, t167 * t316 + t221 * t171 + t107 + 0.2e1 * t156 + 0.2e1 * t158, -t43 * pkin(3) + g(3) * t305 - t42 * qJ(4) - t82 * qJD(4) - t81 * t106 - t126 * t204 - t208 * t309, -t101 * t197 + t284 (-t35 - t274) * t170 + (-t34 + t288) * t166 (t101 * t167 - t136 * t271) * qJD(1) + t311 (-t136 * t267 - t167 * t99) * qJD(1) + t195, t136 * t258, t19 * t258 + qJ(4) * t35 - t58 * t136 + t264 * t99 + t302 * t170 + ((t68 - t249) * t136 + t183) * t166, qJ(4) * t34 + t291 * t136 - t20 * t258 + t264 * t101 - t302 * t166 + (-t136 * t249 + t183) * t170, -t198 * t6 + t200 * t293, -t6 * t104 + t198 * t7 - t200 * t292 - t293 * t38, t200 * t258 + t213, -t258 * t38 + t192, t125 * t258 -(t111 * t165 - t113 * t169) * t90 + t138 * t7 + t10 * t104 + t4 * t258 + t283 * t38 + t292 * t36 + (t165 * t211 - t169 * t212) * t125 + t184 * t148 (-t111 * t169 - t113 * t165) * t90 + t138 * t6 - t10 * t198 - t5 * t258 + t283 * t200 + t293 * t36 + (t165 * t212 + t169 * t211) * t125 + t184 * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, qJDD(3) - t237, -t163 * t176 - t175, qJD(3) * t82 + t102 + t185 - t275 - t297, 0, 0, 0, 0, 0, -t136 * t197 - t310, -t136 ^ 2 * t170 - t199, 0, 0, 0, 0, 0, t213 - t279, t192 - t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101 * t99, t101 ^ 2 - t99 ^ 2, t34 + t288, t274 - t35, -t98, -g(1) * t86 + g(2) * t84 - g(3) * t270 - t101 * t56 + t136 * t20 + t182, g(1) * t87 - g(2) * t85 + t136 * t19 + t56 * t99 + (t276 + t297) * t166 + t239, t318, t315, t314, t300, -t90 -(-t13 * t165 - t287) * t125 + (-t101 * t38 - t125 * t248 - t169 * t90) * pkin(5) + t301 (-t125 * t14 - t2) * t165 + (t125 * t13 - t225) * t169 + (-t101 * t200 - t125 * t247 + t165 * t90) * pkin(5) + t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t318, t315, t314, t300, -t90, t125 * t5 + t301, t125 * t4 - t165 * t2 - t169 * t225 + t313;];
tau_reg  = t1;
