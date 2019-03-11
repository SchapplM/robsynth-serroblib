% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:14:04
% EndTime: 2019-03-08 23:14:14
% DurationCPUTime: 4.04s
% Computational Cost: add. (4249->432), mult. (9590->565), div. (0->0), fcn. (7602->14), ass. (0->236)
t170 = sin(qJ(4));
t171 = sin(qJ(3));
t174 = cos(qJ(3));
t315 = cos(qJ(4));
t122 = t170 * t174 + t315 * t171;
t115 = t122 * qJD(2);
t329 = qJD(6) + t115;
t258 = t315 * t174;
t230 = qJD(2) * t258;
t273 = qJD(2) * t171;
t253 = t170 * t273;
t113 = -t230 + t253;
t161 = qJD(3) + qJD(4);
t169 = sin(qJ(6));
t173 = cos(qJ(6));
t91 = -t173 * t113 + t161 * t169;
t331 = t329 * t91;
t324 = t329 - qJD(6);
t176 = -pkin(9) - pkin(8);
t160 = qJDD(3) + qJDD(4);
t165 = qJ(3) + qJ(4);
t158 = sin(t165);
t168 = cos(pkin(6));
t175 = cos(qJ(2));
t292 = cos(pkin(11));
t242 = t292 * t175;
t166 = sin(pkin(11));
t172 = sin(qJ(2));
t286 = t166 * t172;
t107 = -t168 * t242 + t286;
t243 = t292 * t172;
t285 = t166 * t175;
t109 = t168 * t285 + t243;
t227 = g(1) * t109 + g(2) * t107;
t167 = sin(pkin(6));
t282 = t167 * t175;
t197 = g(3) * t282 - t227;
t195 = t197 * t158;
t259 = qJD(3) * t176;
t123 = t171 * t259;
t124 = t174 * t259;
t131 = t176 * t171;
t132 = t176 * t174;
t250 = qJD(4) * t315;
t275 = qJD(1) * t167;
t257 = t175 * t275;
t271 = qJD(4) * t170;
t280 = t170 * t171;
t326 = t258 - t280;
t305 = -t315 * t123 - t170 * t124 - t131 * t250 - t132 * t271 + t326 * t257;
t95 = t170 * t131 - t315 * t132;
t330 = -t160 * t95 + t305 * t161 + t195;
t304 = t95 * qJD(4) - t122 * t257 + t170 * t123 - t315 * t124;
t252 = t172 * t275;
t239 = -t176 * qJD(2) + t252;
t274 = qJD(1) * t168;
t88 = t171 * t274 + t174 * t239;
t82 = t170 * t88;
t87 = -t239 * t171 + t174 * t274;
t46 = t315 * t87 - t82;
t293 = pkin(3) * t250 + qJD(5) - t46;
t83 = t315 * t88;
t45 = t170 * t87 + t83;
t228 = pkin(3) * t271 - t45;
t284 = t167 * t172;
t179 = qJD(2) ^ 2;
t201 = (qJDD(2) * t175 - t172 * t179) * t167;
t302 = qJD(3) * pkin(3);
t84 = t87 + t302;
t43 = -t315 * t84 + t82;
t278 = qJD(5) + t43;
t281 = t169 * t175;
t134 = t167 * t281;
t283 = t167 * t174;
t112 = t168 * t171 + t172 * t283;
t204 = -t168 * t174 + t171 * t284;
t63 = t170 * t112 + t204 * t315;
t328 = t173 * t63 + t134;
t154 = t160 * qJ(5);
t327 = -t161 * qJD(5) - t154;
t325 = qJD(3) * t204;
t152 = -t315 * pkin(3) - pkin(4);
t149 = -pkin(10) + t152;
t310 = t113 * pkin(5);
t222 = t161 * t280;
t245 = qJDD(2) * t315;
t263 = t174 * qJDD(2);
t233 = -t161 * t230 - t170 * t263 - t171 * t245;
t56 = qJD(2) * t222 + t233;
t53 = -qJDD(6) + t56;
t323 = t329 * (t310 + t228) - t149 * t53;
t178 = qJD(3) ^ 2;
t267 = qJD(1) * qJD(2);
t249 = t172 * t267;
t221 = -qJDD(1) * t282 + t167 * t249;
t311 = g(3) * t175;
t321 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t178 + t167 * (t249 - t311) - t221 + t227;
t254 = qJD(2) * t282;
t231 = t174 * t254;
t183 = t231 - t325;
t198 = t112 * qJD(3);
t232 = t171 * t254;
t184 = t198 + t232;
t18 = t112 * t271 + t170 * t184 - t315 * t183 + t204 * t250;
t272 = qJD(2) * t172;
t64 = t315 * t112 - t170 * t204;
t320 = t167 * (-t115 * t272 - t175 * t56) + t160 * t64 - t161 * t18;
t19 = t64 * qJD(4) + t170 * t183 + t315 * t184;
t264 = t171 * qJDD(2);
t220 = t170 * t264 - t174 * t245;
t86 = t161 * t122;
t57 = qJD(2) * t86 + t220;
t319 = t167 * (t113 * t272 - t175 * t57) - t160 * t63 - t161 * t19;
t159 = cos(t165);
t94 = -t315 * t131 - t170 * t132;
t318 = -t197 * t159 - t160 * t94 - t304 * t161;
t317 = t115 ^ 2;
t316 = pkin(4) + pkin(10);
t314 = pkin(4) * t160;
t309 = t115 * pkin(5);
t153 = pkin(3) * t174 + pkin(2);
t212 = -qJ(5) * t122 - t153;
t58 = -t316 * t326 + t212;
t308 = t58 * t53;
t307 = -pkin(5) * t86 - t305;
t85 = -qJD(3) * t258 - t174 * t250 + t222;
t306 = -t85 * pkin(5) + t304;
t44 = t170 * t84 + t83;
t303 = qJD(2) * pkin(2);
t298 = t161 * t44;
t93 = t113 * t169 + t161 * t173;
t297 = t169 * t93;
t269 = qJD(6) * t173;
t270 = qJD(6) * t169;
t24 = t113 * t269 + t173 * t160 - t161 * t270 + t169 * t57;
t296 = t173 * t24;
t50 = t173 * t53;
t294 = t309 + t293;
t290 = t329 * t113;
t289 = t115 * t113;
t288 = t326 * t169;
t287 = t166 * t167;
t279 = t309 + t278;
t277 = qJDD(1) - g(3);
t163 = t171 ^ 2;
t276 = -t174 ^ 2 + t163;
t266 = qJD(2) * qJD(3);
t265 = t168 * qJDD(1);
t262 = g(3) * t284;
t156 = t171 * t302;
t260 = t173 * t282;
t255 = t167 * t272;
t20 = -t316 * t161 + t279;
t104 = -qJD(2) * t153 - t257;
t190 = -qJ(5) * t115 + t104;
t42 = t316 * t113 + t190;
t11 = t169 * t20 + t173 * t42;
t138 = t174 * t265;
t98 = qJDD(2) * pkin(8) + (qJDD(1) * t172 + t175 * t267) * t167;
t238 = pkin(9) * qJDD(2) + t98;
t37 = qJDD(3) * pkin(3) - t88 * qJD(3) - t238 * t171 + t138;
t39 = t87 * qJD(3) + t171 * t265 + t238 * t174;
t240 = -t170 * t37 - t84 * t250 + t88 * t271 - t315 * t39;
t7 = t240 + t327;
t6 = -pkin(5) * t57 - t7;
t251 = -t11 * t113 + t6 * t173;
t248 = t171 * t266;
t108 = t168 * t243 + t285;
t244 = t167 * t292;
t74 = t108 * t158 + t159 * t244;
t75 = t108 * t159 - t158 * t244;
t247 = -t74 * pkin(4) + t75 * qJ(5);
t110 = -t168 * t286 + t242;
t76 = t110 * t158 - t159 * t287;
t77 = t110 * t159 + t158 * t287;
t246 = -t76 * pkin(4) + qJ(5) * t77;
t241 = t170 * t39 + t88 * t250 + t84 * t271 - t315 * t37;
t101 = t158 * t284 - t168 * t159;
t102 = t158 * t168 + t159 * t284;
t237 = -t101 * pkin(4) + qJ(5) * t102;
t236 = t160 * t169 - t173 * t57;
t38 = -qJ(5) * t161 - t44;
t26 = -t38 - t310;
t235 = t329 * t26;
t234 = t169 * t329;
t226 = g(1) * t110 + g(2) * t108;
t203 = qJ(5) * t85 - qJD(5) * t122 + t156;
t225 = -t316 * t86 - t203 + t252;
t27 = pkin(4) * t86 + t203;
t224 = -t27 + t252;
t223 = qJDD(5) + t241;
t67 = pkin(4) * t115 + qJ(5) * t113;
t219 = -(t44 - t310) * t329 + t316 * t53;
t216 = t169 * t42 - t173 * t20;
t215 = -t216 * t113 + t6 * t169 + (t115 * t173 + t269) * t26;
t60 = pkin(3) * t273 + t67;
t210 = -g(1) * t166 + t292 * g(2);
t209 = -t169 * t63 + t260;
t208 = -t234 * t329 - t50;
t202 = -g(1) * t77 - g(2) * t75 - g(3) * t102;
t200 = -t252 + t156;
t199 = -t173 * t329 ^ 2 + t169 * t53;
t126 = -t257 - t303;
t196 = -qJD(2) * t126 + t226 - t98;
t194 = g(1) * t76 + g(2) * t74 + g(3) * t101 - t241;
t193 = -t202 + t240;
t73 = pkin(3) * t248 - qJDD(2) * t153 + t221;
t61 = t122 * pkin(5) + t94;
t191 = t26 * t86 - t326 * t6 + t61 * t53 + t226;
t189 = -pkin(8) * qJDD(3) + (t126 + t257 - t303) * qJD(3);
t105 = t115 * pkin(10);
t188 = (-qJD(6) * t149 + t105 + t60) * t329 + t202;
t187 = (qJD(6) * t316 + t105 + t67) * t329 + t202;
t186 = -t104 * t115 + t194;
t185 = t104 * t113 + t193;
t40 = -t233 + (t113 - t253) * t161;
t55 = pkin(4) * t113 + t190;
t182 = t115 * t55 + qJDD(5) - t194;
t181 = -t113 * t55 - t193 - t327;
t180 = qJ(5) * t56 - qJD(5) * t115 + t73;
t150 = pkin(3) * t170 + qJ(5);
t72 = -pkin(4) * t326 + t212;
t62 = pkin(5) * t326 + t95;
t59 = -t113 ^ 2 + t317;
t34 = -pkin(4) * t161 + t278;
t25 = qJD(6) * t93 + t236;
t15 = -t113 * t91 + t199;
t14 = t113 * t93 + t208;
t13 = pkin(4) * t57 + t180;
t12 = -t234 * t93 + t296;
t9 = t316 * t57 + t180;
t8 = t223 - t314;
t5 = -pkin(5) * t56 - t316 * t160 + t223;
t2 = t173 * t5;
t1 = (-t329 * t93 - t25) * t173 + (-t24 + t331) * t169;
t3 = [t277, 0, t201 (-qJDD(2) * t172 - t175 * t179) * t167, 0, 0, 0, 0, 0, -t204 * qJDD(3) + t174 * t201 + (-t198 - 0.2e1 * t232) * qJD(3), -t112 * qJDD(3) - t171 * t201 + (-0.2e1 * t231 + t325) * qJD(3), 0, 0, 0, 0, 0, t319, -t320, t113 * t18 + t115 * t19 - t56 * t63 - t57 * t64, -t319, t320, t18 * t38 + t19 * t34 + t63 * t8 - t64 * t7 - g(3) + (-t13 * t175 + t272 * t55) * t167, 0, 0, 0, 0, 0 (qJD(6) * t209 - t169 * t255 + t173 * t19) * t329 - t328 * t53 - t18 * t91 + t64 * t25 -(t328 * qJD(6) + t169 * t19 + t173 * t255) * t329 - t209 * t53 - t18 * t93 + t64 * t24; 0, qJDD(2), t277 * t282 + t227, -t277 * t284 + t226, qJDD(2) * t163 + 0.2e1 * t174 * t248, 0.2e1 * t171 * t263 - 0.2e1 * t276 * t266, qJDD(3) * t171 + t174 * t178, qJDD(3) * t174 - t171 * t178, 0, t189 * t171 + t321 * t174, -t321 * t171 + t189 * t174, -t115 * t85 - t122 * t56, t113 * t85 - t115 * t86 - t122 * t57 - t326 * t56, t122 * t160 - t161 * t85, t160 * t326 - t161 * t86, 0, t104 * t86 + t200 * t113 - t153 * t57 - t326 * t73 + t318, -t104 * t85 + t200 * t115 + t122 * t73 + t153 * t56 + t330, t305 * t113 + t304 * t115 + t122 * t8 - t326 * t7 - t34 * t85 + t38 * t86 - t56 * t94 - t57 * t95 - t226 - t262, t224 * t113 + t13 * t326 - t55 * t86 - t57 * t72 - t318, t224 * t115 - t122 * t13 + t55 * t85 + t56 * t72 - t330, t13 * t72 + t55 * t27 - t7 * t95 + t8 * t94 + t305 * t38 + t304 * t34 + t226 * t176 + (g(3) * t176 - qJD(1) * t55) * t284 + (-t311 * t167 + t227) * (pkin(4) * t159 + qJ(5) * t158 + t153) t86 * t297 - (t169 * t24 + t269 * t93) * t326 (-t169 * t91 + t173 * t93) * t86 - (-t169 * t25 + t296 + (-t173 * t91 - t297) * qJD(6)) * t326, t53 * t288 + t122 * t24 - t85 * t93 + (t169 * t86 - t269 * t326) * t329, t326 * t50 - t122 * t25 + t85 * t91 + (t173 * t86 + t270 * t326) * t329, -t122 * t53 - t329 * t85, t216 * t85 + t2 * t122 + t62 * t25 + t307 * t91 + (-t9 * t122 + t158 * t227 + t308) * t169 - t191 * t173 - g(3) * (t158 * t281 + t172 * t173) * t167 + (t225 * t169 + t306 * t173) * t329 + ((-t169 * t61 - t173 * t58) * t329 - t11 * t122 - t26 * t288) * qJD(6), t11 * t85 + t62 * t24 + t307 * t93 + (t308 - (qJD(6) * t20 + t9) * t122 - t26 * qJD(6) * t326 + (-qJD(6) * t61 + t225) * t329 - t195) * t173 + (-(-qJD(6) * t42 + t5) * t122 + t262 + (qJD(6) * t58 - t306) * t329 + t191) * t169; 0, 0, 0, 0, -t171 * t179 * t174, t276 * t179, t264, t263, qJDD(3), g(3) * t204 + t171 * t196 + t210 * t283 + t138, g(3) * t112 + (-t167 * t210 - t265) * t171 + t196 * t174, t289, t59, t40, -t220, t160, t45 * t161 + (-t113 * t273 + t315 * t160 - t161 * t271) * pkin(3) + t186, t46 * t161 + (-t115 * t273 - t160 * t170 - t161 * t250) * pkin(3) + t185, -t150 * t57 - t152 * t56 + (t228 - t38) * t115 + (t34 - t293) * t113, t113 * t60 + t228 * t161 + (-pkin(4) + t152) * t160 + t182, t115 * t60 + t150 * t160 + t293 * t161 + t181, -t7 * t150 + t8 * t152 - t55 * t60 - g(1) * ((-t110 * t171 + t166 * t283) * pkin(3) + t246) - g(2) * ((-t108 * t171 - t174 * t244) * pkin(3) + t247) - g(3) * (-pkin(3) * t204 + t237) - t293 * t38 + t228 * t34, t12, t1, t14, t15, t290, t150 * t25 + t188 * t169 + t323 * t173 + t294 * t91 + t215, t150 * t24 + t294 * t93 + t188 * t173 + (-t235 - t323) * t169 + t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, t59, t40, -t220, t160, t186 + t298, -t161 * t43 + t185, pkin(4) * t56 - qJ(5) * t57 + (-t38 - t44) * t115 + (t34 - t278) * t113, t113 * t67 + t182 - t298 - 0.2e1 * t314, t115 * t67 + t278 * t161 + t154 + t181, -t8 * pkin(4) - g(1) * t246 - g(2) * t247 - g(3) * t237 - t7 * qJ(5) - t278 * t38 - t34 * t44 - t55 * t67, t12, t1, t14, t15, t290, qJ(5) * t25 + t169 * t187 + t173 * t219 + t279 * t91 + t215, qJ(5) * t24 + t279 * t93 + (-t235 - t219) * t169 + t187 * t173 + t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t160 - t289, -t161 ^ 2 - t317, t161 * t38 + t182 - t314, 0, 0, 0, 0, 0, -t161 * t91 + t208, -t161 * t93 + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93 * t91, -t91 ^ 2 + t93 ^ 2, t24 + t331, t324 * t93 - t236, -t53, -t169 * t9 + t2 - t26 * t93 - g(1) * (-t109 * t169 + t173 * t76) - g(2) * (-t107 * t169 + t173 * t74) - g(3) * (t101 * t173 + t134) + t324 * t11, -t173 * t9 - t169 * t5 + t26 * t91 - g(1) * (-t109 * t173 - t169 * t76) - g(2) * (-t107 * t173 - t169 * t74) - g(3) * (-t101 * t169 + t260) - t324 * t216;];
tau_reg  = t3;
