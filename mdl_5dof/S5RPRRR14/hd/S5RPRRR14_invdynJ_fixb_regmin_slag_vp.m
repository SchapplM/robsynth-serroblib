% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR14
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR14_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:45
% EndTime: 2019-12-31 19:20:02
% DurationCPUTime: 6.62s
% Computational Cost: add. (8152->478), mult. (26286->696), div. (0->0), fcn. (23167->14), ass. (0->234)
t203 = sin(pkin(5));
t202 = sin(pkin(6));
t314 = cos(pkin(5));
t271 = t314 * t202;
t327 = cos(qJ(3));
t240 = t327 * t271;
t205 = cos(pkin(6));
t204 = cos(pkin(11));
t285 = t327 * t204;
t262 = t205 * t285;
t342 = t203 * t262 + t240;
t206 = sin(qJ(5));
t210 = cos(qJ(5));
t201 = sin(pkin(11));
t209 = sin(qJ(1));
t212 = cos(qJ(1));
t273 = t212 * t314;
t161 = t209 * t201 - t204 * t273;
t307 = t203 * t212;
t136 = -t161 * t202 + t205 * t307;
t207 = sin(qJ(4));
t211 = cos(qJ(4));
t162 = t201 * t273 + t209 * t204;
t208 = sin(qJ(3));
t306 = t205 * t208;
t310 = t202 * t208;
t90 = -t161 * t306 + t162 * t327 - t307 * t310;
t65 = t136 * t207 - t211 * t90;
t288 = t202 * t327;
t264 = t203 * t288;
t287 = t205 * t327;
t89 = t161 * t287 + t162 * t208 + t212 * t264;
t341 = t206 * t65 + t210 * t89;
t340 = -t206 * t89 + t210 * t65;
t222 = t203 * (-t201 * t306 + t285);
t149 = qJD(1) * t222;
t281 = qJD(3) * t327;
t257 = t202 * t281;
t339 = t149 - t257;
t255 = t208 * t271;
t286 = t327 * t201;
t134 = (t204 * t306 + t286) * t203 + t255;
t127 = t134 * qJD(1);
t268 = qJD(1) * t314;
t304 = qJD(1) * t203;
t283 = t204 * t304;
t296 = -t202 * t283 + qJD(3);
t227 = -t205 * t268 - t296;
t152 = t211 * t227;
t83 = t127 * t207 + t152;
t82 = qJD(5) + t83;
t284 = t201 * t304;
t335 = t342 * qJD(1) - t208 * t284;
t117 = qJD(4) - t335;
t334 = t136 * t211 + t207 * t90;
t199 = t203 ^ 2;
t333 = t199 * (t201 ^ 2 + t204 ^ 2);
t166 = -t205 * t211 + t207 * t310;
t261 = t202 * t284;
t332 = t166 * qJD(4) + t207 * t261 + t339 * t211;
t289 = pkin(1) * t314;
t309 = t203 * t204;
t165 = qJ(2) * t309 + t201 * t289;
t220 = (t205 * t309 + t271) * pkin(8);
t130 = t220 + t165;
t195 = t204 * t289;
t312 = t201 * t203;
t216 = t314 * pkin(2) + (-pkin(8) * t205 - qJ(2)) * t312;
t135 = t195 + t216;
t153 = (-pkin(8) * t201 * t202 - pkin(2) * t204 - pkin(1)) * t203;
t215 = -t208 * t130 + t135 * t287 + t153 * t288;
t167 = t205 * t207 + t211 * t310;
t331 = t167 * qJD(4) - t339 * t207 + t211 * t261;
t221 = t203 * (t204 * t208 + t205 * t286);
t148 = qJD(1) * t221;
t303 = qJD(3) * t208;
t282 = t202 * t303;
t330 = t148 - t282;
t126 = t134 * qJD(3);
t274 = t209 * t314;
t225 = t212 * t201 + t204 * t274;
t308 = t203 * t209;
t138 = t202 * t225 + t205 * t308;
t163 = -t201 * t274 + t212 * t204;
t219 = t225 * t205;
t94 = t163 * t327 + (t202 * t308 - t219) * t208;
t66 = t138 * t211 - t207 * t94;
t270 = t314 * t205;
t160 = t202 * t309 - t270;
t87 = t134 * t207 + t160 * t211;
t235 = g(1) * t66 - g(2) * t334 - g(3) * t87;
t258 = pkin(1) * t268;
t190 = t204 * t258;
t118 = qJD(1) * t216 + t190;
t143 = qJD(1) * t153 + qJD(2);
t81 = -t118 * t202 + t205 * t143;
t40 = -pkin(3) * t335 - pkin(9) * t127 + t81;
t104 = t118 * t306;
t159 = qJ(2) * t283 + t201 * t258;
t111 = qJD(1) * t220 + t159;
t55 = t327 * t111 + t143 * t310 + t104;
t42 = -pkin(9) * t227 + t55;
t17 = t207 * t40 + t211 * t42;
t266 = qJDD(1) * t314;
t293 = qJDD(1) * t204;
t278 = t203 * t293;
t156 = t202 * t278 - t205 * t266 - qJDD(3);
t256 = pkin(1) * t266;
t188 = t204 * t256;
t295 = qJD(1) * qJD(2);
t280 = t203 * t295;
t100 = qJDD(1) * t216 - t201 * t280 + t188;
t139 = qJDD(1) * t153 + qJDD(2);
t263 = t118 * t287;
t145 = qJ(2) * t278 + t201 * t256 + t204 * t280;
t99 = qJDD(1) * t220 + t145;
t223 = -qJD(3) * t263 - t100 * t306 + t111 * t303 - t139 * t310 - t143 * t257 - t327 * t99;
t21 = -pkin(9) * t156 - t223;
t70 = -t100 * t202 + t205 * t139;
t294 = qJDD(1) * t203;
t279 = t201 * t294;
t78 = t335 * qJD(3) + qJDD(1) * t255 + t278 * t306 + t327 * t279;
t232 = -qJDD(1) * t342 + t208 * t279;
t79 = qJD(1) * t126 + t232;
t29 = pkin(3) * t79 - pkin(9) * t78 + t70;
t275 = t207 * t21 - t211 * t29;
t77 = qJDD(4) + t79;
t4 = -pkin(4) * t77 + qJD(4) * t17 + t275;
t85 = t211 * t127 - t207 * t227;
t329 = (pkin(4) * t85 + t82 * pkin(10)) * t82 + t235 + t4;
t182 = -pkin(4) * t211 - pkin(10) * t207 - pkin(3);
t36 = qJD(4) * t85 + t211 * t156 + t207 * t78;
t34 = qJDD(5) + t36;
t328 = (-t55 + t117 * (pkin(4) * t207 - pkin(10) * t211)) * t82 + t182 * t34;
t301 = qJD(4) * t207;
t35 = -qJD(4) * t152 - t127 * t301 - t207 * t156 + t211 * t78;
t58 = t117 * t206 + t210 * t85;
t12 = qJD(5) * t58 + t206 * t35 - t210 * t77;
t14 = pkin(10) * t117 + t17;
t54 = -t208 * t111 + t143 * t288 + t263;
t41 = pkin(3) * t227 - t54;
t23 = t83 * pkin(4) - t85 * pkin(10) + t41;
t249 = t14 * t206 - t210 * t23;
t299 = qJD(4) * t211;
t238 = t207 * t29 + t211 * t21 + t40 * t299 - t42 * t301;
t3 = pkin(10) * t77 + t238;
t239 = -qJD(3) * t104 + t100 * t287 - t111 * t281 + t139 * t288 - t143 * t282 - t208 * t99;
t22 = pkin(3) * t156 - t239;
t6 = pkin(4) * t36 - pkin(10) * t35 + t22;
t1 = -t249 * qJD(5) + t206 * t6 + t210 * t3;
t326 = pkin(1) * t199;
t56 = -t210 * t117 + t206 * t85;
t325 = t56 * t82;
t324 = t58 * t82;
t311 = t201 * t208;
t133 = t203 * t311 - t342;
t86 = -t135 * t202 + t205 * t153;
t49 = pkin(3) * t133 - pkin(9) * t134 + t86;
t116 = t327 * t130;
t290 = t135 * t306 + t153 * t310 + t116;
t53 = -pkin(9) * t160 + t290;
t244 = t207 * t49 + t211 * t53;
t80 = pkin(3) * t127 - pkin(9) * t335;
t323 = t207 * t80 + t211 * t54;
t322 = pkin(9) * qJD(4);
t297 = qJD(5) * t210;
t298 = qJD(5) * t206;
t11 = t117 * t297 + t206 * t77 + t210 * t35 - t85 * t298;
t321 = t11 * t206;
t320 = t117 * t83;
t318 = t206 * t34;
t317 = t210 * t34;
t316 = t85 * t117;
t313 = t335 * t211;
t302 = qJD(4) * t206;
t300 = qJD(4) * t210;
t292 = t82 * t302;
t291 = t82 * t300;
t213 = qJD(1) ^ 2;
t272 = t213 * t314;
t267 = t210 * t82;
t265 = t117 * t211;
t259 = qJD(2) * t202 * t312;
t253 = g(1) * t212 + g(2) * t209;
t252 = g(1) * t209 - g(2) * t212;
t250 = qJD(2) * t268;
t8 = t14 * t210 + t206 * t23;
t20 = pkin(10) * t133 + t244;
t52 = t160 * pkin(3) - t215;
t88 = t134 * t211 - t160 * t207;
t28 = t87 * pkin(4) - t88 * pkin(10) + t52;
t248 = t20 * t210 + t206 * t28;
t247 = -t20 * t206 + t210 * t28;
t16 = -t207 * t42 + t211 * t40;
t44 = qJD(2) * t222 + t215 * qJD(3);
t125 = (t240 + (t262 - t311) * t203) * qJD(3);
t69 = pkin(3) * t126 - pkin(9) * t125 + t259;
t246 = -t207 * t44 + t211 * t69;
t245 = -t207 * t53 + t211 * t49;
t243 = t133 * t210 - t206 * t88;
t62 = t133 * t206 + t210 * t88;
t241 = (-qJ(2) * t284 + t190) * t201 - t159 * t204;
t237 = t207 * t69 + t211 * t44 + t49 * t299 - t53 * t301;
t236 = -pkin(9) * t77 + t117 * t41;
t93 = t163 * t208 - t209 * t264 + t327 * t219;
t234 = g(1) * t93 + g(2) * t89 + g(3) * t133;
t233 = -g(1) * t94 - g(2) * t90 - g(3) * t134;
t230 = -t206 * t167 - t210 * t288;
t228 = -t210 * t167 + t206 * t288;
t226 = -t22 + t234;
t13 = -pkin(4) * t117 - t16;
t224 = -pkin(10) * t34 + (t13 + t16) * t82;
t217 = pkin(9) * qJD(5) * t82 - t234;
t2 = -t8 * qJD(5) - t206 * t3 + t210 * t6;
t214 = (pkin(10) * t127 - qJD(5) * t182 + t323) * t82 + t233;
t45 = qJD(2) * t221 + (t116 + (t135 * t205 + t153 * t202) * t208) * qJD(3);
t191 = -pkin(1) * t294 + qJDD(2);
t164 = -qJ(2) * t312 + t195;
t144 = t188 + (-qJ(2) * qJDD(1) - t295) * t312;
t72 = t127 * t206 + t210 * t313;
t71 = -t210 * t127 + t206 * t313;
t67 = t138 * t207 + t211 * t94;
t60 = -qJD(4) * t87 + t125 * t211;
t59 = qJD(4) * t88 + t125 * t207;
t38 = t206 * t93 + t210 * t67;
t37 = -t206 * t67 + t210 * t93;
t31 = qJD(5) * t243 + t126 * t206 + t210 * t60;
t30 = qJD(5) * t62 - t126 * t210 + t206 * t60;
t24 = -pkin(4) * t127 + t207 * t54 - t211 * t80;
t19 = -pkin(4) * t133 - t245;
t15 = t59 * pkin(4) - t60 * pkin(10) + t45;
t10 = -pkin(4) * t126 + qJD(4) * t244 - t246;
t9 = pkin(10) * t126 + t237;
t5 = [qJDD(1), t252, t253, t144 * t314 + g(1) * t162 - g(2) * t163 + (-t191 * t204 - t201 * t250) * t203 + (t314 * t164 + t204 * t326) * qJDD(1), -t145 * t314 - g(1) * t161 + g(2) * t225 + (t191 * t201 - t204 * t250) * t203 + (-t165 * t314 - t201 * t326) * qJDD(1), t295 * t333 + (-t144 * t201 + t145 * t204 + (-t164 * t201 + t165 * t204) * qJDD(1) - t253) * t203, t144 * t164 + t145 * t165 + t252 * pkin(1) + (-t191 * pkin(1) - qJ(2) * t253 - qJD(2) * t241) * t203, t125 * t127 + t134 * t78, t125 * t335 - t126 * t127 - t133 * t78 - t134 * t79, -t125 * t227 - t134 * t156 - t78 * t160, t126 * t227 + t133 * t156 + t79 * t160, t156 * t160, g(1) * t90 - g(2) * t94 + t81 * t126 + t70 * t133 - t156 * t215 - t160 * t239 + t227 * t45 - t259 * t335 + t86 * t79, -g(1) * t89 + g(2) * t93 + t81 * t125 + t127 * t259 + t70 * t134 + t156 * t290 - t160 * t223 + t227 * t44 + t86 * t78, t35 * t88 + t60 * t85, -t35 * t87 - t36 * t88 - t59 * t85 - t60 * t83, t117 * t60 + t126 * t85 + t133 * t35 + t77 * t88, -t117 * t59 - t126 * t83 - t133 * t36 - t77 * t87, t117 * t126 + t133 * t77, t246 * t117 + t245 * t77 - t275 * t133 + t16 * t126 + t45 * t83 + t52 * t36 + t22 * t87 + t41 * t59 - g(1) * t65 - g(2) * t67 + (-t117 * t244 - t133 * t17) * qJD(4), -g(1) * t334 - g(2) * t66 - t237 * t117 - t17 * t126 - t238 * t133 + t22 * t88 - t244 * t77 + t52 * t35 + t41 * t60 + t45 * t85, t11 * t62 + t31 * t58, t11 * t243 - t12 * t62 - t30 * t58 - t31 * t56, t11 * t87 + t31 * t82 + t34 * t62 + t58 * t59, -t12 * t87 + t243 * t34 - t30 * t82 - t56 * t59, t34 * t87 + t59 * t82, (-qJD(5) * t248 + t15 * t210 - t206 * t9) * t82 + t247 * t34 + t2 * t87 - t249 * t59 + t10 * t56 + t19 * t12 - t4 * t243 + t13 * t30 - g(1) * t340 - g(2) * t38, -(qJD(5) * t247 + t15 * t206 + t210 * t9) * t82 - t248 * t34 - t1 * t87 - t8 * t59 + t10 * t58 + t19 * t11 + t4 * t62 + t13 * t31 + g(1) * t341 - g(2) * t37; 0, 0, 0, (t201 * t272 - t293) * t203, (qJDD(1) * t201 + t204 * t272) * t203, -t213 * t333, -g(3) * t314 + qJDD(2) + (-pkin(1) * qJDD(1) + qJD(1) * t241 - t252) * t203, 0, 0, 0, 0, 0, t205 * t79 - t148 * t227 + (-t327 * t156 + t227 * t303 + t284 * t335) * t202, t205 * t78 - t149 * t227 + (-t127 * t284 + t208 * t156 + t227 * t281) * t202, 0, 0, 0, 0, 0, -t148 * t83 - t166 * t77 + (t83 * t303 - t327 * t36) * t202 - t331 * t117, -t148 * t85 - t167 * t77 + (t85 * t303 - t327 * t35) * t202 + t332 * t117, 0, 0, 0, 0, 0, t166 * t12 + t230 * t34 + (qJD(5) * t228 + t332 * t206 - t330 * t210) * t82 + t331 * t56, t166 * t11 + t228 * t34 + (-qJD(5) * t230 + t330 * t206 + t332 * t210) * t82 + t331 * t58; 0, 0, 0, 0, 0, 0, 0, -t127 * t335, t127 ^ 2 - t335 ^ 2, t227 * t335 + t78, t127 * t296 + (t127 * t270 - t126) * qJD(1) - t232, -t156, -t81 * t127 - t227 * t55 + t234 + t239, -t227 * t54 - t335 * t81 + t223 - t233, t207 * t35 + t265 * t85, (t35 - t320) * t211 + (-t36 - t316) * t207, t117 * t265 - t127 * t85 + t207 * t77, -t117 ^ 2 * t207 + t127 * t83 + t211 * t77, -t117 * t127, -pkin(3) * t36 - t16 * t127 - t55 * t83 + (t54 * t117 + t236) * t207 + ((-t80 - t322) * t117 + t226) * t211, -pkin(3) * t35 + t323 * t117 + t17 * t127 - t55 * t85 + t236 * t211 + (t117 * t322 - t226) * t207, t11 * t207 * t210 + (-t207 * t298 + t210 * t299 - t72) * t58, t56 * t72 + t58 * t71 + (-t206 * t58 - t210 * t56) * t299 + (-t321 - t12 * t210 + (t206 * t56 - t210 * t58) * qJD(5)) * t207, -t72 * t82 + (-t11 + t291) * t211 + (t117 * t58 - t298 * t82 + t317) * t207, t71 * t82 + (t12 - t292) * t211 + (-t117 * t56 - t297 * t82 - t318) * t207, t117 * t207 * t82 - t211 * t34, -t13 * t71 - t24 * t56 + t328 * t210 + t214 * t206 + (t13 * t302 - t2 + (qJD(4) * t56 - t318) * pkin(9) - t217 * t210) * t211 + (t13 * t297 + t4 * t206 - t117 * t249 + (t12 + t292) * pkin(9)) * t207, -t13 * t72 - t24 * t58 - t328 * t206 + t214 * t210 + (t13 * t300 + t1 + (qJD(4) * t58 - t317) * pkin(9) + t217 * t206) * t211 + (-t13 * t298 + t4 * t210 - t117 * t8 + (t11 + t291) * pkin(9)) * t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t83, -t83 ^ 2 + t85 ^ 2, t35 + t320, t316 - t36, t77, -t41 * t85 - t235 - t275 + (-qJD(4) + t117) * t17, g(1) * t67 - g(2) * t65 + g(3) * t88 + t117 * t16 + t41 * t83 - t238, t267 * t58 + t321, (t11 - t325) * t210 + (-t12 - t324) * t206, t267 * t82 - t58 * t85 + t318, -t206 * t82 ^ 2 + t56 * t85 + t317, -t82 * t85, -pkin(4) * t12 - t17 * t56 + t224 * t206 - t329 * t210 + t249 * t85, -pkin(4) * t11 - t17 * t58 + t329 * t206 + t224 * t210 + t8 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t56, -t56 ^ 2 + t58 ^ 2, t11 + t325, -t12 + t324, t34, -g(1) * t37 - g(2) * t341 - g(3) * t243 - t13 * t58 + t8 * t82 + t2, g(1) * t38 - g(2) * t340 + g(3) * t62 + t13 * t56 - t249 * t82 - t1;];
tau_reg = t5;
