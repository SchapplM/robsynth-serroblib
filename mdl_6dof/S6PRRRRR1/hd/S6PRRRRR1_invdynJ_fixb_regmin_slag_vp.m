% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:39:55
% EndTime: 2019-03-09 00:40:09
% DurationCPUTime: 6.03s
% Computational Cost: add. (7520->443), mult. (17501->624), div. (0->0), fcn. (14381->18), ass. (0->250)
t214 = cos(qJ(6));
t302 = qJD(6) * t214;
t211 = sin(qJ(4));
t212 = sin(qJ(3));
t310 = qJD(2) * t212;
t286 = t211 * t310;
t216 = cos(qJ(4));
t217 = cos(qJ(3));
t308 = qJD(2) * t217;
t287 = t216 * t308;
t150 = -t286 + t287;
t151 = -t211 * t308 - t216 * t310;
t210 = sin(qJ(5));
t215 = cos(qJ(5));
t95 = t215 * t150 + t151 * t210;
t369 = t214 * t95;
t373 = t302 - t369;
t157 = t211 * t217 + t212 * t216;
t355 = pkin(8) + pkin(9);
t293 = qJD(3) * t355;
t162 = t212 * t293;
t163 = t217 * t293;
t218 = cos(qJ(2));
t207 = sin(pkin(6));
t312 = qJD(1) * t207;
t291 = t218 * t312;
t168 = t355 * t212;
t169 = t355 * t217;
t315 = -t211 * t168 + t216 * t169;
t372 = -qJD(4) * t315 + t157 * t291 + t211 * t162 - t216 * t163;
t156 = t211 * t212 - t216 * t217;
t306 = qJD(4) * t216;
t307 = qJD(4) * t211;
t371 = -t156 * t291 + t216 * t162 + t211 * t163 + t168 * t306 + t169 * t307;
t318 = -qJD(6) + t95;
t370 = t318 + qJD(6);
t209 = sin(qJ(6));
t201 = qJDD(3) + qJDD(4);
t196 = qJDD(5) + t201;
t202 = qJD(3) + qJD(4);
t197 = qJD(5) + t202;
t250 = t150 * t210 - t215 * t151;
t303 = qJD(6) * t209;
t304 = qJD(5) * t215;
t305 = qJD(5) * t210;
t300 = qJD(2) * qJD(3);
t283 = t217 * t300;
t297 = t217 * qJDD(2);
t298 = t212 * qJDD(2);
t78 = qJD(4) * t287 - t202 * t286 + t211 * t297 + (t283 + t298) * t216;
t120 = t202 * t157;
t254 = t211 * t298 - t216 * t297;
t79 = qJD(2) * t120 + t254;
t42 = t150 * t304 + t151 * t305 - t210 * t79 + t215 * t78;
t23 = t209 * t196 + t197 * t302 + t214 * t42 - t250 * t303;
t21 = t23 * t209;
t82 = t197 * t209 + t214 * t250;
t9 = t373 * t82 + t21;
t43 = qJD(5) * t250 + t210 * t78 + t215 * t79;
t41 = qJDD(6) + t43;
t342 = t209 * t41 - t302 * t318;
t8 = -t250 * t82 + t318 * t369 + t342;
t366 = t209 * t318;
t39 = t214 * t41;
t80 = -t214 * t197 + t209 * t250;
t7 = t250 * t80 - t318 * t366 + t39;
t22 = t23 * t214;
t24 = qJD(6) * t82 - t214 * t196 + t209 * t42;
t4 = -t209 * t24 + t366 * t82 - t373 * t80 + t22;
t213 = sin(qJ(2));
t292 = t213 * t312;
t274 = t355 * qJD(2) + t292;
t208 = cos(pkin(6));
t311 = qJD(1) * t208;
t124 = t212 * t311 + t217 * t274;
t117 = t216 * t124;
t123 = -t274 * t212 + t217 * t311;
t338 = qJD(3) * pkin(3);
t118 = t123 + t338;
t251 = -t211 * t118 - t117;
t353 = pkin(10) * t150;
t62 = -t251 + t353;
t336 = t210 * t62;
t142 = t151 * pkin(10);
t115 = t211 * t124;
t264 = t216 * t118 - t115;
t61 = t142 + t264;
t53 = pkin(4) * t202 + t61;
t31 = t215 * t53 - t336;
t29 = -pkin(5) * t197 - t31;
t349 = t29 * t95;
t330 = cos(pkin(12));
t276 = t330 * t213;
t206 = sin(pkin(12));
t324 = t206 * t218;
t144 = t208 * t276 + t324;
t275 = t330 * t218;
t325 = t206 * t213;
t146 = -t208 * t325 + t275;
t205 = qJ(3) + qJ(4);
t200 = qJ(5) + t205;
t188 = sin(t200);
t189 = cos(t200);
t277 = t207 * t330;
t323 = t207 * t213;
t326 = t206 * t207;
t240 = -g(3) * (-t188 * t323 + t189 * t208) - g(2) * (-t144 * t188 - t189 * t277) - g(1) * (-t146 * t188 + t189 * t326);
t299 = t208 * qJDD(1);
t177 = t217 * t299;
t301 = qJD(1) * qJD(2);
t132 = qJDD(2) * pkin(8) + (qJDD(1) * t213 + t218 * t301) * t207;
t268 = pkin(9) * qJDD(2) + t132;
t67 = qJDD(3) * pkin(3) - qJD(3) * t124 - t268 * t212 + t177;
t68 = qJD(3) * t123 + t212 * t299 + t268 * t217;
t232 = qJD(4) * t251 - t211 * t68 + t216 * t67;
t15 = pkin(4) * t201 - pkin(10) * t78 + t232;
t356 = -(qJD(4) * t118 + t68) * t216 + t124 * t307 - t211 * t67;
t18 = -pkin(10) * t79 - t356;
t332 = t215 * t62;
t32 = t210 * t53 + t332;
t360 = -qJD(5) * t32 + t215 * t15 - t210 * t18;
t3 = -pkin(5) * t196 - t360;
t238 = t240 - t3;
t346 = t250 * t95;
t119 = t202 * t156;
t368 = -pkin(10) * t119 - t372;
t367 = -pkin(10) * t120 - t371;
t239 = -t212 * t338 + t292;
t143 = -t208 * t275 + t325;
t145 = t208 * t324 + t276;
t259 = g(1) * t145 + g(2) * t143;
t321 = t207 * t218;
t237 = g(3) * t321 - t259;
t234 = t237 * t189;
t113 = t215 * t156 + t157 * t210;
t114 = -t156 * t210 + t157 * t215;
t193 = -pkin(3) * t217 - pkin(2);
t134 = pkin(4) * t156 + t193;
t55 = pkin(5) * t113 - pkin(11) * t114 + t134;
t365 = t55 * t41 - t234;
t37 = t250 ^ 2 - t95 ^ 2;
t65 = pkin(5) * t250 - pkin(11) * t95;
t27 = -t197 * t95 + t42;
t101 = t144 * t189 - t188 * t277;
t103 = t146 * t189 + t188 * t326;
t141 = qJD(2) * t193 - t291;
t105 = -pkin(4) * t150 + t141;
t130 = t188 * t208 + t189 * t323;
t357 = (qJD(5) * t53 + t18) * t215 + t210 * t15 - t62 * t305;
t225 = g(1) * t103 + g(2) * t101 + g(3) * t130 - t105 * t95 - t357;
t347 = t318 * t250;
t256 = pkin(4) * t120 - t239;
t30 = pkin(11) * t197 + t32;
t44 = -pkin(5) * t95 - pkin(11) * t250 + t105;
t253 = t209 * t30 - t214 * t44;
t285 = t250 * t253 + t29 * t303;
t11 = t209 * t44 + t214 * t30;
t242 = t11 * t250 - t238 * t209 + t29 * t302;
t223 = -t105 * t250 + t240 + t360;
t2 = pkin(11) * t196 + t357;
t258 = g(1) * t146 + g(2) * t144;
t262 = -t216 * t168 - t169 * t211;
t87 = -pkin(10) * t157 + t262;
t88 = -pkin(10) * t156 + t315;
t56 = t210 * t88 - t215 * t87;
t345 = qJD(5) * t56 + t368 * t210 - t367 * t215;
t48 = -qJD(5) * t113 - t119 * t215 - t120 * t210;
t57 = t210 * t87 + t215 * t88;
t359 = -(qJD(6) * t44 + t2) * t113 + t29 * t48 + t3 * t114 - (-qJD(6) * t55 + t345) * t318 - t57 * t41 - g(3) * t323 - t258;
t28 = t197 * t250 - t43;
t219 = qJD(3) ^ 2;
t284 = t213 * t301;
t255 = -qJDD(1) * t321 + t207 * t284;
t358 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t219 + t207 * (-g(3) * t218 + t284) - t255 + t259;
t354 = pkin(4) * t151;
t344 = qJD(5) * t57 + t367 * t210 + t368 * t215;
t192 = pkin(3) * t216 + pkin(4);
t320 = t210 * t211;
t263 = -t123 * t211 - t117;
t63 = t263 - t353;
t316 = t216 * t123 - t115;
t64 = t142 + t316;
t341 = t210 * t63 + t215 * t64 - t192 * t304 - (-t211 * t305 + (t215 * t216 - t320) * qJD(4)) * pkin(3);
t319 = t211 * t215;
t340 = -t210 * t64 + t215 * t63 + t192 * t305 + (t211 * t304 + (t210 * t216 + t319) * qJD(4)) * pkin(3);
t339 = qJD(2) * pkin(2);
t335 = t214 * t82;
t331 = t29 * t114;
t327 = t151 * t150;
t322 = t207 * t217;
t317 = qJDD(1) - g(3);
t314 = pkin(3) * t319 + t210 * t192;
t203 = t212 ^ 2;
t313 = -t217 ^ 2 + t203;
t309 = qJD(2) * t213;
t194 = pkin(3) * t310;
t295 = t209 * t321;
t294 = t214 * t321;
t289 = t207 * t309;
t288 = qJD(2) * t321;
t282 = t218 * t300;
t140 = pkin(11) + t314;
t50 = -t354 + t65;
t266 = qJD(6) * t140 + t194 + t50;
t190 = pkin(4) * t210 + pkin(11);
t265 = qJD(6) * t190 + t50;
t33 = t210 * t61 + t332;
t260 = pkin(4) * t305 - t33;
t49 = qJD(5) * t114 - t119 * t210 + t215 * t120;
t257 = pkin(5) * t49 - pkin(11) * t48 + t256;
t147 = t208 * t217 - t212 * t323;
t148 = t208 * t212 + t213 * t322;
t90 = t147 * t216 - t148 * t211;
t91 = t147 * t211 + t148 * t216;
t59 = t210 * t91 - t215 * t90;
t60 = t210 * t90 + t215 * t91;
t220 = qJD(2) ^ 2;
t249 = qJDD(2) * t218 - t213 * t220;
t248 = -pkin(3) * t320 + t192 * t215;
t247 = -g(1) * t206 + g(2) * t330;
t246 = -t303 * t318 - t39;
t245 = -t209 * t60 - t294;
t244 = -t214 * t60 + t295;
t165 = -t291 - t339;
t236 = -qJD(2) * t165 - t132 + t258;
t235 = -t140 * t41 - t318 * t341 - t349;
t34 = t215 * t61 - t336;
t229 = -t190 * t41 - t349 - (-pkin(4) * t304 + t34) * t318;
t106 = qJD(3) * t194 + qJDD(2) * t193 + t255;
t228 = -pkin(8) * qJDD(3) + (t165 + t291 - t339) * qJD(3);
t58 = pkin(4) * t79 + t106;
t198 = sin(t205);
t199 = cos(t205);
t224 = -g(1) * (-t146 * t199 - t198 * t326) - g(2) * (-t144 * t199 + t198 * t277) - g(3) * (-t198 * t208 - t199 * t323) - t141 * t150 + t356;
t222 = -g(1) * (-t146 * t198 + t199 * t326) - g(2) * (-t144 * t198 - t199 * t277) - g(3) * (-t198 * t323 + t199 * t208) + t141 * t151 + t232;
t191 = -pkin(4) * t215 - pkin(5);
t139 = -pkin(5) - t248;
t128 = t194 - t354;
t122 = -qJD(3) * t148 - t212 * t288;
t121 = qJD(3) * t147 + t217 * t288;
t83 = -t150 ^ 2 + t151 ^ 2;
t70 = -t254 + (-qJD(2) * t157 - t151) * t202;
t69 = -t150 * t202 + t78;
t46 = -qJD(4) * t91 - t121 * t211 + t122 * t216;
t45 = qJD(4) * t90 + t121 * t216 + t122 * t211;
t13 = qJD(5) * t60 + t210 * t45 - t215 * t46;
t12 = -qJD(5) * t59 + t210 * t46 + t215 * t45;
t6 = pkin(5) * t43 - pkin(11) * t42 + t58;
t5 = t214 * t6;
t1 = [t317, 0, t249 * t207 (-qJDD(2) * t213 - t218 * t220) * t207, 0, 0, 0, 0, 0, qJD(3) * t122 + qJDD(3) * t147 + (-t212 * t282 + t217 * t249) * t207, -qJD(3) * t121 - qJDD(3) * t148 + (-t212 * t249 - t217 * t282) * t207, 0, 0, 0, 0, 0, t201 * t90 + t202 * t46 + (-t150 * t309 - t218 * t79) * t207, -t201 * t91 - t202 * t45 + (-t151 * t309 - t218 * t78) * t207, 0, 0, 0, 0, 0, -t13 * t197 - t196 * t59 + (-t218 * t43 - t309 * t95) * t207, -t12 * t197 - t196 * t60 + (-t218 * t42 + t250 * t309) * t207, 0, 0, 0, 0, 0 -(qJD(6) * t244 - t12 * t209 + t214 * t289) * t318 + t245 * t41 + t13 * t80 + t59 * t24 (qJD(6) * t245 + t12 * t214 + t209 * t289) * t318 + t244 * t41 + t13 * t82 + t59 * t23; 0, qJDD(2), t317 * t321 + t259, -t317 * t323 + t258, qJDD(2) * t203 + 0.2e1 * t212 * t283, 0.2e1 * t212 * t297 - 0.2e1 * t300 * t313, qJDD(3) * t212 + t217 * t219, qJDD(3) * t217 - t212 * t219, 0, t228 * t212 + t217 * t358, -t212 * t358 + t228 * t217, t119 * t151 + t157 * t78, -t119 * t150 + t120 * t151 - t156 * t78 - t157 * t79, -t119 * t202 + t157 * t201, -t120 * t202 - t156 * t201, 0, t106 * t156 + t141 * t120 + t239 * t150 + t193 * t79 - t237 * t199 + t262 * t201 + t372 * t202, t106 * t157 - t141 * t119 + t239 * t151 + t193 * t78 + t237 * t198 - t315 * t201 + t371 * t202, t114 * t42 + t250 * t48, -t113 * t42 - t114 * t43 - t250 * t49 + t48 * t95, t114 * t196 + t197 * t48, -t113 * t196 - t197 * t49, 0, t105 * t49 + t113 * t58 + t134 * t43 - t196 * t56 - t197 * t344 - t256 * t95 - t234, t105 * t48 + t114 * t58 + t134 * t42 + t188 * t237 - t196 * t57 + t197 * t345 + t250 * t256, t48 * t335 + (-t303 * t82 + t22) * t114 (-t209 * t82 - t214 * t80) * t48 + (-t21 - t214 * t24 + (t209 * t80 - t335) * qJD(6)) * t114, -t214 * t318 * t48 + t113 * t23 - t114 * t246 + t49 * t82, -t113 * t24 - t114 * t342 + t366 * t48 - t49 * t80, t113 * t41 - t318 * t49, -t253 * t49 + t5 * t113 + t56 * t24 + t344 * t80 + (-t257 * t318 + (-t30 * t113 + t318 * t57 + t331) * qJD(6) + t365) * t214 + t359 * t209, -t11 * t49 + t56 * t23 + t344 * t82 + (-(-qJD(6) * t30 + t6) * t113 - qJD(6) * t331 - (qJD(6) * t57 - t257) * t318 - t365) * t209 + t359 * t214; 0, 0, 0, 0, -t212 * t220 * t217, t313 * t220, t298, t297, qJDD(3), -g(3) * t147 + t212 * t236 + t247 * t322 + t177, g(3) * t148 + (-t207 * t247 - t299) * t212 + t236 * t217, t327, t83, t69, t70, t201, -t263 * t202 + (t150 * t310 + t216 * t201 - t202 * t307) * pkin(3) + t222, t316 * t202 + (t151 * t310 - t211 * t201 - t202 * t306) * pkin(3) + t224, -t346, t37, t27, t28, t196, t128 * t95 + t196 * t248 - t197 * t340 + t223, -t128 * t250 - t196 * t314 + t197 * t341 + t225, t9, t4, t8, t7, t347, t139 * t24 + t340 * t80 + t235 * t209 + (t266 * t318 + t238) * t214 + t285, t139 * t23 + t214 * t235 - t266 * t366 + t340 * t82 + t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t327, t83, t69, t70, t201, -t202 * t251 + t222, t202 * t264 + t224, -t346, t37, t27, t28, t196, t197 * t33 + (-t151 * t95 + t196 * t215 - t197 * t305) * pkin(4) + t223, t197 * t34 + (t151 * t250 - t196 * t210 - t197 * t304) * pkin(4) + t225, t9, t4, t8, t7, t347, t191 * t24 + t260 * t80 + t229 * t209 + (t265 * t318 + t238) * t214 + t285, t191 * t23 + t214 * t229 + t260 * t82 - t265 * t366 + t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t346, t37, t27, t28, t196, t197 * t32 + t223, t197 * t31 + t225, t9, t4, t8, t7, t347, -pkin(5) * t24 - t32 * t80 + (-pkin(11) * t41 - t31 * t318 - t349) * t209 + (-(-pkin(11) * qJD(6) - t65) * t318 + t238) * t214 + t285, -pkin(5) * t23 - (t209 * t65 + t214 * t31) * t318 - t32 * t82 - t29 * t369 + t246 * pkin(11) + t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * t80, -t80 ^ 2 + t82 ^ 2, -t318 * t80 + t23, -t318 * t82 - t24, t41, -t209 * t2 + t5 - t29 * t82 - g(1) * (-t103 * t209 + t145 * t214) - g(2) * (-t101 * t209 + t143 * t214) - g(3) * (-t130 * t209 - t294) - t370 * t11, -t214 * t2 - t209 * t6 + t29 * t80 - g(1) * (-t103 * t214 - t145 * t209) - g(2) * (-t101 * t214 - t143 * t209) - g(3) * (-t130 * t214 + t295) + t370 * t253;];
tau_reg  = t1;
