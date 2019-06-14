% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:07:26
% EndTime: 2019-05-05 15:07:42
% DurationCPUTime: 6.28s
% Computational Cost: add. (11867->360), mult. (27081->457), div. (0->0), fcn. (19213->8), ass. (0->241)
t204 = sin(pkin(9));
t205 = cos(pkin(9));
t208 = sin(qJ(4));
t211 = cos(qJ(4));
t270 = t205 * t208;
t235 = t204 * t211 + t270;
t186 = t235 * qJD(1);
t181 = qJD(5) + t186;
t292 = t181 ^ 2;
t271 = t204 * t208;
t188 = (t205 * t211 - t271) * qJD(1);
t207 = sin(qJ(5));
t210 = cos(qJ(5));
t169 = -t210 * qJD(4) + t207 * t188;
t293 = t169 ^ 2;
t140 = t293 - t292;
t255 = t188 * qJD(4);
t299 = t235 * qJDD(1);
t161 = -t299 - t255;
t151 = qJDD(5) - t161;
t171 = t207 * qJD(4) + t210 * t188;
t275 = t171 * t169;
t307 = t151 + t275;
t269 = t207 * t307;
t79 = -t210 * t140 + t269;
t144 = t181 * t171;
t252 = t205 * qJDD(1);
t253 = t204 * qJDD(1);
t185 = -t208 * t253 + t211 * t252;
t256 = t186 * qJD(4);
t163 = t185 - t256;
t245 = t210 * qJDD(4) - t207 * t163;
t230 = t171 * qJD(5) - t245;
t90 = -t144 + t230;
t361 = t204 * (t208 * t79 - t211 * t90) - t205 * (t208 * t90 + t211 * t79);
t263 = t210 * t307;
t358 = t207 * t140 + t263;
t168 = t171 ^ 2;
t305 = t168 - t293;
t234 = -t207 * qJDD(4) - t210 * t163;
t225 = -t169 * qJD(5) - t234;
t276 = t169 * t181;
t315 = t276 - t225;
t283 = t207 * t315;
t308 = t144 + t230;
t55 = t210 * t308 - t283;
t357 = t204 * (t208 * t55 + t211 * t305) - t205 * (-t208 * t305 + t211 * t55);
t306 = -t168 - t292;
t65 = -t210 * t306 + t269;
t356 = pkin(3) * t65;
t355 = pkin(4) * t65;
t354 = pkin(8) * t65;
t67 = t207 * t306 + t263;
t353 = pkin(8) * t67;
t352 = qJ(2) * t65;
t351 = t208 * t67;
t350 = t211 * t67;
t288 = pkin(1) + qJ(3);
t303 = -t275 + t151;
t103 = t210 * t303;
t300 = -t292 - t293;
t314 = t207 * t300 + t103;
t268 = t207 * t303;
t312 = t210 * t300 - t268;
t328 = t208 * t308 + t211 * t312;
t329 = t208 * t312 - t211 * t308;
t341 = t204 * t328 + t205 * t329;
t347 = qJ(2) * t314 - t288 * t341;
t346 = pkin(7) * t329;
t272 = t188 * t186;
t330 = qJDD(4) - t272;
t345 = t208 * t330;
t344 = t211 * t330;
t343 = t315 * qJ(6);
t342 = -pkin(3) * t314 + pkin(7) * t328;
t301 = t276 + t225;
t141 = -t168 + t292;
t332 = -t207 * t141 + t103;
t340 = t205 * (t208 * t301 + t211 * t332) - t204 * (t208 * t332 - t211 * t301);
t338 = pkin(4) * t314;
t337 = pkin(8) * t312;
t336 = pkin(8) * t314;
t331 = t210 * t141 + t268;
t304 = t168 + t293;
t327 = pkin(4) * t304;
t324 = t208 * t304;
t319 = t211 * t304;
t213 = qJD(1) ^ 2;
t209 = sin(qJ(1));
t212 = cos(qJ(1));
t239 = t209 * g(1) - t212 * g(2);
t232 = qJDD(2) - t239;
t224 = -t213 * qJ(2) + t232;
t251 = -0.2e1 * qJD(3) * qJD(1);
t316 = -t288 * qJDD(1) + t224 + t251;
t200 = t204 ^ 2;
t201 = t205 ^ 2;
t258 = t200 + t201;
t313 = pkin(3) * t253 - (t258 * pkin(7) + t288) * t213;
t311 = -t207 * t308 - t210 * t315;
t310 = -pkin(7) - t288;
t309 = t258 * t213;
t128 = t169 * pkin(5) - t171 * qJ(6);
t154 = t186 * pkin(4) - t188 * pkin(8);
t217 = (t251 + (-pkin(3) * t204 - qJ(2)) * t213 + t310 * qJDD(1) + t232) * t205;
t289 = t204 * g(3);
t216 = t217 + t289;
t153 = -t205 * g(3) + t316 * t204;
t136 = -t200 * t213 * pkin(3) - pkin(7) * t253 + t153;
t262 = t211 * t136;
t291 = qJD(4) ^ 2;
t63 = -t291 * pkin(4) + qJDD(4) * pkin(8) - t186 * t154 + t208 * t216 + t262;
t254 = qJD(2) * qJD(1);
t199 = 0.2e1 * t254;
t202 = qJDD(1) * qJ(2);
t240 = t212 * g(1) + t209 * g(2);
t233 = -t202 + t240;
t227 = -qJDD(3) + t233;
t76 = t199 + (-t163 + t256) * pkin(8) + (-t161 + t255) * pkin(4) - t227 + t313;
t40 = t207 * t76 + t210 * t63;
t244 = t151 * qJ(6) - t169 * t128 + t40;
t298 = -(t306 + t292) * pkin(5) + qJ(6) * t307 + t244;
t273 = t181 * t210;
t249 = t169 * t273;
t231 = t207 * t230 + t249;
t248 = t211 * t275;
t250 = t208 * t275;
t295 = t205 * (t211 * t231 - t250) - t204 * (t208 * t231 + t248);
t274 = t181 * t207;
t138 = t171 * t274;
t237 = t138 - t249;
t294 = t205 * (t208 * t151 + t211 * t237) - t204 * (-t211 * t151 + t208 * t237);
t183 = t186 ^ 2;
t184 = t188 ^ 2;
t290 = pkin(5) * t210;
t100 = g(3) * t271 + t208 * t217 + t262;
t99 = t208 * t136 - t211 * t216;
t59 = t208 * t100 - t211 * t99;
t287 = t205 * t59;
t62 = -qJDD(4) * pkin(4) - t291 * pkin(8) + t188 * t154 + t99;
t286 = t207 * t62;
t284 = t207 * t301;
t281 = t210 * t62;
t278 = qJ(6) * t210;
t277 = qJDD(1) * pkin(1);
t223 = t227 - 0.2e1 * t254;
t147 = t223 - t313;
t266 = t208 * t147;
t158 = qJDD(4) + t272;
t264 = t208 * t158;
t261 = t211 * t147;
t260 = t211 * t158;
t257 = qJD(6) * t181;
t247 = -pkin(4) * t211 - pkin(3);
t246 = -qJ(6) * t207 - pkin(4);
t39 = t207 * t63 - t210 * t76;
t18 = t207 * t39 + t210 * t40;
t60 = t211 * t100 + t208 * t99;
t173 = t288 * t213 + t223;
t243 = -t173 + t202;
t174 = 0.2e1 * t257;
t236 = t174 + t244;
t31 = -pkin(5) * t292 + t236;
t33 = -t151 * pkin(5) - qJ(6) * t292 + t171 * t128 + qJDD(6) + t39;
t242 = -pkin(5) * t33 + qJ(6) * t31;
t241 = -pkin(5) * t301 - qJ(6) * t90;
t238 = t169 * t274 - t210 * t230;
t4 = t205 * (t208 * t18 - t211 * t62) + t204 * (t211 * t18 + t208 * t62);
t17 = t207 * t40 - t210 * t39;
t109 = t205 * (t316 * t205 + t289) + t204 * t153;
t228 = (-t169 * t207 - t171 * t210) * t181;
t222 = t230 * pkin(5) + t343 + t62;
t221 = 0.2e1 * qJD(6) * t171 - t222;
t87 = t210 * t225 - t138;
t219 = t205 * (t211 * t87 + t250) - t204 * (t208 * t87 - t248);
t218 = pkin(5) * t303 + qJ(6) * t300 - t33;
t191 = t258 * qJDD(1);
t190 = t204 * t309;
t189 = t205 * t309;
t182 = -t224 + t277;
t177 = -t184 - t291;
t176 = -t184 + t291;
t175 = t183 - t291;
t162 = t185 - 0.2e1 * t256;
t160 = t299 + 0.2e1 * t255;
t155 = -t291 - t183;
t132 = -t183 - t184;
t124 = -t208 * t177 - t260;
t123 = t211 * t177 - t264;
t114 = t208 * t185 - t211 * t299;
t113 = -t211 * t185 - t208 * t299;
t112 = t211 * t155 - t345;
t111 = t208 * t155 + t344;
t97 = (qJD(5) + t181) * t169 + t234;
t92 = (-qJD(5) + t181) * t171 + t245;
t88 = t210 * t301;
t86 = t171 * t273 + t207 * t225;
t81 = t205 * t123 + t204 * t124;
t69 = t205 * t113 + t204 * t114;
t64 = t111 * t205 + t112 * t204;
t58 = t210 * t92 + t284;
t56 = -t210 * t90 + t284;
t54 = t207 * t92 - t88;
t53 = -t207 * t90 - t88;
t51 = -t208 * t97 - t350;
t49 = t211 * t97 - t351;
t47 = t208 * t315 + t350;
t45 = -t211 * t315 + t351;
t44 = t211 * t58 - t324;
t43 = t211 * t56 - t324;
t42 = t208 * t58 + t319;
t41 = t208 * t56 + t319;
t38 = t281 + t354;
t36 = t286 - t336;
t35 = (pkin(5) * t181 - 0.2e1 * qJD(6)) * t171 + t222;
t34 = t204 * t60 + t287;
t32 = -pkin(4) * t53 - t241;
t30 = t40 + t355;
t29 = t39 - t338;
t28 = (-t308 - t144) * pkin(5) + t221;
t27 = -pkin(5) * t144 + t221 - t343;
t26 = qJ(6) * t304 + t33;
t25 = (t304 - t292) * pkin(5) + t236;
t23 = t204 * t51 + t205 * t49;
t21 = t204 * t47 + t205 * t45;
t20 = t204 * t44 + t205 * t42;
t19 = t204 * t43 + t205 * t41;
t16 = -t218 - t338;
t15 = -0.2e1 * t257 - t298 - t355;
t14 = -t207 * t28 - t278 * t308 - t336;
t13 = pkin(5) * t283 + t210 * t27 - t354;
t10 = -pkin(8) * t54 - t17;
t9 = t207 * t33 + t210 * t31;
t8 = t207 * t31 - t210 * t33;
t7 = -pkin(8) * t53 - t207 * t25 + t210 * t26;
t6 = t208 * t35 + t211 * t9;
t5 = t208 * t9 - t211 * t35;
t3 = -pkin(8) * t8 + (pkin(5) * t207 - t278) * t35;
t2 = -pkin(4) * t8 - t242;
t1 = t204 * t6 + t205 * t5;
t11 = [0, 0, 0, 0, 0, qJDD(1), t239, t240, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t232 - 0.2e1 * t277, t199 + 0.2e1 * t202 - t240, pkin(1) * t182 + qJ(2) * (-t213 * pkin(1) + t199 - t233), t201 * qJDD(1), -0.2e1 * t204 * t252, 0, t200 * qJDD(1), 0, 0, t190 * t288 + t204 * t243, t189 * t288 + t205 * t243, -qJ(2) * t309 + t191 * t288 - t109, -qJ(2) * t173 - t109 * t288, t205 * (t163 * t211 - t208 * t255) - t204 * (t163 * t208 + t211 * t255), t205 * (-t160 * t211 - t162 * t208) - t204 * (-t160 * t208 + t162 * t211), t205 * (-t176 * t208 + t344) - t204 * (t176 * t211 + t345), t205 * (-t161 * t208 + t211 * t256) - t204 * (t161 * t211 + t208 * t256), t205 * (t175 * t211 - t264) - t204 * (t175 * t208 + t260), (t205 * (-t186 * t211 + t188 * t208) - t204 * (-t186 * t208 - t188 * t211)) * qJD(4), t205 * (-pkin(7) * t111 - t266) - t204 * (-pkin(3) * t160 + pkin(7) * t112 + t261) + qJ(2) * t160 - t288 * t64, t205 * (-pkin(7) * t123 - t261) - t204 * (-pkin(3) * t162 + pkin(7) * t124 - t266) + qJ(2) * t162 - t288 * t81, t205 * (-pkin(7) * t113 - t59) - t204 * (-pkin(3) * t132 + pkin(7) * t114 + t60) + qJ(2) * t132 - t288 * t69, -pkin(7) * t287 - t204 * (pkin(3) * t147 + pkin(7) * t60) - qJ(2) * t147 - t288 * t34, t219, t357, t340, t295, t361, t294, t205 * (-t208 * t29 + t211 * t36 - t346) - t204 * (t208 * t36 + t211 * t29 + t342) + t347, t205 * (-pkin(7) * t49 - t208 * t30 + t211 * t38) - t204 * (pkin(7) * t51 + t208 * t38 + t211 * t30 + t356) - t352 - t288 * t23, t205 * (-pkin(7) * t42 + t10 * t211) - t204 * (pkin(7) * t44 + t208 * t10) + (pkin(4) * t270 - t204 * t247 + qJ(2)) * t54 - t288 * t20, (t205 * (pkin(4) * t208 - pkin(8) * t211) - t204 * (-pkin(8) * t208 + t247) + qJ(2)) * t17 + t310 * t4, t219, t340, -t357, t294, -t361, t295, t205 * (t14 * t211 - t16 * t208 - t346) - t204 * (t14 * t208 + t16 * t211 + t342) + t347, t205 * (-pkin(7) * t41 - t208 * t32 + t211 * t7) - t204 * (-pkin(3) * t53 + pkin(7) * t43 + t208 * t7 + t211 * t32) + qJ(2) * t53 - t288 * t19, t205 * (-pkin(7) * t45 + t13 * t211 - t15 * t208) - t204 * (pkin(7) * t47 + t13 * t208 + t15 * t211 - t356) + t352 - t288 * t21, t205 * (-pkin(7) * t5 - t2 * t208 + t211 * t3) - t204 * (-pkin(3) * t8 + pkin(7) * t6 + t2 * t211 + t208 * t3) + qJ(2) * t8 - t288 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t213, -t182, 0, 0, 0, 0, 0, 0, -t190, -t189, -t191, t109, 0, 0, 0, 0, 0, 0, t64, t81, t69, t34, 0, 0, 0, 0, 0, 0, t341, t23, t20, t4, 0, 0, 0, 0, 0, 0, t341, t19, t21, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t252, -t309, -t173, 0, 0, 0, 0, 0, 0, t160, t162, t132, -t147, 0, 0, 0, 0, 0, 0, t314, -t65, t54, t17, 0, 0, 0, 0, 0, 0, t314, t53, t65, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t184 - t183, t185, -t272, -t299, qJDD(4), -t99, -t100, 0, 0, t86, t311, t331, t238, t358, t228, -pkin(4) * t308 - t281 + t337, pkin(4) * t97 + t286 - t353, pkin(8) * t58 + t18 + t327, -pkin(4) * t62 + pkin(8) * t18, t86, t331, -t311, t228, -t358, t238, t210 * t28 + t246 * t308 + t337, pkin(8) * t56 + t207 * t26 + t210 * t25 + t327, t353 + t207 * t27 - (pkin(4) + t290) * t315, pkin(8) * t9 + (t246 - t290) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, t305, t301, -t275, -t90, t151, -t39, -t40, 0, 0, t275, t301, -t305, t151, t90, -t275, t218, t241, t174 + t298, t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t303, t301, t306, t33;];
tauJ_reg  = t11;
