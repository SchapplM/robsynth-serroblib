% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 16:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:34:15
% EndTime: 2019-05-05 16:34:28
% DurationCPUTime: 5.66s
% Computational Cost: add. (11830->393), mult. (26335->520), div. (0->0), fcn. (17367->10), ass. (0->218)
t217 = qJD(3) ^ 2;
t205 = sin(pkin(10));
t215 = cos(qJ(3));
t207 = cos(pkin(10));
t212 = sin(qJ(3));
t258 = t212 * t207;
t172 = (t215 * t205 + t258) * qJD(1);
t281 = t172 ^ 2;
t153 = t281 + t217;
t255 = qJD(1) * t212;
t170 = -t207 * t215 * qJD(1) + t205 * t255;
t266 = t172 * t170;
t302 = qJDD(3) + t266;
t318 = t205 * t302;
t94 = t207 * t153 + t318;
t336 = pkin(3) * t94;
t335 = qJ(4) * t94;
t316 = t207 * t302;
t96 = -t205 * t153 + t316;
t334 = qJ(4) * t96;
t333 = t212 * t96 + t215 * t94;
t60 = t212 * t94 - t215 * t96;
t154 = t281 - t217;
t303 = qJDD(3) - t266;
t315 = t207 * t303;
t317 = t205 * t303;
t332 = t212 * (t205 * t154 + t315) - t215 * (t207 * t154 - t317);
t206 = sin(pkin(9));
t208 = cos(pkin(9));
t282 = t170 ^ 2;
t283 = -t282 - t281;
t251 = qJD(1) * qJD(3);
t241 = t215 * t251;
t249 = t212 * qJDD(1);
t181 = t241 + t249;
t242 = t212 * t251;
t248 = t215 * qJDD(1);
t228 = -t242 + t248;
t138 = t205 * t181 - t207 * t228;
t254 = qJD(3) * t172;
t104 = t138 - t254;
t139 = t207 * t181 + t205 * t228;
t161 = t170 * qJD(3);
t285 = t161 + t139;
t319 = -t207 * t104 + t205 * t285;
t320 = -t205 * t104 - t207 * t285;
t325 = -t212 * t320 + t215 * t319;
t331 = pkin(1) * (t206 * t325 - t208 * t283) + pkin(7) * t325 - pkin(2) * t283;
t330 = pkin(3) * t320;
t329 = qJ(4) * t320;
t326 = -pkin(3) * t283 + qJ(4) * t319;
t324 = t212 * t319 + t215 * t320;
t149 = t282 - t217;
t323 = t212 * (-t207 * t149 + t318) - t215 * (t205 * t149 + t316);
t286 = -t161 + t139;
t304 = t138 + t254;
t314 = t212 * (-t205 * t286 - t207 * t304) + t215 * (-t205 * t304 + t207 * t286);
t119 = -t217 - t282;
t84 = t205 * t119 + t315;
t313 = pkin(3) * t84;
t312 = qJ(4) * t84;
t87 = -t207 * t119 + t317;
t311 = qJ(4) * t87;
t301 = t212 * t87 - t215 * t84;
t300 = t212 * t84 + t215 * t87;
t297 = qJ(5) * t286;
t211 = sin(qJ(6));
t214 = cos(qJ(6));
t143 = t211 * qJD(3) - t214 * t170;
t145 = t214 * qJD(3) + t211 * t170;
t102 = t145 * t143;
t132 = qJDD(6) + t139;
t287 = -t102 + t132;
t291 = t211 * t287;
t290 = t214 * t287;
t289 = 2 * qJD(4);
t165 = qJD(6) + t172;
t115 = t165 * t143;
t92 = -t143 * qJD(6) + t214 * qJDD(3) + t211 * t138;
t288 = -t115 + t92;
t218 = qJD(1) ^ 2;
t213 = sin(qJ(1));
t216 = cos(qJ(1));
t234 = t216 * g(1) + t213 * g(2);
t179 = -t218 * pkin(1) - t234;
t233 = t213 * g(1) - t216 * g(2);
t227 = qJDD(1) * pkin(1) + t233;
t256 = t208 * t179 + t206 * t227;
t117 = -t218 * pkin(2) + qJDD(1) * pkin(7) + t256;
t202 = -g(3) + qJDD(2);
t111 = t212 * t117 - t215 * t202;
t191 = t215 * t218 * t212;
t186 = qJDD(3) + t191;
t82 = (-t181 + t241) * qJ(4) + t186 * pkin(3) - t111;
t112 = t215 * t117 + t212 * t202;
t185 = qJD(3) * pkin(3) - qJ(4) * t255;
t201 = t215 ^ 2;
t197 = t201 * t218;
t83 = -pkin(3) * t197 + qJ(4) * t228 - qJD(3) * t185 + t112;
t50 = -0.2e1 * qJD(4) * t170 + t205 * t82 + t207 * t83;
t284 = t281 - t282;
t141 = t143 ^ 2;
t142 = t145 ^ 2;
t163 = t165 ^ 2;
t280 = 2 * qJD(5);
t279 = -pkin(4) - pkin(8);
t278 = pkin(4) * t205;
t277 = pkin(4) * t207;
t237 = -t206 * t179 + t208 * t227;
t116 = -qJDD(1) * pkin(2) - t218 * pkin(7) - t237;
t89 = -t228 * pkin(3) - qJ(4) * t197 + t185 * t255 + qJDD(4) + t116;
t275 = t205 * t89;
t274 = t207 * t89;
t148 = t172 * pkin(5) - qJD(3) * pkin(8);
t118 = t170 * pkin(4) - t172 * qJ(5);
t230 = -t217 * pkin(4) - t170 * t118 + t50;
t250 = qJDD(3) * qJ(5);
t28 = t250 - t138 * pkin(5) - t282 * pkin(8) + (t280 + t148) * qJD(3) + t230;
t273 = t211 * t28;
t81 = t102 + t132;
t272 = t211 * t81;
t239 = t205 * t83 - t207 * t82;
t49 = t172 * t289 + t239;
t24 = t205 * t50 - t207 * t49;
t271 = t212 * t24;
t270 = t214 * t28;
t269 = t214 * t81;
t268 = t165 * t211;
t267 = t165 * t214;
t259 = t212 * t186;
t187 = qJDD(3) - t191;
t257 = t215 * t187;
t252 = t289 + t118;
t247 = -t142 - t163;
t246 = t205 * t102;
t245 = t207 * t102;
t244 = pkin(1) * t208 + pkin(2);
t243 = pkin(1) * t206 + pkin(7);
t240 = qJ(5) * t205 + pkin(3);
t25 = t205 * t49 + t207 * t50;
t229 = -qJDD(3) * pkin(4) - t217 * qJ(5) + qJDD(5) + t239;
t27 = -qJDD(3) * pkin(8) + t285 * pkin(5) + (pkin(8) * t170 + t252) * t172 + t229;
t221 = t138 * pkin(4) - t297 + t89;
t238 = pkin(4) * qJD(3) - (2 * qJD(5));
t33 = t221 + (-t148 + t238) * t172 + t138 * pkin(8) - t282 * pkin(5);
t17 = t211 * t33 - t214 * t27;
t74 = t212 * t111 + t215 * t112;
t236 = t211 * qJDD(3) - t214 * t138;
t18 = t211 * t27 + t214 * t33;
t7 = -t214 * t17 + t211 * t18;
t8 = t211 * t17 + t214 * t18;
t182 = -0.2e1 * t242 + t248;
t226 = (-qJD(6) + t165) * t145 - t236;
t225 = qJD(3) * t280 + t230;
t39 = t252 * t172 + t229;
t38 = t225 + t250;
t224 = t212 * (t207 * t139 - t205 * t254) + t215 * (t205 * t139 + t207 * t254);
t223 = t212 * (t205 * t138 + t207 * t161) + t215 * (-t207 * t138 + t205 * t161);
t220 = (t212 * (-t170 * t207 + t172 * t205) + t215 * (-t170 * t205 - t172 * t207)) * qJD(3);
t219 = -t172 * t280 + t221;
t200 = t212 ^ 2;
t196 = t200 * t218;
t190 = -t197 - t217;
t189 = -t196 - t217;
t184 = t196 + t197;
t183 = (t200 + t201) * qJDD(1);
t180 = 0.2e1 * t241 + t249;
t147 = -t212 * t189 - t257;
t146 = t215 * t190 - t259;
t114 = -t142 + t163;
t113 = t141 - t163;
t99 = t142 - t141;
t91 = -t145 * qJD(6) - t236;
t90 = -t163 - t141;
t88 = -t141 - t142;
t75 = (t143 * t211 + t145 * t214) * t165;
t68 = t115 + t92;
t64 = (qJD(6) + t165) * t145 + t236;
t62 = -t145 * t267 - t211 * t92;
t61 = -t143 * t268 - t214 * t91;
t58 = -t211 * t113 - t269;
t57 = -t214 * t114 - t291;
t56 = -t211 * t247 - t269;
t55 = t214 * t247 - t272;
t54 = t214 * t90 - t291;
t53 = t211 * t90 + t290;
t47 = t238 * t172 + t221;
t44 = t211 * t68 + t214 * t226;
t43 = t211 * t226 - t214 * t68;
t42 = t211 * t64 - t214 * t288;
t41 = t219 + (t304 + t254) * pkin(4);
t40 = -pkin(4) * t254 - t219 + t297;
t37 = t205 * t55 + t207 * t288;
t36 = t205 * t288 - t207 * t55;
t35 = t205 * t53 + t207 * t64;
t34 = t205 * t64 - t207 * t53;
t32 = -qJ(5) * t283 + t39;
t31 = -pkin(4) * t283 + t38;
t30 = t205 * t43 + t207 * t88;
t29 = t205 * t88 - t207 * t43;
t23 = pkin(5) * t43 - qJ(5) * t44;
t22 = t205 * t39 + t207 * t38;
t21 = t205 * t38 - t207 * t39;
t20 = -t212 * t36 + t215 * t37;
t19 = -t212 * t34 + t215 * t35;
t16 = -t212 * t29 + t215 * t30;
t14 = pkin(5) * t288 + t279 * t56 - t273;
t13 = pkin(5) * t64 + t279 * t54 + t270;
t12 = t215 * t25 - t271;
t11 = pkin(5) * t55 - qJ(5) * t56 - t18;
t10 = pkin(5) * t53 - qJ(5) * t54 - t17;
t6 = t205 * t7 + t207 * t28;
t5 = t205 * t28 - t207 * t7;
t4 = pkin(5) * t88 + t279 * t44 - t8;
t3 = pkin(5) * t7 - qJ(5) * t8;
t2 = pkin(5) * t28 + t279 * t8;
t1 = -t212 * t5 + t215 * t6;
t9 = [0, 0, 0, 0, 0, qJDD(1), t233, t234, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t208 * qJDD(1) - t206 * t218) + t237, pkin(1) * (-t206 * qJDD(1) - t208 * t218) - t256, 0, pkin(1) * (t206 * t256 + t208 * t237), (t181 + t241) * t212, t215 * t180 + t212 * t182, t259 + t215 * (-t196 + t217), t182 * t215, t212 * (t197 - t217) + t257, 0, -t215 * t116 + pkin(2) * t182 + pkin(7) * t146 + pkin(1) * (t206 * t146 + t208 * t182), t212 * t116 - pkin(2) * t180 + pkin(7) * t147 + pkin(1) * (t206 * t147 - t208 * t180), pkin(2) * t184 + pkin(7) * t183 + pkin(1) * (t206 * t183 + t208 * t184) + t74, -pkin(2) * t116 + pkin(7) * t74 + pkin(1) * (-t208 * t116 + t206 * t74), t224, t314, t332, t223, -t323, t220, t212 * (t275 - t312) + t215 * (-pkin(3) * t304 - t274 - t311) - pkin(2) * t304 - pkin(7) * t300 + pkin(1) * (-t206 * t300 - t208 * t304), t212 * (t274 + t335) + t215 * (-pkin(3) * t286 + t275 - t334) - pkin(2) * t286 + pkin(7) * t60 + pkin(1) * (t206 * t60 - t208 * t286), t212 * (-t24 - t329) + t215 * (t25 + t326) + t331, -qJ(4) * t271 + t215 * (-pkin(3) * t89 + qJ(4) * t25) - pkin(2) * t89 + pkin(7) * t12 + pkin(1) * (t206 * t12 - t208 * t89), t220, -t332, t323, t224, t314, t223, t212 * (-t205 * t31 + t207 * t32 - t329) + t215 * (t205 * t32 + t207 * t31 + t326) + t331, t212 * (-t205 * t41 + t312) + t215 * (t207 * t41 + t311) + t243 * t300 + (qJ(5) * t258 + t215 * t240 + t244) * t304, t212 * (t207 * t40 - t335) + t215 * (t205 * t40 + t334) - t243 * t60 + (-t212 * t278 + t215 * (pkin(3) + t277) + t244) * t286, (t212 * (-qJ(5) * t207 + t278) + t215 * (-t240 - t277) - t244) * t47 + (t243 + qJ(4)) * (-t212 * t21 + t215 * t22), t212 * (-t205 * t62 + t245) + t215 * (t207 * t62 + t246), t212 * (-t205 * t42 + t207 * t99) + t215 * (t205 * t99 + t207 * t42), t212 * (-t205 * t57 + t207 * t68) + t215 * (t205 * t68 + t207 * t57), t212 * (-t205 * t61 - t245) + t215 * (t207 * t61 - t246), t212 * (-t205 * t58 + t207 * t226) + t215 * (t205 * t226 + t207 * t58), t212 * (t207 * t132 - t205 * t75) + t215 * (t205 * t132 + t207 * t75), t212 * (-qJ(4) * t34 + t207 * t10 - t205 * t13) + t215 * (-pkin(3) * t54 + qJ(4) * t35 + t205 * t10 + t207 * t13) - pkin(2) * t54 + pkin(7) * t19 + pkin(1) * (t206 * t19 - t208 * t54), t212 * (-qJ(4) * t36 + t207 * t11 - t205 * t14) + t215 * (-pkin(3) * t56 + qJ(4) * t37 + t205 * t11 + t207 * t14) - pkin(2) * t56 + pkin(7) * t20 + pkin(1) * (t206 * t20 - t208 * t56), t212 * (-qJ(4) * t29 - t205 * t4 + t207 * t23) + t215 * (-pkin(3) * t44 + qJ(4) * t30 + t205 * t23 + t207 * t4) - pkin(2) * t44 + pkin(7) * t16 + pkin(1) * (t206 * t16 - t208 * t44), t212 * (-qJ(4) * t5 - t205 * t2 + t207 * t3) + t215 * (-pkin(3) * t8 + qJ(4) * t6 + t207 * t2 + t205 * t3) - pkin(2) * t8 + pkin(7) * t1 + pkin(1) * (t206 * t1 - t208 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, 0, 0, 0, 0, 0, 0, t215 * t186 + t212 * t190, -t212 * t187 + t215 * t189, 0, -t215 * t111 + t212 * t112, 0, 0, 0, 0, 0, 0, -t301, -t333, t324, t212 * t25 + t215 * t24, 0, 0, 0, 0, 0, 0, t324, t301, t333, t215 * t21 + t212 * t22, 0, 0, 0, 0, 0, 0, t212 * t35 + t215 * t34, t212 * t37 + t215 * t36, t212 * t30 + t215 * t29, t212 * t6 + t215 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, -t197 + t196, t249, t191, t248, qJDD(3), -t111, -t112, 0, 0, t266, t284, t285, -t266, -t104, qJDD(3), -t49 + t313, -t50 - t336, t330, pkin(3) * t24, qJDD(3), -t285, t104, t266, t284, -t266, -pkin(4) * t285 - qJ(5) * t104 + t330, -pkin(4) * t303 - qJ(5) * t119 - t313 + t39, t336 + pkin(4) * t153 + (qJDD(3) + t302) * qJ(5) + t225, pkin(3) * t21 - pkin(4) * t39 + qJ(5) * t38, -t145 * t268 + t214 * t92, -t211 * t288 - t214 * t64, -t211 * t114 + t290, t143 * t267 - t211 * t91, t214 * t113 - t272, (-t143 * t214 + t145 * t211) * t165, pkin(3) * t34 + qJ(5) * t64 + t279 * t53 + t273, pkin(3) * t36 + qJ(5) * t288 + t279 * t55 + t270, pkin(3) * t29 + qJ(5) * t88 + t279 * t43 - t7, pkin(3) * t5 + qJ(5) * t28 + t279 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t304, t286, t283, t89, 0, 0, 0, 0, 0, 0, t283, -t304, -t286, t47, 0, 0, 0, 0, 0, 0, t54, t56, t44, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, t303, -t153, t39, 0, 0, 0, 0, 0, 0, t53, t55, t43, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t99, t68, -t102, t226, t132, -t17, -t18, 0, 0;];
tauJ_reg  = t9;
