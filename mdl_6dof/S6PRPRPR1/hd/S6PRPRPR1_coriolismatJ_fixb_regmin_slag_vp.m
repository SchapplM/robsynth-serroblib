% Calculate minimal parameter regressor of coriolis matrix for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRPRPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:27:59
% EndTime: 2019-03-08 19:28:06
% DurationCPUTime: 2.75s
% Computational Cost: add. (3536->232), mult. (8661->396), div. (0->0), fcn. (10216->12), ass. (0->211)
t206 = sin(qJ(6));
t210 = cos(qJ(4));
t203 = sin(pkin(11));
t196 = pkin(2) * t203 + pkin(8);
t294 = qJ(5) + t196;
t179 = t294 * t210;
t202 = sin(pkin(12));
t207 = sin(qJ(4));
t250 = t294 * t207;
t311 = cos(pkin(12));
t217 = t179 * t311 - t202 * t250;
t310 = t217 * t206;
t209 = cos(qJ(6));
t309 = t217 * t209;
t182 = t202 * t207 - t210 * t311;
t127 = t206 * t182;
t120 = t127 * qJD(6);
t251 = t311 * t207;
t304 = t202 * t210;
t184 = t251 + t304;
t287 = qJD(4) * t209;
t174 = t184 * t287;
t338 = -t174 + t120;
t180 = t182 ^ 2;
t181 = t184 ^ 2;
t337 = -t181 - t180;
t200 = t206 ^ 2;
t201 = t209 ^ 2;
t192 = t201 - t200;
t129 = t206 * t184;
t249 = 0.2e1 * t209 * t129;
t226 = qJD(2) * t249 - qJD(4) * t192;
t204 = sin(pkin(6));
t208 = sin(qJ(2));
t312 = cos(pkin(11));
t253 = t312 * t208;
t326 = cos(qJ(2));
t265 = t326 * t203;
t159 = (t265 + t253) * t204;
t205 = cos(pkin(6));
t141 = t159 * t210 + t205 * t207;
t236 = t159 * t207 - t205 * t210;
t75 = t141 * t202 + t236 * t311;
t336 = -t75 / 0.2e1;
t334 = -t217 / 0.2e1;
t333 = t182 / 0.2e1;
t332 = -t184 / 0.2e1;
t331 = t184 / 0.2e1;
t328 = -t209 / 0.2e1;
t327 = t209 / 0.2e1;
t325 = t207 * pkin(4);
t20 = (t334 + t217 / 0.2e1) * t182;
t324 = t20 * qJD(4);
t215 = t141 * t311 - t202 * t236;
t316 = t215 * t182;
t317 = t75 * t184;
t323 = -t316 / 0.2e1 + t317 / 0.2e1;
t322 = qJD(4) * pkin(4);
t266 = t204 * t326;
t303 = t204 * t208;
t158 = t203 * t303 - t266 * t312;
t95 = t182 * t158;
t321 = t206 * t95;
t320 = t209 * t95;
t47 = -t158 * t209 + t206 * t215;
t319 = t47 * t182;
t48 = t158 * t206 + t209 * t215;
t318 = t48 * t182;
t271 = t75 / 0.2e1 + t336;
t247 = t271 * t182;
t315 = t247 * qJD(2);
t94 = t184 * t158;
t314 = t94 * t206;
t313 = t94 * t209;
t100 = t158 * t207;
t308 = t159 * t206;
t307 = t159 * t209;
t306 = t182 * t202;
t305 = t184 * t209;
t132 = pkin(5) * t184 + pkin(9) * t182 + t325;
t302 = t206 * t132;
t301 = t209 * t132;
t36 = t331 * t95 - t333 * t94;
t300 = t36 * qJD(1);
t273 = t181 - t180;
t96 = t273 * t206;
t297 = t96 * qJD(2);
t97 = t337 * t206;
t296 = t97 * qJD(2);
t98 = t273 * t209;
t295 = t98 * qJD(2);
t292 = qJD(2) * t184;
t291 = qJD(2) * t207;
t290 = qJD(2) * t209;
t289 = qJD(2) * t210;
t288 = qJD(4) * t206;
t286 = qJD(5) * t209;
t285 = qJD(6) * t206;
t284 = qJD(6) * t209;
t252 = t311 * t184;
t222 = -t306 / 0.2e1 - t252 / 0.2e1;
t106 = (-t207 / 0.2e1 + t222) * pkin(4);
t283 = t106 * qJD(2);
t119 = t127 * qJD(2);
t282 = t129 * qJD(2);
t131 = t209 * t182;
t281 = t131 * qJD(2);
t280 = t131 * qJD(4);
t134 = t337 * t209;
t279 = t134 * qJD(2);
t278 = t337 * qJD(2);
t178 = t251 / 0.2e1 + t304 / 0.2e1;
t277 = t178 * qJD(2);
t193 = -t207 ^ 2 + t210 ^ 2;
t276 = t193 * qJD(2);
t275 = t207 * qJD(4);
t274 = t210 * qJD(4);
t272 = t325 / 0.2e1;
t270 = -t321 / 0.2e1;
t269 = -t320 / 0.2e1;
t268 = -t317 / 0.2e1;
t267 = t75 * t327;
t264 = t184 * t285;
t263 = t184 * t284;
t262 = t182 * t292;
t261 = t182 * t184 * qJD(4);
t198 = -pkin(2) * t312 - pkin(3);
t260 = t198 * t291;
t259 = t198 * t289;
t258 = t206 * t284;
t257 = t206 * t287;
t256 = t207 * t289;
t255 = t184 * t290;
t114 = t179 * t202 + t250 * t311;
t248 = -qJD(2) * t182 - qJD(6);
t246 = t159 / 0.2e1 + t268;
t244 = qJD(4) * t249;
t188 = -pkin(4) * t210 + t198;
t214 = pkin(5) * t182 - pkin(9) * t184 + t188;
t59 = -t209 * t214 + t310;
t11 = (-t59 + t310) * t184 + t301 * t182;
t220 = t215 * t331 + t247;
t212 = t206 * t220 + t332 * t47;
t4 = -t313 / 0.2e1 + t212;
t243 = t4 * qJD(1) + t11 * qJD(2);
t10 = t215 * t333 - t331 * t75 + t323;
t242 = -t10 * qJD(1) - t20 * qJD(2);
t14 = t270 + t318 / 0.2e1 + t246 * t209;
t60 = t206 * t214 + t309;
t31 = t114 * t305 - t182 * t60;
t241 = qJD(1) * t14 - qJD(2) * t31;
t15 = t269 - t319 / 0.2e1 - t246 * t206;
t30 = -t114 * t129 + t182 * t59;
t240 = -qJD(1) * t15 + qJD(2) * t30;
t216 = (t253 / 0.2e1 + t265 / 0.2e1) * t204;
t22 = t316 / 0.2e1 + t268 + t216;
t45 = t114 * t184 - t182 * t217;
t239 = -qJD(1) * t22 + qJD(2) * t45;
t238 = t10 * qJD(3);
t13 = t158 * t159 + t215 * t95 - t75 * t94;
t237 = qJD(1) * t13 + qJD(3) * t36;
t195 = pkin(4) * t202 + pkin(9);
t197 = -pkin(4) * t311 - pkin(5);
t235 = -t182 * t197 - t184 * t195;
t234 = t248 * t209;
t233 = t195 * t333 + t197 * t332;
t232 = t184 * t234;
t125 = (t200 / 0.2e1 - t201 / 0.2e1) * t184;
t231 = -qJD(2) * t125 + t257;
t230 = qJD(6) * t178 + t262;
t229 = t95 * t202 / 0.2e1 + t94 * t311 / 0.2e1;
t228 = t181 * t206 * t290 + qJD(4) * t125;
t133 = t192 * t181;
t227 = qJD(2) * t133 + t244;
t213 = t217 * t336 - t334 * t75;
t1 = (-t100 / 0.2e1 + t229) * pkin(4) + t213;
t28 = t188 * t325;
t225 = -t1 * qJD(1) + t28 * qJD(2) + t20 * qJD(3);
t12 = (-t60 + t309) * t184 - t302 * t182;
t211 = t209 * t220 + t332 * t48;
t7 = t314 / 0.2e1 + t211;
t224 = t7 * qJD(1) + t12 * qJD(2);
t223 = qJD(1) * t247;
t221 = t132 / 0.2e1 + t233;
t26 = t221 * t209;
t32 = t271 * t206;
t219 = qJD(1) * t32 + qJD(2) * t26 - t197 * t288;
t24 = t221 * t206;
t33 = t271 * t209;
t218 = qJD(1) * t33 - qJD(2) * t24 - t197 * t287;
t175 = t178 * qJD(4);
t173 = t184 * t288;
t124 = t131 * qJD(6);
t118 = t127 * qJD(4);
t117 = t125 * qJD(6);
t105 = pkin(4) * t222 + t272;
t104 = -t119 - t285;
t102 = t158 * t210;
t35 = -t328 * t75 + t267;
t34 = t75 * t206;
t27 = t301 / 0.2e1 - t233 * t209 + t114 * t206;
t25 = -t302 / 0.2e1 + t233 * t206 + (t327 - t328) * t114;
t23 = t216 + t323;
t17 = -t318 / 0.2e1 + t184 * t267 + t270 + t307 / 0.2e1;
t16 = t319 / 0.2e1 + t206 * t268 + t269 - t308 / 0.2e1;
t8 = t247 * qJD(4);
t6 = -t314 / 0.2e1 + t211;
t5 = t313 / 0.2e1 + t212;
t3 = qJD(2) * t36 + qJD(4) * t10;
t2 = pkin(4) * t229 + t158 * t272 - t213;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, -qJD(2) * t303, -qJD(2) * t266 (-t158 * t203 - t159 * t312) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, qJD(4) * t100 - t159 * t289, qJD(4) * t102 + t159 * t291 (-t182 * t95 - t184 * t94) * qJD(2) + t8 (-t114 * t94 + t159 * t188 + t217 * t95) * qJD(2) + t2 * qJD(4) + t23 * qJD(5) + t237, 0, 0, 0, 0, 0 ((t307 - t321) * t182 - t94 * t129) * qJD(2) + t5 * qJD(4) + t17 * qJD(6) (-(t308 + t320) * t182 - t94 * t305) * qJD(2) + t6 * qJD(4) + t16 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t100 - qJD(4) * t141, qJD(2) * t102 + qJD(4) * t236, t315, t2 * qJD(2) + (-t202 * t75 - t215 * t311) * t322 + t238, 0, 0, 0, 0, 0, qJD(2) * t5 + qJD(6) * t34 - t215 * t287, qJD(2) * t6 + qJD(6) * t35 + t215 * t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t17 + qJD(4) * t34 - qJD(6) * t48, qJD(2) * t16 + qJD(4) * t35 + qJD(6) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -qJD(4) * t1 - qJD(5) * t22 - t237, 0, 0, 0, 0, 0, qJD(4) * t4 - qJD(6) * t14, qJD(4) * t7 - qJD(6) * t15; 0, 0, 0, 0, 0, t207 * t274, t193 * qJD(4), 0, 0, 0, t198 * t275, t198 * t274, -t337 * qJD(5), qJD(4) * t28 + qJD(5) * t45, -t181 * t258 - t201 * t261, -qJD(6) * t133 + t182 * t244, qJD(4) * t98 - t182 * t264, -qJD(4) * t96 - t182 * t263, t261, qJD(4) * t11 - qJD(5) * t97 + qJD(6) * t31, qJD(4) * t12 - qJD(5) * t134 + qJD(6) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t300 + t324, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t256, t276, t274, -t275, 0, -t196 * t274 + t260, t196 * t275 + t259 (t182 * t311 - t184 * t202) * t322 + t223 (-t114 * t202 - t217 * t311) * t322 + t105 * qJD(5) + t225, -t117 + (-t201 * t292 - t257) * t182, t182 * t226 - 0.2e1 * t184 * t258, t173 + t295, t174 - t297, t230 (t206 * t235 - t309) * qJD(4) + t27 * qJD(6) + t243 (t209 * t235 + t310) * qJD(4) + t25 * qJD(6) + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, qJD(4) * t105 + t239, 0, 0, 0, 0, 0, -t296, -t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228, -t227, t248 * t129, t232, t175, qJD(4) * t27 - qJD(6) * t60 - t241, qJD(4) * t25 + qJD(6) * t59 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t300 + t324, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t275, -t274, 0 (-t252 - t306) * t322 - t242, 0, 0, 0, 0, 0, t338, t124 + t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118 - t263, t264 + t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t315, qJD(2) * t1 - t238, 0, 0, 0, 0, 0, -qJD(2) * t4 - qJD(6) * t32, -qJD(2) * t7 - qJD(6) * t33; 0, 0, 0, 0, 0, -t256, -t276, 0, 0, 0, -t260, -t259, -t223, qJD(5) * t106 - t225, t201 * t262 - t117, 0.2e1 * t206 * t232, t124 - t295, -t120 + t297, -t230, -qJD(6) * t26 - t184 * t286 - t243, qJD(5) * t129 + qJD(6) * t24 - t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t258, t192 * qJD(6), 0, 0, 0, t197 * t285, t197 * t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, 0, 0, 0, 0, 0, -t255, t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, -t226, t281 + t284, t104, -t277, -t195 * t284 - t219, t195 * t285 - t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, -qJD(4) * t106 - t239, 0, 0, 0, 0, 0, t296 - t338, -qJD(4) * t129 - t182 * t284 + t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t283, 0, 0, 0, 0, 0, t255, -t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t14 + qJD(4) * t32, qJD(2) * t15 + qJD(4) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t227, t206 * t262 - t280, t182 * t255 + t118, t175, qJD(4) * t26 + qJD(5) * t127 + t241, -qJD(4) * t24 + t182 * t286 - t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231, t226, -t281, t119, t277, t219, t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t182 * t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t9;