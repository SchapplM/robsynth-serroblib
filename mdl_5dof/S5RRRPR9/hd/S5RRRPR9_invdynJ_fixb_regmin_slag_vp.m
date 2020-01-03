% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:41
% EndTime: 2019-12-31 21:24:52
% DurationCPUTime: 4.50s
% Computational Cost: add. (4517->414), mult. (10416->578), div. (0->0), fcn. (7519->12), ass. (0->223)
t209 = cos(qJ(2));
t264 = t209 * qJD(1);
t180 = -qJD(3) + t264;
t172 = -qJD(5) + t180;
t203 = sin(qJ(5));
t207 = cos(qJ(5));
t204 = sin(qJ(3));
t205 = sin(qJ(2));
t277 = qJD(1) * t205;
t252 = t204 * t277;
t208 = cos(qJ(3));
t266 = t208 * qJD(2);
t158 = t252 - t266;
t267 = t204 * qJD(2);
t160 = t208 * t277 + t267;
t200 = sin(pkin(9));
t201 = cos(pkin(9));
t233 = -t201 * t158 - t200 * t160;
t301 = t207 * t233;
t96 = t200 * t158 - t201 * t160;
t51 = t203 * t96 + t301;
t299 = t51 * t172;
t269 = qJD(5) * t203;
t263 = qJD(1) * qJD(2);
t250 = t209 * t263;
t262 = t205 * qJDD(1);
t271 = qJD(3) * t205;
t326 = -qJD(1) * t271 + qJDD(2);
t89 = qJD(3) * t266 + (t250 + t262) * t208 + t326 * t204;
t90 = t204 * (qJD(2) * (qJD(3) + t264) + t262) - t326 * t208;
t38 = -t200 * t89 - t201 * t90;
t39 = -t200 * t90 + t201 * t89;
t8 = qJD(5) * t301 + t203 * t38 + t207 * t39 + t269 * t96;
t339 = t8 + t299;
t325 = t203 * t233 - t207 * t96;
t338 = t325 * t51;
t337 = t325 ^ 2 - t51 ^ 2;
t195 = qJ(3) + pkin(9) + qJ(5);
t184 = cos(t195);
t183 = sin(t195);
t210 = cos(qJ(1));
t288 = t210 * t183;
t206 = sin(qJ(1));
t290 = t206 * t209;
t116 = -t184 * t290 + t288;
t287 = t210 * t184;
t118 = t206 * t183 + t209 * t287;
t166 = -t209 * pkin(2) - t205 * pkin(7) - pkin(1);
t148 = t166 * qJD(1);
t190 = pkin(6) * t264;
t170 = qJD(2) * pkin(7) + t190;
t107 = t204 * t148 + t208 * t170;
t72 = -t158 * qJ(4) + t107;
t303 = t201 * t72;
t106 = t208 * t148 - t204 * t170;
t71 = -t160 * qJ(4) + t106;
t63 = -t180 * pkin(3) + t71;
t28 = t200 * t63 + t303;
t324 = pkin(8) * t233;
t20 = t28 + t324;
t18 = t20 * t269;
t194 = t209 * qJDD(1);
t248 = t205 * t263;
t322 = -t248 + t194;
t150 = qJDD(3) - t322;
t240 = pkin(2) * t205 - pkin(7) * t209;
t162 = t240 * qJD(2);
t108 = qJD(1) * t162 + qJDD(1) * t166;
t102 = t208 * t108;
t134 = pkin(6) * t322 + qJDD(2) * pkin(7);
t16 = t150 * pkin(3) - t89 * qJ(4) - qJD(3) * t107 - t160 * qJD(4) - t204 * t134 + t102;
t270 = qJD(3) * t208;
t272 = qJD(3) * t204;
t220 = t204 * t108 + t208 * t134 + t148 * t270 - t170 * t272;
t21 = -t90 * qJ(4) - t158 * qJD(4) + t220;
t4 = t201 * t16 - t200 * t21;
t2 = t150 * pkin(4) - t39 * pkin(8) + t4;
t310 = g(3) * t205;
t169 = -qJD(2) * pkin(2) + pkin(6) * t277;
t112 = t158 * pkin(3) + qJD(4) + t169;
t57 = -pkin(4) * t233 + t112;
t336 = g(1) * t118 - g(2) * t116 + t184 * t310 - t203 * t2 - t57 * t51 + t18;
t300 = t325 * t172;
t9 = qJD(5) * t325 + t203 * t39 - t207 * t38;
t334 = -t9 - t300;
t289 = t208 * t209;
t315 = pkin(3) * t205;
t229 = -qJ(4) * t289 + t315;
t308 = qJ(4) + pkin(7);
t245 = qJD(3) * t308;
t161 = t240 * qJD(1);
t282 = pkin(6) * t252 + t208 * t161;
t333 = -qJD(1) * t229 - t204 * qJD(4) - t208 * t245 - t282;
t143 = t204 * t161;
t265 = t208 * qJD(4);
t291 = t205 * t208;
t292 = t204 * t209;
t332 = t143 + (-pkin(6) * t291 - qJ(4) * t292) * qJD(1) + t204 * t245 - t265;
t188 = pkin(6) * t262;
t135 = -qJDD(2) * pkin(2) + pkin(6) * t250 + t188;
t239 = g(1) * t210 + g(2) * t206;
t309 = g(3) * t209;
t216 = t205 * t239 - t309;
t331 = qJD(3) * pkin(7) * t180 - t135 + t216;
t115 = t183 * t290 + t287;
t117 = t206 * t184 - t209 * t288;
t5 = t200 * t16 + t201 * t21;
t3 = pkin(8) * t38 + t5;
t257 = t207 * t2 - t203 * t3;
t330 = -g(1) * t117 + g(2) * t115 + t183 * t310 - t57 * t325 + t257;
t329 = pkin(8) * t96;
t152 = t200 * t208 + t201 * t204;
t221 = t152 * t209;
t328 = qJD(1) * t221 - t152 * qJD(3);
t230 = t200 * t204 - t201 * t208;
t327 = t180 * t230;
t305 = t332 * t200 + t333 * t201;
t304 = t333 * t200 - t332 * t201;
t316 = pkin(3) * t204;
t321 = pkin(3) * t272 - t264 * t316 - t190;
t285 = t210 * t208;
t139 = t204 * t290 + t285;
t286 = t210 * t204;
t141 = t206 * t208 - t209 * t286;
t319 = -g(1) * t141 + g(2) * t139;
t317 = pkin(3) * t200;
t197 = t205 * pkin(6);
t234 = -t203 * t152 - t207 * t230;
t307 = qJD(5) * t234 + t328 * t203 + t327 * t207;
t95 = t207 * t152 - t203 * t230;
t306 = qJD(5) * t95 + t327 * t203 - t328 * t207;
t182 = pkin(6) * t289;
t283 = t208 * t162 + t267 * t197;
t44 = -t205 * t265 + t229 * qJD(2) + (-t182 + (qJ(4) * t205 - t166) * t204) * qJD(3) + t283;
t284 = t204 * t162 + t166 * t270;
t53 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t291 + (-qJD(4) * t205 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t209) * t204 + t284;
t15 = t200 * t44 + t201 * t53;
t67 = t200 * t72;
t32 = t201 * t71 - t67;
t27 = t201 * t63 - t67;
t17 = -t180 * pkin(4) + t27 + t329;
t302 = t207 * t17;
t298 = t89 * t204;
t154 = t208 * t166;
t103 = -qJ(4) * t291 + t154 + (-pkin(6) * t204 - pkin(3)) * t209;
t281 = t204 * t166 + t182;
t293 = t204 * t205;
t109 = -qJ(4) * t293 + t281;
t55 = t200 * t103 + t201 * t109;
t297 = -t328 * pkin(4) + t321;
t296 = t158 * t180;
t295 = t160 * t180;
t294 = t160 * t208;
t167 = t308 * t204;
t168 = t308 * t208;
t111 = -t200 * t167 + t201 * t168;
t279 = pkin(3) * t293 + t197;
t198 = t205 ^ 2;
t278 = -t209 ^ 2 + t198;
t276 = qJD(2) * t158;
t275 = qJD(2) * t160;
t274 = qJD(2) * t205;
t273 = qJD(2) * t209;
t268 = t169 * qJD(3);
t255 = t209 * t267;
t259 = pkin(3) * t255 + pkin(6) * t273 + t270 * t315;
t187 = t208 * pkin(3) + pkin(2);
t258 = pkin(6) + t316;
t256 = t180 * t266;
t254 = t180 * t272;
t253 = t180 * t270;
t251 = t209 * t266;
t247 = qJD(5) * t17 + t3;
t14 = -t200 * t53 + t201 * t44;
t31 = -t200 * t71 - t303;
t54 = t201 * t103 - t200 * t109;
t110 = -t201 * t167 - t200 * t168;
t244 = -qJD(3) * t148 - t134;
t83 = -pkin(8) * t230 + t111;
t242 = pkin(4) * t277 + t327 * pkin(8) + qJD(5) * t83 - t305;
t82 = -t152 * pkin(8) + t110;
t241 = t328 * pkin(8) + qJD(5) * t82 + t304;
t238 = g(1) * t206 - g(2) * t210;
t237 = t170 * t270 - t102;
t7 = t203 * t17 + t207 * t20;
t236 = -pkin(7) * t150 + t268;
t129 = t152 * t205;
t130 = t230 * t205;
t235 = -t207 * t129 + t203 * t130;
t75 = -t203 * t129 - t207 * t130;
t232 = t209 * t187 + t205 * t308;
t185 = t201 * pkin(3) + pkin(4);
t228 = t203 * t185 + t207 * t317;
t227 = t207 * t185 - t203 * t317;
t225 = pkin(1) + t232;
t224 = -0.2e1 * pkin(1) * t263 - pkin(6) * qJDD(2);
t223 = t204 * t150 - t253;
t222 = t208 * t150 + t254;
t212 = qJD(1) ^ 2;
t219 = pkin(1) * t212 + t239;
t62 = t90 * pkin(3) + qJDD(4) + t135;
t211 = qJD(2) ^ 2;
t214 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t211 + t238;
t147 = qJDD(5) + t150;
t142 = t206 * t204 + t209 * t285;
t140 = -t206 * t289 + t286;
t122 = pkin(4) * t230 - t187;
t104 = t129 * pkin(4) + t279;
t77 = t152 * t271 + t200 * t255 - t201 * t251;
t76 = -qJD(2) * t221 + t230 * t271;
t64 = t160 * pkin(3) - t96 * pkin(4);
t56 = -t76 * pkin(4) + t259;
t34 = -t129 * pkin(8) + t55;
t33 = -t209 * pkin(4) + t130 * pkin(8) + t54;
t26 = t32 + t329;
t25 = t31 - t324;
t24 = qJD(5) * t75 - t203 * t77 - t207 * t76;
t23 = qJD(5) * t235 + t203 * t76 - t207 * t77;
t22 = -t38 * pkin(4) + t62;
t11 = pkin(8) * t76 + t15;
t10 = pkin(4) * t274 + t77 * pkin(8) + t14;
t6 = -t203 * t20 + t302;
t1 = [qJDD(1), t238, t239, t198 * qJDD(1) + 0.2e1 * t209 * t248, 0.2e1 * t194 * t205 - 0.2e1 * t263 * t278, qJDD(2) * t205 + t211 * t209, qJDD(2) * t209 - t211 * t205, 0, t205 * t224 + t209 * t214, -t205 * t214 + t209 * t224, t89 * t291 + (-t204 * t271 + t251) * t160, (-t158 * t208 - t160 * t204) * t273 + (-t298 - t208 * t90 + (t158 * t204 - t294) * qJD(3)) * t205, (-t89 - t256) * t209 + (t222 + t275) * t205, (t180 * t267 + t90) * t209 + (-t223 - t276) * t205, -t150 * t209 - t180 * t274, -(-t166 * t272 + t283) * t180 + t154 * t150 - g(1) * t140 - g(2) * t142 + ((t253 + t276) * pkin(6) + (-pkin(6) * t150 + qJD(2) * t169 - t244) * t204 + t237) * t209 + (pkin(6) * t90 + t106 * qJD(2) + t135 * t204 + t208 * t268) * t205, t284 * t180 - t281 * t150 - g(1) * t139 - g(2) * t141 + (t169 * t266 + (-t254 + t275) * pkin(6) + t220) * t209 + (-t204 * t268 - t107 * qJD(2) + t135 * t208 + (t89 - t256) * pkin(6)) * t205, -t5 * t129 + t4 * t130 + t14 * t96 + t15 * t233 + t205 * t238 + t27 * t77 + t28 * t76 + t55 * t38 - t54 * t39, t5 * t55 + t28 * t15 + t4 * t54 + t27 * t14 + t62 * t279 + t112 * t259 + (-g(1) * t258 - g(2) * t225) * t210 + (g(1) * t225 - g(2) * t258) * t206, t23 * t325 + t75 * t8, t23 * t51 + t235 * t8 - t24 * t325 - t75 * t9, t75 * t147 - t23 * t172 - t8 * t209 + t274 * t325, t147 * t235 + t24 * t172 + t9 * t209 + t274 * t51, -t147 * t209 - t172 * t274, -(t207 * t10 - t203 * t11) * t172 + (-t203 * t34 + t207 * t33) * t147 - t257 * t209 + t6 * t274 - t56 * t51 + t104 * t9 - t22 * t235 + t57 * t24 - g(1) * t116 - g(2) * t118 + (-(-t203 * t33 - t207 * t34) * t172 + t7 * t209) * qJD(5), -t7 * t274 - g(1) * t115 - g(2) * t117 + t104 * t8 - t18 * t209 + t22 * t75 + t57 * t23 + t56 * t325 + ((-qJD(5) * t34 + t10) * t172 - t33 * t147 + t2 * t209) * t203 + ((qJD(5) * t33 + t11) * t172 - t34 * t147 + t247 * t209) * t207; 0, 0, 0, -t205 * t212 * t209, t278 * t212, t262, t194, qJDD(2), t205 * t219 - t188 - t309, t310 + (-pkin(6) * qJDD(1) + t219) * t209, -t180 * t294 + t298, (t89 + t296) * t208 + (-t90 + t295) * t204, (-t160 * t205 + t180 * t289) * qJD(1) + t223, (t158 * t205 - t180 * t292) * qJD(1) + t222, t180 * t277, -pkin(2) * t90 + t282 * t180 + t236 * t204 + (-t106 * t205 + (-pkin(6) * t158 - t169 * t204) * t209) * qJD(1) + t331 * t208, -pkin(2) * t89 - t143 * t180 + t236 * t208 + (-t169 * t289 + t107 * t205 + (-t160 * t209 + t180 * t291) * pkin(6)) * qJD(1) - t331 * t204, -t110 * t39 + t111 * t38 - t4 * t152 - t239 * t209 - t5 * t230 + t304 * t233 - t327 * t27 + t328 * t28 + t305 * t96 - t310, t5 * t111 + t4 * t110 - t62 * t187 - g(3) * t232 + t304 * t28 + t305 * t27 + t321 * t112 + t239 * (t187 * t205 - t209 * t308), t307 * t325 + t8 * t95, t234 * t8 - t306 * t325 + t307 * t51 - t9 * t95, t95 * t147 - t307 * t172 - t277 * t325, t147 * t234 + t306 * t172 - t277 * t51, t172 * t277, (-t203 * t83 + t207 * t82) * t147 + t122 * t9 - t22 * t234 - t6 * t277 + t306 * t57 - t297 * t51 + (t203 * t241 + t207 * t242) * t172 + t216 * t184, -(t203 * t82 + t207 * t83) * t147 + t122 * t8 + t22 * t95 + t7 * t277 + t307 * t57 + t297 * t325 + (-t203 * t242 + t207 * t241) * t172 - t216 * t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160 * t158, -t158 ^ 2 + t160 ^ 2, t89 - t296, -t295 - t90, t150, -t107 * t180 - t169 * t160 + (t244 + t310) * t204 - t237 + t319, g(1) * t142 - g(2) * t140 + g(3) * t291 - t106 * t180 + t169 * t158 - t220, (t200 * t38 - t201 * t39) * pkin(3) + (t27 - t32) * t233 + (-t28 - t31) * t96, -t27 * t31 - t28 * t32 + (g(3) * t293 - t112 * t160 + t5 * t200 + t4 * t201 + t319) * pkin(3), -t338, t337, t339, t334, t147, t227 * t147 + (-t203 * t26 + t207 * t25) * t172 + t64 * t51 + (t172 * t228 - t7) * qJD(5) + t330, -t228 * t147 - t207 * t3 - (t203 * t25 + t207 * t26) * t172 - t64 * t325 + (t172 * t227 - t302) * qJD(5) + t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t233 ^ 2 - t96 ^ 2, -t233 * t28 - t27 * t96 - t216 + t62, 0, 0, 0, 0, 0, t9 - t300, t8 - t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t338, t337, t339, t334, t147, (-qJD(5) - t172) * t7 + t330, -t6 * t172 - t207 * t247 + t336;];
tau_reg = t1;
