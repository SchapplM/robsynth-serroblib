% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:11
% EndTime: 2021-01-16 00:50:32
% DurationCPUTime: 6.30s
% Computational Cost: add. (4512->454), mult. (11335->643), div. (0->0), fcn. (10404->14), ass. (0->238)
t183 = sin(pkin(6));
t184 = cos(pkin(12));
t186 = cos(pkin(7));
t277 = t184 * t186;
t182 = sin(pkin(7));
t187 = cos(pkin(6));
t281 = t182 * t187;
t124 = t183 * t277 + t281;
t191 = sin(qJ(3));
t180 = sin(pkin(12));
t194 = cos(qJ(3));
t282 = t180 * t194;
t103 = t124 * t191 + t183 * t282;
t190 = sin(qJ(4));
t276 = t187 * t186;
t278 = t183 * t182;
t122 = t184 * t278 - t276;
t193 = cos(qJ(4));
t294 = t122 * t193;
t344 = t103 * t190 + t294;
t123 = t184 * t276 - t278;
t181 = sin(pkin(11));
t185 = cos(pkin(11));
t285 = t180 * t186;
t102 = -t123 * t185 + t181 * t285;
t284 = t180 * t187;
t125 = t181 * t184 + t185 * t284;
t75 = t102 * t191 - t125 * t194;
t121 = t183 * t186 + t184 * t281;
t286 = t180 * t182;
t98 = t121 * t185 - t181 * t286;
t38 = t190 * t75 - t193 * t98;
t126 = -t181 * t284 + t184 * t185;
t99 = t123 * t181 + t185 * t285;
t69 = -t126 * t194 + t191 * t99;
t97 = t121 * t181 + t185 * t286;
t40 = t190 * t69 + t193 * t97;
t345 = g(1) * t40 + g(2) * t38 - g(3) * t344;
t161 = t187 * qJDD(1) + qJDD(2);
t250 = qJDD(1) * t183;
t233 = t184 * t250;
t114 = t161 * t186 - t182 * t233;
t164 = qJD(1) * t187 + qJD(2);
t266 = qJD(1) * t183;
t243 = t184 * t266;
t116 = t164 * t186 - t182 * t243;
t257 = qJD(4) * t193;
t258 = qJD(4) * t190;
t216 = t191 * t277 + t282;
t246 = t194 * t277;
t283 = t180 * t191;
t217 = t246 - t283;
t262 = qJD(3) * t194;
t49 = qJDD(3) * pkin(9) + (t161 * t191 + t164 * t262) * t182 + (qJD(1) * qJD(3) * t217 + qJDD(1) * t216) * t183;
t205 = t216 * t183;
t280 = t182 * t191;
t87 = qJD(1) * t205 + t164 * t280;
t85 = qJD(3) * pkin(9) + t87;
t224 = -t193 * t114 + t116 * t258 + t190 * t49 + t85 * t257;
t12 = -qJDD(4) * pkin(4) + t224;
t189 = sin(qJ(5));
t192 = cos(qJ(5));
t249 = t190 * qJDD(3);
t263 = qJD(3) * t193;
t255 = qJD(5) * t190;
t346 = qJD(3) * t255 - qJDD(4);
t89 = (qJD(4) * (qJD(5) + t263) + t249) * t189 + t346 * t192;
t7 = t89 * pkin(5) + qJDD(6) + t12;
t355 = -t345 - t7;
t195 = qJD(4) ^ 2;
t247 = t183 * t283;
t106 = -t124 * t194 + t247;
t70 = t126 * t191 + t194 * t99;
t76 = t102 * t194 + t125 * t191;
t338 = g(1) * t70 + g(2) * t76;
t209 = g(3) * t106 + t338;
t227 = t186 * t243;
t264 = qJD(3) * t191;
t242 = t182 * t264;
t244 = t180 * t266;
t198 = -t194 * (t161 * t182 + t186 * t233) + qJDD(1) * t247 + t164 * t242 + t227 * t264 + t244 * t262;
t347 = qJD(3) * t87 - t198;
t354 = -0.2e1 * qJDD(3) * pkin(3) + pkin(9) * t195 - t209 - t347;
t170 = pkin(5) * t192 + pkin(4);
t323 = qJ(6) + pkin(10);
t353 = -t170 * t193 - t323 * t190;
t37 = -t97 * t190 + t193 * t69;
t42 = t190 * t98 + t193 * t75;
t226 = pkin(4) * t190 - pkin(10) * t193;
t151 = t226 * qJD(4);
t349 = t151 - t87;
t295 = t122 * t190;
t78 = t103 * t193 - t295;
t274 = t189 * t193;
t325 = pkin(9) * t189;
t86 = -t191 * t244 + t194 * (t164 * t182 + t227);
t343 = t349 * t192 + t258 * t325 + t86 * t274;
t155 = -pkin(4) * t193 - pkin(10) * t190 - pkin(3);
t254 = qJD(5) * t192;
t271 = t192 * t193;
t342 = t155 * t254 + t349 * t189 - t86 * t271;
t341 = t116 * t193 - t190 * t85;
t253 = t192 * qJD(4);
t265 = qJD(3) * t190;
t140 = t189 * t265 - t253;
t165 = -qJD(5) + t263;
t291 = t140 * t165;
t251 = qJD(3) * qJD(4);
t236 = t193 * t251;
t88 = -qJD(5) * t253 + (-t236 - t249) * t192 + t346 * t189;
t340 = -t88 + t291;
t259 = qJD(4) * t189;
t142 = t192 * t265 + t259;
t290 = t142 * t165;
t339 = t89 - t290;
t96 = t187 * t280 + t205;
t68 = t193 * t96 - t295;
t279 = t182 * t194;
t245 = t187 * t279;
t95 = -t183 * t246 - t245 + t247;
t92 = t95 * t192;
t35 = -t189 * t68 + t92;
t336 = -g(2) * (t189 * t42 + t192 * t76) - g(1) * (t189 * t37 + t192 * t70) - g(3) * (-t189 * t78 + t92);
t335 = t142 ^ 2;
t174 = t193 * qJDD(3);
t137 = t190 * t251 + qJDD(5) - t174;
t327 = pkin(5) * t137;
t326 = pkin(5) * t140;
t324 = t189 * t7;
t63 = t190 * t116 + t193 * t85;
t58 = qJD(4) * pkin(10) + t63;
t80 = qJD(3) * t155 - t86;
t20 = -t189 * t58 + t192 * t80;
t18 = -qJ(6) * t142 + t20;
t13 = -pkin(5) * t165 + t18;
t322 = -t18 + t13;
t166 = pkin(9) * t271;
t221 = pkin(5) * t190 - qJ(6) * t271;
t252 = t192 * qJD(6);
t299 = qJ(6) * t190;
t321 = -t190 * t252 + t221 * qJD(4) + (-t166 + (-t155 + t299) * t189) * qJD(5) + t343;
t272 = t190 * t192;
t320 = (-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t272 + (-qJD(6) * t190 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t193) * t189 + t342;
t319 = qJ(6) * t88;
t318 = qJ(6) * t89;
t317 = qJD(3) * pkin(3);
t316 = t12 * t189;
t315 = t140 * t86;
t21 = t189 * t80 + t192 * t58;
t19 = -qJ(6) * t140 + t21;
t314 = t165 * t19;
t312 = t189 * t69;
t311 = t189 * t75;
t309 = t189 * t95;
t308 = t189 * t96;
t304 = t86 * t142;
t303 = t88 * t189;
t232 = qJD(5) * t323;
t240 = t189 * t263;
t150 = t226 * qJD(3);
t300 = -t189 * t150 - t192 * t341;
t302 = qJ(6) * t240 - t189 * t232 + t252 + t300;
t131 = t192 * t150;
t301 = -qJD(3) * t221 - t192 * t232 - t131 + (-qJD(6) + t341) * t189;
t289 = t142 * t189;
t288 = t142 * t192;
t273 = t190 * t114;
t268 = t189 * t155 + t166;
t178 = t190 ^ 2;
t267 = -t193 ^ 2 + t178;
t261 = qJD(4) * t140;
t260 = qJD(4) * t142;
t256 = qJD(5) * t189;
t248 = t345 * t189;
t241 = t182 * t262;
t239 = t165 * t253;
t238 = t165 * t256;
t237 = t165 * t254;
t235 = t194 * t251;
t11 = qJDD(4) * pkin(10) + qJD(4) * t341 + t193 * t49 + t273;
t28 = qJD(3) * t151 + qJDD(3) * t155 + t198;
t231 = t192 * t11 + t189 * t28 + t80 * t254 - t58 * t256;
t230 = -qJD(6) - t326;
t36 = t192 * t68 + t309;
t67 = t190 * t96 + t294;
t196 = qJD(3) ^ 2;
t223 = qJDD(3) * t194 - t191 * t196;
t147 = pkin(3) * t277 + pkin(9) * t180;
t222 = pkin(3) * t278 - t147 * t187;
t146 = -t180 * pkin(3) + pkin(9) * t277;
t220 = pkin(9) * t278 - t146 * t187;
t57 = -qJD(4) * pkin(4) - t341;
t128 = t186 * t190 + t193 * t280;
t109 = -t128 * t189 - t192 * t279;
t218 = -t128 * t192 + t189 * t279;
t127 = -t186 * t193 + t190 * t280;
t215 = t189 * t137 - t237;
t214 = t192 * t137 + t238;
t213 = -g(1) * (-t192 * t69 + t274 * t70) - g(2) * (-t192 * t75 + t274 * t76) - g(3) * (t106 * t274 + t192 * t96);
t212 = -g(1) * (-t271 * t70 - t312) - g(2) * (-t271 * t76 - t311) - g(3) * (-t106 * t271 + t308);
t211 = -g(1) * t37 - g(2) * t42 + g(3) * t78;
t208 = -g(3) * t187 + (-g(1) * t181 + g(2) * t185) * t183;
t204 = -pkin(10) * t137 - t165 * t57;
t84 = -t86 - t317;
t202 = -pkin(9) * qJDD(4) + (t84 + t86 - t317) * qJD(4);
t201 = -g(1) * (-t189 * t70 + t192 * t37) - g(2) * (-t189 * t76 + t192 * t42) - g(3) * (-t192 * t78 - t309) - t231;
t27 = t192 * t28;
t199 = -qJD(5) * t21 - t189 * t11 + t27;
t197 = t199 + t336;
t157 = t323 * t192;
t156 = t323 * t189;
t152 = (pkin(5) * t189 + pkin(9)) * t190;
t145 = pkin(3) * t285 - pkin(9) * t184;
t144 = pkin(3) * t184 + pkin(9) * t285;
t139 = t192 * t155;
t136 = t140 ^ 2;
t115 = pkin(9) * t257 + (t189 * t257 + t190 * t254) * pkin(5);
t111 = -t189 * t299 + t268;
t108 = qJD(4) * t128 + t190 * t241;
t107 = -qJD(4) * t127 + t193 * t241;
t94 = -qJ(6) * t272 + t139 + (-pkin(5) - t325) * t193;
t91 = t96 * qJD(3);
t90 = (t183 * t217 + t245) * qJD(3);
t56 = qJD(5) * t218 - t189 * t107 + t192 * t242;
t55 = qJD(5) * t109 + t192 * t107 + t189 * t242;
t51 = pkin(5) * t240 + t63;
t45 = -t230 + t57;
t34 = qJD(4) * t68 + t190 * t90;
t33 = -qJD(4) * t67 + t193 * t90;
t9 = t108 * t140 + t109 * t137 + t127 * t89 - t165 * t56;
t8 = t108 * t142 - t127 * t88 + t137 * t218 + t165 * t55;
t6 = -qJD(5) * t36 - t189 * t33 + t91 * t192;
t5 = t35 * qJD(5) + t91 * t189 + t192 * t33;
t4 = -qJD(6) * t140 + t231 - t318;
t3 = -qJD(6) * t142 + t199 + t319 + t327;
t2 = t137 * t35 + t140 * t34 - t165 * t6 + t67 * t89;
t1 = -t137 * t36 + t142 * t34 + t165 * t5 - t67 * t88;
t10 = [qJDD(1) - g(3), t161 * t187 - g(3) + (t180 ^ 2 + t184 ^ 2) * t183 ^ 2 * qJDD(1), 0, -qJD(3) * t91 - qJDD(3) * t95, -qJD(3) * t90 - qJDD(3) * t96, 0, 0, 0, 0, 0, -t95 * t174 - t34 * qJD(4) - t67 * qJDD(4) + (-t193 * t91 + t95 * t258) * qJD(3), t95 * t249 - t33 * qJD(4) - t68 * qJDD(4) + (t190 * t91 + t95 * t257) * qJD(3), 0, 0, 0, 0, 0, t2, t1, t2, t1, -t140 * t5 - t142 * t6 + t35 * t88 - t36 * t89, t13 * t6 + t19 * t5 + t3 * t35 + t34 * t45 + t36 * t4 + t67 * t7 - g(3); 0, t208 + t161, 0, t223 * t182, (-qJDD(3) * t191 - t194 * t196) * t182, 0, 0, 0, 0, 0, -t108 * qJD(4) - t127 * qJDD(4) + (-t190 * t235 + t193 * t223) * t182, -t107 * qJD(4) - t128 * qJDD(4) + (-t190 * t223 - t193 * t235) * t182, 0, 0, 0, 0, 0, t9, t8, t9, t8, t109 * t88 - t140 * t55 - t142 * t56 + t218 * t89, t108 * t45 + t109 * t3 + t127 * t7 + t13 * t56 + t19 * t55 - t218 * t4 + t208; 0, 0, qJDD(3), g(3) * t95 + t338 + t347, -t161 * t280 - g(1) * t69 - g(2) * t75 + g(3) * t96 - t216 * t250 + (-t164 * t279 - t217 * t266 + t86) * qJD(3), qJDD(3) * t178 + 0.2e1 * t190 * t236, 0.2e1 * t190 * t174 - 0.2e1 * t267 * t251, qJDD(4) * t190 + t193 * t195, qJDD(4) * t193 - t190 * t195, 0, t202 * t190 - t354 * t193, t354 * t190 + t202 * t193, -t88 * t272 + (-t189 * t255 + t193 * t253) * t142, (-t140 * t192 - t289) * t257 + (t303 - t192 * t89 + (t140 * t189 - t288) * qJD(5)) * t190, (t88 - t239) * t193 + (t214 + t260) * t190, (t165 * t259 + t89) * t193 + (-t215 - t261) * t190, -t137 * t193 - t165 * t258, t139 * t137 + (t155 * t256 - t343) * t165 + (t58 * t254 - t27 + (t237 + t261) * pkin(9) + (-pkin(9) * t137 + qJD(4) * t57 + qJD(5) * t80 + t11) * t189) * t193 + (pkin(9) * t89 + qJD(4) * t20 + t254 * t57 - t315 + t316) * t190 + t212, -t268 * t137 + t342 * t165 + (t57 * t253 + (-t238 + t260) * pkin(9) + t231) * t193 + (-t57 * t256 - t21 * qJD(4) + t12 * t192 - t304 + (-t88 - t239) * pkin(9)) * t190 + t213, t115 * t140 + t137 * t94 + t152 * t89 + (t45 * t259 - t3) * t193 - t321 * t165 + (qJD(4) * t13 + t254 * t45 - t315 + t324) * t190 + t212, -t111 * t137 + t115 * t142 - t152 * t88 + (t253 * t45 + t4) * t193 + t320 * t165 + (-qJD(4) * t19 + t192 * t7 - t45 * t256 - t304) * t190 + t213, -t111 * t89 + t88 * t94 - t321 * t142 - t320 * t140 + (-t13 * t192 - t189 * t19) * t257 + (-t189 * t4 - t192 * t3 + (t13 * t189 - t19 * t192) * qJD(5) + t209) * t190, t4 * t111 + t3 * t94 + t7 * t152 - g(1) * (-pkin(5) * t312 - (t185 * t144 - t181 * t220) * t191 + (-t185 * t145 + t181 * t222) * t194 + t353 * t70) - g(2) * (-pkin(5) * t311 - (t181 * t144 + t185 * t220) * t191 + (-t181 * t145 - t185 * t222) * t194 + t353 * t76) - g(3) * (pkin(5) * t308 - (-pkin(9) * t281 - t146 * t183) * t191 + (pkin(3) * t281 + t147 * t183) * t194 + t353 * t106) + (-t190 * t86 + t115) * t45 + t320 * t19 + t321 * t13; 0, 0, 0, 0, 0, -t190 * t196 * t193, t267 * t196, t249, t174, qJDD(4), qJD(4) * t63 - t84 * t265 - t224 - t345, -t273 + (-qJD(3) * t84 - t49) * t193 + t211, -t165 * t288 - t303, -t339 * t189 + t340 * t192, (-t142 * t190 + t165 * t271) * qJD(3) + t215, (t140 * t190 - t165 * t274) * qJD(3) + t214, t165 * t265, -t20 * t265 - pkin(4) * t89 + t131 * t165 - t63 * t140 + (-t165 * t341 + t204) * t189 + (pkin(10) * qJD(5) * t165 - t12 - t345) * t192, t21 * t265 + pkin(4) * t88 + t316 - t63 * t142 + (-pkin(10) * t256 + t300) * t165 + t204 * t192 + t248, -t13 * t265 - t137 * t156 - t140 * t51 - t170 * t89 - t301 * t165 + (-t45 * t263 + (t45 + t326) * qJD(5)) * t189 + t355 * t192, -t137 * t157 - t142 * t51 + t170 * t88 + t324 + t302 * t165 + (pkin(5) * t289 + t192 * t45) * qJD(5) + (t19 * t190 - t45 * t271) * qJD(3) + t248, -t156 * t88 - t157 * t89 - t301 * t142 - t302 * t140 + (t13 * t165 + t4) * t192 + (-t3 + t314) * t189 - t211, t4 * t157 - t3 * t156 - t7 * t170 - g(1) * (t170 * t40 - t323 * t37) - g(2) * (t170 * t38 - t323 * t42) - g(3) * (-t170 * t344 + t323 * t78) + (pkin(5) * t256 - t51) * t45 + t302 * t19 + t301 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142 * t140, -t136 + t335, -t88 - t291, -t290 - t89, t137, -t142 * t57 - t165 * t21 + t197, t140 * t57 - t165 * t20 + t201, 0.2e1 * t327 + t319 - t314 + (t230 - t45) * t142 + t197, -pkin(5) * t335 + t318 - t165 * t18 + (qJD(6) + t45) * t140 + t201, t88 * pkin(5) - t322 * t140, t322 * t19 + (-t45 * t142 + t3 + t336) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t339, t340, -t136 - t335, t13 * t142 + t19 * t140 - t355;];
tau_reg = t10;
