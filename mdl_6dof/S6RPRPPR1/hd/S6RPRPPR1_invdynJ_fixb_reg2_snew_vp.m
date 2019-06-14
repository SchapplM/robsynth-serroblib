% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 16:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:25:28
% EndTime: 2019-05-05 16:25:40
% DurationCPUTime: 6.05s
% Computational Cost: add. (32675->493), mult. (73939->729), div. (0->0), fcn. (51321->12), ass. (0->282)
t250 = sin(pkin(11));
t251 = sin(pkin(10));
t254 = cos(pkin(10));
t260 = cos(qJ(3));
t294 = qJD(1) * t260;
t258 = sin(qJ(3));
t295 = qJD(1) * t258;
t219 = t251 * t294 + t254 * t295;
t253 = cos(pkin(11));
t202 = -t253 * qJD(3) + t219 * t250;
t204 = qJD(3) * t250 + t219 * t253;
t169 = t204 * t202;
t287 = qJD(1) * qJD(3);
t277 = t260 * t287;
t286 = t258 * qJDD(1);
t225 = t277 + t286;
t242 = t260 * qJDD(1);
t278 = t258 * t287;
t269 = t242 - t278;
t195 = t225 * t251 - t254 * t269;
t332 = -t169 + t195;
t340 = t250 * t332;
t217 = t251 * t295 - t254 * t294;
t194 = t219 * t217;
t330 = qJDD(3) - t194;
t339 = t251 * t330;
t338 = t253 * t332;
t337 = t254 * t330;
t257 = sin(qJ(6));
t259 = cos(qJ(6));
t164 = t259 * t202 + t204 * t257;
t166 = -t202 * t257 + t204 * t259;
t129 = t166 * t164;
t192 = qJDD(6) + t195;
t334 = -t129 + t192;
t336 = t257 * t334;
t335 = t259 * t334;
t196 = t254 * t225 + t251 * t269;
t184 = qJDD(3) * t250 + t196 * t253;
t272 = -t253 * qJDD(3) + t196 * t250;
t123 = -t164 * qJD(6) + t259 * t184 - t257 * t272;
t213 = qJD(6) + t217;
t155 = t213 * t164;
t333 = -t155 + t123;
t182 = t217 * t202;
t145 = -t182 - t184;
t331 = -t182 + t184;
t292 = qJD(3) * t219;
t171 = t195 + t292;
t229 = qJD(3) * pkin(3) - qJ(4) * t295;
t247 = t260 ^ 2;
t252 = sin(pkin(9));
t324 = sin(qJ(1));
t325 = cos(qJ(1));
t266 = t324 * g(1) - t325 * g(2);
t265 = qJDD(1) * pkin(1) + t266;
t267 = t325 * g(1) + t324 * g(2);
t327 = qJD(1) ^ 2;
t223 = -t327 * pkin(1) - t267;
t255 = cos(pkin(9));
t301 = t255 * t223;
t264 = qJDD(1) * pkin(7) + t252 * t265 + t301;
t298 = -g(3) + qJDD(2);
t276 = t258 * t298;
t296 = qJ(4) * qJD(3);
t149 = t260 * t264 + t276 - qJD(3) * t229 + t242 * qJ(4) + (-t258 * t296 + (-pkin(2) * t260 - pkin(3) * t247) * qJD(1)) * qJD(1);
t275 = t260 * t298;
t262 = -t258 * t264 + t275 - t225 * qJ(4) + qJDD(3) * pkin(3) + (t260 * t296 + (t260 * pkin(3) + pkin(2)) * t295) * qJD(1);
t104 = -0.2e1 * qJD(4) * t217 + t254 * t149 + t251 * t262;
t273 = t257 * t184 + t259 * t272;
t96 = (qJD(6) - t213) * t166 + t273;
t162 = t164 ^ 2;
t163 = t166 ^ 2;
t329 = t202 ^ 2;
t201 = t204 ^ 2;
t212 = t213 ^ 2;
t328 = t217 ^ 2;
t216 = t219 ^ 2;
t326 = 0.2e1 * qJD(4);
t323 = pkin(4) * t251;
t220 = t255 * t265;
t271 = -t252 * t223 + t220;
t185 = -qJDD(1) * pkin(2) - t327 * pkin(7) - t271;
t244 = t247 * t327;
t157 = -t269 * pkin(3) - qJ(4) * t244 + t229 * t295 + qJDD(4) + t185;
t293 = qJD(3) * t217;
t270 = -t196 + t293;
t107 = pkin(4) * t171 + t270 * qJ(5) + t157;
t186 = pkin(4) * t217 - qJ(5) * t219;
t261 = qJD(3) ^ 2;
t84 = -pkin(4) * t261 + qJDD(3) * qJ(5) - t186 * t217 + t104;
t62 = 0.2e1 * qJD(5) * t204 - t253 * t107 + t250 * t84;
t46 = pkin(5) * t332 + pkin(8) * t145 - t62;
t177 = pkin(5) * t217 - pkin(8) * t204;
t63 = -0.2e1 * qJD(5) * t202 + t250 * t107 + t253 * t84;
t51 = -t329 * pkin(5) - t272 * pkin(8) - t217 * t177 + t63;
t24 = t257 * t51 - t259 * t46;
t25 = t257 * t46 + t259 * t51;
t12 = -t24 * t259 + t25 * t257;
t322 = t12 * t250;
t321 = t12 * t253;
t274 = t251 * t149 - t254 * t262;
t83 = qJDD(5) - t261 * qJ(5) - qJDD(3) * pkin(4) + (t326 + t186) * t219 + t274;
t320 = t250 * t83;
t319 = t253 * t83;
t70 = t272 * pkin(5) - t329 * pkin(8) + t204 * t177 + t83;
t318 = t257 * t70;
t103 = t219 * t326 + t274;
t71 = -t103 * t254 + t104 * t251;
t317 = t258 * t71;
t316 = t259 * t70;
t120 = t129 + t192;
t315 = t120 * t257;
t314 = t120 * t259;
t147 = t169 + t195;
t313 = t147 * t250;
t312 = t147 * t253;
t311 = t157 * t251;
t310 = t157 * t254;
t189 = qJDD(3) + t194;
t309 = t189 * t251;
t308 = t189 * t254;
t307 = t195 * t251;
t306 = t204 * t217;
t305 = t213 * t257;
t304 = t213 * t259;
t303 = t217 * t250;
t302 = t217 * t253;
t237 = t260 * t327 * t258;
t230 = qJDD(3) + t237;
t300 = t258 * t230;
t231 = qJDD(3) - t237;
t299 = t260 * t231;
t291 = qJD(3) * t251;
t290 = qJD(3) * t254;
t285 = t251 * t129;
t284 = t254 * t129;
t283 = t251 * t169;
t282 = t254 * t169;
t281 = -pkin(1) * t255 - pkin(2);
t280 = pkin(1) * t252 + pkin(7);
t279 = -pkin(4) * t254 - pkin(3);
t13 = t24 * t257 + t259 * t25;
t32 = t250 * t62 + t253 * t63;
t72 = t103 * t251 + t254 * t104;
t263 = -t327 * pkin(2) + t264;
t175 = t258 * t263 - t275;
t176 = t260 * t263 + t276;
t132 = t258 * t175 + t260 * t176;
t226 = t242 - 0.2e1 * t278;
t31 = t250 * t63 - t253 * t62;
t141 = t272 - t306;
t172 = -t195 + t292;
t246 = t258 ^ 2;
t243 = t246 * t327;
t236 = -t244 - t261;
t235 = -t243 - t261;
t228 = t243 + t244;
t227 = (t246 + t247) * qJDD(1);
t224 = 0.2e1 * t277 + t286;
t209 = -t216 - t261;
t208 = -t216 + t261;
t207 = t328 - t261;
t206 = -t235 * t258 - t299;
t205 = t236 * t260 - t300;
t191 = t254 * t195;
t187 = -t328 - t261;
t180 = -t201 + t328;
t179 = -t328 + t329;
t174 = t196 + t293;
t170 = -t328 - t216;
t167 = -t201 + t329;
t161 = -t201 - t328;
t160 = -t209 * t251 - t308;
t159 = t209 * t254 - t309;
t158 = -t328 - t329;
t154 = -t163 + t212;
t153 = t162 - t212;
t152 = -t201 - t329;
t151 = t187 * t254 - t339;
t150 = t187 * t251 + t337;
t140 = t272 + t306;
t136 = t184 * t253 - t204 * t303;
t135 = t202 * t302 + t250 * t272;
t134 = (-t202 * t253 + t204 * t250) * t217;
t133 = -t163 - t212;
t131 = t172 * t254 + t174 * t251;
t130 = t172 * t251 - t174 * t254;
t128 = t163 - t162;
t127 = -t159 * t258 + t160 * t260;
t126 = t179 * t253 - t313;
t125 = -t180 * t250 + t338;
t124 = -t212 - t162;
t122 = -qJD(6) * t166 - t273;
t118 = -t161 * t250 - t312;
t117 = t161 * t253 - t313;
t116 = (-t164 * t259 + t166 * t257) * t213;
t115 = (-t164 * t257 - t166 * t259) * t213;
t114 = t158 * t253 - t340;
t113 = t158 * t250 + t338;
t112 = -t150 * t258 + t151 * t260;
t111 = -t141 * t253 - t145 * t250;
t110 = -t140 * t253 - t250 * t331;
t109 = -t141 * t250 + t145 * t253;
t108 = -t162 - t163;
t101 = -t130 * t258 + t131 * t260;
t99 = t155 + t123;
t95 = (qJD(6) + t213) * t166 + t273;
t94 = t123 * t259 - t166 * t305;
t93 = t123 * t257 + t166 * t304;
t92 = -t122 * t257 + t164 * t304;
t91 = t122 * t259 + t164 * t305;
t90 = t153 * t259 - t315;
t89 = -t154 * t257 + t335;
t88 = t153 * t257 + t314;
t87 = t154 * t259 + t336;
t86 = t118 * t254 + t251 * t331;
t85 = t118 * t251 - t254 * t331;
t81 = t114 * t254 + t140 * t251;
t80 = t114 * t251 - t140 * t254;
t79 = -t133 * t257 - t314;
t78 = t133 * t259 - t315;
t77 = t111 * t254 + t152 * t251;
t76 = t111 * t251 - t152 * t254;
t75 = t124 * t259 - t336;
t74 = t124 * t257 + t335;
t73 = -t115 * t250 + t116 * t253;
t69 = -qJ(5) * t117 + t319;
t68 = -qJ(5) * t113 + t320;
t67 = t257 * t99 - t259 * t96;
t66 = -t257 * t333 - t259 * t95;
t65 = -t257 * t96 - t259 * t99;
t64 = -t257 * t95 + t259 * t333;
t60 = -t250 * t93 + t253 * t94;
t59 = -t250 * t91 + t253 * t92;
t58 = -t250 * t88 + t253 * t90;
t57 = -t250 * t87 + t253 * t89;
t56 = -t258 * t85 + t260 * t86;
t55 = -t258 * t80 + t260 * t81;
t54 = -t250 * t78 + t253 * t79;
t53 = t250 * t79 + t253 * t78;
t50 = -t250 * t74 + t253 * t75;
t49 = t250 * t75 + t253 * t74;
t48 = -pkin(4) * t117 + t63;
t47 = -pkin(4) * t113 + t62;
t44 = -pkin(8) * t78 + t316;
t43 = -pkin(8) * t74 + t318;
t42 = t251 * t333 + t254 * t54;
t41 = t251 * t54 - t254 * t333;
t40 = t260 * t72 - t317;
t39 = t251 * t95 + t254 * t50;
t38 = t251 * t50 - t254 * t95;
t37 = -pkin(5) * t333 + pkin(8) * t79 + t318;
t36 = -pkin(5) * t95 + pkin(8) * t75 - t316;
t35 = -t250 * t65 + t253 * t67;
t34 = -t250 * t64 + t253 * t66;
t33 = t250 * t67 + t253 * t65;
t30 = -qJ(5) * t109 - t31;
t29 = t108 * t251 + t254 * t35;
t28 = -t108 * t254 + t251 * t35;
t27 = t251 * t83 + t254 * t32;
t26 = t251 * t32 - t254 * t83;
t22 = -pkin(4) * t33 - pkin(5) * t65;
t21 = -t258 * t41 + t260 * t42;
t20 = -t258 * t38 + t260 * t39;
t19 = -pkin(4) * t53 - pkin(5) * t78 + t25;
t18 = -qJ(5) * t53 - t250 * t37 + t253 * t44;
t17 = -pkin(4) * t49 - pkin(5) * t74 + t24;
t16 = -qJ(5) * t49 - t250 * t36 + t253 * t43;
t15 = -t258 * t28 + t260 * t29;
t11 = -pkin(5) * t70 + pkin(8) * t13;
t10 = -pkin(8) * t65 - t12;
t9 = -pkin(5) * t108 + pkin(8) * t67 + t13;
t8 = t13 * t253 - t322;
t7 = t13 * t250 + t321;
t6 = t251 * t70 + t254 * t8;
t5 = t251 * t8 - t254 * t70;
t4 = -qJ(5) * t33 + t10 * t253 - t250 * t9;
t3 = -pkin(4) * t7 - pkin(5) * t12;
t2 = -pkin(8) * t321 - qJ(5) * t7 - t11 * t250;
t1 = -t258 * t5 + t260 * t6;
t14 = [0, 0, 0, 0, 0, qJDD(1), t266, t267, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t255 - t327 * t252) + t271, -t301 - t252 * t266 + (-0.2e1 * qJDD(1) * t252 - t327 * t255) * pkin(1), 0, pkin(1) * (t252 ^ 2 * t265 + t255 * t220), (t225 + t277) * t258, t224 * t260 + t226 * t258, t300 + t260 * (-t243 + t261), t226 * t260, t258 * (t244 - t261) + t299, 0, -t260 * t185 + pkin(2) * t226 + pkin(7) * t205 + pkin(1) * (t205 * t252 + t226 * t255), t258 * t185 - pkin(2) * t224 + pkin(7) * t206 + pkin(1) * (t206 * t252 - t224 * t255), pkin(2) * t228 + pkin(7) * t227 + pkin(1) * (t227 * t252 + t228 * t255) + t132, -pkin(2) * t185 + pkin(7) * t132 + pkin(1) * (t132 * t252 - t185 * t255), t258 * (t196 * t254 - t219 * t291) + t260 * (t196 * t251 + t219 * t290), t258 * (-t171 * t254 + t251 * t270) + t260 * (-t171 * t251 - t254 * t270), t258 * (-t208 * t251 + t337) + t260 * (t208 * t254 + t339), t258 * (t217 * t290 + t307) + t260 * (t217 * t291 - t191), t258 * (t207 * t254 - t309) + t260 * (t207 * t251 + t308), (t258 * (-t217 * t254 + t219 * t251) + t260 * (-t217 * t251 - t219 * t254)) * qJD(3), t258 * (-qJ(4) * t150 + t311) + t260 * (-pkin(3) * t171 + qJ(4) * t151 - t310) - pkin(2) * t171 + pkin(7) * t112 + pkin(1) * (t112 * t252 - t171 * t255), t258 * (-qJ(4) * t159 + t310) + t260 * (pkin(3) * t270 + qJ(4) * t160 + t311) + pkin(2) * t270 + pkin(7) * t127 + pkin(1) * (t127 * t252 + t255 * t270), t258 * (-qJ(4) * t130 - t71) + t260 * (-pkin(3) * t170 + qJ(4) * t131 + t72) - pkin(2) * t170 + pkin(7) * t101 + pkin(1) * (t101 * t252 - t170 * t255), -qJ(4) * t317 + t260 * (-pkin(3) * t157 + qJ(4) * t72) - pkin(2) * t157 + pkin(7) * t40 + pkin(1) * (-t157 * t255 + t252 * t40), t258 * (t136 * t254 + t283) + t260 * (t136 * t251 - t282), t258 * (t110 * t254 - t167 * t251) + t260 * (t110 * t251 + t167 * t254), t258 * (t125 * t254 - t145 * t251) + t260 * (t125 * t251 + t145 * t254), t258 * (t135 * t254 - t283) + t260 * (t135 * t251 + t282), t258 * (t126 * t254 - t141 * t251) + t260 * (t126 * t251 + t141 * t254), t258 * (t134 * t254 + t307) + t260 * (t134 * t251 - t191), t258 * (-qJ(4) * t80 - t251 * t47 + t254 * t68) + t260 * (-pkin(3) * t113 + qJ(4) * t81 + t251 * t68 + t254 * t47) - pkin(2) * t113 + pkin(7) * t55 + pkin(1) * (-t113 * t255 + t252 * t55), t258 * (-qJ(4) * t85 - t251 * t48 + t254 * t69) + t260 * (-pkin(3) * t117 + qJ(4) * t86 + t251 * t69 + t254 * t48) - pkin(2) * t117 + pkin(7) * t56 + pkin(1) * (-t117 * t255 + t252 * t56), t258 * (-qJ(4) * t76 + t254 * t30) + t260 * (qJ(4) * t77 + t251 * t30) + t280 * (-t258 * t76 + t260 * t77) + (t258 * t323 + t260 * t279 + t281) * t109, (t258 * (-qJ(5) * t254 + t323) + t260 * (-qJ(5) * t251 + t279) + t281) * t31 + (t280 + qJ(4)) * (-t258 * t26 + t260 * t27), t258 * (t254 * t60 + t285) + t260 * (t251 * t60 - t284), t258 * (t128 * t251 + t254 * t34) + t260 * (-t128 * t254 + t251 * t34), t258 * (t251 * t99 + t254 * t57) + t260 * (t251 * t57 - t254 * t99), t258 * (t254 * t59 - t285) + t260 * (t251 * t59 + t284), t258 * (-t251 * t96 + t254 * t58) + t260 * (t251 * t58 + t254 * t96), t258 * (t192 * t251 + t254 * t73) + t260 * (-t192 * t254 + t251 * t73), t258 * (-qJ(4) * t38 + t16 * t254 - t17 * t251) + t260 * (-pkin(3) * t49 + qJ(4) * t39 + t16 * t251 + t17 * t254) - pkin(2) * t49 + pkin(7) * t20 + pkin(1) * (t20 * t252 - t255 * t49), t258 * (-qJ(4) * t41 + t18 * t254 - t19 * t251) + t260 * (-pkin(3) * t53 + qJ(4) * t42 + t18 * t251 + t19 * t254) - pkin(2) * t53 + pkin(7) * t21 + pkin(1) * (t21 * t252 - t255 * t53), t258 * (-qJ(4) * t28 - t22 * t251 + t254 * t4) + t260 * (-pkin(3) * t33 + qJ(4) * t29 + t22 * t254 + t251 * t4) - pkin(2) * t33 + pkin(7) * t15 + pkin(1) * (t15 * t252 - t255 * t33), t258 * (-qJ(4) * t5 + t2 * t254 - t251 * t3) + t260 * (-pkin(3) * t7 + qJ(4) * t6 + t2 * t251 + t254 * t3) - pkin(2) * t7 + pkin(7) * t1 + pkin(1) * (t1 * t252 - t255 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t298, 0, 0, 0, 0, 0, 0, t230 * t260 + t236 * t258, -t231 * t258 + t235 * t260, 0, -t175 * t260 + t176 * t258, 0, 0, 0, 0, 0, 0, t150 * t260 + t151 * t258, t159 * t260 + t160 * t258, t130 * t260 + t131 * t258, t258 * t72 + t260 * t71, 0, 0, 0, 0, 0, 0, t258 * t81 + t260 * t80, t258 * t86 + t260 * t85, t258 * t77 + t260 * t76, t258 * t27 + t26 * t260, 0, 0, 0, 0, 0, 0, t258 * t39 + t260 * t38, t258 * t42 + t260 * t41, t258 * t29 + t260 * t28, t258 * t6 + t260 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, t243 - t244, t286, t237, t242, qJDD(3), -t175, -t176, 0, 0, t194, t216 - t328, t174, -t194, t172, qJDD(3), pkin(3) * t150 - t103, pkin(3) * t159 - t104, pkin(3) * t130, pkin(3) * t71, t184 * t250 + t204 * t302, -t140 * t250 + t253 * t331, t180 * t253 + t340, t202 * t303 - t253 * t272, t179 * t250 + t312, (-t202 * t250 - t204 * t253) * t217, pkin(3) * t80 - pkin(4) * t140 + qJ(5) * t114 - t319, pkin(3) * t85 - pkin(4) * t331 + qJ(5) * t118 + t320, pkin(3) * t76 - pkin(4) * t152 + qJ(5) * t111 + t32, pkin(3) * t26 - pkin(4) * t83 + qJ(5) * t32, t250 * t94 + t253 * t93, t250 * t66 + t253 * t64, t250 * t89 + t253 * t87, t250 * t92 + t253 * t91, t250 * t90 + t253 * t88, t115 * t253 + t116 * t250, pkin(3) * t38 - pkin(4) * t95 + qJ(5) * t50 + t250 * t43 + t253 * t36, pkin(3) * t41 - pkin(4) * t333 + qJ(5) * t54 + t250 * t44 + t253 * t37, pkin(3) * t28 - pkin(4) * t108 + qJ(5) * t35 + t10 * t250 + t253 * t9, pkin(3) * t5 - pkin(4) * t70 - pkin(8) * t322 + qJ(5) * t8 + t11 * t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, -t270, t170, t157, 0, 0, 0, 0, 0, 0, t113, t117, t109, t31, 0, 0, 0, 0, 0, 0, t49, t53, t33, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, t331, t152, t83, 0, 0, 0, 0, 0, 0, t95, t333, t108, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t128, t99, -t129, -t96, t192, -t24, -t25, 0, 0;];
tauJ_reg  = t14;
