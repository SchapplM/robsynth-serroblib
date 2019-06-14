% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 22:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:18:01
% EndTime: 2019-05-05 22:18:15
% DurationCPUTime: 5.40s
% Computational Cost: add. (15049->412), mult. (29064->530), div. (0->0), fcn. (19002->10), ass. (0->255)
t224 = cos(qJ(4));
t220 = sin(qJ(4));
t221 = sin(qJ(3));
t269 = qJD(1) * qJD(3);
t204 = t221 * t269;
t225 = cos(qJ(3));
t267 = t225 * qJDD(1);
t185 = -t204 + t267;
t176 = -qJDD(4) + t185;
t272 = qJD(1) * t221;
t177 = -t224 * qJD(3) + t220 * t272;
t179 = t220 * qJD(3) + t224 * t272;
t283 = t179 * t177;
t232 = t176 - t283;
t289 = t232 * t220;
t175 = t179 ^ 2;
t200 = qJD(1) * t225 - qJD(4);
t307 = t200 ^ 2;
t318 = -t175 - t307;
t94 = t224 * t318 + t289;
t360 = pkin(2) * t94;
t359 = pkin(3) * t94;
t358 = pkin(8) * t94;
t288 = t232 * t224;
t96 = -t220 * t318 + t288;
t357 = pkin(8) * t96;
t216 = cos(pkin(10));
t356 = t216 * t94;
t355 = t221 * t96;
t354 = t225 * t96;
t162 = t179 * t200;
t261 = t225 * t269;
t268 = t221 * qJDD(1);
t184 = t261 + t268;
t254 = -t224 * qJDD(3) + t220 * t184;
t238 = qJD(4) * t179 + t254;
t107 = t162 + t238;
t308 = t177 ^ 2;
t155 = t308 - t307;
t353 = t225 * t107 + t221 * (t155 * t224 + t289);
t240 = -t220 * qJDD(3) - t224 * t184;
t137 = -qJD(4) * t177 - t240;
t284 = t177 * t200;
t320 = t137 + t284;
t292 = t320 * t220;
t317 = t175 - t308;
t319 = -t162 + t238;
t352 = t221 * (t319 * t224 + t292) + t225 * t317;
t215 = sin(pkin(10));
t126 = t176 + t283;
t286 = t126 * t224;
t315 = -t307 - t308;
t326 = t220 * t315 - t286;
t287 = t126 * t220;
t325 = t224 * t315 + t287;
t339 = t221 * t319 + t225 * t325;
t351 = pkin(1) * (t215 * t339 - t216 * t326) + pkin(7) * t339 - pkin(2) * t326;
t349 = pkin(3) * t326;
t348 = pkin(8) * t325;
t347 = pkin(8) * t326;
t346 = qJ(5) * t320;
t156 = -t175 + t307;
t344 = t224 * t156 - t287;
t342 = t155 * t220 - t288;
t340 = t221 * t325 - t225 * t319;
t321 = t137 - t284;
t338 = t221 * (-t156 * t220 - t286) - t225 * t321;
t316 = t175 + t308;
t337 = pkin(3) * t316;
t275 = -g(3) + qJDD(2);
t203 = t225 * t275;
t227 = qJD(1) ^ 2;
t222 = sin(qJ(1));
t226 = cos(qJ(1));
t259 = t222 * g(1) - g(2) * t226;
t180 = qJDD(1) * pkin(1) + t259;
t248 = g(1) * t226 + g(2) * t222;
t181 = -pkin(1) * t227 - t248;
t273 = t215 * t180 + t216 * t181;
t133 = -pkin(2) * t227 + qJDD(1) * pkin(7) + t273;
t251 = -pkin(3) * t225 - pkin(8) * t221;
t252 = t227 * t251 + t133;
t306 = qJD(3) ^ 2;
t99 = -qJDD(3) * pkin(3) - t306 * pkin(8) + t252 * t221 - t203;
t231 = t238 * pkin(4) - t346 + t99;
t219 = sin(qJ(6));
t223 = cos(qJ(6));
t142 = -t223 * t177 + t179 * t219;
t144 = t177 * t219 + t179 * t223;
t102 = t144 * t142;
t170 = qJDD(6) + t176;
t322 = -t102 + t170;
t336 = t219 * t322;
t334 = t221 * t316;
t332 = t223 * t322;
t330 = t225 * t316;
t324 = -t319 * t220 + t224 * t320;
t194 = qJD(6) + t200;
t121 = t194 * t142;
t78 = -t142 * qJD(6) + t223 * t137 + t219 * t238;
t323 = -t121 + t78;
t145 = pkin(4) * t177 - qJ(5) * t179;
t257 = t221 * t275;
t100 = -t306 * pkin(3) + qJDD(3) * pkin(8) + t225 * t252 + t257;
t255 = t216 * t180 - t215 * t181;
t132 = -qJDD(1) * pkin(2) - t227 * pkin(7) - t255;
t242 = -t185 + t204;
t243 = t184 + t261;
t92 = pkin(3) * t242 - pkin(8) * t243 + t132;
t56 = t224 * t100 + t220 * t92;
t253 = -t176 * qJ(5) - t177 * t145 + t56;
t313 = -pkin(4) * (t318 + t307) - qJ(5) * t232 + t253;
t305 = pkin(4) + pkin(5);
t55 = t220 * t100 - t224 * t92;
t44 = t176 * pkin(4) - qJ(5) * t307 + t179 * t145 + qJDD(5) + t55;
t28 = t126 * pkin(5) - t321 * pkin(9) + t44;
t249 = pkin(5) * t200 - pkin(9) * t179;
t271 = qJD(5) * t200;
t191 = -0.2e1 * t271;
t245 = t191 + t253;
t43 = -pkin(4) * t307 + t245;
t35 = -pkin(5) * t308 + pkin(9) * t238 - t200 * t249 + t43;
t14 = t219 * t35 - t223 * t28;
t15 = t219 * t28 + t223 * t35;
t8 = -t14 * t223 + t15 * t219;
t9 = t219 * t14 + t223 * t15;
t312 = qJ(5) * t9 - t305 * t8;
t141 = t144 ^ 2;
t192 = t194 ^ 2;
t115 = -t141 - t192;
t86 = t102 + t170;
t300 = t219 * t86;
t65 = t115 * t223 - t300;
t297 = t223 * t86;
t66 = -t115 * t219 - t297;
t311 = qJ(5) * t66 - t305 * t65 + t15;
t140 = t142 ^ 2;
t93 = -t192 - t140;
t48 = t219 * t93 + t332;
t49 = t223 * t93 - t336;
t310 = qJ(5) * t49 - t305 * t48 + t14;
t256 = t219 * t137 - t223 * t238;
t58 = (qJD(6) - t194) * t144 + t256;
t61 = t121 + t78;
t32 = -t219 * t58 - t223 * t61;
t34 = t219 * t61 - t223 * t58;
t309 = qJ(5) * t34 - t305 * t32;
t304 = pkin(4) * t224;
t258 = -pkin(4) * t200 - 0.2e1 * qJD(5);
t36 = t254 * pkin(5) + t308 * pkin(9) + t231 + (pkin(5) * qJD(4) - t249 + t258) * t179;
t301 = t219 * t36;
t299 = t220 * t99;
t298 = t223 * t36;
t296 = t224 * t99;
t295 = qJ(5) * t224;
t291 = t321 * t220;
t290 = t321 * t224;
t282 = t194 * t219;
t281 = t194 * t223;
t280 = t200 * t220;
t279 = t200 * t224;
t199 = t225 * t227 * t221;
t190 = qJDD(3) + t199;
t278 = t221 * t190;
t189 = -t199 + qJDD(3);
t276 = t225 * t189;
t266 = t177 * t279;
t265 = t225 * t102;
t264 = t225 * t283;
t262 = pkin(1) * t215 + pkin(7);
t260 = -qJ(5) * t220 - pkin(3);
t30 = t220 * t55 + t224 * t56;
t116 = t221 * t133 - t203;
t117 = t225 * t133 + t257;
t79 = t221 * t116 + t225 * t117;
t250 = -pkin(4) * t44 + qJ(5) * t43;
t154 = t179 * t280;
t247 = t221 * (t137 * t224 + t154) - t264;
t246 = -t177 * t280 - t224 * t238;
t244 = -pkin(4) * t321 - qJ(5) * t107;
t241 = t220 * t56 - t224 * t55;
t236 = (t177 * t220 + t179 * t224) * t200;
t234 = -pkin(1) * t216 - pkin(2) + t251;
t233 = t221 * (-t154 + t266) + t225 * t176;
t230 = 0.2e1 * qJD(5) * t179 - t231;
t229 = -pkin(4) * t126 + qJ(5) * t315 - t44;
t228 = t221 * (t220 * t238 - t266) + t264;
t212 = t225 ^ 2;
t211 = t221 ^ 2;
t209 = t212 * t227;
t207 = t211 * t227;
t196 = -t209 - t306;
t195 = -t207 - t306;
t188 = t207 + t209;
t187 = (t211 + t212) * qJDD(1);
t186 = -0.2e1 * t204 + t267;
t183 = 0.2e1 * t261 + t268;
t153 = -t195 * t221 - t276;
t152 = t196 * t225 - t278;
t119 = -t141 + t192;
t118 = t140 - t192;
t113 = (qJD(4) - t200) * t177 + t240;
t108 = (-qJD(4) - t200) * t179 - t254;
t104 = t137 * t220 - t179 * t279;
t101 = t141 - t140;
t82 = (-t142 * t223 + t144 * t219) * t194;
t81 = (t142 * t219 + t144 * t223) * t194;
t80 = -t140 - t141;
t77 = -qJD(6) * t144 - t256;
t76 = t108 * t224 + t291;
t75 = -t107 * t224 + t291;
t73 = -t107 * t220 - t290;
t72 = t118 * t223 - t300;
t71 = -t119 * t219 + t332;
t70 = -t118 * t219 - t297;
t69 = -t119 * t223 - t336;
t67 = -t113 * t221 + t354;
t63 = -t221 * t320 - t354;
t57 = (qJD(6) + t194) * t144 + t256;
t53 = -t144 * t282 + t223 * t78;
t52 = -t144 * t281 - t219 * t78;
t51 = t142 * t281 - t219 * t77;
t50 = -t142 * t282 - t223 * t77;
t46 = t225 * t75 - t334;
t45 = t179 * t258 + t231;
t42 = (-t319 + t162) * pkin(4) + t230;
t41 = pkin(4) * t162 + t230 + t346;
t40 = qJ(5) * t316 + t44;
t39 = (t316 - t307) * pkin(4) + t245;
t38 = t220 * t65 + t224 * t66;
t37 = t220 * t66 - t224 * t65;
t33 = -t219 * t323 - t223 * t57;
t31 = t219 * t57 - t223 * t323;
t26 = t220 * t48 + t224 * t49;
t25 = t220 * t49 - t224 * t48;
t23 = -t221 * t323 + t225 * t38;
t22 = t220 * t44 + t224 * t43;
t21 = t220 * t43 - t224 * t44;
t20 = -t221 * t57 + t225 * t26;
t19 = -pkin(9) * t65 + qJ(5) * t323 - t298;
t18 = -pkin(9) * t48 + qJ(5) * t57 - t301;
t17 = t220 * t32 + t224 * t34;
t16 = t220 * t34 - t224 * t32;
t13 = t22 * t225 + t221 * t45;
t12 = t17 * t225 - t221 * t80;
t11 = -pkin(9) * t66 + t305 * t323 + t301;
t10 = -pkin(9) * t49 + t305 * t57 - t298;
t7 = -pkin(9) * t8 - qJ(5) * t36;
t6 = -pkin(9) * t32 + qJ(5) * t80 - t8;
t5 = -pkin(9) * t34 + t305 * t80 - t9;
t4 = -pkin(9) * t9 - t305 * t36;
t3 = t220 * t8 + t224 * t9;
t2 = t220 * t9 - t224 * t8;
t1 = t221 * t36 + t225 * t3;
t24 = [0, 0, 0, 0, 0, qJDD(1), t259, t248, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t216 - t215 * t227) + t255, pkin(1) * (-qJDD(1) * t215 - t216 * t227) - t273, 0, pkin(1) * (t215 * t273 + t216 * t255), t243 * t221, t183 * t225 + t186 * t221, t278 + t225 * (-t207 + t306), -t242 * t225, t221 * (t209 - t306) + t276, 0, -t225 * t132 + pkin(2) * t186 + pkin(7) * t152 + pkin(1) * (t152 * t215 + t186 * t216), t221 * t132 - pkin(2) * t183 + pkin(7) * t153 + pkin(1) * (t153 * t215 - t183 * t216), pkin(2) * t188 + pkin(7) * t187 + pkin(1) * (t187 * t215 + t188 * t216) + t79, -pkin(2) * t132 + pkin(7) * t79 + pkin(1) * (-t132 * t216 + t215 * t79), t247, -t352, t338, t228, t353, t233, t221 * (t299 - t347) + t225 * (t55 - t349) + t351, t221 * (t296 - t358) + t225 * (t56 - t359) - t360 + pkin(7) * t67 + pkin(1) * (t215 * t67 - t356), -t221 * t241 + t262 * (t225 * t76 - t334) + t234 * (t108 * t220 - t290), t262 * (t221 * t99 + t225 * t30) + t234 * t241, t247, t338, t352, t233, -t353, t228, t221 * (-t220 * t42 - t295 * t319 - t347) + t225 * (-t229 - t349) + t351, t221 * (-pkin(8) * t73 - t220 * t39 + t224 * t40) + t225 * (-pkin(3) * t73 - t244) - pkin(2) * t73 + pkin(7) * t46 + pkin(1) * (t215 * t46 - t216 * t73), t221 * (-pkin(4) * t292 + t224 * t41 + t358) + t225 * (0.2e1 * t271 - t313 + t359) + t360 + pkin(7) * t63 + pkin(1) * (t215 * t63 + t356), t221 * (-pkin(8) * t21 + (pkin(4) * t220 - t295) * t45) + t225 * (-pkin(3) * t21 - t250) - pkin(2) * t21 + pkin(7) * t13 + pkin(1) * (t13 * t215 - t21 * t216), t221 * (-t220 * t52 + t224 * t53) + t265, t221 * (-t220 * t31 + t224 * t33) + t225 * t101, t221 * (-t220 * t69 + t224 * t71) + t225 * t61, t221 * (-t220 * t50 + t224 * t51) - t265, t221 * (-t220 * t70 + t224 * t72) - t225 * t58, t221 * (-t220 * t81 + t224 * t82) + t225 * t170, t221 * (-pkin(8) * t25 - t10 * t220 + t18 * t224) + t225 * (-pkin(3) * t25 - t310) - pkin(2) * t25 + pkin(7) * t20 + pkin(1) * (t20 * t215 - t216 * t25), t221 * (-pkin(8) * t37 - t11 * t220 + t19 * t224) + t225 * (-pkin(3) * t37 - t311) - pkin(2) * t37 + pkin(7) * t23 + pkin(1) * (t215 * t23 - t216 * t37), t221 * (-pkin(8) * t16 - t220 * t5 + t224 * t6) + t225 * (-pkin(3) * t16 - t309) - pkin(2) * t16 + pkin(7) * t12 + pkin(1) * (t12 * t215 - t16 * t216), t221 * (-pkin(8) * t2 - t220 * t4 + t224 * t7) + t225 * (-pkin(3) * t2 - t312) - pkin(2) * t2 + pkin(7) * t1 + pkin(1) * (t1 * t215 - t2 * t216); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, 0, 0, 0, 0, 0, 0, t190 * t225 + t196 * t221, -t189 * t221 + t195 * t225, 0, -t116 * t225 + t117 * t221, 0, 0, 0, 0, 0, 0, t340, t113 * t225 + t355, t221 * t76 + t330, t221 * t30 - t225 * t99, 0, 0, 0, 0, 0, 0, t340, t221 * t75 + t330, t225 * t320 - t355, t22 * t221 - t225 * t45, 0, 0, 0, 0, 0, 0, t221 * t26 + t225 * t57, t221 * t38 + t225 * t323, t17 * t221 + t225 * t80, t221 * t3 - t225 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, t207 - t209, t268, t199, t267, qJDD(3), -t116, -t117, 0, 0, t104, t324, t344, t246, t342, t236, -pkin(3) * t319 - t296 + t348, pkin(3) * t113 + t299 + t357, pkin(8) * t76 + t30 + t337, -pkin(3) * t99 + pkin(8) * t30, t104, t344, -t324, t236, -t342, t246, t224 * t42 + t260 * t319 + t348, pkin(8) * t75 + t220 * t40 + t224 * t39 + t337, -t357 + t220 * t41 + (pkin(3) + t304) * t320, pkin(8) * t22 + (t260 - t304) * t45, t220 * t53 + t224 * t52, t220 * t33 + t224 * t31, t220 * t71 + t224 * t69, t220 * t51 + t224 * t50, t220 * t72 + t224 * t70, t220 * t82 + t224 * t81, pkin(3) * t57 + pkin(8) * t26 + t10 * t224 + t18 * t220, pkin(3) * t323 + pkin(8) * t38 + t11 * t224 + t19 * t220, pkin(3) * t80 + pkin(8) * t17 + t220 * t6 + t224 * t5, -pkin(3) * t36 + pkin(8) * t3 + t220 * t7 + t224 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, t317, t321, -t283, -t107, -t176, -t55, -t56, 0, 0, t283, t321, -t317, -t176, t107, -t283, t229, t244, t191 + t313, t250, -t102, -t101, -t61, t102, t58, -t170, t310, t311, t309, t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, t321, t318, t44, 0, 0, 0, 0, 0, 0, t48, t65, t32, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t101, t61, -t102, -t58, t170, -t14, -t15, 0, 0;];
tauJ_reg  = t24;
