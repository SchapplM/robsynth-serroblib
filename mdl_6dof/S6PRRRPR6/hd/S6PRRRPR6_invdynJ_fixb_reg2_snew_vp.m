% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRPR6
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
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 08:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRPR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:28:13
% EndTime: 2019-05-05 08:28:29
% DurationCPUTime: 5.76s
% Computational Cost: add. (15522->415), mult. (30047->552), div. (0->0), fcn. (21306->12), ass. (0->260)
t223 = cos(qJ(4));
t219 = sin(qJ(4));
t220 = sin(qJ(3));
t268 = qJD(2) * qJD(3);
t204 = t220 * t268;
t224 = cos(qJ(3));
t266 = t224 * qJDD(2);
t184 = -t204 + t266;
t177 = -qJDD(4) + t184;
t271 = qJD(2) * t220;
t178 = -t223 * qJD(3) + t219 * t271;
t180 = t219 * qJD(3) + t223 * t271;
t281 = t180 * t178;
t234 = t177 - t281;
t287 = t234 * t219;
t176 = t180 ^ 2;
t201 = qJD(2) * t224 - qJD(4);
t307 = t201 ^ 2;
t318 = -t176 - t307;
t98 = t223 * t318 + t287;
t362 = pkin(2) * t98;
t361 = pkin(3) * t98;
t360 = pkin(9) * t98;
t286 = t234 * t223;
t100 = -t219 * t318 + t286;
t359 = pkin(9) * t100;
t225 = cos(qJ(2));
t358 = t225 * t98;
t357 = t100 * t220;
t356 = t100 * t224;
t128 = t177 + t281;
t284 = t128 * t223;
t308 = t178 ^ 2;
t315 = -t307 - t308;
t326 = t219 * t315 - t284;
t164 = t180 * t201;
t261 = t224 * t268;
t267 = t220 * qJDD(2);
t183 = t261 + t267;
t256 = -t223 * qJDD(3) + t219 * t183;
t240 = qJD(4) * t180 + t256;
t319 = -t164 + t240;
t285 = t128 * t219;
t325 = t223 * t315 + t285;
t341 = t220 * t319 + t224 * t325;
t355 = -pkin(2) * t326 + pkin(8) * t341;
t110 = t164 + t240;
t158 = t308 - t307;
t354 = t224 * t110 + t220 * (t158 * t223 + t287);
t242 = -t219 * qJDD(3) - t223 * t183;
t139 = -qJD(4) * t178 - t242;
t282 = t178 * t201;
t320 = t139 + t282;
t290 = t320 * t219;
t317 = t176 - t308;
t353 = t220 * (t319 * t223 + t290) + t224 * t317;
t214 = sin(pkin(6));
t215 = cos(pkin(6));
t221 = sin(qJ(2));
t352 = t215 * (t220 * t325 - t224 * t319) + (t221 * t341 - t225 * t326) * t214;
t350 = pkin(3) * t326;
t349 = pkin(9) * t325;
t348 = pkin(9) * t326;
t347 = qJ(5) * t320;
t159 = -t176 + t307;
t345 = t223 * t159 - t285;
t343 = t158 * t219 - t286;
t321 = t139 - t282;
t339 = t220 * (-t159 * t219 - t284) - t224 * t321;
t316 = t176 + t308;
t338 = pkin(3) * t316;
t294 = sin(pkin(11));
t295 = cos(pkin(11));
t236 = t294 * g(1) - t295 * g(2);
t273 = -g(3) + qJDD(1);
t230 = -t214 * t236 + t215 * t273;
t156 = t224 * t230;
t188 = -t295 * g(1) - t294 * g(2);
t233 = t215 * t236;
t327 = t214 * t273 + t233;
t133 = t225 * t188 + t327 * t221;
t226 = qJD(2) ^ 2;
t119 = -t226 * pkin(2) + qJDD(2) * pkin(8) + t133;
t253 = -pkin(3) * t224 - pkin(9) * t220;
t254 = t226 * t253 + t119;
t306 = qJD(3) ^ 2;
t82 = -qJDD(3) * pkin(3) - t306 * pkin(9) + t220 * t254 - t156;
t232 = t240 * pkin(4) - t347 + t82;
t218 = sin(qJ(6));
t222 = cos(qJ(6));
t144 = -t222 * t178 + t180 * t218;
t146 = t178 * t218 + t180 * t222;
t105 = t146 * t144;
t171 = qJDD(6) + t177;
t322 = -t105 + t171;
t337 = t218 * t322;
t335 = t220 * t316;
t333 = t222 * t322;
t331 = t224 * t316;
t324 = -t319 * t219 + t223 * t320;
t195 = qJD(6) + t201;
t123 = t195 * t144;
t79 = -t144 * qJD(6) + t222 * t139 + t218 * t240;
t323 = -t123 + t79;
t147 = pkin(4) * t178 - qJ(5) * t180;
t228 = t220 * t230;
t83 = -t306 * pkin(3) + qJDD(3) * pkin(9) + t224 * t254 + t228;
t248 = t221 * t188 - t327 * t225;
t118 = -qJDD(2) * pkin(2) - t226 * pkin(8) + t248;
t244 = -t184 + t204;
t245 = t183 + t261;
t88 = pkin(3) * t244 - pkin(9) * t245 + t118;
t48 = t219 * t88 + t223 * t83;
t255 = -t177 * qJ(5) - t178 * t147 + t48;
t313 = -(t318 + t307) * pkin(4) - qJ(5) * t234 + t255;
t305 = pkin(4) + pkin(5);
t47 = t219 * t83 - t223 * t88;
t44 = t177 * pkin(4) - qJ(5) * t307 + t180 * t147 + qJDD(5) + t47;
t28 = t128 * pkin(5) - pkin(10) * t321 + t44;
t251 = pkin(5) * t201 - pkin(10) * t180;
t270 = qJD(5) * t201;
t191 = -0.2e1 * t270;
t247 = t191 + t255;
t43 = -pkin(4) * t307 + t247;
t29 = -pkin(5) * t308 + pkin(10) * t240 - t201 * t251 + t43;
t14 = t218 * t29 - t222 * t28;
t15 = t218 * t28 + t222 * t29;
t8 = -t14 * t222 + t15 * t218;
t9 = t218 * t14 + t222 * t15;
t312 = qJ(5) * t9 - t305 * t8;
t143 = t146 ^ 2;
t193 = t195 ^ 2;
t117 = -t143 - t193;
t90 = t105 + t171;
t300 = t218 * t90;
t66 = t117 * t222 - t300;
t297 = t222 * t90;
t67 = -t117 * t218 - t297;
t311 = qJ(5) * t67 - t305 * t66 + t15;
t142 = t144 ^ 2;
t97 = -t193 - t142;
t51 = t218 * t97 + t333;
t52 = t222 * t97 - t337;
t310 = qJ(5) * t52 - t305 * t51 + t14;
t257 = t218 * t139 - t222 * t240;
t58 = (qJD(6) - t195) * t146 + t257;
t61 = t123 + t79;
t32 = -t218 * t58 - t222 * t61;
t34 = t218 * t61 - t222 * t58;
t309 = qJ(5) * t34 - t305 * t32;
t304 = pkin(4) * t223;
t259 = -pkin(4) * t201 - 0.2e1 * qJD(5);
t36 = t256 * pkin(5) + t308 * pkin(10) + t232 + (pkin(5) * qJD(4) - t251 + t259) * t180;
t301 = t218 * t36;
t299 = t219 * t82;
t298 = t222 * t36;
t296 = t223 * t82;
t293 = qJ(5) * t223;
t289 = t321 * t219;
t288 = t321 * t223;
t280 = t195 * t218;
t279 = t195 * t222;
t278 = t201 * t219;
t277 = t201 * t223;
t199 = t220 * t226 * t224;
t190 = qJDD(3) + t199;
t276 = t220 * t190;
t189 = -t199 + qJDD(3);
t274 = t224 * t189;
t265 = t178 * t277;
t264 = t224 * t105;
t263 = t224 * t281;
t260 = -qJ(5) * t219 - pkin(3);
t26 = t219 * t47 + t223 * t48;
t102 = t220 * t119 - t156;
t103 = t224 * t119 + t228;
t63 = t102 * t220 + t224 * t103;
t252 = -pkin(4) * t44 + qJ(5) * t43;
t157 = t180 * t278;
t250 = t220 * (t139 * t223 + t157) - t263;
t249 = -t178 * t278 - t223 * t240;
t246 = -pkin(4) * t321 - qJ(5) * t110;
t25 = t219 * t48 - t223 * t47;
t243 = -pkin(2) + t253;
t238 = (t178 * t219 + t180 * t223) * t201;
t235 = t220 * (-t157 + t265) + t224 * t177;
t231 = 0.2e1 * qJD(5) * t180 - t232;
t229 = -pkin(4) * t128 + qJ(5) * t315 - t44;
t227 = t220 * (t219 * t240 - t265) + t263;
t211 = t224 ^ 2;
t210 = t220 ^ 2;
t209 = t211 * t226;
t207 = t210 * t226;
t197 = -t209 - t306;
t196 = -t207 - t306;
t187 = t207 + t209;
t186 = (t210 + t211) * qJDD(2);
t185 = -0.2e1 * t204 + t266;
t182 = 0.2e1 * t261 + t267;
t155 = -t196 * t220 - t274;
t154 = t197 * t224 - t276;
t121 = -t143 + t193;
t120 = t142 - t193;
t116 = (qJD(4) - t201) * t178 + t242;
t111 = (-qJD(4) - t201) * t180 - t256;
t107 = t139 * t219 - t180 * t277;
t104 = t143 - t142;
t85 = (-t144 * t222 + t146 * t218) * t195;
t84 = (t144 * t218 + t146 * t222) * t195;
t80 = -t142 - t143;
t78 = -qJD(6) * t146 - t257;
t77 = t111 * t223 + t289;
t76 = -t110 * t223 + t289;
t75 = t111 * t219 - t288;
t74 = -t110 * t219 - t288;
t73 = t120 * t222 - t300;
t72 = -t121 * t218 + t333;
t71 = -t120 * t218 - t297;
t70 = -t121 * t222 - t337;
t68 = -t116 * t220 + t356;
t64 = -t220 * t320 - t356;
t57 = (qJD(6) + t195) * t146 + t257;
t56 = -t146 * t280 + t222 * t79;
t55 = -t146 * t279 - t218 * t79;
t54 = t144 * t279 - t218 * t78;
t53 = -t144 * t280 - t222 * t78;
t50 = t224 * t77 - t335;
t49 = t224 * t76 - t335;
t45 = t180 * t259 + t232;
t42 = (-t319 + t164) * pkin(4) + t231;
t41 = pkin(4) * t164 + t231 + t347;
t40 = t219 * t66 + t223 * t67;
t39 = t219 * t67 - t223 * t66;
t38 = qJ(5) * t316 + t44;
t37 = (t316 - t307) * pkin(4) + t247;
t35 = -t218 * t323 - t222 * t57;
t33 = t218 * t57 - t222 * t323;
t31 = t219 * t51 + t223 * t52;
t30 = t219 * t52 - t223 * t51;
t24 = -t220 * t323 + t224 * t40;
t23 = t220 * t82 + t224 * t26;
t22 = -t220 * t57 + t224 * t31;
t21 = t219 * t44 + t223 * t43;
t20 = t219 * t43 - t223 * t44;
t19 = -pkin(10) * t66 + qJ(5) * t323 - t298;
t18 = t219 * t32 + t223 * t34;
t17 = t219 * t34 - t223 * t32;
t16 = -pkin(10) * t51 + qJ(5) * t57 - t301;
t13 = t18 * t224 - t220 * t80;
t12 = t21 * t224 + t220 * t45;
t11 = -pkin(10) * t67 + t305 * t323 + t301;
t10 = -pkin(10) * t52 + t305 * t57 - t298;
t7 = -pkin(10) * t8 - qJ(5) * t36;
t6 = -pkin(10) * t32 + qJ(5) * t80 - t8;
t5 = -pkin(10) * t34 + t305 * t80 - t9;
t4 = -pkin(10) * t9 - t305 * t36;
t3 = t219 * t8 + t223 * t9;
t2 = t219 * t9 - t223 * t8;
t1 = t220 * t36 + t224 * t3;
t27 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t273, 0, 0, 0, 0, 0, 0, (qJDD(2) * t225 - t221 * t226) * t214, (-qJDD(2) * t221 - t225 * t226) * t214, 0, t215 ^ 2 * t273 + (t221 * t133 - t225 * t248 - t233) * t214, 0, 0, 0, 0, 0, 0, t215 * (t190 * t224 + t197 * t220) + (t154 * t221 + t185 * t225) * t214, t215 * (-t189 * t220 + t196 * t224) + (t155 * t221 - t182 * t225) * t214, (t186 * t221 + t187 * t225) * t214, t215 * (-t102 * t224 + t103 * t220) + (-t118 * t225 + t221 * t63) * t214, 0, 0, 0, 0, 0, 0, t352, t215 * (t116 * t224 + t357) + (t221 * t68 - t358) * t214, t215 * (t220 * t77 + t331) + (t221 * t50 - t225 * t75) * t214, t215 * (t220 * t26 - t224 * t82) + (t221 * t23 - t225 * t25) * t214, 0, 0, 0, 0, 0, 0, t352, t215 * (t220 * t76 + t331) + (t221 * t49 - t225 * t74) * t214, t215 * (t224 * t320 - t357) + (t221 * t64 + t358) * t214, t215 * (t21 * t220 - t224 * t45) + (t12 * t221 - t20 * t225) * t214, 0, 0, 0, 0, 0, 0, t215 * (t220 * t31 + t224 * t57) + (t22 * t221 - t225 * t30) * t214, t215 * (t220 * t40 + t224 * t323) + (t221 * t24 - t225 * t39) * t214, t215 * (t18 * t220 + t224 * t80) + (t13 * t221 - t17 * t225) * t214, t215 * (t220 * t3 - t224 * t36) + (t1 * t221 - t2 * t225) * t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t248, -t133, 0, 0, t245 * t220, t182 * t224 + t185 * t220, t276 + t224 * (-t207 + t306), -t244 * t224, t220 * (t209 - t306) + t274, 0, pkin(2) * t185 + pkin(8) * t154 - t118 * t224, -pkin(2) * t182 + pkin(8) * t155 + t118 * t220, pkin(2) * t187 + pkin(8) * t186 + t63, -pkin(2) * t118 + pkin(8) * t63, t250, -t353, t339, t227, t354, t235, t220 * (t299 - t348) + t224 * (t47 - t350) + t355, t220 * (t296 - t360) + t224 * (t48 - t361) - t362 + pkin(8) * t68, pkin(8) * t50 - t220 * t25 + t243 * t75, pkin(8) * t23 + t243 * t25, t250, t339, t353, t235, -t354, t227, t220 * (-t219 * t42 - t293 * t319 - t348) + t224 * (-t229 - t350) + t355, t220 * (-pkin(9) * t74 - t219 * t37 + t223 * t38) + t224 * (-pkin(3) * t74 - t246) - pkin(2) * t74 + pkin(8) * t49, t220 * (-pkin(4) * t290 + t223 * t41 + t360) + t224 * (0.2e1 * t270 - t313 + t361) + t362 + pkin(8) * t64, t220 * (-pkin(9) * t20 + (pkin(4) * t219 - t293) * t45) + t224 * (-pkin(3) * t20 - t252) - pkin(2) * t20 + pkin(8) * t12, t220 * (-t219 * t55 + t223 * t56) + t264, t220 * (-t219 * t33 + t223 * t35) + t224 * t104, t220 * (-t219 * t70 + t223 * t72) + t224 * t61, t220 * (-t219 * t53 + t223 * t54) - t264, t220 * (-t219 * t71 + t223 * t73) - t224 * t58, t220 * (-t219 * t84 + t223 * t85) + t224 * t171, t220 * (-pkin(9) * t30 - t10 * t219 + t16 * t223) + t224 * (-pkin(3) * t30 - t310) - pkin(2) * t30 + pkin(8) * t22, t220 * (-pkin(9) * t39 - t11 * t219 + t19 * t223) + t224 * (-pkin(3) * t39 - t311) - pkin(2) * t39 + pkin(8) * t24, t220 * (-pkin(9) * t17 - t219 * t5 + t223 * t6) + t224 * (-pkin(3) * t17 - t309) - pkin(2) * t17 + pkin(8) * t13, t220 * (-pkin(9) * t2 - t219 * t4 + t223 * t7) + t224 * (-pkin(3) * t2 - t312) - pkin(2) * t2 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, t207 - t209, t267, t199, t266, qJDD(3), -t102, -t103, 0, 0, t107, t324, t345, t249, t343, t238, -pkin(3) * t319 - t296 + t349, pkin(3) * t116 + t299 + t359, pkin(9) * t77 + t26 + t338, -pkin(3) * t82 + pkin(9) * t26, t107, t345, -t324, t238, -t343, t249, t223 * t42 + t260 * t319 + t349, pkin(9) * t76 + t219 * t38 + t223 * t37 + t338, -t359 + t219 * t41 + (pkin(3) + t304) * t320, pkin(9) * t21 + (t260 - t304) * t45, t219 * t56 + t223 * t55, t219 * t35 + t223 * t33, t219 * t72 + t223 * t70, t219 * t54 + t223 * t53, t219 * t73 + t223 * t71, t219 * t85 + t223 * t84, pkin(3) * t57 + pkin(9) * t31 + t10 * t223 + t16 * t219, pkin(3) * t323 + pkin(9) * t40 + t11 * t223 + t19 * t219, pkin(3) * t80 + pkin(9) * t18 + t219 * t6 + t223 * t5, -pkin(3) * t36 + pkin(9) * t3 + t219 * t7 + t223 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t281, t317, t321, -t281, -t110, -t177, -t47, -t48, 0, 0, t281, t321, -t317, -t177, t110, -t281, t229, t246, t191 + t313, t252, -t105, -t104, -t61, t105, t58, -t171, t310, t311, t309, t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t321, t318, t44, 0, 0, 0, 0, 0, 0, t51, t66, t32, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t104, t61, -t105, -t58, t171, -t14, -t15, 0, 0;];
tauJ_reg  = t27;
