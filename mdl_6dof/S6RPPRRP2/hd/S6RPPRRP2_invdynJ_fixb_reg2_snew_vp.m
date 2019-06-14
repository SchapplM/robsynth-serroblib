% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:48:08
% EndTime: 2019-05-05 14:48:24
% DurationCPUTime: 6.09s
% Computational Cost: add. (13451->390), mult. (30244->525), div. (0->0), fcn. (21596->10), ass. (0->241)
t210 = sin(pkin(10));
t212 = cos(pkin(10));
t216 = sin(qJ(4));
t218 = cos(qJ(4));
t264 = qJD(1) * t212;
t271 = t210 * t216;
t185 = qJD(1) * t271 - t218 * t264;
t181 = qJD(5) + t185;
t301 = t181 ^ 2;
t238 = t210 * t218 + t212 * t216;
t187 = t238 * qJD(1);
t215 = sin(qJ(5));
t217 = cos(qJ(5));
t170 = -t217 * qJD(4) + t187 * t215;
t302 = t170 ^ 2;
t142 = t302 - t301;
t258 = t212 * qJDD(1);
t259 = t210 * qJDD(1);
t241 = t216 * t259 - t218 * t258;
t261 = t187 * qJD(4);
t162 = -t241 - t261;
t154 = qJDD(5) - t162;
t172 = qJD(4) * t215 + t187 * t217;
t274 = t172 * t170;
t107 = -t274 - t154;
t283 = t107 * t215;
t82 = -t142 * t217 - t283;
t146 = t181 * t172;
t184 = t238 * qJDD(1);
t262 = t185 * qJD(4);
t164 = t184 - t262;
t249 = t217 * qJDD(4) - t215 * t164;
t235 = qJD(5) * t172 - t249;
t93 = -t146 + t235;
t369 = t210 * (t216 * t93 + t218 * t82) + t212 * (t216 * t82 - t218 * t93);
t169 = t172 ^ 2;
t315 = -t169 - t301;
t73 = t217 * t315 + t283;
t368 = pkin(2) * t73;
t367 = pkin(3) * t73;
t366 = pkin(4) * t73;
t365 = pkin(8) * t73;
t282 = t107 * t217;
t75 = -t215 * t315 + t282;
t364 = pkin(8) * t75;
t213 = cos(pkin(9));
t363 = t213 * t73;
t362 = t216 * t75;
t361 = t218 * t75;
t314 = t169 - t302;
t237 = -t215 * qJDD(4) - t217 * t164;
t230 = -qJD(5) * t170 - t237;
t275 = t170 * t181;
t310 = -t275 + t230;
t290 = t215 * t310;
t316 = t146 + t235;
t57 = t217 * t316 + t290;
t358 = t210 * (-t216 * t314 + t218 * t57) + t212 * (t216 * t57 + t218 * t314);
t211 = sin(pkin(9));
t311 = -t274 + t154;
t280 = t311 * t217;
t308 = -t301 - t302;
t322 = t215 * t308 + t280;
t281 = t311 * t215;
t321 = t217 * t308 - t281;
t337 = t216 * t316 + t218 * t321;
t338 = t216 * t321 - t218 * t316;
t351 = -t210 * t338 + t212 * t337;
t355 = pkin(1) * (t211 * t351 - t213 * t322) + qJ(3) * t351 - pkin(2) * t322;
t354 = pkin(7) * t338;
t353 = -t142 * t215 + t282;
t352 = -pkin(3) * t322 + pkin(7) * t337;
t350 = t210 * t337 + t212 * t338;
t309 = t275 + t230;
t143 = -t169 + t301;
t339 = -t143 * t215 + t280;
t349 = t210 * (t216 * t309 + t218 * t339) + t212 * (t216 * t339 - t218 * t309);
t346 = pkin(4) * t322;
t345 = pkin(8) * t321;
t344 = pkin(8) * t322;
t205 = t210 ^ 2;
t206 = t212 ^ 2;
t300 = qJD(1) ^ 2;
t193 = (t205 + t206) * t300;
t340 = t217 * t143 + t281;
t313 = t169 + t302;
t336 = pkin(4) * t313;
t335 = qJ(6) * t310;
t165 = t187 * t185;
t307 = qJDD(4) - t165;
t333 = t216 * t307;
t331 = t216 * t313;
t327 = t218 * t307;
t325 = t218 * t313;
t296 = sin(qJ(1));
t297 = cos(qJ(1));
t231 = t296 * g(1) - t297 * g(2);
t229 = qJDD(1) * pkin(1) + t231;
t232 = t297 * g(1) + t296 * g(2);
t191 = -t300 * pkin(1) - t232;
t270 = t213 * t191;
t223 = qJDD(1) * qJ(3) + t211 * t229 + t270;
t266 = -g(3) + qJDD(2);
t298 = 2 * qJD(3);
t130 = t212 * (-t300 * pkin(2) + t223) + t210 * t266 + t264 * t298;
t260 = t206 * t300;
t122 = -pkin(3) * t260 + pkin(7) * t258 + t130;
t227 = -t211 * t231 - t270;
t250 = t212 * t266;
t252 = pkin(1) * t211 + qJ(3);
t318 = t252 + pkin(7);
t221 = t250 + (-t318 * qJDD(1) + (-(2 * qJD(3)) + (t212 * pkin(3) + pkin(2)) * qJD(1)) * qJD(1) + t227) * t210;
t79 = t218 * t122 + t216 * t221;
t320 = -t215 * t316 + t217 * t310;
t188 = t213 * t229;
t248 = -t211 * t191 + t188;
t150 = -qJDD(1) * pkin(2) - t300 * qJ(3) + qJDD(3) - t248;
t254 = pkin(1) * t213 + pkin(2);
t319 = -t254 * qJDD(1) + t252 * t193 + t150;
t131 = pkin(5) * t170 - qJ(6) * t172;
t155 = pkin(4) * t185 - pkin(8) * t187;
t299 = qJD(4) ^ 2;
t63 = -t299 * pkin(4) + qJDD(4) * pkin(8) - t185 * t155 + t79;
t128 = -pkin(3) * t258 + t150 + (-t205 * t300 - t260) * pkin(7);
t66 = (-t164 + t262) * pkin(8) + (-t162 + t261) * pkin(4) + t128;
t38 = t215 * t66 + t217 * t63;
t247 = t154 * qJ(6) - t170 * t131 + t38;
t306 = -(t315 + t301) * pkin(5) - qJ(6) * t107 + t247;
t272 = t181 * t217;
t256 = t170 * t272;
t236 = t215 * t235 + t256;
t255 = t216 * t274;
t257 = t218 * t274;
t304 = t210 * (t218 * t236 - t255) + t212 * (t216 * t236 + t257);
t273 = t181 * t215;
t140 = t172 * t273;
t243 = t140 - t256;
t303 = t210 * (t216 * t154 + t218 * t243) + t212 * (-t218 * t154 + t216 * t243);
t182 = t185 ^ 2;
t183 = t187 ^ 2;
t295 = pkin(5) * t217;
t78 = t216 * t122 - t218 * t221;
t45 = t216 * t79 - t218 * t78;
t294 = t210 * t45;
t62 = -qJDD(4) * pkin(4) - t299 * pkin(8) + t187 * t155 + t78;
t293 = t215 * t62;
t291 = t215 * t309;
t288 = t217 * t62;
t286 = t217 * t309;
t284 = qJ(6) * t217;
t279 = t128 * t216;
t278 = t128 * t218;
t159 = qJDD(4) + t165;
t276 = t159 * t218;
t268 = t216 * t159;
t263 = qJD(6) * t181;
t253 = -pkin(4) * t218 - pkin(3);
t251 = -qJ(6) * t215 - pkin(4);
t37 = t215 * t63 - t217 * t66;
t18 = t215 * t37 + t217 * t38;
t46 = t216 * t78 + t218 * t79;
t129 = -t250 + ((-pkin(2) * qJD(1) + t298) * qJD(1) + t223) * t210;
t87 = t210 * t129 + t212 * t130;
t173 = 0.2e1 * t263;
t242 = t173 + t247;
t32 = -pkin(5) * t301 + t242;
t33 = -t154 * pkin(5) - qJ(6) * t301 + t131 * t172 + qJDD(6) + t37;
t246 = -pkin(5) * t33 + qJ(6) * t32;
t245 = -pkin(5) * t309 - qJ(6) * t93;
t244 = t170 * t273 - t217 * t235;
t17 = t215 * t38 - t217 * t37;
t233 = (-t170 * t215 - t172 * t217) * t181;
t228 = t235 * pkin(5) - t335 + t62;
t226 = 0.2e1 * qJD(6) * t172 - t228;
t91 = t217 * t230 - t140;
t225 = t210 * (t218 * t91 + t255) + t212 * (t216 * t91 - t257);
t224 = pkin(5) * t311 + qJ(6) * t308 - t33;
t202 = t206 * qJDD(1);
t201 = t205 * qJDD(1);
t192 = t202 + t201;
t176 = -t183 - t299;
t175 = -t183 + t299;
t174 = t182 - t299;
t163 = t184 - 0.2e1 * t262;
t161 = t241 + 0.2e1 * t261;
t156 = -t299 - t182;
t135 = -t182 - t183;
t124 = -t216 * t176 - t276;
t123 = t176 * t218 - t268;
t112 = t216 * t184 - t218 * t241;
t111 = -t184 * t218 - t216 * t241;
t110 = t156 * t218 - t333;
t109 = t216 * t156 + t327;
t100 = (qJD(5) + t181) * t170 + t237;
t95 = (-qJD(5) + t181) * t172 + t249;
t90 = t172 * t272 + t215 * t230;
t84 = -t123 * t210 + t124 * t212;
t72 = -t111 * t210 + t112 * t212;
t67 = -t109 * t210 + t110 * t212;
t60 = t217 * t95 + t291;
t58 = -t217 * t93 + t291;
t56 = t215 * t95 - t286;
t55 = -t215 * t93 - t286;
t53 = -t216 * t100 + t361;
t51 = t100 * t218 + t362;
t49 = -t216 * t310 - t361;
t47 = t218 * t310 - t362;
t44 = t218 * t60 - t331;
t43 = t218 * t58 - t331;
t42 = t216 * t60 + t325;
t41 = t216 * t58 + t325;
t40 = t288 - t365;
t39 = t293 - t344;
t35 = (pkin(5) * t181 - 0.2e1 * qJD(6)) * t172 + t228;
t34 = -pkin(4) * t55 - t245;
t31 = t38 - t366;
t30 = t37 - t346;
t28 = -t210 * t51 + t212 * t53;
t26 = -t210 * t47 + t212 * t49;
t25 = (-t316 - t146) * pkin(5) + t226;
t24 = -pkin(5) * t146 + t226 + t335;
t23 = qJ(6) * t313 + t33;
t22 = t212 * t46 - t294;
t21 = (t313 - t301) * pkin(5) + t242;
t19 = -t210 * t41 + t212 * t43;
t16 = -t224 - t346;
t15 = -t215 * t25 - t284 * t316 - t344;
t14 = -0.2e1 * t263 - t306 + t366;
t13 = -pkin(5) * t290 + t217 * t24 + t365;
t12 = t18 * t218 + t216 * t62;
t11 = t216 * t18 - t218 * t62;
t10 = -pkin(8) * t56 - t17;
t9 = t215 * t33 + t217 * t32;
t8 = t215 * t32 - t217 * t33;
t7 = -pkin(8) * t55 - t21 * t215 + t217 * t23;
t6 = t216 * t35 + t218 * t9;
t5 = t216 * t9 - t218 * t35;
t3 = -pkin(8) * t8 + (pkin(5) * t215 - t284) * t35;
t2 = -pkin(4) * t8 - t246;
t1 = -t210 * t5 + t212 * t6;
t4 = [0, 0, 0, 0, 0, qJDD(1), t231, t232, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t213 - t300 * t211) + t248, (-0.2e1 * qJDD(1) * t211 - t300 * t213) * pkin(1) + t227, 0, pkin(1) * (t211 ^ 2 * t229 + t213 * t188), t201, 0.2e1 * t210 * t258, 0, t202, 0, 0, -t319 * t212, t319 * t210, pkin(2) * t193 + qJ(3) * t192 + pkin(1) * (t192 * t211 + t193 * t213) + t87, -pkin(2) * t150 + qJ(3) * t87 + pkin(1) * (-t150 * t213 + t211 * t87), t210 * (t164 * t218 - t216 * t261) + t212 * (t216 * t164 + t218 * t261), t210 * (-t161 * t218 - t216 * t163) + t212 * (-t216 * t161 + t163 * t218), t210 * (-t216 * t175 + t327) + t212 * (t175 * t218 + t333), t210 * (-t216 * t162 + t218 * t262) + t212 * (t162 * t218 + t216 * t262), t210 * (t174 * t218 - t268) + t212 * (t216 * t174 + t276), (t210 * (-t185 * t218 + t187 * t216) + t212 * (-t185 * t216 - t187 * t218)) * qJD(4), t210 * (-pkin(7) * t109 + t279) + t212 * (-pkin(3) * t161 + pkin(7) * t110 - t278) - pkin(2) * t161 + qJ(3) * t67 + pkin(1) * (-t161 * t213 + t211 * t67), t210 * (-pkin(7) * t123 + t278) + t212 * (-pkin(3) * t163 + pkin(7) * t124 + t279) - pkin(2) * t163 + qJ(3) * t84 + pkin(1) * (-t163 * t213 + t211 * t84), t210 * (-pkin(7) * t111 - t45) + t212 * (-pkin(3) * t135 + pkin(7) * t112 + t46) - pkin(2) * t135 + qJ(3) * t72 + pkin(1) * (-t135 * t213 + t211 * t72), -pkin(7) * t294 + t212 * (-pkin(3) * t128 + pkin(7) * t46) - pkin(2) * t128 + qJ(3) * t22 + pkin(1) * (-t128 * t213 + t211 * t22), t225, -t358, t349, t304, -t369, t303, t210 * (-t216 * t30 + t218 * t39 - t354) + t212 * (t216 * t39 + t218 * t30 + t352) + t355, t210 * (-pkin(7) * t51 - t216 * t31 + t218 * t40) + t212 * (pkin(7) * t53 + t216 * t40 + t218 * t31 - t367) - t368 + qJ(3) * t28 + pkin(1) * (t211 * t28 - t363), t210 * (-pkin(7) * t42 + t10 * t218) + t212 * (pkin(7) * t44 + t216 * t10) + t252 * (-t210 * t42 + t212 * t44) + (pkin(4) * t271 + t212 * t253 - t254) * t56, (t210 * (pkin(4) * t216 - pkin(8) * t218) + t212 * (-pkin(8) * t216 + t253) - t254) * t17 + t318 * (-t11 * t210 + t12 * t212), t225, t349, t358, t303, t369, t304, t210 * (t15 * t218 - t216 * t16 - t354) + t212 * (t216 * t15 + t16 * t218 + t352) + t355, t210 * (-pkin(7) * t41 - t216 * t34 + t218 * t7) + t212 * (-pkin(3) * t55 + pkin(7) * t43 + t216 * t7 + t218 * t34) - pkin(2) * t55 + qJ(3) * t19 + pkin(1) * (t19 * t211 - t213 * t55), t210 * (-pkin(7) * t47 + t13 * t218 - t216 * t14) + t212 * (pkin(7) * t49 + t216 * t13 + t14 * t218 + t367) + t368 + qJ(3) * t26 + pkin(1) * (t211 * t26 + t363), t210 * (-pkin(7) * t5 - t216 * t2 + t218 * t3) + t212 * (-pkin(3) * t8 + pkin(7) * t6 + t2 * t218 + t216 * t3) - pkin(2) * t8 + qJ(3) * t1 + pkin(1) * (t1 * t211 - t213 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129 * t212 + t130 * t210, 0, 0, 0, 0, 0, 0, t109 * t212 + t110 * t210, t123 * t212 + t124 * t210, t111 * t212 + t112 * t210, t210 * t46 + t212 * t45, 0, 0, 0, 0, 0, 0, t350, t210 * t53 + t212 * t51, t210 * t44 + t212 * t42, t11 * t212 + t12 * t210, 0, 0, 0, 0, 0, 0, t350, t210 * t43 + t212 * t41, t210 * t49 + t212 * t47, t210 * t6 + t212 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258, t259, -t193, t150, 0, 0, 0, 0, 0, 0, t161, t163, t135, t128, 0, 0, 0, 0, 0, 0, t322, t73, t56, t17, 0, 0, 0, 0, 0, 0, t322, t55, -t73, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, t183 - t182, t184, -t165, -t241, qJDD(4), -t78, -t79, 0, 0, t90, t320, t340, t244, -t353, t233, -pkin(4) * t316 - t288 + t345, pkin(4) * t100 + t293 + t364, pkin(8) * t60 + t18 + t336, -pkin(4) * t62 + pkin(8) * t18, t90, t340, -t320, t233, t353, t244, t217 * t25 + t251 * t316 + t345, pkin(8) * t58 + t21 * t217 + t215 * t23 + t336, -t364 + t215 * t24 + (pkin(4) + t295) * t310, pkin(8) * t9 + (t251 - t295) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t274, t314, t309, -t274, -t93, t154, -t37, -t38, 0, 0, t274, t309, -t314, t154, t93, -t274, t224, t245, t173 + t306, t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t311, t309, t315, t33;];
tauJ_reg  = t4;
