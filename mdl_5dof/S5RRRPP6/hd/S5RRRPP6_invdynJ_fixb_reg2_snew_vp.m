% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:01
% EndTime: 2019-12-31 21:02:15
% DurationCPUTime: 5.85s
% Computational Cost: add. (13469->367), mult. (28103->455), div. (0->0), fcn. (19349->8), ass. (0->246)
t220 = sin(qJ(2));
t223 = cos(qJ(2));
t219 = sin(qJ(3));
t222 = cos(qJ(3));
t272 = qJD(1) * t220;
t191 = -t222 * qJD(2) + t219 * t272;
t208 = t220 * qJDD(1);
t266 = qJD(1) * qJD(2);
t261 = t223 * t266;
t196 = t208 + t261;
t239 = -t219 * qJDD(2) - t222 * t196;
t164 = -t191 * qJD(3) - t239;
t216 = sin(pkin(8));
t217 = cos(pkin(8));
t193 = t219 * qJD(2) + t222 * t272;
t240 = t222 * qJDD(2) - t219 * t196;
t230 = -t193 * qJD(3) + t240;
t120 = t217 * t164 + t216 * t230;
t169 = t217 * t191 + t216 * t193;
t205 = t223 * qJD(1) - qJD(3);
t290 = t169 * t205;
t325 = t120 + t290;
t207 = t220 * t266;
t265 = t223 * qJDD(1);
t197 = -t207 + t265;
t190 = -qJDD(3) + t197;
t171 = -t216 * t191 + t217 * t193;
t289 = t171 * t169;
t232 = t190 - t289;
t283 = t216 * t232;
t168 = t171 ^ 2;
t310 = t205 ^ 2;
t322 = -t168 - t310;
t63 = -t217 * t322 - t283;
t280 = t217 * t232;
t69 = -t216 * t322 + t280;
t41 = t219 * t69 - t222 * t63;
t43 = t219 * t63 + t222 * t69;
t381 = -pkin(6) * (t220 * t325 + t223 * t43) + pkin(1) * t41;
t379 = pkin(2) * t41;
t378 = pkin(7) * t41;
t377 = pkin(2) * t325 - pkin(7) * t43;
t312 = t169 ^ 2;
t145 = t312 - t310;
t74 = t216 * t145 - t280;
t78 = t217 * t145 + t283;
t153 = t171 * t205;
t257 = t216 * t164 - t217 * t230;
t90 = -t257 - t153;
t375 = t220 * (t219 * t74 - t222 * t78) + t223 * t90;
t321 = t168 - t312;
t324 = -t153 + t257;
t50 = -t216 * t324 + t217 * t325;
t299 = t216 * t325;
t52 = t217 * t324 + t299;
t374 = t220 * (t219 * t50 + t222 * t52) + t223 * t321;
t373 = pkin(3) * t63;
t320 = -t290 + t120;
t340 = t216 * t320 + t217 * t90;
t341 = t216 * t90 - t217 * t320;
t354 = t219 * t340 + t222 * t341;
t372 = pkin(2) * t354;
t371 = pkin(7) * t354;
t370 = qJ(4) * t63;
t369 = qJ(4) * t69;
t102 = -t312 - t168;
t355 = -t219 * t341 + t222 * t340;
t368 = -pkin(2) * t102 + pkin(7) * t355;
t367 = t219 * t52 - t222 * t50;
t365 = t219 * t78 + t222 * t74;
t363 = pkin(6) * (t220 * t102 + t223 * t355) - pkin(1) * t354;
t113 = t190 + t289;
t282 = t216 * t113;
t319 = -t310 - t312;
t328 = t217 * t319 + t282;
t108 = t217 * t113;
t329 = t216 * t319 - t108;
t338 = t219 * t328 + t222 * t329;
t362 = pkin(2) * t338;
t307 = pkin(3) * t341;
t361 = pkin(7) * t338;
t358 = qJ(4) * t341;
t339 = -t219 * t329 + t222 * t328;
t357 = -pkin(2) * t324 + pkin(7) * t339;
t356 = -pkin(3) * t102 + qJ(4) * t340;
t146 = -t168 + t310;
t342 = t217 * t146 - t282;
t343 = -t216 * t146 - t108;
t353 = t219 * t343 + t222 * t342;
t352 = pkin(6) * (t220 * t324 + t223 * t339) - pkin(1) * t338;
t351 = t220 * (-t219 * t342 + t222 * t343) - t223 * t320;
t306 = pkin(3) * t329;
t345 = qJ(5) * t325;
t225 = qJD(1) ^ 2;
t221 = sin(qJ(1));
t224 = cos(qJ(1));
t248 = t224 * g(1) + t221 * g(2);
t291 = qJDD(1) * pkin(6);
t184 = -t225 * pkin(1) - t248 + t291;
t249 = -t223 * pkin(2) - t220 * pkin(7);
t253 = t225 * t249 + t184;
t303 = t223 * g(3);
t309 = qJD(2) ^ 2;
t141 = -qJDD(2) * pkin(2) - t309 * pkin(7) + t253 * t220 + t303;
t175 = -t205 * pkin(3) - t193 * qJ(4);
t311 = t191 ^ 2;
t71 = -t230 * pkin(3) - t311 * qJ(4) + t193 * t175 + qJDD(4) + t141;
t348 = pkin(4) * t257 - t345 + t71;
t347 = qJ(4) * t328;
t346 = qJ(4) * t329;
t288 = t193 * t191;
t231 = -t190 - t288;
t327 = t219 * t231;
t326 = t222 * t231;
t271 = qJD(4) * t169;
t161 = -0.2e1 * t271;
t269 = qJD(5) * t205;
t323 = t161 - 0.2e1 * t269;
t180 = t191 * t205;
t136 = t164 - t180;
t125 = t169 * pkin(4) - t171 * qJ(5);
t259 = t221 * g(1) - t224 * g(2);
t183 = qJDD(1) * pkin(1) + t225 * pkin(6) + t259;
t244 = -t197 + t207;
t245 = t196 + t261;
t131 = pkin(2) * t244 - pkin(7) * t245 - t183;
t304 = t220 * g(3);
t142 = -t309 * pkin(2) + qJDD(2) * pkin(7) + t223 * t253 - t304;
t99 = -t222 * t131 + t219 * t142;
t60 = t231 * pkin(3) - t136 * qJ(4) - t99;
t100 = t219 * t131 + t222 * t142;
t62 = -t311 * pkin(3) + qJ(4) * t230 + t205 * t175 + t100;
t302 = t216 * t60 + t217 * t62;
t254 = -t190 * qJ(5) - t169 * t125 + t302;
t317 = t373 - pkin(4) * (t322 + t310) - qJ(5) * t232 + t254;
t132 = (qJD(3) + t205) * t193 - t240;
t235 = (t169 * t216 + t171 * t217) * t205;
t287 = t205 * t216;
t144 = t171 * t287;
t286 = t205 * t217;
t264 = t169 * t286;
t246 = -t144 + t264;
t316 = t219 * t246 + t222 * t235;
t237 = t216 * t257 - t264;
t247 = -t169 * t287 - t217 * t257;
t315 = t219 * t237 + t222 * t247;
t181 = t223 * t190;
t314 = t220 * (-t219 * t235 + t222 * t246) + t181;
t263 = t223 * t289;
t313 = t220 * (-t219 * t247 + t222 * t237) + t263;
t189 = t193 ^ 2;
t258 = t216 * t62 - t217 * t60;
t270 = qJD(4) * t171;
t31 = t258 + 0.2e1 * t270;
t32 = t161 + t302;
t16 = t216 * t32 - t217 * t31;
t308 = pkin(3) * t16;
t305 = pkin(4) * t217;
t301 = t216 * t71;
t297 = t217 * t71;
t295 = t219 * t16;
t294 = t222 * t16;
t292 = qJ(5) * t217;
t285 = t205 * t219;
t284 = t205 * t222;
t279 = t219 * t141;
t155 = t190 - t288;
t278 = t219 * t155;
t204 = t223 * t225 * t220;
t277 = t220 * (qJDD(2) + t204);
t276 = t222 * t141;
t275 = t222 * t155;
t274 = t223 * (-t204 + qJDD(2));
t268 = qJD(3) - t205;
t262 = t223 * t288;
t260 = -qJ(5) * t216 - pkin(3);
t17 = t216 * t31 + t217 * t32;
t57 = t222 * t100 + t219 * t99;
t176 = t220 * t184 + t303;
t177 = t223 * t184 - t304;
t256 = t220 * t176 + t223 * t177;
t255 = (0.2e1 * qJD(4) + t125) * t171;
t252 = -t302 - t373;
t84 = t216 * t120 - t171 * t286;
t85 = t217 * t120 + t144;
t250 = t220 * (-t219 * t84 + t222 * t85) - t263;
t243 = t219 * t100 - t222 * t99;
t241 = -pkin(1) + t249;
t238 = t254 + t323;
t236 = t190 * pkin(4) - qJ(5) * t310 + qJDD(5) + t258;
t22 = -pkin(4) * t310 + t238;
t23 = t255 + t236;
t10 = t216 * t22 - t217 * t23;
t234 = pkin(3) * t10 - pkin(4) * t23 + qJ(5) * t22;
t233 = -pkin(4) * t320 + qJ(5) * t90 + t307;
t227 = pkin(4) * t113 - qJ(5) * t319 + t236 - t306;
t226 = 0.2e1 * qJD(5) * t171 - t348;
t214 = t223 ^ 2;
t213 = t220 ^ 2;
t211 = t214 * t225;
t209 = t213 * t225;
t198 = -0.2e1 * t207 + t265;
t195 = t208 + 0.2e1 * t261;
t179 = -t189 + t310;
t178 = -t310 + t311;
t173 = t189 - t311;
t172 = -t189 - t310;
t165 = -t310 - t311;
t163 = -0.2e1 * t270;
t162 = 0.2e1 * t271;
t154 = t189 + t311;
t137 = t191 * t268 + t239;
t135 = t164 + t180;
t133 = -t193 * t268 + t240;
t124 = -t219 * t172 + t275;
t123 = t222 * t172 + t278;
t111 = t222 * t165 - t327;
t110 = t219 * t165 + t326;
t97 = -t132 * t222 + t219 * t136;
t47 = t219 * t85 + t222 * t84;
t45 = t297 + t370;
t40 = t301 - t346;
t35 = -pkin(3) * t325 + t301 + t369;
t34 = (-pkin(4) * t205 - 0.2e1 * qJD(5)) * t171 + t348;
t33 = -pkin(3) * t324 - t297 + t347;
t25 = (-t324 + t153) * pkin(4) + t226;
t24 = pkin(4) * t153 + t226 + t345;
t21 = -qJ(5) * t102 + t23;
t20 = (-t102 - t310) * pkin(4) + t238;
t19 = -t216 * t25 - t292 * t324 - t346;
t18 = -pkin(4) * t299 + t217 * t24 - t370;
t15 = t217 * t25 + t260 * t324 + t347;
t14 = -t369 + t216 * t24 + (pkin(3) + t305) * t325;
t13 = -pkin(3) * t71 + qJ(4) * t17;
t12 = -t16 - t358;
t11 = t216 * t23 + t217 * t22;
t9 = t17 + t356;
t8 = -t216 * t20 + t217 * t21 - t358;
t7 = t217 * t20 + t216 * t21 + t356;
t6 = t222 * t17 - t295;
t5 = t219 * t17 + t294;
t4 = -qJ(4) * t10 + (pkin(4) * t216 - t292) * t34;
t3 = -t219 * t10 + t222 * t11;
t2 = t222 * t10 + t219 * t11;
t1 = qJ(4) * t11 + (t260 - t305) * t34;
t26 = [0, 0, 0, 0, 0, qJDD(1), t259, t248, 0, 0, t245 * t220, t223 * t195 + t220 * t198, t277 + t223 * (-t209 + t309), -t244 * t223, t220 * (t211 - t309) + t274, 0, t223 * t183 + pkin(1) * t198 + pkin(6) * (t223 * (-t211 - t309) - t277), -t220 * t183 - pkin(1) * t195 + pkin(6) * (-t274 - t220 * (-t209 - t309)), pkin(1) * (t209 + t211) + (t213 + t214) * t291 + t256, pkin(1) * t183 + pkin(6) * t256, t220 * (t222 * t164 + t193 * t285) - t262, t220 * (t222 * t133 - t219 * t135) - t223 * t173, t220 * (-t219 * t179 + t326) - t223 * t136, t220 * (-t191 * t284 - t219 * t230) + t262, t220 * (t222 * t178 + t278) + t223 * t132, t181 + t220 * (t191 * t222 - t193 * t219) * t205, t220 * (-pkin(7) * t110 + t279) + t223 * (-pkin(2) * t110 + t99) - pkin(1) * t110 + pkin(6) * (t223 * t111 - t220 * t133), t220 * (-pkin(7) * t123 + t276) + t223 * (-pkin(2) * t123 + t100) - pkin(1) * t123 + pkin(6) * (t223 * t124 - t220 * t137), -t220 * t243 + pkin(6) * (-t220 * t154 + t223 * t97) + t241 * (-t132 * t219 - t222 * t136), pkin(6) * (t220 * t141 + t223 * t57) + t241 * t243, t250, -t374, t351, t313, -t375, t314, t220 * (-t219 * t33 + t222 * t40 - t361) + t223 * (-t306 + t31 - t362) + t352, t220 * (-t219 * t35 + t222 * t45 - t378) + t223 * (t161 - t252 - t379) - t381, t220 * (t222 * t12 - t219 * t9 - t371) + t223 * (-t307 - t372) + t363, t220 * (-pkin(7) * t5 - qJ(4) * t294 - t219 * t13) + t223 * (-pkin(2) * t5 - t308) - pkin(1) * t5 + pkin(6) * (t220 * t71 + t223 * t6), t250, t351, t374, t314, t375, t313, t220 * (-t219 * t15 + t222 * t19 - t361) + t223 * (t227 + t255 - t362) + t352, t220 * (-t219 * t7 + t222 * t8 - t371) + t223 * (-t233 - t372) + t363, t220 * (-t219 * t14 + t222 * t18 + t378) + t223 * (t162 + 0.2e1 * t269 - t317 + t379) + t381, t220 * (-pkin(7) * t2 - t219 * t1 + t222 * t4) + t223 * (-pkin(2) * t2 - t234) - pkin(1) * t2 + pkin(6) * (t220 * t34 + t223 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, t209 - t211, t208, t204, t265, qJDD(2), -t176, -t177, 0, 0, t219 * t164 - t193 * t284, t219 * t133 + t222 * t135, t222 * t179 + t327, -t191 * t285 + t222 * t230, t219 * t178 - t275, (t191 * t219 + t193 * t222) * t205, pkin(2) * t133 + pkin(7) * t111 - t276, pkin(2) * t137 + pkin(7) * t124 + t279, pkin(2) * t154 + pkin(7) * t97 + t57, -pkin(2) * t141 + pkin(7) * t57, t47, -t367, t353, t315, t365, t316, t219 * t40 + t222 * t33 + t357, t219 * t45 + t222 * t35 - t377, t219 * t12 + t222 * t9 + t368, -pkin(2) * t71 + pkin(7) * t6 - qJ(4) * t295 + t222 * t13, t47, t353, t367, t316, -t365, t315, t222 * t15 + t219 * t19 + t357, t219 * t8 + t222 * t7 + t368, t222 * t14 + t219 * t18 + t377, -pkin(2) * t34 + pkin(7) * t3 + t222 * t1 + t219 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, t173, t136, -t288, -t132, -t190, -t99, -t100, 0, 0, t289, t321, t320, -t289, t90, -t190, t163 - t258 + t306, t162 + t252, t307, t308, t289, t320, -t321, -t190, -t90, -t289, -t171 * t125 + t163 - t227, t233, t317 + t323, t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t325, t102, t71, 0, 0, 0, 0, 0, 0, t324, t102, -t325, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t320, t322, t23;];
tauJ_reg = t26;
