% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:20:08
% EndTime: 2018-11-23 15:20:19
% DurationCPUTime: 11.31s
% Computational Cost: add. (5890->611), mult. (15370->817), div. (0->0), fcn. (10733->10), ass. (0->276)
t223 = cos(qJ(4));
t271 = qJD(3) * qJD(4);
t220 = sin(qJ(4));
t221 = sin(qJ(3));
t275 = qJD(4) * t221;
t224 = cos(qJ(3));
t277 = qJD(3) * t224;
t354 = t220 * t275 - t223 * t277;
t142 = -qJD(2) * t354 + t223 * t271;
t274 = qJD(4) * t223;
t362 = t220 * t277 + t221 * t274;
t143 = -qJD(2) * t362 - t220 * t271;
t217 = sin(pkin(11));
t300 = cos(pkin(11));
t81 = t142 * t217 - t143 * t300;
t345 = -t81 / 0.2e1;
t82 = t142 * t300 + t217 * t143;
t343 = t82 / 0.2e1;
t369 = Ifges(6,1) + Ifges(7,1);
t368 = Ifges(6,4) - Ifges(7,5);
t357 = Ifges(7,4) + Ifges(6,5);
t250 = t300 * t220;
t179 = t217 * t223 + t250;
t280 = qJD(2) * t224;
t153 = t179 * t280;
t163 = t179 * qJD(4);
t289 = t153 - t163;
t249 = t300 * t223;
t296 = t217 * t220;
t235 = t249 - t296;
t231 = t224 * t235;
t154 = qJD(2) * t231;
t164 = t235 * qJD(4);
t288 = t154 - t164;
t372 = -Ifges(4,1) / 0.2e1;
t214 = Ifges(4,4) * t280;
t371 = -t214 / 0.2e1;
t272 = qJD(2) * qJD(3);
t252 = t221 * t272;
t370 = t369 * t343 + t368 * t345 + t357 * t252 / 0.2e1;
t278 = qJD(3) * t223;
t282 = qJD(2) * t221;
t183 = -t220 * t282 + t278;
t184 = qJD(3) * t220 + t223 * t282;
t116 = -t300 * t183 + t184 * t217;
t340 = t116 / 0.2e1;
t257 = Ifges(4,5) * qJD(3) / 0.2e1;
t367 = Ifges(6,6) - Ifges(7,6);
t206 = qJD(4) - t280;
t236 = t217 * t183 + t184 * t300;
t365 = -t368 * t116 + t357 * t206 + t369 * t236;
t222 = sin(qJ(2));
t218 = sin(pkin(6));
t284 = qJD(1) * t218;
t266 = t222 * t284;
t189 = qJD(2) * pkin(8) + t266;
t219 = cos(pkin(6));
t283 = qJD(1) * t219;
t150 = t224 * t189 + t221 * t283;
t262 = t220 * t280;
t122 = pkin(4) * t262 + t150;
t276 = qJD(4) * t220;
t364 = pkin(4) * t276 - t289 * pkin(5) + t288 * qJ(6) - qJD(6) * t179 - t122;
t149 = -t221 * t189 + t224 * t283;
t225 = cos(qJ(2));
t265 = t225 * t284;
t190 = -qJD(2) * pkin(2) - t265;
t306 = t184 * Ifges(5,4);
t107 = t183 * Ifges(5,2) + t206 * Ifges(5,6) + t306;
t177 = Ifges(5,4) * t183;
t108 = t184 * Ifges(5,1) + t206 * Ifges(5,5) + t177;
t137 = -qJD(3) * pkin(3) - t149;
t138 = qJD(3) * pkin(9) + t150;
t194 = -pkin(3) * t224 - pkin(9) * t221 - pkin(2);
t152 = qJD(2) * t194 - t265;
t85 = -t138 * t220 + t223 * t152;
t87 = t138 * t223 + t152 * t220;
t239 = t220 * t87 + t223 * t85;
t311 = Ifges(5,4) * t223;
t242 = -Ifges(5,2) * t220 + t311;
t312 = Ifges(5,4) * t220;
t244 = Ifges(5,1) * t223 - t312;
t245 = mrSges(5,1) * t220 + mrSges(5,2) * t223;
t309 = Ifges(5,6) * t220;
t310 = Ifges(5,5) * t223;
t326 = t223 / 0.2e1;
t327 = -t220 / 0.2e1;
t328 = t206 / 0.2e1;
t330 = t184 / 0.2e1;
t227 = -t239 * mrSges(5,3) + t244 * t330 + t137 * t245 + t108 * t326 + t107 * t327 + (-t309 + t310) * t328 + t183 * t242 / 0.2e1;
t363 = -t190 * mrSges(4,2) + t149 * mrSges(4,3) + t282 * t372 - t227 - t257 + t371;
t60 = qJ(5) * t183 + t87;
t303 = t217 * t60;
t59 = -qJ(5) * t184 + t85;
t45 = pkin(4) * t206 + t59;
t14 = t300 * t45 - t303;
t10 = -t206 * pkin(5) + qJD(6) - t14;
t104 = -pkin(4) * t183 + qJD(5) + t137;
t35 = pkin(5) * t116 - qJ(6) * t236 + t104;
t361 = -t104 * mrSges(6,2) - mrSges(7,2) * t10 + mrSges(6,3) * t14 + t35 * mrSges(7,3);
t52 = t300 * t60;
t15 = t217 * t45 + t52;
t12 = qJ(6) * t206 + t15;
t39 = Ifges(7,5) * t236 + t206 * Ifges(7,6) + t116 * Ifges(7,3);
t42 = Ifges(6,4) * t236 - t116 * Ifges(6,2) + t206 * Ifges(6,6);
t360 = -t104 * mrSges(6,1) - t35 * mrSges(7,1) + mrSges(7,2) * t12 + mrSges(6,3) * t15 + t42 / 0.2e1 - t39 / 0.2e1;
t359 = -Ifges(6,6) / 0.2e1;
t358 = Ifges(7,6) / 0.2e1;
t344 = t81 / 0.2e1;
t341 = -t116 / 0.2e1;
t337 = t236 / 0.2e1;
t256 = -Ifges(4,6) * qJD(3) / 0.2e1;
t31 = t81 * mrSges(7,1) - t82 * mrSges(7,3);
t32 = t81 * mrSges(6,1) + t82 * mrSges(6,2);
t355 = t31 + t32;
t94 = -mrSges(6,2) * t206 - mrSges(6,3) * t116;
t97 = -mrSges(7,2) * t116 + mrSges(7,3) * t206;
t317 = t94 + t97;
t95 = mrSges(6,1) * t206 - mrSges(6,3) * t236;
t96 = -mrSges(7,1) * t206 + mrSges(7,2) * t236;
t316 = t95 - t96;
t291 = t223 * t224;
t208 = pkin(8) * t291;
t156 = t220 * t194 + t208;
t294 = t218 * t225;
t263 = qJD(2) * t294;
t232 = qJD(1) * (qJD(3) * t219 + t263);
t279 = qJD(3) * t221;
t109 = -t189 * t279 + t224 * t232;
t247 = pkin(3) * t221 - pkin(9) * t224;
t187 = t247 * qJD(3);
t146 = (t187 + t266) * qJD(2);
t23 = t223 * t109 - t138 * t276 + t220 * t146 + t152 * t274;
t24 = -qJD(4) * t87 - t109 * t220 + t223 * t146;
t353 = -t220 * t24 + t223 * t23;
t7 = pkin(4) * t252 - qJ(5) * t142 - qJD(5) * t184 + t24;
t9 = qJ(5) * t143 + qJD(5) * t183 + t23;
t4 = t217 * t7 + t300 * t9;
t1 = qJ(6) * t252 + qJD(6) * t206 + t4;
t3 = -t217 * t9 + t300 * t7;
t2 = -pkin(5) * t252 - t3;
t352 = -t24 * mrSges(5,1) - t3 * mrSges(6,1) + t2 * mrSges(7,1) + t23 * mrSges(5,2) + t4 * mrSges(6,2) - t1 * mrSges(7,3) - Ifges(5,5) * t142 - Ifges(5,6) * t143;
t248 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t269 = t358 + t359;
t270 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t313 = Ifges(4,4) * t221;
t350 = -t269 * t116 - t248 * t206 - t270 * t236 - t12 * mrSges(7,3) - t14 * mrSges(6,1) - t190 * mrSges(4,1) - t85 * mrSges(5,1) - t184 * Ifges(5,5) - t183 * Ifges(5,6) - t256 + (t224 * Ifges(4,2) + t313) * qJD(2) / 0.2e1 - Ifges(6,6) * t341 - Ifges(7,6) * t340 + t10 * mrSges(7,1) + t15 * mrSges(6,2) + t150 * mrSges(4,3) + t87 * mrSges(5,2) - t357 * t337 - (Ifges(5,3) + Ifges(6,3) + Ifges(7,2)) * t328;
t349 = Ifges(7,5) * t343 + Ifges(7,3) * t344 + t252 * t358;
t348 = -t82 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t344 + t252 * t359;
t338 = -t236 / 0.2e1;
t336 = t142 / 0.2e1;
t335 = t143 / 0.2e1;
t332 = -t183 / 0.2e1;
t331 = -t184 / 0.2e1;
t329 = -t206 / 0.2e1;
t325 = pkin(4) * t184;
t324 = pkin(4) * t217;
t323 = pkin(8) * t220;
t322 = pkin(8) * t224;
t320 = -qJ(5) - pkin(9);
t238 = pkin(4) * t221 - qJ(5) * t291;
t273 = qJD(5) * t223;
t286 = t223 * t187 + t279 * t323;
t48 = -t221 * t273 + t238 * qJD(3) + (-t208 + (qJ(5) * t221 - t194) * t220) * qJD(4) + t286;
t287 = t220 * t187 + t194 * t274;
t292 = t221 * t223;
t56 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t292 + (-qJD(5) * t221 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t224) * t220 + t287;
t18 = t217 * t48 + t300 * t56;
t63 = -mrSges(7,2) * t81 + mrSges(7,3) * t252;
t64 = -mrSges(6,2) * t252 - mrSges(6,3) * t81;
t319 = t63 + t64;
t65 = mrSges(6,1) * t252 - mrSges(6,3) * t82;
t66 = -mrSges(7,1) * t252 + t82 * mrSges(7,2);
t318 = t66 - t65;
t186 = t247 * qJD(2);
t98 = -t149 * t220 + t223 * t186;
t83 = qJD(2) * t238 + t98;
t99 = t223 * t149 + t220 * t186;
t91 = -qJ(5) * t262 + t99;
t34 = t217 * t83 + t300 * t91;
t315 = mrSges(5,3) * t183;
t314 = mrSges(5,3) * t184;
t110 = t189 * t277 + t221 * t232;
t295 = t218 * t222;
t165 = -t219 * t224 + t221 * t295;
t297 = t110 * t165;
t293 = t220 * t221;
t290 = t224 * t225;
t181 = t223 * t194;
t120 = -qJ(5) * t292 + t181 + (-pkin(4) - t323) * t224;
t128 = -qJ(5) * t293 + t156;
t62 = t217 * t120 + t300 * t128;
t285 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t183 - mrSges(5,2) * t184 - mrSges(4,3) * t282;
t188 = pkin(4) * t293 + t221 * pkin(8);
t281 = qJD(2) * t222;
t268 = mrSges(4,3) * t280;
t151 = pkin(4) * t362 + pkin(8) * t277;
t213 = -pkin(4) * t223 - pkin(3);
t267 = t300 * pkin(4);
t264 = t218 * t281;
t58 = mrSges(6,1) * t116 + mrSges(6,2) * t236;
t253 = m(6) * t104 + t58;
t251 = qJD(4) * t320;
t246 = mrSges(5,1) * t223 - mrSges(5,2) * t220;
t243 = Ifges(5,1) * t220 + t311;
t241 = Ifges(5,2) * t223 + t312;
t240 = Ifges(5,5) * t220 + Ifges(5,6) * t223;
t166 = t219 * t221 + t224 * t295;
t126 = -t166 * t220 - t223 * t294;
t237 = -t166 * t223 + t220 * t294;
t16 = -t217 * t56 + t300 * t48;
t33 = -t217 * t91 + t300 * t83;
t61 = t120 * t300 - t217 * t128;
t234 = -m(4) * t149 + m(5) * t137 - t285;
t233 = -qJD(5) * t220 + t223 * t251;
t69 = -pkin(4) * t143 + t110;
t226 = qJD(2) ^ 2;
t212 = -t267 - pkin(5);
t209 = qJ(6) + t324;
t205 = Ifges(7,2) * t252;
t204 = Ifges(5,3) * t252;
t203 = Ifges(6,3) * t252;
t197 = t320 * t223;
t196 = -qJD(3) * mrSges(4,2) + t268;
t185 = (-mrSges(4,1) * t224 + mrSges(4,2) * t221) * qJD(2);
t173 = (mrSges(4,1) * t221 + mrSges(4,2) * t224) * t272;
t160 = t220 * t251 + t273;
t159 = t235 * t221;
t158 = t179 * t221;
t155 = -t220 * t322 + t181;
t148 = mrSges(5,1) * t206 - t314;
t147 = -mrSges(5,2) * t206 + t315;
t145 = (t220 * t222 + t223 * t290) * t284;
t144 = (-t220 * t290 + t222 * t223) * t284;
t135 = -t197 * t300 + t296 * t320;
t134 = -t197 * t217 - t250 * t320;
t125 = qJD(3) * t166 + t221 * t263;
t124 = -qJD(3) * t165 + t224 * t263;
t114 = -mrSges(5,2) * t252 + mrSges(5,3) * t143;
t113 = mrSges(5,1) * t252 - mrSges(5,3) * t142;
t111 = -pkin(5) * t235 - qJ(6) * t179 + t213;
t103 = t160 * t300 + t217 * t233;
t102 = t160 * t217 - t233 * t300;
t101 = qJD(3) * t231 - t179 * t275;
t100 = t217 * t354 - t249 * t275 - t250 * t277;
t93 = -qJD(4) * t156 + t286;
t92 = (-t221 * t278 - t224 * t276) * pkin(8) + t287;
t90 = pkin(5) * t158 - qJ(6) * t159 + t188;
t89 = -mrSges(5,1) * t143 + mrSges(5,2) * t142;
t86 = t217 * t144 + t145 * t300;
t84 = -t144 * t300 + t145 * t217;
t79 = Ifges(7,4) * t82;
t78 = Ifges(6,5) * t82;
t77 = Ifges(6,6) * t81;
t76 = Ifges(7,6) * t81;
t72 = t142 * Ifges(5,1) + t143 * Ifges(5,4) + Ifges(5,5) * t252;
t71 = t142 * Ifges(5,4) + t143 * Ifges(5,2) + Ifges(5,6) * t252;
t68 = t217 * t126 - t237 * t300;
t67 = -t126 * t300 - t217 * t237;
t57 = mrSges(7,1) * t116 - mrSges(7,3) * t236;
t55 = t224 * pkin(5) - t61;
t54 = -qJ(6) * t224 + t62;
t51 = qJD(4) * t237 - t124 * t220 + t223 * t264;
t50 = qJD(4) * t126 + t124 * t223 + t220 * t264;
t36 = pkin(5) * t236 + qJ(6) * t116 + t325;
t30 = -pkin(5) * t282 - t33;
t29 = qJ(6) * t282 + t34;
t22 = -pkin(5) * t100 - qJ(6) * t101 - qJD(6) * t159 + t151;
t21 = t300 * t59 - t303;
t20 = t217 * t59 + t52;
t19 = t217 * t51 + t300 * t50;
t17 = t217 * t50 - t300 * t51;
t13 = -pkin(5) * t279 - t16;
t11 = qJ(6) * t279 - qJD(6) * t224 + t18;
t5 = pkin(5) * t81 - qJ(6) * t82 - qJD(6) * t236 + t69;
t6 = [-t166 * mrSges(4,3) * t252 + t126 * t113 - t237 * t114 + t124 * t196 + t50 * t147 + t51 * t148 + t319 * t68 + t318 * t67 + t317 * t19 - t316 * t17 + ((-mrSges(3,2) * t226 - t173) * t225 + (-mrSges(3,1) * t226 + qJD(2) * t185) * t222) * t218 + (qJD(3) * t268 + t355 + t89) * t165 + (t57 + t58 - t285) * t125 + m(4) * (t109 * t166 + t297 + t124 * t150 - t125 * t149 + (t190 - t265) * t264) + m(7) * (t1 * t68 + t10 * t17 + t12 * t19 + t125 * t35 + t165 * t5 + t2 * t67) + m(6) * (t104 * t125 - t14 * t17 + t15 * t19 + t165 * t69 - t3 * t67 + t4 * t68) + m(5) * (t125 * t137 + t126 * t24 - t23 * t237 + t50 * t87 + t51 * t85 + t297); (-t204 / 0.2e1 - t205 / 0.2e1 - t203 / 0.2e1 - t79 / 0.2e1 - t78 / 0.2e1 + t77 / 0.2e1 - t76 / 0.2e1 + t109 * mrSges(4,3) - t270 * t82 - t269 * t81 + t352) * t224 + (t244 * t336 + t242 * t335 + t72 * t326 + t71 * t327 + (-t220 * t23 - t223 * t24) * mrSges(5,3) + (mrSges(4,3) + t245) * t110 + (mrSges(4,2) * t281 + (-m(7) * t35 - t234 - t253 - t57) * t225) * t284 + (t108 * t327 - t223 * t107 / 0.2e1 + t240 * t329 + t241 * t332 + t243 * t331 + t137 * t246 + (t220 * t85 - t223 * t87) * mrSges(5,3)) * qJD(4) + (((-0.3e1 / 0.2e1 * Ifges(4,4) + t310 / 0.2e1 - t309 / 0.2e1) * t221 + t270 * t159 + t269 * t158 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - t248) * t224) * qJD(2) + t256 - t350) * qJD(3) + (t89 + (m(4) + m(5)) * t110 + (-m(4) * t150 - t196) * qJD(3)) * pkin(8)) * t221 + (-t158 * t4 - t159 * t3) * mrSges(6,3) + (Ifges(6,4) * t341 + Ifges(7,5) * t340 + t357 * t328 + t369 * t337 - t361 + t365 / 0.2e1) * t101 + (-t1 * t158 + t159 * t2) * mrSges(7,2) + (-t145 + t92) * t147 + (t257 + 0.3e1 / 0.2e1 * t214 + t234 * pkin(8) - t363) * t277 + (-t158 * t368 + t159 * t369) * t343 + (t104 * t151 + t188 * t69 + t3 * t61 + t4 * t62 + (t18 - t86) * t15 + (t16 + t84) * t14) * m(6) + (Ifges(7,5) * t159 + Ifges(7,3) * t158) * t344 + (Ifges(6,4) * t159 - Ifges(6,2) * t158) * t345 + t158 * t348 + t158 * t349 + (Ifges(6,2) * t341 - Ifges(7,3) * t340 + t328 * t367 + t337 * t368 + t360) * t100 + m(5) * (t155 * t24 + t156 * t23 + t85 * t93 + t87 * t92) + m(7) * (t1 * t54 + t10 * t13 + t11 * t12 + t2 * t55 + t22 * t35 + t5 * t90) - m(7) * (t10 * t84 + t12 * t86) - m(5) * (t144 * t85 + t145 * t87) + t188 * t32 - pkin(2) * t173 + t69 * (mrSges(6,1) * t158 + mrSges(6,2) * t159) + t5 * (mrSges(7,1) * t158 - mrSges(7,3) * t159) + t155 * t113 + t156 * t114 + t151 * t58 + t18 * t94 + t16 * t95 + t13 * t96 + t11 * t97 + t90 * t31 + t54 * t63 + t62 * t64 + t61 * t65 + t55 * t66 + t22 * t57 + m(4) * (-pkin(2) * qJD(1) * t264 + t109 * t322) + t316 * t84 - t317 * t86 + t159 * t370 + (-t144 + t93) * t148 - m(4) * (t150 * t224 * t265 + t190 * t266) + (-t196 * t290 + (-mrSges(4,1) * t280 - t185) * t222) * t284; m(5) * (-pkin(3) * t110 + pkin(9) * t353) + t353 * mrSges(5,3) + (Ifges(6,4) * t164 + Ifges(7,5) * t154 - Ifges(6,2) * t163 + Ifges(7,3) * t153) * t341 + (Ifges(6,4) * t154 + Ifges(7,5) * t164 - Ifges(6,2) * t153 + Ifges(7,3) * t163) * t340 + (-t104 * t122 - t134 * t3 + t135 * t4 + t213 * t69 + (t103 - t34) * t15 + (-t102 - t33) * t14) * m(6) + (t39 - t42) * (t163 / 0.2e1 - t153 / 0.2e1) + (Ifges(7,5) * t179 - Ifges(7,3) * t235) * t344 + (Ifges(6,4) * t179 + Ifges(6,2) * t235) * t345 + t69 * (-mrSges(6,1) * t235 + mrSges(6,2) * t179) + t5 * (-mrSges(7,1) * t235 - mrSges(7,3) * t179) + (t14 * t288 + t15 * t289 - t179 * t3 + t235 * t4) * mrSges(6,3) + (t1 * t235 - t10 * t288 + t12 * t289 + t179 * t2) * mrSges(7,2) - t235 * t348 - t235 * t349 + (t1 * t135 + t111 * t5 + t134 * t2 + t364 * t35 + (t103 - t29) * t12 + (t102 - t30) * t10) * m(7) + t364 * t57 + ((t257 + t371 + t363) * t224 + (t256 + (t313 / 0.2e1 + (t372 + Ifges(4,2) / 0.2e1) * t224) * qJD(2) + t350) * t221 + (t179 * t357 + t235 * t367 + t240) * t279 / 0.2e1) * qJD(2) + (t179 * t369 + t235 * t368) * t343 + (-t153 * t368 + t154 * t369) * t338 + (-t163 * t368 + t164 * t369) * t337 + t365 * (t164 / 0.2e1 - t154 / 0.2e1) + (-t163 * t367 + t164 * t357) * t328 + (-t153 * t367 + t154 * t357) * t329 + (-mrSges(4,1) - t246) * t110 - m(5) * (t137 * t150 + t85 * t98 + t87 * t99) + t220 * t72 / 0.2e1 + t213 * t32 - t149 * t196 - t99 * t147 - t98 * t148 - t122 * t58 - t109 * mrSges(4,2) + t111 * t31 - t34 * t94 - t33 * t95 - t30 * t96 - t29 * t97 - pkin(3) * t89 - t316 * t102 + t317 * t103 + t318 * t134 + t319 * t135 + t241 * t335 + t243 * t336 + t71 * t326 + t179 * t370 + (-t113 * t220 + t114 * t223) * pkin(9) + (t253 * t220 * pkin(4) + (-m(5) * t239 - t220 * t147 - t223 * t148) * pkin(9) + t227) * qJD(4) + t285 * t150 + (-mrSges(7,1) * t289 + mrSges(7,3) * t288) * t35 + (-mrSges(6,1) * t289 - mrSges(6,2) * t288) * t104; -t352 + (t1 * t209 - t10 * t20 + t2 * t212 - t35 * t36 + (qJD(6) - t21) * t12) * m(7) + (t148 + t314) * t87 + (-Ifges(6,2) * t340 + Ifges(7,3) * t341 - t367 * t329 - t368 * t338 + t360) * t236 + t204 + t205 + t203 - (Ifges(6,4) * t340 + Ifges(7,5) * t341 + t357 * t329 + t369 * t338 + t361) * t116 + (-t147 + t315) * t85 + t79 + t78 - t77 + t76 + (-t104 * t325 + t14 * t20 - t15 * t21 + (t217 * t4 + t3 * t300) * pkin(4)) * m(6) + t316 * t20 - t317 * t21 + t209 * t63 + t212 * t66 - t137 * (mrSges(5,1) * t184 + mrSges(5,2) * t183) + qJD(6) * t97 - t36 * t57 - t58 * t325 + t107 * t330 + (Ifges(5,1) * t183 - t306) * t331 + t64 * t324 + (Ifges(5,5) * t183 - Ifges(5,6) * t184) * t329 + t65 * t267 + (-Ifges(5,2) * t184 + t108 + t177) * t332 + t365 * t340; t316 * t236 + t317 * t116 + (-t10 * t236 + t116 * t12 + t5) * m(7) + (t116 * t15 + t14 * t236 + t69) * m(6) + t355; t236 * t57 - t206 * t97 + 0.2e1 * (t2 / 0.2e1 + t35 * t337 + t12 * t329) * m(7) + t66;];
tauc  = t6(:);
