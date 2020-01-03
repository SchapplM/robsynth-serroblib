% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:05
% EndTime: 2019-12-31 21:37:34
% DurationCPUTime: 13.15s
% Computational Cost: add. (9461->689), mult. (25589->999), div. (0->0), fcn. (19654->10), ass. (0->297)
t260 = sin(qJ(3));
t263 = cos(qJ(3));
t274 = pkin(3) * t260 - qJ(4) * t263;
t224 = qJD(3) * t274 - qJD(4) * t260;
t255 = sin(pkin(10));
t257 = cos(pkin(10));
t297 = qJD(3) * t260;
t294 = pkin(8) * t297;
t166 = t257 * t224 + t255 * t294;
t261 = sin(qJ(2));
t256 = sin(pkin(5));
t301 = qJD(1) * t256;
t290 = t261 * t301;
t258 = cos(pkin(5));
t264 = cos(qJ(2));
t347 = pkin(1) * t264;
t295 = t258 * t347;
t218 = -pkin(7) * t290 + qJD(1) * t295;
t271 = (pkin(2) * t261 - pkin(8) * t264) * t256;
t219 = qJD(1) * t271;
t148 = t263 * t218 + t260 * t219;
t130 = qJ(4) * t290 + t148;
t348 = pkin(1) * t261;
t252 = t258 * t348;
t311 = t256 * t264;
t151 = (t252 + (pkin(7) + t274) * t311) * qJD(1);
t78 = -t255 * t130 + t257 * t151;
t391 = -t78 + t166;
t308 = t263 * t264;
t183 = (t255 * t261 + t257 * t308) * t301;
t289 = t264 * t301;
t309 = t257 * t263;
t346 = pkin(4) * t260;
t390 = t183 * pkin(9) - t289 * t346 + (-pkin(9) * t309 + t346) * qJD(3) + t391;
t182 = (-t255 * t308 + t257 * t261) * t301;
t209 = t255 * t224;
t310 = t257 * t260;
t313 = t255 * t263;
t79 = t257 * t130 + t255 * t151;
t389 = -pkin(9) * t182 + t209 + (-pkin(8) * t310 - pkin(9) * t313) * qJD(3) - t79;
t239 = -pkin(3) * t263 - qJ(4) * t260 - pkin(2);
t233 = t257 * t239;
t164 = -pkin(9) * t310 + t233 + (-pkin(8) * t255 - pkin(4)) * t263;
t197 = pkin(8) * t309 + t255 * t239;
t173 = -pkin(9) * t255 * t260 + t197;
t259 = sin(qJ(5));
t262 = cos(qJ(5));
t104 = t164 * t259 + t173 * t262;
t388 = -qJD(5) * t104 - t389 * t259 + t390 * t262;
t103 = t164 * t262 - t173 * t259;
t387 = qJD(5) * t103 + t390 * t259 + t389 * t262;
t147 = -t260 * t218 + t219 * t263;
t131 = -pkin(3) * t290 - t147;
t291 = pkin(4) * t255 + pkin(8);
t296 = qJD(3) * t263;
t386 = pkin(4) * t182 + t291 * t296 - t131;
t283 = qJD(1) * t258 + qJD(2);
t379 = -t260 * t290 + t263 * t283;
t191 = Ifges(4,4) * t379;
t244 = qJD(3) - t289;
t319 = t244 * Ifges(4,5);
t199 = t260 * t283 + t263 * t290;
t321 = t199 * Ifges(4,1);
t113 = t191 + t319 + t321;
t251 = pkin(7) * t311;
t181 = qJD(2) * pkin(8) + (t251 + (pkin(8) + t348) * t258) * qJD(1);
t213 = (-pkin(2) * t264 - pkin(8) * t261 - pkin(1)) * t256;
t190 = qJD(1) * t213;
t116 = -t260 * t181 + t190 * t263;
t180 = -pkin(2) * t283 - t218;
t385 = -t191 / 0.2e1 + t116 * mrSges(4,3) - t113 / 0.2e1 - t180 * mrSges(4,2) - t319 / 0.2e1;
t154 = t199 * t257 + t244 * t255;
t284 = -t199 * t255 + t257 * t244;
t384 = -t154 * t259 + t262 * t284;
t91 = t154 * t262 + t259 * t284;
t300 = qJD(2) * t256;
t285 = qJD(1) * t300;
t280 = t264 * t285;
t160 = qJD(3) * t379 + t263 * t280;
t281 = t261 * t285;
t128 = -t160 * t255 + t257 * t281;
t129 = t160 * t257 + t255 * t281;
t34 = qJD(5) * t384 + t128 * t259 + t129 * t262;
t374 = t34 / 0.2e1;
t35 = -qJD(5) * t91 + t128 * t262 - t129 * t259;
t373 = t35 / 0.2e1;
t161 = qJD(3) * t199 + t260 * t280;
t361 = t161 / 0.2e1;
t286 = -Ifges(3,6) * qJD(2) / 0.2e1;
t339 = pkin(9) + qJ(4);
t240 = t339 * t255;
t241 = t339 * t257;
t177 = -t240 * t259 + t241 * t262;
t236 = t255 * t262 + t257 * t259;
t136 = pkin(3) * t199 - qJ(4) * t379;
t65 = -t116 * t255 + t257 * t136;
t42 = -pkin(9) * t257 * t379 + pkin(4) * t199 + t65;
t314 = t379 * t255;
t66 = t257 * t116 + t255 * t136;
t52 = -pkin(9) * t314 + t66;
t383 = -qJD(4) * t236 - qJD(5) * t177 + t259 * t52 - t262 * t42;
t176 = -t240 * t262 - t241 * t259;
t272 = t255 * t259 - t257 * t262;
t382 = -qJD(4) * t272 + qJD(5) * t176 - t259 * t42 - t262 * t52;
t312 = t256 * t261;
t248 = pkin(7) * t312;
t229 = -t248 + t295;
t225 = t272 * qJD(5);
t318 = t244 * Ifges(4,6);
t337 = Ifges(4,4) * t199;
t112 = Ifges(4,2) * t379 + t318 + t337;
t117 = t181 * t263 + t190 * t260;
t192 = qJD(5) - t379;
t323 = t192 * Ifges(6,3);
t341 = t91 * Ifges(6,5);
t342 = t384 * Ifges(6,6);
t31 = t323 + t341 + t342;
t327 = t154 * Ifges(5,5);
t328 = t284 * Ifges(5,6);
t100 = -pkin(3) * t379 - t199 * qJ(4) + t180;
t105 = qJ(4) * t244 + t117;
t45 = t257 * t100 - t105 * t255;
t46 = t255 * t100 + t257 * t105;
t73 = -Ifges(5,3) * t379 + t327 + t328;
t25 = -pkin(4) * t379 - pkin(9) * t154 + t45;
t36 = pkin(9) * t284 + t46;
t8 = t25 * t262 - t259 * t36;
t9 = t25 * t259 + t262 * t36;
t265 = -t180 * mrSges(4,1) + t9 * mrSges(6,2) - t8 * mrSges(6,1) - t45 * mrSges(5,1) + t46 * mrSges(5,2) + t117 * mrSges(4,3) - t328 / 0.2e1 - t327 / 0.2e1 + t318 / 0.2e1 - t31 / 0.2e1 - t73 / 0.2e1 - t342 / 0.2e1 - t341 / 0.2e1 - t323 / 0.2e1 + t112 / 0.2e1 + t337 / 0.2e1;
t293 = Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1;
t378 = t293 * t379 + t265;
t102 = -pkin(3) * t244 + qJD(4) - t116;
t335 = Ifges(5,4) * t257;
t276 = -Ifges(5,2) * t255 + t335;
t336 = Ifges(5,4) * t255;
t277 = Ifges(5,1) * t257 - t336;
t278 = mrSges(5,1) * t255 + mrSges(5,2) * t257;
t350 = t257 / 0.2e1;
t351 = -t255 / 0.2e1;
t362 = t154 / 0.2e1;
t363 = t284 / 0.2e1;
t74 = t154 * Ifges(5,4) + Ifges(5,2) * t284 - Ifges(5,6) * t379;
t75 = t154 * Ifges(5,1) + Ifges(5,4) * t284 - Ifges(5,5) * t379;
t377 = t102 * t278 + t276 * t363 + t277 * t362 + (-t46 * t255 - t45 * t257) * mrSges(5,3) + t74 * t351 + t75 * t350 - t385;
t376 = Ifges(6,4) * t374 + Ifges(6,2) * t373 + Ifges(6,6) * t361;
t375 = Ifges(6,1) * t374 + Ifges(6,4) * t373 + Ifges(6,5) * t361;
t50 = t129 * Ifges(5,1) + t128 * Ifges(5,4) + t161 * Ifges(5,5);
t372 = t50 / 0.2e1;
t371 = -t384 / 0.2e1;
t370 = t384 / 0.2e1;
t369 = -t91 / 0.2e1;
t368 = t91 / 0.2e1;
t367 = pkin(1) * mrSges(3,1);
t366 = pkin(1) * mrSges(3,2);
t365 = t128 / 0.2e1;
t364 = t129 / 0.2e1;
t358 = -t192 / 0.2e1;
t357 = t192 / 0.2e1;
t356 = t379 / 0.2e1;
t355 = -t379 / 0.2e1;
t227 = -t258 * t263 + t260 * t312;
t354 = -t227 / 0.2e1;
t228 = t258 * t260 + t263 * t312;
t352 = t228 / 0.2e1;
t349 = Ifges(6,4) * t91;
t30 = Ifges(6,5) * t34;
t29 = Ifges(6,6) * t35;
t345 = pkin(8) * t260;
t220 = qJD(2) * t271;
t206 = qJD(1) * t220;
t222 = t229 * qJD(2);
t207 = qJD(1) * t222;
t70 = -t181 * t297 + t190 * t296 + t260 * t206 + t263 * t207;
t344 = t70 * mrSges(4,2);
t71 = -t181 * t296 - t190 * t297 + t206 * t263 - t260 * t207;
t343 = t71 * mrSges(4,1);
t340 = qJD(2) / 0.2e1;
t55 = qJ(4) * t281 + qJD(4) * t244 + t70;
t302 = t251 + t252;
t223 = t302 * qJD(2);
t208 = qJD(1) * t223;
t59 = t161 * pkin(3) - t160 * qJ(4) - t199 * qJD(4) + t208;
t22 = t255 * t59 + t257 * t55;
t299 = qJD(2) * t261;
t212 = pkin(8) * t258 + t302;
t84 = -t212 * t297 + t213 * t296 + t260 * t220 + t263 * t222;
t72 = (qJ(4) * t299 - qJD(4) * t264) * t256 + t84;
t298 = qJD(2) * t264;
t287 = t256 * t298;
t171 = qJD(3) * t228 + t260 * t287;
t172 = -qJD(3) * t227 + t263 * t287;
t81 = t171 * pkin(3) - t172 * qJ(4) - t228 * qJD(4) + t223;
t27 = t255 * t81 + t257 * t72;
t338 = Ifges(3,4) * t261;
t334 = Ifges(3,5) * t258;
t333 = Ifges(5,5) * t257;
t332 = Ifges(3,6) * t258;
t331 = Ifges(5,6) * t255;
t330 = t128 * Ifges(5,6);
t329 = t129 * Ifges(5,5);
t326 = t160 * Ifges(4,1);
t325 = t160 * Ifges(4,4);
t324 = t161 * Ifges(4,4);
t163 = mrSges(4,1) * t244 - mrSges(4,3) * t199;
t95 = -mrSges(5,1) * t284 + mrSges(5,2) * t154;
t316 = -t95 + t163;
t211 = t248 + (-pkin(2) - t347) * t258;
t124 = t227 * pkin(3) - t228 * qJ(4) + t211;
t139 = t263 * t212 + t260 * t213;
t126 = -qJ(4) * t311 + t139;
t62 = t255 * t124 + t257 * t126;
t119 = t182 * t262 - t183 * t259;
t150 = t225 * t260 - t236 * t296;
t307 = t119 - t150;
t120 = t182 * t259 + t183 * t262;
t226 = t236 * qJD(5);
t149 = -t226 * t260 - t272 * t296;
t306 = t120 - t149;
t121 = t236 * t379;
t305 = -t121 + t226;
t122 = t272 * t379;
t304 = -t122 + t225;
t303 = -mrSges(3,1) * t283 - mrSges(4,1) * t379 + mrSges(4,2) * t199 + mrSges(3,3) * t290;
t5 = Ifges(6,3) * t161 + t29 + t30;
t292 = Ifges(4,5) * t160 - Ifges(4,6) * t161 + Ifges(4,3) * t281;
t288 = t256 * t299;
t10 = -t35 * mrSges(6,1) + t34 * mrSges(6,2);
t21 = -t255 * t55 + t257 * t59;
t26 = -t255 * t72 + t257 * t81;
t69 = -t128 * mrSges(5,1) + t129 * mrSges(5,2);
t61 = t257 * t124 - t126 * t255;
t138 = -t260 * t212 + t213 * t263;
t11 = pkin(4) * t161 - pkin(9) * t129 + t21;
t14 = pkin(9) * t128 + t22;
t1 = qJD(5) * t8 + t11 * t259 + t14 * t262;
t2 = -qJD(5) * t9 + t11 * t262 - t14 * t259;
t279 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t127 = pkin(3) * t311 - t138;
t275 = -t331 + t333;
t273 = -t21 * t255 + t22 * t257;
t170 = t228 * t257 - t255 * t311;
t39 = pkin(4) * t227 - pkin(9) * t170 + t61;
t169 = -t228 * t255 - t257 * t311;
t47 = pkin(9) * t169 + t62;
t12 = -t259 * t47 + t262 * t39;
t13 = t259 * t39 + t262 * t47;
t106 = t169 * t262 - t170 * t259;
t107 = t169 * t259 + t170 * t262;
t85 = -t212 * t296 - t213 * t297 + t220 * t263 - t260 * t222;
t245 = Ifges(3,4) * t289;
t269 = -t218 * mrSges(3,3) + Ifges(3,1) * t290 / 0.2e1 + t245 / 0.2e1 + (t283 / 0.2e1 + t340) * Ifges(3,5);
t80 = -pkin(3) * t288 - t85;
t60 = -pkin(3) * t281 - t71;
t221 = t302 * qJD(1);
t267 = t116 * mrSges(4,1) + t244 * Ifges(4,3) + t199 * Ifges(4,5) + t379 * Ifges(4,6) + t286 - (t332 + (Ifges(3,2) * t264 + t338) * t256) * qJD(1) / 0.2e1 - t117 * mrSges(4,2) - t221 * mrSges(3,3);
t254 = -pkin(4) * t257 - pkin(3);
t243 = Ifges(3,5) * t280;
t237 = t291 * t260;
t217 = -mrSges(3,2) * t283 + mrSges(3,3) * t289;
t215 = t272 * t260;
t214 = t236 * t260;
t196 = -pkin(8) * t313 + t233;
t167 = -t257 * t294 + t209;
t162 = -mrSges(4,2) * t244 + mrSges(4,3) * t379;
t146 = t172 * t257 + t255 * t288;
t145 = -t172 * t255 + t257 * t288;
t135 = -mrSges(4,2) * t281 - mrSges(4,3) * t161;
t134 = mrSges(4,1) * t281 - mrSges(4,3) * t160;
t109 = -mrSges(5,1) * t379 - mrSges(5,3) * t154;
t108 = mrSges(5,2) * t379 + mrSges(5,3) * t284;
t99 = mrSges(4,1) * t161 + mrSges(4,2) * t160;
t94 = pkin(4) * t314 + t117;
t93 = -pkin(4) * t169 + t127;
t88 = Ifges(4,5) * t281 - t324 + t326;
t87 = -t161 * Ifges(4,2) + Ifges(4,6) * t281 + t325;
t86 = Ifges(6,4) * t384;
t83 = mrSges(5,1) * t161 - mrSges(5,3) * t129;
t82 = -mrSges(5,2) * t161 + mrSges(5,3) * t128;
t68 = -pkin(4) * t284 + t102;
t64 = mrSges(6,1) * t192 - mrSges(6,3) * t91;
t63 = -mrSges(6,2) * t192 + mrSges(6,3) * t384;
t51 = -pkin(4) * t145 + t80;
t49 = t129 * Ifges(5,4) + t128 * Ifges(5,2) + t161 * Ifges(5,6);
t48 = t161 * Ifges(5,3) + t329 + t330;
t41 = -qJD(5) * t107 + t145 * t262 - t146 * t259;
t40 = qJD(5) * t106 + t145 * t259 + t146 * t262;
t38 = -pkin(4) * t128 + t60;
t37 = -mrSges(6,1) * t384 + mrSges(6,2) * t91;
t33 = Ifges(6,1) * t91 + Ifges(6,5) * t192 + t86;
t32 = Ifges(6,2) * t384 + Ifges(6,6) * t192 + t349;
t24 = -mrSges(6,2) * t161 + mrSges(6,3) * t35;
t23 = mrSges(6,1) * t161 - mrSges(6,3) * t34;
t20 = pkin(9) * t145 + t27;
t17 = pkin(4) * t171 - pkin(9) * t146 + t26;
t4 = -qJD(5) * t13 + t17 * t262 - t20 * t259;
t3 = qJD(5) * t12 + t17 * t259 + t20 * t262;
t6 = [(t48 + t5) * t227 / 0.2e1 + (t73 + t31) * t171 / 0.2e1 + m(3) * (t207 * t302 - t208 * t229 - t218 * t223 + t221 * t222) + ((-t229 * mrSges(3,3) + t334 + (-0.2e1 * t366 + 0.3e1 / 0.2e1 * Ifges(3,4) * t264) * t256) * t298 + (-t302 * mrSges(3,3) - 0.3e1 / 0.2e1 * t332 + Ifges(4,5) * t352 + Ifges(4,6) * t354 + (-0.2e1 * t367 - 0.3e1 / 0.2e1 * t338 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1) * t264) * t256) * t299) * t301 + (Ifges(5,5) * t170 + Ifges(6,5) * t107 + Ifges(5,6) * t169 + Ifges(6,6) * t106 + (Ifges(5,3) + Ifges(6,3)) * t227) * t361 + t311 * t344 + (t269 * t264 + (t286 + t267) * t261) * t300 + (-t116 * t172 - t117 * t171 - t227 * t70 - t228 * t71) * mrSges(4,3) + m(6) * (t1 * t13 + t12 * t2 + t3 * t9 + t38 * t93 + t4 * t8 + t51 * t68) + m(5) * (t102 * t80 + t127 * t60 + t21 * t61 + t22 * t62 + t26 * t45 + t27 * t46) + m(4) * (t116 * t85 + t117 * t84 + t138 * t71 + t139 * t70 + t180 * t223 + t208 * t211) + t244 * (Ifges(4,5) * t172 - Ifges(4,6) * t171) / 0.2e1 + t199 * (Ifges(4,1) * t172 - Ifges(4,4) * t171) / 0.2e1 - t311 * t343 + t88 * t352 + t87 * t354 + (Ifges(5,5) * t146 + Ifges(5,6) * t145 + Ifges(5,3) * t171) * t355 + (Ifges(4,4) * t172 - Ifges(4,2) * t171) * t356 + (Ifges(6,5) * t40 + Ifges(6,6) * t41 + Ifges(6,3) * t171) * t357 + (Ifges(5,1) * t146 + Ifges(5,4) * t145 + Ifges(5,5) * t171) * t362 + (Ifges(5,4) * t146 + Ifges(5,2) * t145 + Ifges(5,6) * t171) * t363 + (Ifges(5,1) * t170 + Ifges(5,4) * t169 + Ifges(5,5) * t227) * t364 + (Ifges(5,4) * t170 + Ifges(5,2) * t169 + Ifges(5,6) * t227) * t365 + (Ifges(6,1) * t40 + Ifges(6,4) * t41 + Ifges(6,5) * t171) * t368 + (Ifges(6,4) * t40 + Ifges(6,2) * t41 + Ifges(6,6) * t171) * t370 + t170 * t372 + (Ifges(6,4) * t107 + Ifges(6,2) * t106 + Ifges(6,6) * t227) * t373 + (Ifges(6,1) * t107 + Ifges(6,4) * t106 + Ifges(6,5) * t227) * t374 + t107 * t375 + t106 * t376 + t258 * t243 / 0.2e1 + t12 * t23 + t13 * t24 + t40 * t33 / 0.2e1 + t41 * t32 / 0.2e1 + t51 * t37 + t3 * t63 + t4 * t64 + t68 * (-mrSges(6,1) * t41 + mrSges(6,2) * t40) + t62 * t82 + t61 * t83 + t93 * t10 + t80 * t95 + t38 * (-mrSges(6,1) * t106 + mrSges(6,2) * t107) + t27 * t108 + t26 * t109 + t127 * t69 + t138 * t134 + t139 * t135 + t145 * t74 / 0.2e1 + t146 * t75 / 0.2e1 + t102 * (-mrSges(5,1) * t145 + mrSges(5,2) * t146) + t84 * t162 + t85 * t163 + t169 * t49 / 0.2e1 + t60 * (-mrSges(5,1) * t169 + mrSges(5,2) * t170) + t8 * (mrSges(6,1) * t171 - mrSges(6,3) * t40) + t9 * (-mrSges(6,2) * t171 + mrSges(6,3) * t41) + t46 * (-mrSges(5,2) * t171 + mrSges(5,3) * t145) + t45 * (mrSges(5,1) * t171 - mrSges(5,3) * t146) - t171 * t112 / 0.2e1 + t172 * t113 / 0.2e1 + t180 * (mrSges(4,1) * t171 + mrSges(4,2) * t172) + t211 * t99 + t303 * t223 + t222 * t217 + t1 * (-mrSges(6,2) * t227 + mrSges(6,3) * t106) + t2 * (mrSges(6,1) * t227 - mrSges(6,3) * t107) + t22 * (-mrSges(5,2) * t227 + mrSges(5,3) * t169) + t21 * (mrSges(5,1) * t227 - mrSges(5,3) * t170) - t292 * t311 / 0.2e1 - t161 * (Ifges(4,4) * t228 - Ifges(4,2) * t227 - Ifges(4,6) * t311) / 0.2e1 + t160 * (Ifges(4,1) * t228 - Ifges(4,4) * t227 - Ifges(4,5) * t311) / 0.2e1 + t207 * (-t258 * mrSges(3,2) + mrSges(3,3) * t311) + (-mrSges(3,1) * t258 + mrSges(4,1) * t227 + mrSges(4,2) * t228 + mrSges(3,3) * t312) * t208; -t284 * (Ifges(5,4) * t183 + Ifges(5,2) * t182) / 0.2e1 + (-Ifges(6,5) * t215 - Ifges(6,6) * t214) * t361 + (-Ifges(6,4) * t215 - Ifges(6,2) * t214) * t373 + (-Ifges(6,1) * t215 - Ifges(6,4) * t214) * t374 + t38 * (mrSges(6,1) * t214 - mrSges(6,2) * t215) + (-t1 * t214 + t2 * t215 + t306 * t8 - t307 * t9) * mrSges(6,3) + t391 * t109 + t386 * t37 + ((t286 + (Ifges(4,5) * t260 + Ifges(4,6) * t263) * t340 + (t332 / 0.2e1 + (t367 + t338 / 0.2e1) * t256) * qJD(1) - t267) * t261 + (-t245 / 0.2e1 + (-t334 / 0.2e1 + (t366 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t261) * t256) * qJD(1) + (-t321 / 0.2e1 + t385) * t263 + t378 * t260 - t269) * t264) * t301 - t154 * (Ifges(5,1) * t183 + Ifges(5,4) * t182) / 0.2e1 + (((-m(4) * t117 - t162) * pkin(8) - t378) * t260 + (t275 * t355 + t321 / 0.2e1 + (-m(4) * t116 + m(5) * t102 - t316) * pkin(8) + t377) * t263) * qJD(3) + t243 + (-t119 / 0.2e1 + t150 / 0.2e1) * t32 + (-t120 / 0.2e1 + t149 / 0.2e1) * t33 + (-t182 * t46 + t183 * t45) * mrSges(5,3) + (-t79 + t167) * t108 - m(5) * (t102 * t131 + t45 * t78 + t46 * t79) - m(4) * (t116 * t147 + t117 * t148 + t180 * t221) + (t275 * t361 + t276 * t365 + t277 * t364 + t326 / 0.2e1 - t324 / 0.2e1 + t208 * mrSges(4,2) + t88 / 0.2e1 + t60 * t278 - t71 * mrSges(4,3) + t49 * t351 + t50 * t350 + (t69 - t134) * pkin(8) + (-t21 * t257 - t22 * t255) * mrSges(5,3)) * t260 + (t22 * mrSges(5,2) - t21 * mrSges(5,1) - t330 / 0.2e1 - t329 / 0.2e1 - t5 / 0.2e1 - t48 / 0.2e1 + t87 / 0.2e1 - t30 / 0.2e1 - t29 / 0.2e1 + t325 / 0.2e1 - t208 * mrSges(4,1) + t70 * mrSges(4,3) + pkin(8) * t135 + (-Ifges(6,3) / 0.2e1 - t293) * t161 - t279) * t263 + m(4) * (pkin(8) * t263 * t70 - pkin(2) * t208 - t345 * t71) + m(5) * (t166 * t45 + t167 * t46 + t196 * t21 + t197 * t22 + t345 * t60) + t387 * t63 + t388 * t64 + (t1 * t104 + t103 * t2 + t237 * t38 + t386 * t68 + t387 * t9 + t388 * t8) * m(6) + (Ifges(5,5) * t183 + Ifges(5,6) * t182) * t356 + (Ifges(6,5) * t149 + Ifges(6,6) * t150) * t357 + (Ifges(6,5) * t120 + Ifges(6,6) * t119) * t358 + (Ifges(6,1) * t149 + Ifges(6,4) * t150) * t368 + (Ifges(6,1) * t120 + Ifges(6,4) * t119) * t369 + (Ifges(6,4) * t149 + Ifges(6,2) * t150) * t370 + (Ifges(6,4) * t120 + Ifges(6,2) * t119) * t371 - t215 * t375 - t214 * t376 - pkin(2) * t99 + t103 * t23 + t104 * t24 - t131 * t95 - t148 * t162 - t147 * t163 - t182 * t74 / 0.2e1 - t183 * t75 / 0.2e1 - t102 * (-mrSges(5,1) * t182 + mrSges(5,2) * t183) + t196 * t83 + t197 * t82 - t207 * mrSges(3,2) - t208 * mrSges(3,1) - t303 * t221 - t218 * t217 + (mrSges(6,1) * t307 - mrSges(6,2) * t306) * t68 + t237 * t10; (-Ifges(6,5) * t225 - Ifges(6,6) * t226) * t357 + (-Ifges(6,1) * t225 - Ifges(6,4) * t226) * t368 + (-Ifges(6,4) * t225 - Ifges(6,2) * t226) * t370 - (-(t333 / 0.2e1 - t331 / 0.2e1) * t379 + (Ifges(4,1) / 0.2e1 - t293) * t199 + t377) * t379 + (t121 / 0.2e1 - t226 / 0.2e1) * t32 + (t122 / 0.2e1 - t225 / 0.2e1) * t33 + (-Ifges(6,5) * t122 - Ifges(6,6) * t121) * t358 + (-Ifges(6,1) * t122 - Ifges(6,4) * t121) * t369 + (-Ifges(6,4) * t122 - Ifges(6,2) * t121) * t371 + (-t255 * t83 + t257 * t82) * qJ(4) - t344 + t343 + t382 * t63 + (t1 * t177 + t176 * t2 + t254 * t38 + t382 * t9 + t383 * t8 - t68 * t94) * m(6) + t383 * t64 + (-t102 * t117 - t45 * t65 - t46 * t66 - pkin(3) * t60 + (-t255 * t45 + t257 * t46) * qJD(4) + t273 * qJ(4)) * m(5) + (t108 * t257 - t109 * t255) * qJD(4) + t60 * (-mrSges(5,1) * t257 + mrSges(5,2) * t255) + t273 * mrSges(5,3) + t49 * t350 + (Ifges(5,1) * t255 + t335) * t364 + (Ifges(5,2) * t257 + t336) * t365 + t255 * t372 + t236 * t375 + t265 * t199 + t292 - pkin(3) * t69 - t94 * t37 + (Ifges(5,5) * t255 + Ifges(6,5) * t236 + Ifges(5,6) * t257 - Ifges(6,6) * t272) * t361 + (Ifges(6,4) * t236 - Ifges(6,2) * t272) * t373 + (Ifges(6,1) * t236 - Ifges(6,4) * t272) * t374 + (-t1 * t272 - t2 * t236 + t304 * t8 - t305 * t9) * mrSges(6,3) + t38 * (mrSges(6,1) * t272 + mrSges(6,2) * t236) - t272 * t376 - t66 * t108 - t65 * t109 - t116 * t162 + t176 * t23 + t177 * t24 + (mrSges(6,1) * t305 - mrSges(6,2) * t304) * t68 + t254 * t10 + t316 * t117; -t284 * t108 + t154 * t109 - t384 * t63 + t91 * t64 + t10 + t69 + (-t384 * t9 + t8 * t91 + t38) * m(6) + (t154 * t45 - t284 * t46 + t60) * m(5); -t68 * (mrSges(6,1) * t91 + mrSges(6,2) * t384) + (Ifges(6,1) * t384 - t349) * t369 + t32 * t368 + (Ifges(6,5) * t384 - Ifges(6,6) * t91) * t358 + t9 * t64 - t8 * t63 + (t384 * t8 + t9 * t91) * mrSges(6,3) + t279 + t5 + (-Ifges(6,2) * t91 + t33 + t86) * t371;];
tauc = t6(:);
