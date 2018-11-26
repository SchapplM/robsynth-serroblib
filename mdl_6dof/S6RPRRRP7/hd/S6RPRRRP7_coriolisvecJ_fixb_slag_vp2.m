% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:27:45
% EndTime: 2018-11-23 16:27:57
% DurationCPUTime: 12.28s
% Computational Cost: add. (12085->604), mult. (31250->779), div. (0->0), fcn. (23503->8), ass. (0->265)
t369 = Ifges(7,5) - Ifges(6,4);
t379 = t369 + Ifges(7,5);
t218 = sin(pkin(10));
t219 = cos(pkin(10));
t222 = sin(qJ(3));
t225 = cos(qJ(3));
t197 = t218 * t225 + t219 * t222;
t191 = t197 * qJD(1);
t221 = sin(qJ(4));
t224 = cos(qJ(4));
t169 = qJD(3) * t221 + t191 * t224;
t220 = sin(qJ(5));
t223 = cos(qJ(5));
t232 = qJD(3) * t224 - t191 * t221;
t113 = t220 * t169 - t223 * t232;
t377 = -t218 * t222 + t225 * t219;
t192 = t377 * qJD(3);
t181 = qJD(1) * t192;
t124 = qJD(4) * t232 + t181 * t224;
t125 = -qJD(4) * t169 - t181 * t221;
t53 = -qJD(5) * t113 + t223 * t124 + t220 * t125;
t338 = t53 / 0.2e1;
t228 = t223 * t169 + t220 * t232;
t54 = qJD(5) * t228 + t220 * t124 - t223 * t125;
t336 = t54 / 0.2e1;
t193 = t197 * qJD(3);
t182 = qJD(1) * t193;
t321 = t182 / 0.2e1;
t371 = Ifges(6,1) + Ifges(7,1);
t370 = Ifges(7,4) + Ifges(6,5);
t368 = -Ifges(6,6) + Ifges(7,6);
t367 = Ifges(6,3) + Ifges(7,2);
t190 = t377 * qJD(1);
t276 = t220 * t224;
t201 = t221 * t223 + t276;
t136 = t201 * t190;
t352 = qJD(4) + qJD(5);
t161 = t352 * t201;
t271 = t136 - t161;
t277 = t220 * t221;
t200 = -t223 * t224 + t277;
t137 = t200 * t190;
t160 = t352 * t200;
t358 = t137 - t160;
t266 = qJD(4) * t224;
t283 = t192 * t221;
t378 = t197 * t266 + t283;
t306 = pkin(7) + qJ(2);
t205 = t306 * t218;
t198 = qJD(1) * t205;
t206 = t306 * t219;
t199 = qJD(1) * t206;
t157 = -t225 * t198 - t222 * t199;
t151 = -qJD(3) * pkin(3) - t157;
t106 = -pkin(4) * t232 + t151;
t185 = qJD(4) - t190;
t180 = qJD(5) + t185;
t259 = -pkin(2) * t219 - pkin(1);
t204 = qJD(1) * t259 + qJD(2);
t130 = -pkin(3) * t190 - pkin(8) * t191 + t204;
t158 = -t222 * t198 + t225 * t199;
t152 = qJD(3) * pkin(8) + t158;
t89 = t221 * t130 + t224 * t152;
t78 = pkin(9) * t232 + t89;
t290 = t220 * t78;
t88 = t224 * t130 - t152 * t221;
t77 = -pkin(9) * t169 + t88;
t66 = pkin(4) * t185 + t77;
t23 = t223 * t66 - t290;
t360 = qJD(6) - t23;
t19 = -pkin(5) * t180 + t360;
t110 = Ifges(6,4) * t113;
t291 = Ifges(7,5) * t113;
t365 = t180 * t370 + t228 * t371 - t110 + t291;
t47 = t113 * pkin(5) - qJ(6) * t228 + t106;
t376 = t106 * mrSges(6,2) + t19 * mrSges(7,2) - t23 * mrSges(6,3) - t47 * mrSges(7,3) + t365 / 0.2e1;
t73 = pkin(5) * t228 + qJ(6) * t113;
t230 = t197 * qJD(2);
t117 = qJD(1) * t230 + qJD(3) * t158;
t82 = -pkin(4) * t125 + t117;
t11 = pkin(5) * t54 - qJ(6) * t53 - qJD(6) * t228 + t82;
t229 = t377 * qJD(2);
t116 = qJD(1) * t229 + qJD(3) * t157;
t138 = pkin(3) * t182 - pkin(8) * t181;
t42 = -qJD(4) * t89 - t116 * t221 + t224 * t138;
t22 = pkin(4) * t182 - pkin(9) * t124 + t42;
t288 = t223 * t78;
t24 = t220 * t66 + t288;
t267 = qJD(4) * t221;
t41 = t224 * t116 + t130 * t266 + t221 * t138 - t152 * t267;
t26 = pkin(9) * t125 + t41;
t6 = -qJD(5) * t24 + t22 * t223 - t220 * t26;
t3 = -pkin(5) * t182 - t6;
t337 = -t54 / 0.2e1;
t375 = mrSges(6,2) * t82 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t11 + Ifges(6,4) * t337 + 0.2e1 * t370 * t321 + t336 * t379 + 0.2e1 * t338 * t371;
t374 = -t232 / 0.2e1;
t373 = t232 / 0.2e1;
t284 = t190 * t221;
t119 = pkin(4) * t284 + t158;
t364 = pkin(4) * t267 - pkin(5) * t271 - qJ(6) * t358 - qJD(6) * t201 - t119;
t363 = t232 * Ifges(5,6);
t362 = Ifges(4,5) * qJD(3);
t361 = Ifges(4,6) * qJD(3);
t359 = Ifges(5,5) * t124 + Ifges(5,6) * t125;
t155 = -pkin(3) * t377 - pkin(8) * t197 + t259;
t165 = -t205 * t222 + t206 * t225;
t159 = t224 * t165;
t105 = t221 * t155 + t159;
t357 = -t225 * t205 - t206 * t222;
t166 = Ifges(5,4) * t232;
t101 = t169 * Ifges(5,1) + t185 * Ifges(5,5) + t166;
t274 = t224 * t101;
t295 = Ifges(5,4) * t169;
t100 = t232 * Ifges(5,2) + Ifges(5,6) * t185 + t295;
t275 = t221 * t100;
t356 = t274 / 0.2e1 - t275 / 0.2e1;
t355 = t182 * t367 + t368 * t54 + t370 * t53;
t354 = -t221 * t42 + t224 * t41;
t353 = t42 * mrSges(5,1) - t41 * mrSges(5,2);
t351 = Ifges(5,5) * t169 + Ifges(5,3) * t185 + t113 * t368 + t180 * t367 + t228 * t370 + t363;
t350 = (m(3) * qJ(2) + mrSges(3,3)) * (t218 ^ 2 + t219 ^ 2);
t104 = t224 * t155 - t165 * t221;
t280 = t197 * t224;
t84 = -pkin(4) * t377 - pkin(9) * t280 + t104;
t281 = t197 * t221;
t90 = -pkin(9) * t281 + t105;
t303 = t220 * t84 + t223 * t90;
t131 = qJD(3) * t357 + t229;
t154 = pkin(3) * t193 - pkin(8) * t192;
t251 = -t131 * t221 + t224 * t154;
t282 = t192 * t224;
t32 = -pkin(9) * t282 + pkin(4) * t193 + (-t159 + (pkin(9) * t197 - t155) * t221) * qJD(4) + t251;
t55 = t224 * t131 + t221 * t154 + t155 * t266 - t165 * t267;
t38 = -pkin(9) * t378 + t55;
t10 = -qJD(5) * t303 - t220 * t38 + t223 * t32;
t349 = t197 * t352;
t348 = -t362 / 0.2e1 - t204 * mrSges(4,2);
t236 = t221 * t89 + t224 * t88;
t293 = Ifges(5,4) * t224;
t244 = -Ifges(5,2) * t221 + t293;
t247 = mrSges(5,1) * t221 + mrSges(5,2) * t224;
t347 = -t236 * mrSges(5,3) + t151 * t247 + t244 * t373;
t264 = qJD(5) * t223;
t265 = qJD(5) * t220;
t5 = t220 * t22 + t223 * t26 + t66 * t264 - t265 * t78;
t2 = qJ(6) * t182 + qJD(6) * t180 + t5;
t346 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3);
t20 = qJ(6) * t180 + t24;
t109 = Ifges(7,5) * t228;
t60 = t180 * Ifges(7,6) + t113 * Ifges(7,3) + t109;
t292 = Ifges(6,4) * t228;
t63 = -t113 * Ifges(6,2) + t180 * Ifges(6,6) + t292;
t344 = t106 * mrSges(6,1) + t47 * mrSges(7,1) + t60 / 0.2e1 - t63 / 0.2e1 - t20 * mrSges(7,2) - t24 * mrSges(6,3);
t342 = mrSges(6,1) * t82 + mrSges(7,1) * t11 - mrSges(7,2) * t2 - mrSges(6,3) * t5 + 0.2e1 * Ifges(7,3) * t336 - t53 * Ifges(6,4) / 0.2e1 - t182 * Ifges(6,6) / 0.2e1 + Ifges(7,6) * t321 + t379 * t338 + (-t337 + t336) * Ifges(6,2);
t341 = t19 * mrSges(7,1) + t24 * mrSges(6,2) + t89 * mrSges(5,2) + t361 / 0.2e1 - t20 * mrSges(7,3) - t204 * mrSges(4,1) - t23 * mrSges(6,1) - t363 / 0.2e1 - t88 * mrSges(5,1);
t334 = -pkin(9) - pkin(8);
t333 = -t113 / 0.2e1;
t332 = t113 / 0.2e1;
t330 = -t228 / 0.2e1;
t329 = t228 / 0.2e1;
t328 = t124 / 0.2e1;
t327 = t125 / 0.2e1;
t325 = -t169 / 0.2e1;
t324 = t169 / 0.2e1;
t323 = -t180 / 0.2e1;
t322 = t180 / 0.2e1;
t320 = -t185 / 0.2e1;
t319 = t185 / 0.2e1;
t318 = -t190 / 0.2e1;
t317 = t190 / 0.2e1;
t316 = -t191 / 0.2e1;
t315 = t191 / 0.2e1;
t311 = -t221 / 0.2e1;
t310 = t224 / 0.2e1;
t309 = m(6) * t106;
t43 = mrSges(6,1) * t182 - mrSges(6,3) * t53;
t44 = -t182 * mrSges(7,1) + t53 * mrSges(7,2);
t305 = -t43 + t44;
t45 = -mrSges(6,2) * t182 - mrSges(6,3) * t54;
t46 = -mrSges(7,2) * t54 + mrSges(7,3) * t182;
t304 = t45 + t46;
t153 = pkin(3) * t191 - pkin(8) * t190;
t97 = t224 * t153 - t157 * t221;
t79 = -pkin(9) * t190 * t224 + pkin(4) * t191 + t97;
t98 = t221 * t153 + t224 * t157;
t86 = -pkin(9) * t284 + t98;
t34 = t220 * t79 + t223 * t86;
t93 = -mrSges(7,2) * t113 + mrSges(7,3) * t180;
t298 = mrSges(6,3) * t113;
t94 = -mrSges(6,2) * t180 - t298;
t302 = t93 + t94;
t297 = mrSges(6,3) * t228;
t95 = mrSges(6,1) * t180 - t297;
t96 = -mrSges(7,1) * t180 + mrSges(7,2) * t228;
t301 = t96 - t95;
t300 = mrSges(4,3) * t190;
t299 = mrSges(4,3) * t191;
t296 = Ifges(4,4) * t191;
t294 = Ifges(5,4) * t221;
t286 = t117 * t357;
t269 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t232 - t169 * mrSges(5,2) - t299;
t260 = Ifges(5,3) * t182 + t359;
t215 = -pkin(4) * t224 - pkin(3);
t258 = qJD(4) * t334;
t252 = t182 * mrSges(4,1) + t181 * mrSges(4,2);
t250 = t224 * t258;
t133 = pkin(4) * t281 - t357;
t248 = mrSges(5,1) * t224 - mrSges(5,2) * t221;
t246 = Ifges(5,1) * t224 - t294;
t245 = Ifges(5,1) * t221 + t293;
t243 = Ifges(5,2) * t224 + t294;
t242 = Ifges(5,5) * t224 - Ifges(5,6) * t221;
t241 = Ifges(5,5) * t221 + Ifges(5,6) * t224;
t33 = -t220 * t86 + t223 * t79;
t39 = -t220 * t90 + t223 * t84;
t237 = -t221 * t41 - t224 * t42;
t235 = t221 * t88 - t224 * t89;
t128 = -t185 * mrSges(5,2) + mrSges(5,3) * t232;
t129 = mrSges(5,1) * t185 - mrSges(5,3) * t169;
t234 = t128 * t224 - t129 * t221;
t207 = t334 * t221;
t208 = t334 * t224;
t233 = t223 * t207 + t208 * t220;
t168 = t207 * t220 - t208 * t223;
t9 = t220 * t32 + t223 * t38 + t84 * t264 - t265 * t90;
t227 = t346 + t355;
t132 = qJD(3) * t165 + t230;
t91 = pkin(4) * t378 + t132;
t214 = -pkin(4) * t223 - pkin(5);
t212 = pkin(4) * t220 + qJ(6);
t209 = pkin(4) * t264 + qJD(6);
t203 = t221 * t258;
t184 = Ifges(4,4) * t190;
t170 = -qJD(3) * mrSges(4,2) + t300;
t156 = pkin(5) * t200 - qJ(6) * t201 + t215;
t147 = t200 * t197;
t146 = t201 * t197;
t141 = t191 * Ifges(4,1) + t184 + t362;
t140 = t190 * Ifges(4,2) + t296 + t361;
t108 = qJD(5) * t168 + t203 * t220 - t223 * t250;
t107 = qJD(5) * t233 + t223 * t203 + t220 * t250;
t103 = -mrSges(5,2) * t182 + mrSges(5,3) * t125;
t102 = mrSges(5,1) * t182 - mrSges(5,3) * t124;
t83 = -mrSges(5,1) * t125 + mrSges(5,2) * t124;
t75 = mrSges(6,1) * t113 + mrSges(6,2) * t228;
t74 = mrSges(7,1) * t113 - mrSges(7,3) * t228;
t71 = t124 * Ifges(5,1) + t125 * Ifges(5,4) + t182 * Ifges(5,5);
t70 = t124 * Ifges(5,4) + t125 * Ifges(5,2) + t182 * Ifges(5,6);
t69 = pkin(5) * t146 + qJ(6) * t147 + t133;
t68 = t192 * t276 + (t280 * t352 + t283) * t223 - t277 * t349;
t67 = -t200 * t192 - t201 * t349;
t58 = pkin(4) * t169 + t73;
t56 = -qJD(4) * t105 + t251;
t36 = pkin(5) * t377 - t39;
t35 = -qJ(6) * t377 + t303;
t30 = -pkin(5) * t191 - t33;
t29 = qJ(6) * t191 + t34;
t28 = t223 * t77 - t290;
t27 = t220 * t77 + t288;
t18 = mrSges(6,1) * t54 + mrSges(6,2) * t53;
t17 = mrSges(7,1) * t54 - mrSges(7,3) * t53;
t12 = pkin(5) * t68 - qJ(6) * t67 + qJD(6) * t147 + t91;
t8 = -pkin(5) * t193 - t10;
t7 = qJ(6) * t193 - qJD(6) * t377 + t9;
t1 = [-(mrSges(4,3) * t181 + t83) * t357 - (-mrSges(4,3) * t116 - Ifges(4,4) * t181 + Ifges(6,6) * t337 + Ifges(7,6) * t336 + t367 * t321 + t370 * t338 + t346 + t353) * t377 + m(5) * (t104 * t42 + t105 * t41 + t132 * t151 + t55 * t89 + t56 * t88 - t286) + m(4) * (t116 * t165 + t131 * t158 - t132 * t157 - t286) + (Ifges(5,3) * t319 + Ifges(5,5) * t324 + Ifges(7,6) * t332 + Ifges(6,6) * t333 - t158 * mrSges(4,3) - Ifges(4,4) * t315 - Ifges(4,2) * t317 - t140 / 0.2e1 - t341 + t370 * t329 + t367 * t322 + t351 / 0.2e1) * t193 + m(6) * (t10 * t23 + t106 * t91 + t133 * t82 + t24 * t9 + t303 * t5 + t39 * t6) + t303 * t45 + (t368 * t321 + t342) * t146 + (-Ifges(6,2) * t333 + Ifges(7,3) * t332 + t322 * t368 + t329 * t369 + t344) * t68 + (Ifges(4,1) * t181 + t71 * t310 + t70 * t311 + t246 * t328 + t244 * t327 + t242 * t321 + (mrSges(4,3) + t247) * t117 + (t151 * t248 + t243 * t374 + t241 * t320 + t245 * t325 + t101 * t311 - t224 * t100 / 0.2e1) * qJD(4)) * t197 + (t141 / 0.2e1 - t157 * mrSges(4,3) + Ifges(4,1) * t315 + Ifges(4,4) * t317 - t348 + t356) * t192 + (Ifges(5,5) * t282 - Ifges(5,6) * t283) * t319 + m(7) * (t11 * t69 + t12 * t47 + t19 * t8 + t2 * t35 + t20 * t7 + t3 * t36) + 0.2e1 * t350 * qJD(2) * qJD(1) + (Ifges(6,4) * t333 + Ifges(7,5) * t332 + t370 * t322 + t371 * t329 + t376) * t67 + (Ifges(5,4) * t282 - Ifges(5,2) * t283) * t373 - (t260 + t355 + t359) * t377 / 0.2e1 + (-(Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t377 - Ifges(4,4) * t197 - t165 * mrSges(4,3)) * t182 + ((qJD(4) * t235 + t237) * t197 - t89 * t283 - t88 * t282) * mrSges(5,3) + (Ifges(5,1) * t282 - Ifges(5,4) * t283) * t324 + t131 * t170 + t56 * t129 + t133 * t18 + t55 * t128 + t8 * t96 + t104 * t102 + t105 * t103 + t91 * t75 + t7 * t93 + t9 * t94 + t10 * t95 + t12 * t74 + t69 * t17 + t36 * t44 + t35 * t46 + t39 * t43 + t151 * (mrSges(5,1) * t283 + mrSges(5,2) * t282) - t269 * t132 - t375 * t147 + t259 * t252; t252 + t304 * t201 + t305 * t200 + (-t75 - t74 + t269) * t191 - m(4) * (-t157 * t191 + t158 * t190) + t221 * t103 + t224 * t102 + (-t170 - t234) * t190 + t234 * qJD(4) + t358 * t302 - t271 * t301 - t350 * qJD(1) ^ 2 + (-t19 * t271 - t191 * t47 + t2 * t201 + t20 * t358 + t200 * t3) * m(7) + (-t106 * t191 - t200 * t6 + t201 * t5 + t23 * t271 + t24 * t358) * m(6) + (-t151 * t191 - t185 * t235 - t237) * m(5); (t60 - t63) * (t161 / 0.2e1 - t136 / 0.2e1) + (Ifges(5,5) * t325 - Ifges(4,2) * t318 + Ifges(6,6) * t332 + Ifges(7,6) * t333 + Ifges(5,3) * t320 + t323 * t367 + t330 * t370 + t341) * t191 + (t136 * t368 - t137 * t370) * t323 + (t136 * t369 - t137 * t371) * t330 + (-t102 * t221 + t103 * t224 + (-m(5) * t236 - t221 * t128 - t224 * t129) * qJD(4) + m(5) * t354) * pkin(8) + t354 * mrSges(5,3) + (t242 * t319 + t246 * t324 + (t75 + t309) * t221 * pkin(4) + t347 + t356) * qJD(4) + (-pkin(3) * t117 - t151 * t158 - t88 * t97 - t89 * t98) * m(5) + (-t106 * t119 + t233 * t6 + t168 * t5 + t215 * t82 + (t107 - t34) * t24 + (-t108 - t33) * t23) * m(6) + (t11 * t156 - t233 * t3 + t168 * t2 + t364 * t47 + (t107 - t29) * t20 + (t108 - t30) * t19) * m(7) - t305 * t233 + (-Ifges(6,4) * t137 - Ifges(7,5) * t160 - Ifges(6,2) * t136 + Ifges(7,3) * t161) * t332 + (-Ifges(6,4) * t160 - Ifges(7,5) * t137 - Ifges(6,2) * t161 + Ifges(7,3) * t136) * t333 + (-t160 * t370 + t161 * t368) * t322 + (-t160 * t371 + t161 * t369) * t329 + t365 * (-t160 / 0.2e1 + t137 / 0.2e1) + t364 * t74 + (-t296 + t351) * t316 + (-t23 * t358 + t24 * t271) * mrSges(6,3) + (t19 * t358 + t20 * t271) * mrSges(7,2) + (-mrSges(7,1) * t271 - mrSges(7,3) * t358) * t47 + (-mrSges(6,1) * t271 + mrSges(6,2) * t358) * t106 + (t200 * t368 + t241) * t321 + t215 * t18 + (t184 + t274 + t141) * t318 + (Ifges(4,1) * t316 + t242 * t320 + t246 * t325 - t347 + t348) * t190 + t342 * t200 - Ifges(4,6) * t182 + Ifges(4,5) * t181 + t156 * t17 - t97 * t129 - t98 * t128 - t119 * t75 - t116 * mrSges(4,2) - t30 * t96 - t29 * t93 - t34 * t94 - t33 * t95 - pkin(3) * t83 + t304 * t168 + t301 * t108 + t302 * t107 + (t269 + t299) * t158 + (-t170 + t300) * t157 + t221 * t71 / 0.2e1 + t375 * t201 + t70 * t310 + t140 * t315 + t275 * t317 + (-mrSges(4,1) - t248) * t117 + t243 * t327 + t245 * t328; (-t19 * t27 + t2 * t212 + t214 * t3 - t47 * t58 + (t209 - t28) * t20) * m(7) + t227 + (-Ifges(6,2) * t332 + Ifges(7,3) * t333 + t368 * t323 + t369 * t330 - t344) * t228 + (m(6) * (t220 * t5 + t223 * t6 - t23 * t265 + t24 * t264) + 0.2e1 * t309 * t325 + m(7) * t19 * t265 - t169 * t75 + t220 * t45 + t223 * t43 + (t220 * t301 + t223 * t94) * qJD(5)) * pkin(4) + t209 * t93 + t212 * t46 + t214 * t44 + t260 - m(6) * (-t23 * t27 + t24 * t28) + t89 * t129 - t88 * t128 - t58 * t74 - t301 * t27 - t302 * t28 + t353 + (t89 * t169 + t232 * t88) * mrSges(5,3) - t151 * (t169 * mrSges(5,1) + mrSges(5,2) * t232) + (Ifges(5,5) * t232 - Ifges(5,6) * t169) * t320 + t100 * t324 + (Ifges(5,1) * t232 - t295) * t325 + (-Ifges(5,2) * t169 + t101 + t166) * t374 + (-Ifges(6,4) * t332 - Ifges(7,5) * t333 - t370 * t323 - t371 * t330 + t376) * t113; t227 + (t297 - t301) * t24 + (-t298 - t302) * t23 + (t113 * t19 + t20 * t228) * mrSges(7,2) - t47 * (mrSges(7,1) * t228 + mrSges(7,3) * t113) - t106 * (mrSges(6,1) * t228 - mrSges(6,2) * t113) + t63 * t329 + (Ifges(7,3) * t228 - t291) * t333 + qJD(6) * t93 - t73 * t74 - pkin(5) * t44 + qJ(6) * t46 + (-t113 * t370 + t228 * t368) * t323 + (-pkin(5) * t3 + qJ(6) * t2 - t19 * t24 + t20 * t360 - t47 * t73) * m(7) + (-Ifges(6,2) * t228 - t110 + t365) * t332 + (-t113 * t371 + t109 - t292 + t60) * t330; t228 * t74 - t180 * t93 + 0.2e1 * (t3 / 0.2e1 + t47 * t329 + t20 * t323) * m(7) + t44;];
tauc  = t1(:);
