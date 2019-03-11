% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:22
% EndTime: 2019-03-09 04:46:42
% DurationCPUTime: 10.54s
% Computational Cost: add. (5399->544), mult. (12368->733), div. (0->0), fcn. (7592->6), ass. (0->243)
t197 = sin(qJ(3));
t199 = cos(qJ(3));
t225 = pkin(3) * t199 + pkin(8) * t197;
t169 = t225 * qJD(1);
t200 = -pkin(1) - pkin(7);
t182 = qJD(1) * t200 + qJD(2);
t196 = sin(qJ(4));
t198 = cos(qJ(4));
t264 = t198 * t199;
t111 = t196 * t169 + t182 * t264;
t292 = -qJ(5) - pkin(8);
t227 = qJD(4) * t292;
t260 = qJD(1) * t197;
t243 = t196 * t260;
t251 = qJD(5) * t198;
t364 = -qJ(5) * t243 + t196 * t227 - t111 + t251;
t266 = t196 * t199;
t110 = t198 * t169 - t182 * t266;
t363 = -qJD(5) * t196 + t198 * t227 - (qJ(5) * t197 * t198 + pkin(4) * t199) * qJD(1) - t110;
t348 = Ifges(6,1) + Ifges(7,1);
t347 = Ifges(6,4) - Ifges(7,5);
t346 = Ifges(7,4) + Ifges(6,5);
t194 = sin(pkin(9));
t195 = cos(pkin(9));
t162 = t194 * t198 + t195 * t196;
t145 = t162 * qJD(1);
t127 = t197 * t145;
t337 = t162 * qJD(4);
t263 = t127 + t337;
t161 = t194 * t196 - t195 * t198;
t339 = t161 * t197;
t128 = qJD(1) * t339;
t146 = t161 * qJD(4);
t262 = t128 + t146;
t345 = Ifges(6,6) - Ifges(7,6);
t356 = -t364 * t194 + t363 * t195;
t355 = t363 * t194 + t364 * t195;
t269 = t182 * t199;
t154 = -qJD(3) * pkin(3) - t269;
t257 = qJD(3) * t198;
t259 = qJD(1) * t199;
t165 = -t196 * t259 + t257;
t171 = pkin(3) * t197 - pkin(8) * t199 + qJ(2);
t150 = t171 * qJD(1);
t170 = t197 * t182;
t153 = qJD(3) * pkin(8) + t170;
t97 = t198 * t150 - t153 * t196;
t98 = t150 * t196 + t153 * t198;
t215 = t196 * t98 + t198 * t97;
t282 = Ifges(5,4) * t198;
t220 = -Ifges(5,2) * t196 + t282;
t283 = Ifges(5,4) * t196;
t222 = Ifges(5,1) * t198 - t283;
t223 = mrSges(5,1) * t196 + mrSges(5,2) * t198;
t280 = Ifges(5,6) * t196;
t281 = Ifges(5,5) * t198;
t298 = t198 / 0.2e1;
t300 = -t196 / 0.2e1;
t188 = qJD(4) + t260;
t301 = t188 / 0.2e1;
t166 = qJD(3) * t196 + t198 * t259;
t303 = t166 / 0.2e1;
t284 = Ifges(5,4) * t166;
t94 = t165 * Ifges(5,2) + t188 * Ifges(5,6) + t284;
t156 = Ifges(5,4) * t165;
t95 = t166 * Ifges(5,1) + t188 * Ifges(5,5) + t156;
t362 = -t215 * mrSges(5,3) + t154 * t223 + (-t280 + t281) * t301 + t220 * t165 / 0.2e1 + t222 * t303 + t298 * t95 + t300 * t94;
t361 = Ifges(6,3) + Ifges(7,2) + Ifges(5,3);
t360 = qJD(1) / 0.2e1;
t256 = qJD(3) * t199;
t229 = qJD(1) * t256;
t252 = qJD(4) * t199;
t236 = t196 * t252;
t240 = t197 * t257;
t250 = qJD(3) * qJD(4);
t118 = t198 * t250 + (-t236 - t240) * qJD(1);
t235 = t198 * t252;
t258 = qJD(3) * t197;
t352 = -t196 * t258 + t235;
t119 = -qJD(1) * t352 - t196 * t250;
t68 = t118 * t194 - t195 * t119;
t69 = t118 * t195 + t119 * t194;
t359 = t346 * t229 - t347 * t68 + t348 * t69;
t358 = -qJ(6) * t259 + t355;
t357 = pkin(5) * t259 - t356;
t104 = -t195 * t165 + t166 * t194;
t211 = t165 * t194 + t195 * t166;
t343 = -t347 * t104 + t346 * t188 + t348 * t211;
t130 = -pkin(4) * t243 + t170;
t254 = qJD(4) * t196;
t354 = pkin(4) * t254 + t263 * pkin(5) + t262 * qJ(6) - qJD(6) * t162 - t130;
t272 = Ifges(4,5) * qJD(3);
t286 = Ifges(4,4) * t197;
t353 = t272 / 0.2e1 + (t199 * Ifges(4,1) - t286) * t360 + t362;
t108 = -pkin(4) * t165 + qJD(5) + t154;
t75 = qJ(5) * t165 + t98;
t275 = t194 * t75;
t74 = -qJ(5) * t166 + t97;
t57 = pkin(4) * t188 + t74;
t20 = t195 * t57 - t275;
t11 = -pkin(5) * t188 + qJD(6) - t20;
t27 = pkin(5) * t104 - qJ(6) * t211 + t108;
t351 = t108 * mrSges(6,2) + t11 * mrSges(7,2) - t20 * mrSges(6,3) - t27 * mrSges(7,3) + t343 / 0.2e1;
t350 = -Ifges(6,6) / 0.2e1;
t349 = Ifges(7,6) / 0.2e1;
t323 = t68 / 0.2e1;
t322 = t69 / 0.2e1;
t315 = -t104 / 0.2e1;
t314 = t104 / 0.2e1;
t311 = t211 / 0.2e1;
t233 = -Ifges(4,6) * qJD(3) / 0.2e1;
t23 = t68 * mrSges(7,1) - t69 * mrSges(7,3);
t24 = t68 * mrSges(6,1) + t69 * mrSges(6,2);
t344 = -t23 - t24;
t237 = t197 * t254;
t239 = t198 * t256;
t241 = t196 * t256;
t253 = qJD(4) * t198;
t342 = (-t197 * t253 - t241) * t195 + (t237 - t239) * t194 + t161 * qJD(1);
t135 = t161 * t199;
t341 = -qJD(3) * t135 - t197 * t337 - t145;
t340 = qJ(2) * (m(3) + m(4));
t338 = Ifges(5,5) * t118 + Ifges(5,6) * t119;
t265 = t197 * t200;
t124 = t196 * t171 + t198 * t265;
t160 = qJD(3) * t225 + qJD(2);
t137 = t160 * qJD(1);
t42 = t196 * t137 + t150 * t253 - t153 * t254 + t182 * t239;
t205 = -qJD(4) * t98 + t198 * t137;
t43 = -t182 * t241 + t205;
t216 = -t196 * t43 + t198 * t42;
t336 = t361 * t229 - t345 * t68 + t346 * t69 + t338;
t121 = -mrSges(5,2) * t188 + mrSges(5,3) * t165;
t122 = mrSges(5,1) * t188 - mrSges(5,3) * t166;
t212 = -t196 * t121 - t198 * t122;
t335 = -m(5) * t215 + t212;
t14 = -qJ(5) * t118 - qJD(5) * t166 + (pkin(4) * qJD(1) - t182 * t196) * t256 + t205;
t22 = qJ(5) * t119 + qJD(5) * t165 + t42;
t4 = t194 * t14 + t195 * t22;
t1 = qJ(6) * t229 + qJD(6) * t188 + t4;
t3 = t14 * t195 - t194 * t22;
t2 = -pkin(5) * t229 - t3;
t334 = t43 * mrSges(5,1) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t42 * mrSges(5,2) - t4 * mrSges(6,2) + t1 * mrSges(7,3);
t70 = t195 * t75;
t21 = t194 * t57 + t70;
t12 = qJ(6) * t188 + t21;
t226 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t246 = t349 + t350;
t247 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t285 = Ifges(4,4) * t199;
t333 = -t246 * t104 - t226 * t188 - t247 * t211 - t12 * mrSges(7,3) - t20 * mrSges(6,1) - t97 * mrSges(5,1) - t233 + (-Ifges(4,2) * t197 + t285) * t360 - Ifges(6,6) * t315 - Ifges(7,6) * t314 - t166 * Ifges(5,5) - t165 * Ifges(5,6) + t11 * mrSges(7,1) + t21 * mrSges(6,2) + t98 * mrSges(5,2) - t346 * t311 - t361 * t301;
t36 = Ifges(7,5) * t211 + t188 * Ifges(7,6) + t104 * Ifges(7,3);
t39 = Ifges(6,4) * t211 - t104 * Ifges(6,2) + t188 * Ifges(6,6);
t330 = -t108 * mrSges(6,1) - t27 * mrSges(7,1) + t12 * mrSges(7,2) + t21 * mrSges(6,3) + t39 / 0.2e1 - t36 / 0.2e1;
t329 = m(5) / 0.2e1;
t328 = Ifges(7,5) * t322 + Ifges(7,3) * t323 + t229 * t349;
t327 = -t69 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t323 + t229 * t350;
t324 = -t68 / 0.2e1;
t312 = -t211 / 0.2e1;
t310 = t118 / 0.2e1;
t309 = t119 / 0.2e1;
t306 = -t165 / 0.2e1;
t304 = -t166 / 0.2e1;
t302 = -t188 / 0.2e1;
t297 = m(6) * t108;
t296 = pkin(4) * t166;
t204 = -qJD(4) * t124 + t198 * t160;
t228 = -t196 * t200 + pkin(4);
t34 = qJ(5) * t240 + (qJ(5) * t254 + qJD(3) * t228 - t251) * t199 + t204;
t255 = qJD(3) * t200;
t238 = t199 * t255;
t244 = t196 * t160 + t171 * t253 + t198 * t238;
t45 = -qJ(5) * t235 + (-qJD(5) * t199 + (qJ(5) * qJD(3) - qJD(4) * t200) * t197) * t196 + t244;
t9 = t194 * t34 + t195 * t45;
t52 = -mrSges(7,2) * t68 + mrSges(7,3) * t229;
t53 = -mrSges(6,2) * t229 - mrSges(6,3) * t68;
t291 = t52 + t53;
t54 = mrSges(6,1) * t229 - mrSges(6,3) * t69;
t55 = -mrSges(7,1) * t229 + t69 * mrSges(7,2);
t290 = t55 - t54;
t79 = -mrSges(6,2) * t188 - mrSges(6,3) * t104;
t80 = -mrSges(7,2) * t104 + mrSges(7,3) * t188;
t289 = t79 + t80;
t81 = mrSges(6,1) * t188 - mrSges(6,3) * t211;
t82 = -mrSges(7,1) * t188 + mrSges(7,2) * t211;
t288 = t81 - t82;
t287 = mrSges(4,1) * t197;
t279 = qJ(2) * mrSges(4,1);
t278 = qJ(2) * mrSges(4,2);
t159 = t198 * t171;
t100 = -qJ(5) * t264 + t197 * t228 + t159;
t109 = -qJ(5) * t266 + t124;
t48 = t194 * t100 + t195 * t109;
t270 = qJD(3) * mrSges(4,2);
t261 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t165 - mrSges(5,2) * t166 - mrSges(4,3) * t259;
t193 = -pkin(4) * t198 - pkin(3);
t167 = t182 * t258;
t234 = -t272 / 0.2e1;
t230 = t292 * t196;
t91 = -pkin(4) * t119 + t167;
t163 = pkin(4) * t266 - t199 * t200;
t224 = mrSges(5,1) * t198 - mrSges(5,2) * t196;
t221 = Ifges(5,1) * t196 + t282;
t219 = Ifges(5,2) * t198 + t283;
t217 = Ifges(5,5) * t196 + Ifges(5,6) * t198;
t8 = -t194 * t45 + t195 * t34;
t214 = t196 * t97 - t198 * t98;
t47 = t100 * t195 - t109 * t194;
t101 = mrSges(5,1) * t229 - mrSges(5,3) * t118;
t102 = -mrSges(5,2) * t229 + mrSges(5,3) * t119;
t213 = -t196 * t101 + t198 * t102;
t120 = pkin(4) * t352 + t197 * t255;
t132 = t162 * t197;
t192 = -pkin(4) * t195 - pkin(5);
t190 = pkin(4) * t194 + qJ(6);
t179 = t292 * t198;
t177 = -mrSges(4,3) * t260 - t270;
t168 = (t199 * mrSges(4,2) + t287) * qJD(1);
t133 = t162 * t199;
t123 = -t196 * t265 + t159;
t115 = -t195 * t179 + t194 * t230;
t114 = -t179 * t194 - t195 * t230;
t99 = pkin(5) * t161 - qJ(6) * t162 + t193;
t88 = qJD(3) * t339 - t199 * t337;
t86 = qJD(3) * t132 + t194 * t236 - t195 * t235;
t77 = -t196 * t238 + t204;
t76 = -t200 * t237 + t244;
t73 = -mrSges(5,1) * t119 + mrSges(5,2) * t118;
t72 = pkin(5) * t133 + qJ(6) * t135 + t163;
t59 = t118 * Ifges(5,1) + t119 * Ifges(5,4) + Ifges(5,5) * t229;
t58 = t118 * Ifges(5,4) + t119 * Ifges(5,2) + Ifges(5,6) * t229;
t50 = mrSges(6,1) * t104 + mrSges(6,2) * t211;
t49 = mrSges(7,1) * t104 - mrSges(7,3) * t211;
t46 = -pkin(5) * t197 - t47;
t44 = qJ(6) * t197 + t48;
t32 = pkin(5) * t211 + qJ(6) * t104 + t296;
t26 = t195 * t74 - t275;
t25 = t194 * t74 + t70;
t10 = -pkin(5) * t86 - qJ(6) * t88 + qJD(6) * t135 + t120;
t7 = -pkin(5) * t256 - t8;
t6 = qJ(6) * t256 + qJD(6) * t197 + t9;
t5 = pkin(5) * t68 - qJ(6) * t69 - qJD(6) * t211 + t91;
t13 = [(t168 + ((2 * mrSges(3,3)) + t287 + 0.2e1 * t340) * qJD(1)) * qJD(2) + m(5) * (t123 * t43 + t124 * t42 + t76 * t98 + t77 * t97) + (Ifges(6,2) * t315 - Ifges(7,3) * t314 + t301 * t345 + t311 * t347 + t330) * t86 + (t222 * t310 + t220 * t309 + t58 * t300 + t59 * t298 - t200 * t73 + qJD(1) * qJD(2) * mrSges(4,2) + (-t196 * t42 - t198 * t43) * mrSges(5,3) + (t154 * t224 + t219 * t306 + t221 * t304 + t217 * t302 + t95 * t300 - t198 * t94 / 0.2e1 + t214 * mrSges(5,3)) * qJD(4) + ((-m(5) * t200 + t223) * t170 + t200 * t177 + ((-0.3e1 / 0.2e1 * Ifges(4,4) + t281 / 0.2e1 - t280 / 0.2e1) * t199 + 0.2e1 * t279 - t247 * t135 + t246 * t133) * qJD(1) + t233 - t333) * qJD(3)) * t199 + (t336 + t338) * t197 / 0.2e1 + (-t133 * t4 + t135 * t3) * mrSges(6,3) + (-Ifges(6,4) * t135 - Ifges(6,2) * t133) * t324 + (-Ifges(7,5) * t135 + Ifges(7,3) * t133) * t323 + (-t1 * t133 - t135 * t2) * mrSges(7,2) + t5 * (mrSges(7,1) * t133 + mrSges(7,3) * t135) + t91 * (mrSges(6,1) * t133 - mrSges(6,2) * t135) + t133 * t327 + t133 * t328 + ((0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,1) + t226) * t229 + Ifges(6,6) * t324 + Ifges(7,6) * t323 + t346 * t322 + t334) * t197 + (-t133 * t347 - t135 * t348) * t322 + t10 * t49 + t44 * t52 + t48 * t53 + t47 * t54 + t46 * t55 + t72 * t23 + t9 * t79 + t6 * t80 + t8 * t81 + t7 * t82 + t120 * t50 + t76 * t121 + t77 * t122 + t123 * t101 + t124 * t102 + (Ifges(6,4) * t315 + Ifges(7,5) * t314 + t346 * t301 + t348 * t311 + t351) * t88 + (t234 + (-0.2e1 * t278 + 0.3e1 / 0.2e1 * t286) * qJD(1) + (m(5) * t154 - t261) * t200 - t353) * t258 + t163 * t24 - t359 * t135 / 0.2e1 + m(6) * (t108 * t120 + t163 * t91 + t20 * t8 + t21 * t9 + t3 * t47 + t4 * t48) + m(7) * (t1 * t44 + t10 * t27 + t11 * t7 + t12 * t6 + t2 * t46 + t5 * t72); -t291 * t339 + t290 * t132 + (-t73 + (t121 * t198 - t122 * t196 + t177) * qJD(3) + t344) * t199 + (t212 * qJD(4) + (t49 + t50 - t261) * qJD(3) + t213) * t197 - m(5) * t214 * t256 + 0.2e1 * ((-t215 * qJD(4) + t216) * t329 + (t297 / 0.2e1 + (t154 - t269) * t329) * qJD(3)) * t197 + t341 * t289 + t342 * t288 + (-t1 * t339 - t11 * t342 + t12 * t341 + t132 * t2 - t199 * t5 + t258 * t27) * m(7) + (-t132 * t3 - t199 * t91 + t20 * t342 + t21 * t341 - t339 * t4) * m(6) + (-t168 + (-mrSges(3,3) - t340) * qJD(1) + t335) * qJD(1); ((t50 + t297) * t196 * pkin(4) + t335 * pkin(8) + t362) * qJD(4) + (t36 - t39) * (t127 / 0.2e1 + t337 / 0.2e1) + (Ifges(6,4) * t128 - Ifges(7,5) * t146 + Ifges(6,2) * t127 + Ifges(7,3) * t337) * t314 + (-Ifges(6,4) * t146 + Ifges(7,5) * t128 - Ifges(6,2) * t337 - Ifges(7,3) * t127) * t315 + (-pkin(3) * t167 + pkin(8) * t216 - t110 * t97 - t111 * t98 - t154 * t170) * m(5) + t290 * t114 + t291 * t115 + ((-t177 - t270) * t199 + ((-mrSges(4,1) - t224) * qJD(3) + t261) * t197) * t182 + (mrSges(7,1) * t263 + mrSges(7,3) * t262) * t27 + (mrSges(6,1) * t263 - mrSges(6,2) * t262) * t108 + (-t161 * t4 - t162 * t3 + t20 * t262 - t21 * t263) * mrSges(6,3) + (-t1 * t161 - t11 * t262 - t12 * t263 + t162 * t2) * mrSges(7,2) + (Ifges(6,4) * t162 - Ifges(6,2) * t161) * t324 + t161 * t327 + t161 * t328 + t221 * t310 + (Ifges(7,5) * t162 + Ifges(7,3) * t161) * t323 + t58 * t298 + t219 * t309 - pkin(3) * t73 + t99 * t23 + t213 * pkin(8) + t216 * mrSges(5,3) - t111 * t121 - t110 * t122 - t130 * t50 + ((t234 + (t278 - t286 / 0.2e1) * qJD(1) + t353) * t197 + ((-t279 + t285 / 0.2e1 + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t197) * qJD(1) + t233 + t333) * t199 + (-t161 * t345 + t162 * t346 + t217) * t256 / 0.2e1) * qJD(1) + t5 * (mrSges(7,1) * t161 - mrSges(7,3) * t162) + t91 * (mrSges(6,1) * t161 + mrSges(6,2) * t162) + t354 * t49 + t343 * (-t128 / 0.2e1 - t146 / 0.2e1) + t355 * t79 + t356 * t81 + (-t108 * t130 - t114 * t3 + t115 * t4 + t193 * t91 + t20 * t356 + t21 * t355) * m(6) + t357 * t82 + t358 * t80 + (t1 * t115 + t11 * t357 + t114 * t2 + t12 * t358 + t27 * t354 + t5 * t99) * m(7) + t359 * t162 / 0.2e1 + t193 * t24 + (t127 * t345 + t128 * t346) * t302 + (-t146 * t346 - t337 * t345) * t301 + (t127 * t347 + t128 * t348) * t312 + (-t161 * t347 + t162 * t348) * t322 + (-t146 * t348 - t337 * t347) * t311 + t196 * t59 / 0.2e1; (-Ifges(6,2) * t314 + Ifges(7,3) * t315 - t312 * t347 + t330) * t211 + (-t108 * t296 + t20 * t25 - t21 * t26 + (t194 * t4 + t195 * t3) * pkin(4)) * m(6) + (t1 * t190 - t11 * t25 + t192 * t2 - t27 * t32 + (qJD(6) - t26) * t12) * m(7) + t336 + (t165 * t97 + t166 * t98) * mrSges(5,3) + t334 + (-t166 * t50 + t194 * t53 + t195 * t54) * pkin(4) + (Ifges(5,5) * t165 - Ifges(5,6) * t166 - t211 * t345) * t302 + t288 * t25 - t289 * t26 + (-Ifges(5,2) * t166 + t156 + t95) * t306 - t32 * t49 + t94 * t303 + (Ifges(5,1) * t165 - t284) * t304 + qJD(6) * t80 - t97 * t121 + t98 * t122 - t154 * (mrSges(5,1) * t166 + mrSges(5,2) * t165) + t190 * t52 + t192 * t55 + (-Ifges(6,4) * t314 - Ifges(7,5) * t315 - t346 * t302 - t348 * t312 + t351) * t104; t288 * t211 + t289 * t104 + (t104 * t12 - t11 * t211 + t5) * m(7) + (t104 * t21 + t20 * t211 + t91) * m(6) - t344; t211 * t49 - t188 * t80 + 0.2e1 * (t2 / 0.2e1 + t27 * t311 + t12 * t302) * m(7) + t55;];
tauc  = t13(:);
