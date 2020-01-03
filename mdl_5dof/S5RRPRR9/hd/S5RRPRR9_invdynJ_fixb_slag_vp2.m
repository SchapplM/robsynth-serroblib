% Calculate vector of inverse dynamics joint torques for
% S5RRPRR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:03
% EndTime: 2019-12-31 20:20:33
% DurationCPUTime: 16.25s
% Computational Cost: add. (8185->653), mult. (18678->872), div. (0->0), fcn. (13466->14), ass. (0->311)
t216 = sin(qJ(4));
t307 = sin(pkin(9));
t265 = t307 * pkin(2);
t202 = t265 + pkin(7);
t327 = pkin(8) + t202;
t258 = qJD(4) * t327;
t217 = sin(qJ(2));
t221 = cos(qJ(2));
t253 = qJD(1) * t307;
t308 = cos(pkin(9));
t254 = qJD(1) * t308;
t166 = -t217 * t253 + t221 * t254;
t304 = t166 * t216;
t214 = -qJ(3) - pkin(6);
t196 = t214 * t221;
t187 = qJD(1) * t196;
t170 = t307 * t187;
t194 = t214 * t217;
t186 = qJD(1) * t194;
t130 = t186 * t308 + t170;
t220 = cos(qJ(4));
t167 = -t217 * t254 - t221 * t253;
t285 = qJD(1) * t217;
t272 = pkin(2) * t285;
t99 = -pkin(3) * t167 - pkin(7) * t166 + t272;
t57 = t220 * t130 + t216 * t99;
t410 = -pkin(8) * t304 + t216 * t258 + t57;
t303 = t166 * t220;
t56 = -t130 * t216 + t220 * t99;
t409 = pkin(4) * t167 + pkin(8) * t303 - t220 * t258 - t56;
t331 = t220 * pkin(4);
t204 = pkin(3) + t331;
t213 = qJ(4) + qJ(5);
t210 = sin(t213);
t211 = cos(t213);
t247 = -mrSges(5,1) * t220 + mrSges(5,2) * t216;
t408 = -m(5) * pkin(3) - m(6) * t204 - t211 * mrSges(6,1) + t210 * mrSges(6,2) + t247;
t280 = qJD(4) * t216;
t407 = t280 - t304;
t386 = m(4) + m(5) + m(6);
t406 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t371 = m(6) * pkin(4);
t195 = -mrSges(3,1) * t221 + mrSges(3,2) * t217;
t212 = qJ(2) + pkin(9);
t208 = sin(t212);
t209 = cos(t212);
t405 = -t209 * mrSges(4,1) + t406 * t208 + t195;
t215 = sin(qJ(5));
t219 = cos(qJ(5));
t140 = qJD(2) * t220 + t167 * t216;
t141 = qJD(2) * t216 - t167 * t220;
t252 = t219 * t140 - t141 * t215;
t277 = qJD(1) * qJD(2);
t263 = t217 * t277;
t276 = qJDD(1) * t221;
t188 = -t263 + t276;
t189 = qJDD(1) * t217 + t221 * t277;
t132 = t188 * t307 + t189 * t308;
t74 = qJD(4) * t140 + qJDD(2) * t216 + t132 * t220;
t75 = -qJD(4) * t141 + qJDD(2) * t220 - t132 * t216;
t25 = qJD(5) * t252 + t215 * t75 + t219 * t74;
t370 = t25 / 0.2e1;
t81 = t140 * t215 + t141 * t219;
t26 = -qJD(5) * t81 - t215 * t74 + t219 * t75;
t369 = t26 / 0.2e1;
t362 = t74 / 0.2e1;
t361 = t75 / 0.2e1;
t131 = t188 * t308 - t307 * t189;
t128 = qJDD(4) - t131;
t125 = qJDD(5) + t128;
t356 = t125 / 0.2e1;
t355 = t128 / 0.2e1;
t404 = t188 / 0.2e1;
t403 = t216 * t371;
t329 = qJD(2) / 0.2e1;
t177 = t327 * t216;
t178 = t327 * t220;
t118 = -t177 * t219 - t178 * t215;
t402 = qJD(5) * t118 + t409 * t215 - t410 * t219;
t119 = -t177 * t215 + t178 * t219;
t401 = -qJD(5) * t119 + t410 * t215 + t409 * t219;
t157 = qJD(4) - t166;
t153 = qJD(5) + t157;
t398 = t157 * Ifges(5,3);
t399 = t140 * Ifges(5,6);
t400 = t141 * Ifges(5,5) + t81 * Ifges(6,5) + Ifges(6,6) * t252 + t153 * Ifges(6,3) + t398 + t399;
t156 = Ifges(4,4) * t166;
t397 = t166 * Ifges(4,2);
t330 = t221 * pkin(2);
t205 = pkin(1) + t330;
t190 = -qJD(1) * t205 + qJD(3);
t396 = t190 * mrSges(4,2);
t309 = qJDD(2) / 0.2e1;
t395 = mrSges(5,1) + t371;
t257 = t308 * t187;
t129 = t186 * t307 - t257;
t394 = t407 * pkin(4) - t129;
t315 = t167 * mrSges(4,3);
t393 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t140 + mrSges(5,2) * t141 - t315;
t185 = t215 * t220 + t216 * t219;
t385 = qJD(4) + qJD(5);
t135 = t385 * t185;
t94 = t185 * t166;
t392 = t94 - t135;
t236 = t215 * t216 - t219 * t220;
t134 = t385 * t236;
t95 = t236 * t166;
t391 = t95 - t134;
t390 = Ifges(4,5) * qJD(2);
t389 = Ifges(4,6) * qJD(2);
t182 = t217 * t308 + t221 * t307;
t106 = t236 * t182;
t181 = t217 * t307 - t221 * t308;
t121 = pkin(3) * t181 - pkin(7) * t182 - t205;
t138 = t194 * t307 - t196 * t308;
t133 = t220 * t138;
t71 = t216 * t121 + t133;
t169 = t181 * qJD(2);
t279 = qJD(4) * t220;
t233 = -t169 * t216 + t182 * t279;
t206 = pkin(6) * t276;
t179 = -pkin(6) * t263 + t206;
t180 = t189 * pkin(6);
t388 = t179 * t221 + t180 * t217;
t176 = qJD(2) * pkin(2) + t186;
t123 = t307 * t176 - t257;
t110 = qJD(2) * pkin(7) + t123;
t306 = qJDD(1) * pkin(1);
t155 = -pkin(2) * t188 + qJDD(3) - t306;
t54 = -pkin(3) * t131 - pkin(7) * t132 + t155;
t282 = qJD(3) * t217;
t115 = qJDD(2) * pkin(2) - qJ(3) * t189 - qJD(1) * t282 - t180;
t283 = qJD(2) * t217;
t269 = pkin(6) * t283;
t281 = qJD(3) * t221;
t126 = qJ(3) * t188 + t206 + (-t269 + t281) * qJD(1);
t66 = t307 * t115 + t308 * t126;
t60 = qJDD(2) * pkin(7) + t66;
t88 = -pkin(3) * t166 + pkin(7) * t167 + t190;
t11 = -t110 * t280 + t216 * t54 + t220 * t60 + t88 * t279;
t50 = t110 * t220 + t216 * t88;
t12 = -t50 * qJD(4) - t216 * t60 + t220 * t54;
t387 = t11 * t220 - t12 * t216;
t384 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t383 = 0.2e1 * t309;
t43 = pkin(8) * t140 + t50;
t310 = t219 * t43;
t49 = -t110 * t216 + t220 * t88;
t42 = -pkin(8) * t141 + t49;
t38 = pkin(4) * t157 + t42;
t14 = t215 * t38 + t310;
t122 = t176 * t308 + t170;
t109 = -qJD(2) * pkin(3) - t122;
t76 = -t140 * pkin(4) + t109;
t381 = -t76 * mrSges(6,1) + t14 * mrSges(6,3);
t311 = t215 * t43;
t13 = t219 * t38 - t311;
t380 = t76 * mrSges(6,2) - t13 * mrSges(6,3);
t10 = pkin(8) * t75 + t11;
t9 = pkin(4) * t128 - pkin(8) * t74 + t12;
t2 = qJD(5) * t13 + t10 * t219 + t215 * t9;
t3 = -qJD(5) * t14 - t10 * t215 + t219 * t9;
t378 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t377 = t12 * mrSges(5,1) - t11 * mrSges(5,2);
t376 = m(3) * pkin(1) + mrSges(2,1) - t405;
t375 = t190 * mrSges(4,1) + t49 * mrSges(5,1) + t13 * mrSges(6,1) - t50 * mrSges(5,2) - t14 * mrSges(6,2);
t373 = Ifges(6,4) * t370 + Ifges(6,2) * t369 + Ifges(6,6) * t356;
t372 = Ifges(6,1) * t370 + Ifges(6,4) * t369 + Ifges(6,5) * t356;
t368 = Ifges(5,1) * t362 + Ifges(5,4) * t361 + Ifges(5,5) * t355;
t342 = Ifges(6,4) * t81;
t35 = Ifges(6,2) * t252 + Ifges(6,6) * t153 + t342;
t367 = -t35 / 0.2e1;
t366 = t35 / 0.2e1;
t77 = Ifges(6,4) * t252;
t36 = Ifges(6,1) * t81 + Ifges(6,5) * t153 + t77;
t365 = -t36 / 0.2e1;
t364 = t36 / 0.2e1;
t360 = -t252 / 0.2e1;
t359 = t252 / 0.2e1;
t358 = -t81 / 0.2e1;
t357 = t81 / 0.2e1;
t354 = -t140 / 0.2e1;
t353 = -t141 / 0.2e1;
t352 = t141 / 0.2e1;
t351 = -t153 / 0.2e1;
t350 = t153 / 0.2e1;
t349 = -t157 / 0.2e1;
t348 = -t166 / 0.2e1;
t347 = -t167 / 0.2e1;
t343 = t220 / 0.2e1;
t340 = pkin(4) * t141;
t338 = pkin(6) * t221;
t337 = pkin(7) * t208;
t334 = g(3) * t208;
t326 = mrSges(5,3) * t140;
t325 = mrSges(5,3) * t141;
t324 = Ifges(3,4) * t217;
t323 = Ifges(3,4) * t221;
t322 = Ifges(5,4) * t141;
t321 = Ifges(5,4) * t216;
t320 = Ifges(5,4) * t220;
t316 = t166 * mrSges(4,3);
t314 = t167 * Ifges(4,4);
t301 = t169 * t220;
t298 = t182 * t216;
t297 = t182 * t220;
t223 = -pkin(8) - pkin(7);
t296 = t208 * t223;
t218 = sin(qJ(1));
t295 = t210 * t218;
t222 = cos(qJ(1));
t294 = t210 * t222;
t293 = t211 * t218;
t292 = t211 * t222;
t291 = t216 * t218;
t290 = t216 * t222;
t289 = t218 * t220;
t288 = t220 * t222;
t148 = t209 * t295 + t292;
t149 = -t209 * t293 + t294;
t287 = -t148 * mrSges(6,1) + t149 * mrSges(6,2);
t150 = -t209 * t294 + t293;
t151 = t209 * t292 + t295;
t286 = t150 * mrSges(6,1) - t151 * mrSges(6,2);
t284 = qJD(1) * t221;
t275 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t285) * t338;
t274 = Ifges(6,5) * t25 + Ifges(6,6) * t26 + Ifges(6,3) * t125;
t273 = Ifges(5,5) * t74 + Ifges(5,6) * t75 + Ifges(5,3) * t128;
t271 = pkin(2) * t283;
t139 = Ifges(5,4) * t140;
t69 = t141 * Ifges(5,1) + t157 * Ifges(5,5) + t139;
t267 = t69 * t343;
t266 = t308 * pkin(2);
t261 = -t280 / 0.2e1;
t168 = t182 * qJD(2);
t100 = pkin(3) * t168 + pkin(7) * t169 + t271;
t259 = qJD(2) * t214;
t164 = t217 * t259 + t281;
t165 = t221 * t259 - t282;
t98 = t164 * t308 + t165 * t307;
t260 = t220 * t100 - t216 * t98;
t255 = -t131 * mrSges(4,1) + t132 * mrSges(4,2);
t70 = t220 * t121 - t138 * t216;
t203 = -t266 - pkin(3);
t250 = pkin(3) * t209 + t337;
t249 = mrSges(3,1) * t217 + mrSges(3,2) * t221;
t246 = mrSges(5,1) * t216 + mrSges(5,2) * t220;
t245 = -mrSges(6,1) * t210 - mrSges(6,2) * t211;
t244 = Ifges(5,1) * t220 - t321;
t243 = t221 * Ifges(3,2) + t324;
t242 = -Ifges(5,2) * t216 + t320;
t241 = Ifges(3,5) * t221 - Ifges(3,6) * t217;
t240 = Ifges(5,5) * t220 - Ifges(5,6) * t216;
t47 = pkin(4) * t181 - pkin(8) * t297 + t70;
t53 = -pkin(8) * t298 + t71;
t27 = -t215 * t53 + t219 * t47;
t28 = t215 * t47 + t219 * t53;
t86 = -mrSges(5,2) * t157 + t326;
t87 = mrSges(5,1) * t157 - t325;
t238 = -t216 * t87 + t220 * t86;
t65 = t115 * t308 - t307 * t126;
t97 = t164 * t307 - t308 * t165;
t137 = -t308 * t194 - t196 * t307;
t237 = t209 * t204 - t296;
t235 = t274 + t378;
t234 = pkin(1) * t249;
t160 = -t209 * t290 + t289;
t158 = t209 * t291 + t288;
t232 = t182 * t280 + t301;
t231 = t109 * t246;
t230 = t217 * (Ifges(3,1) * t221 - t324);
t31 = t216 * t100 + t121 * t279 - t138 * t280 + t220 * t98;
t59 = -qJDD(2) * pkin(3) - t65;
t207 = Ifges(3,4) * t284;
t193 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t284;
t191 = t203 - t331;
t174 = Ifges(3,1) * t285 + Ifges(3,5) * qJD(2) + t207;
t173 = Ifges(3,6) * qJD(2) + qJD(1) * t243;
t161 = t209 * t288 + t291;
t159 = -t209 * t289 + t290;
t146 = -qJD(2) * mrSges(4,2) + t316;
t116 = -mrSges(4,1) * t166 - mrSges(4,2) * t167;
t112 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t132;
t111 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t131;
t105 = t185 * t182;
t104 = -t167 * Ifges(4,1) + t156 + t390;
t103 = -t314 + t389 + t397;
t96 = pkin(4) * t298 + t137;
t68 = t140 * Ifges(5,2) + t157 * Ifges(5,6) + t322;
t62 = mrSges(6,1) * t153 - mrSges(6,3) * t81;
t61 = -mrSges(6,2) * t153 + mrSges(6,3) * t252;
t58 = pkin(4) * t233 + t97;
t46 = -mrSges(5,2) * t128 + mrSges(5,3) * t75;
t45 = mrSges(5,1) * t128 - mrSges(5,3) * t74;
t44 = -mrSges(6,1) * t252 + mrSges(6,2) * t81;
t40 = t106 * t385 + t185 * t169;
t39 = -t135 * t182 + t169 * t236;
t37 = -mrSges(5,1) * t75 + mrSges(5,2) * t74;
t33 = -t75 * pkin(4) + t59;
t32 = -qJD(4) * t71 + t260;
t29 = t74 * Ifges(5,4) + t75 * Ifges(5,2) + t128 * Ifges(5,6);
t24 = -pkin(8) * t233 + t31;
t21 = pkin(8) * t301 + pkin(4) * t168 + (-t133 + (pkin(8) * t182 - t121) * t216) * qJD(4) + t260;
t18 = -mrSges(6,2) * t125 + mrSges(6,3) * t26;
t17 = mrSges(6,1) * t125 - mrSges(6,3) * t25;
t16 = t219 * t42 - t311;
t15 = -t215 * t42 - t310;
t8 = -mrSges(6,1) * t26 + mrSges(6,2) * t25;
t5 = -qJD(5) * t28 + t21 * t219 - t215 * t24;
t4 = qJD(5) * t27 + t21 * t215 + t219 * t24;
t1 = [(t155 * mrSges(4,2) - t65 * mrSges(4,3) + Ifges(4,1) * t132 + Ifges(4,4) * t131 + Ifges(4,5) * t383 + t240 * t355 + t242 * t361 + t244 * t362 + t59 * t246 + t69 * t261) * t182 + (t155 * mrSges(4,1) - t66 * mrSges(4,3) - Ifges(4,4) * t132 + Ifges(5,5) * t362 + Ifges(6,5) * t370 - Ifges(4,2) * t131 - t383 * Ifges(4,6) + Ifges(5,6) * t361 + Ifges(6,6) * t369 + Ifges(5,3) * t355 + Ifges(6,3) * t356 + t377 + t378) * t181 + (t188 * t338 + t388) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t388) - t233 * t68 / 0.2e1 + (-m(4) * t122 + m(5) * t109 + t393) * t97 + t189 * t323 / 0.2e1 + t109 * (mrSges(5,1) * t233 - mrSges(5,2) * t232) + (Ifges(6,4) * t39 + Ifges(6,2) * t40) * t359 + (-t267 - t104 / 0.2e1 - t156 / 0.2e1 - Ifges(4,5) * t329 - Ifges(4,1) * t347 - t396 + t122 * mrSges(4,3)) * t169 + m(6) * (t13 * t5 + t14 * t4 + t2 * t28 + t27 * t3 + t33 * t96 + t58 * t76) + (-m(4) * t65 + m(5) * t59 - t112 + t37) * t137 + t140 * (-Ifges(5,4) * t232 - Ifges(5,2) * t233) / 0.2e1 + (-t159 * mrSges(5,1) - t149 * mrSges(6,1) - t158 * mrSges(5,2) - t148 * mrSges(6,2) + (t214 * t386 + t384 - t403) * t222 + (-m(6) * (-t205 - t237) - m(5) * (-t205 - t250) + m(4) * t205 + t376) * t218) * g(1) + (Ifges(3,1) * t189 + Ifges(3,4) * t404 - pkin(6) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t189) + t383 * Ifges(3,5)) * t217 + (-t291 * t371 - t161 * mrSges(5,1) - t151 * mrSges(6,1) - t160 * mrSges(5,2) - t150 * mrSges(6,2) - t386 * (t222 * t205 - t218 * t214) + t384 * t218 + (-m(5) * t250 - m(6) * t237 - t376) * t222) * g(2) + (t241 * t329 - t275) * qJD(2) + (-t11 * t298 - t12 * t297 + t232 * t49 - t233 * t50) * mrSges(5,3) + m(5) * (t11 * t71 + t12 * t70 + t31 * t50 + t32 * t49) + (Ifges(6,5) * t39 + Ifges(6,6) * t40) * t350 + t243 * t404 + t157 * (-Ifges(5,5) * t232 - Ifges(5,6) * t233) / 0.2e1 + (t400 / 0.2e1 + t399 / 0.2e1 + t398 / 0.2e1 - t103 / 0.2e1 - t397 / 0.2e1 - Ifges(4,6) * t329 - Ifges(4,4) * t347 + Ifges(6,3) * t350 + Ifges(5,5) * t352 + Ifges(6,5) * t357 + Ifges(6,6) * t359 + t375 - t123 * mrSges(4,3)) * t168 + m(4) * (t123 * t98 + t138 * t66 - t155 * t205 + t190 * t271) + t116 * t271 + t96 * t8 + (-Ifges(6,4) * t106 - Ifges(6,2) * t105) * t369 + (-Ifges(6,5) * t106 - Ifges(6,6) * t105) * t356 + (-t105 * t2 + t106 * t3 - t13 * t39 + t14 * t40) * mrSges(6,3) + t33 * (mrSges(6,1) * t105 - mrSges(6,2) * t106) + (-Ifges(6,1) * t106 - Ifges(6,4) * t105) * t370 + t31 * t86 + t32 * t87 + (t221 * (-Ifges(3,2) * t217 + t323) + t230) * t277 / 0.2e1 + (t274 + t273) * t181 / 0.2e1 + t76 * (-mrSges(6,1) * t40 + mrSges(6,2) * t39) + t70 * t45 + t71 * t46 + t58 * t44 + t4 * t61 + t5 * t62 + t27 * t17 + t28 * t18 - qJDD(2) * mrSges(3,2) * t338 + Ifges(2,3) * qJDD(1) - t195 * t306 - t29 * t298 / 0.2e1 + t138 * t111 + t98 * t146 + Ifges(3,6) * t221 * t309 + t221 * t174 * t329 + t39 * t364 + t40 * t366 + t297 * t368 - t106 * t372 - t105 * t373 - pkin(1) * (-mrSges(3,1) * t188 + mrSges(3,2) * t189) + t221 * (Ifges(3,4) * t189 + Ifges(3,2) * t188 + Ifges(3,6) * qJDD(2)) / 0.2e1 + (-Ifges(5,1) * t232 - Ifges(5,4) * t233) * t352 - t205 * t255 + (Ifges(6,1) * t39 + Ifges(6,4) * t40) * t357 - t193 * t269 - t234 * t277 - t173 * t283 / 0.2e1; ((t307 * t66 + t308 * t65) * pkin(2) + t122 * t129 - t123 * t130 - t190 * t272) * m(4) + (-t109 * t129 + t203 * t59 - t49 * t56 - t50 * t57) * m(5) - (Ifges(6,1) * t357 + Ifges(6,4) * t359 + Ifges(6,5) * t350 + t364 + t380) * t134 - (Ifges(6,4) * t357 + Ifges(6,2) * t359 + Ifges(6,6) * t350 + t366 + t381) * t135 + (m(5) * ((-t50 * t216 - t49 * t220) * qJD(4) + t387) - t87 * t279 - t86 * t280 - t216 * t45 + t220 * t46) * t202 + (-t389 / 0.2e1 + Ifges(4,2) * t348 - Ifges(5,3) * t349 - Ifges(6,3) * t351 - Ifges(5,5) * t353 - Ifges(5,6) * t354 - Ifges(6,5) * t358 - Ifges(6,6) * t360 + t375) * t167 - t393 * t129 + t394 * t44 + (-t231 - t390 / 0.2e1 + t240 * t349 + t244 * t353 + t242 * t354 - t396) * t166 + (-t407 * t50 + (-t279 + t303) * t49 + t387) * mrSges(5,3) + (g(1) * t222 + g(2) * t218) * (t249 + t386 * pkin(2) * t217 + (-m(5) * pkin(7) + m(6) * t223 + t406) * t209 + (mrSges(4,1) - t408) * t208) + (-m(5) * (t330 + t337) - m(4) * t330 - m(6) * (-t296 + t330) + t408 * t209 + t405) * g(3) + (pkin(6) * t193 + t173 / 0.2e1) * t285 + (Ifges(4,1) * t166 + t314 + t400) * t167 / 0.2e1 + t401 * t62 + t402 * t61 + (t118 * t3 + t119 * t2 + t13 * t401 + t14 * t402 + t191 * t33 + t394 * t76) * m(6) + (t156 + t104) * t348 + (t267 + t231) * qJD(4) + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (t304 / 0.2e1 + t261) * t68 - t57 * t86 - t56 * t87 - (-Ifges(3,2) * t285 + t174 + t207) * t284 / 0.2e1 + (t140 * t242 + t141 * t244 + t157 * t240) * qJD(4) / 0.2e1 + (-Ifges(6,4) * t95 - Ifges(6,2) * t94) * t360 + (-Ifges(6,1) * t95 - Ifges(6,4) * t94) * t358 + (-Ifges(6,5) * t95 - Ifges(6,6) * t94) * t351 - t76 * (mrSges(6,1) * t94 - mrSges(6,2) * t95) + t65 * mrSges(4,1) - t66 * mrSges(4,2) + (t275 + (-t230 / 0.2e1 + t234) * qJD(1)) * qJD(1) - t123 * t315 - t69 * t303 / 0.2e1 + t118 * t17 + t119 * t18 + Ifges(4,6) * t131 + Ifges(4,5) * t132 + t59 * t247 - t130 * t146 - t179 * mrSges(3,2) - t180 * mrSges(3,1) + t111 * t265 + t112 * t266 + t122 * t316 + t29 * t343 + t103 * t347 + (Ifges(5,5) * t216 + Ifges(5,6) * t220) * t355 + (Ifges(5,2) * t220 + t321) * t361 + (Ifges(5,1) * t216 + t320) * t362 - t95 * t365 - t94 * t367 + t216 * t368 + t185 * t372 + Ifges(3,6) * t188 + Ifges(3,5) * t189 + t191 * t8 + t203 * t37 + (-t13 * t95 + t14 * t94 - t185 * t3 - t2 * t236) * mrSges(6,3) + t33 * (mrSges(6,1) * t236 + mrSges(6,2) * t185) + (Ifges(6,5) * t185 - Ifges(6,6) * t236) * t356 + (Ifges(6,4) * t185 - Ifges(6,2) * t236) * t369 + (Ifges(6,1) * t185 - Ifges(6,4) * t236) * t370 - t236 * t373 - t116 * t272 - t241 * t277 / 0.2e1; -t236 * t17 + t185 * t18 + t216 * t46 + t220 * t45 + t392 * t62 + t391 * t61 + t238 * qJD(4) + (t44 + t393) * t167 + (-t146 - t238) * t166 + t255 + (-g(1) * t218 + g(2) * t222) * t386 + (t13 * t392 + t14 * t391 + t167 * t76 + t185 * t2 - t236 * t3) * m(6) + (t109 * t167 + t11 * t216 + t12 * t220 + t157 * (-t216 * t49 + t220 * t50)) * m(5) + (-t122 * t167 - t123 * t166 + t155) * m(4); (mrSges(5,2) * t161 - t160 * t395 - t286) * g(1) + (-mrSges(5,2) * t159 + t158 * t395 - t287) * g(2) + (-Ifges(5,2) * t141 + t139 + t69) * t354 + (-t86 + t326) * t49 + (t87 + t325) * t50 + t235 + (-t245 + t246 + t403) * t334 + t377 - (Ifges(6,4) * t358 + Ifges(6,2) * t360 + Ifges(6,6) * t351 + t367 - t381) * t81 + (Ifges(6,1) * t358 + Ifges(6,4) * t360 + Ifges(6,5) * t351 + t365 - t380) * t252 - t16 * t61 - t15 * t62 + t273 - t44 * t340 - m(6) * (t13 * t15 + t14 * t16 + t340 * t76) - t109 * (mrSges(5,1) * t141 + mrSges(5,2) * t140) + (Ifges(5,5) * t140 - Ifges(5,6) * t141) * t349 + t68 * t352 + (Ifges(5,1) * t140 - t322) * t353 + (t2 * t215 + t219 * t3 + (-t13 * t215 + t14 * t219) * qJD(5)) * t371 + ((-t215 * t62 + t219 * t61) * qJD(5) + t17 * t219 + t18 * t215) * pkin(4); -t76 * (mrSges(6,1) * t81 + mrSges(6,2) * t252) + (Ifges(6,1) * t252 - t342) * t358 + t35 * t357 + (Ifges(6,5) * t252 - Ifges(6,6) * t81) * t351 - t13 * t61 + t14 * t62 - g(1) * t286 - g(2) * t287 - t245 * t334 + (t13 * t252 + t14 * t81) * mrSges(6,3) + t235 + (-Ifges(6,2) * t81 + t36 + t77) * t360;];
tau = t1;
