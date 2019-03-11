% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:23
% EndTime: 2019-03-09 07:00:51
% DurationCPUTime: 13.30s
% Computational Cost: add. (13559->647), mult. (31697->906), div. (0->0), fcn. (21529->10), ass. (0->289)
t249 = sin(pkin(11)) * pkin(1) + pkin(7);
t233 = t249 * qJD(1);
t259 = sin(qJ(3));
t263 = cos(qJ(3));
t191 = qJD(2) * t263 - t259 * t233;
t286 = pkin(3) * t259 - pkin(8) * t263;
t227 = t286 * qJD(1);
t258 = sin(qJ(4));
t262 = cos(qJ(4));
t146 = -t258 * t191 + t262 * t227;
t317 = t262 * t263;
t272 = pkin(4) * t259 - pkin(9) * t317;
t380 = -pkin(9) - pkin(8);
t296 = qJD(4) * t380;
t439 = -qJD(1) * t272 + t262 * t296 - t146;
t147 = t262 * t191 + t258 * t227;
t311 = qJD(1) * t263;
t295 = t258 * t311;
t438 = -pkin(9) * t295 - t258 * t296 + t147;
t308 = qJD(3) * t262;
t312 = qJD(1) * t259;
t219 = -t258 * t312 + t308;
t220 = qJD(3) * t258 + t262 * t312;
t257 = sin(qJ(5));
t261 = cos(qJ(5));
t160 = t219 * t257 + t220 * t261;
t256 = sin(qJ(6));
t260 = cos(qJ(6));
t287 = t261 * t219 - t220 * t257;
t101 = t160 * t260 + t256 * t287;
t182 = -qJD(3) * pkin(3) - t191;
t153 = -t219 * pkin(4) + t182;
t103 = -pkin(5) * t287 + t153;
t404 = pkin(10) * t287;
t192 = t259 * qJD(2) + t263 * t233;
t183 = qJD(3) * pkin(8) + t192;
t297 = -cos(pkin(11)) * pkin(1) - pkin(2);
t213 = -pkin(3) * t263 - t259 * pkin(8) + t297;
t186 = t213 * qJD(1);
t122 = -t183 * t258 + t262 * t186;
t108 = -pkin(9) * t220 + t122;
t247 = qJD(4) - t311;
t102 = pkin(4) * t247 + t108;
t123 = t183 * t262 + t186 * t258;
t109 = pkin(9) * t219 + t123;
t107 = t261 * t109;
t51 = t102 * t257 + t107;
t41 = t51 + t404;
t331 = t256 * t41;
t240 = qJD(5) + t247;
t419 = pkin(10) * t160;
t105 = t257 * t109;
t50 = t261 * t102 - t105;
t40 = t50 - t419;
t38 = pkin(5) * t240 + t40;
t14 = t260 * t38 - t331;
t329 = t260 * t41;
t15 = t256 * t38 + t329;
t309 = qJD(3) * t259;
t289 = qJD(1) * t309;
t244 = Ifges(7,3) * t289;
t345 = Ifges(7,4) * t101;
t232 = qJD(6) + t240;
t363 = -t232 / 0.2e1;
t377 = t101 / 0.2e1;
t378 = -t101 / 0.2e1;
t408 = -t160 * t256 + t260 * t287;
t382 = -t408 / 0.2e1;
t303 = qJD(5) * t261;
t304 = qJD(5) * t257;
t300 = qJD(3) * qJD(4);
t306 = qJD(4) * t258;
t307 = qJD(3) * t263;
t176 = t262 * t300 + (-t259 * t306 + t262 * t307) * qJD(1);
t184 = t191 * qJD(3);
t230 = t286 * qJD(3);
t212 = qJD(1) * t230;
t71 = -qJD(4) * t123 - t184 * t258 + t262 * t212;
t58 = pkin(4) * t289 - pkin(9) * t176 + t71;
t305 = qJD(4) * t262;
t411 = t258 * t307 + t259 * t305;
t177 = -qJD(1) * t411 - t258 * t300;
t70 = -t183 * t306 + t262 * t184 + t186 * t305 + t258 * t212;
t61 = pkin(9) * t177 + t70;
t12 = t102 * t303 - t109 * t304 + t257 * t58 + t261 * t61;
t75 = -qJD(5) * t160 - t176 * t257 + t177 * t261;
t10 = pkin(10) * t75 + t12;
t13 = -qJD(5) * t51 - t257 * t61 + t261 * t58;
t74 = qJD(5) * t287 + t176 * t261 + t177 * t257;
t9 = pkin(5) * t289 - pkin(10) * t74 + t13;
t2 = qJD(6) * t14 + t10 * t260 + t256 * t9;
t3 = -qJD(6) * t15 - t10 * t256 + t260 * t9;
t30 = qJD(6) * t408 + t256 * t75 + t260 * t74;
t31 = -qJD(6) * t101 - t256 * t74 + t260 * t75;
t391 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t30 + Ifges(7,6) * t31;
t47 = Ifges(7,2) * t408 + Ifges(7,6) * t232 + t345;
t95 = Ifges(7,4) * t408;
t48 = Ifges(7,1) * t101 + Ifges(7,5) * t232 + t95;
t437 = (Ifges(7,1) * t408 - t345) * t378 + (Ifges(7,5) * t408 - Ifges(7,6) * t101) * t363 + (t101 * t15 + t14 * t408) * mrSges(7,3) - t103 * (mrSges(7,1) * t101 + mrSges(7,2) * t408) + t47 * t377 + t244 + t391 + (-Ifges(7,2) * t101 + t48 + t95) * t382;
t273 = t257 * t258 - t261 * t262;
t393 = qJD(4) + qJD(5);
t165 = t393 * t273;
t269 = t273 * t263;
t188 = qJD(1) * t269;
t436 = t165 - t188;
t222 = t257 * t262 + t258 * t261;
t166 = t393 * t222;
t270 = t222 * t263;
t187 = qJD(1) * t270;
t429 = t166 - t187;
t154 = Ifges(6,4) * t287;
t245 = Ifges(6,3) * t289;
t346 = Ifges(6,4) * t160;
t361 = -t240 / 0.2e1;
t372 = -t160 / 0.2e1;
t374 = -t287 / 0.2e1;
t390 = t13 * mrSges(6,1) - t12 * mrSges(6,2) + Ifges(6,5) * t74 + Ifges(6,6) * t75;
t91 = Ifges(6,1) * t160 + Ifges(6,5) * t240 + t154;
t435 = t245 + t390 + (Ifges(6,5) * t287 - Ifges(6,6) * t160) * t361 + (t160 * t51 + t287 * t50) * mrSges(6,3) + (-Ifges(6,2) * t160 + t154 + t91) * t374 - t153 * (mrSges(6,1) * t160 + mrSges(6,2) * t287) + (Ifges(6,1) * t287 - t346) * t372 + t437;
t238 = t380 * t258;
t239 = t380 * t262;
t173 = t257 * t238 - t261 * t239;
t414 = -qJD(5) * t173 + t438 * t257 + t261 * t439;
t413 = t238 * t303 + t239 * t304 + t257 * t439 - t438 * t261;
t432 = -pkin(5) * t312 + pkin(10) * t436 + t414;
t431 = pkin(10) * t429 - t413;
t421 = -Ifges(4,1) / 0.2e1;
t252 = Ifges(4,4) * t311;
t420 = -t252 / 0.2e1;
t172 = t261 * t238 + t239 * t257;
t140 = -pkin(10) * t222 + t172;
t141 = -pkin(10) * t273 + t173;
t83 = t140 * t260 - t141 * t256;
t418 = qJD(6) * t83 + t256 * t432 - t260 * t431;
t84 = t140 * t256 + t141 * t260;
t417 = -qJD(6) * t84 + t256 * t431 + t260 * t432;
t250 = pkin(4) * t261 + pkin(5);
t301 = qJD(6) * t260;
t302 = qJD(6) * t256;
t319 = t257 * t260;
t53 = -t108 * t257 - t107;
t44 = t53 - t404;
t54 = t261 * t108 - t105;
t45 = t54 - t419;
t416 = t256 * t45 - t260 * t44 - t250 * t302 + (-t257 * t301 + (-t256 * t261 - t319) * qJD(5)) * pkin(4);
t320 = t256 * t257;
t415 = -t256 * t44 - t260 * t45 + t250 * t301 + (-t257 * t302 + (t260 * t261 - t320) * qJD(5)) * pkin(4);
t169 = pkin(4) * t295 + t192;
t412 = pkin(4) * t306 + pkin(5) * t429 - t169;
t337 = t220 * Ifges(5,4);
t143 = t219 * Ifges(5,2) + t247 * Ifges(5,6) + t337;
t214 = Ifges(5,4) * t219;
t144 = t220 * Ifges(5,1) + t247 * Ifges(5,5) + t214;
t277 = t122 * t262 + t123 * t258;
t347 = Ifges(5,4) * t262;
t281 = -Ifges(5,2) * t258 + t347;
t348 = Ifges(5,4) * t258;
t283 = Ifges(5,1) * t262 - t348;
t284 = mrSges(5,1) * t258 + mrSges(5,2) * t262;
t343 = Ifges(5,6) * t258;
t344 = Ifges(5,5) * t262;
t357 = t262 / 0.2e1;
t358 = -t258 / 0.2e1;
t364 = t220 / 0.2e1;
t265 = -t277 * mrSges(5,3) + t182 * t284 + t219 * t281 / 0.2e1 + t283 * t364 + t247 * (-t343 + t344) / 0.2e1 + t143 * t358 + t144 * t357;
t325 = Ifges(4,5) * qJD(3);
t410 = t191 * mrSges(4,3) + t312 * t421 + t420 - t325 / 0.2e1 - t265;
t90 = Ifges(6,2) * t287 + Ifges(6,6) * t240 + t346;
t405 = t90 / 0.2e1;
t291 = -Ifges(4,6) * qJD(3) / 0.2e1;
t196 = t273 * t259;
t104 = -mrSges(6,1) * t287 + mrSges(6,2) * t160;
t394 = m(6) * t153 + t104;
t198 = t262 * t213;
t321 = t249 * t258;
t355 = pkin(9) * t259;
t138 = -t262 * t355 + t198 + (-pkin(4) - t321) * t263;
t226 = t249 * t317;
t164 = t258 * t213 + t226;
t318 = t258 * t259;
t148 = -pkin(9) * t318 + t164;
t86 = t257 * t138 + t261 * t148;
t278 = -t258 * t71 + t262 * t70;
t389 = -t71 * mrSges(5,1) + t70 * mrSges(5,2) - Ifges(5,5) * t176 - Ifges(5,6) * t177;
t388 = t30 / 0.2e1;
t387 = t31 / 0.2e1;
t384 = t74 / 0.2e1;
t383 = t75 / 0.2e1;
t381 = t408 / 0.2e1;
t195 = t222 * t259;
t136 = -t195 * t260 + t196 * t256;
t376 = t136 / 0.2e1;
t137 = -t195 * t256 - t196 * t260;
t375 = t137 / 0.2e1;
t373 = t287 / 0.2e1;
t371 = t160 / 0.2e1;
t370 = t176 / 0.2e1;
t369 = t177 / 0.2e1;
t368 = -t195 / 0.2e1;
t367 = -t196 / 0.2e1;
t366 = -t219 / 0.2e1;
t365 = -t220 / 0.2e1;
t362 = t232 / 0.2e1;
t360 = t240 / 0.2e1;
t359 = -t247 / 0.2e1;
t349 = Ifges(4,4) * t259;
t235 = t297 * qJD(1);
t334 = t235 * mrSges(4,2);
t127 = -t187 * t260 + t188 * t256;
t162 = t222 * t260 - t256 * t273;
t67 = -qJD(6) * t162 + t165 * t256 - t166 * t260;
t327 = t127 - t67;
t128 = -t187 * t256 - t188 * t260;
t161 = -t222 * t256 - t260 * t273;
t66 = qJD(6) * t161 - t165 * t260 - t166 * t256;
t326 = t128 - t66;
t298 = mrSges(4,3) * t312;
t314 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t219 + mrSges(5,2) * t220 + t298;
t313 = t262 * t230 + t309 * t321;
t199 = pkin(4) * t318 + t259 * t249;
t299 = mrSges(4,3) * t311;
t168 = t411 * pkin(4) + t249 * t307;
t251 = -pkin(4) * t262 - pkin(3);
t292 = t325 / 0.2e1;
t290 = m(4) * t249 + mrSges(4,3);
t85 = t261 * t138 - t257 * t148;
t285 = mrSges(5,1) * t262 - mrSges(5,2) * t258;
t282 = Ifges(5,1) * t258 + t347;
t280 = Ifges(5,2) * t262 + t348;
t279 = Ifges(5,5) * t258 + Ifges(5,6) * t262;
t62 = -pkin(5) * t263 + t196 * pkin(10) + t85;
t63 = -pkin(10) * t195 + t86;
t34 = -t256 * t63 + t260 * t62;
t35 = t256 * t62 + t260 * t63;
t276 = t122 * t258 - t123 * t262;
t155 = mrSges(5,1) * t289 - mrSges(5,3) * t176;
t156 = -mrSges(5,2) * t289 + mrSges(5,3) * t177;
t275 = -t258 * t155 + t262 * t156;
t179 = -mrSges(5,2) * t247 + mrSges(5,3) * t219;
t180 = mrSges(5,1) * t247 - mrSges(5,3) * t220;
t274 = -t258 * t179 - t262 * t180;
t87 = t272 * qJD(3) + (-t226 + (-t213 + t355) * t258) * qJD(4) + t313;
t110 = t213 * t305 + t258 * t230 + (-t259 * t308 - t263 * t306) * t249;
t92 = -pkin(9) * t411 + t110;
t32 = t138 * t303 - t148 * t304 + t257 * t87 + t261 * t92;
t185 = t192 * qJD(3);
t129 = -t177 * pkin(4) + t185;
t33 = -qJD(5) * t86 - t257 * t92 + t261 * t87;
t264 = t122 * mrSges(5,1) + t14 * mrSges(7,1) + t235 * mrSges(4,1) + t50 * mrSges(6,1) + t247 * Ifges(5,3) + t220 * Ifges(5,5) + t219 * Ifges(5,6) + t291 - (t263 * Ifges(4,2) + t349) * qJD(1) / 0.2e1 + t232 * Ifges(7,3) + t101 * Ifges(7,5) + t408 * Ifges(7,6) + t240 * Ifges(6,3) + t160 * Ifges(6,5) + t287 * Ifges(6,6) - t123 * mrSges(5,2) - t15 * mrSges(7,2) - t51 * mrSges(6,2);
t246 = Ifges(5,3) * t289;
t236 = -qJD(3) * mrSges(4,2) + t299;
t201 = pkin(4) * t319 + t250 * t256;
t200 = -pkin(4) * t320 + t250 * t260;
t189 = pkin(5) * t273 + t251;
t163 = -t263 * t321 + t198;
t149 = pkin(5) * t195 + t199;
t133 = mrSges(6,1) * t240 - mrSges(6,3) * t160;
t132 = -mrSges(6,2) * t240 + mrSges(6,3) * t287;
t124 = pkin(4) * t220 + pkin(5) * t160;
t120 = -mrSges(5,1) * t177 + mrSges(5,2) * t176;
t115 = t176 * Ifges(5,1) + t177 * Ifges(5,4) + Ifges(5,5) * t289;
t114 = t176 * Ifges(5,4) + t177 * Ifges(5,2) + Ifges(5,6) * t289;
t113 = -qJD(3) * t270 + t196 * t393;
t112 = -qJD(3) * t269 - t166 * t259;
t111 = -qJD(4) * t164 + t313;
t82 = mrSges(7,1) * t232 - mrSges(7,3) * t101;
t81 = -mrSges(7,2) * t232 + mrSges(7,3) * t408;
t78 = -pkin(5) * t113 + t168;
t69 = -mrSges(6,2) * t289 + mrSges(6,3) * t75;
t68 = mrSges(6,1) * t289 - mrSges(6,3) * t74;
t52 = -t75 * pkin(5) + t129;
t49 = -mrSges(7,1) * t408 + mrSges(7,2) * t101;
t43 = -qJD(6) * t137 - t112 * t256 + t113 * t260;
t42 = qJD(6) * t136 + t112 * t260 + t113 * t256;
t39 = -mrSges(6,1) * t75 + mrSges(6,2) * t74;
t37 = t74 * Ifges(6,1) + t75 * Ifges(6,4) + Ifges(6,5) * t289;
t36 = t74 * Ifges(6,4) + t75 * Ifges(6,2) + Ifges(6,6) * t289;
t25 = -mrSges(7,2) * t289 + mrSges(7,3) * t31;
t24 = mrSges(7,1) * t289 - mrSges(7,3) * t30;
t21 = pkin(10) * t113 + t32;
t20 = pkin(5) * t309 - pkin(10) * t112 + t33;
t17 = t260 * t40 - t331;
t16 = -t256 * t40 - t329;
t8 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t7 = t30 * Ifges(7,1) + t31 * Ifges(7,4) + Ifges(7,5) * t289;
t6 = t30 * Ifges(7,4) + t31 * Ifges(7,2) + Ifges(7,6) * t289;
t5 = -qJD(6) * t35 + t20 * t260 - t21 * t256;
t4 = qJD(6) * t34 + t20 * t256 + t21 * t260;
t1 = [(-t244 / 0.2e1 - t245 / 0.2e1 - t246 / 0.2e1 + t290 * t184 + (0.3e1 / 0.2e1 * t252 + t292 + 0.2e1 * t334 + (-m(4) * t191 + m(5) * t182 + t314) * t249 - t410) * qJD(3) + t389 - t390 - t391) * t263 + t129 * (mrSges(6,1) * t195 - mrSges(6,2) * t196) + t37 * t367 + t36 * t368 + (Ifges(6,1) * t112 + Ifges(6,4) * t113) * t371 + (Ifges(6,4) * t112 + Ifges(6,2) * t113) * t373 + t7 * t375 + t6 * t376 + (Ifges(6,5) * t112 + Ifges(6,6) * t113) * t360 + (Ifges(7,5) * t42 + Ifges(7,6) * t43) * t362 + t113 * t405 + (t283 * t370 + t281 * t369 + t114 * t358 + t115 * t357 + (-t258 * t70 - t262 * t71) * mrSges(5,3) + (mrSges(4,3) + t284) * t185 + (t144 * t358 - t262 * t143 / 0.2e1 + t182 * t285 + t280 * t366 + t282 * t365 + t279 * t359 + t276 * mrSges(5,3)) * qJD(4) + (t264 - t290 * t192 + t291 + (t297 * mrSges(4,1) + Ifges(7,5) * t375 + Ifges(7,6) * t376 + Ifges(6,5) * t367 + Ifges(6,6) * t368 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t344 / 0.2e1 - t343 / 0.2e1) * t259 + (-Ifges(5,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1) * t263) * qJD(1)) * qJD(3) + (t120 + (m(4) + m(5)) * t185 - t236 * qJD(3)) * t249) * t259 + (t136 * t2 - t137 * t3 - t14 * t42 + t15 * t43) * mrSges(7,3) + m(7) * (t103 * t78 + t14 * t5 + t149 * t52 + t15 * t4 + t2 * t35 + t3 * t34) + m(6) * (t12 * t86 + t129 * t199 + t13 * t85 + t153 * t168 + t32 * t51 + t33 * t50) + m(5) * (t123 * t110 + t122 * t111 + t71 * t163 + t70 * t164) + (-Ifges(6,4) * t196 - Ifges(6,2) * t195) * t383 + (-Ifges(6,1) * t196 - Ifges(6,4) * t195) * t384 + (-t112 * t50 + t113 * t51 - t12 * t195 + t13 * t196) * mrSges(6,3) + (Ifges(7,1) * t42 + Ifges(7,4) * t43) * t377 + (Ifges(7,4) * t42 + Ifges(7,2) * t43) * t381 + (Ifges(7,4) * t137 + Ifges(7,2) * t136) * t387 + t34 * t24 + t35 * t25 + t43 * t47 / 0.2e1 + t42 * t48 / 0.2e1 + (Ifges(7,1) * t137 + Ifges(7,4) * t136) * t388 + t78 * t49 + t4 * t81 + t5 * t82 + t85 * t68 + t86 * t69 + t103 * (-mrSges(7,1) * t43 + mrSges(7,2) * t42) + t112 * t91 / 0.2e1 + t32 * t132 + t33 * t133 + t52 * (-mrSges(7,1) * t136 + mrSges(7,2) * t137) + t149 * t8 + t153 * (-mrSges(6,1) * t113 + mrSges(6,2) * t112) + t163 * t155 + t164 * t156 + t168 * t104 + t110 * t179 + t111 * t180 + t199 * t39; t112 * t132 + t113 * t133 + t136 * t24 + t137 * t25 - t195 * t68 - t196 * t69 + t42 * t81 + t43 * t82 + (-t120 - t39 - t8) * t263 + (t274 * qJD(4) + t275) * t259 + ((t179 * t262 - t180 * t258 + t236 - t299) * t263 + (t104 + t49 - t298 + t314) * t259) * qJD(3) + m(7) * (t103 * t309 + t3 * t136 + t2 * t137 + t14 * t43 + t15 * t42 - t263 * t52) + m(4) * (t184 * t259 - t185 * t263 + (-t191 * t259 + t192 * t263) * qJD(3)) + m(6) * (t51 * t112 + t50 * t113 - t12 * t196 - t129 * t263 - t13 * t195 + t153 * t309) + m(5) * ((-qJD(3) * t276 - t185) * t263 + (qJD(3) * t182 - qJD(4) * t277 + t278) * t259); (-t12 * t273 - t13 * t222 - t429 * t51 + t436 * t50) * mrSges(6,3) + (mrSges(6,1) * t429 - mrSges(6,2) * t436) * t153 + (-Ifges(6,5) * t165 - Ifges(6,6) * t166) * t360 + (-t166 / 0.2e1 + t187 / 0.2e1) * t90 + (-Ifges(6,1) * t188 - Ifges(6,4) * t187) * t372 + (-Ifges(6,4) * t188 - Ifges(6,2) * t187) * t374 + (-Ifges(6,5) * t188 - Ifges(6,6) * t187) * t361 + (-t165 / 0.2e1 + t188 / 0.2e1) * t91 + (Ifges(6,4) * t222 - Ifges(6,2) * t273) * t383 + (Ifges(6,1) * t222 - Ifges(6,4) * t273) * t384 + t129 * (mrSges(6,1) * t273 + mrSges(6,2) * t222) - t273 * t36 / 0.2e1 + (-Ifges(6,1) * t165 - Ifges(6,4) * t166) * t371 + (-Ifges(6,4) * t165 - Ifges(6,2) * t166) * t373 + ((t292 - t334 + t420 + t410) * t263 + (-t264 + (t349 / 0.2e1 + (t421 + Ifges(4,2) / 0.2e1) * t263) * qJD(1) + t192 * mrSges(4,3) + t291) * t259 + (Ifges(6,5) * t222 + Ifges(7,5) * t162 - Ifges(6,6) * t273 + Ifges(7,6) * t161 + t279) * t309 / 0.2e1) * qJD(1) + t417 * t82 + (t412 * t103 + t417 * t14 + t418 * t15 + t189 * t52 + t2 * t84 + t3 * t83) * m(7) + t418 * t81 + t413 * t132 + (t12 * t173 + t129 * t251 + t13 * t172 - t153 * t169 + t413 * t51 + t414 * t50) * m(6) + t414 * t133 + t412 * t49 + (Ifges(7,5) * t128 + Ifges(7,6) * t127) * t363 + t280 * t369 + t282 * t370 + (Ifges(7,1) * t66 + Ifges(7,4) * t67) * t377 + t114 * t357 + (Ifges(7,5) * t66 + Ifges(7,6) * t67) * t362 + (t66 / 0.2e1 - t128 / 0.2e1) * t48 - m(5) * (t122 * t146 + t123 * t147 + t182 * t192) + m(5) * (-pkin(3) * t185 + pkin(8) * t278) + (t394 * t258 * pkin(4) + (-m(5) * t277 + t274) * pkin(8) + t265) * qJD(4) + (mrSges(7,1) * t327 - mrSges(7,2) * t326) * t103 + (t14 * t326 - t15 * t327 + t161 * t2 - t162 * t3) * mrSges(7,3) - t314 * t192 + (-mrSges(4,1) - t285) * t185 + t278 * mrSges(5,3) + t275 * pkin(8) + (Ifges(7,1) * t128 + Ifges(7,4) * t127) * t378 + (Ifges(7,4) * t66 + Ifges(7,2) * t67) * t381 + (Ifges(7,4) * t128 + Ifges(7,2) * t127) * t382 + (Ifges(7,4) * t162 + Ifges(7,2) * t161) * t387 + (Ifges(7,1) * t162 + Ifges(7,4) * t161) * t388 + t83 * t24 + t84 * t25 - pkin(3) * t120 + t161 * t6 / 0.2e1 + t162 * t7 / 0.2e1 + t52 * (-mrSges(7,1) * t161 + mrSges(7,2) * t162) + (t67 / 0.2e1 - t127 / 0.2e1) * t47 - t169 * t104 + t172 * t68 + t173 * t69 - t147 * t179 - t146 * t180 - t184 * mrSges(4,2) + t189 * t8 + t222 * t37 / 0.2e1 - t191 * t236 + t251 * t39 + t258 * t115 / 0.2e1; t160 * t405 + t415 * t81 + (-t103 * t124 + t416 * t14 + t415 * t15 + t2 * t201 + t200 * t3) * m(7) + t416 * t82 + t435 + t143 * t364 + (Ifges(5,1) * t219 - t337) * t365 + (Ifges(5,5) * t219 - Ifges(5,6) * t220) * t359 - m(6) * (t50 * t53 + t51 * t54) - t389 + t246 + (t257 * t69 + t261 * t68 + (t132 * t261 - t133 * t257) * qJD(5) + m(6) * (t12 * t257 + t13 * t261 + t303 * t51 - t304 * t50) - t394 * t220) * pkin(4) + (t122 * t219 + t123 * t220) * mrSges(5,3) + (-Ifges(5,2) * t220 + t144 + t214) * t366 - t124 * t49 - t54 * t132 - t53 * t133 - t122 * t179 + t123 * t180 + t200 * t24 + t201 * t25 - t182 * (mrSges(5,1) * t220 + mrSges(5,2) * t219); t90 * t371 + (-t160 * t49 + t260 * t24 + t256 * t25 + (-t256 * t82 + t260 * t81) * qJD(6) + (-t103 * t160 - t14 * t302 + t15 * t301 + t2 * t256 + t260 * t3) * m(7)) * pkin(5) - m(7) * (t14 * t16 + t15 * t17) - t17 * t81 - t16 * t82 - t50 * t132 + t51 * t133 + t435; -t14 * t81 + t15 * t82 + t437;];
tauc  = t1(:);
