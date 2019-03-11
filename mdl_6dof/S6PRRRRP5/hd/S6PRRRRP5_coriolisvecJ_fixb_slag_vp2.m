% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:15
% EndTime: 2019-03-09 00:20:53
% DurationCPUTime: 17.70s
% Computational Cost: add. (10831->717), mult. (29429->993), div. (0->0), fcn. (23296->12), ass. (0->313)
t426 = Ifges(7,4) + Ifges(6,4);
t427 = Ifges(7,1) + Ifges(6,1);
t425 = Ifges(7,5) + Ifges(6,5);
t424 = Ifges(7,2) + Ifges(6,2);
t423 = Ifges(7,6) + Ifges(6,6);
t249 = cos(pkin(7));
t239 = qJD(2) * t249 + qJD(3);
t252 = sin(qJ(4));
t256 = cos(qJ(4));
t253 = sin(qJ(3));
t247 = sin(pkin(7));
t330 = qJD(2) * t247;
t313 = t253 * t330;
t190 = t239 * t252 + t256 * t313;
t257 = cos(qJ(3));
t329 = qJD(2) * t257;
t312 = t247 * t329;
t232 = qJD(4) - t312;
t251 = sin(qJ(5));
t255 = cos(qJ(5));
t153 = -t190 * t251 + t232 * t255;
t189 = t239 * t256 - t252 * t313;
t184 = qJD(5) - t189;
t154 = t190 * t255 + t232 * t251;
t438 = t426 * t154;
t418 = t153 * t424 + t184 * t423 + t438;
t440 = -t418 / 0.2e1;
t439 = t426 * t153;
t417 = t154 * t427 + t425 * t184 + t439;
t340 = t256 * t257;
t181 = (t251 * t253 + t255 * t340) * t330;
t231 = -pkin(4) * t256 - pkin(11) * t252 - pkin(3);
t341 = t255 * t256;
t244 = pkin(10) * t341;
t322 = qJD(6) * t255;
t297 = pkin(4) * t252 - pkin(11) * t256;
t225 = t297 * qJD(4);
t326 = qJD(4) * t252;
t379 = pkin(10) * t251;
t333 = t255 * t225 + t326 * t379;
t353 = qJ(6) * t252;
t258 = cos(qJ(2));
t248 = sin(pkin(6));
t332 = qJD(1) * t248;
t228 = qJD(2) * pkin(2) + t258 * t332;
t347 = t249 * t253;
t254 = sin(qJ(2));
t314 = t254 * t332;
t220 = pkin(9) * t330 + t314;
t351 = t220 * t257;
t276 = t228 * t347 + t351;
t250 = cos(pkin(6));
t331 = qJD(1) * t250;
t109 = (t253 * t331 + t297 * t329) * t247 + t276;
t209 = t253 * t220;
t271 = t228 * t249 + t247 * t331;
t142 = t257 * t271 - t209;
t298 = pkin(3) * t253 - pkin(10) * t257;
t205 = t298 * t330;
t107 = t256 * t142 + t252 * t205;
t94 = pkin(11) * t313 + t107;
t38 = t255 * t109 - t251 * t94;
t380 = pkin(5) * t252;
t437 = qJ(6) * t181 - t312 * t380 - t38 - t252 * t322 + (-qJ(6) * t341 + t380) * qJD(4) + (-t244 + (-t231 + t353) * t251) * qJD(5) + t333;
t180 = (-t251 * t340 + t253 * t255) * t330;
t323 = qJD(5) * t255;
t334 = t251 * t225 + t231 * t323;
t345 = t252 * t255;
t39 = t251 * t109 + t255 * t94;
t436 = -qJ(6) * t180 + (-pkin(10) * qJD(4) - qJ(6) * qJD(5)) * t345 + (-qJD(6) * t252 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t256) * t251 + t334 - t39;
t325 = qJD(4) * t256;
t106 = -t252 * t142 + t205 * t256;
t93 = -pkin(4) * t313 - t106;
t435 = pkin(10) * t325 - t93 + (t251 * t325 + t252 * t323 + t180) * pkin(5);
t434 = t426 * t255;
t433 = t426 * t251;
t327 = qJD(3) * t257;
t310 = t252 * t327;
t160 = t239 * t326 + (t253 * t325 + t310) * t330;
t309 = t256 * t327;
t159 = t239 * t325 + (-t253 * t326 + t309) * t330;
t328 = qJD(3) * t247;
t308 = qJD(2) * t328;
t302 = t253 * t308;
t79 = qJD(5) * t153 + t159 * t255 + t251 * t302;
t80 = -qJD(5) * t154 - t159 * t251 + t255 * t302;
t420 = t160 * t423 + t424 * t80 + t426 * t79;
t419 = t160 * t425 + t426 * t80 + t427 * t79;
t432 = -t251 * t423 + t255 * t425;
t431 = -t251 * t424 + t434;
t430 = t255 * t427 - t433;
t429 = t313 / 0.2e1;
t428 = -t239 * Ifges(4,6) / 0.2e1;
t422 = Ifges(6,3) + Ifges(7,3);
t18 = Ifges(7,5) * t79 + Ifges(7,6) * t80 + Ifges(7,3) * t160;
t19 = Ifges(6,5) * t79 + Ifges(6,6) * t80 + Ifges(6,3) * t160;
t421 = t19 + t18;
t339 = t257 * t258;
t344 = t253 * t254;
t273 = -t249 * t344 + t339;
t176 = t273 * t332;
t306 = t247 * t314;
t149 = t176 * t256 + t252 * t306;
t349 = t247 * t257;
t218 = pkin(2) * t347 + pkin(9) * t349;
t200 = pkin(10) * t249 + t218;
t299 = -pkin(3) * t257 - pkin(10) * t253;
t201 = (-pkin(2) + t299) * t247;
t272 = t298 * qJD(3);
t206 = t247 * t272;
t350 = t247 * t253;
t240 = pkin(9) * t350;
t346 = t249 * t257;
t217 = pkin(2) * t346 - t240;
t207 = t217 * qJD(3);
t81 = -t200 * t326 + t201 * t325 + t252 * t206 + t256 * t207;
t416 = -t149 + t81;
t199 = t240 + (-pkin(2) * t257 - pkin(3)) * t249;
t213 = -t256 * t249 + t252 * t350;
t214 = t249 * t252 + t256 * t350;
t129 = pkin(4) * t213 - pkin(11) * t214 + t199;
t145 = t256 * t200 + t252 * t201;
t131 = -pkin(11) * t349 + t145;
t53 = t251 * t129 + t255 * t131;
t336 = mrSges(4,1) * t239 + mrSges(5,1) * t189 - mrSges(5,2) * t190 - mrSges(4,3) * t313;
t342 = t254 * t257;
t343 = t253 * t258;
t275 = t249 * t342 + t343;
t175 = t275 * t332;
t208 = t218 * qJD(3);
t415 = t208 - t175;
t197 = t251 * t231 + t244;
t183 = Ifges(5,4) * t189;
t369 = Ifges(5,5) * t232;
t375 = Ifges(5,1) * t190;
t118 = t183 + t369 + t375;
t125 = -pkin(3) * t239 - t142;
t143 = t253 * t271 + t351;
t126 = pkin(10) * t239 + t143;
t238 = t249 * t331;
t150 = t238 + (qJD(2) * t299 - t228) * t247;
t66 = -t252 * t126 + t150 * t256;
t414 = t66 * mrSges(5,3) - t118 / 0.2e1 - t369 / 0.2e1 - t125 * mrSges(5,2) - t183 / 0.2e1;
t413 = t251 * t425 + t255 * t423;
t412 = t255 * t424 + t433;
t411 = t251 * t427 + t434;
t410 = t249 * t339 - t344;
t164 = (t272 + t314) * t330;
t269 = t273 * qJD(2);
t315 = t250 * t349;
t303 = qJD(3) * t315;
t96 = (t228 * t346 - t209) * qJD(3) + (t248 * t269 + t303) * qJD(1);
t26 = -t126 * t326 + t150 * t325 + t252 * t164 + t256 * t96;
t16 = pkin(11) * t302 + t26;
t324 = qJD(5) * t251;
t270 = t275 * qJD(2);
t311 = t253 * t328;
t304 = t250 * t311;
t97 = t276 * qJD(3) + (t248 * t270 + t304) * qJD(1);
t43 = pkin(4) * t160 - pkin(11) * t159 + t97;
t67 = t126 * t256 + t150 * t252;
t55 = pkin(11) * t232 + t67;
t69 = -pkin(4) * t189 - pkin(11) * t190 + t125;
t3 = t255 * t16 + t251 * t43 + t69 * t323 - t324 * t55;
t25 = t251 * t69 + t255 * t55;
t4 = -qJD(5) * t25 - t16 * t251 + t255 * t43;
t409 = -t251 * t4 + t255 * t3;
t361 = t189 * Ifges(5,2);
t368 = Ifges(5,6) * t232;
t374 = Ifges(5,4) * t190;
t117 = t361 + t368 + t374;
t24 = -t251 * t55 + t255 * t69;
t13 = -qJ(6) * t154 + t24;
t12 = pkin(5) * t184 + t13;
t14 = qJ(6) * t153 + t25;
t317 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t318 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t319 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t58 = Ifges(7,5) * t154 + Ifges(7,6) * t153 + Ifges(7,3) * t184;
t59 = Ifges(6,5) * t154 + Ifges(6,6) * t153 + Ifges(6,3) * t184;
t262 = t318 * t153 - t319 * t154 - t317 * t184 + t14 * mrSges(7,2) + t25 * mrSges(6,2) + t67 * mrSges(5,3) + t117 / 0.2e1 - t58 / 0.2e1 - t59 / 0.2e1 + t374 / 0.2e1 + t368 / 0.2e1 - t12 * mrSges(7,1) - t125 * mrSges(5,1) - t24 * mrSges(6,1);
t408 = t262 + t361 / 0.2e1;
t27 = -t126 * t325 - t150 * t326 + t164 * t256 - t252 * t96;
t407 = -t27 * mrSges(5,1) + t26 * mrSges(5,2) - Ifges(5,5) * t159 + Ifges(5,6) * t160;
t278 = t24 * t255 + t25 * t251;
t292 = mrSges(7,1) * t251 + mrSges(7,2) * t255;
t294 = mrSges(6,1) * t251 + mrSges(6,2) * t255;
t382 = t255 / 0.2e1;
t385 = -t251 / 0.2e1;
t389 = t184 / 0.2e1;
t398 = t154 / 0.2e1;
t54 = -pkin(4) * t232 - t66;
t40 = -pkin(5) * t153 + qJD(6) + t54;
t400 = t153 / 0.2e1;
t406 = t278 * mrSges(6,3) + (t12 * t255 + t14 * t251) * mrSges(7,3) - t292 * t40 - t294 * t54 - t431 * t400 - t430 * t398 - t432 * t389 - t418 * t385 - t417 * t382;
t405 = t79 / 0.2e1;
t404 = t80 / 0.2e1;
t401 = -t153 / 0.2e1;
t399 = -t154 / 0.2e1;
t397 = t160 / 0.2e1;
t390 = -t184 / 0.2e1;
t388 = -t213 / 0.2e1;
t386 = t214 / 0.2e1;
t381 = pkin(5) * t251;
t376 = -qJ(6) - pkin(11);
t365 = t159 * Ifges(5,1);
t364 = t159 * Ifges(5,4);
t363 = t160 * Ifges(5,4);
t166 = -t248 * t410 - t315;
t362 = t166 * t97;
t136 = mrSges(5,1) * t302 - mrSges(5,3) * t159;
t31 = -mrSges(6,1) * t80 + mrSges(6,2) * t79;
t355 = t31 - t136;
t141 = pkin(4) * t190 - pkin(11) * t189;
t37 = t251 * t141 + t255 * t66;
t163 = mrSges(5,1) * t232 - mrSges(5,3) * t190;
t90 = -mrSges(6,1) * t153 + mrSges(6,2) * t154;
t354 = t90 - t163;
t352 = t189 * t251;
t348 = t248 * t254;
t110 = -mrSges(7,2) * t184 + mrSges(7,3) * t153;
t111 = -mrSges(6,2) * t184 + mrSges(6,3) * t153;
t338 = t110 + t111;
t112 = mrSges(7,1) * t184 - mrSges(7,3) * t154;
t113 = mrSges(6,1) * t184 - mrSges(6,3) * t154;
t337 = t112 + t113;
t89 = -mrSges(7,1) * t153 + mrSges(7,2) * t154;
t320 = -t89 - t354;
t30 = -t80 * mrSges(7,1) + t79 * mrSges(7,2);
t307 = qJD(5) * t376;
t36 = t255 * t141 - t251 * t66;
t52 = t255 * t129 - t131 * t251;
t144 = -t252 * t200 + t201 * t256;
t305 = t330 * t348;
t301 = -t97 * mrSges(4,1) - t96 * mrSges(4,2);
t130 = pkin(4) * t349 - t144;
t296 = -mrSges(4,1) * t257 + mrSges(4,2) * t253;
t295 = mrSges(6,1) * t255 - mrSges(6,2) * t251;
t293 = mrSges(7,1) * t255 - mrSges(7,2) * t251;
t274 = t249 * t343 + t342;
t167 = t248 * t274 + t250 * t350;
t212 = -t247 * t248 * t258 + t249 * t250;
t128 = t167 * t256 + t212 * t252;
t78 = t128 * t255 + t166 * t251;
t77 = -t128 * t251 + t166 * t255;
t127 = t167 * t252 - t212 * t256;
t82 = -t200 * t325 - t201 * t326 + t206 * t256 - t252 * t207;
t170 = -t214 * t251 - t255 * t349;
t277 = -t214 * t255 + t251 * t349;
t70 = pkin(11) * t311 + t81;
t168 = -qJD(4) * t213 + t247 * t309;
t169 = qJD(4) * t214 + t247 * t310;
t95 = pkin(4) * t169 - pkin(11) * t168 + t208;
t8 = t129 * t323 - t131 * t324 + t251 * t95 + t255 * t70;
t1 = pkin(5) * t160 - qJ(6) * t79 - qJD(6) * t154 + t4;
t2 = qJ(6) * t80 + qJD(6) * t153 + t3;
t268 = -t4 * mrSges(6,1) - t1 * mrSges(7,1) + t3 * mrSges(6,2) + t2 * mrSges(7,2);
t182 = -t228 * t247 + t238;
t236 = Ifges(4,4) * t312;
t266 = t182 * mrSges(4,2) + t239 * Ifges(4,5) - t142 * mrSges(4,3) + Ifges(4,1) * t429 + t236 / 0.2e1;
t71 = -pkin(4) * t311 - t82;
t9 = -qJD(5) * t53 - t251 * t70 + t255 * t95;
t17 = -pkin(4) * t302 - t27;
t265 = -t375 / 0.2e1 + t414;
t264 = t182 * mrSges(4,1) + t66 * mrSges(5,1) + t232 * Ifges(5,3) + t190 * Ifges(5,5) + t189 * Ifges(5,6) + t428 - (Ifges(4,4) * t253 + Ifges(4,2) * t257) * t330 / 0.2e1 - t143 * mrSges(4,3) - t67 * mrSges(5,2);
t246 = -pkin(5) * t255 - pkin(4);
t234 = t376 * t255;
t233 = t376 * t251;
t230 = Ifges(4,5) * t257 * t308;
t229 = Ifges(5,3) * t302;
t227 = (pkin(10) + t381) * t252;
t224 = t255 * t231;
t211 = -qJD(6) * t251 + t255 * t307;
t210 = t251 * t307 + t322;
t204 = t296 * t330;
t203 = -mrSges(4,2) * t239 + mrSges(4,3) * t312;
t196 = -t256 * t379 + t224;
t195 = (mrSges(4,1) * t253 + mrSges(4,2) * t257) * t308;
t172 = -t251 * t353 + t197;
t165 = -qJ(6) * t345 + t224 + (-pkin(5) - t379) * t256;
t162 = -mrSges(5,2) * t232 + mrSges(5,3) * t189;
t148 = t176 * t252 - t256 * t306;
t140 = -qJD(5) * t197 + t333;
t139 = (-t255 * t326 - t256 * t324) * pkin(10) + t334;
t137 = -mrSges(5,2) * t302 - mrSges(5,3) * t160;
t121 = t303 + (qJD(3) * t410 + t269) * t248;
t120 = t304 + (qJD(3) * t274 + t270) * t248;
t103 = t149 * t255 + t175 * t251;
t102 = -t149 * t251 + t175 * t255;
t101 = qJD(5) * t277 - t168 * t251 + t255 * t311;
t100 = qJD(5) * t170 + t168 * t255 + t251 * t311;
t99 = mrSges(5,1) * t160 + mrSges(5,2) * t159;
t85 = -pkin(5) * t170 + t130;
t84 = Ifges(5,5) * t302 - t363 + t365;
t83 = -t160 * Ifges(5,2) + Ifges(5,6) * t302 + t364;
t51 = pkin(5) * t352 + t67;
t50 = -mrSges(6,2) * t160 + mrSges(6,3) * t80;
t49 = -mrSges(7,2) * t160 + mrSges(7,3) * t80;
t48 = mrSges(6,1) * t160 - mrSges(6,3) * t79;
t47 = mrSges(7,1) * t160 - mrSges(7,3) * t79;
t46 = -qJD(4) * t127 + t121 * t256 + t252 * t305;
t45 = qJD(4) * t128 + t121 * t252 - t256 * t305;
t44 = qJ(6) * t170 + t53;
t35 = pkin(5) * t213 + qJ(6) * t277 + t52;
t34 = -pkin(5) * t101 + t71;
t32 = -qJ(6) * t352 + t37;
t28 = -qJ(6) * t189 * t255 + pkin(5) * t190 + t36;
t11 = qJD(5) * t77 + t120 * t251 + t255 * t46;
t10 = -qJD(5) * t78 + t120 * t255 - t251 * t46;
t7 = -pkin(5) * t80 + t17;
t6 = qJ(6) * t101 + qJD(6) * t170 + t8;
t5 = pkin(5) * t169 - qJ(6) * t100 + qJD(6) * t277 + t9;
t15 = [t121 * t203 + t128 * t137 + t46 * t162 + t166 * t99 + t212 * t195 + (t49 + t50) * t78 + (t47 + t48) * t77 + (-mrSges(3,1) * t254 - mrSges(3,2) * t258) * qJD(2) ^ 2 * t248 - t336 * t120 + t338 * t11 + t337 * t10 - t320 * t45 + (t30 + t355) * t127 + (t204 * t348 + (t166 * t257 - t167 * t253) * qJD(3) * mrSges(4,3)) * t330 + m(4) * (-t120 * t142 + t121 * t143 + t362 + t167 * t96 + (qJD(1) * t212 + t182) * t305) + m(5) * (t120 * t125 - t127 * t27 + t128 * t26 - t45 * t66 + t46 * t67 + t362) + m(7) * (t1 * t77 + t10 * t12 + t11 * t14 + t127 * t7 + t2 * t78 + t40 * t45) + m(6) * (t10 * t24 + t11 * t25 + t127 * t17 + t3 * t78 + t4 * t77 + t45 * t54); -t336 * t415 + t1 * (mrSges(7,1) * t213 + mrSges(7,3) * t277) + t4 * (mrSges(6,1) * t213 + mrSges(6,3) * t277) + t7 * (-mrSges(7,1) * t170 - mrSges(7,2) * t277) + t17 * (-mrSges(6,1) * t170 - mrSges(6,2) * t277) - t419 * t277 / 0.2e1 + (t170 * t423 + t213 * t422 - t277 * t425) * t397 + (t170 * t424 + t213 * t423 - t277 * t426) * t404 + (t170 * t426 + t213 * t425 - t277 * t427) * t405 - m(4) * (-t142 * t175 + t143 * t176 + t182 * t306) + (-t168 * t66 - t169 * t67 - t213 * t26 - t214 * t27) * mrSges(5,3) + ((-m(4) * pkin(2) + t296) * t306 + ((-t217 * mrSges(4,3) + Ifges(4,5) * t249 / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,4) * t349) * t257 + (-t218 * mrSges(4,3) - Ifges(4,6) * t249 + Ifges(5,5) * t386 + Ifges(5,6) * t388 - 0.3e1 / 0.2e1 * Ifges(4,4) * t350 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - Ifges(5,3) / 0.2e1) * t349) * t253) * qJD(3)) * t330 + (t1 * t35 + t2 * t44 + t7 * t85 + (-t148 + t34) * t40 + (-t103 + t6) * t14 + (-t102 + t5) * t12) * m(7) + t232 * (Ifges(5,5) * t168 - Ifges(5,6) * t169) / 0.2e1 - t160 * (Ifges(5,4) * t214 - Ifges(5,2) * t213) / 0.2e1 + t189 * (Ifges(5,4) * t168 - Ifges(5,2) * t169) / 0.2e1 + t320 * t148 + (t230 / 0.2e1 + t301) * t249 + m(4) * (-t142 * t208 + t143 * t207 - t217 * t97 + t218 * t96) + t159 * (Ifges(5,1) * t214 - Ifges(5,4) * t213) / 0.2e1 + t190 * (Ifges(5,1) * t168 - Ifges(5,4) * t169) / 0.2e1 + t84 * t386 + t83 * t388 + t97 * (mrSges(5,1) * t213 + mrSges(5,2) * t214) + t2 * (-mrSges(7,2) * t213 + mrSges(7,3) * t170) + t3 * (-mrSges(6,2) * t213 + mrSges(6,3) * t170) + t199 * t99 + t168 * t118 / 0.2e1 + t14 * (-mrSges(7,2) * t169 + mrSges(7,3) * t101) + t25 * (-mrSges(6,2) * t169 + mrSges(6,3) * t101) + t12 * (mrSges(7,1) * t169 - mrSges(7,3) * t100) + t24 * (mrSges(6,1) * t169 - mrSges(6,3) * t100) + t125 * (mrSges(5,1) * t169 + mrSges(5,2) * t168) - t169 * t117 / 0.2e1 + t82 * t163 + t144 * t136 + t145 * t137 + t130 * t31 + t6 * t110 + t8 * t111 + t5 * t112 + t9 * t113 + t40 * (-mrSges(7,1) * t101 + mrSges(7,2) * t100) + t54 * (-mrSges(6,1) * t101 + mrSges(6,2) * t100) + t34 * t89 + t71 * t90 + t85 * t30 - t337 * t102 - t338 * t103 + (t144 * t27 + t145 * t26 + t199 * t97 + t416 * t67 + (t148 + t82) * t66 + t415 * t125) * m(5) + t416 * t162 + (t130 * t17 + t3 * t53 + t4 * t52 + (-t148 + t71) * t54 + (-t103 + t8) * t25 + (-t102 + t9) * t24) * m(6) + t417 * t100 / 0.2e1 + t418 * t101 / 0.2e1 + t420 * t170 / 0.2e1 + t421 * t213 / 0.2e1 + (t100 * t425 + t101 * t423 + t169 * t422) * t389 + (t100 * t426 + t101 * t424 + t169 * t423) * t400 + (t100 * t427 + t101 * t426 + t169 * t425) * t398 + (-t204 * t314 + t97 * mrSges(4,3) * t253 - pkin(2) * t195 + (t96 * mrSges(4,3) - t229 / 0.2e1 + t407) * t257 + (t266 * t257 + (t428 + t264) * t253) * qJD(3)) * t247 + (t58 + t59) * t169 / 0.2e1 + (-t176 + t207) * t203 + t35 * t47 + t44 * t49 + t52 * t48 + t53 * t50; t435 * t89 + t436 * t110 + (-m(5) * t97 - t99) * pkin(3) + (t1 * t165 + t12 * t437 + t14 * t436 + t172 * t2 + t227 * t7 + t40 * t435) * m(7) + t437 * t112 - t417 * t181 / 0.2e1 + t230 + t301 + (t12 * t181 - t14 * t180) * mrSges(7,3) + (-t180 * t25 + t181 * t24) * mrSges(6,3) + (t140 - t38) * t113 + (((-m(5) * t67 - t162) * pkin(10) - t408) * t252 + (-t406 + (-m(5) * t66 + m(6) * t54 + t354) * pkin(10) - t265) * t256) * qJD(4) + (t139 - t39) * t111 + (Ifges(7,1) * t181 + Ifges(7,4) * t180) * t399 + (Ifges(6,1) * t181 + Ifges(6,4) * t180) * t399 + (Ifges(7,4) * t181 + Ifges(7,2) * t180) * t401 + (Ifges(6,4) * t181 + Ifges(6,2) * t180) * t401 + (Ifges(7,5) * t181 + Ifges(7,6) * t180) * t390 + (Ifges(6,5) * t181 + Ifges(6,6) * t180) * t390 + t227 * t30 - t142 * t203 + t196 * t48 + t197 * t50 - t40 * (-mrSges(7,1) * t180 + mrSges(7,2) * t181) - t54 * (-mrSges(6,1) * t180 + mrSges(6,2) * t181) + t172 * t49 + t165 * t47 - t107 * t162 - t106 * t163 - t93 * t90 + t336 * t143 - m(5) * (t106 * t66 + t107 * t67 + t125 * t143) - m(6) * (t24 * t38 + t25 * t39 + t54 * t93) + ((qJD(3) * (Ifges(5,5) * t252 + Ifges(5,6) * t256) / 0.2e1 + Ifges(4,4) * t429 + (-qJD(3) + t239 / 0.2e1) * Ifges(4,6) - t264) * t253 + (-t236 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t313 + t265 * t256 + t408 * t252 - t266) * t257) * t330 + (t84 / 0.2e1 - t363 / 0.2e1 + t365 / 0.2e1 - t27 * mrSges(5,3) + t97 * mrSges(5,2) + t7 * t292 + t17 * t294 + (-t1 * t255 - t2 * t251) * mrSges(7,3) + (-t251 * t3 - t255 * t4) * mrSges(6,3) + (-m(5) * t27 + m(6) * t17 + t355) * pkin(10) + (t40 * t293 + t54 * t295 + (t12 * t251 - t14 * t255) * mrSges(7,3) + (t24 * t251 - t25 * t255) * mrSges(6,3) + t412 * t401 + t411 * t399 + t413 * t390 + t255 * t440) * qJD(5) + t430 * t405 + t431 * t404 + t432 * t397 + (qJD(5) * t417 + t420) * t385 + t419 * t382) * t252 + m(6) * (t139 * t25 + t140 * t24 + t196 * t4 + t197 * t3) + (-t97 * mrSges(5,1) + t83 / 0.2e1 - t18 / 0.2e1 - t19 / 0.2e1 + t364 / 0.2e1 + t26 * mrSges(5,3) + t318 * t80 - t319 * t79 + (m(5) * t26 + t137) * pkin(10) + (-Ifges(5,2) / 0.2e1 - t317) * t160 + t268) * t256 + t180 * t440; t411 * t405 + t412 * t404 + t413 * t397 + ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t190 + t406 + t414) * t189 + (m(6) * t409 + (-m(6) * t278 - t251 * t111 - t255 * t113) * qJD(5) - t251 * t48 + t255 * t50) * pkin(11) + t409 * mrSges(6,3) - t407 + t229 + t262 * t190 - m(7) * (t12 * t28 + t14 * t32 + t40 * t51) - t354 * t67 + m(7) * (t1 * t233 + t12 * t211 + t14 * t210 - t2 * t234 + t246 * t7) + (-pkin(4) * t17 - t24 * t36 - t25 * t37 - t54 * t67) * m(6) + (t211 - t28) * t112 + (t210 - t32) * t110 + (-t406 + (m(7) * t40 + t89) * t381) * qJD(5) - t7 * t293 - t17 * t295 + (-t1 * t251 + t2 * t255) * mrSges(7,3) + t246 * t30 + t233 * t47 - t234 * t49 - t66 * t162 - t37 * t111 - t36 * t113 - t51 * t89 + t419 * t251 / 0.2e1 + t420 * t382 - pkin(4) * t31; t421 + (-(-t12 + t13) * t14 + (-t154 * t40 + t1) * pkin(5)) * m(7) - t268 + (-t154 * t89 + t47) * pkin(5) + (t12 * t153 + t14 * t154) * mrSges(7,3) + (t153 * t24 + t154 * t25) * mrSges(6,3) - t40 * (mrSges(7,1) * t154 + mrSges(7,2) * t153) - t54 * (mrSges(6,1) * t154 + mrSges(6,2) * t153) - t13 * t110 - t24 * t111 + t14 * t112 + t25 * t113 + (t153 * t427 - t438) * t399 + t418 * t398 + (t153 * t425 - t154 * t423) * t390 + (-t154 * t424 + t417 + t439) * t401; -t153 * t110 + t154 * t112 + 0.2e1 * (t7 / 0.2e1 + t12 * t398 + t14 * t401) * m(7) + t30;];
tauc  = t15(:);
