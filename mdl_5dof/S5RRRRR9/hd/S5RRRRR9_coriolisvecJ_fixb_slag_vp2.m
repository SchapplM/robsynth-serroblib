% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:44
% EndTime: 2019-12-31 22:28:13
% DurationCPUTime: 12.87s
% Computational Cost: add. (10236->599), mult. (25283->857), div. (0->0), fcn. (17460->8), ass. (0->281)
t250 = sin(qJ(2));
t254 = cos(qJ(2));
t274 = pkin(2) * t250 - pkin(7) * t254;
t213 = t274 * qJD(1);
t253 = cos(qJ(3));
t249 = sin(qJ(3));
t301 = qJD(1) * t250;
t285 = t249 * t301;
t169 = pkin(6) * t285 + t253 * t213;
t305 = t253 * t254;
t263 = pkin(3) * t250 - pkin(8) * t305;
t368 = -pkin(8) - pkin(7);
t286 = qJD(3) * t368;
t428 = -qJD(1) * t263 + t253 * t286 - t169;
t192 = t249 * t213;
t306 = t250 * t253;
t307 = t249 * t254;
t427 = t192 + (-pkin(6) * t306 - pkin(8) * t307) * qJD(1) - t249 * t286;
t298 = qJD(2) * t253;
t206 = -t285 + t298;
t223 = -qJD(2) * pkin(2) + pkin(6) * t301;
t173 = -pkin(3) * t206 + t223;
t284 = t253 * t301;
t207 = qJD(2) * t249 + t284;
t248 = sin(qJ(4));
t252 = cos(qJ(4));
t277 = t252 * t206 - t207 * t248;
t108 = -pkin(4) * t277 + t173;
t251 = cos(qJ(5));
t247 = sin(qJ(5));
t395 = pkin(9) * t277;
t218 = -pkin(2) * t254 - t250 * pkin(7) - pkin(1);
t198 = t218 * qJD(1);
t300 = qJD(1) * t254;
t244 = pkin(6) * t300;
t224 = qJD(2) * pkin(7) + t244;
t155 = t253 * t198 - t224 * t249;
t120 = -pkin(8) * t207 + t155;
t237 = qJD(3) - t300;
t110 = pkin(3) * t237 + t120;
t156 = t198 * t249 + t224 * t253;
t121 = pkin(8) * t206 + t156;
t115 = t252 * t121;
t52 = t110 * t248 + t115;
t47 = t52 + t395;
t319 = t247 * t47;
t227 = qJD(4) + t237;
t150 = t206 * t248 + t207 * t252;
t410 = pkin(9) * t150;
t113 = t248 * t121;
t51 = t252 * t110 - t113;
t46 = t51 - t410;
t44 = pkin(4) * t227 + t46;
t12 = t251 * t44 - t319;
t318 = t251 * t47;
t13 = t247 * t44 + t318;
t299 = qJD(2) * t250;
t279 = qJD(1) * t299;
t399 = -t150 * t247 + t251 * t277;
t290 = qJD(2) * qJD(3);
t296 = qJD(3) * t249;
t297 = qJD(2) * t254;
t167 = t253 * t290 + (-t250 * t296 + t253 * t297) * qJD(1);
t295 = qJD(3) * t253;
t402 = t249 * t297 + t250 * t295;
t168 = -qJD(1) * t402 - t249 * t290;
t67 = qJD(4) * t277 + t167 * t252 + t168 * t248;
t68 = -qJD(4) * t150 - t167 * t248 + t168 * t252;
t28 = qJD(5) * t399 + t247 * t68 + t251 * t67;
t90 = t150 * t251 + t247 * t277;
t29 = -qJD(5) * t90 - t247 * t67 + t251 * t68;
t289 = Ifges(6,5) * t28 + Ifges(6,6) * t29 + Ifges(6,3) * t279;
t342 = Ifges(6,4) * t90;
t219 = qJD(5) + t227;
t351 = -t219 / 0.2e1;
t369 = t90 / 0.2e1;
t370 = -t90 / 0.2e1;
t372 = -t399 / 0.2e1;
t293 = qJD(4) * t252;
t294 = qJD(4) * t248;
t216 = t274 * qJD(2);
t199 = qJD(1) * t216;
t276 = pkin(6) * t279;
t100 = -qJD(3) * t156 + t253 * t199 + t249 * t276;
t63 = pkin(3) * t279 - pkin(8) * t167 + t100;
t99 = t198 * t295 + t249 * t199 - t224 * t296 - t253 * t276;
t72 = pkin(8) * t168 + t99;
t14 = t110 * t293 - t121 * t294 + t248 * t63 + t252 * t72;
t10 = pkin(9) * t68 + t14;
t15 = -qJD(4) * t52 - t248 * t72 + t252 * t63;
t9 = pkin(4) * t279 - pkin(9) * t67 + t15;
t2 = qJD(5) * t12 + t10 * t251 + t247 * t9;
t3 = -qJD(5) * t13 - t10 * t247 + t251 * t9;
t394 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t42 = Ifges(6,2) * t399 + Ifges(6,6) * t219 + t342;
t84 = Ifges(6,4) * t399;
t43 = Ifges(6,1) * t90 + Ifges(6,5) * t219 + t84;
t426 = (Ifges(6,1) * t399 - t342) * t370 + (Ifges(6,5) * t399 - Ifges(6,6) * t90) * t351 + (t12 * t399 + t13 * t90) * mrSges(6,3) - t108 * (mrSges(6,1) * t90 + mrSges(6,2) * t399) + t42 * t369 + t289 + t394 + (-Ifges(6,2) * t90 + t43 + t84) * t372;
t264 = t248 * t249 - t252 * t253;
t379 = qJD(3) + qJD(4);
t159 = t379 * t264;
t260 = t264 * t254;
t178 = qJD(1) * t260;
t425 = t159 - t178;
t209 = t248 * t253 + t249 * t252;
t160 = t379 * t209;
t261 = t209 * t254;
t177 = qJD(1) * t261;
t418 = t160 - t177;
t143 = Ifges(5,4) * t277;
t288 = Ifges(5,5) * t67 + Ifges(5,6) * t68 + Ifges(5,3) * t279;
t331 = Ifges(5,4) * t150;
t349 = -t227 / 0.2e1;
t360 = -t150 / 0.2e1;
t362 = -t277 / 0.2e1;
t385 = t15 * mrSges(5,1) - t14 * mrSges(5,2);
t81 = Ifges(5,1) * t150 + Ifges(5,5) * t227 + t143;
t424 = t288 + t385 + (Ifges(5,5) * t277 - Ifges(5,6) * t150) * t349 + (t150 * t52 + t277 * t51) * mrSges(5,3) + (-Ifges(5,2) * t150 + t143 + t81) * t362 - t173 * (mrSges(5,1) * t150 + mrSges(5,2) * t277) + (Ifges(5,1) * t277 - t331) * t360 + t426;
t225 = t368 * t249;
t226 = t368 * t253;
t407 = t225 * t293 + t226 * t294 + t428 * t248 - t427 * t252;
t166 = t248 * t225 - t252 * t226;
t406 = -qJD(4) * t166 + t427 * t248 + t428 * t252;
t421 = -pkin(4) * t301 + pkin(9) * t425 + t406;
t420 = pkin(9) * t418 - t407;
t281 = Ifges(3,5) * qJD(2) / 0.2e1;
t165 = t252 * t225 + t226 * t248;
t128 = -pkin(9) * t209 + t165;
t129 = -pkin(9) * t264 + t166;
t75 = t128 * t251 - t129 * t247;
t409 = qJD(5) * t75 + t421 * t247 - t420 * t251;
t76 = t128 * t247 + t129 * t251;
t408 = -qJD(5) * t76 + t420 * t247 + t421 * t251;
t240 = pkin(3) * t252 + pkin(4);
t291 = qJD(5) * t251;
t292 = qJD(5) * t247;
t310 = t247 * t248;
t59 = -t120 * t248 - t115;
t48 = t59 - t395;
t60 = t252 * t120 - t113;
t49 = t60 - t410;
t405 = -t247 * t48 - t251 * t49 + t240 * t291 + (-t248 * t292 + (t251 * t252 - t310) * qJD(4)) * pkin(3);
t309 = t248 * t251;
t404 = t247 * t49 - t251 * t48 - t240 * t292 + (-t248 * t291 + (-t247 * t252 - t309) * qJD(4)) * pkin(3);
t341 = pkin(3) * t249;
t201 = t300 * t341 + t244;
t403 = pkin(3) * t296 + t418 * pkin(4) - t201;
t242 = Ifges(3,4) * t300;
t324 = t207 * Ifges(4,4);
t131 = t206 * Ifges(4,2) + t237 * Ifges(4,6) + t324;
t202 = Ifges(4,4) * t206;
t132 = t207 * Ifges(4,1) + t237 * Ifges(4,5) + t202;
t265 = t155 * t253 + t156 * t249;
t332 = Ifges(4,4) * t253;
t269 = -Ifges(4,2) * t249 + t332;
t333 = Ifges(4,4) * t249;
t271 = Ifges(4,1) * t253 - t333;
t272 = mrSges(4,1) * t249 + mrSges(4,2) * t253;
t329 = Ifges(4,6) * t249;
t330 = Ifges(4,5) * t253;
t345 = t253 / 0.2e1;
t346 = -t249 / 0.2e1;
t352 = t207 / 0.2e1;
t256 = -t265 * mrSges(4,3) + t223 * t272 + t206 * t269 / 0.2e1 + t271 * t352 + t237 * (-t329 + t330) / 0.2e1 + t131 * t346 + t132 * t345;
t401 = t256 + Ifges(3,1) * t301 / 0.2e1 + t242 / 0.2e1 + t281;
t80 = Ifges(5,2) * t277 + Ifges(5,6) * t227 + t331;
t396 = t80 / 0.2e1;
t280 = -Ifges(3,6) * qJD(2) / 0.2e1;
t93 = -mrSges(5,1) * t277 + mrSges(5,2) * t150;
t384 = m(5) * t173 + t93;
t185 = t264 * t250;
t205 = t253 * t218;
t340 = pkin(6) * t249;
t154 = -pkin(8) * t306 + t205 + (-pkin(3) - t340) * t254;
t239 = pkin(6) * t305;
t176 = t249 * t218 + t239;
t308 = t249 * t250;
t162 = -pkin(8) * t308 + t176;
t96 = t248 * t154 + t252 * t162;
t382 = Ifges(4,5) * t167 + Ifges(4,6) * t168;
t381 = -t100 * mrSges(4,1) + t99 * mrSges(4,2);
t380 = pkin(1) * mrSges(3,2) * qJD(1);
t378 = t28 / 0.2e1;
t377 = t29 / 0.2e1;
t374 = t67 / 0.2e1;
t373 = t68 / 0.2e1;
t371 = t399 / 0.2e1;
t367 = pkin(1) * mrSges(3,1);
t184 = t209 * t250;
t124 = -t184 * t251 + t185 * t247;
t364 = t124 / 0.2e1;
t125 = -t184 * t247 - t185 * t251;
t363 = t125 / 0.2e1;
t361 = t277 / 0.2e1;
t359 = t150 / 0.2e1;
t358 = t167 / 0.2e1;
t357 = t168 / 0.2e1;
t356 = -t184 / 0.2e1;
t355 = -t185 / 0.2e1;
t354 = -t206 / 0.2e1;
t353 = -t207 / 0.2e1;
t350 = t219 / 0.2e1;
t348 = t227 / 0.2e1;
t347 = -t237 / 0.2e1;
t334 = Ifges(3,4) * t250;
t116 = -t177 * t251 + t178 * t247;
t152 = t209 * t251 - t247 * t264;
t55 = -qJD(5) * t152 + t159 * t247 - t160 * t251;
t317 = t116 - t55;
t117 = -t177 * t247 - t178 * t251;
t151 = -t209 * t247 - t251 * t264;
t54 = qJD(5) * t151 - t159 * t251 - t160 * t247;
t316 = t117 - t54;
t313 = qJD(2) * mrSges(3,2);
t302 = t253 * t216 + t299 * t340;
t217 = pkin(3) * t308 + t250 * pkin(6);
t287 = Ifges(4,3) * t279 + t382;
t174 = pkin(3) * t402 + pkin(6) * t297;
t241 = -pkin(3) * t253 - pkin(2);
t146 = -pkin(3) * t168 + qJD(2) * t244;
t95 = t252 * t154 - t248 * t162;
t275 = m(4) * t223 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t206 + mrSges(4,2) * t207 + mrSges(3,3) * t301;
t273 = mrSges(4,1) * t253 - mrSges(4,2) * t249;
t270 = Ifges(4,1) * t249 + t332;
t268 = Ifges(4,2) * t253 + t333;
t267 = Ifges(4,5) * t249 + Ifges(4,6) * t253;
t71 = -pkin(4) * t254 + t185 * pkin(9) + t95;
t77 = -pkin(9) * t184 + t96;
t37 = -t247 * t77 + t251 * t71;
t38 = t247 * t71 + t251 * t77;
t266 = -t100 * t249 + t253 * t99;
t94 = t263 * qJD(2) + (-t239 + (pkin(8) * t250 - t218) * t249) * qJD(3) + t302;
t118 = t249 * t216 + t218 * t295 + (-t250 * t298 - t254 * t296) * pkin(6);
t98 = -pkin(8) * t402 + t118;
t32 = t154 * t293 - t162 * t294 + t248 * t94 + t252 * t98;
t33 = -qJD(4) * t96 - t248 * t98 + t252 * t94;
t255 = t12 * mrSges(6,1) + t155 * mrSges(4,1) + t51 * mrSges(5,1) + t237 * Ifges(4,3) + t207 * Ifges(4,5) + t206 * Ifges(4,6) + t280 - (t254 * Ifges(3,2) + t334) * qJD(1) / 0.2e1 + t219 * Ifges(6,3) + t90 * Ifges(6,5) + t399 * Ifges(6,6) + t227 * Ifges(5,3) + t150 * Ifges(5,5) + t277 * Ifges(5,6) - t13 * mrSges(6,2) - t156 * mrSges(4,2) - t52 * mrSges(5,2);
t221 = mrSges(3,3) * t300 - t313;
t187 = pkin(3) * t309 + t240 * t247;
t186 = -pkin(3) * t310 + t240 * t251;
t179 = pkin(4) * t264 + t241;
t175 = -pkin(6) * t307 + t205;
t172 = mrSges(4,1) * t237 - mrSges(4,3) * t207;
t171 = -mrSges(4,2) * t237 + mrSges(4,3) * t206;
t170 = -pkin(6) * t284 + t192;
t157 = pkin(4) * t184 + t217;
t145 = -mrSges(4,2) * t279 + mrSges(4,3) * t168;
t144 = mrSges(4,1) * t279 - mrSges(4,3) * t167;
t123 = mrSges(5,1) * t227 - mrSges(5,3) * t150;
t122 = -mrSges(5,2) * t227 + mrSges(5,3) * t277;
t119 = -qJD(3) * t176 + t302;
t112 = pkin(3) * t207 + pkin(4) * t150;
t111 = -mrSges(4,1) * t168 + mrSges(4,2) * t167;
t104 = t167 * Ifges(4,1) + t168 * Ifges(4,4) + Ifges(4,5) * t279;
t103 = t167 * Ifges(4,4) + t168 * Ifges(4,2) + Ifges(4,6) * t279;
t102 = -qJD(2) * t261 + t185 * t379;
t101 = -qJD(2) * t260 - t160 * t250;
t78 = -pkin(4) * t102 + t174;
t74 = mrSges(6,1) * t219 - mrSges(6,3) * t90;
t73 = -mrSges(6,2) * t219 + mrSges(6,3) * t399;
t57 = -mrSges(5,2) * t279 + mrSges(5,3) * t68;
t56 = mrSges(5,1) * t279 - mrSges(5,3) * t67;
t50 = -pkin(4) * t68 + t146;
t45 = -mrSges(6,1) * t399 + mrSges(6,2) * t90;
t40 = -qJD(5) * t125 - t101 * t247 + t102 * t251;
t39 = qJD(5) * t124 + t101 * t251 + t102 * t247;
t36 = -mrSges(5,1) * t68 + mrSges(5,2) * t67;
t35 = t67 * Ifges(5,1) + t68 * Ifges(5,4) + Ifges(5,5) * t279;
t34 = t67 * Ifges(5,4) + t68 * Ifges(5,2) + Ifges(5,6) * t279;
t25 = pkin(9) * t102 + t32;
t24 = -mrSges(6,2) * t279 + mrSges(6,3) * t29;
t23 = mrSges(6,1) * t279 - mrSges(6,3) * t28;
t20 = pkin(4) * t299 - pkin(9) * t101 + t33;
t17 = t251 * t46 - t319;
t16 = -t247 * t46 - t318;
t8 = -mrSges(6,1) * t29 + mrSges(6,2) * t28;
t7 = t28 * Ifges(6,1) + t29 * Ifges(6,4) + Ifges(6,5) * t279;
t6 = t28 * Ifges(6,4) + t29 * Ifges(6,2) + Ifges(6,6) * t279;
t5 = -qJD(5) * t38 + t20 * t251 - t247 * t25;
t4 = qJD(5) * t37 + t20 * t247 + t25 * t251;
t1 = [(-Ifges(5,1) * t185 - Ifges(5,4) * t184) * t374 + (-t101 * t51 + t102 * t52 - t14 * t184 + t15 * t185) * mrSges(5,3) + t146 * (mrSges(5,1) * t184 - mrSges(5,2) * t185) + (-Ifges(5,4) * t185 - Ifges(5,2) * t184) * t373 + t102 * t396 + (pkin(6) * t111 + t103 * t346 + t104 * t345 + t271 * t358 + t269 * t357 + (-t100 * t253 - t249 * t99) * mrSges(4,3) + (t267 * t347 + t223 * t273 + t268 * t354 + t270 * t353 + t132 * t346 - t253 * t131 / 0.2e1 + (t155 * t249 - t156 * t253) * mrSges(4,3)) * qJD(3) + (t255 - pkin(6) * t221 + (-0.2e1 * t367 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t330 / 0.2e1 - t329 / 0.2e1) * t250 + Ifges(5,5) * t355 + Ifges(5,6) * t356 + Ifges(6,5) * t363 + Ifges(6,6) * t364) * qJD(1) + t280) * qJD(2)) * t250 - (t289 + t288 + t287 + t382) * t254 / 0.2e1 + (-Ifges(5,5) * t374 - Ifges(6,5) * t378 - Ifges(5,6) * t373 - Ifges(6,6) * t377 + (0.3e1 / 0.2e1 * Ifges(3,4) * t297 + (-0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + (m(4) * pkin(6) + t272) * pkin(6)) * t299) * qJD(1) + t381 - t385 - t394) * t254 + (t275 * pkin(6) + t281 - 0.2e1 * t380 + t401) * t297 + m(4) * (t100 * t175 + t156 * t118 + t155 * t119 + t99 * t176) + (-t12 * t39 + t124 * t2 - t125 * t3 + t13 * t40) * mrSges(6,3) + (Ifges(6,1) * t125 + Ifges(6,4) * t124) * t378 + t217 * t36 + t175 * t144 + t176 * t145 + t118 * t171 + t119 * t172 + t173 * (-mrSges(5,1) * t102 + mrSges(5,2) * t101) + t174 * t93 + t157 * t8 + t50 * (-mrSges(6,1) * t124 + mrSges(6,2) * t125) + t32 * t122 + t33 * t123 + t108 * (-mrSges(6,1) * t40 + mrSges(6,2) * t39) + t101 * t81 / 0.2e1 + t95 * t56 + t96 * t57 + (Ifges(5,5) * t101 + Ifges(5,6) * t102) * t348 + (Ifges(6,5) * t39 + Ifges(6,6) * t40) * t350 + t35 * t355 + t34 * t356 + (Ifges(5,1) * t101 + Ifges(5,4) * t102) * t359 + (Ifges(5,4) * t101 + Ifges(5,2) * t102) * t361 + t7 * t363 + t6 * t364 + (Ifges(6,1) * t39 + Ifges(6,4) * t40) * t369 + (Ifges(6,4) * t39 + Ifges(6,2) * t40) * t371 + m(5) * (t14 * t96 + t146 * t217 + t15 * t95 + t173 * t174 + t32 * t52 + t33 * t51) + m(6) * (t108 * t78 + t12 * t5 + t13 * t4 + t157 * t50 + t2 * t38 + t3 * t37) + (Ifges(6,4) * t125 + Ifges(6,2) * t124) * t377 + t37 * t23 + t38 * t24 + t40 * t42 / 0.2e1 + t39 * t43 / 0.2e1 + t4 * t73 + t5 * t74 + t78 * t45; t406 * t123 + (t14 * t166 + t146 * t241 + t15 * t165 - t173 * t201 + t406 * t51 + t407 * t52) * m(5) + t407 * t122 + (t384 * t341 + t256) * qJD(3) + t103 * t345 + t266 * mrSges(4,3) + (t12 * t316 - t13 * t317 + t151 * t2 - t152 * t3) * mrSges(6,3) + (mrSges(6,1) * t317 - mrSges(6,2) * t316) * t108 + (-t116 / 0.2e1 + t55 / 0.2e1) * t42 + (m(4) * t266 - t144 * t249 + t145 * t253 + (-m(4) * t265 - t249 * t171 - t253 * t172) * qJD(3)) * pkin(7) + (-t14 * t264 - t15 * t209 - t418 * t52 + t425 * t51) * mrSges(5,3) + (mrSges(5,1) * t418 - mrSges(5,2) * t425) * t173 + ((t281 - t242 / 0.2e1 + t380 + ((-m(4) * pkin(2) - mrSges(3,1) - t273) * qJD(2) - t275) * pkin(6) - t401) * t254 + (-t255 + t280 + (t367 + t334 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t254) * qJD(1) + (t221 + t313) * pkin(6)) * t250 + (Ifges(5,5) * t209 + Ifges(6,5) * t152 - Ifges(5,6) * t264 + Ifges(6,6) * t151 + t267) * t299 / 0.2e1) * qJD(1) + t403 * t45 + t146 * (mrSges(5,1) * t264 + mrSges(5,2) * t209) + (Ifges(5,4) * t209 - Ifges(5,2) * t264) * t373 + (Ifges(5,1) * t209 - Ifges(5,4) * t264) * t374 - t264 * t34 / 0.2e1 + (-Ifges(5,5) * t159 - Ifges(5,6) * t160) * t348 + (-Ifges(5,1) * t159 - Ifges(5,4) * t160) * t359 + (-Ifges(5,4) * t159 - Ifges(5,2) * t160) * t361 + (t177 / 0.2e1 - t160 / 0.2e1) * t80 + (t178 / 0.2e1 - t159 / 0.2e1) * t81 + (-Ifges(5,5) * t178 - Ifges(5,6) * t177) * t349 + (-Ifges(5,1) * t178 - Ifges(5,4) * t177) * t360 + (-Ifges(5,4) * t178 - Ifges(5,2) * t177) * t362 + t408 * t74 + (t108 * t403 + t12 * t408 + t13 * t409 + t179 * t50 + t2 * t76 + t3 * t75) * m(6) + t409 * t73 - m(4) * (t155 * t169 + t156 * t170) + (-t117 / 0.2e1 + t54 / 0.2e1) * t43 + t249 * t104 / 0.2e1 + t241 * t36 + t209 * t35 / 0.2e1 - t201 * t93 + t179 * t8 - t170 * t171 - t169 * t172 + t165 * t56 + t166 * t57 + t151 * t6 / 0.2e1 + t50 * (-mrSges(6,1) * t151 + mrSges(6,2) * t152) + t152 * t7 / 0.2e1 - pkin(2) * t111 + (Ifges(6,5) * t54 + Ifges(6,6) * t55) * t350 + (Ifges(6,5) * t117 + Ifges(6,6) * t116) * t351 + t268 * t357 + t270 * t358 + (Ifges(6,1) * t54 + Ifges(6,4) * t55) * t369 + (Ifges(6,1) * t117 + Ifges(6,4) * t116) * t370 + (Ifges(6,4) * t54 + Ifges(6,2) * t55) * t371 + (Ifges(6,4) * t117 + Ifges(6,2) * t116) * t372 + (Ifges(6,4) * t152 + Ifges(6,2) * t151) * t377 + (Ifges(6,1) * t152 + Ifges(6,4) * t151) * t378 + t75 * t23 + t76 * t24; t150 * t396 + t404 * t74 + (-t108 * t112 + t12 * t404 + t13 * t405 + t186 * t3 + t187 * t2) * m(6) + t405 * t73 + t424 + (Ifges(4,5) * t206 - Ifges(4,6) * t207) * t347 + t287 - t381 + (-Ifges(4,2) * t207 + t132 + t202) * t354 - m(5) * (t51 * t59 + t52 * t60) - t223 * (mrSges(4,1) * t207 + mrSges(4,2) * t206) + t186 * t23 + t187 * t24 - t155 * t171 + t156 * t172 - t60 * t122 - t59 * t123 - t112 * t45 + t131 * t352 + (t155 * t206 + t156 * t207) * mrSges(4,3) + (Ifges(4,1) * t206 - t324) * t353 + (t248 * t57 + t252 * t56 + (t122 * t252 - t123 * t248) * qJD(4) + m(5) * (t14 * t248 + t15 * t252 + t293 * t52 - t294 * t51) - t384 * t207) * pkin(3); (-t150 * t45 + t251 * t23 + t247 * t24 + (-t247 * t74 + t251 * t73) * qJD(5) + (-t108 * t150 - t12 * t292 + t13 * t291 + t2 * t247 + t251 * t3) * m(6)) * pkin(4) - m(6) * (t12 * t16 + t13 * t17) - t51 * t122 + t52 * t123 + t80 * t359 - t17 * t73 - t16 * t74 + t424; -t12 * t73 + t13 * t74 + t426;];
tauc = t1(:);
