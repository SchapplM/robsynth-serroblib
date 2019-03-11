% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:01
% EndTime: 2019-03-09 02:23:20
% DurationCPUTime: 14.83s
% Computational Cost: add. (6470->600), mult. (12251->809), div. (0->0), fcn. (7437->14), ass. (0->278)
t194 = sin(qJ(4));
t184 = t194 * qJD(1);
t168 = t184 + qJD(5);
t166 = qJD(6) + t168;
t329 = -t166 / 0.2e1;
t193 = sin(qJ(5));
t197 = cos(qJ(5));
t282 = qJD(4) * t197;
t198 = cos(qJ(4));
t286 = qJD(1) * t198;
t149 = -t193 * t286 + t282;
t284 = qJD(4) * t193;
t150 = t197 * t286 + t284;
t192 = sin(qJ(6));
t196 = cos(qJ(6));
t80 = t149 * t192 + t150 * t196;
t338 = -t80 / 0.2e1;
t244 = t196 * t149 - t150 * t192;
t340 = -t244 / 0.2e1;
t398 = Ifges(7,5) * t338 + Ifges(7,6) * t340 + Ifges(7,3) * t329;
t397 = -t149 / 0.2e1;
t396 = -t150 / 0.2e1;
t395 = -t168 / 0.2e1;
t177 = pkin(5) * t197 + pkin(4);
t394 = m(7) * t177;
t200 = -pkin(9) - pkin(8);
t393 = -m(7) * t200 + mrSges(6,3) + mrSges(7,3);
t275 = qJD(1) * qJD(4);
t157 = qJDD(1) * t198 - t194 * t275;
t392 = t157 / 0.2e1;
t158 = -t194 * qJDD(1) - t198 * t275;
t391 = t158 / 0.2e1;
t258 = t193 * t184;
t259 = qJD(5) * t200;
t191 = cos(pkin(10));
t174 = -pkin(1) * t191 - pkin(2);
t167 = -pkin(7) + t174;
t144 = qJD(1) * t167 + qJD(3);
t105 = -t194 * qJD(2) + t144 * t198;
t242 = pkin(4) * t198 + pkin(8) * t194;
t154 = t242 * qJD(1);
t65 = t197 * t105 + t193 * t154;
t390 = -pkin(9) * t258 + t193 * t259 - t65;
t288 = t194 * t197;
t271 = pkin(9) * t288;
t64 = -t105 * t193 + t197 * t154;
t389 = t197 * t259 - (pkin(5) * t198 + t271) * qJD(1) - t64;
t188 = qJ(1) + pkin(10);
t178 = sin(t188);
t179 = cos(t188);
t388 = g(1) * t178 - g(2) * t179;
t264 = mrSges(5,3) * t286;
t367 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t149 + mrSges(6,2) * t150 + t264;
t90 = -qJD(4) * pkin(4) - t105;
t387 = -m(6) * t90 - t367;
t190 = sin(pkin(10));
t171 = pkin(1) * t190 + qJ(3);
t241 = pkin(4) * t194 - pkin(8) * t198;
t139 = t171 + t241;
t107 = t139 * qJD(1);
t279 = qJD(5) * t197;
t280 = qJD(5) * t193;
t281 = qJD(4) * t198;
t378 = -qJD(2) * qJD(4) + qJDD(1) * t167 + qJDD(3);
t57 = t198 * qJDD(2) + t144 * t281 + t194 * t378;
t54 = qJDD(4) * pkin(8) + t57;
t145 = qJD(1) * qJD(3) + qJDD(1) * t171;
t66 = -pkin(4) * t158 - pkin(8) * t157 + t145;
t285 = qJD(2) * t198;
t106 = t144 * t194 + t285;
t91 = qJD(4) * pkin(8) + t106;
t11 = t107 * t279 + t193 * t66 + t197 * t54 - t280 * t91;
t75 = -qJD(5) * t150 + qJDD(4) * t197 - t157 * t193;
t10 = pkin(9) * t75 + t11;
t51 = t107 * t193 + t197 * t91;
t39 = pkin(9) * t149 + t51;
t303 = t192 * t39;
t50 = t197 * t107 - t193 * t91;
t38 = -pkin(9) * t150 + t50;
t35 = pkin(5) * t168 + t38;
t13 = t196 * t35 - t303;
t12 = -qJD(5) * t51 - t193 * t54 + t197 * t66;
t146 = qJDD(5) - t158;
t74 = qJD(5) * t149 + qJDD(4) * t193 + t157 * t197;
t9 = pkin(5) * t146 - pkin(9) * t74 + t12;
t2 = qJD(6) * t13 + t10 * t196 + t192 * t9;
t300 = t196 * t39;
t14 = t192 * t35 + t300;
t3 = -qJD(6) * t14 - t10 * t192 + t196 * t9;
t386 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t385 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t310 = Ifges(5,4) * t198;
t231 = -Ifges(5,2) * t194 + t310;
t384 = t14 * mrSges(7,2) + Ifges(5,6) * qJD(4) / 0.2e1 + qJD(1) * t231 / 0.2e1 + Ifges(6,5) * t396 + Ifges(6,6) * t397 + Ifges(6,3) * t395 - t13 * mrSges(7,1) + t398;
t351 = m(7) * pkin(5);
t383 = m(6) + m(7);
t322 = pkin(5) * t193;
t382 = m(7) * t322;
t36 = -mrSges(6,1) * t75 + mrSges(6,2) * t74;
t380 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t157 - t36;
t222 = t192 * t193 - t196 * t197;
t357 = qJD(5) + qJD(6);
t379 = t357 * t222;
t225 = t11 * t197 - t12 * t193;
t377 = -t50 * t279 - t51 * t280 + t225;
t311 = Ifges(5,4) * t194;
t234 = t198 * Ifges(5,1) - t311;
t142 = Ifges(6,4) * t149;
t71 = t150 * Ifges(6,1) + t168 * Ifges(6,5) + t142;
t376 = Ifges(5,5) * qJD(4) + qJD(1) * t234 + t197 * t71;
t283 = qJD(4) * t194;
t58 = -t194 * qJDD(2) - t144 * t283 + t198 * t378;
t375 = t194 * t57 + t198 * t58;
t374 = -m(5) - t383;
t373 = -mrSges(5,3) - mrSges(3,1) + mrSges(4,2);
t238 = mrSges(5,1) * t194 + mrSges(5,2) * t198;
t372 = -m(6) * t241 - t194 * t394 + t393 * t198 + mrSges(3,2) - mrSges(4,3) - t238;
t21 = qJD(6) * t244 + t192 * t75 + t196 * t74;
t350 = t21 / 0.2e1;
t22 = -qJD(6) * t80 - t192 * t74 + t196 * t75;
t349 = t22 / 0.2e1;
t342 = t74 / 0.2e1;
t341 = t75 / 0.2e1;
t138 = qJDD(6) + t146;
t333 = t138 / 0.2e1;
t332 = t146 / 0.2e1;
t164 = t200 * t193;
t165 = t200 * t197;
t93 = t164 * t192 - t165 * t196;
t370 = -qJD(6) * t93 - t192 * t390 + t196 * t389;
t92 = t164 * t196 + t165 * t192;
t369 = qJD(6) * t92 + t192 * t389 + t196 * t390;
t368 = -mrSges(6,1) - t351;
t366 = pkin(5) * t280 - t285 - (-qJD(1) * t322 + t144) * t194;
t123 = t222 * t198;
t152 = t192 * t197 + t193 * t196;
t134 = t152 * qJD(1);
t85 = t357 * t152;
t365 = -qJD(4) * t123 - t194 * t85 - t134;
t121 = t152 * t198;
t364 = t222 * qJD(1) - qJD(4) * t121 + t194 * t379;
t363 = t222 * t194;
t148 = t167 * t288;
t82 = t193 * t139 + t148;
t189 = qJ(5) + qJ(6);
t185 = sin(t189);
t186 = cos(t189);
t237 = -mrSges(6,1) * t197 + mrSges(6,2) * t193;
t362 = m(6) * pkin(4) + mrSges(7,1) * t186 - mrSges(7,2) * t185 - t237 + t394;
t361 = -m(6) * pkin(8) - t393;
t52 = mrSges(6,1) * t146 - mrSges(6,3) * t74;
t53 = -mrSges(6,2) * t146 + mrSges(6,3) * t75;
t359 = -t193 * t52 + t197 * t53;
t126 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t158;
t312 = mrSges(6,3) * t150;
t100 = mrSges(6,1) * t168 - t312;
t313 = mrSges(6,3) * t149;
t99 = -mrSges(6,2) * t168 + t313;
t223 = -t197 * t100 - t193 * t99;
t37 = -mrSges(7,1) * t244 + mrSges(7,2) * t80;
t67 = -pkin(5) * t149 + t90;
t356 = m(5) * t57 + m(6) * t377 + t223 * qJD(5) + t126 + t359 + (-m(5) * t105 + m(7) * t67 + t37 - t387) * qJD(4);
t265 = mrSges(5,3) * t184;
t161 = -qJD(4) * mrSges(5,2) - t265;
t55 = -qJDD(4) * pkin(4) - t58;
t29 = -pkin(5) * t75 + t55;
t8 = -mrSges(7,1) * t22 + mrSges(7,2) * t21;
t355 = m(5) * (qJD(4) * t106 + t58) - m(6) * (-t282 * t51 + t284 * t50 + t55) - qJD(4) * (t100 * t193 - t197 * t99 - t161) - t8 - m(7) * t29 + t380;
t353 = Ifges(7,4) * t350 + Ifges(7,2) * t349 + Ifges(7,6) * t333;
t352 = Ifges(7,1) * t350 + Ifges(7,4) * t349 + Ifges(7,5) * t333;
t348 = Ifges(6,1) * t342 + Ifges(6,4) * t341 + Ifges(6,5) * t332;
t325 = Ifges(7,4) * t80;
t33 = Ifges(7,2) * t244 + Ifges(7,6) * t166 + t325;
t347 = -t33 / 0.2e1;
t346 = t33 / 0.2e1;
t76 = Ifges(7,4) * t244;
t34 = Ifges(7,1) * t80 + Ifges(7,5) * t166 + t76;
t345 = -t34 / 0.2e1;
t344 = t34 / 0.2e1;
t339 = t244 / 0.2e1;
t337 = t80 / 0.2e1;
t336 = -m(4) - m(5);
t330 = t150 / 0.2e1;
t328 = t166 / 0.2e1;
t326 = mrSges(7,3) * t13;
t195 = sin(qJ(1));
t324 = pkin(1) * t195;
t323 = pkin(5) * t150;
t319 = g(3) * t198;
t318 = t14 * mrSges(7,3);
t199 = cos(qJ(1));
t187 = t199 * pkin(1);
t292 = t185 * t194;
t101 = -t178 * t292 + t179 * t186;
t291 = t186 * t194;
t102 = t178 * t291 + t179 * t185;
t315 = t101 * mrSges(7,1) - t102 * mrSges(7,2);
t103 = t178 * t186 + t179 * t292;
t104 = -t178 * t185 + t179 * t291;
t314 = t103 * mrSges(7,1) + t104 * mrSges(7,2);
t309 = Ifges(6,4) * t150;
t308 = Ifges(6,4) * t193;
t307 = Ifges(6,4) * t197;
t293 = t179 * t193;
t290 = t193 * t194;
t289 = t193 * t198;
t287 = t197 * t198;
t278 = qJD(5) * t198;
t277 = qJDD(1) * mrSges(4,2);
t270 = Ifges(7,5) * t21 + Ifges(7,6) * t22 + Ifges(7,3) * t138;
t269 = Ifges(6,5) * t74 + Ifges(6,6) * t75 + Ifges(6,3) * t146;
t267 = -t336 + t383;
t261 = t167 * t290;
t260 = t179 * pkin(2) + t178 * qJ(3) + t187;
t257 = t167 * t283;
t256 = t167 * t281;
t255 = t193 * t283;
t247 = -t167 * t193 + pkin(5);
t246 = -t275 / 0.2e1;
t239 = mrSges(5,1) * t198 - mrSges(5,2) * t194;
t236 = mrSges(6,1) * t193 + mrSges(6,2) * t197;
t235 = -mrSges(7,1) * t185 - mrSges(7,2) * t186;
t233 = Ifges(6,1) * t197 - t308;
t232 = Ifges(6,1) * t193 + t307;
t230 = -Ifges(6,2) * t193 + t307;
t229 = Ifges(6,2) * t197 + t308;
t228 = -Ifges(5,5) * t194 - Ifges(5,6) * t198;
t227 = Ifges(6,5) * t197 - Ifges(6,6) * t193;
t226 = Ifges(6,5) * t193 + Ifges(6,6) * t197;
t118 = t197 * t139;
t63 = -pkin(9) * t287 + t194 * t247 + t118;
t68 = -pkin(9) * t289 + t82;
t27 = -t192 * t68 + t196 * t63;
t28 = t192 * t63 + t196 * t68;
t224 = t193 * t51 + t197 * t50;
t160 = t171 * qJD(1);
t221 = qJD(3) * t160 + t145 * t171;
t220 = t270 + t386;
t112 = t178 * t197 + t179 * t290;
t110 = -t178 * t290 + t179 * t197;
t218 = t160 * t239;
t217 = t194 * (-Ifges(5,2) * t198 - t311);
t216 = t198 * (-Ifges(5,1) * t194 - t310);
t212 = t193 * t278 + t194 * t282;
t211 = -t197 * t278 + t255;
t147 = qJD(4) * t242 + qJD(3);
t40 = -qJD(5) * t261 + t139 * t279 + t193 * t147 + t197 * t256;
t207 = Ifges(6,5) * t198 - t194 * t233;
t206 = Ifges(6,6) * t198 - t194 * t230;
t205 = Ifges(6,3) * t198 - t194 * t227;
t159 = qJDD(1) * t174 + qJDD(3);
t153 = t238 * qJD(1);
t141 = t236 * t198;
t128 = t197 * t147;
t124 = (-t167 + t322) * t198;
t120 = t152 * t194;
t113 = -t178 * t193 + t179 * t288;
t111 = t178 * t288 + t293;
t109 = qJD(1) * t363;
t108 = t194 * t134;
t87 = -pkin(5) * t211 + t257;
t81 = t118 - t261;
t70 = t149 * Ifges(6,2) + t168 * Ifges(6,6) + t309;
t60 = mrSges(7,1) * t166 - mrSges(7,3) * t80;
t59 = -mrSges(7,2) * t166 + mrSges(7,3) * t244;
t45 = t152 * t283 + t198 * t379;
t43 = qJD(4) * t363 - t198 * t85;
t41 = -qJD(5) * t82 - t193 * t256 + t128;
t31 = pkin(9) * t211 + t40;
t30 = t128 + (-t148 + (pkin(9) * t198 - t139) * t193) * qJD(5) + (t198 * t247 + t271) * qJD(4);
t25 = t74 * Ifges(6,4) + t75 * Ifges(6,2) + t146 * Ifges(6,6);
t18 = -mrSges(7,2) * t138 + mrSges(7,3) * t22;
t17 = mrSges(7,1) * t138 - mrSges(7,3) * t21;
t16 = t196 * t38 - t303;
t15 = -t192 * t38 - t300;
t5 = -qJD(6) * t28 - t192 * t31 + t196 * t30;
t4 = qJD(6) * t27 + t192 * t30 + t196 * t31;
t1 = [(t270 / 0.2e1 + t269 / 0.2e1 - Ifges(5,4) * t157 / 0.2e1 - Ifges(5,2) * t158 / 0.2e1 + Ifges(7,3) * t333 + Ifges(7,6) * t349 + Ifges(7,5) * t350 - Ifges(5,6) * qJDD(4) + Ifges(6,3) * t332 + Ifges(6,6) * t341 + Ifges(6,5) * t342 + t385 + t386) * t194 + (-t121 * t2 + t123 * t3 - t13 * t43 + t14 * t45) * mrSges(7,3) + t29 * (mrSges(7,1) * t121 - mrSges(7,2) * t123) + (Ifges(4,1) + Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t191 - 0.2e1 * mrSges(3,2) * t190 + m(3) * (t190 ^ 2 + t191 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) - (t193 * t71 + t197 * t70) * t278 / 0.2e1 + t231 * t391 + t234 * t392 + (Ifges(5,1) * t392 + Ifges(5,4) * t391 + Ifges(5,5) * qJDD(4) + t227 * t332 + t230 * t341 + t233 * t342) * t198 + (-Ifges(7,4) * t123 - Ifges(7,2) * t121) * t349 + (t238 + 0.2e1 * mrSges(4,3)) * t145 + (t50 * mrSges(6,1) - t51 * mrSges(6,2) - t106 * mrSges(5,3) + Ifges(7,5) * t337 + Ifges(7,6) * t339 + Ifges(7,3) * t328 - t384) * t281 + (m(3) * t324 + mrSges(2,1) * t195 - t113 * mrSges(6,1) - t104 * mrSges(7,1) + mrSges(2,2) * t199 + t112 * mrSges(6,2) + t103 * mrSges(7,2) + (-m(4) + t374) * (t179 * qJ(3) - t324) + (m(4) * pkin(2) + t382 + t374 * (-pkin(2) - pkin(7)) - t373) * t178 + t372 * t179) * g(1) + t367 * t257 - t376 * t283 / 0.2e1 + m(5) * t221 + t161 * t256 + (t380 * t198 + m(5) * ((-t105 * t194 + t106 * t198) * qJD(4) + t375) + m(6) * (-t198 * t55 + t283 * t90) + t194 * t126) * t167 + (Ifges(7,1) * t43 + Ifges(7,4) * t45) * t337 + t217 * t246 + t174 * t277 + qJD(4) * t218 + t90 * (-mrSges(6,1) * t211 - mrSges(6,2) * t212) + m(7) * (t124 * t29 + t13 * t5 + t14 * t4 + t2 * t28 + t27 * t3 + t67 * t87) + (-Ifges(7,5) * t123 - Ifges(7,6) * t121) * t333 + m(4) * (t159 * t174 + t221) + (-Ifges(7,1) * t123 - Ifges(7,4) * t121) * t350 + (-t293 * t351 - m(3) * t187 - m(4) * t260 - mrSges(2,1) * t199 - t111 * mrSges(6,1) - t102 * mrSges(7,1) + mrSges(2,2) * t195 - t110 * mrSges(6,2) - t101 * mrSges(7,2) + t374 * (t179 * pkin(7) + t260) + t373 * t179 + t372 * t178) * g(2) + (Ifges(7,5) * t43 + Ifges(7,6) * t45) * t328 + m(6) * (t11 * t82 + t12 * t81 + t40 * t51 + t41 * t50) + qJD(4) ^ 2 * t228 / 0.2e1 + t70 * t255 / 0.2e1 + (-t11 * t289 - t12 * t287 + t211 * t51 + t212 * t50) * mrSges(6,3) + (qJD(4) * t207 - t232 * t278) * t330 + (Ifges(7,4) * t43 + Ifges(7,2) * t45) * t339 + t216 * t275 / 0.2e1 + t43 * t344 + t45 * t346 + t287 * t348 - t123 * t352 - t121 * t353 + t149 * (qJD(4) * t206 - t229 * t278) / 0.2e1 + t168 * (qJD(4) * t205 - t226 * t278) / 0.2e1 + (t105 * t283 - t375) * mrSges(5,3) - t25 * t289 / 0.2e1 + t27 * t17 + t28 * t18 + t4 * t59 + t5 * t60 + t67 * (-mrSges(7,1) * t45 + mrSges(7,2) * t43) + t81 * t52 + t82 * t53 + t87 * t37 + t40 * t99 + t41 * t100 + t124 * t8 + t55 * t141 + qJD(3) * t153 + t159 * mrSges(4,2) + t171 * (-mrSges(5,1) * t158 + mrSges(5,2) * t157); m(7) * (-t121 * t3 - t123 * t2 + t13 * t45 + t14 * t43) + t43 * t59 + t45 * t60 - t121 * t17 - t123 * t18 + (m(4) + m(3)) * qJDD(2) + (-m(3) - t267) * g(3) - t355 * t194 + t356 * t198; t277 - t120 * t17 - t363 * t18 + t364 * t60 + t365 * t59 + m(4) * t159 - t388 * t267 + (-m(6) * t224 - mrSges(4,3) * qJD(1) + t160 * t336 - t153 + t223) * qJD(1) + (-t120 * t3 + t13 * t364 + t14 * t365 - t2 * t363) * m(7) + t355 * t198 + t356 * t194; (-t265 - t161) * t105 + (-t51 * (-mrSges(6,2) * t198 + mrSges(6,3) * t290) - t50 * (mrSges(6,1) * t198 + mrSges(6,3) * t288) - t218) * qJD(1) + (Ifges(7,4) * t109 + Ifges(7,2) * t108) * t340 + (Ifges(7,1) * t109 + Ifges(7,4) * t108) * t338 + t366 * t37 + (-t100 * t279 - t99 * t280 + m(6) * (-qJD(5) * t224 + t225) + t359) * pkin(8) + (-t108 * t14 + t109 * t13 - t152 * t3 - t2 * t222) * mrSges(7,3) + (Ifges(7,5) * t152 - Ifges(7,6) * t222) * t333 + (Ifges(7,4) * t152 - Ifges(7,2) * t222) * t349 + (Ifges(7,1) * t152 - Ifges(7,4) * t222) * t350 + t29 * (mrSges(7,1) * t222 + mrSges(7,2) * t152) - t222 * t353 + (t194 * t362 + t198 * t361 + t238) * g(3) + t168 * t90 * t236 + t376 * t184 / 0.2e1 + t377 * mrSges(6,3) + ((-t109 - t379) * mrSges(7,2) + (t108 + t85) * mrSges(7,1)) * t67 - (Ifges(7,1) * t337 + Ifges(7,4) * t339 + Ifges(7,5) * t328 - t326 + t344) * t379 - (Ifges(7,4) * t337 + Ifges(7,2) * t339 + Ifges(7,6) * t328 + t318 + t346) * t85 + t228 * t246 + t388 * (t194 * t361 - t198 * t362 - t239) + (Ifges(7,5) * t109 + Ifges(7,6) * t108) * t329 + (-pkin(4) * t55 - t50 * t64 - t51 * t65) * m(6) + (t149 * t230 + t150 * t233 + t168 * t227) * qJD(5) / 0.2e1 - (t149 * t206 + t150 * t207 + t168 * t205) * qJD(1) / 0.2e1 - (t280 + t258) * t70 / 0.2e1 + (t217 / 0.2e1 - t216 / 0.2e1) * qJD(1) ^ 2 + t55 * t237 + (t384 + t398) * t286 + Ifges(5,3) * qJDD(4) + t197 * t25 / 0.2e1 + t226 * t332 + t369 * t59 + t370 * t60 + (t13 * t370 + t14 * t369 - t177 * t29 + t2 * t93 + t3 * t92 + t366 * t67) * m(7) + t229 * t341 + t232 * t342 + t109 * t345 + t108 * t347 + t193 * t348 + t152 * t352 + t71 * t279 / 0.2e1 - pkin(4) * t36 - t57 * mrSges(5,2) + t58 * mrSges(5,1) + t92 * t17 + t93 * t18 - t65 * t99 - t64 * t100 + (t264 + t387) * t106 + Ifges(5,5) * t157 + Ifges(5,6) * t158 - t177 * t8; (-mrSges(6,2) * t113 + t112 * t368 - t314) * g(2) + (mrSges(6,2) * t111 + t110 * t368 - t315) * g(1) - (t235 - t382) * t319 - (mrSges(7,1) * t67 + Ifges(7,4) * t338 + Ifges(7,2) * t340 + Ifges(7,6) * t329 - t318 + t347) * t80 + (-mrSges(7,2) * t67 + Ifges(7,1) * t338 + Ifges(7,4) * t340 + Ifges(7,5) * t329 + t326 + t345) * t244 + t385 + (t313 - t99) * t50 + t220 + (t312 + t100) * t51 + t269 + (-Ifges(6,2) * t150 + t142 + t71) * t397 + t70 * t330 + (t192 * t2 + t196 * t3 + (-t13 * t192 + t14 * t196) * qJD(6)) * t351 + (Ifges(6,1) * t149 - t309) * t396 - t16 * t59 - t15 * t60 - t37 * t323 - m(7) * (t13 * t15 + t14 * t16 + t323 * t67) + g(3) * t141 - t90 * (mrSges(6,1) * t150 + mrSges(6,2) * t149) + (Ifges(6,5) * t149 - Ifges(6,6) * t150) * t395 + ((-t192 * t60 + t196 * t59) * qJD(6) + t196 * t17 + t192 * t18) * pkin(5); -t67 * (mrSges(7,1) * t80 + mrSges(7,2) * t244) + (Ifges(7,1) * t244 - t325) * t338 + t33 * t337 + (Ifges(7,5) * t244 - Ifges(7,6) * t80) * t329 - t13 * t59 + t14 * t60 - g(1) * t315 - g(2) * t314 - t235 * t319 + (t13 * t244 + t14 * t80) * mrSges(7,3) + t220 + (-Ifges(7,2) * t80 + t34 + t76) * t340;];
tau  = t1;
