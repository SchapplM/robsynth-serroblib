% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:33
% EndTime: 2019-03-09 02:30:51
% DurationCPUTime: 13.55s
% Computational Cost: add. (5848->615), mult. (10733->834), div. (0->0), fcn. (6288->10), ass. (0->285)
t406 = -Ifges(5,2) / 0.2e1;
t201 = -pkin(9) - pkin(8);
t405 = m(7) * t201 - mrSges(6,3) - mrSges(7,3);
t195 = sin(qJ(4));
t178 = t195 * qJD(1);
t273 = qJD(5) + qJD(6);
t166 = t178 + t273;
t338 = -t166 / 0.2e1;
t194 = sin(qJ(5));
t198 = cos(qJ(5));
t281 = qJD(4) * t198;
t199 = cos(qJ(4));
t285 = qJD(1) * t199;
t141 = -t194 * t285 + t281;
t142 = qJD(4) * t194 + t198 * t285;
t193 = sin(qJ(6));
t197 = cos(qJ(6));
t72 = t141 * t193 + t142 * t197;
t345 = -t72 / 0.2e1;
t249 = t141 * t197 - t142 * t193;
t347 = -t249 / 0.2e1;
t404 = Ifges(7,5) * t345 + Ifges(7,6) * t347 + Ifges(7,3) * t338;
t403 = -t141 / 0.2e1;
t402 = -t142 / 0.2e1;
t171 = t178 + qJD(5);
t401 = -t171 / 0.2e1;
t318 = Ifges(5,4) * t199;
t400 = t195 * t406 + t318 / 0.2e1;
t174 = pkin(5) * t198 + pkin(4);
t399 = m(7) * t174;
t264 = t194 * t178;
t265 = qJD(5) * t201;
t248 = pkin(4) * t199 + pkin(8) * t195;
t148 = t248 * qJD(1);
t175 = qJ(2) * qJD(1) + qJD(3);
t164 = -pkin(7) * qJD(1) + t175;
t288 = t198 * t199;
t77 = t148 * t194 + t164 * t288;
t398 = -pkin(9) * t264 + t194 * t265 - t77;
t292 = t195 * t198;
t223 = pkin(5) * t199 + pkin(9) * t292;
t295 = t194 * t199;
t76 = t148 * t198 - t164 * t295;
t397 = -qJD(1) * t223 + t198 * t265 - t76;
t190 = qJ(5) + qJ(6);
t179 = sin(t190);
t180 = cos(t190);
t243 = -mrSges(6,1) * t198 + mrSges(6,2) * t194;
t396 = m(6) * pkin(4) + mrSges(7,1) * t180 - mrSges(7,2) * t179 - t243 + t399;
t395 = -m(6) * pkin(8) + t405;
t279 = qJD(5) * t194;
t394 = t279 + t264;
t376 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t141 - mrSges(6,2) * t142 - mrSges(5,3) * t285;
t144 = t193 * t198 + t194 * t197;
t393 = t273 * t144;
t191 = qJ(2) - pkin(7);
t300 = t191 * t195;
t392 = qJD(4) * t248 - qJD(5) * t300 + qJD(3);
t280 = qJD(4) * t199;
t370 = qJD(2) * t195 + t191 * t280;
t196 = sin(qJ(1));
t200 = cos(qJ(1));
t391 = g(1) * t200 + g(2) * t196;
t192 = pkin(1) + qJ(3);
t247 = pkin(4) * t195 - pkin(8) * t199;
t154 = t247 + t192;
t100 = qJD(1) * t154 - qJD(2);
t153 = t195 * t164;
t121 = qJD(4) * pkin(8) + t153;
t278 = qJD(5) * t198;
t275 = qJD(1) * qJD(3);
t146 = qJDD(1) * t192 - qJDD(2) + t275;
t274 = qJD(1) * qJD(4);
t151 = qJDD(1) * t199 - t195 * t274;
t152 = -qJDD(1) * t195 - t199 * t274;
t57 = -pkin(4) * t152 - pkin(8) * t151 + t146;
t186 = qJD(1) * qJD(2);
t372 = qJ(2) * qJDD(1) + t186;
t160 = qJDD(3) + t372;
t145 = -pkin(7) * qJDD(1) + t160;
t83 = t145 * t195 + t164 * t280;
t74 = qJDD(4) * pkin(8) + t83;
t15 = t100 * t278 - t121 * t279 + t194 * t57 + t198 * t74;
t66 = -qJD(5) * t142 + qJDD(4) * t198 - t151 * t194;
t10 = pkin(9) * t66 + t15;
t54 = t100 * t194 + t121 * t198;
t47 = pkin(9) * t141 + t54;
t312 = t193 * t47;
t53 = t100 * t198 - t121 * t194;
t46 = -pkin(9) * t142 + t53;
t37 = pkin(5) * t171 + t46;
t13 = t197 * t37 - t312;
t137 = qJDD(5) - t152;
t16 = -qJD(5) * t54 - t194 * t74 + t198 * t57;
t65 = qJD(5) * t141 + qJDD(4) * t194 + t151 * t198;
t9 = pkin(5) * t137 - pkin(9) * t65 + t16;
t2 = qJD(6) * t13 + t10 * t197 + t193 * t9;
t309 = t197 * t47;
t14 = t193 * t37 + t309;
t3 = -qJD(6) * t14 - t10 * t193 + t197 * t9;
t390 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t389 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t388 = t14 * mrSges(7,2) + Ifges(5,6) * qJD(4) / 0.2e1 + qJD(1) * t400 + Ifges(6,5) * t402 + Ifges(6,6) * t403 + Ifges(6,3) * t401 - t13 * mrSges(7,1) + t404;
t130 = qJDD(6) + t137;
t21 = qJD(6) * t249 + t193 * t66 + t197 * t65;
t22 = -qJD(6) * t72 - t193 * t65 + t197 * t66;
t387 = -t21 * Ifges(7,4) / 0.2e1 - t22 * Ifges(7,2) / 0.2e1 - t130 * Ifges(7,6) / 0.2e1;
t357 = m(7) * pkin(5);
t356 = t21 / 0.2e1;
t355 = t22 / 0.2e1;
t349 = t65 / 0.2e1;
t348 = t66 / 0.2e1;
t343 = -m(4) - m(5);
t386 = -m(6) - m(7);
t342 = t130 / 0.2e1;
t341 = t137 / 0.2e1;
t335 = t199 / 0.2e1;
t385 = mrSges(3,2) - mrSges(4,3);
t384 = -mrSges(4,2) - mrSges(3,3);
t161 = t201 * t194;
t162 = t201 * t198;
t86 = t161 * t193 - t162 * t197;
t383 = -qJD(6) * t86 - t193 * t398 + t197 * t397;
t85 = t161 * t197 + t162 * t193;
t382 = qJD(6) * t85 + t193 * t397 + t197 * t398;
t276 = qJD(6) * t197;
t298 = t193 * t194;
t78 = -t197 * t278 - t198 * t276 + t273 * t298;
t289 = t197 * t198;
t98 = t178 * t289 - t193 * t264;
t380 = t78 - t98;
t120 = t144 * qJD(1);
t97 = t195 * t120;
t379 = t393 + t97;
t378 = mrSges(6,1) + t357;
t32 = -mrSges(6,1) * t66 + mrSges(6,2) * t65;
t377 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t151 - t32;
t375 = pkin(5) * t394 - t153;
t143 = -t289 + t298;
t108 = t143 * t199;
t374 = -qJD(4) * t108 - t195 * t393 - t120;
t106 = t144 * t199;
t203 = t273 * t143;
t373 = qJD(1) * t143 - qJD(4) * t106 + t195 * t203;
t282 = qJD(4) * t195;
t371 = qJD(2) * t199 - t191 * t282;
t319 = Ifges(5,4) * t195;
t240 = t199 * Ifges(5,1) - t319;
t134 = Ifges(6,4) * t141;
t61 = t142 * Ifges(6,1) + t171 * Ifges(6,5) + t134;
t369 = Ifges(5,5) * qJD(4) + qJD(1) * t240 + t198 * t61;
t48 = mrSges(6,1) * t137 - mrSges(6,3) * t65;
t49 = -mrSges(6,2) * t137 + mrSges(6,3) * t66;
t368 = -t194 * t48 + t198 * t49;
t82 = t145 * t199 - t164 * t282;
t367 = -t195 * t83 - t199 * t82;
t366 = t15 * t198 - t16 * t194;
t365 = -m(5) + t386;
t301 = t164 * t199;
t122 = -qJD(4) * pkin(4) - t301;
t80 = -pkin(5) * t141 + t122;
t363 = t80 * mrSges(7,1) - t14 * mrSges(7,3);
t362 = mrSges(5,3) + mrSges(2,2) + t384;
t68 = Ifges(7,4) * t249;
t30 = Ifges(7,1) * t72 + Ifges(7,5) * t166 + t68;
t361 = -t80 * mrSges(7,2) + t13 * mrSges(7,3) - t30 / 0.2e1;
t244 = mrSges(5,1) * t195 + mrSges(5,2) * t199;
t360 = t195 * t399 + t199 * t405 + mrSges(2,1) + t244 - t385;
t202 = qJD(1) ^ 2;
t359 = -m(5) / 0.2e1;
t358 = Ifges(7,1) * t356 + Ifges(7,4) * t355 + Ifges(7,5) * t342;
t354 = Ifges(6,1) * t349 + Ifges(6,4) * t348 + Ifges(6,5) * t341;
t334 = Ifges(7,4) * t72;
t29 = Ifges(7,2) * t249 + Ifges(7,6) * t166 + t334;
t353 = -t29 / 0.2e1;
t351 = t30 / 0.2e1;
t346 = t249 / 0.2e1;
t344 = t72 / 0.2e1;
t339 = t142 / 0.2e1;
t337 = t166 / 0.2e1;
t333 = pkin(5) * t142;
t331 = pkin(5) * t194;
t328 = g(3) * t199;
t293 = t195 * t196;
t101 = t179 * t293 - t180 * t200;
t102 = -t179 * t200 - t180 * t293;
t323 = -mrSges(7,1) * t101 + mrSges(7,2) * t102;
t291 = t195 * t200;
t103 = -t179 * t291 - t180 * t196;
t104 = -t179 * t196 + t180 * t291;
t322 = mrSges(7,1) * t103 - mrSges(7,2) * t104;
t321 = mrSges(6,3) * t141;
t320 = mrSges(6,3) * t142;
t317 = Ifges(6,4) * t142;
t316 = Ifges(6,4) * t194;
t315 = Ifges(6,4) * t198;
t303 = qJDD(1) * pkin(1);
t297 = t194 * t195;
t296 = t194 * t196;
t294 = t194 * t200;
t290 = t196 * t198;
t287 = t198 * t200;
t89 = t154 * t194 + t191 * t292;
t286 = pkin(1) * t200 + qJ(2) * t196;
t277 = qJD(5) * t199;
t35 = -mrSges(7,1) * t249 + mrSges(7,2) * t72;
t272 = t35 - t376;
t271 = Ifges(7,5) * t21 + Ifges(7,6) * t22 + Ifges(7,3) * t130;
t270 = Ifges(6,5) * t65 + Ifges(6,6) * t66 + Ifges(6,3) * t137;
t268 = t343 + t386;
t266 = qJ(3) * t200 + t286;
t261 = t194 * t282;
t252 = -t274 / 0.2e1;
t251 = (t195 ^ 2 + t199 ^ 2) * t164;
t246 = -t152 * mrSges(5,1) + t151 * mrSges(5,2);
t245 = mrSges(5,1) * t199 - mrSges(5,2) * t195;
t242 = mrSges(6,1) * t194 + mrSges(6,2) * t198;
t241 = -mrSges(7,1) * t179 - mrSges(7,2) * t180;
t239 = Ifges(6,1) * t198 - t316;
t238 = Ifges(6,1) * t194 + t315;
t236 = -Ifges(6,2) * t194 + t315;
t235 = Ifges(6,2) * t198 + t316;
t234 = -Ifges(5,5) * t195 - Ifges(5,6) * t199;
t233 = Ifges(6,5) * t198 - Ifges(6,6) * t194;
t232 = Ifges(6,5) * t194 + Ifges(6,6) * t198;
t133 = t198 * t154;
t67 = -pkin(9) * t288 + t133 + (-t191 * t194 + pkin(5)) * t195;
t75 = -pkin(9) * t295 + t89;
t33 = -t193 * t75 + t197 * t67;
t34 = t193 * t67 + t197 * t75;
t231 = t194 * t54 + t198 * t53;
t230 = t194 * t53 - t198 * t54;
t90 = -mrSges(6,2) * t171 + t321;
t91 = mrSges(6,1) * t171 - t320;
t229 = t194 * t91 - t198 * t90;
t228 = -t194 * t90 - t198 * t91;
t165 = qJD(1) * t192 - qJD(2);
t225 = qJD(3) * t165 + t146 * t192;
t222 = t271 + t390;
t221 = t392 * t198;
t157 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t178;
t220 = -t157 + t229;
t125 = -t194 * t291 - t290;
t123 = t194 * t293 - t287;
t73 = -qJDD(4) * pkin(4) - t82;
t219 = t165 * t245;
t218 = t195 * (-Ifges(5,2) * t199 - t319);
t217 = t199 * (-Ifges(5,1) * t195 - t318);
t212 = t194 * t277 + t195 * t281;
t211 = -t198 * t277 + t261;
t38 = t154 * t278 + t194 * t392 + t198 * t370;
t210 = Ifges(6,5) * t199 - t195 * t239;
t209 = Ifges(6,6) * t199 - t195 * t236;
t208 = Ifges(6,3) * t199 - t195 * t233;
t207 = -qJD(5) * t231 + t366;
t183 = t200 * qJ(2);
t176 = qJDD(2) - t303;
t147 = t244 * qJD(1);
t136 = (-t191 + t331) * t199;
t131 = t242 * t199;
t126 = t195 * t287 - t296;
t124 = -t195 * t290 - t294;
t110 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t152;
t107 = t143 * t195;
t105 = t144 * t195;
t88 = -t191 * t297 + t133;
t84 = -pkin(5) * t211 - t371;
t60 = Ifges(6,2) * t141 + Ifges(6,6) * t171 + t317;
t51 = mrSges(7,1) * t166 - mrSges(7,3) * t72;
t50 = -mrSges(7,2) * t166 + mrSges(7,3) * t249;
t43 = t144 * t282 + t199 * t203;
t41 = t143 * t282 - t199 * t393;
t39 = (-qJD(5) * t154 - t370) * t194 + t221;
t36 = -pkin(5) * t66 + t73;
t31 = pkin(9) * t211 + t38;
t27 = t223 * qJD(4) + ((pkin(9) * t199 - t154) * qJD(5) - t370) * t194 + t221;
t23 = Ifges(6,4) * t65 + Ifges(6,2) * t66 + Ifges(6,6) * t137;
t18 = t197 * t46 - t312;
t17 = -t193 * t46 - t309;
t12 = -mrSges(7,2) * t130 + mrSges(7,3) * t22;
t11 = mrSges(7,1) * t130 - mrSges(7,3) * t21;
t8 = -mrSges(7,1) * t22 + mrSges(7,2) * t21;
t7 = -qJD(6) * t34 - t193 * t31 + t197 * t27;
t6 = qJD(6) * t33 + t193 * t27 + t197 * t31;
t1 = [(-t303 + t176) * mrSges(3,2) + (t271 / 0.2e1 + t270 / 0.2e1 + Ifges(7,3) * t342 + Ifges(7,6) * t355 + Ifges(7,5) * t356 - Ifges(5,4) * t151 / 0.2e1 + t152 * t406 + Ifges(6,3) * t341 + Ifges(6,6) * t348 + Ifges(6,5) * t349 - Ifges(5,6) * qJDD(4) + t389 + t390) * t195 + (-t106 * t2 + t108 * t3 - t13 * t41 + t14 * t43) * mrSges(7,3) + (mrSges(4,3) * t192 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) - (t194 * t61 + t198 * t60) * t277 / 0.2e1 + t36 * (mrSges(7,1) * t106 - mrSges(7,2) * t108) + 0.2e1 * t372 * mrSges(3,3) + t376 * t371 + (-t15 * t295 - t16 * t288 + t211 * t54 + t212 * t53) * mrSges(6,3) + t146 * t244 + (-Ifges(7,1) * t108 - Ifges(7,4) * t106) * t356 + (-Ifges(7,4) * t108 - Ifges(7,2) * t106) * t355 + m(3) * (-pkin(1) * t176 + (t372 + t186) * qJ(2)) + m(6) * (-t122 * t371 + t15 * t89 + t16 * t88 + t38 * t54 + t39 * t53) + m(5) * (qJD(2) * t251 - t191 * t367 + t225) + t367 * mrSges(5,3) - t369 * t282 / 0.2e1 + t370 * t157 + (t160 + t372) * mrSges(4,2) + (mrSges(6,1) * t53 - mrSges(6,2) * t54 + Ifges(7,5) * t344 + Ifges(7,6) * t346 + Ifges(7,3) * t337 - t388) * t280 + t41 * t351 + t288 * t354 + (qJD(4) * t210 - t238 * t277) * t339 + t110 * t300 + qJD(4) * t219 + m(7) * (t13 * t7 + t136 * t36 + t14 * t6 + t2 * t34 + t3 * t33 + t80 * t84) - t23 * t295 / 0.2e1 + (Ifges(5,1) * t151 + Ifges(5,4) * t152 + 0.2e1 * Ifges(5,5) * qJDD(4)) * t335 + (Ifges(7,1) * t41 + Ifges(7,4) * t43) * t344 + (t296 * t357 - m(3) * t286 - m(4) * t266 - t126 * mrSges(6,1) - t104 * mrSges(7,1) - t125 * mrSges(6,2) - t103 * mrSges(7,2) + t365 * (-pkin(7) * t196 + t266) + t362 * t196 + (-m(6) * t247 - t360) * t200) * g(2) + (t294 * t357 - t124 * mrSges(6,1) - t102 * mrSges(7,1) - t123 * mrSges(6,2) - t101 * mrSges(7,2) + t365 * (-pkin(7) * t200 + t183) + (-m(4) - m(3)) * t183 + t362 * t200 + (m(3) * pkin(1) + m(6) * t154 + (m(7) - t343) * t192 + t360) * t196) * g(1) + t60 * t261 / 0.2e1 + t171 * (qJD(4) * t208 - t232 * t277) / 0.2e1 + t141 * (qJD(4) * t209 - t235 * t277) / 0.2e1 + t217 * t274 / 0.2e1 + (-Ifges(7,5) * t108 - Ifges(7,6) * t106) * t342 - t108 * t358 + (t275 + t146) * mrSges(4,3) + (Ifges(7,5) * t41 + Ifges(7,6) * t43) * t337 + t152 * t400 + t151 * t240 / 0.2e1 + t192 * t246 + qJD(4) ^ 2 * t234 / 0.2e1 + m(4) * (qJ(2) * t160 + qJD(2) * t175 + t225) + t33 * t11 + t34 * t12 + t43 * t29 / 0.2e1 + t6 * t50 + t7 * t51 + t80 * (-mrSges(7,1) * t43 + mrSges(7,2) * t41) + t84 * t35 + t88 * t48 + t89 * t49 + t38 * t90 + t39 * t91 + (t233 * t341 + t236 * t348 + t239 * t349 + (-m(6) * t73 + t377) * t191) * t199 + t218 * t252 + t106 * t387 + t73 * t131 + t136 * t8 + t122 * (-mrSges(6,1) * t211 - mrSges(6,2) * t212) + qJD(3) * t147 + (Ifges(7,4) * t41 + Ifges(7,2) * t43) * t346; t143 * t11 - t144 * t12 - t194 * t49 - t198 * t48 + t379 * t51 + t380 * t50 + t385 * qJDD(1) + t229 * qJD(5) + (-m(3) * qJ(2) + t384) * t202 + m(3) * t176 + m(6) * (qJD(5) * t230 - t15 * t194 - t16 * t198) + 0.2e1 * (-m(4) / 0.2e1 + t359) * t146 - t246 + (-g(1) * t196 + g(2) * t200) * (m(3) - t268) + (t13 * t379 + t14 * t380 + t143 * t3 - t144 * t2) * m(7) + (-m(4) * t175 + t220 * t195 + t272 * t199 + 0.2e1 * m(7) * t80 * t335 - m(6) * (-t122 * t199 + t292 * t54 - t297 * t53) + 0.2e1 * t251 * t359) * qJD(1); qJDD(1) * mrSges(4,2) - t202 * mrSges(4,3) - t105 * t11 - t107 * t12 + t373 * t51 + t374 * t50 + (t165 * t343 - t147 + t228) * qJD(1) + (-qJD(4) * t220 + t377 - t8) * t199 + (qJD(4) * t272 + qJD(5) * t228 + t110 + t368) * t195 + m(4) * t160 - m(5) * t367 + t391 * t268 + (-t105 * t3 - t107 * t2 + t13 * t373 + t14 * t374 - t199 * t36 + t282 * t80) * m(7) + ((-qJD(4) * t230 - t73) * t199 + (qJD(4) * t122 + t207) * t195 - t231 * qJD(1)) * m(6); (-Ifges(7,5) * t98 + Ifges(7,6) * t97) * t338 + (-t217 / 0.2e1 + t218 / 0.2e1) * t202 + (-t53 * (mrSges(6,1) * t199 + mrSges(6,3) * t292) - t54 * (-mrSges(6,2) * t199 + mrSges(6,3) * t297) - t219) * qJD(1) + (-Ifges(7,4) * t98 + Ifges(7,2) * t97) * t347 + (-Ifges(7,1) * t98 + Ifges(7,4) * t97) * t345 + (-t13 * t98 - t14 * t97 - t143 * t2 - t144 * t3) * mrSges(7,3) + (t141 * t236 + t142 * t239 + t171 * t233) * qJD(5) / 0.2e1 - (t141 * t209 + t142 * t210 + t171 * t208) * qJD(1) / 0.2e1 + (-Ifges(7,1) * t344 - Ifges(7,4) * t346 - Ifges(7,5) * t337 + t361) * t78 - t394 * t60 / 0.2e1 + (t195 * t396 + t199 * t395 + t244) * g(3) + t391 * (t195 * t395 - t199 * t396 - t245) + (t388 + t404) * t285 + (-Ifges(7,4) * t344 - Ifges(7,2) * t346 - Ifges(7,6) * t337 + t363) * t393 + t376 * t153 + t375 * t35 + (-t278 * t53 - t279 * t54 + t366) * mrSges(6,3) + (m(6) * t207 - t278 * t91 - t279 * t90 + t368) * pkin(8) + t369 * t178 / 0.2e1 - t157 * t301 + t232 * t341 + (Ifges(7,5) * t144 - Ifges(7,6) * t143) * t342 + t235 * t348 + t238 * t349 + t98 * t351 + t194 * t354 + (Ifges(7,4) * t144 - Ifges(7,2) * t143) * t355 + (Ifges(7,1) * t144 - Ifges(7,4) * t143) * t356 + t61 * t278 / 0.2e1 + t144 * t358 + t379 * t353 + (-pkin(4) * t73 - t122 * t153 - t53 * t76 - t54 * t77) * m(6) + Ifges(5,3) * qJDD(4) + t73 * t243 - pkin(4) * t32 + t82 * mrSges(5,1) - t83 * mrSges(5,2) + t85 * t11 + t86 * t12 - t77 * t90 - t76 * t91 - t80 * (-mrSges(7,1) * t97 - mrSges(7,2) * t98) + t234 * t252 + t143 * t387 + t36 * (mrSges(7,1) * t143 + mrSges(7,2) * t144) + Ifges(5,5) * t151 + Ifges(5,6) * t152 - t174 * t8 + t382 * t50 + t383 * t51 + (t13 * t383 + t14 * t382 - t174 * t36 + t2 * t86 + t3 * t85 + t375 * t80) * m(7) + t171 * t242 * t122 + t198 * t23 / 0.2e1; (-Ifges(6,2) * t142 + t134 + t61) * t403 + (Ifges(7,1) * t345 + Ifges(7,4) * t347 + Ifges(7,5) * t338 + t361) * t249 - (Ifges(7,4) * t345 + Ifges(7,2) * t347 + Ifges(7,6) * t338 + t353 + t363) * t72 + (-mrSges(6,2) * t124 + t123 * t378 - t323) * g(2) + (mrSges(6,2) * t126 - t125 * t378 - t322) * g(1) + (t193 * t2 + t197 * t3 + (-t13 * t193 + t14 * t197) * qJD(6)) * t357 + t60 * t339 + (t321 - t90) * t53 + (Ifges(6,1) * t141 - t317) * t402 + (t320 + t91) * t54 - (-m(7) * t331 + t241) * t328 + t270 + t389 - t35 * t333 - m(7) * (t13 * t17 + t14 * t18 + t333 * t80) + t222 + ((-qJD(6) * t51 + t12) * t193 + t11 * t197 + t276 * t50) * pkin(5) - t18 * t50 - t17 * t51 + g(3) * t131 - t122 * (mrSges(6,1) * t142 + mrSges(6,2) * t141) + (Ifges(6,5) * t141 - Ifges(6,6) * t142) * t401; -t80 * (mrSges(7,1) * t72 + mrSges(7,2) * t249) + (Ifges(7,1) * t249 - t334) * t345 + t29 * t344 + (Ifges(7,5) * t249 - Ifges(7,6) * t72) * t338 - t13 * t50 + t14 * t51 - g(1) * t322 - g(2) * t323 - t241 * t328 + (t13 * t249 + t14 * t72) * mrSges(7,3) + t222 + (-Ifges(7,2) * t72 + t30 + t68) * t347;];
tau  = t1;
