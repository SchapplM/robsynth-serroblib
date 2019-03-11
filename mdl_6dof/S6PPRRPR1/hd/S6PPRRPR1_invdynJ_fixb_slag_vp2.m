% Calculate vector of inverse dynamics joint torques for
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:34
% EndTime: 2019-03-08 18:45:54
% DurationCPUTime: 14.24s
% Computational Cost: add. (7330->601), mult. (18571->866), div. (0->0), fcn. (16420->18), ass. (0->294)
t337 = cos(pkin(6));
t211 = qJD(1) * t337 + qJD(2);
t223 = sin(pkin(12));
t225 = sin(pkin(6));
t231 = sin(qJ(3));
t227 = cos(pkin(12));
t336 = cos(pkin(7));
t289 = t227 * t336;
t280 = t231 * t289;
t353 = cos(qJ(3));
t247 = (t223 * t353 + t280) * t225;
t224 = sin(pkin(7));
t324 = t224 * t231;
t111 = qJD(1) * t247 + t211 * t324;
t230 = sin(qJ(4));
t233 = cos(qJ(4));
t268 = pkin(4) * t230 - qJ(5) * t233;
t314 = qJD(5) * t230;
t419 = qJD(4) * t268 - t111 - t314;
t222 = sin(pkin(13));
t226 = cos(pkin(13));
t264 = t353 * t289;
t259 = t225 * t264;
t250 = qJD(1) * t259;
t326 = t223 * t225;
t302 = qJD(1) * t326;
t303 = t224 * t353;
t239 = -t211 * t303 + t231 * t302 - t250;
t316 = qJD(4) * t230;
t311 = pkin(9) * t316;
t327 = t222 * t233;
t398 = t222 * t311 + t419 * t226 - t239 * t327;
t322 = t226 * t233;
t418 = t419 * t222 + t239 * t322;
t317 = qJD(3) * t233;
t213 = qJD(6) - t317;
t356 = -t213 / 0.2e1;
t319 = qJD(3) * t230;
t186 = qJD(4) * t226 - t222 * t319;
t187 = qJD(4) * t222 + t226 * t319;
t229 = sin(qJ(6));
t232 = cos(qJ(6));
t123 = t186 * t229 + t187 * t232;
t362 = -t123 / 0.2e1;
t286 = t232 * t186 - t187 * t229;
t364 = -t286 / 0.2e1;
t417 = Ifges(7,5) * t362 + Ifges(7,6) * t364 + Ifges(7,3) * t356;
t404 = qJD(4) / 0.2e1;
t263 = pkin(5) * t230 - pkin(10) * t322;
t416 = qJD(4) * t263 + t398;
t323 = t226 * t230;
t415 = (-pkin(9) * t323 - pkin(10) * t327) * qJD(4) + t418;
t267 = t222 * t229 - t226 * t232;
t257 = t267 * t233;
t154 = qJD(3) * t257;
t168 = t267 * qJD(6);
t414 = -t168 + t154;
t190 = t222 * t232 + t226 * t229;
t258 = t190 * t233;
t153 = qJD(3) * t258;
t169 = t190 * qJD(6);
t413 = -t169 + t153;
t215 = pkin(5) * t226 + pkin(4);
t221 = pkin(13) + qJ(6);
t217 = sin(t221);
t218 = cos(t221);
t275 = -mrSges(6,1) * t226 + mrSges(6,2) * t222;
t412 = m(6) * pkin(4) + m(7) * t215 + mrSges(7,1) * t218 - mrSges(7,2) * t217 - t275;
t346 = pkin(10) + qJ(5);
t411 = m(6) * qJ(5) + m(7) * t346 + mrSges(6,3) + mrSges(7,3);
t341 = Ifges(7,4) * t123;
t66 = Ifges(7,2) * t286 + Ifges(7,6) * t213 + t341;
t367 = t66 / 0.2e1;
t117 = Ifges(7,4) * t286;
t67 = Ifges(7,1) * t123 + Ifges(7,5) * t213 + t117;
t366 = t67 / 0.2e1;
t276 = -mrSges(5,1) * t233 + mrSges(5,2) * t230;
t385 = m(5) + m(7) + m(6);
t410 = pkin(3) * t385 + t230 * t411 + t233 * t412 + mrSges(4,1) - t276;
t108 = qJD(3) * pkin(9) + t111;
t306 = t225 * t227 * t224;
t148 = -qJD(1) * t306 + t211 * t336;
t330 = t148 * t230;
t77 = t108 * t233 + t330;
t74 = qJD(4) * qJ(5) + t77;
t203 = -pkin(4) * t233 - qJ(5) * t230 - pkin(3);
t94 = qJD(3) * t203 + t239;
t34 = -t222 * t74 + t226 * t94;
t30 = -pkin(5) * t317 - pkin(10) * t187 + t34;
t35 = t222 * t94 + t226 * t74;
t31 = pkin(10) * t186 + t35;
t10 = t229 * t30 + t232 * t31;
t345 = Ifges(5,4) * t230;
t272 = t233 * Ifges(5,2) + t345;
t9 = -t229 * t31 + t232 * t30;
t409 = -t9 * mrSges(7,1) + Ifges(5,6) * t404 + qJD(3) * t272 / 0.2e1 + t10 * mrSges(7,2) - t187 * Ifges(6,5) / 0.2e1 - t186 * Ifges(6,6) / 0.2e1 + Ifges(6,3) * t317 / 0.2e1 + t417;
t313 = qJD(3) * qJD(4);
t198 = -t233 * qJDD(3) + t230 * t313;
t357 = t198 / 0.2e1;
t199 = qJDD(3) * t230 + t233 * t313;
t156 = qJDD(4) * t222 + t199 * t226;
t359 = t156 / 0.2e1;
t408 = Ifges(6,1) * t359 + Ifges(6,5) * t357;
t407 = m(5) * pkin(9);
t155 = qJDD(4) * t226 - t199 * t222;
t61 = qJD(6) * t286 + t155 * t229 + t156 * t232;
t369 = t61 / 0.2e1;
t62 = -qJD(6) * t123 + t155 * t232 - t156 * t229;
t368 = t62 / 0.2e1;
t360 = t155 / 0.2e1;
t188 = qJDD(6) + t198;
t358 = t188 / 0.2e1;
t406 = -t198 / 0.2e1;
t405 = t199 / 0.2e1;
t183 = t226 * t203;
t126 = -pkin(10) * t323 + t183 + (-pkin(9) * t222 - pkin(5)) * t233;
t150 = pkin(9) * t322 + t222 * t203;
t328 = t222 * t230;
t137 = -pkin(10) * t328 + t150;
t79 = t126 * t229 + t137 * t232;
t403 = -qJD(6) * t79 - t229 * t415 + t232 * t416;
t78 = t126 * t232 - t137 * t229;
t402 = qJD(6) * t78 + t229 * t416 + t232 * t415;
t204 = t346 * t222;
t205 = t346 * t226;
t138 = -t204 * t232 - t205 * t229;
t192 = t268 * qJD(3);
t101 = t230 * t108;
t76 = t148 * t233 - t101;
t51 = t226 * t192 - t222 * t76;
t39 = qJD(3) * t263 + t51;
t301 = t222 * t317;
t52 = t222 * t192 + t226 * t76;
t42 = -pkin(10) * t301 + t52;
t401 = -qJD(5) * t267 + qJD(6) * t138 - t229 * t39 - t232 * t42;
t139 = -t204 * t229 + t205 * t232;
t400 = -qJD(5) * t190 - qJD(6) * t139 + t229 * t42 - t232 * t39;
t397 = -t226 * t311 + t418;
t95 = -t155 * mrSges(6,1) + t156 * mrSges(6,2);
t396 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t199 + t95;
t308 = mrSges(5,3) * t319;
t395 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t186 + mrSges(6,2) * t187 + t308;
t335 = cos(pkin(11));
t278 = t337 * t335;
t334 = sin(pkin(11));
t245 = t223 * t334 - t227 * t278;
t291 = t225 * t335;
t392 = t224 * t291 + t245 * t336;
t277 = t337 * t334;
t246 = t223 * t335 + t227 * t277;
t290 = t225 * t334;
t391 = -t224 * t290 + t246 * t336;
t342 = Ifges(6,4) * t226;
t271 = -Ifges(6,2) * t222 + t342;
t343 = Ifges(6,4) * t222;
t273 = Ifges(6,1) * t226 - t343;
t390 = t186 * (Ifges(6,6) * t230 + t233 * t271) + t187 * (Ifges(6,5) * t230 + t233 * t273);
t140 = mrSges(5,1) * t198 + mrSges(5,2) * t199;
t389 = mrSges(4,1) * qJDD(3) - t140;
t191 = t276 * qJD(3);
t388 = mrSges(4,1) * qJD(3) - t191;
t210 = qJDD(1) * t337 + qJDD(2);
t146 = -qJDD(1) * t306 + t210 * t336;
t315 = qJD(4) * t233;
t266 = t225 * t280;
t282 = qJD(3) * t302;
t283 = qJD(3) * t303;
t297 = qJDD(1) * t326;
t68 = qJD(3) * t250 + qJDD(1) * t266 + t210 * t324 + t211 * t283 - t231 * t282 + t353 * t297;
t63 = qJDD(3) * pkin(9) + t68;
t309 = t230 * t146 + t148 * t315 + t233 * t63;
t21 = -t108 * t316 + t309;
t22 = -t108 * t315 + t146 * t233 - t148 * t316 - t230 * t63;
t387 = t21 * t233 - t22 * t230;
t17 = qJDD(4) * qJ(5) + (qJD(5) - t101) * qJD(4) + t309;
t318 = qJD(3) * t231;
t300 = t224 * t318;
t69 = -qJD(3) * qJD(1) * t266 + qJDD(1) * t259 + t210 * t303 - t211 * t300 - t231 * t297 - t353 * t282;
t64 = -qJDD(3) * pkin(3) - t69;
t38 = t198 * pkin(4) - t199 * qJ(5) - qJD(3) * t314 + t64;
t7 = -t17 * t222 + t226 * t38;
t8 = t226 * t17 + t222 * t38;
t386 = -t222 * t7 + t226 * t8;
t384 = mrSges(5,1) + t412;
t383 = mrSges(5,2) - t411;
t381 = m(4) + m(3) + t385;
t107 = -qJD(3) * pkin(3) + t239;
t274 = t222 * mrSges(6,1) + t226 * mrSges(6,2);
t72 = -qJD(4) * pkin(4) + qJD(5) - t76;
t379 = -t233 * t72 * t274 - t35 * (-mrSges(6,2) * t230 - mrSges(6,3) * t327) - t34 * (mrSges(6,1) * t230 - mrSges(6,3) * t322) - t107 * (mrSges(5,1) * t230 + mrSges(5,2) * t233);
t378 = -m(6) * t72 - t395;
t5 = pkin(5) * t198 - pkin(10) * t156 + t7;
t6 = pkin(10) * t155 + t8;
t1 = qJD(6) * t9 + t229 * t5 + t232 * t6;
t2 = -qJD(6) * t10 - t229 * t6 + t232 * t5;
t377 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t352 = pkin(5) * t222;
t305 = pkin(9) + t352;
t376 = -m(6) * pkin(9) - m(7) * t305 - t217 * mrSges(7,1) - t218 * mrSges(7,2) + mrSges(4,2) - t274 - t407;
t53 = -pkin(5) * t186 + t72;
t75 = -mrSges(7,1) * t286 + mrSges(7,2) * t123;
t374 = -m(7) * t53 + t378 - t75;
t18 = -qJDD(4) * pkin(4) + qJDD(5) - t22;
t11 = -pkin(5) * t155 + t18;
t25 = -t62 * mrSges(7,1) + t61 * mrSges(7,2);
t373 = -m(5) * t22 + m(6) * t18 + m(7) * t11 + t25 + t396;
t372 = -m(5) * t76 - t374;
t234 = qJD(3) ^ 2;
t371 = Ifges(7,4) * t369 + Ifges(7,2) * t368 + Ifges(7,6) * t358;
t370 = Ifges(7,1) * t369 + Ifges(7,4) * t368 + Ifges(7,5) * t358;
t365 = Ifges(6,4) * t360 + t408;
t363 = t286 / 0.2e1;
t361 = t123 / 0.2e1;
t355 = t213 / 0.2e1;
t344 = Ifges(5,4) * t233;
t340 = t18 * t230;
t332 = mrSges(4,2) * qJD(3);
t325 = t223 * t231;
t320 = mrSges(4,2) * qJDD(3);
t312 = Ifges(7,5) * t61 + Ifges(7,6) * t62 + Ifges(7,3) * t188;
t307 = mrSges(5,3) * t317;
t295 = -t317 / 0.2e1;
t292 = t224 * t337;
t288 = -t313 / 0.2e1;
t287 = t313 / 0.2e1;
t270 = Ifges(5,5) * t233 - Ifges(5,6) * t230;
t269 = Ifges(6,5) * t226 - Ifges(6,6) * t222;
t265 = t353 * t292;
t128 = t225 * t325 - t259 - t265;
t129 = t231 * t292 + t247;
t249 = t336 * t337 - t306;
t93 = t129 * t233 + t230 * t249;
t45 = t128 * t226 - t222 * t93;
t46 = t128 * t222 + t226 * t93;
t15 = -t229 * t46 + t232 * t45;
t16 = t229 * t45 + t232 * t46;
t171 = t230 * t336 + t233 * t324;
t133 = -t222 * t171 - t226 * t303;
t134 = t226 * t171 - t222 * t303;
t80 = t133 * t232 - t134 * t229;
t81 = t133 * t229 + t134 * t232;
t261 = t230 * (Ifges(5,1) * t233 - t345);
t235 = t224 * t245 - t291 * t336;
t165 = t223 * t278 + t227 * t334;
t89 = t165 * t353 - t231 * t392;
t47 = t230 * t89 - t233 * t235;
t236 = t224 * t246 + t290 * t336;
t166 = -t223 * t277 + t227 * t335;
t91 = t166 * t353 - t231 * t391;
t49 = t230 * t91 - t233 * t236;
t92 = t129 * t230 - t233 * t249;
t260 = -g(1) * t49 - g(2) * t47 - g(3) * t92;
t170 = t230 * t324 - t233 * t336;
t242 = t233 * (Ifges(6,3) * t230 + t233 * t269);
t216 = Ifges(5,4) * t317;
t207 = -qJD(4) * mrSges(5,2) + t307;
t196 = t305 * t230;
t179 = t305 * t315;
t175 = Ifges(5,1) * t319 + Ifges(5,5) * qJD(4) + t216;
t162 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t198;
t161 = t267 * t230;
t160 = t190 * t230;
t152 = -mrSges(6,1) * t317 - mrSges(6,3) * t187;
t151 = mrSges(6,2) * t317 + mrSges(6,3) * t186;
t149 = -pkin(9) * t327 + t183;
t135 = -qJD(4) * t170 + t233 * t283;
t119 = t129 * qJD(3);
t118 = (t265 + (t264 - t325) * t225) * qJD(3);
t116 = t187 * Ifges(6,1) + t186 * Ifges(6,4) - Ifges(6,5) * t317;
t115 = t187 * Ifges(6,4) + t186 * Ifges(6,2) - Ifges(6,6) * t317;
t113 = mrSges(6,1) * t198 - mrSges(6,3) * t156;
t112 = -mrSges(6,2) * t198 + mrSges(6,3) * t155;
t106 = -qJD(4) * t258 + t168 * t230;
t105 = -qJD(4) * t257 - t169 * t230;
t104 = t135 * t226 + t222 * t300;
t103 = -t135 * t222 + t226 * t300;
t100 = mrSges(7,1) * t213 - mrSges(7,3) * t123;
t99 = -mrSges(7,2) * t213 + mrSges(7,3) * t286;
t90 = t166 * t231 + t353 * t391;
t88 = t165 * t231 + t353 * t392;
t82 = t156 * Ifges(6,4) + t155 * Ifges(6,2) + t198 * Ifges(6,6);
t70 = t330 + (qJD(3) * t352 + t108) * t233;
t50 = t230 * t236 + t91 * t233;
t48 = t230 * t235 + t89 * t233;
t44 = -qJD(4) * t92 + t118 * t233;
t41 = -mrSges(7,2) * t188 + mrSges(7,3) * t62;
t40 = mrSges(7,1) * t188 - mrSges(7,3) * t61;
t33 = t119 * t222 + t226 * t44;
t32 = t119 * t226 - t222 * t44;
t27 = -qJD(6) * t81 + t103 * t232 - t104 * t229;
t26 = qJD(6) * t80 + t103 * t229 + t104 * t232;
t4 = -qJD(6) * t16 - t229 * t33 + t232 * t32;
t3 = qJD(6) * t15 + t229 * t32 + t232 * t33;
t12 = [t93 * t162 + t33 * t151 + t32 * t152 + t46 * t112 + t45 * t113 + t3 * t99 + t4 * t100 + m(5) * (t21 * t93 + t44 * t77) + m(6) * (t32 * t34 + t33 * t35 + t45 * t7 + t46 * t8) + m(7) * (t1 * t16 + t10 * t3 + t15 * t2 + t4 * t9) + m(4) * (t111 * t118 + t68 * t129 + t146 * t249) + m(3) * (t210 * t337 + (t223 ^ 2 + t227 ^ 2) * t225 ^ 2 * qJDD(1)) - t118 * t332 - t129 * t320 + t44 * t207 + t15 * t40 + t16 * t41 + m(2) * qJDD(1) + (-m(4) * t69 + m(5) * t64 - t389) * t128 + (m(4) * t239 + m(5) * t107 - t388) * t119 + t373 * t92 + t372 * (qJD(4) * t93 + t118 * t230) + (-m(2) - t381) * g(3); t104 * t151 + t103 * t152 + t133 * t113 + t134 * t112 + t26 * t99 + t27 * t100 + t80 * t40 + t81 * t41 + m(3) * t210 + m(6) * (t103 * t34 + t104 * t35 + t133 * t7 + t134 * t8) + m(7) * (t1 * t81 + t10 * t26 + t2 * t80 + t27 * t9) + m(4) * (t146 * t336 + (t353 * t69 + t231 * t68 + (t111 * t353 + t231 * t239) * qJD(3)) * t224) + m(5) * (t77 * t135 + t21 * t171 + (t107 * t318 - t353 * t64) * t224) + t171 * t162 + t135 * t207 + t191 * t300 + (-mrSges(4,1) * t234 - t320) * t324 + (-mrSges(4,2) * t234 + t389) * t303 + t373 * t170 + t372 * (qJD(4) * t171 + t230 * t283) + (-g(1) * t290 + g(2) * t291 - g(3) * t337) * t381; (-pkin(3) * t64 - t107 * t111) * m(5) + ((-Ifges(5,2) * t230 + t344) * t287 - Ifges(7,6) * t368 - Ifges(7,5) * t369 - Ifges(6,3) * t357 - Ifges(7,3) * t358 - Ifges(6,5) * t359 - Ifges(6,6) * t360 + t8 * mrSges(6,2) - t7 * mrSges(6,1) + Ifges(5,4) * t405 + Ifges(5,2) * t406 - t377) * t233 + (m(5) * (-t230 * t76 + t233 * t77) - t332 - t374 * t230 + t233 * t207) * t239 + (t376 * t91 + t410 * t90) * g(1) + (t376 * t89 + t410 * t88) * g(2) + (t128 * t410 + t129 * t376) * g(3) + (-t77 * mrSges(5,3) + Ifges(7,5) * t361 + Ifges(7,6) * t363 + Ifges(7,3) * t355 - t409) * t316 + (Ifges(7,1) * t105 + Ifges(7,4) * t106) * t361 + t388 * t111 + t274 * t340 + t261 * t287 + t242 * t288 + (-Ifges(7,1) * t161 - Ifges(7,4) * t160) * t369 + t11 * (mrSges(7,1) * t160 - mrSges(7,2) * t161) + (-Ifges(7,4) * t161 - Ifges(7,2) * t160) * t368 + (-t1 * t160 + t10 * t106 - t105 * t9 + t161 * t2) * mrSges(7,3) + (-Ifges(7,5) * t161 - Ifges(7,6) * t160) * t358 - pkin(3) * t140 + t149 * t113 + t150 * t112 + (-g(1) * t91 - g(2) * t89 - g(3) * t129 - t315 * t76 + t387) * mrSges(5,3) + t53 * (-mrSges(7,1) * t106 + mrSges(7,2) * t105) + t78 * t40 + t79 * t41 - t68 * mrSges(4,2) + t69 * mrSges(4,1) + (t149 * t7 + t150 * t8 + t398 * t34 + t397 * t35) * m(6) + (t387 * m(5) + t395 * t315 + t396 * t230 + m(6) * (t315 * t72 + t340) + t233 * t162) * pkin(9) - (Ifges(6,5) * t156 + Ifges(6,6) * t155 + Ifges(6,3) * t198 + t312) * t233 / 0.2e1 + t64 * t276 + (Ifges(7,4) * t105 + Ifges(7,2) * t106) * t363 + (t226 * t116 + t175) * t315 / 0.2e1 + (-t323 * t7 - t328 * t8) * mrSges(6,3) + t105 * t366 + t106 * t367 - t161 * t370 - t160 * t371 + t323 * t365 - t222 * t115 * t315 / 0.2e1 + Ifges(4,3) * qJDD(3) - t82 * t328 / 0.2e1 - t207 * t311 + t390 * t404 + t344 * t405 + t272 * t406 + t179 * t75 + t196 * t25 + t397 * t151 + t398 * t152 + t402 * t99 + t403 * t100 + (t1 * t79 + t10 * t402 + t11 * t196 + t179 * t53 + t2 * t78 + t403 * t9) * m(7) + (Ifges(5,1) * t199 + Ifges(5,4) * t406 + t269 * t357 + t271 * t360 + t273 * t359) * t230 + ((-t230 * t77 - t233 * t76) * t407 + t270 * t404 - t379) * qJD(4) + qJDD(4) * (Ifges(5,5) * t230 + Ifges(5,6) * t233) + (Ifges(7,5) * t105 + Ifges(7,6) * t106) * t355; (Ifges(7,5) * t190 - Ifges(7,6) * t267) * t358 + t11 * (mrSges(7,1) * t267 + mrSges(7,2) * t190) - t267 * t371 + (-Ifges(7,4) * t154 - Ifges(7,2) * t153) * t364 + (-Ifges(7,1) * t154 - Ifges(7,4) * t153) * t362 + (-Ifges(7,5) * t154 - Ifges(7,6) * t153) * t356 + (Ifges(7,4) * t190 - Ifges(7,2) * t267) * t368 + (Ifges(7,1) * t190 - Ifges(7,4) * t267) * t369 + (-t1 * t267 + t10 * t413 - t190 * t2 - t414 * t9) * mrSges(7,3) + (-mrSges(7,1) * t413 + mrSges(7,2) * t414) * t53 + (-qJ(5) * t113 + t365 + (-m(6) * t34 - t152) * qJD(5) + t408) * t222 + t413 * t367 + t414 * t366 + (-pkin(4) * t18 + qJ(5) * t386 - t34 * t51 - t35 * t52) * m(6) + t386 * mrSges(6,3) + t270 * t288 + (-Ifges(7,5) * t168 - Ifges(7,6) * t169) * t355 + (-Ifges(7,1) * t168 - Ifges(7,4) * t169) * t361 + (-Ifges(7,4) * t168 - Ifges(7,2) * t169) * t363 + t139 * t41 - t52 * t151 - t51 * t152 + t138 * t40 - pkin(4) * t95 - t70 * t75 + (-t261 / 0.2e1 + t242 / 0.2e1) * t234 + t18 * t275 + t343 * t360 + (-t390 / 0.2e1 + t379) * qJD(3) + (Ifges(6,6) * t357 + Ifges(6,2) * t360 + qJ(5) * t112 + t116 * t295 + t82 / 0.2e1 + (m(6) * t35 + t151) * qJD(5)) * t226 + t190 * t370 + t342 * t359 + (t308 + t378) * t77 + (t383 * t93 + t384 * t92) * g(3) + (t383 * t48 + t384 * t47) * g(2) + (t383 * t50 + t384 * t49) * g(1) + Ifges(5,3) * qJDD(4) + (t175 + t216) * t295 + t115 * t301 / 0.2e1 + (-Ifges(5,2) * t295 + t409 + t417) * t319 + (-t207 + t307) * t76 - Ifges(5,6) * t198 + Ifges(5,5) * t199 - t215 * t25 + t400 * t100 + t401 * t99 + (t1 * t139 + t10 * t401 - t11 * t215 + t138 * t2 + t400 * t9 - t53 * t70) * m(7) - t21 * mrSges(5,2) + t22 * mrSges(5,1); t123 * t100 - t286 * t99 - t186 * t151 + t187 * t152 + t25 + t95 + (-t10 * t286 + t123 * t9 + t11 + t260) * m(7) + (-t186 * t35 + t187 * t34 + t18 + t260) * m(6); -t53 * (mrSges(7,1) * t123 + mrSges(7,2) * t286) + (Ifges(7,1) * t286 - t341) * t362 + t66 * t361 + (Ifges(7,5) * t286 - Ifges(7,6) * t123) * t356 - t9 * t99 + t10 * t100 - g(1) * ((-t217 * t50 + t218 * t90) * mrSges(7,1) + (-t217 * t90 - t218 * t50) * mrSges(7,2)) - g(2) * ((-t217 * t48 + t218 * t88) * mrSges(7,1) + (-t217 * t88 - t218 * t48) * mrSges(7,2)) - g(3) * ((t128 * t218 - t217 * t93) * mrSges(7,1) + (-t128 * t217 - t218 * t93) * mrSges(7,2)) + (t10 * t123 + t286 * t9) * mrSges(7,3) + t312 + (-Ifges(7,2) * t123 + t117 + t67) * t364 + t377;];
tau  = t12;
