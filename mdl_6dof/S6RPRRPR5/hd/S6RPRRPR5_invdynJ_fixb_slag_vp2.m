% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:12:14
% EndTime: 2019-03-09 05:12:47
% DurationCPUTime: 21.00s
% Computational Cost: add. (12309->673), mult. (30113->831), div. (0->0), fcn. (23360->14), ass. (0->321)
t471 = mrSges(5,1) - mrSges(6,2);
t239 = sin(pkin(10));
t235 = t239 ^ 2;
t240 = cos(pkin(10));
t236 = t240 ^ 2;
t333 = t235 + t236;
t326 = qJD(1) * qJD(2);
t213 = qJ(2) * qJDD(1) + t326;
t242 = sin(qJ(6));
t246 = cos(qJ(6));
t457 = mrSges(7,1) * t242 + mrSges(7,2) * t246;
t244 = sin(qJ(3));
t396 = cos(qJ(3));
t272 = t244 * t239 - t240 * t396;
t180 = t272 * qJD(1);
t187 = t239 * t396 + t244 * t240;
t181 = t187 * qJD(1);
t243 = sin(qJ(4));
t395 = cos(qJ(4));
t142 = t395 * t180 + t181 * t243;
t238 = qJD(3) + qJD(4);
t124 = t142 * t242 + t238 * t246;
t137 = Ifges(5,4) * t142;
t271 = -t243 * t180 + t181 * t395;
t138 = qJD(6) + t271;
t461 = t138 * Ifges(7,3);
t123 = t142 * t246 - t238 * t242;
t462 = t123 * Ifges(7,6);
t441 = Ifges(5,1) * t271 + t238 * Ifges(5,5) + t124 * Ifges(7,5) - t137 + t461 + t462;
t377 = pkin(7) + qJ(2);
t207 = t377 * t239;
t192 = qJD(1) * t207;
t208 = t377 * t240;
t193 = qJD(1) * t208;
t152 = -t396 * t192 - t244 * t193;
t121 = -pkin(8) * t181 + t152;
t119 = qJD(3) * pkin(3) + t121;
t153 = -t244 * t192 + t193 * t396;
t122 = -t180 * pkin(8) + t153;
t340 = t243 * t122;
t84 = -t395 * t119 + t340;
t76 = -pkin(4) * t238 + qJD(5) + t84;
t136 = Ifges(6,6) * t142;
t97 = t238 * Ifges(6,4) - Ifges(6,2) * t271 + t136;
t470 = t76 * mrSges(6,1) + t84 * mrSges(5,3) - t97 / 0.2e1 + t441 / 0.2e1;
t463 = Ifges(6,5) - Ifges(5,6);
t303 = qJD(4) * t395;
t88 = t121 * t395 - t340;
t469 = pkin(3) * t303 - t88;
t369 = mrSges(5,3) * t142;
t126 = -mrSges(5,2) * t238 - t369;
t376 = mrSges(6,1) * t142;
t128 = -mrSges(6,3) * t238 + t376;
t118 = t395 * t122;
t85 = t243 * t119 + t118;
t82 = -t238 * qJ(5) - t85;
t425 = -m(6) * t82 + t126 - t128;
t468 = -m(5) * t85 - t425;
t368 = mrSges(5,3) * t271;
t375 = mrSges(6,1) * t271;
t436 = t238 * t471 - t368 - t375;
t424 = -m(6) * t76 + t436;
t467 = m(5) + m(6);
t466 = m(5) * t84;
t233 = qJDD(3) + qJDD(4);
t400 = t233 / 0.2e1;
t385 = t142 * pkin(5);
t443 = Ifges(6,4) - Ifges(5,5);
t364 = Ifges(4,4) * t187;
t266 = Ifges(4,4) * t272;
t439 = -qJD(5) - t469;
t353 = qJ(5) * t142;
t284 = mrSges(7,1) * t246 - mrSges(7,2) * t242;
t46 = -t82 - t385;
t361 = t124 * Ifges(7,4);
t55 = t123 * Ifges(7,2) + t138 * Ifges(7,6) + t361;
t460 = t46 * t284 - t246 * t55 / 0.2e1;
t237 = pkin(10) + qJ(3);
t232 = qJ(4) + t237;
t222 = cos(t232);
t458 = t457 * t222;
t221 = sin(t232);
t456 = -t471 * t222 + (mrSges(5,2) - mrSges(6,3)) * t221;
t444 = pkin(5) * t271;
t274 = t84 + t444;
t455 = t274 + qJD(5);
t454 = 0.2e1 * t400;
t399 = -t238 / 0.2e1;
t403 = t271 / 0.2e1;
t404 = -t271 / 0.2e1;
t409 = -t138 / 0.2e1;
t411 = -t124 / 0.2e1;
t412 = -t123 / 0.2e1;
t224 = t240 * pkin(2) + pkin(1);
t196 = -qJD(1) * t224 + qJD(2);
t154 = pkin(3) * t180 + t196;
t414 = pkin(4) + pkin(9);
t45 = -t238 * t414 + t455;
t256 = -qJ(5) * t271 + t154;
t51 = t142 * t414 + t256;
t22 = -t242 * t51 + t246 * t45;
t23 = t242 * t45 + t246 * t51;
t86 = pkin(4) * t142 + t256;
t446 = -t22 * mrSges(7,1) - t154 * mrSges(5,2) + t23 * mrSges(7,2) + t86 * mrSges(6,3);
t452 = -Ifges(5,1) * t404 - Ifges(7,5) * t411 + Ifges(6,2) * t403 - Ifges(7,6) * t412 - Ifges(7,3) * t409 + t443 * t399 - t446;
t330 = qJD(4) * t243;
t299 = pkin(7) * qJDD(1) + t213;
t173 = t299 * t239;
t174 = t299 * t240;
t104 = -qJD(3) * t153 - t396 * t173 - t244 * t174;
t149 = -qJD(3) * t180 + qJDD(1) * t187;
t79 = qJDD(3) * pkin(3) - t149 * pkin(8) + t104;
t304 = qJD(3) * t396;
t331 = qJD(3) * t244;
t103 = -t244 * t173 + t396 * t174 - t192 * t304 - t193 * t331;
t261 = qJD(3) * t187;
t258 = qJD(1) * t261;
t253 = -qJDD(1) * t272 - t258;
t83 = pkin(8) * t253 + t103;
t19 = -t119 * t330 - t122 * t303 - t243 * t83 + t395 * t79;
t265 = qJDD(5) - t19;
t80 = -t395 * t149 + t180 * t303 + t181 * t330 - t243 * t253;
t10 = -t80 * pkin(5) - t233 * t414 + t265;
t229 = -qJDD(1) * pkin(1) + qJDD(2);
t324 = qJDD(1) * t240;
t81 = qJD(4) * t271 + t243 * t149 - t253 * t395;
t21 = -pkin(2) * t324 - pkin(3) * t253 + t81 * pkin(4) + t80 * qJ(5) - qJD(5) * t271 + t229;
t12 = t81 * pkin(9) + t21;
t1 = qJD(6) * t22 + t10 * t242 + t12 * t246;
t2 = -qJD(6) * t23 + t10 * t246 - t12 * t242;
t451 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t450 = -t154 * mrSges(5,1) + t86 * mrSges(6,2);
t359 = t271 * Ifges(6,6);
t96 = t238 * Ifges(6,5) + t142 * Ifges(6,3) - t359;
t360 = t271 * Ifges(5,4);
t98 = -t142 * Ifges(5,2) + t238 * Ifges(5,6) + t360;
t449 = t85 * mrSges(5,3) - t82 * mrSges(6,1) - t96 / 0.2e1 + t98 / 0.2e1;
t230 = sin(t237);
t231 = cos(t237);
t288 = mrSges(4,1) * t231 - mrSges(4,2) * t230;
t289 = -t240 * mrSges(3,1) + mrSges(3,2) * t239;
t448 = -m(3) * pkin(1) - m(4) * t224 - mrSges(2,1) - t288 + t289 + t456;
t234 = -pkin(8) - t377;
t394 = m(3) * qJ(2);
t447 = -m(4) * t377 - m(7) * (pkin(5) - t234) - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - t394;
t42 = qJD(6) * t123 + t233 * t246 + t242 * t81;
t418 = t42 / 0.2e1;
t43 = -qJD(6) * t124 - t233 * t242 + t246 * t81;
t417 = t43 / 0.2e1;
t75 = qJDD(6) - t80;
t416 = t75 / 0.2e1;
t445 = pkin(4) * t271;
t380 = qJD(3) / 0.2e1;
t13 = -mrSges(7,1) * t43 + mrSges(7,2) * t42;
t64 = mrSges(6,1) * t81 - mrSges(6,3) * t233;
t442 = -t64 + t13;
t440 = t444 - t439;
t89 = -mrSges(7,1) * t123 + mrSges(7,2) * t124;
t438 = t89 - t128;
t437 = t271 * t414;
t156 = -t244 * t207 + t396 * t208;
t247 = cos(qJ(1));
t343 = t222 * t247;
t345 = t221 * t247;
t435 = pkin(4) * t343 + qJ(5) * t345;
t199 = qJ(5) * t343;
t432 = -m(7) * t199 - mrSges(6,2) * t345 - mrSges(6,3) * t343 - t247 * t458;
t245 = sin(qJ(1));
t344 = t222 * t245;
t197 = qJ(5) * t344;
t346 = t221 * t245;
t431 = -m(7) * t197 - mrSges(6,2) * t346 - mrSges(6,3) * t344 - t245 * t458;
t26 = mrSges(7,1) * t75 - mrSges(7,3) * t42;
t27 = -mrSges(7,2) * t75 + mrSges(7,3) * t43;
t430 = t242 * t27 + t246 * t26;
t429 = g(1) * t247 + g(2) * t245;
t428 = -t222 * mrSges(7,3) - t221 * t457 + t456;
t386 = g(3) * t222;
t426 = -t221 * t429 + t386;
t277 = t22 * t242 - t23 * t246;
t255 = -qJD(6) * t277 + t1 * t242 + t2 * t246;
t328 = qJD(6) * t246;
t329 = qJD(6) * t242;
t91 = -mrSges(7,2) * t138 + mrSges(7,3) * t123;
t92 = mrSges(7,1) * t138 - mrSges(7,3) * t124;
t423 = m(7) * t255 + t91 * t328 - t92 * t329 + t430;
t278 = t22 * t246 + t23 * t242;
t422 = -m(7) * t278 + t424;
t280 = Ifges(7,5) * t242 + Ifges(7,6) * t246;
t363 = Ifges(7,4) * t242;
t281 = Ifges(7,2) * t246 + t363;
t362 = Ifges(7,4) * t246;
t282 = Ifges(7,1) * t242 + t362;
t397 = -t242 / 0.2e1;
t406 = t142 / 0.2e1;
t407 = -t142 / 0.2e1;
t120 = Ifges(7,4) * t123;
t56 = t124 * Ifges(7,1) + t138 * Ifges(7,5) + t120;
t420 = -Ifges(5,2) * t406 + Ifges(6,3) * t407 + t280 * t409 + t281 * t412 + t282 * t411 + t56 * t397 + t399 * t463 + t450 + t460;
t419 = Ifges(7,1) * t418 + Ifges(7,4) * t417 + Ifges(7,5) * t416;
t410 = t124 / 0.2e1;
t365 = Ifges(4,4) * t181;
t408 = -Ifges(4,2) * t180 / 0.2e1 + Ifges(4,6) * t380 + t365 / 0.2e1;
t398 = t238 / 0.2e1;
t393 = pkin(3) * t181;
t392 = pkin(3) * t230;
t219 = pkin(3) * t231;
t391 = pkin(3) * t243;
t217 = t222 * pkin(4);
t371 = mrSges(4,3) * t180;
t370 = mrSges(4,3) * t181;
t367 = mrSges(7,3) * t242;
t366 = mrSges(7,3) * t246;
t357 = t235 * mrSges(3,3);
t356 = t236 * mrSges(3,3);
t262 = qJD(3) * t272;
t264 = t243 * t272;
t106 = -qJD(4) * t264 + t187 * t303 - t243 * t262 + t261 * t395;
t352 = t106 * t242;
t351 = t106 * t246;
t350 = t271 * t242;
t349 = t271 * t246;
t259 = t395 * t272;
t150 = t187 * t243 + t259;
t348 = t150 * t242;
t347 = t150 * t246;
t209 = t221 * qJ(5);
t342 = t242 * t245;
t341 = t242 * t247;
t339 = t245 * t246;
t338 = t246 * t247;
t305 = qJD(2) * t396;
t337 = -t207 * t304 + t240 * t305;
t334 = t217 + t209;
t332 = qJD(2) * t244;
t325 = qJDD(1) * t239;
t323 = Ifges(7,5) * t42 + Ifges(7,6) * t43 + Ifges(7,3) * t75;
t322 = t395 * pkin(3);
t313 = t219 + t334;
t307 = t81 * mrSges(5,1) - t80 * mrSges(5,2);
t306 = -t81 * mrSges(6,2) + t80 * mrSges(6,3);
t302 = -t329 / 0.2e1;
t65 = -t80 * mrSges(6,1) + t233 * mrSges(6,2);
t194 = t219 + t224;
t298 = -t194 - t209;
t87 = t121 * t243 + t118;
t155 = -t396 * t207 - t244 * t208;
t134 = -pkin(8) * t187 + t155;
t135 = -pkin(8) * t272 + t156;
t94 = -t395 * t134 + t135 * t243;
t182 = t247 * t194;
t297 = -t234 * t245 + t182;
t295 = -m(7) * t414 - mrSges(7,3);
t225 = -t322 - pkin(4);
t292 = -pkin(4) * t221 - t392;
t290 = -mrSges(3,1) * t324 + mrSges(3,2) * t325;
t285 = mrSges(5,1) * t221 + mrSges(5,2) * t222;
t279 = t393 + t353;
t151 = t187 * t395 - t264;
t162 = pkin(3) * t272 - t224;
t254 = -t151 * qJ(5) + t162;
t59 = t150 * t414 + t254;
t60 = pkin(5) * t151 + t94;
t34 = t242 * t60 + t246 * t59;
t33 = -t242 * t59 + t246 * t60;
t276 = -t242 * t92 + t246 * t91;
t275 = t295 * t221;
t95 = t243 * t134 + t135 * t395;
t270 = t150 * t328 + t352;
t269 = t150 * t329 - t351;
t18 = t119 * t303 - t122 * t330 + t243 * t79 + t395 * t83;
t112 = -pkin(8) * t261 - t208 * t331 - t239 * t332 + t337;
t251 = pkin(8) * t262 + t207 * t331 - t208 * t304 - t239 * t305 - t240 * t332;
t35 = -t395 * t112 - t134 * t303 + t135 * t330 - t243 * t251;
t260 = pkin(3) * t261;
t257 = -m(7) * t392 + t275;
t14 = -qJ(5) * t233 - qJD(5) * t238 - t18;
t36 = qJD(4) * t95 + t243 * t112 - t395 * t251;
t105 = t187 * t330 + t238 * t259 + t243 * t261;
t37 = t106 * pkin(4) + t105 * qJ(5) - t151 * qJD(5) + t260;
t252 = -mrSges(4,1) * t253 + t149 * mrSges(4,2);
t11 = -pkin(5) * t81 - t14;
t16 = -t233 * pkin(4) + t265;
t8 = t42 * Ifges(7,4) + t43 * Ifges(7,2) + t75 * Ifges(7,6);
t250 = t8 * t397 + t56 * t302 - t18 * mrSges(5,2) + t19 * mrSges(5,1) + t16 * mrSges(6,2) - t14 * mrSges(6,3) + t11 * t457 - t2 * t366 - t1 * t367 + (Ifges(7,5) * t246 - Ifges(7,6) * t242) * t416 + (-Ifges(7,2) * t242 + t362) * t417 + (Ifges(7,1) * t246 - t363) * t418 + t246 * t419 + t463 * t81 + t443 * t80 + (Ifges(6,1) + Ifges(5,3)) * t233 + t460 * qJD(6) + (t22 * t329 - t23 * t328) * mrSges(7,3) - (t123 * t281 + t124 * t282 + t138 * t280) * qJD(6) / 0.2e1;
t223 = qJ(5) + t391;
t216 = t222 * pkin(9);
t195 = -qJDD(1) * t224 + qJDD(2);
t176 = Ifges(4,4) * t180;
t172 = -t221 * t342 + t338;
t171 = t221 * t339 + t341;
t170 = t221 * t341 + t339;
t169 = t221 * t338 - t342;
t158 = qJD(3) * mrSges(4,1) - t370;
t157 = -qJD(3) * mrSges(4,2) - t371;
t140 = Ifges(4,1) * t181 + Ifges(4,5) * qJD(3) - t176;
t132 = -t187 * qJD(2) - qJD(3) * t156;
t131 = (-qJD(2) * t239 - qJD(3) * t208) * t244 + t337;
t125 = pkin(3) * t258 + qJDD(1) * t162 + qJDD(2);
t102 = -mrSges(6,2) * t142 - mrSges(6,3) * t271;
t101 = mrSges(5,1) * t142 + mrSges(5,2) * t271;
t100 = t353 + t445;
t93 = t150 * pkin(4) + t254;
t90 = t279 + t445;
t66 = t353 + t437;
t63 = -mrSges(5,2) * t233 - mrSges(5,3) * t81;
t62 = mrSges(5,1) * t233 + mrSges(5,3) * t80;
t61 = -t150 * pkin(5) + t95;
t57 = t279 + t437;
t52 = t87 - t385;
t50 = t85 - t385;
t32 = t242 * t50 + t246 * t66;
t31 = -t242 * t66 + t246 * t50;
t30 = t106 * pkin(9) + t37;
t29 = t242 * t52 + t246 * t57;
t28 = -t242 * t57 + t246 * t52;
t25 = -t105 * pkin(5) + t36;
t24 = -pkin(5) * t106 - t35;
t4 = -qJD(6) * t34 - t242 * t30 + t246 * t25;
t3 = qJD(6) * t33 + t242 * t25 + t246 * t30;
t5 = [(t125 * mrSges(5,1) + t14 * mrSges(6,1) - t21 * mrSges(6,2) - t18 * mrSges(5,3) - t11 * t284 + t280 * t416 + t281 * t417 + t282 * t418 + t55 * t302 + (Ifges(5,2) + Ifges(6,3)) * t81 + (Ifges(5,4) + Ifges(6,6)) * t80 + t463 * t454) * t150 + (-Ifges(5,1) * t403 + Ifges(6,2) * t404 - t462 / 0.2e1 - t461 / 0.2e1 + Ifges(6,6) * t406 - Ifges(5,4) * t407 - Ifges(7,5) * t410 + t443 * t398 + t446 - t470) * t105 + m(6) * (t21 * t93 + t37 * t86) + (-Ifges(5,4) * t403 - Ifges(5,2) * t407 + Ifges(6,6) * t404 + Ifges(6,3) * t406 + t398 * t463 - t449 - t450) * t106 + (-m(7) * (pkin(9) * t343 + t182 + t435) - t170 * mrSges(7,1) - t169 * mrSges(7,2) - mrSges(7,3) * t343 - m(6) * (t297 + t435) - m(5) * t297 + t448 * t247 + t447 * t245) * g(2) + (t466 - t424) * t36 + t468 * t35 + (-t103 * t272 - t104 * t187 + t152 * t262 - t153 * t261 + t156 * t253) * mrSges(4,3) + (-mrSges(4,3) * t155 + Ifges(4,1) * t187 - t266) * t149 + (-t172 * mrSges(7,1) + t171 * mrSges(7,2) + (t234 * t467 + t447) * t247 + (-m(6) * (t298 - t217) - m(7) * t298 - t222 * t295 + m(5) * t194 - t448) * t245) * g(1) + t196 * (mrSges(4,1) * t187 - mrSges(4,2) * t272) * qJD(3) + t333 * t213 * mrSges(3,3) - t140 * t262 / 0.2e1 + t37 * t102 + t24 * t89 + t3 * t91 + t4 * t92 + (m(5) * t18 - m(6) * t14 + t63 - t64) * t95 - t80 * Ifges(5,1) * t151 - t224 * t252 + qJD(3) ^ 2 * (-Ifges(4,5) * t272 - Ifges(4,6) * t187) / 0.2e1 + t195 * (mrSges(4,1) * t272 + t187 * mrSges(4,2)) + (Ifges(3,4) * t239 + Ifges(3,2) * t240) * t324 + (Ifges(3,1) * t239 + Ifges(3,4) * t240) * t325 + t213 * t356 + t213 * t357 + t46 * (mrSges(7,1) * t269 + mrSges(7,2) * t270) + (-m(5) * t19 + m(6) * t16 - t62 + t65) * t94 + t123 * (Ifges(7,4) * t270 - Ifges(7,2) * t269) / 0.2e1 + t101 * t260 + (-t180 * (-Ifges(4,2) * t187 - t266) + t181 * (-Ifges(4,1) * t272 - t364)) * t380 + (mrSges(4,1) * t155 - mrSges(4,2) * t156 + Ifges(4,5) * t187 - Ifges(4,6) * t272) * qJDD(3) + t138 * (Ifges(7,5) * t270 - Ifges(7,6) * t269) / 0.2e1 + t61 * t13 + (t16 * mrSges(6,1) + t125 * mrSges(5,2) - t19 * mrSges(5,3) - t21 * mrSges(6,3) + Ifges(5,5) * t400 + Ifges(7,5) * t418 + Ifges(7,6) * t417 + Ifges(7,3) * t416 + (-Ifges(5,4) / 0.2e1 - Ifges(6,6)) * t81 - Ifges(6,2) * t80 - t454 * Ifges(6,4) + t451) * t151 + (-Ifges(5,4) * t81 + Ifges(5,5) * t233 + t323) * t151 / 0.2e1 + t33 * t26 + t34 * t27 + Ifges(2,3) * qJDD(1) + m(3) * (-pkin(1) * t229 + (t213 + t326) * qJ(2) * t333) + t253 * (-Ifges(4,2) * t272 + t364) + t56 * t352 / 0.2e1 + t55 * t351 / 0.2e1 + t93 * t306 + t162 * t307 + m(5) * (t125 * t162 + t154 * t260) + (qJD(6) * t56 + t8) * t347 / 0.2e1 - t261 * t408 + t229 * t289 - pkin(1) * t290 + (Ifges(7,1) * t270 - Ifges(7,4) * t269) * t410 + t131 * t157 + t132 * t158 + t348 * t419 + m(7) * (t1 * t34 + t11 * t61 + t2 * t33 + t22 * t4 + t23 * t3 + t24 * t46) + m(4) * (t103 * t156 + t104 * t155 + t131 * t153 + t132 * t152 - t195 * t224) + (t1 * t347 - t2 * t348 - t22 * t270 - t23 * t269) * mrSges(7,3); t252 + m(7) * (-qJD(6) * t278 + t1 * t246 - t2 * t242) + m(6) * t21 + m(5) * t125 + m(3) * t229 + t290 + t180 * t157 + t181 * t158 + t306 + t307 - t242 * t26 + t246 * t27 + (-t349 - t328) * t92 + (-t350 - t329) * t91 + (t152 * t181 + t153 * t180 + t195) * m(4) + (-t333 * t394 - t356 - t357) * qJD(1) ^ 2 + (t422 - t466) * t271 - (-m(7) * t46 + t468 - t89) * t142 + (-g(1) * t245 + g(2) * t247) * (m(7) + m(4) + m(3) + t467); ((t395 * t19 + t18 * t243 + (t243 * t84 + t395 * t85) * qJD(4)) * pkin(3) - t154 * t393 - t84 * t87 - t85 * t88) * m(5) + (-Ifges(5,4) * t404 + Ifges(6,6) * t403 + t420 + t449) * t271 + (-Ifges(5,4) * t406 + Ifges(6,6) * t407 + t452 + t470) * t142 + (-t245 * t257 + t431) * g(2) + (-t247 * t257 + t432) * g(1) + t436 * t87 + (t242 * t91 + t246 * t92 - t422) * pkin(3) * t330 - t181 * (-Ifges(4,1) * t180 - t365) / 0.2e1 - qJD(3) * (-Ifges(4,5) * t180 - Ifges(4,6) * t181) / 0.2e1 - t196 * (mrSges(4,1) * t181 - mrSges(4,2) * t180) + t442 * t223 + t63 * t391 - t90 * t102 - t29 * t91 - t28 * t92 + t423 * (-pkin(9) + t225) + t62 * t322 + t250 + Ifges(4,6) * t253 + (-m(5) * t219 - m(6) * t313 - m(7) * (t216 + t313) - t288 + t428) * g(3) + (m(5) * t392 + mrSges(4,1) * t230 + mrSges(4,2) * t231 + t285) * t429 + (t370 + t158) * t153 + (-t76 * t87 - t86 * t90 - t14 * t223 + t16 * t225 - g(1) * (t247 * t292 + t199) - g(2) * (t245 * t292 + t197) + t439 * t82) * m(6) + t439 * t128 + (t11 * t223 - t22 * t28 - t23 * t29 + t440 * t46) * m(7) + t440 * t89 + (t22 * t350 - t23 * t349) * mrSges(7,3) - t101 * t393 + Ifges(4,3) * qJDD(3) + (-Ifges(4,2) * t181 + t140 - t176) * t180 / 0.2e1 + t181 * t408 + Ifges(4,5) * t149 + t225 * t65 + t469 * t126 - t103 * mrSges(4,2) + t104 * mrSges(4,1) + (-t371 - t157) * t152; (qJ(5) * t11 - t22 * t31 - t23 * t32 + t455 * t46) * m(7) + t452 * t142 + t429 * t285 + (-t245 * t275 + t431) * g(2) + (-t247 * t275 + t432) * g(1) + (t136 + t97) * t407 + (-pkin(4) * t16 - qJ(5) * t14 - qJD(5) * t82 - t100 * t86 - g(1) * (-pkin(4) * t345 + t199) - g(2) * (-pkin(4) * t346 + t197)) * m(6) + (t22 * t367 - t23 * t366 + t420) * t271 + (-t137 + t441) * t406 + t442 * qJ(5) + t76 * t376 - t100 * t102 - t32 * t91 - t31 * t92 - t423 * t414 + t250 + (-m(7) * (t216 + t334) - m(6) * t334 + t428) * g(3) + t438 * qJD(5) + (-t360 + t96) * t404 - pkin(4) * t65 - t82 * t375 + t274 * t89 + (t368 + t424) * t85 + (t369 + t425) * t84 + (t98 + t359) * t403; -t438 * t238 + t276 * qJD(6) + (t102 + t276) * t271 + t65 + (-t238 * t46 - t271 * t277 + t255 + t426) * m(7) + (t238 * t82 + t271 * t86 + t16 + t426) * m(6) + t430; -t46 * (mrSges(7,1) * t124 + mrSges(7,2) * t123) + (Ifges(7,1) * t123 - t361) * t411 + t55 * t410 + (Ifges(7,5) * t123 - Ifges(7,6) * t124) * t409 - t22 * t91 + t23 * t92 - g(1) * (mrSges(7,1) * t169 - mrSges(7,2) * t170) - g(2) * (mrSges(7,1) * t171 + mrSges(7,2) * t172) + t284 * t386 + (t123 * t22 + t124 * t23) * mrSges(7,3) + t323 + (-Ifges(7,2) * t124 + t120 + t56) * t412 + t451;];
tau  = t5;
