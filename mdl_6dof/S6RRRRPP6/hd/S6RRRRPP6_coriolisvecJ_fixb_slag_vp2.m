% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2018-11-23 18:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:07:56
% EndTime: 2018-11-23 18:08:20
% DurationCPUTime: 24.55s
% Computational Cost: add. (9102->757), mult. (22484->983), div. (0->0), fcn. (14858->6), ass. (0->326)
t262 = cos(qJ(3));
t389 = cos(qJ(4));
t313 = t389 * qJD(3);
t300 = t262 * t313;
t312 = t389 * qJD(4);
t259 = sin(qJ(4));
t260 = sin(qJ(3));
t360 = t259 * t260;
t419 = qJD(3) + qJD(4);
t159 = -t262 * t312 + t360 * t419 - t300;
t263 = cos(qJ(2));
t336 = t263 * qJD(1);
t322 = t260 * t336;
t324 = t389 * t262;
t179 = -t259 * t322 + t324 * t336;
t351 = t159 + t179;
t261 = sin(qJ(2));
t345 = qJD(1) * t261;
t323 = t260 * t345;
t343 = qJD(2) * t262;
t280 = t323 - t343;
t186 = t389 * t280;
t321 = t262 * t345;
t337 = t260 * qJD(2);
t205 = t321 + t337;
t340 = qJD(3) * t261;
t318 = t260 * t340;
t342 = qJD(2) * t263;
t278 = t262 * t342 - t318;
t334 = qJD(2) * qJD(3);
t308 = t262 * t334;
t267 = qJD(1) * t278 + t308;
t320 = t263 * t337;
t339 = qJD(3) * t262;
t277 = t261 * t339 + t320;
t309 = t260 * t334;
t268 = qJD(1) * t277 + t309;
t338 = qJD(4) * t259;
t66 = qJD(4) * t186 + t205 * t338 + t259 * t268 - t267 * t389;
t412 = -t66 / 0.2e1;
t411 = t66 / 0.2e1;
t271 = t205 * t389 - t259 * t280;
t67 = qJD(4) * t271 + t259 * t267 + t268 * t389;
t409 = t67 / 0.2e1;
t335 = qJD(1) * qJD(2);
t311 = t261 * t335;
t465 = t311 / 0.2e1;
t436 = Ifges(5,1) + Ifges(7,3);
t435 = Ifges(7,4) + Ifges(6,5);
t434 = Ifges(5,5) + Ifges(7,5);
t433 = Ifges(7,2) + Ifges(6,3);
t432 = Ifges(6,6) - Ifges(7,6);
t431 = Ifges(7,6) - Ifges(5,4);
t460 = Ifges(6,4) - t434;
t459 = -Ifges(5,6) + t435;
t359 = t259 * t262;
t207 = t260 * t389 + t359;
t160 = t419 * t207;
t178 = t207 * t336;
t350 = t160 - t178;
t253 = pkin(7) * t336;
t200 = pkin(3) * t322 + t253;
t341 = qJD(3) * t260;
t464 = pkin(3) * t341 - t200;
t463 = -Ifges(3,6) * qJD(2) / 0.2e1;
t462 = t431 * t409 + t436 * t412 + t434 * t465;
t461 = t433 * t409 + t411 * t432 + t435 * t465;
t407 = -pkin(9) - pkin(8);
t325 = qJD(3) * t407;
t210 = t260 * t325;
t223 = t407 * t260;
t224 = t407 * t262;
t111 = -t389 * t210 - t223 * t312 - t224 * t338 - t325 * t359;
t296 = pkin(2) * t261 - pkin(8) * t263;
t209 = t296 * qJD(1);
t170 = pkin(7) * t323 + t262 * t209;
t355 = t262 * t263;
t286 = pkin(3) * t261 - pkin(9) * t355;
t135 = qJD(1) * t286 + t170;
t191 = t260 * t209;
t356 = t261 * t262;
t357 = t260 * t263;
t152 = t191 + (-pkin(7) * t356 - pkin(9) * t357) * qJD(1);
t87 = t259 * t135 + t389 * t152;
t74 = -qJ(5) * t345 - t87;
t458 = t111 - t74;
t213 = -pkin(2) * t263 - t261 * pkin(8) - pkin(1);
t197 = t213 * qJD(1);
t222 = qJD(2) * pkin(8) + t253;
t155 = t197 * t262 - t222 * t260;
t121 = -pkin(9) * t205 + t155;
t156 = t260 * t197 + t262 * t222;
t122 = -pkin(9) * t280 + t156;
t361 = t259 * t122;
t52 = t121 * t389 - t361;
t457 = pkin(3) * t312 - t52;
t169 = t259 * t223 - t224 * t389;
t112 = qJD(4) * t169 + t259 * t210 - t407 * t300;
t293 = -t135 * t389 + t259 * t152;
t456 = t112 - t293;
t455 = qJ(5) * t351 - qJD(5) * t207 + t464;
t149 = t205 * t259 + t186;
t446 = t149 * pkin(5) - qJD(6);
t418 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t206 = -t324 + t360;
t384 = pkin(4) + qJ(6);
t452 = qJD(6) * t206 + t350 * t384 + t455;
t314 = t261 * t384;
t451 = -pkin(5) * t351 + qJD(1) * t314 + t456;
t450 = -pkin(5) * t350 - t458;
t141 = Ifges(7,6) * t271;
t244 = qJD(3) - t336;
t226 = -qJD(4) - t244;
t373 = Ifges(6,6) * t271;
t429 = t149 * t433 - t226 * t435 + t141 - t373;
t143 = Ifges(5,4) * t149;
t372 = Ifges(7,6) * t149;
t428 = -t226 * t434 + t271 * t436 - t143 + t372;
t117 = t389 * t122;
t51 = t121 * t259 + t117;
t449 = pkin(3) * t338 - t51;
t423 = -qJD(5) - t457;
t448 = pkin(4) * t350 + t455;
t364 = qJ(5) * t149;
t126 = mrSges(6,1) * t271 - mrSges(6,2) * t226;
t382 = mrSges(5,3) * t271;
t128 = -mrSges(5,1) * t226 - t382;
t352 = t126 - t128;
t273 = pkin(3) * t244 + t121;
t113 = t389 * t273;
t41 = -t113 + t361;
t437 = pkin(5) * t271;
t284 = t41 + t437;
t447 = qJD(5) + t284;
t252 = pkin(7) * t345;
t221 = -qJD(2) * pkin(2) + t252;
t174 = pkin(3) * t280 + t221;
t272 = t259 * t273;
t42 = t117 + t272;
t39 = t226 * qJ(5) - t42;
t30 = -t39 - t446;
t266 = -qJ(5) * t271 + t174;
t36 = t149 * t384 + t266;
t56 = t149 * pkin(4) + t266;
t376 = Ifges(5,4) * t271;
t82 = -Ifges(5,2) * t149 - Ifges(5,6) * t226 + t376;
t443 = t30 * mrSges(7,1) + t42 * mrSges(5,3) + t56 * mrSges(6,2) + t82 / 0.2e1 - t174 * mrSges(5,1) - t36 * mrSges(7,3) - t39 * mrSges(6,1);
t29 = t226 * t384 + t447;
t421 = -qJD(5) - t41;
t38 = pkin(4) * t226 - t421;
t142 = Ifges(6,6) * t149;
t81 = -Ifges(6,4) * t226 - Ifges(6,2) * t271 + t142;
t442 = mrSges(6,1) * t38 + mrSges(7,1) * t29 + mrSges(5,2) * t174 + mrSges(5,3) * t41 - t81 / 0.2e1 + t428 / 0.2e1 - mrSges(7,2) * t36 - mrSges(6,3) * t56;
t380 = Ifges(3,4) * t261;
t441 = t155 * mrSges(4,1) + t30 * mrSges(7,2) + t38 * mrSges(6,2) + t463 - (t263 * Ifges(3,2) + t380) * qJD(1) / 0.2e1 - t156 * mrSges(4,2) - t29 * mrSges(7,3) - t39 * mrSges(6,3) - t41 * mrSges(5,1) - t42 * mrSges(5,2);
t410 = -t67 / 0.2e1;
t399 = -t271 / 0.2e1;
t439 = -t311 / 0.2e1;
t344 = qJD(2) * t261;
t307 = t344 / 0.2e1;
t438 = pkin(4) * t271;
t46 = mrSges(6,1) * t67 - mrSges(6,3) * t311;
t47 = -t67 * mrSges(7,1) + mrSges(7,2) * t311;
t430 = -t46 + t47;
t427 = Ifges(4,5) * t260;
t425 = t449 + t446;
t424 = t437 - t423;
t422 = t271 * t384;
t420 = -t42 + t446;
t246 = pkin(7) * t355;
t177 = t260 * t213 + t246;
t415 = t418 * t311 + t459 * t67 + t460 * t66;
t414 = Ifges(6,4) * t439 + Ifges(6,2) * t412 + Ifges(6,6) * t410;
t413 = Ifges(5,4) * t411 + Ifges(5,2) * t409 + Ifges(5,6) * t439;
t378 = Ifges(4,4) * t260;
t381 = Ifges(4,1) * t262;
t291 = -t378 + t381;
t404 = t291 * t334 / 0.2e1 + (Ifges(4,1) * t278 - Ifges(4,4) * t277 + Ifges(4,5) * t344) * qJD(1) / 0.2e1;
t402 = -t149 / 0.2e1;
t401 = t149 / 0.2e1;
t398 = t271 / 0.2e1;
t395 = t205 / 0.2e1;
t392 = -t226 / 0.2e1;
t391 = t226 / 0.2e1;
t388 = pkin(3) * t205;
t387 = pkin(3) * t259;
t256 = t261 * pkin(7);
t383 = mrSges(5,3) * t149;
t379 = Ifges(4,4) * t205;
t377 = Ifges(4,4) * t262;
t375 = Ifges(4,5) * t205;
t374 = Ifges(4,6) * t244;
t371 = Ifges(4,3) * t244;
t370 = t156 * mrSges(4,3);
t369 = t205 * Ifges(4,1);
t368 = t244 * Ifges(4,5);
t367 = Ifges(3,5) * qJD(2);
t366 = Ifges(4,2) * qJD(3);
t363 = qJD(2) * mrSges(3,1);
t362 = qJD(2) * mrSges(3,2);
t358 = t260 * t261;
t124 = mrSges(6,1) * t149 + mrSges(6,3) * t226;
t125 = -mrSges(7,1) * t149 - mrSges(7,2) * t226;
t354 = -t124 + t125;
t127 = mrSges(5,2) * t226 - t383;
t353 = t124 - t127;
t204 = t262 * t213;
t154 = -pkin(9) * t356 + t204 + (-pkin(7) * t260 - pkin(3)) * t263;
t162 = -pkin(9) * t358 + t177;
t98 = t259 * t154 + t389 * t162;
t349 = mrSges(4,1) * t280 + t205 * mrSges(4,2) + mrSges(3,3) * t345 - t363;
t211 = t296 * qJD(2);
t198 = qJD(1) * t211;
t348 = t197 * t339 + t260 * t198;
t302 = pkin(7) * t311;
t347 = t262 * t198 + t260 * t302;
t346 = t262 * t211 + t337 * t256;
t212 = pkin(3) * t358 + t256;
t330 = t389 * pkin(3);
t255 = pkin(7) * t342;
t328 = t259 * t358;
t327 = -t378 / 0.2e1;
t326 = Ifges(4,5) * t267 - Ifges(4,6) * t268 + Ifges(4,3) * t311;
t175 = pkin(3) * t277 + t255;
t250 = -pkin(3) * t262 - pkin(2);
t319 = t261 * t343;
t310 = t263 * t335;
t306 = t342 / 0.2e1;
t305 = -t340 / 0.2e1;
t48 = -t66 * mrSges(6,1) + mrSges(6,2) * t311;
t168 = -t389 * t223 - t224 * t259;
t243 = pkin(7) * t310;
t249 = -t330 - pkin(4);
t301 = t389 * t342;
t298 = qJD(2) * t314;
t89 = qJ(5) * t263 - t98;
t184 = t261 * t324 - t328;
t294 = -qJ(5) * t184 + t212;
t97 = t154 * t389 - t259 * t162;
t292 = mrSges(4,1) * t260 + mrSges(4,2) * t262;
t290 = -Ifges(4,2) * t260 + t377;
t289 = Ifges(4,5) * t262 - Ifges(4,6) * t260;
t288 = t388 + t364;
t287 = -t155 * t262 - t156 * t260;
t45 = -t66 * mrSges(7,1) - mrSges(7,3) * t311;
t285 = -qJ(5) * t207 + t250;
t90 = t263 * pkin(4) - t97;
t279 = t286 * qJD(2);
t55 = qJD(1) * t279 - qJD(3) * t122 + t347;
t274 = t277 * pkin(9);
t71 = (-pkin(9) * qJD(2) - t222) * t341 + (-pkin(7) * t319 - t274) * qJD(1) + t348;
t282 = qJD(4) * t272 + t122 * t312 + t259 * t71 - t389 * t55;
t119 = t260 * t211 + t213 * t339 + (-t263 * t341 - t319) * pkin(7);
t101 = -t274 + t119;
t96 = t279 + (-t246 + (pkin(9) * t261 - t213) * t260) * qJD(3) + t346;
t281 = t259 * t101 + t154 * t338 + t162 * t312 - t389 * t96;
t9 = qJD(4) * t113 - t122 * t338 + t259 * t55 + t389 * t71;
t17 = t389 * t101 + t154 * t312 - t162 * t338 + t259 * t96;
t276 = t280 * mrSges(4,3);
t104 = t160 * t261 + t259 * t320 - t262 * t301;
t275 = qJ(5) * t104 - qJD(5) * t184 + t175;
t201 = Ifges(4,4) * t280;
t134 = -t201 + t368 + t369;
t270 = t155 * mrSges(4,3) - t134 / 0.2e1 - t369 / 0.2e1 - t221 * mrSges(4,2) - t368 / 0.2e1;
t133 = -Ifges(4,2) * t280 + t374 + t379;
t269 = t370 + t133 / 0.2e1 + t379 / 0.2e1 + t374 / 0.2e1 - t221 * mrSges(4,1);
t15 = -qJ(5) * t344 + qJD(5) * t263 - t17;
t5 = -qJ(5) * t311 + qJD(5) * t226 - t9;
t147 = pkin(3) * t268 + t243;
t1 = -t66 * pkin(5) - qJD(1) * t298 + t226 * qJD(6) + t282;
t3 = -pkin(5) * t67 - t5;
t7 = -pkin(4) * t311 + t282;
t265 = -mrSges(5,1) * t282 - t9 * mrSges(5,2) + t7 * mrSges(6,2) + t3 * mrSges(7,2) - t5 * mrSges(6,3) - t1 * mrSges(7,3) + t415;
t264 = t66 * qJ(5) - qJD(5) * t271 + t147;
t251 = Ifges(3,4) * t336;
t247 = qJ(5) + t387;
t242 = -qJ(6) + t249;
t219 = mrSges(3,3) * t336 - t362;
t189 = Ifges(3,1) * t345 + t251 + t367;
t183 = t207 * t261;
t176 = -pkin(7) * t357 + t204;
t173 = mrSges(4,1) * t244 - mrSges(4,3) * t205;
t172 = -t244 * mrSges(4,2) - t276;
t171 = -pkin(7) * t321 + t191;
t146 = -mrSges(4,3) * t309 + (-mrSges(4,2) * t344 - mrSges(4,3) * t277) * qJD(1);
t145 = -mrSges(4,3) * t308 + (mrSges(4,1) * t344 - mrSges(4,3) * t278) * qJD(1);
t144 = pkin(4) * t206 + t285;
t132 = -Ifges(4,6) * t280 + t371 + t375;
t131 = -t206 * pkin(5) + t169;
t130 = pkin(5) * t207 + t168;
t123 = mrSges(7,1) * t271 + mrSges(7,3) * t226;
t120 = -qJD(3) * t177 + t346;
t118 = t206 * t384 + t285;
t115 = pkin(4) * t183 + t294;
t114 = t292 * t334 + (mrSges(4,1) * t277 + mrSges(4,2) * t278) * qJD(1);
t106 = t290 * t334 + (Ifges(4,4) * t278 - Ifges(4,2) * t277 + Ifges(4,6) * t344) * qJD(1);
t105 = t260 * t301 - t259 * t318 - qJD(4) * t328 + (t259 * t342 + (t313 + t312) * t261) * t262;
t103 = -qJD(3) * t156 + t347;
t102 = -t222 * t341 - t262 * t302 + t348;
t95 = -mrSges(6,2) * t149 - mrSges(6,3) * t271;
t94 = mrSges(5,1) * t149 + mrSges(5,2) * t271;
t93 = t364 + t438;
t92 = -mrSges(7,2) * t271 + mrSges(7,3) * t149;
t91 = t183 * t384 + t294;
t84 = -t226 * Ifges(6,1) - Ifges(6,4) * t271 + t149 * Ifges(6,5);
t83 = -t226 * Ifges(7,1) + t149 * Ifges(7,4) + Ifges(7,5) * t271;
t79 = Ifges(5,5) * t271 - t149 * Ifges(5,6) - t226 * Ifges(5,3);
t75 = -pkin(4) * t345 + t293;
t73 = t288 + t438;
t68 = -t183 * pkin(5) - t89;
t54 = t184 * pkin(5) + t263 * qJ(6) + t90;
t50 = -mrSges(5,2) * t311 - mrSges(5,3) * t67;
t49 = mrSges(5,1) * t311 + mrSges(5,3) * t66;
t44 = t364 + t422;
t37 = t288 + t422;
t28 = pkin(4) * t105 + t275;
t27 = mrSges(7,2) * t66 + mrSges(7,3) * t67;
t26 = mrSges(5,1) * t67 - mrSges(5,2) * t66;
t25 = -mrSges(6,2) * t67 + mrSges(6,3) * t66;
t16 = -pkin(4) * t344 + t281;
t14 = qJD(6) * t183 + t105 * t384 + t275;
t13 = t67 * pkin(4) + t264;
t12 = -t105 * pkin(5) - t15;
t11 = -t104 * pkin(5) + t263 * qJD(6) + t281 - t298;
t4 = t149 * qJD(6) + t384 * t67 + t264;
t2 = [(t132 + t84 + t83 + t79 + (t261 * t289 - t460 * t184 + t459 * t183 + (-Ifges(4,3) - t418) * t263) * qJD(1)) * t307 + (-Ifges(5,4) * t402 + Ifges(6,2) * t399 + t392 * t460 - t436 * t398 + t432 * t401 - t442) * t104 + (t459 * t392 + t429 / 0.2e1 + t433 * t401 + t431 * t398 + Ifges(6,6) * t399 - Ifges(5,2) * t402 - t443) * t105 - t282 * (-mrSges(5,1) * t263 - t184 * mrSges(5,3)) + m(5) * (t147 * t212 + t17 * t42 + t174 * t175 + t281 * t41 - t282 * t97 + t9 * t98) - t281 * t128 + ((-0.3e1 / 0.2e1 * t261 ^ 2 + 0.3e1 / 0.2e1 * t263 ^ 2) * Ifges(3,4) - 0.2e1 * (mrSges(3,1) * t261 + mrSges(3,2) * t263) * pkin(1)) * t335 + (Ifges(6,4) * t399 + Ifges(5,6) * t402 + t418 * t392 + t434 * t398 + t435 * t401 + t441) * t344 + t189 * t306 + (t102 * t263 + t221 * t278) * mrSges(4,2) - t106 * t358 / 0.2e1 + m(7) * (t1 * t54 + t11 * t29 + t12 * t30 + t14 * t36 + t3 * t68 + t4 * t91) + m(6) * (t115 * t13 + t15 * t39 + t16 * t38 + t28 * t56 + t5 * t89 + t7 * t90) + (t287 * t342 + (-t102 * t260 - t103 * t262 + (t155 * t260 - t156 * t262) * qJD(3)) * t261) * mrSges(4,3) - (t326 + t415) * t263 / 0.2e1 + (t260 * t305 + t262 * t306) * t134 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t261 * t310 + (t261 * t114 + (t349 * t263 + (t292 * t336 - t219) * t261) * qJD(2)) * pkin(7) + m(4) * (t102 * t177 + t103 * t176 + t156 * t119 + t155 * t120 + (t221 + t252) * t255) - t280 * ((-Ifges(4,2) * t262 - t378) * t340 + (Ifges(4,6) * t261 + t263 * t290) * qJD(2)) / 0.2e1 - t268 * (-Ifges(4,6) * t263 + t261 * t290) / 0.2e1 + t267 * (-Ifges(4,5) * t263 + t261 * t291) / 0.2e1 + t3 * (-t183 * mrSges(7,1) - mrSges(7,2) * t263) + t9 * (mrSges(5,2) * t263 - t183 * mrSges(5,3)) + t1 * (t184 * mrSges(7,1) + mrSges(7,3) * t263) + t7 * (t184 * mrSges(6,1) - mrSges(6,2) * t263) + t5 * (t183 * mrSges(6,1) + mrSges(6,3) * t263) + qJD(2) ^ 2 * (Ifges(3,5) * t263 - Ifges(3,6) * t261) / 0.2e1 + t244 * ((-Ifges(4,6) * t262 - t427) * t340 + (Ifges(4,3) * t261 + t263 * t289) * qJD(2)) / 0.2e1 + t54 * t45 + (t183 * t433 - t184 * t432 - t263 * t435) * t409 + (t183 * t431 + t184 * t436 - t263 * t434) * t412 + t68 * t47 + t89 * t46 + t90 * t48 + t91 * t27 + t14 * t92 + t28 * t95 + t97 * t49 + t98 * t50 + t115 * t25 + (-t320 / 0.2e1 + t262 * t305) * t133 + (-t103 * t263 + t221 * t277) * mrSges(4,1) + t11 * t123 + t15 * t124 + t12 * t125 + t16 * t126 + t17 * t127 + t119 * t172 + t120 * t173 + t175 * t94 + t176 * t145 + t177 * t146 + t4 * (-mrSges(7,2) * t184 + mrSges(7,3) * t183) + t147 * (mrSges(5,1) * t183 + mrSges(5,2) * t184) + t13 * (-mrSges(6,2) * t183 - mrSges(6,3) * t184) + t212 * t26 + ((-Ifges(4,1) * t260 - t377) * t340 + (Ifges(4,5) * t261 + t263 * t291) * qJD(2)) * t395 + t356 * t404 + (Ifges(5,4) * t184 - Ifges(5,2) * t183 - Ifges(5,6) * t263) * t410 + (-Ifges(6,4) * t263 - Ifges(6,2) * t184 + Ifges(6,6) * t183) * t411 + t183 * t413 + t184 * t414 + t183 * t461 + t184 * t462; t448 * t95 + ((t367 / 0.2e1 - t251 / 0.2e1 + pkin(1) * mrSges(3,2) * qJD(1) - t189 / 0.2e1 + (-t349 - t363) * pkin(7) + (-qJD(2) * pkin(7) * mrSges(4,1) + t270) * t262 + ((t381 / 0.2e1 + pkin(7) * mrSges(4,2) + t327) * qJD(2) + t269) * t260) * t263 + (-t441 + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1) * t226 - t371 / 0.2e1 - t132 / 0.2e1 + (-Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1 + Ifges(5,6) / 0.2e1) * t149 + t463 + (Ifges(4,6) * t345 / 0.2e1 + (-Ifges(4,2) * t336 / 0.2e1 + t366 / 0.2e1 - qJD(3) * Ifges(4,1) / 0.2e1) * t260 + (t336 / 0.2e1 - 0.3e1 / 0.2e1 * qJD(3)) * t377) * t260 + (t219 + t362) * pkin(7) - t262 ^ 2 * t366 / 0.2e1 - t84 / 0.2e1 - t83 / 0.2e1 - t79 / 0.2e1 + (pkin(1) * mrSges(3,1) + t380 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t263) * qJD(1) + (-Ifges(7,5) / 0.2e1 - Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t271 - t375 / 0.2e1) * t261 + (t206 * t459 - t207 * t460 + t427) * t307) * qJD(1) + (t159 * t460 + t160 * t459) * t392 + (t178 * t459 - t179 * t460) * t391 + (t282 * t168 + t147 * t250 + t169 * t9 + (-t111 - t87) * t42 + t456 * t41 + t464 * t174) * m(5) + (-t206 * t9 + t207 * t282 - t350 * t42 - t351 * t41) * mrSges(5,3) + t293 * t128 + (-t103 * mrSges(4,3) - pkin(8) * t145 + t404) * t260 + (mrSges(7,2) * t351 + mrSges(7,3) * t350) * t36 + (-mrSges(6,2) * t350 + mrSges(6,3) * t351) * t56 + (t206 * t5 + t207 * t7 + t350 * t39 - t351 * t38) * mrSges(6,1) + (mrSges(5,1) * t350 - mrSges(5,2) * t351) * t174 + (t1 * t207 - t206 * t3 - t29 * t351 - t30 * t350) * mrSges(7,1) + t352 * t112 + t353 * t111 + (-t155 * t170 - t156 * t171 - t221 * t253 + (qJD(3) * t287 + t102 * t262 - t103 * t260) * pkin(8) - pkin(2) * t243) * m(4) + ((Ifges(4,4) * t343 - pkin(8) * t173 - t270) * t262 + (pkin(3) * t94 - pkin(8) * t172 + (t327 + (-Ifges(4,2) + Ifges(4,1) / 0.2e1) * t262) * qJD(2) - t269) * t260) * qJD(3) + (Ifges(6,2) * t159 + Ifges(6,6) * t160 + t178 * t431 + t179 * t436) * t399 + (t206 * t431 + t207 * t436) * t412 + (-Ifges(6,2) * t179 + Ifges(6,6) * t178 - t159 * t436 + t160 * t431) * t398 + (t13 * t144 + t168 * t7 - t169 * t5 + t448 * t56 + t458 * t39 + (t112 - t75) * t38) * m(6) + (t50 - t46) * t169 + (t106 / 0.2e1 + t102 * mrSges(4,3) + pkin(8) * t146) * t262 + (t81 - t428) * (t159 / 0.2e1 + t179 / 0.2e1) + (-t82 + t429) * (t160 / 0.2e1 - t178 / 0.2e1) + t450 * t125 + t451 * t123 + (t1 * t130 + t118 * t4 + t131 * t3 + t29 * t451 + t30 * t450 + t36 * t452) * m(7) + t452 * t92 + (Ifges(5,4) * t179 - Ifges(5,2) * t178 + t159 * t432 + t160 * t433) * t401 + (-Ifges(5,4) * t159 - Ifges(5,2) * t160 + t178 * t433 - t179 * t432) * t402 + (t206 * t433 - t207 * t432) * t409 + (t48 - t49) * t168 - pkin(2) * t114 + t118 * t27 - t74 * t124 - t75 * t126 - t87 * t127 + t130 * t45 + t131 * t47 + t144 * t25 - t171 * t172 - t170 * t173 - t200 * t94 + t4 * (-mrSges(7,2) * t207 + mrSges(7,3) * t206) + t147 * (mrSges(5,1) * t206 + mrSges(5,2) * t207) + t13 * (-mrSges(6,2) * t206 - mrSges(6,3) * t207) + t250 * t26 + (Ifges(5,4) * t207 - Ifges(5,2) * t206) * t410 + (-Ifges(6,2) * t207 + Ifges(6,6) * t206) * t411 + t206 * t413 + t207 * t414 + t206 * t461 + t207 * t462; (m(6) * t38 + t352) * t449 + (-Ifges(5,4) * t401 + Ifges(6,2) * t398 + t391 * t460 - t436 * t399 + t432 * t402 + t442) * t149 + (-Ifges(5,2) * t401 + Ifges(6,6) * t398 + t391 * t459 + t431 * t399 + t433 * t402 + t443) * t271 + (-Ifges(4,2) * t205 + t134 - t201) * t280 / 0.2e1 + (-t174 * t388 - t41 * t51 - t42 * t52 + (-t389 * t282 + t259 * t9 + (t259 * t41 + t389 * t42) * qJD(4)) * pkin(3)) * m(5) - t94 * t388 + (-t276 - t172) * t155 + t326 + (-t247 * t5 + t249 * t7 + t39 * t423 - t56 * t73) * m(6) + t457 * t127 + t265 - t205 * (-Ifges(4,1) * t280 - t379) / 0.2e1 + t423 * t124 + t424 * t125 + t425 * t123 + (t1 * t242 + t247 * t3 + t29 * t425 + t30 * t424 - t36 * t37) * m(7) + t430 * t247 + t49 * t330 + t205 * t370 + t429 * t399 - t37 * t92 - t73 * t95 - t102 * mrSges(4,2) + t103 * mrSges(4,1) - t221 * (t205 * mrSges(4,1) - mrSges(4,2) * t280) - t244 * (-Ifges(4,5) * t280 - Ifges(4,6) * t205) / 0.2e1 + t156 * t173 + t242 * t45 + t249 * t48 + t50 * t387 + t133 * t395; t354 * qJD(5) + t265 + (t149 * t29 + t271 * t30) * mrSges(7,1) + (t149 * t38 - t271 * t39) * mrSges(6,1) + t430 * qJ(5) + (-t352 + t382) * t42 + (-t353 + t383) * t41 + t420 * t123 - pkin(4) * t48 - t44 * t92 - t93 * t95 + t284 * t125 - t56 * (-mrSges(6,2) * t271 + mrSges(6,3) * t149) - t36 * (mrSges(7,2) * t149 + mrSges(7,3) * t271) - t174 * (mrSges(5,1) * t271 - mrSges(5,2) * t149) - t384 * t45 + (Ifges(6,2) * t149 + t373 + t82) * t398 + (qJ(5) * t3 - t1 * t384 + t420 * t29 + t30 * t447 - t36 * t44) * m(7) + (-pkin(4) * t7 - qJ(5) * t5 - t38 * t42 + t39 * t421 - t56 * t93) * m(6) + (t271 * t433 + t142 - t372 + t81) * t402 + (-Ifges(5,2) * t271 - t143 + t428) * t401 + (t149 * t460 + t271 * t459) * t391 + (-t149 * t436 + t141 - t376 + t429) * t399; t354 * t226 + (t92 + t95) * t271 + t45 + t48 + (t226 * t30 + t271 * t36 + t1) * m(7) + (-t226 * t39 + t271 * t56 + t7) * m(6); -t226 * t123 - t149 * t92 + 0.2e1 * (t3 / 0.2e1 + t36 * t402 + t29 * t392) * m(7) + t47;];
tauc  = t2(:);
