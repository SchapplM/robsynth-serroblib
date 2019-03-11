% Calculate time derivative of joint inertia matrix for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:34
% EndTime: 2019-03-10 01:35:00
% DurationCPUTime: 11.98s
% Computational Cost: add. (14321->781), mult. (36994->1099), div. (0->0), fcn. (36111->10), ass. (0->313)
t475 = Ifges(6,5) + Ifges(7,5);
t474 = Ifges(6,6) + Ifges(7,6);
t323 = cos(qJ(5));
t319 = sin(qJ(5));
t417 = Ifges(7,4) * t319;
t419 = Ifges(6,4) * t319;
t464 = t417 + t419 + (Ifges(6,2) + Ifges(7,2)) * t323;
t416 = Ifges(7,4) * t323;
t418 = Ifges(6,4) * t323;
t463 = t416 + t418 + (Ifges(6,1) + Ifges(7,1)) * t319;
t324 = cos(qJ(4));
t456 = (t319 ^ 2 + t323 ^ 2) * t324;
t318 = cos(pkin(6));
t321 = sin(qJ(3));
t325 = cos(qJ(3));
t317 = sin(pkin(6));
t322 = sin(qJ(2));
t398 = t317 * t322;
t252 = t318 * t325 - t321 * t398;
t253 = t318 * t321 + t325 * t398;
t320 = sin(qJ(4));
t187 = t252 * t320 + t253 * t324;
t326 = cos(qJ(2));
t384 = qJD(2) * t317;
t362 = t326 * t384;
t222 = -qJD(3) * t253 - t321 * t362;
t223 = qJD(3) * t252 + t325 * t362;
t109 = qJD(4) * t187 - t324 * t222 + t223 * t320;
t335 = t324 * t252 - t253 * t320;
t108 = qJD(4) * t335 + t222 * t320 + t223 * t324;
t397 = t317 * t326;
t162 = -t319 * t187 - t323 * t397;
t383 = qJD(2) * t322;
t363 = t317 * t383;
t59 = qJD(5) * t162 + t323 * t108 + t319 * t363;
t332 = -t323 * t187 + t319 * t397;
t60 = qJD(5) * t332 - t319 * t108 + t323 * t363;
t15 = Ifges(7,4) * t59 + Ifges(7,2) * t60 + Ifges(7,6) * t109;
t16 = Ifges(6,4) * t59 + Ifges(6,2) * t60 + Ifges(6,6) * t109;
t473 = t15 + t16;
t17 = Ifges(7,1) * t59 + Ifges(7,4) * t60 + Ifges(7,5) * t109;
t18 = Ifges(6,1) * t59 + Ifges(6,4) * t60 + Ifges(6,5) * t109;
t472 = t17 + t18;
t71 = -Ifges(7,4) * t332 + Ifges(7,2) * t162 - Ifges(7,6) * t335;
t72 = -Ifges(6,4) * t332 + Ifges(6,2) * t162 - Ifges(6,6) * t335;
t471 = t71 + t72;
t73 = -Ifges(7,1) * t332 + Ifges(7,4) * t162 - Ifges(7,5) * t335;
t74 = -Ifges(6,1) * t332 + Ifges(6,4) * t162 - Ifges(6,5) * t335;
t470 = t73 + t74;
t266 = t320 * t325 + t321 * t324;
t455 = qJD(3) + qJD(4);
t214 = t455 * t266;
t378 = qJD(5) * t319;
t361 = t266 * t378;
t265 = t320 * t321 - t324 * t325;
t213 = t455 * t265;
t407 = t213 * t323;
t330 = t361 + t407;
t377 = qJD(5) * t323;
t360 = t266 * t377;
t396 = t319 * t213;
t331 = t360 - t396;
t83 = -Ifges(7,4) * t330 - Ifges(7,2) * t331 + Ifges(7,6) * t214;
t84 = -Ifges(6,4) * t330 - Ifges(6,2) * t331 + Ifges(6,6) * t214;
t469 = t83 + t84;
t85 = -Ifges(7,1) * t330 - Ifges(7,4) * t331 + Ifges(7,5) * t214;
t86 = -Ifges(6,1) * t330 - Ifges(6,4) * t331 + Ifges(6,5) * t214;
t468 = t85 + t86;
t340 = -Ifges(7,2) * t319 + t416;
t166 = Ifges(7,6) * t265 + t266 * t340;
t341 = -Ifges(6,2) * t319 + t418;
t167 = Ifges(6,6) * t265 + t266 * t341;
t467 = t166 + t167;
t342 = Ifges(7,1) * t323 - t417;
t168 = Ifges(7,5) * t265 + t266 * t342;
t343 = Ifges(6,1) * t323 - t419;
t169 = Ifges(6,5) * t265 + t266 * t343;
t388 = t168 + t169;
t466 = (t340 + t341) * qJD(5);
t465 = (t342 + t343) * qJD(5);
t429 = t323 / 0.2e1;
t430 = t319 / 0.2e1;
t460 = t474 * t429 + t430 * t475;
t310 = Ifges(7,5) * t377;
t311 = Ifges(6,5) * t377;
t355 = -t378 / 0.2e1;
t459 = t311 / 0.2e1 + t310 / 0.2e1 + t474 * t355;
t13 = Ifges(7,5) * t59 + Ifges(7,6) * t60 + Ifges(7,3) * t109;
t14 = Ifges(6,5) * t59 + Ifges(6,6) * t60 + Ifges(6,3) * t109;
t458 = t13 + t14;
t307 = -pkin(3) * t325 - pkin(2);
t199 = pkin(4) * t265 - pkin(11) * t266 + t307;
t438 = -pkin(10) - pkin(9);
t294 = t438 * t321;
t295 = t438 * t325;
t234 = t294 * t320 - t295 * t324;
t224 = t323 * t234;
t141 = t319 * t199 + t224;
t260 = t318 * t322 * pkin(1) + pkin(8) * t397;
t243 = pkin(9) * t318 + t260;
t244 = (-pkin(2) * t326 - pkin(9) * t322 - pkin(1)) * t317;
t174 = t325 * t243 + t321 * t244;
t457 = t324 * t294 + t295 * t320;
t299 = pkin(8) * t398;
t428 = pkin(1) * t326;
t259 = t318 * t428 - t299;
t454 = Ifges(4,5) * t223 + Ifges(4,6) * t222 + Ifges(4,3) * t363;
t329 = t319 * t465 + t323 * t466 + t377 * t463 - t378 * t464;
t453 = 2 * m(5);
t452 = 2 * m(6);
t451 = 2 * m(7);
t450 = -2 * mrSges(3,3);
t449 = -2 * mrSges(5,3);
t364 = qJD(3) * t438;
t280 = t321 * t364;
t348 = t325 * t364;
t154 = qJD(4) * t234 + t280 * t320 - t324 * t348;
t448 = 0.2e1 * t154;
t447 = -0.2e1 * t457;
t249 = t260 * qJD(2);
t446 = 0.2e1 * t249;
t268 = mrSges(7,1) * t378 + mrSges(7,2) * t377;
t445 = 0.2e1 * t268;
t283 = -mrSges(7,1) * t323 + mrSges(7,2) * t319;
t444 = 0.2e1 * t283;
t443 = -0.2e1 * t319;
t442 = 0.2e1 * t323;
t441 = m(6) * pkin(4);
t437 = mrSges(7,3) * pkin(5);
t427 = pkin(3) * t324;
t173 = -t321 * t243 + t325 * t244;
t132 = -pkin(3) * t397 - t253 * pkin(10) + t173;
t149 = pkin(10) * t252 + t174;
t379 = qJD(4) * t324;
t380 = qJD(4) * t320;
t247 = (pkin(2) * t322 - pkin(9) * t326) * t384;
t248 = t259 * qJD(2);
t118 = -t174 * qJD(3) + t325 * t247 - t248 * t321;
t88 = pkin(3) * t363 - pkin(10) * t223 + t118;
t381 = qJD(3) * t325;
t382 = qJD(3) * t321;
t117 = -t243 * t382 + t244 * t381 + t321 * t247 + t325 * t248;
t92 = pkin(10) * t222 + t117;
t24 = t132 * t379 - t149 * t380 + t320 * t88 + t324 * t92;
t21 = pkin(11) * t363 + t24;
t172 = -t222 * pkin(3) + t249;
t41 = t109 * pkin(4) - t108 * pkin(11) + t172;
t80 = t320 * t132 + t324 * t149;
t68 = -pkin(11) * t397 + t80;
t242 = t299 + (-pkin(2) - t428) * t318;
t192 = -t252 * pkin(3) + t242;
t96 = -pkin(4) * t335 - t187 * pkin(11) + t192;
t5 = t323 * t21 + t319 * t41 + t96 * t377 - t378 * t68;
t426 = t323 * t5;
t37 = t319 * t96 + t323 * t68;
t6 = -qJD(5) * t37 - t21 * t319 + t323 * t41;
t425 = t6 * t319;
t423 = -qJ(6) - pkin(11);
t422 = mrSges(7,3) * t323;
t421 = Ifges(4,4) * t321;
t420 = Ifges(4,4) * t325;
t415 = pkin(3) * qJD(4);
t414 = t248 * mrSges(3,2);
t413 = t249 * mrSges(3,1);
t412 = t320 * mrSges(5,1);
t373 = pkin(3) * t382;
t135 = pkin(4) * t214 + pkin(11) * t213 + t373;
t153 = qJD(4) * t457 + t324 * t280 + t320 * t348;
t365 = t319 * t135 + t323 * t153 + t199 * t377;
t47 = -t234 * t378 + t365;
t411 = t323 * t47;
t410 = t324 * mrSges(5,2);
t352 = t323 * t135 - t153 * t319;
t48 = -qJD(5) * t141 + t352;
t409 = t48 * t319;
t408 = t154 * t457;
t406 = t457 * t320;
t405 = t266 * t319;
t404 = t266 * t323;
t304 = pkin(3) * t320 + pkin(11);
t402 = t304 * t319;
t401 = t304 * t323;
t395 = t319 * t324;
t284 = -mrSges(6,1) * t323 + mrSges(6,2) * t319;
t394 = t320 * t284;
t314 = t323 * qJ(6);
t124 = -mrSges(6,2) * t214 - mrSges(6,3) * t331;
t393 = t323 * t124;
t392 = t323 * t324;
t391 = -qJ(6) - t304;
t103 = -mrSges(6,1) * t162 - mrSges(6,2) * t332;
t171 = -mrSges(5,1) * t397 - t187 * mrSges(5,3);
t390 = t103 - t171;
t387 = -Ifges(7,5) * t407 + Ifges(7,3) * t214;
t386 = -Ifges(6,5) * t407 + Ifges(6,3) * t214;
t385 = t310 + t311;
t376 = 0.2e1 * t317;
t375 = 0.2e1 * qJD(5);
t372 = pkin(3) * t379;
t371 = pkin(5) * t378;
t370 = Ifges(4,6) * t397;
t369 = -t71 / 0.2e1 - t72 / 0.2e1;
t368 = t73 / 0.2e1 + t74 / 0.2e1;
t367 = m(7) * pkin(5) + mrSges(7,1);
t366 = Ifges(5,5) * t108 - Ifges(5,6) * t109 + Ifges(5,3) * t363;
t306 = -pkin(5) * t323 - pkin(4);
t27 = -t60 * mrSges(7,1) + t59 * mrSges(7,2);
t354 = t377 / 0.2e1;
t36 = -t319 * t68 + t323 * t96;
t353 = qJD(5) * t423;
t79 = t132 * t324 - t320 * t149;
t140 = t323 * t199 - t234 * t319;
t351 = qJD(5) * t391;
t350 = Ifges(5,5) * t363;
t349 = Ifges(5,6) * t363;
t345 = mrSges(6,3) * t456;
t67 = pkin(4) * t397 - t79;
t344 = mrSges(6,1) * t319 + mrSges(6,2) * t323;
t26 = -pkin(5) * t335 + qJ(6) * t332 + t36;
t29 = qJ(6) * t162 + t37;
t339 = -t26 * t323 - t29 * t319;
t338 = -t319 * t37 - t323 * t36;
t112 = pkin(5) * t265 - t266 * t314 + t140;
t125 = -qJ(6) * t405 + t141;
t337 = -t112 * t323 - t125 * t319;
t336 = -t140 * t323 - t141 * t319;
t334 = qJ(6) * t213 - qJD(6) * t266;
t25 = -t132 * t380 - t149 * t379 - t320 * t92 + t324 * t88;
t110 = mrSges(7,1) * t331 - mrSges(7,2) * t330;
t22 = -pkin(4) * t363 - t25;
t269 = t344 * qJD(5);
t3 = qJ(6) * t60 + qJD(6) * t162 + t5;
t49 = -pkin(5) * t162 + t67;
t8 = -pkin(5) * t60 + t22;
t328 = t25 * mrSges(5,1) - t24 * mrSges(5,2) + mrSges(6,3) * t426 + t22 * t284 + t49 * t268 + t67 * t269 + t8 * t283 + t3 * t422 + t366 + t463 * t59 / 0.2e1 + t464 * t60 / 0.2e1 + t466 * t162 / 0.2e1 - t465 * t332 / 0.2e1 + t472 * t430 + t473 * t429 + t471 * t355 + t470 * t354 - t459 * t335 + t460 * t109;
t101 = pkin(5) * t331 + t154;
t181 = pkin(5) * t405 - t457;
t207 = Ifges(5,6) * t214;
t208 = Ifges(5,5) * t213;
t38 = -qJ(6) * t360 + (-qJD(5) * t234 + t334) * t319 + t365;
t327 = -t153 * mrSges(5,2) + mrSges(6,3) * t411 + t101 * t283 + t181 * t268 - t457 * t269 + t38 * t422 - t207 - t208 + t468 * t430 + t469 * t429 - t466 * t405 / 0.2e1 + t465 * t404 / 0.2e1 + t467 * t355 + t388 * t354 + t459 * t265 + t460 * t214 + (t284 - mrSges(5,1)) * t154 + t464 * (-t360 / 0.2e1 + t396 / 0.2e1) + t463 * (t266 * t355 - t407 / 0.2e1);
t313 = t323 * qJD(6);
t312 = Ifges(4,5) * t381;
t305 = -pkin(4) - t427;
t298 = Ifges(3,5) * t362;
t293 = Ifges(4,1) * t321 + t420;
t290 = Ifges(4,2) * t325 + t421;
t285 = pkin(11) * t323 + t314;
t282 = t423 * t319;
t281 = t306 - t427;
t279 = pkin(3) * t380 + t371;
t278 = (Ifges(4,1) * t325 - t421) * qJD(3);
t275 = (-Ifges(4,2) * t321 + t420) * qJD(3);
t270 = (mrSges(4,1) * t321 + mrSges(4,2) * t325) * qJD(3);
t262 = t314 + t401;
t261 = t391 * t319;
t251 = -qJD(6) * t319 + t323 * t353;
t250 = t319 * t353 + t313;
t232 = -mrSges(4,1) * t397 - t253 * mrSges(4,3);
t231 = mrSges(4,2) * t397 + t252 * mrSges(4,3);
t221 = Ifges(5,1) * t266 - Ifges(5,4) * t265;
t220 = Ifges(5,4) * t266 - Ifges(5,2) * t265;
t219 = mrSges(5,1) * t265 + mrSges(5,2) * t266;
t212 = (-qJD(6) - t372) * t319 + t323 * t351;
t211 = t319 * t351 + t323 * t372 + t313;
t204 = mrSges(6,1) * t265 - mrSges(6,3) * t404;
t203 = mrSges(7,1) * t265 - mrSges(7,3) * t404;
t202 = -mrSges(6,2) * t265 - mrSges(6,3) * t405;
t201 = -mrSges(7,2) * t265 - mrSges(7,3) * t405;
t194 = t344 * t266;
t193 = (mrSges(7,1) * t319 + mrSges(7,2) * t323) * t266;
t185 = mrSges(4,1) * t363 - mrSges(4,3) * t223;
t184 = -mrSges(4,2) * t363 + mrSges(4,3) * t222;
t180 = Ifges(4,1) * t253 + Ifges(4,4) * t252 - Ifges(4,5) * t397;
t179 = Ifges(4,4) * t253 + Ifges(4,2) * t252 - t370;
t170 = mrSges(5,2) * t397 + mrSges(5,3) * t335;
t165 = Ifges(6,3) * t265 + (Ifges(6,5) * t323 - Ifges(6,6) * t319) * t266;
t164 = Ifges(7,3) * t265 + (Ifges(7,5) * t323 - Ifges(7,6) * t319) * t266;
t150 = -mrSges(4,1) * t222 + mrSges(4,2) * t223;
t148 = -Ifges(5,1) * t213 - Ifges(5,4) * t214;
t147 = -Ifges(5,4) * t213 - Ifges(5,2) * t214;
t146 = mrSges(5,1) * t214 - mrSges(5,2) * t213;
t134 = Ifges(4,1) * t223 + Ifges(4,4) * t222 + Ifges(4,5) * t363;
t133 = Ifges(4,4) * t223 + Ifges(4,2) * t222 + Ifges(4,6) * t363;
t126 = -mrSges(5,1) * t335 + mrSges(5,2) * t187;
t123 = -mrSges(7,2) * t214 - mrSges(7,3) * t331;
t122 = mrSges(6,1) * t214 + mrSges(6,3) * t330;
t121 = mrSges(7,1) * t214 + mrSges(7,3) * t330;
t120 = Ifges(5,1) * t187 + Ifges(5,4) * t335 - Ifges(5,5) * t397;
t119 = Ifges(5,4) * t187 + Ifges(5,2) * t335 - Ifges(5,6) * t397;
t116 = -mrSges(6,1) * t335 + mrSges(6,3) * t332;
t115 = -mrSges(7,1) * t335 + mrSges(7,3) * t332;
t114 = mrSges(6,2) * t335 + mrSges(6,3) * t162;
t113 = mrSges(7,2) * t335 + mrSges(7,3) * t162;
t111 = mrSges(6,1) * t331 - mrSges(6,2) * t330;
t102 = -mrSges(7,1) * t162 - mrSges(7,2) * t332;
t100 = -mrSges(5,2) * t363 - mrSges(5,3) * t109;
t99 = mrSges(5,1) * t363 - mrSges(5,3) * t108;
t82 = -Ifges(6,5) * t361 - Ifges(6,6) * t331 + t386;
t81 = -Ifges(7,5) * t361 - Ifges(7,6) * t331 + t387;
t70 = -Ifges(6,5) * t332 + Ifges(6,6) * t162 - Ifges(6,3) * t335;
t69 = -Ifges(7,5) * t332 + Ifges(7,6) * t162 - Ifges(7,3) * t335;
t44 = mrSges(5,1) * t109 + mrSges(5,2) * t108;
t43 = Ifges(5,1) * t108 - Ifges(5,4) * t109 + t350;
t42 = Ifges(5,4) * t108 - Ifges(5,2) * t109 + t349;
t34 = -mrSges(6,2) * t109 + mrSges(6,3) * t60;
t33 = -mrSges(7,2) * t109 + mrSges(7,3) * t60;
t32 = mrSges(6,1) * t109 - mrSges(6,3) * t59;
t31 = mrSges(7,1) * t109 - mrSges(7,3) * t59;
t30 = pkin(5) * t214 + t334 * t323 + (-t224 + (qJ(6) * t266 - t199) * t319) * qJD(5) + t352;
t28 = -mrSges(6,1) * t60 + mrSges(6,2) * t59;
t1 = pkin(5) * t109 - qJ(6) * t59 + qJD(6) * t332 + t6;
t2 = [t470 * t59 + t471 * t60 - t472 * t332 + t473 * t162 - (-t42 + t458) * t335 + (mrSges(3,3) * t322 * t446 + (0.2e1 * t248 * mrSges(3,3) - t366 - t454) * t326 + ((t259 * t450 + Ifges(3,5) * t318 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t326) * t376) * t326 + (t260 * t450 + Ifges(4,5) * t253 + Ifges(5,5) * t187 - 0.2e1 * Ifges(3,6) * t318 + Ifges(4,6) * t252 + Ifges(5,6) * t335 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t322) * t376 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3)) * t397) * t322) * qJD(2)) * t317 + (-t119 + t69 + t70) * t109 + (t1 * t26 + t29 * t3 + t49 * t8) * t451 + (t22 * t67 + t36 * t6 + t37 * t5) * t452 + (t172 * t192 + t24 * t80 + t25 * t79) * t453 + (-mrSges(4,1) * t252 + mrSges(4,2) * t253) * t446 + t252 * t133 + t253 * t134 + 0.2e1 * t242 * t150 + t222 * t179 + t223 * t180 + 0.2e1 * t117 * t231 + 0.2e1 * t118 * t232 + t187 * t43 + 0.2e1 * t192 * t44 + 0.2e1 * t174 * t184 + 0.2e1 * t173 * t185 + 0.2e1 * t24 * t170 + 0.2e1 * t25 * t171 + 0.2e1 * t172 * t126 + 0.2e1 * t3 * t113 + 0.2e1 * t5 * t114 + 0.2e1 * t1 * t115 + 0.2e1 * t6 * t116 + t108 * t120 + 0.2e1 * t80 * t100 + 0.2e1 * t8 * t102 + 0.2e1 * t22 * t103 + 0.2e1 * t79 * t99 + 0.2e1 * t67 * t28 + 0.2e1 * t49 * t27 + 0.2e1 * t36 * t32 + 0.2e1 * t37 * t34 + 0.2e1 * t26 * t31 + 0.2e1 * t29 * t33 + (t298 - 0.2e1 * t413 - 0.2e1 * t414) * t318 + 0.2e1 * m(3) * (t248 * t260 - t249 * t259) + 0.2e1 * m(4) * (t117 * t174 + t118 * t173 + t242 * t249); -(t28 - t99) * t457 + m(6) * (t140 * t6 + t141 * t5 + t154 * t67 - t22 * t457 + t36 * t48 + t37 * t47) + m(5) * (t153 * t80 - t154 * t79 + t172 * t307 + t192 * t373 + t234 * t24 + t25 * t457) + ((t208 / 0.2e1 + t207 / 0.2e1 - t312 / 0.2e1) * t326 + (Ifges(4,5) * t321 / 0.2e1 + Ifges(4,6) * t325 / 0.2e1 - Ifges(3,6)) * t383) * t317 + (-pkin(2) * t249 + (t117 * t325 - t118 * t321 + (-t173 * t325 - t174 * t321) * qJD(3)) * pkin(9)) * m(4) + (-t249 * mrSges(4,1) + t133 / 0.2e1 + pkin(9) * t184 + t117 * mrSges(4,3)) * t325 + (t249 * mrSges(4,2) + t134 / 0.2e1 - pkin(9) * t185 - t118 * mrSges(4,3)) * t321 - (t85 / 0.2e1 + t86 / 0.2e1) * t332 - (t82 / 0.2e1 + t81 / 0.2e1 - t147 / 0.2e1) * t335 + (-t220 / 0.2e1 + t164 / 0.2e1 + t165 / 0.2e1) * t109 + t298 + m(7) * (t1 * t112 + t101 * t49 + t125 * t3 + t181 * t8 + t26 * t30 + t29 * t38) - t414 - t413 + (t167 / 0.2e1 + t166 / 0.2e1) * t60 + (t168 / 0.2e1 + t169 / 0.2e1) * t59 + t307 * t44 + t222 * t290 / 0.2e1 + t223 * t293 / 0.2e1 + t253 * t278 / 0.2e1 + t242 * t270 + t252 * t275 / 0.2e1 + t234 * t100 + t172 * t219 + t108 * t221 / 0.2e1 + t1 * t203 + t6 * t204 + t187 * t148 / 0.2e1 + t192 * t146 + t8 * t193 + t22 * t194 + t3 * t201 + t5 * t202 + t153 * t170 + t181 * t27 + t141 * t34 - pkin(2) * t150 + t140 * t32 + t26 * t121 + t36 * t122 + t29 * t123 + t37 * t124 + t125 * t33 + t112 * t31 + t38 * t113 + t47 * t114 + t30 * t115 + t48 * t116 + t49 * t110 + t67 * t111 + t101 * t102 + t390 * t154 + (t43 / 0.2e1 - t25 * mrSges(5,3) + t350 / 0.2e1 + (t17 / 0.2e1 + t18 / 0.2e1) * t323 + (-t15 / 0.2e1 - t16 / 0.2e1) * t319 + (-t319 * t368 + t323 * t369) * qJD(5)) * t266 - (t120 / 0.2e1 - t79 * mrSges(5,3) + t368 * t323 + t369 * t319) * t213 + ((t180 / 0.2e1 - pkin(9) * t232 - t173 * mrSges(4,3)) * t325 + (t370 / 0.2e1 - t179 / 0.2e1 - pkin(9) * t231 - t174 * mrSges(4,3) + pkin(3) * t126) * t321) * qJD(3) + (t13 / 0.2e1 + t14 / 0.2e1 - t42 / 0.2e1 - t24 * mrSges(5,3) - t349 / 0.2e1) * t265 + (t69 / 0.2e1 + t70 / 0.2e1 - t119 / 0.2e1 - t80 * mrSges(5,3)) * t214 + (t83 / 0.2e1 + t84 / 0.2e1) * t162; t325 * t275 + t321 * t278 + 0.2e1 * t307 * t146 - 0.2e1 * pkin(2) * t270 + t111 * t447 + 0.2e1 * t30 * t203 + 0.2e1 * t48 * t204 + 0.2e1 * t101 * t193 + t194 * t448 + 0.2e1 * t38 * t201 + 0.2e1 * t47 * t202 + 0.2e1 * t181 * t110 + 0.2e1 * t141 * t124 + 0.2e1 * t140 * t122 + 0.2e1 * t112 * t121 + 0.2e1 * t125 * t123 + (t325 * t293 + (0.2e1 * pkin(3) * t219 - t290) * t321) * qJD(3) + (t153 * t234 + t307 * t373 - t408) * t453 + (t140 * t48 + t141 * t47 - t408) * t452 + (t101 * t181 + t112 * t30 + t125 * t38) * t451 + (t153 * t449 - t147 + t81 + t82) * t265 + (t234 * t449 + t164 + t165 - t220) * t214 - (mrSges(5,3) * t447 - t319 * t467 + t323 * t388 + t221) * t213 + (mrSges(5,3) * t448 + t148 + t468 * t323 - t469 * t319 + (-t319 * t388 - t323 * t467) * qJD(5)) * t266; (m(5) * (t24 * t320 + t25 * t324) + t324 * t99 + t320 * t100 + (t390 * t320 + (t114 * t323 - t116 * t319 + t170) * t324 + m(6) * (t320 * t67 - t36 * t395 + t37 * t392) + m(5) * (-t320 * t79 + t324 * t80)) * qJD(4)) * pkin(3) + (t339 * mrSges(7,3) + t338 * mrSges(6,3) + (m(6) * t338 - t319 * t114 - t323 * t116) * t304) * qJD(5) + (-mrSges(6,3) * t6 - mrSges(7,3) * t1 - t304 * t32) * t319 + m(6) * (t22 * t305 + t5 * t401 - t6 * t402) + m(7) * (t1 * t261 + t211 * t29 + t212 * t26 + t262 * t3 + t279 * t49 + t281 * t8) + t328 + t305 * t28 + t279 * t102 + t281 * t27 + t261 * t31 + t262 * t33 + t212 * t115 + t211 * t113 - t117 * mrSges(4,2) + t118 * mrSges(4,1) + t34 * t401 + t454; (-Ifges(4,6) * t321 + (-mrSges(4,1) * t325 + mrSges(4,2) * t321) * pkin(9)) * qJD(3) + (-t48 * mrSges(6,3) - t30 * mrSges(7,3) - t304 * t122) * t319 + m(6) * (t154 * t305 + t47 * t401 - t48 * t402) + (m(5) * (t153 * t320 - t154 * t324) + (t213 * t324 - t214 * t320) * mrSges(5,3) + ((mrSges(5,3) * t266 + t194) * t320 + (-t265 * mrSges(5,3) + t202 * t323 - t204 * t319) * t324 + m(6) * (-t140 * t395 + t141 * t392 - t406) + m(5) * (t234 * t324 - t406)) * qJD(4)) * pkin(3) + (t337 * mrSges(7,3) + t336 * mrSges(6,3) + (m(6) * t336 - t319 * t202 - t323 * t204) * t304) * qJD(5) + t312 + m(7) * (t101 * t281 + t112 * t212 + t125 * t211 + t181 * t279 + t261 * t30 + t262 * t38) + t305 * t111 + t279 * t193 + t281 * t110 + t261 * t121 + t262 * t123 + t212 * t203 + t211 * t201 + t304 * t393 + t327; t279 * t444 + t281 * t445 + 0.2e1 * t305 * t269 + (t211 * t262 + t212 * t261 + t279 * t281) * t451 + (t211 * t442 + t212 * t443 + (-t261 * t323 - t262 * t319) * t375) * mrSges(7,3) + (-0.2e1 * t410 - 0.2e1 * t412 + (t304 * t456 + t305 * t320) * t452 + 0.2e1 * t394 + 0.2e1 * t345) * t415 + t329; (m(6) * (-t36 * t377 - t37 * t378 - t425 + t426) + t323 * t34 - t319 * t32 - t116 * t377 - t114 * t378) * pkin(11) + t328 + t306 * t27 + t282 * t31 + t285 * t33 + t250 * t113 + t251 * t115 - pkin(4) * t28 + (qJD(5) * t338 - t425) * mrSges(6,3) + (qJD(5) * t339 - t1 * t319) * mrSges(7,3) + t102 * t371 + m(7) * (t1 * t282 + t250 * t29 + t251 * t26 + t285 * t3 + t306 * t8 + t371 * t49) - t22 * t441; (m(6) * (-t140 * t377 - t141 * t378 - t409 + t411) + t393 - t319 * t122 - t204 * t377 - t202 * t378) * pkin(11) + (qJD(5) * t336 - t409) * mrSges(6,3) + (qJD(5) * t337 - t30 * t319) * mrSges(7,3) + m(7) * (t101 * t306 + t112 * t251 + t125 * t250 + t181 * t371 + t282 * t30 + t285 * t38) + t306 * t110 + t282 * t121 + t285 * t123 + t250 * t201 + t251 * t203 - pkin(4) * t111 - t154 * t441 + t193 * t371 + t327; m(7) * (t211 * t285 + t212 * t282 + t250 * t262 + t251 * t261 + t279 * t306 + t281 * t371) + (t279 + t371) * t283 + (t305 - pkin(4)) * t269 + (t306 + t281) * t268 + (m(6) * (-pkin(4) * t320 + pkin(11) * t456) - t410 + t394 - t412 + t345) * t415 + ((t211 + t250) * t323 + (-t212 - t251) * t319 + ((-t261 - t282) * t323 + (-t262 - t285) * t319) * qJD(5)) * mrSges(7,3) + t329; -0.2e1 * pkin(4) * t269 + t371 * t444 + t306 * t445 + (t250 * t285 + t251 * t282 + t306 * t371) * t451 + (t250 * t442 + t251 * t443 + (-t282 * t323 - t285 * t319) * t375) * mrSges(7,3) + t329; mrSges(6,1) * t6 + mrSges(7,1) * t1 - mrSges(6,2) * t5 - mrSges(7,2) * t3 + (m(7) * t1 + t31) * pkin(5) + t458; mrSges(6,1) * t48 + mrSges(7,1) * t30 - mrSges(6,2) * t47 - mrSges(7,2) * t38 + t474 * t396 + (m(7) * t30 + t121) * pkin(5) + (-t319 * t475 - t323 * t474) * t266 * qJD(5) + t386 + t387; -mrSges(7,2) * t211 + t367 * t212 - t344 * t372 + ((-mrSges(6,1) * t304 - t437) * t323 + (mrSges(6,2) * t304 - t474) * t319) * qJD(5) + t385; -mrSges(7,2) * t250 + t367 * t251 + ((-mrSges(6,1) * pkin(11) - t437) * t323 + (mrSges(6,2) * pkin(11) - t474) * t319) * qJD(5) + t385; 0; m(7) * t8 + t27; m(7) * t101 + t110; m(7) * t279 + t268; m(7) * t371 + t268; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
