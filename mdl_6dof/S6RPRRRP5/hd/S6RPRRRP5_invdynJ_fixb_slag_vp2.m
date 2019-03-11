% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:46
% EndTime: 2019-03-09 06:11:26
% DurationCPUTime: 24.34s
% Computational Cost: add. (15728->700), mult. (38629->885), div. (0->0), fcn. (30126->14), ass. (0->310)
t523 = -mrSges(6,3) - mrSges(7,2);
t518 = Ifges(6,1) + Ifges(7,1);
t527 = -Ifges(6,4) + Ifges(7,5);
t494 = Ifges(6,5) + Ifges(7,4);
t517 = Ifges(6,6) - Ifges(7,6);
t516 = Ifges(6,3) + Ifges(7,2);
t270 = sin(pkin(10));
t271 = cos(pkin(10));
t374 = t270 ^ 2 + t271 ^ 2;
t366 = qJD(1) * qJD(2);
t244 = qJ(2) * qJDD(1) + t366;
t277 = cos(qJ(5));
t368 = qJD(5) * t277;
t275 = sin(qJ(3));
t279 = cos(qJ(3));
t384 = t271 * t279;
t216 = -t270 * t275 + t384;
t206 = t216 * qJD(1);
t217 = t270 * t279 + t271 * t275;
t207 = t217 * qJD(1);
t274 = sin(qJ(4));
t278 = cos(qJ(4));
t328 = t278 * t206 - t207 * t274;
t482 = t328 * t277;
t526 = t368 - t482;
t273 = sin(qJ(5));
t369 = qJD(5) * t273;
t525 = t328 * t273 - t369;
t208 = t216 * qJD(3);
t176 = qJD(1) * t208 + qJDD(1) * t217;
t209 = t217 * qJD(3);
t177 = -qJD(1) * t209 + qJDD(1) * t216;
t101 = qJD(4) * t328 + t176 * t278 + t177 * t274;
t265 = qJDD(3) + qJDD(4);
t299 = t206 * t274 + t278 * t207;
t363 = qJD(3) + qJD(4);
t157 = t273 * t299 - t277 * t363;
t370 = qJD(5) * t157;
t71 = t277 * t101 + t273 * t265 - t370;
t464 = t71 / 0.2e1;
t158 = t273 * t363 + t277 * t299;
t72 = qJD(5) * t158 + t273 * t101 - t277 * t265;
t462 = t72 / 0.2e1;
t102 = -qJD(4) * t299 - t176 * t274 + t177 * t278;
t97 = qJDD(5) - t102;
t460 = t97 / 0.2e1;
t520 = m(6) + m(7);
t524 = m(5) + t520;
t496 = mrSges(6,1) + mrSges(7,1);
t495 = mrSges(6,2) - mrSges(7,3);
t168 = qJD(5) - t328;
t425 = mrSges(7,2) * t157;
t121 = mrSges(7,3) * t168 - t425;
t422 = mrSges(6,3) * t157;
t122 = -mrSges(6,2) * t168 - t422;
t480 = t121 + t122;
t421 = mrSges(6,3) * t158;
t123 = mrSges(6,1) * t168 - t421;
t424 = mrSges(7,2) * t158;
t124 = -mrSges(7,1) * t168 + t424;
t479 = t123 - t124;
t269 = pkin(10) + qJ(3);
t264 = qJ(4) + t269;
t252 = sin(t264);
t522 = t523 * t252;
t427 = pkin(7) + qJ(2);
t239 = t427 * t270;
t220 = qJD(1) * t239;
t240 = t427 * t271;
t221 = qJD(1) * t240;
t181 = -t220 * t275 + t221 * t279;
t156 = pkin(8) * t206 + t181;
t150 = t274 * t156;
t394 = t221 * t275;
t180 = -t279 * t220 - t394;
t155 = -pkin(8) * t207 + t180;
t152 = qJD(3) * pkin(3) + t155;
t109 = t278 * t152 - t150;
t103 = -pkin(4) * t363 - t109;
t410 = t299 * mrSges(5,3);
t478 = mrSges(5,1) * t363 - mrSges(6,1) * t157 - mrSges(6,2) * t158 - t410;
t521 = -m(6) * t103 + t478;
t515 = t494 * t97 + t518 * t71 + t527 * t72;
t491 = -t157 * t517 + t158 * t494 + t516 * t168;
t154 = Ifges(6,4) * t157;
t412 = t157 * Ifges(7,5);
t490 = t158 * t518 + t494 * t168 - t154 + t412;
t151 = t278 * t156;
t112 = t155 * t274 + t151;
t372 = qJD(4) * t274;
t514 = -pkin(3) * t372 + t112;
t306 = pkin(5) * t273 - qJ(6) * t277;
t513 = pkin(5) * t369 - qJ(6) * t368 - qJD(6) * t273 - t306 * t328;
t307 = pkin(5) * t277 + qJ(6) * t273;
t233 = -pkin(4) - t307;
t315 = -t277 * mrSges(7,1) - t273 * mrSges(7,3);
t426 = mrSges(6,1) * t277;
t474 = (-m(7) * t233 - t315 + t426) * t252;
t314 = t273 * mrSges(7,1) - t277 * mrSges(7,3);
t316 = mrSges(6,1) * t273 + mrSges(6,2) * t277;
t53 = t157 * pkin(5) - t158 * qJ(6) + t103;
t512 = t103 * t316 + t53 * t314;
t511 = -t273 * t517 + t277 * t494;
t416 = Ifges(7,5) * t273;
t418 = Ifges(6,4) * t273;
t510 = t277 * t518 + t416 - t418;
t298 = t278 * t216 - t217 * t274;
t136 = qJD(4) * t298 + t208 * t278 - t209 * t274;
t179 = t216 * t274 + t217 * t278;
t293 = t136 * t273 + t179 * t368;
t509 = t494 * t71 + t516 * t97 - t517 * t72;
t110 = t274 * t152 + t151;
t104 = pkin(9) * t363 + t110;
t254 = pkin(2) * t271 + pkin(1);
t224 = -qJD(1) * t254 + qJD(2);
t182 = -pkin(3) * t206 + t224;
t111 = -pkin(4) * t328 - pkin(9) * t299 + t182;
t329 = pkin(7) * qJDD(1) + t244;
t199 = t329 * t270;
t200 = t329 * t271;
t135 = -qJD(3) * t181 - t279 * t199 - t200 * t275;
t100 = qJDD(3) * pkin(3) - pkin(8) * t176 + t135;
t373 = qJD(3) * t279;
t134 = -qJD(3) * t394 - t275 * t199 + t279 * t200 - t220 * t373;
t107 = pkin(8) * t177 + t134;
t371 = qJD(4) * t278;
t30 = t274 * t100 + t278 * t107 + t152 * t371 - t156 * t372;
t27 = pkin(9) * t265 + t30;
t223 = -qJDD(1) * t254 + qJDD(2);
t159 = -pkin(3) * t177 + t223;
t37 = -pkin(4) * t102 - pkin(9) * t101 + t159;
t6 = -t104 * t369 + t111 * t368 + t277 * t27 + t273 * t37;
t46 = t104 * t277 + t111 * t273;
t7 = -qJD(5) * t46 - t27 * t273 + t277 * t37;
t508 = -t273 * t7 + t277 * t6;
t2 = qJ(6) * t97 + qJD(6) * t168 + t6;
t4 = -pkin(5) * t97 + qJDD(6) - t7;
t507 = t2 * t277 + t273 * t4;
t506 = Ifges(7,5) * t464 + Ifges(7,6) * t460 - Ifges(6,4) * t71 / 0.2e1 - Ifges(6,6) * t97 / 0.2e1 + (Ifges(7,3) + Ifges(6,2)) * t462;
t45 = -t104 * t273 + t111 * t277;
t39 = -pkin(5) * t168 + qJD(6) - t45;
t40 = qJ(6) * t168 + t46;
t33 = mrSges(6,1) * t97 - mrSges(6,3) * t71;
t34 = -t97 * mrSges(7,1) + t71 * mrSges(7,2);
t492 = -t33 + t34;
t32 = -mrSges(7,2) * t72 + mrSges(7,3) * t97;
t35 = -mrSges(6,2) * t97 - mrSges(6,3) * t72;
t493 = t32 + t35;
t504 = m(6) * ((-t273 * t46 - t277 * t45) * qJD(5) + t508) + m(7) * ((-t273 * t40 + t277 * t39) * qJD(5) + t507) - t480 * t369 - t479 * t368 + t277 * t493 + t273 * t492;
t341 = m(3) * qJ(2) + mrSges(3,3);
t503 = -m(4) * t427 + mrSges(2,2) - mrSges(4,3) - mrSges(5,3) - t341;
t502 = m(7) * pkin(5) + t496;
t262 = sin(t269);
t263 = cos(t269);
t320 = t263 * mrSges(4,1) - t262 * mrSges(4,2);
t321 = -mrSges(3,1) * t271 + mrSges(3,2) * t270;
t253 = cos(t264);
t476 = t253 * mrSges(5,1) - t252 * mrSges(5,2);
t501 = -m(3) * pkin(1) - m(4) * t254 - mrSges(2,1) - t320 + t321 - t476;
t500 = m(7) * qJ(6) - t495;
t499 = -t7 * mrSges(6,1) + t4 * mrSges(7,1) + t6 * mrSges(6,2) - t2 * mrSges(7,3);
t31 = t100 * t278 - t274 * t107 - t152 * t372 - t156 * t371;
t28 = -pkin(4) * t265 - t31;
t415 = Ifges(7,5) * t277;
t308 = Ifges(7,3) * t273 + t415;
t417 = Ifges(6,4) * t277;
t311 = -Ifges(6,2) * t273 + t417;
t336 = t368 / 0.2e1;
t337 = -t369 / 0.2e1;
t153 = Ifges(7,5) * t158;
t80 = Ifges(7,6) * t168 + Ifges(7,3) * t157 + t153;
t405 = t273 * t80;
t355 = t405 / 0.2e1;
t445 = t273 / 0.2e1;
t463 = -t72 / 0.2e1;
t419 = Ifges(6,4) * t158;
t83 = -Ifges(6,2) * t157 + Ifges(6,6) * t168 + t419;
t9 = pkin(5) * t72 - qJ(6) * t71 - qJD(6) * t158 + t28;
t498 = t31 * mrSges(5,1) - t30 * mrSges(5,2) + Ifges(5,5) * t101 + Ifges(5,6) * t102 + Ifges(5,3) * t265 - t28 * t426 + t9 * t315 + t83 * t337 + t416 * t462 + t418 * t463 + (-t415 + t417) * t464 + t515 * t445 + (-t311 / 0.2e1 + t308 / 0.2e1) * t370 + (t28 * mrSges(6,2) + t460 * t494 + t464 * t518) * t273 + (Ifges(6,2) * t463 - Ifges(7,3) * t462 + t460 * t517 - t506) * t277 + (t355 + t512) * qJD(5) + (t158 * t510 + t168 * t511) * qJD(5) / 0.2e1 + (t336 - t482 / 0.2e1) * t490 + (-t45 * t526 + t46 * t525 + t508) * mrSges(6,3) + (t39 * t526 + t40 * t525 + t507) * mrSges(7,2);
t450 = -t299 / 0.2e1;
t497 = pkin(5) * t299;
t489 = Ifges(5,5) * t363;
t488 = Ifges(5,6) * t363;
t487 = -t110 + t513;
t486 = t513 - t514;
t485 = qJ(6) * t299;
t183 = -t279 * t239 - t240 * t275;
t165 = -pkin(8) * t217 + t183;
t184 = -t275 * t239 + t279 * t240;
t166 = pkin(8) * t216 + t184;
t127 = t165 * t274 + t166 * t278;
t191 = -pkin(3) * t216 - t254;
t128 = -pkin(4) * t298 - pkin(9) * t179 + t191;
t481 = t277 * t127 + t273 * t128;
t477 = t278 * t165 - t166 * t274;
t375 = t253 * pkin(4) + t252 * pkin(9);
t438 = pkin(4) * t252;
t441 = pkin(3) * t262;
t475 = m(7) * t441 - m(6) * (-t438 - t441) + t474;
t280 = cos(qJ(1));
t388 = t253 * t280;
t237 = pkin(9) * t388;
t393 = t252 * t273;
t360 = mrSges(6,2) * t393;
t473 = -m(7) * t237 - t280 * t360 + t388 * t523;
t276 = sin(qJ(1));
t390 = t253 * t276;
t234 = pkin(9) * t390;
t472 = -m(7) * t234 - t276 * t360 + t390 * t523;
t471 = g(1) * t280 + g(2) * t276;
t389 = t253 * t277;
t391 = t253 * t273;
t470 = -t389 * t496 + t391 * t495 - t476 + t522;
t162 = -t239 * t373 + qJD(2) * t384 + (-qJD(2) * t270 - qJD(3) * t240) * t275;
t143 = -pkin(8) * t209 + t162;
t163 = -t217 * qJD(2) - qJD(3) * t184;
t144 = -pkin(8) * t208 + t163;
t49 = qJD(4) * t477 + t143 * t278 + t144 * t274;
t137 = qJD(4) * t179 + t208 * t274 + t278 * t209;
t442 = pkin(3) * t209;
t74 = pkin(4) * t137 - pkin(9) * t136 + t442;
t13 = -qJD(5) * t481 - t273 * t49 + t277 * t74;
t133 = pkin(4) * t299 - pkin(9) * t328;
t454 = -t168 / 0.2e1;
t456 = -t158 / 0.2e1;
t457 = t157 / 0.2e1;
t458 = -t157 / 0.2e1;
t466 = -t182 * mrSges(5,2) + Ifges(5,1) * t450 + t311 * t457 + t308 * t458 + t83 * t445 - t405 / 0.2e1 - t489 / 0.2e1 + t510 * t456 + t511 * t454 - t512;
t451 = -t328 / 0.2e1;
t465 = -t45 * mrSges(6,1) + t39 * mrSges(7,1) + t46 * mrSges(6,2) - t40 * mrSges(7,3) - t182 * mrSges(5,1) - Ifges(5,2) * t451 + Ifges(6,6) * t457 + Ifges(7,6) * t458 + t488 / 0.2e1 + t494 * t456 + t516 * t454;
t455 = t158 / 0.2e1;
t449 = t299 / 0.2e1;
t446 = t207 / 0.2e1;
t443 = pkin(3) * t207;
t250 = pkin(3) * t263;
t440 = pkin(3) * t274;
t439 = pkin(3) * t278;
t433 = g(3) * t252;
t420 = Ifges(5,4) * t299;
t414 = t109 * mrSges(5,3);
t413 = t110 * mrSges(5,3);
t411 = t328 * mrSges(5,3);
t409 = t180 * mrSges(4,3);
t408 = t181 * mrSges(4,3);
t407 = t207 * Ifges(4,4);
t403 = t136 * t277;
t392 = t252 * t280;
t383 = t273 * t276;
t382 = t276 * t277;
t381 = t277 * t280;
t380 = t280 * t273;
t57 = t277 * t109 + t273 * t133;
t113 = t155 * t278 - t150;
t118 = t133 + t443;
t55 = t277 * t113 + t273 * t118;
t365 = qJDD(1) * t270;
t364 = qJDD(1) * t271;
t358 = pkin(3) * t371;
t333 = -t102 * mrSges(5,1) + t101 * mrSges(5,2);
t332 = -t177 * mrSges(4,1) + t176 * mrSges(4,2);
t222 = t250 + t254;
t266 = -pkin(8) - t427;
t327 = t280 * t222 - t266 * t276;
t326 = t273 * t358;
t325 = t277 * t358;
t324 = pkin(5) * t389 + qJ(6) * t391 + t375;
t322 = -mrSges(3,1) * t364 + mrSges(3,2) * t365;
t317 = mrSges(5,1) * t252 + mrSges(5,2) * t253;
t304 = t273 * t39 + t277 * t40;
t303 = -t273 * t45 + t277 * t46;
t56 = -t109 * t273 + t133 * t277;
t54 = -t113 * t273 + t118 * t277;
t64 = -t127 * t273 + t128 * t277;
t292 = t179 * t369 - t403;
t12 = -t127 * t369 + t128 * t368 + t273 * t74 + t277 * t49;
t50 = qJD(4) * t127 + t143 * t274 - t278 * t144;
t261 = -qJDD(1) * pkin(1) + qJDD(2);
t258 = -pkin(4) - t439;
t212 = t233 - t439;
t201 = Ifges(4,4) * t206;
t198 = t253 * t381 + t383;
t197 = t253 * t380 - t382;
t196 = t253 * t382 - t380;
t195 = t253 * t383 + t381;
t187 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t207;
t186 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t206;
t170 = t207 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t201;
t169 = t206 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t407;
t167 = Ifges(5,4) * t328;
t160 = -mrSges(5,2) * t363 + t411;
t132 = -mrSges(5,1) * t328 + mrSges(5,2) * t299;
t130 = Ifges(5,1) * t299 + t167 + t489;
t129 = Ifges(5,2) * t328 + t420 + t488;
t115 = mrSges(7,1) * t157 - mrSges(7,3) * t158;
t114 = pkin(5) * t158 + qJ(6) * t157;
t90 = -mrSges(5,2) * t265 + mrSges(5,3) * t102;
t89 = mrSges(5,1) * t265 - mrSges(5,3) * t101;
t75 = t179 * t306 - t477;
t52 = pkin(5) * t298 - t64;
t51 = -qJ(6) * t298 + t481;
t44 = -t56 - t497;
t43 = t57 + t485;
t42 = -t54 - t497;
t41 = t55 + t485;
t25 = mrSges(6,1) * t72 + mrSges(6,2) * t71;
t24 = mrSges(7,1) * t72 - mrSges(7,3) * t71;
t14 = t306 * t136 + (qJD(5) * t307 - qJD(6) * t277) * t179 + t50;
t11 = -pkin(5) * t137 - t13;
t10 = qJ(6) * t137 - qJD(6) * t298 + t12;
t1 = [-(-m(5) * t31 + m(6) * t28 + t25 - t89) * t477 + ((mrSges(7,2) * t4 - mrSges(6,3) * t7 + t515 / 0.2e1) * t277 + t490 * t337 + (-t2 * mrSges(7,2) - t6 * mrSges(6,3) + t506) * t273 + t159 * mrSges(5,2) - t31 * mrSges(5,3) + Ifges(5,1) * t101 + Ifges(5,4) * t102 + Ifges(5,5) * t265 + t28 * t316 + t308 * t462 + t311 * t463 + t314 * t9 + t336 * t80 + t460 * t511 + t464 * t510) * t179 + (-t509 / 0.2e1 - t159 * mrSges(5,1) + t30 * mrSges(5,3) + Ifges(5,4) * t101 + Ifges(5,2) * t102 + Ifges(5,6) * t265 - Ifges(6,6) * t463 - Ifges(7,6) * t462 - t460 * t516 - t464 * t494 + t499) * t298 + (Ifges(4,1) * t208 - Ifges(4,4) * t209) * t446 + qJD(3) * (Ifges(4,5) * t208 - Ifges(4,6) * t209) / 0.2e1 + t206 * (Ifges(4,4) * t208 - Ifges(4,2) * t209) / 0.2e1 + t224 * (mrSges(4,1) * t209 + mrSges(4,2) * t208) + t328 * (Ifges(5,4) * t136 - Ifges(5,2) * t137) / 0.2e1 + m(5) * (t110 * t49 + t127 * t30 + t159 * t191 + t182 * t442) + (Ifges(3,4) * t270 + Ifges(3,2) * t271) * t364 + (Ifges(3,1) * t270 + Ifges(3,4) * t271) * t365 - t137 * t413 - t136 * t414 + (mrSges(4,2) * t223 - mrSges(4,3) * t135 + Ifges(4,1) * t176 + Ifges(4,4) * t177 + Ifges(4,5) * qJDD(3)) * t217 + (t502 * t196 + t500 * t195 + (m(5) * t222 - t520 * (-t222 - t375) - t501 - t522) * t276 + (t266 * t524 + t503) * t280) * g(1) + (-mrSges(4,1) * t223 + mrSges(4,3) * t134 + Ifges(4,4) * t176 + Ifges(4,2) * t177 + Ifges(4,6) * qJDD(3)) * t216 + (-m(5) * t327 + t523 * t392 - t520 * (pkin(4) * t388 + pkin(9) * t392 + t327) - t502 * t198 - t500 * t197 + t501 * t280 + t503 * t276) * g(2) + t49 * t160 + m(6) * (t12 * t46 + t13 * t45 + t481 * t6 + t64 * t7) + t481 * t35 - t209 * t408 + t132 * t442 + (-m(5) * t109 - t521) * t50 + m(4) * (t134 * t184 + t135 * t183 + t162 * t181 + t163 * t180 - t223 * t254) + m(7) * (t10 * t40 + t11 * t39 + t14 * t53 + t2 * t51 + t4 * t52 + t75 * t9) - t254 * t332 + t191 * t333 - t137 * t129 / 0.2e1 + t136 * t130 / 0.2e1 + t127 * t90 + t10 * t121 + t12 * t122 + t13 * t123 + t11 * t124 + t14 * t115 + Ifges(2,3) * qJDD(1) + t75 * t24 + t64 * t33 + t51 * t32 + t52 * t34 - t208 * t409 - pkin(1) * t322 + t261 * t321 + (t494 * t137 - t518 * t292 + t293 * t527) * t455 + t103 * (mrSges(6,1) * t293 - mrSges(6,2) * t292) + t53 * (mrSges(7,1) * t293 + mrSges(7,3) * t292) + t40 * (-mrSges(7,2) * t293 + mrSges(7,3) * t137) + t46 * (-mrSges(6,2) * t137 - mrSges(6,3) * t293) + t45 * (mrSges(6,1) * t137 + mrSges(6,3) * t292) + t39 * (-mrSges(7,1) * t137 - mrSges(7,2) * t292) + (Ifges(5,1) * t136 - Ifges(5,4) * t137) * t449 - t293 * t83 / 0.2e1 + t182 * (mrSges(5,1) * t137 + mrSges(5,2) * t136) + t183 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t176) + t184 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t177) + t162 * t186 + t163 * t187 + t208 * t170 / 0.2e1 - t209 * t169 / 0.2e1 + t490 * t403 / 0.2e1 + t491 * t137 / 0.2e1 + (t137 * t516 - t292 * t494 - t293 * t517) * t168 / 0.2e1 + (-Ifges(7,5) * t292 + Ifges(7,6) * t137 + Ifges(7,3) * t293) * t457 + (-Ifges(6,4) * t292 - Ifges(6,2) * t293 + Ifges(6,6) * t137) * t458 + m(3) * (-pkin(1) * t261 + (t244 + t366) * qJ(2) * t374) + t136 * t355 + t363 * (Ifges(5,5) * t136 - Ifges(5,6) * t137) / 0.2e1 + 0.2e1 * t374 * t244 * mrSges(3,3); t332 + t333 + (t168 * t480 - t492) * t277 + (-t168 * t479 + t493) * t273 + t322 + (-t115 + t478) * t299 + m(3) * t261 - t328 * t160 - t206 * t186 + t207 * t187 + (-g(1) * t276 + g(2) * t280) * (m(3) + m(4) + t524) - t341 * t374 * qJD(1) ^ 2 + (t168 * t304 + t2 * t273 - t277 * t4 - t299 * t53) * m(7) + (-t103 * t299 + t168 * t303 + t273 * t6 + t277 * t7) * m(6) + (t109 * t299 - t110 * t328 + t159) * m(5) + (t180 * t207 - t181 * t206 + t223) * m(4); t486 * t115 + (t212 * t9 + t304 * t358 - t39 * t42 - t40 * t41 + t486 * t53) * m(7) + ((t274 * t30 + t278 * t31 + (-t109 * t274 + t110 * t278) * qJD(4)) * pkin(3) + t109 * t112 - t110 * t113 - t182 * t443) * m(5) + (-t41 + t325) * t121 + (-Ifges(5,4) * t450 + t465 + t129 / 0.2e1 + t413) * t299 + (Ifges(5,4) * t451 + t466 - t130 / 0.2e1 + t414) * t328 + (t276 * t475 + t472) * g(2) + (t280 * t475 + t473) * g(1) + (-m(5) * t250 - m(7) * (t250 + t324) - t320 - m(6) * (t250 + t375) + t470) * g(3) + (-t113 + t358) * t160 + (-t42 + t326) * t124 + (-t326 - t54) * t123 + (-t55 + t325) * t122 - (-Ifges(4,2) * t207 + t170 + t201) * t206 / 0.2e1 + t498 + (m(5) * t441 + mrSges(4,1) * t262 + mrSges(4,2) * t263 + t317) * t471 + (t258 * t28 + (t103 * t274 + t278 * t303) * qJD(4) * pkin(3) - g(1) * t237 - g(2) * t234 - t103 * t112 - t45 * t54 - t46 * t55) * m(6) + t207 * t408 + t206 * t409 - t132 * t443 - t134 * mrSges(4,2) + t135 * mrSges(4,1) - t207 * (Ifges(4,1) * t206 - t407) / 0.2e1 + t491 * t450 + Ifges(4,3) * qJDD(3) + t89 * t439 + t90 * t440 + t169 * t446 + t504 * (pkin(9) + t440) + Ifges(4,5) * t176 + Ifges(4,6) * t177 + t514 * t478 - t180 * t186 + t181 * t187 - qJD(3) * (Ifges(4,5) * t206 - Ifges(4,6) * t207) / 0.2e1 + t212 * t24 - t224 * (mrSges(4,1) * t207 + mrSges(4,2) * t206) + t258 * t25; (t167 + t130) * t451 + (t410 + t521) * t110 + t487 * t115 + (t233 * t9 - t39 * t44 - t40 * t43 + t487 * t53) * m(7) + t465 * t299 + t466 * t328 + t471 * t317 + (t276 * t474 + t472) * g(2) + (t280 * t474 + t473) * g(1) + (-m(6) * t375 - m(7) * t324 + t470) * g(3) + t504 * pkin(9) + (-t45 * t56 - t46 * t57 - g(1) * (-pkin(4) * t392 + t237) - pkin(4) * t28 - g(2) * (-t276 * t438 + t234)) * m(6) + t498 + (-t420 + t491) * t450 - t43 * t121 - t57 * t122 - t56 * t123 - t44 * t124 - pkin(4) * t25 + t129 * t449 + (-t160 + t411) * t109 + t233 * t24; (-m(7) * t39 + t421 + t479) * t46 + (-m(7) * t40 - t422 - t480) * t45 + (-t114 * t53 + t306 * t433 - pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t40 - g(2) * (-pkin(5) * t195 + qJ(6) * t196) - g(1) * (-pkin(5) * t197 + qJ(6) * t198)) * m(7) + (t195 * t496 + t196 * t495) * g(2) + (t197 * t496 + t198 * t495) * g(1) + t40 * t424 + t39 * t425 - t53 * (mrSges(7,1) * t158 + mrSges(7,3) * t157) - t103 * (mrSges(6,1) * t158 - mrSges(6,2) * t157) - t499 + (t316 + t314) * t433 + t509 + (-Ifges(6,2) * t158 - t154 + t490) * t457 + (-t494 * t157 - t158 * t517) * t454 + (-t157 * t518 + t153 - t419 + t80) * t456 + qJD(6) * t121 - t114 * t115 + qJ(6) * t32 - pkin(5) * t34 + t83 * t455 + (Ifges(7,3) * t158 - t412) * t458; t158 * t115 - t168 * t121 + (-g(1) * t197 - g(2) * t195 - g(3) * t393 + t53 * t158 - t40 * t168 + t4) * m(7) + t34;];
tau  = t1;
