% Calculate vector of inverse dynamics joint torques for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:19
% EndTime: 2019-12-31 21:44:11
% DurationCPUTime: 29.02s
% Computational Cost: add. (7540->778), mult. (18572->1030), div. (0->0), fcn. (14053->10), ass. (0->364)
t252 = cos(qJ(2));
t245 = sin(pkin(5));
t372 = qJD(1) * t245;
t346 = t252 * t372;
t531 = qJD(3) - t346;
t247 = sin(qJ(3));
t251 = cos(qJ(3));
t248 = sin(qJ(2));
t361 = qJD(1) * qJD(2);
t272 = qJDD(1) * t248 + t252 * t361;
t260 = t272 * t245;
t393 = cos(pkin(5));
t326 = t393 * qJD(1);
t283 = t326 + qJD(2);
t267 = qJD(3) * t283;
t324 = t393 * qJDD(1);
t281 = t324 + qJDD(2);
t385 = t245 * t248;
t353 = t247 * t385;
t316 = qJD(3) * t353;
t87 = qJD(1) * t316 - t247 * t281 + (-t260 - t267) * t251;
t451 = -t87 / 0.2e1;
t530 = Ifges(4,1) * t451;
t371 = qJD(1) * t248;
t345 = t251 * t371;
t88 = t247 * t267 - t251 * t281 + (qJD(3) * t345 + t247 * t272) * t245;
t449 = -t88 / 0.2e1;
t529 = Ifges(4,4) * t449;
t315 = pkin(1) * t326;
t185 = pkin(7) * t346 + t248 * t315;
t366 = qJD(3) * t247;
t528 = -qJD(4) * t247 - t185 + (-t247 * t346 + t366) * pkin(3);
t188 = (-qJDD(1) * t252 + t248 * t361) * t245;
t178 = qJDD(3) + t188;
t434 = t178 / 0.2e1;
t504 = Ifges(5,4) - Ifges(4,5);
t503 = Ifges(4,6) - Ifges(5,5);
t502 = Ifges(4,3) + Ifges(5,1);
t527 = t528 + t531 * (pkin(9) * t247 - qJ(4) * t251);
t347 = t245 * t371;
t182 = -pkin(7) * t347 + t252 * t315;
t279 = (pkin(2) * t248 - pkin(8) * t252) * t245;
t183 = qJD(1) * t279;
t108 = -t247 * t182 + t183 * t251;
t446 = pkin(4) + pkin(8);
t222 = t446 * t251;
t375 = t251 * t252;
t447 = pkin(3) + pkin(9);
t526 = qJD(3) * t222 - (pkin(4) * t375 - t248 * t447) * t372 + t108;
t246 = sin(qJ(5));
t250 = cos(qJ(5));
t302 = mrSges(6,1) * t246 + mrSges(6,2) * t250;
t450 = t87 / 0.2e1;
t448 = t88 / 0.2e1;
t525 = -t178 / 0.2e1;
t157 = t247 * t347 - t251 * t283;
t112 = t157 * t250 - t246 * t531;
t36 = qJD(5) * t112 + t178 * t250 + t246 * t88;
t113 = t157 * t246 + t250 * t531;
t37 = -qJD(5) * t113 - t178 * t246 + t250 * t88;
t84 = qJDD(5) - t87;
t5 = Ifges(6,5) * t36 + Ifges(6,6) * t37 + Ifges(6,3) * t84;
t524 = t530 + t529 + Ifges(4,5) * t434 + t5 / 0.2e1;
t360 = qJDD(1) * t245;
t521 = pkin(7) * t360 + qJD(2) * t315;
t520 = -pkin(7) * t245 * t361 + pkin(1) * t324;
t513 = m(6) * pkin(9);
t479 = mrSges(4,1) - mrSges(5,2) + t513;
t512 = m(6) + m(5);
t519 = pkin(3) * t512 + t479;
t142 = -pkin(2) * t283 - t182;
t158 = t245 * t345 + t247 * t283;
t255 = -t158 * qJ(4) + t142;
t63 = t157 * pkin(3) + t255;
t518 = t142 * mrSges(4,2) - t63 * mrSges(5,3);
t143 = pkin(8) * t283 + t185;
t149 = (-pkin(2) * t252 - pkin(8) * t248 - pkin(1)) * t372;
t85 = t143 * t247 - t251 * t149;
t280 = pkin(4) * t158 + t85;
t517 = qJD(4) + t280;
t516 = Ifges(4,6) * t525 + 0.2e1 * Ifges(5,3) * t448 + (Ifges(4,4) + Ifges(5,6)) * t450 + (t448 - t449) * Ifges(4,2) + (-t503 + Ifges(5,5)) * t434;
t16 = mrSges(6,1) * t84 - mrSges(6,3) * t36;
t17 = -mrSges(6,2) * t84 + mrSges(6,3) * t37;
t152 = qJD(5) + t158;
t70 = -mrSges(6,2) * t152 + mrSges(6,3) * t112;
t71 = mrSges(6,1) * t152 - mrSges(6,3) * t113;
t285 = -t246 * t71 + t250 * t70;
t515 = t285 * qJD(5) + t250 * t16 + t246 * t17;
t44 = -t447 * t531 + t517;
t47 = t157 * t447 + t255;
t12 = -t246 * t47 + t250 * t44;
t13 = t246 * t44 + t250 * t47;
t476 = t12 * mrSges(6,1) - t13 * mrSges(6,2);
t151 = Ifges(4,4) * t157;
t498 = t158 * Ifges(4,1) + Ifges(4,5) * t531 + Ifges(6,5) * t113 + Ifges(6,6) * t112 + Ifges(6,3) * t152 - t151;
t67 = -pkin(3) * t531 + qJD(4) + t85;
t150 = Ifges(5,6) * t157;
t76 = Ifges(5,4) * t531 - t158 * Ifges(5,2) + t150;
t514 = -t76 / 0.2e1 + t476 + t498 / 0.2e1 + t67 * mrSges(5,1) + t85 * mrSges(4,3);
t511 = -t188 / 0.2e1;
t510 = t260 / 0.2e1;
t509 = t281 / 0.2e1;
t506 = -mrSges(5,1) - mrSges(4,3);
t505 = mrSges(4,2) - mrSges(5,3);
t391 = qJ(4) * t247;
t331 = -pkin(2) - t391;
t204 = -t251 * t447 + t331;
t221 = t446 * t247;
t121 = -t204 * t246 + t221 * t250;
t501 = qJD(5) * t121 + t246 * t526 + t250 * t527;
t122 = t204 * t250 + t221 * t246;
t500 = -qJD(5) * t122 - t246 * t527 + t250 * t526;
t499 = -t157 * t503 - t158 * t504 + t502 * t531;
t419 = mrSges(5,1) * t157;
t118 = -mrSges(5,3) * t531 + t419;
t51 = -mrSges(6,1) * t112 + mrSges(6,2) * t113;
t496 = -t118 + t51;
t109 = t251 * t182 + t247 * t183;
t378 = t247 * t252;
t495 = -t446 * t366 - (-pkin(4) * t378 + qJ(4) * t248) * t372 - t109;
t317 = t251 * t346;
t365 = qJD(3) * t251;
t494 = (t317 - t365) * qJ(4) + t528;
t376 = t250 * t252;
t278 = -t246 * t248 + t247 * t376;
t147 = t278 * t372;
t363 = qJD(5) * t251;
t342 = t246 * t363;
t493 = -t250 * t366 + t147 - t342;
t277 = t246 * t378 + t248 * t250;
t148 = t277 * t372;
t492 = -t246 * t366 + t250 * t363 + t148;
t249 = sin(qJ(1));
t429 = cos(qJ(1));
t307 = t393 * t429;
t197 = t248 * t249 - t252 * t307;
t387 = t197 * t251;
t491 = -pkin(3) * t387 - t197 * t391;
t327 = t252 * t393;
t199 = t248 * t429 + t249 * t327;
t386 = t199 * t251;
t490 = -pkin(3) * t386 - t199 * t391;
t202 = pkin(1) * t327 - pkin(7) * t385;
t489 = -t63 * (-mrSges(5,2) * t247 - mrSges(5,3) * t251) - t142 * (mrSges(4,1) * t247 + mrSges(4,2) * t251);
t488 = -t247 * t503 - t251 * t504;
t300 = mrSges(5,2) * t251 - mrSges(5,3) * t247;
t305 = mrSges(4,1) * t251 - mrSges(4,2) * t247;
t487 = -t305 + t300;
t486 = t178 * t502 - t503 * t88 + t504 * t87;
t114 = t248 * t520 + t252 * t521;
t100 = pkin(8) * t281 + t114;
t107 = t188 * pkin(2) + (-qJDD(1) * pkin(1) - pkin(8) * t272) * t245;
t19 = t251 * t100 + t247 * t107 - t143 * t366 + t149 * t365;
t20 = -t247 * t100 + t107 * t251 - t143 * t365 - t149 * t366;
t484 = t19 * t251 - t20 * t247;
t14 = -qJ(4) * t178 - qJD(4) * t531 - t19;
t269 = qJDD(4) - t20;
t18 = -pkin(3) * t178 + t269;
t483 = -t14 * t251 + t18 * t247;
t414 = Ifges(3,4) * t248;
t463 = t245 ^ 2;
t482 = (t248 * (t252 * Ifges(3,1) - t414) / 0.2e1 - pkin(1) * (mrSges(3,1) * t248 + mrSges(3,2) * t252)) * t463;
t303 = t250 * mrSges(6,1) - t246 * mrSges(6,2);
t111 = Ifges(6,4) * t112;
t43 = Ifges(6,1) * t113 + Ifges(6,5) * t152 + t111;
t397 = t246 * t43;
t427 = pkin(4) * t157;
t86 = t143 * t251 + t149 * t247;
t69 = -qJ(4) * t531 - t86;
t48 = -t69 - t427;
t480 = t48 * t303 - t397 / 0.2e1;
t328 = t248 * t393;
t382 = t245 * t252;
t203 = pkin(1) * t328 + pkin(7) * t382;
t174 = pkin(8) * t393 + t203;
t374 = pkin(2) * t382 + pkin(8) * t385;
t428 = pkin(1) * t245;
t175 = -t374 - t428;
t184 = qJD(2) * t279;
t186 = t202 * qJD(2);
t270 = t174 * t366 - t175 * t365 - t247 * t184 - t251 * t186;
t370 = qJD(2) * t248;
t40 = -t245 * (qJ(4) * t370 - qJD(4) * t252) + t270;
t115 = -t248 * t521 + t252 * t520;
t101 = -pkin(2) * t281 - t115;
t254 = t87 * qJ(4) - t158 * qJD(4) + t101;
t11 = t447 * t88 + t254;
t8 = -pkin(4) * t87 - t178 * t447 + t269;
t1 = qJD(5) * t12 + t11 * t250 + t246 * t8;
t2 = -qJD(5) * t13 - t11 * t246 + t250 * t8;
t478 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t321 = mrSges(3,3) * t347;
t477 = -m(4) * t142 + mrSges(3,1) * t283 - mrSges(4,1) * t157 - mrSges(4,2) * t158 - t321;
t475 = -mrSges(6,3) - t479;
t474 = -t142 * mrSges(4,1) + t63 * mrSges(5,2);
t473 = m(6) * t446 - mrSges(3,2) - t506;
t265 = m(6) * qJ(4) + t302;
t472 = -m(5) * qJ(4) - t265 + t505;
t198 = t248 * t307 + t249 * t252;
t348 = t245 * t429;
t131 = t198 * t247 + t251 * t348;
t200 = -t249 * t328 + t252 * t429;
t383 = t245 * t251;
t135 = t200 * t247 - t249 * t383;
t195 = -t251 * t393 + t353;
t470 = g(1) * t135 + g(2) * t131 + g(3) * t195;
t468 = t247 * t302 + mrSges(3,1) - t487;
t467 = -t20 * mrSges(4,1) + t19 * mrSges(4,2) - t18 * mrSges(5,2) + t14 * mrSges(5,3);
t466 = -t303 - t473;
t440 = -t152 / 0.2e1;
t443 = -t113 / 0.2e1;
t445 = -t112 / 0.2e1;
t465 = -Ifges(6,5) * t443 - Ifges(6,6) * t445 - Ifges(6,3) * t440 + t476;
t452 = t84 / 0.2e1;
t457 = t37 / 0.2e1;
t458 = t36 / 0.2e1;
t464 = t530 + Ifges(5,4) * t525 + Ifges(6,5) * t458 + Ifges(5,6) * t449 + Ifges(6,6) * t457 + Ifges(6,3) * t452 - t434 * t504 + t478 + (-t450 + t451) * Ifges(5,2);
t6 = Ifges(6,4) * t36 + Ifges(6,2) * t37 + Ifges(6,6) * t84;
t461 = -t6 / 0.2e1;
t7 = Ifges(6,1) * t36 + Ifges(6,4) * t37 + Ifges(6,5) * t84;
t460 = -t7 / 0.2e1;
t410 = Ifges(6,4) * t113;
t42 = Ifges(6,2) * t112 + Ifges(6,6) * t152 + t410;
t456 = -t42 / 0.2e1;
t455 = t42 / 0.2e1;
t403 = t158 * Ifges(4,4);
t77 = -Ifges(4,2) * t157 + Ifges(4,6) * t531 + t403;
t453 = -t77 / 0.2e1;
t444 = t112 / 0.2e1;
t442 = t113 / 0.2e1;
t439 = t152 / 0.2e1;
t438 = -t157 / 0.2e1;
t437 = t157 / 0.2e1;
t436 = -t158 / 0.2e1;
t435 = t158 / 0.2e1;
t432 = t531 / 0.2e1;
t431 = -t531 / 0.2e1;
t426 = pkin(8) * t197;
t425 = pkin(8) * t199;
t422 = t195 * pkin(9);
t421 = t2 * t250;
t418 = mrSges(5,1) * t158;
t417 = mrSges(4,3) * t157;
t416 = mrSges(4,3) * t158;
t415 = mrSges(6,3) * t246;
t413 = Ifges(3,4) * t252;
t412 = Ifges(4,4) * t247;
t411 = Ifges(4,4) * t251;
t409 = Ifges(6,4) * t246;
t408 = Ifges(6,4) * t250;
t407 = Ifges(5,6) * t247;
t406 = Ifges(5,6) * t251;
t405 = t12 * t246;
t402 = t158 * Ifges(5,6);
t392 = qJ(4) * t157;
t390 = t158 * t250;
t389 = t197 * t246;
t388 = t197 * t250;
t384 = t245 * t249;
t380 = t246 * t251;
t377 = t250 * t251;
t106 = t251 * t174 + t247 * t175;
t373 = t429 * pkin(1) + pkin(7) * t384;
t369 = qJD(2) * t252;
t364 = qJD(5) * t250;
t237 = pkin(3) * t382;
t357 = pkin(8) * t366;
t356 = pkin(8) * t365;
t354 = qJ(4) * t382;
t351 = Ifges(3,5) * t260 - Ifges(3,6) * t188 + Ifges(3,3) * t281;
t350 = t200 * pkin(2) + t373;
t344 = t245 * t370;
t343 = t245 * t369;
t341 = t385 / 0.2e1;
t338 = t372 / 0.2e1;
t333 = -t364 / 0.2e1;
t57 = -t87 * mrSges(5,1) + t178 * mrSges(5,2);
t332 = -pkin(1) * t249 + pkin(7) * t348;
t190 = t197 * pkin(2);
t330 = pkin(8) * t198 - t190;
t192 = t199 * pkin(2);
t329 = pkin(8) * t200 - t192;
t189 = t195 * pkin(3);
t196 = t247 * t393 + t248 * t383;
t325 = t196 * qJ(4) - t189;
t105 = -t247 * t174 + t175 * t251;
t132 = t198 * t251 - t247 * t348;
t320 = mrSges(3,3) * t346;
t311 = t252 * t338;
t309 = -t198 * pkin(2) + t332;
t91 = -t105 + t237;
t306 = t195 * mrSges(4,1) + t196 * mrSges(4,2);
t129 = t195 * t250 + t246 * t382;
t130 = -t195 * t246 + t245 * t376;
t304 = mrSges(6,1) * t129 + mrSges(6,2) * t130;
t301 = -t195 * mrSges(5,2) - t196 * mrSges(5,3);
t299 = Ifges(4,1) * t251 - t412;
t298 = Ifges(6,1) * t250 - t409;
t297 = Ifges(6,1) * t246 + t408;
t296 = -Ifges(4,2) * t247 + t411;
t294 = -Ifges(6,2) * t246 + t408;
t293 = Ifges(6,2) * t250 + t409;
t291 = Ifges(6,5) * t250 - Ifges(6,6) * t246;
t290 = Ifges(6,5) * t246 + Ifges(6,6) * t250;
t289 = -Ifges(5,2) * t251 + t407;
t288 = Ifges(5,3) * t247 - t406;
t286 = -t13 * t250 + t405;
t58 = pkin(4) * t196 + pkin(9) * t382 + t91;
t173 = -pkin(2) * t393 - t202;
t89 = t173 - t325;
t64 = t89 + t422;
t23 = -t246 * t64 + t250 * t58;
t24 = t246 * t58 + t250 * t64;
t90 = t354 - t106;
t50 = -t174 * t365 - t175 * t366 + t184 * t251 - t247 * t186;
t136 = t200 * t251 + t247 * t384;
t273 = t136 * pkin(3) + qJ(4) * t135 + t350;
t261 = -pkin(3) * t132 - qJ(4) * t131 + t309;
t187 = t203 * qJD(2);
t259 = Ifges(3,6) * t393 + (t252 * Ifges(3,2) + t414) * t245;
t258 = -qJD(5) * t286 + t1 * t246 + t421;
t257 = t245 * t283 * (Ifges(3,5) * t252 - Ifges(3,6) * t248);
t128 = -t316 + (qJD(3) * t393 + t343) * t251;
t256 = -t128 * qJ(4) - t196 * qJD(4) + t187;
t230 = Ifges(3,4) * t346;
t212 = -pkin(3) * t251 + t331;
t201 = (-mrSges(3,1) * t252 + mrSges(3,2) * t248) * t245;
t181 = -mrSges(3,2) * t283 + t320;
t140 = Ifges(3,1) * t347 + Ifges(3,5) * t283 + t230;
t139 = Ifges(3,6) * qJD(2) + qJD(1) * t259;
t127 = qJD(3) * t196 + t247 * t343;
t119 = mrSges(5,2) * t531 + t418;
t117 = mrSges(4,1) * t531 - t416;
t116 = -mrSges(4,2) * t531 - t417;
t104 = -mrSges(5,2) * t157 - mrSges(5,3) * t158;
t102 = pkin(3) * t158 + t392;
t96 = t135 * t246 + t199 * t250;
t95 = t135 * t250 - t199 * t246;
t94 = -pkin(3) * t347 - t108;
t92 = -qJ(4) * t347 - t109;
t74 = Ifges(5,5) * t531 + Ifges(5,3) * t157 - t402;
t68 = t158 * t447 + t392;
t65 = -pkin(4) * t195 - t90;
t62 = qJD(5) * t129 + t127 * t246 + t250 * t344;
t61 = qJD(5) * t130 + t127 * t250 - t246 * t344;
t56 = mrSges(5,1) * t88 - mrSges(5,3) * t178;
t55 = -mrSges(4,2) * t178 - mrSges(4,3) * t88;
t54 = mrSges(4,1) * t178 + mrSges(4,3) * t87;
t53 = t86 - t427;
t46 = t127 * pkin(3) + t256;
t45 = -pkin(3) * t344 - t50;
t39 = mrSges(4,1) * t88 - mrSges(4,2) * t87;
t38 = -mrSges(5,2) * t88 + mrSges(5,3) * t87;
t33 = t127 * t447 + t256;
t26 = -pkin(4) * t127 - t40;
t25 = pkin(4) * t128 - t344 * t447 - t50;
t22 = t246 * t53 + t250 * t68;
t21 = -t246 * t68 + t250 * t53;
t15 = t88 * pkin(3) + t254;
t10 = -pkin(4) * t88 - t14;
t9 = -mrSges(6,1) * t37 + mrSges(6,2) * t36;
t4 = -qJD(5) * t24 - t246 * t33 + t25 * t250;
t3 = qJD(5) * t23 + t246 * t25 + t250 * t33;
t27 = [(mrSges(5,1) * t18 - mrSges(4,3) * t20 - Ifges(5,6) * t448 + t464 + t524 + t529) * t196 + (-t115 * t385 - t182 * t343 - t185 * t344 - t188 * t203 - t202 * t260) * mrSges(3,3) + (Ifges(3,3) * t393 + (Ifges(3,5) * t248 + Ifges(3,6) * t252) * t245) * t509 + (Ifges(3,5) * t393 + (t248 * Ifges(3,1) + t413) * t245) * t510 + t259 * t511 + (t257 / 0.2e1 + t499 * t341) * qJD(2) + (t1 * t129 - t12 * t62 + t13 * t61 + t130 * t2) * mrSges(6,3) + (-Ifges(6,1) * t130 + Ifges(6,4) * t129) * t458 + (Ifges(6,1) * t62 + Ifges(6,4) * t61) * t442 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t463 + t114 * t203 + t115 * t202 + t185 * t186) + m(4) * (t101 * t173 + t105 * t20 + t106 * t19 - t270 * t86 - t50 * t85) - t270 * t116 + (-m(6) * t261 + t388 * mrSges(6,1) - t389 * mrSges(6,2) - m(3) * t332 + t198 * mrSges(3,1) - mrSges(3,3) * t348 + t249 * mrSges(2,1) + mrSges(2,2) * t429 - m(4) * (t309 - t426) - m(5) * (t261 - t426) - (-t302 + t505) * t131 + t473 * t197 - t475 * t132) * g(1) + t46 * t104 + t105 * t54 + t106 * t55 + t89 * t38 + t90 * t56 + t91 * t57 + (-m(3) * t182 - t477) * t187 + t186 * t181 + t173 * t39 + t129 * t6 / 0.2e1 + (Ifges(3,4) * t510 - Ifges(5,4) * t450 - Ifges(4,5) * t451 - Ifges(5,5) * t448 + Ifges(3,2) * t511 + Ifges(3,6) * t509 - Ifges(4,6) * t449 - t434 * t502 + t467 + t114 * mrSges(3,3) - t486 / 0.2e1) * t382 + m(6) * (t1 * t24 + t10 * t65 + t12 * t4 + t13 * t3 + t2 * t23 + t26 * t48) + m(5) * (t14 * t90 + t15 * t89 + t18 * t91 + t40 * t69 + t45 * t67 + t46 * t63) + t15 * t301 - t10 * t304 + t3 * t70 + t4 * t71 + t61 * t455 + t130 * t460 + (t14 * mrSges(5,1) - t19 * mrSges(4,3) - Ifges(4,4) * t451 + Ifges(5,6) * t450 + t516) * t195 + t393 * t351 / 0.2e1 + (t463 * qJD(1) * (-Ifges(3,2) * t248 + t413) + t245 * t140) * t369 / 0.2e1 + (t115 * t393 - t188 * t428 + t202 * t281) * mrSges(3,1) + (Ifges(3,1) * t260 - Ifges(3,4) * t188 + Ifges(3,5) * t281) * t341 + t101 * t306 + (Ifges(4,1) * t435 + Ifges(4,4) * t438 + Ifges(6,5) * t442 - Ifges(5,2) * t436 - Ifges(5,6) * t437 + Ifges(6,6) * t444 + Ifges(6,3) * t439 - t432 * t504 + t514 + t518) * t128 + (-t114 * t393 - t203 * t281 - t260 * t428) * mrSges(3,2) + t482 * t361 - pkin(1) * t201 * t360 + t48 * (-mrSges(6,1) * t61 + mrSges(6,2) * t62) + t62 * t43 / 0.2e1 + t65 * t9 + t26 * t51 + (-Ifges(6,4) * t130 + Ifges(6,2) * t129) * t457 + (Ifges(6,4) * t62 + Ifges(6,2) * t61) * t444 + t24 * t17 + t23 * t16 + Ifges(2,3) * qJDD(1) + t50 * t117 + t40 * t118 + t45 * t119 + (Ifges(4,6) * t438 + Ifges(4,5) * t435 + Ifges(5,4) * t436 + Ifges(5,5) * t437 - t85 * mrSges(4,1) - t139 / 0.2e1 - t69 * mrSges(5,3) - t86 * mrSges(4,2) + t67 * mrSges(5,2) + t502 * t432) * t344 + (t74 / 0.2e1 - Ifges(4,2) * t438 + t453 - Ifges(4,4) * t435 + Ifges(5,6) * t436 + Ifges(5,3) * t437 + t69 * mrSges(5,1) - t86 * mrSges(4,3) - t503 * t432 - t474) * t127 + (Ifges(6,5) * t62 + Ifges(6,6) * t61) * t439 + (-Ifges(6,5) * t130 + Ifges(6,6) * t129) * t452 + (-m(3) * t373 - t200 * mrSges(3,1) - mrSges(3,3) * t384 - m(6) * t273 - t96 * mrSges(6,1) - t95 * mrSges(6,2) - mrSges(2,1) * t429 + t249 * mrSges(2,2) - m(4) * (t350 + t425) - m(5) * (t273 + t425) + t505 * t135 - t473 * t199 + t475 * t136) * g(2); (-t291 * t439 - t294 * t444 - t298 * t442) * t363 + (t357 - t92) * t118 + (-t356 - t108) * t117 + (-t357 - t109) * t116 + t514 * t365 + ((t299 / 0.2e1 - t289 / 0.2e1) * t158 + (-t296 / 0.2e1 + t288 / 0.2e1) * t157 + (Ifges(6,3) * t251 + t247 * t290) * t439 + (Ifges(6,5) * t251 + t247 * t297) * t442 + (Ifges(6,6) * t251 + t247 * t293) * t444 - t489 + t488 * t432) * qJD(3) + (-t181 + t320) * t182 + (t158 * (Ifges(5,4) * t248 + t252 * t289) + t157 * (Ifges(4,6) * t248 + t252 * t296) + t248 * t139 - (t248 * t502 + t252 * t488) * t531) * t338 + (Ifges(6,1) * t148 + Ifges(6,4) * t147) * t443 + (t366 * t69 + t483) * mrSges(5,1) + (-t366 * t86 + t484) * mrSges(4,3) + t222 * t9 + t212 * t38 + (t321 + t477) * t185 - t148 * t43 / 0.2e1 - t465 * t317 + t500 * t71 + (t1 * t122 + t10 * t222 + t12 * t500 + t121 * t2 + t13 * t501 + t48 * t495) * m(6) + t501 * t70 + t15 * t300 + t342 * t455 + t147 * t456 + t380 * t460 + t377 * t461 + (t10 * t303 - t290 * t452 - t293 * t457 - t297 * t458 + t76 * t311 + t43 * t333 - t516) * t251 + t366 * t453 + (t356 - t94) * t119 + (Ifges(6,5) * t148 + Ifges(6,6) * t147) * t440 - t406 * t450 + t411 * t451 + (-m(4) * t374 + t201 - t512 * (t251 * t237 + t247 * t354 + t374) + (-t277 * mrSges(6,1) - t278 * mrSges(6,2) + (-mrSges(6,3) - t513) * t375 + t487 * t252 + (-m(6) * pkin(4) + t506) * t248) * t245) * g(3) + t495 * t51 + (-t257 / 0.2e1 - t482 * qJD(1)) * qJD(1) + (t250 * t42 + t397 + t74) * t366 / 0.2e1 - t101 * t305 + (-t1 * t377 + t12 * t492 - t13 * t493 + t2 * t380) * mrSges(6,3) + t494 * t104 + t489 * t346 + (-m(5) * (t329 + t490) - m(4) * t329 - m(6) * (-pkin(9) * t386 - t192 + t490) + mrSges(6,3) * t386 + t466 * t200 + t468 * t199) * g(1) + (-m(5) * (t330 + t491) - m(4) * t330 - m(6) * (-pkin(9) * t387 - t190 + t491) + mrSges(6,3) * t387 + t466 * t198 + t468 * t197) * g(2) + (mrSges(6,1) * t493 - mrSges(6,2) * t492) * t48 - t407 * t448 + t412 * t449 + (t15 * t212 + t494 * t63 - t67 * t94 - t69 * t92) * m(5) + (-pkin(2) * t101 + t108 * t85 - t109 * t86) * m(4) + ((t55 - t56) * t251 + (t57 - t54) * t247 + ((t247 * t69 + t251 * t67) * qJD(3) + t483) * m(5) + ((-t247 * t86 + t251 * t85) * qJD(3) + t484) * m(4)) * pkin(8) - (t158 * (Ifges(4,5) * t248 + t252 * t299) + t157 * (Ifges(5,5) * t248 + t252 * t288) + t499 * t248 + (-Ifges(3,2) * t347 + t247 * t74 + t251 * t498 + t140 + t230) * t252) * t372 / 0.2e1 + (t85 * (mrSges(4,1) * t248 - mrSges(4,3) * t375) - t67 * (mrSges(5,1) * t375 + mrSges(5,2) * t248) - t86 * (-mrSges(4,2) * t248 - mrSges(4,3) * t378) - t69 * (mrSges(5,1) * t378 - mrSges(5,3) * t248)) * t372 - pkin(2) * t39 + t351 + (Ifges(6,4) * t148 + Ifges(6,2) * t147) * t445 - t114 * mrSges(3,2) + t115 * mrSges(3,1) + (t311 * t77 + t464) * t247 + t121 * t16 + t122 * t17 + t247 * t524; (-pkin(3) * t18 - qJ(4) * t14 - qJD(4) * t69 - t102 * t63) * m(5) + (-m(5) * t325 + t301 + t306 - m(6) * (-t189 - t422) - t265 * t196) * g(3) + t486 + (t10 * qJ(4) - t12 * t21 - t13 * t22 + t517 * t48) * m(6) + t280 * t51 - t102 * t104 + t480 * qJD(5) + (qJD(5) * t405 - t421 + (-t364 - t390) * t13 + t470) * mrSges(6,3) + (-Ifges(4,2) * t158 - t151 + t498) * t437 + t10 * t302 - t22 * t70 - t21 * t71 + t390 * t456 + t294 * t457 + t298 * t458 + t246 * t461 - t467 + (-m(6) * t258 - t515) * t447 + t291 * t452 + (-m(5) * t67 + t117 - t119 + t416) * t86 + (-t403 + t74) * t436 + t496 * qJD(4) + t67 * t419 - (t112 * t293 + t113 * t297 + t152 * t290) * qJD(5) / 0.2e1 + (-t56 + t9) * qJ(4) + (-Ifges(4,1) * t436 + Ifges(5,2) * t435 + t431 * t504 + t465 + t518) * t157 + (t131 * t519 + t132 * t472) * g(2) + (t135 * t519 + t136 * t472) * g(1) + t42 * t333 + (t402 + t77) * t435 + (-m(5) * t69 + t116 - t118 + t417) * t85 + (t76 + t150) * t438 - t69 * t418 - t1 * t415 - pkin(3) * t57 + t250 * t7 / 0.2e1 + (Ifges(5,3) * t438 + t12 * t415 + t290 * t440 + t293 * t445 + t297 * t443 - t431 * t503 + t474 + t480) * t158; -t496 * t531 + (t104 + t285) * t158 + t57 + (-t158 * t286 - t48 * t531 + t258 - t470) * m(6) + (t158 * t63 + t531 * t69 + t18 - t470) * m(5) + t515; -t48 * (mrSges(6,1) * t113 + mrSges(6,2) * t112) + (Ifges(6,1) * t112 - t410) * t443 + t42 * t442 + (Ifges(6,5) * t112 - Ifges(6,6) * t113) * t440 - t12 * t70 + t13 * t71 - g(1) * (mrSges(6,1) * t95 - mrSges(6,2) * t96) - g(2) * ((t131 * t250 - t389) * mrSges(6,1) + (-t131 * t246 - t388) * mrSges(6,2)) - g(3) * t304 + (t112 * t12 + t113 * t13) * mrSges(6,3) + t5 + (-Ifges(6,2) * t113 + t111 + t43) * t445 + t478;];
tau = t27;
