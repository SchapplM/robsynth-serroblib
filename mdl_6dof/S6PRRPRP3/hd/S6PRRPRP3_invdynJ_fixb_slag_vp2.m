% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:34:14
% EndTime: 2019-03-08 21:34:51
% DurationCPUTime: 23.72s
% Computational Cost: add. (7293->677), mult. (16850->888), div. (0->0), fcn. (12799->14), ass. (0->319)
t284 = cos(qJ(3));
t366 = qJD(2) * t284;
t264 = qJD(5) - t366;
t413 = -t264 / 0.2e1;
t474 = Ifges(6,3) + Ifges(7,2);
t510 = t474 * t413;
t281 = sin(qJ(3));
t313 = pkin(3) * t281 - qJ(4) * t284;
t363 = qJD(4) * t281;
t201 = qJD(3) * t313 - t363;
t275 = sin(pkin(11));
t277 = cos(pkin(11));
t282 = sin(qJ(2));
t365 = qJD(3) * t281;
t354 = pkin(8) * t365;
t276 = sin(pkin(6));
t371 = qJD(1) * t276;
t285 = cos(qJ(2));
t378 = t284 * t285;
t458 = t277 * t201 + t275 * t354 - (-t275 * t378 + t277 * t282) * t371;
t509 = t275 * t201 - (t275 * t282 + t277 * t378) * t371;
t274 = pkin(11) + qJ(5);
t270 = sin(t274);
t271 = cos(t274);
t325 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t499 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t508 = t270 * t499 + t325 * t271;
t368 = qJD(2) * t281;
t225 = qJD(3) * t275 + t277 * t368;
t280 = sin(qJ(5));
t283 = cos(qJ(5));
t304 = -t277 * qJD(3) + t275 * t368;
t127 = t280 * t225 + t283 * t304;
t359 = qJD(2) * qJD(3);
t234 = qJDD(2) * t281 + t284 * t359;
t184 = qJDD(3) * t277 - t234 * t275;
t185 = qJDD(3) * t275 + t234 * t277;
t45 = -qJD(5) * t127 + t280 * t184 + t283 * t185;
t507 = -t45 / 0.2e1;
t290 = t283 * t225 - t280 * t304;
t46 = qJD(5) * t290 - t283 * t184 + t280 * t185;
t431 = -t46 / 0.2e1;
t425 = -t127 / 0.2e1;
t476 = Ifges(7,4) + Ifges(6,5);
t506 = -t476 / 0.2e1;
t481 = qJD(3) / 0.2e1;
t267 = pkin(4) * t277 + pkin(3);
t487 = -m(6) - m(7);
t504 = t487 * t267;
t380 = t277 * t284;
t307 = pkin(4) * t281 - pkin(9) * t380;
t503 = qJD(3) * t307 + t458;
t381 = t277 * t281;
t384 = t275 * t284;
t502 = (-pkin(8) * t381 - pkin(9) * t384) * qJD(3) + t509;
t478 = Ifges(6,1) + Ifges(7,1);
t477 = -Ifges(6,4) + Ifges(7,5);
t475 = -Ifges(6,6) + Ifges(7,6);
t229 = t275 * t283 + t277 * t280;
t303 = t229 * t284;
t182 = qJD(2) * t303;
t209 = t229 * qJD(5);
t456 = t182 - t209;
t386 = t275 * t280;
t228 = -t283 * t277 + t386;
t302 = t228 * t284;
t183 = qJD(2) * t302;
t208 = t228 * qJD(5);
t455 = -t183 + t208;
t348 = t282 * t371;
t235 = qJD(2) * pkin(8) + t348;
t278 = cos(pkin(6));
t370 = qJD(1) * t278;
t170 = t284 * t235 + t281 * t370;
t159 = qJD(3) * qJ(4) + t170;
t239 = -pkin(3) * t284 - qJ(4) * t281 - pkin(2);
t347 = t285 * t371;
t172 = qJD(2) * t239 - t347;
t82 = -t159 * t275 + t277 * t172;
t55 = -pkin(4) * t366 - pkin(9) * t225 + t82;
t83 = t277 * t159 + t275 * t172;
t62 = -pkin(9) * t304 + t83;
t22 = -t280 * t62 + t283 * t55;
t460 = qJD(6) - t22;
t13 = -pkin(5) * t264 + t460;
t23 = t280 * t55 + t283 * t62;
t15 = qJ(6) * t264 + t23;
t404 = Ifges(4,4) * t281;
t317 = t284 * Ifges(4,2) + t404;
t501 = -t15 * mrSges(7,3) - t22 * mrSges(6,1) + Ifges(4,6) * t481 + qJD(2) * t317 / 0.2e1 + t13 * mrSges(7,1) + t23 * mrSges(6,2) - Ifges(5,5) * t225 / 0.2e1 + Ifges(5,6) * t304 / 0.2e1 + Ifges(5,3) * t366 / 0.2e1 + t475 * t425 + t510 + t290 * t506;
t233 = -t284 * qJDD(2) + t281 * t359;
t227 = qJDD(5) + t233;
t221 = t281 * t235;
t357 = qJDD(1) * t278;
t369 = qJD(2) * t276;
t341 = qJD(1) * t369;
t251 = t285 * t341;
t358 = qJDD(1) * t276;
t198 = t282 * t358 + t251;
t495 = qJDD(2) * pkin(8) + qJD(3) * t370 + t198;
t349 = t281 * t357 + t284 * t495;
t68 = qJDD(3) * qJ(4) + (qJD(4) - t221) * qJD(3) + t349;
t250 = t282 * t341;
t197 = t285 * t358 - t250;
t186 = -qJDD(2) * pkin(2) - t197;
t87 = pkin(3) * t233 - qJ(4) * t234 - qJD(2) * t363 + t186;
t26 = -t275 * t68 + t277 * t87;
t19 = pkin(4) * t233 - pkin(9) * t185 + t26;
t27 = t275 * t87 + t277 * t68;
t25 = pkin(9) * t184 + t27;
t4 = -qJD(5) * t23 + t19 * t283 - t25 * t280;
t2 = -pkin(5) * t227 + qJDD(6) - t4;
t416 = t227 / 0.2e1;
t430 = t46 / 0.2e1;
t432 = t45 / 0.2e1;
t364 = qJD(3) * t284;
t74 = -t235 * t364 - t281 * t495 + t284 * t357;
t69 = -qJDD(3) * pkin(3) + qJDD(4) - t74;
t47 = -pkin(4) * t184 + t69;
t5 = pkin(5) * t46 - qJ(6) * t45 - qJD(6) * t290 + t47;
t500 = -mrSges(6,2) * t47 - mrSges(7,2) * t2 + mrSges(6,3) * t4 + mrSges(7,3) * t5 - Ifges(7,5) * t430 + t227 * t506 - t416 * t476 + (-t432 + t507) * t478 + (-Ifges(6,4) + t477) * t431;
t125 = Ifges(6,4) * t127;
t399 = Ifges(7,5) * t127;
t471 = t476 * t264 + t290 * t478 - t125 + t399;
t498 = -t471 / 0.2e1;
t497 = -m(4) + t487;
t60 = mrSges(7,1) * t127 - mrSges(7,3) * t290;
t61 = mrSges(6,1) * t127 + mrSges(6,2) * t290;
t496 = t60 + t61;
t224 = t277 * t239;
t132 = -pkin(9) * t381 + t224 + (-pkin(8) * t275 - pkin(4)) * t284;
t175 = pkin(8) * t380 + t275 * t239;
t385 = t275 * t281;
t150 = -pkin(9) * t385 + t175;
t459 = t280 * t132 + t283 * t150;
t469 = -qJD(5) * t459 - t280 * t502 + t283 * t503;
t360 = qJD(5) * t283;
t362 = qJD(5) * t280;
t467 = t132 * t360 - t150 * t362 + t280 * t503 + t283 * t502;
t407 = pkin(9) + qJ(4);
t241 = t407 * t275;
t242 = t407 * t277;
t308 = -t283 * t241 - t242 * t280;
t169 = t284 * t370 - t221;
t231 = t313 * qJD(2);
t107 = -t169 * t275 + t277 * t231;
t79 = qJD(2) * t307 + t107;
t108 = t277 * t169 + t275 * t231;
t346 = t275 * t366;
t89 = -pkin(9) * t346 + t108;
t464 = -qJD(4) * t228 + qJD(5) * t308 - t280 * t79 - t283 * t89;
t153 = -t241 * t280 + t242 * t283;
t463 = -qJD(4) * t229 - qJD(5) * t153 + t280 * t89 - t283 * t79;
t320 = t275 * mrSges(5,1) + t277 * mrSges(5,2);
t494 = -t320 + mrSges(3,2) - mrSges(4,3);
t155 = -qJD(3) * pkin(3) + qJD(4) - t169;
t109 = pkin(4) * t304 + t155;
t33 = t127 * pkin(5) - qJ(6) * t290 + t109;
t124 = Ifges(7,5) * t290;
t49 = Ifges(7,6) * t264 + Ifges(7,3) * t127 + t124;
t400 = Ifges(6,4) * t290;
t52 = -Ifges(6,2) * t127 + Ifges(6,6) * t264 + t400;
t493 = t33 * mrSges(7,1) + t49 / 0.2e1 - t52 / 0.2e1;
t414 = t233 / 0.2e1;
t419 = t185 / 0.2e1;
t492 = Ifges(5,1) * t419 + Ifges(5,5) * t414;
t410 = pkin(4) * t275;
t490 = -m(5) * pkin(8) - t325 * t270 + t271 * t499 + t487 * t410 + t494;
t321 = -mrSges(5,1) * t277 + mrSges(5,2) * t275;
t301 = m(5) * pkin(3) - t321;
t322 = mrSges(4,1) * t284 - mrSges(4,2) * t281;
t342 = m(5) * qJ(4) + mrSges(5,3);
t446 = t281 * t342 + t284 * t301 + mrSges(3,1) + t322;
t479 = mrSges(6,3) + mrSges(7,2);
t489 = t446 + (-t487 * t407 + t479) * t281 + (-t504 + t508) * t284;
t488 = -m(5) - m(4);
t420 = t184 / 0.2e1;
t485 = -t233 / 0.2e1;
t484 = t234 / 0.2e1;
t483 = mrSges(7,3) * t33;
t34 = mrSges(6,1) * t227 - mrSges(6,3) * t45;
t35 = -t227 * mrSges(7,1) + t45 * mrSges(7,2);
t473 = t35 - t34;
t36 = -mrSges(6,2) * t227 - mrSges(6,3) * t46;
t37 = -mrSges(7,2) * t46 + mrSges(7,3) * t227;
t472 = t36 + t37;
t470 = -pkin(5) * t365 - t469;
t468 = qJ(6) * t365 - qJD(6) * t284 + t467;
t137 = pkin(4) * t346 + t170;
t466 = -pkin(5) * t456 + qJ(6) * t455 - qJD(6) * t229 - t137;
t465 = -qJ(6) * t368 + t464;
t462 = pkin(5) * t368 - t463;
t461 = mrSges(4,2) - t342;
t457 = -t277 * t354 + t509;
t353 = mrSges(4,3) * t368;
t454 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t304 + t225 * mrSges(5,2) + t353;
t401 = Ifges(5,4) * t277;
t316 = -Ifges(5,2) * t275 + t401;
t402 = Ifges(5,4) * t275;
t318 = Ifges(5,1) * t277 - t402;
t453 = t225 * (Ifges(5,5) * t281 + t284 * t318) - (Ifges(5,6) * t281 + t284 * t316) * t304;
t452 = t227 * t474 + t45 * t476 + t46 * t475;
t73 = -t235 * t365 + t349;
t451 = -t281 * t74 + t284 * t73;
t450 = -t26 * t275 + t27 * t277;
t10 = t46 * mrSges(7,1) - t45 * mrSges(7,3);
t11 = t46 * mrSges(6,1) + t45 * mrSges(6,2);
t97 = -t184 * mrSges(5,1) + t185 * mrSges(5,2);
t449 = t10 + t11 + t97;
t444 = mrSges(4,1) + t301 + t508;
t236 = -qJD(2) * pkin(2) - t347;
t443 = -t155 * t284 * t320 - t83 * (-mrSges(5,2) * t281 - mrSges(5,3) * t384) - t82 * (mrSges(5,1) * t281 - mrSges(5,3) * t380) - t236 * (mrSges(4,1) * t281 + mrSges(4,2) * t284);
t442 = -m(5) * t155 - t454;
t392 = sin(pkin(10));
t328 = t392 * t285;
t393 = cos(pkin(10));
t331 = t393 * t282;
t205 = t278 * t331 + t328;
t333 = t276 * t393;
t145 = t205 * t284 - t281 * t333;
t329 = t392 * t282;
t330 = t393 * t285;
t207 = -t278 * t329 + t330;
t332 = t276 * t392;
t147 = t207 * t284 + t281 * t332;
t383 = t276 * t282;
t212 = t278 * t281 + t284 * t383;
t441 = -g(1) * t147 - g(2) * t145 - g(3) * t212;
t3 = t280 * t19 + t283 * t25 + t55 * t360 - t362 * t62;
t1 = qJ(6) * t227 + qJD(6) * t264 + t3;
t440 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t412 = t264 / 0.2e1;
t422 = t290 / 0.2e1;
t424 = t127 / 0.2e1;
t438 = Ifges(6,4) * t425 + Ifges(7,5) * t424 + t412 * t476 + t422 * t478 - t483;
t437 = -Ifges(6,2) * t425 + Ifges(7,3) * t424 + t412 * t475 + t422 * t477 + t493;
t435 = mrSges(6,1) * t47 + mrSges(7,1) * t5 + 0.2e1 * Ifges(7,3) * t430 + Ifges(6,4) * t507 - t227 * Ifges(6,6) / 0.2e1 + (t477 + Ifges(7,5)) * t432 + (t475 + Ifges(7,6)) * t416 + (-t431 + t430) * Ifges(6,2);
t286 = qJD(2) ^ 2;
t427 = Ifges(5,4) * t420 + t492;
t423 = -t290 / 0.2e1;
t273 = t281 * pkin(8);
t406 = mrSges(6,3) * t127;
t405 = mrSges(6,3) * t290;
t403 = Ifges(4,4) * t284;
t396 = t281 * t69;
t382 = t276 * t285;
t100 = -mrSges(6,2) * t264 - t406;
t103 = -mrSges(7,2) * t127 + mrSges(7,3) * t264;
t377 = t100 + t103;
t101 = mrSges(6,1) * t264 - t405;
t102 = -mrSges(7,1) * t264 + mrSges(7,2) * t290;
t376 = -t102 + t101;
t372 = pkin(2) * t382 + pkin(8) * t383;
t269 = pkin(8) * t364;
t343 = t275 * t364;
t220 = pkin(4) * t343 + t269;
t232 = pkin(4) * t385 + t273;
t367 = qJD(2) * t282;
t361 = qJD(5) * t281;
t352 = mrSges(4,3) * t366;
t351 = t281 * t382;
t350 = t270 * t382;
t345 = t276 * t367;
t344 = t285 * t369;
t338 = -t366 / 0.2e1;
t327 = -t359 / 0.2e1;
t326 = t359 / 0.2e1;
t315 = Ifges(4,5) * t284 - Ifges(4,6) * t281;
t314 = Ifges(5,5) * t277 - t275 * Ifges(5,6);
t63 = t132 * t283 - t150 * t280;
t142 = -t212 * t275 - t277 * t382;
t143 = t212 * t277 - t275 * t382;
t309 = t283 * t142 - t143 * t280;
t67 = t142 * t280 + t143 * t283;
t211 = -t278 * t284 + t281 * t383;
t305 = t281 * (Ifges(4,1) * t284 - t404);
t144 = t205 * t281 + t284 * t333;
t146 = t207 * t281 - t284 * t332;
t299 = -g(1) * t146 - g(2) * t144 - g(3) * t211;
t291 = t284 * (Ifges(5,3) * t281 + t284 * t314);
t268 = Ifges(4,4) * t366;
t245 = -qJD(3) * mrSges(4,2) + t352;
t230 = t322 * qJD(2);
t216 = Ifges(4,1) * t368 + Ifges(4,5) * qJD(3) + t268;
t206 = t278 * t328 + t331;
t204 = -t278 * t330 + t329;
t200 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t234;
t199 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t233;
t196 = t206 * pkin(2);
t195 = t204 * pkin(2);
t194 = t228 * t281;
t193 = t229 * t281;
t181 = -mrSges(5,1) * t366 - mrSges(5,3) * t225;
t180 = mrSges(5,2) * t366 - mrSges(5,3) * t304;
t174 = -pkin(8) * t384 + t224;
t154 = mrSges(4,1) * t233 + mrSges(4,2) * t234;
t149 = -qJD(3) * t211 + t284 * t344;
t148 = qJD(3) * t212 + t281 * t344;
t130 = t212 * t270 + t271 * t382;
t123 = pkin(5) * t228 - qJ(6) * t229 - t267;
t119 = Ifges(5,1) * t225 - Ifges(5,4) * t304 - Ifges(5,5) * t366;
t118 = Ifges(5,4) * t225 - Ifges(5,2) * t304 - Ifges(5,6) * t366;
t116 = mrSges(5,1) * t233 - mrSges(5,3) * t185;
t115 = -mrSges(5,2) * t233 + mrSges(5,3) * t184;
t111 = qJD(3) * t303 + t360 * t381 - t361 * t386;
t110 = -qJD(3) * t302 - t229 * t361;
t106 = t149 * t277 + t275 * t345;
t105 = -t149 * t275 + t277 * t345;
t88 = pkin(5) * t193 + qJ(6) * t194 + t232;
t77 = t147 * t270 - t206 * t271;
t75 = t145 * t270 - t204 * t271;
t71 = t185 * Ifges(5,4) + t184 * Ifges(5,2) + t233 * Ifges(5,6);
t59 = pkin(5) * t290 + qJ(6) * t127;
t58 = pkin(5) * t284 - t63;
t57 = -qJ(6) * t284 + t459;
t28 = pkin(5) * t111 - qJ(6) * t110 + qJD(6) * t194 + t220;
t18 = qJD(5) * t67 - t283 * t105 + t106 * t280;
t17 = qJD(5) * t309 + t105 * t280 + t106 * t283;
t6 = [m(2) * qJDD(1) + t105 * t181 + t106 * t180 + t143 * t115 + t142 * t116 + t149 * t245 + t212 * t199 + t472 * t67 - t473 * t309 - t376 * t18 + t377 * t17 + (-t200 + t449) * t211 + (t454 + t496) * t148 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t286 - t154) * t285 + (-mrSges(3,1) * t286 - mrSges(3,2) * qJDD(2) - qJD(2) * t230) * t282) * t276 + (-m(2) - m(3) + t487 + t488) * g(3) + m(3) * (qJDD(1) * t278 ^ 2 + (t197 * t285 + t198 * t282) * t276) + m(4) * (-t148 * t169 + t149 * t170 - t211 * t74 + t212 * t73 + (-t186 * t285 + t236 * t367) * t276) + m(7) * (t1 * t67 + t13 * t18 + t148 * t33 + t15 * t17 - t2 * t309 + t211 * t5) + m(6) * (t109 * t148 + t17 * t23 - t18 * t22 + t211 * t47 + t3 * t67 + t309 * t4) + m(5) * (t105 * t82 + t106 * t83 + t142 * t26 + t143 * t27 + t148 * t155 + t211 * t69); (-m(6) * t109 - m(7) * t33 + t442 - t496) * t281 * t347 + t459 * t36 + (t109 * t220 + t22 * t469 + t23 * t467 + t232 * t47 + t3 * t459 + t4 * t63) * m(6) + (t488 * t372 - t479 * t351 + t487 * (t351 * t407 + t383 * t410 + t372) - t499 * (-t271 * t383 + t284 * t350) + (t378 * t504 - t325 * (t270 * t282 + t271 * t378) - t446 * t285 + t494 * t282) * t276) * g(3) + (-t169 * t364 + t451) * mrSges(4,3) + (t471 / 0.2e1 + t109 * mrSges(6,2) + t13 * mrSges(7,2) - t22 * mrSges(6,3) + t438) * t110 + (m(5) * t196 + t497 * (pkin(8) * t207 - t196) + t490 * t207 + t489 * t206) * g(1) + (m(5) * t195 + t497 * (pkin(8) * t205 - t195) + t490 * t205 + t489 * t204) * g(2) + (-t26 * t381 - t27 * t385) * mrSges(5,3) - t245 * t354 + t381 * t427 + (t277 * t119 + t216) * t364 / 0.2e1 + (t251 - t198) * mrSges(3,2) - t71 * t385 / 0.2e1 + t230 * t348 - t118 * t343 / 0.2e1 + t305 * t326 + t291 * t327 + t453 * t481 + t403 * t484 + t317 * t485 + (t250 + t197) * mrSges(3,1) + (-t200 + t97) * t273 + qJDD(3) * (Ifges(4,5) * t281 + Ifges(4,6) * t284) + t232 * t11 + Ifges(3,3) * qJDD(2) + t220 * t61 + t174 * t116 + t175 * t115 - pkin(2) * t154 + t88 * t10 + t57 * t37 + t58 * t35 - t186 * t322 + t28 * t60 + t63 * t34 + (-Ifges(7,6) * t430 - Ifges(6,6) * t431 + Ifges(4,4) * t484 + Ifges(4,2) * t485 + t27 * mrSges(5,2) + pkin(8) * t199 - t245 * t347 + (-Ifges(4,2) * t281 + t403) * t326 - t26 * mrSges(5,1) - Ifges(5,3) * t414 - Ifges(5,5) * t419 - Ifges(5,6) * t420 - t476 * t432 - t474 * t416 + t440) * t284 + (Ifges(4,1) * t234 + Ifges(4,4) * t485 + t314 * t414 + t316 * t420 + t318 * t419) * t281 + (m(4) * (-t169 * t284 - t170 * t281) * pkin(8) + t315 * t481 - t443) * qJD(3) + t467 * t100 + t468 * t103 + t469 * t101 + (t1 * t57 + t13 * t470 + t15 * t468 + t2 * t58 + t28 * t33 + t5 * t88) * m(7) + t470 * t102 + t457 * t180 + t458 * t181 + (t174 * t26 + t175 * t27 + (t155 * t364 + t396) * pkin(8) + t457 * t83 + t458 * t82) * m(5) + (-pkin(2) * t186 + t451 * pkin(8) - (t236 * t282 + (-t169 * t281 + t170 * t284) * t285) * t371) * m(4) - (Ifges(5,5) * t185 + Ifges(5,6) * t184 + Ifges(5,3) * t233 + t452) * t284 / 0.2e1 + t454 * t269 + (-mrSges(7,2) * t1 - mrSges(6,3) * t3 + t435) * t193 + (mrSges(6,1) * t109 - mrSges(7,2) * t15 - mrSges(6,3) * t23 + t437) * t111 + t500 * t194 + (-t170 * mrSges(4,3) + Ifges(6,6) * t425 + Ifges(7,6) * t424 + t412 * t474 + t422 * t476 - t501) * t365 + t320 * t396; -t473 * t308 + (t1 * t153 + t123 * t5 + t13 * t462 + t15 * t465 - t2 * t308 + t33 * t466) * m(7) + (-t109 * t137 + t153 * t3 + t22 * t463 + t23 * t464 - t267 * t47 + t308 * t4) * m(6) + (t487 * (-t211 * t267 + t212 * t407) + t461 * t212 + t444 * t211) * g(3) + (t487 * (-t146 * t267 + t147 * t407) + t461 * t147 + t444 * t146) * g(1) + (t487 * (-t144 * t267 + t145 * t407) + t461 * t145 + t444 * t144) * g(2) + (-qJ(4) * t116 + t427 + (-m(5) * t82 - t181) * qJD(4) + t492) * t275 + (-Ifges(6,2) * t424 + Ifges(7,3) * t425 + t475 * t413 + t477 * t423 - t493) * t182 - (Ifges(6,4) * t424 + Ifges(7,5) * t425 + t476 * t413 + t478 * t423 + t483 + t498) * t183 + (t498 - t438) * t208 + (-t1 * t228 - t13 * t455 + t15 * t456 + t441) * mrSges(7,2) + (-t453 / 0.2e1 + t443) * qJD(2) + (-Ifges(4,2) * t338 + Ifges(6,6) * t424 + Ifges(7,6) * t425 + t476 * t423 + t501 + t510) * t368 + t118 * t346 / 0.2e1 + t315 * t327 + (t352 - t245) * t169 - t267 * t11 + Ifges(4,5) * t234 - Ifges(4,6) * t233 + Ifges(4,3) * qJDD(3) - t108 * t180 - t107 * t181 - t137 * t61 + t123 * t10 - pkin(3) * t97 + (t119 * t338 + t71 / 0.2e1 + qJ(4) * t115 + Ifges(5,6) * t414 + Ifges(5,2) * t420 + (m(5) * t83 + t180) * qJD(4)) * t277 + t74 * mrSges(4,1) - t73 * mrSges(4,2) + t69 * t321 + (t216 + t268) * t338 + (-t305 / 0.2e1 + t291 / 0.2e1) * t286 + t472 * t153 + t465 * t103 + t466 * t60 + t462 * t102 + t463 * t101 + t464 * t100 + t402 * t420 + (-mrSges(6,1) * t456 - mrSges(6,2) * t455) * t109 + t450 * mrSges(5,3) + (-pkin(3) * t69 + qJ(4) * t450 - t107 * t82 - t108 * t83) * m(5) + (t353 + t442) * t170 + t401 * t419 + t435 * t228 + t437 * t209 + (t22 * t455 - t228 * t3 + t23 * t456 + t441) * mrSges(6,3) - t500 * t229; t376 * t290 + t377 * t127 + t304 * t180 + t225 * t181 + (t127 * t15 - t13 * t290 + t299 + t5) * m(7) + (t127 * t23 + t22 * t290 + t299 + t47) * m(6) + (t82 * t225 + t304 * t83 + t299 + t69) * m(5) + t449; -t440 + (t127 * t13 + t15 * t290) * mrSges(7,2) + (Ifges(7,3) * t290 - t399) * t425 - t33 * (mrSges(7,1) * t290 + mrSges(7,3) * t127) - t109 * (mrSges(6,1) * t290 - mrSges(6,2) * t127) + (-t127 * t476 + t290 * t475) * t413 + (-Ifges(6,2) * t290 - t125 + t471) * t424 + (-t499 * (t145 * t271 + t204 * t270) + t325 * t75) * g(2) + (-t499 * (t147 * t271 + t206 * t270) + t325 * t77) * g(1) + (-t499 * (t212 * t271 - t350) + t325 * t130) * g(3) + (t376 + t405) * t23 + (-t377 - t406) * t22 + (-t127 * t478 + t124 - t400 + t49) * t423 + (-pkin(5) * t2 + qJ(6) * t1 - t13 * t23 + t15 * t460 - t33 * t59) * m(7) + qJD(6) * t103 - t59 * t60 + qJ(6) * t37 - pkin(5) * t35 + t52 * t422 + t452; -t264 * t103 + t290 * t60 + (-g(1) * t77 - g(2) * t75 - g(3) * t130 - t15 * t264 + t290 * t33 + t2) * m(7) + t35;];
tau  = t6;
