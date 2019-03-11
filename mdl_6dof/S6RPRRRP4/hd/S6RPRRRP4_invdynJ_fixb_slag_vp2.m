% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP4
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:06:59
% EndTime: 2019-03-09 06:07:39
% DurationCPUTime: 24.46s
% Computational Cost: add. (15557->672), mult. (38228->840), div. (0->0), fcn. (29875->14), ass. (0->308)
t285 = cos(qJ(5));
t543 = mrSges(6,1) + mrSges(7,1);
t544 = t285 * t543;
t542 = Ifges(6,4) + Ifges(7,4);
t277 = sin(pkin(10));
t283 = sin(qJ(3));
t278 = cos(pkin(10));
t287 = cos(qJ(3));
t393 = t278 * t287;
t223 = -t277 * t283 + t393;
t214 = t223 * qJD(3);
t224 = t277 * t287 + t278 * t283;
t180 = qJD(1) * t214 + qJDD(1) * t224;
t215 = t224 * qJD(3);
t181 = -qJD(1) * t215 + qJDD(1) * t223;
t282 = sin(qJ(4));
t286 = cos(qJ(4));
t212 = t223 * qJD(1);
t213 = t224 * qJD(1);
t333 = t286 * t212 - t213 * t282;
t102 = qJD(4) * t333 + t180 * t286 + t181 * t282;
t276 = qJD(3) + qJD(4);
t281 = sin(qJ(5));
t308 = t212 * t282 + t286 * t213;
t161 = t276 * t285 - t281 * t308;
t271 = qJDD(3) + qJDD(4);
t68 = qJD(5) * t161 + t102 * t285 + t271 * t281;
t467 = t68 / 0.2e1;
t162 = t276 * t281 + t285 * t308;
t69 = -qJD(5) * t162 - t102 * t281 + t271 * t285;
t466 = t69 / 0.2e1;
t103 = -qJD(4) * t308 - t180 * t282 + t181 * t286;
t98 = qJDD(5) - t103;
t465 = t98 / 0.2e1;
t537 = -mrSges(6,3) - mrSges(7,3);
t529 = Ifges(6,1) + Ifges(7,1);
t498 = Ifges(6,5) + Ifges(7,5);
t527 = Ifges(6,2) + Ifges(7,2);
t497 = Ifges(6,6) + Ifges(7,6);
t526 = Ifges(6,3) + Ifges(7,3);
t373 = qJD(1) * qJD(2);
t252 = qJ(2) * qJDD(1) + t373;
t375 = qJD(5) * t285;
t485 = t333 * t285;
t541 = -t375 + t485;
t376 = qJD(5) * t281;
t486 = t333 * t281;
t540 = t376 - t486;
t172 = qJD(5) - t333;
t536 = t542 * t161;
t494 = t162 * t529 + t172 * t498 + t536;
t539 = t494 / 0.2e1;
t538 = t498 * t465 + t542 * t466 + t467 * t529;
t499 = -mrSges(7,2) - mrSges(6,2);
t380 = t277 ^ 2 + t278 ^ 2;
t535 = t542 * t162;
t534 = t542 * t285;
t533 = t542 * t281;
t275 = pkin(10) + qJ(3);
t268 = qJ(4) + t275;
t258 = sin(t268);
t532 = t258 * t499;
t531 = m(7) * pkin(5);
t525 = t497 * t98 + t527 * t69 + t542 * t68;
t496 = t161 * t497 + t162 * t498 + t172 * t526;
t495 = t161 * t527 + t172 * t497 + t535;
t171 = Ifges(5,4) * t333;
t523 = Ifges(5,2) * t333;
t522 = t276 * Ifges(5,5);
t521 = t276 * Ifges(5,6);
t377 = qJD(4) * t286;
t366 = pkin(3) * t377;
t432 = pkin(7) + qJ(2);
t244 = t432 * t277;
t229 = qJD(1) * t244;
t245 = t432 * t278;
t230 = qJD(1) * t245;
t185 = -t229 * t283 + t230 * t287;
t160 = pkin(8) * t212 + t185;
t154 = t282 * t160;
t401 = t230 * t283;
t184 = -t287 * t229 - t401;
t159 = -pkin(8) * t213 + t184;
t115 = t159 * t286 - t154;
t137 = pkin(4) * t308 - pkin(9) * t333;
t446 = pkin(3) * t213;
t120 = t137 + t446;
t53 = t285 * t115 + t281 * t120;
t520 = t285 * t366 - t53;
t155 = t286 * t160;
t114 = t159 * t282 + t155;
t378 = qJD(4) * t282;
t519 = -pkin(3) * t378 + t114;
t484 = mrSges(5,1) * t276 + mrSges(6,1) * t161 - mrSges(6,2) * t162 - mrSges(5,3) * t308;
t156 = qJD(3) * pkin(3) + t159;
t111 = t156 * t286 - t154;
t104 = -pkin(4) * t276 - t111;
t430 = mrSges(7,2) * t285;
t318 = mrSges(7,1) * t281 + t430;
t319 = mrSges(6,1) * t281 + mrSges(6,2) * t285;
t73 = -pkin(5) * t161 + qJD(6) + t104;
t518 = t104 * t319 + t73 * t318;
t517 = -t281 * t497 + t285 * t498;
t516 = -t281 * t527 + t534;
t515 = t285 * t529 - t533;
t514 = t258 * t544;
t513 = t497 * t69 + t498 * t68 + t526 * t98;
t112 = t156 * t282 + t155;
t105 = pkin(9) * t276 + t112;
t260 = pkin(2) * t278 + pkin(1);
t235 = -qJD(1) * t260 + qJD(2);
t188 = -pkin(3) * t212 + t235;
t113 = -pkin(4) * t333 - pkin(9) * t308 + t188;
t334 = pkin(7) * qJDD(1) + t252;
t204 = t334 * t277;
t205 = t334 * t278;
t139 = -qJD(3) * t185 - t287 * t204 - t205 * t283;
t101 = qJDD(3) * pkin(3) - pkin(8) * t180 + t139;
t379 = qJD(3) * t287;
t138 = -qJD(3) * t401 - t283 * t204 + t287 * t205 - t229 * t379;
t108 = pkin(8) * t181 + t138;
t28 = t282 * t101 + t286 * t108 + t156 * t377 - t160 * t378;
t25 = pkin(9) * t271 + t28;
t234 = -qJDD(1) * t260 + qJDD(2);
t163 = -pkin(3) * t181 + t234;
t42 = -pkin(4) * t103 - pkin(9) * t102 + t163;
t5 = -t105 * t376 + t113 * t375 + t285 * t25 + t281 * t42;
t47 = t105 * t285 + t113 * t281;
t6 = -qJD(5) * t47 - t25 * t281 + t285 * t42;
t512 = -t281 * t6 + t285 * t5;
t511 = m(6) + m(7) + m(5);
t259 = cos(t268);
t510 = t259 * mrSges(5,1) + (-mrSges(5,2) - t537) * t258;
t508 = t540 * pkin(5);
t507 = t188 * mrSges(5,2) - t111 * mrSges(5,3);
t506 = qJ(6) * t486 + t285 * qJD(6);
t347 = m(3) * qJ(2) + mrSges(3,3);
t505 = -m(4) * t432 + mrSges(2,2) - mrSges(4,3) - mrSges(5,3) - t347;
t480 = t259 * pkin(4) + t258 * pkin(9);
t262 = pkin(5) * t285 + pkin(4);
t279 = -qJ(6) - pkin(9);
t482 = -t258 * t279 + t259 * t262;
t504 = -m(6) * t480 - m(7) * t482;
t266 = sin(t275);
t267 = cos(t275);
t323 = mrSges(4,1) * t267 - mrSges(4,2) * t266;
t324 = -mrSges(3,1) * t278 + mrSges(3,2) * t277;
t503 = m(3) * pkin(1) + m(4) * t260 + mrSges(2,1) + t323 - t324 + t510;
t3 = qJ(6) * t69 + qJD(6) * t161 + t5;
t502 = -t6 * mrSges(6,1) + t5 * mrSges(6,2) + t3 * mrSges(7,2);
t46 = -t105 * t281 + t285 * t113;
t39 = -qJ(6) * t162 + t46;
t31 = pkin(5) * t172 + t39;
t40 = qJ(6) * t161 + t47;
t501 = -t188 * mrSges(5,1) - t46 * mrSges(6,1) - t31 * mrSges(7,1) + t47 * mrSges(6,2) + t40 * mrSges(7,2) + t112 * mrSges(5,3);
t1 = pkin(5) * t98 - qJ(6) * t68 - qJD(6) * t162 + t6;
t29 = t101 * t286 - t282 * t108 - t156 * t378 - t160 * t377;
t26 = -pkin(4) * t271 - t29;
t12 = -pkin(5) * t69 + qJDD(6) + t26;
t342 = -t376 / 0.2e1;
t500 = Ifges(5,3) * t271 + Ifges(5,6) * t103 + Ifges(5,5) * t102 + t29 * mrSges(5,1) - t28 * mrSges(5,2) + t26 * (-mrSges(6,1) * t285 + mrSges(6,2) * t281) + t12 * (-mrSges(7,1) * t285 + mrSges(7,2) * t281) + (t281 * t529 + t534) * t467 + (t285 * t527 + t533) * t466 + (t281 * t498 + t285 * t497) * t465 + t281 * t538 + t525 * t285 / 0.2e1 + t375 * t539 + t518 * qJD(5) + (t161 * t516 + t162 * t515 + t172 * t517) * qJD(5) / 0.2e1 - t494 * t485 / 0.2e1 + (t342 + t486 / 0.2e1) * t495 + (-t1 * t281 + t285 * t3 + t31 * t541 - t40 * t540) * mrSges(7,3) + (t46 * t541 - t47 * t540 + t512) * mrSges(6,3);
t455 = -t308 / 0.2e1;
t443 = pkin(3) * t282;
t261 = pkin(9) + t443;
t388 = -qJ(6) - t261;
t330 = qJD(5) * t388;
t493 = t281 * t330 + t506 + t520;
t270 = t285 * qJ(6);
t472 = pkin(5) * t308 - t270 * t333;
t52 = -t115 * t281 + t285 * t120;
t492 = (-qJD(6) - t366) * t281 + t285 * t330 - t472 - t52;
t338 = qJD(5) * t279;
t55 = t285 * t111 + t281 * t137;
t491 = t281 * t338 + t506 - t55;
t54 = -t111 * t281 + t285 * t137;
t490 = -qJD(6) * t281 + t285 * t338 - t472 - t54;
t489 = t508 - t519;
t488 = -t112 + t508;
t487 = mrSges(7,1) + t531;
t189 = -t287 * t244 - t245 * t283;
t169 = -pkin(8) * t224 + t189;
t190 = -t283 * t244 + t287 * t245;
t170 = pkin(8) * t223 + t190;
t130 = t169 * t282 + t170 * t286;
t122 = t285 * t130;
t183 = t223 * t282 + t224 * t286;
t196 = -pkin(3) * t223 - t260;
t307 = t286 * t223 - t224 * t282;
t131 = -pkin(4) * t307 - pkin(9) * t183 + t196;
t61 = t281 * t131 + t122;
t483 = t286 * t169 - t170 * t282;
t305 = -t258 * t262 - t259 * t279;
t441 = pkin(4) * t258;
t444 = pkin(3) * t266;
t479 = -m(7) * (t305 - t444) - m(6) * (-t441 - t444) + t514;
t478 = -m(7) * t305 + t514;
t288 = cos(qJ(1));
t391 = t281 * t288;
t395 = t259 * t288;
t477 = t391 * t532 + t395 * t537;
t284 = sin(qJ(1));
t392 = t281 * t284;
t397 = t259 * t284;
t476 = t392 * t532 + t397 * t537;
t475 = g(1) * t288 + g(2) * t284;
t474 = mrSges(6,1) + t487;
t473 = -t510 + (-t281 * t499 - t544) * t259;
t429 = mrSges(6,3) * t161;
t124 = -mrSges(6,2) * t172 + t429;
t428 = mrSges(6,3) * t162;
t126 = mrSges(6,1) * t172 - t428;
t33 = mrSges(6,1) * t98 - mrSges(6,3) * t68;
t471 = m(6) * ((-t281 * t47 - t285 * t46) * qJD(5) + t512) - t126 * t375 - t124 * t376 - t281 * t33;
t450 = -t276 / 0.2e1;
t459 = -t172 / 0.2e1;
t461 = -t162 / 0.2e1;
t463 = -t161 / 0.2e1;
t470 = Ifges(5,1) * t455 + Ifges(5,5) * t450 + t459 * t517 + t461 * t515 + t463 * t516 - t507 - t518;
t456 = -t333 / 0.2e1;
t469 = -Ifges(5,2) * t456 - Ifges(5,6) * t450 + t459 * t526 + t498 * t461 + t497 * t463 + t501;
t462 = t161 / 0.2e1;
t460 = t162 / 0.2e1;
t458 = t172 / 0.2e1;
t454 = t308 / 0.2e1;
t451 = t213 / 0.2e1;
t445 = pkin(3) * t215;
t256 = pkin(3) * t267;
t442 = pkin(3) * t286;
t439 = pkin(9) * t285;
t427 = mrSges(7,3) * t161;
t426 = mrSges(7,3) * t162;
t425 = Ifges(4,4) * t213;
t424 = Ifges(5,4) * t308;
t415 = t184 * mrSges(4,3);
t414 = t185 * mrSges(4,3);
t140 = qJD(4) * t307 + t214 * t286 - t215 * t282;
t411 = t140 * t281;
t410 = t140 * t285;
t404 = t183 * t281;
t403 = t183 * t285;
t394 = t261 * t285;
t390 = t284 * t285;
t389 = t285 * t288;
t372 = qJDD(1) * t277;
t371 = qJDD(1) * t278;
t166 = -t244 * t379 + qJD(2) * t393 + (-qJD(2) * t277 - qJD(3) * t245) * t283;
t147 = -pkin(8) * t215 + t166;
t167 = -t224 * qJD(2) - qJD(3) * t190;
t148 = -pkin(8) * t214 + t167;
t50 = qJD(4) * t483 + t147 * t286 + t148 * t282;
t141 = qJD(4) * t183 + t214 * t282 + t286 * t215;
t72 = pkin(4) * t141 - pkin(9) * t140 + t445;
t368 = t131 * t375 + t281 * t72 + t285 * t50;
t356 = t183 * t375;
t22 = -t69 * mrSges(7,1) + t68 * mrSges(7,2);
t340 = -t103 * mrSges(5,1) + t102 * mrSges(5,2);
t339 = -t281 * t50 + t285 * t72;
t337 = -t181 * mrSges(4,1) + t180 * mrSges(4,2);
t60 = -t130 * t281 + t285 * t131;
t325 = -mrSges(3,1) * t371 + mrSges(3,2) * t372;
t320 = mrSges(5,1) * t258 + mrSges(5,2) * t259;
t310 = -t281 * t46 + t285 * t47;
t304 = -qJ(6) * t140 - qJD(6) * t183;
t202 = -t259 * t391 + t390;
t200 = t259 * t392 + t389;
t300 = t356 + t411;
t299 = t183 * t376 - t410;
t51 = qJD(4) * t130 + t147 * t282 - t286 * t148;
t272 = -pkin(8) - t432;
t265 = -qJDD(1) * pkin(1) + qJDD(2);
t263 = -pkin(4) - t442;
t248 = t270 + t439;
t247 = t279 * t281;
t246 = -t262 - t442;
t243 = pkin(9) * t395;
t242 = pkin(9) * t397;
t231 = t256 + t260;
t219 = t270 + t394;
t218 = t388 * t281;
t206 = Ifges(4,4) * t212;
t203 = t259 * t389 + t392;
t201 = -t259 * t390 + t391;
t192 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t213;
t191 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t212;
t174 = t213 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t206;
t173 = t212 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t425;
t164 = -mrSges(5,2) * t276 + mrSges(5,3) * t333;
t136 = -mrSges(5,1) * t333 + mrSges(5,2) * t308;
t133 = Ifges(5,1) * t308 + t171 + t522;
t132 = t424 + t521 + t523;
t125 = mrSges(7,1) * t172 - t426;
t123 = -mrSges(7,2) * t172 + t427;
t116 = -mrSges(7,1) * t161 + mrSges(7,2) * t162;
t91 = -mrSges(5,2) * t271 + mrSges(5,3) * t103;
t90 = mrSges(5,1) * t271 - mrSges(5,3) * t102;
t88 = pkin(5) * t404 - t483;
t49 = -qJ(6) * t404 + t61;
t43 = -pkin(5) * t307 - t183 * t270 + t60;
t35 = -mrSges(6,2) * t98 + mrSges(6,3) * t69;
t34 = -mrSges(7,2) * t98 + mrSges(7,3) * t69;
t32 = mrSges(7,1) * t98 - mrSges(7,3) * t68;
t30 = pkin(5) * t300 + t51;
t23 = -mrSges(6,1) * t69 + mrSges(6,2) * t68;
t10 = -qJD(5) * t61 + t339;
t9 = -t130 * t376 + t368;
t8 = -qJ(6) * t356 + (-qJD(5) * t130 + t304) * t281 + t368;
t7 = pkin(5) * t141 + t304 * t285 + (-t122 + (qJ(6) * t183 - t131) * t281) * qJD(5) + t339;
t2 = [t403 * t538 + t410 * t539 + t235 * (mrSges(4,1) * t215 + mrSges(4,2) * t214) + t212 * (Ifges(4,4) * t214 - Ifges(4,2) * t215) / 0.2e1 + qJD(3) * (Ifges(4,5) * t214 - Ifges(4,6) * t215) / 0.2e1 + (Ifges(4,1) * t214 - Ifges(4,4) * t215) * t451 + (-t299 * t542 - t300 * t527) * t462 + (-t299 * t529 - t300 * t542) * t460 + (-t392 * t531 - t511 * (t288 * t231 - t272 * t284) - t543 * t203 + t499 * t202 + t505 * t284 + (-t503 + t504) * t288) * g(2) + (-t543 * t201 + t499 * t200 + (t272 * t511 - t281 * t531 + t505) * t288 + (-m(6) * (-t231 - t480) - m(7) * (-t231 - t482) + m(5) * t231 + t503) * t284) * g(1) + (mrSges(4,2) * t234 - mrSges(4,3) * t139 + Ifges(4,1) * t180 + Ifges(4,4) * t181 + Ifges(4,5) * qJDD(3)) * t224 + (-t1 * t403 + t299 * t31 - t3 * t404 - t300 * t40) * mrSges(7,3) - (-m(5) * t29 + m(6) * t26 + t23 - t90) * t483 + m(5) * (t112 * t50 + t130 * t28 + t163 * t196 + t188 * t445) + (-mrSges(4,1) * t234 + mrSges(4,3) * t138 + Ifges(4,4) * t180 + Ifges(4,2) * t181 + Ifges(4,6) * qJDD(3)) * t223 - t214 * t415 + (-t523 / 0.2e1 - t521 / 0.2e1 - t132 / 0.2e1 - Ifges(5,4) * t454 + t497 * t462 + t498 * t460 + t526 * t458 - t501 + t496 / 0.2e1) * t141 + (t299 * t46 - t300 * t47 - t403 * t6 - t404 * t5) * mrSges(6,3) + m(3) * (-pkin(1) * t265 + (t252 + t373) * qJ(2) * t380) - t260 * t337 + 0.2e1 * t380 * t252 * mrSges(3,3) + (t163 * mrSges(5,2) - t29 * mrSges(5,3) + Ifges(5,1) * t102 + Ifges(5,4) * t103 + Ifges(5,5) * t271 + t12 * t318 + t26 * t319 + t342 * t494 + t465 * t517 + t466 * t516 + t467 * t515) * t183 + t196 * t340 + t136 * t445 + m(6) * (t10 * t46 + t47 * t9 + t5 * t61 + t6 * t60) + m(4) * (t138 * t190 + t139 * t189 + t166 * t185 + t167 * t184 - t234 * t260) + m(7) * (t1 * t43 + t12 * t88 + t3 * t49 + t30 * t73 + t31 * t7 + t40 * t8) + (-mrSges(5,1) * t163 - t1 * mrSges(7,1) + mrSges(5,3) * t28 + Ifges(5,4) * t102 + Ifges(5,2) * t103 + Ifges(5,6) * t271 - t465 * t526 - t466 * t497 - t467 * t498 + t502 - t513 / 0.2e1) * t307 + t214 * t174 / 0.2e1 - t215 * t173 / 0.2e1 + t166 * t191 + t167 * t192 + t189 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t180) + t190 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t181) + (-t356 / 0.2e1 - t411 / 0.2e1) * t495 + Ifges(2,3) * qJDD(1) + t50 * t164 + t130 * t91 + t8 * t123 + t9 * t124 + t7 * t125 + t10 * t126 + t30 * t116 + t88 * t22 + t60 * t33 + t61 * t35 + (-t299 * t498 - t300 * t497) * t458 + t49 * t34 - t525 * t404 / 0.2e1 + t43 * t32 + (t171 / 0.2e1 + t522 / 0.2e1 + t133 / 0.2e1 + Ifges(5,1) * t454 + t507) * t140 + (-m(5) * t111 + m(6) * t104 - t484) * t51 + t265 * t324 - pkin(1) * t325 + (Ifges(3,4) * t277 + Ifges(3,2) * t278) * t371 + (Ifges(3,1) * t277 + Ifges(3,4) * t278) * t372 - t215 * t414 + t73 * (mrSges(7,1) * t300 - mrSges(7,2) * t299) + t104 * (mrSges(6,1) * t300 - mrSges(6,2) * t299); t340 + m(3) * t265 + t325 + t337 - t212 * t191 + t213 * t192 - t333 * t164 + (-t116 + t484) * t308 + (t32 + t33 + t172 * (t123 + t124)) * t285 + (t34 + t35 - t172 * (t125 + t126)) * t281 + (-g(1) * t284 + g(2) * t288) * (m(3) + m(4) + t511) - t347 * t380 * qJD(1) ^ 2 + (t1 * t285 - t308 * t73 + t3 * t281 + t172 * (-t281 * t31 + t285 * t40)) * m(7) + (-t104 * t308 + t172 * t310 + t5 * t281 + t6 * t285) * m(6) + (t111 * t308 - t112 * t333 + t163) * m(5) + (t184 * t213 - t185 * t212 + t234) * m(4); (t284 * t479 + t476) * g(2) + (t288 * t479 + t477) * g(1) + (-m(5) * t256 - m(6) * (t256 + t480) - m(7) * (t256 + t482) - t323 + t473) * g(3) + (-Ifges(5,4) * t455 + t469 + t132 / 0.2e1) * t308 + (-t133 / 0.2e1 + Ifges(5,4) * t456 + t470) * t333 + (-t281 * t366 - t52) * t126 + (t366 - t115) * t164 - t235 * (mrSges(4,1) * t213 + mrSges(4,2) * t212) + t246 * t22 + (t26 * t263 + (t104 * t282 + t286 * t310) * qJD(4) * pkin(3) - g(1) * t243 - g(2) * t242 - t104 * t114 - t46 * t52 - t47 * t53) * m(6) + t489 * t116 + t500 - t213 * (Ifges(4,1) * t212 - t425) / 0.2e1 + (t111 * t114 - t112 * t115 - t188 * t446 + (t28 * t282 + t286 * t29 + (-t111 * t282 + t112 * t286) * qJD(4)) * pkin(3)) * m(5) + t471 * t261 - t136 * t446 + t492 * t125 + t493 * t123 + (t1 * t218 + t12 * t246 + t219 * t3 + t31 * t492 + t40 * t493 + t489 * t73) * m(7) - (-Ifges(4,2) * t213 + t174 + t206) * t212 / 0.2e1 + t218 * t32 + t219 * t34 + t263 * t23 - qJD(3) * (Ifges(4,5) * t212 - Ifges(4,6) * t213) / 0.2e1 - t184 * t191 + t185 * t192 + Ifges(4,5) * t180 + Ifges(4,6) * t181 - t138 * mrSges(4,2) + t139 * mrSges(4,1) + t520 * t124 + t519 * t484 + Ifges(4,3) * qJDD(3) + t496 * t455 + t35 * t394 + t213 * t414 + t212 * t415 + t90 * t442 + t91 * t443 + t173 * t451 + (m(5) * t444 + mrSges(4,1) * t266 + mrSges(4,2) * t267 + t320) * t475; (-t104 * t112 - t46 * t54 - t47 * t55 - g(1) * (-t288 * t441 + t243) - g(2) * (-t284 * t441 + t242) - pkin(4) * t26) * m(6) + t475 * t320 + (t284 * t478 + t476) * g(2) + (t288 * t478 + t477) * g(1) + t469 * t308 + t470 * t333 + t484 * t112 + (-t424 + t496) * t455 + t247 * t32 + t248 * t34 + (t171 + t133) * t456 + t488 * t116 + t490 * t125 + t500 + t471 * pkin(9) + t491 * t123 + (t1 * t247 - t12 * t262 + t248 * t3 + t31 * t490 + t40 * t491 + t488 * t73) * m(7) - t262 * t22 - t111 * t164 - t55 * t124 - t54 * t126 - pkin(4) * t23 + (t473 + t504) * g(3) + t35 * t439 + t132 * t454; (t281 * t487 + t319 + t430) * g(3) * t258 + (t161 * t529 - t535) * t461 + (-t162 * t527 + t494 + t536) * t463 + t495 * t460 + (t161 * t498 - t162 * t497) * t459 + (t200 * t474 + t201 * t499) * g(2) + (-t202 * t474 - t203 * t499) * g(1) + t487 * t1 + (t126 + t428) * t47 - t502 + (-m(7) * (-t31 + t39) + t125 + t426) * t40 - t104 * (mrSges(6,1) * t162 + mrSges(6,2) * t161) - t73 * (mrSges(7,1) * t162 + mrSges(7,2) * t161) - t39 * t123 + (t32 + (-m(7) * t73 - t116) * t162) * pkin(5) + (-t124 + t429) * t46 + t31 * t427 + t513; -t161 * t123 + t162 * t125 + (g(3) * t259 - t40 * t161 + t162 * t31 - t258 * t475 + t12) * m(7) + t22;];
tau  = t2;
