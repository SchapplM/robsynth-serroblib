% Calculate vector of inverse dynamics joint torques for
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:53
% EndTime: 2019-03-09 09:59:37
% DurationCPUTime: 30.31s
% Computational Cost: add. (7674->766), mult. (15965->976), div. (0->0), fcn. (9908->10), ass. (0->347)
t268 = sin(qJ(4));
t271 = cos(qJ(4));
t272 = cos(qJ(2));
t375 = qJD(1) * t272;
t187 = -qJD(2) * t268 - t271 * t375;
t266 = sin(pkin(9));
t373 = qJD(2) * t271;
t293 = t268 * t375 - t373;
t402 = cos(pkin(9));
t118 = -t402 * t187 - t266 * t293;
t246 = pkin(7) * t375;
t194 = pkin(3) * t375 + t246;
t265 = qJD(2) * qJ(3);
t164 = t265 + t194;
t121 = -pkin(4) * t187 + qJD(5) + t164;
t292 = t266 * t187 - t293 * t402;
t36 = pkin(5) * t118 - qJ(6) * t292 + t121;
t269 = sin(qJ(2));
t376 = qJD(1) * t269;
t231 = qJD(4) + t376;
t441 = t231 / 0.2e1;
t442 = -t231 / 0.2e1;
t450 = t292 / 0.2e1;
t451 = -t292 / 0.2e1;
t454 = -t118 / 0.2e1;
t535 = -t36 * mrSges(7,1) + Ifges(6,4) * t450 + Ifges(7,5) * t451 + Ifges(6,6) * t441 + Ifges(7,6) * t442 + (Ifges(6,2) + Ifges(7,3)) * t454;
t479 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t475 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t503 = Ifges(7,2) + Ifges(6,3);
t529 = Ifges(5,3) + t503;
t333 = t402 * t271;
t323 = t269 * t333;
t396 = t266 * t268;
t135 = qJD(1) * t323 - t376 * t396;
t370 = qJD(4) * t268;
t157 = -qJD(4) * t333 + t266 * t370;
t526 = t135 - t157;
t291 = -t266 * t271 - t268 * t402;
t286 = t269 * t291;
t136 = qJD(1) * t286;
t158 = t291 * qJD(4);
t525 = t136 + t158;
t453 = t118 / 0.2e1;
t502 = Ifges(6,6) - Ifges(7,6);
t506 = Ifges(6,4) - Ifges(7,5);
t534 = Ifges(6,2) * t453 - Ifges(7,3) * t454 + t502 * t442 + t506 * t451 - t535;
t459 = -m(6) - m(7);
t533 = mrSges(7,2) + mrSges(6,3);
t508 = Ifges(6,1) + Ifges(7,1);
t507 = -Ifges(4,4) + Ifges(3,5);
t505 = Ifges(7,4) + Ifges(6,5);
t504 = Ifges(4,5) - Ifges(3,6);
t316 = t272 * mrSges(4,2) - t269 * mrSges(4,3);
t321 = mrSges(3,1) * t272 - mrSges(3,2) * t269;
t532 = -t321 + t316;
t509 = t36 * mrSges(7,3);
t531 = -Ifges(6,4) * t454 - Ifges(7,5) * t453 - t505 * t441 - t508 * t450 + t509;
t244 = pkin(7) * t376;
t530 = qJD(3) + t244;
t261 = qJ(4) + pkin(9);
t250 = sin(t261);
t251 = cos(t261);
t318 = t268 * mrSges(5,1) + t271 * mrSges(5,2);
t528 = -t250 * t479 + t251 * t475 - t318;
t365 = qJD(1) * qJD(2);
t197 = -t272 * qJDD(1) + t269 * t365;
t113 = qJD(4) * t187 + qJDD(2) * t271 + t197 * t268;
t114 = qJD(4) * t293 - qJDD(2) * t268 + t197 * t271;
t62 = t113 * t266 - t114 * t402;
t462 = t62 / 0.2e1;
t63 = t113 * t402 + t266 * t114;
t461 = t63 / 0.2e1;
t198 = qJDD(1) * t269 + t272 * t365;
t184 = qJDD(4) + t198;
t527 = -t184 / 0.2e1;
t446 = t184 / 0.2e1;
t437 = pkin(4) * t271;
t240 = pkin(3) + t437;
t369 = qJD(4) * t271;
t491 = pkin(4) * t369 + t240 * t376 + t530;
t185 = -t333 + t396;
t253 = t269 * qJ(3);
t334 = -pkin(1) - t253;
t458 = pkin(2) + pkin(8);
t138 = (-t272 * t458 + t334) * qJD(1);
t193 = -pkin(3) * t376 - t244;
t484 = -t193 + qJD(3);
t145 = -qJD(2) * t458 + t484;
t87 = t138 * t271 + t145 * t268;
t76 = qJ(5) * t187 + t87;
t410 = t266 * t76;
t86 = -t138 * t268 + t271 * t145;
t75 = qJ(5) * t293 + t86;
t71 = pkin(4) * t231 + t75;
t24 = t402 * t71 - t410;
t73 = t402 * t76;
t25 = t266 * t71 + t73;
t183 = t198 * pkin(7);
t335 = qJDD(3) + t183;
t115 = pkin(3) * t198 - qJDD(2) * t458 + t335;
t371 = qJD(3) * t269;
t401 = qJDD(1) * pkin(1);
t284 = -qJ(3) * t198 - qJD(1) * t371 - t401;
t84 = t197 * t458 + t284;
t23 = -qJD(4) * t87 + t271 * t115 - t268 * t84;
t11 = pkin(4) * t184 - qJ(5) * t113 + qJD(5) * t293 + t23;
t22 = t268 * t115 - t138 * t370 + t145 * t369 + t271 * t84;
t13 = qJ(5) * t114 + qJD(5) * t187 + t22;
t3 = t11 * t402 - t266 * t13;
t4 = t266 * t11 + t402 * t13;
t523 = t185 * t3 - t24 * t525 - t25 * t526 + t291 * t4;
t1 = qJ(6) * t184 + qJD(6) * t231 + t4;
t2 = -t184 * pkin(5) + qJDD(6) - t3;
t20 = -t231 * pkin(5) + qJD(6) - t24;
t21 = qJ(6) * t231 + t25;
t522 = t1 * t291 - t185 * t2 + t20 * t525 - t21 * t526;
t521 = -t272 * t533 + t532;
t520 = -t121 * mrSges(6,1) + mrSges(7,2) * t21 + mrSges(6,3) * t25;
t463 = -t62 / 0.2e1;
t182 = t197 * pkin(7);
t139 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t182;
t116 = -pkin(3) * t197 - t139;
t68 = -pkin(4) * t114 + qJDD(5) + t116;
t5 = pkin(5) * t62 - qJ(6) * t63 - qJD(6) * t292 + t68;
t519 = 0.2e1 * Ifges(7,3) * t462 + t5 * mrSges(7,1) + t68 * mrSges(6,1) + Ifges(7,6) * t446 - t63 * Ifges(6,4) / 0.2e1 + Ifges(6,6) * t527 + (t462 - t463) * Ifges(6,2) + (-t506 + Ifges(7,5)) * t461;
t424 = mrSges(5,3) * t187;
t130 = -mrSges(5,2) * t231 + t424;
t423 = mrSges(5,3) * t293;
t131 = mrSges(5,1) * t231 + t423;
t302 = t271 * t130 - t268 * t131;
t81 = mrSges(5,1) * t184 - mrSges(5,3) * t113;
t82 = -mrSges(5,2) * t184 + mrSges(5,3) * t114;
t518 = t302 * qJD(4) + t268 * t82 + t271 * t81;
t498 = -t118 * t506 + t231 * t505 + t292 * t508;
t516 = mrSges(7,2) * t20 + t121 * mrSges(6,2) + t498 / 0.2e1 - mrSges(6,3) * t24;
t417 = Ifges(4,6) * t272;
t306 = -t269 * Ifges(4,2) - t417;
t514 = t20 * mrSges(7,1) + t25 * mrSges(6,2) + Ifges(4,4) * qJD(2) / 0.2e1 + qJD(1) * t306 / 0.2e1 - t21 * mrSges(7,3) - t24 * mrSges(6,1);
t513 = -m(5) - m(4);
t501 = t184 * t505 - t506 * t62 + t508 * t63;
t38 = -mrSges(7,2) * t62 + mrSges(7,3) * t184;
t39 = -mrSges(6,2) * t184 - mrSges(6,3) * t62;
t500 = t38 + t39;
t40 = mrSges(6,1) * t184 - mrSges(6,3) * t63;
t41 = -t184 * mrSges(7,1) + t63 * mrSges(7,2);
t499 = t41 - t40;
t497 = pkin(5) * t526 - qJ(6) * t525 + qJD(6) * t185 + t491;
t92 = -mrSges(6,2) * t231 - mrSges(6,3) * t118;
t93 = -mrSges(7,2) * t118 + mrSges(7,3) * t231;
t426 = t92 + t93;
t94 = mrSges(6,1) * t231 - mrSges(6,3) * t292;
t95 = -mrSges(7,1) * t231 + mrSges(7,2) * t292;
t425 = t95 - t94;
t270 = sin(qJ(1));
t273 = cos(qJ(1));
t485 = g(1) * t273 + g(2) * t270;
t495 = t272 * t485;
t258 = t272 * pkin(2);
t378 = t258 + t253;
t324 = pkin(8) * t272 + t378;
t179 = -pkin(1) - t324;
t457 = pkin(3) + pkin(7);
t216 = t457 * t269;
t123 = t271 * t179 + t268 * t216;
t126 = -mrSges(5,1) * t187 - mrSges(5,2) * t293;
t353 = mrSges(4,1) * t375;
t208 = -qJD(2) * mrSges(4,3) - t353;
t494 = -t208 + t126;
t493 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t375 - t208;
t354 = mrSges(4,1) * t376;
t492 = -mrSges(3,3) * t376 - t354 + (mrSges(3,1) - mrSges(4,2)) * qJD(2);
t418 = Ifges(4,6) * t269;
t490 = t269 * (-Ifges(4,2) * t272 + t418) + t272 * (Ifges(4,3) * t269 - t417);
t489 = t269 * t504 + t272 * t507;
t488 = -t182 * t272 + t183 * t269;
t149 = -qJDD(2) * pkin(2) + t335;
t487 = -t139 * t272 + t149 * t269;
t482 = Ifges(5,5) * t113 + Ifges(5,6) * t114 + t184 * t529 - t502 * t62 + t505 * t63;
t412 = t293 * Ifges(5,4);
t102 = t187 * Ifges(5,2) + t231 * Ifges(5,6) - t412;
t180 = Ifges(5,4) * t187;
t103 = -Ifges(5,1) * t293 + t231 * Ifges(5,5) + t180;
t305 = -t272 * Ifges(4,3) - t418;
t481 = Ifges(4,5) * qJD(2) + qJD(1) * t305 + t271 * t102 + t268 * t103;
t243 = Ifges(3,4) * t375;
t480 = Ifges(3,1) * t376 + Ifges(3,5) * qJD(2) - Ifges(5,5) * t293 + t187 * Ifges(5,6) - t118 * t502 + t231 * t529 + t292 * t505 + t243;
t326 = -m(5) * t458 - mrSges(5,3);
t405 = t272 * mrSges(4,3);
t478 = -t405 + t528 * t272 + (-mrSges(4,2) - t326 + (m(4) - t459) * pkin(2) + t533) * t269;
t474 = -m(5) * pkin(3) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t469 = Ifges(6,2) * t454 - Ifges(7,3) * t453 + t450 * t506 + t535;
t468 = Ifges(6,4) * t453 + Ifges(7,5) * t454 + t442 * t505 + t451 * t508 + t509;
t466 = -t113 * Ifges(5,4) / 0.2e1 - t114 * Ifges(5,2) / 0.2e1 + Ifges(5,6) * t527;
t456 = t113 / 0.2e1;
t455 = t114 / 0.2e1;
t443 = -t293 / 0.2e1;
t439 = pkin(4) * t293;
t438 = pkin(4) * t266;
t436 = pkin(7) * t269;
t433 = g(3) * t272;
t256 = t272 * pkin(7);
t304 = pkin(8) * t269 - qJ(3) * t272;
t374 = qJD(2) * t269;
t328 = pkin(2) * t374 - t371;
t132 = qJD(2) * t304 + t328;
t372 = qJD(2) * t272;
t196 = t457 * t372;
t173 = t271 * t196;
t329 = qJ(5) * t272 - t179;
t32 = pkin(4) * t372 + t173 + t329 * t369 + (-qJ(5) * t374 - qJD(4) * t216 + qJD(5) * t272 - t132) * t268;
t368 = qJD(4) * t272;
t343 = t268 * t368;
t289 = t269 * t373 + t343;
t366 = t271 * qJD(5);
t65 = t271 * t132 - t179 * t370 + t268 * t196 + t216 * t369;
t37 = qJ(5) * t289 - t272 * t366 + t65;
t10 = t266 * t32 + t402 * t37;
t245 = pkin(2) * t376;
t150 = qJD(1) * t304 + t245;
t109 = -t150 * t268 + t271 * t194;
t393 = t268 * t269;
t83 = (pkin(4) * t272 - qJ(5) * t393) * qJD(1) + t109;
t110 = t271 * t150 + t268 * t194;
t89 = qJ(5) * t271 * t376 + t110;
t34 = t266 * t83 + t402 * t89;
t383 = t271 * t272;
t104 = -qJ(5) * t383 + t123;
t190 = t271 * t216;
t96 = pkin(4) * t269 + t268 * t329 + t190;
t47 = t402 * t104 + t266 * t96;
t422 = Ifges(3,4) * t269;
t421 = Ifges(3,4) * t272;
t420 = Ifges(5,4) * t268;
t419 = Ifges(5,4) * t271;
t411 = t22 * t268;
t408 = t271 * mrSges(5,3);
t392 = t268 * t270;
t391 = t268 * t272;
t390 = t268 * t273;
t389 = t269 * t270;
t388 = t269 * t273;
t387 = t270 * t271;
t386 = t270 * t272;
t382 = t271 * t273;
t381 = t272 * t273;
t239 = t268 * pkin(4) + qJ(3);
t380 = qJ(5) + t458;
t217 = t272 * pkin(3) + t256;
t377 = t273 * pkin(1) + t270 * pkin(7);
t360 = pkin(4) * t392;
t232 = pkin(4) * t393;
t234 = pkin(4) * t390;
t355 = -t459 - t513;
t352 = t268 * t388;
t351 = t269 * t387;
t350 = t269 * t382;
t233 = pkin(4) * t383;
t159 = t233 + t217;
t345 = t402 * pkin(4);
t344 = t268 * t374;
t19 = t62 * mrSges(6,1) + t63 * mrSges(6,2);
t18 = t62 * mrSges(7,1) - t63 * mrSges(7,3);
t336 = -t369 / 0.2e1;
t332 = -t365 / 0.2e1;
t153 = t198 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t330 = t380 * t271;
t327 = pkin(2) * t381 + qJ(3) * t388 + t377;
t320 = mrSges(3,1) * t269 + mrSges(3,2) * t272;
t319 = mrSges(5,1) * t271 - mrSges(5,2) * t268;
t315 = Ifges(5,1) * t271 - t420;
t314 = Ifges(5,1) * t268 + t419;
t313 = t272 * Ifges(3,2) + t422;
t311 = -Ifges(5,2) * t268 + t419;
t310 = Ifges(5,2) * t271 + t420;
t308 = Ifges(5,5) * t271 - Ifges(5,6) * t268;
t307 = Ifges(5,5) * t268 + Ifges(5,6) * t271;
t303 = t86 * t268 - t87 * t271;
t201 = -qJD(2) * pkin(2) + t530;
t206 = -t246 - t265;
t301 = t201 * t272 + t206 * t269;
t300 = t334 - t258;
t298 = pkin(1) * t320;
t9 = -t266 * t37 + t32 * t402;
t33 = -t266 * t89 + t402 * t83;
t165 = t300 * qJD(1);
t297 = t165 * (-mrSges(4,2) * t269 - t405);
t296 = t269 * (Ifges(3,1) * t272 - t422);
t46 = -t266 * t104 + t402 * t96;
t288 = -t271 * t368 + t344;
t287 = t370 * t380 - t366;
t283 = Ifges(5,5) * t272 + t269 * t314;
t282 = Ifges(5,6) * t272 + t269 * t310;
t281 = Ifges(5,3) * t272 + t269 * t307;
t127 = -pkin(4) * t343 + (-pkin(7) - t240) * t374;
t277 = -qJD(4) * t303 + t23 * t271 + t411;
t267 = -qJ(5) - pkin(8);
t259 = t273 * pkin(7);
t238 = -t345 - pkin(5);
t236 = qJ(6) + t438;
t228 = qJ(3) * t381;
t227 = qJ(3) * t386;
t202 = -pkin(1) - t378;
t199 = t380 * t268;
t195 = t457 * t374;
t192 = -qJ(3) * t375 + t245;
t191 = t316 * qJD(1);
t177 = t319 * t272;
t169 = -t268 * t389 + t382;
t168 = t351 + t390;
t167 = t352 + t387;
t166 = t350 - t392;
t160 = Ifges(3,6) * qJD(2) + qJD(1) * t313;
t154 = -qJ(3) * t372 + t328;
t152 = mrSges(4,1) * t197 - qJDD(2) * mrSges(4,3);
t148 = t291 * t272;
t147 = t266 * t391 - t272 * t333;
t146 = -qJD(4) * t330 - t268 * qJD(5);
t144 = -t250 * t389 + t251 * t273;
t143 = t250 * t273 + t251 * t389;
t142 = t250 * t388 + t251 * t270;
t141 = t250 * t270 - t251 * t388;
t125 = -t199 * t402 - t266 * t330;
t124 = -t199 * t266 + t330 * t402;
t122 = -t179 * t268 + t190;
t112 = -pkin(5) * t291 + qJ(6) * t185 + t239;
t111 = pkin(2) * t197 + t284;
t98 = -qJD(2) * t286 + t185 * t368;
t97 = qJD(2) * t323 - t158 * t272 - t266 * t344;
t91 = t146 * t402 + t266 * t287;
t90 = t146 * t266 - t287 * t402;
t78 = -pkin(5) * t147 - qJ(6) * t148 + t159;
t70 = mrSges(6,1) * t118 + mrSges(6,2) * t292;
t69 = mrSges(7,1) * t118 - mrSges(7,3) * t292;
t66 = -qJD(4) * t123 - t132 * t268 + t173;
t64 = -mrSges(5,1) * t114 + mrSges(5,2) * t113;
t48 = pkin(5) * t292 + qJ(6) * t118 - t439;
t45 = t113 * Ifges(5,1) + t114 * Ifges(5,4) + t184 * Ifges(5,5);
t43 = -t269 * pkin(5) - t46;
t42 = qJ(6) * t269 + t47;
t30 = -pkin(5) * t375 - t33;
t29 = qJ(6) * t375 + t34;
t28 = t402 * t75 - t410;
t27 = t266 * t75 + t73;
t26 = -pkin(5) * t97 - qJ(6) * t98 - qJD(6) * t148 + t127;
t7 = -pkin(5) * t372 - t9;
t6 = qJ(6) * t372 + qJD(6) * t269 + t10;
t8 = [(t1 * t269 - t148 * t5) * mrSges(7,3) + (t148 * t68 - t269 * t4) * mrSges(6,2) + (t516 - t531) * t98 + (-t492 * pkin(7) + Ifges(6,6) * t454 + Ifges(7,6) * t453 + t86 * mrSges(5,1) - t87 * mrSges(5,2) + t480 / 0.2e1 + t201 * mrSges(4,1) + t503 * t441 + t505 * t450 - t514) * t372 + t321 * t401 + (-qJDD(2) * mrSges(3,2) - t152) * t256 + (-t493 * pkin(7) - t160 / 0.2e1 + t481 / 0.2e1 + t206 * mrSges(4,1)) * t374 + (t469 + t520) * t97 + (t272 * (-Ifges(3,2) * t269 + t421) + t296) * t365 / 0.2e1 + (-t169 * mrSges(5,1) + t168 * mrSges(5,2) + t459 * (t273 * t240 + t267 * t386 + t259) - t479 * t144 - t475 * t143 + (-m(3) + t513) * t259 + t474 * t273 + (m(3) * pkin(1) - m(4) * t300 - m(5) * t334 - t272 * t326 + mrSges(2,1) + t459 * (t300 - t232) - t521) * t270) * g(1) - t198 * t306 / 0.2e1 - t197 * t313 / 0.2e1 + t197 * t305 / 0.2e1 + t164 * (-mrSges(5,1) * t289 + mrSges(5,2) * t288) + m(5) * (t116 * t217 + t122 * t23 + t123 * t22 - t164 * t195 + t65 * t87 + t66 * t86) + (Ifges(6,4) * t148 + Ifges(6,6) * t269) * t463 + t383 * t466 + (Ifges(5,6) * t269 - t272 * t310) * t455 + (Ifges(5,5) * t269 - t272 * t314) * t456 + (qJD(2) * t283 - t315 * t368) * t443 + qJD(2) * t297 + t111 * t316 - t86 * mrSges(5,3) * t288 + t487 * mrSges(4,1) + t198 * t269 * Ifges(3,1) + m(7) * (t1 * t42 + t2 * t43 + t20 * t7 + t21 * t6 + t26 * t36 + t5 * t78) + m(6) * (t10 * t25 + t121 * t127 + t159 * t68 + t24 * t9 + t3 * t46 + t4 * t47) + (Ifges(7,5) * t148 + Ifges(7,6) * t269) * t462 + (-qJDD(2) * mrSges(3,1) + t153) * t436 + t198 * t421 / 0.2e1 + (-Ifges(3,4) * t197 + Ifges(3,5) * qJDD(2) + t482) * t269 / 0.2e1 + t87 * mrSges(5,3) * t289 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) - t519) * t147 + t272 * (Ifges(3,4) * t198 - Ifges(3,2) * t197 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t272 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t198 + Ifges(4,3) * t197) / 0.2e1 + t2 * (-t269 * mrSges(7,1) + mrSges(7,2) * t148) + t3 * (t269 * mrSges(6,1) - mrSges(6,3) * t148) - t269 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t198 + Ifges(4,6) * t197) / 0.2e1 + Ifges(2,3) * qJDD(1) + t217 * t64 - pkin(1) * (mrSges(3,1) * t197 + mrSges(3,2) * t198) + t202 * (-mrSges(4,2) * t197 - mrSges(4,3) * t198) - t195 * t126 + t154 * t191 + t116 * t177 + t159 * t19 + t65 * t130 + t66 * t131 + t127 * t70 + t122 * t81 + t123 * t82 - t45 * t391 / 0.2e1 + t23 * (mrSges(5,1) * t269 + mrSges(5,3) * t391) + t7 * t95 + t22 * (-mrSges(5,2) * t269 - mrSges(5,3) * t383) + t10 * t92 + t6 * t93 + t9 * t94 + t78 * t18 + t26 * t69 + t187 * (qJD(2) * t282 - t311 * t368) / 0.2e1 + t42 * t38 + t43 * t41 + t46 * t40 + t47 * t39 - t298 * t365 + (-m(3) * t377 - t167 * mrSges(5,1) - t166 * mrSges(5,2) + t513 * t327 + t459 * (pkin(4) * t352 + t270 * t240 - t267 * t381 + t327) - t479 * t142 - t475 * t141 + (-m(5) * pkin(8) - mrSges(5,3) - t533) * t381 + (-mrSges(2,1) + t532) * t273 + t474 * t270) * g(2) + t272 * t103 * t336 + m(4) * (t111 * t202 + t154 * t165 + (qJD(2) * t301 + t487) * pkin(7)) + (-t197 * t256 + t198 * t436 + t488) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t488) + t489 * qJD(2) ^ 2 / 0.2e1 + t490 * t332 + (t502 * t147 + t505 * t148 + t269 * t529 - t272 * t307) * t446 + t501 * t148 / 0.2e1 + (t269 * t507 - t272 * t504) * qJDD(2) / 0.2e1 + (t148 * t508 + t269 * t505) * t461 + (qJD(2) * t281 - t308 * t368 + t502 * t97) * t441 + t102 * t343 / 0.2e1; -t531 * t158 + (Ifges(6,6) * t453 + Ifges(7,6) * t454 + t503 * t442 + t505 * t451 + t514) * t375 + (-m(5) * t277 - t518) * t458 + (-m(4) * t378 - m(5) * t324 - t272 * mrSges(5,3) + t459 * (-t267 * t272 + t232 + t378) + t528 * t269 + t521) * g(3) + (t298 - t296 / 0.2e1 + t490 / 0.2e1) * qJD(1) ^ 2 + (t291 * t502 + t308) * t446 - (t187 * t310 + t231 * t307 - t293 * t314) * qJD(4) / 0.2e1 - (t187 * t282 + t231 * t281 - t283 * t293) * qJD(1) / 0.2e1 + (t158 / 0.2e1 + t136 / 0.2e1) * t498 + (-t501 / 0.2e1 - t505 * t446 - t68 * mrSges(6,2) + t5 * mrSges(7,3) - Ifges(6,4) * t463 - Ifges(7,5) * t462 - t508 * t461) * t185 + (mrSges(6,1) * t526 + mrSges(6,2) * t525) * t121 + t522 * mrSges(7,2) + t523 * mrSges(6,3) + t231 * t319 * t164 + (-pkin(2) * t149 - qJ(3) * t139 - qJD(3) * t206 - t165 * t192) * m(4) + (t116 * qJ(3) - t109 * t86 - t110 * t87 + t164 * t484) * m(5) + (Ifges(3,3) + Ifges(4,1)) * qJDD(2) + t116 * t318 - t468 * t136 + t268 * t466 + t311 * t455 + t315 * t456 + (-t152 + t64) * qJ(3) - t519 * t291 + (-m(4) * t301 * pkin(7) - t297 - t87 * (-mrSges(5,2) * t272 + t269 * t408) - t86 * (mrSges(5,1) * t272 - mrSges(5,3) * t393)) * qJD(1) + (-t369 * t87 + t370 * t86 - t411) * mrSges(5,3) + t271 * t45 / 0.2e1 + t239 * t19 - t193 * t126 - t192 * t191 - t183 * mrSges(3,1) + t182 * mrSges(3,2) + t149 * mrSges(4,2) - pkin(2) * t153 - t23 * t408 - t139 * mrSges(4,3) - t110 * t130 - t109 * t131 + t112 * t18 - t30 * t95 - t34 * t92 - t29 * t93 - t33 * t94 + t160 * t376 / 0.2e1 - t103 * t370 / 0.2e1 - t201 * t353 - t206 * t354 - (-Ifges(3,2) * t376 + t243 + t480) * t375 / 0.2e1 - t481 * t376 / 0.2e1 + t485 * t320 + t489 * t332 + (-t124 * t3 + t125 * t4 + t239 * t68 + (t91 - t34) * t25 + (-t90 - t33) * t24 + t491 * t121) * m(6) + t491 * t70 + t492 * t246 + t493 * t244 + t494 * qJD(3) + t534 * t135 + t102 * t336 + t425 * t90 + t426 * t91 + t497 * t69 + (t1 * t125 + t112 * t5 + t124 * t2 + t497 * t36 + (t91 - t29) * t21 + (-t30 + t90) * t20) * m(7) + t499 * t124 + t500 * t125 + (t502 * t441 + t469) * t157 + t504 * t197 + t507 * t198 + (t459 * (t272 * t234 + t267 * t388 + t228) + t513 * t228 + t478 * t273) * g(1) + (t459 * (t267 * t389 + t272 * t360 + t227) + t513 * t227 + t478 * t270) * g(2); (-t69 - t70 - t494) * qJD(2) - t500 * t291 + t499 * t185 + t355 * t433 + ((t191 + t302) * qJD(1) - t485 * t355) * t269 + t153 + t526 * t426 - t525 * t425 + (-qJD(2) * t36 - t522) * m(7) + (-qJD(2) * t121 - t523) * m(6) + (-qJD(2) * t164 - t303 * t376 + t277) * m(5) + (qJD(2) * t206 + t165 * t376 + t149) * m(4) + t518; -(Ifges(5,2) * t293 + t103 + t180) * t187 / 0.2e1 + (Ifges(5,5) * t187 + Ifges(5,6) * t293) * t442 - t164 * (-mrSges(5,1) * t293 + mrSges(5,2) * t187) + t293 * (Ifges(5,1) * t187 + t412) / 0.2e1 + (-t20 * t27 - t36 * t48 + t1 * t236 + t2 * t238 - g(3) * (-t233 + (-pkin(5) * t251 - qJ(6) * t250) * t272) + (qJD(6) - t28) * t21) * m(7) + (t424 - t130) * t86 + (-(-t251 * mrSges(7,1) - t250 * mrSges(7,3)) * t272 + t177) * g(3) + (-t468 + t516) * t118 + t39 * t438 + t102 * t443 + t40 * t345 - (-mrSges(6,1) * t251 + mrSges(6,2) * t250) * t433 + t482 + (-t423 + t131) * t87 + ((t266 * t4 + t3 * t402) * pkin(4) + t437 * t433 + t121 * t439 + t24 * t27 - t25 * t28) * m(6) + t236 * t38 + t238 * t41 + t70 * t439 + qJD(6) * t93 - t48 * t69 - t22 * mrSges(5,2) + t23 * mrSges(5,1) - t4 * mrSges(6,2) - t2 * mrSges(7,1) + t3 * mrSges(6,1) + t1 * mrSges(7,3) + (t520 - t534) * t292 - t425 * t27 - t426 * t28 + (-mrSges(5,1) * t166 + mrSges(5,2) * t167 + t459 * (pkin(4) * t350 - t360) - t475 * t142 + t479 * t141) * g(1) + (-mrSges(5,1) * t168 - mrSges(5,2) * t169 + t459 * (pkin(4) * t351 + t234) + t475 * t144 - t479 * t143) * g(2); t459 * t269 * g(3) - t425 * t292 + t426 * t118 + t18 + t19 + (t118 * t21 - t20 * t292 - t495 + t5) * m(7) + (t118 * t25 + t24 * t292 - t495 + t68) * m(6); t292 * t69 - t231 * t93 + (-g(1) * t141 + g(2) * t143 - t21 * t231 - t251 * t433 + t292 * t36 + t2) * m(7) + t41;];
tau  = t8;
