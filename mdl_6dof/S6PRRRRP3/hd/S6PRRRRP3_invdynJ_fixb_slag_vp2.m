% Calculate vector of inverse dynamics joint torques for
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:24
% EndTime: 2019-03-09 00:07:12
% DurationCPUTime: 27.45s
% Computational Cost: add. (10238->788), mult. (23299->1057), div. (0->0), fcn. (17400->14), ass. (0->350)
t536 = Ifges(6,4) + Ifges(7,4);
t306 = sin(qJ(2));
t301 = sin(pkin(6));
t405 = qJD(1) * t301;
t377 = t306 * t405;
t259 = qJD(2) * pkin(8) + t377;
t305 = sin(qJ(3));
t309 = cos(qJ(3));
t302 = cos(pkin(6));
t404 = qJD(1) * t302;
t195 = -t305 * t259 + t309 * t404;
t349 = pkin(3) * t305 - pkin(9) * t309;
t252 = t349 * qJD(2);
t304 = sin(qJ(4));
t308 = cos(qJ(4));
t135 = -t195 * t304 + t308 * t252;
t413 = t308 * t309;
t332 = pkin(4) * t305 - pkin(10) * t413;
t311 = -pkin(10) - pkin(9);
t378 = qJD(4) * t311;
t556 = -qJD(2) * t332 + t308 * t378 - t135;
t136 = t308 * t195 + t304 * t252;
t400 = qJD(2) * t309;
t373 = t304 * t400;
t555 = -pkin(10) * t373 - t304 * t378 + t136;
t300 = qJ(4) + qJ(5);
t297 = cos(t300);
t453 = pkin(4) * t308;
t262 = pkin(5) * t297 + t453;
t293 = pkin(3) + t453;
t347 = -mrSges(5,1) * t308 + mrSges(5,2) * t304;
t554 = -m(6) * t293 - m(7) * (pkin(3) + t262) - m(5) * pkin(3) + t347;
t553 = m(6) * t311 + m(7) * (-qJ(6) + t311) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(9) - mrSges(5,3);
t539 = -mrSges(7,1) - mrSges(6,1);
t538 = -mrSges(7,2) - mrSges(6,2);
t537 = Ifges(6,1) + Ifges(7,1);
t535 = -Ifges(7,5) - Ifges(6,5);
t534 = Ifges(6,2) + Ifges(7,2);
t533 = Ifges(6,6) + Ifges(7,6);
t532 = Ifges(7,3) + Ifges(6,3);
t303 = sin(qJ(5));
t307 = cos(qJ(5));
t333 = t303 * t304 - t307 * t308;
t504 = qJD(4) + qJD(5);
t167 = t504 * t333;
t326 = t333 * t309;
t204 = qJD(2) * t326;
t552 = t167 - t204;
t246 = t303 * t308 + t304 * t307;
t168 = t504 * t246;
t327 = t246 * t309;
t203 = qJD(2) * t327;
t551 = t168 - t203;
t286 = qJD(4) - t400;
t276 = qJD(5) + t286;
t398 = qJD(3) * t308;
t402 = qJD(2) * t305;
t243 = -t304 * t402 + t398;
t244 = qJD(3) * t304 + t308 * t402;
t351 = t307 * t243 - t244 * t303;
t160 = t243 * t303 + t244 * t307;
t544 = t536 * t160;
t527 = t276 * t533 + t351 * t534 + t544;
t550 = t527 / 0.2e1;
t391 = qJD(2) * qJD(3);
t256 = qJDD(2) * t309 - t305 * t391;
t240 = qJDD(4) - t256;
t234 = qJDD(5) + t240;
t257 = qJDD(2) * t305 + t309 * t391;
t147 = qJD(4) * t243 + qJDD(3) * t304 + t257 * t308;
t148 = -qJD(4) * t244 + qJDD(3) * t308 - t257 * t304;
t49 = qJD(5) * t351 + t147 * t307 + t148 * t303;
t50 = -qJD(5) * t160 - t147 * t303 + t148 * t307;
t549 = -t534 * t50 / 0.2e1 - t536 * t49 / 0.2e1 - t533 * t234 / 0.2e1;
t270 = t311 * t304;
t271 = t311 * t308;
t392 = qJD(5) * t307;
t393 = qJD(5) * t303;
t525 = t270 * t392 + t271 * t393 + t303 * t556 - t307 * t555;
t186 = t303 * t270 - t307 * t271;
t524 = -qJD(5) * t186 + t303 * t555 + t307 * t556;
t548 = qJ(6) * t160;
t547 = t536 * t351;
t310 = cos(qJ(2));
t403 = qJD(2) * t301;
t367 = qJD(1) * t403;
t273 = t310 * t367;
t390 = qJDD(1) * t301;
t211 = t306 * t390 + t273;
t546 = qJDD(2) * pkin(8) + qJD(3) * t404 + t211;
t296 = sin(t300);
t455 = pkin(4) * t304;
t261 = pkin(5) * t296 + t455;
t346 = t304 * mrSges(5,1) + t308 * mrSges(5,2);
t545 = -m(7) * t261 + mrSges(3,2) - mrSges(4,3) - t346;
t394 = qJD(4) * t308;
t397 = qJD(3) * t309;
t320 = t304 * t397 + t305 * t394;
t196 = t309 * t259 + t305 * t404;
t396 = qJD(4) * t304;
t513 = -t196 + (-t373 + t396) * pkin(4);
t348 = mrSges(4,1) * t309 - mrSges(4,2) * t305;
t543 = t305 * t553 + t309 * t554 - mrSges(3,1) - t348;
t482 = t147 / 0.2e1;
t481 = t148 / 0.2e1;
t467 = t240 / 0.2e1;
t542 = t256 / 0.2e1;
t541 = t257 / 0.2e1;
t479 = -t351 / 0.2e1;
t531 = -t234 * t535 + t49 * t537 + t50 * t536;
t40 = -mrSges(7,2) * t234 + mrSges(7,3) * t50;
t41 = -mrSges(6,2) * t234 + mrSges(6,3) * t50;
t530 = t40 + t41;
t529 = -qJ(6) * t551 - qJD(6) * t333 + t525;
t528 = -pkin(5) * t402 + qJ(6) * t552 - qJD(6) * t246 + t524;
t526 = t160 * t537 - t276 * t535 + t547;
t69 = -mrSges(5,1) * t148 + mrSges(5,2) * t147;
t522 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t257 + t69;
t489 = m(6) * pkin(4);
t521 = -t489 - mrSges(5,1);
t520 = qJ(6) * t351;
t207 = t333 * t305;
t444 = mrSges(7,3) * t351;
t131 = -mrSges(7,2) * t276 + t444;
t446 = mrSges(6,3) * t351;
t132 = -mrSges(6,2) * t276 + t446;
t519 = t131 + t132;
t443 = mrSges(7,3) * t160;
t133 = mrSges(7,1) * t276 - t443;
t445 = mrSges(6,3) * t160;
t134 = mrSges(6,1) * t276 - t445;
t518 = t133 + t134;
t517 = pkin(5) * t551 + t513;
t263 = -pkin(3) * t309 - pkin(9) * t305 - pkin(2);
t242 = t308 * t263;
t415 = t305 * t308;
t163 = -pkin(10) * t415 + t242 + (-pkin(8) * t304 - pkin(4)) * t309;
t288 = pkin(8) * t413;
t202 = t304 * t263 + t288;
t418 = t304 * t305;
t178 = -pkin(10) * t418 + t202;
t83 = t303 * t163 + t307 * t178;
t412 = t309 * t310;
t191 = (-t304 * t412 + t306 * t308) * t405;
t255 = t349 * qJD(3);
t399 = qJD(3) * t305;
t383 = pkin(8) * t399;
t407 = t308 * t255 + t304 * t383;
t516 = -qJD(4) * t202 - t191 + t407;
t125 = t304 * t255 + t263 * t394 + (-t305 * t398 - t309 * t396) * pkin(8);
t417 = t304 * t306;
t192 = (t308 * t412 + t417) * t405;
t515 = -t192 + t125;
t381 = mrSges(4,3) * t402;
t514 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t243 + mrSges(5,2) * t244 + t381;
t512 = t234 * t532 - t49 * t535 + t50 * t533;
t420 = t301 * t306;
t223 = t302 * t305 + t309 * t420;
t419 = t301 * t310;
t161 = -t223 * t296 - t297 * t419;
t511 = t538 * (-t223 * t297 + t296 * t419) + t539 * t161;
t426 = sin(pkin(11));
t354 = t426 * t306;
t427 = cos(pkin(11));
t355 = t427 * t310;
t219 = -t302 * t354 + t355;
t357 = t301 * t426;
t173 = t219 * t309 + t305 * t357;
t353 = t426 * t310;
t356 = t427 * t306;
t218 = t302 * t353 + t356;
t109 = -t173 * t296 + t218 * t297;
t510 = t538 * (-t173 * t297 - t218 * t296) + t539 * t109;
t238 = Ifges(5,4) * t243;
t142 = t244 * Ifges(5,1) + t286 * Ifges(5,5) + t238;
t294 = Ifges(4,4) * t400;
t509 = Ifges(4,1) * t402 + Ifges(4,5) * qJD(3) + t308 * t142 + t294;
t389 = qJDD(1) * t302;
t105 = -t259 * t399 + t305 * t389 + t309 * t546;
t106 = -t259 * t397 - t305 * t546 + t309 * t389;
t508 = t105 * t309 - t106 * t305;
t217 = t302 * t356 + t353;
t358 = t301 * t427;
t171 = t217 * t309 - t305 * t358;
t216 = -t302 * t355 + t354;
t107 = -t171 * t296 + t216 * t297;
t507 = t538 * (-t171 * t297 - t216 * t296) + t539 * t107;
t272 = t306 * t367;
t210 = t310 * t390 - t272;
t199 = -qJDD(2) * pkin(2) - t210;
t124 = -pkin(3) * t256 - pkin(9) * t257 + t199;
t183 = qJD(3) * pkin(9) + t196;
t376 = t310 * t405;
t198 = qJD(2) * t263 - t376;
t89 = qJDD(3) * pkin(9) + t105;
t30 = t304 * t124 - t183 * t396 + t198 * t394 + t308 * t89;
t113 = t183 * t308 + t198 * t304;
t31 = -qJD(4) * t113 + t308 * t124 - t304 * t89;
t506 = t30 * t308 - t304 * t31;
t505 = -m(7) - m(4) - m(6);
t503 = t244 * Ifges(5,5) + t243 * Ifges(5,6) + t286 * Ifges(5,3) - t160 * t535 + t276 * t532 + t351 * t533;
t502 = -m(5) + t505;
t182 = -qJD(3) * pkin(3) - t195;
t137 = -pkin(4) * t243 + t182;
t68 = -pkin(5) * t351 + qJD(6) + t137;
t79 = -mrSges(7,1) * t351 + mrSges(7,2) * t160;
t501 = -m(7) * t68 - t79;
t500 = t538 * t296 - t297 * t539 + mrSges(4,1) - t554;
t499 = mrSges(4,2) + t553;
t497 = -m(5) * t182 - t514;
t112 = -t183 * t304 + t308 * t198;
t77 = -pkin(10) * t244 + t112;
t67 = pkin(4) * t286 + t77;
t78 = pkin(10) * t243 + t113;
t74 = t307 * t78;
t33 = t303 * t67 + t74;
t27 = t33 + t520;
t496 = -t68 * mrSges(7,1) + t33 * mrSges(6,3) + t27 * mrSges(7,3);
t72 = t303 * t78;
t32 = t307 * t67 - t72;
t26 = t32 - t548;
t25 = pkin(5) * t276 + t26;
t495 = t68 * mrSges(7,2) - mrSges(6,3) * t32 - mrSges(7,3) * t25;
t494 = -t31 * mrSges(5,1) + t30 * mrSges(5,2);
t492 = -m(5) * pkin(8) - m(6) * t455 + t545;
t17 = pkin(4) * t240 - pkin(10) * t147 + t31;
t22 = pkin(10) * t148 + t30;
t6 = -qJD(5) * t33 + t307 * t17 - t22 * t303;
t2 = pkin(5) * t234 - qJ(6) * t49 - qJD(6) * t160 + t6;
t5 = t303 * t17 + t307 * t22 + t67 * t392 - t393 * t78;
t3 = qJ(6) * t50 + qJD(6) * t351 + t5;
t491 = -t6 * mrSges(6,1) - t2 * mrSges(7,1) + t5 * mrSges(6,2) + t3 * mrSges(7,2);
t442 = Ifges(4,4) * t305;
t343 = t309 * Ifges(4,2) + t442;
t490 = t27 * mrSges(7,2) + t33 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(2) * t343 / 0.2e1 - t25 * mrSges(7,1) - t32 * mrSges(6,1);
t312 = qJD(2) ^ 2;
t488 = m(7) * pkin(5);
t487 = t49 / 0.2e1;
t486 = t50 / 0.2e1;
t485 = Ifges(5,1) * t482 + Ifges(5,4) * t481 + Ifges(5,5) * t467;
t478 = t351 / 0.2e1;
t476 = -t160 / 0.2e1;
t475 = t160 / 0.2e1;
t468 = t234 / 0.2e1;
t465 = t244 / 0.2e1;
t462 = -t276 / 0.2e1;
t461 = t276 / 0.2e1;
t457 = pkin(4) * t244;
t454 = pkin(4) * t307;
t298 = t305 * pkin(8);
t35 = t307 * t77 - t72;
t441 = Ifges(4,4) * t309;
t440 = Ifges(5,4) * t304;
t439 = Ifges(5,4) * t308;
t436 = t112 * mrSges(5,3);
t435 = t113 * mrSges(5,3);
t434 = t244 * Ifges(5,4);
t90 = -qJDD(3) * pkin(3) - t106;
t429 = t305 * t90;
t422 = t296 * t309;
t421 = t297 * t309;
t416 = t304 * t309;
t258 = pkin(4) * t418 + t298;
t401 = qJD(2) * t306;
t395 = qJD(4) * t305;
t385 = pkin(4) * t393;
t384 = pkin(4) * t392;
t295 = pkin(8) * t397;
t380 = mrSges(4,3) * t400;
t379 = Ifges(5,5) * t147 + Ifges(5,6) * t148 + Ifges(5,3) * t240;
t197 = pkin(4) * t320 + t295;
t375 = t301 * t401;
t374 = t310 * t403;
t141 = t243 * Ifges(5,2) + t286 * Ifges(5,6) + t434;
t370 = -t304 * t141 / 0.2e1;
t13 = -t50 * mrSges(7,1) + t49 * mrSges(7,2);
t34 = -t303 * t77 - t74;
t352 = t391 / 0.2e1;
t82 = t307 * t163 - t178 * t303;
t185 = t307 * t270 + t271 * t303;
t345 = Ifges(5,1) * t308 - t440;
t344 = Ifges(5,1) * t304 + t439;
t342 = -Ifges(5,2) * t304 + t439;
t341 = Ifges(5,2) * t308 + t440;
t340 = Ifges(4,5) * t309 - Ifges(4,6) * t305;
t339 = Ifges(5,5) * t308 - Ifges(5,6) * t304;
t338 = Ifges(5,5) * t304 + Ifges(5,6) * t308;
t176 = -t223 * t304 - t308 * t419;
t331 = -t223 * t308 + t304 * t419;
t85 = t176 * t307 + t303 * t331;
t86 = t176 * t303 - t307 * t331;
t222 = -t302 * t309 + t305 * t420;
t330 = t182 * t346;
t260 = -qJD(2) * pkin(2) - t376;
t329 = t260 * (mrSges(4,1) * t305 + mrSges(4,2) * t309);
t328 = t305 * (Ifges(4,1) * t309 - t442);
t81 = t332 * qJD(3) + (-t288 + (pkin(10) * t305 - t263) * t304) * qJD(4) + t407;
t88 = -pkin(10) * t320 + t125;
t23 = t163 * t392 - t178 * t393 + t303 * t81 + t307 * t88;
t321 = -t304 * t395 + t308 * t397;
t317 = Ifges(5,5) * t305 + t309 * t345;
t316 = Ifges(5,6) * t305 + t309 * t342;
t315 = Ifges(5,3) * t305 + t309 * t339;
t24 = -qJD(5) * t83 - t303 * t88 + t307 * t81;
t314 = -t491 + t512;
t53 = -pkin(4) * t148 + t90;
t292 = pkin(5) + t454;
t265 = -qJD(3) * mrSges(4,2) + t380;
t250 = t348 * qJD(2);
t212 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t256;
t209 = t218 * pkin(2);
t208 = t216 * pkin(2);
t206 = t246 * t305;
t205 = pkin(5) * t333 - t293;
t201 = -pkin(8) * t416 + t242;
t194 = mrSges(5,1) * t286 - mrSges(5,3) * t244;
t193 = -mrSges(5,2) * t286 + mrSges(5,3) * t243;
t179 = -mrSges(4,1) * t256 + mrSges(4,2) * t257;
t175 = -qJD(3) * t222 + t309 * t374;
t174 = qJD(3) * t223 + t305 * t374;
t172 = t219 * t305 - t309 * t357;
t170 = t217 * t305 + t309 * t358;
t165 = pkin(5) * t206 + t258;
t139 = -qJ(6) * t333 + t186;
t138 = -qJ(6) * t246 + t185;
t123 = pkin(5) * t160 + t457;
t119 = -mrSges(5,2) * t240 + mrSges(5,3) * t148;
t118 = mrSges(5,1) * t240 - mrSges(5,3) * t147;
t115 = t191 * t303 + t192 * t307;
t114 = t191 * t307 - t192 * t303;
t92 = -qJD(3) * t327 + t207 * t504;
t91 = -qJD(3) * t326 - t168 * t305;
t80 = -mrSges(6,1) * t351 + mrSges(6,2) * t160;
t71 = qJD(4) * t176 + t175 * t308 + t304 * t375;
t70 = qJD(4) * t331 - t175 * t304 + t308 * t375;
t58 = -pkin(5) * t92 + t197;
t57 = -qJ(6) * t206 + t83;
t55 = t147 * Ifges(5,4) + t148 * Ifges(5,2) + t240 * Ifges(5,6);
t54 = -pkin(5) * t309 + qJ(6) * t207 + t82;
t39 = mrSges(6,1) * t234 - mrSges(6,3) * t49;
t38 = mrSges(7,1) * t234 - mrSges(7,3) * t49;
t29 = t35 - t548;
t28 = t34 - t520;
t20 = -pkin(5) * t50 + qJDD(6) + t53;
t19 = -qJD(5) * t86 - t303 * t71 + t307 * t70;
t18 = qJD(5) * t85 + t303 * t70 + t307 * t71;
t14 = -mrSges(6,1) * t50 + mrSges(6,2) * t49;
t8 = qJ(6) * t92 - qJD(6) * t206 + t23;
t7 = pkin(5) * t399 - qJ(6) * t91 + qJD(6) * t207 + t24;
t1 = [m(2) * qJDD(1) + t176 * t118 - t331 * t119 + t175 * t265 + t71 * t193 + t70 * t194 + t223 * t212 + t530 * t86 + (t38 + t39) * t85 + t518 * t19 + t519 * t18 + (t13 + t14 + t522) * t222 + (t79 + t80 + t514) * t174 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t312 - t179) * t310 + (-mrSges(3,1) * t312 - mrSges(3,2) * qJDD(2) - qJD(2) * t250) * t306) * t301 + (-m(2) - m(3) + t502) * g(3) + m(4) * (t105 * t223 - t106 * t222 - t174 * t195 + t175 * t196 + (-t199 * t310 + t260 * t401) * t301) + m(3) * (qJDD(1) * t302 ^ 2 + (t210 * t310 + t211 * t306) * t301) + m(7) * (t174 * t68 + t18 * t27 + t19 * t25 + t2 * t85 + t20 * t222 + t3 * t86) + m(6) * (t137 * t174 + t18 * t33 + t19 * t32 + t222 * t53 + t5 * t86 + t6 * t85) + m(5) * (t112 * t70 + t113 * t71 + t174 * t182 + t176 * t31 + t222 * t90 - t30 * t331); (-t195 * t397 + t508) * mrSges(4,3) + (t502 * (pkin(2) * t419 + pkin(8) * t420) + (t539 * (t296 * t306 + t297 * t412) + t538 * (-t296 * t412 + t297 * t306) - t417 * t489 + t545 * t306 + t543 * t310) * t301) * g(3) + (t137 * t197 + t258 * t53 + t5 * t83 + t6 * t82 + (-t115 + t23) * t33 + (-t114 + t24) * t32) * m(6) + (-pkin(2) * t199 + ((-t195 * t309 - t196 * t305) * qJD(3) + t508) * pkin(8) - (t260 * t306 + (-t195 * t305 + t196 * t309) * t310) * t405) * m(4) + t509 * t397 / 0.2e1 - (t379 + t512) * t309 / 0.2e1 + t514 * t295 + t370 * t397 + t515 * t193 + (t201 * t31 + t202 * t30 + (t182 * t397 + t429) * pkin(8) + t515 * t113 + t516 * t112) * m(5) + t516 * t194 - t518 * t114 - t519 * t115 + (-Ifges(5,6) * t481 - Ifges(5,5) * t482 + Ifges(4,4) * t541 + Ifges(4,2) * t542 - Ifges(5,3) * t467 + (-Ifges(4,2) * t305 + t441) * t352 + pkin(8) * t212 - t265 * t376 + t535 * t487 - t533 * t486 - t532 * t468 + t491 + t494) * t309 + (Ifges(4,1) * t257 + Ifges(4,4) * t542 + t339 * t467 + t342 * t481 + t345 * t482) * t305 - (t308 * t141 + t304 * t142) * t395 / 0.2e1 + t415 * t485 + t250 * t377 - t55 * t418 / 0.2e1 - t265 * t383 + t182 * (mrSges(5,1) * t320 + mrSges(5,2) * t321) + (m(5) * t209 + t539 * (-t218 * t421 + t219 * t296) + t538 * (t218 * t422 + t219 * t297) + t505 * (pkin(8) * t219 - t209) + t492 * t219 - t543 * t218) * g(1) + (m(5) * t208 + t539 * (-t216 * t421 + t217 * t296) + t538 * (t216 * t422 + t217 * t297) + t505 * (pkin(8) * t217 - t208) + t492 * t217 - t543 * t216) * g(2) + t346 * t429 + (-t206 * t5 + t207 * t6 - t32 * t91 + t33 * t92) * mrSges(6,3) + (t2 * t207 - t206 * t3 - t25 * t91 + t27 * t92) * mrSges(7,3) + t53 * (mrSges(6,1) * t206 - mrSges(6,2) * t207) + t20 * (mrSges(7,1) * t206 - mrSges(7,2) * t207) + t243 * (qJD(3) * t316 - t341 * t395) / 0.2e1 + t286 * (qJD(3) * t315 - t338 * t395) / 0.2e1 + (-m(6) * t137 + t497 + t501 - t80) * t305 * t376 + (t273 - t211) * mrSges(3,2) + t328 * t352 - t199 * t348 + qJD(3) ^ 2 * t340 / 0.2e1 + t522 * t298 + t526 * t91 / 0.2e1 - t531 * t207 / 0.2e1 + t206 * t549 + t92 * t550 + (t272 + t210) * mrSges(3,1) + qJD(3) * t329 + (qJD(3) * t317 - t344 * t395) * t465 + (-t112 * t321 - t113 * t320 - t30 * t418 - t31 * t415) * mrSges(5,3) + (t533 * t92 - t535 * t91) * t461 + (-t206 * t533 + t207 * t535) * t468 + (t534 * t92 + t536 * t91) * t478 + (-t206 * t534 - t207 * t536) * t486 + (t536 * t92 + t537 * t91) * t475 + (-t206 * t536 - t207 * t537) * t487 + Ifges(3,3) * qJDD(2) + qJDD(3) * (Ifges(4,5) * t305 + Ifges(4,6) * t309) + t258 * t14 + t201 * t118 + t202 * t119 + t197 * t80 - pkin(2) * t179 + t165 * t13 + t8 * t131 + t23 * t132 + t7 * t133 + t24 * t134 + t137 * (-mrSges(6,1) * t92 + mrSges(6,2) * t91) + t68 * (-mrSges(7,1) * t92 + mrSges(7,2) * t91) + t82 * t39 + t83 * t41 + t58 * t79 + t57 * t40 + t54 * t38 + t441 * t541 + t343 * t542 + (-t196 * mrSges(4,3) + t503 / 0.2e1 + t112 * mrSges(5,1) - t113 * mrSges(5,2) + t461 * t532 - t475 * t535 + t478 * t533 - t490) * t399 + (t165 * t20 + t2 * t54 + t3 * t57 + t58 * t68 + (-t115 + t8) * t27 + (-t114 + t7) * t25) * m(7); (m(5) * ((-t112 * t308 - t113 * t304) * qJD(4) + t506) - t194 * t394 - t193 * t396 + t308 * t119 - t304 * t118) * pkin(9) + t506 * mrSges(5,3) - (-Ifges(4,2) * t402 + t294 + t509) * t400 / 0.2e1 + t513 * t80 + t517 * t79 + (t170 * t500 + t171 * t499) * g(2) + (t172 * t500 + t173 * t499) * g(1) + (t222 * t500 + t223 * t499) * g(3) - t340 * t391 / 0.2e1 - t330 * t400 + (t243 * t342 + t244 * t345 + t286 * t339) * qJD(4) / 0.2e1 - (t243 * t316 + t244 * t317 + t286 * t315) * qJD(2) / 0.2e1 + t341 * t481 + t344 * t482 + t304 * t485 + (t142 / 0.2e1 - t436) * t394 + (mrSges(6,1) * t551 - mrSges(6,2) * t552) * t137 + (-t329 - t112 * (mrSges(5,1) * t305 - mrSges(5,3) * t413) - t113 * (-mrSges(5,2) * t305 - mrSges(5,3) * t416)) * qJD(2) + (-pkin(3) * t90 - t112 * t135 - t113 * t136) * m(5) - t312 * t328 / 0.2e1 + t90 * t347 + t141 * t373 / 0.2e1 + t524 * t134 + t525 * t132 + (t137 * t513 + t185 * t6 + t186 * t5 - t293 * t53 + t32 * t524 + t33 * t525) * m(6) + t528 * t133 + (t138 * t2 + t139 * t3 + t20 * t205 + t25 * t528 + t27 * t529 + t517 * t68) * m(7) + t529 * t131 + t531 * t246 / 0.2e1 + (t370 + t330) * qJD(4) + t333 * t549 + (t532 * t462 - t535 * t476 + t533 * t479 + t490) * t402 + t338 * t467 + (t380 - t265) * t195 + (t381 + t497) * t196 - t503 * t402 / 0.2e1 - (t533 * t461 + t536 * t475 + t534 * t478 + t496) * t168 - (-t535 * t461 + t537 * t475 + t536 * t478 + t495) * t167 - t396 * t435 + Ifges(4,3) * qJDD(3) + t308 * t55 / 0.2e1 + (t203 / 0.2e1 - t168 / 0.2e1) * t527 - t293 * t14 + (-t203 * t534 - t204 * t536) * t479 + (-t203 * t536 - t204 * t537) * t476 + (-t2 * t246 + t203 * t27 - t204 * t25 - t3 * t333) * mrSges(7,3) + (t203 * t33 - t204 * t32 - t246 * t6 - t333 * t5) * mrSges(6,3) - t68 * (mrSges(7,1) * t203 - mrSges(7,2) * t204) + (t204 / 0.2e1 - t167 / 0.2e1) * t526 + (-t203 * t533 + t204 * t535) * t462 + (-t246 * t535 - t333 * t533) * t468 + (t246 * t536 - t333 * t534) * t486 + (t246 * t537 - t333 * t536) * t487 + t53 * (mrSges(6,1) * t333 + mrSges(6,2) * t246) + t20 * (mrSges(7,1) * t333 + mrSges(7,2) * t246) + Ifges(4,6) * t256 + Ifges(4,5) * t257 + t205 * t13 - t135 * t194 + t185 * t39 + t186 * t41 - t136 * t193 + t138 * t38 + t139 * t40 - t105 * mrSges(4,2) + t106 * mrSges(4,1) - pkin(3) * t69; -t494 + t530 * pkin(4) * t303 - (-Ifges(5,2) * t244 + t142 + t238) * t243 / 0.2e1 + (-t123 * t68 - t25 * t28 - t27 * t29 + t2 * t292 + (t3 * t303 + (-t25 * t303 + t27 * t307) * qJD(5)) * pkin(4)) * m(7) + t314 + (t303 * t5 + t307 * t6 + (-t303 * t32 + t307 * t33) * qJD(5)) * t489 - t80 * t457 - m(6) * (t137 * t457 + t32 * t34 + t33 * t35) + t379 + t244 * t435 + t243 * t436 + t39 * t454 + (-(-t171 * t308 - t216 * t304) * mrSges(5,2) - m(7) * (-t171 * t261 + t216 * t262) + t521 * (-t171 * t304 + t216 * t308) + t507) * g(2) + (-(-t173 * t308 - t218 * t304) * mrSges(5,2) - m(7) * (-t173 * t261 + t218 * t262) + t521 * (-t173 * t304 + t218 * t308) + t510) * g(1) - t244 * (Ifges(5,1) * t243 - t434) / 0.2e1 + t141 * t465 + (-t137 * mrSges(6,2) - t535 * t462 + t537 * t476 + t536 * t479 - t495) * t351 + t292 * t38 - t286 * (Ifges(5,5) * t243 - Ifges(5,6) * t244) / 0.2e1 + (-m(7) * (-t223 * t261 - t262 * t419) - mrSges(5,2) * t331 + t521 * t176 + t511) * g(3) - t182 * (mrSges(5,1) * t244 + mrSges(5,2) * t243) + t113 * t194 - t112 * t193 - t123 * t79 + (-t385 - t34) * t134 + t526 * t479 + (t384 - t29) * t131 + (-t137 * mrSges(6,1) - t533 * t462 - t536 * t476 - t534 * t479 + t496 + t550) * t160 + (t384 - t35) * t132 + (-t385 - t28) * t133; t314 + t2 * t488 + t25 * t444 - t68 * (mrSges(7,1) * t160 + mrSges(7,2) * t351) - t137 * (mrSges(6,1) * t160 + mrSges(6,2) * t351) - t26 * t131 + (t351 * t537 - t544) * t476 + t527 * t475 + (-t160 * t533 - t351 * t535) * t462 + (t445 + t134) * t33 + (t446 - t132) * t32 + (-m(7) * (-t25 + t26) + t443 + t133) * t27 + (-t161 * t488 + t511) * g(3) + (-t107 * t488 + t507) * g(2) + (-t109 * t488 + t510) * g(1) + (-t160 * t534 + t526 + t547) * t479 + (t160 * t501 + t38) * pkin(5); -t351 * t131 + t160 * t133 + (-g(1) * t172 - g(2) * t170 - g(3) * t222 + t25 * t160 - t27 * t351 + t20) * m(7) + t13;];
tau  = t1;
