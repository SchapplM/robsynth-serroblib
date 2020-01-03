% Calculate vector of inverse dynamics joint torques for
% S5RRRPR12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:05
% EndTime: 2019-12-31 21:38:05
% DurationCPUTime: 33.30s
% Computational Cost: add. (11961->849), mult. (29592->1208), div. (0->0), fcn. (23415->14), ass. (0->354)
t306 = sin(qJ(2));
t408 = cos(pkin(5));
t358 = t306 * t408;
t293 = pkin(1) * t358;
t305 = sin(qJ(3));
t309 = cos(qJ(3));
t332 = pkin(3) * t305 - qJ(4) * t309;
t301 = sin(pkin(5));
t310 = cos(qJ(2));
t398 = t301 * t310;
t528 = -(t293 + (pkin(7) + t332) * t398) * qJD(1) + qJD(3) * t332 - qJD(4) * t305;
t356 = t408 * qJD(1);
t348 = pkin(1) * t356;
t391 = qJD(1) * t301;
t372 = t306 * t391;
t233 = -pkin(7) * t372 + t310 * t348;
t326 = (pkin(2) * t306 - pkin(8) * t310) * t301;
t234 = qJD(1) * t326;
t162 = t309 * t233 + t305 * t234;
t140 = qJ(4) * t372 + t162;
t300 = sin(pkin(10));
t302 = cos(pkin(10));
t387 = qJD(3) * t305;
t383 = pkin(8) * t387;
t493 = t528 * t302 + (t140 + t383) * t300;
t527 = -t302 * t140 + t528 * t300;
t295 = pkin(4) * t302 + pkin(3);
t341 = -mrSges(5,1) * t302 + mrSges(5,2) * t300;
t526 = -m(5) * pkin(3) - m(6) * t295 + t341;
t342 = mrSges(4,1) * t309 - mrSges(4,2) * t305;
t423 = pkin(9) + qJ(4);
t482 = -m(5) * qJ(4) - m(6) * t423 - mrSges(5,3) - mrSges(6,3);
t299 = pkin(10) + qJ(5);
t296 = sin(t299);
t297 = cos(t299);
t513 = -mrSges(6,1) * t297 + mrSges(6,2) * t296;
t484 = -t513 - t526;
t525 = t305 * t482 - t309 * t484 - t342;
t507 = -m(5) - m(4);
t340 = t300 * mrSges(5,1) + t302 * mrSges(5,2);
t524 = mrSges(4,3) + t340;
t395 = t309 * t310;
t200 = (t300 * t306 + t302 * t395) * t391;
t371 = t310 * t391;
t349 = t305 * t371;
t396 = t302 * t309;
t523 = -pkin(4) * t349 + t200 * pkin(9) + (pkin(4) * t305 - pkin(9) * t396) * qJD(3) + t493;
t199 = (-t300 * t395 + t302 * t306) * t391;
t397 = t302 * t305;
t402 = t300 * t309;
t522 = -pkin(9) * t199 + (-pkin(8) * t397 - pkin(9) * t402) * qJD(3) + t527;
t304 = sin(qJ(5));
t308 = cos(qJ(5));
t329 = t356 + qJD(2);
t213 = t305 * t329 + t309 * t372;
t272 = qJD(3) - t371;
t168 = t213 * t302 + t272 * t300;
t212 = t305 * t372 - t309 * t329;
t197 = -pkin(2) * t329 - t233;
t105 = t212 * pkin(3) - t213 * qJ(4) + t197;
t393 = pkin(7) * t398 + t293;
t224 = t408 * pkin(8) + t393;
t198 = qJD(2) * pkin(8) + qJD(1) * t224;
t204 = (-pkin(2) * t310 - pkin(8) * t306 - pkin(1)) * t391;
t126 = t198 * t309 + t204 * t305;
t110 = qJ(4) * t272 + t126;
t57 = t302 * t105 - t110 * t300;
t34 = pkin(4) * t212 - pkin(9) * t168 + t57;
t355 = -t213 * t300 + t302 * t272;
t58 = t300 * t105 + t302 * t110;
t42 = pkin(9) * t355 + t58;
t11 = -t304 * t42 + t308 * t34;
t12 = t304 * t34 + t308 * t42;
t415 = t126 * mrSges(4,3);
t496 = t272 * Ifges(4,6);
t519 = -t415 - t496 / 0.2e1 + t197 * mrSges(4,1) + t57 * mrSges(5,1) + t11 * mrSges(6,1) - t58 * mrSges(5,2) - t12 * mrSges(6,2);
t518 = Ifges(5,6) * t355;
t517 = t168 * Ifges(5,5);
t384 = qJDD(1) * t301;
t516 = pkin(7) * t384 + qJD(2) * t348;
t354 = t408 * qJDD(1);
t385 = qJD(1) * qJD(2);
t515 = -pkin(7) * t301 * t385 + pkin(1) * t354;
t477 = m(6) - t507;
t514 = pkin(2) * t477 + mrSges(3,1) - t525;
t426 = pkin(4) * t300;
t375 = pkin(8) + t426;
t469 = m(6) * t375 - mrSges(3,2) + t524;
t512 = -t168 * t304 + t308 * t355;
t92 = t168 * t308 + t304 * t355;
t239 = (-qJDD(1) * t310 + t306 * t385) * t301;
t228 = qJDD(3) + t239;
t171 = t306 * t515 + t310 * t516;
t284 = t354 + qJDD(2);
t147 = pkin(8) * t284 + t171;
t240 = (qJDD(1) * t306 + t310 * t385) * t301;
t381 = pkin(1) * t384;
t155 = pkin(2) * t239 - pkin(8) * t240 - t381;
t386 = qJD(3) * t309;
t55 = t309 * t147 + t305 * t155 - t198 * t387 + t204 * t386;
t41 = qJ(4) * t228 + qJD(4) * t272 + t55;
t388 = qJD(3) * t212;
t132 = t309 * t240 + t305 * t284 - t388;
t133 = qJD(3) * t213 + t305 * t240 - t309 * t284;
t172 = -t306 * t516 + t310 * t515;
t148 = -t284 * pkin(2) - t172;
t45 = t133 * pkin(3) - t132 * qJ(4) - t213 * qJD(4) + t148;
t14 = t300 * t45 + t302 * t41;
t97 = -t132 * t300 + t228 * t302;
t10 = pkin(9) * t97 + t14;
t13 = -t300 * t41 + t302 * t45;
t98 = t132 * t302 + t228 * t300;
t8 = pkin(4) * t133 - pkin(9) * t98 + t13;
t1 = qJD(5) * t11 + t10 * t308 + t304 * t8;
t2 = -qJD(5) * t12 - t10 * t304 + t308 * t8;
t511 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t442 = t133 / 0.2e1;
t447 = t98 / 0.2e1;
t510 = Ifges(5,1) * t447 + Ifges(5,5) * t442;
t412 = t213 * Ifges(4,4);
t118 = -t212 * Ifges(4,2) + t412 + t496;
t206 = qJD(5) + t212;
t437 = t206 / 0.2e1;
t449 = t92 / 0.2e1;
t451 = t512 / 0.2e1;
t500 = t92 * Ifges(6,5) + Ifges(6,6) * t512 + t212 * Ifges(5,3) + t206 * Ifges(6,3) + t517 + t518;
t509 = -t118 / 0.2e1 + Ifges(6,5) * t449 + Ifges(6,6) * t451 + Ifges(6,3) * t437 + t500 / 0.2e1;
t28 = qJD(5) * t512 + t304 * t97 + t308 * t98;
t463 = t28 / 0.2e1;
t29 = -qJD(5) * t92 - t304 * t98 + t308 * t97;
t462 = t29 / 0.2e1;
t448 = t97 / 0.2e1;
t124 = qJDD(5) + t133;
t445 = t124 / 0.2e1;
t444 = t132 / 0.2e1;
t443 = -t133 / 0.2e1;
t432 = t228 / 0.2e1;
t5 = Ifges(6,5) * t28 + Ifges(6,6) * t29 + Ifges(6,3) * t124;
t506 = t98 * Ifges(5,5) + t97 * Ifges(5,6) + t133 * Ifges(5,3) + t5;
t266 = -pkin(3) * t309 - qJ(4) * t305 - pkin(2);
t259 = t302 * t266;
t175 = -pkin(9) * t397 + t259 + (-pkin(8) * t300 - pkin(4)) * t309;
t211 = pkin(8) * t396 + t300 * t266;
t403 = t300 * t305;
t184 = -pkin(9) * t403 + t211;
t109 = t175 * t304 + t184 * t308;
t502 = -qJD(5) * t109 - t304 * t522 + t308 * t523;
t108 = t175 * t308 - t184 * t304;
t501 = qJD(5) * t108 + t304 * t523 + t308 * t522;
t497 = t272 * Ifges(4,5);
t267 = t423 * t300;
t268 = t423 * t302;
t193 = -t267 * t308 - t268 * t304;
t331 = t300 * t304 - t302 * t308;
t406 = t212 * t302;
t125 = -t305 * t198 + t204 * t309;
t149 = pkin(3) * t213 + qJ(4) * t212;
t73 = -t125 * t300 + t302 * t149;
t52 = pkin(4) * t213 + pkin(9) * t406 + t73;
t407 = t212 * t300;
t74 = t302 * t125 + t300 * t149;
t63 = pkin(9) * t407 + t74;
t495 = -qJD(4) * t331 + qJD(5) * t193 - t304 * t52 - t308 * t63;
t194 = -t267 * t304 + t268 * t308;
t262 = t300 * t308 + t302 * t304;
t494 = -qJD(4) * t262 - qJD(5) * t194 + t304 * t63 - t308 * t52;
t492 = -t302 * t383 + t527;
t128 = t199 * t308 - t200 * t304;
t246 = t331 * qJD(5);
t164 = t246 * t305 - t262 * t386;
t489 = t128 - t164;
t129 = t199 * t304 + t200 * t308;
t247 = t262 * qJD(5);
t163 = -t247 * t305 - t331 * t386;
t488 = t129 - t163;
t130 = t262 * t212;
t487 = t130 + t247;
t131 = t331 * t212;
t486 = t131 + t246;
t161 = -t305 * t233 + t234 * t309;
t141 = -pkin(3) * t372 - t161;
t485 = pkin(4) * t199 + t375 * t386 - t141;
t357 = t310 * t408;
t401 = t301 * t306;
t255 = pkin(1) * t357 - pkin(7) * t401;
t483 = t349 - t387;
t107 = -pkin(3) * t272 + qJD(4) - t125;
t56 = -t305 * t147 + t155 * t309 - t198 * t386 - t204 * t387;
t46 = -pkin(3) * t228 + qJDD(4) - t56;
t481 = t107 * t386 + t305 * t46;
t480 = -t305 * t56 + t309 * t55;
t479 = -t13 * t300 + t14 * t302;
t79 = t168 * Ifges(5,4) + Ifges(5,2) * t355 + Ifges(5,6) * t212;
t454 = -t79 / 0.2e1;
t478 = -t125 * mrSges(4,3) + t300 * t454;
t422 = Ifges(3,4) * t306;
t467 = t301 ^ 2;
t476 = (pkin(1) * (mrSges(3,1) * t306 + mrSges(3,2) * t310) - t306 * (Ifges(3,1) * t310 - t422) / 0.2e1) * t467;
t475 = mrSges(4,1) + t484;
t470 = mrSges(4,2) + t482;
t352 = mrSges(3,3) * t372;
t472 = -m(4) * t197 + mrSges(3,1) * t329 - mrSges(4,1) * t212 - mrSges(4,2) * t213 - t352;
t339 = t296 * mrSges(6,1) + t297 * mrSges(6,2);
t471 = pkin(8) * t507 - t339 - t469;
t468 = -mrSges(4,1) + t526;
t465 = Ifges(6,4) * t463 + Ifges(6,2) * t462 + Ifges(6,6) * t445;
t464 = Ifges(6,1) * t463 + Ifges(6,4) * t462 + Ifges(6,5) * t445;
t32 = t98 * Ifges(5,4) + t97 * Ifges(5,2) + t133 * Ifges(5,6);
t461 = t32 / 0.2e1;
t460 = Ifges(5,4) * t448 + t510;
t428 = Ifges(6,4) * t92;
t38 = Ifges(6,2) * t512 + Ifges(6,6) * t206 + t428;
t459 = -t38 / 0.2e1;
t458 = t38 / 0.2e1;
t89 = Ifges(6,4) * t512;
t39 = Ifges(6,1) * t92 + Ifges(6,5) * t206 + t89;
t457 = -t39 / 0.2e1;
t456 = t39 / 0.2e1;
t455 = Ifges(4,1) * t444 + Ifges(4,4) * t443 + Ifges(4,5) * t432;
t80 = t168 * Ifges(5,1) + Ifges(5,4) * t355 + Ifges(5,5) * t212;
t453 = t80 / 0.2e1;
t452 = -t512 / 0.2e1;
t450 = -t92 / 0.2e1;
t441 = -t355 / 0.2e1;
t440 = -t168 / 0.2e1;
t438 = -t206 / 0.2e1;
t436 = -t212 / 0.2e1;
t435 = t212 / 0.2e1;
t433 = t213 / 0.2e1;
t429 = cos(qJ(1));
t427 = pkin(1) * t301;
t390 = qJD(2) * t306;
t394 = pkin(2) * t398 + pkin(8) * t401;
t225 = -t394 - t427;
t235 = qJD(2) * t326;
t237 = t255 * qJD(2);
t87 = -t224 * t387 + t225 * t386 + t305 * t235 + t309 * t237;
t77 = (qJ(4) * t390 - qJD(4) * t310) * t301 + t87;
t399 = t301 * t309;
t249 = t305 * t408 + t306 * t399;
t389 = qJD(2) * t310;
t369 = t301 * t389;
t182 = qJD(3) * t249 + t305 * t369;
t248 = t305 * t401 - t309 * t408;
t183 = -qJD(3) * t248 + t309 * t369;
t238 = t393 * qJD(2);
t86 = t182 * pkin(3) - t183 * qJ(4) - t249 * qJD(4) + t238;
t36 = t300 * t86 + t302 * t77;
t421 = Ifges(3,4) * t310;
t420 = Ifges(4,4) * t305;
t419 = Ifges(4,4) * t309;
t418 = Ifges(5,4) * t300;
t417 = Ifges(5,4) * t302;
t307 = sin(qJ(1));
t344 = t408 * t429;
t250 = t306 * t307 - t310 * t344;
t405 = t250 * t296;
t404 = t250 * t297;
t400 = t301 * t307;
t223 = -pkin(2) * t408 - t255;
t137 = t248 * pkin(3) - t249 * qJ(4) + t223;
t152 = t309 * t224 + t305 * t225;
t138 = -qJ(4) * t398 + t152;
t69 = t300 * t137 + t302 * t138;
t392 = t429 * pkin(1) + pkin(7) * t400;
t382 = pkin(8) * t386;
t378 = Ifges(4,5) * t132 - Ifges(4,6) * t133 + Ifges(4,3) * t228;
t377 = Ifges(3,5) * t240 - Ifges(3,6) * t239 + Ifges(3,3) * t284;
t253 = -t307 * t358 + t310 * t429;
t376 = t253 * pkin(2) + t392;
t374 = t301 * t429;
t370 = t301 * t390;
t367 = t401 / 0.2e1;
t48 = -t97 * mrSges(5,1) + t98 * mrSges(5,2);
t9 = -t29 * mrSges(6,1) + t28 * mrSges(6,2);
t359 = -pkin(1) * t307 + pkin(7) * t374;
t35 = -t300 * t77 + t302 * t86;
t68 = t302 * t137 - t138 * t300;
t151 = -t305 * t224 + t225 * t309;
t251 = t306 * t344 + t307 * t310;
t186 = t251 * t309 - t305 * t374;
t351 = mrSges(3,3) * t371;
t346 = -t251 * pkin(2) + t359;
t139 = pkin(3) * t398 - t151;
t343 = t248 * mrSges(4,1) + t249 * mrSges(4,2);
t338 = Ifges(4,1) * t309 - t420;
t337 = Ifges(5,1) * t302 - t418;
t336 = -Ifges(4,2) * t305 + t419;
t335 = -Ifges(5,2) * t300 + t417;
t334 = Ifges(4,5) * t309 - Ifges(4,6) * t305;
t333 = Ifges(5,5) * t302 - Ifges(5,6) * t300;
t181 = t249 * t302 - t300 * t398;
t49 = pkin(4) * t248 - pkin(9) * t181 + t68;
t180 = -t249 * t300 - t302 * t398;
t59 = pkin(9) * t180 + t69;
t17 = -t304 * t59 + t308 * t49;
t18 = t304 * t49 + t308 * t59;
t111 = t180 * t308 - t181 * t304;
t112 = t180 * t304 + t181 * t308;
t88 = -t224 * t386 - t225 * t387 + t235 * t309 - t305 * t237;
t325 = t197 * (mrSges(4,1) * t305 + mrSges(4,2) * t309);
t185 = t251 * t305 + t309 * t374;
t189 = t253 * t305 - t307 * t399;
t320 = -g(1) * t189 - g(2) * t185 - g(3) * t248;
t85 = -pkin(3) * t370 - t88;
t314 = Ifges(3,6) * t408 + (t310 * Ifges(3,2) + t422) * t301;
t313 = t301 * t329 * (Ifges(3,5) * t310 - Ifges(3,6) * t306);
t281 = Ifges(3,4) * t371;
t263 = t375 * t305;
t254 = (-mrSges(3,1) * t310 + mrSges(3,2) * t306) * t301;
t252 = t306 * t429 + t307 * t357;
t236 = t393 * qJD(1);
t232 = -mrSges(3,2) * t329 + t351;
t230 = t331 * t305;
t229 = t262 * t305;
t210 = -pkin(8) * t402 + t259;
t205 = Ifges(4,4) * t212;
t196 = Ifges(3,1) * t372 + Ifges(3,5) * t329 + t281;
t195 = Ifges(3,6) * qJD(2) + qJD(1) * t314;
t190 = t253 * t309 + t305 * t400;
t174 = mrSges(4,1) * t272 - mrSges(4,3) * t213;
t173 = -mrSges(4,2) * t272 - mrSges(4,3) * t212;
t160 = t183 * t302 + t300 * t370;
t159 = -t183 * t300 + t302 * t370;
t135 = t190 * t297 + t252 * t296;
t134 = -t190 * t296 + t252 * t297;
t119 = t213 * Ifges(4,1) - t205 + t497;
t117 = t213 * Ifges(4,5) - t212 * Ifges(4,6) + t272 * Ifges(4,3);
t114 = mrSges(5,1) * t212 - mrSges(5,3) * t168;
t113 = -mrSges(5,2) * t212 + mrSges(5,3) * t355;
t101 = -mrSges(4,2) * t228 - mrSges(4,3) * t133;
t100 = mrSges(4,1) * t228 - mrSges(4,3) * t132;
t99 = -mrSges(5,1) * t355 + mrSges(5,2) * t168;
t96 = -pkin(4) * t407 + t126;
t95 = -pkin(4) * t180 + t139;
t76 = -pkin(4) * t355 + t107;
t72 = mrSges(6,1) * t206 - mrSges(6,3) * t92;
t71 = -mrSges(6,2) * t206 + mrSges(6,3) * t512;
t70 = mrSges(4,1) * t133 + mrSges(4,2) * t132;
t65 = t132 * Ifges(4,4) - t133 * Ifges(4,2) + t228 * Ifges(4,6);
t62 = mrSges(5,1) * t133 - mrSges(5,3) * t98;
t61 = -mrSges(5,2) * t133 + mrSges(5,3) * t97;
t60 = -pkin(4) * t159 + t85;
t51 = -qJD(5) * t112 + t159 * t308 - t160 * t304;
t50 = qJD(5) * t111 + t159 * t304 + t160 * t308;
t47 = -mrSges(6,1) * t512 + mrSges(6,2) * t92;
t30 = pkin(9) * t159 + t36;
t24 = -pkin(4) * t97 + t46;
t21 = pkin(4) * t182 - pkin(9) * t160 + t35;
t16 = -mrSges(6,2) * t124 + mrSges(6,3) * t29;
t15 = mrSges(6,1) * t124 - mrSges(6,3) * t28;
t4 = -qJD(5) * t18 + t21 * t308 - t30 * t304;
t3 = qJD(5) * t17 + t21 * t304 + t30 * t308;
t6 = [(-m(3) * t233 - t472) * t238 - t476 * t385 + (-Ifges(4,6) * t432 + Ifges(5,3) * t442 - Ifges(4,2) * t443 - Ifges(4,4) * t444 + Ifges(6,3) * t445 + Ifges(5,5) * t447 + Ifges(5,6) * t448 - t55 * mrSges(4,3) - t14 * mrSges(5,2) + t13 * mrSges(5,1) - t65 / 0.2e1 + t506 / 0.2e1 + Ifges(6,6) * t462 + Ifges(6,5) * t463 + t511) * t248 + (t301 * t196 + t467 * qJD(1) * (-Ifges(3,2) * t306 + t421)) * t389 / 0.2e1 + (-t126 * t370 + t183 * t197 + t398 * t55) * mrSges(4,2) + (-t13 * t181 + t14 * t180 + t159 * t58 - t160 * t57) * mrSges(5,3) + (-m(3) * t359 - m(6) * t346 + t307 * mrSges(2,1) + t251 * mrSges(3,1) + t405 * mrSges(6,1) + mrSges(2,2) * t429 + t404 * mrSges(6,2) - mrSges(3,3) * t374 + t507 * (-pkin(8) * t250 + t346) - (t468 + t513) * t186 + t469 * t250 - t470 * t185) * g(1) + (t1 * t111 - t11 * t50 - t112 * t2 + t12 * t51) * mrSges(6,3) + t240 * (Ifges(3,5) * t408 + (t306 * Ifges(3,1) + t421) * t301) / 0.2e1 + (t313 / 0.2e1 + t117 * t367) * qJD(2) + (-t171 * t408 - t240 * t427 - t284 * t393) * mrSges(3,2) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t467 + t171 * t393 + t172 * t255 + t236 * t237) + (t172 * t408 - t239 * t427 + t255 * t284) * mrSges(3,1) + (t171 * t398 - t172 * t401 - t233 * t369 - t236 * t370 - t239 * t393 - t240 * t255) * mrSges(3,3) + (Ifges(3,1) * t240 - Ifges(3,4) * t239 + Ifges(3,5) * t284) * t367 + (Ifges(3,4) * t240 - Ifges(3,2) * t239 + Ifges(3,6) * t284) * t398 / 0.2e1 - t195 * t370 / 0.2e1 + t125 * (mrSges(4,1) * t370 - mrSges(4,3) * t183) + (Ifges(6,5) * t50 + Ifges(6,6) * t51) * t437 + (Ifges(6,5) * t112 + Ifges(6,6) * t111) * t445 + t183 * t119 / 0.2e1 + t46 * (-mrSges(5,1) * t180 + mrSges(5,2) * t181) + t148 * t343 + t272 * (Ifges(4,5) * t183 + Ifges(4,3) * t370) / 0.2e1 + (Ifges(4,5) * t249 - Ifges(4,3) * t398) * t432 + t284 * (Ifges(3,3) * t408 + (Ifges(3,5) * t306 + Ifges(3,6) * t310) * t301) / 0.2e1 - t254 * t381 + (Ifges(6,4) * t50 + Ifges(6,2) * t51) * t451 + (Ifges(6,4) * t112 + Ifges(6,2) * t111) * t462 + m(4) * (t125 * t88 + t126 * t87 + t148 * t223 + t151 * t56 + t152 * t55) + (Ifges(5,5) * t160 + Ifges(5,6) * t159) * t435 + (Ifges(5,5) * t181 + Ifges(5,6) * t180) * t442 + (t518 / 0.2e1 + t517 / 0.2e1 - Ifges(4,4) * t433 + Ifges(5,3) * t435 - Ifges(4,2) * t436 + t509 + t519) * t182 - t378 * t398 / 0.2e1 + t56 * (-mrSges(4,1) * t398 - t249 * mrSges(4,3)) + t355 * (Ifges(5,4) * t160 + Ifges(5,2) * t159) / 0.2e1 + (Ifges(5,4) * t181 + Ifges(5,2) * t180) * t448 + t408 * t377 / 0.2e1 + (Ifges(4,4) * t183 + Ifges(4,6) * t370) * t436 + (Ifges(4,4) * t249 - Ifges(4,6) * t398) * t443 - t239 * t314 / 0.2e1 + (Ifges(6,1) * t50 + Ifges(6,4) * t51) * t449 + (Ifges(6,1) * t112 + Ifges(6,4) * t111) * t463 + t87 * t173 + t88 * t174 + t159 * t79 / 0.2e1 + t107 * (-mrSges(5,1) * t159 + mrSges(5,2) * t160) + t151 * t100 + t168 * (Ifges(5,1) * t160 + Ifges(5,4) * t159) / 0.2e1 + t152 * t101 + (Ifges(5,1) * t181 + Ifges(5,4) * t180) * t447 + t139 * t48 + t24 * (-mrSges(6,1) * t111 + mrSges(6,2) * t112) + t36 * t113 + t35 * t114 + t160 * t453 + t223 * t70 + t237 * t232 + (-m(3) * t392 - m(6) * t376 - mrSges(2,1) * t429 - t253 * mrSges(3,1) - t135 * mrSges(6,1) + t307 * mrSges(2,2) - t134 * mrSges(6,2) - mrSges(3,3) * t400 + t507 * (pkin(8) * t252 + t376) + t468 * t190 - t469 * t252 + t470 * t189) * g(2) + (Ifges(4,1) * t183 + Ifges(4,5) * t370) * t433 + (Ifges(4,1) * t249 - Ifges(4,5) * t398) * t444 + m(6) * (t1 * t18 + t11 * t4 + t12 * t3 + t17 * t2 + t24 * t95 + t60 * t76) + m(5) * (t107 * t85 + t13 * t68 + t139 * t46 + t14 * t69 + t35 * t57 + t36 * t58) + t249 * t455 + t50 * t456 + t51 * t458 + t181 * t460 + t180 * t461 + t112 * t464 + t111 * t465 + t17 * t15 + t18 * t16 + t60 * t47 + t68 * t62 + t69 * t61 + t3 * t71 + t4 * t72 + t76 * (-mrSges(6,1) * t51 + mrSges(6,2) * t50) + Ifges(2,3) * qJDD(1) + t95 * t9 + t85 * t99; (-t415 + t509) * t387 + (t352 + t472) * t236 + t478 * t386 + t480 * mrSges(4,3) + t481 * t340 + t485 * t47 + (-mrSges(6,1) * t483 + mrSges(6,3) * t488) * t11 + (mrSges(6,1) * t489 - mrSges(6,2) * t488) * t76 + (mrSges(6,2) * t483 - mrSges(6,3) * t489) * t12 + (t302 * t80 + t119) * t386 / 0.2e1 + (t212 * (Ifges(4,6) * t306 + t310 * t336) + t306 * t195) * t391 / 0.2e1 + (-t126 * (-mrSges(4,3) * t305 * t310 - mrSges(4,2) * t306) - t125 * (mrSges(4,1) * t306 - mrSges(4,3) * t395)) * t391 + (-pkin(2) * t148 - t125 * t161 - t126 * t162) * m(4) + (-t107 * t141 + t13 * t210 + t14 * t211 + t492 * t58 + t493 * t57) * m(5) + (((-t125 * t309 - t126 * t305) * qJD(3) + t480) * m(4) + (-t100 + t48) * t305 + t309 * t101 + t481 * m(5)) * pkin(8) + (t250 * t514 + t251 * t471) * g(2) + (t252 * t514 + t253 * t471) * g(1) + (-t141 + t382) * t99 + (-t232 + t351) * t233 + (-t382 - t161) * t174 + (-t383 - t162) * t173 + (t58 * (-mrSges(5,2) * t305 - mrSges(5,3) * t402) + t57 * (mrSges(5,1) * t305 - mrSges(5,3) * t396) + t325) * qJD(3) + (t168 * (Ifges(5,5) * t305 + t309 * t337) + t355 * (Ifges(5,6) * t305 + t309 * t335) + t272 * t334 + t213 * t338) * qJD(3) / 0.2e1 + t1 * (mrSges(6,2) * t309 - mrSges(6,3) * t229) + (-Ifges(6,5) * t230 - Ifges(6,6) * t229 - Ifges(6,3) * t309) * t445 + t24 * (mrSges(6,1) * t229 - mrSges(6,2) * t230) + (-Ifges(6,4) * t230 - Ifges(6,2) * t229 - Ifges(6,6) * t309) * t462 + (-Ifges(6,1) * t230 - Ifges(6,4) * t229 - Ifges(6,5) * t309) * t463 + t2 * (-mrSges(6,1) * t309 + mrSges(6,3) * t230) + (Ifges(6,5) * t163 + Ifges(6,6) * t164) * t437 + (Ifges(6,4) * t163 + Ifges(6,2) * t164) * t451 - t325 * t371 - t107 * (-mrSges(5,1) * t199 + mrSges(5,2) * t200) - t200 * t80 / 0.2e1 - t148 * t342 - t32 * t403 / 0.2e1 + t14 * (mrSges(5,2) * t309 - mrSges(5,3) * t403) + t13 * (-mrSges(5,1) * t309 - mrSges(5,3) * t397) + t210 * t62 + t211 * t61 + (-t313 / 0.2e1 + t476 * qJD(1)) * qJD(1) + (Ifges(6,1) * t163 + Ifges(6,4) * t164) * t449 + t118 * t349 / 0.2e1 - t58 * (-mrSges(5,2) * t349 + t199 * mrSges(5,3)) - t57 * (mrSges(5,1) * t349 - t200 * mrSges(5,3)) + (Ifges(4,5) * t305 + Ifges(4,6) * t309) * t432 + (Ifges(5,5) * t200 + Ifges(5,6) * t199 + Ifges(5,3) * t349) * t436 + (Ifges(6,5) * t129 + Ifges(6,6) * t128 + Ifges(6,3) * t349) * t438 + (Ifges(5,1) * t200 + Ifges(5,4) * t199 + Ifges(5,5) * t349) * t440 + (Ifges(5,4) * t200 + Ifges(5,2) * t199 + Ifges(5,6) * t349) * t441 + (-Ifges(5,3) * t309 + t305 * t333) * t442 + (Ifges(4,2) * t309 + t420) * t443 + (Ifges(4,1) * t305 + t419) * t444 + (-Ifges(5,5) * t309 + t305 * t337) * t447 + (Ifges(5,3) * t305 + t309 * t333) * t388 / 0.2e1 - t336 * t388 / 0.2e1 - t171 * mrSges(3,2) + t172 * mrSges(3,1) - (t272 * (Ifges(4,3) * t306 + t310 * t334) + t213 * (Ifges(4,5) * t306 + t310 * t338) + t306 * t117 + (-Ifges(3,2) * t372 + t309 * t119 + t305 * t500 + t196 + t281) * t310) * t391 / 0.2e1 + t108 * t15 + t109 * t16 + (-Ifges(5,6) * t309 + t305 * t335) * t448 + (Ifges(6,1) * t129 + Ifges(6,4) * t128 + Ifges(6,5) * t349) * t450 + (Ifges(6,4) * t129 + Ifges(6,2) * t128 + Ifges(6,6) * t349) * t452 + t492 * t113 + t493 * t114 + (t254 - t477 * t394 + (t525 * t310 + (-m(6) * t426 - t339 - t524) * t306) * t301) * g(3) + t263 * t9 + t501 * t71 + t502 * t72 + (t1 * t109 + t108 * t2 + t11 * t502 + t12 * t501 + t24 * t263 + t485 * t76) * m(6) - t506 * t309 / 0.2e1 + t199 * t454 + t305 * t455 + t377 + t309 * t65 / 0.2e1 + t163 * t456 + t129 * t457 + t164 * t458 + t128 * t459 + t397 * t460 - t230 * t464 - t229 * t465 - pkin(2) * t70; t418 * t448 + t417 * t447 + (t185 * t475 + t186 * t470) * g(2) + (t189 * t475 + t190 * t470) * g(1) + (-t406 * t57 - t407 * t58 + t479) * mrSges(5,3) + (-pkin(3) * t46 + (-t300 * t57 + t302 * t58) * qJD(4) + t479 * qJ(4) - t107 * t126 - t57 * t73 - t58 * t74) * m(5) + (t248 * t484 + t249 * t482 + t343) * g(3) + (mrSges(6,1) * t487 - mrSges(6,2) * t486) * t76 + (Ifges(6,5) * t131 + Ifges(6,6) * t130) * t438 + (-t205 + t119) * t435 + (-t1 * t331 + t11 * t486 - t12 * t487 - t2 * t262) * mrSges(6,3) + (Ifges(6,5) * t262 - Ifges(6,6) * t331) * t445 + t24 * (mrSges(6,1) * t331 + mrSges(6,2) * t262) + (Ifges(6,4) * t262 - Ifges(6,2) * t331) * t462 + (Ifges(6,1) * t262 - Ifges(6,4) * t331) * t463 - t331 * t465 + (-Ifges(6,5) * t246 - Ifges(6,6) * t247) * t437 + (-Ifges(6,1) * t246 - Ifges(6,4) * t247) * t449 + (-Ifges(6,4) * t246 - Ifges(6,2) * t247) * t451 + t193 * t15 + t194 * t16 + t46 * t341 + (t174 - t99) * t126 + (-qJ(4) * t62 - qJD(4) * t114 + t460 + t510) * t300 + t118 * t433 + (Ifges(5,5) * t440 + Ifges(6,5) * t450 - Ifges(4,2) * t435 + Ifges(5,6) * t441 + Ifges(6,6) * t452 + Ifges(5,3) * t436 + Ifges(6,3) * t438 - t519) * t213 + t378 + (Ifges(5,2) * t448 + Ifges(5,6) * t442 + qJ(4) * t61 + qJD(4) * t113 + t461) * t302 - t125 * t173 - t74 * t113 - t73 * t114 + (Ifges(6,4) * t131 + Ifges(6,2) * t130) * t452 + t406 * t453 + t494 * t72 + t495 * t71 + (t1 * t194 + t11 * t494 + t12 * t495 + t193 * t2 - t24 * t295 - t76 * t96) * m(6) + (t107 * t340 - t333 * t436 - t337 * t440 - t335 * t441 + t197 * mrSges(4,2) + t497 / 0.2e1 + t478) * t212 - (-Ifges(4,1) * t212 - t412 + t500) * t213 / 0.2e1 - t295 * t9 + (Ifges(6,1) * t131 + Ifges(6,4) * t130) * t450 - t246 * t456 + t131 * t457 - t247 * t458 + t130 * t459 + t262 * t464 - pkin(3) * t48 - t55 * mrSges(4,2) + t56 * mrSges(4,1) - t96 * t47; -t355 * t113 + t168 * t114 - t512 * t71 + t92 * t72 + t48 + t9 + (t11 * t92 - t12 * t512 + t24 + t320) * m(6) + (t168 * t57 - t355 * t58 + t320 + t46) * m(5); -t76 * (mrSges(6,1) * t92 + mrSges(6,2) * t512) + (Ifges(6,1) * t512 - t428) * t450 + t38 * t449 + (Ifges(6,5) * t512 - Ifges(6,6) * t92) * t438 + t12 * t72 - t11 * t71 - g(1) * (mrSges(6,1) * t134 - mrSges(6,2) * t135) - g(2) * ((-t186 * t296 + t404) * mrSges(6,1) + (-t186 * t297 - t405) * mrSges(6,2)) - g(3) * ((-t249 * t296 - t297 * t398) * mrSges(6,1) + (-t249 * t297 + t296 * t398) * mrSges(6,2)) + (t11 * t512 + t12 * t92) * mrSges(6,3) + t5 + (-Ifges(6,2) * t92 + t39 + t89) * t452 + t511;];
tau = t6;
