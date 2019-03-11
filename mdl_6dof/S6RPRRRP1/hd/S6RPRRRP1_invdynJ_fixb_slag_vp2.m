% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP1
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:32
% EndTime: 2019-03-09 05:55:59
% DurationCPUTime: 16.75s
% Computational Cost: add. (9777->672), mult. (20187->854), div. (0->0), fcn. (13492->14), ass. (0->294)
t508 = -mrSges(6,3) - mrSges(7,2);
t505 = Ifges(6,1) + Ifges(7,1);
t509 = -Ifges(6,4) + Ifges(7,5);
t481 = Ifges(7,4) + Ifges(6,5);
t480 = Ifges(7,2) + Ifges(6,3);
t504 = Ifges(6,6) - Ifges(7,6);
t274 = sin(qJ(3));
t278 = cos(qJ(3));
t364 = qJD(1) * qJD(3);
t214 = qJDD(1) * t278 - t274 * t364;
t215 = qJDD(1) * t274 + t278 * t364;
t273 = sin(qJ(4));
t277 = cos(qJ(4));
t212 = t273 * t274 - t277 * t278;
t290 = t212 * qJD(4);
t131 = -qJD(1) * t290 + t214 * t273 + t215 * t277;
t267 = qJDD(3) + qJDD(4);
t272 = sin(qJ(5));
t276 = cos(qJ(5));
t213 = t273 * t278 + t274 * t277;
t204 = t213 * qJD(1);
t363 = qJD(3) + qJD(4);
t176 = t272 * t204 - t276 * t363;
t367 = qJD(5) * t176;
t72 = t276 * t131 + t272 * t267 - t367;
t448 = t72 / 0.2e1;
t177 = t276 * t204 + t272 * t363;
t73 = qJD(5) * t177 + t272 * t131 - t276 * t267;
t446 = t73 / 0.2e1;
t291 = t213 * qJD(4);
t132 = -qJD(1) * t291 + t214 * t277 - t215 * t273;
t129 = qJDD(5) - t132;
t444 = t129 / 0.2e1;
t484 = mrSges(6,1) + mrSges(7,1);
t483 = mrSges(6,2) - mrSges(7,3);
t203 = t212 * qJD(1);
t200 = qJD(5) + t203;
t413 = mrSges(7,2) * t176;
t133 = mrSges(7,3) * t200 - t413;
t411 = mrSges(6,3) * t176;
t134 = -mrSges(6,2) * t200 - t411;
t507 = -t134 - t133;
t410 = mrSges(6,3) * t177;
t135 = mrSges(6,1) * t200 - t410;
t412 = mrSges(7,2) * t177;
t136 = -mrSges(7,1) * t200 + t412;
t378 = -t135 + t136;
t270 = sin(pkin(10));
t246 = pkin(1) * t270 + pkin(7);
t229 = t246 * qJD(1);
t328 = pkin(8) * qJD(1) + t229;
t372 = qJD(2) * t274;
t179 = t278 * t328 + t372;
t174 = t277 * t179;
t262 = t278 * qJD(2);
t178 = -t274 * t328 + t262;
t175 = qJD(3) * pkin(3) + t178;
t108 = t273 * t175 + t174;
t101 = pkin(9) * t363 + t108;
t271 = cos(pkin(10));
t247 = -pkin(1) * t271 - pkin(2);
t265 = t278 * pkin(3);
t225 = t247 - t265;
t205 = t225 * qJD(1);
t130 = pkin(4) * t203 - pkin(9) * t204 + t205;
t191 = t229 * t278 + t372;
t227 = t246 * qJDD(1);
t151 = -qJD(3) * t191 + t278 * qJDD(2) - t227 * t274;
t118 = qJDD(3) * pkin(3) - pkin(8) * t215 + t151;
t371 = qJD(3) * t274;
t150 = qJD(3) * t262 + t274 * qJDD(2) + t278 * t227 - t229 * t371;
t122 = pkin(8) * t214 + t150;
t368 = qJD(4) * t277;
t369 = qJD(4) * t273;
t30 = t273 * t118 + t277 * t122 + t175 * t368 - t179 * t369;
t25 = pkin(9) * t267 + t30;
t365 = qJD(5) * t276;
t366 = qJD(5) * t272;
t228 = t247 * qJDD(1);
t180 = -pkin(3) * t214 + t228;
t41 = -pkin(4) * t132 - pkin(9) * t131 + t180;
t6 = -t101 * t366 + t130 * t365 + t276 * t25 + t272 * t41;
t2 = qJ(6) * t129 + qJD(6) * t200 + t6;
t45 = t101 * t276 + t130 * t272;
t7 = -qJD(5) * t45 - t25 * t272 + t276 * t41;
t4 = -pkin(5) * t129 + qJDD(6) - t7;
t322 = t2 * t276 + t272 * t4;
t44 = -t101 * t272 + t130 * t276;
t37 = -pkin(5) * t200 + qJD(6) - t44;
t38 = qJ(6) * t200 + t45;
t506 = t37 * t365 - t38 * t366 + t322;
t503 = t129 * t481 + t505 * t72 + t509 * t73;
t502 = -t176 * t504 + t177 * t481 + t480 * t200;
t172 = Ifges(6,4) * t176;
t404 = Ifges(7,5) * t176;
t478 = t177 * t505 + t481 * t200 - t172 + t404;
t501 = t363 * Ifges(5,5);
t500 = t363 * Ifges(5,6);
t399 = t204 * mrSges(5,3);
t472 = -mrSges(5,1) * t363 + mrSges(6,1) * t176 + mrSges(6,2) * t177 + t399;
t119 = t178 * t273 + t174;
t499 = -pkin(3) * t369 + t119;
t391 = t203 * t276;
t392 = t203 * t272;
t498 = -qJD(6) * t272 + (-t365 - t391) * qJ(6) + (t366 + t392) * pkin(5);
t269 = qJ(3) + qJ(4);
t263 = sin(t269);
t497 = t508 * t263;
t304 = pkin(5) * t276 + qJ(6) * t272;
t226 = -pkin(4) - t304;
t315 = -t276 * mrSges(7,1) - t272 * mrSges(7,3);
t414 = mrSges(6,1) * t276;
t466 = (-m(7) * t226 - t315 + t414) * t263;
t496 = -t272 * t504 + t276 * t481;
t403 = Ifges(7,5) * t272;
t406 = Ifges(6,4) * t272;
t495 = -t276 * t505 - t403 + t406;
t321 = -t7 * t272 + t276 * t6;
t494 = -t44 * t365 - t45 * t366 + t321;
t33 = -mrSges(7,2) * t73 + mrSges(7,3) * t129;
t34 = mrSges(6,1) * t129 - mrSges(6,3) * t72;
t35 = -t129 * mrSges(7,1) + t72 * mrSges(7,2);
t36 = -mrSges(6,2) * t129 - mrSges(6,3) * t73;
t492 = (-t34 + t35) * t272 + (t33 + t36) * t276;
t173 = t273 * t179;
t107 = t277 * t175 - t173;
t100 = -pkin(4) * t363 - t107;
t314 = t272 * mrSges(7,1) - t276 * mrSges(7,3);
t316 = mrSges(6,1) * t272 + mrSges(6,2) * t276;
t46 = t176 * pkin(5) - t177 * qJ(6) + t100;
t491 = t100 * t316 + t314 * t46;
t490 = Ifges(7,5) * t448 + Ifges(7,6) * t444 - t72 * Ifges(6,4) / 0.2e1 - t129 * Ifges(6,6) / 0.2e1 + (Ifges(7,3) + Ifges(6,2)) * t446;
t489 = m(6) * ((-t272 * t45 - t276 * t44) * qJD(5) + t321) + m(7) * ((-t272 * t38 + t276 * t37) * qJD(5) + t322) + t507 * t366 + t378 * t365 + t492;
t401 = t177 * Ifges(6,4);
t94 = -t176 * Ifges(6,2) + t200 * Ifges(6,6) + t401;
t488 = -t94 / 0.2e1;
t487 = m(3) + m(4);
t486 = -m(6) - m(7);
t485 = t214 / 0.2e1;
t114 = mrSges(5,1) * t267 - mrSges(5,3) * t131;
t28 = mrSges(6,1) * t73 + mrSges(6,2) * t72;
t476 = t28 - t114;
t475 = -t108 + t498;
t474 = t498 - t499;
t153 = pkin(4) * t212 - pkin(9) * t213 + t225;
t415 = pkin(8) + t246;
t207 = t415 * t274;
t208 = t415 * t278;
t160 = -t207 * t273 + t208 * t277;
t470 = t272 * t153 + t276 * t160;
t469 = -t277 * t207 - t208 * t273;
t264 = cos(t269);
t468 = t264 * mrSges(5,1) - t263 * mrSges(5,2);
t375 = t264 * pkin(4) + t263 * pkin(9);
t427 = pkin(4) * t263;
t429 = pkin(3) * t274;
t467 = m(7) * t429 - m(6) * (-t427 - t429) + t466;
t268 = qJ(1) + pkin(10);
t259 = cos(t268);
t383 = t259 * t264;
t223 = pkin(9) * t383;
t382 = t263 * t272;
t360 = mrSges(6,2) * t382;
t465 = -m(7) * t223 - t259 * t360 + t383 * t508;
t258 = sin(t268);
t385 = t258 * t264;
t221 = pkin(9) * t385;
t464 = -m(7) * t221 - t258 * t360 + t385 * t508;
t168 = -qJD(3) * t212 - t290;
t394 = t168 * t272;
t296 = t213 * t365 + t394;
t462 = t129 * t480 + t481 * t72 - t504 * t73;
t461 = t150 * t278 - t151 * t274;
t460 = g(1) * t259 + g(2) * t258;
t458 = -m(5) + t486;
t233 = -mrSges(4,1) * t278 + mrSges(4,2) * t274;
t457 = -m(4) * pkin(2) - mrSges(3,1) + t233 - t468;
t380 = t264 * t276;
t381 = t264 * t272;
t456 = -t380 * t484 + t381 * t483 - t468 + t497;
t331 = qJD(3) * t415;
t196 = t274 * t331;
t197 = t278 * t331;
t85 = qJD(4) * t469 - t196 * t277 - t197 * t273;
t169 = qJD(3) * t213 + t291;
t358 = pkin(3) * t371;
t99 = pkin(4) * t169 - pkin(9) * t168 + t358;
t13 = -qJD(5) * t470 - t272 * t85 + t276 * t99;
t453 = m(7) * pkin(5) + t484;
t452 = m(7) * qJ(6) - t483;
t451 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t450 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t447 = -t73 / 0.2e1;
t442 = -t176 / 0.2e1;
t441 = t176 / 0.2e1;
t440 = -t177 / 0.2e1;
t439 = t177 / 0.2e1;
t438 = -t200 / 0.2e1;
t436 = t203 / 0.2e1;
t434 = t204 / 0.2e1;
t432 = t272 / 0.2e1;
t275 = sin(qJ(1));
t431 = pkin(1) * t275;
t430 = pkin(3) * t273;
t428 = pkin(3) * t277;
t426 = pkin(5) * t204;
t421 = g(3) * t263;
t279 = cos(qJ(1));
t266 = t279 * pkin(1);
t409 = Ifges(4,4) * t274;
t408 = Ifges(4,4) * t278;
t407 = Ifges(5,4) * t204;
t405 = Ifges(6,4) * t276;
t402 = Ifges(7,5) * t276;
t400 = t203 * mrSges(5,3);
t393 = t168 * t276;
t384 = t259 * t263;
t163 = pkin(4) * t204 + pkin(9) * t203;
t52 = t276 * t107 + t272 * t163;
t120 = t178 * t277 - t173;
t374 = qJD(1) * t274;
t359 = pkin(3) * t374;
t145 = t163 + t359;
t50 = t276 * t120 + t272 * t145;
t373 = qJD(1) * t278;
t370 = qJD(3) * t278;
t356 = pkin(3) * t368;
t354 = mrSges(4,3) * t374;
t353 = mrSges(4,3) * t373;
t171 = Ifges(7,5) * t177;
t91 = t200 * Ifges(7,6) + t176 * Ifges(7,3) + t171;
t348 = t91 * t432;
t335 = -t366 / 0.2e1;
t334 = t365 / 0.2e1;
t327 = t272 * t356;
t326 = t276 * t356;
t325 = pkin(5) * t380 + qJ(6) * t381 + t375;
t254 = t265 + pkin(2);
t280 = -pkin(8) - pkin(7);
t320 = t259 * t254 - t258 * t280 + t266;
t319 = mrSges(4,1) * t274 + mrSges(4,2) * t278;
t317 = mrSges(5,1) * t263 + mrSges(5,2) * t264;
t311 = t278 * Ifges(4,2) + t409;
t310 = -Ifges(6,2) * t272 + t405;
t308 = Ifges(4,5) * t278 - Ifges(4,6) * t274;
t306 = Ifges(7,3) * t272 + t402;
t303 = pkin(5) * t272 - qJ(6) * t276;
t51 = -t107 * t272 + t163 * t276;
t49 = -t120 * t272 + t145 * t276;
t79 = t153 * t276 - t160 * t272;
t31 = t118 * t277 - t273 * t122 - t175 * t369 - t179 * t368;
t295 = t213 * t366 - t393;
t294 = t247 * qJD(1) * t319;
t293 = t274 * (Ifges(4,1) * t278 - t409);
t12 = t153 * t365 - t160 * t366 + t272 * t99 + t276 * t85;
t26 = -pkin(4) * t267 - t31;
t86 = qJD(4) * t160 - t196 * t273 + t277 * t197;
t146 = -Ifges(5,2) * t203 + t407 + t500;
t198 = Ifges(5,4) * t203;
t147 = Ifges(5,1) * t204 - t198 + t501;
t9 = pkin(5) * t73 - qJ(6) * t72 - qJD(6) * t177 + t26;
t282 = (-t177 * t495 + t200 * t496) * qJD(5) / 0.2e1 + t503 * t432 + (Ifges(6,2) * t447 - Ifges(7,3) * t446 + t444 * t504 - t490) * t276 + (Ifges(6,6) * t441 + Ifges(7,6) * t442 - t44 * mrSges(6,1) + t37 * mrSges(7,1) + t45 * mrSges(6,2) - t38 * mrSges(7,3) - t205 * mrSges(5,1) + t500 / 0.2e1 - Ifges(5,2) * t436 + t481 * t440 + t480 * t438) * t204 + (t26 * mrSges(6,2) + t481 * t444 + t448 * t505) * t272 + (t348 + t491) * qJD(5) + (-t391 * t44 - t392 * t45 + t494) * mrSges(6,3) + (-t310 / 0.2e1 + t306 / 0.2e1) * t367 - t26 * t414 + (-t402 + t405) * t448 + (-t310 * t441 - t306 * t442 + t205 * mrSges(5,2) + t501 / 0.2e1 + t495 * t440 - t496 * t438 + t491) * t203 - (-Ifges(5,1) * t203 - t407 + t502) * t204 / 0.2e1 + t403 * t446 + t406 * t447 + (-t198 + t147) * t436 - t30 * mrSges(5,2) + t31 * mrSges(5,1) + t94 * t335 + (t91 / 0.2e1 + t488) * t392 + (t391 / 0.2e1 + t334) * t478 + t146 * t434 + t9 * t315 + Ifges(5,5) * t131 + Ifges(5,6) * t132 + (t37 * t391 - t38 * t392 + t506) * mrSges(7,2) + t108 * t399 - t107 * t400 + Ifges(5,3) * t267;
t255 = Ifges(4,4) * t373;
t253 = -pkin(4) - t428;
t232 = -qJD(3) * mrSges(4,2) + t353;
t230 = qJD(3) * mrSges(4,1) - t354;
t206 = t226 - t428;
t202 = Ifges(4,1) * t374 + Ifges(4,5) * qJD(3) + t255;
t201 = Ifges(4,6) * qJD(3) + qJD(1) * t311;
t195 = t204 * qJ(6);
t194 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t215;
t193 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t214;
t190 = -t229 * t274 + t262;
t187 = t258 * t272 + t259 * t380;
t186 = -t258 * t276 + t259 * t381;
t185 = t258 * t380 - t259 * t272;
t184 = t258 * t381 + t259 * t276;
t181 = -mrSges(5,2) * t363 - t400;
t162 = mrSges(5,1) * t203 + mrSges(5,2) * t204;
t115 = -mrSges(5,2) * t267 + mrSges(5,3) * t132;
t112 = mrSges(7,1) * t176 - mrSges(7,3) * t177;
t111 = pkin(5) * t177 + qJ(6) * t176;
t98 = t213 * t303 - t469;
t59 = -pkin(5) * t212 - t79;
t58 = qJ(6) * t212 + t470;
t48 = -t51 - t426;
t47 = t195 + t52;
t43 = -t49 - t426;
t42 = t195 + t50;
t27 = mrSges(7,1) * t73 - mrSges(7,3) * t72;
t22 = t303 * t168 + (qJD(5) * t304 - qJD(6) * t276) * t213 + t86;
t11 = -pkin(5) * t169 - t13;
t10 = qJ(6) * t169 + qJD(6) * t212 + t12;
t1 = [(t294 + t308 * qJD(3) / 0.2e1) * qJD(3) + m(5) * (t108 * t85 + t160 * t30 + t180 * t225 + t205 * t358) - (-m(5) * t31 + m(6) * t26 + t476) * t469 + (m(4) * t247 + t233) * t228 + (t180 * mrSges(5,1) - t30 * mrSges(5,3) - Ifges(5,4) * t131 - Ifges(5,2) * t132 - Ifges(5,6) * t267 + Ifges(6,6) * t447 + Ifges(7,6) * t446 + t444 * t480 + t448 * t481 + t450 + t462 / 0.2e1) * t212 + (t169 * t480 - t295 * t481 - t296 * t504) * t200 / 0.2e1 + (Ifges(4,1) * t215 + Ifges(4,4) * t485) * t274 + t478 * t393 / 0.2e1 + t215 * t408 / 0.2e1 + (-m(5) * t107 + m(6) * t100 + t472) * t86 + (-t230 * t370 - t232 * t371 + m(4) * ((-t190 * t278 - t191 * t274) * qJD(3) + t461) - t274 * t194 + t278 * t193) * t246 + (-t190 * t370 - t191 * t371 + t461) * mrSges(4,3) + t278 * (Ifges(4,4) * t215 + Ifges(4,2) * t214) / 0.2e1 + (t293 + t278 * (-Ifges(4,2) * t274 + t408)) * t364 / 0.2e1 + t502 * t169 / 0.2e1 + (-t107 * t168 - t108 * t169) * mrSges(5,3) + t470 * t36 + m(6) * (t12 * t45 + t13 * t44 + t470 * t6 + t7 * t79) + t22 * t112 + t98 * t27 + t79 * t34 + t58 * t33 + t59 * t35 + t202 * t370 / 0.2e1 - t201 * t371 / 0.2e1 + t162 * t358 + (mrSges(2,1) * t275 + mrSges(2,2) * t279 + t487 * t431 + t458 * (-t259 * t280 - t431) + t453 * t185 + t452 * t184 + t451 * t259 + (m(5) * t254 + t486 * (-t254 - t375) - t457 - t497) * t258) * g(1) + m(7) * (t10 * t38 + t11 * t37 + t2 * t58 + t22 * t46 + t4 * t59 + t9 * t98) + t363 * (Ifges(5,5) * t168 - Ifges(5,6) * t169) / 0.2e1 + (t481 * t169 - t505 * t295 + t296 * t509) * t439 + ((mrSges(7,2) * t4 - mrSges(6,3) * t7 + t503 / 0.2e1) * t276 + t478 * t335 + t180 * mrSges(5,2) - t31 * mrSges(5,3) + Ifges(5,1) * t131 + Ifges(5,4) * t132 + Ifges(5,5) * t267 + t26 * t316 + t306 * t446 + t310 * t447 + t314 * t9 + t334 * t91 - t495 * t448 + t496 * t444 + (-t2 * mrSges(7,2) - t6 * mrSges(6,3) + t490) * t272) * t213 + t168 * t348 + t45 * (-mrSges(6,2) * t169 - mrSges(6,3) * t296) + t38 * (-mrSges(7,2) * t296 + mrSges(7,3) * t169) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t271 - 0.2e1 * mrSges(3,2) * t270 + m(3) * (t270 ^ 2 + t271 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (-Ifges(7,5) * t295 + Ifges(7,6) * t169 + Ifges(7,3) * t296) * t441 + (-Ifges(6,4) * t295 - Ifges(6,2) * t296 + Ifges(6,6) * t169) * t442 + (Ifges(5,1) * t168 - Ifges(5,4) * t169) * t434 + t10 * t133 + t12 * t134 + t13 * t135 + t11 * t136 + t160 * t115 + t168 * t147 / 0.2e1 - t169 * t146 / 0.2e1 + t85 * t181 + t311 * t485 + t296 * t488 - t203 * (Ifges(5,4) * t168 - Ifges(5,2) * t169) / 0.2e1 + t205 * (mrSges(5,1) * t169 + mrSges(5,2) * t168) + (-m(5) * t320 - mrSges(2,1) * t279 + mrSges(2,2) * t275 + t508 * t384 + t486 * (pkin(4) * t383 + pkin(9) * t384 + t320) - t487 * t266 - t453 * t187 - t452 * t186 + t457 * t259 + t451 * t258) * g(2) + t225 * (-mrSges(5,1) * t132 + mrSges(5,2) * t131) + t247 * (-mrSges(4,1) * t214 + mrSges(4,2) * t215) + qJDD(3) * (Ifges(4,5) * t274 + Ifges(4,6) * t278) + t44 * (mrSges(6,1) * t169 + mrSges(6,3) * t295) + t37 * (-mrSges(7,1) * t169 - mrSges(7,2) * t295) + t46 * (mrSges(7,1) * t296 + mrSges(7,3) * t295) + t100 * (mrSges(6,1) * t296 - mrSges(6,2) * t295); m(3) * qJDD(2) + t274 * t193 + t278 * t194 + (-t230 * t274 + t232 * t278) * qJD(3) + (t27 + t476) * t212 + (t112 + t472) * t169 + (t272 * t378 - t276 * t507 + t181) * t168 + (t458 - t487) * g(3) + m(5) * (-t107 * t169 + t108 * t168 - t212 * t31) + m(7) * (t169 * t46 + t212 * t9 + t37 * t394 + t38 * t393) + m(6) * (t100 * t169 + t212 * t26 + t393 * t45 - t394 * t44) + m(4) * (t150 * t274 + t151 * t278 + (-t190 * t274 + t191 * t278) * qJD(3)) + (t115 + (t272 * t507 + t276 * t378) * qJD(5) + m(5) * t30 + m(7) * t506 + m(6) * t494 + t492) * t213; (-m(6) * (t265 + t375) - m(5) * t265 - m(7) * (t265 + t325) + t233 + t456) * g(3) + (t206 * t9 + (t272 * t37 + t276 * t38) * t356 - t37 * t43 - t38 * t42 + t474 * t46) * m(7) + t474 * t112 + (m(5) * t429 + t317 + t319) * t460 - t472 * t499 + (t253 * t26 + (t100 * t273 + (-t272 * t44 + t276 * t45) * t277) * qJD(4) * pkin(3) - g(1) * t223 - g(2) * t221 - t100 * t119 - t44 * t49 - t45 * t50) * m(6) + (t258 * t467 + t464) * g(2) + (t259 * t467 + t465) * g(1) + t489 * (pkin(9) + t430) + (t354 + t230) * t191 + (t353 - t232) * t190 - (-Ifges(4,2) * t374 + t202 + t255) * t373 / 0.2e1 + (t356 - t120) * t181 + (t327 - t43) * t136 + (-t327 - t49) * t135 + (t326 - t50) * t134 + (t326 - t42) * t133 + (-t294 - t293 * qJD(1) / 0.2e1) * qJD(1) + (t107 * t119 - t108 * t120 - t205 * t359 + (t273 * t30 + t277 * t31 + (-t107 * t273 + t108 * t277) * qJD(4)) * pkin(3)) * m(5) + t201 * t374 / 0.2e1 - t308 * t364 / 0.2e1 + t282 + t114 * t428 + t115 * t430 + Ifges(4,3) * qJDD(3) - t150 * mrSges(4,2) + t151 * mrSges(4,1) + t206 * t27 + Ifges(4,6) * t214 + Ifges(4,5) * t215 + t253 * t28 - t162 * t359; -pkin(4) * t28 - t107 * t181 - t47 * t133 - t52 * t134 - t51 * t135 - t48 * t136 + t226 * t27 + t282 + t460 * t317 + t475 * t112 - t472 * t108 + (t258 * t466 + t464) * g(2) + (t466 * t259 + t465) * g(1) + (t226 * t9 - t37 * t48 - t38 * t47 + t46 * t475) * m(7) + (-g(1) * (-pkin(4) * t384 + t223) - pkin(4) * t26 - g(2) * (-t258 * t427 + t221) - t100 * t108 - t44 * t51 - t45 * t52) * m(6) + (-m(6) * t375 - m(7) * t325 + t456) * g(3) + t489 * pkin(9); (-t176 * t505 + t171 - t401 + t91) * t440 + (t484 * t186 + t483 * t187) * g(1) + (t184 * t484 + t185 * t483) * g(2) + (-t176 * t481 - t177 * t504) * t438 + (t314 + t316) * t421 + (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t38 + t303 * t421 - g(2) * (-pkin(5) * t184 + qJ(6) * t185) - g(1) * (-pkin(5) * t186 + qJ(6) * t187) - t111 * t46) * m(7) + (-m(7) * t37 - t378 + t410) * t45 + (-Ifges(6,2) * t177 - t172 + t478) * t441 + (-m(7) * t38 - t411 + t507) * t44 - t111 * t112 + qJ(6) * t33 - pkin(5) * t35 + t450 + t94 * t439 + (Ifges(7,3) * t177 - t404) * t442 + qJD(6) * t133 - t46 * (mrSges(7,1) * t177 + mrSges(7,3) * t176) - t100 * (mrSges(6,1) * t177 - mrSges(6,2) * t176) + t38 * t412 + t37 * t413 + t462; t177 * t112 - t200 * t133 + (-g(1) * t186 - g(2) * t184 - g(3) * t382 + t46 * t177 - t38 * t200 + t4) * m(7) + t35;];
tau  = t1;
