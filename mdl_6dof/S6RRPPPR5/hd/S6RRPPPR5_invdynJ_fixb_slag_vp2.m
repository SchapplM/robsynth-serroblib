% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:21
% EndTime: 2019-03-09 08:22:01
% DurationCPUTime: 27.06s
% Computational Cost: add. (5173->760), mult. (11494->945), div. (0->0), fcn. (7439->8), ass. (0->365)
t516 = -Ifges(5,3) - Ifges(6,1);
t515 = Ifges(6,3) + Ifges(4,1);
t248 = cos(qJ(2));
t372 = qJD(1) * t248;
t207 = qJD(6) - t372;
t430 = t207 / 0.2e1;
t242 = sin(pkin(9));
t245 = sin(qJ(2));
t373 = qJD(1) * t245;
t355 = t242 * t373;
t243 = cos(pkin(9));
t364 = t243 * qJD(2);
t171 = t355 - t364;
t354 = t243 * t373;
t172 = qJD(2) * t242 + t354;
t244 = sin(qJ(6));
t247 = cos(qJ(6));
t101 = t171 * t247 + t172 * t244;
t439 = t101 / 0.2e1;
t337 = -t171 * t244 + t247 * t172;
t442 = t337 / 0.2e1;
t526 = Ifges(7,5) * t439 + Ifges(7,6) * t442 + Ifges(7,3) * t430;
t363 = qJD(1) * qJD(2);
t182 = qJDD(1) * t245 + t248 * t363;
t138 = -t243 * qJDD(2) + t182 * t242;
t438 = -t138 / 0.2e1;
t139 = qJDD(2) * t242 + t182 * t243;
t435 = t139 / 0.2e1;
t525 = -t171 / 0.2e1;
t524 = -t172 / 0.2e1;
t228 = t248 * qJDD(1);
t349 = t245 * t363;
t181 = -t228 + t349;
t432 = t181 / 0.2e1;
t523 = -t372 / 0.2e1;
t345 = t372 / 0.2e1;
t522 = t515 * t435;
t496 = mrSges(5,1) - mrSges(6,2);
t521 = -m(7) * pkin(8) - mrSges(4,3) - mrSges(7,3) - t496;
t437 = t138 / 0.2e1;
t433 = -t181 / 0.2e1;
t520 = -Ifges(4,4) + Ifges(6,5);
t493 = -Ifges(5,4) + Ifges(4,5);
t519 = Ifges(5,5) - Ifges(6,4);
t518 = Ifges(6,5) - Ifges(5,6);
t492 = Ifges(4,6) - Ifges(5,5);
t517 = Ifges(6,6) - Ifges(4,5);
t491 = -Ifges(4,3) - Ifges(5,1);
t177 = t242 * t244 - t243 * t247;
t283 = t177 * t248;
t136 = qJD(1) * t283;
t153 = t177 * qJD(6);
t513 = t136 - t153;
t407 = Ifges(5,6) * t243;
t409 = Ifges(6,5) * t243;
t512 = -t242 * t516 - t407 + t409;
t410 = Ifges(6,5) * t242;
t413 = Ifges(4,4) * t242;
t511 = t243 * t515 + t410 - t413;
t229 = t248 * qJD(5);
t218 = pkin(7) * t228;
t220 = pkin(7) * t373;
t142 = qJDD(2) * qJ(3) + t218 + (qJD(3) - t220) * qJD(2);
t367 = qJD(3) * t245;
t398 = qJDD(1) * pkin(1);
t84 = pkin(2) * t181 - qJ(3) * t182 - qJD(1) * t367 - t398;
t37 = -t242 * t142 + t243 * t84;
t322 = qJDD(4) - t37;
t418 = pkin(3) + qJ(5);
t264 = qJD(1) * t229 - t181 * t418 + t322;
t13 = pkin(4) * t139 + t264;
t31 = -pkin(3) * t181 + t322;
t510 = t31 * mrSges(5,1) - t13 * mrSges(6,2) + Ifges(5,2) * t435 + Ifges(5,6) * t438 + t437 * t520 + t522 + (Ifges(5,4) + t517) * t433;
t231 = t245 * qJ(3);
t238 = t248 * pkin(2);
t357 = -pkin(1) - t238;
t291 = t357 - t231;
t163 = t291 * qJD(1);
t221 = pkin(7) * t372;
t190 = qJD(2) * qJ(3) + t221;
t104 = t163 * t243 - t242 * t190;
t105 = t242 * t163 + t243 * t190;
t408 = Ifges(5,6) * t242;
t302 = -Ifges(5,2) * t243 + t408;
t412 = Ifges(4,4) * t243;
t308 = -Ifges(4,2) * t242 + t412;
t312 = -mrSges(5,2) * t242 - mrSges(5,3) * t243;
t317 = mrSges(6,1) * t243 - mrSges(6,3) * t242;
t318 = mrSges(4,1) * t242 + mrSges(4,2) * t243;
t326 = qJD(2) * pkin(2) - qJD(3) - t220;
t387 = t243 * t248;
t389 = t242 * t248;
t222 = pkin(3) * t372;
t72 = qJD(4) - t104 + t222;
t282 = qJ(5) * t372 + t72;
t41 = pkin(4) * t172 + t282;
t281 = qJ(4) * t172 + t326;
t481 = t171 * t418;
t49 = t281 - t481;
t76 = qJ(4) * t372 - t105;
t50 = -pkin(4) * t171 + qJD(5) - t76;
t69 = pkin(3) * t171 - t281;
t509 = t105 * (-mrSges(4,2) * t245 - mrSges(4,3) * t389) + t104 * (mrSges(4,1) * t245 - mrSges(4,3) * t387) + t76 * (mrSges(5,1) * t389 - mrSges(5,3) * t245) + t72 * (mrSges(5,1) * t387 + mrSges(5,2) * t245) + t50 * (mrSges(6,1) * t245 + mrSges(6,2) * t389) + t41 * (-mrSges(6,2) * t387 - mrSges(6,3) * t245) - (-t312 * t69 - t317 * t49 + t318 * t326) * t248 + (Ifges(4,6) * t245 + t248 * t308) * t525 + (Ifges(5,4) * t245 + t248 * t302) * t524;
t365 = qJD(4) * t248;
t38 = t243 * t142 + t242 * t84;
t27 = -t181 * qJ(4) + qJD(1) * t365 - t38;
t267 = qJDD(5) - t27;
t17 = -pkin(4) * t138 + t267;
t436 = -t139 / 0.2e1;
t508 = t27 * mrSges(5,1) + t17 * mrSges(6,2) + Ifges(4,4) * t436 + Ifges(4,2) * t437 + Ifges(4,6) * t433 + t432 * t519 + t435 * t518 + t438 * t516;
t313 = t243 * mrSges(5,2) - t242 * mrSges(5,3);
t176 = t242 * t247 + t243 * t244;
t314 = -mrSges(7,1) * t176 + mrSges(7,2) * t177;
t319 = mrSges(4,1) * t243 - mrSges(4,2) * t242;
t507 = t313 + t314 - t319;
t219 = Ifges(3,4) * t372;
t506 = Ifges(3,1) * t373 + Ifges(3,5) * qJD(2) + t219 + t243 * (t171 * t520 + t172 * t515 + t372 * t517) + t242 * (-t171 * t516 + t172 * t518 - t372 * t519);
t415 = Ifges(3,4) * t245;
t485 = t248 * Ifges(3,2);
t309 = t415 + t485;
t441 = pkin(4) + pkin(8);
t33 = t172 * t441 + t282;
t416 = pkin(5) + qJ(4);
t350 = t416 * t248;
t34 = -qJD(1) * t350 - t171 * t441 + qJD(5) + t105;
t8 = -t244 * t33 + t247 * t34;
t9 = t244 * t34 + t247 * t33;
t505 = t8 * mrSges(7,1) - t9 * mrSges(7,2) + t172 * t493 / 0.2e1 + t491 * t345 + Ifges(6,2) * t523 - Ifges(3,6) * qJD(2) / 0.2e1 + Ifges(6,6) * t524 - qJD(1) * t309 / 0.2e1 + (t492 + Ifges(6,4)) * t525 + t526;
t411 = Ifges(7,4) * t101;
t29 = Ifges(7,2) * t337 + Ifges(7,6) * t207 + t411;
t446 = t29 / 0.2e1;
t503 = t248 * g(3);
t134 = mrSges(4,2) * t372 - mrSges(4,3) * t171;
t131 = mrSges(5,1) * t171 + mrSges(5,3) * t372;
t133 = -mrSges(6,1) * t372 + mrSges(6,2) * t171;
t382 = t133 - t131;
t502 = t134 + t382;
t366 = qJD(4) * t242;
t380 = t242 * t222 + t221;
t501 = -qJD(5) * t243 - t366 - t380;
t500 = m(7) * pkin(5);
t25 = qJD(6) * t337 + t138 * t247 + t139 * t244;
t448 = t25 / 0.2e1;
t26 = -qJD(6) * t101 - t138 * t244 + t139 * t247;
t447 = t26 / 0.2e1;
t499 = -m(4) - m(5);
t498 = -m(6) - m(7);
t175 = qJDD(6) + t181;
t434 = t175 / 0.2e1;
t170 = t182 * pkin(7);
t497 = qJD(2) / 0.2e1;
t495 = mrSges(2,2) - mrSges(3,3);
t494 = mrSges(7,3) - mrSges(6,2);
t356 = -pkin(7) * t242 - pkin(3);
t292 = (-qJ(5) + t356) * t245;
t265 = pkin(8) * t387 + t292;
t353 = t243 * t372;
t300 = pkin(2) * t245 - qJ(3) * t248;
t180 = t300 * qJD(1);
t391 = t180 * t243;
t336 = pkin(4) * t353 - t391;
t54 = qJD(1) * t265 + t336;
t358 = t242 * t441;
t263 = (-pkin(7) * t243 + pkin(5)) * t245 - t248 * t358;
t159 = t242 * t180;
t215 = qJ(4) * t373;
t381 = t159 + t215;
t55 = qJD(1) * t263 + t381;
t226 = t242 * qJ(3);
t187 = t242 * pkin(4) + t226;
t165 = pkin(8) * t242 + t187;
t227 = t243 * qJ(3);
t188 = t243 * pkin(4) + t227;
t166 = pkin(8) * t243 + t188;
t95 = t165 * t247 + t166 * t244;
t490 = -qJD(3) * t177 - qJD(6) * t95 + t244 * t54 - t247 * t55;
t94 = -t165 * t244 + t166 * t247;
t489 = qJD(3) * t176 + qJD(6) * t94 - t244 * t55 - t247 * t54;
t486 = m(6) * qJD(3);
t321 = mrSges(3,1) * t248 - mrSges(3,2) * t245;
t484 = -mrSges(2,1) - t321;
t399 = qJ(5) * t242;
t287 = -t243 * t416 + t399;
t268 = t248 * t287;
t482 = -qJD(1) * t268 + t501;
t480 = t243 * t418;
t284 = t176 * t248;
t137 = qJD(1) * t284;
t152 = t176 * qJD(6);
t478 = t137 - t152;
t296 = qJ(4) * t243 - t399;
t285 = t296 * t248;
t477 = -qJD(1) * t285 - t501;
t476 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t171 + mrSges(4,2) * t172 + mrSges(3,3) * t373;
t233 = t248 * qJ(5);
t388 = t243 * t245;
t475 = pkin(4) * t388 + t233;
t306 = Ifges(6,4) * t242 + Ifges(6,6) * t243;
t473 = t245 * (Ifges(3,1) * t248 - t415) + t248 * (-Ifges(6,2) * t245 + t248 * t306);
t371 = qJD(2) * t245;
t472 = -qJ(4) * t371 + t365;
t169 = -pkin(7) * t349 + t218;
t471 = t169 * t248 + t170 * t245;
t470 = (mrSges(5,1) + mrSges(4,3)) * t245;
t469 = -t242 * t37 + t243 * t38;
t246 = sin(qJ(1));
t249 = cos(qJ(1));
t468 = g(1) * t249 + g(2) * t246;
t466 = mrSges(6,1) + mrSges(5,3) - mrSges(4,2);
t465 = mrSges(6,3) + mrSges(4,1) - mrSges(5,2);
t130 = -mrSges(6,2) * t172 + mrSges(6,3) * t372;
t132 = mrSges(5,1) * t172 - mrSges(5,2) * t372;
t135 = -mrSges(4,1) * t372 - mrSges(4,3) * t172;
t462 = -t135 + t132 + t130;
t370 = qJD(2) * t248;
t223 = pkin(7) * t370;
t379 = t242 * pkin(3) * t370 + t223;
t459 = -t245 * (qJD(4) * t243 - qJD(5) * t242) + t379;
t141 = -t242 * t416 - pkin(2) - t480;
t340 = qJ(4) * t242 + pkin(2);
t183 = -pkin(3) * t243 - t340;
t405 = t242 * mrSges(6,1);
t458 = t521 * t248 + (-m(7) * t141 + m(6) * t340 + t405 - (-m(6) * t418 - mrSges(6,3)) * t243 - m(5) * t183 + m(4) * pkin(2) - t507) * t245;
t457 = (-t245 * t517 + t248 * t511) * t172 + (t245 * t519 + t248 * t512) * t171;
t10 = t139 * t441 + t264;
t11 = pkin(5) * t181 - t138 * t441 + t267;
t1 = qJD(6) * t8 + t10 * t247 + t11 * t244;
t2 = -qJD(6) * t9 - t10 * t244 + t11 * t247;
t456 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t450 = Ifges(7,4) * t448 + Ifges(7,2) * t447 + Ifges(7,6) * t434;
t449 = Ifges(7,1) * t448 + Ifges(7,4) * t447 + Ifges(7,5) * t434;
t96 = Ifges(7,4) * t337;
t30 = Ifges(7,1) * t101 + Ifges(7,5) * t207 + t96;
t445 = t30 / 0.2e1;
t443 = -t337 / 0.2e1;
t440 = -t101 / 0.2e1;
t431 = -t207 / 0.2e1;
t425 = pkin(8) * t245;
t234 = t245 * pkin(7);
t417 = -pkin(4) - qJ(3);
t414 = Ifges(3,4) * t248;
t406 = t139 * mrSges(6,1);
t109 = mrSges(6,1) * t172 - mrSges(6,3) * t171;
t35 = -mrSges(7,1) * t337 + mrSges(7,2) * t101;
t400 = t109 - t35;
t147 = -qJDD(2) * pkin(2) + qJDD(3) + t170;
t397 = t147 * t245;
t150 = qJD(2) * t300 - t367;
t396 = t150 * t243;
t383 = t249 * t242;
t157 = -t246 * t243 + t248 * t383;
t395 = t157 * qJ(4);
t390 = t242 * t245;
t386 = t245 * t249;
t385 = t246 * t248;
t384 = t248 * t249;
t378 = pkin(3) * t390 + t234;
t375 = t238 + t231;
t186 = -pkin(1) - t375;
t129 = pkin(7) * t387 + t242 * t186;
t374 = t249 * pkin(1) + t246 * pkin(7);
t361 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t175;
t360 = pkin(7) * t371;
t7 = -t26 * mrSges(7,1) + t25 * mrSges(7,2);
t339 = -t363 / 0.2e1;
t338 = t363 / 0.2e1;
t81 = t181 * mrSges(6,1) + t138 * mrSges(6,2);
t200 = pkin(7) * t389;
t128 = t186 * t243 - t200;
t335 = pkin(3) * t387 + qJ(4) * t389 + t375;
t334 = pkin(2) * t384 + qJ(3) * t386 + t374;
t329 = t356 * t245;
t155 = t242 * t385 + t243 * t249;
t156 = t243 * t385 - t383;
t239 = t249 * pkin(7);
t325 = -t156 * pkin(3) - qJ(4) * t155 + t239;
t116 = qJ(4) * t248 - t129;
t324 = t248 * pkin(4) * t364 + t229 - t396;
t237 = t248 * pkin(3);
t117 = -t128 + t237;
t140 = t242 * t150;
t323 = t140 - t472;
t320 = mrSges(3,1) * t245 + mrSges(3,2) * t248;
t316 = t243 * mrSges(6,3) + t405;
t145 = t177 * t245;
t146 = t176 * t245;
t315 = -mrSges(7,1) * t145 - mrSges(7,2) * t146;
t307 = -Ifges(5,4) * t243 + Ifges(5,5) * t242;
t305 = Ifges(3,5) * t248 - Ifges(3,6) * t245;
t304 = Ifges(4,5) * t243 - Ifges(4,6) * t242;
t62 = -mrSges(7,2) * t207 + mrSges(7,3) * t337;
t63 = mrSges(7,1) * t207 - mrSges(7,3) * t101;
t299 = -t244 * t63 + t247 * t62;
t298 = -t244 * t62 - t247 * t63;
t64 = t200 + t237 + (-t186 + t425) * t243 + t475;
t65 = -t245 * t358 + t129 - t350;
t21 = t244 * t65 + t247 * t64;
t20 = -t244 * t64 + t247 * t65;
t158 = t246 * t242 + t243 * t384;
t297 = t158 * pkin(3) + t334;
t295 = -t155 * t247 - t156 * t244;
t294 = t155 * t244 - t156 * t247;
t119 = -pkin(7) * t354 + t159;
t114 = -t243 * t360 + t140;
t290 = t245 * pkin(4) + t243 * t233 + t335;
t289 = pkin(1) * t320;
t286 = -pkin(4) * t389 - pkin(7) * t388;
t270 = -g(1) * t157 - g(2) * t155 - g(3) * t390;
t269 = -g(1) * t158 - g(2) * t156 - g(3) * t388;
t266 = pkin(4) * t386 + t158 * qJ(5) + t297;
t256 = t248 * (Ifges(5,1) * t245 + t248 * t307);
t254 = t248 * (Ifges(4,3) * t245 + t248 * t304);
t32 = t138 * pkin(3) - t139 * qJ(4) - t172 * qJD(4) + t147;
t14 = qJ(5) * t138 + qJD(5) * t171 + t32;
t205 = qJ(3) * t384;
t203 = qJ(3) * t385;
t191 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t372;
t174 = t181 * mrSges(6,3);
t173 = t181 * mrSges(5,2);
t154 = t340 + t480;
t143 = -qJ(4) * t388 + t378;
t126 = t139 * mrSges(5,3);
t125 = t139 * mrSges(4,2);
t122 = -qJ(4) * t353 + t380;
t118 = pkin(7) * t355 + t391;
t115 = t245 * t296 - t378;
t113 = t242 * t360 + t396;
t112 = qJD(1) * t329 - t391;
t111 = -mrSges(5,2) * t171 - mrSges(5,3) * t172;
t108 = -t119 - t215;
t107 = (-qJ(4) * t370 - qJD(4) * t245) * t243 + t379;
t103 = -pkin(4) * t390 - t116;
t98 = t245 * t287 + t378;
t97 = qJD(2) * t329 - t396;
t91 = t172 * Ifges(4,4) - t171 * Ifges(4,2) - Ifges(4,6) * t372;
t86 = -Ifges(5,4) * t372 - t172 * Ifges(5,2) + t171 * Ifges(5,6);
t83 = mrSges(4,1) * t181 - mrSges(4,3) * t139;
t82 = -mrSges(4,2) * t181 - mrSges(4,3) * t138;
t80 = -t139 * mrSges(6,2) - t174;
t79 = t139 * mrSges(5,1) + t173;
t78 = mrSges(5,1) * t138 - mrSges(5,3) * t181;
t77 = t117 + t475;
t75 = t157 * t247 + t158 * t244;
t74 = -t157 * t244 + t158 * t247;
t73 = qJD(1) * t286 + t381;
t70 = -t114 + t472;
t67 = qJD(2) * t284 - t153 * t245;
t66 = -qJD(2) * t283 - t152 * t245;
t61 = qJD(1) * t292 + t336;
t60 = -t138 * mrSges(5,2) - t126;
t59 = t138 * mrSges(4,1) + t125;
t58 = -t138 * mrSges(6,3) + t406;
t57 = qJD(2) * t285 - t459;
t56 = qJD(2) * t286 + t323;
t53 = qJD(2) * t292 + t324;
t48 = qJD(2) * t268 + t459;
t40 = qJD(2) * t263 + t323;
t39 = qJD(2) * t265 + t324;
t36 = -t172 * t416 - t326 + t481;
t19 = -mrSges(7,2) * t175 + mrSges(7,3) * t26;
t18 = mrSges(7,1) * t175 - mrSges(7,3) * t25;
t12 = -pkin(5) * t139 + t14;
t4 = -qJD(6) * t21 - t244 * t39 + t247 * t40;
t3 = qJD(6) * t20 + t244 * t40 + t247 * t39;
t5 = [m(4) * (t104 * t113 + t105 * t114 + t128 * t37 + t129 * t38 + (-t326 * t370 + t397) * pkin(7)) + t321 * t398 + (Ifges(3,4) * t182 + Ifges(6,4) * t138 + Ifges(6,6) * t139 + (-Ifges(3,2) - Ifges(6,2)) * t181) * t248 / 0.2e1 + (-m(4) * t334 - m(6) * (t266 + t395) - m(5) * (t297 + t395) - m(7) * t266 - t75 * mrSges(7,1) - t74 * mrSges(7,2) - m(3) * t374 + t484 * t249 + t495 * t246 - t465 * t158 + (-m(7) * t416 - t466) * t157 + t521 * t386) * g(2) + t457 * t497 + (t256 + t254) * t339 + (t505 + t526) * t371 + (Ifges(7,1) * t67 + Ifges(7,4) * t66) * t439 + t182 * t414 / 0.2e1 - t289 * t363 + (Ifges(7,5) * t67 + Ifges(7,6) * t66) * t430 + t318 * t397 + qJDD(2) * Ifges(3,6) * t248 + (Ifges(7,4) * t67 + Ifges(7,2) * t66) * t442 + (Ifges(7,1) * t146 - Ifges(7,4) * t145) * t448 + (Ifges(7,5) * t146 - Ifges(7,6) * t145) * t434 + (Ifges(7,4) * t146 - Ifges(7,2) * t145) * t447 + (-t1 * t145 - t146 * t2 + t66 * t9 - t67 * t8) * mrSges(7,3) - (t242 * t91 + t243 * t86) * t370 / 0.2e1 - t191 * t360 - t12 * t315 + Ifges(2,3) * qJDD(1) - pkin(1) * (mrSges(3,1) * t181 + mrSges(3,2) * t182) + t143 * t60 + t128 * t83 + t129 * t82 + t53 * t130 + t70 * t131 + t97 * t132 + t56 * t133 + t114 * t134 + t113 * t135 + t115 * t58 + t116 * t78 + t117 * t79 + t103 * t81 + t57 * t109 + t107 * t111 + t98 * t7 + t77 * t80 + m(5) * (t107 * t69 + t116 * t27 + t117 * t31 + t143 * t32 + t70 * t76 + t72 * t97) + m(6) * (t103 * t17 - t115 * t14 + t13 * t77 + t41 * t53 + t49 * t57 + t50 * t56) + m(7) * (t1 * t21 + t12 * t98 + t2 * t20 + t3 * t9 + t36 * t48 + t4 * t8) + t506 * t370 / 0.2e1 + (-t38 * mrSges(4,3) + t508) * t390 + (t305 * t497 + t509) * qJD(2) + (t182 * t234 + t471) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t471) + (-mrSges(4,3) * t37 + t510) * t388 + t473 * t338 + (t32 * t312 - t338 * t485 + Ifges(3,1) * t182 + Ifges(3,5) * qJDD(2) + t302 * t436 + t308 * t438 - t14 * t317 + t512 * t437 + t511 * t435 + (t304 + t307) * t432 + (Ifges(3,4) + t306) * t433) * t245 + t476 * t223 + (-qJDD(2) * mrSges(3,1) + t59) * t234 + t309 * t433 + t67 * t445 + t66 * t446 + t146 * t449 - t145 * t450 + (Ifges(6,2) * t433 - Ifges(5,4) * t436 - Ifges(4,6) * t438 + t414 * t338 - t37 * mrSges(4,1) + t13 * mrSges(6,3) - t31 * mrSges(5,2) + t38 * mrSges(4,2) - t17 * mrSges(6,1) + t27 * mrSges(5,3) + pkin(7) * (-qJDD(2) * mrSges(3,2) - mrSges(3,3) * t181) - Ifges(7,3) * t434 - Ifges(7,6) * t447 - Ifges(7,5) * t448 - t519 * t437 + t517 * t435 + t491 * t432 - t456) * t248 - (-t138 * t492 + t139 * t493 - t181 * t491 + t361) * t248 / 0.2e1 + t20 * t18 + t21 * t19 + (-m(5) * t325 - t295 * mrSges(7,1) - t294 * mrSges(7,2) + t498 * (-qJ(5) * t156 + t325) + t495 * t249 + (-m(4) - m(3)) * t239 + t465 * t156 + (t466 + t500) * t155 + (m(3) * pkin(1) + t498 * t357 + t499 * t291 + (-m(7) * (-pkin(8) + t417) - m(6) * t417 + t494) * t245 + t470 - t484) * t246) * g(1) + t48 * t35 + t3 * t62 + t4 * t63 + t36 * (-mrSges(7,1) * t66 + mrSges(7,2) * t67); (-pkin(2) * t147 + t469 * qJ(3) - t104 * t118 - t105 * t119 + t221 * t326) * m(4) + (Ifges(7,5) * t137 - Ifges(7,6) * t136) * t431 + (Ifges(7,1) * t137 - Ifges(7,4) * t136) * t440 + (Ifges(7,4) * t137 - Ifges(7,2) * t136) * t443 + (m(5) * (qJ(3) * t31 - qJD(4) * t69) + t41 * t486 + t91 * t345 + Ifges(6,6) * t433 - Ifges(5,2) * t436 + t522 + t493 * t432 - t500 * t503 + (-m(4) * t104 + m(5) * t72 + t462) * qJD(3) + t510) * t242 + (Ifges(7,5) * t440 + Ifges(7,6) * t443 + Ifges(7,3) * t431 - t505) * t373 + (-mrSges(7,1) * t513 - mrSges(7,2) * t478) * t36 + (t1 * t176 - t177 * t2 + t478 * t8 + t513 * t9) * mrSges(7,3) + (-t409 + t412) * t435 + (-t108 * t76 - t112 * t72 - t122 * t69 + t183 * t32) * m(5) + (-t366 - t122) * t111 + (-t408 + t410) * t437 - t407 * t436 + (t82 - t78) * t227 + (t79 - t83) * t226 + t413 * t438 + (-Ifges(3,2) * t373 + t219 + t506) * t523 + (Ifges(7,5) * t152 - Ifges(7,6) * t153) * t430 + (Ifges(7,1) * t152 - Ifges(7,4) * t153) * t439 + (Ifges(7,4) * t152 - Ifges(7,2) * t153) * t442 - t147 * t319 + t32 * t313 + t12 * t314 - t14 * t316 + Ifges(3,3) * qJDD(2) + Ifges(3,5) * t182 + t183 * t60 + t187 * t80 + t188 * t81 - Ifges(3,6) * t181 - t169 * mrSges(3,2) - t170 * mrSges(3,1) + t154 * t58 + t141 * t7 - t61 * t130 - t108 * t131 - t112 * t132 - t73 * t133 - t119 * t134 - t118 * t135 - t137 * t30 / 0.2e1 + t94 * t18 + t95 * t19 + (-m(7) * (t290 + t425) - m(4) * t375 - m(5) * t335 - m(6) * t290 - t321 - t494 * t245 + (-t316 + t507) * t248 - t470) * g(3) + (-t457 / 0.2e1 + (t254 / 0.2e1 + t256 / 0.2e1 + t289 - t473 / 0.2e1) * qJD(1) - t509) * qJD(1) + t191 * t220 + t468 * t320 + t469 * mrSges(4,3) - t476 * t221 + t477 * t109 + (t13 * t187 - t14 * t154 + t17 * t188 - t41 * t61 + t477 * t49 - t50 * t73) * m(6) + t513 * t446 + t482 * t35 + t305 * t339 + (Ifges(7,5) * t177 + Ifges(7,6) * t176) * t434 + t152 * t445 + (Ifges(7,4) * t177 + Ifges(7,2) * t176) * t447 + (Ifges(7,1) * t177 + Ifges(7,4) * t176) * t448 + t177 * t449 + t176 * t450 + (-m(5) * qJ(3) * t27 + t50 * t486 + t86 * t345 - Ifges(6,4) * t433 + Ifges(4,2) * t438 + t516 * t437 + t492 * t432 + (m(4) * t105 - m(5) * t76 + t502) * qJD(3) - t508) * t243 + t489 * t62 + t490 * t63 + (t1 * t95 + t12 * t141 + t2 * t94 + t36 * t482 + t489 * t9 + t490 * t8) * m(7) + (t498 * (pkin(4) * t384 + t205) + t499 * t205 + t458 * t249) * g(1) + (t498 * (pkin(4) * t385 + t203) + t499 * t203 + t458 * t246) * g(2) - pkin(2) * t59; -t406 + t101 * t63 - t337 * t62 + t125 - t126 - t462 * t172 + t502 * t171 + t465 * t138 + t7 + (t101 * t8 - t337 * t9 + t12) * m(7) + (t171 * t50 - t172 * t41 + t14) * m(6) + (-t171 * t76 - t172 * t72 + t32) * m(5) + (t104 * t172 + t105 * t171 + t147) * m(4) + (-t245 * t468 + t503) * (-t498 - t499); -t244 * t18 + t247 * t19 + t173 - t174 + t496 * t139 + t298 * qJD(6) + (t111 - t400) * t172 + (-t298 + t382) * t372 + (t1 * t247 + t172 * t36 - t2 * t244 + t270 + t207 * (-t244 * t9 - t247 * t8)) * m(7) + (-t172 * t49 + t372 * t50 + t13 + t270) * m(6) + (t172 * t69 - t372 * t76 + t270 + t31) * m(5); t247 * t18 + t244 * t19 + t400 * t171 + t299 * qJD(6) + (-t130 - t299) * t372 + t81 + (t1 * t244 - t171 * t36 + t2 * t247 + t269 + t207 * (-t244 * t8 + t247 * t9)) * m(7) + (t171 * t49 - t372 * t41 + t17 + t269) * m(6); -t36 * (mrSges(7,1) * t101 + mrSges(7,2) * t337) + (Ifges(7,1) * t337 - t411) * t440 + t29 * t439 + (Ifges(7,5) * t337 - Ifges(7,6) * t101) * t431 - t8 * t62 + t9 * t63 - g(1) * (mrSges(7,1) * t74 - mrSges(7,2) * t75) - g(2) * (-mrSges(7,1) * t294 + mrSges(7,2) * t295) - g(3) * t315 + (t101 * t9 + t337 * t8) * mrSges(7,3) + t361 + (-Ifges(7,2) * t101 + t30 + t96) * t443 + t456;];
tau  = t5;
