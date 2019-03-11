% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:51
% EndTime: 2019-03-09 02:47:02
% DurationCPUTime: 5.74s
% Computational Cost: add. (13101->467), mult. (26698->645), div. (0->0), fcn. (29936->8), ass. (0->240)
t304 = sin(pkin(10));
t306 = cos(pkin(10));
t308 = sin(qJ(6));
t310 = cos(qJ(6));
t287 = t304 * t310 - t308 * t306;
t305 = sin(pkin(9));
t307 = cos(pkin(9));
t309 = sin(qJ(3));
t427 = cos(qJ(3));
t288 = t305 * t427 + t309 * t307;
t454 = t287 * t288;
t411 = t454 * mrSges(7,2);
t463 = t411 / 0.2e1;
t337 = t308 * t304 + t310 * t306;
t327 = t337 * t288;
t408 = t327 * mrSges(7,1);
t285 = t305 * t309 - t307 * t427;
t378 = t288 * t304;
t217 = -t285 * mrSges(5,2) - mrSges(5,3) * t378;
t222 = -mrSges(6,2) * t378 + t285 * mrSges(6,3);
t453 = t217 + t222;
t462 = t453 * t306;
t300 = t304 ^ 2;
t302 = t306 ^ 2;
t367 = t300 + t302;
t456 = mrSges(5,3) + mrSges(6,2);
t461 = t408 / 0.2e1;
t438 = t337 / 0.2e1;
t460 = t438 * t454;
t231 = t287 * mrSges(7,1) - mrSges(7,2) * t337;
t459 = t231 * qJD(6);
t420 = Ifges(5,6) - Ifges(6,6);
t455 = Ifges(5,5) + Ifges(6,4);
t458 = t304 * t420 - t306 * t455 + Ifges(4,4);
t190 = Ifges(7,4) * t454;
t91 = Ifges(7,1) * t327 - t285 * Ifges(7,5) + t190;
t457 = t91 / 0.2e1;
t381 = t285 * t304;
t196 = t337 * t285;
t409 = t196 * mrSges(7,2);
t193 = t287 * t285;
t412 = t193 * mrSges(7,1);
t372 = -t412 / 0.2e1 + t409 / 0.2e1;
t403 = t285 * mrSges(7,2);
t410 = t454 * mrSges(7,3);
t154 = t403 + t410;
t404 = t285 * mrSges(7,1);
t407 = t327 * mrSges(7,3);
t156 = -t404 - t407;
t435 = -t287 / 0.2e1;
t437 = -t337 / 0.2e1;
t331 = t154 * t437 + t156 * t435;
t319 = (t327 * t435 + t460) * mrSges(7,3) + t331;
t12 = t319 - t372;
t451 = t12 * qJD(1);
t269 = t310 * t287;
t382 = t337 * t308;
t446 = m(7) / 0.2e1;
t326 = (t269 + t382) * t446;
t450 = m(5) / 0.2e1;
t449 = -m(6) / 0.2e1;
t448 = m(6) / 0.2e1;
t447 = -m(7) / 0.2e1;
t406 = t327 * Ifges(7,4);
t90 = Ifges(7,2) * t454 - t285 * Ifges(7,6) + t406;
t445 = -t90 / 0.2e1;
t444 = pkin(4) + pkin(5);
t442 = -t193 / 0.2e1;
t282 = Ifges(7,4) * t337;
t237 = Ifges(7,1) * t287 - t282;
t440 = t237 / 0.2e1;
t418 = -pkin(8) + qJ(4);
t291 = t418 * t304;
t293 = t418 * t306;
t241 = t308 * t291 + t293 * t310;
t439 = -t241 / 0.2e1;
t436 = -t285 / 0.2e1;
t434 = t287 / 0.2e1;
t433 = t288 / 0.2e1;
t432 = -t304 / 0.2e1;
t431 = t304 / 0.2e1;
t429 = -t308 / 0.2e1;
t428 = -t310 / 0.2e1;
t419 = pkin(7) + qJ(2);
t294 = t419 * t307;
t349 = t419 * t305;
t240 = t294 * t309 + t427 * t349;
t227 = t304 * t240;
t392 = qJ(4) * t285;
t424 = pkin(3) * t288;
t232 = t392 + t424;
t139 = t232 * t306 + t227;
t109 = -pkin(4) * t288 - t139;
t426 = m(6) * t109;
t425 = m(6) * t304;
t423 = Ifges(5,1) + Ifges(6,1);
t422 = Ifges(5,4) - Ifges(6,5);
t421 = Ifges(5,2) + Ifges(6,3);
t415 = Ifges(7,4) * t287;
t279 = t288 * mrSges(4,1);
t140 = t304 * t232 - t306 * t240;
t106 = t288 * qJ(5) + t140;
t153 = mrSges(7,2) * t288 - t193 * mrSges(7,3);
t155 = -mrSges(7,1) * t288 + t196 * mrSges(7,3);
t216 = -t288 * mrSges(5,2) + mrSges(5,3) * t381;
t380 = t285 * t306;
t218 = t288 * mrSges(5,1) + mrSges(5,3) * t380;
t363 = mrSges(6,2) * t380;
t401 = t288 * mrSges(6,1);
t219 = -t363 - t401;
t223 = mrSges(6,2) * t381 + mrSges(6,3) * t288;
t75 = -t227 + (pkin(8) * t285 - t232) * t306 - t444 * t288;
t87 = -pkin(8) * t381 + t106;
t40 = -t308 * t87 + t310 * t75;
t41 = t308 * t75 + t310 * t87;
t312 = (t218 / 0.2e1 - t219 / 0.2e1) * t306 + (t223 / 0.2e1 + t216 / 0.2e1) * t304 + (t139 * t306 + t140 * t304) * t450 + (t106 * t304 - t109 * t306) * t448 + (t287 * t41 - t337 * t40) * t446 + t155 * t437 + t153 * t434;
t239 = t291 * t310 - t308 * t293;
t348 = qJ(5) * t304 + pkin(3);
t276 = t306 * t444 + t348;
t346 = t367 * t392;
t289 = -pkin(4) * t306 - t348;
t379 = t288 * t289;
t315 = -m(5) * (-t346 - t424) / 0.2e1 + (-t346 + t379) * t449 + (-t193 * t239 - t196 * t241 - t276 * t288) * t447;
t233 = mrSges(7,1) * t337 + mrSges(7,2) * t287;
t292 = -mrSges(6,1) * t306 - mrSges(6,3) * t304;
t352 = t233 / 0.2e1 - t292 / 0.2e1;
t399 = t306 * mrSges(5,1);
t400 = t304 * mrSges(5,2);
t8 = t279 + (t399 / 0.2e1 - t400 / 0.2e1 + t352) * t288 + (-t193 * t434 - t196 * t438) * mrSges(7,3) + (-mrSges(4,2) + t456 * (t300 / 0.2e1 + t302 / 0.2e1)) * t285 + t312 + t315;
t414 = qJD(1) * t8;
t103 = -mrSges(7,1) * t454 + mrSges(7,2) * t327;
t360 = t444 * t304;
t390 = qJ(5) * t306;
t107 = (-t360 + t390) * t288 - t240;
t242 = t294 * t427 - t309 * t349;
t321 = qJ(5) * t380 + t242;
t108 = -t285 * t360 + t321;
t355 = -pkin(2) * t307 - pkin(1);
t212 = pkin(3) * t285 - qJ(4) * t288 + t355;
t228 = t304 * t242;
t120 = t212 * t306 - t228;
t121 = t304 * t212 + t306 * t242;
t144 = (pkin(4) * t304 - t390) * t288 + t240;
t145 = -pkin(4) * t381 + t321;
t343 = mrSges(6,1) * t304 - mrSges(6,3) * t306;
t211 = t343 * t288;
t377 = t306 * t288;
t220 = t285 * mrSges(5,1) - mrSges(5,3) * t377;
t221 = -t285 * mrSges(6,1) + mrSges(6,2) * t377;
t332 = -Ifges(7,5) * t196 / 0.2e1 + Ifges(7,6) * t442;
t340 = -Ifges(7,4) * t196 - Ifges(7,2) * t193;
t341 = -Ifges(7,1) * t196 - Ifges(7,4) * t193;
t342 = -t409 + t412;
t344 = mrSges(5,1) * t304 + mrSges(5,2) * t306;
t350 = t422 * t306;
t65 = t228 + (-pkin(8) * t288 - t212) * t306 - t444 * t285;
t95 = t285 * qJ(5) + t121;
t80 = pkin(8) * t378 + t95;
t36 = -t308 * t80 + t310 * t65;
t37 = t308 * t65 + t310 * t80;
t98 = -pkin(4) * t285 - t120;
t1 = -m(5) * (t120 * t139 + t121 * t140 + t240 * t242) - m(6) * (t106 * t95 + t109 * t98 + t144 * t145) - m(7) * (-t107 * t108 + t36 * t40 + t37 * t41) - t355 * t279 - t327 * t341 / 0.2e1 - t454 * t340 / 0.2e1 - t193 * t445 + (t355 * mrSges(4,2) + t144 * t343 + t240 * t344 - t458 * t285 + t332) * t285 + ((t423 * t302 + Ifges(4,1) - Ifges(4,2) - Ifges(6,2) - Ifges(5,3) - Ifges(7,3) + (t304 * t421 - 0.2e1 * t350) * t304) * t285 - t344 * t242 + Ifges(7,5) * t327 + Ifges(7,6) * t454 + t458 * t288) * t288 + t196 * t457 - t107 * t342 + t108 * t103 - t37 * t153 - t41 * t154 - t36 * t155 - t40 * t156 - t145 * t211 - t121 * t216 - t140 * t217 - t120 * t218 - t98 * t219 - t139 * t220 - t109 * t221 - t106 * t222 - t95 * t223;
t413 = t1 * qJD(1);
t405 = t337 * mrSges(7,3);
t402 = t287 * mrSges(7,3);
t398 = t306 * mrSges(6,2);
t397 = t308 * mrSges(7,1);
t396 = t310 * mrSges(7,2);
t102 = t408 + t411;
t104 = -Ifges(7,2) * t327 + t190;
t105 = Ifges(7,1) * t454 - t406;
t371 = Ifges(7,5) * t454 - Ifges(7,6) * t327;
t4 = -t37 * t156 + t36 * t154 + t107 * t102 + t371 * t436 + (-t37 * mrSges(7,3) + t105 / 0.2e1 + t445) * t327 + (-t36 * mrSges(7,3) + t457 + t104 / 0.2e1) * t454;
t395 = t4 * qJD(1);
t338 = t120 * t304 - t121 * t306;
t339 = t304 * t98 + t306 * t95;
t369 = t220 - t221;
t373 = t103 - t211;
t384 = t240 * t288;
t5 = -t196 * t154 - t193 * t156 + (mrSges(4,3) * t285 + t304 * t369 - t462) * t285 + ((mrSges(4,3) + t344) * t288 - t373) * t288 + m(7) * (-t107 * t288 - t193 * t36 - t196 * t37) + m(6) * (t144 * t288 - t285 * t339) + m(5) * (t285 * t338 + t384) + m(4) * (-t242 * t285 + t384) + (m(3) * qJ(2) + mrSges(3,3)) * (t305 ^ 2 + t307 ^ 2);
t394 = t5 * qJD(1);
t391 = qJ(4) * t306;
t10 = m(7) * (t327 * t36 - t37 * t454) + t327 * t156 - t454 * t154 + (-t369 * t306 - t453 * t304 + m(6) * (-t304 * t95 + t306 * t98) + m(5) * (-t120 * t306 - t121 * t304)) * t288;
t389 = qJD(1) * t10;
t374 = t310 * t154;
t376 = t308 * t156;
t15 = t373 * t377 + (t222 + t374 - t376) * t285 + m(7) * (t107 * t377 + (-t308 * t36 + t310 * t37) * t285) + m(6) * (-t144 * t377 + t285 * t95);
t388 = qJD(1) * t15;
t330 = m(7) * (-t193 * t310 - t308 * t196);
t44 = -t330 / 0.2e1 + (t326 + t425) * t285;
t387 = qJD(1) * t44;
t385 = t327 * t310;
t375 = t308 * t454;
t368 = -Ifges(7,5) * t337 - Ifges(7,6) * t287;
t364 = m(6) / 0.4e1 + m(5) / 0.4e1;
t362 = m(6) * t433;
t357 = -t405 / 0.2e1;
t356 = -t402 / 0.2e1;
t345 = t352 * t306;
t101 = (-m(6) * t289 + m(7) * t276 + t233 - t292) * t304;
t314 = (t103 / 0.2e1 - t211 / 0.2e1) * t304 + (-t304 * t144 + (-t379 + t392) * t306) * t448 + (t276 * t377 + t304 * t107 + (-t239 * t308 + t241 * t310) * t285) * t446;
t317 = -t426 / 0.2e1 + (t308 * t41 + t310 * t40) * t447 + t153 * t429 + t155 * t428;
t325 = (t308 * t434 + t310 * t437) * mrSges(7,3);
t14 = (mrSges(6,1) / 0.2e1 + t345) * t288 + (t325 + t398) * t285 + t314 + t317;
t336 = -qJD(1) * t14 - qJD(3) * t101;
t50 = 0.2e1 * t461 + 0.2e1 * t463;
t335 = qJD(1) * t50 + qJD(3) * t231;
t167 = t425 + (t382 / 0.2e1 + t269 / 0.2e1 + t431) * m(7);
t57 = -m(6) * t377 + (-t377 / 0.2e1 - t385 / 0.2e1 + t375 / 0.2e1) * m(7);
t334 = qJD(1) * t57 - qJD(3) * t167;
t19 = (-t403 / 0.2e1 + t410 / 0.2e1 - t154 / 0.2e1) * t310 + (-t404 / 0.2e1 + t156 / 0.2e1 + t407 / 0.2e1) * t308;
t329 = t19 * qJD(1);
t320 = (-t287 * t454 - t327 * t337) * t446 - 0.2e1 * t364 * t367 * t288;
t35 = 0.2e1 * (-m(7) / 0.4e1 - t364) * t288 + t320;
t328 = t35 * qJD(1);
t234 = -Ifges(7,2) * t287 - t282;
t235 = -Ifges(7,2) * t337 + t415;
t236 = -Ifges(7,1) * t337 - t415;
t26 = t276 * t231 + (-t235 / 0.2e1 + t236 / 0.2e1) * t287 - (t234 / 0.2e1 + t440) * t337;
t311 = (t105 / 0.4e1 - t90 / 0.4e1) * t287 - (t91 / 0.4e1 + t104 / 0.4e1) * t337 + (-t239 * mrSges(7,3) / 0.2e1 + t237 / 0.4e1 + t234 / 0.4e1) * t454 + (mrSges(7,3) * t439 + t236 / 0.4e1 - t235 / 0.4e1) * t327 + t107 * t231 / 0.2e1 + t239 * t154 / 0.2e1 + t156 * t439 + t276 * t102 / 0.2e1 - t285 * t368 / 0.4e1;
t318 = Ifges(7,3) * t433 - t40 * mrSges(7,1) / 0.2e1 + t41 * mrSges(7,2) / 0.2e1 - t332;
t3 = t311 + t318;
t324 = -t3 * qJD(1) - t26 * qJD(3);
t46 = (-t287 ^ 2 - t337 ^ 2) * mrSges(7,3) + m(7) * (t239 * t287 + t241 * t337) + (0.4e1 * qJ(4) * t364 + t456) * t367;
t313 = t338 * t450 + t339 * t449 + (t239 * t327 - t241 * t454 + t287 * t36 + t337 * t37) * t447 + t331;
t316 = ((mrSges(6,3) / 0.2e1 - mrSges(5,2) / 0.2e1) * t306 + (-mrSges(6,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t304) * t285 + t242 * t450 + t145 * t448 + t108 * t446 + t372;
t6 = (-t222 / 0.2e1 - t217 / 0.2e1) * t306 + (t220 / 0.2e1 - t221 / 0.2e1) * t304 + (t327 * t434 - t460) * mrSges(7,3) + t313 + t316;
t323 = -t6 * qJD(1) + t46 * qJD(3);
t176 = m(7) * t432 + t326;
t56 = (-t375 + t385) * t446 + t306 * t362 + (t447 + t449) * t377;
t51 = -t408 / 0.2e1 - t411 / 0.2e1 + t461 + t463;
t43 = t381 * t449 + t330 / 0.2e1 + (t326 + t425 / 0.2e1) * t285;
t34 = t320 + t362 + (m(5) + m(7)) * t433;
t20 = -t376 / 0.2e1 + t410 * t428 + t374 / 0.2e1 + t407 * t429 + (-t397 / 0.2e1 - t396 / 0.2e1) * t285;
t13 = -t363 / 0.2e1 - t401 / 0.2e1 + t288 * t345 + (t398 / 0.2e1 + t325) * t285 + t314 - t317;
t11 = t319 + t372;
t9 = t312 - t193 * t356 - t196 * t357 - t288 * t233 / 0.2e1 - t315 + (-t399 + t400 + t292) * t433 + t456 * t367 * t436;
t7 = t327 * t356 - t454 * t357 + t220 * t432 + t221 * t431 - t313 + t316 + t462 / 0.2e1;
t2 = t311 - t318;
t16 = [qJD(2) * t5 - qJD(3) * t1 + qJD(4) * t10 + qJD(5) * t15 + qJD(6) * t4, t394 + (t193 * t337 - t196 * t287) * m(7) * qJD(2) + t9 * qJD(3) + t34 * qJD(4) + t43 * qJD(5) + t11 * qJD(6), t9 * qJD(2) + t7 * qJD(4) + t13 * qJD(5) + t2 * qJD(6) - t413 + ((-t242 * mrSges(5,1) + t106 * mrSges(6,2) + t140 * mrSges(5,3) + (t216 + t223) * qJ(4) + (t421 - t423) * t381) * t306 - t242 * mrSges(4,1) + t240 * mrSges(4,2) - t108 * t233 + t145 * t292 + t241 * t153 + t239 * t155 - t196 * t440 + t235 * t442 + t276 * t342 + t340 * t437 + t341 * t434 - t40 * t402 - t41 * t405 + 0.2e1 * (-t108 * t276 + t239 * t40 + t241 * t41) * t446 + 0.2e1 * (-pkin(3) * t242 + t140 * t391) * t450 + 0.2e1 * (t106 * t391 + t145 * t289) * t448 + (-Ifges(4,5) + (mrSges(5,2) * pkin(3) + mrSges(6,3) * t289 - t350) * t306) * t285 + (-Ifges(7,5) * t287 + Ifges(7,6) * t337 + t420 * t306 - Ifges(4,6)) * t288 + (t242 * mrSges(5,2) + t109 * mrSges(6,2) - t139 * mrSges(5,3) + t455 * t288 + (-m(5) * t139 - t218 + t219 + t426) * qJ(4) + (pkin(3) * mrSges(5,1) - t289 * mrSges(6,1) + t304 * t422) * t285) * t304) * qJD(3), qJD(2) * t34 + qJD(3) * t7 + qJD(5) * t56 + qJD(6) * t51 + t389, qJD(2) * t43 + qJD(3) * t13 + qJD(4) * t56 + qJD(6) * t20 + t388, t395 + t11 * qJD(2) + t2 * qJD(3) + t51 * qJD(4) + t20 * qJD(5) + (-mrSges(7,1) * t37 - mrSges(7,2) * t36 + t371) * qJD(6); qJD(3) * t8 + qJD(4) * t35 + qJD(5) * t44 + qJD(6) * t12 - t394, 0, t414, t328, t387, t451 - t459; -qJD(2) * t8 - qJD(4) * t6 + qJD(5) * t14 + qJD(6) * t3 + t413, -t414, qJD(4) * t46 + qJD(5) * t101 + qJD(6) * t26, t176 * qJD(5) + t323, qJD(4) * t176 - t336 (-mrSges(7,1) * t241 - mrSges(7,2) * t239 + t368) * qJD(6) - t324; -qJD(2) * t35 + qJD(3) * t6 + qJD(5) * t57 - qJD(6) * t50 - t389, -t328, -t167 * qJD(5) - t323 - t459, 0, t334, -t335; -qJD(2) * t44 - qJD(3) * t14 - qJD(4) * t57 - qJD(6) * t19 - t388, -t387, qJD(4) * t167 + t336, -t334, 0 (-t396 - t397) * qJD(6) - t329; -qJD(2) * t12 - qJD(3) * t3 + qJD(4) * t50 + qJD(5) * t19 - t395, -t451, qJD(4) * t231 + t324, t335, t329, 0;];
Cq  = t16;
