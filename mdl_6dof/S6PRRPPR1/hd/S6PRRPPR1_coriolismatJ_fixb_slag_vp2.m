% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:40
% EndTime: 2019-03-08 20:58:54
% DurationCPUTime: 7.26s
% Computational Cost: add. (16219->492), mult. (36337->718), div. (0->0), fcn. (42327->12), ass. (0->266)
t308 = cos(qJ(3));
t306 = sin(qJ(3));
t396 = cos(pkin(11));
t353 = t396 * t306;
t395 = sin(pkin(11));
t273 = -t308 * t395 - t353;
t302 = sin(pkin(12));
t304 = cos(pkin(12));
t427 = sin(qJ(6));
t428 = cos(qJ(6));
t469 = t427 * t302 - t428 * t304;
t484 = t469 * t273;
t487 = -t484 / 0.2e1;
t488 = mrSges(7,1) * t487;
t327 = t302 * t428 + t304 * t427;
t431 = -t327 / 0.2e1;
t486 = t431 * t484;
t352 = t395 * t306;
t271 = -t308 * t396 + t352;
t386 = t271 * t302;
t423 = -qJ(4) - pkin(8);
t282 = t423 * t308;
t471 = -t396 * t282 + t423 * t352;
t173 = -pkin(5) * t386 + t471;
t459 = -m(7) / 0.2e1;
t461 = -m(6) / 0.2e1;
t197 = t469 * t271;
t411 = t197 * mrSges(7,2);
t194 = t327 * t271;
t413 = t194 * mrSges(7,1);
t474 = t413 / 0.2e1 - t411 / 0.2e1;
t485 = t173 * t459 + t461 * t471 + t474;
t458 = m(7) / 0.2e1;
t401 = t304 * mrSges(6,2);
t404 = t302 * mrSges(6,1);
t347 = -t401 - t404;
t217 = t347 * t271;
t384 = t273 * t302;
t223 = -t271 * mrSges(6,2) + mrSges(6,3) * t384;
t377 = t304 * t223;
t383 = t273 * t304;
t225 = t271 * mrSges(6,1) + mrSges(6,3) * t383;
t381 = t302 * t225;
t332 = t377 / 0.2e1 - t381 / 0.2e1;
t483 = t217 / 0.2e1 - t332;
t360 = t395 * pkin(3);
t285 = t360 + qJ(5);
t299 = t304 ^ 2;
t372 = t302 ^ 2 + t299;
t351 = t372 * t285;
t482 = -m(6) * t351 - t372 * mrSges(6,3);
t293 = -pkin(3) * t308 - pkin(2);
t481 = m(5) * t293;
t324 = t327 * t273;
t477 = -t324 / 0.2e1;
t480 = mrSges(7,2) * t477;
t305 = cos(pkin(6));
t303 = sin(pkin(6));
t307 = sin(qJ(2));
t380 = t303 * t307;
t261 = t305 * t308 - t306 * t380;
t379 = t303 * t308;
t262 = t305 * t306 + t307 * t379;
t200 = -t396 * t261 + t262 * t395;
t325 = t261 * t395 + t262 * t396;
t123 = t200 * t325;
t369 = m(7) / 0.4e1 + m(6) / 0.4e1;
t478 = -0.2e1 * t369;
t309 = cos(qJ(2));
t476 = m(4) * t309;
t230 = mrSges(7,1) * t469 + mrSges(7,2) * t327;
t281 = -mrSges(6,1) * t304 + mrSges(6,2) * t302;
t472 = t230 + t281;
t470 = t306 ^ 2 + t308 ^ 2;
t426 = pkin(3) * t306;
t221 = -pkin(4) * t273 + qJ(5) * t271 + t426;
t242 = -t282 * t395 - t423 * t353;
t141 = t304 * t221 + t242 * t302;
t142 = t302 * t221 - t242 * t304;
t468 = -t141 * t302 + t142 * t304;
t283 = t306 * mrSges(4,1) + t308 * mrSges(4,2);
t467 = -mrSges(4,1) * t308 + mrSges(4,2) * t306;
t457 = m(5) * pkin(3);
t465 = t395 * t457 - mrSges(5,2);
t105 = t327 * t200;
t106 = t469 * t200;
t116 = -mrSges(7,1) * t324 + mrSges(7,2) * t484;
t146 = mrSges(7,2) * t273 + t194 * mrSges(7,3);
t148 = -mrSges(7,1) * t273 - t197 * mrSges(7,3);
t149 = mrSges(7,1) * t271 - mrSges(7,3) * t484;
t378 = t303 * t309;
t163 = -t302 * t325 - t304 * t378;
t164 = -t302 * t378 + t304 * t325;
t174 = -pkin(5) * t384 + t242;
t218 = t347 * t273;
t222 = mrSges(6,2) * t273 + mrSges(6,3) * t386;
t385 = t271 * t304;
t224 = -mrSges(6,1) * t273 + mrSges(6,3) * t385;
t220 = pkin(4) * t271 + qJ(5) * t273 + t293;
t139 = t304 * t220 - t302 * t471;
t140 = t302 * t220 + t304 * t471;
t342 = t139 * t302 - t140 * t304;
t147 = -mrSges(7,2) * t271 + mrSges(7,3) * t324;
t452 = t147 / 0.2e1;
t454 = t105 / 0.2e1;
t460 = m(6) / 0.2e1;
t462 = m(5) / 0.2e1;
t110 = pkin(9) * t384 + t140;
t90 = pkin(5) * t271 + pkin(9) * t383 + t139;
t56 = -t110 * t427 + t428 * t90;
t57 = t110 * t428 + t427 * t90;
t111 = pkin(9) * t386 + t142;
t91 = -pkin(5) * t273 + pkin(9) * t385 + t141;
t58 = -t111 * t427 + t428 * t91;
t59 = t111 * t428 + t427 * t91;
t84 = t163 * t428 - t164 * t427;
t85 = t163 * t427 + t164 * t428;
t464 = -t378 * t426 * t462 + (t141 * t163 + t142 * t164 + t242 * t325 + (t342 + t471) * t200) * t460 + (t105 * t56 + t106 * t57 + t173 * t200 + t174 * t325 + t58 * t84 + t59 * t85) * t458 + t149 * t454 + t106 * t452 + t163 * t224 / 0.2e1 + t164 * t222 / 0.2e1 + t84 * t148 / 0.2e1 + t85 * t146 / 0.2e1 + (t116 / 0.2e1 + t218 / 0.2e1) * t325;
t362 = t396 * pkin(3);
t292 = -t362 - pkin(4);
t463 = m(6) * t292 - t396 * t457 - mrSges(5,1) + t281;
t456 = -mrSges(7,2) / 0.2e1;
t455 = -mrSges(6,3) / 0.2e1;
t450 = t194 / 0.2e1;
t447 = t197 / 0.2e1;
t444 = t200 / 0.2e1;
t424 = pkin(9) + t285;
t263 = t424 * t302;
t264 = t424 * t304;
t216 = -t263 * t427 + t264 * t428;
t442 = -t216 / 0.2e1;
t228 = mrSges(7,1) * t327 - mrSges(7,2) * t469;
t439 = t228 / 0.2e1;
t414 = Ifges(7,4) * t327;
t232 = -Ifges(7,2) * t469 + t414;
t438 = t232 / 0.2e1;
t269 = Ifges(7,4) * t469;
t234 = Ifges(7,1) * t327 - t269;
t437 = t234 / 0.2e1;
t435 = -t273 / 0.2e1;
t434 = t273 / 0.2e1;
t433 = -t469 / 0.2e1;
t432 = t327 / 0.2e1;
t430 = t302 / 0.2e1;
t429 = t304 / 0.2e1;
t48 = -t105 * t469 + t106 * t327;
t422 = t48 * qJD(3) * t458;
t420 = m(7) * qJD(4);
t417 = Ifges(6,4) * t302;
t416 = Ifges(6,4) * t304;
t415 = Ifges(7,4) * t484;
t412 = t324 * mrSges(7,2);
t410 = t484 * mrSges(7,1);
t408 = t271 * mrSges(5,3);
t407 = t273 * mrSges(5,3);
t406 = t469 * mrSges(7,3);
t405 = t327 * mrSges(7,3);
t403 = t302 * Ifges(6,2);
t402 = t302 * Ifges(6,6);
t400 = t304 * Ifges(6,5);
t246 = t273 * t378;
t143 = t200 * t246;
t247 = t271 * t378;
t213 = t247 * t302 + t304 * t380;
t392 = t213 * t302;
t215 = -t247 * t304 + t302 * t380;
t391 = t215 * t304;
t389 = t242 * t246;
t388 = t242 * t273;
t119 = t213 * t428 - t215 * t427;
t120 = t213 * t427 + t215 * t428;
t382 = t303 ^ 2 * t307;
t25 = m(7) * (t119 * t84 + t120 * t85 - t143) + m(6) * (t163 * t213 + t164 * t215 - t143) + (-t261 * t303 * t306 + t262 * t379 - t382) * t476 + (-t247 * t325 - t309 * t382 - t143) * m(5);
t387 = t25 * qJD(1);
t162 = t273 * t200;
t331 = m(7) * (t324 * t327 + t469 * t484);
t354 = m(6) * t372;
t62 = -t331 / 0.2e1 + 0.2e1 * (-t354 / 0.4e1 - t369) * t273;
t376 = t62 * qJD(2);
t374 = Ifges(7,5) * t324 - Ifges(7,6) * t484;
t373 = -Ifges(7,5) * t469 - Ifges(7,6) * t327;
t371 = t228 * qJD(6);
t370 = t457 / 0.2e1;
t365 = -t406 / 0.2e1;
t364 = -t405 / 0.2e1;
t359 = t200 * t439;
t358 = -t385 / 0.2e1;
t356 = -t378 / 0.2e1;
t227 = -t273 * mrSges(5,1) - t271 * mrSges(5,2);
t350 = -t116 - t218 + t407;
t349 = 0.2e1 * t369 * t325;
t348 = t246 * t478;
t115 = t411 - t413;
t214 = -t263 * t428 - t264 * t427;
t280 = -t304 * pkin(5) + t292;
t311 = (t391 - t392) * t285 * t460 + (t119 * t214 + t120 * t216) * t458 + t119 * t364 + t120 * t365 + t392 * t455 + mrSges(6,3) * t391 / 0.2e1 + t283 * t356 + (mrSges(5,2) / 0.2e1 - t395 * t370) * t247 + (-t292 * t460 - t280 * t458 + mrSges(5,1) / 0.2e1 + t396 * t370 - t472 / 0.2e1) * t246;
t2 = t115 * t444 - t311 + (t227 + t283) * t356 + t464 + t483 * t200;
t167 = -Ifges(6,6) * t273 + (t403 - t416) * t271;
t168 = -Ifges(6,5) * t273 + (-t304 * Ifges(6,1) + t417) * t271;
t229 = mrSges(5,1) * t271 - mrSges(5,2) * t273;
t336 = Ifges(7,5) * t447 + Ifges(7,6) * t450;
t96 = Ifges(7,4) * t197 + Ifges(7,2) * t194 - Ifges(7,6) * t273;
t97 = Ifges(7,2) * t324 + t271 * Ifges(7,6) + t415;
t98 = Ifges(7,1) * t197 + Ifges(7,4) * t194 - Ifges(7,5) * t273;
t191 = Ifges(7,4) * t324;
t99 = Ifges(7,1) * t484 + t271 * Ifges(7,5) + t191;
t3 = (-Ifges(4,4) * t306 + pkin(3) * t229) * t306 + (Ifges(4,4) * t308 + (Ifges(4,1) - Ifges(4,2)) * t306) * t308 + m(6) * (t139 * t141 + t140 * t142 + t471 * t242) + t484 * t98 / 0.2e1 + m(7) * (t173 * t174 + t56 * t58 + t57 * t59) + t426 * t481 + t99 * t447 + t97 * t450 + (Ifges(7,5) * t487 + Ifges(7,6) * t477 - t304 * t168 / 0.2e1 + t167 * t430 + (-Ifges(5,4) + t400 / 0.2e1 - t402 / 0.2e1) * t273) * t273 + t471 * t218 + ((t299 * Ifges(6,1) / 0.2e1 - Ifges(5,2) - Ifges(7,3) + Ifges(5,1) - Ifges(6,3) + (-t416 + t403 / 0.2e1) * t302) * t273 + t336 + (Ifges(5,4) - t400 + t402) * t271) * t271 + t57 * t146 + t59 * t147 + t56 * t148 + t58 * t149 + t173 * t116 + t174 * t115 + t140 * t222 + t142 * t223 + t139 * t224 + t141 * t225 + t242 * t217 - pkin(2) * t283 + t293 * t227 + t324 * t96 / 0.2e1;
t346 = t2 * qJD(1) + t3 * qJD(2);
t114 = t410 + t412;
t117 = -Ifges(7,2) * t484 + t191;
t118 = Ifges(7,1) * t324 - t415;
t6 = t56 * t147 - t57 * t149 + t174 * t114 + t271 * t374 / 0.2e1 + (-t57 * mrSges(7,3) + t118 / 0.2e1 - t97 / 0.2e1) * t484 + (-t56 * mrSges(7,3) + t99 / 0.2e1 + t117 / 0.2e1) * t324;
t319 = (t477 * t84 + t487 * t85) * mrSges(7,3) + t114 * t444 + t84 * t452 - t85 * t149 / 0.2e1;
t334 = t119 * mrSges(7,1) / 0.2e1 + t120 * t456;
t7 = t319 - t334;
t345 = t7 * qJD(1) + t6 * qJD(2);
t11 = t197 * t147 + t194 * t149 + t350 * t273 + (-t377 + t381 + t408) * t271 + m(7) * (-t174 * t273 + t194 * t56 + t197 * t57) + m(6) * (t271 * t342 - t388) + m(5) * (-t271 * t471 - t388);
t341 = t163 * t302 - t164 * t304;
t316 = (t194 * t84 + t197 * t85 - t162) * t458 + (-t271 * t325 - t162) * t462 + (t271 * t341 - t162) * t460;
t318 = (t213 * t304 + t215 * t302) * t460 + (-t119 * t469 + t120 * t327) * t458 + t380 * t462;
t22 = t316 - t318;
t344 = t22 * qJD(1) + t11 * qJD(2);
t23 = m(7) * (t324 * t57 - t484 * t56) + t324 * t147 - t484 * t149 + (m(6) * (t139 * t304 + t140 * t302) + t304 * t225 + t302 * t223) * t273;
t320 = (t324 * t85 - t484 * t84) * t459 + m(6) * (t163 * t304 + t164 * t302) * t435;
t31 = t348 + t320;
t343 = -qJD(1) * t31 + qJD(2) * t23;
t68 = 0.2e1 * t480 + 0.2e1 * t488;
t339 = qJD(2) * t68 - qJD(3) * t228;
t335 = mrSges(7,1) * t454 + t106 * t456;
t333 = t147 * t433 + t149 * t431;
t312 = (-t271 * t351 - t273 * t292) * t460 + (t194 * t214 + t197 * t216 - t273 * t280) * t458 + (-t271 * t395 + t273 * t396) * t370 + t194 * t364 + t197 * t365 + t472 * t435 + t372 * t271 * t455;
t315 = (t141 * t304 + t142 * t302) * t460 + (t327 * t59 - t469 * t58) * t458 + t148 * t433 + t146 * t432 + t222 * t430 + t224 * t429 + t306 * t370;
t10 = -t227 + t312 - t315;
t330 = qJD(1) * t459 * t48 + t10 * qJD(2);
t317 = (-t469 * t477 + t486) * mrSges(7,3) + t333;
t17 = t317 - t474;
t329 = t17 * qJD(2);
t20 = m(7) * (t105 * t84 + t106 * t85 + t123) + m(6) * (t200 * t341 + t123);
t328 = -t20 * qJD(1) - t420 * t48 / 0.2e1;
t314 = (t324 * t433 - t486) * mrSges(7,3) + t342 * t461 + (-t214 * t484 + t216 * t324 - t327 * t56 - t469 * t57) * t458 + t332 + t333;
t14 = (t401 / 0.2e1 + t404 / 0.2e1) * t271 + t314 + t485;
t321 = t341 * t460 + (-t327 * t84 - t469 * t85) * t459;
t35 = t349 + t321;
t60 = (t327 ^ 2 + t469 ^ 2) * mrSges(7,3) + m(7) * (-t214 * t327 - t216 * t469) - t482;
t326 = -t35 * qJD(1) + t14 * qJD(2) + t60 * qJD(3);
t18 = t359 - t335;
t231 = -Ifges(7,2) * t327 - t269;
t233 = -Ifges(7,1) * t469 - t414;
t34 = t280 * t228 - (-t233 / 0.2e1 + t438) * t327 - (t437 + t231 / 0.2e1) * t469;
t313 = -(t117 / 0.4e1 + t99 / 0.4e1) * t469 - (t97 / 0.4e1 - t118 / 0.4e1) * t327 + (t234 / 0.4e1 + t231 / 0.4e1 - t214 * mrSges(7,3) / 0.2e1) * t324 + (mrSges(7,3) * t442 + t233 / 0.4e1 - t232 / 0.4e1) * t484 + t174 * t439 + t214 * t452 + t149 * t442 + t271 * t373 / 0.4e1 + t280 * t114 / 0.2e1;
t322 = Ifges(7,3) * t434 - t58 * mrSges(7,1) / 0.2e1 + t59 * mrSges(7,2) / 0.2e1 - t336;
t5 = t313 + t322;
t323 = -t18 * qJD(1) - t5 * qJD(2) - t34 * qJD(3);
t69 = t412 / 0.2e1 + t410 / 0.2e1 + t488 + t480;
t63 = t331 / 0.2e1 + t354 * t434 + t273 * t478;
t36 = t349 - t321;
t32 = t348 - t320;
t21 = t316 + t318;
t19 = t359 + t335;
t16 = t317 + t474;
t13 = mrSges(6,2) * t358 - mrSges(6,1) * t386 / 0.2e1 + t314 - t485;
t12 = t312 + t315;
t8 = t319 + t334;
t4 = t313 - t322;
t1 = t311 + (-t283 / 0.2e1 - t227 / 0.2e1) * t378 + t464 + (t115 / 0.2e1 + t483) * t200;
t9 = [t25 * qJD(2) + t20 * qJD(3), t1 * qJD(3) + t21 * qJD(4) + t32 * qJD(5) + t8 * qJD(6) + t387 + (t119 * t149 + t120 * t147 + t213 * t225 + t215 * t223 + t247 * t408 + t350 * t246 + 0.2e1 * (t119 * t56 + t120 * t57 - t174 * t246) * t458 + 0.2e1 * (-t247 * t471 - t389) * t462 + 0.2e1 * (t139 * t213 + t140 * t215 - t389) * t460 + ((t470 * mrSges(4,3) - mrSges(3,2)) * t309 + t470 * pkin(8) * t476 + (-m(4) * pkin(2) - mrSges(3,1) + t229 + t467 + t481) * t307) * t303) * qJD(2), t1 * qJD(2) + t36 * qJD(5) + t19 * qJD(6) - t328 + (m(7) * (t105 * t214 + t106 * t216) - t106 * t406 - t105 * t405 - t261 * mrSges(4,2) - t262 * mrSges(4,1) + (m(7) * t280 + t230 + t463) * t325 + (-t465 + t482) * t200) * qJD(3), qJD(2) * t21 + t422, qJD(2) * t32 + qJD(3) * t36, t8 * qJD(2) + t19 * qJD(3) + (-mrSges(7,1) * t85 - mrSges(7,2) * t84) * qJD(6); qJD(3) * t2 + qJD(4) * t22 - qJD(5) * t31 + qJD(6) * t7 - t387, qJD(3) * t3 + qJD(4) * t11 + qJD(5) * t23 + qJD(6) * t6 ((Ifges(6,5) * t302 + Ifges(7,5) * t327 + Ifges(6,6) * t304 - Ifges(7,6) * t469) * t435 + (m(6) * t468 + t304 * t222 - t302 * t224) * t285 + t468 * mrSges(6,3) + t467 * pkin(8) + t167 * t429 + t168 * t430 + t98 * t432 + t96 * t433 + t197 * t437 + t194 * t438 - t58 * t405 - t59 * t406 + t360 * t407 + t362 * t408 + m(7) * (t173 * t280 + t214 * t58 + t216 * t59) + (Ifges(6,1) * t302 + t416) * t358 + (Ifges(6,2) * t304 + t417) * t386 / 0.2e1 - t465 * t242 + t463 * t471 + t214 * t148 + t216 * t146 + t173 * t230 - Ifges(5,5) * t271 + Ifges(5,6) * t273 + t280 * t115 + t292 * t217 - Ifges(4,6) * t306 + Ifges(4,5) * t308) * qJD(3) + t12 * qJD(4) + t13 * qJD(5) + t4 * qJD(6) + t346, t12 * qJD(3) + (-t194 * t469 + t197 * t327) * t420 + t63 * qJD(5) + t16 * qJD(6) + t344, qJD(3) * t13 + qJD(4) * t63 + qJD(6) * t69 + t343, t4 * qJD(3) + t16 * qJD(4) + t69 * qJD(5) + (-mrSges(7,1) * t57 - mrSges(7,2) * t56 + t374) * qJD(6) + t345; -t2 * qJD(2) - t35 * qJD(5) + t18 * qJD(6) + t328, qJD(4) * t10 + qJD(5) * t14 + qJD(6) * t5 - t346, qJD(5) * t60 + qJD(6) * t34, t330, t326 (-mrSges(7,1) * t216 - mrSges(7,2) * t214 + t373) * qJD(6) - t323; -qJD(2) * t22 + t422, -qJD(3) * t10 - qJD(5) * t62 + qJD(6) * t17 - t344, -t330, 0, -t376, t329 - t371; qJD(2) * t31 + qJD(3) * t35, -qJD(3) * t14 + qJD(4) * t62 - qJD(6) * t68 - t343, -t326 + t371, t376, 0, -t339; -t7 * qJD(2) - t18 * qJD(3), -qJD(3) * t5 - qJD(4) * t17 + qJD(5) * t68 - t345, -t228 * qJD(5) + t323, -t329, t339, 0;];
Cq  = t9;
