% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:57
% EndTime: 2019-03-09 17:04:58
% DurationCPUTime: 32.12s
% Computational Cost: add. (18308->841), mult. (48602->1146), div. (0->0), fcn. (38098->10), ass. (0->352)
t312 = cos(qJ(2));
t309 = sin(qJ(2));
t304 = sin(pkin(6));
t383 = qJD(1) * t304;
t365 = t309 * t383;
t306 = cos(pkin(6));
t382 = qJD(1) * t306;
t374 = pkin(1) * t382;
t259 = -pkin(8) * t365 + t312 * t374;
t320 = (pkin(2) * t309 - pkin(9) * t312) * t304;
t260 = qJD(1) * t320;
t308 = sin(qJ(3));
t311 = cos(qJ(3));
t195 = -t308 * t259 + t311 * t260;
t423 = -qJ(4) - pkin(9);
t357 = qJD(3) * t423;
t507 = -(-qJ(4) * t311 * t312 + pkin(3) * t309) * t383 - t195 - qJD(4) * t308 + t311 * t357;
t196 = t311 * t259 + t308 * t260;
t364 = t312 * t383;
t355 = t308 * t364;
t506 = -qJ(4) * t355 - qJD(4) * t311 - t308 * t357 + t196;
t303 = sin(pkin(11));
t305 = cos(pkin(11));
t276 = t303 * t311 + t305 * t308;
t221 = t276 * t364;
t268 = t276 * qJD(3);
t499 = t221 - t268;
t485 = t303 * t507 - t506 * t305;
t262 = pkin(8) * t364 + t309 * t374;
t219 = pkin(3) * t355 + t262;
t379 = qJD(3) * t308;
t505 = pkin(3) * t379 - t219;
t504 = pkin(10) * t365 - t485;
t275 = t303 * t308 - t305 * t311;
t222 = t275 * t364;
t269 = t275 * qJD(3);
t503 = t505 + (-t222 + t269) * pkin(10) - t499 * pkin(4);
t294 = qJD(2) + t382;
t229 = pkin(9) * t294 + t262;
t255 = (-pkin(2) * t312 - pkin(9) * t309 - pkin(1)) * t304;
t236 = qJD(1) * t255;
t177 = t229 * t311 + t236 * t308;
t261 = qJD(2) * t320;
t249 = qJD(1) * t261;
t390 = t304 * t309;
t295 = pkin(8) * t390;
t426 = pkin(1) * t312;
t272 = t306 * t426 - t295;
t263 = t272 * qJD(2);
t250 = qJD(1) * t263;
t125 = -qJD(3) * t177 + t311 * t249 - t250 * t308;
t380 = qJD(2) * t312;
t361 = t311 * t380;
t378 = qJD(3) * t311;
t210 = t294 * t378 + (-t309 * t379 + t361) * t383;
t242 = t294 * t308 + t311 * t365;
t381 = qJD(2) * t304;
t359 = qJD(1) * t381;
t354 = t309 * t359;
t70 = pkin(3) * t354 - qJ(4) * t210 - qJD(4) * t242 + t125;
t124 = -t229 * t379 + t236 * t378 + t308 * t249 + t311 * t250;
t362 = t308 * t380;
t211 = -t294 * t379 + (-t309 * t378 - t362) * t383;
t241 = t294 * t311 - t308 * t365;
t78 = qJ(4) * t211 + qJD(4) * t241 + t124;
t25 = t303 * t70 + t305 * t78;
t502 = t25 * mrSges(5,3);
t479 = Ifges(7,1) + Ifges(6,1);
t478 = Ifges(7,4) + Ifges(6,5);
t501 = Ifges(6,6) - Ifges(7,6);
t389 = t304 * t312;
t273 = t306 * t309 * pkin(1) + pkin(8) * t389;
t264 = t273 * qJD(2);
t251 = qJD(1) * t264;
t173 = -t211 * pkin(3) + t251;
t500 = t173 * mrSges(5,1);
t487 = t506 * t303 + t305 * t507;
t287 = qJD(3) - t364;
t394 = t287 * Ifges(5,6);
t321 = t241 * t303 + t305 * t242;
t401 = t321 * Ifges(5,4);
t356 = t305 * t241 - t242 * t303;
t404 = t356 * Ifges(5,2);
t122 = t394 + t401 + t404;
t181 = qJD(5) - t356;
t307 = sin(qJ(5));
t310 = cos(qJ(5));
t176 = -t229 * t308 + t311 * t236;
t147 = -qJ(4) * t242 + t176;
t130 = pkin(3) * t287 + t147;
t148 = qJ(4) * t241 + t177;
t388 = t305 * t148;
t68 = t303 * t130 + t388;
t63 = pkin(10) * t287 + t68;
t228 = -t294 * pkin(2) - t259;
t182 = -t241 * pkin(3) + qJD(4) + t228;
t90 = -pkin(4) * t356 - pkin(10) * t321 + t182;
t26 = -t307 * t63 + t310 * t90;
t470 = qJD(6) - t26;
t17 = -pkin(5) * t181 + t470;
t27 = t307 * t90 + t310 * t63;
t18 = qJ(6) * t181 + t27;
t498 = t182 * mrSges(5,1) + t26 * mrSges(6,1) - t17 * mrSges(7,1) - t27 * mrSges(6,2) - t68 * mrSges(5,3) + t18 * mrSges(7,3) - t122 / 0.2e1;
t432 = t287 / 0.2e1;
t497 = t321 / 0.2e1;
t496 = t356 / 0.2e1;
t477 = Ifges(7,5) - Ifges(6,4);
t495 = Ifges(4,3) + Ifges(5,3);
t302 = -pkin(3) * t311 - pkin(2);
t212 = pkin(4) * t275 - pkin(10) * t276 + t302;
t288 = t423 * t308;
t289 = t423 * t311;
t227 = t288 * t303 - t289 * t305;
t376 = qJD(5) * t310;
t377 = qJD(5) * t307;
t490 = t212 * t376 - t227 * t377 + t503 * t307 - t310 * t504;
t469 = t307 * t212 + t310 * t227;
t489 = -qJD(5) * t469 + t307 * t504 + t503 * t310;
t486 = pkin(4) * t365 - t487;
t158 = -t310 * t287 + t307 * t321;
t159 = t287 * t307 + t310 * t321;
t138 = t303 * t148;
t67 = t130 * t305 - t138;
t62 = -pkin(4) * t287 - t67;
t30 = pkin(5) * t158 - qJ(6) * t159 + t62;
t329 = t26 * t310 + t27 * t307;
t332 = t17 * t310 - t18 * t307;
t409 = Ifges(7,5) * t310;
t336 = Ifges(7,3) * t307 + t409;
t413 = Ifges(6,4) * t310;
t342 = -Ifges(6,2) * t307 + t413;
t347 = mrSges(7,1) * t307 - mrSges(7,3) * t310;
t349 = mrSges(6,1) * t307 + mrSges(6,2) * t310;
t428 = t310 / 0.2e1;
t430 = t307 / 0.2e1;
t431 = -t307 / 0.2e1;
t446 = t181 / 0.2e1;
t450 = t159 / 0.2e1;
t452 = t158 / 0.2e1;
t453 = -t158 / 0.2e1;
t410 = Ifges(7,5) * t307;
t414 = Ifges(6,4) * t307;
t466 = t310 * t479 + t410 - t414;
t467 = -t307 * t501 + t310 * t478;
t155 = Ifges(6,4) * t158;
t411 = Ifges(7,5) * t158;
t472 = t159 * t479 + t181 * t478 - t155 + t411;
t154 = Ifges(7,5) * t159;
t56 = Ifges(7,6) * t181 + Ifges(7,3) * t158 + t154;
t415 = Ifges(6,4) * t159;
t59 = -Ifges(6,2) * t158 + Ifges(6,6) * t181 + t415;
t494 = t332 * mrSges(7,2) - t329 * mrSges(6,3) + t30 * t347 + t336 * t452 + t342 * t453 + t349 * t62 + t428 * t472 + t430 * t56 + t431 * t59 + t446 * t467 + t450 * t466;
t476 = Ifges(7,2) + Ifges(6,3);
t493 = -qJ(6) * t499 + qJD(6) * t275 + t490;
t57 = t159 * Ifges(6,5) - t158 * Ifges(6,6) + t181 * Ifges(6,3);
t58 = t159 * Ifges(7,4) + t181 * Ifges(7,2) + t158 * Ifges(7,6);
t492 = t58 + t57;
t491 = pkin(5) * t499 - t489;
t200 = -t222 * t307 - t310 * t365;
t201 = -t222 * t310 + t307 * t365;
t333 = pkin(5) * t307 - qJ(6) * t310;
t334 = pkin(5) * t310 + qJ(6) * t307;
t488 = -pkin(5) * t200 + qJ(6) * t201 - t333 * t269 + (qJD(5) * t334 - qJD(6) * t310) * t276 + t486;
t455 = Ifges(5,1) * t497 + Ifges(5,4) * t496 + Ifges(5,5) * t432;
t484 = t182 * mrSges(5,2) - t67 * mrSges(5,3) + 0.2e1 * t455 + t494;
t483 = t307 * t478 + t310 * t501;
t482 = t307 * t479 - t409 + t413;
t481 = qJD(5) * t276;
t480 = t30 * mrSges(7,1) + t62 * mrSges(6,1) + t56 / 0.2e1 - t59 / 0.2e1 - t18 * mrSges(7,2) - t27 * mrSges(6,3);
t149 = t210 * t303 - t305 * t211;
t150 = t210 * t305 + t211 * t303;
t84 = -qJD(5) * t158 + t150 * t310 + t307 * t354;
t85 = qJD(5) * t159 + t150 * t307 - t310 * t354;
t12 = Ifges(6,5) * t84 - Ifges(6,6) * t85 + Ifges(6,3) * t149;
t13 = Ifges(7,4) * t84 + Ifges(7,2) * t149 + Ifges(7,6) * t85;
t474 = t13 + t12;
t473 = t149 * t478 + t477 * t85 + t479 * t84;
t270 = t306 * t311 - t308 * t390;
t271 = t306 * t308 + t311 * t390;
t203 = -t305 * t270 + t271 * t303;
t204 = t270 * t303 + t271 * t305;
t253 = t295 + (-pkin(2) - t426) * t306;
t209 = -t270 * pkin(3) + t253;
t120 = t203 * pkin(4) - t204 * pkin(10) + t209;
t254 = pkin(9) * t306 + t273;
t191 = -t308 * t254 + t311 * t255;
t153 = -pkin(3) * t389 - t271 * qJ(4) + t191;
t192 = t311 * t254 + t308 * t255;
t167 = qJ(4) * t270 + t192;
t99 = t303 * t153 + t305 * t167;
t94 = -pkin(10) * t389 + t99;
t471 = t307 * t120 + t310 * t94;
t465 = Ifges(4,5) * t210 + Ifges(5,5) * t150 + Ifges(4,6) * t211 - Ifges(5,6) * t149 + t354 * t495;
t398 = t242 * Ifges(4,4);
t171 = t241 * Ifges(4,2) + t287 * Ifges(4,6) + t398;
t237 = Ifges(4,4) * t241;
t172 = t242 * Ifges(4,1) + t287 * Ifges(4,5) + t237;
t324 = t176 * t311 + t177 * t308;
t416 = Ifges(4,4) * t311;
t417 = Ifges(4,4) * t308;
t427 = t311 / 0.2e1;
t436 = t242 / 0.2e1;
t437 = t241 / 0.2e1;
t464 = -t324 * mrSges(4,3) + t228 * (mrSges(4,1) * t308 + mrSges(4,2) * t311) + (-Ifges(4,2) * t308 + t416) * t437 + (Ifges(4,1) * t311 - t417) * t436 + (Ifges(4,5) * t311 - Ifges(4,6) * t308) * t432 - t308 * t171 / 0.2e1 + t172 * t427;
t21 = pkin(10) * t354 + t25;
t53 = t149 * pkin(4) - t150 * pkin(10) + t173;
t4 = -qJD(5) * t27 - t21 * t307 + t310 * t53;
t363 = t309 * t381;
t131 = -t254 * t379 + t255 * t378 + t308 * t261 + t311 * t263;
t217 = -qJD(3) * t271 - t304 * t362;
t100 = qJ(4) * t217 + qJD(4) * t270 + t131;
t132 = -qJD(3) * t192 + t311 * t261 - t263 * t308;
t218 = qJD(3) * t270 + t304 * t361;
t95 = pkin(3) * t363 - qJ(4) * t218 - qJD(4) * t271 + t132;
t39 = t305 * t100 + t303 * t95;
t37 = pkin(10) * t363 + t39;
t163 = -t305 * t217 + t218 * t303;
t164 = t217 * t303 + t218 * t305;
t189 = -t217 * pkin(3) + t264;
t65 = t163 * pkin(4) - t164 * pkin(10) + t189;
t8 = -qJD(5) * t471 - t307 * t37 + t310 * t65;
t370 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t371 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t372 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t463 = t158 * t370 - t159 * t372 - t181 * t371 - t57 / 0.2e1 - t58 / 0.2e1 + t404 / 0.2e1 + t401 / 0.2e1 + t394 / 0.2e1 - t498;
t462 = t84 / 0.2e1;
t461 = -t85 / 0.2e1;
t460 = t85 / 0.2e1;
t459 = pkin(1) * mrSges(3,1);
t458 = pkin(1) * mrSges(3,2);
t454 = t149 / 0.2e1;
t451 = -t159 / 0.2e1;
t447 = -t181 / 0.2e1;
t444 = -t203 / 0.2e1;
t442 = t204 / 0.2e1;
t441 = t210 / 0.2e1;
t440 = t211 / 0.2e1;
t435 = t270 / 0.2e1;
t434 = t271 / 0.2e1;
t433 = -t287 / 0.2e1;
t429 = -t310 / 0.2e1;
t425 = pkin(3) * t242;
t46 = mrSges(6,1) * t149 - mrSges(6,3) * t84;
t47 = -t149 * mrSges(7,1) + t84 * mrSges(7,2);
t422 = -t46 + t47;
t48 = -mrSges(6,2) * t149 - mrSges(6,3) * t85;
t49 = -mrSges(7,2) * t85 + mrSges(7,3) * t149;
t421 = t48 + t49;
t420 = mrSges(6,3) * t158;
t419 = mrSges(6,3) * t159;
t418 = Ifges(3,4) * t309;
t412 = Ifges(3,5) * t312;
t408 = t149 * Ifges(5,4);
t407 = t150 * Ifges(5,1);
t406 = t150 * Ifges(5,4);
t403 = t356 * Ifges(5,6);
t400 = t321 * Ifges(5,5);
t399 = t241 * Ifges(4,6);
t397 = t242 * Ifges(4,5);
t396 = t250 * mrSges(3,2);
t393 = t294 * Ifges(3,5);
t111 = pkin(4) * t321 - pkin(10) * t356 + t425;
t77 = t147 * t305 - t138;
t35 = t307 * t111 + t310 * t77;
t300 = pkin(3) * t303 + pkin(10);
t392 = t300 * t307;
t391 = t300 * t310;
t114 = -mrSges(7,2) * t158 + mrSges(7,3) * t181;
t115 = -mrSges(6,2) * t181 - t420;
t387 = t114 + t115;
t116 = mrSges(6,1) * t181 - t419;
t117 = -mrSges(7,1) * t181 + mrSges(7,2) * t159;
t386 = t116 - t117;
t103 = mrSges(6,1) * t158 + mrSges(6,2) * t159;
t166 = mrSges(5,1) * t287 - mrSges(5,3) * t321;
t385 = t166 - t103;
t384 = -mrSges(3,1) * t294 - mrSges(4,1) * t241 + mrSges(4,2) * t242 + mrSges(3,3) * t365;
t369 = -Ifges(4,3) / 0.2e1 - Ifges(5,3) / 0.2e1;
t368 = t307 * t389;
t301 = -pkin(3) * t305 - pkin(4);
t24 = -t303 * t78 + t305 * t70;
t38 = -t303 * t100 + t305 * t95;
t86 = t149 * mrSges(5,1) + t150 * mrSges(5,2);
t76 = t147 * t303 + t388;
t98 = t153 * t305 - t303 * t167;
t226 = -t305 * t288 - t289 * t303;
t3 = t310 * t21 + t307 * t53 + t90 * t376 - t377 * t63;
t1 = qJ(6) * t149 + qJD(6) * t181 + t3;
t2 = -pkin(5) * t149 - t4;
t352 = -t1 * t307 + t2 * t310;
t351 = -t3 * t307 - t310 * t4;
t93 = pkin(4) * t389 - t98;
t350 = mrSges(6,1) * t310 - mrSges(6,2) * t307;
t348 = mrSges(7,1) * t310 + mrSges(7,3) * t307;
t341 = Ifges(6,2) * t310 + t414;
t335 = -Ifges(7,3) * t310 + t410;
t331 = t17 * t307 + t18 * t310;
t328 = t26 * t307 - t27 * t310;
t34 = t111 * t310 - t307 * t77;
t41 = t120 * t310 - t307 * t94;
t325 = t124 * t311 - t125 * t308;
t156 = t212 * t310 - t227 * t307;
t179 = t307 * t204 + t310 * t389;
t7 = t120 * t376 + t307 * t65 + t310 * t37 - t377 * t94;
t36 = -pkin(4) * t363 - t38;
t317 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t20 = -pkin(4) * t354 - t24;
t290 = Ifges(3,4) * t364;
t286 = t359 * t412;
t274 = t301 - t334;
t265 = qJD(5) * t333 - qJD(6) * t307;
t258 = -t294 * mrSges(3,2) + mrSges(3,3) * t364;
t224 = Ifges(3,1) * t365 + t290 + t393;
t223 = Ifges(3,6) * t294 + (Ifges(3,2) * t312 + t418) * t383;
t214 = mrSges(4,1) * t287 - mrSges(4,3) * t242;
t213 = -mrSges(4,2) * t287 + mrSges(4,3) * t241;
t188 = -mrSges(4,2) * t354 + mrSges(4,3) * t211;
t187 = mrSges(4,1) * t354 - mrSges(4,3) * t210;
t180 = t310 * t204 - t368;
t170 = t287 * Ifges(4,3) + t397 + t399;
t169 = t276 * t333 + t226;
t165 = -mrSges(5,2) * t287 + mrSges(5,3) * t356;
t151 = -mrSges(4,1) * t211 + mrSges(4,2) * t210;
t141 = -pkin(5) * t275 - t156;
t137 = qJ(6) * t275 + t469;
t136 = t210 * Ifges(4,1) + t211 * Ifges(4,4) + Ifges(4,5) * t354;
t135 = t210 * Ifges(4,4) + t211 * Ifges(4,2) + Ifges(4,6) * t354;
t129 = mrSges(5,1) * t354 - mrSges(5,3) * t150;
t128 = -mrSges(5,2) * t354 - mrSges(5,3) * t149;
t126 = -mrSges(5,1) * t356 + mrSges(5,2) * t321;
t121 = t287 * Ifges(5,3) + t400 + t403;
t107 = -qJD(5) * t368 + t164 * t307 + t204 * t376 - t310 * t363;
t106 = -qJD(5) * t179 + t310 * t164 + t307 * t363;
t102 = mrSges(7,1) * t158 - mrSges(7,3) * t159;
t101 = pkin(5) * t159 + qJ(6) * t158;
t72 = Ifges(5,5) * t354 + t407 - t408;
t71 = -t149 * Ifges(5,2) + Ifges(5,6) * t354 + t406;
t43 = pkin(5) * t179 - qJ(6) * t180 + t93;
t40 = t333 * t356 + t76;
t32 = -pkin(5) * t203 - t41;
t31 = qJ(6) * t203 + t471;
t29 = mrSges(6,1) * t85 + mrSges(6,2) * t84;
t28 = mrSges(7,1) * t85 - mrSges(7,3) * t84;
t23 = -pkin(5) * t321 - t34;
t22 = qJ(6) * t321 + t35;
t14 = t84 * Ifges(6,4) - t85 * Ifges(6,2) + t149 * Ifges(6,6);
t11 = t84 * Ifges(7,5) + t149 * Ifges(7,6) + t85 * Ifges(7,3);
t10 = pkin(5) * t107 - qJ(6) * t106 - qJD(6) * t180 + t36;
t9 = pkin(5) * t85 - qJ(6) * t84 - qJD(6) * t159 + t20;
t6 = -pkin(5) * t163 - t8;
t5 = qJ(6) * t163 + qJD(6) * t203 + t7;
t15 = [t251 * (-mrSges(4,1) * t270 + mrSges(4,2) * t271) + t263 * t258 - t465 * t389 / 0.2e1 + (Ifges(6,4) * t180 + Ifges(6,6) * t203) * t461 + (Ifges(5,1) * t164 + Ifges(5,5) * t363) * t497 + (Ifges(5,4) * t164 + Ifges(5,6) * t363) * t496 + (-Ifges(6,2) * t453 + Ifges(7,3) * t452 - t446 * t501 + t477 * t450 + t480) * t107 + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) + t9 * mrSges(7,1) + t20 * mrSges(6,1) + t11 / 0.2e1 - t14 / 0.2e1 + Ifges(7,3) * t460 - Ifges(6,2) * t461 + t477 * t462 - t501 * t454) * t179 - t203 * t502 + (Ifges(6,4) * t453 - t30 * mrSges(7,3) - t26 * mrSges(6,3) + t17 * mrSges(7,2) + Ifges(7,5) * t452 + t472 / 0.2e1 + t478 * t446 + t479 * t450 + t62 * mrSges(6,2)) * t106 + (t164 * t182 + t173 * t204 + t25 * t389 - t363 * t68) * mrSges(5,2) + (t250 * t312 + t251 * t309 + (-t259 * t312 - t262 * t309) * qJD(2)) * t304 * mrSges(3,3) + t253 * t151 + t228 * (-mrSges(4,1) * t217 + mrSges(4,2) * t218) + t218 * t172 / 0.2e1 + t217 * t171 / 0.2e1 + t131 * t213 + t132 * t214 + t209 * t86 + t4 * (mrSges(6,1) * t203 - mrSges(6,3) * t180) + t2 * (-mrSges(7,1) * t203 + mrSges(7,2) * t180) + t191 * t187 + t192 * t188 + t189 * t126 + t38 * t166 + t39 * t165 + ((t170 + t121) * t309 + t312 * t224 + t294 * (-Ifges(3,6) * t309 + t412)) * t381 / 0.2e1 + t471 * t48 + m(6) * (t20 * t93 + t26 * t8 + t27 * t7 + t3 * t471 + t36 * t62 + t4 * t41) + t99 * t128 + t98 * t129 + t5 * t114 + t7 * t115 + t8 * t116 + t6 * t117 + t10 * t102 + t36 * t103 + t93 * t29 + t41 * t46 + t32 * t47 + t31 * t49 + t43 * t28 + ((-t272 * mrSges(3,3) + Ifges(3,5) * t306 / 0.2e1 + (-0.2e1 * t458 + 0.3e1 / 0.2e1 * Ifges(3,4) * t312) * t304) * t312 + (-t273 * mrSges(3,3) + Ifges(4,5) * t434 + Ifges(4,6) * t435 + Ifges(5,5) * t442 + Ifges(5,6) * t444 - Ifges(3,6) * t306 + (-0.2e1 * t459 - 0.3e1 / 0.2e1 * t418) * t304 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + t369) * t389) * t309) * t359 + (Ifges(7,5) * t180 + Ifges(7,6) * t203) * t460 + (t286 / 0.2e1 - t251 * mrSges(3,1) - t396) * t306 - t149 * (Ifges(5,4) * t204 - Ifges(5,2) * t203 - Ifges(5,6) * t389) / 0.2e1 + t150 * (Ifges(5,1) * t204 - Ifges(5,4) * t203 - Ifges(5,5) * t389) / 0.2e1 + t24 * (-mrSges(5,1) * t389 - t204 * mrSges(5,3)) + t125 * (-mrSges(4,1) * t389 - t271 * mrSges(4,3)) + t124 * (mrSges(4,2) * t389 + t270 * mrSges(4,3)) + t384 * t264 + (t1 * t203 - t180 * t9) * mrSges(7,3) + (t180 * t20 - t203 * t3) * mrSges(6,2) + t177 * (-mrSges(4,2) * t363 + mrSges(4,3) * t217) + t67 * (mrSges(5,1) * t363 - mrSges(5,3) * t164) + t176 * (mrSges(4,1) * t363 - mrSges(4,3) * t218) - t223 * t363 / 0.2e1 + t136 * t434 + t135 * t435 + (Ifges(4,1) * t218 + Ifges(4,4) * t217 + Ifges(4,5) * t363) * t436 + (Ifges(4,4) * t218 + Ifges(4,2) * t217 + Ifges(4,6) * t363) * t437 + t473 * t180 / 0.2e1 + t474 * t203 / 0.2e1 + (t180 * t478 + t203 * t476) * t454 + (t180 * t479 + t203 * t478) * t462 + (Ifges(6,6) * t453 + t492 / 0.2e1 + Ifges(7,6) * t452 + t476 * t446 + t478 * t450 - Ifges(5,6) * t432 - Ifges(5,2) * t496 - Ifges(5,4) * t497 + t498) * t163 + (Ifges(4,4) * t271 + Ifges(4,2) * t270 - Ifges(4,6) * t389) * t440 + (Ifges(4,1) * t271 + Ifges(4,4) * t270 - Ifges(4,5) * t389) * t441 + t72 * t442 + t71 * t444 + t164 * t455 + t203 * t500 + m(3) * (t250 * t273 - t251 * t272 - t259 * t264 + t262 * t263) + m(4) * (t124 * t192 + t125 * t191 + t131 * t177 + t132 * t176 + t228 * t264 + t251 * t253) + m(5) * (t173 * t209 + t182 * t189 + t24 * t98 + t25 * t99 + t38 * t67 + t39 * t68) + m(7) * (t1 * t31 + t10 * t30 + t17 * t6 + t18 * t5 + t2 * t32 + t43 * t9) + (Ifges(4,5) * t218 + Ifges(5,5) * t164 + Ifges(4,6) * t217 + t363 * t495) * t432; (-t201 * t62 + t221 * t27) * mrSges(6,2) + (-t18 * t221 + t201 * t30) * mrSges(7,3) + (t72 / 0.2e1 + t173 * mrSges(5,2) - t408 / 0.2e1 + t407 / 0.2e1 - t24 * mrSges(5,3) + t11 * t430 + t336 * t460 + t342 * t461 + t20 * t349 + t9 * t347 + t351 * mrSges(6,3) + t352 * mrSges(7,2) + (-mrSges(7,2) * t331 + mrSges(6,3) * t328 + t30 * t348 + t335 * t453 + t341 * t452 + t350 * t62 + t429 * t59) * qJD(5) + t466 * t462 + t467 * t454 + (qJD(5) * t472 + t14) * t431 + (qJD(5) * t56 + t473) * t428) * t276 + t488 * t102 + (t1 * t137 + t141 * t2 + t169 * t9 + t17 * t491 + t18 * t493 + t30 * t488) * m(7) + (t201 * t478 + t221 * t476 + t481 * t483) * t447 + (t201 * t479 + t478 * t221 + t482 * t481) * t451 - t356 * (-Ifges(5,4) * t222 - Ifges(5,2) * t221) / 0.2e1 - t321 * (-Ifges(5,1) * t222 - Ifges(5,4) * t221) / 0.2e1 + (t221 * t68 - t222 * t67) * mrSges(5,3) - t182 * (mrSges(5,1) * t221 - mrSges(5,2) * t222) + (-Ifges(5,5) * t222 - Ifges(5,6) * t221) * t433 + t485 * t165 + t486 * t103 + t487 * t166 + t286 + (-Ifges(6,2) * t452 + Ifges(7,3) * t453 - t447 * t501 + t451 * t477 - t480) * t200 + (t308 * pkin(3) * t126 + (-t213 * t308 - t214 * t311) * pkin(9) + t464) * qJD(3) + ((-t393 / 0.2e1 - t224 / 0.2e1 - t290 / 0.2e1 + t259 * mrSges(3,3) + t383 * t458 - t464) * t312 + (t177 * mrSges(4,2) - t176 * mrSges(4,1) - t397 / 0.2e1 - t399 / 0.2e1 + t223 / 0.2e1 - t121 / 0.2e1 - t170 / 0.2e1 + t262 * mrSges(3,3) - t400 / 0.2e1 - t403 / 0.2e1 - t67 * mrSges(5,1) + t68 * mrSges(5,2) + t369 * t287 + (t459 + t418 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t312) * t383 + (t294 / 0.2e1 - qJD(2)) * Ifges(3,6) + (Ifges(4,5) * t308 + Ifges(5,5) * t276 + Ifges(4,6) * t311 - Ifges(5,6) * t275) * qJD(2) / 0.2e1) * t309) * t383 - t396 + (-mrSges(4,1) * t311 + mrSges(4,2) * t308 - mrSges(3,1)) * t251 + t489 * t116 + (t156 * t4 + t20 * t226 + t26 * t489 + t27 * t490 + t469 * t3 + t486 * t62) * m(6) + t490 * t115 + t491 * t117 - t472 * t201 / 0.2e1 - t492 * t221 / 0.2e1 + t493 * t114 + (t13 / 0.2e1 - t71 / 0.2e1 + t500 + t12 / 0.2e1 - t502 - t406 / 0.2e1 - t370 * t85 + t372 * t84 + (Ifges(5,2) / 0.2e1 + t371) * t149 + t317) * t275 + (t29 - t129) * t226 + (Ifges(6,4) * t201 + Ifges(6,6) * t221) * t452 + (Ifges(7,5) * t201 + Ifges(7,6) * t221) * t453 + (-t187 * t308 + t188 * t311) * pkin(9) - t463 * t268 - t259 * t258 + t227 * t128 + t221 * t122 / 0.2e1 - t219 * t126 - t26 * (mrSges(6,1) * t221 - mrSges(6,3) * t201) - t17 * (-mrSges(7,1) * t221 + mrSges(7,2) * t201) - t196 * t213 - t195 * t214 + t169 * t28 + t156 * t46 + t469 * t48 - pkin(2) * t151 + t137 * t49 + t141 * t47 + (t173 * t302 + t182 * t505 - t226 * t24 + t227 * t25 + t485 * t68 + t487 * t67) * m(5) - t384 * t262 + t302 * t86 + t135 * t427 + t308 * t136 / 0.2e1 + (-pkin(2) * t251 - t176 * t195 - t177 * t196 - t228 * t262 + (-qJD(3) * t324 + t325) * pkin(9)) * m(4) + (Ifges(4,2) * t311 + t417) * t440 + (Ifges(4,1) * t308 + t416) * t441 + t325 * mrSges(4,3) + t222 * t455 - t484 * t269; -m(6) * (t26 * t34 + t27 * t35 + t62 * t76) - m(7) * (t17 * t23 + t18 * t22 + t30 * t40) - (-Ifges(4,2) * t242 + t172 + t237) * t241 / 0.2e1 + t274 * t28 + ((t24 * t305 + t25 * t303) * pkin(3) - t182 * t425 + t67 * t76 - t68 * t77) * m(5) + (-t126 * t242 + t128 * t303 + t129 * t305) * pkin(3) + t465 + t463 * t321 + (t3 * t310 - t307 * t4) * mrSges(6,3) + (t265 - t40) * t102 - t228 * (mrSges(4,1) * t242 + mrSges(4,2) * t241) - t176 * t213 + t177 * t214 - t77 * t165 - t124 * mrSges(4,2) + t125 * mrSges(4,1) - t22 * t114 - t35 * t115 - t34 * t116 - t23 * t117 + t24 * mrSges(5,1) - t25 * mrSges(5,2) + (t307 * t422 + t310 * t421) * t300 - t242 * (Ifges(4,1) * t241 - t398) / 0.2e1 + m(6) * (t20 * t301 + t3 * t391 - t392 * t4) + m(7) * (t1 * t391 + t2 * t392 + t265 * t30 + t274 * t9) + t385 * t76 + t14 * t428 + t11 * t429 + (Ifges(4,5) * t241 - Ifges(4,6) * t242) * t433 + t171 * t436 + t473 * t430 + t301 * t29 + t335 * t460 + t341 * t461 + (t1 * t310 + t2 * t307) * mrSges(7,2) + ((-m(6) * t329 + m(7) * t332 - t307 * t387 - t310 * t386) * t300 + t494) * qJD(5) - t9 * t348 - t20 * t350 + (t176 * t241 + t177 * t242) * mrSges(4,3) + t482 * t462 + t483 * t454 - t484 * t356; -t356 * t165 + (-t102 + t385) * t321 + (t181 * t387 - t422) * t310 + (-t181 * t386 + t421) * t307 + t86 + (t181 * t331 - t321 * t30 - t352) * m(7) + (-t181 * t328 - t321 * t62 - t351) * m(6) + (t321 * t67 - t356 * t68 + t173) * m(5); (t158 * t17 + t159 * t18) * mrSges(7,2) + t317 + (t386 + t419) * t27 + (-t387 - t420) * t26 - t30 * (mrSges(7,1) * t159 + mrSges(7,3) * t158) - t62 * (mrSges(6,1) * t159 - mrSges(6,2) * t158) + t59 * t450 + (Ifges(7,3) * t159 - t411) * t453 + qJD(6) * t114 - t101 * t102 - pkin(5) * t47 + qJ(6) * t49 + (-t158 * t478 - t159 * t501) * t447 + (-pkin(5) * t2 + qJ(6) * t1 - t101 * t30 - t17 * t27 + t18 * t470) * m(7) + (-Ifges(6,2) * t159 - t155 + t472) * t452 + (-t158 * t479 + t154 - t415 + t56) * t451 + t474; t159 * t102 - t181 * t114 + 0.2e1 * (t2 / 0.2e1 + t30 * t450 + t18 * t447) * m(7) + t47;];
tauc  = t15(:);
