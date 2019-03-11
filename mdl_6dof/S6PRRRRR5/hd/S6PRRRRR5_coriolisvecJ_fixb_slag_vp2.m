% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:03:56
% EndTime: 2019-03-09 01:04:47
% DurationCPUTime: 23.98s
% Computational Cost: add. (18328->828), mult. (48206->1182), div. (0->0), fcn. (39193->14), ass. (0->378)
t302 = cos(pkin(7));
t312 = cos(qJ(3));
t313 = cos(qJ(2));
t386 = t312 * t313;
t307 = sin(qJ(3));
t308 = sin(qJ(2));
t391 = t307 * t308;
t329 = -t302 * t391 + t386;
t301 = sin(pkin(6));
t377 = qJD(1) * t301;
t224 = t329 * t377;
t300 = sin(pkin(7));
t396 = t300 * t307;
t293 = pkin(9) * t396;
t392 = t302 * t312;
t265 = pkin(2) * t392 - t293;
t255 = t265 * qJD(3);
t506 = -t224 + t255;
t393 = t302 * t307;
t395 = t300 * t312;
t266 = pkin(2) * t393 + pkin(9) * t395;
t248 = pkin(10) * t302 + t266;
t348 = -pkin(3) * t312 - pkin(10) * t307;
t249 = (-pkin(2) + t348) * t300;
t347 = pkin(3) * t307 - pkin(10) * t312;
t328 = t347 * qJD(3);
t254 = t300 * t328;
t306 = sin(qJ(4));
t311 = cos(qJ(4));
t362 = t308 * t377;
t354 = t300 * t362;
t370 = qJD(4) * t311;
t371 = qJD(4) * t306;
t471 = -t248 * t371 + t249 * t370 + t506 * t311 + (t254 - t354) * t306;
t389 = t308 * t312;
t390 = t307 * t313;
t331 = t302 * t389 + t390;
t223 = t331 * t377;
t256 = t266 * qJD(3);
t470 = t223 - t256;
t373 = qJD(3) * t300;
t359 = t307 * t373;
t505 = pkin(11) * t359 + t471;
t261 = -t311 * t302 + t306 * t396;
t372 = qJD(3) * t312;
t357 = t311 * t372;
t213 = -qJD(4) * t261 + t300 * t357;
t262 = t302 * t306 + t311 * t396;
t358 = t306 * t372;
t214 = qJD(4) * t262 + t300 * t358;
t504 = pkin(4) * t214 - pkin(11) * t213 - t470;
t305 = sin(qJ(5));
t310 = cos(qJ(5));
t375 = qJD(2) * t300;
t387 = t311 * t312;
t228 = (-t305 * t387 + t307 * t310) * t375;
t368 = qJD(5) * t310;
t503 = t305 * t370 + t306 * t368 + t228;
t346 = pkin(4) * t306 - pkin(11) * t311;
t277 = t346 * qJD(4);
t284 = -pkin(4) * t311 - pkin(11) * t306 - pkin(3);
t369 = qJD(5) * t305;
t182 = t305 * t277 + t284 * t368 + (-t310 * t371 - t311 * t369) * pkin(10);
t268 = pkin(9) * t375 + t362;
t259 = t307 * t268;
t280 = qJD(2) * pkin(2) + t313 * t377;
t303 = cos(pkin(6));
t376 = qJD(1) * t303;
t326 = t280 * t302 + t300 * t376;
t185 = t312 * t326 - t259;
t253 = t347 * t375;
t139 = t311 * t185 + t306 * t253;
t361 = t307 * t375;
t127 = pkin(11) * t361 + t139;
t397 = t268 * t312;
t332 = t280 * t393 + t397;
t374 = qJD(2) * t312;
t141 = (t307 * t376 + t346 * t374) * t300 + t332;
t72 = t310 * t127 + t305 * t141;
t502 = t182 - t72;
t247 = t293 + (-pkin(2) * t312 - pkin(3)) * t302;
t172 = pkin(4) * t261 - pkin(11) * t262 + t247;
t188 = t311 * t248 + t306 * t249;
t174 = -pkin(11) * t395 + t188;
t86 = t305 * t172 + t310 * t174;
t476 = -qJD(5) * t86 - t305 * t505 + t504 * t310;
t475 = t172 * t368 - t174 * t369 + t504 * t305 + t310 * t505;
t229 = (t305 * t307 + t310 * t387) * t375;
t388 = t310 * t311;
t297 = pkin(10) * t388;
t360 = t300 * t374;
t432 = pkin(10) * t305;
t378 = t310 * t277 + t371 * t432;
t430 = pkin(12) * t306;
t433 = pkin(5) * t306;
t71 = -t127 * t305 + t310 * t141;
t501 = pkin(12) * t229 - t360 * t433 + (-pkin(12) * t388 + t433) * qJD(4) + (-t297 + (-t284 + t430) * t305) * qJD(5) + t378 - t71;
t500 = -pkin(12) * t503 + t502;
t215 = -t262 * t305 - t310 * t395;
t132 = qJD(5) * t215 + t213 * t310 + t305 * t359;
t499 = pkin(5) * t214 - pkin(12) * t132 + t476;
t333 = -t262 * t310 + t305 * t395;
t133 = qJD(5) * t333 - t213 * t305 + t310 * t359;
t498 = -pkin(12) * t133 - t475;
t459 = -pkin(12) - pkin(11);
t363 = qJD(5) * t459;
t292 = qJD(2) * t302 + qJD(3);
t237 = t292 * t311 - t306 * t361;
t398 = t237 * t305;
t238 = t292 * t306 + t311 * t361;
t184 = pkin(4) * t238 - pkin(11) * t237;
t186 = t307 * t326 + t397;
t167 = pkin(10) * t292 + t186;
t291 = t302 * t376;
t194 = t291 + (qJD(2) * t348 - t280) * t300;
t98 = -t306 * t167 + t194 * t311;
t70 = t305 * t184 + t310 * t98;
t497 = pkin(12) * t398 + t305 * t363 - t70;
t429 = pkin(12) * t310;
t69 = t310 * t184 - t305 * t98;
t496 = -pkin(5) * t238 + t237 * t429 + t310 * t363 - t69;
t166 = -pkin(3) * t292 - t185;
t101 = -pkin(4) * t237 - pkin(11) * t238 + t166;
t285 = qJD(4) - t360;
t99 = t167 * t311 + t194 * t306;
t88 = pkin(11) * t285 + t99;
t48 = t310 * t101 - t305 * t88;
t49 = t101 * t305 + t310 * t88;
t335 = t305 * t49 + t310 * t48;
t337 = Ifges(6,5) * t310 - Ifges(6,6) * t305;
t424 = Ifges(6,4) * t310;
t339 = -Ifges(6,2) * t305 + t424;
t425 = Ifges(6,4) * t305;
t341 = Ifges(6,1) * t310 - t425;
t342 = mrSges(6,1) * t305 + mrSges(6,2) * t310;
t435 = t310 / 0.2e1;
t436 = -t305 / 0.2e1;
t232 = qJD(5) - t237;
t440 = t232 / 0.2e1;
t197 = t238 * t310 + t285 * t305;
t447 = t197 / 0.2e1;
t196 = -t238 * t305 + t285 * t310;
t449 = t196 / 0.2e1;
t87 = -pkin(4) * t285 - t98;
t426 = Ifges(6,4) * t197;
t94 = Ifges(6,2) * t196 + Ifges(6,6) * t232 + t426;
t195 = Ifges(6,4) * t196;
t95 = Ifges(6,1) * t197 + Ifges(6,5) * t232 + t195;
t495 = t335 * mrSges(6,3) - t337 * t440 - t339 * t449 - t341 * t447 - t342 * t87 - t435 * t95 - t436 * t94;
t304 = sin(qJ(6));
t309 = cos(qJ(6));
t355 = t309 * t196 - t197 * t304;
t114 = Ifges(7,4) * t355;
t118 = t196 * t304 + t197 * t309;
t40 = -pkin(12) * t197 + t48;
t35 = pkin(5) * t232 + t40;
t41 = pkin(12) * t196 + t49;
t403 = t304 * t41;
t16 = t309 * t35 - t403;
t402 = t309 * t41;
t17 = t304 * t35 + t402;
t423 = Ifges(7,4) * t118;
t222 = qJD(6) + t232;
t444 = -t222 / 0.2e1;
t452 = -t118 / 0.2e1;
t454 = -t355 / 0.2e1;
t59 = Ifges(7,1) * t118 + Ifges(7,5) * t222 + t114;
t75 = -pkin(5) * t196 + t87;
t494 = (Ifges(7,5) * t355 - Ifges(7,6) * t118) * t444 + (t118 * t17 + t16 * t355) * mrSges(7,3) + (-Ifges(7,2) * t118 + t114 + t59) * t454 - t75 * (mrSges(7,1) * t118 + mrSges(7,2) * t355) + (Ifges(7,1) * t355 - t423) * t452;
t202 = t292 * t370 + (-t307 * t371 + t357) * t375;
t356 = qJD(2) * t373;
t350 = t307 * t356;
t108 = qJD(5) * t196 + t202 * t310 + t305 * t350;
t322 = t329 * qJD(2);
t364 = t303 * t395;
t351 = qJD(3) * t364;
t129 = (t280 * t392 - t259) * qJD(3) + (t301 * t322 + t351) * qJD(1);
t207 = (t328 + t362) * t375;
t52 = t311 * t129 - t167 * t371 + t194 * t370 + t306 * t207;
t43 = pkin(11) * t350 + t52;
t323 = t331 * qJD(2);
t352 = t303 * t359;
t130 = t332 * qJD(3) + (t301 * t323 + t352) * qJD(1);
t203 = t292 * t371 + (t307 * t370 + t358) * t375;
t78 = pkin(4) * t203 - pkin(11) * t202 + t130;
t11 = -qJD(5) * t49 - t305 * t43 + t310 * t78;
t8 = pkin(5) * t203 - pkin(12) * t108 + t11;
t10 = t101 * t368 + t305 * t78 + t310 * t43 - t369 * t88;
t109 = -qJD(5) * t197 - t202 * t305 + t310 * t350;
t9 = pkin(12) * t109 + t10;
t2 = qJD(6) * t16 + t304 * t8 + t309 * t9;
t3 = -qJD(6) * t17 - t304 * t9 + t309 * t8;
t493 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t272 = t310 * t284;
t208 = -t306 * t429 + t272 + (-pkin(5) - t432) * t311;
t245 = t305 * t284 + t297;
t217 = -t305 * t430 + t245;
t142 = t208 * t309 - t217 * t304;
t492 = qJD(6) * t142 + t304 * t501 + t309 * t500;
t143 = t208 * t304 + t217 * t309;
t491 = -qJD(6) * t143 - t304 * t500 + t309 * t501;
t138 = -t306 * t185 + t253 * t311;
t126 = -pkin(4) * t361 - t138;
t490 = pkin(5) * t503 + pkin(10) * t370 - t126;
t231 = Ifges(5,4) * t237;
t422 = Ifges(5,5) * t285;
t428 = Ifges(5,1) * t238;
t157 = t231 + t422 + t428;
t321 = t98 * mrSges(5,3) - t157 / 0.2e1 - t422 / 0.2e1 - t166 * mrSges(5,2);
t489 = t321 + t495;
t38 = qJD(6) * t355 + t108 * t309 + t109 * t304;
t463 = t38 / 0.2e1;
t39 = -qJD(6) * t118 - t108 * t304 + t109 * t309;
t462 = t39 / 0.2e1;
t58 = Ifges(7,2) * t355 + Ifges(7,6) * t222 + t423;
t487 = t58 / 0.2e1;
t446 = t203 / 0.2e1;
t486 = t361 / 0.2e1;
t485 = -t292 * Ifges(4,6) / 0.2e1;
t85 = t310 * t172 - t174 * t305;
t68 = pkin(5) * t261 + pkin(12) * t333 + t85;
t79 = pkin(12) * t215 + t86;
t30 = t304 * t68 + t309 * t79;
t484 = -qJD(6) * t30 + t498 * t304 + t309 * t499;
t29 = -t304 * t79 + t309 * t68;
t483 = qJD(6) * t29 + t304 * t499 - t498 * t309;
t482 = -t11 * mrSges(6,1) + t10 * mrSges(6,2) - t493;
t36 = Ifges(7,6) * t39;
t37 = Ifges(7,5) * t38;
t12 = Ifges(7,3) * t203 + t36 + t37;
t104 = Ifges(6,6) * t109;
t105 = Ifges(6,5) * t108;
t45 = Ifges(6,3) * t203 + t104 + t105;
t481 = t45 + t12;
t287 = t459 * t305;
t288 = t459 * t310;
t220 = t287 * t309 + t288 * t304;
t474 = qJD(6) * t220 + t304 * t496 + t309 * t497;
t221 = t287 * t304 - t288 * t309;
t473 = -qJD(6) * t221 - t304 * t497 + t309 * t496;
t67 = -mrSges(7,1) * t355 + mrSges(7,2) * t118;
t472 = m(7) * t75 + t67;
t334 = t304 * t305 - t309 * t310;
t258 = t334 * t306;
t380 = -mrSges(4,1) * t292 - mrSges(5,1) * t237 + mrSges(5,2) * t238 + mrSges(4,3) * t361;
t469 = t302 * t386 - t391;
t468 = t10 * t310 - t11 * t305;
t467 = qJD(5) + qJD(6);
t53 = -t306 * t129 - t167 * t370 - t194 * t371 + t207 * t311;
t466 = -t53 * mrSges(5,1) + t52 * mrSges(5,2) - Ifges(5,5) * t202 + Ifges(5,6) * t203;
t465 = Ifges(7,4) * t463 + Ifges(7,2) * t462 + Ifges(7,6) * t446;
t464 = Ifges(7,1) * t463 + Ifges(7,4) * t462 + Ifges(7,5) * t446;
t47 = t108 * Ifges(6,1) + t109 * Ifges(6,4) + t203 * Ifges(6,5);
t461 = t47 / 0.2e1;
t460 = -t94 / 0.2e1;
t456 = t108 / 0.2e1;
t455 = t109 / 0.2e1;
t453 = t355 / 0.2e1;
t451 = t118 / 0.2e1;
t450 = -t196 / 0.2e1;
t448 = -t197 / 0.2e1;
t443 = t222 / 0.2e1;
t442 = -t231 / 0.2e1;
t441 = -t232 / 0.2e1;
t439 = -t261 / 0.2e1;
t437 = t262 / 0.2e1;
t434 = pkin(5) * t305;
t431 = pkin(10) * t311;
t427 = Ifges(5,4) * t238;
t419 = t355 * Ifges(7,6);
t418 = t118 * Ifges(7,5);
t417 = t196 * Ifges(6,6);
t416 = t197 * Ifges(6,5);
t415 = t202 * Ifges(5,1);
t414 = t202 * Ifges(5,4);
t413 = t203 * Ifges(5,4);
t412 = t222 * Ifges(7,3);
t411 = t232 * Ifges(6,3);
t410 = t237 * Ifges(5,2);
t407 = t285 * Ifges(5,6);
t179 = mrSges(5,1) * t350 - mrSges(5,3) * t202;
t63 = -mrSges(6,1) * t109 + mrSges(6,2) * t108;
t401 = t63 - t179;
t209 = -t301 * t469 - t364;
t400 = t130 * t209;
t394 = t301 * t308;
t123 = -mrSges(6,1) * t196 + mrSges(6,2) * t197;
t206 = mrSges(5,1) * t285 - mrSges(5,3) * t238;
t385 = t123 - t206;
t274 = t304 * t310 + t305 * t309;
t212 = t467 * t274;
t149 = -t212 * t306 - t334 * t370;
t169 = t228 * t304 + t229 * t309;
t384 = t149 - t169;
t150 = t258 * t467 - t274 * t370;
t168 = t228 * t309 - t229 * t304;
t383 = t150 - t168;
t164 = t274 * t237;
t382 = -t164 + t212;
t165 = t334 * t237;
t211 = t467 * t334;
t381 = -t165 + t211;
t366 = t67 + t385;
t187 = -t306 * t248 + t249 * t311;
t353 = t375 * t394;
t173 = pkin(4) * t395 - t187;
t345 = -t130 * mrSges(4,1) - t129 * mrSges(4,2);
t344 = -mrSges(4,1) * t312 + mrSges(4,2) * t307;
t343 = mrSges(6,1) * t310 - mrSges(6,2) * t305;
t340 = Ifges(6,1) * t305 + t424;
t338 = Ifges(6,2) * t310 + t425;
t336 = Ifges(6,5) * t305 + Ifges(6,6) * t310;
t330 = t302 * t390 + t389;
t210 = t301 * t330 + t303 * t396;
t260 = -t300 * t301 * t313 + t302 * t303;
t171 = t210 * t311 + t260 * t306;
t106 = -t171 * t305 + t209 * t310;
t107 = t171 * t310 + t209 * t305;
t61 = t106 * t309 - t107 * t304;
t62 = t106 * t304 + t107 * t309;
t170 = t210 * t306 - t260 * t311;
t146 = t215 * t309 + t304 * t333;
t147 = t215 * t304 - t309 * t333;
t111 = -t248 * t370 - t249 * t371 + t254 * t311 - t306 * t255;
t230 = -t280 * t300 + t291;
t289 = Ifges(4,4) * t360;
t320 = -t185 * mrSges(4,3) + Ifges(4,1) * t486 + t289 / 0.2e1 + t292 * Ifges(4,5) + t230 * mrSges(4,2);
t103 = -pkin(4) * t359 - t111;
t44 = -pkin(4) * t350 - t53;
t319 = t230 * mrSges(4,1) + t98 * mrSges(5,1) + t285 * Ifges(5,3) + t238 * Ifges(5,5) + t237 * Ifges(5,6) + t485 - (Ifges(4,4) * t307 + Ifges(4,2) * t312) * t375 / 0.2e1 - t186 * mrSges(4,3) - t99 * mrSges(5,2);
t156 = t407 + t410 + t427;
t57 = t412 + t418 + t419;
t93 = t411 + t416 + t417;
t317 = t407 / 0.2e1 - t417 / 0.2e1 - t416 / 0.2e1 - t411 / 0.2e1 - t166 * mrSges(5,1) + t99 * mrSges(5,3) + t49 * mrSges(6,2) - t48 * mrSges(6,1) + t17 * mrSges(7,2) - t16 * mrSges(7,1) - t419 / 0.2e1 - t418 / 0.2e1 - t93 / 0.2e1 - t57 / 0.2e1 - t412 / 0.2e1 + t156 / 0.2e1 + t427 / 0.2e1;
t315 = t317 + t410 / 0.2e1;
t299 = -pkin(5) * t310 - pkin(4);
t283 = Ifges(4,5) * t312 * t356;
t282 = Ifges(5,3) * t350;
t279 = (pkin(10) + t434) * t306;
t257 = t274 * t306;
t252 = t344 * t375;
t251 = -mrSges(4,2) * t292 + mrSges(4,3) * t360;
t244 = -t305 * t431 + t272;
t243 = (mrSges(4,1) * t307 + mrSges(4,2) * t312) * t356;
t205 = -mrSges(5,2) * t285 + mrSges(5,3) * t237;
t192 = t224 * t306 - t311 * t354;
t183 = -qJD(5) * t245 + t378;
t180 = -mrSges(5,2) * t350 - mrSges(5,3) * t203;
t160 = t351 + (qJD(3) * t469 + t322) * t301;
t159 = t352 + (qJD(3) * t330 + t323) * t301;
t145 = mrSges(6,1) * t232 - mrSges(6,3) * t197;
t144 = -mrSges(6,2) * t232 + mrSges(6,3) * t196;
t131 = mrSges(5,1) * t203 + mrSges(5,2) * t202;
t119 = -pkin(5) * t215 + t173;
t113 = Ifges(5,5) * t350 - t413 + t415;
t112 = -t203 * Ifges(5,2) + Ifges(5,6) * t350 + t414;
t90 = mrSges(7,1) * t222 - mrSges(7,3) * t118;
t89 = -mrSges(7,2) * t222 + mrSges(7,3) * t355;
t84 = pkin(5) * t398 + t99;
t83 = -mrSges(6,2) * t203 + mrSges(6,3) * t109;
t82 = mrSges(6,1) * t203 - mrSges(6,3) * t108;
t81 = -qJD(4) * t170 + t160 * t311 + t306 * t353;
t80 = qJD(4) * t171 + t160 * t306 - t311 * t353;
t66 = -pkin(5) * t133 + t103;
t51 = -qJD(6) * t147 - t132 * t304 + t133 * t309;
t50 = qJD(6) * t146 + t132 * t309 + t133 * t304;
t46 = t108 * Ifges(6,4) + t109 * Ifges(6,2) + t203 * Ifges(6,6);
t34 = qJD(5) * t106 + t159 * t305 + t310 * t81;
t33 = -qJD(5) * t107 + t159 * t310 - t305 * t81;
t28 = -mrSges(7,2) * t203 + mrSges(7,3) * t39;
t27 = mrSges(7,1) * t203 - mrSges(7,3) * t38;
t26 = -pkin(5) * t109 + t44;
t19 = t309 * t40 - t403;
t18 = -t304 * t40 - t402;
t15 = -mrSges(7,1) * t39 + mrSges(7,2) * t38;
t7 = -qJD(6) * t62 - t304 * t34 + t309 * t33;
t6 = qJD(6) * t61 + t304 * t33 + t309 * t34;
t1 = [t106 * t82 + t107 * t83 + t209 * t131 + t34 * t144 + t33 * t145 + t160 * t251 + t171 * t180 + t81 * t205 + t260 * t243 + t61 * t27 + t62 * t28 + t6 * t89 + t7 * t90 + (-mrSges(3,1) * t308 - mrSges(3,2) * t313) * qJD(2) ^ 2 * t301 + t380 * t159 + t366 * t80 + (t15 + t401) * t170 + (t252 * t394 + (t209 * t312 - t210 * t307) * qJD(3) * mrSges(4,3)) * t375 + m(4) * (t129 * t210 + t400 - t159 * t185 + t160 * t186 + (qJD(1) * t260 + t230) * t353) + m(5) * (t159 * t166 - t170 * t53 + t171 * t52 - t80 * t98 + t81 * t99 + t400) + m(6) * (t10 * t107 + t106 * t11 + t170 * t44 + t33 * t48 + t34 * t49 + t80 * t87) + m(7) * (t16 * t7 + t17 * t6 + t170 * t26 + t2 * t62 + t3 * t61 + t75 * t80); t187 * t179 + t188 * t180 + t173 * t63 + m(4) * (t129 * t266 - t130 * t265 - t185 * t256 + t186 * t255) + (t93 + t57) * t214 / 0.2e1 + (-t213 * t98 - t214 * t99 - t261 * t52 - t262 * t53) * mrSges(5,3) + t113 * t437 + t112 * t439 + (Ifges(6,5) * t132 + Ifges(6,6) * t133 + Ifges(6,3) * t214) * t440 + (Ifges(7,5) * t50 + Ifges(7,6) * t51 + Ifges(7,3) * t214) * t443 + t237 * (Ifges(5,4) * t213 - Ifges(5,2) * t214) / 0.2e1 - t203 * (Ifges(5,4) * t262 - Ifges(5,2) * t261) / 0.2e1 + (-Ifges(6,5) * t333 + Ifges(7,5) * t147 + Ifges(6,6) * t215 + Ifges(7,6) * t146 + (Ifges(6,3) + Ifges(7,3)) * t261) * t446 + t44 * (-mrSges(6,1) * t215 - mrSges(6,2) * t333) + t11 * (mrSges(6,1) * t261 + mrSges(6,3) * t333) + (-Ifges(6,4) * t333 + Ifges(6,2) * t215 + Ifges(6,6) * t261) * t455 + (-Ifges(6,1) * t333 + Ifges(6,4) * t215 + Ifges(6,5) * t261) * t456 - t333 * t461 + ((-m(4) * pkin(2) + t344) * t354 + ((Ifges(4,5) * t302 / 0.2e1 - t265 * mrSges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,4) * t395) * t312 + (-t266 * mrSges(4,3) - Ifges(4,6) * t302 + Ifges(5,5) * t437 + Ifges(5,6) * t439 - 0.3e1 / 0.2e1 * Ifges(4,4) * t396 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(5,3) / 0.2e1) * t395) * t307) * qJD(3)) * t375 + t26 * (-mrSges(7,1) * t146 + mrSges(7,2) * t147) + t132 * t95 / 0.2e1 + t87 * (-mrSges(6,1) * t133 + mrSges(6,2) * t132) + t133 * t94 / 0.2e1 + t119 * t15 + t103 * t123 + t85 * t82 + t86 * t83 + t75 * (-mrSges(7,1) * t51 + mrSges(7,2) * t50) + t66 * t67 + t50 * t59 / 0.2e1 - t366 * t192 + t29 * t27 + t30 * t28 - m(4) * (-t185 * t223 + t186 * t224 + t230 * t354) + t238 * (Ifges(5,1) * t213 - Ifges(5,4) * t214) / 0.2e1 + t202 * (Ifges(5,1) * t262 - Ifges(5,4) * t261) / 0.2e1 + (t283 / 0.2e1 + t345) * t302 + (Ifges(6,1) * t132 + Ifges(6,4) * t133 + Ifges(6,5) * t214) * t447 + (Ifges(6,4) * t132 + Ifges(6,2) * t133 + Ifges(6,6) * t214) * t449 + (Ifges(7,1) * t50 + Ifges(7,4) * t51 + Ifges(7,5) * t214) * t451 + (Ifges(7,4) * t50 + Ifges(7,2) * t51 + Ifges(7,6) * t214) * t453 + t111 * t206 + t213 * t157 / 0.2e1 + t48 * (mrSges(6,1) * t214 - mrSges(6,3) * t132) + t49 * (-mrSges(6,2) * t214 + mrSges(6,3) * t133) + t16 * (mrSges(7,1) * t214 - mrSges(7,3) * t50) + t17 * (-mrSges(7,2) * t214 + mrSges(7,3) * t51) - t214 * t156 / 0.2e1 + t166 * (mrSges(5,1) * t214 + mrSges(5,2) * t213) + t215 * t46 / 0.2e1 + t247 * t131 + t2 * (-mrSges(7,2) * t261 + mrSges(7,3) * t146) + t3 * (mrSges(7,1) * t261 - mrSges(7,3) * t147) + t10 * (-mrSges(6,2) * t261 + mrSges(6,3) * t215) + t130 * (mrSges(5,1) * t261 + mrSges(5,2) * t262) + (Ifges(7,4) * t147 + Ifges(7,2) * t146 + Ifges(7,6) * t261) * t462 + (Ifges(7,1) * t147 + Ifges(7,4) * t146 + Ifges(7,5) * t261) * t463 + t147 * t464 + t146 * t465 + t506 * t251 + t471 * t205 + (t130 * t247 + t187 * t53 + t188 * t52 + t471 * t99 + (t111 + t192) * t98 - t470 * t166) * m(5) - t380 * t470 + t285 * (Ifges(5,5) * t213 - Ifges(5,6) * t214) / 0.2e1 + t475 * t144 + t476 * t145 + (t10 * t86 + t11 * t85 + t173 * t44 + (t103 - t192) * t87 + t475 * t49 + t476 * t48) * m(6) + t481 * t261 / 0.2e1 + t483 * t89 + t484 * t90 + (t119 * t26 + t2 * t30 + t29 * t3 + (-t192 + t66) * t75 + t483 * t17 + t484 * t16) * m(7) + (-t252 * t362 + t130 * mrSges(4,3) * t307 - pkin(2) * t243 + (-t282 / 0.2e1 + t129 * mrSges(4,3) + t466) * t312 + (t320 * t312 + (t485 + t319) * t307) * qJD(3)) * t300 + t51 * t487; t491 * t90 + (t142 * t3 + t143 * t2 + t16 * t491 + t17 * t492 + t26 * t279 + t490 * t75) * m(7) + t492 * t89 + (-t169 / 0.2e1 + t149 / 0.2e1) * t59 - m(5) * (t138 * t98 + t139 * t99 + t166 * t186) - m(6) * (t126 * t87 + t48 * t71 + t49 * t72) + (t113 / 0.2e1 + t130 * mrSges(5,2) + t415 / 0.2e1 - t413 / 0.2e1 + t47 * t435 - t53 * mrSges(5,3) + t46 * t436 + t44 * t342 + t341 * t456 + t339 * t455 + t337 * t446 + (-t10 * t305 - t11 * t310) * mrSges(6,3) + (-m(5) * t53 + m(6) * t44 + t401) * pkin(10) + (t95 * t436 + t310 * t460 + t338 * t450 + t340 * t448 + t336 * t441 + t87 * t343 + (t305 * t48 - t310 * t49) * mrSges(6,3)) * qJD(5)) * t306 + (t183 - t71) * t145 + (-Ifges(7,5) * t258 - Ifges(7,6) * t257) * t446 + (-t16 * t384 + t17 * t383 - t2 * t257 + t258 * t3) * mrSges(7,3) + t26 * (mrSges(7,1) * t257 - mrSges(7,2) * t258) + (-Ifges(7,4) * t258 - Ifges(7,2) * t257) * t462 + (-Ifges(7,1) * t258 - Ifges(7,4) * t257) * t463 + (-t168 / 0.2e1 + t150 / 0.2e1) * t58 + t490 * t67 + (Ifges(6,5) * t229 + Ifges(6,6) * t228) * t441 + (Ifges(7,5) * t149 + Ifges(7,6) * t150) * t443 + (Ifges(7,5) * t169 + Ifges(7,6) * t168) * t444 + t283 + t345 + (-t228 * t49 + t229 * t48) * mrSges(6,3) + m(5) * (-pkin(3) * t130 + t431 * t52) + t142 * t27 + t143 * t28 - t126 * t123 - pkin(3) * t131 + m(6) * (t10 * t245 + t11 * t244 + t182 * t49 + t183 * t48) + (-mrSges(7,1) * t383 + mrSges(7,2) * t384) * t75 - t380 * t186 + t502 * t144 + (Ifges(6,1) * t229 + Ifges(6,4) * t228) * t448 + (Ifges(6,4) * t229 + Ifges(6,2) * t228) * t450 + (Ifges(7,1) * t149 + Ifges(7,4) * t150) * t451 + (Ifges(7,1) * t169 + Ifges(7,4) * t168) * t452 + (Ifges(7,4) * t149 + Ifges(7,2) * t150) * t453 - t139 * t205 - t138 * t206 + ((-t315 + (-m(5) * t99 - t205) * pkin(10)) * t306 + (t231 / 0.2e1 + t428 / 0.2e1 + (-m(5) * t98 + m(6) * t87 + t385) * pkin(10) - t489) * t311) * qJD(4) - t229 * t95 / 0.2e1 - t87 * (-mrSges(6,1) * t228 + mrSges(6,2) * t229) + t244 * t82 + t245 * t83 - t185 * t251 + (Ifges(7,4) * t169 + Ifges(7,2) * t168) * t454 + t228 * t460 - t258 * t464 - t257 * t465 + t279 * t15 + (-t130 * mrSges(5,1) - t105 / 0.2e1 - t104 / 0.2e1 + pkin(10) * t180 + t52 * mrSges(5,3) + t414 / 0.2e1 + t112 / 0.2e1 - t45 / 0.2e1 - t12 / 0.2e1 - t37 / 0.2e1 - t36 / 0.2e1 + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(7,3) / 0.2e1) * t203 + t482) * t311 + ((qJD(3) * (Ifges(5,5) * t306 + Ifges(5,6) * t311) / 0.2e1 + Ifges(4,4) * t486 + (t292 / 0.2e1 - qJD(3)) * Ifges(4,6) - t319) * t307 + (-t289 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t361 + (t442 - t428 / 0.2e1 + t321) * t311 + t315 * t306 - t320) * t312) * t375; (-pkin(4) * t44 - t48 * t69 - t49 * t70 - t87 * t99) * m(6) - t466 + (-Ifges(7,1) * t165 - Ifges(7,4) * t164) * t452 + (-Ifges(7,4) * t165 - Ifges(7,2) * t164) * t454 + (-Ifges(7,5) * t165 - Ifges(7,6) * t164) * t444 + t46 * t435 + t282 + (t16 * t381 - t17 * t382 - t2 * t334 - t274 * t3) * mrSges(7,3) + (Ifges(7,5) * t274 - Ifges(7,6) * t334 + t336) * t446 + t26 * (mrSges(7,1) * t334 + mrSges(7,2) * t274) + (Ifges(7,4) * t274 - Ifges(7,2) * t334) * t462 + (Ifges(7,1) * t274 - Ifges(7,4) * t334) * t463 - t334 * t465 + (t165 / 0.2e1 - t211 / 0.2e1) * t59 + (t164 / 0.2e1 - t212 / 0.2e1) * t58 + (-Ifges(7,5) * t211 - Ifges(7,6) * t212) * t443 + (-Ifges(7,1) * t211 - Ifges(7,4) * t212) * t451 + (-Ifges(7,4) * t211 - Ifges(7,2) * t212) * t453 - t69 * t145 - t70 * t144 - t84 * t67 - pkin(4) * t63 + (mrSges(7,1) * t382 - mrSges(7,2) * t381) * t75 - t385 * t99 - t44 * t343 - t98 * t205 + t220 * t27 + t221 * t28 + (t442 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t238 + t489) * t237 + (t472 * t434 - t495) * qJD(5) + t338 * t455 + t340 * t456 + t305 * t461 + t274 * t464 + ((-m(6) * t335 - t305 * t144 - t310 * t145) * qJD(5) + m(6) * t468 - t305 * t82 + t310 * t83) * pkin(11) + t468 * mrSges(6,3) + t299 * t15 + t473 * t90 + t474 * t89 + (t16 * t473 + t17 * t474 + t2 * t221 + t220 * t3 + t26 * t299 - t75 * t84) * m(7) + t317 * t238; -t87 * (mrSges(6,1) * t197 + mrSges(6,2) * t196) + t481 - t482 + t118 * t487 + (Ifges(6,5) * t196 - Ifges(6,6) * t197) * t441 + (t196 * t48 + t197 * t49) * mrSges(6,3) + t49 * t145 - t48 * t144 - t19 * t89 - t18 * t90 + (-Ifges(6,2) * t197 + t195 + t95) * t450 + t94 * t447 + (Ifges(6,1) * t196 - t426) * t448 - m(7) * (t16 * t18 + t17 * t19) + (t309 * t27 + t304 * t28 + m(7) * (t2 * t304 + t3 * t309) - t472 * t197 + (-t304 * t90 + t309 * t89 + m(7) * (-t16 * t304 + t17 * t309)) * qJD(6)) * pkin(5) + t494; -t16 * t89 + t17 * t90 + t58 * t451 + t12 + t493 + t494;];
tauc  = t1(:);
