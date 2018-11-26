% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:22:38
% EndTime: 2018-11-23 17:23:01
% DurationCPUTime: 22.67s
% Computational Cost: add. (26463->828), mult. (77570->1182), div. (0->0), fcn. (63153->12), ass. (0->377)
t305 = sin(qJ(2));
t301 = cos(pkin(6));
t435 = pkin(1) * t301;
t293 = t305 * t435;
t300 = sin(pkin(6));
t309 = cos(qJ(2));
t391 = t300 * t309;
t428 = pkin(8) + qJ(3);
t488 = t391 * t428 + t293;
t241 = t488 * qJD(1);
t399 = sin(pkin(12));
t230 = t399 * t241;
t294 = t309 * t435;
t289 = qJD(1) * t294;
t365 = t428 * t305;
t352 = t300 * t365;
t240 = -qJD(1) * t352 + t289;
t400 = cos(pkin(12));
t190 = t240 * t400 - t230;
t385 = qJD(1) * t300;
t347 = t399 * t385;
t348 = t400 * t385;
t256 = -t305 * t347 + t309 * t348;
t257 = -t305 * t348 - t309 * t347;
t369 = t305 * t385;
t356 = pkin(2) * t369;
t207 = -pkin(3) * t257 - pkin(9) * t256 + t356;
t304 = sin(qJ(4));
t308 = cos(qJ(4));
t126 = -t190 * t304 + t308 * t207;
t370 = t399 * pkin(2);
t295 = t370 + pkin(9);
t427 = pkin(10) + t295;
t361 = qJD(4) * t427;
t395 = t256 * t308;
t517 = pkin(4) * t257 + pkin(10) * t395 - t308 * t361 - t126;
t127 = t308 * t190 + t304 * t207;
t396 = t256 * t304;
t516 = -pkin(10) * t396 + t304 * t361 + t127;
t382 = qJD(4) * t304;
t515 = t382 - t396;
t303 = sin(qJ(5));
t307 = cos(qJ(5));
t279 = t303 * t308 + t304 * t307;
t205 = t279 * t256;
t490 = qJD(4) + qJD(5);
t234 = t490 * t279;
t495 = t205 - t234;
t278 = t303 * t304 - t307 * t308;
t206 = t278 * t256;
t233 = t490 * t278;
t514 = t206 - t233;
t276 = t427 * t304;
t277 = t427 * t308;
t328 = -t307 * t276 - t277 * t303;
t500 = qJD(5) * t328 + t303 * t517 - t516 * t307;
t360 = t400 * t241;
t189 = t240 * t399 + t360;
t498 = pkin(4) * t515 - t189;
t513 = pkin(11) * t257 + t500;
t512 = -pkin(5) * t495 - t514 * pkin(11) + t498;
t254 = qJD(4) - t256;
t250 = qJD(5) + t254;
t302 = sin(qJ(6));
t306 = cos(qJ(6));
t291 = qJD(1) * t301 + qJD(2);
t219 = t257 * t304 + t291 * t308;
t329 = t257 * t308 - t291 * t304;
t331 = t219 * t303 - t307 * t329;
t132 = t250 * t306 - t302 * t331;
t264 = (t305 * t400 + t309 * t399) * t300;
t258 = qJD(2) * t264;
t251 = qJD(1) * t258;
t383 = qJD(2) * t300;
t493 = -t399 * t305 + t400 * t309;
t259 = t493 * t383;
t252 = qJD(1) * t259;
t174 = qJD(4) * t219 + t252 * t308;
t175 = qJD(4) * t329 - t252 * t304;
t357 = t307 * t219 + t303 * t329;
t90 = qJD(5) * t357 + t174 * t307 + t175 * t303;
t54 = qJD(6) * t132 + t251 * t302 + t306 * t90;
t133 = t250 * t302 + t306 * t331;
t55 = -qJD(6) * t133 + t251 * t306 - t302 * t90;
t21 = -mrSges(7,1) * t55 + mrSges(7,2) * t54;
t223 = pkin(2) * t291 + t240;
t169 = t399 * t223 + t360;
t163 = pkin(9) * t291 + t169;
t280 = (-pkin(2) * t309 - pkin(1)) * t300;
t273 = qJD(1) * t280 + qJD(3);
t185 = -t256 * pkin(3) + t257 * pkin(9) + t273;
t121 = t163 * t308 + t185 * t304;
t105 = pkin(10) * t219 + t121;
t389 = t307 * t105;
t120 = -t163 * t304 + t308 * t185;
t104 = pkin(10) * t329 + t120;
t93 = pkin(4) * t254 + t104;
t44 = t303 * t93 + t389;
t286 = qJD(2) * t289;
t317 = (-qJD(2) * t365 + qJD(3) * t309) * t300;
t212 = qJD(1) * t317 + t286;
t392 = t300 * t305;
t225 = -qJD(2) * t488 - qJD(3) * t392;
t314 = qJD(1) * t225;
t147 = t212 * t400 + t314 * t399;
t364 = qJD(1) * t383;
t351 = t305 * t364;
t335 = pkin(2) * t351;
t186 = pkin(3) * t251 - pkin(9) * t252 + t335;
t76 = -t121 * qJD(4) - t147 * t304 + t308 * t186;
t45 = pkin(4) * t251 - pkin(10) * t174 + t76;
t381 = qJD(4) * t308;
t75 = t308 * t147 - t163 * t382 + t185 * t381 + t304 * t186;
t51 = pkin(10) * t175 + t75;
t11 = -qJD(5) * t44 - t303 * t51 + t307 * t45;
t8 = -pkin(5) * t251 - t11;
t511 = m(7) * t8 + t21;
t168 = t223 * t400 - t230;
t162 = -t291 * pkin(3) - t168;
t130 = -t219 * pkin(4) + t162;
t41 = pkin(11) * t250 + t44;
t77 = -pkin(5) * t357 - pkin(11) * t331 + t130;
t23 = t302 * t77 + t306 * t41;
t507 = t23 * mrSges(7,2);
t22 = -t302 * t41 + t306 * t77;
t508 = t22 * mrSges(7,1);
t510 = -t130 * mrSges(6,1) + t44 * mrSges(6,3) + t507 - t508;
t222 = -t276 * t303 + t277 * t307;
t501 = -qJD(5) * t222 + t516 * t303 + t307 * t517;
t446 = t250 / 0.2e1;
t456 = t331 / 0.2e1;
t457 = t357 / 0.2e1;
t154 = qJD(6) - t357;
t458 = t154 / 0.2e1;
t463 = t133 / 0.2e1;
t465 = t132 / 0.2e1;
t414 = t154 * Ifges(7,3);
t415 = t133 * Ifges(7,5);
t416 = t132 * Ifges(7,6);
t72 = t414 + t415 + t416;
t406 = t250 * Ifges(6,6);
t413 = t357 * Ifges(6,2);
t421 = Ifges(6,4) * t331;
t99 = t406 + t413 + t421;
t509 = -Ifges(6,4) * t456 + Ifges(7,5) * t463 - Ifges(6,2) * t457 - Ifges(6,6) * t446 + Ifges(7,6) * t465 + Ifges(7,3) * t458 - t99 / 0.2e1 + t72 / 0.2e1 - t510;
t80 = mrSges(6,1) * t251 - mrSges(6,3) * t90;
t506 = t21 - t80;
t371 = t400 * pkin(2);
t296 = -t371 - pkin(3);
t284 = -t308 * pkin(4) + t296;
t216 = t278 * pkin(5) - t279 * pkin(11) + t284;
t152 = t216 * t306 - t222 * t302;
t505 = qJD(6) * t152 + t512 * t302 + t513 * t306;
t153 = t216 * t302 + t222 * t306;
t504 = -qJD(6) * t153 - t513 * t302 + t512 * t306;
t503 = -Ifges(5,5) * t329 + Ifges(6,5) * t331 + t219 * Ifges(5,6) + Ifges(6,6) * t357 + t254 * Ifges(5,3) + t250 * Ifges(6,3);
t502 = -pkin(5) * t257 - t501;
t135 = mrSges(6,1) * t250 - mrSges(6,3) * t331;
t86 = -mrSges(7,1) * t132 + mrSges(7,2) * t133;
t401 = t86 - t135;
t237 = pkin(2) * t301 + t294 - t352;
t386 = pkin(8) * t391 + t293;
t255 = qJ(3) * t391 + t386;
t200 = t399 * t237 + t400 * t255;
t188 = pkin(9) * t301 + t200;
t263 = t493 * t300;
t209 = -t263 * pkin(3) - t264 * pkin(9) + t280;
t128 = -t188 * t304 + t308 * t209;
t236 = t264 * t308 + t301 * t304;
t109 = -pkin(4) * t263 - pkin(10) * t236 + t128;
t129 = t308 * t188 + t304 * t209;
t235 = -t264 * t304 + t301 * t308;
t118 = pkin(10) * t235 + t129;
t499 = t303 * t109 + t307 * t118;
t149 = t206 * t302 - t257 * t306;
t377 = qJD(6) * t306;
t325 = -t302 * t233 + t279 * t377;
t497 = t149 + t325;
t150 = -t206 * t306 - t257 * t302;
t378 = qJD(6) * t302;
t324 = t306 * t233 + t279 * t378;
t496 = t150 + t324;
t425 = mrSges(4,3) * t257;
t388 = -mrSges(4,1) * t291 - mrSges(5,1) * t219 - mrSges(5,2) * t329 - t425;
t151 = Ifges(6,4) * t357;
t340 = Ifges(7,5) * t306 - Ifges(7,6) * t302;
t321 = t154 * t340;
t419 = Ifges(7,4) * t302;
t344 = Ifges(7,1) * t306 - t419;
t322 = t133 * t344;
t418 = Ifges(7,4) * t306;
t342 = -Ifges(7,2) * t302 + t418;
t323 = t132 * t342;
t346 = mrSges(7,1) * t302 + mrSges(7,2) * t306;
t398 = t105 * t303;
t43 = t307 * t93 - t398;
t40 = -pkin(5) * t250 - t43;
t326 = t40 * t346;
t131 = Ifges(7,4) * t132;
t74 = Ifges(7,1) * t133 + Ifges(7,5) * t154 + t131;
t403 = t306 * t74;
t407 = t250 * Ifges(6,5);
t437 = t302 / 0.2e1;
t412 = t331 * Ifges(6,1);
t100 = t151 + t407 + t412;
t468 = -t100 / 0.2e1;
t489 = -t130 * mrSges(6,2) + t43 * mrSges(6,3);
t420 = Ifges(7,4) * t133;
t73 = Ifges(7,2) * t132 + Ifges(7,6) * t154 + t420;
t494 = t468 - t323 / 0.2e1 - t322 / 0.2e1 - t321 / 0.2e1 - t407 / 0.2e1 + t73 * t437 - t403 / 0.2e1 - t326 + t489 - t151 / 0.2e1;
t492 = -t304 * t76 + t308 * t75;
t491 = -t22 * t302 + t23 * t306;
t197 = qJD(4) * t235 + t259 * t308;
t290 = qJD(2) * t294;
t224 = t290 + t317;
t166 = t224 * t400 + t225 * t399;
t366 = t305 * t383;
t355 = pkin(2) * t366;
t208 = pkin(3) * t258 - pkin(9) * t259 + t355;
t79 = -qJD(4) * t129 - t166 * t304 + t308 * t208;
t61 = pkin(4) * t258 - pkin(10) * t197 + t79;
t198 = -qJD(4) * t236 - t259 * t304;
t78 = t308 * t166 - t188 * t382 + t304 * t208 + t209 * t381;
t67 = pkin(10) * t198 + t78;
t15 = -qJD(5) * t499 - t303 * t67 + t307 * t61;
t146 = t212 * t399 - t400 * t314;
t116 = -t175 * pkin(4) + t146;
t91 = qJD(5) * t331 + t174 * t303 - t307 * t175;
t24 = t91 * pkin(5) - t90 * pkin(11) + t116;
t379 = qJD(5) * t307;
t380 = qJD(5) * t303;
t10 = -t105 * t380 + t303 * t45 + t307 * t51 + t93 * t379;
t7 = pkin(11) * t251 + t10;
t2 = qJD(6) * t22 + t24 * t302 + t306 * t7;
t3 = -qJD(6) * t23 + t24 * t306 - t302 * t7;
t487 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t114 = pkin(5) * t331 - pkin(11) * t357;
t18 = Ifges(7,5) * t54 + Ifges(7,6) * t55 + Ifges(7,3) * t91;
t471 = t91 / 0.2e1;
t472 = -t91 / 0.2e1;
t473 = t90 / 0.2e1;
t477 = t55 / 0.2e1;
t478 = t54 / 0.2e1;
t485 = t116 * mrSges(6,1) - t10 * mrSges(6,3) - Ifges(6,4) * t473 + Ifges(7,5) * t478 - Ifges(6,2) * t472 + Ifges(7,6) * t477 + Ifges(7,3) * t471 + t18 / 0.2e1 + t487;
t483 = -0.2e1 * pkin(1);
t19 = t54 * Ifges(7,4) + t55 * Ifges(7,2) + t91 * Ifges(7,6);
t480 = t19 / 0.2e1;
t20 = t54 * Ifges(7,1) + t55 * Ifges(7,4) + t91 * Ifges(7,5);
t479 = t20 / 0.2e1;
t476 = -t72 / 0.2e1;
t469 = t99 / 0.2e1;
t467 = t100 / 0.2e1;
t466 = -t132 / 0.2e1;
t464 = -t133 / 0.2e1;
t409 = t329 * Ifges(5,4);
t137 = t219 * Ifges(5,2) + t254 * Ifges(5,6) - t409;
t462 = t137 / 0.2e1;
t217 = Ifges(5,4) * t219;
t138 = -Ifges(5,1) * t329 + t254 * Ifges(5,5) + t217;
t461 = t138 / 0.2e1;
t459 = -t154 / 0.2e1;
t455 = t174 / 0.2e1;
t454 = t175 / 0.2e1;
t330 = t307 * t235 - t236 * t303;
t453 = t330 / 0.2e1;
t180 = t235 * t303 + t236 * t307;
t452 = t180 / 0.2e1;
t451 = -t219 / 0.2e1;
t450 = t329 / 0.2e1;
t449 = -t329 / 0.2e1;
t448 = t235 / 0.2e1;
t447 = t236 / 0.2e1;
t445 = t251 / 0.2e1;
t444 = -t254 / 0.2e1;
t442 = -t257 / 0.2e1;
t438 = t301 / 0.2e1;
t436 = -t305 / 0.2e1;
t434 = t2 * t306;
t433 = t3 * t302;
t426 = mrSges(4,3) * t256;
t424 = Ifges(3,4) * t305;
t423 = Ifges(5,4) * t304;
t422 = Ifges(5,4) * t308;
t417 = Ifges(3,5) * t309;
t410 = t22 * t306;
t405 = t257 * Ifges(4,4);
t394 = t279 * t302;
t393 = t279 * t306;
t375 = Ifges(6,5) * t90 - Ifges(6,6) * t91 + Ifges(6,3) * t251;
t373 = t403 / 0.2e1;
t372 = Ifges(5,5) * t174 + Ifges(5,6) * t175 + Ifges(5,3) * t251;
t368 = t309 * t385;
t362 = -t378 / 0.2e1;
t354 = mrSges(3,3) * t369;
t353 = mrSges(3,3) * t368;
t345 = Ifges(5,1) * t308 - t423;
t343 = -Ifges(5,2) * t304 + t422;
t341 = Ifges(5,5) * t308 - Ifges(5,6) * t304;
t339 = t23 * t302 + t410;
t25 = mrSges(7,1) * t91 - mrSges(7,3) * t54;
t26 = -mrSges(7,2) * t91 + mrSges(7,3) * t55;
t338 = -t302 * t25 + t306 * t26;
t57 = -pkin(11) * t263 + t499;
t199 = t237 * t400 - t399 * t255;
t187 = -t301 * pkin(3) - t199;
t142 = -t235 * pkin(4) + t187;
t85 = -pkin(5) * t330 - t180 * pkin(11) + t142;
t28 = t302 * t85 + t306 * t57;
t27 = -t302 * t57 + t306 * t85;
t165 = t224 * t399 - t400 * t225;
t62 = t109 * t307 - t118 * t303;
t144 = t180 * t306 - t263 * t302;
t143 = -t180 * t302 - t263 * t306;
t183 = -mrSges(5,2) * t254 + mrSges(5,3) * t219;
t184 = mrSges(5,1) * t254 + mrSges(5,3) * t329;
t332 = t183 * t308 - t184 * t304;
t134 = -mrSges(6,2) * t250 + mrSges(6,3) * t357;
t96 = -mrSges(7,2) * t154 + mrSges(7,3) * t132;
t97 = mrSges(7,1) * t154 - mrSges(7,3) * t133;
t327 = -t302 * t97 + t306 * t96 + t134;
t319 = t291 * (-Ifges(3,6) * t305 + t417);
t14 = t109 * t379 - t118 * t380 + t303 * t61 + t307 * t67;
t123 = -t198 * pkin(4) + t165;
t270 = t386 * qJD(2);
t316 = -qJD(6) * t339 - t433;
t260 = -pkin(8) * t351 + t286;
t261 = qJD(1) * t270;
t315 = -t261 * mrSges(3,1) - t146 * mrSges(4,1) - t260 * mrSges(3,2) - t147 * mrSges(4,2);
t313 = -t10 * mrSges(6,2) + mrSges(7,3) * t434 + t20 * t437 + t306 * t480 + t375 + (Ifges(7,1) * t302 + t418) * t478 + (Ifges(7,2) * t306 + t419) * t477 + t8 * (-mrSges(7,1) * t306 + mrSges(7,2) * t302) + t73 * t362 + (Ifges(7,5) * t302 + Ifges(7,6) * t306) * t471 + t11 * mrSges(6,1) + (t326 + t373) * qJD(6) + (t323 + t322 + t321) * qJD(6) / 0.2e1;
t312 = m(7) * (-t22 * t377 - t23 * t378 - t433 + t434) - t96 * t378 - t97 * t377 + t338;
t311 = t476 + t469 + t406 / 0.2e1 - t414 / 0.2e1 - t416 / 0.2e1 - t415 / 0.2e1 + t421 / 0.2e1 + t510;
t288 = Ifges(3,4) * t368;
t285 = t364 * t417;
t274 = -pkin(8) * t392 + t294;
t269 = -pkin(8) * t366 + t290;
t268 = t386 * qJD(1);
t267 = -pkin(8) * t369 + t289;
t266 = -t291 * mrSges(3,2) + t353;
t265 = mrSges(3,1) * t291 - t354;
t253 = Ifges(4,4) * t256;
t246 = Ifges(4,5) * t252;
t245 = Ifges(4,6) * t251;
t242 = t252 * mrSges(4,2);
t239 = Ifges(3,1) * t369 + t291 * Ifges(3,5) + t288;
t238 = Ifges(3,6) * t291 + (Ifges(3,2) * t309 + t424) * t385;
t226 = -mrSges(4,2) * t291 + t426;
t210 = -mrSges(4,1) * t256 - mrSges(4,2) * t257;
t204 = -t257 * Ifges(4,1) + t291 * Ifges(4,5) + t253;
t203 = t256 * Ifges(4,2) + t291 * Ifges(4,6) - t405;
t140 = -mrSges(5,2) * t251 + mrSges(5,3) * t175;
t139 = mrSges(5,1) * t251 - mrSges(5,3) * t174;
t122 = -mrSges(5,1) * t175 + mrSges(5,2) * t174;
t113 = -mrSges(6,1) * t357 + mrSges(6,2) * t331;
t112 = t174 * Ifges(5,1) + t175 * Ifges(5,4) + t251 * Ifges(5,5);
t111 = t174 * Ifges(5,4) + t175 * Ifges(5,2) + t251 * Ifges(5,6);
t103 = qJD(5) * t180 + t197 * t303 - t307 * t198;
t102 = qJD(5) * t330 + t197 * t307 + t198 * t303;
t95 = -pkin(4) * t329 + t114;
t81 = -mrSges(6,2) * t251 - mrSges(6,3) * t91;
t71 = -qJD(6) * t144 - t102 * t302 + t258 * t306;
t70 = qJD(6) * t143 + t102 * t306 + t258 * t302;
t56 = pkin(5) * t263 - t62;
t47 = t104 * t307 - t398;
t46 = t104 * t303 + t389;
t38 = mrSges(6,1) * t91 + mrSges(6,2) * t90;
t37 = t90 * Ifges(6,1) - t91 * Ifges(6,4) + t251 * Ifges(6,5);
t36 = t90 * Ifges(6,4) - t91 * Ifges(6,2) + t251 * Ifges(6,6);
t35 = t103 * pkin(5) - t102 * pkin(11) + t123;
t32 = t114 * t302 + t306 * t43;
t31 = t114 * t306 - t302 * t43;
t30 = t302 * t95 + t306 * t47;
t29 = -t302 * t47 + t306 * t95;
t13 = -pkin(5) * t258 - t15;
t12 = pkin(11) * t258 + t14;
t5 = -qJD(6) * t28 - t12 * t302 + t306 * t35;
t4 = qJD(6) * t27 + t12 * t306 + t302 * t35;
t1 = [m(3) * (t260 * t386 - t261 * t274 - t267 * t270 + t268 * t269) + (-t200 * mrSges(4,3) + Ifges(6,5) * t452 + Ifges(6,6) * t453 + Ifges(5,5) * t447 + Ifges(5,6) * t448 - Ifges(4,4) * t264 + t280 * mrSges(4,1) - Ifges(4,6) * t301 / 0.2e1 - (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t263) * t251 + (-t199 * mrSges(4,3) + Ifges(4,1) * t264 + Ifges(4,4) * t263 + Ifges(4,5) * t438) * t252 + (t10 * t263 + t102 * t130 + t116 * t180 - t258 * t44) * mrSges(6,2) + (t146 * t264 + t147 * t263 - t168 * t259 - t169 * t258) * mrSges(4,3) + (Ifges(6,4) * t180 - Ifges(6,6) * t263) * t472 + ((Ifges(3,5) * t438 - t274 * mrSges(3,3) + (mrSges(3,2) * t483 + 0.3e1 / 0.2e1 * Ifges(3,4) * t309) * t300) * t309 + (-t386 * mrSges(3,3) - Ifges(3,6) * t301 + (mrSges(3,1) * t483 - 0.3e1 / 0.2e1 * t424 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t309) * t300 + (m(4) * t280 - mrSges(4,1) * t263 + mrSges(4,2) * t264) * pkin(2)) * t305) * t364 + t11 * (-mrSges(6,1) * t263 - mrSges(6,3) * t180) + t75 * (mrSges(5,2) * t263 + mrSges(5,3) * t235) + t76 * (-mrSges(5,1) * t263 - mrSges(5,3) * t236) + (Ifges(5,4) * t236 + Ifges(5,2) * t235 - Ifges(5,6) * t263) * t454 + (Ifges(5,1) * t236 + Ifges(5,4) * t235 - Ifges(5,5) * t263) * t455 + (Ifges(6,1) * t180 - Ifges(6,5) * t263) * t473 + t388 * t165 + (Ifges(7,4) * t70 + Ifges(7,2) * t71) * t465 + (Ifges(7,4) * t144 + Ifges(7,2) * t143) * t477 + (t260 * t309 + t261 * t305) * t300 * mrSges(3,3) + t509 * t103 + (t305 * pkin(2) * t210 + t309 * t239 / 0.2e1 + t238 * t436 + t319 / 0.2e1 + (-t267 * t309 - t268 * t305) * mrSges(3,3)) * t383 + t123 * t113 + (Ifges(7,5) * t70 + Ifges(7,6) * t71) * t458 + (Ifges(7,5) * t144 + Ifges(7,6) * t143) * t471 + m(4) * (-t146 * t199 + t147 * t200 - t165 * t168 + t166 * t169 + t273 * t355) + m(6) * (t10 * t499 + t11 * t62 + t116 * t142 + t123 * t130 + t14 * t44 + t15 * t43) + t499 * t81 - t485 * t330 + (t143 * t2 - t144 * t3 - t22 * t70 + t23 * t71) * mrSges(7,3) + t503 * t258 / 0.2e1 + t280 * t242 + (Ifges(6,5) * t102 + Ifges(6,3) * t258) * t446 + (Ifges(6,4) * t102 + Ifges(6,6) * t258) * t457 + t129 * t140 + t142 * t38 + t8 * (-mrSges(7,1) * t143 + mrSges(7,2) * t144) + t14 * t134 + t15 * t135 + t128 * t139 + m(5) * (t120 * t79 + t121 * t78 + t128 * t76 + t129 * t75 + t146 * t187 + t162 * t165) + m(7) * (t13 * t40 + t2 * t28 + t22 * t5 + t23 * t4 + t27 * t3 + t56 * t8) + t4 * t96 + t5 * t97 + t13 * t86 + t62 * t80 + t71 * t73 / 0.2e1 + t70 * t74 / 0.2e1 + t40 * (-mrSges(7,1) * t71 + mrSges(7,2) * t70) + t56 * t21 + t27 * t25 + t28 * t26 - (t375 + t372) * t263 / 0.2e1 + (Ifges(7,1) * t70 + Ifges(7,4) * t71) * t463 + (Ifges(7,1) * t144 + Ifges(7,4) * t143) * t478 + t78 * t183 + t79 * t184 + t187 * t122 + t162 * (-mrSges(5,1) * t198 + mrSges(5,2) * t197) + t166 * t226 + t146 * (-mrSges(5,1) * t235 + mrSges(5,2) * t236) + t43 * (mrSges(6,1) * t258 - mrSges(6,3) * t102) + t120 * (mrSges(5,1) * t258 - mrSges(5,3) * t197) + t121 * (-mrSges(5,2) * t258 + mrSges(5,3) * t198) + t219 * (Ifges(5,4) * t197 + Ifges(5,2) * t198 + Ifges(5,6) * t258) / 0.2e1 + t254 * (Ifges(5,5) * t197 + Ifges(5,6) * t198 + Ifges(5,3) * t258) / 0.2e1 - t258 * t203 / 0.2e1 + t259 * t204 / 0.2e1 + t256 * (Ifges(4,4) * t259 - Ifges(4,2) * t258) / 0.2e1 + t269 * t266 - t270 * t265 + t273 * (mrSges(4,1) * t258 + mrSges(4,2) * t259) + t291 * (Ifges(4,5) * t259 - Ifges(4,6) * t258) / 0.2e1 + (Ifges(4,1) * t259 - Ifges(4,4) * t258) * t442 + t112 * t447 + t111 * t448 + (Ifges(5,1) * t197 + Ifges(5,4) * t198 + Ifges(5,5) * t258) * t449 + t37 * t452 + t36 * t453 + t197 * t461 + t198 * t462 + t102 * t467 + t144 * t479 + t143 * t480 + (Ifges(6,1) * t102 + Ifges(6,5) * t258) * t456 + (t285 / 0.2e1 + t246 / 0.2e1 - t245 / 0.2e1 + t315) * t301; (-t515 * t121 + (-t381 + t395) * t120 + t492) * mrSges(5,3) + (-t120 * t126 - t121 * t127 + t146 * t296 - t162 * t189) * m(5) + t205 * t507 + (pkin(1) * (mrSges(3,1) * t305 + mrSges(3,2) * t309) + (Ifges(3,1) * t309 - t424) * t436) * qJD(1) ^ 2 * t300 ^ 2 + (-t36 / 0.2e1 - Ifges(6,6) * t445 + t485) * t278 + t254 * t162 * (mrSges(5,1) * t304 + mrSges(5,2) * t308) - t137 * t382 / 0.2e1 + (-t11 * mrSges(6,3) + t8 * t346 + t74 * t362 + t37 / 0.2e1 + t116 * mrSges(6,2) + Ifges(6,5) * t445 + t340 * t471 + Ifges(6,4) * t472 + Ifges(6,1) * t473 + t342 * t477 + t344 * t478) * t279 + t285 + t509 * t234 - t121 * mrSges(5,2) * t257 + t120 * mrSges(5,1) * t257 + (-Ifges(7,4) * t324 - Ifges(7,2) * t325) * t465 + (m(5) * ((-t120 * t308 - t121 * t304) * qJD(4) + t492) - t304 * t139 + t308 * t140 - t184 * t381 - t183 * t382) * t295 - t388 * t189 + (-Ifges(7,5) * t324 - Ifges(7,6) * t325) * t458 - t19 * t394 / 0.2e1 + (-Ifges(7,1) * t324 - Ifges(7,4) * t325) * t463 - t245 + t246 - t138 * t395 / 0.2e1 - t210 * t356 - ((-Ifges(3,2) * t369 + t239 + t288) * t309 + t319) * t385 / 0.2e1 + (t219 * t343 + t254 * t341 - t329 * t345) * qJD(4) / 0.2e1 + (t10 * t222 + t11 * t328 + t116 * t284 + t130 * t498 + t43 * t501 + t44 * t500) * m(6) + (t152 * t3 + t153 * t2 + t22 * t504 + t23 * t505 - t328 * t8 + t40 * t502) * m(7) - t506 * t328 + t315 + t500 * t134 + t501 * t135 + t502 * t86 + (Ifges(4,1) * t256 + t405 + t503) * t257 / 0.2e1 - t169 * t425 + t152 * t25 + t153 * t26 - t150 * t74 / 0.2e1 + (t354 + t265) * t268 + (-t251 * t370 - t252 * t371) * mrSges(4,3) + (t353 - t266) * t267 - Ifges(3,6) * t351 - t331 * (-Ifges(6,1) * t206 - Ifges(6,4) * t205 - Ifges(6,5) * t257) / 0.2e1 - t357 * (-Ifges(6,4) * t206 - Ifges(6,2) * t205 - Ifges(6,6) * t257) / 0.2e1 - t130 * (mrSges(6,1) * t205 - mrSges(6,2) * t206) - t43 * (-mrSges(6,1) * t257 + mrSges(6,3) * t206) - t250 * (-Ifges(6,5) * t206 - Ifges(6,6) * t205 - Ifges(6,3) * t257) / 0.2e1 - (Ifges(4,2) * t257 + t204 + t253) * t256 / 0.2e1 + (t168 * t189 - t169 * t190 - t273 * t356 + (-t146 * t400 + t147 * t399) * pkin(2)) * m(4) - t497 * t73 / 0.2e1 + (mrSges(7,1) * t497 - mrSges(7,2) * t496) * t40 + (-t2 * t394 + t22 * t496 - t23 * t497 - t3 * t393) * mrSges(7,3) + t498 * t113 + t504 * t97 + t505 * t96 - t205 * t508 - t127 * t183 - t126 * t184 + t168 * t426 + t222 * t81 - t190 * t226 - t44 * (mrSges(6,2) * t257 - mrSges(6,3) * t205) + t238 * t369 / 0.2e1 - t273 * (-mrSges(4,1) * t257 + mrSges(4,2) * t256) + t284 * t38 - t291 * (Ifges(4,5) * t256 + Ifges(4,6) * t257) / 0.2e1 + t296 * t122 + t203 * t442 + (-Ifges(5,3) * t257 + t256 * t341) * t444 + (Ifges(5,5) * t304 + Ifges(5,6) * t308) * t445 + (-Ifges(5,5) * t257 + t256 * t345) * t450 + (-Ifges(5,6) * t257 + t256 * t343) * t451 + (Ifges(5,2) * t308 + t423) * t454 + (Ifges(5,1) * t304 + t422) * t455 + (Ifges(7,5) * t150 + Ifges(7,6) * t149 + Ifges(7,3) * t205) * t459 + t381 * t461 + t396 * t462 + (Ifges(7,1) * t150 + Ifges(7,4) * t149 + Ifges(7,5) * t205) * t464 + (Ifges(7,4) * t150 + Ifges(7,2) * t149 + Ifges(7,6) * t205) * t466 - t206 * t468 + t205 * t469 + t205 * t476 + t393 * t479 - (Ifges(6,1) * t456 + Ifges(6,4) * t457 + Ifges(6,5) * t446 + t373 + t467 - t489) * t233 + t304 * t112 / 0.2e1 + t308 * t111 / 0.2e1 + t146 * (-mrSges(5,1) * t308 + mrSges(5,2) * t304); t251 * mrSges(4,1) + t206 * t134 + t308 * t139 + t304 * t140 - t149 * t97 - t150 * t96 + t242 + t506 * t278 + t332 * qJD(4) + (t113 + t388) * t257 + (-t226 - t332) * t256 - t327 * t233 + (t81 + (-t302 * t96 - t306 * t97) * qJD(6) + t338) * t279 - t401 * t495 + (-t149 * t22 - t150 * t23 + t278 * t8 - t491 * t233 + (t316 + t434) * t279 - t495 * t40) * m(7) + (t10 * t279 - t11 * t278 + t130 * t257 + t495 * t43 + t44 * t514) * m(6) + (t162 * t257 + t304 * t75 + t308 * t76 + t254 * (-t120 * t304 + t121 * t308)) * m(5) + (-t168 * t257 - t169 * t256 + t335) * m(4); (-t154 * t410 + (-t154 * t23 - t3) * t302) * mrSges(7,3) + (-t412 / 0.2e1 + t494) * t357 + t372 + (t120 * t219 - t121 * t329) * mrSges(5,3) - t162 * (-mrSges(5,1) * t329 + mrSges(5,2) * t219) + (Ifges(5,5) * t219 + Ifges(5,6) * t329) * t444 + (t329 * t113 + t303 * t81 + t307 * t80 + ((m(7) * t40 + t401) * t303 + (m(7) * t491 + t327) * t307) * qJD(5) + (0.2e1 * t130 * t450 + t10 * t303 + t11 * t307 + (-t303 * t43 + t307 * t44) * qJD(5)) * m(6)) * pkin(4) - t47 * t134 - t30 * t96 - t29 * t97 - t75 * mrSges(5,2) + t76 * mrSges(5,1) + (Ifges(5,1) * t219 + t409) * t450 + (t413 / 0.2e1 + t311) * t331 + t313 - t120 * t183 + t121 * t184 - m(6) * (-t43 * t46 + t44 * t47) - t401 * t46 - m(7) * (t22 * t29 + t23 * t30 + t40 * t46) + t137 * t449 + t312 * (pkin(4) * t303 + pkin(11)) + (Ifges(5,2) * t329 + t138 + t217) * t451 + t511 * (-pkin(4) * t307 - pkin(5)); t316 * mrSges(7,3) - m(7) * (t22 * t31 + t23 * t32 + t40 * t44) - t43 * t134 - t401 * t44 - t32 * t96 - t31 * t97 + t313 + ((Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t331 + t339 * mrSges(7,3) + t494) * t357 + t311 * t331 + t312 * pkin(11) - t511 * pkin(5); -t40 * (mrSges(7,1) * t133 + mrSges(7,2) * t132) + (Ifges(7,1) * t132 - t420) * t464 + t73 * t463 + (Ifges(7,5) * t132 - Ifges(7,6) * t133) * t459 - t22 * t96 + t23 * t97 + (t132 * t22 + t133 * t23) * mrSges(7,3) + t18 + (-Ifges(7,2) * t133 + t131 + t74) * t466 + t487;];
tauc  = t1(:);
