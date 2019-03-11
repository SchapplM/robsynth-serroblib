% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:29:01
% EndTime: 2019-03-09 10:30:00
% DurationCPUTime: 31.89s
% Computational Cost: add. (19581->880), mult. (59203->1258), div. (0->0), fcn. (47923->12), ass. (0->364)
t309 = sin(pkin(6));
t308 = sin(pkin(11));
t314 = sin(qJ(2));
t317 = cos(qJ(2));
t392 = cos(pkin(11));
t348 = t392 * t317;
t323 = -t308 * t314 + t348;
t320 = t309 * t323;
t255 = qJD(1) * t320;
t466 = -t255 + qJD(4);
t311 = cos(pkin(6));
t426 = pkin(1) * t311;
t302 = t317 * t426;
t298 = qJD(1) * t302;
t422 = pkin(8) + qJ(3);
t356 = t422 * t314;
t342 = t309 * t356;
t242 = -qJD(1) * t342 + t298;
t301 = t314 * t426;
t383 = t309 * t317;
t465 = t383 * t422 + t301;
t243 = t465 * qJD(1);
t350 = t392 * t243;
t184 = t242 * t308 + t350;
t313 = sin(qJ(4));
t316 = cos(qJ(4));
t492 = -qJD(5) * t313 - t184 + t466 * (pkin(4) * t313 - qJ(5) * t316);
t231 = t308 * t243;
t185 = t242 * t392 - t231;
t349 = t392 * t314;
t373 = qJD(1) * t309;
t361 = t317 * t373;
t256 = -t308 * t361 - t349 * t373;
t362 = t314 * t373;
t346 = pkin(2) * t362;
t199 = -pkin(3) * t256 - pkin(9) * t255 + t346;
t116 = t316 * t185 + t313 * t199;
t101 = -qJ(5) * t256 + t116;
t307 = sin(pkin(12));
t310 = cos(pkin(12));
t491 = -t310 * t101 + t307 * t492;
t425 = pkin(2) * t308;
t304 = pkin(9) + t425;
t369 = qJD(4) * t313;
t359 = t304 * t369;
t476 = t492 * t310 + (t101 + t359) * t307;
t381 = t310 * t316;
t202 = t255 * t381 - t256 * t307;
t389 = t255 * t313;
t490 = -pkin(5) * t389 + pkin(10) * t202 + (pkin(5) * t313 - pkin(10) * t381) * qJD(4) + t476;
t385 = t307 * t316;
t201 = -t255 * t385 - t256 * t310;
t382 = t310 * t313;
t489 = -pkin(10) * t201 + (-pkin(10) * t385 - t304 * t382) * qJD(4) + t491;
t345 = qJD(1) * t311 + qJD(2);
t219 = -t313 * t256 - t316 * t345;
t218 = qJD(6) + t219;
t402 = t218 * Ifges(7,3);
t220 = -t316 * t256 + t313 * t345;
t253 = -t323 * t373 + qJD(4);
t169 = t220 * t310 + t253 * t307;
t312 = sin(qJ(6));
t315 = cos(qJ(6));
t347 = -t220 * t307 + t310 * t253;
t110 = t169 * t315 + t312 * t347;
t406 = t110 * Ifges(7,5);
t483 = -t169 * t312 + t315 * t347;
t407 = t483 * Ifges(7,6);
t44 = t402 + t406 + t407;
t403 = t169 * Ifges(6,5);
t404 = t347 * Ifges(6,6);
t85 = t219 * Ifges(6,3) + t403 + t404;
t480 = t85 + t44;
t488 = t480 / 0.2e1;
t258 = qJD(2) * t320;
t251 = qJD(1) * t258;
t370 = qJD(4) * t219;
t170 = t316 * t251 - t370;
t265 = (t308 * t317 + t349) * t309;
t257 = qJD(2) * t265;
t250 = qJD(1) * t257;
t133 = -t170 * t307 + t250 * t310;
t134 = t170 * t310 + t250 * t307;
t171 = qJD(4) * t220 + t313 * t251;
t38 = qJD(6) * t483 + t133 * t312 + t134 * t315;
t39 = -qJD(6) * t110 + t133 * t315 - t134 * t312;
t8 = Ifges(7,5) * t38 + Ifges(7,6) * t39 + Ifges(7,3) * t171;
t487 = t134 * Ifges(6,5) + t133 * Ifges(6,6) + t171 * Ifges(6,3) + t8;
t486 = Ifges(6,3) + Ifges(7,3);
t485 = Ifges(4,5) * t251;
t355 = qJD(2) * t373;
t412 = Ifges(3,5) * t317;
t484 = -Ifges(4,6) * t250 + t355 * t412 + t485;
t460 = t38 / 0.2e1;
t459 = t39 / 0.2e1;
t447 = t133 / 0.2e1;
t446 = t134 / 0.2e1;
t441 = t171 / 0.2e1;
t363 = t392 * pkin(2);
t305 = -t363 - pkin(3);
t283 = -t316 * pkin(4) - t313 * qJ(5) + t305;
t269 = t310 * t283;
t210 = -pkin(10) * t382 + t269 + (-t304 * t307 - pkin(5)) * t316;
t229 = t307 * t283 + t304 * t381;
t386 = t307 * t313;
t215 = -pkin(10) * t386 + t229;
t139 = t210 * t315 - t215 * t312;
t482 = qJD(6) * t139 + t490 * t312 + t489 * t315;
t140 = t210 * t312 + t215 * t315;
t481 = -qJD(6) * t140 - t489 * t312 + t490 * t315;
t135 = mrSges(5,1) * t250 - mrSges(5,3) * t170;
t77 = -t133 * mrSges(6,1) + t134 * mrSges(6,2);
t479 = -t135 + t77;
t421 = pkin(10) + qJ(5);
t292 = t421 * t307;
t293 = t421 * t310;
t237 = -t292 * t312 + t293 * t315;
t286 = t307 * t315 + t310 * t312;
t155 = pkin(4) * t220 + qJ(5) * t219;
t321 = t311 * pkin(2) - t342;
t221 = qJD(2) * pkin(2) + qJD(1) * t321 + t298;
t160 = t308 * t221 + t350;
t154 = pkin(9) * t345 + t160;
t287 = (-pkin(2) * t317 - pkin(1)) * t309;
t374 = qJD(1) * t287;
t278 = qJD(3) + t374;
t180 = -t255 * pkin(3) + t256 * pkin(9) + t278;
t99 = -t313 * t154 + t180 * t316;
t71 = t310 * t155 - t307 * t99;
t49 = pkin(10) * t219 * t310 + pkin(5) * t220 + t71;
t390 = t219 * t307;
t72 = t307 * t155 + t310 * t99;
t59 = pkin(10) * t390 + t72;
t478 = -qJD(5) * t286 - qJD(6) * t237 + t312 * t59 - t315 * t49;
t236 = -t292 * t315 - t293 * t312;
t325 = t307 * t312 - t310 * t315;
t477 = -qJD(5) * t325 + qJD(6) * t236 - t312 * t49 - t315 * t59;
t475 = -t310 * t359 + t491;
t115 = -t313 * t185 + t199 * t316;
t102 = pkin(4) * t256 - t115;
t351 = pkin(5) * t307 + t304;
t368 = qJD(4) * t316;
t474 = pkin(5) * t201 + t351 * t368 - t102;
t124 = t201 * t315 - t202 * t312;
t279 = t325 * qJD(6);
t212 = t279 * t313 - t286 * t368;
t471 = t124 - t212;
t125 = t201 * t312 + t202 * t315;
t280 = t286 * qJD(6);
t211 = -t280 * t313 - t325 * t368;
t470 = t125 - t211;
t398 = t256 * mrSges(4,3);
t377 = -mrSges(4,1) * t345 + mrSges(5,1) * t219 + mrSges(5,2) * t220 - t398;
t295 = qJD(2) * t298;
t319 = (-qJD(2) * t356 + qJD(3) * t317) * t309;
t213 = qJD(1) * t319 + t295;
t384 = t309 * t314;
t225 = -qJD(2) * t465 - qJD(3) * t384;
t318 = qJD(1) * t225;
t142 = t213 * t392 + t308 * t318;
t341 = t314 * t355;
t328 = pkin(2) * t341;
t181 = pkin(3) * t250 - pkin(9) * t251 + t328;
t52 = -t313 * t142 - t154 * t368 - t180 * t369 + t181 * t316;
t42 = -pkin(4) * t250 - t52;
t88 = -pkin(4) * t253 + qJD(5) - t99;
t469 = t313 * t42 + t88 * t368;
t51 = t316 * t142 - t154 * t369 + t180 * t368 + t313 * t181;
t468 = -t313 * t52 + t316 * t51;
t86 = t169 * Ifges(6,4) + Ifges(6,2) * t347 + Ifges(6,6) * t219;
t455 = -t86 / 0.2e1;
t467 = -t99 * mrSges(5,3) + t307 * t455;
t100 = t154 * t316 + t180 * t313;
t89 = qJ(5) * t253 + t100;
t159 = t221 * t392 - t231;
t153 = -pkin(3) * t345 - t159;
t96 = t219 * pkin(4) - t220 * qJ(5) + t153;
t40 = -t307 * t89 + t310 * t96;
t25 = pkin(5) * t219 - pkin(10) * t169 + t40;
t41 = t307 * t96 + t310 * t89;
t29 = pkin(10) * t347 + t41;
t5 = t25 * t315 - t29 * t312;
t6 = t25 * t312 + t29 * t315;
t464 = -t5 * mrSges(7,1) + t6 * mrSges(7,2) + t100 * mrSges(5,3);
t463 = -0.2e1 * pkin(1);
t462 = Ifges(7,4) * t460 + Ifges(7,2) * t459 + Ifges(7,6) * t441;
t461 = Ifges(7,1) * t460 + Ifges(7,4) * t459 + Ifges(7,5) * t441;
t413 = Ifges(7,4) * t110;
t45 = Ifges(7,2) * t483 + Ifges(7,6) * t218 + t413;
t458 = t45 / 0.2e1;
t104 = Ifges(7,4) * t483;
t46 = Ifges(7,1) * t110 + Ifges(7,5) * t218 + t104;
t457 = t46 / 0.2e1;
t456 = Ifges(6,1) * t446 + Ifges(6,4) * t447 + Ifges(6,5) * t441;
t454 = -t483 / 0.2e1;
t453 = t483 / 0.2e1;
t452 = -t110 / 0.2e1;
t451 = t110 / 0.2e1;
t400 = t253 * Ifges(5,6);
t418 = Ifges(5,4) * t220;
t130 = -t219 * Ifges(5,2) + t400 + t418;
t450 = -t130 / 0.2e1;
t449 = t130 / 0.2e1;
t216 = Ifges(5,4) * t219;
t401 = t253 * Ifges(5,5);
t131 = t220 * Ifges(5,1) - t216 + t401;
t448 = t131 / 0.2e1;
t445 = t347 / 0.2e1;
t444 = t169 / 0.2e1;
t443 = t170 / 0.2e1;
t442 = -t171 / 0.2e1;
t439 = -t218 / 0.2e1;
t438 = t218 / 0.2e1;
t437 = -t219 / 0.2e1;
t436 = t219 / 0.2e1;
t234 = t265 * t313 - t311 * t316;
t435 = -t234 / 0.2e1;
t235 = t265 * t316 + t311 * t313;
t433 = t235 / 0.2e1;
t431 = -t256 / 0.2e1;
t429 = t310 / 0.2e1;
t33 = qJ(5) * t250 + qJD(5) * t253 + t51;
t141 = t213 * t308 - t392 * t318;
t69 = pkin(4) * t171 - qJ(5) * t170 - qJD(5) * t220 + t141;
t19 = t307 * t69 + t310 * t33;
t264 = t308 * t384 - t309 * t348;
t299 = qJD(2) * t302;
t224 = t299 + t319;
t158 = t224 * t392 + t308 * t225;
t238 = t302 + t321;
t375 = pkin(8) * t383 + t301;
t254 = qJ(3) * t383 + t375;
t193 = t308 * t238 + t392 * t254;
t183 = pkin(9) * t311 + t193;
t371 = qJD(2) * t314;
t360 = t309 * t371;
t200 = pkin(2) * t360 + pkin(3) * t257 - pkin(9) * t258;
t203 = t264 * pkin(3) - t265 * pkin(9) + t287;
t67 = t316 * t158 - t183 * t369 + t313 * t200 + t203 * t368;
t50 = qJ(5) * t257 + qJD(5) * t264 + t67;
t157 = t224 * t308 - t392 * t225;
t190 = qJD(4) * t235 + t258 * t313;
t191 = -qJD(4) * t234 + t258 * t316;
t75 = pkin(4) * t190 - qJ(5) * t191 - qJD(5) * t235 + t157;
t24 = t307 * t75 + t310 * t50;
t420 = Ifges(3,4) * t314;
t419 = Ifges(4,4) * t256;
t417 = Ifges(5,4) * t313;
t416 = Ifges(5,4) * t316;
t415 = Ifges(6,4) * t307;
t414 = Ifges(6,4) * t310;
t411 = Ifges(6,5) * t310;
t410 = Ifges(3,6) * t311;
t409 = Ifges(6,6) * t307;
t405 = t142 * mrSges(4,2);
t399 = t255 * mrSges(4,3);
t259 = -pkin(8) * t341 + t295;
t397 = t259 * mrSges(3,2);
t275 = t375 * qJD(2);
t260 = qJD(1) * t275;
t396 = t260 * mrSges(3,1);
t391 = Ifges(3,6) * qJD(2);
t388 = t255 * t316;
t119 = t316 * t183 + t313 * t203;
t105 = qJ(5) * t264 + t119;
t192 = t238 * t392 - t308 * t254;
t182 = -t311 * pkin(3) - t192;
t117 = t234 * pkin(4) - t235 * qJ(5) + t182;
t58 = t310 * t105 + t307 * t117;
t111 = -mrSges(6,1) * t347 + mrSges(6,2) * t169;
t179 = mrSges(5,1) * t253 - mrSges(5,3) * t220;
t380 = -t111 + t179;
t145 = t286 * t219;
t379 = t145 + t280;
t146 = t325 * t219;
t378 = t146 + t279;
t364 = Ifges(5,5) * t170 - Ifges(5,6) * t171 + Ifges(5,3) * t250;
t358 = t304 * t368;
t12 = -t39 * mrSges(7,1) + t38 * mrSges(7,2);
t352 = t368 / 0.2e1;
t18 = -t307 * t33 + t310 * t69;
t23 = -t307 * t50 + t310 * t75;
t57 = -t105 * t307 + t310 * t117;
t118 = -t313 * t183 + t203 * t316;
t344 = mrSges(3,3) * t362;
t343 = mrSges(3,3) * t361;
t339 = mrSges(6,1) * t307 + mrSges(6,2) * t310;
t338 = Ifges(5,1) * t316 - t417;
t337 = Ifges(6,1) * t310 - t415;
t336 = -Ifges(5,2) * t313 + t416;
t335 = -Ifges(6,2) * t307 + t414;
t334 = Ifges(5,5) * t316 - Ifges(5,6) * t313;
t333 = -t409 + t411;
t331 = -t18 * t307 + t19 * t310;
t330 = -t307 * t40 + t310 * t41;
t90 = -mrSges(6,2) * t171 + mrSges(6,3) * t133;
t91 = mrSges(6,1) * t171 - mrSges(6,3) * t134;
t329 = -t307 * t91 + t310 * t90;
t198 = t235 * t310 + t264 * t307;
t31 = pkin(5) * t234 - pkin(10) * t198 + t57;
t197 = -t235 * t307 + t264 * t310;
t43 = pkin(10) * t197 + t58;
t13 = t31 * t315 - t312 * t43;
t14 = t31 * t312 + t315 * t43;
t127 = -mrSges(6,2) * t219 + mrSges(6,3) * t347;
t128 = mrSges(6,1) * t219 - mrSges(6,3) * t169;
t326 = t127 * t310 - t128 * t307;
t120 = t197 * t315 - t198 * t312;
t121 = t197 * t312 + t198 * t315;
t68 = -t313 * t158 - t183 * t368 + t200 * t316 - t203 * t369;
t106 = -pkin(4) * t264 - t118;
t324 = t153 * (mrSges(5,1) * t313 + mrSges(5,2) * t316);
t56 = -pkin(4) * t257 - t68;
t306 = -pkin(5) * t310 - pkin(4);
t297 = Ifges(3,4) * t361;
t281 = -pkin(8) * t384 + t302;
t276 = t351 * t313;
t274 = -pkin(8) * t360 + t299;
t273 = t375 * qJD(1);
t272 = -pkin(8) * t362 + t298;
t271 = -mrSges(3,2) * t345 + t343;
t270 = mrSges(3,1) * t345 - t344;
t267 = t325 * t313;
t266 = t286 * t313;
t252 = Ifges(4,4) * t255;
t244 = t251 * mrSges(4,2);
t240 = Ifges(3,1) * t362 + Ifges(3,5) * t345 + t297;
t239 = t391 + (t410 + (t317 * Ifges(3,2) + t420) * t309) * qJD(1);
t228 = -t304 * t385 + t269;
t226 = -mrSges(4,2) * t345 + t399;
t205 = -mrSges(4,1) * t255 - mrSges(4,2) * t256;
t196 = -Ifges(4,1) * t256 + t345 * Ifges(4,5) + t252;
t195 = Ifges(4,2) * t255 + t345 * Ifges(4,6) - t419;
t178 = -mrSges(5,2) * t253 - mrSges(5,3) * t219;
t144 = t191 * t310 + t257 * t307;
t143 = -t191 * t307 + t257 * t310;
t136 = -mrSges(5,2) * t250 - mrSges(5,3) * t171;
t129 = t220 * Ifges(5,5) - t219 * Ifges(5,6) + t253 * Ifges(5,3);
t112 = mrSges(5,1) * t171 + mrSges(5,2) * t170;
t93 = t170 * Ifges(5,1) - t171 * Ifges(5,4) + t250 * Ifges(5,5);
t92 = t170 * Ifges(5,4) - t171 * Ifges(5,2) + t250 * Ifges(5,6);
t87 = t169 * Ifges(6,1) + Ifges(6,4) * t347 + Ifges(6,5) * t219;
t83 = mrSges(7,1) * t218 - mrSges(7,3) * t110;
t82 = -mrSges(7,2) * t218 + mrSges(7,3) * t483;
t81 = -pkin(5) * t390 + t100;
t78 = -pkin(5) * t197 + t106;
t70 = -pkin(5) * t347 + t88;
t61 = Ifges(6,4) * t134 + Ifges(6,2) * t133 + Ifges(6,6) * t171;
t55 = -mrSges(7,1) * t483 + mrSges(7,2) * t110;
t54 = -qJD(6) * t121 + t143 * t315 - t144 * t312;
t53 = qJD(6) * t120 + t143 * t312 + t144 * t315;
t30 = -pkin(5) * t143 + t56;
t28 = -mrSges(7,2) * t171 + mrSges(7,3) * t39;
t27 = mrSges(7,1) * t171 - mrSges(7,3) * t38;
t26 = -pkin(5) * t133 + t42;
t20 = pkin(10) * t143 + t24;
t15 = pkin(5) * t190 - pkin(10) * t144 + t23;
t11 = pkin(10) * t133 + t19;
t7 = pkin(5) * t171 - pkin(10) * t134 + t18;
t4 = -qJD(6) * t14 + t15 * t315 - t20 * t312;
t3 = qJD(6) * t13 + t15 * t312 + t20 * t315;
t2 = -qJD(6) * t6 - t11 * t312 + t315 * t7;
t1 = qJD(6) * t5 + t11 * t315 + t312 * t7;
t9 = [m(5) * (t100 * t67 + t118 * t52 + t119 * t51 + t141 * t182 + t153 * t157 + t68 * t99) + m(6) * (t106 * t42 + t18 * t57 + t19 * t58 + t23 * t40 + t24 * t41 + t56 * t88) + m(7) * (t1 * t14 + t13 * t2 + t26 * t78 + t3 * t6 + t30 * t70 + t4 * t5) + (t141 * t265 - t142 * t264 - t159 * t258 - t160 * t257) * mrSges(4,3) + (-t192 * mrSges(4,3) + Ifges(4,1) * t265 - Ifges(4,4) * t264) * t251 + (t260 * mrSges(3,3) * t314 + (-t273 * mrSges(3,3) - t239 / 0.2e1 - t391 / 0.2e1 + (-t375 * mrSges(3,3) - 0.3e1 / 0.2e1 * t410 + (mrSges(3,1) * t463 - 0.3e1 / 0.2e1 * t420) * t309) * qJD(1) + (t205 + m(4) * (t278 + t374) + qJD(1) * (mrSges(4,1) * t264 + mrSges(4,2) * t265)) * pkin(2)) * t371) * t309 + m(3) * (t259 * t375 - t260 * t281 - t272 * t275 + t273 * t274) + t487 * t234 / 0.2e1 + (t259 * mrSges(3,3) + (-t272 * mrSges(3,3) + t240 / 0.2e1 + qJD(2) * Ifges(3,5) / 0.2e1 + (Ifges(3,5) * t311 - t281 * mrSges(3,3) + (mrSges(3,2) * t463 + 0.3e1 / 0.2e1 * Ifges(3,4) * t317 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t314) * t309) * qJD(1)) * qJD(2)) * t383 + t287 * t244 + (Ifges(6,5) * t198 + Ifges(7,5) * t121 + Ifges(6,6) * t197 + Ifges(7,6) * t120 + t234 * t486) * t441 + (t484 + t485) * t311 / 0.2e1 + (-t193 * mrSges(4,3) + Ifges(5,5) * t433 + Ifges(5,6) * t435 - Ifges(4,4) * t265 + t287 * mrSges(4,1) - Ifges(4,6) * t311 / 0.2e1 + (Ifges(5,3) / 0.2e1 + Ifges(4,2)) * t264) * t250 - t311 * t405 - t311 * t397 + t377 * t157 + t264 * t364 / 0.2e1 + t58 * t90 + t57 * t91 + t3 * t82 + t4 * t83 + t78 * t12 + t70 * (-mrSges(7,1) * t54 + mrSges(7,2) * t53) + t30 * t55 + t345 * (Ifges(4,5) * t258 - Ifges(4,6) * t257) / 0.2e1 + t54 * t458 + (Ifges(7,4) * t121 + Ifges(7,2) * t120 + Ifges(7,6) * t234) * t459 + (Ifges(7,1) * t121 + Ifges(7,4) * t120 + Ifges(7,5) * t234) * t460 + t121 * t461 + t120 * t462 + (Ifges(6,4) * t144 + Ifges(6,2) * t143 + Ifges(6,6) * t190) * t445 + (Ifges(6,1) * t198 + Ifges(6,4) * t197 + Ifges(6,5) * t234) * t446 + (Ifges(6,4) * t198 + Ifges(6,2) * t197 + Ifges(6,6) * t234) * t447 + t191 * t448 + t190 * t450 + (Ifges(7,1) * t53 + Ifges(7,4) * t54 + Ifges(7,5) * t190) * t451 + (Ifges(7,4) * t53 + Ifges(7,2) * t54 + Ifges(7,6) * t190) * t453 + t198 * t456 + t53 * t457 + t93 * t433 + t92 * t435 + (Ifges(6,5) * t144 + Ifges(6,6) * t143 + Ifges(6,3) * t190) * t436 + (Ifges(5,4) * t191 - Ifges(5,2) * t190 + Ifges(5,6) * t257) * t437 + (Ifges(7,5) * t53 + Ifges(7,6) * t54 + Ifges(7,3) * t190) * t438 + (Ifges(5,4) * t235 - Ifges(5,2) * t234 + Ifges(5,6) * t264) * t442 + (Ifges(5,1) * t235 - Ifges(5,4) * t234 + Ifges(5,5) * t264) * t443 + (Ifges(6,1) * t144 + Ifges(6,4) * t143 + Ifges(6,5) * t190) * t444 + (Ifges(4,1) * t258 - Ifges(4,4) * t257) * t431 + t13 * t27 + t14 * t28 + t106 * t77 + m(4) * (-t141 * t192 + t142 * t193 - t157 * t159 + t158 * t160) + t56 * t111 + t26 * (-mrSges(7,1) * t120 + mrSges(7,2) * t121) + t24 * t127 + t23 * t128 + t118 * t135 + t119 * t136 + t143 * t86 / 0.2e1 + t144 * t87 / 0.2e1 + t88 * (-mrSges(6,1) * t143 + mrSges(6,2) * t144) + t67 * t178 + t68 * t179 + t182 * t112 + t5 * (mrSges(7,1) * t190 - mrSges(7,3) * t53) + t6 * (-mrSges(7,2) * t190 + mrSges(7,3) * t54) + t41 * (-mrSges(6,2) * t190 + mrSges(6,3) * t143) + t40 * (mrSges(6,1) * t190 - mrSges(6,3) * t144) + t153 * (mrSges(5,1) * t190 + mrSges(5,2) * t191) + t197 * t61 / 0.2e1 + t42 * (-mrSges(6,1) * t197 + mrSges(6,2) * t198) + t158 * t226 + t1 * (-mrSges(7,2) * t234 + mrSges(7,3) * t120) + t2 * (mrSges(7,1) * t234 - mrSges(7,3) * t121) + t19 * (-mrSges(6,2) * t234 + mrSges(6,3) * t197) + t18 * (mrSges(6,1) * t234 - mrSges(6,3) * t198) + t190 * t488 - t311 * t396 + (-mrSges(4,1) * t311 + mrSges(5,1) * t234 + mrSges(5,2) * t235) * t141 + t257 * t129 / 0.2e1 + t100 * (-mrSges(5,2) * t257 - mrSges(5,3) * t190) + t99 * (mrSges(5,1) * t257 - mrSges(5,3) * t191) + t220 * (Ifges(5,1) * t191 - Ifges(5,4) * t190 + Ifges(5,5) * t257) / 0.2e1 + t253 * (Ifges(5,5) * t191 - Ifges(5,6) * t190 + Ifges(5,3) * t257) / 0.2e1 - t257 * t195 / 0.2e1 + t258 * t196 / 0.2e1 + t255 * (Ifges(4,4) * t258 - Ifges(4,2) * t257) / 0.2e1 + t51 * (-mrSges(5,2) * t264 - mrSges(5,3) * t234) + t52 * (mrSges(5,1) * t264 - mrSges(5,3) * t235) + t274 * t271 - t275 * t270 + t278 * (mrSges(4,1) * t257 + mrSges(4,2) * t258); (-t314 * (Ifges(3,1) * t317 - t420) / 0.2e1 + pkin(1) * (mrSges(3,1) * t314 + mrSges(3,2) * t317)) * qJD(1) ^ 2 * t309 ^ 2 + t467 * t368 + (t100 * t389 + t388 * t99 + t468) * mrSges(5,3) + t469 * t339 - t377 * t184 + (t1 * t316 - t26 * t267 + t389 * t6 - t470 * t70) * mrSges(7,2) + (-t2 * t316 + t26 * t266 - t389 * t5 + t471 * t70) * mrSges(7,1) + (-t1 * t266 + t2 * t267 + t470 * t5 - t471 * t6) * mrSges(7,3) + t474 * t55 + t475 * t127 + t476 * t128 + (t40 * (mrSges(6,1) * t313 - mrSges(6,3) * t381) + t41 * (-mrSges(6,2) * t313 - mrSges(6,3) * t385) + t324) * qJD(4) + (-t359 - t116) * t178 + (t159 * t184 - t160 * t185 - t278 * t346 + (-t141 * t392 + t142 * t308) * pkin(2)) * m(4) + (t344 + t270) * t273 + (-t250 * t425 - t251 * t363) * mrSges(4,3) - t396 + (-t100 * t116 - t115 * t99 + t141 * t305 - t153 * t184) * m(5) + (-t102 * t88 + t18 * t228 + t19 * t229 + t40 * t476 + t41 * t475) * m(6) + (t479 * t313 + ((-t100 * t313 - t316 * t99) * qJD(4) + t468) * m(5) + t469 * m(6) + t316 * t136) * t304 + (t343 - t271) * t272 + (Ifges(7,5) * t211 + Ifges(7,6) * t212) * t438 - t397 + t484 + (t310 * t352 - t202 / 0.2e1) * t87 + t99 * mrSges(5,1) * t256 + (-t358 - t115) * t179 - ((-Ifges(3,2) * t362 + t240 + t297) * t317 + t345 * (-Ifges(3,6) * t314 + t412)) * t373 / 0.2e1 - t487 * t316 / 0.2e1 + (Ifges(7,5) * t451 + Ifges(7,6) * t453 + Ifges(7,3) * t438 + t450 - t464 + t488) * t369 + (-mrSges(5,1) * t316 + mrSges(5,2) * t313 - mrSges(4,1)) * t141 + (Ifges(7,4) * t211 + Ifges(7,2) * t212) * t453 + (-t388 / 0.2e1 + t352) * t131 + (Ifges(4,1) * t255 + t129 + t419) * t256 / 0.2e1 - (Ifges(4,2) * t256 + t196 + t252) * t255 / 0.2e1 + (t169 * (Ifges(6,5) * t313 + t316 * t337) + t347 * (Ifges(6,6) * t313 + t316 * t335) + t253 * t334 + t220 * t338) * qJD(4) / 0.2e1 - t100 * mrSges(5,2) * t256 + (Ifges(7,1) * t211 + Ifges(7,4) * t212) * t451 + (-Ifges(7,5) * t267 - Ifges(7,6) * t266 + t313 * t333 - t316 * t486) * t441 - t160 * t398 - t41 * (-mrSges(6,2) * t389 + mrSges(6,3) * t201) - t169 * (Ifges(6,1) * t202 + Ifges(6,4) * t201 + Ifges(6,5) * t389) / 0.2e1 - t40 * (mrSges(6,1) * t389 - mrSges(6,3) * t202) - t61 * t386 / 0.2e1 + t19 * (mrSges(6,2) * t316 - mrSges(6,3) * t386) + t18 * (-mrSges(6,1) * t316 - mrSges(6,3) * t382) - t336 * t370 / 0.2e1 + (Ifges(6,3) * t313 + t316 * t333) * t370 / 0.2e1 + t239 * t362 / 0.2e1 + t250 * (Ifges(5,5) * t313 + Ifges(5,6) * t316) / 0.2e1 + t316 * t92 / 0.2e1 + t313 * t93 / 0.2e1 + (t358 - t102) * t111 - t205 * t346 - t345 * (Ifges(4,5) * t255 + Ifges(4,6) * t256) / 0.2e1 + t212 * t458 - t267 * t461 - t266 * t462 + (-Ifges(6,5) * t316 + t313 * t337) * t446 + (-Ifges(6,6) * t316 + t313 * t335) * t447 + t389 * t449 + (Ifges(7,1) * t125 + Ifges(7,4) * t124 + Ifges(7,5) * t389) * t452 + (Ifges(7,4) * t125 + Ifges(7,2) * t124 + Ifges(7,6) * t389) * t454 + t201 * t455 + t382 * t456 + t211 * t457 + t195 * t431 + (-Ifges(5,6) * t256 + t255 * t336) * t436 + (Ifges(6,5) * t202 + Ifges(6,6) * t201 + Ifges(6,3) * t389) * t437 + (Ifges(7,5) * t125 + Ifges(7,6) * t124 + Ifges(7,3) * t389) * t439 + (Ifges(5,2) * t316 + t417) * t442 + (Ifges(5,1) * t313 + t416) * t443 - t405 - t347 * (Ifges(6,4) * t202 + Ifges(6,2) * t201 + Ifges(6,6) * t389) / 0.2e1 + (-Ifges(7,4) * t267 - Ifges(7,2) * t266 - Ifges(7,6) * t316) * t459 + (-Ifges(7,1) * t267 - Ifges(7,4) * t266 - Ifges(7,5) * t316) * t460 - t255 * t324 - t480 * t389 / 0.2e1 + t481 * t83 + t482 * t82 + (t1 * t140 + t139 * t2 + t26 * t276 + t474 * t70 + t481 * t5 + t482 * t6) * m(7) - t124 * t45 / 0.2e1 - t125 * t46 / 0.2e1 + t139 * t27 + t140 * t28 - t88 * (-mrSges(6,1) * t201 + mrSges(6,2) * t202) - t253 * (-Ifges(5,3) * t256 + t255 * t334) / 0.2e1 - t220 * (-Ifges(5,5) * t256 + t255 * t338) / 0.2e1 - t185 * t226 + t228 * t91 + t229 * t90 + t159 * t399 + t276 * t12 - t278 * (-mrSges(4,1) * t256 + mrSges(4,2) * t255) + t305 * t112 - Ifges(3,6) * t341; t250 * mrSges(4,1) - t202 * t127 - t201 * t128 - t255 * t226 - t266 * t27 - t267 * t28 + t244 - t471 * t83 - t470 * t82 + t377 * t256 + (-t255 * t178 - t12 + (t178 + t326) * qJD(4) - t479) * t316 + (t136 + t329 + t466 * (t55 - t380)) * t313 + (-t1 * t267 - t2 * t266 - t26 * t316 + (t369 - t389) * t70 - t470 * t6 - t471 * t5) * m(7) + (-t201 * t40 - t202 * t41 - t88 * t389 - t316 * t42 + t331 * t313 + (t313 * t88 + t316 * t330) * qJD(4)) * m(6) + (t153 * t256 + t313 * t51 + t316 * t52 + t466 * (t100 * t316 - t313 * t99)) * m(5) + (-t159 * t256 - t160 * t255 + t328) * m(4); (-t404 / 0.2e1 - t403 / 0.2e1 + t400 / 0.2e1 + t41 * mrSges(6,2) - t40 * mrSges(6,1) - t85 / 0.2e1 - t44 / 0.2e1 - t407 / 0.2e1 - t406 / 0.2e1 - t402 / 0.2e1 + t449 - t153 * mrSges(5,1) + t418 / 0.2e1 + t464) * t220 + (t87 * t429 + t448 + t153 * mrSges(5,2) - t216 / 0.2e1 + t401 / 0.2e1 + t88 * t339 + t335 * t445 + t337 * t444 + (t411 / 0.2e1 - t409 / 0.2e1) * t219 + (-t307 * t41 - t310 * t40) * mrSges(6,3) + (-Ifges(5,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(6,3) / 0.2e1) * t220 + t467) * t219 + t477 * t82 + (t1 * t237 + t2 * t236 + t26 * t306 + t477 * t6 + t478 * t5 - t70 * t81) * m(7) + t478 * t83 + t364 + (Ifges(6,5) * t307 + Ifges(7,5) * t286 + Ifges(6,6) * t310 - Ifges(7,6) * t325) * t441 + (-t1 * t325 - t2 * t286 + t378 * t5 - t379 * t6) * mrSges(7,3) + (Ifges(7,4) * t286 - Ifges(7,2) * t325) * t459 + (Ifges(7,1) * t286 - Ifges(7,4) * t325) * t460 + t26 * (mrSges(7,1) * t325 + mrSges(7,2) * t286) - t325 * t462 + (-t145 / 0.2e1 - t280 / 0.2e1) * t45 + (-Ifges(7,5) * t279 - Ifges(7,6) * t280) * t438 + (-Ifges(7,1) * t279 - Ifges(7,4) * t280) * t451 + (-Ifges(7,4) * t279 - Ifges(7,2) * t280) * t453 + t380 * t100 + (mrSges(7,1) * t379 - mrSges(7,2) * t378) * t70 - t81 * t55 - pkin(4) * t77 - t51 * mrSges(5,2) + t52 * mrSges(5,1) + t286 * t461 + (Ifges(6,1) * t307 + t414) * t446 + (Ifges(6,2) * t310 + t415) * t447 + (Ifges(7,1) * t146 + Ifges(7,4) * t145) * t452 + (Ifges(7,4) * t146 + Ifges(7,2) * t145) * t454 + t307 * t456 + (Ifges(7,5) * t146 + Ifges(7,6) * t145) * t439 + t61 * t429 + (-t146 / 0.2e1 - t279 / 0.2e1) * t46 - t72 * t127 - t71 * t128 - t99 * t178 + t326 * qJD(5) + (-pkin(4) * t42 + qJ(5) * t331 + qJD(5) * t330 - t100 * t88 - t40 * t71 - t41 * t72) * m(6) + t329 * qJ(5) + t331 * mrSges(6,3) + t236 * t27 + t237 * t28 + t306 * t12 + t42 * (-mrSges(6,1) * t310 + mrSges(6,2) * t307); t110 * t83 - t483 * t82 - t347 * t127 + t169 * t128 + t12 + t77 + (t110 * t5 - t483 * t6 + t26) * m(7) + (t169 * t40 - t347 * t41 + t42) * m(6); -t1 * mrSges(7,2) + t2 * mrSges(7,1) - t70 * (mrSges(7,1) * t110 + mrSges(7,2) * t483) + t45 * t451 + (Ifges(7,5) * t483 - Ifges(7,6) * t110) * t439 + (Ifges(7,1) * t483 - t413) * t452 - t5 * t82 + t6 * t83 + (t110 * t6 + t483 * t5) * mrSges(7,3) + t8 + (-Ifges(7,2) * t110 + t104 + t46) * t454;];
tauc  = t9(:);
