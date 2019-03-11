% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:16:42
% EndTime: 2019-03-09 13:17:08
% DurationCPUTime: 13.01s
% Computational Cost: add. (22359->645), mult. (57354->881), div. (0->0), fcn. (44567->10), ass. (0->309)
t303 = sin(qJ(6));
t304 = sin(qJ(5));
t307 = cos(qJ(6));
t308 = cos(qJ(5));
t275 = t303 * t308 + t304 * t307;
t440 = qJD(5) + qJD(6);
t215 = t440 * t275;
t301 = sin(pkin(11));
t302 = cos(pkin(11));
t306 = sin(qJ(2));
t310 = cos(qJ(2));
t269 = -t301 * t306 + t302 * t310;
t258 = t269 * qJD(1);
t357 = qJD(1) * t310;
t358 = qJD(1) * t306;
t259 = -t301 * t357 - t302 * t358;
t305 = sin(qJ(4));
t309 = cos(qJ(4));
t442 = t309 * t258 + t305 * t259;
t467 = t275 * t442;
t446 = t467 - t215;
t326 = t303 * t304 - t307 * t308;
t214 = t440 * t326;
t466 = t326 * t442;
t445 = t466 - t214;
t293 = pkin(2) * t302 + pkin(3);
t404 = pkin(2) * t301;
t253 = t293 * t309 - t305 * t404;
t240 = t253 * qJD(4);
t398 = -qJ(3) - pkin(7);
t283 = t398 * t306;
t276 = qJD(1) * t283;
t284 = t398 * t310;
t277 = qJD(1) * t284;
t362 = t302 * t277;
t212 = -t276 * t301 + t362;
t402 = pkin(8) * t258;
t189 = t212 - t402;
t262 = t301 * t277;
t213 = t302 * t276 + t262;
t401 = pkin(8) * t259;
t190 = t213 + t401;
t131 = t189 * t305 + t190 * t309;
t328 = t258 * t305 - t309 * t259;
t160 = pkin(4) * t328 - pkin(9) * t442;
t297 = pkin(2) * t358;
t227 = -pkin(3) * t259 + t297;
t132 = t160 + t227;
t76 = t308 * t131 + t304 * t132;
t480 = t240 * t308 - t76;
t75 = -t131 * t304 + t308 * t132;
t479 = -t240 * t304 - t75;
t270 = t301 * t310 + t302 * t306;
t260 = t270 * qJD(2);
t246 = qJD(1) * t260;
t261 = t269 * qJD(2);
t247 = qJD(1) * t261;
t155 = qJD(4) * t442 - t246 * t305 + t247 * t309;
t156 = qJD(4) * t328 + t309 * t246 + t247 * t305;
t300 = qJD(2) + qJD(4);
t186 = t300 * t308 - t304 * t328;
t268 = qJD(2) * pkin(2) + t276;
t208 = t302 * t268 + t262;
t177 = qJD(2) * pkin(3) + t208 + t401;
t209 = t301 * t268 - t362;
t181 = t209 + t402;
t117 = t177 * t305 + t181 * t309;
t113 = pkin(9) * t300 + t117;
t295 = -pkin(2) * t310 - pkin(1);
t359 = qJD(1) * t295;
t280 = qJD(3) + t359;
t216 = -t258 * pkin(3) + t280;
t121 = -pkin(4) * t442 - pkin(9) * t328 + t216;
t69 = t113 * t308 + t121 * t304;
t55 = pkin(10) * t186 + t69;
t377 = t303 * t55;
t200 = qJD(5) - t442;
t188 = t300 * t304 + t308 * t328;
t68 = -t113 * t304 + t308 * t121;
t54 = -pkin(10) * t188 + t68;
t47 = pkin(5) * t200 + t54;
t20 = t307 * t47 - t377;
t376 = t307 * t55;
t21 = t303 * t47 + t376;
t109 = -qJD(5) * t188 - t155 * t304;
t354 = qJD(5) * t308;
t355 = qJD(5) * t304;
t116 = t177 * t309 - t305 * t181;
t346 = qJD(2) * t398;
t255 = qJD(3) * t310 + t306 * t346;
t233 = t255 * qJD(1);
t256 = -t306 * qJD(3) + t310 * t346;
t234 = t256 * qJD(1);
t187 = t302 * t233 + t301 * t234;
t166 = -pkin(8) * t246 + t187;
t185 = -t233 * t301 + t234 * t302;
t320 = -pkin(8) * t247 + t185;
t60 = qJD(4) * t116 + t309 * t166 + t305 * t320;
t292 = qJD(2) * t297;
t217 = pkin(3) * t246 + t292;
t84 = pkin(4) * t156 - pkin(9) * t155 + t217;
t16 = -t113 * t355 + t121 * t354 + t304 * t84 + t308 * t60;
t12 = pkin(10) * t109 + t16;
t108 = qJD(5) * t186 + t155 * t308;
t17 = -qJD(5) * t69 - t304 * t60 + t308 * t84;
t7 = pkin(5) * t156 - pkin(10) * t108 + t17;
t3 = qJD(6) * t20 + t12 * t307 + t303 * t7;
t334 = Ifges(6,5) * t304 + Ifges(6,6) * t308;
t394 = Ifges(6,4) * t304;
t336 = Ifges(6,2) * t308 + t394;
t393 = Ifges(6,4) * t308;
t338 = Ifges(6,1) * t304 + t393;
t341 = mrSges(6,1) * t308 - mrSges(6,2) * t304;
t61 = qJD(4) * t117 + t166 * t305 - t309 * t320;
t37 = -pkin(5) * t109 + t61;
t4 = -qJD(6) * t21 - t12 * t303 + t307 * t7;
t405 = t308 / 0.2e1;
t194 = qJD(6) + t200;
t413 = t194 / 0.2e1;
t414 = -t194 / 0.2e1;
t418 = t156 / 0.2e1;
t129 = t186 * t303 + t188 * t307;
t419 = t129 / 0.2e1;
t420 = -t129 / 0.2e1;
t342 = t307 * t186 - t188 * t303;
t421 = t342 / 0.2e1;
t422 = -t342 / 0.2e1;
t423 = t109 / 0.2e1;
t424 = t108 / 0.2e1;
t122 = Ifges(7,4) * t342;
t65 = Ifges(7,1) * t129 + Ifges(7,5) * t194 + t122;
t429 = t65 / 0.2e1;
t430 = -t65 / 0.2e1;
t392 = Ifges(7,4) * t129;
t64 = Ifges(7,2) * t342 + Ifges(7,6) * t194 + t392;
t431 = t64 / 0.2e1;
t432 = -t64 / 0.2e1;
t41 = -qJD(6) * t129 - t108 * t303 + t109 * t307;
t433 = t41 / 0.2e1;
t40 = qJD(6) * t342 + t108 * t307 + t109 * t303;
t434 = t40 / 0.2e1;
t435 = Ifges(7,1) * t434 + Ifges(7,4) * t433 + Ifges(7,5) * t418;
t436 = Ifges(7,4) * t434 + Ifges(7,2) * t433 + Ifges(7,6) * t418;
t44 = Ifges(6,4) * t108 + Ifges(6,2) * t109 + Ifges(6,6) * t156;
t45 = t108 * Ifges(6,1) + t109 * Ifges(6,4) + t156 * Ifges(6,5);
t379 = t188 * Ifges(6,4);
t104 = t186 * Ifges(6,2) + t200 * Ifges(6,6) + t379;
t184 = Ifges(6,4) * t186;
t105 = Ifges(6,1) * t188 + Ifges(6,5) * t200 + t184;
t112 = -pkin(4) * t300 - t116;
t335 = Ifges(6,5) * t308 - Ifges(6,6) * t304;
t337 = -Ifges(6,2) * t304 + t393;
t339 = Ifges(6,1) * t308 - t394;
t340 = mrSges(6,1) * t304 + mrSges(6,2) * t308;
t406 = -t304 / 0.2e1;
t415 = t188 / 0.2e1;
t459 = t200 * t335 / 0.2e1 + t339 * t415 + t186 * t337 / 0.2e1 + t105 * t405 + t104 * t406 + t112 * t340;
t463 = t16 * t308 - t17 * t304;
t94 = -pkin(5) * t186 + t112;
t478 = (-t446 * mrSges(7,1) + mrSges(7,2) * t445) * t94 - t467 * t432 + (-Ifges(7,1) * t466 - Ifges(7,4) * t467) * t420 + (-Ifges(7,4) * t466 - Ifges(7,2) * t467) * t422 + (-Ifges(7,5) * t466 - Ifges(7,6) * t467) * t414 - t466 * t430 + (-t20 * t445 + t21 * t446 - t4 * t275 - t3 * t326) * mrSges(7,3) + t463 * mrSges(6,3) + t459 * qJD(5) + (-mrSges(5,1) - t341) * t61 + t275 * t435 + t336 * t423 + t338 * t424 - t214 * t429 - t215 * t431 - t60 * mrSges(5,2) + t44 * t405 + (Ifges(7,4) * t275 - Ifges(7,2) * t326) * t433 + (Ifges(7,1) * t275 - Ifges(7,4) * t326) * t434 + t37 * (mrSges(7,1) * t326 + mrSges(7,2) * t275) - t326 * t436 + (Ifges(7,5) * t275 - Ifges(7,6) * t326 + t334) * t418 + (-Ifges(7,4) * t214 - Ifges(7,2) * t215) * t421 + (-Ifges(7,5) * t214 - Ifges(7,6) * t215) * t413 + (-Ifges(7,1) * t214 - Ifges(7,4) * t215) * t419 + Ifges(5,5) * t155 - Ifges(5,6) * t156 + t304 * t45 / 0.2e1;
t465 = t442 * t304;
t477 = pkin(5) * t465;
t476 = pkin(10) * t465;
t254 = t305 * t293 + t309 * t404;
t251 = pkin(9) + t254;
t399 = -pkin(10) - t251;
t345 = qJD(5) * t399;
t474 = t304 * t345 + t476 + t480;
t299 = t308 * pkin(10);
t460 = pkin(5) * t328 - t299 * t442;
t473 = t308 * t345 - t460 + t479;
t428 = -pkin(10) - pkin(9);
t352 = qJD(5) * t428;
t83 = t308 * t116 + t304 * t160;
t472 = t304 * t352 + t476 - t83;
t82 = -t116 * t304 + t308 * t160;
t471 = t308 * t352 - t460 - t82;
t151 = Ifges(7,3) * t156;
t469 = t151 + (Ifges(7,5) * t342 - Ifges(7,6) * t129) * t414 + (t129 * t21 + t20 * t342) * mrSges(7,3) - t94 * (mrSges(7,1) * t129 + mrSges(7,2) * t342) + (Ifges(7,1) * t342 - t392) * t420;
t443 = t254 * qJD(4) + t309 * t189 - t190 * t305;
t332 = t304 * t69 + t308 * t68;
t324 = t332 * mrSges(6,3);
t199 = Ifges(5,4) * t442;
t397 = Ifges(5,1) * t328;
t438 = t199 / 0.2e1 + t397 / 0.2e1;
t464 = -t216 * mrSges(5,2) - Ifges(5,5) * t300 + t324 - t438 - t459;
t462 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + Ifges(7,5) * t40 + Ifges(7,6) * t41;
t461 = -Ifges(7,2) * t129 + t122;
t220 = t399 * t304;
t221 = t251 * t308 + t299;
t171 = t220 * t303 + t221 * t307;
t455 = -qJD(6) * t171 - t303 * t474 + t307 * t473;
t170 = t220 * t307 - t221 * t303;
t454 = qJD(6) * t170 + t303 * t473 + t307 * t474;
t286 = t428 * t304;
t287 = pkin(9) * t308 + t299;
t224 = t286 * t303 + t287 * t307;
t449 = -qJD(6) * t224 - t303 * t472 + t307 * t471;
t223 = t286 * t307 - t287 * t303;
t448 = qJD(6) * t223 + t303 * t471 + t307 * t472;
t211 = t269 * t305 + t270 * t309;
t162 = t326 * t211;
t353 = pkin(5) * t355;
t447 = t353 + t443 - t477;
t218 = t302 * t283 + t284 * t301;
t197 = -pkin(8) * t270 + t218;
t219 = t301 * t283 - t302 * t284;
t198 = pkin(8) * t269 + t219;
t147 = t197 * t305 + t198 * t309;
t140 = t308 * t147;
t235 = -t269 * pkin(3) + t295;
t327 = t309 * t269 - t270 * t305;
t141 = -pkin(4) * t327 - t211 * pkin(9) + t235;
t86 = t304 * t141 + t140;
t444 = t309 * t197 - t198 * t305;
t441 = -t199 / 0.2e1 + t464;
t437 = t17 * mrSges(6,1) - t16 * mrSges(6,2) + Ifges(6,5) * t108 + Ifges(6,6) * t109 + t462;
t427 = pkin(1) * mrSges(3,1);
t426 = pkin(1) * mrSges(3,2);
t417 = -t186 / 0.2e1;
t416 = -t188 / 0.2e1;
t411 = -t200 / 0.2e1;
t409 = -t259 / 0.2e1;
t408 = -t260 / 0.2e1;
t407 = t261 / 0.2e1;
t403 = pkin(5) * t308;
t396 = Ifges(3,4) * t306;
t388 = Ifges(5,2) * t442;
t382 = t444 * t61;
t378 = t259 * Ifges(4,4);
t375 = Ifges(3,5) * qJD(2);
t374 = Ifges(3,6) * qJD(2);
t373 = qJD(2) * mrSges(3,1);
t372 = qJD(2) * mrSges(3,2);
t366 = t211 * t304;
t361 = -mrSges(5,1) * t300 - mrSges(6,1) * t186 + mrSges(6,2) * t188 + mrSges(5,3) * t328;
t196 = t302 * t255 + t301 * t256;
t356 = qJD(2) * t306;
t351 = t375 / 0.2e1;
t350 = -t374 / 0.2e1;
t228 = pkin(2) * t356 + pkin(3) * t260;
t195 = -t255 * t301 + t302 * t256;
t172 = -pkin(8) * t261 + t195;
t173 = -pkin(8) * t260 + t196;
t78 = qJD(4) * t444 + t172 * t305 + t173 * t309;
t163 = qJD(4) * t327 - t260 * t305 + t261 * t309;
t164 = qJD(4) * t211 + t309 * t260 + t261 * t305;
t89 = pkin(4) * t164 - pkin(9) * t163 + t228;
t347 = -t304 * t78 + t308 * t89;
t344 = t246 * mrSges(4,1) + t247 * mrSges(4,2);
t343 = t156 * mrSges(5,1) + t155 * mrSges(5,2);
t85 = t308 * t141 - t147 * t304;
t250 = -pkin(4) - t253;
t333 = -t16 * t304 - t17 * t308;
t66 = -pkin(5) * t327 - t211 * t299 + t85;
t73 = -pkin(10) * t366 + t86;
t30 = -t303 * t73 + t307 * t66;
t31 = t303 * t66 + t307 * t73;
t331 = t304 * t68 - t308 * t69;
t71 = mrSges(6,1) * t156 - mrSges(6,3) * t108;
t72 = -mrSges(6,2) * t156 + mrSges(6,3) * t109;
t330 = -t304 * t71 + t308 * t72;
t136 = -mrSges(6,2) * t200 + mrSges(6,3) * t186;
t137 = mrSges(6,1) * t200 - mrSges(6,3) * t188;
t329 = t136 * t308 - t137 * t304;
t191 = -mrSges(5,2) * t300 + mrSges(5,3) * t442;
t325 = t191 + t329;
t323 = t163 * t304 + t211 * t354;
t22 = t141 * t354 - t147 * t355 + t304 * t89 + t308 * t78;
t79 = qJD(4) * t147 - t309 * t172 + t173 * t305;
t316 = m(6) * (-qJD(5) * t332 + t463);
t314 = t21 * mrSges(7,2) + t69 * mrSges(6,2) - Ifges(6,3) * t200 - Ifges(6,6) * t186 - Ifges(6,5) * t188 + Ifges(5,6) * t300 + t388 / 0.2e1 + Ifges(5,4) * t328 - Ifges(7,3) * t194 - Ifges(7,6) * t342 - Ifges(7,5) * t129 - t20 * mrSges(7,1) - t216 * mrSges(5,1) - t68 * mrSges(6,1);
t312 = -t314 - t388 / 0.2e1;
t296 = Ifges(3,4) * t357;
t294 = -pkin(4) - t403;
t282 = mrSges(3,3) * t357 - t372;
t281 = -mrSges(3,3) * t358 + t373;
t267 = Ifges(3,1) * t358 + t296 + t375;
t266 = t374 + (t310 * Ifges(3,2) + t396) * qJD(1);
t252 = Ifges(4,4) * t258;
t232 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t259;
t231 = -qJD(2) * mrSges(4,2) + t258 * mrSges(4,3);
t229 = t250 - t403;
t207 = -mrSges(4,1) * t258 - mrSges(4,2) * t259;
t202 = -t259 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t252;
t201 = t258 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t378;
t161 = t275 * t211;
t159 = -mrSges(5,1) * t442 + mrSges(5,2) * t328;
t152 = Ifges(6,3) * t156;
t111 = pkin(5) * t366 - t444;
t99 = mrSges(7,1) * t194 - mrSges(7,3) * t129;
t98 = -mrSges(7,2) * t194 + mrSges(7,3) * t342;
t95 = t117 + t477;
t77 = -mrSges(7,1) * t342 + mrSges(7,2) * t129;
t56 = -mrSges(6,1) * t109 + mrSges(6,2) * t108;
t49 = t162 * t440 - t275 * t163;
t48 = -t163 * t326 - t211 * t215;
t46 = pkin(5) * t323 + t79;
t33 = -mrSges(7,2) * t156 + mrSges(7,3) * t41;
t32 = mrSges(7,1) * t156 - mrSges(7,3) * t40;
t25 = t307 * t54 - t377;
t24 = -t303 * t54 - t376;
t23 = -qJD(5) * t86 + t347;
t18 = -pkin(10) * t323 + t22;
t14 = -t163 * t299 + pkin(5) * t164 + (-t140 + (pkin(10) * t211 - t141) * t304) * qJD(5) + t347;
t13 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t6 = -qJD(6) * t31 + t14 * t307 - t18 * t303;
t5 = qJD(6) * t30 + t14 * t303 + t18 * t307;
t1 = [(-t116 * t163 - t117 * t164 - t147 * t156 - t155 * t444) * mrSges(5,3) - t444 * t56 + (t438 - t464) * t163 + m(4) * (t185 * t218 + t187 * t219 + t195 * t208 + t196 * t209) + (-Ifges(7,5) * t162 - Ifges(7,6) * t161) * t418 + (-t161 * t3 + t162 * t4 - t20 * t48 + t21 * t49) * mrSges(7,3) + t37 * (mrSges(7,1) * t161 - mrSges(7,2) * t162) + t280 * (mrSges(4,1) * t260 + mrSges(4,2) * t261) + m(7) * (t111 * t37 + t20 * t6 + t21 * t5 + t3 * t31 + t30 * t4 + t46 * t94) + (-Ifges(7,4) * t162 - Ifges(7,2) * t161) * t433 + (-Ifges(7,1) * t162 - Ifges(7,4) * t161) * t434 + t30 * t32 + t31 * t33 + t312 * t164 - t162 * t435 - t161 * t436 + (Ifges(7,4) * t48 + Ifges(7,2) * t49) * t421 + t48 * t429 + t49 * t431 + t201 * t408 + (Ifges(7,5) * t48 + Ifges(7,6) * t49) * t413 + (Ifges(7,1) * t48 + Ifges(7,4) * t49) * t419 + t202 * t407 + t46 * t77 - (-t60 * mrSges(5,3) + t151 / 0.2e1 + t152 / 0.2e1 - Ifges(5,4) * t155 + t217 * mrSges(5,1) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + Ifges(5,2)) * t156 + t437) * t327 + t85 * t71 + t86 * t72 + (-t185 * t270 + t187 * t269 - t208 * t261 - t209 * t260 - t218 * t247 - t219 * t246) * mrSges(4,3) + (-t270 * t246 + t269 * t247 + t258 * t407 - t260 * t409) * Ifges(4,4) + (-t269 * t246 + t258 * t408) * Ifges(4,2) + m(6) * (t112 * t79 + t16 * t86 + t17 * t85 + t22 * t69 + t23 * t68 - t382) + m(5) * (-t116 * t79 + t117 * t78 + t147 * t60 + t216 * t228 + t217 * t235 - t382) + t94 * (-mrSges(7,1) * t49 + mrSges(7,2) * t48) + t5 * t98 + t6 * t99 + t111 * t13 + t235 * t343 + t295 * t344 + t22 * t136 + t23 * t137 + t361 * t79 + t78 * t191 + t228 * t159 + t196 * t231 + t195 * t232 + (t270 * t247 + t261 * t409) * Ifges(4,1) + (t339 * t424 + t337 * t423 + t335 * t418 + t44 * t406 + t45 * t405 + Ifges(5,1) * t155 - Ifges(5,4) * t156 + t217 * mrSges(5,2) + (mrSges(5,3) + t340) * t61 + t333 * mrSges(6,3) + (t112 * t341 + t336 * t417 + t338 * t416 + t334 * t411 + t105 * t406 - t308 * t104 / 0.2e1 + t331 * mrSges(6,3)) * qJD(5)) * t211 + (Ifges(4,5) * t407 + Ifges(4,6) * t408 + (t267 / 0.2e1 - pkin(7) * t281 + t351 + (-0.2e1 * t426 + 0.3e1 / 0.2e1 * Ifges(3,4) * t310) * qJD(1)) * t310) * qJD(2) + (-t266 / 0.2e1 - pkin(7) * t282 + t350 + (-0.2e1 * t427 - 0.3e1 / 0.2e1 * t396 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t310) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t269 + mrSges(4,2) * t270) + m(4) * (t280 + t359) + t207) * pkin(2)) * t356; (t443 * t112 + t250 * t61 + t479 * t68 + t480 * t69) * m(6) - t324 * qJD(5) + ((-t136 * t304 - t137 * t308) * qJD(5) + t330 + t316) * t251 + t478 - m(4) * (t208 * t212 + t209 * t213) - (Ifges(4,2) * t259 + t202 + t252) * t258 / 0.2e1 - t312 * t328 + t443 * t361 + (-t216 * t227 - t253 * t61 + t254 * t60 + (-t131 + t240) * t117 - t443 * t116) * m(5) + t447 * t77 + t201 * t409 + m(4) * (t185 * t302 + t187 * t301) * pkin(2) + t454 * t98 + t455 * t99 + t325 * t240 + (t170 * t4 + t171 * t3 + t20 * t455 + t21 * t454 + t229 * t37 + t447 * t94) * m(7) + (t116 * t442 + t117 * t328 - t155 * t253 - t156 * t254) * mrSges(5,3) + (-t397 / 0.2e1 + t441) * t442 + (t208 * t258 - t209 * t259 + (-t246 * t301 - t247 * t302) * pkin(2)) * mrSges(4,3) - t76 * t136 - t75 * t137 + t170 * t32 + t171 * t33 + t185 * mrSges(4,1) - t187 * mrSges(4,2) - t131 * t191 + t259 * (Ifges(4,1) * t258 + t378) / 0.2e1 - t227 * t159 + t229 * t13 - t213 * t231 - t212 * t232 - Ifges(4,6) * t246 + Ifges(4,5) * t247 + t250 * t56 - qJD(2) * (Ifges(4,5) * t258 + Ifges(4,6) * t259) / 0.2e1 - t280 * (-mrSges(4,1) * t259 + mrSges(4,2) * t258) + ((t351 - t267 / 0.2e1 - t296 / 0.2e1 + qJD(1) * t426 + (t281 - t373) * pkin(7)) * t310 + (t350 + t266 / 0.2e1 + (t427 + t396 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t310) * qJD(1) + (t282 + t372) * pkin(7) + (-m(4) * t280 - t207) * pkin(2)) * t306) * qJD(1); t343 + t344 + (-t77 - t361) * t328 + t446 * t99 + t445 * t98 + t308 * t71 - t325 * t442 + t329 * qJD(5) - t258 * t231 - t259 * t232 - t326 * t32 + t275 * t33 + t304 * t72 + (t20 * t446 + t21 * t445 + t275 * t3 - t326 * t4 - t328 * t94) * m(7) + (-t112 * t328 - t200 * t331 - t333) * m(6) + (t116 * t328 - t117 * t442 + t217) * m(5) + (-t208 * t259 - t209 * t258 + t292) * m(4); (-pkin(4) * t61 - t112 * t117 - t68 * t82 - t69 * t83) * m(6) + (t117 * mrSges(5,3) + t314) * t328 + ((-t68 * mrSges(6,3) - pkin(9) * t137) * t308 + (-t69 * mrSges(6,3) + pkin(5) * t77 - pkin(9) * t136) * t304) * qJD(5) + t448 * t98 + t449 * t99 + (t223 * t4 + t224 * t3 + t294 * t37 + (t353 - t95) * t94 + t448 * t21 + t449 * t20) * m(7) - pkin(4) * t56 + t330 * pkin(9) + pkin(9) * t316 + (t116 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t328 + t441) * t442 - t95 * t77 - t83 * t136 - t82 * t137 - t361 * t117 - t116 * t191 + t223 * t32 + t224 * t33 + t294 * t13 + t478; -m(7) * (t20 * t24 + t21 * t25) + t342 * t430 + (t186 * t68 + t188 * t69) * mrSges(6,3) + (-Ifges(6,2) * t188 + t105 + t184) * t417 + (-t188 * t77 + t303 * t33 + t307 * t32 + (-t303 * t99 + t307 * t98) * qJD(6) + (-t188 * t94 + t3 * t303 + t307 * t4 + (-t20 * t303 + t21 * t307) * qJD(6)) * m(7)) * pkin(5) + t461 * t422 + (Ifges(6,5) * t186 - Ifges(6,6) * t188) * t411 + t104 * t415 + (Ifges(6,1) * t186 - t379) * t416 + t152 + t437 - t129 * t432 - t25 * t98 - t24 * t99 - t68 * t136 + t69 * t137 - t112 * (mrSges(6,1) * t188 + mrSges(6,2) * t186) + t469; t64 * t419 - t20 * t98 + t21 * t99 + (t461 + t65) * t422 + t462 + t469;];
tauc  = t1(:);
