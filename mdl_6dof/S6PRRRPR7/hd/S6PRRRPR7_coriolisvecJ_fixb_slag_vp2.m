% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:39:02
% EndTime: 2019-03-08 23:39:44
% DurationCPUTime: 20.60s
% Computational Cost: add. (13678->785), mult. (37195->1134), div. (0->0), fcn. (30181->14), ass. (0->350)
t295 = cos(pkin(7));
t303 = cos(qJ(3));
t304 = cos(qJ(2));
t369 = t303 * t304;
t299 = sin(qJ(3));
t300 = sin(qJ(2));
t374 = t299 * t300;
t314 = -t295 * t374 + t369;
t293 = sin(pkin(6));
t361 = qJD(1) * t293;
t218 = t314 * t361;
t292 = sin(pkin(7));
t381 = t292 * t299;
t284 = pkin(9) * t381;
t375 = t295 * t303;
t259 = pkin(2) * t375 - t284;
t250 = t259 * qJD(3);
t446 = -t218 + t250;
t376 = t295 * t299;
t380 = t292 * t303;
t260 = pkin(2) * t376 + pkin(9) * t380;
t241 = pkin(10) * t295 + t260;
t328 = -pkin(3) * t303 - pkin(10) * t299;
t242 = (-pkin(2) + t328) * t292;
t327 = pkin(3) * t299 - pkin(10) * t303;
t313 = t327 * qJD(3);
t249 = t292 * t313;
t298 = sin(qJ(4));
t302 = cos(qJ(4));
t346 = t300 * t361;
t337 = t292 * t346;
t354 = qJD(4) * t302;
t355 = qJD(4) * t298;
t447 = -t241 * t355 + t242 * t354 + t446 * t302 + (t249 - t337) * t298;
t372 = t300 * t303;
t373 = t299 * t304;
t316 = t295 * t372 + t373;
t445 = t260 * qJD(3) - t316 * t361;
t356 = qJD(3) * t299;
t469 = -(qJ(5) * t356 - qJD(5) * t303) * t292 - t447;
t257 = -t302 * t295 + t298 * t381;
t342 = qJD(3) * t380;
t208 = -qJD(4) * t257 + t302 * t342;
t258 = t295 * t298 + t302 * t381;
t209 = qJD(4) * t258 + t298 * t342;
t468 = pkin(4) * t209 - qJ(5) * t208 - qJD(5) * t258 + t445;
t320 = pkin(4) * t298 - qJ(5) * t302;
t253 = qJD(4) * t320 - qJD(5) * t298;
t291 = sin(pkin(13));
t294 = cos(pkin(13));
t351 = pkin(10) * t355;
t204 = t294 * t253 + t291 * t351;
t359 = qJD(2) * t292;
t262 = pkin(9) * t359 + t346;
t252 = t299 * t262;
t272 = qJD(2) * pkin(2) + t304 * t361;
t296 = cos(pkin(6));
t360 = qJD(1) * t296;
t312 = t272 * t295 + t292 * t360;
t171 = t303 * t312 - t252;
t248 = t327 * t359;
t131 = t302 * t171 + t298 * t248;
t345 = t299 * t359;
t119 = qJ(5) * t345 + t131;
t370 = t303 * t262;
t317 = t272 * t376 + t370;
t357 = qJD(2) * t303;
t134 = (t299 * t360 + t320 * t357) * t292 + t317;
t66 = -t119 * t291 + t294 * t134;
t467 = t204 - t66;
t453 = t469 * t291 + t468 * t294;
t452 = t468 * t291 - t469 * t294;
t371 = t302 * t303;
t216 = (t291 * t299 + t294 * t371) * t359;
t344 = t292 * t357;
t377 = t294 * t302;
t414 = pkin(5) * t298;
t466 = pkin(11) * t216 - t344 * t414 + (-pkin(11) * t377 + t414) * qJD(4) + t467;
t215 = (-t291 * t371 + t294 * t299) * t359;
t237 = t291 * t253;
t378 = t294 * t298;
t382 = t291 * t302;
t67 = t294 * t119 + t291 * t134;
t465 = -pkin(11) * t215 + t237 + (-pkin(10) * t378 - pkin(11) * t382) * qJD(4) - t67;
t343 = t292 * t356;
t181 = t208 * t294 + t291 * t343;
t464 = pkin(5) * t209 - pkin(11) * t181 + t453;
t180 = -t208 * t291 + t294 * t343;
t463 = -pkin(11) * t180 - t452;
t412 = -qJD(3) / 0.2e1;
t274 = -pkin(4) * t302 - qJ(5) * t298 - pkin(3);
t264 = t294 * t274;
t200 = -pkin(11) * t378 + t264 + (-pkin(10) * t291 - pkin(5)) * t302;
t229 = pkin(10) * t377 + t291 * t274;
t210 = -pkin(11) * t291 * t298 + t229;
t297 = sin(qJ(6));
t301 = cos(qJ(6));
t132 = t200 * t301 - t210 * t297;
t462 = qJD(6) * t132 + t297 * t466 + t301 * t465;
t133 = t200 * t297 + t210 * t301;
t461 = -qJD(6) * t133 - t297 * t465 + t301 * t466;
t130 = -t298 * t171 + t248 * t302;
t122 = -pkin(4) * t345 - t130;
t347 = pkin(5) * t291 + pkin(10);
t460 = pkin(5) * t215 + t347 * t354 - t122;
t358 = qJD(2) * t295;
t339 = qJD(3) + t358;
t444 = -t298 * t345 + t302 * t339;
t223 = Ifges(5,4) * t444;
t279 = qJD(4) - t344;
t389 = t279 * Ifges(5,5);
t231 = t298 * t339 + t302 * t345;
t391 = t231 * Ifges(5,1);
t142 = t223 + t389 + t391;
t153 = -pkin(3) * t339 - t171;
t172 = t299 * t312 + t370;
t154 = pkin(10) * t339 + t172;
t283 = t295 * t360;
t186 = t283 + (qJD(2) * t328 - t272) * t292;
t96 = -t298 * t154 + t186 * t302;
t459 = -t223 / 0.2e1 + t96 * mrSges(5,3) - t142 / 0.2e1 - t153 * mrSges(5,2) - t389 / 0.2e1;
t189 = t231 * t294 + t279 * t291;
t340 = -t231 * t291 + t294 * t279;
t458 = -t189 * t297 + t301 * t340;
t114 = t189 * t301 + t297 * t340;
t341 = qJD(3) * t359;
t331 = t303 * t341;
t195 = qJD(4) * t444 + t302 * t331;
t332 = t299 * t341;
t161 = -t195 * t291 + t294 * t332;
t162 = t195 * t294 + t291 * t332;
t46 = qJD(6) * t458 + t161 * t297 + t162 * t301;
t437 = t46 / 0.2e1;
t47 = -qJD(6) * t114 + t161 * t301 - t162 * t297;
t436 = t47 / 0.2e1;
t196 = qJD(4) * t231 + t298 * t331;
t426 = t196 / 0.2e1;
t457 = t345 / 0.2e1;
t456 = Ifges(4,6) * t412;
t207 = t258 * t294 - t291 * t380;
t240 = t284 + (-pkin(2) * t303 - pkin(3)) * t295;
t155 = pkin(4) * t257 - qJ(5) * t258 + t240;
t174 = t302 * t241 + t298 * t242;
t159 = -qJ(5) * t380 + t174;
t80 = t294 * t155 - t159 * t291;
t59 = pkin(5) * t257 - pkin(11) * t207 + t80;
t206 = -t258 * t291 - t294 * t380;
t81 = t291 * t155 + t294 * t159;
t71 = pkin(11) * t206 + t81;
t22 = t297 * t59 + t301 * t71;
t455 = -qJD(6) * t22 + t463 * t297 + t464 * t301;
t21 = -t297 * t71 + t301 * t59;
t454 = qJD(6) * t21 + t464 * t297 - t463 * t301;
t410 = pkin(11) + qJ(5);
t275 = t410 * t291;
t276 = t410 * t294;
t211 = -t275 * t301 - t276 * t297;
t318 = t291 * t297 - t294 * t301;
t169 = pkin(4) * t231 - qJ(5) * t444;
t64 = t294 * t169 - t291 * t96;
t39 = -pkin(11) * t294 * t444 + pkin(5) * t231 + t64;
t383 = t444 * t291;
t65 = t291 * t169 + t294 * t96;
t53 = -pkin(11) * t383 + t65;
t451 = -qJD(5) * t318 + qJD(6) * t211 - t297 * t39 - t301 * t53;
t212 = -t275 * t297 + t276 * t301;
t267 = t291 * t301 + t294 * t297;
t450 = -qJD(5) * t267 - qJD(6) * t212 + t297 * t53 - t301 * t39;
t363 = mrSges(4,1) * t339 + mrSges(5,1) * t444 - mrSges(5,2) * t231 - mrSges(4,3) * t345;
t443 = t295 * t369 - t374;
t255 = t318 * qJD(6);
t310 = t314 * qJD(2);
t334 = t296 * t342;
t123 = (t272 * t375 - t252) * qJD(3) + (t293 * t310 + t334) * qJD(1);
t201 = (t313 + t346) * t359;
t35 = t302 * t123 - t154 * t355 + t186 * t354 + t298 * t201;
t36 = -t298 * t123 - t154 * t354 - t186 * t355 + t201 * t302;
t442 = -t36 * mrSges(5,1) + t35 * mrSges(5,2) - Ifges(5,5) * t195 + Ifges(5,6) * t196;
t388 = t279 * Ifges(5,6);
t409 = Ifges(5,4) * t231;
t141 = Ifges(5,2) * t444 + t388 + t409;
t97 = t154 * t302 + t186 * t298;
t83 = qJ(5) * t279 + t97;
t98 = -pkin(4) * t444 - t231 * qJ(5) + t153;
t33 = -t291 * t83 + t294 * t98;
t34 = t291 * t98 + t294 * t83;
t224 = qJD(6) - t444;
t393 = t224 * Ifges(7,3);
t397 = t189 * Ifges(6,5);
t398 = t340 * Ifges(6,6);
t401 = t114 * Ifges(7,5);
t402 = t458 * Ifges(7,6);
t43 = t393 + t401 + t402;
t24 = -pkin(5) * t444 - pkin(11) * t189 + t33;
t26 = pkin(11) * t340 + t34;
t7 = t24 * t301 - t26 * t297;
t8 = t24 * t297 + t26 * t301;
t91 = -Ifges(6,3) * t444 + t397 + t398;
t306 = t388 / 0.2e1 - t398 / 0.2e1 - t397 / 0.2e1 - t393 / 0.2e1 - t402 / 0.2e1 - t401 / 0.2e1 - t91 / 0.2e1 - t43 / 0.2e1 + t141 / 0.2e1 + t97 * mrSges(5,3) + t8 * mrSges(7,2) - t7 * mrSges(7,1) + t34 * mrSges(6,2) - t33 * mrSges(6,1) - t153 * mrSges(5,1) + t409 / 0.2e1;
t350 = -Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1;
t441 = -t350 * t444 + t306;
t407 = Ifges(6,4) * t294;
t322 = -Ifges(6,2) * t291 + t407;
t408 = Ifges(6,4) * t291;
t323 = Ifges(6,1) * t294 - t408;
t324 = mrSges(6,1) * t291 + mrSges(6,2) * t294;
t415 = t294 / 0.2e1;
t416 = -t291 / 0.2e1;
t427 = t189 / 0.2e1;
t428 = t340 / 0.2e1;
t82 = -pkin(4) * t279 + qJD(5) - t96;
t92 = t189 * Ifges(6,4) + Ifges(6,2) * t340 - Ifges(6,6) * t444;
t93 = t189 * Ifges(6,1) + Ifges(6,4) * t340 - Ifges(6,5) * t444;
t440 = t322 * t428 + t323 * t427 + t82 * t324 + (-t291 * t34 - t294 * t33) * mrSges(6,3) + t92 * t416 + t93 * t415 - t459;
t439 = Ifges(7,4) * t437 + Ifges(7,2) * t436 + Ifges(7,6) * t426;
t438 = Ifges(7,1) * t437 + Ifges(7,4) * t436 + Ifges(7,5) * t426;
t74 = t162 * Ifges(6,1) + t161 * Ifges(6,4) + t196 * Ifges(6,5);
t435 = t74 / 0.2e1;
t434 = -t458 / 0.2e1;
t433 = t458 / 0.2e1;
t432 = -t114 / 0.2e1;
t431 = t114 / 0.2e1;
t430 = t161 / 0.2e1;
t429 = t162 / 0.2e1;
t423 = -t224 / 0.2e1;
t422 = t224 / 0.2e1;
t421 = t444 / 0.2e1;
t420 = -t444 / 0.2e1;
t419 = -t257 / 0.2e1;
t417 = t258 / 0.2e1;
t42 = Ifges(7,5) * t46;
t41 = Ifges(7,6) * t47;
t413 = pkin(10) * t298;
t411 = qJD(3) / 0.2e1;
t29 = qJ(5) * t332 + qJD(5) * t279 + t35;
t311 = t316 * qJD(2);
t335 = t296 * t343;
t124 = t317 * qJD(3) + (t293 * t311 + t335) * qJD(1);
t57 = pkin(4) * t196 - qJ(5) * t195 - qJD(5) * t231 + t124;
t15 = t294 * t29 + t291 * t57;
t406 = Ifges(7,4) * t114;
t405 = Ifges(6,5) * t294;
t404 = Ifges(4,6) * t295;
t403 = Ifges(6,6) * t291;
t400 = t161 * Ifges(6,6);
t399 = t162 * Ifges(6,5);
t396 = t195 * Ifges(5,1);
t395 = t195 * Ifges(5,4);
t394 = t196 * Ifges(5,4);
t167 = mrSges(5,1) * t332 - mrSges(5,3) * t195;
t87 = -t161 * mrSges(6,1) + t162 * mrSges(6,2);
t386 = -t167 + t87;
t202 = -t293 * t443 - t296 * t380;
t384 = t124 * t202;
t379 = t293 * t300;
t149 = t215 * t301 - t216 * t297;
t183 = t255 * t298 - t267 * t354;
t368 = t149 - t183;
t150 = t215 * t297 + t216 * t301;
t256 = t267 * qJD(6);
t182 = -t256 * t298 - t318 * t354;
t367 = t150 - t182;
t151 = t267 * t444;
t366 = -t151 + t256;
t152 = t318 * t444;
t365 = -t152 + t255;
t118 = -mrSges(6,1) * t340 + mrSges(6,2) * t189;
t199 = mrSges(5,1) * t279 - mrSges(5,3) * t231;
t364 = t199 - t118;
t11 = Ifges(7,3) * t196 + t41 + t42;
t58 = -mrSges(7,1) * t458 + mrSges(7,2) * t114;
t349 = -t58 + t364;
t16 = -t47 * mrSges(7,1) + t46 * mrSges(7,2);
t14 = -t29 * t291 + t294 * t57;
t173 = -t298 * t241 + t242 * t302;
t222 = -t272 * t292 + t283;
t281 = Ifges(4,4) * t344;
t338 = t222 * mrSges(4,2) + Ifges(4,1) * t457 + Ifges(4,5) * t339 / 0.2e1 + t281 / 0.2e1;
t336 = t359 * t379;
t10 = pkin(11) * t161 + t15;
t9 = pkin(5) * t196 - pkin(11) * t162 + t14;
t1 = qJD(6) * t7 + t10 * t301 + t297 * t9;
t2 = -qJD(6) * t8 - t10 * t297 + t301 * t9;
t330 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t160 = pkin(4) * t380 - t173;
t326 = -t124 * mrSges(4,1) - t123 * mrSges(4,2);
t325 = -mrSges(4,1) * t303 + mrSges(4,2) * t299;
t321 = -t403 + t405;
t319 = -t14 * t291 + t15 * t294;
t315 = t295 * t373 + t372;
t203 = t293 * t315 + t296 * t381;
t254 = -t292 * t293 * t304 + t295 * t296;
t157 = t203 * t302 + t254 * t298;
t103 = -t157 * t291 + t202 * t294;
t104 = t157 * t294 + t202 * t291;
t51 = t103 * t301 - t104 * t297;
t52 = t103 * t297 + t104 * t301;
t156 = t203 * t298 - t254 * t302;
t135 = t206 * t301 - t207 * t297;
t136 = t206 * t297 + t207 * t301;
t108 = -t241 * t354 - t242 * t355 + t249 * t302 - t298 * t250;
t101 = -pkin(4) * t343 - t108;
t32 = -pkin(4) * t332 - t36;
t308 = t222 * mrSges(4,1) + t96 * mrSges(5,1) + t279 * Ifges(5,3) + t231 * Ifges(5,5) + t444 * Ifges(5,6) + t456 - (t404 + (Ifges(4,4) * t299 + Ifges(4,2) * t303) * t292) * qJD(2) / 0.2e1 - t172 * mrSges(4,3) - t97 * mrSges(5,2);
t290 = -pkin(5) * t294 - pkin(4);
t278 = Ifges(4,5) * t331;
t277 = Ifges(5,3) * t332;
t270 = t347 * t298;
t247 = t325 * t359;
t246 = -mrSges(4,2) * t339 + mrSges(4,3) * t344;
t244 = t318 * t298;
t243 = t267 * t298;
t236 = (mrSges(4,1) * t299 + mrSges(4,2) * t303) * t341;
t228 = -pkin(10) * t382 + t264;
t205 = -t294 * t351 + t237;
t198 = -mrSges(5,2) * t279 + mrSges(5,3) * t444;
t184 = t218 * t298 - t302 * t337;
t168 = -mrSges(5,2) * t332 - mrSges(5,3) * t196;
t145 = t334 + (qJD(3) * t443 + t310) * t293;
t144 = t335 + (qJD(3) * t315 + t311) * t293;
t138 = -mrSges(6,1) * t444 - mrSges(6,3) * t189;
t137 = mrSges(6,2) * t444 + mrSges(6,3) * t340;
t125 = mrSges(5,1) * t196 + mrSges(5,2) * t195;
t116 = -pkin(5) * t206 + t160;
t111 = Ifges(5,5) * t332 - t394 + t396;
t110 = -t196 * Ifges(5,2) + Ifges(5,6) * t332 + t395;
t109 = Ifges(7,4) * t458;
t106 = mrSges(6,1) * t196 - mrSges(6,3) * t162;
t105 = -mrSges(6,2) * t196 + mrSges(6,3) * t161;
t85 = mrSges(7,1) * t224 - mrSges(7,3) * t114;
t84 = -mrSges(7,2) * t224 + mrSges(7,3) * t458;
t78 = pkin(5) * t383 + t97;
t77 = -qJD(4) * t156 + t145 * t302 + t298 * t336;
t76 = qJD(4) * t157 + t145 * t298 - t302 * t336;
t75 = -pkin(5) * t180 + t101;
t73 = t162 * Ifges(6,4) + t161 * Ifges(6,2) + t196 * Ifges(6,6);
t72 = t196 * Ifges(6,3) + t399 + t400;
t68 = -pkin(5) * t340 + t82;
t61 = -qJD(6) * t136 + t180 * t301 - t181 * t297;
t60 = qJD(6) * t135 + t180 * t297 + t181 * t301;
t49 = t144 * t291 + t294 * t77;
t48 = t144 * t294 - t291 * t77;
t45 = Ifges(7,1) * t114 + Ifges(7,5) * t224 + t109;
t44 = Ifges(7,2) * t458 + Ifges(7,6) * t224 + t406;
t31 = -mrSges(7,2) * t196 + mrSges(7,3) * t47;
t30 = mrSges(7,1) * t196 - mrSges(7,3) * t46;
t25 = -pkin(5) * t161 + t32;
t6 = -qJD(6) * t52 - t297 * t49 + t301 * t48;
t5 = qJD(6) * t51 + t297 * t48 + t301 * t49;
t3 = [t103 * t106 + t104 * t105 + t202 * t125 + t49 * t137 + t48 * t138 + t145 * t246 + t157 * t168 + t77 * t198 + t254 * t236 + t51 * t30 + t52 * t31 + t5 * t84 + t6 * t85 + (-mrSges(3,1) * t300 - mrSges(3,2) * t304) * qJD(2) ^ 2 * t293 - t363 * t144 - t349 * t76 + (t16 + t386) * t156 + (t247 * t379 + (t202 * t303 - t203 * t299) * qJD(3) * mrSges(4,3)) * t359 + m(4) * (t123 * t203 + t384 - t144 * t171 + t145 * t172 + (qJD(1) * t254 + t222) * t336) + m(5) * (t144 * t153 - t156 * t36 + t157 * t35 - t76 * t96 + t77 * t97 + t384) + m(6) * (t103 * t14 + t104 * t15 + t156 * t32 + t33 * t48 + t34 * t49 + t76 * t82) + m(7) * (t1 * t52 + t156 * t25 + t2 * t51 + t5 * t8 + t6 * t7 + t68 * t76); t279 * (Ifges(5,5) * t208 - Ifges(5,6) * t209) / 0.2e1 - t196 * (Ifges(5,4) * t258 - Ifges(5,2) * t257) / 0.2e1 + t231 * (Ifges(5,1) * t208 - Ifges(5,4) * t209) / 0.2e1 + t195 * (Ifges(5,1) * t258 - Ifges(5,4) * t257) / 0.2e1 + (Ifges(6,1) * t181 + Ifges(6,4) * t180 + Ifges(6,5) * t209) * t427 + (Ifges(6,4) * t181 + Ifges(6,2) * t180 + Ifges(6,6) * t209) * t428 + (Ifges(6,1) * t207 + Ifges(6,4) * t206 + Ifges(6,5) * t257) * t429 + t124 * (mrSges(5,1) * t257 + mrSges(5,2) * t258) + t14 * (mrSges(6,1) * t257 - mrSges(6,3) * t207) + t15 * (-mrSges(6,2) * t257 + mrSges(6,3) * t206) + t1 * (-mrSges(7,2) * t257 + mrSges(7,3) * t135) + t2 * (mrSges(7,1) * t257 - mrSges(7,3) * t136) + t240 * t125 + t32 * (-mrSges(6,1) * t206 + mrSges(6,2) * t207) + t208 * t142 / 0.2e1 + t33 * (mrSges(6,1) * t209 - mrSges(6,3) * t181) + t34 * (-mrSges(6,2) * t209 + mrSges(6,3) * t180) + t7 * (mrSges(7,1) * t209 - mrSges(7,3) * t60) + t8 * (-mrSges(7,2) * t209 + mrSges(7,3) * t61) + t153 * (mrSges(5,1) * t209 + mrSges(5,2) * t208) - t209 * t141 / 0.2e1 + t206 * t73 / 0.2e1 + t108 * t199 + t82 * (-mrSges(6,1) * t180 + mrSges(6,2) * t181) + t181 * t93 / 0.2e1 + t173 * t167 + t174 * t168 + t180 * t92 / 0.2e1 + t160 * t87 + t25 * (-mrSges(7,1) * t135 + mrSges(7,2) * t136) + t116 * t16 + t101 * t118 + t81 * t105 + t80 * t106 + t75 * t58 + t68 * (-mrSges(7,1) * t61 + mrSges(7,2) * t60) + t60 * t45 / 0.2e1 + t61 * t44 / 0.2e1 + t21 * t30 + t22 * t31 + t111 * t417 + t110 * t419 + (Ifges(6,5) * t181 + Ifges(6,6) * t180 + Ifges(6,3) * t209) * t420 + (Ifges(5,4) * t208 - Ifges(5,2) * t209) * t421 + (Ifges(7,5) * t60 + Ifges(7,6) * t61 + Ifges(7,3) * t209) * t422 + (t72 + t11) * t257 / 0.2e1 + (t91 + t43) * t209 / 0.2e1 + (Ifges(6,4) * t207 + Ifges(6,2) * t206 + Ifges(6,6) * t257) * t430 + (Ifges(7,1) * t60 + Ifges(7,4) * t61 + Ifges(7,5) * t209) * t431 + (Ifges(7,4) * t60 + Ifges(7,2) * t61 + Ifges(7,6) * t209) * t433 + t207 * t435 + (Ifges(7,4) * t136 + Ifges(7,2) * t135 + Ifges(7,6) * t257) * t436 + (Ifges(7,1) * t136 + Ifges(7,4) * t135 + Ifges(7,5) * t257) * t437 + (t278 / 0.2e1 + t326) * t295 + (Ifges(6,5) * t207 + Ifges(7,5) * t136 + Ifges(6,6) * t206 + Ifges(7,6) * t135 + (Ifges(6,3) + Ifges(7,3)) * t257) * t426 + t136 * t438 - t363 * t445 + t446 * t246 + (t123 * t260 - t124 * t259 - t171 * t445 + t172 * t446) * m(4) + t447 * t198 + (t124 * t240 + t173 * t36 + t174 * t35 + t447 * t97 + (t108 + t184) * t96 + t445 * t153) * m(5) + t452 * t137 + t453 * t138 + (t14 * t80 + t15 * t81 + t160 * t32 + (t101 - t184) * t82 + t452 * t34 + t453 * t33) * m(6) + t454 * t84 + t455 * t85 + (t1 * t22 + t116 * t25 + t2 * t21 + t454 * t8 + t455 * t7 + (-t184 + t75) * t68) * m(7) + (t124 * mrSges(4,3) * t299 - pkin(2) * t236 + (t123 * mrSges(4,3) - t277 / 0.2e1 + t442) * t303 + (-m(4) * t222 - t247 + (-m(4) * pkin(2) + t325) * t359) * t346 + (((t358 + t411) * Ifges(4,5) + (-qJD(2) * t259 - t171) * mrSges(4,3) + t338) * t303 + (t456 + (-t260 * mrSges(4,3) - 0.3e1 / 0.2e1 * t404 + Ifges(5,5) * t417 + Ifges(5,6) * t419 + (-Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t380) * qJD(2) + t308) * t299) * qJD(3) + (0.3e1 / 0.2e1 * t303 ^ 2 - 0.3e1 / 0.2e1 * t299 ^ 2) * Ifges(4,4) * t341) * t292 + (-t208 * t96 - t209 * t97 - t257 * t35 - t258 * t36) * mrSges(5,3) + t135 * t439 + t349 * t184; t461 * t85 + (t1 * t133 + t132 * t2 + t25 * t270 + t460 * t68 + t461 * t7 + t462 * t8) * m(7) + t462 * t84 + t460 * t58 + (((Ifges(5,5) * t298 + Ifges(5,6) * t302) * t411 + Ifges(4,4) * t457 + (t412 + t358 / 0.2e1) * Ifges(4,6) - t308) * t299 + (-t281 / 0.2e1 + t171 * mrSges(4,3) + (-t295 * Ifges(4,5) / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t381) * qJD(2) + Ifges(4,5) * t412 + (-t391 / 0.2e1 + t459) * t302 + t441 * t298 - t338) * t303) * t359 + t467 * t138 + t278 + (t205 - t67) * t137 - t189 * (Ifges(6,1) * t216 + Ifges(6,4) * t215) / 0.2e1 + (-t215 * t34 + t216 * t33) * mrSges(6,3) + (t182 / 0.2e1 - t150 / 0.2e1) * t45 + t326 + t270 * t16 - t171 * t246 + t228 * t106 + t229 * t105 - t215 * t92 / 0.2e1 - t82 * (-mrSges(6,1) * t215 + mrSges(6,2) * t216) - t216 * t93 / 0.2e1 - t131 * t198 - t130 * t199 + t132 * t30 + t133 * t31 - pkin(3) * t125 - t122 * t118 + (Ifges(6,5) * t216 + Ifges(6,6) * t215) * t421 + (Ifges(7,5) * t182 + Ifges(7,6) * t183) * t422 + (Ifges(7,5) * t150 + Ifges(7,6) * t149) * t423 + (-Ifges(7,1) * t244 - Ifges(7,4) * t243) * t437 + (-t1 * t243 + t2 * t244 + t367 * t7 - t368 * t8) * mrSges(7,3) + t25 * (mrSges(7,1) * t243 - mrSges(7,2) * t244) + (-Ifges(7,5) * t244 - Ifges(7,6) * t243) * t426 + (-Ifges(7,4) * t244 - Ifges(7,2) * t243) * t436 - t340 * (Ifges(6,4) * t216 + Ifges(6,2) * t215) / 0.2e1 + (t73 * t416 - t36 * mrSges(5,3) + t74 * t415 + t111 / 0.2e1 + t32 * t324 + t396 / 0.2e1 - t394 / 0.2e1 + t124 * mrSges(5,2) + t321 * t426 + t323 * t429 + t322 * t430 + t386 * pkin(10) + (-t14 * t294 - t15 * t291) * mrSges(6,3)) * t298 + (t35 * mrSges(5,3) + pkin(10) * t168 - t42 / 0.2e1 - t41 / 0.2e1 + t110 / 0.2e1 - t72 / 0.2e1 - t11 / 0.2e1 + t15 * mrSges(6,2) + t395 / 0.2e1 - t124 * mrSges(5,1) - t399 / 0.2e1 - t400 / 0.2e1 - t14 * mrSges(6,1) + (-Ifges(7,3) / 0.2e1 + t350) * t196 - t330) * t302 + m(5) * (pkin(10) * t302 * t35 - pkin(3) * t124 - t36 * t413) + m(6) * (t14 * t228 + t15 * t229 + t204 * t33 + t205 * t34 + t32 * t413) - m(5) * (t130 * t96 + t131 * t97 + t153 * t172) - m(6) * (t122 * t82 + t33 * t66 + t34 * t67) + (Ifges(7,1) * t182 + Ifges(7,4) * t183) * t431 + (Ifges(7,1) * t150 + Ifges(7,4) * t149) * t432 + (Ifges(7,4) * t182 + Ifges(7,2) * t183) * t433 + (Ifges(7,4) * t150 + Ifges(7,2) * t149) * t434 + (((-m(5) * t97 - t198) * pkin(10) - t441) * t298 + (t391 / 0.2e1 + t321 * t420 + (-m(5) * t96 + m(6) * t82 - t364) * pkin(10) + t440) * t302) * qJD(4) - t244 * t438 + (t183 / 0.2e1 - t149 / 0.2e1) * t44 - t243 * t439 + t363 * t172 + (mrSges(7,1) * t368 - mrSges(7,2) * t367) * t68; (-Ifges(7,4) * t152 - Ifges(7,2) * t151) * t434 + (-t255 / 0.2e1 + t152 / 0.2e1) * t45 - (-(t405 / 0.2e1 - t403 / 0.2e1) * t444 + (Ifges(5,1) / 0.2e1 + t350) * t231 + t440) * t444 + (-t256 / 0.2e1 + t151 / 0.2e1) * t44 + (-Ifges(7,5) * t152 - Ifges(7,6) * t151) * t423 + (-Ifges(7,1) * t152 - Ifges(7,4) * t151) * t432 + (-t1 * t318 - t2 * t267 + t365 * t7 - t366 * t8) * mrSges(7,3) + t25 * (mrSges(7,1) * t318 + mrSges(7,2) * t267) + (Ifges(7,4) * t267 - Ifges(7,2) * t318) * t436 + (Ifges(7,1) * t267 - Ifges(7,4) * t318) * t437 + (Ifges(6,5) * t291 + Ifges(7,5) * t267 + Ifges(6,6) * t294 - Ifges(7,6) * t318) * t426 - t318 * t439 + t277 + (t137 * t294 - t138 * t291) * qJD(5) + (-pkin(4) * t32 + (-t291 * t33 + t294 * t34) * qJD(5) + t319 * qJ(5) - t33 * t64 - t34 * t65 - t82 * t97) * m(6) + (t105 * t294 - t106 * t291) * qJ(5) + t306 * t231 - t442 + t290 * t16 + t211 * t30 + t212 * t31 - t96 * t198 - t65 * t137 - t64 * t138 - pkin(4) * t87 - t78 * t58 + t73 * t415 + (-Ifges(7,5) * t255 - Ifges(7,6) * t256) * t422 + (-Ifges(7,1) * t255 - Ifges(7,4) * t256) * t431 + (-Ifges(7,4) * t255 - Ifges(7,2) * t256) * t433 + (Ifges(6,1) * t291 + t407) * t429 + (Ifges(6,2) * t294 + t408) * t430 + t291 * t435 + t319 * mrSges(6,3) + t267 * t438 + t32 * (-mrSges(6,1) * t294 + mrSges(6,2) * t291) + t450 * t85 + t451 * t84 + (t1 * t212 + t2 * t211 + t25 * t290 + t450 * t7 + t451 * t8 - t68 * t78) * m(7) + t364 * t97 + (mrSges(7,1) * t366 - mrSges(7,2) * t365) * t68; t114 * t85 - t458 * t84 - t340 * t137 + t189 * t138 + t16 + t87 + (t114 * t7 - t458 * t8 + t25) * m(7) + (t189 * t33 - t34 * t340 + t32) * m(6); -t68 * (mrSges(7,1) * t114 + mrSges(7,2) * t458) + (Ifges(7,1) * t458 - t406) * t432 + t44 * t431 + (Ifges(7,5) * t458 - Ifges(7,6) * t114) * t423 - t7 * t84 + t8 * t85 + (t114 * t8 + t458 * t7) * mrSges(7,3) + t330 + t11 + (-Ifges(7,2) * t114 + t109 + t45) * t434;];
tauc  = t3(:);
