% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:47
% EndTime: 2019-03-09 15:48:50
% DurationCPUTime: 33.53s
% Computational Cost: add. (14625->820), mult. (38987->1091), div. (0->0), fcn. (30322->10), ass. (0->343)
t292 = cos(qJ(2));
t289 = sin(qJ(2));
t284 = sin(pkin(6));
t360 = qJD(1) * t284;
t343 = t289 * t360;
t286 = cos(pkin(6));
t359 = qJD(1) * t286;
t352 = pkin(1) * t359;
t238 = -pkin(8) * t343 + t292 * t352;
t301 = (pkin(2) * t289 - pkin(9) * t292) * t284;
t239 = qJD(1) * t301;
t288 = sin(qJ(3));
t291 = cos(qJ(3));
t176 = -t288 * t238 + t291 * t239;
t387 = -qJ(4) - pkin(9);
t333 = qJD(3) * t387;
t470 = -(-qJ(4) * t291 * t292 + pkin(3) * t289) * t360 - t176 - qJD(4) * t288 + t291 * t333;
t177 = t291 * t238 + t288 * t239;
t342 = t292 * t360;
t329 = t288 * t342;
t469 = -qJ(4) * t329 - qJD(4) * t291 - t288 * t333 + t177;
t283 = sin(pkin(11));
t285 = cos(pkin(11));
t368 = t285 * t291;
t201 = -t283 * t329 + t342 * t368;
t355 = qJD(3) * t291;
t356 = qJD(3) * t288;
t247 = -t283 * t356 + t285 * t355;
t361 = t201 - t247;
t241 = pkin(8) * t342 + t289 * t352;
t198 = pkin(3) * t329 + t241;
t468 = pkin(3) * t356 - t198;
t253 = t283 * t291 + t285 * t288;
t200 = t253 * t342;
t246 = t253 * qJD(3);
t362 = t200 - t246;
t460 = t469 * t283 + t285 * t470;
t467 = qJ(5) * t361 - qJD(5) * t253 + t468;
t421 = pkin(4) + pkin(10);
t466 = -t362 * t421 + t467;
t371 = t284 * t289;
t330 = t421 * t371;
t465 = -pkin(5) * t361 + qJD(1) * t330 - t460;
t457 = t283 * t470 - t469 * t285;
t458 = qJ(5) * t343 - t457;
t464 = Ifges(5,5) / 0.2e1;
t272 = qJD(2) + t359;
t357 = qJD(2) * t292;
t339 = t291 * t357;
t190 = t272 * t355 + (-t289 * t356 + t339) * t360;
t340 = t288 * t357;
t191 = -t272 * t356 + (-t289 * t355 - t340) * t360;
t129 = t190 * t283 - t285 * t191;
t444 = -t129 / 0.2e1;
t130 = t190 * t285 + t191 * t283;
t418 = t130 / 0.2e1;
t439 = -Ifges(6,4) + Ifges(5,5);
t438 = Ifges(6,5) - Ifges(5,6);
t252 = t283 * t288 - t368;
t281 = -pkin(3) * t291 - pkin(2);
t303 = -qJ(5) * t253 + t281;
t169 = t252 * t421 + t303;
t263 = t387 * t288;
t264 = t387 * t291;
t205 = -t285 * t263 - t264 * t283;
t180 = pkin(5) * t253 + t205;
t287 = sin(qJ(6));
t290 = cos(qJ(6));
t105 = t169 * t290 + t180 * t287;
t463 = -qJD(6) * t105 - t287 * t466 + t465 * t290;
t104 = -t169 * t287 + t180 * t290;
t462 = qJD(6) * t104 + t465 * t287 + t290 * t466;
t461 = pkin(5) * t362 - t458;
t459 = pkin(4) * t343 - t460;
t456 = -pkin(4) * t362 + t467;
t370 = t284 * t292;
t251 = t286 * t289 * pkin(1) + pkin(8) * t370;
t354 = qJD(1) * qJD(2);
t230 = t251 * t354;
t154 = -t191 * pkin(3) + t230;
t208 = pkin(9) * t272 + t241;
t234 = (-pkin(2) * t292 - pkin(9) * t289 - pkin(1)) * t284;
t215 = qJD(1) * t234;
t158 = t208 * t291 + t215 * t288;
t240 = qJD(2) * t301;
t228 = qJD(1) * t240;
t273 = pkin(8) * t371;
t392 = pkin(1) * t292;
t250 = t286 * t392 - t273;
t242 = t250 * qJD(2);
t229 = qJD(1) * t242;
t101 = -qJD(3) * t158 + t291 * t228 - t229 * t288;
t221 = t272 * t288 + t291 * t343;
t335 = t284 * t354;
t328 = t289 * t335;
t53 = pkin(3) * t328 - qJ(4) * t190 - qJD(4) * t221 + t101;
t100 = -t208 * t356 + t215 * t355 + t288 * t228 + t291 * t229;
t220 = t272 * t291 - t288 * t343;
t63 = qJ(4) * t191 + qJD(4) * t220 + t100;
t20 = -t283 * t63 + t285 * t53;
t19 = -pkin(4) * t328 - t20;
t304 = t220 * t283 + t285 * t221;
t294 = -t130 * qJ(5) - qJD(5) * t304 + t154;
t31 = t129 * pkin(4) + t294;
t22 = t129 * t421 + t294;
t262 = -qJD(3) + t342;
t157 = -t208 * t288 + t291 * t215;
t127 = -qJ(4) * t221 + t157;
t111 = -pkin(3) * t262 + t127;
t128 = qJ(4) * t220 + t158;
t117 = t283 * t128;
t51 = t111 * t285 - t117;
t320 = qJD(5) - t51;
t440 = pkin(5) * t304;
t32 = t262 * t421 + t320 + t440;
t166 = -t285 * t220 + t221 * t283;
t207 = -t272 * pkin(2) - t238;
t164 = -t220 * pkin(3) + qJD(4) + t207;
t295 = -qJ(5) * t304 + t164;
t41 = t166 * t421 + t295;
t5 = -t287 * t41 + t290 * t32;
t307 = qJD(2) * t330;
t8 = pkin(5) * t130 - qJD(1) * t307 - t20;
t1 = qJD(6) * t5 + t22 * t290 + t287 * t8;
t6 = t287 * t32 + t290 * t41;
t2 = -qJD(6) * t6 - t22 * t287 + t290 * t8;
t325 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t441 = -t328 / 0.2e1;
t137 = t166 * t287 - t262 * t290;
t67 = -qJD(6) * t137 + t129 * t290 - t287 * t328;
t64 = Ifges(7,6) * t67;
t136 = t166 * t290 + t262 * t287;
t66 = qJD(6) * t136 + t129 * t287 + t290 * t328;
t65 = Ifges(7,5) * t66;
t9 = Ifges(7,3) * t130 + t64 + t65;
t455 = t325 + Ifges(6,4) * t441 + t154 * mrSges(5,2) + t19 * mrSges(6,1) + t328 * t464 + t9 / 0.2e1 - t20 * mrSges(5,3) - t31 * mrSges(6,3) + (Ifges(6,6) + Ifges(5,4)) * t444 + (Ifges(6,2) + Ifges(5,1)) * t418;
t454 = Ifges(6,1) + Ifges(4,3) + Ifges(5,3);
t453 = pkin(5) * t166;
t163 = qJD(6) + t304;
t43 = Ifges(7,5) * t137 + Ifges(7,6) * t136 + Ifges(7,3) * t163;
t46 = pkin(4) * t262 + t320;
t72 = t166 * pkin(4) + t295;
t95 = -Ifges(6,4) * t262 - Ifges(6,2) * t304 + Ifges(6,6) * t166;
t98 = Ifges(5,1) * t304 - Ifges(5,4) * t166 - Ifges(5,5) * t262;
t452 = t164 * mrSges(5,2) + t46 * mrSges(6,1) + t5 * mrSges(7,1) - t95 / 0.2e1 - t51 * mrSges(5,3) - t6 * mrSges(7,2) - t72 * mrSges(6,3) + t98 / 0.2e1 + t43 / 0.2e1;
t442 = -t130 / 0.2e1;
t443 = t129 / 0.2e1;
t451 = Ifges(6,6) * t442 + Ifges(6,3) * t443;
t21 = t283 * t53 + t285 * t63;
t16 = -qJ(5) * t328 + qJD(5) * t262 - t21;
t450 = -t101 * mrSges(4,1) - t20 * mrSges(5,1) + t100 * mrSges(4,2) + t21 * mrSges(5,2) - t19 * mrSges(6,2) + t16 * mrSges(6,3);
t447 = Ifges(6,5) / 0.2e1;
t449 = t154 * mrSges(5,1) + t16 * mrSges(6,1) - t31 * mrSges(6,2) - t21 * mrSges(5,3) + Ifges(5,4) * t442 + Ifges(5,2) * t443 + Ifges(5,6) * t441 + t328 * t447 + t451;
t233 = pkin(9) * t286 + t251;
t175 = t291 * t233 + t288 * t234;
t434 = t1 * t287 + t2 * t290;
t433 = Ifges(4,5) * t190 + Ifges(4,6) * t191 + t129 * t438 + t130 * t439 + t328 * t454;
t376 = t221 * Ifges(4,4);
t152 = t220 * Ifges(4,2) - t262 * Ifges(4,6) + t376;
t216 = Ifges(4,4) * t220;
t153 = t221 * Ifges(4,1) - t262 * Ifges(4,5) + t216;
t305 = t157 * t291 + t158 * t288;
t384 = Ifges(4,4) * t291;
t385 = Ifges(4,4) * t288;
t393 = t291 / 0.2e1;
t397 = -t262 / 0.2e1;
t400 = t221 / 0.2e1;
t401 = t220 / 0.2e1;
t432 = -t305 * mrSges(4,3) + t207 * (mrSges(4,1) * t288 + mrSges(4,2) * t291) + (-Ifges(4,2) * t288 + t384) * t401 + (Ifges(4,1) * t291 - t385) * t400 + (Ifges(4,5) * t291 - Ifges(4,6) * t288) * t397 - t288 * t152 / 0.2e1 + t153 * t393;
t113 = -qJD(3) * t175 + t291 * t240 - t242 * t288;
t248 = t286 * t291 - t288 * t371;
t197 = qJD(3) * t248 + t284 * t339;
t249 = t286 * t288 + t291 * t371;
t358 = qJD(2) * t289;
t341 = t284 * t358;
t76 = pkin(3) * t341 - qJ(4) * t197 - qJD(4) * t249 + t113;
t112 = -t233 * t356 + t234 * t355 + t288 * t240 + t291 * t242;
t196 = -qJD(3) * t249 - t284 * t340;
t81 = qJ(4) * t196 + qJD(4) * t248 + t112;
t27 = t283 * t76 + t285 * t81;
t24 = -t284 * (qJ(5) * t358 - qJD(5) * t292) - t27;
t369 = t285 * t128;
t52 = t283 * t111 + t369;
t47 = qJ(5) * t262 - t52;
t93 = -Ifges(6,5) * t262 - Ifges(6,6) * t304 + Ifges(6,3) * t166;
t96 = Ifges(5,4) * t304 - Ifges(5,2) * t166 - Ifges(5,6) * t262;
t431 = t164 * mrSges(5,1) + t47 * mrSges(6,1) + t93 / 0.2e1 - t96 / 0.2e1 - t52 * mrSges(5,3) - t72 * mrSges(6,2);
t382 = Ifges(7,4) * t287;
t314 = Ifges(7,2) * t290 + t382;
t381 = Ifges(7,4) * t290;
t316 = Ifges(7,1) * t287 + t381;
t319 = mrSges(7,1) * t290 - mrSges(7,2) * t287;
t321 = t287 * t5 - t290 * t6;
t35 = -t47 - t453;
t378 = Ifges(7,6) * t290;
t379 = Ifges(7,5) * t287;
t395 = -t287 / 0.2e1;
t412 = -t163 / 0.2e1;
t415 = -t137 / 0.2e1;
t417 = -t136 / 0.2e1;
t383 = Ifges(7,4) * t137;
t44 = Ifges(7,2) * t136 + Ifges(7,6) * t163 + t383;
t135 = Ifges(7,4) * t136;
t45 = Ifges(7,1) * t137 + Ifges(7,5) * t163 + t135;
t430 = t321 * mrSges(7,3) + (t378 + t379) * t412 + t314 * t417 + t316 * t415 + t35 * t319 - t290 * t44 / 0.2e1 + t45 * t395;
t386 = Ifges(3,4) * t289;
t428 = t157 * mrSges(4,1) + t46 * mrSges(6,2) + t51 * mrSges(5,1) - Ifges(3,6) * t272 / 0.2e1 - (Ifges(3,2) * t292 + t386) * t360 / 0.2e1 - t158 * mrSges(4,2) - t47 * mrSges(6,3) - t52 * mrSges(5,2);
t11 = t66 * Ifges(7,1) + t67 * Ifges(7,4) + t130 * Ifges(7,5);
t427 = t11 / 0.2e1;
t423 = t66 / 0.2e1;
t422 = t67 / 0.2e1;
t420 = pkin(1) * mrSges(3,1);
t419 = pkin(1) * mrSges(3,2);
t416 = t136 / 0.2e1;
t414 = t137 / 0.2e1;
t411 = t163 / 0.2e1;
t410 = -t166 / 0.2e1;
t409 = t166 / 0.2e1;
t407 = -t304 / 0.2e1;
t406 = t304 / 0.2e1;
t404 = t190 / 0.2e1;
t403 = t191 / 0.2e1;
t399 = t248 / 0.2e1;
t398 = t249 / 0.2e1;
t396 = t262 / 0.2e1;
t394 = t290 / 0.2e1;
t391 = pkin(3) * t221;
t388 = qJD(2) / 0.2e1;
t380 = Ifges(3,5) * t292;
t377 = t220 * Ifges(4,6);
t375 = t221 * Ifges(4,5);
t374 = t229 * mrSges(3,2);
t373 = t272 * Ifges(3,5);
t143 = mrSges(6,1) * t166 + mrSges(6,3) * t262;
t82 = -mrSges(7,1) * t136 + mrSges(7,2) * t137;
t372 = -t143 + t82;
t367 = t287 * t246;
t366 = t290 * t246;
t174 = -t288 * t233 + t291 * t234;
t134 = -pkin(3) * t370 - t249 * qJ(4) + t174;
t147 = qJ(4) * t248 + t175;
t80 = t283 * t134 + t285 * t147;
t145 = mrSges(5,2) * t262 - mrSges(5,3) * t166;
t365 = -t145 + t143;
t144 = mrSges(6,1) * t304 - mrSges(6,2) * t262;
t146 = -mrSges(5,1) * t262 - mrSges(5,3) * t304;
t364 = t146 - t144;
t363 = -mrSges(3,1) * t272 - mrSges(4,1) * t220 + mrSges(4,2) * t221 + mrSges(3,3) * t343;
t243 = t251 * qJD(2);
t351 = t464 - Ifges(6,4) / 0.2e1;
t350 = -Ifges(5,6) / 0.2e1 + t447;
t349 = -Ifges(6,6) / 0.2e1 - Ifges(5,4) / 0.2e1;
t280 = -pkin(3) * t285 - pkin(4);
t26 = -t283 * t81 + t285 * t76;
t108 = t130 * mrSges(6,1) + mrSges(6,2) * t328;
t61 = t127 * t283 + t369;
t79 = t285 * t134 - t283 * t147;
t182 = t200 * t290 - t287 * t343;
t332 = -t182 + t366;
t183 = t200 * t287 + t290 * t343;
t331 = -t183 + t367;
t326 = Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1;
t172 = -pkin(3) * t196 + t243;
t324 = t1 * t290 - t2 * t287;
t322 = t287 * t6 + t290 * t5;
t75 = pkin(4) * t370 - t79;
t318 = mrSges(7,1) * t287 + mrSges(7,2) * t290;
t317 = Ifges(7,1) * t290 - t382;
t315 = -Ifges(7,2) * t287 + t381;
t313 = Ifges(7,5) * t290 - Ifges(7,6) * t287;
t311 = qJ(5) * t166 + t391;
t33 = mrSges(7,1) * t130 - mrSges(7,3) * t66;
t34 = -mrSges(7,2) * t130 + mrSges(7,3) * t67;
t310 = t287 * t34 + t290 * t33;
t185 = t248 * t283 + t249 * t285;
t42 = t185 * pkin(5) + pkin(10) * t370 + t75;
t184 = -t285 * t248 + t249 * t283;
t232 = t273 + (-pkin(2) - t392) * t286;
t189 = -t248 * pkin(3) + t232;
t296 = -t185 * qJ(5) + t189;
t71 = t184 * t421 + t296;
t17 = -t287 * t71 + t290 * t42;
t18 = t287 * t42 + t290 * t71;
t90 = -mrSges(7,2) * t163 + mrSges(7,3) * t136;
t91 = mrSges(7,1) * t163 - mrSges(7,3) * t137;
t309 = -t287 * t91 + t290 * t90;
t308 = -t287 * t90 - t290 * t91;
t306 = t100 * t291 - t101 * t288;
t62 = t127 * t285 - t117;
t206 = t263 * t283 - t264 * t285;
t74 = qJ(5) * t370 - t80;
t300 = -t287 * t184 + t290 * t370;
t161 = t290 * t184 + t287 * t370;
t142 = t196 * t283 + t197 * t285;
t297 = -qJ(5) * t142 - qJD(5) * t185 + t172;
t278 = pkin(3) * t283 + qJ(5);
t265 = Ifges(3,4) * t342;
t260 = t335 * t380;
t237 = -t272 * mrSges(3,2) + mrSges(3,3) * t342;
t203 = Ifges(3,1) * t343 + t265 + t373;
t194 = -mrSges(4,1) * t262 - mrSges(4,3) * t221;
t193 = mrSges(4,2) * t262 + mrSges(4,3) * t220;
t192 = pkin(4) * t252 + t303;
t181 = -pkin(5) * t252 + t206;
t171 = -mrSges(4,2) * t328 + mrSges(4,3) * t191;
t170 = mrSges(4,1) * t328 - mrSges(4,3) * t190;
t151 = -t262 * Ifges(4,3) + t375 + t377;
t141 = -t285 * t196 + t197 * t283;
t131 = -mrSges(4,1) * t191 + mrSges(4,2) * t190;
t121 = t130 * mrSges(6,3);
t120 = t130 * mrSges(5,2);
t116 = t190 * Ifges(4,1) + t191 * Ifges(4,4) + Ifges(4,5) * t328;
t115 = t190 * Ifges(4,4) + t191 * Ifges(4,2) + Ifges(4,6) * t328;
t110 = mrSges(5,1) * t328 - mrSges(5,3) * t130;
t109 = -mrSges(5,2) * t328 - mrSges(5,3) * t129;
t107 = mrSges(6,1) * t129 - mrSges(6,3) * t328;
t103 = -mrSges(6,2) * t166 - mrSges(6,3) * t304;
t102 = mrSges(5,1) * t166 + mrSges(5,2) * t304;
t97 = -t262 * Ifges(6,1) - Ifges(6,4) * t304 + t166 * Ifges(6,5);
t94 = Ifges(5,5) * t304 - t166 * Ifges(5,6) - t262 * Ifges(5,3);
t92 = t184 * pkin(4) + t296;
t87 = pkin(4) * t304 + t311;
t84 = qJD(6) * t161 + t287 * t141 + t290 * t341;
t83 = qJD(6) * t300 + t290 * t141 - t287 * t341;
t69 = -t129 * mrSges(6,2) - t121;
t68 = t129 * mrSges(5,1) + t120;
t50 = t304 * t421 + t311;
t48 = -t184 * pkin(5) - t74;
t40 = t62 - t440;
t39 = t61 - t453;
t38 = pkin(4) * t141 + t297;
t30 = t141 * t421 + t297;
t25 = -pkin(4) * t341 - t26;
t23 = -mrSges(7,1) * t67 + mrSges(7,2) * t66;
t15 = -t141 * pkin(5) - t24;
t14 = pkin(5) * t142 - t26 - t307;
t13 = t287 * t39 + t290 * t50;
t12 = -t287 * t50 + t290 * t39;
t10 = t66 * Ifges(7,4) + t67 * Ifges(7,2) + t130 * Ifges(7,6);
t7 = -pkin(5) * t129 - t16;
t4 = -qJD(6) * t18 + t14 * t290 - t287 * t30;
t3 = qJD(6) * t17 + t14 * t287 + t290 * t30;
t28 = [(t100 * t248 - t101 * t249 - t157 * t197 + t158 * t196) * mrSges(4,3) + (Ifges(4,5) * t197 + Ifges(4,6) * t196 + t438 * t141 + t341 * t454) * t397 + (t351 * t328 - Ifges(6,2) * t442 - Ifges(6,6) * t443 + Ifges(5,4) * t444 + (Ifges(5,1) + Ifges(7,3)) * t418 + Ifges(7,6) * t422 + Ifges(7,5) * t423 + t455) * t185 + (Ifges(7,5) * t84 + Ifges(7,6) * t83) * t411 + m(7) * (t1 * t18 + t15 * t35 + t17 * t2 + t3 * t6 + t4 * t5 + t48 * t7) + m(6) * (t16 * t74 + t19 * t75 + t24 * t47 + t25 * t46 + t31 * t92 + t38 * t72) + m(5) * (t154 * t189 + t164 * t172 + t20 * t79 + t21 * t80 + t26 * t51 + t27 * t52) + m(4) * (t100 * t175 + t101 * t174 + t112 * t158 + t113 * t157 + t207 * t243 + t230 * t232) + m(3) * (t229 * t251 - t230 * t250 - t238 * t243 + t241 * t242) + (-t433 / 0.2e1 - Ifges(4,6) * t403 - Ifges(4,5) * t404 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - t326) * t328 - Ifges(6,4) * t442 - Ifges(6,5) * t443 - Ifges(5,6) * t444 - Ifges(5,5) * t418 + t450) * t370 + (-Ifges(5,4) * t418 - Ifges(5,2) * t444 + t328 * t350 + t449 + t451) * t184 + (Ifges(7,1) * t84 + Ifges(7,4) * t83) * t414 + (Ifges(4,1) * t197 + Ifges(4,4) * t196) * t400 + (Ifges(4,4) * t249 + Ifges(4,2) * t248) * t403 + (Ifges(5,1) * t406 + Ifges(5,4) * t410 + Ifges(7,5) * t414 - Ifges(6,2) * t407 - Ifges(6,6) * t409 + Ifges(7,6) * t416 + Ifges(7,3) * t411 + t439 * t397 + t452) * t142 + (Ifges(4,4) * t197 + Ifges(4,2) * t196) * t401 + (Ifges(7,4) * t84 + Ifges(7,2) * t83) * t416 + (Ifges(4,1) * t249 + Ifges(4,4) * t248) * t404 + (-Ifges(5,4) * t406 - Ifges(5,2) * t410 + Ifges(6,6) * t407 + Ifges(6,3) * t409 + t431) * t141 + (Ifges(6,4) * t407 + Ifges(4,5) * t400 + Ifges(5,5) * t406 + Ifges(6,5) * t409 + Ifges(4,6) * t401 + Ifges(5,6) * t410 + t428) * t341 + ((Ifges(3,5) * t286 / 0.2e1 - t250 * mrSges(3,3) + (-0.2e1 * t419 + 0.3e1 / 0.2e1 * Ifges(3,4) * t292) * t284) * t292 + (-Ifges(3,6) * t286 + Ifges(4,5) * t398 + Ifges(4,6) * t399 - t251 * mrSges(3,3) + (-0.2e1 * t420 - 0.3e1 / 0.2e1 * t386) * t284) * t289) * t335 + (t229 * t292 + t230 * t289 + (-t238 * t292 - t241 * t289) * qJD(2)) * t284 * mrSges(3,3) + (t292 * t203 + t272 * (-Ifges(3,6) * t289 + t380) + (t151 + t97 + t94) * t289) * t284 * t388 + t116 * t398 + t115 * t399 + (-t374 - t230 * mrSges(3,1) + t260 / 0.2e1) * t286 + t363 * t243 + t230 * (-mrSges(4,1) * t248 + mrSges(4,2) * t249) + t242 * t237 + t232 * t131 + t207 * (-mrSges(4,1) * t196 + mrSges(4,2) * t197) + t197 * t153 / 0.2e1 + t113 * t194 + t196 * t152 / 0.2e1 + t112 * t193 + t189 * t68 + t174 * t170 + t175 * t171 + t172 * t102 + t161 * t10 / 0.2e1 + t27 * t145 + t26 * t146 + t24 * t143 + t25 * t144 + (-Ifges(7,1) * t300 + Ifges(7,4) * t161) * t423 + (-Ifges(7,5) * t300 + Ifges(7,6) * t161) * t418 - t300 * t427 + (-Ifges(7,4) * t300 + Ifges(7,2) * t161) * t422 + (t1 * t161 + t2 * t300 - t5 * t84 + t6 * t83) * mrSges(7,3) + t7 * (-mrSges(7,1) * t161 - mrSges(7,2) * t300) + t17 * t33 + t18 * t34 + t48 * t23 + t15 * t82 + t83 * t44 / 0.2e1 + t84 * t45 / 0.2e1 + t35 * (-mrSges(7,1) * t83 + mrSges(7,2) * t84) + t3 * t90 + t4 * t91 + t92 * t69 + t38 * t103 + t74 * t107 + t75 * t108 + t80 * t109 + t79 * t110; (t64 / 0.2e1 + t65 / 0.2e1 + t349 * t129 + (Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(7,3) / 0.2e1) * t130 + t455) * t253 + (Ifges(5,4) * t201 - Ifges(5,2) * t200 - Ifges(6,6) * t247 + Ifges(6,3) * t246) * t409 + t260 + (Ifges(5,1) * t247 - Ifges(5,4) * t246 - Ifges(6,2) * t201 + Ifges(6,6) * t200) * t406 + (t108 - t110) * t205 + (t154 * t281 + t164 * t468 - t20 * t205 + t206 * t21 + t457 * t52 + t460 * t51) * m(5) + t457 * t145 + t458 * t143 + (-t16 * t206 + t19 * t205 + t192 * t31 + t456 * t72 + t458 * t47 + t459 * t46) * m(6) + t459 * t144 + t460 * t146 + t461 * t82 + (Ifges(7,5) * t183 + Ifges(7,6) * t182 + Ifges(7,3) * t201) * t412 + (Ifges(7,1) * t367 + Ifges(7,4) * t366 + Ifges(7,5) * t247) * t414 + (Ifges(7,1) * t183 + Ifges(7,4) * t182 + Ifges(7,5) * t201) * t415 + (Ifges(7,4) * t367 + Ifges(7,2) * t366 + Ifges(7,6) * t247) * t416 + (Ifges(7,4) * t183 + Ifges(7,2) * t182 + Ifges(7,6) * t201) * t417 + (Ifges(5,1) * t201 - Ifges(5,4) * t200 - Ifges(6,2) * t247 + Ifges(6,6) * t246) * t407 + (t109 - t107) * t206 + (t288 * pkin(3) * t102 + (-t193 * t288 - t194 * t291) * pkin(9) + t432) * qJD(3) + (t316 * t423 + t314 * t422 - t7 * t319 + t10 * t394 + t287 * t427 + (t379 / 0.2e1 + t378 / 0.2e1 + t349) * t130 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t129 + t324 * mrSges(7,3) + (-mrSges(7,3) * t322 + t313 * t411 + t315 * t416 + t317 * t414 + t318 * t35 + t394 * t45 + t395 * t44) * qJD(6) + t449) * t252 + ((t238 * mrSges(3,3) + t360 * t419 - t373 / 0.2e1 - t203 / 0.2e1 - t265 / 0.2e1 - t432) * t292 + (-t428 - t377 / 0.2e1 - t375 / 0.2e1 - t97 / 0.2e1 - t94 / 0.2e1 + (t386 / 0.2e1 + t420 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t292) * t360 - t151 / 0.2e1 - t350 * t166 + t241 * mrSges(3,3) - t351 * t304 + t326 * t262 + (t272 / 0.2e1 - qJD(2)) * Ifges(3,6) + (Ifges(4,5) * t288 + Ifges(4,6) * t291 + t252 * t438 + t253 * t439) * t388) * t289) * t360 + t115 * t393 + (-t170 * t288 + t171 * t291) * pkin(9) + (-mrSges(4,1) * t291 + mrSges(4,2) * t288 - mrSges(3,1)) * t230 - t374 + (-pkin(2) * t230 - t157 * t176 - t158 * t177 - t207 * t241 + (-qJD(3) * t305 + t306) * pkin(9)) * m(4) + (Ifges(4,2) * t291 + t385) * t403 + (Ifges(4,1) * t288 + t384) * t404 + (Ifges(7,5) * t367 + Ifges(7,6) * t366 + Ifges(7,3) * t247) * t411 + (t367 / 0.2e1 - t183 / 0.2e1) * t45 + (mrSges(6,2) * t362 + mrSges(6,3) * t361) * t72 + (t361 * t51 + t362 * t52) * mrSges(5,3) + (-t361 * t46 - t362 * t47) * mrSges(6,1) + (-mrSges(5,1) * t362 - mrSges(5,2) * t361) * t164 - t363 * t241 + (t366 / 0.2e1 - t182 / 0.2e1) * t44 + (mrSges(7,2) * t361 + mrSges(7,3) * t332) * t6 + (-mrSges(7,1) * t361 - mrSges(7,3) * t331) * t5 + t288 * t116 / 0.2e1 + t281 * t68 + (t98 + t43 - t95) * (t247 / 0.2e1 - t201 / 0.2e1) + (t93 - t96) * (t246 / 0.2e1 - t200 / 0.2e1) + (-mrSges(7,1) * t332 + mrSges(7,2) * t331) * t35 - t238 * t237 - t198 * t102 - t176 * t194 + t192 * t69 - t177 * t193 + t181 * t23 - pkin(2) * t131 + t462 * t90 + (t1 * t105 + t104 * t2 + t181 * t7 + t35 * t461 + t462 * t6 + t463 * t5) * m(7) + t463 * t91 + t456 * t103 + (t246 * t438 + t247 * t439) * t397 + (t200 * t438 + t201 * t439) * t396 + (Ifges(5,4) * t247 - Ifges(5,2) * t246 - Ifges(6,6) * t201 + Ifges(6,3) * t200) * t410 + t306 * mrSges(4,3) + t104 * t33 + t105 * t34; (t23 - t107) * t278 + (-t102 * t221 + t109 * t283 + t110 * t285) * pkin(3) + ((-m(7) * t321 + t309) * qJD(6) + m(7) * t434 + t310) * (-pkin(10) + t280) + t313 * t418 + t315 * t422 + t317 * t423 + (t157 * t220 + t158 * t221) * mrSges(4,3) + (-Ifges(5,1) * t407 - Ifges(5,4) * t409 - Ifges(7,5) * t415 + Ifges(6,2) * t406 + Ifges(6,6) * t410 - Ifges(7,6) * t417 - Ifges(7,3) * t412 - t439 * t396 + t452) * t166 + (-t16 * t278 + t19 * t280 - t46 * t61 - t72 * t87 + (-qJD(5) + t62) * t47) * m(6) + (-Ifges(5,4) * t407 - Ifges(5,2) * t409 + Ifges(6,6) * t406 + Ifges(6,3) * t410 + t396 * t438 + t430 - t431) * t304 + ((t20 * t285 + t21 * t283) * pkin(3) - t164 * t391 + t51 * t61 - t52 * t62) * m(5) + (-t12 * t5 - t13 * t6 + t278 * t7 + (qJD(5) - t40) * t35) * m(7) + t11 * t394 + t10 * t395 + (Ifges(4,5) * t220 - Ifges(4,6) * t221) * t396 + t430 * qJD(6) + t152 * t400 - t221 * (Ifges(4,1) * t220 - t376) / 0.2e1 + t372 * qJD(5) + t433 + t364 * t61 + t365 * t62 + t280 * t108 - t434 * mrSges(7,3) - t207 * (mrSges(4,1) * t221 + mrSges(4,2) * t220) + t158 * t194 - t157 * t193 - t450 - (-Ifges(4,2) * t221 + t153 + t216) * t220 / 0.2e1 + t7 * t318 - t40 * t82 - t13 * t90 - t12 * t91 - t87 * t103; -t287 * t33 + t290 * t34 + t120 - t121 + (mrSges(5,1) - mrSges(6,2)) * t129 + t308 * qJD(6) - (-t82 + t365) * t166 + (t308 + t364) * t304 + (-t163 * t322 + t166 * t35 + t324) * m(7) + (-t166 * t47 - t304 * t46 + t31) * m(6) + (t166 * t52 + t304 * t51 + t154) * m(5); t372 * t262 + t309 * qJD(6) + (t103 + t309) * t304 + t310 + t108 + (-t163 * t321 + t262 * t35 + t434) * m(7) + (-t262 * t47 + t304 * t72 + t19) * m(6); -t35 * (mrSges(7,1) * t137 + mrSges(7,2) * t136) + (Ifges(7,1) * t136 - t383) * t415 + t44 * t414 + (Ifges(7,5) * t136 - Ifges(7,6) * t137) * t412 - t5 * t90 + t6 * t91 + (t136 * t5 + t137 * t6) * mrSges(7,3) + t325 + t9 + (-Ifges(7,2) * t137 + t135 + t45) * t417;];
tauc  = t28(:);
