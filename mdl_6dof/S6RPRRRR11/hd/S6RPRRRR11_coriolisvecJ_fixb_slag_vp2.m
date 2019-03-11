% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:19
% EndTime: 2019-03-09 07:39:19
% DurationCPUTime: 29.83s
% Computational Cost: add. (40516->867), mult. (132390->1231), div. (0->0), fcn. (113823->14), ass. (0->386)
t320 = sin(pkin(13));
t322 = sin(pkin(6));
t323 = cos(pkin(13));
t333 = cos(qJ(3));
t324 = cos(pkin(7));
t329 = sin(qJ(3));
t404 = t324 * t329;
t341 = (-t320 * t404 + t323 * t333) * t322;
t277 = qJD(1) * t341;
t321 = sin(pkin(7));
t386 = qJD(3) * t333;
t540 = -t321 * t386 + t277;
t325 = cos(pkin(6));
t403 = t324 * t333;
t407 = t321 * t333;
t337 = t325 * t407 + t322 * (-t320 * t329 + t323 * t403);
t254 = t337 * qJD(1);
t408 = t321 * t329;
t262 = t325 * t408 + (t320 * t333 + t323 * t404) * t322;
t257 = t262 * qJD(1);
t331 = cos(qJ(5));
t327 = sin(qJ(5));
t332 = cos(qJ(4));
t401 = t327 * t332;
t206 = -t254 * t401 + t257 * t331;
t328 = sin(qJ(4));
t382 = qJD(5) * t331;
t384 = qJD(4) * t332;
t539 = t327 * t384 + t328 * t382 + t206;
t366 = pkin(4) * t328 - pkin(11) * t332;
t304 = t366 * qJD(4);
t306 = -pkin(4) * t332 - pkin(11) * t328 - pkin(3);
t383 = qJD(5) * t327;
t385 = qJD(4) * t328;
t236 = t327 * t304 + t306 * t382 + (-t331 * t385 - t332 * t383) * pkin(10);
t406 = t322 * t323;
t314 = qJ(2) * t406;
t446 = pkin(1) * t325;
t381 = qJD(1) * t446;
t285 = qJD(1) * t314 + t320 * t381;
t344 = (t321 * t325 + t324 * t406) * pkin(9);
t246 = qJD(1) * t344 + t285;
t312 = t323 * t381;
t410 = t320 * t322;
t340 = pkin(2) * t325 + (-pkin(9) * t324 - qJ(2)) * t410;
t253 = qJD(1) * t340 + t312;
t278 = (-pkin(9) * t320 * t321 - pkin(2) * t323 - pkin(1)) * t322;
t269 = qJD(1) * t278 + qJD(2);
t355 = t253 * t324 + t269 * t321;
t179 = -t329 * t246 + t355 * t333;
t210 = pkin(3) * t257 - pkin(10) * t254;
t132 = t332 * t179 + t328 * t210;
t112 = pkin(11) * t257 + t132;
t180 = t246 * t333 + t329 * t355;
t133 = t254 * t366 + t180;
t70 = t331 * t112 + t327 * t133;
t538 = -t70 + t236;
t291 = -t332 * t324 + t328 * t408;
t389 = qJD(1) * t322;
t375 = t320 * t389;
t369 = t321 * t375;
t537 = qJD(4) * t291 + t328 * t369 + t332 * t540;
t342 = (t320 * t403 + t323 * t329) * t322;
t276 = qJD(1) * t342;
t387 = qJD(3) * t329;
t536 = -t321 * t387 + t276;
t399 = t331 * t332;
t207 = t254 * t399 + t257 * t327;
t317 = pkin(10) * t399;
t444 = pkin(10) * t327;
t391 = t331 * t304 + t385 * t444;
t412 = t254 * t328;
t69 = -t112 * t327 + t331 * t133;
t535 = -pkin(5) * t412 + pkin(12) * t207 + (pkin(5) * t328 - pkin(12) * t399) * qJD(4) + (-t317 + (pkin(12) * t328 - t306) * t327) * qJD(5) + t391 - t69;
t534 = -pkin(12) * t539 + t538;
t288 = -t321 * t406 + t324 * t325;
t280 = qJD(1) * t288 + qJD(3);
t224 = -t328 * t257 + t280 * t332;
t221 = qJD(5) - t224;
t225 = t257 * t332 + t280 * t328;
t252 = qJD(4) - t254;
t190 = t225 * t331 + t252 * t327;
t156 = -t280 * pkin(3) - t179;
t100 = -t224 * pkin(4) - t225 * pkin(11) + t156;
t215 = -t253 * t321 + t324 * t269;
t154 = -pkin(3) * t254 - pkin(10) * t257 + t215;
t157 = t280 * pkin(10) + t180;
t89 = t154 * t328 + t157 * t332;
t82 = pkin(11) * t252 + t89;
t46 = t331 * t100 - t327 * t82;
t37 = -pkin(12) * t190 + t46;
t31 = pkin(5) * t221 + t37;
t330 = cos(qJ(6));
t326 = sin(qJ(6));
t189 = -t225 * t327 + t252 * t331;
t47 = t100 * t327 + t331 * t82;
t38 = pkin(12) * t189 + t47;
t421 = t326 * t38;
t16 = t31 * t330 - t421;
t418 = t330 * t38;
t17 = t31 * t326 + t418;
t533 = -t46 * mrSges(6,1) - t16 * mrSges(7,1) + t47 * mrSges(6,2) + t17 * mrSges(7,2);
t423 = t252 * Ifges(5,6);
t441 = Ifges(5,4) * t225;
t532 = -t156 * mrSges(5,1) + t89 * mrSges(5,3) + t533 + t423 / 0.2e1 + t441 / 0.2e1;
t426 = t224 * Ifges(5,2);
t531 = -t426 / 0.2e1;
t88 = t154 * t332 - t328 * t157;
t526 = t88 * mrSges(5,3);
t525 = t156 * mrSges(5,2);
t477 = -pkin(12) - pkin(11);
t376 = qJD(5) * t477;
t413 = t224 * t327;
t169 = pkin(4) * t225 - pkin(11) * t224;
t73 = t327 * t169 + t331 * t88;
t524 = pkin(12) * t413 + t327 * t376 - t73;
t72 = t331 * t169 - t327 * t88;
t523 = -pkin(5) * t225 - t72 + (pkin(12) * t224 + t376) * t331;
t292 = t324 * t328 + t332 * t408;
t349 = -t331 * t292 + t327 * t407;
t492 = qJD(5) * t349 + t327 * t537 - t331 * t536;
t270 = -t327 * t292 - t331 * t407;
t522 = -qJD(5) * t270 + t327 * t536 + t331 * t537;
t357 = t327 * t47 + t331 * t46;
t359 = Ifges(6,5) * t331 - Ifges(6,6) * t327;
t436 = Ifges(6,4) * t331;
t361 = -Ifges(6,2) * t327 + t436;
t437 = Ifges(6,4) * t327;
t363 = Ifges(6,1) * t331 - t437;
t364 = mrSges(6,1) * t327 + mrSges(6,2) * t331;
t448 = t331 / 0.2e1;
t449 = -t327 / 0.2e1;
t455 = t221 / 0.2e1;
t460 = t190 / 0.2e1;
t462 = t189 / 0.2e1;
t81 = -pkin(4) * t252 - t88;
t438 = Ifges(6,4) * t190;
t96 = Ifges(6,2) * t189 + Ifges(6,6) * t221 + t438;
t186 = Ifges(6,4) * t189;
t97 = Ifges(6,1) * t190 + Ifges(6,5) * t221 + t186;
t521 = t357 * mrSges(6,3) - t359 * t455 - t361 * t462 - t363 * t460 - t364 * t81 - t448 * t97 - t449 * t96;
t255 = t337 * qJD(3);
t244 = qJD(1) * t255;
t187 = qJD(4) * t224 + t244 * t332;
t256 = t262 * qJD(3);
t245 = qJD(1) * t256;
t103 = qJD(5) * t189 + t187 * t331 + t245 * t327;
t338 = qJD(2) * t341;
t149 = qJD(1) * t338 + qJD(3) * t179;
t388 = qJD(2) * t322;
t374 = t320 * t388;
t367 = qJD(1) * t374;
t351 = t321 * t367;
t200 = pkin(3) * t245 - pkin(10) * t244 + t351;
t53 = t332 * t149 + t154 * t384 - t157 * t385 + t328 * t200;
t49 = pkin(11) * t245 + t53;
t339 = qJD(2) * t342;
t150 = qJD(1) * t339 + qJD(3) * t180;
t188 = qJD(4) * t225 + t244 * t328;
t79 = t188 * pkin(4) - t187 * pkin(11) + t150;
t14 = -qJD(5) * t47 - t327 * t49 + t331 * t79;
t6 = pkin(5) * t188 - pkin(12) * t103 + t14;
t104 = -qJD(5) * t190 - t187 * t327 + t245 * t331;
t13 = t100 * t382 + t327 * t79 + t331 * t49 - t383 * t82;
t7 = pkin(12) * t104 + t13;
t2 = qJD(6) * t16 + t326 * t6 + t330 * t7;
t3 = -qJD(6) * t17 - t326 * t7 + t330 * t6;
t520 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t519 = t14 * mrSges(6,1) - t13 * mrSges(6,2);
t370 = t330 * t189 - t190 * t326;
t118 = Ifges(7,4) * t370;
t122 = t189 * t326 + t190 * t330;
t435 = Ifges(7,4) * t122;
t217 = qJD(6) + t221;
t458 = -t217 / 0.2e1;
t470 = -t122 / 0.2e1;
t472 = -t370 / 0.2e1;
t64 = Ifges(7,1) * t122 + Ifges(7,5) * t217 + t118;
t71 = -pkin(5) * t189 + t81;
t34 = qJD(6) * t370 + t103 * t330 + t104 * t326;
t35 = -qJD(6) * t122 - t103 * t326 + t104 * t330;
t9 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t188;
t518 = t9 + (Ifges(7,5) * t370 - Ifges(7,6) * t122) * t458 + (t122 * t17 + t16 * t370) * mrSges(7,3) + (-Ifges(7,2) * t122 + t118 + t64) * t472 - t71 * (mrSges(7,1) * t122 + mrSges(7,2) * t370) + t520 + (Ifges(7,1) * t370 - t435) * t470;
t453 = t245 / 0.2e1;
t465 = -t188 / 0.2e1;
t466 = t187 / 0.2e1;
t473 = Ifges(5,1) * t466 + Ifges(5,4) * t465 + Ifges(5,5) * t453;
t42 = Ifges(6,5) * t103 + Ifges(6,6) * t104 + Ifges(6,3) * t188;
t516 = t42 + t9;
t515 = Ifges(6,3) + Ifges(7,3);
t428 = t217 * Ifges(7,3);
t433 = t122 * Ifges(7,5);
t434 = t370 * Ifges(7,6);
t62 = t428 + t433 + t434;
t427 = t221 * Ifges(6,3);
t429 = t190 * Ifges(6,5);
t430 = t189 * Ifges(6,6);
t95 = t427 + t429 + t430;
t514 = t95 + t62;
t298 = t331 * t306;
t400 = t328 * t331;
t260 = -pkin(12) * t400 + t298 + (-pkin(5) - t444) * t332;
t282 = t327 * t306 + t317;
t402 = t327 * t328;
t272 = -pkin(12) * t402 + t282;
t219 = t260 * t326 + t272 * t330;
t513 = -qJD(6) * t219 - t534 * t326 + t330 * t535;
t218 = t260 * t330 - t272 * t326;
t512 = qJD(6) * t218 + t326 * t535 + t534 * t330;
t131 = -t328 * t179 + t210 * t332;
t111 = -pkin(4) * t257 - t131;
t511 = pkin(5) * t539 + pkin(10) * t384 - t111;
t390 = t320 * t446 + t314;
t258 = t344 + t390;
t316 = t323 * t446;
t263 = t316 + t340;
t354 = t263 * t324 + t278 * t321;
t193 = -t329 * t258 + t354 * t333;
t220 = Ifges(5,4) * t224;
t424 = t252 * Ifges(5,5);
t425 = t225 * Ifges(5,1);
t140 = t220 + t424 + t425;
t467 = t140 / 0.2e1;
t510 = -t220 / 0.2e1 - t525 - t467 - t424 / 0.2e1 + t526 + t521;
t481 = t34 / 0.2e1;
t480 = t35 / 0.2e1;
t63 = Ifges(7,2) * t370 + Ifges(7,6) * t217 + t435;
t508 = t63 / 0.2e1;
t475 = t103 / 0.2e1;
t474 = t104 / 0.2e1;
t464 = t188 / 0.2e1;
t507 = t88 * mrSges(5,1);
t506 = t89 * mrSges(5,2);
t503 = Ifges(4,5) * t244;
t501 = Ifges(4,6) * t245;
t500 = t149 * mrSges(4,2);
t309 = t477 * t327;
t310 = t477 * t331;
t275 = t309 * t326 - t310 * t330;
t498 = -qJD(6) * t275 - t326 * t524 + t330 * t523;
t274 = t309 * t330 + t310 * t326;
t497 = qJD(6) * t274 + t326 * t523 + t330 * t524;
t68 = -mrSges(7,1) * t370 + mrSges(7,2) * t122;
t496 = m(7) * t71 + t68;
t177 = -t288 * pkin(3) - t193;
t234 = t262 * t328 - t288 * t332;
t235 = t262 * t332 + t288 * t328;
t115 = t234 * pkin(4) - t235 * pkin(11) + t177;
t227 = -t263 * t321 + t324 * t278;
t171 = -pkin(3) * t337 - pkin(10) * t262 + t227;
t250 = t333 * t258;
t194 = t263 * t404 + t278 * t408 + t250;
t178 = pkin(10) * t288 + t194;
t108 = t328 * t171 + t332 * t178;
t94 = -pkin(11) * t337 + t108;
t57 = t327 * t115 + t331 * t94;
t352 = t326 * t327 - t330 * t331;
t287 = t352 * t328;
t222 = t270 * t330 + t326 * t349;
t494 = qJD(6) * t222 + t326 * t492 - t330 * t522;
t223 = t270 * t326 - t330 * t349;
t493 = -qJD(6) * t223 + t326 * t522 + t330 * t492;
t442 = mrSges(4,3) * t257;
t393 = -mrSges(4,1) * t280 - mrSges(5,1) * t224 + mrSges(5,2) * t225 + t442;
t489 = qJD(4) * t292 - t328 * t540 + t332 * t369;
t54 = -t328 * t149 - t154 * t385 - t157 * t384 + t200 * t332;
t488 = -t328 * t54 + t332 * t53;
t487 = t13 * t331 - t14 * t327;
t486 = qJD(5) + qJD(6);
t484 = t54 * mrSges(5,1) - t53 * mrSges(5,2);
t483 = Ifges(7,4) * t481 + Ifges(7,2) * t480 + Ifges(7,6) * t464;
t482 = Ifges(7,1) * t481 + Ifges(7,4) * t480 + Ifges(7,5) * t464;
t479 = Ifges(6,1) * t475 + Ifges(6,4) * t474 + Ifges(6,5) * t464;
t478 = -t96 / 0.2e1;
t471 = t370 / 0.2e1;
t469 = t122 / 0.2e1;
t139 = t423 + t426 + t441;
t468 = t139 / 0.2e1;
t463 = -t189 / 0.2e1;
t461 = -t190 / 0.2e1;
t457 = t217 / 0.2e1;
t456 = -t221 / 0.2e1;
t450 = t257 / 0.2e1;
t445 = pkin(5) * t327;
t443 = mrSges(4,3) * t254;
t440 = Ifges(5,4) * t328;
t439 = Ifges(5,4) * t332;
t422 = t257 * Ifges(4,4);
t50 = -pkin(4) * t245 - t54;
t420 = t328 * t50;
t147 = mrSges(5,1) * t245 - mrSges(5,3) * t187;
t55 = -mrSges(6,1) * t104 + mrSges(6,2) * t103;
t416 = -t147 + t55;
t415 = t150 * t333;
t405 = t323 * (-mrSges(3,2) * t325 + mrSges(3,3) * t406) * qJD(1);
t126 = -mrSges(6,1) * t189 + mrSges(6,2) * t190;
t192 = mrSges(5,1) * t252 - mrSges(5,3) * t225;
t398 = t126 - t192;
t141 = t206 * t330 - t207 * t326;
t300 = t326 * t331 + t327 * t330;
t231 = t287 * t486 - t300 * t384;
t397 = t141 - t231;
t142 = t206 * t326 + t207 * t330;
t265 = t486 * t300;
t230 = -t265 * t328 - t352 * t384;
t396 = t142 - t230;
t158 = t300 * t224;
t395 = -t158 + t265;
t159 = t352 * t224;
t264 = t486 * t352;
t394 = -t159 + t264;
t392 = -t501 + t503;
t377 = Ifges(5,5) * t187 - Ifges(5,6) * t188 + Ifges(5,3) * t245;
t56 = t331 * t115 - t327 * t94;
t107 = t171 * t332 - t328 * t178;
t368 = t321 * t374;
t365 = mrSges(6,1) * t331 - mrSges(6,2) * t327;
t362 = Ifges(6,1) * t327 + t436;
t360 = Ifges(6,2) * t331 + t437;
t358 = Ifges(6,5) * t327 + Ifges(6,6) * t331;
t202 = t235 * t331 - t327 * t337;
t39 = pkin(5) * t234 - pkin(12) * t202 + t56;
t201 = -t235 * t327 - t331 * t337;
t45 = pkin(12) * t201 + t57;
t22 = -t326 * t45 + t330 * t39;
t23 = t326 * t39 + t330 * t45;
t356 = t328 * t89 + t332 * t88;
t134 = t201 * t330 - t202 * t326;
t135 = t201 * t326 + t202 * t330;
t353 = -(-qJ(2) * t375 + t312) * t320 + t285 * t323;
t164 = qJD(3) * t193 + t338;
t204 = pkin(3) * t256 - pkin(10) * t255 + t368;
t67 = -t328 * t164 - t171 * t385 - t178 * t384 + t204 * t332;
t93 = pkin(4) * t337 - t107;
t66 = t332 * t164 + t171 * t384 - t178 * t385 + t328 * t204;
t60 = pkin(11) * t256 + t66;
t165 = t339 + (t329 * t354 + t250) * qJD(3);
t198 = qJD(4) * t235 + t255 * t328;
t199 = -qJD(4) * t234 + t255 * t332;
t85 = t198 * pkin(4) - t199 * pkin(11) + t165;
t20 = t115 * t382 + t327 * t85 + t331 * t60 - t383 * t94;
t289 = (mrSges(3,1) * t325 - mrSges(3,3) * t410) * qJD(1);
t61 = -pkin(4) * t256 - t67;
t21 = -qJD(5) * t57 - t327 * t60 + t331 * t85;
t335 = -t430 / 0.2e1 - t429 / 0.2e1 - t427 / 0.2e1 - t62 / 0.2e1 - t95 / 0.2e1 - t434 / 0.2e1 - t433 / 0.2e1 - t428 / 0.2e1 + t468 + t532;
t319 = -pkin(5) * t331 - pkin(4);
t305 = (pkin(10) + t445) * t328;
t286 = t300 * t328;
t281 = -pkin(10) * t401 + t298;
t251 = Ifges(4,4) * t254;
t237 = -qJD(5) * t282 + t391;
t228 = -mrSges(4,2) * t280 + t443;
t209 = -mrSges(4,1) * t254 + mrSges(4,2) * t257;
t205 = mrSges(4,1) * t245 + mrSges(4,2) * t244;
t197 = t257 * Ifges(4,1) + t280 * Ifges(4,5) + t251;
t196 = t254 * Ifges(4,2) + t280 * Ifges(4,6) + t422;
t191 = -mrSges(5,2) * t252 + mrSges(5,3) * t224;
t148 = -mrSges(5,2) * t245 - mrSges(5,3) * t188;
t138 = t225 * Ifges(5,5) + t224 * Ifges(5,6) + t252 * Ifges(5,3);
t137 = mrSges(6,1) * t221 - mrSges(6,3) * t190;
t136 = -mrSges(6,2) * t221 + mrSges(6,3) * t189;
t123 = mrSges(5,1) * t188 + mrSges(5,2) * t187;
t117 = qJD(5) * t201 + t199 * t331 + t256 * t327;
t116 = -qJD(5) * t202 - t199 * t327 + t256 * t331;
t105 = t187 * Ifges(5,4) - t188 * Ifges(5,2) + t245 * Ifges(5,6);
t92 = mrSges(7,1) * t217 - mrSges(7,3) * t122;
t91 = -mrSges(7,2) * t217 + mrSges(7,3) * t370;
t80 = pkin(5) * t413 + t89;
t76 = -mrSges(6,2) * t188 + mrSges(6,3) * t104;
t75 = mrSges(6,1) * t188 - mrSges(6,3) * t103;
t74 = -pkin(5) * t201 + t93;
t43 = Ifges(6,4) * t103 + Ifges(6,2) * t104 + Ifges(6,6) * t188;
t41 = -qJD(6) * t135 + t116 * t330 - t117 * t326;
t40 = qJD(6) * t134 + t116 * t326 + t117 * t330;
t36 = -pkin(5) * t116 + t61;
t30 = -pkin(5) * t104 + t50;
t29 = -mrSges(7,2) * t188 + mrSges(7,3) * t35;
t28 = mrSges(7,1) * t188 - mrSges(7,3) * t34;
t19 = t330 * t37 - t421;
t18 = -t326 * t37 - t418;
t15 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t12 = pkin(12) * t116 + t20;
t8 = pkin(5) * t198 - pkin(12) * t117 + t21;
t5 = -qJD(6) * t23 - t12 * t326 + t330 * t8;
t4 = qJD(6) * t22 + t12 * t330 + t326 * t8;
t1 = [(Ifges(6,5) * t117 + Ifges(6,6) * t116) * t455 + t225 * (Ifges(5,1) * t199 + Ifges(5,5) * t256) / 0.2e1 + (t116 * t47 - t117 * t46 + t13 * t201 - t14 * t202) * mrSges(6,3) - (t377 / 0.2e1 + mrSges(4,1) * t351 + Ifges(5,3) * t453 + Ifges(5,6) * t465 + Ifges(5,5) * t466 - Ifges(4,4) * t244 + Ifges(4,2) * t245 - t149 * mrSges(4,3) + t484) * t337 + t224 * (Ifges(5,4) * t199 + Ifges(5,6) * t256) / 0.2e1 + (t514 / 0.2e1 - t139 / 0.2e1 + t531 + Ifges(6,3) * t455 + Ifges(7,3) * t457 + Ifges(6,5) * t460 + Ifges(6,6) * t462 + Ifges(7,5) * t469 + Ifges(7,6) * t471 - t532) * t198 + (Ifges(6,1) * t117 + Ifges(6,4) * t116) * t460 + (t134 * t2 - t135 * t3 - t16 * t40 + t17 * t41) * mrSges(7,3) + (Ifges(6,1) * t202 + Ifges(6,4) * t201) * t475 + (Ifges(7,5) * t40 + Ifges(7,6) * t41) * t457 + t209 * t368 - 0.2e1 * t289 * t374 + (mrSges(4,2) * t351 + mrSges(4,3) * t150 + Ifges(4,1) * t244 - Ifges(4,4) * t245) * t262 + m(7) * (t16 * t5 + t17 * t4 + t2 * t23 + t22 * t3 + t30 * t74 + t36 * t71) + m(6) * (t13 * t57 + t14 * t56 + t20 * t47 + t21 * t46 + t50 * t93 + t61 * t81) + t252 * (Ifges(5,5) * t199 + Ifges(5,3) * t256) / 0.2e1 + (Ifges(6,5) * t202 + Ifges(7,5) * t135 + Ifges(6,6) * t201 + Ifges(7,6) * t134) * t464 + (m(3) * ((t323 * t390 + (qJ(2) * t410 - t316) * t320) * qJD(1) + t353) + 0.2e1 * t405) * t388 + m(5) * (t107 * t54 + t108 * t53 + t150 * t177 + t66 * t89 + t67 * t88) + (-t179 * t255 - t180 * t256 - t193 * t244 - t194 * t245) * mrSges(4,3) + m(4) * (t149 * t194 - t150 * t193 + t164 * t180 + (qJD(1) * t227 + t215) * t368) + (Ifges(7,4) * t40 + Ifges(7,2) * t41) * t471 + (Ifges(7,4) * t135 + Ifges(7,2) * t134) * t480 - t199 * t526 + t71 * (-mrSges(7,1) * t41 + mrSges(7,2) * t40) + t36 * t68 + (Ifges(6,4) * t117 + Ifges(6,2) * t116) * t462 + (Ifges(6,4) * t202 + Ifges(6,2) * t201) * t474 + t40 * t64 / 0.2e1 + (-Ifges(5,4) * t466 - t53 * mrSges(5,3) - Ifges(5,6) * t453 - Ifges(5,2) * t465 + t515 * t464 + t516 / 0.2e1 - t105 / 0.2e1 + mrSges(5,1) * t150 + Ifges(6,6) * t474 + Ifges(6,5) * t475 + Ifges(7,6) * t480 + Ifges(7,5) * t481 + t519 + t520) * t234 + (mrSges(5,2) * t150 - mrSges(5,3) * t54 + 0.2e1 * t473) * t235 + t23 * t29 + t22 * t28 + t74 * t15 + t56 * t75 + t57 * t76 + (-m(4) * t179 + m(5) * t156 + t393) * t165 + t4 * t91 + t5 * t92 + t93 * t55 + t199 * t525 + t256 * t507 + t41 * t508 + (Ifges(7,1) * t40 + Ifges(7,4) * t41) * t469 + (Ifges(7,1) * t135 + Ifges(7,4) * t134) * t481 + t116 * t96 / 0.2e1 + t117 * t97 / 0.2e1 + t81 * (-mrSges(6,1) * t116 + mrSges(6,2) * t117) + t61 * t126 + (-t500 - t150 * mrSges(4,1) + t503 / 0.2e1 - t501 / 0.2e1 + t392 / 0.2e1) * t288 + t30 * (-mrSges(7,1) * t134 + mrSges(7,2) * t135) + t20 * t136 + t21 * t137 + t107 * t147 + t108 * t148 - t256 * t506 + t177 * t123 + t66 * t191 + t67 * t192 + t201 * t43 / 0.2e1 + t50 * (-mrSges(6,1) * t201 + mrSges(6,2) * t202) + t227 * t205 + t164 * t228 + t255 * t197 / 0.2e1 + t256 * t138 / 0.2e1 + (Ifges(4,1) * t255 - Ifges(4,4) * t256) * t450 + t199 * t467 + t202 * t479 + t135 * t482 + t134 * t483 - t256 * t196 / 0.2e1 + t215 * (mrSges(4,1) * t256 + mrSges(4,2) * t255) + t254 * (Ifges(4,4) * t255 - Ifges(4,2) * t256) / 0.2e1 + t280 * (Ifges(4,5) * t255 - Ifges(4,6) * t256) / 0.2e1; t493 * t92 - m(4) * (-t179 * t276 + t180 * t277) + t494 * t91 - t393 * t276 + t492 * t137 - t522 * t136 - t537 * t191 + (t15 + t416) * t291 + t222 * t28 + t223 * t29 + t270 * t75 - t349 * t76 - t277 * t228 + t292 * t148 + t324 * t205 + (-m(3) * t353 + t289 * t320 - t405) * t389 + t489 * (t68 + t398) + (t16 * t493 + t17 * t494 + t2 * t223 + t222 * t3 + t291 * t30 + t489 * t71) * m(7) + (-t13 * t349 + t14 * t270 + t291 * t50 + t46 * t492 - t47 * t522 + t489 * t81) * m(6) + (-t156 * t276 - t54 * t291 + t53 * t292 - t489 * t88 - t537 * t89) * m(5) + (-t209 * t375 - t333 * t123 + (-t244 * t333 - t245 * t329) * mrSges(4,3) + (t228 * t333 + t329 * t393) * qJD(3) + m(5) * (t156 * t387 - t415) + (t149 * t329 - t179 * t387 + t180 * t386 - t215 * t375 + t324 * t367 - t415) * m(4)) * t321; (t2 * t332 - t287 * t30 - t396 * t71) * mrSges(7,2) + (t46 * mrSges(6,3) - t97 / 0.2e1 - t81 * mrSges(6,2) + Ifges(6,5) * t456 + Ifges(6,1) * t461 + Ifges(6,4) * t463) * t207 + (-t514 / 0.2e1 + Ifges(6,3) * t456 + Ifges(7,3) * t458 + Ifges(6,5) * t461 + Ifges(6,6) * t463 + t468 + Ifges(7,5) * t470 + Ifges(7,6) * t472 + t533) * t412 + (t81 * mrSges(6,1) - t47 * mrSges(6,3) + Ifges(6,4) * t461 + Ifges(6,2) * t463 + Ifges(6,6) * t456 + t478) * t206 + (-t142 / 0.2e1 + t230 / 0.2e1) * t64 + (Ifges(7,4) * t142 + Ifges(7,2) * t141) * t472 - t500 + m(6) * (pkin(10) * t420 + t13 * t282 + t14 * t281 + t236 * t47 + t237 * t46) + (Ifges(7,1) * t142 + Ifges(7,4) * t141) * t470 - t224 * (Ifges(5,6) * t257 + (-Ifges(5,2) * t328 + t439) * t254) / 0.2e1 - t225 * (Ifges(5,5) * t257 + (Ifges(5,1) * t332 - t440) * t254) / 0.2e1 + (Ifges(7,5) * t142 + Ifges(7,6) * t141) * t458 + (-t141 / 0.2e1 + t231 / 0.2e1) * t63 - (Ifges(4,1) * t254 + t138 - t422) * t257 / 0.2e1 - (-Ifges(4,2) * t257 + t197 + t251) * t254 / 0.2e1 + (-t69 + t237) * t137 - t43 * t402 / 0.2e1 + t13 * (mrSges(6,2) * t332 - mrSges(6,3) * t402) - m(6) * (t111 * t81 + t46 * t69 + t47 * t70) - m(5) * (t131 * t88 + t132 * t89 + t156 * t180) + (t16 * t396 - t17 * t397 - t2 * t286 + t287 * t3) * mrSges(7,3) + (-Ifges(7,4) * t287 - Ifges(7,2) * t286 - Ifges(7,6) * t332) * t480 + (-Ifges(7,1) * t287 - Ifges(7,4) * t286 - Ifges(7,5) * t332) * t481 + t538 * t136 + (t81 * t365 + t360 * t463 + t362 * t461 + t358 * t456 + t97 * t449 + t331 * t478 + (t327 * t46 - t331 * t47) * mrSges(6,3)) * t328 * qJD(5) + t14 * (-mrSges(6,1) * t332 - mrSges(6,3) * t400) + (-t393 + t442) * t180 + (-t228 + t443) * t179 - t252 * (Ifges(5,3) * t257 + (Ifges(5,5) * t332 - Ifges(5,6) * t328) * t254) / 0.2e1 + (-mrSges(5,1) * t332 + mrSges(5,2) * t328 - mrSges(4,1)) * t150 + t392 + ((-t335 + t531) * t328 + (-m(5) * t356 - t328 * t191) * pkin(10) + (t425 / 0.2e1 + (m(6) * t81 + t398) * pkin(10) - t510) * t332) * qJD(4) + (t148 * t332 + t328 * t416) * pkin(10) + t332 * t105 / 0.2e1 + (t286 * t30 - t3 * t332 + t397 * t71) * mrSges(7,1) + (t254 * t356 + t488) * mrSges(5,3) + m(5) * (-pkin(3) * t150 + pkin(10) * t488) + t511 * t68 + t512 * t91 + t513 * t92 + (t16 * t513 + t17 * t512 + t2 * t219 + t218 * t3 + t30 * t305 + t511 * t71) * m(7) + t257 * t506 - pkin(3) * t123 - t111 * t126 + (-Ifges(7,5) * t287 - Ifges(7,6) * t286 + t328 * t359 - t332 * t515) * t464 - t156 * (mrSges(5,1) * t328 + mrSges(5,2) * t332) * t254 - t257 * t507 - (t140 * t254 + t516) * t332 / 0.2e1 - t132 * t191 - t131 * t192 + t218 * t28 + t219 * t29 + t364 * t420 + t196 * t450 + (Ifges(5,5) * t328 + Ifges(5,6) * t332) * t453 + (Ifges(7,5) * t230 + Ifges(7,6) * t231) * t457 + (Ifges(5,2) * t332 + t440) * t465 + (Ifges(5,1) * t328 + t439) * t466 + (Ifges(7,1) * t230 + Ifges(7,4) * t231) * t469 + (Ifges(7,4) * t230 + Ifges(7,2) * t231) * t471 + t328 * t473 + (-Ifges(6,6) * t332 + t328 * t361) * t474 + (-Ifges(6,5) * t332 + t328 * t363) * t475 + t400 * t479 - t287 * t482 - t286 * t483 - t215 * (mrSges(4,1) * t257 + mrSges(4,2) * t254) - t280 * (Ifges(4,5) * t254 - Ifges(4,6) * t257) / 0.2e1 + t281 * t75 + t282 * t76 + t305 * t15; -t50 * t365 + t484 + (Ifges(7,5) * t300 - Ifges(7,6) * t352 + t358) * t464 + (t16 * t394 - t17 * t395 - t2 * t352 - t3 * t300) * mrSges(7,3) + (Ifges(7,4) * t300 - Ifges(7,2) * t352) * t480 + (Ifges(7,1) * t300 - Ifges(7,4) * t352) * t481 + t30 * (mrSges(7,1) * t352 + mrSges(7,2) * t300) - t352 * t483 + (t159 / 0.2e1 - t264 / 0.2e1) * t64 + (-Ifges(7,5) * t264 - Ifges(7,6) * t265) * t457 + (-Ifges(7,1) * t264 - Ifges(7,4) * t265) * t469 + (-Ifges(7,4) * t264 - Ifges(7,2) * t265) * t471 + (t158 / 0.2e1 - t265 / 0.2e1) * t63 + t377 + (-Ifges(7,5) * t159 - Ifges(7,6) * t158) * t458 + (-Ifges(7,1) * t159 - Ifges(7,4) * t158) * t470 + (-Ifges(7,4) * t159 - Ifges(7,2) * t158) * t472 + t335 * t225 + (-pkin(4) * t50 - t46 * t72 - t47 * t73 - t81 * t89) * m(6) - t398 * t89 + (mrSges(7,1) * t395 - mrSges(7,2) * t394) * t71 + (t445 * t496 - t521) * qJD(5) - pkin(4) * t55 - t80 * t68 + (-t327 * t75 + t331 * t76 + m(6) * t487 + (-m(6) * t357 - t327 * t136 - t331 * t137) * qJD(5)) * pkin(11) + t487 * mrSges(6,3) + ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t225 + t510) * t224 + t497 * t91 + t498 * t92 + (t16 * t498 + t17 * t497 + t2 * t275 + t274 * t3 + t30 * t319 - t71 * t80) * m(7) - t73 * t136 - t72 * t137 - t88 * t191 + t43 * t448 + t360 * t474 + t362 * t475 + t327 * t479 + t300 * t482 + t274 * t28 + t275 * t29 + t319 * t15; t42 + (-Ifges(6,2) * t190 + t186 + t97) * t463 + t122 * t508 + (t189 * t46 + t190 * t47) * mrSges(6,3) - m(7) * (t16 * t18 + t17 * t19) - t19 * t91 - t18 * t92 + (t330 * t28 + t326 * t29 + m(7) * (t2 * t326 + t3 * t330) - t496 * t190 + (-t326 * t92 + t330 * t91 + m(7) * (-t16 * t326 + t17 * t330)) * qJD(6)) * pkin(5) - t46 * t136 + t47 * t137 - t81 * (mrSges(6,1) * t190 + mrSges(6,2) * t189) + (Ifges(6,5) * t189 - Ifges(6,6) * t190) * t456 + t96 * t460 + (Ifges(6,1) * t189 - t438) * t461 + t518 + t519; -t16 * t91 + t17 * t92 + t63 * t469 + t518;];
tauc  = t1(:);
