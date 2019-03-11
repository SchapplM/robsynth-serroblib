% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:53
% EndTime: 2019-03-09 19:25:47
% DurationCPUTime: 28.25s
% Computational Cost: add. (17561->878), mult. (44641->1198), div. (0->0), fcn. (34009->10), ass. (0->367)
t320 = cos(pkin(6));
t328 = cos(qJ(2));
t309 = t320 * t328 * pkin(1);
t324 = sin(qJ(2));
t319 = sin(pkin(6));
t394 = qJD(1) * t319;
t375 = t324 * t394;
t252 = -pkin(8) * t375 + qJD(1) * t309;
t339 = (pkin(2) * t324 - pkin(9) * t328) * t319;
t253 = qJD(1) * t339;
t323 = sin(qJ(3));
t327 = cos(qJ(3));
t167 = -t323 * t252 + t253 * t327;
t476 = pkin(9) - pkin(10);
t296 = t476 * t327;
t329 = -pkin(3) - pkin(4);
t530 = -(-pkin(10) * t327 * t328 + t324 * t329) * t394 + t167 + qJD(3) * t296;
t168 = t327 * t252 + t323 * t253;
t150 = qJ(4) * t375 + t168;
t374 = t328 * t394;
t366 = t323 * t374;
t390 = qJD(3) * t323;
t529 = pkin(10) * t366 + t476 * t390 + t150;
t322 = sin(qJ(5));
t326 = cos(qJ(5));
t268 = t322 * t323 + t326 * t327;
t195 = (qJD(3) - qJD(5)) * t268;
t210 = t268 * t374;
t399 = t195 - t210;
t387 = qJD(5) * t326;
t388 = qJD(5) * t322;
t389 = qJD(3) * t327;
t196 = t322 * t389 + t323 * t387 - t326 * t390 - t327 * t388;
t365 = t327 * t374;
t209 = t322 * t365 - t326 * t366;
t398 = t196 - t209;
t295 = t476 * t323;
t341 = t326 * t295 - t296 * t322;
t520 = qJD(5) * t341 + t322 * t530 - t529 * t326;
t280 = qJ(4) * t365;
t450 = pkin(1) * t324;
t308 = t320 * t450;
t380 = t329 * t323;
t395 = qJ(4) * t389 + t323 * qJD(4);
t402 = t319 * t328;
t516 = qJD(3) * t380 + t395 - t280 - (-t308 + (-pkin(8) + t380) * t402) * qJD(1);
t512 = Ifges(4,1) + Ifges(5,1);
t511 = -Ifges(4,4) + Ifges(5,5);
t510 = Ifges(5,4) + Ifges(4,5);
t528 = t398 * pkin(5) - pkin(11) * t399 + t516;
t527 = pkin(11) * t375 + t520;
t526 = Ifges(5,2) + Ifges(4,3);
t509 = Ifges(4,6) - Ifges(5,6);
t212 = t295 * t322 + t296 * t326;
t519 = -qJD(5) * t212 + t529 * t322 + t326 * t530;
t307 = pkin(8) * t402;
t208 = qJD(2) * pkin(9) + (t307 + (pkin(9) + t450) * t320) * qJD(1);
t245 = (-pkin(2) * t328 - pkin(9) * t324 - pkin(1)) * t319;
t220 = qJD(1) * t245;
t142 = -t323 * t208 + t327 * t220;
t368 = qJD(1) * t320 + qJD(2);
t233 = t323 * t368 + t327 * t375;
t112 = pkin(10) * t233 + t142;
t525 = qJD(4) - t112;
t232 = t323 * t375 - t327 * t368;
t155 = t232 * t326 - t322 * t233;
t393 = qJD(2) * t319;
t370 = qJD(1) * t393;
t363 = t328 * t370;
t185 = -qJD(3) * t232 + t327 * t363;
t186 = qJD(3) * t233 + t323 * t363;
t77 = qJD(5) * t155 + t185 * t326 + t186 * t322;
t524 = t77 / 0.2e1;
t157 = t232 * t322 + t233 * t326;
t78 = qJD(5) * t157 + t185 * t322 - t326 * t186;
t523 = -t78 / 0.2e1;
t285 = -t327 * pkin(3) - t323 * qJ(4) - pkin(2);
t266 = t327 * pkin(4) - t285;
t269 = -t322 * t327 + t323 * t326;
t169 = pkin(5) * t268 - pkin(11) * t269 + t266;
t321 = sin(qJ(6));
t325 = cos(qJ(6));
t116 = t169 * t325 - t212 * t321;
t522 = qJD(6) * t116 + t321 * t528 + t325 * t527;
t117 = t169 * t321 + t212 * t325;
t521 = -qJD(6) * t117 - t321 * t527 + t325 * t528;
t518 = -pkin(5) * t375 - t519;
t364 = t324 * t370;
t517 = t185 * t512 + t511 * t186 + t510 * t364;
t396 = t307 + t308;
t495 = t396 * qJD(2);
t242 = qJD(1) * t495;
t83 = t186 * pkin(3) - t185 * qJ(4) - t233 * qJD(4) + t242;
t403 = t319 * t324;
t262 = t320 * t323 + t327 * t403;
t391 = qJD(2) * t328;
t372 = t319 * t391;
t199 = qJD(3) * t262 + t323 * t372;
t261 = -t320 * t327 + t323 * t403;
t200 = -qJD(3) * t261 + t327 * t372;
t99 = t199 * pkin(3) - t200 * qJ(4) - t262 * qJD(4) + t495;
t290 = qJD(3) - t374;
t493 = t290 - qJD(5);
t353 = Ifges(7,5) * t325 - Ifges(7,6) * t321;
t433 = Ifges(7,4) * t325;
t354 = -Ifges(7,2) * t321 + t433;
t434 = Ifges(7,4) * t321;
t355 = Ifges(7,1) * t325 - t434;
t356 = mrSges(7,1) * t321 + mrSges(7,2) * t325;
t279 = t290 * qJ(4);
t143 = t327 * t208 + t323 * t220;
t359 = pkin(10) * t232 + t143;
t101 = t279 + t359;
t95 = t290 * t329 + t525;
t42 = -t101 * t322 + t326 * t95;
t40 = pkin(5) * t493 - t42;
t384 = -qJD(6) + t155;
t467 = -t384 / 0.2e1;
t129 = t157 * t325 - t321 * t493;
t470 = t129 / 0.2e1;
t128 = -t157 * t321 - t325 * t493;
t472 = t128 / 0.2e1;
t515 = t353 * t467 + t354 * t472 + t355 * t470 + t40 * t356;
t97 = pkin(5) * t157 - pkin(11) * t155;
t514 = -Ifges(6,5) / 0.2e1;
t479 = t78 / 0.2e1;
t513 = -t77 * Ifges(6,4) / 0.2e1;
t132 = -mrSges(6,1) * t493 - mrSges(6,3) * t157;
t69 = -mrSges(7,1) * t128 + mrSges(7,2) * t129;
t408 = t132 - t69;
t282 = -t322 * qJ(4) + t326 * t329;
t247 = t326 * qJD(4) + qJD(5) * t282;
t54 = t326 * t112 + t322 * t359;
t507 = t247 - t54;
t283 = t326 * qJ(4) + t322 * t329;
t506 = qJD(5) * t283 + t322 * t525 + t326 * t359;
t244 = pkin(9) * t320 + t396;
t165 = -t323 * t244 + t245 * t327;
t146 = pkin(3) * t402 - t165;
t113 = pkin(4) * t402 - pkin(10) * t262 + t146;
t166 = t327 * t244 + t323 * t245;
t145 = -qJ(4) * t402 + t166;
t122 = pkin(10) * t261 + t145;
t502 = t322 * t113 + t326 * t122;
t225 = Ifges(4,4) * t232;
t418 = t232 * Ifges(5,5);
t501 = t233 * t512 + t510 * t290 - t225 + t418;
t500 = t185 * t510 - t186 * t509 + t364 * t526;
t254 = qJD(2) * t339;
t240 = qJD(1) * t254;
t264 = -pkin(8) * t403 + t309;
t256 = t264 * qJD(2);
t241 = qJD(1) * t256;
t89 = -t208 * t390 + t220 * t389 + t323 * t240 + t327 * t241;
t90 = -t208 * t389 - t220 * t390 + t240 * t327 - t323 * t241;
t499 = -t323 * t90 + t327 * t89;
t79 = qJ(4) * t364 + t290 * qJD(4) + t89;
t84 = -pkin(3) * t364 - t90;
t498 = t323 * t84 + t327 * t79;
t43 = t101 * t326 + t322 * t95;
t41 = -pkin(11) * t493 + t43;
t207 = -t368 * pkin(2) - t252;
t120 = t232 * pkin(3) - t233 * qJ(4) + t207;
t100 = -pkin(4) * t232 - t120;
t44 = -pkin(5) * t155 - pkin(11) * t157 + t100;
t17 = -t321 * t41 + t325 * t44;
t18 = t321 * t44 + t325 * t41;
t497 = -t17 * t321 + t18 * t325;
t496 = -t100 * mrSges(6,2) + t42 * mrSges(6,3);
t494 = qJD(4) - t142;
t149 = Ifges(6,4) * t155;
t452 = t325 / 0.2e1;
t454 = -t321 / 0.2e1;
t435 = Ifges(7,4) * t129;
t46 = Ifges(7,2) * t128 - Ifges(7,6) * t384 + t435;
t124 = Ifges(7,4) * t128;
t47 = Ifges(7,1) * t129 - Ifges(7,5) * t384 + t124;
t340 = t452 * t47 + t454 * t46;
t415 = t493 * Ifges(6,5);
t422 = t157 * Ifges(6,1);
t87 = t149 - t415 + t422;
t477 = t87 / 0.2e1;
t492 = t477 - t415 / 0.2e1 + t340 + t149 / 0.2e1 + t515;
t491 = Ifges(6,1) * t524 + Ifges(6,4) * t523;
t55 = pkin(10) * t186 + t79;
t392 = qJD(2) * t324;
t373 = t319 * t392;
t346 = t329 * t373;
t56 = -pkin(10) * t185 + qJD(1) * t346 - t90;
t11 = -qJD(5) * t43 - t322 * t55 + t326 * t56;
t103 = -t244 * t389 - t245 * t390 + t254 * t327 - t323 * t256;
t72 = -pkin(10) * t200 - t103 + t346;
t102 = -t244 * t390 + t245 * t389 + t323 * t254 + t327 * t256;
t91 = qJ(4) * t373 - qJD(4) * t402 + t102;
t73 = pkin(10) * t199 + t91;
t15 = -qJD(5) * t502 - t322 * t73 + t326 * t72;
t224 = Ifges(5,5) * t233;
t136 = t290 * Ifges(5,6) + t232 * Ifges(5,3) + t224;
t343 = t142 * t327 + t143 * t323;
t125 = -pkin(3) * t290 + t494;
t127 = t279 + t143;
t344 = t125 * t327 - t127 * t323;
t430 = Ifges(5,5) * t327;
t431 = Ifges(5,5) * t323;
t437 = Ifges(4,4) * t327;
t438 = Ifges(4,4) * t323;
t451 = t327 / 0.2e1;
t453 = t323 / 0.2e1;
t455 = t290 / 0.2e1;
t458 = t233 / 0.2e1;
t460 = t232 / 0.2e1;
t461 = -t232 / 0.2e1;
t417 = t233 * Ifges(4,4);
t139 = -t232 * Ifges(4,2) + t290 * Ifges(4,6) + t417;
t469 = -t139 / 0.2e1;
t490 = t344 * mrSges(5,2) - t343 * mrSges(4,3) + (-Ifges(4,2) * t323 + t437) * t461 + (Ifges(5,3) * t323 + t430) * t460 + t207 * (mrSges(4,1) * t323 + mrSges(4,2) * t327) + t120 * (mrSges(5,1) * t323 - mrSges(5,3) * t327) + t136 * t453 + t323 * t469 + (t327 * t512 + t431 - t438) * t458 + (-t323 * t509 + t327 * t510) * t455 + t501 * t451;
t38 = qJD(6) * t128 - t321 * t364 + t325 * t77;
t39 = -qJD(6) * t129 - t321 * t77 - t325 * t364;
t7 = t38 * Ifges(7,1) + t39 * Ifges(7,4) + t78 * Ifges(7,5);
t489 = t7 / 0.2e1;
t488 = Ifges(6,2) / 0.2e1;
t487 = Ifges(6,2) * t479 + Ifges(6,6) * t364 / 0.2e1 + t513;
t486 = t364 * t514 + t491;
t485 = t38 / 0.2e1;
t484 = t39 / 0.2e1;
t425 = t384 * Ifges(7,3);
t426 = t129 * Ifges(7,5);
t427 = t128 * Ifges(7,6);
t45 = -t425 + t426 + t427;
t483 = t45 / 0.2e1;
t482 = -t46 / 0.2e1;
t481 = t46 / 0.2e1;
t480 = -t47 / 0.2e1;
t414 = t493 * Ifges(6,6);
t424 = t155 * Ifges(6,2);
t436 = Ifges(6,4) * t157;
t86 = -t414 + t424 + t436;
t478 = -t86 / 0.2e1;
t475 = pkin(1) * mrSges(3,1);
t474 = pkin(1) * mrSges(3,2);
t473 = -t128 / 0.2e1;
t471 = -t129 / 0.2e1;
t468 = t384 / 0.2e1;
t342 = t261 * t326 - t262 * t322;
t466 = -t342 / 0.2e1;
t465 = t185 / 0.2e1;
t464 = -t186 / 0.2e1;
t463 = t186 / 0.2e1;
t37 = Ifges(7,5) * t38;
t36 = Ifges(7,6) * t39;
t445 = -qJD(2) / 0.2e1;
t16 = -mrSges(7,1) * t39 + mrSges(7,2) * t38;
t65 = -mrSges(6,1) * t364 - mrSges(6,3) * t77;
t443 = t16 - t65;
t442 = Ifges(6,5) * t77 - Ifges(6,6) * t78;
t441 = mrSges(4,3) * t232;
t440 = mrSges(4,3) * t233;
t439 = Ifges(3,4) * t324;
t432 = Ifges(3,5) * t320;
t429 = Ifges(3,6) * t320;
t423 = t155 * Ifges(6,6);
t421 = t157 * Ifges(6,5);
t416 = t241 * mrSges(3,2);
t413 = t493 * Ifges(6,3);
t406 = Ifges(3,6) * qJD(2);
t404 = t290 * t326;
t187 = -mrSges(4,2) * t290 - t441;
t190 = -mrSges(5,2) * t232 + mrSges(5,3) * t290;
t401 = -t187 - t190;
t188 = mrSges(4,1) * t290 - t440;
t189 = -mrSges(5,1) * t290 + mrSges(5,2) * t233;
t400 = -t188 + t189;
t162 = t233 * pkin(3) + t232 * qJ(4);
t397 = -mrSges(3,1) * t368 + mrSges(4,1) * t232 + mrSges(4,2) * t233 + mrSges(3,3) * t375;
t386 = qJD(6) * t321;
t385 = qJD(6) * t325;
t5 = Ifges(7,3) * t78 + t36 + t37;
t383 = Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t382 = Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t381 = -Ifges(4,3) / 0.2e1 - Ifges(5,2) / 0.2e1;
t243 = -t320 * pkin(2) - t264;
t59 = -pkin(4) * t186 - t83;
t19 = pkin(5) * t78 - pkin(11) * t77 + t59;
t10 = -t101 * t388 + t322 * t56 + t326 * t55 + t95 * t387;
t8 = -pkin(11) * t364 + t10;
t1 = qJD(6) * t17 + t19 * t321 + t325 * t8;
t2 = -qJD(6) * t18 + t19 * t325 - t321 * t8;
t360 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t126 = -pkin(4) * t233 - t162;
t358 = t1 * t325 - t2 * t321;
t357 = mrSges(7,1) * t325 - mrSges(7,2) * t321;
t294 = Ifges(7,1) * t321 + t433;
t293 = Ifges(7,2) * t325 + t434;
t292 = Ifges(7,5) * t321 + Ifges(7,6) * t325;
t352 = t17 * t325 + t18 * t321;
t20 = mrSges(7,1) * t78 - mrSges(7,3) * t38;
t21 = -mrSges(7,2) * t78 + mrSges(7,3) * t39;
t350 = -t321 * t20 + t325 * t21;
t58 = pkin(11) * t402 + t502;
t144 = t261 * pkin(3) - t262 * qJ(4) + t243;
t121 = -pkin(4) * t261 - t144;
t177 = t261 * t322 + t262 * t326;
t62 = -pkin(5) * t342 - pkin(11) * t177 + t121;
t26 = t321 * t62 + t325 * t58;
t25 = -t321 * t58 + t325 * t62;
t60 = t113 * t326 - t122 * t322;
t147 = -t177 * t321 + t325 * t402;
t148 = t177 * t325 + t321 * t402;
t14 = t113 * t387 - t122 * t388 + t322 * t72 + t326 * t73;
t160 = -mrSges(5,1) * t364 + t185 * mrSges(5,2);
t336 = -qJD(6) * t352 + t358;
t80 = -pkin(4) * t199 - t99;
t335 = t483 + t478 - t436 / 0.2e1 + t427 / 0.2e1 + t426 / 0.2e1 - t425 / 0.2e1 + t414 / 0.2e1;
t334 = -t100 * mrSges(6,1) - t17 * mrSges(7,1) + t18 * mrSges(7,2) + t43 * mrSges(6,3) - t335;
t9 = pkin(5) * t364 - t11;
t333 = t11 * mrSges(6,1) - t10 * mrSges(6,2) - Ifges(6,3) * t364 + qJD(6) * t515 + t292 * t479 + t293 * t484 + t294 * t485 - t357 * t9 + t442;
t331 = t422 / 0.2e1 + t492;
t299 = Ifges(3,4) * t374;
t289 = Ifges(3,5) * t363;
t276 = -pkin(11) + t283;
t275 = pkin(5) - t282;
t260 = pkin(3) * t390 - t395;
t255 = t396 * qJD(1);
t250 = -mrSges(3,2) * t368 + mrSges(3,3) * t374;
t204 = Ifges(3,1) * t375 + Ifges(3,5) * t368 + t299;
t203 = t406 + (t429 + (Ifges(3,2) * t328 + t439) * t319) * qJD(1);
t174 = -t280 + (t308 + (pkin(3) * t323 + pkin(8)) * t402) * qJD(1);
t173 = t210 * t325 - t321 * t375;
t172 = -t210 * t321 - t325 * t375;
t171 = t233 * t321 + t325 * t404;
t170 = t233 * t325 - t321 * t404;
t163 = mrSges(5,1) * t232 - mrSges(5,3) * t233;
t161 = -mrSges(4,2) * t364 - mrSges(4,3) * t186;
t159 = mrSges(4,1) * t364 - mrSges(4,3) * t185;
t158 = -mrSges(5,2) * t186 + mrSges(5,3) * t364;
t153 = -pkin(3) * t375 - t167;
t138 = t233 * Ifges(5,4) + t290 * Ifges(5,2) + t232 * Ifges(5,6);
t137 = t233 * Ifges(4,5) - t232 * Ifges(4,6) + t290 * Ifges(4,3);
t131 = mrSges(6,2) * t493 + mrSges(6,3) * t155;
t115 = mrSges(4,1) * t186 + mrSges(4,2) * t185;
t114 = mrSges(5,1) * t186 - mrSges(5,3) * t185;
t105 = t185 * Ifges(4,4) - t186 * Ifges(4,2) + Ifges(4,6) * t364;
t104 = t185 * Ifges(5,5) + Ifges(5,6) * t364 + t186 * Ifges(5,3);
t98 = -pkin(3) * t373 - t103;
t96 = -mrSges(6,1) * t155 + mrSges(6,2) * t157;
t94 = qJD(5) * t342 + t199 * t322 + t200 * t326;
t93 = qJD(5) * t177 - t199 * t326 + t200 * t322;
t85 = -t413 + t421 + t423;
t82 = -mrSges(7,1) * t384 - mrSges(7,3) * t129;
t81 = mrSges(7,2) * t384 + mrSges(7,3) * t128;
t66 = mrSges(6,2) * t364 - mrSges(6,3) * t78;
t57 = -pkin(5) * t402 - t60;
t50 = t126 - t97;
t49 = qJD(6) * t147 - t321 * t373 + t325 * t94;
t48 = -qJD(6) * t148 - t321 * t94 - t325 * t373;
t31 = mrSges(6,1) * t78 + mrSges(6,2) * t77;
t28 = t321 * t97 + t325 * t42;
t27 = -t321 * t42 + t325 * t97;
t24 = pkin(5) * t93 - pkin(11) * t94 + t80;
t23 = t321 * t50 + t325 * t54;
t22 = -t321 * t54 + t325 * t50;
t13 = pkin(5) * t373 - t15;
t12 = -pkin(11) * t373 + t14;
t6 = Ifges(7,4) * t38 + Ifges(7,2) * t39 + Ifges(7,6) * t78;
t4 = -qJD(6) * t26 - t12 * t321 + t24 * t325;
t3 = qJD(6) * t25 + t12 * t325 + t24 * t321;
t29 = [-t493 * (Ifges(6,5) * t94 - Ifges(6,6) * t93 - Ifges(6,3) * t373) / 0.2e1 + (t328 * t204 + (t138 + t137) * t324) * t393 / 0.2e1 + m(6) * (t10 * t502 + t100 * t80 + t11 * t60 + t121 * t59 + t14 * t43 + t15 * t42) + t502 * t66 - t342 * t487 + t10 * (-mrSges(6,2) * t402 + mrSges(6,3) * t342) + t1 * (mrSges(7,2) * t342 + mrSges(7,3) * t147) + t2 * (-mrSges(7,1) * t342 - mrSges(7,3) * t148) + t59 * (-mrSges(6,1) * t342 + mrSges(6,2) * t177) + (Ifges(7,4) * t148 + Ifges(7,2) * t147 - Ifges(7,6) * t342) * t484 + (Ifges(7,1) * t148 + Ifges(7,4) * t147 - Ifges(7,5) * t342) * t485 + (Ifges(7,5) * t148 + Ifges(7,6) * t147 - Ifges(7,3) * t342) * t479 + t397 * t495 + m(3) * (t241 * t396 - t242 * t264 - t252 * t495 + t255 * t256) + m(4) * (t102 * t143 + t103 * t142 + t165 * t90 + t166 * t89 + t207 * t495 + t242 * t243) + (t510 * t200 + t373 * t526) * t455 + (-t120 * t200 + t127 * t373 - t262 * t83 - t402 * t79) * mrSges(5,3) + ((-t264 * mrSges(3,3) + t432 + (-0.2e1 * t474 + 0.3e1 / 0.2e1 * Ifges(3,4) * t328) * t319) * t391 + (-0.3e1 / 0.2e1 * t429 + t177 * t514 + Ifges(6,6) * t466 - t396 * mrSges(3,3) + (-0.2e1 * t475 - 0.3e1 / 0.2e1 * t439) * t319 + t383 * t262 + (-Ifges(6,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + t381) * t402) * t392) * t394 + t90 * (-mrSges(4,1) * t402 - mrSges(4,3) * t262) + t155 * (Ifges(6,4) * t94 - Ifges(6,2) * t93 - Ifges(6,6) * t373) / 0.2e1 + t157 * (Ifges(6,1) * t94 - Ifges(6,4) * t93 - Ifges(6,5) * t373) / 0.2e1 + (-t416 + t289 / 0.2e1 - t242 * mrSges(3,1)) * t320 + (-t143 * mrSges(4,3) - t127 * mrSges(5,2) + t207 * mrSges(4,1) + t120 * mrSges(5,1) + t136 / 0.2e1 + t469 + Ifges(5,3) * t460 - Ifges(4,2) * t461 + t511 * t458 - t509 * t455) * t199 + (-t79 * mrSges(5,2) - t89 * mrSges(4,3) + t382 * t392 * t394 - t105 / 0.2e1 + mrSges(4,1) * t242 + t83 * mrSges(5,1) + t104 / 0.2e1 + Ifges(5,3) * t463 - Ifges(4,2) * t464 + t511 * t465) * t261 + (t200 * t512 + t373 * t510) * t458 + (t262 * t512 - t402 * t510) * t465 + (t241 * t328 + t242 * t324 + (-t252 * t328 - t255 * t324) * qJD(2)) * t319 * mrSges(3,3) - (t203 + t85) * t373 / 0.2e1 + (Ifges(5,5) * t200 + Ifges(5,6) * t373) * t460 + (Ifges(5,5) * t262 - Ifges(5,6) * t402) * t463 + t11 * (mrSges(6,1) * t402 - mrSges(6,3) * t177) + t84 * (mrSges(5,1) * t402 + mrSges(5,2) * t262) + m(5) * (t120 * t99 + t125 * t98 + t127 * t91 + t144 * t83 + t145 * t79 + t146 * t84) + m(7) * (t1 * t26 + t13 * t40 + t17 * t4 + t18 * t3 + t2 * t25 + t57 * t9) + t42 * (-mrSges(6,1) * t373 - mrSges(6,3) * t94) + t125 * (-mrSges(5,1) * t373 + mrSges(5,2) * t200) + t43 * (mrSges(6,2) * t373 - mrSges(6,3) * t93) + t142 * (mrSges(4,1) * t373 - mrSges(4,3) * t200) - t500 * t402 / 0.2e1 + t501 * t200 / 0.2e1 + (-t143 * t373 + t200 * t207 + t242 * t262 + t402 * t89) * mrSges(4,2) + t442 * t402 / 0.2e1 + t256 * t250 + t243 * t115 + t102 * t187 + t103 * t188 + t98 * t189 + t91 * t190 + t146 * t160 + t99 * t163 + t165 * t159 + t166 * t161 + t145 * t158 + t144 * t114 + t147 * t6 / 0.2e1 + t9 * (-mrSges(7,1) * t147 + mrSges(7,2) * t148) + t14 * t131 + t15 * t132 + t121 * t31 + t80 * t96 + t100 * (mrSges(6,1) * t93 + mrSges(6,2) * t94) + t17 * (mrSges(7,1) * t93 - mrSges(7,3) * t49) + t18 * (-mrSges(7,2) * t93 + mrSges(7,3) * t48) + t3 * t81 + t4 * t82 + t60 * t65 + t13 * t69 + t57 * t16 + t49 * t47 / 0.2e1 + t40 * (-mrSges(7,1) * t48 + mrSges(7,2) * t49) + (Ifges(6,4) * t177 + Ifges(6,2) * t342 + Ifges(6,6) * t402) * t523 + (Ifges(6,1) * t177 + Ifges(6,4) * t342 + Ifges(6,5) * t402) * t524 + t48 * t481 + t93 * t483 + t177 * t486 + t148 * t489 + t5 * t466 + (Ifges(7,5) * t49 + Ifges(7,6) * t48 + Ifges(7,3) * t93) * t467 + (Ifges(7,1) * t49 + Ifges(7,4) * t48 + Ifges(7,5) * t93) * t470 + (Ifges(7,4) * t49 + Ifges(7,2) * t48 + Ifges(7,6) * t93) * t472 + t94 * t477 + t93 * t478 + t25 * t20 + t26 * t21 + qJD(2) ^ 2 * (Ifges(3,5) * t328 - Ifges(3,6) * t324) * t319 / 0.2e1 + (Ifges(4,4) * t200 + Ifges(4,6) * t373) * t461 + (Ifges(4,4) * t262 - Ifges(4,6) * t402) * t464 + t517 * t262 / 0.2e1; t493 * (Ifges(6,5) * t210 - Ifges(6,6) * t209) / 0.2e1 + (t353 * t479 + t354 * t484 + t355 * t485 + t9 * t356 + t59 * mrSges(6,2) + t486 + t7 * t452 + t6 * t454 - t11 * mrSges(6,3) + (-t1 * t321 - t2 * t325) * mrSges(7,3) + (-mrSges(7,3) * t497 + t292 * t468 + t293 * t473 + t294 * t471 + t325 * t482 + t357 * t40 + t454 * t47) * qJD(6) + t491) * t269 - t443 * t341 + ((t203 / 0.2e1 - t138 / 0.2e1 - t137 / 0.2e1 + t125 * mrSges(5,1) - t142 * mrSges(4,1) - t127 * mrSges(5,3) + t143 * mrSges(4,2) - t413 / 0.2e1 + t421 / 0.2e1 + t423 / 0.2e1 + t42 * mrSges(6,1) - t43 * mrSges(6,2) + (Ifges(6,5) * t269 - Ifges(6,6) * t268) * t445 + t255 * mrSges(3,3) + t85 / 0.2e1 - t406 / 0.2e1 + t381 * t290 - t383 * t233 - t382 * t232 + (t429 / 0.2e1 + (t475 + t439 / 0.2e1) * t319) * qJD(1) + (t323 * t510 + t327 * t509) * qJD(2) / 0.2e1) * t324 + (-t299 / 0.2e1 - t204 / 0.2e1 + Ifges(3,5) * t445 + ((t474 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t324) * t319 - t432 / 0.2e1) * qJD(1) + t252 * mrSges(3,3) - t490) * t328) * t394 + (-mrSges(4,1) * t327 + mrSges(4,2) * t323 - mrSges(3,1)) * t242 - t416 + (-t424 / 0.2e1 + t335) * t196 + (-t10 * mrSges(6,3) + t513 + t36 / 0.2e1 + t37 / 0.2e1 + t5 / 0.2e1 + t487 + t59 * mrSges(6,1) + (t488 + Ifges(7,3) / 0.2e1) * t78 + t360) * t268 - m(4) * (t142 * t167 + t143 * t168 + t207 * t255) - m(5) * (t120 * t174 + t125 * t153 + t127 * t150) + t289 - t397 * t255 + ((-t325 * t195 + t173) * mrSges(7,3) + t398 * mrSges(7,1)) * t17 + ((-t321 * t195 - t172) * mrSges(7,3) - t398 * mrSges(7,2)) * t18 + (mrSges(6,1) * t398 + mrSges(6,2) * t399) * t100 + (-t398 * t43 - t399 * t42) * mrSges(6,3) + (t260 - t174) * t163 - t155 * (Ifges(6,4) * t210 - Ifges(6,2) * t209) / 0.2e1 + ((-m(4) * t343 + m(5) * t344 + t323 * t401 + t327 * t400) * pkin(9) + t490) * qJD(3) + ((t158 + t161) * t327 + (-t159 + t160) * t323) * pkin(9) + m(5) * (pkin(9) * t498 + t120 * t260 + t285 * t83) + t498 * mrSges(5,2) + t499 * mrSges(4,3) + m(4) * (-pkin(2) * t242 + pkin(9) * t499) + t83 * (-mrSges(5,1) * t327 - mrSges(5,3) * t323) - t327 * t104 / 0.2e1 + t331 * t195 + t285 * t114 + t266 * t31 - t252 * t250 + t212 * t66 - t210 * t87 / 0.2e1 - t209 * t45 / 0.2e1 + t209 * t86 / 0.2e1 - t157 * (Ifges(6,1) * t210 - Ifges(6,4) * t209) / 0.2e1 - t168 * t187 - t167 * t188 - t153 * t189 - t150 * t190 - t40 * (-mrSges(7,1) * t172 + mrSges(7,2) * t173) - pkin(2) * t115 + t116 * t20 + t117 * t21 + t173 * t480 + t172 * t482 + (Ifges(7,5) * t173 + Ifges(7,6) * t172 + Ifges(7,3) * t209) * t468 + (Ifges(7,1) * t173 + Ifges(7,4) * t172 + Ifges(7,5) * t209) * t471 + (Ifges(7,4) * t173 + Ifges(7,2) * t172 + Ifges(7,6) * t209) * t473 + (-Ifges(5,3) * t327 + t431) * t463 + (Ifges(4,2) * t327 + t438) * t464 + t105 * t451 + t516 * t96 + t517 * t453 + t518 * t69 + t519 * t132 + t520 * t131 + (t10 * t212 + t100 * t516 + t11 * t341 + t266 * t59 + t42 * t519 + t43 * t520) * m(6) + t521 * t82 + t522 * t81 + (t1 * t117 + t116 * t2 + t17 * t521 + t18 * t522 - t341 * t9 + t40 * t518) * m(7) + (t323 * t512 - t430 + t437) * t465; t500 - (t424 / 0.2e1 + t334) * t157 - (mrSges(7,3) * t352 - t331 + t496) * t155 - t408 * t506 + (-t17 * t22 - t18 * t23 + t247 * t497 + t275 * t9 + t276 * t336 + t40 * t506) * m(7) + t507 * t131 + (t10 * t283 - t100 * t126 + t11 * t282 - t42 * t506 + t43 * t507) * m(6) + ((t17 * mrSges(7,3) - t276 * t82 + t480) * t325 + (t18 * mrSges(7,3) - t276 * t81 + t481) * t321) * qJD(6) + (-t6 / 0.2e1 - t1 * mrSges(7,3) + t276 * t21 + t247 * t81) * t325 - (-t232 * t510 - t233 * t509) * t290 / 0.2e1 - (-t512 * t232 + t136 + t224 - t417) * t233 / 0.2e1 + (-t7 / 0.2e1 + t2 * mrSges(7,3) - t276 * t20 - t247 * t82) * t321 + (t125 * t232 + t127 * t233) * mrSges(5,2) + (-Ifges(4,2) * t233 - t225 + t501) * t460 + (-pkin(3) * t84 + qJ(4) * t79 - t120 * t162 - t125 * t143 + t127 * t494) * m(5) + (-t400 + t440) * t143 + (t401 - t441) * t142 + t275 * t16 + t282 * t65 + t283 * t66 - t120 * (mrSges(5,1) * t233 + mrSges(5,3) * t232) - t207 * (mrSges(4,1) * t233 - mrSges(4,2) * t232) + qJD(4) * t190 - t162 * t163 + qJ(4) * t158 - pkin(3) * t160 - t126 * t96 - t89 * mrSges(4,2) + t90 * mrSges(4,1) - t23 * t81 - t22 * t82 - t84 * mrSges(5,1) + t79 * mrSges(5,3) + t139 * t458 + (Ifges(5,3) * t233 - t418) * t461 - t333; -t170 * t82 - t171 * t81 - t290 * t190 + (t163 - t96) * t233 + (-t290 * t131 + (-t321 * t82 + t325 * t81 + t131) * qJD(5) - t443) * t326 + (t66 + (-t321 * t81 - t325 * t82) * qJD(6) + t350 + t493 * t408) * t322 + t160 + ((qJD(5) * t497 - t9) * t326 - t17 * t170 - t171 * t18 + (-t40 * t493 + t336) * t322) * m(7) + (t10 * t322 - t100 * t233 + t11 * t326 - t493 * (-t322 * t42 + t326 * t43)) * m(6) + (t120 * t233 - t127 * t290 + t84) * m(5); ((t17 * t384 + t1) * t325 + (t18 * t384 - t2) * t321) * mrSges(7,3) + t334 * t157 + t340 * qJD(6) + t6 * t452 + t408 * t43 + ((-Ifges(6,1) / 0.2e1 + t488) * t157 - t492 + t496) * t155 + t321 * t489 - t42 * t131 - t28 * t81 - t27 * t82 + (-t82 * t385 - t81 * t386 + t350) * pkin(11) - pkin(5) * t16 + t333 + (-t17 * t27 - t18 * t28 - t40 * t43 - pkin(5) * t9 + (-t17 * t385 - t18 * t386 + t358) * pkin(11)) * m(7); -t40 * (mrSges(7,1) * t129 + mrSges(7,2) * t128) + (Ifges(7,1) * t128 - t435) * t471 + t46 * t470 + (Ifges(7,5) * t128 - Ifges(7,6) * t129) * t468 - t17 * t81 + t18 * t82 + (t128 * t17 + t129 * t18) * mrSges(7,3) + t360 + t5 + (-Ifges(7,2) * t129 + t124 + t47) * t473;];
tauc  = t29(:);
