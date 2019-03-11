% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:02:12
% EndTime: 2019-03-09 18:02:49
% DurationCPUTime: 21.22s
% Computational Cost: add. (28884->775), mult. (69202->1020), div. (0->0), fcn. (53239->18), ass. (0->377)
t356 = qJ(2) + qJ(3);
t348 = pkin(11) + t356;
t337 = qJ(5) + t348;
t328 = sin(t337);
t329 = cos(t337);
t359 = sin(qJ(6));
t508 = mrSges(7,2) * t359;
t609 = t328 * t508 + t329 * (m(7) * pkin(10) + mrSges(7,3));
t355 = qJD(2) + qJD(3);
t347 = qJD(5) + t355;
t364 = cos(qJ(6));
t361 = sin(qJ(3));
t362 = sin(qJ(2));
t366 = cos(qJ(3));
t367 = cos(qJ(2));
t296 = -t361 * t362 + t366 * t367;
t281 = t296 * qJD(1);
t297 = t361 * t367 + t362 * t366;
t282 = t297 * qJD(1);
t357 = sin(pkin(11));
t358 = cos(pkin(11));
t227 = t281 * t357 + t282 * t358;
t360 = sin(qJ(5));
t365 = cos(qJ(5));
t421 = t358 * t281 - t282 * t357;
t588 = t227 * t365 + t360 * t421;
t152 = t347 * t364 - t359 * t588;
t153 = t347 * t359 + t364 * t588;
t480 = mrSges(6,1) * t347 + mrSges(7,1) * t152 - mrSges(7,2) * t153 - mrSges(6,3) * t588;
t369 = -pkin(8) - pkin(7);
t321 = t369 * t367;
t303 = qJD(1) * t321;
t283 = t361 * t303;
t320 = t369 * t362;
t302 = qJD(1) * t320;
t289 = qJD(2) * pkin(2) + t302;
t232 = t366 * t289 + t283;
t276 = t282 * qJ(4);
t199 = t232 - t276;
t186 = pkin(3) * t355 + t199;
t286 = t366 * t303;
t233 = t289 * t361 - t286;
t479 = qJ(4) * t281;
t200 = t233 + t479;
t189 = t357 * t200;
t126 = t358 * t186 - t189;
t219 = pkin(9) * t227;
t112 = pkin(4) * t355 + t126 - t219;
t468 = t358 * t200;
t127 = t357 * t186 + t468;
t519 = pkin(9) * t421;
t114 = t127 + t519;
t56 = t112 * t365 - t114 * t360;
t54 = -pkin(5) * t347 - t56;
t585 = m(6) * t56 - m(7) * t54 + t480;
t589 = -t227 * t360 + t365 * t421;
t164 = qJD(6) - t589;
t163 = Ifges(6,4) * t589;
t352 = t367 * pkin(2);
t339 = t352 + pkin(1);
t319 = t339 * qJD(1);
t244 = -pkin(3) * t281 + qJD(4) - t319;
t177 = -pkin(4) * t421 + t244;
t411 = mrSges(7,1) * t359 + mrSges(7,2) * t364;
t391 = t54 * t411;
t151 = Ifges(7,4) * t152;
t78 = t153 * Ifges(7,1) + t164 * Ifges(7,5) + t151;
t482 = t364 * t78;
t489 = t347 * Ifges(6,5);
t514 = t56 * mrSges(6,3);
t530 = t359 / 0.2e1;
t498 = Ifges(7,4) * t153;
t77 = Ifges(7,2) * t152 + Ifges(7,6) * t164 + t498;
t597 = t588 * Ifges(6,1);
t98 = t163 + t489 + t597;
t608 = -t177 * mrSges(6,2) - t391 - t482 / 0.2e1 + t77 * t530 + t514 - t98 / 0.2e1 - t489 / 0.2e1 - t163 / 0.2e1;
t107 = pkin(5) * t588 - pkin(10) * t589;
t488 = t347 * Ifges(6,6);
t499 = Ifges(6,4) * t588;
t57 = t112 * t360 + t114 * t365;
t513 = t57 * mrSges(6,3);
t55 = pkin(10) * t347 + t57;
t85 = -pkin(5) * t589 - pkin(10) * t588 + t177;
t23 = t359 * t85 + t364 * t55;
t577 = t23 * mrSges(7,2);
t22 = -t359 * t55 + t364 * t85;
t578 = t22 * mrSges(7,1);
t492 = t164 * Ifges(7,3);
t493 = t153 * Ifges(7,5);
t494 = t152 * Ifges(7,6);
t76 = t492 + t493 + t494;
t598 = t589 * Ifges(6,2);
t97 = t488 + t499 + t598;
t607 = t577 - t177 * mrSges(6,1) + t513 + t97 / 0.2e1 - t76 / 0.2e1 + t499 / 0.2e1 + t488 / 0.2e1 - t578;
t605 = t355 / 0.2e1;
t604 = t421 / 0.2e1;
t562 = t329 * pkin(5) + t328 * pkin(10);
t600 = m(7) * t562;
t599 = Ifges(5,4) * t227;
t596 = -t329 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t328;
t333 = sin(t348);
t334 = cos(t348);
t349 = sin(t356);
t350 = cos(t356);
t595 = -t350 * mrSges(4,1) - t334 * mrSges(5,1) + mrSges(4,2) * t349 + t333 * mrSges(5,2);
t451 = qJD(1) * qJD(2);
t306 = qJDD(1) * t367 - t362 * t451;
t511 = mrSges(7,1) * t364;
t594 = t508 - t511;
t363 = sin(qJ(1));
t368 = cos(qJ(1));
t593 = g(1) * t368 + g(2) * t363;
t509 = mrSges(6,2) * t329;
t592 = mrSges(4,1) * t349 + mrSges(5,1) * t333 + mrSges(6,1) * t328 + mrSges(4,2) * t350 + mrSges(5,2) * t334 + t509;
t307 = qJDD(1) * t362 + t367 * t451;
t385 = t296 * qJD(3);
t206 = qJD(1) * t385 + t306 * t361 + t307 * t366;
t386 = t297 * qJD(3);
t207 = -qJD(1) * t386 + t306 * t366 - t307 * t361;
t146 = t206 * t358 + t207 * t357;
t353 = qJDD(2) + qJDD(3);
t295 = t307 * pkin(7);
t243 = qJDD(2) * pkin(2) - pkin(8) * t307 - t295;
t294 = t306 * pkin(7);
t245 = pkin(8) * t306 + t294;
t150 = -qJD(3) * t233 + t366 * t243 - t245 * t361;
t108 = pkin(3) * t353 - qJ(4) * t206 - qJD(4) * t282 + t150;
t455 = qJD(3) * t366;
t456 = qJD(3) * t361;
t149 = t361 * t243 + t366 * t245 + t289 * t455 + t303 * t456;
t110 = qJ(4) * t207 + qJD(4) * t281 + t149;
t51 = t358 * t108 - t110 * t357;
t34 = pkin(4) * t353 - pkin(9) * t146 + t51;
t145 = -t206 * t357 + t207 * t358;
t52 = t357 * t108 + t358 * t110;
t36 = pkin(9) * t145 + t52;
t11 = -qJD(5) * t57 + t34 * t365 - t36 * t360;
t346 = qJDD(5) + t353;
t68 = qJD(5) * t589 + t145 * t360 + t146 * t365;
t59 = mrSges(6,1) * t346 - mrSges(6,3) * t68;
t591 = m(6) * t11 + t59;
t582 = -m(7) - m(6);
t581 = t599 / 0.2e1 + Ifges(5,2) * t604 + Ifges(5,6) * t605;
t580 = t306 / 0.2e1;
t529 = t367 / 0.2e1;
t576 = Ifges(5,4) * t421;
t571 = t367 * Ifges(3,2);
t238 = -t302 * t361 + t286;
t204 = t238 - t479;
t239 = t366 * t302 + t283;
t205 = -t276 + t239;
t143 = t358 * t204 - t205 * t357;
t467 = t358 * t361;
t495 = pkin(2) * qJD(3);
t270 = (-t357 * t366 - t467) * t495;
t566 = -t143 + t270;
t144 = t357 * t204 + t358 * t205;
t469 = t357 * t361;
t271 = (t358 * t366 - t469) * t495;
t565 = -t144 + t271;
t154 = -mrSges(6,2) * t347 + mrSges(6,3) * t589;
t94 = -mrSges(7,2) * t164 + mrSges(7,3) * t152;
t95 = mrSges(7,1) * t164 - mrSges(7,3) * t153;
t402 = -t359 * t95 + t364 * t94;
t564 = -t154 - t402;
t527 = pkin(2) * t366;
t338 = pkin(3) + t527;
t272 = -pkin(2) * t469 + t358 * t338;
t266 = pkin(4) + t272;
t274 = pkin(2) * t467 + t338 * t357;
t216 = t360 * t266 + t365 * t274;
t247 = t361 * t320 - t366 * t321;
t523 = pkin(3) * t358;
t335 = pkin(4) + t523;
t524 = pkin(3) * t357;
t275 = t360 * t335 + t365 * t524;
t561 = t329 * t594 + t596;
t234 = t296 * t358 - t297 * t357;
t235 = t296 * t357 + t297 * t358;
t172 = t234 * t360 + t235 * t365;
t453 = qJD(6) * t364;
t240 = qJD(2) * t296 + t385;
t241 = -qJD(2) * t297 - t386;
t175 = -t240 * t357 + t241 * t358;
t176 = t240 * t358 + t241 * t357;
t399 = t365 * t234 - t235 * t360;
t88 = qJD(5) * t399 + t175 * t360 + t176 * t365;
t393 = t172 * t453 + t359 * t88;
t454 = qJD(6) * t359;
t560 = -t22 * t453 - t23 * t454;
t458 = qJD(1) * t367;
t459 = qJD(1) * t362;
t520 = pkin(7) * t367;
t521 = pkin(7) * t362;
t559 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t459) * t520 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t458) * t521;
t558 = t294 * t367 + t295 * t362;
t69 = -qJD(5) * t588 + t145 * t365 - t146 * t360;
t556 = t561 + t595;
t41 = qJD(6) * t152 + t346 * t359 + t364 * t68;
t67 = qJDD(6) - t69;
t20 = mrSges(7,1) * t67 - mrSges(7,3) * t41;
t42 = -qJD(6) * t153 + t346 * t364 - t359 * t68;
t21 = -mrSges(7,2) * t67 + mrSges(7,3) * t42;
t554 = -t359 * t20 + t364 * t21 - t95 * t453 - t94 * t454;
t478 = qJDD(1) * pkin(1);
t277 = -pkin(2) * t306 - t478;
t174 = -pkin(3) * t207 + qJDD(4) + t277;
t104 = -pkin(4) * t145 + t174;
t19 = -pkin(5) * t69 - pkin(10) * t68 + t104;
t10 = qJD(5) * t56 + t360 * t34 + t365 * t36;
t7 = pkin(10) * t346 + t10;
t2 = qJD(6) * t22 + t19 * t359 + t364 * t7;
t3 = -qJD(6) * t23 + t19 * t364 - t359 * t7;
t553 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t354 = -qJ(4) + t369;
t552 = -m(3) * pkin(7) + m(4) * t369 + m(5) * t354 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t316 = -mrSges(3,1) * t367 + mrSges(3,2) * t362;
t336 = pkin(3) * t350;
t460 = t336 + t352;
t551 = mrSges(2,1) + m(5) * (pkin(1) + t460) + m(4) * t339 + m(3) * pkin(1) - t316 - t595 - t596;
t405 = t22 * t364 + t23 * t359;
t515 = t3 * t359;
t374 = -qJD(6) * t405 - t515;
t516 = t2 * t364;
t550 = m(7) * (t374 + t516) + t554;
t404 = -t22 * t359 + t23 * t364;
t549 = m(6) * t57 + m(7) * t404 - t564;
t547 = m(7) * pkin(5);
t546 = t41 / 0.2e1;
t545 = t42 / 0.2e1;
t544 = t67 / 0.2e1;
t539 = -t152 / 0.2e1;
t538 = -t153 / 0.2e1;
t537 = t153 / 0.2e1;
t536 = -t164 / 0.2e1;
t534 = t282 / 0.2e1;
t528 = pkin(2) * t362;
t526 = pkin(3) * t282;
t525 = pkin(3) * t349;
t522 = pkin(5) * t328;
t507 = mrSges(4,3) * t281;
t506 = mrSges(4,3) * t282;
t505 = mrSges(5,3) * t126;
t504 = mrSges(5,3) * t127;
t503 = mrSges(7,3) * t364;
t502 = Ifges(3,4) * t362;
t501 = Ifges(3,4) * t367;
t500 = Ifges(4,4) * t282;
t497 = Ifges(7,4) * t359;
t496 = Ifges(7,4) * t364;
t487 = t359 * mrSges(7,3);
t475 = t172 * t359;
t474 = t172 * t364;
t466 = t359 * t363;
t465 = t359 * t368;
t464 = t363 * t364;
t463 = t364 * t368;
t435 = qJD(2) * t369;
t304 = t362 * t435;
t305 = t367 * t435;
t181 = t366 * t304 + t361 * t305 + t320 * t455 + t321 * t456;
t141 = qJ(4) * t241 + qJD(4) * t296 + t181;
t182 = -qJD(3) * t247 - t304 * t361 + t366 * t305;
t142 = -qJ(4) * t240 - qJD(4) * t297 + t182;
t91 = t358 * t141 + t357 * t142;
t133 = t358 * t199 - t189;
t246 = t366 * t320 + t321 * t361;
t222 = -qJ(4) * t297 + t246;
t223 = qJ(4) * t296 + t247;
t157 = t357 * t222 + t358 * t223;
t327 = pkin(4) * t334;
t461 = t327 + t336;
t457 = qJD(2) * t362;
t449 = Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t67;
t343 = pkin(2) * t457;
t445 = t328 * t511;
t342 = pkin(2) * t459;
t437 = t482 / 0.2e1;
t436 = t327 + t460;
t431 = -t69 * mrSges(6,1) + t68 * mrSges(6,2);
t428 = -t454 / 0.2e1;
t228 = -pkin(3) * t241 + t343;
t427 = t451 / 0.2e1;
t426 = -t145 * mrSges(5,1) + t146 * mrSges(5,2);
t90 = -t141 * t357 + t358 * t142;
t132 = -t199 * t357 - t468;
t156 = t358 * t222 - t223 * t357;
t420 = t461 + t562;
t257 = -pkin(3) * t296 - t339;
t419 = t609 * t363;
t418 = t609 * t368;
t183 = pkin(4) * t227 + t526;
t290 = -pkin(4) * t333 - t525;
t416 = mrSges(3,1) * t362 + mrSges(3,2) * t367;
t410 = Ifges(7,1) * t364 - t497;
t409 = t502 + t571;
t408 = -Ifges(7,2) * t359 + t496;
t407 = Ifges(3,5) * t367 - Ifges(3,6) * t362;
t406 = Ifges(7,5) * t364 - Ifges(7,6) * t359;
t122 = -pkin(9) * t235 + t156;
t123 = pkin(9) * t234 + t157;
t84 = t122 * t360 + t123 * t365;
t193 = -pkin(4) * t234 + t257;
t92 = -pkin(5) * t399 - pkin(10) * t172 + t193;
t32 = t359 * t92 + t364 * t84;
t31 = -t359 * t84 + t364 * t92;
t401 = t365 * t122 - t123 * t360;
t215 = t266 * t365 - t274 * t360;
t137 = -pkin(4) * t175 + t228;
t273 = t335 * t365 - t360 * t524;
t397 = -pkin(9) * t176 + t90;
t396 = t132 - t519;
t395 = t143 - t519;
t394 = pkin(1) * t416;
t392 = t172 * t454 - t364 * t88;
t390 = t152 * t408;
t389 = t153 * t410;
t388 = t164 * t406;
t387 = t362 * (Ifges(3,1) * t367 - t502);
t265 = t290 - t528;
t376 = m(7) * (t265 - t522) - t445;
t375 = m(7) * (t290 - t522) - t445;
t87 = t107 + t183;
t373 = t509 + (mrSges(6,1) + t511 + t547) * t328;
t14 = Ifges(7,4) * t41 + Ifges(7,2) * t42 + Ifges(7,6) * t67;
t15 = t41 * Ifges(7,1) + t42 * Ifges(7,4) + t67 * Ifges(7,5);
t8 = -pkin(5) * t346 - t11;
t371 = -t10 * mrSges(6,2) + t2 * t503 + t15 * t530 + t364 * t14 / 0.2e1 + Ifges(6,3) * t346 + (Ifges(7,1) * t359 + t496) * t546 + (Ifges(7,2) * t364 + t497) * t545 + (Ifges(7,5) * t359 + Ifges(7,6) * t364) * t544 + t8 * t594 + Ifges(6,6) * t69 + Ifges(6,5) * t68 + t77 * t428 + t11 * mrSges(6,1) + (t391 + t437) * qJD(6) + (t390 + t389 + t388) * qJD(6) / 0.2e1;
t162 = Ifges(5,1) * t227 + Ifges(5,5) * t355 + t576;
t220 = Ifges(4,2) * t281 + t355 * Ifges(4,6) + t500;
t278 = Ifges(4,4) * t281;
t221 = Ifges(4,1) * t282 + t355 * Ifges(4,5) + t278;
t370 = (Ifges(7,3) * t536 + Ifges(7,5) * t538 + Ifges(7,6) * t539 + t598 / 0.2e1 + t607) * t588 - (-Ifges(4,2) * t282 + t221 + t278) * t281 / 0.2e1 - (Ifges(4,5) * t281 + Ifges(5,5) * t421 - Ifges(4,6) * t282 - Ifges(5,6) * t227) * t355 / 0.2e1 - (-Ifges(5,2) * t227 + t162 + t576) * t421 / 0.2e1 - t244 * (mrSges(5,1) * t227 + mrSges(5,2) * t421) + t227 * t581 + t227 * t504 + t421 * t505 + t233 * t506 + t232 * t507 + (Ifges(4,3) + Ifges(5,3)) * t353 + t319 * (mrSges(4,1) * t282 + mrSges(4,2) * t281) + t220 * t534 - t227 * (Ifges(5,1) * t421 - t599) / 0.2e1 + (t23 * t487 + t22 * t503 + t406 * t536 + t410 * t538 + t408 * t539 - t597 / 0.2e1 + t608) * t589 + Ifges(4,5) * t206 + Ifges(4,6) * t207 + t560 * mrSges(7,3) - t149 * mrSges(4,2) + t150 * mrSges(4,1) + t371 + Ifges(5,6) * t145 + Ifges(5,5) * t146 - t282 * (Ifges(4,1) * t281 - t500) / 0.2e1 - t3 * t487 + t51 * mrSges(5,1) - t52 * mrSges(5,2);
t351 = -pkin(9) + t354;
t341 = Ifges(3,4) * t458;
t280 = Ifges(3,1) * t459 + Ifges(3,5) * qJD(2) + t341;
t279 = Ifges(3,6) * qJD(2) + qJD(1) * t409;
t267 = -pkin(5) - t273;
t262 = pkin(1) + t436;
t261 = t329 * t463 + t466;
t260 = -t329 * t465 + t464;
t259 = -t329 * t464 + t465;
t258 = t329 * t466 + t463;
t250 = mrSges(4,1) * t355 - t506;
t249 = -mrSges(4,2) * t355 + t507;
t248 = t342 + t526;
t231 = -mrSges(4,1) * t281 + mrSges(4,2) * t282;
t210 = -pkin(5) - t215;
t209 = mrSges(5,1) * t355 - mrSges(5,3) * t227;
t208 = -mrSges(5,2) * t355 + mrSges(5,3) * t421;
t188 = -mrSges(4,2) * t353 + mrSges(4,3) * t207;
t187 = mrSges(4,1) * t353 - mrSges(4,3) * t206;
t178 = t183 + t342;
t170 = -mrSges(5,1) * t421 + mrSges(5,2) * t227;
t129 = mrSges(5,1) * t353 - mrSges(5,3) * t146;
t128 = -mrSges(5,2) * t353 + mrSges(5,3) * t145;
t118 = -t219 + t144;
t116 = -t219 + t133;
t106 = -mrSges(6,1) * t589 + mrSges(6,2) * t588;
t89 = qJD(5) * t172 - t365 * t175 + t176 * t360;
t86 = t342 + t87;
t75 = pkin(9) * t175 + t91;
t71 = t365 * t118 + t360 * t395;
t62 = t365 * t116 + t360 * t396;
t60 = -mrSges(6,2) * t346 + mrSges(6,3) * t69;
t30 = t107 * t359 + t364 * t56;
t29 = t107 * t364 - t359 * t56;
t28 = t359 * t86 + t364 * t71;
t27 = -t359 * t71 + t364 * t86;
t26 = t359 * t87 + t364 * t62;
t25 = -t359 * t62 + t364 * t87;
t24 = pkin(5) * t89 - pkin(10) * t88 + t137;
t18 = -mrSges(7,1) * t42 + mrSges(7,2) * t41;
t16 = qJD(5) * t401 + t360 * t397 + t365 * t75;
t5 = -qJD(6) * t32 - t16 * t359 + t24 * t364;
t4 = qJD(6) * t31 + t16 * t364 + t24 * t359;
t1 = [(t367 * t501 + t387) * t427 - (Ifges(7,3) * t544 + Ifges(7,6) * t545 + Ifges(7,5) * t546 + t104 * mrSges(6,1) - Ifges(6,2) * t69 - Ifges(6,4) * t68 + t449 / 0.2e1 - Ifges(6,6) * t346 - t10 * mrSges(6,3) + t553) * t399 + m(5) * (t126 * t90 + t127 * t91 + t156 * t51 + t157 * t52 + t174 * t257 + t228 * t244) + m(7) * (t2 * t32 + t22 * t5 + t23 * t4 + t3 * t31) + m(6) * (t10 * t84 + t104 * t193 + t137 * t177 + t16 * t57) + (mrSges(4,1) * t339 + Ifges(4,4) * t297 + Ifges(4,2) * t296) * t207 + t175 * t504 - (m(7) * t8 + t18 - t591) * t401 + (t104 * mrSges(6,2) - t11 * mrSges(6,3) + Ifges(6,1) * t68 + Ifges(6,4) * t69 + Ifges(6,5) * t346 + t406 * t544 + t408 * t545 + t410 * t546 + t411 * t8 + t428 * t78) * t172 + t54 * (mrSges(7,1) * t393 - mrSges(7,2) * t392) + t164 * (-Ifges(7,5) * t392 - Ifges(7,6) * t393 + Ifges(7,3) * t89) / 0.2e1 + t152 * (-Ifges(7,4) * t392 - Ifges(7,2) * t393 + Ifges(7,6) * t89) / 0.2e1 + (t149 * t296 - t150 * t297 - t232 * t240 + t233 * t241) * mrSges(4,3) + (t234 * t52 - t235 * t51) * mrSges(5,3) + (-mrSges(4,2) * t339 + Ifges(4,1) * t297 + Ifges(4,4) * t296) * t206 + (Ifges(4,5) * t297 + Ifges(5,5) * t235 + Ifges(4,6) * t296 + Ifges(5,6) * t234) * t353 + (-Ifges(7,1) * t392 - Ifges(7,4) * t393 + Ifges(7,5) * t89) * t537 + t89 * t578 + t409 * t580 + t175 * t581 + (Ifges(4,1) * t240 + Ifges(4,4) * t241) * t534 + (t306 * t520 + t307 * t521 + t558) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t558) + t16 * t154 + t156 * t129 + t157 * t128 + (-t259 * mrSges(7,1) - t258 * mrSges(7,2) + (-t351 * t582 + t552) * t368 + (m(6) * t262 - m(7) * (-t262 - t562) + t551) * t363) * g(1) + m(4) * (t149 * t247 + t150 * t246 + t181 * t233 + t182 * t232 - t277 * t339 - t319 * t343) - t585 * (qJD(5) * t84 + t360 * t75 - t365 * t397) + t176 * t162 / 0.2e1 + t177 * (mrSges(6,1) * t89 + mrSges(6,2) * t88) + (t280 * t529 + t407 * qJD(2) / 0.2e1 - t559) * qJD(2) - t393 * t77 / 0.2e1 + t307 * t501 / 0.2e1 + (-t261 * mrSges(7,1) - t260 * mrSges(7,2) + t582 * (t368 * t262 - t351 * t363) + t552 * t363 + (-t551 - t600) * t368) * g(2) + (-t2 * t475 + t22 * t392 - t23 * t393 - t3 * t474) * mrSges(7,3) + t88 * t437 + t137 * t106 + t5 * t95 - t89 * t97 / 0.2e1 + t88 * t98 / 0.2e1 + t4 * t94 + t89 * t76 / 0.2e1 + Ifges(2,3) * qJDD(1) + t84 * t60 - t89 * t513 - t88 * t514 - t176 * t505 - t316 * t478 + t15 * t474 / 0.2e1 - t14 * t475 / 0.2e1 + t31 * t20 + t32 * t21 - t279 * t457 / 0.2e1 - t394 * t451 + t193 * t431 + t257 * t426 + t91 * t208 + t90 * t209 + t227 * (Ifges(5,1) * t176 + Ifges(5,4) * t175) / 0.2e1 + t228 * t170 + t174 * (-mrSges(5,1) * t234 + mrSges(5,2) * t235) + t146 * (Ifges(5,1) * t235 + Ifges(5,4) * t234) + t145 * (Ifges(5,4) * t235 + Ifges(5,2) * t234) + t240 * t221 / 0.2e1 + t241 * t220 / 0.2e1 + t231 * t343 + t244 * (-mrSges(5,1) * t175 + mrSges(5,2) * t176) + t246 * t187 + t247 * t188 + t181 * t249 + t182 * t250 + t281 * (Ifges(4,4) * t240 + Ifges(4,2) * t241) / 0.2e1 - t89 * t577 + (-mrSges(3,1) * t521 - mrSges(3,2) * t520 + 0.2e1 * Ifges(3,6) * t529) * qJDD(2) + (Ifges(3,1) * t307 + Ifges(3,4) * t580 + Ifges(3,5) * qJDD(2) - t427 * t571) * t362 + t277 * (-mrSges(4,1) * t296 + mrSges(4,2) * t297) + (Ifges(5,4) * t176 + Ifges(5,2) * t175) * t604 + (Ifges(4,5) * t240 + Ifges(5,5) * t176 + Ifges(4,6) * t241 + Ifges(5,6) * t175) * t605 - pkin(1) * (-mrSges(3,1) * t306 + mrSges(3,2) * t307) - t319 * (-mrSges(4,1) * t241 + mrSges(4,2) * t240) + t347 * (Ifges(6,5) * t88 - Ifges(6,6) * t89) / 0.2e1 + t589 * (Ifges(6,4) * t88 - Ifges(6,2) * t89) / 0.2e1 + t588 * (Ifges(6,1) * t88 - Ifges(6,4) * t89) / 0.2e1 + (Ifges(3,4) * t307 + Ifges(3,2) * t306) * t529; t593 * (m(4) * t528 - m(6) * t265 - m(5) * (-t525 - t528) + t416 + t592) + t585 * (-qJD(5) * t216 + (t270 - t395) * t365 + (t118 - t271) * t360) + t370 + (t10 * t216 + t11 * t215 - t177 * t178 - t57 * t71) * m(6) + (t210 * t8 - t22 * t27 - t23 * t28) * m(7) + (m(4) * (t149 * t361 + t150 * t366 + (-t232 * t361 + t233 * t366) * qJD(3)) + t361 * t188 - t250 * t456 + t249 * t455) * pkin(2) - (-Ifges(3,2) * t459 + t280 + t341) * t458 / 0.2e1 + t187 * t527 + (t559 + (t394 - t387 / 0.2e1) * qJD(1)) * qJD(1) + t549 * (qJD(5) * t215 + t270 * t360 + t271 * t365) + t550 * (pkin(10) + t216) - t71 * t154 - m(4) * (t232 * t238 + t233 * t239 - t319 * t342) - t178 * t106 + t565 * t208 + t566 * t209 + (t126 * t566 + t127 * t565 - t244 * t248 + t272 * t51 + t274 * t52) * m(5) - t27 * t95 - t28 * t94 + Ifges(3,3) * qJDD(2) + t279 * t459 / 0.2e1 - t407 * t451 / 0.2e1 - t231 * t342 + t210 * t18 + t215 * t59 + t216 * t60 - g(1) * (t368 * t376 + t418) - g(2) * (t363 * t376 + t419) - t248 * t170 - t239 * t249 - t238 * t250 + t272 * t129 + t274 * t128 - t294 * mrSges(3,2) - t295 * mrSges(3,1) + (-m(4) * t352 - m(5) * t460 - m(6) * t436 + t316 - m(7) * (t352 + t420) + t556) * g(3) + Ifges(3,6) * t306 + Ifges(3,5) * t307; t370 + ((t357 * t52 + t358 * t51) * pkin(3) - t126 * t132 - t127 * t133 - t244 * t526) * m(5) + t129 * t523 + t128 * t524 + t550 * (pkin(10) + t275) - t62 * t154 + (t10 * t275 - t177 * t183 - t57 * t62) * m(6) - t183 * t106 - t25 * t95 - t26 * t94 - t170 * t526 - t133 * t208 - t132 * t209 - g(1) * (t368 * t375 + t418) - g(2) * (t363 * t375 + t419) + (-t22 * t25 - t23 * t26 + t267 * t8) * m(7) - t232 * t249 + t233 * t250 + t267 * t18 + t275 * t60 + (-m(5) * t336 - m(6) * t461 - m(7) * t420 + t556) * g(3) + t585 * (-t275 * qJD(5) + t116 * t360 - t365 * t396) + t593 * (m(5) * t525 - m(6) * t290 + t592) + (qJD(5) * t549 + t591) * t273; t364 * t20 - t421 * t208 + t227 * t209 + t359 * t21 + t480 * t588 + t402 * qJD(6) + t564 * t589 + t426 + t431 + (-g(1) * t363 + g(2) * t368) * (m(5) - t582) + (t164 * t404 + t2 * t359 + t3 * t364 - t588 * t54) * m(7) + (t56 * t588 - t57 * t589 + t104) * m(6) + (t126 * t227 - t127 * t421 + t174) * m(5); (-t390 / 0.2e1 - t389 / 0.2e1 - t388 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t588 + t405 * mrSges(7,3) + t608) * t589 + t374 * mrSges(7,3) + (t368 * t373 - t418) * g(1) + (t363 * t373 - t419) * g(2) + t480 * t57 + (-t494 / 0.2e1 - t493 / 0.2e1 - t492 / 0.2e1 + t607) * t588 - t8 * t547 - m(7) * (t22 * t29 + t23 * t30 + t54 * t57) - t56 * t154 + (t561 - t600) * g(3) + t371 - t29 * t95 - t30 * t94 - pkin(5) * t18 + (m(7) * (-t515 + t516 + t560) + t554) * pkin(10); -t54 * (mrSges(7,1) * t153 + mrSges(7,2) * t152) + (Ifges(7,1) * t152 - t498) * t538 + t77 * t537 + (Ifges(7,5) * t152 - Ifges(7,6) * t153) * t536 - t22 * t94 + t23 * t95 - g(1) * (mrSges(7,1) * t260 - mrSges(7,2) * t261) - g(2) * (-mrSges(7,1) * t258 + mrSges(7,2) * t259) + g(3) * t411 * t328 + (t152 * t22 + t153 * t23) * mrSges(7,3) + t449 + (-Ifges(7,2) * t153 + t151 + t78) * t539 + t553;];
tau  = t1;
