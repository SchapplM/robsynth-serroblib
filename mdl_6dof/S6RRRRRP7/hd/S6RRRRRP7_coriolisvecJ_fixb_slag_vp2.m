% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:42
% EndTime: 2019-03-10 01:35:53
% DurationCPUTime: 37.85s
% Computational Cost: add. (25149->907), mult. (64527->1248), div. (0->0), fcn. (50860->10), ass. (0->379)
t367 = cos(qJ(2));
t363 = sin(qJ(2));
t358 = sin(pkin(6));
t451 = qJD(1) * t358;
t429 = t363 * t451;
t359 = cos(pkin(6));
t450 = qJD(1) * t359;
t440 = pkin(1) * t450;
t310 = -pkin(8) * t429 + t367 * t440;
t385 = (pkin(2) * t363 - pkin(9) * t367) * t358;
t311 = qJD(1) * t385;
t362 = sin(qJ(3));
t366 = cos(qJ(3));
t244 = -t362 * t310 + t366 * t311;
t537 = -pkin(10) - pkin(9);
t430 = qJD(3) * t537;
t596 = -(-pkin(10) * t366 * t367 + pkin(3) * t363) * t451 - t244 + t366 * t430;
t245 = t366 * t310 + t362 * t311;
t428 = t367 * t451;
t414 = t362 * t428;
t595 = -pkin(10) * t414 - t362 * t430 + t245;
t361 = sin(qJ(4));
t365 = cos(qJ(4));
t326 = t361 * t366 + t362 * t365;
t544 = qJD(3) + qJD(4);
t267 = t544 * t326;
t278 = t326 * t428;
t594 = t267 - t278;
t340 = t537 * t362;
t341 = t537 * t366;
t547 = t365 * t340 + t341 * t361;
t572 = qJD(4) * t547 + t361 * t596 - t595 * t365;
t313 = pkin(8) * t428 + t363 * t440;
t270 = pkin(3) * t414 + t313;
t447 = qJD(3) * t362;
t593 = pkin(3) * t447 - t270;
t558 = Ifges(7,4) + Ifges(6,4);
t592 = pkin(11) * t429 - t572;
t325 = t361 * t362 - t365 * t366;
t266 = t544 * t325;
t279 = t325 * t428;
t591 = t593 + (t266 - t279) * pkin(11) + t594 * pkin(4);
t559 = Ifges(7,1) + Ifges(6,1);
t557 = Ifges(7,5) + Ifges(6,5);
t556 = Ifges(7,2) + Ifges(6,2);
t555 = Ifges(7,6) + Ifges(6,6);
t281 = t340 * t361 - t341 * t365;
t571 = -qJD(4) * t281 + t595 * t361 + t365 * t596;
t360 = sin(qJ(5));
t364 = cos(qJ(5));
t247 = t279 * t360 + t364 * t429;
t442 = qJD(5) * t364;
t590 = -t326 * t442 - t247;
t337 = qJD(3) - t428;
t331 = qJD(4) + t337;
t346 = qJD(2) + t450;
t293 = t346 * t366 - t362 * t429;
t294 = t346 * t362 + t366 * t429;
t388 = t293 * t361 + t365 * t294;
t203 = t331 * t364 - t360 * t388;
t416 = t365 * t293 - t294 * t361;
t232 = qJD(5) - t416;
t204 = t331 * t360 + t364 * t388;
t584 = t558 * t204;
t550 = t556 * t203 + t555 * t232 + t584;
t589 = -t550 / 0.2e1;
t588 = t360 * t592 + t364 * t591;
t355 = -pkin(3) * t366 - pkin(2);
t262 = pkin(4) * t325 - pkin(11) * t326 + t355;
t587 = t262 * t442 + t360 * t591 - t364 * t592;
t586 = t558 * t203;
t570 = pkin(4) * t429 - t571;
t549 = t416 * t360;
t585 = qJ(6) * t549 + t364 * qJD(6);
t548 = t204 * t559 + t557 * t232 + t586;
t583 = t548 / 0.2e1;
t248 = -t279 * t364 + t360 * t429;
t271 = t364 * t281;
t386 = qJ(6) * t266 - qJD(6) * t326;
t582 = t386 * t364 + (-t271 + (qJ(6) * t326 - t262) * t360) * qJD(5) + qJ(6) * t248 + t588 + t594 * pkin(5);
t581 = (-qJD(5) * t281 + t386) * t360 + t587 + t590 * qJ(6);
t209 = t360 * t262 + t271;
t580 = -qJD(5) * t209 + t588;
t443 = qJD(5) * t360;
t579 = -t281 * t443 + t587;
t97 = t204 * Ifges(7,5) + t203 * Ifges(7,6) + t232 * Ifges(7,3);
t98 = t204 * Ifges(6,5) + t203 * Ifges(6,6) + t232 * Ifges(6,3);
t578 = t98 + t97;
t352 = pkin(3) * t361 + pkin(11);
t454 = -qJ(6) - t352;
t415 = qJD(5) * t454;
t444 = qJD(4) * t365;
t438 = pkin(3) * t444;
t277 = pkin(9) * t346 + t313;
t306 = (-pkin(2) * t367 - pkin(9) * t363 - pkin(1)) * t358;
t288 = qJD(1) * t306;
t225 = t277 * t366 + t288 * t362;
t192 = pkin(10) * t293 + t225;
t185 = t361 * t192;
t224 = -t277 * t362 + t366 * t288;
t191 = -pkin(10) * t294 + t224;
t113 = t191 * t365 - t185;
t171 = pkin(4) * t388 - pkin(11) * t416;
t146 = pkin(3) * t294 + t171;
t58 = t364 * t113 + t360 * t146;
t577 = t360 * t415 + t364 * t438 - t58 + t585;
t357 = t364 * qJ(6);
t542 = pkin(5) * t388 - t357 * t416;
t57 = -t113 * t360 + t364 * t146;
t576 = (-qJD(6) - t438) * t360 + t364 * t415 - t542 - t57;
t186 = t365 * t192;
t112 = t191 * t361 + t186;
t437 = pkin(5) * t443;
t445 = qJD(4) * t361;
t567 = pkin(5) * t549;
t575 = pkin(3) * t445 - t112 + t437 - t567;
t312 = qJD(2) * t385;
t301 = qJD(1) * t312;
t456 = t358 * t363;
t347 = pkin(8) * t456;
t498 = pkin(1) * t367;
t321 = t359 * t498 - t347;
t314 = t321 * qJD(2);
t302 = qJD(1) * t314;
t165 = -qJD(3) * t225 + t366 * t301 - t302 * t362;
t448 = qJD(2) * t367;
t425 = t366 * t448;
t446 = qJD(3) * t366;
t258 = t346 * t446 + (-t363 * t447 + t425) * t451;
t449 = qJD(2) * t358;
t422 = qJD(1) * t449;
t413 = t363 * t422;
t117 = pkin(3) * t413 - pkin(10) * t258 + t165;
t164 = -t277 * t447 + t288 * t446 + t362 * t301 + t366 * t302;
t426 = t362 * t448;
t259 = -t346 * t447 + (-t363 * t446 - t426) * t451;
t124 = pkin(10) * t259 + t164;
t175 = pkin(3) * t337 + t191;
t30 = t117 * t365 - t361 * t124 - t175 * t445 - t192 * t444;
t27 = -pkin(4) * t413 - t30;
t141 = qJD(4) * t416 + t258 * t365 + t259 * t361;
t87 = qJD(5) * t203 + t141 * t364 + t360 * t413;
t88 = -qJD(5) * t204 - t141 * t360 + t364 * t413;
t33 = -mrSges(6,1) * t88 + mrSges(6,2) * t87;
t574 = -m(6) * t27 - t33;
t573 = t570 + (-t266 * t360 - t590) * pkin(5);
t475 = t388 * Ifges(5,1);
t106 = t175 * t365 - t185;
t231 = Ifges(5,4) * t416;
t276 = -t346 * pkin(2) - t310;
t237 = -t293 * pkin(3) + t276;
t107 = t175 * t361 + t186;
t104 = pkin(11) * t331 + t107;
t118 = -pkin(4) * t416 - pkin(11) * t388 + t237;
t49 = -t104 * t360 + t364 * t118;
t34 = -qJ(6) * t204 + t49;
t31 = pkin(5) * t232 + t34;
t50 = t104 * t364 + t118 * t360;
t35 = qJ(6) * t203 + t50;
t394 = Ifges(7,5) * t364 - Ifges(7,6) * t360;
t376 = t232 * t394;
t396 = Ifges(6,5) * t364 - Ifges(6,6) * t360;
t377 = t232 * t396;
t483 = Ifges(7,4) * t360;
t402 = Ifges(7,1) * t364 - t483;
t378 = t204 * t402;
t486 = Ifges(6,4) * t360;
t404 = Ifges(6,1) * t364 - t486;
t379 = t204 * t404;
t482 = Ifges(7,4) * t364;
t398 = -Ifges(7,2) * t360 + t482;
t380 = t203 * t398;
t485 = Ifges(6,4) * t364;
t400 = -Ifges(6,2) * t360 + t485;
t381 = t203 * t400;
t391 = t360 * t50 + t364 * t49;
t469 = t331 * Ifges(5,5);
t501 = -t364 / 0.2e1;
t502 = t360 / 0.2e1;
t162 = t231 + t469 + t475;
t529 = -t162 / 0.2e1;
t103 = -pkin(4) * t331 - t106;
t405 = mrSges(7,1) * t360 + mrSges(7,2) * t364;
t407 = mrSges(6,1) * t360 + mrSges(6,2) * t364;
t73 = -pkin(5) * t203 + qJD(6) + t103;
t563 = t103 * t407 + t73 * t405;
t546 = t391 * mrSges(6,3) + (t31 * t364 + t35 * t360) * mrSges(7,3) + t106 * mrSges(5,3) + t529 - t381 / 0.2e1 - t380 / 0.2e1 - t379 / 0.2e1 - t378 / 0.2e1 - t377 / 0.2e1 - t376 / 0.2e1 - t237 * mrSges(5,2) - t469 / 0.2e1 + t550 * t502 + t548 * t501 - t563 - t231 / 0.2e1;
t569 = t475 / 0.2e1 - t546;
t477 = t416 * Ifges(5,2);
t434 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t435 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t436 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t468 = t331 * Ifges(5,6);
t488 = Ifges(5,4) * t388;
t161 = t468 + t477 + t488;
t530 = t161 / 0.2e1;
t541 = t435 * t203 + t436 * t204 + t434 * t232 - t107 * mrSges(5,3) - t35 * mrSges(7,2) - t50 * mrSges(6,2) - t530 + t97 / 0.2e1 + t98 / 0.2e1 - t488 / 0.2e1 + t237 * mrSges(5,1) + t31 * mrSges(7,1) - t468 / 0.2e1 + t49 * mrSges(6,1);
t568 = -t541 + t477 / 0.2e1;
t142 = qJD(4) * t388 + t258 * t361 - t365 * t259;
t566 = t142 * t555 + t556 * t88 + t558 * t87;
t565 = t142 * t557 + t558 * t88 + t559 * t87;
t562 = t360 * t557 + t364 * t555;
t561 = t364 * t556 + t483 + t486;
t560 = t360 * t559 + t482 + t485;
t554 = Ifges(7,3) + Ifges(6,3);
t18 = Ifges(7,5) * t87 + Ifges(7,6) * t88 + Ifges(7,3) * t142;
t19 = Ifges(6,5) * t87 + Ifges(6,6) * t88 + Ifges(6,3) * t142;
t553 = t19 + t18;
t492 = -qJ(6) - pkin(11);
t418 = qJD(5) * t492;
t61 = t364 * t106 + t360 * t171;
t552 = t360 * t418 + t585 - t61;
t60 = -t106 * t360 + t364 * t171;
t551 = -qJD(6) * t360 + t364 * t418 - t542 - t60;
t455 = t358 * t367;
t322 = t359 * t363 * pkin(1) + pkin(8) * t455;
t305 = pkin(9) * t359 + t322;
t242 = -t362 * t305 + t366 * t306;
t319 = t359 * t362 + t366 * t456;
t199 = -pkin(3) * t455 - t319 * pkin(10) + t242;
t243 = t366 * t305 + t362 * t306;
t318 = t359 * t366 - t362 * t456;
t215 = pkin(10) * t318 + t243;
t130 = t361 * t199 + t365 * t215;
t121 = -pkin(11) * t455 + t130;
t251 = t318 * t361 + t319 * t365;
t304 = t347 + (-pkin(2) - t498) * t359;
t257 = -t318 * pkin(3) + t304;
t387 = t365 * t318 - t319 * t361;
t159 = -pkin(4) * t387 - t251 * pkin(11) + t257;
t66 = t364 * t121 + t360 * t159;
t545 = -t360 * t49 + t364 * t50;
t472 = t294 * Ifges(4,4);
t221 = t293 * Ifges(4,2) + t337 * Ifges(4,6) + t472;
t289 = Ifges(4,4) * t293;
t222 = t294 * Ifges(4,1) + t337 * Ifges(4,5) + t289;
t389 = t224 * t366 + t225 * t362;
t489 = Ifges(4,4) * t366;
t490 = Ifges(4,4) * t362;
t499 = t366 / 0.2e1;
t504 = t337 / 0.2e1;
t507 = t294 / 0.2e1;
t509 = t293 / 0.2e1;
t543 = -t389 * mrSges(4,3) + t276 * (mrSges(4,1) * t362 + mrSges(4,2) * t366) + (-Ifges(4,2) * t362 + t489) * t509 + (Ifges(4,1) * t366 - t490) * t507 + (Ifges(4,5) * t366 - Ifges(4,6) * t362) * t504 - t362 * t221 / 0.2e1 + t222 * t499;
t540 = Ifges(5,2) / 0.2e1;
t539 = t87 / 0.2e1;
t538 = t88 / 0.2e1;
t535 = pkin(1) * mrSges(3,1);
t534 = pkin(1) * mrSges(3,2);
t531 = t142 / 0.2e1;
t527 = -t203 / 0.2e1;
t525 = -t204 / 0.2e1;
t524 = t204 / 0.2e1;
t520 = -t232 / 0.2e1;
t516 = t387 / 0.2e1;
t514 = t251 / 0.2e1;
t513 = t258 / 0.2e1;
t512 = t259 / 0.2e1;
t508 = -t294 / 0.2e1;
t506 = t318 / 0.2e1;
t505 = t319 / 0.2e1;
t500 = t364 / 0.2e1;
t497 = pkin(3) * t365;
t496 = pkin(11) * t364;
t29 = t361 * t117 + t365 * t124 + t175 * t444 - t192 * t445;
t26 = pkin(11) * t413 + t29;
t315 = t322 * qJD(2);
t303 = qJD(1) * t315;
t223 = -t259 * pkin(3) + t303;
t55 = t142 * pkin(4) - t141 * pkin(11) + t223;
t7 = -t104 * t443 + t118 * t442 + t364 * t26 + t360 * t55;
t495 = t364 * t7;
t491 = Ifges(3,4) * t363;
t481 = Ifges(3,5) * t367;
t480 = t141 * Ifges(5,1);
t479 = t141 * Ifges(5,4);
t478 = t142 * Ifges(5,4);
t476 = t416 * Ifges(5,6);
t474 = t388 * Ifges(5,5);
t473 = t293 * Ifges(4,6);
t471 = t294 * Ifges(4,5);
t470 = t302 * mrSges(3,2);
t467 = t331 * Ifges(5,3);
t466 = t337 * Ifges(4,3);
t465 = t346 * Ifges(3,5);
t459 = t326 * t360;
t457 = t352 * t364;
t132 = -mrSges(6,1) * t203 + mrSges(6,2) * t204;
t214 = mrSges(5,1) * t331 - mrSges(5,3) * t388;
t453 = t214 - t132;
t452 = -mrSges(3,1) * t346 - mrSges(4,1) * t293 + mrSges(4,2) * t294 + mrSges(3,3) * t429;
t433 = Ifges(5,5) * t141 - Ifges(5,6) * t142 + Ifges(5,3) * t413;
t431 = Ifges(4,5) * t258 + Ifges(4,6) * t259 + Ifges(4,3) * t413;
t354 = -pkin(5) * t364 - pkin(4);
t427 = t363 * t449;
t32 = -t88 * mrSges(7,1) + t87 * mrSges(7,2);
t65 = -t121 * t360 + t364 * t159;
t129 = t199 * t365 - t361 * t215;
t208 = t364 * t262 - t281 * t360;
t8 = -qJD(5) * t50 - t26 * t360 + t364 * t55;
t1 = pkin(5) * t142 - qJ(6) * t87 - qJD(6) * t204 + t8;
t411 = -mrSges(6,3) * t8 - mrSges(7,3) * t1;
t410 = -t49 * mrSges(6,3) - t31 * mrSges(7,3);
t409 = -t50 * mrSges(6,3) - t35 * mrSges(7,3);
t120 = pkin(4) * t455 - t129;
t408 = mrSges(6,1) * t364 - mrSges(6,2) * t360;
t406 = mrSges(7,1) * t364 - mrSges(7,2) * t360;
t390 = t164 * t366 - t165 * t362;
t177 = -qJD(3) * t243 + t366 * t312 - t314 * t362;
t269 = qJD(3) * t318 + t358 * t425;
t145 = pkin(3) * t427 - pkin(10) * t269 + t177;
t176 = -t305 * t447 + t306 * t446 + t362 * t312 + t366 * t314;
t268 = -qJD(3) * t319 - t358 * t426;
t156 = pkin(10) * t268 + t176;
t41 = t145 * t365 - t361 * t156 - t199 * t445 - t215 * t444;
t229 = -t360 * t251 - t364 * t455;
t384 = -t364 * t251 + t360 * t455;
t40 = t361 * t145 + t365 * t156 + t199 * t444 - t215 * t445;
t37 = pkin(11) * t427 + t40;
t168 = qJD(4) * t387 + t268 * t361 + t269 * t365;
t169 = qJD(4) * t251 - t365 * t268 + t269 * t361;
t240 = -t268 * pkin(3) + t315;
t72 = t169 * pkin(4) - t168 * pkin(11) + t240;
t9 = -t121 * t443 + t159 * t442 + t360 * t72 + t364 * t37;
t3 = qJ(6) * t88 + qJD(6) * t203 + t7;
t375 = t8 * mrSges(6,1) + t1 * mrSges(7,1) - t7 * mrSges(6,2) - t3 * mrSges(7,2);
t38 = -pkin(4) * t427 - t41;
t10 = -qJD(5) * t66 - t360 * t37 + t364 * t72;
t374 = m(6) * (-qJD(5) * t391 - t8 * t360 + t495);
t12 = -pkin(5) * t88 + t27;
t372 = t3 * t364 * mrSges(7,3) + t30 * mrSges(5,1) - t29 * mrSges(5,2) + mrSges(6,3) * t495 - t12 * t406 - t27 * t408 + t433 + t560 * t539 + t561 * t538 + t562 * t531 + t565 * t502 + t566 * t500 + t443 * t589 + t442 * t583 + t563 * qJD(5) + (t381 + t380 + t379 + t378 + t377 + t376) * qJD(5) / 0.2e1;
t342 = Ifges(3,4) * t428;
t339 = t357 + t496;
t338 = t492 * t360;
t336 = t354 - t497;
t335 = t422 * t481;
t324 = t357 + t457;
t323 = t454 * t360;
t309 = -t346 * mrSges(3,2) + mrSges(3,3) * t428;
t274 = Ifges(3,1) * t429 + t342 + t465;
t273 = Ifges(3,6) * t346 + (Ifges(3,2) * t367 + t491) * t451;
t261 = mrSges(4,1) * t337 - mrSges(4,3) * t294;
t260 = -mrSges(4,2) * t337 + mrSges(4,3) * t293;
t246 = pkin(5) * t459 - t547;
t239 = -mrSges(4,2) * t413 + mrSges(4,3) * t259;
t238 = mrSges(4,1) * t413 - mrSges(4,3) * t258;
t220 = t466 + t471 + t473;
t213 = -mrSges(5,2) * t331 + mrSges(5,3) * t416;
t193 = -mrSges(4,1) * t259 + mrSges(4,2) * t258;
t184 = -qJ(6) * t459 + t209;
t181 = t258 * Ifges(4,1) + t259 * Ifges(4,4) + Ifges(4,5) * t413;
t180 = t258 * Ifges(4,4) + t259 * Ifges(4,2) + Ifges(4,6) * t413;
t172 = pkin(5) * t325 - t326 * t357 + t208;
t170 = -mrSges(5,1) * t416 + mrSges(5,2) * t388;
t160 = t467 + t474 + t476;
t152 = mrSges(6,1) * t232 - mrSges(6,3) * t204;
t151 = mrSges(7,1) * t232 - mrSges(7,3) * t204;
t150 = -mrSges(6,2) * t232 + mrSges(6,3) * t203;
t149 = -mrSges(7,2) * t232 + mrSges(7,3) * t203;
t131 = -mrSges(7,1) * t203 + mrSges(7,2) * t204;
t126 = -mrSges(5,2) * t413 - mrSges(5,3) * t142;
t125 = mrSges(5,1) * t413 - mrSges(5,3) * t141;
t110 = qJD(5) * t384 - t360 * t168 + t364 * t427;
t109 = qJD(5) * t229 + t364 * t168 + t360 * t427;
t92 = -pkin(5) * t229 + t120;
t74 = t107 + t567;
t71 = mrSges(5,1) * t142 + mrSges(5,2) * t141;
t63 = Ifges(5,5) * t413 - t478 + t480;
t62 = -t142 * Ifges(5,2) + Ifges(5,6) * t413 + t479;
t52 = qJ(6) * t229 + t66;
t48 = -mrSges(6,2) * t142 + mrSges(6,3) * t88;
t47 = -mrSges(7,2) * t142 + mrSges(7,3) * t88;
t46 = mrSges(6,1) * t142 - mrSges(6,3) * t87;
t45 = mrSges(7,1) * t142 - mrSges(7,3) * t87;
t43 = -pkin(5) * t387 + qJ(6) * t384 + t65;
t13 = -pkin(5) * t110 + t38;
t5 = qJ(6) * t110 + qJD(6) * t229 + t9;
t4 = pkin(5) * t169 - qJ(6) * t109 + qJD(6) * t384 + t10;
t2 = [m(6) * (t10 * t49 + t103 * t38 + t120 * t27 + t50 * t9 + t65 * t8 + t66 * t7) + m(7) * (t1 * t43 + t12 * t92 + t13 * t73 + t3 * t52 + t31 * t4 + t35 * t5) + m(5) * (t106 * t41 + t107 * t40 + t129 * t30 + t130 * t29 + t223 * t257 + t237 * t240) + m(4) * (t164 * t243 + t165 * t242 + t176 * t225 + t177 * t224 + t276 * t315 + t303 * t304) + m(3) * (t302 * t322 - t303 * t321 - t310 * t315 + t313 * t314) + (t367 * t274 + t346 * (-Ifges(3,6) * t363 + t481) + (t220 + t160) * t363) * t449 / 0.2e1 + t12 * (-mrSges(7,1) * t229 - mrSges(7,2) * t384) + t27 * (-mrSges(6,1) * t229 - mrSges(6,2) * t384) + t1 * (-mrSges(7,1) * t387 + mrSges(7,3) * t384) + t8 * (-mrSges(6,1) * t387 + mrSges(6,3) * t384) + (t229 * t555 - t384 * t557 - t387 * t554) * t531 + (t229 * t556 - t384 * t558 - t387 * t555) * t538 + (t558 * t229 - t384 * t559 - t387 * t557) * t539 + t29 * (mrSges(5,2) * t455 + mrSges(5,3) * t387) - t142 * (Ifges(5,4) * t251 + Ifges(5,2) * t387 - Ifges(5,6) * t455) / 0.2e1 + t141 * (Ifges(5,1) * t251 + Ifges(5,4) * t387 - Ifges(5,5) * t455) / 0.2e1 + t223 * (-mrSges(5,1) * t387 + mrSges(5,2) * t251) + t3 * (mrSges(7,2) * t387 + mrSges(7,3) * t229) + t7 * (mrSges(6,2) * t387 + mrSges(6,3) * t229) - t553 * t387 / 0.2e1 + t164 * (mrSges(4,2) * t455 + t318 * mrSges(4,3)) + t452 * t315 + ((Ifges(3,5) * t359 / 0.2e1 - t321 * mrSges(3,3) + (-0.2e1 * t534 + 0.3e1 / 0.2e1 * Ifges(3,4) * t367) * t358) * t367 + (Ifges(4,5) * t505 + Ifges(4,6) * t506 + Ifges(5,5) * t514 + Ifges(5,6) * t516 - Ifges(3,6) * t359 - t322 * mrSges(3,3) + (-0.2e1 * t535 - 0.3e1 / 0.2e1 * t491) * t358 + (-Ifges(4,3) / 0.2e1 - Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t455) * t363) * t422 + (-t470 - t303 * mrSges(3,1) + t335 / 0.2e1) * t359 - (t431 + t433) * t455 / 0.2e1 + t416 * (Ifges(5,4) * t168 - Ifges(5,2) * t169 + Ifges(5,6) * t427) / 0.2e1 + t388 * (Ifges(5,1) * t168 - Ifges(5,4) * t169 + Ifges(5,5) * t427) / 0.2e1 + t550 * t110 / 0.2e1 + t106 * (mrSges(5,1) * t427 - mrSges(5,3) * t168) + t224 * (mrSges(4,1) * t427 - mrSges(4,3) * t269) + (t302 * t367 + t303 * t363 + (-t310 * t367 - t313 * t363) * qJD(2)) * t358 * mrSges(3,3) + t30 * (-mrSges(5,1) * t455 - t251 * mrSges(5,3)) + t165 * (-mrSges(4,1) * t455 - t319 * mrSges(4,3)) + t276 * (-mrSges(4,1) * t268 + mrSges(4,2) * t269) + t268 * t221 / 0.2e1 + t269 * t222 / 0.2e1 + t177 * t261 + t257 * t71 + t176 * t260 + t237 * (mrSges(5,1) * t169 + mrSges(5,2) * t168) + t240 * t170 + t242 * t238 + t243 * t239 + t40 * t213 + t41 * t214 + t49 * (mrSges(6,1) * t169 - mrSges(6,3) * t109) + t35 * (-mrSges(7,2) * t169 + mrSges(7,3) * t110) + t50 * (-mrSges(6,2) * t169 + mrSges(6,3) * t110) - t169 * t161 / 0.2e1 + t168 * t162 / 0.2e1 + t31 * (mrSges(7,1) * t169 - mrSges(7,3) * t109) - t273 * t427 / 0.2e1 + t107 * (-mrSges(5,2) * t427 - mrSges(5,3) * t169) + t225 * (-mrSges(4,2) * t427 + mrSges(4,3) * t268) + t331 * (Ifges(5,5) * t168 - Ifges(5,6) * t169 + Ifges(5,3) * t427) / 0.2e1 + t304 * t193 + (t109 * t557 + t110 * t555 + t169 * t554) * t232 / 0.2e1 + (t109 * t558 + t110 * t556 + t169 * t555) * t203 / 0.2e1 + t314 * t309 + (t109 * t559 + t558 * t110 + t557 * t169) * t524 + t303 * (-mrSges(4,1) * t318 + mrSges(4,2) * t319) + (Ifges(4,5) * t269 + Ifges(4,6) * t268 + Ifges(4,3) * t427) * t504 + t181 * t505 + t180 * t506 + (Ifges(4,1) * t269 + Ifges(4,4) * t268 + Ifges(4,5) * t427) * t507 + (Ifges(4,4) * t269 + Ifges(4,2) * t268 + Ifges(4,6) * t427) * t509 + (Ifges(4,4) * t319 + Ifges(4,2) * t318 - Ifges(4,6) * t455) * t512 + t578 * t169 / 0.2e1 + t43 * t45 + t52 * t47 + (Ifges(4,1) * t319 + Ifges(4,4) * t318 - Ifges(4,5) * t455) * t513 + t63 * t514 + t62 * t516 + t109 * t583 + t65 * t46 + t66 * t48 + t92 * t32 + t73 * (-mrSges(7,1) * t110 + mrSges(7,2) * t109) + t103 * (-mrSges(6,1) * t110 + mrSges(6,2) * t109) + t120 * t33 + t129 * t125 + t130 * t126 + t13 * t131 + t38 * t132 - t565 * t384 / 0.2e1 + t566 * t229 / 0.2e1 + t5 * t149 + t9 * t150 + t4 * t151 + t10 * t152; (-t29 * mrSges(5,3) - t479 / 0.2e1 + t18 / 0.2e1 + t19 / 0.2e1 - t62 / 0.2e1 + t223 * mrSges(5,1) + t435 * t88 + t436 * t87 + (t540 + t434) * t142 + t375) * t325 + (-t238 * t362 + t239 * t366) * pkin(9) + t390 * mrSges(4,3) - t470 + (-t224 * t244 - t225 * t245 - t276 * t313 - pkin(2) * t303 + (-qJD(3) * t389 + t390) * pkin(9)) * m(4) - t452 * t313 + (-t106 * t279 + t107 * t278) * mrSges(5,3) - t416 * (-Ifges(5,4) * t279 - Ifges(5,2) * t278) / 0.2e1 - t388 * (-Ifges(5,1) * t279 - Ifges(5,4) * t278) / 0.2e1 - t331 * (-Ifges(5,5) * t279 - Ifges(5,6) * t278) / 0.2e1 - t237 * (mrSges(5,1) * t278 - mrSges(5,2) * t279) + (t362 * pkin(3) * t170 + (-t260 * t362 - t261 * t366) * pkin(9) + t543) * qJD(3) + ((-t465 / 0.2e1 - t274 / 0.2e1 - t342 / 0.2e1 + t310 * mrSges(3,3) + t451 * t534 - t543) * t367 + (-t473 / 0.2e1 - t471 / 0.2e1 - t466 / 0.2e1 - t224 * mrSges(4,1) + t225 * mrSges(4,2) + t273 / 0.2e1 - t220 / 0.2e1 - t160 / 0.2e1 - t474 / 0.2e1 - t476 / 0.2e1 + t107 * mrSges(5,2) - t106 * mrSges(5,1) - t467 / 0.2e1 + t313 * mrSges(3,3) + (t535 + t491 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t367) * t451 + (t346 / 0.2e1 - qJD(2)) * Ifges(3,6) + (Ifges(4,5) * t362 + Ifges(5,5) * t326 + Ifges(4,6) * t366 - Ifges(5,6) * t325) * qJD(2) / 0.2e1) * t363) * t451 + t281 * t126 + (t571 * t106 + t572 * t107 + t223 * t355 + t237 * t593 + t281 * t29 + t547 * t30) * m(5) + t335 + t570 * t132 + t571 * t214 - (t33 - t125) * t547 - t568 * t267 - t569 * t266 + t362 * t181 / 0.2e1 + (-mrSges(4,1) * t366 + mrSges(4,2) * t362 - mrSges(3,1)) * t303 - t31 * (mrSges(7,1) * t278 - mrSges(7,3) * t248) - t49 * (mrSges(6,1) * t278 - mrSges(6,3) * t248) - t35 * (-mrSges(7,2) * t278 + mrSges(7,3) * t247) - t50 * (-mrSges(6,2) * t278 + mrSges(6,3) * t247) - t270 * t170 - t244 * t261 - t245 * t260 + t246 * t32 - t73 * (-mrSges(7,1) * t247 + mrSges(7,2) * t248) - t103 * (-mrSges(6,1) * t247 + mrSges(6,2) * t248) + t208 * t46 + t209 * t48 - pkin(2) * t193 + t184 * t47 + t172 * t45 + t572 * t213 + t573 * t131 - t548 * t248 / 0.2e1 - t310 * t309 + t355 * t71 + t247 * t589 + t180 * t499 + (Ifges(4,2) * t366 + t490) * t512 - t578 * t278 / 0.2e1 + t579 * t150 + t580 * t152 + (t103 * t570 + t208 * t8 + t209 * t7 - t27 * t547 + t49 * t580 + t50 * t579) * m(6) + t581 * t149 + t582 * t151 + (t1 * t172 + t12 * t246 + t184 * t3 + t31 * t582 + t35 * t581 + t573 * t73) * m(7) + (Ifges(4,1) * t362 + t489) * t513 + (Ifges(7,5) * t248 + Ifges(7,6) * t247 + Ifges(7,3) * t278) * t520 + (Ifges(6,5) * t248 + Ifges(6,6) * t247 + Ifges(6,3) * t278) * t520 + (Ifges(7,1) * t248 + Ifges(7,4) * t247 + Ifges(7,5) * t278) * t525 + (Ifges(6,1) * t248 + Ifges(6,4) * t247 + Ifges(6,5) * t278) * t525 + (Ifges(7,4) * t248 + Ifges(7,2) * t247 + Ifges(7,6) * t278) * t527 + (Ifges(6,4) * t248 + Ifges(6,2) * t247 + Ifges(6,6) * t278) * t527 - t279 * t529 + t278 * t530 + (t480 / 0.2e1 - t478 / 0.2e1 - t30 * mrSges(5,3) + t223 * mrSges(5,2) + t12 * t405 + t27 * t407 + t63 / 0.2e1 + (-t1 * t364 - t3 * t360) * mrSges(7,3) + (-t360 * t7 - t364 * t8) * mrSges(6,3) + (t73 * t406 + t103 * t408 + (t31 * t360 - t35 * t364) * mrSges(7,3) - t545 * mrSges(6,3) + t561 * t527 + t560 * t525 + t562 * t520 + t550 * t501) * qJD(5) + (t402 + t404) * t539 + (t400 + t398) * t538 + (t394 + t396) * t531 - (qJD(5) * t548 + t566) * t360 / 0.2e1 + t565 * t500) * t326; (-t352 * t46 + t411) * t360 + ((-t152 * t352 + t410) * t364 + (-t150 * t352 + t409) * t360) * qJD(5) + t431 - m(5) * (-t106 * t112 + t107 * t113) + t352 * t374 - t276 * (mrSges(4,1) * t294 + mrSges(4,2) * t293) + t372 - (-Ifges(4,2) * t294 + t222 + t289) * t293 / 0.2e1 + t453 * t112 + t568 * t388 - t569 * t416 + (t224 * t293 + t225 * t294) * mrSges(4,3) + t225 * t261 - t224 * t260 - t113 * t213 + t48 * t457 + t323 * t45 + t324 * t47 + t336 * t32 - t337 * (Ifges(4,5) * t293 - Ifges(4,6) * t294) / 0.2e1 - t574 * (-pkin(4) - t497) + t575 * t131 + t576 * t151 + t577 * t149 + (t1 * t323 + t12 * t336 + t3 * t324 + t31 * t576 + t35 * t577 + t575 * t73) * m(7) + t221 * t507 + (Ifges(4,1) * t293 - t472) * t508 + (t365 * t125 + t361 * t126 - t294 * t170 + ((m(6) * t103 - t453) * t361 + (m(6) * t545 + t150 * t364 - t152 * t360 + t213) * t365) * qJD(4) + (0.2e1 * t237 * t508 + t29 * t361 + t30 * t365 + (-t106 * t361 + t107 * t365) * qJD(4)) * m(5)) * pkin(3) - m(6) * (t103 * t112 + t49 * t57 + t50 * t58) - t58 * t150 - t57 * t152 - t164 * mrSges(4,2) + t165 * mrSges(4,1); t552 * t149 + pkin(11) * t374 + (-pkin(11) * t46 + t411) * t360 + ((-Ifges(5,1) / 0.2e1 + t540) * t388 + t546) * t416 + t453 * t107 - m(6) * (t103 * t107 + t49 * t60 + t50 * t61) + ((-pkin(11) * t152 + t410) * t364 + (pkin(5) * t131 - pkin(11) * t150 + t409) * t360) * qJD(5) + t372 - t541 * t388 + t551 * t151 - t106 * t213 + t338 * t45 + t339 * t47 + t354 * t32 + t48 * t496 - t74 * t131 - t61 * t150 - t60 * t152 + t574 * pkin(4) + (t1 * t338 + t12 * t354 + t3 * t339 + (t437 - t74) * t73 + t552 * t35 + t551 * t31) * m(7); t553 + t375 + (-t131 * t204 + t45) * pkin(5) + (t203 * t31 + t204 * t35) * mrSges(7,3) + (t203 * t49 + t204 * t50) * mrSges(6,3) - t73 * (mrSges(7,1) * t204 + mrSges(7,2) * t203) - t103 * (mrSges(6,1) * t204 + mrSges(6,2) * t203) + (-(-t31 + t34) * t35 + (-t204 * t73 + t1) * pkin(5)) * m(7) - t34 * t149 - t49 * t150 + t35 * t151 + t50 * t152 + (t203 * t559 - t584) * t525 + t550 * t524 + (t203 * t557 - t204 * t555) * t520 + (-t204 * t556 + t548 + t586) * t527; -t203 * t149 + t204 * t151 + 0.2e1 * (t12 / 0.2e1 + t35 * t527 + t31 * t524) * m(7) + t32;];
tauc  = t2(:);
