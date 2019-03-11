% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:50
% EndTime: 2019-03-08 21:23:03
% DurationCPUTime: 7.48s
% Computational Cost: add. (11997->517), mult. (27287->708), div. (0->0), fcn. (29907->10), ass. (0->269)
t356 = sin(qJ(5));
t350 = t356 ^ 2;
t359 = cos(qJ(5));
t352 = t359 ^ 2;
t440 = t350 + t352;
t565 = mrSges(6,3) + mrSges(7,3);
t571 = t565 * t440;
t570 = -m(7) / 0.2e1;
t357 = sin(qJ(3));
t360 = cos(qJ(3));
t474 = sin(pkin(11));
t475 = cos(pkin(11));
t317 = t357 * t474 - t360 * t475;
t377 = t357 * t475 + t360 * t474;
t454 = t377 * t359;
t242 = mrSges(7,1) * t317 - mrSges(7,3) * t454;
t510 = -t242 / 0.2e1;
t501 = t359 / 0.2e1;
t422 = t474 * pkin(3);
t340 = t422 + pkin(9);
t444 = qJ(6) + t340;
t301 = t444 * t356;
t302 = t444 * t359;
t552 = t301 * t356 + t302 * t359;
t569 = m(7) * t552;
t481 = t317 * mrSges(5,3);
t479 = t377 * mrSges(5,3);
t424 = -pkin(3) * t360 - pkin(2);
t568 = -m(5) * t424 - mrSges(5,1) * t317 - mrSges(5,2) * t377;
t493 = t317 * pkin(5);
t231 = pkin(4) * t317 - pkin(9) * t377 + t424;
t491 = -qJ(4) - pkin(8);
t327 = t491 * t357;
t328 = t491 * t360;
t549 = t474 * t327 - t475 * t328;
t114 = t359 * t231 - t356 * t549;
t87 = -qJ(6) * t454 + t114;
t72 = t87 + t493;
t115 = t231 * t356 + t359 * t549;
t455 = t377 * t356;
t88 = -qJ(6) * t455 + t115;
t392 = t356 * t72 - t359 * t88;
t238 = -mrSges(7,2) * t317 - mrSges(7,3) * t455;
t443 = t238 * t501 + t356 * t510;
t567 = (-(-t301 * t359 + t302 * t356) * t377 - t392) * t570 - t443;
t492 = Ifges(6,6) + Ifges(7,6);
t564 = t492 * t359;
t482 = t317 * mrSges(6,2);
t239 = -mrSges(6,3) * t455 - t482;
t483 = t317 * mrSges(6,1);
t243 = -mrSges(6,3) * t454 + t483;
t509 = -t243 / 0.2e1;
t562 = t239 * t501 + t356 * t509 + t443;
t436 = m(7) / 0.4e1 + m(6) / 0.4e1;
t561 = 0.4e1 * t436;
t528 = mrSges(7,3) / 0.2e1;
t361 = cos(qJ(2));
t560 = m(4) * t361;
t423 = t475 * pkin(3);
t341 = -t423 - pkin(4);
t324 = -t359 * pkin(5) + t341;
t497 = m(7) * t324;
t559 = mrSges(6,2) + mrSges(7,2);
t558 = Ifges(6,5) + Ifges(7,5);
t557 = Ifges(6,3) + Ifges(7,3);
t532 = m(7) * pkin(5);
t556 = -mrSges(7,1) - t532;
t555 = t340 * t440;
t554 = t238 + t239;
t553 = t242 + t243;
t343 = t356 * mrSges(7,1);
t344 = t359 * mrSges(7,2);
t441 = t344 + t343;
t345 = Ifges(7,5) * t359;
t346 = Ifges(6,5) * t359;
t551 = -t345 - t346;
t347 = Ifges(7,4) * t359;
t332 = Ifges(7,1) * t356 + t347;
t348 = Ifges(6,4) * t359;
t333 = Ifges(6,1) * t356 + t348;
t351 = t357 ^ 2;
t353 = t360 ^ 2;
t550 = t351 + t353;
t496 = pkin(3) * t357;
t233 = pkin(4) * t377 + pkin(9) * t317 + t496;
t258 = -t475 * t327 - t328 * t474;
t124 = t359 * t233 + t258 * t356;
t125 = t356 * t233 - t258 * t359;
t545 = -t124 * t356 + t125 * t359;
t329 = t357 * mrSges(4,1) + t360 * mrSges(4,2);
t544 = -mrSges(4,1) * t360 + mrSges(4,2) * t357;
t221 = t441 * t317;
t477 = t359 * mrSges(6,2);
t400 = t356 * mrSges(6,1) + t477;
t222 = t400 * t317;
t543 = -t221 / 0.2e1 - t222 / 0.2e1;
t533 = m(5) * pkin(3);
t541 = t474 * t533 - mrSges(5,2);
t395 = Ifges(7,2) * t356 - t347;
t153 = t317 * Ifges(7,6) - t377 * t395;
t396 = Ifges(6,2) * t356 - t348;
t154 = t317 * Ifges(6,6) - t377 * t396;
t540 = -(Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1) * t317 - t153 / 0.2e1 - t154 / 0.2e1;
t437 = t532 / 0.2e1;
t405 = t437 + mrSges(6,1) / 0.2e1;
t539 = t477 / 0.2e1 + t405 * t356;
t326 = -mrSges(6,1) * t359 + mrSges(6,2) * t356;
t538 = m(6) * t341 - t475 * t533 - mrSges(5,1) + t326;
t355 = cos(pkin(6));
t354 = sin(pkin(6));
t358 = sin(qJ(2));
t451 = t354 * t358;
t298 = t355 * t360 - t357 * t451;
t450 = t354 * t360;
t299 = t355 * t357 + t358 * t450;
t378 = t298 * t474 + t299 * t475;
t449 = t354 * t361;
t149 = -t356 * t378 - t359 * t449;
t150 = -t356 * t449 + t359 * t378;
t457 = t317 * t356;
t173 = -pkin(5) * t457 + t549;
t174 = pkin(5) * t455 + t258;
t197 = -t475 * t298 + t299 * t474;
t224 = t400 * t377;
t236 = -mrSges(7,2) * t377 + mrSges(7,3) * t457;
t237 = -mrSges(6,2) * t377 + mrSges(6,3) * t457;
t456 = t317 * t359;
t241 = mrSges(6,1) * t377 + mrSges(6,3) * t456;
t389 = t114 * t356 - t115 * t359;
t240 = mrSges(7,1) * t377 + mrSges(7,3) * t456;
t512 = t240 / 0.2e1;
t223 = t441 * t377;
t517 = t223 / 0.2e1;
t534 = m(7) / 0.2e1;
t535 = m(6) / 0.2e1;
t536 = m(5) / 0.2e1;
t73 = pkin(5) * t377 + qJ(6) * t456 + t124;
t91 = qJ(6) * t457 + t125;
t537 = -t449 * t496 * t536 + (t124 * t149 + t125 * t150 + t258 * t378 + (t389 + t549) * t197) * t535 + (t149 * t73 + t150 * t91 + t174 * t378 + (t173 + t392) * t197) * t534 + (t241 / 0.2e1 + t512) * t149 + (t236 / 0.2e1 + t237 / 0.2e1) * t150 + (t517 + t224 / 0.2e1) * t378;
t531 = -mrSges(7,1) / 0.2e1;
t530 = mrSges(6,2) / 0.2e1;
t529 = mrSges(7,2) / 0.2e1;
t527 = -t73 / 0.2e1;
t219 = mrSges(7,1) * t454 - mrSges(7,2) * t455;
t521 = t219 / 0.2e1;
t220 = t377 * t326;
t520 = -t220 / 0.2e1;
t513 = t238 / 0.2e1;
t508 = t258 / 0.2e1;
t267 = t377 * t449;
t507 = t267 / 0.2e1;
t505 = t377 / 0.2e1;
t325 = -mrSges(7,1) * t359 + mrSges(7,2) * t356;
t504 = -t325 / 0.2e1;
t503 = -t340 / 0.2e1;
t502 = t356 / 0.2e1;
t500 = m(7) * t174;
t498 = m(7) * t377;
t495 = pkin(5) * t356;
t488 = mrSges(7,3) * t359;
t485 = Ifges(6,4) * t356;
t484 = Ifges(7,4) * t356;
t480 = t377 * mrSges(5,1);
t438 = t533 / 0.2e1;
t364 = (-t317 * t555 + t341 * t377) * t535 + (-t317 * t552 + t324 * t377) * t534 - t377 * t504 + t326 * t505 + (-t317 * t474 - t377 * t475) * t438 - t317 * t571 / 0.2e1;
t368 = (t124 * t359 + t125 * t356) * t535 + (t356 * t91 + t359 * t73) * t534 + t357 * t438 + (t236 + t237) * t502 + (t240 + t241) * t501;
t300 = t317 * mrSges(5,2);
t406 = -t300 + t480;
t12 = t364 - t368 - t406;
t473 = qJD(2) * t12;
t446 = t359 * t150;
t448 = t356 * t149;
t388 = -t446 + t448;
t14 = (t378 + t388) * t197 * t561;
t470 = t14 * qJD(1);
t268 = t317 * t449;
t214 = t268 * t356 + t359 * t451;
t215 = -t268 * t359 + t356 * t451;
t452 = t354 ^ 2 * t358;
t467 = t197 * t267;
t18 = (t149 * t214 + t150 * t215 + t467) * t561 + (-t298 * t354 * t357 + t299 * t450 - t452) * t560 + (-t268 * t378 - t361 * t452 + t467) * m(5);
t469 = t18 * qJD(1);
t466 = t197 * t377;
t463 = t214 * t356;
t462 = t215 * t359;
t460 = t258 * t267;
t459 = t258 * t377;
t453 = t340 * t359;
t417 = t457 / 0.2e1;
t442 = mrSges(7,1) * t417 + t456 * t529;
t138 = -0.2e1 * (t350 / 0.4e1 + t352 / 0.4e1 + 0.1e1 / 0.4e1) * t498;
t439 = t138 * qJD(2);
t434 = m(7) * t507;
t433 = t530 + t529;
t432 = mrSges(6,3) / 0.2e1 + t528;
t431 = -Ifges(6,1) / 0.2e1 - Ifges(7,1) / 0.2e1;
t430 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t429 = -Ifges(6,2) / 0.2e1 - Ifges(7,2) / 0.2e1;
t427 = t87 / 0.2e1 - t72 / 0.2e1;
t426 = m(7) * (-t72 + t87);
t303 = -m(7) * t495 - t441;
t416 = -t456 / 0.2e1;
t415 = -t449 / 0.2e1;
t411 = t239 / 0.2e1 + t513;
t410 = t509 + t510;
t409 = t346 / 0.4e1 + t345 / 0.4e1;
t407 = t454 * t532;
t404 = t426 / 0.2e1;
t403 = -t223 - t224 - t479;
t402 = mrSges(7,2) * t416 + t173 * t534 + t457 * t531;
t398 = -Ifges(6,1) * t359 + t485;
t397 = -Ifges(7,1) * t359 + t484;
t331 = Ifges(6,2) * t359 + t485;
t330 = Ifges(7,2) * t359 + t484;
t362 = (t267 * t341 + (t462 - t463) * t340) * t535 + (-t214 * t301 + t215 * t302 + t267 * t324) * t534 - t267 * mrSges(5,1) / 0.2e1 + t268 * mrSges(5,2) / 0.2e1 + (-t267 * t475 - t268 * t474) * t438 + (t325 + t326) * t507 + t329 * t415 + t565 * (-t463 / 0.2e1 + t462 / 0.2e1);
t2 = -t362 + (t329 + t406) * t415 + t537 + (t543 - t562) * t197;
t155 = t317 * Ifges(7,5) - t377 * t397;
t156 = t317 * Ifges(6,5) - t377 * t398;
t376 = t156 / 0.2e1 + t155 / 0.2e1 + t430 * t317;
t3 = pkin(2) * t329 - t424 * t406 - t72 * t240 - t73 * t242 - t88 * t236 - t91 * t238 - t114 * t241 - t115 * t237 - t124 * t243 - t125 * t239 - t173 * t223 + t174 * t221 + (-Ifges(4,1) + Ifges(4,2)) * t360 * t357 - m(6) * (t114 * t124 + t115 * t125 + t258 * t549) - m(7) * (t173 * t174 + t72 * t73 + t88 * t91) + (-Ifges(5,4) * t317 + t540 * t356 + t376 * t359) * t317 + (-t353 + t351) * Ifges(4,4) - (-(t492 * t356 + Ifges(5,4) + t551) * t377 + (Ifges(5,2) - Ifges(5,1) + t431 * t352 + (t429 * t356 + (Ifges(6,4) + Ifges(7,4)) * t359) * t356 + t557) * t317) * t377 - t549 * t224 + t568 * t496 + t258 * t222;
t394 = t2 * qJD(1) - t3 * qJD(2);
t225 = t377 * t330;
t226 = t377 * t331;
t227 = t332 * t377;
t228 = t333 * t377;
t6 = t114 * t239 - t115 * t243 + t174 * t219 - t258 * t220 + t87 * t238 + (t426 - t242) * t88 - ((-t72 * mrSges(7,3) - t114 * mrSges(6,3) - t226 / 0.2e1 - t225 / 0.2e1 + t376) * t356 + (t88 * mrSges(7,3) + t115 * mrSges(6,3) + t228 / 0.2e1 + t227 / 0.2e1 + (-t223 - t500) * pkin(5) - t540) * t359) * t377;
t382 = t510 + t404;
t367 = (t520 + t521 + t407 / 0.2e1) * t197 + (t432 * t455 + t411) * t149 + (-t432 * t454 + t382 + t509) * t150;
t384 = t531 - t405;
t7 = t214 * t384 + t215 * t433 + t367;
t393 = t7 * qJD(1) + t6 * qJD(2);
t11 = -t403 * t377 + (t553 * t356 - t554 * t359 + t481) * t317 + m(7) * (t174 * t377 + t317 * t392) + m(6) * (t317 * t389 + t459) + m(5) * (-t317 * t549 + t459);
t369 = t451 * t536 + (t534 + t535) * (t359 * t214 + t356 * t215);
t373 = (-t317 * t378 + t466) * t536 + 0.2e1 * t436 * (t317 * t388 + t466);
t20 = -t369 + t373;
t391 = -qJD(1) * t20 - qJD(2) * t11;
t31 = (m(7) * (t356 * t88 + t72 * t359) + t356 * t238 + t359 * t242) * t377;
t380 = (t149 * t359 + t150 * t356) * t498;
t45 = t380 / 0.2e1 + t434;
t390 = -qJD(1) * t45 - qJD(2) * t31;
t148 = -t219 - t407;
t383 = qJD(2) * t148 + qJD(3) * t303;
t15 = (t482 / 0.2e1 - t411) * t359 + (t483 / 0.2e1 + (t493 / 0.2e1 - t427) * m(7) - t410) * t356 + t442;
t381 = t15 * qJD(2);
t157 = t440 * mrSges(7,3) + t569;
t27 = t402 + t567;
t48 = 0.2e1 * (t378 / 0.4e1 + t448 / 0.4e1 - t446 / 0.4e1) * m(7);
t379 = -qJD(1) * t48 - qJD(2) * t27 + qJD(3) * t157;
t34 = -t325 * t495 - t324 * t441 - t341 * t400 + (-t332 / 0.2e1 - t347 / 0.2e1 - t333 / 0.2e1 - t348 / 0.2e1) * t359 + (-pkin(5) * t497 + t330 / 0.2e1 + t331 / 0.2e1 + (Ifges(7,4) / 0.2e1 + Ifges(6,4) / 0.2e1) * t356 + (-t429 + t431) * t359) * t356;
t365 = (t350 / 0.2e1 + t352 / 0.2e1) * t340 * mrSges(6,3) + (t348 / 0.4e1 + t347 / 0.4e1 + t333 / 0.4e1 + t332 / 0.4e1 + t301 * t528 + (-Ifges(6,2) / 0.4e1 - Ifges(7,2) / 0.4e1) * t356) * t356 + (t331 / 0.4e1 + t330 / 0.4e1 + t302 * t528 + (-Ifges(6,1) / 0.4e1 - Ifges(7,1) / 0.4e1) * t359 + (Ifges(6,4) / 0.4e1 + Ifges(7,4) / 0.4e1) * t356 + (-t497 / 0.2e1 + t504) * pkin(5)) * t359;
t366 = t382 * t302 + (t243 * t503 + mrSges(6,2) * t508 - t226 / 0.4e1 - t225 / 0.4e1 + t156 / 0.4e1 + t155 / 0.4e1 + t427 * mrSges(7,3)) * t359 + t174 * t441 / 0.2e1 - t301 * t513 + t324 * t521 + t341 * t520;
t370 = t239 * t503 + mrSges(6,1) * t508 - t228 / 0.4e1 - t227 / 0.4e1 - t154 / 0.4e1 - t153 / 0.4e1 + (t517 + t500 / 0.2e1) * pkin(5);
t374 = -t124 * mrSges(6,1) / 0.2e1 + t125 * t530 + mrSges(7,1) * t527 + t91 * t529;
t5 = (m(7) * t527 - t240 / 0.2e1) * pkin(5) + t370 * t356 + (t430 * t359 + (-0.3e1 / 0.4e1 * Ifges(6,6) - 0.3e1 / 0.4e1 * Ifges(7,6)) * t356 + t409) * t317 - (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + t365) * t377 + t366 + t374;
t371 = (t356 * t384 - t359 * t433) * t197;
t9 = (-t344 / 0.2e1 - t343 / 0.2e1 - t539) * t197 - t371;
t375 = t9 * qJD(1) - t5 * qJD(2) + t34 * qJD(3);
t137 = (-0.1e1 / 0.2e1 + t440 / 0.2e1) * t498;
t49 = t378 * t534 + t388 * t570;
t46 = t434 - t380 / 0.2e1;
t28 = t402 - t567;
t19 = t369 + t373;
t16 = t539 * t317 + t356 * t404 + t442 + t562;
t13 = t364 + t368;
t10 = t356 * t197 * t437 - t371 + (t400 + t441) * t197 / 0.2e1;
t8 = t214 * t437 + t367 + (mrSges(6,1) + mrSges(7,1)) * t214 / 0.2e1 - t559 * t215 / 0.2e1;
t4 = t366 - t365 * t377 + t73 * t437 + t409 * t317 + pkin(5) * t512 + ((-Ifges(6,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t317 + t370) * t356 - t374 + t557 * t505 + t492 * t417 + t558 * t416;
t1 = t362 + (-t480 / 0.2e1 + t300 / 0.2e1 - t329 / 0.2e1) * t449 + t537 + (-t356 * t410 - t359 * t411 + t543) * t197;
t17 = [qJD(2) * t18 + qJD(3) * t14, t1 * qJD(3) + t19 * qJD(4) + t8 * qJD(5) + t46 * qJD(6) + t469 + (t553 * t214 + t554 * t215 + t268 * t481 - t403 * t267 + 0.2e1 * (t174 * t267 + t214 * t72 + t215 * t88) * t534 + 0.2e1 * (t114 * t214 + t115 * t215 + t460) * t535 + 0.2e1 * (-t268 * t549 + t460) * t536 + ((t550 * mrSges(4,3) - mrSges(3,2)) * t361 + t550 * pkin(8) * t560 + (-m(4) * pkin(2) - mrSges(3,1) + t544 - t568) * t358) * t354) * qJD(2), t1 * qJD(2) + t10 * qJD(5) + t49 * qJD(6) + t470 + (-t299 * mrSges(4,1) - t298 * mrSges(4,2) + (t325 + t538 + t497) * t378 + (-m(6) * t555 - t541 - t569 - t571) * t197) * qJD(3), qJD(2) * t19, t8 * qJD(2) + t10 * qJD(3) + (-t559 * t149 + (-mrSges(6,1) + t556) * t150) * qJD(5), qJD(2) * t46 + qJD(3) * t49; qJD(3) * t2 + qJD(4) * t20 + qJD(5) * t7 - qJD(6) * t45 - t469, -qJD(3) * t3 + qJD(4) * t11 + qJD(5) * t6 - qJD(6) * t31 (m(7) * (t173 * t324 - t301 * t73 + t302 * t91) - t73 * t356 * mrSges(7,3) + t544 * pkin(8) + (m(6) * t545 - t356 * t241) * t340 + t545 * mrSges(6,3) - t422 * t479 + t423 * t481 + t91 * t488 + (t492 * t377 + (t395 + t396) * t317) * t501 + (t558 * t377 + (t397 + t398) * t317) * t502 - Ifges(5,6) * t377 + t538 * t549 - t541 * t258 + Ifges(4,5) * t360 - Ifges(4,6) * t357 - t341 * t222 - t324 * t221 + t173 * t325 - Ifges(5,5) * t317 - t301 * t240 + t302 * t236 + (t558 * t356 + t564) * t505 + (t333 + t332) * t416 + (t331 + t330) * t417 + t237 * t453) * qJD(3) + t13 * qJD(4) + t4 * qJD(5) + t28 * qJD(6) + t394, qJD(3) * t13 + qJD(5) * t16 - qJD(6) * t137 - t391, t4 * qJD(3) + t16 * qJD(4) + t393 + (-mrSges(6,1) * t115 - mrSges(6,2) * t114 - mrSges(7,2) * t87 - (t564 + (-mrSges(7,3) * pkin(5) + t558) * t356) * t377 + t556 * t88) * qJD(5), qJD(3) * t28 - qJD(4) * t137 + t390; -qJD(2) * t2 - qJD(5) * t9 - qJD(6) * t48 - t470, qJD(4) * t12 + qJD(5) * t5 - qJD(6) * t27 - t394, -qJD(5) * t34 + qJD(6) * t157, t473, -t375 + (-mrSges(6,1) * t453 + mrSges(7,2) * t301 - pkin(5) * t488 + (mrSges(6,2) * t340 - t492) * t356 + t556 * t302 - t551) * qJD(5), t379; -qJD(2) * t20, -qJD(3) * t12 - qJD(5) * t15 + qJD(6) * t138 + t391, -t473, 0 (-t400 + t303) * qJD(5) - t381, t439; -qJD(2) * t7 + qJD(3) * t9, -qJD(3) * t5 + qJD(4) * t15 + qJD(6) * t148 - t393, t303 * qJD(6) + t375, t381, 0, t383; qJD(2) * t45 + qJD(3) * t48, qJD(3) * t27 - qJD(4) * t138 - qJD(5) * t148 - t390, -qJD(5) * t303 - t379, -t439, -t383, 0;];
Cq  = t17;
