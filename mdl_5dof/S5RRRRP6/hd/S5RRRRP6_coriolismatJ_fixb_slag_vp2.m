% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:53
% EndTime: 2019-12-31 21:53:10
% DurationCPUTime: 7.69s
% Computational Cost: add. (11925->540), mult. (24376->728), div. (0->0), fcn. (24378->6), ass. (0->305)
t335 = sin(qJ(4));
t333 = t335 ^ 2;
t336 = cos(qJ(4));
t334 = t336 ^ 2;
t578 = t333 + t334;
t530 = sin(qJ(3));
t531 = sin(qJ(2));
t532 = cos(qJ(3));
t533 = cos(qJ(2));
t302 = t530 * t531 - t532 * t533;
t493 = t336 * mrSges(6,2);
t497 = t335 * mrSges(6,1);
t557 = m(6) / 0.2e1;
t474 = t302 * t335;
t452 = t533 * pkin(6);
t321 = pkin(7) * t533 + t452;
t449 = t531 * pkin(6);
t362 = -pkin(7) * t531 - t449;
t571 = t532 * t321 + t530 * t362;
t579 = -pkin(4) * t474 + t571;
t608 = (-t493 / 0.2e1 - t497 / 0.2e1) * t302 + t579 * t557;
t303 = -t530 * t533 - t531 * t532;
t472 = t303 * t335;
t446 = mrSges(6,3) * t472;
t216 = -t302 * mrSges(6,2) + t446;
t471 = t303 * t336;
t220 = t302 * mrSges(6,1) + mrSges(6,3) * t471;
t537 = t335 / 0.2e1;
t593 = -t336 / 0.2e1;
t566 = t216 * t593 + t220 * t537;
t607 = t566 + t608;
t606 = -t566 + t608;
t509 = Ifges(6,4) * t336;
t387 = Ifges(6,1) * t335 + t509;
t365 = t336 * t387;
t511 = Ifges(5,4) * t336;
t389 = Ifges(5,1) * t335 + t511;
t366 = t336 * t389;
t510 = Ifges(6,4) * t335;
t383 = Ifges(6,2) * t336 + t510;
t367 = t335 * t383;
t512 = Ifges(5,4) * t335;
t385 = Ifges(5,2) * t336 + t512;
t368 = t335 * t385;
t535 = t336 / 0.2e1;
t542 = -t303 / 0.2e1;
t544 = t302 / 0.2e1;
t545 = -t302 / 0.2e1;
t388 = Ifges(6,1) * t336 - t510;
t139 = -Ifges(6,5) * t303 - t302 * t388;
t390 = Ifges(5,1) * t336 - t512;
t141 = -Ifges(5,5) * t303 - t302 * t390;
t581 = t141 + t139;
t384 = -Ifges(6,2) * t335 + t509;
t135 = -Ifges(6,6) * t303 - t302 * t384;
t386 = -Ifges(5,2) * t335 + t511;
t137 = -Ifges(5,6) * t303 - t302 * t386;
t583 = t137 + t135;
t591 = Ifges(6,5) + Ifges(5,5);
t518 = Ifges(5,6) + Ifges(6,6);
t599 = t518 * t336;
t346 = -Ifges(4,5) * t302 + Ifges(4,6) * t303 + (t365 + t366) * t545 + (t367 + t368) * t544 + (t591 * t335 + t599) * t542 + t581 * t537 + t583 * t535;
t515 = mrSges(5,1) * t336;
t311 = t335 * mrSges(5,2) - t515;
t584 = t571 * t311;
t589 = t571 * mrSges(4,1);
t245 = t321 * t530 - t532 * t362;
t590 = t245 * mrSges(4,2);
t496 = t335 * mrSges(6,2);
t310 = -mrSges(6,1) * t336 + t496;
t598 = t579 * t310;
t605 = t346 + t584 + t590 - t589 + t598;
t604 = t590 / 0.2e1 + t598 / 0.2e1;
t168 = -pkin(4) * t472 + t245;
t603 = t168 * t579;
t331 = Ifges(6,5) * t336;
t381 = Ifges(6,6) * t335 - t331;
t332 = Ifges(5,5) * t336;
t382 = Ifges(5,6) * t335 - t332;
t570 = t381 + t382;
t602 = t302 * t570;
t451 = t532 * pkin(2);
t327 = -t451 - pkin(3);
t520 = t336 * pkin(4);
t307 = t327 - t520;
t601 = t307 * t579;
t328 = -pkin(3) - t520;
t600 = t328 * t579;
t359 = t578 * t532;
t596 = t584 / 0.2e1 - t589 / 0.2e1;
t313 = t493 + t497;
t200 = t313 * t302;
t314 = t335 * mrSges(5,1) + mrSges(5,2) * t336;
t201 = t314 * t302;
t202 = t313 * t303;
t214 = mrSges(6,2) * t303 + mrSges(6,3) * t474;
t215 = mrSges(5,2) * t303 + mrSges(5,3) * t474;
t473 = t302 * t336;
t218 = -t303 * mrSges(6,1) + mrSges(6,3) * t473;
t219 = -t303 * mrSges(5,1) + mrSges(5,3) * t473;
t396 = -pkin(2) * t533 - pkin(1);
t480 = t571 * t303;
t521 = t303 * pkin(8);
t207 = t302 * pkin(3) + t396 + t521;
t85 = t336 * t207 - t335 * t571;
t65 = qJ(5) * t471 + t85;
t54 = t302 * pkin(4) + t65;
t86 = t335 * t207 + t336 * t571;
t66 = qJ(5) * t472 + t86;
t595 = t579 * t202 - mrSges(4,3) * t480 + t168 * t200 + t245 * t201 - t66 * t214 - t86 * t215 - t54 * t218 - t85 * t219 - t396 * (-mrSges(4,1) * t303 - mrSges(4,2) * t302);
t592 = pkin(3) * t571;
t588 = t245 * t530;
t482 = t245 * t571;
t587 = t327 * t571;
t586 = t335 * t245;
t585 = t336 * t245;
t136 = Ifges(6,6) * t302 - t303 * t384;
t138 = Ifges(5,6) * t302 - t303 * t386;
t582 = t138 + t136;
t140 = Ifges(6,5) * t302 - t303 * t388;
t142 = Ifges(5,5) * t302 - t303 * t390;
t580 = t142 + t140;
t577 = t384 + t386;
t576 = t390 + t388;
t552 = m(6) * pkin(4);
t573 = mrSges(6,1) + t552;
t519 = Ifges(6,4) + Ifges(5,4);
t572 = t519 * t335;
t447 = mrSges(5,3) * t472;
t217 = -mrSges(5,2) * t302 + t447;
t221 = t302 * mrSges(5,1) + mrSges(5,3) * t471;
t369 = t221 * t593 - t335 * t217 / 0.2e1;
t523 = pkin(8) * t215;
t222 = -t303 * pkin(3) + t302 * pkin(8);
t450 = t531 * pkin(2);
t209 = t450 + t222;
t101 = t335 * t209 - t585;
t547 = t101 / 0.2e1;
t568 = mrSges(5,3) * t547 + t523 / 0.2e1;
t100 = t336 * t209 + t586;
t375 = -t100 * t335 + t101 * t336;
t567 = -Ifges(4,1) + Ifges(5,3) + Ifges(6,3);
t565 = -t245 * t314 / 0.2e1 - t168 * t313 / 0.2e1;
t272 = mrSges(6,1) * t471;
t394 = mrSges(6,2) * t472 - t272;
t454 = pkin(4) * t471;
t126 = m(6) * t454 - t394;
t282 = t335 * t573 + t493;
t564 = -qJD(1) * t126 + (qJD(2) + qJD(3)) * t282;
t439 = -Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1;
t563 = -t439 * t302 + t140 / 0.2e1 + t142 / 0.2e1;
t109 = t336 * t222 + t586;
t110 = t335 * t222 - t585;
t455 = t552 / 0.2e1;
t551 = -mrSges(6,2) / 0.2e1;
t392 = -t303 * pkin(4) + qJ(5) * t473;
t58 = t109 + t392;
t436 = qJ(5) * t474;
t72 = t436 + t110;
t562 = (t455 + mrSges(6,1) / 0.2e1) * t58 - t110 * mrSges(5,2) / 0.2e1 + t109 * mrSges(5,1) / 0.2e1 + t72 * t551;
t448 = t530 * pkin(2);
t326 = t448 + pkin(8);
t405 = t333 / 0.2e1 + t334 / 0.2e1;
t461 = qJ(5) + t326;
t284 = t461 * t336;
t475 = t284 * t336;
t283 = t461 * t335;
t477 = t283 * t335;
t541 = t311 / 0.2e1;
t561 = t326 * mrSges(5,3) * t405 + (t477 / 0.2e1 + t475 / 0.2e1) * mrSges(6,3) + t327 * t541;
t546 = -t284 / 0.2e1;
t560 = -t283 * t216 / 0.2e1 + t220 * t546 + t369 * t326;
t559 = m(5) / 0.2e1;
t558 = -m(6) / 0.2e1;
t556 = -pkin(3) / 0.2e1;
t555 = pkin(8) / 0.2e1;
t554 = m(5) * pkin(2);
t553 = m(6) * pkin(2);
t55 = t100 + t392;
t550 = -t55 / 0.2e1;
t68 = t436 + t101;
t549 = t68 / 0.2e1;
t548 = -t100 / 0.2e1;
t517 = -qJ(5) - pkin(8);
t312 = t517 * t336;
t540 = t312 / 0.2e1;
t539 = t328 / 0.2e1;
t534 = t336 / 0.4e1;
t529 = m(6) * t168;
t528 = m(6) * t328;
t527 = pkin(3) * t201;
t526 = pkin(3) * t314;
t525 = pkin(4) * t218;
t524 = pkin(4) * t335;
t522 = pkin(8) * t219;
t516 = -t54 + t65;
t514 = mrSges(6,3) * t336;
t513 = Ifges(4,4) * t303;
t506 = Ifges(5,3) * t303;
t505 = Ifges(6,3) * t303;
t287 = Ifges(4,4) * t302;
t438 = -Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t353 = -t138 / 0.2e1 - t136 / 0.2e1 + t438 * t302;
t363 = t314 * t303;
t498 = t302 * mrSges(4,3);
t2 = (t531 ^ 2 - t533 ^ 2) * Ifges(3,4) + (t571 * mrSges(4,3) + t513 + (t141 / 0.2e1 + t139 / 0.2e1 + t439 * t303) * t336 + (-t137 / 0.2e1 - t135 / 0.2e1 - t438 * t303) * t335 + (Ifges(4,2) + t567) * t302) * t303 + pkin(1) * (mrSges(3,1) * t531 + mrSges(3,2) * t533) + t571 * t363 - t245 * t498 - m(5) * (t100 * t85 + t101 * t86 + t482) - m(6) * (t54 * t55 + t66 * t68 + t603) + (t245 * mrSges(4,3) + t353 * t335 + t563 * t336 - t287) * t302 - t68 * t216 - t101 * t217 - t55 * t220 - t100 * t221 + (Ifges(3,2) - Ifges(3,1)) * t533 * t531 + (-m(4) * t396 - mrSges(4,1) * t302 + mrSges(4,2) * t303) * t450 + t595;
t504 = t2 * qJD(1);
t414 = -t471 / 0.2e1;
t417 = -t473 / 0.2e1;
t418 = t474 / 0.2e1;
t3 = t303 * (-Ifges(4,2) * t302 - t513) / 0.2e1 + m(5) * (t109 * t85 + t110 * t86 + t482) + m(6) * (t54 * t58 + t66 * t72 + t603) + (-mrSges(4,3) - t314) * t480 + t72 * t216 + t110 * t217 + t58 * t220 + t109 * t221 + (-0.2e1 * t287 + (-Ifges(4,1) + Ifges(4,2)) * t303) * t545 + (-t505 - t506 + t602) * t544 + t582 * t418 + t580 * t417 + t583 * t472 / 0.2e1 + t581 * t414 + (t302 * t567 + t303 * t570 + t513) * t542 - t595;
t499 = t3 * qJD(1);
t495 = t335 * mrSges(5,3);
t494 = t335 * mrSges(6,3);
t203 = t383 * t303;
t204 = t385 * t303;
t205 = t387 * t303;
t206 = t389 * t303;
t364 = t311 * t303;
t437 = m(6) * t516;
t9 = -t202 * t454 - t168 * t394 - t245 * t364 - t65 * t216 + t66 * t220 + t54 * t446 + t86 * t221 - t66 * t437 + ((-t203 / 0.2e1 - t204 / 0.2e1 - t563) * t335 + (t206 / 0.2e1 + t205 / 0.2e1 + pkin(4) * t529 - t66 * mrSges(6,3) - t86 * mrSges(5,3) + t353) * t336) * t303 + (-t217 + t447) * t85;
t492 = t9 * qJD(1);
t23 = (m(6) * (t335 * t66 + t54 * t336) + t335 * t216 + t336 * t220) * t303;
t491 = qJD(1) * t23;
t488 = t110 * t336;
t485 = t168 * t335;
t478 = t283 * t218;
t476 = t284 * t214;
t470 = t307 * t200;
t309 = t517 * t335;
t469 = t309 * t216;
t468 = t309 * t218;
t467 = t312 * t214;
t466 = t326 * t336;
t465 = t327 * t201;
t464 = t328 * t200;
t460 = t578 * mrSges(6,3);
t458 = qJD(4) * t335;
t456 = mrSges(5,3) * t521;
t453 = -t532 / 0.2e1;
t443 = t524 / 0.2e1;
t441 = -t522 / 0.2e1;
t435 = t335 * t326 * t219;
t434 = t215 * t466;
t433 = -t514 / 0.2e1;
t432 = t514 / 0.2e1;
t428 = -t495 / 0.2e1;
t427 = t495 / 0.2e1;
t426 = -t494 / 0.2e1;
t425 = t494 / 0.2e1;
t424 = t335 * t532;
t423 = t336 * t532;
t404 = (t283 - t309) * t335;
t403 = (t284 - t312) * t336;
t399 = (t307 + t328) * t552;
t397 = t448 / 0.2e1;
t301 = t310 * t524;
t391 = -t334 * t519 - t301;
t380 = -t335 * t54 + t336 * t66;
t257 = t307 * t313;
t276 = t327 * t314;
t350 = t572 + (-Ifges(6,1) + Ifges(6,2) - Ifges(5,1) + Ifges(5,2)) * t336;
t25 = -t257 - t276 + (-t307 * t552 + t350) * t335 + t391;
t341 = pkin(4) * t310 * t414 - t202 * t443 + t65 * t432 + t54 * t433 - t565 + (t427 + t428) * t86 + (t425 + t426) * t66 + 0.2e1 * (t204 + t203) * t534 + (t206 + t205) * t537 - t602 / 0.4e1 - t582 * t335 / 0.4e1 + t580 * t534 + t577 * t472 / 0.4e1 - t576 * t471 / 0.4e1;
t347 = mrSges(5,1) * t548 + mrSges(6,1) * t550 + mrSges(5,2) * t547 + mrSges(6,2) * t549;
t6 = t341 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + (t520 * t558 + t496 / 0.2e1) * t307 + t561) * t303 + (t168 * t443 + t54 * t546 + t284 * t65 / 0.2e1 + pkin(4) * t550) * m(6) + (t335 * t438 - t336 * t439) * t302 - t307 * t272 / 0.2e1 - t525 / 0.2e1 + t347 + t560;
t377 = t6 * qJD(1) - t25 * qJD(2);
t339 = -mrSges(4,2) * t451 + (-mrSges(4,1) + t310 + t311) * t448 + (mrSges(5,3) + mrSges(6,3)) * t359 * pkin(2);
t31 = (t326 * t359 + t327 * t530) * t554 + (t283 * t424 + t284 * t423 + t307 * t530) * t553 + t339;
t374 = -t109 * t335 + t488;
t337 = -m(5) * (t587 + t374 * t326 + (t423 * t86 - t424 * t85 + t588) * pkin(2)) / 0.2e1 + (t601 - t283 * t58 + t284 * t72 + (t168 * t530 + t423 * t66 - t424 * t54) * pkin(2)) * t558 + t478 / 0.2e1 - t476 / 0.2e1 + t470 / 0.2e1 + t465 / 0.2e1 + t109 * t427 - mrSges(5,3) * t488 / 0.2e1 + t435 / 0.2e1 - t434 / 0.2e1 + t58 * t425 + t72 * t433 - (-t202 - t363) * t448 / 0.2e1 + (t220 + t221) * pkin(2) * t424 / 0.2e1 - (t216 + t217) * pkin(2) * t423 / 0.2e1 - t596 - t604;
t338 = (pkin(8) * t375 - t592) * t559 + (t309 * t55 - t312 * t68 + t600) * t557 + t527 / 0.2e1 + t468 / 0.2e1 - t467 / 0.2e1 - t464 / 0.2e1 + t604;
t4 = t100 * t428 + t335 * t441 + t568 * t336 + t55 * t426 + t68 * t432 + t337 + t338 + t596;
t376 = -t4 * qJD(1) + t31 * qJD(2);
t373 = -pkin(4) * t514 + t331 + t332;
t132 = m(6) * (t475 + t477) + t460;
t344 = m(6) * ((-t283 * t336 + t284 * t335) * t303 + t380);
t17 = -t344 / 0.2e1 + t607;
t372 = -qJD(1) * t17 + qJD(2) * t132;
t370 = Ifges(6,2) / 0.4e1 - Ifges(5,1) / 0.4e1 - Ifges(6,1) / 0.4e1 + Ifges(5,2) / 0.4e1;
t361 = t525 / 0.2e1 - t505 / 0.2e1 - t506 / 0.2e1 + t518 * t418 + t591 * t417;
t345 = (mrSges(5,1) + t573) * t453;
t277 = t328 * t313;
t351 = -t257 / 0.2e1 - t276 / 0.2e1 - t277 / 0.2e1 - t301 + t526 / 0.2e1;
t352 = (mrSges(5,2) + mrSges(6,2)) * t453;
t14 = (pkin(2) * t352 - t336 * t519) * t336 + (-t399 / 0.2e1 + t345 * pkin(2) + t350) * t335 + t351;
t29 = t526 - t277 + (-pkin(4) * t528 + t350) * t335 + t391;
t7 = t272 * t539 + t361 - (t220 / 0.2e1 - t437 / 0.2e1) * t312 + (t217 * t555 - t206 / 0.4e1 - t205 / 0.4e1 + t138 / 0.4e1 + t136 / 0.4e1 + (t202 / 0.2e1 - t529 / 0.2e1) * pkin(4) + (pkin(3) * mrSges(5,2) / 0.2e1 + t328 * t551 + t309 * mrSges(6,3) / 0.2e1 + t370 * t335) * t303) * t335 + (t221 * t555 - t204 / 0.4e1 - t203 / 0.4e1 - t142 / 0.4e1 - t140 / 0.4e1 + (-t65 / 0.2e1 + t54 / 0.2e1) * mrSges(6,3) + (mrSges(5,1) * t556 + mrSges(6,3) * t540 - t370 * t336 - t572 + (t528 / 0.2e1 + t310 / 0.2e1) * pkin(4)) * t303) * t336 + (t381 / 0.4e1 + t382 / 0.4e1) * t302 - t405 * t456 - t469 / 0.2e1 + t562 + t565;
t356 = -t7 * qJD(1) - t14 * qJD(2) - t29 * qJD(3);
t177 = m(6) * (-t309 * t335 - t312 * t336) + t460;
t343 = m(6) * ((t309 * t336 - t312 * t335) * t303 + t380);
t19 = -t343 / 0.2e1 + t607;
t70 = (t397 - t403 / 0.2e1 - t404 / 0.2e1) * m(6) - t460;
t354 = -qJD(1) * t19 - qJD(2) * t70 + qJD(3) * t177;
t340 = t341 + t361;
t275 = t282 * qJD(5);
t274 = t282 * qJD(4);
t71 = (t403 + t404) * t557 + m(6) * t397 + t460;
t20 = t343 / 0.2e1 + t606;
t18 = t344 / 0.2e1 + t606;
t15 = t365 / 0.2e1 - t367 / 0.2e1 + t366 / 0.2e1 - t368 / 0.2e1 + (t335 * t345 + t336 * t352) * pkin(2) - t351 + t577 * t535 + (t399 + t576) * t537;
t8 = t340 + t364 * t556 + (-t516 * t312 + (-t328 * t471 + t485) * pkin(4)) * t557 + t394 * t539 + t469 / 0.2e1 + t220 * t540 + t578 * t456 / 0.2e1 + (t309 * t426 - t312 * t432) * t303 + t369 * pkin(8) + t562;
t5 = t340 + (t516 * t284 + (-t307 * t471 + t485) * pkin(4)) * t557 + t307 * t394 / 0.2e1 + t55 * t455 - t347 + t561 * t303 + t560;
t1 = t338 - t337 + t346 + (-mrSges(4,1) / 0.2e1 + t541) * t571 + (mrSges(6,3) * t549 + t568) * t336 + (mrSges(5,3) * t548 + mrSges(6,3) * t550 + t441) * t335;
t10 = [-qJD(2) * t2 + qJD(3) * t3 - qJD(4) * t9 + qJD(5) * t23, -t504 + (m(6) * (-t283 * t55 + t284 * t68 + t601) + m(4) * (-t532 * t571 - t588) * pkin(2) - t465 + m(5) * (t326 * t375 + t587) + Ifges(3,5) * t533 - Ifges(3,6) * t531 - t55 * t494 - mrSges(3,1) * t452 + t68 * t514 + t451 * t498 + mrSges(3,2) * t449 + t375 * mrSges(5,3) + t434 - t435 - t470 + t476 - t478 + t303 * mrSges(4,3) * t448 + t605) * qJD(2) + t1 * qJD(3) + t5 * qJD(4) + t18 * qJD(5), t1 * qJD(2) + t8 * qJD(4) + t20 * qJD(5) + t499 + (-t464 - t467 + t468 + t527 + (t110 * mrSges(5,3) + t72 * mrSges(6,3) + t523) * t336 + (-t109 * mrSges(5,3) - t58 * mrSges(6,3) - t522) * t335 + 0.2e1 * (pkin(8) * t374 - t592) * t559 + 0.2e1 * (t309 * t58 - t312 * t72 + t600) * t557 + t605) * qJD(3), t5 * qJD(2) + t8 * qJD(3) - t492 + (-t86 * mrSges(5,1) - t85 * mrSges(5,2) - t65 * mrSges(6,2) + (t599 + (-mrSges(6,3) * pkin(4) + t591) * t335) * t303 - t573 * t66) * qJD(4), qJD(2) * t18 + qJD(3) * t20 + t491; -qJD(3) * t4 + qJD(4) * t6 - qJD(5) * t17 + t504, qJD(3) * t31 - qJD(4) * t25 + qJD(5) * t132, ((-pkin(3) * t530 + pkin(8) * t359) * t554 + (-t309 * t424 - t312 * t423 + t328 * t530) * t553 + t339) * qJD(3) + t15 * qJD(4) + t71 * qJD(5) + t376, t15 * qJD(3) + (-mrSges(5,1) * t466 + t283 * mrSges(6,2) - t284 * t573 + t373) * qJD(4) + (mrSges(5,2) * t326 - t518) * t458 + t377, qJD(3) * t71 + t372; qJD(2) * t4 - qJD(4) * t7 - qJD(5) * t19 - t499, -qJD(4) * t14 - qJD(5) * t70 - t376, -qJD(4) * t29 + qJD(5) * t177, (-t309 * mrSges(6,2) - pkin(8) * t515 + t312 * t573 + t373) * qJD(4) + (mrSges(5,2) * pkin(8) - t518) * t458 + t356, t354; -qJD(2) * t6 + qJD(3) * t7 + qJD(5) * t126 + t492, qJD(3) * t14 - t275 - t377, -t275 - t356, 0, -t564; qJD(2) * t17 + qJD(3) * t19 - qJD(4) * t126 - t491, qJD(3) * t70 + t274 - t372, t274 - t354, t564, 0;];
Cq = t10;
