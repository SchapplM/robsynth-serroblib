% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR9
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:38
% EndTime: 2019-03-09 10:57:07
% DurationCPUTime: 55.74s
% Computational Cost: add. (27436->978), mult. (67011->1350), div. (0->0), fcn. (55847->18), ass. (0->427)
t394 = pkin(12) + qJ(6);
t389 = sin(t394);
t391 = cos(t394);
t399 = cos(pkin(12));
t387 = pkin(5) * t399 + pkin(4);
t396 = sin(pkin(12));
t445 = -mrSges(6,1) * t399 + mrSges(6,2) * t396;
t415 = -m(6) * pkin(4) - m(7) * t387 - mrSges(5,1) + t445;
t670 = -mrSges(7,1) * t391 + mrSges(7,2) * t389 + t415;
t397 = sin(pkin(11));
t400 = cos(pkin(11));
t405 = sin(qJ(4));
t409 = cos(qJ(4));
t350 = t397 * t409 + t400 * t405;
t398 = sin(pkin(6));
t410 = cos(qJ(2));
t498 = t398 * t410;
t420 = t350 * t498;
t287 = qJD(1) * t420;
t337 = t350 * qJD(4);
t669 = -t287 + t337;
t406 = sin(qJ(2));
t436 = pkin(2) * t406 - qJ(3) * t410;
t488 = qJD(1) * t398;
t325 = t436 * t488;
t469 = t406 * t488;
t401 = cos(pkin(6));
t487 = qJD(1) * t401;
t479 = pkin(1) * t487;
t326 = -pkin(8) * t469 + t410 * t479;
t248 = t400 * t325 - t397 * t326;
t494 = t400 * t410;
t424 = (pkin(3) * t406 - pkin(9) * t494) * t398;
t206 = qJD(1) * t424 + t248;
t249 = t397 * t325 + t400 * t326;
t468 = t410 * t488;
t455 = t397 * t468;
t222 = -pkin(9) * t455 + t249;
t348 = t397 * t405 - t409 * t400;
t536 = pkin(9) + qJ(3);
t357 = t536 * t397;
t359 = t536 * t400;
t613 = -t409 * t357 - t359 * t405;
t618 = -t348 * qJD(3) + qJD(4) * t613 - t405 * t206 - t409 * t222;
t535 = pkin(10) + qJ(5);
t422 = -m(6) * qJ(5) - m(7) * t535 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t668 = qJ(5) * t469 - t618;
t374 = pkin(8) * t468;
t327 = t406 * t479 + t374;
t278 = pkin(3) * t455 + t327;
t419 = t348 * t498;
t288 = qJD(1) * t419;
t336 = t348 * qJD(4);
t667 = -qJD(5) * t350 - t278 + (-t288 + t336) * qJ(5) + t669 * pkin(4);
t365 = qJD(4) - t468;
t546 = t365 / 0.2e1;
t459 = qJD(2) + t487;
t305 = t397 * t459 + t400 * t469;
t414 = t397 * t469 - t400 * t459;
t411 = t409 * t305 - t405 * t414;
t556 = t411 / 0.2e1;
t235 = t405 * t305 + t409 * t414;
t558 = t235 / 0.2e1;
t559 = -t235 / 0.2e1;
t231 = qJD(6) + t235;
t560 = t231 / 0.2e1;
t197 = t365 * t396 + t399 * t411;
t562 = t197 / 0.2e1;
t460 = t399 * t365 - t396 * t411;
t564 = t460 / 0.2e1;
t404 = sin(qJ(6));
t408 = cos(qJ(6));
t122 = t197 * t408 + t404 * t460;
t574 = t122 / 0.2e1;
t640 = -t197 * t404 + t408 * t460;
t576 = t640 / 0.2e1;
t280 = -pkin(2) * t459 + qJD(3) - t326;
t232 = pkin(3) * t414 + t280;
t104 = t235 * pkin(4) - qJ(5) * t411 + t232;
t291 = qJ(3) * t459 + t327;
t313 = (-pkin(2) * t410 - qJ(3) * t406 - pkin(1)) * t398;
t296 = qJD(1) * t313;
t210 = -t397 * t291 + t400 * t296;
t165 = -pkin(3) * t468 - t305 * pkin(9) + t210;
t211 = t400 * t291 + t397 * t296;
t179 = -pkin(9) * t414 + t211;
t98 = t165 * t405 + t179 * t409;
t95 = qJ(5) * t365 + t98;
t51 = t399 * t104 - t396 * t95;
t33 = pkin(5) * t235 - pkin(10) * t197 + t51;
t52 = t396 * t104 + t399 * t95;
t44 = pkin(10) * t460 + t52;
t13 = t33 * t408 - t404 * t44;
t14 = t33 * t404 + t408 * t44;
t653 = t232 * mrSges(5,1) + t51 * mrSges(6,1) + t13 * mrSges(7,1) - t52 * mrSges(6,2) - t14 * mrSges(7,2);
t525 = Ifges(5,4) * t411;
t146 = -t235 * Ifges(5,2) + t365 * Ifges(5,6) + t525;
t665 = mrSges(5,3) * t98 + t146 / 0.2e1;
t666 = -Ifges(5,4) * t556 + Ifges(6,5) * t562 + Ifges(7,5) * t574 - Ifges(5,2) * t559 - Ifges(5,6) * t546 + Ifges(6,6) * t564 + Ifges(7,6) * t576 + Ifges(6,3) * t558 + Ifges(7,3) * t560 + t653 - t665;
t484 = qJD(2) * t410;
t331 = (qJD(1) * t484 + qJDD(1) * t406) * t398;
t480 = qJDD(1) * t401;
t379 = qJDD(2) + t480;
t266 = -t331 * t397 + t379 * t400;
t267 = t331 * t400 + t379 * t397;
t133 = -qJD(4) * t235 + t405 * t266 + t409 * t267;
t485 = qJD(2) * t406;
t467 = t398 * t485;
t642 = -qJD(1) * t467 + qJDD(1) * t498;
t318 = qJDD(4) - t642;
t109 = -t133 * t396 + t318 * t399;
t579 = t109 / 0.2e1;
t110 = t133 * t399 + t318 * t396;
t578 = t110 / 0.2e1;
t134 = qJD(4) * t411 - t409 * t266 + t405 * t267;
t570 = t134 / 0.2e1;
t250 = t288 * t396 + t399 * t469;
t646 = -t396 * t336 + t250;
t251 = -t288 * t399 + t396 * t469;
t664 = t399 * t336 + t251;
t539 = pkin(1) * t401;
t478 = qJD(2) * t539;
t458 = qJD(1) * t478;
t475 = pkin(1) * t480;
t258 = pkin(8) * t642 + t406 * t475 + t410 * t458;
t219 = t379 * qJ(3) + qJD(3) * t459 + t258;
t483 = qJD(3) * t406;
t227 = -pkin(2) * t642 - qJ(3) * t331 + (-pkin(1) * qJDD(1) - qJD(1) * t483) * t398;
t143 = -t219 * t397 + t400 * t227;
t101 = -pkin(3) * t642 - pkin(9) * t267 + t143;
t144 = t400 * t219 + t397 * t227;
t111 = pkin(9) * t266 + t144;
t481 = qJD(4) * t409;
t482 = qJD(4) * t405;
t31 = t405 * t101 + t409 * t111 + t165 * t481 - t179 * t482;
t27 = qJ(5) * t318 + qJD(5) * t365 + t31;
t500 = t398 * t406;
t380 = pkin(8) * t500;
t259 = -qJD(2) * t374 - qJDD(1) * t380 - t406 * t458 + t410 * t475;
t239 = -t379 * pkin(2) + qJDD(3) - t259;
t178 = -t266 * pkin(3) + t239;
t49 = t134 * pkin(4) - t133 * qJ(5) - qJD(5) * t411 + t178;
t10 = -t27 * t396 + t399 * t49;
t11 = t399 * t27 + t396 * t49;
t551 = t318 / 0.2e1;
t571 = -t134 / 0.2e1;
t572 = t133 / 0.2e1;
t132 = qJDD(6) + t134;
t573 = t132 / 0.2e1;
t40 = -qJD(6) * t122 + t109 * t408 - t110 * t404;
t588 = t40 / 0.2e1;
t39 = qJD(6) * t640 + t109 * t404 + t110 * t408;
t589 = t39 / 0.2e1;
t5 = pkin(5) * t134 - pkin(10) * t110 + t10;
t6 = pkin(10) * t109 + t11;
t1 = qJD(6) * t13 + t404 * t5 + t408 * t6;
t2 = -qJD(6) * t14 - t404 * t6 + t408 * t5;
t597 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t7 = Ifges(7,5) * t39 + Ifges(7,6) * t40 + Ifges(7,3) * t132;
t663 = t597 + mrSges(5,1) * t178 + mrSges(6,1) * t10 + 0.2e1 * Ifges(6,5) * t578 + Ifges(7,5) * t589 + 0.2e1 * Ifges(6,6) * t579 + Ifges(7,6) * t588 + 0.2e1 * Ifges(6,3) * t570 + Ifges(7,3) * t573 + t7 / 0.2e1 - mrSges(6,2) * t11 - mrSges(5,3) * t31 + (-t318 / 0.2e1 - t551) * Ifges(5,6) + (t570 - t571) * Ifges(5,2) + (-t133 / 0.2e1 - t572) * Ifges(5,4);
t32 = t101 * t409 - t405 * t111 - t165 * t482 - t179 * t481;
t661 = mrSges(5,2) * t178 - mrSges(5,3) * t32 + 0.2e1 * Ifges(5,1) * t572 + 0.2e1 * Ifges(5,4) * t571 + 0.2e1 * Ifges(5,5) * t551;
t395 = pkin(11) + qJ(4);
t390 = sin(t395);
t392 = cos(t395);
t448 = -mrSges(4,1) * t400 + mrSges(4,2) * t397;
t427 = m(4) * pkin(2) - t448;
t660 = t390 * t422 + t392 * t670 - t427;
t658 = t232 * mrSges(5,2);
t631 = t396 * t668 + t399 * t667;
t630 = t396 * t667 - t399 * t668;
t230 = Ifges(5,4) * t235;
t147 = Ifges(5,1) * t411 + t365 * Ifges(5,5) - t230;
t568 = t147 / 0.2e1;
t97 = t165 * t409 - t405 * t179;
t654 = -t97 * mrSges(5,3) + Ifges(5,1) * t556 + Ifges(5,4) * t559 + Ifges(5,5) * t546 + t568;
t643 = t197 * Ifges(6,5) + t122 * Ifges(7,5) + Ifges(6,6) * t460 + Ifges(7,6) * t640 + t235 * Ifges(6,3) + t231 * Ifges(7,3);
t649 = pkin(5) * t669 + t664 * pkin(10) + t631;
t648 = -pkin(10) * t646 + t630;
t349 = t396 * t408 + t399 * t404;
t148 = t349 * t235;
t335 = t349 * qJD(6);
t624 = t148 + t335;
t432 = t396 * t404 - t399 * t408;
t150 = t432 * t235;
t334 = t432 * qJD(6);
t623 = t150 + t334;
t286 = -t357 * t405 + t359 * t409;
t617 = -qJD(3) * t350 - qJD(4) * t286 - t206 * t409 + t405 * t222;
t644 = m(7) * pkin(5);
t616 = pkin(4) * t469 - t617;
t443 = t389 * mrSges(7,1) + t391 * mrSges(7,2);
t513 = t399 * mrSges(6,2);
t444 = t396 * mrSges(6,1) + t513;
t632 = -m(4) * qJ(3) - mrSges(4,3) - mrSges(5,3);
t641 = -t396 * t644 - t443 - t444 + t632;
t550 = -t642 / 0.2e1;
t553 = t267 / 0.2e1;
t639 = Ifges(4,1) * t553 + Ifges(4,5) * t550;
t638 = m(3) + m(4);
t554 = t266 / 0.2e1;
t636 = pkin(3) * t397;
t388 = pkin(3) * t400 + pkin(2);
t262 = pkin(4) * t348 - qJ(5) * t350 - t388;
t186 = t399 * t262 - t286 * t396;
t504 = t350 * t399;
t160 = pkin(5) * t348 - pkin(10) * t504 + t186;
t187 = t396 * t262 + t399 * t286;
t505 = t350 * t396;
t170 = -pkin(10) * t505 + t187;
t88 = t160 * t408 - t170 * t404;
t635 = qJD(6) * t88 + t404 * t649 + t408 * t648;
t89 = t160 * t404 + t170 * t408;
t634 = -qJD(6) * t89 - t404 * t648 + t408 * t649;
t629 = pkin(5) * t646 + t616;
t356 = t535 * t396;
t358 = t535 * t399;
t283 = -t356 * t408 - t358 * t404;
t506 = t235 * t399;
t156 = pkin(4) * t411 + qJ(5) * t235;
t65 = t399 * t156 - t396 * t97;
t45 = pkin(5) * t411 + pkin(10) * t506 + t65;
t507 = t235 * t396;
t66 = t396 * t156 + t399 * t97;
t53 = pkin(10) * t507 + t66;
t628 = -qJD(5) * t432 + qJD(6) * t283 - t404 * t45 - t408 * t53;
t285 = -t356 * t404 + t358 * t408;
t627 = -qJD(5) * t349 - qJD(6) * t285 + t404 * t53 - t408 * t45;
t163 = t250 * t408 - t251 * t404;
t174 = t334 * t350 + t336 * t349;
t622 = t163 - t174;
t164 = t250 * t404 + t251 * t408;
t173 = -t335 * t350 + t336 * t432;
t621 = t164 - t173;
t124 = -mrSges(6,1) * t460 + mrSges(6,2) * t197;
t529 = mrSges(5,3) * t411;
t205 = mrSges(5,1) * t365 - t529;
t620 = t205 - t124;
t457 = mrSges(3,3) * t469;
t619 = -mrSges(3,1) * t459 + mrSges(4,1) * t414 + t305 * mrSges(4,2) + t457;
t609 = -t143 * t397 + t144 * t400;
t140 = -mrSges(6,2) * t235 + mrSges(6,3) * t460;
t141 = mrSges(6,1) * t235 - mrSges(6,3) * t197;
t608 = t140 * t399 - t141 * t396;
t607 = -t10 * t396 + t11 * t399;
t447 = mrSges(4,1) * t397 + mrSges(4,2) * t400;
t604 = m(6) + m(5) + m(7);
t528 = Ifges(3,4) * t406;
t603 = pkin(1) * (mrSges(3,1) * t406 + mrSges(3,2) * t410) - t406 * (Ifges(3,1) * t410 - t528) / 0.2e1;
t600 = -t396 * (mrSges(6,1) + t644) + mrSges(3,2) - t513 + t632;
t599 = -t32 * mrSges(5,1) + t31 * mrSges(5,2) - Ifges(5,5) * t133 + Ifges(5,6) * t134 - Ifges(5,3) * t318;
t598 = t259 * mrSges(3,1) - t258 * mrSges(3,2) + Ifges(3,5) * t331 + Ifges(3,6) * t642 + Ifges(3,3) * t379;
t596 = mrSges(3,2) + t641;
t595 = mrSges(3,1) - t660;
t544 = t399 / 0.2e1;
t92 = t197 * Ifges(6,4) + Ifges(6,2) * t460 + Ifges(6,6) * t235;
t580 = -t92 / 0.2e1;
t93 = t197 * Ifges(6,1) + Ifges(6,4) * t460 + Ifges(6,5) * t235;
t94 = -pkin(4) * t365 + qJD(5) - t97;
t594 = t396 * t580 + t444 * t94 + t93 * t544 + t658;
t547 = -t365 / 0.2e1;
t561 = -t231 / 0.2e1;
t563 = -t197 / 0.2e1;
t565 = -t460 / 0.2e1;
t575 = -t122 / 0.2e1;
t577 = -t640 / 0.2e1;
t592 = Ifges(6,5) * t563 + Ifges(7,5) * t575 - Ifges(5,2) * t558 - Ifges(5,6) * t547 + Ifges(6,6) * t565 + Ifges(7,6) * t577 + Ifges(6,3) * t559 + Ifges(7,3) * t561 - t653;
t591 = Ifges(7,4) * t589 + Ifges(7,2) * t588 + Ifges(7,6) * t573;
t590 = Ifges(7,1) * t589 + Ifges(7,4) * t588 + Ifges(7,5) * t573;
t587 = Ifges(6,1) * t578 + Ifges(6,4) * t579 + Ifges(6,5) * t570;
t522 = Ifges(7,4) * t122;
t55 = Ifges(7,2) * t640 + Ifges(7,6) * t231 + t522;
t586 = -t55 / 0.2e1;
t585 = t55 / 0.2e1;
t116 = Ifges(7,4) * t640;
t56 = Ifges(7,1) * t122 + Ifges(7,5) * t231 + t116;
t584 = -t56 / 0.2e1;
t583 = t56 / 0.2e1;
t567 = Ifges(4,4) * t554 + t639;
t557 = -t411 / 0.2e1;
t545 = -t397 / 0.2e1;
t543 = t400 / 0.2e1;
t542 = t410 / 0.2e1;
t541 = cos(qJ(1));
t538 = pkin(1) * t410;
t293 = (qJD(2) * t436 - t483) * t398;
t328 = -pkin(8) * t467 + t410 * t478;
t300 = qJD(3) * t401 + t328;
t220 = t400 * t293 - t397 * t300;
t183 = qJD(2) * t424 + t220;
t344 = pkin(8) * t498 + t406 * t539;
t312 = qJ(3) * t401 + t344;
t241 = -t397 * t312 + t400 * t313;
t333 = t397 * t401 + t400 * t500;
t185 = -pkin(3) * t498 - t333 * pkin(9) + t241;
t221 = t397 * t293 + t400 * t300;
t454 = t397 * t398 * t484;
t200 = -pkin(9) * t454 + t221;
t242 = t400 * t312 + t397 * t313;
t495 = t400 * t401;
t332 = -t397 * t500 + t495;
t203 = pkin(9) * t332 + t242;
t63 = t405 * t183 + t185 * t481 + t409 * t200 - t203 * t482;
t59 = (qJ(5) * t485 - qJD(5) * t410) * t398 + t63;
t433 = t409 * t332 - t333 * t405;
t192 = -qJD(2) * t419 + qJD(4) * t433;
t253 = t332 * t405 + t333 * t409;
t193 = qJD(2) * t420 + qJD(4) * t253;
t329 = t344 * qJD(2);
t279 = pkin(3) * t454 + t329;
t84 = t193 * pkin(4) - t192 * qJ(5) - t253 * qJD(5) + t279;
t30 = t396 * t84 + t399 * t59;
t530 = mrSges(5,3) * t235;
t527 = Ifges(4,4) * t397;
t526 = Ifges(4,4) * t400;
t524 = Ifges(6,4) * t396;
t523 = Ifges(6,4) * t399;
t521 = Ifges(4,5) * t305;
t520 = Ifges(4,6) * t406;
t519 = Ifges(4,3) * t406;
t516 = t235 * Ifges(5,6);
t515 = t411 * Ifges(5,5);
t514 = t365 * Ifges(5,3);
t512 = t400 * Ifges(4,6);
t501 = t397 * t410;
t407 = sin(qJ(1));
t499 = t398 * t407;
t496 = t400 * (Ifges(4,1) * t305 - Ifges(4,4) * t414 - Ifges(4,5) * t468);
t493 = t406 * t407;
t492 = t407 * t410;
t118 = t405 * t185 + t409 * t203;
t112 = -qJ(5) * t498 + t118;
t315 = t380 + (-pkin(2) - t538) * t401;
t260 = -t332 * pkin(3) + t315;
t142 = -pkin(4) * t433 - t253 * qJ(5) + t260;
t73 = t399 * t112 + t396 * t142;
t489 = t541 * pkin(1) + pkin(8) * t499;
t472 = t398 * t541;
t471 = t541 * t406;
t470 = t541 * t410;
t12 = -t40 * mrSges(7,1) + t39 * mrSges(7,2);
t461 = pkin(1) * t407 - pkin(8) * t472;
t29 = -t396 * t59 + t399 * t84;
t191 = -t266 * mrSges(4,1) + t267 * mrSges(4,2);
t78 = t134 * mrSges(5,1) + t133 * mrSges(5,2);
t60 = -t109 * mrSges(6,1) + t110 * mrSges(6,2);
t72 = -t112 * t396 + t399 * t142;
t117 = t185 * t409 - t405 * t203;
t339 = t401 * t471 + t492;
t269 = t339 * t392 - t390 * t472;
t456 = mrSges(3,3) * t468;
t113 = pkin(4) * t498 - t117;
t442 = Ifges(4,1) * t400 - t527;
t441 = Ifges(6,1) * t399 - t524;
t440 = -Ifges(4,2) * t397 + t526;
t439 = -Ifges(6,2) * t396 + t523;
t438 = Ifges(3,5) * t410 - Ifges(3,6) * t406;
t437 = Ifges(6,5) * t399 - Ifges(6,6) * t396;
t435 = -t396 * t51 + t399 * t52;
t229 = t399 * t253 - t396 * t498;
t50 = -pkin(5) * t433 - pkin(10) * t229 + t72;
t228 = -t396 * t253 - t399 * t498;
t57 = pkin(10) * t228 + t73;
t17 = -t404 * t57 + t408 * t50;
t18 = t404 * t50 + t408 * t57;
t152 = t228 * t408 - t229 * t404;
t153 = t228 * t404 + t229 * t408;
t64 = t183 * t409 - t185 * t482 - t405 * t200 - t203 * t481;
t429 = t280 * t447;
t268 = t339 * t390 + t392 * t472;
t341 = -t401 * t493 + t470;
t272 = t341 * t390 - t392 * t499;
t310 = t390 * t500 - t401 * t392;
t425 = -g(1) * t272 - g(2) * t268 - g(3) * t310;
t418 = mrSges(3,1) + t427;
t28 = -pkin(4) * t318 + qJDD(5) - t32;
t416 = Ifges(4,4) * t494 - Ifges(4,2) * t501 + t520;
t61 = -pkin(4) * t467 - t64;
t372 = Ifges(3,4) * t468;
t343 = t401 * t538 - t380;
t342 = (-mrSges(3,1) * t410 + mrSges(3,2) * t406) * t398;
t340 = t401 * t492 + t471;
t338 = -t401 * t470 + t493;
t324 = -mrSges(3,2) * t459 + t456;
t311 = t390 * t401 + t392 * t500;
t290 = Ifges(3,1) * t469 + Ifges(3,5) * t459 + t372;
t289 = Ifges(3,6) * qJD(2) + (Ifges(3,6) * t401 + (t410 * Ifges(3,2) + t528) * t398) * qJD(1);
t273 = t341 * t392 + t390 * t499;
t265 = -mrSges(4,1) * t468 - t305 * mrSges(4,3);
t264 = mrSges(4,2) * t468 - mrSges(4,3) * t414;
t255 = t432 * t350;
t254 = t349 * t350;
t247 = pkin(5) * t505 - t613;
t226 = -mrSges(4,1) * t642 - mrSges(4,3) * t267;
t225 = mrSges(4,2) * t642 + mrSges(4,3) * t266;
t217 = Ifges(4,4) * t305 - Ifges(4,2) * t414 - Ifges(4,6) * t468;
t216 = -Ifges(4,6) * t414 - Ifges(4,3) * t468 + t521;
t208 = t273 * t391 + t340 * t389;
t207 = -t273 * t389 + t340 * t391;
t204 = -mrSges(5,2) * t365 - t530;
t177 = t192 * t399 + t396 * t467;
t176 = -t192 * t396 + t399 * t467;
t168 = t267 * Ifges(4,4) + t266 * Ifges(4,2) - Ifges(4,6) * t642;
t157 = mrSges(5,1) * t235 + mrSges(5,2) * t411;
t145 = t514 + t515 - t516;
t115 = -mrSges(5,2) * t318 - mrSges(5,3) * t134;
t114 = mrSges(5,1) * t318 - mrSges(5,3) * t133;
t87 = mrSges(7,1) * t231 - mrSges(7,3) * t122;
t86 = -mrSges(7,2) * t231 + mrSges(7,3) * t640;
t85 = -pkin(5) * t228 + t113;
t79 = -pkin(5) * t507 + t98;
t77 = -pkin(5) * t460 + t94;
t71 = mrSges(6,1) * t134 - mrSges(6,3) * t110;
t70 = -mrSges(6,2) * t134 + mrSges(6,3) * t109;
t69 = -qJD(6) * t153 + t176 * t408 - t177 * t404;
t68 = qJD(6) * t152 + t176 * t404 + t177 * t408;
t67 = -mrSges(7,1) * t640 + mrSges(7,2) * t122;
t48 = -pkin(5) * t176 + t61;
t42 = t110 * Ifges(6,4) + t109 * Ifges(6,2) + t134 * Ifges(6,6);
t23 = -mrSges(7,2) * t132 + mrSges(7,3) * t40;
t22 = mrSges(7,1) * t132 - mrSges(7,3) * t39;
t21 = pkin(10) * t176 + t30;
t20 = -pkin(5) * t109 + t28;
t19 = pkin(5) * t193 - pkin(10) * t177 + t29;
t4 = -qJD(6) * t18 + t19 * t408 - t21 * t404;
t3 = qJD(6) * t17 + t19 * t404 + t21 * t408;
t8 = [(-t208 * mrSges(7,1) - t207 * mrSges(7,2) - t541 * mrSges(2,1) + t407 * mrSges(2,2) - t418 * t341 + t415 * t273 + t600 * t340 + t422 * t272 - t638 * t489 - t604 * (t340 * t536 + t341 * t388 + t499 * t636 + t489)) * g(2) + (Ifges(7,4) * t153 + Ifges(7,2) * t152) * t588 + (Ifges(7,4) * t68 + Ifges(7,2) * t69) * t576 + t619 * t329 + (Ifges(7,1) * t153 + Ifges(7,4) * t152) * t589 + (Ifges(7,1) * t68 + Ifges(7,4) * t69) * t574 + (Ifges(6,4) * t229 + Ifges(6,2) * t228) * t579 + (Ifges(6,4) * t177 + Ifges(6,2) * t176) * t564 + t598 * t401 + ((mrSges(3,1) * t642 - mrSges(3,2) * t331 + (m(3) * t398 * pkin(1) - t342) * qJDD(1)) * pkin(1) + (-mrSges(3,3) * t259 + Ifges(3,1) * t331 + Ifges(3,4) * t642 + Ifges(3,5) * t379) * t406 + (-t266 * Ifges(4,6) - t267 * Ifges(4,5) - t143 * mrSges(4,1) + t379 * Ifges(3,6) + t331 * Ifges(3,4) + t258 * mrSges(3,3) + t144 * mrSges(4,2) - (-Ifges(4,3) - Ifges(3,2)) * t642 + t599) * t410 + ((t290 / 0.2e1 + t305 * t442 / 0.2e1 + t429 + t496 / 0.2e1 + t217 * t545 - t326 * mrSges(3,3) + (Ifges(3,5) / 0.2e1 + t440 * t543) * qJD(2) + (-t210 * t400 - t211 * t397) * mrSges(4,3)) * t410 + (t145 / 0.2e1 + t216 / 0.2e1 - t289 / 0.2e1 + t521 / 0.2e1 + t210 * mrSges(4,1) - t211 * mrSges(4,2) + t515 / 0.2e1 - t516 / 0.2e1 + t97 * mrSges(5,1) - t98 * mrSges(5,2) + t514 / 0.2e1 - t327 * mrSges(3,3) + (-Ifges(3,6) / 0.2e1 + t512 / 0.2e1) * qJD(2)) * t406 + (t401 * t438 / 0.2e1 + (t406 * t416 * t545 + (Ifges(3,4) * t410 - Ifges(3,2) * t406) * t542 - t410 * (Ifges(4,5) * t494 - Ifges(4,6) * t501 + t519) / 0.2e1 - t603) * t398 + t416 * t495 / 0.2e1) * qJD(1)) * qJD(2) + (g(1) * t541 + g(2) * t407) * (-mrSges(3,3) - t447)) * t398 + t344 * (-mrSges(3,2) * t379 + mrSges(3,3) * t642) + (Ifges(6,1) * t229 + Ifges(6,4) * t228) * t578 + (Ifges(6,1) * t177 + Ifges(6,4) * t176) * t562 + (t643 / 0.2e1 + t666) * t193 + (Ifges(4,5) * t333 + Ifges(4,6) * t332) * t550 + (Ifges(4,1) * t333 + Ifges(4,4) * t332) * t553 + (Ifges(4,4) * t333 + Ifges(4,2) * t332) * t554 + (t654 + t658) * t192 + (t1 * t152 - t13 * t68 + t14 * t69 - t153 * t2) * mrSges(7,3) + t333 * t567 + (Ifges(7,5) * t153 + Ifges(7,6) * t152) * t573 + m(3) * (t258 * t344 + t259 * t343 - t326 * t329 + t327 * t328) + (Ifges(7,5) * t68 + Ifges(7,6) * t69) * t560 + (Ifges(6,5) * t177 + Ifges(6,6) * t176) * t558 + (Ifges(6,5) * t229 + Ifges(6,6) * t228) * t570 + m(4) * (t143 * t241 + t144 * t242 + t210 * t220 + t211 * t221 + t239 * t315 + t280 * t329) + m(5) * (t117 * t32 + t118 * t31 + t178 * t260 + t232 * t279 + t63 * t98 + t64 * t97) + m(6) * (t10 * t72 + t11 * t73 + t113 * t28 + t29 * t51 + t30 * t52 + t61 * t94) + m(7) * (t1 * t18 + t13 * t4 + t14 * t3 + t17 * t2 + t20 * t85 + t48 * t77) + (-t10 * t229 + t11 * t228 + t176 * t52 - t177 * t51) * mrSges(6,3) + Ifges(2,3) * qJDD(1) + t343 * (mrSges(3,1) * t379 - mrSges(3,3) * t331) + t332 * t168 / 0.2e1 + t239 * (-mrSges(4,1) * t332 + mrSges(4,2) * t333) + t328 * t324 + t315 * t191 + t279 * t157 + t220 * t265 + t221 * t264 + t260 * t78 + t241 * t226 + t242 * t225 + t28 * (-mrSges(6,1) * t228 + mrSges(6,2) * t229) + t228 * t42 / 0.2e1 + t63 * t204 + t64 * t205 + t176 * t92 / 0.2e1 + t94 * (-mrSges(6,1) * t176 + mrSges(6,2) * t177) + t177 * t93 / 0.2e1 + t20 * (-mrSges(7,1) * t152 + mrSges(7,2) * t153) + t30 * t140 + t29 * t141 + t61 * t124 + t117 * t114 + t118 * t115 + t113 * t60 + t85 * t12 + t3 * t86 + t4 * t87 + t17 * t22 + t18 * t23 + t661 * t253 - t663 * t433 + t48 * t67 + t72 * t71 + t73 * t70 + t77 * (-mrSges(7,1) * t69 + mrSges(7,2) * t68) + t68 * t583 + t69 * t585 + t229 * t587 + t153 * t590 + t152 * t591 + (t407 * mrSges(2,1) + t541 * mrSges(2,2) + t418 * t339 - t670 * t269 + (t443 - t600) * t338 - t422 * t268 + t638 * t461 + t604 * (t338 * t536 + t339 * t388 - t472 * t636 + t461)) * g(1) + (-t143 * t333 + t144 * t332) * mrSges(4,3); -((t216 + t145) * t406 + (-Ifges(3,2) * t469 + t290 + t372 + t496) * t410 + t305 * (Ifges(4,5) * t406 + t410 * t442) + t459 * t438) * t488 / 0.2e1 + (-t604 * (-t340 * t388 + t341 * t536) + t596 * t341 + t595 * t340) * g(1) + t598 + (-t178 * t388 - t232 * t278 + t286 * t31 + t32 * t613 + t617 * t97 + t618 * t98) * m(5) + (t10 * t186 + t11 * t187 - t28 * t613 + t51 * t631 + t52 * t630 + t616 * t94) * m(6) - (t60 - t114) * t613 + (-t604 * (-t338 * t388 + t339 * t536) + t596 * t339 + t595 * t338) * g(2) + (-Ifges(5,1) * t288 + Ifges(5,5) * t469) * t557 - t97 * (mrSges(5,1) * t469 + mrSges(5,3) * t288) + (t232 * t288 + t469 * t98) * mrSges(5,2) + (-Ifges(5,5) * t288 + Ifges(5,3) * t469) * t547 + (-Ifges(5,4) * t288 + Ifges(5,6) * t469) * t558 + t616 * t124 + t617 * t205 + t618 * t204 + (-m(4) * t280 + t457 - t619) * t327 + (mrSges(7,1) * t622 - mrSges(7,2) * t621) * t77 + (-t1 * t254 + t13 * t621 - t14 * t622 + t2 * t255) * mrSges(7,3) + t609 * mrSges(4,3) + (-t210 * t248 - t211 * t249 - pkin(2) * t239 + (-t210 * t397 + t211 * t400) * qJD(3) + t609 * qJ(3)) * m(4) + t666 * t337 + t512 * t550 + (Ifges(4,2) * t400 + t527) * t554 + t288 * t568 - (t437 * t558 + t439 * t564 + t441 * t562 + t594 + t654) * t336 + (Ifges(7,4) * t173 + Ifges(7,2) * t174) * t576 + (Ifges(7,4) * t164 + Ifges(7,2) * t163) * t577 + (t414 * (t410 * t440 + t520) + t406 * t289) * t488 / 0.2e1 + (t342 - t604 * t388 * t498 + (t660 * t410 + (-t536 * t604 + t641) * t406) * t398) * g(3) + (qJD(3) * t400 - t249) * t264 + (t456 - t324) * t326 + (Ifges(7,5) * t173 + Ifges(7,6) * t174) * t560 + (Ifges(7,5) * t164 + Ifges(7,6) * t163) * t561 + (-Ifges(5,4) * t557 + t592 + t665) * t287 - t429 * t468 + (-qJ(3) * t226 - qJD(3) * t265 + t567 + t639) * t397 + (-Ifges(7,4) * t255 - Ifges(7,2) * t254) * t588 + (-Ifges(7,5) * t255 - Ifges(7,6) * t254) * t573 + t20 * (mrSges(7,1) * t254 - mrSges(7,2) * t255) + (-Ifges(7,1) * t255 - Ifges(7,4) * t254) * t589 + (Ifges(6,4) * t251 + Ifges(6,2) * t250) * t565 + t643 * (t337 / 0.2e1 - t287 / 0.2e1) + ((t519 + (Ifges(4,5) * t400 - t397 * Ifges(4,6)) * t410) * t542 + t603) * qJD(1) ^ 2 * t398 ^ 2 + t526 * t553 - t42 * t505 / 0.2e1 - t388 * t78 + (-t210 * (mrSges(4,1) * t406 - mrSges(4,3) * t494) - t211 * (-mrSges(4,2) * t406 - mrSges(4,3) * t501)) * t488 + t286 * t115 - t278 * t157 - t248 * t265 - t94 * (-mrSges(6,1) * t250 + mrSges(6,2) * t251) - t251 * t93 / 0.2e1 + t247 * t12 + t187 * t70 - pkin(2) * t191 + t217 * t455 / 0.2e1 + t186 * t71 + t88 * t22 + t89 * t23 + t239 * t448 + (Ifges(7,1) * t173 + Ifges(7,4) * t174) * t574 + (Ifges(7,1) * t164 + Ifges(7,4) * t163) * t575 + t250 * t580 + (t28 * t444 + t437 * t570 + t439 * t579 + t441 * t578 + t661) * t350 + t400 * qJ(3) * t225 + t663 * t348 + (-t10 * t504 - t11 * t505 + t51 * t664 - t646 * t52) * mrSges(6,3) + t629 * t67 + t630 * t140 + t631 * t141 + t173 * t583 + t164 * t584 + t174 * t585 + t163 * t586 + t504 * t587 - t255 * t590 - t254 * t591 + t634 * t87 + t635 * t86 + (t1 * t89 + t13 * t634 + t14 * t635 + t2 * t88 + t20 * t247 + t629 * t77) * m(7) + t168 * t543 + (Ifges(6,5) * t251 + Ifges(6,6) * t250) * t559 + (Ifges(6,1) * t251 + Ifges(6,4) * t250) * t563; t414 * t264 - t623 * t86 - t624 * t87 - (-t204 - t608) * t235 + t78 + t191 + (-t67 + t620) * t411 + t399 * t71 + t396 * t70 + t349 * t23 - t432 * t22 + t305 * t265 + (t1 * t349 - t13 * t624 - t14 * t623 - t2 * t432 - t411 * t77) * m(7) + (t10 * t399 + t11 * t396 + t235 * t435 - t411 * t94) * m(6) + (t235 * t98 + t411 * t97 + t178) * m(5) + (t210 * t305 + t211 * t414 + t239) * m(4) + (-g(1) * t340 - g(2) * t338 + g(3) * t498) * (m(4) + t604); (-Ifges(5,1) * t557 - Ifges(5,5) * t547 - t437 * t559 - t439 * t565 - t441 * t563 + t594) * t235 + (Ifges(7,5) * t349 - Ifges(7,6) * t432) * t573 + t20 * (mrSges(7,1) * t432 + mrSges(7,2) * t349) + (Ifges(7,4) * t349 - Ifges(7,2) * t432) * t588 + (Ifges(7,1) * t349 - Ifges(7,4) * t432) * t589 + (-t1 * t432 + t13 * t623 - t14 * t624 - t2 * t349) * mrSges(7,3) - t432 * t591 - t599 + (-m(6) * t94 + t529 + t620) * t98 + (-t272 * t670 + t273 * t422) * g(1) + (-t268 * t670 + t269 * t422) * g(2) + (-t310 * t670 + t311 * t422) * g(3) + (-t506 * t51 - t507 * t52 + t607) * mrSges(6,3) + (-pkin(4) * t28 + qJ(5) * t607 + t435 * qJD(5) - t51 * t65 - t52 * t66) * m(6) + t608 * qJD(5) + t146 * t556 + (-t396 * t71 + t399 * t70) * qJ(5) + (Ifges(7,5) * t150 + Ifges(7,6) * t148) * t561 + (Ifges(7,4) * t150 + Ifges(7,2) * t148) * t577 + (-Ifges(7,1) * t334 - Ifges(7,4) * t335) * t574 + (-Ifges(7,4) * t334 - Ifges(7,2) * t335) * t576 + (Ifges(6,5) * t396 + Ifges(6,6) * t399) * t570 + t592 * t411 + (-Ifges(7,5) * t334 - Ifges(7,6) * t335) * t560 + (-t525 + t643) * t557 + (Ifges(7,1) * t150 + Ifges(7,4) * t148) * t575 - t387 * t12 + t283 * t22 + t285 * t23 + (t147 - t230) * t558 - t66 * t140 - t65 * t141 + t28 * t445 + (-t530 - t204) * t97 + (Ifges(6,1) * t396 + t523) * t578 + (Ifges(6,2) * t399 + t524) * t579 + (mrSges(7,1) * t624 - mrSges(7,2) * t623) * t77 - pkin(4) * t60 + t627 * t87 + t628 * t86 + (t1 * t285 + t13 * t627 + t14 * t628 + t2 * t283 - t20 * t387 - t77 * t79) * m(7) - t79 * t67 - t334 * t583 + t150 * t584 - t335 * t585 + t148 * t586 + t396 * t587 + t349 * t590 + t42 * t544; t122 * t87 - t640 * t86 - t460 * t140 + t197 * t141 + t12 + t60 + (t122 * t13 - t14 * t640 + t20 + t425) * m(7) + (t197 * t51 - t460 * t52 + t28 + t425) * m(6); -t77 * (mrSges(7,1) * t122 + mrSges(7,2) * t640) + (Ifges(7,1) * t640 - t522) * t575 + t55 * t574 + (Ifges(7,5) * t640 - Ifges(7,6) * t122) * t561 - t13 * t86 + t14 * t87 - g(1) * (mrSges(7,1) * t207 - mrSges(7,2) * t208) - g(2) * ((-t269 * t389 + t338 * t391) * mrSges(7,1) + (-t269 * t391 - t338 * t389) * mrSges(7,2)) - g(3) * ((-t311 * t389 - t391 * t498) * mrSges(7,1) + (-t311 * t391 + t389 * t498) * mrSges(7,2)) + (t122 * t14 + t13 * t640) * mrSges(7,3) + t7 + (-Ifges(7,2) * t122 + t116 + t56) * t577 + t597;];
tau  = t8;
