% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PPRRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:35
% EndTime: 2019-03-08 19:03:52
% DurationCPUTime: 11.09s
% Computational Cost: add. (23109->693), mult. (62332->990), div. (0->0), fcn. (72665->14), ass. (0->368)
t423 = sin(qJ(3));
t570 = sin(pkin(6));
t501 = sin(pkin(13)) * t570;
t631 = cos(qJ(3));
t419 = sin(pkin(7));
t488 = cos(pkin(13)) * t570;
t571 = cos(pkin(7));
t572 = cos(pkin(6));
t683 = t572 * t419 + t571 * t488;
t277 = t423 * t501 - t631 * t683;
t278 = t423 * t683 + t631 * t501;
t425 = cos(qJ(5));
t421 = sin(qJ(5));
t426 = cos(qJ(4));
t547 = t421 * t426;
t162 = t277 * t547 + t278 * t425;
t542 = t425 * t426;
t163 = -t277 * t542 + t278 * t421;
t420 = sin(qJ(6));
t424 = cos(qJ(6));
t107 = t162 * t424 - t163 * t420;
t108 = t162 * t420 + t163 * t424;
t668 = m(7) * pkin(5);
t712 = -mrSges(7,2) / 0.2e1;
t713 = mrSges(7,1) / 0.2e1;
t687 = t107 * t713 + t108 * t712;
t721 = -t162 * mrSges(6,1) / 0.2e1 + t163 * mrSges(6,2) / 0.2e1 - (t107 * t424 + t108 * t420) * t668 / 0.2e1 - t687;
t414 = t421 ^ 2;
t416 = t425 ^ 2;
t684 = t416 + t414;
t502 = pkin(10) * t684;
t422 = sin(qJ(4));
t443 = -t419 * t488 + t571 * t572;
t212 = t278 * t426 + t422 * t443;
t134 = -t212 * t421 + t277 * t425;
t135 = t212 * t425 + t277 * t421;
t100 = t134 * t420 + t135 * t424;
t498 = t424 * t134 - t135 * t420;
t28 = -t100 * mrSges(7,1) - t498 * mrSges(7,2);
t720 = t28 * qJD(6);
t412 = Ifges(6,4) * t425;
t485 = -t421 * Ifges(6,1) - t412;
t632 = t425 / 0.2e1;
t719 = t485 * t632;
t552 = t419 * t423;
t360 = t422 * t571 + t426 * t552;
t518 = t419 * t631;
t291 = -t421 * t360 - t425 * t518;
t292 = t360 * t425 - t421 * t518;
t194 = t291 * t420 + t292 * t424;
t497 = t424 * t291 - t292 * t420;
t61 = -t194 * mrSges(7,1) - t497 * mrSges(7,2);
t718 = t61 * qJD(6);
t544 = t424 * t425;
t470 = t420 * t421 - t544;
t471 = t420 * t425 + t424 * t421;
t286 = -Ifges(7,5) * t470 - Ifges(7,6) * t471;
t664 = -pkin(11) - pkin(10);
t397 = t664 * t421;
t398 = t664 * t425;
t308 = t397 * t420 - t398 * t424;
t496 = t424 * t397 + t398 * t420;
t75 = -t308 * mrSges(7,1) - t496 * mrSges(7,2) + t286;
t717 = t75 * qJD(6);
t386 = -pkin(4) * t426 - t422 * pkin(10) - pkin(3);
t369 = t425 * t386;
t546 = t422 * t425;
t490 = -pkin(11) * t546 + t369;
t276 = (-pkin(9) * t421 - pkin(5)) * t426 + t490;
t331 = pkin(9) * t542 + t421 * t386;
t548 = t421 * t422;
t298 = -pkin(11) * t548 + t331;
t559 = t298 * t420;
t181 = t276 * t424 - t559;
t536 = pkin(9) * t547;
t297 = t490 - t536;
t199 = t297 * t424 - t559;
t716 = -t181 + t199;
t348 = t471 * t422;
t347 = t420 * t548 - t422 * t544;
t617 = Ifges(7,4) * t347;
t240 = -Ifges(7,2) * t348 - t426 * Ifges(7,6) - t617;
t258 = -Ifges(7,1) * t348 + t617;
t616 = Ifges(7,4) * t471;
t288 = -Ifges(7,2) * t470 + t616;
t289 = -Ifges(7,1) * t470 - t616;
t630 = pkin(5) * t421;
t519 = pkin(9) + t630;
t381 = t519 * t422;
t629 = pkin(5) * t425;
t409 = -pkin(4) - t629;
t284 = mrSges(7,1) * t471 - mrSges(7,2) * t470;
t660 = t284 / 0.2e1;
t254 = -mrSges(7,1) * t347 - mrSges(7,2) * t348;
t710 = t254 / 0.2e1;
t715 = (t288 / 0.4e1 - t289 / 0.4e1) * t347 - t426 * t286 / 0.4e1 + t409 * t710 + t381 * t660 + (-t240 / 0.4e1 + t258 / 0.4e1) * t471;
t714 = -mrSges(7,1) / 0.2e1;
t711 = mrSges(7,2) / 0.2e1;
t635 = t422 / 0.2e1;
t585 = t347 * mrSges(7,3);
t312 = -mrSges(7,1) * t426 + t585;
t655 = -t312 / 0.2e1;
t706 = t194 * t655;
t705 = t308 * t655;
t669 = m(7) / 0.2e1;
t704 = t630 * t669;
t558 = t298 * t424;
t182 = t276 * t420 + t558;
t198 = -t297 * t420 - t558;
t703 = t182 + t198;
t399 = t422 * pkin(4) - pkin(10) * t426;
t333 = pkin(9) * t548 + t425 * t399;
t279 = t422 * pkin(5) - pkin(11) * t542 + t333;
t334 = -pkin(9) * t546 + t421 * t399;
t302 = -pkin(11) * t547 + t334;
t189 = t279 * t424 - t302 * t420;
t192 = t279 * t420 + t302 * t424;
t359 = t422 * t552 - t426 * t571;
t238 = t471 * t359;
t239 = t470 * t359;
t382 = t519 * t426;
t330 = t369 - t536;
t472 = t330 * t421 - t331 * t425;
t457 = pkin(9) * t426 + t472;
t255 = mrSges(7,1) * t348 - mrSges(7,2) * t347;
t577 = t425 * mrSges(6,2);
t581 = t421 * mrSges(6,1);
t390 = t577 + t581;
t363 = t390 * t422;
t503 = t363 / 0.2e1 + t255 / 0.2e1;
t627 = pkin(9) * t422;
t380 = t422 * mrSges(6,1) - mrSges(6,3) * t542;
t640 = t380 / 0.2e1;
t378 = -t422 * mrSges(6,2) - mrSges(6,3) * t547;
t641 = t378 / 0.2e1;
t350 = t470 * t426;
t313 = mrSges(7,1) * t422 + mrSges(7,3) * t350;
t653 = t313 / 0.2e1;
t654 = t312 / 0.2e1;
t349 = t471 * t426;
t311 = -mrSges(7,2) * t422 - mrSges(7,3) * t349;
t656 = t311 / 0.2e1;
t584 = t348 * mrSges(7,3);
t310 = mrSges(7,2) * t426 - t584;
t657 = t310 / 0.2e1;
t671 = m(6) / 0.2e1;
t576 = t426 * mrSges(5,2);
t391 = t422 * mrSges(5,1) + t576;
t688 = -t391 * t518 / 0.2e1;
t701 = (t333 * t291 + t334 * t292 + t359 * t457 + t360 * t627) * t671 + (t181 * t238 + t182 * t239 + t189 * t497 + t192 * t194 + t359 * t382 + t360 * t381) * t669 + t497 * t653 + t194 * t656 + t238 * t654 + t239 * t657 + t291 * t640 + t292 * t641 + t688 + t360 * t503;
t467 = Ifges(7,5) * t350 / 0.2e1 + Ifges(7,6) * t349 / 0.2e1;
t481 = Ifges(7,3) * t635 + t189 * t713 + t192 * t712 - t467;
t699 = t577 / 0.2e1 + t581 / 0.2e1;
t697 = t181 / 0.2e1;
t582 = t471 * mrSges(7,3);
t578 = t422 * mrSges(5,2);
t691 = -t426 * mrSges(5,1) - mrSges(4,1) + t578;
t612 = Ifges(6,6) * t421;
t615 = Ifges(6,5) * t425;
t468 = t615 / 0.2e1 - t612 / 0.2e1;
t690 = Ifges(5,4) - t468;
t285 = mrSges(7,1) * t470 + mrSges(7,2) * t471;
t680 = -mrSges(6,1) * t425 + mrSges(6,2) * t421;
t686 = t285 + t680;
t685 = -Ifges(6,2) * t421 + t412;
t681 = -t333 * t421 + t334 * t425;
t364 = t390 * t426;
t535 = mrSges(6,3) * t548;
t377 = mrSges(6,2) * t426 - t535;
t379 = -mrSges(6,1) * t426 - mrSges(6,3) * t546;
t633 = -t425 / 0.2e1;
t636 = t421 / 0.2e1;
t679 = t377 * t633 + t379 * t636 + t364 / 0.2e1;
t639 = t390 / 0.2e1;
t678 = t639 + t699;
t465 = t238 * t714 + t239 * t711;
t516 = t426 * t631;
t326 = (-t421 * t516 + t423 * t425) * t419;
t327 = (t421 * t423 + t425 * t516) * t419;
t222 = t326 * t424 - t327 * t420;
t223 = t326 * t420 + t327 * t424;
t540 = t222 * t713 + t223 * t712;
t211 = t278 * t422 - t426 * t443;
t117 = t471 * t211;
t118 = t470 * t211;
t466 = t117 * t714 + t118 * t711;
t677 = m(7) * t381 + t255;
t663 = t100 / 0.2e1;
t665 = -t498 / 0.2e1;
t676 = t498 * t657 + (t347 * t663 - t348 * t665) * mrSges(7,3);
t673 = 2 * qJD(4);
t672 = m(5) / 0.2e1;
t670 = -m(7) / 0.2e1;
t667 = mrSges(6,1) / 0.2e1;
t666 = -mrSges(6,2) / 0.2e1;
t662 = -t182 / 0.2e1;
t659 = t285 / 0.2e1;
t658 = -t310 / 0.2e1;
t652 = -t347 / 0.2e1;
t651 = -t348 / 0.2e1;
t649 = -t349 / 0.2e1;
t648 = -t350 / 0.2e1;
t647 = t359 / 0.2e1;
t362 = t680 * t422;
t646 = -t362 / 0.2e1;
t642 = t471 / 0.2e1;
t638 = t391 / 0.2e1;
t637 = -t421 / 0.2e1;
t634 = t424 / 0.2e1;
t417 = t426 ^ 2;
t628 = pkin(9) * t417;
t626 = pkin(10) * t377;
t625 = pkin(10) * t379;
t619 = mrSges(6,3) * t422;
t618 = Ifges(6,4) * t421;
t609 = pkin(5) * qJD(5);
t603 = t181 * mrSges(7,2);
t602 = t182 * mrSges(7,1);
t595 = t198 * mrSges(7,1);
t594 = t199 * mrSges(7,2);
t583 = t470 * mrSges(7,3);
t575 = t426 * Ifges(6,5);
t574 = t426 * Ifges(6,6);
t573 = t680 - mrSges(5,1);
t569 = mrSges(7,3) * qJD(4);
t568 = t134 * t421;
t567 = t135 * t425;
t566 = t162 * t421;
t565 = t163 * t425;
t564 = t189 * t424;
t563 = t192 * t420;
t562 = t277 * t422;
t561 = t291 * t421;
t560 = t292 * t425;
t554 = t381 * t421;
t551 = t420 * t347;
t550 = t421 * t255;
t343 = t422 * t685 - t574;
t549 = t421 * t343;
t545 = t424 * t348;
t486 = Ifges(6,1) * t425 - t618;
t462 = t486 * t422;
t345 = t462 - t575;
t543 = t425 * t345;
t539 = -Ifges(7,5) * t348 + Ifges(7,6) * t347;
t537 = pkin(5) * t546;
t534 = Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1;
t533 = t665 + t498 / 0.2e1;
t532 = -t619 / 0.2e1;
t529 = t585 / 0.2e1;
t528 = t584 / 0.2e1;
t527 = -t583 / 0.2e1;
t526 = -t582 / 0.2e1;
t524 = mrSges(6,3) * t637;
t522 = mrSges(6,3) * t632;
t521 = t663 - t100 / 0.2e1;
t520 = t194 * t529 + t254 * t647 + t497 * t528;
t517 = t422 * t631;
t515 = t211 * t710;
t514 = t211 * t660;
t512 = t284 * t647;
t508 = t662 - t198 / 0.2e1;
t505 = -t199 / 0.2e1 + t697;
t504 = t362 / 0.2e1 - t254 / 0.2e1;
t415 = t422 ^ 2;
t494 = t415 * t518;
t493 = t419 * t517;
t492 = t359 * t517;
t489 = qJD(4) * (t285 + t573);
t392 = t425 * Ifges(6,2) + t618;
t483 = -t612 + t615;
t482 = Ifges(6,5) * t421 + Ifges(6,6) * t425;
t17 = t514 - t466;
t46 = t512 - t465;
t452 = (t211 * t518 - t277 * t359) * t422;
t14 = (t223 * t100 + t107 * t497 + t194 * t108 + t222 * t498 + t452) * t669 + (t326 * t134 + t327 * t135 + t291 * t162 + t292 * t163 + t452) * t671 + ((-t359 * t422 - t360 * t426) * t277 + (t211 * t517 + t212 * t516 + t277 * t423 - t278 * t631) * t419) * t672;
t142 = t211 * t562;
t20 = m(7) * (t100 * t108 + t107 * t498 - t142) + m(6) * (t134 * t162 + t135 * t163 - t142) + m(5) * (-t211 * t422 - t212 * t426 + t278) * t277;
t480 = t20 * qJD(1) + t14 * qJD(2);
t309 = t419 * t492;
t59 = m(7) * (t194 * t223 + t222 * t497 + t309) + m(6) * (t291 * t326 + t292 * t327 + t309) + m(5) * (t360 * t516 - t423 * t518 + t492) * t419;
t479 = t14 * qJD(1) + t59 * qJD(2);
t473 = -t560 + t561;
t474 = -t567 + t568;
t15 = (t100 * t239 + t117 * t497 + t118 * t194 + t211 * t360 + t212 * t359 + t238 * t498) * t669 + ((t212 + t474) * t359 + (t360 + t473) * t211) * t671;
t120 = t211 * t212;
t25 = m(7) * (t100 * t118 + t117 * t498 + t120) + m(6) * (t211 * t474 + t120);
t478 = t25 * qJD(1) + t15 * qJD(2);
t274 = t359 * t360;
t58 = m(7) * (t194 * t239 + t238 * t497 + t274) + m(6) * (t359 * t473 + t274);
t477 = t15 * qJD(1) + t58 * qJD(2);
t464 = t313 * t634 + t420 * t656;
t463 = t392 * t637 - t719;
t461 = qJD(4) * (-mrSges(6,3) * t684 + mrSges(5,2));
t460 = (t117 * t424 + t118 * t420) * t668;
t459 = (t238 * t424 + t239 * t420) * t668;
t256 = mrSges(7,1) * t349 - mrSges(7,2) * t350;
t428 = (-pkin(4) * t493 + (-t326 * t421 + t327 * t425) * pkin(10)) * t671 + (t222 * t496 + t308 * t223 + t409 * t493) * t669 + t222 * t526 + t223 * t527 + t326 * t524 + t327 * t522 + t688 + t686 * t493 / 0.2e1;
t13 = t256 * t647 + t359 * t679 - t428 + t701;
t241 = -Ifges(7,4) * t350 - Ifges(7,2) * t349 + t422 * Ifges(7,6);
t338 = Ifges(7,4) * t348;
t242 = -Ifges(7,1) * t347 - t426 * Ifges(7,5) - t338;
t243 = -Ifges(7,1) * t350 - Ifges(7,4) * t349 + t422 * Ifges(7,5);
t344 = Ifges(6,6) * t422 + t426 * t685;
t346 = Ifges(6,5) * t422 + t426 * t486;
t21 = t381 * t256 + t382 * t255 - pkin(3) * t391 + t334 * t377 + t331 * t378 + t333 * t379 + t330 * t380 + t240 * t649 + t242 * t648 + t243 * t652 + t241 * t651 + t181 * t313 + t192 * t310 + t182 * t311 + t189 * t312 + m(6) * (t330 * t333 + t331 * t334) + m(7) * (t181 * t189 + t182 * t192 + t381 * t382) + (t543 / 0.2e1 - t549 / 0.2e1 + pkin(9) * t363 + t467 + t690 * t426) * t426 + (Ifges(7,5) * t652 + Ifges(7,6) * t651 + t346 * t632 + t344 * t637 + pkin(9) * t364 - t690 * t422 + (m(6) * pkin(9) ^ 2 + Ifges(5,1) - Ifges(5,2) - Ifges(6,3) - Ifges(7,3)) * t426) * t422;
t448 = t256 / 0.2e1 + t679;
t427 = t503 * t212 + t448 * t211 + (t333 * t134 + t334 * t135 + t211 * t457 + t212 * t627) * t671 + (t100 * t192 + t117 * t181 + t118 * t182 + t189 * t498 + t211 * t382 + t212 * t381) * t669 + t100 * t656 + t117 * t654 + t118 * t657 + t134 * t640 + t135 * t641 + t498 * t653;
t435 = -m(6) * (pkin(4) * t562 + (t565 - t566) * pkin(10)) / 0.2e1 + (t107 * t496 + t108 * t308 - t409 * t562) * t670;
t3 = t427 + (t638 - t576 / 0.2e1 + (t659 - mrSges(5,1) / 0.2e1 + t680 / 0.2e1) * t422) * t277 + (t107 * t642 + t108 * t470 / 0.2e1) * mrSges(7,3) + (-t565 / 0.2e1 + t566 / 0.2e1) * mrSges(6,3) + t435;
t455 = t3 * qJD(1) + t13 * qJD(2) + t21 * qJD(3);
t431 = (t716 * t194 + t359 * t537) * t670 - t706 - t291 * t377 / 0.2e1 + t292 * t379 / 0.2e1 + (t703 * t670 + t658) * t497;
t438 = t326 * t667 + t327 * t666 + (t222 * t424 + t223 * t420) * t668 / 0.2e1 + t540;
t18 = t504 * t359 + (t194 * t652 + t497 * t651) * mrSges(7,3) + (t560 / 0.2e1 - t561 / 0.2e1) * t619 + t431 + t438;
t257 = Ifges(7,2) * t347 - t338;
t436 = t381 * t254 - (t242 / 0.2e1 + t257 / 0.2e1) * t348 + (-t258 / 0.2e1 + t240 / 0.2e1 + t182 * mrSges(7,3)) * t347 + t181 * t584 - t426 * t539 / 0.2e1;
t26 = t199 * t310 + t198 * t312 + m(7) * (t181 * t198 + t182 * t199) + t330 * t377 - t331 * t379 + (t426 * t482 / 0.2e1 - pkin(9) * t362 + t345 * t637 + t343 * t633 + (t392 * t636 + t719) * t422 + t677 * t629 + t472 * mrSges(6,3)) * t422 + t436;
t430 = (-t567 / 0.2e1 + t568 / 0.2e1) * t619 + (t211 * t537 + t703 * t498) * t669 + t134 * t377 / 0.2e1 - t135 * t379 / 0.2e1 + t676 + (t716 * t669 - t654) * t100;
t5 = -t211 * t504 + t430 + t721;
t454 = t5 * qJD(1) - t18 * qJD(2) + t26 * qJD(3);
t449 = t497 * t657 + t520 + t706;
t36 = t449 - t540;
t41 = t181 * t310 - t182 * t312 + t436;
t437 = t100 * t655 + t515 + t676;
t8 = t437 - t687;
t453 = t8 * qJD(1) + t36 * qJD(2) + t41 * qJD(3);
t451 = pkin(4) * t646 + t496 * t658 - t705;
t450 = t660 + t639 - t699;
t367 = Ifges(7,4) * t470;
t287 = -Ifges(7,2) * t471 - t367;
t290 = Ifges(7,1) * t471 - t367;
t444 = t333 * t667 + t334 * t666 + t481;
t445 = t716 * t308 + t703 * t496;
t10 = t444 + t445 * t670 + (-t550 / 0.2e1 + 0.2e1 * (t564 / 0.4e1 + t563 / 0.4e1 - t554 / 0.4e1) * m(7) + t464) * pkin(5) + (0.3e1 / 0.4e1 * t575 - t345 / 0.4e1 + t625 / 0.2e1 + (0.3e1 / 0.4e1 * t618 + t392 / 0.4e1 + t534 * t425 + (-t285 / 0.2e1 + t409 * t670) * pkin(5)) * t422) * t425 + (Ifges(6,3) / 0.2e1 - pkin(9) * t390 / 0.2e1 + (t416 / 0.2e1 + t414 / 0.2e1) * pkin(10) * mrSges(6,3) + (-t485 / 0.4e1 + t412 / 0.4e1 - t534 * t421) * t421) * t422 + (t308 * t652 - t470 * t505 - t471 * t508 + t496 * t651) * mrSges(7,3) + (t343 / 0.4e1 - 0.3e1 / 0.4e1 * t574 + t626 / 0.2e1) * t421 - (-t290 / 0.4e1 - t287 / 0.4e1) * t348 - (-t257 / 0.4e1 - t242 / 0.4e1) * t470 + t451 - t715;
t432 = t359 * t704;
t30 = t450 * t359 - t459 / 0.2e1 + t432 + t465;
t64 = t409 * t284 - (-t289 / 0.2e1 + t288 / 0.2e1) * t471 - (t290 / 0.2e1 + t287 / 0.2e1) * t470;
t47 = -pkin(4) * t390 + t486 * t636 + t632 * t685 + t463 + t64 + (m(7) * t409 + t285) * t630;
t433 = (-t470 * t533 - t471 * t521) * mrSges(7,3) + t211 * t704;
t7 = t450 * t211 - t460 / 0.2e1 + t433 + t466;
t447 = t7 * qJD(1) + t30 * qJD(2) - t10 * qJD(3) + t47 * qJD(4);
t16 = t514 + t466;
t441 = t182 * t526 + t308 * t529 + t496 * t528 - (t287 + t290) * t348 / 0.4e1 - (t242 + t257) * t470 / 0.4e1 + t715;
t434 = t496 * t657 - t582 * t662 + t441 + t705;
t24 = t434 - t481;
t45 = t512 + t465;
t446 = t16 * qJD(1) + t45 * qJD(2) + t24 * qJD(3) + t64 * qJD(4);
t27 = -mrSges(7,1) * t521 + mrSges(7,2) * t533;
t385 = (mrSges(7,1) * t420 + mrSges(7,2) * t424) * pkin(5);
t439 = (t310 * t634 + t420 * t655 + (t551 / 0.2e1 + t545 / 0.2e1) * mrSges(7,3)) * pkin(5);
t54 = mrSges(7,1) * t508 - mrSges(7,2) * t505 + t439;
t442 = -qJD(1) * t27 - qJD(3) * t54 + qJD(5) * t385;
t387 = pkin(9) * t494;
t370 = t385 * qJD(6);
t265 = t415 * pkin(9) * t277;
t48 = -t603 / 0.2e1 - t602 / 0.2e1 - t594 / 0.2e1 + t595 / 0.2e1 + t439 + t539;
t37 = t449 + t540;
t31 = t459 / 0.2e1 + t432 + t46 + t678 * t359;
t23 = t434 + t481;
t19 = t359 * t646 + t532 * t560 + t291 * t535 / 0.2e1 - t431 + t438 + t520;
t12 = t359 * t448 + t428 + t701;
t11 = t444 + t441 + t625 * t633 + t583 * t697 - t392 * t546 / 0.2e1 + Ifges(6,3) * t635 + t626 * t637 + t627 * t639 - t549 / 0.4e1 - t451 + t537 * t659 + t198 * t526 + t199 * t527 + t445 * t669 + t543 / 0.4e1 + t425 * t462 / 0.4e1 + t532 * t502 + (-t483 / 0.4e1 + t468) * t426 + (t464 + t550 / 0.2e1 + (t409 * t546 + t554 + t563 + t564) * t669) * pkin(5) + (t485 / 0.2e1 - t685 / 0.4e1) * t548;
t9 = t437 + t687;
t6 = t460 / 0.2e1 + t433 + t17 + t678 * t211;
t4 = t211 * t646 + t430 + t515 - t721;
t2 = t427 + t163 * t522 + t162 * t524 + t107 * t526 + t108 * t527 - t435 + (t638 + t576 / 0.2e1) * t277 + (mrSges(5,1) / 0.2e1 - t686 / 0.2e1) * t562;
t1 = t14 * qJD(3) + t15 * qJD(4);
t22 = [t20 * qJD(3) + t25 * qJD(4), t1 (t107 * t312 + t108 * t310 + t162 * t379 + t163 * t377 + t691 * t278) * qJD(3) + t2 * qJD(4) + t4 * qJD(5) + t9 * qJD(6) + (mrSges(4,2) + (-t255 - t363) * t422 + (-t415 - t417) * mrSges(5,3)) * qJD(3) * t277 + 0.2e1 * ((t107 * t181 + t108 * t182 - t381 * t562) * t669 + (t162 * t330 + t163 * t331 - t265) * t671 + (-pkin(3) * t278 - t277 * t628 - t265) * t672) * qJD(3) + t480, t2 * qJD(3) + t6 * qJD(5) + t17 * qJD(6) + (-t117 * t471 - t118 * t470) * t569 + t212 * t489 + t211 * t461 + ((t117 * t496 + t118 * t308 + t212 * t409) * t669 + (-pkin(4) * t212 - t211 * t502) * t671) * t673 + t478, t4 * qJD(3) + t6 * qJD(4) + (-t135 * mrSges(6,1) - t134 * mrSges(6,2) + (-t100 * t424 + t420 * t498) * t668 + t28) * qJD(5) + t720, t9 * qJD(3) + t17 * qJD(4) + t28 * qJD(5) + t720; t1, t59 * qJD(3) + t58 * qJD(4) (t223 * t310 + t222 * t312 + m(7) * (t181 * t222 + t182 * t223) + m(6) * (t326 * t330 + t327 * t331 + t387) + t327 * t377 + t326 * t379 + m(5) * (t387 + (-pkin(3) * t423 + t628 * t631) * t419) - mrSges(4,2) * t518 + t691 * t552 + (t363 + t677) * t493 + (t417 * t518 + t494) * mrSges(5,3)) * qJD(3) + t12 * qJD(4) + t19 * qJD(5) + t37 * qJD(6) + t479, t12 * qJD(3) + t31 * qJD(5) + t46 * qJD(6) + (-t238 * t471 - t239 * t470) * t569 + t360 * t489 + t359 * t461 + ((t238 * t496 + t239 * t308 + t360 * t409) * t669 + (-pkin(4) * t360 - t359 * t502) * t671) * t673 + t477, t19 * qJD(3) + t31 * qJD(4) + (-t292 * mrSges(6,1) - t291 * mrSges(6,2) + (-t194 * t424 + t420 * t497) * t668 + t61) * qJD(5) + t718, t37 * qJD(3) + t46 * qJD(4) + t61 * qJD(5) + t718; qJD(4) * t3 + qJD(5) * t5 + qJD(6) * t8 - t480, qJD(4) * t13 - qJD(5) * t18 + qJD(6) * t36 - t479, qJD(4) * t21 + qJD(5) * t26 + qJD(6) * t41, t11 * qJD(5) + t23 * qJD(6) + t455 + (t344 * t632 + t346 * t636 - Ifges(5,6) * t422 + t409 * t256 + t382 * t285 + t243 * t642 - t470 * t241 / 0.2e1 - pkin(4) * t364 + t288 * t649 + t290 * t648 + t496 * t313 + t308 * t311 + pkin(9) * t578 + m(7) * (t189 * t496 + t192 * t308 + t382 * t409) - t189 * t582 - t192 * t583 + (m(6) * t681 + t425 * t378 - t421 * t380) * pkin(10) + (Ifges(5,5) + (-m(6) * pkin(4) + t573) * pkin(9) + t463) * t426 + (Ifges(7,5) * t471 - Ifges(7,6) * t470 + t482) * t635 + t681 * mrSges(6,3)) * qJD(4), t11 * qJD(4) + (-t331 * mrSges(6,1) - t330 * mrSges(6,2) - Ifges(6,5) * t548 - Ifges(6,6) * t546 + t539 - t594 + t595) * qJD(5) + t48 * qJD(6) + (m(7) * (t198 * t424 + t199 * t420) + (t545 + t551) * mrSges(7,3)) * t609 + t454, t23 * qJD(4) + t48 * qJD(5) + (t539 - t602 - t603) * qJD(6) + t453; -qJD(3) * t3 + qJD(5) * t7 + qJD(6) * t16 - t478, -qJD(3) * t13 + qJD(5) * t30 + qJD(6) * t45 - t477, -qJD(5) * t10 + qJD(6) * t24 - t455, qJD(5) * t47 + qJD(6) * t64 (pkin(10) * t680 + t483 + t75) * qJD(5) + t717 + (m(7) * (-t308 * t424 + t420 * t496) + (-t420 * t471 + t424 * t470) * mrSges(7,3)) * t609 + t447, t75 * qJD(5) + t446 + t717; -qJD(3) * t5 - qJD(4) * t7 + qJD(6) * t27, qJD(3) * t18 - qJD(4) * t30, qJD(4) * t10 + qJD(6) * t54 - t454, -t447, -t370, -t370 - t442; -t8 * qJD(3) - t16 * qJD(4) - t27 * qJD(5), -t36 * qJD(3) - t45 * qJD(4), -qJD(4) * t24 - qJD(5) * t54 - t453, -t446, t442, 0;];
Cq  = t22;
