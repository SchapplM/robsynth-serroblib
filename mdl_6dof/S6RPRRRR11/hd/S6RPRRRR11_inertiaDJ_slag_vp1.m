% Calculate time derivative of joint inertia matrix for
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
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR11_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR11_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:24
% EndTime: 2019-03-09 07:39:33
% DurationCPUTime: 37.96s
% Computational Cost: add. (276240->1364), mult. (774269->1728), div. (0->0), fcn. (934857->16), ass. (0->551)
t549 = cos(qJ(1));
t547 = sin(qJ(1));
t703 = cos(pkin(13));
t705 = cos(pkin(6));
t601 = t705 * t703;
t701 = sin(pkin(13));
t570 = t547 * t701 - t549 * t601;
t543 = sin(pkin(6));
t704 = cos(pkin(7));
t625 = t543 * t704;
t702 = sin(pkin(7));
t507 = -t549 * t625 + t570 * t702;
t599 = t705 * t701;
t522 = -t547 * t599 + t549 * t703;
t518 = t522 * qJD(1);
t546 = sin(qJ(3));
t571 = t547 * t601 + t549 * t701;
t563 = t571 * t704;
t714 = cos(qJ(3));
t557 = t714 * t563;
t565 = t570 * t704;
t624 = t543 * t702;
t591 = t714 * t624;
t586 = qJD(1) * t591;
t602 = t549 * t624;
t569 = -t547 * t703 - t549 * t599;
t631 = t569 * t714;
t649 = qJD(3) * t546;
t439 = t518 * t546 + qJD(1) * t557 - t602 * t649 - t547 * t586 + (-t546 * t565 - t631) * qJD(3);
t544 = sin(qJ(5));
t492 = -t631 + (-t565 - t602) * t546;
t545 = sin(qJ(4));
t713 = cos(qJ(4));
t465 = t492 * t713 + t507 * t545;
t558 = t714 * t565;
t491 = -t546 * t569 + t549 * t591 + t558;
t548 = cos(qJ(5));
t593 = t465 * t544 - t491 * t548;
t732 = qJD(5) * t593 - t439 * t544;
t738 = pkin(5) * t732;
t494 = t522 * t714 + (t547 * t624 - t563) * t546;
t562 = t571 * t702;
t604 = t547 * t625;
t508 = t562 + t604;
t467 = t494 * t713 + t508 * t545;
t493 = t522 * t546 - t547 * t591 + t557;
t410 = -t467 * t544 + t493 * t548;
t517 = t569 * qJD(1);
t566 = qJD(1) * t570;
t559 = t704 * t566;
t437 = qJD(3) * t494 + t517 * t546 - t549 * t586 - t559 * t714;
t737 = qJD(5) * t410 + t437 * t544;
t736 = t507 * pkin(9);
t598 = t704 * t703;
t600 = t705 * t702;
t735 = t543 * (-t701 * t546 + t714 * t598) + t714 * t600;
t585 = qJD(3) * t591;
t653 = qJD(1) * t546;
t590 = t624 * t653;
t438 = -qJD(3) * t557 + t517 * t714 - t522 * t649 + t546 * t559 + t547 * t585 + t549 * t590;
t500 = t507 * qJD(1);
t356 = qJD(4) * t467 + t438 * t545 + t500 * t713;
t731 = t356 / 0.2e1;
t440 = -qJD(3) * t558 + t518 * t714 + t547 * t590 - t549 * t585 - t563 * t653 + t569 * t649;
t501 = t508 * qJD(1);
t358 = qJD(4) * t465 + t440 * t545 - t501 * t713;
t730 = t358 / 0.2e1;
t729 = t437 / 0.2e1;
t728 = t439 / 0.2e1;
t506 = t546 * t600 + (t546 * t598 + t701 * t714) * t543;
t520 = -t624 * t703 + t704 * t705;
t490 = t506 * t713 + t520 * t545;
t498 = t735 * qJD(3);
t451 = qJD(4) * t490 + t498 * t545;
t727 = t451 / 0.2e1;
t589 = -t492 * t545 + t507 * t713;
t726 = -t589 / 0.2e1;
t588 = -t494 * t545 + t508 * t713;
t725 = -t588 / 0.2e1;
t587 = -t506 * t545 + t520 * t713;
t724 = -t587 / 0.2e1;
t723 = t491 / 0.2e1;
t722 = t493 / 0.2e1;
t499 = t506 * qJD(3);
t721 = t499 / 0.2e1;
t720 = -t500 / 0.2e1;
t719 = t501 / 0.2e1;
t718 = -t735 / 0.2e1;
t717 = t507 / 0.2e1;
t716 = t508 / 0.2e1;
t715 = t520 / 0.2e1;
t338 = t438 * rSges(4,1) - t437 * rSges(4,2) - t500 * rSges(4,3);
t339 = t440 * rSges(4,1) - t439 * rSges(4,2) + t501 * rSges(4,3);
t425 = rSges(4,1) * t492 - rSges(4,2) * t491 + rSges(4,3) * t507;
t426 = t494 * rSges(4,1) - t493 * rSges(4,2) + t508 * rSges(4,3);
t712 = m(4) * (-t338 * t507 + t339 * t508 - t425 * t500 - t426 * t501);
t540 = t549 * pkin(1);
t537 = pkin(5) * t548 + pkin(4);
t711 = -pkin(4) + t537;
t550 = -pkin(12) - pkin(11);
t710 = -pkin(11) - t550;
t542 = qJ(5) + qJ(6);
t538 = sin(t542);
t539 = cos(t542);
t541 = qJD(5) + qJD(6);
t611 = -t465 * t541 + t439;
t359 = qJD(4) * t589 + t440 * t713 + t501 * t545;
t613 = t491 * t541 + t359;
t272 = -t538 * t613 + t539 * t611;
t273 = t538 * t611 + t539 * t613;
t176 = Icges(7,5) * t273 + Icges(7,6) * t272 + Icges(7,3) * t358;
t178 = Icges(7,4) * t273 + Icges(7,2) * t272 + Icges(7,6) * t358;
t180 = Icges(7,1) * t273 + Icges(7,4) * t272 + Icges(7,5) * t358;
t402 = -t465 * t538 + t491 * t539;
t403 = t465 * t539 + t491 * t538;
t304 = Icges(7,5) * t403 + Icges(7,6) * t402 - Icges(7,3) * t589;
t306 = Icges(7,4) * t403 + Icges(7,2) * t402 - Icges(7,6) * t589;
t308 = Icges(7,1) * t403 + Icges(7,4) * t402 - Icges(7,5) * t589;
t609 = -t490 * t541 + t499;
t452 = qJD(4) * t587 + t498 * t713;
t610 = -t541 * t735 + t452;
t365 = -t538 * t610 + t539 * t609;
t366 = t538 * t609 + t539 * t610;
t448 = -t490 * t538 - t539 * t735;
t449 = t490 * t539 - t538 * t735;
t57 = -t176 * t587 + t178 * t448 + t180 * t449 + t304 * t451 + t306 * t365 + t308 * t366;
t708 = t57 * t589;
t612 = -t467 * t541 + t437;
t357 = qJD(4) * t588 + t438 * t713 - t500 * t545;
t614 = t493 * t541 + t357;
t270 = -t538 * t614 + t539 * t612;
t271 = t538 * t612 + t539 * t614;
t175 = Icges(7,5) * t271 + Icges(7,6) * t270 + Icges(7,3) * t356;
t177 = Icges(7,4) * t271 + Icges(7,2) * t270 + Icges(7,6) * t356;
t179 = Icges(7,1) * t271 + Icges(7,4) * t270 + Icges(7,5) * t356;
t404 = -t467 * t538 + t493 * t539;
t405 = t467 * t539 + t493 * t538;
t305 = Icges(7,5) * t405 + Icges(7,6) * t404 - Icges(7,3) * t588;
t307 = Icges(7,4) * t405 + Icges(7,2) * t404 - Icges(7,6) * t588;
t309 = Icges(7,1) * t405 + Icges(7,4) * t404 - Icges(7,5) * t588;
t58 = -t175 * t587 + t177 * t448 + t179 * t449 + t305 * t451 + t307 * t365 + t309 * t366;
t707 = t58 * t588;
t706 = t550 - rSges(7,3);
t166 = -t304 * t587 + t306 * t448 + t308 * t449;
t700 = t166 * t358;
t167 = -t305 * t587 + t307 * t448 + t309 * t449;
t699 = t167 * t356;
t413 = Icges(5,5) * t490 + Icges(5,6) * t587 - Icges(5,3) * t735;
t414 = Icges(5,4) * t490 + Icges(5,2) * t587 - Icges(5,6) * t735;
t415 = Icges(5,1) * t490 + Icges(5,4) * t587 - Icges(5,5) * t735;
t257 = -t413 * t735 + t414 * t587 + t415 * t490;
t251 = t257 * t499;
t419 = Icges(4,5) * t492 - Icges(4,6) * t491 + Icges(4,3) * t507;
t698 = t419 * t508;
t420 = Icges(4,5) * t494 - Icges(4,6) * t493 + Icges(4,3) * t508;
t697 = t420 * t507;
t693 = t491 * t544;
t692 = t493 * t544;
t369 = Icges(5,5) * t465 + Icges(5,6) * t589 + Icges(5,3) * t491;
t371 = Icges(5,4) * t465 + Icges(5,2) * t589 + Icges(5,6) * t491;
t373 = Icges(5,1) * t465 + Icges(5,4) * t589 + Icges(5,5) * t491;
t220 = -t369 * t735 + t371 * t587 + t373 * t490;
t690 = t499 * t220;
t689 = t735 * t544;
t688 = t543 * t547;
t687 = t543 * t549;
t252 = Icges(7,5) * t366 + Icges(7,6) * t365 + Icges(7,3) * t451;
t253 = Icges(7,4) * t366 + Icges(7,2) * t365 + Icges(7,6) * t451;
t254 = Icges(7,1) * t366 + Icges(7,4) * t365 + Icges(7,5) * t451;
t345 = Icges(7,5) * t449 + Icges(7,6) * t448 - Icges(7,3) * t587;
t346 = Icges(7,4) * t449 + Icges(7,2) * t448 - Icges(7,6) * t587;
t347 = Icges(7,1) * t449 + Icges(7,4) * t448 - Icges(7,5) * t587;
t116 = -t252 * t587 + t448 * t253 + t449 * t254 + t451 * t345 + t365 * t346 + t366 * t347;
t207 = -t345 * t587 + t346 * t448 + t347 * t449;
t686 = -t116 * t587 + t207 * t451;
t685 = -t116 * t735 + t207 * t499;
t456 = t490 * t548 - t689;
t380 = -qJD(5) * t456 - t452 * t544 + t499 * t548;
t455 = -t490 * t544 - t548 * t735;
t572 = qJD(5) * t455 + t499 * t544;
t381 = t452 * t548 + t572;
t262 = Icges(6,5) * t381 + Icges(6,6) * t380 + Icges(6,3) * t451;
t263 = Icges(6,4) * t381 + Icges(6,2) * t380 + Icges(6,6) * t451;
t264 = Icges(6,1) * t381 + Icges(6,4) * t380 + Icges(6,5) * t451;
t361 = Icges(6,5) * t456 + Icges(6,6) * t455 - Icges(6,3) * t587;
t362 = Icges(6,4) * t456 + Icges(6,2) * t455 - Icges(6,6) * t587;
t363 = Icges(6,1) * t456 + Icges(6,4) * t455 - Icges(6,5) * t587;
t120 = -t262 * t587 + t455 * t263 + t456 * t264 + t451 * t361 + t380 * t362 + t381 * t363;
t209 = -t361 * t587 + t362 * t455 + t363 * t456;
t684 = -t120 * t587 + t209 * t451;
t683 = -t120 * t735 + t209 * t499;
t382 = Icges(5,5) * t452 - Icges(5,6) * t451 + Icges(5,3) * t499;
t383 = Icges(5,4) * t452 - Icges(5,2) * t451 + Icges(5,6) * t499;
t384 = Icges(5,1) * t452 - Icges(5,4) * t451 + Icges(5,5) * t499;
t163 = -t382 * t735 + t383 * t587 + t490 * t384 + t499 * t413 - t451 * t414 + t452 * t415;
t682 = -t163 * t735 + t251;
t596 = -t273 * rSges(7,1) - t272 * rSges(7,2);
t182 = t358 * rSges(7,3) - t596;
t595 = -t403 * rSges(7,1) - t402 * rSges(7,2);
t310 = -rSges(7,3) * t589 - t595;
t681 = -t182 * t588 + t356 * t310;
t181 = t271 * rSges(7,1) + t270 * rSges(7,2) + t356 * rSges(7,3);
t311 = t405 * rSges(7,1) + t404 * rSges(7,2) - rSges(7,3) * t588;
t680 = -t181 * t587 + t451 * t311;
t288 = t357 * pkin(4) + t356 * pkin(11);
t581 = pkin(5) * t737 - t356 * t550 + t357 * t537;
t191 = -t288 + t581;
t679 = t181 + t191;
t355 = t358 * pkin(11);
t192 = -t358 * t550 + t359 * t711 - t355 - t738;
t678 = t182 + t192;
t411 = t467 * t548 + t692;
t279 = -qJD(5) * t411 - t357 * t544 + t437 * t548;
t280 = t357 * t548 + t737;
t189 = t280 * rSges(6,1) + t279 * rSges(6,2) + t356 * rSges(6,3);
t677 = -t189 - t288;
t409 = t465 * t548 + t693;
t281 = -qJD(5) * t409 - t359 * t544 + t439 * t548;
t282 = t359 * t548 - t732;
t190 = t282 * rSges(6,1) + t281 * rSges(6,2) + t358 * rSges(6,3);
t289 = t359 * pkin(4) + t355;
t676 = -t190 - t289;
t255 = rSges(7,1) * t366 + rSges(7,2) * t365 + rSges(7,3) * t451;
t348 = rSges(7,1) * t449 + rSges(7,2) * t448 - rSges(7,3) * t587;
t675 = -t255 * t589 + t358 * t348;
t276 = pkin(5) * t572 + t451 * t710 + t452 * t711;
t674 = t255 + t276;
t421 = Icges(4,4) * t492 - Icges(4,2) * t491 + Icges(4,6) * t507;
t422 = Icges(4,4) * t494 - Icges(4,2) * t493 + Icges(4,6) * t508;
t423 = Icges(4,1) * t492 - Icges(4,4) * t491 + Icges(4,5) * t507;
t424 = Icges(4,1) * t494 - Icges(4,4) * t493 + Icges(4,5) * t508;
t673 = -t421 * t493 - t422 * t491 + t423 * t494 + t424 * t492 + t697 + t698;
t265 = rSges(6,1) * t381 + rSges(6,2) * t380 + rSges(6,3) * t451;
t395 = t452 * pkin(4) + t451 * pkin(11);
t672 = -t265 - t395;
t462 = t589 * pkin(11);
t396 = t465 * pkin(4) - t462;
t671 = t493 * t289 + t437 * t396;
t397 = t467 * pkin(4) - pkin(11) * t588;
t670 = -t288 * t735 + t499 * t397;
t377 = t438 * pkin(3) + t437 * pkin(10);
t344 = t520 * t377;
t669 = t520 * t288 + t344;
t646 = pkin(5) * t693;
t300 = t465 * t711 + t550 * t589 + t462 + t646;
t668 = t300 + t310;
t635 = pkin(5) * t692 + t467 * t537 + t550 * t588;
t301 = -t397 + t635;
t667 = t301 + t311;
t318 = rSges(6,1) * t409 - rSges(6,2) * t593 - rSges(6,3) * t589;
t666 = -t318 - t396;
t319 = t411 * rSges(6,1) + t410 * rSges(6,2) - rSges(6,3) * t588;
t665 = -t319 - t397;
t378 = t440 * pkin(3) + t439 * pkin(10);
t442 = t492 * pkin(3) + t491 * pkin(10);
t664 = t508 * t378 - t500 * t442;
t341 = -pkin(5) * t689 + t490 * t711 - t587 * t710;
t663 = t341 + t348;
t441 = t490 * pkin(4) - pkin(11) * t587;
t662 = t491 * t395 + t439 * t441;
t364 = rSges(6,1) * t456 + rSges(6,2) * t455 - rSges(6,3) * t587;
t661 = -t364 - t441;
t376 = t467 * rSges(5,1) + rSges(5,2) * t588 + t493 * rSges(5,3);
t443 = t494 * pkin(3) + t493 * pkin(10);
t660 = -t376 - t443;
t427 = t508 * t442;
t659 = t508 * t396 + t427;
t428 = t520 * t443;
t658 = t520 * t397 + t428;
t416 = rSges(5,1) * t490 + rSges(5,2) * t587 - rSges(5,3) * t735;
t479 = pkin(3) * t506 - pkin(10) * t735;
t657 = -t416 - t479;
t453 = t507 * t479;
t656 = t507 * t441 + t453;
t476 = pkin(3) * t498 + pkin(10) * t499;
t655 = t507 * t476 + t501 * t479;
t654 = qJ(2) * t688 + t540;
t652 = qJD(1) * t547;
t651 = qJD(1) * t549;
t650 = qJD(2) * t543;
t184 = Icges(6,5) * t282 + Icges(6,6) * t281 + Icges(6,3) * t358;
t186 = Icges(6,4) * t282 + Icges(6,2) * t281 + Icges(6,6) * t358;
t188 = Icges(6,1) * t282 + Icges(6,4) * t281 + Icges(6,5) * t358;
t312 = Icges(6,5) * t409 - Icges(6,6) * t593 - Icges(6,3) * t589;
t314 = Icges(6,4) * t409 - Icges(6,2) * t593 - Icges(6,6) * t589;
t316 = Icges(6,1) * t409 - Icges(6,4) * t593 - Icges(6,5) * t589;
t59 = -t184 * t587 + t186 * t455 + t188 * t456 + t312 * t451 + t314 * t380 + t316 * t381;
t97 = -t262 * t589 - t263 * t593 + t264 * t409 + t281 * t362 + t282 * t363 + t358 * t361;
t645 = t59 / 0.2e1 + t97 / 0.2e1;
t183 = Icges(6,5) * t280 + Icges(6,6) * t279 + Icges(6,3) * t356;
t185 = Icges(6,4) * t280 + Icges(6,2) * t279 + Icges(6,6) * t356;
t187 = Icges(6,1) * t280 + Icges(6,4) * t279 + Icges(6,5) * t356;
t313 = Icges(6,5) * t411 + Icges(6,6) * t410 - Icges(6,3) * t588;
t315 = Icges(6,4) * t411 + Icges(6,2) * t410 - Icges(6,6) * t588;
t317 = Icges(6,1) * t411 + Icges(6,4) * t410 - Icges(6,5) * t588;
t60 = -t183 * t587 + t185 * t455 + t187 * t456 + t313 * t451 + t315 * t380 + t317 * t381;
t96 = -t262 * t588 + t263 * t410 + t264 * t411 + t279 * t362 + t280 * t363 + t356 * t361;
t644 = t60 / 0.2e1 + t96 / 0.2e1;
t643 = -t288 - t679;
t642 = -t289 - t678;
t641 = -t395 - t674;
t640 = -t396 - t668;
t639 = -t397 - t667;
t638 = -t443 + t665;
t637 = -t441 - t663;
t248 = t357 * rSges(5,1) - t356 * rSges(5,2) + t437 * rSges(5,3);
t636 = -t479 + t661;
t634 = m(5) * t705;
t633 = m(6) * t705;
t632 = m(7) * t705;
t629 = t543 * t651;
t168 = -t312 * t587 + t314 * t455 + t316 * t456;
t199 = -t361 * t589 - t362 * t593 + t363 * t409;
t628 = t168 / 0.2e1 + t199 / 0.2e1;
t169 = -t313 * t587 + t315 * t455 + t317 * t456;
t200 = -t361 * t588 + t362 * t410 + t363 * t411;
t627 = t169 / 0.2e1 + t200 / 0.2e1;
t626 = -t547 * pkin(1) + qJ(2) * t687;
t621 = 0.2e1 * m(4);
t619 = 0.2e1 * m(5);
t617 = 0.2e1 * m(6);
t615 = 0.2e1 * m(7);
t608 = t508 * t289 - t500 * t396 + t664;
t607 = -t443 + t639;
t606 = -t479 + t637;
t605 = t507 * t395 + t501 * t441 + t655;
t26 = t686 + t699 + t700 - t707 - t708;
t152 = -t304 * t588 + t306 * t404 + t308 * t405;
t153 = -t305 * t588 + t307 * t404 + t309 * t405;
t198 = -t345 * t588 + t346 * t404 + t347 * t405;
t43 = -t176 * t588 + t178 * t404 + t180 * t405 + t270 * t306 + t271 * t308 + t304 * t356;
t44 = -t175 * t588 + t177 * t404 + t179 * t405 + t270 * t307 + t271 * t309 + t305 * t356;
t90 = -t252 * t588 + t253 * t404 + t254 * t405 + t270 * t346 + t271 * t347 + t345 * t356;
t7 = t152 * t358 + t153 * t356 + t198 * t451 - t43 * t589 - t44 * t588 - t587 * t90;
t150 = -t304 * t589 + t306 * t402 + t308 * t403;
t151 = -t305 * t589 + t307 * t402 + t309 * t403;
t197 = -t345 * t589 + t346 * t402 + t347 * t403;
t75 = -t150 * t589 - t151 * t588 - t197 * t587;
t76 = -t152 * t589 - t153 * t588 - t198 * t587;
t45 = -t176 * t589 + t178 * t402 + t180 * t403 + t272 * t306 + t273 * t308 + t304 * t358;
t46 = -t175 * t589 + t177 * t402 + t179 * t403 + t272 * t307 + t273 * t309 + t305 * t358;
t91 = -t252 * t589 + t253 * t402 + t254 * t403 + t272 * t346 + t273 * t347 + t345 * t358;
t8 = t150 * t358 + t151 * t356 + t197 * t451 - t45 * t589 - t46 * t588 - t587 * t91;
t94 = -t166 * t589 - t167 * t588 - t207 * t587;
t597 = -t26 * t587 + t356 * t76 + t358 * t75 + t451 * t94 - t588 * t7 - t589 * t8;
t469 = Icges(4,4) * t506 + Icges(4,2) * t735 + Icges(4,6) * t520;
t470 = Icges(4,1) * t506 + Icges(4,4) * t735 + Icges(4,5) * t520;
t472 = Icges(4,5) * t498 - Icges(4,6) * t499;
t473 = Icges(4,4) * t498 - Icges(4,2) * t499;
t474 = Icges(4,1) * t498 - Icges(4,4) * t499;
t594 = -t499 * t469 + t498 * t470 + t520 * t472 + t473 * t735 + t506 * t474;
t592 = -pkin(1) * t652 + qJ(2) * t629 + t547 * t650;
t249 = t359 * rSges(5,1) - t358 * rSges(5,2) + t439 * rSges(5,3);
t375 = rSges(5,1) * t465 + rSges(5,2) * t589 + rSges(5,3) * t491;
t582 = pkin(2) * t569 + t626 - t736;
t243 = Icges(5,5) * t359 - Icges(5,6) * t358 + Icges(5,3) * t439;
t245 = Icges(5,4) * t359 - Icges(5,2) * t358 + Icges(5,6) * t439;
t247 = Icges(5,1) * t359 - Icges(5,4) * t358 + Icges(5,5) * t439;
t125 = -t243 * t735 + t245 * t587 + t247 * t490 + t369 * t499 - t371 * t451 + t373 * t452;
t143 = -t358 * t414 + t359 * t415 + t382 * t491 + t383 * t589 + t384 * t465 + t413 * t439;
t580 = t125 / 0.2e1 + t57 / 0.2e1 + t143 / 0.2e1 + t91 / 0.2e1 + t645;
t242 = Icges(5,5) * t357 - Icges(5,6) * t356 + Icges(5,3) * t437;
t244 = Icges(5,4) * t357 - Icges(5,2) * t356 + Icges(5,6) * t437;
t246 = Icges(5,1) * t357 - Icges(5,4) * t356 + Icges(5,5) * t437;
t370 = Icges(5,5) * t467 + Icges(5,6) * t588 + Icges(5,3) * t493;
t372 = Icges(5,4) * t467 + Icges(5,2) * t588 + Icges(5,6) * t493;
t374 = Icges(5,1) * t467 + Icges(5,4) * t588 + Icges(5,5) * t493;
t126 = -t242 * t735 + t244 * t587 + t246 * t490 + t370 * t499 - t372 * t451 + t374 * t452;
t142 = -t356 * t414 + t357 * t415 + t382 * t493 + t383 * t588 + t384 * t467 + t413 * t437;
t579 = t126 / 0.2e1 + t58 / 0.2e1 + t142 / 0.2e1 + t90 / 0.2e1 + t644;
t102 = t166 * t507 + t167 * t508 + t207 * t520;
t19 = t152 * t501 - t153 * t500 + t43 * t507 + t44 * t508 + t520 * t90;
t20 = t150 * t501 - t151 * t500 + t45 * t507 + t46 * t508 + t520 * t91;
t115 = t116 * t520;
t32 = t166 * t501 - t167 * t500 + t57 * t507 + t58 * t508 + t115;
t81 = t150 * t507 + t151 * t508 + t197 * t520;
t82 = t152 * t507 + t153 * t508 + t198 * t520;
t578 = t102 * t727 + t19 * t725 + t20 * t726 + t26 * t715 + t32 * t724 + t7 * t716 + t8 * t717 + t75 * t719 + t76 * t720 + t81 * t730 + t731 * t82;
t234 = t413 * t491 + t414 * t589 + t415 * t465;
t577 = t220 / 0.2e1 + t166 / 0.2e1 + t234 / 0.2e1 + t197 / 0.2e1 + t628;
t221 = -t370 * t735 + t372 * t587 + t374 * t490;
t235 = t413 * t493 + t414 * t588 + t415 * t467;
t576 = t221 / 0.2e1 + t167 / 0.2e1 + t235 / 0.2e1 + t198 / 0.2e1 + t627;
t575 = t700 / 0.2e1 + t699 / 0.2e1 + t198 * t731 + t197 * t730 - t708 / 0.2e1 - t707 / 0.2e1 + t91 * t726 + t90 * t725 + t686;
t101 = t166 * t491 + t167 * t493 - t207 * t735;
t11 = t152 * t439 + t153 * t437 + t198 * t499 + t43 * t491 + t44 * t493 - t735 * t90;
t12 = t150 * t439 + t151 * t437 + t197 * t499 + t45 * t491 + t46 * t493 - t735 * t91;
t28 = t166 * t439 + t167 * t437 + t57 * t491 + t58 * t493 + t685;
t79 = t150 * t491 + t151 * t493 - t197 * t735;
t80 = t152 * t491 + t153 * t493 - t198 * t735;
t574 = t101 * t727 + t11 * t725 + t12 * t726 + t26 * t718 + t28 * t724 + t7 * t722 + t94 * t721 + t8 * t723 + t75 * t728 + t76 * t729 + t79 * t730 + t731 * t80;
t568 = -t442 + t582;
t567 = t571 * rSges(3,2);
t560 = pkin(9) * t562;
t534 = t549 * t650;
t556 = -t518 * pkin(2) + t534 + (-t560 - t540 + (-pkin(9) * t704 - qJ(2)) * t688) * qJD(1);
t555 = t522 * pkin(2) + pkin(9) * t604 + t560 + t654;
t554 = -t378 + t556;
t553 = t443 + t555;
t552 = t517 * pkin(2) - qJD(1) * t736 + t592;
t551 = t377 + t552;
t488 = t522 * rSges(3,1) + rSges(3,3) * t688 - t567 + t654;
t487 = rSges(3,1) * t569 + rSges(3,2) * t570 + rSges(3,3) * t687 + t626;
t478 = -t518 * rSges(3,1) - pkin(1) * t651 + qJD(1) * t567 + t534 + (-rSges(3,3) - qJ(2)) * t543 * t652;
t477 = t517 * rSges(3,1) + rSges(3,2) * t566 + rSges(3,3) * t629 + t592;
t475 = rSges(4,1) * t498 - rSges(4,2) * t499;
t471 = rSges(4,1) * t506 + rSges(4,2) * t735 + rSges(4,3) * t520;
t468 = Icges(4,5) * t506 + Icges(4,6) * t735 + Icges(4,3) * t520;
t399 = t491 * t441;
t393 = t555 + t426;
t392 = -t425 + t582;
t390 = t735 * t397;
t385 = rSges(5,1) * t452 - rSges(5,2) * t451 + rSges(5,3) * t499;
t379 = t493 * t396;
t337 = t426 * t520 - t471 * t508;
t336 = -t425 * t520 + t471 * t507;
t335 = Icges(4,1) * t440 - Icges(4,4) * t439 + Icges(4,5) * t501;
t334 = Icges(4,1) * t438 - Icges(4,4) * t437 - Icges(4,5) * t500;
t333 = Icges(4,4) * t440 - Icges(4,2) * t439 + Icges(4,6) * t501;
t332 = Icges(4,4) * t438 - Icges(4,2) * t437 - Icges(4,6) * t500;
t331 = Icges(4,5) * t440 - Icges(4,6) * t439 + Icges(4,3) * t501;
t330 = Icges(4,5) * t438 - Icges(4,6) * t437 - Icges(4,3) * t500;
t325 = t589 * t348;
t323 = -t339 + t556;
t322 = t552 + t338;
t321 = t468 * t508 - t469 * t493 + t470 * t494;
t320 = t468 * t507 - t469 * t491 + t470 * t492;
t303 = t553 + t376;
t302 = -t375 + t568;
t299 = t587 * t311;
t298 = -t376 * t735 - t416 * t493;
t297 = t375 * t735 + t416 * t491;
t294 = t588 * t310;
t287 = t420 * t520 + t422 * t735 + t424 * t506;
t286 = t419 * t520 + t421 * t735 + t423 * t506;
t285 = -t339 * t520 + t471 * t501 + t475 * t507;
t284 = t338 * t520 + t471 * t500 - t475 * t508;
t266 = t594 * t520;
t256 = t375 * t493 - t376 * t491;
t239 = t376 * t520 + t508 * t657 + t428;
t238 = t416 * t507 + t453 + (-t375 - t442) * t520;
t230 = t553 - t665;
t229 = t568 + t666;
t228 = -t319 * t587 + t364 * t588;
t227 = t318 * t587 - t364 * t589;
t226 = t348 * t588 - t299;
t225 = t310 * t587 - t325;
t224 = t375 * t508 + t507 * t660 + t427;
t223 = t553 + t635 + t311;
t222 = -t465 * t537 - t589 * t706 + t568 + t595 - t646;
t219 = -t439 * t469 + t440 * t470 + t468 * t501 + t472 * t507 - t473 * t491 + t474 * t492;
t218 = -t437 * t469 + t438 * t470 - t468 * t500 + t472 * t508 - t473 * t493 + t474 * t494;
t217 = -t249 + t554;
t216 = t551 + t248;
t215 = t370 * t493 + t372 * t588 + t374 * t467;
t214 = t369 * t493 + t371 * t588 + t373 * t467;
t213 = t370 * t491 + t372 * t589 + t374 * t465;
t212 = t369 * t491 + t371 * t589 + t373 * t465;
t210 = -t318 * t588 + t319 * t589;
t208 = t311 * t589 - t294;
t202 = -t319 * t735 + t493 * t661 - t390;
t201 = t364 * t491 - t666 * t735 + t399;
t196 = t319 * t520 + t508 * t636 + t658;
t195 = t364 * t507 + (-t442 + t666) * t520 + t656;
t194 = t330 * t520 + t332 * t735 + t334 * t506 - t422 * t499 + t424 * t498;
t193 = t331 * t520 + t333 * t735 + t335 * t506 - t421 * t499 + t423 * t498;
t174 = t318 * t493 + t491 * t665 + t379;
t165 = t249 * t735 - t375 * t499 + t385 * t491 + t416 * t439;
t164 = -t248 * t735 + t376 * t499 - t385 * t493 - t416 * t437;
t162 = t163 * t520;
t160 = t385 * t507 + t416 * t501 + (-t249 - t378) * t520 + t655;
t159 = t248 * t520 + t344 + (-t385 - t476) * t508 - t657 * t500;
t158 = t318 * t508 + t507 * t638 + t659;
t157 = -t313 * t588 + t315 * t410 + t317 * t411;
t156 = -t312 * t588 + t314 * t410 + t316 * t411;
t155 = -t313 * t589 - t315 * t593 + t317 * t409;
t154 = -t312 * t589 - t314 * t593 + t316 * t409;
t149 = -t301 * t587 + t588 * t663 - t299;
t148 = -t341 * t589 + t587 * t668 - t325;
t145 = t493 * t637 - t667 * t735 - t390;
t144 = t491 * t663 - t640 * t735 + t399;
t141 = t508 * t606 + t520 * t667 + t658;
t140 = t663 * t507 + (-t442 + t640) * t520 + t656;
t139 = -t300 * t588 + t589 * t667 - t294;
t138 = -t248 * t491 + t249 * t493 + t375 * t437 - t376 * t439;
t137 = t554 + t676;
t136 = t551 - t677;
t135 = t491 * t639 + t493 * t668 + t379;
t134 = t706 * t358 - t359 * t537 + t554 + t596 + t738;
t133 = t551 + t581 + t181;
t132 = t507 * t607 + t508 * t668 + t659;
t131 = t214 * t507 + t215 * t508 + t235 * t520;
t130 = t212 * t507 + t213 * t508 + t234 * t520;
t129 = t249 * t508 - t375 * t500 + (-t248 - t377) * t507 + t660 * t501 + t664;
t128 = t214 * t491 + t215 * t493 - t235 * t735;
t127 = t212 * t491 + t213 * t493 - t234 * t735;
t124 = t190 * t587 - t265 * t589 - t318 * t451 + t358 * t364;
t123 = -t189 * t587 + t265 * t588 + t319 * t451 - t356 * t364;
t122 = t182 * t587 - t310 * t451 + t675;
t121 = t255 * t588 - t348 * t356 + t680;
t119 = t120 * t520;
t112 = t242 * t491 + t244 * t589 + t246 * t465 - t358 * t372 + t359 * t374 + t370 * t439;
t111 = t243 * t491 + t245 * t589 + t247 * t465 - t358 * t371 + t359 * t373 + t369 * t439;
t110 = t242 * t493 + t244 * t588 + t246 * t467 - t356 * t372 + t357 * t374 + t370 * t437;
t109 = t243 * t493 + t245 * t588 + t247 * t467 - t356 * t371 + t357 * t373 + t369 * t437;
t108 = t265 * t507 + t364 * t501 + (-t378 + t676) * t520 + t605;
t107 = t189 * t520 + (-t476 + t672) * t508 - t636 * t500 + t669;
t106 = t168 * t507 + t169 * t508 + t209 * t520;
t105 = t168 * t491 + t169 * t493 - t209 * t735;
t104 = t265 * t491 + t364 * t439 + t499 * t666 - t676 * t735 + t662;
t103 = -t189 * t735 + t319 * t499 + t437 * t661 + t493 * t672 + t670;
t100 = -t168 * t589 - t169 * t588 - t209 * t587;
t99 = t189 * t589 - t190 * t588 + t318 * t356 - t319 * t358;
t95 = t181 * t589 - t311 * t358 + t681;
t86 = t156 * t507 + t157 * t508 + t200 * t520;
t85 = t154 * t507 + t155 * t508 + t199 * t520;
t84 = t156 * t491 + t157 * t493 - t200 * t735;
t83 = t154 * t491 + t155 * t493 - t199 * t735;
t78 = -t156 * t589 - t157 * t588 - t200 * t587;
t77 = -t154 * t589 - t155 * t588 - t199 * t587;
t64 = t190 * t493 + t318 * t437 + t439 * t665 + t491 * t677 + t671;
t63 = t674 * t507 + t663 * t501 + (-t378 + t642) * t520 + t605;
t62 = t679 * t520 + (-t476 + t641) * t508 - t606 * t500 + t669;
t61 = t190 * t508 - t318 * t500 + (-t377 + t677) * t507 + t638 * t501 + t608;
t56 = -t276 * t589 + t341 * t358 - t451 * t668 + t587 * t678 + t675;
t55 = -t191 * t587 + t301 * t451 - t356 * t663 + t588 * t674 + t680;
t52 = t439 * t663 + t491 * t674 + t499 * t640 - t642 * t735 + t662;
t51 = t437 * t637 + t493 * t641 + t499 * t667 - t679 * t735 + t670;
t50 = -t183 * t589 - t185 * t593 + t187 * t409 + t281 * t315 + t282 * t317 + t313 * t358;
t49 = -t184 * t589 - t186 * t593 + t188 * t409 + t281 * t314 + t282 * t316 + t312 * t358;
t48 = -t183 * t588 + t185 * t410 + t187 * t411 + t279 * t315 + t280 * t317 + t313 * t356;
t47 = -t184 * t588 + t186 * t410 + t188 * t411 + t279 * t314 + t280 * t316 + t312 * t356;
t42 = -t192 * t588 + t300 * t356 - t358 * t667 + t589 * t679 + t681;
t41 = t125 * t507 + t126 * t508 + t220 * t501 - t221 * t500 + t162;
t40 = t437 * t668 + t439 * t639 + t491 * t643 + t493 * t678 + t671;
t39 = t678 * t508 - t668 * t500 + (-t377 + t643) * t507 + t607 * t501 + t608;
t38 = t125 * t491 + t126 * t493 + t220 * t439 + t221 * t437 + t682;
t37 = t111 * t507 + t112 * t508 + t143 * t520 + t212 * t501 - t213 * t500;
t36 = t109 * t507 + t110 * t508 + t142 * t520 + t214 * t501 - t215 * t500;
t35 = t111 * t491 + t112 * t493 - t143 * t735 + t212 * t439 + t213 * t437 + t234 * t499;
t34 = t109 * t491 + t110 * t493 - t142 * t735 + t214 * t439 + t215 * t437 + t235 * t499;
t33 = t168 * t501 - t169 * t500 + t59 * t507 + t60 * t508 + t119;
t30 = t168 * t439 + t169 * t437 + t59 * t491 + t60 * t493 + t683;
t29 = t168 * t358 + t169 * t356 - t588 * t60 - t589 * t59 + t684;
t22 = t154 * t501 - t155 * t500 + t49 * t507 + t50 * t508 + t520 * t97;
t21 = t156 * t501 - t157 * t500 + t47 * t507 + t48 * t508 + t520 * t96;
t16 = t154 * t439 + t155 * t437 + t199 * t499 + t49 * t491 + t493 * t50 - t735 * t97;
t15 = t156 * t439 + t157 * t437 + t200 * t499 + t47 * t491 + t48 * t493 - t735 * t96;
t14 = t154 * t358 + t155 * t356 + t199 * t451 - t49 * t589 - t50 * t588 - t587 * t97;
t13 = t156 * t358 + t157 * t356 + t200 * t451 - t47 * t589 - t48 * t588 - t587 * t96;
t1 = [t120 + t116 + (t133 * t223 + t134 * t222) * t615 + (t136 * t230 + t137 * t229) * t617 + (t216 * t303 + t217 * t302) * t619 + (t322 * t393 + t323 * t392) * t621 + 0.2e1 * m(3) * (t477 * t488 + t478 * t487) + t163 + t594; ((-t133 * t549 + t134 * t547 + t222 * t651 + t223 * t652) * m(7) + (-t136 * t549 + t137 * t547 + t229 * t651 + t230 * t652) * m(6) + (-t216 * t549 + t217 * t547 + t302 * t651 + t303 * t652) * m(5) + (-t322 * t549 + t323 * t547 + t392 * t651 + t393 * t652) * m(4) + m(3) * (-t477 * t549 + t478 * t547 + t487 * t651 + t488 * t652)) * t543; 0; t115 + (t284 * t393 + t285 * t392 + t322 * t337 + t323 * t336) * m(4) + (t159 * t303 + t160 * t302 + t216 * t239 + t217 * t238) * m(5) + (t107 * t230 + t108 * t229 + t136 * t196 + t137 * t195) * m(6) + (t133 * t141 + t134 * t140 + t222 * t63 + t223 * t62) * m(7) + t162 + t119 + t266 + (t194 / 0.2e1 + t218 / 0.2e1 + t579) * t508 + (t193 / 0.2e1 + t219 / 0.2e1 + t580) * t507 + (t286 / 0.2e1 + t320 / 0.2e1 + t577) * t501 - (t287 / 0.2e1 + t321 / 0.2e1 + t576) * t500; t705 * t712 + t129 * t634 + t61 * t633 + t39 * t632 + ((-t284 * t549 + t285 * t547 + t336 * t651 + t337 * t652) * m(4) + (-t159 * t549 + t160 * t547 + t238 * t651 + t239 * t652) * m(5) + (-t107 * t549 + t108 * t547 + t195 * t651 + t196 * t652) * m(6) + (t140 * t651 + t141 * t652 + t547 * t63 - t549 * t62) * m(7)) * t543; (t132 * t39 + t140 * t63 + t141 * t62) * t615 + (t107 * t196 + t108 * t195 + t158 * t61) * t617 + (t129 * t224 + t159 * t239 + t160 * t238) * t619 + (t337 * t284 + t336 * t285) * t621 + (0.2e1 * t425 * t712 + t19 + t21 + t36 + (t330 * t508 - t332 * t493 + t334 * t494 - t422 * t437 + t424 * t438) * t508) * t508 + (t20 + t22 + t37 - 0.2e1 * t426 * t712 + (t331 * t507 - t333 * t491 + t335 * t492 - t421 * t439 + t423 * t440) * t507 + (t330 * t507 + t331 * t508 - t332 * t491 - t333 * t493 + t334 * t492 + t335 * t494 - t421 * t437 - t422 * t439 + t423 * t438 + t424 * t440) * t508) * t507 + (t266 + t32 + t33 + t41 + (t194 + t218) * t508 + (t193 + t219) * t507) * t520 + (t81 + t85 + t130 + (0.3e1 * t419 * t507 - 0.2e1 * t421 * t491 + 0.2e1 * t423 * t492) * t507 + (t286 + t320) * t520 + (t673 + t697) * t508) * t501 - (t82 + t86 + t131 + (0.3e1 * t420 * t508 - 0.2e1 * t422 * t493 + 0.2e1 * t424 * t494) * t508 + (t321 + t287) * t520 + (t673 + t698) * t507) * t500; (t164 * t303 + t165 * t302 + t216 * t298 + t217 * t297) * m(5) + (t103 * t230 + t104 * t229 + t136 * t202 + t137 * t201) * m(6) + (t133 * t145 + t134 * t144 + t222 * t52 + t223 * t51) * m(7) + t579 * t493 + t580 * t491 + t577 * t439 + t576 * t437 + t682 + t683 + t685; t138 * t634 + t64 * t633 + t40 * t632 + ((-t164 * t549 + t165 * t547 + t297 * t651 + t298 * t652) * m(5) + (-t103 * t549 + t104 * t547 + t201 * t651 + t202 * t652) * m(6) + (t144 * t651 + t145 * t652 - t51 * t549 + t52 * t547) * m(7)) * t543; (t129 * t256 + t138 * t224 + t159 * t298 + t160 * t297 + t164 * t239 + t165 * t238) * m(5) + (t103 * t196 + t104 * t195 + t107 * t202 + t108 * t201 + t158 * t64 + t174 * t61) * m(6) + (t132 * t40 + t135 * t39 + t140 * t52 + t141 * t51 + t144 * t63 + t145 * t62) * m(7) + (t30 / 0.2e1 + t28 / 0.2e1 + t38 / 0.2e1) * t520 + (t15 / 0.2e1 + t11 / 0.2e1 + t34 / 0.2e1) * t508 + (t221 * t716 + t257 * t715 + t106 / 0.2e1 + t102 / 0.2e1) * t499 + (t690 / 0.2e1 + t16 / 0.2e1 + t12 / 0.2e1 + t35 / 0.2e1) * t507 - (t41 / 0.2e1 + t32 / 0.2e1 + t33 / 0.2e1) * t735 + (t127 / 0.2e1 + t83 / 0.2e1 + t79 / 0.2e1) * t501 - (t128 / 0.2e1 + t84 / 0.2e1 + t80 / 0.2e1) * t500 + (t21 / 0.2e1 + t19 / 0.2e1 + t36 / 0.2e1) * t493 + (t22 / 0.2e1 + t20 / 0.2e1 + t37 / 0.2e1) * t491 + (t130 / 0.2e1 + t85 / 0.2e1 + t81 / 0.2e1) * t439 + (t131 / 0.2e1 + t86 / 0.2e1 + t82 / 0.2e1) * t437; (t135 * t40 + t144 * t52 + t145 * t51) * t615 + (t103 * t202 + t104 * t201 + t174 * t64) * t617 + (t138 * t256 + t164 * t298 + t165 * t297) * t619 + (t101 + t105) * t499 - (t28 + t30 + t38 + t251) * t735 + (t499 * t221 + t11 + t15 + t34) * t493 + (t12 + t16 + t35 + t690) * t491 + (t79 + t83 + t127) * t439 + (t80 + t84 + t128) * t437; -t644 * t588 - t645 * t589 + t627 * t356 + t628 * t358 + t575 + (t123 * t230 + t124 * t229 + t136 * t228 + t137 * t227) * m(6) + (t133 * t149 + t134 * t148 + t222 * t56 + t223 * t55) * m(7) + t684; t99 * t633 + t42 * t632 + ((-t123 * t549 + t124 * t547 + t227 * t651 + t228 * t652) * m(6) + (t148 * t651 + t149 * t652 + t547 * t56 - t549 * t55) * m(7)) * t543; t29 * t715 + t13 * t716 + t14 * t717 + t77 * t719 + t78 * t720 + t33 * t724 + t21 * t725 + t22 * t726 + t106 * t727 + t85 * t730 + t86 * t731 + t578 + (t107 * t228 + t108 * t227 + t123 * t196 + t124 * t195 + t158 * t99 + t210 * t61) * m(6) + (t132 * t42 + t139 * t39 + t140 * t56 + t141 * t55 + t148 * t63 + t149 * t62) * m(7); t29 * t718 + t100 * t721 + t13 * t722 + t14 * t723 + t30 * t724 + t15 * t725 + t16 * t726 + t105 * t727 + t77 * t728 + t78 * t729 + t83 * t730 + t84 * t731 + t574 + (t103 * t228 + t104 * t227 + t123 * t202 + t124 * t201 + t174 * t99 + t210 * t64) * m(6) + (t135 * t42 + t139 * t40 + t144 * t56 + t145 * t55 + t148 * t52 + t149 * t51) * m(7); (t139 * t42 + t148 * t56 + t149 * t55) * t615 + (t123 * t228 + t124 * t227 + t210 * t99) * t617 + t358 * t77 - t589 * t14 + t451 * t100 - t587 * t29 + t356 * t78 - t588 * t13 + t597; (t121 * t223 + t122 * t222 + t133 * t226 + t134 * t225) * m(7) + t575; m(7) * (t95 * t705 + (-t121 * t549 + t122 * t547 + (t225 * t549 + t226 * t547) * qJD(1)) * t543); (t121 * t141 + t122 * t140 + t132 * t95 + t208 * t39 + t225 * t63 + t226 * t62) * m(7) + t578; t574 + (t121 * t145 + t122 * t144 + t135 * t95 + t208 * t40 + t225 * t52 + t226 * t51) * m(7); (t121 * t149 + t122 * t148 + t139 * t95 + t208 * t42 + t225 * t56 + t226 * t55) * m(7) + t597; (t121 * t226 + t122 * t225 + t208 * t95) * t615 + t597;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
