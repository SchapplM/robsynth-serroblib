% Calculate time derivative of joint inertia matrix for
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR6_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:49
% EndTime: 2019-03-09 15:49:03
% DurationCPUTime: 45.55s
% Computational Cost: add. (92715->1389), mult. (181676->1799), div. (0->0), fcn. (199078->12), ass. (0->600)
t551 = cos(pkin(6));
t555 = sin(qJ(2));
t550 = sin(pkin(6));
t706 = qJ(3) + pkin(11);
t640 = sin(t706);
t601 = t550 * t640;
t594 = t555 * t601;
t641 = cos(t706);
t502 = -t551 * t641 + t594;
t602 = t550 * t641;
t503 = t551 * t640 + t555 * t602;
t756 = cos(qJ(2));
t673 = t550 * t756;
t417 = -Icges(6,5) * t673 - Icges(6,6) * t503 + Icges(6,3) * t502;
t421 = Icges(5,4) * t503 - Icges(5,2) * t502 - Icges(5,6) * t673;
t822 = t417 - t421;
t418 = -Icges(6,4) * t673 - Icges(6,2) * t503 + Icges(6,6) * t502;
t422 = Icges(5,1) * t503 - Icges(5,4) * t502 - Icges(5,5) * t673;
t821 = -t418 + t422;
t419 = -Icges(6,1) * t673 - Icges(6,4) * t503 + Icges(6,5) * t502;
t420 = Icges(5,5) * t503 - Icges(5,6) * t502 - Icges(5,3) * t673;
t554 = sin(qJ(3));
t558 = cos(qJ(3));
t742 = t551 * t558;
t746 = t550 * t555;
t520 = -t554 * t746 + t742;
t743 = t551 * t554;
t744 = t550 * t558;
t521 = t555 * t744 + t743;
t432 = Icges(4,5) * t521 + Icges(4,6) * t520 - Icges(4,3) * t673;
t820 = t419 + t420 + t432;
t757 = cos(qJ(1));
t663 = t757 * t555;
t556 = sin(qJ(1));
t672 = t556 * t756;
t575 = -t551 * t663 - t672;
t576 = -t551 * t672 - t663;
t460 = qJD(1) * t575 + qJD(2) * t576;
t629 = t757 * t756;
t741 = t555 * t556;
t525 = -t551 * t741 + t629;
t592 = t556 * t601;
t480 = t525 * t641 + t592;
t586 = t757 * t602;
t318 = -qJD(1) * t586 + qJD(3) * t480 + t460 * t640;
t588 = t757 * t640;
t598 = qJD(3) * t640;
t599 = qJD(3) * t641;
t319 = t460 * t641 - t525 * t598 + (qJD(1) * t588 + t556 * t599) * t550;
t597 = t551 * t629;
t459 = -qJD(1) * t597 - qJD(2) * t629 + (qJD(2) * t551 + qJD(1)) * t741;
t209 = -Icges(6,5) * t459 - Icges(6,6) * t319 + Icges(6,3) * t318;
t215 = Icges(5,4) * t319 - Icges(5,2) * t318 - Icges(5,6) * t459;
t819 = -t209 + t215;
t574 = t597 - t741;
t462 = qJD(1) * t525 + qJD(2) * t574;
t585 = t550 * t588;
t593 = t556 * t602;
t320 = -qJD(1) * t593 - qJD(3) * t585 + t462 * t640 - t575 * t599;
t321 = qJD(1) * t592 - qJD(3) * t586 + t462 * t641 + t575 * t598;
t461 = -qJD(1) * t576 - qJD(2) * t575;
t210 = Icges(6,5) * t461 - Icges(6,6) * t321 + Icges(6,3) * t320;
t216 = Icges(5,4) * t321 - Icges(5,2) * t320 + Icges(5,6) * t461;
t818 = -t210 + t216;
t213 = -Icges(6,4) * t459 - Icges(6,2) * t319 + Icges(6,6) * t318;
t219 = Icges(5,1) * t319 - Icges(5,4) * t318 - Icges(5,5) * t459;
t817 = t213 - t219;
t214 = Icges(6,4) * t461 - Icges(6,2) * t321 + Icges(6,6) * t320;
t220 = Icges(5,1) * t321 - Icges(5,4) * t320 + Icges(5,5) * t461;
t816 = t214 - t220;
t477 = -t575 * t640 + t586;
t478 = -t575 * t641 - t585;
t327 = -Icges(6,5) * t574 - Icges(6,6) * t478 + Icges(6,3) * t477;
t333 = Icges(5,4) * t478 - Icges(5,2) * t477 - Icges(5,6) * t574;
t815 = t327 - t333;
t479 = t525 * t640 - t593;
t328 = -Icges(6,5) * t576 - Icges(6,6) * t480 + Icges(6,3) * t479;
t334 = Icges(5,4) * t480 - Icges(5,2) * t479 - Icges(5,6) * t576;
t814 = -t328 + t334;
t331 = -Icges(6,4) * t574 - Icges(6,2) * t478 + Icges(6,6) * t477;
t337 = Icges(5,1) * t478 - Icges(5,4) * t477 - Icges(5,5) * t574;
t813 = t331 - t337;
t332 = -Icges(6,4) * t576 - Icges(6,2) * t480 + Icges(6,6) * t479;
t338 = Icges(5,1) * t480 - Icges(5,4) * t479 - Icges(5,5) * t576;
t812 = t332 - t338;
t211 = Icges(5,5) * t319 - Icges(5,6) * t318 - Icges(5,3) * t459;
t217 = -Icges(6,1) * t459 - Icges(6,4) * t319 + Icges(6,5) * t318;
t745 = t550 * t556;
t687 = t554 * t745;
t489 = t525 * t558 + t687;
t654 = qJD(1) * t757;
t628 = t550 * t654;
t372 = -qJD(3) * t489 - t460 * t554 + t558 * t628;
t688 = t556 * t744;
t488 = -t525 * t554 + t688;
t596 = t554 * t628;
t373 = qJD(3) * t488 + t460 * t558 + t596;
t241 = Icges(4,5) * t373 + Icges(4,6) * t372 - Icges(4,3) * t459;
t811 = t211 + t217 + t241;
t212 = Icges(5,5) * t321 - Icges(5,6) * t320 + Icges(5,3) * t461;
t218 = Icges(6,1) * t461 - Icges(6,4) * t321 + Icges(6,5) * t320;
t674 = t550 * t757;
t632 = t554 * t674;
t583 = t558 * t575 + t632;
t711 = qJD(1) * t556;
t659 = t550 * t711;
t374 = qJD(3) * t583 - t462 * t554 + t558 * t659;
t631 = t558 * t674;
t486 = t554 * t575 - t631;
t630 = t554 * t659;
t375 = qJD(3) * t486 + t462 * t558 + t630;
t242 = Icges(4,5) * t375 + Icges(4,6) * t374 + Icges(4,3) * t461;
t810 = t212 + t218 + t242;
t329 = Icges(5,5) * t478 - Icges(5,6) * t477 - Icges(5,3) * t574;
t335 = -Icges(6,1) * t574 - Icges(6,4) * t478 + Icges(6,5) * t477;
t378 = -Icges(4,5) * t583 + Icges(4,6) * t486 - Icges(4,3) * t574;
t809 = t329 + t335 + t378;
t330 = Icges(5,5) * t480 - Icges(5,6) * t479 - Icges(5,3) * t576;
t336 = -Icges(6,1) * t576 - Icges(6,4) * t480 + Icges(6,5) * t479;
t379 = Icges(4,5) * t489 + Icges(4,6) * t488 - Icges(4,3) * t576;
t808 = t330 + t336 + t379;
t653 = qJD(2) * t756;
t463 = qJD(3) * t503 + t601 * t653;
t464 = -qJD(3) * t594 + t551 * t599 + t602 * t653;
t709 = qJD(2) * t555;
t658 = t550 * t709;
t357 = Icges(6,5) * t658 - Icges(6,6) * t464 + Icges(6,3) * t463;
t358 = Icges(6,4) * t658 - Icges(6,2) * t464 + Icges(6,6) * t463;
t359 = Icges(6,1) * t658 - Icges(6,4) * t464 + Icges(6,5) * t463;
t360 = Icges(5,5) * t464 - Icges(5,6) * t463 + Icges(5,3) * t658;
t361 = Icges(5,4) * t464 - Icges(5,2) * t463 + Icges(5,6) * t658;
t362 = Icges(5,1) * t464 - Icges(5,4) * t463 + Icges(5,5) * t658;
t627 = t550 * t653;
t484 = -qJD(3) * t521 - t554 * t627;
t485 = qJD(3) * t520 + t558 * t627;
t398 = Icges(4,5) * t485 + Icges(4,6) * t484 + Icges(4,3) * t658;
t399 = Icges(4,4) * t485 + Icges(4,2) * t484 + Icges(4,6) * t658;
t400 = Icges(4,1) * t485 + Icges(4,4) * t484 + Icges(4,5) * t658;
t433 = Icges(4,4) * t521 + Icges(4,2) * t520 - Icges(4,6) * t673;
t434 = Icges(4,1) * t521 + Icges(4,4) * t520 - Icges(4,5) * t673;
t807 = t520 * t399 + t521 * t400 + t484 * t433 + t485 * t434 + (-t358 + t362) * t503 + (t357 - t361) * t502 + t821 * t464 + t822 * t463 + (-t359 - t360 - t398) * t673 + t820 * t658;
t788 = t520 * t433 + t521 * t434 + t502 * t822 + t821 * t503 - t820 * t673;
t244 = Icges(4,4) * t375 + Icges(4,2) * t374 + Icges(4,6) * t461;
t246 = Icges(4,1) * t375 + Icges(4,4) * t374 + Icges(4,5) * t461;
t380 = -Icges(4,4) * t583 + Icges(4,2) * t486 - Icges(4,6) * t574;
t382 = -Icges(4,1) * t583 + Icges(4,4) * t486 - Icges(4,5) * t574;
t806 = -t244 * t488 - t246 * t489 - t318 * t815 + t319 * t813 - t372 * t380 - t373 * t382 + t459 * t809 + t479 * t818 + t480 * t816 + t576 * t810;
t805 = -t244 * t486 + t246 * t583 - t320 * t815 + t321 * t813 - t374 * t380 - t375 * t382 - t461 * t809 + t477 * t818 + t478 * t816 + t574 * t810;
t243 = Icges(4,4) * t373 + Icges(4,2) * t372 - Icges(4,6) * t459;
t245 = Icges(4,1) * t373 + Icges(4,4) * t372 - Icges(4,5) * t459;
t381 = Icges(4,4) * t489 + Icges(4,2) * t488 - Icges(4,6) * t576;
t383 = Icges(4,1) * t489 + Icges(4,4) * t488 - Icges(4,5) * t576;
t804 = -t243 * t488 - t245 * t489 + t814 * t318 + t812 * t319 - t372 * t381 - t373 * t383 + t808 * t459 + t479 * t819 + t817 * t480 + t811 * t576;
t803 = -t243 * t486 + t245 * t583 + t814 * t320 + t812 * t321 - t374 * t381 - t375 * t383 - t808 * t461 + t477 * t819 + t817 * t478 + t811 * t574;
t78 = t502 * t209 - t503 * t213 + t463 * t328 - t464 * t332 + (-t217 * t756 + t336 * t709) * t550;
t80 = -t502 * t215 + t503 * t219 - t463 * t334 + t464 * t338 + (-t211 * t756 + t330 * t709) * t550;
t94 = t520 * t243 + t521 * t245 + t484 * t381 + t485 * t383 + (-t241 * t756 + t379 * t709) * t550;
t802 = t78 + t80 + t94;
t104 = t318 * t417 - t319 * t418 + t357 * t479 - t358 * t480 - t359 * t576 - t419 * t459;
t106 = -t318 * t421 + t319 * t422 - t360 * t576 - t361 * t479 + t362 * t480 - t420 * t459;
t113 = t372 * t433 + t373 * t434 - t398 * t576 + t399 * t488 + t400 * t489 - t432 * t459;
t801 = t113 + t106 + t104;
t105 = t320 * t417 - t321 * t418 + t357 * t477 - t358 * t478 - t359 * t574 + t419 * t461;
t107 = -t320 * t421 + t321 * t422 - t360 * t574 - t361 * t477 + t362 * t478 + t420 * t461;
t114 = t374 * t433 + t375 * t434 - t398 * t574 + t399 * t486 - t400 * t583 + t432 * t461;
t800 = t114 + t107 + t105;
t798 = t380 * t486 - t382 * t583 + t477 * t815 - t478 * t813 - t574 * t809;
t797 = t380 * t488 + t382 * t489 + t479 * t815 - t480 * t813 - t576 * t809;
t174 = t502 * t327 - t503 * t331 - t335 * t673;
t176 = -t329 * t673 - t502 * t333 + t503 * t337;
t187 = -t378 * t673 + t520 * t380 + t521 * t382;
t790 = t174 + t176 + t187;
t175 = t502 * t328 - t503 * t332 - t336 * t673;
t177 = -t330 * t673 - t502 * t334 + t503 * t338;
t188 = -t379 * t673 + t520 * t381 + t521 * t383;
t789 = t175 + t177 + t188;
t796 = -t381 * t486 + t383 * t583 + t477 * t814 + t478 * t812 + t574 * t808;
t795 = -t381 * t488 - t383 * t489 + t479 * t814 + t480 * t812 + t576 * t808;
t192 = t417 * t477 - t418 * t478 - t419 * t574;
t194 = -t420 * t574 - t421 * t477 + t422 * t478;
t239 = -t432 * t574 + t433 * t486 - t434 * t583;
t794 = t239 + t194 + t192;
t193 = t417 * t479 - t418 * t480 - t419 * t576;
t195 = -t420 * t576 - t421 * t479 + t422 * t480;
t240 = -t432 * t576 + t433 * t488 + t434 * t489;
t793 = t240 + t195 + t193;
t792 = t807 * t551;
t791 = t788 * t658;
t553 = sin(qJ(6));
t557 = cos(qJ(6));
t414 = t479 * t553 - t557 * t576;
t253 = -qJD(6) * t414 + t318 * t557 + t459 * t553;
t413 = t479 * t557 + t553 * t576;
t254 = qJD(6) * t413 + t318 * t553 - t459 * t557;
t135 = t254 * rSges(7,1) + t253 * rSges(7,2) + t319 * rSges(7,3);
t740 = -t459 * pkin(5) + pkin(10) * t319 + t135;
t412 = t477 * t553 - t557 * t574;
t255 = -qJD(6) * t412 + t320 * t557 - t461 * t553;
t411 = t477 * t557 + t553 * t574;
t256 = qJD(6) * t411 + t320 * t553 + t461 * t557;
t604 = -t256 * rSges(7,1) - t255 * rSges(7,2);
t136 = t321 * rSges(7,3) - t604;
t739 = t461 * pkin(5) + t321 * pkin(10) + t136;
t584 = -t502 * t553 + t557 * t673;
t348 = qJD(6) * t584 + t463 * t557 - t553 * t658;
t475 = t502 * t557 + t553 * t673;
t349 = qJD(6) * t475 + t463 * t553 + t557 * t658;
t236 = rSges(7,1) * t349 + rSges(7,2) * t348 + rSges(7,3) * t464;
t734 = pkin(5) * t658 + pkin(10) * t464 + t236;
t603 = -rSges(7,1) * t412 - rSges(7,2) * t411;
t272 = rSges(7,3) * t478 - t603;
t732 = -pkin(5) * t574 + pkin(10) * t478 + t272;
t273 = t414 * rSges(7,1) + t413 * rSges(7,2) + t480 * rSges(7,3);
t731 = -pkin(5) * t576 + pkin(10) * t480 + t273;
t307 = -rSges(7,1) * t584 + rSges(7,2) * t475 + rSges(7,3) * t503;
t727 = -pkin(5) * t673 + t503 * pkin(10) + t307;
t266 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t478;
t268 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t478;
t270 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t478;
t118 = t266 * t480 + t268 * t413 + t270 * t414;
t267 = Icges(7,5) * t414 + Icges(7,6) * t413 + Icges(7,3) * t480;
t269 = Icges(7,4) * t414 + Icges(7,2) * t413 + Icges(7,6) * t480;
t271 = Icges(7,1) * t414 + Icges(7,4) * t413 + Icges(7,5) * t480;
t119 = t267 * t480 + t269 * t413 + t271 * t414;
t302 = -Icges(7,5) * t584 + Icges(7,6) * t475 + Icges(7,3) * t503;
t303 = -Icges(7,4) * t584 + Icges(7,2) * t475 + Icges(7,6) * t503;
t304 = -Icges(7,1) * t584 + Icges(7,4) * t475 + Icges(7,5) * t503;
t149 = t302 * t480 + t303 * t413 + t304 * t414;
t130 = Icges(7,5) * t256 + Icges(7,6) * t255 + Icges(7,3) * t321;
t132 = Icges(7,4) * t256 + Icges(7,2) * t255 + Icges(7,6) * t321;
t134 = Icges(7,1) * t256 + Icges(7,4) * t255 + Icges(7,5) * t321;
t30 = t130 * t480 + t132 * t413 + t134 * t414 + t253 * t268 + t254 * t270 + t266 * t319;
t129 = Icges(7,5) * t254 + Icges(7,6) * t253 + Icges(7,3) * t319;
t131 = Icges(7,4) * t254 + Icges(7,2) * t253 + Icges(7,6) * t319;
t133 = Icges(7,1) * t254 + Icges(7,4) * t253 + Icges(7,5) * t319;
t31 = t129 * t480 + t131 * t413 + t133 * t414 + t253 * t269 + t254 * t271 + t267 * t319;
t233 = Icges(7,5) * t349 + Icges(7,6) * t348 + Icges(7,3) * t464;
t234 = Icges(7,4) * t349 + Icges(7,2) * t348 + Icges(7,6) * t464;
t235 = Icges(7,1) * t349 + Icges(7,4) * t348 + Icges(7,5) * t464;
t49 = t233 * t480 + t234 * t413 + t235 * t414 + t253 * t303 + t254 * t304 + t302 * t319;
t3 = t118 * t461 - t119 * t459 - t30 * t574 - t31 * t576 + (t149 * t709 - t49 * t756) * t550;
t787 = t3 + t804 * t576 + t806 * t574 + (t709 * t793 - t756 * t801) * t550 + t797 * t461 + t795 * t459;
t116 = t266 * t478 + t268 * t411 + t270 * t412;
t117 = t267 * t478 + t269 * t411 + t271 * t412;
t148 = t302 * t478 + t303 * t411 + t304 * t412;
t32 = t130 * t478 + t132 * t411 + t134 * t412 + t255 * t268 + t256 * t270 + t266 * t321;
t33 = t129 * t478 + t131 * t411 + t133 * t412 + t255 * t269 + t256 * t271 + t267 * t321;
t50 = t233 * t478 + t234 * t411 + t235 * t412 + t255 * t303 + t256 * t304 + t302 * t321;
t4 = t116 * t461 - t117 * t459 - t32 * t574 - t33 * t576 + (t148 * t709 - t50 * t756) * t550;
t786 = t4 + t803 * t576 + t805 * t574 + (t709 * t794 - t756 * t800) * t550 + t798 * t461 + t796 * t459;
t5 = t49 * t551 + (-t757 * t30 + t31 * t556 + (t118 * t556 + t119 * t757) * qJD(1)) * t550;
t785 = t5 + t801 * t551 + (t806 * t757 - t804 * t556 + (t556 * t797 - t757 * t795) * qJD(1)) * t550;
t6 = t50 * t551 + (-t757 * t32 + t33 * t556 + (t116 * t556 + t117 * t757) * qJD(1)) * t550;
t784 = t6 + t800 * t551 + (t805 * t757 - t803 * t556 + (t556 * t798 - t757 * t796) * qJD(1)) * t550;
t77 = t502 * t210 - t503 * t214 + t463 * t327 - t464 * t331 + (-t218 * t756 + t335 * t709) * t550;
t79 = -t502 * t216 + t503 * t220 - t463 * t333 + t464 * t337 + (-t212 * t756 + t329 * t709) * t550;
t127 = t266 * t503 + t268 * t475 - t270 * t584;
t128 = t267 * t503 + t269 * t475 - t271 * t584;
t158 = t302 * t503 + t303 * t475 - t304 * t584;
t155 = t158 * t658;
t36 = t130 * t503 + t132 * t475 - t134 * t584 + t266 * t464 + t268 * t348 + t270 * t349;
t37 = t129 * t503 + t131 * t475 - t133 * t584 + t267 * t464 + t269 * t348 + t271 * t349;
t70 = t503 * t233 + t475 * t234 - t235 * t584 + t464 * t302 + t348 * t303 + t349 * t304;
t8 = t127 * t461 - t128 * t459 - t36 * t574 - t37 * t576 - t673 * t70 + t155;
t93 = t520 * t244 + t521 * t246 + t484 * t380 + t485 * t382 + (-t242 * t756 + t378 * t709) * t550;
t783 = t791 + t8 - t807 * t673 - t802 * t576 + (-t77 - t79 - t93) * t574 + t790 * t461 - t789 * t459;
t667 = t757 * t188;
t668 = t757 * t177;
t669 = t757 * t175;
t689 = t757 * t93;
t690 = t757 * t79;
t691 = t757 * t77;
t671 = t757 * t128;
t69 = t70 * t551;
t692 = t757 * t36;
t9 = t69 + (-t692 + t37 * t556 + (t127 * t556 + t671) * qJD(1)) * t550;
t782 = t9 + t792 + (-t689 - t690 - t691 + t802 * t556 + (t556 * t790 + t667 + t668 + t669) * qJD(1)) * t550;
t44 = -t116 * t574 - t117 * t576 - t148 * t673;
t781 = -t574 * t798 + t576 * t796 - t673 * t794 + t44;
t45 = -t118 * t574 - t119 * t576 - t149 * t673;
t780 = -t574 * t797 + t576 * t795 - t673 * t793 + t45;
t46 = t148 * t551 + (-t116 * t757 + t117 * t556) * t550;
t779 = t46 + t794 * t551 + (-t556 * t796 - t757 * t798) * t550;
t47 = t149 * t551 + (-t118 * t757 + t119 * t556) * t550;
t778 = t47 + t793 * t551 + (-t556 * t795 - t757 * t797) * t550;
t777 = -t70 - t807;
t776 = t550 ^ 2;
t775 = -0.2e1 * t459;
t774 = 0.2e1 * t461;
t773 = -0.2e1 * t574;
t772 = -0.2e1 * t576;
t771 = m(7) / 0.2e1;
t770 = t319 / 0.2e1;
t769 = t321 / 0.2e1;
t768 = -t459 / 0.2e1;
t767 = t461 / 0.2e1;
t766 = t464 / 0.2e1;
t765 = t478 / 0.2e1;
t764 = t480 / 0.2e1;
t763 = t503 / 0.2e1;
t762 = -t574 / 0.2e1;
t761 = -t576 / 0.2e1;
t760 = t551 / 0.2e1;
t759 = t556 / 0.2e1;
t758 = rSges(6,2) - pkin(4);
t548 = pkin(3) * t558 + pkin(2);
t755 = -pkin(2) + t548;
t552 = -qJ(4) - pkin(9);
t754 = -pkin(5) + t552;
t753 = rSges(6,3) * t477;
t752 = pkin(3) * qJD(3);
t751 = t320 * rSges(6,3);
t750 = -rSges(6,1) + t552;
t749 = -rSges(5,3) + t552;
t748 = t158 * t464 + t70 * t503;
t747 = Icges(3,4) * t555;
t458 = t461 * pkin(9);
t695 = t554 * t752;
t675 = -qJD(4) * t574 + t575 * t695 - t631 * t752;
t230 = pkin(3) * t630 - t461 * t552 + t462 * t755 - t458 + t675;
t517 = t574 * pkin(9);
t539 = pkin(3) * t632;
t346 = t552 * t574 - t575 * t755 + t517 - t539;
t738 = -t230 * t576 - t459 * t346;
t391 = t460 * pkin(2) - pkin(9) * t459;
t567 = pkin(3) * t596 - qJD(4) * t576 + t459 * t552 + t460 * t548 - t525 * t695 + t688 * t752;
t229 = -t391 + t567;
t377 = t551 * t391;
t737 = t551 * t229 + t377;
t227 = t319 * pkin(4) + qJ(5) * t318 + qJD(5) * t479;
t736 = -t227 - t229;
t392 = t462 * pkin(2) + t458;
t735 = -t230 - t392;
t248 = t375 * rSges(4,1) + t374 * rSges(4,2) + t461 * rSges(4,3);
t733 = -t248 - t392;
t296 = t576 * t346;
t468 = t477 * qJ(5);
t402 = pkin(4) * t478 + t468;
t730 = -t402 * t576 - t296;
t301 = pkin(4) * t464 + qJ(5) * t463 + qJD(5) * t502;
t697 = t756 * pkin(2);
t404 = t742 * t752 + (-t555 * t695 - t756 * qJD(4) + (-t697 + t756 * t548 + (-pkin(9) - t552) * t555) * qJD(2)) * t550;
t729 = -t301 - t404;
t467 = t525 * pkin(2) - pkin(9) * t576;
t676 = pkin(3) * t687 + t525 * t548 + t552 * t576;
t347 = -t467 + t676;
t305 = t347 * t658;
t403 = t480 * pkin(4) + qJ(5) * t479;
t728 = t403 * t658 + t305;
t696 = t756 * pkin(9);
t435 = pkin(3) * t743 + (t552 * t756 + t555 * t755 + t696) * t550;
t726 = t346 * t673 - t435 * t574;
t725 = -t320 * qJ(5) - t477 * qJD(5);
t450 = t551 * t467;
t724 = t551 * t347 + t450;
t605 = -rSges(5,1) * t478 + rSges(5,2) * t477;
t343 = -rSges(5,3) * t574 - t605;
t723 = -t343 - t346;
t344 = t480 * rSges(5,1) - t479 * rSges(5,2) - rSges(5,3) * t576;
t722 = -t344 - t347;
t721 = -t346 - t402;
t720 = -t347 - t403;
t366 = rSges(5,1) * t464 - rSges(5,2) * t463 + rSges(5,3) * t658;
t719 = -t366 - t404;
t384 = -rSges(4,1) * t583 + rSges(4,2) * t486 - rSges(4,3) * t574;
t466 = -pkin(2) * t575 - t517;
t718 = -t384 - t466;
t385 = t489 * rSges(4,1) + t488 * rSges(4,2) - rSges(4,3) * t576;
t717 = -t385 - t467;
t425 = t503 * rSges(5,1) - t502 * rSges(5,2) - rSges(5,3) * t673;
t716 = t425 + t435;
t526 = (pkin(2) * t555 - t696) * t550;
t493 = t526 * t659;
t715 = t435 * t659 + t493;
t431 = pkin(4) * t503 + qJ(5) * t502;
t714 = t431 + t435;
t713 = t466 * t745 + t467 * t674;
t549 = t757 * pkin(1);
t712 = pkin(8) * t745 + t549;
t710 = qJD(2) * t550;
t705 = -rSges(7,3) - pkin(4) - pkin(10);
t704 = 0.2e1 * t550;
t703 = m(6) / 0.2e1 + t771;
t702 = -t757 / 0.2e1;
t701 = t757 / 0.2e1;
t700 = m(5) * t756;
t699 = m(6) * t756;
t698 = m(7) * t756;
t694 = t36 / 0.2e1 + t50 / 0.2e1;
t693 = t37 / 0.2e1 + t49 / 0.2e1;
t686 = t230 * t673 - t404 * t574 + t461 * t435;
t685 = t551 * t227 + t737;
t228 = t321 * pkin(4) - t725;
t684 = -t228 + t735;
t365 = rSges(6,1) * t658 - rSges(6,2) * t464 + rSges(6,3) * t463;
t683 = -t365 + t729;
t223 = t319 * rSges(5,1) - t318 * rSges(5,2) - t459 * rSges(5,3);
t682 = t551 * t403 + t724;
t341 = -rSges(6,1) * t574 - rSges(6,2) * t478 + t753;
t681 = -t341 + t721;
t342 = -rSges(6,1) * t576 - t480 * rSges(6,2) + t479 * rSges(6,3);
t680 = -t342 + t720;
t247 = t373 * rSges(4,1) + t372 * rSges(4,2) - t459 * rSges(4,3);
t679 = t391 * t674 + t392 * t745 + t466 * t628;
t678 = t431 * t659 + t715;
t424 = -rSges(6,1) * t673 - t503 * rSges(6,2) + t502 * rSges(6,3);
t677 = t424 + t714;
t221 = -t459 * rSges(6,1) - t319 * rSges(6,2) + t318 * rSges(6,3);
t363 = t460 * rSges(3,1) + t459 * rSges(3,2) + rSges(3,3) * t628;
t444 = t525 * rSges(3,1) + rSges(3,2) * t576 + rSges(3,3) * t745;
t351 = Icges(3,5) * t462 - Icges(3,6) * t461 + Icges(3,3) * t659;
t353 = Icges(3,4) * t462 - Icges(3,2) * t461 + Icges(3,6) * t659;
t355 = Icges(3,1) * t462 - Icges(3,4) * t461 + Icges(3,5) * t659;
t439 = -Icges(3,4) * t575 + Icges(3,2) * t574 - Icges(3,6) * t674;
t441 = -Icges(3,1) * t575 + Icges(3,4) * t574 - Icges(3,5) * t674;
t670 = t757 * (t551 * t351 + (t756 * t353 + t355 * t555 + (-t439 * t555 + t441 * t756) * qJD(2)) * t550);
t438 = Icges(3,5) * t525 + Icges(3,6) * t576 + Icges(3,3) * t745;
t440 = Icges(3,4) * t525 + Icges(3,2) * t576 + Icges(3,6) * t745;
t442 = Icges(3,1) * t525 + Icges(3,4) * t576 + Icges(3,5) * t745;
t666 = t757 * (t551 * t438 + (t440 * t756 + t442 * t555) * t550);
t510 = (pkin(9) * t555 + t697) * t710;
t665 = t757 * t510;
t664 = t757 * t526;
t662 = t756 * Icges(3,4);
t661 = t756 * t229;
t660 = t756 * t347;
t656 = t127 / 0.2e1 + t148 / 0.2e1;
t655 = t128 / 0.2e1 + t149 / 0.2e1;
t652 = t709 / 0.2e1;
t651 = -t556 * pkin(1) + pkin(8) * t674;
t650 = 2 * m(3);
t648 = 2 * m(4);
t646 = 0.2e1 * m(5);
t644 = 0.2e1 * m(6);
t642 = 0.2e1 * m(7);
t639 = -t228 * t576 - t459 * t402 + t738;
t638 = t729 - t734;
t637 = t721 - t732;
t636 = t720 - t731;
t635 = t714 + t727;
t634 = t346 * t745 + t347 * t674 + t713;
t633 = t402 * t673 - t431 * t574 + t726;
t579 = -t227 * t756 - t661;
t34 = -t638 * t576 + t635 * t459 + (t731 * t709 - t740 * t756 + t579) * t550 + t728;
t54 = -t683 * t576 + t677 * t459 + (-t221 * t756 + t342 * t709 + t579) * t550 + t728;
t622 = m(6) * t54 + m(7) * t34;
t591 = t228 * t673 - t301 * t574 + t461 * t431 + t686;
t35 = -t734 * t574 + t727 * t461 + (t637 * t709 + t739 * t756) * t550 + t591;
t222 = t461 * rSges(6,1) - t321 * rSges(6,2) + t751;
t55 = -t574 * t365 + t461 * t424 + (t222 * t756 + t681 * t709) * t550 + t591;
t621 = m(6) * t55 + m(7) * t35;
t582 = -t404 * t757 - t665;
t570 = -t301 * t757 + t582;
t40 = (t684 - t739) * t551 + (t727 * t711 - t734 * t757 + t570) * t550 + t678;
t67 = (-t222 + t684) * t551 + (-t365 * t757 + t424 * t711 + t570) * t550 + t678;
t620 = m(6) * t67 + m(7) * t40;
t581 = -t435 * t757 - t664;
t568 = -t431 * t757 + t581;
t561 = -t727 * t757 + t568;
t41 = t740 * t551 + ((-t510 + t638) * t556 + t561 * qJD(1)) * t550 + t685;
t564 = -t424 * t757 + t568;
t68 = t551 * t221 + ((-t510 + t683) * t556 + t564 * qJD(1)) * t550 + t685;
t619 = m(6) * t68 + m(7) * t41;
t589 = t402 * t745 + t403 * t674 + t634;
t115 = (t341 * t556 + t342 * t757) * t550 + t589;
t99 = (t732 * t556 + t731 * t757) * t550 + t589;
t618 = m(6) * t115 + m(7) * t99;
t615 = -pkin(1) * t711 + pkin(8) * t628;
t563 = t567 + t615;
t559 = t227 + t563;
t120 = t559 + t221;
t83 = t559 + t740;
t617 = m(6) * t120 + m(7) * t83;
t562 = (-t549 + (-pkin(3) * t554 - pkin(8)) * t745) * qJD(1) - t462 * t548 - t675;
t560 = t562 + t725;
t121 = t321 * t758 + t461 * t750 + t560 - t751;
t84 = t321 * t705 + t461 * t754 + t560 + t604;
t616 = m(6) * t121 + m(7) * t84;
t100 = -t574 * t636 - t576 * t732 + t730;
t122 = -t341 * t576 - t574 * t680 + t730;
t614 = m(6) * t122 + m(7) * t100;
t108 = t731 * t551 + (-t526 - t635) * t745 + t682;
t137 = t342 * t551 + (-t526 - t677) * t745 + t682;
t613 = m(6) * t137 + m(7) * t108;
t109 = (-t466 + t637) * t551 + t561 * t550;
t138 = (-t466 + t681) * t551 + t564 * t550;
t612 = m(6) * t138 + m(7) * t109;
t111 = -t574 * t727 + t673 * t732 + t633;
t146 = t341 * t673 - t424 * t574 + t633;
t611 = m(6) * t146 + m(7) * t111;
t578 = -t403 * t756 - t660;
t112 = (-t731 * t756 + t578) * t550 + t635 * t576;
t147 = (-t342 * t756 + t578) * t550 + t677 * t576;
t610 = m(6) * t147 + m(7) * t112;
t587 = t548 * t575 + t539 + t651;
t577 = -t468 + t587;
t170 = t478 * t705 - t574 * t754 + t577 + t603;
t225 = t478 * t758 - t574 * t750 + t577 - t753;
t609 = m(6) * t225 + m(7) * t170;
t600 = t676 + t712;
t573 = t403 + t600;
t171 = t573 + t731;
t226 = t573 + t342;
t608 = m(6) * t226 + m(7) * t171;
t607 = -t462 * rSges(3,1) + t461 * rSges(3,2);
t606 = -t321 * rSges(5,1) + t320 * rSges(5,2);
t590 = t229 * t674 + t230 * t745 + t346 * t628 + t679;
t436 = t521 * rSges(4,1) + t520 * rSges(4,2) - rSges(4,3) * t673;
t580 = -t436 * t757 - t664;
t500 = Icges(3,6) * t551 + (Icges(3,2) * t756 + t747) * t550;
t501 = Icges(3,5) * t551 + (Icges(3,1) * t555 + t662) * t550;
t506 = (Icges(3,5) * t756 - Icges(3,6) * t555) * t710;
t507 = (-Icges(3,2) * t555 + t662) * t710;
t508 = (Icges(3,1) * t756 - t747) * t710;
t572 = -t500 * t658 + t501 * t627 + t551 * t506 + t507 * t673 + t508 * t746;
t571 = t227 * t674 + t228 * t745 + t402 * t628 + t590;
t443 = -rSges(3,1) * t575 + rSges(3,2) * t574 - rSges(3,3) * t674;
t569 = -t425 * t757 + t581;
t566 = t113 / 0.2e1 + t106 / 0.2e1 + t104 / 0.2e1 + t94 / 0.2e1 + t80 / 0.2e1 + t78 / 0.2e1 + t693;
t565 = t187 / 0.2e1 + t176 / 0.2e1 + t174 / 0.2e1 + t239 / 0.2e1 + t192 / 0.2e1 + t194 / 0.2e1 + t656;
t509 = (rSges(3,1) * t756 - rSges(3,2) * t555) * t710;
t504 = t551 * rSges(3,3) + (rSges(3,1) * t555 + rSges(3,2) * t756) * t550;
t499 = Icges(3,3) * t551 + (Icges(3,5) * t555 + Icges(3,6) * t756) * t550;
t437 = -Icges(3,5) * t575 + Icges(3,6) * t574 - Icges(3,3) * t674;
t416 = t444 + t712;
t415 = -t443 + t651;
t401 = rSges(4,1) * t485 + rSges(4,2) * t484 + rSges(4,3) * t658;
t397 = -t551 * t443 - t504 * t674;
t396 = t444 * t551 - t504 * t745;
t364 = rSges(3,3) * t659 - t607;
t354 = Icges(3,1) * t460 + Icges(3,4) * t459 + Icges(3,5) * t628;
t352 = Icges(3,4) * t460 + Icges(3,2) * t459 + Icges(3,6) * t628;
t350 = Icges(3,5) * t460 + Icges(3,6) * t459 + Icges(3,3) * t628;
t300 = (-t549 + (-rSges(3,3) - pkin(8)) * t745) * qJD(1) + t607;
t299 = t615 + t363;
t298 = t499 * t745 + t500 * t576 + t501 * t525;
t297 = -t499 * t674 + t500 * t574 - t501 * t575;
t287 = t712 - t717;
t286 = t651 + t718;
t282 = t572 * t551;
t281 = t551 * t363 + (-t504 * t654 - t509 * t556) * t550;
t280 = -t551 * t364 + (t504 * t711 - t509 * t757) * t550;
t279 = -t385 * t673 + t436 * t576;
t278 = t384 * t673 - t436 * t574;
t276 = t551 * t437 + (t439 * t756 + t441 * t555) * t550;
t275 = t600 + t344;
t274 = -t574 * t749 + t587 + t605;
t265 = t438 * t745 + t440 * t576 + t442 * t525;
t264 = t437 * t745 + t439 * t576 + t441 * t525;
t263 = -t438 * t674 + t440 * t574 - t442 * t575;
t262 = -t437 * t674 + t439 * t574 - t441 * t575;
t257 = -t384 * t576 + t385 * t574;
t252 = t550 * t580 + t551 * t718;
t251 = t385 * t551 + t450 + (-t436 - t526) * t745;
t224 = t461 * rSges(5,3) - t606;
t197 = -t461 * t500 + t462 * t501 + t574 * t507 - t575 * t508 + (t499 * t711 - t506 * t757) * t550;
t196 = t459 * t500 + t460 * t501 + t576 * t507 + t525 * t508 + (t499 * t654 + t506 * t556) * t550;
t189 = (t384 * t556 + t385 * t757) * t550 + t713;
t186 = -qJD(1) * t712 + t733;
t185 = t391 + t615 + t247;
t184 = t273 * t503 - t307 * t480;
t183 = -t272 * t503 + t307 * t478;
t173 = (-t344 * t756 - t660) * t550 + t716 * t576;
t172 = t343 * t673 - t425 * t574 + t726;
t161 = (-t466 + t723) * t551 + t569 * t550;
t160 = t344 * t551 + (-t526 - t716) * t745 + t724;
t159 = t272 * t480 - t273 * t478;
t157 = t551 * t350 + (t756 * t352 + t354 * t555 + (-t440 * t555 + t442 * t756) * qJD(2)) * t550;
t154 = t461 * t749 + t562 + t606;
t153 = t563 + t223;
t151 = -t343 * t576 - t574 * t722 - t296;
t144 = (t343 * t556 + t344 * t757) * t550 + t634;
t142 = t551 * t247 + t377 + ((-t401 - t510) * t556 + t580 * qJD(1)) * t550;
t141 = t493 + t733 * t551 + (-t401 * t757 + t436 * t711 - t665) * t550;
t140 = -t574 * t401 + t461 * t436 + (t248 * t756 - t384 * t709) * t550;
t139 = t576 * t401 + t459 * t436 + (-t247 * t756 + t385 * t709) * t550;
t110 = t247 * t574 - t248 * t576 - t384 * t459 - t385 * t461;
t103 = (t757 * t247 + t248 * t556 + (t384 * t757 + t556 * t717) * qJD(1)) * t550 + t679;
t98 = t551 * t223 + ((-t510 + t719) * t556 + t569 * qJD(1)) * t550 + t737;
t97 = (-t224 + t735) * t551 + (-t366 * t757 + t425 * t711 + t582) * t550 + t715;
t82 = -t574 * t366 + t461 * t425 + (t224 * t756 + t709 * t723) * t550 + t686;
t81 = t305 - t719 * t576 + t716 * t459 + (-t223 * t756 + t344 * t709 - t661) * t550;
t76 = -t136 * t503 + t236 * t478 - t272 * t464 + t307 * t321;
t75 = t135 * t503 - t236 * t480 + t273 * t464 - t307 * t319;
t57 = -t224 * t576 - t343 * t459 - (-t223 - t229) * t574 + t722 * t461 + t738;
t56 = t158 * t551 + (-t127 * t757 + t128 * t556) * t550;
t53 = -t127 * t574 - t128 * t576 - t158 * t673;
t52 = -t135 * t478 + t136 * t480 + t272 * t319 - t273 * t321;
t51 = (t757 * t223 + t224 * t556 + (t757 * t343 + (-t467 + t722) * t556) * qJD(1)) * t550 + t590;
t48 = t127 * t478 + t128 * t480 + t158 * t503;
t43 = t118 * t478 + t119 * t480 + t149 * t503;
t42 = t116 * t478 + t117 * t480 + t148 * t503;
t39 = -t222 * t576 - t341 * t459 - (-t221 + t736) * t574 + t680 * t461 + t639;
t38 = (t757 * t221 + t222 * t556 + (t757 * t341 + (-t467 + t680) * t556) * qJD(1)) * t550 + t571;
t28 = t571 + (t739 + (-t467 + t636) * qJD(1)) * t745 + (qJD(1) * t732 + t740) * t674;
t26 = -t739 * t576 - t732 * t459 - (t736 - t740) * t574 + t636 * t461 + t639;
t7 = t127 * t321 + t128 * t319 + t36 * t478 + t37 * t480 + t748;
t2 = t116 * t321 + t117 * t319 + t148 * t464 + t32 * t478 + t33 * t480 + t50 * t503;
t1 = t118 * t321 + t119 * t319 + t149 * t464 + t30 * t478 + t31 * t480 + t49 * t503;
t10 = [(t120 * t226 + t121 * t225) * t644 + (t153 * t275 + t154 * t274) * t646 + (t170 * t84 + t171 * t83) * t642 + (t185 * t287 + t186 * t286) * t648 + (t299 * t416 + t300 * t415) * t650 + t572 - t777; m(3) * (t280 * t415 + t281 * t416 + t299 * t396 + t300 * t397) + (t141 * t286 + t142 * t287 + t185 * t251 + t186 * t252) * m(4) + (t153 * t160 + t154 * t161 + t274 * t97 + t275 * t98) * m(5) + (t120 * t137 + t121 * t138 + t225 * t67 + t226 * t68) * m(6) + t282 + (t108 * t83 + t109 * t84 + t170 * t40 + t171 * t41) * m(7) + t69 + t792 + (-t691 / 0.2e1 - t692 / 0.2e1 - t689 / 0.2e1 - t670 / 0.2e1 - t690 / 0.2e1 + (t196 / 0.2e1 + t157 / 0.2e1 + t566) * t556 + (t667 / 0.2e1 + t671 / 0.2e1 + t668 / 0.2e1 + t669 / 0.2e1 + t666 / 0.2e1 + (t297 / 0.2e1 + t276 / 0.2e1 + t565) * t556 + (t298 + t149 + t793) * t701) * qJD(1) + (t50 + t197 + t800) * t702) * t550; (t397 * t280 + t396 * t281 + (t443 * t556 + t444 * t757) * (t757 * t363 + t364 * t556 + (t443 * t757 - t444 * t556) * qJD(1)) * t776) * t650 + (t108 * t41 + t109 * t40 + t28 * t99) * t642 + (t115 * t38 + t137 * t68 + t138 * t67) * t644 + (t144 * t51 + t160 * t98 + t161 * t97) * t646 + (t103 * t189 + t141 * t252 + t142 * t251) * t648 + ((t576 * t352 + t525 * t354 + t459 * t440 + t460 * t442 + (t350 * t556 + t438 * t654) * t550) * t745 + t785) * t745 + ((t574 * t353 - t575 * t355 - t461 * t439 + t462 * t441 + (-t351 * t757 + t437 * t711) * t550) * t674 + (-t574 * t352 - t576 * t353 + t575 * t354 - t525 * t355 - t459 * t439 + t461 * t440 - t460 * t441 - t462 * t442 + (t350 * t757 - t351 * t556 - t437 * t654 - t438 * t711) * t550) * t745 - t784) * t674 + (-t262 * t674 + t264 * t745 + (-t262 * t757 + t263 * t556) * t550 + t779) * t659 + (-t263 * t674 + t265 * t745 + (-t264 * t757 + t265 * t556) * t550 + t778) * t628 + (t298 * t628 - t197 * t674 + t297 * t659 + t196 * t745 + t282 + (-t670 + t157 * t556 + (t276 * t556 + t666) * qJD(1)) * t550 + t782) * t551; (t139 * t287 + t140 * t286 + t185 * t279 + t186 * t278) * m(4) + (t153 * t173 + t154 * t172 + t274 * t82 + t275 * t81) * m(5) + (t120 * t147 + t121 * t146 + t225 * t55 + t226 * t54) * m(6) + (t111 * t84 + t112 * t83 + t170 * t35 + t171 * t34) * m(7) + t155 - t566 * t576 - (t114 / 0.2e1 + t107 / 0.2e1 + t105 / 0.2e1 + t93 / 0.2e1 + t79 / 0.2e1 + t77 / 0.2e1 + t694) * t574 + t565 * t461 + (-t240 / 0.2e1 - t193 / 0.2e1 - t195 / 0.2e1 - t188 / 0.2e1 - t177 / 0.2e1 - t175 / 0.2e1 - t655) * t459 + t777 * t673 + t791; (t100 * t28 + t108 * t34 + t109 * t35 + t111 * t40 + t112 * t41 + t26 * t99) * m(7) + (t115 * t39 + t122 * t38 + t137 * t54 + t138 * t55 + t146 * t67 + t147 * t68) * m(6) + (t144 * t57 + t151 * t51 + t160 * t81 + t161 * t82 + t172 * t97 + t173 * t98) * m(5) + (t103 * t257 + t110 * t189 + t139 * t251 + t140 * t252 + t141 * t278 + t142 * t279) * m(4) + t778 * t768 + t779 * t767 + t784 * t762 + t785 * t761 + t783 * t760 - t786 * t674 / 0.2e1 - t782 * t673 / 0.2e1 + (t56 + t788 * t551 + (t556 * t789 - t757 * t790) * t550) * t550 * t652 + t780 * t628 / 0.2e1 + (qJD(1) * t781 + t787) * t745 / 0.2e1; (t110 * t257 + t139 * t279 + t140 * t278) * t648 + (t151 * t57 + t172 * t82 + t173 * t81) * t646 + (t122 * t39 + t146 * t55 + t147 * t54) * t644 + (t100 * t26 + t111 * t35 + t112 * t34) * t642 - t783 * t673 - t787 * t576 - t786 * t574 + t781 * t461 - t780 * t459 + (-t574 * t790 - t576 * t789 - t673 * t788 + t53) * t658; -(m(5) * t154 + t616) * t576 - (m(5) * t153 + t617) * t574 + (m(5) * t275 + t608) * t461 + (-m(5) * t274 - t609) * t459; (-t28 * t698 - t38 * t699 - t51 * t700) * t550 - (m(5) * t97 + t620) * t576 - (m(5) * t98 + t619) * t574 + (m(5) * t160 + t613) * t461 + (-m(5) * t161 - t612) * t459 + (m(5) * t144 + t618) * t658; (-t26 * t698 - t39 * t699 - t57 * t700) * t550 - (m(5) * t82 + t621) * t576 - (m(5) * t81 + t622) * t574 + (m(5) * t173 + t610) * t461 + (-m(5) * t172 - t611) * t459 + (m(5) * t151 + t614) * t658; 0.4e1 * (m(5) / 0.2e1 + t703) * (-t555 * t653 * t776 + t459 * t576 - t461 * t574); t318 * t609 + t320 * t608 + t477 * t617 + t479 * t616; (m(6) * t38 + m(7) * t28) * t502 + t620 * t479 + t619 * t477 + t618 * t463 + t613 * t320 + t612 * t318; (m(6) * t39 + m(7) * t26) * t502 + t621 * t479 + t622 * t477 + t614 * t463 + t610 * t320 + t611 * t318; t703 * (t318 * t772 + t320 * t773 + t479 * t775 + t477 * t774 + (-t463 * t756 + t502 * t709) * t704); 0.4e1 * t703 * (t318 * t479 + t320 * t477 + t463 * t502); (t170 * t76 + t171 * t75 + t183 * t84 + t184 * t83) * m(7) + t693 * t480 + t694 * t478 + t656 * t321 + t655 * t319 + t748; (t108 * t75 + t109 * t76 + t159 * t28 + t183 * t40 + t184 * t41 + t52 * t99) * m(7) + t46 * t769 + t6 * t765 + t56 * t766 + t9 * t763 + t47 * t770 + t5 * t764 + t7 * t760 + (t1 * t759 + t2 * t702 + (t42 * t759 + t43 * t701) * qJD(1)) * t550; (t100 * t52 + t111 * t76 + t112 * t75 + t159 * t26 + t183 * t35 + t184 * t34) * m(7) + t53 * t766 + t8 * t763 + t43 * t768 + t1 * t761 + t42 * t767 + t2 * t762 + t45 * t770 + t3 * t764 + t44 * t769 + t4 * t765 + (t48 * t652 - t756 * t7 / 0.2e1) * t550; (t183 * t775 + t184 * t774 + t75 * t773 + t76 * t772 + (t159 * t709 - t52 * t756) * t704) * t771; (t159 * t463 + t183 * t318 + t184 * t320 + t477 * t75 + t479 * t76 + t502 * t52) * m(7); t319 * t43 + t480 * t1 + t321 * t42 + t478 * t2 + t464 * t48 + t503 * t7 + (t159 * t52 + t183 * t76 + t184 * t75) * t642;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;