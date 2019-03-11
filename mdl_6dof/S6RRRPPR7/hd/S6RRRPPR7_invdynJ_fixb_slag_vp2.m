% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:34
% EndTime: 2019-03-09 15:57:40
% DurationCPUTime: 43.17s
% Computational Cost: add. (12608->978), mult. (27097->1255), div. (0->0), fcn. (18686->12), ass. (0->433)
t364 = cos(qJ(3));
t543 = pkin(8) - qJ(5);
t304 = t543 * t364;
t360 = sin(qJ(3));
t361 = sin(qJ(2));
t470 = -pkin(7) * t360 - pkin(3);
t438 = -pkin(4) + t470;
t365 = cos(qJ(2));
t503 = t364 * t365;
t430 = pkin(2) * t361 - pkin(8) * t365;
t292 = t430 * qJD(1);
t512 = t292 * t364;
t690 = t512 - (-qJ(5) * t503 + t361 * t438) * qJD(1) + qJD(3) * t304 - qJD(5) * t360;
t257 = t360 * t292;
t492 = qJD(1) * t361;
t328 = qJ(4) * t492;
t483 = qJD(5) * t364;
t487 = qJD(3) * t360;
t506 = t361 * t364;
t507 = t360 * t365;
t689 = t257 + t328 + (-pkin(7) * t506 + qJ(5) * t507) * qJD(1) + t487 * t543 + t483;
t680 = Ifges(4,1) + Ifges(5,1);
t679 = Ifges(5,4) + Ifges(4,5);
t678 = Ifges(4,6) - Ifges(5,6);
t356 = cos(pkin(10));
t355 = sin(pkin(10));
t510 = t355 * t364;
t395 = -t356 * t360 + t510;
t385 = t395 * t365;
t215 = qJD(1) * t385;
t247 = t395 * qJD(3);
t688 = t215 - t247;
t511 = t355 * t360;
t394 = t356 * t364 + t511;
t386 = t365 * t394;
t216 = qJD(1) * t386;
t248 = t394 * qJD(3);
t687 = t216 - t248;
t686 = -mrSges(6,3) - mrSges(7,3);
t641 = t689 * t355 + t356 * t690;
t640 = t355 * t690 - t689 * t356;
t358 = -pkin(9) - qJ(5);
t602 = -m(7) * t358 - t686;
t681 = mrSges(5,2) + mrSges(4,3);
t685 = m(6) * qJ(5) + t602 - t681;
t359 = sin(qJ(6));
t363 = cos(qJ(6));
t469 = t360 * t492;
t482 = t364 * qJD(2);
t287 = t469 - t482;
t467 = t364 * t492;
t288 = qJD(2) * t360 + t467;
t439 = t287 * t355 + t356 * t288;
t440 = t356 * t287 - t288 * t355;
t684 = -t359 * t439 + t363 * t440;
t104 = t359 * t440 + t363 * t439;
t491 = qJD(1) * t365;
t318 = -qJD(3) + t491;
t556 = t318 / 0.2e1;
t312 = qJD(6) + t318;
t558 = t312 / 0.2e1;
t567 = t439 / 0.2e1;
t569 = t440 / 0.2e1;
t574 = t104 / 0.2e1;
t576 = t684 / 0.2e1;
t683 = Ifges(6,5) * t567 + Ifges(7,5) * t574 + Ifges(6,6) * t569 + Ifges(7,6) * t576 + Ifges(6,3) * t556 + Ifges(7,3) * t558;
t481 = qJD(1) * qJD(2);
t295 = qJDD(1) * t361 + t365 * t481;
t488 = qJD(3) * t287;
t174 = qJDD(2) * t360 + t295 * t364 - t488;
t573 = t174 / 0.2e1;
t175 = qJD(3) * t288 - t364 * qJDD(2) + t295 * t360;
t571 = t175 / 0.2e1;
t294 = t365 * qJDD(1) - t361 * t481;
t278 = qJDD(3) - t294;
t564 = t278 / 0.2e1;
t562 = t287 / 0.2e1;
t561 = -t288 / 0.2e1;
t557 = -t318 / 0.2e1;
t682 = pkin(9) * t440;
t677 = -Ifges(4,3) - Ifges(5,2);
t676 = pkin(5) * t492 + pkin(9) * t687 + t641;
t675 = pkin(9) * t688 + t640;
t622 = -t360 * t678 + t364 * t679;
t537 = Ifges(5,5) * t360;
t540 = Ifges(4,4) * t360;
t619 = t364 * t680 + t537 - t540;
t337 = pkin(7) * t491;
t485 = qJD(3) * t364;
t673 = -t360 * qJD(4) - t337 + (t364 * t491 - t485) * qJ(4);
t649 = -mrSges(4,2) + mrSges(5,3);
t657 = -t649 - m(7) * (pkin(5) * t355 + qJ(4));
t366 = cos(qJ(1));
t362 = sin(qJ(1));
t504 = t362 * t365;
t253 = t360 * t504 + t364 * t366;
t501 = t366 * t360;
t254 = t362 * t503 - t501;
t353 = pkin(10) + qJ(6);
t338 = sin(t353);
t339 = cos(t353);
t401 = t253 * t338 + t254 * t339;
t614 = t253 * t339 - t254 * t338;
t671 = mrSges(7,1) * t614 - t401 * mrSges(7,2);
t397 = t338 * t360 + t339 * t364;
t398 = t338 * t364 - t339 * t360;
t422 = t364 * mrSges(5,1) + t360 * mrSges(5,3);
t424 = mrSges(4,1) * t364 - mrSges(4,2) * t360;
t655 = t394 * mrSges(6,1) - t395 * mrSges(6,2);
t670 = t397 * mrSges(7,1) - t398 * mrSges(7,2) + t422 + t424 + t655;
t669 = -m(5) - m(6);
t572 = -t175 / 0.2e1;
t565 = -t278 / 0.2e1;
t652 = t294 / 0.2e1;
t668 = t295 / 0.2e1;
t667 = t679 * t564 + (-Ifges(4,4) + Ifges(5,5)) * t571 + t680 * t573;
t578 = pkin(3) + pkin(4);
t296 = -qJ(4) * t355 - t356 * t578;
t289 = -pkin(5) + t296;
t297 = t356 * qJ(4) - t355 * t578;
t186 = t289 * t363 - t297 * t359;
t396 = t355 * t359 - t356 * t363;
t345 = t361 * pkin(8);
t350 = t365 * pkin(2);
t471 = -pkin(1) - t350;
t393 = t471 - t345;
t263 = t393 * qJD(1);
t306 = qJD(2) * pkin(8) + t337;
t188 = t364 * t263 - t360 * t306;
t134 = qJ(5) * t288 + t188;
t189 = t360 * t263 + t364 * t306;
t135 = qJ(5) * t287 + t189;
t62 = -t134 * t355 + t356 * t135;
t51 = t62 + t682;
t63 = t356 * t134 + t355 * t135;
t650 = pkin(9) * t439;
t52 = t63 + t650;
t664 = -qJD(4) * t396 + qJD(6) * t186 - t359 * t51 - t363 * t52;
t187 = t289 * t359 + t297 * t363;
t282 = t355 * t363 + t356 * t359;
t663 = -qJD(4) * t282 - qJD(6) * t187 + t359 * t52 - t363 * t51;
t272 = Ifges(4,4) * t287;
t531 = t287 * Ifges(5,5);
t662 = t288 * t680 - t318 * t679 - t272 + t531;
t468 = t360 * t491;
t474 = t578 * t360;
t632 = -qJD(3) * t474 + t468 * t578 - t673;
t309 = t318 * qJ(4);
t154 = -t309 + t189;
t661 = -t154 * mrSges(5,2) - t189 * mrSges(4,3);
t152 = pkin(3) * t318 + qJD(4) - t188;
t660 = t152 * mrSges(5,2) - t188 * mrSges(4,3);
t347 = t364 * pkin(4);
t344 = t360 * qJ(4);
t607 = -t364 * pkin(3) - pkin(2) - t344;
t659 = t607 - t347;
t137 = -mrSges(6,2) * t318 + mrSges(6,3) * t440;
t138 = mrSges(6,1) * t318 - mrSges(6,3) * t439;
t210 = -mrSges(5,2) * t287 - mrSges(5,3) * t318;
t658 = t137 * t356 - t138 * t355 + t210;
t108 = t318 * t578 + qJD(4) - t134;
t114 = t135 - t309;
t53 = t356 * t108 - t114 * t355;
t39 = pkin(5) * t318 + t53 - t650;
t54 = t355 * t108 + t356 * t114;
t45 = t54 + t682;
t12 = -t359 * t45 + t363 * t39;
t13 = t359 * t39 + t363 * t45;
t542 = Ifges(3,4) * t361;
t643 = t365 * Ifges(3,2);
t414 = t542 + t643;
t656 = t53 * mrSges(6,1) + t12 * mrSges(7,1) - t54 * mrSges(6,2) - t13 * mrSges(7,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t414 / 0.2e1 + t677 * t557 + t679 * t561 + t678 * t562 + t683;
t326 = pkin(5) * t356 + pkin(4);
t597 = m(6) * pkin(4) + m(7) * t326 + mrSges(4,1) + mrSges(5,1);
t90 = -t174 * t355 + t175 * t356;
t91 = t174 * t356 + t175 * t355;
t21 = qJD(6) * t684 + t359 * t90 + t363 * t91;
t591 = t21 / 0.2e1;
t22 = -qJD(6) * t104 - t359 * t91 + t363 * t90;
t590 = t22 / 0.2e1;
t654 = -t91 * Ifges(6,4) / 0.2e1 - t90 * Ifges(6,2) / 0.2e1 + Ifges(6,6) * t564;
t580 = t90 / 0.2e1;
t579 = t91 / 0.2e1;
t653 = m(6) + m(7);
t262 = qJDD(6) - t278;
t566 = t262 / 0.2e1;
t489 = qJD(2) * t365;
t477 = pkin(7) * t489;
t648 = -mrSges(3,3) + mrSges(2,2);
t302 = t543 * t360;
t200 = t356 * t302 - t304 * t355;
t150 = pkin(9) * t395 + t200;
t201 = t355 * t302 + t356 * t304;
t151 = -pkin(9) * t394 + t201;
t67 = t150 * t363 - t151 * t359;
t647 = qJD(6) * t67 + t359 * t676 + t363 * t675;
t68 = t150 * t359 + t151 * t363;
t646 = -qJD(6) * t68 - t359 * t675 + t363 * t676;
t644 = Ifges(6,4) * t439;
t426 = mrSges(3,1) * t365 - mrSges(3,2) * t361;
t639 = -t426 - mrSges(2,1);
t634 = -pkin(5) * t688 + t632;
t631 = (-t468 + t487) * pkin(3) + t673;
t630 = t312 * t396;
t629 = t312 * t282;
t271 = Ifges(5,5) * t288;
t159 = -t318 * Ifges(5,6) + t287 * Ifges(5,3) + t271;
t335 = Ifges(3,4) * t491;
t628 = Ifges(3,1) * t492 + Ifges(3,5) * qJD(2) + t360 * t159 + t335;
t627 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t287 + mrSges(4,2) * t288 + mrSges(3,3) * t492;
t490 = qJD(2) * t361;
t626 = qJ(4) * t490 - qJD(4) * t365;
t625 = -t361 * t677 + t365 * t622;
t624 = t361 * t679 + t365 * t619;
t336 = pkin(7) * t492;
t305 = -qJD(2) * pkin(2) + t336;
t158 = t287 * pkin(3) - t288 * qJ(4) + t305;
t525 = t364 * mrSges(5,3);
t421 = t360 * mrSges(5,1) - t525;
t423 = mrSges(4,1) * t360 + mrSges(4,2) * t364;
t623 = -t158 * t421 - t305 * t423;
t621 = t360 * t679 + t364 * t678;
t536 = Ifges(5,5) * t364;
t539 = Ifges(4,4) * t364;
t620 = t360 * t680 - t536 + t539;
t616 = t174 * t679 - t175 * t678 - t278 * t677;
t276 = t294 * pkin(7);
t277 = t295 * pkin(7);
t615 = t276 * t365 + t277 * t361;
t613 = t681 * t361;
t521 = qJDD(1) * pkin(1);
t191 = -pkin(2) * t294 - pkin(8) * t295 - t521;
t242 = qJDD(2) * pkin(8) + t276;
t76 = t360 * t191 + t364 * t242 + t263 * t485 - t306 * t487;
t77 = t191 * t364 - t360 * t242 - t263 * t487 - t306 * t485;
t612 = -t360 * t77 + t364 * t76;
t57 = t278 * qJ(4) - t318 * qJD(4) + t76;
t384 = qJDD(4) - t77;
t60 = -pkin(3) * t278 + t384;
t611 = t360 * t60 + t364 * t57;
t610 = g(1) * t366 + g(2) * t362;
t608 = -m(5) - t653;
t605 = -Ifges(5,5) * t174 / 0.2e1 + Ifges(5,6) * t565 + Ifges(4,4) * t573 + Ifges(4,6) * t564 + (Ifges(5,3) + Ifges(4,2)) * t572;
t391 = pkin(5) * t511 + t326 * t364;
t603 = t685 * t365 + (-m(7) * (-t391 + t607) - m(6) * t659 - m(5) * t607 + m(4) * pkin(2) + t670) * t361;
t601 = pkin(8) * (-m(4) + t608);
t38 = -qJ(5) * t174 - qJD(5) * t288 - t278 * t578 + t384;
t41 = qJ(5) * t175 + qJD(5) * t287 + t57;
t10 = -t355 * t41 + t356 * t38;
t8 = -pkin(5) * t278 - pkin(9) * t91 + t10;
t11 = t355 * t38 + t356 * t41;
t9 = pkin(9) * t90 + t11;
t1 = qJD(6) * t12 + t359 * t8 + t363 * t9;
t2 = -qJD(6) * t13 - t359 * t9 + t363 * t8;
t600 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t595 = -t77 * mrSges(4,1) + t60 * mrSges(5,1) + t10 * mrSges(6,1) + t76 * mrSges(4,2) - t11 * mrSges(6,2) - t57 * mrSges(5,3);
t593 = Ifges(7,4) * t591 + Ifges(7,2) * t590 + Ifges(7,6) * t566;
t592 = Ifges(7,1) * t591 + Ifges(7,4) * t590 + Ifges(7,5) * t566;
t589 = Ifges(6,1) * t579 + Ifges(6,4) * t580 + Ifges(6,5) * t565;
t538 = Ifges(7,4) * t104;
t43 = Ifges(7,2) * t684 + Ifges(7,6) * t312 + t538;
t588 = -t43 / 0.2e1;
t587 = t43 / 0.2e1;
t94 = Ifges(7,4) * t684;
t44 = Ifges(7,1) * t104 + Ifges(7,5) * t312 + t94;
t586 = -t44 / 0.2e1;
t585 = t44 / 0.2e1;
t87 = Ifges(6,2) * t440 + t318 * Ifges(6,6) + t644;
t584 = -t87 / 0.2e1;
t583 = t87 / 0.2e1;
t88 = Ifges(6,1) * t439 + Ifges(6,4) * t440 + t318 * Ifges(6,5);
t582 = -t88 / 0.2e1;
t581 = t88 / 0.2e1;
t577 = -t684 / 0.2e1;
t575 = -t104 / 0.2e1;
t570 = -t440 / 0.2e1;
t568 = -t439 / 0.2e1;
t563 = -t287 / 0.2e1;
t560 = t288 / 0.2e1;
t559 = -t312 / 0.2e1;
t552 = mrSges(6,3) * t53;
t551 = mrSges(6,3) * t54;
t550 = mrSges(7,3) * t12;
t549 = mrSges(7,3) * t13;
t548 = pkin(7) * t361;
t293 = t430 * qJD(2);
t495 = t350 + t345;
t299 = -pkin(1) - t495;
t320 = pkin(7) * t503;
t499 = qJD(3) * t320 + t299 * t487;
t83 = (-qJ(5) * t489 - t293) * t364 + (qJ(5) * t487 + qJD(2) * t438 - t483) * t361 + t499;
t500 = t360 * t293 + t299 * t485;
t84 = (-pkin(7) * qJD(2) + qJ(5) * qJD(3)) * t506 + (qJD(5) * t361 + (-pkin(7) * qJD(3) + qJ(5) * qJD(2)) * t365) * t360 + t500 + t626;
t37 = t355 * t83 + t356 * t84;
t541 = Ifges(3,4) * t365;
t530 = t288 * Ifges(4,4);
t522 = qJ(5) * t361;
t518 = t253 * t355;
t508 = t360 * t361;
t505 = t361 * t366;
t502 = t365 * t366;
t319 = pkin(7) * t507;
t349 = t365 * pkin(3);
t172 = pkin(4) * t365 + t319 + t349 + (-t299 - t522) * t364;
t218 = t360 * t299 + t320;
t203 = -qJ(4) * t365 + t218;
t185 = qJ(5) * t508 + t203;
t93 = t355 * t172 + t356 * t185;
t194 = t288 * pkin(3) + t287 * qJ(4);
t466 = t365 * t482;
t497 = qJ(4) * t466 + qJD(4) * t506;
t494 = t366 * pkin(1) + t362 * pkin(7);
t486 = qJD(3) * t361;
t479 = Ifges(7,5) * t21 + Ifges(7,6) * t22 + Ifges(7,3) * t262;
t478 = pkin(7) * t490;
t476 = pkin(8) * t487;
t475 = pkin(8) * t485;
t243 = -qJDD(2) * pkin(2) + t277;
t465 = t360 * t486;
t162 = -t287 * Ifges(4,2) - t318 * Ifges(4,6) + t530;
t464 = -t360 * t162 / 0.2e1;
t46 = -t90 * mrSges(6,1) + t91 * mrSges(6,2);
t7 = -t22 * mrSges(7,1) + t21 * mrSges(7,2);
t456 = -t491 / 0.2e1;
t453 = t489 / 0.2e1;
t450 = -t486 / 0.2e1;
t449 = t485 / 0.2e1;
t36 = -t355 * t84 + t356 * t83;
t447 = t481 / 0.2e1;
t118 = -t278 * mrSges(5,1) + t174 * mrSges(5,2);
t231 = t355 * t506 - t356 * t508;
t232 = t394 * t361;
t446 = t231 * mrSges(6,1) + t232 * mrSges(6,2);
t92 = t356 * t172 - t185 * t355;
t442 = -t253 * t356 + t254 * t355;
t255 = -t362 * t364 + t365 * t501;
t256 = t360 * t362 + t364 * t502;
t441 = t255 * t355 + t256 * t356;
t217 = t299 * t364 - t319;
t436 = pkin(2) * t502 + pkin(8) * t505 + t494;
t432 = -pkin(7) - t474;
t141 = -pkin(4) * t288 - t194;
t431 = t470 * t361;
t427 = -t293 * t364 + t499;
t425 = mrSges(3,1) * t361 + mrSges(3,2) * t365;
t143 = t255 * t339 - t256 * t338;
t144 = t255 * t338 + t256 * t339;
t420 = mrSges(7,1) * t143 - mrSges(7,2) * t144;
t419 = (-mrSges(7,1) * t398 - mrSges(7,2) * t397) * t361;
t413 = -Ifges(4,2) * t360 + t539;
t412 = Ifges(4,2) * t364 + t540;
t409 = Ifges(3,5) * t365 - Ifges(3,6) * t361;
t406 = Ifges(5,3) * t360 + t536;
t405 = -Ifges(5,3) * t364 + t537;
t404 = t355 * t53 - t356 * t54;
t61 = pkin(5) * t365 - pkin(9) * t232 + t92;
t64 = -pkin(9) * t231 + t93;
t30 = -t359 * t64 + t363 * t61;
t31 = t359 * t61 + t363 * t64;
t403 = t256 * pkin(3) + t436;
t139 = -t231 * t363 - t232 * t359;
t402 = t231 * t359 - t232 * t363;
t400 = t254 * t356 + t518;
t399 = t255 * t356 - t256 * t355;
t176 = t359 * t395 - t363 * t394;
t177 = -t359 * t394 - t363 * t395;
t206 = -pkin(7) * t467 + t257;
t59 = t175 * pkin(3) - t174 * qJ(4) - t288 * qJD(4) + t243;
t390 = pkin(1) * t425;
t387 = t361 * (Ifges(3,1) * t365 - t542);
t383 = Ifges(6,5) * t91 + Ifges(6,6) * t90 - Ifges(6,3) * t278;
t382 = t479 + t600;
t315 = qJ(4) * t506;
t202 = t361 * t432 + t315;
t381 = -t465 + t466;
t380 = t360 * t489 + t361 * t485;
t116 = -pkin(4) * t287 + qJD(5) - t158;
t378 = -g(1) * t255 - g(2) * t253 - g(3) * t508;
t374 = Ifges(4,6) * t361 + t365 * t413;
t373 = Ifges(5,6) * t361 + t365 * t406;
t48 = -pkin(4) * t175 + qJDD(5) - t59;
t132 = (-t361 * t482 - t365 * t487) * pkin(7) + t500;
t109 = (-t364 * t578 - t344) * t486 + t432 * t489 + t497;
t351 = t366 * pkin(7);
t301 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t491;
t265 = t423 * t361;
t240 = t255 * pkin(3);
t238 = t253 * pkin(3);
t228 = -t315 + (pkin(3) * t360 + pkin(7)) * t361;
t209 = mrSges(5,1) * t318 + mrSges(5,2) * t288;
t208 = -mrSges(4,1) * t318 - mrSges(4,3) * t288;
t207 = mrSges(4,2) * t318 - mrSges(4,3) * t287;
t205 = pkin(7) * t469 + t512;
t204 = -t217 + t349;
t197 = pkin(5) * t394 - t659;
t195 = mrSges(5,1) * t287 - mrSges(5,3) * t288;
t193 = qJD(1) * t431 - t512;
t192 = t206 + t328;
t147 = qJD(2) * t386 + t395 * t486;
t146 = -qJD(2) * t385 + t248 * t361;
t136 = pkin(5) * t231 + t202;
t133 = t360 * t478 - t427;
t131 = pkin(3) * t380 + qJ(4) * t465 + t477 - t497;
t123 = -t215 * t359 + t216 * t363;
t122 = -t215 * t363 - t216 * t359;
t121 = qJD(2) * t431 + t427;
t120 = -mrSges(5,2) * t175 + mrSges(5,3) * t278;
t119 = -mrSges(4,2) * t278 - mrSges(4,3) * t175;
t117 = mrSges(4,1) * t278 - mrSges(4,3) * t174;
t113 = t132 + t626;
t105 = -mrSges(6,1) * t440 + mrSges(6,2) * t439;
t98 = mrSges(4,1) * t175 + mrSges(4,2) * t174;
t97 = mrSges(5,1) * t175 - mrSges(5,3) * t174;
t96 = -qJD(6) * t177 - t247 * t363 - t248 * t359;
t95 = qJD(6) * t176 - t247 * t359 + t248 * t363;
t89 = -pkin(5) * t439 + t141;
t79 = mrSges(7,1) * t312 - mrSges(7,3) * t104;
t78 = -mrSges(7,2) * t312 + mrSges(7,3) * t684;
t75 = -pkin(5) * t440 + t116;
t66 = -mrSges(6,1) * t278 - mrSges(6,3) * t91;
t65 = mrSges(6,2) * t278 + mrSges(6,3) * t90;
t58 = -pkin(5) * t146 + t109;
t50 = qJD(6) * t402 + t146 * t363 - t147 * t359;
t49 = qJD(6) * t139 + t146 * t359 + t147 * t363;
t47 = -mrSges(7,1) * t684 + mrSges(7,2) * t104;
t29 = pkin(9) * t146 + t37;
t26 = -pkin(5) * t490 - pkin(9) * t147 + t36;
t25 = -pkin(5) * t90 + t48;
t17 = -mrSges(7,2) * t262 + mrSges(7,3) * t22;
t16 = mrSges(7,1) * t262 - mrSges(7,3) * t21;
t4 = -qJD(6) * t31 + t26 * t363 - t29 * t359;
t3 = qJD(6) * t30 + t26 * t359 + t29 * t363;
t5 = [(-Ifges(7,5) * t402 + Ifges(7,6) * t139) * t566 + (-Ifges(7,4) * t402 + Ifges(7,2) * t139) * t590 + (-Ifges(7,1) * t402 + Ifges(7,4) * t139) * t591 + (t1 * t139 - t12 * t49 + t13 * t50 + t2 * t402) * mrSges(7,3) + (-t188 * t381 - t189 * t380 - t506 * t77) * mrSges(4,3) + m(4) * (t132 * t189 + t133 * t188 + t217 * t77 + t218 * t76 + t305 * t477) + (Ifges(7,5) * t49 + Ifges(7,6) * t50) * t558 + (Ifges(6,5) * t147 + Ifges(6,6) * t146) * t556 + (Ifges(6,5) * t232 - Ifges(6,6) * t231) * t565 + (-t57 * mrSges(5,2) - t76 * mrSges(4,3) - t605) * t508 + (Ifges(7,4) * t49 + Ifges(7,2) * t50) * t576 + (Ifges(6,4) * t147 + Ifges(6,2) * t146) * t569 + (Ifges(6,4) * t232 - Ifges(6,2) * t231) * t580 + (t152 * t381 - t154 * t380 + t506 * t60) * mrSges(5,2) + (Ifges(7,1) * t49 + Ifges(7,4) * t50) * t574 + (Ifges(6,1) * t147 + Ifges(6,4) * t146) * t567 + (Ifges(6,1) * t232 - Ifges(6,4) * t231) * t579 + t25 * (-mrSges(7,1) * t139 - mrSges(7,2) * t402) - t402 * t592 + t426 * t521 + (-qJDD(2) * mrSges(3,1) + t98) * t548 + (qJD(2) * t624 - t486 * t620) * t560 + (qJD(2) * t625 - t486 * t621) * t557 + t627 * t477 + t628 * t453 + (-t10 * t232 - t11 * t231 + t146 * t54 - t147 * t53) * mrSges(6,3) + m(7) * (t1 * t31 + t12 * t4 + t13 * t3 + t136 * t25 + t2 * t30 + t58 * t75) + m(6) * (t10 * t92 + t109 * t116 + t11 * t93 + t202 * t48 + t36 * t53 + t37 * t54) + m(5) * (t113 * t154 + t121 * t152 + t131 * t158 + t203 * t57 + t204 * t60 + t228 * t59) - t390 * t481 + (t295 * t548 + t615) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t615) + (-t679 * t573 + t541 * t447 + Ifges(6,3) * t565 + Ifges(7,3) * t566 - Ifges(5,6) * t571 - Ifges(4,6) * t572 + Ifges(6,5) * t579 + Ifges(6,6) * t580 + Ifges(7,6) * t590 + Ifges(7,5) * t591 + t677 * t564 + pkin(7) * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t294) + t595 + Ifges(3,6) * qJDD(2) + t600 + Ifges(3,4) * t668 + Ifges(3,2) * t652 + t383 / 0.2e1 + t479 / 0.2e1 - t616 / 0.2e1) * t365 - t301 * t478 + qJD(2) ^ 2 * t409 / 0.2e1 + (m(7) * pkin(5) * t518 + t400 * mrSges(6,1) + t401 * mrSges(7,1) - t442 * mrSges(6,2) + t614 * mrSges(7,2) + t608 * (-t254 * pkin(3) - qJ(4) * t253 + t351) + t648 * t366 + (-m(3) - m(4)) * t351 + t597 * t254 + t649 * t253 + (m(3) * pkin(1) - t653 * t471 + (-m(4) - m(5)) * t393 + (m(6) * t543 - m(7) * (-pkin(8) - t358) + t686) * t361 + t613 - t639) * t362) * g(1) + t231 * t654 + (m(4) * t243 * pkin(7) + Ifges(3,1) * t295 + Ifges(3,4) * t652 + Ifges(3,5) * qJDD(2) + t159 * t449 + t406 * t571 + t413 * t572 + t59 * t421 - t447 * t643 + t564 * t622 + t573 * t619) * t361 + t662 * (t360 * t450 + t364 * t453) + t48 * t446 + (-m(3) * t494 - m(4) * t436 - m(7) * t403 - t441 * mrSges(6,1) - t144 * mrSges(7,1) - t399 * mrSges(6,2) - t143 * mrSges(7,2) + t669 * (qJ(4) * t255 + t403) + t639 * t366 + t648 * t362 - t597 * t256 + t657 * t255 + t685 * t505) * g(2) + (t188 * mrSges(4,1) - t152 * mrSges(5,1) - t189 * mrSges(4,2) + t154 * mrSges(5,3) - t656 - t683) * t490 + Ifges(2,3) * qJDD(1) - pkin(1) * (-mrSges(3,1) * t294 + mrSges(3,2) * t295) + t243 * t265 + t228 * t97 + t217 * t117 + t218 * t119 + t121 * t209 + t113 * t210 + t305 * (mrSges(4,1) * t380 + mrSges(4,2) * t381) + t158 * (mrSges(5,1) * t380 - mrSges(5,3) * t381) + t132 * t207 + t133 * t208 + t202 * t46 + t203 * t120 + t204 * t118 + t131 * t195 + t116 * (-mrSges(6,1) * t146 + mrSges(6,2) * t147) + t37 * t137 + t36 * t138 + t136 * t7 + t109 * t105 + t92 * t66 + t93 * t65 + t75 * (-mrSges(7,1) * t50 + mrSges(7,2) * t49) + t3 * t78 + t4 * t79 + t58 * t47 + t30 * t16 + t31 * t17 + t364 * t162 * t450 + t387 * t447 + t464 * t489 + t506 * t667 + t541 * t668 + (qJD(2) * t373 - t405 * t486) * t562 + (qJD(2) * t374 - t412 * t486) * t563 + t147 * t581 + t146 * t583 + t49 * t585 + t50 * t587 + t232 * t589 + t139 * t593 + t414 * t652; (-mrSges(6,1) * t688 - mrSges(6,2) * t687) * t116 + (Ifges(6,4) * t216 - Ifges(6,2) * t215) * t570 + (t10 * t395 - t11 * t394 + t215 * t54 + t216 * t53) * mrSges(6,3) + (Ifges(6,5) * t216 - Ifges(6,6) * t215) * t557 + (Ifges(6,1) * t216 - Ifges(6,4) * t215) * t568 + t631 * t195 + t632 * t105 + t634 * t47 + (-pkin(2) * t243 - t188 * t205 - t189 * t206 - t305 * t337) * m(4) + (((t152 * t364 - t154 * t360) * qJD(3) + t611) * m(5) + (t120 + t119) * t364 + (-t117 + t118) * t360 + ((-t188 * t364 - t189 * t360) * qJD(3) + t612) * m(4)) * pkin(8) + t605 * t364 + (-t188 * (mrSges(4,1) * t361 - mrSges(4,3) * t503) - t152 * (-mrSges(5,1) * t361 + mrSges(5,2) * t503) - t189 * (-mrSges(4,2) * t361 - mrSges(4,3) * t507) - t154 * (-mrSges(5,2) * t507 + mrSges(5,3) * t361) + (t374 / 0.2e1 - t373 / 0.2e1) * t287 + t624 * t561 + t625 * t556 + (-t387 / 0.2e1 + t390) * qJD(1)) * qJD(1) + (-m(4) * t495 + m(6) * t522 - t426 + t608 * (pkin(3) * t503 + qJ(4) * t507 + t495) + t602 * t361 + (-m(6) * t347 - m(7) * t391 - t670) * t365 - t613) * g(3) - (Ifges(6,4) * t567 + Ifges(6,2) * t569 + Ifges(6,6) * t556 + t551 + t583) * t247 + t394 * t654 + (-Ifges(6,5) * t395 - Ifges(6,6) * t394) * t565 + (-Ifges(6,1) * t395 - Ifges(6,4) * t394) * t579 + (-Ifges(6,4) * t395 - Ifges(6,2) * t394) * t580 - t395 * t589 + (-t152 * t193 - t154 * t192 + t631 * t158 + t59 * t607) * m(5) + t607 * t97 + t646 * t79 + (t1 * t68 + t12 * t646 + t13 * t647 + t197 * t25 + t2 * t67 + t634 * t75) * m(7) + t647 * t78 + (-Ifges(6,5) * t568 - Ifges(7,5) * t575 - Ifges(3,2) * t456 - Ifges(6,6) * t570 - Ifges(7,6) * t577 - Ifges(6,3) * t557 - Ifges(7,3) * t559 + t656) * t492 + (t335 + t628) * t456 + t48 * t655 + (t362 * t603 + t504 * t601) * g(2) + (t366 * t603 + t502 * t601) * g(1) + (Ifges(7,4) * t574 + Ifges(7,2) * t576 + Ifges(7,6) * t558 + t549 + t587) * t96 + t620 * t573 + t621 * t564 + t623 * t491 - t627 * t337 - t59 * t422 - t243 * t424 - t409 * t481 / 0.2e1 + (Ifges(7,5) * t123 + Ifges(7,6) * t122) * t559 + t610 * t425 + t611 * mrSges(5,2) + t612 * mrSges(4,3) + (t1 * t176 + t12 * t123 - t122 * t13 - t177 * t2) * mrSges(7,3) + t162 * t468 / 0.2e1 + (-t476 - t206) * t207 + (-t476 - t192) * t210 + t662 * (t364 * t456 + t449) + (t10 * t200 + t11 * t201 + t116 * t632 - t48 * t659 + t53 * t641 + t54 * t640) * m(6) - t659 * t46 + (Ifges(6,1) * t567 + Ifges(6,4) * t569 + Ifges(6,5) * t556 - t552 + t581) * t248 + (Ifges(7,4) * t123 + Ifges(7,2) * t122) * t577 + t640 * t137 + t641 * t138 + (Ifges(7,1) * t574 + Ifges(7,4) * t576 + Ifges(7,5) * t558 - t550 + t585) * t95 + (-t413 / 0.2e1 + t406 / 0.2e1) * t488 + (Ifges(7,1) * t123 + Ifges(7,4) * t122) * t575 + (t622 * t557 + t619 * t560 + t464 - t623) * qJD(3) + (-t193 + t475) * t209 + t660 * t485 + (t159 / 0.2e1 + t661) * t487 + Ifges(3,6) * t294 + Ifges(3,5) * t295 - t276 * mrSges(3,2) - t277 * mrSges(3,1) + t197 * t7 + t200 * t66 + t201 * t65 + t25 * (-mrSges(7,1) * t176 + mrSges(7,2) * t177) - pkin(2) * t98 + Ifges(3,3) * qJDD(2) + t67 * t16 + t68 * t17 + ((-t123 + t95) * mrSges(7,2) + (t122 - t96) * mrSges(7,1)) * t75 + t301 * t336 + t360 * t667 + (Ifges(7,5) * t177 + Ifges(7,6) * t176) * t566 + t405 * t571 + t412 * t572 + t216 * t582 - t215 * t584 + t123 * t586 + t122 * t588 + (Ifges(7,4) * t177 + Ifges(7,2) * t176) * t590 + (Ifges(7,1) * t177 + Ifges(7,4) * t176) * t591 + t177 * t592 + t176 * t593 + (-t475 - t205) * t208; (-t210 - t207) * t188 + (m(5) * t154 - m(6) * t404 + t658) * qJD(4) + t663 * t79 + (t1 * t187 + t663 * t12 + t664 * t13 + t186 * t2 - t75 * t89) * m(7) + t664 * t78 + (-t305 * mrSges(4,1) - t158 * mrSges(5,1) + Ifges(5,3) * t563 - t556 * t678 - t661) * t288 + t616 - t531 * t563 + (mrSges(6,1) * t116 + Ifges(6,2) * t570 + Ifges(6,6) * t557 - t551 + t584) * t439 + (-t287 * t680 + t159 + t271 - t530) * t561 + (m(7) * t240 + t399 * mrSges(6,1) - t441 * mrSges(6,2) + t420 + t669 * (qJ(4) * t256 - t240) + t657 * t256 + t597 * t255) * g(1) + (-m(5) * t315 - (t525 + (-m(5) * pkin(3) - mrSges(5,1)) * t360) * t361 - m(6) * (-t361 * t474 + t315) - t446 + t419 - (t315 + (pkin(5) * t510 + (-pkin(3) - t326) * t360) * t361) * m(7) + t265) * g(3) + (t10 * t296 + t11 * t297 - t116 * t141 - t53 * t62 - t54 * t63) * m(6) - t595 + t644 * t568 + (-t209 + t208) * t189 + (-Ifges(4,2) * t288 - t272 + t662) * t562 + (t305 * mrSges(4,2) - t158 * mrSges(5,3) - t556 * t679 + t660) * t287 + (-pkin(3) * t60 + qJ(4) * t57 - t152 * t189 - t154 * t188 - t158 * t194) * m(5) + (mrSges(7,1) * t75 + Ifges(7,4) * t575 + Ifges(7,2) * t577 + Ifges(7,6) * t559 - t549 + t588) * t104 - (-mrSges(7,2) * t75 + Ifges(7,1) * t575 + Ifges(7,4) * t577 + Ifges(7,5) * t559 + t550 + t586) * t684 - (-mrSges(6,2) * t116 + Ifges(6,1) * t568 + Ifges(6,4) * t570 + Ifges(6,5) * t557 + t552 + t582) * t440 + t296 * t66 + t297 * t65 - t194 * t195 + t186 * t16 + t187 * t17 - t63 * t137 - t62 * t138 - t141 * t105 + qJ(4) * t120 - pkin(3) * t118 - t89 * t47 + t162 * t560 + (m(7) * t238 - t442 * mrSges(6,1) - t400 * mrSges(6,2) + t669 * (qJ(4) * t254 - t238) + t657 * t254 + t597 * t253 + t671) * g(2) - t382 - t383; -t396 * t16 + t282 * t17 + t355 * t65 + t356 * t66 - t629 * t79 - t630 * t78 + t658 * t318 + (t195 - t105 - t47) * t288 + t118 + (t1 * t282 - t12 * t629 - t13 * t630 - t2 * t396 - t288 * t75 + t378) * m(7) + (t10 * t356 + t11 * t355 - t116 * t288 - t318 * t404 + t378) * m(6) + (t154 * t318 + t158 * t288 + t378 + t60) * m(5); -t684 * t78 - t440 * t137 + t439 * t138 + t104 * t79 + t46 + t7 + (-g(3) * t365 + t361 * t610) * t653 + (t104 * t12 - t13 * t684 + t25) * m(7) + (t439 * t53 - t440 * t54 + t48) * m(6); -t75 * (mrSges(7,1) * t104 + mrSges(7,2) * t684) + (Ifges(7,1) * t684 - t538) * t575 + t43 * t574 + (Ifges(7,5) * t684 - Ifges(7,6) * t104) * t559 - t12 * t78 + t13 * t79 - g(1) * t420 - g(2) * t671 - g(3) * t419 + (t104 * t13 + t12 * t684) * mrSges(7,3) + t382 + (-Ifges(7,2) * t104 + t44 + t94) * t577;];
tau  = t5;
