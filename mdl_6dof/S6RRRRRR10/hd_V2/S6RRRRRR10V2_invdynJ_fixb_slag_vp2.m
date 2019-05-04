% Calculate vector of inverse dynamics joint torques for
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRR10V2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:39
% EndTime: 2019-04-11 14:43:17
% DurationCPUTime: 36.21s
% Computational Cost: add. (20850->982), mult. (44080->1363), div. (0->0), fcn. (34308->14), ass. (0->454)
t375 = sin(qJ(4));
t526 = qJD(4) * t375;
t376 = sin(qJ(3));
t377 = sin(qJ(2));
t382 = cos(qJ(3));
t383 = cos(qJ(2));
t316 = t376 * t377 - t382 * t383;
t306 = t316 * qJD(1);
t559 = t306 * t375;
t734 = -t526 - t559;
t373 = sin(qJ(6));
t379 = cos(qJ(6));
t449 = -mrSges(7,1) * t379 + mrSges(7,2) * t373;
t418 = mrSges(6,1) - t449;
t317 = t376 * t383 + t377 * t382;
t307 = t317 * qJD(1);
t235 = pkin(3) * t307 + pkin(5) * t306;
t381 = cos(qJ(4));
t583 = pkin(2) * qJD(2);
t511 = t382 * t583;
t206 = t235 * t375 + t381 * t511;
t374 = sin(qJ(5));
t380 = cos(qJ(5));
t520 = qJD(5) * t381;
t525 = qJD(4) * t380;
t399 = -t374 * t520 - t375 * t525;
t513 = t376 * t583;
t521 = qJD(5) * t380;
t697 = -pkin(3) * t521 + pkin(5) * t399 - t380 * t206 - t374 * t513;
t752 = t734 * pkin(6);
t536 = t380 * t381;
t215 = -t306 * t536 + t307 * t374;
t522 = qJD(5) * t375;
t484 = t374 * t522;
t751 = (t215 + t484) * pkin(6);
t518 = qJD(1) * qJD(2);
t320 = qJDD(1) * t383 - t377 * t518;
t321 = qJDD(1) * t377 + t383 * t518;
t404 = t316 * qJD(3);
t201 = -qJD(1) * t404 + t320 * t376 + t321 * t382;
t371 = qJD(2) + qJD(3);
t254 = t307 * t381 + t371 * t375;
t370 = qJDD(2) + qJDD(3);
t143 = -qJD(4) * t254 - t201 * t375 + t370 * t381;
t140 = qJDD(5) - t143;
t253 = -t307 * t375 + t371 * t381;
t142 = qJD(4) * t253 + t201 * t381 + t370 * t375;
t299 = qJD(4) + t306;
t193 = -t254 * t374 + t299 * t380;
t405 = t317 * qJD(3);
t202 = -qJD(1) * t405 + t320 * t382 - t321 * t376;
t198 = qJDD(4) - t202;
t77 = qJD(5) * t193 + t142 * t380 + t198 * t374;
t194 = t254 * t380 + t299 * t374;
t78 = -qJD(5) * t194 - t142 * t374 + t198 * t380;
t20 = Ifges(6,4) * t77 + Ifges(6,2) * t78 + Ifges(6,6) * t140;
t326 = pkin(5) * t371 + t513;
t621 = pkin(2) * t383;
t362 = pkin(1) + t621;
t331 = t362 * qJD(1);
t391 = t306 * pkin(3) - t307 * pkin(5) - t331;
t174 = t381 * t326 + t375 * t391;
t517 = qJDD(2) * t382;
t528 = qJD(3) * t376;
t312 = (qJD(2) * t528 - t517) * pkin(2);
t608 = t370 * pkin(3);
t390 = t608 - t312;
t607 = t371 * pkin(3);
t419 = t511 + t607;
t401 = t380 * t419;
t523 = qJD(5) * t374;
t582 = pkin(2) * qJD(3);
t509 = t382 * t582;
t624 = pkin(2) * t376;
t313 = qJD(2) * t509 + qJDD(2) * t624;
t288 = pkin(5) * t370 + t313;
t574 = qJDD(1) * pkin(1);
t297 = -pkin(2) * t320 - t574;
t393 = -pkin(3) * t202 - pkin(5) * t201 + t297;
t173 = t375 * t326 - t381 * t391;
t527 = qJD(4) * t173;
t69 = t381 * t288 + t375 * t393 - t527;
t34 = -qJD(5) * t401 - t174 * t523 - t374 * t390 + t380 * t69;
t74 = qJDD(6) - t78;
t655 = t74 / 0.2e1;
t250 = qJD(5) - t253;
t150 = t194 * t379 + t250 * t373;
t30 = -qJD(6) * t150 + t140 * t379 - t373 * t77;
t662 = t30 / 0.2e1;
t149 = -t194 * t373 + t250 * t379;
t29 = qJD(6) * t149 + t140 * t373 + t379 * t77;
t663 = t29 / 0.2e1;
t24 = pkin(6) * t140 + t34;
t70 = qJD(4) * t174 + t375 * t288 - t381 * t393;
t25 = -t77 * pkin(6) + t70;
t117 = -pkin(6) * t194 + t173;
t156 = t380 * t174 - t374 * t419;
t121 = t250 * pkin(6) + t156;
t59 = t117 * t373 + t121 * t379;
t4 = -qJD(6) * t59 - t24 * t373 + t25 * t379;
t720 = t4 * mrSges(7,1);
t58 = t117 * t379 - t121 * t373;
t3 = qJD(6) * t58 + t24 * t379 + t25 * t373;
t721 = t3 * mrSges(7,2);
t676 = t720 - t721;
t8 = Ifges(7,5) * t29 + Ifges(7,6) * t30 + Ifges(7,3) * t74;
t750 = -mrSges(6,3) * t34 + Ifges(7,5) * t663 + Ifges(7,6) * t662 + Ifges(7,3) * t655 - t20 / 0.2e1 + t8 / 0.2e1 + t676;
t372 = qJ(2) + qJ(3);
t368 = sin(t372);
t369 = cos(t372);
t285 = t368 * t536 - t369 * t374;
t544 = t374 * t381;
t416 = t368 * t544 + t369 * t380;
t604 = mrSges(5,1) * t381;
t667 = m(7) * pkin(6);
t507 = mrSges(7,3) + t667;
t737 = mrSges(6,2) - t507;
t749 = t285 * t418 + t368 * t604 - t416 * t737;
t713 = mrSges(6,3) - mrSges(5,2);
t748 = mrSges(7,1) * t373;
t448 = mrSges(7,2) * t379 + t748;
t708 = -mrSges(6,3) - t448;
t681 = mrSges(5,2) + t708;
t601 = mrSges(5,3) * t254;
t701 = -mrSges(5,1) * t299 - mrSges(6,1) * t193 + mrSges(6,2) * t194 + t601;
t205 = -t381 * t235 + t375 * t511;
t614 = pkin(6) * t380;
t493 = pkin(5) - t614;
t524 = qJD(4) * t381;
t747 = t493 * t524 - t205 + t751;
t531 = qJD(1) * t377;
t514 = pkin(2) * t531;
t219 = t235 + t514;
t464 = t375 * t509;
t360 = pkin(5) + t624;
t470 = t360 - t614;
t746 = t219 * t381 + t470 * t524 + t464 + t751;
t745 = t752 - t697;
t622 = pkin(2) * t382;
t361 = -pkin(3) - t622;
t463 = t381 * t509;
t510 = pkin(2) * t528;
t184 = t360 * t399 + t361 * t521 + t374 * t510 + t380 * t463;
t541 = t375 * t380;
t744 = t219 * t541 - t184 + t752;
t352 = pkin(5) * t536;
t488 = t374 * t526;
t743 = -pkin(3) * t523 - pkin(5) * t488 + qJD(5) * t352 - t206 * t374 + t380 * t513;
t398 = t374 * t524 + t375 * t521;
t19 = Ifges(6,5) * t77 + Ifges(6,6) * t78 + Ifges(6,3) * t140;
t35 = qJD(5) * t156 + t374 * t69 + t380 * t390;
t650 = t140 / 0.2e1;
t653 = t78 / 0.2e1;
t654 = t77 / 0.2e1;
t715 = t34 * mrSges(6,2);
t742 = -t715 + Ifges(6,5) * t654 + Ifges(6,6) * t653 + Ifges(6,3) * t650 + t19 / 0.2e1 - Ifges(5,4) * t142 / 0.2e1 - Ifges(5,2) * t143 / 0.2e1 - Ifges(5,6) * t198 / 0.2e1 - t35 * mrSges(6,1);
t739 = t156 * mrSges(6,2);
t581 = t156 * mrSges(6,3);
t738 = t173 * mrSges(6,1);
t152 = mrSges(6,1) * t250 - mrSges(6,3) * t194;
t87 = -mrSges(7,1) * t149 + mrSges(7,2) * t150;
t707 = -t152 + t87;
t542 = t375 * t379;
t166 = -t215 * t373 - t306 * t542;
t461 = -qJD(6) + t525;
t462 = -qJD(6) * t380 + qJD(4);
t209 = t462 * t542 + (-t381 * t461 + t484) * t373;
t736 = t166 - t209;
t547 = t373 * t375;
t167 = t215 * t379 - t306 * t547;
t537 = t379 * t381;
t208 = t461 * t537 + (t373 * t462 - t379 * t523) * t375;
t735 = t167 - t208;
t586 = Ifges(6,6) * t250;
t594 = Ifges(6,4) * t194;
t112 = Ifges(6,2) * t193 + t586 + t594;
t191 = qJD(6) - t193;
t584 = Ifges(7,3) * t191;
t585 = Ifges(7,6) * t149;
t587 = Ifges(7,5) * t150;
t60 = t584 + t585 + t587;
t733 = t112 / 0.2e1 - t60 / 0.2e1;
t732 = t60 / 0.2e1 - t581;
t101 = -mrSges(5,2) * t198 + mrSges(5,3) * t143;
t602 = mrSges(5,3) * t253;
t203 = -mrSges(5,2) * t299 + t602;
t567 = t173 * t381;
t424 = -t174 * t375 + t567;
t575 = t70 * t375;
t576 = t69 * t381;
t731 = m(5) * (qJD(4) * t424 + t575 + t576) + t101 * t381 - t203 * t526;
t190 = Ifges(6,4) * t193;
t588 = Ifges(6,5) * t250;
t113 = Ifges(6,1) * t194 + t190 + t588;
t651 = -t113 / 0.2e1;
t730 = -t173 * mrSges(6,2) + t651;
t384 = cos(qJ(1));
t532 = t384 * t375;
t378 = sin(qJ(1));
t539 = t378 * t381;
t294 = t369 * t539 - t532;
t554 = t368 * t378;
t239 = t294 * t380 + t374 * t554;
t533 = t381 * t384;
t540 = t378 * t375;
t293 = t369 * t540 + t533;
t729 = t239 * t373 - t293 * t379;
t728 = -t239 * t379 - t293 * t373;
t675 = t58 * mrSges(7,1) - t59 * mrSges(7,2);
t727 = t675 - t581;
t664 = Ifges(6,1) * t654 + Ifges(6,4) * t653 + Ifges(6,5) * t650;
t642 = t191 / 0.2e1;
t644 = t150 / 0.2e1;
t646 = t149 / 0.2e1;
t726 = Ifges(7,5) * t644 + Ifges(7,6) * t646 + Ifges(7,3) * t642;
t643 = -t191 / 0.2e1;
t645 = -t150 / 0.2e1;
t647 = -t149 / 0.2e1;
t725 = Ifges(7,5) * t645 + Ifges(7,6) * t647 + Ifges(7,3) * t643 + t733;
t722 = -m(5) - m(6);
t597 = Ifges(5,4) * t254;
t160 = Ifges(5,2) * t253 + Ifges(5,6) * t299 + t597;
t719 = -t160 / 0.2e1;
t249 = Ifges(5,4) * t253;
t161 = Ifges(5,1) * t254 + Ifges(5,5) * t299 + t249;
t718 = t161 / 0.2e1;
t717 = t320 / 0.2e1;
t155 = t174 * t374 + t401;
t716 = m(7) * t155;
t714 = t69 * mrSges(5,2);
t712 = -mrSges(6,1) * t140 - mrSges(7,1) * t30 + mrSges(7,2) * t29 + mrSges(6,3) * t77;
t278 = t360 * t536 + t374 * t361;
t613 = pkin(6) * t381;
t261 = t278 - t613;
t300 = t470 * t375;
t200 = t261 * t379 + t300 * t373;
t711 = -qJD(6) * t200 + t373 * t744 + t379 * t746;
t199 = -t261 * t373 + t300 * t379;
t710 = qJD(6) * t199 + t373 * t746 - t379 * t744;
t709 = t383 * Ifges(3,2);
t706 = -mrSges(5,1) * t198 - mrSges(6,1) * t78 + mrSges(6,2) * t77 + mrSges(5,3) * t142;
t451 = -mrSges(6,1) * t380 + mrSges(6,2) * t374;
t705 = t451 - mrSges(5,1);
t319 = -pkin(3) * t374 + t352;
t302 = t319 - t613;
t322 = t493 * t375;
t232 = t302 * t379 + t322 * t373;
t704 = -qJD(6) * t232 + t373 * t745 + t379 * t747;
t231 = -t302 * t373 + t322 * t379;
t703 = qJD(6) * t231 + t373 * t747 - t379 * t745;
t452 = mrSges(5,1) * t375 + mrSges(5,2) * t381;
t702 = t419 * t452;
t562 = t253 * t380;
t170 = t254 * t379 - t373 * t562;
t519 = qJD(6) * t374;
t692 = -t373 * t521 - t379 * t519;
t700 = t170 - t692;
t538 = t379 * t380;
t171 = t253 * t538 + t254 * t373;
t693 = t373 * t519 - t379 * t521;
t699 = t171 + t693;
t698 = mrSges(4,1) * t371 + mrSges(5,1) * t253 - mrSges(5,2) * t254 - mrSges(4,3) * t307;
t467 = t369 * mrSges(4,1) - mrSges(4,2) * t368;
t315 = t369 * pkin(3) + t368 * pkin(5);
t318 = pkin(3) * t380 + pkin(5) * t544;
t695 = t155 * t743 + t318 * t35;
t236 = pkin(3) * t316 - pkin(5) * t317 - t362;
t694 = t398 * t236;
t400 = t380 * t524 - t484;
t690 = t173 * t524 + t575;
t579 = t35 * t374;
t687 = t34 * t380 + t579;
t686 = g(1) * t384 + g(2) * t378;
t685 = Ifges(3,6) * qJDD(2);
t683 = -m(7) + t722;
t682 = -mrSges(3,3) - mrSges(4,3) + mrSges(2,2);
t457 = t507 * t374;
t680 = -t380 * t449 + t457 - t705;
t330 = -mrSges(3,1) * t383 + mrSges(3,2) * t377;
t625 = m(4) * t362;
t679 = m(3) * pkin(1) + mrSges(2,1) - t330 + t467 + t625;
t429 = t373 * t59 + t379 * t58;
t432 = Ifges(7,5) * t379 - Ifges(7,6) * t373;
t589 = Ifges(7,4) * t379;
t438 = -Ifges(7,2) * t373 + t589;
t590 = Ifges(7,4) * t373;
t444 = Ifges(7,1) * t379 - t590;
t144 = Ifges(7,4) * t149;
t62 = Ifges(7,1) * t150 + Ifges(7,5) * t191 + t144;
t577 = t379 * t62;
t591 = Ifges(7,4) * t150;
t61 = Ifges(7,2) * t149 + Ifges(7,6) * t191 + t591;
t578 = t373 * t61;
t678 = -t429 * mrSges(7,3) + t438 * t646 + t444 * t644 + t432 * t642 - t578 / 0.2e1 + t577 / 0.2e1;
t553 = t368 * t380;
t286 = t369 * t544 - t553;
t677 = -t368 * mrSges(5,3) - t467 + (mrSges(6,2) - mrSges(7,3)) * t286 - t418 * (t368 * t374 + t369 * t536) + (-mrSges(7,1) * t547 - mrSges(7,2) * t542 - t375 * t713 - t604) * t369;
t549 = t369 * t384;
t674 = -t368 * t532 * t681 - mrSges(5,3) * t549 + t384 * t749;
t502 = t368 * t542;
t550 = t369 * t378;
t673 = -mrSges(5,3) * t550 + (t748 + t713) * t368 * t540 + (mrSges(7,2) * t502 + t749) * t378;
t672 = -t675 - t738;
t151 = -mrSges(6,2) * t250 + mrSges(6,3) * t193;
t425 = t155 * t374 + t156 * t380;
t546 = t374 * t375;
t671 = m(5) * t424 - m(6) * (t375 * t425 - t567) - t151 * t541 + t152 * t546 - t375 * t203 + t701 * t381;
t9 = Ifges(7,4) * t29 + Ifges(7,2) * t30 + Ifges(7,6) * t74;
t668 = t9 / 0.2e1;
t666 = Ifges(7,1) * t663 + Ifges(7,4) * t662 + Ifges(7,5) * t655;
t659 = -t61 / 0.2e1;
t658 = t61 / 0.2e1;
t657 = -t62 / 0.2e1;
t656 = t62 / 0.2e1;
t649 = t142 / 0.2e1;
t648 = t143 / 0.2e1;
t641 = -t193 / 0.2e1;
t640 = t193 / 0.2e1;
t639 = -t194 / 0.2e1;
t638 = t194 / 0.2e1;
t637 = t198 / 0.2e1;
t636 = -t250 / 0.2e1;
t635 = t250 / 0.2e1;
t634 = -t253 / 0.2e1;
t633 = -t254 / 0.2e1;
t632 = t254 / 0.2e1;
t631 = -t299 / 0.2e1;
t628 = t307 / 0.2e1;
t627 = t375 / 0.2e1;
t626 = t383 / 0.2e1;
t623 = pkin(2) * t377;
t617 = pkin(6) * t193;
t600 = Ifges(3,4) * t377;
t599 = Ifges(3,4) * t383;
t598 = Ifges(4,4) * t307;
t596 = Ifges(5,4) * t375;
t595 = Ifges(5,4) * t381;
t593 = Ifges(6,4) * t374;
t592 = Ifges(6,4) * t380;
t570 = t173 * t205;
t569 = t173 * t374;
t568 = t173 * t380;
t566 = t236 * t381;
t245 = -qJD(2) * t316 - t404;
t565 = t245 * t375;
t564 = t245 * t381;
t563 = t253 * t374;
t558 = t306 * t381;
t557 = t317 * t375;
t552 = t368 * t384;
t548 = t373 * t374;
t545 = t374 * t379;
t530 = qJD(1) * t383;
t529 = qJD(2) * t377;
t516 = m(4) * t623;
t512 = pkin(2) * t529;
t246 = qJD(2) * t317 + t405;
t164 = pkin(3) * t246 - pkin(5) * t245 + t512;
t506 = t164 * t546;
t505 = t219 * t546;
t498 = Ifges(5,5) * t142 + Ifges(5,6) * t143 + Ifges(5,3) * t198;
t490 = t236 * t526;
t111 = Ifges(6,5) * t194 + Ifges(6,6) * t193 + Ifges(6,3) * t250;
t479 = t111 * t627;
t477 = -t526 / 0.2e1;
t476 = -t524 / 0.2e1;
t475 = t524 / 0.2e1;
t474 = -t523 / 0.2e1;
t473 = -t521 / 0.2e1;
t472 = t521 / 0.2e1;
t469 = t331 * t516;
t468 = t518 / 0.2e1;
t466 = mrSges(4,3) * t513;
t465 = mrSges(4,3) * t511;
t458 = pkin(6) * t286 + t315;
t456 = -pkin(3) * t368 - t623;
t413 = t317 * t526 - t564;
t114 = (-t317 * t520 + t246) * t374 + (qJD(5) * t316 - t413) * t380;
t455 = qJD(6) * t557 + t114;
t454 = mrSges(3,1) * t377 + mrSges(3,2) * t383;
t453 = mrSges(4,1) * t368 + mrSges(4,2) * t369;
t450 = mrSges(6,1) * t374 + mrSges(6,2) * t380;
t447 = Ifges(5,1) * t381 - t596;
t446 = Ifges(6,1) * t380 - t593;
t445 = Ifges(6,1) * t374 + t592;
t443 = Ifges(7,1) * t373 + t589;
t442 = t600 + t709;
t441 = -Ifges(5,2) * t375 + t595;
t440 = -Ifges(6,2) * t374 + t592;
t439 = Ifges(6,2) * t380 + t593;
t437 = Ifges(7,2) * t379 + t590;
t436 = Ifges(3,5) * t383 - Ifges(3,6) * t377;
t435 = Ifges(5,5) * t381 - Ifges(5,6) * t375;
t434 = Ifges(6,5) * t380 - Ifges(6,6) * t374;
t433 = Ifges(6,5) * t374 + Ifges(6,6) * t380;
t431 = Ifges(7,5) * t373 + Ifges(7,6) * t379;
t430 = pkin(6) * t317 + t236 * t380;
t185 = qJD(5) * t278 - t360 * t488 + t374 * t463 - t380 * t510;
t277 = t360 * t544 - t380 * t361;
t428 = t155 * t185 + t277 * t35;
t230 = t316 * t374 + t317 * t536;
t162 = -pkin(6) * t230 - t566;
t182 = t430 * t375;
t105 = t162 * t379 - t182 * t373;
t106 = t162 * t373 + t182 * t379;
t423 = t173 * t375 + t174 * t381;
t417 = pkin(1) * t454;
t240 = -t294 * t374 + t378 * t553;
t309 = -t373 * t381 + t375 * t538;
t415 = t373 * t541 + t537;
t414 = t317 * t524 + t565;
t411 = t377 * (Ifges(3,1) * t383 - t600);
t394 = -qJD(6) * t230 + t414;
t296 = t369 * t533 + t540;
t243 = t296 * t380 + t374 * t552;
t389 = -g(1) * t243 - g(2) * t239 - g(3) * t285 + t3 * t379 - t373 * t4;
t159 = Ifges(5,5) * t254 + Ifges(5,6) * t253 + Ifges(5,3) * t299;
t220 = -Ifges(4,2) * t306 + Ifges(4,6) * t371 + t598;
t298 = Ifges(4,4) * t306;
t221 = Ifges(4,1) * t307 + Ifges(4,5) * t371 - t298;
t57 = Ifges(5,1) * t142 + Ifges(5,4) * t143 + Ifges(5,5) * t198;
t385 = (Ifges(7,5) * t208 + Ifges(7,6) * t209) * t642 + (Ifges(7,5) * t167 + Ifges(7,6) * t166) * t643 + (Ifges(5,2) * t648 + Ifges(5,6) * t637 - t742) * t381 + ((Ifges(6,3) * t375 + t381 * t434) * t635 + (Ifges(6,5) * t375 + t381 * t446) * t638 + (Ifges(6,6) * t375 + t381 * t440) * t640 - t702 + t479) * qJD(4) + t596 * t648 + ((-t215 + t400) * mrSges(6,3) - t735 * mrSges(7,2) + t736 * mrSges(7,1) + t734 * mrSges(6,1)) * t155 + (-t3 * t415 - t4 * t309 + t58 * t735 - t59 * t736) * mrSges(7,3) + t174 * mrSges(5,2) * t307 + t558 * t718 + (t726 + t727) * t398 + t595 * t649 + t374 * t112 * t476 + (Ifges(6,1) * t639 + Ifges(6,4) * t641 + Ifges(6,5) * t636 + t730) * t215 + (Ifges(5,3) * t307 - t306 * t435) * t631 + (Ifges(5,5) * t307 - t306 * t447) * t633 + (Ifges(5,6) * t307 - t306 * t441) * t634 - t371 * (-Ifges(4,5) * t306 - Ifges(4,6) * t307) / 0.2e1 + t331 * (mrSges(4,1) * t307 - mrSges(4,2) * t306) + (t173 * t558 + t174 * t734 + t576 + t690) * mrSges(5,3) + t750 * t546 + t57 * t627 + t220 * t628 + t208 * t656 + t167 * t657 + t209 * t658 + t166 * t659 + t541 * t664 + t390 * t604 + (Ifges(7,1) * t208 + Ifges(7,4) * t209) * t644 + (Ifges(7,1) * t167 + Ifges(7,4) * t166) * t645 + t309 * t666 + (t253 * t441 + t254 * t447 + t299 * t435) * qJD(4) / 0.2e1 - (-Ifges(4,1) * t306 + t159 - t598) * t307 / 0.2e1 + (Ifges(7,5) * t309 - Ifges(7,6) * t415) * t655 + (mrSges(7,1) * t415 + mrSges(7,2) * t309 + mrSges(6,3) * t541) * t35 + (Ifges(7,1) * t309 - Ifges(7,4) * t415) * t663 - t415 * t668 + (Ifges(7,4) * t309 - Ifges(7,2) * t415) * t662 + (-Ifges(4,2) * t307 + t221 - t298) * t306 / 0.2e1 + t450 * t575 - t306 * t465 + t307 * t466 + (mrSges(5,1) * t307 + mrSges(6,1) * t398 + mrSges(6,2) * t400) * t173 + (-mrSges(5,2) * t390 + Ifges(5,1) * t649 + Ifges(5,5) * t637 + t112 * t473 + t113 * t474 + t434 * t650 + t440 * t653 + t446 * t654 + t472 * t60) * t375 + (Ifges(7,4) * t208 + Ifges(7,2) * t209) * t646 + (Ifges(7,4) * t167 + Ifges(7,2) * t166) * t647 + (-t433 * t635 - t439 * t640 - t445 * t638) * t522 + (t113 * t380 + t374 * t60 + t161) * t475 + Ifges(4,3) * t370 - t312 * mrSges(4,1) - t313 * mrSges(4,2) + Ifges(4,6) * t202 + Ifges(4,5) * t201 + (-Ifges(6,4) * t639 - Ifges(6,2) * t641 - Ifges(6,6) * t636 + t725 - t727 - t738) * (-t306 * t544 - t380 * t307) - t526 * t739 + (t719 - Ifges(6,3) * t636 - Ifges(6,5) * t639 - Ifges(6,6) * t641 - t739 + t111 / 0.2e1) * t559 - t306 * t702 + t160 * t477;
t364 = Ifges(3,4) * t530;
t341 = pkin(5) * t549;
t339 = pkin(5) * t550;
t305 = Ifges(3,1) * t531 + Ifges(3,5) * qJD(2) + t364;
t304 = Ifges(3,6) * qJD(2) + qJD(1) * t442;
t295 = t369 * t532 - t539;
t269 = -mrSges(4,2) * t371 - mrSges(4,3) * t306;
t242 = t296 * t374 - t380 * t552;
t234 = mrSges(4,1) * t306 + mrSges(4,2) * t307;
t187 = t243 * t379 + t295 * t373;
t186 = -t243 * t373 + t295 * t379;
t179 = t230 * t379 + t317 * t547;
t178 = -t230 * t373 + t317 * t542;
t136 = -pkin(6) * t562 + t174;
t132 = pkin(6) * t254 - t568;
t108 = -t155 * t379 - t373 * t617;
t107 = t155 * t373 - t379 * t617;
t103 = mrSges(7,1) * t191 - mrSges(7,3) * t150;
t102 = -mrSges(7,2) * t191 + mrSges(7,3) * t149;
t81 = -mrSges(5,1) * t143 + mrSges(5,2) * t142;
t80 = t132 * t379 + t136 * t373;
t79 = -t132 * t373 + t136 * t379;
t66 = t430 * t524 + (pkin(6) * t245 + t164 * t380 - t236 * t523) * t375;
t50 = -pkin(6) * t114 - t164 * t381 + t490;
t48 = -t373 * t455 + t379 * t394;
t47 = t373 * t394 + t379 * t455;
t41 = -mrSges(6,2) * t140 + mrSges(6,3) * t78;
t15 = -mrSges(7,2) * t74 + mrSges(7,3) * t30;
t14 = mrSges(7,1) * t74 - mrSges(7,3) * t29;
t13 = -qJD(6) * t106 - t373 * t66 + t379 * t50;
t12 = qJD(6) * t105 + t373 * t50 + t379 * t66;
t1 = [(Ifges(7,4) * t47 + Ifges(7,2) * t48) * t646 + (Ifges(7,4) * t179 + Ifges(7,2) * t178) * t662 + (Ifges(6,5) * t114 + Ifges(6,3) * t414) * t635 - t671 * t164 + (Ifges(6,4) * t114 + Ifges(6,6) * t414) * t640 + (t57 / 0.2e1 + t70 * mrSges(5,3)) * t317 * t381 + t321 * t599 / 0.2e1 + ((m(7) * t35 + t712) * t546 + t400 * t151 + m(5) * (qJD(4) * t423 + t375 * t69 - t381 * t70) + t203 * t524 + t398 * t716 + m(6) * ((qJD(4) * t425 - t70) * t381 + (t527 + (t155 * t380 - t156 * t374) * qJD(5) + t687) * t375) + t41 * t541 + t375 * t101) * t236 + (-t69 * mrSges(5,3) + t742) * t557 + (t506 + t694) * t87 + t442 * t717 + t564 * t718 + t565 * t719 + (mrSges(6,2) * t70 + t35 * mrSges(6,3) + 0.2e1 * t664) * t230 + (mrSges(4,2) * t297 + mrSges(4,3) * t312 + Ifges(4,1) * t201 + Ifges(4,4) * t202 + Ifges(4,5) * t370 + t111 * t475 + t160 * t476 + t161 * t477 - t390 * t452 + t435 * t637 + t441 * t648 + t447 * t649) * t317 + (t305 * t626 - t469 + t436 * qJD(2) / 0.2e1) * qJD(2) + t234 * t512 + (Ifges(7,5) * t47 + Ifges(7,6) * t48) * t642 + (Ifges(7,5) * t179 + Ifges(7,6) * t178) * t655 + (Ifges(7,1) * t47 + Ifges(7,4) * t48) * t644 + (Ifges(7,1) * t179 + Ifges(7,4) * t178) * t663 + (m(7) * t506 - mrSges(6,1) * t414 - mrSges(7,1) * t48 + mrSges(7,2) * t47 + mrSges(6,3) * t114) * t155 + (-Ifges(6,6) * t635 - Ifges(6,4) * t638 - Ifges(6,2) * t640 - t112 / 0.2e1 - t672 + t726 + t732) * (qJD(5) * t230 + t245 * t544 - t380 * t246 - t317 * t488) + (Ifges(6,1) * t114 + Ifges(6,5) * t414) * t638 - t417 * t518 + (mrSges(6,1) * t70 - Ifges(6,4) * t654 - Ifges(6,2) * t653 - Ifges(6,6) * t650 + t750) * (-t380 * t316 + t317 * t544) - t304 * t529 / 0.2e1 + (Ifges(4,1) * t245 - Ifges(4,4) * t246) * t628 + (-Ifges(5,1) * t413 - Ifges(5,4) * t414 + Ifges(5,5) * t246) * t632 + t47 * t656 + t48 * t658 - t245 * t465 - t246 * t466 + t179 * t666 + t178 * t668 - t419 * (mrSges(5,1) * t414 - mrSges(5,2) * t413) + m(7) * (t105 * t4 + t106 * t3 + t12 * t59 + t13 * t58) - t173 * (mrSges(5,1) * t246 + mrSges(5,3) * t413) + t253 * (-Ifges(5,4) * t413 - Ifges(5,2) * t414 + Ifges(5,6) * t246) / 0.2e1 + t299 * (-Ifges(5,5) * t413 - Ifges(5,6) * t414 + Ifges(5,3) * t246) / 0.2e1 + t174 * (-mrSges(5,2) * t246 - mrSges(5,3) * t414) + (t178 * t3 - t179 * t4 - t47 * t58 + t48 * t59) * mrSges(7,3) + Ifges(2,3) * qJDD(1) + t371 * (Ifges(4,5) * t245 - Ifges(4,6) * t246) / 0.2e1 - t362 * (-mrSges(4,1) * t202 + mrSges(4,2) * t201) - t331 * (mrSges(4,1) * t246 + mrSges(4,2) * t245) - t306 * (Ifges(4,4) * t245 - Ifges(4,2) * t246) / 0.2e1 + (t114 * t173 - t156 * t414) * mrSges(6,2) + t245 * t221 / 0.2e1 + t246 * t159 / 0.2e1 - t246 * t220 / 0.2e1 + t35 * (-mrSges(7,1) * t178 + mrSges(7,2) * t179) + (t599 * t468 + t685 / 0.2e1) * t383 + (Ifges(3,4) * t321 + Ifges(3,2) * t320 + t685) * t626 + (mrSges(6,1) * t239 + mrSges(5,1) * t294 + mrSges(5,3) * t554 - t728 * mrSges(7,1) - t729 * mrSges(7,2) + t713 * t293 + t737 * t240 + t682 * t384 + (t679 + t683 * (-t362 - t315)) * t378) * g(1) + (-t296 * mrSges(5,1) - mrSges(6,1) * t243 - t187 * mrSges(7,1) - t186 * mrSges(7,2) - mrSges(5,3) * t552 - t713 * t295 + t683 * (pkin(3) * t549 + pkin(5) * t552 + t384 * t362) + t737 * t242 - t679 * t384 + t682 * t378) * g(2) - t694 * t152 + t701 * t490 + (m(3) * t574 + mrSges(3,1) * t320 - mrSges(3,2) * t321) * pkin(1) - t706 * t566 - t330 * t574 + t411 * t468 + t245 * t479 + (Ifges(5,3) * t637 + Ifges(5,6) * t648 + Ifges(5,5) * t649 + t498 / 0.2e1 - Ifges(4,4) * t201 - Ifges(4,2) * t202 - Ifges(4,6) * t370 + t297 * mrSges(4,1) - t70 * mrSges(5,1) - t714 - t313 * mrSges(4,3)) * t316 + (Ifges(3,1) * t321 + Ifges(3,4) * t717 + Ifges(3,5) * qJDD(2) - t468 * t709) * t377 + t12 * t102 + t13 * t103 + t105 * t14 + t106 * t15 + t114 * t113 / 0.2e1 - t297 * t625; m(6) * (t156 * t184 + t278 * t34 + t428) + m(5) * ((-pkin(2) * t517 - t608) * t361 + (t423 * t382 + (-t607 + (t361 - t622) * qJD(2)) * t376) * t582) + t671 * t219 + (t454 + t453 + t516) * t686 + (t683 * (t378 * t456 + t339) + t673) * g(2) + (t683 * (t384 * t456 + t341) + t674) * g(1) + (-m(6) * t173 - t701) * (-t360 * t524 - t464) + m(4) * (-t312 * t382 + t313 * t376) * pkin(2) + t269 * t509 + ((m(6) * t70 + t706) * t375 + t731) * t360 - t234 * t514 + (t469 + (t417 - t411 / 0.2e1) * qJD(1)) * qJD(1) - t436 * t518 / 0.2e1 + t385 + t304 * t531 / 0.2e1 + (mrSges(4,1) * t370 - mrSges(4,3) * t201) * t622 + (-mrSges(4,2) * t370 + mrSges(4,3) * t202) * t624 - t87 * t505 + t203 * t463 + Ifges(3,3) * qJDD(2) + t361 * t81 + Ifges(3,6) * t320 + Ifges(3,5) * t321 + t278 * t41 + t199 * t14 + t200 * t15 + t184 * t151 - t698 * t510 - (-Ifges(3,2) * t531 + t305 + t364) * t530 / 0.2e1 + t707 * t185 + t710 * t102 + t711 * t103 + (-t155 * t505 + t199 * t4 + t200 * t3 + t711 * t58 + t710 * t59 + t428) * m(7) + t712 * t277 + (t330 - m(4) * t621 - m(7) * (t458 + t621) + t722 * (t315 + t621) + t677) * g(3); (t683 * (-pkin(3) * t554 + t339) + t673) * g(2) + (t683 * (-pkin(3) * t552 + t341) + t674) * g(1) + (pkin(3) * t390 - t174 * t206 + t419 * t513 - t570) * m(5) + t385 - t269 * t511 + t319 * t41 + t231 * t14 + t232 * t15 - t206 * t203 + t686 * t453 + (t156 * t697 + t319 * t34 - t570 + t695) * m(6) + t697 * t151 + t698 * t513 + t703 * t102 + t704 * t103 + (t231 * t4 + t232 * t3 + t58 * t704 + t59 * t703 + t695) * m(7) + t712 * t318 - pkin(3) * t81 + (-m(7) * t458 + t315 * t722 + t677) * g(3) + t707 * t743 + t701 * (pkin(5) * t524 - t205) + (m(6) * t690 + t375 * t706 + t731) * pkin(5); (-m(6) * (t174 - t425) + t203 - t602 + t250 * t450) * t173 + (t295 * t680 + t296 * t681) * g(1) + (t293 * t680 + t294 * t681) * g(2) + (t102 * t692 + t103 * t693 - t14 * t545 - t15 * t548) * pkin(6) - m(7) * (t58 * t79 + t59 * t80) + t380 * t721 + (-t675 + t725) * t563 + (Ifges(7,5) * t171 + Ifges(7,6) * t170) * t643 + (Ifges(7,1) * t171 + Ifges(7,4) * t170) * t645 + (mrSges(7,1) * t309 - mrSges(7,2) * t415 - mrSges(6,3) * t381 + t452 + (-t451 + t457) * t375) * g(3) * t368 - t714 + (t675 + t732) * t523 + (Ifges(5,1) * t253 + t111 - t597) * t633 + (t577 + t113) * t472 + (Ifges(5,5) * t253 - Ifges(5,6) * t254) * t631 + t160 * t632 + (Ifges(6,3) * t254 + t253 * t434) * t636 + (Ifges(6,5) * t254 + t253 * t446) * t639 + (Ifges(6,6) * t254 + t253 * t440) * t641 + (-t431 * t519 + (Ifges(7,3) * t374 + t380 * t432) * qJD(5)) * t642 + (-t443 * t519 + (Ifges(7,5) * t374 + t380 * t444) * qJD(5)) * t644 + (-t437 * t519 + (Ifges(7,6) * t374 + t380 * t438) * qJD(5)) * t646 + t433 * t650 + t562 * t651 + t439 * t653 + t445 * t654 + (-Ifges(7,3) * t380 + t374 * t432) * t655 + t171 * t657 + t170 * t659 + (-Ifges(7,6) * t380 + t374 * t438) * t662 + (-Ifges(7,5) * t380 + t374 * t444) * t663 + t374 * t664 + t254 * t739 + (-Ifges(5,2) * t254 + t161 + t249) * t634 + t498 + t545 * t666 + (-t429 * t521 + (-t3 * t373 - t379 * t4 + (t373 * t58 - t379 * t59) * qJD(6)) * t374) * t667 + t151 * t568 + t473 * t578 + t448 * t579 + t419 * (mrSges(5,1) * t254 + mrSges(5,2) * t253) + t380 * t20 / 0.2e1 - t380 * t8 / 0.2e1 - t9 * t548 / 0.2e1 + (t156 * t563 + t687) * mrSges(6,3) + (m(7) * t569 + mrSges(6,1) * t254 + (t521 - t562) * mrSges(6,3) - t699 * mrSges(7,2) + t700 * mrSges(7,1)) * t155 + (-t3 * t548 - t4 * t545 + t58 * t699 - t59 * t700) * mrSges(7,3) + (t601 - t701) * t174 - (t373 * t62 + t379 * t61) * t519 / 0.2e1 + (t193 * t440 + t194 * t446 + t250 * t434) * qJD(5) / 0.2e1 + t705 * t70 + t707 * t569 + t112 * t474 - t80 * t102 - t79 * t103 - t380 * t720 + (Ifges(7,4) * t171 + Ifges(7,2) * t170) * t647; (m(7) * t389 - t373 * t14 + t379 * t15) * pkin(6) + t389 * mrSges(7,3) + t19 + t443 * t663 + t431 * t655 + t437 * t662 + (-t584 / 0.2e1 - t585 / 0.2e1 - t587 / 0.2e1 + t586 / 0.2e1 + t581 + t594 / 0.2e1 + t672 + t733) * t194 + (mrSges(6,2) * t243 + t242 * t418) * g(1) - t418 * t35 + (mrSges(6,2) * t285 + t416 * t418) * g(3) + (mrSges(6,2) * t239 - t240 * t418) * g(2) + t379 * t668 + t373 * t666 + (-t588 / 0.2e1 - t190 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t194 + t708 * t155 - t678 + t730) * t193 + (t155 * t448 + (-m(7) * t429 - t373 * t102 - t379 * t103) * pkin(6) + t678) * qJD(6) - m(7) * (t107 * t58 + t108 * t59) - t715 - t107 * t103 - t108 * t102 + t155 * t151 + (-t707 - t716) * t156; -t155 * (mrSges(7,1) * t150 + mrSges(7,2) * t149) + (Ifges(7,1) * t149 - t591) * t645 + t61 * t644 + (Ifges(7,5) * t149 - Ifges(7,6) * t150) * t643 - t58 * t102 + t59 * t103 - g(1) * (mrSges(7,1) * t186 - mrSges(7,2) * t187) - g(2) * (-mrSges(7,1) * t729 + mrSges(7,2) * t728) - g(3) * ((-t285 * t373 + t502) * mrSges(7,1) + (-t285 * t379 - t368 * t547) * mrSges(7,2)) + (t149 * t58 + t150 * t59) * mrSges(7,3) + t8 + (-Ifges(7,2) * t150 + t144 + t62) * t647 + t676;];
tau  = t1;
