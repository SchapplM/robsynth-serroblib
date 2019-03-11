% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:40
% EndTime: 2019-03-08 22:17:02
% DurationCPUTime: 14.15s
% Computational Cost: add. (27868->739), mult. (64780->1043), div. (0->0), fcn. (73159->12), ass. (0->364)
t460 = sin(qJ(6));
t464 = cos(qJ(6));
t457 = sin(pkin(12));
t459 = cos(pkin(12));
t461 = sin(qJ(5));
t622 = cos(qJ(5));
t490 = t461 * t457 - t622 * t459;
t491 = t457 * t622 + t461 * t459;
t336 = -t460 * t490 + t464 * t491;
t512 = -t460 * t491 - t464 * t490;
t536 = Ifges(7,5) * t512 - Ifges(7,6) * t336;
t613 = pkin(9) + qJ(4);
t430 = t613 * t457;
t432 = t613 * t459;
t360 = -t461 * t430 + t432 * t622;
t301 = -pkin(10) * t490 + t360;
t359 = -t622 * t430 - t432 * t461;
t499 = -pkin(10) * t491 + t359;
t185 = t301 * t464 + t460 * t499;
t683 = -t301 * t460 + t464 * t499;
t736 = -t185 * mrSges(7,1) - t683 * mrSges(7,2);
t38 = t536 + t736;
t741 = t38 * qJD(6);
t462 = sin(qJ(3));
t391 = t490 * t462;
t392 = t491 * t462;
t513 = t391 * t460 - t464 * t392;
t291 = t391 * t464 + t392 * t460;
t714 = t291 * mrSges(7,1);
t725 = t513 * mrSges(7,2) - t714;
t740 = t725 / 0.2e1;
t713 = t336 * mrSges(7,1);
t726 = t512 * mrSges(7,2) + t713;
t739 = t726 / 0.2e1;
t465 = cos(qJ(3));
t429 = -pkin(3) * t465 - qJ(4) * t462 - pkin(2);
t410 = t459 * t429;
t546 = t459 * t462;
t341 = -pkin(9) * t546 + t410 + (-pkin(8) * t457 - pkin(4)) * t465;
t545 = t459 * t465;
t370 = pkin(8) * t545 + t457 * t429;
t550 = t457 * t462;
t355 = -pkin(9) * t550 + t370;
t222 = t622 * t341 - t355 * t461;
t189 = pkin(10) * t391 + t222;
t168 = -pkin(5) * t465 + t189;
t223 = t461 * t341 + t355 * t622;
t190 = -t392 * pkin(10) + t223;
t564 = t190 * t464;
t78 = t168 * t460 + t564;
t91 = -t189 * t460 - t564;
t733 = t78 + t91;
t625 = t462 / 0.2e1;
t565 = t190 * t460;
t77 = t168 * t464 - t565;
t92 = t189 * t464 - t565;
t738 = -t77 + t92;
t451 = t462 * pkin(8);
t421 = pkin(4) * t550 + t451;
t344 = pkin(5) * t392 + t421;
t737 = t344 * t725;
t607 = Ifges(7,4) * t336;
t214 = Ifges(7,2) * t512 + t607;
t728 = Ifges(7,1) * t512 - t607;
t735 = t214 / 0.4e1 - t728 / 0.4e1;
t608 = Ifges(7,4) * t291;
t160 = Ifges(7,2) * t513 - t465 * Ifges(7,6) - t608;
t282 = Ifges(7,4) * t513;
t162 = -Ifges(7,1) * t291 - t465 * Ifges(7,5) + t282;
t447 = -pkin(4) * t459 - pkin(3);
t378 = pkin(5) * t490 + t447;
t251 = mrSges(7,2) * t465 + mrSges(7,3) * t513;
t665 = t251 / 0.2e1;
t708 = Ifges(7,2) * t291 + t282;
t727 = Ifges(7,1) * t513 + t608;
t734 = t683 * t665 + (t162 / 0.4e1 + t708 / 0.4e1) * t512 + t378 * t740 + t344 * t739 + (t727 / 0.4e1 - t160 / 0.4e1) * t336;
t526 = t713 / 0.2e1;
t528 = -t714 / 0.2e1;
t458 = sin(pkin(6));
t463 = sin(qJ(2));
t548 = t458 * t463;
t566 = cos(pkin(6));
t402 = t462 * t566 + t465 * t548;
t466 = cos(qJ(2));
t547 = t458 * t466;
t353 = -t402 * t457 - t459 * t547;
t354 = t402 * t459 - t457 * t547;
t229 = t353 * t622 - t461 * t354;
t230 = t461 * t353 + t354 * t622;
t119 = t229 * t460 + t230 * t464;
t514 = t464 * t229 - t230 * t460;
t35 = -t119 * mrSges(7,1) - t514 * mrSges(7,2);
t731 = t35 * qJD(6);
t401 = t462 * t548 - t465 * t566;
t288 = t491 * t401;
t289 = t490 * t401;
t158 = t288 * t464 - t289 * t460;
t159 = t288 * t460 + t289 * t464;
t672 = m(7) * pkin(5);
t531 = -t672 / 0.2e1;
t720 = -mrSges(7,2) / 0.2e1;
t721 = mrSges(7,1) / 0.2e1;
t539 = t158 * t721 + t159 * t720;
t670 = mrSges(6,2) / 0.2e1;
t671 = -mrSges(6,1) / 0.2e1;
t730 = t288 * t671 + t289 * t670 + (t158 * t464 + t159 * t460) * t531 - t539;
t541 = t465 * t466;
t367 = (-t457 * t541 + t459 * t463) * t458;
t368 = (t457 * t463 + t459 * t541) * t458;
t255 = t367 * t622 - t461 * t368;
t256 = t461 * t367 + t368 * t622;
t149 = t255 * t464 - t256 * t460;
t150 = t255 * t460 + t256 * t464;
t540 = t149 * t721 + t150 * t720;
t729 = t255 * t671 + t256 * t670 + (t149 * t464 + t150 * t460) * t531 - t540;
t697 = Ifges(7,5) * t513;
t717 = Ifges(7,6) * t291;
t538 = t697 + t717;
t393 = t491 * t465;
t394 = t490 * t465;
t293 = -t393 * t464 + t394 * t460;
t296 = -t393 * t460 - t394 * t464;
t496 = -Ifges(7,5) * t296 / 0.2e1 - Ifges(7,6) * t293 / 0.2e1;
t434 = pkin(3) * t462 - qJ(4) * t465;
t379 = pkin(8) * t550 + t459 * t434;
t343 = pkin(4) * t462 - pkin(9) * t545 + t379;
t380 = -pkin(8) * t546 + t457 * t434;
t549 = t457 * t465;
t356 = -pkin(9) * t549 + t380;
t225 = t622 * t343 - t356 * t461;
t186 = pkin(5) * t462 + pkin(10) * t394 + t225;
t226 = t461 * t343 + t622 * t356;
t193 = -pkin(10) * t393 + t226;
t89 = t186 * t464 - t193 * t460;
t90 = t186 * t460 + t193 * t464;
t482 = t496 + t90 * mrSges(7,2) / 0.2e1 - t89 * mrSges(7,1) / 0.2e1;
t510 = Ifges(7,3) * t625 - t482;
t329 = Ifges(7,4) * t512;
t217 = Ifges(7,1) * t336 + t329;
t707 = -Ifges(7,2) * t336 + t329;
t724 = t707 / 0.4e1 + t217 / 0.4e1;
t517 = -t717 / 0.2e1 - t697 / 0.2e1;
t722 = 0.2e1 * mrSges(7,2);
t719 = t119 / 0.2e1;
t718 = -t336 / 0.2e1;
t702 = -t513 / 0.2e1;
t711 = t514 * t702;
t529 = t462 * t547;
t553 = t368 * t459;
t554 = t367 * t457;
t675 = -m(7) / 0.2e1;
t677 = -m(6) / 0.2e1;
t679 = -m(5) / 0.2e1;
t709 = -(t553 / 0.2e1 - t554 / 0.2e1) * mrSges(5,3) + (-pkin(3) * t529 + (t553 - t554) * qJ(4)) * t679 + (t359 * t255 + t360 * t256 + t447 * t529) * t677 + (t149 * t683 + t185 * t150 + t378 * t529) * t675;
t705 = t229 / 0.2e1;
t704 = -t512 / 0.2e1;
t703 = t512 / 0.2e1;
t701 = t513 / 0.2e1;
t693 = t514 * t665;
t302 = -t391 * mrSges(6,1) - t392 * mrSges(6,2);
t692 = t725 + t302;
t347 = mrSges(6,1) * t491 - mrSges(6,2) * t490;
t691 = t726 + t347;
t211 = -mrSges(7,1) * t512 + mrSges(7,2) * t336;
t348 = mrSges(6,1) * t490 + mrSges(6,2) * t491;
t690 = t211 + t348;
t689 = -Ifges(6,5) * t490 - Ifges(6,6) * t491 + t536;
t688 = -Ifges(6,5) * t392 + Ifges(6,6) * t391 + t538;
t687 = -t379 * t457 + t380 * t459;
t568 = t465 * mrSges(4,2);
t435 = t462 * mrSges(4,1) + t568;
t609 = Ifges(6,4) * t491;
t350 = -Ifges(6,2) * t490 + t609;
t351 = -Ifges(6,1) * t490 - t609;
t686 = -t351 / 0.4e1 + t350 / 0.4e1;
t176 = -mrSges(7,1) * t513 - mrSges(7,2) * t291;
t630 = t491 / 0.2e1;
t637 = -t391 / 0.2e1;
t674 = m(7) / 0.2e1;
t682 = (t344 * t491 - t378 * t391) * t674 + t176 * t630 + t211 * t637;
t610 = Ifges(6,4) * t391;
t284 = -Ifges(6,2) * t392 - t465 * Ifges(6,6) - t610;
t388 = Ifges(6,4) * t392;
t286 = -Ifges(6,1) * t391 - t465 * Ifges(6,5) - t388;
t305 = Ifges(6,2) * t391 - t388;
t306 = -Ifges(6,1) * t392 + t610;
t408 = Ifges(6,4) * t490;
t349 = -Ifges(6,2) * t491 - t408;
t352 = Ifges(6,1) * t491 - t408;
t623 = -t465 / 0.4e1;
t361 = mrSges(6,2) * t465 - t392 * mrSges(6,3);
t640 = t361 / 0.2e1;
t643 = t347 / 0.2e1;
t653 = t302 / 0.2e1;
t661 = t291 / 0.2e1;
t253 = -mrSges(7,1) * t465 + t291 * mrSges(7,3);
t663 = t253 / 0.2e1;
t681 = -(t305 / 0.4e1 + t286 / 0.4e1) * t490 - (t352 / 0.4e1 + t349 / 0.4e1 - t359 * mrSges(6,3) / 0.2e1) * t392 + (t683 * t702 + t703 * t92 + t704 * t77 + t718 * t733) * mrSges(7,3) + t359 * t640 + t421 * t643 + t447 * t653 + t689 * t623 + t724 * t513 + t733 * t683 * t674 - (-t306 / 0.4e1 + t284 / 0.4e1) * t491 + t735 * t291 + (t661 * mrSges(7,3) + t674 * t738 - t663) * t185 + t734;
t454 = t459 ^ 2;
t680 = 2 * qJD(3);
t678 = m(5) / 0.2e1;
t676 = m(6) / 0.2e1;
t669 = -mrSges(7,3) / 0.2e1;
t668 = -t119 / 0.2e1;
t664 = -t253 / 0.2e1;
t658 = t293 / 0.2e1;
t655 = -t291 / 0.2e1;
t654 = t296 / 0.2e1;
t645 = t336 / 0.2e1;
t578 = t391 * mrSges(6,3);
t363 = -mrSges(6,1) * t465 + t578;
t639 = -t363 / 0.2e1;
t636 = t391 / 0.2e1;
t635 = -t392 / 0.2e1;
t634 = -t393 / 0.2e1;
t633 = -t394 / 0.2e1;
t632 = t401 / 0.2e1;
t631 = t490 / 0.2e1;
t629 = -t457 / 0.2e1;
t628 = t457 / 0.2e1;
t627 = t459 / 0.2e1;
t626 = -t460 / 0.2e1;
t624 = -t465 / 0.2e1;
t621 = pkin(5) * t391;
t620 = pkin(5) * t491;
t452 = t465 * pkin(8);
t619 = t77 * mrSges(7,2);
t618 = t78 * mrSges(7,1);
t615 = t91 * mrSges(7,1);
t614 = t92 * mrSges(7,2);
t612 = Ifges(5,4) * t457;
t611 = Ifges(5,4) * t459;
t606 = Ifges(5,5) * t459;
t604 = Ifges(5,6) * t457;
t601 = pkin(5) * qJD(5);
t587 = t293 * mrSges(7,1);
t584 = t296 * mrSges(7,2);
t582 = t512 * mrSges(7,3);
t579 = t336 * mrSges(7,3);
t577 = t393 * mrSges(6,1);
t576 = t394 * mrSges(6,2);
t575 = t490 * mrSges(6,3);
t574 = t491 * mrSges(6,3);
t573 = t457 * mrSges(5,1);
t572 = t457 * Ifges(5,2);
t571 = t459 * mrSges(5,2);
t569 = t462 * mrSges(4,2);
t431 = -mrSges(5,1) * t459 + mrSges(5,2) * t457;
t567 = -mrSges(4,1) + t431;
t563 = t291 * t464;
t562 = t513 * t460;
t365 = t401 * t529;
t30 = m(7) * (t119 * t150 + t149 * t514 + t365) + m(6) * (t229 * t255 + t230 * t256 + t365) + m(5) * (t353 * t367 + t354 * t368 + t365) + m(4) * (t401 * t462 + t402 * t465 - t548) * t547;
t561 = t30 * qJD(1);
t342 = t401 * t402;
t502 = t457 * t353 - t459 * t354;
t31 = m(7) * (t119 * t159 + t158 * t514 + t342) + m(6) * (t229 * t288 + t230 * t289 + t342) + m(5) * (t401 * t502 + t342);
t558 = t31 * qJD(1);
t557 = t336 * t464;
t556 = t512 * t460;
t252 = -mrSges(7,2) * t462 + mrSges(7,3) * t293;
t544 = t460 * t252;
t254 = mrSges(7,1) * t462 - mrSges(7,3) * t296;
t542 = t464 * t254;
t422 = pkin(4) * t549 + t452;
t532 = t457 ^ 2 + t454;
t530 = t672 / 0.2e1;
t527 = t582 / 0.2e1;
t525 = -t579 / 0.2e1;
t524 = -t575 / 0.2e1;
t523 = -t574 / 0.2e1;
t520 = t726 * t632;
t509 = t571 + t573;
t177 = t584 - t587;
t303 = mrSges(6,1) * t392 - mrSges(6,2) * t391;
t304 = -t576 + t577;
t345 = pkin(5) * t393 + t422;
t362 = -mrSges(6,2) * t462 - mrSges(6,3) * t393;
t364 = mrSges(6,1) * t462 + mrSges(6,3) * t394;
t403 = t509 * t462;
t404 = t509 * t465;
t418 = -mrSges(5,2) * t462 - mrSges(5,3) * t549;
t420 = mrSges(5,1) * t462 - mrSges(5,3) * t545;
t417 = t465 * mrSges(5,2) - mrSges(5,3) * t550;
t419 = -t465 * mrSges(5,1) - mrSges(5,3) * t546;
t492 = t419 * t628 - t459 * t417 / 0.2e1;
t369 = -pkin(8) * t549 + t410;
t501 = -t369 * t457 + t370 * t459;
t467 = (t303 / 0.2e1 + t176 / 0.2e1 + t403 / 0.2e1) * t402 + (t404 / 0.2e1 + t304 / 0.2e1 + t177 / 0.2e1 + t492) * t401 + (t402 * t451 + t353 * t379 + t354 * t380 + (-t501 + t452) * t401) * t678 + (t222 * t288 + t223 * t289 + t225 * t229 + t226 * t230 + t401 * t422 + t402 * t421) * t676 + (t119 * t90 + t158 * t77 + t159 * t78 + t344 * t402 + t345 * t401 + t514 * t89) * t674 + t514 * t254 / 0.2e1 + t252 * t719 + t158 * t663 + t159 * t665 + t364 * t705 + t230 * t362 / 0.2e1 + t288 * t363 / 0.2e1 + t289 * t640 + t353 * t420 / 0.2e1 + t354 * t418 / 0.2e1;
t3 = t467 + (-t435 / 0.2e1 + t568 / 0.2e1 + (mrSges(4,1) / 0.2e1 - t431 / 0.2e1 - t348 / 0.2e1 - t211 / 0.2e1) * t462) * t547 + (t149 * t645 + t150 * t704) * mrSges(7,3) + (t255 * t630 + t256 * t631) * mrSges(6,3) + t709;
t161 = Ifges(7,4) * t296 + Ifges(7,2) * t293 + Ifges(7,6) * t462;
t163 = Ifges(7,1) * t296 + Ifges(7,4) * t293 + Ifges(7,5) * t462;
t285 = -Ifges(6,4) * t394 - Ifges(6,2) * t393 + Ifges(6,6) * t462;
t287 = -Ifges(6,1) * t394 - Ifges(6,4) * t393 + Ifges(6,5) * t462;
t389 = Ifges(5,6) * t462 + (-t572 + t611) * t465;
t390 = Ifges(5,5) * t462 + (t459 * Ifges(5,1) - t612) * t465;
t497 = Ifges(6,5) * t394 / 0.2e1 + Ifges(6,6) * t393 / 0.2e1;
t5 = t286 * t633 + t284 * t634 + t285 * t635 + t287 * t637 + (Ifges(7,5) * t655 + Ifges(7,6) * t701 + Ifges(6,5) * t637 + Ifges(6,6) * t635 + t390 * t627 + t389 * t629 + pkin(8) * t404 + (-Ifges(4,4) + t606 / 0.2e1 - t604 / 0.2e1) * t462) * t462 + t161 * t701 + t421 * t304 + t422 * t303 - pkin(2) * t435 + t380 * t417 + t370 * t418 + t379 * t419 + t369 * t420 + t226 * t361 + t223 * t362 + t225 * t363 + t222 * t364 + t344 * t177 + t345 * t176 + t77 * t254 + t90 * t251 + t78 * t252 + t89 * t253 + m(5) * (t369 * t379 + t370 * t380) + (pkin(8) * t403 + (Ifges(4,1) - Ifges(4,2) - Ifges(5,3) - Ifges(6,3) - Ifges(7,3) + t454 * Ifges(5,1) / 0.2e1 + m(5) * pkin(8) ^ 2 + (-t611 + t572 / 0.2e1) * t457) * t462 + t496 + t497 + (Ifges(4,4) + t604 - t606) * t465) * t465 + m(6) * (t222 * t225 + t223 * t226 + t421 * t422) + m(7) * (t344 * t345 + t77 * t89 + t78 * t90) + t162 * t654 + t163 * t655 + t160 * t658;
t508 = t3 * qJD(1) + t5 * qJD(2);
t10 = t92 * t251 - t176 * t621 + t737 + t727 * t655 + t160 * t661 + t91 * t253 + m(7) * (-t344 * t621 + t77 * t91 + t78 * t92) + t306 * t637 + t421 * t302 + t222 * t361 + t284 * t636 + (-t363 + t578) * t223 + (t291 * t78 - t513 * t77) * mrSges(7,3) - (t286 / 0.2e1 - t222 * mrSges(6,3) + t305 / 0.2e1) * t392 + t688 * t624 + (t162 + t708) * t701;
t470 = (t119 * t661 + t711) * mrSges(7,3) + (t230 * t636 + t392 * t705) * mrSges(6,3) + (t119 * t738 - t401 * t621 + t733 * t514) * t674 - t119 * t663 + t693 + t229 * t640 + t230 * t639;
t8 = (t740 + t653) * t401 + t470 + t729;
t507 = t8 * qJD(1) + t10 * qJD(2);
t13 = t77 * t251 + t737 + t538 * t624 - t78 * t253 - (-t78 * mrSges(7,3) + t727 / 0.2e1 - t160 / 0.2e1) * t291 + (t162 / 0.2e1 + t708 / 0.2e1 - t77 * mrSges(7,3)) * t513;
t477 = (-t291 * t668 + t711) * mrSges(7,3) + t693 + t119 * t664 + t725 * t632;
t14 = t477 - t540;
t506 = t14 * qJD(1) + t13 * qJD(2);
t29 = t513 * t251 + t291 * t253 - t392 * t361 + t391 * t363 + m(7) * (t291 * t77 + t513 * t78) + m(6) * (t222 * t391 - t223 * t392) + (m(5) * (-t369 * t459 - t370 * t457) - t457 * t417 - t459 * t419) * t462;
t478 = (t119 * t513 + t291 * t514) * t675 + (t229 * t391 - t230 * t392) * t677;
t488 = (t674 + t676 + t678) * t547;
t489 = m(5) * (-t353 * t459 - t354 * t457);
t40 = (-t489 / 0.2e1 + t488) * t462 + t478;
t505 = -qJD(1) * t40 + qJD(2) * t29;
t57 = (t637 - t563 / 0.2e1 - t562 / 0.2e1) * t672 + t692;
t67 = (t630 + t557 / 0.2e1 - t556 / 0.2e1) * t672 + t691;
t504 = qJD(2) * t57 + qJD(3) * t67;
t106 = t703 * t722 + 0.2e1 * t526;
t83 = t701 * t722 + 0.2e1 * t528;
t500 = qJD(2) * t83 + qJD(3) * t106;
t498 = 0.2e1 * (m(7) / 0.4e1 + m(6) / 0.4e1 + m(5) / 0.4e1) * t402;
t493 = m(7) * (t460 * t90 + t464 * t89);
t473 = t401 * t620 * t674 + (t704 + t703) * mrSges(7,3) * t514;
t11 = (t739 + t643) * t401 + t473 + t730;
t16 = t214 * t718 + t728 * t645 + t211 * t620 + t351 * t630 - t491 * t350 / 0.2e1 + t447 * t347 - (t352 / 0.2e1 + t349 / 0.2e1) * t490 + (t707 + t217) * t703 + (m(7) * t620 + t726) * t378;
t481 = t225 * t671 + t226 * t670 + t497;
t2 = (-t542 / 0.2e1 - t544 / 0.2e1 - t493 / 0.2e1 + t682) * pkin(5) + (t360 * mrSges(6,3) / 0.2e1 + t686) * t391 + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1) * t462 + t360 * t639 + t481 + t482 + t681;
t487 = t11 * qJD(1) + t2 * qJD(2) + t16 * qJD(3);
t19 = t520 - t539;
t27 = t378 * t726 + (-t214 / 0.2e1 + t728 / 0.2e1) * t336 + (t707 / 0.2e1 + t217 / 0.2e1) * t512;
t471 = (t683 * t669 + t724) * t513 - (t185 * t669 - t735) * t291 + t185 * t664 + t536 * t623 + t734;
t7 = t471 - t510;
t486 = t19 * qJD(1) + t7 * qJD(2) + t27 * qJD(3);
t469 = t501 * t679 + (-t222 * t491 - t223 * t490 + t359 * t391 - t360 * t392) * t677 + (t185 * t513 + t291 * t683 - t336 * t77 + t512 * t78) * t675 - t336 * t664 + t251 * t704 + t361 * t631 + t363 * t630 + t492;
t474 = (pkin(8) * t678 + t571 / 0.2e1 + t573 / 0.2e1) * t465 + t422 * t676 + t345 * t674 - t587 / 0.2e1 + t584 / 0.2e1 + t577 / 0.2e1 - t576 / 0.2e1;
t17 = (t291 * t645 + t513 * t704) * mrSges(7,3) + (t391 * t630 - t392 * t631) * mrSges(6,3) + t469 + t474;
t475 = t502 * t678 + (-t229 * t491 - t230 * t490) * t677 + (t119 * t512 - t336 * t514) * t675;
t42 = t498 + t475;
t46 = (t336 ^ 2 + t512 ^ 2) * mrSges(7,3) + (t490 ^ 2 + t491 ^ 2) * mrSges(6,3) + m(7) * (t185 * t512 - t336 * t683) + m(6) * (-t359 * t491 - t360 * t490) + (m(5) * qJ(4) + mrSges(5,3)) * t532;
t485 = qJD(1) * t42 + qJD(2) * t17 - qJD(3) * t46;
t476 = (t464 * t665 + t253 * t626 + (-t291 * t626 + t464 * t702) * mrSges(7,3)) * pkin(5) - t517;
t22 = (-t77 / 0.2e1 + t92 / 0.2e1) * mrSges(7,2) + (-t78 / 0.2e1 - t91 / 0.2e1) * mrSges(7,1) + t476 + t517;
t34 = (t668 + t719) * mrSges(7,1);
t425 = (mrSges(7,1) * t460 + mrSges(7,2) * t464) * pkin(5);
t483 = -qJD(1) * t34 - qJD(2) * t22 + qJD(5) * t425;
t456 = t465 ^ 2;
t455 = t462 ^ 2;
t433 = t455 * pkin(8) * t547;
t414 = t425 * qJD(6);
t157 = -t491 * t531 + (t556 - t557) * t530;
t132 = t391 * t531 + (t562 + t563) * t530;
t107 = t526 - t713 / 0.2e1;
t84 = t528 + t714 / 0.2e1;
t43 = t498 - t475;
t41 = t462 * t488 + t489 * t625 - t478;
t21 = -t619 / 0.2e1 - t618 / 0.2e1 - t614 / 0.2e1 + t615 / 0.2e1 + t476 - t517;
t20 = t520 + t539;
t18 = t291 * t525 + t391 * t523 - t392 * t524 + t513 * t527 - t469 + t474;
t15 = t477 + t540;
t12 = t632 * t691 + t473 - t730;
t9 = t632 * t692 + t470 - t729;
t6 = t471 + t510;
t4 = t149 * t525 + t150 * t527 + t255 * t523 + t256 * t524 + t467 - t435 * t547 + (t431 + t690) * t529 / 0.2e1 - t709;
t1 = -t481 + t510 + Ifges(6,3) * t625 + (t578 / 0.2e1 + t639) * t360 + t686 * t391 + t682 * pkin(5) + (t493 + t542 + t544) * pkin(5) / 0.2e1 + t681;
t23 = [t30 * qJD(2) + t31 * qJD(3), t4 * qJD(3) + t41 * qJD(4) + t9 * qJD(5) + t15 * qJD(6) + t561 + (t149 * t253 + t150 * t251 + t255 * t363 + t256 * t361 + t367 * t419 + t368 * t417 + ((-t465 * mrSges(4,1) - mrSges(3,1) + t569) * t463 + (-mrSges(3,2) + (t455 + t456) * mrSges(4,3) + (t176 + t303 + t403) * t462) * t466) * t458 + 0.2e1 * (t77 * t149 + t78 * t150 + t344 * t529) * t674 + 0.2e1 * (t222 * t255 + t223 * t256 + t421 * t529) * t676 + 0.2e1 * (t367 * t369 + t368 * t370 + t433) * t678 + m(4) * (t433 + (pkin(8) * t456 * t466 - pkin(2) * t463) * t458)) * qJD(2), t558 + t4 * qJD(2) + t43 * qJD(4) + t12 * qJD(5) + t20 * qJD(6) + ((t158 * t683 + t159 * t185 + t378 * t402) * t674 + (t288 * t359 + t289 * t360 + t402 * t447) * t676 + (-qJ(4) * t401 * t532 - pkin(3) * t402) * t678) * t680 + (-t158 * t579 + t159 * t582 - t288 * t574 - t289 * t575 + (-t532 * mrSges(5,3) + mrSges(4,2)) * t401 + (t567 + t690) * t402) * qJD(3), qJD(2) * t41 + qJD(3) * t43, t9 * qJD(2) + t12 * qJD(3) + (-t230 * mrSges(6,1) - t229 * mrSges(6,2) + (-t119 * t464 + t460 * t514) * t672 + t35) * qJD(5) + t731, t15 * qJD(2) + t20 * qJD(3) + t35 * qJD(5) + t731; qJD(3) * t3 - qJD(4) * t40 + qJD(5) * t8 + qJD(6) * t14 - t561, qJD(3) * t5 + qJD(4) * t29 + qJD(5) * t10 + qJD(6) * t13, t18 * qJD(4) + t1 * qJD(5) + t6 * qJD(6) + ((-pkin(3) * t452 + qJ(4) * t687) * t678 + (t225 * t359 + t226 * t360 + t422 * t447) * t676 + (t185 * t90 + t345 * t378 + t683 * t89) * t674) * t680 + t508 + (t163 * t645 + t352 * t633 + t350 * t634 + t389 * t627 + t390 * t628 + t287 * t630 - t225 * t574 - t226 * t575 + t90 * t582 + pkin(8) * t569 + (t459 * t418 - t457 * t420) * qJ(4) + (Ifges(5,5) * t457 + Ifges(6,5) * t491 + Ifges(7,5) * t336 + Ifges(5,6) * t459 - Ifges(6,6) * t490 + Ifges(7,6) * t512) * t625 + t161 * t703 + t683 * t254 - t490 * t285 / 0.2e1 - Ifges(4,6) * t462 + t447 * t304 + t422 * t348 - pkin(3) * t404 + t378 * t177 + t360 * t362 + t359 * t364 + t345 * t211 + t185 * t252 + t687 * mrSges(5,3) + (Ifges(4,5) + (Ifges(5,1) * t457 + t611) * t627 + (Ifges(5,2) * t459 + t612) * t629 + t567 * pkin(8)) * t465 - t89 * t579 + t217 * t654 + t214 * t658) * qJD(3), qJD(3) * t18 + qJD(5) * t132 + qJD(6) * t84 + t505, t1 * qJD(3) + t132 * qJD(4) + (-t223 * mrSges(6,1) - t222 * mrSges(6,2) - t614 + t615 + t688) * qJD(5) + t21 * qJD(6) + (m(7) * (t460 * t92 + t464 * t91) + (t291 * t460 - t464 * t513) * mrSges(7,3)) * t601 + t507, t6 * qJD(3) + t84 * qJD(4) + t21 * qJD(5) + (t538 - t618 - t619) * qJD(6) + t506; -qJD(2) * t3 - qJD(4) * t42 + qJD(5) * t11 + qJD(6) * t19 - t558, -qJD(4) * t17 + qJD(5) * t2 + qJD(6) * t7 - t508, qJD(4) * t46 + qJD(5) * t16 + qJD(6) * t27, qJD(5) * t157 + qJD(6) * t107 - t485, t157 * qJD(4) + (-t360 * mrSges(6,1) - t359 * mrSges(6,2) + t689 + t736) * qJD(5) + t741 + (m(7) * (-t185 * t464 + t460 * t683) + (-t336 * t460 - t464 * t512) * mrSges(7,3)) * t601 + t487, t107 * qJD(4) + t38 * qJD(5) + t486 + t741; qJD(2) * t40 + qJD(3) * t42, qJD(3) * t17 + qJD(5) * t57 + qJD(6) * t83 - t505, qJD(5) * t67 + qJD(6) * t106 + t485, 0, t504, t500; -qJD(2) * t8 - qJD(3) * t11 + qJD(6) * t34, -qJD(3) * t2 - qJD(4) * t57 + qJD(6) * t22 - t507, -qJD(4) * t67 - t487, -t504, -t414, -t414 - t483; -t14 * qJD(2) - t19 * qJD(3) - t34 * qJD(5), -qJD(3) * t7 - qJD(4) * t83 - qJD(5) * t22 - t506, -qJD(4) * t106 - t486, -t500, t483, 0;];
Cq  = t23;
