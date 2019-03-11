% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR12_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR12_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:35:17
% EndTime: 2019-03-09 14:36:28
% DurationCPUTime: 44.40s
% Computational Cost: add. (21985->986), mult. (50898->1294), div. (0->0), fcn. (39131->14), ass. (0->466)
t703 = mrSges(6,2) - mrSges(7,3);
t641 = m(7) * pkin(11);
t461 = -t641 + t703;
t386 = sin(qJ(2));
t381 = sin(pkin(6));
t513 = qJD(1) * t381;
t351 = t386 * t513;
t342 = pkin(2) * t351;
t391 = cos(qJ(2));
t434 = pkin(9) * t386 - qJ(3) * t391;
t238 = t434 * t513 + t342;
t475 = t391 * t513;
t345 = pkin(8) * t475;
t382 = cos(pkin(6));
t512 = qJD(1) * t382;
t497 = pkin(1) * t512;
t291 = t386 * t497 + t345;
t242 = pkin(3) * t475 + t291;
t385 = sin(qJ(4));
t390 = cos(qJ(4));
t148 = -t238 * t385 + t390 * t242;
t507 = qJD(4) * t385;
t530 = t385 * t386;
t626 = pkin(2) + pkin(9);
t586 = pkin(10) + t626;
t715 = -(pkin(4) * t391 - pkin(10) * t530) * t513 - t148 + t586 * t507;
t149 = t390 * t238 + t385 * t242;
t320 = t586 * t390;
t454 = t390 * t351;
t714 = pkin(10) * t454 + qJD(4) * t320 + t149;
t681 = m(7) + m(6);
t682 = m(4) + m(5);
t706 = -t682 - t681;
t653 = qJ(3) * t706 + mrSges(3,2) - mrSges(4,3);
t380 = qJ(4) + qJ(5);
t376 = sin(t380);
t377 = cos(t380);
t445 = t385 * mrSges(5,1) + t390 * mrSges(5,2);
t642 = m(7) * pkin(5);
t383 = sin(qJ(6));
t388 = cos(qJ(6));
t443 = -mrSges(7,1) * t388 + mrSges(7,2) * t383;
t702 = mrSges(6,1) - t443;
t708 = -t642 - t702;
t713 = t376 * t708 - t461 * t377 - t445;
t499 = qJD(4) + qJD(5);
t317 = t351 + t499;
t604 = -t317 / 0.2e1;
t360 = qJD(2) + t512;
t252 = -t360 * t385 - t390 * t475;
t384 = sin(qJ(5));
t389 = cos(qJ(5));
t407 = -t360 * t390 + t385 * t475;
t425 = t252 * t384 - t389 * t407;
t610 = -t425 / 0.2e1;
t712 = Ifges(6,4) * t610 + Ifges(6,6) * t604;
t463 = t389 * t252 + t384 * t407;
t165 = qJD(6) - t463;
t613 = t165 / 0.2e1;
t139 = t317 * t383 + t388 * t425;
t617 = t139 / 0.2e1;
t138 = t317 * t388 - t383 * t425;
t619 = t138 / 0.2e1;
t711 = Ifges(7,5) * t617 + Ifges(7,6) * t619 + Ifges(7,3) * t613;
t612 = -t463 / 0.2e1;
t710 = Ifges(6,2) * t612;
t503 = qJD(5) * t389;
t504 = qJD(5) * t384;
t506 = qJD(4) * t390;
t215 = -t384 * t506 - t385 * t503 - t389 * t507 - t390 * t504;
t314 = t384 * t390 + t385 * t389;
t240 = t314 * t351;
t663 = t215 - t240;
t524 = t389 * t390;
t216 = -t384 * t507 - t385 * t504 + t499 * t524;
t533 = t384 * t385;
t239 = -t351 * t533 + t389 * t454;
t698 = t239 + t216;
t347 = t391 * t497;
t709 = qJD(3) - t347;
t319 = t586 * t385;
t226 = -t319 * t384 + t389 * t320;
t669 = -qJD(5) * t226 + t384 * t715 - t714 * t389;
t372 = pkin(4) * t390 + pkin(3);
t664 = -(-pkin(8) - t372) * t351 + pkin(4) * t506 + t709;
t705 = t711 + t712;
t535 = t381 * t391;
t294 = qJD(2) * t351 - qJDD(1) * t535;
t500 = qJDD(1) * t382;
t352 = qJDD(2) + t500;
t157 = qJD(4) * t252 + t294 * t385 + t352 * t390;
t158 = qJD(4) * t407 + t294 * t390 - t352 * t385;
t80 = qJD(5) * t463 + t157 * t389 + t158 * t384;
t629 = t80 / 0.2e1;
t81 = -qJD(5) * t425 - t157 * t384 + t158 * t389;
t628 = t81 / 0.2e1;
t509 = qJD(2) * t391;
t295 = (qJD(1) * t509 + qJDD(1) * t386) * t381;
t278 = qJDD(4) + t295;
t259 = qJDD(5) + t278;
t606 = t259 / 0.2e1;
t704 = -mrSges(3,1) + mrSges(4,2);
t701 = -m(5) * pkin(9) - mrSges(5,3) - mrSges(6,3);
t700 = -pkin(11) * t475 + t669;
t699 = pkin(5) * t698 - pkin(11) * t663 + t664;
t333 = t360 * qJ(3);
t189 = t333 + t242;
t153 = -pkin(4) * t252 + t189;
t625 = pkin(3) + pkin(8);
t174 = t351 * t625 - t360 * t626 + t709;
t466 = -qJ(3) * t386 - pkin(1);
t204 = (-t391 * t626 + t466) * t513;
t117 = t390 * t174 - t204 * t385;
t106 = pkin(10) * t407 + t117;
t325 = t351 + qJD(4);
t100 = pkin(4) * t325 + t106;
t118 = t174 * t385 + t204 * t390;
t107 = pkin(10) * t252 + t118;
t525 = t389 * t107;
t55 = t100 * t384 + t525;
t697 = -t153 * mrSges(6,1) + t55 * mrSges(6,3);
t603 = t317 / 0.2e1;
t609 = t425 / 0.2e1;
t611 = t463 / 0.2e1;
t51 = pkin(11) * t317 + t55;
t86 = -pkin(5) * t463 - pkin(11) * t425 + t153;
t25 = -t383 * t51 + t388 * t86;
t26 = t383 * t86 + t388 * t51;
t649 = t25 * mrSges(7,1) - t26 * mrSges(7,2);
t691 = t705 + t710;
t696 = -Ifges(6,4) * t609 - Ifges(6,2) * t611 - Ifges(6,6) * t603 + t649 + t691 + t711;
t45 = qJD(6) * t138 + t259 * t383 + t388 * t80;
t46 = -qJD(6) * t139 + t259 * t388 - t383 * t80;
t79 = qJDD(6) - t81;
t14 = Ifges(7,5) * t45 + Ifges(7,6) * t46 + Ifges(7,3) * t79;
t630 = t79 / 0.2e1;
t635 = t46 / 0.2e1;
t636 = t45 / 0.2e1;
t601 = pkin(1) * t382;
t496 = qJD(2) * t601;
t455 = qJD(1) * t496;
t487 = pkin(1) * t500;
t182 = -pkin(8) * t294 + t386 * t487 + t391 * t455;
t156 = -t352 * qJ(3) - t360 * qJD(3) - t182;
t130 = -pkin(3) * t294 - t156;
t94 = -pkin(4) * t158 + t130;
t22 = -pkin(5) * t81 - pkin(11) * t80 + t94;
t537 = t381 * t386;
t361 = pkin(8) * t537;
t183 = -qJD(2) * t345 - qJDD(1) * t361 - t386 * t455 + t391 * t487;
t403 = qJDD(3) - t183;
t129 = pkin(3) * t295 - t352 * t626 + t403;
t508 = qJD(3) * t386;
t400 = -qJ(3) * t295 + (-pkin(1) * qJDD(1) - qJD(1) * t508) * t381;
t132 = t294 * t626 + t400;
t53 = -qJD(4) * t118 + t390 * t129 - t132 * t385;
t35 = pkin(4) * t278 - pkin(10) * t157 + t53;
t52 = t385 * t129 + t390 * t132 + t174 * t506 - t204 * t507;
t39 = pkin(10) * t158 + t52;
t10 = t100 * t503 - t107 * t504 + t384 * t35 + t389 * t39;
t7 = pkin(11) * t259 + t10;
t2 = qJD(6) * t25 + t22 * t383 + t388 * t7;
t3 = -qJD(6) * t26 + t22 * t388 - t383 * t7;
t652 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t695 = t652 + mrSges(6,1) * t94 + Ifges(7,5) * t636 + Ifges(7,6) * t635 + Ifges(7,3) * t630 + t14 / 0.2e1 + (-t606 - t259 / 0.2e1) * Ifges(6,6) + (-t628 - t81 / 0.2e1) * Ifges(6,2) + (-t629 - t80 / 0.2e1) * Ifges(6,4);
t546 = t107 * t384;
t54 = t100 * t389 - t546;
t694 = t153 * mrSges(6,2) - t54 * mrSges(6,3);
t227 = -t319 * t389 - t320 * t384;
t670 = -qJD(5) * t227 + t714 * t384 + t389 * t715;
t144 = mrSges(6,1) * t317 - mrSges(6,3) * t425;
t89 = -mrSges(7,1) * t138 + mrSges(7,2) * t139;
t693 = t89 - t144;
t419 = t701 + t704;
t442 = mrSges(7,1) * t383 + mrSges(7,2) * t388;
t647 = -t419 + t442;
t392 = cos(qJ(1));
t521 = t391 * t392;
t387 = sin(qJ(1));
t529 = t386 * t387;
t303 = -t382 * t521 + t529;
t534 = t381 * t392;
t692 = t303 * t390 + t385 * t534;
t527 = t387 * t391;
t528 = t386 * t392;
t305 = t382 * t527 + t528;
t536 = t381 * t387;
t228 = t305 * t390 - t385 * t536;
t126 = mrSges(5,1) * t278 - mrSges(5,3) * t157;
t127 = -mrSges(5,2) * t278 + mrSges(5,3) * t158;
t190 = -mrSges(5,2) * t325 + mrSges(5,3) * t252;
t191 = mrSges(5,1) * t325 + mrSges(5,3) * t407;
t426 = t390 * t190 - t385 * t191;
t689 = t426 * qJD(4) + t390 * t126 + t385 * t127;
t688 = -t691 + t697;
t50 = -pkin(5) * t317 - t54;
t416 = t50 * t442;
t576 = Ifges(7,4) * t139;
t65 = Ifges(7,2) * t138 + Ifges(7,6) * t165 + t576;
t552 = t383 * t65;
t135 = Ifges(7,4) * t138;
t66 = Ifges(7,1) * t139 + Ifges(7,5) * t165 + t135;
t632 = -t66 / 0.2e1;
t687 = t388 * t632 - t416 + t552 / 0.2e1 - t694;
t164 = Ifges(6,4) * t463;
t557 = t317 * Ifges(6,5);
t103 = Ifges(6,1) * t425 + t164 + t557;
t686 = Ifges(6,1) * t609 + Ifges(6,4) * t611 + Ifges(6,5) * t603 + t103 / 0.2e1;
t685 = t94 * mrSges(6,2) + 0.2e1 * Ifges(6,1) * t629 + 0.2e1 * Ifges(6,4) * t628 + 0.2e1 * Ifges(6,5) * t606;
t616 = t157 / 0.2e1;
t615 = t158 / 0.2e1;
t605 = t278 / 0.2e1;
t680 = -Ifges(4,4) + Ifges(3,5);
t679 = Ifges(4,5) - Ifges(3,6);
t17 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t62 = mrSges(6,1) * t259 - mrSges(6,3) * t80;
t678 = t17 - t62;
t313 = -t524 + t533;
t378 = t385 * pkin(4);
t368 = qJ(3) + t378;
t195 = pkin(5) * t314 + pkin(11) * t313 + t368;
t136 = t195 * t388 - t227 * t383;
t677 = qJD(6) * t136 + t383 * t699 + t388 * t700;
t137 = t195 * t383 + t227 * t388;
t676 = -qJD(6) * t137 - t383 * t700 + t388 * t699;
t672 = pkin(5) * t475 - t670;
t600 = pkin(1) * t391;
t476 = -pkin(2) - t600;
t213 = pkin(3) * t537 + t361 + (-pkin(9) + t476) * t382;
t515 = pkin(2) * t535 + qJ(3) * t537;
t232 = (-pkin(9) * t391 - pkin(1)) * t381 - t515;
t141 = t390 * t213 - t232 * t385;
t417 = -t382 * t390 + t385 * t535;
t498 = pkin(4) * t537;
t115 = pkin(10) * t417 + t141 + t498;
t142 = t385 * t213 + t390 * t232;
t301 = -t382 * t385 - t390 * t535;
t120 = pkin(10) * t301 + t142;
t668 = t384 * t115 + t389 * t120;
t173 = -mrSges(5,1) * t252 - mrSges(5,2) * t407;
t458 = mrSges(4,1) * t475;
t286 = -mrSges(4,3) * t360 - t458;
t667 = t173 - t286;
t178 = -t240 * t383 + t388 * t475;
t501 = qJD(6) * t388;
t415 = -t215 * t383 + t313 * t501;
t666 = t178 - t415;
t179 = t240 * t388 + t383 * t475;
t502 = qJD(6) * t383;
t526 = t388 * t215;
t414 = t313 * t502 + t526;
t665 = t179 - t414;
t457 = mrSges(3,3) * t351;
t459 = mrSges(4,1) * t351;
t662 = t360 * t704 + t457 + t459;
t181 = t301 * t384 - t389 * t417;
t335 = t388 * t537;
t162 = -t181 * t383 + t335;
t265 = -t376 * t382 - t377 * t535;
t266 = -t376 * t535 + t377 * t382;
t661 = -t265 * t702 + t266 * t703;
t207 = t303 * t377 + t376 * t534;
t209 = -t303 * t376 + t377 * t534;
t660 = -t207 * t702 - t209 * t703;
t205 = -t305 * t377 + t376 * t536;
t206 = t305 * t376 + t377 * t536;
t659 = t205 * t702 + t206 * t703;
t23 = mrSges(7,1) * t79 - mrSges(7,3) * t45;
t24 = -mrSges(7,2) * t79 + mrSges(7,3) * t46;
t657 = -t383 * t23 + t388 * t24;
t656 = -t385 * t52 - t390 * t53;
t143 = -mrSges(6,2) * t317 + mrSges(6,3) * t463;
t97 = -mrSges(7,2) * t165 + mrSges(7,3) * t138;
t98 = mrSges(7,1) * t165 - mrSges(7,3) * t139;
t654 = -t383 * t98 + t388 * t97 + t143;
t11 = -qJD(5) * t55 + t35 * t389 - t384 * t39;
t474 = qJD(2) * t537;
t224 = qJD(4) * t301 + t385 * t474;
t473 = t381 * t509;
t344 = pkin(2) * t474;
t194 = t344 + (qJD(2) * t434 - t508) * t381;
t367 = t386 * t601;
t243 = (t535 * t625 + t367) * qJD(2);
t93 = -qJD(4) * t142 - t194 * t385 + t390 * t243;
t74 = pkin(4) * t473 - pkin(10) * t224 + t93;
t225 = qJD(4) * t417 + t390 * t474;
t92 = t390 * t194 + t213 * t506 - t232 * t507 + t385 * t243;
t85 = pkin(10) * t225 + t92;
t21 = -qJD(5) * t668 - t384 * t85 + t389 * t74;
t651 = t11 * mrSges(6,1) - t10 * mrSges(6,2) + Ifges(6,5) * t80 + Ifges(6,6) * t81 + Ifges(6,3) * t259;
t650 = t54 * mrSges(6,1) - t55 * mrSges(6,2);
t111 = pkin(5) * t425 - pkin(11) * t463;
t648 = t53 * mrSges(5,1) - t52 * mrSges(5,2) + Ifges(5,5) * t157 + Ifges(5,6) * t158 + Ifges(5,3) * t278;
t646 = t653 + t713;
t622 = -t103 / 0.2e1;
t645 = Ifges(6,1) * t610 + Ifges(6,4) * t612 + Ifges(6,5) * t604 + t622;
t644 = -t10 * t314 + t11 * t313 - t54 * t663 - t55 * t698;
t614 = -t165 / 0.2e1;
t618 = -t139 / 0.2e1;
t620 = -t138 / 0.2e1;
t643 = -Ifges(7,5) * t618 - Ifges(7,6) * t620 - Ifges(7,3) * t614 + t649 + t710 + t712;
t15 = t45 * Ifges(7,4) + t46 * Ifges(7,2) + t79 * Ifges(7,6);
t639 = t15 / 0.2e1;
t16 = t45 * Ifges(7,1) + t46 * Ifges(7,4) + t79 * Ifges(7,5);
t638 = t16 / 0.2e1;
t631 = t66 / 0.2e1;
t627 = Ifges(5,1) * t616 + Ifges(5,4) * t615 + Ifges(5,5) * t605;
t607 = -t407 / 0.2e1;
t602 = pkin(1) * t381;
t599 = pkin(4) * t407;
t405 = pkin(4) * t301;
t598 = pkin(4) * t384;
t597 = pkin(4) * t389;
t8 = -pkin(5) * t259 - t11;
t595 = t313 * t8;
t587 = -Ifges(4,6) - Ifges(3,4);
t584 = mrSges(5,2) * t385;
t582 = mrSges(7,3) * t383;
t581 = mrSges(7,3) * t388;
t580 = Ifges(3,4) * t386;
t579 = Ifges(5,4) * t385;
t578 = Ifges(5,4) * t390;
t575 = Ifges(7,4) * t383;
t574 = Ifges(7,4) * t388;
t573 = Ifges(4,6) * t386;
t572 = Ifges(4,6) * t391;
t569 = t117 * mrSges(5,3);
t568 = t118 * mrSges(5,3);
t564 = t463 * Ifges(6,6);
t563 = t425 * Ifges(6,5);
t560 = t252 * Ifges(5,6);
t559 = t407 * Ifges(5,4);
t558 = t407 * Ifges(5,5);
t555 = t317 * Ifges(6,3);
t554 = t325 * Ifges(5,3);
t543 = t303 * t385;
t541 = t305 * t385;
t539 = t313 * t383;
t538 = t313 * t388;
t516 = t692 * pkin(4);
t309 = pkin(8) * t535 + t367;
t514 = t392 * pkin(1) + pkin(8) * t536;
t510 = qJD(1) ^ 2 * t381 ^ 2;
t493 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t492 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t491 = Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1;
t486 = t383 * t537;
t483 = -t552 / 0.2e1;
t306 = -t382 * t529 + t521;
t478 = t306 * pkin(2) + t514;
t263 = -t382 * qJ(3) - t309;
t467 = t501 / 0.2e1;
t465 = -pkin(1) * t387 + pkin(8) * t534;
t464 = -t205 * pkin(5) + pkin(11) * t206;
t223 = t295 * mrSges(4,1) + t352 * mrSges(4,2);
t462 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(3,3);
t460 = t625 * t537;
t456 = mrSges(3,3) * t475;
t231 = pkin(3) * t535 - t263;
t304 = t382 * t528 + t527;
t449 = t304 * pkin(2) - t465;
t447 = mrSges(5,1) * t301 + mrSges(5,2) * t417;
t441 = mrSges(4,2) * t391 - mrSges(4,3) * t386;
t440 = Ifges(5,1) * t385 + t578;
t439 = Ifges(7,1) * t388 - t575;
t438 = Ifges(5,2) * t390 + t579;
t437 = -Ifges(7,2) * t383 + t574;
t436 = Ifges(5,5) * t385 + Ifges(5,6) * t390;
t435 = Ifges(7,5) * t388 - Ifges(7,6) * t383;
t433 = t25 * t388 + t26 * t383;
t432 = -t25 * t383 + t26 * t388;
t68 = pkin(11) * t537 + t668;
t175 = t231 - t405;
t424 = t389 * t301 + t384 * t417;
t96 = -pkin(5) * t424 - pkin(11) * t181 + t175;
t37 = t383 * t96 + t388 * t68;
t36 = -t383 * t68 + t388 * t96;
t70 = t115 * t389 - t120 * t384;
t428 = t117 * t385 - t118 * t390;
t423 = t228 * pkin(4);
t290 = pkin(8) * t351 - t347;
t348 = t391 * t496;
t292 = -pkin(8) * t474 + t348;
t393 = -pkin(10) - pkin(9);
t418 = pkin(4) * t541 - t306 * t393 + t372 * t536 + t478;
t163 = t181 * t388 + t486;
t413 = t138 * t437;
t412 = t139 * t439;
t411 = t165 * t435;
t20 = t115 * t503 - t120 * t504 + t384 * t74 + t389 * t85;
t374 = t382 * qJD(3);
t212 = -qJD(2) * t460 + t348 + t374;
t402 = -qJD(6) * t433 - t3 * t383;
t171 = -pkin(2) * t352 + t403;
t401 = t183 * mrSges(3,1) - t182 * mrSges(3,2) + t171 * mrSges(4,2) - t156 * mrSges(4,3);
t147 = -pkin(4) * t225 + t212;
t399 = t2 * t388 + t402;
t398 = (-t383 * t97 - t388 * t98) * qJD(6) + t657;
t397 = -qJD(4) * t428 - t656;
t395 = t2 * t581 + t383 * t638 + t388 * t639 + t8 * t443 + (Ifges(7,1) * t383 + t574) * t636 + (Ifges(7,2) * t388 + t575) * t635 + t66 * t467 + (Ifges(7,5) * t383 + Ifges(7,6) * t388) * t630 + t651 + (t416 + t483) * qJD(6) + (t413 + t412 + t411) * qJD(6) / 0.2e1;
t371 = -pkin(5) - t597;
t341 = Ifges(3,4) * t475;
t336 = t390 * t534;
t332 = Ifges(4,1) * t352;
t331 = Ifges(3,3) * t352;
t308 = t382 * t600 - t361;
t307 = (-mrSges(3,1) * t391 + mrSges(3,2) * t386) * t381;
t298 = t305 * pkin(2);
t296 = t303 * pkin(2);
t293 = t309 * qJD(2);
t289 = -qJ(3) * t475 + t342;
t288 = t441 * t513;
t285 = -mrSges(3,2) * t360 + t456;
t277 = Ifges(4,4) * t295;
t276 = Ifges(3,5) * t295;
t275 = Ifges(4,5) * t294;
t274 = Ifges(3,6) * t294;
t273 = t382 * t476 + t361;
t264 = -t515 - t602;
t257 = t265 * pkin(5);
t250 = -t292 - t374;
t249 = Ifges(5,4) * t252;
t247 = (-pkin(2) * t391 + t466) * t513;
t244 = t344 + (-qJ(3) * t509 - t508) * t381;
t241 = -qJD(1) * t460 + t347;
t237 = -t333 - t291;
t236 = t360 * Ifges(4,4) + (-Ifges(4,2) * t386 - t572) * t513;
t235 = t360 * Ifges(4,5) + (-t391 * Ifges(4,3) - t573) * t513;
t234 = Ifges(3,1) * t351 + t360 * Ifges(3,5) + t341;
t233 = t360 * Ifges(3,6) + (t391 * Ifges(3,2) + t580) * t513;
t230 = -pkin(2) * t360 + qJD(3) + t290;
t229 = t390 * t536 + t541;
t222 = mrSges(4,1) * t294 - mrSges(4,3) * t352;
t203 = t207 * pkin(5);
t177 = -t239 * t388 + t360 * t383;
t176 = t239 * t383 + t360 * t388;
t161 = pkin(2) * t294 + t400;
t160 = t206 * t388 + t306 * t383;
t159 = -t206 * t383 + t306 * t388;
t152 = -Ifges(5,1) * t407 + t325 * Ifges(5,5) + t249;
t151 = t252 * Ifges(5,2) + t325 * Ifges(5,6) - t559;
t150 = t554 - t558 + t560;
t110 = -mrSges(6,1) * t463 + mrSges(6,2) * t425;
t109 = qJD(5) * t181 + t224 * t384 - t389 * t225;
t108 = qJD(5) * t424 + t224 * t389 + t225 * t384;
t104 = -mrSges(5,1) * t158 + mrSges(5,2) * t157;
t101 = t555 + t563 + t564;
t95 = t111 - t599;
t90 = t157 * Ifges(5,4) + t158 * Ifges(5,2) + t278 * Ifges(5,6);
t76 = -qJD(6) * t163 - t108 * t383 + t388 * t473;
t75 = qJD(6) * t162 + t108 * t388 + t383 * t473;
t67 = -pkin(5) * t537 - t70;
t63 = -mrSges(6,2) * t259 + mrSges(6,3) * t81;
t57 = t106 * t389 - t546;
t56 = t106 * t384 + t525;
t42 = pkin(5) * t109 - pkin(11) * t108 + t147;
t33 = t111 * t383 + t388 * t54;
t32 = t111 * t388 - t383 * t54;
t31 = -mrSges(6,1) * t81 + mrSges(6,2) * t80;
t30 = t383 * t95 + t388 * t57;
t29 = -t383 * t57 + t388 * t95;
t19 = -pkin(5) * t473 - t21;
t18 = pkin(11) * t473 + t20;
t5 = -qJD(6) * t37 - t18 * t383 + t388 * t42;
t4 = qJD(6) * t36 + t18 * t388 + t383 * t42;
t1 = [(-t11 * mrSges(6,3) + t685) * t181 + (-m(3) * t465 - t336 * mrSges(5,1) + mrSges(2,1) * t387 + mrSges(2,2) * t392 + t461 * t207 + (t445 - t653) * t303 + t708 * t209 + t647 * t304 + t682 * t449 + t681 * (pkin(4) * t543 - t304 * t393 - t372 * t534 + t449)) * g(1) + m(6) * (t10 * t668 + t11 * t70 + t147 * t153 + t175 * t94 + t20 * t55 + t21 * t54) + t668 * t63 - t417 * t627 + (-Ifges(5,5) * t417 + Ifges(5,6) * t301) * t605 + (-Ifges(5,4) * t417 + Ifges(5,2) * t301) * t615 + (-Ifges(5,1) * t417 + Ifges(5,4) * t301) * t616 + (-t117 * t224 + t118 * t225 + t301 * t52 + t417 * t53) * mrSges(5,3) + (Ifges(7,5) * t75 + Ifges(7,6) * t76) * t613 + (Ifges(7,5) * t163 + Ifges(7,6) * t162) * t630 + (Ifges(7,1) * t75 + Ifges(7,4) * t76) * t617 + (Ifges(7,1) * t163 + Ifges(7,4) * t162) * t636 + (t162 * t2 - t163 * t3 - t25 * t75 + t26 * t76) * mrSges(7,3) + (Ifges(5,1) * t224 + Ifges(5,4) * t225) * t607 + t75 * t631 + (t276 / 0.2e1 - t274 / 0.2e1 + t331 / 0.2e1 + t332 / 0.2e1 - t277 / 0.2e1 + t275 / 0.2e1 + (Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * t352 + t493 * t295 + t492 * t294 + t401) * t382 + t252 * (Ifges(5,4) * t224 + Ifges(5,2) * t225) / 0.2e1 + m(3) * (t182 * t309 + t183 * t308 + t290 * t293 + t291 * t292) + (Ifges(7,4) * t75 + Ifges(7,2) * t76) * t619 + (Ifges(7,4) * t163 + Ifges(7,2) * t162) * t635 + t175 * t31 + m(5) * (t117 * t93 + t118 * t92 + t130 * t231 + t141 * t53 + t142 * t52 + t189 * t212) + m(7) * (t19 * t50 + t2 * t37 + t25 * t5 + t26 * t4 + t3 * t36 + t67 * t8) + m(4) * (t156 * t263 + t161 * t264 + t171 * t273 + t230 * t293 + t237 * t250 + t244 * t247) + t163 * t638 + t162 * t639 + (t686 + t694) * t108 + (mrSges(6,3) * t10 - t695) * t424 + (t696 - t697) * t109 - t130 * t447 + Ifges(2,3) * qJDD(1) + t8 * (-mrSges(7,1) * t162 + mrSges(7,2) * t163) + t20 * t143 + t21 * t144 + t147 * t110 + t141 * t126 + t142 * t127 + t4 * t97 + t5 * t98 + t19 * t89 + t50 * (-mrSges(7,1) * t76 + mrSges(7,2) * t75) + t76 * t65 / 0.2e1 + t70 * t62 + t67 * t17 + t37 * t24 + t36 * t23 + (-m(3) * t514 - t229 * mrSges(5,1) - t228 * mrSges(5,2) - mrSges(2,1) * t392 + mrSges(2,2) * t387 - m(6) * t418 - t206 * mrSges(6,1) - m(7) * (pkin(5) * t206 + t418) - t160 * mrSges(7,1) - t159 * mrSges(7,2) + t653 * t305 + t461 * t205 + t419 * t306 - t682 * t478) * g(2) + (t462 * t387 * g(2) + (-mrSges(3,1) * t294 - mrSges(3,2) * t295 + (m(3) * t602 - t307) * qJDD(1)) * pkin(1) + (t462 + t584) * g(1) * t392 + (-t156 * mrSges(4,1) + t161 * mrSges(4,2) + t182 * mrSges(3,3) - t679 * t352 - t587 * t295 + (-Ifges(3,2) - Ifges(4,3)) * t294) * t391 + (-t183 * mrSges(3,3) + t171 * mrSges(4,1) - t161 * mrSges(4,3) + t680 * t352 + (Ifges(4,2) + Ifges(3,1)) * t295 + t587 * t294 + t648 + t651) * t386 + ((t237 * mrSges(4,1) - t291 * mrSges(3,3) - t233 / 0.2e1 + t235 / 0.2e1 - t247 * mrSges(4,2) + t492 * t360 + (-pkin(1) * mrSges(3,1) - t386 * t491) * t513) * t386 + (t230 * mrSges(4,1) + t290 * mrSges(3,3) + t564 / 0.2e1 + t563 / 0.2e1 + t117 * mrSges(5,1) - t118 * mrSges(5,2) + t560 / 0.2e1 - t558 / 0.2e1 + t555 / 0.2e1 + t554 / 0.2e1 + t234 / 0.2e1 - t236 / 0.2e1 + t101 / 0.2e1 + t150 / 0.2e1 - t247 * mrSges(4,3) + t493 * t360 + ((-pkin(1) * mrSges(3,2) + t391 * t491) * t381 + (-Ifges(4,3) / 0.2e1 + Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1 + Ifges(4,2) / 0.2e1) * t537) * qJD(1) + t650) * t391) * qJD(2)) * t381 + t308 * (mrSges(3,1) * t352 - mrSges(3,3) * t295) + t309 * (-mrSges(3,2) * t352 - mrSges(3,3) * t294) + t92 * t190 + t93 * t191 + t212 * t173 + t224 * t152 / 0.2e1 + t225 * t151 / 0.2e1 + t189 * (-mrSges(5,1) * t225 + mrSges(5,2) * t224) + t231 * t104 + t263 * t222 + t273 * t223 + t250 * t286 + t244 * t288 + t292 * t285 + t264 * (-mrSges(4,2) * t294 - mrSges(4,3) * t295) + t301 * t90 / 0.2e1 + t662 * t293 + t325 * (Ifges(5,5) * t224 + Ifges(5,6) * t225) / 0.2e1; (-t117 * (mrSges(5,1) * t391 - mrSges(5,3) * t530) - t118 * (mrSges(5,3) * t386 * t390 - mrSges(5,2) * t391) - t247 * (-mrSges(4,2) * t386 - mrSges(4,3) * t391)) * t513 + (t643 + t691) * t239 + (Ifges(7,5) * t414 + Ifges(7,6) * t415) * t613 + (Ifges(7,5) * t179 + Ifges(7,6) * t178) * t614 + (t483 + t686) * t215 + (Ifges(7,4) * t414 + Ifges(7,2) * t415) * t619 + (Ifges(7,4) * t179 + Ifges(7,2) * t178) * t620 + (Ifges(7,1) * t414 + Ifges(7,4) * t415) * t617 + (Ifges(7,1) * t179 + Ifges(7,4) * t178) * t618 + (-pkin(2) * t171 - qJ(3) * t156 - qJD(3) * t237 - t247 * t289) * m(4) + (-t435 * t630 - t437 * t635 - t439 * t636 + t65 * t467 - t685) * t313 + (qJD(6) * t66 + t15) * t539 / 0.2e1 + (t391 * (Ifges(4,3) * t386 - t572) + t386 * (-Ifges(4,2) * t391 + t573)) * t510 / 0.2e1 + (t386 * t233 + t391 * t236) * t513 / 0.2e1 + (-t568 - t151 / 0.2e1) * t506 + (-m(4) * t237 + t285 - t286 - t456) * t290 + (t130 * qJ(3) - t117 * t148 - t118 * t149 + (qJD(3) - t241) * t189) * m(5) + t325 * (mrSges(5,1) * t390 - t584) * t189 + (-t386 * (Ifges(3,1) * t391 - t580) / 0.2e1 + pkin(1) * (mrSges(3,1) * t386 + mrSges(3,2) * t391)) * t510 + t645 * t240 + t644 * mrSges(6,3) - ((t390 * t151 + t385 * t152 + t235) * t386 + (-Ifges(3,2) * t351 + t101 + t150 + t234 + t341) * t391 + t325 * (Ifges(5,3) * t391 + t386 * t436) - t407 * (Ifges(5,5) * t391 + t386 * t440) + t252 * (Ifges(5,6) * t391 + t386 * t438) + (t386 * t679 + t391 * t680) * t360) * t513 / 0.2e1 - (t252 * t438 + t325 * t436 - t407 * t440) * qJD(4) / 0.2e1 + (t569 - t152 / 0.2e1) * t507 + (Ifges(5,5) * t390 - Ifges(5,6) * t385) * t605 - t230 * t458 + t401 + (-Ifges(5,2) * t385 + t578) * t615 + (Ifges(5,1) * t390 - t579) * t616 + t390 * t627 + t526 * t631 + t179 * t632 + (t104 - t222) * qJ(3) - t442 * t595 - t16 * t538 / 0.2e1 - t277 - t274 + t275 + t276 + (t307 - t682 * t515 - t681 * (t385 * t498 + t515) + (t441 + (t393 * t681 - t442 + t701) * t391 + t713 * t386) * t381) * g(3) - t178 * t65 / 0.2e1 + (-m(5) * t397 - t689) * t626 + t332 + t695 * t314 + t696 * t216 + (mrSges(6,1) * t698 + t663 * mrSges(6,2)) * t153 + t331 + t130 * t445 + t136 * t23 + t137 * t24 + (-t681 * (t305 * t393 + t306 * t378 - t298) + t682 * t298 + t646 * t306 + t647 * t305) * g(1) + (-t681 * (t303 * t393 + t304 * t378 - t296) + t682 * t296 + t646 * t304 + t647 * t303) * g(2) - t385 * t90 / 0.2e1 + t678 * t226 + t676 * t98 + t677 * t97 + (t136 * t3 + t137 * t2 + t226 * t8 + t25 * t676 + t26 * t677 + t50 * t672) * m(7) + t368 * t31 - t149 * t190 - t148 * t191 + (Ifges(6,5) * t610 + Ifges(6,6) * t612 + Ifges(6,3) * t604 - t650) * t475 - pkin(2) * t223 + t227 * t63 - t241 * t173 - t289 * t288 + t656 * mrSges(5,3) + (-m(4) * t230 + t457 - t662) * t291 + t664 * t110 + (mrSges(7,1) * t666 - mrSges(7,2) * t665) * t50 + (t2 * t539 + t25 * t665 - t26 * t666 + t3 * t538) * mrSges(7,3) + t667 * qJD(3) + t669 * t143 + t670 * t144 + (t10 * t227 - t11 * t226 + t153 * t664 + t368 * t94 + t54 * t670 + t55 * t669) * m(6) + t672 * t89 - t237 * t459; t223 + (t288 + t426) * t351 - t176 * t98 - t177 * t97 + (t63 + t398) * t314 + (-t110 - t667) * t360 + t678 * t313 + t654 * t216 + t239 * t143 - t693 * t663 - (-g(1) * t305 - g(2) * t303 + g(3) * t535) * t706 + (-t176 * t25 - t177 * t26 + t216 * t432 + t314 * t399 - t50 * t663 + t595) * m(7) + (-t153 * t360 - t644) * m(6) + (-t189 * t360 - t351 * t428 + t397) * m(5) + (t237 * t360 + t247 * t351 + t171) * m(4) + t689; (t25 * t581 + t26 * t582 + t435 * t614 + t437 * t620 + t439 * t618 + t645 + t687) * t463 + (t371 * t8 + (t384 * t50 + t389 * t432) * qJD(5) * pkin(4) - t25 * t29 - t26 * t30 - t50 * t56) * m(7) + t654 * pkin(4) * t503 + t407 * (Ifges(5,1) * t252 + t559) / 0.2e1 - t407 * t568 - (Ifges(5,2) * t407 + t152 + t249) * t252 / 0.2e1 - t189 * (-mrSges(5,1) * t407 + mrSges(5,2) * t252) - t325 * (Ifges(5,5) * t252 + Ifges(5,6) * t407) / 0.2e1 + ((t10 * t384 + t11 * t389 + (-t384 * t54 + t389 * t55) * qJD(5)) * pkin(4) + t153 * t599 + t54 * t56 - t55 * t57) * m(6) + t62 * t597 + t63 * t598 + t151 * t607 + (-m(7) * (-pkin(11) * t209 + t203 + t516) - m(6) * t516 - t692 * mrSges(5,1) - (t336 - t543) * mrSges(5,2) + t660) * g(2) + t252 * t569 + (-t25 * t501 - t26 * t502) * mrSges(7,3) - t3 * t582 + (-t643 + t688) * t425 + t110 * t599 + t648 + t693 * (pkin(4) * t504 - t56) - t57 * t143 - t30 * t97 - t29 * t98 + t371 * t17 - t117 * t190 + t118 * t191 + (m(7) * t399 - t501 * t98 - t502 * t97 + t657) * (pkin(11) + t598) + (-m(6) * t423 - m(7) * (t423 + t464) - mrSges(5,1) * t228 + mrSges(5,2) * t229 + t659) * g(1) + (-m(6) * t405 - m(7) * (pkin(11) * t266 + t257 + t405) - t447 + t661) * g(3) + t395; -t8 * t642 - m(7) * (t25 * t32 + t26 * t33 + t50 * t55) + ((g(2) * t209 - g(3) * t266) * m(7) + t398) * pkin(11) + (-t649 + t688 - t705) * t425 + t402 * mrSges(7,3) + (t622 - t164 / 0.2e1 - t411 / 0.2e1 - t413 / 0.2e1 - t412 / 0.2e1 - t557 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t425 + t433 * mrSges(7,3) + t687) * t463 + t399 * t641 + (-m(7) * t464 + t659) * g(1) + (-m(7) * t257 + t661) * g(3) - t54 * t143 - t33 * t97 - t32 * t98 - pkin(5) * t17 - t693 * t55 + (-m(7) * t203 + t660) * g(2) + t395; -t50 * (mrSges(7,1) * t139 + mrSges(7,2) * t138) + (Ifges(7,1) * t138 - t576) * t618 + t65 * t617 + (Ifges(7,5) * t138 - Ifges(7,6) * t139) * t614 - t25 * t97 + t26 * t98 - g(1) * (mrSges(7,1) * t159 - mrSges(7,2) * t160) - g(2) * ((t209 * t383 + t304 * t388) * mrSges(7,1) + (t209 * t388 - t304 * t383) * mrSges(7,2)) - g(3) * ((-t266 * t383 + t335) * mrSges(7,1) + (-t266 * t388 - t486) * mrSges(7,2)) + (t138 * t25 + t139 * t26) * mrSges(7,3) + t14 + (-Ifges(7,2) * t139 + t135 + t66) * t620 + t652;];
tau  = t1;
