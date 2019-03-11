% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:53
% EndTime: 2019-03-09 19:26:25
% DurationCPUTime: 53.67s
% Computational Cost: add. (21625->1041), mult. (51277->1387), div. (0->0), fcn. (40540->12), ass. (0->446)
t414 = cos(qJ(2));
t567 = cos(pkin(6));
t496 = t567 * qJD(1);
t482 = pkin(1) * t496;
t410 = sin(qJ(2));
t406 = sin(pkin(6));
t547 = qJD(1) * t406;
t514 = t410 * t547;
t306 = -pkin(8) * t514 + t414 * t482;
t436 = (pkin(2) * t410 - pkin(9) * t414) * t406;
t307 = qJD(1) * t436;
t409 = sin(qJ(3));
t413 = cos(qJ(3));
t199 = -t409 * t306 + t307 * t413;
t637 = pkin(9) - pkin(10);
t367 = t637 * t413;
t415 = -pkin(3) - pkin(4);
t552 = t413 * t414;
t744 = -(-pkin(10) * t552 + t410 * t415) * t547 + t199 + qJD(3) * t367;
t200 = t413 * t306 + t409 * t307;
t176 = qJ(4) * t514 + t200;
t513 = t414 * t547;
t485 = t409 * t513;
t544 = qJD(3) * t409;
t743 = pkin(10) * t485 + t637 * t544 + t176;
t715 = Ifges(4,1) + Ifges(5,1);
t714 = Ifges(5,4) + Ifges(4,5);
t484 = t413 * t513;
t543 = qJD(3) * t413;
t742 = t409 * qJD(4) + (-t484 + t543) * qJ(4);
t408 = sin(qJ(5));
t412 = cos(qJ(5));
t366 = t637 * t409;
t441 = t412 * t366 - t367 * t408;
t702 = qJD(5) * t441 + t408 * t744 - t743 * t412;
t518 = pkin(1) * t567;
t394 = t410 * t518;
t532 = t415 * t409;
t556 = t406 * t414;
t694 = -(-t394 + (-pkin(8) + t532) * t556) * qJD(1) + qJD(3) * t532 + t742;
t607 = cos(qJ(1));
t470 = t567 * t607;
t606 = sin(qJ(1));
t325 = t410 * t470 + t414 * t606;
t517 = t406 * t607;
t237 = t325 * t409 + t413 * t517;
t238 = t325 * t413 - t409 * t517;
t143 = t237 * t408 + t238 * t412;
t324 = t410 * t606 - t414 * t470;
t407 = sin(qJ(6));
t411 = cos(qJ(6));
t741 = t143 * t407 + t324 * t411;
t740 = -t143 * t411 + t324 * t407;
t717 = -mrSges(5,2) - mrSges(4,3);
t719 = -m(7) - m(6);
t668 = t637 * t719 + mrSges(3,2) + mrSges(6,3) + t717;
t557 = t406 * t410;
t512 = qJD(2) * t557;
t537 = qJDD(1) * t406;
t313 = -qJD(1) * t512 + t414 * t537;
t299 = qJDD(3) - t313;
t615 = t299 / 0.2e1;
t538 = qJD(1) * qJD(2);
t506 = t414 * t538;
t314 = (qJDD(1) * t410 + t506) * t406;
t492 = t567 * qJDD(1);
t384 = t492 + qJDD(2);
t439 = t496 + qJD(2);
t274 = t409 * t514 - t413 * t439;
t545 = qJD(3) * t274;
t168 = t413 * t314 + t409 * t384 - t545;
t630 = t168 / 0.2e1;
t739 = t615 * t714 + t630 * t715;
t275 = t409 * t439 + t413 * t514;
t169 = qJD(3) * t275 + t409 * t314 - t413 * t384;
t628 = t169 / 0.2e1;
t713 = Ifges(5,5) - Ifges(4,4);
t712 = Ifges(5,2) + Ifges(4,3);
t711 = Ifges(4,6) - Ifges(5,6);
t334 = t408 * t409 + t412 * t413;
t223 = (qJD(3) - qJD(5)) * t334;
t541 = qJD(5) * t412;
t542 = qJD(5) * t408;
t224 = t408 * t543 + t409 * t541 - t412 * t544 - t413 * t542;
t251 = t408 * t484 - t412 * t485;
t261 = t334 * t556;
t252 = qJD(1) * t261;
t738 = t694 + (-t223 + t252) * pkin(11) + (t224 - t251) * pkin(5);
t737 = pkin(11) * t514 + t702;
t266 = Ifges(4,4) * t274;
t579 = t274 * Ifges(5,5);
t733 = t513 - qJD(3);
t695 = t715 * t275 - t714 * t733 - t266 + t579;
t550 = pkin(8) * t556 + t394;
t295 = t567 * pkin(9) + t550;
t250 = qJD(2) * pkin(9) + qJD(1) * t295;
t259 = (-pkin(2) * t414 - pkin(9) * t410 - pkin(1)) * t547;
t165 = -t409 * t250 + t413 * t259;
t117 = pkin(10) * t275 + t165;
t734 = -qJD(4) + t117;
t102 = -t415 * t733 - t734;
t345 = t733 * qJ(4);
t166 = t413 * t250 + t409 * t259;
t473 = pkin(10) * t274 + t166;
t108 = -t345 + t473;
t48 = t102 * t408 + t108 * t412;
t185 = t274 * t408 + t275 * t412;
t678 = -t733 - qJD(5);
t135 = -t185 * t407 - t411 * t678;
t136 = t185 * t411 - t407 * t678;
t183 = t274 * t412 - t408 * t275;
t178 = qJD(6) - t183;
t52 = t136 * Ifges(7,5) + t135 * Ifges(7,6) + t178 * Ifges(7,3);
t587 = Ifges(6,4) * t185;
t93 = t183 * Ifges(6,2) - Ifges(6,6) * t678 + t587;
t736 = mrSges(6,3) * t48 - t52 / 0.2e1 + t93 / 0.2e1;
t47 = t102 * t412 - t108 * t408;
t45 = pkin(5) * t678 - t47;
t594 = mrSges(6,3) * t185;
t139 = -mrSges(6,1) * t678 - t594;
t80 = -mrSges(7,1) * t135 + mrSges(7,2) * t136;
t703 = t139 - t80;
t735 = -m(7) * t45 + t703;
t533 = m(7) * pkin(11) + mrSges(7,3);
t704 = mrSges(6,2) - t533;
t361 = -mrSges(7,1) * t411 + mrSges(7,2) * t407;
t429 = m(7) * pkin(5) - t361;
t705 = mrSges(6,1) + t429;
t720 = t237 * t412 - t238 * t408;
t732 = t143 * t704 - t705 * t720;
t469 = t567 * t606;
t327 = -t410 * t469 + t414 * t607;
t516 = t406 * t606;
t241 = t327 * t409 - t413 * t516;
t242 = t327 * t413 + t409 * t516;
t149 = t241 * t408 + t242 * t412;
t443 = -t241 * t412 + t242 * t408;
t731 = t149 * t704 + t705 * t443;
t629 = -t169 / 0.2e1;
t729 = t628 * t713 + t739;
t254 = t366 * t408 + t367 * t412;
t701 = -qJD(5) * t254 + t743 * t408 + t412 * t744;
t134 = -t345 + t166;
t265 = Ifges(5,5) * t275;
t153 = -Ifges(5,6) * t733 + t274 * Ifges(5,3) + t265;
t578 = t275 * Ifges(4,4);
t156 = -t274 * Ifges(4,2) - Ifges(4,6) * t733 + t578;
t580 = t166 * mrSges(4,3);
t726 = -t580 - mrSges(5,2) * t134 - t156 / 0.2e1 + t153 / 0.2e1;
t612 = -t678 / 0.2e1;
t622 = t185 / 0.2e1;
t624 = t183 / 0.2e1;
t626 = t178 / 0.2e1;
t633 = t136 / 0.2e1;
t635 = t135 / 0.2e1;
t249 = -t439 * pkin(2) - t306;
t127 = t274 * pkin(3) - t275 * qJ(4) + t249;
t107 = -pkin(4) * t274 - t127;
t46 = -pkin(11) * t678 + t48;
t49 = -pkin(5) * t183 - pkin(11) * t185 + t107;
t20 = -t407 * t46 + t411 * t49;
t21 = t407 * t49 + t411 * t46;
t665 = mrSges(6,1) * t107 + mrSges(7,1) * t20 - mrSges(7,2) * t21;
t725 = -Ifges(6,4) * t622 + Ifges(7,5) * t633 - Ifges(6,2) * t624 - Ifges(6,6) * t612 + Ifges(7,6) * t635 + Ifges(7,3) * t626 + t665 - t736;
t623 = -t185 / 0.2e1;
t613 = t678 / 0.2e1;
t625 = -t183 / 0.2e1;
t627 = -t178 / 0.2e1;
t634 = -t136 / 0.2e1;
t636 = -t135 / 0.2e1;
t661 = Ifges(7,5) * t634 - Ifges(6,2) * t625 - Ifges(6,6) * t613 + Ifges(7,6) * t636 + Ifges(7,3) * t627 - t665;
t724 = -Ifges(6,4) * t623 + t661 + t736;
t335 = -t408 * t413 + t409 * t412;
t462 = mrSges(5,1) * t413 + mrSges(5,3) * t409;
t464 = mrSges(4,1) * t413 - mrSges(4,2) * t409;
t685 = -t462 - t464;
t722 = t334 * t705 + t335 * t704 + mrSges(3,1) - t685;
t721 = 0.2e1 * Ifges(4,2) * t629 - t168 * Ifges(5,5) / 0.2e1 - t299 * Ifges(5,6) / 0.2e1 + Ifges(4,4) * t630 + (t629 - t628) * Ifges(5,3) + (t711 + Ifges(4,6)) * t615;
t104 = pkin(5) * t185 - pkin(11) * t183;
t287 = qJDD(5) - t299;
t616 = t287 / 0.2e1;
t66 = -qJD(5) * t185 - t168 * t408 + t169 * t412;
t642 = t66 / 0.2e1;
t65 = qJD(5) * t183 + t168 * t412 + t169 * t408;
t643 = t65 / 0.2e1;
t650 = Ifges(6,1) * t643 + Ifges(6,4) * t642 + Ifges(6,5) * t616;
t718 = mrSges(5,1) + mrSges(4,1);
t716 = -mrSges(5,3) + mrSges(4,2);
t34 = qJD(6) * t135 + t287 * t407 + t411 * t65;
t35 = -qJD(6) * t136 + t287 * t411 - t407 * t65;
t13 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t50 = mrSges(6,1) * t287 - mrSges(6,3) * t65;
t710 = t13 - t50;
t400 = t409 * qJ(4);
t353 = -t413 * pkin(3) - pkin(2) - t400;
t332 = t413 * pkin(4) - t353;
t201 = pkin(5) * t334 - pkin(11) * t335 + t332;
t123 = t201 * t411 - t254 * t407;
t709 = qJD(6) * t123 + t407 * t738 + t411 * t737;
t124 = t201 * t407 + t254 * t411;
t708 = -qJD(6) * t124 - t407 * t737 + t411 * t738;
t706 = t107 * mrSges(6,2);
t459 = t407 * mrSges(7,1) + t411 * mrSges(7,2);
t435 = t45 * t459;
t700 = -pkin(5) * t514 - t701;
t591 = Ifges(3,4) * t410;
t419 = Ifges(3,6) * t567 + (Ifges(3,2) * t414 + t591) * t406;
t699 = t185 * Ifges(6,5) + Ifges(3,6) * qJD(2) + t183 * Ifges(6,6) - Ifges(6,3) * t678 + qJD(1) * t419;
t551 = pkin(2) * t556 + pkin(9) * t557;
t604 = pkin(1) * t406;
t296 = -t551 - t604;
t192 = -t409 * t295 + t296 * t413;
t172 = pkin(3) * t556 - t192;
t323 = t409 * t567 + t413 * t557;
t122 = pkin(4) * t556 - pkin(10) * t323 + t172;
t193 = t413 * t295 + t409 * t296;
t171 = -qJ(4) * t556 + t193;
t322 = t409 * t557 - t413 * t567;
t129 = pkin(10) * t322 + t171;
t697 = t408 * t122 + t412 * t129;
t696 = -t274 * t711 + t275 * t714 - t712 * t733;
t204 = -t252 * t407 - t411 * t514;
t539 = qJD(6) * t411;
t434 = t223 * t407 + t335 * t539;
t693 = t204 + t434;
t205 = t252 * t411 - t407 * t514;
t540 = qJD(6) * t407;
t433 = -t223 * t411 + t335 * t540;
t692 = t205 + t433;
t562 = t324 * t413;
t691 = -pkin(3) * t562 - t324 * t400;
t326 = t410 * t607 + t414 * t469;
t561 = t326 * t413;
t690 = -pkin(3) * t561 - t326 * t400;
t546 = qJD(2) * t414;
t511 = t406 * t546;
t235 = qJD(3) * t323 + t409 * t511;
t236 = -qJD(3) * t322 + t413 * t511;
t311 = t550 * qJD(2);
t106 = t235 * pkin(3) - t236 * qJ(4) - t323 * qJD(4) + t311;
t493 = -t322 * pkin(3) + t323 * qJ(4);
t689 = pkin(3) * t544 - (t394 + (pkin(3) * t409 + pkin(8)) * t556) * qJD(1) - t742;
t350 = t412 * qJ(4) + t408 * t415;
t687 = -t409 * t711 + t413 * t714;
t583 = Ifges(5,5) * t409;
t589 = Ifges(4,4) * t409;
t686 = t413 * t715 + t583 - t589;
t683 = t168 * t714 - t169 * t711 + t299 * t712;
t438 = qJD(2) * t482;
t472 = pkin(1) * t492;
t211 = pkin(8) * t313 + t410 * t472 + t414 * t438;
t187 = pkin(9) * t384 + t211;
t531 = pkin(1) * t537;
t194 = -pkin(2) * t313 - pkin(9) * t314 - t531;
t71 = t413 * t187 + t409 * t194 - t250 * t544 + t259 * t543;
t55 = t299 * qJ(4) - qJD(4) * t733 + t71;
t72 = -t409 * t187 + t194 * t413 - t250 * t543 - t259 * t544;
t427 = qJDD(4) - t72;
t57 = -pkin(3) * t299 + t427;
t682 = t409 * t57 + t413 * t55;
t64 = qJDD(6) - t66;
t16 = mrSges(7,1) * t64 - mrSges(7,3) * t34;
t17 = -mrSges(7,2) * t64 + mrSges(7,3) * t35;
t681 = -t407 * t16 + t411 * t17;
t680 = -t409 * t72 + t413 * t71;
t679 = qJD(4) - t165;
t656 = t406 ^ 2;
t677 = (pkin(1) * (mrSges(3,1) * t410 + mrSges(3,2) * t414) - t410 * (Ifges(3,1) * t414 - t591) / 0.2e1) * t656;
t595 = mrSges(6,3) * t183;
t138 = mrSges(6,2) * t678 + t595;
t90 = -mrSges(7,2) * t178 + mrSges(7,3) * t135;
t91 = mrSges(7,1) * t178 - mrSges(7,3) * t136;
t674 = -t407 * t91 + t411 * t90 + t138;
t39 = -pkin(10) * t168 + t299 * t415 + t427;
t41 = pkin(10) * t169 + t55;
t8 = -qJD(5) * t48 + t39 * t412 - t408 * t41;
t308 = qJD(2) * t436;
t387 = pkin(8) * t557;
t330 = t414 * t518 - t387;
t310 = t330 * qJD(2);
t110 = -t295 * t543 - t296 * t544 + t308 * t413 - t409 * t310;
t83 = -pkin(10) * t236 + t415 * t512 - t110;
t109 = -t295 * t544 + t296 * t543 + t409 * t308 + t413 * t310;
t98 = qJ(4) * t512 - qJD(4) * t556 + t109;
t84 = pkin(10) * t235 + t98;
t19 = -qJD(5) * t697 - t408 * t84 + t412 * t83;
t212 = -t406 * pkin(8) * t506 - qJDD(1) * t387 - t410 * t438 + t414 * t472;
t188 = -t384 * pkin(2) - t212;
t56 = t169 * pkin(3) - t168 * qJ(4) - t275 * qJD(4) + t188;
t42 = -pkin(4) * t169 - t56;
t12 = -pkin(5) * t66 - pkin(11) * t65 + t42;
t7 = t102 * t541 - t108 * t542 + t408 * t39 + t412 * t41;
t5 = pkin(11) * t287 + t7;
t1 = qJD(6) * t20 + t12 * t407 + t411 * t5;
t2 = -qJD(6) * t21 + t12 * t411 - t407 * t5;
t671 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t491 = mrSges(3,3) * t514;
t670 = -m(4) * t249 + mrSges(3,1) * t439 - mrSges(4,1) * t274 - mrSges(4,2) * t275 - t491;
t418 = t1 * t411 - t2 * t407 + (-t20 * t411 - t21 * t407) * qJD(6);
t667 = m(7) * t418 - t91 * t539 - t90 * t540 + t681;
t662 = t459 + t668;
t450 = Ifges(7,5) * t411 - Ifges(7,6) * t407;
t584 = Ifges(7,4) * t411;
t453 = -Ifges(7,2) * t407 + t584;
t585 = Ifges(7,4) * t407;
t456 = Ifges(7,1) * t411 - t585;
t586 = Ifges(7,4) * t136;
t53 = Ifges(7,2) * t135 + Ifges(7,6) * t178 + t586;
t131 = Ifges(7,4) * t135;
t54 = Ifges(7,1) * t136 + Ifges(7,5) * t178 + t131;
t572 = t411 * t54;
t592 = mrSges(7,3) * t411;
t593 = mrSges(7,3) * t407;
t609 = t407 / 0.2e1;
t660 = Ifges(6,1) * t623 + Ifges(6,5) * t613 + t20 * t592 + t21 * t593 + t450 * t627 + t453 * t636 + t456 * t634 - t435 - t572 / 0.2e1 + t53 * t609 - t706;
t644 = t64 / 0.2e1;
t648 = t35 / 0.2e1;
t649 = t34 / 0.2e1;
t9 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t64;
t659 = -mrSges(6,3) * t7 + Ifges(7,5) * t649 + Ifges(7,6) * t648 + Ifges(7,3) * t644 + t9 / 0.2e1 + t671 + (-t616 - t287 / 0.2e1) * Ifges(6,6) + (-t642 - t66 / 0.2e1) * Ifges(6,2) + (-t643 - t65 / 0.2e1) * Ifges(6,4);
t10 = t34 * Ifges(7,4) + t35 * Ifges(7,2) + t64 * Ifges(7,6);
t653 = t10 / 0.2e1;
t11 = t34 * Ifges(7,1) + t35 * Ifges(7,4) + t64 * Ifges(7,5);
t652 = t11 / 0.2e1;
t175 = Ifges(6,4) * t183;
t94 = t185 * Ifges(6,1) - Ifges(6,5) * t678 + t175;
t639 = -t94 / 0.2e1;
t638 = t94 / 0.2e1;
t620 = -t274 / 0.2e1;
t619 = t274 / 0.2e1;
t617 = t275 / 0.2e1;
t610 = -t733 / 0.2e1;
t603 = pkin(9) * t326;
t315 = t322 * pkin(4);
t600 = t324 * pkin(9);
t599 = t47 * mrSges(6,3);
t597 = mrSges(5,2) * t274;
t596 = mrSges(5,2) * t275;
t590 = Ifges(3,4) * t414;
t588 = Ifges(4,4) * t413;
t582 = Ifges(5,5) * t413;
t581 = t165 * mrSges(4,3);
t560 = t335 * t407;
t559 = t335 * t411;
t558 = t733 * t412;
t553 = t409 * t414;
t189 = t275 * pkin(3) + t274 * qJ(4);
t548 = t607 * pkin(1) + pkin(8) * t516;
t536 = Ifges(6,5) * t65 + Ifges(6,6) * t66 + Ifges(6,3) * t287;
t535 = pkin(9) * t544;
t534 = pkin(9) * t543;
t528 = t406 * t553;
t527 = t406 * t552;
t524 = t572 / 0.2e1;
t520 = Ifges(3,5) * t314 + Ifges(3,6) * t313 + Ifges(3,3) * t384;
t519 = t327 * pkin(2) + t548;
t294 = -t567 * pkin(2) - t330;
t509 = t557 / 0.2e1;
t505 = -t547 / 0.2e1;
t499 = -t540 / 0.2e1;
t317 = t324 * pkin(2);
t498 = pkin(9) * t325 - t317;
t319 = t326 * pkin(2);
t497 = pkin(9) * t327 - t319;
t119 = -t299 * mrSges(5,1) + t168 * mrSges(5,2);
t495 = -t237 * pkin(3) + qJ(4) * t238;
t494 = -t241 * pkin(3) + qJ(4) * t242;
t490 = mrSges(3,3) * t513;
t487 = pkin(3) * t527 + qJ(4) * t528 + t551;
t476 = t414 * t505;
t133 = -pkin(4) * t275 - t189;
t471 = -pkin(1) * t606 + pkin(8) * t517;
t465 = mrSges(4,1) * t322 + mrSges(4,2) * t323;
t463 = t322 * mrSges(5,1) - t323 * mrSges(5,3);
t210 = t322 * t408 + t323 * t412;
t442 = t322 * t412 - t323 * t408;
t461 = -mrSges(6,1) * t442 + mrSges(6,2) * t210;
t173 = -t210 * t407 + t411 * t556;
t174 = t210 * t411 + t407 * t556;
t460 = mrSges(7,1) * t173 - mrSges(7,2) * t174;
t455 = -Ifges(4,2) * t409 + t588;
t451 = Ifges(5,3) * t409 + t582;
t449 = -t20 * t407 + t21 * t411;
t70 = pkin(11) * t556 + t697;
t170 = t294 - t493;
t128 = -t170 - t315;
t75 = -pkin(5) * t442 - pkin(11) * t210 + t128;
t29 = t407 * t75 + t411 * t70;
t28 = -t407 * t70 + t411 * t75;
t349 = -t408 * qJ(4) + t412 * t415;
t73 = t122 * t412 - t129 * t408;
t437 = -t325 * pkin(2) + t471;
t18 = t122 * t541 - t129 * t542 + t408 * t83 + t412 * t84;
t430 = t242 * pkin(3) + qJ(4) * t241 + t519;
t428 = -g(1) * t241 - g(2) * t237 - g(3) * t322;
t424 = t242 * pkin(4) + t430;
t422 = pkin(4) * t527 - pkin(10) * t557 + t487;
t421 = -pkin(3) * t238 - qJ(4) * t237 + t437;
t420 = -pkin(4) * t238 + t421;
t417 = t406 * t439 * (Ifges(3,5) * t414 - Ifges(3,6) * t410);
t89 = -t235 * pkin(4) - t106;
t6 = -pkin(5) * t287 - t8;
t416 = t53 * t499 + t536 + t8 * mrSges(6,1) - t7 * mrSges(6,2) + t411 * t653 + t11 * t609 + t6 * t361 + (Ifges(7,5) * t407 + Ifges(7,6) * t411) * t644 + (Ifges(7,2) * t411 + t585) * t648 + (Ifges(7,1) * t407 + t584) * t649 + t1 * t592 - t2 * t593 + (-t20 * t539 - t21 * t540) * mrSges(7,3) + (t450 * t626 + t453 * t635 + t456 * t633 + t435 + t524) * qJD(6);
t379 = Ifges(3,4) * t513;
t341 = pkin(5) - t349;
t328 = (-mrSges(3,1) * t414 + mrSges(3,2) * t410) * t406;
t309 = t550 * qJD(1);
t304 = -mrSges(3,2) * t439 + t490;
t246 = Ifges(3,1) * t514 + Ifges(3,5) * t439 + t379;
t216 = -mrSges(5,3) * t733 - t597;
t215 = mrSges(5,1) * t733 + t596;
t214 = -mrSges(4,1) * t733 - mrSges(4,3) * t275;
t213 = mrSges(4,2) * t733 - mrSges(4,3) * t274;
t203 = t275 * t407 - t411 * t558;
t202 = t275 * t411 + t407 * t558;
t190 = mrSges(5,1) * t274 - mrSges(5,3) * t275;
t179 = -pkin(3) * t514 - t199;
t132 = pkin(3) * t733 + t679;
t121 = -mrSges(5,2) * t169 + mrSges(5,3) * t299;
t120 = -mrSges(4,2) * t299 - mrSges(4,3) * t169;
t118 = mrSges(4,1) * t299 - mrSges(4,3) * t168;
t112 = t149 * t411 - t326 * t407;
t111 = -t149 * t407 - t326 * t411;
t105 = -pkin(3) * t512 - t110;
t103 = -mrSges(6,1) * t183 + mrSges(6,2) * t185;
t101 = qJD(5) * t442 + t235 * t408 + t236 * t412;
t100 = qJD(5) * t210 - t235 * t412 + t236 * t408;
t96 = mrSges(4,1) * t169 + mrSges(4,2) * t168;
t95 = mrSges(5,1) * t169 - mrSges(5,3) * t168;
t69 = -pkin(5) * t556 - t73;
t68 = t412 * t117 + t408 * t473;
t61 = -t104 + t133;
t60 = qJD(6) * t173 + t101 * t411 - t407 * t512;
t59 = -qJD(6) * t174 - t101 * t407 - t411 * t512;
t51 = -mrSges(6,2) * t287 + mrSges(6,3) * t66;
t31 = t104 * t407 + t411 * t47;
t30 = t104 * t411 - t407 * t47;
t27 = t100 * pkin(5) - t101 * pkin(11) + t89;
t26 = -mrSges(6,1) * t66 + mrSges(6,2) * t65;
t25 = t407 * t61 + t411 * t68;
t24 = -t407 * t68 + t411 * t61;
t15 = pkin(5) * t512 - t19;
t14 = -pkin(11) * t512 + t18;
t4 = -qJD(6) * t29 - t14 * t407 + t27 * t411;
t3 = qJD(6) * t28 + t14 * t411 + t27 * t407;
t22 = [(Ifges(7,5) * t60 + Ifges(7,6) * t59) * t626 + (Ifges(7,5) * t174 + Ifges(7,6) * t173) * t644 + (t406 * t246 + t656 * qJD(1) * (-Ifges(3,2) * t410 + t590)) * t546 / 0.2e1 + (Ifges(3,4) * t314 + Ifges(3,2) * t313 + Ifges(3,6) * t384 + t536) * t556 / 0.2e1 + (t211 * t556 - t212 * t557 - t306 * t511 - t309 * t512 + t313 * t550 - t314 * t330) * mrSges(3,3) + (-t211 * t567 - t314 * t604 - t384 * t550) * mrSges(3,2) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t656 + t211 * t550 + t212 * t330 + t309 * t310) + (-t127 * t236 + t134 * t512 - t55 * t556) * mrSges(5,3) + (Ifges(7,4) * t174 + Ifges(7,2) * t173) * t648 + (Ifges(7,4) * t60 + Ifges(7,2) * t59) * t635 - t677 * t538 - t6 * t460 + t42 * t461 + t314 * (Ifges(3,5) * t567 + (t410 * Ifges(3,1) + t590) * t406) / 0.2e1 + (Ifges(4,4) * t323 - Ifges(4,6) * t556) * t629 + (Ifges(4,4) * t236 + Ifges(4,6) * t512) * t620 + (Ifges(6,1) * t101 - Ifges(6,5) * t512) * t622 + (Ifges(6,1) * t210 + Ifges(6,5) * t556) * t643 + t188 * t465 + (t236 * t714 + t512 * t712) * t610 + (t323 * t714 - t556 * t712) * t615 + (t236 * t715 + t512 * t714) * t617 + (t323 * t715 - t556 * t714) * t630 + (t101 * t107 + t48 * t512 - t556 * t7) * mrSges(6,2) + m(6) * (t107 * t89 + t128 * t42 + t18 * t48 + t19 * t47 + t697 * t7 + t73 * t8) + t697 * t51 - t659 * t442 + (Ifges(7,1) * t174 + Ifges(7,4) * t173) * t649 + (Ifges(7,1) * t60 + Ifges(7,4) * t59) * t633 + (-t166 * t512 + t236 * t249 + t556 * t71) * mrSges(4,2) + (t1 * t173 - t174 * t2 - t20 * t60 + t21 * t59) * mrSges(7,3) + (-m(3) * t306 - t670) * t311 - t328 * t531 + (Ifges(3,1) * t314 + Ifges(3,4) * t313 + Ifges(3,5) * t384) * t509 + t47 * (-mrSges(6,1) * t512 - mrSges(6,3) * t101) + t132 * (-mrSges(5,1) * t512 + mrSges(5,2) * t236) + t165 * (mrSges(4,1) * t512 - mrSges(4,3) * t236) + t695 * t236 / 0.2e1 + (Ifges(6,4) * t101 - Ifges(6,6) * t512) * t624 + (Ifges(6,4) * t210 + Ifges(6,6) * t556) * t642 + (Ifges(5,5) * t323 - Ifges(5,6) * t556) * t628 + (Ifges(5,5) * t236 + Ifges(5,6) * t512) * t619 + t313 * t419 / 0.2e1 + t567 * t520 / 0.2e1 + t384 * (Ifges(3,3) * t567 + (Ifges(3,5) * t410 + Ifges(3,6) * t414) * t406) / 0.2e1 + (t417 / 0.2e1 + t696 * t509) * qJD(2) + m(4) * (t109 * t166 + t110 * t165 + t188 * t294 + t192 * t72 + t193 * t71) + (Ifges(6,5) * t210 + Ifges(6,3) * t556) * t616 + (Ifges(6,5) * t101 - Ifges(6,3) * t512) * t612 + (t212 * t567 + t313 * t604 + t330 * t384) * mrSges(3,1) + (mrSges(2,1) * t606 + mrSges(2,2) * t607 - m(6) * t420 + t143 * mrSges(6,1) - m(7) * (-pkin(5) * t143 + t420) - t740 * mrSges(7,1) - t741 * mrSges(7,2) - m(3) * t471 + t325 * mrSges(3,1) - mrSges(3,3) * t517 - m(4) * (t437 - t600) - m(5) * (t421 - t600) + t718 * t238 - t716 * t237 + t704 * t720 - t668 * t324) * g(1) + m(5) * (t105 * t132 + t106 * t127 + t134 * t98 + t170 * t56 + t171 * t55 + t172 * t57) + m(7) * (t1 * t29 + t15 * t45 + t2 * t28 + t20 * t4 + t21 * t3 + t6 * t69) + Ifges(2,3) * qJDD(1) + t18 * t138 + t19 * t139 + t128 * t26 + (-mrSges(2,1) * t607 + mrSges(2,2) * t606 - m(3) * t548 - t327 * mrSges(3,1) - mrSges(3,3) * t516 - m(7) * (pkin(5) * t149 + t424) - t112 * mrSges(7,1) - t111 * mrSges(7,2) - m(4) * (t519 + t603) - m(5) * (t430 + t603) - m(6) * t424 - t149 * mrSges(6,1) - t718 * t242 + t716 * t241 + t704 * t443 + t668 * t326) * g(2) + t89 * t103 + t3 * t90 + t4 * t91 + t15 * t80 + t73 * t50 + t69 * t13 + t45 * (-mrSges(7,1) * t59 + mrSges(7,2) * t60) + t60 * t54 / 0.2e1 + t59 * t53 / 0.2e1 + t725 * t100 + (mrSges(4,1) * t249 + mrSges(5,1) * t127 - Ifges(4,2) * t620 + Ifges(5,3) * t619 - t610 * t711 + t617 * t713 + t726) * t235 + t29 * t17 + t28 * t16 + (-t55 * mrSges(5,2) - t71 * mrSges(4,3) + t630 * t713 - t721) * t322 + t210 * t650 + t174 * t652 + t173 * t653 + t310 * t304 - t683 * t556 / 0.2e1 + t101 * t638 + t170 * t95 + t171 * t121 + t172 * t119 - t699 * t512 / 0.2e1 + t106 * t190 + t192 * t118 + t193 * t120 + t109 * t213 + t110 * t214 + t105 * t215 + t98 * t216 + t294 * t96 + t56 * t463 + t72 * (-mrSges(4,1) * t556 - mrSges(4,3) * t323) + t8 * (mrSges(6,1) * t556 - mrSges(6,3) * t210) + t57 * (mrSges(5,1) * t556 + mrSges(5,2) * t323) + t323 * t729; (t534 - t179) * t215 + t583 * t628 + t589 * t629 + t733 * (-t127 * (mrSges(5,1) * t409 - mrSges(5,3) * t413) - t249 * (mrSges(4,1) * t409 + mrSges(4,2) * t413)) + (-t582 + t588) * t630 + (-t535 - t200) * t213 + (Ifges(6,5) * t252 - Ifges(6,3) * t514) * t613 + (-t107 * t252 - t48 * t514) * mrSges(6,2) + (-t534 - t199) * t214 + (-t304 + t490) * t306 + (Ifges(6,1) * t252 - Ifges(6,5) * t514) * t623 + (-Ifges(7,1) * t433 - Ifges(7,4) * t434) * t633 + (Ifges(7,1) * t205 + Ifges(7,4) * t204) * t634 + (Ifges(6,4) * t252 - Ifges(6,6) * t514) * t625 + (-Ifges(7,4) * t433 - Ifges(7,2) * t434) * t635 + (Ifges(7,4) * t205 + Ifges(7,2) * t204) * t636 + (-t417 / 0.2e1 + t677 * qJD(1)) * qJD(1) - t543 * t581 + t520 + (-t455 / 0.2e1 + t451 / 0.2e1) * t545 + (t689 * t127 - t132 * t179 - t134 * t176 + t353 * t56) * m(5) + (-pkin(2) * t188 - t165 * t199 - t166 * t200) * m(4) + (-t535 - t176) * t216 - t188 * t464 + t708 * t91 + t709 * t90 + t689 * t190 + (mrSges(7,1) * t693 - mrSges(7,2) * t692) * t45 - t693 * t53 / 0.2e1 + (-t1 * t560 - t2 * t559 + t20 * t692 - t21 * t693) * mrSges(7,3) + t694 * t103 + (t1 * t124 + t123 * t2 + t708 * t20 + t709 * t21 - t441 * t6 + t700 * t45) * m(7) - t710 * t441 + (t107 * t694 + t254 * t7 + t332 * t42 + t441 * t8 + t47 * t701 + t48 * t702) * m(6) + (t491 + t670) * t309 + (-Ifges(7,5) * t433 - Ifges(7,6) * t434) * t626 + (Ifges(7,5) * t205 + Ifges(7,6) * t204) * t627 + ((Ifges(5,6) * t410 + t414 * t451) * t274 + (t410 * t714 + t414 * t686) * t275 - (t410 * t712 + t414 * t687) * t733 + t696 * t410) * t505 + (t275 * t686 - t687 * t733) * qJD(3) / 0.2e1 - t47 * (-mrSges(6,1) * t514 - mrSges(6,3) * t252) + (mrSges(6,1) * t42 + t659) * t334 + (-Ifges(3,2) * t514 + t246 + t379) * t476 + (-t166 * (-mrSges(4,2) * t410 - mrSges(4,3) * t553) - t134 * (-mrSges(5,2) * t553 + mrSges(5,3) * t410) - t165 * (mrSges(4,1) * t410 - mrSges(4,3) * t552) - t132 * (-mrSges(5,1) * t410 + mrSges(5,2) * t552)) * t547 + t123 * t16 + t124 * t17 + ((Ifges(4,6) * t410 + t414 * t455) * t274 + t699 * t410) * t547 / 0.2e1 + (-m(4) * t551 - m(5) * t487 + t328 - m(7) * (pkin(5) * t261 + t422) - (t261 * t411 - t407 * t557) * mrSges(7,1) - (-t261 * t407 - t411 * t557) * mrSges(7,2) - m(6) * t422 - t261 * mrSges(6,1) + mrSges(6,3) * t557 + (t410 * t717 + t414 * t685) * t406 + t704 * (t408 * t527 - t412 * t528)) * g(3) - pkin(2) * t96 + t726 * t544 + ((-t118 + t119) * t409 + ((t132 * t413 - t134 * t409) * qJD(3) + t682) * m(5) + ((-t165 * t413 - t166 * t409) * qJD(3) + t680) * m(4) + (t121 + t120) * t413) * pkin(9) + t724 * t251 + t725 * t224 + (-m(5) * (t498 + t691) - m(4) * t498 + t719 * (-pkin(4) * t562 - t317 + t691) + t662 * t325 + t722 * t324) * g(2) + (-m(5) * (t497 + t690) - m(4) * t497 + t719 * (-pkin(4) * t561 - t319 + t690) + t662 * t327 + t722 * t326) * g(1) + t721 * t413 + t156 * t485 / 0.2e1 + (mrSges(6,2) * t42 - mrSges(6,3) * t8 + t450 * t644 + t453 * t648 + t456 * t649 + t459 * t6 + t499 * t54 + 0.2e1 * t650) * t335 + t353 * t95 + t559 * t652 + t332 * t26 + t680 * mrSges(4,3) + (t132 * t543 + t682) * mrSges(5,2) + t252 * t639 + t700 * t80 + t701 * t139 + t702 * t138 - t205 * t54 / 0.2e1 - t211 * mrSges(3,2) + t212 * mrSges(3,1) + t695 * (t413 * t476 + t543 / 0.2e1) + t254 * t51 + (Ifges(6,1) * t622 + Ifges(6,4) * t624 + Ifges(6,5) * t612 + t524 - t599 + t638 + t706) * t223 - t10 * t560 / 0.2e1 - t56 * t462 + (t153 * t476 + t729 + t739) * t409; (-m(5) * t494 + t719 * (-t241 * pkin(4) + t494) + t716 * t242 + t718 * t241 - t731) * g(1) + (-m(5) * t495 + t719 * (-t237 * pkin(4) + t495) + t716 * t238 + t718 * t237 - t732) * g(2) + t683 + (m(6) * t47 + t735) * (-qJD(5) * t350 + t408 * t734 - t412 * t473) + (-pkin(3) * t57 + qJ(4) * t55 - t127 * t189 - t132 * t166 + t134 * t679) * m(5) - t274 * t581 + (-t107 * t133 + t349 * t8 + t350 * t7 - t48 * t68) * m(6) + (-t20 * t24 - t21 * t25 + t341 * t6) * m(7) + (m(6) * t48 + m(7) * t449 + t674) * (t412 * qJD(4) + qJD(5) * t349) - (-t274 * t715 + t153 + t265 - t578) * t275 / 0.2e1 + t667 * (-pkin(11) + t350) - (Ifges(6,4) * t625 + t599 + t639 + t660) * t183 - t416 + (-t274 * t714 - t275 * t711) * t733 / 0.2e1 + (-Ifges(4,2) * t275 - t266 + t695) * t619 - t68 * t138 - t133 * t103 - pkin(3) * t119 + qJ(4) * t121 + (-m(5) * t493 + t463 + t465 + t719 * (-t315 + t493) + t705 * t442 - t704 * t210) * g(3) - t25 * t90 - t24 * t91 + t156 * t617 + (Ifges(5,3) * t275 - t579) * t620 + t275 * t580 + t134 * t596 + t132 * t597 - t71 * mrSges(4,2) + t72 * mrSges(4,1) + t55 * mrSges(5,3) - t57 * mrSges(5,1) + (t214 - t215) * t166 - t724 * t185 + (-t213 - t216) * t165 + t349 * t50 + t350 * t51 + t341 * t13 - t189 * t190 + qJD(4) * t216 - t127 * (mrSges(5,1) * t275 + mrSges(5,3) * t274) - t249 * (mrSges(4,1) * t275 - mrSges(4,2) * t274); -t202 * t91 - t203 * t90 + t733 * t216 + (-t103 + t190) * t275 + (qJD(5) * t674 + t138 * t733 - t710) * t412 + (t51 + (-t407 * t90 - t411 * t91) * qJD(6) + t678 * t703 + t681) * t408 + t119 + (t428 - t20 * t202 - t203 * t21 + (qJD(5) * t449 - t6) * t412 + (-t45 * t678 + t418) * t408) * m(7) + (-t107 * t275 + t408 * t7 + t412 * t8 + t428 - t678 * (-t408 * t47 + t412 * t48)) * m(6) + (t127 * t275 + t134 * t733 + t428 + t57) * m(5); -pkin(5) * t13 - t30 * t91 - t31 * t90 + t93 * t622 + t416 + (t94 + t175) * t625 + (-t587 + t52) * t623 + (-t138 + t595) * t47 + (-t210 * t533 - t429 * t442 + t461) * g(3) + t732 * g(2) + t731 * g(1) + (-pkin(5) * t6 - t20 * t30 - t21 * t31) * m(7) + (t594 + t735) * t48 + t661 * t185 + t660 * t183 + t667 * pkin(11); -t45 * (mrSges(7,1) * t136 + mrSges(7,2) * t135) + (Ifges(7,1) * t135 - t586) * t634 + t53 * t633 + (Ifges(7,5) * t135 - Ifges(7,6) * t136) * t627 - t20 * t90 + t21 * t91 - g(1) * (mrSges(7,1) * t111 - mrSges(7,2) * t112) - g(2) * (-mrSges(7,1) * t741 + mrSges(7,2) * t740) - g(3) * t460 + (t135 * t20 + t136 * t21) * mrSges(7,3) + t9 + (-Ifges(7,2) * t136 + t131 + t54) * t636 + t671;];
tau  = t22;
