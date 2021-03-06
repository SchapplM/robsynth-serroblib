% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:46
% EndTime: 2019-03-10 02:35:11
% DurationCPUTime: 89.49s
% Computational Cost: add. (45921->1230), mult. (127002->1672), div. (0->0), fcn. (107699->14), ass. (0->564)
t485 = sin(pkin(6));
t491 = sin(qJ(2));
t495 = cos(qJ(2));
t649 = qJD(1) * qJD(2);
t420 = (qJDD(1) * t495 - t491 * t649) * t485;
t699 = cos(pkin(6));
t597 = t699 * qJDD(1);
t467 = t597 + qJDD(2);
t484 = sin(pkin(7));
t486 = cos(pkin(7));
t342 = -t420 * t484 + t467 * t486 + qJDD(3);
t887 = t342 / 0.2e1;
t923 = 0.2e1 * t887;
t494 = cos(qJ(3));
t676 = t491 * t494;
t490 = sin(qJ(3));
t677 = t490 * t495;
t536 = t486 * t677 + t676;
t468 = qJD(1) * t699 + qJD(2);
t690 = t484 * t490;
t640 = t468 * t690;
t661 = qJD(1) * t485;
t421 = (qJDD(1) * t491 + t495 * t649) * t485;
t543 = t420 * t486 + t467 * t484;
t809 = -t490 * t421 + t543 * t494;
t230 = (-t536 * t661 - t640) * qJD(3) + t809;
t888 = t230 / 0.2e1;
t922 = 0.2e1 * t888;
t674 = t494 * t495;
t678 = t490 * t491;
t538 = t486 * t674 - t678;
t687 = t484 * t494;
t329 = t468 * t687 + t538 * t661;
t229 = qJD(3) * t329 + t494 * t421 + t490 * t543;
t889 = t229 / 0.2e1;
t921 = 0.2e1 * t889;
t633 = pkin(1) * t699;
t478 = t495 * t633;
t464 = qJD(1) * t478;
t749 = pkin(10) * t486;
t634 = pkin(9) + t749;
t579 = t485 * t634;
t541 = t491 * t579;
t361 = -qJD(1) * t541 + t464;
t477 = t491 * t633;
t507 = -t495 * t579 - t477;
t362 = t507 * qJD(1);
t750 = pkin(10) * t484;
t525 = t485 * (pkin(2) * t491 - t495 * t750);
t401 = qJD(1) * t525;
t469 = pkin(10) * t690;
t683 = t486 * t494;
t435 = pkin(2) * t683 - t469;
t684 = t486 * t490;
t834 = t435 * qJD(3) - t494 * t361 - t362 * t684 - t401 * t690;
t284 = -t362 * t484 + t486 * t401;
t537 = t486 * t676 + t677;
t396 = t537 * t485;
t379 = qJD(1) * t396;
t535 = -t486 * t678 + t674;
t397 = t535 * t485;
t380 = qJD(1) * t397;
t659 = qJD(3) * t484;
t920 = -pkin(3) * t379 + pkin(11) * t380 - t284 + (pkin(3) * t490 - pkin(11) * t494) * t659;
t630 = t491 * t661;
t589 = t484 * t630;
t919 = pkin(11) * t589 - t834;
t519 = t536 * t485;
t330 = qJD(1) * t519 + t640;
t629 = t495 * t661;
t384 = t486 * t468 - t484 * t629 + qJD(3);
t489 = sin(qJ(4));
t493 = cos(qJ(4));
t271 = t330 * t493 + t384 * t489;
t136 = -qJD(4) * t271 - t229 * t489 + t342 * t493;
t132 = Ifges(5,6) * t136;
t270 = -t330 * t489 + t384 * t493;
t135 = qJD(4) * t270 + t229 * t493 + t342 * t489;
t133 = Ifges(5,5) * t135;
t226 = qJD(3) * t330 + qJDD(4) - t809;
t222 = Ifges(5,3) * t226;
t47 = t133 + t132 + t222;
t328 = pkin(2) * t468 + t361;
t689 = t484 * t491;
t375 = (-pkin(2) * t495 - pkin(10) * t689 - pkin(1)) * t661;
t262 = -t328 * t484 + t486 * t375;
t164 = -pkin(3) * t329 - pkin(11) * t330 + t262;
t685 = t485 * t495;
t323 = t468 * t750 + (t634 * t685 + t477) * qJD(1);
t545 = t328 * t486 + t375 * t484;
t194 = t323 * t494 + t490 * t545;
t173 = pkin(11) * t384 + t194;
t655 = qJD(4) * t493;
t656 = qJD(4) * t489;
t474 = pkin(9) * t685;
t587 = qJD(2) * t633;
t542 = qJD(1) * t587;
t582 = pkin(1) * t597;
t660 = qJD(2) * t485;
t628 = t491 * t660;
t596 = pkin(9) * t628;
t326 = -qJD(1) * t596 + qJDD(1) * t474 + t491 * t582 + t495 * t542;
t263 = pkin(10) * t543 + t326;
t327 = -pkin(9) * t421 - t491 * t542 + t495 * t582;
t269 = t467 * pkin(2) - t421 * t749 + t327;
t754 = pkin(1) * t485;
t642 = qJDD(1) * t754;
t306 = -pkin(2) * t420 - t421 * t750 - t642;
t657 = qJD(3) * t494;
t624 = t486 * t657;
t626 = t484 * t657;
t658 = qJD(3) * t490;
t84 = t494 * t263 + t269 * t684 + t306 * t690 - t323 * t658 + t328 * t624 + t375 * t626;
t78 = pkin(11) * t342 + t84;
t188 = -t269 * t484 + t486 * t306;
t86 = -pkin(3) * t230 - pkin(11) * t229 + t188;
t19 = t164 * t655 - t173 * t656 + t489 * t86 + t493 * t78;
t20 = -t164 * t656 - t173 * t655 - t489 * t78 + t493 * t86;
t574 = -t20 * mrSges(5,1) + t19 * mrSges(5,2);
t918 = Ifges(4,4) * t921 + Ifges(4,2) * t922 + Ifges(4,6) * t923 + t574 + t84 * mrSges(4,3) - t47 / 0.2e1;
t317 = t380 * t489 - t493 * t589;
t688 = t484 * t493;
t429 = t486 * t489 + t490 * t688;
t366 = qJD(4) * t429 + t489 * t626;
t668 = t317 - t366;
t318 = t380 * t493 + t489 * t589;
t428 = -t493 * t486 + t489 * t690;
t365 = -qJD(4) * t428 + t493 * t626;
t667 = t318 - t365;
t437 = pkin(2) * t684 + pkin(10) * t687;
t917 = -t437 * qJD(3) + t490 * t361;
t627 = t484 * t658;
t916 = t379 - t627;
t325 = qJD(4) - t329;
t864 = Ifges(6,4) + Ifges(7,4);
t792 = t135 / 0.2e1;
t915 = Ifges(5,4) * t792;
t837 = t362 * t683 - (-pkin(3) * t630 - t401 * t494) * t484 - t917;
t407 = pkin(11) * t486 + t437;
t408 = (-pkin(3) * t494 - pkin(11) * t490 - pkin(2)) * t484;
t833 = -t407 * t656 + t408 * t655 + t489 * t920 - t919 * t493;
t791 = t136 / 0.2e1;
t781 = t226 / 0.2e1;
t865 = Ifges(6,1) + Ifges(7,1);
t863 = -Ifges(7,5) - Ifges(6,5);
t862 = Ifges(6,2) + Ifges(7,2);
t861 = Ifges(6,6) + Ifges(7,6);
t756 = cos(qJ(1));
t578 = t699 * t756;
t755 = sin(qJ(1));
t524 = t755 * t491 - t495 * t578;
t632 = t485 * t756;
t504 = t524 * t484 - t486 * t632;
t914 = t489 * t504;
t913 = t493 * t504;
t912 = pkin(12) * t916 - t833;
t911 = -t194 + t325 * (pkin(4) * t489 - pkin(12) * t493);
t910 = -pkin(4) * t668 + t667 * pkin(12) + t837;
t134 = qJDD(5) - t136;
t488 = sin(qJ(5));
t492 = cos(qJ(5));
t210 = t271 * t492 + t325 * t488;
t11 = pkin(12) * t226 + t19;
t625 = t486 * t658;
t85 = t494 * (t269 * t486 + t306 * t484) - t490 * t263 - t323 * t657 - t328 * t625 - t375 * t627;
t79 = -pkin(3) * t342 - t85;
t27 = -pkin(4) * t136 - pkin(12) * t135 + t79;
t193 = -t490 * t323 + t494 * t545;
t172 = -pkin(3) * t384 - t193;
t105 = -pkin(4) * t270 - pkin(12) * t271 + t172;
t92 = t164 * t489 + t173 * t493;
t81 = pkin(12) * t325 + t92;
t41 = t105 * t488 + t492 * t81;
t4 = -qJD(5) * t41 - t11 * t488 + t492 * t27;
t209 = -t271 * t488 + t325 * t492;
t70 = qJD(5) * t209 + t135 * t492 + t226 * t488;
t1 = pkin(5) * t134 - qJ(6) * t70 - qJD(6) * t210 + t4;
t793 = t134 / 0.2e1;
t71 = -qJD(5) * t210 - t135 * t488 + t226 * t492;
t797 = t71 / 0.2e1;
t798 = t70 / 0.2e1;
t860 = -Ifges(7,3) - Ifges(6,3);
t652 = qJD(5) * t492;
t654 = qJD(5) * t488;
t3 = t105 * t652 + t492 * t11 + t488 * t27 - t654 * t81;
t2 = qJ(6) * t71 + qJD(6) * t209 + t3;
t874 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t908 = t1 * mrSges(7,1) - 0.2e1 * Ifges(5,2) * t791 - 0.2e1 * Ifges(5,6) * t781 - t860 * t793 + t861 * t797 - t863 * t798 + t874 - t915;
t875 = Ifges(5,1) * t792 + Ifges(5,5) * t781;
t799 = Ifges(5,4) * t791 + t875;
t907 = -mrSges(5,3) * t20 + 0.2e1 * t799;
t906 = -t85 * mrSges(4,3) + Ifges(4,1) * t921 + Ifges(4,4) * t922 + Ifges(4,5) * t923;
t224 = Ifges(4,6) * t230;
t225 = Ifges(4,5) * t229;
t336 = Ifges(4,3) * t342;
t117 = t225 + t224 + t336;
t575 = t85 * mrSges(4,1) - t84 * mrSges(4,2);
t905 = t575 + t117 / 0.2e1;
t903 = t864 * t797 + t865 * t798;
t482 = pkin(5) * t492 + pkin(4);
t900 = m(7) * t482;
t897 = t864 * t209;
t832 = -t407 * t655 - t408 * t656 + t919 * t489 + t493 * t920;
t681 = t488 * t493;
t245 = -t329 * t681 + t330 * t492;
t623 = t488 * t655;
t831 = t489 * t652 + t245 + t623;
t367 = -t429 * t488 - t492 * t687;
t258 = qJD(5) * t367 + t365 * t492 + t488 * t627;
t261 = t318 * t492 + t379 * t488;
t672 = t258 - t261;
t539 = -t429 * t492 + t488 * t687;
t259 = qJD(5) * t539 - t365 * t488 + t492 * t627;
t260 = -t318 * t488 + t379 * t492;
t671 = t259 - t260;
t698 = t270 * t488;
t896 = t654 - t698;
t895 = t864 * t210;
t894 = t262 * mrSges(4,2) - t193 * mrSges(4,3);
t91 = t164 * t493 - t489 * t173;
t893 = t262 * mrSges(4,1) + t91 * mrSges(5,1) - t92 * mrSges(5,2) - t194 * mrSges(4,3);
t731 = Ifges(5,4) * t271;
t155 = t270 * Ifges(5,2) + t325 * Ifges(5,6) + t731;
t268 = qJD(5) - t270;
t853 = t209 * t861 - t210 * t863 - t268 * t860;
t878 = t853 / 0.2e1 - t155 / 0.2e1;
t487 = -qJ(6) - pkin(12);
t891 = -m(7) * t487 + mrSges(6,3) + mrSges(7,3);
t886 = -t863 * t134 / 0.2e1 + t903;
t858 = t134 * t861 + t70 * t864 + t71 * t862;
t885 = t858 / 0.2e1;
t406 = t469 + (-pkin(2) * t494 - pkin(3)) * t486;
t295 = pkin(4) * t428 - pkin(12) * t429 + t406;
t310 = t493 * t407 + t489 * t408;
t298 = -pkin(12) * t687 + t310;
t854 = t295 * t652 - t298 * t654 + t488 * t910 - t492 * t912;
t852 = t209 * t862 + t268 * t861 + t895;
t851 = t210 * t865 - t268 * t863 + t897;
t208 = t488 * t295 + t492 * t298;
t839 = -qJD(5) * t208 + t488 * t912 + t492 * t910;
t252 = pkin(3) * t330 - pkin(11) * t329;
t142 = t493 * t193 + t489 * t252;
t116 = pkin(12) * t330 + t142;
t646 = pkin(11) * t656;
t884 = t911 * t492 + (t116 + t646) * t488;
t573 = -pkin(4) * t493 - pkin(12) * t489;
t453 = -pkin(3) + t573;
t883 = -t492 * t116 + t453 * t652 + t488 * t911;
t141 = -t489 * t193 + t252 * t493;
t115 = -pkin(4) * t330 - t141;
t645 = pkin(11) * t655;
t882 = -t115 + t645;
t836 = pkin(4) * t916 - t832;
t577 = t699 * t755;
t508 = t491 * t756 + t495 * t577;
t631 = t485 * t755;
t501 = t508 * t484 + t486 * t631;
t762 = -t384 / 0.2e1;
t766 = -t330 / 0.2e1;
t881 = Ifges(4,1) * t766 + Ifges(4,5) * t762 - t894;
t880 = t864 * t492;
t879 = t864 * t488;
t767 = -t329 / 0.2e1;
t769 = -t325 / 0.2e1;
t772 = -t271 / 0.2e1;
t774 = -t270 / 0.2e1;
t876 = Ifges(5,5) * t772 - Ifges(4,2) * t767 - Ifges(4,6) * t762 + Ifges(5,6) * t774 + Ifges(5,3) * t769 - t893;
t776 = -t268 / 0.2e1;
t785 = -t210 / 0.2e1;
t787 = -t209 / 0.2e1;
t873 = -t860 * t776 - t863 * t785 + t861 * t787;
t40 = t492 * t105 - t488 * t81;
t28 = -qJ(6) * t210 + t40;
t24 = pkin(5) * t268 + t28;
t29 = qJ(6) * t209 + t41;
t872 = t172 * mrSges(5,1) + t40 * mrSges(6,1) + t24 * mrSges(7,1) - t41 * mrSges(6,2) - t29 * mrSges(7,2);
t871 = t3 * mrSges(6,3) + t2 * mrSges(7,3) + t861 * t793 + t862 * t797 + t864 * t798 + t885;
t870 = -t4 * mrSges(6,3) - t1 * mrSges(7,3) - t863 * t793 + t886 + t903;
t859 = -t134 * t860 - t70 * t863 + t71 * t861;
t868 = t859 / 0.2e1 - t19 * mrSges(5,3) - t915 + t908;
t801 = m(7) * pkin(5);
t741 = -mrSges(6,1) - mrSges(7,1);
t867 = mrSges(4,2) - mrSges(5,3);
t866 = mrSges(6,2) + mrSges(7,2);
t856 = -pkin(5) * t668 - qJ(6) * t672 + qJD(6) * t539 + t839;
t855 = qJ(6) * t671 + qJD(6) * t367 + t854;
t848 = mrSges(7,1) + t801;
t847 = -pkin(5) * t671 + t836;
t675 = t492 * t493;
t246 = t329 * t675 + t330 * t488;
t479 = pkin(11) * t675;
t651 = qJD(6) * t492;
t694 = t329 * t489;
t846 = -pkin(5) * t694 + qJ(6) * t246 - t489 * t651 + (pkin(5) * t489 - qJ(6) * t675) * qJD(4) + (-t479 + (qJ(6) * t489 - t453) * t488) * qJD(5) + t884;
t680 = t489 * t492;
t845 = -qJ(6) * t245 + (-pkin(11) * qJD(4) - qJ(6) * qJD(5)) * t680 + (-qJD(6) * t489 + (-pkin(11) * qJD(5) - qJ(6) * qJD(4)) * t493) * t488 + t883;
t844 = (-t492 * t656 - t493 * t654) * pkin(11) + t883;
t404 = t488 * t453 + t479;
t843 = -qJD(5) * t404 + t884;
t842 = pkin(5) * t831 + t882;
t178 = pkin(4) * t271 - pkin(12) * t270;
t60 = t492 * t178 - t488 * t91;
t600 = qJD(5) * t487;
t697 = t270 * t492;
t841 = -pkin(5) * t271 + qJ(6) * t697 - qJD(6) * t488 + t492 * t600 - t60;
t61 = t488 * t178 + t492 * t91;
t840 = qJ(6) * t698 + t488 * t600 - t61 + t651;
t838 = pkin(5) * t896 - t92;
t360 = pkin(2) * t699 + t478 - t541;
t686 = t485 * t491;
t641 = t484 * t686;
t664 = pkin(2) * t685 + pkin(10) * t641;
t390 = -t664 - t754;
t280 = -t360 * t484 + t486 * t390;
t599 = t699 * t484;
t355 = -t485 * t538 - t494 * t599;
t356 = t490 * t599 + t519;
t601 = -t355 * pkin(3) + pkin(11) * t356;
t189 = t280 - t601;
t663 = t474 + t477;
t345 = (t486 * t685 + t599) * pkin(10) + t663;
t228 = t494 * t345 + t360 * t684 + t390 * t690;
t520 = t484 * t685 - t486 * t699;
t196 = -pkin(11) * t520 + t228;
t111 = t489 * t189 + t493 * t196;
t108 = pkin(12) * t355 + t111;
t227 = -t490 * t345 + t494 * (t360 * t486 + t390 * t484);
t195 = pkin(3) * t520 - t227;
t293 = t356 * t489 + t493 * t520;
t294 = t356 * t493 - t489 * t520;
t125 = t293 * pkin(4) - t294 * pkin(12) + t195;
t53 = t492 * t108 + t488 * t125;
t835 = -(t362 * t486 + t401 * t484) * t494 + t917;
t653 = qJD(5) * t489;
t830 = t488 * t653 - t492 * t655 + t246;
t564 = -mrSges(7,1) * t492 + mrSges(7,2) * t488;
t567 = -mrSges(6,1) * t492 + mrSges(6,2) * t488;
t829 = m(6) * pkin(4) - t564 - t567 + t900;
t563 = mrSges(7,1) * t488 + mrSges(7,2) * t492;
t566 = mrSges(6,1) * t488 + mrSges(6,2) * t492;
t80 = -pkin(4) * t325 - t91;
t57 = -pkin(5) * t209 + qJD(6) + t80;
t828 = t57 * t563 + t80 * t566;
t827 = -t488 * t863 + t492 * t861;
t826 = -t488 * t861 - t492 * t863;
t825 = t492 * t862 + t879;
t824 = -t488 * t862 + t880;
t823 = t488 * t865 + t880;
t822 = t492 * t865 - t879;
t821 = -m(6) * pkin(12) - t891;
t820 = -t652 + t697;
t818 = t656 - t694;
t817 = t19 * t493 - t20 * t489;
t816 = t3 * t492 - t4 * t488;
t650 = -m(5) - m(6) - m(7);
t815 = -mrSges(6,1) - t848;
t814 = mrSges(5,1) + t829;
t813 = mrSges(5,2) + t821;
t569 = -mrSges(5,1) * t493 + mrSges(5,2) * t489;
t812 = -m(6) * t573 + t489 * t891 + t493 * t900 - t569;
t808 = mrSges(4,1) + t812;
t751 = pkin(5) * t488;
t635 = pkin(11) + t751;
t807 = -m(7) * t635 + t867;
t802 = t485 ^ 2;
t711 = t325 * Ifges(5,3);
t712 = t271 * Ifges(5,5);
t713 = t270 * Ifges(5,6);
t154 = t711 + t712 + t713;
t790 = t154 / 0.2e1;
t786 = t209 / 0.2e1;
t784 = t210 / 0.2e1;
t702 = t384 * Ifges(4,6);
t707 = t330 * Ifges(4,4);
t710 = t329 * Ifges(4,2);
t220 = t702 + t707 + t710;
t783 = -t220 / 0.2e1;
t324 = Ifges(4,4) * t329;
t703 = t384 * Ifges(4,5);
t708 = t330 * Ifges(4,1);
t221 = t324 + t703 + t708;
t782 = t221 / 0.2e1;
t775 = t268 / 0.2e1;
t773 = t270 / 0.2e1;
t771 = t271 / 0.2e1;
t768 = t325 / 0.2e1;
t765 = t330 / 0.2e1;
t759 = t491 / 0.2e1;
t753 = pkin(5) * t210;
t235 = -t294 * t488 + t355 * t492;
t752 = pkin(5) * t235;
t739 = mrSges(5,3) * t270;
t738 = mrSges(5,3) * t271;
t737 = mrSges(6,3) * t209;
t736 = mrSges(6,3) * t210;
t735 = mrSges(7,3) * t209;
t734 = mrSges(7,3) * t210;
t733 = Ifges(3,4) * t491;
t732 = Ifges(3,4) * t495;
t730 = Ifges(5,4) * t489;
t729 = Ifges(5,4) * t493;
t12 = -pkin(4) * t226 - t20;
t724 = t12 * t489;
t709 = t329 * Ifges(4,6);
t706 = t330 * Ifges(4,5);
t701 = t384 * Ifges(4,3);
t700 = t495 * Ifges(3,2);
t514 = t524 * t490;
t430 = t491 * t578 + t495 * t755;
t593 = t484 * t632;
t598 = -t430 * t494 + t490 * t593;
t300 = -t486 * t514 - t598;
t696 = t300 * t488;
t431 = -t491 * t577 + t495 * t756;
t591 = t484 * t631;
t304 = t431 * t494 + (-t486 * t508 + t591) * t490;
t695 = t304 * t488;
t693 = t329 * t493;
t692 = t356 * t488;
t691 = t484 * t489;
t682 = t488 * t489;
t463 = Ifges(3,4) * t629;
t673 = t495 * (Ifges(3,1) * t630 + t468 * Ifges(3,5) + t463);
t513 = t524 * t494;
t320 = -t430 * t684 - t513;
t422 = t524 * pkin(2);
t670 = t320 * pkin(3) - t422;
t505 = t508 * t494;
t322 = -t431 * t684 - t505;
t424 = t508 * pkin(2);
t669 = t322 * pkin(3) - t424;
t662 = t756 * pkin(1) + pkin(9) * t631;
t644 = mrSges(3,3) * t685;
t637 = t397 * pkin(3) + t664;
t636 = Ifges(3,5) * t421 + Ifges(3,6) * t420 + Ifges(3,3) * t467;
t22 = -t71 * mrSges(7,1) + t70 * mrSges(7,2);
t613 = t802 * t649;
t610 = t655 / 0.2e1;
t605 = -t653 / 0.2e1;
t52 = -t108 * t488 + t492 * t125;
t110 = t189 * t493 - t489 * t196;
t207 = t492 * t295 - t298 * t488;
t364 = t507 * qJD(2);
t402 = qJD(2) * t525;
t285 = -t364 * t484 + t486 * t402;
t309 = -t489 * t407 + t408 * t493;
t465 = t495 * t587;
t363 = -qJD(2) * t541 + t465;
t594 = -t345 * t657 - t360 * t625 - t490 * t363 - t390 * t627;
t588 = t484 * t628;
t580 = -pkin(1) * t755 + pkin(9) * t632;
t576 = qJD(3) * t599;
t297 = pkin(4) * t687 - t309;
t571 = mrSges(4,1) * t355 + mrSges(4,2) * t356;
t570 = mrSges(5,1) * t293 + mrSges(5,2) * t294;
t236 = t294 * t492 + t355 * t488;
t568 = mrSges(6,1) * t235 - mrSges(6,2) * t236;
t565 = -t235 * mrSges(7,1) + t236 * mrSges(7,2);
t562 = Ifges(5,1) * t493 - t730;
t557 = -Ifges(5,2) * t489 + t729;
t552 = Ifges(5,5) * t493 - Ifges(5,6) * t489;
t244 = t304 * t493 + t489 * t501;
t303 = t431 * t490 + t486 * t505 - t494 * t591;
t170 = -t244 * t488 + t303 * t492;
t145 = -t345 * t658 + t360 * t624 + t494 * t363 + t364 * t684 + t390 * t626 + t402 * t690;
t137 = pkin(11) * t588 + t145;
t282 = t490 * t576 + (qJD(2) * t537 + qJD(3) * t536) * t485;
t283 = t494 * t576 + (qJD(2) * t535 + qJD(3) * t538) * t485;
t158 = pkin(3) * t282 - pkin(11) * t283 + t285;
t39 = -t489 * t137 + t158 * t493 - t189 * t656 - t196 * t655;
t107 = -pkin(4) * t355 - t110;
t540 = pkin(1) * (mrSges(3,1) * t491 + mrSges(3,2) * t495);
t532 = t172 * (mrSges(5,1) * t489 + mrSges(5,2) * t493);
t531 = t468 * (Ifges(3,5) * t495 - Ifges(3,6) * t491);
t530 = t491 * (Ifges(3,1) * t495 - t733);
t529 = (t700 + t733) * t485;
t38 = t493 * t137 + t489 * t158 + t189 * t655 - t196 * t656;
t31 = pkin(12) * t282 + t38;
t138 = -t364 * t683 + (-pkin(3) * t628 - t402 * t494) * t484 - t594;
t181 = qJD(4) * t294 + t283 * t489 - t493 * t588;
t182 = -qJD(4) * t293 + t283 * t493 + t489 * t588;
t56 = pkin(4) * t181 - pkin(12) * t182 + t138;
t7 = -t108 * t654 + t125 * t652 + t492 * t31 + t488 * t56;
t32 = -pkin(4) * t282 - t39;
t517 = -t430 * t490 - t494 * t593;
t516 = mrSges(3,2) + (-mrSges(4,3) + (-m(4) + t650) * pkin(10)) * t484;
t515 = t486 * t524;
t511 = -m(7) * t751 + pkin(11) * t650 + t867;
t8 = -qJD(5) * t53 - t31 * t488 + t492 * t56;
t502 = -t430 * pkin(2) - pkin(10) * t504 + t580;
t500 = t431 * pkin(2) + pkin(10) * t501 + t662;
t302 = t490 * t515 + t598;
t499 = t302 * pkin(3) + t502;
t498 = t304 * pkin(3) + t500;
t301 = -t494 * t515 + t517;
t497 = t301 * pkin(11) + t499;
t496 = t303 * pkin(11) + t498;
t456 = t487 * t492;
t455 = t487 * t488;
t449 = t635 * t489;
t443 = t492 * t453;
t436 = -pkin(9) * t686 + t478;
t432 = (-mrSges(3,1) * t495 + mrSges(3,2) * t491) * t485;
t419 = t663 * qJD(2);
t418 = t465 - t596;
t415 = t663 * qJD(1);
t413 = -pkin(9) * t630 + t464;
t412 = -t468 * mrSges(3,2) + mrSges(3,3) * t629;
t411 = mrSges(3,1) * t468 - mrSges(3,3) * t630;
t403 = -pkin(11) * t681 + t443;
t371 = Ifges(3,6) * t468 + qJD(1) * t529;
t370 = -qJ(6) * t682 + t404;
t354 = -qJ(6) * t680 + t443 + (-pkin(11) * t488 - pkin(5)) * t493;
t332 = t397 * t493 + t489 * t641;
t321 = t431 * t683 - t490 * t508;
t319 = t430 * t683 - t514;
t299 = t486 * t513 - t517;
t279 = mrSges(4,1) * t384 - mrSges(4,3) * t330;
t278 = -mrSges(4,2) * t384 + mrSges(4,3) * t329;
t275 = t322 * t493 + t431 * t691;
t273 = t320 * t493 + t430 * t691;
t267 = Ifges(5,4) * t270;
t251 = -mrSges(4,1) * t329 + mrSges(4,2) * t330;
t247 = -pkin(5) * t367 + t297;
t243 = t304 * t489 - t493 * t501;
t242 = t302 * t493 - t914;
t240 = t300 * t493 + t914;
t239 = t300 * t489 - t913;
t219 = t701 + t706 + t709;
t213 = mrSges(5,1) * t325 - t738;
t212 = -mrSges(5,2) * t325 + t739;
t180 = -mrSges(4,2) * t342 + mrSges(4,3) * t230;
t179 = mrSges(4,1) * t342 - mrSges(4,3) * t229;
t177 = -mrSges(5,1) * t270 + mrSges(5,2) * t271;
t174 = qJ(6) * t367 + t208;
t171 = t244 * t492 + t303 * t488;
t160 = pkin(5) * t428 + qJ(6) * t539 + t207;
t156 = t271 * Ifges(5,1) + t325 * Ifges(5,5) + t267;
t153 = mrSges(6,1) * t268 - t736;
t152 = mrSges(7,1) * t268 - t734;
t151 = -mrSges(6,2) * t268 + t737;
t150 = -mrSges(7,2) * t268 + t735;
t146 = (t364 * t486 + t402 * t484) * t494 + t594;
t144 = -mrSges(4,1) * t230 + mrSges(4,2) * t229;
t124 = -mrSges(6,1) * t209 + mrSges(6,2) * t210;
t123 = -mrSges(7,1) * t209 + mrSges(7,2) * t210;
t102 = qJD(5) * t235 + t182 * t492 + t282 * t488;
t101 = -qJD(5) * t236 - t182 * t488 + t282 * t492;
t89 = -mrSges(5,2) * t226 + mrSges(5,3) * t136;
t88 = mrSges(5,1) * t226 - mrSges(5,3) * t135;
t75 = t107 - t752;
t72 = -mrSges(5,1) * t136 + mrSges(5,2) * t135;
t42 = qJ(6) * t235 + t53;
t37 = -mrSges(6,2) * t134 + mrSges(6,3) * t71;
t36 = -mrSges(7,2) * t134 + mrSges(7,3) * t71;
t35 = mrSges(6,1) * t134 - mrSges(6,3) * t70;
t34 = mrSges(7,1) * t134 - mrSges(7,3) * t70;
t33 = pkin(5) * t293 - qJ(6) * t236 + t52;
t23 = -mrSges(6,1) * t71 + mrSges(6,2) * t70;
t21 = -pkin(5) * t101 + t32;
t9 = -pkin(5) * t71 + qJDD(6) + t12;
t6 = qJ(6) * t101 + qJD(6) * t235 + t7;
t5 = pkin(5) * t181 - qJ(6) * t102 - qJD(6) * t236 + t8;
t10 = [t870 * t236 + t871 * t235 + t868 * t293 + (-t863 * t775 + t864 * t786 + t865 * t784 + t851 / 0.2e1 - mrSges(7,3) * t24 - mrSges(6,3) * t40 + mrSges(7,2) * t57 + mrSges(6,2) * t80) * t102 + (Ifges(5,1) * t182 + Ifges(5,5) * t282) * t771 + (t861 * t775 + t862 * t786 + t864 * t784 + t852 / 0.2e1 + mrSges(7,3) * t29 + mrSges(6,3) * t41 - mrSges(7,1) * t57 - mrSges(6,1) * t80) * t101 + (Ifges(5,4) * t182 + Ifges(5,6) * t282) * t773 + (t495 * (-Ifges(3,2) * t491 + t732) + t530) * t613 / 0.2e1 + (t673 + t531) * t660 / 0.2e1 + (-t92 * mrSges(5,3) - Ifges(5,4) * t771 - Ifges(5,2) * t773 - Ifges(5,6) * t768 - t860 * t775 - t863 * t784 + t861 * t786 + t872 + t878) * t181 + (t172 * t182 - t282 * t92) * mrSges(5,2) + t663 * (-mrSges(3,2) * t467 + mrSges(3,3) * t420) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t802 + t326 * t663 + t327 * t436 - t413 * t419 + t415 * t418) + (-mrSges(2,1) * t756 + mrSges(2,2) * t755 - m(6) * (t244 * pkin(4) + t496) - m(4) * t500 - t304 * mrSges(4,1) - mrSges(4,3) * t501 - m(5) * t496 - t244 * mrSges(5,1) - m(7) * (t244 * t482 + t498) - m(3) * t662 - t431 * mrSges(3,1) + mrSges(3,2) * t508 - mrSges(3,3) * t631 + t807 * t303 + t741 * t171 - t866 * t170 + t813 * t243) * g(2) + (Ifges(3,4) * t421 + Ifges(3,2) * t420 + Ifges(3,6) * t467) * t685 / 0.2e1 + (Ifges(3,1) * t421 + Ifges(3,4) * t420 + Ifges(3,5) * t467) * t686 / 0.2e1 + t9 * t565 - t12 * t568 + t79 * t570 + t188 * t571 + m(6) * (t107 * t12 + t3 * t53 + t32 * t80 + t4 * t52 + t40 * t8 + t41 * t7) + m(7) * (t1 * t33 + t2 * t42 + t21 * t57 + t24 * t5 + t29 * t6 + t75 * t9) + m(4) * (t145 * t194 + t146 * t193 + t188 * t280 + t227 * t85 + t228 * t84 + t262 * t285) + m(5) * (t110 * t20 + t111 * t19 + t138 * t172 + t195 * t79 + t38 * t92 + t39 * t91) - (-mrSges(3,1) * t420 + mrSges(3,2) * t421) * t754 + t421 * (Ifges(3,5) * t699 + (t491 * Ifges(3,1) + t732) * t485) / 0.2e1 - t432 * t642 + t327 * (mrSges(3,1) * t699 - mrSges(3,3) * t686) + t326 * (-mrSges(3,2) * t699 + t644) + t699 * t636 / 0.2e1 + t420 * (Ifges(3,6) * t699 + t529) / 0.2e1 + t467 * (Ifges(3,3) * t699 + (Ifges(3,5) * t491 + Ifges(3,6) * t495) * t485) / 0.2e1 + (-m(4) * t502 - t302 * mrSges(4,1) + mrSges(4,3) * t504 + mrSges(2,1) * t755 + mrSges(2,2) * t756 - m(3) * t580 + t430 * mrSges(3,1) - mrSges(3,2) * t524 - mrSges(3,3) * t632 - m(6) * (t242 * pkin(4) + t497) - m(5) * t497 - t242 * mrSges(5,1) - m(7) * (t242 * t482 + t499) + t807 * t301 + t741 * (t242 * t492 + t301 * t488) - t866 * (-t242 * t488 + t301 * t492) + t813 * (t302 * t489 + t913)) * g(1) + t906 * t356 + t907 * t294 + (-t371 / 0.2e1 - t415 * mrSges(3,3)) * t628 + (Ifges(5,5) * t182 + Ifges(5,3) * t282) * t768 - t540 * t613 + (-Ifges(4,5) * t889 - Ifges(4,6) * t888 - Ifges(4,3) * t887 - t905) * t520 + Ifges(2,3) * qJDD(1) + t262 * (mrSges(4,1) * t282 + mrSges(4,2) * t283) + t285 * t251 + t280 * t144 + t91 * (mrSges(5,1) * t282 - mrSges(5,3) * t182) + t145 * t278 + t146 * t279 + t227 * t179 + t228 * t180 - t413 * qJD(2) * t644 + t38 * t212 + t39 * t213 + t219 * t588 / 0.2e1 + t194 * (-mrSges(4,2) * t588 - mrSges(4,3) * t282) + t384 * (Ifges(4,5) * t283 - Ifges(4,6) * t282 + Ifges(4,3) * t588) / 0.2e1 + t329 * (Ifges(4,4) * t283 - Ifges(4,2) * t282 + Ifges(4,6) * t588) / 0.2e1 + t195 * t72 + t193 * (mrSges(4,1) * t588 - mrSges(4,3) * t283) + t182 * t156 / 0.2e1 + t138 * t177 + t7 * t151 + t5 * t152 + t8 * t153 + t6 * t150 + t21 * t123 + t32 * t124 + t111 * t89 + t418 * t412 - t419 * t411 + t42 * t36 + t52 * t35 + t53 * t37 + t436 * (mrSges(3,1) * t467 - mrSges(3,3) * t421) + (Ifges(5,5) * t792 + Ifges(5,6) * t791 + Ifges(5,3) * t781 - t918) * t355 + (Ifges(4,1) * t283 - Ifges(4,4) * t282 + Ifges(4,5) * t588) * t765 + t283 * t782 + t282 * t783 + t282 * t790 + t75 * t22 + t107 * t23 + t110 * t88 + t33 * t34; (-t12 * mrSges(6,2) - t9 * mrSges(7,2) - t870) * t539 + (-t12 * mrSges(6,1) - t9 * mrSges(7,1) + t871) * t367 + (t260 * t861 - t261 * t863 - t317 * t860) * t776 + (-t258 * t863 + t259 * t861 - t366 * t860) * t775 + (t260 * t862 + t261 * t864 + t317 * t861) * t787 + (t258 * t864 + t259 * t862 + t366 * t861) * t786 + (t258 * t865 + t259 * t864 - t366 * t863) * t784 + (t260 * t864 + t261 * t865 - t317 * t863) * t785 + t851 * (-t261 / 0.2e1 + t258 / 0.2e1) + t852 * (-t260 / 0.2e1 + t259 / 0.2e1) + (-t155 + t853) * (t366 / 0.2e1 - t317 / 0.2e1) + (t12 * t297 + t207 * t4 + t208 * t3 + (-pkin(4) * t275 - t669) * g(1) + (-pkin(4) * t273 - t670) * g(2) + (-pkin(4) * t332 - t637) * g(3) + t836 * t80 + t854 * t41 + t839 * t40) * m(6) + t854 * t151 + t855 * t150 + (t1 * t160 + t174 * t2 + t247 * t9 + (-t275 * t482 - t669) * g(1) + (-t273 * t482 - t670) * g(2) + (-t332 * t482 - t637) * g(3) + t847 * t57 + t855 * t29 + t856 * t24) * m(7) + t856 * t152 + (t79 * mrSges(5,1) + t868) * t428 + (Ifges(5,5) * t318 - Ifges(5,6) * t317) * t769 + t847 * t123 + t839 * t153 + (-t154 / 0.2e1 + t220 / 0.2e1 - Ifges(4,4) * t766 + t876) * t379 + (Ifges(5,1) * t318 - Ifges(5,4) * t317) * t772 + ((t782 + t703 / 0.2e1 + t324 / 0.2e1 + t708 / 0.2e1 + t894) * t494 + (t783 + t790 - t702 / 0.2e1 - t710 / 0.2e1 - t707 / 0.2e1 + t711 / 0.2e1 + t712 / 0.2e1 + t713 / 0.2e1 + t893) * t490) * t659 + (-t221 / 0.2e1 + Ifges(4,4) * t767 + t881) * t380 + t832 * t213 + t833 * t212 + t834 * t278 + (-pkin(2) * t188 * t484 + g(1) * t424 + g(2) * t422 - g(3) * t664 + t193 * t835 + t194 * t834 - t262 * t284 + t435 * t85 + t437 * t84) * m(4) + t835 * t279 + t836 * t124 + (-t669 * g(1) - t670 * g(2) - t637 * g(3) + t172 * t837 + t19 * t310 + t20 * t309 + t406 * t79 + t832 * t91 + t833 * t92) * m(5) + t837 * t177 + (-mrSges(7,1) * t671 + mrSges(7,2) * t672) * t57 + (-mrSges(6,1) * t671 + mrSges(6,2) * t672) * t80 + (-mrSges(7,1) * t668 - mrSges(7,3) * t672) * t24 + (-mrSges(6,1) * t668 - mrSges(6,3) * t672) * t40 + (t365 / 0.2e1 - t318 / 0.2e1) * t156 + (-mrSges(5,1) * t668 - mrSges(5,2) * t667) * t172 + (mrSges(7,2) * t668 + mrSges(7,3) * t671) * t29 + (mrSges(6,2) * t668 + mrSges(6,3) * t671) * t41 + (t508 * mrSges(3,1) - t275 * mrSges(5,1) - t322 * mrSges(4,1) + t516 * t431 + t741 * (t275 * t492 + t321 * t488) - t866 * (-t275 * t488 + t321 * t492) + t511 * t321 + t813 * (t322 * t489 - t431 * t688)) * g(1) + (t524 * mrSges(3,1) - t273 * mrSges(5,1) - t320 * mrSges(4,1) + t516 * t430 + t741 * (t273 * t492 + t319 * t488) - t866 * (-t273 * t488 + t319 * t492) + t511 * t319 + t813 * (t320 * t489 - t430 * t688)) * g(2) + (t432 - t332 * mrSges(5,1) - t397 * mrSges(4,1) - mrSges(4,3) * t641 + t741 * (t332 * t492 + t396 * t488) - t866 * (-t332 * t488 + t396 * t492) + t511 * t396 + t813 * (t397 * t489 - t493 * t641)) * g(3) + (-t495 * t463 / 0.2e1 - t673 / 0.2e1 + t371 * t759 - t531 / 0.2e1 + (t540 + t700 * t759 - t530 / 0.2e1) * t661 + (t413 * t495 + t415 * t491) * mrSges(3,3) + (-t219 / 0.2e1 - t709 / 0.2e1 - t706 / 0.2e1 - t193 * mrSges(4,1) - t701 / 0.2e1 + t194 * mrSges(4,2)) * t689) * t661 + (mrSges(5,2) * t79 + t907) * t429 + (t336 / 0.2e1 + t225 / 0.2e1 + t224 / 0.2e1 + t905) * t486 - t326 * mrSges(3,2) + t327 * mrSges(3,1) + t309 * t88 + t310 * t89 + t297 * t23 - t284 * t251 + t247 * t22 + (Ifges(5,4) * t318 - Ifges(5,2) * t317) * t774 + t207 * t35 + t208 * t37 + t174 * t36 + (t667 * t91 + t668 * t92) * mrSges(5,3) + t160 * t34 + t406 * t72 - t413 * t412 + t415 * t411 + t435 * t179 + t437 * t180 + (-pkin(2) * t144 + (t188 * mrSges(4,2) + t906) * t490 + (-t222 / 0.2e1 - t132 / 0.2e1 - t133 / 0.2e1 - t188 * mrSges(4,1) + t918) * t494) * t484 + (Ifges(5,5) * t365 - Ifges(5,6) * t366) * t768 + (Ifges(5,1) * t365 - Ifges(5,4) * t366) * t771 + (Ifges(5,4) * t365 - Ifges(5,2) * t366) * t773 + t636; (-t825 * t653 + (t489 * t861 + t493 * t824) * qJD(4)) * t786 + (-t823 * t653 + (-t489 * t863 + t493 * t822) * qJD(4)) * t784 - t858 * t682 / 0.2e1 - t859 * t493 / 0.2e1 + (-t827 * t653 + (-t489 * t860 + t493 * t826) * qJD(4)) * t775 + t730 * t791 + t680 * t886 + (t245 * t862 + t246 * t864) * t787 + (t245 * t864 + t246 * t865) * t785 + (t245 * t861 - t246 * t863) * t776 + t843 * t153 + (t3 * t404 + t4 * t403 + (t655 * t80 + t724) * pkin(11) - t115 * t80 + t844 * t41 + t843 * t40) * m(6) + t844 * t151 + t845 * t150 + (t1 * t354 + t2 * t370 + t24 * t846 + t29 * t845 + t449 * t9 + t57 * t842) * m(7) + t846 * t152 + t842 * t123 + (-t1 * t680 - t2 * t682) * mrSges(7,3) + (t324 + t221) * t767 + t876 * t330 + t729 * t792 + (-t646 - t142) * t212 + (t270 * t557 + t271 * t562 + t325 * t552) * qJD(4) / 0.2e1 + (-t693 / 0.2e1 + t610) * t156 + (t873 - t878) * t694 + t882 * t124 + (-t696 * t801 + mrSges(4,2) * t300 + t650 * (-t299 * pkin(3) + pkin(11) * t300) + t741 * (-t299 * t675 + t696) - t866 * (t299 * t681 + t300 * t492) + t808 * t299) * g(2) + (-t695 * t801 + mrSges(4,2) * t304 + t650 * (-t303 * pkin(3) + pkin(11) * t304) + t741 * (-t303 * t675 + t695) - t866 * (t303 * t681 + t304 * t492) + t808 * t303) * g(1) + (-t692 * t801 + t571 + t650 * t601 + t741 * (-t355 * t675 + t692) - t866 * (t355 * t681 + t356 * t492) + t812 * t355) * g(3) + (t552 * t769 + t557 * t774 + t562 * t772 - t532 + t881) * t329 + t878 * t656 + (-t707 + t154) * t766 + qJD(4) * t532 + (t9 * t563 + t793 * t826 + t797 * t824 + t798 * t822 + t799 + t875) * t489 + (mrSges(7,1) * t818 + mrSges(7,3) * t830) * t24 + (mrSges(6,1) * t818 + mrSges(6,3) * t830) * t40 + (mrSges(7,1) * t831 - mrSges(7,2) * t830) * t57 + (mrSges(6,1) * t831 - mrSges(6,2) * t830) * t80 + (-mrSges(7,2) * t818 - mrSges(7,3) * t831) * t29 + (-mrSges(6,2) * t818 - mrSges(6,3) * t831) * t41 + (-t88 + t23) * pkin(11) * t489 + (-t141 * t91 - t142 * t92 - t172 * t194 - pkin(3) * t79 + ((-t489 * t92 - t493 * t91) * qJD(4) + t817) * pkin(11)) * m(5) + (-g(1) * t304 - g(2) * t300 - g(3) * t356 - t818 * t92 + (-t655 + t693) * t91 + t817) * mrSges(5,3) + (-t645 - t141) * t213 + (-t3 * t682 - t4 * t680) * mrSges(6,3) + t79 * t569 + t117 + t575 + t370 * t36 + (pkin(11) * t89 - t908) * t493 + t354 * t34 - t193 * t278 + (t279 - t177) * t194 + t403 * t35 + t404 * t37 + t449 * t22 + t566 * t724 + t220 * t765 - pkin(3) * t72 + t851 * (t488 * t605 + t492 * t610 - t246 / 0.2e1) + t852 * (t492 * t605 - t623 / 0.2e1 - t245 / 0.2e1); (-t731 + t853) * t772 + (t267 + t156) * t774 + t488 * t886 + t838 * t123 + t840 * t150 + (t1 * t455 - t2 * t456 + t24 * t841 + t29 * t840 - t482 * t9 + t57 * t838) * m(7) + t841 * t152 + (-Ifges(5,2) * t774 - Ifges(5,6) * t769 - t872 + t873) * t271 + (-pkin(4) * t12 - t40 * t60 - t41 * t61) * m(6) + (-m(6) * t80 - t124 + t213 + t738) * t92 + (-t212 + t739) * t91 + (-t1 * t488 + t2 * t492 + t24 * t820 - t29 * t896) * mrSges(7,3) + (t40 * t820 - t41 * t896 + t816) * mrSges(6,3) + (t243 * t814 + t244 * t813) * g(1) + (t239 * t814 + t240 * t813) * g(2) + (m(6) * ((-t40 * t492 - t41 * t488) * qJD(5) + t816) - t151 * t654 - t153 * t652 - t488 * t35 + t492 * t37) * pkin(12) + (t209 * t824 + t210 * t822 + t268 * t826) * qJD(5) / 0.2e1 + t823 * t798 + t9 * t564 + t12 * t567 + t492 * t885 + t825 * t797 + t827 * t793 + (-t172 * mrSges(5,2) + Ifges(5,1) * t772 + Ifges(5,5) * t769 + t776 * t826 + t785 * t822 + t787 * t824 - t828) * t270 + t828 * qJD(5) + (t293 * t829 + t294 * t821 + t570) * g(3) - t574 + t47 - pkin(4) * t23 + (-t654 / 0.2e1 + t698 / 0.2e1) * t852 + (t652 / 0.2e1 - t697 / 0.2e1) * t851 - t61 * t151 - t60 * t153 + t455 * t34 - t456 * t36 - t482 * t22 + t155 * t771; (-t866 * (-t240 * t492 - t299 * t488) + t815 * (-t240 * t488 + t299 * t492)) * g(2) + t874 - t123 * t753 + (-m(7) * t753 - mrSges(7,1) * t210 - mrSges(7,2) * t209) * t57 + (-t151 + t737) * t40 + (t209 * t865 - t895) * t785 + (-t209 * t863 - t210 * t861) * t776 + (t153 + t736) * t41 + t852 * t784 + t848 * t1 + (-m(7) * (-t24 + t28) + t152 + t734) * t29 + (-m(7) * t752 + t565 - t568) * g(3) - t80 * (mrSges(6,1) * t210 + mrSges(6,2) * t209) - t28 * t150 + (-t210 * t862 + t851 + t897) * t787 + (t170 * t815 + t171 * t866) * g(1) + t24 * t735 + pkin(5) * t34 + t859; -t209 * t150 + t210 * t152 + (-g(1) * t243 - g(2) * t239 - g(3) * t293 - t29 * t209 + t24 * t210 + t9) * m(7) + t22;];
tau  = t10;
