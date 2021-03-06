% Calculate vector of inverse dynamics joint torques for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR6_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:36
% EndTime: 2019-12-31 19:33:03
% DurationCPUTime: 76.74s
% Computational Cost: add. (26924->1239), mult. (32907->1566), div. (0->0), fcn. (29893->10), ass. (0->582)
t997 = Icges(3,3) + Icges(4,3);
t519 = qJ(2) + pkin(8);
t496 = sin(t519);
t498 = cos(t519);
t526 = sin(qJ(2));
t528 = cos(qJ(2));
t984 = Icges(3,5) * t528 + Icges(4,5) * t498 - Icges(3,6) * t526 - Icges(4,6) * t496;
t529 = cos(qJ(1));
t996 = t997 * t529;
t527 = sin(qJ(1));
t801 = t527 * t528;
t803 = t526 * t527;
t808 = t498 * t527;
t811 = t496 * t527;
t985 = -Icges(3,5) * t801 - Icges(4,5) * t808 + Icges(3,6) * t803 + Icges(4,6) * t811 + t996;
t992 = t997 * t527 + t984 * t529;
t836 = Icges(4,6) * t529;
t302 = Icges(4,4) * t808 - Icges(4,2) * t811 - t836;
t837 = Icges(3,6) * t529;
t347 = Icges(3,4) * t801 - Icges(3,2) * t803 - t837;
t995 = t302 * t496 + t347 * t526;
t847 = Icges(4,4) * t496;
t392 = Icges(4,1) * t498 - t847;
t305 = Icges(4,5) * t527 + t392 * t529;
t848 = Icges(3,4) * t526;
t431 = Icges(3,1) * t528 - t848;
t350 = Icges(3,5) * t527 + t431 * t529;
t994 = -t305 * t808 - t350 * t801;
t993 = Icges(3,5) * t526 + Icges(4,5) * t496 + Icges(3,6) * t528 + Icges(4,6) * t498;
t448 = Icges(4,4) * t811;
t842 = Icges(4,5) * t529;
t304 = Icges(4,1) * t808 - t448 - t842;
t474 = Icges(3,4) * t803;
t843 = Icges(3,5) * t529;
t349 = Icges(3,1) * t801 - t474 - t843;
t979 = -t304 * t498 - t349 * t528 + t995;
t389 = Icges(4,2) * t498 + t847;
t479 = Icges(4,4) * t498;
t391 = Icges(4,1) * t496 + t479;
t428 = Icges(3,2) * t528 + t848;
t509 = Icges(3,4) * t528;
t430 = Icges(3,1) * t526 + t509;
t982 = t389 * t496 - t391 * t498 + t428 * t526 - t430 * t528;
t991 = -t529 * t992 - t994;
t800 = t528 * t529;
t807 = t498 * t529;
t949 = -t305 * t807 - t350 * t800 - t527 * t992;
t990 = -t304 * t807 - t349 * t800 + t527 * t985;
t632 = -Icges(4,2) * t496 + t479;
t303 = Icges(4,6) * t527 + t529 * t632;
t633 = -Icges(3,2) * t526 + t509;
t348 = Icges(3,6) * t527 + t529 * t633;
t989 = t303 * t496 + t348 * t526;
t522 = sin(pkin(9));
t797 = t529 * t522;
t523 = cos(pkin(9));
t805 = t523 * t527;
t369 = -t498 * t797 + t805;
t804 = t523 * t529;
t806 = t522 * t527;
t370 = t498 * t804 + t806;
t810 = t496 * t529;
t175 = Icges(5,5) * t370 + Icges(5,6) * t369 + Icges(5,3) * t810;
t178 = Icges(5,4) * t370 + Icges(5,2) * t369 + Icges(5,6) * t810;
t181 = Icges(5,1) * t370 + Icges(5,4) * t369 + Icges(5,5) * t810;
t367 = t498 * t806 + t804;
t368 = t498 * t805 - t797;
t59 = t175 * t811 - t367 * t178 + t368 * t181;
t956 = -t303 * t811 - t348 * t803 + t991;
t939 = t59 + t956;
t61 = t175 * t810 + t369 * t178 + t370 * t181;
t802 = t526 * t529;
t955 = -t303 * t810 - t348 * t802 - t949;
t937 = t61 + t955;
t988 = t302 * t810 + t347 * t802 + t990;
t987 = t302 * t498 + t304 * t496 + t347 * t528 + t349 * t526;
t986 = -t303 * t498 - t305 * t496 - t348 * t528 - t350 * t526;
t983 = t993 * qJD(2);
t981 = -t389 * t498 - t391 * t496 - t428 * t528 - t430 * t526;
t980 = -t305 * t498 - t350 * t528 + t989;
t948 = t993 * t529;
t947 = t993 * t527;
t173 = Icges(5,5) * t368 - Icges(5,6) * t367 + Icges(5,3) * t811;
t176 = Icges(5,4) * t368 - Icges(5,2) * t367 + Icges(5,6) * t811;
t179 = Icges(5,1) * t368 - Icges(5,4) * t367 + Icges(5,5) * t811;
t620 = -t176 * t367 + t179 * t368;
t940 = t173 * t811 - t527 * t979 + t529 * t985 + t620;
t978 = t992 * qJD(1);
t629 = Icges(5,5) * t523 - Icges(5,6) * t522;
t283 = -Icges(5,3) * t498 + t496 * t629;
t631 = Icges(5,4) * t523 - Icges(5,2) * t522;
t285 = -Icges(5,6) * t498 + t496 * t631;
t635 = Icges(5,1) * t523 - Icges(5,4) * t522;
t287 = -Icges(5,5) * t498 + t496 * t635;
t936 = t283 * t811 - t285 * t367 + t287 * t368 - t527 * t982 - t948;
t935 = t283 * t810 + t285 * t369 + t287 * t370 - t529 * t982 + t947;
t518 = pkin(9) + qJ(5);
t495 = sin(t518);
t497 = cos(t518);
t331 = t495 * t808 + t497 * t529;
t798 = t529 * t495;
t332 = t497 * t808 - t798;
t152 = Icges(6,5) * t332 - Icges(6,6) * t331 + Icges(6,3) * t811;
t317 = Icges(6,4) * t332;
t155 = -Icges(6,2) * t331 + Icges(6,6) * t811 + t317;
t316 = Icges(6,4) * t331;
t159 = -Icges(6,1) * t332 - Icges(6,5) * t811 + t316;
t963 = t155 * t495 + t159 * t497;
t56 = -t152 * t498 - t496 * t963;
t373 = t632 * qJD(2);
t374 = t392 * qJD(2);
t403 = t633 * qJD(2);
t404 = t431 * qJD(2);
t977 = qJD(1) * t993 + t981 * qJD(2) - t373 * t496 + t374 * t498 - t403 * t526 + t404 * t528;
t575 = qJD(2) * t389;
t184 = -t529 * t575 + (-t527 * t632 + t836) * qJD(1);
t577 = qJD(2) * t391;
t186 = -t529 * t577 + (-t392 * t527 + t842) * qJD(1);
t576 = qJD(2) * t428;
t222 = -t529 * t576 + (-t527 * t633 + t837) * qJD(1);
t578 = qJD(2) * t430;
t224 = -t529 * t578 + (-t431 * t527 + t843) * qJD(1);
t976 = t986 * qJD(2) - t184 * t496 + t186 * t498 - t222 * t526 + t224 * t528 + t978;
t185 = qJD(1) * t303 - t527 * t575;
t187 = qJD(1) * t305 - t527 * t577;
t223 = qJD(1) * t348 - t527 * t576;
t225 = qJD(1) * t350 - t527 * t578;
t975 = t985 * qJD(1) + t987 * qJD(2) + t185 * t496 - t187 * t498 + t223 * t526 - t225 * t528;
t974 = t982 * qJD(1) + t984 * qJD(2);
t973 = t979 * qJD(1) - t983 * t527 + t978;
t972 = -t983 * t529 + (-t984 * t527 + t980 + t996) * qJD(1);
t927 = -t369 * t176 - t370 * t179;
t60 = t173 * t810 - t927;
t855 = t529 * t60;
t971 = t937 * t527 + t988 * t529 - t855;
t970 = t939 * t527 - t940 * t529;
t917 = -t332 * rSges(6,1) + t331 * rSges(6,2);
t163 = rSges(6,3) * t811 - t917;
t863 = rSges(6,2) * t495;
t866 = rSges(6,1) * t497;
t655 = -t863 + t866;
t270 = -rSges(6,3) * t498 + t496 * t655;
t737 = qJD(5) * t496;
t741 = qJD(2) * t529;
t384 = -t527 * t737 + t741;
t736 = qJD(5) * t498;
t460 = qJD(1) - t736;
t525 = -pkin(7) - qJ(4);
t796 = qJ(4) + t525;
t921 = t498 * t796;
t485 = pkin(4) * t523 + pkin(3);
t870 = pkin(3) - t485;
t922 = t496 * t870;
t250 = t921 - t922;
t397 = pkin(3) * t496 - qJ(4) * t498;
t879 = pkin(2) * t526;
t683 = -t397 - t879;
t667 = -t250 + t683;
t969 = -t163 * t460 - t270 * t384 + t667 * t741;
t398 = rSges(4,1) * t496 + rSges(4,2) * t498;
t358 = t398 * t527;
t510 = t527 * rSges(4,3);
t482 = t498 * rSges(4,1);
t915 = -rSges(4,2) * t496 + t482;
t192 = -qJD(2) * t358 + (t529 * t915 + t510) * qJD(1);
t375 = t915 * qJD(2);
t733 = qJD(1) * qJD(2);
t413 = -qJDD(2) * t529 + t527 * t733;
t732 = qJD(1) * qJD(3);
t706 = qJDD(3) * t527 + t413 * t879 + t529 * t732;
t325 = rSges(4,1) * t808 - rSges(4,2) * t811 - t529 * rSges(4,3);
t516 = t529 * pkin(6);
t457 = pkin(1) * t527 - t516;
t524 = -qJ(3) - pkin(6);
t490 = t529 * t524;
t515 = t528 * pkin(2);
t491 = t515 + pkin(1);
t756 = -t527 * t491 - t490;
t298 = t457 + t756;
t783 = t298 - t457;
t711 = -t325 + t783;
t799 = t528 * qJD(2) ^ 2;
t725 = pkin(2) * t799;
t514 = t527 * pkin(6);
t746 = qJD(1) * t527;
t728 = pkin(2) * t803;
t754 = qJD(2) * t728 + qJD(3) * t529;
t705 = t524 * t746 + t754;
t871 = pkin(1) - t491;
t219 = (-t529 * t871 - t514) * qJD(1) - t705;
t458 = t529 * pkin(1) + t514;
t411 = t458 * qJD(1);
t792 = -t219 - t411;
t47 = t398 * t413 + (-qJD(2) * t375 - t725) * t529 + t711 * qJDD(1) + (-t192 + t792) * qJD(1) + t706;
t968 = t47 - g(1);
t480 = t496 * rSges(6,3);
t271 = t498 * t655 + t480;
t481 = t496 * rSges(5,3);
t864 = rSges(5,2) * t522;
t867 = rSges(5,1) * t523;
t658 = -t864 + t867;
t967 = t498 * t658 + t481;
t659 = t368 * rSges(5,1) - t367 * rSges(5,2);
t966 = -t659 + t756;
t189 = -rSges(5,3) * t811 - t659;
t965 = t935 * qJD(1);
t964 = t936 * qJD(1);
t742 = qJD(2) * t527;
t383 = t529 * t737 + t742;
t49 = t152 * t811 - t155 * t331 - t159 * t332;
t333 = t497 * t527 - t498 * t798;
t334 = t495 * t527 + t497 * t807;
t154 = Icges(6,5) * t334 + Icges(6,6) * t333 + Icges(6,3) * t810;
t846 = Icges(6,4) * t334;
t157 = Icges(6,2) * t333 + Icges(6,6) * t810 + t846;
t318 = Icges(6,4) * t333;
t160 = Icges(6,1) * t334 + Icges(6,5) * t810 + t318;
t50 = t154 * t811 - t331 * t157 + t332 * t160;
t628 = Icges(6,5) * t497 - Icges(6,6) * t495;
t263 = -Icges(6,3) * t498 + t496 * t628;
t844 = Icges(6,4) * t497;
t630 = -Icges(6,2) * t495 + t844;
t265 = -Icges(6,6) * t498 + t496 * t630;
t845 = Icges(6,4) * t495;
t634 = Icges(6,1) * t497 - t845;
t267 = -Icges(6,5) * t498 + t496 * t634;
t81 = t263 * t811 - t265 * t331 + t267 * t332;
t13 = t383 * t50 - t384 * t49 + t460 * t81;
t51 = t152 * t810 + t333 * t155 - t159 * t334;
t52 = t154 * t810 + t333 * t157 + t334 * t160;
t82 = t263 * t810 + t265 * t333 + t267 * t334;
t14 = t383 * t52 - t384 * t51 + t82 * t460;
t700 = t496 * t741;
t569 = -t498 * t746 - t700;
t704 = t496 * t746;
t699 = t498 * t741;
t421 = qJ(4) * t699;
t738 = qJD(4) * t529;
t437 = t496 * t738;
t762 = t421 + t437;
t169 = pkin(3) * t569 - qJ(4) * t704 + t762;
t745 = qJD(1) * t529;
t339 = t496 * t745 + t498 * t742;
t701 = t496 * t742;
t425 = pkin(3) * t701;
t461 = pkin(3) * t807;
t739 = qJD(4) * t527;
t694 = t496 * t739;
t170 = qJ(4) * t339 + qJD(1) * t461 - t425 + t694;
t443 = qJ(4) * t808;
t357 = -pkin(3) * t811 + t443;
t478 = t496 * qJ(4);
t483 = t498 * pkin(3);
t914 = t483 + t478;
t359 = t914 * t527;
t445 = qJ(4) * t807;
t360 = -pkin(3) * t810 + t445;
t477 = qJD(4) * t496;
t494 = pkin(6) * t745;
t501 = qJD(3) * t527;
t698 = t526 * t741;
t599 = -pkin(2) * t698 + t501;
t218 = -t494 + (t527 * t871 - t490) * qJD(1) + t599;
t713 = t529 * t218 + t527 * t219 - t298 * t745;
t953 = t529 * t169 + t527 * t170 - t357 * t742 + t359 * t745 - t360 * t741 - t477 + t713;
t952 = t985 + t989;
t731 = qJD(2) * qJD(4);
t951 = qJDD(4) * t496 + t498 * t731;
t950 = t984 * qJD(1);
t946 = qJD(2) * t970 + t964;
t945 = qJD(2) * t971 + t965;
t254 = qJD(1) * t369 + t522 * t701;
t255 = qJD(1) * t370 - t523 * t701;
t103 = Icges(5,5) * t255 + Icges(5,6) * t254 + Icges(5,3) * t339;
t105 = Icges(5,4) * t255 + Icges(5,2) * t254 + Icges(5,6) * t339;
t107 = Icges(5,1) * t255 + Icges(5,4) * t254 + Icges(5,5) * t339;
t926 = t176 * t522 - t179 * t523;
t944 = -t223 * t528 - t225 * t526 + (t103 - t185) * t498 + (t105 * t522 - t107 * t523 - t187) * t496 + (-t173 * t496 + t498 * t926 + t979) * qJD(2);
t252 = qJD(1) * t367 + t522 * t700;
t253 = -qJD(1) * t368 - t523 * t700;
t340 = t699 - t704;
t102 = Icges(5,5) * t253 + Icges(5,6) * t252 + Icges(5,3) * t340;
t104 = Icges(5,4) * t253 + Icges(5,2) * t252 + Icges(5,6) * t340;
t106 = Icges(5,1) * t253 + Icges(5,4) * t252 + Icges(5,5) * t340;
t618 = -t178 * t522 + t181 * t523;
t943 = t222 * t528 + t224 * t526 + (-t102 + t184) * t498 + (-t104 * t522 + t106 * t523 + t186) * t496 + (t175 * t496 + t498 * t618 - t980) * qJD(2);
t284 = Icges(5,3) * t496 + t498 * t629;
t260 = t284 * qJD(2);
t286 = Icges(5,6) * t496 + t498 * t631;
t261 = t286 * qJD(2);
t288 = Icges(5,5) * t496 + t498 * t635;
t262 = t288 * qJD(2);
t942 = t252 * t285 + t253 * t287 + t260 * t810 + t261 * t369 + t262 * t370 + t283 * t340 + t974 * t527 + t529 * t977;
t941 = t254 * t285 + t255 * t287 + t260 * t811 - t261 * t367 + t262 * t368 + t283 * t339 + t527 * t977 - t974 * t529;
t938 = t60 - t988;
t934 = -t173 * t498 - t496 * t926 + t987;
t933 = -t175 * t498 + t496 * t618 - t986;
t747 = qJD(1) * t498;
t759 = t430 + t633;
t760 = -t428 + t431;
t769 = t391 + t632;
t770 = -t389 + t392;
t900 = -t527 * t618 - t529 * t926;
t613 = -t285 * t522 + t287 * t523;
t910 = (t284 - t613) * qJD(1);
t932 = (qJD(2) * t900 + t910) * t496 + t283 * t747 + (-t496 * t769 + t498 * t770 - t526 * t759 + t528 * t760) * qJD(1);
t559 = t347 * t529 - t348 * t527;
t560 = t302 * t529 - t303 * t527;
t898 = t527 * (-t389 * t529 + t305) - t529 * (-Icges(4,2) * t808 + t304 - t448);
t899 = t527 * (-t428 * t529 + t350) - t529 * (-Icges(3,2) * t801 + t349 - t474);
t931 = -t496 * t898 - t526 * t899 + t559 * t528 + (-t173 * t529 + t175 * t527 + t560) * t498;
t930 = (-t103 * t811 + t105 * t367 - t107 * t368 - t173 * t339 - t176 * t254 - t179 * t255 + t973 * t529) * t529 + (t102 * t811 - t104 * t367 + t106 * t368 + t175 * t339 + t178 * t254 + t181 * t255 + t976 * t527 + (-t972 + t975) * t529) * t527;
t929 = (-t103 * t810 - t105 * t369 - t107 * t370 - t173 * t340 - t176 * t252 - t179 * t253 + t975 * t529) * t529 + (t102 * t810 + t104 * t369 + t106 * t370 + t175 * t340 + t178 * t252 + t181 * t253 + t972 * t527 + (-t973 + t976) * t529) * t527;
t726 = pkin(4) * t797;
t668 = -t485 * t808 + t726;
t923 = t496 * t796;
t171 = (t483 + t923) * t527 + t668;
t710 = -t359 + t783;
t670 = t171 + t710;
t758 = t437 + t501;
t39 = qJD(1) * t670 + t758 + t969;
t165 = t334 * rSges(6,1) + t333 * rSges(6,2) + rSges(6,3) * t810;
t362 = qJ(4) * t810 + t461;
t471 = pkin(4) * t806;
t597 = t485 * t807 - t525 * t810 + t471;
t172 = t597 - t362;
t707 = -t397 * t742 - t754;
t675 = t529 * t491 - t524 * t527;
t299 = t675 - t458;
t780 = t299 + t458;
t709 = t362 + t780;
t40 = t165 * t460 - t270 * t383 + (-qJD(2) * t250 + t477) * t527 + (t172 + t709) * qJD(1) + t707;
t928 = t39 * t529 + t40 * t527;
t925 = -rSges(5,3) - qJ(4);
t920 = t498 * t870;
t571 = -t398 - t879;
t919 = t529 * t571;
t417 = qJD(1) * t457;
t918 = qJD(1) * t298 - t417;
t682 = -t914 - t515;
t916 = t915 + t515;
t876 = g(2) * t527;
t911 = g(1) * t529 + t876;
t909 = t911 * t496;
t674 = -qJD(5) + t747;
t908 = t527 * t674 + t700;
t907 = t529 * t674 - t701;
t264 = Icges(6,3) * t496 + t498 * t628;
t614 = -t265 * t495 + t267 * t497;
t622 = -t157 * t495 + t160 * t497;
t897 = t383 * (-t263 * t529 - t622) - t384 * (-t263 * t527 + t963) + t460 * (t264 - t614);
t307 = (-Icges(6,2) * t497 - t845) * t496;
t896 = t383 * (-Icges(6,2) * t334 + t160 + t318) - t384 * (-Icges(6,2) * t332 - t159 - t316) + t460 * (t267 + t307);
t895 = m(5) / 0.2e1;
t894 = m(6) / 0.2e1;
t893 = -m(5) - m(6);
t412 = qJDD(2) * t527 + t529 * t733;
t729 = qJDD(5) * t496;
t216 = qJD(5) * t340 + t529 * t729 + t412;
t892 = t216 / 0.2e1;
t217 = qJD(5) * t339 + t527 * t729 + t413;
t891 = t217 / 0.2e1;
t371 = qJD(2) * t737 - qJDD(5) * t498 + qJDD(1);
t890 = t371 / 0.2e1;
t889 = -t383 / 0.2e1;
t888 = t383 / 0.2e1;
t887 = -t384 / 0.2e1;
t886 = t384 / 0.2e1;
t885 = t412 / 0.2e1;
t884 = t413 / 0.2e1;
t883 = -t460 / 0.2e1;
t882 = t460 / 0.2e1;
t881 = t527 / 0.2e1;
t880 = -t529 / 0.2e1;
t878 = g(1) * t527;
t320 = (-rSges(6,1) * t495 - rSges(6,2) * t497) * t496;
t146 = qJD(2) * t271 + qJD(5) * t320;
t766 = qJD(1) * (-pkin(1) * t746 + t494) + qJDD(1) * t458;
t563 = qJD(1) * t218 + qJDD(1) * t299 - qJDD(3) * t529 + t527 * t732 + t766;
t543 = qJDD(1) * t362 + t563 + t951 * t527 + (t169 + t437) * qJD(1);
t561 = -t920 - t923;
t228 = t561 * qJD(2);
t740 = qJD(4) * t498;
t319 = qJD(2) * t914 - t740;
t567 = -t725 + (-t228 - t319) * qJD(2);
t600 = t529 * t460;
t137 = t495 * t908 + t497 * t600;
t138 = t495 * t600 - t497 * t908;
t717 = t138 * rSges(6,1) + t137 * rSges(6,2) + rSges(6,3) * t699;
t79 = -rSges(6,3) * t704 + t717;
t809 = t498 * t525;
t583 = -t809 + t922;
t562 = t583 * t529;
t765 = qJD(1) * t726 + t525 * t704;
t95 = -t421 + qJD(2) * t562 + (t478 + t920) * t746 + t765;
t4 = qJD(1) * t95 + qJDD(1) * t172 - t383 * t146 + t371 * t165 - t216 * t270 + t412 * t667 + t460 * t79 + t527 * t567 + t543;
t875 = t4 * t529;
t580 = -t170 - t694 + t792;
t594 = t413 * t397 + t529 * t951 + t706;
t601 = t527 * t460;
t139 = -t495 * t907 + t497 * t601;
t140 = t495 * t601 + t497 * t907;
t657 = rSges(6,1) * t140 + rSges(6,2) * t139;
t80 = rSges(6,3) * t339 + t657;
t813 = t485 * t496;
t96 = t425 + (-t813 - t921) * t742 + (t529 * t561 + t471) * qJD(1);
t5 = -t384 * t146 - t371 * t163 + t217 * t270 + t413 * t250 - t460 * t80 + t567 * t529 + t670 * qJDD(1) + (t580 - t96) * qJD(1) + t594;
t874 = t5 * t527;
t306 = (-Icges(6,5) * t495 - Icges(6,6) * t497) * t496;
t143 = qJD(2) * t264 + qJD(5) * t306;
t266 = Icges(6,6) * t496 + t498 * t630;
t144 = qJD(2) * t266 + qJD(5) * t307;
t268 = Icges(6,5) * t496 + t498 * t634;
t308 = (-Icges(6,1) * t495 - t844) * t496;
t145 = qJD(2) * t268 + qJD(5) * t308;
t29 = (qJD(2) * t614 - t143) * t498 + (qJD(2) * t263 - t144 * t495 + t145 * t497 + (-t265 * t497 - t267 * t495) * qJD(5)) * t496;
t93 = -t263 * t498 + t496 * t614;
t869 = t29 * t460 + t93 * t371;
t868 = rSges(3,1) * t528;
t74 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t339;
t76 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t339;
t78 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t339;
t10 = (-qJD(2) * t963 - t74) * t498 + (qJD(2) * t152 - t495 * t76 + t497 * t78 + (-t155 * t497 + t159 * t495) * qJD(5)) * t496;
t862 = t10 * t384;
t73 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t340;
t75 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t340;
t77 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t340;
t11 = (qJD(2) * t622 - t73) * t498 + (qJD(2) * t154 - t495 * t75 + t497 * t77 + (-t157 * t497 - t160 * t495) * qJD(5)) * t496;
t861 = t11 * t383;
t858 = t39 * t527;
t857 = t40 * t529;
t511 = t527 * rSges(3,3);
t854 = t56 * t217;
t57 = -t154 * t498 + t496 * t622;
t853 = t57 * t216;
t851 = -rSges(6,3) + t525;
t568 = qJD(2) * t919 + t501;
t100 = qJD(1) * t711 + t568;
t828 = t100 * t398;
t439 = rSges(3,1) * t526 + rSges(3,2) * t528;
t702 = t439 * t741;
t753 = rSges(3,2) * t803 + t529 * rSges(3,3);
t363 = rSges(3,1) * t801 - t753;
t772 = -t363 - t457;
t205 = qJD(1) * t772 - t702;
t827 = t205 * t527;
t826 = t205 * t529;
t364 = rSges(3,1) * t800 - rSges(3,2) * t802 + t511;
t280 = t364 + t458;
t206 = qJD(1) * t280 - t439 * t742;
t386 = t439 * t529;
t825 = t206 * t386;
t812 = t496 * t525;
t414 = t498 * t485;
t793 = t163 - t171;
t787 = -t298 * t742 + t299 * t741;
t786 = -t527 * t298 + t529 * t299;
t326 = rSges(4,1) * t807 - rSges(4,2) * t810 + t510;
t782 = -t299 - t326;
t781 = -t299 - t362;
t776 = qJD(1) * t360 + t498 * t739;
t703 = t526 * t746;
t468 = pkin(2) * t703;
t771 = t397 * t746 + t468;
t721 = t496 * t863;
t768 = rSges(6,3) * t808 + t527 * t721;
t767 = rSges(6,3) * t807 + t529 * t721;
t722 = t496 * t864;
t764 = rSges(5,3) * t808 + t527 * t722;
t763 = rSges(5,3) * t807 + t529 * t722;
t761 = rSges(4,2) * t704 + rSges(4,3) * t745;
t757 = t498 * t738 + t468;
t755 = rSges(3,2) * t703 + rSges(3,3) * t745;
t752 = t527 ^ 2 + t529 ^ 2;
t744 = qJD(2) * t496;
t727 = pkin(2) * t802;
t724 = t496 * t867;
t723 = t496 * t866;
t718 = qJD(2) * t515;
t716 = -t172 + t781;
t190 = t370 * rSges(5,1) + t369 * rSges(5,2) + rSges(5,3) * t810;
t715 = -t190 + t781;
t714 = t218 * t741 + t219 * t742 - t412 * t298;
t712 = t253 * rSges(5,1) + t252 * rSges(5,2) + rSges(5,3) * t699;
t697 = t528 * t741;
t693 = -pkin(1) - t868;
t689 = t745 / 0.2e1;
t688 = -t742 / 0.2e1;
t687 = t742 / 0.2e1;
t686 = -t741 / 0.2e1;
t685 = t741 / 0.2e1;
t676 = t414 - t812;
t673 = t189 + t710;
t671 = t527 * t359 + t529 * t362 + t786;
t669 = t752 * t879;
t289 = -rSges(5,3) * t498 + t496 * t658;
t665 = -t289 + t683;
t663 = -t319 - t718;
t442 = rSges(2,1) * t529 - rSges(2,2) * t527;
t440 = rSges(2,1) * t527 + rSges(2,2) * t529;
t441 = -rSges(3,2) * t526 + t868;
t660 = rSges(5,1) * t255 + rSges(5,2) * t254;
t654 = -t809 - t879;
t646 = t49 * t529 - t50 * t527;
t645 = t49 * t527 + t50 * t529;
t644 = t51 * t529 - t52 * t527;
t643 = t51 * t527 + t52 * t529;
t642 = t527 * t57 - t529 * t56;
t641 = t527 * t56 + t529 * t57;
t638 = -qJD(1) * t359 + t758 + t918;
t621 = t163 * t529 - t165 * t527;
t191 = rSges(4,1) * t569 - rSges(4,2) * t699 + t761;
t617 = t191 * t529 + t192 * t527;
t616 = -t206 * t527 - t826;
t226 = -rSges(3,2) * t697 + (-t528 * t746 - t698) * rSges(3,1) + t755;
t385 = t439 * t527;
t227 = -qJD(2) * t385 + (t441 * t529 + t511) * qJD(1);
t615 = t226 * t529 + t227 * t527;
t610 = t325 * t527 + t326 * t529;
t607 = t363 * t527 + t364 * t529;
t602 = -t270 + t667;
t269 = t967 * qJD(2);
t598 = -t269 + t663;
t67 = (-qJD(2) * t289 + t477) * t527 + (t190 + t709) * qJD(1) + t707;
t596 = t67 * t665;
t593 = -t414 - t491 - t480;
t584 = t496 * t925 - t483;
t582 = -t146 - t228 + t663;
t579 = t359 * t742 + t362 * t741 - t740 + t787;
t570 = -t914 - t481;
t566 = -t725 + (-t269 - t319) * qJD(2);
t565 = -t152 * t384 + t154 * t383 + t263 * t460;
t564 = (-Icges(6,5) * t331 - Icges(6,6) * t332) * t384 - (Icges(6,5) * t333 - Icges(6,6) * t334) * t383 - t306 * t460;
t557 = t496 * t564;
t550 = -qJDD(4) * t498 + t169 * t741 + t170 * t742 + t412 * t359 + t496 * t731 + t714;
t545 = (Icges(6,1) * t333 - t157 - t846) * t383 - (-Icges(6,1) * t331 - t155 - t317) * t384 + (-t265 + t308) * t460;
t30 = t163 * t383 + t165 * t384 + (-t171 * t527 + t172 * t529) * qJD(2) + t579;
t542 = t30 * t621 + (-t857 + t858) * t270;
t533 = t897 * t496;
t406 = t441 * qJD(2);
t361 = t398 * t529;
t338 = t752 * t744;
t249 = -t529 * t724 + t763;
t248 = -t527 * t724 + t764;
t247 = t287 * t529;
t246 = t287 * t527;
t245 = t285 * t529;
t244 = t285 * t527;
t238 = -t529 * t723 + t767;
t237 = -t527 * t723 + t768;
t236 = t267 * t529;
t235 = t267 * t527;
t234 = t265 * t529;
t233 = t265 * t527;
t215 = -t445 + t562;
t214 = t527 * t583 - t443;
t201 = rSges(6,1) * t333 - rSges(6,2) * t334;
t200 = -rSges(6,1) * t331 - rSges(6,2) * t332;
t199 = t607 * qJD(2);
t113 = rSges(5,3) * t339 + t660;
t112 = -rSges(5,3) * t704 + t712;
t101 = -t398 * t742 + (t326 + t780) * qJD(1) - t754;
t94 = qJD(2) * t610 + t787;
t90 = qJD(1) * t226 + qJDD(1) * t364 - t406 * t742 - t412 * t439 + t766;
t89 = -t406 * t741 + t413 * t439 + t772 * qJDD(1) + (-t227 - t411) * qJD(1);
t66 = qJD(1) * t673 + t665 * t741 + t758;
t53 = (-t189 * t527 + t190 * t529) * qJD(2) + t579;
t48 = -t375 * t742 + qJD(1) * t191 + qJDD(1) * t326 - t398 * t412 + (-t412 * t526 - t527 * t799) * pkin(2) + t563;
t25 = t139 * t265 + t140 * t267 + t143 * t811 - t144 * t331 + t145 * t332 + t263 * t339;
t24 = t137 * t265 + t138 * t267 + t143 * t810 + t144 * t333 + t145 * t334 + t263 * t340;
t23 = qJD(1) * t112 + qJDD(1) * t190 + t412 * t665 + t527 * t566 + t543;
t22 = t289 * t413 + t566 * t529 + t673 * qJDD(1) + (-t113 + t580) * qJD(1) + t594;
t19 = t383 * t57 - t384 * t56 + t460 * t93;
t12 = -t189 * t412 + (t112 * t529 + t113 * t527) * qJD(2) + t715 * t413 + t550;
t9 = t139 * t157 + t140 * t160 + t154 * t339 - t331 * t75 + t332 * t77 + t73 * t811;
t8 = t139 * t155 - t140 * t159 + t152 * t339 - t331 * t76 + t332 * t78 + t74 * t811;
t7 = t137 * t157 + t138 * t160 + t154 * t340 + t333 * t75 + t334 * t77 + t73 * t810;
t6 = t137 * t155 - t138 * t159 + t152 * t340 + t333 * t76 + t334 * t78 + t74 * t810;
t3 = t163 * t216 - t165 * t217 - t171 * t412 + t383 * t80 + t384 * t79 + (t527 * t96 + t529 * t95) * qJD(2) + t716 * t413 + t550;
t2 = t216 * t50 + t217 * t49 + t25 * t460 + t371 * t81 + t383 * t9 - t384 * t8;
t1 = t216 * t52 + t217 * t51 + t24 * t460 + t371 * t82 + t383 * t7 - t384 * t6;
t15 = [(t933 + t935) * t885 + (((t529 * t952 + t620 + t949 + t955) * t529 + (t527 * t952 - t59 + t927 + t938 - t991) * t527) * qJD(2) + t946 - t964) * t688 + (-t596 * t741 - g(1) * t966 - t584 * t878 + (t570 * t527 + t966) * t22 + (t425 - t660 + t705 + (qJD(2) * t498 * t925 - t477) * t527 + (-t491 + t570) * t745) * t66 + (-pkin(3) * t700 + t599 - t638 + t66 + t712 + t762 + ((-t491 + t584) * t527 - t490 - t189) * qJD(1)) * t67 + (t23 - g(2)) * (t675 + t190 + t362)) * m(5) + (-(qJD(1) * t171 - t39 + t638 + t969) * t40 + t39 * (-t657 - t694 + t705) + t40 * (t717 + t758 + t765) + ((t654 - t813) * t857 + (t498 * t851 + t813) * t858) * qJD(2) + ((-t39 * pkin(4) * t522 + t40 * t593) * t527 + (t39 * (t593 + t812) - t40 * t524) * t529) * qJD(1) + (t4 - g(2)) * (t597 + t675 + t165) + (t5 - g(1)) * (t811 * t851 + t668 + t756 + t917)) * m(6) + t14 * t886 + (t403 * t528 + t404 * t526 + (t373 - t260) * t498 + (-t261 * t522 + t262 * t523 + t374) * t496 + (t283 * t496 + t498 * t613 - t982) * qJD(2)) * qJD(1) + (t934 + t936) * t884 + (t942 + t943) * t687 + (t941 - t944 + t945) * t686 + (-(-qJD(1) * t325 - t100 + t568 + t918) * t101 + t100 * t705 + t101 * (t501 + t761) + (t101 * t919 + t527 * t828) * qJD(2) + ((-t100 * rSges(4,3) + t101 * (-t491 - t482)) * t527 + (t100 * (-t491 - t915) - t101 * t524) * t529) * qJD(1) + (t48 - g(2)) * (t326 + t675) + t968 * (-t325 + t756)) * m(4) + (-t283 * t498 + t496 * t613 + m(2) * (t440 ^ 2 + t442 ^ 2) + Icges(2,3) - t981) * qJDD(1) + (t25 + t14) * t887 + t869 + t861 / 0.2e1 - t862 / 0.2e1 + t853 / 0.2e1 + t854 / 0.2e1 + ((-t855 + ((t992 + t995) * t529 + t956 + t990 + t994) * t529 + (t61 - t949) * t527) * qJD(2) + t965) * t685 - m(2) * (-g(1) * t440 + g(2) * t442) + t24 * t888 + t81 * t891 + t82 * t892 + (t206 * (t494 + t755) + (t439 * t827 - t825) * qJD(2) + ((-pkin(1) - t441) * t826 + (t205 * (-rSges(3,3) - pkin(6)) + t206 * t693) * t527) * qJD(1) - (-qJD(1) * t363 - t205 - t417 - t702) * t206 + (t90 - g(2)) * t280 + (t89 - g(1)) * (t693 * t527 + t516 + t753)) * m(3); (t527 * t933 - t529 * t934) * qJDD(1) / 0.2e1 + (-g(3) * t916 - t571 * t876 + t100 * (-pkin(2) * t697 - t375 * t529) + (qJD(2) * t617 + t325 * t412 + t413 * t782 + t714) * (t610 + t786) + t94 * (t617 + t713) + (t101 * t571 + t94 * t325) * t745 + (t48 * t571 + t101 * (-t375 - t718) + (t782 * t94 + t828) * qJD(1)) * t527 - (t100 * t358 + t101 * (-t361 - t727)) * qJD(1) - (-t94 * t669 + (-t100 * t916 - t94 * t361) * t529 + (-t101 * t916 - t94 * t358) * t527) * qJD(2) + t968 * t919) * m(4) + ((-t245 * t369 - t247 * t370) * t742 + (t286 * t369 + t288 * t370) * qJD(1) + (-t742 * t948 + t950) * t527 + ((t369 * t244 + t370 * t246 + t527 * t947 + t931) * qJD(2) + t932) * t529) * t688 + (-(t244 * t367 - t246 * t368) * t741 + (-t286 * t367 + t288 * t368) * qJD(1) + (-t741 * t947 - t950) * t529 + ((t367 * t245 - t368 * t247 + t529 * t948 + t931) * qJD(2) + t932) * t527) * t685 + (g(1) * t386 + g(2) * t385 - g(3) * t441 - (t205 * t385 - t825) * qJD(1) - (t199 * (-t385 * t527 - t386 * t529) + t616 * t441) * qJD(2) + (qJD(2) * t615 + t363 * t412 - t364 * t413) * t607 + t199 * ((t363 * t529 - t364 * t527) * qJD(1) + t615) + t616 * t406 + (-t90 * t527 - t89 * t529 + (-t206 * t529 + t827) * qJD(1)) * t439) * m(3) - t217 * t646 / 0.2e1 - t216 * t644 / 0.2e1 + ((t527 * t938 + t529 * t937) * qJD(1) + t929) * t687 + ((t527 * t940 + t529 * t939) * qJD(1) + t930) * t686 + (qJD(1) * t941 + qJD(2) * t930 + qJDD(1) * t936 + t412 * t939 + t413 * t940 + t2) * t880 + (qJD(1) * t942 + qJD(2) * t929 + qJDD(1) * t935 + t412 * t937 + t413 * t938 + t1) * t881 + (t944 * t529 + t943 * t527 + (t527 * t934 + t529 * t933) * qJD(1)) * qJD(1) / 0.2e1 + (t14 + t945) * t689 + (t13 + t946) * t746 / 0.2e1 + (((t234 * t495 - t236 * t497 + t154) * t383 - (t233 * t495 - t235 * t497 + t152) * t384 + (-t266 * t495 + t268 * t497 + t263) * t460 + t93 * qJD(5)) * t496 + (qJD(5) * t641 - t897) * t498) * t883 + ((t234 * t331 - t236 * t332) * t383 - (t233 * t331 - t235 * t332) * t384 + (-t266 * t331 + t268 * t332) * t460 + (t496 * t81 + t50 * t807) * qJD(5) + ((qJD(5) * t49 + t565) * t498 + t533) * t527) * t886 + ((-t234 * t333 - t236 * t334) * t383 - (-t233 * t333 - t235 * t334) * t384 + (t266 * t333 + t268 * t334) * t460 + (t496 * t82 + t51 * t808) * qJD(5) + ((qJD(5) * t52 + t565) * t498 + t533) * t529) * t889 - (-t910 * t498 + (t559 * t526 + t528 * t899 + (t898 - t900) * t498 + ((t245 * t522 - t247 * t523 + t175) * t527 - (t244 * t522 - t246 * t523 + t173) * t529 + t560) * t496) * qJD(2) + (t498 * t769 + t526 * t760 + t528 * t759 + (-t286 * t522 + t288 * t523 + t283 + t770) * t496) * qJD(1)) * qJD(1) / 0.2e1 - (t13 * t527 + t14 * t529) * t736 / 0.2e1 + t970 * t884 + t971 * t885 - t19 * t737 / 0.2e1 + (t12 * t671 + (qJD(1) * t596 + t12 * t190 + t22 * t665) * t529 + (-t12 * t189 + t23 * t665) * t527 - g(1) * (t445 - t727 + t763) - g(2) * (t443 - t728 + t764) - (-pkin(3) - t867) * t909 + (t598 * t527 - t776 - (t249 - t727) * qJD(1)) * t67 + (t598 * t529 - t757 + t771 + (t289 * t527 + t248 + t357) * qJD(1)) * t66 + (g(3) - (t67 * t527 + t66 * t529) * qJD(2)) * (t682 - t967) + ((-qJD(1) * t189 + t112) * t529 + (qJD(1) * t715 + t113) * t527 - (t248 * t527 + t249 * t529 - t669) * qJD(2) + t953) * t53) * m(5) + (qJD(1) * t641 - t10 * t529 + t11 * t527) * t882 + (qJD(1) * t645 + t527 * t9 - t529 * t8) * t887 + (qJD(1) * t643 + t527 * t7 - t529 * t6) * t888 + t642 * t890 + (-g(1) * t767 - g(2) * t768 - g(3) * (t271 + t515 + t676) - t911 * ((-t485 - t866) * t496 + t654) + t39 * t771 + t3 * t671 + (t39 * t582 + t3 * (t165 + t172) + (t40 * qJD(1) + t5) * t602) * t529 + (t4 * t602 + t40 * t582 + t3 * t793 + t39 * (t250 + t270) * qJD(1)) * t527 - t39 * (-t237 * t460 - t271 * t384 + t757) - t40 * (t238 * t460 - t271 * t383 + t776) - (t39 * (-t214 - t357) + t40 * (t215 - t727)) * qJD(1) - ((-t163 * t39 + t165 * t40) * t496 + t542 * t498) * qJD(5) - t928 * qJD(2) * (-t676 + t914 + t682) + ((qJD(1) * t793 + t79 + t95) * t529 + (t80 + t96 + (-t165 + t716) * qJD(1)) * t527 - t237 * t383 - t238 * t384 - (t214 * t527 + t215 * t529 - t669) * qJD(2) + t953) * t30) * m(6); (-m(4) + t893) * (-g(2) * t529 + t878) + 0.2e1 * (-t875 / 0.2e1 + t874 / 0.2e1) * m(6) + 0.2e1 * (t22 * t881 + t23 * t880) * m(5) + 0.2e1 * (t47 * t881 + t48 * t880) * m(4); t893 * (-g(3) * t498 + t909) - m(5) * (t338 * t53 + t339 * t67 + t340 * t66) - m(6) * (t30 * t338 + t339 * t40 + t340 * t39) + 0.2e1 * ((t66 * t741 + t67 * t742 - t12) * t895 + (t39 * t741 + t40 * t742 - t3) * t894) * t498 + 0.2e1 * ((qJD(2) * t53 + t22 * t529 + t23 * t527 - t66 * t746 + t67 * t745) * t895 + (qJD(2) * t30 - t39 * t746 + t4 * t527 + t40 * t745 + t5 * t529) * t894) * t496; t1 * t810 / 0.2e1 + (t496 * t643 - t498 * t82) * t892 + ((qJD(2) * t643 - t24) * t498 + (qJD(1) * t644 + qJD(2) * t82 + t527 * t6 + t529 * t7) * t496) * t888 + t2 * t811 / 0.2e1 + (t496 * t645 - t498 * t81) * t891 + ((qJD(2) * t645 - t25) * t498 + (qJD(1) * t646 + qJD(2) * t81 + t527 * t8 + t529 * t9) * t496) * t887 + t19 * t744 / 0.2e1 - t498 * (t853 + t854 + t861 - t862 + t869) / 0.2e1 + (t496 * t641 - t498 * t93) * t890 + ((qJD(2) * t641 - t29) * t498 + (-qJD(1) * t642 + qJD(2) * t93 + t10 * t527 + t11 * t529) * t496) * t882 + (t333 * t896 + t545 * t334 - t529 * t557) * t889 + (-t331 * t896 + t332 * t545 - t527 * t557) * t886 + (t564 * t498 + (-t495 * t896 + t497 * t545) * t496) * t883 + (-t704 / 0.2e1 + t498 * t685) * t14 + (t496 * t689 + t498 * t687) * t13 + ((qJD(2) * t542 + t5 * t163 - t4 * t165 + t39 * t80 - t40 * t79) * t498 + (t39 * (-qJD(2) * t163 + t146 * t527) + t40 * (qJD(2) * t165 - t146 * t529) + t3 * t621 + t30 * (-t163 * t746 - t165 * t745 - t527 * t79 + t529 * t80) + (qJD(1) * t928 + t874 - t875) * t270) * t496 - t39 * (-t200 * t460 - t320 * t384) - t40 * (t201 * t460 - t320 * t383) - t30 * (t200 * t383 + t201 * t384) - g(1) * t201 - g(2) * t200 - g(3) * t320) * m(6);];
tau = t15;
