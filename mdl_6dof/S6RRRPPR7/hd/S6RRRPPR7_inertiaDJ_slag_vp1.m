% Calculate time derivative of joint inertia matrix for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR7_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:35
% EndTime: 2019-03-09 15:57:35
% DurationCPUTime: 35.92s
% Computational Cost: add. (54228->1405), mult. (120513->1858), div. (0->0), fcn. (127097->10), ass. (0->631)
t501 = sin(qJ(2));
t498 = cos(pkin(10));
t500 = sin(qJ(3));
t497 = sin(pkin(10));
t503 = cos(qJ(3));
t725 = t497 * t503;
t546 = -t498 * t500 + t725;
t426 = t546 * t501;
t726 = t497 * t500;
t545 = t498 * t503 + t726;
t427 = t545 * t501;
t504 = cos(qJ(2));
t310 = Icges(6,5) * t427 - Icges(6,6) * t426 + Icges(6,3) * t504;
t572 = Icges(4,5) * t503 - Icges(4,6) * t500;
t406 = -Icges(4,3) * t504 + t501 * t572;
t574 = Icges(5,4) * t503 + Icges(5,6) * t500;
t409 = -Icges(5,2) * t504 + t501 * t574;
t800 = -t310 + t406 + t409;
t761 = -qJD(1) / 0.2e1;
t799 = t501 * t761;
t673 = qJD(3) * t503;
t627 = t501 * t673;
t677 = qJD(2) * t504;
t635 = t500 * t677;
t524 = t627 + t635;
t674 = qJD(3) * t501;
t628 = t500 * t674;
t630 = t503 * t677;
t798 = -t628 + t630;
t502 = sin(qJ(1));
t631 = t502 * t677;
t505 = cos(qJ(1));
t680 = qJD(1) * t505;
t526 = t501 * t680 + t631;
t311 = Icges(6,4) * t427 - Icges(6,2) * t426 + Icges(6,6) * t504;
t312 = Icges(6,1) * t427 - Icges(6,4) * t426 + Icges(6,5) * t504;
t744 = Icges(5,5) * t503;
t571 = Icges(5,3) * t500 + t744;
t405 = -Icges(5,6) * t504 + t501 * t571;
t747 = Icges(4,4) * t503;
t575 = -Icges(4,2) * t500 + t747;
t410 = -Icges(4,6) * t504 + t501 * t575;
t745 = Icges(5,5) * t500;
t578 = Icges(5,1) * t503 + t745;
t413 = -Icges(5,4) * t504 + t501 * t578;
t748 = Icges(4,4) * t500;
t579 = Icges(4,1) * t503 - t748;
t414 = -Icges(4,5) * t504 + t501 * t579;
t797 = t311 * t426 - t312 * t427 + ((-t413 - t414) * t503 + (-t405 + t410) * t500) * t501 + t800 * t504;
t493 = pkin(10) + qJ(6);
t487 = sin(t493);
t488 = cos(t493);
t548 = t487 * t503 - t488 * t500;
t397 = t548 * t501;
t547 = t487 * t500 + t488 * t503;
t398 = t547 * t501;
t289 = Icges(7,5) * t398 - Icges(7,6) * t397 + Icges(7,3) * t504;
t290 = Icges(7,4) * t398 - Icges(7,2) * t397 + Icges(7,6) * t504;
t291 = Icges(7,1) * t398 - Icges(7,4) * t397 + Icges(7,5) * t504;
t718 = t504 * t505;
t442 = t500 * t718 - t502 * t503;
t443 = t502 * t500 + t503 * t718;
t322 = t442 * t488 - t443 * t487;
t323 = t442 * t487 + t443 * t488;
t721 = t501 * t505;
t136 = -t289 * t721 + t322 * t290 + t323 * t291;
t678 = qJD(2) * t502;
t634 = t501 * t678;
t675 = qJD(3) * t500;
t682 = qJD(1) * t502;
t316 = -t500 * t634 - t505 * t675 - t503 * t682 + (t500 * t680 + t502 * t673) * t504;
t672 = qJD(3) * t504;
t544 = (qJD(1) - t672) * t500;
t681 = qJD(1) * t504;
t617 = -qJD(3) + t681;
t679 = qJD(2) * t501;
t719 = t503 * t505;
t317 = t617 * t719 + (-t503 * t679 + t544) * t502;
t720 = t502 * t504;
t440 = t500 * t720 + t719;
t441 = -t500 * t505 + t503 * t720;
t321 = t440 * t487 + t441 * t488;
t188 = -qJD(6) * t321 + t316 * t488 - t317 * t487;
t320 = t440 * t488 - t441 * t487;
t189 = qJD(6) * t320 + t316 * t487 + t317 * t488;
t116 = Icges(7,5) * t189 + Icges(7,6) * t188 - Icges(7,3) * t526;
t118 = Icges(7,4) * t189 + Icges(7,2) * t188 - Icges(7,6) * t526;
t120 = Icges(7,1) * t189 + Icges(7,4) * t188 - Icges(7,5) * t526;
t626 = t503 * t672;
t676 = qJD(2) * t505;
t632 = t501 * t676;
t314 = qJD(1) * t440 + t500 * t632 - t502 * t675 - t505 * t626;
t315 = t505 * t544 + (-t502 * t617 - t632) * t503;
t186 = -qJD(6) * t323 - t314 * t488 - t315 * t487;
t187 = qJD(6) * t322 - t314 * t487 + t315 * t488;
t723 = t501 * t502;
t241 = Icges(7,5) * t321 + Icges(7,6) * t320 - Icges(7,3) * t723;
t243 = Icges(7,4) * t321 + Icges(7,2) * t320 - Icges(7,6) * t723;
t245 = Icges(7,1) * t321 + Icges(7,4) * t320 - Icges(7,5) * t723;
t629 = t504 * t676;
t640 = t501 * t682;
t525 = -t629 + t640;
t17 = -t116 * t721 + t322 * t118 + t323 * t120 + t186 * t243 + t187 * t245 + t241 * t525;
t115 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t525;
t117 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t525;
t119 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t525;
t242 = Icges(7,5) * t323 + Icges(7,6) * t322 - Icges(7,3) * t721;
t244 = Icges(7,4) * t323 + Icges(7,2) * t322 - Icges(7,6) * t721;
t246 = Icges(7,1) * t323 + Icges(7,4) * t322 - Icges(7,5) * t721;
t18 = -t115 * t721 + t322 * t117 + t323 * t119 + t186 * t244 + t187 * t246 + t242 * t525;
t783 = t501 * (qJD(3) - qJD(6));
t281 = t547 * t783 - t548 * t677;
t282 = t547 * t677 + t548 * t783;
t178 = Icges(7,5) * t282 + Icges(7,6) * t281 - Icges(7,3) * t679;
t179 = Icges(7,4) * t282 + Icges(7,2) * t281 - Icges(7,6) * t679;
t180 = Icges(7,1) * t282 + Icges(7,4) * t281 - Icges(7,5) * t679;
t37 = -t178 * t721 + t322 * t179 + t323 * t180 + t186 * t290 + t187 * t291 + t289 * t525;
t96 = -t241 * t721 + t322 * t243 + t323 * t245;
t97 = -t242 * t721 + t322 * t244 + t323 * t246;
t585 = t502 * t96 + t505 * t97;
t1 = -t37 * t504 - t97 * t640 + t585 * t677 + (qJD(2) * t136 + t17 * t502 + (qJD(1) * t96 + t18) * t505) * t501;
t135 = -t289 * t723 + t290 * t320 + t291 * t321;
t19 = -t116 * t723 + t320 * t118 + t321 * t120 + t188 * t243 + t189 * t245 - t241 * t526;
t20 = -t115 * t723 + t320 * t117 + t321 * t119 + t188 * t244 + t189 * t246 - t242 * t526;
t38 = -t178 * t723 + t320 * t179 + t321 * t180 + t188 * t290 + t189 * t291 - t289 * t526;
t94 = -t241 * t723 + t243 * t320 + t245 * t321;
t95 = -t242 * t723 + t244 * t320 + t246 * t321;
t586 = t502 * t94 + t505 * t95;
t2 = -t38 * t504 - t95 * t640 + t586 * t677 + (qJD(2) * t135 + t19 * t502 + (qJD(1) * t94 + t20) * t505) * t501;
t109 = t242 * t504 - t244 * t397 + t246 * t398;
t103 = t109 * t640;
t108 = t241 * t504 - t243 * t397 + t245 * t398;
t28 = t116 * t504 - t118 * t397 + t120 * t398 - t241 * t679 + t243 * t281 + t245 * t282;
t29 = t115 * t504 - t117 * t397 + t119 * t398 - t242 * t679 + t244 * t281 + t246 * t282;
t515 = t504 * t178 - t397 * t179 + t398 * t180 + t281 * t290 + t282 * t291 - t289 * t679;
t49 = t515 * t504;
t568 = t108 * t502 + t109 * t505;
t151 = t289 * t504 - t290 * t397 + t291 * t398;
t670 = t151 * qJD(2);
t3 = -t103 - t49 + t568 * t677 + (t670 + t28 * t502 + (qJD(1) * t108 + t29) * t505) * t501;
t42 = -t135 * t504 + t501 * t586;
t43 = -t136 * t504 + t501 * t585;
t48 = -t151 * t504 + t501 * t568;
t796 = -(qJD(1) * (t42 * t505 - t43 * t502) + qJD(2) * t48 + t505 * t1 + t502 * t2) * t501 - (qJD(2) * (t42 * t502 + t43 * t505) - t3) * t504;
t750 = Icges(3,4) * t501;
t581 = Icges(3,1) * t504 - t750;
t416 = Icges(3,5) * t502 + t505 * t581;
t729 = t416 * t504;
t749 = Icges(3,4) * t504;
t577 = -Icges(3,2) * t501 + t749;
t412 = Icges(3,6) * t502 + t505 * t577;
t734 = t412 * t501;
t549 = -t729 + t734;
t795 = t502 * t549;
t415 = -Icges(3,5) * t505 + t502 * t581;
t731 = t415 * t504;
t411 = -Icges(3,6) * t505 + t502 * t577;
t736 = t411 * t501;
t550 = -t731 + t736;
t794 = t505 * t550;
t248 = t323 * rSges(7,1) + t322 * rSges(7,2) - rSges(7,3) * t721;
t483 = pkin(5) * t498 + pkin(4);
t499 = -pkin(9) - qJ(5);
t727 = t442 * t497;
t793 = pkin(5) * t727 + t443 * t483 + t499 * t721 + t248;
t792 = -rSges(3,2) * t721 + t502 * rSges(3,3);
t230 = t316 * t498 - t317 * t497;
t739 = t316 * t497;
t231 = t317 * t498 + t739;
t142 = Icges(6,5) * t231 + Icges(6,6) * t230 - Icges(6,3) * t526;
t144 = Icges(6,4) * t231 + Icges(6,2) * t230 - Icges(6,6) * t526;
t146 = Icges(6,1) * t231 + Icges(6,4) * t230 - Icges(6,5) * t526;
t364 = t440 * t498 - t441 * t497;
t728 = t440 * t497;
t365 = t441 * t498 + t728;
t258 = Icges(6,5) * t365 + Icges(6,6) * t364 - Icges(6,3) * t723;
t260 = Icges(6,4) * t365 + Icges(6,2) * t364 - Icges(6,6) * t723;
t262 = Icges(6,1) * t365 + Icges(6,4) * t364 - Icges(6,5) * t723;
t336 = qJD(3) * t427 - t546 * t677;
t337 = t545 * t677 + t546 * t674;
t39 = t142 * t504 - t144 * t426 + t146 * t427 - t258 * t679 + t260 * t336 + t262 * t337;
t209 = Icges(5,5) * t317 + Icges(5,6) * t526 + Icges(5,3) * t316;
t213 = Icges(5,4) * t317 + Icges(5,2) * t526 + Icges(5,6) * t316;
t217 = Icges(5,1) * t317 + Icges(5,4) * t526 + Icges(5,5) * t316;
t324 = Icges(5,5) * t441 + Icges(5,6) * t723 + Icges(5,3) * t440;
t328 = Icges(5,4) * t441 + Icges(5,2) * t723 + Icges(5,6) * t440;
t332 = Icges(5,1) * t441 + Icges(5,4) * t723 + Icges(5,5) * t440;
t558 = t324 * t500 + t332 * t503;
t71 = (qJD(2) * t558 - t213) * t504 + (qJD(2) * t328 + t209 * t500 + t217 * t503 + (t324 * t503 - t332 * t500) * qJD(3)) * t501;
t211 = Icges(4,5) * t317 - Icges(4,6) * t316 + Icges(4,3) * t526;
t215 = Icges(4,4) * t317 - Icges(4,2) * t316 + Icges(4,6) * t526;
t219 = Icges(4,1) * t317 - Icges(4,4) * t316 + Icges(4,5) * t526;
t326 = Icges(4,5) * t441 - Icges(4,6) * t440 + Icges(4,3) * t723;
t330 = Icges(4,4) * t441 - Icges(4,2) * t440 + Icges(4,6) * t723;
t334 = Icges(4,1) * t441 - Icges(4,4) * t440 + Icges(4,5) * t723;
t556 = -t330 * t500 + t334 * t503;
t73 = (qJD(2) * t556 - t211) * t504 + (qJD(2) * t326 - t215 * t500 + t219 * t503 + (-t330 * t503 - t334 * t500) * qJD(3)) * t501;
t791 = -t39 - t71 - t73;
t228 = -t314 * t498 - t315 * t497;
t740 = t314 * t497;
t229 = t315 * t498 - t740;
t141 = Icges(6,5) * t229 + Icges(6,6) * t228 + Icges(6,3) * t525;
t143 = Icges(6,4) * t229 + Icges(6,2) * t228 + Icges(6,6) * t525;
t145 = Icges(6,1) * t229 + Icges(6,4) * t228 + Icges(6,5) * t525;
t366 = t442 * t498 - t443 * t497;
t367 = t443 * t498 + t727;
t259 = Icges(6,5) * t367 + Icges(6,6) * t366 - Icges(6,3) * t721;
t261 = Icges(6,4) * t367 + Icges(6,2) * t366 - Icges(6,6) * t721;
t263 = Icges(6,1) * t367 + Icges(6,4) * t366 - Icges(6,5) * t721;
t40 = t141 * t504 - t143 * t426 + t145 * t427 - t259 * t679 + t261 * t336 + t263 * t337;
t208 = Icges(5,5) * t315 - Icges(5,6) * t525 - Icges(5,3) * t314;
t212 = Icges(5,4) * t315 - Icges(5,2) * t525 - Icges(5,6) * t314;
t216 = Icges(5,1) * t315 - Icges(5,4) * t525 - Icges(5,5) * t314;
t325 = Icges(5,5) * t443 + Icges(5,6) * t721 + Icges(5,3) * t442;
t329 = Icges(5,4) * t443 + Icges(5,2) * t721 + Icges(5,6) * t442;
t333 = Icges(5,1) * t443 + Icges(5,4) * t721 + Icges(5,5) * t442;
t557 = t325 * t500 + t333 * t503;
t72 = (qJD(2) * t557 - t212) * t504 + (qJD(2) * t329 + t208 * t500 + t216 * t503 + (t325 * t503 - t333 * t500) * qJD(3)) * t501;
t210 = Icges(4,5) * t315 + Icges(4,6) * t314 - Icges(4,3) * t525;
t214 = Icges(4,4) * t315 + Icges(4,2) * t314 - Icges(4,6) * t525;
t218 = Icges(4,1) * t315 + Icges(4,4) * t314 - Icges(4,5) * t525;
t327 = Icges(4,5) * t443 - Icges(4,6) * t442 + Icges(4,3) * t721;
t331 = Icges(4,4) * t443 - Icges(4,2) * t442 + Icges(4,6) * t721;
t335 = Icges(4,1) * t443 - Icges(4,4) * t442 + Icges(4,5) * t721;
t555 = -t331 * t500 + t335 * t503;
t74 = (qJD(2) * t555 - t210) * t504 + (qJD(2) * t327 - t214 * t500 + t218 * t503 + (-t331 * t503 - t335 * t500) * qJD(3)) * t501;
t790 = t40 + t72 + t74;
t139 = -t310 * t723 + t311 * t364 + t312 * t365;
t235 = t405 * t440 + t409 * t723 + t413 * t441;
t236 = t406 * t723 - t410 * t440 + t414 * t441;
t167 = t326 * t723 - t330 * t440 + t334 * t441;
t168 = t327 * t723 - t331 * t440 + t335 * t441;
t564 = t167 * t502 + t168 * t505;
t165 = t324 * t440 + t328 * t723 + t332 * t441;
t166 = t325 * t440 + t329 * t723 + t333 * t441;
t566 = t165 * t502 + t166 * t505;
t98 = -t258 * t723 + t260 * t364 + t262 * t365;
t99 = -t259 * t723 + t261 * t364 + t263 * t365;
t584 = t502 * t98 + t505 * t99;
t789 = (-t139 - t235 - t236) * t504 + (t564 + t566 + t584) * t501;
t140 = -t310 * t721 + t366 * t311 + t367 * t312;
t237 = t442 * t405 + t409 * t721 + t443 * t413;
t238 = t406 * t721 - t442 * t410 + t443 * t414;
t171 = t326 * t721 - t442 * t330 + t443 * t334;
t172 = t327 * t721 - t442 * t331 + t443 * t335;
t560 = t171 * t502 + t172 * t505;
t169 = t442 * t324 + t328 * t721 + t443 * t332;
t170 = t442 * t325 + t329 * t721 + t443 * t333;
t562 = t169 * t502 + t170 * t505;
t100 = -t258 * t721 + t366 * t260 + t367 * t262;
t101 = -t259 * t721 + t366 * t261 + t367 * t263;
t569 = t100 * t502 + t101 * t505;
t788 = (-t140 - t237 - t238) * t504 + (t560 + t562 + t569) * t501;
t252 = Icges(6,5) * t337 + Icges(6,6) * t336 - Icges(6,3) * t679;
t253 = Icges(6,4) * t337 + Icges(6,2) * t336 - Icges(6,6) * t679;
t254 = Icges(6,1) * t337 + Icges(6,4) * t336 - Icges(6,5) * t679;
t344 = (Icges(5,3) * t503 - t745) * t674 + (Icges(5,6) * t501 + t504 * t571) * qJD(2);
t345 = (-Icges(4,5) * t500 - Icges(4,6) * t503) * t674 + (Icges(4,3) * t501 + t504 * t572) * qJD(2);
t348 = (-Icges(5,4) * t500 + Icges(5,6) * t503) * t674 + (Icges(5,2) * t501 + t504 * t574) * qJD(2);
t352 = (-Icges(5,1) * t500 + t744) * t674 + (Icges(5,4) * t501 + t504 * t578) * qJD(2);
t353 = (-Icges(4,1) * t500 - t747) * t674 + (Icges(4,5) * t501 + t504 * t579) * qJD(2);
t722 = t501 * t503;
t724 = t500 * t501;
t787 = -t426 * t253 + t427 * t254 + t336 * t311 + t337 * t312 + t344 * t724 + t524 * t405 - t410 * t627 + t798 * t413 + t414 * t630 + (t352 + t353) * t722 + t800 * t679 + (t252 - t345 - t348) * t504;
t123 = t258 * t504 - t260 * t426 + t262 * t427;
t193 = -t328 * t504 + t501 * t558;
t195 = -t326 * t504 + t501 * t556;
t660 = t123 + t193 + t195;
t124 = t259 * t504 - t261 * t426 + t263 * t427;
t194 = -t329 * t504 + t501 * t557;
t196 = -t327 * t504 + t501 * t555;
t659 = t124 + t194 + t196;
t573 = Icges(3,5) * t504 - Icges(3,6) * t501;
t407 = -Icges(3,3) * t505 + t502 * t573;
t785 = t504 * t659 - t43 - t788;
t782 = 2 * m(3);
t781 = 2 * m(4);
t780 = 2 * m(5);
t779 = 2 * m(6);
t778 = 2 * m(7);
t777 = m(6) / 0.2e1;
t776 = m(7) / 0.2e1;
t775 = -pkin(3) - pkin(4);
t774 = -t136 / 0.2e1;
t773 = t502 / 0.2e1;
t772 = -t504 / 0.2e1;
t771 = -t505 / 0.2e1;
t770 = t505 / 0.2e1;
t769 = -rSges(5,1) - pkin(3);
t768 = -rSges(5,2) - pkin(8);
t767 = -rSges(4,3) - pkin(8);
t766 = rSges(6,3) - pkin(8);
t460 = rSges(3,1) * t501 + rSges(3,2) * t504;
t765 = m(3) * t460;
t764 = pkin(2) * t501;
t763 = pkin(2) * t504;
t490 = t502 * pkin(7);
t762 = t504 * (qJD(1) * t568 - t28 * t505 + t29 * t502);
t760 = qJD(1) / 0.2e1;
t758 = -pkin(3) - t483;
t757 = pkin(4) - t483;
t756 = -pkin(8) - t499;
t755 = rSges(3,3) * t505;
t754 = rSges(5,3) * t440;
t753 = rSges(7,3) * t504;
t752 = t316 * rSges(5,3);
t751 = -rSges(6,3) - qJ(5);
t591 = -rSges(6,1) * t365 - rSges(6,2) * t364;
t264 = -rSges(6,3) * t723 - t591;
t741 = t264 * t505;
t340 = rSges(5,1) * t441 + rSges(5,2) * t723 + t754;
t738 = t340 * t505;
t595 = -rSges(4,1) * t441 + rSges(4,2) * t440;
t341 = rSges(4,3) * t723 - t595;
t737 = t341 * t505;
t735 = t411 * t504;
t733 = t412 * t504;
t732 = t415 * t501;
t730 = t416 * t501;
t717 = qJ(5) + t499;
t656 = t187 * rSges(7,1) + t186 * rSges(7,2) + rSges(7,3) * t640;
t121 = -rSges(7,3) * t629 + t656;
t651 = -pkin(5) * t740 + t315 * t483 + t499 * t629;
t706 = t315 * pkin(4) + qJ(5) * t640;
t716 = qJ(5) * t629 - t499 * t640 + t121 + t651 - t706;
t590 = t189 * rSges(7,1) + t188 * rSges(7,2);
t122 = -rSges(7,3) * t526 + t590;
t665 = pkin(5) * t739;
t689 = t526 * qJ(5);
t715 = -t317 * t757 + t499 * t526 + t122 + t665 + t689;
t181 = rSges(7,1) * t282 + rSges(7,2) * t281 - rSges(7,3) * t679;
t529 = pkin(5) * t726 - t503 * t757;
t714 = t181 + (pkin(5) * t725 + t500 * t757) * t674 + (t501 * t717 + t504 * t529) * qJD(2);
t707 = -t316 * qJ(4) - t440 * qJD(4);
t221 = t317 * pkin(3) - t707;
t432 = t440 * qJ(4);
t375 = pkin(3) * t441 + t432;
t713 = t221 * t721 + t375 * t629;
t220 = t315 * pkin(3) - t314 * qJ(4) + t442 * qJD(4);
t649 = t315 * rSges(5,1) + rSges(5,2) * t629 - t314 * rSges(5,3);
t222 = -rSges(5,2) * t640 + t649;
t712 = -t220 - t222;
t671 = qJD(5) * t501;
t283 = (-qJ(5) * t677 - t671) * t505 + t706;
t711 = -t220 - t283;
t589 = -rSges(7,1) * t321 - rSges(7,2) * t320;
t247 = -rSges(7,3) * t723 - t589;
t477 = qJ(5) * t723;
t666 = pkin(5) * t728;
t710 = -t441 * t757 + t499 * t723 + t247 + t477 + t666;
t438 = t443 * pkin(4);
t394 = -qJ(5) * t721 + t438;
t709 = -t394 + t793;
t294 = rSges(7,1) * t398 - rSges(7,2) * t397 + t753;
t708 = t501 * t529 - t504 * t717 + t294;
t339 = (pkin(3) * t677 + qJ(4) * t674) * t503 + (qJ(4) * t677 + (-pkin(3) * qJD(3) + qJD(4)) * t501) * t500;
t593 = rSges(5,1) * t503 + rSges(5,3) * t500;
t359 = (-rSges(5,1) * t500 + rSges(5,3) * t503) * t674 + (rSges(5,2) * t501 + t504 * t593) * qJD(2);
t705 = -t339 - t359;
t386 = t798 * pkin(4) - qJ(5) * t679 + qJD(5) * t504;
t704 = -t339 - t386;
t703 = -t340 - t375;
t342 = t443 * rSges(5,1) + rSges(5,2) * t721 + t442 * rSges(5,3);
t376 = t443 * pkin(3) + t442 * qJ(4);
t702 = -t342 - t376;
t356 = t375 * t721;
t393 = pkin(4) * t441 - t477;
t701 = t393 * t721 + t356;
t700 = t367 * rSges(6,1) + t366 * rSges(6,2);
t594 = rSges(4,1) * t503 - rSges(4,2) * t500;
t360 = (-rSges(4,1) * t500 - rSges(4,2) * t503) * t674 + (rSges(4,3) * t501 + t504 * t594) * qJD(2);
t599 = pkin(8) * t501 + t763;
t453 = t599 * qJD(2);
t699 = -t360 - t453;
t444 = (pkin(3) * t503 + qJ(4) * t500) * t501;
t420 = t444 * t682;
t698 = t376 * t679 + t501 * t420;
t697 = t504 * t375 + t444 * t723;
t696 = -t375 - t393;
t695 = -t376 - t394;
t461 = -pkin(8) * t504 + t764;
t447 = t461 * t682;
t694 = t420 + t447;
t421 = -rSges(5,2) * t504 + t501 * t593;
t693 = -t421 - t444;
t422 = -rSges(4,3) * t504 + t501 * t594;
t692 = -t422 - t461;
t445 = t599 * t502;
t482 = pkin(2) * t718;
t446 = pkin(8) * t721 + t482;
t691 = t502 * t445 + t505 * t446;
t452 = pkin(4) * t722 + qJ(5) * t504;
t690 = -t444 - t452;
t688 = rSges(3,2) * t640 + rSges(3,3) * t680;
t472 = pkin(8) * t629;
t486 = pkin(7) * t680;
t687 = t472 + t486;
t686 = t505 * pkin(1) + t490;
t491 = t505 * pkin(7);
t685 = t491 - t432;
t684 = t502 ^ 2 + t505 ^ 2;
t408 = Icges(3,3) * t502 + t505 * t573;
t683 = qJD(1) * t408;
t669 = t49 + t103 / 0.2e1;
t668 = t777 + t776;
t667 = rSges(7,3) + t756;
t663 = -t28 / 0.2e1 - t38 / 0.2e1;
t662 = -t29 / 0.2e1 - t37 / 0.2e1;
t349 = (-Icges(4,2) * t503 - t748) * t674 + (Icges(4,6) * t501 + t504 * t575) * qJD(2);
t661 = (-t410 * t677 + (-qJD(3) * t414 - t349) * t501) * t500 + t787;
t655 = t229 * rSges(6,1) + t228 * rSges(6,2) + rSges(6,3) * t640;
t147 = -rSges(6,3) * t629 + t655;
t658 = -t147 + t711;
t657 = t797 * t679;
t255 = rSges(6,1) * t337 + rSges(6,2) * t336 - rSges(6,3) * t679;
t654 = -t255 + t704;
t653 = -t264 + t696;
t265 = -rSges(6,3) * t721 + t700;
t652 = -t265 + t695;
t650 = t315 * rSges(4,1) + t314 * rSges(4,2) + rSges(4,3) * t629;
t313 = rSges(6,1) * t427 - rSges(6,2) * t426 + rSges(6,3) * t504;
t648 = -t313 + t690;
t647 = -t453 + t705;
t471 = pkin(2) * t634;
t527 = -t502 * t681 - t632;
t646 = t502 * (pkin(8) * t526 + qJD(1) * t482 - t471) + t505 * (pkin(2) * t527 - pkin(8) * t640 + t472) + t445 * t680;
t437 = t452 * t682;
t644 = t437 + t694;
t643 = -t461 + t693;
t343 = t443 * rSges(4,1) - t442 * rSges(4,2) + rSges(4,3) * t721;
t642 = t471 + t707;
t641 = -pkin(1) - t763;
t638 = t313 * t682;
t637 = t421 * t682;
t636 = t422 * t682;
t633 = t501 * t677;
t625 = -t108 / 0.2e1 - t135 / 0.2e1;
t624 = -t677 / 0.2e1;
t622 = t502 * t708;
t621 = t710 * t505;
t620 = t693 * t505;
t374 = t692 * t505;
t619 = qJD(1) * t708;
t618 = pkin(2) * t632;
t616 = t711 - t716;
t615 = t704 - t714;
t475 = t502 * t671;
t284 = t317 * pkin(4) - t475 - t689;
t614 = t284 * t721 + t393 * t629 + t713;
t613 = t504 * t221 + t339 * t723 + t526 * t444;
t612 = t696 - t710;
t611 = t695 - t709;
t610 = -t453 + t654;
t609 = t690 - t708;
t608 = -t461 + t648;
t607 = t394 * t679 + t501 * t437 + t698;
t606 = t502 * t375 + t505 * t376 + t691;
t605 = t504 * t393 + t452 * t723 + t697;
t604 = t475 + t642;
t603 = t686 + t446;
t602 = -t1 / 0.2e1 + t42 * t761;
t601 = t2 / 0.2e1 + t43 * t761;
t600 = t648 * t505;
t293 = t643 * t505;
t598 = t502 * t619;
t597 = rSges(3,1) * t504 - rSges(3,2) * t501;
t596 = t317 * rSges(4,1) - t316 * rSges(4,2);
t592 = t231 * rSges(6,1) + t230 * rSges(6,2);
t56 = t99 * t502 - t505 * t98;
t583 = -t453 + t615;
t582 = -t461 + t609;
t576 = Icges(3,2) * t504 + t750;
t570 = t100 * t505 - t101 * t502;
t567 = t165 * t505 - t166 * t502;
t565 = t167 * t505 - t168 * t502;
t563 = t169 * t505 - t170 * t502;
t561 = t171 * t505 - t172 * t502;
t559 = -t247 * t505 + t248 * t502;
t554 = -t343 * t502 + t737;
t553 = -t502 * t341 - t343 * t505;
t424 = rSges(3,1) * t718 + t792;
t543 = t609 * t505;
t251 = t608 * t505;
t542 = -pkin(1) - t597;
t51 = t95 * t502 - t505 * t94;
t541 = -t567 / 0.2e1 - t565 / 0.2e1 + t56 / 0.2e1 + t51 / 0.2e1;
t52 = t97 * t502 - t505 * t96;
t540 = t52 / 0.2e1 - t563 / 0.2e1 - t561 / 0.2e1 - t570 / 0.2e1;
t539 = t505 * t220 + t502 * t221 + t375 * t680 + t646;
t538 = t502 * t393 + t505 * t394 + t606;
t537 = qJD(2) * t460;
t191 = t582 * t505;
t536 = t502 * t702 + t738;
t535 = t501 * t768 + t641;
t534 = t501 * t767 + t641;
t533 = t501 * t766 + t641;
t531 = qJD(2) * t576;
t530 = qJD(2) * (-Icges(3,5) * t501 - Icges(3,6) * t504);
t523 = t376 + t603;
t522 = t502 * t652 + t741;
t521 = t501 * t667 + t641;
t520 = t504 * t284 + t386 * t723 + t526 * t452 + t613;
t519 = t535 * t502;
t518 = t534 * t502;
t517 = t220 + t687;
t516 = -t504 * t660 + t42 + t789;
t6 = qJD(1) * t585 - t17 * t505 + t18 * t502;
t7 = qJD(1) * t586 - t19 * t505 + t20 * t502;
t514 = -qJD(2) * (-t108 * t505 + t109 * t502) / 0.2e1 - t502 * t7 / 0.2e1 + t6 * t771;
t513 = t505 * t283 + t502 * t284 + t393 * t680 + t539;
t512 = t502 * t611 + t621;
t510 = t774 - t194 / 0.2e1 - t196 / 0.2e1 - t124 / 0.2e1 - t237 / 0.2e1 - t238 / 0.2e1 - t140 / 0.2e1;
t104 = -t314 * t405 + t315 * t413 + t442 * t344 + t348 * t721 + t443 * t352 - t409 * t525;
t105 = t314 * t410 + t315 * t414 + t345 * t721 - t442 * t349 + t443 * t353 - t406 * t525;
t53 = t228 * t311 + t229 * t312 - t252 * t721 + t366 * t253 + t367 * t254 + t310 * t525;
t509 = t104 / 0.2e1 + t105 / 0.2e1 + t40 / 0.2e1 + t53 / 0.2e1 + t72 / 0.2e1 + t74 / 0.2e1 - t662;
t106 = t316 * t405 + t317 * t413 + t440 * t344 + t348 * t723 + t441 * t352 + t409 * t526;
t107 = -t316 * t410 + t317 * t414 + t345 * t723 - t440 * t349 + t441 * t353 + t406 * t526;
t54 = t230 * t311 + t231 * t312 - t252 * t723 + t364 * t253 + t365 * t254 - t310 * t526;
t508 = t106 / 0.2e1 + t107 / 0.2e1 + t39 / 0.2e1 + t54 / 0.2e1 + t71 / 0.2e1 + t73 / 0.2e1 - t663;
t507 = t193 / 0.2e1 + t195 / 0.2e1 + t123 / 0.2e1 + t236 / 0.2e1 + t235 / 0.2e1 + t139 / 0.2e1 - t625;
t506 = t109 / 0.2e1 - t510;
t494 = t501 ^ 2;
t451 = t597 * qJD(2);
t448 = t573 * qJD(2);
t423 = t502 * t597 - t755;
t388 = t424 + t686;
t387 = t502 * t542 + t491 + t755;
t373 = t692 * t502;
t347 = t502 * t530 + t683;
t346 = -qJD(1) * t407 + t505 * t530;
t297 = t460 * t678 + ((-rSges(3,3) - pkin(7)) * t502 + t542 * t505) * qJD(1);
t296 = rSges(3,1) * t527 - rSges(3,2) * t629 - pkin(1) * t682 + t486 + t688;
t292 = t643 * t502;
t288 = t603 + t343;
t287 = t491 + t518 + t595;
t286 = -t504 * t343 - t422 * t721;
t285 = t341 * t504 + t422 * t723;
t280 = t502 * t408 - t505 * t549;
t279 = t502 * t407 - t794;
t278 = -t408 * t505 - t795;
t277 = -t407 * t505 - t502 * t550;
t250 = t608 * t502;
t249 = t554 * t501;
t240 = qJD(1) * t374 + t502 * t699;
t239 = t505 * t699 + t447 + t636;
t234 = t523 + t342;
t233 = t441 * t769 + t519 + t685 - t754;
t225 = rSges(4,3) * t526 + t596;
t224 = t317 * rSges(5,1) + rSges(5,2) * t526 + t752;
t223 = -rSges(4,3) * t640 + t650;
t201 = -t553 + t691;
t198 = t501 * t620 + t504 * t702;
t197 = t340 * t504 + t421 * t723 + t697;
t190 = t582 * t502;
t183 = t471 + t767 * t631 + (t505 * t534 - t490) * qJD(1) - t596;
t182 = qJD(1) * t518 - t618 + t650 + t687;
t177 = t504 * t248 + t294 * t721;
t176 = -t247 * t504 - t294 * t723;
t175 = qJD(1) * t293 + t502 * t647;
t174 = t505 * t647 + t637 + t694;
t164 = t721 * t751 + t438 + t523 + t700;
t163 = t441 * t775 + t502 * t533 + t477 + t591 + t685;
t162 = t501 * t536 + t356;
t155 = t502 * t340 + t342 * t505 + t606;
t154 = t523 + t793;
t153 = t441 * t758 + t502 * t521 + t589 - t666 + t685;
t152 = t559 * t501;
t148 = -rSges(6,3) * t526 + t592;
t138 = t501 * t600 + t504 * t652;
t137 = t264 * t504 + t313 * t723 + t605;
t132 = (t422 * t678 + t225) * t504 + (-qJD(2) * t341 + t502 * t360 + t422 * t680) * t501;
t131 = (-t422 * t676 - t223) * t504 + (qJD(2) * t343 - t360 * t505 + t636) * t501;
t130 = -t752 + t769 * t317 + t768 * t631 + (t505 * t535 - t490) * qJD(1) + t642;
t129 = qJD(1) * t519 + t517 - t618 + t649;
t126 = qJD(1) * t251 + t502 * t610;
t125 = t505 * t610 + t638 + t644;
t114 = t501 * t522 + t701;
t102 = t502 * t264 + t265 * t505 + t538;
t93 = t501 * t543 + t504 * t611;
t92 = t501 * t622 + t504 * t710 + t605;
t89 = t554 * t677 + (qJD(1) * t553 - t223 * t502 + t225 * t505) * t501;
t88 = t223 * t505 + t502 * t225 + (t737 + (-t343 - t446) * t502) * qJD(1) + t646;
t83 = t775 * t317 + t766 * t631 + (t505 * t533 - t490) * qJD(1) - t592 + t604 + t689;
t82 = (-pkin(1) - t599) * t682 + (-t671 + (t504 * t751 - t764) * qJD(2)) * t505 + t517 + t655 + t706;
t81 = qJD(1) * t191 + t502 * t583;
t80 = t505 * t583 + t598 + t644;
t79 = (t421 * t678 + t224) * t504 + (qJD(2) * t703 + t502 * t359 + t421 * t680) * t501 + t613;
t78 = (qJD(2) * t620 + t712) * t504 + (qJD(2) * t342 + t505 * t705 + t637) * t501 + t698;
t77 = t501 * t512 + t701;
t76 = t502 * t710 + t505 * t709 + t538;
t70 = -t665 + t758 * t317 + t667 * t631 + (t505 * t521 - t490) * qJD(1) - t590 + t604;
t69 = (-t671 + (-t753 - t764) * qJD(2)) * t505 + (t501 * t756 + t641) * t682 + t517 + t651 + t656;
t68 = t210 * t723 - t440 * t214 + t441 * t218 - t316 * t331 + t317 * t335 + t327 * t526;
t67 = t211 * t723 - t440 * t215 + t441 * t219 - t316 * t330 + t317 * t334 + t326 * t526;
t66 = t440 * t208 + t212 * t723 + t441 * t216 + t316 * t325 + t317 * t333 + t329 * t526;
t65 = t440 * t209 + t213 * t723 + t441 * t217 + t316 * t324 + t317 * t332 + t328 * t526;
t64 = t210 * t721 - t442 * t214 + t443 * t218 + t314 * t331 + t315 * t335 - t327 * t525;
t63 = t211 * t721 - t442 * t215 + t443 * t219 + t314 * t330 + t315 * t334 - t326 * t525;
t62 = t442 * t208 + t212 * t721 + t443 * t216 - t314 * t325 + t315 * t333 - t329 * t525;
t61 = t442 * t209 + t213 * t721 + t443 * t217 - t314 * t324 + t315 * t332 - t328 * t525;
t59 = (-t294 * t678 - t122) * t504 + (qJD(2) * t247 - t502 * t181 - t294 * t680) * t501;
t58 = (t294 * t676 + t121) * t504 + (-qJD(2) * t248 + t181 * t505 - t294 * t682) * t501;
t55 = t222 * t505 + t502 * t224 + (t738 + (-t446 + t702) * t502) * qJD(1) + t539;
t50 = t536 * t677 + (t224 * t505 + t712 * t502 + (t502 * t703 + t505 * t702) * qJD(1)) * t501 + t713;
t45 = (t313 * t678 + t148) * t504 + (qJD(2) * t653 + t502 * t255 + t313 * t680) * t501 + t520;
t44 = (qJD(2) * t265 + t505 * t654 + t638) * t501 + (qJD(2) * t600 + t658) * t504 + t607;
t41 = t559 * t677 + (t121 * t502 - t122 * t505 + (t502 * t247 + t248 * t505) * qJD(1)) * t501;
t34 = t147 * t505 + t502 * t148 + (t741 + (-t446 + t652) * t502) * qJD(1) + t513;
t33 = -t141 * t723 + t364 * t143 + t365 * t145 + t230 * t261 + t231 * t263 - t259 * t526;
t32 = -t142 * t723 + t364 * t144 + t365 * t146 + t230 * t260 + t231 * t262 - t258 * t526;
t31 = -t141 * t721 + t366 * t143 + t367 * t145 + t228 * t261 + t229 * t263 + t259 * t525;
t30 = -t142 * t721 + t366 * t144 + t367 * t146 + t228 * t260 + t229 * t262 + t258 * t525;
t27 = t522 * t677 + (t148 * t505 + t658 * t502 + (t502 * t653 + t505 * t652) * qJD(1)) * t501 + t614;
t26 = (qJD(2) * t622 + t715) * t504 + (qJD(2) * t612 + t502 * t714 + t505 * t619) * t501 + t520;
t25 = (qJD(2) * t543 + t616) * t504 + (qJD(2) * t709 + t505 * t615 + t598) * t501 + t607;
t24 = qJD(1) * t564 + t68 * t502 - t505 * t67;
t23 = qJD(1) * t566 + t66 * t502 - t505 * t65;
t22 = qJD(1) * t560 + t64 * t502 - t505 * t63;
t21 = qJD(1) * t562 + t62 * t502 - t505 * t61;
t16 = t716 * t505 + t715 * t502 + (t621 + (-t446 + t611) * t502) * qJD(1) + t513;
t15 = (qJD(2) * t564 - t107) * t504 + (qJD(1) * t565 + qJD(2) * t236 + t502 * t67 + t505 * t68) * t501;
t14 = (qJD(2) * t566 - t106) * t504 + (qJD(1) * t567 + qJD(2) * t235 + t502 * t65 + t505 * t66) * t501;
t13 = (qJD(2) * t560 - t105) * t504 + (qJD(1) * t561 + qJD(2) * t238 + t502 * t63 + t505 * t64) * t501;
t12 = (qJD(2) * t562 - t104) * t504 + (qJD(1) * t563 + qJD(2) * t237 + t502 * t61 + t505 * t62) * t501;
t11 = t512 * t677 + (t715 * t505 + t616 * t502 + (t502 * t612 + t505 * t611) * qJD(1)) * t501 + t614;
t10 = qJD(1) * t584 - t32 * t505 + t33 * t502;
t9 = qJD(1) * t569 - t30 * t505 + t31 * t502;
t5 = (qJD(2) * t584 - t54) * t504 + (-qJD(1) * t56 + qJD(2) * t139 + t32 * t502 + t33 * t505) * t501;
t4 = (qJD(2) * t569 - t53) * t504 + (qJD(1) * t570 + qJD(2) * t140 + t30 * t502 + t31 * t505) * t501;
t8 = [t515 + (t129 * t234 + t130 * t233) * t780 + (t153 * t70 + t154 * t69) * t778 + (t163 * t83 + t164 * t82) * t779 + (t182 * t288 + t183 * t287) * t781 + (t296 * t388 + t297 * t387) * t782 - t410 * t635 - t414 * t628 - t349 * t724 + (-t576 + t581) * t679 + (Icges(3,1) * t501 + t577 + t749) * t677 + t787; m(6) * (t125 * t163 + t126 * t164 + t250 * t82 + t251 * t83) + m(7) * (t153 * t80 + t154 * t81 + t190 * t69 + t191 * t70) + m(5) * (t129 * t292 + t130 * t293 + t174 * t233 + t175 * t234) + m(4) * (t182 * t373 + t183 * t374 + t239 * t287 + t240 * t288) + ((qJD(1) * t412 - t502 * t531) * t772 + t416 * t799 + t448 * t770 + (t736 / 0.2e1 - t731 / 0.2e1) * qJD(2) - t508 + m(3) * (-t297 * t460 - t387 * t451) + (t733 / 0.2e1 + t730 / 0.2e1 - t388 * t765 + t506) * qJD(1)) * t505 + ((-qJD(1) * t411 - t505 * t531) * t504 / 0.2e1 + t415 * t799 + t448 * t773 + (-t734 / 0.2e1 + t729 / 0.2e1) * qJD(2) + t509 + m(3) * (-t296 * t460 - t388 * t451) + (t387 * t765 + t735 / 0.2e1 + t732 / 0.2e1 + t507) * qJD(1)) * t502; (t155 * t55 + t174 * t293 + t175 * t292) * t780 + (t201 * t88 + t239 * t374 + t240 * t373) * t781 + ((t502 * t423 + t424 * t505) * ((qJD(1) * t423 - t505 * t537 + t688) * t505 + (-t502 * t537 + (-t424 + t792) * qJD(1)) * t502) + t684 * t460 * t451) * t782 + (t16 * t76 + t190 * t81 + t191 * t80) * t778 + (t102 * t34 + t125 * t251 + t126 * t250) * t779 - t505 * t24 - t505 * t23 - t505 * t10 - t505 * t7 + t502 * t21 + t502 * t22 + t502 * t9 + t502 * t6 + t502 * ((t502 * t346 + (t279 + t795) * qJD(1)) * t502 + (t280 * qJD(1) + (t411 * t677 + t415 * t679) * t505 + (-t347 + (-t730 - t733) * qJD(2) + (t408 - t550) * qJD(1)) * t502) * t505) - t505 * ((t505 * t347 + (t278 + t794) * qJD(1)) * t505 + (t277 * qJD(1) + (-t412 * t677 - t416 * t679 + t683) * t502 + (-t346 + (t732 + t735) * qJD(2) - t549 * qJD(1)) * t505) * t502) + (-t277 * t505 + t278 * t502 + t51 + t56 - t565 - t567) * t682 + (-t279 * t505 + t280 * t502 + t52 - t561 - t563 - t570) * t680; m(5) * (t129 * t198 + t130 * t197 + t233 * t79 + t234 * t78) + m(7) * (t153 * t26 + t154 * t25 + t69 * t93 + t70 * t92) + m(6) * (t137 * t83 + t138 * t82 + t163 * t45 + t164 * t44) + m(4) * (t131 * t288 + t132 * t287 + t182 * t286 + t183 * t285) + ((t502 * t507 + t505 * t506) * qJD(2) - t661) * t504 + (t670 + (qJD(1) * t510 + t508) * t502 + (qJD(1) * t507 + t509) * t505) * t501 - t657 - t669; -t762 / 0.2e1 + m(6) * (t102 * t27 + t114 * t34 + t125 * t137 + t126 * t138 + t250 * t44 + t251 * t45) + m(7) * (t11 * t76 + t16 * t77 + t190 * t25 + t191 * t26 + t80 * t92 + t81 * t93) + m(5) * (t155 * t50 + t162 * t55 + t174 * t197 + t175 * t198 + t292 * t78 + t293 * t79) + m(4) * (t131 * t373 + t132 * t374 + t201 * t89 + t239 * t285 + t240 * t286 + t249 * t88) + (-t15 / 0.2e1 - t14 / 0.2e1 - t5 / 0.2e1 + t540 * t677 - t601) * t505 + (t12 / 0.2e1 + t13 / 0.2e1 + t4 / 0.2e1 + t541 * t677 - t602) * t502 + (t502 * t659 - t505 * t660) * t679 / 0.2e1 + (t791 * t505 + t790 * t502 + (t502 * t660 + t505 * t659) * qJD(1)) * t772 + (t502 * t789 + t505 * t788) * t760 + ((-t502 * t540 + t505 * t541) * qJD(1) - t514 + (t24 + t10 + t23) * t773 + (t9 + t21 + t22) * t770) * t501; (t131 * t286 + t132 * t285 + t249 * t89) * t781 + (t162 * t50 + t197 * t79 + t198 * t78) * t780 + (t114 * t27 + t137 * t45 + t138 * t44) * t779 + (t11 * t77 + t25 * t93 + t26 * t92) * t778 + (-t3 + t661 * t504 + (t516 * t502 - t505 * t785) * qJD(2) + t657) * t504 + ((-t504 * t790 + t1 + t12 + t13 + t4) * t505 + (t504 * t791 + t14 + t15 + t2 + t5) * t502 + (t797 * t504 + t659 * t721 + t660 * t723 + t48) * qJD(2) + (t502 * t785 + t516 * t505) * qJD(1)) * t501; m(7) * (-t153 * t314 + t154 * t316 + t440 * t69 + t442 * t70) + m(6) * (-t163 * t314 + t164 * t316 + t440 * t82 + t442 * t83) + m(5) * (t129 * t440 + t130 * t442 - t233 * t314 + t234 * t316); m(7) * (t76 * t627 + t190 * t316 - t191 * t314 + t440 * t81 + t442 * t80 + (t16 * t501 + t677 * t76) * t500) + m(6) * (t102 * t524 + t125 * t442 + t126 * t440 + t250 * t316 - t251 * t314 + t34 * t724) + m(5) * (t155 * t524 + t174 * t442 + t175 * t440 + t292 * t316 - t293 * t314 + t55 * t724); m(7) * (t77 * t627 + t25 * t440 + t26 * t442 - t314 * t92 + t316 * t93 + (t11 * t501 + t677 * t77) * t500) + m(6) * (t114 * t524 - t137 * t314 + t138 * t316 + t27 * t724 + t44 * t440 + t442 * t45) + m(5) * (t162 * t524 - t197 * t314 + t198 * t316 + t440 * t78 + t442 * t79 + t50 * t724); 0.4e1 * (m(5) / 0.2e1 + t668) * (-t314 * t442 + t316 * t440 + (t494 * t673 + t500 * t633) * t500); 0.2e1 * ((-t153 * t505 - t154 * t502) * t776 + (-t163 * t505 - t164 * t502) * t777) * t677 + 0.2e1 * ((t153 * t682 - t154 * t680 - t502 * t69 - t505 * t70) * t776 + (t163 * t682 - t164 * t680 - t502 * t82 - t505 * t83) * t777) * t501; 0.2e1 * ((-t190 * t678 - t191 * t676 + t16) * t776 + (-t250 * t678 - t251 * t676 + t34) * t777) * t504 + 0.2e1 * ((-qJD(2) * t76 - t190 * t680 + t191 * t682 - t502 * t81 - t505 * t80) * t776 + (-qJD(2) * t102 - t125 * t505 - t126 * t502 - t250 * t680 + t251 * t682) * t777) * t501; 0.2e1 * ((-t676 * t92 - t678 * t93 + t11) * t776 + (-t137 * t676 - t138 * t678 + t27) * t777) * t504 + 0.2e1 * ((-qJD(2) * t77 - t25 * t502 - t26 * t505 - t680 * t93 + t682 * t92) * t776 + (-qJD(2) * t114 + t137 * t682 - t138 * t680 - t44 * t502 - t45 * t505) * t777) * t501; 0.2e1 * t668 * ((-t494 * t500 + (-t440 * t502 - t442 * t505 + t500 * t504) * t504) * qJD(2) + (t626 + t314 * t505 - t316 * t502 + (-t440 * t505 + t442 * t502) * qJD(1)) * t501); 0.4e1 * t668 * (-0.1e1 + t684) * t633; m(7) * (t153 * t59 + t154 * t58 + t176 * t70 + t177 * t69) + ((-t109 / 0.2e1 + t774) * t505 + t625 * t502) * t677 + (-t670 + (t136 * t760 + t663) * t502 + (qJD(1) * t625 + t662) * t505) * t501 + t669; m(7) * (t152 * t16 + t176 * t80 + t177 * t81 + t190 * t58 + t191 * t59 + t41 * t76) + t762 / 0.2e1 + (t52 * t624 + t601) * t505 + (t51 * t624 + t602) * t502 + ((t51 * t771 + t52 * t773) * qJD(1) + t514) * t501; m(7) * (t11 * t152 + t176 * t26 + t177 * t25 + t41 * t77 + t58 * t93 + t59 * t92) + t796; m(7) * (t152 * t524 - t176 * t314 + t177 * t316 + t41 * t724 + t440 * t58 + t442 * t59); m(7) * ((t41 + (-t176 * t505 - t177 * t502) * qJD(2)) * t504 + (-qJD(2) * t152 - t502 * t58 - t505 * t59 + (t176 * t502 - t177 * t505) * qJD(1)) * t501); (t152 * t41 + t176 * t59 + t177 * t58) * t778 - t796;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t8(1) t8(2) t8(4) t8(7) t8(11) t8(16); t8(2) t8(3) t8(5) t8(8) t8(12) t8(17); t8(4) t8(5) t8(6) t8(9) t8(13) t8(18); t8(7) t8(8) t8(9) t8(10) t8(14) t8(19); t8(11) t8(12) t8(13) t8(14) t8(15) t8(20); t8(16) t8(17) t8(18) t8(19) t8(20) t8(21);];
Mq  = res;
