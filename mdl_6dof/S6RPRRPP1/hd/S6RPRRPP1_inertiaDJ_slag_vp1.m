% Calculate time derivative of joint inertia matrix for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:55
% EndTime: 2019-03-09 04:28:25
% DurationCPUTime: 17.88s
% Computational Cost: add. (41626->994), mult. (44757->1362), div. (0->0), fcn. (42552->10), ass. (0->477)
t365 = sin(qJ(3));
t368 = cos(qJ(3));
t360 = qJ(4) + pkin(10);
t355 = sin(t360);
t357 = cos(t360);
t418 = Icges(6,5) * t357 - Icges(6,6) * t355;
t284 = -Icges(6,3) * t368 + t365 * t418;
t421 = Icges(7,4) * t357 + Icges(7,6) * t355;
t285 = -Icges(7,2) * t368 + t365 * t421;
t364 = sin(qJ(4));
t367 = cos(qJ(4));
t419 = Icges(5,5) * t367 - Icges(5,6) * t364;
t297 = -Icges(5,3) * t368 + t365 * t419;
t626 = t297 + t284 + t285;
t361 = qJ(1) + pkin(9);
t358 = cos(t361);
t356 = sin(t361);
t546 = t356 * t368;
t278 = t355 * t546 + t357 * t358;
t279 = -t355 * t358 + t357 * t546;
t619 = rSges(7,3) + qJ(6);
t620 = rSges(7,1) + pkin(5);
t625 = t619 * t278 + t279 * t620;
t509 = qJD(3) * t365;
t474 = t356 * t509;
t506 = qJD(4) * t357;
t507 = qJD(4) * t355;
t514 = qJD(1) * t358;
t515 = qJD(1) * t356;
t174 = -t355 * t474 - t358 * t507 - t357 * t515 + (t355 * t514 + t356 * t506) * t368;
t503 = qJD(4) * t368;
t458 = qJD(1) - t503;
t403 = t458 * t355;
t512 = qJD(1) * t368;
t457 = -qJD(4) + t512;
t404 = t358 * t457;
t175 = t357 * t404 + (-t357 * t509 + t403) * t356;
t624 = -t278 * qJD(6) - t619 * t174 - t175 * t620;
t513 = qJD(1) * t365;
t539 = t364 * t368;
t293 = -t356 * t539 - t358 * t367;
t538 = t367 * t368;
t544 = t358 * t364;
t294 = t356 * t538 - t544;
t547 = t356 * t365;
t208 = Icges(5,5) * t294 + Icges(5,6) * t293 + Icges(5,3) * t547;
t210 = Icges(5,4) * t294 + Icges(5,2) * t293 + Icges(5,6) * t547;
t212 = Icges(5,1) * t294 + Icges(5,4) * t293 + Icges(5,5) * t547;
t295 = t356 * t367 - t358 * t539;
t548 = t356 * t364;
t296 = t358 * t538 + t548;
t543 = t358 * t365;
t86 = t208 * t543 + t210 * t295 + t212 * t296;
t209 = Icges(5,5) * t296 + Icges(5,6) * t295 + Icges(5,3) * t543;
t211 = Icges(5,4) * t296 + Icges(5,2) * t295 + Icges(5,6) * t543;
t213 = Icges(5,1) * t296 + Icges(5,4) * t295 + Icges(5,5) * t543;
t87 = t209 * t543 + t211 * t295 + t213 * t296;
t432 = t356 * t86 + t358 * t87;
t178 = Icges(6,5) * t279 - Icges(6,6) * t278 + Icges(6,3) * t547;
t182 = Icges(6,4) * t279 - Icges(6,2) * t278 + Icges(6,6) * t547;
t186 = Icges(6,1) * t279 - Icges(6,4) * t278 + Icges(6,5) * t547;
t542 = t358 * t368;
t280 = t355 * t542 - t356 * t357;
t281 = t355 * t356 + t357 * t542;
t82 = t178 * t543 - t182 * t280 + t186 * t281;
t179 = Icges(6,5) * t281 - Icges(6,6) * t280 + Icges(6,3) * t543;
t183 = Icges(6,4) * t281 - Icges(6,2) * t280 + Icges(6,6) * t543;
t187 = Icges(6,1) * t281 - Icges(6,4) * t280 + Icges(6,5) * t543;
t83 = t179 * t543 - t183 * t280 + t187 * t281;
t434 = t356 * t82 + t358 * t83;
t176 = Icges(7,5) * t279 + Icges(7,6) * t547 + Icges(7,3) * t278;
t180 = Icges(7,4) * t279 + Icges(7,2) * t547 + Icges(7,6) * t278;
t184 = Icges(7,1) * t279 + Icges(7,4) * t547 + Icges(7,5) * t278;
t80 = t176 * t280 + t180 * t543 + t184 * t281;
t177 = Icges(7,5) * t281 + Icges(7,6) * t543 + Icges(7,3) * t280;
t181 = Icges(7,4) * t281 + Icges(7,2) * t543 + Icges(7,6) * t280;
t185 = Icges(7,1) * t281 + Icges(7,4) * t543 + Icges(7,5) * t280;
t81 = t177 * t280 + t181 * t543 + t185 * t281;
t435 = t356 * t80 + t358 * t81;
t623 = t432 + t434 + t435;
t84 = t208 * t547 + t210 * t293 + t212 * t294;
t85 = t209 * t547 + t211 * t293 + t213 * t294;
t433 = t356 * t84 + t358 * t85;
t78 = t178 * t547 - t182 * t278 + t186 * t279;
t79 = t179 * t547 - t183 * t278 + t187 * t279;
t436 = t356 * t78 + t358 * t79;
t76 = t176 * t278 + t180 * t547 + t184 * t279;
t77 = t177 * t278 + t181 * t547 + t185 * t279;
t437 = t356 * t76 + t358 * t77;
t622 = t433 + t436 + t437;
t621 = -t513 / 0.2e1;
t352 = pkin(4) * t367 + pkin(3);
t578 = pkin(3) - t352;
t464 = t578 * t368;
t581 = pkin(8) * t365;
t618 = t464 + t581;
t505 = qJD(4) * t365;
t467 = t357 * t505;
t508 = qJD(3) * t368;
t475 = t355 * t508;
t617 = t467 + t475;
t473 = t356 * t508;
t386 = t358 * t513 + t473;
t363 = -qJ(5) - pkin(8);
t504 = qJD(4) * t367;
t616 = pkin(4) * t504 + t363 * t513;
t411 = -t211 * t364 + t213 * t367;
t100 = -t209 * t368 + t365 * t411;
t415 = t177 * t355 + t185 * t357;
t89 = -t181 * t368 + t365 * t415;
t413 = -t183 * t355 + t187 * t357;
t91 = -t179 * t368 + t365 * t413;
t497 = t100 + t89 + t91;
t416 = t176 * t355 + t184 * t357;
t88 = -t180 * t368 + t365 * t416;
t414 = -t182 * t355 + t186 * t357;
t90 = -t178 * t368 + t365 * t414;
t412 = -t210 * t364 + t212 * t367;
t99 = -t208 * t368 + t365 * t412;
t498 = t88 + t90 + t99;
t615 = t356 * t498 + t358 * t497;
t563 = Icges(7,5) * t357;
t417 = Icges(7,3) * t355 + t563;
t283 = -Icges(7,6) * t368 + t365 * t417;
t566 = Icges(6,4) * t357;
t422 = -Icges(6,2) * t355 + t566;
t286 = -Icges(6,6) * t368 + t365 * t422;
t564 = Icges(7,5) * t355;
t426 = Icges(7,1) * t357 + t564;
t287 = -Icges(7,4) * t368 + t365 * t426;
t567 = Icges(6,4) * t355;
t427 = Icges(6,1) * t357 - t567;
t288 = -Icges(6,5) * t368 + t365 * t427;
t568 = Icges(5,4) * t367;
t423 = -Icges(5,2) * t364 + t568;
t298 = -Icges(5,6) * t368 + t365 * t423;
t569 = Icges(5,4) * t364;
t428 = Icges(5,1) * t367 - t569;
t299 = -Icges(5,5) * t368 + t365 * t428;
t550 = t299 * t367;
t614 = t626 * t368 + (t298 * t364 - t550 + (-t287 - t288) * t357 + (-t283 + t286) * t355) * t365;
t596 = 2 * m(4);
t576 = rSges(4,1) * t368;
t445 = -rSges(4,2) * t365 + t576;
t574 = rSges(4,3) * t358;
t271 = t356 * t445 - t574;
t607 = -rSges(4,2) * t543 + t356 * rSges(4,3);
t272 = rSges(4,1) * t542 + t607;
t336 = rSges(4,1) * t365 + rSges(4,2) * t368;
t395 = qJD(3) * t336;
t479 = t356 * t513;
t374 = rSges(4,2) * t479 + rSges(4,3) * t514 - t358 * t395;
t98 = (qJD(1) * t271 + t374) * t358 + (-t356 * t395 + (-t272 + t607) * qJD(1)) * t356;
t613 = t596 * t98;
t571 = Icges(4,4) * t365;
t430 = Icges(4,1) * t368 - t571;
t266 = Icges(4,5) * t356 + t358 * t430;
t551 = t266 * t368;
t570 = Icges(4,4) * t368;
t425 = -Icges(4,2) * t365 + t570;
t264 = Icges(4,6) * t356 + t358 * t425;
t556 = t264 * t365;
t405 = -t551 + t556;
t611 = t358 * t405;
t609 = t365 * t578;
t466 = t357 * t503;
t471 = t358 * t509;
t172 = qJD(1) * t278 + t355 * t471 - t356 * t507 - t358 * t466;
t402 = t457 * t356;
t173 = t358 * t403 + (-t402 - t471) * t357;
t470 = t358 * t508;
t608 = rSges(7,2) * t470 + qJD(6) * t280 - t619 * t172 + t173 * t620;
t359 = cos(qJ(1)) * pkin(1);
t606 = t356 * pkin(7) + t359;
t350 = t358 * pkin(7);
t343 = pkin(4) * t544;
t541 = t363 * t365;
t518 = -t356 * t541 - t343;
t605 = t350 - t518;
t107 = Icges(7,5) * t175 + Icges(7,6) * t386 + Icges(7,3) * t174;
t111 = Icges(7,4) * t175 + Icges(7,2) * t386 + Icges(7,6) * t174;
t115 = Icges(7,1) * t175 + Icges(7,4) * t386 + Icges(7,5) * t174;
t385 = t470 - t479;
t19 = t107 * t280 + t111 * t543 + t115 * t281 - t172 * t176 + t173 * t184 + t180 * t385;
t106 = Icges(7,5) * t173 + Icges(7,6) * t385 - Icges(7,3) * t172;
t110 = Icges(7,4) * t173 + Icges(7,2) * t385 - Icges(7,6) * t172;
t114 = Icges(7,1) * t173 + Icges(7,4) * t385 - Icges(7,5) * t172;
t20 = t106 * t280 + t110 * t543 + t114 * t281 - t172 * t177 + t173 * t185 + t181 * t385;
t109 = Icges(6,5) * t175 - Icges(6,6) * t174 + Icges(6,3) * t386;
t113 = Icges(6,4) * t175 - Icges(6,2) * t174 + Icges(6,6) * t386;
t117 = Icges(6,1) * t175 - Icges(6,4) * t174 + Icges(6,5) * t386;
t21 = t109 * t543 - t113 * t280 + t117 * t281 + t172 * t182 + t173 * t186 + t178 * t385;
t108 = Icges(6,5) * t173 + Icges(6,6) * t172 + Icges(6,3) * t385;
t112 = Icges(6,4) * t173 + Icges(6,2) * t172 + Icges(6,6) * t385;
t116 = Icges(6,1) * t173 + Icges(6,4) * t172 + Icges(6,5) * t385;
t22 = t108 * t543 - t112 * t280 + t116 * t281 + t172 * t183 + t173 * t187 + t179 * t385;
t378 = t364 * t509 + t367 * t458;
t205 = t356 * t378 - t457 * t544;
t377 = t364 * t458 - t367 * t509;
t206 = t356 * t377 + t367 * t404;
t129 = Icges(5,5) * t206 + Icges(5,6) * t205 + Icges(5,3) * t386;
t131 = Icges(5,4) * t206 + Icges(5,2) * t205 + Icges(5,6) * t386;
t133 = Icges(5,1) * t206 + Icges(5,4) * t205 + Icges(5,5) * t386;
t203 = t358 * t378 + t457 * t548;
t204 = t358 * t377 - t367 * t402;
t31 = t129 * t543 + t131 * t295 + t133 * t296 + t203 * t210 + t204 * t212 + t208 * t385;
t128 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t385;
t130 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t385;
t132 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t385;
t32 = t128 * t543 + t130 * t295 + t132 * t296 + t203 * t211 + t204 * t213 + t209 * t385;
t604 = (-t19 - t21 - t31) * t358 + (t20 + t22 + t32) * t356 + t623 * qJD(1);
t23 = t107 * t278 + t111 * t547 + t115 * t279 + t174 * t176 + t175 * t184 + t180 * t386;
t24 = t106 * t278 + t110 * t547 + t114 * t279 + t174 * t177 + t175 * t185 + t181 * t386;
t25 = t109 * t547 - t113 * t278 + t117 * t279 - t174 * t182 + t175 * t186 + t178 * t386;
t26 = t108 * t547 - t112 * t278 + t116 * t279 - t174 * t183 + t175 * t187 + t179 * t386;
t33 = t129 * t547 + t131 * t293 + t133 * t294 + t205 * t210 + t206 * t212 + t208 * t386;
t34 = t128 * t547 + t130 * t293 + t132 * t294 + t205 * t211 + t206 * t213 + t209 * t386;
t603 = (-t23 - t25 - t33) * t358 + (t24 + t26 + t34) * t356 + t622 * qJD(1);
t27 = (qJD(3) * t416 - t111) * t368 + (qJD(3) * t180 + t107 * t355 + t115 * t357 + (t176 * t357 - t184 * t355) * qJD(4)) * t365;
t29 = (qJD(3) * t414 - t109) * t368 + (qJD(3) * t178 - t113 * t355 + t117 * t357 + (-t182 * t357 - t186 * t355) * qJD(4)) * t365;
t35 = (qJD(3) * t412 - t129) * t368 + (qJD(3) * t208 - t131 * t364 + t133 * t367 + (-t210 * t367 - t212 * t364) * qJD(4)) * t365;
t602 = -t27 - t29 - t35;
t28 = (qJD(3) * t415 - t110) * t368 + (qJD(3) * t181 + t106 * t355 + t114 * t357 + (t177 * t357 - t185 * t355) * qJD(4)) * t365;
t30 = (qJD(3) * t413 - t108) * t368 + (qJD(3) * t179 - t112 * t355 + t116 * t357 + (-t183 * t357 - t187 * t355) * qJD(4)) * t365;
t36 = (qJD(3) * t411 - t128) * t368 + (qJD(3) * t209 - t130 * t364 + t132 * t367 + (-t211 * t367 - t213 * t364) * qJD(4)) * t365;
t601 = t28 + t30 + t36;
t136 = t278 * t283 + t279 * t287 + t285 * t547;
t137 = -t278 * t286 + t279 * t288 + t284 * t547;
t143 = t293 * t298 + t294 * t299 + t297 * t547;
t600 = (-t136 - t137 - t143) * t368 + t622 * t365;
t138 = t280 * t283 + t281 * t287 + t285 * t543;
t139 = -t280 * t286 + t281 * t288 + t284 * t543;
t144 = t295 * t298 + t296 * t299 + t297 * t543;
t599 = (-t138 - t139 - t144) * t368 + t623 * t365;
t228 = (Icges(7,3) * t357 - t564) * t505 + (Icges(7,6) * t365 + t368 * t417) * qJD(3);
t229 = (-Icges(6,5) * t355 - Icges(6,6) * t357) * t505 + (Icges(6,3) * t365 + t368 * t418) * qJD(3);
t230 = (-Icges(7,4) * t355 + Icges(7,6) * t357) * t505 + (Icges(7,2) * t365 + t368 * t421) * qJD(3);
t232 = (-Icges(7,1) * t355 + t563) * t505 + (Icges(7,4) * t365 + t368 * t426) * qJD(3);
t233 = (-Icges(6,1) * t355 - t566) * t505 + (Icges(6,5) * t365 + t368 * t427) * qJD(3);
t242 = (-Icges(5,5) * t364 - Icges(5,6) * t367) * t505 + (Icges(5,3) * t365 + t368 * t419) * qJD(3);
t244 = (-Icges(5,1) * t364 - t568) * t505 + (Icges(5,5) * t365 + t368 * t428) * qJD(3);
t468 = t355 * t505;
t472 = t357 * t508;
t549 = t355 * t365;
t598 = t508 * t550 - t286 * t467 + t288 * t472 + t228 * t549 + (-t468 + t472) * t287 + t617 * t283 + t626 * t509 + (-t242 - t229 - t230) * t368 + (t244 * t367 - t298 * t504 + (t233 + t232) * t357) * t365;
t420 = Icges(4,5) * t368 - Icges(4,6) * t365;
t261 = -Icges(4,3) * t358 + t356 * t420;
t263 = -Icges(4,6) * t358 + t356 * t425;
t265 = -Icges(4,5) * t358 + t356 * t430;
t597 = t368 * t497 - t599;
t595 = 2 * m(5);
t594 = 2 * m(6);
t593 = 2 * m(7);
t592 = m(6) / 0.2e1;
t591 = m(7) / 0.2e1;
t590 = t356 / 0.2e1;
t589 = t358 / 0.2e1;
t588 = -t368 / 0.2e1;
t586 = -rSges(5,3) - pkin(8);
t585 = m(4) * t336;
t584 = sin(qJ(1)) * pkin(1);
t583 = pkin(3) * t368;
t582 = pkin(4) * t364;
t577 = pkin(8) + t363;
t575 = rSges(7,2) * t365;
t573 = rSges(6,3) * t365;
t440 = -rSges(6,1) * t279 + rSges(6,2) * t278;
t190 = rSges(6,3) * t547 - t440;
t560 = t190 * t358;
t443 = -rSges(5,1) * t294 - rSges(5,2) * t293;
t217 = rSges(5,3) * t547 - t443;
t559 = t217 * t358;
t558 = t263 * t365;
t557 = t263 * t368;
t555 = t264 * t368;
t554 = t265 * t365;
t553 = t265 * t368;
t552 = t266 * t365;
t540 = t363 * t368;
t537 = -rSges(7,2) * t479 + t608;
t536 = rSges(7,2) * t386 - t624;
t492 = t173 * rSges(6,1) + t172 * rSges(6,2) + rSges(6,3) * t470;
t119 = -rSges(6,3) * t479 + t492;
t329 = pkin(8) * t470;
t502 = qJD(5) * t365;
t452 = qJD(1) * t343 + t356 * t616 + t358 * t502;
t459 = t503 * t582;
t126 = -t329 + t618 * t515 + (-t459 + (-t540 + t609) * qJD(3)) * t358 + t452;
t535 = -t119 - t126;
t327 = pkin(3) * t474;
t342 = pkin(4) * t548;
t431 = t352 * t474 + t356 * t459 + t358 * t616 + t363 * t473;
t127 = t327 + (-pkin(8) * t508 + t502) * t356 + (-t358 * t618 + t342) * qJD(1) - t431;
t215 = -t356 * t618 + t518;
t534 = t127 * t543 + t215 * t470;
t533 = rSges(7,2) * t547 + t625;
t532 = -t190 - t215;
t531 = rSges(7,2) * t543 + t619 * t280 + t281 * t620;
t192 = t281 * rSges(6,1) - t280 * rSges(6,2) + rSges(6,3) * t543;
t345 = pkin(3) * t542;
t303 = pkin(8) * t543 + t345;
t400 = t352 * t542 - t358 * t541 + t342;
t216 = t400 - t303;
t530 = -t192 - t216;
t277 = t368 * t577 - t609;
t529 = t216 * t509 + t277 * t479;
t528 = t368 * t215 + t277 * t547;
t438 = rSges(7,1) * t357 + rSges(7,3) * t355;
t527 = (pkin(5) * t508 + qJ(6) * t505) * t357 + (qJ(6) * t508 + (-pkin(5) * qJD(4) + qJD(6)) * t365) * t355 + (-rSges(7,1) * t355 + rSges(7,3) * t357) * t505 + (t368 * t438 + t575) * qJD(3);
t218 = t296 * rSges(5,1) + t295 * rSges(5,2) + rSges(5,3) * t543;
t526 = -t218 - t303;
t439 = rSges(6,1) * t357 - rSges(6,2) * t355;
t235 = (-rSges(6,1) * t355 - rSges(6,2) * t357) * t505 + (t368 * t439 + t573) * qJD(3);
t465 = t364 * t505;
t241 = -pkin(4) * t465 - qJD(5) * t368 + (-t365 * t577 - t464) * qJD(3);
t525 = -t235 - t241;
t442 = rSges(5,1) * t367 - rSges(5,2) * t364;
t245 = (-rSges(5,1) * t364 - rSges(5,2) * t367) * t505 + (rSges(5,3) * t365 + t368 * t442) * qJD(3);
t448 = t581 + t583;
t317 = t448 * qJD(3);
t524 = -t245 - t317;
t341 = pkin(3) * t365 - pkin(8) * t368;
t304 = t341 * t515;
t523 = t277 * t515 + t304;
t302 = t448 * t356;
t522 = t356 * t302 + t358 * t303;
t290 = -rSges(6,3) * t368 + t365 * t439;
t521 = -t277 - t290;
t520 = -rSges(7,2) * t368 + (pkin(5) * t357 + qJ(6) * t355 + t438) * t365;
t301 = -rSges(5,3) * t368 + t365 * t442;
t519 = -t301 - t341;
t517 = t356 ^ 2 + t358 ^ 2;
t262 = Icges(4,3) * t356 + t358 * t420;
t516 = qJD(1) * t262;
t511 = qJD(3) * t356;
t510 = qJD(3) * t358;
t231 = (-Icges(6,2) * t357 - t567) * t505 + (Icges(6,6) * t365 + t368 * t422) * qJD(3);
t243 = (-Icges(5,2) * t367 - t569) * t505 + (Icges(5,6) * t365 + t368 * t423) * qJD(3);
t476 = t298 * t508;
t499 = (-t286 * t508 + (-qJD(4) * t288 - t231) * t365) * t355 + (-t476 + (-qJD(4) * t299 - t243) * t365) * t364 + t598;
t494 = -t126 - t537;
t493 = t614 * t509;
t490 = -t215 - t533;
t489 = -t216 - t531;
t488 = t204 * rSges(5,1) + t203 * rSges(5,2) + rSges(5,3) * t470;
t487 = -t241 - t527;
t486 = -t317 + t525;
t485 = t356 * (pkin(8) * t386 + qJD(1) * t345 - t327) + t358 * (-pkin(8) * t479 + t329 + (-t356 * t512 - t471) * pkin(3)) + t302 * t514;
t484 = -t277 - t520;
t483 = -t341 + t521;
t482 = t358 * pkin(2) + t606;
t481 = t290 * t515;
t480 = t301 * t515;
t469 = t365 * t508;
t463 = -t352 * t368 - pkin(2);
t462 = t533 * t358;
t461 = t520 * t356;
t247 = t519 * t358;
t460 = qJD(1) * t520;
t456 = t368 * t127 + t241 * t547 + t277 * t386;
t455 = t356 * t215 + t358 * t216 + t522;
t454 = -t317 + t487;
t453 = -t341 + t484;
t47 = t356 * t77 - t358 * t76;
t48 = t356 * t79 - t358 * t78;
t57 = t356 * t85 - t358 * t84;
t451 = t47 / 0.2e1 + t48 / 0.2e1 + t57 / 0.2e1;
t49 = t356 * t81 - t358 * t80;
t50 = t356 * t83 - t358 * t82;
t58 = t356 * t87 - t358 * t86;
t450 = t49 / 0.2e1 + t50 / 0.2e1 + t58 / 0.2e1;
t449 = t484 * t358;
t164 = t483 * t358;
t447 = t356 * t460;
t444 = t206 * rSges(5,1) + t205 * rSges(5,2);
t441 = t175 * rSges(6,1) - t174 * rSges(6,2);
t424 = Icges(4,2) * t368 + t571;
t410 = -t218 * t356 + t559;
t409 = -t217 * t356 - t218 * t358;
t406 = -t553 + t558;
t152 = t453 * t358;
t401 = -pkin(2) - t445;
t399 = t463 - t575;
t398 = t463 - t573;
t397 = t358 * t126 + t356 * t127 + t215 * t514 + t485;
t394 = t356 * t530 + t560;
t393 = t365 * t586 - pkin(2) - t583;
t391 = qJD(3) * t424;
t390 = qJD(3) * (-Icges(4,5) * t365 - Icges(4,6) * t368);
t388 = -t368 * t498 + t600;
t387 = -t359 + (-pkin(7) - t582) * t356;
t51 = -t172 * t283 + t173 * t287 + t228 * t280 + t230 * t543 + t232 * t281 + t285 * t385;
t52 = t172 * t286 + t173 * t288 + t229 * t543 - t231 * t280 + t233 * t281 + t284 * t385;
t61 = t203 * t298 + t204 * t299 + t242 * t543 + t243 * t295 + t244 * t296 + t297 * t385;
t384 = t51 / 0.2e1 + t52 / 0.2e1 + t61 / 0.2e1 + t30 / 0.2e1 + t36 / 0.2e1 + t28 / 0.2e1;
t53 = t174 * t283 + t175 * t287 + t228 * t278 + t230 * t547 + t232 * t279 + t285 * t386;
t54 = -t174 * t286 + t175 * t288 + t229 * t547 - t231 * t278 + t233 * t279 + t284 * t386;
t62 = t205 * t298 + t206 * t299 + t242 * t547 + t243 * t293 + t244 * t294 + t297 * t386;
t383 = t53 / 0.2e1 + t54 / 0.2e1 + t62 / 0.2e1 + t29 / 0.2e1 + t35 / 0.2e1 + t27 / 0.2e1;
t382 = t400 + t482;
t381 = t90 / 0.2e1 + t99 / 0.2e1 + t88 / 0.2e1 + t136 / 0.2e1 + t137 / 0.2e1 + t143 / 0.2e1;
t380 = t91 / 0.2e1 + t100 / 0.2e1 + t89 / 0.2e1 + t138 / 0.2e1 + t139 / 0.2e1 + t144 / 0.2e1;
t379 = t356 * t489 + t462;
t376 = t356 * t399 - t584;
t375 = t356 * t398 - t584;
t371 = t356 * t393 - t584;
t347 = pkin(7) * t514;
t370 = t347 + (-t459 + (-t352 * t365 - t540) * qJD(3)) * t358 + t452;
t362 = t365 ^ 2;
t316 = t445 * qJD(3);
t310 = t420 * qJD(3);
t246 = t519 * t356;
t240 = t272 + t482;
t239 = t356 * t401 + t350 + t574 - t584;
t222 = t356 * t390 + t516;
t221 = -qJD(1) * t261 + t358 * t390;
t193 = t215 * t543;
t166 = t336 * t511 + (-t359 + (-rSges(4,3) - pkin(7)) * t356 + t401 * t358) * qJD(1);
t165 = t347 + (-t584 + (-pkin(2) - t576) * t356) * qJD(1) + t374;
t163 = t483 * t356;
t160 = -t218 * t368 - t301 * t543;
t159 = t217 * t368 + t301 * t547;
t158 = t482 - t526;
t157 = t350 + t371 + t443;
t151 = t453 * t356;
t150 = t262 * t356 - t611;
t149 = t261 * t356 - t358 * t406;
t148 = -t262 * t358 - t356 * t405;
t147 = -t261 * t358 - t356 * t406;
t146 = qJD(1) * t247 + t356 * t524;
t145 = t358 * t524 + t304 + t480;
t142 = t382 + t192;
t141 = t375 + t440 + t605;
t140 = t410 * t365;
t135 = rSges(5,3) * t386 + t444;
t134 = -rSges(5,3) * t479 + t488;
t121 = rSges(6,3) * t386 + t441;
t105 = t382 + t531;
t104 = t376 + t605 - t625;
t101 = -t409 + t522;
t97 = t368 * t530 + t521 * t543;
t96 = t190 * t368 + t290 * t547 + t528;
t95 = qJD(1) * t164 + t356 * t486;
t94 = t358 * t486 + t481 + t523;
t93 = t327 + t586 * t473 + (t358 * t393 - t606) * qJD(1) - t444;
t92 = -pkin(3) * t471 + qJD(1) * t371 + t329 + t347 + t488;
t75 = t365 * t394 + t193;
t74 = t365 * t449 + t368 * t489;
t73 = t365 * t461 + t368 * t533 + t528;
t71 = (-rSges(6,3) * t508 - t502) * t356 + (t358 * t398 + t387) * qJD(1) + t431 - t441;
t70 = qJD(1) * t375 + t370 + t492;
t69 = qJD(1) * t152 + t356 * t454;
t68 = t358 * t454 + t447 + t523;
t67 = t190 * t356 + t192 * t358 + t455;
t66 = (t301 * t511 + t135) * t368 + (-qJD(3) * t217 + t245 * t356 + t301 * t514) * t365;
t65 = (-t301 * t510 - t134) * t368 + (qJD(3) * t218 - t245 * t358 + t480) * t365;
t60 = t365 * t379 + t193;
t59 = t356 * t533 + t358 * t531 + t455;
t56 = (-rSges(7,2) * t508 - t502) * t356 + (t358 * t399 + t387) * qJD(1) + t431 + t624;
t55 = qJD(1) * t376 + t370 + t608;
t46 = t410 * t508 + (qJD(1) * t409 - t134 * t356 + t135 * t358) * t365;
t43 = t134 * t358 + t135 * t356 + (t356 * t526 + t559) * qJD(1) + t485;
t38 = (t290 * t511 + t121) * t368 + (qJD(3) * t532 + t235 * t356 + t290 * t514) * t365 + t456;
t37 = t535 * t368 + (qJD(3) * t192 + t481) * t365 + (t365 * t525 + t508 * t521) * t358 + t529;
t18 = (qJD(3) * t461 + t536) * t368 + (qJD(3) * t490 + t356 * t527 + t358 * t460) * t365 + t456;
t17 = (qJD(3) * t449 + t494) * t368 + (qJD(3) * t531 + t358 * t487 + t447) * t365 + t529;
t16 = t119 * t358 + t121 * t356 + (t560 + (-t303 + t530) * t356) * qJD(1) + t397;
t15 = t394 * t508 + (t121 * t358 + t535 * t356 + (t356 * t532 + t358 * t530) * qJD(1)) * t365 + t534;
t14 = t537 * t358 + t536 * t356 + (t462 + (-t303 + t489) * t356) * qJD(1) + t397;
t13 = t379 * t508 + (t536 * t358 + t494 * t356 + (t356 * t490 + t358 * t489) * qJD(1)) * t365 + t534;
t6 = (qJD(3) * t433 - t62) * t368 + (-qJD(1) * t57 + qJD(3) * t143 + t33 * t356 + t34 * t358) * t365;
t5 = (qJD(3) * t432 - t61) * t368 + (-qJD(1) * t58 + qJD(3) * t144 + t31 * t356 + t32 * t358) * t365;
t4 = (qJD(3) * t436 - t54) * t368 + (-qJD(1) * t48 + qJD(3) * t137 + t25 * t356 + t26 * t358) * t365;
t3 = (qJD(3) * t437 - t53) * t368 + (-qJD(1) * t47 + qJD(3) * t136 + t23 * t356 + t24 * t358) * t365;
t2 = (qJD(3) * t434 - t52) * t368 + (-qJD(1) * t50 + qJD(3) * t139 + t21 * t356 + t22 * t358) * t365;
t1 = (qJD(3) * t435 - t51) * t368 + (-qJD(1) * t49 + qJD(3) * t138 + t19 * t356 + t20 * t358) * t365;
t7 = [(t104 * t56 + t105 * t55) * t593 + (t141 * t71 + t142 * t70) * t594 + (t157 * t93 + t158 * t92) * t595 + (t165 * t240 + t166 * t239) * t596 - t286 * t475 - t299 * t465 - t288 * t468 - t231 * t549 + t598 + (-t424 + t430) * t509 + (Icges(4,1) * t365 + t425 + t570) * t508 + (-t243 * t365 - t476) * t364; 0; 0; m(7) * (t104 * t68 + t105 * t69 + t151 * t55 + t152 * t56) + m(6) * (t141 * t94 + t142 * t95 + t163 * t70 + t164 * t71) + m(5) * (t145 * t157 + t146 * t158 + t246 * t92 + t247 * t93) + ((qJD(1) * t264 - t356 * t391) * t588 + t266 * t621 + m(4) * (-t166 * t336 - t239 * t316) + t310 * t589 + (t558 / 0.2e1 - t553 / 0.2e1) * qJD(3) - t383) * t358 + ((-qJD(1) * t263 - t358 * t391) * t368 / 0.2e1 + t265 * t621 + m(4) * (-t165 * t336 - t240 * t316) + t310 * t590 + (-t556 / 0.2e1 + t551 / 0.2e1) * qJD(3) + t384) * t356 + ((t555 / 0.2e1 + t552 / 0.2e1 - t240 * t585 + t380) * t358 + (t557 / 0.2e1 + t554 / 0.2e1 + t239 * t585 + t381) * t356) * qJD(1); m(4) * t98 + m(5) * t43 + m(6) * t16 + m(7) * t14; (t16 * t67 + t163 * t95 + t164 * t94) * t594 + (t14 * t59 + t151 * t69 + t152 * t68) * t593 + (t101 * t43 + t145 * t247 + t146 * t246) * t595 + t517 * t336 * t316 * t596 + (t272 * t613 + (-t148 * qJD(1) + (-qJD(1) * t406 - t222) * t358) * t358 - t603) * t358 + (t271 * t613 + (t149 * qJD(1) + (t405 * qJD(1) + t221) * t356) * t356 + ((t264 * t508 + t266 * t509 - t516 - t222 + (-t552 - t555) * qJD(3)) * t356 + (t221 - (t554 + t557) * qJD(3) + t263 * t508 + t265 * t509) * t358 + (-t147 + t150 + (t262 - t406) * t356 + t611) * qJD(1)) * t358 + t604) * t356 + (-t147 * t358 + t148 * t356 + t47 + t48 + t57) * t515 + (-t149 * t358 + t150 * t356 + t49 + t50 + t58) * t514; m(7) * (t104 * t18 + t105 * t17 + t55 * t74 + t56 * t73) + m(6) * (t141 * t38 + t142 * t37 + t70 * t97 + t71 * t96) + m(5) * (t157 * t66 + t158 * t65 + t159 * t93 + t160 * t92) + ((t356 * t381 + t358 * t380) * qJD(3) - t499) * t368 + (t384 * t358 + t383 * t356 + (-t356 * t380 + t358 * t381) * qJD(1)) * t365 - t493; m(5) * t46 + m(6) * t15 + m(7) * t13; m(7) * (t13 * t59 + t14 * t60 + t151 * t17 + t152 * t18 + t68 * t73 + t69 * t74) + m(6) * (t15 * t67 + t16 * t75 + t163 * t37 + t164 * t38 + t94 * t96 + t95 * t97) + m(5) * (t101 * t46 + t140 * t43 + t145 * t159 + t146 * t160 + t246 * t65 + t247 * t66) + (-t6 / 0.2e1 - t4 / 0.2e1 - t3 / 0.2e1 + t450 * t508) * t358 + (t5 / 0.2e1 + t2 / 0.2e1 + t1 / 0.2e1 + t451 * t508) * t356 + (-t356 * t450 + t358 * t451) * t513 + (t356 * t497 - t358 * t498) * t509 / 0.2e1 + (qJD(1) * t615 + t601 * t356 + t602 * t358) * t588 + (t600 * t356 + t599 * t358) * qJD(1) / 0.2e1 + (t589 * t604 + t590 * t603) * t365; (t13 * t60 + t17 * t74 + t18 * t73) * t593 + (t15 * t75 + t37 * t97 + t38 * t96) * t594 + (t140 * t46 + t159 * t66 + t160 * t65) * t595 + (t499 * t368 + (t388 * t356 - t358 * t597) * qJD(3) + t493) * t368 + ((-t368 * t601 + t1 + t2 + t5) * t358 + (t368 * t602 + t3 + t4 + t6) * t356 + (t365 * t615 + t368 * t614) * qJD(3) + (t356 * t597 + t388 * t358) * qJD(1)) * t365; 0.2e1 * ((t141 * t358 + t142 * t356) * t592 + (t104 * t358 + t105 * t356) * t591) * t508 + 0.2e1 * ((-t141 * t515 + t142 * t514 + t356 * t70 + t358 * t71) * t592 + (-t104 * t515 + t105 * t514 + t356 * t55 + t358 * t56) * t591) * t365; (m(6) + m(7)) * t509; 0.2e1 * ((t163 * t511 + t164 * t510 - t16) * t592 + (t151 * t511 + t152 * t510 - t14) * t591) * t368 + 0.2e1 * ((qJD(3) * t67 + t163 * t514 - t164 * t515 + t356 * t95 + t358 * t94) * t592 + (qJD(3) * t59 + t151 * t514 - t152 * t515 + t356 * t69 + t358 * t68) * t591) * t365; 0.2e1 * ((t510 * t73 + t511 * t74 - t13) * t591 + (t510 * t96 + t511 * t97 - t15) * t592) * t368 + 0.2e1 * ((qJD(3) * t60 + t17 * t356 + t18 * t358 + t514 * t74 - t515 * t73) * t591 + (qJD(3) * t75 + t356 * t37 + t358 * t38 + t514 * t97 - t515 * t96) * t592) * t365; 0.4e1 * (t592 + t591) * (-0.1e1 + t517) * t469; m(7) * (-t104 * t172 + t105 * t174 + t278 * t55 + t280 * t56); t617 * m(7); m(7) * (t59 * t467 + t151 * t174 - t152 * t172 + t278 * t69 + t280 * t68 + (t14 * t365 + t508 * t59) * t355); m(7) * (t60 * t467 + t17 * t278 - t172 * t73 + t174 * t74 + t18 * t280 + (t13 * t365 + t508 * t60) * t355); m(7) * ((t355 * t362 + (t278 * t356 + t280 * t358 - t355 * t368) * t368) * qJD(3) + (-t466 - t172 * t358 + t174 * t356 + (t278 * t358 - t280 * t356) * qJD(1)) * t365); (-t172 * t280 + t174 * t278 + (t355 * t469 + t362 * t506) * t355) * t593;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;