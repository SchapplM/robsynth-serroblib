% Calculate time derivative of joint inertia matrix for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP6_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:44
% EndTime: 2019-03-09 00:29:49
% DurationCPUTime: 38.95s
% Computational Cost: add. (292402->1485), mult. (876078->2002), div. (0->0), fcn. (1057210->14), ass. (0->572)
t678 = rSges(7,1) + pkin(5);
t677 = rSges(7,3) + qJ(6);
t494 = sin(pkin(6));
t676 = t494 / 0.2e1;
t496 = cos(pkin(6));
t675 = t496 / 0.2e1;
t492 = sin(pkin(12));
t495 = cos(pkin(12));
t500 = sin(qJ(2));
t501 = cos(qJ(2));
t649 = t496 * t501;
t483 = -t492 * t500 + t495 * t649;
t474 = t483 * qJD(2);
t650 = t496 * t500;
t484 = t492 * t501 + t495 * t650;
t475 = t484 * qJD(2);
t499 = sin(qJ(3));
t661 = cos(pkin(7));
t665 = cos(qJ(3));
t529 = t661 * t665;
t493 = sin(pkin(7));
t576 = t493 * t665;
t544 = t494 * t576;
t502 = t483 * t529 - t484 * t499 - t495 * t544;
t571 = t499 * t661;
t404 = qJD(3) * t502 + t474 * t665 - t475 * t571;
t570 = t661 * t483;
t577 = t484 * t665;
t651 = t494 * t495;
t428 = t577 + (-t493 * t651 + t570) * t499;
t572 = t494 * t661;
t462 = -t483 * t493 - t495 * t572;
t498 = sin(qJ(4));
t664 = cos(qJ(4));
t411 = t428 * t664 + t462 * t498;
t575 = t493 * t664;
t337 = qJD(4) * t411 + t404 * t498 - t475 * t575;
t514 = -t428 * t498 + t462 * t664;
t653 = t493 * t498;
t338 = qJD(4) * t514 + t404 * t664 + t475 * t653;
t264 = pkin(4) * t338 + pkin(11) * t337;
t360 = pkin(4) * t411 - pkin(11) * t514;
t485 = -t492 * t649 - t495 * t500;
t518 = t492 * t650 - t495 * t501;
t654 = t493 * t494;
t430 = -t518 * t665 + (t485 * t661 + t492 * t654) * t499;
t476 = t485 * qJD(2);
t477 = t518 * qJD(2);
t405 = qJD(3) * t430 + t476 * t499 - t477 * t529;
t503 = t485 * t529 + t492 * t544 + t499 * t518;
t640 = -t264 * t503 + t405 * t360;
t657 = t475 * t493;
t445 = pkin(2) * t474 + pkin(9) * t657;
t656 = t477 * t493;
t446 = pkin(2) * t476 - pkin(9) * t656;
t655 = t492 * t494;
t619 = t445 * t655 + t446 * t651;
t674 = -t499 * t500 + t501 * t529;
t669 = m(7) / 0.2e1;
t673 = m(6) / 0.2e1 + t669;
t543 = t496 * t576;
t422 = qJD(3) * t543 + (t674 * qJD(3) + (-t500 * t571 + t501 * t665) * qJD(2)) * t494;
t507 = t500 * t665 + t501 * t571;
t652 = t493 * t499;
t461 = t494 * t507 + t496 * t652;
t482 = t496 * t661 - t501 * t654;
t512 = -t461 * t498 + t482 * t664;
t618 = qJD(2) * t494;
t574 = t500 * t618;
t542 = t493 * t574;
t391 = qJD(4) * t512 + t422 * t664 + t498 * t542;
t432 = t461 * t664 + t482 * t498;
t460 = -t494 * t674 - t543;
t497 = sin(qJ(5));
t663 = cos(qJ(5));
t409 = t432 * t663 + t460 * t497;
t573 = qJD(3) * t652;
t421 = t496 * t573 + (t507 * qJD(3) + (t499 * t501 + t500 * t529) * qJD(2)) * t494;
t302 = t409 * qJD(5) + t391 * t497 - t421 * t663;
t672 = -0.2e1 * t462;
t671 = m(5) / 0.2e1;
t668 = -t492 / 0.2e1;
t667 = t495 / 0.2e1;
t666 = -t496 / 0.2e1;
t406 = qJD(3) * t503 + t476 * t665 + t477 * t571;
t335 = rSges(4,1) * t406 - rSges(4,2) * t405 - rSges(4,3) * t656;
t662 = m(4) * t335;
t660 = Icges(3,4) * t500;
t659 = Icges(3,4) * t501;
t364 = t411 * t663 - t497 * t502;
t403 = t474 * t499 + t475 * t529 - t573 * t651 + (t499 * t570 + t577) * qJD(3);
t248 = qJD(5) * t364 + t338 * t497 - t403 * t663;
t517 = -t411 * t497 - t502 * t663;
t249 = qJD(5) * t517 + t338 * t663 + t403 * t497;
t647 = rSges(7,2) * t337 - qJD(6) * t517 + t677 * t248 + t249 * t678;
t183 = rSges(6,1) * t249 - rSges(6,2) * t248 + rSges(6,3) * t337;
t646 = -t183 - t264;
t463 = -t485 * t493 + t492 * t572;
t513 = -t430 * t498 + t463 * t664;
t340 = qJD(4) * t513 + t406 * t664 - t477 * t653;
t413 = t430 * t664 + t463 * t498;
t366 = t413 * t663 - t497 * t503;
t250 = qJD(5) * t366 + t340 * t497 - t405 * t663;
t516 = -t413 * t497 - t503 * t663;
t251 = qJD(5) * t516 + t340 * t663 + t405 * t497;
t339 = qJD(4) * t413 + t406 * t498 + t477 * t575;
t645 = rSges(7,2) * t339 - qJD(6) * t516 + t677 * t250 + t251 * t678;
t185 = rSges(6,1) * t251 - rSges(6,2) * t250 + rSges(6,3) * t339;
t265 = pkin(4) * t340 + pkin(11) * t339;
t644 = -t185 - t265;
t515 = -t432 * t497 + t460 * t663;
t303 = qJD(5) * t515 + t391 * t663 + t421 * t497;
t390 = qJD(4) * t432 + t422 * t498 - t542 * t664;
t643 = rSges(7,2) * t390 - qJD(6) * t515 + t677 * t302 + t303 * t678;
t218 = rSges(6,1) * t303 - rSges(6,2) * t302 + rSges(6,3) * t390;
t306 = pkin(4) * t391 + pkin(11) * t390;
t642 = -t218 - t306;
t241 = rSges(5,1) * t338 - rSges(5,2) * t337 + rSges(5,3) * t403;
t349 = pkin(3) * t404 + pkin(10) * t403;
t641 = -t241 - t349;
t361 = pkin(4) * t413 - pkin(11) * t513;
t639 = t460 * t265 + t421 * t361;
t638 = -rSges(7,2) * t514 + t364 * t678 - t677 * t517;
t282 = rSges(6,1) * t364 + rSges(6,2) * t517 - rSges(6,3) * t514;
t637 = -t282 - t360;
t636 = -rSges(7,2) * t513 + t366 * t678 - t677 * t516;
t284 = rSges(6,1) * t366 + rSges(6,2) * t516 - rSges(6,3) * t513;
t635 = -t284 - t361;
t292 = rSges(5,1) * t391 - rSges(5,2) * t390 + rSges(5,3) * t421;
t394 = pkin(3) * t422 + pkin(10) * t421;
t634 = -t292 - t394;
t400 = pkin(4) * t432 - pkin(11) * t512;
t633 = -t306 * t502 + t403 * t400;
t367 = t462 * t394;
t632 = t462 * t306 + t367;
t309 = t463 * t349;
t399 = pkin(3) * t430 - pkin(10) * t503;
t384 = t399 * t657;
t631 = t309 - t384;
t323 = rSges(5,1) * t411 + rSges(5,2) * t514 - rSges(5,3) * t502;
t398 = pkin(3) * t428 - pkin(10) * t502;
t630 = -t323 - t398;
t629 = -rSges(7,2) * t512 + t409 * t678 - t677 * t515;
t326 = rSges(6,1) * t409 + rSges(6,2) * t515 - rSges(6,3) * t512;
t628 = -t326 - t400;
t350 = pkin(3) * t406 + pkin(10) * t405;
t440 = t496 * t446;
t627 = t496 * t350 + t440;
t380 = t463 * t398;
t626 = t463 * t360 + t380;
t392 = t482 * t399;
t625 = t482 * t361 + t392;
t379 = rSges(5,1) * t432 + rSges(5,2) * t512 + rSges(5,3) * t460;
t419 = pkin(3) * t461 + pkin(10) * t460;
t624 = -t379 - t419;
t407 = t462 * t419;
t623 = t462 * t400 + t407;
t439 = -pkin(2) * t518 + pkin(9) * t463;
t433 = t496 * t439;
t622 = t496 * t399 + t433;
t438 = t484 * pkin(2) + pkin(9) * t462;
t621 = t438 * t655 + t439 * t651;
t620 = 0.2e1 * t619;
t615 = -0.2e1 * t656;
t269 = Icges(7,5) * t364 - Icges(7,6) * t514 - Icges(7,3) * t517;
t273 = Icges(7,4) * t364 - Icges(7,2) * t514 - Icges(7,6) * t517;
t277 = Icges(7,1) * t364 - Icges(7,4) * t514 - Icges(7,5) * t517;
t139 = -t269 * t517 - t273 * t514 + t277 * t364;
t270 = Icges(7,5) * t366 - Icges(7,6) * t513 - Icges(7,3) * t516;
t274 = Icges(7,4) * t366 - Icges(7,2) * t513 - Icges(7,6) * t516;
t278 = Icges(7,1) * t366 - Icges(7,4) * t513 - Icges(7,5) * t516;
t140 = -t270 * t517 - t274 * t514 + t278 * t364;
t311 = Icges(7,5) * t409 - Icges(7,6) * t512 - Icges(7,3) * t515;
t315 = Icges(7,4) * t409 - Icges(7,2) * t512 - Icges(7,6) * t515;
t319 = Icges(7,1) * t409 - Icges(7,4) * t512 - Icges(7,5) * t515;
t188 = -t311 * t517 - t315 * t514 + t319 * t364;
t170 = Icges(7,5) * t249 + Icges(7,6) * t337 + Icges(7,3) * t248;
t174 = Icges(7,4) * t249 + Icges(7,2) * t337 + Icges(7,6) * t248;
t178 = Icges(7,1) * t249 + Icges(7,4) * t337 + Icges(7,5) * t248;
t41 = -t170 * t517 - t174 * t514 + t178 * t364 + t248 * t269 + t249 * t277 + t273 * t337;
t171 = Icges(7,5) * t251 + Icges(7,6) * t339 + Icges(7,3) * t250;
t175 = Icges(7,4) * t251 + Icges(7,2) * t339 + Icges(7,6) * t250;
t179 = Icges(7,1) * t251 + Icges(7,4) * t339 + Icges(7,5) * t250;
t42 = -t171 * t517 - t175 * t514 + t179 * t364 + t248 * t270 + t249 * t278 + t274 * t337;
t211 = Icges(7,5) * t303 + Icges(7,6) * t390 + Icges(7,3) * t302;
t213 = Icges(7,4) * t303 + Icges(7,2) * t390 + Icges(7,6) * t302;
t215 = Icges(7,1) * t303 + Icges(7,4) * t390 + Icges(7,5) * t302;
t67 = -t211 * t517 - t213 * t514 + t215 * t364 + t248 * t311 + t249 * t319 + t315 * t337;
t1 = t139 * t337 + t140 * t339 + t188 * t390 - t41 * t514 - t42 * t513 - t512 * t67;
t271 = Icges(6,5) * t364 + Icges(6,6) * t517 - Icges(6,3) * t514;
t275 = Icges(6,4) * t364 + Icges(6,2) * t517 - Icges(6,6) * t514;
t279 = Icges(6,1) * t364 + Icges(6,4) * t517 - Icges(6,5) * t514;
t141 = -t271 * t514 + t275 * t517 + t279 * t364;
t272 = Icges(6,5) * t366 + Icges(6,6) * t516 - Icges(6,3) * t513;
t276 = Icges(6,4) * t366 + Icges(6,2) * t516 - Icges(6,6) * t513;
t280 = Icges(6,1) * t366 + Icges(6,4) * t516 - Icges(6,5) * t513;
t142 = -t272 * t514 + t276 * t517 + t280 * t364;
t312 = Icges(6,5) * t409 + Icges(6,6) * t515 - Icges(6,3) * t512;
t316 = Icges(6,4) * t409 + Icges(6,2) * t515 - Icges(6,6) * t512;
t320 = Icges(6,1) * t409 + Icges(6,4) * t515 - Icges(6,5) * t512;
t189 = -t312 * t514 + t316 * t517 + t320 * t364;
t172 = Icges(6,5) * t249 - Icges(6,6) * t248 + Icges(6,3) * t337;
t176 = Icges(6,4) * t249 - Icges(6,2) * t248 + Icges(6,6) * t337;
t180 = Icges(6,1) * t249 - Icges(6,4) * t248 + Icges(6,5) * t337;
t43 = -t172 * t514 + t176 * t517 + t180 * t364 - t248 * t275 + t249 * t279 + t271 * t337;
t173 = Icges(6,5) * t251 - Icges(6,6) * t250 + Icges(6,3) * t339;
t177 = Icges(6,4) * t251 - Icges(6,2) * t250 + Icges(6,6) * t339;
t181 = Icges(6,1) * t251 - Icges(6,4) * t250 + Icges(6,5) * t339;
t44 = -t173 * t514 + t177 * t517 + t181 * t364 - t248 * t276 + t249 * t280 + t272 * t337;
t212 = Icges(6,5) * t303 - Icges(6,6) * t302 + Icges(6,3) * t390;
t214 = Icges(6,4) * t303 - Icges(6,2) * t302 + Icges(6,6) * t390;
t216 = Icges(6,1) * t303 - Icges(6,4) * t302 + Icges(6,5) * t390;
t68 = -t212 * t514 + t214 * t517 + t216 * t364 - t248 * t316 + t249 * t320 + t312 * t337;
t2 = t141 * t337 + t142 * t339 + t189 * t390 - t43 * t514 - t44 * t513 - t512 * t68;
t614 = t2 / 0.2e1 + t1 / 0.2e1;
t143 = -t269 * t516 - t273 * t513 + t277 * t366;
t144 = -t270 * t516 - t274 * t513 + t278 * t366;
t190 = -t311 * t516 - t315 * t513 + t319 * t366;
t45 = -t170 * t516 - t174 * t513 + t178 * t366 + t250 * t269 + t251 * t277 + t273 * t339;
t46 = -t171 * t516 - t175 * t513 + t179 * t366 + t250 * t270 + t251 * t278 + t274 * t339;
t69 = -t211 * t516 - t213 * t513 + t215 * t366 + t250 * t311 + t251 * t319 + t315 * t339;
t3 = t143 * t337 + t144 * t339 + t190 * t390 - t45 * t514 - t46 * t513 - t512 * t69;
t145 = -t271 * t513 + t275 * t516 + t279 * t366;
t146 = -t272 * t513 + t276 * t516 + t280 * t366;
t191 = -t312 * t513 + t316 * t516 + t320 * t366;
t47 = -t172 * t513 + t176 * t516 + t180 * t366 - t250 * t275 + t251 * t279 + t271 * t339;
t48 = -t173 * t513 + t177 * t516 + t181 * t366 - t250 * t276 + t251 * t280 + t272 * t339;
t70 = -t212 * t513 + t214 * t516 + t216 * t366 - t250 * t316 + t251 * t320 + t312 * t339;
t4 = t145 * t337 + t146 * t339 + t191 * t390 - t47 * t514 - t48 * t513 - t512 * t70;
t613 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t139 * t403 + t140 * t405 + t188 * t421 - t41 * t502 - t42 * t503 + t460 * t67;
t6 = t141 * t403 + t142 * t405 + t189 * t421 - t43 * t502 - t44 * t503 + t460 * t68;
t612 = t5 / 0.2e1 + t6 / 0.2e1;
t7 = t143 * t403 + t144 * t405 + t190 * t421 - t45 * t502 - t46 * t503 + t460 * t69;
t8 = t145 * t403 + t146 * t405 + t191 * t421 + t460 * t70 - t47 * t502 - t48 * t503;
t611 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t43 * t462 + t44 * t463 + t482 * t68 + (t141 * t475 - t142 * t477 + t189 * t574) * t493;
t9 = t41 * t462 + t42 * t463 + t482 * t67 + (t139 * t475 - t140 * t477 + t188 * t574) * t493;
t610 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t45 * t462 + t46 * t463 + t482 * t69 + (t143 * t475 - t144 * t477 + t190 * t574) * t493;
t12 = t462 * t47 + t463 * t48 + t482 * t70 + (t145 * t475 - t146 * t477 + t191 * t574) * t493;
t609 = t12 / 0.2e1 + t11 / 0.2e1;
t155 = -t269 * t515 - t273 * t512 + t277 * t409;
t156 = -t270 * t515 - t274 * t512 + t278 * t409;
t201 = -t311 * t515 - t315 * t512 + t319 * t409;
t52 = -t170 * t515 - t174 * t512 + t178 * t409 + t269 * t302 + t273 * t390 + t277 * t303;
t53 = -t171 * t515 - t175 * t512 + t179 * t409 + t270 * t302 + t274 * t390 + t278 * t303;
t99 = -t211 * t515 - t213 * t512 + t215 * t409 + t302 * t311 + t303 * t319 + t315 * t390;
t13 = t155 * t337 + t156 * t339 + t201 * t390 - t512 * t99 - t513 * t53 - t514 * t52;
t100 = -t212 * t512 + t214 * t515 + t216 * t409 - t302 * t316 + t303 * t320 + t312 * t390;
t157 = -t271 * t512 + t275 * t515 + t279 * t409;
t158 = -t272 * t512 + t276 * t515 + t280 * t409;
t202 = -t312 * t512 + t316 * t515 + t320 * t409;
t54 = -t172 * t512 + t176 * t515 + t180 * t409 + t271 * t390 - t275 * t302 + t279 * t303;
t55 = -t173 * t512 + t177 * t515 + t181 * t409 + t272 * t390 - t276 * t302 + t280 * t303;
t14 = -t100 * t512 + t157 * t337 + t158 * t339 + t202 * t390 - t513 * t55 - t514 * t54;
t608 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = t155 * t403 + t156 * t405 + t201 * t421 + t460 * t99 - t502 * t52 - t503 * t53;
t16 = t100 * t460 + t157 * t403 + t158 * t405 + t202 * t421 - t502 * t54 - t503 * t55;
t607 = t16 / 0.2e1 + t15 / 0.2e1;
t17 = t462 * t52 + t463 * t53 + t482 * t99 + (t155 * t475 - t156 * t477 + t201 * t574) * t493;
t18 = t100 * t482 + t462 * t54 + t463 * t55 + (t157 * t475 - t158 * t477 + t202 * t574) * t493;
t606 = t18 / 0.2e1 + t17 / 0.2e1;
t19 = t496 * t67 + (-t41 * t495 + t42 * t492) * t494;
t20 = t496 * t68 + (-t43 * t495 + t44 * t492) * t494;
t605 = t20 / 0.2e1 + t19 / 0.2e1;
t21 = t496 * t69 + (-t45 * t495 + t46 * t492) * t494;
t22 = t496 * t70 + (-t47 * t495 + t48 * t492) * t494;
t604 = t21 / 0.2e1 + t22 / 0.2e1;
t23 = t496 * t99 + (t492 * t53 - t495 * t52) * t494;
t24 = t100 * t496 + (t492 * t55 - t495 * t54) * t494;
t603 = t23 / 0.2e1 + t24 / 0.2e1;
t63 = -t139 * t514 - t140 * t513 - t188 * t512;
t64 = -t141 * t514 - t142 * t513 - t189 * t512;
t602 = t64 / 0.2e1 + t63 / 0.2e1;
t65 = -t143 * t514 - t144 * t513 - t190 * t512;
t66 = -t145 * t514 - t146 * t513 - t191 * t512;
t601 = t66 / 0.2e1 + t65 / 0.2e1;
t71 = -t139 * t502 - t140 * t503 + t188 * t460;
t72 = -t141 * t502 - t142 * t503 + t189 * t460;
t600 = t72 / 0.2e1 + t71 / 0.2e1;
t73 = -t143 * t502 - t144 * t503 + t190 * t460;
t74 = -t145 * t502 - t146 * t503 + t191 * t460;
t599 = t74 / 0.2e1 + t73 / 0.2e1;
t75 = t139 * t462 + t140 * t463 + t188 * t482;
t76 = t141 * t462 + t142 * t463 + t189 * t482;
t598 = t76 / 0.2e1 + t75 / 0.2e1;
t77 = t143 * t462 + t144 * t463 + t190 * t482;
t78 = t145 * t462 + t146 * t463 + t191 * t482;
t597 = t78 / 0.2e1 + t77 / 0.2e1;
t596 = ((-t139 - t141) * t495 + (t140 + t142) * t492) * t676 + (t189 + t188) * t675;
t595 = ((-t143 - t145) * t495 + (t144 + t146) * t492) * t676 + (t190 + t191) * t675;
t87 = -t155 * t514 - t156 * t513 - t201 * t512;
t88 = -t157 * t514 - t158 * t513 - t202 * t512;
t594 = t88 / 0.2e1 + t87 / 0.2e1;
t91 = -t155 * t502 - t156 * t503 + t201 * t460;
t92 = -t157 * t502 - t158 * t503 + t202 * t460;
t593 = t92 / 0.2e1 + t91 / 0.2e1;
t93 = t155 * t462 + t156 * t463 + t201 * t482;
t94 = t157 * t462 + t158 * t463 + t202 * t482;
t592 = t94 / 0.2e1 + t93 / 0.2e1;
t591 = ((-t155 - t157) * t495 + (t156 + t158) * t492) * t676 + (t201 + t202) * t675;
t590 = -t264 - t647;
t589 = -t349 + t646;
t588 = -t265 - t645;
t587 = -t306 - t643;
t586 = -t394 + t642;
t585 = t496 * t265 + t627;
t584 = -t360 - t638;
t583 = -t398 + t637;
t582 = -t361 - t636;
t581 = -t400 - t629;
t580 = -t419 + t628;
t579 = t482 * t350 + t399 * t542 + t419 * t656;
t578 = t496 * t361 + t622;
t569 = 0.2e1 * m(4);
t567 = 0.2e1 * m(5);
t565 = 0.2e1 * m(6);
t563 = 0.2e1 * m(7);
t389 = rSges(4,1) * t422 - rSges(4,2) * t421 + rSges(4,3) * t542;
t469 = (pkin(9) * t493 * t500 + pkin(2) * t501) * t618;
t562 = t494 * (-t389 - t469);
t418 = rSges(4,1) * t461 - rSges(4,2) * t460 + rSges(4,3) * t482;
t464 = t494 * t500 * pkin(2) + pkin(9) * t482;
t561 = (-t418 - t464) * t494;
t235 = Icges(5,5) * t338 - Icges(5,6) * t337 + Icges(5,3) * t403;
t237 = Icges(5,4) * t338 - Icges(5,2) * t337 + Icges(5,6) * t403;
t239 = Icges(5,1) * t338 - Icges(5,4) * t337 + Icges(5,5) * t403;
t313 = Icges(5,5) * t411 + Icges(5,6) * t514 - Icges(5,3) * t502;
t317 = Icges(5,4) * t411 + Icges(5,2) * t514 - Icges(5,6) * t502;
t321 = Icges(5,1) * t411 + Icges(5,4) * t514 - Icges(5,5) * t502;
t104 = -t235 * t503 + t237 * t513 + t239 * t413 + t313 * t405 - t317 * t339 + t321 * t340;
t236 = Icges(5,5) * t340 - Icges(5,6) * t339 + Icges(5,3) * t405;
t238 = Icges(5,4) * t340 - Icges(5,2) * t339 + Icges(5,6) * t405;
t240 = Icges(5,1) * t340 - Icges(5,4) * t339 + Icges(5,5) * t405;
t314 = Icges(5,5) * t413 + Icges(5,6) * t513 - Icges(5,3) * t503;
t318 = Icges(5,4) * t413 + Icges(5,2) * t513 - Icges(5,6) * t503;
t322 = Icges(5,1) * t413 + Icges(5,4) * t513 - Icges(5,5) * t503;
t105 = -t236 * t503 + t238 * t513 + t240 * t413 + t314 * t405 - t318 * t339 + t322 * t340;
t289 = Icges(5,5) * t391 - Icges(5,6) * t390 + Icges(5,3) * t421;
t290 = Icges(5,4) * t391 - Icges(5,2) * t390 + Icges(5,6) * t421;
t291 = Icges(5,1) * t391 - Icges(5,4) * t390 + Icges(5,5) * t421;
t374 = Icges(5,5) * t432 + Icges(5,6) * t512 + Icges(5,3) * t460;
t375 = Icges(5,4) * t432 + Icges(5,2) * t512 + Icges(5,6) * t460;
t376 = Icges(5,1) * t432 + Icges(5,4) * t512 + Icges(5,5) * t460;
t124 = -t289 * t503 + t290 * t513 + t291 * t413 - t339 * t375 + t340 * t376 + t374 * t405;
t199 = -t313 * t503 + t317 * t513 + t321 * t413;
t200 = -t314 * t503 + t318 * t513 + t322 * t413;
t223 = -t374 * t503 + t375 * t513 + t376 * t413;
t26 = -t104 * t502 - t105 * t503 + t124 * t460 + t199 * t403 + t200 * t405 + t223 * t421;
t560 = t26 / 0.2e1 + t611;
t102 = -t235 * t502 + t237 * t514 + t239 * t411 + t313 * t403 - t317 * t337 + t321 * t338;
t103 = -t236 * t502 + t238 * t514 + t240 * t411 + t314 * t403 - t318 * t337 + t322 * t338;
t123 = -t289 * t502 + t290 * t514 + t291 * t411 - t337 * t375 + t338 * t376 + t374 * t403;
t197 = -t313 * t502 + t317 * t514 + t321 * t411;
t198 = -t314 * t502 + t318 * t514 + t322 * t411;
t222 = -t374 * t502 + t375 * t514 + t376 * t411;
t25 = -t102 * t502 - t103 * t503 + t123 * t460 + t197 * t403 + t198 * t405 + t222 * t421;
t559 = -t25 / 0.2e1 - t612;
t558 = -t349 + t590;
t557 = -t394 + t587;
t255 = t463 * t264;
t353 = t361 * t657;
t555 = t255 - t353 + t631;
t554 = -t398 + t584;
t553 = t350 * t672 + t398 * t615 + 0.2e1 * t309 - 0.2e1 * t384;
t552 = -t419 + t581;
t343 = t349 * t655;
t344 = t350 * t651;
t551 = 0.2e1 * t343 + 0.2e1 * t344 + t620;
t550 = t343 + t344 + t619;
t549 = t398 * t655 + t399 * t651 + t621;
t548 = 0.2e1 * t647;
t547 = 0.2e1 * t645;
t546 = 0.2e1 * t638;
t545 = 0.2e1 * t636;
t27 = t102 * t462 + t103 * t463 + t123 * t482 + (t197 * t475 - t198 * t477 + t222 * t574) * t493;
t541 = t27 / 0.2e1 + t610;
t28 = t104 * t462 + t105 * t463 + t124 * t482 + (t199 * t475 - t200 * t477 + t223 * t574) * t493;
t540 = t28 / 0.2e1 + t609;
t110 = t235 * t460 + t237 * t512 + t239 * t432 + t313 * t421 - t317 * t390 + t321 * t391;
t111 = t236 * t460 + t238 * t512 + t240 * t432 + t314 * t421 - t318 * t390 + t322 * t391;
t130 = t289 * t460 + t290 * t512 + t291 * t432 + t374 * t421 - t375 * t390 + t376 * t391;
t204 = t313 * t460 + t317 * t512 + t321 * t432;
t205 = t314 * t460 + t318 * t512 + t322 * t432;
t234 = t374 * t460 + t375 * t512 + t376 * t432;
t29 = -t110 * t502 - t111 * t503 + t130 * t460 + t204 * t403 + t205 * t405 + t234 * t421;
t539 = t29 / 0.2e1 + t607;
t30 = t110 * t462 + t111 * t463 + t130 * t482 + (t204 * t475 - t205 * t477 + t234 * t574) * t493;
t538 = t30 / 0.2e1 + t606;
t31 = t123 * t496 + (-t102 * t495 + t103 * t492) * t494;
t537 = t31 / 0.2e1 + t605;
t32 = t124 * t496 + (-t104 * t495 + t105 * t492) * t494;
t536 = t32 / 0.2e1 + t604;
t35 = t130 * t496 + (-t110 * t495 + t111 * t492) * t494;
t535 = t35 / 0.2e1 + t603;
t534 = t222 * t675 + (-t197 * t495 + t198 * t492) * t676 + t596;
t533 = t223 * t675 + (-t199 * t495 + t200 * t492) * t676 + t595;
t532 = t234 * t675 + (-t204 * t495 + t205 * t492) * t676 + t591;
t531 = t494 * (-t464 + t624);
t530 = (-t469 + t634) * t494;
t242 = rSges(5,1) * t340 - rSges(5,2) * t339 + rSges(5,3) * t405;
t528 = -m(5) * t242 - m(6) * t185;
t527 = m(5) * t323 + m(6) * t282;
t324 = rSges(5,1) * t413 + rSges(5,2) * t513 - rSges(5,3) * t503;
t526 = -m(5) * t324 - m(6) * t284;
t525 = t494 * (-t469 + t586);
t524 = (-t464 + t580) * t494;
t522 = t482 * t265 + t361 * t542 + t400 * t656 + t579;
t261 = t264 * t655;
t262 = t265 * t651;
t520 = t261 + t262 + t550;
t519 = t360 * t655 + t361 * t651 + t549;
t511 = (-t469 + t557) * t494;
t510 = (-t464 + t552) * t494;
t96 = -t183 * t513 + t185 * t514 + t282 * t339 - t284 * t337;
t506 = m(5) * t241 + m(6) * t183 + t548 * t669;
t505 = -t547 * t669 + t528;
t334 = rSges(4,1) * t404 - rSges(4,2) * t403 + rSges(4,3) * t657;
t504 = m(4) * t334 + t506;
t481 = (rSges(3,1) * t501 - rSges(3,2) * t500) * t618;
t480 = (Icges(3,1) * t501 - t660) * t618;
t479 = (-Icges(3,2) * t500 + t659) * t618;
t478 = (Icges(3,5) * t501 - Icges(3,6) * t500) * t618;
t472 = t496 * rSges(3,3) + (rSges(3,1) * t500 + rSges(3,2) * t501) * t494;
t471 = Icges(3,5) * t496 + (Icges(3,1) * t500 + t659) * t494;
t470 = Icges(3,6) * t496 + (Icges(3,2) * t501 + t660) * t494;
t456 = rSges(3,1) * t476 + rSges(3,2) * t477;
t455 = rSges(3,1) * t474 - rSges(3,2) * t475;
t454 = Icges(3,1) * t476 + Icges(3,4) * t477;
t453 = Icges(3,1) * t474 - Icges(3,4) * t475;
t452 = Icges(3,4) * t476 + Icges(3,2) * t477;
t451 = Icges(3,4) * t474 - Icges(3,2) * t475;
t450 = Icges(3,5) * t476 + Icges(3,6) * t477;
t449 = Icges(3,5) * t474 - Icges(3,6) * t475;
t448 = -rSges(3,1) * t518 + rSges(3,2) * t485 + rSges(3,3) * t655;
t447 = rSges(3,1) * t484 + rSges(3,2) * t483 - rSges(3,3) * t651;
t444 = -Icges(3,1) * t518 + Icges(3,4) * t485 + Icges(3,5) * t655;
t443 = Icges(3,1) * t484 + Icges(3,4) * t483 - Icges(3,5) * t651;
t442 = -Icges(3,4) * t518 + Icges(3,2) * t485 + Icges(3,6) * t655;
t441 = Icges(3,4) * t484 + Icges(3,2) * t483 - Icges(3,6) * t651;
t417 = Icges(4,1) * t461 - Icges(4,4) * t460 + Icges(4,5) * t482;
t416 = Icges(4,4) * t461 - Icges(4,2) * t460 + Icges(4,6) * t482;
t415 = Icges(4,5) * t461 - Icges(4,6) * t460 + Icges(4,3) * t482;
t388 = Icges(4,1) * t422 - Icges(4,4) * t421 + Icges(4,5) * t542;
t387 = Icges(4,4) * t422 - Icges(4,2) * t421 + Icges(4,6) * t542;
t386 = Icges(4,5) * t422 - Icges(4,6) * t421 + Icges(4,3) * t542;
t378 = rSges(4,1) * t430 + rSges(4,2) * t503 + rSges(4,3) * t463;
t377 = rSges(4,1) * t428 + rSges(4,2) * t502 + rSges(4,3) * t462;
t373 = Icges(4,1) * t430 + Icges(4,4) * t503 + Icges(4,5) * t463;
t372 = Icges(4,1) * t428 + Icges(4,4) * t502 + Icges(4,5) * t462;
t371 = Icges(4,4) * t430 + Icges(4,2) * t503 + Icges(4,6) * t463;
t370 = Icges(4,4) * t428 + Icges(4,2) * t502 + Icges(4,6) * t462;
t369 = Icges(4,5) * t430 + Icges(4,6) * t503 + Icges(4,3) * t463;
t368 = Icges(4,5) * t428 + Icges(4,6) * t502 + Icges(4,3) * t462;
t362 = t502 * t400;
t347 = t460 * t361;
t333 = Icges(4,1) * t406 - Icges(4,4) * t405 - Icges(4,5) * t656;
t332 = Icges(4,1) * t404 - Icges(4,4) * t403 + Icges(4,5) * t657;
t331 = Icges(4,4) * t406 - Icges(4,2) * t405 - Icges(4,6) * t656;
t330 = Icges(4,4) * t404 - Icges(4,2) * t403 + Icges(4,6) * t657;
t329 = Icges(4,5) * t406 - Icges(4,6) * t405 - Icges(4,3) * t656;
t328 = Icges(4,5) * t404 - Icges(4,6) * t403 + Icges(4,3) * t657;
t327 = t503 * t360;
t305 = t378 * t482 - t418 * t463;
t304 = -t377 * t482 + t418 * t462;
t294 = (-t377 - t438) * t496 + t495 * t561;
t293 = t378 * t496 + t492 * t561 + t433;
t288 = t415 * t482 - t416 * t460 + t417 * t461;
t287 = t377 * t463 - t378 * t462;
t286 = t415 * t463 + t416 * t503 + t417 * t430;
t285 = t415 * t462 + t416 * t502 + t417 * t428;
t268 = (t377 * t492 + t378 * t495) * t494 + t621;
t267 = (-t334 - t445) * t496 + t495 * t562;
t266 = t335 * t496 + t492 * t562 + t440;
t258 = t324 * t460 + t379 * t503;
t257 = -t323 * t460 - t379 * t502;
t247 = t369 * t482 - t371 * t460 + t373 * t461;
t246 = t368 * t482 - t370 * t460 + t372 * t461;
t233 = t369 * t463 + t371 * t503 + t373 * t430;
t232 = t368 * t463 + t370 * t503 + t372 * t430;
t231 = t369 * t462 + t371 * t502 + t373 * t428;
t230 = t368 * t462 + t370 * t502 + t372 * t428;
t229 = (t334 * t492 + t335 * t495) * t494 + t619;
t228 = -t323 * t503 + t324 * t502;
t227 = t324 * t482 + t463 * t624 + t392;
t226 = t379 * t462 + t482 * t630 + t407;
t225 = (-t438 + t630) * t496 + t495 * t531;
t224 = t324 * t496 + t492 * t531 + t622;
t220 = t335 * t482 - t389 * t463 + (t378 * t574 + t418 * t477) * t493;
t219 = -t334 * t482 + t389 * t462 + (-t377 * t574 + t418 * t475) * t493;
t210 = -t284 * t512 + t326 * t513;
t209 = t282 * t512 - t326 * t514;
t208 = t323 * t463 + t380 + (-t324 - t399) * t462;
t207 = (t323 * t492 + t324 * t495) * t494 + t549;
t206 = t386 * t482 - t387 * t460 + t388 * t461 + t415 * t542 - t416 * t421 + t417 * t422;
t203 = t334 * t463 - t335 * t462 + (-t377 * t477 - t378 * t475) * t493;
t196 = -t282 * t513 + t284 * t514;
t195 = t386 * t463 + t387 * t503 + t388 * t430 - t405 * t416 + t406 * t417 - t415 * t656;
t194 = t386 * t462 + t387 * t502 + t388 * t428 - t403 * t416 + t404 * t417 + t415 * t657;
t193 = t284 * t460 - t503 * t628 + t347;
t192 = -t326 * t502 + t460 * t637 - t362;
t169 = (-t438 + t583) * t496 + t495 * t524;
t168 = t284 * t496 + t492 * t524 + t578;
t167 = t284 * t482 + t463 * t580 + t625;
t166 = t326 * t462 + t482 * t583 + t623;
t165 = (-t445 + t641) * t496 + t495 * t530;
t164 = t242 * t496 + t492 * t530 + t627;
t163 = t329 * t482 - t331 * t460 + t333 * t461 + t369 * t542 - t371 * t421 + t373 * t422;
t162 = t328 * t482 - t330 * t460 + t332 * t461 + t368 * t542 - t370 * t421 + t372 * t422;
t161 = -t282 * t503 - t502 * t635 - t327;
t160 = -t512 * t636 + t513 * t629;
t159 = t512 * t638 - t514 * t629;
t154 = t329 * t463 + t331 * t503 + t333 * t430 - t369 * t656 - t371 * t405 + t373 * t406;
t153 = t328 * t463 + t330 * t503 + t332 * t430 - t368 * t656 - t370 * t405 + t372 * t406;
t152 = t329 * t462 + t331 * t502 + t333 * t428 + t369 * t657 - t371 * t403 + t373 * t404;
t151 = t328 * t462 + t330 * t502 + t332 * t428 + t368 * t657 - t370 * t403 + t372 * t404;
t150 = t242 * t460 + t292 * t503 + t324 * t421 - t379 * t405;
t149 = -t241 * t460 - t292 * t502 - t323 * t421 + t379 * t403;
t148 = (t282 * t492 + t284 * t495) * t494 + t519;
t147 = t282 * t463 + (-t399 + t635) * t462 + t626;
t138 = (t241 * t492 + t242 * t495) * t494 + t550;
t137 = t460 * t636 - t503 * t581 + t347;
t136 = t460 * t584 - t502 * t629 - t362;
t135 = (-t438 + t554) * t496 + t495 * t510;
t134 = t492 * t510 + t496 * t636 + t578;
t133 = -t513 * t638 + t514 * t636;
t132 = t463 * t552 + t482 * t636 + t625;
t131 = t462 * t629 + t482 * t554 + t623;
t129 = -t241 * t503 + t242 * t502 + t323 * t405 - t324 * t403;
t128 = t242 * t482 + (t324 * t574 + t379 * t477) * t493 + t634 * t463 + t579;
t127 = t292 * t462 + t367 + t641 * t482 + (-t475 * t624 + t574 * t630) * t493;
t126 = -t502 * t582 - t503 * t638 - t327;
t125 = (t492 * t638 + t495 * t636) * t494 + t519;
t122 = t638 * t463 + (-t399 + t582) * t462 + t626;
t120 = t204 * t462 + t205 * t463 + t234 * t482;
t119 = -t204 * t502 - t205 * t503 + t234 * t460;
t118 = t241 * t463 + (-t242 - t350) * t462 + (-t324 * t475 + t477 * t630) * t493 + t631;
t115 = t199 * t462 + t200 * t463 + t223 * t482;
t114 = t197 * t462 + t198 * t463 + t222 * t482;
t113 = -t199 * t502 - t200 * t503 + t223 * t460;
t112 = -t197 * t502 - t198 * t503 + t222 * t460;
t109 = -t185 * t512 + t218 * t513 + t284 * t390 - t326 * t339;
t108 = t183 * t512 - t218 * t514 - t282 * t390 + t326 * t337;
t107 = (-t445 + t589) * t496 + t495 * t525;
t106 = t185 * t496 + t492 * t525 + t585;
t101 = t206 * t496 + (-t162 * t495 + t163 * t492) * t494;
t95 = (t183 * t492 + t185 * t495) * t494 + t520;
t90 = t195 * t496 + (-t153 * t495 + t154 * t492) * t494;
t89 = t194 * t496 + (-t151 * t495 + t152 * t492) * t494;
t86 = t185 * t460 + t284 * t421 + t405 * t628 - t503 * t642 + t639;
t85 = -t218 * t502 + t326 * t403 + t421 * t637 + t460 * t646 + t633;
t84 = t185 * t482 + (t284 * t574 + t326 * t477) * t493 + t586 * t463 + t522;
t83 = t218 * t462 + t589 * t482 + (-t475 * t580 + t574 * t583) * t493 + t632;
t62 = (-t445 + t558) * t496 + t495 * t511;
t61 = t492 * t511 + t496 * t645 + t585;
t60 = t162 * t462 + t163 * t463 + t206 * t482 + (t246 * t475 - t247 * t477 + t288 * t574) * t493;
t59 = -t183 * t503 + t282 * t405 + t403 * t635 - t502 * t644 + t640;
t58 = t183 * t463 + (-t350 + t644) * t462 + (-t284 * t475 + t477 * t583) * t493 + t555;
t57 = t153 * t462 + t154 * t463 + t195 * t482 + (t232 * t475 - t233 * t477 + t286 * t574) * t493;
t56 = t151 * t462 + t152 * t463 + t194 * t482 + (t230 * t475 - t231 * t477 + t285 * t574) * t493;
t51 = -t339 * t629 + t390 * t636 - t512 * t645 + t513 * t643;
t50 = t337 * t629 - t390 * t638 + t512 * t647 - t514 * t643;
t49 = (t492 * t647 + t495 * t645) * t494 + t520;
t40 = t405 * t581 + t421 * t636 + t460 * t645 - t503 * t587 + t639;
t39 = t403 * t629 + t421 * t584 + t460 * t590 - t502 * t643 + t633;
t38 = t645 * t482 + (t477 * t629 + t574 * t636) * t493 + t557 * t463 + t522;
t37 = t643 * t462 + t558 * t482 + (-t475 * t552 + t554 * t574) * t493 + t632;
t36 = -t337 * t636 + t339 * t638 - t513 * t647 + t514 * t645;
t34 = t403 * t582 + t405 * t638 - t502 * t588 - t503 * t647 + t640;
t33 = t647 * t463 + (-t350 + t588) * t462 + (-t475 * t636 + t477 * t554) * t493 + t555;
t79 = [0; m(4) * t620 / 0.2e1 + t551 * t671 + (m(3) * t456 + m(7) * t645 - t528 + t662) * t651 + (m(3) * t455 + t504) * t655 + t673 * (0.2e1 * t261 + 0.2e1 * t262 + t551); t21 * t655 - t19 * t651 + t22 * t655 - t20 * t651 + t32 * t655 - t31 * t651 + t90 * t655 - t89 * t651 + ((t442 * t477 + t444 * t476 + t450 * t655 + t452 * t485 - t454 * t518) * t655 - (t441 * t477 + t443 * t476 + t449 * t655 + t451 * t485 - t453 * t518) * t651 + (t470 * t477 + t471 * t476 + t478 * t655 + t479 * t485 - t480 * t518) * t496) * t655 - ((-t442 * t475 + t444 * t474 - t450 * t651 + t452 * t483 + t454 * t484) * t655 - (-t441 * t475 + t443 * t474 - t449 * t651 + t451 * t483 + t453 * t484) * t651 + (-t470 * t475 + t471 * t474 - t478 * t651 + t479 * t483 + t480 * t484) * t496) * t651 + (t125 * t49 + t134 * t61 + t135 * t62) * t563 + (t106 * t168 + t107 * t169 + t148 * t95) * t565 + t496 * t23 + t496 * t24 + (t138 * t207 + t164 * t224 + t165 * t225) * t567 + t496 * t35 + (t229 * t268 + t266 * t293 + t267 * t294) * t569 + t496 * t101 + t496 * (t496 ^ 2 * t478 + (((t452 * t501 + t454 * t500) * t492 - (t451 * t501 + t453 * t500) * t495 + ((-t442 * t500 + t444 * t501) * t492 - (-t441 * t500 + t443 * t501) * t495) * qJD(2)) * t494 + (-t449 * t495 + t450 * t492 + t479 * t501 + t480 * t500 + (-t470 * t500 + t471 * t501) * qJD(2)) * t496) * t494) + 0.2e1 * m(3) * ((-t447 * t496 - t472 * t651) * (-t455 * t496 - t481 * t651) + (t448 * t496 - t472 * t655) * (t456 * t496 - t481 * t655) + (t447 * t492 + t448 * t495) * t494 ^ 2 * (t455 * t492 + t456 * t495)); t553 * t671 + t504 * t463 + ((-m(4) * t377 - m(7) * t638 - t527) * t477 - (m(4) * t378 + m(7) * t636 - t526) * t475) * t493 + (t505 - t662) * t462 + t673 * (t265 * t672 + t360 * t615 + 0.2e1 * t255 - 0.2e1 * t353 + t553); (t122 * t49 + t125 * t33 + t131 * t62 + t132 * t61 + t134 * t38 + t135 * t37) * m(7) + m(4) * (t203 * t268 + t219 * t294 + t220 * t293 + t229 * t287 + t266 * t305 + t267 * t304) + (t118 * t207 + t127 * t225 + t128 * t224 + t138 * t208 + t164 * t227 + t165 * t226) * m(5) + (t106 * t167 + t107 * t166 + t147 * t95 + t148 * t58 + t168 * t84 + t169 * t83) * m(6) + (t60 / 0.2e1 + t538) * t496 + (t101 / 0.2e1 + t535) * t482 + (t90 / 0.2e1 + t536) * t463 + (t89 / 0.2e1 + t537) * t462 + ((-t56 / 0.2e1 - t541) * t495 + (t57 / 0.2e1 + t540) * t492) * t494 + ((t286 * t666 + (t232 * t667 + t233 * t668) * t494 - t533) * t477 - (t285 * t666 + (t230 * t667 + t231 * t668) * t494 - t534) * t475 + (t288 * t675 + (-t246 * t495 / 0.2e1 + t247 * t492 / 0.2e1) * t494 + t532) * t574) * t493; (t122 * t33 + t131 * t37 + t132 * t38) * t563 + (t203 * t287 + t219 * t304 + t220 * t305) * t569 + (t118 * t208 + t127 * t226 + t128 * t227) * t567 + (t147 * t58 + t166 * t83 + t167 * t84) * t565 + (t60 + t30 + t18 + t17) * t482 + (t57 + t28 + t12 + t11) * t463 + (t56 + t27 + t10 + t9) * t462 + ((-t232 * t462 - t233 * t463 - t286 * t482 - t115 - t77 - t78) * t477 - (-t230 * t462 - t231 * t463 - t285 * t482 - t114 - t75 - t76) * t475 + (t246 * t462 + t247 * t463 + t288 * t482 + t120 + t93 + t94) * t574) * t493; -t506 * t503 - t505 * t502 + (t546 * t669 + t527) * t405 + (-t545 * t669 + t526) * t403 + 0.2e1 * t673 * (t265 * t502 - t403 * t361 + t640); (t125 * t34 + t126 * t49 + t134 * t40 + t135 * t39 + t136 * t62 + t137 * t61) * m(7) + (t129 * t207 + t138 * t228 + t149 * t225 + t150 * t224 + t164 * t258 + t165 * t257) * m(5) + (t106 * t193 + t107 * t192 + t148 * t59 + t161 * t95 + t168 * t86 + t169 * t85) * m(6) + t539 * t496 + t535 * t460 - t536 * t503 - t537 * t502 + t532 * t421 + t533 * t405 + t534 * t403 + (t492 * t560 + t495 * t559) * t494; (t122 * t34 + t126 * t33 + t131 * t39 + t132 * t40 + t136 * t37 + t137 * t38) * m(7) + (t118 * t228 + t127 * t257 + t128 * t258 + t129 * t208 + t149 * t226 + t150 * t227) * m(5) + (t147 * t59 + t161 * t58 + t166 * t85 + t167 * t86 + t192 * t83 + t193 * t84) * m(6) + t539 * t482 + t560 * t463 - t559 * t462 + t538 * t460 - t540 * t503 - t541 * t502 + (t120 / 0.2e1 + t592) * t421 + (t115 / 0.2e1 + t597) * t405 + (t114 / 0.2e1 + t598) * t403 + ((-t113 / 0.2e1 - t599) * t477 - (-t112 / 0.2e1 - t600) * t475 + (t119 / 0.2e1 + t593) * t574) * t493; (t126 * t34 + t136 * t39 + t137 * t40) * t563 + (t161 * t59 + t192 * t85 + t193 * t86) * t565 + (t129 * t228 + t149 * t257 + t150 * t258) * t567 + (t16 + t15 + t29) * t460 - (t7 + t8 + t26) * t503 - (t6 + t5 + t25) * t502 + (t92 + t91 + t119) * t421 + (t73 + t74 + t113) * t405 + (t72 + t71 + t112) * t403; t96 * m(6) + (-t337 * t545 + t339 * t546 - t513 * t548 + t514 * t547) * t669; (t125 * t36 + t133 * t49 + t134 * t51 + t135 * t50 + t159 * t62 + t160 * t61) * m(7) + (t106 * t210 + t107 * t209 + t108 * t169 + t109 * t168 + t148 * t96 + t196 * t95) * m(6) + t608 * t496 - t603 * t512 - t604 * t513 - t605 * t514 + t591 * t390 + t595 * t339 + t596 * t337 + (t492 * t613 - t495 * t614) * t494; (t122 * t36 + t131 * t50 + t132 * t51 + t133 * t33 + t159 * t37 + t160 * t38) * m(7) + (t108 * t166 + t109 * t167 + t147 * t96 + t196 * t58 + t209 * t83 + t210 * t84) * m(6) + t608 * t482 + t613 * t463 + t614 * t462 - t606 * t512 - t609 * t513 - t610 * t514 + t592 * t390 + t597 * t339 + t598 * t337 + (t475 * t602 - t477 * t601 + t574 * t594) * t493; (t126 * t36 + t133 * t34 + t136 * t50 + t137 * t51 + t159 * t39 + t160 * t40) * m(7) + (t108 * t192 + t109 * t193 + t161 * t96 + t196 * t59 + t209 * t85 + t210 * t86) * m(6) + t608 * t460 - t607 * t512 - t613 * t503 - t614 * t502 + t594 * t421 - t611 * t513 - t612 * t514 + t601 * t405 + t602 * t403 + t593 * t390 + t599 * t339 + t600 * t337; (t133 * t36 + t159 * t50 + t160 * t51) * t563 + (t108 * t209 + t109 * t210 + t196 * t96) * t565 - (t14 + t13) * t512 - (t4 + t3) * t513 - (t1 + t2) * t514 + (t88 + t87) * t390 + (t66 + t65) * t339 + (t63 + t64) * t337; 0.2e1 * t302 * t669; (t125 * t302 + t134 * t248 + t135 * t250 - t49 * t515 - t516 * t62 - t517 * t61) * m(7); (t122 * t302 + t131 * t250 + t132 * t248 - t33 * t515 - t37 * t516 - t38 * t517) * m(7); (t126 * t302 + t136 * t250 + t137 * t248 - t34 * t515 - t39 * t516 - t40 * t517) * m(7); (t133 * t302 + t159 * t250 + t160 * t248 - t36 * t515 - t50 * t516 - t51 * t517) * m(7); (-t248 * t517 - t250 * t516 - t302 * t515) * t563;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t79(1) t79(2) t79(4) t79(7) t79(11) t79(16); t79(2) t79(3) t79(5) t79(8) t79(12) t79(17); t79(4) t79(5) t79(6) t79(9) t79(13) t79(18); t79(7) t79(8) t79(9) t79(10) t79(14) t79(19); t79(11) t79(12) t79(13) t79(14) t79(15) t79(20); t79(16) t79(17) t79(18) t79(19) t79(20) t79(21);];
Mq  = res;
