% Calculate vector of inverse dynamics joint torques for
% S5PPRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:38
% EndTime: 2019-12-05 15:19:42
% DurationCPUTime: 40.24s
% Computational Cost: add. (77274->1278), mult. (223318->1783), div. (0->0), fcn. (281150->14), ass. (0->425)
t572 = 2 * qJD(3);
t571 = 2 * qJDD(3);
t408 = sin(pkin(10));
t410 = cos(pkin(10));
t519 = sin(pkin(11));
t523 = cos(pkin(5));
t466 = t523 * t519;
t521 = cos(pkin(11));
t396 = -t408 * t466 + t410 * t521;
t413 = sin(qJ(3));
t468 = t523 * t521;
t432 = t408 * t468 + t410 * t519;
t522 = cos(pkin(6));
t422 = t432 * t522;
t409 = sin(pkin(5));
t520 = sin(pkin(6));
t478 = t409 * t520;
t533 = cos(qJ(3));
t362 = t396 * t533 + (t408 * t478 - t422) * t413;
t479 = t409 * t522;
t387 = t408 * t479 + t432 * t520;
t412 = sin(qJ(4));
t532 = cos(qJ(4));
t304 = t362 * t532 + t387 * t412;
t454 = t533 * t478;
t361 = t396 * t413 - t408 * t454 + t422 * t533;
t339 = t361 * qJD(3);
t201 = qJD(4) * t304 - t339 * t412;
t340 = t362 * qJD(3);
t380 = qJDD(3) * t387;
t214 = qJD(4) * t340 + qJDD(4) * t361 + t380;
t451 = -t362 * t412 + t387 * t532;
t100 = qJD(5) * t201 - qJDD(5) * t451 + t214;
t465 = t522 * t521;
t467 = t523 * t520;
t385 = t413 * t467 + (t413 * t465 + t519 * t533) * t409;
t394 = -t478 * t521 + t522 * t523;
t364 = t385 * t532 + t394 * t412;
t566 = t409 * (-t519 * t413 + t533 * t465) + t533 * t467;
t369 = t566 * qJD(3);
t292 = qJD(4) * t364 + t369 * t412;
t370 = t385 * qJD(3);
t392 = qJDD(3) * t394;
t309 = qJD(4) * t370 - qJDD(4) * t566 + t392;
t450 = -t385 * t412 + t394 * t532;
t157 = qJD(5) * t292 - qJDD(5) * t450 + t309;
t395 = t408 * t521 + t410 * t466;
t433 = t408 * t519 - t410 * t468;
t423 = t433 * t522;
t359 = t395 * t413 + t410 * t454 + t423 * t533;
t386 = -t410 * t479 + t433 * t520;
t382 = qJD(3) * t386;
t307 = qJD(4) * t359 + t382;
t360 = t395 * t533 + (-t410 * t478 - t423) * t413;
t452 = -t360 * t412 + t386 * t532;
t181 = -qJD(5) * t452 + t307;
t383 = qJD(3) * t387;
t308 = qJD(4) * t361 + t383;
t182 = -qJD(5) * t451 + t308;
t337 = t359 * qJD(3);
t200 = qJD(4) * t452 - t337 * t532;
t302 = t360 * t532 + t386 * t412;
t411 = sin(qJ(5));
t414 = cos(qJ(5));
t208 = t302 * t414 + t359 * t411;
t338 = t360 * qJD(3);
t103 = -qJD(5) * t208 - t200 * t411 + t338 * t414;
t207 = -t302 * t411 + t359 * t414;
t104 = qJD(5) * t207 + t200 * t414 + t338 * t411;
t305 = -t364 * t411 - t414 * t566;
t306 = t364 * t414 - t411 * t566;
t153 = Icges(6,5) * t306 + Icges(6,6) * t305 - Icges(6,3) * t450;
t510 = Icges(6,4) * t306;
t154 = Icges(6,2) * t305 - Icges(6,6) * t450 + t510;
t298 = Icges(6,4) * t305;
t155 = Icges(6,1) * t306 - Icges(6,5) * t450 + t298;
t199 = qJD(4) * t302 - t337 * t412;
t293 = qJD(4) * t450 + t369 * t532;
t158 = -qJD(5) * t306 - t293 * t411 + t370 * t414;
t159 = qJD(5) * t305 + t293 * t414 + t370 * t411;
t82 = Icges(6,5) * t159 + Icges(6,6) * t158 + Icges(6,3) * t292;
t83 = Icges(6,4) * t159 + Icges(6,2) * t158 + Icges(6,6) * t292;
t84 = Icges(6,1) * t159 + Icges(6,4) * t158 + Icges(6,5) * t292;
t19 = t103 * t154 + t104 * t155 + t153 * t199 + t207 * t83 + t208 * t84 - t452 * t82;
t393 = qJD(3) * t394;
t365 = -qJD(4) * t566 + t393;
t279 = -qJD(5) * t450 + t365;
t91 = Icges(6,5) * t208 + Icges(6,6) * t207 - Icges(6,3) * t452;
t512 = Icges(6,4) * t208;
t93 = Icges(6,2) * t207 - Icges(6,6) * t452 + t512;
t203 = Icges(6,4) * t207;
t95 = Icges(6,1) * t208 - Icges(6,5) * t452 + t203;
t39 = t207 * t93 + t208 * t95 - t452 * t91;
t209 = -t304 * t411 + t361 * t414;
t210 = t304 * t414 + t361 * t411;
t92 = Icges(6,5) * t210 + Icges(6,6) * t209 - Icges(6,3) * t451;
t511 = Icges(6,4) * t210;
t94 = Icges(6,2) * t209 - Icges(6,6) * t451 + t511;
t204 = Icges(6,4) * t209;
t96 = Icges(6,1) * t210 - Icges(6,5) * t451 + t204;
t40 = t207 * t94 + t208 * t96 - t452 * t92;
t46 = -t153 * t452 + t154 * t207 + t155 * t208;
t54 = Icges(6,5) * t104 + Icges(6,6) * t103 + Icges(6,3) * t199;
t56 = Icges(6,4) * t104 + Icges(6,2) * t103 + Icges(6,6) * t199;
t58 = Icges(6,1) * t104 + Icges(6,4) * t103 + Icges(6,5) * t199;
t8 = t103 * t93 + t104 * t95 + t199 * t91 + t207 * t56 + t208 * t58 - t452 * t54;
t202 = qJD(4) * t451 - t339 * t532;
t105 = -qJD(5) * t210 - t202 * t411 + t340 * t414;
t106 = qJD(5) * t209 + t202 * t414 + t340 * t411;
t55 = Icges(6,5) * t106 + Icges(6,6) * t105 + Icges(6,3) * t201;
t57 = Icges(6,4) * t106 + Icges(6,2) * t105 + Icges(6,6) * t201;
t59 = Icges(6,1) * t106 + Icges(6,4) * t105 + Icges(6,5) * t201;
t9 = t103 * t94 + t104 * t96 + t199 * t92 + t207 * t57 + t208 * t59 - t452 * t55;
t379 = qJDD(3) * t386;
t213 = qJD(4) * t338 + qJDD(4) * t359 + t379;
t99 = qJD(5) * t199 - qJDD(5) * t452 + t213;
t1 = t100 * t40 + t157 * t46 + t181 * t8 + t182 * t9 + t19 * t279 + t39 * t99;
t107 = Icges(5,5) * t200 - Icges(5,6) * t199 + Icges(5,3) * t338;
t109 = Icges(5,4) * t200 - Icges(5,2) * t199 + Icges(5,6) * t338;
t111 = Icges(5,1) * t200 - Icges(5,4) * t199 + Icges(5,5) * t338;
t145 = Icges(5,5) * t302 + Icges(5,6) * t452 + Icges(5,3) * t359;
t515 = Icges(5,4) * t302;
t147 = Icges(5,2) * t452 + Icges(5,6) * t359 + t515;
t294 = Icges(5,4) * t452;
t149 = Icges(5,1) * t302 + Icges(5,5) * t359 + t294;
t25 = t107 * t359 + t109 * t452 + t111 * t302 + t145 * t338 - t147 * t199 + t149 * t200;
t108 = Icges(5,5) * t202 - Icges(5,6) * t201 + Icges(5,3) * t340;
t110 = Icges(5,4) * t202 - Icges(5,2) * t201 + Icges(5,6) * t340;
t112 = Icges(5,1) * t202 - Icges(5,4) * t201 + Icges(5,5) * t340;
t146 = Icges(5,5) * t304 + Icges(5,6) * t451 + Icges(5,3) * t361;
t514 = Icges(5,4) * t304;
t148 = Icges(5,2) * t451 + Icges(5,6) * t361 + t514;
t295 = Icges(5,4) * t451;
t150 = Icges(5,1) * t304 + Icges(5,5) * t361 + t295;
t26 = t108 * t359 + t110 * t452 + t112 * t302 + t146 * t338 - t148 * t199 + t150 * t200;
t160 = Icges(5,5) * t293 - Icges(5,6) * t292 + Icges(5,3) * t370;
t161 = Icges(5,4) * t293 - Icges(5,2) * t292 + Icges(5,6) * t370;
t162 = Icges(5,1) * t293 - Icges(5,4) * t292 + Icges(5,5) * t370;
t228 = Icges(5,5) * t364 + Icges(5,6) * t450 - Icges(5,3) * t566;
t513 = Icges(5,4) * t364;
t229 = Icges(5,2) * t450 - Icges(5,6) * t566 + t513;
t356 = Icges(5,4) * t450;
t230 = Icges(5,1) * t364 - Icges(5,5) * t566 + t356;
t35 = t160 * t359 + t161 * t452 + t162 * t302 - t199 * t229 + t200 * t230 + t228 * t338;
t48 = t145 * t359 + t147 * t452 + t149 * t302;
t49 = t146 * t359 + t148 * t452 + t150 * t302;
t72 = t228 * t359 + t229 * t452 + t230 * t302;
t570 = t213 * t48 + t214 * t49 + t25 * t307 + t26 * t308 + t309 * t72 + t35 * t365 + t1;
t10 = t105 * t93 + t106 * t95 + t201 * t91 + t209 * t56 + t210 * t58 - t451 * t54;
t11 = t105 * t94 + t106 * t96 + t201 * t92 + t209 * t57 + t210 * t59 - t451 * t55;
t20 = t105 * t154 + t106 * t155 + t153 * t201 + t209 * t83 + t210 * t84 - t451 * t82;
t41 = t209 * t93 + t210 * t95 - t451 * t91;
t42 = t209 * t94 + t210 * t96 - t451 * t92;
t47 = -t153 * t451 + t154 * t209 + t155 * t210;
t2 = t10 * t181 + t100 * t42 + t182 * t11 + t157 * t47 + t20 * t279 + t41 * t99;
t27 = t107 * t361 + t109 * t451 + t111 * t304 + t145 * t340 - t147 * t201 + t149 * t202;
t28 = t108 * t361 + t110 * t451 + t112 * t304 + t146 * t340 - t148 * t201 + t150 * t202;
t36 = t160 * t361 + t161 * t451 + t162 * t304 - t201 * t229 + t202 * t230 + t228 * t340;
t50 = t145 * t361 + t147 * t451 + t149 * t304;
t51 = t146 * t361 + t148 * t451 + t150 * t304;
t73 = t228 * t361 + t229 * t451 + t230 * t304;
t569 = t213 * t50 + t214 * t51 + t27 * t307 + t308 * t28 + t309 * t73 + t36 * t365 + t2;
t14 = t158 * t93 + t159 * t95 + t292 * t91 + t305 * t56 + t306 * t58 - t450 * t54;
t15 = t158 * t94 + t159 * t96 + t292 * t92 + t305 * t57 + t306 * t59 - t450 * t55;
t21 = t153 * t292 + t154 * t158 + t155 * t159 + t305 * t83 + t306 * t84 - t450 * t82;
t44 = t305 * t93 + t306 * t95 - t450 * t91;
t45 = t305 * t94 + t306 * t96 - t450 * t92;
t52 = -t153 * t450 + t154 * t305 + t155 * t306;
t3 = t100 * t45 + t14 * t181 + t15 * t182 + t157 * t52 + t21 * t279 + t44 * t99;
t30 = -t107 * t566 + t109 * t450 + t111 * t364 + t145 * t370 - t147 * t292 + t149 * t293;
t31 = -t108 * t566 + t110 * t450 + t112 * t364 + t146 * t370 - t148 * t292 + t150 * t293;
t43 = -t160 * t566 + t161 * t450 + t162 * t364 + t228 * t370 - t229 * t292 + t230 * t293;
t66 = -t145 * t566 + t147 * t450 + t149 * t364;
t67 = -t146 * t566 + t148 * t450 + t150 * t364;
t81 = -t228 * t566 + t229 * t450 + t230 * t364;
t568 = t213 * t66 + t214 * t67 + t30 * t307 + t308 * t31 + t309 * t81 + t365 * t43 + t3;
t567 = rSges(6,1) * t414 - rSges(6,2) * t411;
t565 = (-m(3) - m(5));
t563 = m(3) / 0.2e1;
t562 = m(5) / 0.2e1;
t561 = m(6) / 0.2e1;
t560 = t99 / 0.2e1;
t559 = t100 / 0.2e1;
t558 = t157 / 0.2e1;
t557 = -t181 / 0.2e1;
t556 = t181 / 0.2e1;
t555 = -t182 / 0.2e1;
t554 = t182 / 0.2e1;
t553 = t213 / 0.2e1;
t552 = t214 / 0.2e1;
t551 = -t279 / 0.2e1;
t550 = t279 / 0.2e1;
t549 = -t307 / 0.2e1;
t548 = t307 / 0.2e1;
t547 = -t308 / 0.2e1;
t546 = t308 / 0.2e1;
t545 = t309 / 0.2e1;
t540 = -t365 / 0.2e1;
t539 = t365 / 0.2e1;
t117 = pkin(4) * t200 + pkin(9) * t199;
t60 = rSges(6,1) * t104 + rSges(6,2) * t103 + rSges(6,3) * t199;
t528 = t117 + t60;
t118 = pkin(4) * t202 + pkin(9) * t201;
t61 = rSges(6,1) * t106 + rSges(6,2) * t105 + rSges(6,3) * t201;
t527 = t118 + t61;
t164 = pkin(4) * t293 + pkin(9) * t292;
t85 = rSges(6,1) * t159 + rSges(6,2) * t158 + rSges(6,3) * t292;
t526 = t164 + t85;
t177 = pkin(4) * t302 - pkin(9) * t452;
t97 = rSges(6,1) * t208 + rSges(6,2) * t207 - rSges(6,3) * t452;
t525 = t177 + t97;
t179 = pkin(4) * t304 - pkin(9) * t451;
t98 = rSges(6,1) * t210 + rSges(6,2) * t209 - rSges(6,3) * t451;
t524 = t179 + t98;
t518 = Icges(4,4) * t360;
t517 = Icges(4,4) * t362;
t516 = Icges(4,4) * t385;
t113 = rSges(5,1) * t200 - rSges(5,2) * t199 + rSges(5,3) * t338;
t509 = t113 * t308;
t114 = rSges(5,1) * t202 - rSges(5,2) * t201 + rSges(5,3) * t340;
t508 = t114 * t307;
t151 = rSges(5,1) * t302 + rSges(5,2) * t452 + rSges(5,3) * t359;
t507 = t151 * t214;
t152 = rSges(5,1) * t304 + rSges(5,2) * t451 + rSges(5,3) * t361;
t506 = t152 * t213;
t505 = t359 * t412;
t504 = t361 * t412;
t503 = t566 * t412;
t156 = rSges(6,1) * t306 + rSges(6,2) * t305 - rSges(6,3) * t450;
t278 = pkin(4) * t364 - pkin(9) * t450;
t502 = t156 + t278;
t486 = t359 * t532;
t257 = -pkin(4) * t486 - pkin(9) * t505;
t485 = t361 * t532;
t258 = -pkin(4) * t485 - pkin(9) * t504;
t484 = t566 * t532;
t323 = pkin(4) * t484 + pkin(9) * t503;
t405 = qJDD(2) * t523;
t404 = 0.2e1 * t405;
t415 = 2 * qJDD(1);
t501 = t415 + t404;
t500 = qJD(2) * t409;
t499 = qJD(4) * t360;
t498 = qJD(4) * t362;
t497 = qJD(4) * t385;
t496 = qJD(5) * t412;
t495 = qJD(2) * t523 + qJD(1);
t494 = 2 * m(4);
t493 = 0.2e1 * m(5);
t492 = 0.2e1 * m(6);
t491 = qJDD(2) * t409;
t489 = qJDD(1) + t405;
t488 = 0.2e1 * t523;
t487 = 0.2e1 * t409;
t483 = t411 * t532;
t253 = t359 * t483 + t360 * t414;
t482 = t414 * t532;
t254 = -t359 * t482 + t360 * t411;
t133 = t254 * rSges(6,1) + t253 * rSges(6,2) - rSges(6,3) * t505;
t255 = t361 * t483 + t362 * t414;
t256 = -t361 * t482 + t362 * t411;
t134 = t256 * rSges(6,1) + t255 * rSges(6,2) - rSges(6,3) * t504;
t321 = t385 * t414 - t483 * t566;
t322 = t385 * t411 + t482 * t566;
t186 = t322 * rSges(6,1) + t321 * rSges(6,2) + rSges(6,3) * t503;
t481 = t410 * t500;
t480 = t410 * t491;
t269 = -t359 * pkin(3) + pkin(8) * t360;
t271 = -t361 * pkin(3) + pkin(8) * t362;
t329 = pkin(3) * t566 + pkin(8) * t385;
t176 = pkin(4) * t452 + pkin(9) * t302;
t178 = pkin(4) * t451 + pkin(9) * t304;
t277 = pkin(4) * t450 + pkin(9) * t364;
t476 = -t493 / 0.2e1;
t475 = t493 / 0.2e1;
t474 = -t492 / 0.2e1;
t473 = t492 / 0.2e1;
t472 = -g(3) + t404 / 0.2e1 + t415 / 0.2e1;
t471 = t269 * t383 - t271 * t382;
t470 = t271 * t393 - t329 * t383;
t469 = -t269 * t393 + t329 * t382;
t464 = -Icges(6,1) * t414 + Icges(6,4) * t411;
t463 = -Icges(6,4) * t414 + Icges(6,2) * t411;
t462 = -Icges(6,5) * t414 + Icges(6,6) * t411;
t224 = rSges(4,1) * t360 - rSges(4,2) * t359 + rSges(4,3) * t386;
t225 = rSges(4,1) * t362 - rSges(4,2) * t361 + rSges(4,3) * t387;
t461 = t224 * t387 - t225 * t386;
t313 = rSges(4,1) * t385 + rSges(4,2) * t566 + rSges(4,3) * t394;
t460 = -t224 * t394 + t313 * t386;
t459 = t225 * t394 - t313 * t387;
t245 = -rSges(4,1) * t337 - rSges(4,2) * t338;
t246 = -rSges(4,1) * t339 - rSges(4,2) * t340;
t458 = t245 * t387 - t246 * t386;
t317 = rSges(4,1) * t369 - rSges(4,2) * t370;
t457 = -t245 * t394 + t317 * t386;
t456 = t246 * t394 - t317 * t387;
t193 = -rSges(5,1) * t486 + rSges(5,2) * t505 + t360 * rSges(5,3);
t194 = -rSges(5,1) * t485 + rSges(5,2) * t504 + t362 * rSges(5,3);
t285 = rSges(5,1) * t484 - rSges(5,2) * t503 + t385 * rSges(5,3);
t143 = t302 * rSges(6,3) + t452 * t567;
t144 = t304 * rSges(6,3) + t451 * t567;
t198 = t364 * rSges(6,3) + t450 * t567;
t270 = pkin(3) * t360 + pkin(8) * t359;
t330 = pkin(3) * t385 - pkin(8) * t566;
t403 = t408 * t500;
t453 = -t270 * t393 + t330 * t382 + t403;
t449 = -Icges(5,1) * t532 + Icges(5,4) * t412;
t448 = -Icges(5,4) * t532 + Icges(5,2) * t412;
t447 = -Icges(5,5) * t532 + Icges(5,6) * t412;
t272 = pkin(3) * t362 + pkin(8) * t361;
t446 = t270 * t383 - t272 * t382 + t495;
t443 = (Icges(6,5) * t207 - Icges(6,6) * t208) * t181 + (Icges(6,5) * t209 - Icges(6,6) * t210) * t182 + (Icges(6,5) * t305 - Icges(6,6) * t306) * t279;
t442 = (Icges(5,5) * t452 - Icges(5,6) * t302) * t307 + (Icges(5,5) * t451 - Icges(5,6) * t304) * t308 + (Icges(5,5) * t450 - Icges(5,6) * t364) * t365;
t441 = (-Icges(4,5) * t359 - Icges(4,6) * t360) * t386 + (-Icges(4,5) * t361 - Icges(4,6) * t362) * t387 + (Icges(4,5) * t566 - Icges(4,6) * t385) * t394;
t440 = t272 * t393 - t330 * t383 - t481;
t247 = -pkin(3) * t337 + pkin(8) * t338;
t318 = pkin(3) * t369 + pkin(8) * t370;
t402 = t408 * t491;
t431 = t318 * t382 + t330 * t379 + t402 + (-qJD(3) * t247 - qJDD(3) * t270) * t394;
t211 = t247 * t383;
t223 = t270 * t380;
t248 = -pkin(3) * t339 + pkin(8) * t340;
t430 = t211 + t223 + (-qJD(3) * t248 - qJDD(3) * t272) * t386 + t489;
t429 = (-Icges(6,2) * t210 + t204 + t96) * t182 + (-Icges(6,2) * t208 + t203 + t95) * t181 + (-Icges(6,2) * t306 + t155 + t298) * t279;
t428 = (Icges(6,1) * t209 - t511 - t94) * t182 + (Icges(6,1) * t207 - t512 - t93) * t181 + (Icges(6,1) * t305 - t154 - t510) * t279;
t427 = (Icges(5,1) * t451 - t148 - t514) * t308 + (Icges(5,1) * t452 - t147 - t515) * t307 + (Icges(5,1) * t450 - t229 - t513) * t365;
t426 = (Icges(5,2) * t304 - t150 - t295) * t308 + (Icges(5,2) * t302 - t149 - t294) * t307 + (Icges(5,2) * t364 - t230 - t356) * t365;
t219 = -Icges(4,2) * t359 + Icges(4,6) * t386 + t518;
t220 = -Icges(4,2) * t361 + Icges(4,6) * t387 + t517;
t311 = Icges(4,2) * t566 + Icges(4,6) * t394 + t516;
t425 = (-Icges(4,1) * t361 - t220 - t517) * t387 + (-Icges(4,1) * t359 - t219 - t518) * t386 + (Icges(4,1) * t566 - t311 - t516) * t394;
t350 = Icges(4,4) * t359;
t221 = Icges(4,1) * t360 + Icges(4,5) * t386 - t350;
t351 = Icges(4,4) * t361;
t222 = Icges(4,1) * t362 + Icges(4,5) * t387 - t351;
t377 = Icges(4,4) * t566;
t312 = Icges(4,1) * t385 + Icges(4,5) * t394 + t377;
t424 = (Icges(4,2) * t362 - t222 + t351) * t387 + (Icges(4,2) * t360 - t221 + t350) * t386 + (Icges(4,2) * t385 - t312 - t377) * t394;
t421 = t248 * t393 + t272 * t392 + (-qJD(3) * t318 - qJDD(3) * t330) * t387 - t480;
t420 = qJD(3) * t458 + qJDD(3) * t461;
t418 = t100 * t97 + t117 * t308 - t118 * t307 + t177 * t214 - t179 * t213 - t181 * t61 + t182 * t60 - t98 * t99;
t417 = (Icges(6,3) * t304 + t411 * t94 - t414 * t96 - t451 * t462) * t182 + (Icges(6,3) * t302 + t411 * t93 - t414 * t95 - t452 * t462) * t181 + (Icges(6,3) * t364 + t154 * t411 - t155 * t414 - t450 * t462) * t279;
t416 = (Icges(5,3) * t362 + t148 * t412 - t150 * t532 + t361 * t447) * t308 + (Icges(5,3) * t360 + t147 * t412 - t149 * t532 + t359 * t447) * t307 + (Icges(5,3) * t385 + t229 * t412 - t230 * t532 - t447 * t566) * t365;
t328 = rSges(4,1) * t566 - rSges(4,2) * t385;
t324 = t496 * t566 + t497;
t316 = Icges(4,1) * t369 - Icges(4,4) * t370;
t315 = Icges(4,4) * t369 - Icges(4,2) * t370;
t314 = Icges(4,5) * t369 - Icges(4,6) * t370;
t310 = Icges(4,5) * t385 + Icges(4,6) * t566 + Icges(4,3) * t394;
t291 = t386 * t330;
t284 = Icges(5,5) * t385 - t449 * t566;
t283 = Icges(5,6) * t385 - t448 * t566;
t281 = t386 * t318;
t276 = rSges(5,1) * t450 - rSges(5,2) * t364;
t268 = -rSges(4,1) * t361 - rSges(4,2) * t362;
t267 = -rSges(4,1) * t359 - rSges(4,2) * t360;
t260 = -t361 * t496 + t498;
t259 = -t359 * t496 + t499;
t244 = -Icges(4,1) * t339 - Icges(4,4) * t340;
t243 = -Icges(4,1) * t337 - Icges(4,4) * t338;
t242 = -Icges(4,4) * t339 - Icges(4,2) * t340;
t241 = -Icges(4,4) * t337 - Icges(4,2) * t338;
t240 = -Icges(4,5) * t339 - Icges(4,6) * t340;
t239 = -Icges(4,5) * t337 - Icges(4,6) * t338;
t238 = t394 * t272;
t234 = t394 * t248;
t232 = t387 * t270;
t231 = rSges(5,1) * t364 + rSges(5,2) * t450 - rSges(5,3) * t566;
t218 = Icges(4,5) * t362 - Icges(4,6) * t361 + Icges(4,3) * t387;
t217 = Icges(4,5) * t360 - Icges(4,6) * t359 + Icges(4,3) * t386;
t216 = 0.2e1 * t223;
t215 = -0.2e1 * t272 * t379;
t212 = t387 * t247;
t206 = 0.2e1 * t211;
t205 = -0.2e1 * t248 * t382;
t197 = Icges(6,5) * t364 - t450 * t464;
t196 = Icges(6,6) * t364 - t450 * t463;
t192 = Icges(5,5) * t362 + t361 * t449;
t191 = Icges(5,5) * t360 + t359 * t449;
t190 = Icges(5,6) * t362 + t361 * t448;
t189 = Icges(5,6) * t360 + t359 * t448;
t185 = Icges(6,1) * t322 + Icges(6,4) * t321 + Icges(6,5) * t503;
t184 = Icges(6,4) * t322 + Icges(6,2) * t321 + Icges(6,6) * t503;
t183 = Icges(6,5) * t322 + Icges(6,6) * t321 + Icges(6,3) * t503;
t180 = rSges(6,1) * t305 - rSges(6,2) * t306;
t172 = rSges(5,1) * t451 - rSges(5,2) * t304;
t171 = rSges(5,1) * t452 - rSges(5,2) * t302;
t163 = rSges(5,1) * t293 - rSges(5,2) * t292 + rSges(5,3) * t370;
t142 = Icges(6,5) * t304 - t451 * t464;
t141 = Icges(6,5) * t302 - t452 * t464;
t140 = Icges(6,6) * t304 - t451 * t463;
t139 = Icges(6,6) * t302 - t452 * t463;
t136 = qJD(3) * t459 - t481;
t135 = qJD(3) * t460 + t403;
t132 = Icges(6,1) * t256 + Icges(6,4) * t255 - Icges(6,5) * t504;
t131 = Icges(6,1) * t254 + Icges(6,4) * t253 - Icges(6,5) * t505;
t130 = Icges(6,4) * t256 + Icges(6,2) * t255 - Icges(6,6) * t504;
t129 = Icges(6,4) * t254 + Icges(6,2) * t253 - Icges(6,6) * t505;
t128 = Icges(6,5) * t256 + Icges(6,6) * t255 - Icges(6,3) * t504;
t127 = Icges(6,5) * t254 + Icges(6,6) * t253 - Icges(6,3) * t505;
t126 = rSges(6,1) * t209 - rSges(6,2) * t210;
t125 = rSges(6,1) * t207 - rSges(6,2) * t208;
t115 = qJD(3) * t461 + t495;
t88 = qJD(3) * t456 + qJDD(3) * t459 - t480;
t87 = qJD(3) * t457 + qJDD(3) * t460 + t402;
t74 = t420 + t489;
t71 = t152 * t365 - t231 * t308 + t440;
t70 = -t151 * t365 + t231 * t307 + t453;
t53 = t151 * t308 - t152 * t307 + t446;
t38 = -t156 * t182 + t179 * t365 - t278 * t308 + t279 * t98 + t440;
t37 = t156 * t181 - t177 * t365 + t278 * t307 - t279 * t97 + t453;
t34 = t114 * t365 + t152 * t309 - t163 * t308 - t214 * t231 + t421;
t33 = -t113 * t365 - t151 * t309 + t163 * t307 + t213 * t231 + t431;
t32 = t177 * t308 - t179 * t307 - t181 * t98 + t182 * t97 + t446;
t29 = t430 - t506 + t507 - t508 + t509;
t24 = t307 * t66 + t308 * t67 + t365 * t81;
t23 = t307 * t50 + t308 * t51 + t365 * t73;
t22 = t307 * t48 + t308 * t49 + t365 * t72;
t18 = t181 * t44 + t182 * t45 + t279 * t52;
t17 = -t100 * t156 + t118 * t365 + t157 * t98 - t164 * t308 + t179 * t309 - t182 * t85 - t214 * t278 + t279 * t61 + t421;
t16 = -t117 * t365 + t156 * t99 - t157 * t97 + t164 * t307 - t177 * t309 + t181 * t85 + t213 * t278 - t279 * t60 + t431;
t13 = t181 * t41 + t182 * t42 + t279 * t47;
t12 = t181 * t39 + t182 * t40 + t279 * t46;
t7 = t418 + t430;
t4 = [(m(2) * qJDD(1)) + t501 * t563 + (t205 + t206 + t215 + t216 + t501 - 0.2e1 * t506 + 0.2e1 * t507 - 0.2e1 * t508 + 0.2e1 * t509) * t562 + (-m(2) + t565) * g(3) + (t205 / 0.2e1 + t206 / 0.2e1 + t215 / 0.2e1 + t216 / 0.2e1 + t418 + t472) * m(6) + (t420 + t472) * m(4); (t489 * t488 + 0.2e1 * (t408 ^ 2 + t410 ^ 2) * t409 ^ 2 * qJDD(2)) * t563 + m(4) * (t74 * t488 + (t408 * t87 - t410 * t88) * t487) / 0.2e1 + (t29 * t488 + (t33 * t408 - t34 * t410) * t487) * t562 + (t7 * t488 + (t16 * t408 - t17 * t410) * t487) * t561 + (-m(4) - m(6) + t565) * (g(3) * t523 + (g(1) * t408 - g(2) * t410) * t409); -(t394 * (t385 * t425 + t394 * t441 - t424 * t566) + t387 * (t361 * t424 + t362 * t425 + t387 * t441) + t386 * (t359 * t424 + t360 * t425 + t386 * t441)) * (qJD(3) ^ 2) / 0.2e1 + (((t310 * t394 + t311 * t566 + t312 * t385) * t394 + t386 * (t217 * t394 + t219 * t566 + t221 * t385) + t387 * (t218 * t394 + t220 * t566 + t222 * t385)) * t571 + (t386 * (-t219 * t370 + t221 * t369 + t239 * t394 + t241 * t566 + t243 * t385) + t387 * (-t220 * t370 + t222 * t369 + t240 * t394 + t242 * t566 + t244 * t385) + t394 * (-t311 * t370 + t312 * t369 + t314 * t394 + t315 * t566 + t316 * t385)) * t572 + t568) * t394 / 0.2e1 + ((t385 * t146 + t190 * t450 + t364 * t192) * t308 + (t385 * t145 + t189 * t450 + t364 * t191) * t307 + (t385 * t228 + t283 * t450 + t364 * t284) * t365 + (t360 * t66 + t362 * t67 + t385 * t81) * qJD(4) - t416 * t566) * t540 - m(5) * (g(1) * (t194 + t271) + g(2) * (t193 + t269) + g(3) * (t285 + t329)) - t259 * t12 / 0.2e1 - t260 * t13 / 0.2e1 + (t115 * t458 + t135 * t457 + t136 * t456 + t459 * t88 + t460 * t87 + t461 * t74) * t494 / 0.2e1 + ((-t128 * t450 + t130 * t305 + t132 * t306 + t321 * t94 + t322 * t96 + t503 * t92) * t182 + t45 * t260 + (-t127 * t450 + t129 * t305 + t131 * t306 + t321 * t93 + t322 * t95 + t503 * t91) * t181 + t44 * t259 + (t153 * t503 + t154 * t321 + t155 * t322 - t183 * t450 + t184 * t305 + t185 * t306) * t279 + t52 * t324) * t551 + ((-t128 * t452 + t130 * t207 + t132 * t208 + t253 * t94 + t254 * t96 - t505 * t92) * t182 + t40 * t260 + (-t127 * t452 + t129 * t207 + t131 * t208 + t253 * t93 + t254 * t95 - t505 * t91) * t181 + t39 * t259 + (-t153 * t505 + t154 * t253 + t155 * t254 - t183 * t452 + t184 * t207 + t185 * t208) * t279 + t46 * t324) * t557 + ((t360 * t146 + t190 * t452 + t302 * t192) * t308 + (t360 * t145 + t189 * t452 + t302 * t191) * t307 + (t360 * t228 + t283 * t452 + t302 * t284) * t365 + (t360 * t48 + t362 * t49 + t385 * t72) * qJD(4) + t416 * t359) * t549 - t22 * t499 / 0.2e1 - t324 * t18 / 0.2e1 - m(4) * (g(1) * t268 + g(2) * t267 + g(3) * t328) + (t19 * t394 + t386 * t8 + t387 * t9) * t556 + (t386 * t44 + t387 * t45 + t394 * t52) * t558 + (t386 * t41 + t387 * t42 + t394 * t47) * t559 + (t386 * t39 + t387 * t40 + t394 * t46) * t560 + (t37 * (-t133 * t279 + t156 * t259 + t181 * t186 - t257 * t365 + t307 * t323 - t324 * t97 + (-t177 * t385 + t278 * t360) * qJD(4) + t469) + t38 * (t134 * t279 - t156 * t260 - t182 * t186 + t258 * t365 - t308 * t323 + t324 * t98 + (t179 * t385 - t278 * t362) * qJD(4) + t470) + t32 * (t133 * t182 - t134 * t181 + t257 * t308 - t258 * t307 - t259 * t98 + t260 * t97 + (t177 * t362 - t179 * t360) * qJD(4) + t471)) * t474 + (t30 * t386 + t31 * t387 + t394 * t43) * t539 + (t386 * t66 + t387 * t67 + t394 * t81) * t545 + (t27 * t386 + t28 * t387 + t36 * t394) * t546 + (t25 * t386 + t26 * t387 + t35 * t394) * t548 + (t14 * t386 + t15 * t387 + t21 * t394) * t550 + (t386 * t50 + t387 * t51 + t394 * t73) * t552 + (t386 * t48 + t387 * t49 + t394 * t72) * t553 + (t10 * t386 + t11 * t387 + t20 * t394) * t554 - t24 * t497 / 0.2e1 - t23 * t498 / 0.2e1 + (((t310 * t387 - t311 * t361 + t312 * t362) * t394 + t386 * (t217 * t387 - t219 * t361 + t221 * t362) + t387 * (t218 * t387 - t220 * t361 + t222 * t362)) * t571 + (t386 * (-t219 * t340 - t221 * t339 + t239 * t387 - t241 * t361 + t243 * t362) + t387 * (-t220 * t340 - t222 * t339 + t240 * t387 - t242 * t361 + t244 * t362) + t394 * (-t311 * t340 - t312 * t339 + t314 * t387 - t315 * t361 + t316 * t362)) * t572 + t569) * t387 / 0.2e1 + (((t310 * t386 - t311 * t359 + t312 * t360) * t394 + t386 * (t217 * t386 - t219 * t359 + t221 * t360) + t387 * (t218 * t386 - t220 * t359 + t222 * t360)) * t571 + (t386 * (-t219 * t338 - t221 * t337 + t239 * t386 - t241 * t359 + t243 * t360) + t387 * (-t220 * t338 - t222 * t337 + t240 * t386 - t242 * t359 + t244 * t360) + t394 * (-t311 * t338 - t312 * t337 + t314 * t386 - t315 * t359 + t316 * t360)) * t572 + t570) * t386 / 0.2e1 - (t135 * (-t267 * t394 + t328 * t386) + t136 * (t268 * t394 - t328 * t387) + t115 * (t267 * t387 - t268 * t386)) * qJD(3) * t494 / 0.2e1 - m(6) * (g(1) * (t271 + t134 + t258) + g(2) * (t269 + t133 + t257) + g(3) * (t329 + t186 + t323)) + ((t362 * t146 + t190 * t451 + t304 * t192) * t308 + (t362 * t145 + t189 * t451 + t304 * t191) * t307 + (t362 * t228 + t283 * t451 + t304 * t284) * t365 + (t360 * t50 + t362 * t51 + t385 * t73) * qJD(4) + t416 * t361) * t547 + ((-t128 * t451 + t130 * t209 + t132 * t210 + t255 * t94 + t256 * t96 - t504 * t92) * t182 + t42 * t260 + (-t127 * t451 + t129 * t209 + t131 * t210 + t255 * t93 + t256 * t95 - t504 * t91) * t181 + t41 * t259 + (-t153 * t504 + t154 * t255 + t155 * t256 - t183 * t451 + t184 * t209 + t185 * t210) * t279 + t47 * t324) * t555 + 0.2e1 * (t16 * t291 + t17 * t238 + t32 * t212 + t7 * t232 + t38 * t234 + t37 * t281 + (t16 * (-t270 - t525) + t37 * (-t247 - t528) + t17 * t524 + t38 * t527) * t394 + (t17 * (-t330 - t502) + t38 * (-t318 - t526) + t7 * t525 + t32 * t528) * t387 + (t16 * t502 + t37 * t526 + t7 * (-t272 - t524) + t32 * (-t248 - t527)) * t386) * t561 + (t33 * (t231 * t386 + t291 + (-t151 - t270) * t394) + t70 * (t163 * t386 + t281 + (-t113 - t247) * t394) + t34 * (t152 * t394 + t238 + (-t231 - t330) * t387) + t71 * (t114 * t394 + t234 + (-t163 - t318) * t387) + t29 * (t151 * t387 + t232 + (-t152 - t272) * t386) + t53 * (t113 * t387 + t212 + (-t114 - t248) * t386)) * t475 + (t70 * (-t193 * t365 + t285 * t307 + (-t151 * t385 + t231 * t360) * qJD(4) + t469) + t71 * (t194 * t365 - t285 * t308 + (t152 * t385 - t231 * t362) * qJD(4) + t470) + t53 * (t193 * t308 - t194 * t307 + (t151 * t362 - t152 * t360) * qJD(4) + t471)) * t476; (t16 * (t359 * t502 + t525 * t566) + t37 * (t338 * t502 + t359 * t526 - t370 * t525 + t528 * t566) + t17 * (-t361 * t502 - t524 * t566) + t38 * (-t340 * t502 - t361 * t526 + t370 * t524 - t527 * t566) + t7 * (-t359 * t524 + t361 * t525) + t32 * (-t338 * t524 + t340 * t525 - t359 * t527 + t361 * t528)) * t473 + (t33 * (t151 * t566 + t231 * t359) + t70 * (t113 * t566 - t151 * t370 + t163 * t359 + t231 * t338) + t34 * (-t152 * t566 - t231 * t361) + t71 * (-t114 * t566 + t152 * t370 - t163 * t361 - t231 * t340) + t29 * (t151 * t361 - t152 * t359) + t53 * (t113 * t361 - t114 * t359 + t151 * t340 - t152 * t338)) * t475 - t568 * t566 / 0.2e1 + (-t19 * t566 + t338 * t39 + t340 * t40 + t359 * t8 + t361 * t9 + t370 * t46) * t556 + (t359 * t44 + t361 * t45 - t52 * t566) * t558 + (t359 * t41 + t361 * t42 - t47 * t566) * t559 + (t359 * t39 + t361 * t40 - t46 * t566) * t560 + (t30 * t359 + t31 * t361 + t338 * t66 + t340 * t67 + t370 * t81 - t43 * t566) * t539 + (t359 * t66 + t361 * t67 - t566 * t81) * t545 + (t27 * t359 + t28 * t361 + t338 * t50 + t340 * t51 - t36 * t566 + t370 * t73) * t546 + (t25 * t359 + t26 * t361 + t338 * t48 + t340 * t49 - t35 * t566 + t370 * t72) * t548 + (t14 * t359 + t15 * t361 - t21 * t566 + t338 * t44 + t340 * t45 + t370 * t52) * t550 + (t359 * t50 + t361 * t51 - t566 * t73) * t552 + (t359 * t48 + t361 * t49 - t566 * t72) * t553 + (t10 * t359 + t11 * t361 - t20 * t566 + t338 * t41 + t340 * t42 + t370 * t47) * t554 + (t364 * t427 - t426 * t450 - t442 * t566) * t540 - m(5) * (g(1) * t172 + g(2) * t171 + g(3) * t276) - (t12 * t302 + t13 * t304 + t18 * t364) * qJD(5) / 0.2e1 + (t18 + t24) * t370 / 0.2e1 + (t23 + t13) * t340 / 0.2e1 + (t22 + t12) * t338 / 0.2e1 + ((t140 * t207 + t142 * t208 + t302 * t92) * t182 + (t139 * t207 + t141 * t208 + t302 * t91) * t181 + (t153 * t302 + t196 * t207 + t197 * t208) * t279 + (t302 * t39 + t304 * t40 + t364 * t46) * qJD(5) - t417 * t452) * t557 + (t302 * t427 + t359 * t442 - t426 * t452) * t549 + t569 * t361 / 0.2e1 + t570 * t359 / 0.2e1 + (t37 * (-t143 * t279 - t176 * t365 + t181 * t198 + t277 * t307 + (t156 * t302 - t364 * t97) * qJD(5)) + t38 * (t144 * t279 + t178 * t365 - t182 * t198 - t277 * t308 + (-t156 * t304 + t364 * t98) * qJD(5)) + t32 * (t143 * t182 - t144 * t181 + t176 * t308 - t178 * t307 + (-t302 * t98 + t304 * t97) * qJD(5))) * t474 - m(6) * (g(1) * (t144 + t178) + g(2) * (t143 + t176) + g(3) * (t198 + t277)) + ((t140 * t209 + t142 * t210 + t304 * t92) * t182 + (t139 * t209 + t141 * t210 + t304 * t91) * t181 + (t153 * t304 + t196 * t209 + t197 * t210) * t279 + (t302 * t41 + t304 * t42 + t364 * t47) * qJD(5) - t417 * t451) * t555 + (t304 * t427 + t361 * t442 - t426 * t451) * t547 + ((t140 * t305 + t142 * t306 + t364 * t92) * t182 + (t139 * t305 + t141 * t306 + t364 * t91) * t181 + (t153 * t364 + t196 * t305 + t197 * t306) * t279 + (t302 * t44 + t304 * t45 + t364 * t52) * qJD(5) - t417 * t450) * t551 + (t70 * (-t171 * t365 + t276 * t307) + t71 * (t172 * t365 - t276 * t308) + t53 * (t171 * t308 - t172 * t307)) * t476; (t16 * (-t156 * t452 + t450 * t97) + t37 * (t156 * t199 - t292 * t97 + t450 * t60 - t452 * t85) + t17 * (t156 * t451 - t450 * t98) + t38 * (-t156 * t201 + t292 * t98 - t450 * t61 + t451 * t85) + t7 * (-t451 * t97 + t452 * t98) + t32 * (-t199 * t98 + t201 * t97 - t451 * t60 + t452 * t61)) * t473 + t201 * t13 / 0.2e1 - t451 * t2 / 0.2e1 + (-t41 * t452 - t42 * t451 - t450 * t47) * t559 + (-t10 * t452 - t11 * t451 + t199 * t41 - t20 * t450 + t201 * t42 + t292 * t47) * t554 + t199 * t12 / 0.2e1 - t452 * t1 / 0.2e1 + (-t39 * t452 - t40 * t451 - t450 * t46) * t560 + (-t19 * t450 + t199 * t39 + t201 * t40 + t292 * t46 - t451 * t9 - t452 * t8) * t556 + t292 * t18 / 0.2e1 - t450 * t3 / 0.2e1 + (-t44 * t452 - t45 * t451 - t450 * t52) * t558 + (-t14 * t452 - t15 * t451 + t199 * t44 + t201 * t45 - t21 * t450 + t292 * t52) * t550 + (t37 * (-t125 * t279 + t180 * t181) + t38 * (t126 * t279 - t180 * t182) + t32 * (t125 * t182 - t126 * t181)) * t474 + (t209 * t429 + t210 * t428 - t443 * t451) * t555 + (t207 * t429 + t208 * t428 - t443 * t452) * t557 + (t305 * t429 + t306 * t428 - t443 * t450) * t551 - m(6) * (g(1) * t126 + g(2) * t125 + g(3) * t180);];
tau = t4;
