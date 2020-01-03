% Calculate vector of inverse dynamics joint torques for
% S5RPPRR12
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR12_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR12_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:57
% EndTime: 2019-12-31 18:07:20
% DurationCPUTime: 18.51s
% Computational Cost: add. (13835->818), mult. (20089->1077), div. (0->0), fcn. (18388->8), ass. (0->396)
t308 = sin(qJ(1));
t557 = pkin(1) * t308;
t310 = cos(qJ(1));
t556 = t310 * pkin(1);
t303 = pkin(8) + qJ(4);
t282 = sin(t303);
t309 = cos(qJ(5));
t462 = t308 * t309;
t307 = sin(qJ(5));
t464 = t307 * t310;
t209 = t282 * t464 + t462;
t461 = t309 * t310;
t463 = t308 * t307;
t210 = t282 * t461 - t463;
t443 = t210 * rSges(6,1) - t209 * rSges(6,2);
t283 = cos(t303);
t466 = t283 * t310;
t117 = rSges(6,3) * t466 - t443;
t416 = qJD(1) * qJD(4);
t233 = qJDD(4) * t308 + t310 * t416;
t424 = qJD(4) * t310;
t398 = t282 * t424;
t428 = qJD(1) * t308;
t329 = -t283 * t428 - t398;
t414 = qJDD(5) * t283;
t120 = qJD(5) * t329 + t310 * t414 + t233;
t396 = t283 * t424;
t122 = t329 * pkin(7) + (t282 * t428 - t396) * pkin(4);
t494 = rSges(6,2) * t307;
t498 = rSges(6,1) * t309;
t374 = -t494 + t498;
t492 = rSges(6,3) * t282;
t160 = t283 * t374 + t492;
t422 = qJD(5) * t283;
t205 = qJD(4) * t422 + qJDD(5) * t282 + qJDD(1);
t273 = t283 * pkin(7);
t507 = pkin(4) * t282;
t227 = t273 - t507;
t216 = t227 * qJD(4);
t420 = qJD(5) * t310;
t395 = t283 * t420;
t425 = qJD(4) * t308;
t217 = t395 + t425;
t504 = t283 * pkin(4);
t228 = pkin(7) * t282 + t504;
t423 = qJD(5) * t282;
t257 = qJD(1) + t423;
t417 = qJD(1) * qJD(3);
t418 = qJD(1) * qJD(2);
t436 = qJDD(2) * t308 + t310 * t418;
t542 = -qJ(3) * qJD(1) ^ 2 + qJDD(3);
t328 = -0.2e1 * t308 * t417 + t310 * t542 + t436;
t261 = pkin(7) * t466;
t468 = t282 * t310;
t197 = pkin(4) * t468 - t261;
t306 = -pkin(6) - qJ(3);
t459 = qJ(3) + t306;
t304 = sin(pkin(8));
t508 = pkin(3) * t304;
t200 = t308 * t459 + t310 * t508;
t293 = t310 * qJ(2);
t244 = -t293 + t557;
t474 = qJ(3) * t308;
t387 = -t244 - t474;
t380 = t200 + t387;
t340 = t197 + t380;
t247 = t308 * qJ(2) + t556;
t291 = qJD(2) * t310;
t214 = qJD(1) * t247 - t291;
t427 = qJD(1) * t310;
t266 = t306 * t427;
t465 = t304 * t308;
t269 = pkin(3) * t465;
t473 = qJ(3) * t310;
t445 = t266 - (t269 - t473) * qJD(1) - t214;
t342 = t257 * t310;
t384 = qJD(1) * t282 + qJD(5);
t532 = t308 * t384 - t396;
t91 = -t307 * t532 + t309 * t342;
t92 = t307 * t342 + t309 * t532;
t382 = rSges(6,1) * t92 + rSges(6,2) * t91;
t50 = rSges(6,3) * t329 + t382;
t188 = (-rSges(6,1) * t307 - rSges(6,2) * t309) * t283;
t272 = t283 * rSges(6,3);
t95 = qJD(5) * t188 + (-t282 * t374 + t272) * qJD(4);
t10 = t216 * t425 - t117 * t205 + t120 * t160 + t217 * t95 + t228 * t233 - t257 * t50 + (-t122 + t445) * qJD(1) + t340 * qJDD(1) + t328;
t555 = t10 * t310;
t554 = t117 * t257 - t160 * t217 - t228 * t425;
t207 = -t282 * t463 + t461;
t208 = t282 * t462 + t464;
t467 = t283 * t308;
t106 = Icges(6,5) * t208 + Icges(6,6) * t207 - Icges(6,3) * t467;
t108 = -Icges(6,5) * t210 + Icges(6,6) * t209 + Icges(6,3) * t466;
t186 = Icges(6,4) * t210;
t111 = Icges(6,2) * t209 + Icges(6,6) * t466 - t186;
t185 = Icges(6,4) * t209;
t113 = Icges(6,1) * t210 - Icges(6,5) * t466 - t185;
t352 = -t111 * t209 - t113 * t210;
t481 = Icges(6,4) * t208;
t109 = Icges(6,2) * t207 - Icges(6,6) * t467 + t481;
t184 = Icges(6,4) * t207;
t112 = Icges(6,1) * t208 - Icges(6,5) * t467 + t184;
t499 = t207 * t109 + t208 * t112;
t551 = t352 + t499 + (-t106 * t308 - t108 * t310) * t283;
t421 = qJD(5) * t308;
t218 = -t283 * t421 + t424;
t29 = t106 * t466 + t209 * t109 - t210 * t112;
t355 = Icges(6,5) * t309 - Icges(6,6) * t307;
t154 = Icges(6,3) * t282 + t283 * t355;
t479 = Icges(6,4) * t309;
t357 = -Icges(6,2) * t307 + t479;
t156 = Icges(6,6) * t282 + t283 * t357;
t480 = Icges(6,4) * t307;
t359 = Icges(6,1) * t309 - t480;
t158 = Icges(6,5) * t282 + t283 * t359;
t53 = t154 * t466 + t156 * t209 - t158 * t210;
t550 = -t218 * t29 - t257 * t53;
t351 = t111 * t307 + t113 * t309;
t42 = t108 * t282 - t283 * t351;
t28 = -t108 * t467 + t207 * t111 - t113 * t208;
t397 = t283 * t425;
t545 = t282 * t427 + t397;
t250 = Icges(5,4) * t467;
t469 = t282 * t308;
t478 = Icges(5,5) * t310;
t166 = Icges(5,1) * t469 + t250 + t478;
t482 = Icges(5,4) * t283;
t360 = Icges(5,1) * t282 + t482;
t167 = -Icges(5,5) * t308 + t310 * t360;
t222 = -Icges(5,2) * t282 + t482;
t180 = t222 * t310;
t319 = t308 * (t167 + t180) - t310 * (-Icges(5,2) * t469 + t166 + t250);
t483 = Icges(5,4) * t282;
t358 = Icges(5,2) * t283 + t483;
t164 = Icges(5,6) * t310 + t308 * t358;
t165 = -Icges(5,6) * t308 + t310 * t358;
t224 = Icges(5,1) * t283 - t483;
t182 = t224 * t308;
t183 = t224 * t310;
t320 = t308 * (t165 - t183) - t310 * (t164 - t182);
t544 = -t320 * t282 + t319 * t283;
t441 = t222 + t360;
t442 = -t358 + t224;
t543 = (t282 * t441 - t283 * t442) * qJD(1);
t52 = -t154 * t467 + t156 * t207 + t158 * t208;
t540 = t217 * t28 + t52 * t257;
t509 = rSges(6,3) + pkin(7);
t539 = t282 * t509;
t153 = Icges(6,3) * t283 - t282 * t355;
t349 = t156 * t307 - t158 * t309;
t353 = t109 * t307 - t112 * t309;
t312 = t217 * (-t154 * t310 + t351) + t218 * (t154 * t308 + t353) + t257 * (t153 + t349);
t538 = t312 * t283;
t346 = t165 * t283 + t167 * t282;
t535 = t346 * t310;
t305 = cos(pkin(8));
t496 = rSges(4,2) * t305;
t173 = rSges(4,1) * t465 + t310 * rSges(4,3) + t308 * t496;
t386 = t247 + t473;
t534 = t173 + t386;
t248 = -rSges(3,2) * t310 + t308 * rSges(3,3);
t100 = -qJD(4) * t183 + (t308 * t360 + t478) * qJD(1);
t356 = Icges(5,5) * t282 + Icges(5,6) * t283;
t163 = -Icges(5,3) * t308 + t310 * t356;
t430 = qJD(1) * t163;
t80 = t165 * t282 - t167 * t283;
t98 = qJD(1) * t164 - qJD(4) * t180;
t533 = qJD(4) * t80 + t100 * t282 + t283 * t98 + t430;
t212 = t358 * qJD(4);
t213 = t360 * qJD(4);
t220 = Icges(5,5) * t283 - Icges(5,6) * t282;
t343 = t222 * t282 - t224 * t283;
t531 = qJD(1) * t220 + qJD(4) * t343 + t212 * t283 + t213 * t282;
t101 = qJD(1) * t167 + qJD(4) * t182;
t347 = t164 * t282 - t166 * t283;
t162 = Icges(5,3) * t310 + t308 * t356;
t431 = qJD(1) * t162;
t99 = qJD(1) * t165 + t222 * t425;
t530 = qJD(4) * t347 - t101 * t282 - t283 * t99 + t431;
t178 = (-Icges(6,2) * t309 - t480) * t283;
t314 = t217 * (Icges(6,2) * t210 - t113 + t185) + t218 * (-Icges(6,2) * t208 + t112 + t184) + t257 * (t158 + t178);
t181 = (-Icges(6,1) * t307 - t479) * t283;
t528 = t217 * (-Icges(6,1) * t209 + t111 - t186) + t218 * (-Icges(6,1) * t207 + t109 + t481) + t257 * (t156 - t181);
t27 = -t106 * t467 + t499;
t12 = t218 * t27 + t540;
t527 = -t12 / 0.2e1;
t526 = t120 / 0.2e1;
t285 = qJDD(4) * t310;
t121 = -qJD(1) * t395 + t285 + (-t414 + (-qJD(1) + t423) * qJD(4)) * t308;
t525 = t121 / 0.2e1;
t524 = t205 / 0.2e1;
t523 = -t217 / 0.2e1;
t522 = t217 / 0.2e1;
t521 = -t218 / 0.2e1;
t520 = t218 / 0.2e1;
t519 = t233 / 0.2e1;
t234 = -t308 * t416 + t285;
t518 = t234 / 0.2e1;
t517 = -t257 / 0.2e1;
t516 = t257 / 0.2e1;
t515 = t282 / 0.2e1;
t514 = t308 / 0.2e1;
t513 = -t310 / 0.2e1;
t511 = rSges(3,2) - pkin(1);
t510 = -rSges(5,3) - pkin(1);
t506 = g(1) * t310;
t505 = g(2) * t310;
t399 = t282 * t425;
t401 = t283 * t427;
t330 = t399 - t401;
t93 = -t257 * t462 + (-t310 * t384 - t397) * t307;
t426 = qJD(4) * t283;
t94 = t384 * t461 + (-t257 * t307 + t309 * t426) * t308;
t45 = Icges(6,5) * t94 + Icges(6,6) * t93 + Icges(6,3) * t330;
t47 = Icges(6,4) * t94 + Icges(6,2) * t93 + Icges(6,6) * t330;
t49 = Icges(6,1) * t94 + Icges(6,4) * t93 + Icges(6,5) * t330;
t7 = (qJD(4) * t353 + t45) * t282 + (qJD(4) * t106 - t307 * t47 + t309 * t49 + (-t109 * t309 - t112 * t307) * qJD(5)) * t283;
t503 = t7 * t218;
t44 = Icges(6,5) * t92 + Icges(6,6) * t91 + Icges(6,3) * t329;
t46 = Icges(6,4) * t92 + Icges(6,2) * t91 + Icges(6,6) * t329;
t48 = Icges(6,1) * t92 + Icges(6,4) * t91 + Icges(6,5) * t329;
t8 = (qJD(4) * t351 + t44) * t282 + (qJD(4) * t108 - t307 * t46 + t309 * t48 + (-t111 * t309 + t113 * t307) * qJD(5)) * t283;
t502 = t8 * t217;
t501 = -pkin(1) - qJ(3);
t175 = (-Icges(6,5) * t307 - Icges(6,6) * t309) * t283;
t88 = qJD(4) * t153 + qJD(5) * t175;
t155 = Icges(6,6) * t283 - t282 * t357;
t89 = qJD(4) * t155 + qJD(5) * t178;
t157 = Icges(6,5) * t283 - t282 * t359;
t90 = qJD(4) * t157 + qJD(5) * t181;
t18 = (qJD(4) * t349 + t88) * t282 + (qJD(4) * t154 - t307 * t89 + t309 * t90 + (-t156 * t309 - t158 * t307) * qJD(5)) * t283;
t57 = t154 * t282 - t283 * t349;
t500 = t18 * t257 + t57 * t205;
t495 = rSges(5,2) * t282;
t493 = rSges(3,3) * t310;
t226 = rSges(5,1) * t283 - t495;
t189 = t226 * t310;
t300 = t310 * rSges(5,3);
t376 = rSges(5,1) * t282 + rSges(5,2) * t283;
t104 = -qJD(4) * t189 + (t308 * t376 + t300) * qJD(1);
t215 = t376 * qJD(4);
t297 = t308 * rSges(5,3);
t169 = t310 * t376 - t297;
t341 = t169 + t380;
t25 = -t215 * t425 + t226 * t233 + (-t104 + t445) * qJD(1) + t341 * qJDD(1) + t328;
t490 = t25 * t308;
t407 = rSges(5,1) * t545 + rSges(5,2) * t401;
t105 = (-rSges(5,3) * qJD(1) - qJD(4) * t495) * t308 + t407;
t168 = rSges(5,1) * t469 + rSges(5,2) * t467 + t300;
t201 = -t310 * t459 + t269;
t278 = qJ(2) * t427;
t290 = qJD(2) * t308;
t435 = t278 + t290;
t408 = qJD(1) * (-pkin(1) * t428 + t435) + qJDD(1) * t247 + t308 * t418;
t317 = qJDD(1) * t473 + t308 * t542 + 0.2e1 * t310 * t417 + t408;
t400 = t304 * t427;
t437 = pkin(3) * t400 + t306 * t428;
t316 = qJD(1) * (qJ(3) * t428 + t437) + qJDD(1) * t201 + t317;
t26 = qJD(1) * t105 + qJDD(1) * t168 - t226 * t234 + (qJD(4) * t215 - qJDD(2)) * t310 + t316;
t488 = t26 * t310;
t199 = t226 * t425;
t289 = qJD(3) * t310;
t433 = t289 + t290;
t65 = qJD(1) * t341 + t199 + t433;
t487 = t310 * t65;
t41 = t106 * t282 - t283 * t353;
t486 = t41 * t121;
t485 = t42 * t120;
t484 = t216 + t95;
t470 = t220 * t310;
t176 = t308 * t220;
t344 = t222 * t283 + t224 * t282;
t78 = t310 * t344 - t176;
t460 = t78 * qJD(1);
t444 = t208 * rSges(6,1) + t207 * rSges(6,2);
t115 = -rSges(6,3) * t467 + t444;
t259 = pkin(4) * t469;
t195 = -pkin(7) * t467 + t259;
t454 = -t115 - t195;
t453 = t117 - t197;
t450 = t160 + t228;
t245 = rSges(3,2) * t308 + t493;
t440 = -t244 + t245;
t204 = t247 + t248;
t439 = t282 * t494 + t272;
t196 = pkin(4) * t467 + pkin(7) * t469;
t438 = rSges(4,1) * t400 + t427 * t496;
t434 = rSges(3,2) * t428 + rSges(3,3) * t427;
t432 = -qJD(1) * t244 + t290;
t429 = qJD(1) * t356;
t419 = -m(4) - m(5) - m(6);
t415 = qJDD(2) * t310;
t413 = -rSges(4,3) + t501;
t412 = t94 * rSges(6,1) + t93 * rSges(6,2) + rSges(6,3) * t399;
t411 = t283 * t498;
t410 = t283 * t494;
t61 = t310 * t162 + t164 * t467 + t166 * t469;
t62 = -t310 * t163 - t165 * t467 - t167 * t469;
t406 = pkin(4) * t545 + pkin(7) * t399;
t405 = t278 + t433;
t404 = t289 + t432;
t403 = t283 * t509;
t394 = -pkin(4) - t498;
t393 = -t428 / 0.2e1;
t392 = t427 / 0.2e1;
t391 = -t425 / 0.2e1;
t390 = t425 / 0.2e1;
t389 = -t424 / 0.2e1;
t388 = t424 / 0.2e1;
t385 = qJD(3) * t308 - t291;
t383 = qJD(1) * t200 + t404;
t377 = rSges(4,1) * t304 + t496;
t174 = -t308 * rSges(4,3) + t310 * t377;
t381 = t174 + t387;
t379 = t201 + t386;
t378 = t266 - t385;
t249 = rSges(2,1) * t310 - rSges(2,2) * t308;
t246 = rSges(2,1) * t308 + rSges(2,2) * t310;
t123 = -pkin(7) * t401 + t406;
t51 = -rSges(6,3) * t401 + t412;
t11 = qJD(1) * t123 + qJDD(1) * t195 + t115 * t205 - t121 * t160 - t218 * t95 - t228 * t234 + t257 * t51 + (-qJD(4) * t216 - qJDD(2)) * t310 + t316;
t373 = t11 * t308 + t555;
t348 = t164 * t283 + t166 * t282;
t321 = qJD(1) * t348 + qJD(4) * t176 + t430;
t322 = -qJD(1) * t346 - qJD(4) * t470 + t431;
t372 = (t321 * t308 + t310 * t530) * t310 + (t322 * t308 - t310 * t533) * t308;
t371 = (-t308 * t530 + t321 * t310) * t310 + (t308 * t533 + t322 * t310) * t308;
t370 = t27 * t310 + t28 * t308;
t369 = t27 * t308 - t28 * t310;
t30 = t108 * t466 - t352;
t368 = t29 * t310 + t30 * t308;
t367 = t29 * t308 - t30 * t310;
t366 = t308 * t42 + t310 * t41;
t365 = t308 * t41 - t310 * t42;
t364 = t308 * t62 + t310 * t61;
t150 = t308 * t162;
t63 = -t348 * t310 + t150;
t64 = -t163 * t308 + t535;
t363 = t308 * t64 + t310 * t63;
t66 = -t226 * t424 + (t168 + t379) * qJD(1) + t385;
t362 = t308 * t65 - t310 * t66;
t361 = t405 + t437;
t354 = t104 * t310 - t105 * t308;
t350 = t115 * t310 + t117 * t308;
t345 = -t168 * t308 - t169 * t310;
t339 = -t306 * t310 + t247 + t269;
t138 = rSges(6,3) * t469 + (-t410 + t411) * t308;
t332 = -t272 + t507 + t508;
t331 = t376 + t508;
t325 = t106 * t218 + t108 * t217 + t154 * t257;
t324 = (Icges(6,5) * t207 - Icges(6,6) * t208) * t218 + (Icges(6,5) * t209 + Icges(6,6) * t210) * t217 + t175 * t257;
t318 = t344 * qJD(1) - t356 * qJD(4);
t34 = qJD(1) * t340 + t433 - t554;
t35 = -t228 * t424 + t115 * t257 - t160 * t218 + (t195 + t379) * qJD(1) + t385;
t38 = -t115 * t217 + t117 * t218 + (-t195 * t308 - t197 * t310) * qJD(4);
t313 = t38 * t350 + (-t308 * t35 - t310 * t34) * t160;
t232 = t310 * t410;
t198 = t228 * t310;
t187 = t226 * t308;
t159 = -t282 * t498 + t439;
t149 = qJD(1) * t204 - t291;
t148 = qJD(1) * t440 + t290;
t139 = t232 + (-t411 - t492) * t310;
t137 = t158 * t310;
t136 = t158 * t308;
t135 = t156 * t310;
t134 = t156 * t308;
t131 = rSges(6,1) * t209 + rSges(6,2) * t210;
t130 = rSges(6,1) * t207 - rSges(6,2) * t208;
t103 = qJD(1) * t534 + t385;
t102 = qJD(1) * t381 + t433;
t81 = t345 * qJD(4);
t77 = t308 * t344 + t470;
t76 = t77 * qJD(1);
t69 = qJD(1) * t434 + qJDD(1) * t248 + t408 - t415;
t68 = t440 * qJDD(1) + (-qJD(1) * t248 - t214) * qJD(1) + t436;
t56 = -t415 + qJDD(1) * t173 + qJD(1) * (-rSges(4,3) * t428 + t438) + t317;
t55 = t381 * qJDD(1) + (-qJD(1) * t173 - t214) * qJD(1) + t328;
t40 = t346 * qJD(4) + t100 * t283 - t282 * t98;
t39 = -qJD(4) * t348 + t101 * t283 - t282 * t99;
t37 = -t308 * t531 + t318 * t310;
t36 = t318 * t308 + t310 * t531;
t24 = qJD(4) * t363 - t460;
t23 = qJD(4) * t364 + t76;
t16 = t154 * t330 + t156 * t93 + t158 * t94 + t207 * t89 + t208 * t90 - t467 * t88;
t15 = t154 * t329 + t156 * t91 + t158 * t92 + t209 * t89 - t210 * t90 + t466 * t88;
t14 = t217 * t42 + t218 * t41 + t257 * t57;
t13 = t217 * t30 - t550;
t9 = -t115 * t120 + t117 * t121 - t195 * t233 - t197 * t234 - t217 * t51 + t218 * t50 + (t122 * t310 - t123 * t308) * qJD(4);
t6 = t108 * t330 + t111 * t93 - t113 * t94 + t207 * t46 + t208 * t48 - t44 * t467;
t5 = t106 * t330 + t109 * t93 + t112 * t94 + t207 * t47 + t208 * t49 - t45 * t467;
t4 = t108 * t329 + t111 * t91 - t113 * t92 + t209 * t46 - t210 * t48 + t44 * t466;
t3 = t106 * t329 + t109 * t91 + t112 * t92 + t209 * t47 - t210 * t49 + t45 * t466;
t2 = t120 * t28 + t121 * t27 + t16 * t257 + t205 * t52 + t217 * t6 + t218 * t5;
t1 = t120 * t30 + t121 * t29 + t15 * t257 + t205 * t53 + t217 * t4 + t218 * t3;
t17 = [(-(-t102 + (t174 - t474) * qJD(1) + t404) * t103 - t102 * t385 + t103 * (t405 + t438) + (t102 * t413 * t310 + (t102 * (-qJ(2) - t377) + t103 * t413) * t308) * qJD(1) + (-g(2) + t56) * t534 + (-g(1) + t55) * (t308 * t501 + t174 + t293)) * m(4) + (-t347 + t77) * t518 + t486 / 0.2e1 + (t39 + t37) * t388 - m(2) * (-g(1) * t246 + g(2) * t249) - t233 * t78 / 0.2e1 + t502 / 0.2e1 + (-(t199 - t65 + (t169 - t474) * qJD(1) + t383) * t66 + t65 * (rSges(5,1) * t396 - rSges(5,2) * t398 + t378) + t66 * (-rSges(5,2) * t399 + t361 + t407) + (t510 * t487 + (t65 * (-qJ(2) - t331) + t66 * t510) * t308) * qJD(1) + (-g(2) + t26) * (t339 + t168) + (-g(1) + t25) * (t293 - t297 + (-pkin(1) + t306) * t308 + t331 * t310)) * m(5) + (Icges(4,1) * t305 ^ 2 + (-0.2e1 * Icges(4,4) * t305 + Icges(4,2) * t304) * t304 + m(2) * (t246 ^ 2 + t249 ^ 2) - t343 + Icges(2,3) + Icges(3,1)) * qJDD(1) + (-(qJD(1) * t245 - t148 + t432) * t149 + t148 * t291 + t149 * (t434 + t435) + (t148 * t511 * t310 + (t148 * (-rSges(3,3) - qJ(2)) - t149 * pkin(1)) * t308) * qJD(1) + (-g(2) + t69) * t204 + (-g(1) + t68) * (t308 * t511 + t293 + t493)) * m(3) + (-qJD(4) * t344 + t212 * t282 - t213 * t283) * qJD(1) + t485 / 0.2e1 + t80 * t519 + t16 * t520 + t52 * t525 + t53 * t526 + t500 + (t24 + t460 + (t163 * t308 ^ 2 + (-t150 + t62 + (t163 + t348) * t310) * t310) * qJD(4)) * t389 + t503 / 0.2e1 + (t15 + t12) * t522 + (t76 + ((-t63 + t150 + t62) * t308 + (t64 - t535 + (t163 - t348) * t308 + t61) * t310) * qJD(4)) * t391 + ((t30 + t551) * t218 + t540) * t523 + ((-t27 + t551) * t217 + t13 + t550) * t521 + ((t11 - g(2)) * (-t308 * t403 + t259 + t339 + t444) + (t10 - g(1)) * (t306 * t308 - t244 - t261 + t443) + (-t506 + t555) * t332 + (t378 - t382 + (t504 + t539) * t424 + (-t556 + (-qJ(2) - t332 + t273) * t308) * qJD(1)) * t34 + (t34 - t383 + t361 + t406 + t412 + (-t403 * t310 - t197 + t474 - t557) * qJD(1) + t554) * t35) * m(6) + (t36 + t23 + t40) * t390; (-m(3) + t419) * (g(1) * t308 - t505) + 0.2e1 * (t10 * t514 + t11 * t513) * m(6) + 0.2e1 * (t490 / 0.2e1 - t488 / 0.2e1) * m(5) + 0.2e1 * (t513 * t56 + t514 * t55) * m(4) + 0.2e1 * (t513 * t69 + t514 * t68) * m(3); t419 * (g(2) * t308 + t506) + m(4) * (t308 * t56 + t310 * t55) + m(5) * (t25 * t310 + t26 * t308) + m(6) * t373; qJDD(1) * (t308 * t80 - t310 * t347) / 0.2e1 + qJD(1) * (t308 * t40 + t310 * t39 + (t308 * t347 + t80 * t310) * qJD(1)) / 0.2e1 + (-g(1) * t187 + g(2) * t189 + g(3) * t376 + (qJD(4) * t354 - t168 * t233 - t169 * t234) * t345 + t81 * ((-t168 * t310 + t169 * t308) * qJD(1) + t354) - t362 * t215 + (t490 - t488 + (t308 * t66 + t487) * qJD(1)) * t226 - (t187 * t66 + t189 * t65) * qJD(1) - (t81 * (-t187 * t308 - t189 * t310) - t362 * t376) * qJD(4)) * m(5) + ((t176 * t424 - t429) * t310 + (-t543 + (-t310 * t470 - t544) * qJD(4)) * t308) * t389 + ((-t425 * t470 - t429) * t308 + (t543 + (t308 * t176 + t544) * qJD(4)) * t310) * t391 + (t23 + t12) * t393 - t14 * t422 / 0.2e1 + t24 * t392 + ((-t63 * t308 + t64 * t310) * qJD(1) + t372) * t390 + ((-t61 * t308 + t62 * t310) * qJD(1) + t371) * t388 + (t420 * t515 + t392) * t13 + (qJD(1) * t37 + qJD(4) * t371 + qJDD(1) * t77 + t233 * t62 + t234 * t61 + t2) * t310 / 0.2e1 + (-qJD(1) * t365 + t308 * t8 + t310 * t7) * t516 + t364 * t518 + t363 * t519 + (-qJD(1) * t369 + t308 * t6 + t310 * t5) * t520 + (-qJD(1) * t367 + t3 * t310 + t308 * t4) * t522 + t366 * t524 + t370 * t525 + t368 * t526 + (((-t134 * t307 + t136 * t309 + t106) * t218 + (t135 * t307 - t137 * t309 + t108) * t217 + (-t155 * t307 + t157 * t309 + t154) * t257 + t57 * qJD(5)) * t283 + (qJD(5) * t365 + t312) * t282) * t517 + ((t134 * t207 + t136 * t208) * t218 + (-t135 * t207 - t137 * t208) * t217 + (t155 * t207 + t157 * t208) * t257 + (-t28 * t468 + t283 * t52) * qJD(5) + ((qJD(5) * t27 + t325) * t282 - t538) * t308) * t521 + ((t134 * t209 - t136 * t210) * t218 + (-t135 * t209 + t137 * t210) * t217 + (t155 * t209 - t157 * t210) * t257 + (t283 * t53 + t29 * t469) * qJD(5) + ((-qJD(5) * t30 - t325) * t282 + t538) * t310) * t523 + (-t34 * (qJD(1) * t198 - t139 * t257 + t159 * t217 + t227 * t425) - t35 * (qJD(1) * t196 + t138 * t257 - t159 * t218 - t227 * t424) - t38 * (-t138 * t217 + t139 * t218 - t196 * t425 - t198 * t424) - ((t115 * t35 - t117 * t34) * t283 + t313 * t282) * qJD(5) + (-t11 * t450 - t35 * t484 + t9 * t453 + t38 * (t122 + t50) + (t34 * t450 + t38 * t454) * qJD(1)) * t310 + (t10 * t450 + t34 * t484 + t9 * t454 + t38 * (-t123 - t51) + (t35 * t450 - t38 * t453) * qJD(1)) * t308 - g(1) * (t138 + t196) - g(2) * t232 - g(3) * (t282 * t394 + t273 + t439) - (t283 * t394 - t539) * t505) * m(6) - qJD(1) * ((-t282 * t442 - t283 * t441) * qJD(1) + (t282 * t319 + t283 * t320) * qJD(4)) / 0.2e1 + t282 * t421 * t527 + (qJD(1) * t36 + qJD(4) * t372 - qJDD(1) * t78 + t233 * t64 + t234 * t63 + t1) * t514; t401 * t527 + t282 * t12 * t390 - t2 * t467 / 0.2e1 + (t282 * t52 - t283 * t369) * t525 + ((qJD(4) * t369 + t16) * t282 + (-qJD(1) * t370 + qJD(4) * t52 - t308 * t5 + t310 * t6) * t283) * t520 + t1 * t466 / 0.2e1 + (t282 * t53 - t283 * t367) * t526 + ((qJD(4) * t367 + t15) * t282 + (-qJD(1) * t368 + qJD(4) * t53 - t3 * t308 + t310 * t4) * t283) * t522 + t14 * t426 / 0.2e1 + (t485 + t486 + t500 + t502 + t503) * t515 + (t282 * t57 - t283 * t365) * t524 + ((qJD(4) * t365 + t18) * t282 + (-qJD(1) * t366 + qJD(4) * t57 - t308 * t7 + t310 * t8) * t283) * t516 + (t207 * t314 - t208 * t528 - t324 * t467) * t521 + (t314 * t209 + t210 * t528 + t324 * t466) * t523 + (t324 * t282 + (-t307 * t314 - t309 * t528) * t283) * t517 + (t282 * t389 + t283 * t393) * t13 + ((qJD(4) * t313 - t10 * t117 + t11 * t115 - t34 * t50 + t35 * t51) * t282 + (t34 * (-qJD(4) * t117 + t310 * t95) + t35 * (qJD(4) * t115 + t308 * t95) - t9 * t350 + t38 * (t115 * t428 - t117 * t427 - t308 * t50 - t310 * t51) + ((-t308 * t34 + t310 * t35) * qJD(1) + t373) * t160) * t283 - t34 * (-t131 * t257 + t188 * t217) - t35 * (t130 * t257 - t188 * t218) - t38 * (-t130 * t217 + t131 * t218) - g(1) * t130 - g(2) * t131 - g(3) * t188) * m(6);];
tau = t17;
