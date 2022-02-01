% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m [6x1]
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:35
% EndTime: 2022-01-20 09:48:46
% DurationCPUTime: 4.56s
% Computational Cost: add. (42067->334), mult. (23472->445), div. (0->0), fcn. (21350->10), ass. (0->220)
t314 = qJ(4) + qJ(5);
t311 = cos(t314);
t303 = Icges(6,4) * t311;
t310 = sin(t314);
t264 = -Icges(6,2) * t310 + t303;
t265 = Icges(6,1) * t310 + t303;
t482 = t264 + t265;
t317 = cos(qJ(4));
t312 = Icges(5,4) * t317;
t315 = sin(qJ(4));
t282 = -Icges(5,2) * t315 + t312;
t283 = Icges(5,1) * t315 + t312;
t481 = t282 + t283;
t313 = qJ(1) + pkin(9);
t309 = qJ(3) + t313;
t305 = cos(t309);
t300 = t305 * pkin(7);
t304 = sin(t309);
t425 = rSges(5,1) * t317;
t358 = pkin(3) + t425;
t396 = t304 * t315;
t361 = rSges(5,2) * t396 + t305 * rSges(5,3);
t182 = -t304 * t358 + t300 + t361;
t391 = t305 * t315;
t279 = rSges(5,2) * t391;
t183 = -t279 + t358 * t305 + (rSges(5,3) + pkin(7)) * t304;
t285 = rSges(5,1) * t315 + rSges(5,2) * t317;
t240 = t285 * t304;
t241 = t285 * t305;
t103 = t182 * t240 - t183 * t241;
t397 = t304 * t311;
t398 = t304 * t310;
t208 = rSges(6,1) * t397 - rSges(6,2) * t398 - t305 * rSges(6,3);
t429 = pkin(4) * t317;
t306 = pkin(3) + t429;
t463 = -pkin(8) - pkin(7);
t362 = -t304 * t306 - t305 * t463;
t177 = -t208 + t362;
t289 = t304 * t463;
t393 = t305 * t310;
t355 = -rSges(6,2) * t393 + t304 * rSges(6,3);
t424 = rSges(6,1) * t311;
t178 = -t289 + (t306 + t424) * t305 + t355;
t267 = rSges(6,1) * t310 + rSges(6,2) * t311;
t430 = pkin(4) * t315;
t333 = t267 + t430;
t473 = t333 * t305;
t474 = t333 * t304;
t91 = t177 * t474 - t178 * t473;
t480 = m(5) * t103 + m(6) * t91;
t465 = m(5) / 0.2e1;
t464 = m(6) / 0.2e1;
t438 = t304 / 0.2e1;
t437 = -t305 / 0.2e1;
t479 = t305 / 0.2e1;
t347 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t313);
t348 = -sin(qJ(1)) * pkin(1) - pkin(2) * sin(t313);
t435 = m(4) * (t347 * (-rSges(4,1) * t304 - rSges(4,2) * t305) - (t305 * rSges(4,1) - t304 * rSges(4,2)) * t348);
t301 = t304 ^ 2;
t302 = t305 ^ 2;
t360 = t301 + t302;
t439 = -t304 / 0.2e1;
t476 = t438 + t439;
t166 = t177 + t348;
t167 = t178 + t347;
t88 = t166 * t474 - t167 * t473;
t232 = t267 * t304;
t233 = t267 * t305;
t168 = -t304 * t232 - t305 * t233;
t268 = -rSges(6,2) * t310 + t424;
t342 = Icges(6,5) * t310 + Icges(6,6) * t311;
t226 = t342 * t304;
t227 = t305 * t342;
t414 = Icges(6,4) * t310;
t266 = Icges(6,1) * t311 - t414;
t205 = Icges(6,5) * t304 + t266 * t305;
t263 = Icges(6,2) * t311 + t414;
t367 = -t263 * t305 + t205;
t203 = Icges(6,6) * t304 + t264 * t305;
t369 = -t265 * t305 - t203;
t328 = -t310 * t367 + t311 * t369;
t271 = Icges(6,4) * t398;
t204 = Icges(6,1) * t397 - Icges(6,5) * t305 - t271;
t368 = -Icges(6,2) * t397 + t204 - t271;
t202 = Icges(6,4) * t397 - Icges(6,2) * t398 - Icges(6,6) * t305;
t370 = t265 * t304 + t202;
t329 = t310 * t368 + t311 * t370;
t428 = (-t301 * t227 + (t329 * t305 + (t226 + t328) * t304) * t305) * t438 + (-t302 * t226 + (t328 * t304 + (t227 + t329) * t305) * t304) * t437;
t392 = t305 * t311;
t148 = t304 * t208 + t305 * (rSges(6,1) * t392 + t355);
t94 = -t304 * (pkin(3) * t304 - t300 + t362) + (-t304 * pkin(7) - t289 + (-pkin(3) + t306) * t305) * t305 + t148;
t15 = t428 + m(6) * (t94 * t168 + (t304 * t474 + t305 * t473) * t268);
t475 = t15 * qJD(5);
t84 = -t178 * t166 + t167 * t177;
t175 = t182 + t348;
t176 = t183 + t347;
t89 = -t183 * t175 + t176 * t182;
t472 = qJD(1) + qJD(3);
t427 = (t88 - t91) * t464 + ((-t176 + t183) * t305 + (t175 - t182) * t304) * t285 * t465;
t97 = t175 * t240 - t176 * t241;
t471 = (t91 + t88) * t464 + (t103 + t97) * t465;
t415 = Icges(5,4) * t315;
t281 = Icges(5,2) * t317 + t415;
t284 = Icges(5,1) * t317 - t415;
t470 = t481 * t317 / 0.2e1 + (t284 / 0.2e1 - t281 / 0.2e1) * t315;
t350 = t482 * t311 / 0.2e1 + (-t263 / 0.2e1 + t266 / 0.2e1) * t310;
t179 = t205 * t397;
t262 = Icges(6,5) * t311 - Icges(6,6) * t310;
t403 = t262 * t305;
t201 = Icges(6,3) * t304 + t403;
t354 = t201 * t305 - t179;
t105 = -t203 * t398 - t354;
t200 = Icges(6,5) * t397 - Icges(6,6) * t398 - Icges(6,3) * t305;
t375 = -t304 * t200 - t204 * t392;
t106 = -t202 * t393 - t375;
t374 = t304 * t201 + t205 * t392;
t107 = -t203 * t393 + t374;
t352 = t203 * t310 - t200;
t407 = t202 * t310;
t356 = ((t105 - t179 + (t201 + t407) * t305 + t375) * t305 + t374 * t304) * t437 + (-t106 * t305 + t107 * t304) * t479 + ((t304 * t352 + t105 + t106 + t354) * t304 + (t107 - t374 + t304 * (-t204 * t311 + t407) + (t352 + t200) * t305) * t305) * t438;
t469 = 4 * qJD(1);
t467 = 2 * qJD(4);
t460 = m(5) * t89;
t458 = m(5) * t97;
t90 = t166 * t232 - t167 * t233;
t96 = t177 * t232 - t178 * t233;
t452 = m(6) * (t96 + t90);
t450 = m(6) * ((-t167 + t178) * t305 + (t166 - t177) * t304) * t267;
t337 = (-t166 * t305 - t167 * t304) * t268;
t376 = -t232 * t473 + t233 * t474;
t449 = m(6) * (t337 + t376);
t335 = (t304 * t473 - t305 * t474) * t267;
t448 = m(6) * (t337 + t335);
t336 = (-t177 * t305 - t178 * t304) * t268;
t447 = m(6) * (t336 + t376);
t446 = m(6) * (t336 + t335);
t445 = m(6) * t84;
t443 = m(6) * t88;
t442 = m(6) * t90;
t440 = m(6) * t96;
t395 = t304 * t317;
t219 = Icges(5,4) * t395 - Icges(5,2) * t396 - Icges(5,6) * t305;
t217 = Icges(5,5) * t395 - Icges(5,6) * t396 - Icges(5,3) * t305;
t277 = Icges(5,4) * t396;
t221 = Icges(5,1) * t395 - Icges(5,5) * t305 - t277;
t390 = t305 * t317;
t372 = -t304 * t217 - t221 * t390;
t118 = -t219 * t391 - t372;
t220 = Icges(5,6) * t304 + t282 * t305;
t280 = Icges(5,5) * t317 - Icges(5,6) * t315;
t401 = t280 * t305;
t218 = Icges(5,3) * t304 + t401;
t222 = Icges(5,5) * t304 + t284 * t305;
t371 = t304 * t218 + t222 * t390;
t119 = -t220 * t391 + t371;
t351 = t220 * t315 - t217;
t192 = t222 * t395;
t353 = t218 * t305 - t192;
t28 = (t305 * t351 + t119 - t371) * t305 + (t304 * t351 + t118 + t353) * t304;
t117 = -t220 * t396 - t353;
t405 = t219 * t315;
t29 = (t117 - t192 + (t218 + t405) * t305 + t372) * t305 + t371 * t304;
t81 = -(-t304 * (-t221 * t317 + t405) - t217 * t305) * t305 + t117 * t304;
t82 = -t118 * t305 + t119 * t304;
t2 = (t82 / 0.2e1 - t29 / 0.2e1) * t305 + (t28 / 0.2e1 + t81 / 0.2e1) * t304 + t356;
t436 = t2 * qJD(4);
t431 = m(6) * t168;
t419 = qJD(5) * t356;
t366 = t283 * t304 + t219;
t365 = -t283 * t305 - t220;
t364 = -Icges(5,2) * t395 + t221 - t277;
t363 = -t281 * t305 + t222;
t357 = -t268 - t429;
t346 = t452 / 0.2e1 + t350;
t343 = Icges(5,5) * t315 + Icges(5,6) * t317;
t334 = (-t240 * t305 + t241 * t304) * t285;
t325 = (-t263 + t266) * t311 - t482 * t310;
t332 = -t356 + (t262 * t304 + t305 * t325 + t310 * t369 + t311 * t367) * t438 + (t304 * t325 - t310 * t370 + t311 * t368 - t403) * t437;
t331 = -t350 + t476 * (t311 * t202 + t310 * t204);
t330 = t350 + t470;
t327 = t315 * t364 + t317 * t366;
t326 = -t315 * t363 + t317 * t365;
t324 = (-t281 + t284) * t317 - t481 * t315;
t322 = t330 + t471;
t321 = (-t232 * t305 + t233 * t304) * t267;
t320 = t331 - t470 + t476 * (t317 * t219 + t315 * t221);
t319 = (t29 * t479 + t332 + (t304 * t324 - t315 * t366 + t317 * t364 - t401 + t82) * t437 + (t28 + t81) * t439 + (t280 * t304 + t305 * t324 + t315 * t365 + t317 * t363) * t438) * qJD(4);
t287 = -rSges(5,2) * t315 + t425;
t235 = t305 * t343;
t234 = t343 * t304;
t216 = t357 * t305;
t214 = t357 * t304;
t171 = -t240 * t304 - t241 * t305;
t159 = qJD(5) * t431;
t145 = -t360 * t430 + t168;
t83 = t350 + t440;
t70 = t446 / 0.2e1;
t69 = t350 + t442;
t68 = t447 / 0.2e1;
t58 = t448 / 0.2e1;
t56 = t449 / 0.2e1;
t53 = t450 / 0.2e1;
t39 = t330 + t480;
t34 = t330 + t443 + t458;
t33 = t435 + t445 + t460;
t19 = -t450 / 0.2e1 + t346;
t18 = t53 + t346;
t17 = m(6) * (t360 * t267 * t268 + t148 * t168) + t428;
t16 = t17 * qJD(5);
t14 = t53 - t452 / 0.2e1 + t331;
t13 = t322 - t427;
t12 = t322 + t427;
t9 = t320 + t427 - t471;
t8 = t68 - t446 / 0.2e1 + t356;
t7 = t70 - t447 / 0.2e1 + t356;
t6 = t56 - t448 / 0.2e1 + t356;
t5 = t58 - t449 / 0.2e1 + t356;
t4 = t68 + t70 + t332;
t3 = t56 + t58 + t332;
t1 = [qJD(3) * t33 + qJD(4) * t34 + qJD(5) * t69, 0, t33 * qJD(1) + t12 * qJD(4) + t18 * qJD(5) + 0.2e1 * (t435 / 0.2e1 + t84 * t464 + t89 * t465) * qJD(3), t34 * qJD(1) + t12 * qJD(3) + t3 * qJD(5) + ((t166 * t216 + t167 * t214) * t464 + ((-t175 * t305 - t176 * t304) * t287 + t334) * t465) * t467 + t319, t69 * qJD(1) + t18 * qJD(3) + t3 * qJD(4) + ((t337 + t321) * m(6) + t332) * qJD(5); 0, 0, 0, (t145 * t464 + t171 * t465) * t467 + t159, qJD(4) * t431 + t159; t13 * qJD(4) + t19 * qJD(5) + (-t445 / 0.4e1 - t460 / 0.4e1 - t435 / 0.4e1) * t469, 0, qJD(4) * t39 + qJD(5) * t83, t13 * qJD(1) + t39 * qJD(3) + t4 * qJD(5) + ((t177 * t216 + t178 * t214) * t464 + ((-t182 * t305 - t183 * t304) * t287 + t334) * t465) * t467 + t319, t19 * qJD(1) + t83 * qJD(3) + t4 * qJD(4) + ((t336 + t321) * m(6) + t332) * qJD(5); t320 * qJD(1) + t9 * qJD(3) + t6 * qJD(5) + (-t443 / 0.4e1 - t458 / 0.4e1) * t469 + t436, 0, t9 * qJD(1) + t8 * qJD(5) + t436 + (t320 - t480) * qJD(3), (m(5) * (t285 * t287 * t360 + (t304 * (rSges(5,1) * t395 - t361) + t305 * (rSges(5,1) * t390 + t304 * rSges(5,3) - t279)) * t171) + (-t301 * t235 + (t327 * t305 + (t234 + t326) * t304) * t305) * t438 + (-t302 * t234 + (t326 * t304 + (t235 + t327) * t305) * t304) * t437 + m(6) * (t145 * t94 - t214 * t474 - t216 * t473) + t428) * qJD(4) + t475 + t472 * t2, t6 * qJD(1) + t8 * qJD(3) + t15 * qJD(4) + t475; (t331 - t442) * qJD(1) + t14 * qJD(3) + t5 * qJD(4) + t419, 0, t14 * qJD(1) + (t331 - t440) * qJD(3) + t7 * qJD(4) + t419, t5 * qJD(1) + t7 * qJD(3) + ((t145 * t148 + (-t214 * t304 - t216 * t305) * t267) * m(6) + t428) * qJD(4) + t16, qJD(4) * t17 + t356 * t472 + t16;];
Cq = t1;
