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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:53:49
% EndTime: 2020-01-03 11:53:58
% DurationCPUTime: 4.80s
% Computational Cost: add. (42067->346), mult. (23472->450), div. (0->0), fcn. (21350->10), ass. (0->233)
t310 = qJ(4) + qJ(5);
t305 = sin(t310);
t306 = cos(t310);
t267 = rSges(6,1) * t305 + rSges(6,2) * t306;
t309 = qJ(1) + pkin(9);
t304 = qJ(3) + t309;
t300 = cos(t304);
t233 = t267 * t300;
t299 = sin(t304);
t311 = sin(qJ(4));
t425 = pkin(4) * t311;
t342 = t267 + t425;
t478 = t342 * t299;
t176 = t478 * t233;
t409 = Icges(6,4) * t305;
t266 = Icges(6,1) * t306 - t409;
t205 = -Icges(6,5) * t300 + t266 * t299;
t385 = t299 * t306;
t184 = t205 * t385;
t379 = t300 * t306;
t380 = t300 * t305;
t202 = Icges(6,5) * t379 - Icges(6,6) * t380 + Icges(6,3) * t299;
t271 = Icges(6,4) * t380;
t206 = Icges(6,1) * t379 + Icges(6,5) * t299 - t271;
t204 = Icges(6,4) * t379 - Icges(6,2) * t380 + Icges(6,6) * t299;
t399 = t204 * t305;
t327 = t206 * t306 - t399;
t489 = -t202 * t299 - t327 * t300 - t184;
t313 = cos(qJ(4));
t424 = pkin(4) * t313;
t301 = pkin(3) + t424;
t386 = t299 * t305;
t350 = -rSges(6,2) * t386 - t300 * rSges(6,3);
t315 = -pkin(8) - pkin(7);
t376 = t300 * t315;
t419 = rSges(6,1) * t306;
t180 = t376 + (t301 + t419) * t299 + t350;
t335 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t309);
t168 = t180 + t335;
t488 = -t168 + t180;
t272 = t300 * t301;
t336 = rSges(6,1) * t379 - rSges(6,2) * t380;
t181 = t272 + (rSges(6,3) - t315) * t299 + t336;
t344 = pkin(2) * cos(t309) + cos(qJ(1)) * pkin(1);
t169 = t181 + t344;
t487 = -t169 + t181;
t262 = Icges(6,5) * t306 - Icges(6,6) * t305;
t392 = t262 * t299;
t201 = -Icges(6,3) * t300 + t392;
t486 = t201 * t299 + t205 * t379;
t298 = Icges(6,4) * t306;
t264 = -Icges(6,2) * t305 + t298;
t474 = Icges(6,1) * t305 + t298;
t485 = t264 + t474;
t307 = Icges(5,4) * t313;
t282 = -Icges(5,2) * t311 + t307;
t473 = Icges(5,1) * t311 + t307;
t484 = t282 + t473;
t293 = t300 * pkin(7);
t384 = t299 * t311;
t349 = -rSges(5,2) * t384 - t300 * rSges(5,3);
t420 = rSges(5,1) * t313;
t187 = -t293 + (pkin(3) + t420) * t299 + t349;
t377 = t300 * t313;
t378 = t300 * t311;
t324 = rSges(5,1) * t377 - rSges(5,2) * t378 + rSges(5,3) * t299;
t346 = -t300 * pkin(3) - t299 * pkin(7);
t188 = t324 - t346;
t285 = rSges(5,1) * t311 + rSges(5,2) * t313;
t240 = t285 * t299;
t241 = t285 * t300;
t103 = -t187 * t240 - t188 * t241;
t477 = t342 * t300;
t91 = -t180 * t478 - t181 * t477;
t483 = m(5) * t103 + m(6) * t91;
t459 = m(5) / 0.2e1;
t458 = m(6) / 0.2e1;
t434 = -t299 / 0.2e1;
t482 = t299 / 0.2e1;
t433 = -t300 / 0.2e1;
t430 = m(4) * (-(rSges(4,1) * t299 + rSges(4,2) * t300) * t344 + (t300 * rSges(4,1) - rSges(4,2) * t299) * t335);
t287 = -rSges(5,2) * t311 + t420;
t479 = t287 * t459;
t295 = t299 ^ 2;
t296 = t300 ^ 2;
t345 = t295 + t296;
t432 = t300 / 0.2e1;
t476 = t432 + t433;
t232 = t267 * t299;
t170 = -t299 * t232 - t233 * t300;
t268 = -rSges(6,2) * t305 + t419;
t389 = t268 * t300;
t390 = t268 * t299;
t329 = Icges(6,5) * t305 + Icges(6,6) * t306;
t226 = t299 * t329;
t227 = t329 * t300;
t263 = Icges(6,2) * t306 + t409;
t358 = -t263 * t299 + t205;
t203 = -Icges(6,6) * t300 + t264 * t299;
t360 = t299 * t474 + t203;
t468 = -t358 * t305 - t360 * t306;
t357 = -Icges(6,2) * t379 + t206 - t271;
t359 = t300 * t474 + t204;
t469 = -t357 * t305 - t359 * t306;
t423 = (-t296 * t226 + (t469 * t299 + (t227 - t468) * t300) * t299) * t433 + (t295 * t227 + (t468 * t300 + (-t226 - t469) * t299) * t300) * t434;
t191 = t299 * (rSges(6,1) * t385 + t350);
t209 = t299 * rSges(6,3) + t336;
t95 = t191 + (-t299 * t315 + t209 + t272 + t346) * t300 + (t376 + t293 + (-pkin(3) + t301) * t299) * t299;
t15 = t423 + m(6) * (t95 * t170 + t389 * t477 + t390 * t478);
t475 = t15 * qJD(5);
t84 = t168 * t181 - t180 * t169;
t178 = t187 + t335;
t179 = t188 + t344;
t89 = t178 * t188 - t187 * t179;
t472 = qJD(1) + qJD(3);
t422 = (t487 * t477 + t488 * t478) * t458 + ((-t179 + t188) * t300 + (-t178 + t187) * t299) * t285 * t459;
t88 = -t168 * t478 - t169 * t477;
t97 = -t178 * t240 - t179 * t241;
t471 = (t91 + t88) * t458 + (t103 + t97) * t459;
t470 = t485 * t305 + (t263 - t266) * t306;
t410 = Icges(5,4) * t311;
t281 = Icges(5,2) * t313 + t410;
t284 = Icges(5,1) * t313 - t410;
t467 = t484 * t311 + (t281 - t284) * t313;
t277 = Icges(5,4) * t378;
t222 = Icges(5,1) * t377 + Icges(5,5) * t299 - t277;
t353 = -Icges(5,2) * t377 + t222 - t277;
t220 = Icges(5,4) * t377 - Icges(5,2) * t378 + Icges(5,6) * t299;
t355 = t300 * t473 + t220;
t466 = -t311 * t353 - t313 * t355;
t221 = -Icges(5,5) * t300 + t284 * t299;
t354 = -t281 * t299 + t221;
t219 = -Icges(5,6) * t300 + t282 * t299;
t356 = t299 * t473 + t219;
t465 = -t311 * t354 - t313 * t356;
t464 = t484 * t313 / 0.2e1 + (t284 / 0.2e1 - t281 / 0.2e1) * t311;
t338 = t485 * t306 / 0.2e1 + (-t263 / 0.2e1 + t266 / 0.2e1) * t305;
t104 = -t201 * t300 - t203 * t386 + t184;
t185 = t206 * t385;
t105 = t202 * t300 + t204 * t386 - t185;
t398 = t205 * t306;
t400 = t203 * t305;
t340 = ((t185 + (t201 - t399) * t299 - t486) * t299 + ((t201 + t327) * t300 + (t398 + t400) * t299 + t489) * t300) * t434 + (-t104 * t300 - t105 * t299) * t482 + ((-t105 + (t202 - t398) * t300 + t486) * t300 + (t104 + (t202 + t400) * t299 + t489) * t299) * t433;
t463 = 4 * qJD(1);
t461 = 2 * qJD(4);
t455 = m(5) * t89;
t453 = m(5) * t97;
t90 = -t168 * t232 - t169 * t233;
t96 = -t180 * t232 - t181 * t233;
t447 = m(6) * (t96 + t90);
t445 = m(6) * (t488 * t299 + t487 * t300) * t267;
t146 = t168 * t389;
t362 = -t232 * t477 + t176;
t444 = m(6) * (-t169 * t390 + t146 + t362);
t401 = t477 * t267;
t405 = t169 * t268;
t443 = m(6) * (t146 - t176 + (t401 - t405) * t299);
t158 = t180 * t389;
t442 = m(6) * (-t181 * t390 + t158 + t362);
t403 = t181 * t268;
t441 = m(6) * (t158 - t176 + (t401 - t403) * t299);
t440 = m(6) * t84;
t438 = m(6) * t88;
t437 = m(6) * t90;
t435 = m(6) * t96;
t196 = t219 * t378;
t280 = Icges(5,5) * t313 - Icges(5,6) * t311;
t388 = t280 * t299;
t217 = -Icges(5,3) * t300 + t388;
t118 = -t217 * t299 - t221 * t377 + t196;
t218 = Icges(5,5) * t377 - Icges(5,6) * t378 + Icges(5,3) * t299;
t396 = t220 * t311;
t325 = t222 * t313 - t396;
t119 = t218 * t299 + t325 * t300;
t383 = t299 * t313;
t194 = t221 * t383;
t195 = t222 * t383;
t395 = t221 * t313;
t397 = t219 * t311;
t28 = (t118 + t195 - t196 + (t217 - t396) * t299) * t299 + (-t194 - t119 + (t217 + t325) * t300 + (t395 + t397) * t299) * t300;
t116 = -t217 * t300 - t219 * t384 + t194;
t117 = t218 * t300 + t220 * t384 - t195;
t29 = (-t117 + t196 + (t218 - t395) * t300) * t300 + (t116 - t194 + (t218 + t397) * t299) * t299;
t81 = -t116 * t300 - t117 * t299;
t82 = -t118 * t300 - t119 * t299;
t2 = (-t29 / 0.2e1 - t82 / 0.2e1) * t300 + (t81 / 0.2e1 - t28 / 0.2e1) * t299 + t340;
t431 = t2 * qJD(4);
t426 = m(6) * t170;
t414 = qJD(5) * t340;
t394 = t233 * t267;
t341 = t268 + t424;
t339 = t232 * t233;
t333 = t447 / 0.2e1 + t338;
t330 = Icges(5,5) * t311 + Icges(5,6) * t313;
t322 = -t340 + (t470 * t300 + t359 * t305 - t357 * t306 - t392) * t434 + (-t262 * t300 - t470 * t299 - t360 * t305 + t358 * t306) * t433;
t321 = -t338 + t476 * (t306 * t204 + t305 * t206);
t320 = t338 + t464;
t318 = t320 + t471;
t317 = t321 - t464 + t476 * (t313 * t220 + t311 * t222);
t316 = (t322 + t28 * t482 + (t467 * t300 + t355 * t311 - t353 * t313 - t388 + t81) * t434 + (-t280 * t300 - t467 * t299 - t356 * t311 + t354 * t313) * t433 + (t82 + t29) * t432) * qJD(4);
t235 = t330 * t300;
t234 = t299 * t330;
t216 = t341 * t300;
t214 = t341 * t299;
t173 = -t240 * t299 - t241 * t300;
t161 = qJD(5) * t426;
t150 = t300 * t209 + t191;
t145 = -t345 * t425 + t170;
t83 = t338 + t435;
t70 = t441 / 0.2e1;
t69 = t338 + t437;
t68 = t442 / 0.2e1;
t58 = t443 / 0.2e1;
t56 = t444 / 0.2e1;
t53 = t445 / 0.2e1;
t39 = t320 + t483;
t34 = t320 + t438 + t453;
t33 = t430 + t440 + t455;
t19 = -t445 / 0.2e1 + t333;
t18 = t53 + t333;
t17 = m(6) * (t267 * t268 * t345 + t150 * t170) + t423;
t16 = t17 * qJD(5);
t14 = t53 - t447 / 0.2e1 + t321;
t13 = t318 + t422;
t12 = t318 - t422;
t9 = t317 + t422 - t471;
t8 = t68 - t441 / 0.2e1 + t340;
t7 = t70 - t442 / 0.2e1 + t340;
t6 = t56 - t443 / 0.2e1 + t340;
t5 = t58 - t444 / 0.2e1 + t340;
t4 = t68 + t70 + t322;
t3 = t56 + t58 + t322;
t1 = [qJD(3) * t33 + qJD(4) * t34 + qJD(5) * t69, 0, t33 * qJD(1) + t13 * qJD(4) + t18 * qJD(5) + 0.2e1 * (t430 / 0.2e1 + t84 * t458 + t89 * t459) * qJD(3), t34 * qJD(1) + t13 * qJD(3) + t3 * qJD(5) + ((t178 * t300 - t179 * t299) * t479 + (t168 * t216 - t169 * t214) * t458) * t461 + t316, t69 * qJD(1) + t18 * qJD(3) + t3 * qJD(4) + ((t146 + (t394 - t405) * t299 - t339) * m(6) + t322) * qJD(5); 0, 0, 0, (t145 * t458 + t173 * t459) * t461 + t161, qJD(4) * t426 + t161; t12 * qJD(4) + t19 * qJD(5) + (-t430 / 0.4e1 - t455 / 0.4e1 - t440 / 0.4e1) * t463, 0, qJD(4) * t39 + qJD(5) * t83, t12 * qJD(1) + t39 * qJD(3) + t4 * qJD(5) + ((t180 * t216 - t181 * t214) * t458 + (t187 * t300 - t188 * t299) * t479) * t461 + t316, t19 * qJD(1) + t83 * qJD(3) + t4 * qJD(4) + ((t158 + (t394 - t403) * t299 - t339) * m(6) + t322) * qJD(5); t317 * qJD(1) + t9 * qJD(3) + t6 * qJD(5) + (-t453 / 0.4e1 - t438 / 0.4e1) * t463 + t431, 0, t9 * qJD(1) + t8 * qJD(5) + t431 + (t317 - t483) * qJD(3), (m(5) * (t285 * t287 * t345 + (t300 * t324 + t299 * (rSges(5,1) * t383 + t349)) * t173) + (-t296 * t234 + (t466 * t299 + (t235 - t465) * t300) * t299) * t433 + (t295 * t235 + (t465 * t300 + (-t234 - t466) * t299) * t300) * t434 + m(6) * (t145 * t95 + t214 * t478 + t216 * t477) + t423) * qJD(4) + t475 + t472 * t2, t6 * qJD(1) + t8 * qJD(3) + t15 * qJD(4) + t475; (t321 - t437) * qJD(1) + t14 * qJD(3) + t5 * qJD(4) + t414, 0, t14 * qJD(1) + (t321 - t435) * qJD(3) + t7 * qJD(4) + t414, t5 * qJD(1) + t7 * qJD(3) + ((t145 * t150 + (t214 * t299 + t216 * t300) * t267) * m(6) + t423) * qJD(4) + t16, qJD(4) * t17 + t340 * t472 + t16;];
Cq = t1;
