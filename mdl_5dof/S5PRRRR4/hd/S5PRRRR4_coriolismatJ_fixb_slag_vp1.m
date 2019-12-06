% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:41
% EndTime: 2019-12-05 17:07:50
% DurationCPUTime: 4.27s
% Computational Cost: add. (41881->332), mult. (23278->442), div. (0->0), fcn. (21162->8), ass. (0->221)
t312 = qJ(4) + qJ(5);
t309 = cos(t312);
t301 = Icges(6,4) * t309;
t308 = sin(t312);
t264 = -Icges(6,2) * t308 + t301;
t265 = Icges(6,1) * t308 + t301;
t479 = t264 + t265;
t314 = cos(qJ(4));
t310 = Icges(5,4) * t314;
t313 = sin(qJ(4));
t282 = -Icges(5,2) * t313 + t310;
t283 = Icges(5,1) * t313 + t310;
t478 = t282 + t283;
t311 = pkin(9) + qJ(2);
t307 = qJ(3) + t311;
t303 = cos(t307);
t298 = t303 * pkin(7);
t302 = sin(t307);
t420 = rSges(5,1) * t314;
t352 = pkin(3) + t420;
t393 = t302 * t313;
t355 = rSges(5,2) * t393 + t303 * rSges(5,3);
t182 = -t352 * t302 + t298 + t355;
t414 = t313 * rSges(5,2);
t279 = t303 * t414;
t183 = -t279 + t352 * t303 + (rSges(5,3) + pkin(7)) * t302;
t285 = t313 * rSges(5,1) + rSges(5,2) * t314;
t240 = t285 * t302;
t241 = t285 * t303;
t103 = t182 * t240 - t183 * t241;
t394 = t302 * t309;
t395 = t302 * t308;
t208 = rSges(6,1) * t394 - rSges(6,2) * t395 - t303 * rSges(6,3);
t424 = pkin(4) * t314;
t304 = pkin(3) + t424;
t460 = -pkin(8) - pkin(7);
t356 = -t302 * t304 - t303 * t460;
t175 = -t208 + t356;
t287 = t302 * t460;
t388 = t303 * t308;
t349 = -rSges(6,2) * t388 + rSges(6,3) * t302;
t419 = rSges(6,1) * t309;
t176 = -t287 + (t304 + t419) * t303 + t349;
t267 = rSges(6,1) * t308 + rSges(6,2) * t309;
t425 = pkin(4) * t313;
t329 = t267 + t425;
t470 = t329 * t303;
t471 = t329 * t302;
t90 = t175 * t471 - t176 * t470;
t477 = m(5) * t103 + m(6) * t90;
t462 = m(5) / 0.2e1;
t461 = m(6) / 0.2e1;
t436 = t302 / 0.2e1;
t435 = -t303 / 0.2e1;
t476 = t303 / 0.2e1;
t426 = pkin(2) * cos(t311);
t427 = pkin(2) * sin(t311);
t433 = m(4) * (t426 * (-rSges(4,1) * t302 - rSges(4,2) * t303) + (rSges(4,1) * t303 - t302 * rSges(4,2)) * t427);
t299 = t302 ^ 2;
t300 = t303 ^ 2;
t354 = t299 + t300;
t437 = -t302 / 0.2e1;
t473 = t436 + t437;
t169 = t175 - t427;
t170 = t176 + t426;
t88 = t169 * t471 - t170 * t470;
t230 = t267 * t302;
t231 = t267 * t303;
t166 = -t302 * t230 - t303 * t231;
t268 = -rSges(6,2) * t308 + t419;
t338 = Icges(6,5) * t308 + Icges(6,6) * t309;
t224 = t338 * t302;
t225 = t303 * t338;
t408 = Icges(6,4) * t308;
t266 = Icges(6,1) * t309 - t408;
t205 = Icges(6,5) * t302 + t266 * t303;
t263 = Icges(6,2) * t309 + t408;
t361 = -t263 * t303 + t205;
t203 = Icges(6,6) * t302 + t264 * t303;
t363 = -t265 * t303 - t203;
t323 = -t361 * t308 + t363 * t309;
t271 = Icges(6,4) * t395;
t204 = Icges(6,1) * t394 - Icges(6,5) * t303 - t271;
t362 = -Icges(6,2) * t394 + t204 - t271;
t202 = Icges(6,4) * t394 - Icges(6,2) * t395 - Icges(6,6) * t303;
t364 = t265 * t302 + t202;
t324 = t362 * t308 + t364 * t309;
t423 = (-t299 * t225 + (t324 * t303 + (t224 + t323) * t302) * t303) * t436 + (-t300 * t224 + (t323 * t302 + (t225 + t324) * t303) * t302) * t435;
t387 = t303 * t309;
t147 = t302 * t208 + t303 * (rSges(6,1) * t387 + t349);
t94 = -t302 * (t302 * pkin(3) - t298 + t356) + (-t302 * pkin(7) - t287 + (-pkin(3) + t304) * t303) * t303 + t147;
t15 = t423 + m(6) * (t94 * t166 + (t302 * t471 + t303 * t470) * t268);
t472 = t15 * qJD(5);
t85 = -t176 * t169 + t170 * t175;
t177 = t182 - t427;
t178 = t183 + t426;
t89 = -t183 * t177 + t178 * t182;
t469 = qJD(2) + qJD(3);
t422 = (t88 - t90) * t461 + ((-t178 + t183) * t303 + (t177 - t182) * t302) * t285 * t462;
t101 = t177 * t240 - t178 * t241;
t468 = (t90 + t88) * t461 + (t103 + t101) * t462;
t409 = Icges(5,4) * t313;
t281 = Icges(5,2) * t314 + t409;
t284 = Icges(5,1) * t314 - t409;
t467 = t478 * t314 / 0.2e1 + (t284 / 0.2e1 - t281 / 0.2e1) * t313;
t344 = t479 * t309 / 0.2e1 + (-t263 / 0.2e1 + t266 / 0.2e1) * t308;
t179 = t205 * t394;
t262 = Icges(6,5) * t309 - Icges(6,6) * t308;
t399 = t262 * t303;
t201 = Icges(6,3) * t302 + t399;
t348 = t201 * t303 - t179;
t105 = -t203 * t395 - t348;
t200 = Icges(6,5) * t394 - Icges(6,6) * t395 - Icges(6,3) * t303;
t369 = -t302 * t200 - t204 * t387;
t106 = -t202 * t388 - t369;
t368 = t302 * t201 + t205 * t387;
t107 = -t203 * t388 + t368;
t346 = t203 * t308 - t200;
t401 = t202 * t308;
t350 = ((t105 - t179 + (t201 + t401) * t303 + t369) * t303 + t368 * t302) * t435 + (-t106 * t303 + t107 * t302) * t476 + ((t346 * t302 + t105 + t106 + t348) * t302 + (t107 - t368 + (-t204 * t309 + t401) * t302 + (t346 + t200) * t303) * t303) * t436;
t466 = 4 * qJD(2);
t464 = 2 * qJD(4);
t457 = m(5) * t89;
t92 = t169 * t230 - t170 * t231;
t96 = t175 * t230 - t176 * t231;
t449 = m(6) * (t96 + t92);
t448 = m(6) * ((-t170 + t176) * t303 + (t169 - t175) * t302) * t267;
t333 = (-t169 * t303 - t170 * t302) * t268;
t370 = -t230 * t470 + t231 * t471;
t447 = m(6) * (t333 + t370);
t331 = (t302 * t470 - t303 * t471) * t267;
t446 = m(6) * (t333 + t331);
t332 = (-t175 * t303 - t176 * t302) * t268;
t445 = m(6) * (t332 + t370);
t444 = m(6) * (t332 + t331);
t443 = m(6) * t85;
t441 = m(6) * t88;
t439 = m(6) * t92;
t438 = m(6) * t96;
t392 = t302 * t314;
t217 = Icges(5,5) * t392 - Icges(5,6) * t393 - Icges(5,3) * t303;
t277 = Icges(5,4) * t393;
t221 = Icges(5,1) * t392 - Icges(5,5) * t303 - t277;
t386 = t303 * t314;
t366 = -t302 * t217 - t221 * t386;
t219 = Icges(5,4) * t392 - Icges(5,2) * t393 - Icges(5,6) * t303;
t381 = t313 * t219;
t118 = -t303 * t381 - t366;
t280 = Icges(5,5) * t314 - Icges(5,6) * t313;
t389 = t303 * t280;
t218 = Icges(5,3) * t302 + t389;
t222 = Icges(5,5) * t302 + t284 * t303;
t365 = t302 * t218 + t222 * t386;
t220 = Icges(5,6) * t302 + t282 * t303;
t380 = t313 * t220;
t119 = -t303 * t380 + t365;
t345 = -t217 + t380;
t192 = t222 * t392;
t347 = t303 * t218 - t192;
t28 = (t345 * t303 + t119 - t365) * t303 + (t345 * t302 + t118 + t347) * t302;
t117 = -t302 * t380 - t347;
t29 = (t117 - t192 + (t218 + t381) * t303 + t366) * t303 + t365 * t302;
t81 = -(-(-t221 * t314 + t381) * t302 - t303 * t217) * t303 + t117 * t302;
t82 = -t118 * t303 + t119 * t302;
t2 = (t82 / 0.2e1 - t29 / 0.2e1) * t303 + (t28 / 0.2e1 + t81 / 0.2e1) * t302 + t350;
t434 = t2 * qJD(4);
t431 = m(5) * t101;
t428 = m(6) * t166;
t413 = qJD(5) * t350;
t360 = t283 * t302 + t219;
t359 = -t283 * t303 - t220;
t358 = -Icges(5,2) * t392 + t221 - t277;
t357 = -t281 * t303 + t222;
t351 = -t268 - t424;
t342 = t449 / 0.2e1 + t344;
t339 = Icges(5,5) * t313 + Icges(5,6) * t314;
t330 = (-t240 * t303 + t241 * t302) * t285;
t320 = (-t263 + t266) * t309 - t479 * t308;
t328 = -t350 + (t262 * t302 + t320 * t303 + t363 * t308 + t361 * t309) * t436 + (t320 * t302 - t364 * t308 + t362 * t309 - t399) * t435;
t327 = -t344 + t473 * (t309 * t202 + t308 * t204);
t326 = t344 + t467;
t322 = t358 * t313 + t360 * t314;
t321 = -t357 * t313 + t359 * t314;
t319 = (-t281 + t284) * t314 - t478 * t313;
t318 = t326 + t468;
t317 = (-t230 * t303 + t231 * t302) * t267;
t316 = t327 - t467 + t473 * (t314 * t219 + t313 * t221);
t315 = (t29 * t476 + t328 + (t319 * t302 - t360 * t313 + t358 * t314 - t389 + t82) * t435 + (t28 + t81) * t437 + (t302 * t280 + t319 * t303 + t359 * t313 + t357 * t314) * t436) * qJD(4);
t286 = -t414 + t420;
t235 = t303 * t339;
t234 = t339 * t302;
t216 = t351 * t303;
t214 = t351 * t302;
t171 = -t240 * t302 - t241 * t303;
t159 = qJD(5) * t428;
t144 = -t354 * t425 + t166;
t83 = t344 + t438;
t72 = t344 + t439;
t68 = t444 / 0.2e1;
t67 = t445 / 0.2e1;
t62 = t446 / 0.2e1;
t60 = t447 / 0.2e1;
t53 = t448 / 0.2e1;
t39 = t326 + t477;
t34 = t326 + t431 + t441;
t33 = t433 + t443 + t457;
t19 = -t448 / 0.2e1 + t342;
t18 = t53 + t342;
t17 = m(6) * (t354 * t267 * t268 + t147 * t166) + t423;
t16 = t17 * qJD(5);
t14 = t53 - t449 / 0.2e1 + t327;
t13 = t318 - t422;
t12 = t318 + t422;
t9 = t316 + t422 - t468;
t8 = t67 - t444 / 0.2e1 + t350;
t7 = t68 - t445 / 0.2e1 + t350;
t6 = t60 - t446 / 0.2e1 + t350;
t5 = t62 - t447 / 0.2e1 + t350;
t4 = t67 + t68 + t328;
t3 = t60 + t62 + t328;
t1 = [0, 0, 0, (t144 * t461 + t171 * t462) * t464 + t159, qJD(4) * t428 + t159; 0, qJD(3) * t33 + qJD(4) * t34 + qJD(5) * t72, t33 * qJD(2) + t12 * qJD(4) + t18 * qJD(5) + 0.2e1 * (t433 / 0.2e1 + t85 * t461 + t89 * t462) * qJD(3), t34 * qJD(2) + t12 * qJD(3) + t3 * qJD(5) + ((t169 * t216 + t170 * t214) * t461 + ((-t177 * t303 - t178 * t302) * t286 + t330) * t462) * t464 + t315, t72 * qJD(2) + t18 * qJD(3) + t3 * qJD(4) + ((t333 + t317) * m(6) + t328) * qJD(5); 0, t13 * qJD(4) + t19 * qJD(5) + (-t443 / 0.4e1 - t457 / 0.4e1 - t433 / 0.4e1) * t466, qJD(4) * t39 + qJD(5) * t83, t13 * qJD(2) + t39 * qJD(3) + t4 * qJD(5) + ((t175 * t216 + t176 * t214) * t461 + ((-t182 * t303 - t183 * t302) * t286 + t330) * t462) * t464 + t315, t19 * qJD(2) + t83 * qJD(3) + t4 * qJD(4) + ((t332 + t317) * m(6) + t328) * qJD(5); 0, t316 * qJD(2) + t9 * qJD(3) + t6 * qJD(5) + (-t441 / 0.4e1 - t431 / 0.4e1) * t466 + t434, t9 * qJD(2) + t8 * qJD(5) + t434 + (t316 - t477) * qJD(3), (m(5) * (t285 * t286 * t354 + (t302 * (rSges(5,1) * t392 - t355) + t303 * (rSges(5,1) * t386 + t302 * rSges(5,3) - t279)) * t171) + (-t299 * t235 + (t322 * t303 + (t234 + t321) * t302) * t303) * t436 + (-t300 * t234 + (t321 * t302 + (t235 + t322) * t303) * t302) * t435 + m(6) * (t144 * t94 - t214 * t471 - t216 * t470) + t423) * qJD(4) + t472 + t469 * t2, t6 * qJD(2) + t8 * qJD(3) + t15 * qJD(4) + t472; 0, (t327 - t439) * qJD(2) + t14 * qJD(3) + t5 * qJD(4) + t413, t14 * qJD(2) + (t327 - t438) * qJD(3) + t7 * qJD(4) + t413, t5 * qJD(2) + t7 * qJD(3) + ((t144 * t147 + (-t214 * t302 - t216 * t303) * t267) * m(6) + t423) * qJD(4) + t16, qJD(4) * t17 + t350 * t469 + t16;];
Cq = t1;
