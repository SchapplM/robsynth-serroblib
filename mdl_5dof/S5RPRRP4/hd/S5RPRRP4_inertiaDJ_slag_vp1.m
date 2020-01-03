% Calculate time derivative of joint inertia matrix for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:16
% EndTime: 2020-01-03 11:49:54
% DurationCPUTime: 18.96s
% Computational Cost: add. (15761->586), mult. (22702->817), div. (0->0), fcn. (21282->8), ass. (0->266)
t468 = Icges(5,4) + Icges(6,4);
t466 = Icges(5,2) + Icges(6,2);
t465 = Icges(5,6) + Icges(6,6);
t297 = qJ(3) + qJ(4);
t289 = cos(t297);
t470 = t468 * t289;
t469 = Icges(5,1) + Icges(6,1);
t467 = Icges(5,5) + Icges(6,5);
t464 = -Icges(6,3) - Icges(5,3);
t288 = sin(t297);
t298 = sin(pkin(8));
t299 = cos(pkin(8));
t457 = -t465 * t299 + (-t466 * t288 + t470) * t298;
t463 = t468 * t288;
t296 = qJD(3) + qJD(4);
t388 = t296 * t298;
t459 = (t467 * t288 + t465 * t289) * t388;
t462 = (-t469 * t288 - t470) * t388;
t451 = t462 * t289 * t298 + t459 * t299;
t452 = t467 * t299 + (-t469 * t289 + t463) * t298;
t455 = (t466 * t289 + t463) * t388;
t390 = t289 * t296;
t456 = t457 * t390;
t443 = t451 + ((t452 * t296 + t455) * t288 - t456) * t298;
t461 = t443 * t299;
t303 = cos(qJ(1));
t301 = sin(qJ(1));
t384 = t299 * t301;
t244 = -t288 * t384 - t289 * t303;
t383 = t299 * t303;
t311 = t288 * t301 + t289 * t383;
t183 = qJD(1) * t244 + t296 * t311;
t245 = -t288 * t303 + t289 * t384;
t246 = t288 * t383 - t289 * t301;
t184 = qJD(1) * t245 + t246 * t296;
t360 = qJD(1) * t301;
t342 = t298 * t360;
t460 = -t465 * t183 - t467 * t184 + t464 * t342;
t185 = -qJD(1) * t246 - t245 * t296;
t186 = qJD(1) * t311 + t244 * t296;
t359 = qJD(1) * t303;
t341 = t298 * t359;
t442 = t465 * t185 + t467 * t186 - t464 * t341;
t441 = t466 * t183 + t468 * t184 + t465 * t342;
t440 = t466 * t185 + t468 * t186 + t465 * t341;
t439 = -t468 * t183 - t469 * t184 - t467 * t342;
t438 = t468 * t185 + t469 * t186 + t467 * t341;
t387 = t298 * t301;
t436 = t465 * t244 + t467 * t245 - t464 * t387;
t386 = t298 * t303;
t435 = t465 * t246 - t467 * t311 + t464 * t386;
t434 = t466 * t244 + t468 * t245 + t465 * t387;
t433 = t466 * t246 - t468 * t311 - t465 * t386;
t432 = t468 * t244 + t469 * t245 + t467 * t387;
t431 = t468 * t246 - t469 * t311 - t467 * t386;
t454 = t464 * t299 + (-t465 * t288 + t467 * t289) * t298;
t450 = t442 * t299 + ((t434 * t296 - t438) * t289 + (t432 * t296 + t440) * t288) * t298;
t449 = t460 * t299 + ((-t433 * t296 - t439) * t289 + (-t431 * t296 - t441) * t288) * t298;
t448 = -t462 * t311 + (t459 * t303 + t454 * t360) * t298 - t455 * t246 - t452 * t184 + t457 * t183;
t447 = t434 * t246 - t432 * t311 - t436 * t386;
t446 = -t433 * t246 + t431 * t311 + t435 * t386;
t445 = t436 * t299 + (t434 * t288 - t432 * t289) * t298;
t444 = t435 * t299 + (t433 * t288 - t431 * t289) * t298;
t300 = sin(qJ(3));
t357 = qJD(3) * t300;
t353 = pkin(3) * t357;
t407 = pkin(4) * t288;
t255 = -t296 * t407 - t353;
t302 = cos(qJ(3));
t285 = t302 * pkin(3) + pkin(2);
t265 = pkin(4) * t289 + t285;
t408 = pkin(3) * t300;
t268 = t407 + t408;
t339 = t299 * t359;
t355 = qJD(5) * t298;
t437 = t186 * rSges(6,1) + t185 * rSges(6,2) + rSges(6,3) * t341 + t255 * t384 + t265 * t339 + t268 * t360 + t301 * t355;
t430 = t245 * rSges(6,1) + t244 * rSges(6,2) + rSges(6,3) * t387 + t265 * t384;
t429 = -rSges(6,1) * t311 + rSges(6,2) * t246 - t301 * t268;
t356 = qJD(3) * t302;
t352 = pkin(3) * t356;
t256 = pkin(4) * t390 + t352;
t428 = -rSges(6,1) * t184 - rSges(6,2) * t183 + t256 * t301;
t427 = t460 * t303;
t426 = (t459 * t301 - t454 * t359) * t298 - t462 * t245 + t455 * t244 + t452 * t186 - t457 * t185;
t425 = t434 * t244 + t432 * t245 + t387 * t436;
t424 = t244 * t433 + t245 * t431 + t387 * t435;
t423 = -t244 * t457 + t245 * t452 - t387 * t454;
t422 = -t246 * t457 - t311 * t452 + t386 * t454;
t304 = -pkin(7) - pkin(6);
t294 = -qJ(5) + t304;
t323 = t255 + t353;
t354 = qJD(1) * t408;
t361 = qJD(1) * t298;
t416 = t304 * t361 + t352;
t343 = t301 * t416 + t303 * t354;
t367 = t265 - t285;
t376 = t303 * t268;
t421 = rSges(6,3) * t342 + (-t299 * t323 - t355) * t303 + (-t376 + (-t294 * t298 + t299 * t367) * t301) * qJD(1) + t343 - t428;
t328 = t301 * t353;
t309 = -t285 * t339 + t299 * t328 - t301 * t354;
t362 = t294 - t304;
t332 = t362 * t298;
t420 = (-qJD(1) * t332 - t256 + t352) * t303 + t309 + t437;
t263 = t285 * t384;
t419 = -t263 + (-t268 + t408) * t303 - t301 * t332 + t430;
t381 = t300 * t301;
t368 = pkin(3) * t381 + t285 * t383;
t392 = t265 * t299;
t418 = (-t392 + t332) * t303 + t368 - rSges(6,3) * t386 + t429;
t417 = (t362 - rSges(6,3)) * t299 + (rSges(6,1) * t289 - rSges(6,2) * t288 + t367) * t298;
t415 = t461 + (t449 * t303 + t450 * t301 + (t444 * t301 + t445 * t303) * qJD(1)) * t298;
t414 = (rSges(6,3) - t294) * t298 + t392;
t413 = (t426 * t299 + (t424 * t360 + t425 * t359 + (-t433 * t185 - t431 * t186 - t441 * t244 + t439 * t245 - t435 * t341) * t303 + ((t442 * t301 + t436 * t359 + t427) * t298 + t438 * t245 + t440 * t244 + t432 * t186 + t434 * t185) * t301) * t298) * t387 + (t422 * t299 + (t447 * t301 + t446 * t303) * t298) * t342 + (t423 * t299 + (t301 * t425 - t303 * t424) * t298) * t341 + (t448 * t299 + (t446 * t360 - t447 * t359 + (t439 * t311 + (t435 * t360 + t427) * t298 + t441 * t246 + t431 * t184 + t433 * t183) * t303 + ((t442 * t303 - t436 * t360) * t298 + t438 * t311 - t440 * t246 - t432 * t184 - t434 * t183) * t301) * t298) * t386;
t412 = 2 * m(4);
t411 = 2 * m(5);
t410 = 2 * m(6);
t295 = t298 ^ 2;
t409 = pkin(2) * t299;
t406 = pkin(6) * t298;
t405 = -pkin(2) + t285;
t404 = pkin(6) + t304;
t400 = Icges(4,4) * t300;
t399 = Icges(4,4) * t302;
t223 = (-rSges(5,1) * t288 - rSges(5,2) * t289) * t388;
t394 = t223 * t301;
t385 = t298 * t304;
t358 = qJD(3) * t298;
t249 = (-Icges(4,2) * t302 - t400) * t358;
t382 = t300 * t249;
t380 = t300 * t303;
t378 = t301 * t302;
t377 = t302 * t303;
t375 = t418 * t299;
t173 = t245 * rSges(5,1) + t244 * rSges(5,2) + rSges(5,3) * t387;
t316 = rSges(5,1) * t311 - t246 * rSges(5,2);
t175 = -rSges(5,3) * t386 - t316;
t102 = t173 * t386 + t175 * t387;
t312 = -t298 * t404 - t409;
t195 = -pkin(3) * t380 + t301 * t312 + t263;
t373 = -t173 - t195;
t366 = pkin(2) * t383 + pkin(6) * t386;
t196 = t303 * t385 + t366 - t368;
t372 = t195 * t386 + t196 * t387;
t370 = qJD(5) * t299 - t298 * t323 - (-rSges(6,1) * t288 - rSges(6,2) * t289) * t388;
t224 = t298 * t405 + t299 * t404;
t233 = -rSges(5,3) * t299 + (rSges(5,1) * t289 - rSges(5,2) * t288) * t298;
t369 = -t224 - t233;
t365 = pkin(1) * t359 + qJ(2) * t360;
t364 = qJ(2) * t359 + qJD(2) * t301;
t363 = t303 * pkin(1) + t301 * qJ(2);
t317 = -rSges(5,1) * t184 - rSges(5,2) * t183;
t118 = rSges(5,3) * t342 - t317;
t120 = t186 * rSges(5,1) + t185 * rSges(5,2) + rSges(5,3) * t341;
t351 = t118 * t387 + t120 * t386 + t175 * t341;
t327 = t303 * t353;
t308 = t299 * t327 - t343;
t148 = (t299 * t405 - t406) * t360 + t308;
t149 = (qJD(1) * t312 - t352) * t303 - t309;
t347 = t148 * t387 + t149 * t386 + t196 * t341;
t346 = t299 * t148 + t224 * t342 + t295 * t327;
t345 = -t195 - t419;
t344 = -t224 - t417;
t259 = t299 * t378 - t380;
t260 = t299 * t380 - t378;
t206 = -qJD(1) * t260 - qJD(3) * t259;
t258 = -t299 * t381 - t377;
t310 = t299 * t377 + t381;
t207 = qJD(1) * t310 + qJD(3) * t258;
t134 = t207 * rSges(4,1) + t206 * rSges(4,2) + rSges(4,3) * t341;
t197 = t259 * rSges(4,1) + t258 * rSges(4,2) + rSges(4,3) * t387;
t335 = t417 * t303;
t334 = t370 * t301;
t333 = t369 * t303;
t248 = (-Icges(4,5) * t300 - Icges(4,6) * t302) * t358;
t250 = (-Icges(4,1) * t300 - t399) * t358;
t329 = t298 * t302 * t250 - t299 * t248;
t42 = t386 * t419 + t387 * t418;
t322 = t344 * t303;
t321 = t406 + t409;
t319 = rSges(3,1) * t299 - rSges(3,2) * t298;
t204 = qJD(1) * t258 + qJD(3) * t310;
t205 = qJD(1) * t259 + qJD(3) * t260;
t318 = -rSges(4,1) * t205 - rSges(4,2) * t204;
t313 = t341 * t418 + t386 * t420 + t387 * t421;
t62 = t299 * t118 - t223 * t386 + t233 * t342;
t198 = -rSges(4,1) * t310 + rSges(4,2) * t260 - rSges(4,3) * t386;
t24 = t299 * t421 + t342 * t417 + t370 * t386;
t307 = rSges(3,3) * t301 + t303 * t319;
t306 = t415 * t299 + t413;
t305 = (-t426 - t450) * t387 / 0.2e1 - (t448 + t449) * t386 / 0.2e1 + ((-t422 - t444) * t301 + (-t423 - t445) * t303) * t361 / 0.2e1;
t291 = t301 * pkin(1);
t266 = t295 * t328;
t251 = (-rSges(4,1) * t300 - rSges(4,2) * t302) * t358;
t242 = -rSges(4,3) * t299 + (rSges(4,1) * t302 - rSges(4,2) * t300) * t298;
t241 = -Icges(4,5) * t299 + (Icges(4,1) * t302 - t400) * t298;
t240 = -Icges(4,6) * t299 + (-Icges(4,2) * t300 + t399) * t298;
t239 = -Icges(4,3) * t299 + (Icges(4,5) * t302 - Icges(4,6) * t300) * t298;
t215 = t307 + t363;
t214 = t291 + (-rSges(3,3) - qJ(2)) * t303 + t319 * t301;
t200 = qJD(1) * t307 - qJD(2) * t303 + t365;
t199 = (rSges(3,3) * t303 + (-pkin(1) - t319) * t301) * qJD(1) + t364;
t194 = -Icges(4,1) * t310 + Icges(4,4) * t260 - Icges(4,5) * t386;
t193 = Icges(4,1) * t259 + Icges(4,4) * t258 + Icges(4,5) * t387;
t192 = -Icges(4,4) * t310 + Icges(4,2) * t260 - Icges(4,6) * t386;
t191 = Icges(4,4) * t259 + Icges(4,2) * t258 + Icges(4,6) * t387;
t190 = -Icges(4,5) * t310 + Icges(4,6) * t260 - Icges(4,3) * t386;
t189 = Icges(4,5) * t259 + Icges(4,6) * t258 + Icges(4,3) * t387;
t188 = t299 * t196;
t159 = t299 * t175;
t151 = -t198 + t363 + t366;
t150 = -t303 * qJ(2) + t301 * t321 + t197 + t291;
t144 = t198 * t299 - t242 * t386;
t143 = -t197 * t299 - t242 * t387;
t136 = -t233 * t386 + t159;
t135 = -t173 * t299 - t233 * t387;
t133 = rSges(4,3) * t342 - t318;
t132 = Icges(4,1) * t207 + Icges(4,4) * t206 + Icges(4,5) * t341;
t131 = Icges(4,1) * t205 + Icges(4,4) * t204 + Icges(4,5) * t342;
t130 = Icges(4,4) * t207 + Icges(4,2) * t206 + Icges(4,6) * t341;
t129 = Icges(4,4) * t205 + Icges(4,2) * t204 + Icges(4,6) * t342;
t128 = Icges(4,5) * t207 + Icges(4,6) * t206 + Icges(4,3) * t341;
t127 = Icges(4,5) * t205 + Icges(4,6) * t204 + Icges(4,3) * t342;
t126 = (rSges(5,3) - t304) * t386 + t316 + t363 + t368;
t125 = -t301 * t385 + t263 + t291 + (-qJ(2) - t408) * t303 + t173;
t124 = -t239 * t386 + t240 * t260 - t241 * t310;
t123 = t239 * t387 + t240 * t258 + t241 * t259;
t122 = t414 * t303 + t363 - t429;
t121 = -t294 * t387 + t291 + (-qJ(2) - t268) * t303 + t430;
t97 = (qJD(1) * t321 - qJD(2)) * t303 + t134 + t365;
t96 = (-t409 - pkin(1) + (-rSges(4,3) - pkin(6)) * t298) * t360 + t318 + t364;
t91 = (-t382 + (-t240 * t302 - t241 * t300) * qJD(3)) * t298 + t329;
t86 = -t134 * t299 + (-t242 * t359 - t251 * t301) * t298;
t85 = t133 * t299 + (t242 * t360 - t251 * t303) * t298;
t84 = -t190 * t299 + (-t192 * t300 + t194 * t302) * t298;
t83 = -t189 * t299 + (-t191 * t300 + t193 * t302) * t298;
t77 = t298 * t333 + t159 + t188;
t76 = t299 * t373 + t369 * t387;
t75 = -t190 * t386 + t192 * t260 - t194 * t310;
t74 = -t189 * t386 + t191 * t260 - t193 * t310;
t73 = t190 * t387 + t192 * t258 + t194 * t259;
t72 = t189 * t387 + t191 * t258 + t193 * t259;
t65 = (-qJD(2) - t416) * t303 - t309 + t120 + t365;
t64 = (-rSges(5,3) * t298 - t285 * t299 - pkin(1)) * t360 - t308 + t317 + t364;
t63 = -t120 * t299 + (-t233 * t359 - t394) * t298;
t49 = -t298 * t335 + t375;
t48 = -t299 * t419 - t387 * t417;
t47 = t372 + t102;
t46 = (-t294 * t361 - qJD(2) - t256) * t303 + t365 + t437;
t45 = (t255 * t299 + t355) * t303 + (t376 + (-pkin(1) - t414) * t301) * qJD(1) + t364 + t428;
t44 = t206 * t240 + t207 * t241 + t249 * t258 + t250 * t259 + (t239 * t359 + t248 * t301) * t298;
t43 = t204 * t240 + t205 * t241 + t249 * t260 - t250 * t310 + (t239 * t360 - t248 * t303) * t298;
t41 = t298 * t322 + t188 + t375;
t40 = t299 * t345 + t344 * t387;
t39 = t266 + (-t120 - t149) * t299 + (qJD(1) * t333 - t394) * t298;
t38 = t62 + t346;
t29 = -t173 * t342 + t351;
t28 = t42 + t372;
t27 = -t127 * t299 + (-t129 * t300 + t131 * t302 + (-t192 * t302 - t194 * t300) * qJD(3)) * t298;
t26 = -t128 * t299 + (-t130 * t300 + t132 * t302 + (-t191 * t302 - t193 * t300) * qJD(3)) * t298;
t25 = -t420 * t299 + (-qJD(1) * t335 + t334) * t298;
t15 = t266 + (-t149 - t420) * t299 + (qJD(1) * t322 + t334) * t298;
t14 = t24 + t346;
t9 = t342 * t373 + t347 + t351;
t8 = -t342 * t419 + t313;
t7 = t342 * t345 + t313 + t347;
t1 = [0.2e1 * m(3) * (t199 * t215 + t200 * t214) + (t150 * t97 + t151 * t96) * t412 + (t125 * t65 + t126 * t64) * t411 + (t121 * t46 + t122 * t45) * t410 + t329 + (-t240 * t356 - t241 * t357 - t382 - t456) * t298 + (t298 * t455 + t388 * t452) * t288 + t451; m(3) * (-t199 * t303 - t200 * t301 + (-t214 * t303 + t215 * t301) * qJD(1)) + m(4) * (-t301 * t97 - t303 * t96 + (-t150 * t303 + t151 * t301) * qJD(1)) + m(5) * (-t301 * t65 - t303 * t64 + (-t125 * t303 + t126 * t301) * qJD(1)) + m(6) * (-t301 * t46 - t303 * t45 + (-t121 * t303 + t122 * t301) * qJD(1)); 0; m(6) * (t121 * t15 + t122 * t14 + t40 * t46 + t41 * t45) + m(5) * (t125 * t39 + t126 * t38 + t64 * t77 + t65 * t76) + m(4) * (t143 * t97 + t144 * t96 + t150 * t86 + t151 * t85) + ((-t27 / 0.2e1 - t43 / 0.2e1) * t303 + (t26 / 0.2e1 + t44 / 0.2e1) * t301 + ((t83 / 0.2e1 + t123 / 0.2e1) * t303 + (t84 / 0.2e1 + t124 / 0.2e1) * t301) * qJD(1)) * t298 + t305 + (-t91 - t443) * t299; m(4) * (-t301 * t86 - t303 * t85 + (-t143 * t303 + t144 * t301) * qJD(1)) + m(5) * (-t301 * t39 - t303 * t38 + (t301 * t77 - t303 * t76) * qJD(1)) + m(6) * (-t14 * t303 - t15 * t301 + (t301 * t41 - t303 * t40) * qJD(1)); (t143 * t86 + t144 * t85 + (t197 * t303 + t198 * t301) * (t133 * t301 + t134 * t303 + (-t197 * t301 + t198 * t303) * qJD(1)) * t295) * t412 + (t14 * t41 + t15 * t40 + t28 * t7) * t410 + (t38 * t77 + t39 * t76 + t47 * t9) * t411 + ((t301 * t72 - t303 * t73) * t341 + (t301 * t74 - t303 * t75) * t342 + ((t258 * t130 + t259 * t132 + t206 * t191 + t207 * t193) * t301 + t72 * t359 - (t258 * t129 + t259 * t131 + t206 * t192 + t207 * t194) * t303 + t73 * t360) * t387 - ((t260 * t130 - t132 * t310 + t204 * t191 + t205 * t193) * t301 + t74 * t359 - (t260 * t129 - t131 * t310 + t204 * t192 + t205 * t194) * t303 + t75 * t360) * t386 + (((t128 * t301 + t189 * t359) * t301 - (t127 * t301 + t190 * t359) * t303) * t387 - ((-t128 * t303 + t189 * t360) * t301 - (-t127 * t303 + t190 * t360) * t303) * t386) * t298) * t298 + (-(t26 * t301 - t27 * t303 + (t301 * t84 + t303 * t83) * qJD(1)) * t298 - t123 * t341 - t124 * t342 - t44 * t387 + t43 * t386 + t91 * t299 + t415) * t299 + t413; t305 + m(5) * (t125 * t63 + t126 * t62 + t135 * t65 + t136 * t64) + m(6) * (t121 * t25 + t122 * t24 + t49 * t45 + t48 * t46) - t461; m(5) * (-t301 * t63 - t303 * t62 + (-t135 * t303 + t136 * t301) * qJD(1)) + m(6) * (-t24 * t303 - t25 * t301 + (t301 * t49 - t303 * t48) * qJD(1)); m(6) * (t14 * t49 + t15 * t48 + t24 * t41 + t25 * t40 + t28 * t8 + t42 * t7) + m(5) * (t102 * t9 + t135 * t39 + t136 * t38 + t29 * t47 + t62 * t77 + t63 * t76) + t306; (t102 * t29 + t135 * t63 + t136 * t62) * t411 + (t24 * t49 + t25 * t48 + t42 * t8) * t410 + t306; m(6) * (t301 * t45 - t303 * t46 + (t121 * t301 + t122 * t303) * qJD(1)) * t298; 0; m(6) * (-t299 * t7 + (t14 * t301 - t15 * t303 + (t301 * t40 + t303 * t41) * qJD(1)) * t298); m(6) * (-t299 * t8 + (t24 * t301 - t25 * t303 + (t301 * t48 + t303 * t49) * qJD(1)) * t298); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
