% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:38
% EndTime: 2019-12-31 19:52:46
% DurationCPUTime: 5.27s
% Computational Cost: add. (18239->349), mult. (18026->449), div. (0->0), fcn. (15941->6), ass. (0->228)
t342 = qJ(1) + qJ(2);
t339 = sin(t342);
t340 = cos(t342);
t343 = sin(qJ(4));
t345 = cos(qJ(4));
t455 = rSges(6,1) + pkin(4);
t518 = rSges(6,3) + qJ(5);
t525 = t518 * t343 + t455 * t345;
t533 = t525 * t340;
t535 = t339 * t533;
t414 = t339 * t345;
t374 = t533 * t414;
t366 = rSges(5,1) * t343 + rSges(5,2) * t345;
t352 = -t339 * rSges(5,3) + t366 * t340;
t324 = t340 * qJ(3);
t484 = pkin(2) + pkin(7);
t369 = -t484 * t339 + t324;
t204 = t352 + t369;
t437 = sin(qJ(1)) * pkin(1);
t201 = t204 - t437;
t205 = (rSges(5,3) + t484) * t340 + (qJ(3) + t366) * t339;
t436 = cos(qJ(1)) * pkin(1);
t202 = t205 + t436;
t117 = t201 * t340 + t202 * t339;
t125 = t204 * t340 + t205 * t339;
t519 = pkin(2) - rSges(4,2);
t240 = t340 * rSges(4,3) - t519 * t339 + t324;
t231 = t240 - t437;
t241 = (rSges(4,3) + qJ(3)) * t339 + t519 * t340;
t232 = t241 + t436;
t152 = t231 * t340 + t232 * t339;
t174 = t240 * t340 + t241 * t339;
t485 = m(6) / 0.2e1;
t486 = m(5) / 0.2e1;
t372 = t455 * t343;
t410 = t340 * t345;
t373 = -t339 * rSges(6,2) - t518 * t410;
t178 = t340 * t372 + t369 + t373;
t381 = t518 * t414;
t179 = (rSges(6,2) + t484) * t340 + (qJ(3) + t372) * t339 - t381;
t499 = -t178 * t340 - t179 * t339;
t166 = t178 - t437;
t167 = t179 + t436;
t506 = -t166 * t340 - t167 * t339;
t524 = m(4) / 0.2e1;
t378 = (-t499 - t506) * t485 + (t125 + t117) * t486 + (t174 + t152) * t524;
t379 = (-t506 + t499) * t485 + (t117 - t125) * t486 + (t152 - t174) * t524;
t6 = t379 - t378;
t534 = t6 * qJD(1);
t430 = Icges(5,4) * t343;
t364 = Icges(5,2) * t345 + t430;
t428 = Icges(6,5) * t343;
t532 = Icges(6,3) * t345 + t364 - t428;
t302 = Icges(5,1) * t345 - t430;
t531 = Icges(6,1) * t345 + t302 + t428;
t530 = t166 - t178;
t529 = t167 - t179;
t341 = Icges(6,5) * t345;
t294 = Icges(6,3) * t343 + t341;
t299 = -Icges(6,1) * t343 + t341;
t528 = t294 + t299;
t429 = Icges(5,4) * t345;
t298 = -Icges(5,2) * t343 + t429;
t365 = Icges(5,1) * t343 + t429;
t527 = t298 + t365;
t526 = (Icges(6,4) + Icges(5,5)) * t345 + (-Icges(5,6) + Icges(6,6)) * t343;
t453 = m(3) * (t436 * (-rSges(3,1) * t339 - rSges(3,2) * t340) + (t340 * rSges(3,1) - t339 * rSges(3,2)) * t437);
t449 = m(4) * (-t241 * t231 + t232 * t240);
t337 = t339 ^ 2;
t338 = t340 ^ 2;
t380 = t337 + t338;
t237 = (0.1e1 - t380) * t345 * t343;
t433 = m(6) * qJD(5);
t517 = t237 * t433;
t360 = Icges(5,5) * t343 + Icges(5,6) * t345;
t362 = Icges(6,4) * t343 - Icges(6,6) * t345;
t516 = -t360 - t362;
t317 = Icges(6,5) * t410;
t411 = t340 * t343;
t251 = Icges(6,1) * t411 - Icges(6,4) * t339 - t317;
t253 = -Icges(5,5) * t339 + t365 * t340;
t503 = -Icges(6,3) * t411 + t298 * t340 + t251 + t253 - t317;
t250 = Icges(6,4) * t340 - t299 * t339;
t318 = Icges(5,4) * t414;
t415 = t339 * t343;
t252 = Icges(5,1) * t415 + Icges(5,5) * t340 + t318;
t502 = -Icges(5,2) * t415 - t294 * t339 + t250 + t252 + t318;
t243 = Icges(6,5) * t411 - Icges(6,6) * t339 - Icges(6,3) * t410;
t249 = -Icges(5,6) * t339 + t364 * t340;
t501 = -t340 * t531 - t243 + t249;
t316 = Icges(6,5) * t415;
t242 = Icges(6,6) * t340 - Icges(6,3) * t414 + t316;
t248 = Icges(5,6) * t340 + t364 * t339;
t500 = -Icges(6,1) * t414 - t302 * t339 - t242 + t248 - t316;
t514 = (t531 - t532) * t345 + (t528 - t527) * t343;
t512 = t343 * t500 - t345 * t502;
t511 = t343 * t501 - t345 * t503;
t510 = m(6) * t345;
t509 = t366 * t486;
t508 = (t249 * t345 + t253 * t343) * t340;
t507 = (t243 * t345 - t251 * t343) * t340;
t58 = -t179 * t166 + t167 * t178;
t88 = -t205 * t201 + t202 * t204;
t505 = t526 * t339;
t504 = t526 * t340;
t498 = qJD(1) + qJD(2);
t207 = t455 * t414 + t518 * t415;
t309 = rSges(5,1) * t345 - rSges(5,2) * t343;
t268 = t309 * t339;
t270 = t309 * t340;
t399 = (-t207 * t340 + t535) * t485 + (-t268 * t340 + t270 * t339) * t486;
t234 = t525 * t339;
t422 = t205 * t268;
t424 = t202 * t268;
t435 = (t529 * t234 + t530 * t533) * t485 + (t424 - t422 + (t201 - t204) * t270) * t486;
t105 = t204 * t270 + t422;
t80 = t166 * t533 + t167 * t207;
t84 = t178 * t533 + t179 * t207;
t98 = t201 * t270 + t424;
t496 = (t84 + t80) * t485 + (t105 + t98) * t486;
t351 = (-t527 / 0.2e1 + t528 / 0.2e1) * t345 + (-t531 / 0.2e1 + t532 / 0.2e1) * t343;
t492 = 0.4e1 * qJD(1);
t490 = 0.4e1 * qJD(2);
t489 = 2 * qJD(4);
t479 = m(5) * t88;
t477 = m(5) * t98;
t473 = (-t529 * t339 - t530 * t340) * t510;
t472 = ((-t166 - t178) * t340 + (-t167 - t179) * t339) * t510;
t469 = m(6) * t58;
t151 = t166 * t415;
t368 = t207 * t410 - t374;
t467 = m(6) * (-t167 * t411 + t151 + t368);
t157 = t178 * t415;
t466 = m(6) * (-t179 * t411 + t157 + t368);
t420 = t234 * t345;
t465 = m(6) * (t374 + t151 + (-t167 * t343 - t420) * t340);
t464 = m(6) * (t374 + t157 + (-t179 * t343 - t420) * t340);
t463 = m(6) * t80;
t462 = m(6) * t84;
t460 = m(6) * t506;
t458 = t339 / 0.2e1;
t457 = -t340 / 0.2e1;
t456 = t340 / 0.2e1;
t247 = Icges(6,4) * t411 - Icges(6,2) * t339 - Icges(6,6) * t410;
t127 = t243 * t414 - t340 * t247 - t251 * t415;
t246 = Icges(6,2) * t340 + t362 * t339;
t396 = t242 * t410 + t339 * t246;
t130 = -t250 * t411 + t396;
t131 = -t247 * t339 - t507;
t370 = -t242 * t345 + t247;
t397 = t340 * t246 + t250 * t415;
t417 = t250 * t343;
t23 = (t131 + t397 + t507) * t340 + (-t130 + (t370 - t417) * t340 + t127 + t396) * t339;
t244 = Icges(5,3) * t340 + t360 * t339;
t128 = t340 * t244 + t248 * t414 + t252 * t415;
t245 = -Icges(5,3) * t339 + t360 * t340;
t129 = -t340 * t245 - t249 * t414 - t253 * t415;
t225 = t339 * t244;
t356 = -t248 * t345 - t252 * t343;
t132 = t356 * t340 + t225;
t133 = -t245 * t339 + t508;
t24 = (-t132 + t225 + t129) * t339 + (t133 - t508 + (t245 + t356) * t339 + t128) * t340;
t126 = -t242 * t414 + t397;
t25 = (t127 + (t247 + t417) * t340 - t396) * t340 + (t370 * t339 - t126 + t397) * t339;
t26 = t245 * t337 + (t129 - t225 + (t245 - t356) * t340) * t340;
t76 = t126 * t340 + t127 * t339;
t77 = t128 * t340 + t129 * t339;
t78 = t130 * t340 + t131 * t339;
t79 = t132 * t340 + t133 * t339;
t2 = (t26 / 0.2e1 + t79 / 0.2e1 + t25 / 0.2e1 + t78 / 0.2e1) * t340 + (-t77 / 0.2e1 + t24 / 0.2e1 - t76 / 0.2e1 + t23 / 0.2e1) * t339;
t439 = m(6) * (t234 * t340 - t535);
t149 = t439 / 0.2e1;
t47 = t149 - t399;
t454 = t47 * qJD(3) + t2 * qJD(4);
t447 = m(4) * t152;
t446 = m(4) * t174;
t445 = m(5) * t105;
t444 = m(5) * t117;
t443 = m(5) * t125;
t441 = m(6) * t499;
t434 = m(6) * qJD(4);
t398 = (-t207 + t234) * t533;
t383 = t518 * t345 - t372;
t90 = t506 * t345;
t377 = m(6) * t90 * qJD(1);
t94 = t499 * t345;
t376 = m(6) * t94 * qJD(2);
t371 = t380 * t366;
t233 = t383 * t339;
t235 = t383 * t340;
t359 = t233 * t339 + t235 * t340;
t349 = t351 + t496;
t348 = -t351 + ((t251 - t253) * t345 + (t243 + t249) * t343) * (t456 + t457);
t45 = t149 + t399;
t347 = t45 * qJD(3) + (-(t23 + t24) * t339 / 0.2e1 + (t25 + t26 + t78 + t79) * t457 + (t339 * t514 + t340 * t516 - t343 * t502 - t345 * t500) * t456 + (t339 * t516 - t340 * t514 + t343 * t503 + t345 * t501 + t76 + t77) * t458) * qJD(4);
t273 = t380 * t343;
t206 = t234 * t415;
t136 = -t207 * t339 - t338 * t525;
t116 = (-t411 * t455 - t373) * t340 + (-t340 * rSges(6,2) - t455 * t415 + t381) * t339;
t83 = t116 * t345 * t380 + t411 * t533 + t206;
t68 = t464 / 0.2e1;
t65 = t465 / 0.2e1;
t64 = t466 / 0.2e1;
t62 = t467 / 0.2e1;
t49 = t472 / 0.2e1;
t48 = t473 / 0.2e1;
t46 = -t439 / 0.2e1 + t399;
t43 = t46 * qJD(4);
t41 = t45 * qJD(4);
t40 = -t441 + t443 + t446;
t39 = t444 + t447 - t460;
t30 = t351 + t445 + t462;
t29 = t351 + t463 + t477;
t18 = t449 + t453 + t469 + t479;
t17 = t68 - t466 / 0.2e1;
t16 = t68 + t64;
t15 = t64 - t464 / 0.2e1;
t14 = t65 - t467 / 0.2e1;
t13 = t65 + t62;
t12 = t62 - t465 / 0.2e1;
t11 = t49 - t473 / 0.2e1;
t10 = t49 + t48;
t9 = t48 - t472 / 0.2e1;
t7 = t378 + t379;
t5 = t349 + t435;
t4 = t349 - t435;
t3 = t348 + t435 - t496;
t1 = [t18 * qJD(2) + t39 * qJD(3) + t29 * qJD(4) + t90 * t433, t18 * qJD(1) + t7 * qJD(3) + t5 * qJD(4) + t10 * qJD(5) + 0.2e1 * (t449 / 0.2e1 + t453 / 0.2e1 + t58 * t485 + t88 * t486) * qJD(2), qJD(1) * t39 + qJD(2) * t7 + t41, t29 * qJD(1) + t5 * qJD(2) + t13 * qJD(5) + ((t166 * t233 - t167 * t235 + t398) * t485 - (t201 * t339 - t202 * t340) * t509) * t489 + t347, t10 * qJD(2) + t13 * qJD(4) + t377; -t6 * qJD(3) + t4 * qJD(4) + t11 * qJD(5) + (-t469 / 0.4e1 - t479 / 0.4e1 - t449 / 0.4e1 - t453 / 0.4e1) * t492, t40 * qJD(3) + t30 * qJD(4) + t94 * t433, qJD(2) * t40 + t41 - t534, t4 * qJD(1) + t30 * qJD(2) + t16 * qJD(5) + ((t178 * t233 - t179 * t235 + t398) * t485 - (t204 * t339 - t205 * t340) * t509) * t489 + t347, t11 * qJD(1) + t16 * qJD(4) + t376; t6 * qJD(2) + t43 + (t460 / 0.4e1 - t444 / 0.4e1 - t447 / 0.4e1) * t492, t534 + t43 + (t441 / 0.4e1 - t443 / 0.4e1 - t446 / 0.4e1) * t490, 0, (t359 * t485 - t371 * t486) * t489 + t273 * t433 + t498 * t46, t273 * t434; t3 * qJD(2) + t14 * qJD(5) + (-t463 / 0.4e1 - t477 / 0.4e1) * t492 + t454 + t348 * qJD(1), t3 * qJD(1) + t348 * qJD(2) + t17 * qJD(5) + (-t462 / 0.4e1 - t445 / 0.4e1) * t490 + t454, t498 * t47, (m(5) * ((-t340 * t352 + (-t340 * rSges(5,3) - t366 * t339) * t339) * (-t268 * t339 - t270 * t340) - t309 * t371) + m(6) * (t116 * t136 + t233 * t234 + t235 * t533) + (-t504 * t337 + (t512 * t340 + (t505 - t511) * t339) * t340) * t458 + (t505 * t338 + (t511 * t339 + (-t504 - t512) * t340) * t339) * t456) * qJD(4) + t83 * t433 + t498 * t2, t14 * qJD(1) + t17 * qJD(2) + t83 * t434 - t517; t9 * qJD(2) + t12 * qJD(4) - t377, t9 * qJD(1) + t15 * qJD(4) - t376, 0, t12 * qJD(1) + t15 * qJD(2) + (t206 + (t340 * t533 + t136) * t343 + (t116 - t359) * t345 - t83) * t434 + t517, t237 * t434;];
Cq = t1;
