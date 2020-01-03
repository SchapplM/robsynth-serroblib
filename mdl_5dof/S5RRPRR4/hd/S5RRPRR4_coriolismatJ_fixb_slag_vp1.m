% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:01:49
% EndTime: 2020-01-03 12:01:56
% DurationCPUTime: 4.48s
% Computational Cost: add. (42777->354), mult. (23906->459), div. (0->0), fcn. (21722->10), ass. (0->238)
t318 = qJ(4) + qJ(5);
t312 = sin(t318);
t314 = cos(t318);
t273 = rSges(6,1) * t312 + rSges(6,2) * t314;
t319 = qJ(1) + qJ(2);
t311 = pkin(9) + t319;
t308 = cos(t311);
t239 = t273 * t308;
t307 = sin(t311);
t320 = sin(qJ(4));
t434 = pkin(4) * t320;
t350 = t273 + t434;
t491 = t350 * t307;
t179 = t491 * t239;
t417 = Icges(6,4) * t312;
t272 = Icges(6,1) * t314 - t417;
t209 = -Icges(6,5) * t308 + t272 * t307;
t392 = t307 * t314;
t192 = t209 * t392;
t386 = t308 * t314;
t387 = t308 * t312;
t206 = Icges(6,5) * t386 - Icges(6,6) * t387 + Icges(6,3) * t307;
t279 = Icges(6,4) * t387;
t210 = Icges(6,1) * t386 + Icges(6,5) * t307 - t279;
t208 = Icges(6,4) * t386 - Icges(6,2) * t387 + Icges(6,6) * t307;
t407 = t208 * t312;
t336 = t210 * t314 - t407;
t502 = -t206 * t307 - t336 * t308 - t192;
t322 = cos(qJ(4));
t433 = pkin(4) * t322;
t310 = pkin(3) + t433;
t393 = t307 * t312;
t357 = -rSges(6,2) * t393 - t308 * rSges(6,3);
t324 = -pkin(8) - pkin(7);
t383 = t308 * t324;
t427 = rSges(6,1) * t314;
t313 = sin(t319);
t435 = pkin(2) * t313;
t174 = t435 + t383 + (t310 + t427) * t307 + t357;
t432 = sin(qJ(1)) * pkin(1);
t168 = t174 + t432;
t501 = -t168 + t174;
t280 = t308 * t310;
t315 = cos(t319);
t309 = pkin(2) * t315;
t344 = rSges(6,1) * t386 - rSges(6,2) * t387;
t175 = t280 + t309 + (rSges(6,3) - t324) * t307 + t344;
t317 = cos(qJ(1)) * pkin(1);
t169 = t317 + t175;
t500 = -t169 + t175;
t268 = Icges(6,5) * t314 - Icges(6,6) * t312;
t400 = t268 * t307;
t205 = -Icges(6,3) * t308 + t400;
t499 = t205 * t307 + t209 * t386;
t305 = Icges(6,4) * t314;
t270 = -Icges(6,2) * t312 + t305;
t487 = Icges(6,1) * t312 + t305;
t498 = t270 + t487;
t316 = Icges(5,4) * t322;
t290 = -Icges(5,2) * t320 + t316;
t486 = Icges(5,1) * t320 + t316;
t497 = t290 + t486;
t301 = t308 * pkin(7);
t391 = t307 * t320;
t356 = -rSges(5,2) * t391 - t308 * rSges(5,3);
t428 = rSges(5,1) * t322;
t185 = t435 - t301 + (pkin(3) + t428) * t307 + t356;
t384 = t308 * t322;
t385 = t308 * t320;
t333 = rSges(5,1) * t384 - rSges(5,2) * t385 + rSges(5,3) * t307;
t353 = -t308 * pkin(3) - t307 * pkin(7);
t186 = t309 + t333 - t353;
t293 = rSges(5,1) * t320 + rSges(5,2) * t322;
t248 = t293 * t307;
t249 = t293 * t308;
t102 = -t185 * t248 - t186 * t249;
t490 = t350 * t308;
t91 = -t174 * t491 - t175 * t490;
t496 = m(5) * t102 + m(6) * t91;
t471 = m(5) / 0.2e1;
t470 = m(6) / 0.2e1;
t446 = -t307 / 0.2e1;
t495 = t307 / 0.2e1;
t445 = -t308 / 0.2e1;
t442 = m(3) * (-t317 * (rSges(3,1) * t313 + rSges(3,2) * t315) + t432 * (t315 * rSges(3,1) - rSges(3,2) * t313));
t440 = m(4) * (-t317 * (rSges(4,1) * t307 + rSges(4,2) * t308 + t435) + t432 * (t308 * rSges(4,1) - rSges(4,2) * t307 + t309));
t295 = -rSges(5,2) * t320 + t428;
t492 = t295 * t471;
t303 = t307 ^ 2;
t304 = t308 ^ 2;
t352 = t303 + t304;
t444 = t308 / 0.2e1;
t489 = t444 + t445;
t238 = t273 * t307;
t170 = -t307 * t238 - t239 * t308;
t275 = -rSges(6,2) * t312 + t427;
t396 = t275 * t308;
t397 = t275 * t307;
t338 = Icges(6,5) * t312 + Icges(6,6) * t314;
t232 = t307 * t338;
t233 = t338 * t308;
t269 = Icges(6,2) * t314 + t417;
t365 = -t269 * t307 + t209;
t207 = -Icges(6,6) * t308 + t270 * t307;
t367 = t307 * t487 + t207;
t481 = -t365 * t312 - t367 * t314;
t364 = -Icges(6,2) * t386 + t210 - t279;
t366 = t308 * t487 + t208;
t482 = -t364 * t312 - t366 * t314;
t431 = (-t304 * t232 + (t482 * t307 + (t233 - t481) * t308) * t307) * t445 + (t303 * t233 + (t481 * t308 + (-t232 - t482) * t307) * t308) * t446;
t195 = t307 * (rSges(6,1) * t392 + t357);
t211 = t307 * rSges(6,3) + t344;
t96 = t195 + (-t307 * t324 + t211 + t280 + t353) * t308 + (t383 + t301 + (-pkin(3) + t310) * t307) * t307;
t15 = t431 + m(6) * (t96 * t170 + t396 * t490 + t397 * t491);
t488 = t15 * qJD(5);
t75 = t168 * t175 - t174 * t169;
t181 = t185 + t432;
t182 = t186 + t317;
t88 = t181 * t186 - t185 * t182;
t485 = qJD(1) + qJD(2);
t430 = (t500 * t490 + t501 * t491) * t470 + ((-t182 + t186) * t308 + (-t181 + t185) * t307) * t293 * t471;
t90 = -t168 * t491 - t169 * t490;
t99 = -t181 * t248 - t182 * t249;
t484 = (t91 + t90) * t470 + (t102 + t99) * t471;
t483 = t498 * t312 + (t269 - t272) * t314;
t418 = Icges(5,4) * t320;
t289 = Icges(5,2) * t322 + t418;
t292 = Icges(5,1) * t322 - t418;
t480 = t497 * t320 + (t289 - t292) * t322;
t285 = Icges(5,4) * t385;
t226 = Icges(5,1) * t384 + Icges(5,5) * t307 - t285;
t360 = -Icges(5,2) * t384 + t226 - t285;
t224 = Icges(5,4) * t384 - Icges(5,2) * t385 + Icges(5,6) * t307;
t362 = t308 * t486 + t224;
t479 = -t360 * t320 - t362 * t322;
t225 = -Icges(5,5) * t308 + t292 * t307;
t361 = -t289 * t307 + t225;
t223 = -Icges(5,6) * t308 + t290 * t307;
t363 = t307 * t486 + t223;
t478 = -t361 * t320 - t363 * t322;
t477 = t497 * t322 / 0.2e1 + (t292 / 0.2e1 - t289 / 0.2e1) * t320;
t346 = t498 * t314 / 0.2e1 + (-t269 / 0.2e1 + t272 / 0.2e1) * t312;
t106 = -t205 * t308 - t207 * t393 + t192;
t193 = t210 * t392;
t107 = t206 * t308 + t208 * t393 - t193;
t406 = t209 * t314;
t408 = t207 * t312;
t348 = ((t193 + (t205 - t407) * t307 - t499) * t307 + ((t205 + t336) * t308 + (t406 + t408) * t307 + t502) * t308) * t446 + (-t106 * t308 - t107 * t307) * t495 + ((-t107 + (t206 - t406) * t308 + t499) * t308 + (t106 + (t206 + t408) * t307 + t502) * t307) * t445;
t476 = 4 * qJD(1);
t474 = 2 * qJD(4);
t467 = m(5) * t88;
t465 = m(5) * t99;
t92 = -t168 * t238 - t169 * t239;
t93 = -t174 * t238 - t175 * t239;
t459 = m(6) * (t93 + t92);
t457 = m(6) * (t501 * t307 + t500 * t308) * t273;
t149 = t168 * t396;
t369 = -t238 * t490 + t179;
t456 = m(6) * (-t169 * t397 + t149 + t369);
t409 = t490 * t273;
t413 = t169 * t275;
t455 = m(6) * (t149 - t179 + (t409 - t413) * t307);
t150 = t174 * t396;
t454 = m(6) * (-t175 * t397 + t150 + t369);
t411 = t175 * t275;
t453 = m(6) * (t150 - t179 + (t409 - t411) * t307);
t452 = m(6) * t75;
t450 = m(6) * t90;
t448 = m(6) * t92;
t447 = m(6) * t93;
t200 = t223 * t385;
t288 = Icges(5,5) * t322 - Icges(5,6) * t320;
t395 = t288 * t307;
t221 = -Icges(5,3) * t308 + t395;
t120 = -t221 * t307 - t225 * t384 + t200;
t222 = Icges(5,5) * t384 - Icges(5,6) * t385 + Icges(5,3) * t307;
t404 = t224 * t320;
t334 = t226 * t322 - t404;
t121 = t222 * t307 + t334 * t308;
t390 = t307 * t322;
t198 = t225 * t390;
t199 = t226 * t390;
t403 = t225 * t322;
t405 = t223 * t320;
t28 = (t120 + t199 - t200 + (t221 - t404) * t307) * t307 + (-t198 - t121 + (t221 + t334) * t308 + (t403 + t405) * t307) * t308;
t118 = -t221 * t308 - t223 * t391 + t198;
t119 = t222 * t308 + t224 * t391 - t199;
t29 = (-t119 + t200 + (t222 - t403) * t308) * t308 + (t118 - t198 + (t222 + t405) * t307) * t307;
t84 = -t118 * t308 - t119 * t307;
t85 = -t120 * t308 - t121 * t307;
t2 = (-t29 / 0.2e1 - t85 / 0.2e1) * t308 + (t84 / 0.2e1 - t28 / 0.2e1) * t307 + t348;
t443 = t2 * qJD(4);
t436 = m(6) * t170;
t422 = qJD(5) * t348;
t402 = t239 * t273;
t349 = t275 + t433;
t347 = t238 * t239;
t342 = t459 / 0.2e1 + t346;
t339 = Icges(5,5) * t320 + Icges(5,6) * t322;
t331 = -t348 + (t483 * t308 + t366 * t312 - t364 * t314 - t400) * t446 + (-t268 * t308 - t483 * t307 - t367 * t312 + t365 * t314) * t445;
t330 = -t346 + t489 * (t314 * t208 + t312 * t210);
t329 = t346 + t477;
t327 = t329 + t484;
t326 = t330 - t477 + t489 * (t224 * t322 + t226 * t320);
t325 = (t331 + t28 * t495 + (t480 * t308 + t362 * t320 - t360 * t322 - t395 + t84) * t446 + (-t288 * t308 - t480 * t307 - t363 * t320 + t361 * t322) * t445 + (t29 + t85) * t444) * qJD(4);
t243 = t339 * t308;
t242 = t307 * t339;
t220 = t349 * t308;
t218 = t349 * t307;
t176 = -t248 * t307 - t249 * t308;
t164 = qJD(5) * t436;
t152 = t308 * t211 + t195;
t148 = -t352 * t434 + t170;
t74 = t346 + t447;
t73 = t346 + t448;
t67 = t453 / 0.2e1;
t65 = t454 / 0.2e1;
t58 = t455 / 0.2e1;
t57 = t456 / 0.2e1;
t53 = t457 / 0.2e1;
t37 = t329 + t496;
t36 = t329 + t450 + t465;
t30 = t440 + t442 + t452 + t467;
t19 = -t457 / 0.2e1 + t342;
t18 = t53 + t342;
t17 = m(6) * (t352 * t273 * t275 + t152 * t170) + t431;
t16 = t17 * qJD(5);
t14 = t53 - t459 / 0.2e1 + t330;
t13 = t327 + t430;
t12 = t327 - t430;
t9 = t326 + t430 - t484;
t8 = t65 - t453 / 0.2e1 + t348;
t7 = t67 - t454 / 0.2e1 + t348;
t6 = t57 - t455 / 0.2e1 + t348;
t5 = t58 - t456 / 0.2e1 + t348;
t4 = t65 + t67 + t331;
t3 = t57 + t58 + t331;
t1 = [qJD(2) * t30 + qJD(4) * t36 + qJD(5) * t73, t30 * qJD(1) + t13 * qJD(4) + t18 * qJD(5) + 0.2e1 * (t440 / 0.2e1 + t442 / 0.2e1 + t75 * t470 + t88 * t471) * qJD(2), 0, t36 * qJD(1) + t13 * qJD(2) + t3 * qJD(5) + ((t181 * t308 - t182 * t307) * t492 + (t168 * t220 - t169 * t218) * t470) * t474 + t325, t73 * qJD(1) + t18 * qJD(2) + t3 * qJD(4) + ((t149 + (t402 - t413) * t307 - t347) * m(6) + t331) * qJD(5); t12 * qJD(4) + t19 * qJD(5) + (-t442 / 0.4e1 - t440 / 0.4e1 - t467 / 0.4e1 - t452 / 0.4e1) * t476, qJD(4) * t37 + qJD(5) * t74, 0, t12 * qJD(1) + t37 * qJD(2) + t4 * qJD(5) + ((t174 * t220 - t175 * t218) * t470 + (t185 * t308 - t186 * t307) * t492) * t474 + t325, t19 * qJD(1) + t74 * qJD(2) + t4 * qJD(4) + ((t150 + (t402 - t411) * t307 - t347) * m(6) + t331) * qJD(5); 0, 0, 0, (t148 * t470 + t176 * t471) * t474 + t164, qJD(4) * t436 + t164; t326 * qJD(1) + t9 * qJD(2) + t6 * qJD(5) + (-t465 / 0.4e1 - t450 / 0.4e1) * t476 + t443, t9 * qJD(1) + t8 * qJD(5) + t443 + (t326 - t496) * qJD(2), 0, (m(5) * (t352 * t295 * t293 + (t308 * t333 + t307 * (rSges(5,1) * t390 + t356)) * t176) + (-t242 * t304 + (t479 * t307 + (t243 - t478) * t308) * t307) * t445 + (t243 * t303 + (t478 * t308 + (-t242 - t479) * t307) * t308) * t446 + m(6) * (t148 * t96 + t218 * t491 + t220 * t490) + t431) * qJD(4) + t488 + t485 * t2, t6 * qJD(1) + t8 * qJD(2) + t15 * qJD(4) + t488; (t330 - t448) * qJD(1) + t14 * qJD(2) + t5 * qJD(4) + t422, t14 * qJD(1) + (t330 - t447) * qJD(2) + t7 * qJD(4) + t422, 0, t5 * qJD(1) + t7 * qJD(2) + ((t148 * t152 + (t218 * t307 + t220 * t308) * t273) * m(6) + t431) * qJD(4) + t16, qJD(4) * t17 + t348 * t485 + t16;];
Cq = t1;
