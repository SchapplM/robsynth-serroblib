% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:05
% EndTime: 2019-12-31 17:49:21
% DurationCPUTime: 12.90s
% Computational Cost: add. (9237->429), mult. (8113->528), div. (0->0), fcn. (6153->8), ass. (0->252)
t225 = pkin(8) + qJ(4);
t222 = cos(t225);
t220 = sin(t225);
t207 = Icges(6,5) * t220;
t277 = Icges(6,3) * t222 - t207;
t376 = Icges(5,4) * t220;
t481 = Icges(5,2) * t222 + t277 + t376;
t373 = Icges(6,5) * t222;
t156 = Icges(6,1) * t220 - t373;
t210 = Icges(5,4) * t222;
t480 = Icges(5,1) * t220 + t156 + t210;
t151 = Icges(5,5) * t222 - Icges(5,6) * t220;
t226 = qJ(1) + pkin(7);
t221 = sin(t226);
t223 = cos(t226);
t100 = Icges(5,3) * t221 + t151 * t223;
t153 = Icges(6,4) * t222 + Icges(6,6) * t220;
t102 = Icges(6,2) * t221 + t153 * t223;
t483 = t100 + t102;
t279 = Icges(6,1) * t222 + t207;
t106 = Icges(6,4) * t221 + t223 * t279;
t159 = Icges(5,1) * t222 - t376;
t108 = Icges(5,5) * t221 + t159 * t223;
t457 = t106 + t108;
t149 = Icges(6,3) * t220 + t373;
t278 = -Icges(5,2) * t220 + t210;
t482 = -t149 + t278;
t475 = (Icges(5,6) - Icges(6,6)) * t222 + (Icges(6,4) + Icges(5,5)) * t220;
t479 = t159 + t279;
t470 = -t481 * t220 + t480 * t222;
t356 = t222 * t223;
t358 = t220 * t223;
t177 = Icges(6,5) * t356;
t369 = Icges(6,6) * t221;
t98 = Icges(6,3) * t358 + t177 + t369;
t478 = t483 * t221 + t457 * t356 + t98 * t358;
t477 = t482 * qJD(4);
t476 = t479 * qJD(4);
t474 = -t151 - t153;
t472 = t481 * qJD(4);
t471 = t480 * qJD(4);
t413 = t475 * t223;
t414 = t475 * t221;
t359 = t220 * t221;
t178 = Icges(5,4) * t359;
t357 = t221 * t222;
t374 = Icges(5,5) * t223;
t107 = Icges(5,1) * t357 - t178 - t374;
t370 = Icges(5,6) * t223;
t103 = Icges(5,4) * t357 - Icges(5,2) * t359 - t370;
t365 = t103 * t220;
t276 = -t107 * t222 + t365;
t101 = -Icges(6,2) * t223 + t153 * t221;
t355 = t223 * t101;
t105 = -Icges(6,4) * t223 + t221 * t279;
t97 = -Icges(6,6) * t223 + t149 * t221;
t282 = t105 * t222 + t220 * t97;
t416 = t221 * t282;
t34 = -t355 + t416;
t367 = Icges(5,3) * t223;
t99 = Icges(5,5) * t357 - Icges(5,6) * t359 - t367;
t469 = -t221 * t276 - t223 * t99 + t34;
t104 = Icges(5,6) * t221 + t223 * t278;
t439 = -t104 * t358 + t478;
t468 = t470 * t221 - t413;
t467 = t470 * t223 + t414;
t290 = t102 * t223 - t106 * t357 - t98 * t359;
t82 = t108 * t357;
t297 = t223 * t100 - t82;
t37 = -t104 * t359 - t297;
t466 = -t290 + t37;
t229 = -pkin(6) - qJ(3);
t194 = t223 * t229;
t228 = cos(pkin(8));
t219 = pkin(3) * t228 + pkin(2);
t334 = -t221 * t219 - t194;
t230 = sin(qJ(1));
t394 = pkin(1) * t230;
t465 = t334 - t394;
t464 = t472 * t223 + (t221 * t278 - t370 - t97) * qJD(1);
t463 = t472 * t221 + (t149 * t223 - t104 + t369) * qJD(1);
t462 = -t471 * t223 + (-t159 * t221 - t105 + t374) * qJD(1);
t461 = -t457 * qJD(1) + t471 * t221;
t460 = t103 - t97;
t459 = t104 - t98;
t458 = t105 + t107;
t456 = t476 * t222 - t477 * t220 + (-t220 * t480 - t222 * t481) * qJD(4) + t475 * qJD(1);
t315 = rSges(5,1) * t357;
t455 = -t315 + t465;
t364 = t104 * t220;
t454 = t220 * t98 + t457 * t222 - t364;
t453 = t276 - t282;
t452 = t470 * qJD(1) + t474 * qJD(4);
t451 = t467 * qJD(1);
t388 = -t107 * t356 - t221 * t99;
t40 = -t103 * t358 - t388;
t92 = t221 * t101;
t38 = t105 * t356 + t97 * t358 + t92;
t427 = t223 * t38;
t449 = (t439 * t221 - t223 * t40 - t427) * qJD(4);
t448 = (t466 * t221 - t469 * t223) * qJD(4);
t447 = t468 * qJD(1);
t446 = t447 + t448;
t445 = t449 + t451;
t444 = t453 * qJD(4) + t461 * t220 + t463 * t222;
t443 = t454 * qJD(4) + t462 * t220 - t464 * t222;
t442 = -t452 * t221 + t456 * t223;
t441 = t456 * t221 + t452 * t223;
t440 = t38 + t40;
t438 = t458 * t220 + t460 * t222;
t437 = t457 * t220 + t459 * t222;
t436 = t475 * qJD(4);
t435 = t355 + t478;
t434 = t221 * (-t223 * t481 + t457) - t223 * (-Icges(5,2) * t357 - t277 * t221 - t178 + t458);
t426 = rSges(6,3) + qJ(5);
t433 = t460 * t223 + (-Icges(6,1) * t358 + t156 * t223 + t177 - t459) * t221;
t432 = t483 * qJD(1);
t431 = t480 + t482;
t430 = -t481 + t479;
t232 = qJD(1) ^ 2;
t395 = rSges(6,1) + pkin(4);
t160 = pkin(4) * t220 - qJ(5) * t222;
t161 = rSges(6,1) * t220 - rSges(6,3) * t222;
t339 = t160 + t161;
t166 = rSges(6,1) * t222 + rSges(6,3) * t220;
t425 = pkin(4) * t222 + qJ(5) * t220 + t166;
t424 = -t437 * qJD(4) + t464 * t220 + t462 * t222 + t432;
t406 = qJD(1) * t101;
t423 = -qJD(1) * t99 + t438 * qJD(4) - t463 * t220 + t461 * t222 - t406;
t384 = rSges(4,2) * sin(pkin(8));
t386 = rSges(4,1) * t228;
t269 = t221 * rSges(4,3) + (-t384 + t386) * t223;
t201 = t221 * qJ(3);
t168 = t223 * pkin(2) + t201;
t231 = cos(qJ(1));
t224 = t231 * pkin(1);
t300 = t168 + t224;
t422 = t269 + t300;
t421 = t453 * qJD(1) - t436 * t221 + t432;
t420 = -t406 - t436 * t223 + (-t151 * t221 + t367 - t454) * qJD(1);
t419 = 0.2e1 * qJD(4);
t202 = t223 * qJ(3);
t163 = pkin(2) * t221 - t202;
t147 = qJD(1) * t163;
t95 = t163 + t334;
t418 = qJD(1) * t95 - t147;
t320 = qJD(5) * t222;
t213 = t221 * rSges(6,2);
t347 = t356 * t395 + t358 * t426 + t213;
t216 = t223 * rSges(6,2);
t348 = t221 * t425 - t216;
t31 = -t320 + qJD(2) + (t221 * t348 + t223 * t347) * qJD(4);
t417 = qJD(4) * t31;
t294 = t339 * qJD(4);
t321 = qJD(5) * t220;
t254 = -t294 + t321;
t415 = t254 * t221;
t322 = qJD(4) * t223;
t310 = t222 * t322;
t324 = qJD(1) * t223;
t412 = rSges(6,2) * t324 + t310 * t426;
t411 = -rSges(5,2) * t359 - t223 * rSges(5,3);
t298 = t223 * rSges(3,1) - rSges(3,2) * t221;
t410 = t224 + t298;
t409 = -t220 * t434 + t433 * t222;
t408 = (-t220 * t431 + t222 * t430) * qJD(1);
t407 = t474 * qJD(1);
t130 = t161 * t221;
t323 = qJD(4) * t221;
t389 = (pkin(4) * t324 + qJ(5) * t323) * t222 + (qJ(5) * t324 + (-pkin(4) * qJD(4) + qJD(5)) * t221) * t220 - qJD(4) * t130 + (t166 * t223 + t213) * qJD(1);
t173 = t223 * t321;
t311 = t220 * t322;
t325 = qJD(1) * t221;
t251 = -t222 * t325 - t311;
t313 = t220 * t325;
t390 = t251 * t395 - t313 * t426 + t173 + t412;
t1 = (t321 + t390 * t223 + t389 * t221 + (-t221 * t347 + t223 * t348) * qJD(1)) * qJD(4);
t398 = m(6) * t1;
t397 = t221 / 0.2e1;
t396 = -t223 / 0.2e1;
t392 = qJD(1) / 0.2e1;
t391 = pkin(2) - t219;
t162 = rSges(5,1) * t220 + rSges(5,2) * t222;
t135 = t162 * t223;
t211 = t221 * rSges(5,3);
t112 = rSges(5,1) * t356 - rSges(5,2) * t358 + t211;
t200 = qJD(3) * t223;
t180 = t223 * t219;
t295 = -t221 * t229 + t180;
t293 = t300 + t295 - t168;
t312 = t162 * t323;
t43 = -t312 - t200 + (t112 + t293) * qJD(1);
t383 = t135 * t43;
t110 = t315 + t411;
t199 = qJD(3) * t221;
t309 = t162 * t322;
t287 = t199 - t309;
t301 = -t163 - t394;
t42 = (-t110 + t301 + t95) * qJD(1) + t287;
t382 = t221 * t42;
t128 = qJD(1) * t168 - t200;
t187 = t229 * t325;
t379 = -t128 + t187 - (-t223 * t391 - t201) * qJD(1);
t346 = -qJD(4) * t425 + t320;
t345 = -t160 * t221 - t130;
t344 = t339 * t223;
t336 = rSges(5,2) * t313 + rSges(5,3) * t324;
t335 = t173 + t199;
t188 = t221 * t384;
t333 = rSges(4,3) * t324 + qJD(1) * t188;
t332 = t187 + t200;
t331 = t223 * rSges(4,3) + t188;
t195 = qJ(3) * t324;
t330 = t195 + t199;
t319 = qJD(1) * qJD(3);
t318 = t232 * t394;
t317 = t232 * t224;
t316 = t221 * t386;
t308 = -pkin(2) - t386;
t305 = -t323 / 0.2e1;
t302 = t322 / 0.2e1;
t296 = -t99 + t364;
t292 = t223 * t319 - t317;
t291 = -t348 - t394;
t164 = rSges(3,1) * t221 + rSges(3,2) * t223;
t284 = rSges(5,1) * t222 - rSges(5,2) * t220;
t283 = -t221 * t43 - t223 * t42;
t274 = t110 * t221 + t112 * t223;
t270 = qJD(1) * (-pkin(2) * t325 + t330) + t221 * t319 - t318;
t267 = t270 + qJD(1) * (-t195 + (t221 * t391 - t194) * qJD(1));
t266 = qJD(4) * (t320 + t346);
t131 = t162 * t221;
t252 = -t223 * t294 + t335;
t246 = -t219 - t425;
t74 = rSges(5,1) * t251 - rSges(5,2) * t310 + t336;
t76 = -qJD(4) * t131 + (t223 * t284 + t211) * qJD(1);
t245 = t221 * t76 + t223 * t74 + (t110 * t223 - t112 * t221) * qJD(1);
t146 = t284 * qJD(4);
t113 = t316 - t331;
t78 = qJD(1) * t422 - t200;
t77 = t199 + (-t113 + t301) * qJD(1);
t50 = (-qJD(1) * t269 - t128) * qJD(1) + t292;
t49 = qJD(1) * (-qJD(1) * t316 + t333) + t270;
t44 = qJD(4) * t274 + qJD(2);
t30 = -t200 + t415 + (t293 + t347) * qJD(1);
t29 = (-t163 + t291 + t95) * qJD(1) + t252;
t26 = -t146 * t322 + (-t76 + t312 + t379) * qJD(1) + t292;
t25 = -t146 * t323 + (t74 - t309) * qJD(1) + t267;
t16 = t245 * qJD(4);
t11 = t223 * t266 + (t379 - t389 - t415) * qJD(1) + t292;
t10 = t221 * t266 + (t223 * t254 + t390) * qJD(1) + t267;
t2 = [m(3) * ((-t164 * t232 - t318) * t410 + (-t317 + (-0.2e1 * t298 - t224 + t410) * t232) * (-t164 - t394)) + (((t37 - t82 + (t100 + t365) * t223 + t388) * t223 - t427 + (t34 - t416 + t435) * t221) * qJD(4) + t451) * t302 + (t470 * qJD(4) + t476 * t220 + t477 * t222) * qJD(1) + (t11 * (t216 + t465) + t29 * t332 + t10 * (t180 + t224 + t347) + t30 * (-t395 * t311 + t335 + t412) + (-t10 * t229 - t11 * t395 * t222 + (-t29 * qJD(5) - t11 * t426) * t220 + t29 * (t220 * t395 - t222 * t426) * qJD(4)) * t221 + ((-t230 * t30 - t231 * t29) * pkin(1) + (-t30 * t229 + t246 * t29) * t223 + (-t29 * rSges(6,2) + t246 * t30) * t221) * qJD(1) - (qJD(1) * t291 + t252 - t29 + t418) * t30) * m(6) + (t26 * (-t411 + t455) + t25 * (t112 + t224 + t295) + (t162 * t382 - t383) * qJD(4) + (t332 + (-t211 - t224 + (-t219 - t284) * t223) * qJD(1)) * t42 + (t42 - t287 - t418 + t199 + t336 + (t110 + t394 + t455) * qJD(1)) * t43) * m(5) + (t50 * (t221 * t308 + t202 + t331 - t394) + t77 * t200 + t49 * t422 + t78 * (t330 + t333) + ((-t230 * t78 - t231 * t77) * pkin(1) + t77 * (t308 + t384) * t223 + (t77 * (-rSges(4,3) - qJ(3)) + t78 * t308) * t221) * qJD(1) - (-t147 + t199 - t77 + (-t113 - t394) * qJD(1)) * t78) * m(4) + (((t223 * t296 - t435 + t439) * t223 + (t221 * t296 + t290 + t297 + t440 - t92) * t221) * qJD(4) + t446 - t447) * t305 + (t442 + t443) * t323 / 0.2e1 - (t441 - t444 + t445) * t322 / 0.2e1 + ((t438 + t468) * t221 + (t437 + t467) * t223) * qJD(4) * t392; m(5) * t16 + t398; 0.2e1 * (t10 * t396 + t11 * t397) * m(6) + 0.2e1 * (t25 * t396 + t26 * t397) * m(5) + 0.2e1 * (t396 * t49 + t397 * t50) * m(4); -((t433 * t220 + t222 * t434) * qJD(4) + (t430 * t220 + t431 * t222) * qJD(1)) * qJD(1) / 0.2e1 + (t444 * t223 + t443 * t221 + (t221 * t438 + t223 * t437) * qJD(1)) * t392 + ((-t323 * t413 - t407) * t221 + ((t221 * t414 + t409) * qJD(4) + t408) * t223) * t305 + ((-t322 * t414 + t407) * t223 + ((t223 * t413 + t409) * qJD(4) + t408) * t221) * t302 + ((-t11 * t339 + t29 * t346 + t1 * t347 + t31 * t390 + (-t30 * t339 + t31 * t348) * qJD(1)) * t223 + (-t10 * t339 + t30 * t346 + t1 * t348 + t31 * t389 + (t29 * t339 - t31 * t347) * qJD(1)) * t221 - (t220 * t31 + (t221 * t30 + t223 * t29) * t222) * qJD(5) - (-t29 * t345 - t30 * t344) * qJD(1) - ((-t29 * t425 - t31 * t344) * t223 + (-t30 * t425 + t31 * t345) * t221) * qJD(4)) * m(6) + (t16 * t274 + t44 * t245 + t283 * t146 + (-t25 * t221 - t26 * t223 + (-t223 * t43 + t382) * qJD(1)) * t162 - (t131 * t42 - t383) * qJD(1) - (t44 * (-t131 * t221 - t135 * t223) + t283 * t284) * qJD(4)) * m(5) + (t442 * qJD(1) + ((t439 * qJD(1) + t423 * t223) * t223 + (t420 * t221 + t440 * qJD(1) + (-t421 + t424) * t223) * t221) * t419) * t397 + (t441 * qJD(1) + ((t466 * qJD(1) + t421 * t223) * t223 + (t424 * t221 + t469 * qJD(1) + (-t420 + t423) * t223) * t221) * t419) * t396 + (t446 + t448) * t325 / 0.2e1 + (t445 + t449) * t324 / 0.2e1; -t222 * t398 + 0.2e1 * (m(6) * (t10 * t221 + t11 * t223 + t417) / 0.2e1 - m(6) * (t221 ^ 2 + t223 ^ 2) * t417 / 0.2e1) * t220;];
tauc = t2(:);
