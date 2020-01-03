% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:20
% EndTime: 2019-12-31 17:42:28
% DurationCPUTime: 6.40s
% Computational Cost: add. (39067->419), mult. (36829->628), div. (0->0), fcn. (38965->10), ass. (0->283)
t302 = sin(pkin(8));
t299 = t302 ^ 2;
t303 = cos(pkin(8));
t300 = t303 ^ 2;
t488 = t299 + t300;
t301 = qJ(2) + qJ(3);
t294 = pkin(9) + t301;
t291 = sin(t294);
t292 = cos(t294);
t295 = sin(t301);
t296 = cos(t301);
t493 = -Icges(4,5) * t295 - Icges(5,5) * t291 - Icges(4,6) * t296 - Icges(5,6) * t292;
t304 = sin(qJ(5));
t306 = cos(qJ(5));
t353 = rSges(6,1) * t306 - rSges(6,2) * t304;
t220 = -rSges(6,3) * t292 + t291 * t353;
t200 = t220 * t302;
t201 = t220 * t303;
t284 = rSges(4,1) * t295 + rSges(4,2) * t296;
t180 = t488 * t284;
t305 = sin(qJ(2));
t307 = cos(qJ(2));
t486 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t307 * t305 + (-0.2e1 * t305 ^ 2 + 0.2e1 * t307 ^ 2) * Icges(3,4);
t451 = t302 / 0.2e1;
t450 = -t303 / 0.2e1;
t337 = -Icges(3,5) * t305 - Icges(3,6) * t307;
t278 = t337 * t302;
t279 = t337 * t303;
t485 = t493 * t302;
t484 = t493 * t303;
t477 = t493 * t488;
t379 = qJD(2) + qJD(3);
t414 = (-Icges(6,5) * t304 - Icges(6,6) * t306) * t291 * t292;
t423 = Icges(6,4) * t306;
t338 = -Icges(6,2) * t304 + t423;
t212 = -Icges(6,6) * t292 + t291 * t338;
t391 = -t212 + (-Icges(6,1) * t304 - t423) * t291;
t424 = Icges(6,4) * t304;
t345 = Icges(6,1) * t306 - t424;
t214 = -Icges(6,5) * t292 + t291 * t345;
t390 = t214 + (-Icges(6,2) * t306 - t424) * t291;
t276 = pkin(4) * t291 - pkin(7) * t292;
t475 = -t220 - t276;
t221 = rSges(6,3) * t291 + t292 * t353;
t277 = pkin(4) * t292 + pkin(7) * t291;
t474 = -t221 - t277;
t285 = rSges(4,1) * t296 - rSges(4,2) * t295;
t170 = t488 * t285;
t297 = t307 * pkin(2);
t389 = t488 * t297;
t142 = t170 + t389;
t445 = pkin(2) * t305;
t369 = -t284 - t445;
t238 = t369 * t302;
t240 = t369 * t303;
t106 = -t142 * t180 + (-t238 * t302 - t240 * t303) * t285;
t444 = pkin(3) * t295;
t359 = t488 * t444;
t360 = -t302 * t200 - t303 * t201 - t276 * t488;
t131 = -t359 + t360;
t287 = -t444 - t445;
t323 = t287 + t475;
t155 = t323 * t302;
t157 = t323 * t303;
t443 = pkin(3) * t296;
t356 = -t443 + t474;
t162 = t356 * t302;
t164 = t356 * t303;
t400 = t155 * t162 + t157 * t164;
t404 = t303 * t306;
t407 = t302 * t304;
t270 = -t292 * t407 - t404;
t405 = t303 * t304;
t406 = t302 * t306;
t271 = t292 * t406 - t405;
t411 = t291 * t302;
t177 = rSges(6,1) * t271 + rSges(6,2) * t270 + rSges(6,3) * t411;
t272 = -t292 * t405 + t406;
t273 = t292 * t404 + t407;
t410 = t291 * t303;
t178 = rSges(6,1) * t273 + rSges(6,2) * t272 + rSges(6,3) * t410;
t392 = t488 * t443;
t103 = t302 * t177 + t303 * t178 + t277 * t488 + t392;
t94 = t103 + t389;
t36 = t94 * t131 + t400;
t275 = rSges(5,1) * t292 - rSges(5,2) * t291;
t132 = t275 * t488 + t392;
t117 = t132 + t389;
t274 = rSges(5,1) * t291 + rSges(5,2) * t292;
t385 = t488 * t274;
t154 = -t359 - t385;
t328 = -t274 + t287;
t202 = t328 * t302;
t204 = t328 * t303;
t365 = -t275 - t443;
t217 = t365 * t302;
t219 = t365 * t303;
t398 = t202 * t217 + t204 * t219;
t63 = t117 * t154 + t398;
t472 = m(4) * t106 + m(5) * t63 + m(6) * t36;
t181 = Icges(6,5) * t270 - Icges(6,6) * t271;
t258 = Icges(6,4) * t270;
t175 = Icges(6,1) * t271 + Icges(6,5) * t411 + t258;
t394 = -Icges(6,2) * t271 + t175 + t258;
t426 = Icges(6,4) * t271;
t173 = Icges(6,2) * t270 + Icges(6,6) * t411 + t426;
t396 = Icges(6,1) * t270 - t173 - t426;
t84 = t181 * t411 + t270 * t394 + t271 * t396;
t182 = Icges(6,5) * t272 - Icges(6,6) * t273;
t259 = Icges(6,4) * t272;
t176 = Icges(6,1) * t273 + Icges(6,5) * t410 + t259;
t393 = -Icges(6,2) * t273 + t176 + t259;
t425 = Icges(6,4) * t273;
t174 = Icges(6,2) * t272 + Icges(6,6) * t410 + t425;
t395 = Icges(6,1) * t272 - t174 - t425;
t85 = t182 * t411 + t270 * t393 + t271 * t395;
t45 = t302 * t85 - t303 * t84;
t86 = t181 * t410 + t272 * t394 + t273 * t396;
t87 = t182 * t410 + t272 * t393 + t273 * t395;
t46 = t302 * t87 - t303 * t86;
t440 = t45 * t450 + t46 * t451;
t334 = Icges(6,5) * t306 - Icges(6,6) * t304;
t210 = -Icges(6,3) * t292 + t291 * t334;
t466 = 2 * qJD(2);
t465 = 4 * qJD(3);
t464 = m(4) / 0.2e1;
t463 = m(5) / 0.2e1;
t462 = m(6) / 0.2e1;
t386 = t488 * (t287 + t445);
t145 = -t385 + t386;
t366 = -t274 - t444;
t216 = t366 * t302;
t218 = t366 * t303;
t397 = t216 * t217 + t218 * t219;
t460 = m(5) * (t132 * t145 + t397);
t330 = t177 * t303 - t178 * t302;
t108 = t330 * t292 + (-t200 * t303 + t201 * t302) * t291;
t123 = (t221 * t302 - t177) * t291;
t124 = (-t221 * t303 + t178) * t291;
t374 = t108 * t94 + t123 * t157 + t124 * t155;
t147 = t177 * t292 + t220 * t411;
t148 = -t178 * t292 - t220 * t410;
t401 = t147 * t164 + t148 * t162;
t137 = t330 * t291;
t99 = t137 * t131;
t459 = m(6) * (t99 + t374 + t401);
t122 = t360 + t386;
t357 = -t444 + t475;
t161 = t357 * t302;
t163 = t357 * t303;
t354 = t108 * t103 + t123 * t163 + t124 * t161 + t401;
t458 = m(6) * (t122 * t137 + t354);
t257 = (-rSges(6,1) * t304 - rSges(6,2) * t306) * t291;
t187 = rSges(6,1) * t270 - rSges(6,2) * t271;
t188 = rSges(6,1) * t272 - rSges(6,2) * t273;
t150 = t187 * t302 + t188 * t303;
t65 = t94 * t150;
t79 = t103 * t150;
t457 = m(6) * (t65 + t79 + ((-t157 - t163) * t303 + (-t155 - t161) * t302) * t257);
t456 = m(6) * (t108 * t137 + t123 * t147 + t124 * t148);
t399 = t161 * t162 + t163 * t164;
t454 = m(6) * (t103 * t122 + t399);
t453 = m(6) * (t123 * t302 - t124 * t303);
t452 = -t292 / 0.2e1;
t449 = t303 / 0.2e1;
t126 = (-t170 + t285) * t180;
t125 = m(4) * t126;
t447 = m(5) * t154;
t446 = m(6) * t131;
t376 = -t453 / 0.2e1;
t92 = 0.2e1 * (t108 / 0.4e1 - t150 / 0.4e1) * m(6);
t402 = t92 * qJD(1);
t439 = qJD(4) * t376 - t402;
t375 = t453 / 0.2e1;
t438 = t379 * t375;
t437 = m(6) * qJD(5);
t116 = m(6) * (-t162 * t303 + t164 * t302) + m(5) * (-t217 * t303 + t219 * t302);
t433 = t116 * qJD(3) + qJD(5) * t375;
t172 = Icges(6,5) * t273 + Icges(6,6) * t272 + Icges(6,3) * t410;
t112 = t172 * t411 + t174 * t270 + t176 * t271;
t419 = t112 * t303;
t171 = Icges(6,5) * t271 + Icges(6,6) * t270 + Icges(6,3) * t411;
t113 = t171 * t410 + t173 * t272 + t175 * t273;
t418 = t113 * t302;
t417 = t171 * t292;
t416 = t172 * t292;
t415 = t210 * t292;
t409 = t292 * t302;
t408 = t292 * t303;
t54 = 0.2e1 * (t131 / 0.4e1 - t122 / 0.4e1) * m(6) + 0.2e1 * (t154 / 0.4e1 - t145 / 0.4e1) * m(5);
t403 = t54 * qJD(1);
t111 = t171 * t411 + t173 * t270 + t175 * t271;
t133 = t210 * t411 + t212 * t270 + t214 * t271;
t213 = Icges(6,6) * t291 + t292 * t338;
t215 = Icges(6,5) * t291 + t292 * t345;
t329 = -t212 * t304 + t214 * t306;
t320 = (Icges(6,3) * t291 + t292 * t334 - t329) * t292;
t195 = t212 * t302;
t197 = t214 * t302;
t332 = -t173 * t304 + t175 * t306;
t325 = -t210 * t302 - t332;
t310 = t291 * t325 + t417;
t75 = -t195 * t270 - t197 * t271 + t302 * t310;
t196 = t212 * t303;
t198 = t214 * t303;
t331 = -t174 * t304 + t176 * t306;
t324 = -t210 * t303 - t331;
t309 = t291 * t324 + t416;
t76 = -t196 * t270 - t198 * t271 + t302 * t309;
t18 = (t419 - t213 * t270 - t215 * t271 + (t111 - t415) * t302) * t292 + (t76 * t303 + t133 + (t75 - t320) * t302) * t291;
t114 = t172 * t410 + t174 * t272 + t176 * t273;
t134 = t210 * t410 + t212 * t272 + t214 * t273;
t77 = -t195 * t272 - t197 * t273 + t303 * t310;
t78 = -t196 * t272 - t198 * t273 + t303 * t309;
t19 = (t418 - t213 * t272 - t215 * t273 + (t114 - t415) * t303) * t292 + (t77 * t302 + t134 + (t78 - t320) * t303) * t291;
t141 = t291 * t329 - t415;
t120 = t291 * t332 - t417;
t121 = t291 * t331 - t416;
t333 = t120 * t302 + t121 * t303;
t82 = -t325 * t292 + (t195 * t304 - t197 * t306 + t171) * t291;
t83 = -t324 * t292 + (t196 * t304 - t198 * t306 + t172) * t291;
t21 = (t320 + t333) * t292 + (t83 * t303 + t82 * t302 - (-t213 * t304 + t215 * t306 + t210) * t292 + t141) * t291;
t50 = -t133 * t292 + (t111 * t302 + t419) * t291;
t51 = -t134 * t292 + (t114 * t303 + t418) * t291;
t58 = -t141 * t292 + t291 * t333;
t3 = t456 + (t51 * t449 + t50 * t451 - t21 / 0.2e1) * t292 + (t19 * t449 + t18 * t451 + t58 / 0.2e1) * t291;
t378 = qJD(4) * t375 + t3 * qJD(5) + t402;
t377 = t457 / 0.2e1 + t440;
t373 = t103 * t131 + t399;
t372 = t132 * t154 + t397;
t371 = t411 / 0.2e1;
t370 = t410 / 0.2e1;
t368 = -t285 - t297;
t358 = t488 * t445;
t355 = -t297 - t443;
t288 = rSges(3,1) * t305 + rSges(3,2) * t307;
t34 = t302 * t76 - t303 * t75;
t35 = t302 * t78 - t303 * t77;
t352 = (t35 + (-t485 * t302 + t477) * t303 + t484 * t299) * t451 + (t34 + (-t484 * t303 + t477) * t302 + t485 * t300) * t450;
t327 = -t275 + t355;
t326 = t125 + t352;
t322 = t355 + t474;
t321 = -t358 + t386;
t319 = t137 * t150 + (-t147 * t303 - t148 * t302) * t257;
t318 = t486 * t302 + t279;
t317 = -t486 * t303 + t278;
t316 = t18 * t450 + t19 * t451 + t34 * t371 + t35 * t370 + (t302 * t83 - t303 * t82) * t452 + (-t111 * t303 + t112 * t302) * t409 / 0.2e1 + (-t113 * t303 + t114 * t302) * t408 / 0.2e1 + t291 * (-t120 * t303 + t121 * t302) / 0.2e1 - t440;
t28 = -(t270 * t390 + t271 * t391) * t292 + (t85 * t303 + (t84 - t414) * t302) * t291;
t29 = -(t272 * t390 + t273 * t391) * t292 + (t86 * t302 + (t87 - t414) * t303) * t291;
t97 = -t181 * t292 + (-t304 * t394 + t306 * t396) * t291;
t98 = -t182 * t292 + (-t304 * t393 + t306 * t395) * t291;
t311 = -t18 * t411 / 0.2e1 - t19 * t410 / 0.2e1 + t292 * t21 / 0.2e1 - t456 + t29 * t451 + t28 * t450 + t45 * t371 + t46 * t370 - t50 * t409 / 0.2e1 - t51 * t408 / 0.2e1 + (t302 * t98 - t303 * t97) * t452 - t291 * t58 / 0.2e1;
t241 = t368 * t303;
t239 = t368 * t302;
t205 = t327 * t303;
t203 = t327 * t302;
t199 = t488 * t288;
t179 = m(4) * t180;
t165 = -t358 - t180;
t158 = t322 * t303;
t156 = t322 * t302;
t152 = -t188 * t292 - t257 * t410;
t151 = t187 * t292 + t257 * t411;
t144 = (t187 * t303 - t188 * t302) * t291;
t140 = t321 - t385;
t118 = t321 + t360;
t93 = (t108 + t150) * t462;
t90 = t92 * qJD(5);
t89 = t93 * qJD(5);
t68 = qJD(5) * t376;
t57 = t79 + (-t161 * t302 - t163 * t303) * t257;
t52 = -t179 + t145 * t463 + t122 * t462 + t447 / 0.2e1 + t446 / 0.2e1;
t49 = t65 + (-t155 * t302 - t157 * t303) * t257;
t14 = t458 / 0.2e1;
t10 = t459 / 0.2e1;
t9 = m(6) * t57 + t440;
t8 = m(6) * t49 + t440;
t7 = t326 + t454 + t460;
t6 = t352 + t472;
t5 = t14 - t459 / 0.2e1 + t377;
t4 = t10 - t458 / 0.2e1 + t377;
t1 = t10 + t14 - t457 / 0.2e1 + t316;
t2 = [0, t52 * qJD(3) + t89 + (-m(3) * t199 / 0.2e1 + t165 * t464 + t140 * t463 + t118 * t462) * t466, t52 * qJD(2) + (-t179 + t446 + t447) * qJD(3) + t89, 0, t144 * t437 + t379 * t93; qJD(3) * t54 - t90, (m(6) * (t118 * t94 + t155 * t156 + t157 * t158) + m(5) * (t117 * t140 + t202 * t203 + t204 * t205) + m(4) * (t142 * t165 + t238 * t239 + t240 * t241) + (t299 * t279 + (t318 * t303 + (-t278 + t317) * t302) * t303) * t451 + (t300 * t278 + (t317 * t302 + (-t279 + t318) * t303) * t302) * t450 + t352 + m(3) * (-t199 + t288) * t488 * (rSges(3,1) * t307 - rSges(3,2) * t305)) * qJD(2) + t6 * qJD(3) + t8 * qJD(5), t403 + t6 * qJD(2) + t352 * qJD(3) + t4 * qJD(5) + (-t454 / 0.4e1 - t460 / 0.4e1 - t125 / 0.4e1) * t465 + 0.2e1 * ((t373 + t36) * t462 + (t372 + t63) * t463 + (t106 + t126) * t464) * qJD(3), t68, t8 * qJD(2) + t4 * qJD(3) + (t311 + m(6) * (t144 * t94 + t151 * t157 + t152 * t155 + t319)) * qJD(5) + t439; -qJD(2) * t54 - t90, -t403 + t7 * qJD(3) + t5 * qJD(5) + ((t103 * t118 + t122 * t94 + t156 * t161 + t158 * t163 + t400) * t462 + (t117 * t145 + t132 * t140 + t203 * t216 + t205 * t218 + t398) * t463 + (t165 * t170 + (-t239 * t302 - t241 * t303) * t284 + t106) * t464) * t466 + (t352 - t472) * qJD(2), t7 * qJD(2) + t326 * qJD(3) + t9 * qJD(5) + (m(6) * t373 / 0.4e1 + m(5) * t372 / 0.4e1) * t465, t68, t5 * qJD(2) + t9 * qJD(3) + (t311 + m(6) * (t103 * t144 + t151 * t163 + t152 * t161 + t319)) * qJD(5) + t439; 0, ((-t156 * t303 + t158 * t302) * t462 + (-t203 * t303 + t205 * t302) * t463) * t466 + t433, qJD(2) * t116 + t433, 0, (t151 * t302 - t152 * t303) * t437 + t438; t379 * t92, ((t118 * t137 + t147 * t158 + t148 * t156 + t374 - t49) * m(6) + t316) * qJD(2) + t1 * qJD(3) + t378, t1 * qJD(2) + ((t99 - t57 + t354) * m(6) + t316) * qJD(3) + t378, t438, t379 * t3 + (m(6) * (t137 * t144 + t147 * t151 + t148 * t152) - t292 ^ 2 * t414 / 0.2e1 + (t29 * t449 + t28 * t451 + (t98 * t303 + t97 * t302 - (-t390 * t304 + t391 * t306) * t292) * t452) * t291) * qJD(5);];
Cq = t2;
