% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:47
% EndTime: 2019-12-31 17:40:55
% DurationCPUTime: 5.12s
% Computational Cost: add. (9305->318), mult. (10010->433), div. (0->0), fcn. (8988->4), ass. (0->203)
t249 = sin(qJ(3));
t245 = Icges(5,5) * t249;
t246 = Icges(6,4) * t249;
t408 = t245 + t246;
t248 = pkin(7) + qJ(2);
t243 = sin(t248);
t244 = cos(t248);
t250 = cos(qJ(3));
t336 = Icges(6,1) * t250;
t270 = t246 + t336;
t149 = Icges(6,5) * t244 + t243 * t270;
t326 = t243 * t249;
t222 = Icges(4,4) * t326;
t325 = t243 * t250;
t153 = Icges(4,1) * t325 - Icges(4,5) * t244 - t222;
t407 = -t149 - t153;
t337 = Icges(5,1) * t250;
t271 = t245 + t337;
t152 = Icges(5,4) * t243 + t244 * t271;
t335 = Icges(4,4) * t249;
t208 = Icges(4,1) * t250 - t335;
t154 = Icges(4,5) * t243 + t208 * t244;
t406 = t152 + t154;
t330 = Icges(5,3) * t250;
t331 = Icges(6,2) * t250;
t405 = -t330 - t331 + t408;
t209 = t249 * pkin(3) - qJ(4) * t250;
t210 = t249 * rSges(6,1) - rSges(6,2) * t250;
t350 = pkin(4) * t249;
t277 = t209 + t210 + t350;
t101 = t277 * t243;
t103 = t277 * t244;
t341 = t249 * rSges(5,1);
t211 = -rSges(5,3) * t250 + t341;
t289 = t209 + t211;
t132 = t289 * t243;
t134 = t289 * t244;
t321 = t244 * t249;
t320 = t244 * t250;
t237 = t244 * rSges(5,2);
t238 = t244 * pkin(6);
t339 = rSges(5,3) + qJ(4);
t352 = rSges(5,1) + pkin(3);
t379 = t339 * t249 + t352 * t250 + pkin(2);
t85 = -t379 * t243 + t237 + t238;
t86 = (rSges(5,2) + pkin(6)) * t243 + t379 * t244;
t344 = t85 * t320 + t86 * t325;
t351 = rSges(6,1) + pkin(4);
t284 = pkin(3) + t351;
t340 = rSges(6,2) + qJ(4);
t380 = t340 * t249 + t284 * t250 + pkin(2);
t389 = rSges(6,3) + qJ(5);
t76 = -t380 * t243 - t389 * t244 + t238;
t77 = t380 * t244 + (pkin(6) - t389) * t243;
t345 = t76 * t320 + t77 * t325;
t374 = m(6) / 0.2e1;
t375 = m(5) / 0.2e1;
t347 = (-t101 * t321 + t103 * t326 + t345) * t374 + (-t132 * t321 + t134 * t326 + t344) * t375;
t227 = pkin(3) * t326;
t110 = t227 + (-t339 * t250 + t341) * t243;
t217 = qJ(4) * t320;
t225 = rSges(5,3) * t320;
t111 = -t352 * t321 + t217 + t225;
t96 = t227 + (t351 * t249 - t340 * t250) * t243;
t226 = rSges(6,2) * t320;
t97 = -t284 * t321 + t217 + t226;
t348 = ((t243 * t97 + t244 * t96) * t249 + t345) * t374 + ((t110 * t244 + t111 * t243) * t249 + t344) * t375;
t2 = t348 - t347;
t241 = t243 ^ 2;
t242 = t244 ^ 2;
t285 = t241 + t242;
t188 = t285 * t249;
t283 = m(5) / 0.4e1 + m(6) / 0.4e1;
t274 = 0.2e1 * t283 * t188;
t282 = t374 + t375;
t88 = -t249 * t282 + t274;
t404 = t88 * qJD(1) - t2 * qJD(2);
t197 = Icges(4,5) * t250 - Icges(4,6) * t249;
t142 = Icges(4,3) * t243 + t244 * t197;
t200 = Icges(5,4) * t250 + Icges(5,6) * t249;
t146 = Icges(5,2) * t243 + t244 * t200;
t220 = Icges(5,5) * t320;
t140 = Icges(5,6) * t243 + Icges(5,3) * t321 + t220;
t318 = t249 * t140;
t403 = t244 * t318 + t406 * t320 + (t142 + t146) * t243;
t116 = t154 * t325;
t194 = Icges(6,5) * t250 + Icges(6,6) * t249;
t137 = Icges(6,3) * t244 + t243 * t194;
t402 = -t243 * t137 - t244 * t142 + t116;
t333 = Icges(5,5) * t250;
t196 = Icges(5,3) * t249 + t333;
t139 = -Icges(5,6) * t244 + t196 * t243;
t334 = Icges(6,4) * t250;
t199 = Icges(6,2) * t249 + t334;
t143 = Icges(6,6) * t244 + t199 * t243;
t203 = Icges(6,1) * t249 - t334;
t205 = Icges(5,1) * t249 - t333;
t401 = t139 + t143 + (-t203 - t205) * t243;
t141 = Icges(4,5) * t325 - Icges(4,6) * t326 - Icges(4,3) * t244;
t317 = t249 * t143;
t400 = -t243 * t141 - t244 * t317 + t407 * t320;
t399 = (-Icges(4,6) + Icges(5,6) - Icges(6,6)) * t250 + (-Icges(5,4) - Icges(4,5) + Icges(6,5)) * t249;
t150 = -Icges(6,5) * t243 + t244 * t270;
t201 = Icges(4,2) * t250 + t335;
t398 = t150 + (-t201 + t405) * t244 + t406;
t151 = -Icges(5,4) * t244 + t243 * t271;
t397 = -Icges(4,2) * t325 + t405 * t243 + t151 - t222 - t407;
t221 = Icges(6,4) * t320;
t144 = Icges(6,2) * t321 - Icges(6,6) * t243 + t221;
t247 = Icges(4,4) * t250;
t332 = Icges(4,2) * t249;
t148 = Icges(4,6) * t243 + (t247 - t332) * t244;
t338 = Icges(4,1) * t249;
t272 = -t247 - t338;
t396 = t272 * t244 + t140 + t144 - t148 + t220 + t221 + (-Icges(5,1) - Icges(6,1)) * t321;
t263 = t149 * t250 + t317;
t323 = t244 * (-Icges(5,2) * t244 + t243 * t200);
t319 = t249 * t139;
t388 = t243 * (t151 * t250 + t319);
t395 = t244 * t137 + t243 * t263 - t323 + t388;
t314 = t249 * t148;
t394 = -t244 * t314 + t403;
t393 = t323 + t403;
t147 = Icges(4,4) * t325 - Icges(4,2) * t326 - Icges(4,6) * t244;
t315 = t249 * t147;
t392 = -t243 * t314 - t244 * t315 - t400 + t402;
t391 = -t243 / 0.2e1;
t355 = t243 / 0.2e1;
t354 = -t244 / 0.2e1;
t390 = t244 / 0.2e1;
t313 = t249 * t250;
t290 = t285 * t313;
t387 = t283 * (t290 - t313);
t214 = rSges(6,1) * t250 + t249 * rSges(6,2);
t386 = pkin(4) * t250 + t214;
t385 = t399 * t243;
t384 = t399 * t244;
t383 = -t398 * t249 + t396 * t250;
t301 = -t272 * t243 + t147;
t382 = (t301 - t401) * t250 + t397 * t249;
t381 = -t249 * (t208 / 0.2e1 - t201 / 0.2e1 + t337 / 0.2e1 - t330 / 0.2e1 + t336 / 0.2e1 - t331 / 0.2e1 + t408) - t250 * (t247 + t338 / 0.2e1 - t332 / 0.2e1 + t205 / 0.2e1 - t196 / 0.2e1 + t203 / 0.2e1 - t199 / 0.2e1);
t378 = 0.4e1 * qJD(2);
t377 = 2 * qJD(3);
t342 = rSges(4,1) * t250;
t280 = pkin(2) + t342;
t287 = rSges(4,2) * t326 + t244 * rSges(4,3);
t108 = -t243 * t280 + t238 + t287;
t224 = rSges(4,2) * t321;
t109 = -t224 + t280 * t244 + (rSges(4,3) + pkin(6)) * t243;
t212 = t249 * rSges(4,1) + rSges(4,2) * t250;
t185 = t212 * t243;
t187 = t212 * t244;
t373 = m(4) * (t108 * t185 - t109 * t187);
t310 = -t132 * t325 - t134 * t320;
t215 = rSges(5,1) * t250 + t249 * rSges(5,3);
t213 = pkin(3) * t250 + t249 * qJ(4);
t292 = t285 * t213;
t52 = t243 * (t215 * t243 - t237) + (t243 * rSges(5,2) + t215 * t244) * t244 + t292;
t369 = m(5) * (t188 * t52 + t310);
t367 = m(5) * (t110 * t85 + t111 * t86);
t366 = m(5) * (t86 * t321 - t85 * t326);
t343 = -t101 * t325 - t103 * t320;
t46 = (pkin(4) * t325 + t214 * t243) * t243 + t292 + t386 * t242;
t363 = m(6) * (t188 * t46 + t343);
t360 = m(6) * (t76 * t96 + t77 * t97);
t359 = m(6) * (t77 * t321 - t76 * t326);
t358 = m(6) * (-t243 * t77 - t244 * t76);
t357 = m(6) * (-t243 * t96 + t244 * t97);
t356 = m(6) * (t101 * t243 + t103 * t244);
t98 = m(6) * t188;
t311 = t98 * qJD(2);
t308 = t144 * t321 + t150 * t320;
t293 = t243 * (qJ(4) * t325 - t227) + t244 * (-pkin(3) * t321 + t217);
t288 = -t213 - t215;
t281 = qJD(2) * t358;
t278 = -t141 + t314;
t276 = -t213 - t386;
t275 = t244 * t146 - t152 * t325 - t243 * t318;
t273 = -t194 / 0.2e1 + t200 / 0.2e1 + t197 / 0.2e1;
t260 = t242 * t137 / 0.2e1 + t275 * t391 + (t244 * t278 - t393 + t394) * t390 + (-t243 * (-t153 * t250 + t315) - t244 * t141 + t395) * t354 + (t151 * t320 + t243 * t278 + t244 * t319 + t275 + t392 - t402) * t355;
t138 = -Icges(6,3) * t243 + t244 * t194;
t253 = ((-t138 - t263) * t243 + t308 - t388 + t393 + t395) * t391 + (-t243 * t138 + t308 + t394) * t355 + (-t116 + (t142 + t315) * t244 + t392 + t400) * t354;
t216 = -t249 * rSges(4,2) + t342;
t135 = t288 * t244;
t133 = t288 * t243;
t104 = t276 * t244;
t102 = t276 * t243;
t92 = -t185 * t243 - t187 * t244;
t87 = t274 + (m(5) + m(6)) * t249 / 0.2e1;
t79 = 0.4e1 * t387;
t70 = t356 / 0.2e1;
t68 = t244 * (-rSges(5,1) * t321 + t225) - t211 * t241 + t293;
t54 = t357 / 0.2e1;
t49 = t244 * (-rSges(6,1) * t321 + t226) - t210 * t241 - t285 * t350 + t293;
t20 = t70 - t357 / 0.2e1;
t19 = t70 + t54;
t18 = t54 - t356 / 0.2e1;
t17 = t359 + t366;
t12 = t363 + t369;
t5 = t360 + t367 + t373 - t381;
t4 = t347 + t348;
t1 = t243 * t260 + t244 * t253;
t3 = [0, 0, t87 * qJD(4) + (m(4) * t92 / 0.2e1 + t68 * t375 + t49 * t374) * t377, t87 * qJD(3), 0; 0, t5 * qJD(3) + t17 * qJD(4) + qJD(5) * t358, t5 * qJD(2) + t4 * qJD(4) + t19 * qJD(5) + ((-t101 * t97 + t102 * t77 - t103 * t96 + t104 * t76) * t374 + (-t110 * t134 - t111 * t132 + t133 * t86 + t135 * t85) * t375) * t377 + ((m(4) * (-t108 * t216 - t185 * t212) + t273 * t244 - t253) * t244 + (m(4) * (-t109 * t216 + t187 * t212) + t273 * t243 - t260) * t243 + (t301 * t390 + t401 * t354 + t396 * t355) * t249 + (t397 * t354 + t398 * t355) * t250) * qJD(3), qJD(2) * t17 + qJD(3) * t4, t19 * qJD(3) + t281; qJD(4) * t88, t1 * qJD(3) - t2 * qJD(4) + t20 * qJD(5) + (-t360 / 0.4e1 - t367 / 0.4e1 - t373 / 0.4e1) * t378 + t381 * qJD(2), t1 * qJD(2) + t12 * qJD(4) + (m(6) * (-t101 * t102 - t103 * t104 + t46 * t49) + m(5) * (-t132 * t133 - t134 * t135 + t52 * t68) + m(4) * (t212 * t216 * t285 + (t243 * (rSges(4,1) * t325 - t287) + t244 * (rSges(4,1) * t320 + t243 * rSges(4,3) - t224)) * t92) + (t384 * t241 + (t382 * t244 + (t383 - t385) * t243) * t244) * t355 + (t385 * t242 + (t383 * t243 + (t382 - t384) * t244) * t243) * t354) * qJD(3), t12 * qJD(3) + (-0.4e1 * t387 + 0.2e1 * t282 * (-t188 * t250 + t290)) * qJD(4) + t404, t20 * qJD(2); -t88 * qJD(3), t2 * qJD(3) - t98 * qJD(5) + (-t366 / 0.4e1 - t359 / 0.4e1) * t378, t79 * qJD(4) + 0.4e1 * (-t363 / 0.4e1 - t369 / 0.4e1) * qJD(3) + ((-t250 * t49 + t343) * t374 + (-t250 * t68 + t310) * t375 + ((t102 * t243 + t104 * t244 + t46) * t374 + (t133 * t243 + t135 * t244 + t52) * t375) * t249) * t377 - t404, t79 * qJD(3), -t311; 0, t18 * qJD(3) + t98 * qJD(4) - t281, t18 * qJD(2) + m(6) * (t102 * t244 - t104 * t243) * qJD(3), t311, 0;];
Cq = t3;
