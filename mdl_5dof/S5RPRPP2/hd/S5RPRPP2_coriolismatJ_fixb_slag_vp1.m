% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:41
% EndTime: 2019-12-31 18:10:49
% DurationCPUTime: 5.15s
% Computational Cost: add. (9433->322), mult. (10146->434), div. (0->0), fcn. (9118->6), ass. (0->200)
t251 = sin(qJ(3));
t247 = Icges(5,5) * t251;
t248 = Icges(6,4) * t251;
t409 = t247 + t248;
t250 = qJ(1) + pkin(7);
t245 = sin(t250);
t246 = cos(t250);
t253 = cos(qJ(3));
t336 = Icges(6,1) * t253;
t274 = t248 + t336;
t149 = Icges(6,5) * t246 + t245 * t274;
t323 = t245 * t251;
t224 = Icges(4,4) * t323;
t322 = t245 * t253;
t153 = Icges(4,1) * t322 - Icges(4,5) * t246 - t224;
t408 = -t149 - t153;
t337 = Icges(5,1) * t253;
t275 = t247 + t337;
t152 = Icges(5,4) * t245 + t246 * t275;
t335 = Icges(4,4) * t251;
t208 = Icges(4,1) * t253 - t335;
t154 = Icges(4,5) * t245 + t208 * t246;
t407 = t152 + t154;
t330 = Icges(5,3) * t253;
t331 = Icges(6,2) * t253;
t406 = -t330 - t331 + t409;
t209 = pkin(3) * t251 - qJ(4) * t253;
t210 = rSges(6,1) * t251 - rSges(6,2) * t253;
t350 = pkin(4) * t251;
t281 = t209 + t210 + t350;
t103 = t281 * t245;
t105 = t281 * t246;
t341 = rSges(5,1) * t251;
t211 = -rSges(5,3) * t253 + t341;
t294 = t209 + t211;
t132 = t294 * t245;
t134 = t294 * t246;
t320 = t246 * t251;
t319 = t246 * t253;
t239 = t246 * rSges(5,2);
t284 = -sin(qJ(1)) * pkin(1) + t246 * pkin(6);
t339 = rSges(5,3) + qJ(4);
t353 = rSges(5,1) + pkin(3);
t380 = t251 * t339 + t253 * t353 + pkin(2);
t85 = -t380 * t245 + t239 + t284;
t351 = cos(qJ(1)) * pkin(1);
t86 = t351 + (rSges(5,2) + pkin(6)) * t245 + t380 * t246;
t344 = t85 * t319 + t86 * t322;
t352 = rSges(6,1) + pkin(4);
t289 = pkin(3) + t352;
t340 = rSges(6,2) + qJ(4);
t381 = t251 * t340 + t253 * t289 + pkin(2);
t390 = rSges(6,3) + qJ(5);
t76 = -t381 * t245 - t390 * t246 + t284;
t77 = t381 * t246 + t351 + (pkin(6) - t390) * t245;
t345 = t76 * t319 + t77 * t322;
t375 = m(6) / 0.2e1;
t376 = m(5) / 0.2e1;
t347 = (-t103 * t320 + t105 * t323 + t345) * t375 + (-t132 * t320 + t134 * t323 + t344) * t376;
t229 = pkin(3) * t323;
t110 = t229 + (-t253 * t339 + t341) * t245;
t219 = qJ(4) * t319;
t227 = rSges(5,3) * t319;
t111 = -t320 * t353 + t219 + t227;
t96 = t229 + (t251 * t352 - t253 * t340) * t245;
t228 = rSges(6,2) * t319;
t97 = -t289 * t320 + t219 + t228;
t348 = ((t245 * t97 + t246 * t96) * t251 + t345) * t375 + ((t110 * t246 + t111 * t245) * t251 + t344) * t376;
t2 = t348 - t347;
t243 = t245 ^ 2;
t244 = t246 ^ 2;
t290 = t243 + t244;
t190 = t290 * t251;
t287 = m(5) / 0.4e1 + m(6) / 0.4e1;
t278 = 0.2e1 * t287 * t190;
t288 = t376 + t375;
t88 = -t251 * t288 + t278;
t405 = -t2 * qJD(1) + t88 * qJD(2);
t222 = Icges(5,5) * t319;
t140 = Icges(5,6) * t245 + Icges(5,3) * t320 + t222;
t197 = Icges(4,5) * t253 - Icges(4,6) * t251;
t142 = Icges(4,3) * t245 + t197 * t246;
t200 = Icges(5,4) * t253 + Icges(5,6) * t251;
t146 = Icges(5,2) * t245 + t200 * t246;
t404 = t140 * t320 + t407 * t319 + (t142 + t146) * t245;
t116 = t154 * t322;
t194 = Icges(6,5) * t253 + Icges(6,6) * t251;
t137 = Icges(6,3) * t246 + t194 * t245;
t403 = -t137 * t245 - t142 * t246 + t116;
t333 = Icges(5,5) * t253;
t196 = Icges(5,3) * t251 + t333;
t139 = -Icges(5,6) * t246 + t196 * t245;
t334 = Icges(6,4) * t253;
t199 = Icges(6,2) * t251 + t334;
t143 = Icges(6,6) * t246 + t199 * t245;
t203 = Icges(6,1) * t251 - t334;
t205 = Icges(5,1) * t251 - t333;
t402 = t139 + t143 + (-t203 - t205) * t245;
t141 = Icges(4,5) * t322 - Icges(4,6) * t323 - Icges(4,3) * t246;
t401 = -t245 * t141 - t143 * t320 + t408 * t319;
t400 = (-Icges(4,6) + Icges(5,6) - Icges(6,6)) * t253 + (-Icges(5,4) - Icges(4,5) + Icges(6,5)) * t251;
t150 = -Icges(6,5) * t245 + t246 * t274;
t201 = Icges(4,2) * t253 + t335;
t399 = t150 + (-t201 + t406) * t246 + t407;
t151 = -Icges(5,4) * t246 + t245 * t275;
t398 = -Icges(4,2) * t322 + t406 * t245 + t151 - t224 - t408;
t223 = Icges(6,4) * t319;
t144 = Icges(6,2) * t320 - Icges(6,6) * t245 + t223;
t249 = Icges(4,4) * t253;
t332 = Icges(4,2) * t251;
t148 = Icges(4,6) * t245 + (t249 - t332) * t246;
t338 = Icges(4,1) * t251;
t276 = -t249 - t338;
t397 = t276 * t246 + t140 + t144 - t148 + t222 + t223 + (-Icges(5,1) - Icges(6,1)) * t320;
t267 = t143 * t251 + t149 * t253;
t327 = (-Icges(5,2) * t246 + t200 * t245) * t246;
t389 = (t139 * t251 + t151 * t253) * t245;
t396 = t137 * t246 + t245 * t267 - t327 + t389;
t395 = -t148 * t320 + t404;
t394 = t327 + t404;
t147 = Icges(4,4) * t322 - Icges(4,2) * t323 - Icges(4,6) * t246;
t393 = -t147 * t320 - t148 * t323 - t401 + t403;
t392 = -t245 / 0.2e1;
t356 = t245 / 0.2e1;
t355 = -t246 / 0.2e1;
t391 = t246 / 0.2e1;
t318 = t251 * t253;
t295 = t290 * t318;
t388 = t287 * (t295 - t318);
t215 = rSges(6,1) * t253 + rSges(6,2) * t251;
t387 = pkin(4) * t253 + t215;
t386 = t400 * t245;
t385 = t400 * t246;
t384 = -t399 * t251 + t397 * t253;
t306 = -t276 * t245 + t147;
t383 = (t306 - t402) * t253 + t398 * t251;
t382 = -t251 * (t208 / 0.2e1 - t201 / 0.2e1 + t337 / 0.2e1 - t330 / 0.2e1 + t336 / 0.2e1 - t331 / 0.2e1 + t409) - t253 * (t249 + t338 / 0.2e1 - t332 / 0.2e1 + t205 / 0.2e1 - t196 / 0.2e1 + t203 / 0.2e1 - t199 / 0.2e1);
t379 = 0.4e1 * qJD(1);
t378 = 2 * qJD(3);
t342 = rSges(4,1) * t253;
t285 = pkin(2) + t342;
t292 = rSges(4,2) * t323 + t246 * rSges(4,3);
t100 = -t245 * t285 + t284 + t292;
t226 = rSges(4,2) * t320;
t101 = t351 - t226 + t285 * t246 + (rSges(4,3) + pkin(6)) * t245;
t212 = rSges(4,1) * t251 + rSges(4,2) * t253;
t185 = t212 * t245;
t187 = t212 * t246;
t374 = m(4) * (t100 * t185 - t101 * t187);
t315 = -t132 * t322 - t134 * t319;
t216 = rSges(5,1) * t253 + rSges(5,3) * t251;
t214 = pkin(3) * t253 + qJ(4) * t251;
t297 = t290 * t214;
t52 = t245 * (t216 * t245 - t239) + (t245 * rSges(5,2) + t216 * t246) * t246 + t297;
t370 = m(5) * (t190 * t52 + t315);
t369 = m(5) * (t110 * t85 + t111 * t86);
t367 = m(5) * (t86 * t320 - t323 * t85);
t343 = -t103 * t322 - t105 * t319;
t46 = (pkin(4) * t322 + t215 * t245) * t245 + t297 + t387 * t244;
t364 = m(6) * (t190 * t46 + t343);
t361 = m(6) * (t76 * t96 + t77 * t97);
t360 = m(6) * (t77 * t320 - t323 * t76);
t359 = m(6) * (-t245 * t77 - t246 * t76);
t358 = m(6) * (-t245 * t96 + t246 * t97);
t357 = m(6) * (t103 * t245 + t105 * t246);
t326 = t147 * t251;
t98 = m(6) * t190;
t316 = t98 * qJD(1);
t313 = t144 * t320 + t150 * t319;
t298 = t245 * (qJ(4) * t322 - t229) + t246 * (-pkin(3) * t320 + t219);
t293 = -t214 - t216;
t286 = qJD(1) * t359;
t282 = t148 * t251 - t141;
t280 = -t214 - t387;
t279 = -t140 * t323 + t146 * t246 - t152 * t322;
t277 = t197 / 0.2e1 - t194 / 0.2e1 + t200 / 0.2e1;
t264 = t279 * t392 + t137 * t244 / 0.2e1 + (t246 * t282 - t394 + t395) * t391 + (-t245 * (-t153 * t253 + t326) - t141 * t246 + t396) * t355 + (t139 * t320 + t151 * t319 + t245 * t282 + t279 + t393 - t403) * t356;
t138 = -Icges(6,3) * t245 + t194 * t246;
t257 = (-t389 + (-t138 - t267) * t245 + t313 + t394 + t396) * t392 + (-t138 * t245 + t313 + t395) * t356 + (-t116 + (t142 + t326) * t246 + t393 + t401) * t355;
t217 = -rSges(4,2) * t251 + t342;
t135 = t293 * t246;
t133 = t293 * t245;
t106 = t280 * t246;
t104 = t280 * t245;
t92 = -t185 * t245 - t187 * t246;
t87 = t278 + (m(5) + m(6)) * t251 / 0.2e1;
t82 = 0.4e1 * t388;
t70 = t357 / 0.2e1;
t68 = t246 * (-rSges(5,1) * t320 + t227) - t211 * t243 + t298;
t54 = t358 / 0.2e1;
t50 = t246 * (-rSges(6,1) * t320 + t228) - t210 * t243 - t290 * t350 + t298;
t20 = t70 - t358 / 0.2e1;
t19 = t70 + t54;
t18 = t54 - t357 / 0.2e1;
t17 = t360 + t367;
t12 = t364 + t370;
t5 = t361 + t369 + t374 - t382;
t3 = t347 + t348;
t1 = t245 * t264 + t246 * t257;
t4 = [t5 * qJD(3) + t17 * qJD(4) + qJD(5) * t359, 0, t5 * qJD(1) + t3 * qJD(4) + t19 * qJD(5) + ((-t103 * t97 + t104 * t77 - t105 * t96 + t106 * t76) * t375 + (-t110 * t134 - t111 * t132 + t133 * t86 + t135 * t85) * t376) * t378 + ((m(4) * (-t100 * t217 - t185 * t212) + t277 * t246 - t257) * t246 + (m(4) * (-t101 * t217 + t187 * t212) + t277 * t245 - t264) * t245 + (t306 * t391 + t402 * t355 + t397 * t356) * t251 + (t398 * t355 + t399 * t356) * t253) * qJD(3), qJD(1) * t17 + qJD(3) * t3, t19 * qJD(3) + t286; 0, 0, t87 * qJD(4) + (m(4) * t92 / 0.2e1 + t68 * t376 + t50 * t375) * t378, t87 * qJD(3), 0; t1 * qJD(3) - t2 * qJD(4) + t20 * qJD(5) + (-t361 / 0.4e1 - t369 / 0.4e1 - t374 / 0.4e1) * t379 + t382 * qJD(1), qJD(4) * t88, t1 * qJD(1) + t12 * qJD(4) + (m(5) * (-t132 * t133 - t134 * t135 + t52 * t68) + m(6) * (-t103 * t104 - t105 * t106 + t46 * t50) + m(4) * (t290 * t217 * t212 + (t245 * (rSges(4,1) * t322 - t292) + t246 * (rSges(4,1) * t319 + t245 * rSges(4,3) - t226)) * t92) + (t385 * t243 + (t383 * t246 + (t384 - t386) * t245) * t246) * t356 + (t386 * t244 + (t384 * t245 + (t383 - t385) * t246) * t245) * t355) * qJD(3), t12 * qJD(3) + (-0.4e1 * t388 + 0.2e1 * t288 * (-t190 * t253 + t295)) * qJD(4) + t405, t20 * qJD(1); t2 * qJD(3) - t98 * qJD(5) + (-t360 / 0.4e1 - t367 / 0.4e1) * t379, -t88 * qJD(3), t82 * qJD(4) + 0.4e1 * (-t370 / 0.4e1 - t364 / 0.4e1) * qJD(3) + ((-t253 * t68 + t315) * t376 + (-t253 * t50 + t343) * t375 + ((t133 * t245 + t135 * t246 + t52) * t376 + (t104 * t245 + t106 * t246 + t46) * t375) * t251) * t378 - t405, t82 * qJD(3), -t316; t18 * qJD(3) + t98 * qJD(4) - t286, 0, t18 * qJD(1) + m(6) * (t104 * t246 - t106 * t245) * qJD(3), t316, 0;];
Cq = t4;
