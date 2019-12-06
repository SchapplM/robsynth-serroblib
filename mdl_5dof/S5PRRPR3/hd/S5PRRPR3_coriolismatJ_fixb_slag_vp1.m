% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR3
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:04
% EndTime: 2019-12-05 16:19:17
% DurationCPUTime: 5.82s
% Computational Cost: add. (21913->373), mult. (14267->506), div. (0->0), fcn. (13084->8), ass. (0->243)
t286 = pkin(8) + qJ(2);
t281 = cos(t286);
t287 = qJ(3) + pkin(9);
t283 = qJ(5) + t287;
t274 = sin(t283);
t275 = cos(t283);
t221 = rSges(6,1) * t274 + rSges(6,2) * t275;
t280 = sin(t287);
t289 = sin(qJ(3));
t414 = pkin(3) * t289;
t243 = -pkin(4) * t280 - t414;
t302 = t221 - t243;
t125 = t302 * t281;
t279 = sin(t286);
t455 = t302 * t279;
t462 = t279 * t455;
t453 = t125 * t281 + t462;
t282 = cos(t287);
t295 = rSges(5,1) * t280 + rSges(5,2) * t282 + t414;
t456 = t295 * t281;
t457 = t295 * t279;
t461 = t279 * t457 + t281 * t456;
t407 = -m(6) * t453 / 0.2e1 - m(5) * t461 / 0.2e1;
t190 = t221 * t281;
t215 = t281 * t243;
t128 = t215 - t190;
t440 = m(6) / 0.2e1;
t441 = m(5) / 0.2e1;
t408 = (-t128 * t281 + t462) * t440 + t461 * t441;
t23 = t408 - t407;
t468 = t23 * qJD(2);
t277 = t279 ^ 2;
t278 = t281 ^ 2;
t334 = t277 + t278;
t467 = t221 * t334;
t259 = Icges(6,4) * t275;
t218 = -Icges(6,2) * t274 + t259;
t219 = Icges(6,1) * t274 + t259;
t466 = t218 + t219;
t290 = cos(qJ(3));
t465 = -Icges(4,5) * t289 - Icges(5,5) * t280 - Icges(4,6) * t290 - Icges(5,6) * t282;
t422 = t279 / 0.2e1;
t421 = -t281 / 0.2e1;
t463 = t281 / 0.2e1;
t393 = Icges(5,4) * t280;
t224 = Icges(5,2) * t282 + t393;
t227 = Icges(5,1) * t282 - t393;
t394 = Icges(4,4) * t289;
t245 = Icges(4,2) * t290 + t394;
t248 = Icges(4,1) * t290 - t394;
t460 = -(t227 / 0.2e1 - t224 / 0.2e1) * t280 - (t248 / 0.2e1 - t245 / 0.2e1) * t289;
t404 = rSges(6,1) * t275;
t222 = -rSges(6,2) * t274 + t404;
t307 = Icges(6,5) * t274 + Icges(6,6) * t275;
t177 = t307 * t279;
t178 = t281 * t307;
t392 = Icges(6,4) * t274;
t220 = Icges(6,1) * t275 - t392;
t158 = Icges(6,5) * t279 + t220 * t281;
t217 = Icges(6,2) * t275 + t392;
t346 = -t217 * t281 + t158;
t377 = t274 * t279;
t234 = Icges(6,4) * t377;
t373 = t275 * t279;
t157 = Icges(6,1) * t373 - Icges(6,5) * t281 - t234;
t347 = -Icges(6,2) * t373 + t157 - t234;
t156 = Icges(6,6) * t279 + t218 * t281;
t348 = -t219 * t281 - t156;
t155 = Icges(6,4) * t373 - Icges(6,2) * t377 - Icges(6,6) * t281;
t349 = t219 * t279 + t155;
t445 = (-t279 * t346 + t281 * t347) * t274 + (t279 * t348 + t281 * t349) * t275;
t410 = (-t278 * t177 + (t281 * t178 + t445) * t279) * t421 + (-t277 * t178 + (t279 * t177 + t445) * t281) * t422;
t285 = t290 * pkin(3);
t276 = t285 + pkin(2);
t413 = pkin(4) * t282;
t237 = t276 + t413;
t288 = -qJ(4) - pkin(6);
t333 = -pkin(7) + t288;
t256 = t279 * t333;
t258 = t279 * t288;
t337 = -t237 * t279 - t281 * t333;
t273 = t281 * pkin(6);
t362 = t281 * t288;
t411 = pkin(2) - t276;
t352 = -t279 * (t279 * t411 - t273 - t362) + t281 * (-t279 * pkin(6) - t281 * t411 - t258);
t160 = rSges(6,1) * t373 - rSges(6,2) * t377 - rSges(6,3) * t281;
t376 = t274 * t281;
t326 = -rSges(6,2) * t376 + rSges(6,3) * t279;
t372 = t275 * t281;
t96 = t279 * t160 + t281 * (rSges(6,1) * t372 + t326);
t45 = -t279 * (t276 * t279 + t337 + t362) + (-t256 + t258 + (t237 - t276) * t281) * t281 + t96 + t352;
t189 = t221 * t279;
t452 = t189 * t279 + t190 * t281;
t6 = t410 + m(6) * (t222 * t453 - t45 * t452);
t458 = t6 * qJD(5);
t110 = -t160 + t337;
t111 = -t256 + (t237 + t404) * t281 + t326;
t454 = t110 * t281 + t111 * t279;
t451 = t465 * t279;
t450 = t465 * t281;
t269 = Icges(5,4) * t282;
t225 = -Icges(5,2) * t280 + t269;
t226 = Icges(5,1) * t280 + t269;
t284 = Icges(4,4) * t290;
t246 = -Icges(4,2) * t289 + t284;
t247 = Icges(4,1) * t289 + t284;
t188 = Icges(4,5) * t279 + t248 * t281;
t338 = -t245 * t281 + t188;
t186 = Icges(4,6) * t279 + t246 * t281;
t340 = -t247 * t281 - t186;
t171 = Icges(5,5) * t279 + t227 * t281;
t342 = -t224 * t281 + t171;
t169 = Icges(5,6) * t279 + t225 * t281;
t344 = -t226 * t281 - t169;
t447 = -t280 * t342 + t282 * t344 - t289 * t338 + t290 * t340;
t369 = t279 * t289;
t253 = Icges(4,4) * t369;
t368 = t279 * t290;
t187 = Icges(4,1) * t368 - Icges(4,5) * t281 - t253;
t339 = -Icges(4,2) * t368 + t187 - t253;
t185 = Icges(4,4) * t368 - Icges(4,2) * t369 - Icges(4,6) * t281;
t341 = t247 * t279 + t185;
t371 = t279 * t280;
t240 = Icges(5,4) * t371;
t370 = t279 * t282;
t170 = Icges(5,1) * t370 - Icges(5,5) * t281 - t240;
t343 = -Icges(5,2) * t370 + t170 - t240;
t168 = Icges(5,4) * t370 - Icges(5,2) * t371 - Icges(5,6) * t281;
t345 = t226 * t279 + t168;
t446 = t280 * t343 + t282 * t345 + t289 * t339 + t290 * t341;
t317 = t466 * t275 / 0.2e1 + (-t217 / 0.2e1 + t220 / 0.2e1) * t274;
t120 = t158 * t373;
t153 = Icges(6,5) * t373 - Icges(6,6) * t377 - Icges(6,3) * t281;
t216 = Icges(6,5) * t275 - Icges(6,6) * t274;
t382 = t216 * t281;
t154 = Icges(6,3) * t279 + t382;
t320 = t156 * t274 - t153;
t323 = t154 * t281 - t120;
t356 = t154 * t279 + t158 * t372;
t357 = -t279 * t153 - t157 * t372;
t385 = t155 * t274;
t66 = -t156 * t377 - t323;
t67 = -t155 * t376 - t357;
t68 = -t156 * t376 + t356;
t329 = ((t66 - t120 + (t154 + t385) * t281 + t357) * t281 + t356 * t279) * t421 + (t279 * t68 - t281 * t67) * t463 + ((t320 * t279 + t323 + t66 + t67) * t279 + (-t356 + t68 + (-t157 * t275 + t385) * t279 + (t320 + t153) * t281) * t281) * t422;
t444 = 0.4e1 * qJD(2);
t443 = 2 * qJD(3);
t442 = m(4) / 0.2e1;
t406 = rSges(4,1) * t290;
t331 = pkin(2) + t406;
t335 = rSges(4,2) * t369 + rSges(4,3) * t281;
t138 = -t279 * t331 + t273 + t335;
t398 = t289 * rSges(4,2);
t255 = t281 * t398;
t139 = -t255 + t331 * t281 + (rSges(4,3) + pkin(6)) * t279;
t249 = rSges(4,1) * t289 + rSges(4,2) * t290;
t212 = t249 * t279;
t213 = t249 * t281;
t439 = m(4) * (t138 * t212 - t139 * t213);
t405 = rSges(5,1) * t282;
t328 = t276 + t405;
t336 = rSges(5,2) * t371 + rSges(5,3) * t281;
t115 = -t279 * t328 + t336 - t362;
t367 = t280 * t281;
t327 = -rSges(5,2) * t367 + rSges(5,3) * t279;
t116 = t281 * t328 - t258 + t327;
t437 = m(5) * (t115 * t457 - t116 * t456);
t436 = m(5) * (t115 * t281 + t116 * t279);
t297 = t454 * t222;
t431 = m(6) * (-t125 * t189 + t190 * t455 - t297);
t430 = m(6) * (-t297 + (-t128 * t279 - t281 * t455) * t221);
t429 = m(6) * (t110 * t455 + t111 * t128);
t428 = m(6) * (t110 * t189 - t111 * t190);
t427 = m(6) * t454;
t423 = -t279 / 0.2e1;
t415 = m(6) * t452;
t383 = t168 * t280;
t363 = t281 * t282;
t361 = t281 * t290;
t360 = t289 * t185;
t359 = t289 * t186;
t299 = t452 * t440;
t315 = m(6) * t467;
t74 = t299 + t315 / 0.2e1;
t358 = t74 * qJD(2);
t166 = Icges(5,5) * t370 - Icges(5,6) * t371 - Icges(5,3) * t281;
t355 = -t279 * t166 - t170 * t363;
t309 = Icges(5,5) * t282 - Icges(5,6) * t280;
t167 = Icges(5,3) * t279 + t281 * t309;
t354 = t167 * t279 + t171 * t363;
t183 = Icges(4,5) * t368 - Icges(4,6) * t369 - Icges(4,3) * t281;
t351 = -t279 * t183 - t187 * t361;
t311 = Icges(4,5) * t290 - Icges(4,6) * t289;
t184 = Icges(4,3) * t279 + t281 * t311;
t350 = t184 * t279 + t188 * t361;
t330 = rSges(5,2) * t280 - t285 - t405;
t129 = t171 * t370;
t322 = t167 * t281 - t129;
t149 = t188 * t368;
t321 = t281 * t184 - t149;
t319 = t169 * t280 - t166;
t318 = -t183 + t359;
t301 = -t222 - t285 - t413;
t291 = (-t217 + t220) * t275 - t466 * t274;
t296 = -t329 + (t216 * t279 + t274 * t348 + t275 * t346 + t281 * t291) * t422 + (-t274 * t349 + t275 * t347 + t279 * t291 - t382) * t421;
t294 = -t317 + (t422 + t423) * (t155 * t275 + t157 * t274);
t250 = -t398 + t406;
t176 = t330 * t281;
t174 = t330 * t279;
t126 = t301 * t281;
t124 = t301 * t279;
t114 = -t212 * t279 - t213 * t281;
t102 = qJD(5) * t415;
t94 = t295 * t334;
t83 = -t281 * t359 + t350;
t82 = -t281 * t360 - t351;
t81 = -t279 * t359 - t321;
t73 = t299 - t315 / 0.2e1;
t72 = -t169 * t367 + t354;
t71 = -t168 * t367 - t355;
t70 = -t169 * t371 - t322;
t60 = t281 * t215 + t243 * t277 - t452;
t52 = t279 * t83 - t281 * t82;
t51 = t279 * t81 - t281 * (-(-t187 * t290 + t360) * t279 - t281 * t183);
t49 = t279 * t72 - t281 * t71;
t48 = t279 * t70 - t281 * (-(-t170 * t282 + t383) * t279 - t166 * t281);
t47 = t427 + t436;
t41 = t317 + t428;
t35 = t430 / 0.2e1;
t33 = t431 / 0.2e1;
t25 = t407 + t408;
t22 = (t81 - t149 + (t184 + t360) * t281 + t351) * t281 + t350 * t279;
t21 = (t281 * t318 - t350 + t83) * t281 + (t279 * t318 + t321 + t82) * t279;
t19 = (t70 - t129 + (t167 + t383) * t281 + t355) * t281 + t354 * t279;
t18 = (t281 * t319 - t354 + t72) * t281 + (t279 * t319 + t322 + t71) * t279;
t17 = (t247 / 0.2e1 + t246 / 0.2e1) * t290 + (t226 / 0.2e1 + t225 / 0.2e1) * t282 + t439 + t437 + t429 + t317 - t460;
t8 = m(6) * (t222 * t467 - t452 * t96) + t410;
t7 = t8 * qJD(5);
t4 = t33 - t430 / 0.2e1 + t329;
t3 = t35 - t431 / 0.2e1 + t329;
t2 = t33 + t35 + t296;
t1 = (-t19 / 0.2e1 - t22 / 0.2e1 + t49 / 0.2e1 + t52 / 0.2e1) * t281 + (t48 / 0.2e1 + t51 / 0.2e1 + t18 / 0.2e1 + t21 / 0.2e1) * t279 + t329;
t5 = [0, 0, -t102 + (t114 * t442 + t440 * t60 - t441 * t94) * t443, 0, -qJD(3) * t415 - t102; 0, qJD(3) * t17 + qJD(4) * t47 + qJD(5) * t41, t17 * qJD(2) + t25 * qJD(4) + t2 * qJD(5) + (((-t138 * t281 - t139 * t279) * t250 + (-t212 * t281 + t213 * t279) * t249) * t442 + (t115 * t176 + t116 * t174) * t441 + (t110 * t126 + t111 * t124 + (-t125 - t128) * t455) * t440) * t443 + ((t311 + t309) * (t277 / 0.2e1 + t278 / 0.2e1) + t296 + (t280 * t344 + t282 * t342 + t289 * t340 + t290 * t338) * t422 + (t19 + t22) * t463 + (t48 + t51 + t18 + t21) * t423 + (-t280 * t345 + t282 * t343 - t289 * t341 + t290 * t339 + t49 + t52) * t421) * qJD(3), qJD(2) * t47 + qJD(3) * t25 + qJD(5) * t73, t41 * qJD(2) + t2 * qJD(3) + t73 * qJD(4) + (t296 + (-t297 + (-t189 * t281 + t190 * t279) * t221) * m(6)) * qJD(5); 0, (t294 - (t226 + t225) * t282 / 0.2e1 - (t247 + t246) * t290 / 0.2e1 + t460) * qJD(2) + t1 * qJD(3) - t23 * qJD(4) + t4 * qJD(5) + (-t439 / 0.4e1 - t429 / 0.4e1 - t437 / 0.4e1) * t444, t1 * qJD(2) + (m(6) * (-t124 * t455 - t125 * t126 + t45 * t60) + m(5) * (-t457 * t174 - t456 * t176 - (t279 * (rSges(5,1) * t370 - t336) + t281 * (rSges(5,1) * t363 + t327) + t352) * t94) + m(4) * (t249 * t250 * t334 + (t279 * (rSges(4,1) * t368 - t335) + t281 * (rSges(4,1) * t361 + t279 * rSges(4,3) - t255)) * t114) + t410 + (t450 * t277 + (t446 * t281 + (t447 - t451) * t279) * t281) * t422 + ((t447 * t279 + (t446 - t450) * t281) * t279 + t451 * t278) * t421) * qJD(3) + t458, -t468, t4 * qJD(2) + t6 * qJD(3) + t458; 0, t23 * qJD(3) + t74 * qJD(5) + (-t427 / 0.4e1 - t436 / 0.4e1) * t444, t468 + ((-t124 * t281 + t126 * t279) * t440 + (-t174 * t281 + t176 * t279) * t441) * t443, 0, t358; 0, (t294 - t428) * qJD(2) + t3 * qJD(3) - t74 * qJD(4) + t329 * qJD(5), t3 * qJD(2) + ((t60 * t96 + (-t124 * t279 - t126 * t281) * t221) * m(6) + t410) * qJD(3) + t7, -t358, qJD(2) * t329 + qJD(3) * t8 + t7;];
Cq = t5;
