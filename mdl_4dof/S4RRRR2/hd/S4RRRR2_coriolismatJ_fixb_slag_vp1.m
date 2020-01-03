% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:08
% EndTime: 2019-12-31 17:23:14
% DurationCPUTime: 3.86s
% Computational Cost: add. (28527->327), mult. (22428->437), div. (0->0), fcn. (20532->8), ass. (0->215)
t293 = qJ(3) + qJ(4);
t290 = cos(t293);
t280 = Icges(5,4) * t290;
t288 = sin(t293);
t245 = -Icges(5,2) * t288 + t280;
t246 = Icges(5,1) * t288 + t280;
t455 = t245 + t246;
t297 = cos(qJ(3));
t292 = Icges(4,4) * t297;
t295 = sin(qJ(3));
t260 = -Icges(4,2) * t295 + t292;
t261 = Icges(4,1) * t295 + t292;
t454 = t260 + t261;
t294 = qJ(1) + qJ(2);
t291 = cos(t294);
t289 = sin(t294);
t371 = t289 * t290;
t375 = t288 * t289;
t193 = rSges(5,1) * t371 - rSges(5,2) * t375 - t291 * rSges(5,3);
t402 = pkin(3) * t297;
t285 = pkin(2) + t402;
t435 = -pkin(7) - pkin(6);
t339 = -t289 * t285 - t291 * t435;
t160 = -t193 + t339;
t272 = t289 * t435;
t374 = t288 * t291;
t333 = -rSges(5,2) * t374 + t289 * rSges(5,3);
t398 = rSges(5,1) * t290;
t161 = -t272 + (t285 + t398) * t291 + t333;
t248 = rSges(5,1) * t288 + rSges(5,2) * t290;
t403 = pkin(3) * t295;
t310 = t248 + t403;
t446 = t310 * t291;
t447 = t310 * t289;
t82 = t160 * t447 - t161 * t446;
t284 = t291 * pkin(6);
t399 = rSges(4,1) * t297;
t336 = pkin(2) + t399;
t370 = t289 * t295;
t338 = rSges(4,2) * t370 + t291 * rSges(4,3);
t167 = -t336 * t289 + t284 + t338;
t362 = t291 * t295;
t271 = rSges(4,2) * t362;
t168 = -t271 + t336 * t291 + (rSges(4,3) + pkin(6)) * t289;
t263 = rSges(4,1) * t295 + rSges(4,2) * t297;
t231 = t263 * t289;
t232 = t263 * t291;
t97 = t167 * t231 - t168 * t232;
t453 = m(4) * t97 + m(5) * t82;
t437 = m(4) / 0.2e1;
t436 = m(5) / 0.2e1;
t409 = t289 / 0.2e1;
t408 = -t291 / 0.2e1;
t451 = t291 / 0.2e1;
t404 = cos(qJ(1)) * pkin(1);
t405 = sin(qJ(1)) * pkin(1);
t407 = m(3) * (t404 * (-rSges(3,1) * t289 - rSges(3,2) * t291) + (t291 * rSges(3,1) - t289 * rSges(3,2)) * t405);
t286 = t289 ^ 2;
t287 = t291 ^ 2;
t337 = t286 + t287;
t410 = -t289 / 0.2e1;
t449 = t409 + t410;
t155 = t160 - t405;
t156 = t161 + t404;
t81 = t155 * t447 - t156 * t446;
t215 = t248 * t289;
t216 = t248 * t291;
t153 = -t289 * t215 - t291 * t216;
t250 = -rSges(5,2) * t288 + t398;
t320 = Icges(5,5) * t288 + Icges(5,6) * t290;
t209 = t320 * t289;
t210 = t291 * t320;
t389 = Icges(5,4) * t288;
t247 = Icges(5,1) * t290 - t389;
t192 = Icges(5,5) * t289 + t247 * t291;
t244 = Icges(5,2) * t290 + t389;
t344 = -t244 * t291 + t192;
t254 = Icges(5,4) * t375;
t191 = Icges(5,1) * t371 - Icges(5,5) * t291 - t254;
t345 = -Icges(5,2) * t371 + t191 - t254;
t190 = Icges(5,6) * t289 + t245 * t291;
t346 = -t246 * t291 - t190;
t189 = Icges(5,4) * t371 - Icges(5,2) * t375 - Icges(5,6) * t291;
t347 = t246 * t289 + t189;
t442 = (-t344 * t289 + t345 * t291) * t288 + (t346 * t289 + t347 * t291) * t290;
t401 = (-t286 * t210 + (t289 * t209 + t442) * t291) * t409 + (-t287 * t209 + (t291 * t210 + t442) * t289) * t408;
t366 = t290 * t291;
t137 = t289 * t193 + t291 * (rSges(5,1) * t366 + t333);
t86 = -t289 * (pkin(2) * t289 - t284 + t339) + (-t289 * pkin(6) - t272 + (-pkin(2) + t285) * t291) * t291 + t137;
t15 = t401 + m(5) * (t86 * t153 + (t289 * t447 + t291 * t446) * t250);
t448 = t15 * qJD(4);
t77 = -t161 * t155 + t156 * t160;
t162 = t167 - t405;
t163 = t168 + t404;
t83 = -t168 * t162 + t163 * t167;
t445 = qJD(1) + qJD(2);
t400 = (t81 - t82) * t436 + ((-t163 + t168) * t291 + (t162 - t167) * t289) * t263 * t437;
t94 = t162 * t231 - t163 * t232;
t444 = (t82 + t81) * t436 + (t97 + t94) * t437;
t390 = Icges(4,4) * t295;
t259 = Icges(4,2) * t297 + t390;
t262 = Icges(4,1) * t297 - t390;
t443 = t454 * t297 / 0.2e1 + (t262 / 0.2e1 - t259 / 0.2e1) * t295;
t326 = t455 * t290 / 0.2e1 + (-t244 / 0.2e1 + t247 / 0.2e1) * t288;
t187 = Icges(5,5) * t371 - Icges(5,6) * t375 - Icges(5,3) * t291;
t352 = -t289 * t187 - t191 * t366;
t100 = -t189 * t374 - t352;
t243 = Icges(5,5) * t290 - Icges(5,6) * t288;
t381 = t243 * t291;
t188 = Icges(5,3) * t289 + t381;
t351 = t289 * t188 + t192 * t366;
t101 = -t190 * t374 + t351;
t164 = t192 * t371;
t328 = t190 * t288 - t187;
t330 = t291 * t188 - t164;
t383 = t189 * t288;
t99 = -t190 * t375 - t330;
t334 = ((t99 - t164 + (t188 + t383) * t291 + t352) * t291 + t351 * t289) * t408 + (-t100 * t291 + t101 * t289) * t451 + ((t328 * t289 + t100 + t330 + t99) * t289 + (t101 - t351 + (-t191 * t290 + t383) * t289 + (t328 + t187) * t291) * t291) * t409;
t441 = 4 * qJD(1);
t439 = 2 * qJD(3);
t432 = m(4) * t83;
t430 = m(4) * t94;
t85 = t155 * t215 - t156 * t216;
t88 = t160 * t215 - t161 * t216;
t422 = m(5) * (t88 + t85);
t421 = m(5) * ((-t156 + t161) * t291 + (t155 - t160) * t289) * t248;
t315 = (-t155 * t291 - t156 * t289) * t250;
t353 = -t215 * t446 + t216 * t447;
t420 = m(5) * (t315 + t353);
t313 = (t289 * t446 - t291 * t447) * t248;
t419 = m(5) * (t315 + t313);
t314 = (-t160 * t291 - t161 * t289) * t250;
t418 = m(5) * (t314 + t353);
t417 = m(5) * (t314 + t313);
t416 = m(5) * t77;
t414 = m(5) * t81;
t412 = m(5) * t85;
t411 = m(5) * t88;
t369 = t289 * t297;
t204 = Icges(4,4) * t369 - Icges(4,2) * t370 - Icges(4,6) * t291;
t382 = t204 * t295;
t258 = Icges(4,5) * t297 - Icges(4,6) * t295;
t378 = t258 * t291;
t361 = t291 * t297;
t202 = Icges(4,5) * t369 - Icges(4,6) * t370 - Icges(4,3) * t291;
t269 = Icges(4,4) * t370;
t206 = Icges(4,1) * t369 - Icges(4,5) * t291 - t269;
t349 = -t289 * t202 - t206 * t361;
t203 = Icges(4,3) * t289 + t378;
t207 = Icges(4,5) * t289 + t262 * t291;
t348 = t289 * t203 + t207 * t361;
t343 = t261 * t289 + t204;
t205 = Icges(4,6) * t289 + t260 * t291;
t342 = -t261 * t291 - t205;
t341 = -Icges(4,2) * t369 + t206 - t269;
t340 = -t259 * t291 + t207;
t335 = -t250 - t402;
t175 = t207 * t369;
t329 = t291 * t203 - t175;
t327 = t205 * t295 - t202;
t324 = t422 / 0.2e1 + t326;
t321 = Icges(4,5) * t295 + Icges(4,6) * t297;
t312 = (-t231 * t291 + t232 * t289) * t263;
t304 = (-t244 + t247) * t290 - t455 * t288;
t311 = -t334 + (t243 * t289 + t346 * t288 + t344 * t290 + t304 * t291) * t409 + (-t347 * t288 + t304 * t289 + t345 * t290 - t381) * t408;
t308 = -t326 + t449 * (t189 * t290 + t191 * t288);
t307 = t326 + t443;
t306 = t341 * t295 + t343 * t297;
t305 = -t340 * t295 + t342 * t297;
t303 = (-t259 + t262) * t297 - t454 * t295;
t302 = t307 + t444;
t301 = (-t215 * t291 + t216 * t289) * t248;
t300 = t308 - t443 + t449 * (t297 * t204 + t295 * t206);
t110 = -t204 * t362 - t349;
t111 = -t205 * t362 + t348;
t28 = (t327 * t291 + t111 - t348) * t291 + (t327 * t289 + t110 + t329) * t289;
t109 = -t205 * t370 - t329;
t29 = (t109 - t175 + (t203 + t382) * t291 + t349) * t291 + t348 * t289;
t75 = -(-(-t206 * t297 + t382) * t289 - t291 * t202) * t291 + t109 * t289;
t76 = -t110 * t291 + t111 * t289;
t299 = (t29 * t451 + t311 + (t28 + t75) * t410 + (t258 * t289 + t303 * t291 + t342 * t295 + t340 * t297) * t409 + (t303 * t289 - t343 * t295 + t341 * t297 - t378 + t76) * t408) * qJD(3);
t265 = -rSges(4,2) * t295 + t399;
t226 = t291 * t321;
t225 = t321 * t289;
t197 = t335 * t291;
t195 = t335 * t289;
t132 = -t337 * t403 + t153;
t70 = t326 + t411;
t68 = t326 + t412;
t65 = t417 / 0.2e1;
t61 = t418 / 0.2e1;
t54 = t419 / 0.2e1;
t53 = t420 / 0.2e1;
t51 = t421 / 0.2e1;
t35 = t307 + t453;
t34 = t307 + t414 + t430;
t33 = t407 + t416 + t432;
t19 = -t421 / 0.2e1 + t324;
t18 = t51 + t324;
t17 = m(5) * (t337 * t248 * t250 + t137 * t153) + t401;
t16 = t17 * qJD(4);
t14 = t51 - t422 / 0.2e1 + t308;
t13 = t302 + t400;
t12 = t302 - t400;
t10 = t334 * qJD(4);
t9 = t300 + t400 - t444;
t8 = t61 - t417 / 0.2e1 + t334;
t7 = t65 - t418 / 0.2e1 + t334;
t6 = t53 - t419 / 0.2e1 + t334;
t5 = t54 - t420 / 0.2e1 + t334;
t4 = t61 + t65 + t311;
t3 = t53 + t54 + t311;
t2 = (t76 / 0.2e1 - t29 / 0.2e1) * t291 + (t28 / 0.2e1 + t75 / 0.2e1) * t289 + t334;
t1 = t2 * qJD(3);
t11 = [qJD(2) * t33 + qJD(3) * t34 + qJD(4) * t68, t33 * qJD(1) + t13 * qJD(3) + t18 * qJD(4) + 0.2e1 * (t407 / 0.2e1 + t77 * t436 + t83 * t437) * qJD(2), t34 * qJD(1) + t13 * qJD(2) + t299 + t3 * qJD(4) + (((-t162 * t291 - t163 * t289) * t265 + t312) * t437 + (t155 * t197 + t156 * t195) * t436) * t439, t68 * qJD(1) + t18 * qJD(2) + t3 * qJD(3) + ((t315 + t301) * m(5) + t311) * qJD(4); t12 * qJD(3) + t19 * qJD(4) + (-t407 / 0.4e1 - t432 / 0.4e1 - t416 / 0.4e1) * t441, qJD(3) * t35 + qJD(4) * t70, t12 * qJD(1) + t35 * qJD(2) + t299 + t4 * qJD(4) + ((t160 * t197 + t161 * t195) * t436 + ((-t167 * t291 - t168 * t289) * t265 + t312) * t437) * t439, t19 * qJD(1) + t70 * qJD(2) + t4 * qJD(3) + ((t314 + t301) * m(5) + t311) * qJD(4); t300 * qJD(1) + t9 * qJD(2) + t1 + t6 * qJD(4) + (-t430 / 0.4e1 - t414 / 0.4e1) * t441, t9 * qJD(1) + t8 * qJD(4) + t1 + (t300 - t453) * qJD(2), (m(4) * ((t289 * (rSges(4,1) * t369 - t338) + t291 * (rSges(4,1) * t361 + t289 * rSges(4,3) - t271)) * (-t231 * t289 - t232 * t291) + t337 * t265 * t263) + (-t286 * t226 + (t306 * t291 + (t225 + t305) * t289) * t291) * t409 + (-t287 * t225 + (t305 * t289 + (t226 + t306) * t291) * t289) * t408 + m(5) * (t132 * t86 - t195 * t447 - t197 * t446) + t401) * qJD(3) + t448 + t445 * t2, t6 * qJD(1) + t8 * qJD(2) + t15 * qJD(3) + t448; (t308 - t412) * qJD(1) + t14 * qJD(2) + t5 * qJD(3) + t10, t14 * qJD(1) + (t308 - t411) * qJD(2) + t7 * qJD(3) + t10, t5 * qJD(1) + t7 * qJD(2) + ((t132 * t137 + (-t195 * t289 - t197 * t291) * t248) * m(5) + t401) * qJD(3) + t16, qJD(3) * t17 + t334 * t445 + t16;];
Cq = t11;
