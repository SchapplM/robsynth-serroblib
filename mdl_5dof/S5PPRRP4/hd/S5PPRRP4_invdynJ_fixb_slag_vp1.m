% Calculate vector of inverse dynamics joint torques for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:27
% EndTime: 2019-12-31 17:34:40
% DurationCPUTime: 11.34s
% Computational Cost: add. (6232->437), mult. (14242->540), div. (0->0), fcn. (15263->6), ass. (0->227)
t211 = sin(pkin(7));
t212 = cos(pkin(7));
t370 = sin(qJ(3));
t371 = cos(qJ(3));
t176 = -t211 * t371 + t212 * t370;
t228 = t211 * t370 + t212 * t371;
t440 = Icges(5,3) + Icges(6,3);
t214 = sin(qJ(4));
t215 = cos(qJ(4));
t243 = Icges(6,5) * t215 - Icges(6,6) * t214;
t245 = Icges(5,5) * t215 - Icges(5,6) * t214;
t460 = t243 + t245;
t400 = t440 * t176 + t460 * t228;
t339 = Icges(6,4) * t215;
t247 = -Icges(6,2) * t214 + t339;
t96 = Icges(6,6) * t176 + t228 * t247;
t341 = Icges(5,4) * t215;
t249 = -Icges(5,2) * t214 + t341;
t99 = Icges(5,6) * t176 + t228 * t249;
t421 = t96 + t99;
t340 = Icges(6,4) * t214;
t251 = Icges(6,1) * t215 - t340;
t102 = Icges(6,5) * t176 + t228 * t251;
t342 = Icges(5,4) * t214;
t253 = Icges(5,1) * t215 - t342;
t105 = Icges(5,5) * t176 + t228 * t253;
t449 = t102 + t105;
t425 = t214 * t421 - t215 * t449;
t405 = t400 * t176 - t425 * t228;
t432 = -t460 * t176 + t440 * t228;
t461 = t432 * t176;
t250 = Icges(6,1) * t214 + t339;
t252 = Icges(5,1) * t214 + t341;
t455 = t250 + t252;
t246 = Icges(6,2) * t215 + t340;
t248 = Icges(5,2) * t215 + t342;
t456 = -t246 - t248;
t391 = t456 * t214 + t455 * t215;
t244 = Icges(5,5) * t214 + Icges(5,6) * t215;
t419 = t244 * t176;
t242 = Icges(6,5) * t214 + Icges(6,6) * t215;
t420 = t242 * t176;
t401 = t228 * t391 + t419 + t420;
t459 = qJD(3) * t401;
t403 = t449 * t214 + t215 * t421;
t101 = Icges(6,5) * t228 - t176 * t251;
t104 = Icges(5,5) * t228 - t176 * t253;
t450 = t101 + t104;
t95 = Icges(6,6) * t228 - t176 * t247;
t98 = Icges(5,6) * t228 - t176 * t249;
t454 = t95 + t98;
t453 = Icges(5,1) + Icges(6,1);
t452 = Icges(5,4) + Icges(6,4);
t442 = Icges(5,5) + Icges(6,5);
t451 = Icges(5,2) + Icges(6,2);
t441 = Icges(5,6) + Icges(6,6);
t329 = t228 * t215;
t448 = t450 * t329 + t461;
t325 = t176 * t215;
t386 = t432 * t228 - t450 * t325;
t326 = t176 * t214;
t408 = t454 * t326 + t386;
t407 = t400 * t228 - t449 * t325 + t421 * t326;
t330 = t228 * t214;
t406 = -t454 * t330 + t448;
t165 = t228 * qJD(3);
t164 = t176 * qJD(3);
t301 = qJD(4) * t214;
t291 = t228 * t301;
t232 = -t164 * t215 - t291;
t300 = qJD(4) * t215;
t233 = t164 * t214 - t228 * t300;
t447 = t441 * t165 + t452 * t232 + t451 * t233;
t290 = t176 * t301;
t331 = t165 * t215;
t230 = t290 - t331;
t332 = t165 * t214;
t231 = t176 * t300 + t332;
t446 = -t441 * t164 + t452 * t230 + t451 * t231;
t445 = -t442 * t165 - t453 * t232 - t452 * t233;
t444 = -t442 * t164 + t453 * t230 + t452 * t231;
t431 = rSges(6,3) + qJ(5) + pkin(6);
t418 = t176 * t431;
t439 = (t247 + t249) * qJD(4);
t438 = (t251 + t253) * qJD(4);
t437 = -t455 * t214 + t456 * t215;
t404 = t450 * t214 + t454 * t215;
t116 = t242 * t228;
t118 = t244 * t228;
t402 = -t391 * t176 + t116 + t118;
t436 = t454 * t214;
t434 = t165 * t440 + t232 * t442 + t233 * t441;
t433 = -t164 * t440 + t230 * t442 + t231 * t441;
t430 = t460 * qJD(4);
t429 = qJD(4) * t437 - t214 * t439 + t215 * t438;
t428 = -t404 * qJD(4) - t446 * t214 + t444 * t215;
t427 = t403 * qJD(4) + t447 * t214 + t445 * t215;
t426 = t242 + t244;
t424 = -t450 * t215 + t436;
t423 = t405 * t176 + t406 * t228;
t422 = t407 * t176 + t408 * t228;
t372 = rSges(6,1) + pkin(4);
t62 = rSges(5,1) * t230 + rSges(5,2) * t231 - t164 * rSges(5,3);
t417 = t402 * qJD(3);
t139 = -rSges(4,1) * t228 + rSges(4,2) * t176;
t416 = t139 * qJD(3);
t135 = -qJD(4) * t164 + qJDD(4) * t228;
t352 = t215 * rSges(6,2);
t268 = rSges(6,1) * t214 + t352;
t208 = qJDD(2) * t211;
t281 = t176 * pkin(3) - pkin(6) * t228;
t283 = -t165 * pkin(3) - pkin(6) * t164;
t294 = qJD(3) * t283 - qJDD(3) * t281 + t208;
t210 = t214 * rSges(6,2);
t392 = -rSges(6,1) * t215 + t210;
t184 = t392 * qJD(4);
t304 = qJD(4) * t184;
t324 = t215 * qJD(4) ^ 2;
t369 = pkin(4) * t215;
t207 = pkin(3) + t369;
t394 = -rSges(6,1) * t325 + rSges(6,2) * t326 - t176 * t207 + t228 * t431;
t348 = t281 + t394;
t385 = t228 * qJD(5);
t387 = -rSges(6,1) * t331 - t164 * t431 - t165 * t207;
t364 = rSges(6,2) * t231 + t290 * t372 - t283 + t385 + t387;
t10 = t228 * t304 + qJD(5) * t165 + qJDD(5) * t176 - t135 * t268 + t348 * qJDD(3) + t364 * qJD(3) + (-t135 * t214 - t228 * t324) * pkin(4) + t294;
t415 = -g(1) + t10;
t414 = qJD(4) * t422 + t417;
t413 = qJD(4) * t423 + t459;
t412 = -t424 * qJD(4) + t444 * t214 + t446 * t215;
t411 = t425 * qJD(4) + t445 * t214 - t447 * t215;
t410 = -t164 * t391 + t165 * t426 + t176 * t430 + t228 * t429;
t409 = -t164 * t426 - t165 * t391 - t176 * t429 + t228 * t430;
t270 = rSges(5,1) * t214 + rSges(5,2) * t215;
t302 = qJD(4) * t270;
t397 = t176 * t302;
t395 = -rSges(6,1) * t329 + rSges(6,2) * t330 - t207 * t228 - t418;
t173 = t176 * pkin(6);
t393 = pkin(3) * t228 + t173;
t313 = pkin(4) * t291 - qJD(5) * t176;
t390 = -t165 * t431 + t313;
t389 = (-t432 * t164 + t424 * t165 + t433 * t228) * t228 + (t427 * t176 + t425 * t165 - t400 * t164 + (-t428 + t434) * t228) * t176;
t388 = (t424 * t164 + t432 * t165 + t428 * t228) * t228 + (t434 * t176 + t400 * t165 + t425 * t164 + (-t427 + t433) * t228) * t176;
t306 = qJD(4) * t228;
t305 = qJD(4) * t176;
t138 = qJD(3) * t281;
t384 = t348 * qJD(3) - t138 - t313;
t319 = t248 * t176 + t104;
t344 = t252 * t176 - t98;
t383 = t214 * t319 - t215 * t344;
t321 = t246 * t176 + t101;
t346 = t250 * t176 - t95;
t382 = t214 * t321 - t215 * t346;
t381 = m(3) + m(4);
t366 = pkin(3) - t207;
t160 = t165 * pkin(6);
t365 = rSges(6,1) * t232 + rSges(6,2) * t233 + t164 * t366 - t160 - t390;
t358 = rSges(5,3) * t165;
t284 = pkin(4) * t214 + t268;
t347 = -t173 + (-t392 - t366) * t228 + t418;
t296 = -t393 - t347;
t307 = qJD(2) * t212;
t29 = qJD(3) * t296 + t284 * t305 - t307 + t385;
t356 = t164 * t29;
t169 = t176 * rSges(5,3);
t349 = t393 + t395;
t345 = t228 * t250 + t96;
t343 = t228 * t252 + t99;
t320 = -t228 * t246 + t102;
t318 = -t228 * t248 + t105;
t202 = -rSges(5,1) * t215 + rSges(5,2) * t214;
t111 = -t202 * t228 + t169;
t317 = -t111 - t393;
t311 = -t246 + t251;
t310 = -t247 - t250;
t309 = -t248 + t253;
t308 = -t249 - t252;
t185 = t202 * qJD(4);
t303 = qJD(4) * t185;
t299 = t243 * qJD(3);
t298 = t245 * qJD(3);
t297 = qJDD(2) * t212;
t295 = -m(5) - m(6) - t381;
t287 = -t306 / 0.2e1;
t286 = -t305 / 0.2e1;
t285 = t305 / 0.2e1;
t109 = -rSges(5,1) * t325 + rSges(5,2) * t326 + rSges(5,3) * t228;
t280 = -qJD(3) * t109 + t138;
t279 = pkin(4) * t326 + t268 * t176;
t278 = -pkin(4) * t330 - t228 * t268;
t209 = qJD(2) * t211;
t28 = -t268 * t306 + t209 + t384;
t271 = t28 * t284;
t39 = -t228 * t302 + t209 - t280;
t40 = qJD(3) * t317 - t307 + t397;
t263 = -t176 * t40 + t228 * t39;
t60 = rSges(5,1) * t232 + rSges(5,2) * t233 + t358;
t262 = t176 * t62 - t228 * t60;
t241 = t109 * t176 - t111 * t228;
t107 = -rSges(5,1) * t329 + rSges(5,2) * t330 - t169;
t236 = pkin(3) - t202;
t235 = t207 - t392;
t234 = t214 * t372 + t352;
t227 = t234 * t176;
t226 = t214 * t320 + t215 * t345;
t225 = t214 * t318 + t215 * t343;
t224 = (t214 * t310 + t215 * t311) * qJD(3);
t223 = (t214 * t308 + t215 * t309) * qJD(3);
t140 = -rSges(4,1) * t176 - rSges(4,2) * t228;
t137 = qJD(3) * t393;
t134 = qJD(4) * t165 + qJDD(4) * t176;
t133 = -pkin(3) * t164 + t160;
t132 = t270 * t228;
t130 = t270 * t176;
t128 = -rSges(4,1) * t165 + rSges(4,2) * t164;
t127 = -rSges(4,1) * t164 - rSges(4,2) * t165;
t113 = -t307 + t416;
t66 = -qJD(3) * t127 + qJDD(3) * t139 - t297;
t65 = qJD(3) * t128 + qJDD(3) * t140 + t208;
t38 = qJD(4) * t241 + qJD(1);
t23 = qJD(1) + (t176 * t348 - t228 * t347) * qJD(4);
t22 = -t176 * t303 - t297 + t134 * t270 + t317 * qJDD(3) + (-t133 - t60) * qJD(3);
t21 = qJD(3) * t62 + qJDD(3) * t109 - t135 * t270 + t228 * t303 + t294;
t16 = qJD(4) * t262 + t109 * t134 - t111 * t135 + qJDD(1);
t11 = -t176 * t304 - qJD(5) * t164 - t297 + qJDD(5) * t228 + t134 * t268 + (t134 * t214 + t176 * t324) * pkin(4) + t296 * qJDD(3) + (-t133 - t365) * qJD(3);
t1 = qJDD(1) - t347 * t135 + t348 * t134 + (t364 * t176 - t228 * t365) * qJD(4);
t2 = [m(5) * t16 + m(6) * t1 + (m(2) + t381) * qJDD(1) + (-m(2) + t295) * g(3); t295 * (g(1) * t211 - g(2) * t212) + m(4) * (t211 * t65 - t212 * t66) + m(5) * (t21 * t211 - t212 * t22) + m(6) * (t10 * t211 - t11 * t212) + m(3) * (t211 ^ 2 + t212 ^ 2) * qJDD(2); t414 * t285 + (t391 * qJD(4) + t214 * t438 + t215 * t439) * qJD(3) + (Icges(4,3) - t437) * qJDD(3) + (-g(2) * t395 + t235 * t356 + t415 * t394 + (-t235 * t228 - t418) * t11 + (t384 + t390) * t29 + (rSges(6,2) * t332 - t349 * qJD(3) + t137 + t387) * t28 + (t227 * t28 - (t271 + t23 * (t347 + t349)) * t176) * qJD(4)) * m(6) + (-g(2) * (t107 - t393) - t38 * (t107 + t111) * t305 + (t21 - g(1)) * (t109 - t281) + (-qJD(3) * t107 + t137 + t283 - t397 + t62) * t39 + (-t228 * t236 - t169 - t173) * t22 + (t236 * t164 - t160 - t280 - t358) * t40) * m(5) + (-t113 * t127 + (t128 - t416) * (qJD(3) * t140 + t209) + (t113 * qJD(3) - g(1) + t65) * t140 + (-g(2) + t66) * t139) * m(4) - (-t401 - t403) * t134 / 0.2e1 - (-t402 - t404) * t135 / 0.2e1 + (((-t406 + t407 + t448) * t176 + t386 * t228) * qJD(4) - t410 + t411 + t417) * t286 + ((((-t400 + t436) * t176 + t386 - t408) * t176 - (-(-t400 + t424) * t228 + t461 - t407) * t228) * qJD(4) - t409 - t412 + t413 - t459) * t287; t423 * t134 / 0.2e1 + t422 * t135 / 0.2e1 - t414 * t164 / 0.2e1 + t413 * t165 / 0.2e1 + (t409 * qJD(3) + t389 * qJD(4) + t402 * qJDD(3) + t407 * t134 + t408 * t135) * t228 / 0.2e1 + (t410 * qJD(3) + t388 * qJD(4) + t401 * qJDD(3) + t405 * t134 + t406 * t135) * t176 / 0.2e1 - (t404 * t164 - t403 * t165 + t411 * t176 - t228 * t412) * qJD(3) / 0.2e1 + ((((-t318 - t320) * t176 - (t319 + t321) * t228) * t215 + ((t343 + t345) * t176 - (t344 + t346) * t228) * t214) * qJD(4) + ((t308 + t310) * t215 + (-t309 - t311) * t214) * qJD(3)) * qJD(3) / 0.2e1 - (-t403 * t176 - t228 * t404) * qJDD(3) / 0.2e1 + (-t164 * t408 + t165 * t407 + t389) * t306 / 0.2e1 + (-(-t306 * t419 - t298) * t228 + (-t223 + (t225 * t176 - (t118 - t383) * t228) * qJD(4)) * t176 - (-t306 * t420 - t299) * t228 + (-t224 + (t226 * t176 - (t116 - t382) * t228) * qJD(4)) * t176) * t287 + ((-t118 * t305 + t298) * t176 - (-t223 + (t383 * t228 + (-t419 + t225) * t176) * qJD(4)) * t228 + (-t116 * t305 + t299) * t176 - (-t224 + (t382 * t228 + (-t420 + t226) * t176) * qJD(4)) * t228) * t286 + (-t164 * t406 + t165 * t405 + t388) * t285 + (g(1) * t132 - g(2) * t130 - g(3) * t202 - (t130 * t39 + t132 * t40) * qJD(3) - (t38 * (t130 * t176 + t132 * t228) + t263 * t202) * qJD(4) + t16 * t241 + t38 * (t109 * t165 + t111 * t164 + t262) + t263 * t185 - (-t164 * t39 - t165 * t40 - t176 * t22 + t21 * t228) * t270) * m(5) + (-(-t278 * t29 + t279 * t28) * qJD(3) - g(3) * (-t215 * t372 + t210) - g(2) * t227 + (t23 * t348 + t284 * t29) * t165 + (t23 * t347 + t271) * t164 + (-(t23 * t279 - t29 * t392) * qJD(4) + t1 * t348 + t11 * t284 - t29 * t184 + t23 * t364) * t176 - (-(t28 * (-t392 + t369) + t278 * t23) * qJD(4) - g(1) * t234 + t10 * t284 + t28 * (pkin(4) * t300 - t184) + t1 * t347 + t23 * t365) * t228) * m(6); (t165 * t28 - t356 + (qJD(3) * t29 + t415) * t176 - (t28 * qJD(3) + g(2) - t11) * t228) * m(6);];
tau = t2;
