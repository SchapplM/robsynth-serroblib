% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR9_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:09
% EndTime: 2019-12-31 16:56:20
% DurationCPUTime: 8.87s
% Computational Cost: add. (14311->488), mult. (35787->729), div. (0->0), fcn. (39289->6), ass. (0->281)
t284 = sin(qJ(1));
t283 = sin(qJ(3));
t285 = cos(qJ(4));
t287 = cos(qJ(1));
t365 = t285 * t287;
t282 = sin(qJ(4));
t370 = t284 * t282;
t237 = -t283 * t370 + t365;
t369 = t284 * t285;
t238 = t282 * t287 + t283 * t369;
t286 = cos(qJ(3));
t368 = t284 * t286;
t165 = Icges(5,5) * t238 + Icges(5,6) * t237 - Icges(5,3) * t368;
t385 = Icges(5,4) * t238;
t168 = Icges(5,2) * t237 - Icges(5,6) * t368 + t385;
t232 = Icges(5,4) * t237;
t171 = Icges(5,1) * t238 - Icges(5,5) * t368 + t232;
t371 = t283 * t287;
t239 = t282 * t371 + t369;
t240 = t283 * t365 - t370;
t364 = t286 * t287;
t91 = t165 * t364 + t239 * t168 - t240 * t171;
t459 = t284 * t91;
t458 = t287 * t91;
t167 = -Icges(5,5) * t240 + Icges(5,6) * t239 + Icges(5,3) * t364;
t377 = t167 * t283;
t234 = Icges(5,4) * t240;
t170 = Icges(5,2) * t239 + Icges(5,6) * t364 - t234;
t233 = Icges(5,4) * t239;
t172 = Icges(5,1) * t240 - Icges(5,5) * t364 - t233;
t455 = t170 * t282 + t172 * t285;
t108 = t286 * t455 - t377;
t310 = -t239 * t170 - t240 * t172;
t358 = t237 * t168 + t238 * t171;
t457 = t310 + t358 + (-t165 * t284 - t167 * t287) * t286;
t439 = t240 * rSges(5,1) - t239 * rSges(5,2);
t176 = rSges(5,3) * t364 - t439;
t323 = rSges(5,1) * t285 - rSges(5,2) * t282;
t228 = rSges(5,3) * t283 + t286 * t323;
t456 = -t176 * t283 + t228 * t364;
t90 = -t167 * t368 + t237 * t170 - t172 * t238;
t280 = t284 ^ 2;
t430 = m(5) / 0.2e1;
t313 = Icges(5,5) * t285 - Icges(5,6) * t282;
t216 = Icges(5,3) * t283 + t286 * t313;
t383 = Icges(5,4) * t285;
t316 = -Icges(5,2) * t282 + t383;
t220 = Icges(5,6) * t283 + t286 * t316;
t384 = Icges(5,4) * t282;
t318 = Icges(5,1) * t285 - t384;
t224 = Icges(5,5) * t283 + t286 * t318;
t120 = t216 * t364 + t239 * t220 - t240 * t224;
t453 = t120 * t283;
t387 = Icges(4,4) * t283;
t317 = Icges(4,2) * t286 + t387;
t221 = Icges(4,6) * t287 + t284 * t317;
t266 = Icges(4,4) * t368;
t372 = t283 * t284;
t225 = Icges(4,1) * t372 + Icges(4,5) * t287 + t266;
t306 = -t221 * t286 - t225 * t283;
t451 = t287 * t306;
t272 = pkin(6) * t364;
t274 = t287 * qJ(2);
t425 = pkin(1) + pkin(5);
t326 = -t425 * t284 + t274;
t391 = rSges(5,3) * t286;
t400 = pkin(3) * t283;
t450 = (-t391 + t400) * t287 - t272 + t326 + t439;
t222 = -Icges(4,6) * t284 + t287 * t317;
t386 = Icges(4,4) * t286;
t319 = Icges(4,1) * t283 + t386;
t226 = -Icges(4,5) * t284 + t287 * t319;
t245 = -Icges(4,2) * t372 + t266;
t256 = -Icges(4,2) * t283 + t386;
t246 = t256 * t287;
t258 = Icges(4,1) * t286 - t387;
t248 = t258 * t284;
t249 = t258 * t287;
t449 = -t283 * ((t222 - t249) * t284 - (t221 - t248) * t287) + t286 * ((t226 + t246) * t284 - (t225 + t245) * t287);
t415 = -t284 / 0.2e1;
t414 = t284 / 0.2e1;
t412 = t287 / 0.2e1;
t446 = t287 * t90;
t445 = t90 * t284;
t444 = t228 * t287;
t215 = Icges(5,3) * t286 - t283 * t313;
t366 = t285 * t224;
t374 = t282 * t220;
t443 = t283 * (t366 / 0.2e1 - t374 / 0.2e1 + t258 / 0.2e1 - t317 / 0.2e1 - t215 / 0.2e1);
t442 = (t222 * t286 + t226 * t283) * t287;
t174 = t238 * rSges(5,1) + t237 * rSges(5,2) - rSges(5,3) * t368;
t440 = -pkin(6) * t368 + t174;
t180 = ((rSges(5,3) + pkin(6)) * t283 + (pkin(3) + t323) * t286) * t287;
t262 = rSges(4,1) * t286 - rSges(4,2) * t283;
t250 = t262 * t284;
t252 = t262 * t287;
t204 = rSges(5,3) * t372 + t323 * t368;
t335 = -pkin(3) * t368 - pkin(6) * t372 - t204;
t361 = (t284 * t180 + t287 * t335) * t430 + m(4) * (-t250 * t287 + t284 * t252) / 0.2e1;
t307 = t366 - t374;
t296 = t215 - t307;
t376 = t216 * t283;
t437 = t286 * t296 - t376;
t297 = -t216 * t287 + t455;
t436 = t286 * t297 - t377;
t311 = -t168 * t282 + t171 * t285;
t298 = t216 * t284 - t311;
t379 = t165 * t283;
t435 = t286 * t298 - t379;
t244 = (-Icges(5,2) * t285 - t384) * t286;
t247 = (-Icges(5,1) * t282 - t383) * t286;
t434 = -(t224 / 0.2e1 + t244 / 0.2e1) * t282 + (t247 / 0.2e1 - t220 / 0.2e1) * t285;
t281 = t287 ^ 2;
t431 = 4 * qJD(1);
t118 = -t216 * t368 + t220 * t237 + t224 * t238;
t117 = t118 * t283;
t89 = -t165 * t368 + t358;
t322 = -t89 * t284 + t446;
t39 = t286 * t322 + t117;
t429 = -t39 / 0.2e1;
t200 = t220 * t284;
t202 = t224 * t284;
t77 = (-t200 * t282 + t202 * t285 + t165) * t286 + t298 * t283;
t428 = t77 / 0.2e1;
t201 = t220 * t287;
t203 = t224 * t287;
t78 = (t201 * t282 - t203 * t285 + t167) * t286 + t297 * t283;
t427 = t78 / 0.2e1;
t189 = Icges(5,5) * t239 + Icges(5,6) * t240;
t354 = Icges(5,2) * t240 - t172 + t233;
t356 = -Icges(5,1) * t239 + t170 - t234;
t83 = t189 * t283 + (-t354 * t282 - t356 * t285) * t286;
t426 = t83 / 0.2e1;
t227 = -t283 * t323 + t391;
t110 = (t227 * t284 + t174) * t286 + (-t228 * t284 + t204) * t283;
t111 = (t227 * t287 - t176) * t286;
t308 = t174 * t287 + t176 * t284;
t122 = t308 * t286;
t146 = t174 * t283 + t228 * t368;
t93 = (-t204 * t287 + t284 * t444) * t286 + t308 * t283;
t423 = m(5) * (t110 * t146 + t111 * t456 - t122 * t93);
t151 = t425 * t287 + (qJ(2) + t400) * t284 + t440;
t422 = m(5) * (t110 * t151 + t111 * t450 - t146 * t335 + t180 * t456);
t265 = pkin(3) * t286 + pkin(6) * t283;
t345 = t228 + t265;
t184 = t345 * t284;
t186 = t345 * t287;
t194 = rSges(5,1) * t237 - rSges(5,2) * t238;
t195 = rSges(5,1) * t239 + rSges(5,2) * t240;
t251 = (-rSges(5,1) * t282 - rSges(5,2) * t285) * t286;
t420 = m(5) * (-t184 * t195 - t186 * t194 + (-t151 * t287 + t284 * t450) * t251);
t418 = m(5) * (-t151 * t335 + t180 * t450);
t417 = -t108 / 0.2e1;
t416 = t283 / 0.2e1;
t413 = t284 / 0.4e1;
t411 = t287 / 0.4e1;
t410 = m(3) * ((rSges(3,3) * t287 + t274) * t287 + (rSges(3,3) + qJ(2)) * t280);
t325 = rSges(4,1) * t283 + rSges(4,2) * t286;
t293 = -t284 * rSges(4,3) + t287 * t325;
t196 = t293 + t326;
t197 = (rSges(4,3) + t425) * t287 + (qJ(2) + t325) * t284;
t409 = m(4) * (t196 * t252 + t197 * t250);
t408 = m(4) * (t196 * t287 + t197 * t284);
t406 = m(5) * (t146 * t284 + t287 * t456);
t405 = m(5) * (t151 * t284 + t287 * t450);
t403 = m(5) * (t184 * t287 - t186 * t284);
t134 = -t194 * t287 - t284 * t195;
t401 = m(5) * t134;
t92 = t167 * t364 - t310;
t399 = t92 + t457;
t396 = -t89 + t457;
t389 = t284 * t39;
t321 = t287 * t92 - t459;
t40 = t286 * t321 + t453;
t388 = t287 * t40;
t219 = Icges(5,6) * t286 - t283 * t316;
t375 = t282 * t219;
t241 = (-Icges(5,5) * t282 - Icges(5,6) * t285) * t286;
t373 = t283 * t241;
t223 = Icges(5,5) * t286 - t283 * t318;
t367 = t285 * t223;
t294 = m(5) * (-t110 * t287 + t111 * t284);
t343 = t280 + t281;
t303 = t343 * t251 * t430;
t56 = -t294 / 0.2e1 + t303;
t363 = t56 * qJD(2);
t357 = -Icges(5,1) * t237 + t168 + t385;
t355 = -Icges(5,2) * t238 + t171 + t232;
t352 = t220 - t247;
t349 = t224 + t244;
t346 = pkin(6) * t286 + t227 - t400;
t342 = qJD(1) * t286;
t341 = qJD(4) * t286;
t188 = Icges(5,5) * t237 - Icges(5,6) * t238;
t71 = -t188 * t368 + t355 * t237 - t357 * t238;
t72 = -t189 * t368 + t354 * t237 - t356 * t238;
t32 = t72 * t284 + t287 * t71;
t66 = t200 * t237 + t202 * t238 - t284 * t435;
t67 = -t201 * t237 - t203 * t238 - t284 * t436;
t87 = t219 * t237 + t223 * t238 - t284 * t437;
t8 = (-t284 * t66 + t287 * t67 + t118) * t286 + (-t322 + t87) * t283;
t340 = -t32 / 0.2e1 + t8 / 0.2e1;
t73 = t188 * t364 + t355 * t239 + t357 * t240;
t74 = t189 * t364 + t354 * t239 + t356 * t240;
t33 = t74 * t284 + t287 * t73;
t68 = t239 * t200 - t240 * t202 + t287 * t435;
t69 = -t239 * t201 + t240 * t203 + t287 * t436;
t88 = t239 * t219 - t240 * t223 + t287 * t437;
t9 = (-t284 * t68 + t287 * t69 + t120) * t286 + (-t321 + t88) * t283;
t339 = t33 / 0.2e1 - t9 / 0.2e1;
t13 = -t453 + (t396 * t287 + t459) * t286;
t337 = -t13 / 0.2e1 - t40 / 0.2e1;
t12 = t117 + (-t399 * t284 + t446) * t286;
t336 = t429 + t12 / 0.2e1;
t314 = Icges(4,5) * t283 + Icges(4,6) * t286;
t136 = t287 * (Icges(4,3) * t287 + t284 * t314) + t221 * t368 + t225 * t372;
t218 = -Icges(4,3) * t284 + t287 * t314;
t137 = -t287 * t218 - t222 * t368 - t226 * t372;
t332 = -t368 / 0.4e1;
t331 = t364 / 0.4e1;
t328 = t325 * t343;
t82 = t188 * t283 + (-t355 * t282 - t357 * t285) * t286;
t98 = t349 * t237 - t352 * t238 - t241 * t368;
t99 = t349 * t239 + t352 * t240 + t241 * t364;
t327 = t420 / 0.2e1 + (t83 + t99) * t413 + (t82 + t98) * t411;
t315 = Icges(4,5) * t286 - Icges(4,6) * t283;
t107 = t286 * t311 + t379;
t312 = -t107 * t284 - t108 * t287;
t304 = -m(5) * (t151 * t194 - t195 * t450) - t373 / 0.2e1;
t140 = -t284 * t218 + t442;
t18 = t396 * t284 - t458;
t46 = t92 * t284 + t458;
t302 = t280 * t218 / 0.2e1 + t140 * t414 + t18 / 0.2e1 + t46 / 0.2e1 + (t137 + (t218 - t306) * t287 + t451) * t412;
t17 = t399 * t287 + t445;
t45 = t287 * t89 + t445;
t301 = -t136 * t287 / 0.2e1 + t137 * t415 + (t137 - t451) * t414 + (t140 - t442 + (t218 + t306) * t284 + t136) * t412 - t45 / 0.2e1 + t17 / 0.2e1;
t295 = t108 / 0.2e1 + t417;
t292 = t12 * t413 + t13 * t411 + t17 * t331 - t389 / 0.4e1 + t388 / 0.4e1 - t45 * t364 / 0.4e1 + (t18 + t46) * t332;
t104 = (t216 + t367 - t375) * t286 + t296 * t283;
t138 = t286 * t307 + t376;
t291 = t104 * t416 + t138 * t286 / 0.2e1 + t422 / 0.2e1 + (t107 + t118) * t372 / 0.4e1 - (-t108 + t120) * t371 / 0.4e1 + (t77 + t87) * t332 + (t78 + t88) * t331;
t290 = -t367 / 0.2e1 + t375 / 0.2e1 + t319 / 0.2e1 + t256 / 0.2e1 - t216 / 0.2e1;
t243 = t315 * t287;
t242 = t284 * t315;
t185 = t346 * t287;
t183 = t346 * t284;
t153 = -t283 * t195 + t251 * t364;
t152 = t194 * t283 + t251 * t368;
t135 = -t284 * t194 + t195 * t287;
t130 = t401 / 0.2e1;
t129 = t134 * t286;
t127 = t403 / 0.2e1;
t121 = t335 * t284 + (-t265 * t287 - t444) * t287;
t114 = (t373 + (-t349 * t282 - t352 * t285) * t286) * t283;
t113 = (-pkin(3) * t371 + t176 + t272) * t287 + (-pkin(3) * t372 - t440) * t284;
t100 = t406 / 0.2e1;
t65 = t405 + t408 + t410;
t60 = -t403 / 0.2e1 + t361;
t59 = t127 + t361;
t58 = t127 - t361;
t55 = t294 / 0.2e1 + t303;
t54 = t113 * t135 + (t184 * t284 + t186 * t287) * t251;
t51 = t286 * t434 - t304;
t42 = t138 * t283 + t286 * t312;
t31 = -t286 * t290 + t409 + t418 - t443;
t30 = t100 - t401 / 0.2e1;
t29 = t130 + t100;
t28 = t130 - t406 / 0.2e1;
t27 = t69 * t284 + t287 * t68;
t26 = t67 * t284 + t287 * t66;
t22 = t99 * t283 + (-t284 * t73 + t287 * t74) * t286;
t21 = t98 * t283 + (-t284 * t71 + t287 * t72) * t286;
t14 = (-t77 * t284 + t78 * t287 + t138) * t286 + (t104 - t312) * t283;
t7 = m(5) * t54 + t32 * t412 + t33 * t414;
t6 = (t337 * t284 + t336 * t287) * t286;
t5 = t284 * t301 + t287 * t302;
t4 = t423 + (t8 * t415 + t9 * t412 + t42 / 0.2e1) * t286 + (t389 / 0.2e1 - t388 / 0.2e1 + t14 / 0.2e1) * t283;
t3 = t291 + (-t98 / 0.4e1 - t82 / 0.4e1) * t287 + (-t99 / 0.4e1 - t83 / 0.4e1) * t284 + t292 - t420 / 0.2e1;
t2 = (-t138 / 0.2e1 + (-t88 / 0.4e1 - t78 / 0.4e1) * t287 + (t87 / 0.4e1 + t77 / 0.4e1) * t284) * t286 + (-t104 / 0.2e1 + (t120 / 0.4e1 - t108 / 0.4e1) * t287 + (-t118 / 0.4e1 - t107 / 0.4e1) * t284) * t283 + t292 - t422 / 0.2e1 + t327;
t1 = t291 + (-t40 / 0.4e1 - t13 / 0.4e1 + (t45 / 0.4e1 - t17 / 0.4e1) * t286) * t287 + (t39 / 0.4e1 - t12 / 0.4e1 + (t46 / 0.4e1 + t18 / 0.4e1) * t286) * t284 + t327;
t10 = [t65 * qJD(2) + t31 * qJD(3) + t51 * qJD(4), qJD(1) * t65 + qJD(3) * t59 + qJD(4) * t29, t31 * qJD(1) + t59 * qJD(2) + t1 * qJD(4) + (m(5) * (-t151 * t185 + t180 * t184 + t183 * t450 + t186 * t335) + (m(4) * (t197 * t325 - t250 * t262) + t428 + t87 / 0.2e1 - t314 * t412 + (-t221 / 0.2e1 + t248 / 0.2e1) * t286 + (-t225 / 0.2e1 - t245 / 0.2e1) * t283 - t302) * t287 + (m(4) * (-t196 * t325 + t252 * t262) + t427 + t88 / 0.2e1 - t314 * t414 + (t222 / 0.2e1 - t249 / 0.2e1) * t286 + (t226 / 0.2e1 + t246 / 0.2e1) * t283 - t301) * t284) * qJD(3), t51 * qJD(1) + t29 * qJD(2) + t1 * qJD(3) + (t114 + m(5) * (t146 * t194 + t151 * t152 + t153 * t450 - t195 * t456)) * qJD(4) + ((t426 + t99 / 0.2e1 - t336) * t287 + (-t82 / 0.2e1 - t98 / 0.2e1 - t337) * t284) * t341; t60 * qJD(3) + t28 * qJD(4) + (-t410 / 0.4e1 - t408 / 0.4e1 - t405 / 0.4e1) * t431, 0, t60 * qJD(1) + 0.2e1 * ((t183 * t284 + t185 * t287) * t430 - m(4) * t328 / 0.2e1) * qJD(3) + t55 * qJD(4), t28 * qJD(1) + t55 * qJD(3) + m(5) * (-t152 * t287 + t153 * t284) * qJD(4); t58 * qJD(2) + t5 * qJD(3) + t2 * qJD(4) + (-t409 / 0.4e1 - t418 / 0.4e1) * t431 + t290 * t342 + (t295 * t287 + t443) * qJD(1), qJD(1) * t58 + qJD(4) * t56, t5 * qJD(1) + (m(4) * ((-t287 * t293 + (-t287 * rSges(4,3) - t284 * t325) * t284) * (-t284 * t250 - t252 * t287) - t262 * t328) + m(5) * (t113 * t121 + t183 * t184 + t185 * t186) + (-t280 * t243 + (t284 * t242 + t449) * t287 + t27) * t414 + (t281 * t242 + (-t287 * t243 - t449) * t284 + t26) * t412) * qJD(3) + t7 * qJD(4), t2 * qJD(1) + t363 + t7 * qJD(3) + (-t42 / 0.2e1 + t339 * t287 + t340 * t284) * t341 + (m(5) * (t129 * t113 - t122 * t135 - t152 * t186 + t153 * t184 + (-t146 * t287 + t284 * t456) * t251) + t21 * t412 + t22 * t414 - t423 + (-t14 / 0.2e1 + (t82 / 0.2e1 + t40 / 0.2e1) * t287 + (t426 + t429) * t284) * t283) * qJD(4); t304 * qJD(1) + t30 * qJD(2) + t3 * qJD(3) + t6 * qJD(4) + (-t295 * t284 - t434) * t342, qJD(1) * t30 - qJD(3) * t56, t3 * qJD(1) - t363 + (((t27 / 0.2e1 + t107 / 0.2e1) * t286 + t340) * t287 + ((-t26 / 0.2e1 + t417) * t286 - t339) * t284 + ((-t46 / 0.2e1 + t428) * t287 + (t45 / 0.2e1 + t427) * t284) * t283 + (-t110 * t186 + t111 * t184 + t113 * t93 - t121 * t122 - t146 * t185 + t183 * t456 - t54) * m(5)) * qJD(3) + t4 * qJD(4), t6 * qJD(1) + t4 * qJD(3) + (m(5) * (-t122 * t129 + t146 * t152 + t153 * t456) + t114 * t416 + (t21 * t415 + t22 * t412 + (-t284 * t82 + t287 * t83) * t416) * t286) * qJD(4);];
Cq = t10;
