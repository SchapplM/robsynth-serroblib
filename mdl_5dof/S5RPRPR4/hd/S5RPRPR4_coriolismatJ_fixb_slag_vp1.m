% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:25
% EndTime: 2022-01-23 09:22:36
% DurationCPUTime: 5.95s
% Computational Cost: add. (22593->391), mult. (14687->521), div. (0->0), fcn. (13480->10), ass. (0->256)
t297 = qJ(1) + pkin(8);
t291 = cos(t297);
t296 = qJ(3) + pkin(9);
t292 = qJ(5) + t296;
t283 = sin(t292);
t284 = cos(t292);
t227 = rSges(6,1) * t283 + rSges(6,2) * t284;
t288 = sin(t296);
t299 = sin(qJ(3));
t429 = pkin(3) * t299;
t247 = -pkin(4) * t288 - t429;
t316 = t227 - t247;
t129 = t316 * t291;
t289 = sin(t297);
t475 = t316 * t289;
t484 = t289 * t475;
t473 = t129 * t291 + t484;
t290 = cos(t296);
t309 = rSges(5,1) * t288 + rSges(5,2) * t290 + t429;
t476 = t309 * t291;
t477 = t309 * t289;
t480 = t289 * t477 + t291 * t476;
t423 = -m(6) * t473 / 0.2e1 - m(5) * t480 / 0.2e1;
t194 = t227 * t291;
t221 = t291 * t247;
t132 = t221 - t194;
t457 = m(6) / 0.2e1;
t458 = m(5) / 0.2e1;
t424 = (-t132 * t291 + t484) * t457 + t480 * t458;
t25 = t424 - t423;
t490 = t25 * qJD(1);
t286 = t289 ^ 2;
t287 = t291 ^ 2;
t349 = t286 + t287;
t489 = t227 * t349;
t267 = Icges(6,4) * t284;
t224 = -Icges(6,2) * t283 + t267;
t225 = Icges(6,1) * t283 + t267;
t488 = t224 + t225;
t301 = cos(qJ(3));
t487 = -Icges(4,5) * t299 - Icges(5,5) * t288 - Icges(4,6) * t301 - Icges(5,6) * t290;
t387 = t284 * t291;
t391 = t283 * t291;
t317 = rSges(6,1) * t387 - rSges(6,2) * t391 + t289 * rSges(6,3);
t430 = cos(qJ(1)) * pkin(1);
t294 = t301 * pkin(3);
t285 = t294 + pkin(2);
t428 = pkin(4) * t290;
t241 = t285 + t428;
t298 = -qJ(4) - pkin(6);
t295 = pkin(7) - t298;
t482 = t291 * t241 + t295 * t289;
t112 = -t317 - t430 - t482;
t264 = t295 * t291;
t392 = t283 * t289;
t352 = rSges(6,2) * t392 + t291 * rSges(6,3);
t419 = rSges(6,1) * t284;
t431 = sin(qJ(1)) * pkin(1);
t468 = t264 + (-t241 - t419) * t289 - t431 + t352;
t474 = t112 * t289 - t291 * t468;
t439 = t289 / 0.2e1;
t437 = -t291 / 0.2e1;
t485 = t291 / 0.2e1;
t386 = t288 * t289;
t384 = t289 * t290;
t409 = Icges(5,4) * t288;
t230 = Icges(5,2) * t290 + t409;
t233 = Icges(5,1) * t290 - t409;
t410 = Icges(4,4) * t299;
t250 = Icges(4,2) * t301 + t410;
t253 = Icges(4,1) * t301 - t410;
t481 = -(t233 / 0.2e1 - t230 / 0.2e1) * t288 - (t253 / 0.2e1 - t250 / 0.2e1) * t299;
t351 = t289 * t285 + t291 * t298;
t467 = t264 + t351;
t122 = t241 * t289 - t467;
t265 = t289 * t298;
t469 = t265 + t482;
t123 = -t291 * t285 + t469;
t33 = (t122 + t467) * t291 + (-t123 + (-t241 - t285) * t291 + t469) * t289;
t479 = m(6) * t33;
t228 = -rSges(6,2) * t283 + t419;
t323 = Icges(6,5) * t283 + Icges(6,6) * t284;
t181 = t323 * t289;
t182 = t291 * t323;
t408 = Icges(6,4) * t283;
t226 = Icges(6,1) * t284 - t408;
t162 = Icges(6,5) * t289 + t226 * t291;
t223 = Icges(6,2) * t284 + t408;
t362 = -t223 * t291 + t162;
t238 = Icges(6,4) * t392;
t388 = t284 * t289;
t161 = Icges(6,1) * t388 - Icges(6,5) * t291 - t238;
t363 = -Icges(6,2) * t388 + t161 - t238;
t160 = Icges(6,6) * t289 + t224 * t291;
t364 = -t225 * t291 - t160;
t159 = Icges(6,4) * t388 - Icges(6,2) * t392 - Icges(6,6) * t291;
t365 = t225 * t289 + t159;
t462 = (-t362 * t289 + t363 * t291) * t283 + (t364 * t289 + t365 * t291) * t284;
t426 = (-t286 * t182 + (t289 * t181 + t462) * t291) * t439 + (-t181 * t287 + (t182 * t291 + t462) * t289) * t437;
t282 = t291 * pkin(6);
t368 = -t289 * (pkin(2) * t289 - t282 - t351) + t291 * (-t289 * pkin(6) - t265 + (-pkin(2) + t285) * t291);
t96 = t289 * (rSges(6,1) * t388 - t352) + t291 * t317;
t47 = t122 * t289 + t123 * t291 + t368 + t96;
t193 = t227 * t289;
t472 = t289 * t193 + t291 * t194;
t6 = t426 + m(6) * (t228 * t473 - t47 * t472);
t478 = t6 * qJD(5);
t471 = t487 * t289;
t470 = t487 * t291;
t277 = Icges(5,4) * t290;
t231 = -Icges(5,2) * t288 + t277;
t232 = Icges(5,1) * t288 + t277;
t293 = Icges(4,4) * t301;
t251 = -Icges(4,2) * t299 + t293;
t252 = Icges(4,1) * t299 + t293;
t175 = Icges(5,5) * t289 + t233 * t291;
t358 = -t230 * t291 + t175;
t173 = Icges(5,6) * t289 + t231 * t291;
t360 = -t232 * t291 - t173;
t464 = -t358 * t386 + t360 * t384;
t383 = t289 * t299;
t260 = Icges(4,4) * t383;
t382 = t289 * t301;
t191 = Icges(4,1) * t382 - Icges(4,5) * t291 - t260;
t355 = -Icges(4,2) * t382 + t191 - t260;
t189 = Icges(4,4) * t382 - Icges(4,2) * t383 - Icges(4,6) * t291;
t357 = t252 * t289 + t189;
t244 = Icges(5,4) * t386;
t174 = Icges(5,1) * t384 - Icges(5,5) * t291 - t244;
t359 = -Icges(5,2) * t384 + t174 - t244;
t172 = Icges(5,4) * t384 - Icges(5,2) * t386 - Icges(5,6) * t291;
t361 = t232 * t289 + t172;
t463 = t359 * t288 + t361 * t290 + t355 * t299 + t357 * t301;
t333 = t488 * t284 / 0.2e1 + (-t223 / 0.2e1 + t226 / 0.2e1) * t283;
t124 = t162 * t388;
t157 = Icges(6,5) * t388 - Icges(6,6) * t392 - Icges(6,3) * t291;
t222 = Icges(6,5) * t284 - Icges(6,6) * t283;
t397 = t222 * t291;
t158 = Icges(6,3) * t289 + t397;
t336 = t160 * t283 - t157;
t339 = t158 * t291 - t124;
t372 = t289 * t158 + t162 * t387;
t373 = -t289 * t157 - t161 * t387;
t402 = t159 * t283;
t66 = -t160 * t392 - t339;
t67 = -t159 * t391 - t373;
t68 = -t160 * t391 + t372;
t345 = ((t66 - t124 + (t158 + t402) * t291 + t373) * t291 + t372 * t289) * t437 + (t289 * t68 - t291 * t67) * t485 + ((t336 * t289 + t339 + t66 + t67) * t289 + (-t372 + t68 + (-t161 * t284 + t402) * t289 + (t336 + t157) * t291) * t291) * t439;
t461 = 0.4e1 * qJD(1);
t460 = 2 * qJD(3);
t459 = m(4) / 0.2e1;
t421 = rSges(4,1) * t301;
t347 = pkin(2) + t421;
t350 = rSges(4,2) * t383 + t291 * rSges(4,3);
t133 = -t347 * t289 + t282 + t350 - t431;
t378 = t291 * t299;
t262 = rSges(4,2) * t378;
t134 = t430 - t262 + t347 * t291 + (rSges(4,3) + pkin(6)) * t289;
t254 = rSges(4,1) * t299 + rSges(4,2) * t301;
t216 = t254 * t289;
t217 = t254 * t291;
t456 = m(4) * (t133 * t216 - t134 * t217);
t318 = rSges(5,1) * t384 - rSges(5,2) * t386 - t291 * rSges(5,3);
t116 = -t318 - t351 - t431;
t385 = t288 * t291;
t344 = -rSges(5,2) * t385 + t289 * rSges(5,3);
t420 = rSges(5,1) * t290;
t117 = t430 - t265 + (t285 + t420) * t291 + t344;
t454 = m(5) * (t116 * t477 - t117 * t476);
t453 = m(5) * (t116 * t291 + t117 * t289);
t451 = t96 * t479;
t450 = t47 * t479;
t311 = t474 * t228;
t448 = m(6) * (-t129 * t193 + t194 * t475 + t311);
t447 = m(6) * (t311 + (-t132 * t289 - t291 * t475) * t227);
t446 = m(6) * (-t112 * t132 + t468 * t475);
t445 = m(6) * (t112 * t194 + t193 * t468);
t444 = m(6) * t474;
t440 = -t289 / 0.2e1;
t432 = m(6) * t472;
t400 = t172 * t288;
t398 = t189 * t299;
t24 = t33 * t457;
t395 = t24 * qJD(2);
t381 = t290 * t291;
t377 = t291 * t301;
t313 = t472 * t457;
t331 = m(6) * t489;
t76 = t313 + t331 / 0.2e1;
t376 = t76 * qJD(1);
t170 = Icges(5,5) * t384 - Icges(5,6) * t386 - Icges(5,3) * t291;
t371 = -t289 * t170 - t174 * t381;
t325 = Icges(5,5) * t290 - Icges(5,6) * t288;
t171 = Icges(5,3) * t289 + t325 * t291;
t370 = t289 * t171 + t175 * t381;
t187 = Icges(4,5) * t382 - Icges(4,6) * t383 - Icges(4,3) * t291;
t367 = -t289 * t187 - t191 * t377;
t327 = Icges(4,5) * t301 - Icges(4,6) * t299;
t188 = Icges(4,3) * t289 + t327 * t291;
t192 = Icges(4,5) * t289 + t253 * t291;
t366 = t289 * t188 + t192 * t377;
t190 = Icges(4,6) * t289 + t251 * t291;
t356 = -t252 * t291 - t190;
t354 = -t250 * t291 + t192;
t346 = rSges(5,2) * t288 - t294 - t420;
t135 = t175 * t384;
t338 = t171 * t291 - t135;
t153 = t192 * t382;
t337 = t188 * t291 - t153;
t335 = t173 * t288 - t170;
t334 = t190 * t299 - t187;
t332 = t451 / 0.2e1 + t345;
t315 = -t228 - t294 - t428;
t312 = t24 * qJD(3);
t305 = (-t223 + t226) * t284 - t488 * t283;
t310 = -t345 + (t222 * t289 + t364 * t283 + t362 * t284 + t305 * t291) * t439 + (-t365 * t283 + t363 * t284 + t305 * t289 - t397) * t437;
t308 = -t333 + (t439 + t440) * (t159 * t284 + t161 * t283);
t306 = -t354 * t299 + t356 * t301;
t256 = -rSges(4,2) * t299 + t421;
t180 = t346 * t291;
t178 = t346 * t289;
t130 = t315 * t291;
t128 = t315 * t289;
t119 = -t216 * t289 - t217 * t291;
t102 = qJD(5) * t432;
t94 = t309 * t349;
t83 = -t190 * t378 + t366;
t82 = -t189 * t378 - t367;
t81 = -t190 * t383 - t337;
t75 = t313 - t331 / 0.2e1;
t74 = -t173 * t385 + t370;
t73 = -t172 * t385 - t371;
t72 = -t173 * t386 - t338;
t60 = t291 * t221 + t247 * t286 - t472;
t53 = t289 * t83 - t291 * t82;
t52 = t289 * t81 - t291 * (-(-t191 * t301 + t398) * t289 - t187 * t291);
t49 = t289 * t74 - t291 * t73;
t48 = t289 * t72 - t291 * (-(-t174 * t290 + t400) * t289 - t170 * t291);
t46 = -t444 + t453;
t39 = t333 + t445;
t37 = t447 / 0.2e1;
t35 = t448 / 0.2e1;
t27 = t423 + t424;
t23 = t24 * qJD(1);
t22 = (t81 - t153 + (t188 + t398) * t291 + t367) * t291 + t366 * t289;
t21 = (t334 * t291 - t366 + t83) * t291 + (t334 * t289 + t337 + t82) * t289;
t20 = (t72 - t135 + (t171 + t400) * t291 + t371) * t291 + t370 * t289;
t19 = (t335 * t291 - t370 + t74) * t291 + (t335 * t289 + t338 + t73) * t289;
t18 = (t252 / 0.2e1 + t251 / 0.2e1) * t301 + (t232 / 0.2e1 + t231 / 0.2e1) * t290 + t456 + t454 + t446 + t333 - t481;
t8 = m(6) * (t228 * t489 - t472 * t96) + t426;
t7 = t8 * qJD(5);
t4 = t35 - t447 / 0.2e1 + t332;
t3 = t37 - t448 / 0.2e1 + t332;
t2 = t35 + t37 - t451 / 0.2e1 + t310;
t1 = t450 + (-t20 / 0.2e1 - t22 / 0.2e1 + t53 / 0.2e1 + t49 / 0.2e1) * t291 + (t48 / 0.2e1 + t52 / 0.2e1 + t21 / 0.2e1 + t19 / 0.2e1) * t289 + t345;
t5 = [t18 * qJD(3) + t46 * qJD(4) + t39 * qJD(5), -t312, t18 * qJD(1) - t395 + t27 * qJD(4) + t2 * qJD(5) + ((-t112 * t128 + t130 * t468 + (-t129 - t132) * t475) * t457 + (t116 * t180 + t117 * t178) * t458 + ((-t133 * t291 - t134 * t289) * t256 + (-t216 * t291 + t217 * t289) * t254) * t459) * t460 + ((t327 + t325) * (t286 / 0.2e1 + t287 / 0.2e1) + t310 - t450 + (t288 * t360 + t290 * t358 + t299 * t356 + t301 * t354) * t439 + (t20 + t22) * t485 + (t48 + t52 + t21 + t19) * t440 + (-t288 * t361 + t290 * t359 - t299 * t357 + t301 * t355 + t49 + t53) * t437) * qJD(3), qJD(1) * t46 + qJD(3) * t27 + qJD(5) * t75, t39 * qJD(1) + t2 * qJD(3) + t75 * qJD(4) + (t310 + (t311 + (-t193 * t291 + t194 * t289) * t227) * m(6)) * qJD(5); t312, 0, t23 - t102 + (t119 * t459 + t60 * t457 - t94 * t458) * t460, 0, -qJD(3) * t432 - t102; (t308 - (t232 + t231) * t290 / 0.2e1 - (t252 + t251) * t301 / 0.2e1 + t481) * qJD(1) + t395 + t1 * qJD(3) - t25 * qJD(4) + t4 * qJD(5) + (-t456 / 0.4e1 - t446 / 0.4e1 - t454 / 0.4e1) * t461, t23, t1 * qJD(1) + (m(6) * (-t128 * t475 - t129 * t130 + t47 * t60) + m(5) * (-t477 * t178 - t476 * t180 - (t289 * t318 + t291 * (rSges(5,1) * t381 + t344) + t368) * t94) + m(4) * (t254 * t256 * t349 + (t289 * (rSges(4,1) * t382 - t350) + t291 * (rSges(4,1) * t377 + t289 * rSges(4,3) - t262)) * t119) + t426 + ((t463 * t291 + (t306 - t471) * t289 + t464) * t291 + t470 * t286) * t439 + ((t306 * t289 + (t463 - t470) * t291 + t464) * t289 + t471 * t287) * t437) * qJD(3) + t478, -t490, t4 * qJD(1) + t6 * qJD(3) + t478; t25 * qJD(3) + t76 * qJD(5) + (t444 / 0.4e1 - t453 / 0.4e1) * t461, 0, t490 + ((-t128 * t291 + t130 * t289) * t457 + (-t178 * t291 + t180 * t289) * t458) * t460, 0, t376; (t308 - t445) * qJD(1) + t3 * qJD(3) - t76 * qJD(4) + t345 * qJD(5), 0, t3 * qJD(1) + ((t60 * t96 + (-t128 * t289 - t130 * t291) * t227) * m(6) + t426) * qJD(3) + t7, -t376, qJD(1) * t345 + qJD(3) * t8 + t7;];
Cq = t5;
