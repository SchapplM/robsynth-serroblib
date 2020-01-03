% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP7_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:16
% DurationCPUTime: 10.93s
% Computational Cost: add. (2850->385), mult. (7410->500), div. (0->0), fcn. (5690->4), ass. (0->218)
t467 = -Icges(5,4) - Icges(4,5);
t466 = Icges(4,6) - Icges(5,6);
t214 = sin(qJ(1));
t216 = cos(qJ(1));
t215 = cos(qJ(3));
t213 = sin(qJ(3));
t342 = Icges(4,4) * t213;
t255 = Icges(4,2) * t215 + t342;
t465 = t466 * t214 - t216 * t255;
t341 = Icges(4,4) * t215;
t257 = Icges(4,1) * t213 + t341;
t464 = t467 * t214 + t216 * t257;
t325 = t215 * t216;
t327 = t213 * t216;
t459 = Icges(5,5) * t327 - Icges(5,3) * t325 + t465;
t182 = Icges(5,5) * t325;
t454 = -Icges(5,1) * t327 + t182 - t464;
t253 = Icges(4,5) * t213 + Icges(4,6) * t215;
t93 = Icges(4,3) * t216 + t214 * t253;
t254 = Icges(5,4) * t213 - Icges(5,6) * t215;
t95 = Icges(5,2) * t216 + t214 * t254;
t463 = t93 + t95;
t462 = t466 * t213 + t467 * t215;
t338 = Icges(5,5) * t213;
t155 = Icges(5,1) * t215 + t338;
t157 = Icges(4,1) * t215 - t342;
t449 = t155 + t157;
t251 = -Icges(5,3) * t215 + t338;
t461 = t251 - t255;
t203 = Icges(5,5) * t215;
t256 = Icges(5,1) * t213 - t203;
t460 = -t256 - t257;
t433 = t454 * t213 + t459 * t215;
t94 = -Icges(4,3) * t214 + t216 * t253;
t337 = Icges(5,2) * t214;
t96 = Icges(5,4) * t327 - Icges(5,6) * t325 - t337;
t458 = -t94 - t96;
t153 = -Icges(4,2) * t213 + t341;
t248 = t153 * t215 + t157 * t213;
t385 = Icges(5,3) * t213 + t203;
t448 = t155 * t213 - t215 * t385 + t248;
t328 = t213 * t214;
t181 = Icges(5,5) * t328;
t326 = t214 * t215;
t335 = Icges(5,6) * t216;
t91 = -Icges(5,3) * t326 + t181 + t335;
t97 = Icges(4,6) * t216 + t214 * t255;
t457 = -t91 + t97;
t183 = Icges(4,4) * t326;
t339 = Icges(4,5) * t216;
t101 = Icges(4,1) * t328 + t183 + t339;
t99 = Icges(5,4) * t216 + t214 * t256;
t455 = -t101 - t99;
t411 = -t463 * t214 - t91 * t325;
t409 = t462 * t214;
t453 = t449 * t216;
t452 = t461 * qJD(3);
t451 = t460 * qJD(3);
t450 = -t153 + t385;
t447 = -t253 - t254;
t402 = t462 * t216;
t446 = t433 * t216;
t362 = t216 * t95 + t99 * t328;
t31 = -t326 * t91 + t362;
t33 = t101 * t328 + t216 * t93 + t97 * t326;
t445 = t31 + t33;
t415 = t458 * t216 + t459 * t326 + t454 * t328;
t260 = t101 * t213 + t215 * t97;
t444 = -t216 * t260 - t327 * t99 - t411;
t413 = t458 * t214 - t446;
t443 = t448 * t214 - t402;
t442 = t155 * t327 + t216 * t248 - t325 * t385 + t409;
t406 = t213 * t99 + t260;
t412 = t457 * t213 + t455 * t215;
t118 = t153 * t216;
t297 = qJD(3) * t216;
t441 = qJD(3) * t118 - t385 * t297 + (t214 * t251 + t335 - t97) * qJD(1);
t111 = t385 * t214;
t298 = qJD(3) * t214;
t440 = qJD(3) * t111 - t153 * t298 + (t216 * t251 + t465) * qJD(1);
t439 = -t453 * qJD(3) + (t214 * t257 + t339 + t99) * qJD(1);
t121 = t157 * t214;
t438 = qJD(3) * t121 + t155 * t298 + (t216 * t256 + t464) * qJD(1);
t437 = t448 * qJD(1) + t447 * qJD(3);
t436 = (Icges(4,2) * t328 + t111 - t183 + t455) * t216 + (-Icges(5,3) * t327 + t118 - t182 - t454) * t214;
t435 = (Icges(5,1) * t326 + t121 + t181 - t457) * t216 + (-t453 - t459) * t214;
t353 = t215 * t91;
t434 = t353 - t406;
t405 = -t213 * t459 + t454 * t215;
t198 = t216 * qJ(2);
t161 = pkin(1) * t214 - t198;
t145 = qJD(1) * t161;
t195 = qJD(2) * t214;
t299 = qJD(1) * t216;
t305 = qJ(2) * t299 + t195;
t432 = t305 - t195 + t145;
t431 = -t450 - t460;
t430 = -t449 - t461;
t429 = t452 * t215 + t451 * t213 + (t450 * t213 + t449 * t215) * qJD(3) + t462 * qJD(1);
t369 = rSges(5,1) + pkin(3);
t428 = t443 * qJD(1);
t349 = rSges(5,3) + qJ(4);
t427 = (t413 * t214 + t444 * t216) * qJD(3);
t426 = (t415 * t214 + t445 * t216) * qJD(3);
t425 = t442 * qJD(1);
t424 = t463 * qJD(1);
t423 = (t213 * t431 + t215 * t430) * qJD(1);
t371 = t214 / 0.2e1;
t422 = -t216 / 0.2e1;
t421 = t426 + t428;
t420 = -t425 + t427;
t419 = qJD(3) * t434 + t213 * t440 + t215 * t438;
t418 = -qJD(3) * t433 + t213 * t441 + t215 * t439;
t417 = t214 * t437 - t216 * t429;
t416 = t214 * t429 + t216 * t437;
t209 = t216 * rSges(4,3);
t104 = rSges(4,1) * t328 + rSges(4,2) * t326 + t209;
t167 = t216 * pkin(1) + t214 * qJ(2);
t211 = t216 * pkin(5);
t384 = t211 + t167;
t410 = t104 + t384;
t265 = rSges(5,1) * t213 - rSges(5,3) * t215;
t309 = -pkin(3) * t213 + qJ(4) * t215 - t265;
t164 = pkin(3) * t215 + qJ(4) * t213;
t165 = rSges(5,1) * t215 + rSges(5,3) * t213;
t308 = t164 + t165;
t332 = qJD(1) * t94;
t408 = t332 - t409 * qJD(3) + (t216 * t254 - t337 - t434) * qJD(1);
t407 = qJD(1) * t433 + t402 * qJD(3) + t424;
t404 = t213 * t299 + t215 * t298;
t403 = t447 * qJD(1);
t400 = qJD(1) * t96 + qJD(3) * t405 + t213 * t439 - t215 * t441 + t332;
t399 = qJD(3) * t412 - t438 * t213 + t440 * t215 + t424;
t398 = -t213 * t435 + t215 * t436;
t397 = 0.2e1 * qJD(3);
t194 = qJD(4) * t213;
t306 = -t214 * rSges(5,2) - rSges(5,3) * t325;
t387 = -qJ(4) * t325 + t306;
t318 = t369 * t327 + t387;
t210 = t216 * rSges(5,2);
t386 = t369 * t328 + t210;
t319 = t349 * t326 - t386;
t30 = t194 + (t214 * t319 - t216 * t318) * qJD(3);
t396 = qJD(3) * t30;
t296 = qJD(4) * t215;
t239 = qJD(3) * t308 - t296;
t395 = t216 * t239;
t274 = -rSges(3,2) * t216 + t214 * rSges(3,3);
t388 = t167 + t274;
t286 = t213 * t298;
t383 = t349 * t286 + t369 * t404;
t217 = qJD(1) ^ 2;
t295 = qJD(1) * qJD(2);
t300 = qJD(1) * t214;
t317 = qJD(1) * (-pkin(1) * t300 + t305) + t214 * t295;
t367 = pkin(5) * t214;
t246 = -t217 * t367 + t317;
t316 = t309 * qJD(3) + t194;
t268 = t194 + t316;
t363 = -(-qJ(4) * t299 - qJD(4) * t214) * t215 - t306 * qJD(1) - t383;
t14 = -t268 * t297 + (t214 * t239 - t363) * qJD(1) + t246;
t196 = qJD(2) * t216;
t109 = qJD(1) * t167 - t196;
t190 = t216 * t295;
t273 = -t211 * t217 + t190;
t129 = t165 * t216;
t364 = (pkin(3) * t300 - qJ(4) * t297) * t213 + (-qJ(4) * t300 + (-pkin(3) * qJD(3) + qJD(4)) * t216) * t215 - qJD(3) * t129 + (t214 * t265 + t210) * qJD(1);
t15 = t268 * t298 + (-t109 - t364 + t395) * qJD(1) + t273;
t382 = t14 * t422 + t15 * t371;
t373 = t214 ^ 2;
t372 = -pkin(1) - pkin(5);
t368 = rSges(3,2) - pkin(1);
t366 = -qJD(1) / 0.2e1;
t359 = rSges(4,2) * t213;
t358 = rSges(4,2) * t215;
t357 = rSges(3,3) * t216;
t166 = rSges(4,1) * t215 - t359;
t133 = t166 * t298;
t266 = rSges(4,1) * t213 + t358;
t143 = t266 * qJD(3);
t290 = t404 * rSges(4,1) + t299 * t358;
t73 = (-rSges(4,3) * qJD(1) - qJD(3) * t359) * t214 + t290;
t24 = t143 * t297 + (t73 + t133) * qJD(1) + t246;
t352 = t216 * t24;
t206 = t214 * rSges(4,3);
t106 = t216 * t266 - t206;
t278 = t106 - t367;
t41 = t133 + t195 + (-t161 + t278) * qJD(1);
t351 = t216 * t41;
t287 = t166 * t297;
t130 = t166 * t216;
t71 = -qJD(3) * t130 + (t214 * t266 + t209) * qJD(1);
t25 = -t143 * t298 + (-t109 - t71 + t287) * qJD(1) + t273;
t350 = t25 * t214;
t315 = t308 * t214;
t314 = -t164 * t216 - t129;
t304 = rSges(3,2) * t300 + rSges(3,3) * t299;
t294 = -rSges(5,2) + t372;
t293 = -rSges(4,3) + t372;
t282 = -t298 / 0.2e1;
t280 = -t297 / 0.2e1;
t279 = t297 / 0.2e1;
t276 = t96 - t353;
t272 = t318 - t367;
t42 = t410 * qJD(1) - t196 - t287;
t263 = t214 * t41 - t216 * t42;
t240 = -t214 * t296 + t308 * t298 + t195;
t43 = (-t104 * t214 - t106 * t216) * qJD(3);
t162 = rSges(3,2) * t214 + t357;
t126 = t166 * t214;
t77 = qJD(1) * t388 - t196;
t76 = t195 + (-t161 + t162) * qJD(1);
t69 = t190 + (-qJD(1) * t274 - t109) * qJD(1);
t68 = qJD(1) * t304 + t317;
t29 = -t196 - t395 + (t384 - t319) * qJD(1);
t28 = (-t161 + t272) * qJD(1) + t240;
t1 = (t296 + t364 * t216 + t363 * t214 + (t214 * t318 + t216 * t319) * qJD(1)) * qJD(3);
t2 = [(((t362 + t33 + t413 + t446) * t216 + ((t276 + t94 - t406) * t216 - t411 - t444 + t415) * t214) * qJD(3) + t428) * t282 + (-(-t145 + t240 - t28) * t29 + t15 * (t198 + t387) + t28 * t196 + t14 * (t384 + t386) + t29 * (t305 + t383) + (t15 * t372 + (-t29 * qJD(4) - t14 * t349) * t215) * t214 + (t15 * t369 * t213 + (-t296 + (t213 * t349 + t215 * t369) * qJD(3)) * t28) * t216) * m(5) + (t25 * (-t161 - t206 - t367) + t41 * t196 + t24 * t410 + (qJD(3) * t166 * t41 + t25 * t266) * t216 + (-rSges(4,2) * t286 - t133 + t290 + t41 + t432) * t42) * m(4) + (t69 * (t214 * t368 + t198 + t357) + t76 * t196 + t68 * t388 + (t76 + t304 + t432) * t77) * m(3) + ((t373 * t94 + (t214 * t276 - t31 + t362) * t214 + ((t406 - t458) * t216 + t411 + t415) * t216) * qJD(3) + t420 + t425) * t280 + (t416 + t419) * t279 + (t417 + t418 + t421) * t298 / 0.2e1 + ((-t412 + t443) * t214 + t442 * t216) * qJD(3) * t366 + (t405 * t279 + t451 * t215 - t452 * t213 - t448 * qJD(3) + (-t272 * t29 + (-t215 * t29 * t349 + t28 * t294) * t216 + (t28 * (-qJ(2) + t309) + t29 * t294) * t214) * m(5) + (-t278 * t42 + t293 * t351 + (t41 * (-qJ(2) - t266) + t42 * t293) * t214) * m(4) + (-t162 * t77 + t76 * t368 * t216 + (t76 * (-rSges(3,3) - qJ(2)) - t77 * pkin(1)) * t214) * m(3)) * qJD(1); 0.2e1 * t382 * m(5) + 0.2e1 * (t350 / 0.2e1 - t352 / 0.2e1) * m(4) + 0.2e1 * (t69 * t371 + t68 * t422) * m(3); ((t213 * t436 + t215 * t435) * qJD(3) + (t213 * t430 - t215 * t431) * qJD(1)) * t366 + (t419 * t216 + t418 * t214 + (t412 * t214 + t405 * t216) * qJD(1)) * qJD(1) / 0.2e1 + ((t402 * t298 + t403) * t214 + ((-t409 * t214 + t398) * qJD(3) + t423) * t216) * t282 + ((-t409 * t297 + t403) * t216 + ((t402 * t216 - t398) * qJD(3) - t423) * t214) * t280 + ((-t14 * t308 - t29 * t316 - t1 * t318 + t30 * t364 + (t28 * t308 + t30 * t319) * qJD(1)) * t216 + (t15 * t308 + t28 * t316 + t1 * t319 + t30 * t363 + (t29 * t308 + t30 * t318) * qJD(1)) * t214 - (t30 * t215 + (t214 * t28 - t216 * t29) * t213) * qJD(4) - (-t28 * t314 + t29 * t315) * qJD(1) - ((-t29 * t309 + t30 * t314) * t216 + (t28 * t309 - t30 * t315) * t214) * qJD(3)) * m(5) + (0.2e1 * t43 * (-t214 * t73 + t216 * t71 + (-t104 * t216 + t106 * t214) * qJD(1)) - t263 * t143 + (t350 - t352 + (t214 * t42 + t351) * qJD(1)) * t166 - (t126 * t42 + t130 * t41) * qJD(1) - (t43 * (-t126 * t214 - t130 * t216) - t263 * t266) * qJD(3)) * m(4) + (t417 * qJD(1) + ((t413 * qJD(1) + t399 * t216) * t216 + (t407 * t214 - t444 * qJD(1) + (-t400 + t408) * t216) * t214) * t397) * t371 + (t416 * qJD(1) + ((t415 * qJD(1) + t408 * t216) * t216 + (t400 * t214 - t445 * qJD(1) + (-t399 + t407) * t216) * t214) * t397) * t216 / 0.2e1 - (t421 + t426) * t300 / 0.2e1 + (t420 + t427) * t299 / 0.2e1; (t1 * t213 + 0.2e1 * (t396 / 0.2e1 - (t216 ^ 2 + t373) * t396 / 0.2e1 - t382) * t215) * m(5);];
tauc = t2(:);
