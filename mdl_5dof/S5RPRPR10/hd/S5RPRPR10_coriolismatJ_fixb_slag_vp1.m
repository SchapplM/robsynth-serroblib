% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR10_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR10_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:55
% EndTime: 2019-12-31 18:26:01
% DurationCPUTime: 4.37s
% Computational Cost: add. (15780->296), mult. (17329->392), div. (0->0), fcn. (19862->8), ass. (0->213)
t337 = m(6) * qJD(5);
t294 = qJ(3) + pkin(8);
t276 = sin(t294);
t277 = cos(t294);
t352 = sin(qJ(1));
t354 = cos(qJ(1));
t195 = t354 * t276 - t352 * t277;
t194 = -t352 * t276 - t354 * t277;
t228 = sin(qJ(5));
t230 = cos(qJ(5));
t315 = t230 * Icges(6,4);
t260 = -Icges(6,2) * t228 + t315;
t139 = Icges(6,6) * t194 + t260 * t195;
t332 = Icges(6,4) * t228;
t262 = Icges(6,1) * t230 - t332;
t143 = Icges(6,5) * t194 + t262 * t195;
t421 = t139 * t230 + t228 * t143;
t334 = t421 * t195;
t138 = -Icges(6,6) * t195 + t260 * t194;
t142 = -Icges(6,5) * t195 + t262 * t194;
t422 = t138 * t230 + t228 * t142;
t432 = t422 * t194;
t433 = t334 / 0.2e1 - t432 / 0.2e1;
t431 = -t195 / 0.4e1;
t213 = -t228 * rSges(6,1) - rSges(6,2) * t230;
t165 = t195 * t213;
t166 = t213 * t194;
t427 = m(6) * (t352 * t165 + t166 * t354);
t286 = -t427 / 0.2e1;
t104 = t427 / 0.2e1;
t424 = t104 + t286;
t429 = qJD(5) * t424;
t310 = t424 * qJD(2);
t415 = -rSges(6,1) * t230 + t228 * rSges(6,2);
t428 = t195 * rSges(6,3) + t415 * t194;
t146 = -t194 * rSges(6,3) + t415 * t195;
t381 = -m(6) / 0.2e1;
t317 = t138 * t228;
t318 = t139 * t228;
t313 = t142 * t230;
t314 = t143 * t230;
t258 = Icges(6,5) * t230 - Icges(6,6) * t228;
t136 = -Icges(6,3) * t194 - t258 * t195;
t426 = t195 * t136;
t137 = Icges(6,3) * t195 - t258 * t194;
t425 = t195 * t137;
t309 = t194 * t136 - t195 * t314;
t401 = -t137 - t318;
t423 = -t401 * t195 + t309;
t383 = m(5) / 0.2e1;
t380 = m(6) / 0.2e1;
t409 = -t228 / 0.2e1;
t420 = -t230 / 0.2e1;
t419 = t230 / 0.2e1;
t353 = cos(qJ(3));
t270 = t352 * t353;
t229 = sin(qJ(3));
t281 = t354 * t229;
t206 = -t270 + t281;
t280 = t352 * t229;
t243 = t354 * t353 + t280;
t177 = -rSges(4,1) * t243 + t206 * rSges(4,2);
t265 = -t206 * rSges(4,1) - rSges(4,2) * t243;
t397 = t177 * t352 + t354 * t265;
t349 = m(4) * t397;
t418 = t349 / 0.2e1;
t257 = Icges(6,5) * t228 + Icges(6,6) * t230;
t160 = t257 * t194;
t417 = t257 * t195;
t226 = t353 * pkin(3) + pkin(2);
t217 = t352 * t226;
t273 = pkin(3) * t281;
t255 = t273 - t217;
t287 = t352 * pkin(2);
t241 = -t287 - t255;
t400 = t195 * rSges(5,1) - t194 * rSges(5,2);
t152 = t241 - t400;
t147 = t152 * t354;
t256 = -pkin(3) * t280 - t354 * t226;
t288 = t354 * pkin(2);
t240 = t288 + t256;
t264 = -t194 * rSges(5,1) - t195 * rSges(5,2);
t153 = t240 - t264;
t232 = t195 * pkin(4) + t194 * pkin(7) - t146 + t273;
t231 = t217 - t232;
t89 = -t287 + t231;
t86 = t89 * t354;
t253 = pkin(4) - t415;
t355 = rSges(6,3) + pkin(7);
t392 = -t253 * t194 + t355 * t195;
t90 = t240 - t392;
t289 = t418 + (t90 * t352 + t86) * t380 + (t153 * t352 + t147) * t383;
t236 = -t255 - t400;
t151 = -t287 + t236;
t393 = t355 * t194 + t253 * t195;
t88 = t241 - t393;
t293 = (-t88 * t354 + t86) * t380 + (-t151 * t354 + t147) * t383;
t416 = t289 - t293;
t413 = 0.2e1 * t286;
t238 = t243 * pkin(3);
t154 = t238 + t264;
t222 = pkin(3) * t270;
t155 = t273 - t222 + t400;
t93 = -t194 * pkin(4) + t195 * pkin(7) + t238 + t428;
t94 = -t222 + t232;
t396 = (t352 * t93 + t94 * t354) * t381 - m(5) * (t352 * t154 + t155 * t354) / 0.2e1;
t247 = t354 * pkin(1) + t352 * qJ(2);
t170 = t288 + t247 - t177;
t267 = -t352 * pkin(1) + t354 * qJ(2);
t233 = -t287 - t265 + t267;
t412 = t170 * t265 - t233 * t177;
t329 = Icges(6,2) * t230;
t259 = t329 + t332;
t248 = t228 * t259;
t261 = Icges(6,1) * t228 + t315;
t311 = t230 * t261;
t399 = t248 - t311;
t108 = t399 * t194 + t417;
t411 = -t194 / 0.2e1;
t410 = t195 / 0.2e1;
t408 = (-t262 + t329) * t419 + (t260 + t261 + t315) * t228 / 0.2e1;
t339 = t88 - t89;
t407 = m(6) * t339;
t81 = -t231 + t267;
t237 = t255 + t267;
t83 = t237 + t393;
t341 = t81 - t83;
t406 = m(6) * t341;
t404 = t195 * t194;
t402 = -t136 + t317;
t109 = t399 * t195 - t160;
t320 = t195 * t109;
t321 = t194 * t108;
t398 = t320 / 0.4e1 + t321 / 0.4e1;
t395 = t432 / 0.4e1 + t165 * t407 / 0.2e1;
t234 = t247 - t256;
t82 = t234 + t392;
t394 = t334 / 0.4e1 + ((-t81 + t94) * t195 + (t82 - t93) * t194) * t213 * t381;
t271 = t248 / 0.2e1 + t262 * t409 + t260 * t420 - t311 / 0.2e1;
t391 = t194 ^ 2;
t390 = t195 ^ 2;
t389 = 2 * qJD(1);
t388 = 4 * qJD(1);
t387 = 2 * qJD(3);
t386 = 4 * qJD(3);
t385 = m(4) / 0.2e1;
t382 = m(5) / 0.4e1;
t38 = -t195 * t318 - t309;
t308 = -t194 * t137 + t195 * t313;
t39 = -t195 * t317 + t308;
t22 = -t194 * t38 + t195 * t39;
t379 = t22 / 0.2e1;
t307 = t194 * t314 + t426;
t40 = -t194 * t318 + t307;
t306 = t194 * t313 + t425;
t41 = -t194 * t317 + t306;
t23 = -t194 * t40 + t195 * t41;
t378 = -t23 / 0.2e1;
t377 = m(4) * t412;
t300 = t151 - t152;
t373 = m(5) * t300 * t153;
t132 = t237 + t400;
t133 = t234 + t264;
t372 = m(5) * (t132 * t154 - t133 * t155);
t124 = t132 * t354;
t371 = m(5) * (t133 * t352 + t124);
t340 = t82 - t90;
t342 = -t81 + t89;
t368 = m(6) * (t342 * t165 + t340 * t166);
t364 = t90 * t407;
t363 = m(6) * (t81 * t93 - t82 * t94);
t362 = m(6) * (t165 * t81 - t166 * t82);
t361 = m(6) * (-t165 * t89 + t166 * t90);
t79 = t81 * t354;
t360 = m(6) * (t82 * t352 + t79);
t351 = m(3) * ((t354 * rSges(3,3) + t267) * t354 + (t352 * rSges(3,3) + t247) * t352);
t350 = m(4) * (t170 * t352 + t233 * t354);
t324 = t165 * t213;
t323 = t166 * t213;
t131 = -t236 + t267;
t305 = t131 - t132;
t296 = qJD(5) * t194;
t295 = qJD(5) * t195;
t290 = -t349 / 0.2e1 - t396;
t67 = t432 / 0.2e1;
t269 = t379 + (-t261 * t194 - t138) * t409 + (t259 * t194 - t142) * t419 + t194 * t408 + t258 * t410;
t268 = t378 + (-t261 * t195 - t139) * t409 + (-t259 * t195 + t143) * t420 + t258 * t411 + t195 * t408;
t252 = t271 + t433;
t251 = -t334 / 0.2e1 + t67 + t271;
t250 = -t321 / 0.4e1 - t395 + t398 + (-t421 + t109) * t431;
t249 = -t320 / 0.4e1 - t394 + t398 - (t108 - t422) * t194 / 0.4e1;
t85 = -t165 * t195 - t166 * t194;
t49 = 0.2e1 * t104;
t32 = -t271 + t361;
t31 = -t271 + t362;
t25 = (-t138 * t195 - t139 * t194) * t228 + t307 + t308;
t21 = t350 + t351 + t360 + t371;
t15 = t368 / 0.2e1;
t14 = t364 + t373;
t13 = t363 + t372 + t377;
t12 = (t38 + (-t314 + t318) * t195 + t306) * t195 + (-t195 * t402 - t25 + t39) * t194;
t11 = (t25 - t40) * t195 + (-t41 + (t313 - t317) * t194 + t423) * t194;
t10 = (t40 + (t402 - t313) * t195) * t195 + (t402 * t194 - t306 + t41 + t425) * t194;
t9 = (-t38 - t423) * t195 + (-t39 + (t401 + t314) * t194 + t426) * t194;
t8 = t418 + t289 + t293 + t396;
t7 = t290 - t416;
t6 = t290 + t416;
t5 = -t368 / 0.2e1 + t249 + t250 - t271;
t4 = t250 - t432 / 0.4e1 + t15 + t251 + t394;
t3 = t421 * t431 + t15 + t249 + t252 + t395;
t2 = (t11 / 0.2e1 - t22 / 0.2e1) * t195 + (t378 - t9 / 0.2e1) * t194;
t1 = (t10 / 0.2e1 + t379) * t195 + (t23 / 0.2e1 - t12 / 0.2e1) * t194;
t16 = [(-t82 * t406 / 0.4e1 + t305 * t133 * t382) * t388 + t21 * qJD(2) + t13 * qJD(3) + t31 * qJD(5), qJD(1) * t21 + qJD(3) * t6 + t429, t13 * qJD(1) + t6 * qJD(2) + t3 * qJD(5) + (-t364 / 0.4e1 - t373 / 0.4e1) * t386 + ((t340 * t94 + t342 * t93) * t380 + ((t133 - t153) * t155 + (-t132 + t152) * t154) * t383 - t412 * t385) * t387, 0, t31 * qJD(1) + t310 + t3 * qJD(3) + (-t10 / 0.2e1 - t269) * t295 + (t12 / 0.2e1 + t268) * t296 + ((-t415 * t82 + t323) * t195 + (-t415 * t81 - t324) * t194) * t337; t7 * qJD(3) + t49 * qJD(5) + (-t360 / 0.4e1 - t371 / 0.4e1 - t350 / 0.4e1 - t351 / 0.4e1) * t388 + ((-t354 * t83 + t79) * t380 + (-t354 * t131 + t124) * t383) * t389, 0, t7 * qJD(1) + t413 * qJD(5) + (t397 * t385 + t396) * t387, 0, t49 * qJD(1) + t413 * qJD(3) + (-t352 * t194 + t354 * t195) * t415 * t337; t8 * qJD(2) + t14 * qJD(3) + t4 * qJD(5) + (-t363 / 0.4e1 - t372 / 0.4e1 - t377 / 0.4e1) * t388 + ((t339 * t82 - t341 * t90) * t380 + (t300 * t133 + t305 * t153) * t383) * t389, qJD(1) * t8 - t429, t14 * qJD(1) + t32 * qJD(5) + (m(6) * (-t89 * t93 + t90 * t94) / 0.4e1 + (-t152 * t154 + t153 * t155) * t382) * t386, 0, t4 * qJD(1) - t310 + t32 * qJD(3) + (-t11 / 0.2e1 + t269) * t295 + (t9 / 0.2e1 - t268) * t296 + ((-t415 * t90 - t323) * t195 + (-t415 * t89 + t324) * t194) * t337; 0, 0, 0, 0, -t85 * t337; t413 * qJD(2) + t5 * qJD(3) + t1 * qJD(5) + (t252 + t67 - t362 + (-t421 / 0.2e1 + t213 * t406) * t195) * qJD(1), qJD(1) * t413 + qJD(3) * t424, t5 * qJD(1) + t2 * qJD(5) + t310 + (((-t89 - t94) * t195 + (t90 + t93) * t194) * t213 * m(6) + t251 - t361 + t433) * qJD(3), 0, t1 * qJD(1) + t2 * qJD(3) + (m(6) * ((t195 * t146 + t194 * t428) * t85 + (t390 + t391) * t415 * t213) + (t390 * t160 - t404 * t417) * t410 + (-t160 * t404 + t391 * t417) * t411) * qJD(5);];
Cq = t16;
