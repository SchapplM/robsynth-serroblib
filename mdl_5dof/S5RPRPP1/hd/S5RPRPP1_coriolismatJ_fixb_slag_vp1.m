% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:39
% EndTime: 2019-12-31 18:08:49
% DurationCPUTime: 5.87s
% Computational Cost: add. (11256->320), mult. (9112->430), div. (0->0), fcn. (8179->8), ass. (0->205)
t419 = Icges(4,3) + Icges(5,3);
t248 = qJ(1) + pkin(7);
t243 = sin(t248);
t249 = sin(qJ(3));
t353 = pkin(3) * t249;
t247 = qJ(3) + pkin(8);
t244 = cos(t247);
t242 = sin(t247);
t356 = rSges(6,1) + pkin(4);
t294 = t356 * t242;
t340 = rSges(6,3) + qJ(5);
t402 = t340 * t244 - t294;
t275 = t353 - t402;
t416 = t275 * t243;
t418 = t243 * t416;
t245 = cos(t248);
t387 = t340 * t242 + t356 * t244;
t417 = -t245 * rSges(6,2) + t387 * t243;
t184 = Icges(5,5) * t244 - Icges(5,6) * t242;
t251 = cos(qJ(3));
t207 = Icges(4,5) * t251 - Icges(4,6) * t249;
t415 = (t184 + t207) * t245 + t419 * t243;
t230 = Icges(6,5) * t242;
t337 = Icges(6,1) * t244;
t269 = t230 + t337;
t138 = Icges(6,4) * t243 + t269 * t245;
t335 = Icges(5,4) * t242;
t191 = Icges(5,1) * t244 - t335;
t140 = Icges(5,5) * t243 + t191 * t245;
t414 = t138 + t140;
t318 = t243 * t251;
t319 = t243 * t249;
t320 = t243 * t244;
t323 = t242 * t243;
t410 = -Icges(4,5) * t318 - Icges(5,5) * t320 + Icges(4,6) * t319 + Icges(5,6) * t323 + t419 * t245;
t376 = m(6) / 0.2e1;
t194 = rSges(5,1) * t242 + rSges(5,2) * t244;
t257 = t194 + t353;
t392 = t257 * t245;
t393 = t257 * t243;
t397 = t243 * t393 + t245 * t392;
t96 = t275 * t245;
t348 = (-t245 * t96 - t418) * t376 - m(5) * t397 / 0.2e1;
t377 = m(5) / 0.2e1;
t317 = t244 * t245;
t300 = t340 * t317;
t86 = (-t294 - t353) * t245 + t300;
t350 = (-t245 * t86 + t418) * t376 + t397 * t377;
t13 = t350 - t348;
t413 = t13 * qJD(1);
t336 = Icges(4,4) * t249;
t211 = Icges(4,1) * t251 - t336;
t154 = Icges(4,5) * t243 + t211 * t245;
t412 = -t140 * t320 - t154 * t318;
t409 = t415 * t245 + t412;
t186 = Icges(5,2) * t244 + t335;
t331 = Icges(6,3) * t244;
t265 = t331 - t230;
t408 = (-t186 - t265) * t245 + t414;
t202 = Icges(5,4) * t323;
t139 = Icges(5,1) * t320 - Icges(5,5) * t245 - t202;
t218 = Icges(4,4) * t319;
t153 = Icges(4,1) * t318 - Icges(4,5) * t245 - t218;
t314 = t245 * t251;
t407 = -t139 * t317 - t153 * t314 + t410 * t243;
t135 = Icges(5,4) * t320 - Icges(5,2) * t323 - Icges(5,6) * t245;
t151 = Icges(4,4) * t318 - Icges(4,2) * t319 - Icges(4,6) * t245;
t406 = t135 * t242 + t151 * t249;
t201 = Icges(6,5) * t317;
t321 = t242 * t245;
t130 = Icges(6,6) * t243 + Icges(6,3) * t321 + t201;
t185 = Icges(6,4) * t244 + Icges(6,6) * t242;
t134 = Icges(6,2) * t243 + t185 * t245;
t405 = t130 * t321 + t154 * t314 + t414 * t317 + (t134 + t415) * t243;
t404 = -Icges(4,5) * t249 - Icges(4,6) * t251 + (-Icges(5,6) + Icges(6,6)) * t244 + (-Icges(6,4) - Icges(5,5)) * t242;
t240 = t243 ^ 2;
t241 = t245 ^ 2;
t297 = t240 + t241;
t234 = Icges(5,4) * t244;
t332 = Icges(5,2) * t242;
t136 = Icges(5,6) * t243 + (t234 - t332) * t245;
t246 = Icges(4,4) * t251;
t209 = -Icges(4,2) * t249 + t246;
t152 = Icges(4,6) * t243 + t209 * t245;
t403 = t136 * t242 + t152 * t249 + t410;
t315 = t245 * t249;
t401 = -t136 * t321 - t152 * t315 + t405;
t329 = (-Icges(6,2) * t245 + t185 * t243) * t245;
t400 = t329 + t405;
t399 = -t135 * t321 - t136 * t323 - t151 * t315 - t152 * t319 - t407 - t409;
t334 = Icges(6,5) * t244;
t183 = Icges(6,3) * t242 + t334;
t188 = Icges(6,1) * t242 - t334;
t208 = Icges(4,2) * t251 + t336;
t338 = Icges(5,1) * t242;
t398 = -(t211 / 0.2e1 - t208 / 0.2e1) * t249 - (t234 + t338 / 0.2e1 - t332 / 0.2e1 + t188 / 0.2e1 - t183 / 0.2e1) * t244 - (t191 / 0.2e1 - t186 / 0.2e1 + t230 + t337 / 0.2e1 - t331 / 0.2e1) * t242;
t396 = -t243 / 0.2e1;
t359 = t243 / 0.2e1;
t395 = -t245 / 0.2e1;
t129 = -Icges(6,6) * t245 + t183 * t243;
t137 = -Icges(6,4) * t245 + t269 * t243;
t389 = (t129 * t242 + t137 * t244) * t243;
t210 = Icges(4,1) * t249 + t246;
t386 = t408 * t243;
t270 = -t234 - t338;
t286 = (t270 * t245 - t136) * t243;
t288 = (-Icges(6,1) * t321 + t130 + t201) * t243;
t385 = t286 + t288;
t384 = t245 * qJD(3);
t383 = t404 * t243;
t382 = t404 * t245;
t379 = 0.4e1 * qJD(1);
t378 = 0.2e1 * qJD(3);
t212 = rSges(4,1) * t249 + rSges(4,2) * t251;
t176 = t212 * t243;
t177 = t212 * t245;
t238 = t245 * pkin(6);
t343 = rSges(4,1) * t251;
t293 = pkin(2) + t343;
t298 = rSges(4,2) * t319 + t245 * rSges(4,3);
t355 = sin(qJ(1)) * pkin(1);
t98 = -t293 * t243 + t238 + t298 - t355;
t220 = rSges(4,2) * t315;
t354 = cos(qJ(1)) * pkin(1);
t99 = t354 - t220 + t293 * t245 + (rSges(4,3) + pkin(6)) * t243;
t375 = m(4) * (t176 * t98 - t177 * t99);
t261 = rSges(5,1) * t320 - rSges(5,2) * t323 - t245 * rSges(5,3);
t352 = pkin(3) * t251;
t239 = pkin(2) + t352;
t351 = -qJ(4) - pkin(6);
t299 = -t243 * t239 - t245 * t351;
t276 = t299 - t355;
t88 = -t261 + t276;
t290 = -rSges(5,2) * t321 + t243 * rSges(5,3);
t221 = t243 * t351;
t292 = -t221 + t354;
t342 = rSges(5,1) * t244;
t89 = (t239 + t342) * t245 + t290 + t292;
t373 = m(5) * (-t392 * t89 + t393 * t88);
t372 = m(5) * (t89 * t243 + t245 * t88);
t70 = t276 - t417;
t341 = t243 * rSges(6,2);
t71 = t341 + (t239 + t387) * t245 + t292;
t347 = t70 * t317 + t71 * t320;
t366 = m(6) * ((t243 * t86 + t245 * t416) * t242 + t347);
t365 = m(6) * (-t321 * t416 + t96 * t323 + t347);
t364 = m(6) * (t416 * t70 + t71 * t86);
t363 = m(6) * (t71 * t243 + t245 * t70);
t169 = t297 * t242;
t90 = m(6) * t169;
t346 = -t96 * t317 - t320 * t416;
t345 = m(6) * qJD(3);
t344 = m(6) * qJD(5);
t322 = t242 * t244;
t313 = t90 * qJD(1);
t310 = -t243 * (pkin(2) * t243 - t238 + t299) + t245 * (-t243 * pkin(6) - t221 + (-pkin(2) + t239) * t245);
t307 = -t188 * t243 + t129;
t306 = -t270 * t243 + t135;
t305 = -t265 * t243 + t137;
t303 = -Icges(5,2) * t320 + t139 - t202;
t301 = t297 * t322;
t123 = (t169 / 0.2e1 - t242 / 0.2e1) * m(6);
t296 = t123 * qJD(2);
t37 = t71 * t321 - t70 * t323;
t295 = m(6) * t37 * qJD(1);
t291 = rSges(5,2) * t242 - t342 - t352;
t289 = t307 * t245;
t287 = t306 * t245;
t285 = t305 * t245;
t283 = t303 * t245;
t277 = t297 * t353;
t274 = -t352 - t387;
t273 = -t130 * t323 + t134 * t245 - t138 * t320;
t272 = t207 / 0.2e1 + t184 / 0.2e1 + t185 / 0.2e1;
t172 = -Icges(4,2) * t318 - t218;
t174 = t210 * t243;
t256 = (t151 + t174) * t251 + (t153 + t172) * t249;
t173 = t208 * t245;
t175 = t210 * t245;
t255 = (-t152 - t175) * t251 + (-t154 + t173) * t249;
t47 = -t329 + t389;
t254 = (t47 - t389 + t400) * t396 + t401 * t359 + ((t406 + t415) * t245 + t399 + t407 + t412) * t395;
t253 = t273 * t396 + (t47 + t410 * t245 + (t139 * t244 + t153 * t251 - t406) * t243) * t395 + (t403 * t245 - t400 + t401) * t245 / 0.2e1 + (t129 * t321 + t137 * t317 + t403 * t243 + t273 + t399 + t409) * t359;
t214 = -rSges(4,2) * t249 + t343;
t146 = t291 * t245;
t144 = t291 * t243;
t122 = t90 / 0.2e1 + t242 * t376;
t100 = t301 - t322;
t97 = t274 * t245;
t95 = t274 * t243;
t93 = -t176 * t243 - t177 * t245;
t68 = -t297 * t194 - t277;
t42 = -t277 + (-t245 * t294 + t300) * t245 + t402 * t240;
t35 = (t387 * t245 + t341) * t245 + t310 + t417 * t243;
t25 = t363 + t372;
t24 = t365 / 0.2e1;
t22 = t366 / 0.2e1;
t19 = t169 * t35 + t346;
t15 = t348 + t350;
t5 = (t210 / 0.2e1 + t209 / 0.2e1) * t251 + t375 + t373 + t364 - t398;
t4 = t24 - t366 / 0.2e1;
t3 = t24 + t22;
t2 = t22 - t365 / 0.2e1;
t1 = t253 * t243 + t254 * t245;
t6 = [t5 * qJD(3) + t25 * qJD(4) + t37 * t344, 0, t5 * qJD(1) + t15 * qJD(4) + t3 * qJD(5) + ((t70 * t97 + t71 * t95 + (-t86 - t96) * t416) * t376 + (t144 * t89 + t146 * t88) * t377) * t378 + (m(4) * (-t176 * t212 - t214 * t98) + t272 * t245 + (-t153 / 0.2e1 - t172 / 0.2e1) * t251 + (t151 / 0.2e1 + t174 / 0.2e1) * t249 - t254) * t384 + ((m(4) * (t177 * t212 - t214 * t99) + t272 * t243 + (t154 / 0.2e1 - t173 / 0.2e1) * t251 + (-t152 / 0.2e1 - t175 / 0.2e1) * t249 - t253) * t243 + (-t289 / 0.2e1 + t287 / 0.2e1 + t288 / 0.2e1 + t286 / 0.2e1) * t242 + (-t285 / 0.2e1 - t283 / 0.2e1 + t408 * t359) * t244) * qJD(3), qJD(1) * t25 + qJD(3) * t15, t3 * qJD(3) + t295; 0, 0, t122 * qJD(5) + (m(4) * t93 / 0.2e1 + t68 * t377 + t42 * t376) * t378, 0, t122 * qJD(3); t1 * qJD(3) - t13 * qJD(4) + t4 * qJD(5) + (-t364 / 0.4e1 - t373 / 0.4e1 - t375 / 0.4e1) * t379 + (-(t210 + t209) * t251 / 0.2e1 + t398) * qJD(1), qJD(5) * t123, t1 * qJD(1) + t19 * t344 - ((t255 * t243 + (-t289 + t287 + t385) * t244 - t386 * t242 + (t256 + (t303 + t305) * t242 - t382) * t245) * t243 + t383 * t241) * t384 / 0.2e1 + (m(6) * (t35 * t42 - t416 * t95 - t96 * t97) + m(5) * (-t393 * t144 - t392 * t146 + (t243 * t261 + t245 * (rSges(5,1) * t317 + t290) + t310) * t68) + m(4) * (t297 * t214 * t212 + (t243 * (rSges(4,1) * t318 - t298) + t245 * (rSges(4,1) * t314 + t243 * rSges(4,3) - t220)) * t93) + ((t256 * t245 + (t285 + t283 - t386) * t242 + ((t306 - t307) * t245 + t385) * t244 + (t255 - t383) * t243) * t245 + t382 * t240) * t359) * qJD(3), -t413, t4 * qJD(1) + t296 + t19 * t345 + (-t169 * t244 - t100 + t301) * t344; t13 * qJD(3) - t90 * qJD(5) + (-t372 / 0.4e1 - t363 / 0.4e1) * t379, 0, t413 + ((t243 * t97 - t245 * t95) * t376 + (-t144 * t245 + t146 * t243) * t377) * t378, 0, -t313; t2 * qJD(3) + t90 * qJD(4) - t295, -t123 * qJD(3), t2 * qJD(1) - t296 + (-t244 * t42 + (t243 * t95 + t245 * t97 + t35) * t242 - t19 + t346) * t345 + t100 * t344, t313, t100 * t345;];
Cq = t6;
