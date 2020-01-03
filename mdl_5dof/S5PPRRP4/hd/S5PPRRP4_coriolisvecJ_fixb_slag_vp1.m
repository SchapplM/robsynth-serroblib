% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:27
% EndTime: 2019-12-31 17:34:39
% DurationCPUTime: 9.95s
% Computational Cost: add. (5646->364), mult. (13234->472), div. (0->0), fcn. (14102->6), ass. (0->204)
t400 = -Icges(5,3) - Icges(6,3);
t192 = sin(qJ(4));
t193 = cos(qJ(4));
t225 = Icges(6,5) * t193 - Icges(6,6) * t192;
t227 = Icges(5,5) * t193 - Icges(5,6) * t192;
t422 = t225 + t227;
t309 = sin(pkin(7));
t310 = cos(pkin(7));
t343 = sin(qJ(3));
t344 = cos(qJ(3));
t158 = -t309 * t343 - t310 * t344;
t159 = -t309 * t344 + t310 * t343;
t306 = Icges(6,4) * t192;
t233 = Icges(6,1) * t193 - t306;
t96 = -Icges(6,5) * t159 + t158 * t233;
t308 = Icges(5,4) * t192;
t235 = Icges(5,1) * t193 - t308;
t99 = -Icges(5,5) * t159 + t158 * t235;
t361 = t96 + t99;
t305 = Icges(6,4) * t193;
t229 = -Icges(6,2) * t192 + t305;
t90 = -Icges(6,6) * t159 + t158 * t229;
t307 = Icges(5,4) * t193;
t231 = -Icges(5,2) * t192 + t307;
t93 = -Icges(5,6) * t159 + t158 * t231;
t362 = t90 + t93;
t377 = -t192 * t362 + t193 * t361;
t384 = t422 * t158 + t400 * t159;
t404 = t377 * t158 - t384 * t159;
t385 = t400 * t158 - t422 * t159;
t424 = t385 * t159;
t100 = -Icges(5,5) * t158 - t159 * t235;
t97 = -Icges(6,5) * t158 - t159 * t233;
t413 = t100 + t97;
t423 = t413 * t193;
t228 = Icges(6,2) * t193 + t306;
t230 = Icges(5,2) * t193 + t308;
t421 = t228 + t230;
t232 = Icges(6,1) * t192 + t305;
t234 = Icges(5,1) * t192 + t307;
t420 = t232 + t234;
t358 = -t421 * t192 + t420 * t193;
t226 = Icges(5,5) * t192 + Icges(5,6) * t193;
t382 = t226 * t159;
t224 = Icges(6,5) * t192 + Icges(6,6) * t193;
t383 = t224 * t159;
t397 = -t158 * t358 + t382 + t383;
t418 = t397 * qJD(3);
t91 = -Icges(6,6) * t158 - t159 * t229;
t94 = -Icges(5,6) * t158 - t159 * t231;
t417 = t91 + t94;
t364 = -t192 * t361 - t193 * t362;
t401 = Icges(5,6) + Icges(6,6);
t402 = Icges(5,5) + Icges(6,5);
t416 = t402 * t192 + t401 * t193;
t410 = Icges(5,4) + Icges(6,4);
t415 = t410 * t192 + (Icges(5,2) + Icges(6,2)) * t193;
t414 = (Icges(5,1) + Icges(6,1)) * t192 + t410 * t193;
t408 = (t229 + t231) * qJD(4);
t407 = (t233 + t235) * qJD(4);
t406 = -t158 * t423 + t424;
t291 = t159 * t193;
t354 = -t385 * t158 - t413 * t291;
t292 = t159 * t192;
t390 = -t417 * t292 - t354;
t405 = t384 * t158 + t361 * t291 - t362 * t292;
t296 = t158 * t192;
t389 = -t417 * t296 - t406;
t112 = t224 * t158;
t114 = t226 * t158;
t403 = -t358 * t159 - t112 - t114;
t399 = t417 * t192;
t191 = -qJ(5) - pkin(6);
t396 = rSges(6,3) - t191;
t151 = qJD(3) * t158;
t150 = t159 * qJD(3);
t274 = qJD(4) * t192;
t269 = t158 * t274;
t217 = -t150 * t193 + t269;
t273 = qJD(4) * t193;
t218 = t150 * t192 + t158 * t273;
t395 = (t400 * t151 + t402 * t217 + t401 * t218) * t159;
t394 = t422 * qJD(4);
t393 = -t407 * t193 + t408 * t192 + (t420 * t192 + t421 * t193) * qJD(4);
t392 = -t224 - t226;
t376 = t399 - t423;
t391 = rSges(6,1) + pkin(4);
t268 = t159 * t274;
t297 = t151 * t193;
t215 = t268 + t297;
t298 = t151 * t192;
t216 = t159 * t273 - t298;
t386 = t400 * t150 + t402 * t215 + t401 * t216;
t62 = rSges(5,1) * t215 + rSges(5,2) * t216 - t150 * rSges(5,3);
t379 = (t389 * t158 + t404 * t159) * qJD(4);
t378 = (t390 * t158 + t405 * t159) * qJD(4);
t375 = t403 * qJD(3);
t276 = qJD(4) * t159;
t374 = 0.2e1 * qJD(4);
t372 = t375 + t378;
t371 = t379 + t418;
t370 = -t376 * qJD(4) - t416 * t150 + t414 * t215 + t415 * t216;
t369 = t377 * qJD(4) + t416 * t151 - t414 * t217 - t415 * t218;
t368 = -t150 * t358 + t151 * t392 + t158 * t393 + t159 * t394;
t367 = t150 * t392 + t151 * t358 - t158 * t394 + t159 * t393;
t365 = t413 * t192 + t417 * t193;
t247 = rSges(5,1) * t192 + rSges(5,2) * t193;
t275 = qJD(4) * t247;
t359 = t159 * t275;
t156 = t159 * pkin(6);
t135 = -pkin(3) * t158 + t156;
t341 = pkin(4) * t193;
t189 = pkin(3) + t341;
t357 = -rSges(6,1) * t291 - t158 * t396 - t159 * t189;
t356 = rSges(6,1) * t297 - t150 * t396 + t151 * t189;
t284 = -pkin(4) * t269 - qJD(5) * t159;
t355 = t151 * t396 + t284;
t277 = qJD(4) * t158;
t153 = t158 * qJD(5);
t260 = t151 * pkin(3) - pkin(6) * t150;
t147 = t151 * pkin(6);
t129 = -pkin(3) * t150 - t147;
t330 = rSges(6,2) * t192;
t246 = rSges(6,1) * t193 - t330;
t167 = t246 * qJD(4);
t329 = rSges(6,2) * t193;
t245 = rSges(6,1) * t192 + t329;
t261 = pkin(4) * t192 + t245;
t270 = qJD(4) ^ 2 * t341;
t338 = pkin(3) - t189;
t337 = rSges(6,1) * t217 + rSges(6,2) * t218 + t150 * t338 + t147 - t355;
t14 = t159 * t270 - qJD(5) * t150 + (-t151 * t261 + t159 * t167) * qJD(4) + (-t129 - t337) * qJD(3);
t258 = qJD(2) * t310;
t349 = -rSges(6,3) * t159 + (t246 - t338) * t158;
t317 = -t159 * t191 - t156 - t349;
t27 = -t258 - t153 + t261 * t276 + (-t135 - t317) * qJD(3);
t353 = t14 * t158 + t150 * t27;
t248 = rSges(5,1) * t193 - rSges(5,2) * t192;
t327 = rSges(5,3) * t159;
t103 = t158 * t248 - t327;
t157 = t159 * pkin(3);
t259 = -pkin(6) * t158 - t157;
t131 = qJD(3) * t259;
t318 = rSges(6,2) * t292 - t259 + t357;
t352 = t318 * qJD(3) + t131 - t284;
t288 = t230 * t159 + t100;
t312 = t234 * t159 - t94;
t351 = t192 * t288 - t193 * t312;
t314 = t232 * t159 - t91;
t316 = t228 * t159 + t97;
t350 = t192 * t316 - t193 * t314;
t336 = rSges(6,2) * t216 + t268 * t391 - t153 - t260 + t356;
t331 = m(4) * qJD(3);
t328 = rSges(5,3) * t151;
t319 = (pkin(6) + t191) * t159 + t349;
t315 = t228 * t158 - t96;
t313 = -t232 * t158 - t90;
t311 = -t234 * t158 - t93;
t287 = t230 * t158 - t99;
t281 = -t228 + t233;
t280 = -t229 - t232;
t279 = -t230 + t235;
t278 = -t231 - t234;
t272 = t225 * qJD(3);
t271 = t227 * qJD(3);
t264 = t277 / 0.2e1;
t263 = -t276 / 0.2e1;
t105 = -rSges(5,1) * t291 + rSges(5,2) * t292 - t158 * rSges(5,3);
t257 = -qJD(3) * t105 - t131;
t256 = pkin(4) * t292 + t245 * t159;
t255 = pkin(4) * t296 + t245 * t158;
t190 = qJD(2) * t309;
t26 = t245 * t277 + t190 + t352;
t250 = t26 * t261;
t249 = rSges(4,1) * t158 + rSges(4,2) * t159;
t39 = t158 * t275 + t190 - t257;
t40 = t359 - t258 + (t103 - t135) * qJD(3);
t244 = -t158 * t39 - t159 * t40;
t223 = -t103 * t158 + t105 * t159;
t220 = pkin(3) + t248;
t209 = t192 * t315 + t193 * t313;
t208 = t192 * t287 + t193 * t311;
t207 = (t192 * t280 + t193 * t281) * qJD(3);
t206 = (t192 * t278 + t193 * t279) * qJD(3);
t60 = rSges(5,1) * t217 + rSges(5,2) * t218 - t328;
t205 = -t103 * t150 - t105 * t151 + t158 * t60 + t159 * t62;
t168 = t248 * qJD(4);
t132 = -rSges(4,1) * t159 + rSges(4,2) * t158;
t130 = qJD(3) * t135;
t128 = t247 * t158;
t126 = t247 * t159;
t124 = rSges(4,1) * t151 + rSges(4,2) * t150;
t123 = -rSges(4,1) * t150 + rSges(4,2) * t151;
t110 = qJD(3) * t260;
t109 = qJD(3) * t249 - t258;
t108 = qJD(3) * t132 + t190;
t38 = qJD(4) * t223 + qJD(1);
t37 = qJD(3) * t62 + t110 + (t150 * t247 + t158 * t168) * qJD(4);
t36 = (-t151 * t247 + t159 * t168) * qJD(4) + (-t129 - t60) * qJD(3);
t21 = qJD(1) + (t158 * t317 + t159 * t318) * qJD(4);
t16 = t205 * qJD(4);
t15 = t158 * t270 - qJD(5) * t151 + t110 + t336 * qJD(3) + (t150 * t261 + t158 * t167) * qJD(4);
t1 = (t317 * t150 - t318 * t151 + t337 * t158 + t336 * t159) * qJD(4);
t2 = [m(5) * t16 + m(6) * t1; (t123 * t310 + t124 * t309) * t331 + m(5) * (t309 * t37 - t310 * t36) + m(6) * (-t14 * t310 + t15 * t309); m(4) * (t108 * t124 - t109 * t123) - (t108 * t249 - t109 * t132) * t331 + t372 * t276 / 0.2e1 + (m(4) * (-t123 * t249 + t124 * t132) + t408 * t193 + t407 * t192 + t358 * qJD(4)) * qJD(3) + (t15 * t357 + (-t14 * t396 + t15 * t330) * t159 + t353 * (t189 + t246) + (t352 + t355) * t27 + (-rSges(6,2) * t298 - t319 * qJD(3) + t130 + t356) * t26 + (t26 * (t192 * t391 + t329) - t250 - t21 * (t317 + t319)) * t276) * m(6) + (t37 * (t105 - t157) + t36 * (-t156 - t327) + (-t37 * pkin(6) + t220 * t36) * t158 + (-qJD(3) * t103 + t130 + t260 - t359 + t62) * t39 + (t150 * t220 + t147 - t257 + t328) * t40) * m(5) + (((t405 + t389 + t406) * t159 - t354 * t158) * qJD(4) - t368 + t369 + t375) * t263 + ((-t365 - t403) * t150 + (-t364 - t397) * t151) * qJD(4) / 0.2e1 + ((((t384 + t399) * t159 + t354 + t390) * t159 + ((t384 + t376) * t158 + t424 - t405) * t158) * qJD(4) - t367 - t370 + t371 - t418) * t264; -(t365 * t150 + t364 * t151 + t370 * t158 + t369 * t159) * qJD(3) / 0.2e1 + ((((-t287 - t315) * t159 + (t288 + t316) * t158) * t193 + ((t311 + t313) * t159 + (t312 + t314) * t158) * t192) * qJD(4) + ((t278 + t280) * t193 + (-t279 - t281) * t192) * qJD(3)) * qJD(3) / 0.2e1 + ((t277 * t383 - t272) * t158 + (-t207 + (t209 * t159 + (-t112 - t350) * t158) * qJD(4)) * t159 + (t277 * t382 - t271) * t158 + (-t206 + (t208 * t159 + (-t114 - t351) * t158) * qJD(4)) * t159) * t264 + ((t112 * t276 + t272) * t159 + (-t207 + (-t350 * t158 + (-t383 + t209) * t159) * qJD(4)) * t158 + (t114 * t276 + t271) * t159 + (-t206 + (-t351 * t158 + (-t382 + t208) * t159) * qJD(4)) * t158) * t263 + (-(t126 * t39 - t128 * t40) * qJD(3) - (t38 * (t126 * t159 + t128 * t158) - t244 * t248) * qJD(4) + t16 * t223 + t38 * t205 - t244 * t168 - (-t150 * t39 + t151 * t40 - t158 * t37 - t159 * t36) * t247) * m(5) - (t367 * qJD(3) + (-t405 * t151 + t390 * t150 + (t384 * t150 - t377 * t151) * t159 + (t385 * t150 + t376 * t151 + t386 * t158 - t395) * t158) * t374) * t158 / 0.2e1 + (t368 * qJD(3) + (-t404 * t151 + t389 * t150 + (t377 * t150 + t384 * t151 + t395) * t159 + (-t376 * t150 + t385 * t151 - t159 * t386) * t158) * t374) * t159 / 0.2e1 + (-(t21 * t318 + t261 * t27) * t151 + (t21 * t317 + t250) * t150 + (t1 * t318 + t14 * t261 + t27 * t167 + t21 * t336) * t159 + (t15 * t261 + t26 * (pkin(4) * t273 + t167) + t1 * t317 + t21 * t337) * t158 - (-t255 * t27 + t256 * t26) * qJD(3) - ((t21 * t256 + t246 * t27) * t159 + (t26 * (t246 + t341) + t255 * t21) * t158) * qJD(4)) * m(6) - (t372 + t378) * t150 / 0.2e1 - (t371 + t379) * t151 / 0.2e1; (t15 * t159 - t151 * t26 - (-t158 * t26 - t159 * t27) * qJD(3) - t353) * m(6);];
tauc = t2(:);
