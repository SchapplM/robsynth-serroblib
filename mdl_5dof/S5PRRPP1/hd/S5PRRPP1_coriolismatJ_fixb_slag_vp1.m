% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:13
% EndTime: 2019-12-05 16:06:29
% DurationCPUTime: 6.12s
% Computational Cost: add. (11152->313), mult. (9000->420), div. (0->0), fcn. (8073->6), ass. (0->199)
t409 = Icges(4,3) + Icges(5,3);
t245 = pkin(7) + qJ(2);
t240 = sin(t245);
t247 = sin(qJ(3));
t346 = pkin(3) * t247;
t246 = qJ(3) + pkin(8);
t243 = cos(t246);
t241 = sin(t246);
t347 = rSges(6,1) + pkin(4);
t284 = t347 * t241;
t332 = rSges(6,3) + qJ(5);
t393 = t332 * t243 - t284;
t271 = t346 - t393;
t406 = t271 * t240;
t408 = t240 * t406;
t242 = cos(t245);
t379 = t332 * t241 + t347 * t243;
t407 = -t242 * rSges(6,2) + t240 * t379;
t182 = Icges(5,5) * t243 - Icges(5,6) * t241;
t248 = cos(qJ(3));
t207 = Icges(4,5) * t248 - Icges(4,6) * t247;
t405 = (t182 + t207) * t242 + t409 * t240;
t228 = Icges(6,5) * t241;
t329 = Icges(6,1) * t243;
t265 = t228 + t329;
t138 = Icges(6,4) * t240 + t265 * t242;
t327 = Icges(5,4) * t241;
t189 = Icges(5,1) * t243 - t327;
t140 = Icges(5,5) * t240 + t189 * t242;
t404 = t138 + t140;
t314 = t240 * t248;
t315 = t240 * t247;
t316 = t240 * t243;
t317 = t240 * t241;
t400 = -Icges(4,5) * t314 - Icges(5,5) * t316 + Icges(4,6) * t315 + Icges(5,6) * t317 + t409 * t242;
t367 = m(6) / 0.2e1;
t193 = rSges(5,1) * t241 + rSges(5,2) * t243;
t253 = t193 + t346;
t384 = t253 * t242;
t385 = t253 * t240;
t388 = t240 * t385 + t242 * t384;
t96 = t271 * t242;
t341 = (-t242 * t96 - t408) * t367 - m(5) * t388 / 0.2e1;
t368 = m(5) / 0.2e1;
t309 = t242 * t243;
t290 = t332 * t309;
t85 = (-t284 - t346) * t242 + t290;
t343 = (-t242 * t85 + t408) * t367 + t388 * t368;
t13 = t343 - t341;
t403 = t13 * qJD(2);
t328 = Icges(4,4) * t247;
t211 = Icges(4,1) * t248 - t328;
t154 = Icges(4,5) * t240 + t211 * t242;
t402 = -t140 * t316 - t154 * t314;
t399 = t405 * t242 + t402;
t202 = Icges(5,4) * t317;
t139 = Icges(5,1) * t316 - Icges(5,5) * t242 - t202;
t216 = Icges(4,4) * t315;
t153 = Icges(4,1) * t314 - Icges(4,5) * t242 - t216;
t308 = t242 * t248;
t398 = -t139 * t309 - t153 * t308 + t400 * t240;
t135 = Icges(5,4) * t316 - Icges(5,2) * t317 - Icges(5,6) * t242;
t151 = Icges(4,4) * t314 - Icges(4,2) * t315 - Icges(4,6) * t242;
t307 = t247 * t151;
t397 = t135 * t241 + t307;
t201 = Icges(6,5) * t309;
t313 = t241 * t242;
t130 = Icges(6,6) * t240 + Icges(6,3) * t313 + t201;
t183 = Icges(6,4) * t243 + Icges(6,6) * t241;
t134 = Icges(6,2) * t240 + t183 * t242;
t396 = t130 * t313 + t154 * t308 + t404 * t309 + (t134 + t405) * t240;
t395 = -Icges(4,5) * t247 - Icges(4,6) * t248 + (-Icges(5,6) + Icges(6,6)) * t243 + (-Icges(6,4) - Icges(5,5)) * t241;
t238 = t240 ^ 2;
t239 = t242 ^ 2;
t287 = t238 + t239;
t232 = Icges(5,4) * t243;
t324 = Icges(5,2) * t241;
t136 = Icges(5,6) * t240 + (t232 - t324) * t242;
t244 = Icges(4,4) * t248;
t209 = -Icges(4,2) * t247 + t244;
t152 = Icges(4,6) * t240 + t209 * t242;
t306 = t247 * t152;
t394 = t136 * t241 + t306 + t400;
t392 = -t136 * t313 - t242 * t306 + t396;
t321 = (-Icges(6,2) * t242 + t183 * t240) * t242;
t391 = t321 + t396;
t390 = -t135 * t313 - t136 * t317 - t240 * t306 - t242 * t307 - t398 - t399;
t326 = Icges(6,5) * t243;
t181 = Icges(6,3) * t241 + t326;
t184 = Icges(5,2) * t243 + t327;
t186 = Icges(6,1) * t241 - t326;
t208 = Icges(4,2) * t248 + t328;
t323 = Icges(6,3) * t243;
t330 = Icges(5,1) * t241;
t389 = -(t211 / 0.2e1 - t208 / 0.2e1) * t247 - (t232 + t330 / 0.2e1 - t324 / 0.2e1 + t186 / 0.2e1 - t181 / 0.2e1) * t243 - (t189 / 0.2e1 - t184 / 0.2e1 + t228 + t329 / 0.2e1 - t323 / 0.2e1) * t241;
t387 = -t240 / 0.2e1;
t350 = t240 / 0.2e1;
t349 = -t242 / 0.2e1;
t129 = -Icges(6,6) * t242 + t181 * t240;
t137 = -Icges(6,4) * t242 + t265 * t240;
t381 = (t129 * t241 + t137 * t243) * t240;
t210 = Icges(4,1) * t247 + t244;
t261 = t323 - t228;
t378 = (-t184 - t261) * t242 + t404;
t377 = -Icges(5,2) * t316 - t261 * t240 + t137 + t139 - t202;
t376 = t395 * t240;
t375 = t395 * t242;
t173 = t208 * t242;
t175 = t210 * t242;
t373 = -t241 * t378 + (-t152 - t175) * t248 + (-t154 + t173) * t247;
t172 = -Icges(4,2) * t314 - t216;
t174 = t210 * t240;
t372 = t241 * t377 + (t151 + t174) * t248 + (t153 + t172) * t247;
t370 = 0.4e1 * qJD(2);
t369 = 2 * qJD(3);
t236 = t242 * pkin(6);
t336 = rSges(4,1) * t248;
t283 = pkin(2) + t336;
t288 = rSges(4,2) * t315 + t242 * rSges(4,3);
t112 = -t283 * t240 + t236 + t288;
t333 = t247 * rSges(4,2);
t218 = t242 * t333;
t113 = -t218 + t283 * t242 + (rSges(4,3) + pkin(6)) * t240;
t212 = t247 * rSges(4,1) + rSges(4,2) * t248;
t176 = t212 * t240;
t177 = t212 * t242;
t366 = m(4) * (t112 * t176 - t113 * t177);
t257 = rSges(5,1) * t316 - rSges(5,2) * t317 - t242 * rSges(5,3);
t345 = pkin(3) * t248;
t237 = pkin(2) + t345;
t344 = -qJ(4) - pkin(6);
t289 = -t240 * t237 - t242 * t344;
t92 = -t257 + t289;
t219 = t240 * t344;
t281 = -rSges(5,2) * t313 + rSges(5,3) * t240;
t335 = rSges(5,1) * t243;
t93 = -t219 + (t237 + t335) * t242 + t281;
t364 = m(5) * (-t384 * t93 + t385 * t92);
t363 = m(5) * (t93 * t240 + t242 * t92);
t72 = t289 - t407;
t334 = rSges(6,2) * t240;
t73 = t334 - t219 + (t237 + t379) * t242;
t340 = t72 * t309 + t73 * t316;
t357 = m(6) * ((t240 * t85 + t242 * t406) * t241 + t340);
t356 = m(6) * (-t313 * t406 + t96 * t317 + t340);
t355 = m(6) * (t406 * t72 + t73 * t85);
t354 = m(6) * (t73 * t240 + t242 * t72);
t169 = t287 * t241;
t88 = m(6) * t169;
t339 = -t96 * t309 - t316 * t406;
t338 = m(6) * qJD(3);
t337 = m(6) * qJD(5);
t312 = t241 * t243;
t305 = t88 * qJD(2);
t302 = -t240 * (t240 * pkin(2) - t236 + t289) + t242 * (-t240 * pkin(6) - t219 + (-pkin(2) + t237) * t242);
t299 = -t186 * t240 + t129;
t298 = -Icges(6,1) * t313 + t130 + t201;
t266 = -t232 - t330;
t297 = -t266 * t240 + t135;
t296 = t266 * t242 - t136;
t291 = t287 * t312;
t123 = (t169 / 0.2e1 - t241 / 0.2e1) * m(6);
t286 = t123 * qJD(1);
t37 = t73 * t313 - t72 * t317;
t285 = m(6) * t37 * qJD(2);
t282 = rSges(5,2) * t241 - t335 - t345;
t280 = t299 * t242;
t279 = t298 * t240;
t278 = t297 * t242;
t277 = t296 * t240;
t272 = t287 * t346;
t270 = -t345 - t379;
t269 = -t130 * t317 + t134 * t242 - t138 * t316;
t268 = t207 / 0.2e1 + t183 / 0.2e1 + t182 / 0.2e1;
t45 = -t321 + t381;
t250 = (t45 - t381 + t391) * t387 + t392 * t350 + ((t397 + t405) * t242 + t390 + t398 + t402) * t349;
t249 = t269 * t387 + (t394 * t242 - t391 + t392) * t242 / 0.2e1 + (t45 + t400 * t242 + (t139 * t243 + t153 * t248 - t397) * t240) * t349 + (t129 * t313 + t137 * t309 + t394 * t240 + t269 + t390 + t399) * t350;
t213 = -t333 + t336;
t146 = t282 * t242;
t144 = t282 * t240;
t122 = t88 / 0.2e1 + t241 * t367;
t98 = t291 - t312;
t97 = t270 * t242;
t95 = t270 * t240;
t91 = -t176 * t240 - t177 * t242;
t68 = -t193 * t287 - t272;
t42 = -t272 + (-t242 * t284 + t290) * t242 + t393 * t238;
t35 = (t242 * t379 + t334) * t242 + t302 + t407 * t240;
t26 = t354 + t363;
t24 = t356 / 0.2e1;
t22 = t357 / 0.2e1;
t19 = t169 * t35 + t339;
t15 = t341 + t343;
t5 = (t210 / 0.2e1 + t209 / 0.2e1) * t248 + t366 + t364 + t355 - t389;
t4 = t24 - t357 / 0.2e1;
t3 = t24 + t22;
t2 = t22 - t356 / 0.2e1;
t1 = t249 * t240 + t250 * t242;
t6 = [0, 0, t122 * qJD(5) + (m(4) * t91 / 0.2e1 + t68 * t368 + t42 * t367) * t369, 0, t122 * qJD(3); 0, t5 * qJD(3) + t26 * qJD(4) + t37 * t337, t5 * qJD(2) + t15 * qJD(4) + t3 * qJD(5) + ((t72 * t97 + t73 * t95 + (-t85 - t96) * t406) * t367 + (t144 * t93 + t146 * t92) * t368) * t369 + ((m(4) * (-t112 * t213 - t176 * t212) + t268 * t242 + (-t153 / 0.2e1 - t172 / 0.2e1) * t248 + (t151 / 0.2e1 + t174 / 0.2e1) * t247 - t250) * t242 + (m(4) * (-t113 * t213 + t177 * t212) + t268 * t240 + (t154 / 0.2e1 - t173 / 0.2e1) * t248 + (-t152 / 0.2e1 - t175 / 0.2e1) * t247 - t249) * t240 + (-t280 / 0.2e1 + t278 / 0.2e1 + t277 / 0.2e1 + t279 / 0.2e1) * t241 + (t349 * t377 + t350 * t378) * t243) * qJD(3), qJD(2) * t26 + qJD(3) * t15, t3 * qJD(3) + t285; qJD(5) * t123, t1 * qJD(3) - t13 * qJD(4) + t4 * qJD(5) + (-t355 / 0.4e1 - t364 / 0.4e1 - t366 / 0.4e1) * t370 + (-(t210 + t209) * t248 / 0.2e1 + t389) * qJD(2), t1 * qJD(2) + t19 * t337 + (m(4) * (t287 * t213 * t212 + (t240 * (rSges(4,1) * t314 - t288) + t242 * (rSges(4,1) * t308 + t240 * rSges(4,3) - t218)) * t91) + m(6) * (t35 * t42 - t406 * t95 - t96 * t97) + m(5) * (-t385 * t144 - t384 * t146 + (t240 * t257 + t242 * (rSges(5,1) * t309 + t281) + t302) * t68) + (t375 * t238 + (((t297 - t299) * t243 + t372) * t242 + ((t296 + t298) * t243 + t373 - t376) * t240) * t242) * t350 + ((t373 * t240 + (t278 + t277 - t280 + t279) * t243 + (t372 - t375) * t242) * t240 + t376 * t239) * t349) * qJD(3), -t403, t286 + t4 * qJD(2) + t19 * t338 + (-t169 * t243 + t291 - t98) * t337; 0, t13 * qJD(3) - t88 * qJD(5) + (-t363 / 0.4e1 - t354 / 0.4e1) * t370, t403 + ((t240 * t97 - t242 * t95) * t367 + (-t144 * t242 + t146 * t240) * t368) * t369, 0, -t305; -t123 * qJD(3), t2 * qJD(3) + t88 * qJD(4) - t285, -t286 + t2 * qJD(2) + (-t243 * t42 + (t240 * t95 + t242 * t97 + t35) * t241 - t19 + t339) * t338 + t98 * t337, t305, t98 * t338;];
Cq = t6;
