% Calculate time derivative of joint inertia matrix for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:22
% EndTime: 2019-12-05 16:26:09
% DurationCPUTime: 17.77s
% Computational Cost: add. (45676->886), mult. (88918->1271), div. (0->0), fcn. (99783->12), ass. (0->364)
t312 = sin(pkin(9));
t314 = cos(pkin(9));
t322 = cos(qJ(2));
t315 = cos(pkin(5));
t319 = sin(qJ(2));
t378 = t315 * t319;
t299 = t312 * t322 + t314 * t378;
t353 = qJ(3) + pkin(10);
t311 = sin(t353);
t340 = cos(t353);
t313 = sin(pkin(5));
t383 = t313 * t314;
t273 = t299 * t340 - t311 * t383;
t332 = t313 * t340;
t323 = -t299 * t311 - t314 * t332;
t376 = t315 * t322;
t327 = -t312 * t319 + t314 * t376;
t188 = Icges(5,5) * t273 + Icges(5,6) * t323 - Icges(5,3) * t327;
t318 = sin(qJ(3));
t321 = cos(qJ(3));
t381 = t313 * t321;
t278 = -t299 * t318 - t314 * t381;
t382 = t313 * t318;
t350 = t314 * t382;
t328 = -t299 * t321 + t350;
t206 = -Icges(4,5) * t328 + Icges(4,6) * t278 - Icges(4,3) * t327;
t414 = t188 + t206;
t351 = t312 * t378;
t301 = t314 * t322 - t351;
t384 = t312 * t313;
t275 = t301 * t340 + t311 * t384;
t300 = t312 * t376 + t314 * t319;
t324 = -t301 * t311 + t312 * t332;
t189 = Icges(5,5) * t275 + Icges(5,6) * t324 + Icges(5,3) * t300;
t280 = -t301 * t318 + t312 * t381;
t352 = t312 * t382;
t281 = t301 * t321 + t352;
t207 = Icges(4,5) * t281 + Icges(4,6) * t280 + Icges(4,3) * t300;
t413 = t189 + t207;
t190 = Icges(5,4) * t273 + Icges(5,2) * t323 - Icges(5,6) * t327;
t192 = Icges(5,1) * t273 + Icges(5,4) * t323 - Icges(5,5) * t327;
t208 = -Icges(4,4) * t328 + Icges(4,2) * t278 - Icges(4,6) * t327;
t210 = -Icges(4,1) * t328 + Icges(4,4) * t278 - Icges(4,5) * t327;
t410 = t190 * t323 + t192 * t273 + t208 * t278 - t210 * t328 - t327 * t414;
t191 = Icges(5,4) * t275 + Icges(5,2) * t324 + Icges(5,6) * t300;
t193 = Icges(5,1) * t275 + Icges(5,4) * t324 + Icges(5,5) * t300;
t209 = Icges(4,4) * t281 + Icges(4,2) * t280 + Icges(4,6) * t300;
t211 = Icges(4,1) * t281 + Icges(4,4) * t280 + Icges(4,5) * t300;
t409 = t191 * t323 + t193 * t273 + t209 * t278 - t211 * t328 - t327 * t413;
t408 = t190 * t324 + t192 * t275 + t208 * t280 + t210 * t281 + t300 * t414;
t407 = t191 * t324 + t193 * t275 + t209 * t280 + t211 * t281 + t300 * t413;
t286 = t311 * t313 * t319 - t315 * t340;
t287 = t311 * t315 + t319 * t332;
t380 = t313 * t322;
t236 = Icges(5,5) * t287 - Icges(5,6) * t286 - Icges(5,3) * t380;
t237 = Icges(5,4) * t287 - Icges(5,2) * t286 - Icges(5,6) * t380;
t238 = Icges(5,1) * t287 - Icges(5,4) * t286 - Icges(5,5) * t380;
t124 = -t236 * t327 + t237 * t323 + t238 * t273;
t375 = t318 * t319;
t377 = t315 * t321;
t302 = -t313 * t375 + t377;
t379 = t315 * t318;
t303 = t319 * t381 + t379;
t250 = Icges(4,5) * t303 + Icges(4,6) * t302 - Icges(4,3) * t380;
t251 = Icges(4,4) * t303 + Icges(4,2) * t302 - Icges(4,6) * t380;
t252 = Icges(4,1) * t303 + Icges(4,4) * t302 - Icges(4,5) * t380;
t130 = -t250 * t327 + t251 * t278 - t252 * t328;
t412 = t124 + t130;
t125 = t236 * t300 + t237 * t324 + t238 * t275;
t131 = t250 * t300 + t251 * t280 + t252 * t281;
t411 = t125 + t131;
t355 = qJD(2) * t319;
t342 = t313 * t355;
t115 = -t188 * t380 - t190 * t286 + t192 * t287;
t121 = -t206 * t380 + t208 * t302 + t210 * t303;
t406 = t115 + t121;
t116 = -t189 * t380 - t191 * t286 + t193 * t287;
t122 = -t207 * t380 + t209 * t302 + t211 * t303;
t405 = t116 + t122;
t290 = t299 * qJD(2);
t354 = qJD(2) * t322;
t292 = -qJD(2) * t351 + t314 * t354;
t317 = sin(qJ(5));
t320 = cos(qJ(5));
t228 = -t273 * t317 - t320 * t327;
t229 = t273 * t320 - t317 * t327;
t138 = Icges(6,5) * t229 + Icges(6,6) * t228 - Icges(6,3) * t323;
t140 = Icges(6,4) * t229 + Icges(6,2) * t228 - Icges(6,6) * t323;
t142 = Icges(6,1) * t229 + Icges(6,4) * t228 - Icges(6,5) * t323;
t289 = t327 * qJD(2);
t224 = qJD(3) * t323 + t289 * t340;
t159 = -qJD(5) * t229 - t224 * t317 + t290 * t320;
t160 = qJD(5) * t228 + t224 * t320 + t290 * t317;
t223 = qJD(3) * t273 + t289 * t311;
t95 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t223;
t97 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t223;
t99 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t223;
t22 = t138 * t223 + t140 * t159 + t142 * t160 + t228 * t97 + t229 * t99 - t323 * t95;
t291 = t300 * qJD(2);
t226 = qJD(3) * t324 - t291 * t340;
t231 = t275 * t320 + t300 * t317;
t161 = -qJD(5) * t231 - t226 * t317 + t292 * t320;
t230 = -t275 * t317 + t300 * t320;
t162 = qJD(5) * t230 + t226 * t320 + t292 * t317;
t225 = qJD(3) * t275 - t291 * t311;
t100 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t225;
t139 = Icges(6,5) * t231 + Icges(6,6) * t230 - Icges(6,3) * t324;
t141 = Icges(6,4) * t231 + Icges(6,2) * t230 - Icges(6,6) * t324;
t143 = Icges(6,1) * t231 + Icges(6,4) * t230 - Icges(6,5) * t324;
t96 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t225;
t98 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t225;
t23 = t100 * t229 + t139 * t223 + t141 * t159 + t143 * t160 + t228 * t98 - t323 * t96;
t271 = -qJD(3) * t286 + t332 * t354;
t329 = -t287 * t320 + t317 * t380;
t199 = qJD(5) * t329 - t271 * t317 + t320 * t342;
t276 = -t287 * t317 - t320 * t380;
t200 = qJD(5) * t276 + t271 * t320 + t317 * t342;
t341 = t313 * t354;
t270 = qJD(3) * t287 + t311 * t341;
t126 = Icges(6,5) * t200 + Icges(6,6) * t199 + Icges(6,3) * t270;
t127 = Icges(6,4) * t200 + Icges(6,2) * t199 + Icges(6,6) * t270;
t128 = Icges(6,1) * t200 + Icges(6,4) * t199 + Icges(6,5) * t270;
t183 = -Icges(6,5) * t329 + Icges(6,6) * t276 + Icges(6,3) * t286;
t184 = -Icges(6,4) * t329 + Icges(6,2) * t276 + Icges(6,6) * t286;
t185 = -Icges(6,1) * t329 + Icges(6,4) * t276 + Icges(6,5) * t286;
t40 = -t126 * t323 + t127 * t228 + t128 * t229 + t159 * t184 + t160 * t185 + t183 * t223;
t69 = -t138 * t323 + t140 * t228 + t142 * t229;
t70 = -t139 * t323 + t141 * t228 + t143 * t229;
t87 = -t183 * t323 + t184 * t228 + t185 * t229;
t3 = -t22 * t327 + t23 * t300 + t69 * t290 + t70 * t292 + (-t322 * t40 + t355 * t87) * t313;
t146 = Icges(5,5) * t224 - Icges(5,6) * t223 + Icges(5,3) * t290;
t148 = Icges(5,4) * t224 - Icges(5,2) * t223 + Icges(5,6) * t290;
t150 = Icges(5,1) * t224 - Icges(5,4) * t223 + Icges(5,5) * t290;
t50 = -t146 * t327 + t148 * t323 + t150 * t273 + t188 * t290 - t190 * t223 + t192 * t224;
t147 = Icges(5,5) * t226 - Icges(5,6) * t225 + Icges(5,3) * t292;
t149 = Icges(5,4) * t226 - Icges(5,2) * t225 + Icges(5,6) * t292;
t151 = Icges(5,1) * t226 - Icges(5,4) * t225 + Icges(5,5) * t292;
t51 = -t147 * t327 + t149 * t323 + t151 * t273 + t189 * t290 - t191 * t223 + t193 * t224;
t232 = qJD(3) * t328 - t289 * t318;
t326 = t278 * qJD(3);
t233 = t289 * t321 + t326;
t166 = Icges(4,5) * t233 + Icges(4,6) * t232 + Icges(4,3) * t290;
t168 = Icges(4,4) * t233 + Icges(4,2) * t232 + Icges(4,6) * t290;
t170 = Icges(4,1) * t233 + Icges(4,4) * t232 + Icges(4,5) * t290;
t54 = -t166 * t327 + t168 * t278 - t170 * t328 + t206 * t290 + t208 * t232 + t210 * t233;
t234 = -qJD(3) * t281 + t291 * t318;
t325 = t280 * qJD(3);
t235 = -t291 * t321 + t325;
t167 = Icges(4,5) * t235 + Icges(4,6) * t234 + Icges(4,3) * t292;
t169 = Icges(4,4) * t235 + Icges(4,2) * t234 + Icges(4,6) * t292;
t171 = Icges(4,1) * t235 + Icges(4,4) * t234 + Icges(4,5) * t292;
t55 = -t167 * t327 + t169 * t278 - t171 * t328 + t207 * t290 + t209 * t232 + t211 * t233;
t201 = Icges(5,5) * t271 - Icges(5,6) * t270 + Icges(5,3) * t342;
t202 = Icges(5,4) * t271 - Icges(5,2) * t270 + Icges(5,6) * t342;
t203 = Icges(5,1) * t271 - Icges(5,4) * t270 + Icges(5,5) * t342;
t67 = -t201 * t327 + t202 * t323 + t203 * t273 - t223 * t237 + t224 * t238 + t236 * t290;
t282 = -qJD(3) * t303 - t318 * t341;
t283 = qJD(3) * t302 + t321 * t341;
t217 = Icges(4,5) * t283 + Icges(4,6) * t282 + Icges(4,3) * t342;
t218 = Icges(4,4) * t283 + Icges(4,2) * t282 + Icges(4,6) * t342;
t219 = Icges(4,1) * t283 + Icges(4,4) * t282 + Icges(4,5) * t342;
t79 = -t217 * t327 + t218 * t278 - t219 * t328 + t232 * t251 + t233 * t252 + t250 * t290;
t404 = t3 + (-t50 - t54) * t327 + (t412 * t355 + (-t67 - t79) * t322) * t313 + (t51 + t55) * t300 + t409 * t292 + t410 * t290;
t24 = t138 * t225 + t140 * t161 + t142 * t162 + t230 * t97 + t231 * t99 - t324 * t95;
t25 = t100 * t231 + t139 * t225 + t141 * t161 + t143 * t162 + t230 * t98 - t324 * t96;
t41 = -t126 * t324 + t127 * t230 + t128 * t231 + t161 * t184 + t162 * t185 + t183 * t225;
t71 = -t138 * t324 + t140 * t230 + t142 * t231;
t72 = -t139 * t324 + t141 * t230 + t143 * t231;
t88 = -t183 * t324 + t184 * t230 + t185 * t231;
t4 = -t24 * t327 + t25 * t300 + t71 * t290 + t72 * t292 + (-t322 * t41 + t355 * t88) * t313;
t52 = t146 * t300 + t148 * t324 + t150 * t275 + t188 * t292 - t190 * t225 + t192 * t226;
t53 = t147 * t300 + t149 * t324 + t151 * t275 + t189 * t292 - t191 * t225 + t193 * t226;
t56 = t166 * t300 + t168 * t280 + t170 * t281 + t206 * t292 + t208 * t234 + t210 * t235;
t57 = t167 * t300 + t169 * t280 + t171 * t281 + t207 * t292 + t209 * t234 + t211 * t235;
t68 = t201 * t300 + t202 * t324 + t203 * t275 - t225 * t237 + t226 * t238 + t236 * t292;
t80 = t217 * t300 + t218 * t280 + t219 * t281 + t234 * t251 + t235 * t252 + t250 * t292;
t403 = t4 + (-t52 - t56) * t327 + (t411 * t355 + (-t68 - t80) * t322) * t313 + (t53 + t57) * t300 + t407 * t292 + t408 * t290;
t402 = 2 * m(4);
t401 = 2 * m(5);
t400 = 2 * m(6);
t399 = t223 / 0.2e1;
t398 = t225 / 0.2e1;
t397 = t270 / 0.2e1;
t396 = -t323 / 0.2e1;
t395 = -t324 / 0.2e1;
t394 = t286 / 0.2e1;
t393 = t290 / 0.2e1;
t392 = t292 / 0.2e1;
t391 = t312 / 0.2e1;
t390 = -t314 / 0.2e1;
t389 = pkin(3) * t321;
t387 = pkin(3) * qJD(3);
t386 = Icges(3,4) * t319;
t385 = Icges(3,4) * t322;
t101 = rSges(6,1) * t160 + rSges(6,2) * t159 + rSges(6,3) * t223;
t374 = pkin(4) * t224 + pkin(8) * t223 + t101;
t102 = rSges(6,1) * t162 + rSges(6,2) * t161 + rSges(6,3) * t225;
t373 = pkin(4) * t226 + pkin(8) * t225 + t102;
t129 = rSges(6,1) * t200 + rSges(6,2) * t199 + rSges(6,3) * t270;
t372 = pkin(4) * t271 + pkin(8) * t270 + t129;
t164 = pkin(3) * t326 + qJ(4) * t290 - qJD(4) * t327 + t289 * t389;
t196 = -pkin(3) * t350 - qJ(4) * t327 + t299 * t389;
t371 = t164 * t300 + t196 * t292;
t144 = rSges(6,1) * t229 + rSges(6,2) * t228 - rSges(6,3) * t323;
t370 = pkin(4) * t273 - pkin(8) * t323 + t144;
t145 = rSges(6,1) * t231 + rSges(6,2) * t230 - rSges(6,3) * t324;
t369 = pkin(4) * t275 - pkin(8) * t324 + t145;
t157 = rSges(5,1) * t226 - rSges(5,2) * t225 + rSges(5,3) * t292;
t165 = pkin(3) * t325 + qJ(4) * t292 + qJD(4) * t300 - t291 * t389;
t368 = -t157 - t165;
t266 = -pkin(2) * t291 + pkin(7) * t292;
t242 = t315 * t266;
t367 = t165 * t315 + t242;
t265 = pkin(2) * t289 + pkin(7) * t290;
t366 = -t164 - t265;
t253 = pkin(3) * t379 + (-qJ(4) * t322 + t319 * t389) * t313;
t365 = t196 * t380 - t253 * t327;
t186 = -rSges(6,1) * t329 + rSges(6,2) * t276 + rSges(6,3) * t286;
t364 = pkin(4) * t287 + pkin(8) * t286 + t186;
t197 = pkin(3) * t352 + qJ(4) * t300 + t301 * t389;
t269 = pkin(2) * t301 + pkin(7) * t300;
t267 = t315 * t269;
t363 = t197 * t315 + t267;
t194 = rSges(5,1) * t273 + rSges(5,2) * t323 - rSges(5,3) * t327;
t362 = -t194 - t196;
t195 = rSges(5,1) * t275 + rSges(5,2) * t324 + rSges(5,3) * t300;
t361 = -t195 - t197;
t204 = rSges(5,1) * t271 - rSges(5,2) * t270 + rSges(5,3) * t342;
t221 = t377 * t387 + (-t375 * t387 - qJD(4) * t322 + (qJ(4) * t319 + t322 * t389) * qJD(2)) * t313;
t360 = -t204 - t221;
t239 = rSges(5,1) * t287 - rSges(5,2) * t286 - rSges(5,3) * t380;
t359 = -t239 - t253;
t358 = t265 * t384 + t266 * t383;
t268 = pkin(2) * t299 - pkin(7) * t327;
t357 = t268 * t384 + t269 * t383;
t356 = qJD(2) * t313;
t348 = -t165 - t373;
t347 = -t221 - t372;
t346 = -t196 - t370;
t345 = -t197 - t369;
t344 = t164 * t380 - t221 * t327 + t253 * t290;
t343 = -t253 - t364;
t220 = rSges(4,1) * t283 + rSges(4,2) * t282 + rSges(4,3) * t342;
t297 = (pkin(2) * t322 + pkin(7) * t319) * t356;
t339 = (-t220 - t297) * t313;
t254 = rSges(4,1) * t303 + rSges(4,2) * t302 - rSges(4,3) * t380;
t304 = (pkin(2) * t319 - pkin(7) * t322) * t313;
t338 = (-t254 - t304) * t313;
t337 = t164 * t384 + t165 * t383 + t358;
t336 = t196 * t384 + t197 * t383 + t357;
t335 = (-t297 + t360) * t313;
t334 = (-t304 + t359) * t313;
t331 = (-t297 + t347) * t313;
t330 = (-t304 + t343) * t313;
t296 = (rSges(3,1) * t322 - rSges(3,2) * t319) * t356;
t295 = (Icges(3,1) * t322 - t386) * t356;
t294 = (-Icges(3,2) * t319 + t385) * t356;
t293 = (Icges(3,5) * t322 - Icges(3,6) * t319) * t356;
t288 = t315 * rSges(3,3) + (rSges(3,1) * t319 + rSges(3,2) * t322) * t313;
t285 = Icges(3,5) * t315 + (Icges(3,1) * t319 + t385) * t313;
t284 = Icges(3,6) * t315 + (Icges(3,2) * t322 + t386) * t313;
t264 = -rSges(3,1) * t291 - rSges(3,2) * t292;
t263 = rSges(3,1) * t289 - rSges(3,2) * t290;
t262 = -Icges(3,1) * t291 - Icges(3,4) * t292;
t261 = Icges(3,1) * t289 - Icges(3,4) * t290;
t260 = -Icges(3,4) * t291 - Icges(3,2) * t292;
t259 = Icges(3,4) * t289 - Icges(3,2) * t290;
t258 = -Icges(3,5) * t291 - Icges(3,6) * t292;
t257 = Icges(3,5) * t289 - Icges(3,6) * t290;
t249 = rSges(3,1) * t301 - rSges(3,2) * t300 + rSges(3,3) * t384;
t248 = rSges(3,1) * t299 + rSges(3,2) * t327 - rSges(3,3) * t383;
t246 = Icges(3,1) * t301 - Icges(3,4) * t300 + Icges(3,5) * t384;
t245 = Icges(3,1) * t299 + Icges(3,4) * t327 - Icges(3,5) * t383;
t244 = Icges(3,4) * t301 - Icges(3,2) * t300 + Icges(3,6) * t384;
t243 = Icges(3,4) * t299 + Icges(3,2) * t327 - Icges(3,6) * t383;
t213 = rSges(4,1) * t281 + rSges(4,2) * t280 + rSges(4,3) * t300;
t212 = -rSges(4,1) * t328 + rSges(4,2) * t278 - rSges(4,3) * t327;
t179 = t197 * t342;
t178 = (t263 * t312 + t264 * t314) * t313;
t177 = t300 * t196;
t173 = rSges(4,1) * t235 + rSges(4,2) * t234 + rSges(4,3) * t292;
t172 = rSges(4,1) * t233 + rSges(4,2) * t232 + rSges(4,3) * t290;
t156 = rSges(5,1) * t224 - rSges(5,2) * t223 + rSges(5,3) * t290;
t153 = -t213 * t380 - t254 * t300;
t152 = t212 * t380 - t254 * t327;
t136 = -t250 * t380 + t251 * t302 + t252 * t303;
t135 = t212 * t300 + t213 * t327;
t134 = (-t212 - t268) * t315 + t314 * t338;
t133 = t213 * t315 + t312 * t338 + t267;
t132 = -t236 * t380 - t237 * t286 + t238 * t287;
t123 = (t212 * t312 + t213 * t314) * t313 + t357;
t120 = (-t172 - t265) * t315 + t314 * t339;
t119 = t173 * t315 + t312 * t339 + t242;
t118 = t145 * t286 + t186 * t324;
t117 = -t144 * t286 - t186 * t323;
t114 = t300 * t359 + t361 * t380;
t113 = t194 * t380 - t239 * t327 + t365;
t108 = (-t268 + t362) * t315 + t314 * t334;
t107 = t195 * t315 + t312 * t334 + t363;
t94 = t183 * t286 + t184 * t276 - t185 * t329;
t93 = (t172 * t312 + t173 * t314) * t313 + t358;
t92 = -t144 * t324 + t145 * t323;
t91 = -t300 * t220 - t292 * t254 + (-t173 * t322 + t213 * t355) * t313;
t90 = -t327 * t220 + t290 * t254 + (t172 * t322 - t212 * t355) * t313;
t89 = t194 * t300 - t327 * t361 + t177;
t86 = t302 * t218 + t303 * t219 + t282 * t251 + t283 * t252 + (-t217 * t322 + t250 * t355) * t313;
t85 = (t194 * t312 + t195 * t314) * t313 + t336;
t84 = t172 * t300 + t173 * t327 + t212 * t292 - t213 * t290;
t83 = (-t156 + t366) * t315 + t314 * t335;
t82 = t157 * t315 + t312 * t335 + t367;
t81 = -t286 * t202 + t287 * t203 - t270 * t237 + t271 * t238 + (-t201 * t322 + t236 * t355) * t313;
t78 = t139 * t286 + t141 * t276 - t143 * t329;
t77 = t138 * t286 + t140 * t276 - t142 * t329;
t76 = t300 * t343 + t345 * t380;
t75 = -t327 * t364 + t370 * t380 + t365;
t74 = (-t268 + t346) * t315 + t314 * t330;
t73 = t312 * t330 + t315 * t369 + t363;
t66 = (t156 * t312 + t157 * t314) * t313 + t337;
t65 = t300 * t370 - t327 * t345 + t177;
t64 = t302 * t169 + t303 * t171 + t282 * t209 + t283 * t211 + (-t167 * t322 + t207 * t355) * t313;
t63 = t302 * t168 + t303 * t170 + t282 * t208 + t283 * t210 + (-t166 * t322 + t206 * t355) * t313;
t62 = (t312 * t370 + t314 * t369) * t313 + t336;
t61 = t179 + t360 * t300 + t359 * t292 + (t195 * t355 + t322 * t368) * t313;
t60 = -t327 * t204 + t290 * t239 + (t156 * t322 + t355 * t362) * t313 + t344;
t59 = -t286 * t149 + t287 * t151 - t270 * t191 + t271 * t193 + (-t147 * t322 + t189 * t355) * t313;
t58 = -t286 * t148 + t287 * t150 - t270 * t190 + t271 * t192 + (-t146 * t322 + t188 * t355) * t313;
t49 = t102 * t286 + t129 * t324 + t145 * t270 - t186 * t225;
t48 = -t101 * t286 - t129 * t323 - t144 * t270 + t186 * t223;
t47 = t156 * t300 + t194 * t292 + t290 * t361 - t327 * t368 + t371;
t46 = (t366 - t374) * t315 + t314 * t331;
t45 = t312 * t331 + t315 * t373 + t367;
t44 = t126 * t286 + t127 * t276 - t128 * t329 + t183 * t270 + t184 * t199 + t185 * t200;
t43 = -t101 * t324 + t102 * t323 + t144 * t225 - t145 * t223;
t42 = (t312 * t374 + t314 * t373) * t313 + t337;
t39 = t315 * t94 + (t312 * t78 - t314 * t77) * t313;
t38 = t300 * t78 - t327 * t77 - t380 * t94;
t37 = t286 * t94 - t323 * t77 - t324 * t78;
t36 = t179 + t347 * t300 + t343 * t292 + (t322 * t348 + t355 * t369) * t313;
t35 = -t372 * t327 + t364 * t290 + (t322 * t374 + t346 * t355) * t313 + t344;
t34 = t315 * t88 + (t312 * t72 - t314 * t71) * t313;
t33 = t315 * t87 + (t312 * t70 - t314 * t69) * t313;
t32 = t300 * t72 - t327 * t71 - t380 * t88;
t31 = t300 * t70 - t327 * t69 - t380 * t87;
t30 = t286 * t88 - t323 * t71 - t324 * t72;
t29 = t286 * t87 - t323 * t69 - t324 * t70;
t28 = -t100 * t329 + t139 * t270 + t141 * t199 + t143 * t200 + t276 * t98 + t286 * t96;
t27 = t138 * t270 + t140 * t199 + t142 * t200 + t276 * t97 + t286 * t95 - t329 * t99;
t26 = t315 * t86 + (t312 * t64 - t314 * t63) * t313;
t21 = t290 * t345 + t292 * t370 + t300 * t374 - t327 * t348 + t371;
t20 = t315 * t81 + (t312 * t59 - t314 * t58) * t313;
t19 = t315 * t80 + (t312 * t57 - t314 * t56) * t313;
t18 = t315 * t79 + (t312 * t55 - t314 * t54) * t313;
t17 = t315 * t68 + (t312 * t53 - t314 * t52) * t313;
t16 = t315 * t67 + (t312 * t51 - t314 * t50) * t313;
t15 = t121 * t290 + t122 * t292 - t63 * t327 + t64 * t300 + (t136 * t355 - t322 * t86) * t313;
t14 = t115 * t290 + t116 * t292 - t58 * t327 + t59 * t300 + (t132 * t355 - t322 * t81) * t313;
t9 = t315 * t44 + (-t27 * t314 + t28 * t312) * t313;
t8 = t315 * t41 + (-t24 * t314 + t25 * t312) * t313;
t7 = t315 * t40 + (-t22 * t314 + t23 * t312) * t313;
t6 = -t27 * t327 + t28 * t300 + t77 * t290 + t78 * t292 + (-t322 * t44 + t355 * t94) * t313;
t5 = t223 * t77 + t225 * t78 - t27 * t323 + t270 * t94 - t28 * t324 + t286 * t44;
t2 = t223 * t71 + t225 * t72 - t24 * t323 - t25 * t324 + t270 * t88 + t286 * t41;
t1 = -t22 * t323 + t223 * t69 + t225 * t70 - t23 * t324 + t270 * t87 + t286 * t40;
t10 = [0; m(3) * t178 + m(4) * t93 + m(5) * t66 + m(6) * t42; -t7 * t383 + t8 * t384 - t18 * t383 - t16 * t383 + t19 * t384 + t17 * t384 - ((-t244 * t290 + t246 * t289 - t258 * t383 + t260 * t327 + t262 * t299) * t384 - (-t243 * t290 + t245 * t289 - t257 * t383 + t259 * t327 + t261 * t299) * t383 + (-t284 * t290 + t285 * t289 - t293 * t383 + t294 * t327 + t295 * t299) * t315) * t383 + ((-t244 * t292 - t246 * t291 + t258 * t384 - t260 * t300 + t262 * t301) * t384 - (-t243 * t292 - t245 * t291 + t257 * t384 - t259 * t300 + t261 * t301) * t383 + (-t284 * t292 - t285 * t291 + t293 * t384 - t294 * t300 + t295 * t301) * t315) * t384 + (t42 * t62 + t45 * t73 + t46 * t74) * t400 + t315 * t9 + t315 * t26 + t315 * t20 + (t107 * t82 + t108 * t83 + t66 * t85) * t401 + (t119 * t133 + t120 * t134 + t123 * t93) * t402 + t315 * (t315 ^ 2 * t293 + (((t260 * t322 + t262 * t319) * t312 - (t259 * t322 + t261 * t319) * t314 + ((-t244 * t319 + t246 * t322) * t312 - (-t243 * t319 + t245 * t322) * t314) * qJD(2)) * t313 + (-t257 * t314 + t258 * t312 + t294 * t322 + t295 * t319 + (-t284 * t319 + t285 * t322) * qJD(2)) * t315) * t313) + 0.2e1 * m(3) * ((-t248 * t315 - t288 * t383) * (-t263 * t315 - t296 * t383) + (t249 * t315 - t288 * t384) * (t264 * t315 - t296 * t384) + (t248 * t312 + t249 * t314) * t313 * t178); m(4) * t84 + m(5) * t47 + m(6) * t21; (t19 / 0.2e1 + t17 / 0.2e1 + t8 / 0.2e1) * t300 - (t18 / 0.2e1 + t16 / 0.2e1 + t7 / 0.2e1) * t327 + m(4) * (t119 * t153 + t120 * t152 + t123 * t84 + t133 * t91 + t134 * t90 + t135 * t93) + m(5) * (t107 * t61 + t108 * t60 + t113 * t83 + t114 * t82 + t47 * t85 + t66 * t89) + m(6) * (t21 * t62 + t35 * t74 + t36 * t73 + t42 * t65 + t45 * t76 + t46 * t75) + (t15 / 0.2e1 + t14 / 0.2e1 + t6 / 0.2e1 + (t131 / 0.2e1 + t125 / 0.2e1) * t292 + (t130 / 0.2e1 + t124 / 0.2e1) * t290) * t315 + t33 * t393 + t34 * t392 + ((-t26 / 0.2e1 - t20 / 0.2e1 - t9 / 0.2e1) * t322 + (t39 / 0.2e1 + (t132 / 0.2e1 + t136 / 0.2e1) * t315) * t355 + (t342 * t405 + t403) * t391 + (t406 * t342 + t404) * t390 + (-t392 * t408 - t393 * t410) * t314 + (t392 * t407 + t393 * t409) * t312) * t313; (-t14 - t15 - t6) * t380 + t403 * t300 - t404 * t327 + (t407 * t300 - t408 * t327 - t380 * t411 + t32) * t292 + (t409 * t300 - t410 * t327 - t380 * t412 + t31) * t290 + (t38 + (-t132 - t136) * t380 + t405 * t300 - t406 * t327) * t342 + (t21 * t65 + t35 * t75 + t36 * t76) * t400 + (t113 * t60 + t114 * t61 + t47 * t89) * t401 + (t135 * t84 + t152 * t90 + t153 * t91) * t402; (m(5) + m(6)) * t342; m(6) * (t290 * t73 + t292 * t74 - t327 * t45 + t300 * t46 + (-t322 * t42 + t355 * t62) * t313) + m(5) * (t290 * t107 + t292 * t108 - t327 * t82 + t300 * t83 + (-t322 * t66 + t355 * t85) * t313); m(6) * (t290 * t76 + t292 * t75 - t327 * t36 + t300 * t35 + (-t21 * t322 + t355 * t65) * t313) + m(5) * (t292 * t113 + t290 * t114 - t327 * t61 + t300 * t60 + (-t322 * t47 + t355 * t89) * t313); 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-t313 ^ 2 * t319 * t354 - t290 * t327 + t292 * t300); m(6) * t43; m(6) * (t117 * t46 + t118 * t45 + t42 * t92 + t43 * t62 + t48 * t74 + t49 * t73) + t315 * t5 / 0.2e1 + t33 * t399 + t7 * t396 + t34 * t398 + t8 * t395 + t39 * t397 + t9 * t394 + (t1 * t390 + t2 * t391) * t313; m(6) * (t117 * t35 + t118 * t36 + t21 * t92 + t43 * t65 + t48 * t75 + t49 * t76) + t31 * t399 + t3 * t396 + t30 * t392 + t300 * t2 / 0.2e1 + t32 * t398 + t4 * t395 + t29 * t393 - t327 * t1 / 0.2e1 + t38 * t397 + t6 * t394 + (t37 * t355 / 0.2e1 - t322 * t5 / 0.2e1) * t313; m(6) * (t117 * t292 + t118 * t290 - t49 * t327 + t48 * t300 + (-t322 * t43 + t355 * t92) * t313); (t117 * t48 + t118 * t49 + t43 * t92) * t400 + t225 * t30 - t324 * t2 + t223 * t29 - t323 * t1 + t270 * t37 + t286 * t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
