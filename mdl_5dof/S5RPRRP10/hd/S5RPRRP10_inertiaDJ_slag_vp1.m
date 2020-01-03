% Calculate time derivative of joint inertia matrix for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP10_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP10_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:56
% EndTime: 2019-12-31 18:51:14
% DurationCPUTime: 10.61s
% Computational Cost: add. (18523->702), mult. (26625->1000), div. (0->0), fcn. (25328->8), ass. (0->343)
t233 = pkin(8) + qJ(3);
t228 = sin(t233);
t229 = cos(t233);
t240 = sin(qJ(4));
t242 = cos(qJ(4));
t276 = Icges(6,5) * t242 - Icges(6,6) * t240;
t167 = -Icges(6,3) * t229 + t228 * t276;
t277 = Icges(5,5) * t242 - Icges(5,6) * t240;
t168 = -Icges(5,3) * t229 + t228 * t277;
t438 = t167 + t168;
t383 = Icges(6,4) * t242;
t279 = -Icges(6,2) * t240 + t383;
t169 = -Icges(6,6) * t229 + t228 * t279;
t385 = Icges(5,4) * t242;
t280 = -Icges(5,2) * t240 + t385;
t170 = -Icges(5,6) * t229 + t228 * t280;
t433 = t170 + t169;
t384 = Icges(6,4) * t240;
t283 = Icges(6,1) * t242 - t384;
t171 = -Icges(6,5) * t229 + t228 * t283;
t386 = Icges(5,4) * t240;
t284 = Icges(5,1) * t242 - t386;
t172 = -Icges(5,5) * t229 + t228 * t284;
t432 = -t172 - t171;
t238 = -qJ(5) - pkin(7);
t437 = rSges(6,3) - t238;
t436 = -qJD(1) * t228 / 0.2e1;
t435 = t438 * t229 + (t433 * t240 + t432 * t242) * t228;
t337 = qJD(4) * t228;
t116 = (-Icges(6,5) * t240 - Icges(6,6) * t242) * t337 + (Icges(6,3) * t228 + t229 * t276) * qJD(3);
t117 = (-Icges(5,5) * t240 - Icges(5,6) * t242) * t337 + (Icges(5,3) * t228 + t229 * t277) * qJD(3);
t434 = -t117 - t116;
t118 = (-Icges(6,2) * t242 - t384) * t337 + (Icges(6,6) * t228 + t229 * t279) * qJD(3);
t119 = (-Icges(5,2) * t242 - t386) * t337 + (Icges(5,6) * t228 + t229 * t280) * qJD(3);
t431 = (-t118 - t119) * t240;
t241 = sin(qJ(1));
t243 = cos(qJ(1));
t358 = t241 * t242;
t360 = t240 * t243;
t195 = -t229 * t360 + t358;
t357 = t242 * t243;
t359 = t241 * t240;
t196 = t229 * t357 + t359;
t365 = t228 * t243;
t137 = Icges(6,5) * t196 + Icges(6,6) * t195 + Icges(6,3) * t365;
t141 = Icges(6,4) * t196 + Icges(6,2) * t195 + Icges(6,6) * t365;
t145 = Icges(6,1) * t196 + Icges(6,4) * t195 + Icges(6,5) * t365;
t271 = -t141 * t240 + t145 * t242;
t62 = -t137 * t229 + t228 * t271;
t139 = Icges(5,5) * t196 + Icges(5,6) * t195 + Icges(5,3) * t365;
t143 = Icges(5,4) * t196 + Icges(5,2) * t195 + Icges(5,6) * t365;
t147 = Icges(5,1) * t196 + Icges(5,4) * t195 + Icges(5,5) * t365;
t269 = -t143 * t240 + t147 * t242;
t64 = -t139 * t229 + t228 * t269;
t396 = t62 + t64;
t193 = -t229 * t359 - t357;
t194 = t229 * t358 - t360;
t367 = t228 * t241;
t136 = Icges(6,5) * t194 + Icges(6,6) * t193 + Icges(6,3) * t367;
t140 = Icges(6,4) * t194 + Icges(6,2) * t193 + Icges(6,6) * t367;
t144 = Icges(6,1) * t194 + Icges(6,4) * t193 + Icges(6,5) * t367;
t272 = -t140 * t240 + t144 * t242;
t61 = -t136 * t229 + t228 * t272;
t138 = Icges(5,5) * t194 + Icges(5,6) * t193 + Icges(5,3) * t367;
t142 = Icges(5,4) * t194 + Icges(5,2) * t193 + Icges(5,6) * t367;
t146 = Icges(5,1) * t194 + Icges(5,4) * t193 + Icges(5,5) * t367;
t270 = -t142 * t240 + t146 * t242;
t63 = -t138 * t229 + t228 * t270;
t397 = t61 + t63;
t430 = t241 * t397 + t243 * t396;
t429 = t241 * t396 - t243 * t397;
t304 = qJD(4) * t229 - qJD(1);
t264 = t242 * t304;
t303 = qJD(1) * t229 - qJD(4);
t340 = qJD(3) * t241;
t316 = t228 * t340;
t124 = -t241 * t264 + (-t243 * t303 + t316) * t240;
t263 = t304 * t240;
t339 = qJD(3) * t242;
t125 = t303 * t357 + (-t228 * t339 - t263) * t241;
t312 = t229 * t340;
t344 = qJD(1) * t243;
t250 = t228 * t344 + t312;
t69 = Icges(6,5) * t125 + Icges(6,6) * t124 + Icges(6,3) * t250;
t73 = Icges(6,4) * t125 + Icges(6,2) * t124 + Icges(6,6) * t250;
t77 = Icges(6,1) * t125 + Icges(6,4) * t124 + Icges(6,5) * t250;
t19 = (qJD(3) * t272 - t69) * t229 + (qJD(3) * t136 - t240 * t73 + t242 * t77 + (-t140 * t242 - t144 * t240) * qJD(4)) * t228;
t71 = Icges(5,5) * t125 + Icges(5,6) * t124 + Icges(5,3) * t250;
t75 = Icges(5,4) * t125 + Icges(5,2) * t124 + Icges(5,6) * t250;
t79 = Icges(5,1) * t125 + Icges(5,4) * t124 + Icges(5,5) * t250;
t21 = (qJD(3) * t270 - t71) * t229 + (qJD(3) * t138 - t240 * t75 + t242 * t79 + (-t142 * t242 - t146 * t240) * qJD(4)) * t228;
t428 = -t19 - t21;
t338 = qJD(3) * t243;
t315 = t228 * t338;
t416 = t241 * t303 + t315;
t122 = t240 * t416 - t243 * t264;
t123 = -t242 * t416 - t243 * t263;
t311 = t229 * t338;
t345 = qJD(1) * t241;
t319 = t228 * t345;
t249 = t311 - t319;
t68 = Icges(6,5) * t123 + Icges(6,6) * t122 + Icges(6,3) * t249;
t72 = Icges(6,4) * t123 + Icges(6,2) * t122 + Icges(6,6) * t249;
t76 = Icges(6,1) * t123 + Icges(6,4) * t122 + Icges(6,5) * t249;
t20 = (qJD(3) * t271 - t68) * t229 + (qJD(3) * t137 - t240 * t72 + t242 * t76 + (-t141 * t242 - t145 * t240) * qJD(4)) * t228;
t70 = Icges(5,5) * t123 + Icges(5,6) * t122 + Icges(5,3) * t249;
t74 = Icges(5,4) * t123 + Icges(5,2) * t122 + Icges(5,6) * t249;
t78 = Icges(5,1) * t123 + Icges(5,4) * t122 + Icges(5,5) * t249;
t22 = (qJD(3) * t269 - t70) * t229 + (qJD(3) * t139 - t240 * t74 + t242 * t78 + (-t143 * t242 - t147 * t240) * qJD(4)) * t228;
t427 = t20 + t22;
t53 = t138 * t367 + t142 * t193 + t146 * t194;
t54 = t139 * t367 + t143 * t193 + t147 * t194;
t291 = t241 * t53 + t243 * t54;
t51 = t136 * t367 + t140 * t193 + t144 * t194;
t52 = t137 * t367 + t141 * t193 + t145 * t194;
t292 = t241 * t51 + t243 * t52;
t84 = t167 * t367 + t169 * t193 + t171 * t194;
t85 = t168 * t367 + t170 * t193 + t172 * t194;
t400 = (-t84 - t85) * t229 + (t291 + t292) * t228;
t57 = t138 * t365 + t195 * t142 + t196 * t146;
t58 = t139 * t365 + t195 * t143 + t196 * t147;
t289 = t241 * t57 + t243 * t58;
t55 = t136 * t365 + t195 * t140 + t196 * t144;
t56 = t137 * t365 + t195 * t141 + t196 * t145;
t290 = t241 * t55 + t243 * t56;
t86 = t167 * t365 + t195 * t169 + t196 * t171;
t87 = t168 * t365 + t195 * t170 + t196 * t172;
t426 = (-t86 - t87) * t229 + (t289 + t290) * t228;
t401 = pkin(7) + t238;
t425 = t229 * t401;
t388 = Icges(4,4) * t228;
t286 = Icges(4,1) * t229 - t388;
t180 = Icges(4,5) * t241 + t243 * t286;
t369 = t180 * t229;
t387 = Icges(4,4) * t229;
t282 = -Icges(4,2) * t228 + t387;
t178 = Icges(4,6) * t241 + t243 * t282;
t374 = t178 * t228;
t265 = -t369 + t374;
t424 = t241 * t265;
t179 = -Icges(4,5) * t243 + t241 * t286;
t371 = t179 * t229;
t177 = -Icges(4,6) * t243 + t241 * t282;
t376 = t177 * t228;
t266 = -t371 + t376;
t423 = t243 * t266;
t226 = pkin(4) * t242 + pkin(3);
t364 = t229 * t243;
t422 = t196 * rSges(6,1) + t195 * rSges(6,2) + pkin(4) * t359 + t226 * t364 + t437 * t365;
t335 = qJD(4) * t242;
t330 = pkin(4) * t335;
t333 = pkin(4) * t360;
t334 = qJD(5) * t228;
t421 = t123 * rSges(6,1) + t122 * rSges(6,2) + rSges(6,3) * t311 + qJD(1) * t333 + t238 * t319 + t241 * t330 + t243 * t334;
t120 = (-Icges(6,1) * t240 - t383) * t337 + (Icges(6,5) * t228 + t229 * t283) * qJD(3);
t121 = (-Icges(5,1) * t240 - t385) * t337 + (Icges(5,5) * t228 + t229 * t284) * qJD(3);
t343 = qJD(3) * t228;
t419 = (t120 + t121) * t228 * t242 + t438 * t343 - t432 * t229 * t339;
t232 = t241 * rSges(4,3);
t418 = -rSges(4,2) * t365 + t232;
t239 = -pkin(6) - qJ(2);
t237 = cos(pkin(8));
t225 = pkin(2) * t237 + pkin(1);
t392 = rSges(4,1) * t229;
t299 = -rSges(4,2) * t228 + t392;
t256 = -t225 - t299;
t158 = (rSges(4,3) - t239) * t243 + t256 * t241;
t182 = rSges(4,1) * t364 + t418;
t307 = t243 * t225 - t241 * t239;
t159 = t182 + t307;
t417 = t158 * t243 + t159 * t241;
t278 = Icges(4,5) * t229 - Icges(4,6) * t228;
t175 = -Icges(4,3) * t243 + t241 * t278;
t262 = rSges(3,1) * t237 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t390 = rSges(3,3) + qJ(2);
t166 = t241 * t390 + t243 * t262;
t220 = pkin(3) * t364;
t191 = pkin(7) * t365 + t220;
t353 = -t191 + t422;
t402 = pkin(3) - t226;
t247 = -t228 * t401 - t229 * t402;
t294 = -t194 * rSges(6,1) - t193 * rSges(6,2);
t354 = rSges(6,3) * t367 + t241 * t247 - t294 - t333;
t415 = -t241 * t354 - t243 * t353;
t414 = 2 * m(4);
t413 = 2 * m(5);
t412 = 2 * m(6);
t234 = t241 ^ 2;
t235 = t243 ^ 2;
t411 = -t229 / 0.2e1;
t408 = -rSges(5,3) - pkin(7);
t205 = rSges(4,1) * t228 + rSges(4,2) * t229;
t407 = m(4) * t205;
t406 = pkin(3) * t229;
t405 = pkin(4) * t240;
t341 = qJD(3) * t240;
t398 = t419 + (-t433 * t341 + t434) * t229 + ((t432 * t240 - t433 * t242) * qJD(4) + t431) * t228;
t213 = pkin(7) * t311;
t336 = qJD(4) * t240;
t331 = pkin(4) * t336;
t395 = -t213 + (pkin(7) * t345 + t338 * t402) * t228 + ((-qJD(3) * t238 - t331) * t243 + t402 * t345) * t229 - rSges(6,3) * t319 + t421;
t212 = pkin(3) * t316;
t295 = t125 * rSges(6,1) + t124 * rSges(6,2);
t368 = t226 * t228;
t394 = t212 + (qJD(1) * t247 - t330) * t243 + (t334 - pkin(4) * t263 + (-t368 - t425) * qJD(3)) * t241 + rSges(6,3) * t250 + t295;
t393 = t435 * t343;
t391 = rSges(6,3) * t228;
t297 = -t194 * rSges(5,1) - t193 * rSges(5,2);
t149 = rSges(5,3) * t367 - t297;
t379 = t149 * t243;
t375 = t177 * t229;
t373 = t178 * t229;
t372 = t179 * t228;
t370 = t180 * t228;
t363 = t239 * t243;
t293 = rSges(6,1) * t242 - rSges(6,2) * t240;
t314 = t228 * t336;
t356 = -pkin(4) * t314 - qJD(5) * t229 + (-rSges(6,1) * t240 - rSges(6,2) * t242) * t337 + (t229 * t293 + t247 + t391) * qJD(3);
t296 = rSges(5,1) * t242 - rSges(5,2) * t240;
t127 = (-rSges(5,1) * t240 - rSges(5,2) * t242) * t337 + (rSges(5,3) * t228 + t229 * t296) * qJD(3);
t301 = pkin(7) * t228 + t406;
t201 = t301 * qJD(3);
t355 = -t127 - t201;
t151 = t196 * rSges(5,1) + t195 * rSges(5,2) + rSges(5,3) * t365;
t352 = -t151 - t191;
t351 = -rSges(6,3) * t229 + t425 + (t293 - t402) * t228;
t174 = -rSges(5,3) * t229 + t228 * t296;
t207 = pkin(3) * t228 - pkin(7) * t229;
t350 = -t174 - t207;
t190 = t301 * t241;
t349 = t241 * t190 + t243 * t191;
t231 = qJD(2) * t243;
t348 = t239 * t345 + t231;
t347 = t234 + t235;
t176 = Icges(4,3) * t241 + t243 * t278;
t346 = qJD(1) * t176;
t342 = qJD(3) * t229;
t35 = t52 * t241 - t243 * t51;
t36 = t54 * t241 - t243 * t53;
t329 = t36 / 0.2e1 + t35 / 0.2e1;
t37 = t56 * t241 - t243 * t55;
t38 = t58 * t241 - t243 * t57;
t328 = t37 / 0.2e1 + t38 / 0.2e1;
t325 = -t201 - t356;
t323 = t123 * rSges(5,1) + t122 * rSges(5,2) + rSges(5,3) * t311;
t322 = t241 * (pkin(7) * t250 + qJD(1) * t220 - t212) + t243 * (-pkin(7) * t319 + t213 + (-t229 * t345 - t315) * pkin(3)) + t190 * t344;
t321 = -t207 - t351;
t320 = t174 * t345;
t310 = t241 * t351;
t309 = t243 * t351;
t308 = t354 * t243;
t153 = t350 * t243;
t306 = -t226 * t229 - t225;
t305 = qJD(1) * t351;
t104 = t321 * t243;
t300 = t241 * t305;
t298 = t125 * rSges(5,1) + t124 * rSges(5,2);
t59 = t228 * t310 + t229 * t354;
t60 = -t228 * t309 - t229 * t353;
t288 = t241 * t60 + t243 * t59;
t248 = -t228 * t437 + t306;
t95 = (-t239 + t405) * t243 + t248 * t241 + t294;
t96 = t307 + t422;
t287 = t241 * t96 + t243 * t95;
t281 = Icges(4,2) * t229 + t388;
t103 = t321 * t241;
t275 = t103 * t241 + t104 * t243;
t268 = -t151 * t241 + t379;
t267 = -t241 * t149 - t151 * t243;
t33 = t116 * t367 + t193 * t118 + t194 * t120 + t124 * t169 + t125 * t171 + t167 * t250;
t34 = t117 * t367 + t193 * t119 + t194 * t121 + t124 * t170 + t125 * t172 + t168 * t250;
t260 = t21 / 0.2e1 + t19 / 0.2e1 + t33 / 0.2e1 + t34 / 0.2e1;
t31 = t116 * t365 + t195 * t118 + t196 * t120 + t122 * t169 + t123 * t171 + t167 * t249;
t32 = t117 * t365 + t195 * t119 + t196 * t121 + t122 * t170 + t123 * t172 + t168 * t249;
t259 = t22 / 0.2e1 + t20 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1;
t258 = t84 / 0.2e1 + t85 / 0.2e1 + t63 / 0.2e1 + t61 / 0.2e1;
t257 = t86 / 0.2e1 + t87 / 0.2e1 + t64 / 0.2e1 + t62 / 0.2e1;
t255 = qJD(3) * t205;
t253 = qJD(3) * t281;
t252 = qJD(3) * (-Icges(4,5) * t228 - Icges(4,6) * t229);
t251 = t228 * t408 - t225 - t406;
t246 = -t241 * t353 + t308;
t245 = rSges(4,2) * t319 + rSges(4,3) * t344 - t243 * t255;
t244 = t241 * t251 - t363;
t165 = -t241 * t262 + t243 * t390;
t230 = qJD(2) * t241;
t200 = t299 * qJD(3);
t192 = t207 * t345;
t181 = -rSges(4,3) * t243 + t241 * t299;
t157 = -qJD(1) * t166 + t231;
t156 = qJD(1) * t165 + t230;
t152 = t350 * t241;
t131 = t241 * t252 + t346;
t130 = -qJD(1) * t175 + t243 * t252;
t108 = t205 * t340 + (t243 * t256 - t232) * qJD(1) + t348;
t107 = t230 + (-t363 + (-t225 - t392) * t241) * qJD(1) + t245;
t106 = t307 - t352;
t105 = t244 + t297;
t102 = -t229 * t151 - t174 * t365;
t101 = t149 * t229 + t174 * t367;
t100 = t241 * t176 - t243 * t265;
t99 = t241 * t175 - t423;
t98 = -t176 * t243 - t424;
t97 = -t175 * t243 - t241 * t266;
t90 = t268 * t228;
t89 = qJD(1) * t153 + t241 * t355;
t88 = t243 * t355 + t192 + t320;
t83 = rSges(5,3) * t250 + t298;
t81 = -rSges(5,3) * t319 + t323;
t65 = -t267 + t349;
t50 = t251 * t344 + t312 * t408 + t212 - t298 + t348;
t49 = -pkin(3) * t315 + qJD(1) * t244 + t213 + t230 + t323;
t48 = qJD(1) * t104 + t241 * t325;
t47 = t243 * t325 + t192 + t300;
t46 = t246 * t228;
t45 = (qJD(1) * t248 + t330) * t243 + (-t334 + t304 * t405 + (-t229 * t437 + t368) * qJD(3)) * t241 - t295 + t348;
t44 = t230 + (-t229 * t331 + (-t229 * t238 - t368) * qJD(3)) * t243 + (-t363 + (t306 - t391) * t241) * qJD(1) + t421;
t43 = t349 - t415;
t42 = (t174 * t340 + t83) * t229 + (-qJD(3) * t149 + t241 * t127 + t174 * t344) * t228;
t41 = (-t174 * t338 - t81) * t229 + (qJD(3) * t151 - t127 * t243 + t320) * t228;
t30 = t268 * t342 + (qJD(1) * t267 - t241 * t81 + t243 * t83) * t228;
t29 = t241 * t83 + t243 * t81 + (t241 * t352 + t379) * qJD(1) + t322;
t24 = (qJD(3) * t310 + t394) * t229 + (-qJD(3) * t354 + t241 * t356 + t243 * t305) * t228;
t23 = (-qJD(3) * t309 - t395) * t229 + (qJD(3) * t353 - t243 * t356 + t300) * t228;
t18 = t124 * t143 + t125 * t147 + t139 * t250 + t193 * t74 + t194 * t78 + t367 * t70;
t17 = t124 * t142 + t125 * t146 + t138 * t250 + t193 * t75 + t194 * t79 + t367 * t71;
t16 = t124 * t141 + t125 * t145 + t137 * t250 + t193 * t72 + t194 * t76 + t367 * t68;
t15 = t124 * t140 + t125 * t144 + t136 * t250 + t193 * t73 + t194 * t77 + t367 * t69;
t14 = t122 * t143 + t123 * t147 + t139 * t249 + t195 * t74 + t196 * t78 + t365 * t70;
t13 = t122 * t142 + t123 * t146 + t138 * t249 + t195 * t75 + t196 * t79 + t365 * t71;
t12 = t122 * t141 + t123 * t145 + t137 * t249 + t195 * t72 + t196 * t76 + t365 * t68;
t11 = t122 * t140 + t123 * t144 + t136 * t249 + t195 * t73 + t196 * t77 + t365 * t69;
t10 = t395 * t243 + t394 * t241 + (t308 + (-t191 - t353) * t241) * qJD(1) + t322;
t9 = t246 * t342 + (qJD(1) * t415 - t395 * t241 + t394 * t243) * t228;
t8 = qJD(1) * t291 - t17 * t243 + t18 * t241;
t7 = qJD(1) * t292 - t15 * t243 + t16 * t241;
t6 = qJD(1) * t289 - t13 * t243 + t14 * t241;
t5 = qJD(1) * t290 - t11 * t243 + t12 * t241;
t4 = (qJD(3) * t291 - t34) * t229 + (-qJD(1) * t36 + qJD(3) * t85 + t17 * t241 + t18 * t243) * t228;
t3 = (qJD(3) * t292 - t33) * t229 + (-qJD(1) * t35 + qJD(3) * t84 + t15 * t241 + t16 * t243) * t228;
t2 = (qJD(3) * t289 - t32) * t229 + (-qJD(1) * t38 + qJD(3) * t87 + t13 * t241 + t14 * t243) * t228;
t1 = (qJD(3) * t290 - t31) * t229 + (-qJD(1) * t37 + qJD(3) * t86 + t11 * t241 + t12 * t243) * t228;
t25 = [0.2e1 * m(3) * (t156 * t166 + t157 * t165) + (t107 * t159 + t108 * t158) * t414 + (t44 * t96 + t45 * t95) * t412 + (t105 * t50 + t106 * t49) * t413 + (-t281 + t286) * t343 + (Icges(4,1) * t228 + t282 + t387) * t342 + t432 * t314 + t434 * t229 + t431 * t228 + t419 + t433 * (-t228 * t335 - t229 * t341); m(6) * (qJD(1) * t287 + t241 * t45 - t243 * t44) + m(5) * (t241 * t50 - t243 * t49 + (t105 * t243 + t106 * t241) * qJD(1)) + m(4) * (qJD(1) * t417 - t107 * t243 + t241 * t108) + m(3) * (-t156 * t243 + t241 * t157 + (t165 * t243 + t166 * t241) * qJD(1)); 0; ((qJD(1) * t178 - t241 * t253) * t411 + t180 * t436 + (t376 / 0.2e1 - t371 / 0.2e1) * qJD(3) - t260) * t243 + ((-qJD(1) * t177 - t243 * t253) * t229 / 0.2e1 + t179 * t436 + (-t374 / 0.2e1 + t369 / 0.2e1) * qJD(3) + t259) * t241 + m(4) * ((-t107 * t241 - t108 * t243) * t205 - t417 * t200) + m(5) * (t105 * t88 + t106 * t89 + t152 * t49 + t153 * t50) + m(6) * (t103 * t44 + t104 * t45 + t47 * t95 + t48 * t96) + (t235 / 0.2e1 + t234 / 0.2e1) * t278 * qJD(3) + ((-t159 * t407 + t373 / 0.2e1 + t370 / 0.2e1 + t257) * t243 + (t158 * t407 + t375 / 0.2e1 + t372 / 0.2e1 + t258) * t241) * qJD(1); m(5) * (t88 * t241 - t243 * t89 + (t152 * t241 + t153 * t243) * qJD(1)) + m(6) * (qJD(1) * t275 + t47 * t241 - t243 * t48); (t10 * t43 + t103 * t48 + t104 * t47) * t412 - t243 * t8 - t243 * t7 + t241 * t6 + (t152 * t89 + t153 * t88 + t29 * t65) * t413 + t241 * t5 + ((t241 * t181 + t182 * t243) * ((qJD(1) * t181 + t245) * t243 + (-t241 * t255 + (-t182 + t418) * qJD(1)) * t241) + t347 * t205 * t200) * t414 + t241 * ((t241 * t130 + (t99 + t424) * qJD(1)) * t241 + (t100 * qJD(1) + (t177 * t342 + t179 * t343) * t243 + (-t131 + (-t370 - t373) * qJD(3) + (t176 - t266) * qJD(1)) * t241) * t243) - t243 * ((t131 * t243 + (t98 + t423) * qJD(1)) * t243 + (t97 * qJD(1) + (-t178 * t342 - t180 * t343 + t346) * t241 + (-t130 + (t372 + t375) * qJD(3) - t265 * qJD(1)) * t243) * t241) + (t98 * t241 - t243 * t97 + t35 + t36) * t345 + (t100 * t241 - t243 * t99 + t37 + t38) * t344; m(6) * (t23 * t96 + t24 * t95 + t44 * t60 + t45 * t59) + m(5) * (t101 * t50 + t102 * t49 + t105 * t42 + t106 * t41) + ((t241 * t258 + t243 * t257) * qJD(3) - t398) * t229 + (t259 * t243 + t260 * t241 + (-t241 * t257 + t243 * t258) * qJD(1)) * t228 - t393; m(5) * (t42 * t241 - t243 * t41 + (t101 * t243 + t102 * t241) * qJD(1)) + m(6) * (qJD(1) * t288 - t23 * t243 + t24 * t241); m(5) * (t101 * t88 + t102 * t89 + t152 * t41 + t153 * t42 + t29 * t90 + t30 * t65) + m(6) * (t10 * t46 + t103 * t23 + t104 * t24 + t43 * t9 + t47 * t59 + t48 * t60) + (-t4 / 0.2e1 - t3 / 0.2e1 + t328 * t342) * t243 + (t2 / 0.2e1 + t1 / 0.2e1 + t329 * t342) * t241 + ((-t241 * t328 + t243 * t329) * qJD(1) + (t7 + t8) * t241 / 0.2e1 + (t5 + t6) * t243 / 0.2e1 + t429 * qJD(3) / 0.2e1) * t228 + (t430 * qJD(1) + t427 * t241 + t428 * t243) * t411 + (t400 * t241 + t426 * t243) * qJD(1) / 0.2e1; (t23 * t60 + t24 * t59 + t46 * t9) * t412 + (t101 * t42 + t102 * t41 + t30 * t90) * t413 + (t398 * t229 + ((-t229 * t396 + t426) * t243 + (-t229 * t397 + t400) * t241) * qJD(3) + t393) * t229 + ((t1 + t2) * t243 + (t3 + t4) * t241 + (t428 * t241 - t427 * t243) * t229 + (t430 * t228 + t435 * t229) * qJD(3) + (t229 * t429 - t241 * t426 + t400 * t243) * qJD(1)) * t228; m(6) * (t287 * t342 + (t241 * t44 + t243 * t45 + (-t241 * t95 + t243 * t96) * qJD(1)) * t228); 0; m(6) * ((qJD(3) * t275 - t10) * t229 + (qJD(3) * t43 + t241 * t48 + t243 * t47 + (t103 * t243 - t104 * t241) * qJD(1)) * t228); m(6) * ((qJD(3) * t288 - t9) * t229 + (qJD(3) * t46 + t23 * t241 + t24 * t243 + (-t241 * t59 + t243 * t60) * qJD(1)) * t228); (-0.1e1 + t347) * t228 * t342 * t412;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t25(1), t25(2), t25(4), t25(7), t25(11); t25(2), t25(3), t25(5), t25(8), t25(12); t25(4), t25(5), t25(6), t25(9), t25(13); t25(7), t25(8), t25(9), t25(10), t25(14); t25(11), t25(12), t25(13), t25(14), t25(15);];
Mq = res;
