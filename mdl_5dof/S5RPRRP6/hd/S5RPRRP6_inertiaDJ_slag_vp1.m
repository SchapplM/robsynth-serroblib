% Calculate time derivative of joint inertia matrix for
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:10
% EndTime: 2019-12-31 18:42:28
% DurationCPUTime: 10.00s
% Computational Cost: add. (19246->664), mult. (26145->926), div. (0->0), fcn. (24923->8), ass. (0->340)
t237 = cos(qJ(4));
t235 = sin(qJ(3));
t238 = cos(qJ(3));
t234 = sin(qJ(4));
t379 = Icges(6,4) * t234;
t282 = Icges(6,1) * t237 - t379;
t187 = -Icges(6,5) * t238 + t235 * t282;
t381 = Icges(5,4) * t234;
t283 = Icges(5,1) * t237 - t381;
t188 = -Icges(5,5) * t238 + t235 * t283;
t431 = -t188 - t187;
t440 = t431 * t237;
t275 = Icges(6,5) * t237 - Icges(6,6) * t234;
t183 = -Icges(6,3) * t238 + t235 * t275;
t276 = Icges(5,5) * t237 - Icges(5,6) * t234;
t184 = -Icges(5,3) * t238 + t235 * t276;
t439 = t183 + t184;
t378 = Icges(6,4) * t237;
t278 = -Icges(6,2) * t234 + t378;
t185 = -Icges(6,6) * t238 + t235 * t278;
t380 = Icges(5,4) * t237;
t279 = -Icges(5,2) * t234 + t380;
t186 = -Icges(5,6) * t238 + t235 * t279;
t437 = -t186 - t185;
t233 = -qJ(5) - pkin(7);
t438 = rSges(6,3) - t233;
t232 = qJ(1) + pkin(8);
t229 = sin(t232);
t230 = cos(t232);
t355 = t234 * t238;
t359 = t230 * t237;
t179 = -t229 * t355 - t359;
t353 = t237 * t238;
t361 = t230 * t234;
t180 = t229 * t353 - t361;
t363 = t229 * t235;
t120 = Icges(5,5) * t180 + Icges(5,6) * t179 + Icges(5,3) * t363;
t124 = Icges(5,4) * t180 + Icges(5,2) * t179 + Icges(5,6) * t363;
t128 = Icges(5,1) * t180 + Icges(5,4) * t179 + Icges(5,5) * t363;
t362 = t229 * t237;
t181 = -t230 * t355 + t362;
t364 = t229 * t234;
t182 = t230 * t353 + t364;
t360 = t230 * t235;
t53 = t120 * t360 + t124 * t181 + t128 * t182;
t121 = Icges(5,5) * t182 + Icges(5,6) * t181 + Icges(5,3) * t360;
t125 = Icges(5,4) * t182 + Icges(5,2) * t181 + Icges(5,6) * t360;
t129 = Icges(5,1) * t182 + Icges(5,4) * t181 + Icges(5,5) * t360;
t54 = t121 * t360 + t125 * t181 + t129 * t182;
t286 = t229 * t53 + t230 * t54;
t118 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t363;
t122 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t363;
t126 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t363;
t51 = t118 * t360 + t122 * t181 + t126 * t182;
t119 = Icges(6,5) * t182 + Icges(6,6) * t181 + Icges(6,3) * t360;
t123 = Icges(6,4) * t182 + Icges(6,2) * t181 + Icges(6,6) * t360;
t127 = Icges(6,1) * t182 + Icges(6,4) * t181 + Icges(6,5) * t360;
t52 = t119 * t360 + t123 * t181 + t127 * t182;
t287 = t229 * t51 + t230 * t52;
t436 = t286 + t287;
t49 = t120 * t363 + t124 * t179 + t128 * t180;
t50 = t121 * t363 + t125 * t179 + t129 * t180;
t288 = t229 * t49 + t230 * t50;
t47 = t118 * t363 + t122 * t179 + t126 * t180;
t48 = t119 * t363 + t123 * t179 + t127 * t180;
t289 = t229 * t47 + t230 * t48;
t435 = t288 + t289;
t339 = qJD(1) * t235;
t434 = -t339 / 0.2e1;
t433 = t439 * t238 + (-t437 * t234 + t440) * t235;
t334 = qJD(4) * t235;
t151 = (-Icges(6,2) * t237 - t379) * t334 + (Icges(6,6) * t235 + t238 * t278) * qJD(3);
t152 = (-Icges(5,2) * t237 - t381) * t334 + (Icges(5,6) * t235 + t238 * t279) * qJD(3);
t432 = -t152 - t151;
t335 = qJD(3) * t238;
t430 = t437 * t335;
t273 = -t123 * t234 + t127 * t237;
t61 = -t119 * t238 + t235 * t273;
t271 = -t125 * t234 + t129 * t237;
t63 = -t121 * t238 + t235 * t271;
t390 = t61 + t63;
t274 = -t122 * t234 + t126 * t237;
t60 = -t118 * t238 + t235 * t274;
t272 = -t124 * t234 + t128 * t237;
t62 = -t120 * t238 + t235 * t272;
t391 = t60 + t62;
t429 = t391 * t229 + t390 * t230;
t332 = qJD(4) * t238;
t302 = qJD(1) - t332;
t336 = qJD(3) * t235;
t245 = t234 * t336 + t237 * t302;
t338 = qJD(1) * t238;
t301 = -qJD(4) + t338;
t114 = t230 * t245 + t301 * t364;
t264 = t302 * t234;
t244 = -t237 * t336 + t264;
t115 = t230 * t244 - t301 * t362;
t315 = t230 * t335;
t317 = t229 * t339;
t248 = t315 - t317;
t116 = t229 * t245 - t301 * t361;
t117 = t229 * t244 + t301 * t359;
t312 = t229 * t335;
t249 = t230 * t339 + t312;
t70 = Icges(6,5) * t117 + Icges(6,6) * t116 + Icges(6,3) * t249;
t74 = Icges(6,4) * t117 + Icges(6,2) * t116 + Icges(6,6) * t249;
t78 = Icges(6,1) * t117 + Icges(6,4) * t116 + Icges(6,5) * t249;
t11 = t114 * t122 + t115 * t126 + t118 * t248 + t181 * t74 + t182 * t78 + t70 * t360;
t69 = Icges(6,5) * t115 + Icges(6,6) * t114 + Icges(6,3) * t248;
t73 = Icges(6,4) * t115 + Icges(6,2) * t114 + Icges(6,6) * t248;
t77 = Icges(6,1) * t115 + Icges(6,4) * t114 + Icges(6,5) * t248;
t12 = t114 * t123 + t115 * t127 + t119 * t248 + t181 * t73 + t182 * t77 + t69 * t360;
t72 = Icges(5,5) * t117 + Icges(5,6) * t116 + Icges(5,3) * t249;
t76 = Icges(5,4) * t117 + Icges(5,2) * t116 + Icges(5,6) * t249;
t80 = Icges(5,1) * t117 + Icges(5,4) * t116 + Icges(5,5) * t249;
t13 = t114 * t124 + t115 * t128 + t120 * t248 + t181 * t76 + t182 * t80 + t72 * t360;
t71 = Icges(5,5) * t115 + Icges(5,6) * t114 + Icges(5,3) * t248;
t75 = Icges(5,4) * t115 + Icges(5,2) * t114 + Icges(5,6) * t248;
t79 = Icges(5,1) * t115 + Icges(5,4) * t114 + Icges(5,5) * t248;
t14 = t114 * t125 + t115 * t129 + t121 * t248 + t181 * t75 + t182 * t79 + t71 * t360;
t428 = (-t11 - t13) * t230 + (t12 + t14) * t229 + t436 * qJD(1);
t15 = t116 * t122 + t117 * t126 + t118 * t249 + t179 * t74 + t180 * t78 + t70 * t363;
t16 = t116 * t123 + t117 * t127 + t119 * t249 + t179 * t73 + t180 * t77 + t69 * t363;
t17 = t116 * t124 + t117 * t128 + t120 * t249 + t179 * t76 + t180 * t80 + t72 * t363;
t18 = t116 * t125 + t117 * t129 + t121 * t249 + t179 * t75 + t180 * t79 + t71 * t363;
t427 = (-t15 - t17) * t230 + (t16 + t18) * t229 + t435 * qJD(1);
t19 = (t274 * qJD(3) - t70) * t238 + (qJD(3) * t118 - t234 * t74 + t237 * t78 + (-t122 * t237 - t126 * t234) * qJD(4)) * t235;
t21 = (t272 * qJD(3) - t72) * t238 + (qJD(3) * t120 - t234 * t76 + t237 * t80 + (-t124 * t237 - t128 * t234) * qJD(4)) * t235;
t426 = -t19 - t21;
t20 = (t273 * qJD(3) - t69) * t238 + (qJD(3) * t119 - t234 * t73 + t237 * t77 + (-t123 * t237 - t127 * t234) * qJD(4)) * t235;
t22 = (t271 * qJD(3) - t71) * t238 + (qJD(3) * t121 - t234 * t75 + t237 * t79 + (-t125 * t237 - t129 * t234) * qJD(4)) * t235;
t425 = t20 + t22;
t86 = t179 * t185 + t180 * t187 + t183 * t363;
t87 = t179 * t186 + t180 * t188 + t184 * t363;
t424 = (-t86 - t87) * t238 + t435 * t235;
t88 = t181 * t185 + t182 * t187 + t183 * t360;
t89 = t181 * t186 + t182 * t188 + t184 * t360;
t423 = (-t88 - t89) * t238 + t436 * t235;
t409 = 2 * m(4);
t387 = rSges(4,1) * t238;
t296 = -rSges(4,2) * t235 + t387;
t386 = rSges(4,3) * t230;
t167 = t229 * t296 - t386;
t358 = t230 * t238;
t416 = -rSges(4,2) * t360 + t229 * rSges(4,3);
t168 = rSges(4,1) * t358 + t416;
t212 = rSges(4,1) * t235 + rSges(4,2) * t238;
t256 = qJD(3) * t212;
t340 = qJD(1) * t230;
t243 = rSges(4,2) * t317 + rSges(4,3) * t340 - t230 * t256;
t59 = (qJD(1) * t167 + t243) * t230 + (-t229 * t256 + (-t168 + t416) * qJD(1)) * t229;
t422 = t409 * t59;
t383 = Icges(4,4) * t235;
t285 = Icges(4,1) * t238 - t383;
t166 = Icges(4,5) * t229 + t230 * t285;
t366 = t166 * t238;
t382 = Icges(4,4) * t238;
t281 = -Icges(4,2) * t235 + t382;
t164 = Icges(4,6) * t229 + t230 * t281;
t371 = t164 * t235;
t265 = -t366 + t371;
t421 = t230 * t265;
t394 = pkin(7) + t233;
t420 = t394 * t238;
t226 = pkin(4) * t237 + pkin(3);
t395 = pkin(3) - t226;
t419 = t395 * t235;
t418 = t395 * t238;
t417 = t182 * rSges(6,1) + t181 * rSges(6,2) + pkin(4) * t364 + t226 * t358 + t438 * t360;
t231 = cos(qJ(1)) * pkin(1);
t415 = t229 * pkin(6) + t231;
t149 = (-Icges(6,5) * t234 - Icges(6,6) * t237) * t334 + (Icges(6,3) * t235 + t238 * t275) * qJD(3);
t150 = (-Icges(5,5) * t234 - Icges(5,6) * t237) * t334 + (Icges(5,3) * t235 + t238 * t276) * qJD(3);
t153 = (-Icges(6,1) * t234 - t378) * t334 + (Icges(6,5) * t235 + t238 * t282) * qJD(3);
t154 = (-Icges(5,1) * t234 - t380) * t334 + (Icges(5,5) * t235 + t238 * t283) * qJD(3);
t333 = qJD(4) * t237;
t414 = t439 * t336 - t335 * t440 + (-t149 - t150) * t238 + ((t153 + t154) * t237 + t437 * t333) * t235;
t330 = pkin(4) * t361;
t413 = -rSges(6,1) * t180 - rSges(6,2) * t179 + t330;
t327 = pkin(4) * t333;
t331 = qJD(5) * t235;
t412 = t115 * rSges(6,1) + t114 * rSges(6,2) + rSges(6,3) * t315 + qJD(1) * t330 + t229 * t327 + t230 * t331 + t233 * t317;
t277 = Icges(4,5) * t238 - Icges(4,6) * t235;
t161 = -Icges(4,3) * t230 + t229 * t277;
t163 = -Icges(4,6) * t230 + t229 * t281;
t165 = -Icges(4,5) * t230 + t229 * t285;
t411 = t390 * t238 - t423;
t219 = pkin(3) * t358;
t192 = pkin(7) * t360 + t219;
t350 = -t192 + t417;
t247 = -t394 * t235 - t418;
t351 = rSges(6,3) * t363 + t229 * t247 - t413;
t410 = -t351 * t229 - t350 * t230;
t408 = 2 * m(5);
t407 = 2 * m(6);
t227 = t229 ^ 2;
t228 = t230 ^ 2;
t404 = -t238 / 0.2e1;
t403 = -rSges(5,3) - pkin(7);
t402 = m(4) * t212;
t401 = sin(qJ(1)) * pkin(1);
t400 = pkin(3) * t238;
t399 = pkin(4) * t234;
t398 = pkin(7) * t235;
t392 = t414 + (t430 + (t431 * qJD(4) + t432) * t235) * t234;
t206 = pkin(7) * t315;
t303 = t332 * t399;
t341 = qJD(1) * t229;
t357 = t233 * t238;
t389 = -t206 + (t398 + t418) * t341 + (-t303 + (-t357 + t419) * qJD(3)) * t230 - rSges(6,3) * t317 + t412;
t205 = t229 * pkin(3) * t336;
t292 = t117 * rSges(6,1) + t116 * rSges(6,2);
t365 = t226 * t235;
t388 = t205 + (qJD(1) * t247 - t327) * t230 + (t331 + pkin(4) * t264 + (-t365 - t420) * qJD(3)) * t229 + rSges(6,3) * t249 + t292;
t385 = rSges(6,3) * t235;
t294 = -rSges(5,1) * t180 - rSges(5,2) * t179;
t133 = rSges(5,3) * t363 - t294;
t374 = t133 * t230;
t373 = t163 * t235;
t372 = t163 * t238;
t370 = t164 * t238;
t369 = t165 * t235;
t368 = t165 * t238;
t367 = t166 * t235;
t352 = t433 * t336;
t135 = t182 * rSges(5,1) + t181 * rSges(5,2) + rSges(5,3) * t360;
t349 = -t135 - t192;
t290 = rSges(6,1) * t237 - rSges(6,2) * t234;
t311 = t234 * t334;
t348 = -pkin(4) * t311 - qJD(5) * t238 + (-rSges(6,1) * t234 - rSges(6,2) * t237) * t334 + (t238 * t290 + t247 + t385) * qJD(3);
t293 = rSges(5,1) * t237 - rSges(5,2) * t234;
t156 = (-rSges(5,1) * t234 - rSges(5,2) * t237) * t334 + (rSges(5,3) * t235 + t238 * t293) * qJD(3);
t299 = t398 + t400;
t200 = t299 * qJD(3);
t347 = -t156 - t200;
t191 = t299 * t229;
t346 = t229 * t191 + t230 * t192;
t345 = -rSges(6,3) * t238 + t235 * t290 - t419 + t420;
t190 = -rSges(5,3) * t238 + t235 * t293;
t216 = pkin(3) * t235 - pkin(7) * t238;
t344 = -t190 - t216;
t343 = t227 + t228;
t162 = Icges(4,3) * t229 + t230 * t277;
t342 = qJD(1) * t162;
t337 = qJD(3) * t229;
t329 = m(6) * t336;
t31 = t229 * t48 - t230 * t47;
t32 = t229 * t50 - t230 * t49;
t326 = t31 / 0.2e1 + t32 / 0.2e1;
t33 = t229 * t52 - t230 * t51;
t34 = t229 * t54 - t230 * t53;
t325 = t33 / 0.2e1 + t34 / 0.2e1;
t323 = t115 * rSges(5,1) + t114 * rSges(5,2) + rSges(5,3) * t315;
t316 = t230 * t336;
t322 = t229 * (pkin(7) * t249 + qJD(1) * t219 - t205) + t230 * (-pkin(7) * t317 + t206 + (-t229 * t338 - t316) * pkin(3)) + t191 * t340;
t321 = -t200 - t348;
t320 = -t216 - t345;
t319 = t230 * pkin(2) + t415;
t318 = t190 * t341;
t224 = t230 * pkin(6);
t309 = t224 - t401;
t308 = -t226 * t238 - pkin(2);
t307 = t351 * t230;
t306 = t345 * t230;
t305 = t345 * t229;
t158 = t344 * t230;
t304 = qJD(1) * t345;
t109 = t320 * t230;
t298 = t229 * t304;
t295 = t117 * rSges(5,1) + t116 * rSges(5,2);
t280 = Icges(4,2) * t238 + t383;
t270 = -t135 * t229 + t374;
t269 = -t133 * t229 - t135 * t230;
t266 = -t368 + t373;
t263 = -t391 * t238 + t424;
t262 = -pkin(2) - t296;
t37 = t116 * t185 + t117 * t187 + t149 * t363 + t151 * t179 + t153 * t180 + t183 * t249;
t38 = t116 * t186 + t117 * t188 + t150 * t363 + t152 * t179 + t154 * t180 + t184 * t249;
t260 = t19 / 0.2e1 + t21 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1;
t35 = t114 * t185 + t115 * t187 + t149 * t360 + t151 * t181 + t153 * t182 + t183 * t248;
t36 = t114 * t186 + t115 * t188 + t150 * t360 + t152 * t181 + t154 * t182 + t184 * t248;
t259 = t20 / 0.2e1 + t22 / 0.2e1 + t35 / 0.2e1 + t36 / 0.2e1;
t258 = t60 / 0.2e1 + t62 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1;
t257 = t61 / 0.2e1 + t63 / 0.2e1 + t88 / 0.2e1 + t89 / 0.2e1;
t255 = t403 * t235 - pkin(2) - t400;
t253 = qJD(3) * t280;
t252 = qJD(3) * (-Icges(4,5) * t235 - Icges(4,6) * t238);
t250 = -t235 * t438 + t308;
t246 = -t350 * t229 + t307;
t240 = t229 * t255 - t401;
t221 = pkin(6) * t340;
t199 = t296 * qJD(3);
t193 = t216 * t341;
t157 = t344 * t229;
t147 = t168 + t319;
t146 = t229 * t262 + t309 + t386;
t137 = t229 * t252 + t342;
t136 = -qJD(1) * t161 + t230 * t252;
t108 = t320 * t229;
t107 = t212 * t337 + (-t231 + (-rSges(4,3) - pkin(6)) * t229 + t262 * t230) * qJD(1);
t106 = t221 + (-t401 + (-pkin(2) - t387) * t229) * qJD(1) + t243;
t101 = -t135 * t238 - t190 * t360;
t100 = t133 * t238 + t190 * t363;
t99 = t319 - t349;
t98 = t224 + t240 + t294;
t97 = t162 * t229 - t421;
t96 = t161 * t229 - t230 * t266;
t95 = -t162 * t230 - t229 * t265;
t94 = -t161 * t230 - t229 * t266;
t93 = t319 + t417;
t92 = t229 * t250 + t309 + t413;
t91 = qJD(1) * t158 + t347 * t229;
t90 = t347 * t230 + t193 + t318;
t85 = t270 * t235;
t84 = rSges(5,3) * t249 + t295;
t82 = -rSges(5,3) * t317 + t323;
t66 = -t269 + t346;
t65 = -t235 * t306 - t350 * t238;
t64 = t235 * t305 + t351 * t238;
t58 = qJD(1) * t109 + t229 * t321;
t57 = t230 * t321 + t193 + t298;
t56 = t205 + t403 * t312 + (t230 * t255 - t415) * qJD(1) - t295;
t55 = -pkin(3) * t316 + qJD(1) * t240 + t206 + t221 + t323;
t46 = t246 * t235;
t45 = t230 * t327 + (t303 - t331 + (-t238 * t438 + t365) * qJD(3)) * t229 + (-t231 + (-pkin(6) - t399) * t229 + t250 * t230) * qJD(1) - t292;
t44 = t221 + (-t303 + (-t357 - t365) * qJD(3)) * t230 + (-t401 + (t308 - t385) * t229) * qJD(1) + t412;
t41 = t346 - t410;
t40 = (t190 * t337 + t84) * t238 + (-qJD(3) * t133 + t156 * t229 + t190 * t340) * t235;
t39 = (-qJD(3) * t190 * t230 - t82) * t238 + (qJD(3) * t135 - t156 * t230 + t318) * t235;
t30 = t270 * t335 + (qJD(1) * t269 - t229 * t82 + t230 * t84) * t235;
t25 = t229 * t84 + t230 * t82 + (t349 * t229 + t374) * qJD(1) + t322;
t24 = (qJD(3) * t305 + t388) * t238 + (-t351 * qJD(3) + t348 * t229 + t230 * t304) * t235;
t23 = (-qJD(3) * t306 - t389) * t238 + (t350 * qJD(3) - t348 * t230 + t298) * t235;
t10 = t389 * t230 + t388 * t229 + (t307 + (-t192 - t350) * t229) * qJD(1) + t322;
t9 = t246 * t335 + (t410 * qJD(1) - t389 * t229 + t388 * t230) * t235;
t4 = (t288 * qJD(3) - t38) * t238 + (-qJD(1) * t32 + qJD(3) * t87 + t17 * t229 + t18 * t230) * t235;
t3 = (t289 * qJD(3) - t37) * t238 + (-qJD(1) * t31 + qJD(3) * t86 + t15 * t229 + t16 * t230) * t235;
t2 = (t286 * qJD(3) - t36) * t238 + (-qJD(1) * t34 + qJD(3) * t89 + t13 * t229 + t14 * t230) * t235;
t1 = (t287 * qJD(3) - t35) * t238 + (-qJD(1) * t33 + qJD(3) * t88 + t11 * t229 + t12 * t230) * t235;
t5 = [(t44 * t93 + t45 * t92) * t407 + (t55 * t99 + t56 * t98) * t408 + (t106 * t147 + t107 * t146) * t409 + t414 + (-t280 + t285) * t336 + (Icges(4,1) * t235 + t281 + t382) * t335 + t431 * t311 + (t432 * t235 + t430) * t234; 0; 0; m(6) * (t108 * t44 + t109 * t45 + t57 * t92 + t58 * t93) + m(5) * (t157 * t55 + t158 * t56 + t90 * t98 + t91 * t99) + (t227 / 0.2e1 + t228 / 0.2e1) * t277 * qJD(3) + ((qJD(1) * t164 - t229 * t253) * t404 + t166 * t434 + (t373 / 0.2e1 - t368 / 0.2e1) * qJD(3) - t260 + m(4) * (-t107 * t212 - t146 * t199) + (t370 / 0.2e1 + t367 / 0.2e1 - t147 * t402 + t257) * qJD(1)) * t230 + ((-t163 * qJD(1) - t230 * t253) * t238 / 0.2e1 + t165 * t434 + (-t371 / 0.2e1 + t366 / 0.2e1) * qJD(3) + t259 + m(4) * (-t106 * t212 - t147 * t199) + (t372 / 0.2e1 + t369 / 0.2e1 + t146 * t402 + t258) * qJD(1)) * t229; m(4) * t59 + m(5) * t25 + m(6) * t10; (t10 * t41 + t108 * t58 + t109 * t57) * t407 + (t157 * t91 + t158 * t90 + t25 * t66) * t408 + t343 * t212 * t199 * t409 + (t168 * t422 + (-t95 * qJD(1) + (-t266 * qJD(1) - t137) * t230) * t230 - t427) * t230 + (t167 * t422 + (t96 * qJD(1) + (t265 * qJD(1) + t136) * t229) * t229 + ((-t137 + (-t367 - t370) * qJD(3) + t164 * t335 + t166 * t336 - t342) * t229 + (t163 * t335 + t165 * t336 + t136 - (t369 + t372) * qJD(3)) * t230 + (t97 - t94 + (t162 - t266) * t229 + t421) * qJD(1)) * t230 + t428) * t229 + (t229 * t95 - t230 * t94 + t31 + t32) * t341 + (t229 * t97 - t230 * t96 + t33 + t34) * t340; m(6) * (t23 * t93 + t24 * t92 + t44 * t65 + t45 * t64) + m(5) * (t100 * t56 + t101 * t55 + t39 * t99 + t40 * t98) + ((t229 * t258 + t230 * t257) * qJD(3) - t392) * t238 + (t259 * t230 + t260 * t229 + (-t229 * t257 + t230 * t258) * qJD(1)) * t235 - t352; m(5) * t30 + m(6) * t9; m(6) * (t10 * t46 + t108 * t23 + t109 * t24 + t41 * t9 + t57 * t64 + t58 * t65) + m(5) * (t100 * t90 + t101 * t91 + t157 * t39 + t158 * t40 + t25 * t85 + t30 * t66) + (-t4 / 0.2e1 - t3 / 0.2e1 + t325 * t335) * t230 + (t2 / 0.2e1 + t1 / 0.2e1 + t326 * t335) * t229 + ((-t229 * t325 + t230 * t326) * qJD(1) + t427 * t229 / 0.2e1 + t428 * t230 / 0.2e1 + (t390 * t229 - t391 * t230) * qJD(3) / 0.2e1) * t235 + (qJD(1) * t429 + t425 * t229 + t426 * t230) * t404 + (t424 * t229 + t423 * t230) * qJD(1) / 0.2e1; (t23 * t65 + t24 * t64 + t46 * t9) * t407 + (t100 * t40 + t101 * t39 + t30 * t85) * t408 + (t392 * t238 + (t263 * t229 - t411 * t230) * qJD(3) + t352) * t238 + ((-t425 * t238 + t1 + t2) * t230 + (t426 * t238 + t3 + t4) * t229 + (t429 * t235 + t433 * t238) * qJD(3) + (t411 * t229 + t263 * t230) * qJD(1)) * t235; m(6) * ((t229 * t93 + t230 * t92) * t335 + (t229 * t44 + t230 * t45 + (-t229 * t92 + t230 * t93) * qJD(1)) * t235); t329; m(6) * ((-t10 + (t108 * t229 + t109 * t230) * qJD(3)) * t238 + (qJD(3) * t41 + t229 * t58 + t230 * t57 + (t108 * t230 - t109 * t229) * qJD(1)) * t235); m(6) * ((-t9 + (t229 * t65 + t230 * t64) * qJD(3)) * t238 + (qJD(3) * t46 + t229 * t23 + t230 * t24 + (-t229 * t64 + t230 * t65) * qJD(1)) * t235); 0.2e1 * (-0.1e1 + t343) * t238 * t329;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
