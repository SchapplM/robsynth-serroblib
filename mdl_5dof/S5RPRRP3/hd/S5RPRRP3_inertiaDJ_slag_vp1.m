% Calculate time derivative of joint inertia matrix for
% S5RPRRP3
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:35
% EndTime: 2022-01-23 09:29:51
% DurationCPUTime: 9.37s
% Computational Cost: add. (9208->377), mult. (7890->511), div. (0->0), fcn. (6082->8), ass. (0->227)
t415 = Icges(5,4) + Icges(6,4);
t414 = Icges(5,1) + Icges(6,1);
t413 = Icges(5,2) + Icges(6,2);
t196 = qJ(3) + qJ(4);
t190 = cos(t196);
t412 = t415 * t190;
t189 = sin(t196);
t411 = t415 * t189;
t410 = Icges(5,5) + Icges(6,5);
t409 = Icges(5,6) + Icges(6,6);
t408 = -t189 * t413 + t412;
t407 = t190 * t414 - t411;
t406 = Icges(5,3) + Icges(6,3);
t195 = qJ(1) + pkin(8);
t187 = sin(t195);
t188 = cos(t195);
t378 = t187 * t409 + t188 * t408;
t376 = t187 * t410 + t188 * t407;
t404 = t190 * t413 + t411;
t403 = t189 * t414 + t412;
t405 = -t189 * t409 + t190 * t410;
t379 = t187 * t408 - t188 * t409;
t377 = t187 * t407 - t188 * t410;
t364 = t189 * t410 + t190 * t409;
t401 = t187 * t405 - t188 * t406;
t393 = t187 * t406 + t188 * t405;
t402 = t189 * t378 - t190 * t376;
t194 = qJD(3) + qJD(4);
t398 = t403 * t194;
t397 = t404 * t194;
t400 = t364 * t194;
t399 = t189 * t379 - t190 * t377;
t396 = qJD(1) * t378 - t187 * t397;
t395 = -qJD(1) * t376 + t187 * t398;
t392 = t408 * t194;
t391 = t407 * t194;
t390 = t189 * t404 - t190 * t403;
t388 = t393 * qJD(1);
t387 = t401 * qJD(1);
t386 = t402 * t187;
t385 = t187 * t399 + t188 * t401;
t384 = -t188 * t393 - t386;
t383 = -t188 * t400 - t387;
t382 = t187 * t400 - t388;
t381 = -qJD(1) * t379 - t188 * t397;
t380 = qJD(1) * t377 + t188 * t398;
t199 = cos(qJ(3));
t184 = pkin(3) * t199 + pkin(2);
t163 = pkin(4) * t190 + t184;
t180 = t187 * rSges(6,3);
t297 = t188 * t190;
t298 = t188 * t189;
t375 = rSges(6,1) * t297 - rSges(6,2) * t298 + t163 * t188 + t180;
t374 = t194 * t377 + t396;
t373 = -t194 * t379 - t395;
t372 = t399 * t188;
t197 = sin(qJ(3));
t277 = qJD(3) * t197;
t271 = pkin(3) * t277;
t295 = t189 * t194;
t147 = -pkin(4) * t295 - t271;
t154 = rSges(6,1) * t189 + rSges(6,2) * t190;
t220 = t154 * t194;
t371 = qJD(5) * t188 - (t147 - t220) * t187;
t370 = -t187 * t401 + t372;
t369 = t187 * t393 - t188 * t402;
t328 = rSges(6,1) * t190;
t248 = -rSges(6,2) * t189 + t328;
t201 = -pkin(7) - pkin(6);
t193 = -qJ(5) + t201;
t285 = t193 - t201;
t289 = t163 - t184;
t91 = t187 * t289 + t188 * t285;
t368 = -rSges(6,3) * t188 + t187 * t248 + t91;
t164 = t188 * t184;
t367 = -t187 * t285 - t164 + t375;
t155 = rSges(5,1) * t189 + rSges(5,2) * t190;
t221 = t155 * t194;
t366 = t187 * t221;
t363 = (t391 - t397) * t190 + (-t392 - t398) * t189 + t364 * qJD(1);
t361 = -t194 * t376 - t381;
t360 = -t194 * t378 - t380;
t279 = qJD(1) * t189;
t268 = t187 * t279;
t280 = qJD(1) * t188;
t359 = rSges(6,2) * t268 + rSges(6,3) * t280 + qJD(5) * t187 + t147 * t188;
t358 = qJD(1) * t390 + t194 * t405;
t294 = t190 * t194;
t357 = (t382 * t188 + (-t372 - t384) * qJD(1)) * t188 + ((t189 * t381 + t190 * t380 + t294 * t378 + t295 * t376 - t388) * t187 + (-t374 * t189 + t373 * t190 + t383) * t188 + ((t402 + t401) * t188 + t385) * qJD(1)) * t187;
t327 = rSges(4,2) * t197;
t330 = rSges(4,1) * t199;
t250 = -t327 + t330;
t326 = rSges(4,3) * t188;
t134 = t187 * t250 - t326;
t272 = t188 * t327;
t182 = t187 * rSges(4,3);
t287 = t188 * t330 + t182;
t135 = -t272 + t287;
t173 = rSges(4,1) * t197 + rSges(4,2) * t199;
t219 = qJD(3) * t173;
t278 = qJD(1) * t197;
t267 = t187 * t278;
t205 = rSges(4,2) * t267 + rSges(4,3) * t280 - t188 * t219;
t355 = t187 * t219;
t22 = (qJD(1) * t134 + t205) * t188 + (-t355 + (-t135 - t272 + t182) * qJD(1)) * t187;
t343 = 2 * m(4);
t356 = t22 * t343;
t333 = pkin(2) - t184;
t352 = t187 * t333;
t322 = Icges(4,4) * t199;
t241 = -Icges(4,2) * t197 + t322;
t131 = Icges(4,6) * t187 + t188 * t241;
t323 = Icges(4,4) * t197;
t245 = Icges(4,1) * t199 - t323;
t133 = Icges(4,5) * t187 + t188 * t245;
t229 = t131 * t197 - t133 * t199;
t349 = t229 * t188;
t183 = t188 * pkin(6);
t334 = sin(qJ(1)) * pkin(1);
t348 = t183 - t334;
t237 = Icges(4,5) * t199 - Icges(4,6) * t197;
t128 = -Icges(4,3) * t188 + t187 * t237;
t130 = -Icges(4,6) * t188 + t187 * t241;
t132 = -Icges(4,5) * t188 + t187 * t245;
t342 = 2 * m(5);
t341 = 2 * m(6);
t185 = t187 ^ 2;
t186 = t188 ^ 2;
t340 = t187 / 0.2e1;
t339 = -t188 / 0.2e1;
t338 = m(4) * t173;
t337 = m(5) * t155;
t336 = pkin(3) * t197;
t335 = t187 * pkin(6);
t192 = cos(qJ(1)) * pkin(1);
t332 = -pkin(6) - t201;
t296 = t188 * t201;
t106 = t183 + t296 - t352;
t107 = -pkin(2) * t188 + t187 * t332 + t164;
t331 = t106 * t187 + t107 * t188;
t329 = rSges(5,1) * t190;
t181 = t187 * rSges(5,3);
t325 = rSges(6,3) - t193;
t308 = t130 * t199;
t307 = t131 * t199;
t306 = t132 * t197;
t305 = t133 * t197;
t249 = -rSges(5,2) * t189 + t329;
t144 = t249 * t194;
t304 = t144 * t187;
t299 = t187 * t193;
t123 = -rSges(5,3) * t188 + t187 * t249;
t125 = rSges(5,1) * t297 - rSges(5,2) * t298 + t181;
t57 = t123 * t187 + t125 * t188;
t281 = qJD(1) * t187;
t292 = pkin(4) * t268 + t154 * t281;
t290 = rSges(5,2) * t268 + rSges(5,3) * t280;
t288 = t187 * t271 + t201 * t281;
t286 = t185 + t186;
t129 = Icges(4,3) * t187 + t188 * t237;
t282 = qJD(1) * t129;
t276 = qJD(3) * t199;
t254 = t188 * t271;
t274 = t187 * ((-t188 * t333 - t335) * qJD(1) - t288) + t188 * (-t254 + (t188 * t332 + t352) * qJD(1)) + t106 * t280;
t208 = -t188 * t295 - t190 * t281;
t269 = t188 * t294;
t273 = t123 * t280 + t187 * (-t366 + (t188 * t249 + t181) * qJD(1)) + t188 * (rSges(5,1) * t208 - rSges(5,2) * t269 + t290);
t270 = pkin(3) * t276;
t264 = -t155 - t336;
t263 = -pkin(4) * t189 - t154;
t23 = t187 * t368 + t188 * t367;
t143 = t248 * t194;
t253 = -pkin(4) * t294 - t143;
t252 = (t370 * t280 + t385 * t281) * t188 + ((t383 * t187 + (-t370 + t386) * qJD(1)) * t187 + t384 * t281 + t369 * t280 + (t382 * t187 + (t379 * t294 + t377 * t295 - t387) * t188 + (t360 * t187 + t188 * t395) * t190 + (t361 * t187 + t188 * t396) * t189 + ((t393 - t399) * t187 + t369) * qJD(1)) * t188) * t187;
t244 = Icges(4,1) * t197 + t322;
t240 = Icges(4,2) * t199 + t323;
t230 = t130 * t197 - t132 * t199;
t226 = t368 * t280 + (rSges(6,1) * t208 - rSges(6,2) * t269 - qJD(1) * t91 + t254 + t359) * t188 + (t288 + (t180 - t299 + (t248 + t289) * t188) * qJD(1) - t371) * t187;
t225 = -pkin(2) - t250;
t224 = t263 - t336;
t223 = -t184 - t249;
t222 = -t163 - t248;
t212 = qJD(3) * t244;
t211 = qJD(3) * t240;
t210 = qJD(3) * (-Icges(4,5) * t197 - Icges(4,6) * t199);
t97 = t224 * t188;
t207 = t253 - t270;
t206 = t188 * t357 + t252;
t204 = (t187 * t358 + t188 * t363 + t189 * t360 - t190 * t361) * t340 + (t363 * t187 - t358 * t188 + t189 * t373 + t190 * t374) * t339 + (-t187 * t390 - t364 * t188 + t189 * t377 + t190 * t379) * t281 / 0.2e1 + (t364 * t187 - t188 * t390 + t189 * t376 + t190 * t378) * t280 / 0.2e1;
t169 = pkin(3) * t267;
t162 = t250 * qJD(3);
t127 = t264 * t188;
t126 = t264 * t187;
t109 = t263 * t188;
t108 = t263 * t187;
t96 = t224 * t187;
t94 = t335 + t192 + (pkin(2) - t327) * t188 + t287;
t93 = t187 * t225 + t326 + t348;
t86 = t187 * t210 + t282;
t85 = -qJD(1) * t128 + t188 * t210;
t82 = -t187 * t201 + t125 + t164 + t192;
t81 = -t334 + (rSges(5,3) - t201) * t188 + t223 * t187;
t77 = t192 - t299 + t375;
t76 = t187 * t222 + t188 * t325 - t334;
t63 = -t155 * t280 - t304 + (-t187 * t276 - t188 * t278) * pkin(3);
t62 = t155 * t281 + t169 + (-t144 - t270) * t188;
t56 = t355 + (-t192 + (-rSges(4,3) - pkin(6)) * t187 + t225 * t188) * qJD(1);
t55 = ((-pkin(2) - t330) * t187 + t348) * qJD(1) + t205;
t54 = -t154 * t280 - t143 * t187 + (-t187 * t294 - t188 * t279) * pkin(4);
t53 = t188 * t253 + t292;
t44 = qJD(1) * t97 + t187 * t207;
t43 = t188 * t207 + t169 + t292;
t40 = t366 + (t188 * t223 - t181 - t192) * qJD(1) + t288;
t39 = (-t221 - t271) * t188 + (-t334 - t296 + (-t184 - t329) * t187) * qJD(1) + t290;
t38 = t129 * t187 - t349;
t37 = t128 * t187 - t188 * t230;
t36 = -t129 * t188 - t187 * t229;
t35 = -t128 * t188 - t187 * t230;
t26 = (-t187 * t325 + t188 * t222 - t192) * qJD(1) + t371;
t25 = -t188 * t220 + (-t334 - t188 * t193 + (-t163 - t328) * t187) * qJD(1) + t359;
t24 = t57 + t331;
t21 = -t125 * t281 + t273;
t8 = t23 + t331;
t7 = (-t107 - t125) * t281 + t273 + t274;
t6 = -t281 * t367 + t226;
t5 = (-t107 - t367) * t281 + t226 + t274;
t1 = [(t25 * t77 + t26 * t76) * t341 + (t39 * t82 + t40 * t81) * t342 + (t55 * t94 + t56 * t93) * t343 - t404 * t295 + t403 * t294 + (t245 - t240) * t277 + (t244 + t241) * t276 + t392 * t190 + t391 * t189; 0; 0; (-qJD(3) * t229 + t197 * (-qJD(1) * t132 - t188 * t212) + t199 * (-qJD(1) * t130 - t188 * t211)) * t340 + (-qJD(3) * t230 + t197 * (qJD(1) * t133 - t187 * t212) + t199 * (qJD(1) * t131 - t187 * t211)) * t339 + ((-t94 * t338 + t307 / 0.2e1 + t305 / 0.2e1) * t188 + (t93 * t338 + t308 / 0.2e1 + t306 / 0.2e1) * t187) * qJD(1) + m(4) * ((-t187 * t55 - t188 * t56) * t173 + (-t187 * t94 - t188 * t93) * t162) + t204 + m(6) * (t25 * t96 + t26 * t97 + t43 * t76 + t44 * t77) + m(5) * (t126 * t39 + t127 * t40 + t62 * t81 + t63 * t82) + (t186 / 0.2e1 + t185 / 0.2e1) * t237 * qJD(3); m(4) * t22 + m(5) * t7 + m(6) * t5; (t43 * t97 + t44 * t96 + t5 * t8) * t341 + (t126 * t63 + t127 * t62 + t24 * t7) * t342 + t286 * t173 * t162 * t343 + t252 + (t134 * t356 + t38 * t280 + t36 * t281 + (t37 * qJD(1) + (t229 * qJD(1) + t85) * t187) * t187) * t187 + (-t37 * t280 - t35 * t281 + t135 * t356 + (-t36 * qJD(1) + (-qJD(1) * t230 - t86) * t188) * t188 + ((t130 * t276 + t132 * t277 + t85 - (t306 + t308) * qJD(3)) * t188 + (-t86 + (-t305 - t307) * qJD(3) + t131 * t276 + t133 * t277 - t282) * t187 + (t38 - t35 + t349 + (t129 - t230) * t187) * qJD(1)) * t187 + t357) * t188; m(6) * (t108 * t25 + t109 * t26 + t53 * t76 + t54 * t77) + t204 + m(5) * (-t187 * t82 - t188 * t81) * t144 + (-t187 * t39 - t188 * t40 + (t187 * t81 - t188 * t82) * qJD(1)) * t337; m(5) * t21 + m(6) * t6; m(6) * (t108 * t44 + t109 * t43 + t23 * t5 + t53 * t97 + t54 * t96 + t6 * t8) + m(5) * (-t127 * t144 * t188 - t126 * t304 + t21 * t24 + t57 * t7) + (-t187 * t63 - t188 * t62 + (-t126 * t188 + t127 * t187) * qJD(1)) * t337 + t206; (t144 * t155 * t286 + t21 * t57) * t342 + (t108 * t54 + t109 * t53 + t23 * t6) * t341 + t206; m(6) * (t187 * t26 - t188 * t25 + (t187 * t77 + t188 * t76) * qJD(1)); 0; m(6) * (t187 * t43 - t188 * t44 + (t187 * t96 + t188 * t97) * qJD(1)); m(6) * (t187 * t53 - t188 * t54 + (t108 * t187 + t109 * t188) * qJD(1)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
