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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:46:55
% EndTime: 2020-01-03 11:47:17
% DurationCPUTime: 10.69s
% Computational Cost: add. (9208->386), mult. (7890->515), div. (0->0), fcn. (6082->8), ass. (0->227)
t419 = Icges(5,4) + Icges(6,4);
t418 = Icges(5,1) + Icges(6,1);
t417 = -Icges(5,2) - Icges(6,2);
t202 = qJ(3) + qJ(4);
t195 = cos(t202);
t416 = t419 * t195;
t194 = sin(t202);
t415 = t419 * t194;
t414 = Icges(5,5) + Icges(6,5);
t413 = Icges(5,6) + Icges(6,6);
t412 = t417 * t194 + t416;
t411 = t418 * t195 - t415;
t410 = Icges(5,3) + Icges(6,3);
t201 = qJ(1) + pkin(8);
t192 = sin(t201);
t193 = cos(t201);
t384 = t412 * t192 - t413 * t193;
t383 = t411 * t192 - t414 * t193;
t408 = -t413 * t194 + t414 * t195;
t377 = t414 * t192 + t411 * t193;
t376 = t413 * t192 + t412 * t193;
t403 = t410 * t192 + t408 * t193;
t409 = t192 * t403;
t407 = t408 * t192 - t410 * t193;
t395 = t417 * t195 - t415;
t394 = t418 * t194 + t416;
t405 = -t376 * t194 + t377 * t195;
t404 = t384 * t194 - t383 * t195;
t406 = t414 * t194 + t413 * t195;
t200 = qJD(3) + qJD(4);
t301 = t193 * t200;
t402 = -t384 * qJD(1) + t395 * t301;
t302 = t192 * t200;
t401 = t376 * qJD(1) + t395 * t302;
t400 = t383 * qJD(1) + t394 * t301;
t399 = -t377 * qJD(1) + t394 * t302;
t205 = cos(qJ(3));
t188 = t205 * pkin(3) + pkin(2);
t166 = pkin(4) * t195 + t188;
t249 = rSges(6,1) * t195 - rSges(6,2) * t194;
t223 = t166 + t249;
t398 = t223 * t193;
t397 = t412 * t200;
t396 = t411 * t200;
t393 = t404 * t193;
t368 = -t395 * t194 - t394 * t195;
t392 = t407 * qJD(1);
t391 = t405 * t192;
t390 = t404 * t192 + t193 * t407;
t389 = -t193 * t403 + t391;
t388 = -t192 * t407 + t393;
t387 = t193 * t405 + t409;
t386 = -qJD(1) * t403 + t302 * t406;
t385 = -t301 * t406 - t392;
t207 = -pkin(7) - pkin(6);
t199 = -qJ(5) + t207;
t303 = t192 * t195;
t304 = t192 * t194;
t329 = rSges(6,3) * t193;
t382 = rSges(6,1) * t303 - rSges(6,2) * t304 + t192 * t166 + t193 * t199 - t329;
t381 = -t200 * t377 - t402;
t380 = t200 * t383 + t401;
t379 = t200 * t376 + t400;
t378 = -t200 * t384 - t399;
t184 = t193 * t207;
t293 = t192 * t188 + t184;
t374 = -t293 + t382;
t168 = t193 * t188;
t290 = t199 - t207;
t373 = t168 - t398 + (-rSges(6,3) + t290) * t192;
t144 = t249 * t200;
t299 = t195 * t200;
t372 = pkin(4) * t299 + t144;
t157 = rSges(6,1) * t194 + rSges(6,2) * t195;
t370 = pkin(4) * t194 + t157;
t369 = -t368 * qJD(1) - t408 * t200;
t203 = sin(qJ(3));
t280 = qJD(3) * t203;
t275 = pkin(3) * t280;
t300 = t194 * t200;
t149 = -pkin(4) * t300 - t275;
t283 = qJD(1) * t193;
t270 = t195 * t283;
t284 = qJD(1) * t192;
t367 = rSges(6,1) * t270 + rSges(6,3) * t284 + t192 * t149 + t166 * t283;
t366 = (-t395 * t200 - t396) * t195 + (t394 * t200 + t397) * t194 - t406 * qJD(1);
t198 = cos(qJ(1)) * pkin(1);
t189 = qJD(1) * t198;
t365 = ((t401 * t194 + t399 * t195 + t384 * t299 + t383 * t300 - t392) * t193 + t387 * qJD(1)) * t193 + (t385 * t192 + (t381 * t194 - t379 * t195 + t386) * t193 + ((t403 - t404) * t193 - t388 - t391) * qJD(1)) * t192;
t342 = (t386 * t193 + (-t389 - t393) * qJD(1)) * t193 + ((t402 * t194 + t400 * t195 + t376 * t299 + t377 * t300) * t192 + (-t380 * t194 + t378 * t195 + t385) * t193 + (-t409 + (-t405 + t407) * t193 + t390) * qJD(1)) * t192;
t364 = t387 * t192 + t388 * t193;
t333 = rSges(4,2) * t203;
t334 = rSges(4,1) * t205;
t178 = t192 * t334;
t331 = rSges(4,3) * t193;
t357 = t331 - t178;
t134 = -t192 * t333 - t357;
t251 = -t333 + t334;
t135 = -rSges(4,3) * t192 - t193 * t251;
t177 = rSges(4,1) * t203 + rSges(4,2) * t205;
t221 = qJD(3) * t177;
t210 = rSges(4,3) * t284 - t192 * t221 + t283 * t334;
t213 = t193 * t221;
t22 = (qJD(1) * t135 + t210) * t192 + (-t213 + (t134 + t357) * qJD(1)) * t193;
t345 = 2 * m(4);
t363 = t22 * t345;
t323 = Icges(4,4) * t205;
t241 = -Icges(4,2) * t203 + t323;
t130 = -Icges(4,6) * t193 + t192 * t241;
t324 = Icges(4,4) * t203;
t245 = Icges(4,1) * t205 - t324;
t132 = -Icges(4,5) * t193 + t192 * t245;
t229 = t130 * t203 - t132 * t205;
t362 = t192 * t229;
t250 = rSges(5,1) * t195 - rSges(5,2) * t194;
t359 = t250 * t193;
t237 = Icges(4,5) * t205 - Icges(4,6) * t203;
t354 = Icges(4,3) * t192 + t193 * t237;
t351 = Icges(4,6) * t192 + t193 * t241;
t348 = Icges(4,5) * t192 + t193 * t245;
t344 = 2 * m(5);
t343 = 2 * m(6);
t190 = t192 ^ 2;
t191 = t193 ^ 2;
t341 = -t192 / 0.2e1;
t340 = -t193 / 0.2e1;
t339 = (t389 * t192 + t390 * t193) * t284;
t338 = m(4) * t177;
t337 = pkin(3) * t203;
t186 = t192 * pkin(6);
t187 = t193 * pkin(2);
t196 = sin(qJ(1)) * pkin(1);
t278 = qJD(5) * t192;
t335 = t278 - (-t149 - t275) * t193 - t157 * t301 + (-t290 * t193 + t329 + (t188 - t223) * t192) * qJD(1);
t330 = rSges(5,3) * t193;
t328 = rSges(5,3) - t207;
t327 = rSges(6,3) - t199;
t326 = t374 * t192;
t309 = t130 * t205;
t308 = t205 * t351;
t307 = t132 * t203;
t306 = t203 * t348;
t145 = t250 * t200;
t305 = t145 * t192;
t298 = t372 * t193;
t294 = rSges(5,1) * t270 + rSges(5,3) * t284;
t109 = t370 * t193;
t292 = t186 + t187;
t291 = t190 + t191;
t128 = -Icges(4,3) * t193 + t192 * t237;
t287 = qJD(1) * t128;
t282 = qJD(1) * t194;
t281 = qJD(3) * t192;
t279 = qJD(3) * t205;
t106 = -pkin(2) * t192 + pkin(6) * t193 + t293;
t107 = t192 * t207 - t168 + t292;
t159 = t188 * t283;
t277 = t192 * (-t192 * t275 + t159 + (-t187 + (-pkin(6) - t207) * t192) * qJD(1)) + t107 * t284 + t106 * t283;
t274 = pkin(3) * t279;
t123 = rSges(5,1) * t303 - rSges(5,2) * t304 - t330;
t125 = -t192 * rSges(5,3) - t359;
t271 = t193 * t282;
t214 = -t192 * t299 - t271;
t272 = t192 * t300;
t273 = t125 * t284 + t123 * t283 + t192 * (-rSges(5,1) * t272 + rSges(5,2) * t214 + t294);
t269 = pkin(2) - t333;
t158 = rSges(5,1) * t194 + rSges(5,2) * t195;
t126 = (-t158 - t337) * t192;
t108 = t370 * t192;
t246 = t192 * t365 + t339;
t244 = Icges(4,1) * t203 + t323;
t240 = Icges(4,2) * t205 + t324;
t236 = Icges(4,5) * t203 + Icges(4,6) * t205;
t228 = -t203 * t351 + t205 * t348;
t225 = t373 * t284 + t374 * t283 + (-rSges(6,1) * t272 + rSges(6,2) * t214 - qJD(5) * t193 - t159 + (-qJD(1) * t290 + t275) * t192 + t367) * t192;
t224 = -t370 - t337;
t222 = t157 * t200;
t218 = t229 * t193;
t217 = qJD(3) * t244;
t216 = qJD(3) * t240;
t215 = qJD(1) * t224;
t209 = -t158 * t200 - t275;
t208 = (t369 * t192 + t366 * t193 + t379 * t194 + t381 * t195) * t341 + (-t366 * t192 + t369 * t193 + t378 * t194 + t380 * t195) * t340 + (-t368 * t192 - t193 * t406 + t383 * t194 + t384 * t195) * t284 / 0.2e1 - (-t192 * t406 + t368 * t193 - t377 * t194 - t376 * t195) * t283 / 0.2e1;
t185 = pkin(6) * t283;
t180 = t193 * t337;
t173 = t193 * t274;
t165 = t251 * qJD(3);
t127 = t158 * t193 + t180;
t105 = t192 * t123;
t99 = t192 * t106;
t98 = t180 + t109;
t97 = t224 * t192;
t94 = -t135 + t198 + t292;
t93 = t178 + t196 + (-rSges(4,3) - pkin(6)) * t193 + t269 * t192;
t90 = t193 * t275 + t185 + (t184 + (-pkin(2) + t188) * t192) * qJD(1);
t85 = qJD(1) * t354 - t236 * t281;
t84 = qJD(3) * t193 * t236 + t287;
t82 = t328 * t192 + t168 + t198 + t359;
t81 = t123 + t196 + t293;
t77 = t327 * t192 + t198 + t398;
t76 = t196 + t382;
t75 = t158 * t301 + (t192 * t250 - t330) * qJD(1);
t61 = -t158 * t283 - t305 + (-t192 * t279 - t203 * t283) * pkin(3);
t60 = qJD(1) * t126 + t145 * t193 + t173;
t57 = -t125 * t193 + t105;
t56 = t189 + (t193 * t269 + t186) * qJD(1) + t210;
t55 = t185 - t213 + (t331 - t196 + (-pkin(2) - t251) * t192) * qJD(1);
t54 = pkin(4) * t214 - t144 * t192 - t157 * t283;
t53 = -qJD(1) * t108 + t298;
t43 = (-t274 - t372) * t192 + t193 * t215;
t42 = t192 * t215 + t173 + t298;
t40 = -rSges(5,2) * t271 + t159 + t189 + (-qJD(1) * t207 + t209) * t192 + t294;
t39 = t209 * t193 + (-t196 + t328 * t193 + (-t188 - t250) * t192) * qJD(1);
t38 = t192 * t354 + t193 * t228;
t37 = -t192 * t128 + t218;
t36 = -t192 * t228 + t193 * t354;
t35 = -t128 * t193 - t362;
t26 = t189 + (-rSges(6,2) * t282 - qJD(5)) * t193 + (-qJD(1) * t199 - t222) * t192 + t367;
t25 = t278 + (t149 - t222) * t193 + (-t223 * t192 + t327 * t193 - t196) * qJD(1);
t24 = t105 + t99 + (-t107 - t125) * t193;
t23 = -t193 * t373 + t326;
t21 = -t193 * t75 + t273;
t8 = t99 + (-t107 - t373) * t193 + t326;
t7 = (-t75 - t90) * t193 + t273 + t277;
t6 = t335 * t193 + t225;
t5 = (-t90 + t335) * t193 + t225 + t277;
t1 = [t205 * t217 + t245 * t280 - t203 * t216 + t241 * t279 + (t55 * t94 + t56 * t93) * t345 + (t39 * t82 + t40 * t81) * t344 + (t25 * t77 + t26 * t76) * t343 + t395 * t300 + t394 * t299 + t397 * t195 + t396 * t194; 0; 0; t208 + ((-t93 * t338 + t309 / 0.2e1 + t307 / 0.2e1) * t192 + (-t94 * t338 + t308 / 0.2e1 + t306 / 0.2e1) * t193) * qJD(1) + m(4) * ((-t192 * t55 + t193 * t56) * t177 + (-t192 * t94 + t193 * t93) * t165) + m(6) * (t25 * t97 + t26 * t98 + t42 * t76 + t43 * t77) + m(5) * (t126 * t39 + t127 * t40 + t60 * t81 + t61 * t82) + (-qJD(3) * t228 + t203 * (qJD(1) * t132 + t193 * t217) + t205 * (qJD(1) * t130 + t193 * t216)) * t341 + (-qJD(3) * t229 + t203 * (qJD(1) * t348 - t244 * t281) + t205 * (qJD(1) * t351 - t240 * t281)) * t340 + (t191 / 0.2e1 + t190 / 0.2e1) * t237 * qJD(3); m(4) * t22 + m(5) * t7 + m(6) * t5; (t42 * t98 + t43 * t97 + t5 * t8) * t343 + (t126 * t61 + t127 * t60 + t24 * t7) * t344 + t291 * t177 * t165 * t345 + t339 + (-t135 * t363 - t35 * t284 + (-(-t36 + t218) * qJD(1) - t193 * t85) * t193 + t342) * t193 + (-t36 * t284 + t134 * t363 + (-t37 * qJD(1) + (-t228 * qJD(1) - t84) * t192) * t192 + ((-t85 - (t306 + t308) * qJD(3) + t351 * t279 + t348 * t280) * t192 + (t130 * t279 + t132 * t280 - t287 - t84 - (t307 + t309) * qJD(3)) * t193 + (t38 - t35 - t362 + (t128 - t228) * t193) * qJD(1)) * t193 + t365) * t192 + (t192 * t38 + t193 * t37 + t364) * t283; t208 + m(6) * (-t108 * t25 + t109 * t26 + t53 * t76 + t54 * t77) + ((-t192 * t39 + t193 * t40 + (-t192 * t81 - t193 * t82) * qJD(1)) * t158 + (-t192 * t82 + t193 * t81) * t145) * m(5); m(5) * t21 + m(6) * t6; m(6) * (-t108 * t43 + t109 * t42 + t23 * t5 + t53 * t98 + t54 * t97 + t6 * t8) + m(5) * (-t126 * t305 + t21 * t24 + t57 * t7 + (-t127 * t284 - t192 * t61) * t158) + (m(5) * (t127 * t145 + t158 * t60) + t342 + (-m(5) * t126 * t158 + t364) * qJD(1)) * t193 + t246; (t145 * t158 * t291 + t21 * t57) * t344 + (-t108 * t54 + t109 * t53 + t23 * t6) * t343 + (t364 * qJD(1) + t342) * t193 + t246; m(6) * (-t192 * t26 - t193 * t25 + (t192 * t77 - t193 * t76) * qJD(1)); 0; m(6) * (-t192 * t42 - t193 * t43 + (t192 * t97 - t193 * t98) * qJD(1)); m(6) * (-t192 * t53 - t193 * t54 + (-t108 * t192 - t109 * t193) * qJD(1)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
