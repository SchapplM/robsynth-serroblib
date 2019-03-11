% Calculate time derivative of joint inertia matrix for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:14:06
% EndTime: 2019-03-08 19:14:25
% DurationCPUTime: 13.89s
% Computational Cost: add. (66486->920), mult. (148045->1363), div. (0->0), fcn. (174731->14), ass. (0->373)
t372 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t437 = 0.2e1 * t372;
t410 = sin(pkin(11));
t411 = cos(pkin(11));
t418 = sin(qJ(2));
t419 = cos(qJ(2));
t339 = t418 * t410 - t419 * t411;
t312 = t339 * qJD(2);
t327 = sin(pkin(10));
t330 = cos(pkin(10));
t412 = cos(pkin(6));
t361 = t412 * t410;
t362 = t412 * t411;
t298 = t361 * t419 + t362 * t418;
t334 = qJD(2) * t298;
t254 = t327 * t312 - t330 * t334;
t256 = t312 * t330 + t327 * t334;
t335 = -t361 * t418 + t362 * t419;
t340 = t410 * t419 + t411 * t418;
t275 = -t327 * t340 + t330 * t335;
t277 = -t327 * t335 - t330 * t340;
t328 = sin(pkin(6));
t296 = t339 * t328;
t297 = t340 * t328;
t393 = pkin(12) + qJ(5);
t325 = sin(t393);
t374 = cos(t393);
t280 = t297 * t374 + t325 * t412;
t343 = -t297 * t325 + t374 * t412;
t207 = Icges(6,5) * t280 + Icges(6,6) * t343 + Icges(6,3) * t296;
t208 = Icges(6,4) * t280 + Icges(6,2) * t343 + Icges(6,6) * t296;
t209 = Icges(6,1) * t280 + Icges(6,4) * t343 + Icges(6,5) * t296;
t106 = t207 * t296 + t208 * t343 + t209 * t280;
t290 = qJD(2) * t297;
t409 = t106 * t290;
t276 = t298 * t330 - t327 * t339;
t407 = t328 * t330;
t237 = t276 * t374 - t325 * t407;
t292 = t335 * qJD(2);
t313 = t340 * qJD(2);
t255 = t292 * t330 - t313 * t327;
t188 = qJD(5) * t237 + t255 * t325;
t358 = t328 * t374;
t344 = -t276 * t325 - t330 * t358;
t189 = qJD(5) * t344 + t255 * t374;
t109 = Icges(6,5) * t189 - Icges(6,6) * t188 - Icges(6,3) * t254;
t111 = Icges(6,4) * t189 - Icges(6,2) * t188 - Icges(6,6) * t254;
t113 = Icges(6,1) * t189 - Icges(6,4) * t188 - Icges(6,5) * t254;
t145 = Icges(6,5) * t237 + Icges(6,6) * t344 - Icges(6,3) * t275;
t147 = Icges(6,4) * t237 + Icges(6,2) * t344 - Icges(6,6) * t275;
t149 = Icges(6,1) * t237 + Icges(6,4) * t344 - Icges(6,5) * t275;
t394 = qJD(2) * t328;
t291 = t339 * t394;
t232 = qJD(5) * t280 - t291 * t325;
t233 = qJD(5) * t343 - t291 * t374;
t47 = t109 * t296 + t111 * t343 + t113 * t280 + t145 * t290 - t147 * t232 + t149 * t233;
t278 = -t298 * t327 - t330 * t339;
t408 = t327 * t328;
t239 = t278 * t374 + t325 * t408;
t257 = -t292 * t327 - t313 * t330;
t190 = qJD(5) * t239 + t257 * t325;
t345 = -t278 * t325 + t327 * t358;
t191 = qJD(5) * t345 + t257 * t374;
t110 = Icges(6,5) * t191 - Icges(6,6) * t190 - Icges(6,3) * t256;
t112 = Icges(6,4) * t191 - Icges(6,2) * t190 - Icges(6,6) * t256;
t114 = Icges(6,1) * t191 - Icges(6,4) * t190 - Icges(6,5) * t256;
t146 = Icges(6,5) * t239 + Icges(6,6) * t345 - Icges(6,3) * t277;
t148 = Icges(6,4) * t239 + Icges(6,2) * t345 - Icges(6,6) * t277;
t150 = Icges(6,1) * t239 + Icges(6,4) * t345 - Icges(6,5) * t277;
t48 = t110 * t296 + t112 * t343 + t114 * t280 + t146 * t290 - t148 * t232 + t150 * t233;
t332 = sin(qJ(6));
t333 = cos(qJ(6));
t192 = -t237 * t332 - t275 * t333;
t193 = t237 * t333 - t275 * t332;
t117 = Icges(7,5) * t193 + Icges(7,6) * t192 - Icges(7,3) * t344;
t119 = Icges(7,4) * t193 + Icges(7,2) * t192 - Icges(7,6) * t344;
t121 = Icges(7,1) * t193 + Icges(7,4) * t192 - Icges(7,5) * t344;
t235 = t280 * t333 + t296 * t332;
t161 = -qJD(6) * t235 - t233 * t332 + t290 * t333;
t234 = -t280 * t332 + t296 * t333;
t162 = qJD(6) * t234 + t233 * t333 + t290 * t332;
t126 = -qJD(6) * t193 - t189 * t332 - t254 * t333;
t127 = qJD(6) * t192 + t189 * t333 - t254 * t332;
t77 = Icges(7,5) * t127 + Icges(7,6) * t126 + Icges(7,3) * t188;
t79 = Icges(7,4) * t127 + Icges(7,2) * t126 + Icges(7,6) * t188;
t81 = Icges(7,1) * t127 + Icges(7,4) * t126 + Icges(7,5) * t188;
t20 = t117 * t232 + t119 * t161 + t121 * t162 + t234 * t79 + t235 * t81 - t343 * t77;
t194 = -t239 * t332 - t277 * t333;
t195 = t239 * t333 - t277 * t332;
t118 = Icges(7,5) * t195 + Icges(7,6) * t194 - Icges(7,3) * t345;
t120 = Icges(7,4) * t195 + Icges(7,2) * t194 - Icges(7,6) * t345;
t122 = Icges(7,1) * t195 + Icges(7,4) * t194 - Icges(7,5) * t345;
t128 = -qJD(6) * t195 - t191 * t332 - t256 * t333;
t129 = qJD(6) * t194 + t191 * t333 - t256 * t332;
t78 = Icges(7,5) * t129 + Icges(7,6) * t128 + Icges(7,3) * t190;
t80 = Icges(7,4) * t129 + Icges(7,2) * t128 + Icges(7,6) * t190;
t82 = Icges(7,1) * t129 + Icges(7,4) * t128 + Icges(7,5) * t190;
t21 = t118 * t232 + t120 * t161 + t122 * t162 + t234 * t80 + t235 * t82 - t343 * t78;
t102 = Icges(7,5) * t162 + Icges(7,6) * t161 + Icges(7,3) * t232;
t103 = Icges(7,4) * t162 + Icges(7,2) * t161 + Icges(7,6) * t232;
t104 = Icges(7,1) * t162 + Icges(7,4) * t161 + Icges(7,5) * t232;
t153 = Icges(7,5) * t235 + Icges(7,6) * t234 - Icges(7,3) * t343;
t154 = Icges(7,4) * t235 + Icges(7,2) * t234 - Icges(7,6) * t343;
t155 = Icges(7,1) * t235 + Icges(7,4) * t234 - Icges(7,5) * t343;
t38 = -t102 * t343 + t103 * t234 + t104 * t235 + t153 * t232 + t154 * t161 + t155 * t162;
t64 = -t117 * t343 + t119 * t234 + t121 * t235;
t65 = -t118 * t343 + t120 * t234 + t122 * t235;
t89 = -t153 * t343 + t154 * t234 + t155 * t235;
t6 = -t20 * t275 - t21 * t277 - t254 * t64 - t256 * t65 + t290 * t89 + t296 * t38;
t157 = Icges(6,5) * t233 - Icges(6,6) * t232 + Icges(6,3) * t290;
t158 = Icges(6,4) * t233 - Icges(6,2) * t232 + Icges(6,6) * t290;
t159 = Icges(6,1) * t233 - Icges(6,4) * t232 + Icges(6,5) * t290;
t63 = t157 * t296 + t158 * t343 + t159 * t280 + t207 * t290 - t208 * t232 + t209 * t233;
t91 = t145 * t296 + t147 * t343 + t149 * t280;
t92 = t146 * t296 + t148 * t343 + t150 * t280;
t436 = -t254 * t91 - t256 * t92 - t275 * t47 - t277 * t48 + t296 * t63 + t409 + t6;
t16 = t117 * t188 + t119 * t126 + t121 * t127 + t192 * t79 + t193 * t81 - t344 * t77;
t17 = t118 * t188 + t120 * t126 + t122 * t127 + t192 * t80 + t193 * t82 - t344 * t78;
t33 = -t102 * t344 + t103 * t192 + t104 * t193 + t126 * t154 + t127 * t155 + t153 * t188;
t55 = -t117 * t344 + t119 * t192 + t121 * t193;
t56 = -t118 * t344 + t120 * t192 + t122 * t193;
t71 = -t153 * t344 + t154 * t192 + t155 * t193;
t3 = -t16 * t275 - t17 * t277 - t254 * t55 - t256 * t56 + t290 * t71 + t296 * t33;
t43 = -t109 * t275 + t111 * t344 + t113 * t237 - t145 * t254 - t147 * t188 + t149 * t189;
t44 = -t110 * t275 + t112 * t344 + t114 * t237 - t146 * t254 - t148 * t188 + t150 * t189;
t51 = -t157 * t275 + t158 * t344 + t159 * t237 - t188 * t208 + t189 * t209 - t207 * t254;
t85 = -t145 * t275 + t147 * t344 + t149 * t237;
t86 = -t146 * t275 + t148 * t344 + t150 * t237;
t99 = -t207 * t275 + t208 * t344 + t209 * t237;
t435 = -t254 * t85 - t256 * t86 - t275 * t43 - t277 * t44 + t290 * t99 + t296 * t51 + t3;
t100 = -t207 * t277 + t208 * t345 + t209 * t239;
t18 = t117 * t190 + t119 * t128 + t121 * t129 + t194 * t79 + t195 * t81 - t345 * t77;
t19 = t118 * t190 + t120 * t128 + t122 * t129 + t194 * t80 + t195 * t82 - t345 * t78;
t34 = -t102 * t345 + t103 * t194 + t104 * t195 + t128 * t154 + t129 * t155 + t153 * t190;
t57 = -t117 * t345 + t119 * t194 + t121 * t195;
t58 = -t118 * t345 + t120 * t194 + t122 * t195;
t72 = -t153 * t345 + t154 * t194 + t155 * t195;
t4 = -t18 * t275 - t19 * t277 - t254 * t57 - t256 * t58 + t290 * t72 + t296 * t34;
t45 = -t109 * t277 + t111 * t345 + t113 * t239 - t145 * t256 - t147 * t190 + t149 * t191;
t46 = -t110 * t277 + t112 * t345 + t114 * t239 - t146 * t256 - t148 * t190 + t150 * t191;
t52 = -t157 * t277 + t158 * t345 + t159 * t239 - t190 * t208 + t191 * t209 - t207 * t256;
t87 = -t145 * t277 + t147 * t345 + t149 * t239;
t88 = -t146 * t277 + t148 * t345 + t150 * t239;
t434 = t100 * t290 - t254 * t87 - t256 * t88 - t275 * t45 - t277 * t46 + t296 * t52 + t4;
t83 = rSges(7,1) * t127 + rSges(7,2) * t126 + rSges(7,3) * t188;
t414 = pkin(5) * t189 + pkin(9) * t188 + t83;
t84 = rSges(7,1) * t129 + rSges(7,2) * t128 + rSges(7,3) * t190;
t413 = pkin(5) * t191 + pkin(9) * t190 + t84;
t123 = rSges(7,1) * t193 + rSges(7,2) * t192 - rSges(7,3) * t344;
t404 = pkin(5) * t237 - pkin(9) * t344 + t123;
t124 = rSges(7,1) * t195 + rSges(7,2) * t194 - rSges(7,3) * t345;
t403 = pkin(5) * t239 - pkin(9) * t345 + t124;
t433 = 0.2e1 * m(6);
t432 = 0.2e1 * m(7);
t431 = t188 / 0.2e1;
t430 = t190 / 0.2e1;
t429 = t232 / 0.2e1;
t428 = -t344 / 0.2e1;
t427 = -t345 / 0.2e1;
t426 = -t254 / 0.2e1;
t425 = -t256 / 0.2e1;
t424 = -t275 / 0.2e1;
t423 = -t277 / 0.2e1;
t422 = -t343 / 0.2e1;
t421 = t290 / 0.2e1;
t420 = t296 / 0.2e1;
t417 = pkin(2) * t328;
t416 = pkin(2) * t419;
t329 = cos(pkin(12));
t415 = pkin(4) * t329;
t105 = rSges(7,1) * t162 + rSges(7,2) * t161 + rSges(7,3) * t232;
t405 = pkin(5) * t233 + pkin(9) * t232 + t105;
t156 = rSges(7,1) * t235 + rSges(7,2) * t234 - rSges(7,3) * t343;
t402 = pkin(5) * t280 - pkin(9) * t343 + t156;
t183 = pkin(3) * t257 - qJ(4) * t256 - qJD(4) * t277;
t366 = t412 * t419;
t315 = pkin(2) * qJD(2) * t366 - t328 * qJD(3);
t379 = qJD(2) * t418;
t371 = pkin(2) * t379;
t287 = -t327 * t315 - t330 * t371;
t283 = t412 * t287;
t401 = t412 * t183 + t283;
t225 = pkin(3) * t278 - qJ(4) * t277;
t365 = t412 * t418;
t378 = pkin(2) * t365 - qJ(3) * t328;
t265 = -t327 * t378 + t330 * t416;
t253 = t412 * t265;
t400 = t412 * t225 + t253;
t380 = qJD(2) * t419;
t314 = qJD(3) * t412 + t380 * t417;
t399 = pkin(3) * t291 - qJ(4) * t290 - qJD(4) * t296 - t314;
t264 = t327 * t416 + t330 * t378;
t398 = t264 * t408 + t265 * t407;
t319 = qJ(3) * t412 + t417 * t418;
t397 = -t297 * rSges(4,1) + t296 * rSges(4,2) - rSges(4,3) * t412 - t319;
t396 = -t297 * pkin(3) - t296 * qJ(4) - t319;
t286 = t330 * t315 - t327 * t371;
t395 = t286 * t408 + t287 * t407;
t326 = sin(pkin(12));
t392 = t326 * t408;
t391 = t326 * t407;
t389 = t412 / 0.2e1;
t139 = -pkin(8) * t256 + t257 * t415;
t388 = t412 * t139 + t401;
t144 = pkin(4) * t392 - pkin(8) * t277 + t278 * t415;
t387 = t412 * t144 + t400;
t386 = -pkin(8) * t290 + t291 * t415 + t399;
t375 = t412 * t326;
t385 = -pkin(4) * t375 - pkin(8) * t296 - t297 * t415 + t396;
t384 = t419 * Icges(3,4);
t383 = t418 * Icges(3,4);
t377 = t412 * t264;
t376 = t412 * t286;
t373 = t328 * (rSges(4,1) * t291 + rSges(4,2) * t290 - t314);
t182 = pkin(3) * t255 - qJ(4) * t254 - qJD(4) * t275;
t370 = t182 * t408 + t183 * t407 + t395;
t224 = pkin(3) * t276 - qJ(4) * t275;
t369 = t224 * t408 + t225 * t407 + t398;
t284 = -t297 * t326 + t329 * t412;
t285 = t297 * t329 + t375;
t368 = (-rSges(5,1) * t285 - rSges(5,2) * t284 - rSges(5,3) * t296 + t396) * t328;
t360 = rSges(5,1) * t329 - rSges(5,2) * t326;
t367 = (-rSges(5,3) * t290 + t291 * t360 + t399) * t328;
t357 = Icges(5,1) * t329 - Icges(5,4) * t326;
t356 = Icges(5,4) * t329 - Icges(5,2) * t326;
t355 = Icges(5,5) * t329 - Icges(5,6) * t326;
t354 = -(Icges(5,4) * t285 + Icges(5,2) * t284 + Icges(5,6) * t296) * t326 + (Icges(5,1) * t285 + Icges(5,4) * t284 + Icges(5,5) * t296) * t329;
t214 = rSges(6,1) * t280 + rSges(6,2) * t343 + rSges(6,3) * t296;
t353 = t328 * (-t214 + t385);
t160 = rSges(6,1) * t233 - rSges(6,2) * t232 + rSges(6,3) * t290;
t352 = (-t160 + t386) * t328;
t138 = -pkin(8) * t254 + t255 * t415;
t351 = t138 * t408 + t139 * t407 + t370;
t143 = -pkin(4) * t391 - pkin(8) * t275 + t276 * t415;
t350 = t143 * t408 + t144 * t407 + t369;
t349 = t328 * (t386 - t405);
t348 = (t385 - t402) * t328;
t347 = -t182 * t412 - t376;
t346 = -t224 * t412 - t377;
t308 = -t418 * t327 + t330 * t366;
t310 = -t327 * t366 - t418 * t330;
t342 = t327 * t365 - t419 * t330;
t341 = -t419 * t327 - t330 * t365;
t338 = -t138 * t412 + t347;
t337 = -t143 * t412 + t346;
t240 = -t276 * t326 - t329 * t407;
t241 = t276 * t329 - t391;
t242 = -t278 * t326 + t329 * t408;
t243 = t278 * t329 + t392;
t336 = (-(Icges(5,4) * t243 + Icges(5,2) * t242 - Icges(5,6) * t277) * t326 + (Icges(5,1) * t243 + Icges(5,4) * t242 - Icges(5,5) * t277) * t329) * t327 - (-(Icges(5,4) * t241 + Icges(5,2) * t240 - Icges(5,6) * t275) * t326 + (Icges(5,1) * t241 + Icges(5,4) * t240 - Icges(5,5) * t275) * t329) * t330;
t306 = (rSges(3,1) * t419 - rSges(3,2) * t418) * t394;
t305 = (Icges(3,1) * t419 - t383) * t394;
t304 = (-Icges(3,2) * t418 + t384) * t394;
t303 = (Icges(3,5) * t419 - Icges(3,6) * t418) * t394;
t302 = t342 * qJD(2);
t301 = t310 * qJD(2);
t300 = t341 * qJD(2);
t299 = t308 * qJD(2);
t295 = t412 * rSges(3,3) + (rSges(3,1) * t418 + rSges(3,2) * t419) * t328;
t294 = Icges(3,5) * t412 + (Icges(3,1) * t418 + t384) * t328;
t293 = Icges(3,6) * t412 + (Icges(3,2) * t419 + t383) * t328;
t274 = rSges(3,1) * t301 + rSges(3,2) * t302;
t273 = rSges(3,1) * t299 + rSges(3,2) * t300;
t272 = Icges(3,1) * t301 + Icges(3,4) * t302;
t271 = Icges(3,1) * t299 + Icges(3,4) * t300;
t270 = Icges(3,4) * t301 + Icges(3,2) * t302;
t269 = Icges(3,4) * t299 + Icges(3,2) * t300;
t268 = Icges(3,5) * t301 + Icges(3,6) * t302;
t267 = Icges(3,5) * t299 + Icges(3,6) * t300;
t263 = -rSges(3,1) * t342 + rSges(3,2) * t310 + rSges(3,3) * t408;
t262 = -rSges(3,1) * t341 + rSges(3,2) * t308 - rSges(3,3) * t407;
t261 = -Icges(3,1) * t342 + Icges(3,4) * t310 + Icges(3,5) * t408;
t260 = -Icges(3,1) * t341 + Icges(3,4) * t308 - Icges(3,5) * t407;
t259 = -Icges(3,4) * t342 + Icges(3,2) * t310 + Icges(3,6) * t408;
t258 = -Icges(3,4) * t341 + Icges(3,2) * t308 - Icges(3,6) * t407;
t251 = Icges(4,1) * t297 - Icges(4,4) * t296 + Icges(4,5) * t412;
t250 = Icges(4,4) * t297 - Icges(4,2) * t296 + Icges(4,6) * t412;
t248 = -Icges(4,1) * t291 - Icges(4,4) * t290;
t247 = -Icges(4,4) * t291 - Icges(4,2) * t290;
t246 = -Icges(4,5) * t291 - Icges(4,6) * t290;
t228 = Icges(5,5) * t290 - t291 * t357;
t227 = Icges(5,6) * t290 - t291 * t356;
t226 = Icges(5,3) * t290 - t291 * t355;
t220 = Icges(5,5) * t285 + Icges(5,6) * t284 + Icges(5,3) * t296;
t216 = rSges(4,1) * t278 + rSges(4,2) * t277 + rSges(4,3) * t408;
t215 = rSges(4,1) * t276 + rSges(4,2) * t275 - rSges(4,3) * t407;
t213 = Icges(4,1) * t278 + Icges(4,4) * t277 + Icges(4,5) * t408;
t212 = Icges(4,1) * t276 + Icges(4,4) * t275 - Icges(4,5) * t407;
t211 = Icges(4,4) * t278 + Icges(4,2) * t277 + Icges(4,6) * t408;
t210 = Icges(4,4) * t276 + Icges(4,2) * t275 - Icges(4,6) * t407;
t206 = rSges(4,1) * t257 + rSges(4,2) * t256;
t205 = rSges(4,1) * t255 + rSges(4,2) * t254;
t204 = Icges(4,1) * t257 + Icges(4,4) * t256;
t203 = Icges(4,1) * t255 + Icges(4,4) * t254;
t202 = Icges(4,4) * t257 + Icges(4,2) * t256;
t201 = Icges(4,4) * t255 + Icges(4,2) * t254;
t200 = Icges(4,5) * t257 + Icges(4,6) * t256;
t199 = Icges(4,5) * t255 + Icges(4,6) * t254;
t198 = (t273 * t327 + t274 * t330) * t328;
t178 = -rSges(5,3) * t256 + t257 * t360;
t177 = -rSges(5,3) * t254 + t255 * t360;
t176 = -Icges(5,5) * t256 + t257 * t357;
t175 = -Icges(5,5) * t254 + t255 * t357;
t174 = -Icges(5,6) * t256 + t257 * t356;
t173 = -Icges(5,6) * t254 + t255 * t356;
t172 = -Icges(5,3) * t256 + t257 * t355;
t171 = -Icges(5,3) * t254 + t255 * t355;
t170 = rSges(5,1) * t243 + rSges(5,2) * t242 - rSges(5,3) * t277;
t169 = rSges(5,1) * t241 + rSges(5,2) * t240 - rSges(5,3) * t275;
t164 = Icges(5,5) * t243 + Icges(5,6) * t242 - Icges(5,3) * t277;
t163 = Icges(5,5) * t241 + Icges(5,6) * t240 - Icges(5,3) * t275;
t152 = rSges(6,1) * t239 + rSges(6,2) * t345 - rSges(6,3) * t277;
t151 = rSges(6,1) * t237 + rSges(6,2) * t344 - rSges(6,3) * t275;
t134 = -t205 * t412 + t330 * t373 - t376;
t133 = t206 * t412 + t327 * t373 + t283;
t125 = (t205 * t327 + t206 * t330) * t328 + t395;
t116 = rSges(6,1) * t191 - rSges(6,2) * t190 - rSges(6,3) * t256;
t115 = rSges(6,1) * t189 - rSges(6,2) * t188 - rSges(6,3) * t254;
t108 = t152 * t296 + t214 * t277;
t107 = -t151 * t296 - t214 * t275;
t101 = -t151 * t277 + t152 * t275;
t98 = -t169 * t412 + t330 * t368 + t346;
t97 = t170 * t412 + t327 * t368 + t400;
t96 = -t177 * t412 + t330 * t367 + t347;
t95 = t178 * t412 + t327 * t367 + t401;
t94 = -t124 * t343 + t156 * t345;
t93 = t123 * t343 - t156 * t344;
t90 = (t169 * t327 + t170 * t330) * t328 + t369;
t76 = -t123 * t345 + t124 * t344;
t75 = (t177 * t327 + t178 * t330) * t328 + t370;
t74 = t277 * t402 + t296 * t403;
t73 = -t275 * t402 - t296 * t404;
t70 = -t151 * t412 + t330 * t353 + t337;
t69 = t152 * t412 + t327 * t353 + t387;
t68 = t116 * t296 + t152 * t290 + t160 * t277 + t214 * t256;
t67 = -t115 * t296 - t151 * t290 - t160 * t275 - t214 * t254;
t66 = t275 * t403 - t277 * t404;
t62 = (t151 * t327 + t152 * t330) * t328 + t350;
t61 = -t115 * t412 + t330 * t352 + t338;
t60 = t116 * t412 + t327 * t352 + t388;
t59 = -t115 * t277 + t116 * t275 - t151 * t256 + t152 * t254;
t54 = t330 * t348 - t404 * t412 + t337;
t53 = t327 * t348 + t403 * t412 + t387;
t50 = (t115 * t327 + t116 * t330) * t328 + t351;
t49 = (t327 * t404 + t330 * t403) * t328 + t350;
t42 = t105 * t345 + t124 * t232 - t156 * t190 - t343 * t84;
t41 = -t105 * t344 - t123 * t232 + t156 * t188 + t343 * t83;
t40 = t330 * t349 - t412 * t414 + t338;
t39 = t327 * t349 + t412 * t413 + t388;
t37 = t123 * t190 - t124 * t188 + t344 * t84 - t345 * t83;
t36 = t256 * t402 + t277 * t405 + t290 * t403 + t296 * t413;
t35 = -t254 * t402 - t275 * t405 - t290 * t404 - t296 * t414;
t32 = (t327 * t414 + t330 * t413) * t328 + t351;
t31 = t89 * t412 + (t327 * t65 - t330 * t64) * t328;
t30 = -t275 * t64 - t277 * t65 + t296 * t89;
t29 = -t343 * t89 - t344 * t64 - t345 * t65;
t28 = t254 * t403 - t256 * t404 + t275 * t413 - t277 * t414;
t27 = t72 * t412 + (t327 * t58 - t330 * t57) * t328;
t26 = t71 * t412 + (t327 * t56 - t330 * t55) * t328;
t25 = -t275 * t57 - t277 * t58 + t296 * t72;
t24 = -t275 * t55 - t277 * t56 + t296 * t71;
t23 = -t343 * t72 - t344 * t57 - t345 * t58;
t22 = -t343 * t71 - t344 * t55 - t345 * t56;
t15 = t63 * t412 + (t327 * t48 - t330 * t47) * t328;
t14 = t52 * t412 + (t327 * t46 - t330 * t45) * t328;
t13 = t51 * t412 + (t327 * t44 - t330 * t43) * t328;
t9 = t38 * t412 + (-t20 * t330 + t21 * t327) * t328;
t8 = t34 * t412 + (-t18 * t330 + t19 * t327) * t328;
t7 = t33 * t412 + (-t16 * t330 + t17 * t327) * t328;
t5 = t188 * t64 + t190 * t65 - t20 * t344 - t21 * t345 + t232 * t89 - t343 * t38;
t2 = -t18 * t344 + t188 * t57 - t19 * t345 + t190 * t58 + t232 * t72 - t34 * t343;
t1 = -t16 * t344 - t17 * t345 + t188 * t55 + t190 * t56 + t232 * t71 - t33 * t343;
t10 = [0; m(3) * t198 + m(4) * t125 + m(5) * t75 + m(6) * t50 + m(7) * t32; -t7 * t407 + t8 * t408 - t13 * t407 + t14 * t408 + ((t259 * t302 + t261 * t301 + t268 * t408 + t270 * t310 - t272 * t342) * t408 - (t258 * t302 + t260 * t301 + t267 * t408 + t269 * t310 - t271 * t342) * t407 + (t293 * t302 + t294 * t301 + t303 * t408 + t304 * t310 - t305 * t342) * t412) * t408 + ((t200 * t408 + t202 * t277 + t204 * t278 + t211 * t256 + t213 * t257) * t408 - (t199 * t408 + t201 * t277 + t203 * t278 + t210 * t256 + t212 * t257) * t407 + (t246 * t408 + t247 * t277 + t248 * t278 + t250 * t256 + t251 * t257) * t412) * t408 - ((t259 * t300 + t261 * t299 - t268 * t407 + t270 * t308 - t272 * t341) * t408 - (t258 * t300 + t260 * t299 - t267 * t407 + t269 * t308 - t271 * t341) * t407 + (t293 * t300 + t294 * t299 - t303 * t407 + t304 * t308 - t305 * t341) * t412) * t407 - ((-t200 * t407 + t202 * t275 + t204 * t276 + t211 * t254 + t213 * t255) * t408 - (-t199 * t407 + t201 * t275 + t203 * t276 + t210 * t254 + t212 * t255) * t407 + (-t246 * t407 + t247 * t275 + t248 * t276 + t250 * t254 + t251 * t255) * t412) * t407 - ((-t254 * t220 - t275 * t226 + t240 * t227 + t241 * t228 + t255 * t354) * t412 + ((-t164 * t254 - t172 * t275 + t174 * t240 + t176 * t241) * t327 - (-t163 * t254 - t171 * t275 + t173 * t240 + t175 * t241) * t330 + t336 * t255) * t328) * t407 + ((-t256 * t220 - t277 * t226 + t242 * t227 + t243 * t228 + t257 * t354) * t412 + ((-t164 * t256 - t172 * t277 + t174 * t242 + t176 * t243) * t327 - (-t163 * t256 - t171 * t277 + t173 * t242 + t175 * t243) * t330 + t336 * t257) * t328) * t408 + t412 * t9 + (t32 * t49 + t39 * t53 + t40 * t54) * t432 + t412 * t15 + (t50 * t62 + t60 * t69 + t61 * t70) * t433 + 0.2e1 * m(5) * (t75 * t90 + t95 * t97 + t96 * t98) + 0.2e1 * m(4) * ((-t215 * t412 - t377) * t134 + (t216 * t412 + t253) * t133 + t398 * t125 + ((t216 * t125 + t134 * t397) * t330 + (t215 * t125 + t133 * t397) * t327) * t328) + t412 * (((-t259 * t379 + t261 * t380 + t270 * t419 + t272 * t418) * t327 - (-t258 * t379 + t260 * t380 + t269 * t419 + t271 * t418) * t330) * t328 ^ 2 + ((-t267 * t330 + t268 * t327 - t293 * t379 + t294 * t380 + t304 * t419 + t305 * t418) * t328 + t303 * t412) * t412) + t412 * ((t412 * t246 - t296 * t247 + t297 * t248 - t290 * t250 - t291 * t251) * t412 + ((t200 * t412 - t296 * t202 + t297 * t204 - t290 * t211 - t291 * t213) * t327 - (t199 * t412 - t296 * t201 + t297 * t203 - t290 * t210 - t291 * t212) * t330) * t328) + t412 * ((t220 * t290 + t226 * t296 + t227 * t284 + t228 * t285 - t291 * t354) * t412 + ((t164 * t290 + t172 * t296 + t174 * t284 + t176 * t285) * t327 - (t163 * t290 + t171 * t296 + t173 * t284 + t175 * t285) * t330 - t336 * t291) * t328) + 0.2e1 * m(3) * ((-t262 * t412 - t295 * t407) * (-t273 * t412 - t306 * t407) + (t263 * t412 - t295 * t408) * (t274 * t412 - t306 * t408) + (t262 * t327 + t263 * t330) * t328 * t198); 0; m(7) * (t412 * t32 + (t327 * t40 - t330 * t39) * t328) + m(6) * (t412 * t50 + (t327 * t61 - t330 * t60) * t328) + m(5) * (t412 * t75 + (t327 * t96 - t330 * t95) * t328) + m(4) * (t412 * t125 + (-t133 * t330 + t134 * t327) * t328); 0; t290 * t437; m(7) * (-t254 * t53 - t256 * t54 - t275 * t39 - t277 * t40 + t290 * t49 + t296 * t32) + m(6) * (-t254 * t69 - t256 * t70 - t275 * t60 - t277 * t61 + t290 * t62 + t296 * t50) + m(5) * (-t254 * t97 - t256 * t98 - t275 * t95 - t277 * t96 + t290 * t90 + t296 * t75); (t290 * t412 + (t254 * t330 - t256 * t327) * t328) * t437; 0.4e1 * t372 * (t254 * t275 + t256 * t277 + t290 * t296); m(6) * t59 + m(7) * t28; m(7) * (t28 * t49 + t32 * t66 + t35 * t54 + t36 * t53 + t39 * t74 + t40 * t73) + m(6) * (t101 * t50 + t107 * t61 + t108 * t60 + t59 * t62 + t67 * t70 + t68 * t69) + (t26 + t99 * t412 + (t327 * t86 - t330 * t85) * t328) * t426 + (t27 + t100 * t412 + (t327 * t88 - t330 * t87) * t328) * t425 + (t7 + t13) * t424 + (t8 + t14) * t423 + (t31 + t106 * t412 + (t327 * t92 - t330 * t91) * t328) * t421 + (t9 + t15) * t420 + t436 * t389 + t434 * t408 / 0.2e1 - t435 * t407 / 0.2e1; m(6) * (t59 * t412 + (t327 * t67 - t330 * t68) * t328) + m(7) * (t28 * t412 + (t327 * t35 - t330 * t36) * t328); m(6) * (t101 * t290 - t107 * t256 - t108 * t254 - t275 * t68 - t277 * t67 + t296 * t59) + m(7) * (-t254 * t74 - t256 * t73 - t275 * t36 - t277 * t35 + t28 * t296 + t290 * t66); t290 * t30 + (t409 + t436) * t296 + (-t290 * t92 - t434) * t277 + (-t290 * t91 - t435) * t275 + (-t100 * t296 + t275 * t87 + t277 * t88 - t25) * t256 + (t275 * t85 + t277 * t86 - t296 * t99 - t24) * t254 + (t28 * t66 + t35 * t73 + t36 * t74) * t432 + (t101 * t59 + t107 * t67 + t108 * t68) * t433; m(7) * t37; t5 * t389 + t31 * t429 + t9 * t422 + t26 * t431 + t7 * t428 + t27 * t430 + t8 * t427 + m(7) * (t32 * t76 + t37 * t49 + t39 * t94 + t40 * t93 + t41 * t54 + t42 * t53) + (t327 * t2 / 0.2e1 - t330 * t1 / 0.2e1) * t328; m(7) * (t37 * t412 + (t327 * t41 - t330 * t42) * t328); m(7) * (-t254 * t94 - t256 * t93 - t275 * t42 - t277 * t41 + t290 * t76 + t296 * t37); t22 * t426 + t1 * t424 + t29 * t421 + t5 * t420 + t25 * t430 + t4 * t427 + m(7) * (t28 * t76 + t35 * t93 + t36 * t94 + t37 * t66 + t41 * t73 + t42 * t74) + t30 * t429 + t6 * t422 + t23 * t425 + t2 * t423 + t24 * t431 + t3 * t428; t190 * t23 - t345 * t2 + t188 * t22 - t344 * t1 + t232 * t29 - t343 * t5 + (t37 * t76 + t41 * t93 + t42 * t94) * t432;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
