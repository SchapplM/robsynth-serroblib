% Calculate time derivative of joint inertia matrix for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR14_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR14_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR14_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:57
% EndTime: 2019-12-31 19:17:32
% DurationCPUTime: 17.23s
% Computational Cost: add. (107505->888), mult. (322494->1178), div. (0->0), fcn. (388054->14), ass. (0->365)
t349 = cos(qJ(1));
t347 = sin(qJ(1));
t440 = cos(pkin(11));
t442 = cos(pkin(5));
t392 = t442 * t440;
t438 = sin(pkin(11));
t369 = t347 * t438 - t349 * t392;
t343 = sin(pkin(5));
t441 = cos(pkin(6));
t404 = t343 * t441;
t439 = sin(pkin(6));
t312 = -t349 * t404 + t369 * t439;
t472 = t312 * pkin(8);
t346 = sin(qJ(3));
t389 = t441 * t440;
t391 = t442 * t439;
t451 = cos(qJ(3));
t471 = (-t438 * t346 + t451 * t389) * t343 + t451 * t391;
t364 = t369 * t441;
t403 = t343 * t439;
t393 = t349 * t403;
t390 = t442 * t438;
t368 = -t347 * t440 - t349 * t390;
t408 = t368 * t451;
t297 = -t408 + (-t364 - t393) * t346;
t345 = sin(qJ(4));
t450 = cos(qJ(4));
t272 = t297 * t450 + t312 * t345;
t357 = t451 * t364;
t384 = t451 * t403;
t296 = -t346 * t368 + t349 * t384 + t357;
t378 = -t297 * t345 + t312 * t450;
t200 = Icges(5,5) * t272 + Icges(5,6) * t378 + Icges(5,3) * t296;
t202 = Icges(5,4) * t272 + Icges(5,2) * t378 + Icges(5,6) * t296;
t204 = Icges(5,1) * t272 + Icges(5,4) * t378 + Icges(5,5) * t296;
t327 = -t347 * t390 + t349 * t440;
t370 = t347 * t392 + t349 * t438;
t362 = t370 * t441;
t299 = t327 * t451 + (t347 * t403 - t362) * t346;
t361 = t370 * t439;
t395 = t347 * t404;
t313 = t361 + t395;
t274 = t299 * t450 + t313 * t345;
t356 = t451 * t362;
t298 = t327 * t346 - t347 * t384 + t356;
t377 = -t299 * t345 + t313 * t450;
t100 = t200 * t298 + t202 * t377 + t204 * t274;
t201 = Icges(5,5) * t274 + Icges(5,6) * t377 + Icges(5,3) * t298;
t203 = Icges(5,4) * t274 + Icges(5,2) * t377 + Icges(5,6) * t298;
t205 = Icges(5,1) * t274 + Icges(5,4) * t377 + Icges(5,5) * t298;
t101 = t201 * t298 + t203 * t377 + t205 * t274;
t311 = t346 * t391 + (t346 * t389 + t438 * t451) * t343;
t325 = -t403 * t440 + t441 * t442;
t295 = t311 * t450 + t325 * t345;
t376 = -t311 * t345 + t325 * t450;
t228 = Icges(5,5) * t295 + Icges(5,6) * t376 - Icges(5,3) * t471;
t229 = Icges(5,4) * t295 + Icges(5,2) * t376 - Icges(5,6) * t471;
t230 = Icges(5,1) * t295 + Icges(5,4) * t376 - Icges(5,5) * t471;
t114 = t228 * t298 + t229 * t377 + t230 * t274;
t322 = t368 * qJD(1);
t365 = qJD(1) * t369;
t358 = t441 * t365;
t375 = qJD(1) * t384;
t250 = qJD(3) * t299 + t322 * t346 - t349 * t375 - t358 * t451;
t323 = t327 * qJD(1);
t419 = qJD(3) * t346;
t252 = t323 * t346 + qJD(1) * t356 - t393 * t419 - t347 * t375 + (-t346 * t364 - t408) * qJD(3);
t374 = qJD(3) * t384;
t423 = qJD(1) * t346;
t383 = t403 * t423;
t251 = -qJD(3) * t356 + t322 * t451 - t327 * t419 + t346 * t358 + t347 * t374 + t349 * t383;
t305 = t312 * qJD(1);
t191 = qJD(4) * t377 + t251 * t450 - t305 * t345;
t344 = sin(qJ(5));
t348 = cos(qJ(5));
t227 = t274 * t348 + t298 * t344;
t142 = -qJD(5) * t227 - t191 * t344 + t250 * t348;
t226 = -t274 * t344 + t298 * t348;
t143 = qJD(5) * t226 + t191 * t348 + t250 * t344;
t224 = -t272 * t344 + t296 * t348;
t225 = t272 * t348 + t296 * t344;
t158 = Icges(6,5) * t225 + Icges(6,6) * t224 - Icges(6,3) * t378;
t160 = Icges(6,4) * t225 + Icges(6,2) * t224 - Icges(6,6) * t378;
t162 = Icges(6,1) * t225 + Icges(6,4) * t224 - Icges(6,5) * t378;
t190 = qJD(4) * t274 + t251 * t345 + t305 * t450;
t253 = -qJD(3) * t357 + t323 * t451 + t347 * t383 - t349 * t374 - t362 * t423 + t368 * t419;
t306 = t313 * qJD(1);
t193 = qJD(4) * t378 + t253 * t450 + t306 * t345;
t144 = -qJD(5) * t225 - t193 * t344 + t252 * t348;
t145 = qJD(5) * t224 + t193 * t348 + t252 * t344;
t192 = qJD(4) * t272 + t253 * t345 - t306 * t450;
t78 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t192;
t80 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t192;
t82 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t192;
t16 = t142 * t160 + t143 * t162 + t158 * t190 + t226 * t80 + t227 * t82 - t377 * t78;
t159 = Icges(6,5) * t227 + Icges(6,6) * t226 - Icges(6,3) * t377;
t161 = Icges(6,4) * t227 + Icges(6,2) * t226 - Icges(6,6) * t377;
t163 = Icges(6,1) * t227 + Icges(6,4) * t226 - Icges(6,5) * t377;
t77 = Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t190;
t79 = Icges(6,4) * t143 + Icges(6,2) * t142 + Icges(6,6) * t190;
t81 = Icges(6,1) * t143 + Icges(6,4) * t142 + Icges(6,5) * t190;
t17 = t142 * t161 + t143 * t163 + t159 * t190 + t226 * t79 + t227 * t81 - t377 * t77;
t303 = t471 * qJD(3);
t262 = qJD(4) * t376 + t303 * t450;
t265 = t295 * t348 - t344 * t471;
t304 = t311 * qJD(3);
t210 = -qJD(5) * t265 - t262 * t344 + t304 * t348;
t264 = -t295 * t344 - t348 * t471;
t211 = qJD(5) * t264 + t262 * t348 + t304 * t344;
t261 = qJD(4) * t295 + t303 * t345;
t135 = Icges(6,5) * t211 + Icges(6,6) * t210 + Icges(6,3) * t261;
t136 = Icges(6,4) * t211 + Icges(6,2) * t210 + Icges(6,6) * t261;
t137 = Icges(6,1) * t211 + Icges(6,4) * t210 + Icges(6,5) * t261;
t194 = Icges(6,5) * t265 + Icges(6,6) * t264 - Icges(6,3) * t376;
t195 = Icges(6,4) * t265 + Icges(6,2) * t264 - Icges(6,6) * t376;
t196 = Icges(6,1) * t265 + Icges(6,4) * t264 - Icges(6,5) * t376;
t30 = -t135 * t377 + t136 * t226 + t137 * t227 + t142 * t195 + t143 * t196 + t190 * t194;
t64 = -t158 * t377 + t160 * t226 + t162 * t227;
t65 = -t159 * t377 + t161 * t226 + t163 * t227;
t90 = -t194 * t377 + t195 * t226 + t196 * t227;
t3 = t16 * t296 + t17 * t298 + t250 * t65 + t252 * t64 - t30 * t471 + t304 * t90;
t120 = Icges(5,5) * t193 - Icges(5,6) * t192 + Icges(5,3) * t252;
t122 = Icges(5,4) * t193 - Icges(5,2) * t192 + Icges(5,6) * t252;
t124 = Icges(5,1) * t193 - Icges(5,4) * t192 + Icges(5,5) * t252;
t40 = t120 * t298 + t122 * t377 + t124 * t274 - t190 * t202 + t191 * t204 + t200 * t250;
t119 = Icges(5,5) * t191 - Icges(5,6) * t190 + Icges(5,3) * t250;
t121 = Icges(5,4) * t191 - Icges(5,2) * t190 + Icges(5,6) * t250;
t123 = Icges(5,1) * t191 - Icges(5,4) * t190 + Icges(5,5) * t250;
t41 = t119 * t298 + t121 * t377 + t123 * t274 - t190 * t203 + t191 * t205 + t201 * t250;
t212 = Icges(5,5) * t262 - Icges(5,6) * t261 + Icges(5,3) * t304;
t213 = Icges(5,4) * t262 - Icges(5,2) * t261 + Icges(5,6) * t304;
t214 = Icges(5,1) * t262 - Icges(5,4) * t261 + Icges(5,5) * t304;
t60 = -t190 * t229 + t191 * t230 + t212 * t298 + t213 * t377 + t214 * t274 + t228 * t250;
t470 = t100 * t252 + t101 * t250 + t114 * t304 + t296 * t40 + t298 * t41 - t471 * t60 + t3;
t113 = t228 * t296 + t229 * t378 + t230 * t272;
t18 = t144 * t160 + t145 * t162 + t158 * t192 + t224 * t80 + t225 * t82 - t378 * t78;
t19 = t144 * t161 + t145 * t163 + t159 * t192 + t224 * t79 + t225 * t81 - t378 * t77;
t31 = -t135 * t378 + t136 * t224 + t137 * t225 + t144 * t195 + t145 * t196 + t192 * t194;
t62 = -t158 * t378 + t160 * t224 + t162 * t225;
t63 = -t159 * t378 + t161 * t224 + t163 * t225;
t89 = -t194 * t378 + t195 * t224 + t196 * t225;
t4 = t18 * t296 + t19 * t298 + t250 * t63 + t252 * t62 + t304 * t89 - t31 * t471;
t42 = t120 * t296 + t122 * t378 + t124 * t272 - t192 * t202 + t193 * t204 + t200 * t252;
t43 = t119 * t296 + t121 * t378 + t123 * t272 - t192 * t203 + t193 * t205 + t201 * t252;
t61 = -t192 * t229 + t193 * t230 + t212 * t296 + t213 * t378 + t214 * t272 + t228 * t252;
t98 = t200 * t296 + t202 * t378 + t204 * t272;
t99 = t201 * t296 + t203 * t378 + t205 * t272;
t469 = t113 * t304 + t250 * t99 + t252 * t98 + t296 * t42 + t298 * t43 - t471 * t61 + t4;
t106 = -t200 * t471 + t202 * t376 + t204 * t295;
t107 = -t201 * t471 + t203 * t376 + t205 * t295;
t130 = -t228 * t471 + t229 * t376 + t230 * t295;
t128 = t130 * t304;
t71 = -t212 * t471 + t213 * t376 + t295 * t214 + t304 * t228 - t261 * t229 + t262 * t230;
t445 = -t471 * t71 + t128;
t50 = -t120 * t471 + t122 * t376 + t124 * t295 + t200 * t304 - t202 * t261 + t204 * t262;
t51 = -t119 * t471 + t121 * t376 + t123 * t295 + t201 * t304 - t203 * t261 + t205 * t262;
t20 = t158 * t261 + t160 * t210 + t162 * t211 + t264 * t80 + t265 * t82 - t376 * t78;
t21 = t159 * t261 + t161 * t210 + t163 * t211 + t264 * t79 + t265 * t81 - t376 * t77;
t47 = -t135 * t376 + t264 * t136 + t265 * t137 + t261 * t194 + t210 * t195 + t211 * t196;
t95 = -t194 * t376 + t195 * t264 + t196 * t265;
t446 = t95 * t304 - t47 * t471;
t74 = -t158 * t376 + t160 * t264 + t162 * t265;
t75 = -t159 * t376 + t161 * t264 + t163 * t265;
t8 = t20 * t296 + t21 * t298 + t75 * t250 + t74 * t252 + t446;
t468 = t106 * t252 + t107 * t250 + t50 * t296 + t51 * t298 + t445 + t8;
t400 = 2 * m(4);
t181 = t251 * rSges(4,1) - t250 * rSges(4,2) - t305 * rSges(4,3);
t182 = t253 * rSges(4,1) - t252 * rSges(4,2) + t306 * rSges(4,3);
t239 = rSges(4,1) * t297 - rSges(4,2) * t296 + rSges(4,3) * t312;
t240 = t299 * rSges(4,1) - t298 * rSges(4,2) + t313 * rSges(4,3);
t97 = -t181 * t312 + t182 * t313 - t239 * t305 - t240 * t306;
t467 = t400 * t97;
t24 = -t376 * t89 - t377 * t63 - t378 * t62;
t464 = t24 / 0.2e1;
t25 = -t376 * t90 - t377 * t65 - t378 * t64;
t463 = t25 / 0.2e1;
t462 = t190 / 0.2e1;
t461 = t192 / 0.2e1;
t460 = t261 / 0.2e1;
t459 = -t378 / 0.2e1;
t458 = -t377 / 0.2e1;
t457 = -t376 / 0.2e1;
t456 = t304 / 0.2e1;
t455 = t312 / 0.2e1;
t454 = t313 / 0.2e1;
t453 = t325 / 0.2e1;
t452 = -rSges(6,3) - pkin(10);
t449 = pkin(4) * t272;
t448 = t193 * pkin(4);
t342 = t349 * pkin(1);
t447 = t95 * t261 - t376 * t47;
t83 = t143 * rSges(6,1) + t142 * rSges(6,2) + t190 * rSges(6,3);
t444 = t191 * pkin(4) + t190 * pkin(10) + t83;
t388 = -t145 * rSges(6,1) - t144 * rSges(6,2);
t84 = t192 * rSges(6,3) - t388;
t443 = t192 * pkin(10) + t448 + t84;
t233 = Icges(4,5) * t297 - Icges(4,6) * t296 + Icges(4,3) * t312;
t437 = t233 * t313;
t234 = Icges(4,5) * t299 - Icges(4,6) * t298 + Icges(4,3) * t313;
t436 = t234 * t312;
t435 = t343 * t347;
t434 = t343 * t349;
t235 = Icges(4,4) * t297 - Icges(4,2) * t296 + Icges(4,6) * t312;
t236 = Icges(4,4) * t299 - Icges(4,2) * t298 + Icges(4,6) * t313;
t237 = Icges(4,1) * t297 - Icges(4,4) * t296 + Icges(4,5) * t312;
t238 = Icges(4,1) * t299 - Icges(4,4) * t298 + Icges(4,5) * t313;
t433 = -t235 * t298 - t236 * t296 + t237 * t299 + t238 * t297 + t436 + t437;
t138 = rSges(6,1) * t211 + rSges(6,2) * t210 + rSges(6,3) * t261;
t432 = pkin(4) * t262 + pkin(10) * t261 + t138;
t387 = -rSges(6,1) * t225 - rSges(6,2) * t224;
t164 = -rSges(6,3) * t378 - t387;
t431 = -pkin(10) * t378 + t164 + t449;
t165 = t227 * rSges(6,1) + t226 * rSges(6,2) - rSges(6,3) * t377;
t430 = t274 * pkin(4) - pkin(10) * t377 + t165;
t209 = t253 * pkin(3) + t252 * pkin(9);
t255 = pkin(3) * t297 + t296 * pkin(9);
t429 = t313 * t209 - t305 * t255;
t197 = rSges(6,1) * t265 + rSges(6,2) * t264 - rSges(6,3) * t376;
t428 = pkin(4) * t295 - pkin(10) * t376 + t197;
t207 = t274 * rSges(5,1) + rSges(5,2) * t377 + t298 * rSges(5,3);
t256 = t299 * pkin(3) + t298 * pkin(9);
t427 = -t207 - t256;
t231 = rSges(5,1) * t295 + rSges(5,2) * t376 - rSges(5,3) * t471;
t286 = pkin(3) * t311 - pkin(9) * t471;
t426 = -t231 - t286;
t283 = pkin(3) * t303 + pkin(9) * t304;
t425 = t312 * t283 + t306 * t286;
t424 = qJ(2) * t435 + t342;
t422 = qJD(1) * t347;
t421 = qJD(1) * t349;
t420 = qJD(2) * t343;
t416 = t20 / 0.2e1 + t31 / 0.2e1;
t415 = t30 / 0.2e1 + t21 / 0.2e1;
t414 = t74 / 0.2e1 + t89 / 0.2e1;
t413 = t90 / 0.2e1 + t75 / 0.2e1;
t412 = -t256 - t430;
t125 = t191 * rSges(5,1) - t190 * rSges(5,2) + t250 * rSges(5,3);
t411 = -t286 - t428;
t410 = m(5) * t442;
t409 = m(6) * t442;
t406 = t343 * t421;
t405 = -t347 * pkin(1) + qJ(2) * t434;
t208 = t251 * pkin(3) + t250 * pkin(9);
t398 = 0.2e1 * m(5);
t396 = 0.2e1 * m(6);
t276 = Icges(4,4) * t311 + Icges(4,2) * t471 + Icges(4,6) * t325;
t277 = Icges(4,1) * t311 + Icges(4,4) * t471 + Icges(4,5) * t325;
t279 = Icges(4,5) * t303 - Icges(4,6) * t304;
t280 = Icges(4,4) * t303 - Icges(4,2) * t304;
t281 = Icges(4,1) * t303 - Icges(4,4) * t304;
t386 = -t304 * t276 + t303 * t277 + t325 * t279 + t280 * t471 + t311 * t281;
t385 = -pkin(1) * t422 + qJ(2) * t406 + t347 * t420;
t382 = t61 / 0.2e1 + t50 / 0.2e1 + t416;
t381 = t51 / 0.2e1 + t60 / 0.2e1 + t415;
t380 = t113 / 0.2e1 + t106 / 0.2e1 + t414;
t379 = t107 / 0.2e1 + t114 / 0.2e1 + t413;
t126 = t193 * rSges(5,1) - t192 * rSges(5,2) + t252 * rSges(5,3);
t206 = rSges(5,1) * t272 + rSges(5,2) * t378 + rSges(5,3) * t296;
t371 = pkin(2) * t368 + t405 - t472;
t367 = -t255 + t371;
t366 = t370 * rSges(3,2);
t359 = pkin(8) * t361;
t339 = t349 * t420;
t355 = -t323 * pkin(2) + t339 + (-t359 - t342 + (-pkin(8) * t441 - qJ(2)) * t435) * qJD(1);
t354 = t327 * pkin(2) + pkin(8) * t395 + t359 + t424;
t353 = -t209 + t355;
t352 = t256 + t354;
t351 = t322 * pkin(2) - t472 * qJD(1) + t385;
t350 = t208 + t351;
t293 = t327 * rSges(3,1) + rSges(3,3) * t435 - t366 + t424;
t292 = rSges(3,1) * t368 + rSges(3,2) * t369 + rSges(3,3) * t434 + t405;
t285 = -t323 * rSges(3,1) - pkin(1) * t421 + qJD(1) * t366 + t339 + (-rSges(3,3) - qJ(2)) * t343 * t422;
t284 = t322 * rSges(3,1) + rSges(3,2) * t365 + rSges(3,3) * t406 + t385;
t282 = rSges(4,1) * t303 - rSges(4,2) * t304;
t278 = rSges(4,1) * t311 + rSges(4,2) * t471 + rSges(4,3) * t325;
t275 = Icges(4,5) * t311 + Icges(4,6) * t471 + Icges(4,3) * t325;
t263 = t312 * t286;
t242 = t325 * t256;
t241 = t313 * t255;
t217 = t354 + t240;
t216 = -t239 + t371;
t215 = rSges(5,1) * t262 - rSges(5,2) * t261 + rSges(5,3) * t304;
t185 = t325 * t208;
t180 = t240 * t325 - t278 * t313;
t179 = -t239 * t325 + t278 * t312;
t178 = Icges(4,1) * t253 - Icges(4,4) * t252 + Icges(4,5) * t306;
t177 = Icges(4,1) * t251 - Icges(4,4) * t250 - Icges(4,5) * t305;
t176 = Icges(4,4) * t253 - Icges(4,2) * t252 + Icges(4,6) * t306;
t175 = Icges(4,4) * t251 - Icges(4,2) * t250 - Icges(4,6) * t305;
t174 = Icges(4,5) * t253 - Icges(4,6) * t252 + Icges(4,3) * t306;
t173 = Icges(4,5) * t251 - Icges(4,6) * t250 - Icges(4,3) * t305;
t169 = -t182 + t355;
t168 = t351 + t181;
t167 = t275 * t313 - t276 * t298 + t277 * t299;
t166 = t275 * t312 - t276 * t296 + t277 * t297;
t157 = t352 + t207;
t156 = -t206 + t367;
t155 = -t207 * t471 - t231 * t298;
t154 = t206 * t471 + t231 * t296;
t149 = t234 * t325 + t236 * t471 + t238 * t311;
t148 = t233 * t325 + t235 * t471 + t237 * t311;
t147 = -t182 * t325 + t278 * t306 + t282 * t312;
t146 = t181 * t325 + t278 * t305 - t282 * t313;
t139 = t386 * t325;
t129 = t206 * t298 - t207 * t296;
t118 = t207 * t325 + t313 * t426 + t242;
t117 = t231 * t312 + t263 + (-t206 - t255) * t325;
t112 = t352 + t430;
t111 = -t378 * t452 + t367 + t387 - t449;
t110 = -t165 * t376 + t197 * t377;
t109 = t164 * t376 - t197 * t378;
t108 = t206 * t313 + t312 * t427 + t241;
t105 = -t252 * t276 + t253 * t277 + t275 * t306 + t279 * t312 - t280 * t296 + t281 * t297;
t104 = -t250 * t276 + t251 * t277 - t275 * t305 + t279 * t313 - t280 * t298 + t281 * t299;
t103 = -t126 + t353;
t102 = t350 + t125;
t96 = -t164 * t377 + t165 * t378;
t92 = -t298 * t428 - t430 * t471;
t91 = t296 * t428 + t431 * t471;
t88 = t313 * t411 + t325 * t430 + t242;
t87 = t263 + t428 * t312 + (-t255 - t431) * t325;
t86 = t173 * t325 + t175 * t471 + t177 * t311 - t236 * t304 + t238 * t303;
t85 = t174 * t325 + t176 * t471 + t178 * t311 - t235 * t304 + t237 * t303;
t76 = -t296 * t430 + t298 * t431;
t73 = t126 * t471 - t206 * t304 + t215 * t296 + t231 * t252;
t72 = -t125 * t471 + t207 * t304 - t215 * t298 - t231 * t250;
t70 = t71 * t325;
t68 = t215 * t312 + t231 * t306 + (-t126 - t209) * t325 + t425;
t67 = t125 * t325 + t185 + (-t215 - t283) * t313 - t426 * t305;
t66 = t312 * t412 + t313 * t431 + t241;
t59 = -t125 * t296 + t126 * t298 + t206 * t250 - t207 * t252;
t58 = t192 * t452 + t353 + t388 - t448;
t57 = t350 + t444;
t56 = t100 * t312 + t101 * t313 + t114 * t325;
t55 = t113 * t325 + t312 * t98 + t313 * t99;
t54 = t126 * t313 - t206 * t305 + (-t125 - t208) * t312 + t427 * t306 + t429;
t53 = t100 * t296 + t101 * t298 - t114 * t471;
t52 = -t113 * t471 + t296 * t98 + t298 * t99;
t49 = -t138 * t378 - t164 * t261 + t192 * t197 + t376 * t84;
t48 = t138 * t377 + t165 * t261 - t190 * t197 - t376 * t83;
t46 = t47 * t325;
t39 = t432 * t312 + t428 * t306 + (-t209 - t443) * t325 + t425;
t38 = t185 + t444 * t325 + (-t283 - t432) * t313 - t411 * t305;
t37 = t312 * t74 + t313 * t75 + t325 * t95;
t36 = t296 * t74 + t298 * t75 - t471 * t95;
t35 = t252 * t428 + t296 * t432 - t304 * t431 + t443 * t471;
t34 = -t250 * t428 - t298 * t432 + t304 * t430 - t444 * t471;
t33 = -t376 * t95 - t377 * t75 - t378 * t74;
t32 = t164 * t190 - t165 * t192 - t377 * t84 + t378 * t83;
t29 = t312 * t64 + t313 * t65 + t325 * t90;
t28 = t312 * t62 + t313 * t63 + t325 * t89;
t27 = t296 * t64 + t298 * t65 - t471 * t90;
t26 = t296 * t62 + t298 * t63 - t471 * t89;
t23 = t250 * t431 - t252 * t430 - t296 * t444 + t298 * t443;
t22 = t443 * t313 - t431 * t305 + (-t208 - t444) * t312 + t412 * t306 + t429;
t15 = t106 * t306 - t107 * t305 + t50 * t312 + t51 * t313 + t70;
t13 = -t305 * t99 + t306 * t98 + t312 * t42 + t313 * t43 + t325 * t61;
t12 = t100 * t306 - t101 * t305 + t312 * t40 + t313 * t41 + t325 * t60;
t9 = t20 * t312 + t21 * t313 - t75 * t305 + t74 * t306 + t46;
t7 = t75 * t190 + t74 * t192 - t20 * t378 - t21 * t377 + t447;
t6 = t18 * t312 + t19 * t313 - t305 * t63 + t306 * t62 + t31 * t325;
t5 = t16 * t312 + t17 * t313 + t30 * t325 - t305 * t65 + t306 * t64;
t2 = -t18 * t378 - t19 * t377 + t190 * t63 + t192 * t62 + t261 * t89 - t31 * t376;
t1 = -t16 * t378 - t17 * t377 + t190 * t65 + t192 * t64 + t261 * t90 - t30 * t376;
t10 = [(t168 * t217 + t169 * t216) * t400 + 0.2e1 * m(3) * (t284 * t293 + t285 * t292) + (t102 * t157 + t103 * t156) * t398 + (t111 * t58 + t112 * t57) * t396 + t71 + t47 + t386; ((t111 * t421 + t112 * t422 + t347 * t58 - t349 * t57) * m(6) + (-t102 * t349 + t103 * t347 + t156 * t421 + t157 * t422) * m(5) + (-t168 * t349 + t169 * t347 + t216 * t421 + t217 * t422) * m(4) + m(3) * (-t284 * t349 + t285 * t347 + t292 * t421 + t293 * t422)) * t343; 0; t46 + t139 + t70 + (t146 * t217 + t147 * t216 + t168 * t180 + t169 * t179) * m(4) + (t102 * t118 + t103 * t117 + t156 * t68 + t157 * t67) * m(5) + (t111 * t39 + t112 * t38 + t57 * t88 + t58 * t87) * m(6) + (t86 / 0.2e1 + t104 / 0.2e1 + t381) * t313 + (t105 / 0.2e1 + t85 / 0.2e1 + t382) * t312 + (t166 / 0.2e1 + t148 / 0.2e1 + t380) * t306 - (t149 / 0.2e1 + t167 / 0.2e1 + t379) * t305; m(4) * t97 * t442 + t54 * t410 + t22 * t409 + ((-t146 * t349 + t147 * t347 + t179 * t421 + t180 * t422) * m(4) + (t117 * t421 + t118 * t422 + t347 * t68 - t349 * t67) * m(5) + (t347 * t39 - t349 * t38 + t421 * t87 + t422 * t88) * m(6)) * t343; (t22 * t66 + t38 * t88 + t39 * t87) * t396 + (t108 * t54 + t117 * t68 + t118 * t67) * t398 + (t180 * t146 + t179 * t147) * t400 + (t139 + t15 + t9) * t325 + (t28 + t55 + (t148 + t166) * t325) * t306 - (t29 + t56 + (t149 + t167) * t325) * t305 + (t5 + t12 + t239 * t467 + (t104 + t86) * t325 + (t433 + t436) * t306 - 0.2e1 * (-t236 * t298 + t238 * t299) * t305 + (t173 * t313 - t175 * t298 + t177 * t299 - 0.3e1 * t234 * t305 - t236 * t250 + t238 * t251) * t313) * t313 + (t6 + t13 - t240 * t467 + (t105 + t85) * t325 + 0.2e1 * (-t235 * t296 + t237 * t297) * t306 - (t433 + t437) * t305 + (t174 * t312 - t176 * t296 + t178 * t297 + 0.3e1 * t233 * t306 - t235 * t252 + t237 * t253) * t312 + (t173 * t312 + t174 * t313 - t175 * t296 - t176 * t298 + t177 * t297 + t178 * t299 - t235 * t250 - t236 * t252 + t237 * t251 + t238 * t253) * t313) * t312; (t111 * t35 + t112 * t34 + t57 * t92 + t58 * t91) * m(6) + (t102 * t155 + t103 * t154 + t156 * t73 + t157 * t72) * m(5) + t381 * t298 + t382 * t296 + t380 * t252 + t379 * t250 + t445 + t446; t59 * t410 + t23 * t409 + ((t154 * t421 + t155 * t422 + t347 * t73 - t349 * t72) * m(5) + (-t34 * t349 + t347 * t35 + t421 * t91 + t422 * t92) * m(6)) * t343; (t22 * t76 + t23 * t66 + t34 * t88 + t35 * t87 + t38 * t92 + t39 * t91) * m(6) + (t108 * t59 + t117 * t73 + t118 * t72 + t129 * t54 + t154 * t68 + t155 * t67) * m(5) - (t9 / 0.2e1 + t15 / 0.2e1) * t471 + (t26 / 0.2e1 + t52 / 0.2e1) * t306 - (t27 / 0.2e1 + t53 / 0.2e1) * t305 + (t5 / 0.2e1 + t12 / 0.2e1) * t298 + (t6 / 0.2e1 + t13 / 0.2e1) * t296 + (t28 / 0.2e1 + t55 / 0.2e1) * t252 + (t29 / 0.2e1 + t56 / 0.2e1) * t250 + (t106 * t312 + t107 * t313 + t130 * t325 + t37) * t456 + t469 * t455 + t470 * t454 + t468 * t453; (t23 * t76 + t34 * t92 + t35 * t91) * t396 + t304 * t36 + (t129 * t59 + t154 * t73 + t155 * t72) * t398 - (t128 + t468) * t471 + (t304 * t107 + t470) * t298 + (t304 * t106 + t469) * t296 + (t26 + t52) * t252 + (t27 + t53) * t250; (t109 * t58 + t110 * t57 + t111 * t49 + t112 * t48) * m(6) - t415 * t377 - t416 * t378 + t414 * t192 + t413 * t190 + t447; m(6) * (t32 * t442 + (t347 * t49 - t349 * t48 + (t109 * t349 + t110 * t347) * qJD(1)) * t343); t28 * t461 + t6 * t459 + t37 * t460 + t9 * t457 - t305 * t463 + t1 * t454 + t306 * t464 + t2 * t455 + t29 * t462 + t5 * t458 + (t109 * t39 + t110 * t38 + t22 * t96 + t32 * t66 + t48 * t88 + t49 * t87) * m(6) + t7 * t453; (t109 * t35 + t110 * t34 + t23 * t96 + t32 * t76 + t48 * t92 + t49 * t91) * m(6) + t26 * t461 + t4 * t459 + t27 * t462 + t3 * t458 + t250 * t463 + t298 * t1 / 0.2e1 + t252 * t464 + t296 * t2 / 0.2e1 + t33 * t456 - t471 * t7 / 0.2e1 + t36 * t460 + t8 * t457; (t109 * t49 + t110 * t48 + t32 * t96) * t396 + t190 * t25 - t377 * t1 + t192 * t24 - t378 * t2 + t261 * t33 - t376 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
