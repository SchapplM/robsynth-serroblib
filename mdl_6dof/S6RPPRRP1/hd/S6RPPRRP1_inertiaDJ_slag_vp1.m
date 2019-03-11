% Calculate time derivative of joint inertia matrix for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:37
% EndTime: 2019-03-09 01:57:54
% DurationCPUTime: 10.20s
% Computational Cost: add. (28850->710), mult. (27587->992), div. (0->0), fcn. (26095->10), ass. (0->351)
t239 = pkin(10) + qJ(4);
t234 = sin(t239);
t236 = cos(t239);
t245 = sin(qJ(5));
t247 = cos(qJ(5));
t284 = Icges(7,5) * t247 - Icges(7,6) * t245;
t179 = -Icges(7,3) * t236 + t234 * t284;
t285 = Icges(6,5) * t247 - Icges(6,6) * t245;
t180 = -Icges(6,3) * t236 + t234 * t285;
t454 = t179 + t180;
t394 = Icges(7,4) * t247;
t287 = -Icges(7,2) * t245 + t394;
t181 = -Icges(7,6) * t236 + t234 * t287;
t396 = Icges(6,4) * t247;
t288 = -Icges(6,2) * t245 + t396;
t182 = -Icges(6,6) * t236 + t234 * t288;
t452 = t182 + t181;
t395 = Icges(7,4) * t245;
t291 = Icges(7,1) * t247 - t395;
t183 = -Icges(7,5) * t236 + t234 * t291;
t397 = Icges(6,4) * t245;
t292 = Icges(6,1) * t247 - t397;
t184 = -Icges(6,5) * t236 + t234 * t292;
t445 = -t184 - t183;
t243 = -qJ(6) - pkin(8);
t453 = rSges(7,3) - t243;
t240 = qJ(1) + pkin(9);
t235 = sin(t240);
t237 = cos(t240);
t371 = t237 * t247;
t375 = t235 * t245;
t193 = -t236 * t375 - t371;
t372 = t237 * t245;
t374 = t235 * t247;
t194 = t236 * t374 - t372;
t378 = t234 * t235;
t130 = Icges(6,5) * t194 + Icges(6,6) * t193 + Icges(6,3) * t378;
t134 = Icges(6,4) * t194 + Icges(6,2) * t193 + Icges(6,6) * t378;
t138 = Icges(6,1) * t194 + Icges(6,4) * t193 + Icges(6,5) * t378;
t195 = -t236 * t372 + t374;
t196 = t236 * t371 + t375;
t377 = t234 * t237;
t55 = t130 * t377 + t134 * t195 + t138 * t196;
t131 = Icges(6,5) * t196 + Icges(6,6) * t195 + Icges(6,3) * t377;
t135 = Icges(6,4) * t196 + Icges(6,2) * t195 + Icges(6,6) * t377;
t139 = Icges(6,1) * t196 + Icges(6,4) * t195 + Icges(6,5) * t377;
t56 = t131 * t377 + t135 * t195 + t139 * t196;
t297 = t235 * t55 + t237 * t56;
t128 = Icges(7,5) * t194 + Icges(7,6) * t193 + Icges(7,3) * t378;
t132 = Icges(7,4) * t194 + Icges(7,2) * t193 + Icges(7,6) * t378;
t136 = Icges(7,1) * t194 + Icges(7,4) * t193 + Icges(7,5) * t378;
t53 = t128 * t377 + t132 * t195 + t136 * t196;
t129 = Icges(7,5) * t196 + Icges(7,6) * t195 + Icges(7,3) * t377;
t133 = Icges(7,4) * t196 + Icges(7,2) * t195 + Icges(7,6) * t377;
t137 = Icges(7,1) * t196 + Icges(7,4) * t195 + Icges(7,5) * t377;
t54 = t129 * t377 + t133 * t195 + t137 * t196;
t298 = t235 * t53 + t237 * t54;
t451 = t297 + t298;
t51 = t130 * t378 + t134 * t193 + t138 * t194;
t52 = t131 * t378 + t135 * t193 + t139 * t194;
t299 = t235 * t51 + t237 * t52;
t49 = t128 * t378 + t132 * t193 + t136 * t194;
t50 = t129 * t378 + t133 * t193 + t137 * t194;
t300 = t235 * t49 + t237 * t50;
t450 = t299 + t300;
t449 = -qJD(1) * t234 / 0.2e1;
t448 = t454 * t236 + (t245 * t452 + t247 * t445) * t234;
t348 = qJD(5) * t234;
t151 = (-Icges(7,5) * t245 - Icges(7,6) * t247) * t348 + (Icges(7,3) * t234 + t236 * t284) * qJD(4);
t152 = (-Icges(6,5) * t245 - Icges(6,6) * t247) * t348 + (Icges(6,3) * t234 + t236 * t285) * qJD(4);
t447 = -t152 - t151;
t153 = (-Icges(7,2) * t247 - t395) * t348 + (Icges(7,6) * t234 + t236 * t287) * qJD(4);
t154 = (-Icges(6,2) * t247 - t397) * t348 + (Icges(6,6) * t234 + t236 * t288) * qJD(4);
t443 = (-t153 - t154) * t245;
t279 = -t133 * t245 + t137 * t247;
t63 = -t129 * t236 + t234 * t279;
t277 = -t135 * t245 + t139 * t247;
t65 = -t131 * t236 + t234 * t277;
t407 = t63 + t65;
t280 = -t132 * t245 + t136 * t247;
t62 = -t128 * t236 + t234 * t280;
t278 = -t134 * t245 + t138 * t247;
t64 = -t130 * t236 + t234 * t278;
t408 = t62 + t64;
t442 = t235 * t408 + t237 * t407;
t315 = -qJD(5) * t236 + qJD(1);
t350 = qJD(4) * t245;
t252 = t234 * t350 + t247 * t315;
t356 = qJD(1) * t236;
t314 = -qJD(5) + t356;
t114 = t237 * t252 + t314 * t375;
t272 = t315 * t245;
t349 = qJD(4) * t247;
t251 = -t234 * t349 + t272;
t115 = t237 * t251 - t314 * t374;
t351 = qJD(4) * t237;
t327 = t236 * t351;
t357 = qJD(1) * t235;
t329 = t234 * t357;
t256 = t327 - t329;
t116 = t235 * t252 - t314 * t372;
t117 = t235 * t251 + t314 * t371;
t352 = qJD(4) * t236;
t322 = t235 * t352;
t355 = qJD(1) * t237;
t257 = t234 * t355 + t322;
t70 = Icges(7,5) * t117 + Icges(7,6) * t116 + Icges(7,3) * t257;
t74 = Icges(7,4) * t117 + Icges(7,2) * t116 + Icges(7,6) * t257;
t78 = Icges(7,1) * t117 + Icges(7,4) * t116 + Icges(7,5) * t257;
t11 = t114 * t132 + t115 * t136 + t128 * t256 + t195 * t74 + t196 * t78 + t377 * t70;
t69 = Icges(7,5) * t115 + Icges(7,6) * t114 + Icges(7,3) * t256;
t73 = Icges(7,4) * t115 + Icges(7,2) * t114 + Icges(7,6) * t256;
t77 = Icges(7,1) * t115 + Icges(7,4) * t114 + Icges(7,5) * t256;
t12 = t114 * t133 + t115 * t137 + t129 * t256 + t195 * t73 + t196 * t77 + t377 * t69;
t72 = Icges(6,5) * t117 + Icges(6,6) * t116 + Icges(6,3) * t257;
t76 = Icges(6,4) * t117 + Icges(6,2) * t116 + Icges(6,6) * t257;
t80 = Icges(6,1) * t117 + Icges(6,4) * t116 + Icges(6,5) * t257;
t13 = t114 * t134 + t115 * t138 + t130 * t256 + t195 * t76 + t196 * t80 + t377 * t72;
t71 = Icges(6,5) * t115 + Icges(6,6) * t114 + Icges(6,3) * t256;
t75 = Icges(6,4) * t115 + Icges(6,2) * t114 + Icges(6,6) * t256;
t79 = Icges(6,1) * t115 + Icges(6,4) * t114 + Icges(6,5) * t256;
t14 = t114 * t135 + t115 * t139 + t131 * t256 + t195 * t75 + t196 * t79 + t377 * t71;
t441 = (-t11 - t13) * t237 + (t12 + t14) * t235 + t451 * qJD(1);
t15 = t116 * t132 + t117 * t136 + t128 * t257 + t193 * t74 + t194 * t78 + t378 * t70;
t16 = t116 * t133 + t117 * t137 + t129 * t257 + t193 * t73 + t194 * t77 + t378 * t69;
t17 = t116 * t134 + t117 * t138 + t130 * t257 + t193 * t76 + t194 * t80 + t378 * t72;
t18 = t116 * t135 + t117 * t139 + t131 * t257 + t193 * t75 + t194 * t79 + t378 * t71;
t440 = (-t15 - t17) * t237 + (t16 + t18) * t235 + t450 * qJD(1);
t19 = (qJD(4) * t280 - t70) * t236 + (qJD(4) * t128 - t245 * t74 + t247 * t78 + (-t132 * t247 - t136 * t245) * qJD(5)) * t234;
t21 = (qJD(4) * t278 - t72) * t236 + (qJD(4) * t130 - t245 * t76 + t247 * t80 + (-t134 * t247 - t138 * t245) * qJD(5)) * t234;
t439 = -t19 - t21;
t20 = (qJD(4) * t279 - t69) * t236 + (qJD(4) * t129 - t245 * t73 + t247 * t77 + (-t133 * t247 - t137 * t245) * qJD(5)) * t234;
t22 = (qJD(4) * t277 - t71) * t236 + (qJD(4) * t131 - t245 * t75 + t247 * t79 + (-t135 * t247 - t139 * t245) * qJD(5)) * t234;
t438 = t20 + t22;
t86 = t179 * t378 + t181 * t193 + t183 * t194;
t87 = t180 * t378 + t182 * t193 + t184 * t194;
t437 = (-t86 - t87) * t236 + t450 * t234;
t88 = t179 * t377 + t181 * t195 + t183 * t196;
t89 = t180 * t377 + t182 * t195 + t184 * t196;
t436 = (-t88 - t89) * t236 + t451 * t234;
t425 = 2 * m(5);
t403 = rSges(5,1) * t236;
t308 = -rSges(5,2) * t234 + t403;
t174 = -rSges(5,3) * t237 + t235 * t308;
t373 = t236 * t237;
t229 = t235 * rSges(5,3);
t431 = -rSges(5,2) * t377 + t229;
t175 = rSges(5,1) * t373 + t431;
t208 = rSges(5,1) * t234 + rSges(5,2) * t236;
t263 = qJD(4) * t208;
t250 = rSges(5,2) * t329 + rSges(5,3) * t355 - t237 * t263;
t61 = (qJD(1) * t174 + t250) * t237 + (-t235 * t263 + (-t175 + t431) * qJD(1)) * t235;
t435 = t425 * t61;
t411 = pkin(8) + t243;
t434 = t236 * t411;
t399 = Icges(5,4) * t234;
t294 = Icges(5,1) * t236 - t399;
t173 = Icges(5,5) * t235 + t237 * t294;
t380 = t173 * t236;
t398 = Icges(5,4) * t236;
t290 = -Icges(5,2) * t234 + t398;
t171 = Icges(5,6) * t235 + t237 * t290;
t385 = t171 * t234;
t273 = -t380 + t385;
t433 = t273 * t237;
t225 = pkin(5) * t375;
t231 = pkin(5) * t247 + pkin(4);
t432 = t196 * rSges(7,1) + t195 * rSges(7,2) + t231 * t373 + t377 * t453 + t225;
t346 = qJD(5) * t247;
t340 = pkin(5) * t346;
t344 = pkin(5) * t372;
t345 = qJD(6) * t234;
t430 = t115 * rSges(7,1) + t114 * rSges(7,2) + rSges(7,3) * t327 + qJD(1) * t344 + t235 * t340 + t237 * t345 + t243 * t329;
t155 = (-Icges(7,1) * t245 - t394) * t348 + (Icges(7,5) * t234 + t236 * t291) * qJD(4);
t156 = (-Icges(6,1) * t245 - t396) * t348 + (Icges(6,5) * t234 + t236 * t292) * qJD(4);
t354 = qJD(4) * t234;
t429 = (t155 + t156) * t234 * t247 + t454 * t354 - t445 * t236 * t349;
t244 = -pkin(7) - qJ(3);
t242 = cos(pkin(10));
t230 = pkin(3) * t242 + pkin(2);
t264 = -t230 - t308;
t417 = sin(qJ(1)) * pkin(1);
t149 = -t417 + (rSges(5,3) - t244) * t237 + t264 * t235;
t238 = cos(qJ(1)) * pkin(1);
t311 = t237 * t230 - t235 * t244 + t238;
t150 = t175 + t311;
t428 = t149 * t237 + t150 * t235;
t286 = Icges(5,5) * t236 - Icges(5,6) * t234;
t168 = -Icges(5,3) * t237 + t235 * t286;
t170 = -Icges(5,6) * t237 + t235 * t290;
t172 = -Icges(5,5) * t237 + t235 * t294;
t427 = t236 * t407 - t436;
t221 = pkin(4) * t373;
t192 = pkin(8) * t377 + t221;
t367 = -t192 + t432;
t412 = pkin(4) - t231;
t254 = -t234 * t411 - t236 * t412;
t303 = -rSges(7,1) * t194 - rSges(7,2) * t193;
t368 = rSges(7,3) * t378 + t235 * t254 - t303 - t344;
t426 = -t235 * t368 - t237 * t367;
t270 = rSges(4,1) * t242 - rSges(4,2) * sin(pkin(10)) + pkin(2);
t401 = rSges(4,3) + qJ(3);
t162 = t235 * t401 + t237 * t270 + t238;
t424 = 2 * m(6);
t423 = 2 * m(7);
t232 = t235 ^ 2;
t233 = t237 ^ 2;
t421 = -t236 / 0.2e1;
t419 = -rSges(6,3) - pkin(8);
t418 = m(5) * t208;
t416 = pkin(4) * t234;
t415 = pkin(4) * t236;
t409 = t429 + (-t350 * t452 + t447) * t236 + ((t245 * t445 - t247 * t452) * qJD(5) + t443) * t234;
t214 = pkin(8) * t327;
t347 = qJD(5) * t245;
t341 = pkin(5) * t347;
t406 = -t214 + (pkin(8) * t357 + t351 * t412) * t234 + ((-qJD(4) * t243 - t341) * t237 + t412 * t357) * t236 - rSges(7,3) * t329 + t430;
t353 = qJD(4) * t235;
t213 = t353 * t416;
t304 = t117 * rSges(7,1) + t116 * rSges(7,2);
t379 = t231 * t234;
t405 = t213 + (qJD(1) * t254 - t340) * t237 + (t345 + pkin(5) * t272 + (-t379 - t434) * qJD(4)) * t235 + rSges(7,3) * t257 + t304;
t404 = t448 * t354;
t402 = rSges(7,3) * t234;
t306 = -rSges(6,1) * t194 - rSges(6,2) * t193;
t141 = rSges(6,3) * t378 - t306;
t390 = t141 * t237;
t387 = t170 * t234;
t386 = t170 * t236;
t384 = t171 * t236;
t383 = t172 * t234;
t382 = t172 * t236;
t381 = t173 * t234;
t143 = t196 * rSges(6,1) + t195 * rSges(6,2) + rSges(6,3) * t377;
t366 = -t143 - t192;
t302 = rSges(7,1) * t247 - rSges(7,2) * t245;
t324 = t234 * t347;
t365 = -pkin(5) * t324 - qJD(6) * t236 + (-rSges(7,1) * t245 - rSges(7,2) * t247) * t348 + (t236 * t302 + t254 + t402) * qJD(4);
t305 = rSges(6,1) * t247 - rSges(6,2) * t245;
t158 = (-rSges(6,1) * t245 - rSges(6,2) * t247) * t348 + (rSges(6,3) * t234 + t236 * t305) * qJD(4);
t312 = pkin(8) * t234 + t415;
t202 = t312 * qJD(4);
t364 = -t158 - t202;
t363 = -rSges(7,3) * t236 + t434 + (t302 - t412) * t234;
t191 = t312 * t235;
t362 = t235 * t191 + t237 * t192;
t186 = -rSges(6,3) * t236 + t234 * t305;
t209 = -pkin(8) * t236 + t416;
t361 = -t186 - t209;
t228 = qJD(3) * t237;
t360 = t244 * t357 + t228;
t359 = t232 + t233;
t169 = Icges(5,3) * t235 + t237 * t286;
t358 = qJD(1) * t169;
t343 = m(7) * t354;
t31 = t235 * t50 - t237 * t49;
t32 = t235 * t52 - t237 * t51;
t339 = t31 / 0.2e1 + t32 / 0.2e1;
t33 = t235 * t54 - t237 * t53;
t34 = t235 * t56 - t237 * t55;
t338 = t33 / 0.2e1 + t34 / 0.2e1;
t336 = t115 * rSges(6,1) + t114 * rSges(6,2) + rSges(6,3) * t327;
t328 = t234 * t351;
t333 = t235 * (pkin(8) * t257 + qJD(1) * t221 - t213) + t237 * (-pkin(8) * t329 + t214 + (-t235 * t356 - t328) * pkin(4)) + t191 * t355;
t332 = -t202 - t365;
t331 = -t209 - t363;
t330 = t186 * t357;
t321 = t235 * t363;
t320 = t237 * t363;
t319 = t368 * t237;
t160 = t361 * t237;
t318 = -t231 * t236 - t230;
t317 = qJD(1) * t363;
t316 = t236 * t341;
t109 = t331 * t237;
t310 = t235 * t317;
t307 = t117 * rSges(6,1) + t116 * rSges(6,2);
t301 = -t237 * t244 - t417;
t59 = t234 * t321 + t236 * t368;
t60 = -t234 * t320 - t236 * t367;
t296 = t235 * t60 + t237 * t59;
t255 = -t234 * t453 + t318;
t92 = -t417 + (pkin(5) * t245 - t244) * t237 + t255 * t235 + t303;
t93 = t311 + t432;
t295 = t235 * t93 + t237 * t92;
t289 = Icges(5,2) * t236 + t399;
t108 = t331 * t235;
t283 = t108 * t235 + t109 * t237;
t276 = -t143 * t235 + t390;
t275 = -t141 * t235 - t143 * t237;
t274 = -t382 + t387;
t271 = -t236 * t408 + t437;
t37 = t116 * t181 + t117 * t183 + t151 * t378 + t153 * t193 + t155 * t194 + t179 * t257;
t38 = t116 * t182 + t117 * t184 + t152 * t378 + t154 * t193 + t156 * t194 + t180 * t257;
t268 = t19 / 0.2e1 + t21 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1;
t35 = t114 * t181 + t115 * t183 + t151 * t377 + t153 * t195 + t155 * t196 + t179 * t256;
t36 = t114 * t182 + t115 * t184 + t152 * t377 + t154 * t195 + t156 * t196 + t180 * t256;
t267 = t20 / 0.2e1 + t22 / 0.2e1 + t35 / 0.2e1 + t36 / 0.2e1;
t266 = t63 / 0.2e1 + t65 / 0.2e1 + t88 / 0.2e1 + t89 / 0.2e1;
t265 = t87 / 0.2e1 + t62 / 0.2e1 + t64 / 0.2e1 + t86 / 0.2e1;
t261 = qJD(4) * t289;
t260 = qJD(4) * (-Icges(5,5) * t234 - Icges(5,6) * t236);
t259 = t234 * t419 - t230 - t415;
t253 = -t235 * t367 + t319;
t249 = t235 * t259 + t301;
t161 = -t235 * t270 + t237 * t401 - t417;
t227 = qJD(3) * t235;
t201 = t308 * qJD(4);
t197 = t209 * t357;
t159 = t361 * t235;
t145 = -qJD(1) * t162 + t228;
t144 = qJD(1) * t161 + t227;
t121 = t235 * t260 + t358;
t120 = -qJD(1) * t168 + t237 * t260;
t107 = t208 * t353 + (t237 * t264 - t229 - t238) * qJD(1) + t360;
t106 = t227 + ((-t230 - t403) * t235 + t301) * qJD(1) + t250;
t103 = -t143 * t236 - t186 * t377;
t102 = t141 * t236 + t186 * t378;
t101 = t311 - t366;
t100 = t249 + t306;
t97 = t169 * t235 - t433;
t96 = t168 * t235 - t237 * t274;
t95 = -t169 * t237 - t273 * t235;
t94 = -t168 * t237 - t235 * t274;
t91 = qJD(1) * t160 + t235 * t364;
t90 = t237 * t364 + t197 + t330;
t85 = t276 * t234;
t84 = rSges(6,3) * t257 + t307;
t82 = -rSges(6,3) * t329 + t336;
t66 = -t275 + t362;
t58 = qJD(1) * t109 + t235 * t332;
t57 = t237 * t332 + t197 + t310;
t48 = t213 + t419 * t322 + (t237 * t259 - t238) * qJD(1) - t307 + t360;
t47 = -pkin(4) * t328 + qJD(1) * t249 + t214 + t227 + t336;
t46 = t253 * t234;
t45 = t237 * t340 + (t316 - t345 + (-t236 * t453 + t379) * qJD(4)) * t235 + (t237 * t255 - t225 - t238) * qJD(1) - t304 + t360;
t44 = t227 + (-t316 + (-t236 * t243 - t379) * qJD(4)) * t237 + ((t318 - t402) * t235 + t301) * qJD(1) + t430;
t43 = t362 - t426;
t40 = (t186 * t353 + t84) * t236 + (-qJD(4) * t141 + t158 * t235 + t186 * t355) * t234;
t39 = (-t186 * t351 - t82) * t236 + (qJD(4) * t143 - t158 * t237 + t330) * t234;
t30 = t276 * t352 + (qJD(1) * t275 - t235 * t82 + t237 * t84) * t234;
t29 = t235 * t84 + t237 * t82 + (t235 * t366 + t390) * qJD(1) + t333;
t24 = (qJD(4) * t321 + t405) * t236 + (-qJD(4) * t368 + t235 * t365 + t237 * t317) * t234;
t23 = (-qJD(4) * t320 - t406) * t236 + (qJD(4) * t367 - t237 * t365 + t310) * t234;
t10 = t406 * t237 + t405 * t235 + (t319 + (-t192 - t367) * t235) * qJD(1) + t333;
t9 = t253 * t352 + (qJD(1) * t426 - t406 * t235 + t405 * t237) * t234;
t4 = (qJD(4) * t299 - t38) * t236 + (-qJD(1) * t32 + qJD(4) * t87 + t17 * t235 + t18 * t237) * t234;
t3 = (qJD(4) * t300 - t37) * t236 + (-qJD(1) * t31 + qJD(4) * t86 + t15 * t235 + t16 * t237) * t234;
t2 = (qJD(4) * t297 - t36) * t236 + (-qJD(1) * t34 + qJD(4) * t89 + t13 * t235 + t14 * t237) * t234;
t1 = (qJD(4) * t298 - t35) * t236 + (-qJD(1) * t33 + qJD(4) * t88 + t11 * t235 + t12 * t237) * t234;
t5 = [(t106 * t150 + t107 * t149) * t425 + 0.2e1 * m(4) * (t144 * t162 + t145 * t161) + (t44 * t93 + t45 * t92) * t423 + (t100 * t48 + t101 * t47) * t424 + t429 + (-t289 + t294) * t354 + (Icges(5,1) * t234 + t290 + t398) * t352 + t445 * t324 + t447 * t236 + t443 * t234 + t452 * (-t234 * t346 - t236 * t350); 0; 0; m(7) * (qJD(1) * t295 + t235 * t45 - t237 * t44) + m(6) * (t235 * t48 - t237 * t47 + (t100 * t237 + t101 * t235) * qJD(1)) + m(5) * (qJD(1) * t428 - t106 * t237 + t107 * t235) + m(4) * (-t144 * t237 + t145 * t235 + (t161 * t237 + t162 * t235) * qJD(1)); 0; 0; ((qJD(1) * t171 - t235 * t261) * t421 + t173 * t449 + (t387 / 0.2e1 - t382 / 0.2e1) * qJD(4) - t268) * t237 + ((-qJD(1) * t170 - t237 * t261) * t236 / 0.2e1 + t172 * t449 + (-t385 / 0.2e1 + t380 / 0.2e1) * qJD(4) + t267) * t235 + m(5) * ((-t106 * t235 - t107 * t237) * t208 - t428 * t201) + m(6) * (t100 * t90 + t101 * t91 + t159 * t47 + t160 * t48) + m(7) * (t108 * t44 + t109 * t45 + t57 * t92 + t58 * t93) + (t233 / 0.2e1 + t232 / 0.2e1) * t286 * qJD(4) + ((t384 / 0.2e1 + t381 / 0.2e1 - t150 * t418 + t266) * t237 + (t149 * t418 + t386 / 0.2e1 + t383 / 0.2e1 + t265) * t235) * qJD(1); m(5) * t61 + m(6) * t29 + m(7) * t10; m(6) * (t235 * t90 - t237 * t91 + (t159 * t235 + t160 * t237) * qJD(1)) + m(7) * (qJD(1) * t283 + t235 * t57 - t237 * t58); (t43 * t10 + t108 * t58 + t109 * t57) * t423 + (t159 * t91 + t160 * t90 + t66 * t29) * t424 + t359 * t208 * t201 * t425 + (t175 * t435 + (-t95 * qJD(1) + (-qJD(1) * t274 - t121) * t237) * t237 - t440) * t237 + (t174 * t435 + (t96 * qJD(1) + (t273 * qJD(1) + t120) * t235) * t235 + ((-t121 + (-t381 - t384) * qJD(4) + t171 * t352 + t173 * t354 - t358) * t235 + (t170 * t352 + t172 * t354 + t120 - (t383 + t386) * qJD(4)) * t237 + (t97 - t94 + (t169 - t274) * t235 + t433) * qJD(1)) * t237 + t441) * t235 + (t235 * t95 - t237 * t94 + t31 + t32) * t357 + (t235 * t97 - t237 * t96 + t33 + t34) * t355; m(6) * (t100 * t40 + t101 * t39 + t102 * t48 + t103 * t47) + m(7) * (t23 * t93 + t24 * t92 + t60 * t44 + t59 * t45) + ((t235 * t265 + t237 * t266) * qJD(4) - t409) * t236 + (t267 * t237 + t268 * t235 + (-t235 * t266 + t237 * t265) * qJD(1)) * t234 - t404; m(6) * t30 + m(7) * t9; m(6) * (t235 * t40 - t237 * t39 + (t102 * t237 + t103 * t235) * qJD(1)) + m(7) * (qJD(1) * t296 - t23 * t237 + t235 * t24); m(6) * (t102 * t90 + t103 * t91 + t159 * t39 + t160 * t40 + t29 * t85 + t30 * t66) + m(7) * (t46 * t10 + t108 * t23 + t109 * t24 + t9 * t43 + t59 * t57 + t60 * t58) + (-t4 / 0.2e1 - t3 / 0.2e1 + t338 * t352) * t237 + (t2 / 0.2e1 + t1 / 0.2e1 + t339 * t352) * t235 + ((-t235 * t338 + t237 * t339) * qJD(1) + t440 * t235 / 0.2e1 + t441 * t237 / 0.2e1 + (t235 * t407 - t237 * t408) * qJD(4) / 0.2e1) * t234 + (qJD(1) * t442 + t438 * t235 + t439 * t237) * t421 + (t437 * t235 + t436 * t237) * qJD(1) / 0.2e1; (t23 * t60 + t24 * t59 + t46 * t9) * t423 + (t102 * t40 + t103 * t39 + t30 * t85) * t424 + (t409 * t236 + (t271 * t235 - t237 * t427) * qJD(4) + t404) * t236 + ((-t236 * t438 + t1 + t2) * t237 + (t236 * t439 + t3 + t4) * t235 + (t442 * t234 + t236 * t448) * qJD(4) + (t235 * t427 + t271 * t237) * qJD(1)) * t234; m(7) * (t295 * t352 + (t235 * t44 + t237 * t45 + (-t235 * t92 + t237 * t93) * qJD(1)) * t234); t343; 0; m(7) * ((qJD(4) * t283 - t10) * t236 + (qJD(4) * t43 + t235 * t58 + t237 * t57 + (t108 * t237 - t109 * t235) * qJD(1)) * t234); m(7) * ((qJD(4) * t296 - t9) * t236 + (qJD(4) * t46 + t23 * t235 + t237 * t24 + (-t235 * t59 + t237 * t60) * qJD(1)) * t234); 0.2e1 * (-0.1e1 + t359) * t236 * t343;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
