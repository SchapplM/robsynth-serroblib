% Calculate time derivative of joint inertia matrix for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:43
% EndTime: 2019-03-09 02:02:59
% DurationCPUTime: 9.84s
% Computational Cost: add. (19996->700), mult. (27570->985), div. (0->0), fcn. (26498->8), ass. (0->340)
t248 = sin(qJ(4));
t251 = cos(qJ(4));
t247 = sin(qJ(5));
t250 = cos(qJ(5));
t289 = Icges(6,5) * t250 - Icges(6,6) * t247;
t190 = Icges(6,3) * t248 + t251 * t289;
t292 = Icges(7,4) * t250 + Icges(7,6) * t247;
t191 = Icges(7,2) * t248 + t251 * t292;
t453 = t191 + t190;
t246 = qJ(1) + pkin(9);
t243 = sin(t246);
t244 = cos(t246);
t371 = t247 * t248;
t185 = t243 * t371 - t244 * t250;
t369 = t248 * t250;
t186 = t243 * t369 + t244 * t247;
t373 = t243 * t251;
t121 = Icges(6,5) * t186 - Icges(6,6) * t185 - Icges(6,3) * t373;
t125 = Icges(6,4) * t186 - Icges(6,2) * t185 - Icges(6,6) * t373;
t129 = Icges(6,1) * t186 - Icges(6,4) * t185 - Icges(6,5) * t373;
t187 = t243 * t250 + t244 * t371;
t344 = t244 * t369;
t188 = t243 * t247 - t344;
t372 = t244 * t251;
t55 = t121 * t372 + t125 * t187 + t129 * t188;
t122 = Icges(6,5) * t188 + Icges(6,6) * t187 + Icges(6,3) * t372;
t126 = Icges(6,4) * t188 + Icges(6,2) * t187 + Icges(6,6) * t372;
t130 = Icges(6,1) * t188 + Icges(6,4) * t187 + Icges(6,5) * t372;
t56 = t122 * t372 + t126 * t187 + t130 * t188;
t302 = t243 * t55 - t244 * t56;
t119 = Icges(7,5) * t186 - Icges(7,6) * t373 + Icges(7,3) * t185;
t123 = Icges(7,4) * t186 - Icges(7,2) * t373 + Icges(7,6) * t185;
t127 = Icges(7,1) * t186 - Icges(7,4) * t373 + Icges(7,5) * t185;
t53 = -t119 * t187 + t123 * t372 + t127 * t188;
t120 = Icges(7,5) * t188 + Icges(7,6) * t372 - Icges(7,3) * t187;
t124 = Icges(7,4) * t188 + Icges(7,2) * t372 - Icges(7,6) * t187;
t128 = Icges(7,1) * t188 + Icges(7,4) * t372 - Icges(7,5) * t187;
t54 = -t120 * t187 + t124 * t372 + t128 * t188;
t303 = t243 * t53 - t244 * t54;
t452 = -t302 - t303;
t51 = -t121 * t373 - t125 * t185 + t129 * t186;
t52 = -t122 * t373 - t126 * t185 + t130 * t186;
t304 = t243 * t51 - t244 * t52;
t49 = t119 * t185 - t123 * t373 + t127 * t186;
t50 = t120 * t185 - t124 * t373 + t128 * t186;
t305 = t243 * t49 - t244 * t50;
t451 = -t304 - t305;
t355 = qJD(1) * t248;
t314 = qJD(5) + t355;
t277 = t247 * t314;
t350 = qJD(4) * t251;
t326 = t247 * t350;
t356 = qJD(1) * t244;
t115 = -qJD(5) * t344 + t243 * t277 - t244 * t326 - t250 * t356;
t276 = t314 * t250;
t315 = qJD(5) * t248 + qJD(1);
t424 = t247 * t315 - t250 * t350;
t116 = t243 * t276 + t244 * t424;
t395 = rSges(7,3) + qJ(6);
t446 = rSges(7,1) + pkin(5);
t450 = t187 * qJD(6) - t115 * t395 - t116 * t446;
t449 = -t187 * t395 + t188 * t446;
t388 = Icges(7,5) * t250;
t288 = Icges(7,3) * t247 + t388;
t189 = Icges(7,6) * t248 + t251 * t288;
t389 = Icges(7,5) * t247;
t296 = Icges(7,1) * t250 + t389;
t193 = Icges(7,4) * t248 + t251 * t296;
t88 = t185 * t189 + t186 * t193 - t191 * t373;
t391 = Icges(6,4) * t250;
t293 = -Icges(6,2) * t247 + t391;
t192 = Icges(6,6) * t248 + t251 * t293;
t392 = Icges(6,4) * t247;
t297 = Icges(6,1) * t250 - t392;
t194 = Icges(6,5) * t248 + t251 * t297;
t89 = -t185 * t192 + t186 * t194 - t190 * t373;
t445 = t451 * t251 + (t88 + t89) * t248;
t90 = -t187 * t189 + t188 * t193 + t191 * t372;
t91 = t187 * t192 + t188 * t194 + t190 * t372;
t444 = t452 * t251 + (t90 + t91) * t248;
t278 = t189 * t247 + t193 * t250;
t443 = (-t192 * t247 + t194 * t250 + t278) * t251 + t453 * t248;
t349 = qJD(5) * t251;
t155 = (Icges(7,3) * t250 - t389) * t349 + (Icges(7,6) * t251 - t248 * t288) * qJD(4);
t156 = (-Icges(6,5) * t247 - Icges(6,6) * t250) * t349 + (Icges(6,3) * t251 - t248 * t289) * qJD(4);
t157 = (-Icges(7,4) * t247 + Icges(7,6) * t250) * t349 + (Icges(7,2) * t251 - t248 * t292) * qJD(4);
t159 = (-Icges(7,1) * t247 + t388) * t349 + (Icges(7,4) * t251 - t248 * t296) * qJD(4);
t160 = (-Icges(6,1) * t247 - t391) * t349 + (Icges(6,5) * t251 - t248 * t297) * qJD(4);
t322 = t250 * t349;
t323 = t247 * t349;
t351 = qJD(4) * t248;
t325 = t250 * t351;
t327 = t247 * t351;
t370 = t247 * t251;
t442 = t155 * t370 + t189 * t322 + t192 * t327 - t193 * t323 - t194 * t325 + (t159 + t160) * t250 * t251 + t453 * t350 + (t156 + t157) * t248;
t328 = t244 * t351;
t354 = qJD(1) * t251;
t441 = t243 * t354 + t328;
t440 = t243 * t350 + t244 * t355;
t286 = t120 * t247 + t128 * t250;
t59 = t124 * t248 + t251 * t286;
t284 = t126 * t247 - t130 * t250;
t61 = t122 * t248 - t251 * t284;
t402 = t59 + t61;
t287 = t119 * t247 + t127 * t250;
t58 = t123 * t248 + t251 * t287;
t285 = t125 * t247 - t129 * t250;
t60 = t121 * t248 - t251 * t285;
t403 = t58 + t60;
t439 = -t243 * t403 + t244 * t402;
t117 = t244 * t277 + (t250 * t315 + t326) * t243;
t118 = -t243 * t424 + t244 * t276;
t330 = t243 * t351;
t331 = t244 * t354;
t261 = t330 - t331;
t68 = Icges(7,5) * t118 + Icges(7,6) * t261 + Icges(7,3) * t117;
t72 = Icges(7,4) * t118 + Icges(7,2) * t261 + Icges(7,6) * t117;
t76 = Icges(7,1) * t118 + Icges(7,4) * t261 + Icges(7,5) * t117;
t11 = t115 * t119 + t116 * t127 - t123 * t441 - t187 * t68 + t188 * t76 + t372 * t72;
t67 = Icges(7,5) * t116 - Icges(7,6) * t441 + Icges(7,3) * t115;
t71 = Icges(7,4) * t116 - Icges(7,2) * t441 + Icges(7,6) * t115;
t75 = Icges(7,1) * t116 - Icges(7,4) * t441 + Icges(7,5) * t115;
t12 = t115 * t120 + t116 * t128 - t124 * t441 - t187 * t67 + t188 * t75 + t372 * t71;
t70 = Icges(6,5) * t118 - Icges(6,6) * t117 + Icges(6,3) * t261;
t74 = Icges(6,4) * t118 - Icges(6,2) * t117 + Icges(6,6) * t261;
t78 = Icges(6,1) * t118 - Icges(6,4) * t117 + Icges(6,5) * t261;
t13 = -t115 * t125 + t116 * t129 - t121 * t441 + t187 * t74 + t188 * t78 + t372 * t70;
t69 = Icges(6,5) * t116 - Icges(6,6) * t115 - Icges(6,3) * t441;
t73 = Icges(6,4) * t116 - Icges(6,2) * t115 - Icges(6,6) * t441;
t77 = Icges(6,1) * t116 - Icges(6,4) * t115 - Icges(6,5) * t441;
t14 = -t115 * t126 + t116 * t130 - t122 * t441 + t187 * t73 + t188 * t77 + t372 * t69;
t438 = (t11 + t13) * t244 + (t12 + t14) * t243 + t452 * qJD(1);
t15 = t117 * t119 + t118 * t127 + t123 * t261 + t185 * t68 + t186 * t76 - t373 * t72;
t16 = t117 * t120 + t118 * t128 + t124 * t261 + t185 * t67 + t186 * t75 - t373 * t71;
t17 = -t117 * t125 + t118 * t129 + t121 * t261 - t185 * t74 + t186 * t78 - t373 * t70;
t18 = -t117 * t126 + t118 * t130 + t122 * t261 - t185 * t73 + t186 * t77 - t373 * t69;
t437 = (t15 + t17) * t244 + (t16 + t18) * t243 + t451 * qJD(1);
t19 = (-qJD(4) * t287 + t72) * t248 + (qJD(4) * t123 + t247 * t68 + t250 * t76 + (t119 * t250 - t127 * t247) * qJD(5)) * t251;
t21 = (qJD(4) * t285 + t70) * t248 + (qJD(4) * t121 - t247 * t74 + t250 * t78 + (-t125 * t250 - t129 * t247) * qJD(5)) * t251;
t436 = t19 + t21;
t20 = (-qJD(4) * t286 + t71) * t248 + (qJD(4) * t124 + t247 * t67 + t250 * t75 + (t120 * t250 - t128 * t247) * qJD(5)) * t251;
t22 = (qJD(4) * t284 + t69) * t248 + (qJD(4) * t122 - t247 * t73 + t250 * t77 + (-t126 * t250 - t130 * t247) * qJD(5)) * t251;
t435 = t20 + t22;
t421 = 2 * m(5);
t238 = t244 * rSges(5,3);
t374 = t243 * t248;
t175 = rSges(5,1) * t374 + rSges(5,2) * t373 + t238;
t310 = rSges(5,1) * t248 + rSges(5,2) * t251;
t267 = t244 * t310;
t396 = t243 * rSges(5,3);
t176 = t396 - t267;
t399 = rSges(5,2) * t248;
t226 = rSges(5,1) * t251 - t399;
t241 = t243 ^ 2;
t242 = t244 ^ 2;
t340 = rSges(5,1) * t440 + rSges(5,2) * t331;
t57 = -t243 * t340 + (-t226 * t242 + t241 * t399) * qJD(4) + ((-t175 + t238) * t244 + (-t176 + t267 + t396) * t243) * qJD(1);
t434 = t421 * t57;
t433 = -rSges(7,2) * t330 - qJD(6) * t185 - t117 * t395 - t118 * t446;
t393 = Icges(5,4) * t251;
t298 = Icges(5,1) * t248 + t393;
t173 = Icges(5,5) * t244 + t243 * t298;
t378 = t173 * t248;
t394 = Icges(5,4) * t248;
t294 = Icges(5,2) * t251 + t394;
t171 = Icges(5,6) * t244 + t243 * t294;
t381 = t171 * t251;
t280 = t378 + t381;
t432 = t243 * t280;
t365 = rSges(7,2) * t372 + t449;
t431 = t244 * t365;
t430 = t185 * t395 + t186 * t446;
t408 = pkin(4) * t248;
t231 = t244 * t408;
t237 = t244 * qJ(3);
t429 = t231 + t237;
t418 = -pkin(2) - pkin(7);
t347 = -rSges(5,3) + t418;
t409 = sin(qJ(1)) * pkin(1);
t265 = t243 * t347 - t409;
t359 = qJ(3) * t356 + qJD(3) * t243;
t100 = -rSges(5,2) * t330 + qJD(1) * t265 + t340 + t359;
t235 = qJD(3) * t244;
t245 = cos(qJ(1)) * pkin(1);
t352 = qJD(4) * t244;
t101 = t235 + t226 * t352 + (-t245 + t347 * t244 + (-qJ(3) - t310) * t243) * qJD(1);
t428 = -t100 * t244 + t101 * t243;
t290 = Icges(5,5) * t248 + Icges(5,6) * t251;
t427 = -Icges(5,3) * t243 + t244 * t290;
t426 = -Icges(5,6) * t243 + t244 * t294;
t425 = -Icges(5,5) * t243 + t244 * t298;
t423 = t248 * t403 + t445;
t367 = -rSges(7,2) * t373 + t430;
t422 = t243 * t365 + t244 * t367;
t420 = 2 * m(6);
t419 = 2 * m(7);
t415 = -t248 / 0.2e1;
t413 = t251 / 0.2e1;
t411 = rSges(4,2) - pkin(2);
t410 = m(5) * t226;
t401 = -rSges(7,2) * t441 - t450;
t400 = rSges(7,2) * t331 + t433;
t398 = rSges(7,2) * t251;
t397 = rSges(6,3) * t251;
t308 = -rSges(6,1) * t188 - rSges(6,2) * t187;
t134 = rSges(6,3) * t372 - t308;
t383 = t134 * t244;
t382 = t171 * t248;
t380 = t426 * t248;
t379 = t426 * t251;
t377 = t173 * t251;
t376 = t425 * t248;
t375 = t425 * t251;
t362 = rSges(6,1) * t186 - rSges(6,2) * t185;
t132 = -rSges(6,3) * t373 + t362;
t230 = pkin(4) * t374;
t197 = -pkin(8) * t373 + t230;
t366 = -t132 - t197;
t306 = rSges(7,1) * t250 + rSges(7,3) * t247;
t364 = (-pkin(5) * t351 + qJ(6) * t349) * t250 + (-qJ(6) * t351 + (-pkin(5) * qJD(5) + qJD(6)) * t251) * t247 + (-rSges(7,1) * t247 + rSges(7,3) * t250) * t349 + (-t248 * t306 + t398) * qJD(4);
t361 = rSges(7,2) * t248 + (pkin(5) * t250 + qJ(6) * t247 + t306) * t251;
t210 = (pkin(8) * t251 - t408) * qJD(4);
t229 = pkin(4) * t251 + pkin(8) * t248;
t360 = t210 * t243 + t229 * t356;
t169 = Icges(5,3) * t244 + t243 * t290;
t358 = qJD(1) * t169;
t357 = qJD(1) * t243;
t353 = qJD(4) * t243;
t31 = t243 * t50 + t244 * t49;
t32 = t243 * t52 + t244 * t51;
t346 = t31 / 0.2e1 + t32 / 0.2e1;
t33 = t243 * t54 + t244 * t53;
t34 = t243 * t56 + t244 * t55;
t345 = -t33 / 0.2e1 - t34 / 0.2e1;
t342 = rSges(6,1) * t118 - rSges(6,2) * t117 + rSges(6,3) * t330;
t341 = -t197 - t367;
t339 = pkin(4) * t440 + pkin(8) * t330;
t338 = t244 * pkin(4) * t350 + pkin(8) * t441;
t337 = t244 * pkin(2) + t243 * qJ(3) + t245;
t336 = (-rSges(7,2) - pkin(8)) * t251;
t335 = (-rSges(6,3) - pkin(8)) * t251;
t307 = rSges(6,1) * t250 - rSges(6,2) * t247;
t196 = rSges(6,3) * t248 + t251 * t307;
t334 = t196 * t357;
t321 = -qJ(3) - t408;
t207 = t310 * qJD(4);
t320 = t207 * (t241 + t242);
t319 = t361 * t243;
t318 = qJD(1) * t361;
t317 = qJD(4) * t361;
t158 = (-Icges(6,2) * t250 - t392) * t349 + (Icges(6,6) * t251 - t248 * t293) * qJD(4);
t316 = t443 * t350 + (-t278 * t351 + (-t247 * t158 + (-t192 * t250 - t194 * t247) * qJD(5)) * t251 + t442) * t248;
t313 = t235 + t338;
t312 = t244 * pkin(7) + t337;
t309 = t116 * rSges(6,1) - t115 * rSges(6,2);
t301 = t339 + t359;
t300 = t230 + t312;
t299 = Icges(5,1) * t251 - t394;
t295 = -Icges(5,2) * t248 + t393;
t291 = Icges(5,5) * t251 - Icges(5,6) * t248;
t283 = t132 * t244 + t134 * t243;
t279 = -t376 - t379;
t275 = -t248 * t402 - t444;
t35 = t115 * t189 + t116 * t193 - t155 * t187 + t157 * t372 + t159 * t188 - t191 * t441;
t36 = -t115 * t192 + t116 * t194 + t156 * t372 + t158 * t187 + t160 * t188 - t190 * t441;
t274 = t35 / 0.2e1 + t36 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1;
t37 = t117 * t189 + t118 * t193 + t155 * t185 - t157 * t373 + t159 * t186 + t191 * t261;
t38 = -t117 * t192 + t118 * t194 - t156 * t373 - t158 * t185 + t160 * t186 + t190 * t261;
t273 = -t37 / 0.2e1 - t38 / 0.2e1 - t21 / 0.2e1 - t19 / 0.2e1;
t272 = t60 / 0.2e1 + t58 / 0.2e1 + t88 / 0.2e1 + t89 / 0.2e1;
t271 = -t61 / 0.2e1 - t59 / 0.2e1 - t90 / 0.2e1 - t91 / 0.2e1;
t270 = t244 * t418 - t245;
t269 = t243 * t418 - t409;
t162 = (-rSges(6,1) * t247 - rSges(6,2) * t250) * t349 + (-t248 * t307 + t397) * qJD(4);
t268 = t162 * t243 + t196 * t356;
t266 = t279 * t243;
t264 = qJD(4) * t299;
t263 = qJD(4) * t295;
t259 = t322 - t327;
t258 = rSges(4,3) * t244 + t243 * t411 - t409;
t256 = t243 * t364 + t244 * t318;
t254 = t244 * t336 + t269;
t253 = t244 * t335 + t269;
t203 = t243 * t229;
t200 = t229 * t357;
t198 = pkin(8) * t372 - t231;
t177 = t244 * t198;
t166 = -rSges(4,2) * t244 + rSges(4,3) * t243 + t337;
t165 = t237 + t258;
t164 = (-t196 - t229) * t244;
t163 = t196 * t243 + t203;
t153 = -pkin(8) * t331 + t339;
t148 = t235 + (-t245 + t411 * t244 + (-rSges(4,3) - qJ(3)) * t243) * qJD(1);
t147 = qJD(1) * t258 + t359;
t143 = t244 * (qJD(1) * t230 - t338);
t142 = t312 + t175;
t141 = t237 + t267 + t265;
t136 = qJD(1) * t427 + t291 * t353;
t135 = -t291 * t352 + t358;
t109 = (-t229 - t361) * t244;
t108 = t203 + t319;
t103 = -t134 * t248 + t196 * t372;
t102 = t132 * t248 + t196 * t373;
t99 = t243 * t335 + t300 + t362;
t98 = t253 + t308 + t429;
t97 = -t243 * t427 - t244 * t279;
t96 = t169 * t243 - t244 * t280;
t95 = -t244 * t427 + t266;
t94 = t169 * t244 + t432;
t93 = t268 + t360;
t92 = t334 + t200 + (-t162 - t210) * t244;
t87 = t283 * t251;
t86 = t243 * t336 + t300 + t430;
t85 = t254 + t429 - t449;
t84 = -rSges(6,3) * t331 + t342;
t82 = -rSges(6,3) * t441 + t309;
t66 = -t248 * t365 + t361 * t372;
t65 = t248 * t367 + t251 * t319;
t64 = t243 * t366 + t177 + t383;
t63 = t256 + t360;
t62 = t200 + t243 * t318 + (-t210 - t364) * t244;
t48 = rSges(6,3) * t328 + ((t321 + t397) * t243 + t270) * qJD(1) - t309 + t313;
t47 = qJD(1) * t253 + t301 + t342;
t46 = t422 * t251;
t45 = t243 * t341 + t177 + t431;
t42 = (-t196 * t353 + t84) * t248 + (qJD(4) * t132 + t268) * t251;
t41 = (-t196 * t352 - t82) * t248 + (-qJD(4) * t134 + t162 * t244 - t334) * t251;
t40 = rSges(7,2) * t328 + ((t321 + t398) * t243 + t270) * qJD(1) + t313 + t450;
t39 = qJD(1) * t254 + t301 - t433;
t30 = t283 * t351 + (-t243 * t82 - t244 * t84 + (t132 * t243 - t383) * qJD(1)) * t251;
t25 = t244 * t82 + t143 + (-t153 - t84) * t243 + (t366 * t244 + (-t134 - t198) * t243) * qJD(1);
t24 = (-t243 * t317 - t400) * t248 + (qJD(4) * t367 + t256) * t251;
t23 = (-t244 * t317 - t401) * t248 + (-qJD(4) * t365 + t244 * t364 - t357 * t361) * t251;
t10 = t143 + t401 * t244 + (-t153 + t400) * t243 + (t341 * t244 + (-t198 - t365) * t243) * qJD(1);
t9 = t422 * t351 + (t400 * t244 - t401 * t243 + (t243 * t367 - t431) * qJD(1)) * t251;
t4 = (qJD(4) * t304 + t38) * t248 + (-qJD(1) * t32 + qJD(4) * t89 - t17 * t243 + t18 * t244) * t251;
t3 = (qJD(4) * t305 + t37) * t248 + (-qJD(1) * t31 + qJD(4) * t88 - t15 * t243 + t16 * t244) * t251;
t2 = (qJD(4) * t302 + t36) * t248 + (-qJD(1) * t34 + qJD(4) * t91 - t13 * t243 + t14 * t244) * t251;
t1 = (qJD(4) * t303 + t35) * t248 + (-qJD(1) * t33 + qJD(4) * t90 - t11 * t243 + t12 * t244) * t251;
t5 = [0.2e1 * m(4) * (t147 * t166 + t148 * t165) + (t100 * t142 + t101 * t141) * t421 + (t47 * t99 + t48 * t98) * t420 + (t39 * t86 + t40 * t85) * t419 - t193 * t325 - t298 * t350 - t194 * t323 - t192 * t322 - t248 * t264 - t251 * t263 - t189 * t327 + t294 * t351 - t158 * t370 + t442; 0; 0; m(7) * (t243 * t40 - t244 * t39 + (t243 * t86 + t244 * t85) * qJD(1)) + m(6) * (t243 * t48 - t244 * t47 + (t243 * t99 + t244 * t98) * qJD(1)) + m(5) * ((t141 * t244 + t142 * t243) * qJD(1) + t428) + m(4) * (-t147 * t244 + t148 * t243 + (t165 * t244 + t166 * t243) * qJD(1)); 0; 0; ((qJD(1) * t426 + t243 * t263) * t415 + (qJD(1) * t425 + t243 * t264) * t413 + (-t381 / 0.2e1 - t378 / 0.2e1) * qJD(4) - t273) * t244 + ((qJD(1) * t171 - t295 * t352) * t415 + (qJD(1) * t173 - t299 * t352) * t413 + (t379 / 0.2e1 + t376 / 0.2e1) * qJD(4) + t274) * t243 + m(5) * (t428 * t226 - (t141 * t243 - t142 * t244) * t207) + m(6) * (t163 * t48 + t164 * t47 + t92 * t99 + t93 * t98) + m(7) * (t108 * t40 + t109 * t39 + t62 * t86 + t63 * t85) - (t241 / 0.2e1 + t242 / 0.2e1) * t290 * qJD(4) + ((t382 / 0.2e1 - t377 / 0.2e1 + t142 * t410 - t272) * t243 + (t380 / 0.2e1 - t375 / 0.2e1 + t141 * t410 - t271) * t244) * qJD(1); m(5) * t57 + m(6) * t25 + m(7) * t10; m(6) * (t243 * t93 - t244 * t92 + (t163 * t244 + t164 * t243) * qJD(1)) + m(7) * (t243 * t63 - t244 * t62 + (t108 * t244 + t109 * t243) * qJD(1)) - m(5) * t320; (t10 * t45 + t108 * t63 + t109 * t62) * t419 + (t163 * t93 + t164 * t92 + t25 * t64) * t420 - t226 * t320 * t421 + (t176 * t434 + (t95 * qJD(1) + (t280 * qJD(1) + t136) * t244) * t244 + t437) * t244 + (-t175 * t434 + (t243 * t135 + (-t96 + t266) * qJD(1)) * t243 + ((t426 * t351 - t425 * t350 + t136 + (t375 - t380) * qJD(4)) * t243 + (t135 + (t377 - t382) * qJD(4) + t171 * t351 - t173 * t350 + t358) * t244 + (-t94 + t97 + t432 + (-t169 + t279) * t244) * qJD(1)) * t244 + t438) * t243 + (-t243 * t95 - t244 * t94 - t31 - t32) * t357 + (t243 * t97 + t244 * t96 + t33 + t34) * t356; m(6) * (t102 * t47 + t103 * t48 + t41 * t98 + t42 * t99) + m(7) * (t23 * t85 + t24 * t86 + t39 * t65 + t40 * t66) + (t243 * t272 + t244 * t271) * t351 + (t274 * t244 + t273 * t243 + (t243 * t271 - t244 * t272) * qJD(1)) * t251 + t316; m(6) * t30 + m(7) * t9; m(6) * (t243 * t41 - t244 * t42 + (t102 * t243 + t103 * t244) * qJD(1)) + m(7) * (t23 * t243 - t24 * t244 + (t243 * t65 + t244 * t66) * qJD(1)); m(6) * (t102 * t92 + t103 * t93 + t163 * t41 + t164 * t42 - t25 * t87 + t30 * t64) + m(7) * (-t10 * t46 + t108 * t23 + t109 * t24 + t45 * t9 + t62 * t65 + t63 * t66) + (t3 / 0.2e1 + t4 / 0.2e1 + t345 * t351 + t444 * qJD(1) / 0.2e1) * t244 + (t1 / 0.2e1 + t2 / 0.2e1 + t346 * t351 - t445 * qJD(1) / 0.2e1) * t243 + ((t243 * t345 - t244 * t346) * qJD(1) - t437 * t243 / 0.2e1 + t438 * t244 / 0.2e1 + (t243 * t402 + t244 * t403) * qJD(4) / 0.2e1) * t251 + (qJD(1) * t439 + t435 * t243 + t436 * t244) * t248 / 0.2e1; (t23 * t66 + t24 * t65 - t46 * t9) * t419 + (t102 * t42 + t103 * t41 - t30 * t87) * t420 + ((t243 * t423 + t275 * t244) * qJD(4) + t316) * t248 + ((t248 * t435 + t1 + t2) * t244 + (-t248 * t436 - t3 - t4) * t243 + (t248 * t443 + t251 * t439) * qJD(4) + (t275 * t243 - t244 * t423) * qJD(1)) * t251; m(7) * (t115 * t86 + t117 * t85 + t185 * t40 - t187 * t39); t259 * m(7); m(7) * (-t115 * t244 + t117 * t243 + (t185 * t244 - t187 * t243) * qJD(1)); m(7) * (t45 * t322 + t108 * t117 + t109 * t115 + t185 * t63 - t187 * t62 + (t10 * t251 - t351 * t45) * t247); m(7) * (-t46 * t322 + t115 * t65 + t117 * t66 + t185 * t23 - t187 * t24 + (t251 * t9 + t351 * t46) * t247); (-t115 * t187 + t117 * t185 + t259 * t370) * t419;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
