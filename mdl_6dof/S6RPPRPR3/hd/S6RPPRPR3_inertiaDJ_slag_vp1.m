% Calculate time derivative of joint inertia matrix for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:00
% EndTime: 2019-03-09 01:44:11
% DurationCPUTime: 9.13s
% Computational Cost: add. (15152->564), mult. (14969->799), div. (0->0), fcn. (13793->10), ass. (0->293)
t415 = Icges(5,3) + Icges(6,3);
t206 = qJ(4) + pkin(10);
t201 = sin(t206);
t203 = cos(t206);
t210 = sin(qJ(4));
t213 = cos(qJ(4));
t414 = -Icges(5,5) * t210 - Icges(6,5) * t201 - Icges(5,6) * t213 - Icges(6,6) * t203;
t207 = qJ(1) + pkin(9);
t204 = cos(t207);
t202 = sin(t207);
t351 = Icges(6,4) * t201;
t246 = Icges(6,2) * t203 + t351;
t114 = Icges(6,6) * t204 + t202 * t246;
t350 = Icges(6,4) * t203;
t251 = Icges(6,1) * t201 + t350;
t116 = Icges(6,5) * t204 + t202 * t251;
t353 = Icges(5,4) * t210;
t248 = Icges(5,2) * t213 + t353;
t127 = Icges(5,6) * t204 + t202 * t248;
t352 = Icges(5,4) * t213;
t253 = Icges(5,1) * t210 + t352;
t130 = Icges(5,5) * t204 + t202 * t253;
t389 = t114 * t203 + t116 * t201 + t127 * t213 + t130 * t210;
t413 = t204 * t389;
t412 = -t414 * t202 + t204 * t415;
t382 = -Icges(5,5) * t202 + t204 * t253;
t383 = -Icges(6,5) * t202 + t204 * t251;
t384 = -Icges(5,6) * t202 + t204 * t248;
t385 = -Icges(6,6) * t202 + t204 * t246;
t411 = -t201 * t383 - t203 * t385 - t210 * t382 - t213 * t384;
t410 = t389 * t202;
t408 = Icges(5,5) * t213 + Icges(6,5) * t203 - Icges(5,6) * t210 - Icges(6,6) * t201;
t407 = t202 * t415 + t414 * t204;
t406 = -t202 * t412 + t413;
t405 = t411 * t202;
t303 = qJD(4) * t213;
t312 = qJD(1) * t204;
t404 = t202 * t303 + t210 * t312;
t308 = qJD(4) * t203;
t403 = -t201 * t312 - t202 * t308;
t370 = rSges(7,3) + pkin(8);
t291 = t370 * t203;
t365 = pkin(5) * t201;
t402 = -t291 + t365;
t401 = t412 * qJD(1);
t209 = sin(qJ(6));
t212 = cos(qJ(6));
t348 = Icges(7,4) * t212;
t245 = -Icges(7,2) * t209 + t348;
t126 = Icges(7,6) * t201 + t203 * t245;
t349 = Icges(7,4) * t209;
t250 = Icges(7,1) * t212 - t349;
t129 = Icges(7,5) * t201 + t203 * t250;
t400 = -t126 * t212 - t129 * t209;
t200 = t204 ^ 2;
t399 = -t204 * t412 - t410;
t398 = -t204 * t407 - t405;
t397 = t202 * t407 - t204 * t411;
t307 = qJD(4) * t204;
t396 = -t307 * t408 + t401;
t309 = qJD(4) * t202;
t395 = -qJD(1) * t407 + t309 * t408;
t328 = t202 * t210;
t186 = pkin(4) * t328;
t208 = -qJ(5) - pkin(7);
t311 = qJD(1) * t208;
t271 = qJD(5) * t202 + (-pkin(4) * t303 - t311) * t204;
t363 = pkin(7) * t204;
t101 = t204 * ((t186 - t363) * qJD(1) + t271);
t275 = pkin(4) * t404 + qJD(5) * t204 + t202 * t311;
t313 = qJD(1) * t202;
t105 = pkin(7) * t313 + t275;
t268 = rSges(6,1) * t201 + rSges(6,2) * t203;
t228 = t204 * t268;
t358 = t202 * rSges(6,3);
t119 = t358 - t228;
t367 = pkin(4) * t210;
t317 = t202 * t208 + t204 * t367;
t135 = -pkin(7) * t202 - t317;
t360 = rSges(6,2) * t201;
t153 = rSges(6,1) * t203 - t360;
t196 = t204 * rSges(6,3);
t199 = t202 ^ 2;
t289 = t203 * t312;
t295 = rSges(6,1) * t403 - rSges(6,2) * t289;
t330 = t202 * t203;
t331 = t201 * t202;
t118 = rSges(6,1) * t331 + rSges(6,2) * t330 + t196;
t136 = t186 + (-pkin(7) - t208) * t204;
t321 = -t118 - t136;
t19 = t101 + (-t105 + t295) * t202 + (-t153 * t200 + t199 * t360) * qJD(4) + ((t321 + t196) * t204 + (-t119 - t135 + t228 + t358) * t202) * qJD(1);
t378 = 2 * m(6);
t394 = t19 * t378;
t197 = t204 * rSges(5,3);
t326 = t202 * t213;
t132 = rSges(5,1) * t328 + rSges(5,2) * t326 + t197;
t269 = rSges(5,1) * t210 + rSges(5,2) * t213;
t229 = t204 * t269;
t359 = t202 * rSges(5,3);
t134 = t359 - t229;
t361 = rSges(5,2) * t210;
t183 = rSges(5,1) * t213 - t361;
t287 = t213 * t312;
t293 = rSges(5,1) * t404 + rSges(5,2) * t287;
t32 = -t202 * t293 + (-t183 * t200 + t199 * t361) * qJD(4) + ((-t132 + t197) * t204 + (-t134 + t229 + t359) * t202) * qJD(1);
t379 = 2 * m(5);
t393 = t32 * t379;
t195 = t204 * qJ(3);
t392 = t195 + t317;
t205 = cos(qJ(1)) * pkin(1);
t391 = -t204 * pkin(2) - t205;
t301 = -rSges(5,3) - pkin(2) - pkin(7);
t368 = sin(qJ(1)) * pkin(1);
t224 = t202 * t301 - t368;
t305 = qJD(4) * t210;
t316 = qJ(3) * t312 + qJD(3) * t202;
t62 = -rSges(5,2) * t202 * t305 + qJD(1) * t224 + t293 + t316;
t193 = qJD(3) * t204;
t63 = t193 + t183 * t307 + (-t205 + t301 * t204 + (-qJ(3) - t269) * t202) * qJD(1);
t388 = t202 * t63 - t204 * t62;
t277 = qJD(6) * t201 + qJD(1);
t304 = qJD(4) * t212;
t381 = -t203 * t304 + t209 * t277;
t306 = qJD(4) * t209;
t380 = t203 * t306 + t212 * t277;
t377 = 2 * m(7);
t376 = t201 / 0.2e1;
t375 = -t202 / 0.2e1;
t373 = t204 / 0.2e1;
t372 = rSges(4,2) - pkin(2);
t371 = -rSges(6,3) - pkin(2);
t369 = m(5) * t183;
t366 = pkin(4) * t213;
t364 = pkin(5) * t203;
t240 = Icges(7,5) * t212 - Icges(7,6) * t209;
t123 = Icges(7,3) * t201 + t203 * t240;
t302 = qJD(6) * t203;
t92 = (-Icges(7,5) * t209 - Icges(7,6) * t212) * t302 + (Icges(7,3) * t203 - t201 * t240) * qJD(4);
t98 = (-Icges(7,1) * t209 - t348) * t302 + (Icges(7,5) * t203 - t201 * t250) * qJD(4);
t215 = t203 * t212 * t98 + t123 * t308 + (t126 * t306 - t129 * t304 + t92) * t201;
t95 = (-Icges(7,2) * t212 - t349) * t302 + (Icges(7,6) * t203 - t201 * t245) * qJD(4);
t355 = t209 * t95;
t59 = t123 * t201 + (-t126 * t209 + t129 * t212) * t203;
t362 = ((qJD(6) * t400 - t355) * t203 + t215) * t201 + t59 * t308;
t272 = pkin(8) * t203 - t365;
t324 = t204 * t209;
t327 = t202 * t212;
t143 = t201 * t324 + t327;
t323 = t204 * t212;
t329 = t202 * t209;
t144 = -t201 * t323 + t329;
t267 = -rSges(7,1) * t144 - rSges(7,2) * t143;
t325 = t203 * t204;
t90 = rSges(7,3) * t325 - t267;
t354 = t272 * t204 + t90;
t341 = t114 * t201;
t340 = t201 * t385;
t339 = t116 * t203;
t338 = t203 * t383;
t336 = t127 * t210;
t335 = t210 * t384;
t333 = t130 * t213;
t332 = t382 * t213;
t266 = rSges(7,1) * t212 - rSges(7,2) * t209;
t102 = (-rSges(7,1) * t209 - rSges(7,2) * t212) * t302 + (rSges(7,3) * t203 - t201 * t266) * qJD(4);
t322 = -t272 * qJD(4) - t102;
t133 = rSges(7,3) * t201 + t203 * t266;
t154 = pkin(8) * t201 + t364;
t320 = t133 + t154;
t141 = -t201 * t329 + t323;
t142 = t201 * t327 + t324;
t319 = t142 * rSges(7,1) + t141 * rSges(7,2);
t187 = pkin(4) * t326;
t299 = pkin(4) * t305;
t318 = qJD(1) * t187 + t204 * t299;
t310 = qJD(4) * t201;
t168 = pkin(5) * t331;
t286 = t201 * t309;
t276 = qJD(1) * t201 + qJD(6);
t234 = t276 * t209;
t72 = -t202 * t380 - t204 * t234;
t235 = t212 * t276;
t73 = -t202 * t381 + t204 * t235;
t300 = t73 * rSges(7,1) + t72 * rSges(7,2) + rSges(7,3) * t286;
t85 = Icges(7,4) * t142 + Icges(7,2) * t141 - Icges(7,6) * t330;
t87 = Icges(7,1) * t142 + Icges(7,4) * t141 - Icges(7,5) * t330;
t258 = t209 * t85 - t212 * t87;
t83 = Icges(7,5) * t142 + Icges(7,6) * t141 - Icges(7,3) * t330;
t30 = t201 * t83 - t203 * t258;
t44 = -t123 * t330 + t126 * t141 + t129 * t142;
t298 = t44 / 0.2e1 + t30 / 0.2e1;
t86 = Icges(7,4) * t144 + Icges(7,2) * t143 + Icges(7,6) * t325;
t88 = Icges(7,1) * t144 + Icges(7,4) * t143 + Icges(7,5) * t325;
t257 = t209 * t86 - t212 * t88;
t84 = Icges(7,5) * t144 + Icges(7,6) * t143 + Icges(7,3) * t325;
t31 = t201 * t84 - t203 * t257;
t45 = t123 * t325 + t126 * t143 + t129 * t144;
t297 = -t45 / 0.2e1 - t31 / 0.2e1;
t89 = -rSges(7,3) * t330 + t319;
t296 = pkin(8) * t330 - t136 - t168 - t89;
t294 = pkin(5) * t403 - pkin(8) * t286;
t292 = t202 * qJ(3) - t391;
t285 = t201 * t307;
t280 = -qJ(3) - t367;
t165 = t269 * qJD(4);
t279 = (t199 + t200) * t165;
t278 = qJD(1) * t320;
t70 = -t202 * t234 + t204 * t380;
t71 = t202 * t235 + t204 * t381;
t274 = t71 * rSges(7,1) + t70 * rSges(7,2);
t273 = -pkin(2) * t202 - t368;
t26 = t141 * t85 + t142 * t87 - t330 * t83;
t27 = t141 * t86 + t142 * t88 - t330 * t84;
t15 = t202 * t27 + t204 * t26;
t263 = t202 * t26 - t204 * t27;
t28 = t143 * t85 + t144 * t87 + t325 * t83;
t29 = t143 * t86 + t144 * t88 + t325 * t84;
t16 = t202 * t29 + t204 * t28;
t262 = t202 * t28 - t204 * t29;
t261 = t202 * t31 + t204 * t30;
t260 = t202 * t30 - t204 * t31;
t259 = t202 * t90 + t204 * t89;
t254 = Icges(5,1) * t213 - t353;
t252 = Icges(6,1) * t203 - t351;
t249 = -Icges(5,2) * t210 + t352;
t247 = -Icges(6,2) * t201 + t350;
t233 = t193 - t271;
t231 = t202 * t371 - t368;
t230 = t275 + t316;
t227 = -t204 * t208 + t186 + t292;
t223 = qJD(4) * t254;
t222 = qJD(4) * t252;
t221 = qJD(4) * t249;
t220 = qJD(4) * t247;
t218 = t286 - t289;
t217 = -t203 * t313 - t285;
t216 = rSges(4,3) * t204 + t202 * t372 - t368;
t177 = pkin(4) * t287;
t148 = t268 * qJD(4);
t122 = (-t153 - t366) * t204;
t121 = t153 * t202 + t187;
t120 = t204 * t135;
t111 = -rSges(4,2) * t204 + rSges(4,3) * t202 + t292;
t110 = t195 + t216;
t107 = t193 + (-t205 + t372 * t204 + (-rSges(4,3) - qJ(3)) * t202) * qJD(1);
t106 = qJD(1) * t216 + t316;
t104 = t132 + t292 + t363;
t103 = t195 + t229 + t224;
t81 = (-t320 - t366) * t204;
t80 = t202 * t320 + t187;
t69 = t153 * t312 + t177 + (-t148 - t299) * t202;
t68 = t148 * t204 + t153 * t313 + t318;
t65 = t227 + t118;
t64 = t228 + t231 + t392;
t56 = t133 * t325 - t201 * t90;
t55 = t133 * t330 + t201 * t89;
t53 = -t202 * t291 + t168 + t227 + t319;
t52 = t204 * t402 + t267 + t273 + t392;
t47 = t153 * t307 + (-t205 + t371 * t204 + (-t268 + t280) * t202) * qJD(1) + t233;
t46 = -rSges(6,2) * t286 + qJD(1) * t231 + t230 - t295;
t43 = t259 * t203;
t42 = t177 + t204 * t278 + (-t299 - t322) * t202;
t41 = t202 * t278 + t204 * t322 + t318;
t40 = -rSges(7,3) * t289 + t300;
t39 = rSges(7,3) * t217 + t274;
t38 = Icges(7,1) * t73 + Icges(7,4) * t72 + Icges(7,5) * t218;
t37 = Icges(7,1) * t71 + Icges(7,4) * t70 + Icges(7,5) * t217;
t36 = Icges(7,4) * t73 + Icges(7,2) * t72 + Icges(7,6) * t218;
t35 = Icges(7,4) * t71 + Icges(7,2) * t70 + Icges(7,6) * t217;
t34 = Icges(7,5) * t73 + Icges(7,6) * t72 + Icges(7,3) * t218;
t33 = Icges(7,5) * t71 + Icges(7,6) * t70 + Icges(7,3) * t217;
t25 = t202 * t296 + t204 * t354 + t120;
t24 = (t201 * t370 + t364) * t307 + ((t280 - t402) * t202 + t391) * qJD(1) + t233 - t274;
t23 = (-t204 * t291 + t273) * qJD(1) + t230 - t294 + t300;
t21 = (-t133 * t309 + t40) * t201 + (qJD(4) * t89 + t102 * t202 + t133 * t312) * t203;
t20 = (-t133 * t307 - t39) * t201 + (-qJD(4) * t90 + t102 * t204 - t133 * t313) * t203;
t18 = t123 * t218 + t126 * t72 + t129 * t73 + t141 * t95 + t142 * t98 - t330 * t92;
t17 = t123 * t217 + t126 * t70 + t129 * t71 + t143 * t95 + t144 * t98 + t325 * t92;
t14 = t259 * t310 + (-t202 * t39 - t204 * t40 + (t202 * t89 - t204 * t90) * qJD(1)) * t203;
t13 = t201 * t45 - t203 * t262;
t12 = t201 * t44 - t203 * t263;
t11 = t101 + (-t105 - t40 + (-t135 - t354) * qJD(1) + t294) * t202 + (t39 - t154 * t307 + (t296 + t168) * qJD(1)) * t204;
t10 = (qJD(4) * t257 + t33) * t201 + (qJD(4) * t84 - t209 * t35 + t212 * t37 + (-t209 * t88 - t212 * t86) * qJD(6)) * t203;
t9 = (qJD(4) * t258 + t34) * t201 + (qJD(4) * t83 - t209 * t36 + t212 * t38 + (-t209 * t87 - t212 * t85) * qJD(6)) * t203;
t8 = -t84 * t289 + t141 * t35 + t142 * t37 + t72 * t86 + t73 * t88 + (-t203 * t33 + t310 * t84) * t202;
t7 = -t83 * t289 + t141 * t36 + t142 * t38 + t72 * t85 + t73 * t87 + (-t203 * t34 + t310 * t83) * t202;
t6 = -t84 * t285 + t143 * t35 + t144 * t37 + t70 * t86 + t71 * t88 + (t204 * t33 - t313 * t84) * t203;
t5 = -t83 * t285 + t143 * t36 + t144 * t38 + t70 * t85 + t71 * t87 + (t204 * t34 - t313 * t83) * t203;
t4 = -qJD(1) * t263 + t202 * t8 + t204 * t7;
t3 = -qJD(1) * t262 + t202 * t6 + t204 * t5;
t2 = (qJD(4) * t263 + t18) * t201 + (-qJD(1) * t15 + qJD(4) * t44 - t202 * t7 + t204 * t8) * t203;
t1 = (qJD(4) * t262 + t17) * t201 + (-qJD(1) * t16 + qJD(4) * t45 - t202 * t5 + t204 * t6) * t203;
t22 = [(t23 * t53 + t24 * t52) * t377 + (t46 * t65 + t47 * t64) * t378 + (t103 * t63 + t104 * t62) * t379 + 0.2e1 * m(4) * (t106 * t111 + t107 * t110) + t215 - t210 * t223 - t213 * t221 - t253 * t303 + t248 * t305 - t201 * t222 - t251 * t308 + t246 * t310 + t400 * t302 + (-t220 - t355) * t203; 0; 0; m(7) * (t202 * t24 - t204 * t23 + (t202 * t53 + t204 * t52) * qJD(1)) + m(6) * (t202 * t47 - t204 * t46 + (t202 * t65 + t204 * t64) * qJD(1)) + m(5) * ((t103 * t204 + t104 * t202) * qJD(1) + t388) + m(4) * (-t106 * t204 + t107 * t202 + (t110 * t204 + t111 * t202) * qJD(1)); 0; 0; m(5) * (t388 * t183 - (t103 * t202 - t104 * t204) * t165) + m(6) * (t121 * t47 + t122 * t46 + t64 * t69 + t65 * t68) + m(7) * (t23 * t81 + t24 * t80 + t41 * t53 + t42 * t52) + ((t336 / 0.2e1 - t333 / 0.2e1 + t104 * t369 + t341 / 0.2e1 - t339 / 0.2e1 - t298) * t202 + (t335 / 0.2e1 - t332 / 0.2e1 + t340 / 0.2e1 - t338 / 0.2e1 + t103 * t369 - t297) * t204) * qJD(1) + t414 * qJD(4) * (t199 / 0.2e1 + t200 / 0.2e1) + (-t411 * qJD(4) - t201 * (qJD(1) * t114 - t247 * t307) + t203 * (qJD(1) * t116 - t252 * t307) - t210 * (qJD(1) * t127 - t249 * t307) + t213 * (qJD(1) * t130 - t254 * t307) + t10 + t17) * t202 / 0.2e1 + (-qJD(4) * t389 + (qJD(1) * t382 + t202 * t223) * t213 - t201 * (qJD(1) * t385 + t202 * t220) + t203 * (qJD(1) * t383 + t202 * t222) - t210 * (qJD(1) * t384 + t202 * t221) + t18 + t9) * t373; m(5) * t32 + m(6) * t19 + m(7) * t11; m(6) * (t202 * t69 - t204 * t68 + (t121 * t204 + t122 * t202) * qJD(1)) + m(7) * (t202 * t42 - t204 * t41 + (t202 * t81 + t204 * t80) * qJD(1)) - m(5) * t279; (t11 * t25 + t41 * t81 + t42 * t80) * t377 + t16 * t312 - t15 * t313 - t183 * t279 * t379 + (t120 * t19 + t121 * t69 + t122 * t68) * t378 + (t119 * t394 + t134 * t393 + t395 * t200 + t399 * t313 + t4 + (-t398 - t406 + t413) * t312) * t204 + (t3 - t132 * t393 + t321 * t394 + (t396 * t202 + (t405 + t406) * qJD(1)) * t202 + t398 * t313 + t397 * t312 + ((-t303 * t382 + t305 * t384 - t308 * t383 + t310 * t385 + t395) * t202 + (t114 * t310 - t116 * t308 + t127 * t305 - t130 * t303 + t396 + t401) * t204 + ((t338 - t340 + t332 - t335) * t202 + (t339 - t341 + t333 - t336) * t204) * qJD(4) + (t410 + (t411 - t412) * t204 + t397 + t399) * qJD(1)) * t204) * t202; m(7) * (t202 * t23 + t204 * t24 + (-t202 * t52 + t204 * t53) * qJD(1)) + m(6) * (t202 * t46 + t204 * t47 + (-t202 * t64 + t204 * t65) * qJD(1)); 0; 0; m(7) * (t202 * t41 + t204 * t42 + (-t202 * t80 + t204 * t81) * qJD(1)) + m(6) * (t202 * t68 + t204 * t69 + (-t121 * t202 + t122 * t204) * qJD(1)); 0; m(7) * (t20 * t52 + t21 * t53 + t23 * t55 + t24 * t56) + (t202 * t298 + t204 * t297) * t310 + ((t17 / 0.2e1 + t10 / 0.2e1) * t204 + (-t18 / 0.2e1 - t9 / 0.2e1) * t202 + (t202 * t297 - t204 * t298) * qJD(1)) * t203 + t362; m(7) * t14; m(7) * (t20 * t202 - t204 * t21 + (t202 * t55 + t204 * t56) * qJD(1)); m(7) * (-t11 * t43 + t14 * t25 + t20 * t80 + t21 * t81 + t41 * t55 + t42 * t56) + (-t16 * t310 / 0.2e1 + (qJD(1) * t31 + t9) * t376 + qJD(1) * t13 / 0.2e1 + t2 / 0.2e1) * t204 + ((-qJD(1) * t30 + t10) * t376 + t1 / 0.2e1 - qJD(1) * t12 / 0.2e1 + t15 * t310 / 0.2e1) * t202 + (t3 * t373 + qJD(4) * t261 / 0.2e1 + t4 * t375 + (t16 * t375 - t204 * t15 / 0.2e1) * qJD(1)) * t203; m(7) * (t20 * t204 + t202 * t21 + (-t202 * t56 + t204 * t55) * qJD(1)); (-t14 * t43 + t20 * t56 + t21 * t55) * t377 + ((t202 * t12 - t204 * t13 + t201 * t260) * qJD(4) + t362) * t201 + (-t202 * t2 + t204 * t1 + t201 * (t10 * t204 - t202 * t9) + (t59 * t201 - t203 * t260) * qJD(4) + (-t204 * t12 - t202 * t13 - t201 * t261) * qJD(1)) * t203;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t22(1) t22(2) t22(4) t22(7) t22(11) t22(16); t22(2) t22(3) t22(5) t22(8) t22(12) t22(17); t22(4) t22(5) t22(6) t22(9) t22(13) t22(18); t22(7) t22(8) t22(9) t22(10) t22(14) t22(19); t22(11) t22(12) t22(13) t22(14) t22(15) t22(20); t22(16) t22(17) t22(18) t22(19) t22(20) t22(21);];
Mq  = res;
