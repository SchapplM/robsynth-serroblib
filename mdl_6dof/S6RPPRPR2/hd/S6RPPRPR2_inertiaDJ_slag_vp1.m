% Calculate time derivative of joint inertia matrix for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:34
% EndTime: 2019-03-09 01:41:48
% DurationCPUTime: 9.60s
% Computational Cost: add. (16594->555), mult. (15822->788), div. (0->0), fcn. (14569->10), ass. (0->283)
t400 = Icges(6,1) + Icges(5,3);
t200 = pkin(10) + qJ(4);
t197 = cos(t200);
t195 = sin(t200);
t330 = Icges(6,6) * t195;
t338 = Icges(5,4) * t195;
t399 = -t330 - t338 + (-Icges(5,2) - Icges(6,3)) * t197;
t329 = Icges(6,6) * t197;
t337 = Icges(5,4) * t197;
t398 = -t329 - t337 + (-Icges(5,1) - Icges(6,2)) * t195;
t240 = Icges(5,5) * t197 - Icges(5,6) * t195;
t242 = Icges(6,4) * t197 - Icges(6,5) * t195;
t397 = t240 - t242;
t201 = qJ(1) + pkin(9);
t198 = cos(t201);
t196 = sin(t201);
t247 = Icges(5,1) * t197 - t338;
t122 = Icges(5,5) * t196 + t198 * t247;
t238 = Icges(6,2) * t197 - t330;
t125 = Icges(6,4) * t196 - t198 * t238;
t370 = t122 - t125;
t244 = -Icges(5,2) * t195 + t337;
t120 = Icges(5,6) * t196 + t198 * t244;
t236 = -Icges(6,3) * t195 + t329;
t123 = Icges(6,5) * t196 - t198 * t236;
t372 = t120 - t123;
t389 = -t195 * t372 + t197 * t370;
t396 = t198 * t389;
t119 = -Icges(5,6) * t198 + t196 * t244;
t364 = Icges(6,5) * t198 + t196 * t236;
t373 = t119 + t364;
t121 = -Icges(5,5) * t198 + t196 * t247;
t363 = Icges(6,4) * t198 + t196 * t238;
t371 = t121 + t363;
t391 = t196 * t400 + t397 * t198;
t395 = t399 * qJD(4);
t394 = t398 * qJD(4);
t367 = -t195 * t373 + t197 * t371;
t393 = t198 * t367;
t392 = t397 * t196 - t198 * t400;
t390 = ((Icges(6,5) - Icges(5,6)) * t197 + (Icges(6,4) - Icges(5,5)) * t195) * qJD(4);
t301 = qJD(4) * t198;
t282 = t195 * t301;
t306 = qJD(1) * t196;
t388 = t197 * t306 + t282;
t387 = t391 * qJD(1);
t386 = -t391 * t198 + t393 + (t389 + t392) * t196;
t385 = t196 * t370 + t198 * t371;
t384 = t196 * t372 + t198 * t373;
t383 = t196 ^ 2;
t382 = t198 ^ 2;
t381 = qJD(4) / 0.2e1;
t380 = -t196 * t367 + t198 * t392;
t377 = t196 * t391 + t396;
t376 = -qJD(1) * t392 + t198 * t390;
t375 = -t196 * t390 - t387;
t345 = rSges(5,1) * t197;
t270 = -rSges(5,2) * t195 + t345;
t130 = -rSges(5,3) * t198 + t196 * t270;
t319 = t197 * t198;
t189 = t196 * rSges(5,3);
t323 = t195 * t198;
t369 = -rSges(5,2) * t323 + t189;
t131 = rSges(5,1) * t319 + t369;
t168 = rSges(5,1) * t195 + rSges(5,2) * t197;
t225 = qJD(4) * t168;
t286 = t195 * t306;
t305 = qJD(1) * t198;
t210 = rSges(5,2) * t286 + rSges(5,3) * t305 - t198 * t225;
t30 = (qJD(1) * t130 + t210) * t198 + (-t196 * t225 + (-t131 + t369) * qJD(1)) * t196;
t361 = 2 * m(5);
t374 = t30 * t361;
t190 = t196 * rSges(6,1);
t368 = -rSges(6,2) * t319 + t190;
t204 = -pkin(7) - qJ(3);
t203 = cos(pkin(10));
t192 = pkin(3) * t203 + pkin(2);
t226 = -t192 - t270;
t350 = sin(qJ(1)) * pkin(1);
t102 = -t350 + (rSges(5,3) - t204) * t198 + t226 * t196;
t199 = cos(qJ(1)) * pkin(1);
t272 = t198 * t192 - t196 * t204 + t199;
t103 = t131 + t272;
t365 = t102 * t198 + t103 * t196;
t227 = rSges(4,1) * t203 - rSges(4,2) * sin(pkin(10)) + pkin(2);
t343 = rSges(4,3) + qJ(3);
t110 = t196 * t343 + t198 * t227 + t199;
t360 = 2 * m(6);
t359 = 2 * m(7);
t358 = m(6) / 0.2e1;
t357 = m(7) / 0.2e1;
t356 = t195 / 0.2e1;
t355 = t196 / 0.2e1;
t354 = -t198 / 0.2e1;
t353 = t198 / 0.2e1;
t352 = rSges(7,3) + pkin(8);
t351 = m(5) * t168;
t349 = pkin(4) * t197;
t191 = t196 * pkin(5);
t348 = qJD(1) / 0.2e1;
t205 = sin(qJ(6));
t207 = cos(qJ(6));
t335 = Icges(7,4) * t205;
t241 = Icges(7,2) * t207 + t335;
t297 = qJD(6) * t197;
t334 = Icges(7,4) * t207;
t105 = (Icges(7,2) * t205 - t334) * t297 + (Icges(7,6) * t197 + t195 * t241) * qJD(4);
t245 = Icges(7,1) * t205 + t334;
t138 = Icges(7,5) * t195 - t197 * t245;
t239 = Icges(7,5) * t205 + Icges(7,6) * t207;
t104 = (-Icges(7,5) * t207 + Icges(7,6) * t205) * t297 + (Icges(7,3) * t197 + t195 * t239) * qJD(4);
t136 = Icges(7,3) * t195 - t197 * t239;
t137 = Icges(7,6) * t195 - t197 * t241;
t299 = qJD(4) * t207;
t300 = qJD(4) * t205;
t302 = qJD(4) * t197;
t248 = t205 * t137 * t297 + t136 * t302 + (t137 * t299 + t138 * t300 + t104) * t195;
t106 = (-Icges(7,1) * t207 + t335) * t297 + (Icges(7,5) * t197 + t195 * t245) * qJD(4);
t315 = t205 * t106;
t62 = t136 * t195 + (-t137 * t207 - t138 * t205) * t197;
t347 = ((-t315 + (-qJD(6) * t138 - t105) * t207) * t197 + t248) * t195 + t62 * t302;
t277 = qJD(6) * t195 + qJD(1);
t211 = t197 * t299 - t205 * t277;
t276 = qJD(1) * t195 + qJD(6);
t230 = t196 * t276;
t75 = t198 * t211 - t207 * t230;
t212 = t197 * t300 + t207 * t277;
t76 = t198 * t212 - t205 * t230;
t346 = t76 * rSges(7,1) + t75 * rSges(7,2);
t344 = rSges(6,2) * t195;
t342 = -rSges(6,3) - qJ(5);
t316 = t198 * t207;
t321 = t196 * t205;
t145 = t195 * t316 - t321;
t317 = t198 * t205;
t320 = t196 * t207;
t146 = t195 * t317 + t320;
t97 = t146 * rSges(7,1) + t145 * rSges(7,2) + rSges(7,3) * t319;
t341 = pkin(8) * t319 + t191 + t97;
t322 = t196 * t197;
t147 = t195 * t320 + t317;
t148 = t195 * t321 - t316;
t269 = -rSges(7,1) * t148 - rSges(7,2) * t147;
t98 = rSges(7,3) * t322 - t269;
t340 = -pkin(5) * t198 + pkin(8) * t322 + t98;
t327 = qJ(5) * t195;
t326 = qJ(5) * t197;
t318 = t198 * t204;
t264 = t327 + t349;
t141 = t264 * t196;
t142 = pkin(4) * t319 + qJ(5) * t323;
t314 = t196 * t141 + t198 * t142;
t140 = qJD(4) * t264 - qJD(5) * t197;
t267 = -rSges(6,2) * t197 + rSges(6,3) * t195;
t313 = -t267 * qJD(4) - t140;
t166 = pkin(4) * t195 - t326;
t266 = rSges(6,3) * t197 + t344;
t312 = -t166 + t266;
t281 = t197 * t301;
t298 = qJD(5) * t195;
t311 = qJ(5) * t281 + t198 * t298;
t188 = qJD(3) * t198;
t310 = t204 * t306 + t188;
t309 = t382 + t383;
t304 = qJD(4) * t195;
t303 = qJD(4) * t196;
t296 = -pkin(4) - t352;
t18 = t104 * t319 + t105 * t145 + t106 * t146 - t136 * t388 + t137 * t75 + t138 * t76;
t93 = Icges(7,4) * t146 + Icges(7,2) * t145 + Icges(7,6) * t319;
t95 = Icges(7,1) * t146 + Icges(7,4) * t145 + Icges(7,5) * t319;
t251 = t205 * t95 + t207 * t93;
t38 = Icges(7,5) * t76 + Icges(7,6) * t75 - Icges(7,3) * t388;
t40 = Icges(7,4) * t76 + Icges(7,2) * t75 - Icges(7,6) * t388;
t42 = Icges(7,1) * t76 + Icges(7,4) * t75 - Icges(7,5) * t388;
t91 = Icges(7,5) * t146 + Icges(7,6) * t145 + Icges(7,3) * t319;
t9 = (qJD(4) * t251 + t38) * t195 + (qJD(4) * t91 - t205 * t42 - t207 * t40 + (t205 * t93 - t207 * t95) * qJD(6)) * t197;
t295 = t18 / 0.2e1 + t9 / 0.2e1;
t283 = t195 * t303;
t174 = pkin(4) * t283;
t284 = t197 * t305;
t294 = t141 * t305 + t196 * (pkin(4) * t284 + t196 * t298 - t174 + (t195 * t305 + t196 * t302) * qJ(5)) + t198 * (-pkin(4) * t388 - qJ(5) * t286 + t311);
t94 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t322;
t96 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t322;
t250 = t205 * t96 + t207 * t94;
t216 = -t283 + t284;
t229 = t198 * t276;
t73 = t196 * t211 + t207 * t229;
t74 = t196 * t212 + t205 * t229;
t37 = Icges(7,5) * t74 + Icges(7,6) * t73 + Icges(7,3) * t216;
t39 = Icges(7,4) * t74 + Icges(7,2) * t73 + Icges(7,6) * t216;
t41 = Icges(7,1) * t74 + Icges(7,4) * t73 + Icges(7,5) * t216;
t92 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t322;
t10 = (qJD(4) * t250 + t37) * t195 + (qJD(4) * t92 - t205 * t41 - t207 * t39 + (t205 * t94 - t207 * t96) * qJD(6)) * t197;
t17 = t104 * t322 + t105 * t147 + t106 * t148 + t136 * t216 + t137 * t73 + t138 * t74;
t291 = t10 / 0.2e1 + t17 / 0.2e1;
t32 = t195 * t92 - t197 * t250;
t47 = t136 * t322 + t137 * t147 + t138 * t148;
t290 = t32 / 0.2e1 + t47 / 0.2e1;
t31 = t195 * t91 - t197 * t251;
t46 = t136 * t319 + t137 * t145 + t138 * t146;
t289 = -t46 / 0.2e1 - t31 / 0.2e1;
t187 = qJD(3) * t196;
t288 = t187 + t311;
t287 = t174 + t310;
t280 = -t242 * qJD(4) / 0.2e1 + t240 * t381;
t279 = -t304 / 0.2e1;
t268 = rSges(7,1) * t205 + rSges(7,2) * t207;
t139 = rSges(7,3) * t195 - t197 * t268;
t278 = pkin(8) * t195 + t139;
t112 = t312 * t198;
t275 = rSges(6,1) * t305 + rSges(6,2) * t388 + rSges(6,3) * t281;
t274 = -t166 - t278;
t273 = t74 * rSges(7,1) + t73 * rSges(7,2);
t265 = -t318 - t350;
t26 = t145 * t93 + t146 * t95 + t319 * t91;
t27 = t145 * t94 + t146 * t96 + t319 * t92;
t259 = t196 * t27 + t198 * t26;
t15 = t196 * t26 - t198 * t27;
t28 = t147 * t93 + t148 * t95 + t322 * t91;
t29 = t147 * t94 + t148 * t96 + t322 * t92;
t258 = t196 * t29 + t198 * t28;
t16 = t196 * t28 - t198 * t29;
t257 = t196 * t32 + t198 * t31;
t256 = t196 * t31 - t198 * t32;
t214 = t197 * t296 - t192 - t327;
t209 = t196 * t214 - t350;
t49 = (pkin(5) - t204) * t198 + t209 + t269;
t224 = t272 + t142;
t50 = t224 + t341;
t255 = t196 * t50 + t198 * t49;
t60 = t139 * t322 - t195 * t98;
t61 = -t139 * t319 + t195 * t97;
t254 = t196 * t61 + t198 * t60;
t249 = t195 * t342 - t192;
t213 = (rSges(6,2) - pkin(4)) * t197 + t249;
t67 = -t350 + (rSges(6,1) - t204) * t198 + t213 * t196;
t132 = rSges(6,3) * t323 + t368;
t68 = t132 + t224;
t253 = t196 * t68 + t198 * t67;
t252 = t196 * t97 - t198 * t98;
t107 = (-rSges(7,1) * t207 + rSges(7,2) * t205) * t297 + (rSges(7,3) * t197 + t195 * t268) * qJD(4);
t228 = -pkin(8) * t302 - t107 - t140;
t78 = t274 * t198;
t109 = -t196 * t227 + t198 * t343 - t350;
t186 = pkin(5) * t305;
t159 = t270 * qJD(4);
t149 = t166 * t306;
t133 = -rSges(6,1) * t198 + t196 * t267;
t111 = t312 * t196;
t100 = -qJD(1) * t110 + t188;
t99 = qJD(1) * t109 + t187;
t77 = t274 * t196;
t66 = qJD(1) * t112 + t196 * t313;
t65 = t198 * t313 - t266 * t306 + t149;
t64 = t168 * t303 + (t198 * t226 - t189 - t199) * qJD(1) + t310;
t63 = t187 + ((-t192 - t345) * t196 + t265) * qJD(1) + t210;
t48 = t132 * t198 + t133 * t196 + t314;
t45 = t252 * t197;
t44 = -rSges(7,3) * t388 + t346;
t43 = rSges(7,3) * t216 + t273;
t36 = (-t298 + (t197 * t342 - t344) * qJD(4)) * t196 + (t198 * t213 - t190 - t199) * qJD(1) + t287;
t35 = -pkin(4) * t282 + ((t249 - t349) * t196 + t265) * qJD(1) + t275 + t288;
t34 = qJD(1) * t78 + t196 * t228;
t33 = t198 * t228 + t278 * t306 + t149;
t25 = t196 * t340 + t198 * t341 + t314;
t24 = (-t298 + (t195 * t352 - t326) * qJD(4)) * t196 + (t198 * t214 - t191 - t199) * qJD(1) - t273 + t287;
t23 = t186 + t296 * t282 + (t209 - t318) * qJD(1) + t288 + t346;
t21 = (-t139 * t303 - t43) * t195 + (-qJD(4) * t98 + t107 * t196 + t139 * t305) * t197;
t20 = (t139 * t301 + t44) * t195 + (qJD(4) * t97 - t107 * t198 + t139 * t306) * t197;
t19 = (qJD(1) * t133 + t275) * t198 + (t266 * t303 + (-t132 - t142 + t368) * qJD(1)) * t196 + t294;
t14 = t252 * t304 + (-t196 * t44 + t198 * t43 + (-t196 * t98 - t198 * t97) * qJD(1)) * t197;
t13 = t195 * t47 + t197 * t258;
t12 = t195 * t46 + t197 * t259;
t11 = (-pkin(8) * t282 + qJD(1) * t340 + t186 + t44) * t198 + (-pkin(8) * t283 + t43 + (-t142 - t341 + t191) * qJD(1)) * t196 + t294;
t8 = -t92 * t282 + t145 * t39 + t146 * t41 + t75 * t94 + t76 * t96 + (t198 * t37 - t306 * t92) * t197;
t7 = -t91 * t282 + t145 * t40 + t146 * t42 + t75 * t93 + t76 * t95 + (t198 * t38 - t306 * t91) * t197;
t6 = t92 * t284 + t147 * t39 + t148 * t41 + t73 * t94 + t74 * t96 + (t197 * t37 - t304 * t92) * t196;
t5 = t91 * t284 + t147 * t40 + t148 * t42 + t73 * t93 + t74 * t95 + (t197 * t38 - t304 * t91) * t196;
t4 = qJD(1) * t259 + t196 * t7 - t198 * t8;
t3 = qJD(1) * t258 + t196 * t5 - t198 * t6;
t2 = (-qJD(4) * t259 + t18) * t195 + (-qJD(1) * t15 + qJD(4) * t46 + t196 * t8 + t198 * t7) * t197;
t1 = (-qJD(4) * t258 + t17) * t195 + (-qJD(1) * t16 + qJD(4) * t47 + t196 * t6 + t198 * t5) * t197;
t22 = [(t23 * t50 + t24 * t49) * t359 + (t35 * t68 + t36 * t67) * t360 + (t102 * t64 + t103 * t63) * t361 + 0.2e1 * m(4) * (t100 * t109 + t110 * t99) + t248 - t207 * t138 * t297 + (-t105 * t207 - t315) * t197 + (t238 + t247 + t399) * t304 + (t244 + t236 - t398) * t302; 0; 0; m(7) * (qJD(1) * t255 + t196 * t24 - t198 * t23) + m(6) * (qJD(1) * t253 + t196 * t36 - t198 * t35) + m(5) * (qJD(1) * t365 + t196 * t64 - t198 * t63) + m(4) * (t100 * t196 - t198 * t99 + (t109 * t198 + t110 * t196) * qJD(1)); 0; 0; (t198 * t280 - t291) * t198 + (t196 * t280 + t295) * t196 + m(5) * ((-t196 * t63 - t198 * t64) * t168 - t365 * t159) + m(6) * (t111 * t35 + t112 * t36 + t65 * t67 + t66 * t68) + m(7) * (t23 * t77 + t24 * t78 + t33 * t49 + t34 * t50) + ((-qJD(4) * t372 + t198 * t394) * t355 + (-qJD(4) * t373 + t196 * t394) * t354 + (t354 * t370 - t355 * t371) * qJD(1)) * t195 + ((qJD(4) * t370 + t198 * t395) * t355 + (qJD(4) * t371 + t196 * t395) * t354 + (t354 * t372 - t355 * t373) * qJD(1)) * t197 + ((-t103 * t351 + (t120 / 0.2e1 - t123 / 0.2e1) * t197 + (t122 / 0.2e1 - t125 / 0.2e1) * t195 - t289) * t198 + (t102 * t351 + (t119 / 0.2e1 + t364 / 0.2e1) * t197 + (t121 / 0.2e1 + t363 / 0.2e1) * t195 + t290) * t196) * qJD(1); m(5) * t30 + m(6) * t19 + m(7) * t11; m(6) * (t196 * t65 - t198 * t66 + (t111 * t196 + t112 * t198) * qJD(1)) + m(7) * (t196 * t33 - t198 * t34 + (t196 * t77 + t198 * t78) * qJD(1)); t16 * t306 + t15 * t305 + (t11 * t25 + t33 * t78 + t34 * t77) * t359 + (t111 * t66 + t112 * t65 + t19 * t48) * t360 + t309 * t168 * t159 * t361 + (t131 * t374 + t380 * t306 + t375 * t382 - t3 + (-t386 + t393) * t305) * t198 + (t4 + t130 * t374 + t376 * t383 + t377 * t305 + ((t375 - t387) * t196 + t376 * t198 + t385 * t304 + t384 * t302 + (-t195 * t385 - t197 * t384) * qJD(4) + ((t367 + t391) * t196 - t396 + t377 + t380) * qJD(1)) * t198 + (-t196 * t389 + t386) * t306) * t196; 0.2e1 * (t253 * t358 + t255 * t357) * t302 + 0.2e1 * ((t196 * t23 + t198 * t24 + t305 * t50 - t306 * t49) * t357 + (t196 * t35 + t198 * t36 + t305 * t68 - t306 * t67) * t358) * t195; (m(6) + m(7)) * t304; 0; 0.2e1 * ((t301 * t78 + t303 * t77 - t11) * t357 + (t111 * t303 + t112 * t301 - t19) * t358) * t197 + 0.2e1 * ((qJD(4) * t25 + t196 * t34 + t198 * t33 + t305 * t77 - t306 * t78) * t357 + (qJD(4) * t48 + t111 * t305 - t112 * t306 + t196 * t66 + t198 * t65) * t358) * t195; 0.4e1 * (t358 + t357) * (-0.1e1 + t309) * t195 * t302; m(7) * (t20 * t50 + t21 * t49 + t23 * t61 + t24 * t60) + (-t196 * t290 + t198 * t289) * t304 + (t295 * t198 + t291 * t196 + (t196 * t289 + t198 * t290) * qJD(1)) * t197 + t347; m(7) * t14; m(7) * (qJD(1) * t254 + t196 * t21 - t198 * t20); m(7) * (-t11 * t45 + t14 * t25 + t20 * t77 + t21 * t78 + t33 * t60 + t34 * t61) + (-t1 / 0.2e1 + t15 * t279 + (qJD(1) * t31 - t10) * t356 + t12 * t348) * t198 + (t13 * t348 + (qJD(1) * t32 + t9) * t356 + t16 * t279 + t2 / 0.2e1) * t196 + (t4 * t353 + t256 * t381 + t3 * t355 + (-t196 * t15 / 0.2e1 + t16 * t353) * qJD(1)) * t197; m(7) * ((qJD(4) * t254 - t14) * t197 + (-qJD(4) * t45 + t196 * t20 + t198 * t21 + (-t196 * t60 + t198 * t61) * qJD(1)) * t195); (-t14 * t45 + t20 * t61 + t21 * t60) * t359 + ((-t198 * t12 - t196 * t13 - t195 * t257) * qJD(4) + t347) * t195 + (t198 * t2 + t196 * t1 + t195 * (t10 * t196 + t198 * t9) + (t195 * t62 + t197 * t257) * qJD(4) + (-t196 * t12 + t198 * t13 - t195 * t256) * qJD(1)) * t197;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t22(1) t22(2) t22(4) t22(7) t22(11) t22(16); t22(2) t22(3) t22(5) t22(8) t22(12) t22(17); t22(4) t22(5) t22(6) t22(9) t22(13) t22(18); t22(7) t22(8) t22(9) t22(10) t22(14) t22(19); t22(11) t22(12) t22(13) t22(14) t22(15) t22(20); t22(16) t22(17) t22(18) t22(19) t22(20) t22(21);];
Mq  = res;
