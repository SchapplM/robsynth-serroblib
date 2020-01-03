% Calculate time derivative of joint inertia matrix for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:22
% EndTime: 2019-12-31 19:23:46
% DurationCPUTime: 11.90s
% Computational Cost: add. (6833->668), mult. (20931->886), div. (0->0), fcn. (21317->8), ass. (0->276)
t391 = -Icges(4,4) - Icges(5,6) + Icges(6,6);
t390 = Icges(5,4) - Icges(4,5) - Icges(6,5);
t389 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t241 = sin(qJ(1));
t243 = cos(qJ(1));
t240 = sin(qJ(2));
t242 = cos(qJ(2));
t348 = Icges(3,4) * t242;
t279 = -Icges(3,2) * t240 + t348;
t182 = Icges(3,6) * t241 + t243 * t279;
t388 = t182 * t242;
t239 = sin(pkin(5));
t340 = t239 * t242;
t351 = cos(pkin(8));
t352 = cos(pkin(5));
t288 = t352 * t351;
t350 = sin(pkin(8));
t194 = t240 * t350 - t242 * t288;
t387 = qJD(1) * t240;
t251 = t240 * t288 + t242 * t350;
t300 = t243 * t351;
t160 = t239 * t300 + t241 * t251;
t339 = t240 * t241;
t196 = t239 * t339 - t243 * t352;
t361 = rSges(6,1) + pkin(4);
t386 = rSges(6,2) * t160 + t196 * t361;
t356 = rSges(6,3) + qJ(5);
t248 = t251 * t243;
t301 = t241 * t351;
t162 = -t239 * t301 + t248;
t303 = t239 * t350;
t252 = t241 * t303 + t242 * t300;
t287 = t352 * t350;
t269 = t240 * t287;
t163 = -t243 * t269 + t252;
t338 = t240 * t243;
t313 = t239 * t338;
t197 = t241 * t352 + t313;
t65 = Icges(6,5) * t197 + Icges(6,6) * t162 + Icges(6,3) * t163;
t73 = Icges(5,4) * t197 - Icges(5,2) * t163 + Icges(5,6) * t162;
t81 = Icges(4,1) * t163 - Icges(4,4) * t162 + Icges(4,5) * t197;
t385 = t65 - t73 + t81;
t67 = Icges(5,5) * t197 - Icges(5,6) * t163 + Icges(5,3) * t162;
t71 = Icges(6,4) * t197 + Icges(6,2) * t162 + Icges(6,6) * t163;
t75 = Icges(4,4) * t163 - Icges(4,2) * t162 + Icges(4,6) * t197;
t384 = t67 + t71 - t75;
t69 = Icges(4,5) * t163 - Icges(4,6) * t162 + Icges(4,3) * t197;
t77 = Icges(6,1) * t197 + Icges(6,4) * t162 + Icges(6,5) * t163;
t79 = Icges(5,1) * t197 - Icges(5,4) * t163 + Icges(5,5) * t162;
t383 = t69 + t77 + t79;
t187 = t251 * qJD(2);
t266 = qJD(2) * t287;
t298 = qJD(2) * t351;
t188 = -t240 * t266 + t242 * t298;
t322 = qJD(2) * t240;
t307 = t239 * t322;
t382 = -t390 * t307 + (Icges(4,1) + Icges(5,2) + Icges(6,3)) * t188 + t391 * t187;
t381 = t389 * t307 + t391 * t188 + (Icges(4,2) + Icges(6,2) + Icges(5,3)) * t187;
t380 = (-Icges(5,1) - Icges(6,1) - Icges(4,3)) * t307 + t390 * t188 - t389 * t187;
t195 = t240 * t351 + t242 * t287;
t131 = -Icges(6,5) * t340 + Icges(6,6) * t194 + Icges(6,3) * t195;
t134 = -Icges(5,4) * t340 - Icges(5,2) * t195 + Icges(5,6) * t194;
t139 = Icges(4,1) * t195 - Icges(4,4) * t194 - Icges(4,5) * t340;
t374 = t131 - t134 + t139;
t132 = -Icges(5,5) * t340 - Icges(5,6) * t195 + Icges(5,3) * t194;
t133 = -Icges(6,4) * t340 + Icges(6,2) * t194 + Icges(6,6) * t195;
t138 = Icges(4,4) * t195 - Icges(4,2) * t194 - Icges(4,6) * t340;
t379 = -t133 + t138 - t132;
t135 = -Icges(6,1) * t340 + Icges(6,4) * t194 + Icges(6,5) * t195;
t136 = -Icges(5,1) * t340 - Icges(5,4) * t195 + Icges(5,5) * t194;
t137 = Icges(4,5) * t195 - Icges(4,6) * t194 - Icges(4,3) * t340;
t372 = t135 + t136 + t137;
t349 = Icges(3,4) * t240;
t281 = Icges(3,1) * t242 - t349;
t184 = Icges(3,5) * t241 + t243 * t281;
t272 = t182 * t240 - t184 * t242;
t378 = t241 * t272;
t181 = -Icges(3,6) * t243 + t241 * t279;
t183 = -Icges(3,5) * t243 + t241 * t281;
t274 = t181 * t240 - t183 * t242;
t377 = t243 * t274;
t376 = t194 * qJD(2);
t375 = -rSges(3,2) * t338 + t241 * rSges(3,3);
t290 = qJD(1) * t301;
t112 = qJD(1) * t248 - t239 * t290 - t241 * t376;
t291 = t243 * t303;
t161 = -t241 * t269 + t242 * t301 - t291;
t320 = qJD(2) * t242;
t341 = t239 * t241;
t170 = qJD(1) * t197 + t320 * t341;
t371 = t112 * rSges(6,2) + t161 * qJD(5) + t361 * t170;
t277 = Icges(3,5) * t242 - Icges(3,6) * t240;
t179 = -Icges(3,3) * t243 + t241 * t277;
t370 = 2 * m(3);
t369 = 2 * m(4);
t368 = 2 * m(5);
t367 = 2 * m(6);
t366 = t241 ^ 2;
t365 = t243 ^ 2;
t360 = rSges(5,2) - pkin(3);
t359 = pkin(2) * t240;
t357 = rSges(3,3) * t243;
t355 = t356 * t161 + t386;
t354 = t162 * rSges(6,2) + t356 * t163 + t361 * t197;
t337 = t242 * t243;
t232 = pkin(2) * t337;
t299 = t352 * qJ(3);
t177 = qJ(3) * t313 + t241 * t299 + t232;
t98 = t163 * pkin(3) + qJ(4) * t162;
t353 = -t177 - t98;
t344 = qJ(3) * t239;
t343 = t183 * t240;
t220 = rSges(3,1) * t240 + rSges(3,2) * t242;
t342 = t220 * t243;
t267 = pkin(2) * t242 + t240 * t344;
t178 = qJD(2) * t267 - qJD(3) * t340;
t336 = -pkin(3) * t188 - qJ(4) * t187 - qJD(4) * t194 - t178;
t335 = -t112 * qJ(4) - t160 * qJD(4);
t334 = -rSges(4,1) * t188 + rSges(4,2) * t187 - rSges(4,3) * t307 - t178;
t333 = rSges(6,2) * t194 + t356 * t195 - t361 * t340;
t142 = rSges(4,1) * t195 - rSges(4,2) * t194 - rSges(4,3) * t340;
t205 = -qJ(3) * t340 + t359;
t332 = -t142 - t205;
t150 = pkin(3) * t195 + qJ(4) * t194;
t324 = qJD(1) * t241;
t193 = t205 * t324;
t331 = t150 * t324 + t193;
t330 = -t150 - t205;
t230 = t243 * t299;
t176 = t241 * t267 - t230;
t329 = t241 * t176 + t243 * t177;
t308 = t240 * t324;
t323 = qJD(1) * t243;
t328 = rSges(3,2) * t308 + rSges(3,3) * t323;
t297 = qJD(3) * t352;
t321 = qJD(2) * t241;
t327 = -t243 * t297 - t321 * t359;
t326 = t243 * pkin(1) + t241 * pkin(7);
t180 = Icges(3,3) * t241 + t243 * t277;
t325 = qJD(1) * t180;
t319 = qJD(2) * t243;
t318 = qJD(3) * t240;
t317 = m(5) / 0.2e1 + m(6) / 0.2e1;
t316 = -pkin(3) - t356;
t306 = t240 * t319;
t256 = -t242 * t324 - t306;
t305 = t242 * t319;
t292 = t239 * t305;
t304 = t239 * t318;
t293 = qJ(3) * t292 + qJD(1) * t230 + t241 * t297 + t243 * t304;
t315 = t176 * t323 + t241 * (qJ(3) * t170 + qJD(1) * t232 + t241 * t304 + t327) + t243 * (pkin(2) * t256 - t308 * t344 + t293);
t312 = -rSges(5,1) * t307 + rSges(5,2) * t188 - rSges(5,3) * t187 + t336;
t110 = qJD(1) * t160 + t243 * t376;
t259 = qJD(1) * t269;
t111 = -qJD(1) * t291 + t195 * t319 - t241 * t259 + t242 * t290;
t169 = -qJD(1) * t196 + t292;
t311 = -t111 * rSges(4,1) + t110 * rSges(4,2) + t169 * rSges(4,3);
t141 = -rSges(5,1) * t340 - rSges(5,2) * t195 + rSges(5,3) * t194;
t310 = -t141 + t330;
t83 = t163 * rSges(4,1) - t162 * rSges(4,2) + t197 * rSges(4,3);
t309 = t169 * rSges(5,1) + t111 * rSges(5,2) - t110 * rSges(5,3);
t87 = t197 * rSges(5,1) - t163 * rSges(5,2) + t162 * rSges(5,3);
t100 = t332 * t243;
t152 = t160 * qJ(4);
t97 = pkin(3) * t161 + t152;
t296 = t241 * t97 + t243 * t98 + t329;
t295 = -rSges(6,2) * t187 - qJD(5) * t195 - t356 * t188 - t361 * t307 + t336;
t294 = t330 - t333;
t62 = t310 * t243;
t286 = rSges(3,1) * t242 - rSges(3,2) * t240;
t285 = -t170 * rSges(5,1) - t112 * rSges(5,3);
t284 = -rSges(5,1) * t196 - rSges(5,3) * t160;
t283 = t177 + t326;
t273 = t184 * t240 + t388;
t113 = -t241 * t242 * t266 + qJD(1) * t252 - t243 * t259 - t298 * t339;
t265 = -t111 * pkin(3) - t110 * qJ(4) + t162 * qJD(4);
t271 = t241 * (t113 * pkin(3) - t335) + t243 * t265 + t97 * t323 + t315;
t186 = rSges(3,1) * t337 + t375;
t58 = t294 * t243;
t270 = -pkin(1) - t286;
t264 = qJD(2) * t220;
t263 = -pkin(1) - t267;
t260 = qJD(2) * (-Icges(3,5) * t240 - Icges(3,6) * t242);
t257 = t113 * rSges(4,1) - t112 * rSges(4,2) + t170 * rSges(4,3);
t82 = rSges(4,1) * t161 - rSges(4,2) * t160 + rSges(4,3) * t196;
t255 = t263 * t241;
t254 = t283 + t98;
t253 = -t110 * rSges(6,2) + t163 * qJD(5) - t356 * t111 + t361 * t169;
t237 = t243 * pkin(7);
t250 = t230 + t237 + t255;
t249 = -t152 + t250;
t234 = pkin(7) * t323;
t247 = -pkin(2) * t306 + qJD(1) * t255 + t234 + t293;
t246 = t247 + t265;
t245 = ((-t299 - pkin(7)) * t241 + t263 * t243) * qJD(1) + (-qJ(3) * t320 - t318) * t341 - t327;
t244 = t245 + t335;
t206 = t286 * qJD(2);
t185 = t241 * t286 - t357;
t174 = t186 + t326;
t173 = t241 * t270 + t237 + t357;
t145 = t241 * t260 + t325;
t144 = -qJD(1) * t179 + t243 * t260;
t124 = t220 * t321 + ((-rSges(3,3) - pkin(7)) * t241 + t270 * t243) * qJD(1);
t123 = rSges(3,1) * t256 - rSges(3,2) * t305 - pkin(1) * t324 + t234 + t328;
t99 = t332 * t241;
t91 = t241 * t180 - t243 * t272;
t90 = t241 * t179 - t377;
t89 = -t180 * t243 - t378;
t88 = -t179 * t243 - t241 * t274;
t85 = -rSges(5,2) * t161 - t284;
t80 = Icges(4,1) * t161 - Icges(4,4) * t160 + Icges(4,5) * t196;
t78 = Icges(5,1) * t196 - Icges(5,4) * t161 + Icges(5,5) * t160;
t76 = Icges(6,1) * t196 + Icges(6,4) * t160 + Icges(6,5) * t161;
t74 = Icges(4,4) * t161 - Icges(4,2) * t160 + Icges(4,6) * t196;
t72 = Icges(5,4) * t196 - Icges(5,2) * t161 + Icges(5,6) * t160;
t70 = Icges(6,4) * t196 + Icges(6,2) * t160 + Icges(6,6) * t161;
t68 = Icges(4,5) * t161 - Icges(4,6) * t160 + Icges(4,3) * t196;
t66 = Icges(5,5) * t196 - Icges(5,6) * t161 + Icges(5,3) * t160;
t64 = Icges(6,5) * t196 + Icges(6,6) * t160 + Icges(6,3) * t161;
t61 = t310 * t241;
t60 = t283 + t83;
t59 = t250 - t82;
t57 = t294 * t241;
t56 = qJD(1) * t100 + t241 * t334;
t55 = t142 * t324 + t243 * t334 + t193;
t54 = Icges(4,1) * t113 - Icges(4,4) * t112 + Icges(4,5) * t170;
t53 = -Icges(4,1) * t111 + Icges(4,4) * t110 + Icges(4,5) * t169;
t52 = Icges(5,1) * t170 - Icges(5,4) * t113 + Icges(5,5) * t112;
t51 = Icges(5,1) * t169 + Icges(5,4) * t111 - Icges(5,5) * t110;
t50 = Icges(6,1) * t170 + Icges(6,4) * t112 + Icges(6,5) * t113;
t49 = Icges(6,1) * t169 - Icges(6,4) * t110 - Icges(6,5) * t111;
t48 = Icges(4,4) * t113 - Icges(4,2) * t112 + Icges(4,6) * t170;
t47 = -Icges(4,4) * t111 + Icges(4,2) * t110 + Icges(4,6) * t169;
t46 = Icges(5,4) * t170 - Icges(5,2) * t113 + Icges(5,6) * t112;
t45 = Icges(5,4) * t169 + Icges(5,2) * t111 - Icges(5,6) * t110;
t44 = Icges(6,4) * t170 + Icges(6,2) * t112 + Icges(6,6) * t113;
t43 = Icges(6,4) * t169 - Icges(6,2) * t110 - Icges(6,6) * t111;
t42 = Icges(4,5) * t113 - Icges(4,6) * t112 + Icges(4,3) * t170;
t41 = -Icges(4,5) * t111 + Icges(4,6) * t110 + Icges(4,3) * t169;
t40 = Icges(5,5) * t170 - Icges(5,6) * t113 + Icges(5,3) * t112;
t39 = Icges(5,5) * t169 + Icges(5,6) * t111 - Icges(5,3) * t110;
t38 = Icges(6,5) * t170 + Icges(6,6) * t112 + Icges(6,3) * t113;
t37 = Icges(6,5) * t169 - Icges(6,6) * t110 - Icges(6,3) * t111;
t34 = t254 + t87;
t33 = t161 * t360 + t249 + t284;
t31 = t245 - t257;
t30 = t247 + t311;
t29 = t241 * t82 + t243 * t83 + t329;
t28 = t254 + t354;
t27 = t161 * t316 + t249 - t386;
t26 = qJD(1) * t62 + t241 * t312;
t25 = t141 * t324 + t243 * t312 + t331;
t23 = t162 * t67 - t163 * t73 + t197 * t79;
t22 = t162 * t66 - t163 * t72 + t197 * t78;
t21 = t162 * t71 + t163 * t65 + t197 * t77;
t20 = t162 * t70 + t163 * t64 + t197 * t76;
t19 = t160 * t67 - t161 * t73 + t196 * t79;
t18 = t160 * t66 - t161 * t72 + t196 * t78;
t17 = t160 * t71 + t161 * t65 + t196 * t77;
t16 = t160 * t70 + t161 * t64 + t196 * t76;
t15 = -t162 * t75 + t163 * t81 + t197 * t69;
t14 = -t162 * t74 + t163 * t80 + t197 * t68;
t13 = -t160 * t75 + t161 * t81 + t196 * t69;
t12 = -t160 * t74 + t161 * t80 + t196 * t68;
t11 = qJD(1) * t58 + t241 * t295;
t10 = t243 * t295 + t324 * t333 + t331;
t9 = t241 * t85 + t243 * t87 + t296;
t8 = t113 * t360 + t244 + t285;
t7 = t246 + t309;
t6 = t241 * t355 + t243 * t354 + t296;
t5 = t113 * t316 + t244 - t371;
t4 = t246 + t253;
t3 = t241 * t257 + t243 * t311 + (t243 * t82 + (-t177 - t83) * t241) * qJD(1) + t315;
t2 = t241 * (-t113 * rSges(5,2) - t285) + t243 * t309 + (t243 * t85 + (-t87 + t353) * t241) * qJD(1) + t271;
t1 = t253 * t243 + (t356 * t113 + t371) * t241 + (t355 * t243 + (t353 - t354) * t241) * qJD(1) + t271;
t24 = [(t27 * t5 + t28 * t4) * t367 + (t33 * t8 + t34 * t7) * t368 + (t30 * t60 + t31 * t59) * t369 + (t123 * t174 + t124 * t173) * t370 + (-Icges(3,2) * t242 + t281 - t349) * t322 + (Icges(3,1) * t240 + t279 + t348) * t320 + t380 * t340 + t372 * t307 + t382 * t195 + t381 * t194 + t374 * t188 - t379 * t187; m(3) * (-t241 * t220 * t123 - t124 * t342 + (-t173 * t243 - t174 * t241) * t206) + m(4) * (t100 * t31 + t30 * t99 + t55 * t59 + t56 * t60) + m(6) * (t10 * t27 + t11 * t28 + t4 * t57 + t5 * t58) + m(5) * (t25 * t33 + t26 * t34 + t61 * t7 + t62 * t8) + (t366 / 0.2e1 + t365 / 0.2e1) * t277 * qJD(2) + (-m(3) * t174 * t342 + (m(3) * t173 * t220 + t343 / 0.2e1 + (t136 / 0.2e1 + t135 / 0.2e1 + t137 / 0.2e1) * t196 + (t64 / 0.2e1 + t80 / 0.2e1 - t72 / 0.2e1) * t195 + (t70 / 0.2e1 - t74 / 0.2e1 + t66 / 0.2e1) * t194 + (-t134 / 0.2e1 + t131 / 0.2e1 + t139 / 0.2e1) * t161 + (t132 / 0.2e1 + t133 / 0.2e1 - t138 / 0.2e1) * t160 + (-t76 / 0.2e1 - t68 / 0.2e1 - t78 / 0.2e1) * t340) * t241 + (-t162 * t379 + t374 * t163 + t384 * t194 + t385 * t195 + t372 * t197 - t383 * t340 + t273) * t243 / 0.2e1) * qJD(1) + (-t272 * qJD(2) - t183 * t387 + (t383 * t322 + (-t41 - t49 - t51) * t242) * t239 - t380 * t197 + (t53 + t37 - t45) * t195 + (-t47 + t43 + t39) * t194 + t385 * t188 + t384 * t187 + t372 * t169 + t382 * t163 + t381 * t162 - t374 * t111 + t379 * t110) * t241 / 0.2e1 - (-t274 * qJD(2) + qJD(1) * t388 + t184 * t387 + ((t68 + t76 + t78) * t322 + (-t42 - t50 - t52) * t242) * t239 - t380 * t196 + (t54 + t38 - t46) * t195 + (-t48 + t44 + t40) * t194 + (t80 + t64 - t72) * t188 + (-t74 + t70 + t66) * t187 + t372 * t170 + t382 * t161 + t381 * t160 + t374 * t113 - t379 * t112) * t243 / 0.2e1; (t1 * t6 + t10 * t58 + t11 * t57) * t367 + (t2 * t9 + t25 * t62 + t26 * t61) * t368 + (t100 * t55 + t29 * t3 + t56 * t99) * t369 - t243 * ((t112 * t67 - t113 * t73 + t160 * t39 - t161 * t45 + t170 * t79 + t196 * t51) * t241 - (t112 * t66 - t113 * t72 + t160 * t40 - t161 * t46 + t170 * t78 + t196 * t52) * t243 + (t18 * t241 + t19 * t243) * qJD(1)) + t241 * ((-t110 * t67 + t111 * t73 + t162 * t39 - t163 * t45 + t169 * t79 + t197 * t51) * t241 - (-t110 * t66 + t111 * t72 + t162 * t40 - t163 * t46 + t169 * t78 + t197 * t52) * t243 + (t22 * t241 + t23 * t243) * qJD(1)) + t241 * ((-t110 * t71 - t111 * t65 + t162 * t43 + t163 * t37 + t169 * t77 + t197 * t49) * t241 - (-t110 * t70 - t111 * t64 + t162 * t44 + t163 * t38 + t169 * t76 + t197 * t50) * t243 + (t20 * t241 + t21 * t243) * qJD(1)) + t241 * ((t110 * t75 - t111 * t81 - t162 * t47 + t163 * t53 + t169 * t69 + t197 * t41) * t241 - (t110 * t74 - t111 * t80 - t162 * t48 + t163 * t54 + t169 * t68 + t197 * t42) * t243 + (t14 * t241 + t15 * t243) * qJD(1)) + t241 * ((t241 * t144 + (t90 + t378) * qJD(1)) * t241 + (t91 * qJD(1) + (t181 * t320 + t183 * t322) * t243 + (-t145 - t273 * qJD(2) + (t180 - t274) * qJD(1)) * t241) * t243) - t243 * ((t112 * t71 + t113 * t65 + t160 * t43 + t161 * t37 + t170 * t77 + t196 * t49) * t241 - (t112 * t70 + t113 * t64 + t160 * t44 + t161 * t38 + t170 * t76 + t196 * t50) * t243 + (t16 * t241 + t17 * t243) * qJD(1)) - t243 * ((-t112 * t75 + t113 * t81 - t160 * t47 + t161 * t53 + t170 * t69 + t196 * t41) * t241 - (-t112 * t74 + t113 * t80 - t160 * t48 + t161 * t54 + t170 * t68 + t196 * t42) * t243 + (t12 * t241 + t13 * t243) * qJD(1)) - t243 * ((t243 * t145 + (t89 + t377) * qJD(1)) * t243 + (t88 * qJD(1) + (-t182 * t320 - t184 * t322 + t325) * t241 + (-t144 + (t181 * t242 + t343) * qJD(2) - t272 * qJD(1)) * t243) * t241) + ((t241 * t185 + t186 * t243) * ((qJD(1) * t185 - t243 * t264 + t328) * t243 + (-t241 * t264 + (-t186 + t375) * qJD(1)) * t241) + (t365 + t366) * t220 * t206) * t370 + ((-t12 - t16 - t18 - t88) * t243 + (t13 + t17 + t19 + t89) * t241) * t324 + ((-t14 - t20 - t22 - t90) * t243 + (t15 + t21 + t23 + t91) * t241) * t323; m(6) * (t169 * t27 + t170 * t28 + t196 * t4 + t197 * t5) + m(4) * (t169 * t59 + t170 * t60 + t196 * t30 + t197 * t31) + m(5) * (t169 * t33 + t170 * t34 + t196 * t7 + t197 * t8); m(6) * (t10 * t197 + t11 * t196 + t169 * t58 + t170 * t57 + (-t1 * t242 + t322 * t6) * t239) + m(5) * (t169 * t62 + t170 * t61 + t196 * t26 + t197 * t25 + (-t2 * t242 + t322 * t9) * t239) + m(4) * (t100 * t169 + t170 * t99 + t196 * t56 + t197 * t55 + (-t242 * t3 + t29 * t322) * t239); 0.4e1 * (m(4) / 0.2e1 + t317) * (-t239 ^ 2 * t240 * t320 + t169 * t197 + t170 * t196); m(6) * (-t110 * t27 + t112 * t28 + t160 * t4 + t162 * t5) + m(5) * (-t110 * t33 + t112 * t34 + t160 * t7 + t162 * t8); m(6) * (t1 * t194 + t10 * t162 + t11 * t160 - t110 * t58 + t112 * t57 + t187 * t6) + m(5) * (-t110 * t62 + t112 * t61 + t160 * t26 + t162 * t25 + t187 * t9 + t194 * t2); 0.2e1 * t317 * (-t110 * t197 + t112 * t196 + t160 * t170 + t162 * t169 + (-t187 * t242 + t194 * t322) * t239); 0.4e1 * t317 * (-t110 * t162 + t112 * t160 + t187 * t194); m(6) * (-t111 * t27 + t113 * t28 + t161 * t4 + t163 * t5); m(6) * (t1 * t195 + t10 * t163 + t11 * t161 - t111 * t58 + t113 * t57 + t188 * t6); m(6) * (-t111 * t197 + t113 * t196 + t161 * t170 + t163 * t169 + (-t188 * t242 + t195 * t322) * t239); m(6) * (-t110 * t163 - t111 * t162 + t112 * t161 + t113 * t160 + t187 * t195 + t188 * t194); (-t111 * t163 + t113 * t161 + t188 * t195) * t367;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t24(1), t24(2), t24(4), t24(7), t24(11); t24(2), t24(3), t24(5), t24(8), t24(12); t24(4), t24(5), t24(6), t24(9), t24(13); t24(7), t24(8), t24(9), t24(10), t24(14); t24(11), t24(12), t24(13), t24(14), t24(15);];
Mq = res;
