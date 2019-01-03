% Calculate joint inertia matrix for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (leg links until cut joint, platform)
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR14_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_inertiaJ_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR14_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:08:01
% EndTime: 2019-01-03 10:08:33
% DurationCPUTime: 12.13s
% Computational Cost: add. (106803->675), mult. (303819->943), div. (0->0), fcn. (404827->18), ass. (0->303)
t290 = sin(pkin(14));
t292 = sin(pkin(6));
t293 = cos(pkin(14));
t299 = sin(qJ(2));
t294 = cos(pkin(7));
t302 = cos(qJ(2));
t364 = t294 * t302;
t291 = sin(pkin(7));
t295 = cos(pkin(6));
t369 = t291 * t295;
t259 = t293 * t369 + (-t290 * t299 + t293 * t364) * t292;
t368 = t292 * t299;
t260 = t293 * t368 + (t292 * t364 + t369) * t290;
t366 = t292 * t302;
t277 = -t291 * t366 + t294 * t295;
t220 = Icges(4,5) * t260 + Icges(4,6) * t259 + Icges(4,3) * t277;
t221 = Icges(4,4) * t260 + Icges(4,2) * t259 + Icges(4,6) * t277;
t222 = Icges(4,1) * t260 + Icges(4,4) * t259 + Icges(4,5) * t277;
t269 = Icges(3,3) * t295 + (Icges(3,5) * t299 + Icges(3,6) * t302) * t292;
t270 = Icges(3,6) * t295 + (Icges(3,4) * t299 + Icges(3,2) * t302) * t292;
t271 = Icges(3,5) * t295 + (Icges(3,1) * t299 + Icges(3,4) * t302) * t292;
t386 = t277 * t220 + t259 * t221 + t260 * t222 + t295 * t269 + t270 * t366 + t271 * t368;
t303 = cos(qJ(1));
t360 = t302 * t303;
t300 = sin(qJ(1));
t362 = t300 * t299;
t281 = -t295 * t362 + t360;
t361 = t300 * t302;
t363 = t299 * t303;
t280 = -t295 * t361 - t363;
t367 = t292 * t300;
t313 = t280 * t294 + t291 * t367;
t240 = -t281 * t290 + t293 * t313;
t263 = -t280 * t291 + t294 * t367;
t372 = sin(pkin(8));
t373 = cos(pkin(8));
t219 = -t240 * t372 + t263 * t373;
t279 = t295 * t363 + t361;
t278 = t295 * t360 - t362;
t365 = t292 * t303;
t314 = t278 * t294 - t291 * t365;
t238 = -t279 * t290 + t293 * t314;
t262 = -t278 * t291 - t294 * t365;
t218 = -t238 * t372 + t262 * t373;
t385 = t386 * t295;
t239 = t279 * t293 + t290 * t314;
t298 = sin(qJ(4));
t378 = cos(qJ(4));
t190 = t239 * t378 + (t238 * t373 + t262 * t372) * t298;
t297 = sin(qJ(5));
t377 = cos(qJ(5));
t173 = t190 * t297 - t218 * t377;
t241 = t281 * t293 + t290 * t313;
t192 = t241 * t378 + (t240 * t373 + t263 * t372) * t298;
t175 = t192 * t297 - t219 * t377;
t213 = t260 * t378 + (t259 * t373 + t277 * t372) * t298;
t232 = -t259 * t372 + t277 * t373;
t180 = t213 * t297 - t232 * t377;
t174 = t190 * t377 + t218 * t297;
t320 = t378 * t372;
t321 = t373 * t378;
t189 = -t238 * t321 + t239 * t298 - t262 * t320;
t296 = sin(qJ(6));
t301 = cos(qJ(6));
t133 = -t174 * t296 + t189 * t301;
t134 = t174 * t301 + t189 * t296;
t87 = Icges(7,5) * t134 + Icges(7,6) * t133 + Icges(7,3) * t173;
t89 = Icges(7,4) * t134 + Icges(7,2) * t133 + Icges(7,6) * t173;
t91 = Icges(7,1) * t134 + Icges(7,4) * t133 + Icges(7,5) * t173;
t27 = t133 * t89 + t134 * t91 + t173 * t87;
t176 = t192 * t377 + t219 * t297;
t191 = -t240 * t321 + t241 * t298 - t263 * t320;
t135 = -t176 * t296 + t191 * t301;
t136 = t176 * t301 + t191 * t296;
t88 = Icges(7,5) * t136 + Icges(7,6) * t135 + Icges(7,3) * t175;
t90 = Icges(7,4) * t136 + Icges(7,2) * t135 + Icges(7,6) * t175;
t92 = Icges(7,1) * t136 + Icges(7,4) * t135 + Icges(7,5) * t175;
t28 = t133 * t90 + t134 * t92 + t173 * t88;
t181 = t213 * t377 + t232 * t297;
t212 = -t259 * t321 + t260 * t298 - t277 * t320;
t163 = -t181 * t296 + t212 * t301;
t164 = t181 * t301 + t212 * t296;
t106 = Icges(7,5) * t164 + Icges(7,6) * t163 + Icges(7,3) * t180;
t107 = Icges(7,4) * t164 + Icges(7,2) * t163 + Icges(7,6) * t180;
t108 = Icges(7,1) * t164 + Icges(7,4) * t163 + Icges(7,5) * t180;
t41 = t106 * t173 + t107 * t133 + t108 * t134;
t1 = t173 * t27 + t175 * t28 + t180 * t41;
t384 = t1 / 0.2e1;
t29 = t135 * t89 + t136 * t91 + t175 * t87;
t30 = t135 * t90 + t136 * t92 + t175 * t88;
t42 = t106 * t175 + t107 * t135 + t108 * t136;
t2 = t173 * t29 + t175 * t30 + t180 * t42;
t383 = t2 / 0.2e1;
t34 = t163 * t89 + t164 * t91 + t180 * t87;
t35 = t163 * t90 + t164 * t92 + t180 * t88;
t49 = t180 * t106 + t163 * t107 + t164 * t108;
t45 = t49 * t180;
t9 = t34 * t173 + t35 * t175 + t45;
t382 = t9 / 0.2e1;
t381 = t173 / 0.2e1;
t380 = t175 / 0.2e1;
t379 = t180 / 0.2e1;
t376 = t174 * pkin(5);
t319 = -t134 * rSges(7,1) - t133 * rSges(7,2);
t93 = t173 * rSges(7,3) - t319;
t375 = t173 * pkin(13) + t376 + t93;
t94 = t136 * rSges(7,1) + t135 * rSges(7,2) + t175 * rSges(7,3);
t374 = t176 * pkin(5) + t175 * pkin(13) + t94;
t109 = rSges(7,1) * t164 + rSges(7,2) * t163 + rSges(7,3) * t180;
t359 = pkin(5) * t181 + pkin(13) * t180 + t109;
t195 = t241 * pkin(3) + pkin(11) * t219;
t246 = t281 * pkin(2) + qJ(3) * t263;
t244 = t295 * t246;
t358 = t295 * t195 + t244;
t194 = t239 * pkin(3) + pkin(11) * t218;
t245 = t279 * pkin(2) + qJ(3) * t262;
t357 = -t194 - t245;
t266 = pkin(2) * t368 + qJ(3) * t277;
t356 = -t260 * pkin(3) - pkin(11) * t232 - t266;
t355 = t245 * t367 + t246 * t365;
t354 = t303 * pkin(1) + pkin(10) * t367;
t110 = Icges(6,5) * t174 - Icges(6,6) * t173 + Icges(6,3) * t189;
t112 = Icges(6,4) * t174 - Icges(6,2) * t173 + Icges(6,6) * t189;
t114 = Icges(6,1) * t174 - Icges(6,4) * t173 + Icges(6,5) * t189;
t51 = t110 * t189 - t112 * t173 + t114 * t174;
t111 = Icges(6,5) * t176 - Icges(6,6) * t175 + Icges(6,3) * t191;
t113 = Icges(6,4) * t176 - Icges(6,2) * t175 + Icges(6,6) * t191;
t115 = Icges(6,1) * t176 - Icges(6,4) * t175 + Icges(6,5) * t191;
t52 = t111 * t189 - t113 * t173 + t115 * t174;
t127 = Icges(6,5) * t181 - Icges(6,6) * t180 + Icges(6,3) * t212;
t128 = Icges(6,4) * t181 - Icges(6,2) * t180 + Icges(6,6) * t212;
t129 = Icges(6,1) * t181 - Icges(6,4) * t180 + Icges(6,5) * t212;
t65 = t127 * t189 - t128 * t173 + t129 * t174;
t13 = t189 * t51 + t191 * t52 + t212 * t65;
t3 = t189 * t27 + t191 * t28 + t212 * t41;
t352 = t3 / 0.2e1 + t13 / 0.2e1;
t53 = t110 * t191 - t112 * t175 + t114 * t176;
t54 = t111 * t191 - t113 * t175 + t115 * t176;
t66 = t127 * t191 - t128 * t175 + t129 * t176;
t14 = t189 * t53 + t191 * t54 + t212 * t66;
t4 = t189 * t29 + t191 * t30 + t212 * t42;
t351 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t218 * t51 + t219 * t52 + t232 * t65;
t5 = t218 * t27 + t219 * t28 + t232 * t41;
t350 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t218 * t53 + t219 * t54 + t232 * t66;
t6 = t218 * t29 + t219 * t30 + t232 * t42;
t349 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t65 * t295 + (t300 * t52 - t303 * t51) * t292;
t7 = t41 * t295 + (-t27 * t303 + t28 * t300) * t292;
t348 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t66 * t295 + (t300 * t54 - t303 * t53) * t292;
t8 = t42 * t295 + (-t29 * t303 + t30 * t300) * t292;
t347 = t8 / 0.2e1 + t18 / 0.2e1;
t46 = t49 * t212;
t10 = t34 * t189 + t35 * t191 + t46;
t55 = t110 * t212 - t112 * t180 + t114 * t181;
t56 = t111 * t212 - t113 * t180 + t115 * t181;
t72 = t212 * t127 - t180 * t128 + t181 * t129;
t69 = t72 * t212;
t19 = t55 * t189 + t56 * t191 + t69;
t346 = t10 / 0.2e1 + t19 / 0.2e1;
t47 = t49 * t232;
t11 = t34 * t218 + t35 * t219 + t47;
t70 = t72 * t232;
t20 = t55 * t218 + t56 * t219 + t70;
t345 = t11 / 0.2e1 + t20 / 0.2e1;
t48 = t49 * t295;
t12 = t48 + (t35 * t300 - t34 * t303) * t292;
t71 = t72 * t295;
t21 = t71 + (t56 * t300 - t55 * t303) * t292;
t344 = t12 / 0.2e1 + t21 / 0.2e1;
t343 = t34 / 0.2e1 + t41 / 0.2e1;
t342 = t35 / 0.2e1 + t42 / 0.2e1;
t165 = Icges(5,5) * t213 - Icges(5,6) * t212 + Icges(5,3) * t232;
t166 = Icges(5,4) * t213 - Icges(5,2) * t212 + Icges(5,6) * t232;
t167 = Icges(5,1) * t213 - Icges(5,4) * t212 + Icges(5,5) * t232;
t97 = t232 * t165 - t212 * t166 + t213 * t167;
t161 = t192 * pkin(4) + t191 * pkin(12);
t339 = t295 * t161 + t358;
t160 = t190 * pkin(4) + t189 * pkin(12);
t338 = -t160 + t357;
t117 = t176 * rSges(6,1) - t175 * rSges(6,2) + t191 * rSges(6,3);
t179 = pkin(4) * t213 + pkin(12) * t212;
t337 = -t179 + t356;
t145 = t192 * rSges(5,1) - t191 * rSges(5,2) + t219 * rSges(5,3);
t205 = t241 * rSges(4,1) + t240 * rSges(4,2) + t263 * rSges(4,3);
t254 = t281 * rSges(3,1) + t280 * rSges(3,2) + rSges(3,3) * t367;
t334 = -t300 * pkin(1) + pkin(10) * t365;
t324 = t292 * (-rSges(4,1) * t260 - rSges(4,2) * t259 - rSges(4,3) * t277 - t266);
t323 = t194 * t367 + t195 * t365 + t355;
t168 = rSges(5,1) * t213 - rSges(5,2) * t212 + rSges(5,3) * t232;
t322 = t292 * (-t168 + t356);
t130 = rSges(6,1) * t181 - rSges(6,2) * t180 + rSges(6,3) * t212;
t318 = t292 * (-t130 + t337);
t317 = t65 / 0.2e1 + t55 / 0.2e1 + t343;
t316 = t66 / 0.2e1 + t56 / 0.2e1 + t342;
t315 = t160 * t367 + t161 * t365 + t323;
t312 = t292 * (t337 - t359);
t204 = rSges(4,1) * t239 + rSges(4,2) * t238 + rSges(4,3) * t262;
t144 = t190 * rSges(5,1) - t189 * rSges(5,2) + t218 * rSges(5,3);
t116 = t174 * rSges(6,1) - t173 * rSges(6,2) + t189 * rSges(6,3);
t311 = -t245 + t334;
t310 = t246 + t354;
t253 = t279 * rSges(3,1) + t278 * rSges(3,2) - rSges(3,3) * t365;
t138 = Icges(5,5) * t190 - Icges(5,6) * t189 + Icges(5,3) * t218;
t140 = Icges(5,4) * t190 - Icges(5,2) * t189 + Icges(5,6) * t218;
t142 = Icges(5,1) * t190 - Icges(5,4) * t189 + Icges(5,5) * t218;
t78 = t138 * t232 - t140 * t212 + t142 * t213;
t85 = t165 * t218 - t166 * t189 + t167 * t190;
t309 = t78 / 0.2e1 + t85 / 0.2e1 + t317;
t139 = Icges(5,5) * t192 - Icges(5,6) * t191 + Icges(5,3) * t219;
t141 = Icges(5,4) * t192 - Icges(5,2) * t191 + Icges(5,6) * t219;
t143 = Icges(5,1) * t192 - Icges(5,4) * t191 + Icges(5,5) * t219;
t79 = t139 * t232 - t141 * t212 + t143 * t213;
t86 = t165 * t219 - t166 * t191 + t167 * t192;
t308 = t79 / 0.2e1 + t86 / 0.2e1 + t316;
t307 = t310 + t195;
t306 = -t194 + t311;
t305 = t161 + t307;
t304 = -t160 + t306;
t285 = rSges(2,1) * t303 - t300 * rSges(2,2);
t284 = -t300 * rSges(2,1) - rSges(2,2) * t303;
t273 = t295 * rSges(3,3) + (rSges(3,1) * t299 + rSges(3,2) * t302) * t292;
t252 = Icges(3,1) * t281 + Icges(3,4) * t280 + Icges(3,5) * t367;
t251 = Icges(3,1) * t279 + Icges(3,4) * t278 - Icges(3,5) * t365;
t250 = Icges(3,4) * t281 + Icges(3,2) * t280 + Icges(3,6) * t367;
t249 = Icges(3,4) * t279 + Icges(3,2) * t278 - Icges(3,6) * t365;
t248 = Icges(3,5) * t281 + Icges(3,6) * t280 + Icges(3,3) * t367;
t247 = Icges(3,5) * t279 + Icges(3,6) * t278 - Icges(3,3) * t365;
t243 = t254 + t354;
t242 = -t253 + t334;
t227 = -t295 * t253 - t273 * t365;
t226 = t254 * t295 - t273 * t367;
t211 = (t253 * t300 + t254 * t303) * t292;
t210 = t269 * t367 + t270 * t280 + t271 * t281;
t209 = -t269 * t365 + t278 * t270 + t279 * t271;
t203 = Icges(4,1) * t241 + Icges(4,4) * t240 + Icges(4,5) * t263;
t202 = Icges(4,1) * t239 + Icges(4,4) * t238 + Icges(4,5) * t262;
t201 = Icges(4,4) * t241 + Icges(4,2) * t240 + Icges(4,6) * t263;
t200 = Icges(4,4) * t239 + Icges(4,2) * t238 + Icges(4,6) * t262;
t199 = Icges(4,5) * t241 + Icges(4,6) * t240 + Icges(4,3) * t263;
t198 = Icges(4,5) * t239 + Icges(4,6) * t238 + Icges(4,3) * t262;
t197 = t295 * t248 + (t250 * t302 + t252 * t299) * t292;
t196 = t295 * t247 + (t249 * t302 + t251 * t299) * t292;
t178 = t310 + t205;
t177 = -t204 + t311;
t162 = t218 * t179;
t154 = (-t204 - t245) * t295 + t303 * t324;
t153 = t205 * t295 + t300 * t324 + t244;
t149 = t232 * t161;
t148 = t220 * t263 + t221 * t240 + t222 * t241;
t147 = t220 * t262 + t221 * t238 + t222 * t239;
t146 = t219 * t160;
t137 = (t204 * t300 + t205 * t303) * t292 + t355;
t124 = t199 * t277 + t201 * t259 + t203 * t260;
t123 = t198 * t277 + t200 * t259 + t202 * t260;
t119 = t307 + t145;
t118 = -t144 + t306;
t105 = t145 * t232 - t168 * t219;
t104 = -t144 * t232 + t168 * t218;
t100 = t144 * t219 - t145 * t218;
t99 = (-t144 + t357) * t295 + t303 * t322;
t98 = t145 * t295 + t300 * t322 + t358;
t96 = t97 * t295;
t95 = t97 * t232;
t84 = t305 + t117;
t83 = -t116 + t304;
t82 = (t144 * t300 + t145 * t303) * t292 + t323;
t81 = t117 * t212 - t130 * t191;
t80 = -t116 * t212 + t130 * t189;
t77 = t139 * t219 - t141 * t191 + t143 * t192;
t76 = t138 * t219 - t140 * t191 + t142 * t192;
t75 = t139 * t218 - t141 * t189 + t143 * t190;
t74 = t138 * t218 - t140 * t189 + t142 * t190;
t73 = t116 * t191 - t117 * t189;
t68 = t117 * t232 + t149 + (-t130 - t179) * t219;
t67 = t130 * t218 + t162 + (-t116 - t160) * t232;
t64 = (-t116 + t338) * t295 + t303 * t318;
t63 = t117 * t295 + t300 * t318 + t339;
t62 = t305 + t374;
t61 = -t376 + (-rSges(7,3) - pkin(13)) * t173 + t304 + t319;
t60 = -t109 * t175 + t180 * t94;
t59 = t109 * t173 - t180 * t93;
t58 = t116 * t219 + t146 + (-t117 - t161) * t218;
t57 = (t116 * t300 + t117 * t303) * t292 + t315;
t50 = -t173 * t94 + t175 * t93;
t44 = -t191 * t359 + t212 * t374;
t43 = t189 * t359 - t212 * t375;
t40 = (t338 - t375) * t295 + t303 * t312;
t39 = t295 * t374 + t300 * t312 + t339;
t38 = t149 + t374 * t232 + (-t179 - t359) * t219;
t37 = t162 + t359 * t218 + (-t160 - t375) * t232;
t36 = -t189 * t374 + t191 * t375;
t33 = t96 + (t79 * t300 - t78 * t303) * t292;
t32 = t78 * t218 + t79 * t219 + t95;
t31 = (t300 * t375 + t303 * t374) * t292 + t315;
t26 = t146 + t375 * t219 + (-t161 - t374) * t218;
t25 = t86 * t295 + (t300 * t77 - t303 * t76) * t292;
t24 = t85 * t295 + (t300 * t75 - t303 * t74) * t292;
t23 = t218 * t76 + t219 * t77 + t232 * t86;
t22 = t218 * t74 + t219 * t75 + t232 * t85;
t101 = [(t61 ^ 2 + t62 ^ 2) * m(7) + (t118 ^ 2 + t119 ^ 2) * m(5) + (t177 ^ 2 + t178 ^ 2) * m(4) + (t242 ^ 2 + t243 ^ 2) * m(3) + m(2) * (t284 ^ 2 + t285 ^ 2) + (t83 ^ 2 + t84 ^ 2) * m(6) + Icges(2,3) + t97 + t72 + t49 + t386; t71 + t96 + (t63 * t84 + t64 * t83) * m(6) + (t153 * t178 + t154 * t177) * m(4) + t48 + (t118 * t99 + t119 * t98) * m(5) + (t226 * t243 + t227 * t242) * m(3) + (t39 * t62 + t40 * t61) * m(7) + ((-t123 / 0.2e1 - t196 / 0.2e1 - t147 / 0.2e1 - t209 / 0.2e1 - t309) * t303 + (t124 / 0.2e1 + t197 / 0.2e1 + t148 / 0.2e1 + t210 / 0.2e1 + t308) * t300) * t292 + t385; (t31 ^ 2 + t39 ^ 2 + t40 ^ 2) * m(7) + (t57 ^ 2 + t63 ^ 2 + t64 ^ 2) * m(6) + (t82 ^ 2 + t98 ^ 2 + t99 ^ 2) * m(5) + (t137 ^ 2 + t153 ^ 2 + t154 ^ 2) * m(4) + (t211 ^ 2 + t226 ^ 2 + t227 ^ 2) * m(3) + (t8 + t18 + t25 + ((t199 * t263 + t201 * t240 + t203 * t241) * t300 - (t263 * t198 + t240 * t200 + t241 * t202) * t303) * t292 + (t248 * t367 + t250 * t280 + t252 * t281) * t367) * t367 + (-t7 - t17 - t24 - ((t262 * t199 + t238 * t201 + t239 * t203) * t300 - (t198 * t262 + t200 * t238 + t202 * t239) * t303) * t292 + (-t247 * t365 + t278 * t249 + t279 * t251) * t365 + (-t247 * t367 + t248 * t365 - t280 * t249 - t278 * t250 - t281 * t251 - t279 * t252) * t367) * t365 + (t12 + t21 + t33 + (t148 + t210) * t367 + (-t147 - t209) * t365 + ((-t123 - t196) * t303 + (t124 + t197) * t300) * t292 + t385) * t295; (m(4) * t177 + m(5) * t118 + m(6) * t83 + m(7) * t61) * t263 + (m(4) * t178 + m(5) * t119 + m(6) * t84 + m(7) * t62) * t262; (m(4) * t137 + m(5) * t82 + m(6) * t57 + m(7) * t31) * t277 + (m(4) * t154 + m(5) * t99 + m(6) * t64 + m(7) * t40) * t263 + (m(4) * t153 + m(5) * t98 + m(6) * t63 + m(7) * t39) * t262; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1 + m(5) / 0.2e1 + m(4) / 0.2e1) * (t262 ^ 2 + t263 ^ 2 + t277 ^ 2); t47 + (t37 * t61 + t38 * t62) * m(7) + (t67 * t83 + t68 * t84) * m(6) + t70 + t95 + (t104 * t118 + t105 * t119) * m(5) + t308 * t219 + t309 * t218; (t26 * t31 + t37 * t40 + t38 * t39) * m(7) + (t58 * t57 + t63 * t68 + t64 * t67) * m(6) + (t100 * t82 + t104 * t99 + t105 * t98) * m(5) + (t32 / 0.2e1 + t345) * t295 + (t33 / 0.2e1 + t344) * t232 + (t25 / 0.2e1 + t347) * t219 + (t24 / 0.2e1 + t348) * t218 + ((-t22 / 0.2e1 - t350) * t303 + (t23 / 0.2e1 + t349) * t300) * t292; (m(5) * t100 + m(6) * t58 + m(7) * t26) * t277 + (m(5) * t104 + m(6) * t67 + m(7) * t37) * t263 + (m(5) * t105 + m(6) * t68 + m(7) * t38) * t262; (t26 ^ 2 + t37 ^ 2 + t38 ^ 2) * m(7) + (t58 ^ 2 + t67 ^ 2 + t68 ^ 2) * m(6) + (t100 ^ 2 + t104 ^ 2 + t105 ^ 2) * m(5) + (t11 + t20 + t32) * t232 + (t6 + t16 + t23) * t219 + (t5 + t15 + t22) * t218; t46 + (t43 * t61 + t44 * t62) * m(7) + (t80 * t83 + t81 * t84) * m(6) + t69 + t316 * t191 + t317 * t189; (t31 * t36 + t39 * t44 + t40 * t43) * m(7) + (t57 * t73 + t63 * t81 + t64 * t80) * m(6) + t346 * t295 + t344 * t212 + t347 * t191 + t348 * t189 + (t300 * t351 - t303 * t352) * t292; (m(6) * t73 + m(7) * t36) * t277 + (m(6) * t80 + m(7) * t43) * t263 + (m(6) * t81 + m(7) * t44) * t262; (t26 * t36 + t37 * t43 + t38 * t44) * m(7) + (t58 * t73 + t67 * t80 + t68 * t81) * m(6) + t346 * t232 + t351 * t219 + t352 * t218 + t345 * t212 + t349 * t191 + t350 * t189; (t36 ^ 2 + t43 ^ 2 + t44 ^ 2) * m(7) + (t73 ^ 2 + t80 ^ 2 + t81 ^ 2) * m(6) + (t10 + t19) * t212 + (t4 + t14) * t191 + (t3 + t13) * t189; t45 + (t59 * t61 + t60 * t62) * m(7) + t342 * t175 + t343 * t173; (t31 * t50 + t39 * t60 + t40 * t59) * m(7) + t295 * t382 + t7 * t381 + t8 * t380 + t12 * t379 + (t300 * t383 - t303 * t1 / 0.2e1) * t292; (t262 * t60 + t263 * t59 + t277 * t50) * m(7); t218 * t384 + t232 * t382 + t11 * t379 + t5 * t381 + t219 * t383 + t6 * t380 + (t26 * t50 + t37 * t59 + t38 * t60) * m(7); (t36 * t50 + t43 * t59 + t44 * t60) * m(7) + t10 * t379 + t189 * t384 + t191 * t383 + t3 * t381 + t4 * t380 + t212 * t382; t180 * t9 + t173 * t1 + t175 * t2 + (t50 ^ 2 + t59 ^ 2 + t60 ^ 2) * m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t101(1) t101(2) t101(4) t101(7) t101(11) t101(16); t101(2) t101(3) t101(5) t101(8) t101(12) t101(17); t101(4) t101(5) t101(6) t101(9) t101(13) t101(18); t101(7) t101(8) t101(9) t101(10) t101(14) t101(19); t101(11) t101(12) t101(13) t101(14) t101(15) t101(20); t101(16) t101(17) t101(18) t101(19) t101(20) t101(21);];
Mq  = res;
