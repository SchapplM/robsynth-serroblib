% Calculate time derivative of joint inertia matrix for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:24
% EndTime: 2019-12-31 20:49:39
% DurationCPUTime: 9.31s
% Computational Cost: add. (7509->429), mult. (6566->578), div. (0->0), fcn. (4884->8), ass. (0->248)
t216 = qJ(3) + pkin(8);
t210 = sin(t216);
t211 = cos(t216);
t333 = Icges(6,5) * t211;
t335 = Icges(5,4) * t211;
t387 = -t333 + t335 + (Icges(5,1) + Icges(6,1)) * t210;
t259 = Icges(6,3) * t210 + t333;
t263 = -Icges(5,2) * t210 + t335;
t334 = Icges(6,5) * t210;
t265 = Icges(6,1) * t211 + t334;
t336 = Icges(5,4) * t210;
t266 = Icges(5,1) * t211 - t336;
t386 = (t265 + t266) * t210 - (t259 - t263) * t211;
t151 = Icges(5,2) * t211 + t336;
t221 = cos(qJ(3));
t219 = sin(qJ(3));
t338 = Icges(4,4) * t219;
t185 = Icges(4,2) * t221 + t338;
t383 = -t151 * t210 - t185 * t219;
t339 = rSges(6,3) + qJ(5);
t380 = rSges(6,1) + pkin(4);
t379 = t210 * t339 + t211 * t380;
t295 = t380 * t210;
t308 = qJD(3) * t219;
t307 = qJD(3) * t221;
t148 = -Icges(6,3) * t211 + t334;
t337 = Icges(4,4) * t221;
t186 = Icges(4,1) * t219 + t337;
t382 = -t148 * t210 - t186 * t221 - t387 * t211 - t383;
t215 = qJD(1) + qJD(2);
t323 = t215 * t219;
t381 = t339 * t211;
t217 = qJ(1) + qJ(2);
t212 = sin(t217);
t213 = cos(t217);
t378 = -t212 * t307 - t213 * t323;
t347 = rSges(4,2) * t219;
t349 = rSges(4,1) * t221;
t377 = -t347 + t349;
t149 = Icges(5,5) * t210 + Icges(5,6) * t211;
t150 = Icges(6,4) * t210 - Icges(6,6) * t211;
t184 = Icges(4,5) * t219 + Icges(4,6) * t221;
t376 = t149 + t150 + t184;
t260 = Icges(5,5) * t211 - Icges(5,6) * t210;
t261 = Icges(4,5) * t221 - Icges(4,6) * t219;
t262 = Icges(6,4) * t211 + Icges(6,6) * t210;
t374 = t382 * t215 + (t260 + t261 + t262) * qJD(3);
t264 = -Icges(4,2) * t219 + t337;
t245 = t264 * t213;
t117 = Icges(4,6) * t212 + t245;
t267 = Icges(4,1) * t221 - t338;
t248 = t267 * t213;
t119 = Icges(4,5) * t212 + t248;
t253 = t117 * t219 - t119 * t221;
t372 = t212 * t253;
t244 = t263 * t213;
t101 = Icges(5,6) * t212 + t244;
t247 = t266 * t213;
t105 = Icges(5,5) * t212 + t247;
t257 = t101 * t210 - t105 * t211;
t371 = t212 * t257;
t246 = t265 * t213;
t103 = Icges(6,4) * t212 + t246;
t240 = t259 * t213;
t95 = Icges(6,6) * t212 + t240;
t268 = t103 * t211 + t210 * t95;
t370 = t212 * t268;
t116 = -Icges(4,6) * t213 + t212 * t264;
t118 = -Icges(4,5) * t213 + t212 * t267;
t255 = t116 * t219 - t118 * t221;
t369 = t213 * t255;
t100 = -Icges(5,6) * t213 + t212 * t263;
t104 = -Icges(5,5) * t213 + t212 * t266;
t258 = t100 * t210 - t104 * t211;
t368 = t213 * t258;
t102 = -Icges(6,4) * t213 + t212 * t265;
t94 = -Icges(6,6) * t213 + t212 * t259;
t269 = t102 * t211 + t210 * t94;
t367 = t213 * t269;
t218 = -qJ(4) - pkin(7);
t325 = t213 * t218;
t207 = pkin(3) * t221 + pkin(2);
t351 = pkin(2) - t207;
t363 = t212 * t351 - t325;
t362 = 2 * m(3);
t361 = 2 * m(4);
t360 = 2 * m(5);
t359 = 2 * m(6);
t170 = t377 * qJD(3);
t355 = m(4) * t170;
t192 = rSges(4,1) * t219 + rSges(4,2) * t221;
t354 = m(4) * t192;
t220 = sin(qJ(1));
t353 = pkin(1) * t220;
t352 = pkin(3) * t219;
t204 = t212 * pkin(7);
t205 = t213 * pkin(7);
t92 = t205 - t363;
t277 = t213 * t207 - t212 * t218;
t314 = -t213 * pkin(2) - t204;
t93 = t277 + t314;
t350 = t212 * t92 + t213 * t93;
t348 = rSges(5,1) * t211;
t346 = rSges(5,2) * t210;
t345 = pkin(1) * qJD(1);
t199 = t212 * rSges(6,2);
t198 = t212 * rSges(4,3);
t197 = t212 * rSges(5,3);
t319 = -t381 + t295;
t275 = -t319 - t352;
t84 = t275 * t213;
t344 = t215 * t84;
t330 = t210 * t213;
t329 = t211 * t213;
t328 = t212 * t215;
t327 = t212 * t221;
t326 = t213 * t215;
t324 = t215 * t218;
t202 = t213 * rSges(6,2);
t322 = t379 * t212 - t202;
t321 = t380 * t329 + t339 * t330 + t199;
t175 = t212 * t346;
t320 = rSges(5,3) * t326 + t215 * t175;
t310 = qJD(3) * t212;
t294 = t210 * t310;
t318 = t380 * t294;
t300 = t212 * t323;
t317 = rSges(4,2) * t300 + rSges(4,3) * t326;
t316 = t213 * rSges(5,3) + t175;
t315 = t213 * rSges(4,3) + t212 * t347;
t313 = t212 ^ 2 + t213 ^ 2;
t312 = qJD(3) * t210;
t311 = qJD(3) * t211;
t309 = qJD(3) * t213;
t306 = qJD(5) * t210;
t191 = pkin(7) * t326;
t195 = qJD(4) * t212;
t289 = t213 * t308;
t291 = t212 * t308;
t296 = pkin(3) * t291 + qJD(4) * t213 + t212 * t324;
t305 = t212 * ((-t213 * t351 - t204) * t215 - t296) + t213 * (-pkin(3) * t289 + t215 * t363 - t191 + t195) + t92 * t326;
t304 = rSges(5,2) * t330;
t303 = t220 * t345;
t222 = cos(qJ(1));
t302 = t222 * t345;
t301 = pkin(3) * t307;
t293 = t211 * t310;
t298 = -rSges(5,1) * t294 - rSges(5,2) * t293 - t215 * t304;
t297 = -rSges(4,1) * t291 + t378 * rSges(4,2);
t286 = -pkin(2) - t349;
t156 = rSges(5,1) * t210 + rSges(5,2) * t211;
t285 = -t156 - t352;
t239 = -t207 - t379;
t228 = t212 * t239 - t325;
t276 = rSges(6,2) * t326 + t213 * t306 + t309 * t381;
t15 = t195 + (-t295 - t352) * t309 + t228 * t215 + t276;
t13 = t15 - t303;
t51 = t202 + t228;
t49 = t51 - t353;
t284 = t215 * t49 - t13;
t16 = (-t311 * t339 - t306) * t212 + (t213 * t239 - t199) * t215 + t296 + t318;
t14 = t16 - t302;
t214 = t222 * pkin(1);
t52 = t277 + t321;
t50 = t214 + t52;
t283 = t215 * t50 + t14;
t282 = t215 * t51 - t15;
t281 = t215 * t52 + t16;
t171 = pkin(3) * t300;
t249 = -t379 * qJD(3) + qJD(5) * t211 - t301;
t19 = t213 * t249 + t319 * t328 + t171;
t83 = t275 * t212;
t280 = t215 * t83 + t19;
t20 = t212 * t249 + t344;
t279 = -t20 + t344;
t278 = -t207 - t348;
t158 = t213 * rSges(3,1) - rSges(3,2) * t212;
t132 = -rSges(3,1) * t326 + rSges(3,2) * t328;
t157 = -rSges(3,1) * t212 - rSges(3,2) * t213;
t270 = t278 * t212;
t256 = t116 * t221 + t118 * t219;
t254 = t117 * t221 + t119 * t219;
t121 = t377 * t213 + t198;
t109 = rSges(5,1) * t329 + t197 - t304;
t131 = t157 * t215;
t243 = t262 * t213;
t242 = t261 * t213;
t241 = t260 * t213;
t91 = t121 - t314;
t90 = t212 * t286 + t205 + t315;
t82 = t109 + t277;
t235 = Icges(6,2) * t215 - qJD(3) * t150;
t231 = Icges(4,3) * t215 - qJD(3) * t184;
t230 = Icges(5,3) * t215 - qJD(3) * t149;
t81 = t270 + t316 - t325;
t229 = t148 * t312 + t267 * t308 + t387 * t311 + (t264 + t186) * t307 + t386 * qJD(3);
t48 = (t286 * t213 + (-rSges(4,3) - pkin(7)) * t212) * t215 - t297;
t32 = (t213 * t278 - t197) * t215 + t296 - t298;
t224 = t383 * qJD(3) + t229;
t223 = (t374 * t212 + (-t253 + t268 - t257) * qJD(3) + (-t219 * t267 - t221 * t264 - t386) * t328) * t212 / 0.2e1 - (t248 * t323 - t374 * t213 + (-t255 + t269 - t258) * qJD(3) + (t221 * t245 + (-t240 + t244) * t211 + (t246 + t247) * t210) * t215) * t213 / 0.2e1 + (t256 + (t100 - t94) * t211 + (t102 + t104) * t210 - t376 * t213 - t382 * t212) * t328 / 0.2e1 + (t254 + (t101 - t95) * t211 + (t103 + t105) * t210 - t382 * t213 + t376 * t212) * t326 / 0.2e1;
t47 = -rSges(4,2) * t213 * t307 - pkin(2) * t328 + t191 + (-t215 * t327 - t289) * rSges(4,1) + t317;
t31 = t195 + t215 * t270 + (qJD(3) * t285 - t324) * t213 + t320;
t142 = (-t346 + t348) * qJD(3);
t134 = t158 + t214;
t133 = t157 - t353;
t120 = rSges(4,1) * t327 - t315;
t115 = Icges(4,3) * t212 + t242;
t114 = -Icges(4,3) * t213 + t212 * t261;
t113 = t132 - t302;
t112 = t131 - t303;
t111 = t285 * t213;
t110 = t285 * t212;
t107 = t212 * t348 - t316;
t99 = Icges(6,2) * t212 + t243;
t98 = -Icges(6,2) * t213 + t212 * t262;
t97 = Icges(5,3) * t212 + t241;
t96 = -Icges(5,3) * t213 + t212 * t260;
t86 = t214 + t91;
t85 = t90 - t353;
t80 = t214 + t82;
t79 = t81 - t353;
t74 = t212 * t231 + t215 * t242;
t73 = t213 * t231 - t261 * t328;
t64 = t212 * t235 + t215 * t243;
t63 = t213 * t235 - t262 * t328;
t62 = t212 * t230 + t215 * t241;
t61 = t213 * t230 - t260 * t328;
t58 = t378 * pkin(3) - t142 * t212 - t156 * t326;
t57 = t156 * t328 + t171 + (-t142 - t301) * t213;
t42 = t48 - t302;
t41 = t47 - t303;
t36 = t115 * t212 - t213 * t253;
t35 = t114 * t212 - t369;
t34 = -t115 * t213 - t372;
t33 = -t114 * t213 - t212 * t255;
t30 = t32 - t302;
t29 = t31 - t303;
t28 = t212 * t97 - t213 * t257;
t27 = t212 * t96 - t368;
t26 = t212 * t99 + t213 * t268;
t25 = t212 * t98 + t367;
t24 = -t213 * t97 - t371;
t23 = -t212 * t258 - t213 * t96;
t22 = -t213 * t99 + t370;
t21 = t212 * t269 - t213 * t98;
t10 = t212 * t322 + t213 * t321 + t350;
t1 = (t215 * t322 - t295 * t309 + t276) * t213 + (t212 * t306 + (-t321 - t93 + t199) * t215 + t339 * t293 - t318) * t212 + t305;
t2 = [t229 + (t13 * t50 + t14 * t49) * t359 + (t29 * t80 + t30 * t79) * t360 + (t41 * t86 + t42 * t85) * t361 + (t112 * t134 + t113 * t133) * t362 - t151 * t312 - t185 * t308; m(6) * (t13 * t52 + t14 * t51 + t15 * t50 + t16 * t49) + m(5) * (t29 * t82 + t30 * t81 + t31 * t80 + t32 * t79) + m(4) * (t41 * t91 + t42 * t90 + t47 * t86 + t48 * t85) + m(3) * (t112 * t158 + t113 * t157 + t131 * t134 + t132 * t133) + t224; (t15 * t52 + t16 * t51) * t359 + (t31 * t82 + t32 * t81) * t360 + (t47 * t91 + t48 * t90) * t361 + (t131 * t158 + t132 * t157) * t362 + t224; ((-t215 * t86 - t42) * t213 + (t215 * t85 - t41) * t212) * t354 + m(6) * (t13 * t83 + t14 * t84 + t19 * t49 + t20 * t50) + m(5) * (t110 * t29 + t111 * t30 + t57 * t79 + t58 * t80) + (-t212 * t86 - t213 * t85) * t355 + t223; (-t212 * t91 - t213 * t90) * t355 + ((-t215 * t91 - t48) * t213 + (t215 * t90 - t47) * t212) * t354 + t223 + m(6) * (t15 * t83 + t16 * t84 + t19 * t51 + t20 * t52) + m(5) * (t110 * t31 + t111 * t32 + t57 * t81 + t58 * t82); (t111 * t57 + t110 * t58 + (t107 * t212 + t109 * t213 + t350) * ((t215 * t107 - t156 * t309 + t320) * t213 + ((-t109 - t93 + t197) * t215 + t298) * t212 + t305)) * t360 + (t1 * t10 + t19 * t84 + t20 * t83) * t359 - t213 * ((t213 * t74 + (t34 + t369) * t215) * t213 + (t33 * t215 + (-t117 * t307 - t119 * t308) * t212 + (t256 * qJD(3) - t215 * t253 - t73) * t213) * t212) - t213 * ((t213 * t62 + (t24 + t368) * t215) * t213 + (t23 * t215 + (-t101 * t311 - t105 * t312) * t212 + (-t61 + (qJD(3) * t100 + t105 * t215) * t211 + (qJD(3) * t104 - t101 * t215) * t210) * t213) * t212) - t213 * ((t213 * t64 + (t22 - t367) * t215) * t213 + (t21 * t215 + (-t103 * t312 + t311 * t95) * t212 + (-t63 + (-qJD(3) * t94 + t103 * t215) * t211 + (qJD(3) * t102 + t215 * t95) * t210) * t213) * t212) + t212 * ((t212 * t61 + (t27 + t371) * t215) * t212 + (t28 * t215 + (t100 * t311 + t104 * t312) * t213 + (-t62 + (-qJD(3) * t101 + t104 * t215) * t211 + (-qJD(3) * t105 - t100 * t215) * t210) * t212) * t213) + t212 * ((t212 * t63 + (t25 - t370) * t215) * t212 + (t26 * t215 + (t102 * t312 - t311 * t94) * t213 + (-t64 + (qJD(3) * t95 + t102 * t215) * t211 + (-qJD(3) * t103 + t215 * t94) * t210) * t212) * t213) + t212 * ((t212 * t73 + (t35 + t372) * t215) * t212 + (t36 * t215 + (t116 * t307 + t118 * t308) * t213 + (-t254 * qJD(3) - t215 * t255 - t74) * t212) * t213) + ((t120 * t212 + t121 * t213) * (((-t121 + t198) * t215 + t297) * t212 + (t215 * t120 - t192 * t309 + t317) * t213) + t313 * t192 * t170) * t361 + ((-t21 - t23 - t33) * t213 + (t22 + t24 + t34) * t212) * t328 + ((-t25 - t27 - t35) * t213 + (t26 + t28 + t36) * t212) * t326; m(6) * (t212 * t283 + t213 * t284) + m(5) * ((t215 * t79 - t29) * t213 + (t215 * t80 + t30) * t212); m(6) * (t212 * t281 + t213 * t282) + m(5) * ((t215 * t81 - t31) * t213 + (t215 * t82 + t32) * t212); m(5) * ((t111 * t215 - t58) * t213 + (t110 * t215 + t57) * t212) + m(6) * (t212 * t280 + t213 * t279); 0; m(6) * ((t212 * t50 + t213 * t49) * t311 + (-t212 * t284 + t213 * t283) * t210); m(6) * ((t212 * t52 + t213 * t51) * t311 + (-t212 * t282 + t213 * t281) * t210); m(6) * ((-t1 + (t212 * t83 + t213 * t84) * qJD(3)) * t211 + (qJD(3) * t10 - t212 * t279 + t213 * t280) * t210); 0; (-0.1e1 + t313) * t210 * t311 * t359;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
