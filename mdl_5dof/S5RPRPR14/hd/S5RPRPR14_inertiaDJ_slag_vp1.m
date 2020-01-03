% Calculate time derivative of joint inertia matrix for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR14_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR14_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR14_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:31
% EndTime: 2019-12-31 18:34:44
% DurationCPUTime: 7.14s
% Computational Cost: add. (9292->580), mult. (14405->836), div. (0->0), fcn. (13371->8), ass. (0->296)
t204 = sin(qJ(3));
t205 = sin(qJ(1));
t207 = cos(qJ(3));
t292 = qJD(3) * t207;
t208 = cos(qJ(1));
t298 = qJD(1) * t208;
t381 = t204 * t298 + t205 * t292;
t199 = qJ(3) + pkin(8);
t189 = sin(t199);
t190 = cos(t199);
t294 = qJD(3) * t205;
t272 = t190 * t294;
t380 = -t189 * t298 - t272;
t358 = rSges(6,3) + pkin(7);
t278 = t358 * t190;
t354 = pkin(4) * t189;
t379 = -t278 + t354;
t340 = Icges(4,4) * t204;
t238 = Icges(4,2) * t207 + t340;
t125 = Icges(4,6) * t208 + t205 * t238;
t339 = Icges(4,4) * t207;
t243 = Icges(4,1) * t204 + t339;
t127 = Icges(4,5) * t208 + t205 * t243;
t227 = t125 * t207 + t127 * t204;
t378 = t208 * t227;
t338 = Icges(5,4) * t189;
t236 = Icges(5,2) * t190 + t338;
t116 = Icges(5,6) * t208 + t205 * t236;
t337 = Icges(5,4) * t190;
t241 = Icges(5,1) * t189 + t337;
t118 = Icges(5,5) * t208 + t205 * t241;
t229 = t116 * t190 + t118 * t189;
t377 = t208 * t229;
t258 = rSges(5,1) * t189 + rSges(5,2) * t190;
t218 = t208 * t258;
t259 = rSges(4,1) * t204 + rSges(4,2) * t207;
t219 = t208 * t259;
t203 = sin(qJ(5));
t206 = cos(qJ(5));
t335 = Icges(6,4) * t206;
t235 = -Icges(6,2) * t203 + t335;
t109 = Icges(6,6) * t189 + t190 * t235;
t336 = Icges(6,4) * t203;
t240 = Icges(6,1) * t206 - t336;
t110 = Icges(6,5) * t189 + t190 * t240;
t376 = -t109 * t206 - t110 * t203;
t275 = t207 * t298;
t280 = rSges(4,1) * t381 + rSges(4,2) * t275;
t289 = -rSges(4,3) - pkin(1) - pkin(6);
t295 = qJD(3) * t204;
t303 = qJ(2) * t298 + qJD(2) * t205;
t60 = (-rSges(4,2) * t295 + qJD(1) * t289) * t205 + t280 + t303;
t349 = rSges(4,2) * t204;
t167 = rSges(4,1) * t207 - t349;
t193 = qJD(2) * t208;
t291 = qJD(3) * t208;
t61 = t193 + t167 * t291 + (t289 * t208 + (-qJ(2) - t259) * t205) * qJD(1);
t375 = t205 * t61 - t208 * t60;
t264 = qJD(1) * t189 + qJD(5);
t374 = -t190 * t291 + t205 * t264;
t231 = Icges(5,5) * t189 + Icges(5,6) * t190;
t373 = -Icges(5,3) * t205 + t208 * t231;
t233 = Icges(4,5) * t204 + Icges(4,6) * t207;
t372 = -Icges(4,3) * t205 + t208 * t233;
t371 = -Icges(5,6) * t205 + t208 * t236;
t370 = -Icges(4,6) * t205 + t208 * t238;
t369 = -Icges(5,5) * t205 + t208 * t241;
t368 = -Icges(4,5) * t205 + t208 * t243;
t367 = 2 * m(4);
t366 = 2 * m(5);
t365 = 2 * m(6);
t200 = t205 ^ 2;
t201 = t208 ^ 2;
t364 = t189 / 0.2e1;
t363 = -t205 / 0.2e1;
t361 = t208 / 0.2e1;
t360 = rSges(3,2) - pkin(1);
t359 = -rSges(5,3) - pkin(1);
t357 = m(4) * t167;
t356 = pkin(3) * t204;
t355 = pkin(3) * t207;
t353 = pkin(4) * t190;
t352 = pkin(6) * t208;
t351 = t205 * pkin(1);
t198 = t208 * pkin(1);
t230 = Icges(6,5) * t206 - Icges(6,6) * t203;
t108 = Icges(6,3) * t189 + t190 * t230;
t293 = qJD(3) * t206;
t296 = qJD(3) * t190;
t297 = qJD(3) * t189;
t328 = t109 * t203;
t290 = qJD(5) * t190;
t70 = (-Icges(6,5) * t203 - Icges(6,6) * t206) * t290 + (Icges(6,3) * t190 - t189 * t230) * qJD(3);
t72 = (-Icges(6,1) * t203 - t335) * t290 + (Icges(6,5) * t190 - t189 * t240) * qJD(3);
t209 = t190 * t206 * t72 + t108 * t296 + t297 * t328 + (-t110 * t293 + t70) * t189;
t71 = (-Icges(6,2) * t206 - t336) * t290 + (Icges(6,6) * t190 - t189 * t235) * qJD(3);
t347 = t203 * t71;
t45 = t108 * t189 + (t110 * t206 - t328) * t190;
t350 = ((qJD(5) * t376 - t347) * t190 + t209) * t189 + t45 * t296;
t348 = rSges(5,2) * t189;
t346 = t205 * rSges(4,3);
t345 = t205 * rSges(5,3);
t197 = t208 * rSges(4,3);
t196 = t208 * rSges(5,3);
t261 = pkin(7) * t190 - t354;
t311 = t205 * t206;
t314 = t203 * t208;
t141 = t189 * t314 + t311;
t309 = t206 * t208;
t312 = t205 * t203;
t142 = -t189 * t309 + t312;
t257 = -rSges(6,1) * t142 - rSges(6,2) * t141;
t315 = t190 * t208;
t91 = rSges(6,3) * t315 - t257;
t342 = t208 * t261 + t91;
t256 = rSges(6,1) * t206 - rSges(6,2) * t203;
t77 = (-rSges(6,1) * t203 - rSges(6,2) * t206) * t290 + (rSges(6,3) * t190 - t189 * t256) * qJD(3);
t341 = -qJD(3) * t261 - t77;
t325 = t116 * t189;
t324 = t189 * t371;
t323 = t118 * t190;
t322 = t190 * t369;
t321 = t125 * t204;
t320 = t204 * t370;
t319 = t127 * t207;
t318 = t207 * t368;
t317 = t189 * t205;
t316 = t190 * t205;
t313 = t204 * t205;
t310 = t205 * t207;
t111 = rSges(6,3) * t189 + t190 * t256;
t152 = pkin(7) * t189 + t353;
t308 = t111 + t152;
t120 = rSges(5,1) * t317 + rSges(5,2) * t316 + t196;
t184 = pkin(3) * t313;
t202 = -qJ(4) - pkin(6);
t138 = t184 + (-pkin(6) - t202) * t208;
t307 = -t120 - t138;
t139 = -t189 * t312 + t309;
t140 = t189 * t311 + t314;
t306 = rSges(6,1) * t140 + rSges(6,2) * t139;
t185 = pkin(3) * t310;
t286 = pkin(3) * t291;
t305 = qJD(1) * t185 + t204 * t286;
t304 = t202 * t205 + t208 * t356;
t302 = t205 * qJ(2) + t198;
t114 = Icges(5,3) * t208 + t205 * t231;
t301 = qJD(1) * t114;
t123 = Icges(4,3) * t208 + t205 * t233;
t300 = qJD(1) * t123;
t299 = qJD(1) * t205;
t170 = pkin(4) * t317;
t274 = t189 * t294;
t265 = qJD(5) * t189 + qJD(1);
t75 = -t265 * t311 + (-t208 * t264 - t272) * t203;
t76 = t264 * t309 + (t190 * t293 - t203 * t265) * t205;
t288 = rSges(6,1) * t76 + rSges(6,2) * t75 + rSges(6,3) * t274;
t287 = pkin(3) * t295;
t86 = Icges(6,4) * t140 + Icges(6,2) * t139 - Icges(6,6) * t316;
t88 = Icges(6,1) * t140 + Icges(6,4) * t139 - Icges(6,5) * t316;
t253 = t203 * t86 - t206 * t88;
t84 = Icges(6,5) * t140 + Icges(6,6) * t139 - Icges(6,3) * t316;
t29 = t189 * t84 - t190 * t253;
t41 = -t108 * t316 + t109 * t139 + t110 * t140;
t285 = t41 / 0.2e1 + t29 / 0.2e1;
t87 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t315;
t89 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t315;
t252 = t203 * t87 - t206 * t89;
t85 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t315;
t30 = t189 * t85 - t190 * t252;
t42 = t108 * t315 + t109 * t141 + t110 * t142;
t284 = -t42 / 0.2e1 - t30 / 0.2e1;
t90 = -rSges(6,3) * t316 + t306;
t283 = pkin(7) * t316 - t138 - t170 - t90;
t269 = t190 * t298;
t282 = rSges(5,1) * t380 - rSges(5,2) * t269;
t281 = pkin(4) * t380 - pkin(7) * t274;
t131 = rSges(4,1) * t313 + rSges(4,2) * t310 + t197;
t195 = t208 * qJ(2);
t279 = t195 + t304;
t273 = t189 * t291;
t268 = -qJ(2) - t356;
t156 = t259 * qJD(3);
t267 = (t200 + t201) * t156;
t266 = qJD(1) * t308;
t263 = pkin(3) * t381 + qJD(4) * t208 + t202 * t299;
t225 = t208 * t265;
t73 = -t203 * t374 + t206 * t225;
t74 = t203 * t225 + t206 * t374;
t262 = t74 * rSges(6,1) + t73 * rSges(6,2);
t260 = qJD(4) * t205 - t202 * t298 - t207 * t286;
t151 = rSges(5,1) * t190 - t348;
t24 = t139 * t86 + t140 * t88 - t316 * t84;
t25 = t139 * t87 + t140 * t89 - t316 * t85;
t17 = t205 * t25 + t208 * t24;
t249 = t205 * t24 - t208 * t25;
t26 = t141 * t86 + t142 * t88 + t315 * t84;
t27 = t141 * t87 + t142 * t89 + t315 * t85;
t18 = t205 * t27 + t208 * t26;
t248 = t205 * t26 - t208 * t27;
t247 = t205 * t30 + t208 * t29;
t246 = t205 * t29 - t208 * t30;
t245 = t205 * t91 + t208 * t90;
t244 = Icges(4,1) * t207 - t340;
t242 = Icges(5,1) * t190 - t338;
t239 = -Icges(4,2) * t204 + t339;
t237 = -Icges(5,2) * t189 + t337;
t234 = Icges(4,5) * t207 - Icges(4,6) * t204;
t232 = Icges(5,5) * t190 - Icges(5,6) * t189;
t228 = -t189 * t369 - t190 * t371;
t226 = -t204 * t368 - t207 * t370;
t224 = -t202 * t208 + t184 + t302;
t223 = t193 - t260;
t221 = t263 + t303;
t220 = rSges(3,3) * t208 + t205 * t360;
t217 = t228 * t205;
t216 = t226 * t205;
t215 = qJD(3) * t244;
t214 = qJD(3) * t242;
t213 = qJD(3) * t239;
t212 = qJD(3) * t237;
t211 = -t269 + t274;
t210 = -t190 * t299 - t273;
t179 = pkin(3) * t275;
t146 = t258 * qJD(3);
t137 = -pkin(6) * t205 - t304;
t136 = -rSges(3,2) * t208 + rSges(3,3) * t205 + t302;
t135 = t195 + t220;
t132 = t346 - t219;
t122 = t208 * t137;
t121 = t345 - t218;
t113 = (-t151 - t355) * t208;
t112 = t151 * t205 + t185;
t107 = t193 + (t360 * t208 + (-rSges(3,3) - qJ(2)) * t205) * qJD(1);
t106 = qJD(1) * t220 + t303;
t105 = t131 + t302 + t352;
t104 = t205 * t289 + t195 + t219;
t102 = pkin(6) * t299 + t263;
t100 = t208 * ((t184 - t352) * qJD(1) + t260);
t95 = qJD(1) * t372 + t234 * t294;
t94 = -t234 * t291 + t300;
t93 = t224 + t120;
t92 = t205 * t359 + t218 + t279;
t79 = qJD(1) * t373 + t232 * t294;
t78 = -t232 * t291 + t301;
t69 = (-t308 - t355) * t208;
t68 = t205 * t308 + t185;
t65 = t151 * t298 + t179 + (-t146 - t287) * t205;
t64 = t146 * t208 + t151 * t299 + t305;
t59 = -t205 * t372 - t208 * t226;
t58 = t123 * t205 - t378;
t57 = -t208 * t372 + t216;
t56 = t123 * t208 + t205 * t227;
t55 = -t205 * t278 + t170 + t224 + t306;
t54 = t208 * t379 + t257 + t279 - t351;
t53 = t111 * t315 - t189 * t91;
t52 = t111 * t316 + t189 * t90;
t51 = -t205 * t373 - t208 * t228;
t50 = t114 * t205 - t377;
t49 = -t208 * t373 + t217;
t48 = t114 * t208 + t205 * t229;
t47 = t151 * t291 + (t359 * t208 + (-t258 + t268) * t205) * qJD(1) + t223;
t46 = (-rSges(5,2) * t297 + qJD(1) * t359) * t205 + t221 - t282;
t43 = t245 * t190;
t40 = -rSges(6,3) * t269 + t288;
t39 = rSges(6,3) * t210 + t262;
t38 = Icges(6,1) * t76 + Icges(6,4) * t75 + Icges(6,5) * t211;
t37 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t210;
t36 = Icges(6,4) * t76 + Icges(6,2) * t75 + Icges(6,6) * t211;
t35 = Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t210;
t34 = Icges(6,5) * t76 + Icges(6,6) * t75 + Icges(6,3) * t211;
t33 = Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t210;
t32 = t179 + t208 * t266 + (-t287 - t341) * t205;
t31 = t205 * t266 + t208 * t341 + t305;
t28 = t205 * t283 + t208 * t342 + t122;
t23 = (t189 * t358 + t353) * t291 + (-t198 + (t268 - t379) * t205) * qJD(1) + t223 - t262;
t22 = (-t208 * t278 - t351) * qJD(1) + t221 - t281 + t288;
t21 = (-t111 * t294 + t40) * t189 + (qJD(3) * t90 + t111 * t298 + t205 * t77) * t190;
t20 = (-t111 * t291 - t39) * t189 + (-qJD(3) * t91 - t111 * t299 + t208 * t77) * t190;
t16 = t108 * t211 + t109 * t75 + t110 * t76 + t139 * t71 + t140 * t72 - t316 * t70;
t15 = t108 * t210 + t109 * t73 + t110 * t74 + t141 * t71 + t142 * t72 + t315 * t70;
t14 = t245 * t297 + (-t205 * t39 - t208 * t40 + (t205 * t90 - t208 * t91) * qJD(1)) * t190;
t13 = t189 * t42 - t190 * t248;
t12 = t189 * t41 - t190 * t249;
t11 = t100 + (-t102 - t40 + (-t137 - t342) * qJD(1) + t281) * t205 + (t39 - t152 * t291 + (t283 + t170) * qJD(1)) * t208;
t10 = (qJD(3) * t252 + t33) * t189 + (qJD(3) * t85 - t203 * t35 + t206 * t37 + (-t203 * t89 - t206 * t87) * qJD(5)) * t190;
t9 = (qJD(3) * t253 + t34) * t189 + (qJD(3) * t84 - t203 * t36 + t206 * t38 + (-t203 * t88 - t206 * t86) * qJD(5)) * t190;
t8 = t85 * t274 + t139 * t35 + t140 * t37 + t75 * t87 + t76 * t89 + (-t205 * t33 - t298 * t85) * t190;
t7 = t84 * t274 + t139 * t36 + t140 * t38 + t75 * t86 + t76 * t88 + (-t205 * t34 - t298 * t84) * t190;
t6 = -t85 * t273 + t141 * t35 + t142 * t37 + t73 * t87 + t74 * t89 + (t208 * t33 - t299 * t85) * t190;
t5 = -t84 * t273 + t141 * t36 + t142 * t38 + t73 * t86 + t74 * t88 + (t208 * t34 - t299 * t84) * t190;
t4 = -qJD(1) * t249 + t205 * t8 + t208 * t7;
t3 = -qJD(1) * t248 + t6 * t205 + t208 * t5;
t2 = (qJD(3) * t249 + t16) * t189 + (-qJD(1) * t17 + qJD(3) * t41 - t205 * t7 + t208 * t8) * t190;
t1 = (qJD(3) * t248 + t15) * t189 + (-qJD(1) * t18 + qJD(3) * t42 - t205 * t5 + t208 * t6) * t190;
t19 = [-t189 * t214 + t236 * t297 - t241 * t296 + t209 - t207 * t213 - t204 * t215 + t238 * t295 - t243 * t292 + (t22 * t55 + t23 * t54) * t365 + (t46 * t93 + t47 * t92) * t366 + (t104 * t61 + t105 * t60) * t367 + 0.2e1 * m(3) * (t106 * t136 + t107 * t135) + t376 * t290 + (-t212 - t347) * t190; m(6) * (t205 * t23 - t208 * t22 + (t205 * t55 + t208 * t54) * qJD(1)) + m(5) * (t205 * t47 - t208 * t46 + (t205 * t93 + t208 * t92) * qJD(1)) + m(4) * ((t104 * t208 + t105 * t205) * qJD(1) + t375) + m(3) * (-t106 * t208 + t205 * t107 + (t135 * t208 + t136 * t205) * qJD(1)); 0; m(4) * (t375 * t167 - (t104 * t205 - t105 * t208) * t156) + m(6) * (t22 * t69 + t23 * t68 + t31 * t55 + t32 * t54) + m(5) * (t112 * t47 + t113 * t46 + t64 * t93 + t65 * t92) + ((t105 * t357 + t321 / 0.2e1 - t319 / 0.2e1 + t325 / 0.2e1 - t323 / 0.2e1 - t285) * t205 + (t320 / 0.2e1 - t318 / 0.2e1 + t104 * t357 + t324 / 0.2e1 - t322 / 0.2e1 - t284) * t208) * qJD(1) + (-t233 - t231) * qJD(3) * (t201 / 0.2e1 + t200 / 0.2e1) + (-t189 * (qJD(1) * t116 - t237 * t291) + t190 * (qJD(1) * t118 - t242 * t291) - t204 * (qJD(1) * t125 - t239 * t291) + t207 * (qJD(1) * t127 - t244 * t291) + t10 + t15 + (-t226 - t228) * qJD(3)) * t205 / 0.2e1 + (-t189 * (qJD(1) * t371 + t205 * t212) + t190 * (qJD(1) * t369 + t205 * t214) - t204 * (qJD(1) * t370 + t205 * t213) + t207 * (qJD(1) * t368 + t205 * t215) + t16 + t9 + (-t227 - t229) * qJD(3)) * t361; m(5) * (t65 * t205 - t208 * t64 + (t112 * t208 + t113 * t205) * qJD(1)) + m(6) * (t32 * t205 - t208 * t31 + (t205 * t69 + t208 * t68) * qJD(1)) - m(4) * t267; (t11 * t28 + t31 * t69 + t32 * t68) * t365 + t208 * t4 + t205 * t3 + (t112 * t65 + t113 * t64 + (t121 * t208 + t205 * t307 + t122) * (t100 + (-t102 + t282) * t205 + (-t151 * t201 + t200 * t348) * qJD(3) + ((t307 + t196) * t208 + (-t121 - t137 + t218 + t345) * t205) * qJD(1))) * t366 + t205 * ((t205 * t78 + (-t50 + t217) * qJD(1)) * t205 + (t51 * qJD(1) + (t116 * t297 - t118 * t296 + t301) * t208 + (t79 + (t322 - t324) * qJD(3) + t229 * qJD(1)) * t205) * t208) + ((-t131 * t205 + t132 * t208) * (-t205 * t280 + (-t167 * t201 + t200 * t349) * qJD(3) + ((-t131 + t197) * t208 + (-t132 + t219 + t346) * t205) * qJD(1)) - t167 * t267) * t367 + t208 * ((t208 * t95 + (t57 + t378) * qJD(1)) * t208 + (-t56 * qJD(1) + (-t292 * t368 + t295 * t370) * t205 + (t94 + (t319 - t321) * qJD(3) + (-t123 + t226) * qJD(1)) * t208) * t205) + t208 * ((t208 * t79 + (t49 + t377) * qJD(1)) * t208 + (-t48 * qJD(1) + (-t296 * t369 + t297 * t371) * t205 + (t78 + (t323 - t325) * qJD(3) + (-t114 + t228) * qJD(1)) * t208) * t205) + t205 * ((t205 * t94 + (-t58 + t216) * qJD(1)) * t205 + (t59 * qJD(1) + (t125 * t295 - t127 * t292 + t300) * t208 + (t95 + (t318 - t320) * qJD(3) + t227 * qJD(1)) * t205) * t208) + (-t17 + (-t48 - t56) * t208 + (-t49 - t57) * t205) * t299 + (t18 + (t50 + t58) * t208 + (t51 + t59) * t205) * t298; m(6) * (t205 * t22 + t208 * t23 + (-t205 * t54 + t208 * t55) * qJD(1)) + m(5) * (t205 * t46 + t208 * t47 + (-t205 * t92 + t208 * t93) * qJD(1)); 0; m(6) * (t205 * t31 + t208 * t32 + (-t205 * t68 + t208 * t69) * qJD(1)) + m(5) * (t205 * t64 + t208 * t65 + (-t112 * t205 + t113 * t208) * qJD(1)); 0; m(6) * (t20 * t54 + t21 * t55 + t22 * t52 + t23 * t53) + (t205 * t285 + t208 * t284) * t297 + ((t15 / 0.2e1 + t10 / 0.2e1) * t208 + (-t16 / 0.2e1 - t9 / 0.2e1) * t205 + (t205 * t284 - t208 * t285) * qJD(1)) * t190 + t350; m(6) * (t20 * t205 - t208 * t21 + (t205 * t52 + t208 * t53) * qJD(1)); m(6) * (-t11 * t43 + t14 * t28 + t20 * t68 + t21 * t69 + t31 * t52 + t32 * t53) + ((qJD(1) * t30 + t9) * t364 + t2 / 0.2e1 + qJD(1) * t13 / 0.2e1 - t18 * t297 / 0.2e1) * t208 + ((-qJD(1) * t29 + t10) * t364 - qJD(1) * t12 / 0.2e1 + t1 / 0.2e1 + t17 * t297 / 0.2e1) * t205 + (qJD(3) * t247 / 0.2e1 + t3 * t361 + t4 * t363 + (t18 * t363 - t208 * t17 / 0.2e1) * qJD(1)) * t190; m(6) * (t20 * t208 + t21 * t205 + (-t205 * t53 + t208 * t52) * qJD(1)); (-t14 * t43 + t20 * t53 + t21 * t52) * t365 + ((t12 * t205 - t13 * t208 + t189 * t246) * qJD(3) + t350) * t189 + (-t205 * t2 + t208 * t1 + t189 * (t10 * t208 - t205 * t9) + (t189 * t45 - t190 * t246) * qJD(3) + (-t12 * t208 - t13 * t205 - t189 * t247) * qJD(1)) * t190;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
