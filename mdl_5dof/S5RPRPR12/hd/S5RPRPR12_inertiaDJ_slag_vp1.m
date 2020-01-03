% Calculate time derivative of joint inertia matrix for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR12_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR12_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR12_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:27
% EndTime: 2019-12-31 18:29:41
% DurationCPUTime: 7.95s
% Computational Cost: add. (15175->646), mult. (16618->937), div. (0->0), fcn. (15637->10), ass. (0->312)
t217 = sin(qJ(1));
t355 = t217 / 0.2e1;
t218 = cos(qJ(1));
t375 = -t218 / 0.2e1;
t374 = -qJD(1) / 0.2e1;
t213 = cos(pkin(9));
t197 = pkin(4) * t213 + pkin(3);
t208 = pkin(8) + qJ(3);
t203 = cos(t208);
t325 = t203 * t218;
t211 = sin(pkin(9));
t319 = t217 * t211;
t196 = pkin(4) * t319;
t215 = -pkin(7) - qJ(4);
t201 = sin(t208);
t329 = t201 * t218;
t370 = -t215 * t329 + t196;
t207 = pkin(9) + qJ(5);
t200 = sin(t207);
t202 = cos(t207);
t320 = t217 * t202;
t160 = -t200 * t325 + t320;
t321 = t217 * t200;
t161 = t202 * t325 + t321;
t94 = t161 * rSges(6,1) + t160 * rSges(6,2) + rSges(6,3) * t329;
t373 = t197 * t325 + t370 + t94;
t349 = pkin(3) - t197;
t372 = t201 * t349;
t317 = qJ(4) + t215;
t371 = t203 * t317;
t340 = Icges(4,4) * t203;
t251 = -Icges(4,2) * t201 + t340;
t147 = Icges(4,6) * t217 + t218 * t251;
t341 = Icges(4,4) * t201;
t255 = Icges(4,1) * t203 - t341;
t149 = Icges(4,5) * t217 + t218 * t255;
t241 = t147 * t201 - t149 * t203;
t231 = t241 * t217;
t146 = -Icges(4,6) * t218 + t217 * t251;
t148 = -Icges(4,5) * t218 + t217 * t255;
t242 = t146 * t201 - t148 * t203;
t232 = t242 * t218;
t206 = t217 * rSges(4,3);
t369 = -rSges(4,2) * t329 + t206;
t216 = -pkin(6) - qJ(2);
t214 = cos(pkin(8));
t198 = pkin(2) * t214 + pkin(1);
t346 = rSges(4,1) * t203;
t272 = -rSges(4,2) * t201 + t346;
t236 = -t198 - t272;
t120 = (rSges(4,3) - t216) * t218 + t236 * t217;
t155 = rSges(4,1) * t325 + t369;
t277 = t218 * t198 - t217 * t216;
t121 = t155 + t277;
t368 = t120 * t218 + t121 * t217;
t247 = Icges(4,5) * t203 - Icges(4,6) * t201;
t144 = -Icges(4,3) * t218 + t217 * t247;
t274 = qJD(1) * t203 - qJD(5);
t301 = qJD(3) * t218;
t282 = t201 * t301;
t367 = t217 * t274 + t282;
t302 = qJD(3) * t217;
t280 = t201 * t302;
t366 = t218 * t274 - t280;
t238 = rSges(3,1) * t214 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t344 = rSges(3,3) + qJ(2);
t139 = t217 * t344 + t218 * t238;
t365 = 2 * m(4);
t364 = 2 * m(5);
t363 = 2 * m(6);
t209 = t217 ^ 2;
t210 = t218 ^ 2;
t362 = m(5) / 0.2e1;
t361 = m(6) / 0.2e1;
t249 = Icges(5,4) * t213 - Icges(5,2) * t211;
t141 = -Icges(5,6) * t203 + t201 * t249;
t360 = t141 / 0.2e1;
t253 = Icges(5,1) * t213 - Icges(5,4) * t211;
t142 = -Icges(5,5) * t203 + t201 * t253;
t359 = t142 / 0.2e1;
t358 = -t203 / 0.2e1;
t357 = -t211 / 0.2e1;
t356 = t213 / 0.2e1;
t353 = t218 / 0.2e1;
t180 = rSges(4,1) * t201 + rSges(4,2) * t203;
t352 = m(4) * t180;
t351 = pkin(3) * t203;
t350 = qJD(1) / 0.2e1;
t222 = -t201 * t317 - t203 * t349;
t324 = t211 * t218;
t298 = pkin(4) * t324;
t158 = -t202 * t218 - t203 * t321;
t159 = -t200 * t218 + t203 * t320;
t268 = -t159 * rSges(6,1) - t158 * rSges(6,2);
t330 = t201 * t217;
t93 = rSges(6,3) * t330 - t268;
t348 = t217 * t222 - t298 + t93;
t193 = pkin(3) * t325;
t165 = qJ(4) * t329 + t193;
t347 = -t165 + t373;
t345 = rSges(6,3) * t201;
t343 = -rSges(5,3) - qJ(4);
t342 = -rSges(6,3) + t215;
t339 = Icges(6,4) * t200;
t338 = Icges(6,4) * t202;
t252 = Icges(6,1) * t202 - t339;
t135 = -Icges(6,5) * t203 + t201 * t252;
t332 = t135 * t202;
t331 = t197 * t201;
t328 = t203 * t211;
t327 = t203 * t213;
t326 = t203 * t215;
t323 = t213 * t218;
t322 = t216 * t218;
t318 = t217 * t213;
t169 = -t203 * t324 + t318;
t170 = t203 * t323 + t319;
t112 = t170 * rSges(5,1) + t169 * rSges(5,2) + rSges(5,3) * t329;
t316 = -t112 - t165;
t267 = rSges(6,1) * t202 - rSges(6,2) * t200;
t137 = -rSges(6,3) * t203 + t201 * t267;
t315 = t137 + t371 - t372;
t266 = qJ(4) * t201 + t351;
t151 = qJD(3) * t266 - qJD(4) * t203;
t269 = rSges(5,1) * t213 - rSges(5,2) * t211;
t314 = -(rSges(5,3) * t201 + t203 * t269) * qJD(3) - t151;
t143 = -rSges(5,3) * t203 + t201 * t269;
t179 = pkin(3) * t201 - qJ(4) * t203;
t313 = -t143 - t179;
t164 = t266 * t217;
t312 = t217 * t164 + t218 * t165;
t306 = qJD(1) * t217;
t284 = t201 * t306;
t311 = qJD(1) * t298 + t215 * t284;
t300 = qJD(4) * t201;
t187 = t218 * t300;
t204 = qJD(2) * t217;
t310 = t187 + t204;
t205 = qJD(2) * t218;
t309 = t216 * t306 + t205;
t308 = t209 + t210;
t145 = Icges(4,3) * t217 + t218 * t247;
t307 = qJD(1) * t145;
t305 = qJD(1) * t218;
t304 = qJD(3) * t201;
t303 = qJD(3) * t203;
t299 = qJD(5) * t201;
t279 = t203 * t301;
t182 = qJ(4) * t279;
t186 = pkin(3) * t280;
t281 = t203 * t302;
t226 = t201 * t305 + t281;
t297 = t164 * t305 + t217 * (qJ(4) * t226 + qJD(1) * t193 + t217 * t300 - t186) + t218 * (-qJ(4) * t284 + t182 + t187 + (-t203 * t306 - t282) * pkin(3));
t275 = -qJD(5) * t203 + qJD(1);
t239 = t275 * t218;
t79 = t200 * t367 + t202 * t239;
t80 = t200 * t239 - t202 * t367;
t296 = t80 * rSges(6,1) + t79 * rSges(6,2) + rSges(6,3) * t279;
t168 = t203 * t318 - t324;
t234 = t203 * t319 + t323;
t99 = Icges(5,5) * t168 - Icges(5,6) * t234 + Icges(5,3) * t330;
t294 = t99 * t330;
t293 = t99 * t329;
t89 = Icges(6,4) * t159 + Icges(6,2) * t158 + Icges(6,6) * t330;
t91 = Icges(6,1) * t159 + Icges(6,4) * t158 + Icges(6,5) * t330;
t265 = -t200 * t89 + t202 * t91;
t87 = Icges(6,5) * t159 + Icges(6,6) * t158 + Icges(6,3) * t330;
t32 = t201 * t265 - t203 * t87;
t245 = Icges(6,5) * t202 - Icges(6,6) * t200;
t133 = -Icges(6,3) * t203 + t201 * t245;
t248 = -Icges(6,2) * t200 + t338;
t134 = -Icges(6,6) * t203 + t201 * t248;
t49 = t133 * t330 + t134 * t158 + t135 * t159;
t292 = t49 / 0.2e1 + t32 / 0.2e1;
t90 = Icges(6,4) * t161 + Icges(6,2) * t160 + Icges(6,6) * t329;
t92 = Icges(6,1) * t161 + Icges(6,4) * t160 + Icges(6,5) * t329;
t264 = -t200 * t90 + t202 * t92;
t88 = Icges(6,5) * t161 + Icges(6,6) * t160 + Icges(6,3) * t329;
t33 = t201 * t264 - t203 * t88;
t50 = t133 * t329 + t160 * t134 + t161 * t135;
t291 = t50 / 0.2e1 + t33 / 0.2e1;
t86 = (-rSges(6,1) * t200 - rSges(6,2) * t202) * t299 + (t203 * t267 + t345) * qJD(3);
t290 = -t222 * qJD(3) - t151 - t86;
t100 = Icges(5,5) * t170 + Icges(5,6) * t169 + Icges(5,3) * t329;
t289 = t100 * t330;
t288 = t100 * t329;
t126 = qJD(1) * t234 + t211 * t282;
t127 = -qJD(1) * t168 - t213 * t282;
t286 = t127 * rSges(5,1) + t126 * rSges(5,2) + rSges(5,3) * t279;
t285 = -t179 - t315;
t283 = t134 * t303;
t278 = t303 / 0.2e1;
t114 = t313 * t218;
t276 = -t197 * t203 - t198;
t71 = t285 * t218;
t240 = t217 * t275;
t81 = -t200 * t366 + t202 * t240;
t82 = t200 * t240 + t202 * t366;
t273 = t82 * rSges(6,1) + t81 * rSges(6,2);
t128 = qJD(1) * t169 + t211 * t280;
t129 = qJD(1) * t170 - t213 * t280;
t271 = -t129 * rSges(5,1) - t128 * rSges(5,2);
t270 = -t168 * rSges(5,1) + rSges(5,2) * t234;
t26 = t158 * t89 + t159 * t91 + t330 * t87;
t27 = t158 * t90 + t159 * t92 + t330 * t88;
t17 = t27 * t217 - t218 * t26;
t263 = t217 * t26 + t218 * t27;
t28 = t160 * t89 + t161 * t91 + t329 * t87;
t29 = t160 * t90 + t161 * t92 + t329 * t88;
t18 = t29 * t217 - t218 * t28;
t262 = t217 * t28 + t218 * t29;
t261 = t33 * t217 - t218 * t32;
t260 = t217 * t32 + t218 * t33;
t224 = t201 * t342 + t276;
t56 = (pkin(4) * t211 - t216) * t218 + t224 * t217 + t268;
t57 = t277 + t373;
t259 = t217 * t57 + t218 * t56;
t58 = t137 * t330 + t203 * t93;
t59 = -t137 * t329 - t203 * t94;
t258 = t217 * t59 + t218 * t58;
t227 = t201 * t343 - t198 - t351;
t219 = t217 * t227 - t322;
t72 = t219 + t270;
t73 = t277 - t316;
t257 = t217 * t73 + t218 * t72;
t256 = -t217 * t94 + t218 * t93;
t254 = Icges(4,1) * t201 + t340;
t250 = Icges(4,2) * t203 + t341;
t246 = Icges(5,5) * t213 - Icges(5,6) * t211;
t233 = qJD(3) * t180;
t230 = qJD(3) * t254;
t229 = qJD(3) * t250;
t228 = qJD(3) * (-Icges(4,5) * t201 - Icges(4,6) * t203);
t225 = t279 - t284;
t83 = (-Icges(6,5) * t200 - Icges(6,6) * t202) * t299 + (Icges(6,3) * t201 + t203 * t245) * qJD(3);
t85 = (-Icges(6,1) * t200 - t338) * t299 + (Icges(6,5) * t201 + t203 * t252) * qJD(3);
t221 = t133 * t304 - t203 * t83 + t303 * t332 + (-t134 * t299 + t201 * t85) * t202;
t220 = rSges(4,2) * t284 + rSges(4,3) * t305 - t218 * t233;
t138 = -t217 * t238 + t218 * t344;
t174 = t272 * qJD(3);
t166 = t179 * t306;
t154 = -rSges(4,3) * t218 + t217 * t272;
t132 = (Icges(5,5) * t201 + t203 * t253) * qJD(3);
t131 = (Icges(5,6) * t201 + t203 * t249) * qJD(3);
t118 = -qJD(1) * t139 + t205;
t117 = qJD(1) * t138 + t204;
t113 = t313 * t217;
t111 = rSges(5,3) * t330 - t270;
t106 = t217 * t228 + t307;
t105 = -qJD(1) * t144 + t218 * t228;
t104 = Icges(5,1) * t170 + Icges(5,4) * t169 + Icges(5,5) * t329;
t103 = Icges(5,1) * t168 - Icges(5,4) * t234 + Icges(5,5) * t330;
t102 = Icges(5,4) * t170 + Icges(5,2) * t169 + Icges(5,6) * t329;
t101 = Icges(5,4) * t168 - Icges(5,2) * t234 + Icges(5,6) * t330;
t84 = (-Icges(6,2) * t202 - t339) * t299 + (Icges(6,6) * t201 + t203 * t248) * qJD(3);
t75 = t180 * t302 + (t218 * t236 - t206) * qJD(1) + t309;
t74 = t204 + (-t322 + (-t198 - t346) * t217) * qJD(1) + t220;
t70 = t285 * t217;
t69 = t217 * t145 - t218 * t241;
t68 = t217 * t144 - t232;
t67 = -t145 * t218 - t231;
t66 = -t144 * t218 - t217 * t242;
t65 = Icges(5,1) * t129 + Icges(5,4) * t128 + Icges(5,5) * t226;
t64 = Icges(5,1) * t127 + Icges(5,4) * t126 + Icges(5,5) * t225;
t63 = Icges(5,4) * t129 + Icges(5,2) * t128 + Icges(5,6) * t226;
t62 = Icges(5,4) * t127 + Icges(5,2) * t126 + Icges(5,6) * t225;
t55 = qJD(1) * t114 + t217 * t314;
t54 = t143 * t306 + t218 * t314 + t166;
t53 = -t133 * t203 + (-t134 * t200 + t332) * t201;
t52 = t53 * t304;
t51 = t256 * t201;
t48 = t217 * t111 + t112 * t218 + t312;
t47 = rSges(6,3) * t226 + t273;
t46 = -rSges(6,3) * t284 + t296;
t45 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t226;
t44 = Icges(6,1) * t80 + Icges(6,4) * t79 + Icges(6,5) * t225;
t43 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t226;
t42 = Icges(6,4) * t80 + Icges(6,2) * t79 + Icges(6,6) * t225;
t41 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t226;
t40 = Icges(6,5) * t80 + Icges(6,6) * t79 + Icges(6,3) * t225;
t39 = t186 + (t303 * t343 - t300) * t217 + t227 * t305 + t271 + t309;
t38 = -pkin(3) * t282 + qJD(1) * t219 + t182 + t286 + t310;
t37 = t169 * t102 + t170 * t104 + t288;
t36 = t169 * t101 + t170 * t103 + t293;
t35 = -t102 * t234 + t104 * t168 + t289;
t34 = -t101 * t234 + t103 * t168 + t294;
t31 = qJD(1) * t71 + t217 * t290;
t30 = t218 * t290 + t306 * t315 + t166;
t25 = (-t300 + (t203 * t342 + t331) * qJD(3)) * t217 + (t218 * t224 - t196) * qJD(1) - t273 + t309;
t24 = (-t326 - t331) * t301 + (-t322 + (t276 - t345) * t217) * qJD(1) + t296 + t310 + t311;
t23 = t217 * t348 + t218 * t347 + t312;
t22 = (t137 * t302 + t47) * t203 + (-qJD(3) * t93 + t137 * t305 + t217 * t86) * t201;
t21 = (-t137 * t301 - t46) * t203 + (qJD(3) * t94 + t137 * t306 - t218 * t86) * t201;
t20 = (-t283 + (-qJD(5) * t135 - t84) * t201) * t200 + t221;
t19 = t217 * (rSges(5,3) * t281 - t271) + t218 * t286 + (t218 * t111 + t217 * t316) * qJD(1) + t297;
t16 = t133 * t226 + t81 * t134 + t82 * t135 + t158 * t84 + t159 * t85 + t330 * t83;
t15 = t133 * t225 + t79 * t134 + t80 * t135 + t160 * t84 + t161 * t85 + t329 * t83;
t14 = t256 * t303 + (-t217 * t46 + t218 * t47 + (-t217 * t93 - t218 * t94) * qJD(1)) * t201;
t13 = t201 * t262 - t50 * t203;
t12 = t201 * t263 - t49 * t203;
t11 = (qJD(3) * t264 - t40) * t203 + (qJD(3) * t88 - t200 * t42 + t202 * t44 + (-t200 * t92 - t202 * t90) * qJD(5)) * t201;
t10 = (qJD(3) * t265 - t41) * t203 + (qJD(3) * t87 - t200 * t43 + t202 * t45 + (-t200 * t91 - t202 * t89) * qJD(5)) * t201;
t9 = t88 * t281 + t158 * t42 + t159 * t44 + t81 * t90 + t82 * t92 + (t217 * t40 + t305 * t88) * t201;
t8 = t87 * t281 + t158 * t43 + t159 * t45 + t81 * t89 + t82 * t91 + (t217 * t41 + t305 * t87) * t201;
t7 = t88 * t279 + t160 * t42 + t161 * t44 + t79 * t90 + t80 * t92 + (t218 * t40 - t306 * t88) * t201;
t6 = t87 * t279 + t160 * t43 + t161 * t45 + t79 * t89 + t80 * t91 + (t218 * t41 - t306 * t87) * t201;
t5 = (-t182 + t46 + t311) * t218 + (t186 + t47) * t217 + (t210 * (-t326 + t372) + (-t331 - t371) * t209) * qJD(3) + (t348 * t218 + (-t165 - t347 + t370) * t217) * qJD(1) + t297;
t4 = qJD(1) * t263 + t9 * t217 - t218 * t8;
t3 = qJD(1) * t262 + t7 * t217 - t218 * t6;
t2 = (qJD(3) * t263 - t16) * t203 + (-qJD(1) * t17 + qJD(3) * t49 + t217 * t8 + t218 * t9) * t201;
t1 = (qJD(3) * t262 - t15) * t203 + (-qJD(1) * t18 + qJD(3) * t50 + t217 * t6 + t218 * t7) * t201;
t60 = [(t24 * t57 + t25 * t56) * t363 + (t38 * t73 + t39 * t72) * t364 + (t120 * t75 + t121 * t74) * t365 + 0.2e1 * m(3) * (t117 * t139 + t118 * t138) + t221 + (-t211 * t131 + t213 * t132) * t201 + (-Icges(5,3) * t203 + t201 * t246 - t250 + t255) * t304 + (-t135 * t299 - t201 * t84 - t283) * t200 + (-Icges(5,3) * t201 - t211 * t141 + t213 * t142 - t203 * t246 + t251 + t254) * t303; m(6) * (qJD(1) * t259 + t217 * t25 - t218 * t24) + m(5) * (qJD(1) * t257 + t217 * t39 - t218 * t38) + m(4) * (qJD(1) * t368 + t217 * t75 - t218 * t74) + m(3) * (-t117 * t218 + t217 * t118 + (t138 * t218 + t139 * t217) * qJD(1)); 0; m(4) * ((-t217 * t74 - t218 * t75) * t180 - t368 * t174) + m(5) * (t113 * t38 + t114 * t39 + t54 * t72 + t55 * t73) + m(6) * (t24 * t70 + t25 * t71 + t30 * t56 + t31 * t57) + ((Icges(5,5) * t129 / 0.2e1 + Icges(5,6) * t128 / 0.2e1 + Icges(5,3) * t226 / 0.2e1 + t147 * t374 + t229 * t355) * t218 + (-Icges(5,5) * t127 / 0.2e1 - Icges(5,6) * t126 / 0.2e1 - Icges(5,3) * t225 / 0.2e1 + t146 * t374 + t229 * t375) * t217) * t203 + ((-t121 * t352 + t169 * t360 + t170 * t359 + (-t100 / 0.2e1 + t147 / 0.2e1) * t203 + (t102 * t357 + t104 * t356 + t149 / 0.2e1) * t201 + t291) * t218 + (t120 * t352 - t234 * t360 + t168 * t359 + (-t99 / 0.2e1 + t146 / 0.2e1) * t203 + (t101 * t357 + t103 * t356 + t148 / 0.2e1) * t201 + t292) * t217) * qJD(1) + ((t210 / 0.2e1 + t209 / 0.2e1) * t247 + t232 / 0.2e1 - t231 / 0.2e1) * qJD(3) + (t126 * t141 + t127 * t142 + t169 * t131 + t170 * t132 + t15 + t11 + (-qJD(1) * t148 - t211 * t62 + t213 * t64 - t218 * t230) * t201 + (t100 * t201 - t102 * t328 + t104 * t327) * qJD(3)) * t355 + (t10 + t128 * t141 + t129 * t142 - t234 * t131 + t168 * t132 + t16 + (qJD(1) * t149 - t211 * t63 + t213 * t65 - t217 * t230) * t201 + (-t101 * t328 + t103 * t327 + t201 * t99) * qJD(3)) * t375; m(5) * (t54 * t217 - t218 * t55 + (t113 * t217 + t114 * t218) * qJD(1)) + m(6) * (t30 * t217 - t218 * t31 + (t217 * t70 + t218 * t71) * qJD(1)); -t218 * t4 + t217 * t3 + (t23 * t5 + t30 * t71 + t31 * t70) * t363 + (t113 * t55 + t114 * t54 + t19 * t48) * t364 - t218 * ((t218 * t106 + (t67 + t232) * qJD(1)) * t218 + (t66 * qJD(1) + (-t147 * t303 - t149 * t304 + t307) * t217 + (-t105 + (t146 * t203 + t148 * t201) * qJD(3) - t241 * qJD(1)) * t218) * t217) + ((t217 * t154 + t155 * t218) * ((qJD(1) * t154 + t220) * t218 + (-t217 * t233 + (-t155 + t369) * qJD(1)) * t217) + t308 * t180 * t174) * t365 + t217 * ((t126 * t102 + t127 * t104 + t169 * t62 + t170 * t64 + (t36 - t289) * qJD(1)) * t217 + (-t126 * t101 - t127 * t103 - t169 * t63 - t170 * t65 + (t37 + t294) * qJD(1)) * t218) + t217 * ((t217 * t105 + (t68 + t231) * qJD(1)) * t217 + (t69 * qJD(1) + (t146 * t303 + t148 * t304) * t218 + (-t106 + (-t147 * t203 - t149 * t201) * qJD(3) + (t145 - t242) * qJD(1)) * t217) * t218) - t218 * ((-t128 * t101 - t129 * t103 + t234 * t63 - t168 * t65 + (t35 - t293) * qJD(1)) * t218 + (t128 * t102 + t129 * t104 - t234 * t62 + t168 * t64 + (t34 + t288) * qJD(1)) * t217) + (t17 + (-t34 - t66) * t218 + (t35 + t67) * t217) * t306 + (t18 + (-t36 - t68) * t218 + (t37 + t69) * t217) * t305; 0.2e1 * (t257 * t362 + t259 * t361) * t303 + 0.2e1 * ((t217 * t24 + t218 * t25 + t305 * t57 - t306 * t56) * t361 + (t217 * t38 + t218 * t39 + t305 * t73 - t306 * t72) * t362) * t201; 0; 0.2e1 * ((t301 * t71 + t302 * t70 - t5) * t361 + (t113 * t302 + t114 * t301 - t19) * t362) * t203 + 0.2e1 * ((qJD(3) * t23 + t217 * t31 + t218 * t30 + t305 * t70 - t306 * t71) * t361 + (qJD(3) * t48 + t113 * t305 - t114 * t306 + t217 * t55 + t218 * t54) * t362) * t201; 0.4e1 * (t362 + t361) * (-0.1e1 + t308) * t201 * t303; t52 + m(6) * (t21 * t57 + t22 * t56 + t24 * t59 + t25 * t58) + (-t20 + (t217 * t292 + t218 * t291) * qJD(3)) * t203 + ((t15 / 0.2e1 + t11 / 0.2e1) * t218 + (t16 / 0.2e1 + t10 / 0.2e1) * t217 + (-t217 * t291 + t218 * t292) * qJD(1)) * t201; m(6) * (qJD(1) * t258 - t21 * t218 + t22 * t217); m(6) * (t14 * t23 + t21 * t70 + t22 * t71 + t30 * t58 + t31 * t59 + t5 * t51) + (t18 * t278 + (qJD(1) * t33 - t10) * t358 + t13 * t350 - t2 / 0.2e1) * t218 + (t17 * t278 + (qJD(1) * t32 + t11) * t358 + t1 / 0.2e1 + t12 * t350) * t217 + (t4 * t355 + t3 * t353 + qJD(3) * t261 / 0.2e1 + (t17 * t353 - t217 * t18 / 0.2e1) * qJD(1)) * t201; m(6) * ((qJD(3) * t258 - t14) * t203 + (qJD(3) * t51 + t21 * t217 + t218 * t22 + (-t217 * t58 + t218 * t59) * qJD(1)) * t201); (t14 * t51 + t21 * t59 + t22 * t58) * t363 + (t20 * t203 - t52 + (t217 * t12 + t218 * t13 - t203 * t260) * qJD(3)) * t203 + (t218 * t1 + t217 * t2 - t203 * (t10 * t217 + t11 * t218) + (t201 * t260 - t53 * t203) * qJD(3) + (t218 * t12 - t217 * t13 + t203 * t261) * qJD(1)) * t201;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t60(1), t60(2), t60(4), t60(7), t60(11); t60(2), t60(3), t60(5), t60(8), t60(12); t60(4), t60(5), t60(6), t60(9), t60(13); t60(7), t60(8), t60(9), t60(10), t60(14); t60(11), t60(12), t60(13), t60(14), t60(15);];
Mq = res;
