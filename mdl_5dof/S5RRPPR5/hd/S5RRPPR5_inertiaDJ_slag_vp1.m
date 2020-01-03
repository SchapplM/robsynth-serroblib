% Calculate time derivative of joint inertia matrix for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:54
% EndTime: 2019-12-31 19:29:15
% DurationCPUTime: 12.96s
% Computational Cost: add. (9215->598), mult. (13625->842), div. (0->0), fcn. (12450->8), ass. (0->303)
t217 = qJ(2) + pkin(8);
t209 = sin(t217);
t210 = cos(t217);
t359 = Icges(5,5) * t210;
t364 = Icges(4,4) * t210;
t411 = -t359 + t364 + (Icges(4,1) + Icges(5,1)) * t209;
t410 = t411 * qJD(2);
t221 = sin(qJ(5));
t224 = cos(qJ(5));
t166 = t209 * t224 - t210 * t221;
t399 = qJD(2) - qJD(5);
t403 = t399 * t166;
t409 = t403 / 0.2e1;
t262 = t209 * t221 + t210 * t224;
t106 = t399 * t262;
t408 = t106 / 0.2e1;
t223 = sin(qJ(1));
t407 = t223 / 0.2e1;
t226 = cos(qJ(1));
t406 = -t226 / 0.2e1;
t405 = -qJD(1) / 0.2e1;
t404 = qJD(1) / 0.2e1;
t402 = Icges(6,5) * t408 + Icges(6,6) * t409;
t225 = cos(qJ(2));
t387 = pkin(2) * t225;
t205 = pkin(1) + t387;
t380 = pkin(1) - t205;
t401 = t223 * t380;
t222 = sin(qJ(2));
t366 = Icges(3,4) * t225;
t279 = -Icges(3,2) * t222 + t366;
t151 = Icges(3,6) * t223 + t226 * t279;
t367 = Icges(3,4) * t222;
t285 = Icges(3,1) * t225 - t367;
t153 = Icges(3,5) * t223 + t226 * t285;
t263 = t151 * t222 - t153 * t225;
t248 = t263 * t223;
t150 = -Icges(3,6) * t226 + t223 * t279;
t152 = -Icges(3,5) * t226 + t223 * t285;
t264 = t150 * t222 - t152 * t225;
t249 = t264 * t226;
t277 = -Icges(4,2) * t209 + t364;
t129 = Icges(4,6) * t223 + t226 * t277;
t365 = Icges(4,4) * t209;
t283 = Icges(4,1) * t210 - t365;
t133 = Icges(4,5) * t223 + t226 * t283;
t265 = t129 * t209 - t133 * t210;
t250 = t265 * t223;
t128 = -Icges(4,6) * t226 + t223 * t277;
t132 = -Icges(4,5) * t226 + t223 * t283;
t266 = t128 * t209 - t132 * t210;
t251 = t266 * t226;
t272 = Icges(5,3) * t209 + t359;
t123 = Icges(5,6) * t223 + t226 * t272;
t360 = Icges(5,5) * t209;
t281 = Icges(5,1) * t210 + t360;
t131 = Icges(5,4) * t223 + t226 * t281;
t267 = t123 * t209 + t131 * t210;
t252 = t267 * t223;
t122 = -Icges(5,6) * t226 + t223 * t272;
t130 = -Icges(5,4) * t226 + t223 * t281;
t268 = t122 * t209 + t130 * t210;
t253 = t268 * t226;
t213 = t223 * rSges(4,3);
t346 = t209 * t226;
t400 = -rSges(4,2) * t346 + t213;
t273 = Icges(4,5) * t210 - Icges(4,6) * t209;
t124 = -Icges(4,3) * t226 + t223 * t273;
t274 = Icges(3,5) * t225 - Icges(3,6) * t222;
t148 = -Icges(3,3) * t226 + t223 * t274;
t275 = Icges(5,4) * t210 + Icges(5,6) * t209;
t126 = -Icges(5,2) * t226 + t223 * t275;
t398 = 2 * m(3);
t397 = 2 * m(4);
t396 = 2 * m(5);
t395 = 2 * m(6);
t218 = t223 ^ 2;
t219 = t226 ^ 2;
t394 = m(5) / 0.2e1;
t393 = m(6) / 0.2e1;
t392 = -pkin(3) - pkin(4);
t391 = -rSges(5,1) - pkin(3);
t390 = rSges(6,3) + pkin(7);
t192 = rSges(3,1) * t222 + rSges(3,2) * t225;
t389 = m(3) * t192;
t388 = pkin(2) * t222;
t386 = pkin(4) * t209;
t385 = pkin(4) * t210;
t384 = pkin(7) * t223;
t383 = pkin(7) * t226;
t382 = t223 * pkin(6);
t216 = t226 * pkin(6);
t220 = -qJ(3) - pkin(6);
t379 = -pkin(6) - t220;
t330 = qJD(1) * t223;
t62 = t106 * t226 - t166 * t330;
t63 = -t226 * t403 - t262 * t330;
t378 = t63 * rSges(6,1) + t62 * rSges(6,2);
t377 = rSges(3,1) * t225;
t376 = rSges(4,1) * t210;
t375 = rSges(5,1) * t209;
t374 = rSges(3,2) * t222;
t373 = rSges(3,3) * t226;
t215 = t223 * rSges(5,2);
t214 = t223 * rSges(3,3);
t146 = t166 * t226;
t147 = t262 * t226;
t73 = Icges(6,5) * t147 + Icges(6,6) * t146 - Icges(6,3) * t223;
t372 = t223 * t73;
t371 = t226 * t73;
t370 = -rSges(5,3) - qJ(4);
t144 = t166 * t223;
t145 = t262 * t223;
t293 = -t145 * rSges(6,1) - t144 * rSges(6,2);
t78 = rSges(6,3) * t226 - t293;
t369 = t223 * t385 + t383 + t78;
t345 = t210 * t226;
t197 = pkin(4) * t345;
t341 = t147 * rSges(6,1) + t146 * rSges(6,2);
t79 = -rSges(6,3) * t223 + t341;
t368 = t197 - t384 + t79;
t352 = qJ(4) * t209;
t351 = qJ(4) * t210;
t344 = t220 * t226;
t120 = t216 + t344 - t401;
t199 = t226 * t205;
t121 = -pkin(1) * t226 + t223 * t379 + t199;
t343 = t223 * t120 + t226 * t121;
t198 = pkin(3) * t345;
t155 = qJ(4) * t346 + t198;
t342 = -t121 - t155;
t177 = pkin(3) * t209 - t351;
t310 = t222 * t330;
t202 = pkin(2) * t310;
t340 = t177 * t330 + t202;
t323 = qJD(2) * t226;
t307 = t210 * t323;
t322 = qJD(4) * t209;
t339 = qJ(4) * t307 + t226 * t322;
t329 = qJD(1) * t226;
t338 = rSges(5,2) * t329 + rSges(5,3) * t307;
t311 = t209 * t330;
t337 = rSges(4,2) * t311 + rSges(4,3) * t329;
t336 = t226 * t377 + t214;
t335 = t218 + t219;
t112 = rSges(6,1) * t166 - rSges(6,2) * t262;
t334 = qJD(1) * t112;
t125 = Icges(4,3) * t223 + t226 * t273;
t333 = qJD(1) * t125;
t127 = Icges(5,2) * t223 + t226 * t275;
t332 = qJD(1) * t127;
t149 = Icges(3,3) * t223 + t226 * t274;
t331 = qJD(1) * t149;
t328 = qJD(2) * t209;
t327 = qJD(2) * t210;
t326 = qJD(2) * t222;
t325 = qJD(2) * t223;
t324 = qJD(2) * t225;
t109 = Icges(6,5) * t166 - Icges(6,6) * t262;
t110 = Icges(6,4) * t166 - Icges(6,2) * t262;
t111 = Icges(6,1) * t166 - Icges(6,4) * t262;
t26 = Icges(6,4) * t63 + Icges(6,2) * t62 - Icges(6,6) * t329;
t27 = Icges(6,1) * t63 + Icges(6,4) * t62 - Icges(6,5) * t329;
t49 = Icges(6,4) * t106 + Icges(6,2) * t403;
t50 = Icges(6,1) * t106 + Icges(6,4) * t403;
t75 = Icges(6,4) * t147 + Icges(6,2) * t146 - Icges(6,6) * t223;
t77 = Icges(6,1) * t147 + Icges(6,4) * t146 - Icges(6,5) * t223;
t321 = t109 * t329 / 0.2e1 - t62 * t110 / 0.2e1 - t63 * t111 / 0.2e1 - t146 * t49 / 0.2e1 - t147 * t50 / 0.2e1 + t223 * t402 - t403 * t75 / 0.2e1 - t106 * t77 / 0.2e1 + t262 * t26 / 0.2e1 - t166 * t27 / 0.2e1;
t64 = qJD(1) * t146 + t106 * t223;
t65 = qJD(1) * t147 - t223 * t403;
t234 = Icges(6,4) * t65 + Icges(6,2) * t64 - Icges(6,6) * t330;
t235 = Icges(6,1) * t65 + Icges(6,4) * t64 - Icges(6,5) * t330;
t74 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t226;
t76 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t226;
t320 = -t109 * t330 / 0.2e1 + t64 * t110 / 0.2e1 + t65 * t111 / 0.2e1 + t144 * t49 / 0.2e1 + t145 * t50 / 0.2e1 + t226 * t402 + t74 * t409 + t76 * t408 - t262 * t234 / 0.2e1 + t166 * t235 / 0.2e1;
t319 = -t220 - t390;
t211 = qJD(3) * t223;
t315 = pkin(2) * t326;
t312 = qJD(3) * t226 + t220 * t330 + t223 * t315;
t318 = t120 * t329 + t223 * ((-t226 * t380 - t382) * qJD(1) - t312) + t226 * (-t226 * t315 + t211 + (t226 * t379 + t401) * qJD(1));
t317 = t226 * t374;
t314 = pkin(2) * t324;
t313 = t211 + t339;
t140 = rSges(5,1) * t345 + rSges(5,3) * t346 + t215;
t309 = t209 * t325;
t308 = t209 * t323;
t306 = -t177 - t388;
t179 = rSges(4,1) * t209 + rSges(4,2) * t210;
t305 = -t179 - t388;
t304 = t112 + t386;
t303 = t112 * t335;
t302 = -t223 * t220 + t199;
t292 = pkin(3) * t210 + t352;
t154 = t292 * t223;
t301 = t223 * t154 + t226 * t155 + t343;
t187 = pkin(3) * t309;
t300 = t187 + t312;
t178 = -rSges(5,3) * t210 + t375;
t299 = -t178 + t306;
t298 = -t65 * rSges(6,1) - t64 * rSges(6,2);
t22 = t146 * t75 + t147 * t77 - t372;
t233 = Icges(6,5) * t65 + Icges(6,6) * t64 - Icges(6,3) * t330;
t72 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t226;
t236 = t146 * t74 + t147 * t76 - t223 * t72;
t25 = Icges(6,5) * t63 + Icges(6,6) * t62 - Icges(6,3) * t329;
t1 = (t146 * t26 + t147 * t27 - t223 * t25 + t62 * t75 + t63 * t77) * t223 - (t146 * t234 + t147 * t235 - t223 * t233 - t329 * t72 + t62 * t74 + t63 * t76) * t226 + (t22 * t226 + (t236 - t371) * t223) * qJD(1);
t21 = t144 * t75 + t145 * t77 + t371;
t237 = t144 * t74 + t145 * t76 + t226 * t72;
t2 = (t144 * t26 + t145 * t27 + t226 * t25 + t64 * t75 + t65 * t77) * t223 - (t144 * t234 + t145 * t235 + t226 * t233 - t330 * t72 + t64 * t74 + t65 * t76) * t226 + (t21 * t226 + (t237 - t372) * t223) * qJD(1);
t297 = t223 * t1 - t226 * t2;
t296 = -t374 + t377;
t295 = -rSges(4,2) * t209 + t376;
t294 = rSges(5,1) * t210 + rSges(5,3) * t209;
t238 = t210 * t392 - t205 - t352;
t228 = t223 * t238 + t226 * t319;
t38 = t228 + t293;
t39 = t223 * t319 + t155 + t197 + t199 + t341;
t287 = t223 * t39 + t226 * t38;
t231 = t209 * t370 + t210 * t391 - t205;
t229 = t231 * t223;
t66 = (rSges(5,2) - t220) * t226 + t229;
t67 = t302 + t140 + t155;
t286 = t223 * t67 + t226 * t66;
t276 = Icges(4,2) * t210 + t365;
t261 = t154 * t329 + t223 * (qJD(1) * t198 + t223 * t322 - t187 + (t209 * t329 + t210 * t325) * qJ(4)) + t226 * (-qJ(4) * t311 + (-t210 * t330 - t308) * pkin(3) + t339) + t318;
t141 = rSges(4,1) * t345 + t400;
t137 = qJD(2) * t292 - qJD(4) * t210;
t260 = -t294 * qJD(2) - t137 - t314;
t259 = -pkin(1) - t296;
t97 = t299 * t226;
t258 = -t205 - t295;
t256 = qJD(2) * t192;
t255 = qJD(2) * t179;
t254 = -t304 + t306;
t243 = qJD(2) * t276;
t242 = qJD(2) * (-Icges(5,4) * t209 + Icges(5,6) * t210);
t241 = qJD(2) * (-Icges(3,5) * t222 - Icges(3,6) * t225);
t240 = qJD(2) * (-Icges(4,5) * t209 - Icges(4,6) * t210);
t57 = t254 * t226;
t51 = rSges(6,1) * t106 + rSges(6,2) * t403;
t232 = -t137 - t51 + (-t385 - t387) * qJD(2);
t230 = rSges(3,2) * t310 + rSges(3,3) * t329 - t226 * t256;
t183 = t296 * qJD(2);
t168 = t295 * qJD(2);
t157 = -t317 + t336;
t156 = t223 * t296 - t373;
t139 = -rSges(4,3) * t226 + t223 * t295;
t138 = -rSges(5,2) * t226 + t223 * t294;
t119 = t305 * t226;
t118 = t305 * t223;
t117 = t382 + (pkin(1) - t374) * t226 + t336;
t116 = t223 * t259 + t216 + t373;
t108 = t141 + t302;
t107 = (rSges(4,3) - t220) * t226 + t258 * t223;
t100 = t223 * t241 + t331;
t99 = -qJD(1) * t148 + t226 * t241;
t96 = t299 * t223;
t87 = t223 * t242 + t332;
t86 = -qJD(1) * t126 + t226 * t242;
t85 = t223 * t240 + t333;
t84 = -qJD(1) * t124 + t226 * t240;
t81 = t192 * t325 + ((-rSges(3,3) - pkin(6)) * t223 + t259 * t226) * qJD(1);
t80 = (t216 + (-pkin(1) - t377) * t223) * qJD(1) + t230;
t71 = -t179 * t329 - t223 * t168 + (-t222 * t329 - t223 * t324) * pkin(2);
t70 = t179 * t330 + t202 + (-t168 - t314) * t226;
t61 = t223 * t149 - t226 * t263;
t60 = t223 * t148 - t249;
t59 = -t149 * t226 - t248;
t58 = -t148 * t226 - t223 * t264;
t56 = t254 * t223;
t53 = t179 * t325 + (t226 * t258 - t213) * qJD(1) + t312;
t52 = t211 + (-t205 - t376) * t330 + (-qJD(1) * t220 + qJD(2) * t305) * t226 + t337;
t47 = t223 * t125 - t226 * t265;
t46 = t223 * t124 - t251;
t45 = t223 * t127 + t226 * t267;
t44 = t223 * t126 + t253;
t43 = -t125 * t226 - t250;
t42 = -t124 * t226 - t223 * t266;
t41 = -t127 * t226 + t252;
t40 = -t126 * t226 + t223 * t268;
t37 = qJD(1) * t97 + t223 * t260;
t36 = t178 * t330 + t226 * t260 + t340;
t35 = -t223 * t78 - t226 * t79;
t34 = (-t322 + (t210 * t370 + t375) * qJD(2)) * t223 + (t226 * t231 - t215) * qJD(1) + t300;
t33 = (t209 * t391 - t388) * t323 + (t229 - t344) * qJD(1) + t313 + t338;
t32 = t166 * t77 - t262 * t75;
t31 = t166 * t76 - t262 * t74;
t30 = t223 * t138 + t140 * t226 + t301;
t29 = -rSges(6,3) * t330 - t298;
t28 = -rSges(6,3) * t329 + t378;
t24 = -t109 * t223 + t110 * t146 + t111 * t147;
t23 = t109 * t226 + t144 * t110 + t145 * t111;
t20 = qJD(1) * t57 + t223 * t232;
t19 = t226 * t232 + t304 * t330 + t340;
t16 = t223 * t369 + t226 * t368 + t301;
t15 = (-t322 + (-t351 + t386) * qJD(2)) * t223 + (t223 * t390 + t226 * t238) * qJD(1) + t298 + t300;
t14 = (t209 * t392 - t388) * t323 + t228 * qJD(1) + t313 + t378;
t13 = t22 * t223 - t226 * t236;
t12 = t21 * t223 - t226 * t237;
t11 = -t223 * t29 - t226 * t28 + (t223 * t79 - t226 * t78) * qJD(1);
t10 = t226 * t338 + (-t178 * t218 - t219 * t375) * qJD(2) + (t226 * t138 + (-t140 + t342 + t215) * t223) * qJD(1) + t261;
t3 = (-pkin(4) * t308 + t28 + (t369 - t383) * qJD(1)) * t226 + (-pkin(4) * t309 + t29 + (t342 - t368 - t384) * qJD(1)) * t223 + t261;
t4 = [t106 * t111 + t166 * t50 + t403 * t110 - t262 * t49 + (t14 * t39 + t15 * t38) * t395 + (t33 * t67 + t34 * t66) * t396 + (t107 * t53 + t108 * t52) * t397 + (t116 * t81 + t117 * t80) * t398 + (-Icges(3,2) * t225 + t285 - t367) * t326 + (Icges(3,1) * t222 + t279 + t366) * t324 + (-Icges(5,3) * t210 - t276 + t281 + t283 + t360) * t328 + (t277 - t272 + t411) * t327; m(4) * (t107 * t70 + t108 * t71 + t118 * t52 + t119 * t53) + m(5) * (t33 * t96 + t34 * t97 + t36 * t66 + t37 * t67) + m(6) * (t14 * t56 + t15 * t57 + t19 * t38 + t20 * t39) + ((t123 * t404 + t129 * t405 + t243 * t407) * t210 + (t410 * t407 + (t131 + t133) * t405) * t209 - t320 + m(3) * (-t116 * t183 - t192 * t81)) * t226 + ((t122 * t404 + t128 * t405 + t243 * t406) * t210 + (t410 * t406 + (t130 + t132) * t405) * t209 - t321 + m(3) * (-t117 * t183 - t192 * t80)) * t223 + (t249 / 0.2e1 - t253 / 0.2e1 + t251 / 0.2e1 - t248 / 0.2e1 + t252 / 0.2e1 - t250 / 0.2e1 + (t274 + t273 + t275) * (t218 / 0.2e1 + t219 / 0.2e1)) * qJD(2) + ((-t117 * t389 + t24 / 0.2e1 + t32 / 0.2e1 + (-t123 / 0.2e1 + t129 / 0.2e1) * t210 + (t131 / 0.2e1 + t133 / 0.2e1) * t209) * t226 + (t116 * t389 + t23 / 0.2e1 + t31 / 0.2e1 + (-t122 / 0.2e1 + t128 / 0.2e1) * t210 + (t130 / 0.2e1 + t132 / 0.2e1) * t209) * t223) * qJD(1); (t16 * t3 + t19 * t57 + t20 * t56) * t395 + (t10 * t30 + t36 * t97 + t37 * t96) * t396 + (t119 * t70 + t118 * t71 + (t223 * t139 + t141 * t226 + t343) * ((qJD(1) * t139 - t226 * t255 + t337) * t226 + (-t223 * t255 + (-t121 - t141 + t400) * qJD(1)) * t223 + t318)) * t397 - t226 * ((t226 * t85 + (t43 + t251) * qJD(1)) * t226 + (t42 * qJD(1) + (-t129 * t327 - t133 * t328 + t333) * t223 + (-t84 + (t128 * t210 + t132 * t209) * qJD(2) - t265 * qJD(1)) * t226) * t223) - t226 * ((t226 * t87 + (t41 - t253) * qJD(1)) * t226 + (t40 * qJD(1) + (t123 * t327 - t131 * t328 + t332) * t223 + (-t86 + (-t122 * t210 + t130 * t209) * qJD(2) + t267 * qJD(1)) * t226) * t223) - t226 * ((t226 * t100 + (t59 + t249) * qJD(1)) * t226 + (t58 * qJD(1) + (-t151 * t324 - t153 * t326 + t331) * t223 + (-t99 + (t150 * t225 + t152 * t222) * qJD(2) - t263 * qJD(1)) * t226) * t223) + t223 * ((t223 * t84 + (t46 + t250) * qJD(1)) * t223 + (t47 * qJD(1) + (t128 * t327 + t132 * t328) * t226 + (-t85 + (-t129 * t210 - t133 * t209) * qJD(2) + (t125 - t266) * qJD(1)) * t223) * t226) + t223 * ((t223 * t86 + (t44 - t252) * qJD(1)) * t223 + (t45 * qJD(1) + (-t122 * t327 + t130 * t328) * t226 + (-t87 + (t123 * t210 - t131 * t209) * qJD(2) + (t127 + t268) * qJD(1)) * t223) * t226) + ((t223 * t156 + t157 * t226) * ((qJD(1) * t156 + t230) * t226 + (-t223 * t256 + (-t157 - t317 + t214) * qJD(1)) * t223) + t335 * t192 * t183) * t398 + t223 * ((t223 * t99 + (t60 + t248) * qJD(1)) * t223 + (t61 * qJD(1) + (t150 * t324 + t152 * t326) * t226 + (-t100 + (-t151 * t225 - t153 * t222) * qJD(2) + (t149 - t264) * qJD(1)) * t223) * t226) + t297 + (t12 + (-t40 - t42 - t58) * t226 + (t41 + t43 + t59) * t223) * t330 + (t13 + (-t44 - t46 - t60) * t226 + (t45 + t47 + t61) * t223) * t329; m(6) * (qJD(1) * t287 - t14 * t226 + t223 * t15) + m(5) * (qJD(1) * t286 + t223 * t34 - t226 * t33) + m(4) * (t223 * t53 - t226 * t52 + (t107 * t226 + t108 * t223) * qJD(1)); m(6) * (t223 * t19 - t20 * t226 + (t223 * t56 + t226 * t57) * qJD(1)) + m(5) * (t223 * t36 - t226 * t37 + (t223 * t96 + t226 * t97) * qJD(1)) + m(4) * (t223 * t70 - t226 * t71 + (t118 * t223 + t119 * t226) * qJD(1)); 0; 0.2e1 * (t286 * t394 + t287 * t393) * t327 + 0.2e1 * ((t14 * t223 + t15 * t226 + t329 * t39 - t330 * t38) * t393 + (t223 * t33 + t226 * t34 + t329 * t67 - t330 * t66) * t394) * t209; 0.2e1 * ((t323 * t57 + t325 * t56 - t3) * t393 + (t323 * t97 + t325 * t96 - t10) * t394) * t210 + 0.2e1 * ((qJD(2) * t16 + t19 * t226 + t20 * t223 + t329 * t56 - t330 * t57) * t393 + (qJD(2) * t30 + t223 * t37 + t226 * t36 + t329 * t96 - t330 * t97) * t394) * t209; 0; 0.4e1 * (t394 + t393) * (-0.1e1 + t335) * t209 * t327; (m(6) * (t112 * t15 + t334 * t39 + t38 * t51) + t320) * t226 + (m(6) * (t112 * t14 - t334 * t38 + t39 * t51) + t321) * t223 + ((t24 + t32) * t226 + (t23 + t31) * t223) * t405; m(6) * (t11 * t16 + t35 * t3) + (m(6) * (t112 * t19 + t334 * t56 + t51 * t57) + t2 - qJD(1) * t13) * t226 + (m(6) * (t112 * t20 - t334 * t57 + t51 * t56) - qJD(1) * t12 - t1) * t223; 0; m(6) * (-t11 * t210 + t335 * t51 * t209 + (t209 * t35 + t210 * t303) * qJD(2)); (t11 * t35 + t303 * t51) * t395 + (t223 * t12 + t226 * t13) * qJD(1) + t297;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
