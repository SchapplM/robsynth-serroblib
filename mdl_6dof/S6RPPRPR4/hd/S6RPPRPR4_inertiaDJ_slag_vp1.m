% Calculate time derivative of joint inertia matrix for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:16
% EndTime: 2019-03-09 01:46:31
% DurationCPUTime: 9.41s
% Computational Cost: add. (15184->599), mult. (27163->844), div. (0->0), fcn. (30219->10), ass. (0->280)
t398 = Icges(5,3) + Icges(6,3);
t217 = qJ(4) + pkin(10);
t209 = sin(t217);
t210 = cos(t217);
t220 = sin(qJ(4));
t222 = cos(qJ(4));
t399 = -Icges(5,5) * t222 - Icges(6,5) * t210 + Icges(5,6) * t220 + Icges(6,6) * t209;
t347 = sin(pkin(9));
t348 = cos(pkin(9));
t366 = sin(qJ(1));
t367 = cos(qJ(1));
t191 = -t347 * t366 - t348 * t367;
t192 = t367 * t347 - t366 * t348;
t382 = -t191 * t398 + t399 * t192;
t397 = t399 * t191 + t192 * t398;
t343 = Icges(6,4) * t210;
t260 = Icges(6,2) * t209 - t343;
t112 = Icges(6,6) * t192 + t191 * t260;
t344 = Icges(6,4) * t209;
t263 = -Icges(6,1) * t210 + t344;
t114 = Icges(6,5) * t192 + t191 * t263;
t345 = Icges(5,4) * t222;
t261 = Icges(5,2) * t220 - t345;
t123 = Icges(5,6) * t192 + t191 * t261;
t346 = Icges(5,4) * t220;
t264 = -Icges(5,1) * t222 + t346;
t125 = Icges(5,5) * t192 + t191 * t264;
t380 = t112 * t209 - t114 * t210 + t123 * t220 - t125 * t222;
t111 = -Icges(6,6) * t191 + t192 * t260;
t113 = -Icges(6,5) * t191 + t192 * t263;
t122 = -Icges(5,6) * t191 + t192 * t261;
t124 = -Icges(5,5) * t191 + t192 * t264;
t379 = t111 * t209 - t113 * t210 + t122 * t220 - t124 * t222;
t316 = t367 * pkin(1) + t366 * qJ(2);
t301 = t367 * pkin(2) + t316;
t278 = rSges(6,1) * t209 + rSges(6,2) * t210;
t363 = pkin(4) * t220;
t396 = t278 + t363;
t180 = t191 * qJD(1);
t310 = qJD(4) * t222;
t243 = t180 * t220 + t192 * t310;
t181 = t192 * qJD(1);
t311 = qJD(4) * t220;
t298 = t191 * t311;
t238 = t181 * t222 + t298;
t312 = qJD(4) * t210;
t336 = t181 * t209;
t241 = t191 * t312 - t336;
t395 = -rSges(3,1) * t367 - rSges(3,3) * t366 - t316;
t392 = -0.2e1 * t192;
t239 = -t181 * t220 + t191 * t310;
t313 = qJD(4) * t209;
t300 = t191 * t313;
t335 = t181 * t210;
t240 = t300 + t335;
t391 = Icges(5,5) * t238 + Icges(6,5) * t240 + Icges(5,6) * t239 + Icges(6,6) * t241 + t180 * t398;
t297 = t192 * t311;
t242 = -t180 * t222 + t297;
t244 = -t180 * t210 + t192 * t313;
t339 = t180 * t209;
t245 = t192 * t312 + t339;
t390 = -Icges(5,5) * t242 - Icges(6,5) * t244 - Icges(5,6) * t243 - Icges(6,6) * t245 - t181 * t398;
t185 = t191 * pkin(7);
t218 = -qJ(5) - pkin(7);
t330 = t191 * t218;
t208 = pkin(4) * t222 + pkin(3);
t359 = pkin(3) - t208;
t102 = t192 * t359 + t185 + t330;
t279 = rSges(6,1) * t210 - rSges(6,2) * t209;
t386 = t192 * t279;
t115 = -t191 * rSges(6,3) - t386;
t331 = t191 * t210;
t332 = t191 * t209;
t116 = -rSges(6,1) * t331 + rSges(6,2) * t332 + t192 * rSges(6,3);
t237 = qJD(4) * t278;
t251 = rSges(6,1) * t335 - rSges(6,2) * t336 + t180 * rSges(6,3);
t322 = -t191 * t208 - t192 * t218;
t360 = pkin(7) * t192;
t103 = pkin(3) * t191 + t322 - t360;
t281 = t192 * qJD(5) - t180 * t218 + t181 * t208;
t286 = t181 * pkin(3) + t180 * pkin(7);
t358 = t191 * (pkin(4) * t298 + t281 - t286) - t181 * t103;
t172 = t181 * pkin(7);
t320 = pkin(4) * t297 - t191 * qJD(5);
t334 = t181 * t218;
t63 = t180 * t359 - t172 + t320 - t334;
t17 = -t181 * t116 + (t191 * t237 + t251) * t191 + (t181 * rSges(6,3) + t192 * t237 + t63) * t192 - (-t102 - t115 + t386) * t180 + t358;
t375 = 2 * m(6);
t389 = t17 * t375;
t356 = rSges(5,2) * t220;
t357 = rSges(5,1) * t222;
t280 = t356 - t357;
t352 = t191 * rSges(5,3);
t127 = t192 * t280 - t352;
t319 = t192 * rSges(5,3) - t191 * t357;
t128 = t191 * t356 + t319;
t236 = rSges(5,1) * t238 + t180 * rSges(5,3);
t353 = t181 * rSges(5,3);
t30 = t180 * t127 + t192 * (rSges(5,1) * t242 + t353) - t181 * t128 + t191 * t236 + (t191 * t239 + t192 * t243) * rSges(5,2);
t376 = 2 * m(5);
t388 = t30 * t376;
t219 = sin(qJ(6));
t221 = cos(qJ(6));
t228 = qJD(6) * t192 + t240;
t308 = qJD(6) * t210;
t283 = t191 * t308 + t180;
t68 = -t219 * t228 + t221 * t283;
t69 = t219 * t283 + t221 * t228;
t40 = t69 * rSges(7,1) + t68 * rSges(7,2) - rSges(7,3) * t241;
t387 = pkin(5) * t335 - pkin(8) * t241 + t40;
t361 = pkin(5) * t210;
t285 = pkin(8) * t209 + t361;
t385 = t192 * t285;
t202 = -t220 * rSges(5,1) - rSges(5,2) * t222;
t384 = t202 ^ 2 * t376;
t197 = t280 * qJD(4);
t383 = t197 * t202 * t376;
t378 = (t380 + t382) * t192 + (t379 - t397) * t191;
t368 = rSges(7,3) + pkin(8);
t377 = t209 * t368 + t208 + t361;
t374 = 2 * m(7);
t373 = t181 / 0.2e1;
t372 = -t191 / 0.2e1;
t371 = -t192 / 0.2e1;
t370 = t192 / 0.2e1;
t369 = t210 / 0.2e1;
t365 = m(5) * t197;
t364 = m(5) * t202;
t362 = pkin(5) * t209;
t326 = t210 * t219;
t137 = -t191 * t221 + t192 * t326;
t325 = t210 * t221;
t138 = -t191 * t219 - t192 * t325;
t329 = t192 * t209;
t80 = Icges(7,4) * t138 + Icges(7,2) * t137 - Icges(7,6) * t329;
t82 = Icges(7,1) * t138 + Icges(7,4) * t137 - Icges(7,5) * t329;
t271 = t219 * t80 - t221 * t82;
t229 = -qJD(6) * t191 + t244;
t282 = t192 * t308 + t181;
t66 = -t219 * t229 + t221 * t282;
t67 = t219 * t282 + t221 * t229;
t33 = Icges(7,5) * t67 + Icges(7,6) * t66 - Icges(7,3) * t245;
t35 = Icges(7,4) * t67 + Icges(7,2) * t66 - Icges(7,6) * t245;
t37 = Icges(7,1) * t67 + Icges(7,4) * t66 - Icges(7,5) * t245;
t78 = Icges(7,5) * t138 + Icges(7,6) * t137 - Icges(7,3) * t329;
t10 = (qJD(4) * t271 + t33) * t210 + (-qJD(4) * t78 + t219 * t35 - t221 * t37 + (t219 * t82 + t221 * t80) * qJD(6)) * t209;
t355 = t10 * t191;
t139 = t191 * t326 + t192 * t221;
t140 = -t191 * t325 + t192 * t219;
t81 = Icges(7,4) * t140 + Icges(7,2) * t139 - Icges(7,6) * t332;
t83 = Icges(7,1) * t140 + Icges(7,4) * t139 - Icges(7,5) * t332;
t270 = t219 * t81 - t221 * t83;
t34 = Icges(7,5) * t69 + Icges(7,6) * t68 - Icges(7,3) * t241;
t36 = Icges(7,4) * t69 + Icges(7,2) * t68 - Icges(7,6) * t241;
t38 = Icges(7,1) * t69 + Icges(7,4) * t68 - Icges(7,5) * t241;
t79 = Icges(7,5) * t140 + Icges(7,6) * t139 - Icges(7,3) * t332;
t11 = (qJD(4) * t270 + t34) * t210 + (-qJD(4) * t79 + t219 * t36 - t221 * t38 + (t219 * t83 + t221 * t81) * qJD(6)) * t209;
t354 = t11 * t192;
t351 = rSges(6,3) - t218;
t277 = -t138 * rSges(7,1) - t137 * rSges(7,2);
t85 = -rSges(7,3) * t329 - t277;
t350 = t85 - t385;
t86 = t140 * rSges(7,1) + t139 * rSges(7,2) - rSges(7,3) * t332;
t349 = pkin(5) * t331 + pkin(8) * t332 - t86;
t342 = Icges(7,4) * t219;
t341 = Icges(7,4) * t221;
t259 = Icges(7,2) * t219 - t341;
t147 = Icges(7,6) * t210 + t209 * t259;
t340 = t147 * t219;
t256 = -Icges(7,5) * t221 + Icges(7,6) * t219;
t146 = Icges(7,3) * t210 + t209 * t256;
t327 = t209 * t146;
t276 = -rSges(7,1) * t221 + rSges(7,2) * t219;
t309 = qJD(6) * t209;
t126 = (rSges(7,1) * t219 + rSges(7,2) * t221) * t309 + (-rSges(7,3) * t209 + t210 * t276) * qJD(4);
t324 = t285 * qJD(4) - t126;
t149 = rSges(7,3) * t210 + t209 * t276;
t284 = -pkin(8) * t210 + t362;
t323 = t149 - t284;
t321 = t243 * pkin(4);
t318 = Icges(5,5) * t220 + Icges(6,5) * t209 + Icges(5,6) * t222 + Icges(6,6) * t210;
t214 = t367 * qJ(2);
t317 = qJD(1) * t214 + qJD(2) * t366;
t315 = qJD(4) * t191;
t314 = qJD(4) * t192;
t307 = m(6) / 0.2e1 + m(7) / 0.2e1;
t306 = t366 * pkin(1);
t31 = t209 * t271 + t210 * t78;
t262 = -Icges(7,1) * t221 + t342;
t148 = Icges(7,5) * t210 + t209 * t262;
t56 = t137 * t147 + t138 * t148 - t192 * t327;
t305 = t31 / 0.2e1 + t56 / 0.2e1;
t32 = t209 * t270 + t210 * t79;
t57 = t139 * t147 + t140 * t148 - t191 * t327;
t304 = t57 / 0.2e1 + t32 / 0.2e1;
t303 = t148 * t325;
t290 = -t323 + t363;
t287 = t67 * rSges(7,1) + t66 * rSges(7,2);
t25 = t137 * t80 + t138 * t82 - t329 * t78;
t26 = t137 * t81 + t138 * t83 - t329 * t79;
t275 = -t191 * t26 - t192 * t25;
t27 = t139 * t80 + t140 * t82 - t332 * t78;
t28 = t139 * t81 + t140 * t83 - t332 * t79;
t274 = -t191 * t28 - t192 * t27;
t273 = -t191 * t32 - t192 * t31;
t272 = -t191 * t85 + t192 * t86;
t269 = t301 + t322;
t250 = pkin(3) - t280;
t249 = t208 + t279;
t247 = -t209 * t33 - t312 * t78;
t246 = -t209 * t34 - t312 * t79;
t235 = -pkin(2) * t366 - t306;
t232 = -t191 * t366 + t192 * t367;
t231 = t214 + t235;
t117 = (Icges(7,5) * t219 + Icges(7,6) * t221) * t309 + (-Icges(7,3) * t209 + t210 * t256) * qJD(4);
t118 = (Icges(7,2) * t221 + t342) * t309 + (-Icges(7,6) * t209 + t210 * t259) * qJD(4);
t119 = (Icges(7,1) * t219 + t341) * t309 + (-Icges(7,5) * t209 + t210 * t262) * qJD(4);
t230 = t210 * t117 + t312 * t340 + (t147 * t221 + t148 * t219) * t309 + (t118 * t219 - t119 * t221) * t209;
t227 = -rSges(3,1) * t366 + rSges(3,3) * t367 - t306;
t226 = qJD(1) * t235 + t317;
t212 = qJD(2) * t367;
t225 = -qJD(1) * t301 + t212;
t224 = t225 - t320;
t223 = t226 + t281;
t177 = t279 * qJD(4);
t173 = t192 * t363;
t167 = t214 + t227;
t145 = qJD(1) * t395 + t212;
t144 = qJD(1) * t227 + t317;
t132 = -rSges(4,1) * t191 - rSges(4,2) * t192 + t301;
t131 = t192 * rSges(4,1) - t191 * rSges(4,2) + t231;
t130 = t396 * t191;
t129 = t192 * t278 + t173;
t107 = t180 * rSges(4,1) + t181 * rSges(4,2) + t225;
t106 = t181 * rSges(4,1) - t180 * rSges(4,2) + t226;
t100 = t192 * t102;
t98 = t290 * t191;
t97 = -t192 * t323 + t173;
t96 = t360 + (-pkin(3) + t356) * t191 + t301 + t319;
t95 = t192 * t250 + t185 + t231 + t352;
t88 = pkin(4) * t239 + t191 * t177 - t181 * t278;
t87 = t177 * t192 + t180 * t278 + t321;
t84 = t146 * t210 + (-t148 * t221 + t340) * t209;
t77 = t116 + t269;
t76 = t191 * t351 + t192 * t249 + t231;
t61 = t149 * t332 + t210 * t86;
t60 = -t149 * t329 - t210 * t85;
t59 = t180 * t250 + t202 * t314 - t172 + t225 - t353;
t58 = rSges(5,2) * t239 + t226 + t236 + t286;
t55 = t269 - t349;
t54 = t192 * t377 + t231 + t277 - t330;
t49 = t180 * t249 - t181 * t351 - t278 * t314 + t224;
t48 = t315 * t396 + t223 + t251;
t47 = (pkin(4) * t310 + t324) * t191 - t290 * t181;
t46 = -t180 * t323 + t192 * t324 + t321;
t41 = t272 * t209;
t39 = -rSges(7,3) * t245 + t287;
t29 = ((-t303 - t327) * qJD(4) + t230) * t210;
t24 = t334 + (t210 * t368 - t362) * t314 + t377 * t180 + t224 - t287;
t23 = (t362 + t363) * t315 + t223 + t387;
t22 = t100 + t350 * t192 + (t103 - t349) * t191;
t21 = (t149 * t315 + t40) * t210 + (-qJD(4) * t86 + t191 * t126 - t181 * t149) * t209;
t20 = (-t149 * t314 - t39) * t210 + (qJD(4) * t85 - t126 * t192 - t149 * t180) * t209;
t19 = -t117 * t332 + t118 * t139 + t119 * t140 - t146 * t241 + t147 * t68 + t148 * t69;
t18 = -t117 * t329 + t118 * t137 + t119 * t138 - t146 * t245 + t147 * t66 + t148 * t67;
t16 = -t191 * t27 + t192 * t28;
t15 = -t191 * t25 + t192 * t26;
t14 = t209 * t274 + t210 * t57;
t13 = t209 * t275 + t210 * t56;
t12 = t272 * t312 + (t180 * t86 + t181 * t85 - t191 * t39 + t192 * t40) * t209;
t9 = (pkin(5) * t300 + t387) * t191 + t349 * t181 + (t284 * t314 + t39 + t63) * t192 - (-t102 - t350 + t385) * t180 + t358;
t8 = t139 * t36 + t140 * t38 + t191 * t246 + t336 * t79 + t68 * t81 + t69 * t83;
t7 = t139 * t35 + t140 * t37 + t191 * t247 + t336 * t78 + t68 * t80 + t69 * t82;
t6 = t137 * t36 + t138 * t38 + t192 * t246 - t339 * t79 + t66 * t81 + t67 * t83;
t5 = t137 * t35 + t138 * t37 + t192 * t247 - t339 * t78 + t66 * t80 + t67 * t82;
t4 = t180 * t28 + t181 * t27 - t191 * t7 + t192 * t8;
t3 = t180 * t26 + t181 * t25 - t191 * t5 + t192 * t6;
t2 = (qJD(4) * t274 + t19) * t210 + (-qJD(4) * t57 - t180 * t27 + t181 * t28 - t191 * t8 - t192 * t7) * t209;
t1 = (qJD(4) * t275 + t18) * t210 + (-qJD(4) * t56 - t180 * t25 + t181 * t26 - t191 * t6 - t192 * t5) * t209;
t42 = [t230 + (t48 * t77 + t49 * t76) * t375 + 0.2e1 * m(4) * (t106 * t132 + t107 * t131) + 0.2e1 * m(3) * (-t144 * t395 + t145 * t167) + (t23 * t55 + t24 * t54) * t374 + (t58 * t96 + t59 * t95) * t376 - qJD(4) * t303 + (Icges(6,1) * t209 - t260 + t343) * t312 + (-Icges(5,2) * t222 - t264 - t346) * t311 + (Icges(5,1) * t220 - t261 + t345) * t310 + (-Icges(6,2) * t210 - t146 - t263 - t344) * t313; m(7) * (t366 * t24 - t367 * t23 + (t366 * t55 + t367 * t54) * qJD(1)) + m(5) * (t366 * t59 - t367 * t58 + (t366 * t96 + t367 * t95) * qJD(1)) + m(6) * (t366 * t49 - t367 * t48 + (t366 * t77 + t367 * t76) * qJD(1)) + m(4) * (t366 * t107 - t367 * t106 + (t131 * t367 + t132 * t366) * qJD(1)) + m(3) * (t366 * t145 - t367 * t144 + (t167 * t367 - t366 * t395) * qJD(1)); 0; 0; 0; 0; t354 / 0.2e1 - t355 / 0.2e1 + (-t122 * t222 / 0.2e1 - t220 * t124 / 0.2e1 - t111 * t210 / 0.2e1 - t113 * t209 / 0.2e1 + t318 * t191 + t305) * t181 - (t123 * t222 / 0.2e1 + t220 * t125 / 0.2e1 + t112 * t369 + t114 * t209 / 0.2e1 + t318 * t192 - t304) * t180 + m(7) * (t23 * t97 + t24 * t98 + t46 * t55 + t47 * t54) + m(6) * (t129 * t48 + t130 * t49 + t76 * t88 + t77 * t87) + (-t191 * t95 - t192 * t96) * t365 + (-t180 * t96 + t181 * t95 - t191 * t59 - t192 * t58) * t364 + t399 * qJD(4) * (t191 ^ 2 / 0.2e1 + t192 ^ 2 / 0.2e1) + (t379 * qJD(4) - t209 * (Icges(6,1) * t244 + Icges(6,4) * t245 + Icges(6,5) * t181) - t210 * (Icges(6,4) * t244 + Icges(6,2) * t245 + Icges(6,6) * t181) - t220 * (Icges(5,1) * t242 + Icges(5,4) * t243 + Icges(5,5) * t181) - t222 * (Icges(5,4) * t242 + Icges(5,2) * t243 + Icges(5,6) * t181) + t18) * t372 + (t380 * qJD(4) - t209 * (Icges(6,1) * t240 + Icges(6,4) * t241 + Icges(6,5) * t180) - t210 * (Icges(6,4) * t240 + Icges(6,2) * t241 + Icges(6,6) * t180) - t220 * (Icges(5,1) * t238 + Icges(5,4) * t239 + Icges(5,5) * t180) - t222 * (Icges(5,4) * t238 + Icges(5,2) * t239 + Icges(5,6) * t180) + t19) * t370; m(6) * (t88 * t366 - t87 * t367 + (t129 * t366 + t130 * t367) * qJD(1)) + m(7) * (t47 * t366 - t46 * t367 + (t366 * t97 + t367 * t98) * qJD(1)) + t232 * t365 + (t366 * t181 + t367 * t180 + (-t191 * t367 - t192 * t366) * qJD(1)) * t364; -m(5) * t30 - m(6) * t17 - m(7) * t9; t180 * t16 + t181 * t15 + (t22 * t9 + t46 * t97 + t47 * t98) * t374 + (t100 * t17 + t129 * t87 + t130 * t88) * t375 + (-t3 + (t103 + t116) * t389 + t128 * t388 + (t191 * t390 + t383) * t191 + (0.3e1 * t382 * t191 + t379 * t392 - t384) * t181 - (-t191 * t379 + t378) * t180) * t191 + (t4 + t115 * t389 + t127 * t388 + (t192 * t391 + t383) * t192 - (-0.2e1 * t191 * t380 - t384 + (-t192 + t392) * t397) * t180 + (-t380 * t192 + t378) * t181 + (t390 * t192 + t391 * t191 - t397 * t181 - t382 * t180 + (-t124 * t181 + t125 * t180) * t222 + (t122 * t181 - t123 * t180) * t220 + (-t113 * t181 + t114 * t180) * t210 + (t111 * t181 - t112 * t180) * t209) * t191) * t192; m(7) * (t180 * t54 + t181 * t55 - t191 * t23 + t192 * t24) + m(6) * (t180 * t76 + t181 * t77 - t191 * t48 + t192 * t49); 0.2e1 * t307 * (qJD(1) * t232 + t180 * t366 - t181 * t367); 0; m(7) * (t180 * t98 + t181 * t97 - t191 * t46 + t192 * t47) + m(6) * (t129 * t181 + t130 * t180 - t191 * t87 + t192 * t88); 0.4e1 * t307 * (t180 * t192 - t181 * t191); m(7) * (t20 * t54 + t21 * t55 + t23 * t61 + t24 * t60) + t29 + (-t191 * t304 - t192 * t305) * t312 + (-t84 * qJD(4) + (-t10 / 0.2e1 - t18 / 0.2e1) * t192 + (-t19 / 0.2e1 - t11 / 0.2e1) * t191 + t304 * t181 - t305 * t180) * t209; m(7) * (t20 * t366 - t21 * t367 + (t366 * t61 + t367 * t60) * qJD(1)); -m(7) * t12; t180 * t14 / 0.2e1 + t2 * t370 + t13 * t373 + t1 * t372 + (t32 * t180 + t31 * t181 + t354 - t355) * t369 + m(7) * (t12 * t22 + t20 * t98 + t21 * t97 + t41 * t9 + t46 * t61 + t47 * t60) + (t15 * t371 + t16 * t372) * t312 + (t16 * t373 + t4 * t372 - t180 * t15 / 0.2e1 + t3 * t371 - qJD(4) * (-t191 * t31 + t192 * t32) / 0.2e1) * t209; m(7) * (t180 * t60 + t181 * t61 - t191 * t21 + t192 * t20); (t12 * t41 + t20 * t60 + t21 * t61) * t374 + (t29 + (-t192 * t13 - t191 * t14 + t210 * t273) * qJD(4)) * t210 + (t181 * t14 - t191 * t2 - t180 * t13 - t192 * t1 + t210 * (-t10 * t192 - t11 * t191 - t31 * t180 + t32 * t181) + (-t209 * t273 - 0.2e1 * t210 * t84) * qJD(4)) * t209;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t42(1) t42(2) t42(4) t42(7) t42(11) t42(16); t42(2) t42(3) t42(5) t42(8) t42(12) t42(17); t42(4) t42(5) t42(6) t42(9) t42(13) t42(18); t42(7) t42(8) t42(9) t42(10) t42(14) t42(19); t42(11) t42(12) t42(13) t42(14) t42(15) t42(20); t42(16) t42(17) t42(18) t42(19) t42(20) t42(21);];
Mq  = res;
