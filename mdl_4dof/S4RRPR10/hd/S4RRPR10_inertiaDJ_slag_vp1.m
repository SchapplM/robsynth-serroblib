% Calculate time derivative of joint inertia matrix for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR10_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR10_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:10
% EndTime: 2019-12-31 17:11:25
% DurationCPUTime: 9.71s
% Computational Cost: add. (5279->540), mult. (14364->791), div. (0->0), fcn. (13443->6), ass. (0->268)
t193 = cos(qJ(2));
t190 = sin(qJ(2));
t303 = Icges(4,6) * t190;
t311 = Icges(3,4) * t190;
t355 = -t303 - t311 + (-Icges(3,2) - Icges(4,3)) * t193;
t302 = Icges(4,6) * t193;
t310 = Icges(3,4) * t193;
t354 = -t302 - t310 + (-Icges(3,1) - Icges(4,2)) * t190;
t353 = t355 * qJD(2);
t352 = t354 * qJD(2);
t191 = sin(qJ(1));
t194 = cos(qJ(1));
t226 = -Icges(3,2) * t190 + t310;
t116 = -Icges(3,6) * t194 + t191 * t226;
t218 = -Icges(4,3) * t190 + t302;
t338 = Icges(4,5) * t194 + t191 * t218;
t351 = -t116 - t338;
t117 = Icges(3,6) * t191 + t194 * t226;
t121 = Icges(4,5) * t191 - t194 * t218;
t350 = t117 - t121;
t229 = Icges(3,1) * t193 - t311;
t119 = -Icges(3,5) * t194 + t191 * t229;
t220 = Icges(4,2) * t193 - t303;
t337 = Icges(4,4) * t194 + t191 * t220;
t349 = t119 + t337;
t120 = Icges(3,5) * t191 + t194 * t229;
t123 = Icges(4,4) * t191 - t194 * t220;
t348 = t120 - t123;
t276 = qJD(2) * t194;
t260 = t190 * t276;
t281 = qJD(1) * t191;
t347 = t193 * t281 + t260;
t192 = cos(qJ(4));
t189 = sin(qJ(4));
t308 = Icges(5,4) * t189;
t223 = Icges(5,2) * t192 + t308;
t115 = Icges(5,6) * t190 - t193 * t223;
t307 = Icges(5,4) * t192;
t227 = Icges(5,1) * t189 + t307;
t118 = Icges(5,5) * t190 - t193 * t227;
t346 = t115 * t192 + t118 * t189;
t345 = qJD(2) / 0.2e1;
t214 = t121 * t190 - t123 * t193;
t344 = t191 * t214;
t215 = t117 * t190 - t120 * t193;
t343 = t191 * t215;
t213 = -t190 * t338 + t193 * t337;
t342 = t194 * t213;
t216 = t116 * t190 - t119 * t193;
t341 = t194 * t216;
t291 = t193 * t194;
t340 = t191 * rSges(4,1) - rSges(4,2) * t291;
t296 = t190 * t194;
t339 = -rSges(3,2) * t296 + t191 * rSges(3,3);
t222 = Icges(3,5) * t193 - Icges(3,6) * t190;
t113 = -Icges(3,3) * t194 + t191 * t222;
t224 = Icges(4,4) * t193 - Icges(4,5) * t190;
t336 = Icges(4,1) * t194 + t191 * t224;
t335 = 2 * m(3);
t334 = 2 * m(4);
t333 = 2 * m(5);
t332 = m(4) / 0.2e1;
t331 = m(5) / 0.2e1;
t330 = t190 / 0.2e1;
t329 = t191 / 0.2e1;
t328 = -t194 / 0.2e1;
t327 = t194 / 0.2e1;
t326 = rSges(5,3) + pkin(6);
t163 = rSges(3,1) * t190 + rSges(3,2) * t193;
t325 = m(3) * t163;
t324 = pkin(2) * t193;
t323 = pkin(3) * t194;
t184 = t191 * pkin(3);
t322 = qJD(1) / 0.2e1;
t221 = Icges(5,5) * t189 + Icges(5,6) * t192;
t112 = Icges(5,3) * t190 - t193 * t221;
t274 = qJD(4) * t193;
t277 = qJD(2) * t193;
t279 = qJD(2) * t190;
t87 = (-Icges(5,5) * t192 + Icges(5,6) * t189) * t274 + (Icges(5,3) * t193 + t190 * t221) * qJD(2);
t247 = t189 * t115 * t274 + t112 * t277 + t190 * t87 + t279 * t346;
t93 = (-Icges(5,1) * t192 + t308) * t274 + (Icges(5,5) * t193 + t190 * t227) * qJD(2);
t316 = t189 * t93;
t53 = t112 * t190 - t193 * t346;
t90 = (Icges(5,2) * t189 - t307) * t274 + (Icges(5,6) * t193 + t190 * t223) * qJD(2);
t321 = ((-t316 + (-qJD(4) * t118 - t90) * t192) * t193 + t247) * t190 + t53 * t277;
t253 = qJD(1) * t190 + qJD(4);
t258 = t193 * t276;
t195 = -t191 * t253 + t258;
t254 = qJD(4) * t190 + qJD(1);
t211 = t254 * t189;
t70 = t192 * t195 - t194 * t211;
t212 = t192 * t254;
t71 = t189 * t195 + t194 * t212;
t320 = t71 * rSges(5,1) + t70 * rSges(5,2);
t319 = rSges(4,1) * t194;
t318 = rSges(4,2) * t190;
t317 = rSges(3,3) * t194;
t315 = -rSges(4,3) - qJ(3);
t292 = t192 * t194;
t295 = t191 * t189;
t138 = t190 * t292 - t295;
t294 = t191 * t192;
t139 = t189 * t296 + t294;
t83 = t139 * rSges(5,1) + t138 * rSges(5,2) + rSges(5,3) * t291;
t314 = pkin(6) * t291 + t184 + t83;
t293 = t191 * t193;
t140 = t189 * t194 + t190 * t294;
t141 = t190 * t295 - t292;
t245 = -t141 * rSges(5,1) - t140 * rSges(5,2);
t84 = rSges(5,3) * t293 - t245;
t313 = pkin(6) * t293 - t323 + t84;
t300 = qJ(3) * t190;
t299 = qJ(3) * t193;
t241 = t300 + t324;
t142 = t241 * t191;
t143 = pkin(2) * t291 + qJ(3) * t296;
t290 = t191 * t142 + t194 * t143;
t135 = qJD(2) * t241 - qJD(3) * t193;
t243 = -rSges(4,2) * t193 + rSges(4,3) * t190;
t289 = -t243 * qJD(2) - t135;
t161 = pkin(2) * t190 - t299;
t242 = rSges(4,3) * t193 + t318;
t288 = -t161 + t242;
t275 = qJD(3) * t190;
t287 = qJ(3) * t258 + t194 * t275;
t264 = t190 * t281;
t280 = qJD(1) * t194;
t286 = rSges(3,2) * t264 + rSges(3,3) * t280;
t285 = t194 * pkin(1) + t191 * pkin(5);
t284 = t191 ^ 2 + t194 ^ 2;
t114 = Icges(3,3) * t191 + t194 * t222;
t283 = qJD(1) * t114;
t125 = Icges(4,1) * t191 - t194 * t224;
t282 = qJD(1) * t125;
t278 = qJD(2) * t191;
t273 = -pkin(2) - t326;
t16 = -t112 * t347 + t70 * t115 + t71 * t118 + t138 * t90 + t139 * t93 + t291 * t87;
t79 = Icges(5,4) * t139 + Icges(5,2) * t138 + Icges(5,6) * t291;
t81 = Icges(5,1) * t139 + Icges(5,4) * t138 + Icges(5,5) * t291;
t240 = t189 * t81 + t192 * t79;
t35 = Icges(5,5) * t71 + Icges(5,6) * t70 - Icges(5,3) * t347;
t37 = Icges(5,4) * t71 + Icges(5,2) * t70 - Icges(5,6) * t347;
t39 = Icges(5,1) * t71 + Icges(5,4) * t70 - Icges(5,5) * t347;
t77 = Icges(5,5) * t139 + Icges(5,6) * t138 + Icges(5,3) * t291;
t9 = (qJD(2) * t240 + t35) * t190 + (qJD(2) * t77 - t189 * t39 - t192 * t37 + (t189 * t79 - t192 * t81) * qJD(4)) * t193;
t272 = t9 / 0.2e1 + t16 / 0.2e1;
t261 = t190 * t278;
t169 = pkin(2) * t261;
t259 = t191 * t277;
t262 = t193 * t280;
t271 = t142 * t280 + t191 * (pkin(2) * t262 + t191 * t275 - t169 + (t190 * t280 + t259) * qJ(3)) + t194 * (-pkin(2) * t347 - qJ(3) * t264 + t287);
t80 = Icges(5,4) * t141 + Icges(5,2) * t140 + Icges(5,6) * t293;
t82 = Icges(5,1) * t141 + Icges(5,4) * t140 + Icges(5,5) * t293;
t239 = t189 * t82 + t192 * t80;
t200 = -t261 + t262;
t210 = t253 * t194;
t68 = t192 * t210 + (t192 * t277 - t211) * t191;
t69 = t191 * t212 + (t210 + t259) * t189;
t34 = Icges(5,5) * t69 + Icges(5,6) * t68 + Icges(5,3) * t200;
t36 = Icges(5,4) * t69 + Icges(5,2) * t68 + Icges(5,6) * t200;
t38 = Icges(5,1) * t69 + Icges(5,4) * t68 + Icges(5,5) * t200;
t78 = Icges(5,5) * t141 + Icges(5,6) * t140 + Icges(5,3) * t293;
t10 = (qJD(2) * t239 + t34) * t190 + (qJD(2) * t78 - t189 * t38 - t192 * t36 + (t189 * t80 - t192 * t82) * qJD(4)) * t193;
t15 = t112 * t200 + t68 * t115 + t69 * t118 + t140 * t90 + t141 * t93 + t293 * t87;
t268 = t10 / 0.2e1 + t15 / 0.2e1;
t30 = t190 * t77 - t193 * t240;
t42 = t112 * t291 + t138 * t115 + t139 * t118;
t267 = -t30 / 0.2e1 - t42 / 0.2e1;
t31 = t190 * t78 - t193 * t239;
t43 = t112 * t293 + t115 * t140 + t118 * t141;
t266 = t31 / 0.2e1 + t43 / 0.2e1;
t179 = pkin(5) * t280;
t265 = t179 + t287;
t257 = -t224 * qJD(2) / 0.2e1 + t222 * t345;
t256 = -t279 / 0.2e1;
t244 = rSges(5,1) * t189 + rSges(5,2) * t192;
t129 = rSges(5,3) * t190 - t193 * t244;
t255 = pkin(6) * t190 + t129;
t110 = t288 * t194;
t252 = rSges(4,1) * t280 + rSges(4,2) * t347 + rSges(4,3) * t258;
t251 = t285 + t143;
t250 = -t161 - t255;
t249 = t69 * rSges(5,1) + t68 * rSges(5,2);
t248 = t190 * t315 - pkin(1);
t246 = rSges(3,1) * t193 - rSges(3,2) * t190;
t25 = t138 * t79 + t139 * t81 + t291 * t77;
t26 = t138 * t80 + t139 * t82 + t291 * t78;
t234 = t191 * t26 + t194 * t25;
t17 = t25 * t191 - t194 * t26;
t27 = t140 * t79 + t141 * t81 + t293 * t77;
t28 = t140 * t80 + t141 * t82 + t293 * t78;
t233 = t191 * t28 + t194 * t27;
t18 = t27 * t191 - t194 * t28;
t232 = t191 * t31 + t194 * t30;
t231 = t30 * t191 - t194 * t31;
t230 = t191 * t83 - t194 * t84;
t130 = rSges(3,1) * t291 + t339;
t131 = rSges(4,3) * t296 + t340;
t102 = (-rSges(5,1) * t192 + rSges(5,2) * t189) * t274 + (rSges(5,3) * t193 + t190 * t244) * qJD(2);
t209 = -pkin(6) * t277 - t102 - t135;
t208 = -pkin(1) - t246;
t76 = t250 * t194;
t207 = qJD(2) * t163;
t204 = qJD(2) * (Icges(4,4) * t190 + Icges(4,5) * t193);
t203 = qJD(2) * (-Icges(3,5) * t190 - Icges(3,6) * t193);
t198 = t193 * t273 - pkin(1) - t300;
t197 = (rSges(4,2) - pkin(2)) * t193 + t248;
t196 = t198 * t191;
t185 = t194 * pkin(5);
t180 = pkin(3) * t280;
t152 = t246 * qJD(2);
t144 = t161 * t281;
t132 = t191 * t243 - t319;
t128 = t191 * t246 - t317;
t109 = t288 * t191;
t105 = t130 + t285;
t104 = t191 * t208 + t185 + t317;
t101 = qJD(1) * t336 + t194 * t204;
t100 = t191 * t204 + t282;
t89 = t191 * t203 + t283;
t88 = -qJD(1) * t113 + t194 * t203;
t86 = t131 + t251;
t85 = t191 * t197 + t185 + t319;
t75 = t250 * t191;
t65 = t163 * t278 + ((-rSges(3,3) - pkin(5)) * t191 + t208 * t194) * qJD(1);
t64 = -rSges(3,1) * t347 - rSges(3,2) * t258 - pkin(1) * t281 + t179 + t286;
t63 = qJD(1) * t110 + t191 * t289;
t62 = t194 * t289 - t242 * t281 + t144;
t61 = -t129 * t291 + t190 * t83;
t60 = t129 * t293 - t190 * t84;
t59 = t191 * t213 + t194 * t336;
t58 = -t125 * t194 + t344;
t57 = -t191 * t336 + t342;
t56 = t191 * t125 + t214 * t194;
t55 = t191 * t114 - t215 * t194;
t54 = t191 * t113 - t341;
t52 = -t114 * t194 - t343;
t51 = -t113 * t194 - t191 * t216;
t50 = t251 + t314;
t49 = t185 + t196 + t245 + t323;
t47 = t131 * t194 + t191 * t132 + t290;
t46 = t169 + (-t275 + (t193 * t315 - t318) * qJD(2)) * t191 + ((-rSges(4,1) - pkin(5)) * t191 + t197 * t194) * qJD(1);
t45 = -pkin(2) * t260 + (t248 - t324) * t281 + t252 + t265;
t44 = t230 * t193;
t41 = -rSges(5,3) * t347 + t320;
t40 = rSges(5,3) * t200 + t249;
t33 = qJD(1) * t76 + t191 * t209;
t32 = t194 * t209 + t255 * t281 + t144;
t29 = t191 * t313 + t194 * t314 + t290;
t24 = t169 + (-t275 + (t190 * t326 - t299) * qJD(2)) * t191 + ((-pkin(3) - pkin(5)) * t191 + t198 * t194) * qJD(1) - t249;
t23 = qJD(1) * t196 + t260 * t273 + t180 + t265 + t320;
t22 = (-t129 * t278 - t40) * t190 + (-qJD(2) * t84 + t191 * t102 + t129 * t280) * t193;
t21 = (t129 * t276 + t41) * t190 + (qJD(2) * t83 - t102 * t194 + t129 * t281) * t193;
t20 = (qJD(1) * t132 + t252) * t194 + (t242 * t278 + (-t131 - t143 + t340) * qJD(1)) * t191 + t271;
t14 = t230 * t279 + (-t191 * t41 + t194 * t40 + (-t191 * t84 - t194 * t83) * qJD(1)) * t193;
t13 = t43 * t190 + t193 * t233;
t12 = t42 * t190 + t193 * t234;
t11 = (-pkin(6) * t260 + qJD(1) * t313 + t180 + t41) * t194 + (-pkin(6) * t261 + t40 + (-t143 - t314 + t184) * qJD(1)) * t191 + t271;
t8 = -t78 * t260 + t138 * t36 + t139 * t38 + t70 * t80 + t71 * t82 + (t194 * t34 - t281 * t78) * t193;
t7 = -t77 * t260 + t138 * t37 + t139 * t39 + t70 * t79 + t71 * t81 + (t194 * t35 - t281 * t77) * t193;
t6 = t78 * t262 + t140 * t36 + t141 * t38 + t68 * t80 + t69 * t82 + (t193 * t34 - t279 * t78) * t191;
t5 = t77 * t262 + t140 * t37 + t141 * t39 + t68 * t79 + t69 * t81 + (t193 * t35 - t279 * t77) * t191;
t4 = qJD(1) * t234 + t7 * t191 - t194 * t8;
t3 = qJD(1) * t233 + t5 * t191 - t194 * t6;
t2 = (-qJD(2) * t234 + t16) * t190 + (-qJD(1) * t17 + qJD(2) * t42 + t191 * t8 + t194 * t7) * t193;
t1 = (-qJD(2) * t233 + t15) * t190 + (-qJD(1) * t18 + qJD(2) * t43 + t191 * t6 + t194 * t5) * t193;
t19 = [t247 + (t104 * t65 + t105 * t64) * t335 + (t45 * t86 + t46 * t85) * t334 + (t23 * t50 + t24 * t49) * t333 - t193 * t316 + (-t118 * t274 - t193 * t90) * t192 + (t220 + t229 + t355) * t279 + (t218 + t226 - t354) * t277; (t194 * t257 - t268) * t194 + (t191 * t257 + t272) * t191 + m(5) * (t23 * t75 + t24 * t76 + t32 * t49 + t33 * t50) + m(4) * (t109 * t45 + t110 * t46 + t62 * t85 + t63 * t86) + m(3) * ((-t191 * t64 - t194 * t65) * t163 + (-t104 * t194 - t105 * t191) * t152) + ((-t350 * qJD(2) + t352 * t194) * t329 + (t351 * qJD(2) + t352 * t191) * t328 + (t348 * t328 - t349 * t329) * qJD(1)) * t190 + ((t348 * qJD(2) + t353 * t194) * t329 + (t349 * qJD(2) + t353 * t191) * t328 + (t350 * t328 + t351 * t329) * qJD(1)) * t193 + ((-t105 * t325 + (t117 / 0.2e1 - t121 / 0.2e1) * t193 + (t120 / 0.2e1 - t123 / 0.2e1) * t190 - t267) * t194 + (t104 * t325 + (t116 / 0.2e1 + t338 / 0.2e1) * t193 + (t119 / 0.2e1 + t337 / 0.2e1) * t190 + t266) * t191) * qJD(1); (t29 * t11 + t32 * t76 + t33 * t75) * t333 + t191 * t4 - t194 * t3 + (t109 * t63 + t110 * t62 + t47 * t20) * t334 + ((t191 * t128 + t130 * t194) * ((qJD(1) * t128 - t194 * t207 + t286) * t194 + (-t191 * t207 + (-t130 + t339) * qJD(1)) * t191) + t284 * t163 * t152) * t335 - t194 * ((t194 * t89 + (t52 + t341) * qJD(1)) * t194 + (t51 * qJD(1) + (-t117 * t277 - t120 * t279 + t283) * t191 + (-t88 + (t116 * t193 + t119 * t190) * qJD(2) - t215 * qJD(1)) * t194) * t191) + t191 * ((t191 * t101 + (t57 - t344) * qJD(1)) * t191 + (t56 * qJD(1) + (t277 * t338 + t279 * t337) * t194 + (-t100 + (t121 * t193 + t123 * t190) * qJD(2) + (t125 + t213) * qJD(1)) * t191) * t194) - t194 * ((t194 * t100 + (t58 - t342) * qJD(1)) * t194 + (t59 * qJD(1) + (t121 * t277 + t123 * t279 + t282) * t191 + (-t101 + (t190 * t337 + t193 * t338) * qJD(2) + t214 * qJD(1)) * t194) * t191) + t191 * ((t191 * t88 + (t54 + t343) * qJD(1)) * t191 + (t55 * qJD(1) + (t116 * t277 + t119 * t279) * t194 + (-t89 + (-t117 * t193 - t120 * t190) * qJD(2) + (t114 - t216) * qJD(1)) * t191) * t194) + (t18 + (-t51 - t59) * t194 + (t52 + t58) * t191) * t281 + (t17 + (-t54 - t57) * t194 + (t55 + t56) * t191) * t280; 0.2e1 * ((t191 * t86 + t194 * t85) * t332 + (t191 * t50 + t194 * t49) * t331) * t277 + 0.2e1 * ((t191 * t45 + t194 * t46 + t280 * t86 - t281 * t85) * t332 + (t191 * t23 + t194 * t24 + t280 * t50 - t281 * t49) * t331) * t190; 0.2e1 * ((t276 * t76 + t278 * t75 - t11) * t331 + (t109 * t278 + t110 * t276 - t20) * t332) * t193 + 0.2e1 * ((qJD(2) * t29 + t191 * t33 + t194 * t32 + t280 * t75 - t281 * t76) * t331 + (qJD(2) * t47 + t109 * t280 - t110 * t281 + t191 * t63 + t194 * t62) * t332) * t190; 0.4e1 * (t332 + t331) * (-0.1e1 + t284) * t190 * t277; m(5) * (t21 * t50 + t22 * t49 + t23 * t61 + t24 * t60) + (-t191 * t266 + t194 * t267) * t279 + (t272 * t194 + t268 * t191 + (t191 * t267 + t194 * t266) * qJD(1)) * t193 + t321; m(5) * (-t44 * t11 + t14 * t29 + t21 * t75 + t22 * t76 + t32 * t60 + t33 * t61) + (t12 * t322 + t17 * t256 + (qJD(1) * t30 - t10) * t330 - t1 / 0.2e1) * t194 + (t2 / 0.2e1 + t18 * t256 + (qJD(1) * t31 + t9) * t330 + t13 * t322) * t191 + (t3 * t329 + t4 * t327 + t231 * t345 + (t18 * t327 - t191 * t17 / 0.2e1) * qJD(1)) * t193; m(5) * ((-t14 + (t191 * t61 + t194 * t60) * qJD(2)) * t193 + (-qJD(2) * t44 + t191 * t21 + t194 * t22 + (-t191 * t60 + t194 * t61) * qJD(1)) * t190); (-t44 * t14 + t21 * t61 + t22 * t60) * t333 + ((-t194 * t12 - t191 * t13 - t190 * t232) * qJD(2) + t321) * t190 + (t194 * t2 + t191 * t1 + t190 * (t10 * t191 + t194 * t9) + (t53 * t190 + t193 * t232) * qJD(2) + (-t191 * t12 + t194 * t13 - t190 * t231) * qJD(1)) * t193;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t19(1), t19(2), t19(4), t19(7); t19(2), t19(3), t19(5), t19(8); t19(4), t19(5), t19(6), t19(9); t19(7), t19(8), t19(9), t19(10);];
Mq = res;
