% Calculate time derivative of joint inertia matrix for
% S5RPRPR5
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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:55:53
% EndTime: 2019-12-05 17:56:25
% DurationCPUTime: 12.41s
% Computational Cost: add. (13742->695), mult. (17234->974), div. (0->0), fcn. (16113->10), ass. (0->318)
t259 = qJ(3) + pkin(9);
t252 = cos(t259);
t261 = cos(pkin(8));
t266 = cos(qJ(1));
t251 = sin(t259);
t264 = sin(qJ(1));
t351 = t264 * t251;
t214 = t252 * t266 + t261 * t351;
t350 = t264 * t252;
t276 = -t251 * t266 + t261 * t350;
t260 = sin(pkin(8));
t358 = t260 * t264;
t134 = -Icges(5,5) * t276 + Icges(5,6) * t214 - Icges(5,3) * t358;
t265 = cos(qJ(3));
t347 = t265 * t266;
t263 = sin(qJ(3));
t349 = t264 * t263;
t225 = t261 * t349 + t347;
t348 = t264 * t265;
t354 = t263 * t266;
t274 = t261 * t348 - t354;
t150 = -Icges(4,5) * t274 + Icges(4,6) * t225 - Icges(4,3) * t358;
t398 = t134 + t150;
t356 = t261 * t266;
t217 = t252 * t356 + t351;
t275 = t251 * t356 - t350;
t357 = t260 * t266;
t135 = Icges(5,5) * t217 - Icges(5,6) * t275 + Icges(5,3) * t357;
t228 = t261 * t347 + t349;
t273 = t261 * t354 - t348;
t151 = Icges(4,5) * t228 - Icges(4,6) * t273 + Icges(4,3) * t357;
t387 = t135 + t151;
t324 = qJD(3) * t260;
t368 = Icges(5,4) * t251;
t199 = (-Icges(5,2) * t252 - t368) * t324;
t370 = Icges(4,4) * t263;
t219 = (-Icges(4,2) * t265 - t370) * t324;
t400 = -t251 * t199 - t263 * t219;
t168 = qJD(1) * t273 + qJD(3) * t274;
t169 = -qJD(1) * t228 + qJD(3) * t225;
t325 = qJD(1) * t266;
t308 = t260 * t325;
t102 = Icges(4,5) * t169 + Icges(4,6) * t168 - Icges(4,3) * t308;
t148 = qJD(1) * t275 + qJD(3) * t276;
t149 = -qJD(1) * t217 + qJD(3) * t214;
t88 = Icges(5,5) * t149 + Icges(5,6) * t148 - Icges(5,3) * t308;
t399 = -t88 - t102;
t198 = (-Icges(5,5) * t251 - Icges(5,6) * t252) * t324;
t367 = Icges(5,4) * t252;
t200 = (-Icges(5,1) * t251 - t367) * t324;
t218 = (-Icges(4,5) * t263 - Icges(4,6) * t265) * t324;
t369 = Icges(4,4) * t265;
t220 = (-Icges(4,1) * t263 - t369) * t324;
t397 = (-t198 - t218) * t261 + (t200 * t252 + t220 * t265) * t260;
t256 = t265 * pkin(3);
t249 = t256 + pkin(2);
t376 = pkin(2) - t249;
t279 = pkin(6) * t260 + t261 * t376;
t396 = -pkin(3) * t349 + t266 * t279;
t136 = -Icges(5,4) * t276 + Icges(5,2) * t214 - Icges(5,6) * t358;
t138 = -Icges(5,1) * t276 + Icges(5,4) * t214 - Icges(5,5) * t358;
t152 = -Icges(4,4) * t274 + Icges(4,2) * t225 - Icges(4,6) * t358;
t154 = -Icges(4,1) * t274 + Icges(4,4) * t225 - Icges(4,5) * t358;
t395 = t136 * t214 - t138 * t276 + t152 * t225 - t154 * t274 - t398 * t358;
t137 = Icges(5,4) * t217 - Icges(5,2) * t275 + Icges(5,6) * t357;
t139 = Icges(5,1) * t217 - Icges(5,4) * t275 + Icges(5,5) * t357;
t153 = Icges(4,4) * t228 - Icges(4,2) * t273 + Icges(4,6) * t357;
t155 = Icges(4,1) * t228 - Icges(4,4) * t273 + Icges(4,5) * t357;
t394 = t137 * t214 - t139 * t276 + t153 * t225 - t155 * t274 - t387 * t358;
t393 = t137 * t275 - t217 * t139 + t153 * t273 - t228 * t155 - t387 * t357;
t392 = t136 * t275 - t217 * t138 + t152 * t273 - t228 * t154 - t398 * t357;
t209 = -Icges(4,6) * t261 + (-Icges(4,2) * t263 + t369) * t260;
t210 = -Icges(4,5) * t261 + (Icges(4,1) * t265 - t370) * t260;
t192 = -Icges(5,6) * t261 + (-Icges(5,2) * t251 + t367) * t260;
t193 = -Icges(5,5) * t261 + (Icges(5,1) * t252 - t368) * t260;
t386 = -t192 * t252 - t193 * t251;
t391 = t397 + ((-t209 * t265 - t210 * t263 + t386) * qJD(3) + t400) * t260;
t262 = -qJ(4) - pkin(6);
t257 = -pkin(7) + t262;
t390 = t260 * (rSges(6,3) - t257);
t379 = pkin(3) * t263;
t235 = pkin(4) * t251 + t379;
t300 = -t235 + t379;
t388 = t300 * t264;
t378 = pkin(4) * t252;
t234 = t249 + t378;
t330 = t234 - t249;
t299 = t330 * t261;
t166 = qJD(1) * t225 - qJD(3) * t228;
t167 = -qJD(1) * t274 - qJD(3) * t273;
t326 = qJD(1) * t264;
t309 = t260 * t326;
t101 = Icges(4,5) * t167 + Icges(4,6) * t166 - Icges(4,3) * t309;
t146 = qJD(1) * t214 - qJD(3) * t217;
t147 = -qJD(1) * t276 - qJD(3) * t275;
t87 = Icges(5,5) * t147 + Icges(5,6) * t146 - Icges(5,3) * t309;
t385 = (t101 + t87) * t266;
t384 = 2 * m(4);
t383 = 2 * m(5);
t382 = 2 * m(6);
t381 = m(5) / 0.2e1;
t380 = m(6) / 0.2e1;
t253 = qJ(5) + t259;
t248 = cos(t253);
t258 = qJD(3) + qJD(5);
t359 = t258 * t260;
t247 = sin(t253);
t366 = Icges(6,4) * t247;
t180 = (-Icges(6,2) * t248 - t366) * t359;
t187 = -Icges(6,5) * t261 + (Icges(6,1) * t248 - t366) * t260;
t179 = (-Icges(6,5) * t247 - Icges(6,6) * t248) * t359;
t365 = Icges(6,4) * t248;
t181 = (-Icges(6,1) * t247 - t365) * t359;
t297 = t260 * t248 * t181 - t261 * t179;
t186 = -Icges(6,6) * t261 + (-Icges(6,2) * t247 + t365) * t260;
t314 = t258 * t248 * t186;
t48 = (-t314 + (-t187 * t258 - t180) * t247) * t260 + t297;
t375 = t48 * t261;
t374 = -rSges(3,3) - qJ(2);
t323 = qJD(3) * t263;
t316 = pkin(3) * t323;
t293 = t261 * t316;
t322 = qJD(3) * t265;
t315 = pkin(3) * t322;
t311 = t262 * t308 + t264 * t293 + t266 * t315;
t321 = qJD(4) * t260;
t271 = -t264 * t321 + t311;
t111 = qJD(1) * t396 + t271;
t341 = t149 * rSges(5,1) + t148 * rSges(5,2);
t372 = rSges(5,3) * t308 - t111 - t341;
t182 = (-rSges(6,1) * t247 - rSges(6,2) * t248) * t359;
t353 = t264 * t247;
t202 = t248 * t266 + t261 * t353;
t205 = t248 * t356 + t353;
t130 = qJD(1) * t202 - t205 * t258;
t352 = t264 * t248;
t277 = t247 * t356 - t352;
t278 = -t247 * t266 + t261 * t352;
t131 = -qJD(1) * t278 - t258 * t277;
t285 = -t131 * rSges(6,1) - t130 * rSges(6,2);
t81 = -rSges(6,3) * t309 - t285;
t371 = t182 * t357 + t261 * t81;
t231 = t235 * qJD(3);
t362 = t231 * t261;
t361 = t234 * t261;
t255 = t266 * qJ(2);
t346 = -qJ(2) - t235;
t320 = qJD(4) * t266;
t307 = t261 * t326;
t328 = pkin(2) * t307 + pkin(6) * t309;
t331 = t249 * t307 + t266 * t293;
t110 = (t262 * t326 + t320) * t260 + (t263 * t325 + t264 * t322) * pkin(3) + t328 - t331;
t306 = t260 * t323;
t230 = -pkin(3) * t306 - qJD(4) * t261;
t345 = t261 * t110 + t230 * t357;
t284 = -t205 * rSges(6,1) + rSges(6,2) * t277;
t127 = rSges(6,3) * t357 - t284;
t188 = -rSges(6,3) * t261 + (rSges(6,1) * t248 - rSges(6,2) * t247) * t260;
t98 = t261 * t127 + t188 * t357;
t132 = qJD(1) * t277 + t258 * t278;
t133 = -qJD(1) * t205 + t202 * t258;
t344 = t133 * rSges(6,1) + t132 * rSges(6,2);
t335 = -rSges(5,1) * t276 + t214 * rSges(5,2);
t140 = -rSges(5,3) * t358 + t335;
t329 = pkin(3) * t354 + t262 * t358;
t156 = t264 * t279 + t329;
t343 = -t140 - t156;
t286 = -t217 * rSges(5,1) + rSges(5,2) * t275;
t141 = rSges(5,3) * t357 - t286;
t242 = t262 * t357;
t157 = -t242 - t396;
t342 = -t141 - t157;
t190 = (pkin(6) + t262) * t261 - t376 * t260;
t340 = t261 * t157 + t190 * t357;
t339 = t169 * rSges(4,1) + t168 * rSges(4,2);
t338 = t182 * t358 + t188 * t308;
t337 = t190 * t308 + t230 * t358;
t336 = -rSges(6,1) * t278 + t202 * rSges(6,2);
t232 = (t256 + t378) * qJD(3);
t334 = t266 * t232 + t257 * t308;
t333 = -rSges(4,1) * t274 + t225 * rSges(4,2);
t332 = t266 * t235 + t257 * t358;
t327 = t257 - t262;
t82 = -rSges(6,3) * t308 + t344;
t318 = -t111 - t264 * t362 - (-t266 * t299 + t388) * qJD(1) + t311 - t334 - t82;
t317 = rSges(4,3) * t358;
t114 = -t264 * t299 - t329 + t332;
t126 = -rSges(6,3) * t358 + t336;
t313 = -t114 - t126 - t156;
t115 = t242 - t388 + (-t257 * t260 + t299) * t266;
t312 = -t115 - t127 - t157;
t122 = -Icges(6,4) * t278 + Icges(6,2) * t202 - Icges(6,6) * t358;
t124 = -Icges(6,1) * t278 + Icges(6,4) * t202 - Icges(6,5) * t358;
t76 = Icges(6,5) * t133 + Icges(6,6) * t132 - Icges(6,3) * t308;
t78 = Icges(6,4) * t133 + Icges(6,2) * t132 - Icges(6,6) * t308;
t80 = Icges(6,1) * t133 + Icges(6,4) * t132 - Icges(6,5) * t308;
t12 = -t261 * t76 + ((-t122 * t258 + t80) * t248 + (-t124 * t258 - t78) * t247) * t260;
t120 = -Icges(6,5) * t278 + Icges(6,6) * t202 - Icges(6,3) * t358;
t121 = Icges(6,5) * t205 - Icges(6,6) * t277 + Icges(6,3) * t357;
t123 = Icges(6,4) * t205 - Icges(6,2) * t277 + Icges(6,6) * t357;
t125 = Icges(6,1) * t205 - Icges(6,4) * t277 + Icges(6,5) * t357;
t75 = Icges(6,5) * t131 + Icges(6,6) * t130 - Icges(6,3) * t309;
t77 = Icges(6,4) * t131 + Icges(6,2) * t130 - Icges(6,6) * t309;
t79 = Icges(6,1) * t131 + Icges(6,4) * t130 - Icges(6,5) * t309;
t13 = -t261 * t75 + ((-t123 * t258 + t79) * t248 + (-t125 * t258 - t77) * t247) * t260;
t185 = -Icges(6,3) * t261 + (Icges(6,5) * t248 - Icges(6,6) * t247) * t260;
t22 = t130 * t186 + t131 * t187 - t277 * t180 + t205 * t181 + (t179 * t266 - t185 * t326) * t260;
t36 = t120 * t357 - t122 * t277 + t205 * t124;
t37 = t121 * t357 - t123 * t277 + t205 * t125;
t46 = -t120 * t261 + (-t122 * t247 + t124 * t248) * t260;
t47 = -t121 * t261 + (-t123 * t247 + t125 * t248) * t260;
t310 = -t261 * (-t375 + (-t12 * t264 + t13 * t266 + (-t264 * t47 - t266 * t46) * qJD(1)) * t260) + (-t22 * t261 + (-(t130 * t122 + t131 * t124 + t205 * t80 - t277 * t78) * t264 - t36 * t325 + (t130 * t123 + t131 * t125 + t205 * t79 - t277 * t77) * t266 - t37 * t326 + (-(-t120 * t326 + t266 * t76) * t264 + (-t121 * t326 + t266 * t75) * t266) * t260) * t260) * t357;
t302 = -pkin(1) - t361;
t301 = -qJ(2) - t379;
t298 = t346 * t264;
t250 = pkin(1) * t326;
t294 = -qJD(2) * t264 + t250;
t290 = rSges(3,1) * t261 - rSges(3,2) * t260;
t289 = -t167 * rSges(4,1) - t166 * rSges(4,2);
t288 = -t228 * rSges(4,1) + rSges(4,2) * t273;
t287 = -t147 * rSges(5,1) - t146 * rSges(5,2);
t283 = -t321 + t362;
t282 = -pkin(1) - t290;
t281 = -rSges(5,3) * t260 - t249 * t261 - pkin(1);
t280 = -rSges(6,3) * t260 + t302;
t272 = -pkin(2) * t261 - pkin(1) + (-rSges(4,3) - pkin(6)) * t260;
t23 = t132 * t186 + t133 * t187 + t202 * t180 - t278 * t181 + (-t179 * t264 - t185 * t325) * t260;
t68 = -t185 * t358 + t186 * t202 - t187 * t278;
t69 = t185 * t357 - t186 * t277 + t205 * t187;
t270 = -(t12 + t23) * t358 / 0.2e1 + (t13 + t22) * t357 / 0.2e1 - ((t47 + t69) * t264 + (t46 + t68) * t266) * qJD(1) * t260 / 0.2e1;
t269 = -t264 * qJ(2) + t266 * t272;
t184 = t264 * t374 + t266 * t282;
t268 = t264 * t301 + t266 * t281;
t34 = -t120 * t358 + t122 * t202 - t124 * t278;
t35 = -t121 * t358 + t123 * t202 - t125 * t278;
t2 = -t23 * t261 + (-(t132 * t122 + t133 * t124 + t202 * t78 - t278 * t80) * t264 - t34 * t325 + (t132 * t123 + t133 * t125 + t202 * t77 - t278 * t79) * t266 - t35 * t326 + (-(-t120 * t325 - t264 * t76) * t264 + (-t121 * t325 - t264 * t75) * t266) * t260) * t260;
t6 = -t68 * t261 + (-t264 * t34 + t266 * t35) * t260;
t7 = -t69 * t261 + (-t264 * t36 + t266 * t37) * t260;
t267 = (-t264 * t2 + (-t264 * t7 - t266 * t6) * qJD(1)) * t260 + t310;
t254 = qJD(2) * t266;
t221 = (-rSges(4,1) * t263 - rSges(4,2) * t265) * t324;
t213 = -rSges(4,3) * t261 + (rSges(4,1) * t265 - rSges(4,2) * t263) * t260;
t208 = -Icges(4,3) * t261 + (Icges(4,5) * t265 - Icges(4,6) * t263) * t260;
t201 = (-rSges(5,1) * t251 - rSges(5,2) * t252) * t324;
t195 = -rSges(5,3) * t261 + (rSges(5,1) * t252 - rSges(5,2) * t251) * t260;
t191 = -Icges(5,3) * t261 + (Icges(5,5) * t252 - Icges(5,6) * t251) * t260;
t189 = (-t231 + t316) * t260;
t183 = rSges(3,3) * t266 + t264 * t282 + t255;
t176 = t190 * t358;
t174 = t188 * t358;
t162 = t260 * t330 + t261 * t327;
t161 = qJD(1) * t184 + t254;
t160 = (t264 * t290 + t266 * t374) * qJD(1) + t294;
t159 = rSges(4,3) * t357 - t288;
t158 = -t317 + t333;
t142 = t156 * t309;
t118 = t126 * t309;
t117 = t269 + t288;
t116 = t264 * t272 + t255 + t333;
t113 = t261 * t159 + t213 * t357;
t112 = -t158 * t261 + t213 * t358;
t108 = -rSges(4,3) * t308 + t339;
t107 = -rSges(4,3) * t309 - t289;
t106 = Icges(4,1) * t169 + Icges(4,4) * t168 - Icges(4,5) * t308;
t105 = Icges(4,1) * t167 + Icges(4,4) * t166 - Icges(4,5) * t309;
t104 = Icges(4,4) * t169 + Icges(4,2) * t168 - Icges(4,6) * t308;
t103 = Icges(4,4) * t167 + Icges(4,2) * t166 - Icges(4,6) * t309;
t100 = t242 + t268 + t286;
t99 = t264 * t281 + t255 + t329 + t335;
t97 = -t126 * t261 + t174;
t96 = t208 * t357 - t209 * t273 + t228 * t210;
t95 = -t208 * t358 + t209 * t225 - t210 * t274;
t93 = -rSges(5,3) * t309 - t287;
t92 = Icges(5,1) * t149 + Icges(5,4) * t148 - Icges(5,5) * t308;
t91 = Icges(5,1) * t147 + Icges(5,4) * t146 - Icges(5,5) * t309;
t90 = Icges(5,4) * t149 + Icges(5,2) * t148 - Icges(5,6) * t308;
t89 = Icges(5,4) * t147 + Icges(5,2) * t146 - Icges(5,6) * t309;
t86 = t298 + (t302 - t390) * t266 + t284;
t85 = t264 * t280 + t255 + t332 + t336;
t84 = qJD(1) * t269 + t254 + t339;
t83 = (t317 - t255) * qJD(1) + t289 + t294 + t328;
t74 = t191 * t357 - t192 * t275 + t217 * t193;
t73 = -t191 * t358 + t192 * t214 - t193 * t276;
t71 = (-t126 * t266 - t127 * t264) * t260;
t66 = -t231 * t356 + (t232 - t315) * t264 + (-t300 * t266 + (t260 * t327 - t361) * t264) * qJD(1) + t331;
t65 = -t261 * t108 + (t213 * t325 + t221 * t264) * t260;
t64 = t261 * t107 + (-t213 * t326 + t221 * t266) * t260;
t61 = -t151 * t261 + (-t153 * t263 + t155 * t265) * t260;
t60 = -t150 * t261 + (-t152 * t263 + t154 * t265) * t260;
t58 = t261 * t141 + t195 * t357 + t340;
t57 = t195 * t358 + t261 * t343 + t176;
t52 = -t135 * t261 + (-t137 * t251 + t139 * t252) * t260;
t51 = -t134 * t261 + (-t136 * t251 + t138 * t252) * t260;
t50 = qJD(1) * t268 + t254 + t271 + t341;
t49 = -t260 * t320 + t250 + (-qJD(2) - t315) * t264 + (t301 * t266 + (rSges(5,3) - t262) * t358) * qJD(1) + t287 + t331;
t45 = -t261 * t82 + t338;
t44 = -t188 * t309 + t371;
t33 = t254 + t283 * t264 + (t266 * t280 + t298) * qJD(1) + t334 + t344;
t32 = t250 + t283 * t266 + (-qJD(2) - t232) * t264 + (t346 * t266 + (t361 + t390) * t264) * qJD(1) + t285;
t31 = t168 * t209 + t169 * t210 + t225 * t219 - t274 * t220 + (-t208 * t325 - t218 * t264) * t260;
t30 = t166 * t209 + t167 * t210 - t273 * t219 + t228 * t220 + (-t208 * t326 + t218 * t266) * t260;
t29 = t148 * t192 + t149 * t193 + t214 * t199 - t276 * t200 + (-t191 * t325 - t198 * t264) * t260;
t28 = t146 * t192 + t147 * t193 - t275 * t199 + t217 * t200 + (-t191 * t326 + t198 * t266) * t260;
t27 = t261 * t115 + t162 * t357 + t340 + t98;
t26 = t162 * t358 + t261 * t313 + t174 + t176;
t25 = t372 * t261 + (t195 * t325 + t201 * t264) * t260 + t337;
t24 = t261 * t93 + (t201 * t266 + (-t190 - t195) * t326) * t260 + t345;
t19 = -t101 * t261 + (-t103 * t263 + t105 * t265 + (-t153 * t265 - t155 * t263) * qJD(3)) * t260;
t18 = -t102 * t261 + (-t104 * t263 + t106 * t265 + (-t152 * t265 - t154 * t263) * qJD(3)) * t260;
t17 = (t264 * t312 + t266 * t313) * t260;
t16 = t118 + (-t264 * t81 + (-qJD(1) * t127 - t82) * t266) * t260;
t15 = -t261 * t87 + (-t251 * t89 + t252 * t91 + (-t137 * t252 - t139 * t251) * qJD(3)) * t260;
t14 = -t261 * t88 + (-t251 * t90 + t252 * t92 + (-t136 * t252 - t138 * t251) * qJD(3)) * t260;
t9 = (t162 * t325 + t189 * t264) * t260 + t318 * t261 + t337 + t338;
t8 = t261 * t66 + (t189 * t266 + (-t162 - t188 - t190) * t326) * t260 + t345 + t371;
t5 = t142 + ((qJD(1) * t140 - t110 - t93) * t264 + (qJD(1) * t342 + t372) * t266) * t260;
t4 = t118 + t142 + ((qJD(1) * t114 - t110 - t66 - t81) * t264 + (qJD(1) * t312 + t318) * t266) * t260;
t1 = [0.2e1 * m(3) * (t160 * t184 + t161 * t183) - t210 * t306 + (t116 * t84 + t117 * t83) * t384 + (t100 * t49 + t50 * t99) * t383 - t247 * t187 * t359 + (t32 * t86 + t33 * t85) * t382 + t297 + t386 * t324 + (-t247 * t180 - t209 * t322 - t314 + t400) * t260 + t397; m(3) * (t160 * t266 + t264 * t161 + (t183 * t266 - t184 * t264) * qJD(1)) + m(4) * (t264 * t84 + t266 * t83 + (t116 * t266 - t117 * t264) * qJD(1)) + m(5) * (t264 * t50 + t266 * t49 + (-t100 * t264 + t266 * t99) * qJD(1)) + m(6) * (t264 * t33 + t266 * t32 + (-t264 * t86 + t266 * t85) * qJD(1)); 0; (-t48 - t391) * t261 + m(4) * (t112 * t84 + t113 * t83 + t116 * t65 + t117 * t64) + m(5) * (t100 * t24 + t25 * t99 + t49 * t58 + t50 * t57) + m(6) * (t26 * t33 + t27 * t32 + t8 * t86 + t85 * t9) + ((t19 / 0.2e1 + t15 / 0.2e1 + t28 / 0.2e1 + t30 / 0.2e1) * t266 + (-t18 / 0.2e1 - t14 / 0.2e1 - t29 / 0.2e1 - t31 / 0.2e1) * t264 + ((-t60 / 0.2e1 - t51 / 0.2e1 - t73 / 0.2e1 - t95 / 0.2e1) * t266 + (-t61 / 0.2e1 - t52 / 0.2e1 - t74 / 0.2e1 - t96 / 0.2e1) * t264) * qJD(1)) * t260 + t270; m(4) * (t65 * t264 + t266 * t64 + (t112 * t266 - t113 * t264) * qJD(1)) + m(5) * (t24 * t266 + t25 * t264 + (-t264 * t58 + t266 * t57) * qJD(1)) + m(6) * (t9 * t264 + t266 * t8 + (t26 * t266 - t264 * t27) * qJD(1)); -t7 * t309 - t6 * t308 - t2 * t358 + (t17 * t4 + t26 * t9 + t27 * t8) * t382 + (t58 * t24 + t57 * t25) * t383 + (t112 * t65 + t113 * t64) * t384 + t310 + ((t264 * t342 + t266 * t343) * t5 * t383 + (-t158 * t266 - t159 * t264) * (-t107 * t264 - t108 * t266 + (t158 * t264 - t159 * t266) * qJD(1)) * t384 * t260 + (-t264 * t392 + t266 * t393) * t309 + (t264 * t395 - t394 * t266) * t308 + (t394 * t326 + t395 * t325 + (-t225 * t103 + t105 * t274 - t148 * t137 - t149 * t139 - t168 * t153 - t169 * t155 - t214 * t89 + t276 * t91 + t308 * t387) * t266 + (t225 * t104 - t274 * t106 + t148 * t136 + t149 * t138 + t168 * t152 + t169 * t154 + t214 * t90 - t276 * t92 + (t264 * t399 - t398 * t325 + t385) * t260) * t264) * t358 + (t393 * t326 + t392 * t325 + (-t273 * t103 + t228 * t105 + t146 * t137 + t147 * t139 + t166 * t153 + t167 * t155 - t275 * t89 + t217 * t91 + (-t326 * t387 + t385) * t260) * t266 + (t273 * t104 - t228 * t106 - t166 * t152 - t167 * t154 - t146 * t136 - t147 * t138 + t275 * t90 - t217 * t92 + (t266 * t399 + t398 * t326) * t260) * t264) * t357) * t260 + (t391 * t261 + (t29 + t31) * t358 + (-t30 - t28) * t357 + (t96 + t74) * t309 + (t73 + t95) * t308 + ((-t15 - t19) * t266 + (t14 + t18) * t264 + ((t51 + t60) * t266 + (t52 + t61) * t264) * qJD(1)) * t260) * t261; 0.2e1 * ((-t100 * t325 - t264 * t49 + t266 * t50 - t326 * t99) * t381 + (-t264 * t32 + t266 * t33 - t325 * t86 - t326 * t85) * t380) * t260; 0; 0.2e1 * (-m(6) * t4 / 0.2e1 - m(5) * t5 / 0.2e1) * t261 + 0.2e1 * ((-t26 * t326 - t264 * t8 + t266 * t9 - t27 * t325) * t380 + (-t24 * t264 + t25 * t266 - t325 * t58 - t326 * t57) * t381) * t260; 0; -t375 + m(6) * (t32 * t98 + t33 * t97 + t44 * t86 + t45 * t85) + t270; m(6) * (t45 * t264 + t266 * t44 + (-t264 * t98 + t266 * t97) * qJD(1)); m(6) * (t16 * t17 + t26 * t45 + t27 * t44 + t4 * t71 + t8 * t98 + t9 * t97) + t267; m(6) * (-t16 * t261 + (-t264 * t44 + t266 * t45 + (-t264 * t97 - t266 * t98) * qJD(1)) * t260); (t16 * t71 + t44 * t98 + t45 * t97) * t382 + t267;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
