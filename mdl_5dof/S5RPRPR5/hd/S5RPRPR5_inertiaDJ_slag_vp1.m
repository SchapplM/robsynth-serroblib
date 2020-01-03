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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:41:27
% EndTime: 2020-01-03 11:42:08
% DurationCPUTime: 12.45s
% Computational Cost: add. (13742->693), mult. (17234->970), div. (0->0), fcn. (16113->10), ass. (0->304)
t271 = qJ(3) + pkin(9);
t262 = cos(t271);
t273 = cos(pkin(8));
t278 = cos(qJ(1));
t261 = sin(t271);
t276 = sin(qJ(1));
t347 = t276 * t261;
t215 = -t262 * t278 - t273 * t347;
t346 = t276 * t262;
t216 = -t261 * t278 + t273 * t346;
t272 = sin(pkin(8));
t357 = t272 * t276;
t140 = Icges(5,5) * t216 + Icges(5,6) * t215 + Icges(5,3) * t357;
t277 = cos(qJ(3));
t343 = t277 * t278;
t275 = sin(qJ(3));
t345 = t276 * t275;
t228 = -t273 * t345 - t343;
t344 = t276 * t277;
t352 = t275 * t278;
t229 = t273 * t344 - t352;
t158 = Icges(4,5) * t229 + Icges(4,6) * t228 + Icges(4,3) * t357;
t394 = t140 + t158;
t354 = t273 * t278;
t217 = t261 * t354 - t346;
t286 = t262 * t354 + t347;
t356 = t272 * t278;
t141 = -Icges(5,5) * t286 + Icges(5,6) * t217 - Icges(5,3) * t356;
t230 = t273 * t352 - t344;
t285 = t273 * t343 + t345;
t159 = -Icges(4,5) * t285 + Icges(4,6) * t230 - Icges(4,3) * t356;
t387 = t141 + t159;
t327 = qJD(3) * t272;
t367 = Icges(5,4) * t261;
t201 = (-Icges(5,2) * t262 - t367) * t327;
t369 = Icges(4,4) * t275;
t221 = (-Icges(4,2) * t277 - t369) * t327;
t396 = -t261 * t201 - t275 * t221;
t176 = -qJD(1) * t230 - qJD(3) * t229;
t177 = qJD(1) * t285 + qJD(3) * t228;
t328 = qJD(1) * t278;
t311 = t272 * t328;
t104 = Icges(4,5) * t177 + Icges(4,6) * t176 + Icges(4,3) * t311;
t156 = -qJD(1) * t217 - qJD(3) * t216;
t157 = qJD(1) * t286 + qJD(3) * t215;
t90 = Icges(5,5) * t157 + Icges(5,6) * t156 + Icges(5,3) * t311;
t395 = t90 + t104;
t200 = (-Icges(5,5) * t261 - Icges(5,6) * t262) * t327;
t366 = Icges(5,4) * t262;
t202 = (-Icges(5,1) * t261 - t366) * t327;
t220 = (-Icges(4,5) * t275 - Icges(4,6) * t277) * t327;
t368 = Icges(4,4) * t277;
t222 = (-Icges(4,1) * t275 - t368) * t327;
t393 = (-t200 - t220) * t273 + (t202 * t262 + t222 * t277) * t272;
t142 = Icges(5,4) * t216 + Icges(5,2) * t215 + Icges(5,6) * t357;
t144 = Icges(5,1) * t216 + Icges(5,4) * t215 + Icges(5,5) * t357;
t160 = Icges(4,4) * t229 + Icges(4,2) * t228 + Icges(4,6) * t357;
t162 = Icges(4,1) * t229 + Icges(4,4) * t228 + Icges(4,5) * t357;
t392 = t142 * t215 + t144 * t216 + t160 * t228 + t162 * t229 + t394 * t357;
t143 = -Icges(5,4) * t286 + Icges(5,2) * t217 - Icges(5,6) * t356;
t145 = -Icges(5,1) * t286 + Icges(5,4) * t217 - Icges(5,5) * t356;
t161 = -Icges(4,4) * t285 + Icges(4,2) * t230 - Icges(4,6) * t356;
t163 = -Icges(4,1) * t285 + Icges(4,4) * t230 - Icges(4,5) * t356;
t391 = t143 * t215 + t145 * t216 + t161 * t228 + t163 * t229 + t387 * t357;
t390 = -t217 * t143 + t145 * t286 - t230 * t161 + t163 * t285 + t387 * t356;
t389 = -t217 * t142 + t144 * t286 - t230 * t160 + t162 * t285 + t394 * t356;
t211 = -Icges(4,6) * t273 + (-Icges(4,2) * t275 + t368) * t272;
t212 = -Icges(4,5) * t273 + (Icges(4,1) * t277 - t369) * t272;
t194 = -Icges(5,6) * t273 + (-Icges(5,2) * t261 + t366) * t272;
t195 = -Icges(5,5) * t273 + (Icges(5,1) * t262 - t367) * t272;
t386 = -t194 * t262 - t195 * t261;
t388 = t393 + ((-t211 * t277 - t212 * t275 + t386) * qJD(3) + t396) * t272;
t274 = -qJ(4) - pkin(6);
t269 = -pkin(7) + t274;
t331 = t269 - t274;
t304 = t331 * t272;
t174 = qJD(1) * t228 + qJD(3) * t285;
t175 = qJD(1) * t229 + qJD(3) * t230;
t329 = qJD(1) * t276;
t312 = t272 * t329;
t103 = Icges(4,5) * t175 + Icges(4,6) * t174 + Icges(4,3) * t312;
t154 = qJD(1) * t215 + qJD(3) * t286;
t155 = qJD(1) * t216 + qJD(3) * t217;
t89 = Icges(5,5) * t155 + Icges(5,6) * t154 + Icges(5,3) * t312;
t385 = (-t103 - t89) * t278;
t267 = t277 * pkin(3);
t259 = t267 + pkin(2);
t376 = pkin(4) * t262;
t236 = t259 + t376;
t361 = t236 * t273;
t384 = t272 * (rSges(6,3) - t269) + t361;
t383 = 2 * m(4);
t382 = 2 * m(5);
t381 = 2 * m(6);
t380 = m(5) / 0.2e1;
t379 = m(6) / 0.2e1;
t378 = pkin(2) * t273;
t377 = pkin(3) * t275;
t375 = pkin(6) * t272;
t374 = -pkin(2) + t259;
t373 = pkin(6) + t274;
t263 = qJ(5) + t271;
t256 = cos(t263);
t270 = qJD(3) + qJD(5);
t358 = t270 * t272;
t255 = sin(t263);
t365 = Icges(6,4) * t255;
t182 = (-Icges(6,2) * t256 - t365) * t358;
t189 = -Icges(6,5) * t273 + (Icges(6,1) * t256 - t365) * t272;
t181 = (-Icges(6,5) * t255 - Icges(6,6) * t256) * t358;
t364 = Icges(6,4) * t256;
t183 = (-Icges(6,1) * t255 - t364) * t358;
t303 = t272 * t256 * t183 - t273 * t181;
t188 = -Icges(6,6) * t273 + (-Icges(6,2) * t255 + t364) * t272;
t319 = t270 * t256 * t188;
t48 = (-t319 + (-t189 * t270 - t182) * t255) * t272 + t303;
t372 = t48 * t273;
t190 = -rSges(6,3) * t273 + (rSges(6,1) * t256 - rSges(6,2) * t255) * t272;
t349 = t276 * t255;
t204 = -t256 * t278 - t273 * t349;
t287 = t256 * t354 + t349;
t136 = qJD(1) * t204 + t270 * t287;
t348 = t276 * t256;
t205 = -t255 * t278 + t273 * t348;
t206 = t255 * t354 - t348;
t137 = qJD(1) * t205 + t206 * t270;
t290 = -t137 * rSges(6,1) - t136 * rSges(6,2);
t83 = rSges(6,3) * t312 - t290;
t370 = t190 * t312 + t273 * t83;
t239 = pkin(4) * t261 + t377;
t360 = t239 * t278;
t355 = t273 * t276;
t234 = (t267 + t376) * qJD(3);
t351 = t276 * t234;
t350 = t276 * t239;
t324 = qJD(4) * t272;
t248 = t278 * t324;
t326 = qJD(3) * t275;
t321 = pkin(3) * t326;
t300 = t273 * t321;
t325 = qJD(3) * t277;
t320 = pkin(3) * t325;
t322 = qJD(1) * t377;
t313 = t274 * t312 + t276 * t320 + t278 * t322;
t279 = t278 * t300 - t248 - t313;
t114 = (t273 * t374 - t375) * t329 + t279;
t192 = t272 * t374 + t273 * t373;
t342 = t273 * t114 + t192 * t312;
t132 = t205 * rSges(6,1) + t204 * rSges(6,2) + rSges(6,3) * t357;
t289 = rSges(6,1) * t287 - t206 * rSges(6,2);
t133 = -rSges(6,3) * t356 - t289;
t73 = t132 * t356 + t133 * t357;
t146 = t216 * rSges(5,1) + t215 * rSges(5,2) + rSges(5,3) * t357;
t237 = t259 * t355;
t288 = -t272 * t373 - t378;
t164 = -pkin(3) * t352 + t276 * t288 + t237;
t341 = -t146 - t164;
t335 = pkin(2) * t354 + pkin(6) * t356;
t336 = pkin(3) * t345 + t259 * t354;
t165 = t274 * t356 + t335 - t336;
t340 = t164 * t356 + t165 * t357;
t197 = -rSges(5,3) * t273 + (rSges(5,1) * t262 - rSges(5,2) * t261) * t272;
t339 = -t192 - t197;
t309 = t272 * t326;
t232 = -pkin(3) * t309 - qJD(4) * t273;
t338 = -(-rSges(5,1) * t261 - rSges(5,2) * t262) * t327 - t232;
t337 = t236 - t259;
t334 = pkin(1) * t328 + qJ(2) * t329;
t333 = qJ(2) * t328 + qJD(2) * t276;
t332 = t278 * pkin(1) + t276 * qJ(2);
t330 = qJD(1) * t272;
t138 = -qJD(1) * t206 - t205 * t270;
t139 = qJD(1) * t287 + t204 * t270;
t84 = t139 * rSges(6,1) + t138 * rSges(6,2) + rSges(6,3) * t311;
t323 = t133 * t311 + t84 * t356 + t83 * t357;
t247 = t276 * t324;
t310 = t273 * t328;
t283 = -t259 * t310 + (t300 - t322) * t276;
t281 = t247 - t283;
t115 = (qJD(1) * t288 - t320) * t278 + t281;
t318 = t114 * t357 + t115 * t356 + t165 * t311;
t223 = t236 * t355;
t118 = t223 - t237 + (-t239 + t377) * t278 - t276 * t304;
t317 = -t118 - t132 - t164;
t96 = t157 * rSges(5,1) + t156 * rSges(5,2) + rSges(5,3) * t311;
t170 = t272 * t337 + t273 * t331;
t316 = -t170 - t190 - t192;
t110 = t177 * rSges(4,1) + t176 * rSges(4,2) + rSges(4,3) * t311;
t184 = (-rSges(6,1) * t255 - rSges(6,2) * t256) * t358;
t233 = t239 * qJD(3);
t297 = -t233 + t321;
t315 = -t297 * t272 - t184 - t232;
t314 = -t233 * t355 + t236 * t310 + t239 * t329;
t166 = t229 * rSges(4,1) + t228 * rSges(4,2) + rSges(4,3) * t357;
t305 = t278 * t339;
t296 = t316 * t278;
t295 = t375 + t378;
t294 = rSges(3,1) * t273 - rSges(3,2) * t272;
t293 = -rSges(4,1) * t175 - rSges(4,2) * t174;
t292 = -t155 * rSges(5,1) - t154 * rSges(5,2);
t291 = rSges(5,1) * t286 - t217 * rSges(5,2);
t128 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t357;
t130 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t357;
t78 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t311;
t80 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t311;
t82 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t311;
t12 = -t273 * t78 + ((-t128 * t270 + t82) * t256 + (-t130 * t270 - t80) * t255) * t272;
t129 = -Icges(6,4) * t287 + Icges(6,2) * t206 - Icges(6,6) * t356;
t131 = -Icges(6,1) * t287 + Icges(6,4) * t206 - Icges(6,5) * t356;
t77 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t312;
t79 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t312;
t81 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t312;
t13 = -t273 * t77 + ((-t129 * t270 + t81) * t256 + (-t131 * t270 - t79) * t255) * t272;
t187 = -Icges(6,3) * t273 + (Icges(6,5) * t256 - Icges(6,6) * t255) * t272;
t22 = t136 * t188 + t137 * t189 + t206 * t182 - t287 * t183 + (-t181 * t278 + t187 * t329) * t272;
t23 = t138 * t188 + t139 * t189 + t204 * t182 + t205 * t183 + (t181 * t276 + t187 * t328) * t272;
t126 = Icges(6,5) * t205 + Icges(6,6) * t204 + Icges(6,3) * t357;
t46 = -t126 * t273 + (-t128 * t255 + t130 * t256) * t272;
t127 = -Icges(6,5) * t287 + Icges(6,6) * t206 - Icges(6,3) * t356;
t47 = -t127 * t273 + (-t129 * t255 + t131 * t256) * t272;
t68 = t187 * t357 + t188 * t204 + t189 * t205;
t69 = -t187 * t356 + t206 * t188 - t189 * t287;
t284 = (t12 + t23) * t357 / 0.2e1 - (t13 + t22) * t356 / 0.2e1 + ((t47 + t69) * t276 + (t46 + t68) * t278) * t330 / 0.2e1;
t34 = t126 * t357 + t128 * t204 + t130 * t205;
t35 = t127 * t357 + t129 * t204 + t131 * t205;
t36 = -t126 * t356 + t206 * t128 - t130 * t287;
t37 = -t127 * t356 + t206 * t129 - t131 * t287;
t282 = -(-t22 * t273 + ((t136 * t128 + t137 * t130 + t206 * t80 - t287 * t82) * t276 + t36 * t328 - (t136 * t129 + t137 * t131 + t206 * t79 - t287 * t81) * t278 + t37 * t329 + ((t126 * t329 - t278 * t78) * t276 - (t127 * t329 - t278 * t77) * t278) * t272) * t272) * t356 - t273 * (-t372 + (t12 * t276 - t13 * t278 + (t276 * t47 + t278 * t46) * qJD(1)) * t272) + (-t23 * t273 + ((t138 * t128 + t139 * t130 + t204 * t80 + t205 * t82) * t276 + t34 * t328 - (t138 * t129 + t139 * t131 + t204 * t79 + t205 * t81) * t278 + t35 * t329 + ((t126 * t328 + t276 * t78) * t276 - (t127 * t328 + t276 * t77) * t278) * t272) * t272) * t357 + (-t69 * t273 + (t276 * t36 - t278 * t37) * t272) * t312 + (-t68 * t273 + (t276 * t34 - t278 * t35) * t272) * t311;
t167 = -rSges(4,1) * t285 + t230 * rSges(4,2) - rSges(4,3) * t356;
t280 = t276 * rSges(3,3) + t278 * t294;
t266 = t276 * pkin(1);
t224 = (-rSges(4,1) * t275 - rSges(4,2) * t277) * t327;
t214 = -rSges(4,3) * t273 + (rSges(4,1) * t277 - rSges(4,2) * t275) * t272;
t210 = -Icges(4,3) * t273 + (Icges(4,5) * t277 - Icges(4,6) * t275) * t272;
t193 = -Icges(5,3) * t273 + (Icges(5,5) * t262 - Icges(5,6) * t261) * t272;
t186 = t280 + t332;
t185 = t266 + (-rSges(3,3) - qJ(2)) * t278 + t294 * t276;
t169 = qJD(1) * t280 - qJD(2) * t278 + t334;
t168 = (rSges(3,3) * t278 + (-pkin(1) - t294) * t276) * qJD(1) + t333;
t153 = t273 * t165;
t147 = -rSges(5,3) * t356 - t291;
t125 = t273 * t133;
t121 = -t167 + t332 + t335;
t120 = -qJ(2) * t278 + t276 * t295 + t166 + t266;
t119 = -t350 + (-t361 + t304) * t278 + t336;
t117 = t273 * t167 - t214 * t356;
t116 = -t166 * t273 - t214 * t357;
t109 = rSges(4,3) * t312 - t293;
t108 = Icges(4,1) * t177 + Icges(4,4) * t176 + Icges(4,5) * t311;
t107 = Icges(4,1) * t175 + Icges(4,4) * t174 + Icges(4,5) * t312;
t106 = Icges(4,4) * t177 + Icges(4,2) * t176 + Icges(4,6) * t311;
t105 = Icges(4,4) * t175 + Icges(4,2) * t174 + Icges(4,6) * t312;
t102 = (rSges(5,3) - t274) * t356 + t291 + t332 + t336;
t101 = -t274 * t357 + t237 + t266 + (-qJ(2) - t377) * t278 + t146;
t100 = -t190 * t356 + t125;
t99 = -t132 * t273 - t190 * t357;
t98 = -t210 * t356 + t230 * t211 - t212 * t285;
t97 = t210 * t357 + t211 * t228 + t212 * t229;
t95 = rSges(5,3) * t312 - t292;
t94 = Icges(5,1) * t157 + Icges(5,4) * t156 + Icges(5,5) * t311;
t93 = Icges(5,1) * t155 + Icges(5,4) * t154 + Icges(5,5) * t312;
t92 = Icges(5,4) * t157 + Icges(5,2) * t156 + Icges(5,6) * t311;
t91 = Icges(5,4) * t155 + Icges(5,2) * t154 + Icges(5,6) * t312;
t88 = t384 * t278 + t289 + t332 + t350;
t87 = -t269 * t357 + t223 + t266 + (-qJ(2) - t239) * t278 + t132;
t86 = (qJD(1) * t295 - qJD(2)) * t278 + t110 + t334;
t85 = (-t378 - pkin(1) + (-rSges(4,3) - pkin(6)) * t272) * t329 + t293 + t333;
t76 = -t193 * t356 + t217 * t194 - t195 * t286;
t75 = t193 * t357 + t194 * t215 + t195 * t216;
t67 = (-qJD(1) * t304 - t234 + t320) * t278 + t283 + t314;
t66 = -t351 - t297 * t354 + (-t360 + (-t269 * t272 + t273 * t337) * t276) * qJD(1) + t313;
t65 = -t273 * t110 + (-t214 * t328 - t224 * t276) * t272;
t64 = t273 * t109 + (t214 * t329 - t224 * t278) * t272;
t61 = -t159 * t273 + (-t161 * t275 + t163 * t277) * t272;
t60 = -t158 * t273 + (-t160 * t275 + t162 * t277) * t272;
t58 = t273 * t147 + t272 * t305 + t153;
t57 = t273 * t341 + t339 * t357;
t52 = -t141 * t273 + (-t143 * t261 + t145 * t262) * t272;
t51 = -t140 * t273 + (-t142 * t261 + t144 * t262) * t272;
t50 = (-t274 * t330 - qJD(2) - t320) * t278 + t281 + t96 + t334;
t49 = (-rSges(5,3) * t272 - t259 * t273 - pkin(1)) * t329 - t279 + t292 + t333;
t45 = -t273 * t84 + (-t184 * t276 - t190 * t328) * t272;
t44 = -t184 * t356 + t370;
t33 = t247 + (-t269 * t330 - qJD(2) - t234) * t278 + t314 + t84 + t334;
t32 = -t233 * t354 + t351 + t248 + (t360 + (-pkin(1) - t384) * t276) * qJD(1) + t290 + t333;
t31 = t176 * t211 + t177 * t212 + t228 * t221 + t229 * t222 + (t210 * t328 + t220 * t276) * t272;
t30 = t174 * t211 + t175 * t212 + t230 * t221 - t285 * t222 + (t210 * t329 - t220 * t278) * t272;
t29 = t156 * t194 + t157 * t195 + t215 * t201 + t216 * t202 + (t193 * t328 + t200 * t276) * t272;
t28 = t154 * t194 + t155 * t195 + t217 * t201 - t286 * t202 + (t193 * t329 - t200 * t278) * t272;
t27 = t273 * t119 + t272 * t296 + t125 + t153;
t26 = t273 * t317 + t316 * t357;
t25 = (-t115 - t96) * t273 + (qJD(1) * t305 + t276 * t338) * t272;
t24 = t273 * t95 + (t197 * t329 + t278 * t338) * t272 + t342;
t19 = -t103 * t273 + (-t105 * t275 + t107 * t277 + (-t161 * t277 - t163 * t275) * qJD(3)) * t272;
t18 = -t104 * t273 + (-t106 * t275 + t108 * t277 + (-t160 * t277 - t162 * t275) * qJD(3)) * t272;
t17 = (t118 * t278 + t119 * t276) * t272 + t340 + t73;
t16 = -t132 * t312 + t323;
t15 = -t273 * t89 + (-t261 * t91 + t262 * t93 + (-t143 * t262 - t145 * t261) * qJD(3)) * t272;
t14 = -t273 * t90 + (-t261 * t92 + t262 * t94 + (-t142 * t262 - t144 * t261) * qJD(3)) * t272;
t9 = (-t115 - t67 - t84) * t273 + (qJD(1) * t296 + t276 * t315) * t272;
t8 = t273 * t66 + (t170 * t329 + t278 * t315) * t272 + t342 + t370;
t7 = (t276 * t95 + t278 * t96 + (t147 * t278 + t276 * t341) * qJD(1)) * t272 + t318;
t4 = (t276 * t66 + t278 * t67 + (t119 * t278 + t276 * t317) * qJD(1)) * t272 + t318 + t323;
t1 = [0.2e1 * m(3) * (t168 * t186 + t169 * t185) - t212 * t309 + (t120 * t86 + t121 * t85) * t383 + (t101 * t50 + t102 * t49) * t382 - t255 * t189 * t358 + (t32 * t88 + t33 * t87) * t381 + t303 + t386 * t327 + (-t255 * t182 - t211 * t325 - t319 + t396) * t272 + t393; m(3) * (-t168 * t278 - t276 * t169 + (-t185 * t278 + t186 * t276) * qJD(1)) + m(4) * (-t276 * t86 - t278 * t85 + (-t120 * t278 + t121 * t276) * qJD(1)) + m(5) * (-t276 * t50 - t278 * t49 + (-t101 * t278 + t102 * t276) * qJD(1)) + m(6) * (-t276 * t33 - t278 * t32 + (t276 * t88 - t278 * t87) * qJD(1)); 0; (-t48 - t388) * t273 + m(4) * (t116 * t86 + t117 * t85 + t120 * t65 + t121 * t64) + m(5) * (t101 * t25 + t102 * t24 + t49 * t58 + t50 * t57) + m(6) * (t26 * t33 + t27 * t32 + t8 * t88 + t87 * t9) + ((-t19 / 0.2e1 - t15 / 0.2e1 - t28 / 0.2e1 - t30 / 0.2e1) * t278 + (t18 / 0.2e1 + t14 / 0.2e1 + t29 / 0.2e1 + t31 / 0.2e1) * t276 + ((t60 / 0.2e1 + t51 / 0.2e1 + t97 / 0.2e1 + t75 / 0.2e1) * t278 + (t61 / 0.2e1 + t52 / 0.2e1 + t98 / 0.2e1 + t76 / 0.2e1) * t276) * qJD(1)) * t272 + t284; m(4) * (-t65 * t276 - t278 * t64 + (-t116 * t278 + t117 * t276) * qJD(1)) + m(5) * (-t24 * t278 - t25 * t276 + (t276 * t58 - t278 * t57) * qJD(1)) + m(6) * (-t9 * t276 - t278 * t8 + (-t26 * t278 + t27 * t276) * qJD(1)); (t17 * t4 + t26 * t9 + t27 * t8) * t381 + (t58 * t24 + t57 * t25 + t340 * t7) * t382 + (t116 * t65 + t117 * t64) * t383 + t282 + ((t146 * t278 + t147 * t276) * t7 * t382 + (t166 * t278 + t167 * t276) * (t109 * t276 + t110 * t278 + (-t166 * t276 + t167 * t278) * qJD(1)) * t383 * t272 + (-t389 * t276 + t390 * t278) * t312 + (t392 * t276 - t391 * t278) * t311 + (t391 * t329 + t392 * t328 + (-t228 * t105 - t229 * t107 - t156 * t143 - t157 * t145 - t176 * t161 - t177 * t163 - t215 * t91 - t216 * t93 - t387 * t311) * t278 + (t228 * t106 + t229 * t108 + t156 * t142 + t157 * t144 + t176 * t160 + t177 * t162 + t215 * t92 + t216 * t94 + (t276 * t395 + t394 * t328 + t385) * t272) * t276) * t357 + (t390 * t329 + t389 * t328 + (t230 * t105 - t285 * t107 + t154 * t143 + t155 * t145 + t174 * t161 + t175 * t163 + t217 * t91 - t286 * t93 + (t387 * t329 + t385) * t272) * t278 + (-t230 * t106 + t285 * t108 - t174 * t160 - t175 * t162 - t154 * t142 - t155 * t144 - t217 * t92 + t286 * t94 + (t278 * t395 - t394 * t329) * t272) * t276) * t356) * t272 + (t388 * t273 + (-t29 - t31) * t357 + (t30 + t28) * t356 + (-t98 - t76) * t312 + (-t75 - t97) * t311 + ((t15 + t19) * t278 + (-t14 - t18) * t276 + ((-t51 - t60) * t278 + (-t52 - t61) * t276) * qJD(1)) * t272) * t273; 0.2e1 * ((t101 * t329 + t102 * t328 + t276 * t49 - t278 * t50) * t380 + (t276 * t32 - t278 * t33 + t328 * t88 + t329 * t87) * t379) * t272; 0; 0.2e1 * (-m(6) * t4 / 0.2e1 - m(5) * t7 / 0.2e1) * t273 + 0.2e1 * ((t26 * t329 + t27 * t328 + t276 * t8 - t278 * t9) * t379 + (t24 * t276 - t25 * t278 + t328 * t58 + t329 * t57) * t380) * t272; 0; -t372 + m(6) * (t100 * t32 + t33 * t99 + t44 * t88 + t45 * t87) + t284; m(6) * (-t45 * t276 - t278 * t44 + (t100 * t276 - t278 * t99) * qJD(1)); m(6) * (t100 * t8 + t16 * t17 + t26 * t45 + t27 * t44 + t4 * t73 + t9 * t99) + t282; m(6) * (-t16 * t273 + (t276 * t44 - t278 * t45 + (t100 * t278 + t276 * t99) * qJD(1)) * t272); (t100 * t44 + t16 * t73 + t45 * t99) * t381 + t282;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
