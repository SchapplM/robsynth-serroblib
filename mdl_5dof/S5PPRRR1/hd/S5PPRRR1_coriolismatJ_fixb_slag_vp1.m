% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:38
% EndTime: 2019-12-05 15:12:47
% DurationCPUTime: 5.08s
% Computational Cost: add. (35645->352), mult. (32695->541), div. (0->0), fcn. (35201->8), ass. (0->231)
t249 = sin(pkin(8));
t246 = t249 ^ 2;
t250 = cos(pkin(8));
t247 = t250 ^ 2;
t388 = t246 + t247;
t248 = pkin(9) + qJ(3);
t243 = qJ(4) + t248;
t238 = sin(t243);
t239 = cos(t243);
t252 = sin(qJ(5));
t253 = cos(qJ(5));
t284 = rSges(6,1) * t253 - rSges(6,2) * t252;
t192 = -t239 * rSges(6,3) + t238 * t284;
t176 = t192 * t249;
t177 = t192 * t250;
t230 = rSges(5,1) * t238 + rSges(5,2) * t239;
t156 = t388 * t230;
t241 = sin(t248);
t242 = cos(t248);
t273 = -Icges(4,5) * t241 - Icges(4,6) * t242;
t272 = -Icges(5,5) * t238 - Icges(5,6) * t239;
t362 = t249 / 0.2e1;
t361 = -t250 / 0.2e1;
t303 = qJD(3) + qJD(4);
t326 = t239 * (-Icges(6,5) * t252 - Icges(6,6) * t253) * t238;
t339 = Icges(6,4) * t253;
t274 = -Icges(6,2) * t252 + t339;
t188 = -Icges(6,6) * t239 + t238 * t274;
t313 = -t188 + (-Icges(6,1) * t252 - t339) * t238;
t340 = Icges(6,4) * t252;
t279 = Icges(6,1) * t253 - t340;
t190 = -Icges(6,5) * t239 + t238 * t279;
t312 = t190 + (-Icges(6,2) * t253 - t340) * t238;
t231 = rSges(5,1) * t239 - rSges(5,2) * t238;
t141 = t388 * t231;
t109 = (-t141 + t231) * t156;
t320 = t250 * t253;
t323 = t249 * t252;
t226 = -t239 * t323 - t320;
t321 = t250 * t252;
t322 = t249 * t253;
t227 = t239 * t322 - t321;
t331 = t238 * t249;
t150 = rSges(6,1) * t227 + rSges(6,2) * t226 + rSges(6,3) * t331;
t228 = -t239 * t321 + t322;
t229 = t239 * t320 + t323;
t330 = t238 * t250;
t151 = rSges(6,1) * t229 + rSges(6,2) * t228 + rSges(6,3) * t330;
t233 = pkin(4) * t239 + pkin(7) * t238;
t112 = t249 * t150 + t250 * t151 + t233 * t388;
t232 = pkin(4) * t238 - pkin(7) * t239;
t119 = -t249 * t176 - t250 * t177 - t232 * t388;
t311 = -t192 - t232;
t152 = t311 * t249;
t193 = t238 * rSges(6,3) + t239 * t284;
t310 = -t193 - t233;
t153 = t310 * t249;
t154 = t311 * t250;
t155 = t310 * t250;
t55 = t112 * t119 + t152 * t153 + t154 * t155;
t381 = m(5) * t109 + m(6) * t55;
t357 = pkin(3) * t241;
t286 = t311 - t357;
t137 = t286 * t249;
t139 = t286 * t250;
t356 = pkin(3) * t242;
t314 = t388 * t356;
t94 = t112 + t314;
t46 = t119 * t94 + t137 * t153 + t139 * t155;
t116 = t141 + t314;
t294 = -t230 - t357;
t182 = t294 * t249;
t184 = t294 * t250;
t88 = -t116 * t156 + (-t182 * t249 - t184 * t250) * t231;
t380 = -m(5) * t88 - m(6) * t46;
t160 = Icges(6,5) * t226 - Icges(6,6) * t227;
t218 = Icges(6,4) * t226;
t148 = Icges(6,1) * t227 + Icges(6,5) * t331 + t218;
t316 = -Icges(6,2) * t227 + t148 + t218;
t342 = Icges(6,4) * t227;
t146 = Icges(6,2) * t226 + Icges(6,6) * t331 + t342;
t318 = Icges(6,1) * t226 - t146 - t342;
t78 = t160 * t331 + t226 * t316 + t227 * t318;
t161 = Icges(6,5) * t228 - Icges(6,6) * t229;
t219 = Icges(6,4) * t228;
t149 = Icges(6,1) * t229 + Icges(6,5) * t330 + t219;
t315 = -Icges(6,2) * t229 + t149 + t219;
t341 = Icges(6,4) * t229;
t147 = Icges(6,2) * t228 + Icges(6,6) * t330 + t341;
t317 = Icges(6,1) * t228 - t147 - t341;
t79 = t161 * t331 + t226 * t315 + t227 * t317;
t44 = t249 * t79 - t250 * t78;
t80 = t160 * t330 + t228 * t316 + t229 * t318;
t81 = t161 * t330 + t228 * t315 + t229 * t317;
t45 = t249 * t81 - t250 * t80;
t355 = t361 * t44 + t362 * t45;
t271 = Icges(6,5) * t253 - Icges(6,6) * t252;
t186 = -Icges(6,3) * t239 + t238 * t271;
t378 = t388 * t273;
t377 = t388 * t272;
t374 = 2 * qJD(3);
t373 = m(5) / 0.2e1;
t372 = m(6) / 0.2e1;
t267 = t150 * t250 - t151 * t249;
t120 = t267 * t238;
t130 = t150 * t239 + t192 * t331;
t131 = -t151 * t239 - t192 * t330;
t297 = t119 * t120 + t130 * t155 + t131 * t153;
t110 = (t193 * t249 - t150) * t238;
t111 = (-t193 * t250 + t151) * t238;
t98 = t267 * t239 + (-t176 * t250 + t177 * t249) * t238;
t300 = t110 * t139 + t111 * t137 + t94 * t98;
t370 = m(6) * (t297 + t300);
t21 = t110 * t154 + t111 * t152 + t112 * t98 + t297;
t369 = m(6) * t21;
t217 = (-rSges(6,1) * t252 - rSges(6,2) * t253) * t238;
t166 = rSges(6,1) * t226 - rSges(6,2) * t227;
t167 = rSges(6,1) * t228 - rSges(6,2) * t229;
t133 = t166 * t249 + t167 * t250;
t73 = t94 * t133;
t93 = t112 * t133;
t368 = m(6) * (t73 + t93 + ((-t139 - t154) * t250 + (-t137 - t152) * t249) * t217);
t367 = m(6) * (t110 * t130 + t111 * t131 + t120 * t98);
t364 = m(6) * (t110 * t249 - t111 * t250);
t363 = -t239 / 0.2e1;
t360 = t250 / 0.2e1;
t358 = m(6) * (-t153 * t250 + t155 * t249);
t298 = t364 / 0.2e1;
t354 = t303 * t298;
t299 = -t364 / 0.2e1;
t86 = 0.2e1 * (t98 / 0.4e1 - t133 / 0.4e1) * m(6);
t319 = t86 * qJD(1);
t353 = qJD(2) * t299 - t319;
t352 = m(6) * qJD(5);
t106 = -m(5) * t156 + m(6) * t119;
t87 = (t133 + t98) * t372;
t348 = qJD(4) * t106 + qJD(5) * t87;
t347 = qJD(4) * t358 + qJD(5) * t298;
t145 = Icges(6,5) * t229 + Icges(6,6) * t228 + Icges(6,3) * t330;
t101 = t145 * t331 + t147 * t226 + t149 * t227;
t335 = t101 * t250;
t144 = Icges(6,5) * t227 + Icges(6,6) * t226 + Icges(6,3) * t331;
t102 = t144 * t330 + t146 * t228 + t148 * t229;
t334 = t102 * t249;
t329 = t239 * t144;
t328 = t239 * t145;
t327 = t239 * t186;
t325 = t239 * t249;
t324 = t239 * t250;
t100 = t144 * t331 + t146 * t226 + t148 * t227;
t117 = t186 * t331 + t188 * t226 + t190 * t227;
t189 = Icges(6,6) * t238 + t239 * t274;
t191 = Icges(6,5) * t238 + t239 * t279;
t266 = -t188 * t252 + t190 * t253;
t263 = (Icges(6,3) * t238 + t239 * t271 - t266) * t239;
t172 = t188 * t249;
t174 = t190 * t249;
t269 = -t146 * t252 + t148 * t253;
t265 = -t186 * t249 - t269;
t255 = t238 * t265 + t329;
t69 = -t226 * t172 - t227 * t174 + t249 * t255;
t173 = t188 * t250;
t175 = t190 * t250;
t268 = -t147 * t252 + t149 * t253;
t264 = -t186 * t250 - t268;
t254 = t238 * t264 + t328;
t70 = -t226 * t173 - t227 * t175 + t249 * t254;
t15 = (t335 - t226 * t189 - t227 * t191 + (t100 - t327) * t249) * t239 + (t70 * t250 + t117 + (t69 - t263) * t249) * t238;
t103 = t145 * t330 + t147 * t228 + t149 * t229;
t118 = t186 * t330 + t188 * t228 + t190 * t229;
t71 = -t228 * t172 - t229 * t174 + t250 * t255;
t72 = -t228 * t173 - t229 * t175 + t250 * t254;
t16 = (t334 - t228 * t189 - t229 * t191 + (t103 - t327) * t250) * t239 + (t71 * t249 + t118 + (t72 - t263) * t250) * t238;
t124 = t238 * t266 - t327;
t107 = t238 * t269 - t329;
t108 = t238 * t268 - t328;
t270 = t107 * t249 + t108 * t250;
t76 = -t265 * t239 + (t172 * t252 - t174 * t253 + t144) * t238;
t77 = -t264 * t239 + (t173 * t252 - t175 * t253 + t145) * t238;
t22 = (t263 + t270) * t239 + (t77 * t250 + t76 * t249 - (-t189 * t252 + t191 * t253 + t186) * t239 + t124) * t238;
t49 = -t117 * t239 + (t100 * t249 + t335) * t238;
t50 = -t118 * t239 + (t103 * t250 + t334) * t238;
t54 = -t124 * t239 + t238 * t270;
t3 = t367 + (t50 * t360 + t49 * t362 - t22 / 0.2e1) * t239 + (t16 * t360 + t15 * t362 + t54 / 0.2e1) * t238;
t302 = qJD(2) * t298 + qJD(5) * t3 + t319;
t301 = t368 / 0.2e1 + t355;
t296 = t331 / 0.2e1;
t295 = t330 / 0.2e1;
t293 = -t231 - t356;
t208 = t272 * t249;
t209 = t272 * t250;
t35 = t249 * t70 - t250 * t69;
t36 = t249 * t72 - t250 * t71;
t292 = (t36 + t246 * t209 + (-t249 * t208 + t377) * t250) * t362 + (t35 + t247 * t208 + (-t250 * t209 + t377) * t249) * t361;
t287 = t388 * t357;
t285 = t310 - t356;
t234 = rSges(4,1) * t241 + rSges(4,2) * t242;
t262 = t120 * t133 + (-t130 * t250 - t131 * t249) * t217;
t261 = t15 * t361 + t16 * t362 + t36 * t295 + t35 * t296 + (t249 * t77 - t250 * t76) * t363 + (-t100 * t250 + t101 * t249) * t325 / 0.2e1 + (-t102 * t250 + t103 * t249) * t324 / 0.2e1 + t238 * (-t107 * t250 + t108 * t249) / 0.2e1 - t355;
t29 = -(t226 * t312 + t227 * t313) * t239 + (t79 * t250 + (t78 - t326) * t249) * t238;
t30 = -(t228 * t312 + t229 * t313) * t239 + (t80 * t249 + (t81 - t326) * t250) * t238;
t90 = -t239 * t160 + (-t252 * t316 + t253 * t318) * t238;
t91 = -t239 * t161 + (-t252 * t315 + t253 * t317) * t238;
t256 = -t15 * t331 / 0.2e1 - t16 * t330 / 0.2e1 + t239 * t22 / 0.2e1 - t367 + t30 * t362 + t29 * t361 + t44 * t296 + t45 * t295 - t49 * t325 / 0.2e1 - t50 * t324 / 0.2e1 + (t249 * t91 - t250 * t90) * t363 - t238 * t54 / 0.2e1;
t221 = t273 * t250;
t220 = t273 * t249;
t185 = t293 * t250;
t183 = t293 * t249;
t157 = t388 * t234;
t140 = t285 * t250;
t138 = t285 * t249;
t136 = -t287 - t156;
t135 = -t167 * t239 - t217 * t330;
t134 = t166 * t239 + t217 * t331;
t129 = (t166 * t250 - t167 * t249) * t238;
t113 = -t287 + t119;
t84 = t86 * qJD(5);
t63 = qJD(5) * t299;
t59 = t93 + (-t152 * t249 - t154 * t250) * t217;
t53 = t73 + (-t137 * t249 - t139 * t250) * t217;
t20 = t369 / 0.2e1;
t17 = t370 / 0.2e1;
t10 = m(6) * t59 + t355;
t9 = m(6) * t53 + t355;
t8 = t292 + t381;
t7 = t8 * qJD(4);
t6 = t292 - t380;
t5 = t20 - t370 / 0.2e1 + t301;
t4 = t17 - t369 / 0.2e1 + t301;
t1 = t17 + t20 - t368 / 0.2e1 + t261;
t2 = [0, 0, (-m(4) * t157 / 0.2e1 + t136 * t373 + t113 * t372) * t374 + t348, qJD(3) * t106 + t348, t129 * t352 + t303 * t87; 0, 0, ((-t183 * t250 + t185 * t249) * t373 + (-t138 * t250 + t140 * t249) * t372) * t374 + t347, qJD(3) * t358 + t347, (t134 * t249 - t135 * t250) * t352 + t354; -t84, t63, (m(6) * (t113 * t94 + t137 * t138 + t139 * t140) + (t246 * t221 + (-t249 * t220 + t378) * t250) * t362 + (t247 * t220 + (-t250 * t221 + t378) * t249) * t361 + m(5) * (t116 * t136 + t182 * t183 + t184 * t185) + t292 + m(4) * (-t157 + t234) * t388 * (rSges(4,1) * t242 - rSges(4,2) * t241)) * qJD(3) + t6 * qJD(4) + t9 * qJD(5), t6 * qJD(3) + t4 * qJD(5) + (t292 + 0.2e1 * (t55 + t46) * t372 + 0.2e1 * (t88 + t109) * t373 - t381) * qJD(4), t9 * qJD(3) + t4 * qJD(4) + (t256 + m(6) * (t129 * t94 + t134 * t139 + t135 * t137 + t262)) * qJD(5) + t353; -t84, t63, t7 + t5 * qJD(5) + ((t112 * t113 + t138 * t152 + t140 * t154 + t46) * t372 + (t136 * t141 + (-t183 * t249 - t185 * t250) * t230 + t88) * t373) * t374 + (t292 + t380) * qJD(3), qJD(3) * t8 + qJD(5) * t10 + t7, t5 * qJD(3) + t10 * qJD(4) + (t256 + m(6) * (t112 * t129 + t134 * t154 + t135 * t152 + t262)) * qJD(5) + t353; t303 * t86, t354, ((t113 * t120 + t130 * t140 + t131 * t138 + t300 - t53) * m(6) + t261) * qJD(3) + t1 * qJD(4) + t302, t1 * qJD(3) + ((t21 - t59) * m(6) + t261) * qJD(4) + t302, t303 * t3 + (m(6) * (t120 * t129 + t130 * t134 + t131 * t135) - t239 ^ 2 * t326 / 0.2e1 + (t30 * t360 + t29 * t362 + (t91 * t250 + t90 * t249 - (-t252 * t312 + t253 * t313) * t239) * t363) * t238) * qJD(5);];
Cq = t2;
