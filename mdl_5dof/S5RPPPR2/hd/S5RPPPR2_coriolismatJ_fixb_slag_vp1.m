% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:05
% EndTime: 2022-01-23 08:59:14
% DurationCPUTime: 4.52s
% Computational Cost: add. (12042->374), mult. (30213->535), div. (0->0), fcn. (37708->10), ass. (0->241)
t247 = sin(pkin(9));
t249 = sin(pkin(7));
t251 = cos(pkin(8));
t250 = cos(pkin(9));
t252 = cos(pkin(7));
t289 = t252 * t250;
t224 = t249 * t247 + t251 * t289;
t254 = sin(qJ(1));
t248 = sin(pkin(8));
t256 = cos(qJ(1));
t283 = t256 * t248;
t198 = t224 * t254 - t250 * t283;
t282 = t256 * t251;
t286 = t254 * t248;
t227 = t252 * t286 + t282;
t253 = sin(qJ(5));
t255 = cos(qJ(5));
t160 = -t198 * t253 + t227 * t255;
t294 = t248 * t253;
t199 = t224 * t255 + t252 * t294;
t293 = t248 * t255;
t226 = t250 * t293 - t253 * t251;
t360 = -t199 * t254 + t256 * t226;
t127 = Icges(6,5) * t160 + Icges(6,6) * t360;
t285 = t254 * t251;
t228 = t252 * t285 - t283;
t291 = t249 * t254;
t203 = -t228 * t247 + t250 * t291;
t386 = t127 * t203;
t161 = t198 * t255 + t227 * t253;
t200 = -t224 * t253 + t252 * t293;
t225 = t250 * t294 + t255 * t251;
t165 = t200 * t254 + t256 * t225;
t363 = Icges(6,2) * t165 - Icges(6,6) * t203;
t100 = Icges(6,4) * t161 + t363;
t382 = Icges(6,4) * t360;
t131 = Icges(6,1) * t160 + t382;
t385 = -t100 + t131;
t230 = t252 * t282 + t286;
t290 = t249 * t256;
t204 = t230 * t247 - t250 * t290;
t352 = t200 * t256 - t254 * t225;
t279 = Icges(6,2) * t352 + Icges(6,6) * t204;
t163 = t199 * t256 + t254 * t226;
t383 = Icges(6,4) * t163;
t102 = t279 + t383;
t130 = Icges(6,1) * t352 - t383;
t384 = -t102 + t130;
t229 = t252 * t283 - t285;
t175 = t227 * t290 - t229 * t291;
t269 = -m(5) / 0.4e1 - m(6) / 0.4e1;
t357 = 0.2e1 * t269;
t381 = t175 * t357;
t187 = t227 * t254 + t256 * t229;
t380 = t187 * t357;
t374 = m(5) + m(6);
t379 = t374 * t175;
t378 = t374 * t187;
t364 = Icges(6,4) * t165 - Icges(6,5) * t203;
t104 = Icges(6,1) * t161 + t364;
t129 = Icges(6,4) * t160 + Icges(6,2) * t360;
t377 = t104 + t129;
t372 = Icges(6,4) * t352;
t278 = Icges(6,5) * t204 + t372;
t106 = Icges(6,1) * t163 + t278;
t128 = -Icges(6,2) * t163 + t372;
t376 = t106 + t128;
t362 = Icges(6,6) * t165 - Icges(6,3) * t203;
t97 = Icges(6,5) * t360 - t362;
t164 = (t224 * t256 + t250 * t286) * t255 + t229 * t253;
t280 = Icges(6,6) * t352 + Icges(6,3) * t204;
t99 = Icges(6,5) * t164 + t280;
t375 = t203 * t99 - t204 * t97;
t324 = -t203 / 0.2e1;
t292 = t249 * t251;
t223 = -t252 * t247 + t250 * t292;
t195 = t223 * t255 + t249 * t294;
t373 = Icges(6,4) * t195;
t222 = t247 * t292 + t289;
t262 = -t223 * t253 + t249 * t293;
t139 = t195 * rSges(6,1) + rSges(6,2) * t262 + t222 * rSges(6,3);
t369 = t203 * t139;
t318 = t248 * qJ(4) + pkin(2);
t319 = t249 * qJ(3) + pkin(1);
t210 = (t251 * pkin(3) + t318) * t252 + t319;
t263 = qJ(4) * t251 - qJ(2);
t232 = -t248 * pkin(3) + t263;
t368 = -t203 * rSges(5,2) - rSges(5,1) * (t228 * t250 + t247 * t291) - t227 * rSges(5,3) - t210 * t254 - t232 * t256;
t367 = t165 * rSges(6,2) - t203 * rSges(6,3);
t98 = Icges(6,5) * t163 + t280;
t366 = t165 * t102 - t203 * t98;
t96 = Icges(6,5) * t161 + t362;
t365 = t165 * t100 - t203 * t96;
t136 = Icges(6,5) * t195 + Icges(6,6) * t262 + Icges(6,3) * t222;
t137 = Icges(6,2) * t262 + Icges(6,6) * t222 + t373;
t361 = -t203 * t136 + t165 * t137;
t234 = t250 * pkin(4) + t247 * pkin(6) + pkin(3);
t184 = (t234 * t251 + t318) * t252 + pkin(1) + (t247 * pkin(4) - t250 * pkin(6) + qJ(3)) * t249;
t211 = -t234 * t248 + t263;
t277 = rSges(6,2) * t352 + t204 * rSges(6,3);
t359 = t184 * t256 - t211 * t254 + t277;
t358 = -t184 * t254 - t211 * t256;
t356 = t203 / 0.2e1;
t355 = t204 / 0.2e1;
t354 = t222 / 0.2e1;
t353 = Icges(6,4) * t262;
t231 = (-t254 ^ 2 - t256 ^ 2) * t249;
t350 = 0.2e1 * t231;
t349 = 2 * qJD(1);
t348 = 4 * qJD(1);
t347 = m(4) / 0.2e1;
t345 = m(5) / 0.2e1;
t344 = m(6) / 0.2e1;
t123 = (t230 * t250 + t247 * t290) * rSges(5,1) - t204 * rSges(5,2) + t229 * rSges(5,3) + t210 * t256 - t232 * t254;
t118 = t123 * t290;
t343 = m(5) * (-t291 * t368 + t118);
t342 = m(5) * (t123 * t254 + t256 * t368);
t108 = t161 * rSges(6,1) + t367;
t109 = rSges(6,1) * t360 - t367;
t156 = t163 * rSges(6,1);
t110 = t156 + t277;
t314 = t164 * rSges(6,1);
t111 = t277 + t314;
t29 = (-t110 + t111) * t204 + (t108 + t109) * t203;
t295 = t248 * t249;
t66 = -t222 * t108 - t369;
t298 = t204 * t139;
t69 = t222 * t110 - t298;
t34 = -t66 * t227 + t69 * t229;
t67 = t222 * t109 - t369;
t68 = -t222 * t111 + t298;
t341 = m(6) * (t67 * t227 + t68 * t229 + t29 * t295 + t34);
t62 = t69 * t290;
t340 = m(6) * (-t29 * t252 + t62 + (t256 * t68 + (-t66 + t67) * t254) * t249);
t42 = t69 * t254 + t66 * t256;
t339 = m(6) * (t68 * t254 - t67 * t256 + t42);
t338 = m(6) * t34;
t337 = m(6) * (-t66 * t291 + t62);
t336 = m(6) * t42;
t90 = t156 + t359;
t85 = t90 * t290;
t88 = -t108 + t358;
t335 = m(6) * (-t88 * t291 + t85);
t50 = t90 * t254 + t88 * t256;
t334 = m(6) * t50;
t132 = rSges(6,1) * t352 - rSges(6,2) * t163;
t133 = t160 * rSges(6,1) + rSges(6,2) * t360;
t333 = m(6) * (t227 * t132 - t229 * t133);
t332 = m(6) * (t132 * t254 - t133 * t256) * t249;
t331 = m(6) * (-t256 * t132 - t254 * t133);
t246 = t256 * qJ(2);
t322 = m(3) * ((rSges(3,2) * t291 + t256 * rSges(3,3) + t246) * t256 + (-rSges(3,2) * t290 + (rSges(3,3) + qJ(2)) * t254) * t254);
t264 = rSges(4,3) * t249 + pkin(2) * t252 + t319;
t149 = t230 * rSges(4,1) - t229 * rSges(4,2) + t254 * qJ(2) + t256 * t264;
t140 = t149 * t290;
t257 = -t228 * rSges(4,1) + t227 * rSges(4,2) - t254 * t264 + t246;
t321 = m(4) * (-t257 * t291 + t140);
t320 = m(4) * (t149 * t254 + t256 * t257);
t48 = -t88 * t227 + t90 * t229;
t317 = t100 * t352 + t204 * t96;
t316 = t102 * t352 + t204 * t98;
t315 = m(6) * qJD(5);
t312 = (t102 * t262 + t195 * t106 + t222 * t98) * t203;
t103 = Icges(6,4) * t164 + t279;
t107 = Icges(6,1) * t164 + t278;
t311 = (t103 * t262 + t195 * t107 + t222 * t99) * t203;
t138 = Icges(6,1) * t195 + Icges(6,5) * t222 + t353;
t303 = t262 * t138;
t146 = Icges(6,1) * t262 - t373;
t302 = t195 * t146;
t145 = -Icges(6,2) * t195 + t353;
t301 = t262 * t145;
t300 = t195 * t137;
t144 = Icges(6,5) * t262 - Icges(6,6) * t195;
t297 = t222 * t144;
t281 = t204 * t136 + t137 * t352;
t276 = t379 / 0.2e1;
t275 = -t379 / 0.2e1;
t273 = -t378 / 0.2e1;
t271 = t378 / 0.2e1;
t135 = (t344 + t345 + t347) * t350;
t270 = t135 * qJD(1);
t101 = -t363 + t382;
t105 = Icges(6,1) * t360 - t364;
t35 = t161 * t104 + t365;
t36 = t161 * t106 + t366;
t37 = t163 * t104 + t317;
t38 = t163 * t106 + t316;
t51 = t161 * t138 + t361;
t52 = t138 * t360 - t361;
t268 = t38 * t356 + (t101 * t352 + t163 * t105 + t106 * t360 + t36 - t366 + t37 - t375) * t355 + (t52 + t51) * t354 + (t103 * t352 + t104 * t360 + t163 * t107 + t35 - t365) * t324;
t28 = m(5) * (t123 * t229 - t227 * t368) + m(6) * t48;
t44 = t100 * t262 + t195 * t104 + t222 * t96;
t45 = t101 * t262 + t195 * t105 + t222 * t97;
t53 = t163 * t138 + t281;
t55 = t204 * t108 + t110 * t203;
t258 = -m(6) * (t55 * t29 + t66 * t68 + t69 * t67) + (-t203 * t37 + t38 * t204 + t53 * t222) * t324 - t222 * (-t311 + t312 + (t44 + t45) * t204) / 0.2e1;
t21 = t303 / 0.2e1 + t302 / 0.2e1 - t300 / 0.2e1 + t301 / 0.2e1 + t297 / 0.2e1 + m(6) * (t90 * t132 - t88 * t133);
t147 = rSges(6,1) * t262 - rSges(6,2) * t195;
t134 = (m(4) / 0.4e1 - t269) * t350 - (m(4) + t374) * t231 / 0.2e1;
t126 = Icges(6,5) * t352 - Icges(6,6) * t163;
t91 = -t314 - t359;
t89 = t109 + t358;
t84 = t222 * t132 - t204 * t147;
t83 = -t222 * t133 - t147 * t203;
t79 = t331 / 0.2e1;
t77 = t332 / 0.2e1;
t73 = t333 / 0.2e1;
t72 = t271 - t380;
t71 = t273 + t380;
t70 = t271 + t273;
t61 = t132 * t203 + t204 * t133;
t58 = t275 + t381;
t57 = t275 + t276;
t56 = t276 - t381;
t54 = t164 * t138 + t281;
t43 = (t297 - t300 + t301 + t302 + t303) * t222;
t41 = t336 / 0.2e1;
t39 = t337 / 0.2e1;
t33 = t338 / 0.2e1;
t32 = t137 * t360 + t160 * t138 - t144 * t203 + t165 * t145 + t161 * t146;
t31 = t204 * t144 + (t138 + t145) * t352 + (-t137 + t146) * t163;
t27 = t222 * t126 + t384 * t195 + t376 * t262;
t26 = t222 * t127 + t385 * t195 + t377 * t262;
t25 = t321 + t335 + t343;
t24 = t320 + t322 + t334 + t342;
t22 = t339 / 0.2e1;
t17 = t340 / 0.2e1;
t15 = t341 / 0.2e1;
t14 = t41 + t22 - t331 / 0.2e1;
t13 = t79 + t41 - t339 / 0.2e1;
t12 = t79 + t22 - t336 / 0.2e1;
t9 = t39 + t17 - t332 / 0.2e1;
t8 = t77 + t39 - t340 / 0.2e1;
t7 = t77 + t17 - t337 / 0.2e1;
t6 = t33 + t15 - t333 / 0.2e1;
t5 = t73 + t33 - t341 / 0.2e1;
t4 = t73 + t15 - t338 / 0.2e1;
t3 = t54 * t222 + (t165 * t101 + t161 * t105 + t164 * t106 + t316 + t35) * t204 + (-t165 * t103 - t164 * t104 - t161 * t107 - t317 + t36 + t375) * t203;
t1 = t268 * t204 + t3 * t324 - t258;
t2 = [m(6) * (t88 * t91 + t90 * t89) * t348 / 0.4e1 + t24 * qJD(2) + t25 * qJD(3) + t28 * qJD(4) + t21 * qJD(5), t24 * qJD(1) + t134 * qJD(3) + t70 * qJD(4) + t13 * qJD(5), t25 * qJD(1) + t134 * qJD(2) + t57 * qJD(4) + t8 * qJD(5), t28 * qJD(1) + t70 * qJD(2) + t57 * qJD(3) + t5 * qJD(5), t21 * qJD(1) + t13 * qJD(2) + t8 * qJD(3) + t5 * qJD(4) + (t43 + m(6) * (t69 * t132 - t66 * t133 + t83 * t88 + t84 * t90) + t258 - (t32 / 0.2e1 + t26 / 0.2e1 - t3 / 0.2e1) * t203 + (t27 / 0.2e1 + t31 / 0.2e1 - t268) * t204) * qJD(5); t135 * qJD(3) + t71 * qJD(4) + t12 * qJD(5) + (-t334 / 0.4e1 - t342 / 0.4e1 - t320 / 0.4e1 - t322 / 0.4e1) * t348 + (t254 * t91 - t256 * t89 + t50) * t344 * t349, 0, t270, t71 * qJD(1), t12 * qJD(1) + (t83 * t254 - t84 * t256) * t315; -t135 * qJD(2) + t58 * qJD(4) + t7 * qJD(5) + (-t335 / 0.4e1 - t343 / 0.4e1 - t321 / 0.4e1) * t348 + (t85 * t344 + t118 * t345 + t140 * t347 + ((-t123 * t345 - t149 * t347 + t91 * t344) * t256 + (-t88 + t89) * t344 * t254) * t249) * t349, -t270, 0, t58 * qJD(1), t7 * qJD(1) + (-t61 * t252 + (t254 * t84 + t256 * t83) * t249) * t315; (m(6) * (t227 * t89 + t229 * t91 + t48) - t28) * qJD(1) + t72 * qJD(2) + t56 * qJD(3) + t4 * qJD(5), t72 * qJD(1), t56 * qJD(1), 0, t4 * qJD(1) + (t84 * t227 + t83 * t229 + t61 * t295) * t315; t14 * qJD(2) + t9 * qJD(3) + t6 * qJD(4) + t1 * qJD(5) + (t54 * t324 + m(6) * (t66 * t91 + t67 * t90 + t68 * t88 + t69 * t89) + t312 / 0.2e1 + t53 * t356 - t311 / 0.2e1 - t21 + (t51 / 0.2e1 + t45 / 0.2e1 + t52 / 0.2e1 + t44 / 0.2e1) * t204) * qJD(1), t14 * qJD(1), t9 * qJD(1), t6 * qJD(1), t1 * qJD(1) + (m(6) * (t55 * t61 + t66 * t83 + t69 * t84) + (-t377 * t203 * t352 + t31 * t222 + (t204 * t126 + t376 * t352 - t386) * t204 + (-t385 * t203 + t384 * t204) * t163) * t355 + ((t102 * t360 + t160 * t106 - t126 * t203 + t165 * t128 + t161 * t130) * t204 - (t100 * t360 + t160 * t104 + t165 * t129 + t161 * t131 - t386) * t203 + t32 * t222) * t324 + (-t203 * t26 + t27 * t204 + t43) * t354) * qJD(5);];
Cq = t2;
