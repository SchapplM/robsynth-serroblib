% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:30
% EndTime: 2019-12-31 17:56:34
% DurationCPUTime: 3.25s
% Computational Cost: add. (14974->248), mult. (15169->342), div. (0->0), fcn. (17654->8), ass. (0->182)
t246 = qJ(1) + pkin(8);
t234 = sin(t246);
t235 = cos(t246);
t307 = sin(qJ(4));
t308 = cos(qJ(4));
t172 = -t234 * t308 + t235 * t307;
t171 = -t234 * t307 - t235 * t308;
t196 = sin(qJ(5));
t198 = cos(qJ(5));
t264 = t198 * Icges(6,4);
t219 = -t196 * Icges(6,2) + t264;
t127 = Icges(6,6) * t171 + t219 * t172;
t285 = Icges(6,4) * t196;
t221 = Icges(6,1) * t198 - t285;
t131 = Icges(6,5) * t171 + t221 * t172;
t363 = t127 * t198 + t131 * t196;
t287 = t363 * t172;
t126 = -Icges(6,6) * t172 + t219 * t171;
t130 = -Icges(6,5) * t172 + t221 * t171;
t364 = t126 * t198 + t130 * t196;
t384 = t364 * t171;
t385 = t287 / 0.2e1 - t384 / 0.2e1;
t383 = -t172 / 0.4e1;
t188 = -rSges(6,1) * t196 - rSges(6,2) * t198;
t148 = t172 * t188;
t149 = t188 * t171;
t379 = m(6) * (t234 * t148 + t149 * t235);
t244 = -t379 / 0.2e1;
t91 = t379 / 0.2e1;
t329 = m(6) / 0.2e1;
t331 = m(5) / 0.2e1;
t150 = t171 * rSges(5,1) + t172 * rSges(5,2);
t223 = -t172 * rSges(5,1) + t171 * rSges(5,2);
t341 = t150 * t234 + t235 * t223;
t166 = t171 * pkin(7);
t267 = t172 * t196;
t163 = rSges(6,2) * t267;
t165 = t171 * rSges(6,3);
t236 = -t163 + t165;
t290 = rSges(6,1) * t198;
t238 = pkin(4) + t290;
t95 = -t238 * t172 - t166 - t236;
t89 = t235 * t95;
t271 = t171 * t196;
t356 = rSges(6,2) * t271 + t172 * rSges(6,3);
t93 = -t172 * pkin(7) + t238 * t171 - t356;
t348 = t93 * t234 + t89;
t296 = t348 * t329 + t341 * t331;
t377 = t91 + t244;
t381 = qJD(5) * t377;
t266 = t172 * t198;
t162 = rSges(6,1) * t266;
t249 = t162 + t165;
t289 = rSges(6,2) * t196;
t94 = (-pkin(4) + t289) * t172 - t166 - t249;
t292 = t94 - t95;
t349 = m(6) * t292;
t317 = t93 * t349;
t380 = t317 * qJD(4);
t261 = t377 * qJD(3);
t291 = m(6) * qJD(5);
t134 = -t236 - t162;
t56 = (t134 - t163 + t249) * t171;
t378 = t56 * t291;
t376 = t126 * t196;
t375 = t127 * t196;
t374 = t130 * t198;
t373 = t131 * t198;
t217 = Icges(6,5) * t198 - Icges(6,6) * t196;
t358 = t217 * t172;
t124 = -Icges(6,3) * t171 - t358;
t372 = t172 * t124;
t359 = t217 * t171;
t125 = Icges(6,3) * t172 - t359;
t371 = t172 * t125;
t368 = t56 * qJD(2);
t367 = t329 * qJD(1);
t207 = -sin(qJ(1)) * pkin(1) - t234 * pkin(2) + t235 * qJ(3);
t202 = -t234 * pkin(3) + t207;
t77 = t202 - t95;
t203 = cos(qJ(1)) * pkin(1) + t235 * pkin(2) + t234 * qJ(3);
t201 = t235 * pkin(3) + t203;
t78 = t201 - t93;
t366 = -t77 * t93 + t78 * t95;
t260 = t171 * t124 - t131 * t266;
t343 = -t125 - t375;
t365 = -t343 * t172 + t260;
t362 = -t196 / 0.2e1;
t361 = -t198 / 0.2e1;
t216 = Icges(6,5) * t196 + Icges(6,6) * t198;
t141 = t216 * t171;
t360 = t216 * t172;
t357 = 0.2e1 * t367;
t297 = -m(6) * t348 / 0.2e1 - m(5) * t341 / 0.2e1;
t298 = (-t94 * t235 + t89) * t329;
t355 = t297 - t298;
t353 = 0.2e1 * t244;
t138 = -t150 + t201;
t200 = t202 - t223;
t352 = t138 * t223 - t200 * t150;
t282 = Icges(6,2) * t198;
t218 = t282 + t285;
t210 = t196 * t218;
t220 = Icges(6,1) * t196 + t264;
t262 = t198 * t220;
t342 = t210 - t262;
t100 = t342 * t171 + t360;
t351 = t196 / 0.2e1;
t350 = (-t221 + t282) * t361 + (t219 + t220 + t264) * t362;
t301 = m(6) * t188;
t346 = t172 * t171;
t344 = -t124 + t376;
t101 = t342 * t172 - t141;
t269 = t172 * t101;
t272 = t171 * t100;
t340 = t269 / 0.4e1 + t272 / 0.4e1;
t339 = t384 / 0.4e1 + t148 * t349 / 0.2e1;
t338 = t287 / 0.4e1 - ((-t77 - t95) * t172 + (t78 + t93) * t171) * t301 / 0.2e1;
t230 = t210 / 0.2e1 + t221 * t362 + t219 * t361 - t262 / 0.2e1;
t337 = t171 ^ 2;
t336 = t172 ^ 2;
t334 = 0.4e1 * qJD(1);
t333 = 0.2e1 * qJD(4);
t36 = -t127 * t267 - t260;
t259 = -t171 * t125 + t130 * t266;
t37 = -t126 * t267 + t259;
t19 = -t171 * t36 + t172 * t37;
t328 = -t19 / 0.2e1;
t270 = t171 * t198;
t258 = t131 * t270 + t372;
t38 = -t127 * t271 + t258;
t257 = t130 * t270 + t371;
t39 = -t126 * t271 + t257;
t20 = -t171 * t38 + t172 * t39;
t327 = t20 / 0.2e1;
t326 = m(5) * t352;
t323 = m(5) * (t138 * t234 + t200 * t235);
t322 = m(6) * ((t78 - t93) * t149 + (-t77 + t95) * t148);
t319 = m(6) * t366;
t316 = m(6) * (t148 * t77 - t149 * t78);
t315 = m(6) * (-t148 * t95 + t149 * t93);
t75 = t77 * t235;
t314 = m(6) * (t78 * t234 + t75);
t306 = m(4) * ((t235 * rSges(4,3) + t207) * t235 + (t234 * rSges(4,3) + t203) * t234);
t295 = t56 * t367;
t79 = t202 - t94;
t294 = t77 - t79;
t274 = t148 * t188;
t273 = t149 * t188;
t248 = qJD(5) * t171;
t247 = qJD(5) * t172;
t67 = t384 / 0.2e1;
t229 = t328 + (-t220 * t171 - t126) * t351 + (t218 * t171 - t130) * t361 + t171 * t350 - t358 / 0.2e1;
t228 = t327 + (-t220 * t172 - t127) * t351 + (-t218 * t172 + t131) * t198 / 0.2e1 + t359 / 0.2e1 + t172 * t350;
t215 = t230 + t385;
t214 = -t287 / 0.2e1 + t67 + t230;
t212 = -t269 / 0.4e1 - t338 + t340 - (-t364 + t100) * t171 / 0.4e1;
t211 = -t272 / 0.4e1 - t339 + t340 + (-t363 + t101) * t383;
t190 = t289 - t290;
t80 = -t148 * t172 - t149 * t171;
t58 = t134 * t172 + (-rSges(6,1) * t270 + t356) * t171;
t41 = 0.2e1 * t91;
t32 = -t230 + t315;
t31 = -t230 + t316;
t28 = t56 * t58;
t24 = t306 + t314 + t323;
t23 = (-t126 * t172 - t127 * t171) * t196 + t258 + t259;
t14 = t322 / 0.2e1;
t13 = t319 + t326;
t12 = t296 - t355;
t11 = t297 + t298 - t296;
t10 = t296 + t355;
t9 = (t36 + (-t373 + t375) * t172 + t257) * t172 + (-t172 * t344 - t23 + t37) * t171;
t8 = (t23 - t38) * t172 + (-t39 + (t374 - t376) * t171 + t365) * t171;
t7 = (t38 + (t344 - t374) * t172) * t172 + (t344 * t171 - t257 + t371 + t39) * t171;
t6 = (-t36 - t365) * t172 + (-t37 + (t343 + t373) * t171 + t372) * t171;
t5 = t211 - t322 / 0.2e1 + t212 - t230;
t4 = t363 * t383 + t14 + t212 + t215 + t339;
t3 = t211 + t14 - t384 / 0.4e1 + t214 + t338;
t2 = (t8 / 0.2e1 + t328) * t172 + (-t20 / 0.2e1 - t6 / 0.2e1) * t171;
t1 = m(6) * t28 + (t7 / 0.2e1 + t19 / 0.2e1) * t172 + (t327 - t9 / 0.2e1) * t171;
t15 = [-m(6) * t294 * t78 * t334 / 0.4e1 + t24 * qJD(3) + t13 * qJD(4) + t31 * qJD(5), -t378 / 0.2e1, qJD(1) * t24 + qJD(4) * t10 + t381, t13 * qJD(1) + t10 * qJD(3) + t4 * qJD(5) - t380 + (-t366 * t329 - t352 * t331) * t333, t31 * qJD(1) + t261 + t4 * qJD(4) + (-t7 / 0.2e1 + t229) * t247 + (t9 / 0.2e1 - t228) * t248 + (-t368 / 0.2e1 + ((-t190 * t78 + t273) * t172 + (-t190 * t77 - t274) * t171 - t28) * qJD(5)) * m(6); t378 / 0.2e1, 0, 0, 0, t80 * t291 + t295; t11 * qJD(4) + t41 * qJD(5) + (-t314 / 0.4e1 - t323 / 0.4e1 - t306 / 0.4e1) * t334 + (-t235 * t79 + t75) * t357, 0, 0, t11 * qJD(1) + t353 * qJD(5) + t296 * t333, t41 * qJD(1) + t353 * qJD(4) + (-t171 * t234 + t172 * t235) * t190 * t291; t12 * qJD(3) + t380 + t3 * qJD(5) + (-t319 / 0.4e1 - t326 / 0.4e1) * t334 + (t292 * t78 - t294 * t93) * t357, 0, qJD(1) * t12 - t381, qJD(1) * t317 + t32 * qJD(5), t3 * qJD(1) - t261 + t32 * qJD(4) + (-t8 / 0.2e1 - t229) * t247 + (t6 / 0.2e1 + t228) * t248 + ((-t190 * t93 - t273) * t172 + (-t190 * t95 + t274) * t171) * t291; t353 * qJD(3) + t5 * qJD(4) + t1 * qJD(5) + t329 * t368 + (t215 + t67 - t316 + (-t363 / 0.2e1 + t294 * t301) * t172) * qJD(1), t295, qJD(1) * t353 + qJD(4) * t377, t5 * qJD(1) + t2 * qJD(5) + t261 + (t214 - t315 + t385) * qJD(4), t1 * qJD(1) + t2 * qJD(4) + (m(6) * (t58 * t80 + (t336 + t337) * t190 * t188) + t172 * (t336 * t141 - t346 * t360) / 0.2e1 - t171 * (-t141 * t346 + t337 * t360) / 0.2e1) * qJD(5);];
Cq = t15;
