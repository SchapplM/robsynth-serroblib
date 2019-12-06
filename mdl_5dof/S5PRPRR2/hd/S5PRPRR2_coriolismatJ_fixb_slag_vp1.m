% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:48
% EndTime: 2019-12-05 15:45:00
% DurationCPUTime: 6.79s
% Computational Cost: add. (36204->382), mult. (33513->582), div. (0->0), fcn. (35991->10), ass. (0->253)
t269 = sin(pkin(8));
t266 = t269 ^ 2;
t270 = cos(pkin(8));
t267 = t270 ^ 2;
t434 = t266 + t267;
t268 = qJ(2) + pkin(9);
t263 = qJ(4) + t268;
t258 = sin(t263);
t259 = cos(t263);
t272 = sin(qJ(5));
t274 = cos(qJ(5));
t319 = rSges(6,1) * t274 - rSges(6,2) * t272;
t193 = -rSges(6,3) * t259 + t319 * t258;
t177 = t193 * t269;
t178 = t193 * t270;
t241 = rSges(5,1) * t258 + rSges(5,2) * t259;
t157 = t434 * t241;
t261 = sin(t268);
t262 = cos(t268);
t303 = -Icges(4,5) * t261 - Icges(4,6) * t262;
t302 = -Icges(5,5) * t258 - Icges(5,6) * t259;
t273 = sin(qJ(2));
t275 = cos(qJ(2));
t433 = 0.2e1 * t273 * (Icges(3,1) - Icges(3,2)) * t275 + (-0.2e1 * t273 ^ 2 + 0.2e1 * t275 ^ 2) * Icges(3,4);
t402 = t269 / 0.2e1;
t401 = -t270 / 0.2e1;
t304 = -Icges(3,5) * t273 - Icges(3,6) * t275;
t245 = t304 * t269;
t246 = t304 * t270;
t432 = t303 * t269 + t245;
t431 = t303 * t270 + t246;
t339 = qJD(2) + qJD(4);
t366 = (-Icges(6,5) * t272 - Icges(6,6) * t274) * t258 * t259;
t375 = Icges(6,4) * t274;
t305 = -Icges(6,2) * t272 + t375;
t189 = -Icges(6,6) * t259 + t305 * t258;
t350 = -t189 + (-Icges(6,1) * t272 - t375) * t258;
t376 = Icges(6,4) * t272;
t312 = Icges(6,1) * t274 - t376;
t191 = -Icges(6,5) * t259 + t312 * t258;
t349 = t191 + (-Icges(6,2) * t274 - t376) * t258;
t242 = rSges(5,1) * t259 - rSges(5,2) * t258;
t142 = t434 * t242;
t111 = (-t142 + t242) * t157;
t356 = t270 * t274;
t359 = t269 * t272;
t237 = -t259 * t359 - t356;
t357 = t270 * t272;
t358 = t269 * t274;
t238 = t259 * t358 - t357;
t363 = t258 * t269;
t151 = rSges(6,1) * t238 + rSges(6,2) * t237 + rSges(6,3) * t363;
t239 = -t259 * t357 + t358;
t240 = t259 * t356 + t359;
t362 = t258 * t270;
t152 = rSges(6,1) * t240 + rSges(6,2) * t239 + rSges(6,3) * t362;
t244 = pkin(4) * t259 + pkin(7) * t258;
t114 = t269 * t151 + t270 * t152 + t244 * t434;
t243 = pkin(4) * t258 - pkin(7) * t259;
t119 = -t269 * t177 - t270 * t178 - t243 * t434;
t348 = -t193 - t243;
t153 = t348 * t269;
t194 = rSges(6,3) * t258 + t319 * t259;
t347 = -t194 - t244;
t154 = t347 * t269;
t155 = t348 * t270;
t156 = t347 * t270;
t56 = t114 * t119 + t153 * t154 + t155 * t156;
t422 = m(5) * t111 + m(6) * t56;
t397 = pkin(2) * t273;
t254 = -pkin(3) * t261 - t397;
t290 = t254 + t348;
t137 = t290 * t269;
t139 = t290 * t270;
t264 = t275 * pkin(2);
t346 = t434 * t264;
t396 = pkin(3) * t262;
t322 = t396 * t434 + t346;
t89 = t114 + t322;
t37 = t89 * t119 + t137 * t154 + t139 * t156;
t105 = t322 + t142;
t295 = -t241 + t254;
t179 = t295 * t269;
t181 = t295 * t270;
t70 = -t105 * t157 + (-t179 * t269 - t181 * t270) * t242;
t421 = -m(5) * t70 - m(6) * t37;
t160 = Icges(6,5) * t237 - Icges(6,6) * t238;
t225 = Icges(6,4) * t237;
t149 = Icges(6,1) * t238 + Icges(6,5) * t363 + t225;
t352 = -Icges(6,2) * t238 + t149 + t225;
t378 = Icges(6,4) * t238;
t147 = Icges(6,2) * t237 + Icges(6,6) * t363 + t378;
t354 = Icges(6,1) * t237 - t147 - t378;
t79 = t160 * t363 + t352 * t237 + t354 * t238;
t161 = Icges(6,5) * t239 - Icges(6,6) * t240;
t226 = Icges(6,4) * t239;
t150 = Icges(6,1) * t240 + Icges(6,5) * t362 + t226;
t351 = -Icges(6,2) * t240 + t150 + t226;
t377 = Icges(6,4) * t240;
t148 = Icges(6,2) * t239 + Icges(6,6) * t362 + t377;
t353 = Icges(6,1) * t239 - t148 - t377;
t80 = t161 * t363 + t351 * t237 + t353 * t238;
t45 = t269 * t80 - t270 * t79;
t81 = t160 * t362 + t352 * t239 + t354 * t240;
t82 = t161 * t362 + t351 * t239 + t353 * t240;
t46 = t269 * t82 - t270 * t81;
t393 = t401 * t45 + t402 * t46;
t301 = Icges(6,5) * t274 - Icges(6,6) * t272;
t187 = -Icges(6,3) * t259 + t301 * t258;
t419 = t434 * t303;
t418 = t434 * t302;
t415 = 2 * qJD(2);
t414 = m(4) / 0.2e1;
t413 = m(5) / 0.2e1;
t412 = m(6) / 0.2e1;
t297 = t151 * t270 - t152 * t269;
t120 = t297 * t258;
t131 = t151 * t259 + t193 * t363;
t132 = -t152 * t259 - t193 * t362;
t333 = t120 * t119 + t131 * t156 + t132 * t154;
t112 = (t269 * t194 - t151) * t258;
t113 = (-t270 * t194 + t152) * t258;
t99 = t297 * t259 + (-t177 * t270 + t178 * t269) * t258;
t336 = t112 * t139 + t113 * t137 + t99 * t89;
t410 = m(6) * (t333 + t336);
t21 = t112 * t155 + t113 * t153 + t114 * t99 + t333;
t409 = m(6) * t21;
t224 = (-rSges(6,1) * t272 - rSges(6,2) * t274) * t258;
t166 = rSges(6,1) * t237 - rSges(6,2) * t238;
t167 = rSges(6,1) * t239 - rSges(6,2) * t240;
t134 = t166 * t269 + t167 * t270;
t61 = t89 * t134;
t94 = t114 * t134;
t408 = m(6) * (t61 + t94 + ((-t139 - t155) * t270 + (-t137 - t153) * t269) * t224);
t407 = m(6) * (t112 * t131 + t113 * t132 + t120 * t99);
t404 = m(6) * (t112 * t269 - t113 * t270);
t403 = -t259 / 0.2e1;
t400 = t270 / 0.2e1;
t398 = m(6) * (-t154 * t270 + t156 * t269);
t335 = -t404 / 0.2e1;
t87 = 0.2e1 * (t99 / 0.4e1 - t134 / 0.4e1) * m(6);
t355 = t87 * qJD(1);
t392 = qJD(3) * t335 - t355;
t334 = t404 / 0.2e1;
t391 = t339 * t334;
t390 = m(6) * qJD(5);
t108 = -m(5) * t157 + m(6) * t119;
t88 = (t134 + t99) * t412;
t386 = t108 * qJD(4) + t88 * qJD(5);
t385 = qJD(4) * t398 + qJD(5) * t334;
t146 = Icges(6,5) * t240 + Icges(6,6) * t239 + Icges(6,3) * t362;
t102 = t146 * t363 + t148 * t237 + t150 * t238;
t371 = t102 * t270;
t145 = Icges(6,5) * t238 + Icges(6,6) * t237 + Icges(6,3) * t363;
t103 = t145 * t362 + t147 * t239 + t149 * t240;
t370 = t103 * t269;
t369 = t145 * t259;
t368 = t146 * t259;
t367 = t187 * t259;
t361 = t259 * t269;
t360 = t259 * t270;
t101 = t145 * t363 + t147 * t237 + t149 * t238;
t117 = t187 * t363 + t189 * t237 + t191 * t238;
t190 = Icges(6,6) * t258 + t305 * t259;
t192 = Icges(6,5) * t258 + t312 * t259;
t296 = -t189 * t272 + t191 * t274;
t287 = (Icges(6,3) * t258 + t301 * t259 - t296) * t259;
t172 = t189 * t269;
t174 = t191 * t269;
t299 = -t147 * t272 + t149 * t274;
t293 = -t187 * t269 - t299;
t277 = t293 * t258 + t369;
t71 = -t172 * t237 - t174 * t238 + t277 * t269;
t173 = t189 * t270;
t175 = t191 * t270;
t298 = -t148 * t272 + t150 * t274;
t292 = -t187 * t270 - t298;
t276 = t292 * t258 + t368;
t72 = -t173 * t237 - t175 * t238 + t276 * t269;
t17 = (t371 - t190 * t237 - t192 * t238 + (t101 - t367) * t269) * t259 + (t72 * t270 + t117 + (t71 - t287) * t269) * t258;
t104 = t146 * t362 + t148 * t239 + t150 * t240;
t118 = t187 * t362 + t189 * t239 + t191 * t240;
t73 = -t172 * t239 - t174 * t240 + t277 * t270;
t74 = -t173 * t239 - t175 * t240 + t276 * t270;
t18 = (t370 - t190 * t239 - t192 * t240 + (t104 - t367) * t270) * t259 + (t73 * t269 + t118 + (t74 - t287) * t270) * t258;
t125 = t296 * t258 - t367;
t109 = t299 * t258 - t369;
t110 = t298 * t258 - t368;
t300 = t109 * t269 + t110 * t270;
t77 = -t293 * t259 + (t172 * t272 - t174 * t274 + t145) * t258;
t78 = -t292 * t259 + (t173 * t272 - t175 * t274 + t146) * t258;
t22 = (t287 + t300) * t259 + (t78 * t270 + t77 * t269 - (-t190 * t272 + t192 * t274 + t187) * t259 + t125) * t258;
t50 = -t117 * t259 + (t101 * t269 + t371) * t258;
t51 = -t118 * t259 + (t104 * t270 + t370) * t258;
t55 = -t125 * t259 + t300 * t258;
t3 = t407 + (t51 * t400 + t50 * t402 - t22 / 0.2e1) * t259 + (t18 * t400 + t17 * t402 + t55 / 0.2e1) * t258;
t338 = qJD(3) * t334 + t3 * qJD(5) + t355;
t337 = t408 / 0.2e1 + t393;
t332 = t363 / 0.2e1;
t331 = t362 / 0.2e1;
t251 = rSges(4,1) * t261 + rSges(4,2) * t262;
t330 = -t251 - t397;
t252 = rSges(4,1) * t262 - rSges(4,2) * t261;
t329 = -t252 - t264;
t215 = t302 * t269;
t216 = t302 * t270;
t35 = t269 * t72 - t270 * t71;
t36 = t269 * t74 - t270 * t73;
t327 = (t36 + t266 * t216 + (-t269 * t215 + t418) * t270) * t402 + (t35 + t267 * t215 + (-t270 * t216 + t418) * t269) * t401;
t321 = t434 * t397;
t320 = -t264 - t396;
t255 = t273 * rSges(3,1) + rSges(3,2) * t275;
t294 = -t242 + t320;
t141 = -t251 * t434 - t321;
t291 = t141 * t252;
t289 = t320 + t347;
t288 = -t321 + t434 * (t254 + t397);
t286 = t120 * t134 + (-t131 * t270 - t132 * t269) * t224;
t285 = t433 * t269 + t246;
t284 = -t433 * t270 + t245;
t283 = t17 * t401 + t18 * t402 + t36 * t331 + t35 * t332 + (t269 * t78 - t270 * t77) * t403 + (-t101 * t270 + t102 * t269) * t361 / 0.2e1 + (-t103 * t270 + t104 * t269) * t360 / 0.2e1 + t258 * (-t109 * t270 + t110 * t269) / 0.2e1 - t393;
t29 = -(t349 * t237 + t350 * t238) * t259 + (t80 * t270 + (t79 - t366) * t269) * t258;
t30 = -(t349 * t239 + t350 * t240) * t259 + (t81 * t269 + (t82 - t366) * t270) * t258;
t91 = -t160 * t259 + (-t352 * t272 + t354 * t274) * t258;
t92 = -t161 * t259 + (-t351 * t272 + t353 * t274) * t258;
t278 = -t17 * t363 / 0.2e1 - t18 * t362 / 0.2e1 + t259 * t22 / 0.2e1 - t407 + t30 * t402 + t29 * t401 + t45 * t332 + t46 * t331 - t50 * t361 / 0.2e1 - t51 * t360 / 0.2e1 + (t269 * t92 - t270 * t91) * t403 - t258 * t55 / 0.2e1;
t210 = t329 * t270;
t209 = t329 * t269;
t182 = t294 * t270;
t180 = t294 * t269;
t176 = t434 * t255;
t140 = t289 * t270;
t138 = t289 * t269;
t136 = -t167 * t259 - t224 * t362;
t135 = t166 * t259 + t224 * t363;
t130 = (t166 * t270 - t167 * t269) * t258;
t124 = t288 - t157;
t106 = t288 + t119;
t85 = t87 * qJD(5);
t64 = qJD(5) * t335;
t60 = t94 + (-t153 * t269 - t155 * t270) * t224;
t49 = t61 + (-t137 * t269 - t139 * t270) * t224;
t20 = t409 / 0.2e1;
t11 = t410 / 0.2e1;
t10 = m(6) * t60 + t393;
t9 = m(6) * t49 + t393;
t8 = t327 + t422;
t7 = t8 * qJD(4);
t6 = t327 - t421;
t5 = t20 - t410 / 0.2e1 + t337;
t4 = t11 - t409 / 0.2e1 + t337;
t1 = t11 + t20 - t408 / 0.2e1 + t283;
t2 = [0, (-m(3) * t176 / 0.2e1 + t141 * t414 + t124 * t413 + t106 * t412) * t415 + t386, 0, qJD(2) * t108 + t386, t130 * t390 + t339 * t88; -t85, (m(6) * (t106 * t89 + t137 * t138 + t139 * t140) + m(5) * (t105 * t124 + t179 * t180 + t181 * t182) + m(4) * (t346 * t141 + (t330 * t210 + t270 * t291) * t270 + (t330 * t209 + t269 * t291) * t269) + t327 + m(3) * (-t176 + t255) * t434 * (rSges(3,1) * t275 - t273 * rSges(3,2)) + ((t285 * t270 + t419 + (t284 - t432) * t269) * t270 + t431 * t266) * t402 + ((t284 * t269 + t419 + (t285 - t431) * t270) * t269 + t432 * t267) * t401) * qJD(2) + t6 * qJD(4) + t9 * qJD(5), t64, t6 * qJD(2) + t4 * qJD(5) + (t327 + 0.2e1 * (t56 + t37) * t412 + 0.2e1 * (t70 + t111) * t413 - t422) * qJD(4), t9 * qJD(2) + t4 * qJD(4) + (t278 + m(6) * (t130 * t89 + t135 * t139 + t136 * t137 + t286)) * qJD(5) + t392; 0, ((-t138 * t270 + t140 * t269) * t412 + (-t180 * t270 + t182 * t269) * t413 + (-t209 * t270 + t210 * t269) * t414) * t415 + t385, 0, qJD(2) * t398 + t385, (t135 * t269 - t136 * t270) * t390 + t391; -t85, t7 + t5 * qJD(5) + ((t106 * t114 + t138 * t153 + t140 * t155 + t37) * t412 + (t124 * t142 + (-t180 * t269 - t182 * t270) * t241 + t70) * t413) * t415 + (t327 + t421) * qJD(2), t64, qJD(2) * t8 + qJD(5) * t10 + t7, t5 * qJD(2) + t10 * qJD(4) + (t278 + m(6) * (t114 * t130 + t135 * t155 + t136 * t153 + t286)) * qJD(5) + t392; t339 * t87, ((t106 * t120 + t131 * t140 + t132 * t138 + t336 - t49) * m(6) + t283) * qJD(2) + t1 * qJD(4) + t338, t391, t1 * qJD(2) + ((t21 - t60) * m(6) + t283) * qJD(4) + t338, t339 * t3 + (m(6) * (t120 * t130 + t131 * t135 + t132 * t136) - t259 ^ 2 * t366 / 0.2e1 + (t30 * t400 + t29 * t402 + (t92 * t270 + t91 * t269 - (-t349 * t272 + t350 * t274) * t259) * t403) * t258) * qJD(5);];
Cq = t2;
