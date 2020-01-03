% Calculate vector of inverse dynamics joint torques for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:20
% EndTime: 2019-12-31 16:24:40
% DurationCPUTime: 18.45s
% Computational Cost: add. (7861->548), mult. (13383->839), div. (0->0), fcn. (12742->8), ass. (0->279)
t240 = sin(pkin(6));
t242 = cos(pkin(6));
t244 = sin(qJ(2));
t245 = cos(qJ(2));
t304 = Icges(3,5) * t245 - Icges(3,6) * t244;
t448 = (-Icges(3,3) * t242 + t240 * t304) * t242;
t236 = t240 ^ 2;
t237 = t242 ^ 2;
t428 = t236 + t237;
t208 = rSges(3,1) * t245 - rSges(3,2) * t244;
t298 = t428 * t208;
t445 = 2 * qJD(2);
t444 = 2 * qJDD(2);
t238 = pkin(7) + qJ(4);
t229 = sin(t238);
t400 = rSges(5,2) * t229;
t230 = cos(t238);
t402 = rSges(5,1) * t230;
t323 = -t400 + t402;
t141 = t244 * rSges(5,3) + t245 * t323;
t353 = qJD(4) * t244;
t358 = qJD(2) * t240;
t197 = t242 * t353 + t358;
t357 = qJD(2) * t242;
t198 = t240 * t353 - t357;
t381 = t240 * t245;
t154 = -t229 * t381 - t230 * t242;
t155 = -t242 * t229 + t230 * t381;
t382 = t240 * t244;
t63 = Icges(5,5) * t155 + Icges(5,6) * t154 + Icges(5,3) * t382;
t392 = Icges(5,4) * t155;
t65 = Icges(5,2) * t154 + Icges(5,6) * t382 + t392;
t144 = Icges(5,4) * t154;
t67 = Icges(5,1) * t155 + Icges(5,5) * t382 + t144;
t24 = t154 * t65 + t155 * t67 + t382 * t63;
t378 = t242 * t245;
t156 = -t229 * t378 + t230 * t240;
t157 = t229 * t240 + t230 * t378;
t379 = t242 * t244;
t64 = Icges(5,5) * t157 + Icges(5,6) * t156 + Icges(5,3) * t379;
t391 = Icges(5,4) * t157;
t66 = Icges(5,2) * t156 + Icges(5,6) * t379 + t391;
t145 = Icges(5,4) * t156;
t68 = Icges(5,1) * t157 + Icges(5,5) * t379 + t145;
t25 = t154 * t66 + t155 * t68 + t382 * t64;
t26 = t156 * t65 + t157 * t67 + t379 * t63;
t27 = t156 * t66 + t157 * t68 + t379 * t64;
t352 = qJD(4) * t245;
t302 = Icges(5,5) * t230 - Icges(5,6) * t229;
t133 = -Icges(5,3) * t245 + t244 * t302;
t389 = Icges(5,4) * t230;
t305 = -Icges(5,2) * t229 + t389;
t135 = -Icges(5,6) * t245 + t244 * t305;
t390 = Icges(5,4) * t229;
t308 = Icges(5,1) * t230 - t390;
t137 = -Icges(5,5) * t245 + t244 * t308;
t43 = t133 * t382 + t135 * t154 + t137 * t155;
t44 = t133 * t379 + t135 * t156 + t137 * t157;
t443 = t242 * (t197 * t27 + t198 * t26 - t352 * t44) + t240 * (t197 * t25 + t198 * t24 - t352 * t43);
t239 = sin(pkin(7));
t401 = rSges(4,2) * t239;
t241 = cos(pkin(7));
t403 = rSges(4,1) * t241;
t324 = -t401 + t403;
t426 = rSges(4,3) * t245 - t244 * t324;
t442 = t428 * qJD(2) * t426;
t206 = rSges(3,1) * t244 + rSges(3,2) * t245;
t441 = t206 * t428;
t351 = qJD(2) * qJD(3);
t438 = qJDD(3) * t244 + t245 * t351;
t356 = qJD(2) * t244;
t355 = qJD(2) * t245;
t429 = t245 * pkin(2) + t244 * qJ(3);
t427 = g(1) * t242 + g(2) * t240;
t319 = -t229 * t66 + t230 * t68;
t320 = -t229 * t65 + t230 * t67;
t425 = -(-t133 * t242 - t319) * t197 - (-t133 * t240 - t320) * t198;
t243 = -pkin(5) - qJ(3);
t376 = qJ(3) + t243;
t227 = pkin(3) * t241 + pkin(2);
t406 = pkin(2) - t227;
t423 = t244 * t406 - t245 * t376;
t167 = (-Icges(5,2) * t230 - t390) * t244;
t251 = t197 * (-Icges(5,2) * t157 + t145 + t68) + t198 * (-Icges(5,2) * t155 + t144 + t67) - t352 * (t137 + t167);
t421 = -m(4) - m(5);
t350 = qJD(2) * qJD(4);
t278 = qJDD(4) * t244 + t245 * t350;
t349 = qJDD(2) * t240;
t142 = t242 * t278 + t349;
t420 = t142 / 0.2e1;
t348 = qJDD(2) * t242;
t143 = t240 * t278 - t348;
t419 = t143 / 0.2e1;
t418 = -t428 * t356 / 0.2e1;
t417 = -t197 / 0.2e1;
t416 = t197 / 0.2e1;
t415 = -t198 / 0.2e1;
t414 = t198 / 0.2e1;
t202 = -qJDD(4) * t245 + t244 * t350;
t413 = t202 / 0.2e1;
t410 = -t245 / 0.2e1;
t409 = pkin(2) * t244;
t398 = t242 * t25;
t234 = t244 * rSges(4,3);
t397 = t26 * t240;
t396 = t43 * t244;
t395 = t44 * t244;
t385 = t133 * t245;
t384 = t239 * t240;
t383 = t239 * t242;
t380 = t241 * t245;
t377 = t245 * t243;
t354 = qJD(3) * t245;
t181 = qJD(2) * t429 - t354;
t273 = -t244 * t376 - t245 * t406;
t375 = -t273 * qJD(2) - t181;
t231 = qJD(3) * t244;
t216 = t240 * t231;
t205 = -qJ(3) * t245 + t409;
t287 = qJD(2) * t205;
t129 = -t240 * t287 + t216;
t218 = t242 * t231;
t130 = -t242 * t287 + t218;
t374 = t240 * t129 + t242 * t130;
t373 = -t205 + t423;
t370 = -(t245 * t324 + t234) * qJD(2) - t181;
t367 = -t205 + t426;
t194 = t429 * t240;
t196 = t429 * t242;
t366 = t240 * t194 + t242 * t196;
t343 = t244 * t400;
t365 = rSges(5,3) * t381 + t240 * t343;
t364 = rSges(5,3) * t378 + t242 * t343;
t344 = t244 * t401;
t363 = rSges(4,3) * t381 + t240 * t344;
t362 = rSges(4,3) * t378 + t242 * t344;
t361 = t438 * t240;
t360 = t438 * t242;
t346 = t244 * t403;
t345 = t244 * t402;
t170 = (-rSges(5,1) * t229 - rSges(5,2) * t230) * t244;
t75 = qJD(2) * t141 + qJD(4) * t170;
t342 = -t75 + t375;
t140 = -rSges(5,3) * t245 + t244 * t323;
t341 = -t140 + t373;
t220 = qJ(3) * t381;
t221 = qJ(3) * t378;
t340 = (-pkin(2) * t382 + t220) * t358 + (-pkin(2) * t379 + t221) * t357 + t231;
t339 = t240 * t356;
t338 = t240 * t355;
t337 = t242 * t356;
t336 = t242 * t355;
t329 = -t352 / 0.2e1;
t328 = t352 / 0.2e1;
t327 = t245 * t227 - t243 * t244;
t326 = qJD(2) * t373;
t325 = qJD(2) * t367;
t321 = t64 * t197 + t63 * t198;
t316 = t24 * t240 + t398;
t315 = t27 * t242 + t397;
t28 = t244 * t320 - t245 * t63;
t29 = t244 * t319 - t245 * t64;
t314 = t28 * t240 + t29 * t242;
t69 = rSges(5,1) * t155 + rSges(5,2) * t154 + rSges(5,3) * t382;
t70 = rSges(5,1) * t157 + rSges(5,2) * t156 + rSges(5,3) * t379;
t313 = -t240 * t70 + t242 * t69;
t182 = -t239 * t381 - t241 * t242;
t183 = t240 * t380 - t383;
t92 = rSges(4,1) * t183 + rSges(4,2) * t182 + rSges(4,3) * t382;
t184 = -t239 * t378 + t240 * t241;
t185 = t241 * t378 + t384;
t93 = rSges(4,1) * t185 + rSges(4,2) * t184 + rSges(4,3) * t379;
t312 = t240 * t92 + t242 * t93;
t94 = -pkin(3) * t383 + t240 * t273;
t95 = pkin(3) * t384 + t242 * t273;
t311 = t240 * t94 + t242 * t95;
t303 = -Icges(3,5) * t244 - Icges(3,6) * t245;
t301 = -t135 * t229 + t137 * t230;
t297 = qJD(2) * t441;
t294 = t194 * t358 + t196 * t357 + qJD(1) - t354;
t98 = -qJD(4) * t155 + t229 * t339;
t99 = qJD(4) * t154 - t230 * t339;
t49 = Icges(5,5) * t99 + Icges(5,6) * t98 + Icges(5,3) * t338;
t293 = t244 * t49 + t355 * t63;
t100 = -qJD(4) * t157 + t229 * t337;
t101 = qJD(4) * t156 - t230 * t337;
t50 = Icges(5,5) * t101 + Icges(5,6) * t100 + Icges(5,3) * t336;
t292 = t244 * t50 + t355 * t64;
t134 = Icges(5,3) * t244 + t245 * t302;
t166 = (-Icges(5,5) * t229 - Icges(5,6) * t230) * t244;
t72 = qJD(2) * t134 + qJD(4) * t166;
t291 = t133 * t355 + t244 * t72;
t290 = t134 - t301;
t23 = qJD(2) * t311 + t197 * t69 - t198 * t70 + t294;
t289 = t23 * t313;
t168 = (-Icges(5,1) * t229 - t389) * t244;
t275 = qJD(2) * t303;
t274 = -t227 * t244 - t377 + t409;
t270 = t166 * t352 - t197 * (Icges(5,5) * t156 - Icges(5,6) * t157) - t198 * (Icges(5,5) * t154 - Icges(5,6) * t155);
t269 = qJD(2) * t375 + qJDD(2) * t373;
t268 = qJD(2) * t370 + qJDD(2) * t367;
t267 = -qJDD(3) * t245 + t129 * t358 + t130 * t357 + t194 * t349 + t196 * t348 + t244 * t351 + qJDD(1);
t266 = qJD(2) * t423;
t265 = Icges(4,5) * t245 + (-Icges(4,1) * t241 + Icges(4,4) * t239) * t244;
t138 = Icges(5,5) * t244 + t245 * t308;
t263 = Icges(4,6) * t245 + (-Icges(4,4) * t241 + Icges(4,2) * t239) * t244;
t136 = Icges(5,6) * t244 + t245 * t305;
t259 = t244 * t270;
t257 = qJD(2) * t265;
t256 = qJD(2) * t263;
t252 = (Icges(5,1) * t156 - t391 - t66) * t197 + (Icges(5,1) * t154 - t392 - t65) * t198 - (-t135 + t168) * t352;
t247 = (-t290 * t352 - t425) * t244;
t219 = t242 * t354;
t217 = t240 * t354;
t195 = t206 * t242;
t193 = t206 * t240;
t188 = t303 * t242;
t187 = t303 * t240;
t174 = t242 * t275;
t173 = t240 * t275;
t147 = Icges(3,3) * t240 + t242 * t304;
t128 = t265 * t242;
t127 = t265 * t240;
t126 = t263 * t242;
t125 = t263 * t240;
t117 = -t242 * t345 + t364;
t116 = -t240 * t345 + t365;
t113 = t137 * t242;
t112 = t137 * t240;
t111 = t135 * t242;
t110 = t135 * t240;
t107 = t242 * t257;
t106 = t240 * t257;
t105 = t242 * t256;
t104 = t240 * t256;
t97 = t242 * t266;
t96 = t240 * t266;
t91 = rSges(5,1) * t156 - rSges(5,2) * t157;
t90 = rSges(5,1) * t154 - rSges(5,2) * t155;
t89 = Icges(4,1) * t185 + Icges(4,4) * t184 + Icges(4,5) * t379;
t88 = Icges(4,1) * t183 + Icges(4,4) * t182 + Icges(4,5) * t382;
t87 = Icges(4,4) * t185 + Icges(4,2) * t184 + Icges(4,6) * t379;
t86 = Icges(4,4) * t183 + Icges(4,2) * t182 + Icges(4,6) * t382;
t85 = Icges(4,5) * t185 + Icges(4,6) * t184 + Icges(4,3) * t379;
t84 = Icges(4,5) * t183 + Icges(4,6) * t182 + Icges(4,3) * t382;
t77 = t242 * t325 + t218;
t76 = t240 * t325 + t216;
t74 = qJD(2) * t138 + qJD(4) * t168;
t73 = qJD(2) * t136 + qJD(4) * t167;
t56 = rSges(5,1) * t101 + rSges(5,2) * t100 + rSges(5,3) * t336;
t55 = rSges(5,1) * t99 + rSges(5,2) * t98 + rSges(5,3) * t338;
t54 = Icges(5,1) * t101 + Icges(5,4) * t100 + Icges(5,5) * t336;
t53 = Icges(5,1) * t99 + Icges(5,4) * t98 + Icges(5,5) * t338;
t52 = Icges(5,4) * t101 + Icges(5,2) * t100 + Icges(5,6) * t336;
t51 = Icges(5,4) * t99 + Icges(5,2) * t98 + Icges(5,6) * t338;
t48 = t244 * t301 - t385;
t47 = -qJD(2) * t297 + qJDD(2) * t298 + qJDD(1);
t46 = t242 * t268 + t360;
t45 = t240 * t268 + t361;
t38 = qJD(2) * t312 + t294;
t37 = t140 * t198 + t242 * t326 + t352 * t69 + t218;
t36 = -t140 * t197 + t240 * t326 - t352 * t70 + t216;
t22 = qJD(2) * t442 + t312 * qJDD(2) + t267;
t17 = (qJD(2) * t301 - t72) * t245 + (qJD(2) * t133 - t229 * t73 + t230 * t74 + (-t135 * t230 - t137 * t229) * qJD(4)) * t244;
t16 = t100 * t135 + t101 * t137 + t156 * t73 + t157 * t74 + t242 * t291;
t15 = t135 * t98 + t137 * t99 + t154 * t73 + t155 * t74 + t240 * t291;
t14 = t140 * t143 + t198 * t75 - t202 * t69 + t242 * t269 + t352 * t55 + t360;
t13 = -t140 * t142 - t197 * t75 + t202 * t70 + t240 * t269 - t352 * t56 + t361;
t12 = t197 * t29 + t198 * t28 - t352 * t48;
t11 = t100 * t66 + t101 * t68 + t156 * t52 + t157 * t54 + t242 * t292;
t10 = t100 * t65 + t101 * t67 + t156 * t51 + t157 * t53 + t242 * t293;
t9 = t154 * t52 + t155 * t54 + t240 * t292 + t66 * t98 + t68 * t99;
t8 = t154 * t51 + t155 * t53 + t240 * t293 + t65 * t98 + t67 * t99;
t5 = (qJD(2) * t319 - t50) * t245 + (qJD(2) * t64 - t229 * t52 + t230 * t54 + (-t229 * t68 - t230 * t66) * qJD(4)) * t244;
t4 = (qJD(2) * t320 - t49) * t245 + (qJD(2) * t63 - t229 * t51 + t230 * t53 + (-t229 * t67 - t230 * t65) * qJD(4)) * t244;
t3 = t142 * t69 - t143 * t70 + t197 * t55 - t198 * t56 + t311 * qJDD(2) + (t240 * t96 + t242 * t97) * qJD(2) + t267;
t2 = t10 * t198 + t197 * t11 + t142 * t27 + t143 * t26 - t16 * t352 + t202 * t44;
t1 = t142 * t25 + t143 * t24 - t15 * t352 + t197 * t9 + t198 * t8 + t202 * t43;
t6 = [m(2) * qJDD(1) + m(3) * t47 + m(4) * t22 + m(5) * t3 + (-m(2) - m(3) + t421) * g(3); (t240 * t29 - t242 * t28) * t413 + (t240 * t9 - t242 * t8) * t414 + ((-t111 * t154 - t113 * t155) * t197 + (-t110 * t154 - t112 * t155) * t198 + (t396 + (-t136 * t154 - t138 * t155 + t398) * t245) * qJD(4) + (((t24 - t385) * qJD(4) + t321) * t245 + t247) * t240) * t415 + (-t10 * t242 + t11 * t240) * t416 + ((-t111 * t156 - t113 * t157) * t197 + (-t110 * t156 - t112 * t157) * t198 + (t395 + (-t136 * t156 - t138 * t157 + t397) * t245) * qJD(4) + (((t27 - t385) * qJD(4) + t321) * t245 + t247) * t242) * t417 + (-t24 * t242 + t240 * t25) * t419 + (t240 * t27 - t242 * t26) * t420 - t12 * t353 / 0.2e1 + (((t111 * t229 - t113 * t230 + t64) * t197 + (t110 * t229 - t112 * t230 + t63) * t198 + t48 * qJD(4)) * t244 + ((t290 * t245 + (t136 * t229 - t138 * t230 - t133) * t244 + t314) * qJD(4) + t425) * t245) * t328 - ((t126 * t184 + t128 * t185) * t358 + t188 * qJD(2) * t236 + (-t184 * t125 - t185 * t127 - t187 * t240) * t357) * t358 / 0.2e1 + (-(t125 * t182 + t127 * t183) * t357 + t187 * qJD(2) * t237 + (t182 * t126 + t183 * t128 - t188 * t242) * t358) * t357 / 0.2e1 + (t240 * t5 - t242 * t4 + t443) * t329 + (-t37 * (t141 * t198 + t219) - t36 * (-t141 * t197 + t217) - t23 * (t197 * t116 - t198 * t117 + t340) - ((-t37 * t327 + (t242 * t274 - t221) * t23) * t242 + (-t36 * t327 + (t240 * t274 - t220) * t23) * t240) * qJD(2) - ((t36 * t70 - t37 * t69) * t244 + (t37 * (t140 * t240 + t116) + t36 * (-t140 * t242 - t117) + t289) * t245) * qJD(4) + t3 * t366 + t23 * t374 + (t14 * t341 + t37 * t342 + t3 * (t70 + t95) + t23 * (t56 + t97)) * t242 + (t13 * t341 + t36 * t342 + t3 * (t69 + t94) + t23 * (t55 + t96)) * t240 - g(1) * t364 - g(2) * t365 - g(3) * (t141 + t327) - t427 * (-t377 + (-t227 - t402) * t244)) * m(5) + (t22 * t366 + (t22 * t93 + t367 * t46 + t370 * t77) * t242 + (t22 * t92 + t367 * t45 + t370 * t76) * t240 - g(1) * (t221 + t362) - g(2) * (t220 + t363) - t427 * t244 * (-pkin(2) - t403) - t77 * t219 - t76 * t217 + (-(-t76 * t240 - t77 * t242) * qJD(2) - g(3)) * (rSges(4,1) * t380 - t245 * t401 + t234 + t429) + (t374 - t340 - (t240 * (-t240 * t346 + t363) + t242 * (-t242 * t346 + t362)) * qJD(2) + t442) * t38) * m(4) + (g(1) * t195 + g(2) * t193 - g(3) * t208 + t47 * t298 + (-t297 - (-t193 * t240 - t195 * t242) * qJD(2)) * (qJD(2) * t298 + qJD(1)) + ((qJD(2) ^ 2) * t298 + qJDD(2) * t441) * t206) * m(3) + (t2 + ((t105 * t184 + t107 * t185 + t174 * t240) * t240 + (-t104 * t184 - t106 * t185 - t173 * t240) * t242) * t445 + ((-t184 * t86 - t185 * t88 - t379 * t84) * t242 + (t147 * t240 + t184 * t87 + t185 * t89 + t379 * t85 - t448) * t240) * t444) * t240 / 0.2e1 - (t1 + ((-t104 * t182 - t106 * t183 + t173 * t242) * t242 + (t105 * t182 + t107 * t183 - t174 * t242) * t240) * t445 + ((-t182 * t86 - t183 * t88 - t382 * t84 + t448) * t242 + (-t147 * t242 + t182 * t87 + t183 * t89 + t382 * t85) * t240) * t444) * t242 / 0.2e1; -t421 * g(3) * t245 + 0.2e1 * (t23 * t418 + t3 * t410) * m(5) + 0.2e1 * (t22 * t410 + t38 * t418) * m(4) + (t421 * t427 + m(4) * (qJD(2) * t38 + t240 * t45 + t242 * t46) + m(5) * (qJD(2) * t23 + t13 * t240 + t14 * t242)) * t244; t2 * t379 / 0.2e1 + (t244 * t315 - t245 * t44) * t420 + (-t16 * t245 + (t10 * t240 + t11 * t242) * t244 + (t245 * t315 + t395) * qJD(2)) * t416 + t1 * t382 / 0.2e1 + (t244 * t316 - t245 * t43) * t419 + (-t15 * t245 + (t240 * t8 + t242 * t9) * t244 + (t245 * t316 + t396) * qJD(2)) * t414 + t12 * t356 / 0.2e1 + (t142 * t29 + t143 * t28 - t17 * t352 + t197 * t5 + t198 * t4 + t202 * t48) * t410 + (t244 * t314 - t245 * t48) * t413 + (-t17 * t245 + (t240 * t4 + t242 * t5) * t244 + (t48 * t244 + t245 * t314) * qJD(2)) * t329 + (t156 * t251 + t157 * t252 - t242 * t259) * t417 + (t154 * t251 + t155 * t252 - t240 * t259) * t415 + (t270 * t245 + (-t229 * t251 + t252 * t230) * t244) * t328 + t443 * t355 / 0.2e1 + ((-t13 * t70 + t14 * t69 - t36 * t56 + t37 * t55 + (t289 + (t240 * t37 - t242 * t36) * t140) * qJD(2)) * t245 + (t37 * (-qJD(2) * t69 + t240 * t75) + t36 * (qJD(2) * t70 - t242 * t75) + t3 * t313 + t23 * (-t240 * t56 + t242 * t55) + (-t13 * t242 + t14 * t240) * t140) * t244 - t37 * (t170 * t198 + t352 * t90) - t36 * (-t170 * t197 - t352 * t91) - t23 * (t197 * t90 - t198 * t91) - g(1) * t91 - g(2) * t90 - g(3) * t170) * m(5);];
tau = t6;
