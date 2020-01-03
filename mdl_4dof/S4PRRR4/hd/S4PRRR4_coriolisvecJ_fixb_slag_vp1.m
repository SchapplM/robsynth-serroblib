% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:30
% EndTime: 2019-12-31 16:32:37
% DurationCPUTime: 5.40s
% Computational Cost: add. (9192->462), mult. (8462->643), div. (0->0), fcn. (6680->6), ass. (0->278)
t211 = qJ(3) + qJ(4);
t207 = cos(t211);
t202 = Icges(5,4) * t207;
t206 = sin(t211);
t153 = Icges(5,1) * t206 + t202;
t251 = -Icges(5,2) * t206 + t202;
t391 = t153 + t251;
t209 = pkin(7) + qJ(2);
t204 = sin(t209);
t331 = t204 * t206;
t164 = Icges(5,4) * t331;
t330 = t204 * t207;
t205 = cos(t209);
t347 = Icges(5,5) * t205;
t103 = Icges(5,1) * t330 - t164 - t347;
t343 = Icges(5,6) * t205;
t101 = Icges(5,4) * t330 - Icges(5,2) * t331 - t343;
t340 = t101 * t206;
t250 = -t103 * t207 + t340;
t237 = t250 * t204;
t150 = Icges(5,5) * t207 - Icges(5,6) * t206;
t100 = Icges(5,3) * t204 + t150 * t205;
t349 = Icges(5,4) * t206;
t154 = Icges(5,1) * t207 - t349;
t104 = Icges(5,5) * t204 + t154 * t205;
t325 = t205 * t207;
t365 = t100 * t204 + t104 * t325;
t390 = -t237 - t365;
t389 = 2 * qJD(3);
t148 = pkin(2) * t205 + t204 * pkin(5);
t213 = cos(qJ(3));
t208 = Icges(4,4) * t213;
t212 = sin(qJ(3));
t252 = -Icges(4,2) * t212 + t208;
t176 = Icges(4,1) * t212 + t208;
t210 = qJD(3) + qJD(4);
t359 = rSges(5,2) * t207;
t291 = t210 * t359;
t322 = t206 * t210;
t388 = rSges(5,1) * t322 + t291;
t173 = Icges(4,5) * t213 - Icges(4,6) * t212;
t172 = Icges(4,5) * t212 + Icges(4,6) * t213;
t234 = qJD(3) * t172;
t350 = Icges(4,4) * t212;
t177 = Icges(4,1) * t213 - t350;
t113 = Icges(4,5) * t204 + t177 * t205;
t111 = Icges(4,6) * t204 + t205 * t252;
t337 = t111 * t212;
t247 = -t113 * t213 + t337;
t342 = Icges(4,3) * t205;
t387 = -t205 * t234 + (-t173 * t204 + t247 + t342) * qJD(2);
t329 = t204 * t212;
t182 = Icges(4,4) * t329;
t328 = t204 * t213;
t348 = Icges(4,5) * t205;
t112 = Icges(4,1) * t328 - t182 - t348;
t344 = Icges(4,6) * t205;
t110 = Icges(4,4) * t328 - Icges(4,2) * t329 - t344;
t338 = t110 * t212;
t248 = -t112 * t213 + t338;
t109 = Icges(4,3) * t204 + t173 * t205;
t303 = qJD(2) * t109;
t386 = qJD(2) * t248 - t204 * t234 + t303;
t149 = Icges(5,5) * t206 + Icges(5,6) * t207;
t238 = t149 * t210;
t102 = Icges(5,6) * t204 + t205 * t251;
t339 = t102 * t206;
t341 = Icges(5,3) * t205;
t385 = -t205 * t238 + (-t104 * t207 - t150 * t204 + t339 + t341) * qJD(2);
t304 = qJD(2) * t100;
t384 = qJD(2) * t250 - t204 * t238 + t304;
t108 = Icges(4,5) * t328 - Icges(4,6) * t329 - t342;
t46 = -t205 * t108 - t204 * t248;
t151 = Icges(5,2) * t207 + t349;
t245 = t151 * t206 - t153 * t207;
t383 = qJD(2) * t245 + t150 * t210;
t174 = Icges(4,2) * t213 + t350;
t244 = t212 * t174 - t213 * t176;
t382 = t244 * qJD(2) + qJD(3) * t173;
t167 = rSges(5,2) * t331;
t105 = rSges(5,1) * t330 - rSges(5,3) * t205 - t167;
t195 = t204 * rSges(5,3);
t309 = rSges(5,1) * t325 + t195;
t326 = t205 * t206;
t106 = -rSges(5,2) * t326 + t309;
t155 = rSges(5,1) * t206 + t359;
t123 = t155 * t204;
t124 = t155 * t205;
t143 = t204 * t210;
t144 = t205 * t210;
t360 = rSges(5,2) * t206;
t361 = rSges(5,1) * t207;
t156 = -t360 + t361;
t200 = t205 * pkin(5);
t147 = pkin(2) * t204 - t200;
t214 = -pkin(6) - pkin(5);
t187 = t205 * t214;
t370 = pkin(3) * t213;
t203 = pkin(2) + t370;
t310 = -t203 * t204 - t187;
t95 = t147 + t310;
t266 = t203 * t205 - t204 * t214;
t96 = t266 - t148;
t35 = t105 * t143 + t106 * t144 + qJD(1) + (-t204 * t95 + t205 * t96) * qJD(3);
t296 = qJD(3) * t212;
t285 = t205 * t296;
t264 = pkin(3) * t285;
t233 = -t144 * t155 - t264;
t354 = -t105 + t95;
t38 = (-t147 + t354) * qJD(2) + t233;
t290 = pkin(3) * t296;
t170 = t204 * t290;
t353 = -t106 - t96;
t39 = -t143 * t155 - t170 + (t148 - t353) * qJD(2);
t263 = qJD(2) * t210;
t138 = t204 * t263;
t139 = t205 * t263;
t301 = qJD(2) * t204;
t300 = qJD(2) * t205;
t311 = rSges(5,3) * t300 + qJD(2) * t167;
t64 = -t205 * t291 + (-t205 * t322 - t207 * t301) * rSges(5,1) + t311;
t65 = -t210 * t123 + (t156 * t205 + t195) * qJD(2);
t188 = pkin(5) * t300;
t367 = pkin(2) - t203;
t81 = -t264 - t188 + (t204 * t367 - t187) * qJD(2);
t82 = -t170 + (-t367 * t205 + (-pkin(5) - t214) * t204) * qJD(2);
t8 = t105 * t139 - t106 * t138 + t143 * t65 + t144 * t64 + (t204 * t82 + t205 * t81 + (-t204 * t96 - t205 * t95) * qJD(2)) * qJD(3);
t381 = -t38 * (qJD(2) * t123 - t144 * t156) - t35 * (-t123 * t143 - t124 * t144) - t39 * (-qJD(2) * t124 - t143 * t156) + t8 * (t105 * t204 + t106 * t205);
t315 = -Icges(4,2) * t328 + t112 - t182;
t317 = t176 * t204 + t110;
t380 = -t212 * t315 - t213 * t317;
t379 = qJD(2) * t391 + t143 * (-t151 * t205 + t104) - t144 * (-Icges(5,2) * t330 + t103 - t164);
t378 = t138 / 0.2e1;
t377 = t139 / 0.2e1;
t376 = -t143 / 0.2e1;
t375 = t143 / 0.2e1;
t374 = -t144 / 0.2e1;
t373 = t144 / 0.2e1;
t372 = t204 / 0.2e1;
t371 = -t205 / 0.2e1;
t369 = -qJD(2) / 0.2e1;
t368 = qJD(2) / 0.2e1;
t99 = Icges(5,5) * t330 - Icges(5,6) * t331 - t341;
t366 = -t103 * t325 - t204 * t99;
t323 = t205 * t213;
t364 = -t108 * t204 - t112 * t323;
t363 = t109 * t204 + t113 * t323;
t362 = rSges(4,1) * t213;
t178 = rSges(4,1) * t212 + rSges(4,2) * t213;
t136 = t178 * t205;
t298 = qJD(3) * t204;
t287 = t178 * t298;
t196 = t204 * rSges(4,3);
t324 = t205 * t212;
t115 = rSges(4,1) * t323 - rSges(4,2) * t324 + t196;
t313 = t115 + t148;
t67 = qJD(2) * t313 - t287;
t358 = t136 * t67;
t305 = rSges(4,2) * t329 + rSges(4,3) * t205;
t114 = rSges(4,1) * t328 - t305;
t297 = qJD(3) * t205;
t286 = t178 * t297;
t66 = -t286 + (-t114 - t147) * qJD(2);
t357 = t204 * t66;
t356 = t205 * t39;
t355 = t205 * t66;
t336 = t149 * t204;
t335 = t149 * t205;
t334 = t151 * t210;
t333 = t172 * t204;
t332 = t172 * t205;
t56 = -t204 * t245 - t335;
t321 = t56 * qJD(2);
t77 = -t204 * t244 - t332;
t320 = t77 * qJD(2);
t316 = -t176 * t205 - t111;
t314 = -t174 * t205 + t113;
t299 = qJD(2) * t212;
t308 = rSges(4,2) * t204 * t299 + rSges(4,3) * t300;
t307 = -t174 + t177;
t306 = t176 + t252;
t302 = qJD(2) * t173;
t295 = qJD(3) * t213;
t294 = t105 * t300 + t204 * t65 + t205 * t64;
t293 = (qJD(3) ^ 2) * t370;
t284 = -pkin(2) - t362;
t283 = t301 / 0.2e1;
t282 = t300 / 0.2e1;
t281 = -t298 / 0.2e1;
t278 = t297 / 0.2e1;
t277 = -pkin(3) * t212 - t155;
t240 = t153 * t210;
t275 = qJD(2) * t104 - t101 * t210 - t204 * t240;
t274 = -t99 + t339;
t273 = -t102 * t210 - t205 * t240 + (-t154 * t204 + t347) * qJD(2);
t272 = qJD(2) * t102 + t103 * t210 - t204 * t334;
t271 = t104 * t210 - t205 * t334 + (-t204 * t251 + t343) * qJD(2);
t92 = t113 * t328;
t270 = t205 * t109 - t92;
t269 = -t108 + t337;
t268 = t391 * t210;
t267 = t154 * t210 - t334;
t128 = t156 * t210;
t260 = -pkin(3) * t295 - t128;
t83 = t104 * t330;
t259 = t102 * t331 - t83;
t257 = -rSges(4,2) * t212 + t362;
t256 = -t204 * t39 - t205 * t38;
t255 = -t204 * t67 - t355;
t50 = t101 * t207 + t103 * t206;
t68 = t110 * t213 + t112 * t212;
t69 = t111 * t213 + t113 * t212;
t246 = t114 * t204 + t115 * t205;
t135 = t178 * t204;
t47 = -t111 * t329 - t270;
t242 = (t204 * t47 - t205 * t46) * qJD(3);
t48 = -t110 * t324 - t364;
t49 = -t111 * t324 + t363;
t241 = (t204 * t49 - t205 * t48) * qJD(3);
t236 = qJD(3) * t176;
t235 = qJD(3) * t174;
t232 = qJD(2) * t150 - t143 * t335 + t144 * t336;
t231 = -t212 * t314 + t213 * t316;
t222 = -t206 * t271 + t207 * t273 + t304;
t10 = t204 * t385 + t205 * t222;
t223 = qJD(2) * t99 - t206 * t272 + t207 * t275;
t11 = t204 * t223 - t205 * t384;
t12 = t204 * t222 - t205 * t385;
t40 = -t205 * t99 - t237;
t41 = -t100 * t205 - t259;
t20 = t143 * t41 - t144 * t40 + t321;
t42 = -t101 * t326 - t366;
t43 = -t102 * t326 + t365;
t57 = -t205 * t245 + t336;
t52 = t57 * qJD(2);
t21 = t143 * t43 - t144 * t42 + t52;
t226 = (-t153 * t205 - t102) * t143 - (-t153 * t204 - t101) * t144 + (-t151 + t154) * qJD(2);
t217 = -t206 * t379 + t207 * t226;
t28 = t206 * t275 + t207 * t272;
t29 = t206 * t273 + t207 * t271;
t221 = qJD(2) * t149 - t206 * t268 + t207 * t267;
t30 = t204 * t383 + t205 * t221;
t31 = t204 * t221 - t205 * t383;
t51 = t102 * t207 + t104 * t206;
t9 = t204 * t384 + t205 * t223;
t230 = (qJD(2) * t30 + t10 * t143 + t138 * t42 + t139 * t43 - t144 * t9) * t372 + (t206 * t226 + t207 * t379) * t369 + t20 * t283 + t21 * t282 + (qJD(2) * t31 - t11 * t144 + t12 * t143 + t138 * t40 + t139 * t41) * t371 + (t204 * t41 - t205 * t40) * t378 + (t204 * t43 - t205 * t42) * t377 + (t10 * t204 - t205 * t9 + (t204 * t42 + t205 * t43) * qJD(2)) * t375 + (-t11 * t205 + t12 * t204 + (t204 * t40 + t205 * t41) * qJD(2)) * t374 + (t204 * t29 - t205 * t28 + (t204 * t50 + t205 * t51) * qJD(2)) * t368 + (t204 * t232 + t205 * t217) * t376 + (t204 * t217 - t205 * t232) * t373;
t229 = (-t212 * t306 + t213 * t307) * qJD(2);
t79 = -rSges(4,2) * t205 * t295 + (-t213 * t301 - t285) * rSges(4,1) + t308;
t80 = -qJD(3) * t135 + (t205 * t257 + t196) * qJD(2);
t227 = t204 * t80 + t205 * t79 + (t114 * t205 - t115 * t204) * qJD(2);
t74 = qJD(2) * t111 - t204 * t235;
t76 = qJD(2) * t113 - t204 * t236;
t220 = qJD(2) * t108 - qJD(3) * t68 - t212 * t74 + t213 * t76;
t73 = -t205 * t235 + (-t204 * t252 + t344) * qJD(2);
t75 = -t205 * t236 + (-t177 * t204 + t348) * qJD(2);
t219 = -qJD(3) * t69 - t212 * t73 + t213 * t75 + t303;
t158 = t252 * qJD(3);
t159 = t177 * qJD(3);
t218 = qJD(2) * t172 - t158 * t212 + t159 * t213 + (-t174 * t213 - t176 * t212) * qJD(3);
t161 = t257 * qJD(3);
t142 = qJD(2) * t147;
t141 = t148 * qJD(2);
t137 = qJD(2) * (-pkin(2) * t301 + t188);
t78 = -t205 * t244 + t333;
t70 = t78 * qJD(2);
t55 = qJD(3) * t246 + qJD(1);
t45 = -t161 * t297 + (-t141 - t80 + t287) * qJD(2);
t44 = -t161 * t298 + t137 + (t79 - t286) * qJD(2);
t37 = t204 * t218 - t205 * t382;
t36 = t204 * t382 + t205 * t218;
t34 = -qJD(3) * t247 + t212 * t75 + t213 * t73;
t33 = -qJD(3) * t248 + t212 * t76 + t213 * t74;
t32 = t227 * qJD(3);
t27 = -t205 * t293 - t128 * t144 + t138 * t155 + (-t141 - t65 - t82 + t170) * qJD(2);
t26 = -t204 * t293 - t128 * t143 - t139 * t155 + t137 + (t64 + t81 - t264) * qJD(2);
t23 = t70 + t241;
t22 = t242 + t320;
t1 = [m(4) * t32 + m(5) * t8; (t70 + ((t47 - t92 + (t109 + t338) * t205 + t364) * t205 + t363 * t204) * qJD(3)) * t278 + (t52 + (t41 + (t100 + t340) * t205 + t259 + t366) * t144 + (-t205 * t274 - t390 + t40) * t143) * t373 + (t50 + t56) * t378 + (t51 + t57) * t377 + (-t321 + (t43 + t390) * t144 + (t274 * t204 + t42 - t83) * t143 + ((t100 + t250) * t143 + t274 * t144) * t205 + t20) * t376 + (t29 + t30) * t375 + (t22 - t320 + ((t205 * t269 - t363 + t49) * t205 + (t204 * t269 + t270 + t48) * t204) * qJD(3)) * t281 + (t34 + t36) * t298 / 0.2e1 + (-qJD(3) * t244 + t158 * t213 + t159 * t212 + t206 * t267 + t207 * t268) * qJD(2) + (-(qJD(2) * t354 - t142 + t233 - t38) * t39 + t27 * (-t105 + t310) + t38 * (t388 * t204 + t170) + t26 * (t266 + t309) + t39 * t311 + (-t26 * t360 + t39 * (-t290 - t388)) * t205 + ((t38 * (-t156 - t203) - t39 * t214) * t205 + (t38 * (-rSges(5,3) + t214) + t39 * (-t203 - t361)) * t204) * qJD(2)) * m(5) + (-(-qJD(2) * t114 - t142 - t286 - t66) * t67 + t45 * (t204 * t284 + t200 + t305) + t44 * t313 + t67 * (t188 + t308) + (t178 * t357 - t358) * qJD(3) + ((-pkin(2) - t257) * t355 + (t66 * (-rSges(4,3) - pkin(5)) + t67 * t284) * t204) * qJD(2)) * m(4) + (t28 + t31 + t21) * t374 - (t33 + t37 + t23) * t297 / 0.2e1 + ((t68 + t77) * t204 + (t69 + t78) * t205) * qJD(3) * t368; t230 + ((t212 * t307 + t213 * t306) * qJD(2) + ((t204 * t314 - t205 * t315) * t213 + (t204 * t316 + t205 * t317) * t212) * qJD(3)) * t369 + ((-t298 * t332 + t302) * t204 + (t229 + (-t380 * t205 + (t333 + t231) * t204) * qJD(3)) * t205) * t281 + ((-t297 * t333 - t302) * t205 + (t229 + (t231 * t204 + (t332 - t380) * t205) * qJD(3)) * t204) * t278 + (t34 * t204 - t33 * t205 + (t204 * t68 + t205 * t69) * qJD(2)) * t368 + (qJD(2) * t36 + (-(t204 * t386 + t205 * t220) * t205 + (t204 * t387 + t205 * t219) * t204 + (t48 * t204 + t49 * t205) * qJD(2)) * t389) * t372 + (qJD(2) * t37 + (-(t204 * t220 - t205 * t386) * t205 + (t204 * t219 - t205 * t387) * t204 + (t46 * t204 + t47 * t205) * qJD(2)) * t389) * t371 + (t242 + t22) * t283 + (t241 + t23) * t282 + (t35 * t294 + (t27 * t277 + t38 * t260 + t8 * t96 + t35 * t81 + (t277 * t39 - t35 * t95) * qJD(2)) * t205 + (t26 * t277 + t39 * t260 - t8 * t95 + t35 * t82 + (t155 * t38 + t35 * t353) * qJD(2)) * t204 - (-t299 * t356 + (t256 * t213 + t35 * (-t204 ^ 2 - t205 ^ 2) * t212) * qJD(3)) * pkin(3) + t381) * m(5) + (-(t135 * t66 - t358) * qJD(2) - (t55 * (-t135 * t204 - t136 * t205) + t255 * t257) * qJD(3) + t32 * t246 + t55 * t227 + t255 * t161 + (-t44 * t204 - t45 * t205 + (-t205 * t67 + t357) * qJD(2)) * t178) * m(4); t230 + (t35 * (-t106 * t301 + t294) + t256 * t128 + (-t26 * t204 - t27 * t205 + (t204 * t38 - t356) * qJD(2)) * t155 + t381) * m(5);];
tauc = t1(:);
