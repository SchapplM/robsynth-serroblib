% Calculate time derivative of joint inertia matrix for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:12
% EndTime: 2019-03-09 02:00:28
% DurationCPUTime: 9.65s
% Computational Cost: add. (28223->686), mult. (27558->973), div. (0->0), fcn. (26516->10), ass. (0->337)
t241 = pkin(10) + qJ(4);
t236 = sin(t241);
t238 = cos(t241);
t246 = sin(qJ(5));
t248 = cos(qJ(5));
t285 = Icges(6,5) * t248 - Icges(6,6) * t246;
t183 = -Icges(6,3) * t238 + t236 * t285;
t287 = Icges(7,4) * t248 + Icges(7,6) * t246;
t184 = -Icges(7,2) * t238 + t236 * t287;
t439 = t183 + t184;
t242 = qJ(1) + pkin(9);
t239 = cos(t242);
t360 = t239 * t248;
t237 = sin(t242);
t364 = t237 * t246;
t197 = t238 * t364 + t360;
t361 = t239 * t246;
t363 = t237 * t248;
t198 = t238 * t363 - t361;
t432 = rSges(7,3) + qJ(6);
t433 = rSges(7,1) + pkin(5);
t438 = -t432 * t197 - t198 * t433;
t342 = qJD(4) * t246;
t322 = t236 * t342;
t338 = qJD(5) * t248;
t339 = qJD(5) * t246;
t347 = qJD(1) * t239;
t349 = qJD(1) * t237;
t117 = -t237 * t322 - t239 * t339 - t248 * t349 + (t237 * t338 + t246 * t347) * t238;
t341 = qJD(4) * t248;
t253 = -t236 * t341 + (-qJD(5) * t238 + qJD(1)) * t246;
t348 = qJD(1) * t238;
t310 = -qJD(5) + t348;
t118 = t237 * t253 + t310 * t360;
t437 = -t197 * qJD(6) - t432 * t117 - t118 * t433;
t368 = t236 * t237;
t130 = Icges(6,5) * t198 - Icges(6,6) * t197 + Icges(6,3) * t368;
t134 = Icges(6,4) * t198 - Icges(6,2) * t197 + Icges(6,6) * t368;
t138 = Icges(6,1) * t198 - Icges(6,4) * t197 + Icges(6,5) * t368;
t199 = t238 * t361 - t363;
t333 = t238 * t360;
t200 = t333 + t364;
t367 = t236 * t239;
t55 = t130 * t367 - t134 * t199 + t138 * t200;
t131 = Icges(6,5) * t200 - Icges(6,6) * t199 + Icges(6,3) * t367;
t135 = Icges(6,4) * t200 - Icges(6,2) * t199 + Icges(6,6) * t367;
t139 = Icges(6,1) * t200 - Icges(6,4) * t199 + Icges(6,5) * t367;
t56 = t131 * t367 - t135 * t199 + t139 * t200;
t295 = t237 * t55 + t239 * t56;
t128 = Icges(7,5) * t198 + Icges(7,6) * t368 + Icges(7,3) * t197;
t132 = Icges(7,4) * t198 + Icges(7,2) * t368 + Icges(7,6) * t197;
t136 = Icges(7,1) * t198 + Icges(7,4) * t368 + Icges(7,5) * t197;
t53 = t128 * t199 + t132 * t367 + t136 * t200;
t129 = Icges(7,5) * t200 + Icges(7,6) * t367 + Icges(7,3) * t199;
t133 = Icges(7,4) * t200 + Icges(7,2) * t367 + Icges(7,6) * t199;
t137 = Icges(7,1) * t200 + Icges(7,4) * t367 + Icges(7,5) * t199;
t54 = t129 * t199 + t133 * t367 + t137 * t200;
t296 = t237 * t53 + t239 * t54;
t436 = t295 + t296;
t51 = t130 * t368 - t134 * t197 + t138 * t198;
t52 = t131 * t368 - t135 * t197 + t139 * t198;
t297 = t237 * t51 + t239 * t52;
t49 = t128 * t197 + t132 * t368 + t136 * t198;
t50 = t129 * t197 + t133 * t368 + t137 * t198;
t298 = t237 * t49 + t239 * t50;
t435 = t297 + t298;
t434 = -qJD(1) * t236 / 0.2e1;
t382 = Icges(7,5) * t248;
t284 = Icges(7,3) * t246 + t382;
t182 = -Icges(7,6) * t238 + t236 * t284;
t385 = Icges(6,4) * t248;
t288 = -Icges(6,2) * t246 + t385;
t185 = -Icges(6,6) * t238 + t236 * t288;
t383 = Icges(7,5) * t246;
t291 = Icges(7,1) * t248 + t383;
t186 = -Icges(7,4) * t238 + t236 * t291;
t386 = Icges(6,4) * t246;
t292 = Icges(6,1) * t248 - t386;
t187 = -Icges(6,5) * t238 + t236 * t292;
t431 = t439 * t238 + ((-t186 - t187) * t248 + (-t182 + t185) * t246) * t236;
t319 = t238 * t342;
t255 = t236 * t338 + t319;
t280 = t129 * t246 + t137 * t248;
t61 = -t133 * t238 + t236 * t280;
t278 = -t135 * t246 + t139 * t248;
t63 = -t131 * t238 + t236 * t278;
t395 = t61 + t63;
t281 = t128 * t246 + t136 * t248;
t60 = -t132 * t238 + t236 * t281;
t279 = -t134 * t246 + t138 * t248;
t62 = -t130 * t238 + t236 * t279;
t396 = t60 + t62;
t430 = t237 * t396 + t239 * t395;
t115 = qJD(1) * t197 - qJD(5) * t333 - t237 * t339 + t239 * t322;
t116 = t239 * t253 - t310 * t363;
t343 = qJD(4) * t239;
t320 = t238 * t343;
t324 = t236 * t349;
t256 = t320 - t324;
t344 = qJD(4) * t238;
t321 = t237 * t344;
t257 = t236 * t347 + t321;
t68 = Icges(7,5) * t118 + Icges(7,6) * t257 + Icges(7,3) * t117;
t72 = Icges(7,4) * t118 + Icges(7,2) * t257 + Icges(7,6) * t117;
t76 = Icges(7,1) * t118 + Icges(7,4) * t257 + Icges(7,5) * t117;
t11 = -t115 * t128 + t116 * t136 + t132 * t256 + t199 * t68 + t200 * t76 + t367 * t72;
t67 = Icges(7,5) * t116 + Icges(7,6) * t256 - Icges(7,3) * t115;
t71 = Icges(7,4) * t116 + Icges(7,2) * t256 - Icges(7,6) * t115;
t75 = Icges(7,1) * t116 + Icges(7,4) * t256 - Icges(7,5) * t115;
t12 = -t115 * t129 + t116 * t137 + t133 * t256 + t199 * t67 + t200 * t75 + t367 * t71;
t70 = Icges(6,5) * t118 - Icges(6,6) * t117 + Icges(6,3) * t257;
t74 = Icges(6,4) * t118 - Icges(6,2) * t117 + Icges(6,6) * t257;
t78 = Icges(6,1) * t118 - Icges(6,4) * t117 + Icges(6,5) * t257;
t13 = t115 * t134 + t116 * t138 + t130 * t256 - t199 * t74 + t200 * t78 + t367 * t70;
t69 = Icges(6,5) * t116 + Icges(6,6) * t115 + Icges(6,3) * t256;
t73 = Icges(6,4) * t116 + Icges(6,2) * t115 + Icges(6,6) * t256;
t77 = Icges(6,1) * t116 + Icges(6,4) * t115 + Icges(6,5) * t256;
t14 = t115 * t135 + t116 * t139 + t131 * t256 - t199 * t73 + t200 * t77 + t367 * t69;
t429 = (-t11 - t13) * t239 + (t12 + t14) * t237 + t436 * qJD(1);
t15 = t117 * t128 + t118 * t136 + t132 * t257 + t197 * t68 + t198 * t76 + t368 * t72;
t16 = t117 * t129 + t118 * t137 + t133 * t257 + t197 * t67 + t198 * t75 + t368 * t71;
t17 = -t117 * t134 + t118 * t138 + t130 * t257 - t197 * t74 + t198 * t78 + t368 * t70;
t18 = -t117 * t135 + t118 * t139 + t131 * t257 - t197 * t73 + t198 * t77 + t368 * t69;
t428 = (-t15 - t17) * t239 + (t16 + t18) * t237 + t435 * qJD(1);
t19 = (qJD(4) * t281 - t72) * t238 + (qJD(4) * t132 + t246 * t68 + t248 * t76 + (t128 * t248 - t136 * t246) * qJD(5)) * t236;
t21 = (qJD(4) * t279 - t70) * t238 + (qJD(4) * t130 - t246 * t74 + t248 * t78 + (-t134 * t248 - t138 * t246) * qJD(5)) * t236;
t427 = -t19 - t21;
t20 = (qJD(4) * t280 - t71) * t238 + (qJD(4) * t133 + t246 * t67 + t248 * t75 + (t129 * t248 - t137 * t246) * qJD(5)) * t236;
t22 = (qJD(4) * t278 - t69) * t238 + (qJD(4) * t131 - t246 * t73 + t248 * t77 + (-t135 * t248 - t139 * t246) * qJD(5)) * t236;
t426 = t20 + t22;
t88 = t182 * t197 + t184 * t368 + t186 * t198;
t89 = t183 * t368 - t185 * t197 + t187 * t198;
t425 = (-t88 - t89) * t238 + t435 * t236;
t90 = t182 * t199 + t184 * t367 + t186 * t200;
t91 = t183 * t367 - t185 * t199 + t187 * t200;
t424 = (-t90 - t91) * t238 + t436 * t236;
t415 = 2 * m(5);
t391 = rSges(5,1) * t238;
t305 = -rSges(5,2) * t236 + t391;
t177 = -rSges(5,3) * t239 + t237 * t305;
t362 = t238 * t239;
t234 = t237 * rSges(5,3);
t419 = -rSges(5,2) * t367 + t234;
t178 = rSges(5,1) * t362 + t419;
t213 = rSges(5,1) * t236 + rSges(5,2) * t238;
t265 = qJD(4) * t213;
t252 = rSges(5,2) * t324 + rSges(5,3) * t347 - t239 * t265;
t59 = (qJD(1) * t177 + t252) * t239 + (-t237 * t265 + (-t178 + t419) * qJD(1)) * t237;
t423 = t415 * t59;
t340 = qJD(5) * t236;
t151 = (Icges(7,3) * t248 - t383) * t340 + (Icges(7,6) * t236 + t238 * t284) * qJD(4);
t153 = (-Icges(7,4) * t246 + Icges(7,6) * t248) * t340 + (Icges(7,2) * t236 + t238 * t287) * qJD(4);
t155 = (-Icges(7,1) * t246 + t382) * t340 + (Icges(7,4) * t236 + t238 * t291) * qJD(4);
t156 = (-Icges(6,1) * t246 - t385) * t340 + (Icges(6,5) * t236 + t238 * t292) * qJD(4);
t317 = t236 * t339;
t318 = t238 * t341;
t346 = qJD(4) * t236;
t366 = t236 * t246;
t422 = t187 * t318 + t151 * t366 - t238 * t153 + (-t317 + t318) * t186 + t255 * t182 + (t156 + t155) * t236 * t248 + t439 * t346;
t421 = rSges(7,2) * t320 + qJD(6) * t199 - t432 * t115 + t116 * t433;
t388 = Icges(5,4) * t236;
t294 = Icges(5,1) * t238 - t388;
t176 = Icges(5,5) * t237 + t239 * t294;
t369 = t176 * t238;
t387 = Icges(5,4) * t238;
t290 = -Icges(5,2) * t236 + t387;
t174 = Icges(5,6) * t237 + t239 * t290;
t374 = t174 * t236;
t274 = -t369 + t374;
t420 = t274 * t239;
t245 = -pkin(7) - qJ(3);
t244 = cos(pkin(10));
t235 = pkin(3) * t244 + pkin(2);
t267 = -t235 - t305;
t403 = sin(qJ(1)) * pkin(1);
t148 = -t403 + (rSges(5,3) - t245) * t239 + t267 * t237;
t240 = cos(qJ(1)) * pkin(1);
t308 = t239 * t235 - t237 * t245 + t240;
t149 = t178 + t308;
t418 = t148 * t239 + t149 * t237;
t286 = Icges(5,5) * t238 - Icges(5,6) * t236;
t171 = -Icges(5,3) * t239 + t237 * t286;
t173 = -Icges(5,6) * t239 + t237 * t290;
t175 = -Icges(5,5) * t239 + t237 * t294;
t417 = t238 * t395 - t424;
t357 = rSges(7,2) * t367 + t432 * t199 + t200 * t433;
t358 = rSges(7,2) * t368 - t438;
t416 = -t237 * t358 - t239 * t357;
t272 = rSges(4,1) * t244 - rSges(4,2) * sin(pkin(10)) + pkin(2);
t390 = rSges(4,3) + qJ(3);
t164 = t237 * t390 + t239 * t272 + t240;
t414 = 2 * m(6);
t413 = 2 * m(7);
t412 = t237 ^ 2;
t411 = t239 ^ 2;
t409 = -t238 / 0.2e1;
t406 = -rSges(7,2) - pkin(8);
t405 = -rSges(6,3) - pkin(8);
t404 = m(5) * t213;
t402 = pkin(4) * t236;
t401 = pkin(4) * t238;
t152 = (-Icges(6,5) * t246 - Icges(6,6) * t248) * t340 + (Icges(6,3) * t236 + t238 * t285) * qJD(4);
t154 = (-Icges(6,2) * t248 - t386) * t340 + (Icges(6,6) * t236 + t238 * t288) * qJD(4);
t359 = t246 * t154;
t397 = (-t185 * t342 - t152) * t238 + (-t359 + (-t185 * t248 - t187 * t246) * qJD(5)) * t236 + t422;
t394 = -rSges(7,2) * t324 + t421;
t393 = rSges(7,2) * t257 - t437;
t392 = t431 * t346;
t303 = -rSges(6,1) * t198 + rSges(6,2) * t197;
t141 = rSges(6,3) * t368 - t303;
t379 = t141 * t239;
t376 = t173 * t236;
t375 = t173 * t238;
t373 = t174 * t238;
t372 = t175 * t236;
t371 = t175 * t238;
t370 = t176 * t236;
t299 = pkin(5) * t248 + qJ(6) * t246;
t301 = rSges(7,1) * t248 + rSges(7,3) * t246;
t356 = t299 * t344 + (qJD(6) * t246 + (-pkin(5) * t246 + qJ(6) * t248) * qJD(5)) * t236 + (-rSges(7,1) * t246 + rSges(7,3) * t248) * t340 + (rSges(7,2) * t236 + t238 * t301) * qJD(4);
t302 = rSges(6,1) * t248 - rSges(6,2) * t246;
t158 = (-rSges(6,1) * t246 - rSges(6,2) * t248) * t340 + (rSges(6,3) * t236 + t238 * t302) * qJD(4);
t309 = pkin(8) * t236 + t401;
t207 = t309 * qJD(4);
t355 = -t158 - t207;
t195 = t309 * t237;
t225 = pkin(4) * t362;
t196 = pkin(8) * t367 + t225;
t354 = t237 * t195 + t239 * t196;
t353 = -rSges(7,2) * t238 + (t299 + t301) * t236;
t189 = -rSges(6,3) * t238 + t236 * t302;
t214 = -pkin(8) * t238 + t402;
t352 = -t189 - t214;
t233 = qJD(3) * t239;
t351 = t245 * t349 + t233;
t172 = Icges(5,3) * t237 + t239 * t286;
t350 = qJD(1) * t172;
t345 = qJD(4) * t237;
t31 = t237 * t50 - t239 * t49;
t32 = t237 * t52 - t239 * t51;
t335 = t31 / 0.2e1 + t32 / 0.2e1;
t33 = t237 * t54 - t239 * t53;
t34 = t237 * t56 - t239 * t55;
t334 = t33 / 0.2e1 + t34 / 0.2e1;
t332 = t116 * rSges(6,1) + t115 * rSges(6,2) + rSges(6,3) * t320;
t218 = t345 * t402;
t219 = pkin(8) * t320;
t323 = t236 * t343;
t329 = t237 * (pkin(8) * t257 + qJD(1) * t225 - t218) + t239 * (-pkin(8) * t324 + t219 + (-t237 * t348 - t323) * pkin(4)) + t195 * t347;
t328 = -t207 - t356;
t327 = -t214 - t353;
t143 = t200 * rSges(6,1) - t199 * rSges(6,2) + rSges(6,3) * t367;
t326 = t218 + t351;
t325 = t189 * t349;
t315 = -t235 - t401;
t314 = t237 * t353;
t313 = t239 * t353;
t312 = t358 * t239;
t160 = t352 * t239;
t311 = qJD(1) * t353;
t109 = t327 * t239;
t307 = t237 * t311;
t304 = t118 * rSges(6,1) - t117 * rSges(6,2);
t300 = -t239 * t245 - t403;
t289 = Icges(5,2) * t238 + t388;
t277 = -t143 * t237 + t379;
t276 = -t141 * t237 - t143 * t239;
t275 = -t371 + t376;
t273 = -t238 * t396 + t425;
t37 = t117 * t182 + t118 * t186 + t151 * t197 + t153 * t368 + t155 * t198 + t184 * t257;
t38 = -t117 * t185 + t118 * t187 + t152 * t368 - t154 * t197 + t156 * t198 + t183 * t257;
t271 = t21 / 0.2e1 + t19 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1;
t35 = -t115 * t182 + t116 * t186 + t151 * t199 + t153 * t367 + t155 * t200 + t184 * t256;
t36 = t115 * t185 + t116 * t187 + t152 * t367 - t154 * t199 + t156 * t200 + t183 * t256;
t270 = t35 / 0.2e1 + t36 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1;
t269 = t63 / 0.2e1 + t61 / 0.2e1 + t90 / 0.2e1 + t91 / 0.2e1;
t268 = t88 / 0.2e1 + t89 / 0.2e1 + t62 / 0.2e1 + t60 / 0.2e1;
t232 = qJD(3) * t237;
t266 = -pkin(4) * t323 + t219 + t232;
t264 = t308 + t196;
t262 = qJD(4) * t289;
t261 = qJD(4) * (-Icges(5,5) * t236 - Icges(5,6) * t238);
t260 = t236 * t406 + t315;
t259 = t236 * t405 + t315;
t254 = -t237 * t357 + t312;
t251 = t237 * t260 + t300;
t250 = t237 * t259 + t300;
t163 = -t237 * t272 + t239 * t390 - t403;
t206 = t305 * qJD(4);
t201 = t214 * t349;
t159 = t352 * t237;
t145 = -qJD(1) * t164 + t233;
t144 = qJD(1) * t163 + t232;
t120 = t237 * t261 + t350;
t119 = -qJD(1) * t171 + t239 * t261;
t108 = t327 * t237;
t107 = t213 * t345 + (t239 * t267 - t234 - t240) * qJD(1) + t351;
t106 = t232 + ((-t235 - t391) * t237 + t300) * qJD(1) + t252;
t103 = -t143 * t238 - t189 * t367;
t102 = t141 * t238 + t189 * t368;
t101 = t264 + t143;
t100 = t250 + t303;
t97 = t172 * t237 - t420;
t96 = t171 * t237 - t239 * t275;
t95 = -t172 * t239 - t274 * t237;
t94 = -t171 * t239 - t237 * t275;
t93 = qJD(1) * t160 + t237 * t355;
t92 = t239 * t355 + t201 + t325;
t87 = t277 * t236;
t86 = t264 + t357;
t85 = t251 + t438;
t82 = rSges(6,3) * t257 + t304;
t80 = -rSges(6,3) * t324 + t332;
t66 = -t236 * t313 - t238 * t357;
t65 = t236 * t314 + t238 * t358;
t64 = -t276 + t354;
t58 = qJD(1) * t109 + t237 * t328;
t57 = t239 * t328 + t201 + t307;
t48 = t405 * t321 + (t239 * t259 - t240) * qJD(1) - t304 + t326;
t47 = qJD(1) * t250 + t266 + t332;
t46 = t254 * t236;
t45 = t354 - t416;
t42 = (t189 * t345 + t82) * t238 + (-qJD(4) * t141 + t158 * t237 + t189 * t347) * t236;
t41 = (-t189 * t343 - t80) * t238 + (qJD(4) * t143 - t158 * t239 + t325) * t236;
t40 = t406 * t321 + (t239 * t260 - t240) * qJD(1) + t326 + t437;
t39 = qJD(1) * t251 + t266 + t421;
t30 = t277 * t344 + (qJD(1) * t276 - t237 * t80 + t239 * t82) * t236;
t29 = t237 * t82 + t239 * t80 + (t379 + (-t143 - t196) * t237) * qJD(1) + t329;
t24 = (qJD(4) * t314 + t393) * t238 + (-qJD(4) * t358 + t237 * t356 + t239 * t311) * t236;
t23 = (-qJD(4) * t313 - t394) * t238 + (qJD(4) * t357 - t239 * t356 + t307) * t236;
t10 = t394 * t239 + t393 * t237 + (t312 + (-t196 - t357) * t237) * qJD(1) + t329;
t9 = t254 * t344 + (qJD(1) * t416 - t394 * t237 + t393 * t239) * t236;
t4 = (qJD(4) * t297 - t38) * t238 + (-qJD(1) * t32 + qJD(4) * t89 + t17 * t237 + t18 * t239) * t236;
t3 = (qJD(4) * t298 - t37) * t238 + (-qJD(1) * t31 + qJD(4) * t88 + t15 * t237 + t16 * t239) * t236;
t2 = (qJD(4) * t295 - t36) * t238 + (-qJD(1) * t34 + qJD(4) * t91 + t13 * t237 + t14 * t239) * t236;
t1 = (qJD(4) * t296 - t35) * t238 + (-qJD(1) * t33 + qJD(4) * t90 + t11 * t237 + t12 * t239) * t236;
t5 = [-t236 * t359 - t187 * t317 - t238 * t152 + (t106 * t149 + t107 * t148) * t415 + 0.2e1 * m(4) * (t144 * t164 + t145 * t163) + (t100 * t48 + t101 * t47) * t414 + (t39 * t86 + t40 * t85) * t413 + t422 + (-t289 + t294) * t346 + (Icges(5,1) * t236 + t290 + t387) * t344 - t255 * t185; 0; 0; m(7) * (t237 * t40 - t239 * t39 + (t237 * t86 + t239 * t85) * qJD(1)) + m(6) * (t237 * t48 - t239 * t47 + (t100 * t239 + t101 * t237) * qJD(1)) + m(5) * (qJD(1) * t418 - t106 * t239 + t107 * t237) + m(4) * (-t144 * t239 + t145 * t237 + (t163 * t239 + t164 * t237) * qJD(1)); 0; 0; ((qJD(1) * t174 - t237 * t262) * t409 + t176 * t434 + (t376 / 0.2e1 - t371 / 0.2e1) * qJD(4) - t271) * t239 + ((-qJD(1) * t173 - t239 * t262) * t238 / 0.2e1 + t175 * t434 + (-t374 / 0.2e1 + t369 / 0.2e1) * qJD(4) + t270) * t237 + m(5) * ((-t106 * t237 - t107 * t239) * t213 - t418 * t206) + m(6) * (t100 * t92 + t101 * t93 + t159 * t47 + t160 * t48) + m(7) * (t108 * t39 + t109 * t40 + t57 * t85 + t58 * t86) + (t411 / 0.2e1 + t412 / 0.2e1) * t286 * qJD(4) + ((t373 / 0.2e1 + t370 / 0.2e1 - t149 * t404 + t269) * t239 + (t148 * t404 + t375 / 0.2e1 + t372 / 0.2e1 + t268) * t237) * qJD(1); m(5) * t59 + m(6) * t29 + m(7) * t10; m(6) * (t237 * t92 - t239 * t93 + (t159 * t237 + t160 * t239) * qJD(1)) + m(7) * (t237 * t57 - t239 * t58 + (t108 * t237 + t109 * t239) * qJD(1)); (t45 * t10 + t108 * t58 + t109 * t57) * t413 + (t159 * t93 + t160 * t92 + t64 * t29) * t414 + (t411 + t412) * t213 * t206 * t415 + (t178 * t423 + (-t95 * qJD(1) + (-qJD(1) * t275 - t120) * t239) * t239 - t428) * t239 + (t177 * t423 + (t96 * qJD(1) + (t274 * qJD(1) + t119) * t237) * t237 + ((-t120 + (-t370 - t373) * qJD(4) + t174 * t344 + t176 * t346 - t350) * t237 + (t173 * t344 + t175 * t346 + t119 - (t372 + t375) * qJD(4)) * t239 + (t97 - t94 + (t172 - t275) * t237 + t420) * qJD(1)) * t239 + t429) * t237 + (t237 * t95 - t239 * t94 + t31 + t32) * t349 + (t237 * t97 - t239 * t96 + t33 + t34) * t347; m(6) * (t100 * t42 + t101 * t41 + t102 * t48 + t103 * t47) + m(7) * (t23 * t86 + t24 * t85 + t66 * t39 + t65 * t40) + ((t237 * t268 + t239 * t269) * qJD(4) - t397) * t238 + (t270 * t239 + t271 * t237 + (-t237 * t269 + t239 * t268) * qJD(1)) * t236 - t392; m(6) * t30 + m(7) * t9; m(6) * (t237 * t42 - t239 * t41 + (t102 * t239 + t103 * t237) * qJD(1)) + m(7) * (-t23 * t239 + t237 * t24 + (t237 * t66 + t239 * t65) * qJD(1)); m(6) * (t102 * t92 + t103 * t93 + t159 * t41 + t160 * t42 + t29 * t87 + t30 * t64) + m(7) * (t46 * t10 + t108 * t23 + t109 * t24 + t9 * t45 + t65 * t57 + t66 * t58) + (-t3 / 0.2e1 - t4 / 0.2e1 + t334 * t344) * t239 + (t1 / 0.2e1 + t2 / 0.2e1 + t335 * t344) * t237 + ((-t237 * t334 + t239 * t335) * qJD(1) + t428 * t237 / 0.2e1 + t429 * t239 / 0.2e1 + (t237 * t395 - t239 * t396) * qJD(4) / 0.2e1) * t236 + (qJD(1) * t430 + t426 * t237 + t427 * t239) * t409 + (t425 * t237 + t424 * t239) * qJD(1) / 0.2e1; (t23 * t66 + t24 * t65 + t46 * t9) * t413 + (t102 * t42 + t103 * t41 + t30 * t87) * t414 + (t397 * t238 + (t273 * t237 - t239 * t417) * qJD(4) + t392) * t238 + ((-t238 * t426 + t1 + t2) * t239 + (t238 * t427 + t3 + t4) * t237 + (t236 * t430 + t238 * t431) * qJD(4) + (t237 * t417 + t273 * t239) * qJD(1)) * t236; m(7) * (-t115 * t85 + t117 * t86 + t197 * t39 + t199 * t40); t255 * m(7); m(7) * (-t115 * t237 - t117 * t239 + (t197 * t237 + t199 * t239) * qJD(1)); m(7) * (t45 * t319 + t108 * t117 - t109 * t115 + t197 * t58 + t199 * t57 + (t10 * t246 + t338 * t45) * t236); m(7) * (t46 * t319 - t115 * t65 + t117 * t66 + t197 * t23 + t199 * t24 + (t246 * t9 + t338 * t46) * t236); (-t115 * t199 + t117 * t197 + t255 * t366) * t413;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
