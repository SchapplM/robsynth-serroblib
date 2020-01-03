% Calculate vector of inverse dynamics joint torques for
% S4PRPR7
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:54
% DurationCPUTime: 16.21s
% Computational Cost: add. (4345->507), mult. (12056->774), div. (0->0), fcn. (11311->6), ass. (0->253)
t452 = -Icges(3,5) + Icges(4,4);
t450 = -Icges(3,6) + Icges(4,5);
t451 = -Icges(3,2) - Icges(4,3);
t230 = sin(qJ(2));
t232 = cos(qJ(2));
t442 = t230 * t452 + t232 * t450;
t449 = (Icges(3,4) + Icges(4,6)) * t232;
t448 = Icges(4,1) + Icges(3,3);
t227 = sin(pkin(6));
t228 = cos(pkin(6));
t447 = t450 * t228 + (t230 * t451 + t449) * t227;
t370 = t227 * t232;
t371 = t227 * t230;
t382 = Icges(3,4) * t230;
t446 = t227 * (Icges(3,1) * t232 - t382) + Icges(4,2) * t370 - Icges(4,6) * t371 + t452 * t228;
t445 = t442 * qJD(2);
t444 = t230 * t450 - t232 * t452;
t443 = t228 * t227;
t229 = sin(qJ(4));
t231 = cos(qJ(4));
t311 = rSges(5,1) * t229 + rSges(5,2) * t231;
t426 = t232 * t311;
t226 = t228 ^ 2;
t196 = rSges(3,1) * t230 + rSges(3,2) * t232;
t225 = t227 ^ 2;
t412 = t225 + t226;
t441 = t196 * t412;
t440 = t227 * t444 - t228 * t448;
t439 = t227 * t448 + t228 * t444;
t438 = t445 * t227;
t437 = t445 * t228;
t343 = qJD(2) * t232;
t344 = qJD(2) * t230;
t435 = ((-Icges(4,6) * t230 + t232 * t451 - t382) * t344 + ((Icges(3,1) + Icges(4,2)) * t230 + t449) * t343) * t227 + (t230 * t446 + t232 * t447) * qJD(2);
t433 = t230 * t447 - t232 * t446;
t432 = 0.2e1 * qJD(2);
t431 = 2 * qJDD(2);
t430 = t442 * t227;
t429 = t442 * t228;
t338 = qJD(2) * qJD(3);
t427 = qJDD(3) * t230 + t232 * t338;
t340 = qJD(4) * t232;
t347 = qJD(2) * t227;
t183 = t228 * t340 + t347;
t345 = qJD(2) * t228;
t184 = t227 * t340 - t345;
t366 = t230 * t231;
t159 = -t227 * t229 + t228 * t366;
t367 = t229 * t230;
t160 = t227 * t231 + t228 * t367;
t368 = t228 * t232;
t64 = Icges(5,5) * t160 + Icges(5,6) * t159 + Icges(5,3) * t368;
t380 = Icges(5,4) * t160;
t66 = Icges(5,2) * t159 + Icges(5,6) * t368 + t380;
t155 = Icges(5,4) * t159;
t68 = Icges(5,1) * t160 + Icges(5,5) * t368 + t155;
t20 = t159 * t66 + t160 * t68 + t368 * t64;
t161 = t227 * t366 + t228 * t229;
t162 = -t227 * t367 + t228 * t231;
t65 = -Icges(5,5) * t162 + Icges(5,6) * t161 + Icges(5,3) * t370;
t379 = Icges(5,4) * t162;
t67 = Icges(5,2) * t161 + Icges(5,6) * t370 - t379;
t156 = Icges(5,4) * t161;
t69 = -Icges(5,1) * t162 + Icges(5,5) * t370 + t156;
t21 = t159 * t67 + t160 * t69 + t368 * t65;
t22 = t161 * t66 - t162 * t68 + t370 * t64;
t23 = t161 * t67 - t162 * t69 + t370 * t65;
t341 = qJD(4) * t230;
t378 = Icges(5,4) * t229;
t291 = Icges(5,2) * t231 + t378;
t245 = -Icges(5,6) * t230 + t232 * t291;
t377 = Icges(5,4) * t231;
t296 = Icges(5,1) * t229 + t377;
t246 = -Icges(5,5) * t230 + t232 * t296;
t288 = Icges(5,5) * t229 + Icges(5,6) * t231;
t244 = -Icges(5,3) * t230 + t232 * t288;
t365 = t232 * t244;
t39 = -t159 * t245 - t160 * t246 - t228 * t365;
t40 = -t161 * t245 + t162 * t246 - t227 * t365;
t422 = t228 * (t183 * t20 + t184 * t21 + t341 * t39) + t227 * (t183 * t22 + t184 * t23 + t341 * t40);
t421 = -t344 / 0.2e1;
t310 = -rSges(4,2) * t232 + rSges(4,3) * t230;
t413 = pkin(2) * t232 + qJ(3) * t230;
t411 = g(1) * t228 + g(2) * t227;
t180 = (Icges(5,2) * t229 - t377) * t232;
t236 = t183 * (-Icges(5,2) * t160 + t155 + t68) + t184 * (Icges(5,2) * t162 + t156 + t69) + t341 * (-t246 + t180);
t181 = (-Icges(5,1) * t231 + t378) * t232;
t237 = t183 * (-Icges(5,1) * t159 + t380 + t66) + t184 * (-Icges(5,1) * t161 - t379 + t67) + t341 * (-t245 - t181);
t233 = qJD(2) ^ 2;
t408 = -m(4) - m(5);
t337 = qJD(2) * qJD(4);
t259 = qJDD(4) * t232 - t230 * t337;
t336 = qJDD(2) * t227;
t104 = t228 * t259 + t336;
t407 = t104 / 0.2e1;
t335 = qJDD(2) * t228;
t105 = t227 * t259 - t335;
t406 = t105 / 0.2e1;
t405 = t412 * t421;
t404 = -t183 / 0.2e1;
t403 = t183 / 0.2e1;
t402 = -t184 / 0.2e1;
t401 = t184 / 0.2e1;
t189 = qJDD(4) * t230 + t232 * t337;
t400 = t189 / 0.2e1;
t397 = -t232 / 0.2e1;
t220 = qJD(3) * t230;
t205 = t227 * t220;
t194 = pkin(2) * t230 - qJ(3) * t232;
t268 = qJD(2) * t194;
t102 = -t227 * t268 + t205;
t207 = t228 * t220;
t103 = -t228 * t268 + t207;
t390 = t102 * t227 + t103 * t228;
t386 = t21 * t227;
t385 = t22 * t228;
t223 = t232 * rSges(5,3);
t384 = t39 * t232;
t383 = t40 * t232;
t372 = t244 * t230;
t369 = t228 * t230;
t176 = t413 * t227;
t178 = t413 * t228;
t358 = t176 * t227 + t178 * t228;
t342 = qJD(3) * t232;
t157 = qJD(2) * t413 - t342;
t357 = -qJD(2) * t310 - t157;
t356 = t426 * t227;
t355 = t426 * t228;
t309 = rSges(4,2) * t230 + rSges(4,3) * t232;
t354 = -t194 + t309;
t353 = -t413 - t310;
t352 = t427 * t227;
t351 = t427 * t228;
t350 = rSges(4,2) * t371 + rSges(4,3) * t370;
t349 = rSges(4,2) * t369 + rSges(4,3) * t368;
t210 = qJ(3) * t368;
t313 = -pkin(2) * t369 + t210;
t209 = qJ(3) * t370;
t314 = -pkin(2) * t371 + t209;
t329 = t313 * t345 + t314 * t347 + t220;
t135 = rSges(5,1) * t367 + rSges(5,2) * t366 + t223;
t328 = t227 * t344;
t327 = t228 * t344;
t326 = t229 * t343;
t325 = t231 * t343;
t319 = -t341 / 0.2e1;
t318 = t341 / 0.2e1;
t317 = -pkin(5) * t230 - t194;
t316 = t412 * t230;
t315 = qJD(2) * t354;
t136 = rSges(5,3) * t230 - t426;
t312 = -t136 + t317;
t199 = rSges(3,1) * t232 - rSges(3,2) * t230;
t307 = -t183 * t64 - t184 * t65;
t306 = t20 * t228 + t386;
t305 = t227 * t23 + t385;
t301 = t229 * t68 + t231 * t66;
t24 = t230 * t64 - t232 * t301;
t300 = t229 * t69 + t231 * t67;
t25 = t230 * t65 - t232 * t300;
t304 = t227 * t25 + t228 * t24;
t299 = qJD(2) * t317;
t70 = rSges(5,1) * t160 + rSges(5,2) * t159 + rSges(5,3) * t368;
t36 = -t136 * t183 + t227 * t299 + t341 * t70 + t205;
t71 = -rSges(5,1) * t162 + rSges(5,2) * t161 + rSges(5,3) * t370;
t37 = t136 * t184 + t228 * t299 - t341 * t71 + t207;
t303 = -t227 * t36 - t228 * t37;
t302 = t227 * t70 - t228 * t71;
t281 = t412 * t199;
t126 = rSges(4,1) * t227 + t228 * t310;
t127 = -rSges(4,1) * t228 + t227 * t310;
t280 = t126 * t228 + t127 * t227;
t279 = -t229 * t246 - t231 * t245;
t278 = qJD(2) * t441;
t185 = pkin(3) * t227 + pkin(5) * t368;
t186 = -pkin(3) * t228 + pkin(5) * t370;
t277 = t185 * t228 + t186 * t227;
t182 = (-rSges(5,1) * t231 + rSges(5,2) * t229) * t232;
t75 = qJD(4) * t182 + (t230 * t311 + t223) * qJD(2);
t276 = -pkin(5) * t343 - t157 - t75;
t275 = t176 * t347 + t178 * t345 + qJD(1) - t342;
t90 = qJD(4) * t159 + t228 * t326;
t91 = -qJD(4) * t160 + t228 * t325;
t53 = Icges(5,5) * t90 + Icges(5,6) * t91 - Icges(5,3) * t327;
t274 = t232 * t53 - t344 * t64;
t92 = qJD(4) * t161 + t227 * t326;
t93 = qJD(4) * t162 + t227 * t325;
t54 = Icges(5,5) * t92 + Icges(5,6) * t93 - Icges(5,3) * t328;
t273 = t232 * t54 - t344 * t65;
t128 = Icges(5,3) * t232 + t230 * t288;
t179 = (-Icges(5,5) * t231 + Icges(5,6) * t229) * t232;
t72 = qJD(2) * t128 + qJD(4) * t179;
t272 = t232 * t72 + t244 * t344;
t19 = qJD(2) * t277 + t183 * t71 - t184 * t70 + t275;
t271 = t19 * t302;
t269 = qJD(2) * t309;
t252 = (t128 + t279) * t230;
t249 = t179 * t341 + (Icges(5,5) * t159 - Icges(5,6) * t160) * t183 + (Icges(5,5) * t161 + Icges(5,6) * t162) * t184;
t248 = qJD(2) * t357 + qJDD(2) * t354;
t247 = -qJDD(3) * t232 + t102 * t347 + t103 * t345 + t176 * t336 + t178 * t335 + t230 * t338 + qJDD(1);
t132 = Icges(5,5) * t232 + t230 * t296;
t130 = Icges(5,6) * t232 + t230 * t291;
t243 = t232 * t249;
t238 = -qJD(2) * t157 - qJDD(2) * t194 + (-qJDD(2) * t230 - t232 * t233) * pkin(5);
t235 = (t228 * t244 + t301) * t183 + (t227 * t244 + t300) * t184;
t234 = (qJD(4) * t252 + t235) * t232;
t208 = t228 * t342;
t206 = t227 * t342;
t177 = t196 * t228;
t175 = t196 * t227;
t153 = t228 * t269;
t151 = t227 * t269;
t101 = -rSges(5,3) * t369 + t355;
t100 = -rSges(5,3) * t371 + t356;
t99 = t246 * t228;
t98 = t246 * t227;
t97 = t245 * t228;
t96 = t245 * t227;
t85 = rSges(5,1) * t161 + rSges(5,2) * t162;
t84 = rSges(5,1) * t159 - rSges(5,2) * t160;
t77 = t228 * t315 + t207;
t76 = t227 * t315 + t205;
t74 = qJD(2) * t132 + qJD(4) * t181;
t73 = qJD(2) * t130 + qJD(4) * t180;
t60 = rSges(5,1) * t92 + rSges(5,2) * t93 - rSges(5,3) * t328;
t59 = rSges(5,1) * t90 + rSges(5,2) * t91 - rSges(5,3) * t327;
t58 = Icges(5,1) * t92 + Icges(5,4) * t93 - Icges(5,5) * t328;
t57 = Icges(5,1) * t90 + Icges(5,4) * t91 - Icges(5,5) * t327;
t56 = Icges(5,4) * t92 + Icges(5,2) * t93 - Icges(5,6) * t328;
t55 = Icges(5,4) * t90 + Icges(5,2) * t91 - Icges(5,6) * t327;
t52 = -t232 * t279 - t372;
t43 = t228 * t248 + t351;
t42 = t227 * t248 + t352;
t41 = qJD(2) * t280 + t275;
t38 = -qJD(2) * t278 + qJDD(2) * t281 + qJDD(1);
t18 = t280 * qJDD(2) + (t151 * t227 + t153 * t228) * qJD(2) + t247;
t17 = (qJD(2) * t279 + t72) * t230 + (-qJD(2) * t244 - t229 * t74 - t231 * t73 + (-t229 * t245 + t231 * t246) * qJD(4)) * t232;
t16 = t161 * t73 - t162 * t74 + t227 * t272 - t245 * t93 - t246 * t92;
t15 = t159 * t73 + t160 * t74 + t228 * t272 - t245 * t91 - t246 * t90;
t14 = t105 * t136 + t184 * t75 - t189 * t71 + t228 * t238 - t341 * t60 + t351;
t13 = -t104 * t136 - t183 * t75 + t189 * t70 + t227 * t238 + t341 * t59 + t352;
t12 = t183 * t24 + t184 * t25 + t341 * t52;
t11 = t161 * t56 - t162 * t58 + t227 * t273 + t67 * t93 + t69 * t92;
t10 = t161 * t55 - t162 * t57 + t227 * t274 + t66 * t93 + t68 * t92;
t9 = t159 * t56 + t160 * t58 + t228 * t273 + t67 * t91 + t69 * t90;
t8 = t159 * t55 + t160 * t57 + t228 * t274 + t66 * t91 + t68 * t90;
t7 = -pkin(5) * t233 * t316 + qJDD(2) * t277 + t104 * t71 - t105 * t70 + t183 * t60 - t184 * t59 + t247;
t4 = (qJD(2) * t300 + t54) * t230 + (qJD(2) * t65 - t229 * t58 - t231 * t56 + (t229 * t67 - t231 * t69) * qJD(4)) * t232;
t3 = (qJD(2) * t301 + t53) * t230 + (qJD(2) * t64 - t229 * t57 - t231 * t55 + (t229 * t66 - t231 * t68) * qJD(4)) * t232;
t2 = t10 * t183 + t104 * t22 + t105 * t23 + t11 * t184 + t16 * t341 + t189 * t40;
t1 = t104 * t20 + t105 * t21 + t15 * t341 + t183 * t8 + t184 * t9 + t189 * t39;
t5 = [m(2) * qJDD(1) + m(3) * t38 + m(4) * t18 + m(5) * t7 + (-m(2) - m(3) + t408) * g(3); (((-t229 * t99 - t231 * t97 + t64) * t183 + (-t229 * t98 - t231 * t96 + t65) * t184 + t52 * qJD(4)) * t232 + ((t252 + (-t130 * t231 - t132 * t229 - t244) * t232 - t304) * qJD(4) + t235) * t230) * t319 + (t227 * t24 - t228 * t25) * t400 + (t10 * t227 - t11 * t228) * t401 + ((t161 * t97 - t162 * t99) * t183 + (t161 * t96 - t162 * t98) * t184 + (t383 + (t130 * t161 - t132 * t162 - t385) * t230) * qJD(4) + (((-t23 + t372) * qJD(4) + t307) * t230 + t234) * t227) * t402 + (t227 * t8 - t228 * t9) * t403 + ((t159 * t97 + t160 * t99) * t183 + (t159 * t96 + t160 * t98) * t184 + (t384 + (t130 * t159 + t132 * t160 - t386) * t230) * qJD(4) + (((-t20 + t372) * qJD(4) + t307) * t230 + t234) * t228) * t404 + (t22 * t227 - t228 * t23) * t406 + (t20 * t227 - t21 * t228) * t407 - t12 * t340 / 0.2e1 - (t429 * qJD(2) * t225 - t227 * t345 * t430) * t347 / 0.2e1 + (t430 * t226 * qJD(2) - t228 * t347 * t429) * t345 / 0.2e1 + (t227 * t3 - t228 * t4 + t422) * t318 + (-t37 * (t135 * t184 + t208) - t36 * (-t135 * t183 + t206) - t19 * (t100 * t183 - t101 * t184 + t329) - (t303 * t413 + (-t19 * t316 + t232 * t303) * pkin(5)) * qJD(2) - ((t36 * t70 - t37 * t71) * t232 + (t37 * (-t136 * t227 - t100) + t36 * (t136 * t228 + t101) + t271) * t230) * qJD(4) - g(1) * (t210 + t355) - g(2) * (t209 + t356) - g(3) * (pkin(5) * t232 + t135 + t413) - t411 * t230 * (-rSges(5,3) - pkin(2) - pkin(5)) + t7 * t358 + t19 * (-pkin(5) * t344 * t412 + t390) + (t14 * t312 + t37 * t276 + t7 * (t185 + t70) + t19 * t59) * t228 + (t13 * t312 + t36 * t276 + t7 * (t186 + t71) + t19 * t60) * t227) * m(5) + (t18 * t358 + t41 * t390 + (t18 * t126 + t41 * t153 + t354 * t43 + t357 * t77) * t228 + (t18 * t127 + t41 * t151 + t354 * t42 + t357 * t76) * t227 - g(1) * (t313 + t349) - g(2) * (t314 + t350) + g(3) * t353 - t77 * t208 - t76 * t206 - t41 * t329 - ((t349 * t41 + t353 * t77) * t228 + (t350 * t41 + t353 * t76) * t227) * qJD(2)) * m(4) + (g(1) * t177 + g(2) * t175 - g(3) * t199 + t38 * t281 + (-t278 - (-t175 * t227 - t177 * t228) * qJD(2)) * (qJD(2) * t281 + qJD(1)) + (qJDD(2) * t196 + t199 * t233) * t441) * m(3) + (t1 + (t435 * t226 + (t437 * t227 - t228 * t438) * t227) * t432 + (t433 * t226 + (t439 * t227 - t228 * t440) * t227) * t431) * t227 / 0.2e1 - (t2 + (t438 * t226 + (t435 - t437) * t443) * t432 + (t440 * t226 + (t433 - t439) * t443) * t431) * t228 / 0.2e1; -t408 * g(3) * t232 + 0.2e1 * (t19 * t405 + t397 * t7) * m(5) + 0.2e1 * (t18 * t397 + t405 * t41) * m(4) + (t408 * t411 + m(4) * (qJD(2) * t41 + t227 * t42 + t228 * t43) + m(5) * (qJD(2) * t19 + t13 * t227 + t14 * t228)) * t230; t1 * t368 / 0.2e1 + (t230 * t39 + t232 * t306) * t407 + (t15 * t230 + (t227 * t9 + t228 * t8) * t232 + (-t230 * t306 + t384) * qJD(2)) * t403 + t2 * t370 / 0.2e1 + (t230 * t40 + t232 * t305) * t406 + (t16 * t230 + (t10 * t228 + t11 * t227) * t232 + (-t230 * t305 + t383) * qJD(2)) * t401 + t12 * t343 / 0.2e1 + t230 * (t104 * t24 + t105 * t25 + t17 * t341 + t183 * t3 + t184 * t4 + t189 * t52) / 0.2e1 + (t230 * t52 + t232 * t304) * t400 + (t17 * t230 + (t227 * t4 + t228 * t3) * t232 + (-t230 * t304 + t232 * t52) * qJD(2)) * t318 + (t236 * t159 - t160 * t237 + t228 * t243) * t404 + (t161 * t236 + t162 * t237 + t227 * t243) * t402 + (t249 * t230 + (t237 * t229 - t231 * t236) * t232) * t319 + t422 * t421 + ((t13 * t70 - t14 * t71 + t36 * t59 - t37 * t60 + (t271 + (-t227 * t37 + t228 * t36) * t136) * qJD(2)) * t230 + (t37 * (-qJD(2) * t71 + t227 * t75) + t36 * (qJD(2) * t70 - t228 * t75) - t7 * t302 + t19 * (-t227 * t59 + t228 * t60) + (-t13 * t228 + t14 * t227) * t136) * t232 - t37 * (t182 * t184 - t341 * t85) - t36 * (-t182 * t183 + t341 * t84) - t19 * (t183 * t85 - t184 * t84) - g(1) * t84 - g(2) * t85 - g(3) * t182) * m(5);];
tau = t5;
