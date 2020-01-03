% Calculate vector of inverse dynamics joint torques for
% S4PRRR5
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:35
% EndTime: 2019-12-31 16:33:45
% DurationCPUTime: 8.32s
% Computational Cost: add. (11817->566), mult. (16405->879), div. (0->0), fcn. (15625->8), ass. (0->298)
t261 = sin(qJ(2));
t263 = cos(qJ(2));
t234 = rSges(3,1) * t261 + rSges(3,2) * t263;
t258 = sin(pkin(7));
t254 = t258 ^ 2;
t259 = cos(pkin(7));
t255 = t259 ^ 2;
t440 = t254 + t255;
t456 = t234 * t440;
t257 = qJ(2) + qJ(3);
t251 = sin(t257);
t252 = cos(t257);
t410 = Icges(4,4) * t252;
t329 = -Icges(4,2) * t251 + t410;
t164 = -Icges(4,6) * t259 + t258 * t329;
t334 = -Icges(4,1) * t251 - t410;
t455 = t334 * t258 - t164;
t165 = Icges(4,6) * t258 + t259 * t329;
t454 = t334 * t259 - t165;
t411 = Icges(4,4) * t251;
t335 = Icges(4,1) * t252 - t411;
t166 = -Icges(4,5) * t259 + t258 * t335;
t328 = -Icges(4,2) * t252 - t411;
t453 = -t328 * t258 - t166;
t167 = Icges(4,5) * t258 + t259 * t335;
t452 = -t328 * t259 - t167;
t260 = sin(qJ(4));
t262 = cos(qJ(4));
t422 = rSges(5,1) * t262;
t346 = -rSges(5,2) * t260 + t422;
t412 = Icges(3,4) * t263;
t413 = Icges(3,4) * t261;
t451 = (-t261 * (-Icges(3,2) * t263 - t413) + t263 * (-Icges(3,1) * t261 - t412)) * qJD(2);
t250 = qJD(2) * t258;
t232 = qJD(3) * t258 + t250;
t380 = qJD(4) * t251;
t179 = t259 * t380 + t232;
t377 = qJD(4) * t258;
t256 = qJD(2) + qJD(3);
t396 = t256 * t259;
t180 = t251 * t377 - t396;
t379 = qJD(4) * t252;
t392 = t259 * t262;
t395 = t258 * t260;
t210 = -t252 * t395 - t392;
t393 = t259 * t260;
t394 = t258 * t262;
t211 = t252 * t394 - t393;
t401 = t251 * t258;
t93 = Icges(5,5) * t211 + Icges(5,6) * t210 + Icges(5,3) * t401;
t409 = Icges(5,4) * t211;
t95 = Icges(5,2) * t210 + Icges(5,6) * t401 + t409;
t185 = Icges(5,4) * t210;
t97 = Icges(5,1) * t211 + Icges(5,5) * t401 + t185;
t42 = t210 * t95 + t211 * t97 + t401 * t93;
t212 = -t252 * t393 + t394;
t213 = t252 * t392 + t395;
t400 = t251 * t259;
t94 = Icges(5,5) * t213 + Icges(5,6) * t212 + Icges(5,3) * t400;
t408 = Icges(5,4) * t213;
t96 = Icges(5,2) * t212 + Icges(5,6) * t400 + t408;
t186 = Icges(5,4) * t212;
t98 = Icges(5,1) * t213 + Icges(5,5) * t400 + t186;
t43 = t210 * t96 + t211 * t98 + t401 * t94;
t44 = t212 * t95 + t213 * t97 + t400 * t93;
t45 = t212 * t96 + t213 * t98 + t400 * t94;
t321 = Icges(5,5) * t262 - Icges(5,6) * t260;
t154 = -Icges(5,3) * t252 + t251 * t321;
t406 = Icges(5,4) * t262;
t327 = -Icges(5,2) * t260 + t406;
t156 = -Icges(5,6) * t252 + t251 * t327;
t407 = Icges(5,4) * t260;
t333 = Icges(5,1) * t262 - t407;
t158 = -Icges(5,5) * t252 + t251 * t333;
t59 = t154 * t401 + t156 * t210 + t158 * t211;
t60 = t154 * t400 + t156 * t212 + t158 * t213;
t450 = (t179 * t45 + t180 * t44 - t379 * t60) * t259 + (t179 * t43 + t180 * t42 - t379 * t59) * t258;
t223 = rSges(4,1) * t251 + rSges(4,2) * t252;
t193 = t223 * t258;
t148 = t256 * t193;
t194 = t223 * t259;
t149 = t256 * t194;
t442 = t252 * rSges(4,1) - rSges(4,2) * t251;
t169 = -rSges(4,3) * t259 + t258 * t442;
t170 = rSges(4,3) * t258 + t259 * t442;
t253 = t263 * pkin(2);
t152 = -pkin(5) * t259 + t253 * t258;
t153 = pkin(5) * t258 + t253 * t259;
t381 = qJD(2) * t259;
t363 = t152 * t250 + t153 * t381 + qJD(1);
t449 = (-t258 * t148 - t259 * t149 + t232 * t193 + t194 * t396) * (t169 * t232 + t170 * t396 + t363);
t100 = rSges(5,1) * t213 + rSges(5,2) * t212 + rSges(5,3) * t400;
t372 = t251 * t394;
t373 = t251 * t395;
t398 = t252 * t258;
t384 = rSges(5,2) * t373 + rSges(5,3) * t398;
t126 = -rSges(5,1) * t372 + t384;
t370 = t251 * t392;
t371 = t251 * t393;
t397 = t252 * t259;
t383 = rSges(5,2) * t371 + rSges(5,3) * t397;
t127 = -rSges(5,1) * t370 + t383;
t225 = pkin(3) * t251 - pkin(6) * t252;
t304 = t256 * t225;
t150 = t258 * t304;
t151 = t259 * t304;
t238 = pkin(6) * t398;
t239 = pkin(6) * t397;
t365 = t252 * t377;
t99 = rSges(5,1) * t211 + rSges(5,2) * t210 + rSges(5,3) * t401;
t374 = t99 * t379;
t441 = t252 * pkin(3) + t251 * pkin(6);
t195 = t441 * t258;
t197 = t441 * t259;
t38 = -t100 * t180 + t179 * t99 + t195 * t232 + t197 * t396 + t363;
t439 = g(1) * t259 + g(2) * t258;
t111 = -qJD(4) * t211 + t256 * t373;
t112 = qJD(4) * t210 - t256 * t372;
t369 = t256 * t398;
t75 = rSges(5,1) * t112 + rSges(5,2) * t111 + rSges(5,3) * t369;
t113 = -qJD(4) * t213 + t256 * t371;
t114 = qJD(4) * t212 - t256 * t370;
t368 = t252 * t396;
t76 = rSges(5,1) * t114 + rSges(5,2) * t113 + rSges(5,3) * t368;
t448 = (t100 * t365 - t179 * t126 + t127 * t180 - (-pkin(3) * t400 + t239) * t396 - t232 * (-pkin(3) * t401 + t238) + (-t150 + t75) * t258 + (-t374 - t151 + t76) * t259) * t38 - t439 * (-pkin(3) - t422) * t251;
t339 = -t260 * t96 + t262 * t98;
t340 = -t260 * t95 + t262 * t97;
t438 = -(-t154 * t259 - t339) * t179 - (-t154 * t258 - t340) * t180;
t326 = -Icges(5,2) * t262 - t407;
t271 = t179 * (-Icges(5,2) * t213 + t186 + t98) + t180 * (-Icges(5,2) * t211 + t185 + t97) - t379 * (t326 * t251 + t158);
t265 = qJD(2) ^ 2;
t248 = qJDD(2) * t258;
t230 = qJDD(3) * t258 + t248;
t378 = qJD(4) * t256;
t299 = qJDD(4) * t251 + t252 * t378;
t122 = t259 * t299 + t230;
t436 = t122 / 0.2e1;
t231 = (-qJDD(2) - qJDD(3)) * t259;
t123 = t258 * t299 + t231;
t435 = t123 / 0.2e1;
t434 = -t179 / 0.2e1;
t433 = t179 / 0.2e1;
t432 = -t180 / 0.2e1;
t431 = t180 / 0.2e1;
t214 = -qJDD(4) * t252 + t251 * t378;
t430 = t214 / 0.2e1;
t429 = t258 / 0.2e1;
t428 = -t259 / 0.2e1;
t427 = pkin(2) * t261;
t419 = pkin(2) * qJD(2);
t418 = t258 * t44;
t417 = t259 * t43;
t184 = t441 * t256;
t345 = -rSges(5,1) * t260 - rSges(5,2) * t262;
t399 = t252 * t256;
t87 = t346 * t399 + (rSges(5,3) * t256 + qJD(4) * t345) * t251;
t414 = -t184 - t87;
t402 = t154 * t252;
t139 = t259 * t153;
t390 = t258 * t152 + t139;
t389 = t258 * t169 + t259 * t170;
t160 = -rSges(5,3) * t252 + t251 * t346;
t387 = -t160 - t225;
t376 = t261 * t419;
t375 = t263 * t419;
t367 = t238 + t384;
t366 = t239 + t383;
t360 = -t379 / 0.2e1;
t359 = t379 / 0.2e1;
t289 = -t223 - t427;
t356 = (t100 + t197) * t259 + (t195 + t99) * t258;
t355 = t258 * t376;
t354 = t258 * t375;
t353 = t259 * t376;
t352 = t259 * t375;
t350 = t387 - t427;
t181 = t442 * t256;
t349 = -t181 - t375;
t235 = rSges(3,1) * t263 - rSges(3,2) * t261;
t344 = t94 * t179 + t93 * t180;
t343 = t258 * t42 + t417;
t342 = t259 * t45 + t418;
t50 = t251 * t340 - t252 * t93;
t51 = t251 * t339 - t252 * t94;
t341 = t50 * t258 + t51 * t259;
t338 = -t100 * t258 + t259 * t99;
t337 = Icges(3,1) * t263 - t413;
t332 = -Icges(5,1) * t260 - t406;
t331 = -Icges(3,2) * t261 + t412;
t325 = Icges(3,5) * t263 - Icges(3,6) * t261;
t324 = -Icges(3,5) * t261 - Icges(3,6) * t263;
t323 = Icges(4,5) * t252 - Icges(4,6) * t251;
t322 = -Icges(4,5) * t251 - Icges(4,6) * t252;
t320 = -Icges(5,5) * t260 - Icges(5,6) * t262;
t319 = -t156 * t260 + t158 * t262;
t318 = -t164 * t251 + t166 * t252;
t317 = -t165 * t251 + t167 * t252;
t175 = -Icges(3,6) * t259 + t258 * t331;
t177 = -Icges(3,5) * t259 + t258 * t337;
t316 = -t175 * t261 + t177 * t263;
t176 = Icges(3,6) * t258 + t259 * t331;
t178 = Icges(3,5) * t258 + t259 * t337;
t315 = -t176 * t261 + t178 * t263;
t314 = t440 * t235;
t187 = t322 * t258;
t188 = t322 * t259;
t313 = -t187 * t396 + t188 * t232;
t312 = qJD(2) * t456;
t311 = -t375 + t414;
t161 = t251 * rSges(5,3) + t346 * t252;
t69 = Icges(5,5) * t112 + Icges(5,6) * t111 + Icges(5,3) * t369;
t309 = t251 * t69 + t399 * t93;
t70 = Icges(5,5) * t114 + Icges(5,6) * t113 + Icges(5,3) * t368;
t308 = t251 * t70 + t399 * t94;
t295 = t321 * t252;
t84 = t256 * t295 + (Icges(5,3) * t256 + qJD(4) * t320) * t251;
t305 = t154 * t399 + t251 * t84;
t303 = Icges(5,3) * t251 + t295 - t319;
t298 = pkin(2) * (-qJDD(2) * t261 - t263 * t265);
t297 = t333 * t252;
t296 = t327 * t252;
t290 = qJD(2) * t324;
t288 = t161 + t441;
t287 = t126 * t379 + t160 * t365 + t180 * t161 - t380 * t99 - t396 * t441;
t286 = t258 * t298;
t285 = t259 * t298;
t283 = -(Icges(5,5) * t210 - Icges(5,6) * t211) * t180 - (Icges(5,5) * t212 - Icges(5,6) * t213) * t179 + t320 * t251 * t379;
t282 = -t265 * t427 * t440 + qJDD(2) * t139 + t152 * t248 + qJDD(1);
t278 = t251 * t283;
t277 = (t453 * t251 + t455 * t252) * t256;
t276 = (t452 * t251 + t454 * t252) * t256;
t275 = (-t175 * t263 - t177 * t261) * qJD(2) + t451 * t258;
t274 = (-t176 * t263 - t178 * t261) * qJD(2) + t451 * t259;
t272 = -t161 * t179 - t441 * t232 + t100 * t380 + (-t160 * t259 - t127) * t379;
t270 = (Icges(5,1) * t212 - t408 - t96) * t179 + (Icges(5,1) * t210 - t409 - t95) * t180 - (t332 * t251 - t156) * t379;
t118 = t156 * t258;
t119 = t156 * t259;
t120 = t158 * t258;
t121 = t158 * t259;
t71 = Icges(5,4) * t112 + Icges(5,2) * t111 + Icges(5,6) * t369;
t73 = Icges(5,1) * t112 + Icges(5,4) * t111 + Icges(5,5) * t369;
t14 = (t256 * t340 - t69) * t252 + (t256 * t93 - t260 * t71 + t262 * t73 + (-t260 * t97 - t262 * t95) * qJD(4)) * t251;
t72 = Icges(5,4) * t114 + Icges(5,2) * t113 + Icges(5,6) * t368;
t74 = Icges(5,1) * t114 + Icges(5,4) * t113 + Icges(5,5) * t368;
t15 = (t256 * t339 - t70) * t252 + (t256 * t94 - t260 * t72 + t262 * t74 + (-t260 * t98 - t262 * t96) * qJD(4)) * t251;
t157 = Icges(5,6) * t251 + t296;
t159 = Icges(5,5) * t251 + t297;
t20 = t111 * t95 + t112 * t97 + t210 * t71 + t211 * t73 + t258 * t309;
t21 = t111 * t96 + t112 * t98 + t210 * t72 + t211 * t74 + t258 * t308;
t22 = t113 * t95 + t114 * t97 + t212 * t71 + t213 * t73 + t259 * t309;
t23 = t113 * t96 + t114 * t98 + t212 * t72 + t213 * t74 + t259 * t308;
t62 = t251 * t319 - t402;
t25 = t179 * t51 + t180 * t50 - t379 * t62;
t266 = (-t303 * t379 - t438) * t251;
t267 = (t454 * t232 - t455 * t396) * t252 + (t452 * t232 - t453 * t396) * t251;
t85 = t256 * t296 + (Icges(5,6) * t256 + qJD(4) * t326) * t251;
t86 = t256 * t297 + (Icges(5,5) * t256 + qJD(4) * t332) * t251;
t28 = t111 * t156 + t112 * t158 + t210 * t85 + t211 * t86 + t258 * t305;
t3 = t122 * t43 + t123 * t42 + t179 * t21 + t180 * t20 + t214 * t59 - t28 * t379;
t29 = t113 * t156 + t114 * t158 + t212 * t85 + t213 * t86 + t259 * t305;
t4 = t122 * t45 + t123 * t44 + t179 * t23 + t180 * t22 + t214 * t60 - t29 * t379;
t142 = t256 * t187;
t46 = -t142 * t259 + t258 * t277;
t143 = t256 * t188;
t47 = -t143 * t259 + t258 * t276;
t48 = t142 * t258 + t259 * t277;
t49 = t143 * t258 + t259 * t276;
t162 = -Icges(4,3) * t259 + t258 * t323;
t63 = -t162 * t259 + t258 * t318;
t163 = Icges(4,3) * t258 + t259 * t323;
t64 = -t163 * t259 + t258 * t317;
t65 = t162 * t258 + t259 * t318;
t66 = t163 * t258 + t259 * t317;
t268 = (-t20 * t259 + t21 * t258) * t431 - t25 * t380 / 0.2e1 + (t258 * t45 - t259 * t44) * t436 + (t258 * t43 - t259 * t42) * t435 + (t258 * t51 - t259 * t50) * t430 + t232 * (t258 * t49 - t259 * t48) / 0.2e1 - t396 * (t258 * t47 - t259 * t46) / 0.2e1 - t232 * (t258 * t313 + t259 * t267) / 0.2e1 + t396 * (t258 * t267 - t259 * t313) / 0.2e1 + t230 * (t258 * t66 - t259 * t65) / 0.2e1 + t231 * (t258 * t64 - t259 * t63) / 0.2e1 + ((-t119 * t212 - t121 * t213) * t179 + (-t118 * t212 - t120 * t213) * t180 + (t60 * t251 + (-t157 * t212 - t159 * t213 + t418) * t252) * qJD(4) + (((t45 - t402) * qJD(4) + t344) * t252 + t266) * t259) * t434 + ((-t119 * t210 - t121 * t211) * t179 + (-t118 * t210 - t120 * t211) * t180 + (t59 * t251 + (-t157 * t210 - t159 * t211 + t417) * t252) * qJD(4) + (((t42 - t402) * qJD(4) + t344) * t252 + t266) * t258) * t432 + (((t119 * t260 - t121 * t262 + t94) * t179 + (t118 * t260 - t120 * t262 + t93) * t180 + t62 * qJD(4)) * t251 + ((t303 * t252 + (t157 * t260 - t159 * t262 - t154) * t251 + t341) * qJD(4) + t438) * t252) * t359 + (-t22 * t259 + t23 * t258) * t433 + (t230 * t66 + t231 * t65 + t232 * t49 - t396 * t48 + t4) * t429 + (t230 * t64 + t231 * t63 + t232 * t47 - t396 * t46 + t3) * t428 + (-t14 * t259 + t15 * t258 + t450) * t360;
t222 = t234 * t259;
t221 = t234 * t258;
t216 = t324 * t259;
t215 = t324 * t258;
t207 = t345 * t251;
t202 = t259 * t290;
t201 = t258 * t290;
t174 = Icges(3,3) * t258 + t259 * t325;
t173 = -Icges(3,3) * t259 + t258 * t325;
t125 = -t223 * t396 - t353;
t124 = -t223 * t232 - t355;
t110 = rSges(5,1) * t212 - rSges(5,2) * t213;
t109 = rSges(5,1) * t210 - rSges(5,2) * t211;
t82 = -t181 * t396 + t223 * t231 + t285;
t81 = -t181 * t232 - t223 * t230 + t286;
t61 = -qJD(2) * t312 + qJDD(2) * t314 + qJDD(1);
t53 = t160 * t180 - t225 * t396 - t353 + t374;
t52 = -t100 * t379 - t160 * t179 - t225 * t232 - t355;
t39 = -t148 * t232 - t149 * t396 + t169 * t230 - t170 * t231 + t282;
t33 = (t256 * t319 - t84) * t252 + (t154 * t256 - t260 * t85 + t262 * t86 + (-t156 * t262 - t158 * t260) * qJD(4)) * t251;
t31 = t123 * t160 + t180 * t87 - t184 * t396 - t214 * t99 + t225 * t231 + t379 * t75 + t285;
t30 = t100 * t214 - t122 * t160 - t179 * t87 - t184 * t232 - t225 * t230 - t379 * t76 + t286;
t11 = -t100 * t123 + t122 * t99 - t150 * t232 - t151 * t396 + t179 * t75 - t180 * t76 + t195 * t230 - t197 * t231 + t282;
t1 = [m(2) * qJDD(1) + m(3) * t61 + m(4) * t39 + m(5) * t11 + (-m(2) - m(3) - m(4) - m(5)) * g(3); (t215 * qJD(2) * t255 - t216 * t259 * t250) * t381 / 0.2e1 - (t216 * qJD(2) * t254 - t258 * t215 * t381) * t250 / 0.2e1 + t268 + (t11 * (t356 + t390) + (t31 * t350 + t311 * t53) * t259 + (t30 * t350 + t311 * t52) * t258 - t53 * (t287 - t352) - t52 * (t272 - t354) - g(1) * (-t259 * t427 + t366) - g(2) * (-t258 * t427 + t367) - g(3) * (t253 + t288) + t448) * m(5) + (-t125 * (-t396 * t442 - t352) - t124 * (-t232 * t442 - t354) + t39 * (t389 + t390) + (t125 * t349 + t289 * t82) * t259 + (t124 * t349 + t289 * t81) * t258 - g(3) * (t442 + t253) - t439 * t289 + t449) * m(4) + (g(1) * t222 + g(2) * t221 - g(3) * t235 + t61 * t314 + (-(-t221 * t258 - t222 * t259) * qJD(2) - t312) * (qJD(2) * t314 + qJD(1)) + (qJDD(2) * t234 + t235 * t265) * t456) * m(3) + 0.2e1 * ((t258 * (t174 * t258 + t259 * t315) - t259 * (t173 * t258 + t259 * t316)) * qJDD(2) + (t258 * (t202 * t258 + t259 * t274) - t259 * (t201 * t258 + t259 * t275)) * qJD(2)) * t429 + 0.2e1 * ((t258 * (-t174 * t259 + t258 * t315) - t259 * (-t173 * t259 + t258 * t316)) * qJDD(2) + (t258 * (-t202 * t259 + t258 * t274) - t259 * (-t201 * t259 + t258 * t275)) * qJD(2)) * t428; t268 + (-t52 * t272 - t53 * t287 + t11 * t356 + (t31 * t387 + t414 * t53) * t259 + (t30 * t387 + t414 * t52) * t258 - g(1) * t366 - g(2) * t367 - g(3) * t288 + t448) * m(5) + (t39 * t389 + (-t258 * t81 - t259 * t82) * t223 + (-t124 * t258 - t125 * t259) * t181 + g(1) * t194 + g(2) * t193 + t449 + (t124 * t232 + t125 * t396 - g(3)) * t442) * m(4); t4 * t400 / 0.2e1 + (t251 * t342 - t60 * t252) * t436 + ((t256 * t342 - t29) * t252 + (t22 * t258 + t23 * t259 + t256 * t60) * t251) * t433 + t3 * t401 / 0.2e1 + (t251 * t343 - t252 * t59) * t435 + ((t256 * t343 - t28) * t252 + (t20 * t258 + t21 * t259 + t256 * t59) * t251) * t431 + t256 * t251 * t25 / 0.2e1 - t252 * (t122 * t51 + t123 * t50 + t14 * t180 + t15 * t179 + t214 * t62 - t33 * t379) / 0.2e1 + (t251 * t341 - t252 * t62) * t430 + ((t256 * t341 - t33) * t252 + (t14 * t258 + t15 * t259 + t256 * t62) * t251) * t360 + (t212 * t271 + t213 * t270 - t259 * t278) * t434 + (t210 * t271 + t211 * t270 - t258 * t278) * t432 + (t283 * t252 + (-t260 * t271 + t270 * t262) * t251) * t359 + t450 * t399 / 0.2e1 + ((-t30 * t100 + t31 * t99 - t52 * t76 + t53 * t75 + (t38 * t338 + (t258 * t53 - t259 * t52) * t160) * t256) * t252 + (t53 * (-t256 * t99 + t258 * t87) + t52 * (t100 * t256 - t259 * t87) + t11 * t338 + t38 * (-t258 * t76 + t259 * t75) + (t258 * t31 - t259 * t30) * t160) * t251 - t53 * (t109 * t379 + t180 * t207) - t52 * (-t110 * t379 - t179 * t207) - t38 * (t109 * t179 - t110 * t180) - g(1) * t110 - g(2) * t109 - g(3) * t207) * m(5);];
tau = t1;
