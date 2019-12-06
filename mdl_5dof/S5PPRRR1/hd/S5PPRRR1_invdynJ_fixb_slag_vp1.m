% Calculate vector of inverse dynamics joint torques for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:37
% EndTime: 2019-12-05 15:12:49
% DurationCPUTime: 9.09s
% Computational Cost: add. (17455->586), mult. (16797->906), div. (0->0), fcn. (15937->8), ass. (0->310)
t271 = sin(pkin(8));
t267 = t271 ^ 2;
t272 = cos(pkin(8));
t268 = t272 ^ 2;
t399 = t267 + t268;
t269 = pkin(9) + qJ(3);
t264 = qJ(4) + t269;
t253 = sin(t264);
t254 = cos(t264);
t427 = Icges(5,4) * t254;
t344 = -Icges(5,2) * t253 + t427;
t166 = -Icges(5,6) * t272 + t271 * t344;
t349 = -Icges(5,1) * t253 - t427;
t466 = t349 * t271 - t166;
t167 = Icges(5,6) * t271 + t272 * t344;
t465 = t349 * t272 - t167;
t428 = Icges(5,4) * t253;
t350 = Icges(5,1) * t254 - t428;
t168 = -Icges(5,5) * t272 + t271 * t350;
t343 = -Icges(5,2) * t254 - t428;
t464 = -t343 * t271 - t168;
t169 = Icges(5,5) * t271 + t272 * t350;
t463 = -t343 * t272 - t169;
t274 = sin(qJ(5));
t275 = cos(qJ(5));
t438 = rSges(6,1) * t275;
t362 = -rSges(6,2) * t274 + t438;
t262 = sin(t269);
t263 = cos(t269);
t429 = Icges(4,4) * t263;
t430 = Icges(4,4) * t262;
t462 = (-t262 * (-Icges(4,2) * t263 - t430) + t263 * (-Icges(4,1) * t262 - t429)) * qJD(3);
t260 = qJD(3) * t271;
t239 = qJD(4) * t271 + t260;
t396 = qJD(5) * t253;
t203 = t272 * t396 + t239;
t393 = qJD(5) * t271;
t270 = qJD(3) + qJD(4);
t413 = t270 * t272;
t204 = t253 * t393 - t413;
t395 = qJD(5) * t254;
t409 = t272 * t275;
t412 = t271 * t274;
t223 = -t254 * t412 - t409;
t410 = t272 * t274;
t411 = t271 * t275;
t224 = t254 * t411 - t410;
t418 = t253 * t271;
t93 = Icges(6,5) * t224 + Icges(6,6) * t223 + Icges(6,3) * t418;
t426 = Icges(6,4) * t224;
t95 = Icges(6,2) * t223 + Icges(6,6) * t418 + t426;
t212 = Icges(6,4) * t223;
t97 = Icges(6,1) * t224 + Icges(6,5) * t418 + t212;
t42 = t223 * t95 + t224 * t97 + t418 * t93;
t225 = -t254 * t410 + t411;
t226 = t254 * t409 + t412;
t417 = t253 * t272;
t94 = Icges(6,5) * t226 + Icges(6,6) * t225 + Icges(6,3) * t417;
t425 = Icges(6,4) * t226;
t96 = Icges(6,2) * t225 + Icges(6,6) * t417 + t425;
t213 = Icges(6,4) * t225;
t98 = Icges(6,1) * t226 + Icges(6,5) * t417 + t213;
t43 = t223 * t96 + t224 * t98 + t418 * t94;
t44 = t225 * t95 + t226 * t97 + t417 * t93;
t45 = t225 * t96 + t226 * t98 + t417 * t94;
t336 = Icges(6,5) * t275 - Icges(6,6) * t274;
t146 = -Icges(6,3) * t254 + t253 * t336;
t423 = Icges(6,4) * t275;
t342 = -Icges(6,2) * t274 + t423;
t148 = -Icges(6,6) * t254 + t253 * t342;
t424 = Icges(6,4) * t274;
t348 = Icges(6,1) * t275 - t424;
t150 = -Icges(6,5) * t254 + t253 * t348;
t59 = t146 * t418 + t148 * t223 + t150 * t224;
t60 = t146 * t417 + t148 * t225 + t150 * t226;
t461 = (t203 * t45 + t204 * t44 - t395 * t60) * t272 + (t203 * t43 + t204 * t42 - t395 * t59) * t271;
t228 = rSges(5,1) * t253 + rSges(5,2) * t254;
t201 = t228 * t271;
t160 = t270 * t201;
t202 = t228 * t272;
t161 = t270 * t202;
t454 = t254 * rSges(5,1) - rSges(5,2) * t253;
t171 = -rSges(5,3) * t272 + t271 * t454;
t172 = rSges(5,3) * t271 + t272 * t454;
t251 = pkin(3) * t263;
t122 = -pkin(6) * t272 + t251 * t271;
t123 = pkin(6) * t271 + t251 * t272;
t397 = qJD(3) * t272;
t376 = t122 * t260 + t123 * t397 + qJD(1);
t460 = (-t271 * t160 - t272 * t161 + t239 * t201 + t202 * t413) * (t171 * t239 + t172 * t413 + t376);
t100 = rSges(6,1) * t226 + rSges(6,2) * t225 + rSges(6,3) * t417;
t385 = t253 * t411;
t386 = t253 * t412;
t415 = t254 * t271;
t403 = rSges(6,2) * t386 + rSges(6,3) * t415;
t130 = -rSges(6,1) * t385 + t403;
t383 = t253 * t409;
t384 = t253 * t410;
t414 = t254 * t272;
t402 = rSges(6,2) * t384 + rSges(6,3) * t414;
t131 = -rSges(6,1) * t383 + t402;
t230 = pkin(4) * t253 - pkin(7) * t254;
t316 = t270 * t230;
t162 = t271 * t316;
t163 = t272 * t316;
t243 = pkin(7) * t415;
t244 = pkin(7) * t414;
t378 = t254 * t393;
t453 = t254 * pkin(4) + t253 * pkin(7);
t205 = t453 * t271;
t207 = t453 * t272;
t99 = rSges(6,1) * t224 + rSges(6,2) * t223 + rSges(6,3) * t418;
t38 = -t100 * t204 + t203 * t99 + t205 * t239 + t207 * t413 + t376;
t387 = t99 * t395;
t452 = g(1) * t272 + g(2) * t271;
t132 = -qJD(5) * t224 + t270 * t386;
t133 = qJD(5) * t223 - t270 * t385;
t382 = t270 * t415;
t79 = rSges(6,1) * t133 + rSges(6,2) * t132 + rSges(6,3) * t382;
t134 = -qJD(5) * t226 + t270 * t384;
t135 = qJD(5) * t225 - t270 * t383;
t381 = t254 * t413;
t80 = rSges(6,1) * t135 + rSges(6,2) * t134 + rSges(6,3) * t381;
t459 = ((-t162 + t79) * t271 + t100 * t378 - t203 * t130 + t131 * t204 - (-pkin(4) * t417 + t244) * t413 - t239 * (-pkin(4) * t418 + t243) + (-t163 + t80 - t387) * t272) * t38 - t452 * (-pkin(4) - t438) * t253;
t354 = -t274 * t96 + t275 * t98;
t355 = -t274 * t95 + t275 * t97;
t451 = -(-t146 * t272 - t354) * t203 - (-t146 * t271 - t355) * t204;
t341 = -Icges(6,2) * t275 - t424;
t283 = t203 * (-Icges(6,2) * t226 + t213 + t98) + t204 * (-Icges(6,2) * t224 + t212 + t97) - t395 * (t341 * t253 + t150);
t257 = qJDD(3) * t271;
t237 = qJDD(4) * t271 + t257;
t394 = qJD(5) * t270;
t311 = qJDD(5) * t253 + t254 * t394;
t136 = t272 * t311 + t237;
t450 = t136 / 0.2e1;
t238 = (-qJDD(3) - qJDD(4)) * t272;
t137 = t271 * t311 + t238;
t449 = t137 / 0.2e1;
t448 = -t203 / 0.2e1;
t447 = t203 / 0.2e1;
t446 = -t204 / 0.2e1;
t445 = t204 / 0.2e1;
t214 = -qJDD(5) * t254 + t253 * t394;
t444 = t214 / 0.2e1;
t443 = t271 / 0.2e1;
t442 = -t272 / 0.2e1;
t441 = pkin(3) * t262;
t435 = t271 * t44;
t434 = t272 * t43;
t192 = t453 * t270;
t361 = -rSges(6,1) * t274 - rSges(6,2) * t275;
t416 = t254 * t270;
t87 = t362 * t416 + (rSges(6,3) * t270 + qJD(5) * t361) * t253;
t431 = -t192 - t87;
t419 = t146 * t254;
t112 = t272 * t123;
t408 = t271 * t122 + t112;
t406 = t271 * t171 + t272 * t172;
t158 = -rSges(6,3) * t254 + t253 * t362;
t404 = -t158 - t230;
t398 = qJD(2) * t272;
t392 = qJDD(2) * t272;
t391 = qJD(3) * t251;
t390 = pkin(3) * t260;
t389 = pkin(3) * t397;
t388 = -m(3) - m(4) - m(5) - m(6);
t380 = t243 + t403;
t379 = t244 + t402;
t374 = -t395 / 0.2e1;
t373 = t395 / 0.2e1;
t300 = -t228 - t441;
t370 = (t100 + t207) * t272 + (t205 + t99) * t271;
t369 = t263 * t390;
t368 = t263 * t389;
t366 = t404 - t441;
t185 = t454 * t270;
t365 = -t185 - t391;
t233 = rSges(4,1) * t263 - rSges(4,2) * t262;
t232 = rSges(4,1) * t262 + rSges(4,2) * t263;
t360 = t94 * t203 + t93 * t204;
t258 = qJDD(2) * t271;
t276 = qJD(3) ^ 2;
t310 = pkin(3) * (-qJDD(3) * t262 - t263 * t276);
t295 = t272 * t310 + t258;
t31 = t137 * t158 - t192 * t413 + t204 * t87 - t214 * t99 + t230 * t238 + t395 * t79 + t295;
t290 = t271 * t310 - t392;
t32 = t100 * t214 - t136 * t158 - t192 * t239 - t203 * t87 - t230 * t237 - t395 * t80 + t290;
t359 = t271 * t31 - t272 * t32;
t358 = t271 * t42 + t434;
t357 = t272 * t45 + t435;
t50 = t253 * t355 - t254 * t93;
t51 = t253 * t354 - t254 * t94;
t356 = t50 * t271 + t51 * t272;
t353 = -t100 * t271 + t272 * t99;
t352 = Icges(4,1) * t263 - t430;
t347 = -Icges(6,1) * t274 - t423;
t346 = -Icges(4,2) * t262 + t429;
t340 = Icges(4,5) * t263 - Icges(4,6) * t262;
t339 = -Icges(4,5) * t262 - Icges(4,6) * t263;
t338 = Icges(5,5) * t254 - Icges(5,6) * t253;
t337 = -Icges(5,5) * t253 - Icges(5,6) * t254;
t335 = -Icges(6,5) * t274 - Icges(6,6) * t275;
t334 = -t148 * t274 + t150 * t275;
t333 = -t166 * t253 + t168 * t254;
t332 = -t167 * t253 + t169 * t254;
t177 = -Icges(4,6) * t272 + t271 * t346;
t179 = -Icges(4,5) * t272 + t271 * t352;
t331 = -t177 * t262 + t179 * t263;
t178 = Icges(4,6) * t271 + t272 * t346;
t180 = Icges(4,5) * t271 + t272 * t352;
t330 = -t178 * t262 + t180 * t263;
t261 = qJD(2) * t271;
t329 = -(-t232 * t397 + t261) * t272 - (-t232 * t260 - t398) * t271;
t328 = t399 * t233;
t327 = t399 * qJD(3) * t232;
t195 = t337 * t271;
t196 = t337 * t272;
t326 = -t195 * t413 + t196 * t239;
t325 = -t391 + t431;
t323 = -t262 * t389 + t261;
t159 = t253 * rSges(6,3) + t362 * t254;
t227 = t233 * qJD(3);
t322 = -qJD(3) * t227 - qJDD(3) * t232;
t73 = Icges(6,5) * t133 + Icges(6,6) * t132 + Icges(6,3) * t382;
t321 = t253 * t73 + t416 * t93;
t74 = Icges(6,5) * t135 + Icges(6,6) * t134 + Icges(6,3) * t381;
t320 = t253 * t74 + t416 * t94;
t307 = t336 * t254;
t84 = t270 * t307 + (Icges(6,3) * t270 + qJD(5) * t335) * t253;
t317 = t146 * t416 + t253 * t84;
t315 = Icges(6,3) * t253 + t307 - t334;
t309 = t348 * t254;
t308 = t342 * t254;
t302 = qJD(3) * t339;
t301 = -t262 * t390 - t398;
t299 = t159 + t453;
t298 = t130 * t395 + t158 * t378 + t204 * t159 - t396 * t99 - t413 * t453;
t297 = -(Icges(6,5) * t223 - Icges(6,6) * t224) * t204 - (Icges(6,5) * t225 - Icges(6,6) * t226) * t203 + t335 * t253 * t395;
t296 = -t276 * t399 * t441 + qJDD(3) * t112 + t122 * t257 + qJDD(1);
t291 = t253 * t297;
t289 = (t464 * t253 + t466 * t254) * t270;
t288 = (t463 * t253 + t465 * t254) * t270;
t287 = (-t177 * t263 - t179 * t262) * qJD(3) + t462 * t271;
t286 = (-t178 * t263 - t180 * t262) * qJD(3) + t462 * t272;
t284 = -t159 * t203 - t453 * t239 + t100 * t396 + (-t158 * t272 - t131) * t395;
t282 = (Icges(6,1) * t225 - t425 - t96) * t203 + (Icges(6,1) * t223 - t426 - t95) * t204 - (t347 * t253 - t148) * t395;
t126 = t148 * t271;
t127 = t148 * t272;
t128 = t150 * t271;
t129 = t150 * t272;
t75 = Icges(6,4) * t133 + Icges(6,2) * t132 + Icges(6,6) * t382;
t77 = Icges(6,1) * t133 + Icges(6,4) * t132 + Icges(6,5) * t382;
t14 = (t270 * t355 - t73) * t254 + (t270 * t93 - t274 * t75 + t275 * t77 + (-t274 * t97 - t275 * t95) * qJD(5)) * t253;
t149 = Icges(6,6) * t253 + t308;
t76 = Icges(6,4) * t135 + Icges(6,2) * t134 + Icges(6,6) * t381;
t78 = Icges(6,1) * t135 + Icges(6,4) * t134 + Icges(6,5) * t381;
t15 = (t270 * t354 - t74) * t254 + (t270 * t94 - t274 * t76 + t275 * t78 + (-t274 * t98 - t275 * t96) * qJD(5)) * t253;
t151 = Icges(6,5) * t253 + t309;
t20 = t132 * t95 + t133 * t97 + t223 * t75 + t224 * t77 + t271 * t321;
t21 = t132 * t96 + t133 * t98 + t223 * t76 + t224 * t78 + t271 * t320;
t22 = t134 * t95 + t135 * t97 + t225 * t75 + t226 * t77 + t272 * t321;
t23 = t134 * t96 + t135 * t98 + t225 * t76 + t226 * t78 + t272 * t320;
t62 = t253 * t334 - t419;
t25 = t203 * t51 + t204 * t50 - t395 * t62;
t277 = (-t315 * t395 - t451) * t253;
t279 = (t465 * t239 - t466 * t413) * t254 + (t463 * t239 - t464 * t413) * t253;
t85 = t270 * t308 + (Icges(6,6) * t270 + qJD(5) * t341) * t253;
t86 = t270 * t309 + (Icges(6,5) * t270 + qJD(5) * t347) * t253;
t29 = t132 * t148 + t133 * t150 + t223 * t85 + t224 * t86 + t271 * t317;
t3 = t136 * t43 + t137 * t42 + t20 * t204 + t203 * t21 + t214 * t59 - t29 * t395;
t30 = t134 * t148 + t135 * t150 + t225 * t85 + t226 * t86 + t272 * t317;
t4 = t136 * t45 + t137 * t44 + t203 * t23 + t204 * t22 + t214 * t60 - t30 * t395;
t152 = t270 * t195;
t46 = -t152 * t272 + t271 * t289;
t153 = t270 * t196;
t47 = -t153 * t272 + t271 * t288;
t48 = t152 * t271 + t272 * t289;
t49 = t153 * t271 + t272 * t288;
t164 = -Icges(5,3) * t272 + t271 * t338;
t63 = -t164 * t272 + t271 * t333;
t165 = Icges(5,3) * t271 + t272 * t338;
t64 = -t165 * t272 + t271 * t332;
t65 = t164 * t271 + t272 * t333;
t66 = t165 * t271 + t272 * t332;
t280 = (-t20 * t272 + t21 * t271) * t445 - t25 * t396 / 0.2e1 + (t271 * t45 - t272 * t44) * t450 + (t271 * t43 - t272 * t42) * t449 + t239 * (t271 * t49 - t272 * t48) / 0.2e1 - t413 * (t271 * t47 - t272 * t46) / 0.2e1 + (t271 * t51 - t272 * t50) * t444 - t239 * (t271 * t326 + t272 * t279) / 0.2e1 + t413 * (t271 * t279 - t272 * t326) / 0.2e1 + t237 * (t271 * t66 - t272 * t65) / 0.2e1 + t238 * (t271 * t64 - t272 * t63) / 0.2e1 + ((-t127 * t225 - t129 * t226) * t203 + (-t126 * t225 - t128 * t226) * t204 + (t60 * t253 + (-t149 * t225 - t151 * t226 + t435) * t254) * qJD(5) + (((t45 - t419) * qJD(5) + t360) * t254 + t277) * t272) * t448 + ((-t127 * t223 - t129 * t224) * t203 + (-t126 * t223 - t128 * t224) * t204 + (t59 * t253 + (-t149 * t223 - t151 * t224 + t434) * t254) * qJD(5) + (((t42 - t419) * qJD(5) + t360) * t254 + t277) * t271) * t446 + (((t127 * t274 - t129 * t275 + t94) * t203 + (t126 * t274 - t128 * t275 + t93) * t204 + t62 * qJD(5)) * t253 + ((t315 * t254 + (t149 * t274 - t151 * t275 - t146) * t253 + t356) * qJD(5) + t451) * t254) * t373 + (-t22 * t272 + t23 * t271) * t447 + (t237 * t66 + t238 * t65 + t239 * t49 - t413 * t48 + t4) * t443 + (t237 * t64 + t238 * t63 + t239 * t47 - t413 * t46 + t3) * t442 + (-t14 * t272 + t15 * t271 + t461) * t374;
t222 = t232 * t272;
t221 = t232 * t271;
t216 = t339 * t272;
t215 = t339 * t271;
t211 = t361 * t253;
t187 = t272 * t302;
t186 = t271 * t302;
t176 = Icges(4,3) * t271 + t272 * t340;
t175 = -Icges(4,3) * t272 + t271 * t340;
t120 = rSges(6,1) * t225 - rSges(6,2) * t226;
t119 = rSges(6,1) * t223 - rSges(6,2) * t224;
t110 = -t228 * t239 + t301;
t109 = -t228 * t413 + t323;
t103 = t271 * t322 - t392;
t102 = t272 * t322 + t258;
t89 = qJD(3) * t328 + qJD(1);
t82 = -t185 * t239 - t228 * t237 + t290;
t81 = -t185 * t413 + t228 * t238 + t295;
t61 = -qJD(3) * t327 + qJDD(3) * t328 + qJDD(1);
t58 = -t100 * t395 - t158 * t203 - t230 * t239 + t301;
t57 = t158 * t204 - t230 * t413 + t323 + t387;
t39 = -t160 * t239 - t161 * t413 + t171 * t237 - t172 * t238 + t296;
t28 = (t270 * t334 - t84) * t254 + (t146 * t270 - t274 * t85 + t275 * t86 + (-t148 * t275 - t150 * t274) * qJD(5)) * t253;
t11 = -t100 * t137 + t136 * t99 - t162 * t239 - t163 * t413 + t203 * t79 - t204 * t80 + t205 * t237 - t207 * t238 + t296;
t1 = [(m(2) + m(3)) * qJDD(1) + m(4) * t61 + m(5) * t39 + m(6) * t11 + (-m(2) + t388) * g(3); t388 * (g(1) * t271 - g(2) * t272) + m(4) * (t102 * t271 - t103 * t272) + m(5) * (t271 * t81 - t272 * t82) + m(6) * t359 + m(3) * t399 * qJDD(2); t280 + (t215 * qJD(3) * t268 - t272 * t216 * t260) * t397 / 0.2e1 - (t216 * qJD(3) * t267 - t271 * t215 * t397) * t260 / 0.2e1 + (-g(1) * (-t272 * t441 + t379) - g(2) * (-t271 * t441 + t380) - g(3) * (t251 + t299) + t11 * (t370 + t408) + (t31 * t366 + t325 * t57) * t272 + (t32 * t366 + t325 * t58) * t271 - t57 * (t298 - t368) - t58 * (t284 - t369) + t459) * m(6) + (-g(3) * (t454 + t251) - t452 * t300 - t109 * (-t413 * t454 - t368) - t110 * (-t239 * t454 - t369) + t39 * (t406 + t408) + (t109 * t365 + t300 * t81) * t272 + (t110 * t365 + t300 * t82) * t271 + t460) * m(5) + (-(t89 * (-t221 * t271 - t222 * t272) + t329 * t233) * qJD(3) + t61 * t328 - t89 * t327 + (-t102 * t272 - t103 * t271) * t232 + t329 * t227 + g(1) * t222 + g(2) * t221 - g(3) * t233) * m(4) + 0.2e1 * ((t271 * (t176 * t271 + t272 * t330) - t272 * (t175 * t271 + t272 * t331)) * qJDD(3) + (t271 * (t187 * t271 + t272 * t286) - t272 * (t186 * t271 + t272 * t287)) * qJD(3)) * t443 + 0.2e1 * ((t271 * (-t176 * t272 + t271 * t330) - t272 * (-t175 * t272 + t271 * t331)) * qJDD(3) + (t271 * (-t187 * t272 + t271 * t286) - t272 * (-t186 * t272 + t271 * t287)) * qJD(3)) * t442; t280 + (-g(1) * t379 - g(2) * t380 - g(3) * t299 + t11 * t370 + (t31 * t404 + t431 * t57) * t272 + (t32 * t404 + t431 * t58) * t271 - t284 * t58 - t298 * t57 + t459) * m(6) + (t39 * t406 + (-t271 * t82 - t272 * t81) * t228 + (-t109 * t272 - t110 * t271) * t185 + g(1) * t202 + g(2) * t201 + t460 + (t109 * t413 + t110 * t239 - g(3)) * t454) * m(5); t4 * t417 / 0.2e1 + (t253 * t357 - t254 * t60) * t450 + ((t270 * t357 - t30) * t254 + (t22 * t271 + t23 * t272 + t270 * t60) * t253) * t447 + t3 * t418 / 0.2e1 + (t253 * t358 - t254 * t59) * t449 + ((t270 * t358 - t29) * t254 + (t20 * t271 + t21 * t272 + t270 * t59) * t253) * t445 + t270 * t253 * t25 / 0.2e1 - t254 * (t136 * t51 + t137 * t50 + t14 * t204 + t15 * t203 + t214 * t62 - t28 * t395) / 0.2e1 + (t253 * t356 - t254 * t62) * t444 + ((t270 * t356 - t28) * t254 + (t14 * t271 + t15 * t272 + t270 * t62) * t253) * t374 + (t225 * t283 + t226 * t282 - t272 * t291) * t448 + (t223 * t283 + t224 * t282 - t271 * t291) * t446 + (t297 * t254 + (-t274 * t283 + t282 * t275) * t253) * t373 + t461 * t416 / 0.2e1 + ((-t32 * t100 + t31 * t99 + t57 * t79 - t58 * t80 + (t38 * t353 + (t271 * t57 - t272 * t58) * t158) * t270) * t254 + (t57 * (-t270 * t99 + t271 * t87) + t58 * (t100 * t270 - t272 * t87) + t11 * t353 + t38 * (-t271 * t80 + t272 * t79) + t359 * t158) * t253 - t57 * (t119 * t395 + t204 * t211) - t58 * (-t120 * t395 - t203 * t211) - t38 * (t119 * t203 - t120 * t204) - g(1) * t120 - g(2) * t119 - g(3) * t211) * m(6);];
tau = t1;
