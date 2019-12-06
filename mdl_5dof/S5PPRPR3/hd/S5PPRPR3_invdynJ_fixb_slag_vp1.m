% Calculate vector of inverse dynamics joint torques for
% S5PPRPR3
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:47
% EndTime: 2019-12-05 15:05:20
% DurationCPUTime: 25.69s
% Computational Cost: add. (16094->868), mult. (26076->1269), div. (0->0), fcn. (27738->10), ass. (0->336)
t306 = qJ(3) + pkin(9);
t303 = sin(t306);
t304 = cos(t306);
t307 = sin(pkin(8));
t313 = sin(qJ(3));
t315 = cos(qJ(3));
t468 = (-Icges(4,5) * t313 - Icges(5,5) * t303 - Icges(4,6) * t315 - Icges(5,6) * t304) * t307;
t479 = t468 * qJD(3);
t485 = -Icges(5,3) - Icges(4,3);
t308 = sin(pkin(7));
t309 = cos(pkin(8));
t422 = t308 * t309;
t310 = cos(pkin(7));
t427 = t304 * t310;
t245 = t303 * t422 + t427;
t374 = t304 * t422;
t419 = t310 * t303;
t246 = t374 - t419;
t417 = t310 * t315;
t421 = t308 * t313;
t275 = -t309 * t421 - t417;
t418 = t310 * t313;
t420 = t308 * t315;
t276 = t309 * t420 - t418;
t426 = t307 * t308;
t484 = -Icges(4,5) * t276 - Icges(5,5) * t246 - Icges(4,6) * t275 + Icges(5,6) * t245 + t426 * t485;
t247 = -t304 * t308 + t309 * t419;
t248 = t303 * t308 + t309 * t427;
t277 = -t309 * t418 + t420;
t278 = t309 * t417 + t421;
t425 = t307 * t310;
t483 = -Icges(4,5) * t278 - Icges(5,5) * t248 - Icges(4,6) * t277 + Icges(5,6) * t247 + t425 * t485;
t223 = t245 * qJD(3);
t386 = qJD(3) * t310;
t224 = qJD(3) * t374 - t303 * t386;
t261 = t275 * qJD(3);
t262 = t276 * qJD(3);
t482 = -Icges(4,5) * t261 + Icges(5,5) * t223 + Icges(4,6) * t262 + Icges(5,6) * t224;
t225 = t247 * qJD(3);
t226 = t248 * qJD(3);
t263 = t277 * qJD(3);
t264 = t278 * qJD(3);
t481 = -Icges(4,5) * t263 + Icges(5,5) * t225 + Icges(4,6) * t264 + Icges(5,6) * t226;
t480 = t485 * t309 + (Icges(4,5) * t315 + Icges(5,5) * t304 - Icges(4,6) * t313 - Icges(5,6) * t303) * t307;
t312 = sin(qJ(5));
t314 = cos(qJ(5));
t444 = rSges(6,1) * t314;
t478 = rSges(6,2) * t312 - t444;
t433 = Icges(5,4) * t304;
t218 = -Icges(5,6) * t309 + (-Icges(5,2) * t303 + t433) * t307;
t434 = Icges(5,4) * t303;
t219 = -Icges(5,5) * t309 + (Icges(5,1) * t304 - t434) * t307;
t254 = (-Icges(5,2) * t304 - t434) * t307;
t228 = qJD(3) * t254;
t255 = (-Icges(5,1) * t303 - t433) * t307;
t229 = qJD(3) * t255;
t437 = Icges(4,4) * t315;
t250 = -Icges(4,6) * t309 + (-Icges(4,2) * t313 + t437) * t307;
t438 = Icges(4,4) * t313;
t251 = -Icges(4,5) * t309 + (Icges(4,1) * t315 - t438) * t307;
t280 = (-Icges(4,2) * t315 - t438) * t307;
t266 = qJD(3) * t280;
t281 = (-Icges(4,1) * t313 - t437) * t307;
t267 = qJD(3) * t281;
t477 = t218 * t224 + t219 * t223 + t228 * t245 - t229 * t246 + t250 * t262 - t251 * t261 - t266 * t275 - t267 * t276 - t426 * t479;
t476 = t218 * t226 + t219 * t225 + t228 * t247 - t229 * t248 + t250 * t264 - t251 * t263 - t266 * t277 - t267 * t278 - t425 * t479;
t475 = t479 * t309 + (t228 * t303 - t229 * t304 + t266 * t313 - t267 * t315 + (t218 * t304 + t219 * t303 + t250 * t315 + t251 * t313) * qJD(3)) * t307;
t474 = t218 * t245 - t219 * t246 - t250 * t275 - t251 * t276 - t426 * t480;
t473 = t218 * t247 - t219 * t248 - t250 * t277 - t251 * t278 - t425 * t480;
t472 = t480 * t309 + (t218 * t303 - t219 * t304 + t250 * t313 - t251 * t315) * t307;
t389 = qJD(3) * t307;
t370 = t308 * t389;
t210 = qJD(5) * t245 + t370;
t369 = t307 * t386;
t211 = qJD(5) * t247 + t369;
t383 = qJD(5) * t307;
t388 = qJD(3) * t309;
t282 = t303 * t383 - t388;
t424 = t307 * t312;
t269 = -t304 * t424 - t309 * t314;
t423 = t307 * t314;
t270 = t304 * t423 - t309 * t312;
t429 = t303 * t307;
t199 = -t246 * t312 + t308 * t423;
t200 = t246 * t314 + t308 * t424;
t84 = Icges(6,5) * t200 + Icges(6,6) * t199 + Icges(6,3) * t245;
t432 = Icges(6,4) * t200;
t86 = Icges(6,2) * t199 + Icges(6,6) * t245 + t432;
t192 = Icges(6,4) * t199;
t88 = Icges(6,1) * t200 + Icges(6,5) * t245 + t192;
t31 = t269 * t86 + t270 * t88 + t429 * t84;
t201 = -t248 * t312 + t310 * t423;
t202 = t248 * t314 + t310 * t424;
t85 = Icges(6,5) * t202 + Icges(6,6) * t201 + Icges(6,3) * t247;
t431 = Icges(6,4) * t202;
t87 = Icges(6,2) * t201 + Icges(6,6) * t247 + t431;
t193 = Icges(6,4) * t201;
t89 = Icges(6,1) * t202 + Icges(6,5) * t247 + t193;
t32 = t269 * t87 + t270 * t89 + t429 * t85;
t142 = Icges(6,5) * t270 + Icges(6,6) * t269 + Icges(6,3) * t429;
t430 = Icges(6,4) * t270;
t143 = Icges(6,2) * t269 + Icges(6,6) * t429 + t430;
t257 = Icges(6,4) * t269;
t144 = Icges(6,1) * t270 + Icges(6,5) * t429 + t257;
t53 = t142 * t429 + t143 * t269 + t144 * t270;
t471 = (t210 * t31 + t211 * t32 + t282 * t53) * t304;
t470 = Icges(4,5) * t275 - Icges(5,5) * t245 - Icges(4,6) * t276 - Icges(5,6) * t246;
t469 = Icges(4,5) * t277 - Icges(5,5) * t247 - Icges(4,6) * t278 - Icges(5,6) * t248;
t436 = Icges(5,4) * t246;
t123 = -Icges(5,2) * t245 + Icges(5,6) * t426 + t436;
t435 = Icges(5,4) * t248;
t124 = -Icges(5,2) * t247 + Icges(5,6) * t425 + t435;
t232 = Icges(5,4) * t245;
t125 = Icges(5,1) * t246 + Icges(5,5) * t426 - t232;
t233 = Icges(5,4) * t247;
t126 = Icges(5,1) * t248 + Icges(5,5) * t425 - t233;
t440 = Icges(4,4) * t276;
t163 = Icges(4,2) * t275 + Icges(4,6) * t426 + t440;
t439 = Icges(4,4) * t278;
t164 = Icges(4,2) * t277 + Icges(4,6) * t425 + t439;
t271 = Icges(4,4) * t275;
t165 = Icges(4,1) * t276 + Icges(4,5) * t426 + t271;
t272 = Icges(4,4) * t277;
t166 = Icges(4,1) * t278 + Icges(4,5) * t425 + t272;
t467 = (t483 * t309 + (-t124 * t303 + t126 * t304 - t164 * t313 + t166 * t315) * t307) * t310 + (t484 * t309 + (-t123 * t303 + t125 * t304 - t163 * t313 + t165 * t315) * t307) * t308;
t466 = (-t124 * t247 + t126 * t248 + t164 * t277 + t166 * t278 - t425 * t483) * t310 + (-t123 * t247 + t125 * t248 + t163 * t277 + t165 * t278 - t425 * t484) * t308;
t465 = (-t124 * t245 + t126 * t246 + t164 * t275 + t166 * t276 - t426 * t483) * t310 + (-t123 * t245 + t125 * t246 + t163 * t275 + t165 * t276 - t426 * t484) * t308;
t133 = -Icges(5,4) * t223 - Icges(5,2) * t224;
t134 = -Icges(5,4) * t225 - Icges(5,2) * t226;
t135 = -Icges(5,1) * t223 - Icges(5,4) * t224;
t136 = -Icges(5,1) * t225 - Icges(5,4) * t226;
t173 = Icges(4,4) * t261 - Icges(4,2) * t262;
t174 = Icges(4,4) * t263 - Icges(4,2) * t264;
t175 = Icges(4,1) * t261 - Icges(4,4) * t262;
t176 = Icges(4,1) * t263 - Icges(4,4) * t264;
t464 = (t481 * t309 + (-t134 * t303 + t136 * t304 - t174 * t313 + t176 * t315 + (-t124 * t304 - t126 * t303 - t164 * t315 - t166 * t313) * qJD(3)) * t307) * t310 + (t482 * t309 + (-t133 * t303 + t135 * t304 - t173 * t313 + t175 * t315 + (-t123 * t304 - t125 * t303 - t163 * t315 - t165 * t313) * qJD(3)) * t307) * t308;
t463 = (-t124 * t226 - t126 * t225 - t134 * t247 + t136 * t248 - t164 * t264 + t166 * t263 + t174 * t277 + t176 * t278 - t425 * t481) * t310 + (-t123 * t226 - t125 * t225 - t133 * t247 + t135 * t248 - t163 * t264 + t165 * t263 + t173 * t277 + t175 * t278 - t425 * t482) * t308;
t462 = (-t124 * t224 - t126 * t223 - t134 * t245 + t136 * t246 - t164 * t262 + t166 * t261 + t174 * t275 + t176 * t276 - t426 * t481) * t310 + (-t123 * t224 - t125 * t223 - t133 * t245 + t135 * t246 - t163 * t262 + t165 * t261 + t173 * t275 + t175 * t276 - t426 * t482) * t308;
t461 = -t309 * t468 + t425 * t469 + t426 * t470;
t460 = -m(5) - m(6);
t382 = qJDD(3) * t307;
t366 = t308 * t382;
t119 = qJD(5) * t224 + qJDD(5) * t245 + t366;
t459 = t119 / 0.2e1;
t365 = t310 * t382;
t120 = qJD(5) * t226 + qJDD(5) * t247 + t365;
t458 = t120 / 0.2e1;
t457 = -t210 / 0.2e1;
t456 = t210 / 0.2e1;
t455 = -t211 / 0.2e1;
t454 = t211 / 0.2e1;
t381 = qJDD(3) * t309;
t390 = qJD(3) * t304;
t214 = -t381 + (qJD(5) * t390 + qJDD(5) * t303) * t307;
t453 = t214 / 0.2e1;
t452 = -t282 / 0.2e1;
t451 = t282 / 0.2e1;
t449 = pkin(3) * t313;
t448 = g(3) * t307;
t446 = pkin(3) * t315;
t140 = -pkin(4) * t223 + pkin(6) * t224;
t115 = -qJD(5) * t200 + t223 * t312;
t116 = qJD(5) * t199 - t223 * t314;
t73 = rSges(6,1) * t116 + rSges(6,2) * t115 + rSges(6,3) * t224;
t442 = t140 + t73;
t158 = pkin(4) * t246 + pkin(6) * t245;
t90 = rSges(6,1) * t200 + rSges(6,2) * t199 + rSges(6,3) * t245;
t441 = t158 + t90;
t428 = t304 * t307;
t416 = -Icges(5,1) * t245 - t123 - t436;
t415 = -Icges(5,1) * t247 - t124 - t435;
t414 = Icges(5,2) * t246 - t125 + t232;
t413 = Icges(5,2) * t248 - t126 + t233;
t322 = qJ(4) * t307 + t309 * t446;
t379 = pkin(3) * t421;
t168 = t310 * t322 + t379;
t412 = -rSges(5,1) * t248 + rSges(5,2) * t247 - rSges(5,3) * t425 - t168;
t385 = qJD(4) * t307;
t296 = t310 * t385;
t213 = pkin(3) * t263 + t296;
t411 = rSges(5,1) * t225 + rSges(5,2) * t226 - t213;
t410 = pkin(4) * t225 - pkin(6) * t226 - t213;
t167 = -pkin(3) * t418 + t308 * t322;
t216 = -qJ(4) * t309 + t307 * t446;
t409 = t167 * t309 + t216 * t426;
t274 = t277 * pkin(3);
t408 = rSges(5,1) * t247 + rSges(5,2) * t248 - t274;
t407 = pkin(4) * t247 - pkin(6) * t248 - t274;
t406 = -pkin(4) * t248 - pkin(6) * t247 - t168;
t405 = Icges(4,1) * t275 - t163 - t440;
t404 = Icges(4,1) * t277 - t164 - t439;
t403 = -Icges(4,2) * t276 + t165 + t271;
t402 = -Icges(4,2) * t278 + t166 + t272;
t295 = t308 * t385;
t212 = pkin(3) * t261 + t295;
t378 = qJD(3) * t449;
t384 = qJD(4) * t309;
t284 = -t307 * t378 - t384;
t401 = t212 * t309 + t284 * t426;
t222 = -rSges(5,3) * t309 + (rSges(5,1) * t304 - rSges(5,2) * t303) * t307;
t400 = -t216 - t222;
t260 = (pkin(4) * t304 + pkin(6) * t303) * t307;
t399 = -t216 - t260;
t398 = t218 - t255;
t397 = t219 + t254;
t349 = -rSges(5,1) * t303 - rSges(5,2) * t304;
t323 = t349 * t307;
t231 = qJD(3) * t323;
t396 = -t231 - t284;
t155 = -rSges(5,1) * t245 - rSges(5,2) * t246;
t240 = (-pkin(4) * t303 + pkin(6) * t304) * t389;
t395 = -t240 - t284;
t394 = t250 - t281;
t393 = t251 + t280;
t392 = rSges(6,2) * t303 * t424 + rSges(6,3) * t428;
t391 = qJD(2) * t310;
t380 = qJDD(4) * t307;
t117 = -qJD(5) * t202 + t225 * t312;
t118 = qJD(5) * t201 - t225 * t314;
t74 = rSges(6,1) * t118 + rSges(6,2) * t117 + rSges(6,3) * t226;
t377 = -t74 + t410;
t91 = rSges(6,1) * t202 + rSges(6,2) * t201 + rSges(6,3) * t247;
t376 = -t91 + t406;
t375 = -m(3) - m(4) + t460;
t372 = t303 * t389;
t371 = t304 * t389;
t364 = t389 / 0.2e1;
t157 = -pkin(4) * t245 + pkin(6) * t246;
t305 = t307 ^ 2;
t358 = t305 * t379;
t357 = t295 - t391;
t302 = qJD(2) * t308;
t356 = t167 * t388 + t216 * t370 + t296 + t302;
t177 = rSges(4,1) * t261 - rSges(4,2) * t262;
t283 = (-rSges(4,1) * t313 - rSges(4,2) * t315) * t307;
t268 = qJD(3) * t283;
t301 = qJDD(2) * t308;
t169 = rSges(4,1) * t276 + rSges(4,2) * t275 + rSges(4,3) * t426;
t256 = -t309 * rSges(4,3) + (rSges(4,1) * t315 - rSges(4,2) * t313) * t307;
t325 = t169 * t309 + t256 * t426;
t71 = t301 + t325 * qJDD(3) + (t177 * t309 + t268 * t426) * qJD(3);
t170 = rSges(4,1) * t278 + rSges(4,2) * t277 + rSges(4,3) * t425;
t178 = rSges(4,1) * t263 - rSges(4,2) * t264;
t72 = (-qJD(3) * t178 - qJDD(3) * t170) * t309 + (-qJDD(2) + (-qJD(3) * t268 - qJDD(3) * t256) * t307) * t310;
t337 = t308 * t71 - t310 * t72;
t92 = qJD(3) * t325 + t302;
t93 = -t391 + (-t170 * t309 - t256 * t425) * qJD(3);
t335 = t308 * t92 - t310 * t93;
t334 = t167 * t369 + qJD(1) - t384;
t333 = -Icges(6,1) * t314 + Icges(6,4) * t312;
t332 = -Icges(6,4) * t314 + Icges(6,2) * t312;
t331 = -Icges(6,5) * t314 + Icges(6,6) * t312;
t330 = t169 * t310 - t170 * t308;
t329 = t177 * t310 - t178 * t308;
t105 = t246 * rSges(6,3) + t245 * t478;
t106 = rSges(6,3) * t248 + t247 * t478;
t328 = t167 * t381 + t212 * t388 + t216 * t366 + t284 * t370 + t310 * t380 + t301;
t127 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t426;
t327 = t127 * t309 + t222 * t426;
t326 = t158 * t309 + t260 * t426;
t324 = -qJDD(4) * t309 + t167 * t365 + t212 * t369 + qJDD(1);
t273 = t275 * pkin(3);
t321 = (Icges(6,5) * t199 - Icges(6,6) * t200) * t210 + (Icges(6,5) * t201 - Icges(6,6) * t202) * t211 + (Icges(6,5) * t269 - Icges(6,6) * t270) * t282;
t320 = qJD(3) * t411 + qJDD(3) * t412;
t319 = qJD(3) * t410 + qJDD(3) * t406;
t318 = (-Icges(6,2) * t202 + t193 + t89) * t211 + (-Icges(6,2) * t200 + t192 + t88) * t210 + (-Icges(6,2) * t270 + t144 + t257) * t282;
t317 = (Icges(6,1) * t201 - t431 - t87) * t211 + (Icges(6,1) * t199 - t432 - t86) * t210 + (Icges(6,1) * t269 - t143 - t430) * t282;
t316 = (Icges(6,3) * t248 + t247 * t331 + t312 * t87 - t314 * t89) * t211 + (Icges(6,3) * t246 + t245 * t331 + t312 * t86 - t314 * t88) * t210 + (t143 * t312 - t144 * t314 + (Icges(6,3) * t304 + t303 * t331) * t307) * t282;
t293 = t308 * t380;
t290 = pkin(6) * t428;
t286 = t310 * t305 * t378;
t230 = t273 * t388;
t215 = t273 * t369;
t207 = qJD(5) * t269 - t314 * t372;
t206 = -qJD(5) * t270 + t312 * t372;
t198 = -rSges(6,1) * t303 * t423 + t392;
t197 = (Icges(6,5) * t304 + t303 * t333) * t307;
t196 = (Icges(6,6) * t304 + t303 * t332) * t307;
t194 = t212 * t425;
t191 = rSges(4,1) * t277 - rSges(4,2) * t278;
t190 = rSges(4,1) * t275 - rSges(4,2) * t276;
t182 = rSges(6,1) * t269 - rSges(6,2) * t270;
t145 = rSges(6,1) * t270 + rSges(6,2) * t269 + rSges(6,3) * t429;
t139 = t167 * t425;
t137 = -rSges(5,1) * t223 - rSges(5,2) * t224;
t114 = rSges(6,1) * t201 - rSges(6,2) * t202;
t113 = rSges(6,1) * t199 - rSges(6,2) * t200;
t104 = Icges(6,5) * t248 + t247 * t333;
t103 = Icges(6,5) * t246 + t245 * t333;
t102 = Icges(6,6) * t248 + t247 * t332;
t101 = Icges(6,6) * t246 + t245 * t332;
t97 = rSges(6,1) * t207 + rSges(6,2) * t206 + rSges(6,3) * t371;
t96 = Icges(6,1) * t207 + Icges(6,4) * t206 + Icges(6,5) * t371;
t95 = Icges(6,4) * t207 + Icges(6,2) * t206 + Icges(6,6) * t371;
t94 = Icges(6,5) * t207 + Icges(6,6) * t206 + Icges(6,3) * t371;
t80 = t330 * t389 + qJD(1);
t70 = Icges(6,1) * t118 + Icges(6,4) * t117 + Icges(6,5) * t226;
t69 = Icges(6,1) * t116 + Icges(6,4) * t115 + Icges(6,5) * t224;
t68 = Icges(6,4) * t118 + Icges(6,2) * t117 + Icges(6,6) * t226;
t67 = Icges(6,4) * t116 + Icges(6,2) * t115 + Icges(6,6) * t224;
t66 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t226;
t65 = Icges(6,5) * t116 + Icges(6,6) * t115 + Icges(6,3) * t224;
t55 = (t309 * t412 + t400 * t425) * qJD(3) + t357;
t54 = qJD(3) * t327 + t356;
t52 = qJDD(1) + (qJD(3) * t329 + qJDD(3) * t330) * t307;
t45 = (t127 * t310 + t308 * t412) * t389 + t334;
t44 = t142 * t247 + t143 * t201 + t144 * t202;
t43 = t142 * t245 + t143 * t199 + t144 * t200;
t34 = t293 + t320 * t309 + (-qJDD(2) + (qJD(3) * t396 + qJDD(3) * t400) * t307) * t310;
t33 = t327 * qJDD(3) + (t137 * t309 + t231 * t426) * qJD(3) + t328;
t26 = -t145 * t211 + t282 * t91 + (t309 * t406 + t399 * t425) * qJD(3) + t357;
t25 = qJD(3) * t326 + t145 * t210 - t282 * t90 + t356;
t24 = t201 * t87 + t202 * t89 + t247 * t85;
t23 = t201 * t86 + t202 * t88 + t247 * t84;
t22 = t199 * t87 + t200 * t89 + t245 * t85;
t21 = t199 * t86 + t200 * t88 + t245 * t84;
t20 = ((qJD(3) * t137 + qJDD(3) * t127) * t310 + t320 * t308) * t307 + t324;
t19 = -t210 * t91 + t211 * t90 + (t158 * t310 + t308 * t406) * t389 + t334;
t18 = t143 * t206 + t144 * t207 + t269 * t95 + t270 * t96 + (t142 * t390 + t303 * t94) * t307;
t17 = t117 * t143 + t118 * t144 + t142 * t226 + t201 * t95 + t202 * t96 + t247 * t94;
t16 = t115 * t143 + t116 * t144 + t142 * t224 + t199 * t95 + t200 * t96 + t245 * t94;
t15 = -t120 * t145 - t211 * t97 + t214 * t91 + t282 * t74 + t293 + t319 * t309 + (-qJDD(2) + (qJD(3) * t395 + qJDD(3) * t399) * t307) * t310;
t14 = t119 * t145 + t210 * t97 - t214 * t90 - t282 * t73 + t326 * qJDD(3) + (t140 * t309 + t240 * t426) * qJD(3) + t328;
t13 = t206 * t87 + t207 * t89 + t269 * t68 + t270 * t70 + (t303 * t66 + t390 * t85) * t307;
t12 = t206 * t86 + t207 * t88 + t269 * t67 + t270 * t69 + (t303 * t65 + t390 * t84) * t307;
t10 = t117 * t87 + t118 * t89 + t201 * t68 + t202 * t70 + t226 * t85 + t247 * t66;
t9 = t117 * t86 + t118 * t88 + t201 * t67 + t202 * t69 + t226 * t84 + t247 * t65;
t8 = t115 * t87 + t116 * t89 + t199 * t68 + t200 * t70 + t224 * t85 + t245 * t66;
t7 = t115 * t86 + t116 * t88 + t199 * t67 + t200 * t69 + t224 * t84 + t245 * t65;
t6 = t210 * t23 + t211 * t24 + t282 * t44;
t5 = t21 * t210 + t211 * t22 + t282 * t43;
t4 = -t119 * t91 + t120 * t90 - t210 * t74 + t211 * t73 + ((qJD(3) * t140 + qJDD(3) * t158) * t310 + t319 * t308) * t307 + t324;
t3 = t119 * t31 + t12 * t210 + t120 * t32 + t13 * t211 + t18 * t282 + t214 * t53;
t2 = t10 * t211 + t119 * t23 + t120 * t24 + t17 * t282 + t210 * t9 + t214 * t44;
t1 = t119 * t21 + t120 * t22 + t16 * t282 + t210 * t7 + t211 * t8 + t214 * t43;
t11 = [(m(2) + m(3)) * qJDD(1) + m(4) * t52 + m(5) * t20 + m(6) * t4 + (-m(2) + t375) * g(3); t375 * (g(1) * t308 - g(2) * t310) + m(4) * t337 + m(5) * (t308 * t33 - t310 * t34) + m(6) * (t14 * t308 - t15 * t310) + m(3) * (t308 ^ 2 + t310 ^ 2) * qJDD(2); (-t25 * (-t105 * t282 + t198 * t210 + t230) - t26 * (t106 * t282 - t198 * t211 + t286) - t19 * (t105 * t211 - t106 * t210 + t215) - (t25 * (t145 * t246 - t428 * t90) + t26 * (-t145 * t248 + t428 * t91) + t19 * (-t246 * t91 + t248 * t90)) * qJD(5) - (-t25 * t358 + (t157 * t25 + t26 * t407) * t309 + (t19 * (t157 * t310 + t308 * t407) + (t25 * t308 - t26 * t310) * (-pkin(4) * t429 + t290)) * t307) * qJD(3) - g(1) * (t106 - t407) - g(2) * (t105 + t157 + t273) - g(3) * (t290 + t392) - (-t449 + (-pkin(4) - t444) * t303) * t448 + t14 * t409 + t25 * t401 + t4 * t139 + t19 * t194 + (t14 * t441 + t15 * t376 + t25 * t442 + t26 * t377) * t309 + ((t15 * (-t145 + t399) + t26 * (-t97 + t395) + t4 * t441 + t19 * t442) * t310 + (t14 * (t145 + t260) + t25 * (t240 + t97) + t4 * t376 + t19 * t377) * t308) * t307) * m(6) - (t246 * t5 + t248 * t6) * qJD(5) / 0.2e1 + (-t18 * t309 + (t12 * t308 + t13 * t310) * t307) * t451 + ((t102 * t269 + t104 * t270) * t211 + (t101 * t269 + t103 * t270) * t210 + (t196 * t269 + t197 * t270) * t282 + (t246 * t31 + t248 * t32) * qJD(5) + ((qJD(5) * t53 + t142 * t282 + t210 * t84 + t211 * t85) * t304 + t316 * t303) * t307) * t452 + (-t309 * t53 + (t308 * t31 + t310 * t32) * t307) * t453 + (-t17 * t309 + (t10 * t310 + t308 * t9) * t307) * t454 + ((t102 * t201 + t104 * t202 + t248 * t85) * t211 + (t101 * t201 + t103 * t202 + t248 * t84) * t210 + (t142 * t248 + t196 * t201 + t197 * t202) * t282 + (t23 * t246 + t24 * t248 + t428 * t44) * qJD(5) + t316 * t247) * t455 + (-t16 * t309 + (t308 * t7 + t310 * t8) * t307) * t456 + ((t102 * t199 + t104 * t200 + t246 * t85) * t211 + (t101 * t199 + t103 * t200 + t246 * t84) * t210 + (t142 * t246 + t196 * t199 + t197 * t200) * t282 + (t21 * t246 + t22 * t248 + t428 * t43) * qJD(5) + t316 * t245) * t457 + (-t309 * t44 + (t23 * t308 + t24 * t310) * t307) * t458 + (-t309 * t43 + (t21 * t308 + t22 * t310) * t307) * t459 - (t307 * t467 + t309 * t472) * t381 / 0.2e1 - (t3 + (qJD(3) * t475 + qJDD(3) * t472) * t309 + (qJD(3) * t464 + qJDD(3) * t467) * t307) * t309 / 0.2e1 - (t307 * t464 + t309 * t475) * t388 / 0.2e1 + (t2 + (qJD(3) * t476 + qJDD(3) * t473) * t309 + (qJD(3) * t463 + qJDD(3) * t466) * t307) * t425 / 0.2e1 + (t1 + (qJD(3) * t477 + qJDD(3) * t474) * t309 + (qJD(3) * t462 + qJDD(3) * t465) * t307) * t426 / 0.2e1 - t383 * t471 / 0.2e1 + (((t303 * t397 + t304 * t398 - t308 * t470 - t310 * t469 + t313 * t393 + t315 * t394) * t309 + ((t308 * t416 + t310 * t415) * t304 + (t308 * t414 + t310 * t413) * t303 + (t308 * t405 + t310 * t404) * t315 + (-t308 * t403 - t310 * t402) * t313) * t307) * t389 + t309 ^ 2 * t479) * t388 / 0.2e1 + (-t45 * t215 - t54 * t230 - t55 * t286 - (-t54 * t358 + (t155 * t54 + t408 * t55) * t309 + (t45 * (t155 * t310 + t308 * t408) + (t308 * t54 - t310 * t55) * t323) * t307) * qJD(3) + t33 * t409 + t54 * t401 + t20 * t139 + t45 * t194 + (t33 * t127 + t54 * t137 + t34 * t412 + t411 * t55) * t309 + ((t127 * t20 + t137 * t45 + t34 * t400 + t396 * t55) * t310 + (t20 * t412 + t222 * t33 + t231 * t54 + t411 * t45) * t308) * t307 + g(1) * t408 - g(2) * (t273 + t155) - (t349 - t449) * t448) * m(5) + ((t169 * t71 - t170 * t72 + t177 * t92 - t178 * t93) * t309 + (t256 * t337 + t268 * t335 + t329 * t80 + t330 * t52) * t307 - ((t190 * t92 - t191 * t93) * t309 + (t80 * (t190 * t310 - t191 * t308) + t335 * t283) * t307) * qJD(3) - g(1) * t191 - g(2) * t190 - g(3) * t283) * m(4) + ((t307 * t463 + t309 * t476) * t310 + (t307 * t462 + t309 * t477) * t308) * t364 - (((t247 * t414 + t248 * t416 + t277 * t403 + t278 * t405) * t426 + (t247 * t397 + t248 * t398 - t277 * t393 + t278 * t394) * t309 + (t247 * t413 + t248 * t415 + t277 * t402 + t278 * t404 + t461) * t425) * t310 + ((t245 * t413 + t246 * t415 + t275 * t402 + t276 * t404) * t425 + (t245 * t397 + t246 * t398 - t275 * t393 + t276 * t394) * t309 + (t245 * t414 + t246 * t416 + t275 * t403 + t276 * t405 + t461) * t426) * t308) * t307 * qJD(3) ^ 2 / 0.2e1 + ((t307 * t466 + t309 * t473) * t310 + (t307 * t465 + t309 * t474) * t308) * t382 / 0.2e1; t460 * (-g(3) * t309 + (g(1) * t310 + g(2) * t308) * t307) + m(5) * (-t20 * t309 + (t308 * t34 + t310 * t33) * t307) + m(6) * (-t309 * t4 + (t14 * t310 + t15 * t308) * t307); t226 * t6 / 0.2e1 + t247 * t2 / 0.2e1 + (t23 * t245 + t24 * t247 + t429 * t44) * t458 + (t10 * t247 + t224 * t23 + t226 * t24 + t245 * t9 + (t17 * t303 + t390 * t44) * t307) * t454 + t224 * t5 / 0.2e1 + t245 * t1 / 0.2e1 + (t21 * t245 + t22 * t247 + t429 * t43) * t459 + (t21 * t224 + t22 * t226 + t245 * t7 + t247 * t8 + (t16 * t303 + t390 * t43) * t307) * t456 + t364 * t471 + t3 * t429 / 0.2e1 + (t245 * t31 + t247 * t32 + t429 * t53) * t453 + (t12 * t245 + t13 * t247 + t224 * t31 + t226 * t32 + (t18 * t303 + t390 * t53) * t307) * t451 + (t201 * t318 + t202 * t317 + t247 * t321) * t455 + (t199 * t318 + t200 * t317 + t245 * t321) * t457 + (t269 * t318 + t270 * t317 + t321 * t429) * t452 + (t4 * (-t245 * t91 + t247 * t90) + (t245 * t25 - t247 * t26) * t97 + (t14 * t245 - t15 * t247 + t224 * t25 - t226 * t26) * t145 + ((-t25 * t90 + t26 * t91) * t390 + (-t14 * t90 + t15 * t91 - t25 * t73 + t26 * t74) * t303) * t307 - t25 * (-t113 * t282 + t182 * t210) - t26 * (t114 * t282 - t182 * t211) - g(1) * t114 - g(2) * t113 - g(3) * t182 + (-t113 * t211 + t114 * t210 - t224 * t91 + t226 * t90 - t245 * t74 + t247 * t73) * t19) * m(6);];
tau = t11;
