% Calculate vector of inverse dynamics joint torques for
% S4RPRR9
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR9_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR9_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:08
% EndTime: 2019-12-31 16:56:27
% DurationCPUTime: 16.12s
% Computational Cost: add. (7068->745), mult. (18347->1009), div. (0->0), fcn. (17182->6), ass. (0->350)
t269 = sin(qJ(1));
t459 = -pkin(1) - pkin(5);
t489 = t459 * t269;
t271 = cos(qJ(3));
t267 = sin(qJ(4));
t431 = rSges(5,2) * t267;
t270 = cos(qJ(4));
t433 = rSges(5,1) * t270;
t327 = -t431 + t433;
t268 = sin(qJ(3));
t429 = rSges(5,3) * t268;
t156 = t271 * t327 + t429;
t272 = cos(qJ(1));
t366 = qJD(4) * t272;
t346 = t271 * t366;
t371 = qJD(3) * t269;
t193 = t346 + t371;
t438 = t271 * pkin(3);
t223 = pkin(6) * t268 + t438;
t368 = qJD(4) * t268;
t237 = qJD(1) + t368;
t401 = t269 * t270;
t403 = t268 * t272;
t171 = t267 * t403 + t401;
t399 = t270 * t272;
t402 = t269 * t267;
t172 = t268 * t399 - t402;
t387 = t172 * rSges(5,1) - t171 * rSges(5,2);
t398 = t271 * t272;
t96 = rSges(5,3) * t398 - t387;
t488 = -t156 * t193 - t223 * t371 + t237 * t96;
t369 = qJD(3) * t272;
t241 = pkin(6) * t398;
t262 = t271 * rSges(5,3);
t441 = pkin(3) * t268;
t330 = -t262 + t441;
t487 = t272 * t330 - t241 + t387;
t163 = Icges(5,4) * t172;
t90 = Icges(5,2) * t171 + Icges(5,6) * t398 - t163;
t162 = Icges(5,4) * t171;
t92 = Icges(5,1) * t172 - Icges(5,5) * t398 - t162;
t326 = -t171 * t90 - t172 * t92;
t169 = -t268 * t402 + t399;
t170 = t267 * t272 + t268 * t401;
t400 = t269 * t271;
t412 = Icges(5,4) * t170;
t88 = Icges(5,2) * t169 - Icges(5,6) * t400 + t412;
t161 = Icges(5,4) * t169;
t91 = Icges(5,1) * t170 - Icges(5,5) * t400 + t161;
t434 = t169 * t88 + t170 * t91;
t85 = Icges(5,5) * t170 + Icges(5,6) * t169 - Icges(5,3) * t400;
t87 = -Icges(5,5) * t172 + Icges(5,6) * t171 + Icges(5,3) * t398;
t484 = t326 + t434 + (-t269 * t85 - t272 * t87) * t271;
t367 = qJD(4) * t271;
t194 = -t269 * t367 + t369;
t27 = t171 * t88 - t172 * t91 + t85 * t398;
t306 = Icges(5,5) * t270 - Icges(5,6) * t267;
t143 = Icges(5,3) * t268 + t271 * t306;
t410 = Icges(5,4) * t270;
t308 = -Icges(5,2) * t267 + t410;
t147 = Icges(5,6) * t268 + t271 * t308;
t411 = Icges(5,4) * t267;
t310 = Icges(5,1) * t270 - t411;
t151 = Icges(5,5) * t268 + t271 * t310;
t52 = t143 * t398 + t147 * t171 - t151 * t172;
t483 = -t194 * t27 - t237 * t52;
t320 = t267 * t90 + t270 * t92;
t33 = t268 * t87 - t271 * t320;
t26 = t169 * t90 - t170 * t92 - t400 * t87;
t263 = t272 * rSges(4,3);
t404 = t268 * t269;
t155 = rSges(4,1) * t404 + rSges(4,2) * t400 + t263;
t219 = t272 * pkin(1) + t269 * qJ(2);
t265 = t272 * pkin(5);
t467 = t265 + t219;
t478 = t155 + t467;
t370 = qJD(3) * t271;
t348 = t269 * t370;
t373 = qJD(1) * t272;
t477 = t268 * t373 + t348;
t232 = Icges(4,4) * t400;
t409 = Icges(4,5) * t272;
t152 = Icges(4,1) * t404 + t232 + t409;
t413 = Icges(4,4) * t271;
t311 = Icges(4,1) * t268 + t413;
t153 = -Icges(4,5) * t269 + t272 * t311;
t209 = -Icges(4,2) * t268 + t413;
t178 = t209 * t272;
t279 = t269 * (t153 + t178) - t272 * (-Icges(4,2) * t404 + t152 + t232);
t414 = Icges(4,4) * t268;
t309 = Icges(4,2) * t271 + t414;
t148 = Icges(4,6) * t272 + t269 * t309;
t149 = -Icges(4,6) * t269 + t272 * t309;
t211 = Icges(4,1) * t271 - t414;
t180 = t211 * t269;
t181 = t211 * t272;
t280 = t269 * (t149 - t181) - t272 * (t148 - t180);
t476 = -t280 * t268 + t279 * t271;
t385 = t209 + t311;
t386 = -t309 + t211;
t475 = (t268 * t385 - t271 * t386) * qJD(1);
t51 = -t143 * t400 + t147 * t169 + t151 * t170;
t473 = t193 * t26 + t51 * t237;
t442 = rSges(5,3) + pkin(6);
t472 = t268 * t442;
t142 = Icges(5,3) * t271 - t268 * t306;
t304 = t147 * t267 - t151 * t270;
t321 = t267 * t88 - t270 * t91;
t274 = t193 * (-t143 * t272 + t320) + t194 * (t143 * t269 + t321) + t237 * (t142 + t304);
t471 = t274 * t271;
t301 = t149 * t271 + t153 * t268;
t468 = t301 * t272;
t256 = t272 * qJ(2);
t215 = pkin(1) * t269 - t256;
t440 = pkin(5) * t269;
t337 = -t215 - t440;
t220 = -rSges(3,2) * t272 + t269 * rSges(3,3);
t264 = t271 * pkin(6);
t222 = t264 - t441;
t218 = rSges(4,1) * t271 - rSges(4,2) * t268;
t184 = t218 * t272;
t329 = rSges(4,1) * t268 + rSges(4,2) * t271;
t107 = -qJD(3) * t184 + (t269 * t329 + t263) * qJD(1);
t254 = qJD(2) * t272;
t164 = qJD(1) * t219 - t254;
t198 = t329 * qJD(3);
t364 = qJD(1) * qJD(3);
t202 = qJDD(3) * t269 + t272 * t364;
t273 = qJD(1) ^ 2;
t365 = qJD(1) * qJD(2);
t382 = qJDD(2) * t269 + t272 * t365;
t297 = -t265 * t273 + t382;
t157 = -t269 * rSges(4,3) + t272 * t329;
t334 = t157 + t337;
t29 = -t198 * t371 + t202 * t218 + (-t107 - t164) * qJD(1) + t334 * qJDD(1) + t297;
t351 = t271 * t373;
t356 = rSges(4,1) * t477 + rSges(4,2) * t351;
t372 = qJD(3) * t268;
t108 = (-rSges(4,2) * t372 - rSges(4,3) * qJD(1)) * t269 + t356;
t251 = qJDD(3) * t272;
t203 = -t269 * t364 + t251;
t375 = qJD(1) * t269;
t253 = qJD(2) * t269;
t381 = qJ(2) * t373 + t253;
t357 = qJD(1) * (-pkin(1) * t375 + t381) + qJDD(1) * t219 + t269 * t365;
t290 = qJDD(1) * t265 - t273 * t440 + t357;
t30 = qJD(1) * t108 + qJDD(1) * t155 - t203 * t218 + (qJD(3) * t198 - qJDD(2)) * t272 + t290;
t466 = t29 * t269 - t30 * t272;
t439 = g(2) * t272;
t465 = g(1) * t269 - t439;
t101 = qJD(1) * t148 - qJD(3) * t178;
t104 = -qJD(3) * t181 + (t269 * t311 + t409) * qJD(1);
t307 = Icges(4,5) * t268 + Icges(4,6) * t271;
t145 = -Icges(4,3) * t269 + t272 * t307;
t377 = qJD(1) * t145;
t78 = t149 * t268 - t153 * t271;
t464 = qJD(3) * t78 + t101 * t271 + t104 * t268 + t377;
t335 = qJD(1) * t268 + qJD(4);
t347 = t271 * t369;
t463 = t269 * t335 - t347;
t196 = t309 * qJD(3);
t197 = t311 * qJD(3);
t207 = Icges(4,5) * t271 - Icges(4,6) * t268;
t298 = t268 * t209 - t211 * t271;
t462 = qJD(1) * t207 + qJD(3) * t298 + t196 * t271 + t197 * t268;
t102 = qJD(1) * t149 + t209 * t371;
t105 = qJD(1) * t153 + qJD(3) * t180;
t302 = t148 * t268 - t152 * t271;
t144 = Icges(4,3) * t272 + t269 * t307;
t378 = qJD(1) * t144;
t461 = qJD(3) * t302 - t102 * t271 - t105 * t268 + t378;
t176 = (-Icges(5,2) * t270 - t411) * t271;
t277 = t193 * (Icges(5,2) * t172 + t162 - t92) + t194 * (-Icges(5,2) * t170 + t161 + t91) + t237 * (t151 + t176);
t179 = (-Icges(5,1) * t267 - t410) * t271;
t276 = t193 * (Icges(5,1) * t171 + t163 - t90) + t194 * (Icges(5,1) * t169 - t412 - t88) - t237 * (t147 - t179);
t349 = t268 * t369;
t374 = qJD(1) * t271;
t288 = -t269 * t374 - t349;
t363 = qJDD(4) * t271;
t109 = qJD(4) * t288 + t272 * t363 + t202;
t458 = t109 / 0.2e1;
t110 = -qJD(1) * t346 + t251 + (-t363 + (-qJD(1) + t368) * qJD(3)) * t269;
t457 = t110 / 0.2e1;
t191 = qJD(3) * t367 + qJDD(4) * t268 + qJDD(1);
t456 = t191 / 0.2e1;
t455 = -t193 / 0.2e1;
t454 = t193 / 0.2e1;
t453 = -t194 / 0.2e1;
t452 = t194 / 0.2e1;
t451 = t202 / 0.2e1;
t450 = t203 / 0.2e1;
t449 = -t237 / 0.2e1;
t448 = t237 / 0.2e1;
t447 = t268 / 0.2e1;
t446 = t269 / 0.2e1;
t445 = -t272 / 0.2e1;
t443 = rSges(3,2) - pkin(1);
t350 = t268 * t371;
t289 = t350 - t351;
t81 = -t237 * t401 + (-t272 * t335 - t348) * t267;
t82 = t335 * t399 + (-t237 * t267 + t270 * t370) * t269;
t43 = Icges(5,5) * t82 + Icges(5,6) * t81 + Icges(5,3) * t289;
t45 = Icges(5,4) * t82 + Icges(5,2) * t81 + Icges(5,6) * t289;
t47 = Icges(5,1) * t82 + Icges(5,4) * t81 + Icges(5,5) * t289;
t7 = (qJD(3) * t321 + t43) * t268 + (qJD(3) * t85 - t267 * t45 + t270 * t47 + (-t267 * t91 - t270 * t88) * qJD(4)) * t271;
t437 = t7 * t194;
t296 = t237 * t272;
t79 = -t267 * t463 + t270 * t296;
t80 = t267 * t296 + t270 * t463;
t42 = Icges(5,5) * t80 + Icges(5,6) * t79 + Icges(5,3) * t288;
t44 = Icges(5,4) * t80 + Icges(5,2) * t79 + Icges(5,6) * t288;
t46 = Icges(5,1) * t80 + Icges(5,4) * t79 + Icges(5,5) * t288;
t8 = (qJD(3) * t320 + t42) * t268 + (qJD(3) * t87 - t267 * t44 + t270 * t46 + (t267 * t92 - t270 * t90) * qJD(4)) * t271;
t436 = t8 * t193;
t146 = Icges(5,6) * t271 - t268 * t308;
t100 = qJD(3) * t146 + qJD(4) * t176;
t150 = Icges(5,5) * t271 - t268 * t310;
t103 = qJD(3) * t150 + qJD(4) * t179;
t173 = (-Icges(5,5) * t267 - Icges(5,6) * t270) * t271;
t97 = qJD(3) * t142 + qJD(4) * t173;
t18 = (qJD(3) * t304 + t97) * t268 + (qJD(3) * t143 - t100 * t267 + t103 * t270 + (-t147 * t270 - t151 * t267) * qJD(4)) * t271;
t56 = t143 * t268 - t271 * t304;
t435 = t18 * t237 + t56 * t191;
t430 = rSges(3,3) * t272;
t189 = t218 * t371;
t70 = qJD(1) * t334 + t189 + t253;
t425 = t272 * t70;
t32 = t268 * t85 - t271 * t321;
t422 = t32 * t110;
t421 = t33 * t109;
t239 = pkin(3) * t404;
t185 = -pkin(6) * t400 + t239;
t388 = t170 * rSges(5,1) + t169 * rSges(5,2);
t94 = -rSges(5,3) * t400 + t388;
t416 = -t185 - t94;
t187 = pkin(3) * t403 - t241;
t415 = t187 - t96;
t405 = t207 * t272;
t174 = t269 * t207;
t299 = t209 * t271 + t211 * t268;
t84 = t272 * t299 - t174;
t397 = t84 * qJD(1);
t183 = (-rSges(5,1) * t267 - rSges(5,2) * t270) * t271;
t106 = qJD(4) * t183 + (-t268 * t327 + t262) * qJD(3);
t201 = t222 * qJD(3);
t396 = t106 + t201;
t389 = t156 + t223;
t216 = rSges(3,2) * t269 + t430;
t384 = -t215 + t216;
t160 = t219 + t220;
t383 = t268 * t431 + t262;
t186 = pkin(3) * t400 + pkin(6) * t404;
t380 = rSges(3,2) * t375 + rSges(3,3) * t373;
t379 = -qJD(1) * t215 + t253;
t376 = qJD(1) * t307;
t362 = -rSges(4,3) + t459;
t361 = t82 * rSges(5,1) + t81 * rSges(5,2) + rSges(5,3) * t350;
t360 = t271 * t433;
t359 = t271 * t431;
t54 = t272 * t144 + t148 * t400 + t152 * t404;
t55 = -t272 * t145 - t149 * t400 - t153 * t404;
t355 = pkin(3) * t477 + pkin(6) * t350;
t353 = t442 * t271;
t345 = -pkin(3) - t433;
t344 = -t375 / 0.2e1;
t343 = -t374 / 0.2e1;
t342 = t373 / 0.2e1;
t341 = -t371 / 0.2e1;
t340 = t371 / 0.2e1;
t339 = -t369 / 0.2e1;
t338 = t369 / 0.2e1;
t333 = t187 + t337;
t332 = t256 + t489;
t331 = rSges(5,1) * t80 + rSges(5,2) * t79;
t221 = rSges(2,1) * t272 - rSges(2,2) * t269;
t217 = rSges(2,1) * t269 + rSges(2,2) * t272;
t303 = t148 * t271 + t152 * t268;
t281 = qJD(1) * t303 + qJD(3) * t174 + t377;
t282 = -qJD(1) * t301 - qJD(3) * t405 + t378;
t325 = (t281 * t269 + t272 * t461) * t272 + (t282 * t269 - t272 * t464) * t269;
t324 = (-t269 * t461 + t281 * t272) * t272 + (t269 * t464 + t282 * t272) * t269;
t25 = -t400 * t85 + t434;
t323 = t25 * t272 + t26 * t269;
t322 = t25 * t269 - t26 * t272;
t28 = t398 * t87 - t326;
t319 = t269 * t28 + t27 * t272;
t318 = t269 * t27 - t272 * t28;
t317 = t269 * t33 + t272 * t32;
t316 = t269 * t32 - t272 * t33;
t315 = t269 * t55 + t272 * t54;
t139 = t269 * t144;
t57 = -t303 * t272 + t139;
t58 = -t269 * t145 + t468;
t314 = t269 * t58 + t272 * t57;
t71 = qJD(1) * t478 - t218 * t369 - t254;
t313 = t269 * t70 - t272 * t71;
t312 = t269 * t96 + t272 * t94;
t305 = t107 * t272 - t108 * t269;
t300 = -t155 * t269 - t157 * t272;
t133 = rSges(5,3) * t404 + (-t359 + t360) * t269;
t287 = t143 * t237 + t193 * t87 + t194 * t85;
t286 = (Icges(5,5) * t169 - Icges(5,6) * t170) * t194 + (Icges(5,5) * t171 + Icges(5,6) * t172) * t193 + t173 * t237;
t278 = t299 * qJD(1) - t307 * qJD(3);
t31 = -t193 * t94 + t194 * t96 + (-t185 * t269 - t187 * t272) * qJD(3);
t40 = qJD(1) * t333 + t253 - t488;
t41 = -t223 * t369 - t156 * t194 + t237 * t94 - t254 + (t185 + t467) * qJD(1);
t275 = t31 * t312 + (-t269 * t41 - t272 * t40) * t156;
t213 = t272 * t359;
t188 = t223 * t272;
t182 = t218 * t269;
t154 = -t268 * t433 + t383;
t134 = t213 + (-t360 - t429) * t272;
t132 = t151 * t272;
t131 = t151 * t269;
t130 = t147 * t272;
t129 = t147 * t269;
t126 = qJD(1) * t160 - t254;
t125 = qJD(1) * t384 + t253;
t120 = -pkin(6) * t351 + t355;
t119 = t288 * pkin(6) + (t268 * t375 - t347) * pkin(3);
t118 = rSges(5,1) * t171 + rSges(5,2) * t172;
t117 = rSges(5,1) * t169 - rSges(5,2) * t170;
t83 = t269 * t299 + t405;
t74 = t83 * qJD(1);
t72 = t300 * qJD(3);
t60 = qJD(1) * t380 + qJDD(1) * t220 - qJDD(2) * t272 + t357;
t59 = t384 * qJDD(1) + (-qJD(1) * t220 - t164) * qJD(1) + t382;
t49 = -rSges(5,3) * t351 + t361;
t48 = rSges(5,3) * t288 + t331;
t39 = -t269 * t462 + t278 * t272;
t38 = t278 * t269 + t272 * t462;
t37 = t301 * qJD(3) - t101 * t268 + t104 * t271;
t36 = -qJD(3) * t303 - t102 * t268 + t105 * t271;
t24 = qJD(3) * t314 - t397;
t23 = qJD(3) * t315 + t74;
t16 = t100 * t169 + t103 * t170 + t143 * t289 + t147 * t81 + t151 * t82 - t400 * t97;
t15 = t100 * t171 - t103 * t172 + t143 * t288 + t147 * t79 + t151 * t80 + t398 * t97;
t14 = t193 * t33 + t194 * t32 + t237 * t56;
t13 = qJD(1) * t120 + qJDD(1) * t185 - t106 * t194 - t110 * t156 + t191 * t94 - t203 * t223 + t237 * t49 + (-qJD(3) * t201 - qJDD(2)) * t272 + t290;
t12 = t201 * t371 + t106 * t193 + t109 * t156 - t191 * t96 + t202 * t223 - t237 * t48 + (-t119 - t164) * qJD(1) + t333 * qJDD(1) + t297;
t11 = t193 * t28 - t483;
t10 = t194 * t25 + t473;
t9 = -t109 * t94 + t110 * t96 - t185 * t202 - t187 * t203 - t193 * t49 + t194 * t48 + (t119 * t272 - t120 * t269) * qJD(3);
t6 = -t87 * t351 + t169 * t44 + t170 * t46 + t81 * t90 - t82 * t92 + (-t271 * t42 + t372 * t87) * t269;
t5 = -t85 * t351 + t169 * t45 + t170 * t47 + t81 * t88 + t82 * t91 + (-t271 * t43 + t372 * t85) * t269;
t4 = -t87 * t349 + t171 * t44 - t172 * t46 + t79 * t90 - t80 * t92 + (t272 * t42 - t375 * t87) * t271;
t3 = -t85 * t349 + t171 * t45 - t172 * t47 + t79 * t88 + t80 * t91 + (t272 * t43 - t375 * t85) * t271;
t2 = t109 * t26 + t110 * t25 + t16 * t237 + t191 * t51 + t193 * t6 + t194 * t5;
t1 = t109 * t28 + t110 * t27 + t15 * t237 + t191 * t52 + t193 * t4 + t194 * t3;
t17 = [(-g(1) * (t332 + t487) + (-g(2) + t13) * (-t269 * t353 + t239 + t388 + t467) + (t337 + t487) * t12 + (t254 - t331 + (t438 + t472) * t369 + (t459 * t272 + (-qJ(2) - t330 + t264) * t269) * qJD(1)) * t40 + (t355 + t361 + t381 + t40 - t379 + (-t353 * t272 - t187 + t440 + t489) * qJD(1) + t488) * t41) * m(5) + (t125 * t254 + t126 * (t380 + t381) + (t125 * t443 * t272 + (t125 * (-rSges(3,3) - qJ(2)) - t126 * pkin(1)) * t269) * qJD(1) - (qJD(1) * t216 - t125 + t379) * t126 + (-g(2) + t60) * t160 + (-g(1) + t59) * (t269 * t443 + t256 + t430)) * m(3) + (-qJD(3) * t299 + t196 * t268 - t197 * t271) * qJD(1) + (-t302 + t83) * t450 + ((-t25 + t484) * t193 + t11 + t483) * t453 + (t397 + (t269 ^ 2 * t145 + (-t139 + t55 + (t145 + t303) * t272) * t272) * qJD(3) + t24) * t339 + (t70 * (rSges(4,1) * t347 - rSges(4,2) * t349 + t254) + t71 * (-rSges(4,2) * t350 + t356 + t381) + (t362 * t425 + (t70 * (-qJ(2) - t329) + t71 * t362) * t269) * qJD(1) - (t189 - t70 + (t157 - t440) * qJD(1) + t379) * t71 + (-g(2) + t30) * t478 + (-g(1) + t29) * (t157 + t332)) * m(4) + (t36 + t39) * t338 + (t15 + t10) * t454 + t421 / 0.2e1 + t422 / 0.2e1 + t435 + t436 / 0.2e1 + t437 / 0.2e1 + ((t28 + t484) * t194 + t473) * t455 + t78 * t451 + t16 * t452 + t51 * t457 + t52 * t458 + (t37 + t38 + t23) * t340 - t202 * t84 / 0.2e1 - m(2) * (-g(1) * t217 + g(2) * t221) + (t74 + ((-t57 + t139 + t55) * t269 + (t58 - t468 + (t145 - t303) * t269 + t54) * t272) * qJD(3)) * t341 + (m(2) * (t217 ^ 2 + t221 ^ 2) - t298 + Icges(2,3) + Icges(3,1)) * qJDD(1); (-m(3) - m(5)) * t465 + 0.2e1 * (t12 * t446 + t13 * t445) * m(5) + 0.2e1 * (t445 * t60 + t446 * t59) * m(3) + (-t465 + t466) * m(4); (qJD(1) * t39 + qJD(3) * t324 + qJDD(1) * t83 + t202 * t55 + t203 * t54 + t2) * t272 / 0.2e1 + (((-t129 * t267 + t131 * t270 + t85) * t194 + (t130 * t267 - t132 * t270 + t87) * t193 + (-t146 * t267 + t150 * t270 + t143) * t237 + t56 * qJD(4)) * t271 + (qJD(4) * t316 + t274) * t268) * t449 + ((t129 * t169 + t131 * t170) * t194 + (-t130 * t169 - t132 * t170) * t193 + (t146 * t169 + t150 * t170) * t237 + (-t26 * t403 + t271 * t51) * qJD(4) + ((qJD(4) * t25 + t287) * t268 - t471) * t269) * t453 + ((t129 * t171 - t131 * t172) * t194 + (-t130 * t171 + t132 * t172) * t193 + (t146 * t171 - t150 * t172) * t237 + (t27 * t404 + t271 * t52) * qJD(4) + ((-qJD(4) * t28 - t287) * t268 + t471) * t272) * t455 + qJD(1) * (t269 * t37 + t272 * t36 + (t269 * t302 + t78 * t272) * qJD(1)) / 0.2e1 + qJDD(1) * (t269 * t78 - t272 * t302) / 0.2e1 + ((qJD(3) * t305 - t155 * t202 - t157 * t203) * t300 + t72 * ((-t155 * t272 + t157 * t269) * qJD(1) + t305) - t313 * t198 + ((t269 * t71 + t425) * qJD(1) + t466) * t218 - (t182 * t71 + t184 * t70) * qJD(1) - (t72 * (-t182 * t269 - t184 * t272) - t313 * t329) * qJD(3) - g(1) * t182 + g(2) * t184 + g(3) * t329) * m(4) + (t366 * t447 + t342) * t11 + ((-t371 * t405 - t376) * t269 + (t475 + (t269 * t174 + t476) * qJD(3)) * t272) * t341 + ((t174 * t369 - t376) * t272 + (-t475 + (-t272 * t405 - t476) * qJD(3)) * t269) * t339 + (qJD(1) * t38 + qJD(3) * t325 - qJDD(1) * t84 + t202 * t58 + t203 * t57 + t1) * t446 - qJD(1) * ((-t268 * t386 - t271 * t385) * qJD(1) + (t268 * t279 + t271 * t280) * qJD(3)) / 0.2e1 - t14 * t367 / 0.2e1 + t23 * t344 + t314 * t451 + (-qJD(1) * t322 + t269 * t6 + t272 * t5) * t452 + (-qJD(1) * t318 + t269 * t4 + t272 * t3) * t454 + t317 * t456 + t323 * t457 + t319 * t458 + (-qJD(1) * t316 + t269 * t8 + t272 * t7) * t448 + t315 * t450 + (t344 - t269 * t368 / 0.2e1) * t10 + ((-t13 * t389 - t41 * t396 - t9 * t415 + t31 * (t119 + t48) + (t31 * t416 + t389 * t40) * qJD(1)) * t272 + (t12 * t389 + t40 * t396 + t9 * t416 + t31 * (-t120 - t49) + (t31 * t415 + t389 * t41) * qJD(1)) * t269 - t40 * (qJD(1) * t188 - t134 * t237 + t154 * t193 + t222 * t371) - t41 * (qJD(1) * t186 + t133 * t237 - t154 * t194 - t222 * t369) - t31 * (-t133 * t193 + t134 * t194 - t186 * t371 - t188 * t369) - ((-t40 * t96 + t41 * t94) * t271 + t275 * t268) * qJD(4) - g(1) * (t133 + t186) - g(2) * t213 - g(3) * (t268 * t345 + t264 + t383) - (t271 * t345 - t472) * t439) * m(5) + ((-t54 * t269 + t55 * t272) * qJD(1) + t324) * t338 + ((-t57 * t269 + t58 * t272) * qJD(1) + t325) * t340 + t24 * t342; -t2 * t400 / 0.2e1 + (t268 * t51 - t271 * t322) * t457 + ((qJD(3) * t322 + t16) * t268 + (-qJD(1) * t323 + qJD(3) * t51 - t269 * t5 + t272 * t6) * t271) * t452 + t1 * t398 / 0.2e1 + (t268 * t52 - t271 * t318) * t458 + ((qJD(3) * t318 + t15) * t268 + (-qJD(1) * t319 + qJD(3) * t52 - t269 * t3 + t272 * t4) * t271) * t454 + t14 * t370 / 0.2e1 + (t421 + t422 + t435 + t436 + t437) * t447 + (t268 * t56 - t271 * t316) * t456 + ((qJD(3) * t316 + t18) * t268 + (-qJD(1) * t317 + qJD(3) * t56 - t269 * t7 + t272 * t8) * t271) * t448 + (t169 * t277 + t170 * t276 - t286 * t400) * t453 + (t277 * t171 - t172 * t276 + t286 * t398) * t455 + (t286 * t268 + (-t267 * t277 + t276 * t270) * t271) * t449 + (t268 * t339 + t269 * t343) * t11 + (t268 * t340 + t272 * t343) * t10 + ((qJD(3) * t275 - t12 * t96 + t13 * t94 - t40 * t48 + t41 * t49) * t268 + (t40 * (-qJD(3) * t96 + t106 * t272) + t41 * (qJD(3) * t94 + t106 * t269) - t9 * t312 + t31 * (-t269 * t48 - t272 * t49 - t373 * t96 + t375 * t94) + (t12 * t272 + t13 * t269 + (-t269 * t40 + t272 * t41) * qJD(1)) * t156) * t271 - t40 * (-t118 * t237 + t183 * t193) - t41 * (t117 * t237 - t183 * t194) - t31 * (-t117 * t193 + t118 * t194) - g(1) * t117 - g(2) * t118 - g(3) * t183) * m(5);];
tau = t17;
