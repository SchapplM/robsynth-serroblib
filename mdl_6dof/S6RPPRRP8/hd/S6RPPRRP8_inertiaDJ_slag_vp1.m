% Calculate time derivative of joint inertia matrix for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP8_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:28
% EndTime: 2019-03-09 02:15:42
% DurationCPUTime: 11.41s
% Computational Cost: add. (18906->752), mult. (28328->1068), div. (0->0), fcn. (27116->8), ass. (0->355)
t252 = pkin(9) + qJ(4);
t243 = sin(t252);
t244 = cos(t252);
t258 = sin(qJ(5));
t260 = cos(qJ(5));
t298 = Icges(6,5) * t260 - Icges(6,6) * t258;
t173 = Icges(6,3) * t243 + t244 * t298;
t301 = Icges(7,4) * t260 + Icges(7,6) * t258;
t174 = Icges(7,2) * t243 + t244 * t301;
t467 = t174 + t173;
t261 = cos(qJ(1));
t378 = t261 * t260;
t353 = t243 * t378;
t365 = qJD(1) * t261;
t259 = sin(qJ(1));
t322 = qJD(1) * t243 + qJD(5);
t360 = qJD(4) * t261;
t334 = t244 * t360;
t446 = t322 * t259 - t334;
t126 = -qJD(5) * t353 + t258 * t446 - t260 * t365;
t323 = qJD(5) * t243 + qJD(1);
t382 = t258 * t261;
t127 = t260 * t446 + t323 * t382;
t379 = t259 * t260;
t204 = t243 * t382 + t379;
t410 = rSges(7,3) + qJ(6);
t462 = rSges(7,1) + pkin(5);
t466 = t204 * qJD(6) - t410 * t126 - t127 * t462;
t380 = t259 * t258;
t205 = -t353 + t380;
t465 = -t204 * t410 + t205 * t462;
t202 = t243 * t380 - t378;
t203 = t243 * t379 + t382;
t386 = t244 * t259;
t142 = Icges(6,5) * t203 - Icges(6,6) * t202 - Icges(6,3) * t386;
t146 = Icges(6,4) * t203 - Icges(6,2) * t202 - Icges(6,6) * t386;
t150 = Icges(6,1) * t203 - Icges(6,4) * t202 - Icges(6,5) * t386;
t53 = -t142 * t386 - t146 * t202 + t150 * t203;
t384 = t244 * t261;
t143 = Icges(6,5) * t205 + Icges(6,6) * t204 + Icges(6,3) * t384;
t147 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t384;
t151 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t384;
t54 = -t143 * t386 - t147 * t202 + t151 * t203;
t312 = t259 * t53 - t261 * t54;
t140 = Icges(7,5) * t203 - Icges(7,6) * t386 + Icges(7,3) * t202;
t144 = Icges(7,4) * t203 - Icges(7,2) * t386 + Icges(7,6) * t202;
t148 = Icges(7,1) * t203 - Icges(7,4) * t386 + Icges(7,5) * t202;
t51 = t140 * t202 - t144 * t386 + t148 * t203;
t141 = Icges(7,5) * t205 + Icges(7,6) * t384 - Icges(7,3) * t204;
t145 = Icges(7,4) * t205 + Icges(7,2) * t384 - Icges(7,6) * t204;
t149 = Icges(7,1) * t205 + Icges(7,4) * t384 - Icges(7,5) * t204;
t52 = t141 * t202 - t145 * t386 + t149 * t203;
t313 = t259 * t51 - t261 * t52;
t403 = Icges(7,5) * t260;
t297 = Icges(7,3) * t258 + t403;
t172 = Icges(7,6) * t243 + t244 * t297;
t404 = Icges(7,5) * t258;
t305 = Icges(7,1) * t260 + t404;
t176 = Icges(7,4) * t243 + t244 * t305;
t86 = t172 * t202 - t174 * t386 + t176 * t203;
t406 = Icges(6,4) * t260;
t302 = -Icges(6,2) * t258 + t406;
t175 = Icges(6,6) * t243 + t244 * t302;
t407 = Icges(6,4) * t258;
t306 = Icges(6,1) * t260 - t407;
t177 = Icges(6,5) * t243 + t244 * t306;
t87 = -t173 * t386 - t175 * t202 + t177 * t203;
t461 = (-t312 - t313) * t244 + (t86 + t87) * t243;
t57 = t142 * t384 + t204 * t146 + t205 * t150;
t58 = t143 * t384 + t204 * t147 + t205 * t151;
t310 = t259 * t57 - t261 * t58;
t55 = -t204 * t140 + t144 * t384 + t205 * t148;
t56 = -t204 * t141 + t145 * t384 + t205 * t149;
t311 = t259 * t55 - t261 * t56;
t88 = -t204 * t172 + t174 * t384 + t205 * t176;
t89 = t173 * t384 + t204 * t175 + t205 * t177;
t460 = (-t310 - t311) * t244 + (t88 + t89) * t243;
t289 = t172 * t258 + t176 * t260;
t459 = (-t175 * t258 + t177 * t260 + t289) * t244 + t467 * t243;
t359 = qJD(5) * t244;
t120 = (Icges(7,3) * t260 - t404) * t359 + (Icges(7,6) * t244 - t243 * t297) * qJD(4);
t121 = (-Icges(6,5) * t258 - Icges(6,6) * t260) * t359 + (Icges(6,3) * t244 - t243 * t298) * qJD(4);
t122 = (-Icges(7,4) * t258 + Icges(7,6) * t260) * t359 + (Icges(7,2) * t244 - t243 * t301) * qJD(4);
t124 = (-Icges(7,1) * t258 + t403) * t359 + (Icges(7,4) * t244 - t243 * t305) * qJD(4);
t125 = (-Icges(6,1) * t258 - t406) * t359 + (Icges(6,5) * t244 - t243 * t306) * qJD(4);
t358 = qJD(5) * t260;
t332 = t244 * t358;
t333 = t258 * t359;
t361 = qJD(4) * t260;
t337 = t243 * t361;
t364 = qJD(4) * t243;
t338 = t258 * t364;
t363 = qJD(4) * t244;
t387 = t244 * t258;
t458 = t120 * t387 + t172 * t332 + t175 * t338 - t176 * t333 - t177 * t337 + (t124 + t125) * t244 * t260 + t467 * t363 + (t121 + t122) * t243;
t362 = qJD(4) * t259;
t335 = t244 * t362;
t457 = t243 * t365 + t335;
t336 = t243 * t360;
t366 = qJD(1) * t259;
t456 = t244 * t366 + t336;
t293 = t141 * t258 + t149 * t260;
t60 = t145 * t243 + t244 * t293;
t291 = t147 * t258 - t151 * t260;
t62 = t143 * t243 - t244 * t291;
t418 = t60 + t62;
t294 = t140 * t258 + t148 * t260;
t59 = t144 * t243 + t244 * t294;
t292 = t146 * t258 - t150 * t260;
t61 = t142 * t243 - t244 * t292;
t419 = t59 + t61;
t455 = -t259 * t419 + t261 * t418;
t454 = t259 * t418 + t261 * t419;
t285 = t322 * t261;
t128 = t323 * t379 + (t285 + t335) * t258;
t129 = t260 * t285 + (t244 * t361 - t258 * t323) * t259;
t331 = t243 * t362;
t339 = t244 * t365;
t269 = t331 - t339;
t67 = Icges(7,5) * t129 + Icges(7,6) * t269 + Icges(7,3) * t128;
t71 = Icges(7,4) * t129 + Icges(7,2) * t269 + Icges(7,6) * t128;
t75 = Icges(7,1) * t129 + Icges(7,4) * t269 + Icges(7,5) * t128;
t19 = (-qJD(4) * t294 + t71) * t243 + (qJD(4) * t144 + t258 * t67 + t260 * t75 + (t140 * t260 - t148 * t258) * qJD(5)) * t244;
t69 = Icges(6,5) * t129 - Icges(6,6) * t128 + Icges(6,3) * t269;
t73 = Icges(6,4) * t129 - Icges(6,2) * t128 + Icges(6,6) * t269;
t77 = Icges(6,1) * t129 - Icges(6,4) * t128 + Icges(6,5) * t269;
t21 = (qJD(4) * t292 + t69) * t243 + (qJD(4) * t142 - t258 * t73 + t260 * t77 + (-t146 * t260 - t150 * t258) * qJD(5)) * t244;
t453 = t19 + t21;
t66 = Icges(7,5) * t127 - Icges(7,6) * t456 + Icges(7,3) * t126;
t70 = Icges(7,4) * t127 - Icges(7,2) * t456 + Icges(7,6) * t126;
t74 = Icges(7,1) * t127 - Icges(7,4) * t456 + Icges(7,5) * t126;
t20 = (-qJD(4) * t293 + t70) * t243 + (qJD(4) * t145 + t258 * t66 + t260 * t74 + (t141 * t260 - t149 * t258) * qJD(5)) * t244;
t68 = Icges(6,5) * t127 - Icges(6,6) * t126 - Icges(6,3) * t456;
t72 = Icges(6,4) * t127 - Icges(6,2) * t126 - Icges(6,6) * t456;
t76 = Icges(6,1) * t127 - Icges(6,4) * t126 - Icges(6,5) * t456;
t22 = (qJD(4) * t291 + t68) * t243 + (qJD(4) * t143 - t258 * t72 + t260 * t76 + (-t147 * t260 - t151 * t258) * qJD(5)) * t244;
t452 = t20 + t22;
t451 = -rSges(7,2) * t331 - t202 * qJD(6) - t410 * t128 - t129 * t462;
t408 = Icges(5,4) * t244;
t307 = Icges(5,1) * t243 + t408;
t184 = Icges(5,5) * t261 + t259 * t307;
t392 = t184 * t243;
t409 = Icges(5,4) * t243;
t303 = Icges(5,2) * t244 + t409;
t182 = Icges(5,6) * t261 + t259 * t303;
t395 = t182 * t244;
t288 = t392 + t395;
t450 = t261 * t288;
t319 = rSges(5,1) * t243 + rSges(5,2) * t244;
t276 = t261 * t319;
t374 = rSges(7,2) * t384 + t465;
t449 = t261 * t374;
t448 = t410 * t202 + t203 * t462;
t426 = pkin(4) * t243;
t233 = t261 * t426;
t255 = sin(pkin(9));
t427 = pkin(3) * t255;
t239 = t261 * t427;
t249 = t261 * qJ(2);
t257 = -pkin(7) - qJ(3);
t346 = t259 * t257 + t239 + t249;
t447 = t233 + t346;
t429 = -rSges(5,3) - pkin(1);
t159 = t259 * t429 + t276 + t346;
t250 = t261 * rSges(5,3);
t388 = t243 * t259;
t186 = rSges(5,1) * t388 + rSges(5,2) * t386 + t250;
t251 = t261 * pkin(1);
t368 = t259 * qJ(2) + t251;
t283 = -t257 * t261 + t259 * t427 + t368;
t160 = t283 + t186;
t445 = -t159 * t259 + t160 * t261;
t369 = qJ(2) * t365 + qJD(2) * t259;
t345 = qJD(3) * t261 + t369;
t309 = qJD(1) * t239 + t257 * t366 + t345;
t349 = rSges(5,1) * t457 + rSges(5,2) * t339;
t105 = (-rSges(5,2) * t364 + qJD(1) * t429) * t259 + t309 + t349;
t415 = rSges(5,2) * t243;
t215 = rSges(5,1) * t244 - t415;
t247 = qJD(2) * t261;
t324 = -qJD(3) * t259 + t247;
t321 = t257 * t365 + t324;
t330 = -qJ(2) - t427;
t106 = t215 * t360 + (t429 * t261 + (-t319 + t330) * t259) * qJD(1) + t321;
t444 = -t105 * t261 + t259 * t106;
t299 = Icges(5,5) * t243 + Icges(5,6) * t244;
t443 = -Icges(5,3) * t259 + t261 * t299;
t442 = -Icges(5,6) * t259 + t261 * t303;
t441 = -Icges(5,5) * t259 + t261 * t307;
t376 = -rSges(7,2) * t386 + t448;
t440 = t259 * t374 + t261 * t376;
t439 = 2 * m(5);
t438 = 2 * m(6);
t437 = 2 * m(7);
t253 = t259 ^ 2;
t254 = t261 ^ 2;
t436 = -t243 / 0.2e1;
t434 = t244 / 0.2e1;
t430 = rSges(3,2) - pkin(1);
t428 = m(5) * t215;
t425 = t259 * pkin(1);
t417 = -rSges(7,2) * t456 - t466;
t416 = rSges(7,2) * t339 + t451;
t414 = rSges(7,2) * t244;
t413 = rSges(6,3) * t244;
t412 = t259 * rSges(5,3);
t411 = rSges(4,3) + qJ(3);
t317 = -t205 * rSges(6,1) - t204 * rSges(6,2);
t155 = rSges(6,3) * t384 - t317;
t399 = t155 * t261;
t396 = t182 * t243;
t394 = t442 * t243;
t393 = t442 * t244;
t391 = t184 * t244;
t390 = t441 * t243;
t389 = t441 * t244;
t123 = (-Icges(6,2) * t260 - t407) * t359 + (Icges(6,6) * t244 - t243 * t302) * qJD(4);
t383 = t258 * t123;
t314 = pkin(5) * t260 + qJ(6) * t258;
t315 = rSges(7,1) * t260 + rSges(7,3) * t258;
t377 = -t314 * t364 + (qJD(6) * t258 + (-pkin(5) * t258 + qJ(6) * t260) * qJD(5)) * t244 + (-rSges(7,1) * t258 + rSges(7,3) * t260) * t359 + (-t243 * t315 + t414) * qJD(4);
t370 = t203 * rSges(6,1) - t202 * rSges(6,2);
t153 = -rSges(6,3) * t386 + t370;
t232 = pkin(4) * t388;
t196 = -pkin(8) * t386 + t232;
t375 = -t153 - t196;
t373 = rSges(7,2) * t243 + (t314 + t315) * t244;
t211 = (pkin(8) * t244 - t426) * qJD(4);
t216 = pkin(4) * t244 + pkin(8) * t243;
t372 = t259 * t211 + t216 * t365;
t180 = Icges(5,3) * t261 + t259 * t299;
t367 = qJD(1) * t180;
t356 = -pkin(1) - t411;
t35 = t52 * t259 + t261 * t51;
t36 = t54 * t259 + t261 * t53;
t355 = t35 / 0.2e1 + t36 / 0.2e1;
t37 = t56 * t259 + t261 * t55;
t38 = t58 * t259 + t261 * t57;
t354 = -t37 / 0.2e1 - t38 / 0.2e1;
t351 = t129 * rSges(6,1) - t128 * rSges(6,2) + rSges(6,3) * t331;
t350 = -t196 - t376;
t348 = pkin(4) * t457 + pkin(8) * t331;
t347 = pkin(4) * t334 + pkin(8) * t456;
t344 = t244 * (-rSges(7,2) - pkin(8));
t343 = t244 * (-rSges(6,3) - pkin(8));
t316 = rSges(6,1) * t260 - rSges(6,2) * t258;
t179 = rSges(6,3) * t243 + t244 * t316;
t340 = t179 * t366;
t329 = t459 * t363 + (-t289 * t364 + (-t383 + (-t175 * t260 - t177 * t258) * qJD(5)) * t244 + t458) * t243;
t210 = t319 * qJD(4);
t328 = t210 * (t253 + t254);
t327 = t373 * t259;
t326 = qJD(1) * t373;
t325 = qJD(4) * t373;
t320 = rSges(4,1) * t255 + rSges(4,2) * cos(pkin(9));
t318 = t127 * rSges(6,1) - t126 * rSges(6,2);
t308 = Icges(5,1) * t244 - t409;
t304 = -Icges(5,2) * t243 + t408;
t300 = Icges(5,5) * t244 - Icges(5,6) * t243;
t290 = t153 * t261 + t155 * t259;
t287 = -t390 - t393;
t284 = t330 - t426;
t33 = t202 * t120 - t122 * t386 + t203 * t124 + t128 * t172 + t129 * t176 + t174 * t269;
t34 = -t121 * t386 - t202 * t123 + t203 * t125 - t128 * t175 + t129 * t177 + t173 * t269;
t282 = -t19 / 0.2e1 - t21 / 0.2e1 - t33 / 0.2e1 - t34 / 0.2e1;
t31 = -t204 * t120 + t122 * t384 + t205 * t124 + t126 * t172 + t127 * t176 - t174 * t456;
t32 = t121 * t384 + t204 * t123 + t205 * t125 - t126 * t175 + t127 * t177 - t173 * t456;
t281 = t20 / 0.2e1 + t22 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1;
t280 = t59 / 0.2e1 + t61 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1;
t279 = -t60 / 0.2e1 - t62 / 0.2e1 - t88 / 0.2e1 - t89 / 0.2e1;
t278 = rSges(3,3) * t261 + t259 * t430;
t131 = (-rSges(6,1) * t258 - rSges(6,2) * t260) * t359 + (-t243 * t316 + t413) * qJD(4);
t277 = t259 * t131 + t179 * t365;
t275 = t232 + t283;
t274 = t287 * t259;
t273 = qJD(4) * t308;
t272 = qJD(4) * t304;
t271 = t261 * t344 - t425;
t270 = t261 * t343 - t425;
t267 = t321 + t347;
t266 = t309 + t348;
t264 = t259 * t377 + t261 * t326;
t262 = t259 * t356 + t261 * t320;
t209 = t259 * t216;
t201 = -rSges(3,2) * t261 + t259 * rSges(3,3) + t368;
t200 = t249 + t278;
t198 = t216 * t366;
t197 = pkin(8) * t384 - t233;
t188 = t261 * t197;
t187 = t412 - t276;
t171 = t247 + (t430 * t261 + (-rSges(3,3) - qJ(2)) * t259) * qJD(1);
t170 = qJD(1) * t278 + t369;
t169 = t259 * t320 + t261 * t411 + t368;
t168 = t249 + t262;
t161 = -pkin(8) * t339 + t348;
t158 = t261 * (qJD(1) * t232 - t347);
t157 = (-t179 - t216) * t261;
t156 = t179 * t259 + t209;
t135 = qJD(1) * t443 + t300 * t362;
t134 = -t300 * t360 + t367;
t133 = (t356 * t261 + (-qJ(2) - t320) * t259) * qJD(1) + t324;
t132 = qJD(1) * t262 + t345;
t108 = (-t216 - t373) * t261;
t107 = t209 + t327;
t104 = t259 * t343 + t275 + t370;
t103 = t270 + t317 + t447;
t102 = -t243 * t155 + t179 * t384;
t101 = t153 * t243 + t179 * t386;
t100 = -t259 * t443 - t261 * t287;
t99 = t259 * t180 - t450;
t98 = -t261 * t443 + t274;
t97 = t180 * t261 + t259 * t288;
t92 = t290 * t244;
t91 = t277 + t372;
t90 = t340 + t198 + (-t131 - t211) * t261;
t85 = t259 * t344 + t275 + t448;
t84 = t271 + t447 - t465;
t81 = -rSges(6,3) * t339 + t351;
t79 = -rSges(6,3) * t456 + t318;
t65 = t259 * t375 + t188 + t399;
t64 = -t243 * t374 + t373 * t384;
t63 = t243 * t376 + t244 * t327;
t50 = t264 + t372;
t49 = t198 + t259 * t326 + (-t211 - t377) * t261;
t48 = t440 * t244;
t47 = rSges(6,3) * t336 + (-t251 + (t284 + t413) * t259) * qJD(1) + t267 - t318;
t46 = qJD(1) * t270 + t266 + t351;
t45 = t259 * t350 + t188 + t449;
t44 = (-t179 * t362 + t81) * t243 + (qJD(4) * t153 + t277) * t244;
t43 = (-t179 * t360 - t79) * t243 + (-qJD(4) * t155 + t131 * t261 - t340) * t244;
t42 = rSges(7,2) * t336 + (-t251 + (t284 + t414) * t259) * qJD(1) + t267 + t466;
t41 = qJD(1) * t271 + t266 - t451;
t30 = t290 * t364 + (-t259 * t79 - t261 * t81 + (t259 * t153 - t399) * qJD(1)) * t244;
t29 = t261 * t79 + t158 + (-t161 - t81) * t259 + (t375 * t261 + (-t155 - t197) * t259) * qJD(1);
t24 = (-t259 * t325 - t416) * t243 + (qJD(4) * t376 + t264) * t244;
t23 = (-t261 * t325 - t417) * t243 + (-qJD(4) * t374 + t261 * t377 - t366 * t373) * t244;
t18 = -t128 * t147 + t129 * t151 + t143 * t269 - t202 * t72 + t203 * t76 - t386 * t68;
t17 = -t128 * t146 + t129 * t150 + t142 * t269 - t202 * t73 + t203 * t77 - t386 * t69;
t16 = t128 * t141 + t129 * t149 + t145 * t269 + t202 * t66 + t203 * t74 - t386 * t70;
t15 = t128 * t140 + t129 * t148 + t144 * t269 + t202 * t67 + t203 * t75 - t386 * t71;
t14 = -t126 * t147 + t127 * t151 - t143 * t456 + t204 * t72 + t205 * t76 + t384 * t68;
t13 = -t126 * t146 + t127 * t150 - t142 * t456 + t204 * t73 + t205 * t77 + t384 * t69;
t12 = t126 * t141 + t127 * t149 - t145 * t456 - t204 * t66 + t205 * t74 + t384 * t70;
t11 = t126 * t140 + t127 * t148 - t144 * t456 - t204 * t67 + t205 * t75 + t384 * t71;
t10 = t158 + t417 * t261 + (-t161 + t416) * t259 + (t350 * t261 + (-t197 - t374) * t259) * qJD(1);
t9 = t440 * t364 + (t416 * t261 - t417 * t259 + (t259 * t376 - t449) * qJD(1)) * t244;
t8 = -qJD(1) * t312 + t17 * t261 + t18 * t259;
t7 = -qJD(1) * t313 + t15 * t261 + t16 * t259;
t6 = -qJD(1) * t310 + t13 * t261 + t14 * t259;
t5 = -qJD(1) * t311 + t11 * t261 + t12 * t259;
t4 = (qJD(4) * t312 + t34) * t243 + (-qJD(1) * t36 + qJD(4) * t87 - t17 * t259 + t18 * t261) * t244;
t3 = (qJD(4) * t313 + t33) * t243 + (-qJD(1) * t35 + qJD(4) * t86 - t15 * t259 + t16 * t261) * t244;
t2 = (qJD(4) * t310 + t32) * t243 + (-qJD(1) * t38 + qJD(4) * t89 - t13 * t259 + t14 * t261) * t244;
t1 = (qJD(4) * t311 + t31) * t243 + (-qJD(1) * t37 + qJD(4) * t88 - t11 * t259 + t12 * t261) * t244;
t25 = [-t177 * t333 - t172 * t338 - t243 * t273 - t175 * t332 - t176 * t337 - t307 * t363 + 0.2e1 * m(3) * (t170 * t201 + t171 * t200) + (t105 * t160 + t106 * t159) * t439 + 0.2e1 * m(4) * (t132 * t169 + t133 * t168) + (t103 * t47 + t104 * t46) * t438 + (t41 * t85 + t42 * t84) * t437 + t303 * t364 + (-t272 - t383) * t244 + t458; m(7) * (t259 * t42 - t261 * t41 + (t259 * t85 + t261 * t84) * qJD(1)) + m(6) * (t259 * t47 - t261 * t46 + (t103 * t261 + t104 * t259) * qJD(1)) + m(5) * ((t159 * t261 + t160 * t259) * qJD(1) + t444) + m(4) * (-t132 * t261 + t259 * t133 + (t168 * t261 + t169 * t259) * qJD(1)) + m(3) * (-t170 * t261 + t259 * t171 + (t200 * t261 + t201 * t259) * qJD(1)); 0; m(7) * (t259 * t41 + t261 * t42 + (-t259 * t84 + t261 * t85) * qJD(1)) + m(6) * (t259 * t46 + t261 * t47 + (-t103 * t259 + t104 * t261) * qJD(1)) + m(5) * (qJD(1) * t445 + t259 * t105 + t106 * t261) + m(4) * (t259 * t132 + t133 * t261 + (-t168 * t259 + t169 * t261) * qJD(1)); 0; 0; ((qJD(1) * t442 + t259 * t272) * t436 + (qJD(1) * t441 + t259 * t273) * t434 + (-t395 / 0.2e1 - t392 / 0.2e1) * qJD(4) - t282) * t261 + ((qJD(1) * t182 - t304 * t360) * t436 + (qJD(1) * t184 - t308 * t360) * t434 + (t393 / 0.2e1 + t390 / 0.2e1) * qJD(4) + t281) * t259 + m(5) * (t210 * t445 + t215 * t444) + m(6) * (t103 * t91 + t104 * t90 + t156 * t47 + t157 * t46) + m(7) * (t107 * t42 + t108 * t41 + t49 * t85 + t50 * t84) - (t254 / 0.2e1 + t253 / 0.2e1) * t299 * qJD(4) + ((t160 * t428 + t396 / 0.2e1 - t391 / 0.2e1 - t280) * t259 + (t159 * t428 + t394 / 0.2e1 - t389 / 0.2e1 - t279) * t261) * qJD(1); m(6) * (t91 * t259 - t261 * t90 + (t156 * t261 + t157 * t259) * qJD(1)) + m(7) * (t50 * t259 - t261 * t49 + (t107 * t261 + t108 * t259) * qJD(1)) - m(5) * t328; m(6) * (t90 * t259 + t261 * t91 + (-t156 * t259 + t157 * t261) * qJD(1)) + m(7) * (t49 * t259 + t261 * t50 + (-t107 * t259 + t108 * t261) * qJD(1)); (t10 * t45 + t107 * t50 + t108 * t49) * t437 + (t156 * t91 + t157 * t90 + t29 * t65) * t438 + t259 * t5 + t259 * t6 + t261 * t8 + t261 * t7 + t261 * ((t261 * t135 + (t98 + t450) * qJD(1)) * t261 + (-t97 * qJD(1) + (-t363 * t441 + t364 * t442) * t259 + (t134 + (t391 - t396) * qJD(4) + (-t180 + t287) * qJD(1)) * t261) * t259) + t259 * ((t259 * t134 + (-t99 + t274) * qJD(1)) * t259 + (t100 * qJD(1) + (t182 * t364 - t184 * t363 + t367) * t261 + (t135 + (t389 - t394) * qJD(4) + t288 * qJD(1)) * t259) * t261) + ((-t259 * t186 + t187 * t261) * (-t259 * t349 + (-t215 * t254 + t253 * t415) * qJD(4) + ((-t186 + t250) * t261 + (-t187 + t276 + t412) * t259) * qJD(1)) - t215 * t328) * t439 + (-t98 * t259 - t261 * t97 - t35 - t36) * t366 + (t100 * t259 + t261 * t99 + t37 + t38) * t365; m(6) * (t101 * t46 + t102 * t47 + t103 * t43 + t104 * t44) + m(7) * (t23 * t84 + t24 * t85 + t41 * t63 + t42 * t64) + (t259 * t280 + t261 * t279) * t364 + (t281 * t261 + t282 * t259 + (t259 * t279 - t261 * t280) * qJD(1)) * t244 + t329; m(6) * (t43 * t259 - t261 * t44 + (t101 * t259 + t102 * t261) * qJD(1)) + m(7) * (t23 * t259 - t24 * t261 + (t259 * t63 + t261 * t64) * qJD(1)); m(6) * (t44 * t259 + t261 * t43 + (t101 * t261 - t102 * t259) * qJD(1)) + m(7) * (t23 * t261 + t24 * t259 + (-t259 * t64 + t261 * t63) * qJD(1)); m(6) * (t101 * t90 + t102 * t91 + t156 * t43 + t157 * t44 - t29 * t92 + t30 * t65) + m(7) * (-t10 * t48 + t107 * t23 + t108 * t24 + t45 * t9 + t49 * t63 + t50 * t64) + (t3 / 0.2e1 + t4 / 0.2e1 + t354 * t364 + t460 * qJD(1) / 0.2e1) * t261 + (t1 / 0.2e1 + t2 / 0.2e1 + t355 * t364 - t461 * qJD(1) / 0.2e1) * t259 + ((t259 * t354 - t261 * t355) * qJD(1) - (t7 + t8) * t259 / 0.2e1 + (t5 + t6) * t261 / 0.2e1 + t454 * qJD(4) / 0.2e1) * t244 + (qJD(1) * t455 + t452 * t259 + t453 * t261) * t243 / 0.2e1; (t23 * t64 + t24 * t63 - t48 * t9) * t437 + (t101 * t44 + t102 * t43 - t30 * t92) * t438 + (((-t243 * t418 - t460) * t261 + (t243 * t419 + t461) * t259) * qJD(4) + t329) * t243 + ((t1 + t2) * t261 + (-t3 - t4) * t259 + (-t259 * t453 + t452 * t261) * t243 + (t243 * t459 + t244 * t455) * qJD(4) + (-t243 * t454 - t259 * t460 - t261 * t461) * qJD(1)) * t244; m(7) * (t126 * t85 + t128 * t84 + t202 * t42 - t204 * t41); m(7) * (-t126 * t261 + t128 * t259 + (t202 * t261 - t204 * t259) * qJD(1)); m(7) * (t126 * t259 + t128 * t261 + (-t202 * t259 - t204 * t261) * qJD(1)); m(7) * (-t45 * t338 + t107 * t128 + t108 * t126 + t202 * t50 - t204 * t49 + (t10 * t258 + t358 * t45) * t244); m(7) * (t48 * t338 + t126 * t63 + t128 * t64 + t202 * t23 - t204 * t24 + (t258 * t9 - t358 * t48) * t244); (-t126 * t204 + t128 * t202 + (t332 - t338) * t387) * t437;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
