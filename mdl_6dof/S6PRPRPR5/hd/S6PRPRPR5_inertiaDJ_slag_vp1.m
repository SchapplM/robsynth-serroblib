% Calculate time derivative of joint inertia matrix for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:43:01
% EndTime: 2019-03-08 19:43:31
% DurationCPUTime: 21.99s
% Computational Cost: add. (53801->1058), mult. (99764->1513), div. (0->0), fcn. (111628->12), ass. (0->403)
t355 = sin(pkin(10));
t358 = cos(pkin(10));
t364 = cos(qJ(2));
t359 = cos(pkin(6));
t362 = sin(qJ(2));
t433 = t359 * t362;
t341 = t355 * t364 + t358 * t433;
t356 = sin(pkin(6));
t408 = pkin(11) + qJ(4);
t392 = cos(t408);
t383 = t356 * t392;
t391 = sin(t408);
t307 = t341 * t391 + t358 * t383;
t382 = t356 * t391;
t308 = t341 * t392 - t358 * t382;
t432 = t359 * t364;
t368 = -t355 * t362 + t358 * t432;
t196 = -Icges(6,5) * t368 - Icges(6,6) * t308 + Icges(6,3) * t307;
t202 = Icges(5,4) * t308 - Icges(5,2) * t307 - Icges(5,6) * t368;
t471 = t196 - t202;
t404 = t355 * t433;
t343 = t358 * t364 - t404;
t309 = t343 * t391 - t355 * t383;
t310 = t343 * t392 + t355 * t382;
t342 = t355 * t432 + t358 * t362;
t197 = Icges(6,5) * t342 - Icges(6,6) * t310 + Icges(6,3) * t309;
t203 = Icges(5,4) * t310 - Icges(5,2) * t309 + Icges(5,6) * t342;
t470 = t197 - t203;
t198 = Icges(5,5) * t308 - Icges(5,6) * t307 - Icges(5,3) * t368;
t204 = -Icges(6,1) * t368 - Icges(6,4) * t308 + Icges(6,5) * t307;
t469 = t198 + t204;
t199 = Icges(5,5) * t310 - Icges(5,6) * t309 + Icges(5,3) * t342;
t205 = Icges(6,1) * t342 - Icges(6,4) * t310 + Icges(6,5) * t309;
t468 = t199 + t205;
t200 = -Icges(6,4) * t368 - Icges(6,2) * t308 + Icges(6,6) * t307;
t206 = Icges(5,1) * t308 - Icges(5,4) * t307 - Icges(5,5) * t368;
t467 = -t200 + t206;
t201 = Icges(6,4) * t342 - Icges(6,2) * t310 + Icges(6,6) * t309;
t207 = Icges(5,1) * t310 - Icges(5,4) * t309 + Icges(5,5) * t342;
t466 = -t201 + t207;
t463 = t307 * t471 + t467 * t308 - t469 * t368;
t462 = t307 * t470 + t308 * t466 - t368 * t468;
t461 = t309 * t471 + t467 * t310 + t469 * t342;
t460 = t309 * t470 + t310 * t466 + t342 * t468;
t327 = -t359 * t392 + t362 * t382;
t328 = t359 * t391 + t362 * t383;
t434 = t356 * t364;
t270 = -Icges(6,5) * t434 - Icges(6,6) * t328 + Icges(6,3) * t327;
t271 = -Icges(6,4) * t434 - Icges(6,2) * t328 + Icges(6,6) * t327;
t272 = -Icges(6,1) * t434 - Icges(6,4) * t328 + Icges(6,5) * t327;
t130 = t270 * t307 - t271 * t308 - t272 * t368;
t273 = Icges(5,5) * t328 - Icges(5,6) * t327 - Icges(5,3) * t434;
t274 = Icges(5,4) * t328 - Icges(5,2) * t327 - Icges(5,6) * t434;
t275 = Icges(5,1) * t328 - Icges(5,4) * t327 - Icges(5,5) * t434;
t132 = -t273 * t368 - t274 * t307 + t275 * t308;
t465 = t130 + t132;
t131 = t270 * t309 - t271 * t310 + t272 * t342;
t133 = t273 * t342 - t274 * t309 + t275 * t310;
t464 = t131 + t133;
t407 = m(6) / 0.2e1 + m(7) / 0.2e1;
t459 = 0.2e1 * t407;
t410 = qJD(2) * t362;
t393 = t356 * t410;
t119 = t327 * t196 - t328 * t200 - t204 * t434;
t121 = -t198 * t434 - t327 * t202 + t328 * t206;
t458 = t119 + t121;
t120 = t327 * t197 - t328 * t201 - t205 * t434;
t122 = -t199 * t434 - t327 * t203 + t328 * t207;
t457 = t120 + t122;
t330 = t368 * qJD(2);
t380 = qJD(4) * t391;
t370 = t356 * t380;
t381 = qJD(4) * t392;
t258 = t330 * t391 + t341 * t381 - t358 * t370;
t361 = sin(qJ(6));
t363 = cos(qJ(6));
t265 = t307 * t361 - t363 * t368;
t331 = t341 * qJD(2);
t178 = -qJD(6) * t265 + t258 * t363 - t331 * t361;
t264 = t307 * t363 + t361 * t368;
t179 = qJD(6) * t264 + t258 * t361 + t331 * t363;
t371 = t356 * t381;
t259 = t330 * t392 - t341 * t380 - t358 * t371;
t100 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t259;
t102 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t259;
t104 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t259;
t154 = Icges(7,5) * t265 + Icges(7,6) * t264 + Icges(7,3) * t308;
t156 = Icges(7,4) * t265 + Icges(7,2) * t264 + Icges(7,6) * t308;
t158 = Icges(7,1) * t265 + Icges(7,4) * t264 + Icges(7,5) * t308;
t22 = t100 * t308 + t102 * t264 + t104 * t265 + t154 * t259 + t156 * t178 + t158 * t179;
t332 = t342 * qJD(2);
t260 = qJD(4) * t310 - t332 * t391;
t267 = t309 * t361 + t342 * t363;
t409 = qJD(2) * t364;
t333 = -qJD(2) * t404 + t358 * t409;
t180 = -qJD(6) * t267 + t260 * t363 - t333 * t361;
t266 = t309 * t363 - t342 * t361;
t181 = qJD(6) * t266 + t260 * t361 + t333 * t363;
t261 = -t332 * t392 - t343 * t380 + t355 * t371;
t101 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t261;
t103 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t261;
t105 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t261;
t155 = Icges(7,5) * t267 + Icges(7,6) * t266 + Icges(7,3) * t310;
t157 = Icges(7,4) * t267 + Icges(7,2) * t266 + Icges(7,6) * t310;
t159 = Icges(7,1) * t267 + Icges(7,4) * t266 + Icges(7,5) * t310;
t23 = t101 * t308 + t103 * t264 + t105 * t265 + t155 * t259 + t157 * t178 + t159 * t179;
t305 = qJD(4) * t328 + t382 * t409;
t369 = -t327 * t361 + t363 * t434;
t217 = qJD(6) * t369 + t305 * t363 - t361 * t393;
t314 = t327 * t363 + t361 * t434;
t218 = qJD(6) * t314 + t305 * t361 + t363 * t393;
t306 = t359 * t381 - t362 * t370 + t383 * t409;
t135 = Icges(7,5) * t218 + Icges(7,6) * t217 + Icges(7,3) * t306;
t136 = Icges(7,4) * t218 + Icges(7,2) * t217 + Icges(7,6) * t306;
t137 = Icges(7,1) * t218 + Icges(7,4) * t217 + Icges(7,5) * t306;
t190 = -Icges(7,5) * t369 + Icges(7,6) * t314 + Icges(7,3) * t328;
t191 = -Icges(7,4) * t369 + Icges(7,2) * t314 + Icges(7,6) * t328;
t192 = -Icges(7,1) * t369 + Icges(7,4) * t314 + Icges(7,5) * t328;
t41 = t135 * t308 + t136 * t264 + t137 * t265 + t178 * t191 + t179 * t192 + t190 * t259;
t77 = t154 * t308 + t156 * t264 + t158 * t265;
t78 = t155 * t308 + t157 * t264 + t159 * t265;
t93 = t190 * t308 + t191 * t264 + t192 * t265;
t3 = -t22 * t368 + t23 * t342 + t77 * t331 + t78 * t333 + (-t364 * t41 + t410 * t93) * t356;
t164 = Icges(5,5) * t259 - Icges(5,6) * t258 + Icges(5,3) * t331;
t168 = Icges(5,4) * t259 - Icges(5,2) * t258 + Icges(5,6) * t331;
t172 = Icges(5,1) * t259 - Icges(5,4) * t258 + Icges(5,5) * t331;
t51 = -t164 * t368 - t168 * t307 + t172 * t308 + t198 * t331 - t202 * t258 + t206 * t259;
t165 = Icges(5,5) * t261 - Icges(5,6) * t260 + Icges(5,3) * t333;
t169 = Icges(5,4) * t261 - Icges(5,2) * t260 + Icges(5,6) * t333;
t173 = Icges(5,1) * t261 - Icges(5,4) * t260 + Icges(5,5) * t333;
t52 = -t165 * t368 - t169 * t307 + t173 * t308 + t199 * t331 - t203 * t258 + t207 * t259;
t162 = Icges(6,5) * t331 - Icges(6,6) * t259 + Icges(6,3) * t258;
t166 = Icges(6,4) * t331 - Icges(6,2) * t259 + Icges(6,6) * t258;
t170 = Icges(6,1) * t331 - Icges(6,4) * t259 + Icges(6,5) * t258;
t55 = t162 * t307 - t166 * t308 - t170 * t368 + t196 * t258 - t200 * t259 + t204 * t331;
t163 = Icges(6,5) * t333 - Icges(6,6) * t261 + Icges(6,3) * t260;
t167 = Icges(6,4) * t333 - Icges(6,2) * t261 + Icges(6,6) * t260;
t171 = Icges(6,1) * t333 - Icges(6,4) * t261 + Icges(6,5) * t260;
t56 = t163 * t307 - t167 * t308 - t171 * t368 + t197 * t258 - t201 * t259 + t205 * t331;
t222 = Icges(5,5) * t306 - Icges(5,6) * t305 + Icges(5,3) * t393;
t223 = Icges(5,4) * t306 - Icges(5,2) * t305 + Icges(5,6) * t393;
t224 = Icges(5,1) * t306 - Icges(5,4) * t305 + Icges(5,5) * t393;
t71 = -t222 * t368 - t223 * t307 + t224 * t308 - t258 * t274 + t259 * t275 + t273 * t331;
t219 = Icges(6,5) * t393 - Icges(6,6) * t306 + Icges(6,3) * t305;
t220 = Icges(6,4) * t393 - Icges(6,2) * t306 + Icges(6,6) * t305;
t221 = Icges(6,1) * t393 - Icges(6,4) * t306 + Icges(6,5) * t305;
t73 = t219 * t307 - t220 * t308 - t221 * t368 + t258 * t270 - t259 * t271 + t272 * t331;
t456 = t3 + (-t51 - t55) * t368 + (t465 * t410 + (-t71 - t73) * t364) * t356 + (t52 + t56) * t342 + t462 * t333 + t463 * t331;
t24 = t100 * t310 + t102 * t266 + t104 * t267 + t154 * t261 + t156 * t180 + t158 * t181;
t25 = t101 * t310 + t103 * t266 + t105 * t267 + t155 * t261 + t157 * t180 + t159 * t181;
t42 = t135 * t310 + t136 * t266 + t137 * t267 + t180 * t191 + t181 * t192 + t190 * t261;
t79 = t154 * t310 + t156 * t266 + t158 * t267;
t80 = t155 * t310 + t157 * t266 + t159 * t267;
t94 = t190 * t310 + t191 * t266 + t192 * t267;
t4 = -t24 * t368 + t25 * t342 + t79 * t331 + t80 * t333 + (-t364 * t42 + t410 * t94) * t356;
t53 = t164 * t342 - t168 * t309 + t172 * t310 + t198 * t333 - t202 * t260 + t206 * t261;
t54 = t165 * t342 - t169 * t309 + t173 * t310 + t199 * t333 - t203 * t260 + t207 * t261;
t57 = t162 * t309 - t166 * t310 + t170 * t342 + t196 * t260 - t200 * t261 + t204 * t333;
t58 = t163 * t309 - t167 * t310 + t171 * t342 + t197 * t260 - t201 * t261 + t205 * t333;
t72 = t222 * t342 - t223 * t309 + t224 * t310 - t260 * t274 + t261 * t275 + t273 * t333;
t74 = t219 * t309 - t220 * t310 + t221 * t342 + t260 * t270 - t261 * t271 + t272 * t333;
t455 = t4 + (-t53 - t57) * t368 + (t464 * t410 + (-t72 - t74) * t364) * t356 + (t54 + t58) * t342 + t460 * t333 + t461 * t331;
t454 = 2 * m(5);
t453 = 0.2e1 * m(6);
t452 = 0.2e1 * m(7);
t451 = t259 / 0.2e1;
t450 = t261 / 0.2e1;
t449 = t306 / 0.2e1;
t448 = t308 / 0.2e1;
t447 = t310 / 0.2e1;
t446 = t328 / 0.2e1;
t445 = t331 / 0.2e1;
t444 = t333 / 0.2e1;
t443 = t355 / 0.2e1;
t442 = -t358 / 0.2e1;
t357 = cos(pkin(11));
t441 = pkin(3) * t357;
t440 = Icges(3,4) * t362;
t439 = Icges(3,4) * t364;
t354 = sin(pkin(11));
t438 = t354 * t359;
t437 = t355 * t356;
t436 = t356 * t358;
t435 = t356 * t362;
t106 = rSges(7,1) * t179 + rSges(7,2) * t178 + rSges(7,3) * t259;
t430 = pkin(5) * t331 + pkin(9) * t259 + t106;
t107 = rSges(7,1) * t181 + rSges(7,2) * t180 + rSges(7,3) * t261;
t429 = pkin(5) * t333 + pkin(9) * t261 + t107;
t138 = rSges(7,1) * t218 + rSges(7,2) * t217 + rSges(7,3) * t306;
t428 = pkin(5) * t393 + pkin(9) * t306 + t138;
t152 = pkin(4) * t259 + qJ(5) * t258 + qJD(5) * t307;
t243 = pkin(4) * t308 + qJ(5) * t307;
t427 = t342 * t152 + t333 * t243;
t153 = pkin(4) * t261 + qJ(5) * t260 + qJD(5) * t309;
t177 = rSges(6,1) * t333 - rSges(6,2) * t261 + rSges(6,3) * t260;
t426 = -t153 - t177;
t160 = rSges(7,1) * t265 + rSges(7,2) * t264 + rSges(7,3) * t308;
t425 = -pkin(5) * t368 + pkin(9) * t308 + t160;
t161 = rSges(7,1) * t267 + rSges(7,2) * t266 + rSges(7,3) * t310;
t424 = pkin(5) * t342 + pkin(9) * t310 + t161;
t209 = pkin(8) * t333 - t332 * t441;
t263 = -pkin(2) * t332 + qJ(3) * t333 + qJD(3) * t342;
t256 = t359 * t263;
t423 = t359 * t209 + t256;
t406 = t354 * t437;
t215 = pkin(3) * t406 + pkin(8) * t342 + t343 * t441;
t304 = pkin(2) * t343 + qJ(3) * t342;
t302 = t359 * t304;
t422 = t359 * t215 + t302;
t194 = -rSges(7,1) * t369 + rSges(7,2) * t314 + rSges(7,3) * t328;
t421 = -pkin(5) * t434 + t328 * pkin(9) + t194;
t208 = pkin(8) * t331 + t330 * t441;
t262 = pkin(2) * t330 + qJ(3) * t331 - qJD(3) * t368;
t420 = -t208 - t262;
t211 = rSges(6,1) * t342 - rSges(6,2) * t310 + rSges(6,3) * t309;
t244 = pkin(4) * t310 + qJ(5) * t309;
t419 = -t211 - t244;
t405 = t354 * t436;
t214 = -pkin(3) * t405 - pkin(8) * t368 + t341 * t441;
t303 = pkin(2) * t341 - qJ(3) * t368;
t418 = -t214 - t303;
t282 = pkin(4) * t328 + qJ(5) * t327;
t417 = t243 * t434 - t282 * t368;
t416 = t262 * t437 + t263 * t436;
t276 = -rSges(6,1) * t434 - t328 * rSges(6,2) + t327 * rSges(6,3);
t415 = -t276 - t282;
t344 = (pkin(2) * t362 - qJ(3) * t364) * t356;
t414 = -pkin(3) * t438 - (-pkin(8) * t364 + t362 * t441) * t356 - t344;
t413 = t303 * t437 + t304 * t436;
t322 = (-qJD(3) * t364 + (pkin(2) * t364 + qJ(3) * t362) * qJD(2)) * t356;
t411 = qJD(2) * t356;
t412 = -(pkin(8) * t362 + t364 * t441) * t411 - t322;
t403 = -t153 - t429;
t184 = pkin(4) * t306 + qJ(5) * t305 + qJD(5) * t327;
t402 = t152 * t434 - t184 * t368 + t331 * t282;
t401 = t359 * t153 + t423;
t400 = -t152 + t420;
t399 = -t244 - t424;
t398 = -t184 + t412;
t397 = t359 * t244 + t422;
t396 = -t282 - t421;
t395 = -t243 + t418;
t394 = -t282 + t414;
t338 = -t354 * t435 + t357 * t359;
t339 = t357 * t435 + t438;
t390 = (-t339 * rSges(4,1) - t338 * rSges(4,2) + rSges(4,3) * t434 - t344) * t356;
t384 = rSges(4,1) * t357 - rSges(4,2) * t354;
t389 = (-(rSges(4,3) * t362 + t364 * t384) * t411 - t322) * t356;
t388 = t208 * t437 + t209 * t436 + t416;
t387 = t214 * t437 + t215 * t436 + t413;
t226 = rSges(5,1) * t306 - rSges(5,2) * t305 + rSges(5,3) * t393;
t386 = (-t226 + t412) * t356;
t277 = t328 * rSges(5,1) - t327 * rSges(5,2) - rSges(5,3) * t434;
t385 = (-t277 + t414) * t356;
t379 = Icges(4,1) * t357 - Icges(4,4) * t354;
t378 = Icges(4,4) * t357 - Icges(4,2) * t354;
t377 = Icges(4,5) * t357 - Icges(4,6) * t354;
t376 = -(Icges(4,4) * t339 + Icges(4,2) * t338 - Icges(4,6) * t434) * t354 + (Icges(4,1) * t339 + Icges(4,4) * t338 - Icges(4,5) * t434) * t357;
t225 = rSges(6,1) * t393 - rSges(6,2) * t306 + rSges(6,3) * t305;
t375 = (-t225 + t398) * t356;
t374 = (-t276 + t394) * t356;
t373 = t152 * t437 + t153 * t436 + t388;
t372 = t243 * t437 + t244 * t436 + t387;
t367 = (t398 - t428) * t356;
t366 = (t394 - t421) * t356;
t317 = -t341 * t354 - t357 * t436;
t318 = t341 * t357 - t405;
t319 = -t343 * t354 + t357 * t437;
t320 = t343 * t357 + t406;
t365 = (-(Icges(4,4) * t320 + Icges(4,2) * t319 + Icges(4,6) * t342) * t354 + (Icges(4,1) * t320 + Icges(4,4) * t319 + Icges(4,5) * t342) * t357) * t355 - (-(Icges(4,4) * t318 + Icges(4,2) * t317 - Icges(4,6) * t368) * t354 + (Icges(4,1) * t318 + Icges(4,4) * t317 - Icges(4,5) * t368) * t357) * t358;
t337 = (rSges(3,1) * t364 - rSges(3,2) * t362) * t411;
t336 = (Icges(3,1) * t364 - t440) * t411;
t335 = (-Icges(3,2) * t362 + t439) * t411;
t334 = (Icges(3,5) * t364 - Icges(3,6) * t362) * t411;
t329 = t359 * rSges(3,3) + (rSges(3,1) * t362 + rSges(3,2) * t364) * t356;
t326 = Icges(3,5) * t359 + (Icges(3,1) * t362 + t439) * t356;
t325 = Icges(3,6) * t359 + (Icges(3,2) * t364 + t440) * t356;
t313 = (Icges(4,5) * t362 + t364 * t379) * t411;
t312 = (Icges(4,6) * t362 + t364 * t378) * t411;
t311 = (Icges(4,3) * t362 + t364 * t377) * t411;
t301 = -rSges(3,1) * t332 - rSges(3,2) * t333;
t300 = rSges(3,1) * t330 - rSges(3,2) * t331;
t298 = -Icges(3,1) * t332 - Icges(3,4) * t333;
t297 = Icges(3,1) * t330 - Icges(3,4) * t331;
t296 = -Icges(3,4) * t332 - Icges(3,2) * t333;
t295 = Icges(3,4) * t330 - Icges(3,2) * t331;
t294 = -Icges(3,5) * t332 - Icges(3,6) * t333;
t293 = Icges(3,5) * t330 - Icges(3,6) * t331;
t290 = rSges(3,1) * t343 - rSges(3,2) * t342 + rSges(3,3) * t437;
t289 = rSges(3,1) * t341 + rSges(3,2) * t368 - rSges(3,3) * t436;
t286 = Icges(3,1) * t343 - Icges(3,4) * t342 + Icges(3,5) * t437;
t285 = Icges(3,1) * t341 + Icges(3,4) * t368 - Icges(3,5) * t436;
t284 = Icges(3,4) * t343 - Icges(3,2) * t342 + Icges(3,6) * t437;
t283 = Icges(3,4) * t341 + Icges(3,2) * t368 - Icges(3,6) * t436;
t279 = Icges(4,5) * t339 + Icges(4,6) * t338 - Icges(4,3) * t434;
t252 = rSges(4,3) * t333 - t332 * t384;
t251 = rSges(4,3) * t331 + t330 * t384;
t250 = Icges(4,5) * t333 - t332 * t379;
t249 = Icges(4,5) * t331 + t330 * t379;
t248 = Icges(4,6) * t333 - t332 * t378;
t247 = Icges(4,6) * t331 + t330 * t378;
t246 = Icges(4,3) * t333 - t332 * t377;
t245 = Icges(4,3) * t331 + t330 * t377;
t238 = t244 * t393;
t237 = rSges(4,1) * t320 + rSges(4,2) * t319 + rSges(4,3) * t342;
t236 = rSges(4,1) * t318 + rSges(4,2) * t317 - rSges(4,3) * t368;
t230 = Icges(4,5) * t320 + Icges(4,6) * t319 + Icges(4,3) * t342;
t229 = Icges(4,5) * t318 + Icges(4,6) * t317 - Icges(4,3) * t368;
t216 = t342 * t243;
t213 = rSges(5,1) * t310 - rSges(5,2) * t309 + rSges(5,3) * t342;
t212 = rSges(5,1) * t308 - rSges(5,2) * t307 - rSges(5,3) * t368;
t210 = -rSges(6,1) * t368 - rSges(6,2) * t308 + rSges(6,3) * t307;
t183 = (t300 * t355 + t301 * t358) * t356;
t176 = rSges(6,1) * t331 - rSges(6,2) * t259 + rSges(6,3) * t258;
t175 = rSges(5,1) * t261 - rSges(5,2) * t260 + rSges(5,3) * t333;
t174 = rSges(5,1) * t259 - rSges(5,2) * t258 + rSges(5,3) * t331;
t151 = -t213 * t434 - t342 * t277;
t150 = t212 * t434 - t277 * t368;
t144 = (-t236 - t303) * t359 + t358 * t390;
t143 = t237 * t359 + t355 * t390 + t302;
t142 = -t273 * t434 - t327 * t274 + t328 * t275;
t141 = t327 * t270 - t328 * t271 - t272 * t434;
t140 = (-t251 - t262) * t359 + t358 * t389;
t139 = t252 * t359 + t355 * t389 + t256;
t134 = t212 * t342 + t213 * t368;
t128 = (t236 * t355 + t237 * t358) * t356 + t413;
t127 = (t251 * t355 + t252 * t358) * t356 + t416;
t126 = t161 * t328 - t194 * t310;
t125 = -t160 * t328 + t194 * t308;
t124 = t342 * t415 + t419 * t434;
t123 = t210 * t434 - t276 * t368 + t417;
t109 = (-t212 + t418) * t359 + t358 * t385;
t108 = t213 * t359 + t355 * t385 + t422;
t99 = t190 * t328 + t191 * t314 - t192 * t369;
t98 = t160 * t310 - t161 * t308;
t97 = t210 * t342 - t368 * t419 + t216;
t96 = -t342 * t226 - t333 * t277 + (-t175 * t364 + t213 * t410) * t356;
t95 = -t368 * t226 + t331 * t277 + (t174 * t364 - t212 * t410) * t356;
t92 = (t212 * t355 + t213 * t358) * t356 + t387;
t91 = (-t174 + t420) * t359 + t358 * t386;
t90 = t175 * t359 + t355 * t386 + t423;
t89 = (-t210 + t395) * t359 + t358 * t374;
t88 = t211 * t359 + t355 * t374 + t397;
t87 = t342 * t396 + t399 * t434;
t86 = -t368 * t421 + t425 * t434 + t417;
t85 = -t327 * t223 + t328 * t224 - t305 * t274 + t306 * t275 + (-t222 * t364 + t273 * t410) * t356;
t84 = t327 * t219 - t328 * t220 + t305 * t270 - t306 * t271 + (-t221 * t364 + t272 * t410) * t356;
t83 = t174 * t342 + t175 * t368 + t212 * t333 - t213 * t331;
t82 = t155 * t328 + t157 * t314 - t159 * t369;
t81 = t154 * t328 + t156 * t314 - t158 * t369;
t76 = (t210 * t355 + t211 * t358) * t356 + t372;
t75 = (t174 * t355 + t175 * t358) * t356 + t388;
t70 = t342 * t425 - t368 * t399 + t216;
t69 = (t395 - t425) * t359 + t358 * t366;
t68 = t355 * t366 + t359 * t424 + t397;
t67 = (-t176 + t400) * t359 + t358 * t375;
t66 = t177 * t359 + t355 * t375 + t401;
t65 = -t327 * t169 + t328 * t173 - t305 * t203 + t306 * t207 + (-t165 * t364 + t199 * t410) * t356;
t64 = -t327 * t168 + t328 * t172 - t305 * t202 + t306 * t206 + (-t164 * t364 + t198 * t410) * t356;
t63 = t327 * t163 - t328 * t167 + t305 * t197 - t306 * t201 + (-t171 * t364 + t205 * t410) * t356;
t62 = t327 * t162 - t328 * t166 + t305 * t196 - t306 * t200 + (-t170 * t364 + t204 * t410) * t356;
t61 = t238 + (-t184 - t225) * t342 + t415 * t333 + (t211 * t410 + t364 * t426) * t356;
t60 = -t368 * t225 + t331 * t276 + (t176 * t364 + (-t210 - t243) * t410) * t356 + t402;
t59 = (t355 * t425 + t358 * t424) * t356 + t372;
t50 = (t176 * t355 + t177 * t358) * t356 + t373;
t49 = t107 * t328 - t138 * t310 + t161 * t306 - t194 * t261;
t48 = -t106 * t328 + t138 * t308 - t160 * t306 + t194 * t259;
t47 = t176 * t342 + t210 * t333 + t331 * t419 - t368 * t426 + t427;
t46 = t135 * t328 + t136 * t314 - t137 * t369 + t190 * t306 + t191 * t217 + t192 * t218;
t45 = t106 * t310 - t107 * t308 + t160 * t261 - t161 * t259;
t44 = (t400 - t430) * t359 + t358 * t367;
t43 = t355 * t367 + t359 * t429 + t401;
t40 = t359 * t99 + (t355 * t82 - t358 * t81) * t356;
t39 = t82 * t342 - t368 * t81 - t434 * t99;
t38 = t238 + (-t184 - t428) * t342 + t396 * t333 + (t364 * t403 + t410 * t424) * t356;
t37 = -t428 * t368 + t421 * t331 + (t430 * t364 + (-t243 - t425) * t410) * t356 + t402;
t36 = t308 * t81 + t310 * t82 + t328 * t99;
t35 = (t355 * t430 + t358 * t429) * t356 + t373;
t34 = t359 * t94 + (t355 * t80 - t358 * t79) * t356;
t33 = t359 * t93 + (t355 * t78 - t358 * t77) * t356;
t32 = t80 * t342 - t368 * t79 - t434 * t94;
t31 = t78 * t342 - t368 * t77 - t434 * t93;
t30 = t308 * t79 + t310 * t80 + t328 * t94;
t29 = t308 * t77 + t310 * t78 + t328 * t93;
t28 = t101 * t328 + t103 * t314 - t105 * t369 + t155 * t306 + t157 * t217 + t159 * t218;
t27 = t100 * t328 + t102 * t314 - t104 * t369 + t154 * t306 + t156 * t217 + t158 * t218;
t26 = t331 * t399 + t333 * t425 + t342 * t430 - t368 * t403 + t427;
t21 = t359 * t85 + (t355 * t65 - t358 * t64) * t356;
t20 = t359 * t84 + (t355 * t63 - t358 * t62) * t356;
t19 = t359 * t74 + (t355 * t58 - t358 * t57) * t356;
t18 = t359 * t73 + (t355 * t56 - t358 * t55) * t356;
t17 = t359 * t72 + (t355 * t54 - t358 * t53) * t356;
t16 = t359 * t71 + (t355 * t52 - t358 * t51) * t356;
t15 = t121 * t331 + t122 * t333 - t64 * t368 + t65 * t342 + (t142 * t410 - t364 * t85) * t356;
t14 = t119 * t331 + t120 * t333 - t62 * t368 + t63 * t342 + (t141 * t410 - t364 * t84) * t356;
t9 = t359 * t46 + (-t27 * t358 + t28 * t355) * t356;
t8 = t359 * t42 + (-t24 * t358 + t25 * t355) * t356;
t7 = t359 * t41 + (-t22 * t358 + t23 * t355) * t356;
t6 = -t27 * t368 + t28 * t342 + t81 * t331 + t82 * t333 + (-t364 * t46 + t410 * t99) * t356;
t5 = t259 * t81 + t261 * t82 + t27 * t308 + t28 * t310 + t306 * t99 + t328 * t46;
t2 = t24 * t308 + t25 * t310 + t259 * t79 + t261 * t80 + t306 * t94 + t328 * t42;
t1 = t22 * t308 + t23 * t310 + t259 * t77 + t261 * t78 + t306 * t93 + t328 * t41;
t10 = [0; m(3) * t183 + m(4) * t127 + m(5) * t75 + m(6) * t50 + m(7) * t35; t8 * t437 - t7 * t436 + t17 * t437 + t19 * t437 - t16 * t436 - t18 * t436 - ((-t284 * t331 + t286 * t330 - t294 * t436 + t296 * t368 + t298 * t341) * t437 - (-t283 * t331 + t285 * t330 - t293 * t436 + t295 * t368 + t297 * t341) * t436 + (-t325 * t331 + t326 * t330 - t334 * t436 + t335 * t368 + t336 * t341) * t359) * t436 - ((t331 * t279 - t311 * t368 + t317 * t312 + t318 * t313 + t330 * t376) * t359 + ((t230 * t331 - t246 * t368 + t248 * t317 + t250 * t318) * t355 - (t229 * t331 - t245 * t368 + t247 * t317 + t249 * t318) * t358 + t365 * t330) * t356) * t436 + ((-t284 * t333 - t286 * t332 + t294 * t437 - t296 * t342 + t298 * t343) * t437 - (-t283 * t333 - t285 * t332 + t293 * t437 - t295 * t342 + t297 * t343) * t436 + (-t325 * t333 - t326 * t332 + t334 * t437 - t335 * t342 + t336 * t343) * t359) * t437 + ((t333 * t279 + t342 * t311 + t319 * t312 + t320 * t313 - t332 * t376) * t359 + ((t230 * t333 + t246 * t342 + t248 * t319 + t250 * t320) * t355 - (t229 * t333 + t245 * t342 + t247 * t319 + t249 * t320) * t358 - t365 * t332) * t356) * t437 + t359 * t9 + (t35 * t59 + t43 * t68 + t44 * t69) * t452 + (t50 * t76 + t66 * t88 + t67 * t89) * t453 + t359 * t21 + (t108 * t90 + t109 * t91 + t75 * t92) * t454 + t359 * t20 + 0.2e1 * m(4) * (t127 * t128 + t139 * t143 + t140 * t144) + t359 * (t359 ^ 2 * t334 + (((t296 * t364 + t298 * t362) * t355 - (t295 * t364 + t297 * t362) * t358 + ((-t284 * t362 + t286 * t364) * t355 - (-t283 * t362 + t285 * t364) * t358) * qJD(2)) * t356 + (-t293 * t358 + t294 * t355 + t335 * t364 + t336 * t362 + (-t325 * t362 + t326 * t364) * qJD(2)) * t359) * t356) + t359 * ((t338 * t312 + t339 * t313) * t359 + ((t338 * t248 + t339 * t250) * t355 - (t338 * t247 + t339 * t249) * t358 + (-t311 * t359 + (t245 * t358 - t246 * t355) * t356) * t364 + ((t279 * t359 + (-t229 * t358 + t230 * t355) * t356) * t362 + (t356 * t365 + t359 * t376) * t364) * qJD(2)) * t356) + 0.2e1 * m(3) * ((-t289 * t359 - t329 * t436) * (-t300 * t359 - t337 * t436) + (t290 * t359 - t329 * t437) * (t301 * t359 - t337 * t437) + (t289 * t355 + t290 * t358) * t356 * t183); (m(4) + m(5) + m(6) + m(7)) * t393; m(7) * (t331 * t68 + t333 * t69 - t368 * t43 + t342 * t44 + (-t35 * t364 + t410 * t59) * t356) + m(6) * (t331 * t88 + t333 * t89 - t368 * t66 + t342 * t67 + (-t364 * t50 + t410 * t76) * t356) + m(5) * (t331 * t108 + t333 * t109 - t368 * t90 + t342 * t91 + (-t364 * t75 + t410 * t92) * t356) + m(4) * (-t368 * t139 + t342 * t140 + t331 * t143 + t333 * t144 + (-t127 * t364 + t128 * t410) * t356); 0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t407) * (-t356 ^ 2 * t362 * t409 - t331 * t368 + t342 * t333); m(5) * t83 + m(6) * t47 + m(7) * t26; (t8 / 0.2e1 + t17 / 0.2e1 + t19 / 0.2e1) * t342 - (t7 / 0.2e1 + t16 / 0.2e1 + t18 / 0.2e1) * t368 + m(5) * (t108 * t96 + t109 * t95 + t134 * t75 + t150 * t91 + t151 * t90 + t83 * t92) + m(6) * (t123 * t67 + t124 * t66 + t47 * t76 + t50 * t97 + t60 * t89 + t61 * t88) + m(7) * (t26 * t59 + t35 * t70 + t37 * t69 + t38 * t68 + t43 * t87 + t44 * t86) + (t6 / 0.2e1 + t14 / 0.2e1 + t15 / 0.2e1 + (t131 / 0.2e1 + t133 / 0.2e1) * t333 + (t132 / 0.2e1 + t130 / 0.2e1) * t331) * t359 + t33 * t445 + t34 * t444 + ((-t20 / 0.2e1 - t21 / 0.2e1 - t9 / 0.2e1) * t364 + (t40 / 0.2e1 + (t141 / 0.2e1 + t142 / 0.2e1) * t359) * t410 + (t457 * t393 + t455) * t443 + (t458 * t393 + t456) * t442 + (-t461 * t444 - t445 * t463) * t358 + (t444 * t460 + t445 * t462) * t355) * t356; m(5) * (t150 * t333 + t151 * t331 - t96 * t368 + t95 * t342 + (t134 * t410 - t364 * t83) * t356) + m(6) * (t123 * t333 + t124 * t331 - t61 * t368 + t60 * t342 + (-t364 * t47 + t410 * t97) * t356) + m(7) * (t87 * t331 + t86 * t333 - t38 * t368 + t37 * t342 + (-t26 * t364 + t410 * t70) * t356); (-t14 - t15 - t6) * t434 + t455 * t342 - t456 * t368 + (t460 * t342 - t461 * t368 - t434 * t464 + t32) * t333 + (t462 * t342 - t463 * t368 - t434 * t465 + t31) * t331 + (t39 + (-t141 - t142) * t434 + t457 * t342 - t458 * t368) * t393 + (t26 * t70 + t37 * t86 + t38 * t87) * t452 + (t123 * t60 + t124 * t61 + t47 * t97) * t453 + (t134 * t83 + t150 * t95 + t151 * t96) * t454; t305 * t459; m(7) * (t258 * t68 + t260 * t69 + t305 * t59 + t307 * t43 + t309 * t44 + t327 * t35) + m(6) * (t258 * t88 + t260 * t89 + t305 * t76 + t307 * t66 + t309 * t67 + t327 * t50); (-t258 * t368 + t260 * t342 + t307 * t331 + t309 * t333 + (-t305 * t364 + t327 * t410) * t356) * t459; m(7) * (t258 * t87 + t26 * t327 + t260 * t86 + t305 * t70 + t307 * t38 + t309 * t37) + m(6) * (t123 * t260 + t124 * t258 + t305 * t97 + t307 * t61 + t309 * t60 + t327 * t47); 0.4e1 * t407 * (t258 * t307 + t260 * t309 + t305 * t327); m(7) * t45; t34 * t450 + t8 * t447 + t359 * t5 / 0.2e1 + t40 * t449 + t9 * t446 + t33 * t451 + t7 * t448 + m(7) * (t125 * t44 + t126 * t43 + t35 * t98 + t45 * t59 + t48 * t69 + t49 * t68) + (t1 * t442 + t2 * t443) * t356; m(7) * (t125 * t333 + t126 * t331 - t49 * t368 + t48 * t342 + (-t364 * t45 + t410 * t98) * t356); m(7) * (t125 * t37 + t126 * t38 + t26 * t98 + t45 * t70 + t48 * t86 + t49 * t87) + t30 * t444 + t342 * t2 / 0.2e1 + t31 * t451 + t3 * t448 + t39 * t449 + t6 * t446 + t32 * t450 + t4 * t447 + t29 * t445 - t368 * t1 / 0.2e1 + (t36 * t410 / 0.2e1 - t364 * t5 / 0.2e1) * t356; m(7) * (t125 * t260 + t126 * t258 + t305 * t98 + t307 * t49 + t309 * t48 + t327 * t45); t261 * t30 + t310 * t2 + t259 * t29 + t308 * t1 + t306 * t36 + t328 * t5 + (t125 * t48 + t126 * t49 + t45 * t98) * t452;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
