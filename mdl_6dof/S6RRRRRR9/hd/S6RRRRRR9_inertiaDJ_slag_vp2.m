% Calculate time derivative of joint inertia matrix for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:53
% EndTime: 2019-03-10 05:17:26
% DurationCPUTime: 15.08s
% Computational Cost: add. (36815->1092), mult. (108177->1572), div. (0->0), fcn. (111402->14), ass. (0->424)
t384 = sin(pkin(6));
t533 = 0.2e1 * t384;
t396 = cos(qJ(2));
t386 = cos(pkin(6));
t474 = pkin(1) * t386;
t375 = t396 * t474;
t366 = qJD(2) * t375;
t391 = sin(qJ(2));
t385 = cos(pkin(7));
t416 = t384 * (-pkin(10) * t385 - pkin(9));
t404 = t391 * t416;
t258 = pkin(2) * t386 + t375 + t404;
t461 = t258 * t385;
t532 = qJD(2) * t404 + qJD(3) * t461 + t366;
t388 = sin(qJ(5));
t393 = cos(qJ(5));
t389 = sin(qJ(4));
t440 = qJD(5) * t389;
t394 = cos(qJ(4));
t442 = qJD(4) * t394;
t398 = -t388 * t440 + t393 * t442;
t351 = -pkin(4) * t394 - pkin(12) * t389 - pkin(3);
t450 = t393 * t394;
t376 = pkin(11) * t450;
t295 = t388 * t351 + t376;
t531 = qJD(5) * t295;
t387 = sin(qJ(6));
t392 = cos(qJ(6));
t406 = t387 * t388 - t392 * t393;
t489 = -t406 / 0.2e1;
t334 = t387 * t393 + t388 * t392;
t488 = t334 / 0.2e1;
t478 = t388 / 0.2e1;
t476 = t393 / 0.2e1;
t395 = cos(qJ(3));
t451 = t391 * t395;
t390 = sin(qJ(3));
t452 = t390 * t396;
t400 = t385 * t452 + t451;
t383 = sin(pkin(7));
t445 = qJD(3) * t390;
t429 = t383 * t445;
t202 = t386 * t429 + (t400 * qJD(3) + (t385 * t451 + t452) * qJD(2)) * t384;
t444 = qJD(3) * t395;
t428 = t383 * t444;
t449 = t395 * t396;
t453 = t390 * t391;
t527 = t385 * t449 - t453;
t203 = t386 * t428 + (t527 * qJD(3) + (-t385 * t453 + t449) * qJD(2)) * t384;
t374 = t391 * t474;
t260 = (t396 * t416 - t374) * qJD(2);
t447 = qJD(2) * t384;
t473 = pkin(10) * t383;
t293 = (pkin(2) * t391 - t396 * t473) * t447;
t204 = -t260 * t383 + t385 * t293;
t103 = pkin(3) * t202 - pkin(11) * t203 + t204;
t289 = (-pkin(2) * t396 - t391 * t473 - pkin(1)) * t384;
t195 = -t258 * t383 + t385 * t289;
t459 = t383 * t395;
t248 = -t384 * t527 - t386 * t459;
t460 = t383 * t390;
t249 = t384 * t400 + t386 * t460;
t133 = pkin(3) * t248 - pkin(11) * t249 + t195;
t458 = t384 * t396;
t324 = pkin(9) * t458 + t374;
t242 = (t383 * t386 + t385 * t458) * pkin(10) + t324;
t457 = t385 * t390;
t147 = t395 * t242 + t258 * t457 + t289 * t460;
t316 = -t383 * t458 + t386 * t385;
t139 = pkin(11) * t316 + t147;
t443 = qJD(4) * t389;
t431 = t391 * t447;
t420 = t383 * t431;
t92 = -t242 * t445 + t260 * t457 + t289 * t428 + t293 * t460 + t532 * t395;
t87 = pkin(11) * t420 + t92;
t32 = t103 * t394 - t133 * t443 - t139 * t442 - t389 * t87;
t29 = -pkin(4) * t202 - t32;
t208 = t249 * t389 - t316 * t394;
t125 = -qJD(4) * t208 + t203 * t394 + t389 * t420;
t209 = t249 * t394 + t316 * t389;
t151 = t209 * t393 + t248 * t388;
t62 = -qJD(5) * t151 - t125 * t388 + t202 * t393;
t150 = -t209 * t388 + t248 * t393;
t63 = qJD(5) * t150 + t125 * t393 + t202 * t388;
t34 = -mrSges(6,1) * t62 + mrSges(6,2) * t63;
t530 = -m(6) * t29 - t34;
t315 = t406 * t389;
t317 = -t394 * t385 + t389 * t460;
t265 = -qJD(4) * t317 + t394 * t428;
t318 = t385 * t389 + t394 * t460;
t267 = -t318 * t388 - t393 * t459;
t169 = qJD(5) * t267 + t265 * t393 + t388 * t429;
t401 = -t318 * t393 + t388 * t459;
t170 = qJD(5) * t401 - t265 * t388 + t393 * t429;
t106 = -mrSges(6,1) * t170 + mrSges(6,2) * t169;
t323 = pkin(2) * t457 + pkin(10) * t459;
t300 = pkin(11) * t385 + t323;
t301 = (-pkin(3) * t395 - pkin(11) * t390 - pkin(2)) * t383;
t446 = qJD(3) * t383;
t309 = (pkin(3) * t390 - pkin(11) * t395) * t446;
t369 = pkin(10) * t460;
t456 = t385 * t395;
t321 = pkin(2) * t456 - t369;
t310 = t321 * qJD(3);
t153 = -t300 * t442 - t301 * t443 + t309 * t394 - t389 * t310;
t149 = -pkin(4) * t429 - t153;
t529 = -m(6) * t149 - t106;
t299 = t369 + (-pkin(2) * t395 - pkin(3)) * t385;
t210 = pkin(4) * t317 - pkin(12) * t318 + t299;
t220 = t394 * t300 + t389 * t301;
t212 = -pkin(12) * t459 + t220;
t141 = t388 * t210 + t393 * t212;
t528 = -m(5) * pkin(3) - mrSges(5,1) * t394 + mrSges(5,2) * t389;
t526 = qJD(5) + qJD(6);
t146 = -t390 * t242 + t395 * (t289 * t383 + t461);
t525 = 2 * m(4);
t524 = 0.2e1 * m(5);
t523 = 0.2e1 * m(6);
t522 = 2 * m(7);
t521 = 0.2e1 * pkin(11);
t520 = -2 * mrSges(3,3);
t519 = -2 * mrSges(4,3);
t124 = qJD(4) * t209 + t203 * t389 - t394 * t420;
t55 = Ifges(5,1) * t125 - Ifges(5,4) * t124 + Ifges(5,5) * t202;
t517 = t55 / 0.2e1;
t96 = t150 * t392 - t151 * t387;
t516 = t96 / 0.2e1;
t97 = t150 * t387 + t151 * t392;
t515 = t97 / 0.2e1;
t514 = -pkin(13) - pkin(12);
t114 = Ifges(5,1) * t209 - Ifges(5,4) * t208 + Ifges(5,5) * t248;
t512 = t114 / 0.2e1;
t266 = qJD(4) * t318 + t389 * t428;
t175 = Ifges(5,1) * t265 - Ifges(5,4) * t266 + Ifges(5,5) * t429;
t511 = t175 / 0.2e1;
t181 = t267 * t392 + t387 * t401;
t510 = t181 / 0.2e1;
t182 = t267 * t387 - t392 * t401;
t509 = t182 / 0.2e1;
t256 = t526 * t406;
t257 = t526 * t334;
t185 = -Ifges(7,4) * t256 - Ifges(7,2) * t257;
t508 = t185 / 0.2e1;
t186 = -Ifges(7,1) * t256 - Ifges(7,4) * t257;
t507 = t186 / 0.2e1;
t192 = -t257 * t389 - t406 * t442;
t506 = t192 / 0.2e1;
t193 = t315 * t526 - t334 * t442;
t505 = t193 / 0.2e1;
t314 = t334 * t389;
t222 = -Ifges(7,4) * t315 - Ifges(7,2) * t314 - Ifges(7,6) * t394;
t504 = t222 / 0.2e1;
t223 = -Ifges(7,1) * t315 - Ifges(7,4) * t314 - Ifges(7,5) * t394;
t503 = t223 / 0.2e1;
t226 = Ifges(5,1) * t318 - Ifges(5,4) * t317 - Ifges(5,5) * t459;
t502 = t226 / 0.2e1;
t466 = Ifges(6,4) * t388;
t356 = Ifges(6,2) * t393 + t466;
t465 = Ifges(6,4) * t393;
t409 = -Ifges(6,2) * t388 + t465;
t230 = -t356 * t440 + (Ifges(6,6) * t389 + t394 * t409) * qJD(4);
t501 = t230 / 0.2e1;
t358 = Ifges(6,1) * t388 + t465;
t410 = Ifges(6,1) * t393 - t466;
t231 = -t358 * t440 + (Ifges(6,5) * t389 + t394 * t410) * qJD(4);
t500 = t231 / 0.2e1;
t499 = -t256 / 0.2e1;
t498 = -t257 / 0.2e1;
t263 = Ifges(7,4) * t334 - Ifges(7,2) * t406;
t497 = t263 / 0.2e1;
t264 = Ifges(7,1) * t334 - Ifges(7,4) * t406;
t496 = t264 / 0.2e1;
t495 = t267 / 0.2e1;
t494 = -t401 / 0.2e1;
t303 = -Ifges(6,6) * t394 + t389 * t409;
t493 = t303 / 0.2e1;
t304 = -Ifges(6,5) * t394 + t389 * t410;
t492 = t304 / 0.2e1;
t491 = -t314 / 0.2e1;
t490 = -t315 / 0.2e1;
t339 = t409 * qJD(5);
t487 = t339 / 0.2e1;
t341 = t410 * qJD(5);
t486 = t341 / 0.2e1;
t468 = Ifges(5,4) * t389;
t342 = (Ifges(5,1) * t394 - t468) * qJD(4);
t485 = t342 / 0.2e1;
t484 = Ifges(5,5) * t389 / 0.2e1 + Ifges(5,6) * t394 / 0.2e1;
t483 = t356 / 0.2e1;
t482 = t358 / 0.2e1;
t467 = Ifges(5,4) * t394;
t359 = Ifges(5,1) * t389 + t467;
t481 = t359 / 0.2e1;
t480 = t385 / 0.2e1;
t479 = -t388 / 0.2e1;
t477 = -t393 / 0.2e1;
t472 = pkin(11) * t388;
t71 = t389 * t133 + t394 * t139;
t66 = pkin(12) * t248 + t71;
t138 = -pkin(3) * t316 - t146;
t80 = pkin(4) * t208 - pkin(12) * t209 + t138;
t36 = t388 * t80 + t393 * t66;
t471 = m(6) * qJD(4);
t470 = Ifges(4,4) * t390;
t469 = Ifges(4,4) * t395;
t464 = Ifges(6,6) * t388;
t312 = -pkin(9) * t431 + t366;
t463 = t312 * mrSges(3,2);
t313 = t324 * qJD(2);
t462 = t313 * mrSges(3,1);
t455 = t388 * t389;
t454 = t389 * t393;
t184 = -Ifges(7,5) * t256 - Ifges(7,6) * t257;
t347 = (pkin(4) * t389 - pkin(12) * t394) * qJD(4);
t448 = t393 * t347 + t443 * t472;
t441 = qJD(5) * t388;
t439 = qJD(5) * t393;
t438 = qJD(6) * t387;
t437 = qJD(6) * t392;
t22 = qJD(6) * t96 + t387 * t62 + t392 * t63;
t23 = -qJD(6) * t97 - t387 * t63 + t392 * t62;
t8 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t124;
t24 = Ifges(6,5) * t63 + Ifges(6,6) * t62 + Ifges(6,3) * t124;
t83 = qJD(6) * t181 + t169 * t392 + t170 * t387;
t84 = -qJD(6) * t182 - t169 * t387 + t170 * t392;
t40 = Ifges(7,5) * t83 + Ifges(7,6) * t84 + Ifges(7,3) * t266;
t435 = pkin(5) * t441;
t76 = Ifges(6,1) * t151 + Ifges(6,4) * t150 + Ifges(6,5) * t208;
t433 = t76 * t476;
t53 = Ifges(5,5) * t125 - Ifges(5,6) * t124 + Ifges(5,3) * t202;
t99 = Ifges(6,5) * t169 + Ifges(6,6) * t170 + Ifges(6,3) * t266;
t126 = Ifges(7,5) * t192 + Ifges(7,6) * t193 + Ifges(7,3) * t443;
t121 = Ifges(4,5) * t203 - Ifges(4,6) * t202 + Ifges(4,3) * t420;
t173 = Ifges(5,5) * t265 - Ifges(5,6) * t266 + Ifges(5,3) * t429;
t432 = qJD(5) * t514;
t165 = -Ifges(6,1) * t401 + Ifges(6,4) * t267 + Ifges(6,5) * t317;
t425 = t165 * t476;
t381 = Ifges(6,5) * t439;
t424 = -Ifges(6,6) * t441 / 0.2e1 + t381 / 0.2e1 + t184 / 0.2e1;
t423 = Ifges(6,5) * t478 + Ifges(7,5) * t488 + Ifges(6,6) * t476 + Ifges(7,6) * t489;
t35 = -t388 * t66 + t393 * t80;
t70 = t133 * t394 - t389 * t139;
t140 = t393 * t210 - t212 * t388;
t219 = -t389 * t300 + t301 * t394;
t217 = t388 * t347 + t351 * t439 + (-t393 * t443 - t394 * t441) * pkin(11);
t330 = t393 * t351;
t294 = -t394 * t472 + t330;
t422 = -qJD(5) * t294 + t217;
t421 = -t242 * t444 - t289 * t429 - t532 * t390;
t54 = Ifges(5,4) * t125 - Ifges(5,2) * t124 + Ifges(5,6) * t202;
t419 = t24 / 0.2e1 + t8 / 0.2e1 - t54 / 0.2e1;
t113 = Ifges(5,4) * t209 - Ifges(5,2) * t208 + Ifges(5,6) * t248;
t44 = Ifges(7,5) * t97 + Ifges(7,6) * t96 + Ifges(7,3) * t208;
t74 = Ifges(6,5) * t151 + Ifges(6,6) * t150 + Ifges(6,3) * t208;
t418 = -t113 / 0.2e1 + t74 / 0.2e1 + t44 / 0.2e1;
t174 = Ifges(5,4) * t265 - Ifges(5,2) * t266 + Ifges(5,6) * t429;
t417 = t99 / 0.2e1 + t40 / 0.2e1 - t174 / 0.2e1;
t31 = t389 * t103 + t133 * t442 - t139 * t443 + t394 * t87;
t28 = pkin(12) * t202 + t31;
t88 = -t260 * t456 + (-pkin(3) * t431 - t293 * t395) * t383 - t421;
t39 = pkin(4) * t124 - pkin(12) * t125 + t88;
t6 = t393 * t28 + t388 * t39 + t80 * t439 - t441 * t66;
t7 = -qJD(5) * t36 - t28 * t388 + t393 * t39;
t415 = -t7 * t388 + t6 * t393;
t211 = pkin(4) * t459 - t219;
t107 = Ifges(7,5) * t182 + Ifges(7,6) * t181 + Ifges(7,3) * t317;
t163 = -Ifges(6,5) * t401 + Ifges(6,6) * t267 + Ifges(6,3) * t317;
t225 = Ifges(5,4) * t318 - Ifges(5,2) * t317 - Ifges(5,6) * t459;
t414 = -t225 / 0.2e1 + t163 / 0.2e1 + t107 / 0.2e1;
t397 = t388 * t442 + t389 * t439;
t229 = t398 * Ifges(6,5) - Ifges(6,6) * t397 + Ifges(6,3) * t443;
t340 = (-Ifges(5,2) * t389 + t467) * qJD(4);
t413 = -t340 / 0.2e1 + t229 / 0.2e1 + t126 / 0.2e1;
t221 = -Ifges(7,5) * t315 - Ifges(7,6) * t314 - Ifges(7,3) * t394;
t302 = -Ifges(6,3) * t394 + (Ifges(6,5) * t393 - t464) * t389;
t357 = Ifges(5,2) * t394 + t468;
t412 = -t357 / 0.2e1 + t302 / 0.2e1 + t221 / 0.2e1;
t352 = -mrSges(6,1) * t393 + mrSges(6,2) * t388;
t411 = mrSges(6,1) * t388 + mrSges(6,2) * t393;
t30 = pkin(5) * t208 - pkin(13) * t151 + t35;
t33 = pkin(13) * t150 + t36;
t12 = t30 * t392 - t33 * t387;
t13 = t30 * t387 + t33 * t392;
t152 = -t300 * t443 + t301 * t442 + t389 * t309 + t394 * t310;
t148 = pkin(12) * t429 + t152;
t311 = t323 * qJD(3);
t166 = pkin(4) * t266 - pkin(12) * t265 + t311;
t58 = t393 * t148 + t388 * t166 + t210 * t439 - t212 * t441;
t59 = -qJD(5) * t141 - t148 * t388 + t393 * t166;
t408 = -t59 * t388 + t58 * t393;
t105 = pkin(5) * t317 + pkin(13) * t401 + t140;
t115 = pkin(13) * t267 + t141;
t56 = t105 * t392 - t115 * t387;
t57 = t105 * t387 + t115 * t392;
t245 = -pkin(13) * t454 + t330 + (-pkin(5) - t472) * t394;
t270 = -pkin(13) * t455 + t295;
t179 = t245 * t392 - t270 * t387;
t180 = t245 * t387 + t270 * t392;
t361 = t514 * t388;
t362 = t514 * t393;
t275 = t361 * t392 + t362 * t387;
t276 = t361 * t387 - t362 * t392;
t345 = t388 * t432;
t346 = t393 * t432;
t200 = qJD(6) * t275 + t345 * t392 + t346 * t387;
t201 = -qJD(6) * t276 - t345 * t387 + t346 * t392;
t405 = t201 * mrSges(7,1) - t200 * mrSges(7,2) + t184;
t4 = pkin(5) * t124 - pkin(13) * t63 + t7;
t5 = pkin(13) * t62 + t6;
t2 = qJD(6) * t12 + t387 * t4 + t392 * t5;
t3 = -qJD(6) * t13 - t387 * t5 + t392 * t4;
t403 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t8;
t65 = -pkin(4) * t248 - t70;
t49 = pkin(5) * t266 - pkin(13) * t169 + t59;
t51 = pkin(13) * t170 + t58;
t15 = qJD(6) * t56 + t387 * t49 + t392 * t51;
t16 = -qJD(6) * t57 - t387 * t51 + t392 * t49;
t402 = t16 * mrSges(7,1) - t15 * mrSges(7,2) + t40;
t176 = (pkin(5) * t389 - pkin(13) * t450) * qJD(4) + (-t376 + (pkin(13) * t389 - t351) * t388) * qJD(5) + t448;
t187 = -pkin(13) * t397 + t217;
t90 = qJD(6) * t179 + t176 * t387 + t187 * t392;
t91 = -qJD(6) * t180 + t176 * t392 - t187 * t387;
t399 = t91 * mrSges(7,1) - t90 * mrSges(7,2) + t126;
t382 = Ifges(5,5) * t442;
t378 = -pkin(5) * t393 - pkin(4);
t365 = Ifges(3,5) * t396 * t447;
t364 = Ifges(4,5) * t428;
t349 = (pkin(5) * t388 + pkin(11)) * t389;
t344 = -mrSges(6,1) * t394 - mrSges(6,3) * t454;
t343 = mrSges(6,2) * t394 - mrSges(6,3) * t455;
t338 = -Ifges(5,6) * t443 + t382;
t336 = (mrSges(5,1) * t389 + mrSges(5,2) * t394) * qJD(4);
t335 = t411 * qJD(5);
t332 = -mrSges(4,2) * t385 + mrSges(4,3) * t459;
t331 = mrSges(4,1) * t385 - mrSges(4,3) * t460;
t327 = (-mrSges(7,1) * t387 - mrSges(7,2) * t392) * qJD(6) * pkin(5);
t325 = t411 * t389;
t322 = -pkin(9) * t384 * t391 + t375;
t308 = (Ifges(4,1) * t395 - t470) * t446;
t307 = (-Ifges(4,2) * t390 + t469) * t446;
t306 = -Ifges(4,6) * t429 + t364;
t305 = (mrSges(4,1) * t390 + mrSges(4,2) * t395) * t446;
t297 = Ifges(4,5) * t385 + (Ifges(4,1) * t390 + t469) * t383;
t296 = Ifges(4,6) * t385 + (Ifges(4,2) * t395 + t470) * t383;
t288 = pkin(5) * t397 + pkin(11) * t442;
t284 = -mrSges(6,2) * t443 - mrSges(6,3) * t397;
t283 = mrSges(6,1) * t443 - mrSges(6,3) * t398;
t278 = -mrSges(7,1) * t394 + mrSges(7,3) * t315;
t277 = mrSges(7,2) * t394 - mrSges(7,3) * t314;
t274 = -mrSges(5,1) * t459 - mrSges(5,3) * t318;
t273 = mrSges(5,2) * t459 - mrSges(5,3) * t317;
t261 = mrSges(7,1) * t406 + mrSges(7,2) * t334;
t241 = mrSges(6,1) * t397 + mrSges(6,2) * t398;
t236 = mrSges(5,1) * t317 + mrSges(5,2) * t318;
t232 = mrSges(7,1) * t314 - mrSges(7,2) * t315;
t228 = -mrSges(5,2) * t429 - mrSges(5,3) * t266;
t227 = mrSges(5,1) * t429 - mrSges(5,3) * t265;
t224 = Ifges(5,5) * t318 - Ifges(5,6) * t317 - Ifges(5,3) * t459;
t218 = t448 - t531;
t216 = mrSges(6,1) * t317 + mrSges(6,3) * t401;
t215 = -mrSges(6,2) * t317 + mrSges(6,3) * t267;
t214 = mrSges(4,1) * t316 - mrSges(4,3) * t249;
t213 = -mrSges(4,2) * t316 - mrSges(4,3) * t248;
t191 = -mrSges(6,1) * t267 - mrSges(6,2) * t401;
t188 = mrSges(5,1) * t266 + mrSges(5,2) * t265;
t183 = mrSges(7,1) * t257 - mrSges(7,2) * t256;
t178 = -mrSges(7,2) * t443 + mrSges(7,3) * t193;
t177 = mrSges(7,1) * t443 - mrSges(7,3) * t192;
t172 = mrSges(4,1) * t420 - mrSges(4,3) * t203;
t171 = -mrSges(4,2) * t420 - mrSges(4,3) * t202;
t164 = -Ifges(6,4) * t401 + Ifges(6,2) * t267 + Ifges(6,6) * t317;
t160 = Ifges(4,1) * t249 - Ifges(4,4) * t248 + Ifges(4,5) * t316;
t159 = Ifges(4,4) * t249 - Ifges(4,2) * t248 + Ifges(4,6) * t316;
t158 = -pkin(5) * t267 + t211;
t157 = mrSges(7,1) * t317 - mrSges(7,3) * t182;
t156 = -mrSges(7,2) * t317 + mrSges(7,3) * t181;
t155 = mrSges(5,1) * t248 - mrSges(5,3) * t209;
t154 = -mrSges(5,2) * t248 - mrSges(5,3) * t208;
t144 = -mrSges(6,2) * t266 + mrSges(6,3) * t170;
t143 = mrSges(6,1) * t266 - mrSges(6,3) * t169;
t142 = mrSges(5,1) * t208 + mrSges(5,2) * t209;
t137 = mrSges(4,1) * t202 + mrSges(4,2) * t203;
t132 = -mrSges(7,1) * t193 + mrSges(7,2) * t192;
t128 = Ifges(7,1) * t192 + Ifges(7,4) * t193 + Ifges(7,5) * t443;
t127 = Ifges(7,4) * t192 + Ifges(7,2) * t193 + Ifges(7,6) * t443;
t123 = Ifges(4,1) * t203 - Ifges(4,4) * t202 + Ifges(4,5) * t420;
t122 = Ifges(4,4) * t203 - Ifges(4,2) * t202 + Ifges(4,6) * t420;
t120 = -mrSges(7,1) * t181 + mrSges(7,2) * t182;
t112 = Ifges(5,5) * t209 - Ifges(5,6) * t208 + Ifges(5,3) * t248;
t111 = mrSges(6,1) * t208 - mrSges(6,3) * t151;
t110 = -mrSges(6,2) * t208 + mrSges(6,3) * t150;
t109 = Ifges(7,1) * t182 + Ifges(7,4) * t181 + Ifges(7,5) * t317;
t108 = Ifges(7,4) * t182 + Ifges(7,2) * t181 + Ifges(7,6) * t317;
t104 = -pkin(5) * t170 + t149;
t101 = Ifges(6,1) * t169 + Ifges(6,4) * t170 + Ifges(6,5) * t266;
t100 = Ifges(6,4) * t169 + Ifges(6,2) * t170 + Ifges(6,6) * t266;
t98 = -mrSges(6,1) * t150 + mrSges(6,2) * t151;
t95 = mrSges(5,1) * t202 - mrSges(5,3) * t125;
t94 = -mrSges(5,2) * t202 - mrSges(5,3) * t124;
t93 = (t260 * t385 + t293 * t383) * t395 + t421;
t75 = Ifges(6,4) * t151 + Ifges(6,2) * t150 + Ifges(6,6) * t208;
t73 = mrSges(7,1) * t208 - mrSges(7,3) * t97;
t72 = -mrSges(7,2) * t208 + mrSges(7,3) * t96;
t69 = -mrSges(7,2) * t266 + mrSges(7,3) * t84;
t68 = mrSges(7,1) * t266 - mrSges(7,3) * t83;
t67 = mrSges(5,1) * t124 + mrSges(5,2) * t125;
t52 = -pkin(5) * t150 + t65;
t50 = -mrSges(7,1) * t96 + mrSges(7,2) * t97;
t48 = mrSges(6,1) * t124 - mrSges(6,3) * t63;
t47 = -mrSges(6,2) * t124 + mrSges(6,3) * t62;
t46 = Ifges(7,1) * t97 + Ifges(7,4) * t96 + Ifges(7,5) * t208;
t45 = Ifges(7,4) * t97 + Ifges(7,2) * t96 + Ifges(7,6) * t208;
t43 = -mrSges(7,1) * t84 + mrSges(7,2) * t83;
t42 = Ifges(7,1) * t83 + Ifges(7,4) * t84 + Ifges(7,5) * t266;
t41 = Ifges(7,4) * t83 + Ifges(7,2) * t84 + Ifges(7,6) * t266;
t26 = Ifges(6,1) * t63 + Ifges(6,4) * t62 + Ifges(6,5) * t124;
t25 = Ifges(6,4) * t63 + Ifges(6,2) * t62 + Ifges(6,6) * t124;
t19 = -mrSges(7,2) * t124 + mrSges(7,3) * t23;
t18 = mrSges(7,1) * t124 - mrSges(7,3) * t22;
t17 = -pkin(5) * t62 + t29;
t11 = -mrSges(7,1) * t23 + mrSges(7,2) * t22;
t10 = Ifges(7,1) * t22 + Ifges(7,4) * t23 + Ifges(7,5) * t124;
t9 = Ifges(7,4) * t22 + Ifges(7,2) * t23 + Ifges(7,6) * t124;
t1 = [t316 * t121 + t248 * t53 - t248 * t122 + 0.2e1 * t204 * (mrSges(4,1) * t248 + mrSges(4,2) * t249) + t249 * t123 + 0.2e1 * t92 * t213 + 0.2e1 * t93 * t214 + t209 * t55 + t203 * t160 + 0.2e1 * t195 * t137 + 0.2e1 * t146 * t172 + 0.2e1 * t147 * t171 + t151 * t26 + 0.2e1 * t31 * t154 + 0.2e1 * t32 * t155 + t150 * t25 + 0.2e1 * t88 * t142 + 0.2e1 * t138 * t67 + t125 * t114 + 0.2e1 * t6 * t110 + 0.2e1 * t7 * t111 + 0.2e1 * t70 * t95 + t96 * t9 + t97 * t10 + 0.2e1 * t29 * t98 + 0.2e1 * t71 * t94 + 0.2e1 * t3 * t73 + t62 * t75 + t63 * t76 + 0.2e1 * t2 * t72 + 0.2e1 * t65 * t34 + 0.2e1 * t52 * t11 + 0.2e1 * t17 * t50 + t23 * t45 + t22 * t46 + 0.2e1 * t36 * t47 + 0.2e1 * t35 * t48 + 0.2e1 * t13 * t19 + 0.2e1 * t12 * t18 + (t12 * t3 + t13 * t2 + t17 * t52) * t522 + (t29 * t65 + t35 * t7 + t36 * t6) * t523 + (t138 * t88 + t31 * t71 + t32 * t70) * t524 + (t146 * t93 + t147 * t92 + t195 * t204) * t525 + (t365 - 0.2e1 * t462 - 0.2e1 * t463) * t386 + 0.2e1 * m(3) * (t312 * t324 - t313 * t322) + (t24 + t8 - t54) * t208 + (-t159 + t112) * t202 + (t74 + t44 - t113) * t124 + (0.2e1 * (t312 * t396 + t313 * t391) * mrSges(3,3) + ((t322 * t520 + Ifges(3,5) * t386 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t396) * t533) * t396 + (-0.2e1 * Ifges(3,6) * t386 + t383 * (Ifges(4,5) * t249 - Ifges(4,6) * t248 + Ifges(4,3) * t316) + t324 * t520 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t391 + (Ifges(3,1) - Ifges(3,2)) * t396) * t533) * t391) * qJD(2)) * t384; t323 * t171 + t93 * t331 + t92 * t332 + t321 * t172 + t316 * t306 / 0.2e1 + t195 * t305 + t249 * t308 / 0.2e1 + t310 * t213 + t203 * t297 / 0.2e1 + t299 * t67 + t32 * t274 + t31 * t273 + t88 * t236 + t220 * t94 + t70 * t227 + t71 * t228 + t211 * t34 + t6 * t215 + t7 * t216 + t219 * t95 + t138 * t188 + t29 * t191 + t169 * t76 / 0.2e1 + t170 * t75 / 0.2e1 + t158 * t11 + t62 * t164 / 0.2e1 + t63 * t165 / 0.2e1 + t151 * t101 / 0.2e1 + t152 * t154 + t153 * t155 + t2 * t156 + t3 * t157 + t149 * t98 + t150 * t100 / 0.2e1 + t35 * t143 + t36 * t144 + t140 * t48 + t141 * t47 + t17 * t120 + t22 * t109 / 0.2e1 + t58 * t110 + t59 * t111 + t104 * t50 + t65 * t106 + t23 * t108 / 0.2e1 + m(4) * (-t146 * t311 + t147 * t310 + t321 * t93 + t323 * t92) + t83 * t46 / 0.2e1 + t84 * t45 / 0.2e1 + t15 * t72 + t16 * t73 + t12 * t68 + t13 * t69 + t56 * t18 + t57 * t19 + t52 * t43 + ((t204 * mrSges(4,2) + t123 / 0.2e1) * t390 + (-t204 * mrSges(4,1) + t122 / 0.2e1 - t53 / 0.2e1) * t395 + (Ifges(4,3) * t480 + (Ifges(4,5) * t390 + Ifges(4,6) * t395) * t383 / 0.2e1) * t431 + (-m(4) * t204 - t137) * pkin(2) + ((t160 / 0.2e1 - t146 * mrSges(4,3)) * t395 + (-t159 / 0.2e1 + t112 / 0.2e1 - t147 * mrSges(4,3)) * t390) * qJD(3)) * t383 + t10 * t509 + t9 * t510 + t209 * t511 + t265 * t512 + t42 * t515 + t41 * t516 + t318 * t517 + t26 * t494 + t25 * t495 + t125 * t502 + t121 * t480 - Ifges(3,6) * t431 + t419 * t317 + t417 * t208 + t418 * t266 + t414 * t124 + m(5) * (t138 * t311 + t152 * t71 + t153 * t70 + t219 * t32 + t220 * t31 + t299 * t88) + m(6) * (t140 * t7 + t141 * t6 + t149 * t65 + t211 * t29 + t35 * t59 + t36 * t58) + m(7) * (t104 * t52 + t12 * t16 + t13 * t15 + t158 * t17 + t2 * t57 + t3 * t56) + t365 + (-t307 / 0.2e1 + t173 / 0.2e1) * t248 + (-t296 / 0.2e1 + t224 / 0.2e1) * t202 + (t142 - t214) * t311 - t462 - t463; -t401 * t101 + t385 * t306 + 0.2e1 * t310 * t332 + t318 * t175 + 0.2e1 * t299 * t188 + 0.2e1 * t153 * t274 + t267 * t100 + 0.2e1 * t152 * t273 + t265 * t226 + 0.2e1 * t219 * t227 + 0.2e1 * t220 * t228 + 0.2e1 * t58 * t215 + 0.2e1 * t59 * t216 + 0.2e1 * t211 * t106 + 0.2e1 * t149 * t191 + t181 * t41 + t182 * t42 + t169 * t165 + t170 * t164 + 0.2e1 * t158 * t43 + 0.2e1 * t15 * t156 + 0.2e1 * t16 * t157 + 0.2e1 * t140 * t143 + 0.2e1 * t141 * t144 + 0.2e1 * t104 * t120 + t83 * t109 + t84 * t108 + 0.2e1 * t56 * t68 + 0.2e1 * t57 * t69 + (-0.2e1 * pkin(2) * t305 + t390 * t308 + (-t173 + t307) * t395 + ((t321 * t519 + t297) * t395 + (t323 * t519 + t224 - t296) * t390) * qJD(3)) * t383 + (t104 * t158 + t15 * t57 + t16 * t56) * t522 + (t140 * t59 + t141 * t58 + t149 * t211) * t523 + (t152 * t220 + t153 * t219 + t299 * t311) * t524 + (t310 * t323 - t311 * t321) * t525 + (t99 + t40 - t174) * t317 + (-t225 + t163 + t107) * t266 + 0.2e1 * (-t331 + t236) * t311; t138 * t336 + t248 * t338 / 0.2e1 + t6 * t343 + t7 * t344 + t349 * t11 + t29 * t325 + t288 * t50 + t294 * t48 + t295 * t47 + t2 * t277 + t3 * t278 + t35 * t283 + t36 * t284 + t65 * t241 + t17 * t232 + t217 * t110 + t218 * t111 + t179 * t18 + t180 * t19 + t12 * t177 + t13 * t178 + t52 * t132 + t91 * t73 - t92 * mrSges(4,2) + t93 * mrSges(4,1) + t90 * t72 - pkin(3) * t67 + (t31 * mrSges(5,3) + (-t70 * mrSges(5,3) + t75 * t479 + t433 + t512) * qJD(4) + (t94 + (-t155 + t98) * qJD(4) + t65 * t471 + m(5) * (-qJD(4) * t70 + t31)) * pkin(11) - t419) * t394 + t128 * t515 + t127 * t516 + t9 * t491 + t63 * t492 + t62 * t493 + t151 * t500 + t150 * t501 + t22 * t503 + t23 * t504 + t45 * t505 + t46 * t506 + t125 * t481 + t202 * t484 + t209 * t485 + t10 * t490 + t412 * t124 + t413 * t208 + m(7) * (t12 * t91 + t13 * t90 + t17 * t349 + t179 * t3 + t180 * t2 + t288 * t52) + (t517 + t26 * t476 + t25 * t479 - t32 * mrSges(5,3) + (t477 * t75 + t479 * t76) * qJD(5) + (-t71 * mrSges(5,3) + t418) * qJD(4) + (-qJD(4) * t154 - t95 + m(5) * (-qJD(4) * t71 - t32) - t530) * pkin(11)) * t389 + t528 * t88 + t121 + m(6) * (t217 * t36 + t218 * t35 + t294 * t7 + t295 * t6); t299 * t336 + t58 * t343 + t59 * t344 + t349 * t43 + t149 * t325 - t310 * mrSges(4,2) + t288 * t120 + t294 * t143 + t295 * t144 + t15 * t277 + t16 * t278 + t140 * t283 + t141 * t284 + t211 * t241 + t104 * t232 + t217 * t215 + t218 * t216 - pkin(3) * t188 + t179 * t68 + t180 * t69 + t56 * t177 + t57 * t178 + t90 * t156 + t91 * t157 + t158 * t132 + (t152 * mrSges(5,3) + (-t219 * mrSges(5,3) + t164 * t479 + t425 + t502) * qJD(4) + (t228 + (t191 - t274) * qJD(4) + t211 * t471 + m(5) * (-qJD(4) * t219 + t152)) * pkin(11) - t417) * t394 + (-t395 * t338 / 0.2e1 + (-Ifges(4,6) + t484) * t445) * t383 + t128 * t509 + t127 * t510 + t41 * t491 + t169 * t492 + t170 * t493 + t231 * t494 + t230 * t495 + t83 * t503 + t84 * t504 + t108 * t505 + t109 * t506 + t265 * t481 + t318 * t485 + t42 * t490 + t412 * t266 + t413 * t317 + t364 + m(6) * (t140 * t218 + t141 * t217 + t294 * t59 + t295 * t58) + (-mrSges(4,1) + t528) * t311 + (t511 + t101 * t476 + t100 * t479 - t153 * mrSges(5,3) + (t164 * t477 + t165 * t479) * qJD(5) + (-t220 * mrSges(5,3) + t414) * qJD(4) + (-qJD(4) * t273 - t227 + m(5) * (-qJD(4) * t220 - t153) - t529) * pkin(11)) * t389 + m(7) * (t104 * t349 + t15 * t180 + t158 * t288 + t16 * t179 + t56 * t91 + t57 * t90); -0.2e1 * pkin(3) * t336 - t314 * t127 - t315 * t128 + 0.2e1 * t349 * t132 + 0.2e1 * t179 * t177 + 0.2e1 * t180 * t178 + t192 * t223 + t193 * t222 + 0.2e1 * t217 * t343 + 0.2e1 * t218 * t344 + 0.2e1 * t288 * t232 + 0.2e1 * t90 * t277 + 0.2e1 * t91 * t278 + 0.2e1 * t294 * t283 + 0.2e1 * t295 * t284 + (t217 * t295 + t218 * t294) * t523 + (t179 * t91 + t180 * t90 + t288 * t349) * t522 + (-t126 - t229 + t340 + (-t303 * t388 + t304 * t393 + t325 * t521 + t359) * qJD(4)) * t394 + (t241 * t521 - t388 * t230 + t393 * t231 + t342 + (-t303 * t393 - t304 * t388) * qJD(5) + (pkin(11) ^ 2 * t394 * t523 + t221 + t302 - t357) * qJD(4)) * t389; (t12 * t256 - t13 * t257 - t2 * t406 - t3 * t334) * mrSges(7,3) + t53 + t378 * t11 + t29 * t352 + t65 * t335 + t275 * t18 + t276 * t19 + t17 * t261 + t200 * t72 + t201 * t73 + t52 * t183 - t31 * mrSges(5,2) + t32 * mrSges(5,1) + (-t111 * t439 - t110 * t441 + m(6) * (-t35 * t439 - t36 * t441 + t415) + t393 * t47 - t388 * t48) * pkin(12) + m(7) * (t12 * t201 + t13 * t200 + t17 * t378 + t2 * t276 + t275 * t3 + t435 * t52) + t97 * t507 + t96 * t508 + t22 * t496 + t23 * t497 + t45 * t498 + t46 * t499 + t63 * t482 + t62 * t483 + t151 * t486 + t150 * t487 + t10 * t488 + t9 * t489 + t25 * t476 + t26 * t478 + (t433 + (pkin(5) * t50 - t75 / 0.2e1) * t388) * qJD(5) + t423 * t124 + t424 * t208 + ((-t35 * t393 - t36 * t388) * qJD(5) + t415) * mrSges(6,3) + t530 * pkin(4); -t401 * t486 + (-t15 * t406 - t16 * t334 + t256 * t56 - t257 * t57) * mrSges(7,3) + t173 + t378 * t43 + t149 * t352 + t211 * t335 + t275 * t68 + t276 * t69 + t104 * t261 + t200 * t156 + t201 * t157 + t158 * t183 - t152 * mrSges(5,2) + t153 * mrSges(5,1) + (-t216 * t439 - t215 * t441 + m(6) * (-t140 * t439 - t141 * t441 + t408) + t393 * t144 - t388 * t143) * pkin(12) + m(7) * (t104 * t378 + t15 * t276 + t158 * t435 + t16 * t275 + t200 * t57 + t201 * t56) + t181 * t508 + t83 * t496 + t84 * t497 + t108 * t498 + t109 * t499 + t182 * t507 + t169 * t482 + t170 * t483 + t267 * t487 + t42 * t488 + t41 * t489 + t100 * t476 + t101 * t478 + (t425 + (pkin(5) * t120 - t164 / 0.2e1) * t388) * qJD(5) + t423 * t266 + t424 * t317 + ((-t140 * t393 - t141 * t388) * qJD(5) + t408) * mrSges(6,3) + t529 * pkin(4); t382 + t378 * t132 + t349 * t183 + t127 * t489 + t128 * t488 + t185 * t491 + t186 * t490 + t275 * t177 + t276 * t178 + t200 * t277 + t201 * t278 + t288 * t261 + t193 * t497 + t192 * t496 + t223 * t499 + t222 * t498 - pkin(4) * t241 + m(7) * (t179 * t201 + t180 * t200 + t275 * t91 + t276 * t90 + t288 * t378) + ((-m(6) * pkin(4) - mrSges(5,1) + t352) * qJD(4) * pkin(11) - t424) * t394 + (t500 - t356 * t442 / 0.2e1 - t218 * mrSges(6,3) + (-t295 * mrSges(6,3) - t303 / 0.2e1 + (m(7) * t349 + t232) * pkin(5)) * qJD(5) + (-qJD(5) * t343 + m(6) * (-t218 - t531) - t283) * pkin(12)) * t388 + (t179 * t256 - t180 * t257 - t334 * t91 - t406 * t90) * mrSges(7,3) + (t501 + t442 * t482 + qJD(5) * t492 + t422 * mrSges(6,3) + (m(6) * t422 - qJD(5) * t344 + t284) * pkin(12)) * t393 + (t341 * t476 + t339 * t479 + pkin(11) * t335 + (t356 * t477 + t358 * t479) * qJD(5) + (pkin(11) * mrSges(5,2) - Ifges(5,6) + t423) * qJD(4)) * t389; (t200 * t276 + t201 * t275 + t378 * t435) * t522 + 0.2e1 * t261 * t435 + 0.2e1 * t378 * t183 - t257 * t263 - t406 * t185 - t256 * t264 + t334 * t186 - 0.2e1 * pkin(4) * t335 + t388 * t341 - t356 * t441 + (qJD(5) * t358 + t339) * t393 + 0.2e1 * (-t200 * t406 - t201 * t334 + t256 * t275 - t257 * t276) * mrSges(7,3); t7 * mrSges(6,1) - t6 * mrSges(6,2) + (-t73 * t438 + t392 * t18 + m(7) * (-t12 * t438 + t13 * t437 + t2 * t387 + t3 * t392) + t72 * t437 + t387 * t19) * pkin(5) + t403 + t24; t59 * mrSges(6,1) - t58 * mrSges(6,2) + (m(7) * (t15 * t387 + t16 * t392 + t437 * t57 - t438 * t56) + t156 * t437 + t387 * t69 - t157 * t438 + t392 * t68) * pkin(5) + t402 + t99; t218 * mrSges(6,1) - t217 * mrSges(6,2) + (-t278 * t438 + t392 * t177 + m(7) * (-t179 * t438 + t180 * t437 + t387 * t90 + t392 * t91) + t277 * t437 + t387 * t178) * pkin(5) + t229 + t399; t381 + (pkin(12) * t352 - t464) * qJD(5) + (m(7) * (t200 * t387 + t201 * t392 + (-t275 * t387 + t276 * t392) * qJD(6)) + (t392 * t256 - t387 * t257 + (t334 * t387 - t392 * t406) * qJD(6)) * mrSges(7,3)) * pkin(5) + t405; 0.2e1 * t327; t403; t402; t399; t405; t327; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
