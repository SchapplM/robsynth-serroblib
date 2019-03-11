% Calculate time derivative of joint inertia matrix for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:47:21
% EndTime: 2019-03-10 01:47:43
% DurationCPUTime: 9.89s
% Computational Cost: add. (14502->724), mult. (37417->1017), div. (0->0), fcn. (36542->10), ass. (0->305)
t483 = Ifges(7,4) + Ifges(6,5);
t318 = sin(qJ(5));
t322 = cos(qJ(5));
t482 = t318 ^ 2 + t322 ^ 2;
t422 = Ifges(7,5) * t318;
t345 = Ifges(7,1) * t322 + t422;
t424 = Ifges(6,4) * t318;
t346 = Ifges(6,1) * t322 - t424;
t476 = (t345 + t346) * qJD(5);
t481 = t476 / 0.2e1;
t421 = Ifges(7,5) * t322;
t423 = Ifges(6,4) * t322;
t475 = -t421 + t423 + (Ifges(6,1) + Ifges(7,1)) * t318;
t317 = cos(pkin(6));
t320 = sin(qJ(3));
t324 = cos(qJ(3));
t316 = sin(pkin(6));
t321 = sin(qJ(2));
t406 = t316 * t321;
t250 = t317 * t324 - t320 * t406;
t251 = t317 * t320 + t324 * t406;
t319 = sin(qJ(4));
t323 = cos(qJ(4));
t188 = t250 * t319 + t251 * t323;
t325 = cos(qJ(2));
t390 = qJD(2) * t316;
t367 = t325 * t390;
t220 = -qJD(3) * t251 - t320 * t367;
t221 = qJD(3) * t250 + t324 * t367;
t110 = qJD(4) * t188 - t323 * t220 + t221 * t319;
t336 = t323 * t250 - t251 * t319;
t109 = qJD(4) * t336 + t220 * t319 + t221 * t323;
t405 = t316 * t325;
t164 = t318 * t188 + t322 * t405;
t389 = qJD(2) * t321;
t368 = t316 * t389;
t63 = -qJD(5) * t164 + t322 * t109 + t318 * t368;
t371 = t318 * t405;
t383 = qJD(5) * t322;
t64 = -qJD(5) * t371 + t109 * t318 + t188 * t383 - t322 * t368;
t18 = Ifges(7,1) * t63 + Ifges(7,4) * t110 + Ifges(7,5) * t64;
t19 = Ifges(6,1) * t63 - Ifges(6,4) * t64 + Ifges(6,5) * t110;
t480 = t18 + t19;
t165 = t322 * t188 - t371;
t77 = Ifges(7,1) * t165 - Ifges(7,4) * t336 + Ifges(7,5) * t164;
t78 = Ifges(6,1) * t165 - Ifges(6,4) * t164 - Ifges(6,5) * t336;
t479 = t77 + t78;
t265 = t319 * t324 + t320 * t323;
t461 = qJD(3) + qJD(4);
t212 = t461 * t265;
t384 = qJD(5) * t318;
t364 = t265 * t384;
t264 = t319 * t320 - t323 * t324;
t211 = t461 * t264;
t411 = t211 * t322;
t333 = t364 + t411;
t363 = t265 * t383;
t404 = t318 * t211;
t334 = t363 - t404;
t89 = -Ifges(7,1) * t333 + Ifges(7,4) * t212 + Ifges(7,5) * t334;
t90 = -Ifges(6,1) * t333 - Ifges(6,4) * t334 + Ifges(6,5) * t212;
t478 = t89 + t90;
t283 = -t322 * mrSges(6,1) + t318 * mrSges(6,2);
t477 = -mrSges(5,1) + t283;
t170 = Ifges(7,4) * t264 + t265 * t345;
t171 = Ifges(6,5) * t264 + t265 * t346;
t397 = t170 + t171;
t385 = qJD(4) * t323;
t376 = pkin(3) * t385;
t474 = t482 * t376;
t388 = qJD(3) * t320;
t378 = pkin(3) * t388;
t137 = pkin(4) * t212 + pkin(11) * t211 + t378;
t447 = -pkin(10) - pkin(9);
t369 = qJD(3) * t447;
t278 = t320 * t369;
t350 = t324 * t369;
t292 = t447 * t320;
t293 = t447 * t324;
t462 = t323 * t292 + t293 * t319;
t156 = qJD(4) * t462 + t323 * t278 + t319 * t350;
t306 = -pkin(3) * t324 - pkin(2);
t199 = pkin(4) * t264 - pkin(11) * t265 + t306;
t232 = t292 * t319 - t293 * t323;
t50 = t318 * t137 + t322 * t156 + t199 * t383 - t232 * t384;
t463 = t318 * t199 + t322 * t232;
t51 = -qJD(5) * t463 + t137 * t322 - t156 * t318;
t472 = -t51 * t318 + t322 * t50;
t257 = t317 * t321 * pkin(1) + pkin(8) * t405;
t242 = pkin(9) * t317 + t257;
t243 = (-pkin(2) * t325 - pkin(9) * t321 - pkin(1)) * t316;
t175 = -t320 * t242 + t324 * t243;
t134 = -pkin(3) * t405 - t251 * pkin(10) + t175;
t176 = t324 * t242 + t320 * t243;
t151 = pkin(10) * t250 + t176;
t386 = qJD(4) * t319;
t246 = (pkin(2) * t321 - pkin(9) * t325) * t390;
t299 = pkin(8) * t406;
t433 = pkin(1) * t325;
t256 = t317 * t433 - t299;
t247 = t256 * qJD(2);
t118 = -t176 * qJD(3) + t324 * t246 - t247 * t320;
t92 = pkin(3) * t368 - pkin(10) * t221 + t118;
t387 = qJD(3) * t324;
t117 = -t242 * t388 + t243 * t387 + t320 * t246 + t324 * t247;
t96 = pkin(10) * t220 + t117;
t25 = t134 * t385 - t151 * t386 + t319 * t92 + t323 * t96;
t22 = pkin(11) * t368 + t25;
t248 = t257 * qJD(2);
t174 = -t220 * pkin(3) + t248;
t40 = t110 * pkin(4) - t109 * pkin(11) + t174;
t84 = t319 * t134 + t323 * t151;
t72 = -pkin(11) * t405 + t84;
t241 = t299 + (-pkin(2) - t433) * t317;
t192 = -t250 * pkin(3) + t241;
t98 = -pkin(4) * t336 - t188 * pkin(11) + t192;
t6 = t322 * t22 + t318 * t40 + t98 * t383 - t384 * t72;
t427 = t318 * t98 + t322 * t72;
t7 = -qJD(5) * t427 - t22 * t318 + t322 * t40;
t471 = -t7 * t318 + t322 * t6;
t284 = -Ifges(7,3) * t322 + t422;
t287 = Ifges(6,2) * t322 + t424;
t438 = -t287 / 0.2e1;
t470 = t438 + t284 / 0.2e1;
t419 = Ifges(6,6) * t322;
t436 = -t322 / 0.2e1;
t437 = t318 / 0.2e1;
t469 = Ifges(7,6) * t436 + t419 / 0.2e1 + t483 * t437;
t272 = Ifges(7,4) * t383 + Ifges(7,6) * t384;
t310 = Ifges(6,5) * t383;
t360 = -t384 / 0.2e1;
t468 = -Ifges(6,6) * t360 - t310 / 0.2e1 - t272 / 0.2e1;
t343 = Ifges(7,3) * t318 + t421;
t270 = t343 * qJD(5);
t344 = -Ifges(6,2) * t318 + t423;
t273 = t344 * qJD(5);
t467 = t270 / 0.2e1 - t273 / 0.2e1;
t2 = qJ(6) * t110 - qJD(6) * t336 + t6;
t30 = -qJ(6) * t336 + t427;
t4 = -pkin(5) * t110 - t7;
t466 = t2 * t322 - t30 * t384 + t318 * t4;
t127 = qJ(6) * t264 + t463;
t142 = t199 * t322 - t232 * t318;
t128 = -pkin(5) * t264 - t142;
t43 = qJ(6) * t212 + qJD(6) * t264 + t50;
t45 = -pkin(5) * t212 - t51;
t465 = -t127 * t384 + t128 * t383 + t318 * t45 + t322 * t43;
t15 = Ifges(6,5) * t63 - Ifges(6,6) * t64 + Ifges(6,3) * t110;
t16 = Ifges(7,4) * t63 + Ifges(7,2) * t110 + Ifges(7,6) * t64;
t464 = t15 + t16;
t460 = Ifges(4,5) * t221 + Ifges(4,6) * t220 + Ifges(4,3) * t368;
t459 = (-mrSges(5,2) + (mrSges(7,2) + mrSges(6,3)) * t482) * t376;
t458 = 2 * m(5);
t457 = 2 * m(6);
t456 = 2 * m(7);
t455 = -2 * mrSges(3,3);
t454 = -2 * mrSges(5,3);
t157 = qJD(4) * t232 + t278 * t319 - t323 * t350;
t453 = 0.2e1 * t157;
t452 = -0.2e1 * t462;
t451 = 0.2e1 * t248;
t347 = t318 * mrSges(7,1) - t322 * mrSges(7,3);
t267 = t347 * qJD(5);
t450 = 0.2e1 * t267;
t449 = m(6) * pkin(4);
t435 = t322 / 0.2e1;
t432 = pkin(3) * t323;
t305 = -pkin(4) - t432;
t434 = m(6) * t305;
t426 = Ifges(4,4) * t320;
t425 = Ifges(4,4) * t324;
t420 = Ifges(6,6) * t318;
t418 = t247 * mrSges(3,2);
t417 = t248 * mrSges(3,1);
t412 = t157 * t462;
t410 = t462 * t319;
t409 = t265 * t318;
t408 = t265 * t322;
t403 = t318 * t323;
t402 = t322 * t323;
t104 = mrSges(6,1) * t164 + mrSges(6,2) * t165;
t173 = -mrSges(5,1) * t405 - t188 * mrSges(5,3);
t401 = t104 - t173;
t113 = mrSges(6,2) * t336 - mrSges(6,3) * t164;
t114 = -mrSges(7,2) * t164 - mrSges(7,3) * t336;
t400 = t113 + t114;
t115 = -mrSges(6,1) * t336 - mrSges(6,3) * t165;
t116 = mrSges(7,1) * t336 + mrSges(7,2) * t165;
t399 = -t115 + t116;
t166 = Ifges(7,6) * t264 + t265 * t343;
t169 = Ifges(6,6) * t264 + t265 * t344;
t398 = t166 - t169;
t396 = -Ifges(6,5) * t411 + Ifges(6,3) * t212;
t201 = -mrSges(6,2) * t264 - mrSges(6,3) * t409;
t204 = -mrSges(7,2) * t409 + mrSges(7,3) * t264;
t395 = t201 + t204;
t202 = mrSges(6,1) * t264 - mrSges(6,3) * t408;
t203 = -mrSges(7,1) * t264 + mrSges(7,2) * t408;
t394 = -t202 + t203;
t304 = pkin(3) * t319 + pkin(11);
t393 = t474 * t304;
t392 = t474 * pkin(11);
t382 = qJD(6) * t322;
t381 = 0.2e1 * t316;
t379 = m(7) * t382;
t377 = pkin(3) * t386;
t375 = Ifges(4,6) * t405;
t73 = Ifges(7,5) * t165 - Ifges(7,6) * t336 + Ifges(7,3) * t164;
t76 = Ifges(6,4) * t165 - Ifges(6,2) * t164 - Ifges(6,6) * t336;
t374 = t73 / 0.2e1 - t76 / 0.2e1;
t373 = -t78 / 0.2e1 - t77 / 0.2e1;
t308 = mrSges(7,2) * t383;
t370 = Ifges(5,5) * t109 - Ifges(5,6) * t110 + Ifges(5,3) * t368;
t359 = t384 / 0.2e1;
t358 = t383 / 0.2e1;
t33 = -t110 * mrSges(7,1) + t63 * mrSges(7,2);
t83 = t134 * t323 - t319 * t151;
t355 = Ifges(5,5) * t368;
t354 = Ifges(5,6) * t368;
t353 = -Ifges(7,4) * t411 + Ifges(7,2) * t212 + t334 * Ifges(7,6);
t352 = mrSges(7,2) * t382 + t272 + t310;
t71 = pkin(4) * t405 - t83;
t348 = t318 * mrSges(6,1) + t322 * mrSges(6,2);
t282 = -t322 * mrSges(7,1) - t318 * mrSges(7,3);
t342 = pkin(5) * t322 + qJ(6) * t318;
t341 = pkin(5) * t318 - qJ(6) * t322;
t36 = -t318 * t72 + t322 * t98;
t335 = t477 * t377;
t26 = -t134 * t386 - t151 * t385 - t319 * t96 + t323 * t92;
t279 = -pkin(4) - t342;
t249 = pkin(5) * t384 - qJ(6) * t383 - qJD(6) * t318;
t123 = -t212 * mrSges(7,1) - t333 * mrSges(7,2);
t332 = -mrSges(7,2) * t342 - t420;
t23 = -pkin(4) * t368 - t26;
t331 = (t284 - t287) * t384 + t475 * t383 + (-t270 + t273) * t322 + t476 * t318;
t330 = -m(7) * t342 + t282 + t283;
t31 = pkin(5) * t336 - t36;
t32 = mrSges(6,1) * t110 - mrSges(6,3) * t63;
t34 = -mrSges(6,2) * t110 - mrSges(6,3) * t64;
t35 = -mrSges(7,2) * t64 + mrSges(7,3) * t110;
t329 = (t34 + t35) * t322 + (-t32 + t33) * t318 + (-t318 * t400 + t322 * t399) * qJD(5) + m(7) * (t31 * t383 + t466) + m(6) * (-t36 * t383 - t384 * t427 + t471);
t14 = Ifges(7,5) * t63 + Ifges(7,6) * t110 + Ifges(7,3) * t64;
t17 = Ifges(6,4) * t63 - Ifges(6,2) * t64 + Ifges(6,6) * t110;
t268 = t348 * qJD(5);
t41 = pkin(5) * t164 - qJ(6) * t165 + t71;
t9 = pkin(5) * t64 - qJ(6) * t63 - qJD(6) * t165 + t23;
t328 = t26 * mrSges(5,1) - t25 * mrSges(5,2) + t14 * t436 + t17 * t435 + t23 * t283 + t41 * t267 + t71 * t268 + t9 * t282 + t31 * t308 + t73 * t359 + t76 * t360 + t370 + t470 * t64 + t475 * t63 / 0.2e1 + t165 * t481 + t480 * t437 + t479 * t358 + t468 * t336 + t467 * t164 + t469 * t110 + ((-t318 * t427 - t322 * t36) * qJD(5) + t471) * mrSges(6,3) + t466 * mrSges(7,2);
t122 = mrSges(6,1) * t212 + mrSges(6,3) * t333;
t124 = -mrSges(6,2) * t212 - mrSges(6,3) * t334;
t125 = -mrSges(7,2) * t334 + mrSges(7,3) * t212;
t327 = (t124 + t125) * t322 + (-t122 + t123) * t318 + (-t318 * t395 + t322 * t394) * qJD(5) + m(7) * t465 + m(6) * (-t142 * t383 - t384 * t463 + t472);
t155 = t265 * t341 - t462;
t206 = Ifges(5,6) * t212;
t208 = Ifges(5,5) * t211;
t57 = -t341 * t211 + (qJD(5) * t342 - t382) * t265 + t157;
t85 = -Ifges(7,5) * t333 + Ifges(7,6) * t212 + Ifges(7,3) * t334;
t88 = -Ifges(6,4) * t333 - Ifges(6,2) * t334 + Ifges(6,6) * t212;
t326 = -t156 * mrSges(5,2) + t155 * t267 + t166 * t359 + t169 * t360 - t462 * t268 + t57 * t282 + t363 * t438 + t88 * t435 + t85 * t436 - t206 - t208 + t478 * t437 + t467 * t409 - t470 * t404 + t408 * t481 - t468 * t264 + t469 * t212 + t477 * t157 + ((-t142 * t322 - t318 * t463) * qJD(5) + t472) * mrSges(6,3) + (t265 * t284 + t397) * t358 + t465 * mrSges(7,2) + t475 * (t265 * t360 - t411 / 0.2e1);
t311 = Ifges(4,5) * t387;
t298 = Ifges(3,5) * t367;
t291 = Ifges(4,1) * t320 + t425;
t288 = Ifges(4,2) * t324 + t426;
t277 = (Ifges(4,1) * t324 - t426) * qJD(3);
t274 = (-Ifges(4,2) * t320 + t425) * qJD(3);
t269 = (mrSges(4,1) * t320 + mrSges(4,2) * t324) * qJD(3);
t258 = t279 - t432;
t237 = t249 + t377;
t230 = -mrSges(4,1) * t405 - t251 * mrSges(4,3);
t229 = mrSges(4,2) * t405 + t250 * mrSges(4,3);
t219 = Ifges(5,1) * t265 - Ifges(5,4) * t264;
t218 = Ifges(5,4) * t265 - Ifges(5,2) * t264;
t217 = mrSges(5,1) * t264 + mrSges(5,2) * t265;
t194 = t348 * t265;
t193 = t347 * t265;
t186 = mrSges(4,1) * t368 - mrSges(4,3) * t221;
t185 = -mrSges(4,2) * t368 + mrSges(4,3) * t220;
t182 = Ifges(4,1) * t251 + Ifges(4,4) * t250 - Ifges(4,5) * t405;
t181 = Ifges(4,4) * t251 + Ifges(4,2) * t250 - t375;
t172 = mrSges(5,2) * t405 + mrSges(5,3) * t336;
t168 = Ifges(7,2) * t264 + (Ifges(7,4) * t322 + Ifges(7,6) * t318) * t265;
t167 = Ifges(6,3) * t264 + (Ifges(6,5) * t322 - t420) * t265;
t152 = -mrSges(4,1) * t220 + mrSges(4,2) * t221;
t150 = -Ifges(5,1) * t211 - Ifges(5,4) * t212;
t149 = -Ifges(5,4) * t211 - Ifges(5,2) * t212;
t148 = mrSges(5,1) * t212 - mrSges(5,2) * t211;
t136 = Ifges(4,1) * t221 + Ifges(4,4) * t220 + Ifges(4,5) * t368;
t135 = Ifges(4,4) * t221 + Ifges(4,2) * t220 + Ifges(4,6) * t368;
t126 = -mrSges(5,1) * t336 + mrSges(5,2) * t188;
t120 = Ifges(5,1) * t188 + Ifges(5,4) * t336 - Ifges(5,5) * t405;
t119 = Ifges(5,4) * t188 + Ifges(5,2) * t336 - Ifges(5,6) * t405;
t112 = mrSges(6,1) * t334 - mrSges(6,2) * t333;
t111 = mrSges(7,1) * t334 + mrSges(7,3) * t333;
t103 = mrSges(7,1) * t164 - mrSges(7,3) * t165;
t102 = -mrSges(5,2) * t368 - mrSges(5,3) * t110;
t101 = mrSges(5,1) * t368 - mrSges(5,3) * t109;
t87 = -Ifges(7,4) * t364 + t353;
t86 = -Ifges(6,5) * t364 - Ifges(6,6) * t334 + t396;
t75 = Ifges(7,4) * t165 - Ifges(7,2) * t336 + Ifges(7,6) * t164;
t74 = Ifges(6,5) * t165 - Ifges(6,6) * t164 - Ifges(6,3) * t336;
t48 = mrSges(5,1) * t110 + mrSges(5,2) * t109;
t47 = Ifges(5,1) * t109 - Ifges(5,4) * t110 + t355;
t46 = Ifges(5,4) * t109 - Ifges(5,2) * t110 + t354;
t28 = mrSges(6,1) * t64 + mrSges(6,2) * t63;
t27 = mrSges(7,1) * t64 - mrSges(7,3) * t63;
t1 = [t479 * t63 + t480 * t165 + 0.2e1 * m(3) * (t247 * t257 - t248 * t256) + 0.2e1 * m(4) * (t117 * t176 + t118 * t175 + t241 * t248) + (t73 - t76) * t64 + (t23 * t71 + t36 * t7 + t427 * t6) * t457 + 0.2e1 * t427 * t34 + (mrSges(3,3) * t321 * t451 + (0.2e1 * mrSges(3,3) * t247 - t370 - t460) * t325 + ((t256 * t455 + Ifges(3,5) * t317 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t325) * t381) * t325 + (t257 * t455 + Ifges(4,5) * t251 + Ifges(5,5) * t188 - 0.2e1 * Ifges(3,6) * t317 + Ifges(4,6) * t250 + Ifges(5,6) * t336 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t321) * t381 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3)) * t405) * t321) * qJD(2)) * t316 - (-t46 + t464) * t336 + (t298 - 0.2e1 * t417 - 0.2e1 * t418) * t317 + (t14 - t17) * t164 + (-mrSges(4,1) * t250 + mrSges(4,2) * t251) * t451 + (t2 * t30 + t31 * t4 + t41 * t9) * t456 + (t174 * t192 + t25 * t84 + t26 * t83) * t458 + (-t119 + t74 + t75) * t110 + t250 * t135 + t251 * t136 + 0.2e1 * t117 * t229 + 0.2e1 * t118 * t230 + 0.2e1 * t241 * t152 + t220 * t181 + t221 * t182 + t188 * t47 + 0.2e1 * t192 * t48 + 0.2e1 * t176 * t185 + 0.2e1 * t175 * t186 + 0.2e1 * t25 * t172 + 0.2e1 * t26 * t173 + 0.2e1 * t174 * t126 + 0.2e1 * t6 * t113 + 0.2e1 * t2 * t114 + 0.2e1 * t7 * t115 + 0.2e1 * t4 * t116 + t109 * t120 + 0.2e1 * t84 * t102 + 0.2e1 * t9 * t103 + 0.2e1 * t23 * t104 + 0.2e1 * t83 * t101 + 0.2e1 * t71 * t28 + 0.2e1 * t41 * t27 + 0.2e1 * t36 * t32 + 0.2e1 * t31 * t33 + 0.2e1 * t30 * t35; m(6) * (t142 * t7 + t157 * t71 - t23 * t462 + t36 * t51 + t427 * t50 + t463 * t6) + m(5) * (t156 * t84 - t157 * t83 + t174 * t306 + t192 * t378 + t232 * t25 + t26 * t462) - (t28 - t101) * t462 + ((t208 / 0.2e1 + t206 / 0.2e1 - t311 / 0.2e1) * t325 + (Ifges(4,5) * t320 / 0.2e1 + Ifges(4,6) * t324 / 0.2e1 - Ifges(3,6)) * t389) * t316 + t298 + (-t119 / 0.2e1 + t74 / 0.2e1 + t75 / 0.2e1 - t84 * mrSges(5,3)) * t212 + (-t248 * mrSges(4,1) + t135 / 0.2e1 + pkin(9) * t185 + t117 * mrSges(4,3)) * t324 + (-t218 / 0.2e1 + t167 / 0.2e1 + t168 / 0.2e1) * t110 + t427 * t124 + t463 * t34 - (t86 / 0.2e1 + t87 / 0.2e1 - t149 / 0.2e1) * t336 + t401 * t157 + (t355 / 0.2e1 + t47 / 0.2e1 - t26 * mrSges(5,3) + (t18 / 0.2e1 + t19 / 0.2e1) * t322 + (t14 / 0.2e1 - t17 / 0.2e1) * t318 + (t318 * t373 + t322 * t374) * qJD(5)) * t265 - (t120 / 0.2e1 - t83 * mrSges(5,3) - t373 * t322 + t374 * t318) * t211 + ((-pkin(9) * t230 - t175 * mrSges(4,3) + t182 / 0.2e1) * t324 + (-pkin(9) * t229 - t176 * mrSges(4,3) + pkin(3) * t126 + t375 / 0.2e1 - t181 / 0.2e1) * t320) * qJD(3) + (-t354 / 0.2e1 + t16 / 0.2e1 + t15 / 0.2e1 - t46 / 0.2e1 - t25 * mrSges(5,3)) * t264 + (t85 / 0.2e1 - t88 / 0.2e1) * t164 - t417 - t418 + t306 * t48 + t220 * t288 / 0.2e1 + t221 * t291 / 0.2e1 + t241 * t269 + t250 * t274 / 0.2e1 + t251 * t277 / 0.2e1 + t232 * t102 + t174 * t217 + t109 * t219 / 0.2e1 + t9 * t193 + t23 * t194 + t6 * t201 + t7 * t202 + t4 * t203 + t2 * t204 + t188 * t150 / 0.2e1 + t192 * t148 + t156 * t172 - pkin(2) * t152 + t155 * t27 + t142 * t32 + t36 * t122 + t31 * t123 + t30 * t125 + t127 * t35 + t128 * t33 + t50 * t113 + t43 * t114 + t51 * t115 + t45 * t116 + t41 * t111 + t71 * t112 + t57 * t103 + m(7) * (t127 * t2 + t128 * t4 + t155 * t9 + t30 * t43 + t31 * t45 + t41 * t57) + ((t117 * t324 - t118 * t320 + (-t175 * t324 - t176 * t320) * qJD(3)) * pkin(9) - pkin(2) * t248) * m(4) + (t248 * mrSges(4,2) + t136 / 0.2e1 - pkin(9) * t186 - t118 * mrSges(4,3)) * t320 + (t166 / 0.2e1 - t169 / 0.2e1) * t64 + (t170 / 0.2e1 + t171 / 0.2e1) * t63 + (t89 / 0.2e1 + t90 / 0.2e1) * t165; t324 * t274 + t320 * t277 + 0.2e1 * t306 * t148 - 0.2e1 * pkin(2) * t269 + t112 * t452 + 0.2e1 * t57 * t193 + t194 * t453 + 0.2e1 * t50 * t201 + 0.2e1 * t51 * t202 + 0.2e1 * t45 * t203 + 0.2e1 * t43 * t204 + 0.2e1 * t155 * t111 + 0.2e1 * t142 * t122 + 0.2e1 * t463 * t124 + 0.2e1 * t127 * t125 + 0.2e1 * t128 * t123 + (t324 * t291 + (0.2e1 * pkin(3) * t217 - t288) * t320) * qJD(3) + (t156 * t232 + t306 * t378 - t412) * t458 + (t142 * t51 + t463 * t50 - t412) * t457 + (t127 * t43 + t128 * t45 + t155 * t57) * t456 + (t156 * t454 - t149 + t86 + t87) * t264 + (t232 * t454 + t167 + t168 - t218) * t212 - (mrSges(5,3) * t452 + t318 * t398 + t322 * t397 + t219) * t211 + (mrSges(5,3) * t453 + t150 + t478 * t322 + (t85 - t88) * t318 + (-t318 * t397 + t322 * t398) * qJD(5)) * t265; t23 * t434 + t329 * t304 + (m(5) * (t25 * t319 + t26 * t323) + t323 * t101 + t319 * t102 + (t401 * t319 + (t318 * t399 + t322 * t400 + t172) * t323 + m(7) * (t30 * t402 + t31 * t403) + m(6) * (t319 * t71 - t36 * t403 + t402 * t427) + m(5) * (-t319 * t83 + t323 * t84)) * qJD(4)) * pkin(3) + t328 + t305 * t28 + t258 * t27 + t237 * t103 - t117 * mrSges(4,2) + t118 * mrSges(4,1) + m(7) * (t237 * t41 + t258 * t9) + t460; t311 + (-Ifges(4,6) * t320 + (-mrSges(4,1) * t324 + mrSges(4,2) * t320) * pkin(9)) * qJD(3) + m(7) * (t155 * t237 + t258 * t57) + t327 * t304 + (m(5) * (t156 * t319 - t157 * t323) + (t211 * t323 - t212 * t319) * mrSges(5,3) + ((mrSges(5,3) * t265 + t194) * t319 + (-t264 * mrSges(5,3) + t318 * t394 + t322 * t395) * t323 + m(7) * (t127 * t402 + t128 * t403) + m(6) * (-t142 * t403 + t402 * t463 - t410) + m(5) * (t232 * t323 - t410)) * qJD(4)) * pkin(3) + t326 + t305 * t112 + t258 * t111 + t237 * t193 + t157 * t434; 0.2e1 * t237 * t282 + t258 * t450 + 0.2e1 * t305 * t268 + 0.2e1 * t335 + (t237 * t258 + t393) * t456 + (t305 * t377 + t393) * t457 + 0.2e1 * t459 + t331; t329 * pkin(11) + t328 + t279 * t27 + t249 * t103 - pkin(4) * t28 + m(7) * (t249 * t41 + t279 * t9) - t23 * t449; t327 * pkin(11) + t326 + t279 * t111 + t249 * t193 + m(7) * (t155 * t249 + t279 * t57) - pkin(4) * t112 - t157 * t449; (t237 + t249) * t282 + (t305 - pkin(4)) * t268 + (t258 + t279) * t267 + t335 + m(7) * (t237 * t279 + t249 * t258 + t392) + m(6) * (-pkin(4) * t377 + t392) + t459 + t331; -0.2e1 * pkin(4) * t268 + t279 * t450 + 0.2e1 * (m(7) * t279 + t282) * t249 + t331; -pkin(5) * t33 + m(7) * (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t30) + qJD(6) * t114 + qJ(6) * t35 + t2 * mrSges(7,3) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t7 * mrSges(6,1) + t464; m(7) * (-pkin(5) * t45 + qJ(6) * t43 + qJD(6) * t127) + t43 * mrSges(7,3) + qJD(6) * t204 + qJ(6) * t125 + Ifges(6,6) * t404 - t50 * mrSges(6,2) - t45 * mrSges(7,1) + t51 * mrSges(6,1) - pkin(5) * t123 + (-t318 * t483 - t419) * t265 * qJD(5) + t353 + t396; t304 * t379 + (-m(7) * t341 - t347 - t348) * t376 + (t304 * t330 + t332) * qJD(5) + t352; t332 * qJD(5) + (qJD(5) * t330 + t379) * pkin(11) + t352; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t4 + t33; m(7) * t45 + t123; t308 + (t304 * t383 + t318 * t376) * m(7); m(7) * pkin(11) * t383 + t308; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
