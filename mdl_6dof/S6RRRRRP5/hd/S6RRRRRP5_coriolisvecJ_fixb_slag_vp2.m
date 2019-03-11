% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:21
% EndTime: 2019-03-10 01:19:07
% DurationCPUTime: 25.91s
% Computational Cost: add. (19831->745), mult. (48784->1006), div. (0->0), fcn. (34932->8), ass. (0->330)
t327 = sin(qJ(2));
t331 = cos(qJ(2));
t351 = pkin(2) * t327 - pkin(8) * t331;
t286 = t351 * qJD(1);
t330 = cos(qJ(3));
t326 = sin(qJ(3));
t384 = qJD(1) * t327;
t364 = t326 * t384;
t241 = pkin(7) * t364 + t330 * t286;
t390 = t330 * t331;
t340 = pkin(3) * t327 - pkin(9) * t390;
t460 = -pkin(9) - pkin(8);
t365 = qJD(3) * t460;
t550 = -qJD(1) * t340 + t330 * t365 - t241;
t265 = t326 * t286;
t391 = t327 * t330;
t392 = t326 * t331;
t549 = t265 + (-pkin(7) * t391 - pkin(9) * t392) * qJD(1) - t326 * t365;
t325 = sin(qJ(4));
t329 = cos(qJ(4));
t341 = t325 * t326 - t329 * t330;
t467 = qJD(3) + qJD(4);
t231 = t467 * t341;
t338 = t341 * t331;
t250 = qJD(1) * t338;
t548 = t231 - t250;
t282 = t325 * t330 + t326 * t329;
t232 = t467 * t282;
t339 = t282 * t331;
t249 = qJD(1) * t339;
t540 = t232 - t249;
t381 = qJD(2) * t330;
t279 = -t364 + t381;
t363 = t330 * t384;
t280 = qJD(2) * t326 + t363;
t222 = t279 * t325 + t280 * t329;
t324 = sin(qJ(5));
t328 = cos(qJ(5));
t354 = t329 * t279 - t280 * t325;
t157 = t222 * t328 + t324 * t354;
t447 = -t157 / 0.2e1;
t382 = qJD(2) * t327;
t356 = qJD(1) * t382;
t373 = qJD(2) * qJD(3);
t379 = qJD(3) * t326;
t380 = qJD(2) * t331;
t239 = t330 * t373 + (-t327 * t379 + t330 * t380) * qJD(1);
t378 = qJD(3) * t330;
t507 = t326 * t380 + t327 * t378;
t240 = -qJD(1) * t507 - t326 * t373;
t128 = qJD(4) * t354 + t239 * t329 + t240 * t325;
t129 = -qJD(4) * t222 - t239 * t325 + t240 * t329;
t504 = -t222 * t324 + t328 * t354;
t58 = qJD(5) * t504 + t128 * t328 + t129 * t324;
t291 = -pkin(2) * t331 - t327 * pkin(8) - pkin(1);
t271 = t291 * qJD(1);
t383 = qJD(1) * t331;
t320 = pkin(7) * t383;
t297 = qJD(2) * pkin(8) + t320;
t227 = t330 * t271 - t297 * t326;
t188 = -pkin(9) * t280 + t227;
t312 = qJD(3) - t383;
t177 = pkin(3) * t312 + t188;
t228 = t271 * t326 + t297 * t330;
t189 = pkin(9) * t279 + t228;
t183 = t329 * t189;
t106 = t177 * t325 + t183;
t289 = t351 * qJD(2);
t272 = qJD(1) * t289;
t353 = pkin(7) * t356;
t167 = -qJD(3) * t228 + t330 * t272 + t326 * t353;
t119 = pkin(3) * t356 - pkin(9) * t239 + t167;
t166 = t271 * t378 + t326 * t272 - t297 * t379 - t330 * t353;
t134 = pkin(9) * t240 + t166;
t36 = -qJD(4) * t106 + t329 * t119 - t134 * t325;
t21 = pkin(4) * t356 - pkin(10) * t128 + t36;
t376 = qJD(4) * t329;
t377 = qJD(4) * t325;
t35 = t325 * t119 + t329 * t134 + t177 * t376 - t189 * t377;
t24 = pkin(10) * t129 + t35;
t301 = qJD(4) + t312;
t181 = t325 * t189;
t105 = t329 * t177 - t181;
t519 = pkin(10) * t222;
t89 = t105 - t519;
t83 = pkin(4) * t301 + t89;
t498 = pkin(10) * t354;
t90 = t106 + t498;
t88 = t328 * t90;
t34 = t324 * t83 + t88;
t6 = -qJD(5) * t34 + t328 * t21 - t24 * t324;
t2 = pkin(5) * t356 - qJ(6) * t58 - qJD(6) * t157 + t6;
t374 = qJD(5) * t328;
t375 = qJD(5) * t324;
t5 = t324 * t21 + t328 * t24 + t83 * t374 - t375 * t90;
t59 = -qJD(5) * t157 - t128 * t324 + t129 * t328;
t3 = qJ(6) * t59 + qJD(6) * t504 + t5;
t466 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t493 = -Ifges(7,6) - Ifges(6,6);
t495 = -Ifges(7,5) - Ifges(6,5);
t529 = Ifges(6,3) + Ifges(7,3);
t470 = t356 * t529 - t493 * t59 - t495 * t58;
t497 = Ifges(6,1) + Ifges(7,1);
t496 = Ifges(6,4) + Ifges(7,4);
t545 = t157 * t496;
t547 = t466 + t470 + (t497 * t504 - t545) * t447;
t296 = -qJD(2) * pkin(2) + pkin(7) * t384;
t245 = -pkin(3) * t279 + t296;
t175 = -pkin(4) * t354 + t245;
t86 = t324 * t90;
t33 = t328 * t83 - t86;
t511 = qJ(6) * t157;
t26 = t33 - t511;
t292 = qJD(5) + t301;
t25 = pkin(5) * t292 + t26;
t474 = qJ(6) * t504;
t27 = t34 + t474;
t430 = -t292 / 0.2e1;
t524 = t496 * t504;
t491 = t157 * t497 - t292 * t495 + t524;
t494 = Ifges(6,2) + Ifges(7,2);
t499 = -t504 / 0.2e1;
t95 = -pkin(5) * t504 + qJD(6) + t175;
t546 = (t157 * t27 + t25 * t504) * mrSges(7,3) + (t157 * t34 + t33 * t504) * mrSges(6,3) - t175 * (mrSges(6,1) * t157 + mrSges(6,2) * t504) - t95 * (mrSges(7,1) * t157 + mrSges(7,2) * t504) + (t157 * t493 - t495 * t504) * t430 + (-t157 * t494 + t491 + t524) * t499;
t298 = t460 * t326;
t299 = t460 * t330;
t510 = t298 * t376 + t299 * t377 + t325 * t550 - t549 * t329;
t238 = t325 * t298 - t329 * t299;
t509 = -qJD(4) * t238 + t549 * t325 + t329 * t550;
t215 = Ifges(5,4) * t354;
t146 = Ifges(5,1) * t222 + Ifges(5,5) * t301 + t215;
t367 = Ifges(5,5) * t128 + Ifges(5,6) * t129 + Ifges(5,3) * t356;
t414 = Ifges(5,4) * t222;
t428 = -t301 / 0.2e1;
t441 = -t222 / 0.2e1;
t443 = -t354 / 0.2e1;
t475 = t36 * mrSges(5,1) - t35 * mrSges(5,2);
t544 = t367 + t475 + (Ifges(5,5) * t354 - Ifges(5,6) * t222) * t428 + (t105 * t354 + t106 * t222) * mrSges(5,3) + (-Ifges(5,2) * t222 + t146 + t215) * t443 - t245 * (mrSges(5,1) * t222 + mrSges(5,2) * t354) + (Ifges(5,1) * t354 - t414) * t441 + t547;
t543 = -t356 / 0.2e1;
t542 = -pkin(4) * t384 + pkin(10) * t548 + t509;
t541 = pkin(10) * t540 - t510;
t464 = t58 / 0.2e1;
t463 = t59 / 0.2e1;
t538 = pkin(5) * t157;
t223 = -t282 * t324 - t328 * t341;
t109 = qJD(5) * t223 - t231 * t328 - t232 * t324;
t185 = -t249 * t324 - t250 * t328;
t389 = t109 - t185;
t224 = t282 * t328 - t324 * t341;
t110 = -qJD(5) * t224 + t231 * t324 - t232 * t328;
t184 = -t249 * t328 + t250 * t324;
t388 = t110 - t184;
t422 = pkin(3) * t326;
t274 = t383 * t422 + t320;
t535 = pkin(3) * t379 - t274;
t360 = Ifges(3,5) * qJD(2) / 0.2e1;
t492 = -t292 * t493 + t494 * t504 + t545;
t532 = t492 / 0.2e1;
t531 = t496 * t463 + t497 * t464 + t495 * t543;
t530 = t494 * t463 + t496 * t464 + t493 * t543;
t237 = t329 * t298 + t299 * t325;
t200 = -pkin(10) * t282 + t237;
t201 = -pkin(10) * t341 + t238;
t141 = t324 * t200 + t328 * t201;
t514 = -qJD(5) * t141 + t324 * t541 + t328 * t542;
t513 = t200 * t374 - t201 * t375 + t324 * t542 - t328 * t541;
t508 = pkin(4) * t540 + t535;
t446 = t157 / 0.2e1;
t521 = t492 * t446 + t546;
t520 = pkin(4) * t222;
t516 = -pkin(5) * t384 - qJ(6) * t389 - qJD(6) * t224 + t514;
t515 = qJ(6) * t388 + qJD(6) * t223 + t513;
t512 = -pkin(5) * t388 + t508;
t316 = pkin(3) * t329 + pkin(4);
t395 = t324 * t325;
t114 = -t188 * t325 - t183;
t96 = t114 - t498;
t115 = t329 * t188 - t181;
t97 = t115 - t519;
t478 = -t324 * t96 - t328 * t97 + t316 * t374 + (-t325 * t375 + (t328 * t329 - t395) * qJD(4)) * pkin(3);
t394 = t325 * t328;
t476 = t324 * t97 - t328 * t96 - t316 * t375 + (-t325 * t374 + (-t324 * t329 - t394) * qJD(4)) * pkin(3);
t318 = Ifges(3,4) * t383;
t404 = t280 * Ifges(4,4);
t203 = t279 * Ifges(4,2) + t312 * Ifges(4,6) + t404;
t275 = Ifges(4,4) * t279;
t204 = t280 * Ifges(4,1) + t312 * Ifges(4,5) + t275;
t342 = t227 * t330 + t228 * t326;
t415 = Ifges(4,4) * t330;
t346 = -Ifges(4,2) * t326 + t415;
t416 = Ifges(4,4) * t326;
t348 = Ifges(4,1) * t330 - t416;
t349 = mrSges(4,1) * t326 + mrSges(4,2) * t330;
t410 = Ifges(4,6) * t326;
t411 = Ifges(4,5) * t330;
t424 = t330 / 0.2e1;
t425 = -t326 / 0.2e1;
t431 = t280 / 0.2e1;
t333 = -t342 * mrSges(4,3) + t296 * t349 + t279 * t346 / 0.2e1 + t348 * t431 + t312 * (-t410 + t411) / 0.2e1 + t203 * t425 + t204 * t424;
t505 = t333 + Ifges(3,1) * t384 / 0.2e1 + t318 / 0.2e1 + t360;
t145 = Ifges(5,2) * t354 + Ifges(5,6) * t301 + t414;
t500 = t145 / 0.2e1;
t429 = t292 / 0.2e1;
t450 = t504 / 0.2e1;
t359 = -Ifges(3,6) * qJD(2) / 0.2e1;
t479 = t478 + t511;
t477 = t476 + t474;
t257 = t341 * t327;
t278 = t330 * t291;
t421 = pkin(7) * t326;
t226 = -pkin(9) * t391 + t278 + (-pkin(3) - t421) * t331;
t314 = pkin(7) * t390;
t248 = t326 * t291 + t314;
t393 = t326 * t327;
t234 = -pkin(9) * t393 + t248;
t162 = t329 * t226 - t325 * t234;
t132 = -pkin(4) * t331 + t257 * pkin(10) + t162;
t163 = t325 * t226 + t329 * t234;
t256 = t282 * t327;
t142 = -pkin(10) * t256 + t163;
t72 = t324 * t132 + t328 * t142;
t471 = Ifges(4,5) * t239 + Ifges(4,6) * t240;
t469 = -t167 * mrSges(4,1) + t166 * mrSges(4,2);
t468 = pkin(1) * mrSges(3,2) * qJD(1);
t368 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t369 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t370 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t417 = Ifges(3,4) * t327;
t465 = t370 * t157 + t368 * t292 + t369 * t504 + t301 * Ifges(5,3) + t280 * Ifges(4,5) + t312 * Ifges(4,3) + t354 * Ifges(5,6) + t227 * mrSges(4,1) - t228 * mrSges(4,2) + t105 * mrSges(5,1) - t106 * mrSges(5,2) - t27 * mrSges(7,2) - t34 * mrSges(6,2) + t25 * mrSges(7,1) + t33 * mrSges(6,1) + t359 - (t331 * Ifges(3,2) + t417) * qJD(1) / 0.2e1 + t222 * Ifges(5,5) + t279 * Ifges(4,6) - t493 * t450 - t495 * t446 + t529 * t429;
t459 = pkin(1) * mrSges(3,1);
t455 = t128 / 0.2e1;
t454 = t129 / 0.2e1;
t442 = t354 / 0.2e1;
t440 = t222 / 0.2e1;
t437 = t239 / 0.2e1;
t436 = t240 / 0.2e1;
t435 = -t256 / 0.2e1;
t434 = -t257 / 0.2e1;
t433 = -t279 / 0.2e1;
t432 = -t280 / 0.2e1;
t427 = t301 / 0.2e1;
t426 = -t312 / 0.2e1;
t49 = -mrSges(7,2) * t356 + mrSges(7,3) * t59;
t50 = -mrSges(6,2) * t356 + mrSges(6,3) * t59;
t418 = t49 + t50;
t38 = t328 * t89 - t86;
t398 = qJD(2) * mrSges(3,2);
t385 = t330 * t289 + t382 * t421;
t290 = pkin(3) * t393 + t327 * pkin(7);
t366 = Ifges(4,3) * t356 + t471;
t246 = pkin(3) * t507 + pkin(7) * t380;
t317 = -pkin(3) * t330 - pkin(2);
t15 = -t59 * mrSges(7,1) + t58 * mrSges(7,2);
t218 = -pkin(3) * t240 + qJD(2) * t320;
t37 = -t324 * t89 - t88;
t71 = t328 * t132 - t324 * t142;
t140 = t328 * t200 - t201 * t324;
t259 = -pkin(3) * t395 + t328 * t316;
t352 = m(4) * t296 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t279 + mrSges(4,2) * t280 + mrSges(3,3) * t384;
t229 = pkin(4) * t256 + t290;
t180 = pkin(3) * t280 + t520;
t350 = mrSges(4,1) * t330 - mrSges(4,2) * t326;
t347 = Ifges(4,1) * t326 + t415;
t345 = Ifges(4,2) * t330 + t416;
t344 = Ifges(4,5) * t326 + Ifges(4,6) * t330;
t343 = t166 * t330 - t167 * t326;
t193 = -t256 * t328 + t257 * t324;
t194 = -t256 * t324 - t257 * t328;
t251 = pkin(4) * t341 + t317;
t169 = -qJD(2) * t339 + t257 * t467;
t143 = -pkin(4) * t169 + t246;
t98 = -pkin(4) * t129 + t218;
t168 = -qJD(2) * t338 - t232 * t327;
t161 = t340 * qJD(2) + (-t314 + (pkin(9) * t327 - t291) * t326) * qJD(3) + t385;
t186 = t326 * t289 + t291 * t378 + (-t327 * t381 - t331 * t379) * pkin(7);
t165 = -pkin(9) * t507 + t186;
t67 = -qJD(4) * t163 + t329 * t161 - t165 * t325;
t43 = pkin(4) * t382 - pkin(10) * t168 + t67;
t66 = t325 * t161 + t329 * t165 + t226 * t376 - t234 * t377;
t51 = pkin(10) * t169 + t66;
t9 = t132 * t374 - t142 * t375 + t324 * t43 + t328 * t51;
t10 = -qJD(5) * t72 - t324 * t51 + t328 * t43;
t315 = pkin(4) * t328 + pkin(5);
t294 = mrSges(3,3) * t383 - t398;
t260 = pkin(3) * t394 + t316 * t324;
t258 = pkin(5) + t259;
t247 = -pkin(7) * t392 + t278;
t244 = mrSges(4,1) * t312 - mrSges(4,3) * t280;
t243 = -mrSges(4,2) * t312 + mrSges(4,3) * t279;
t242 = -pkin(7) * t363 + t265;
t217 = -mrSges(4,2) * t356 + mrSges(4,3) * t240;
t216 = mrSges(4,1) * t356 - mrSges(4,3) * t239;
t191 = mrSges(5,1) * t301 - mrSges(5,3) * t222;
t190 = -mrSges(5,2) * t301 + mrSges(5,3) * t354;
t187 = -qJD(3) * t248 + t385;
t179 = -mrSges(4,1) * t240 + mrSges(4,2) * t239;
t178 = -pkin(5) * t223 + t251;
t171 = t239 * Ifges(4,1) + t240 * Ifges(4,4) + Ifges(4,5) * t356;
t170 = t239 * Ifges(4,4) + t240 * Ifges(4,2) + Ifges(4,6) * t356;
t160 = -mrSges(5,1) * t354 + mrSges(5,2) * t222;
t149 = -pkin(5) * t193 + t229;
t139 = mrSges(6,1) * t292 - mrSges(6,3) * t157;
t138 = mrSges(7,1) * t292 - mrSges(7,3) * t157;
t137 = -mrSges(6,2) * t292 + mrSges(6,3) * t504;
t136 = -mrSges(7,2) * t292 + mrSges(7,3) * t504;
t112 = -mrSges(5,2) * t356 + mrSges(5,3) * t129;
t111 = mrSges(5,1) * t356 - mrSges(5,3) * t128;
t102 = t520 + t538;
t101 = qJ(6) * t223 + t141;
t100 = -qJ(6) * t224 + t140;
t99 = t180 + t538;
t85 = -mrSges(6,1) * t504 + mrSges(6,2) * t157;
t84 = -mrSges(7,1) * t504 + mrSges(7,2) * t157;
t74 = -qJD(5) * t194 - t168 * t324 + t169 * t328;
t73 = qJD(5) * t193 + t168 * t328 + t169 * t324;
t70 = -mrSges(5,1) * t129 + mrSges(5,2) * t128;
t69 = t128 * Ifges(5,1) + t129 * Ifges(5,4) + Ifges(5,5) * t356;
t68 = t128 * Ifges(5,4) + t129 * Ifges(5,2) + Ifges(5,6) * t356;
t63 = qJ(6) * t193 + t72;
t62 = -pkin(5) * t331 - t194 * qJ(6) + t71;
t61 = -pkin(5) * t74 + t143;
t48 = mrSges(6,1) * t356 - mrSges(6,3) * t58;
t47 = mrSges(7,1) * t356 - mrSges(7,3) * t58;
t29 = t38 - t511;
t28 = t37 - t474;
t23 = -pkin(5) * t59 + t98;
t16 = -mrSges(6,1) * t59 + mrSges(6,2) * t58;
t8 = qJ(6) * t74 + qJD(6) * t193 + t9;
t7 = pkin(5) * t382 - qJ(6) * t73 - qJD(6) * t194 + t10;
t1 = [m(5) * (t105 * t67 + t106 * t66 + t162 * t36 + t163 * t35 + t218 * t290 + t245 * t246) + m(6) * (t10 * t33 + t143 * t175 + t229 * t98 + t34 * t9 + t5 * t72 + t6 * t71) + m(7) * (t149 * t23 + t2 * t62 + t25 * t7 + t27 * t8 + t3 * t63 + t61 * t95) + (t193 * t5 - t194 * t6 - t33 * t73 + t34 * t74) * mrSges(6,3) + (-t105 * t168 + t106 * t169 - t256 * t35 + t257 * t36) * mrSges(5,3) + (-Ifges(5,1) * t257 - Ifges(5,4) * t256) * t455 + (-Ifges(5,4) * t257 - Ifges(5,2) * t256) * t454 + t218 * (mrSges(5,1) * t256 - mrSges(5,2) * t257) + (t193 * t3 - t194 * t2 - t25 * t73 + t27 * t74) * mrSges(7,3) + t69 * t434 + t68 * t435 + (Ifges(5,1) * t168 + Ifges(5,4) * t169) * t440 + (Ifges(5,4) * t168 + Ifges(5,2) * t169) * t442 + m(4) * (t166 * t248 + t167 * t247 + t228 * t186 + t227 * t187) + (t193 * t496 + t194 * t497) * t464 + (t496 * t74 + t497 * t73) * t446 + (-t493 * t74 - t495 * t73) * t429 + (t193 * t494 + t194 * t496) * t463 + (t494 * t74 + t496 * t73) * t450 + (-Ifges(5,5) * t455 - Ifges(5,6) * t454 + t495 * t464 + t493 * t463 + (0.3e1 / 0.2e1 * Ifges(3,4) * t380 + (-Ifges(5,3) / 0.2e1 - Ifges(4,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + (m(4) * pkin(7) + t349) * pkin(7) - t368) * t382) * qJD(1) - t466 + t469 - t475) * t331 + t491 * t73 / 0.2e1 + (t352 * pkin(7) + t360 - 0.2e1 * t468 + t505) * t380 - (t367 + t366 + t470 + t471) * t331 / 0.2e1 + t290 * t70 + (t171 * t424 + pkin(7) * t179 + t170 * t425 + t348 * t437 + t346 * t436 + (-t166 * t326 - t167 * t330) * mrSges(4,3) + (t344 * t426 + t296 * t350 + t345 * t433 + t347 * t432 + t204 * t425 - t330 * t203 / 0.2e1 + (t227 * t326 - t228 * t330) * mrSges(4,3)) * qJD(3) + (-pkin(7) * t294 + t359 + ((-0.3e1 / 0.2e1 * Ifges(3,4) + t411 / 0.2e1 - t410 / 0.2e1) * t327 + Ifges(5,5) * t434 + Ifges(5,6) * t435 - 0.2e1 * t459 + t370 * t194 + t369 * t193) * qJD(1) + t465) * qJD(2)) * t327 + t246 * t160 + t247 * t216 + t248 * t217 + t186 * t243 + t187 * t244 + t245 * (-mrSges(5,1) * t169 + mrSges(5,2) * t168) + t229 * t16 + t98 * (-mrSges(6,1) * t193 + mrSges(6,2) * t194) + t23 * (-mrSges(7,1) * t193 + mrSges(7,2) * t194) + t66 * t190 + t67 * t191 + t175 * (-mrSges(6,1) * t74 + mrSges(6,2) * t73) + t168 * t146 / 0.2e1 + t162 * t111 + t163 * t112 + t143 * t85 + t149 * t15 + t8 * t136 + t9 * t137 + t7 * t138 + t10 * t139 + t95 * (-mrSges(7,1) * t74 + mrSges(7,2) * t73) + t61 * t84 + t71 * t48 + t72 * t50 + t62 * t47 + t63 * t49 + (Ifges(5,5) * t168 + Ifges(5,6) * t169) * t427 + t169 * t500 + t193 * t530 + t194 * t531 + t74 * t532; (t509 * t105 + t510 * t106 + t218 * t317 + t237 * t36 + t238 * t35 + t245 * t535) * m(5) + (t105 * t548 - t106 * t540 - t282 * t36 - t341 * t35) * mrSges(5,3) + (mrSges(5,1) * t540 - mrSges(5,2) * t548) * t245 + (t160 * t422 + t333) * qJD(3) - m(4) * (t227 * t241 + t228 * t242) + t218 * (mrSges(5,1) * t341 + mrSges(5,2) * t282) + (Ifges(5,4) * t282 - Ifges(5,2) * t341) * t454 + (Ifges(5,1) * t282 - Ifges(5,4) * t341) * t455 - t341 * t68 / 0.2e1 + (-Ifges(5,1) * t231 - Ifges(5,4) * t232) * t440 + (-Ifges(5,4) * t231 - Ifges(5,2) * t232) * t442 + (-Ifges(5,5) * t231 - Ifges(5,6) * t232) * t427 + t343 * mrSges(4,3) + (t249 / 0.2e1 - t232 / 0.2e1) * t145 + (t250 / 0.2e1 - t231 / 0.2e1) * t146 + (-Ifges(5,1) * t250 - Ifges(5,4) * t249) * t441 + (-Ifges(5,4) * t250 - Ifges(5,2) * t249) * t443 + (-Ifges(5,5) * t250 - Ifges(5,6) * t249) * t428 + (-t216 * t326 + t217 * t330 + m(4) * (-qJD(3) * t342 + t343) + (-t243 * t326 - t244 * t330) * qJD(3)) * pkin(8) + t345 * t436 + t347 * t437 + (t184 * t496 + t185 * t497) * t447 + (t223 * t496 + t224 * t497) * t464 + (t109 * t497 + t110 * t496) * t446 + (-t109 * t495 - t110 * t493) * t429 + (-t184 * t493 - t185 * t495) * t430 + (t184 * t494 + t185 * t496) * t499 + (t223 * t494 + t224 * t496) * t463 + (t109 * t496 + t110 * t494) * t450 + t514 * t139 + (t140 * t6 + t141 * t5 + t175 * t508 + t251 * t98 + t33 * t514 + t34 * t513) * m(6) + t515 * t136 + t516 * t138 + (t100 * t2 + t101 * t3 + t178 * t23 + t25 * t516 + t27 * t515 + t512 * t95) * m(7) + t491 * (-t185 / 0.2e1 + t109 / 0.2e1) + t492 * (-t184 / 0.2e1 + t110 / 0.2e1) + t513 * t137 + t512 * t84 + t510 * t190 + t508 * t85 + t509 * t191 + ((t468 + t360 - t318 / 0.2e1 + ((-m(4) * pkin(2) - mrSges(3,1) - t350) * qJD(2) - t352) * pkin(7) - t505) * t331 + ((t294 + t398) * pkin(7) + t359 + (t417 / 0.2e1 + t459 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t331) * qJD(1) - t465) * t327 + (Ifges(5,5) * t282 - Ifges(5,6) * t341 - t223 * t493 - t224 * t495 + t344) * t382 / 0.2e1) * qJD(1) + t326 * t171 / 0.2e1 + t317 * t70 + t282 * t69 / 0.2e1 - t274 * t160 + t251 * t16 + t237 * t111 + t238 * t112 - t242 * t243 - t241 * t244 + t23 * (-mrSges(7,1) * t223 + mrSges(7,2) * t224) + t98 * (-mrSges(6,1) * t223 + mrSges(6,2) * t224) + t178 * t15 - pkin(2) * t179 + t140 * t48 + t141 * t50 + t100 * t47 + t101 * t49 + (-mrSges(6,1) * t388 + mrSges(6,2) * t389) * t175 + (-mrSges(7,1) * t388 + mrSges(7,2) * t389) * t95 + (t223 * t5 - t224 * t6 - t33 * t389 + t34 * t388) * mrSges(6,3) + (-t2 * t224 + t223 * t3 - t25 * t389 + t27 * t388) * mrSges(7,3) + t170 * t424 + t223 * t530 + t224 * t531; t157 * t532 + (-Ifges(4,2) * t280 + t204 + t275) * t433 - t469 + t222 * t500 - m(5) * (t105 * t114 + t106 * t115) + t366 + (t329 * t111 + t325 * t112 - t280 * t160 + (t190 * t329 - t191 * t325) * qJD(4) + (-t105 * t377 + t106 * t376 + 0.2e1 * t245 * t432 + t325 * t35 + t329 * t36) * m(5)) * pkin(3) + t203 * t431 + (Ifges(4,1) * t279 - t404) * t432 + t479 * t136 + (t2 * t258 + t25 * t477 + t260 * t3 + t27 * t479 - t95 * t99) * m(7) + t476 * t139 + t477 * t138 + (-t175 * t180 + t259 * t6 + t260 * t5 + t33 * t476 + t34 * t478) * m(6) + t478 * t137 + (t227 * t279 + t228 * t280) * mrSges(4,3) - t296 * (mrSges(4,1) * t280 + mrSges(4,2) * t279) + t258 * t47 + t259 * t48 - t227 * t243 + t228 * t244 - t115 * t190 - t114 * t191 - t180 * t85 + t418 * t260 - t99 * t84 + (Ifges(4,5) * t279 - Ifges(4,6) * t280) * t426 + t544 + t546; (-t222 * t85 + t328 * t48 + t418 * t324 + ((t136 + t137) * t328 + (-t138 - t139) * t324) * qJD(5) + (-t175 * t222 + t324 * t5 + t328 * t6 - t33 * t375 + t34 * t374) * m(6)) * pkin(4) + (t2 * t315 - t102 * t95 - t25 * t28 - t27 * t29 + (-t25 * t375 + t27 * t374 + t3 * t324) * pkin(4)) * m(7) - m(6) * (t33 * t37 + t34 * t38) + t145 * t440 + t315 * t47 - t105 * t190 + t106 * t191 - t29 * t136 - t38 * t137 - t28 * t138 - t37 * t139 - t102 * t84 + t521 + t544; (-t157 * t84 + t47) * pkin(5) - t26 * t136 - t33 * t137 + t27 * t138 + t34 * t139 + (-(-t25 + t26) * t27 + (-t157 * t95 + t2) * pkin(5)) * m(7) + t521 + t547; -t504 * t136 + t157 * t138 + 0.2e1 * (t23 / 0.2e1 + t27 * t499 + t25 * t446) * m(7) + t15;];
tauc  = t1(:);
