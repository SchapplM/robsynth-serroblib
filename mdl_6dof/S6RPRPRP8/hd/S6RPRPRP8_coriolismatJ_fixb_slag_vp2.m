% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:30
% EndTime: 2019-03-09 03:24:41
% DurationCPUTime: 6.16s
% Computational Cost: add. (10429->486), mult. (19188->642), div. (0->0), fcn. (19879->6), ass. (0->246)
t248 = sin(pkin(9));
t250 = sin(qJ(3));
t374 = cos(pkin(9));
t412 = cos(qJ(3));
t218 = t248 * t412 + t250 * t374;
t217 = t248 * t250 - t374 * t412;
t249 = sin(qJ(5));
t360 = t217 * t249;
t139 = -mrSges(6,2) * t218 + mrSges(6,3) * t360;
t395 = t218 * mrSges(7,3);
t143 = mrSges(7,2) * t360 + t395;
t337 = t139 + t143;
t246 = t249 ^ 2;
t251 = cos(qJ(5));
t247 = t251 ^ 2;
t331 = t246 + t247;
t418 = -t249 / 0.2e1;
t415 = -t251 / 0.2e1;
t414 = t251 / 0.2e1;
t351 = t247 * t217;
t353 = t246 * t217;
t461 = t351 / 0.2e1 + t353 / 0.2e1;
t460 = t218 * mrSges(5,1) - t217 * mrSges(5,2);
t235 = t250 * pkin(3) + qJ(2);
t459 = m(5) * t235;
t458 = -mrSges(6,3) - mrSges(7,2);
t447 = Ifges(7,4) + Ifges(6,5);
t457 = Ifges(6,6) - Ifges(7,6);
t234 = -pkin(3) * t374 - pkin(4);
t288 = t251 * pkin(5) + t249 * qJ(6);
t213 = -t288 + t234;
t384 = t251 * mrSges(7,1);
t390 = t249 * mrSges(7,3);
t224 = -t384 - t390;
t456 = m(7) * t213 + t224;
t401 = Ifges(7,5) * t251;
t229 = t249 * Ifges(7,1) - t401;
t243 = Ifges(6,4) * t251;
t230 = t249 * Ifges(6,1) + t243;
t455 = t230 + t229;
t166 = t217 ^ 2;
t327 = t166 / 0.2e1;
t436 = t218 ^ 2;
t454 = t327 + t436 / 0.2e1;
t233 = pkin(3) * t248 + pkin(8);
t453 = t331 * t233;
t452 = t217 * t248 + t374 * t218;
t324 = t412 * pkin(3);
t138 = -t217 * pkin(4) + t218 * pkin(8) + t324;
t252 = -pkin(1) - pkin(7);
t222 = (-qJ(4) + t252) * t250;
t318 = t412 * t252;
t223 = -qJ(4) * t412 + t318;
t151 = t222 * t248 - t374 * t223;
t64 = t138 * t251 + t151 * t249;
t65 = t249 * t138 - t151 * t251;
t280 = -t249 * t64 + t251 * t65;
t51 = -qJ(6) * t217 + t65;
t52 = pkin(5) * t217 - t64;
t283 = t249 * t52 + t251 * t51;
t359 = t217 * t251;
t141 = mrSges(6,1) * t218 + mrSges(6,3) * t359;
t142 = -mrSges(7,1) * t218 - mrSges(7,2) * t359;
t451 = t141 * t415 + t142 * t414 + t337 * t418;
t448 = -mrSges(6,1) - mrSges(7,1);
t137 = pkin(4) * t218 + pkin(8) * t217 + t235;
t446 = t374 * t222 + t248 * t223;
t62 = t137 * t251 - t249 * t446;
t48 = -pkin(5) * t218 - t62;
t405 = t48 + t62;
t358 = t218 * qJ(6);
t339 = t251 * t446;
t348 = t249 * t137;
t63 = t339 + t348;
t47 = t63 + t358;
t336 = t141 - t142;
t385 = t251 * mrSges(6,1);
t391 = t249 * mrSges(6,2);
t225 = -t385 + t391;
t445 = t224 + t225;
t242 = Ifges(7,5) * t249;
t227 = -t251 * Ifges(7,3) + t242;
t301 = Ifges(7,1) * t251 + t242;
t402 = Ifges(6,4) * t249;
t292 = Ifges(6,1) * t251 - t402;
t444 = t292 + t301;
t443 = t461 * t233;
t290 = Ifges(6,2) * t249 - t243;
t439 = -t290 + t455;
t241 = m(7) * qJ(6) + mrSges(7,3);
t101 = t218 * Ifges(7,4) - t217 * t301;
t103 = t218 * Ifges(6,5) - t217 * t292;
t322 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t397 = t151 * mrSges(6,2);
t373 = qJ(6) * t251;
t409 = t249 * pkin(5);
t287 = -t373 + t409;
t67 = -t217 * t287 + t151;
t407 = t67 * mrSges(7,3);
t438 = -t322 * t218 + t62 * mrSges(6,3) + t407 - t101 / 0.2e1 - t103 / 0.2e1 - t397;
t406 = t47 - t63;
t430 = m(7) / 0.2e1;
t437 = (-t249 * t406 + t251 * t405) * t430 + t451;
t434 = m(5) / 0.2e1;
t433 = -m(6) / 0.2e1;
t432 = m(6) / 0.2e1;
t431 = -m(7) / 0.2e1;
t429 = m(5) * pkin(3);
t428 = mrSges(6,1) / 0.2e1;
t427 = -mrSges(6,2) / 0.2e1;
t135 = t288 * t217;
t425 = t135 / 0.2e1;
t424 = t213 / 0.2e1;
t423 = -t217 / 0.2e1;
t422 = t217 / 0.2e1;
t421 = -t218 / 0.2e1;
t419 = -t287 / 0.2e1;
t417 = t249 / 0.2e1;
t411 = m(7) * t287;
t410 = m(7) * t251;
t408 = t67 * mrSges(7,1);
t404 = m(7) * qJD(6);
t382 = t251 * mrSges(7,3);
t392 = t249 * mrSges(7,1);
t293 = -t382 + t392;
t136 = t217 * t293;
t281 = t249 * t62 - t251 * t63;
t285 = -t249 * t48 - t251 * t47;
t383 = t251 * mrSges(6,2);
t393 = t249 * mrSges(6,1);
t294 = t383 + t393;
t366 = t151 * t217;
t9 = (mrSges(5,3) * t218 + t249 * t336 - t251 * t337) * t218 + (t136 + (mrSges(5,3) + t294) * t217) * t217 + m(5) * (-t218 * t446 - t366) + m(6) * (t218 * t281 - t366) + m(7) * (-t217 * t67 + t218 * t285);
t399 = qJD(1) * t9;
t398 = t151 * mrSges(6,1);
t396 = t217 * mrSges(7,1);
t394 = t218 * Ifges(6,6);
t386 = t250 * mrSges(4,1);
t204 = Ifges(7,6) * t359;
t205 = Ifges(7,5) * t359;
t97 = Ifges(7,6) * t218 - Ifges(7,3) * t360 - t205;
t99 = t217 * t290 + t394;
t267 = t63 * mrSges(6,3) - t97 / 0.2e1 + t99 / 0.2e1 - t398 - t408;
t4 = t204 * t421 + ((t205 / 0.2e1 + t394 / 0.2e1 + t47 * mrSges(7,2) - Ifges(6,4) * t359 / 0.2e1 + t267) * t251 + (t48 * mrSges(7,2) + (-Ifges(7,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t360 + (Ifges(7,3) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(7,1) / 0.2e1) * t359 - t438) * t249) * t217 - (m(7) * t67 - t136) * t135 + (m(7) * t48 - t336) * t63 + (m(7) * t47 + t337) * t62;
t378 = t4 * qJD(1);
t377 = t47 * t249;
t376 = t63 * t249;
t263 = t385 / 0.2e1 - t391 / 0.2e1 + t384 / 0.2e1 + t390 / 0.2e1;
t260 = t288 * t430 + t263;
t264 = t458 * (t246 / 0.2e1 + t247 / 0.2e1);
t282 = t62 * t251 + t376;
t284 = -t48 * t251 + t377;
t7 = (m(7) * t425 + t217 * t263) * t217 + ((t282 - t284) * t431 + t264 * t217 - t451) * t218 + t260;
t375 = t7 * qJD(1);
t269 = mrSges(4,2) * t412 + t386 + t460;
t15 = mrSges(3,3) + t336 * t251 + t337 * t249 + (m(4) + m(3)) * qJ(2) + t459 + m(6) * t282 + m(7) * t284 + t269;
t372 = qJD(1) * t15;
t22 = m(7) * (t218 * t47 + t359 * t67) - t136 * t359 + t218 * t143;
t371 = qJD(1) * t22;
t70 = (0.1e1 / 0.2e1 + t454) * t410;
t370 = qJD(1) * t70;
t210 = t218 * mrSges(5,2);
t333 = t453 * t218;
t362 = t217 * t234;
t363 = t217 * t224;
t365 = t213 * t217;
t254 = (-t333 - t362) * t432 + (-t333 - t365) * t430 - t363 / 0.2e1 + t225 * t423 + (t217 * t374 - t218 * t248) * t429 / 0.2e1 - t458 * t331 * t421;
t355 = t218 * t251;
t323 = mrSges(7,2) * t355;
t140 = -t323 + t396;
t356 = t218 * t249;
t144 = mrSges(7,2) * t356 - t217 * mrSges(7,3);
t273 = mrSges(6,2) * t217 + mrSges(6,3) * t356;
t274 = -mrSges(6,1) * t217 + mrSges(6,3) * t355;
t255 = (t249 * t65 + t251 * t64) * t432 + (t249 * t51 - t251 * t52) * t430 + t140 * t415 + t274 * t414 + t324 * t434 + (t144 + t273) * t417;
t10 = t217 * mrSges(5,1) + t210 + t254 - t255;
t369 = t10 * qJD(1);
t257 = (t249 * t405 + t251 * t406) * t431 + t141 * t417 + t142 * t418 + t337 * t415;
t259 = t383 / 0.2e1 + t393 / 0.2e1 + t392 / 0.2e1 - t382 / 0.2e1 + t287 * t430;
t258 = t259 * t218;
t12 = t258 + t257;
t368 = t12 * qJD(1);
t364 = t213 * t218;
t357 = t218 * t234;
t266 = -m(5) / 0.2e1 + (t431 + t433) * t331;
t326 = m(6) / 0.4e1 + m(7) / 0.4e1;
t268 = (-t166 - t436) * t434 + 0.2e1 * t326 * (-t331 * t436 - t166);
t24 = t266 + t268;
t354 = t24 * qJD(1);
t346 = t249 * t140;
t344 = t249 * t227;
t228 = t251 * Ifges(6,2) + t402;
t343 = t249 * t228;
t340 = t251 * t144;
t338 = t251 * t224;
t332 = t453 * t217;
t330 = qJD(5) * t249;
t329 = qJD(5) * t251;
t132 = m(7) * t356;
t328 = t132 * qJD(1);
t325 = t52 * t430;
t321 = -Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t320 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t306 = t139 / 0.2e1 + t143 / 0.2e1;
t305 = -t141 / 0.2e1 + t142 / 0.2e1;
t303 = qJ(6) * mrSges(7,2) + Ifges(6,6);
t299 = -mrSges(7,2) * pkin(5) + t447;
t289 = Ifges(7,3) * t249 + t401;
t100 = -t217 * Ifges(7,4) - t218 * t301;
t102 = -t217 * Ifges(6,5) - t218 * t292;
t68 = -t218 * t287 + t446;
t96 = -t217 * Ifges(7,6) - t218 * t289;
t98 = -t217 * Ifges(6,6) + t218 * t290;
t1 = -t68 * t136 + t65 * t139 + t48 * t140 + t64 * t141 + t52 * t142 + t51 * t143 + t47 * t144 - t235 * t210 + m(6) * (t151 * t446 + t62 * t64 + t63 * t65) + m(7) * (t47 * t51 + t48 * t52 + t67 * t68) + (-t235 * mrSges(5,1) - t62 * mrSges(6,1) + t63 * mrSges(6,2) - Ifges(5,4) * t217 + (-t100 / 0.2e1 - t102 / 0.2e1 - t446 * mrSges(6,2) + t322 * t217) * t251 + (t98 / 0.2e1 - t96 / 0.2e1 - t446 * mrSges(6,1) + t320 * t217) * t249) * t217 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t250) * t250 + (Ifges(5,4) * t218 + t438 * t251 + (-t320 * t218 + t267) * t249 + (-Ifges(7,2) - Ifges(6,3) + Ifges(5,1) - Ifges(5,2)) * t217) * t218 + (qJ(2) * mrSges(4,1) + (-Ifges(4,1) + Ifges(4,2)) * t250 - Ifges(4,4) * t412) * t412 + (t459 + t460) * t324;
t6 = (t306 * t251 + t305 * t249 + (t446 + t281) * t433 + (t285 + t68) * t431) * t217 + (t294 * t422 - t346 / 0.2e1 + t136 - t340 / 0.2e1 + (t151 + t280) * t433 + (t283 + t67) * t431) * t218;
t286 = t1 * qJD(1) - t6 * qJD(2);
t25 = 0.4e1 * t326 * (0.1e1 - t331) * t218 * t217;
t279 = -t6 * qJD(1) + t25 * qJD(2);
t278 = -t135 * t213 + t287 * t67;
t146 = t456 * t249;
t262 = (-t249 * t67 + (t218 * t233 + t365) * t251) * t430 - t136 * t418;
t20 = -t323 + (mrSges(7,1) / 0.2e1 - t338 / 0.2e1) * t217 + t325 - t262;
t277 = qJD(1) * t20 + qJD(3) * t146;
t31 = t395 + 0.2e1 * (t358 / 0.2e1 + t348 / 0.4e1 + t339 / 0.4e1 - t63 / 0.4e1) * m(7);
t276 = qJD(1) * t31 + qJD(5) * t241;
t275 = Ifges(7,3) / 0.4e1 + Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1 - Ifges(7,1) / 0.4e1;
t272 = -t136 * t419 + t224 * t425;
t256 = (-pkin(5) * t52 + qJ(6) * t51) * t430 - pkin(5) * t140 / 0.2e1 + qJ(6) * t144 / 0.2e1 + t51 * mrSges(7,3) / 0.2e1 - t52 * mrSges(7,1) / 0.2e1 + t64 * t428 + t65 * t427;
t2 = t278 * t431 + (-t397 / 0.2e1 - t103 / 0.4e1 - t101 / 0.4e1 + t407 / 0.2e1 + (-0.3e1 / 0.4e1 * Ifges(7,4) - 0.3e1 / 0.4e1 * Ifges(6,5)) * t218 + (-t62 / 0.2e1 - t48 / 0.2e1) * mrSges(7,2) + (t405 * t431 - t305) * t233) * t251 + (-t398 / 0.2e1 - t408 / 0.2e1 + t99 / 0.4e1 - t97 / 0.4e1 + t205 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(6,6) - 0.3e1 / 0.4e1 * Ifges(7,6)) * t218 + (-t63 / 0.2e1 + t47 / 0.2e1) * mrSges(7,2) + (t406 * t430 + t306) * t233) * t249 + (t264 * t233 + (t242 / 0.4e1 - t228 / 0.4e1 + t227 / 0.4e1 + t234 * t428 + mrSges(7,1) * t424 - t275 * t251) * t251 + (-t230 / 0.4e1 - t229 / 0.4e1 - t243 / 0.4e1 + t234 * t427 + mrSges(7,3) * t424 + t275 * t249 + (-0.3e1 / 0.4e1 * Ifges(6,4) + Ifges(7,5) / 0.2e1) * t251) * t249 + t321) * t217 + t256 + t272;
t27 = t234 * t294 + t287 * t224 + t344 / 0.2e1 + t289 * t415 - t343 / 0.2e1 + (t293 + t411) * t213 + t444 * t417 + t439 * t414;
t28 = (t409 / 0.2e1 - t373 / 0.2e1 + t419) * m(7) * t217;
t271 = t2 * qJD(1) + t28 * qJD(2) - t27 * qJD(3);
t261 = (-m(7) * t288 + t445) * qJD(5);
t209 = (m(7) * t233 + mrSges(7,2)) * t251;
t129 = m(7) * t360;
t71 = (-0.1e1 / 0.2e1 + t454) * t410;
t29 = t217 * t259 - t411 * t423 + (t293 + t294) * t422;
t26 = 0.2e1 * t430 * t47 + t143;
t23 = -t266 + t268;
t21 = t338 * t422 + t396 / 0.2e1 + t325 + t262;
t13 = t258 - t257;
t11 = t254 + t255;
t8 = -t135 * t217 * t430 + t327 * t445 + t260 + (-t458 * t461 + t437) * t218;
t5 = t6 * qJD(3);
t3 = -t272 + t278 * t430 + t256 - t249 * t99 / 0.4e1 + t225 * t362 / 0.2e1 - t289 * t360 / 0.4e1 + t363 * t424 + (-t249 * t320 - t251 * t322) * t218 + t321 * t217 + t67 * t293 / 0.2e1 + t151 * t294 / 0.2e1 + (-t457 * t249 + t447 * t251) * t218 / 0.4e1 + (Ifges(7,1) * t360 - t205 + t97) * t249 / 0.4e1 + (t103 + t101) * t251 / 0.4e1 + t443 * mrSges(6,3) + (t439 + t230) * t360 / 0.4e1 + t437 * t233 + (t228 / 0.2e1 - t227 / 0.2e1 - t444 / 0.4e1) * t359 + (t405 * t414 - t377 / 0.2e1 + t376 / 0.2e1 + t443) * mrSges(7,2);
t14 = [qJD(2) * t15 + qJD(3) * t1 + qJD(4) * t9 + qJD(5) * t4 + qJD(6) * t22, qJD(4) * t23 + qJD(5) * t8 + qJD(6) * t71 + t372 - t5, t11 * qJD(4) + t3 * qJD(5) + t21 * qJD(6) + t286 + (-Ifges(4,5) * t250 - Ifges(4,6) * t412 - t252 * t386 + Ifges(5,6) * t217 - t293 * t364 - t294 * t357 + t344 * t421 - mrSges(4,2) * t318 + t98 * t414 + t96 * t415 + (m(6) * t234 - t374 * t429 - mrSges(5,1) + t225) * t446 + (m(6) * t280 + m(7) * t283 - t249 * t274 + t251 * t273 + t340 + t346) * t233 + t456 * t68 + t452 * mrSges(5,3) * pkin(3) + (t249 * t447 + t251 * t457) * t423 + (t102 + t100) * t417 - t455 * t355 / 0.2e1 + (-Ifges(5,5) + t343 / 0.2e1) * t218 + (-t248 * t429 + mrSges(5,2)) * t151 + t280 * mrSges(6,3) + t283 * mrSges(7,2)) * qJD(3), qJD(2) * t23 + qJD(3) * t11 + qJD(5) * t13 + t399, t8 * qJD(2) + t3 * qJD(3) + t13 * qJD(4) + t26 * qJD(6) + t378 + (-t204 + (-m(7) * pkin(5) + t448) * t63 + (-mrSges(6,2) + t241) * t62 + (t249 * t299 + t251 * t303) * t217) * qJD(5), qJD(2) * t71 + qJD(3) * t21 + qJD(5) * t26 + t371; qJD(4) * t24 - qJD(5) * t7 + qJD(6) * t70 - t372 - t5, t25 * qJD(3) (m(6) * (-t332 + t357) + m(7) * (-t332 + t364) - t452 * t429 - t269 + t445 * t218 - t458 * (-t353 - t351)) * qJD(3) + t29 * qJD(5) - t129 * qJD(6) + t279, t354, -t375 + t29 * qJD(3) + (t251 * t404 + t261) * t218, m(7) * t218 * t329 - qJD(3) * t129 + t370; qJD(4) * t10 - qJD(5) * t2 - qJD(6) * t20 - t286, -qJD(5) * t28 - t279, qJD(5) * t27 - qJD(6) * t146, t369, t209 * qJD(6) + t299 * t329 + (Ifges(7,6) - t303) * t330 + t233 * t261 - t271, qJD(5) * t209 - t277; -qJD(2) * t24 - qJD(3) * t10 - qJD(5) * t12 + qJD(6) * t132 - t399, -t354, -t369, 0, -t368 + (t382 - t383 - t411) * qJD(5) + (qJD(5) * t448 + t404) * t249, m(7) * t330 + t328; qJD(2) * t7 + qJD(3) * t2 + qJD(4) * t12 + qJD(6) * t31 - t378, qJD(3) * t28 + t375, t271, t368, t241 * qJD(6), t276; -qJD(2) * t70 + qJD(3) * t20 - qJD(4) * t132 - qJD(5) * t31 - t371, -t370, t277, -t328, -t276, 0;];
Cq  = t14;
