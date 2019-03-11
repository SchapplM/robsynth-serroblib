% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:44
% EndTime: 2019-03-09 08:29:56
% DurationCPUTime: 5.82s
% Computational Cost: add. (10599->475), mult. (20369->623), div. (0->0), fcn. (21474->6), ass. (0->244)
t269 = sin(pkin(9));
t271 = sin(qJ(2));
t373 = cos(pkin(9));
t403 = cos(qJ(2));
t242 = t269 * t271 - t373 * t403;
t270 = sin(qJ(5));
t267 = t270 ^ 2;
t272 = cos(qJ(5));
t268 = t272 ^ 2;
t320 = t268 / 0.2e1 + t267 / 0.2e1;
t448 = (mrSges(6,3) + mrSges(7,3)) * t320;
t451 = t242 * t448;
t425 = -m(7) / 0.2e1;
t436 = Ifges(7,6) + Ifges(6,6);
t450 = t270 * t436;
t244 = t269 * t403 + t271 * t373;
t261 = -pkin(2) * t403 - pkin(1);
t286 = -t244 * qJ(4) + t261;
t152 = t242 * pkin(3) + t286;
t449 = m(5) * t152 - mrSges(5,2) * t242 - mrSges(5,3) * t244;
t406 = -t270 / 0.2e1;
t405 = -t272 / 0.2e1;
t266 = t271 * pkin(2);
t316 = qJ(4) * t242 + t266;
t416 = pkin(3) + pkin(8);
t101 = t244 * t416 + t316;
t249 = (-qJ(3) - pkin(7)) * t271;
t341 = t403 * pkin(7);
t251 = qJ(3) * t403 + t341;
t433 = t269 * t249 + t373 * t251;
t442 = -t242 * pkin(4) + t433;
t108 = t272 * t442;
t57 = -t101 * t270 + t108;
t58 = t272 * t101 + t270 * t442;
t301 = t58 * t270 + t57 * t272;
t446 = m(6) * t301;
t445 = -mrSges(4,2) + mrSges(5,3);
t437 = Ifges(7,5) + Ifges(6,5);
t42 = -pkin(5) * t242 + t108 + (-qJ(6) * t244 - t101) * t270;
t363 = t244 * t272;
t50 = qJ(6) * t363 + t58;
t443 = t270 * t50 + t272 * t42;
t330 = t373 * pkin(2);
t260 = -t330 - pkin(3);
t256 = -pkin(8) + t260;
t434 = -qJ(6) + t256;
t224 = t434 * t270;
t225 = t434 * t272;
t295 = t224 * t272 - t225 * t270;
t397 = t244 * pkin(5);
t178 = -t373 * t249 + t251 * t269;
t114 = pkin(4) * t244 + t178;
t109 = t272 * t114;
t100 = t242 * t416 + t286;
t315 = qJ(6) * t242 + t100;
t48 = -t270 * t315 + t109;
t41 = t48 + t397;
t369 = t114 * t270;
t49 = t272 * t315 + t369;
t303 = t270 * t49 + t272 * t41;
t366 = t242 * t270;
t157 = t244 * mrSges(7,1) - mrSges(7,3) * t366;
t365 = t242 * t272;
t161 = -t244 * mrSges(7,2) + mrSges(7,3) * t365;
t353 = t157 * t405 + t161 * t406;
t441 = (t242 * t295 - t303) * t425 - t353;
t440 = m(7) * (t41 * t270 - t272 * t49);
t439 = mrSges(7,1) + mrSges(6,1);
t396 = mrSges(6,2) + mrSges(7,2);
t395 = mrSges(4,3) + mrSges(5,1);
t438 = -Ifges(4,4) - Ifges(5,6);
t435 = Ifges(7,3) + Ifges(6,3);
t422 = m(7) * pkin(5);
t333 = mrSges(7,1) + t422;
t347 = t267 + t268;
t317 = t347 * t244;
t264 = t272 * mrSges(7,1);
t318 = -t270 * mrSges(7,2) + t264;
t337 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t219 = Ifges(7,4) * t365;
t94 = Ifges(7,1) * t366 + t244 * Ifges(7,5) + t219;
t220 = Ifges(6,4) * t365;
t96 = Ifges(6,1) * t366 + t244 * Ifges(6,5) + t220;
t430 = t337 * t244 + t94 / 0.2e1 + t96 / 0.2e1;
t429 = -m(5) / 0.2e1;
t428 = m(5) / 0.2e1;
t427 = m(6) / 0.2e1;
t426 = m(6) / 0.4e1;
t424 = m(7) / 0.2e1;
t423 = m(4) * pkin(2);
t421 = mrSges(6,2) / 0.2e1;
t420 = -mrSges(7,2) / 0.2e1;
t419 = mrSges(7,2) / 0.2e1;
t418 = -mrSges(7,3) / 0.2e1;
t417 = -t42 / 0.2e1;
t399 = pkin(5) * t272;
t331 = -pkin(4) - t399;
t71 = t242 * t331 + t433;
t415 = m(7) * t71;
t414 = t442 / 0.2e1;
t385 = t244 * mrSges(6,1);
t158 = -mrSges(6,3) * t366 + t385;
t413 = t158 / 0.2e1;
t364 = t244 * t270;
t159 = -t242 * mrSges(7,1) - mrSges(7,3) * t364;
t412 = -t159 / 0.2e1;
t411 = -t242 / 0.2e1;
t400 = pkin(2) * t269;
t257 = qJ(4) + t400;
t248 = pkin(5) * t270 + t257;
t410 = t248 / 0.2e1;
t409 = t257 / 0.2e1;
t404 = t272 / 0.2e1;
t402 = m(7) * t242;
t401 = m(7) * t248;
t393 = Ifges(6,4) * t270;
t392 = Ifges(6,4) * t272;
t391 = Ifges(7,4) * t270;
t390 = Ifges(7,4) * t272;
t148 = t318 * t242;
t380 = t272 * mrSges(6,1);
t383 = t270 * mrSges(6,2);
t308 = t380 - t383;
t149 = t308 * t242;
t384 = t244 * mrSges(6,2);
t162 = mrSges(6,3) * t365 - t384;
t351 = t161 + t162;
t352 = t157 + t158;
t55 = -t100 * t270 + t109;
t56 = t100 * t272 + t369;
t5 = (t242 * t395 + t148 + t149) * t242 + (t395 * t244 + t351 * t270 + t352 * t272) * t244 + m(7) * (-t242 * t71 + t244 * t303) + m(6) * (-t442 * t242 + (t270 * t56 + t272 * t55) * t244) + (m(5) + m(4)) * (t178 * t244 - t242 * t433);
t389 = qJD(1) * t5;
t160 = -t242 * mrSges(6,1) - mrSges(6,3) * t364;
t163 = t242 * mrSges(7,2) + mrSges(7,3) * t363;
t164 = t242 * mrSges(6,2) + mrSges(6,3) * t363;
t226 = t244 * mrSges(5,2);
t227 = t244 * mrSges(4,1);
t307 = t270 * mrSges(6,1) + t272 * mrSges(6,2);
t287 = t242 * t307;
t294 = t224 * t270 + t225 * t272;
t343 = t423 / 0.2e1;
t367 = t242 * t257;
t379 = t272 * mrSges(7,2);
t250 = t270 * mrSges(7,1) + t379;
t368 = t242 * t250;
t274 = -t367 * t428 + (t256 * t317 - t367) * t427 - t368 / 0.2e1 - t287 / 0.2e1 + (-t248 * t424 - t269 * t343) * t242 + (t260 * t428 + t294 * t424 - t373 * t343 - t448) * t244;
t153 = pkin(3) * t244 + t316;
t278 = t153 * t428 + (-t270 * t57 + t272 * t58) * t427 + (-t270 * t42 + t272 * t50) * t424 + t271 * t343;
t6 = -t226 + t227 - t274 + t278 + (t159 + t160) * t406 + (t163 + t164) * t404 + t445 * t242;
t388 = qJD(1) * t6;
t150 = t318 * t244;
t151 = t308 * t244;
t336 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t304 = Ifges(7,2) * t272 + t391;
t90 = t244 * Ifges(7,6) + t242 * t304;
t305 = Ifges(6,2) * t272 + t393;
t92 = t244 * Ifges(6,6) + t242 * t305;
t285 = t90 / 0.2e1 + t92 / 0.2e1 + t336 * t244;
t70 = t244 * t331 - t178;
t91 = -Ifges(7,6) * t242 + t244 * t304;
t93 = -Ifges(6,6) * t242 + t244 * t305;
t95 = -Ifges(7,5) * t242 + (Ifges(7,1) * t270 + t390) * t244;
t97 = -Ifges(6,5) * t242 + (Ifges(6,1) * t270 + t392) * t244;
t1 = t261 * t227 - t152 * t226 + t57 * t158 + t41 * t159 + t55 * t160 + t50 * t161 + t58 * t162 + t49 * t163 + t56 * t164 - t70 * t148 + t114 * t149 - t71 * t150 - t442 * t151 + t42 * t157 + m(6) * (-t114 * t442 + t55 * t57 + t56 * t58) + m(7) * (t41 * t42 + t49 * t50 + t70 * t71) + (-t261 * mrSges(4,2) + t152 * mrSges(5,3) + (t91 / 0.2e1 + t93 / 0.2e1) * t272 + (t95 / 0.2e1 + t97 / 0.2e1) * t270 + (-t270 * t337 - t272 * t336 - t438) * t242) * t242 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t271 + (m(4) * t261 + mrSges(4,1) * t242) * pkin(2)) * t271 + (mrSges(4,2) * t266 + t438 * t244 + t285 * t272 + t430 * t270 + (-Ifges(4,1) + Ifges(5,3) + Ifges(4,2) - Ifges(5,2) - t435) * t242) * t244 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) - Ifges(3,2)) * t271 + Ifges(3,4) * t403) * t403 + t449 * t153;
t387 = t1 * qJD(1);
t386 = t242 * mrSges(5,1);
t288 = Ifges(7,2) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(7,1) / 0.2e1;
t332 = m(7) * (-t41 + t48);
t338 = Ifges(7,4) / 0.2e1 + Ifges(6,4) / 0.2e1;
t4 = -t48 * t161 + t49 * t157 - t442 * t287 + t56 * t158 - t55 * t162 + pkin(5) * t148 * t366 - t71 * t368 - t49 * t332 + ((-t219 / 0.2e1 - t220 / 0.2e1 + t41 * mrSges(7,3) + t55 * mrSges(6,3) - t430) * t272 + (t56 * mrSges(6,3) + t49 * mrSges(7,3) - pkin(5) * t415 + t288 * t365 + t338 * t366 + t285) * t270) * t242;
t377 = t4 * qJD(1);
t311 = t332 / 0.2e1;
t290 = -t157 / 0.2e1 + t311;
t342 = t422 / 0.2e1;
t313 = mrSges(6,1) / 0.2e1 + t342;
t292 = mrSges(7,1) / 0.2e1 + t313;
t321 = -t161 / 0.2e1 - t162 / 0.2e1;
t340 = t421 + t419;
t8 = (t244 * t340 + t321) * t272 + (t244 * t292 - t290 + t413) * t270 + t451;
t374 = t8 * qJD(1);
t14 = (-t351 * t272 + t352 * t270 + t440 + m(6) * (t270 * t55 - t272 * t56) - t449) * t244;
t372 = qJD(1) * t14;
t356 = t272 * t161;
t361 = t270 * t157;
t22 = (t356 - t361 - t440) * t242;
t371 = qJD(1) * t22;
t335 = -t48 / 0.2e1 + t41 / 0.2e1;
t324 = t363 / 0.2e1;
t348 = mrSges(7,1) * t324 + t364 * t420;
t10 = (-t384 / 0.2e1 - t321) * t270 + (t385 / 0.2e1 + t157 / 0.2e1 + t413 + (t397 / 0.2e1 + t335) * m(7)) * t272 + t348;
t370 = t10 * qJD(1);
t359 = t270 * t164;
t357 = t272 * t160;
t312 = (m(7) / 0.4e1 + t426) * t317;
t319 = m(7) * t347;
t33 = -0.2e1 * t312 + 0.2e1 * (t429 - m(6) * t347 / 0.4e1 - t319 / 0.4e1) * t244;
t355 = t33 * qJD(1);
t72 = 0.2e1 * (-0.1e1 / 0.4e1 - t267 / 0.4e1 - t268 / 0.4e1) * t402;
t354 = t72 * qJD(1);
t346 = qJD(5) * t270;
t345 = qJD(5) * t272;
t240 = (-0.1e1 / 0.2e1 - t320) * m(7);
t344 = t240 * qJD(2);
t339 = -0.3e1 / 0.4e1 * Ifges(6,4) - 0.3e1 / 0.4e1 * Ifges(7,4);
t325 = t364 / 0.2e1;
t314 = -mrSges(7,3) * pkin(5) + t437;
t309 = mrSges(7,2) * t325 + t424 * t70 - mrSges(7,1) * t363 / 0.2e1;
t252 = -t270 * Ifges(7,2) + t390;
t253 = -t270 * Ifges(6,2) + t392;
t254 = t272 * Ifges(7,1) - t391;
t255 = t272 * Ifges(6,1) - t393;
t24 = -t250 * t399 - t248 * t318 - t257 * t308 + (t254 / 0.2e1 + t255 / 0.2e1 - t338 * t270) * t270 + (-pkin(5) * t401 + t252 / 0.2e1 + t253 / 0.2e1 + t338 * t272 - t288 * t270) * t272;
t289 = Ifges(6,2) / 0.4e1 + Ifges(7,2) / 0.4e1 - Ifges(6,1) / 0.4e1 - Ifges(7,1) / 0.4e1;
t277 = t224 * t418 + mrSges(6,1) * t409 + mrSges(7,1) * t410 - t253 / 0.4e1 - t252 / 0.4e1 + t289 * t270 + (t250 / 0.2e1 + t401 / 0.2e1) * pkin(5);
t279 = t225 * t418 + mrSges(6,2) * t409 + mrSges(7,2) * t410 + t255 / 0.4e1 + t254 / 0.4e1 - t289 * t272;
t280 = t335 * mrSges(7,3) - t219 / 0.4e1 - t220 / 0.4e1 - t94 / 0.4e1 - t96 / 0.4e1 - t442 * mrSges(6,2) / 0.2e1 - t256 * t158 / 0.2e1 + t71 * t420;
t281 = (-t148 / 0.2e1 + t415 / 0.2e1) * pkin(5) - t90 / 0.4e1 - t92 / 0.4e1 + mrSges(6,1) * t414 + t256 * t162 / 0.2e1;
t282 = t290 * t224 + t225 * t161 / 0.2e1 + t71 * t264 / 0.2e1;
t284 = mrSges(7,1) * t417 + t50 * t419 - t57 * mrSges(6,1) / 0.2e1 + t58 * t421;
t291 = t320 * t256 * mrSges(6,3);
t3 = (m(7) * t417 + t412) * pkin(5) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 - t291) * t242 + ((-0.3e1 / 0.4e1 * Ifges(6,5) - 0.3e1 / 0.4e1 * Ifges(7,5)) * t244 + t277 * t242 + t280) * t270 + ((-0.3e1 / 0.4e1 * Ifges(7,6) - 0.3e1 / 0.4e1 * Ifges(6,6)) * t244 + (t270 * t339 + t279) * t242 + t281) * t272 + t282 + t284;
t300 = t3 * qJD(1) - t24 * qJD(2);
t275 = t433 * t428 + m(6) * t414 + (-t244 * t295 + t71) * t424 + t396 * t366 / 0.2e1 - t439 * t365 / 0.2e1;
t276 = t433 * t429 - t446 / 0.2e1 + t443 * t425 + t163 * t406 - t359 / 0.2e1 + t159 * t405 - t357 / 0.2e1;
t13 = t275 + t276;
t78 = mrSges(5,3) + t396 * t272 + t439 * t270 + t401 + 0.4e1 * (t426 + m(5) / 0.4e1) * t257;
t299 = qJD(1) * t13 + qJD(2) * t78;
t19 = t309 + t441;
t87 = -m(7) * t294 + mrSges(7,3) * t347;
t298 = -qJD(1) * t19 + qJD(2) * t87;
t262 = m(7) * t399;
t230 = -t262 - t318;
t79 = (t270 * t333 + t379) * t242;
t293 = qJD(1) * t79 - qJD(2) * t230;
t239 = t347 * t425 + t424;
t73 = t242 * t319 / 0.2e1 - t402 / 0.2e1;
t32 = -0.2e1 * t312 + (m(6) + m(7)) * t317 / 0.2e1;
t20 = t309 - t441;
t12 = t275 - t276 - t386;
t11 = t272 * t311 + t158 * t405 + t162 * t406 + (-t383 / 0.2e1 + t313 * t272) * t244 + t348 + t353;
t9 = t356 / 0.2e1 - t361 / 0.2e1 + t270 * t311 + t158 * t406 + t162 * t404 + (t270 * t292 + t272 * t340) * t244 - t451;
t7 = (t164 / 0.2e1 + t163 / 0.2e1) * t272 + (t412 - t160 / 0.2e1) * t270 + t274 + t278;
t2 = pkin(5) * t159 / 0.2e1 + t42 * t342 + ((-Ifges(7,6) / 0.4e1 - Ifges(6,6) / 0.4e1) * t244 + t281) * t272 + ((-Ifges(7,5) / 0.4e1 - Ifges(6,5) / 0.4e1) * t244 + t280) * t270 + (-t291 + t279 * t272 + (t272 * t339 + t277) * t270) * t242 + t282 - t284 + t435 * t411 + t437 * t325 + t436 * t324;
t15 = [qJD(2) * t1 + qJD(3) * t5 + qJD(4) * t14 - qJD(5) * t4 + qJD(6) * t22, t7 * qJD(3) + t12 * qJD(4) + t2 * qJD(5) + t20 * qJD(6) + t387 + (-t257 * t151 - t248 * t150 + t70 * t250 + t224 * t163 + t225 * t159 + (-t257 * mrSges(5,1) - mrSges(4,3) * t400 + Ifges(5,5) - Ifges(4,6)) * t244 + Ifges(3,5) * t403 - t260 * t386 + (mrSges(4,3) * t330 + Ifges(5,4) - Ifges(4,5)) * t242 - mrSges(3,1) * t341 + m(7) * (t224 * t50 + t225 * t42 + t248 * t70) + (t272 * t437 - t450) * t411 + (t93 + t91) * t406 + (t97 + t95) * t404 + (t255 + t254) * t325 + (t253 + t252) * t324 + (mrSges(3,2) * pkin(7) - Ifges(3,6)) * t271 + (-m(6) * t257 - t307) * t114 + (t357 + t359 + t446) * t256 + (-m(5) * t257 - t269 * t423 - t445) * t178 + (m(5) * t260 - t373 * t423 - mrSges(4,1) + mrSges(5,2)) * t433 - t443 * mrSges(7,3) - t301 * mrSges(6,3)) * qJD(2), qJD(2) * t7 + qJD(4) * t32 + qJD(5) * t11 + qJD(6) * t73 + t389, qJD(2) * t12 + qJD(3) * t32 + qJD(5) * t9 + t372, t2 * qJD(2) + t11 * qJD(3) + t9 * qJD(4) - t377 + (-mrSges(6,1) * t56 - mrSges(6,2) * t55 - mrSges(7,2) * t48 + (t272 * t314 - t450) * t242 - t333 * t49) * qJD(5), qJD(2) * t20 + qJD(3) * t73 + t371; -qJD(3) * t6 + qJD(4) * t13 + qJD(5) * t3 - qJD(6) * t19 - t387, qJD(4) * t78 - qJD(5) * t24 + qJD(6) * t87, -t388, qJD(6) * t239 + t299 (-mrSges(7,2) * t225 - t224 * t333) * qJD(5) + (-mrSges(6,2) * t256 - t436) * t345 + (-mrSges(6,1) * t256 - t314) * t346 + t300, qJD(4) * t239 + t298; qJD(2) * t6 + qJD(4) * t33 - qJD(5) * t10 - qJD(6) * t72 - t389, t388, 0, t355, -t370 + (t270 * t396 - t262 - t264 - t380) * qJD(5), -t354; -qJD(2) * t13 - qJD(3) * t33 - qJD(5) * t8 - t372, qJD(6) * t240 - t299, -t355, 0, -t374 - t396 * t345 + (-mrSges(6,1) - t333) * t346, t344; -qJD(2) * t3 + qJD(3) * t10 + qJD(4) * t8 - qJD(6) * t79 + t377, t230 * qJD(6) - t300, t370, t374, 0, -t293; qJD(2) * t19 + qJD(3) * t72 + qJD(5) * t79 - t371, -qJD(4) * t240 - qJD(5) * t230 - t298, t354, -t344, t293, 0;];
Cq  = t15;
