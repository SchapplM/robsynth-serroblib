% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR15_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:15
% EndTime: 2019-12-31 20:41:24
% DurationCPUTime: 4.50s
% Computational Cost: add. (9537->417), mult. (18097->557), div. (0->0), fcn. (17099->6), ass. (0->218)
t290 = sin(qJ(5));
t291 = sin(qJ(4));
t293 = cos(qJ(5));
t294 = cos(qJ(4));
t248 = t290 * t291 - t293 * t294;
t292 = sin(qJ(2));
t214 = t248 * t292;
t318 = t290 * t294 + t293 * t291;
t217 = t318 * t292;
t447 = m(6) * pkin(4);
t458 = t217 * mrSges(6,1) / 0.2e1 - t214 * mrSges(6,2) / 0.2e1;
t472 = -(t214 * t290 + t217 * t293) * t447 / 0.2e1 - t458;
t295 = cos(qJ(2));
t316 = -Ifges(6,5) * t217 / 0.2e1 + Ifges(6,6) * t214 / 0.2e1;
t380 = qJ(3) * t295;
t444 = pkin(2) + pkin(7);
t246 = t444 * t292 - t380;
t443 = pkin(3) + pkin(6);
t268 = t443 * t295;
t252 = t294 * t268;
t123 = pkin(4) * t295 + t252 + (-pkin(8) * t292 - t246) * t291;
t158 = t294 * t246 + t291 * t268;
t364 = t292 * t294;
t137 = pkin(8) * t364 + t158;
t62 = t123 * t293 - t137 * t290;
t63 = t123 * t290 + t137 * t293;
t454 = t316 + t63 * mrSges(6,2) / 0.2e1 - t62 * mrSges(6,1) / 0.2e1;
t452 = Ifges(6,3) * t295 / 0.2e1 - t454;
t336 = -qJ(3) * t292 - pkin(1);
t230 = -t444 * t295 + t336;
t267 = t443 * t292;
t334 = t230 * t291 - t294 * t267;
t366 = t291 * t295;
t134 = pkin(8) * t366 - t334;
t122 = pkin(4) * t292 + t134;
t152 = t230 * t294 + t267 * t291;
t358 = t294 * t295;
t135 = -pkin(8) * t358 + t152;
t379 = t135 * t290;
t58 = t122 * t293 - t379;
t65 = t134 * t293 - t379;
t471 = -t58 + t65;
t365 = t291 * t444;
t260 = -t291 * pkin(8) - t365;
t261 = (-pkin(8) - t444) * t294;
t167 = t260 * t293 + t261 * t290;
t333 = -t260 * t290 + t293 * t261;
t353 = -Ifges(6,5) * t318 + Ifges(6,6) * t248;
t25 = -t167 * mrSges(6,1) - t333 * mrSges(6,2) + t353;
t470 = t25 * qJD(5);
t335 = -mrSges(6,1) * t318 + t248 * mrSges(6,2);
t469 = qJD(5) * t335;
t279 = pkin(4) * t291 + qJ(3);
t422 = m(6) * t279;
t378 = t135 * t293;
t59 = t122 * t290 + t378;
t64 = -t134 * t290 - t378;
t468 = t59 + t64;
t465 = t248 * t290 + t318 * t293;
t215 = t248 * t295;
t216 = t318 * t295;
t393 = t216 * Ifges(6,4);
t108 = t215 * Ifges(6,2) + t292 * Ifges(6,6) - t393;
t204 = Ifges(6,4) * t215;
t109 = -t216 * Ifges(6,1) + t292 * Ifges(6,5) + t204;
t124 = -mrSges(6,1) * t216 + mrSges(6,2) * t215;
t126 = Ifges(6,2) * t216 + t204;
t127 = Ifges(6,1) * t215 + t393;
t159 = -mrSges(6,1) * t248 - mrSges(6,2) * t318;
t229 = pkin(4) * t358 + t268;
t394 = t216 * mrSges(6,3);
t180 = mrSges(6,1) * t292 + t394;
t435 = -t180 / 0.2e1;
t396 = t215 * mrSges(6,3);
t178 = -mrSges(6,2) * t292 + t396;
t436 = t178 / 0.2e1;
t463 = t167 * t435 + t333 * t436 - (t109 / 0.4e1 + t126 / 0.4e1) * t318 + (-t127 / 0.4e1 + t108 / 0.4e1) * t248 + t229 * t159 / 0.2e1 + t279 * t124 / 0.2e1 + t292 * t353 / 0.4e1;
t411 = Ifges(5,4) * t294;
t264 = -Ifges(5,2) * t291 + t411;
t413 = Ifges(5,1) * t291;
t328 = t411 + t413;
t424 = t294 / 0.2e1;
t462 = -t294 * t264 / 0.2e1 - t328 * t424;
t461 = -Ifges(3,4) - Ifges(4,6);
t317 = pkin(2) * t295 - t336;
t457 = m(4) * t317 - t295 * mrSges(4,2) + t292 * mrSges(4,3);
t157 = -t246 * t291 + t252;
t319 = t294 * t157 + t291 * t158;
t455 = -t248 * t62 + t318 * t63;
t289 = t294 ^ 2;
t451 = m(5) / 0.2e1;
t449 = m(6) / 0.2e1;
t446 = mrSges(5,2) / 0.2e1;
t445 = mrSges(6,3) / 0.2e1;
t442 = -qJ(3) / 0.2e1;
t441 = -t108 / 0.2e1;
t440 = -t109 / 0.2e1;
t410 = Ifges(6,4) * t248;
t162 = -Ifges(6,2) * t318 - t410;
t439 = t162 / 0.2e1;
t234 = Ifges(6,4) * t318;
t164 = -Ifges(6,1) * t248 - t234;
t438 = t164 / 0.2e1;
t437 = -t333 / 0.2e1;
t434 = -t215 / 0.2e1;
t433 = t215 / 0.2e1;
t432 = t216 / 0.2e1;
t431 = -t248 / 0.2e1;
t430 = t318 / 0.2e1;
t429 = -t318 / 0.2e1;
t352 = mrSges(5,3) * t366;
t388 = t292 * mrSges(5,1);
t254 = t352 + t388;
t428 = -t254 / 0.2e1;
t351 = mrSges(5,3) * t358;
t387 = t292 * mrSges(5,2);
t256 = -t351 - t387;
t427 = t256 / 0.2e1;
t426 = -t291 / 0.2e1;
t425 = t291 / 0.2e1;
t421 = pkin(4) * t294;
t420 = t58 * mrSges(6,2);
t419 = t59 * mrSges(6,1);
t416 = t64 * mrSges(6,1);
t415 = t65 * mrSges(6,2);
t412 = Ifges(5,4) * t291;
t409 = Ifges(5,5) * t292;
t407 = Ifges(5,6) * t292;
t406 = Ifges(5,6) * t294;
t403 = pkin(4) * qJD(4);
t397 = t215 * mrSges(6,1);
t395 = t216 * mrSges(6,2);
t389 = t290 * t63;
t386 = t293 * t62;
t385 = t295 * mrSges(5,1);
t384 = t295 * mrSges(5,2);
t125 = -t395 - t397;
t177 = -mrSges(6,2) * t295 - t214 * mrSges(6,3);
t179 = mrSges(6,1) * t295 - t217 * mrSges(6,3);
t228 = (-t421 - t443) * t292;
t367 = t291 * t292;
t253 = -mrSges(5,3) * t367 + t385;
t255 = mrSges(5,3) * t364 - t384;
t262 = t294 * mrSges(5,1) - t291 * mrSges(5,2);
t324 = -Ifges(5,5) * t291 - t406;
t325 = Ifges(6,4) * t217 - Ifges(6,2) * t214;
t327 = Ifges(6,1) * t217 - Ifges(6,4) * t214;
t330 = t214 * mrSges(6,1) + t217 * mrSges(6,2);
t3 = -t229 * t330 - t228 * t125 - t158 * t256 - t157 * t254 - t152 * t255 + t334 * t253 + t325 * t434 - t59 * t177 - t58 * t179 + t317 * (-mrSges(4,2) * t292 - mrSges(4,3) * t295) - t63 * t178 - t62 * t180 + t327 * t432 + t217 * t440 - t214 * t441 - m(5) * (t152 * t158 - t334 * t157 - t268 * t267) - m(6) * (t228 * t229 + t58 * t62 + t59 * t63) + (pkin(1) * mrSges(3,1) + (t324 - t461) * t292 + t316) * t292 + (Ifges(6,5) * t216 - Ifges(6,6) * t215 + pkin(1) * mrSges(3,2) + (t406 + t461) * t295 + (Ifges(5,2) * t289 - Ifges(3,1) + Ifges(3,2) - Ifges(4,2) + Ifges(4,3) - Ifges(5,3) - Ifges(6,3)) * t292 + (Ifges(5,5) * t295 + (0.2e1 * t411 + t413) * t292) * t291) * t295 + t457 * (pkin(2) * t292 - t380) + (t267 * t295 + t268 * t292) * t262;
t383 = t3 * qJD(1);
t354 = Ifges(6,5) * t215 + Ifges(6,6) * t216;
t314 = t59 * t394 + t229 * t124 + t292 * t354 / 0.2e1;
t331 = t291 * mrSges(5,1) + t294 * mrSges(5,2);
t6 = -m(6) * (t58 * t64 + t59 * t65) - t65 * t178 - t64 * t180 + (t127 / 0.2e1 + t441) * t216 + (t58 * mrSges(6,3) + t440 - t126 / 0.2e1) * t215 - t314 + (t351 + t256) * t334 + (t268 * t331 + (-Ifges(5,4) * t358 + t409) * t294 + (Ifges(5,4) * t366 - t407 + (-Ifges(5,1) + Ifges(5,2)) * t358 + (m(6) * t229 + t125) * pkin(4)) * t291) * t295 + (t254 - t352) * t152;
t382 = t6 * qJD(1);
t7 = -t180 * t59 + t108 * t432 - t216 * t127 / 0.2e1 + (t178 - t396) * t58 + t314 + (t109 + t126) * t433;
t381 = t7 * qJD(1);
t307 = (t216 * t430 + t248 * t433) * mrSges(6,3) + t178 * t431 + t180 * t429;
t15 = t307 - t458;
t377 = t15 * qJD(1);
t360 = t294 * t256;
t19 = -m(6) * (t214 * t59 + t217 * t58) - t214 * t178 - t217 * t180 - t254 * t367 + (t360 - m(5) * (-t152 * t294 - t334 * t291) - t457) * t292;
t376 = t19 * qJD(1);
t375 = t229 * t294;
t374 = t248 * t217;
t372 = t318 * t214;
t370 = t290 * t177;
t369 = t290 * t216;
t363 = t293 * t179;
t362 = t293 * t215;
t357 = t294 * t444;
t346 = Ifges(5,1) / 0.4e1 - Ifges(5,2) / 0.4e1;
t345 = -t58 / 0.2e1 + t65 / 0.2e1;
t344 = t59 / 0.2e1 + t64 / 0.2e1;
t341 = t367 / 0.2e1;
t340 = t125 * t424;
t339 = t444 * t427;
t163 = -Ifges(6,1) * t318 + t410;
t338 = t162 / 0.4e1 - t163 / 0.4e1;
t161 = Ifges(6,2) * t248 - t234;
t337 = t164 / 0.4e1 + t161 / 0.4e1;
t332 = mrSges(5,3) * (t291 ^ 2 / 0.2e1 + t289 / 0.2e1);
t326 = Ifges(5,2) * t294 + t412;
t20 = t279 * t159 - (t438 + t161 / 0.2e1) * t318 + (-t163 / 0.2e1 + t439) * t248;
t298 = (t167 * t445 + t338) * t216 + (mrSges(6,3) * t437 + t337) * t215 + t463;
t5 = t298 - t452;
t322 = t5 * qJD(1) + t20 * qJD(2);
t299 = t295 * t332 + (-t468 * t248 + t471 * t318) * t449 + t307;
t9 = (t427 - t387 / 0.2e1) * t294 + (t428 - t388 / 0.2e1) * t291 + t299 + t472;
t321 = t9 * qJD(1);
t302 = -m(5) * t268 / 0.2e1 - (t167 * t214 + t217 * t333 + t229) * m(6) / 0.2e1 + t397 / 0.2e1 + t395 / 0.2e1;
t304 = t177 * t430 + t179 * t431 + t319 * t451 + t455 * t449;
t13 = (t253 / 0.2e1 - t385 / 0.2e1) * t294 + (t255 / 0.2e1 + t384 / 0.2e1) * t291 + (t372 / 0.2e1 - t374 / 0.2e1) * mrSges(6,3) + t302 + t304;
t310 = -t331 + t335;
t91 = t422 + mrSges(4,3) + (m(5) + m(4)) * qJ(3) - t310;
t320 = -qJD(1) * t13 + qJD(2) * t91;
t315 = -t157 * mrSges(5,1) / 0.2e1 + t158 * t446;
t297 = t338 * t216 + t337 * t215 + (t167 * t432 + t344 * t248 - t318 * t345 + t333 * t434) * mrSges(6,3) + t268 * t262 / 0.2e1 + t463;
t300 = (t328 / 0.4e1 + t264 / 0.4e1 + t411 + mrSges(5,1) * t442 + t346 * t291 + (t335 / 0.2e1 - t422 / 0.2e1) * pkin(4)) * t295 - t444 * t428;
t266 = Ifges(5,1) * t294 - t412;
t305 = -t444 * t332 + (mrSges(5,2) * t442 - t266 / 0.4e1 + t326 / 0.4e1 - t346 * t294) * t294;
t309 = t471 * t167 + t468 * t333;
t1 = t297 + (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t305) * t295 + (-t339 - t407) * t294 + t309 * t449 + (t340 - t363 / 0.2e1 - t370 / 0.2e1 + 0.2e1 * (-t389 / 0.4e1 - t386 / 0.4e1 + t375 / 0.4e1) * m(6)) * pkin(4) + (t300 - t409) * t291 + t315 + t454;
t14 = qJ(3) * t262 + (-t266 / 0.2e1 + t326 / 0.2e1) * t291 + t20 + t462 + (t422 - t335) * t421;
t311 = t1 * qJD(1) + t14 * qJD(2);
t306 = (t290 * t435 + t293 * t436 + (-t362 / 0.2e1 + t369 / 0.2e1) * mrSges(6,3)) * pkin(4);
t11 = -t344 * mrSges(6,1) + t345 * mrSges(6,2) + t306;
t257 = (mrSges(6,1) * t290 + mrSges(6,2) * t293) * pkin(4);
t26 = (t437 + t333 / 0.2e1) * mrSges(6,2);
t308 = -t11 * qJD(1) - t26 * qJD(2) + t257 * qJD(4);
t247 = t257 * qJD(5);
t16 = t307 + t458;
t12 = -mrSges(5,2) * t366 / 0.2e1 + mrSges(5,1) * t358 / 0.2e1 - mrSges(6,3) * t372 / 0.2e1 + t374 * t445 + t255 * t425 + t253 * t424 + (m(4) * pkin(6) + mrSges(4,1)) * t295 - t302 + t304;
t10 = -t420 / 0.2e1 - t419 / 0.2e1 - t415 / 0.2e1 + t416 / 0.2e1 + t306 + t354;
t8 = t360 / 0.2e1 + t254 * t426 + mrSges(5,1) * t341 + t364 * t446 + t299 - t472;
t4 = t298 + t452;
t2 = t297 + (-t409 / 0.2e1 + t300) * t291 + (pkin(4) * t375 + t309) * t449 - t294 * t339 + pkin(4) * t340 + Ifges(5,5) * t341 + (t386 + t389) * t447 / 0.2e1 - t315 + (t363 + t370) * pkin(4) / 0.2e1 + (Ifges(5,3) / 0.2e1 + t305) * t295 + t452;
t17 = [-qJD(2) * t3 - qJD(3) * t19 - qJD(4) * t6 + qJD(5) * t7, t12 * qJD(3) + t2 * qJD(4) + t4 * qJD(5) - t383 + (-t228 * t335 + t333 * t179 + t167 * t177 - t214 * t439 + t217 * t438 - t253 * t357 - t255 * t365 - t267 * t331 + t279 * t330 + t325 * t429 + t327 * t431 + 0.2e1 * (-qJ(3) * t267 - t319 * t444) * t451 + 0.2e1 * (t167 * t63 + t228 * t279 + t333 * t62) * t449 + (t266 * t425 + t326 * t426 + Ifges(4,5) - Ifges(3,6) + (-mrSges(4,1) - t262) * qJ(3) + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(6) - t462) * t292 + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t294 - Ifges(6,5) * t248 - Ifges(5,6) * t291 - Ifges(6,6) * t318 - Ifges(4,4) + Ifges(3,5) + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2)) * pkin(6)) * t295 - t455 * mrSges(6,3) - t319 * mrSges(5,3)) * qJD(2), -t376 + t12 * qJD(2) + m(6) * (t372 - t374) * qJD(3) + t8 * qJD(4) + t16 * qJD(5), -t382 + t2 * qJD(2) + t8 * qJD(3) + (-t152 * mrSges(5,1) + t334 * mrSges(5,2) - Ifges(5,5) * t358 + Ifges(5,6) * t366 + t354 - t415 + t416) * qJD(4) + t10 * qJD(5) + (m(6) * (t290 * t65 + t293 * t64) + (-t362 + t369) * mrSges(6,3)) * t403, t381 + t4 * qJD(2) + t16 * qJD(3) + t10 * qJD(4) + (t354 - t419 - t420) * qJD(5); -qJD(3) * t13 + qJD(4) * t1 + qJD(5) * t5 + t383, qJD(3) * t91 + qJD(4) * t14 + qJD(5) * t20, t320, (mrSges(5,1) * t365 + mrSges(5,2) * t357 + t25 + t324) * qJD(4) + t470 + (m(6) * (-t167 * t293 + t290 * t333) + t465 * mrSges(6,3)) * t403 + t311, t25 * qJD(4) + t322 + t470; qJD(2) * t13 + qJD(4) * t9 + qJD(5) * t15 + t376, -t320, 0, (-t465 * t447 + t310) * qJD(4) + t469 + t321, qJD(4) * t335 + t377 + t469; -qJD(2) * t1 - qJD(3) * t9 + qJD(5) * t11 + t382, qJD(5) * t26 - t311, -t321, -t247, -t247 - t308; -qJD(2) * t5 - qJD(3) * t15 - qJD(4) * t11 - t381, -qJD(4) * t26 - t322, -t377, t308, 0;];
Cq = t17;
