% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:11
% EndTime: 2019-12-31 21:16:23
% DurationCPUTime: 5.67s
% Computational Cost: add. (16088->431), mult. (32578->586), div. (0->0), fcn. (35615->8), ass. (0->265)
t418 = sin(qJ(3));
t419 = sin(qJ(2));
t421 = cos(qJ(3));
t422 = cos(qJ(2));
t260 = t418 * t419 - t421 * t422;
t286 = sin(pkin(9));
t287 = cos(pkin(9));
t417 = sin(qJ(5));
t420 = cos(qJ(5));
t453 = t417 * t286 - t420 * t287;
t183 = t453 * t260;
t405 = t183 * mrSges(6,2);
t259 = -t420 * t286 - t417 * t287;
t180 = t259 * t260;
t407 = t180 * mrSges(6,1);
t445 = -m(6) / 0.2e1;
t447 = -m(5) / 0.2e1;
t354 = t422 * pkin(6);
t268 = t422 * pkin(7) + t354;
t351 = t419 * pkin(6);
t312 = -t419 * pkin(7) - t351;
t455 = t421 * t268 + t418 * t312;
t372 = t260 * t286;
t457 = -pkin(4) * t372 + t455;
t484 = -t405 / 0.2e1 - t407 / 0.2e1 + t445 * t457 + t447 * t455;
t261 = -t418 * t422 - t421 * t419;
t474 = t453 * t261;
t435 = t474 / 0.2e1;
t308 = t259 * t261;
t130 = -mrSges(6,2) * t260 - mrSges(6,3) * t308;
t132 = mrSges(6,1) * t260 - mrSges(6,3) * t474;
t428 = t259 / 0.2e1;
t429 = -t453 / 0.2e1;
t446 = m(5) / 0.2e1;
t279 = -t422 * pkin(2) - pkin(1);
t199 = t260 * pkin(3) + t261 * qJ(4) + t279;
t116 = t287 * t199 - t286 * t455;
t117 = t286 * t199 + t287 * t455;
t451 = -t286 * t116 + t117 * t287;
t370 = t261 * t286;
t206 = -t260 * mrSges(5,2) + mrSges(5,3) * t370;
t369 = t261 * t287;
t208 = t260 * mrSges(5,1) + mrSges(5,3) * t369;
t423 = t287 / 0.2e1;
t454 = t206 * t423 - t286 * t208 / 0.2e1;
t297 = (-t308 * t429 - t428 * t474) * mrSges(6,3) + t451 * t446 + t130 * t429 + t132 * t428 + t454;
t392 = t287 * mrSges(5,2);
t394 = t286 * mrSges(5,1);
t483 = (t392 / 0.2e1 + t394 / 0.2e1) * t260 + t297 + t484;
t371 = t260 * t287;
t338 = -t371 / 0.2e1;
t482 = -mrSges(5,1) * t372 / 0.2e1 + mrSges(5,2) * t338 + t297 - t484;
t441 = -mrSges(6,2) / 0.2e1;
t442 = mrSges(6,1) / 0.2e1;
t436 = t183 / 0.2e1;
t438 = -t180 / 0.2e1;
t314 = Ifges(6,5) * t436 + Ifges(6,6) * t438;
t452 = Ifges(6,3) * t261 / 0.2e1 - t314;
t217 = -t261 * pkin(3) + t260 * qJ(4);
t231 = t418 * t268 - t421 * t312;
t462 = t286 * t231;
t127 = t287 * t217 + t462;
t326 = -t261 * pkin(4) + pkin(8) * t371;
t77 = t127 + t326;
t461 = t287 * t231;
t128 = t286 * t217 - t461;
t360 = pkin(8) * t372;
t94 = t360 + t128;
t48 = -t417 * t94 + t420 * t77;
t49 = t417 * t77 + t420 * t94;
t481 = -t49 * t441 - t48 * t442 + t452;
t352 = t419 * pkin(2);
t200 = t352 + t217;
t122 = t287 * t200 + t462;
t66 = t122 + t326;
t123 = t286 * t200 - t461;
t88 = t360 + t123;
t43 = -t417 * t88 + t420 * t66;
t44 = t417 * t66 + t420 * t88;
t480 = -t43 * t442 - t44 * t441 + t452;
t409 = Ifges(5,2) * t286;
t412 = Ifges(5,4) * t287;
t139 = -Ifges(5,6) * t261 + (t409 - t412) * t260;
t413 = Ifges(5,4) * t286;
t140 = -Ifges(5,5) * t261 + (-Ifges(5,1) * t287 + t413) * t260;
t424 = t286 / 0.2e1;
t252 = Ifges(6,4) * t453;
t214 = -Ifges(6,1) * t259 - t252;
t432 = t214 / 0.2e1;
t410 = Ifges(6,4) * t259;
t212 = -Ifges(6,2) * t453 - t410;
t433 = t212 / 0.2e1;
t69 = Ifges(6,4) * t183 - Ifges(6,2) * t180 - t261 * Ifges(6,6);
t71 = Ifges(6,1) * t183 - Ifges(6,4) * t180 - t261 * Ifges(6,5);
t302 = -t259 * t71 / 0.2e1 + t69 * t429 - t180 * t433 + t183 * t432 + t140 * t424 + t139 * t423 + (Ifges(5,2) * t287 + t413) * t372 / 0.2e1 + (Ifges(5,1) * t286 + t412) * t338 + Ifges(4,6) * t261 - Ifges(4,5) * t260 - (Ifges(5,5) * t286 - Ifges(6,5) * t259 + Ifges(5,6) * t287 - Ifges(6,6) * t453) * t261 / 0.2e1;
t266 = -mrSges(5,1) * t287 + t286 * mrSges(5,2);
t458 = t455 * t266;
t464 = t455 * mrSges(4,1);
t466 = t231 * mrSges(4,2);
t210 = mrSges(6,1) * t453 - mrSges(6,2) * t259;
t473 = t457 * t210;
t479 = t302 + t458 + t466 - t464 + t473;
t478 = -t458 / 0.2e1 + t464 / 0.2e1 - t466 / 0.2e1 - t473 / 0.2e1;
t160 = -pkin(4) * t370 + t231;
t477 = t160 * t457;
t353 = t421 * pkin(2);
t278 = -t353 - pkin(3);
t415 = t287 * pkin(4);
t264 = t278 - t415;
t476 = t264 * t457;
t274 = -pkin(3) - t415;
t475 = t274 * t457;
t129 = mrSges(6,2) * t261 - mrSges(6,3) * t180;
t131 = -mrSges(6,1) * t261 - mrSges(6,3) * t183;
t325 = -t392 - t394;
t196 = t325 * t260;
t197 = t325 * t261;
t205 = mrSges(5,2) * t261 + mrSges(5,3) * t372;
t207 = -t261 * mrSges(5,1) + mrSges(5,3) * t371;
t65 = t260 * pkin(4) + pkin(8) * t369 + t116;
t83 = pkin(8) * t370 + t117;
t40 = -t417 * t83 + t420 * t65;
t41 = t417 * t65 + t420 * t83;
t437 = -t308 / 0.2e1;
t411 = Ifges(6,4) * t474;
t70 = -Ifges(6,2) * t308 + t260 * Ifges(6,6) + t411;
t175 = Ifges(6,4) * t308;
t72 = Ifges(6,1) * t474 + t260 * Ifges(6,5) - t175;
t90 = t405 + t407;
t91 = mrSges(6,1) * t308 + mrSges(6,2) * t474;
t472 = t457 * t91 + t455 * t197 + t116 * t207 + t117 * t205 + t41 * t129 + t40 * t131 + t160 * t90 + t231 * t196 + t279 * (-mrSges(4,1) * t261 - mrSges(4,2) * t260) + t71 * t435 + t72 * t436 + t69 * t437 + t70 * t438;
t350 = t418 * pkin(2);
t468 = t350 / 0.2e1;
t467 = pkin(3) * t455;
t285 = t287 ^ 2;
t362 = t286 ^ 2 + t285;
t465 = t362 * mrSges(5,3);
t374 = t231 * t455;
t463 = t278 * t455;
t460 = t418 * t231;
t209 = -t259 * mrSges(6,1) - mrSges(6,2) * t453;
t51 = 0.2e1 * t435 * mrSges(6,1) + 0.2e1 * t437 * mrSges(6,2);
t448 = qJD(1) * t51 + (qJD(2) + qJD(3)) * t209;
t444 = m(6) / 0.2e1;
t443 = m(5) * pkin(2);
t440 = -mrSges(6,3) / 0.2e1;
t439 = t130 / 0.2e1;
t273 = t350 + qJ(4);
t246 = (-pkin(8) - t273) * t286;
t282 = t287 * pkin(8);
t247 = t273 * t287 + t282;
t191 = t417 * t246 + t420 * t247;
t434 = -t191 / 0.2e1;
t265 = (-pkin(8) - qJ(4)) * t286;
t267 = qJ(4) * t287 + t282;
t229 = t417 * t265 + t420 * t267;
t431 = -t229 / 0.2e1;
t327 = t421 * t417;
t328 = t421 * t420;
t237 = (-t286 * t328 - t287 * t327) * pkin(2);
t430 = t237 / 0.2e1;
t426 = t264 / 0.2e1;
t425 = t274 / 0.2e1;
t416 = pkin(3) * t196;
t414 = t259 * t40 - t41 * t453;
t406 = t308 * mrSges(6,2);
t404 = t474 * mrSges(6,1);
t391 = t287 * Ifges(5,5);
t393 = t286 * Ifges(5,6);
t295 = (-Ifges(5,3) + Ifges(4,1) - Ifges(4,2) - Ifges(6,3) + t285 * Ifges(5,1) / 0.2e1 + (-t412 + t409 / 0.2e1) * t286) * t261 + t314 + (Ifges(4,4) - t391 + t393) * t260;
t299 = -Ifges(6,5) * t474 / 0.2e1 + Ifges(6,6) * t308 / 0.2e1 + t139 * t424 - t287 * t140 / 0.2e1 + (t391 / 0.2e1 - t393 / 0.2e1 - Ifges(4,4)) * t261;
t2 = -pkin(1) * (t419 * mrSges(3,1) + t422 * mrSges(3,2)) + m(6) * (t40 * t43 + t41 * t44 + t477) + m(5) * (t116 * t122 + t117 * t123 + t374) + m(4) * t279 * t352 + (mrSges(4,1) * t352 + t295) * t260 + t123 * t206 + t122 * t208 + t44 * t130 + t43 * t132 + (-mrSges(4,2) * t352 + t299) * t261 + (-Ifges(3,2) + Ifges(3,1)) * t422 * t419 + (-t419 ^ 2 + t422 ^ 2) * Ifges(3,4) + t472;
t403 = t2 * qJD(1);
t398 = t453 * mrSges(6,3);
t397 = t259 * mrSges(6,3);
t396 = t264 * t90;
t395 = t274 * t90;
t4 = t299 * t261 + m(6) * (t40 * t48 + t41 * t49 + t477) + m(5) * (t116 * t127 + t117 * t128 + t374) + t295 * t260 + t128 * t206 + t127 * t208 + t49 * t130 + t48 * t132 + t472;
t390 = t4 * qJD(1);
t364 = -Ifges(6,5) * t308 - Ifges(6,6) * t474;
t89 = t404 - t406;
t92 = -Ifges(6,2) * t474 - t175;
t93 = -Ifges(6,1) * t308 - t411;
t9 = t160 * t89 + t260 * t364 / 0.2e1 + t40 * t130 - t41 * t132 + (-t41 * mrSges(6,3) + t93 / 0.2e1 - t70 / 0.2e1) * t474 - (-t40 * mrSges(6,3) + t72 / 0.2e1 + t92 / 0.2e1) * t308;
t389 = t9 * qJD(1);
t15 = m(6) * (-t308 * t41 - t40 * t474) - t308 * t130 - t474 * t132 + (t286 * t206 + t287 * t208 + m(5) * (t116 * t287 + t117 * t286)) * t261;
t388 = qJD(1) * t15;
t385 = t122 * t286;
t384 = t123 * t287;
t383 = t127 * t286;
t382 = t128 * t287;
t190 = t420 * t246 - t417 * t247;
t379 = t190 * t131;
t378 = t191 * t129;
t228 = t420 * t265 - t417 * t267;
t377 = t228 * t131;
t376 = t229 * t129;
t368 = t278 * t196;
t366 = t286 * t207;
t365 = t287 * t205;
t363 = -Ifges(6,5) * t453 + Ifges(6,6) * t259;
t359 = t43 * t397;
t358 = t44 * t398;
t357 = t48 * t397;
t356 = t49 * t398;
t349 = mrSges(5,3) * t385;
t348 = mrSges(5,3) * t384;
t347 = mrSges(5,3) * t383;
t346 = mrSges(5,3) * t382;
t345 = qJ(4) * t366;
t344 = t273 * t366;
t343 = qJ(4) * t365;
t342 = t273 * t365;
t211 = Ifges(6,2) * t259 - t252;
t335 = t211 / 0.4e1 + t214 / 0.4e1;
t213 = -Ifges(6,1) * t453 + t410;
t334 = -t212 / 0.4e1 + t213 / 0.4e1;
t333 = t362 * qJ(4);
t332 = t362 * t273;
t303 = (-t213 / 0.2e1 + t433) * t259 - (t432 + t211 / 0.2e1) * t453;
t23 = t264 * t209 + t303;
t300 = (t70 / 0.4e1 - t93 / 0.4e1) * t259 - (t92 / 0.4e1 + t72 / 0.4e1) * t453 + t160 * t209 / 0.2e1 + t260 * t363 / 0.4e1;
t290 = -(t190 * t440 + t335) * t308 + (mrSges(6,3) * t434 + t334) * t474 + t190 * t439 + t132 * t434 + t89 * t426 + t300;
t6 = t290 + t480;
t324 = t6 * qJD(1) + t23 * qJD(2);
t238 = (-t286 * t327 + t287 * t328) * pkin(2);
t293 = t237 * t397 - t238 * t398 + (-mrSges(4,1) + t210 + t266) * t350 + (-mrSges(4,2) + t465) * t353;
t311 = t362 * t421;
t26 = m(6) * (t190 * t237 + t191 * t238 + t264 * t350) + (t311 * t273 + t418 * t278) * t443 + t293;
t320 = t382 - t383;
t288 = (t463 + t320 * t273 + (t451 * t421 + t460) * pkin(2)) * t446 + (t160 * t350 + t190 * t48 + t191 * t49 + t237 * t40 + t238 * t41 + t476) * t444 + t379 / 0.2e1 + t378 / 0.2e1 + t132 * t430 + t238 * t439 + t396 / 0.2e1 + t368 / 0.2e1 - t347 / 0.2e1 + t346 / 0.2e1 - t344 / 0.2e1 + t342 / 0.2e1 + t357 / 0.2e1 - t356 / 0.2e1 + t454 * t353 + (t197 + t91) * t468 - t478;
t321 = t384 - t385;
t289 = (t321 * qJ(4) - t467) * t447 + (t228 * t43 + t229 * t44 + t475) * t445 + t416 / 0.2e1 - t377 / 0.2e1 - t376 / 0.2e1 - t395 / 0.2e1 + t349 / 0.2e1 - t348 / 0.2e1 + t345 / 0.2e1 - t343 / 0.2e1 - t359 / 0.2e1 + t358 / 0.2e1 + t478;
t3 = t289 + t288;
t323 = t3 * qJD(1) + t26 * qJD(2);
t306 = (-t190 * t474 - t191 * t308 + t414) * t444;
t11 = t306 + t483;
t319 = (t259 ^ 2 + t453 ^ 2) * mrSges(6,3) + t465;
t45 = m(6) * (t190 * t259 - t191 * t453) + m(5) * t332 + t319;
t322 = qJD(1) * t11 + qJD(2) * t45;
t313 = mrSges(6,1) * t430 + t238 * t441;
t301 = (t426 + t425) * t209 + t303;
t18 = t301 - t313;
t24 = t274 * t209 + t303;
t291 = (mrSges(6,3) * t431 + t334) * t474 - (t228 * t440 + t335) * t308 + t228 * t439 + t132 * t431 + t89 * t425 + t300;
t8 = t291 + t481;
t310 = t8 * qJD(1) + t18 * qJD(2) + t24 * qJD(3);
t305 = (-t228 * t474 - t229 * t308 + t414) * t444;
t13 = t305 + t483;
t298 = (t332 + t333) * t447 + ((t190 + t228) * t259 - (t191 + t229) * t453) * t445 - t319;
t307 = (m(5) + m(6)) * t468;
t27 = t307 + t298;
t53 = m(6) * (t228 * t259 - t229 * t453) + m(5) * t333 + t319;
t309 = qJD(1) * t13 - qJD(2) * t27 + qJD(3) * t53;
t204 = t209 * qJD(4);
t203 = t209 * qJD(5);
t52 = -t406 / 0.2e1 + t404 / 0.2e1 - t308 * t441 - t474 * t442;
t28 = t307 - t298;
t19 = t301 + t313;
t12 = t305 + t482;
t10 = t306 + t482;
t7 = t291 - t481;
t5 = t290 - t480;
t1 = -t289 + t288 + t302;
t14 = [qJD(2) * t2 + qJD(3) * t4 + qJD(4) * t15 + qJD(5) * t9, t403 + (t379 + t378 + t348 - t344 - t349 + m(6) * (t190 * t43 + t191 * t44 + t476) + m(5) * (t273 * t321 + t463) + (t260 * t353 + t261 * t350) * mrSges(4,3) + t342 + t396 + m(4) * (-t421 * t455 - t460) * pkin(2) + Ifges(3,5) * t422 - Ifges(3,6) * t419 - mrSges(3,1) * t354 - t358 + t359 + t368 + mrSges(3,2) * t351 + t479) * qJD(2) + t1 * qJD(3) + t10 * qJD(4) + t5 * qJD(5), t390 + t1 * qJD(2) + (t343 - t345 + t346 - t347 - t356 + t357 + t376 + t377 + t395 - t416 + t479) * qJD(3) + t12 * qJD(4) + t7 * qJD(5) + 0.2e1 * ((t228 * t48 + t229 * t49 + t475) * t444 + (qJ(4) * t320 - t467) * t446) * qJD(3), qJD(2) * t10 + qJD(3) * t12 + qJD(5) * t52 + t388, t389 + t5 * qJD(2) + t7 * qJD(3) + t52 * qJD(4) + (-mrSges(6,1) * t41 - mrSges(6,2) * t40 + t364) * qJD(5); qJD(3) * t3 + qJD(4) * t11 + qJD(5) * t6 - t403, qJD(3) * t26 + qJD(4) * t45 + qJD(5) * t23, (m(6) * (t228 * t237 + t229 * t238 + t274 * t350) + (-t418 * pkin(3) + t311 * qJ(4)) * t443 + t293) * qJD(3) + t28 * qJD(4) + t19 * qJD(5) + t323, qJD(3) * t28 + t322, t19 * qJD(3) + (-mrSges(6,1) * t191 - mrSges(6,2) * t190 + t363) * qJD(5) + t324; -qJD(2) * t3 + qJD(4) * t13 + qJD(5) * t8 - t390, -qJD(4) * t27 + qJD(5) * t18 - t323, qJD(4) * t53 + qJD(5) * t24, t309, (-mrSges(6,1) * t229 - mrSges(6,2) * t228 + t363) * qJD(5) + t310; -qJD(2) * t11 - qJD(3) * t13 + qJD(5) * t51 - t388, qJD(3) * t27 + t203 - t322, t203 - t309, 0, t448; -qJD(2) * t6 - qJD(3) * t8 - qJD(4) * t51 - t389, -qJD(3) * t18 - t204 - t324, -t204 - t310, -t448, 0;];
Cq = t14;
