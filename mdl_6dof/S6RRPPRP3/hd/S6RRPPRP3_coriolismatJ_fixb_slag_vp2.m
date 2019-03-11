% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:37
% EndTime: 2019-03-09 08:33:43
% DurationCPUTime: 4.16s
% Computational Cost: add. (5747->457), mult. (10352->566), div. (0->0), fcn. (8043->4), ass. (0->221)
t410 = -m(7) / 0.2e1;
t428 = Ifges(6,5) + Ifges(7,5);
t427 = Ifges(6,6) + Ifges(7,6);
t274 = cos(qJ(2));
t256 = t274 * qJ(4);
t271 = sin(qJ(5));
t386 = pkin(5) * t271;
t332 = -pkin(7) + t386;
t159 = -t274 * t332 - t256;
t259 = t271 * mrSges(7,1);
t273 = cos(qJ(5));
t261 = t273 * mrSges(7,2);
t313 = t261 + t259;
t438 = t159 * t313;
t268 = t271 ^ 2;
t269 = t273 ^ 2;
t343 = -t269 - t268;
t257 = t271 * mrSges(7,2);
t200 = t273 * mrSges(7,1) - t257;
t400 = -t200 / 0.2e1;
t258 = t271 * mrSges(6,2);
t416 = -t273 * mrSges(6,1) + t258;
t437 = t400 + t416 / 0.2e1;
t275 = -pkin(2) - pkin(3);
t267 = -pkin(8) + t275;
t196 = (qJ(6) - t267) * t271;
t356 = t273 * qJ(6);
t197 = t267 * t273 - t356;
t272 = sin(qJ(2));
t358 = t271 * t274;
t337 = mrSges(7,3) * t358;
t192 = -t272 * mrSges(7,2) + t337;
t353 = t273 * t274;
t194 = t272 * mrSges(7,1) + mrSges(7,3) * t353;
t394 = -t273 / 0.2e1;
t433 = t271 / 0.2e1;
t351 = t192 * t394 + t194 * t433;
t385 = t272 * pkin(5);
t338 = t274 * pkin(2) + t272 * qJ(3) + pkin(1);
t323 = t274 * pkin(3) + t338;
t101 = pkin(4) * t272 + pkin(8) * t274 + t323;
t253 = t272 * qJ(4);
t202 = pkin(7) * t272 - t253;
t63 = t273 * t101 - t202 * t271;
t54 = qJ(6) * t353 + t63;
t51 = t54 + t385;
t363 = t202 * t273;
t55 = t363 + (qJ(6) * t274 + t101) * t271;
t436 = ((t196 * t274 - t55) * t273 + (t197 * t274 + t51) * t271) * t410 - t351;
t270 = qJ(3) + pkin(4);
t241 = pkin(5) * t273 + t270;
t376 = t271 * mrSges(6,1);
t314 = t273 * mrSges(6,2) + t376;
t435 = -t241 * t313 - t270 * t314;
t388 = m(7) * t241;
t432 = -t388 / 0.2e1;
t431 = mrSges(7,3) + mrSges(6,3);
t367 = qJ(3) * t274;
t429 = -t272 * t275 - t367;
t426 = Ifges(6,3) + Ifges(7,3);
t408 = m(7) * pkin(5);
t334 = mrSges(7,1) + t408;
t425 = t270 * t274;
t288 = (t196 * t273 + t197 * t271) * t272;
t405 = mrSges(7,2) / 0.2e1;
t406 = mrSges(6,2) / 0.2e1;
t336 = t405 + t406;
t423 = t336 * t273;
t190 = t343 * t272;
t263 = Ifges(7,6) * t271;
t264 = Ifges(6,6) * t271;
t422 = -t263 - t264;
t359 = t271 * t272;
t292 = -t274 * mrSges(7,2) - mrSges(7,3) * t359;
t293 = -t274 * mrSges(6,2) - mrSges(6,3) * t359;
t421 = t293 + t292;
t357 = t272 * t273;
t294 = t274 * mrSges(7,1) - mrSges(7,3) * t357;
t295 = t274 * mrSges(6,1) - mrSges(6,3) * t357;
t420 = t295 + t294;
t418 = t427 * t353 + t358 * t428;
t105 = t267 * t272 + t425;
t205 = pkin(7) * t274 - t256;
t65 = t273 * t105 - t205 * t271;
t66 = t271 * t105 + t273 * t205;
t305 = -t271 * t65 + t273 * t66;
t52 = pkin(5) * t274 - t272 * t356 + t65;
t56 = -qJ(6) * t359 + t66;
t417 = -t271 * t52 + t273 * t56;
t193 = -t272 * mrSges(6,2) + mrSges(6,3) * t358;
t415 = -t193 / 0.2e1 + t336 * t272;
t414 = -0.2e1 * t190;
t413 = m(5) / 0.2e1;
t412 = m(6) / 0.2e1;
t411 = m(6) / 0.4e1;
t409 = m(7) / 0.2e1;
t407 = -mrSges(7,1) / 0.2e1;
t404 = mrSges(7,3) / 0.2e1;
t401 = t192 / 0.2e1;
t399 = t205 / 0.2e1;
t398 = -t267 / 0.2e1;
t397 = -t271 / 0.2e1;
t393 = t273 / 0.2e1;
t391 = t274 / 0.2e1;
t390 = m(6) * t270;
t389 = m(7) * t159;
t387 = m(7) * t274;
t384 = mrSges(5,1) + mrSges(4,3);
t383 = -mrSges(4,2) + mrSges(5,3);
t382 = Ifges(6,4) * t271;
t381 = Ifges(6,4) * t273;
t380 = Ifges(7,4) * t271;
t379 = Ifges(7,4) * t273;
t277 = (t271 * t66 + t273 * t65) * t412 + (t271 * t56 + t273 * t52) * t409 + t421 * t433 + t420 * t393;
t299 = t196 * t271 - t197 * t273;
t289 = m(7) * t299;
t324 = t267 * t343;
t344 = t274 * mrSges(5,1) + t272 * mrSges(5,2);
t5 = 0.2e1 * (t289 / 0.4e1 + t324 * t411) * t272 - t277 - t344 - (t404 + mrSges(6,3) / 0.2e1) * t190 + m(5) * t429 + (t432 - t390 / 0.2e1 + t437) * t274;
t378 = qJD(1) * t5;
t158 = t272 * t332 + t253;
t164 = t313 * t274;
t195 = t272 * mrSges(6,1) + mrSges(6,3) * t353;
t304 = m(5) * t323;
t306 = pkin(2) * t272 - t367;
t316 = (Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1) * t272;
t322 = m(4) * t338;
t247 = Ifges(7,4) * t358;
t123 = -Ifges(7,1) * t353 + Ifges(7,5) * t272 + t247;
t248 = Ifges(6,4) * t358;
t124 = -Ifges(6,1) * t353 + Ifges(6,5) * t272 + t248;
t326 = -t123 / 0.2e1 - t124 / 0.2e1;
t307 = -Ifges(7,2) * t271 + t379;
t121 = Ifges(7,6) * t272 - t274 * t307;
t308 = -Ifges(6,2) * t271 + t381;
t122 = Ifges(6,6) * t272 - t274 * t308;
t327 = t121 / 0.2e1 + t122 / 0.2e1;
t339 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t64 = t101 * t271 + t363;
t1 = -t56 * t192 - t55 * t292 - t52 * t194 - t51 * t294 - t63 * t295 - t64 * t293 - t65 * t195 - t66 * t193 + t158 * t164 - t323 * t344 + t306 * t322 - m(6) * (-t202 * t205 + t63 * t65 + t64 * t66) - m(7) * (t158 * t159 + t51 * t52 + t55 * t56) + (t306 * mrSges(4,1) + pkin(1) * mrSges(3,2) - t429 * mrSges(5,2) - t338 * mrSges(4,3) - qJ(3) * t304 - t202 * t314 + (t273 * t428 - t339 + t422) * t274) * t274 + (-t438 - t205 * t314 + t429 * mrSges(5,1) + t338 * mrSges(4,1) + t306 * mrSges(4,3) - t275 * t304 + pkin(1) * mrSges(3,1) + t326 * t273 + (t316 + t327) * t271 + (-Ifges(3,1) - Ifges(4,1) + Ifges(5,1) + Ifges(3,2) - Ifges(5,2) + Ifges(4,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t269 + ((Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t271 + (-Ifges(6,4) - Ifges(7,4)) * t273) * t271 - t426) * t274 + (t339 + (-Ifges(6,5) / 0.2e1 - Ifges(7,5) / 0.2e1) * t273) * t272) * t272;
t377 = t1 * qJD(1);
t242 = mrSges(7,2) * t358;
t162 = -mrSges(7,1) * t353 + t242;
t163 = t416 * t274;
t165 = Ifges(7,2) * t353 + t247;
t166 = Ifges(6,2) * t353 + t248;
t309 = Ifges(7,1) * t271 + t379;
t167 = t309 * t274;
t311 = Ifges(6,1) * t271 + t381;
t168 = t311 * t274;
t333 = m(7) * (-t51 + t54);
t4 = t159 * t162 + t54 * t192 + t63 * t193 - t64 * t195 + t205 * t163 + (-t194 + t333) * t55 + ((t165 / 0.2e1 + t166 / 0.2e1 - t51 * mrSges(7,3) - t63 * mrSges(6,3) - t326) * t271 + (-t167 / 0.2e1 - t168 / 0.2e1 + t55 * mrSges(7,3) + t64 * mrSges(6,3) + (t164 - t389) * pkin(5) + t327) * t273) * t274 + t418 * t272 / 0.2e1;
t370 = t4 * qJD(1);
t321 = t407 - t408 / 0.2e1;
t298 = mrSges(6,1) / 0.2e1 - t321;
t282 = -t333 / 0.2e1 + t195 / 0.2e1 + t298 * t272;
t325 = t269 / 0.2e1 + t268 / 0.2e1;
t6 = (t401 - t415) * t271 + (t194 / 0.2e1 + t282) * t273 - t431 * t274 * t325;
t369 = t6 * qJD(1);
t9 = t282 * t271 + t273 * t415 + t351;
t368 = t9 * qJD(1);
t348 = t194 + t195;
t349 = t192 + t193;
t362 = t205 * t274;
t11 = (t164 + (mrSges(5,3) + t314) * t274) * t274 + (t272 * mrSges(5,3) + t348 * t271 - t349 * t273) * t272 + m(7) * (-t159 * t274 + (t271 * t51 - t273 * t55) * t272) + m(6) * (-t362 + (t271 * t63 - t273 * t64) * t272) + m(5) * (-t202 * t272 - t362);
t366 = qJD(1) * t11;
t291 = m(7) * (t271 * t55 + t51 * t273);
t14 = (t304 + t322 + t384 * t272 + (-mrSges(5,2) + mrSges(4,1)) * t274 + t348 * t273 + t349 * t271 + t291 + m(6) * (t271 * t64 + t273 * t63)) * t272;
t365 = qJD(1) * t14;
t354 = t273 * t194;
t361 = t271 * t192;
t23 = (t291 + t354 + t361) * t274;
t364 = qJD(1) * t23;
t33 = m(5) * t272 + (t409 + t412) * t414;
t352 = t33 * qJD(1);
t345 = t257 + t258;
t342 = qJD(5) * t273;
t318 = -0.1e1 / 0.2e1 - t325;
t103 = t318 * t387;
t341 = t103 * qJD(1);
t188 = t318 * m(7);
t340 = t188 * qJD(2);
t335 = -t54 / 0.2e1 + t51 / 0.2e1;
t169 = -m(7) * t386 - t313;
t320 = t333 / 0.2e1;
t319 = -mrSges(6,1) * t267 - t428;
t315 = t325 * mrSges(6,3);
t312 = Ifges(6,1) * t273 - t382;
t310 = Ifges(7,1) * t273 - t380;
t206 = -Ifges(7,2) * t273 - t380;
t208 = -Ifges(6,2) * t273 - t382;
t16 = (t309 / 0.2e1 + t307 / 0.2e1 + t311 / 0.2e1 + t308 / 0.2e1) * t273 + (t310 / 0.2e1 + t206 / 0.2e1 + t312 / 0.2e1 + t208 / 0.2e1 + (-t200 - t388) * pkin(5)) * t271 + t435;
t278 = (-t194 / 0.2e1 + t320) * t197 + (t264 / 0.4e1 + t263 / 0.4e1) * t272 - t438 / 0.2e1 + t196 * t401 - t314 * t399 + t241 * t162 / 0.2e1 + t270 * t163 / 0.2e1;
t279 = -t168 / 0.4e1 - t167 / 0.4e1 + t122 / 0.4e1 + t121 / 0.4e1 + t193 * t398 + (-t389 / 0.2e1 + t164 / 0.2e1) * pkin(5);
t280 = t197 * t404 + t312 / 0.4e1 + t310 / 0.4e1 + t208 / 0.4e1 + t206 / 0.4e1 + (t400 + t432) * pkin(5);
t285 = -t166 / 0.4e1 - t165 / 0.4e1 - t124 / 0.4e1 - t123 / 0.4e1 + t195 * t398;
t286 = -t196 * mrSges(7,3) / 0.2e1 - t311 / 0.4e1 - t309 / 0.4e1 - t308 / 0.4e1 - t307 / 0.4e1;
t287 = t56 * t405 - t65 * mrSges(6,1) / 0.2e1 + t66 * t406;
t296 = t267 * t315;
t3 = t321 * t52 + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + pkin(5) * t407 + t296) * t274 + (t274 * t286 + t279 + t316) * t271 + ((-0.3e1 / 0.4e1 * Ifges(6,5) - 0.3e1 / 0.4e1 * Ifges(7,5)) * t272 + (t385 / 0.2e1 + t335) * mrSges(7,3) + t280 * t274 + t285) * t273 + t278 + t287;
t303 = t3 * qJD(1) + t16 * qJD(2);
t281 = -m(6) * t305 / 0.2e1 + t417 * t410;
t13 = t389 / 0.2e1 + m(6) * t399 + t288 * t409 + t281;
t60 = (mrSges(7,1) + mrSges(6,1)) * t273 + (m(5) + m(4)) * qJ(3) + t390 + t388 - t345 + t384;
t302 = qJD(1) * t13 + qJD(2) * t60;
t284 = (t261 / 0.2e1 + t259 / 0.2e1) * t272 + t158 * t409;
t18 = t284 + t436;
t59 = -mrSges(7,3) * t343 + t289;
t301 = -qJD(1) * t18 + qJD(2) * t59;
t102 = t334 * t353 - t242;
t297 = qJD(1) * t102 - qJD(2) * t169;
t187 = (t343 + 0.1e1) * t409;
t104 = (-0.1e1 / 0.2e1 - t343 / 0.2e1) * t387;
t34 = (m(7) / 0.4e1 + t411) * t414 + (m(7) + m(6)) * t190 / 0.2e1;
t19 = t284 - t436;
t12 = (-t256 + t288) * t409 + 0.2e1 * (t411 + t413) * t205 + (-t259 / 0.2e1 + t332 * t410 - t376 / 0.2e1 + m(4) * pkin(7) - t423 - t383) * t274 - t281 + t420 * t397 + t421 * t393;
t10 = t271 * t320 + t193 * t393 + t195 * t397 + (t271 * t298 + t423) * t272 - t351;
t8 = t437 * t274 + (mrSges(7,3) * t325 + t315) * t272 + (-t241 * t274 + t272 * t299) * t409 + (t272 * t324 - t425) * t412 + t277;
t7 = -t354 / 0.2e1 - t361 / 0.2e1 + t273 * t320 + t193 * t397 + t195 * t394 + (-t271 * t336 + t273 * t298) * t272 - t431 * t343 * t391;
t2 = ((-Ifges(6,5) / 0.4e1 - Ifges(7,5) / 0.4e1) * t272 + t335 * mrSges(7,3) + t285) * t273 + t278 + pkin(5) * t294 / 0.2e1 + t279 * t271 + (t286 * t271 + t280 * t273 + t296) * t274 - t287 + t334 * t52 / 0.2e1 + t426 * t391 - t427 * t359 / 0.2e1 + t428 * t357 / 0.2e1;
t15 = [-qJD(2) * t1 + qJD(3) * t14 + qJD(4) * t11 + qJD(5) * t4 + qJD(6) * t23, t12 * qJD(3) + t8 * qJD(4) + t2 * qJD(5) + t19 * qJD(6) - t377 + (-t202 * mrSges(5,1) + t205 * mrSges(5,2) + t158 * t200 + t202 * t416 + 0.2e1 * (-qJ(3) * t202 + t205 * t275) * t413 + 0.2e1 * (-t202 * t270 + t267 * t305) * t412 + 0.2e1 * (t158 * t241 + t196 * t52 + t197 * t56) * t409 + (t196 * mrSges(7,1) - pkin(2) * mrSges(4,2) - t197 * mrSges(7,2) - t275 * mrSges(5,3) + Ifges(4,4) + Ifges(3,5) + Ifges(5,6) + (-mrSges(6,2) * t267 - t427) * t273 + t319 * t271 + (-m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1)) * pkin(7)) * t274 + (-t417 - t288) * mrSges(7,3) - t305 * mrSges(6,3) + (-Ifges(5,5) + Ifges(4,6) - Ifges(3,6) + t383 * qJ(3) + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(7) + (t307 + t308) * t394 + (-t311 - t309) * t393 + (t310 + t312 + t206 + t208) * t397 - t435) * t272) * qJD(2), qJD(2) * t12 + qJD(4) * t34 + qJD(5) * t7 + t365, qJD(2) * t8 + qJD(3) * t34 + qJD(5) * t10 + qJD(6) * t104 + t366, t370 + t2 * qJD(2) + t7 * qJD(3) + t10 * qJD(4) + (-mrSges(6,1) * t64 - mrSges(7,1) * t55 - mrSges(6,2) * t63 - mrSges(7,2) * t54 + (-m(7) * t55 - t337) * pkin(5) + t418) * qJD(5), qJD(2) * t19 + qJD(4) * t104 + t364; qJD(3) * t13 + qJD(4) * t5 + qJD(5) * t3 - qJD(6) * t18 + t377, qJD(3) * t60 + qJD(5) * t16 + qJD(6) * t59, qJD(6) * t187 + t302, t378 (-mrSges(7,2) * t196 - t197 * t334 + t267 * t258 - t422) * qJD(5) + (pkin(5) * mrSges(7,3) + t319) * t342 + t303, qJD(3) * t187 + t301; -qJD(2) * t13 - qJD(4) * t33 - qJD(5) * t6 - t365, qJD(6) * t188 - t302, 0, -t352, -t369 + t345 * qJD(5) + (-mrSges(6,1) - t334) * t342, t340; -qJD(2) * t5 + qJD(3) * t33 - qJD(5) * t9 - qJD(6) * t103 - t366, -t378, t352, 0, -t368 + (-t314 + t169) * qJD(5), -t341; -qJD(2) * t3 + qJD(3) * t6 + qJD(4) * t9 + qJD(6) * t102 - t370, -t169 * qJD(6) - t303, t369, t368, 0, t297; qJD(2) * t18 + qJD(4) * t103 - qJD(5) * t102 - t364, -qJD(3) * t188 + qJD(5) * t169 - t301, -t340, t341, -t297, 0;];
Cq  = t15;
