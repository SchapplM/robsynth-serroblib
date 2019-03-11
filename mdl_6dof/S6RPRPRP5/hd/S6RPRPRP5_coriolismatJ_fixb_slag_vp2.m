% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:35
% EndTime: 2019-03-09 03:14:49
% DurationCPUTime: 6.94s
% Computational Cost: add. (19181->546), mult. (39230->740), div. (0->0), fcn. (44777->8), ass. (0->281)
t304 = sin(pkin(9));
t306 = cos(pkin(9));
t308 = sin(qJ(3));
t458 = cos(qJ(3));
t292 = t304 * t458 + t308 * t306;
t303 = sin(pkin(10));
t305 = cos(pkin(10));
t307 = sin(qJ(5));
t457 = cos(qJ(5));
t504 = -t307 * t303 + t457 * t305;
t319 = t504 * t292;
t517 = -t319 / 0.2e1;
t291 = t303 * t457 + t307 * t305;
t511 = t291 * t292;
t526 = -t511 / 0.2e1;
t514 = Ifges(7,4) + Ifges(6,5);
t513 = -Ifges(6,6) + Ifges(7,6);
t415 = qJ(6) * t511;
t455 = pkin(5) * t319;
t110 = t455 + t415;
t525 = mrSges(6,3) + mrSges(7,2);
t387 = mrSges(6,3) / 0.2e1 + mrSges(7,2) / 0.2e1;
t524 = t387 * t511;
t515 = mrSges(7,1) + mrSges(6,1);
t523 = t515 * t517;
t289 = t304 * t308 - t306 * t458;
t437 = t511 * mrSges(6,3);
t150 = -mrSges(6,2) * t289 - t437;
t279 = t289 * mrSges(7,3);
t155 = -mrSges(7,2) * t511 + t279;
t395 = t150 + t155;
t297 = -pkin(4) * t305 - pkin(3);
t229 = -pkin(5) * t504 - qJ(6) * t291 + t297;
t247 = -mrSges(7,1) * t504 - mrSges(7,3) * t291;
t522 = m(7) * t229 + t247;
t283 = Ifges(6,4) * t504;
t442 = Ifges(7,5) * t504;
t521 = t283 - t442 + (Ifges(6,1) + Ifges(7,1)) * t291;
t416 = qJ(4) * t289;
t456 = pkin(3) * t292;
t246 = t416 + t456;
t452 = pkin(7) + qJ(2);
t294 = t452 * t306;
t362 = t452 * t304;
t255 = t294 * t308 + t458 * t362;
t142 = t305 * t246 + t255 * t303;
t143 = t303 * t246 - t305 * t255;
t520 = -t142 * t303 + t143 * t305;
t451 = pkin(8) + qJ(4);
t293 = t451 * t305;
t361 = t451 * t303;
t254 = t293 * t307 + t361 * t457;
t256 = t293 * t457 - t307 * t361;
t519 = t254 * t526 + t256 * t517;
t404 = t504 * qJ(6);
t454 = pkin(5) * t291;
t245 = -t404 + t454;
t518 = -t245 / 0.2e1;
t475 = t319 / 0.4e1;
t381 = -pkin(2) * t306 - pkin(1);
t230 = pkin(3) * t289 - qJ(4) * t292 + t381;
t257 = t294 * t458 - t308 * t362;
t135 = t305 * t230 - t257 * t303;
t399 = t292 * t305;
t80 = pkin(4) * t289 - pkin(8) * t399 + t135;
t136 = t303 * t230 + t305 * t257;
t400 = t292 * t303;
t98 = -pkin(8) * t400 + t136;
t51 = -t307 * t98 + t457 * t80;
t44 = -t289 * pkin(5) - t51;
t448 = t44 + t51;
t430 = t319 * mrSges(6,3);
t210 = Ifges(7,5) * t319;
t89 = t289 * Ifges(7,6) + Ifges(7,3) * t511 + t210;
t512 = -Ifges(7,1) * t511 + t210 + t89;
t510 = t319 * t387;
t153 = mrSges(6,1) * t289 - t430;
t154 = -mrSges(7,1) * t289 + mrSges(7,2) * t319;
t394 = -t153 + t154;
t282 = Ifges(7,5) * t291;
t248 = -Ifges(7,3) * t504 + t282;
t509 = Ifges(7,1) * t504 + t248 + t282;
t299 = t303 ^ 2;
t301 = t305 ^ 2;
t391 = t301 + t299;
t503 = t319 * t513 - t511 * t514;
t213 = Ifges(6,4) * t511;
t443 = Ifges(7,5) * t511;
t93 = Ifges(7,1) * t319 + t289 * Ifges(7,4) + t443;
t95 = Ifges(6,1) * t319 + t289 * Ifges(6,5) - t213;
t502 = -Ifges(6,2) * t319 - t213 + t93 + t95;
t501 = -Ifges(6,2) * t291 + t283 + t521;
t500 = t154 / 0.2e1 - t153 / 0.2e1;
t298 = m(7) * qJ(6) + mrSges(7,3);
t499 = -m(7) * pkin(5) - t515;
t498 = -mrSges(6,2) + t298;
t218 = t289 * t504;
t433 = t218 * mrSges(7,3);
t435 = t218 * mrSges(6,2);
t215 = t291 * t289;
t439 = t215 * mrSges(7,1);
t440 = t215 * mrSges(6,1);
t497 = -t433 / 0.2e1 + t435 / 0.2e1 + t439 / 0.2e1 + t440 / 0.2e1;
t496 = 0.2e1 * m(7);
t495 = 0.2e1 * t292;
t494 = m(5) / 0.2e1;
t493 = -m(6) / 0.2e1;
t492 = m(6) / 0.2e1;
t491 = -m(7) / 0.2e1;
t490 = m(7) / 0.2e1;
t489 = -mrSges(6,2) / 0.2e1;
t402 = t289 * t303;
t109 = pkin(8) * t402 + t143;
t401 = t289 * t305;
t87 = pkin(4) * t292 + pkin(8) * t401 + t142;
t56 = -t307 * t109 + t457 * t87;
t50 = -t292 * pkin(5) - t56;
t488 = -t50 / 0.2e1;
t384 = t457 * t98;
t419 = t307 * t80;
t52 = t384 + t419;
t487 = t52 / 0.2e1;
t156 = t215 * mrSges(7,2) + mrSges(7,3) * t292;
t486 = t156 / 0.2e1;
t484 = t215 / 0.2e1;
t483 = -t215 / 0.2e1;
t481 = t511 / 0.4e1;
t477 = -t319 / 0.4e1;
t474 = -t247 / 0.2e1;
t472 = t504 / 0.2e1;
t471 = -t504 / 0.2e1;
t470 = -t504 / 0.4e1;
t469 = -t289 / 0.2e1;
t466 = -t291 / 0.2e1;
t465 = -t291 / 0.4e1;
t464 = t291 / 0.2e1;
t463 = t291 / 0.4e1;
t461 = t292 / 0.2e1;
t460 = t303 / 0.2e1;
t459 = t305 / 0.2e1;
t403 = t289 * qJ(6);
t43 = t52 + t403;
t449 = -t43 + t52;
t447 = Ifges(5,4) * t303;
t446 = Ifges(5,4) * t305;
t445 = Ifges(6,4) * t319;
t444 = Ifges(6,4) * t291;
t111 = t433 - t439;
t112 = -t435 - t440;
t113 = mrSges(7,1) * t511 - mrSges(7,3) * t319;
t114 = mrSges(6,1) * t511 + mrSges(6,2) * t319;
t425 = t292 * mrSges(6,2);
t149 = t215 * mrSges(6,3) - t425;
t427 = t292 * mrSges(6,1);
t151 = t218 * mrSges(6,3) + t427;
t426 = t292 * mrSges(7,1);
t434 = t218 * mrSges(7,2);
t152 = -t426 - t434;
t423 = t303 * Ifges(5,2);
t172 = Ifges(5,6) * t292 + (t423 - t446) * t289;
t173 = Ifges(5,5) * t292 + (-t305 * Ifges(5,1) + t447) * t289;
t179 = pkin(4) * t400 + t255;
t180 = -pkin(4) * t402 + t257;
t421 = t305 * mrSges(5,2);
t424 = t303 * mrSges(5,1);
t344 = t421 + t424;
t227 = t344 * t289;
t228 = t344 * t292;
t237 = -t292 * mrSges(5,2) + mrSges(5,3) * t402;
t238 = -t289 * mrSges(5,2) - mrSges(5,3) * t400;
t239 = t292 * mrSges(5,1) + mrSges(5,3) * t401;
t240 = t289 * mrSges(5,1) - mrSges(5,3) * t399;
t281 = t292 * mrSges(4,1);
t385 = -Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t386 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t318 = -t215 * t385 - t218 * t386;
t420 = t305 * Ifges(5,5);
t422 = t303 * Ifges(5,6);
t57 = t457 * t109 + t307 * t87;
t49 = qJ(6) * t292 + t57;
t335 = pkin(5) * t511 - qJ(6) * t319;
t66 = t179 + t335;
t336 = pkin(5) * t215 - t218 * qJ(6);
t67 = t180 - t336;
t88 = -Ifges(7,5) * t218 + Ifges(7,6) * t292 - Ifges(7,3) * t215;
t90 = -Ifges(6,4) * t218 + Ifges(6,2) * t215 + Ifges(6,6) * t292;
t91 = -Ifges(6,2) * t511 + t289 * Ifges(6,6) + t445;
t92 = -Ifges(7,1) * t218 + Ifges(7,4) * t292 - Ifges(7,5) * t215;
t94 = -Ifges(6,1) * t218 + Ifges(6,4) * t215 + Ifges(6,5) * t292;
t1 = (t94 / 0.2e1 + t92 / 0.2e1) * t319 - (t93 / 0.2e1 + t95 / 0.2e1) * t218 + (-t303 * t172 / 0.2e1 + t173 * t459 + (t420 / 0.2e1 - t422 / 0.2e1 - Ifges(4,4)) * t292 + t386 * t319 + t385 * t511) * t292 + (-t90 / 0.2e1 + t88 / 0.2e1) * t511 + m(5) * (t135 * t142 + t136 * t143 + t255 * t257) + m(6) * (t179 * t180 + t51 * t56 + t52 * t57) + m(7) * (t43 * t49 + t44 * t50 + t66 * t67) - (-t91 / 0.2e1 + t89 / 0.2e1) * t215 + (-t381 * mrSges(4,2) + (-t301 * Ifges(5,1) / 0.2e1 + Ifges(5,3) - Ifges(4,1) + Ifges(6,3) + Ifges(7,2) + Ifges(4,2) + (t446 - t423 / 0.2e1) * t303) * t292 + t318 + (Ifges(4,4) - t420 + t422) * t289) * t289 + t257 * t228 - t255 * t227 + t136 * t237 + t143 * t238 + t135 * t239 + t142 * t240 + t179 * t112 + t180 * t114 + t52 * t149 + t57 * t150 + t51 * t151 + t44 * t152 + t56 * t153 + t50 * t154 + t49 * t155 + t43 * t156 + t381 * t281 + t66 * t111 + t67 * t113;
t441 = t1 * qJD(1);
t438 = t511 * mrSges(6,2);
t436 = t511 * mrSges(7,3);
t432 = t319 * mrSges(6,1);
t431 = t319 * mrSges(7,1);
t429 = t504 * mrSges(7,2);
t428 = t291 * mrSges(7,2);
t338 = Ifges(7,3) * t319 - t443;
t340 = -Ifges(6,1) * t511 - t445;
t341 = t431 + t436;
t343 = t432 - t438;
t4 = t179 * t343 + t66 * t341 + t91 * t517 + t511 * t338 / 0.2e1 + (-t319 * t43 - t44 * t511) * mrSges(7,2) + (m(7) * t66 + t113) * t110 + (m(7) * t44 + t394 - t430) * t52 + (m(7) * t43 + t395 + t437) * t51 + t503 * t289 / 0.2e1 + t502 * t526 + (t340 + t512) * t319 / 0.2e1;
t418 = t4 * qJD(1);
t331 = t135 * t303 - t136 * t305;
t397 = t305 * t238;
t398 = t303 * t240;
t406 = t255 * t292;
t9 = -t395 * t218 - t394 * t215 + (mrSges(4,3) * t289 - t397 + t398) * t289 + (mrSges(4,3) * t292 + t113 + t114 + t228) * t292 + m(6) * (t179 * t292 + t215 * t51 - t218 * t52) + m(7) * (-t215 * t44 - t218 * t43 + t292 * t66) + m(5) * (t289 * t331 + t406) + m(4) * (-t257 * t289 + t406) + (m(3) * qJ(2) + mrSges(3,3)) * (t304 ^ 2 + t306 ^ 2);
t417 = t9 * qJD(1);
t328 = -t215 * t254 - t218 * t256;
t345 = -t305 * mrSges(5,1) + t303 * mrSges(5,2);
t310 = -m(5) * (-t391 * t416 - t456) / 0.2e1 + (t292 * t297 + t328) * t493 + (t229 * t292 + t328) * t491 - (t247 + t345) * t292 / 0.2e1;
t313 = (t142 * t305 + t143 * t303) * t494 + (t291 * t57 + t504 * t56) * t492 + (t291 * t49 - t50 * t504) * t490 + t237 * t460 + t239 * t459;
t367 = -t151 / 0.2e1 + t152 / 0.2e1;
t369 = t149 / 0.2e1 + t486;
t10 = t281 + (-mrSges(4,2) + (t301 / 0.2e1 + t299 / 0.2e1) * mrSges(5,3)) * t289 + (-t425 / 0.2e1 + t387 * t215 + t369) * t291 - (-t427 / 0.2e1 - t387 * t218 + t367) * t504 + t310 + t313;
t414 = qJD(1) * t10;
t12 = -t395 * t511 + t394 * t319 + m(6) * (-t319 * t51 - t511 * t52) + m(7) * (t319 * t44 - t43 * t511) + (m(5) * (-t135 * t305 - t136 * t303) - t303 * t238 - t305 * t240) * t292;
t413 = qJD(1) * t12;
t21 = t289 * t155 + m(7) * (t289 * t43 - t319 * t66) - t319 * t113;
t412 = qJD(1) * t21;
t409 = t511 * t504;
t408 = t319 * t291;
t389 = qJD(5) * t291;
t370 = t289 * t464;
t117 = (t370 + t484) * m(7);
t388 = t117 * qJD(1);
t371 = t504 * t469;
t368 = t150 / 0.2e1 + t155 / 0.2e1;
t363 = m(5) * t391;
t360 = t291 * mrSges(7,1) - mrSges(7,3) * t504;
t359 = t291 * mrSges(6,1) + mrSges(6,2) * t504;
t349 = m(6) / 0.4e1 + m(7) / 0.4e1 + m(5) / 0.4e1;
t348 = t153 * t466 + t154 * t464 + t395 * t472;
t342 = -mrSges(6,1) * t504 + mrSges(6,2) * t291;
t339 = Ifges(6,1) * t504 - t444;
t337 = Ifges(7,3) * t291 + t442;
t315 = (-t229 * t319 + t256 * t289 - t66 * t291) * t490 + t319 * t474 + t113 * t466;
t325 = m(7) * t488 + t426 / 0.2e1;
t19 = (-t371 + t218 / 0.2e1) * mrSges(7,2) + t315 + t325;
t96 = t522 * t291;
t334 = -qJD(1) * t19 + qJD(3) * t96;
t25 = (-t455 / 0.4e1 - t415 / 0.4e1 - t110 / 0.4e1) * t496 + 0.2e1 * t523 + 0.2e1 * (mrSges(7,3) - mrSges(6,2)) * t526;
t323 = t359 + t360;
t65 = (t245 / 0.4e1 + t454 / 0.4e1 - t404 / 0.4e1) * t496 + t323;
t333 = qJD(1) * t25 - qJD(3) * t65;
t316 = (t493 + t491) * (-t291 * t511 - t319 * t504);
t32 = (t363 / 0.4e1 + t349) * t495 + t316;
t332 = qJD(1) * t32;
t330 = t254 * t319 - t256 * t511;
t231 = m(7) * t291;
t99 = 0.2e1 * t475 * t496;
t327 = qJD(1) * t99 + qJD(3) * t231;
t23 = -t279 + (t52 / 0.4e1 - t384 / 0.4e1 - t403 / 0.2e1 - t419 / 0.4e1) * t496;
t326 = qJD(1) * t23 - qJD(5) * t298;
t249 = Ifges(6,2) * t504 + t444;
t14 = t337 * t471 + t297 * t359 + t229 * t360 + t249 * t466 + t501 * t472 + (t339 + t509) * t464 + t522 * t245;
t309 = (t110 * t229 + t245 * t66 + t254 * t449 + t256 * t448) * t491 + t110 * t474 - t179 * t359 / 0.2e1 + t249 * t475 - t229 * t341 / 0.2e1 + t113 * t518 + t91 * t463 - t297 * t343 / 0.2e1 - t66 * t360 / 0.2e1 + t509 * t477 - (t291 * t513 + t504 * t514) * t289 / 0.4e1 + t512 * t465 + t501 * t481 + t502 * t470;
t311 = (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t292 + (-pkin(5) * t50 + qJ(6) * t49) * t490 - pkin(5) * t152 / 0.2e1 + qJ(6) * t486 + t49 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t488 + t56 * mrSges(6,1) / 0.2e1 + t57 * t489 + t318;
t2 = t309 + t311 + (t430 / 0.2e1 - t500) * t256 + (t437 / 0.2e1 + t368) * t254 - Ifges(7,5) * t409 / 0.2e1 + Ifges(6,4) * t408 / 0.2e1 + ((-t52 / 0.2e1 + t43 / 0.2e1) * t291 - (t44 / 0.2e1 + t51 / 0.2e1) * t504 - t519) * mrSges(7,2) + (-Ifges(6,1) + Ifges(7,3)) * (t465 * t511 - t477 * t504);
t322 = -t2 * qJD(1) + t14 * qJD(3);
t317 = t336 * t490 + t497;
t5 = (t448 * t491 - t500 + t510) * t291 - (t449 * t491 + t368 + t524) * t504 + t317;
t321 = t5 * qJD(1);
t34 = (t291 ^ 2 + t504 ^ 2) * t525 + (m(5) * qJ(4) + mrSges(5,3)) * t391 + (m(7) + m(6)) * (t254 * t291 + t256 * t504);
t312 = t331 * t494 + (-t291 * t51 + t504 * t52 + t330) * t493 + (t291 * t44 + t43 * t504 + t330) * t491 + t398 / 0.2e1 - t397 / 0.2e1;
t314 = (-t421 / 0.2e1 - t424 / 0.2e1) * t289 + t257 * t494 + t180 * t492 + t67 * t490 - t497;
t7 = (-t500 - t510) * t291 - (t368 - t524) * t504 + t312 + t314;
t320 = -qJD(1) * t7 + qJD(3) * t34;
t144 = m(7) * t256 + t429;
t116 = (t370 + t483) * m(7);
t100 = m(7) * t517 + t319 * t490;
t33 = t349 * t495 - t363 * t461 - t316;
t26 = -t438 / 0.2e1 + t432 / 0.2e1 + t436 / 0.2e1 + t431 / 0.2e1 + mrSges(7,3) * t526 - t511 * t489 + t523;
t24 = (t52 + 0.2e1 * t403) * t490 + m(7) * t487 + t155;
t18 = -mrSges(7,2) * t371 - t434 / 0.2e1 + t315 - t325;
t11 = t391 * mrSges(5,3) * t469 - t367 * t504 + t369 * t291 + t342 * t461 - t310 + t313 + t525 * (-t215 * t464 - t218 * t472);
t8 = -t312 + t314 + t348 + t525 * (t319 * t464 - t472 * t511);
t6 = (t291 * t448 - t449 * t504) * t490 + t317 + t348 + t525 * (t409 / 0.2e1 - t408 / 0.2e1);
t3 = t337 * t481 + t338 * t470 + t339 * t475 + t340 * t463 - t309 + t311 - t395 * t254 / 0.2e1 + (t487 - t43 / 0.2e1) * t428 + t448 * t429 / 0.2e1 + t500 * t256 + t525 * t519;
t13 = [qJD(2) * t9 + qJD(3) * t1 + qJD(4) * t12 + qJD(5) * t4 + qJD(6) * t21, t11 * qJD(3) + t33 * qJD(4) + t6 * qJD(5) + t116 * qJD(6) + t417 + 0.2e1 * (t490 + t492) * qJD(2) * (t215 * t504 - t291 * t218) t11 * qJD(2) + t8 * qJD(4) + t3 * qJD(5) + t18 * qJD(6) + t441 + (t180 * t342 + 0.2e1 * (-pkin(3) * t257 + t520 * qJ(4)) * t494 + t520 * mrSges(5,3) - t521 * t218 / 0.2e1 + 0.2e1 * (t180 * t297 - t254 * t56 + t256 * t57) * t492 + t248 * t483 + t249 * t484 + 0.2e1 * (t229 * t67 + t254 * t50 + t256 * t49) * t490 + t88 * t471 + t90 * t472 + t172 * t459 + t173 * t460 + t305 * qJ(4) * t237 - t303 * qJ(4) * t239 + t257 * t345 + t50 * t428 + (-t305 * (Ifges(5,1) * t303 + t446) / 0.2e1 + (Ifges(5,2) * t305 + t447) * t460 - Ifges(4,5)) * t289 + t297 * t112 - Ifges(4,6) * t292 + t49 * t429 + t255 * mrSges(4,2) + t256 * t149 + t256 * t156 - t257 * mrSges(4,1) + t254 * t152 - t254 * t151 + t67 * t247 + pkin(3) * t227 + t229 * t111 + (t94 + t92) * t464 + (Ifges(5,5) * t303 + Ifges(5,6) * t305 + t291 * t514 - t504 * t513) * t461 + (-t291 * t56 + t504 * t57) * mrSges(6,3)) * qJD(3), qJD(2) * t33 + qJD(3) * t8 + qJD(5) * t26 + qJD(6) * t100 + t413, t418 + t6 * qJD(2) + t3 * qJD(3) + t26 * qJD(4) + (t335 * mrSges(7,2) + t498 * t51 + t499 * t52 + t503) * qJD(5) + t24 * qJD(6), qJD(2) * t116 + qJD(3) * t18 + qJD(4) * t100 + qJD(5) * t24 + t412; qJD(3) * t10 - qJD(4) * t32 - qJD(5) * t5 + qJD(6) * t117 - t417, 0, t414, -t332, -t323 * qJD(5) + (qJD(5) * t518 + qJD(6) * t464) * t496 - t321, m(7) * t389 + t388; -qJD(2) * t10 - qJD(4) * t7 - qJD(5) * t2 + qJD(6) * t19 - t441, -t414, qJD(4) * t34 + qJD(5) * t14 - qJD(6) * t96, t320, t144 * qJD(6) + (-mrSges(7,2) * qJ(6) + t513) * t389 + t322 + (-t498 * t254 + t499 * t256 - (mrSges(7,2) * pkin(5) - t514) * t504) * qJD(5), qJD(5) * t144 - t334; qJD(2) * t32 + qJD(3) * t7 - qJD(5) * t25 - qJD(6) * t99 - t413, t332, qJD(5) * t65 - qJD(6) * t231 - t320, 0, -t333, -t327; qJD(2) * t5 + qJD(3) * t2 + qJD(4) * t25 - qJD(6) * t23 - t418, t321, -qJD(4) * t65 - t322, t333, t298 * qJD(6), -t326; -qJD(2) * t117 - qJD(3) * t19 + qJD(4) * t99 + qJD(5) * t23 - t412, -t388, qJD(4) * t231 + t334, t327, t326, 0;];
Cq  = t13;
