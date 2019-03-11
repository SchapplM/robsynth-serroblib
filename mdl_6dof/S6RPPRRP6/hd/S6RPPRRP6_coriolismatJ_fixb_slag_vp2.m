% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:13
% EndTime: 2019-03-09 02:10:21
% DurationCPUTime: 4.84s
% Computational Cost: add. (4480->416), mult. (8927->545), div. (0->0), fcn. (6558->4), ass. (0->227)
t270 = sin(qJ(5));
t272 = cos(qJ(5));
t422 = Ifges(6,6) - Ifges(7,6);
t423 = Ifges(6,5) + Ifges(7,4);
t475 = -t422 * t270 + t423 * t272;
t474 = -Ifges(5,4) + t475;
t262 = t270 ^ 2;
t264 = t272 ^ 2;
t369 = t262 + t264;
t395 = qJ(6) * t272;
t428 = pkin(5) * t270;
t214 = -t395 + t428;
t266 = qJ(2) - pkin(7);
t273 = cos(qJ(4));
t117 = (t214 - t266) * t273;
t255 = t270 * mrSges(7,1);
t402 = t272 * mrSges(7,3);
t464 = t255 - t402;
t473 = t117 * t464;
t417 = Ifges(6,4) * t270;
t315 = Ifges(6,1) * t272 - t417;
t259 = Ifges(7,5) * t270;
t462 = Ifges(7,1) * t272 + t259;
t471 = t315 + t462;
t470 = t369 * pkin(8);
t469 = mrSges(7,1) + mrSges(6,1);
t468 = Ifges(7,2) + Ifges(6,3);
t271 = sin(qJ(4));
t333 = t266 * t270 - pkin(5);
t267 = pkin(1) + qJ(3);
t426 = pkin(8) * t273;
t429 = pkin(4) * t271;
t205 = t267 - t426 + t429;
t388 = t205 * t272;
t75 = t333 * t271 - t388;
t382 = t270 * t271;
t96 = -t266 * t382 + t388;
t420 = t75 + t96;
t381 = t270 * t273;
t359 = mrSges(6,3) * t381;
t199 = -t271 * mrSges(6,2) - t359;
t358 = mrSges(7,2) * t381;
t405 = t271 * mrSges(7,3);
t203 = -t358 + t405;
t467 = t199 + t203;
t375 = t272 * t273;
t356 = mrSges(6,3) * t375;
t201 = t271 * mrSges(6,1) - t356;
t357 = mrSges(7,2) * t375;
t202 = -t271 * mrSges(7,1) + t357;
t466 = t201 - t202;
t309 = pkin(5) * t272 + qJ(6) * t270;
t210 = -pkin(4) - t309;
t403 = t272 * mrSges(7,1);
t410 = t270 * mrSges(7,3);
t317 = t403 + t410;
t465 = m(7) * t210 - t317;
t256 = t270 * mrSges(6,1);
t258 = t272 * mrSges(6,2);
t463 = t258 + t256;
t260 = Ifges(6,4) * t272;
t461 = -Ifges(6,2) * t270 + t260;
t294 = t315 * t273;
t460 = t423 * t271 + t273 * t462 + t294;
t458 = -t369 * t426 / 0.2e1;
t431 = m(7) * t214;
t299 = t464 + t431;
t216 = t272 * Ifges(6,2) + t417;
t414 = Ifges(7,5) * t272;
t217 = t270 * Ifges(7,1) - t414;
t436 = t272 / 0.2e1;
t441 = -t270 / 0.2e1;
t457 = t216 * t441 + t217 * t436;
t261 = pkin(4) * t273;
t427 = pkin(8) * t271;
t223 = t261 + t427;
t386 = t223 * t272;
t114 = -t266 * t381 + t386;
t115 = t270 * t223 + t266 * t375;
t302 = -t114 * t270 + t115 * t272;
t78 = qJ(6) * t273 + t115;
t93 = t333 * t273 - t386;
t306 = t270 * t93 + t272 * t78;
t310 = t272 * Ifges(7,3) - t259;
t314 = -t270 * Ifges(6,1) - t260;
t455 = -t255 / 0.2e1 - t256 / 0.2e1;
t254 = m(7) * qJ(6) + mrSges(7,3);
t454 = t369 * mrSges(6,3);
t383 = t270 * t205;
t385 = t266 * t272;
t74 = t383 + (qJ(6) + t385) * t271;
t380 = t271 * t272;
t97 = t266 * t380 + t383;
t421 = -t74 + t97;
t437 = -t272 / 0.2e1;
t449 = m(7) / 0.2e1;
t452 = (t421 * t270 + t420 * t272) * t449 + t202 * t436 + t201 * t437;
t263 = t271 ^ 2;
t265 = t273 ^ 2;
t451 = m(6) / 0.2e1;
t450 = -m(7) / 0.2e1;
t448 = -pkin(5) / 0.2e1;
t447 = -mrSges(7,1) / 0.2e1;
t446 = -mrSges(6,2) / 0.2e1;
t445 = mrSges(7,3) / 0.2e1;
t444 = Ifges(6,6) / 0.2e1;
t443 = qJ(6) / 0.2e1;
t442 = t214 / 0.2e1;
t440 = t270 / 0.2e1;
t434 = -t273 / 0.2e1;
t433 = t273 / 0.2e1;
t430 = m(7) * t272;
t419 = m(7) * qJD(6);
t182 = t464 * t273;
t296 = t463 * t273;
t305 = -t270 * t96 + t272 * t97;
t307 = t270 * t75 + t272 * t74;
t368 = t265 + t263;
t334 = m(5) * t368;
t377 = t272 * t203;
t378 = t272 * t199;
t9 = t368 * mrSges(5,3) - mrSges(4,2) - mrSges(3,3) + t466 * t382 + (m(7) * t117 + t182 + t296) * t273 + (-m(6) * t265 - t334) * t266 + (-m(6) * t305 - m(7) * t307 - t377 - t378) * t271 + (-m(4) - m(3)) * qJ(2);
t412 = qJD(1) * t9;
t411 = t270 * mrSges(6,2);
t408 = t270 * t74;
t406 = t270 * t97;
t404 = t272 * mrSges(6,1);
t399 = t273 * mrSges(7,1);
t311 = Ifges(7,3) * t270 + t414;
t293 = t311 * t273;
t285 = Ifges(7,6) * t271 + t293;
t286 = Ifges(6,6) * t271 + t273 * t461;
t292 = t309 * t273;
t295 = t317 * t273;
t340 = -t375 / 0.2e1;
t379 = t271 * t273;
t344 = -t379 / 0.2e1;
t363 = -t265 / 0.2e1;
t319 = t404 - t411;
t180 = t319 * t273;
t374 = t273 * t180;
t4 = t266 * t374 + t286 * t375 / 0.2e1 + t285 * t340 + t75 * t358 + t74 * t357 - t182 * t292 + (t270 * t310 + t272 * t314) * t363 + t460 * t381 / 0.2e1 + (-t423 * t270 - t422 * t272) * t344 + t457 * t265 + (-m(7) * t292 - t295) * t117 + (-m(7) * t75 + t356 + t466) * t97 + (-m(7) * t74 - t359 - t467) * t96;
t398 = t4 * qJD(1);
t348 = -t382 / 0.2e1;
t276 = -t265 * t309 * t449 - t374 / 0.2e1 + t317 * t363 + t467 * t348 + t452 * t271 + (mrSges(7,2) + mrSges(6,3)) * t369 * t344;
t279 = t309 * t450 + t411 / 0.2e1 - t410 / 0.2e1 - t404 / 0.2e1 - t403 / 0.2e1;
t8 = t276 + t279;
t397 = t8 * qJD(1);
t396 = -t319 - mrSges(5,1);
t14 = t271 * mrSges(5,1) + t273 * mrSges(5,2) + mrSges(4,3) + t466 * t272 + t467 * t270 + m(7) * (-t272 * t75 + t408) + m(6) * (t272 * t96 + t406) + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t267;
t394 = qJD(1) * t14;
t27 = m(7) * (-t117 * t375 + t271 * t74) + t271 * t203 - t182 * t375;
t393 = qJD(1) * t27;
t200 = -mrSges(7,2) * t380 - t399;
t204 = mrSges(7,2) * t382 + t273 * mrSges(7,3);
t297 = -mrSges(6,2) * t273 + mrSges(6,3) * t382;
t298 = mrSges(6,1) * t273 + mrSges(6,3) * t380;
t277 = (-t114 * t272 - t115 * t270) * t451 + (-t270 * t78 + t272 * t93) * t449 + t200 * t436 + t298 * t437 + (t204 + t297) * t441;
t372 = t470 * t271;
t387 = t210 * t273;
t283 = (t261 + t372) * t451 + (t372 - t387) * t449;
t335 = t262 / 0.2e1 + t264 / 0.2e1;
t321 = t335 * mrSges(7,2);
t332 = t273 * mrSges(5,1) - t271 * mrSges(5,2);
t338 = t317 * t434;
t10 = -t319 * t434 - t277 + t283 + t332 - t338 + (t321 + t454 / 0.2e1) * t271;
t392 = t10 * qJD(1);
t241 = qJ(6) * t380;
t322 = (t445 + t446) * t271;
t336 = t203 / 0.2e1 + t199 / 0.2e1;
t337 = t202 / 0.2e1 - t201 / 0.2e1;
t353 = -mrSges(6,1) / 0.2e1 + t447;
t12 = t241 * t449 + (t420 * t449 + (m(7) * t448 + t353) * t271 + t337) * t270 + (t421 * t450 + t322 + t336) * t272;
t389 = t12 * qJD(1);
t384 = t270 * t200;
t376 = t272 * t204;
t284 = -m(5) / 0.2e1 + (-m(6) / 0.2e1 + t450) * t369;
t362 = m(7) / 0.4e1 + m(6) / 0.4e1;
t289 = t334 / 0.2e1 + 0.2e1 * t362 * (t369 * t263 + t265);
t31 = m(4) - t284 + t289;
t373 = t31 * qJD(1);
t371 = t470 * t273;
t367 = qJD(4) * t271;
t366 = qJD(4) * t273;
t365 = qJD(5) * t272;
t108 = (t265 / 0.2e1 + t263 / 0.2e1 + 0.1e1 / 0.2e1) * t430;
t364 = t108 * qJD(1);
t361 = m(7) * t382;
t360 = t270 * t419;
t351 = -t258 / 0.2e1;
t326 = pkin(5) * mrSges(7,2) - t423;
t323 = -qJ(6) * mrSges(7,2) - t422;
t320 = t463 * t434;
t116 = t241 + (t266 - t428) * t271;
t181 = t463 * t271;
t3 = -t267 * t332 - t114 * t201 - t115 * t199 - t116 * t182 - t97 * t297 - t96 * t298 - t93 * t202 - t75 * t200 - t78 * t203 - t74 * t204 - m(7) * (t116 * t117 + t74 * t78 + t75 * t93) - m(6) * (t114 * t96 + t115 * t97) + t474 * t263 + (-t266 * t181 - t474 * t273) * t273 + (-t266 * t296 + t473 + (Ifges(5,1) - Ifges(5,2) + m(6) * t266 ^ 2 + (Ifges(6,1) + Ifges(7,1)) * t264 + ((Ifges(6,2) + Ifges(7,3)) * t270 + 0.2e1 * (-Ifges(6,4) + Ifges(7,5)) * t272) * t270 - t468) * t273) * t271;
t6 = (t181 / 0.2e1 + t336 * t272 + t337 * t270 + (-t116 + t307) * t449 + t305 * t451) * t273 + (t384 / 0.2e1 + t376 / 0.2e1 + t182 / 0.2e1 + t464 * t433 + (t117 + t306) * t449 + (-0.2e1 * t266 * t273 + t302) * t451) * t271;
t308 = -t3 * qJD(1) + t6 * qJD(3);
t41 = 0.4e1 * t362 * (-0.1e1 + t369) * t379;
t304 = t6 * qJD(1) + t41 * qJD(3);
t281 = (-t117 * t270 + (-t387 + t427) * t272) * t449 + t182 * t441;
t300 = t93 * t449 - t399 / 0.2e1;
t21 = (t271 * mrSges(7,2) - t338) * t272 + t281 - t300;
t76 = t465 * t270;
t303 = -qJD(1) * t21 + qJD(4) * t76;
t35 = t405 + 0.2e1 * (t383 / 0.4e1 - t97 / 0.4e1 + (t385 / 0.4e1 + t443) * t271) * m(7);
t301 = qJD(1) * t35 + qJD(5) * t254;
t275 = t216 * t340 + t266 * t320 + t182 * t442 + (t214 * t117 + t210 * t292) * t449 + t473 / 0.2e1 + t210 * t295 / 0.2e1 - t317 * t292 / 0.2e1 - pkin(4) * t180 / 0.2e1 - t270 * t286 / 0.4e1 + (t293 + t285) * t270 / 0.4e1 + t475 * t271 / 0.4e1 + t458 * mrSges(6,3) + (t294 + t460) * t272 / 0.4e1 + (t467 * t441 + t452) * pkin(8) + (t420 * t436 + t406 / 0.2e1 - t408 / 0.2e1 + t458) * mrSges(7,2) + (t314 / 0.2e1 - t217 / 0.2e1 - t461 / 0.4e1) * t381 + (-t310 / 0.2e1 + t462 / 0.4e1) * t375;
t278 = (-pkin(5) * t93 + qJ(6) * t78) * t449 + t200 * t448 + t204 * t443 + t114 * mrSges(6,1) / 0.2e1 + t115 * t446 + t78 * t445 + t93 * t447;
t1 = Ifges(7,6) * t348 + t382 * t444 - t275 + t278 + t468 * t433 - t423 * t380 / 0.2e1;
t282 = -t310 * t440 - t314 * t436 + t457;
t17 = -pkin(4) * t463 + t299 * t210 - t214 * t317 + t311 * t437 + t461 * t436 + t471 * t440 + t282;
t29 = (t351 + t258 / 0.2e1 + t353 * t270 + (-t428 / 0.2e1 + t395 / 0.2e1 + t442) * m(7) - t455) * t273;
t291 = t1 * qJD(1) + t29 * qJD(3) - t17 * qJD(4);
t280 = (-m(7) * t309 - t317 - t319) * qJD(5);
t211 = (m(7) * pkin(8) + mrSges(7,2)) * t272;
t109 = (-0.1e1 / 0.2e1 + t368 / 0.2e1) * t430;
t33 = t203 + (t383 + (0.2e1 * qJ(6) + t385) * t271 + t97) * t449;
t32 = t284 + t289;
t30 = t320 + (t351 + t402 / 0.2e1 - t431 / 0.2e1 + t455) * t273 + t299 * t434;
t24 = -t272 * t338 + t281 + t300;
t13 = t202 * t441 - t377 / 0.2e1 - t378 / 0.2e1 + t201 * t440 + t272 * t322 + (-pkin(5) * t382 - t420 * t270 + t421 * t272 + t241) * t449 + t469 * t348;
t11 = (t319 / 0.2e1 + t317 / 0.2e1) * t273 + (t335 * mrSges(6,3) + t321) * t271 + t277 + t283;
t7 = t276 - t279;
t5 = t6 * qJD(4);
t2 = ((-Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1) * t272 + (t444 - Ifges(7,6) / 0.2e1) * t270) * t271 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t273 + t275 + t278;
t15 = [-qJD(2) * t9 + qJD(3) * t14 - qJD(4) * t3 - qJD(5) * t4 + qJD(6) * t27, qJD(3) * t32 + qJD(4) * t11 + qJD(5) * t13 - t412, qJD(2) * t32 + qJD(5) * t7 + qJD(6) * t109 + t394 + t5, t11 * qJD(2) + t2 * qJD(5) + t24 * qJD(6) + (t311 * t436 + t461 * t437 - t210 * t464 - Ifges(5,5) + (-m(6) * pkin(4) + t396) * t266 - t282 + t471 * t441) * t367 + (-mrSges(5,2) * t266 - Ifges(5,6) + (-pkin(8) * mrSges(6,2) + t422) * t272 + (-pkin(8) * mrSges(6,1) + t423) * t270) * t366 + t308 + (pkin(4) * t181 + (m(6) * t302 + m(7) * t306 + t376 + t384) * pkin(8) + t465 * t116 + t302 * mrSges(6,3) + t306 * mrSges(7,2)) * qJD(4), t13 * qJD(2) + t7 * qJD(3) + t2 * qJD(4) + t33 * qJD(6) - t398 + ((-m(7) * pkin(5) - t469) * t97 + (-mrSges(6,2) + t254) * t96 + (t270 * t326 + t272 * t323) * t273) * qJD(5), qJD(3) * t109 + qJD(4) * t24 + qJD(5) * t33 + t393; -qJD(3) * t31 - qJD(4) * t10 - qJD(5) * t12 - t271 * t360 + t412, 0, -t373, -t392, -t389 + (t299 + t463) * qJD(5) - t360 (-qJD(1) * t271 - qJD(5)) * t270 * m(7); qJD(2) * t31 + qJD(5) * t8 + qJD(6) * t108 - t394 + t5, t373, t41 * qJD(4), t30 * qJD(5) + (-t317 + t396) * t367 + 0.2e1 * ((t210 * t271 + t371) * t449 + (t371 - t429) * t451) * qJD(4) + ((t369 * mrSges(7,2) - mrSges(5,2) + t454) * qJD(4) + t360) * t273 + t304, t397 + t30 * qJD(4) + (t272 * t419 + t280) * t271, t364 + (t270 * t366 + t271 * t365) * m(7); qJD(2) * t10 - qJD(5) * t1 + qJD(6) * t21 - t308, t392, -qJD(5) * t29 - t304, qJD(5) * t17 - qJD(6) * t76, qJD(5) * t270 * t323 + pkin(8) * t280 + t211 * qJD(6) - t326 * t365 - t291, qJD(5) * t211 - t303; qJD(2) * t12 - qJD(3) * t8 + qJD(4) * t1 + qJD(6) * t35 + t398, t389, qJD(4) * t29 - t397, t291, t254 * qJD(6), t301; qJD(2) * t361 - qJD(3) * t108 - qJD(4) * t21 - qJD(5) * t35 - t393, qJD(1) * t361, -t364, t303, -t301, 0;];
Cq  = t15;
